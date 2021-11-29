#version 300 es
precision highp float;

out vec4 fragColor;

uniform float iRz;
uniform float iRx;
uniform float iSc;
uniform vec2 iResolution;

#define TREEBUFFER_SIZE 1024
uniform highp usampler2D uTreeBuffer;


#define P0 vec3(-2.0)
#define P1 vec3(2.0)
#define SEARCH_DIF ivec3(4, 4, 4)
#define PLOT_DEPTH 2
#define GRID_SIZE (SEARCH_DIF*(1<<PLOT_DEPTH))
#define EDGE_ROUNDING 128
#define MESH_SIZE (GRID_SIZE*EDGE_ROUNDING)


int getUint8(int i) {
    ivec2 pos = ivec2((i/4)%TREEBUFFER_SIZE, (i/4)/TREEBUFFER_SIZE);
    uvec4 sp = texelFetch(uTreeBuffer, pos, 0);
    return int(i%4==0 ? sp.x : i%4==1 ? sp.y : i%4==2 ? sp.z : sp.w);
}
int getUint16(int i) {
	return getUint8(i) + getUint8(i + 1) * 256;
}
int getUint32(int i) {
	return getUint8(i) + 256 * (getUint8(i + 1) + 256 * (getUint8(i + 2) + 256 * getUint8(i + 3)));
}
ivec3 getUvec3(int i) {
	int x = getUint8(i);
	int y = getUint8(i + 1);
	int z = getUint8(i + 2);
	return ivec3(x, y, z);
}


float intersectBox(vec3 r, vec3 ro, vec3 inv_rd) {  // inv_rd = 1/rd
	vec3 p = -inv_rd * ro;
	vec3 k = abs(inv_rd)*r;
	vec3 t1 = p - k, t2 = p + k;
	float tN = max(max(t1.x, t1.y), t1.z);
	float tF = min(min(t2.x, t2.y), t2.z);
	if (tN > tF || tF < 0.0) return -1.0;
	return tN < 0.0 ? tF : tN;
}

float intersectTriangle(vec3 ro, vec3 rd, vec3 v0, vec3 v1, vec3 v2) {
	ro -= v0;
	vec3 v01 = v1 - v0, v02 = v2 - v0;
	vec3 n = cross(v01, v02);
	vec3 q = cross(ro, rd);
	float d = 1.0 / dot(rd, n);
	float u = -d * dot(q, v02);
	if (u < 0.0 || u > 1.0) return -1.0;
	float v = d * dot(q, v01);
	if (v < 0.0 || u + v > 1.0) return -1.0;
	float t = -d * dot(n, ro);
	return t;
}


ivec3 VERTEX_LIST[8] = ivec3[8](
	ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
	ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
);

// used for tree traversal
struct TreeNode {
	ivec3 pos;  // position of origin of the cell
	int subcell;  // ID of subcell during traversal
	int ptr;  // position in the buffer
};

// ray-object intersection, grid/tree lookup
bool intersectObject(vec3 ro, vec3 rd, out float min_t, float t1, out vec3 min_n, out vec3 col) {
	float t;
	min_t = t1;

	// bounding box
	t = intersectBox(0.5*(P1 - P0), ro - 0.5*(P0 + P1), 1.0 / rd);
	if (t <= 0.0 || t > t1) return false;

	// grid
	for (int xi = 0; xi < SEARCH_DIF.x; xi++) for (int yi = 0; yi < SEARCH_DIF.y; yi++) for (int zi = 0; zi < SEARCH_DIF.z; zi++) {
		int grid_pos = getUint32(4 * ((zi * SEARCH_DIF.y + yi) * SEARCH_DIF.x + xi));
		if (grid_pos == 0) continue;

		// tree traversal
		TreeNode stk[PLOT_DEPTH];  // stack
		int si = -1;  // index of the top of the stack
		int cell_size = 1 << PLOT_DEPTH;
		TreeNode cur;  // current node
		cur.pos = ivec3(xi, yi, zi)*cell_size;
		cur.subcell = 0;
		cur.ptr = grid_pos;
		vec3 p0, p1;

		while (true) {

			// test of current node is none
			if (cur.ptr != 0) {
				p0 = mix(P0, P1, vec3(cur.pos) / vec3(GRID_SIZE));
				p1 = mix(P0, P1, vec3(cur.pos + cell_size) / vec3(GRID_SIZE));
				vec3 r = 0.5*(p1 - p0), c = 0.5*(p0 + p1);
				t = intersectBox(r, ro - c, 1.0 / rd);
				if (!(t > 0.0 && t < min_t)) cur.ptr = 0;
			}

			// go into subtree
			if (cur.ptr != 0) {
				// triangles
				if (cell_size == 1) {
					int n = getUint8(cur.ptr);
					ivec3 po = cur.pos * EDGE_ROUNDING;
					for (int ti = 0; ti < n; ti++) {
						vec3 a = mix(P0, P1, vec3(po + getUvec3(cur.ptr + 12 * ti + 1)) / vec3(MESH_SIZE));
						vec3 b = mix(P0, P1, vec3(po + getUvec3(cur.ptr + 12 * ti + 4)) / vec3(MESH_SIZE));
						vec3 c = mix(P0, P1, vec3(po + getUvec3(cur.ptr + 12 * ti + 7)) / vec3(MESH_SIZE));
						t = intersectTriangle(ro, rd, a, b, c);
						if (t > 0.0 && t < min_t) {
							min_t = t, min_n = cross(b - a, c - a);
							col = vec3(getUvec3(cur.ptr + 12 * ti + 10)) / 255.0;
						}
					}
					cur.ptr = 0;
				}
				// subtree
				else {
					stk[++si] = cur, cell_size /= 2;
					cur.subcell = 0;
					cur.ptr = getUint32(cur.ptr);
				}
			}

			// next node
			else if (si != -1) {
				cur = stk[si--], cell_size *= 2;
				cur.subcell += 1;
				if (cur.subcell >= 8) {
					cur.ptr = 0;
				}
				else {
					stk[++si] = cur, cell_size /= 2;
					cur.pos = cur.pos + VERTEX_LIST[cur.subcell] * cell_size;
					cur.ptr = getUint32(cur.ptr + 4 * cur.subcell);
					cur.subcell = 0;
				}
			}

			else break;
		}
	}

	return min_t < t1;
}

vec3 mainRender(vec3 ro, vec3 rd) {
	const int MAT_BACKGROUND = 0;
	const int MAT_PLANE = 1;
	const int MAT_OBJECT = 2;

	vec3 m_col = vec3(1.0), t_col = vec3(0.0), col;
	bool inside_object = false;

	for (int iter = 0; iter < 64; iter++) {
		rd = normalize(rd);
		ro += 1e-4 * rd;
		float t, min_t = 1e12;
		vec3 n, min_n;
		int material = MAT_BACKGROUND;

		// plane
		t = -(ro.z - P0.z) / rd.z;
		if (t > 0.0) {
			min_t = t, min_n = vec3(0, 0, 1);
			material = MAT_PLANE;
		}

		// object
		t = min_t;
		if (intersectObject(ro, rd, t, min_t, n, col) && t < min_t) {
			min_t = t, min_n = normalize(n);
			material = MAT_OBJECT;
			//return col * max(dot(min_n, -rd), 0.0);
		}
		//return vec3(0, 0, 0);

		// update ray
		if (material == MAT_BACKGROUND) {
			col = vec3(1.0)*max(rd.z, 0.0);
			return m_col * col + t_col;
		}
		min_n = dot(rd, min_n) < 0.0 ? min_n : -min_n;  // ray hits into the surface
		ro += rd * min_t;
		if (material == MAT_PLANE) {
			m_col *= vec3(0.8, 0.9, 1.0);
			rd = rd - 2.0*dot(rd, min_n)*min_n;
		}
		if (material == MAT_OBJECT) {
			m_col *= col;
			rd = rd - 2.0*dot(rd, min_n)*min_n;
		}
		if (inside_object) return 1e12*vec3(1, -1, -1);  // red warning
	}
	return m_col + t_col;
}

void mainImage(out vec4 fragColor, vec2 fragCoord) {
    int i = int(fragCoord.x) + int(fragCoord.y) * int(iResolution.x);
    int pos = getUint16(i);
    //pos = ivec2(fragCoord);
	vec4 s = vec4(pos) / 256.0;
    //s = vec4(i) / 65536.0;
    fragColor = vec4(s.x, s.y, s.z, 1.0);
    //return;

	const float SCALE = 1.0;  // larger = smaller (more view field)
	const vec3 CENTER = vec3(0, 0, 0);
	const float DIST = 20.0;  // larger = smaller
	const float VIEW_FIELD = 0.4;  // larger = larger + more perspective

	// camera
	vec3 w = vec3(cos(iRx)*vec2(cos(iRz), sin(iRz)), sin(iRx));
	vec3 u = vec3(-sin(iRz), cos(iRz), 0);
	vec3 v = cross(w, u);
	vec3 ro = DIST * w + CENTER;
	vec2 uv = iSc * SCALE * (2.0*(fragCoord.xy - 0.5) / iResolution.xy - 1.0);
	vec2 sc = iResolution.xy / length(iResolution.xy);
	vec3 rd = mat3(u, v, -w)*vec3(VIEW_FIELD*uv*sc, 1.0);
	rd = normalize(rd);

	// calculate pixel color
	vec3 col = mainRender(ro, rd);
	fragColor = vec4(col, 1.0);
}

void main(void) {
    mainImage(fragColor, gl_FragCoord.xy);
}
