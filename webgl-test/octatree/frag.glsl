#version 300 es
precision highp float;

out vec4 fragColor;

uniform int iFrame;
uniform float iRz;
uniform float iRx;
uniform float iSc;
uniform vec2 iResolution;

#define TREEBUFFER_SIZE 1024
uniform mediump usampler2D uTreeBuffer;

#define ZERO min(iFrame,0)

#define P0 vec3(-2.0, -2.0, -2.0) /* min coordinates of grid */
#define P1 vec3(2.0, 2.0, 2.0) /* max coordinates of grid */
#define GRID_DIF ivec3(1, 1, 1) /* initial grid size, at least one odd component */
#define PLOT_DEPTH 6 /* depth of the tree */
#define GRID_SIZE (GRID_DIF*(1<<PLOT_DEPTH))
#define EDGE_ROUNDING 255 /* divide edge into # intervals and round to integer coordinate */
#define MESH_SIZE (GRID_SIZE*EDGE_ROUNDING)


int getUint8(int i) {
    ivec2 pos = ivec2((i/4)%TREEBUFFER_SIZE, (i/4)/TREEBUFFER_SIZE);
    uvec4 sp = texelFetch(uTreeBuffer, pos, 0);
    return int(i%4==0 ? sp.x : i%4==1 ? sp.y : i%4==2 ? sp.z : sp.w);
}
int getUint32(int i) {
	//return getUint8(i) + 256 * (getUint8(i + 1) + 256 * (getUint8(i + 2) + 256 * getUint8(i + 3)));
    ivec2 pos = ivec2((i/4)%TREEBUFFER_SIZE, (i/4)/TREEBUFFER_SIZE);
    ivec4 sp = ivec4(texelFetch(uTreeBuffer, pos, 0));
	return sp.x + 256 * (sp.y + 256 * (sp.z + 256 * sp.w));
}
ivec3 getUvec3(int i) {
	int x = getUint8(i);
	int y = getUint8(i + 1);
	int z = getUint8(i + 2);
	return ivec3(x, y, z);
}


const ivec3 VERTEX_LIST[8] = ivec3[8](
	ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
	ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
);

vec2 intersectBox(vec3 r, vec3 ro, vec3 inv_rd) {  // inv_rd = 1/rd
	vec3 p = -inv_rd * ro;
	vec3 k = abs(inv_rd)*r;
	vec3 t1 = p - k, t2 = p + k;
	float tN = max(max(t1.x, t1.y), t1.z);
	float tF = min(min(t2.x, t2.y), t2.z);
	if (tN >= tF || tF < 0.0) return vec2(-1.0);
	return vec2(tN, tF);
}

float intersectTriangle(vec3 ro, vec3 rd, vec3 v01, vec3 v02) {
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


// faster in CPU, faster in Chrome for deep trees, slower in Firefox
#define IN_DISTANCE_ORDER 1

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
	vec3 inv_rd = 1.0 / rd;
	vec2 ib = intersectBox(0.5*(P1 - P0), ro - 0.5*(P0 + P1), inv_rd);
	if (ib.y <= 0.0 || ib.x > t1) return false;

	// calculate the order to traverse subcells
	int subcell_order[8] = int[8](0, 1, 2, 3, 4, 5, 6, 7);
#if IN_DISTANCE_ORDER
	float dist[8];
	for (int i = 0; i < 8; i++) {
		dist[i] = dot(vec3(VERTEX_LIST[i]), rd);
	}
	for (int i = 1; i < 8; i++) {  // sorting
		int soi = subcell_order[i];
		float di = dist[i];
		int j = i - 1;
		while (j >= 0 && dist[j] > di) {
			dist[j + 1] = dist[j];
			subcell_order[j + 1] = subcell_order[j];
			j--;
		}
		dist[j + 1] = di;
		subcell_order[j + 1] = soi;
	}
#endif

	// debug
	int loop_count = 0;
	int trig_int_count = 0;
	int box_int_count = 0;

	// grid
	for (int xi = ZERO; xi < GRID_DIF.x; xi++) for (int yi = ZERO; yi < GRID_DIF.y; yi++) for (int zi = ZERO; zi < GRID_DIF.z; zi++) {
		int grid_pos = getUint32(4 * ((zi * GRID_DIF.y + yi) * GRID_DIF.x + xi));
		if (grid_pos == 0) continue;
		float cur_t1 = min_t;

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
			loop_count++;

			// test if current node is none
			if (cur.ptr != 0) {
				box_int_count++;
				vec3 c = mix(P0, P1, (vec3(cur.pos) + 0.5 * float(cell_size)) / vec3(GRID_SIZE));
				vec3 r = (P1 - P0) / vec3(GRID_SIZE) * 0.5 * float(cell_size);
				ib = intersectBox(r, ro - c, inv_rd);
				if (!(ib.y > 0.0 && ib.x < min_t)) cur.ptr = 0;
			}

			// go into subtree
			if (cur.ptr != 0) {
				// triangles
				if (cell_size == 1) {
					int n = getUint8(cur.ptr);
					ivec3 po = cur.pos * EDGE_ROUNDING;
					for (int ti = 0; ti < n; ti++) {
						trig_int_count++;
						ivec3 vi0 = getUvec3(cur.ptr + 12 * ti + 1);
						ivec3 vi1 = getUvec3(cur.ptr + 12 * ti + 4);
						ivec3 vi2 = getUvec3(cur.ptr + 12 * ti + 7);
						vec3 v0 = mix(P0, P1, vec3(po + vi0) / vec3(MESH_SIZE));
						vec3 v01 = ((P1 - P0) / vec3(MESH_SIZE)) * vec3(vi1 - vi0);
						vec3 v02 = ((P1 - P0) / vec3(MESH_SIZE)) * vec3(vi2 - vi0);
						t = intersectTriangle(ro - v0, rd, v01, v02);
						if (t > 0.0 && t < min_t) {
							min_t = t, min_n = cross(v01, v02);
							col = vec3(getUvec3(cur.ptr + 12 * ti + 10)) / 255.0;
						}
					}
					cur.ptr = 0;
				}
				// subtree
				else {
					stk[++si] = cur; cell_size /= 2;
					cur.subcell = 0;
					int subcell = subcell_order[cur.subcell];
					cur.ptr = getUint32(cur.ptr + 4 * subcell);
					cur.pos += VERTEX_LIST[subcell] * cell_size;
				}
			}

			// next node
			else if (si != -1) {
				cur = stk[si--]; cell_size *= 2;
				cur.subcell += 1;
				if (cur.subcell >= 8) {
					cur.ptr = 0;
				}
				else {
					stk[++si] = cur; cell_size /= 2;
					int subcell = subcell_order[cur.subcell];
					cur.pos = cur.pos + VERTEX_LIST[subcell] * cell_size;
					cur.ptr = getUint32(cur.ptr + 4 * subcell);
					cur.subcell = 0;
				}
			}

			else break;
#if IN_DISTANCE_ORDER
			if (min_t < cur_t1) break;
#endif
		}
	}

	//col = vec3(loop_count, box_int_count, trig_int_count) / 255.0;
	//col = vec3(box_int_count) / 255.0;
	//col = vec3(1.0) * float(box_int_count) / 255.0;

	return min_t < t1;
}


#if 1
// fast preview
vec3 mainRender(vec3 ro, vec3 rd) {
	float t;
	vec3 n, col = vec3(0.0);
	if (intersectObject(ro, rd, t, 1e6, n, col)) {
		//return col;
		return col * abs(dot(normalize(n), rd));
	}
	return col;
}
#else
// ray tracing
vec3 mainRender(vec3 ro, vec3 rd) {
	const int MAT_BACKGROUND = 0;
	const int MAT_PLANE = 1;
	const int MAT_OBJECT = 2;

	vec3 m_col = vec3(1.0), t_col = vec3(0.0);
	bool inside_object = false;

	for (int iter = ZERO; iter < 4; iter++) {
		rd = normalize(rd);
		ro += 1e-4 * rd;
		float t, min_t = 1e12;
		vec3 n, min_n;
		int material = MAT_BACKGROUND;
		vec3 col = vec3(0.0);

		// plane
		t = -(ro.z - (P0.z - 0.01)) / rd.z;
		if (t > 0.0) {
			min_t = t, min_n = vec3(0, 0, 1);
			material = MAT_PLANE;
		}

		// object
		t = min_t;
		if (intersectObject(ro, rd, t, min_t, n, col) && t < min_t) {
			min_t = t, min_n = normalize(n);
			material = MAT_OBJECT;
		}

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
#endif

void mainImage(out vec4 fragColor, vec2 fragCoord) {
	const vec3 CENTER = vec3(0, 0, 0);
	const float DIST = 20.0;  // larger = smaller
	const float VIEW_FIELD = 0.4;  // larger = larger + more perspective

	// camera
	vec3 w = vec3(cos(iRx)*vec2(cos(iRz), sin(iRz)), sin(iRx));
	vec3 u = vec3(-sin(iRz), cos(iRz), 0);
	vec3 v = cross(w, u);
	vec3 ro = DIST * w + CENTER;
	vec2 uv = iSc * (2.0*(fragCoord.xy - 0.5) / iResolution.xy - 1.0);
	vec2 sc = iResolution.xy / min(iResolution.x, iResolution.y);
	vec3 rd = mat3(u, v, -w)*vec3(VIEW_FIELD*uv*sc, 1.0);
	rd = normalize(rd);

	// calculate pixel color
	vec3 col = mainRender(ro, rd);
	fragColor = vec4(col, 1.0);
}

void main(void) {
    mainImage(fragColor, gl_FragCoord.xy);
}
