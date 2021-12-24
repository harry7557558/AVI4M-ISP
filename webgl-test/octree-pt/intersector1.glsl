// stack-based octree traversal


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
bool intersectOctree(usampler2D treeSampler, vec3 P0, vec3 P1, ivec3 GRID_DIF, int PLOT_DEPTH, int EDGE_ROUNDING,
    vec3 ro, vec3 rd, out float min_t, float t1, out vec3 min_n, out vec3 col) {
    ivec3 GRID_SIZE = GRID_DIF*(1<<PLOT_DEPTH);
    ivec3 MESH_SIZE = (GRID_SIZE*EDGE_ROUNDING);

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

	// grid
	for (int xi = int(ZERO); xi < GRID_DIF.x; xi++) for (int yi = int(ZERO); yi < GRID_DIF.y; yi++) for (int zi = int(ZERO); zi < GRID_DIF.z; zi++) {
		int grid_pos = getUint32(treeSampler, 4 * ((zi * GRID_DIF.y + yi) * GRID_DIF.x + xi));
		if (grid_pos == 0) continue;
		float cur_t1 = min_t;

		// tree traversal
		TreeNode stk[12];  // stack
		int si = -1;  // index of the top of the stack
		int cell_size = 1 << PLOT_DEPTH;
		TreeNode cur;  // current node
		cur.pos = ivec3(xi, yi, zi)*cell_size;
		cur.subcell = 0;
		cur.ptr = grid_pos;
		vec3 p0, p1;

		while (true) {

			// test if current node is none
			if (cur.ptr != 0) {
				vec3 c = mix(P0, P1, (vec3(cur.pos) + 0.5 * float(cell_size)) / vec3(GRID_SIZE));
				vec3 r = (P1 - P0) / vec3(GRID_SIZE) * 0.5 * float(cell_size);
				ib = intersectBox(r, ro - c, inv_rd);
				if (!(ib.y > 0.0 && ib.x < min_t)) cur.ptr = 0;
			}

			// go into subtree
			if (cur.ptr != 0) {
				// triangles
				if (cell_size == 1) {
					int n = getUint8(treeSampler, cur.ptr);
					ivec3 po = cur.pos * EDGE_ROUNDING;
					for (int ti = 0; ti < n; ti++) {
						ivec3 vi0 = getUvec3(treeSampler, cur.ptr + 12 * ti + 1);
						ivec3 vi1 = getUvec3(treeSampler, cur.ptr + 12 * ti + 4);
						ivec3 vi2 = getUvec3(treeSampler, cur.ptr + 12 * ti + 7);
						vec3 v0 = mix(P0, P1, vec3(po + vi0) / vec3(MESH_SIZE));
						vec3 v01 = ((P1 - P0) / vec3(MESH_SIZE)) * vec3(vi1 - vi0);
						vec3 v02 = ((P1 - P0) / vec3(MESH_SIZE)) * vec3(vi2 - vi0);
						t = intersectTriangle(ro - v0, rd, v01, v02);
						if (t > 0.0 && t < min_t) {
							min_t = t, min_n = cross(v01, v02);
							col = vec3(getUvec3(treeSampler, cur.ptr + 12 * ti + 10)) / 255.0;
						}
					}
					cur.ptr = 0;
				}
				// subtree
				else {
					stk[++si] = cur; cell_size /= 2;
					cur.subcell = 0;
					int subcell = subcell_order[cur.subcell];
					cur.ptr = getUint32(treeSampler, cur.ptr + 4 * subcell);
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
					cur.ptr = getUint32(treeSampler, cur.ptr + 4 * subcell);
					cur.subcell = 0;
				}
			}

			else break;
#if IN_DISTANCE_ORDER
			if (min_t < cur_t1) break;
#endif
		}
	}

	return min_t < t1;
}
