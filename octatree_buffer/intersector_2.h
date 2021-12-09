// raymarching on octree

#include "test.cpp"


const ivec3 VERTEX_LIST[8] = {
	ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
	ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
};

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



// ray-object intersection, grid/tree lookup
bool intersectObject(vec3 ro, vec3 rd, float &min_t, float t1, vec3 &min_n, vec3 &col) {
	float t;
	min_t = t1;

	// bounding box
	vec3 inv_rd = 1.0 / rd;
	vec2 ib = intersectBox(0.5*(P1 - P0), ro - 0.5*(P0 + P1), inv_rd);
	if (ib.y <= 0.0 || ib.x > t1) return false;

	// calculate the order to traverse subcells
	int subcell_order[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
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

	// debug
	int loop_count = 0;
	int trig_int_count = 0;
	int box_int_count = 0;

	// grid
	for (int xi = ZERO; xi < GRID_DIF.x; xi++) for (int yi = ZERO; yi < GRID_DIF.y; yi++) for (int zi = ZERO; zi < GRID_DIF.z; zi++) {
		int grid_pos = getUint32(4 * ((zi * GRID_DIF.y + yi) * GRID_DIF.x + xi));
		if (grid_pos == 0) continue;

		// intersect grid cell
		box_int_count++;
		vec3 inv_rd = 1.0 / rd;
		vec3 c = mix(P0, P1, (vec3(xi, yi, zi) + 0.5) / vec3(GRID_DIF));
		vec3 r = 0.5 * (P1 - P0) / vec3(GRID_DIF);
		vec2 t01 = intersectBox(r, ro - c, inv_rd);
		if (t01.y <= 0.0 || t01.x > min_t) continue;
		float epsilon = max(0.01f * min(min(r.x, r.y), r.z) * exp2(-float(PLOT_DEPTH)), 1e-4f);

		// raymarching
		float t0 = t01.x;
		float t_max = min(t01.y, min_t) - epsilon;
		while (!(t0 >= t_max)) {
			loop_count++;
			if (loop_count > 1024) { col = vec3(1, 0, 0); return true; }

			// find the next cell
			int cell_size = 1 << PLOT_DEPTH;
			ivec3 cell_pos = ivec3(xi, yi, zi) * cell_size;
			int ptr = grid_pos;
			for (int si = 0; si < PLOT_DEPTH; si++) {
				cell_size /= 2;
				ivec3 cell_pos_1 = cell_pos;
				for (int j = 0; j < 8; j++) {
					int i = subcell_order[j];
					box_int_count++;
					cell_pos_1 = cell_pos + VERTEX_LIST[i] * cell_size;
					vec3 c = mix(P0, P1, (vec3(cell_pos_1) + 0.5 * float(cell_size)) / vec3(GRID_SIZE));
					vec3 r = 0.5 * float(cell_size) * (P1 - P0) / vec3(GRID_SIZE);
					vec2 t01_t = intersectBox(r, ro - c, inv_rd);
					if (t01_t.y > t0 + epsilon) {
						t01 = t01_t;
						ptr = getUint32(ptr + 4 * i);
						break;
					}
				}
				cell_pos = cell_pos_1;
				if (ptr == 0) break;
			}

			// triangles
			if (ptr != 0) {
				int n = getUint8(ptr);
				ivec3 po = cell_pos * EDGE_ROUNDING;
				bool triangle_found = false;
				for (int ti = 0; ti < n; ti++) {
					trig_int_count++;
					ivec3 vi0 = getUvec3(ptr + 12 * ti + 1);
					ivec3 vi1 = getUvec3(ptr + 12 * ti + 4);
					ivec3 vi2 = getUvec3(ptr + 12 * ti + 7);
					vec3 v0 = mix(P0, P1, vec3(po + vi0) / vec3(MESH_SIZE));
					vec3 v01 = ((P1 - P0) / vec3(MESH_SIZE)) * vec3(vi1 - vi0);
					vec3 v02 = ((P1 - P0) / vec3(MESH_SIZE)) * vec3(vi2 - vi0);
					t = intersectTriangle(ro - v0, rd, v01, v02);
					if (t > 0.0 && t < min_t) {
						min_t = t, min_n = cross(v01, v02);
						col = vec3(getUvec3(ptr + 12 * ti + 10)) / 255.0;
						triangle_found = true;
					}
				}
				if (triangle_found) break;
			}

			// march
			t0 = t01.y;
		}
	}

	//col = vec3(loop_count, box_int_count, trig_int_count) / 255.0;
	//col = vec3(box_int_count) / 255.0;
	//col = vec3(1.0) * float(box_int_count) / 255.0;

	return min_t < t1;
}
