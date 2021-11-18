#pragma once

#include <algorithm>

#include "meshtype.h"


class DisjointSet {
	uint8_t *rank;
public:
	int *parent;
	const int inf = 0x7fffffff;
	DisjointSet(int N) {
		parent = new int[N];
		rank = new uint8_t[N];
		for (int i = 0; i < N; i++) {
			parent[i] = -inf;
			rank[i] = 0;
		}
	}
	~DisjointSet() {
		delete parent; parent = 0;
		delete rank; rank = 0;
	}
	int findRepresentative(int i) {
		if (parent[i] < 0) return i;
		else {
			int ans = findRepresentative(parent[i]);
			parent[i] = ans;
			return ans;
		}
	}
	int representative_ID(int i) {
		while (parent[i] >= 0) i = parent[i];
		return -1 - parent[i];
	}
	bool unionSet(int i, int j) {
		int i_rep = findRepresentative(i);
		int j_rep = findRepresentative(j);
		if (i_rep == j_rep) return false;
		if (rank[i_rep] < rank[j_rep])
			parent[i_rep] = parent[i] = j_rep;
		else if (rank[i_rep] > rank[j_rep])
			parent[j_rep] = parent[j] = i_rep;
		else parent[j_rep] = parent[j] = i_rep, rank[i_rep]++;
		return true;
	}
};

void TrigsToMesh(const std::vector<triangle_3d> &trigs, float epsilon,
	std::vector<vec3> &vertice, std::vector<ivec3> &faces) {

#if 0
		{
			int FN = (int)trigs.size();
			vertice.resize(3 * FN);
			faces.resize(FN);
			for (int i = 0; i < FN; i++) {
				for (int u = 0; u < 3; u++) {
					vertice[3 * i + u] = trigs[i][u];
					((int*)&faces[i])[u] = 3 * i + u;
				}
			}
			return;
		}
#endif

		int FN = (int)trigs.size();
		int VN = 0;

		// restore vertice
		struct vec3_id {
			vec3 p;
			int id;
		};
		vec3_id *vtx = new vec3_id[3 * FN];
		for (int i = 0; i < FN; i++) {
			for (int u = 0; u < 3; u++)
				vtx[3 * i + u] = vec3_id{ trigs[i][u], 3 * i + u };
		}
		vertice.clear();
		faces.resize(FN);

		if (!(epsilon > 0.)) {

			std::sort(vtx, vtx + 3 * FN, [](vec3_id a, vec3_id b) {
				return a.p.x < b.p.x ? true : a.p.x > b.p.x ? false : a.p.y < b.p.y ? true : a.p.y > b.p.y ? false : a.p.z < b.p.z;
			});

			vec3 previous_p = vec3(NAN);
			for (int i = 0; i < 3 * FN; i++) {
				if (vtx[i].p != previous_p) {
					previous_p = vtx[i].p;
					vertice.push_back(vtx[i].p);
					VN++;
				}
				((int*)&faces[vtx[i].id / 3])[vtx[i].id % 3] = VN - 1;
			}
		}

		else {
			DisjointSet dsj(3 * FN);

			// three level sorting
			std::sort(vtx, vtx + 3 * FN, [](vec3_id a, vec3_id b) { return a.p.z < b.p.z; });
			for (int i = 0; i < 3 * FN;) {
				int j = i + 1;
				while (j < 3 * FN && vtx[j].p.z - vtx[j - 1].p.z < epsilon) j++;
				std::sort(vtx + i, vtx + j, [](vec3_id a, vec3_id b) { return a.p.y < b.p.y; });
				for (int u = i; u < j;) {
					int v = u + 1;
					while (v < j && vtx[v].p.y - vtx[v - 1].p.y < epsilon) v++;
					std::sort(vtx + u, vtx + v, [](vec3_id a, vec3_id b) { return a.p.x < b.p.x; });
					for (int m = u; m < v;) {
						int n = m + 1;
						while (n < v && vtx[n].p.x - vtx[n - 1].p.x < epsilon) n++;
						if (1) {  // O(N)
							for (int t = m; t + 1 < n; t++) dsj.unionSet(vtx[t].id, vtx[t + 1].id);
						}
						else {  // O(N²), more accurate, slower
							for (int t1 = m; t1 < n; t1++) for (int t2 = m; t2 < t1; t2++) {
								if (length(vtx[t2].p - vtx[t1].p) < epsilon) dsj.unionSet(vtx[t1].id, vtx[t2].id);
							}
						}
						m = n;
					}
					u = v;
				}
				i = j;
			}

			// pull points out from the disjoint set
			int unique_count = 0;
			int *vertice_map = new int[3 * FN];
			for (int i = 0; i < 3 * FN; i++)
				if (dsj.findRepresentative(i) == i) {
					vertice_map[i] = unique_count++;
					vertice.push_back(trigs[i / 3][i % 3]);
				}
			for (int i = 0; i < 3 * FN; i++)
				vertice_map[i] = vertice_map[dsj.findRepresentative(i)];
			for (int i = 0; i < FN; i++) for (int u = 0; u < 3; u++) {
				((int*)&faces[i])[u] = vertice_map[3 * i + u];
			}
			delete vertice_map;

			// remove degenerate faces
			std::vector<ivec3> face_bk = faces;
			faces.clear();
			for (ivec3 f : face_bk) {
				if (f.x == f.y || f.x == f.z || f.y == f.z) continue;
				faces.push_back(f);
			}
		}

		delete vtx;

}
