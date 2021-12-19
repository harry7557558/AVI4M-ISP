// Attempt to find an adaptive marching cube look up table that
// generates a mesh to best represent the isosurface.

// References:
// http://paulbourke.net/geometry/polygonise/ (indices of vertices/edges)
// https://www.cs.upc.edu/~virtual/SGI/docs/1.%20Theory/Unit%2010.%20Volume%20models.%20Marching%20Cubes/Marching%20Cubes.pdf


#include <vector>
#include <initializer_list>
#include <assert.h>

#include "../octatree/trigs2mesh.h"
#include "../octatree/ply_writer.h"

std::vector<std::vector<ivec3>> completeTriangulation(int cubeIndex, int faceIndex[6]);


/* cube elements */

const ivec3 VERTICE_LIST[8] = {
	ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
	ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
};
const int VERTICE_LIST_INV[2][2][2] = {
	{{0, 4}, {1, 5}}, {{3, 7}, {2, 6}}
};

const ivec2 EDGE_LIST[12] = {
	ivec2(0,1), ivec2(1,2), ivec2(2,3), ivec2(3,0),
	ivec2(4,5), ivec2(5,6), ivec2(6,7), ivec2(7,4),
	ivec2(0,4), ivec2(1,5), ivec2(2,6), ivec2(3,7)
};

const int FACE_LIST[6][4] = {  // list vertices for each face
	{0, 4, 5, 1}, {0, 3, 7, 4}, {0, 1, 2, 3}, // ccw
	{3, 7, 6, 2}, {1, 2, 6, 5}, {4, 5, 6, 7}  // cw
};
const ivec3 FACE_DIR[6] = {
	ivec3(-1,0,0), ivec3(0,-1,0), ivec3(0,0,-1),
	ivec3(1,0,0), ivec3(0,1,0), ivec3(0,0,1)
};
const int FACE_EDGE_LIST[6][4] = {  // list edges for each face
	{8, 4, 9, 0}, {3, 11, 7, 8}, {0, 1, 2, 3},
	{11, 6, 10, 2}, {1, 10, 5, 9}, {4, 5, 6, 7}
};

// face elements
const ivec2 VERTICE_LIST_FACE[4] = {
	ivec2(0,0), ivec2(1,0), ivec2(1,1), ivec2(0,1)
};
const int VERTICE_LIST_FACE_INV[2][2] = {
	{0, 3}, {1, 2}
};
const ivec2 EDGE_LIST_FACE[4] = {
	ivec2(0,1), ivec2(1,2), ivec2(2,3), ivec2(3,0)
};
const int LUT_FACE[18][4] = {  // left is negative (1)
	{ -1, }, // 0000
	{ 0, 3, -1 }, // 1000
	{ 1, 0, -1 }, // 0100
	{ 1, 3, -1 }, // 1100
	{ 2, 1, -1 }, // 0010
	{ 0, 3, 2, 1 }, // 1010, #5
	{ 2, 0, -1 }, // 0110
	{ 2, 3, -1 }, // 1110
	{ 3, 2, -1 }, // 0001
	{ 0, 2, -1 }, // 1001
	{ 1, 0, 3, 2 }, // 0101, #10
	{ 1, 2, -1 }, // 1101
	{ 3, 1, -1 }, // 0011
	{ 0, 1, -1 }, // 1011
	{ 3, 0, -1 }, // 0111
	{ -1 }, // 1111
	{ 0, 1, 2, 3 }, // 1010 alternate, #16
	{ 1, 2, 3, 0 }, // 0101 alternate, #17
};


/* value to index */

int calcIndex(const float v[8]) {
	if (isnan(v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7]))
		return 0;
	return int(v[0] <= 0) |
		(int(v[1] <= 0) << 1) |
		(int(v[2] <= 0) << 2) |
		(int(v[3] <= 0) << 3) |
		(int(v[4] <= 0) << 4) |
		(int(v[5] <= 0) << 5) |
		(int(v[6] <= 0) << 6) |
		(int(v[7] <= 0) << 7);
}
void toIndex(const int index, float v[8]) {
	for (int i = 0; i < 8; i++)
		v[i] = (index >> i) & 1 ? -1.0f : 1.0f;
}
int calcIndex(const std::initializer_list<float> c) {
	return calcIndex(&*c.begin());
}
int calcFaceIndex(const float v[4]) {
	if (isnan(v[0] + v[1] + v[2] + v[3]))
		return 0;
	return int(v[0] <= 0) |
		(int(v[1] <= 0) << 1) |
		(int(v[2] <= 0) << 2) |
		(int(v[3] <= 0) << 3);
}


/* cube cases */

struct CubeCase {
	float v[8];
	int index = -1;
	int case_id = -1;  // 0-14
	bool isAmbiguous = false;
	int faceIndex[6] = { -1, -1, -1, -1, -1, -1 };
	bool isAmbiguousFace[6] = { false, false, false, false, false, false };
	bool hasAmbiguousFace = false;
	void calcFaceIndices() {
		hasAmbiguousFace = false;
		for (int i = 0; i < 6; i++) {
			const int* verts = FACE_LIST[i];
			float vf[4] = { v[verts[0]], v[verts[1]], v[verts[2]], v[verts[3]] };
			faceIndex[i] = calcFaceIndex(vf);
			isAmbiguousFace[i] = faceIndex[i] == 5 || faceIndex[i] == 10;
			if (isAmbiguousFace[i]) hasAmbiguousFace = true;
		}
	}
	CubeCase() {};
	CubeCase(int index) {
		this->index = index;
		toIndex(index, this->v);
	}
	CubeCase(const std::initializer_list<float> c) {
		for (int i = 0; i < 8; i++) this->v[i] = *(c.begin() + i);
		this->index = calcIndex(v);
	}
} cubes[256];

void initCubeCases() {

	// 15 basic cases
	{
		int index;
		index = calcIndex({ 1, 1, 1, 1, 1, 1, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 0;
		index = calcIndex({ 0, 1, 1, 1, 1, 1, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 1;
		index = calcIndex({ 0, 1, 1, 0, 1, 1, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 2;
		index = calcIndex({ 0, 1, 1, 1, 1, 1, 1, 0 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 3;
		cubes[index].isAmbiguous = true;
		index = calcIndex({ 0, 1, 1, 1, 1, 1, 0, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 4;
		cubes[index].isAmbiguous = true;  // ambiguous cube
		index = calcIndex({ 1, 0, 0, 0, 1, 1, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 5;
		index = calcIndex({ 0, 1, 1, 0, 1, 1, 0, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 6;
		cubes[index].isAmbiguous = true;
		index = calcIndex({ 1, 1, 1, 0, 0, 1, 0, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 7;
		cubes[index].isAmbiguous = true;
		index = calcIndex({ 0, 0, 0, 0, 1, 1, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 8;
		index = calcIndex({ 0, 0, 0, 1, 1, 0, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 9;
		index = calcIndex({ 0, 1, 0, 1, 0, 1, 0, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 10;
		cubes[index].isAmbiguous = true;
		index = calcIndex({ 0, 0, 0, 1, 1, 1, 0, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 11;
		index = calcIndex({ 1, 0, 0, 0, 0, 1, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 12;
		cubes[index].isAmbiguous = true;
		index = calcIndex({ 0, 1, 0, 1, 1, 0, 1, 0 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 13;
		cubes[index].isAmbiguous = true;
		index = calcIndex({ 1, 0, 0, 0, 1, 0, 1, 1 });
		cubes[index] = CubeCase(index), cubes[index].case_id = 14;
	}
	std::vector<int> cube_ids;
	for (int i = 0; i < 256; i++) if (cubes[i].case_id != -1) {
		printf("%d ", i);
		cube_ids.push_back(i);
		cubes[i].calcFaceIndices();
		if (cubes[i].hasAmbiguousFace) assert(cubes[i].isAmbiguous);
		else if (cubes[i].isAmbiguous) assert(cubes[i].case_id == 4);
	}
	printf("\n");
}


/* brute-force iterate through all possible triangulations */

bool isValidTriangulation(
	const std::vector<int> &cube_edges,
	const std::vector<ivec2> &required_trig_edges,
	std::vector<ivec3> trigs  // ccw
) {
	// check signs
	// check if it contains all required edges
	// check edges of the triangles
	return true;
}

bool isIdenticalTriangle(ivec3 a, ivec3 b) {
	std::sort(&a.x, &a.x + 3);
	std::sort(&b.x, &b.x + 3);
	return a == b;
}
void _completeTriangulation_recurse(
	std::vector<std::vector<ivec3>> &valid_triangulations,
	const std::vector<int> &cube_edges,
	const std::vector<ivec2> &required_trig_edges,
	std::vector<ivec3> &trigs
) {
	// attempt to construct a new triangle
	for (ivec2 ab : required_trig_edges) for (int c : cube_edges) {
		printf("%d %d\n", ab.x, ab.y);
		if (c == ab.x || c == ab.y) continue;
		ivec3 t = ivec3(ab.x, c, ab.y);
		bool is_dup = false;
		for (ivec3 t1 : trigs)
			if (isIdenticalTriangle(t, t1)) { is_dup = true; break; }
		if (is_dup) continue;
		// recursion
		trigs.push_back(t);
		if (isValidTriangulation(cube_edges, required_trig_edges, trigs)) {
			valid_triangulations.push_back(trigs);
		}
		else {
			_completeTriangulation_recurse(
				valid_triangulations,
				cube_edges, required_trig_edges,
				trigs);
		}
		trigs.pop_back();
	}
}
std::vector<std::vector<ivec3>> completeTriangulation(
	int cubeIndex,
	int faceIndex[6]  // 0-17
) {
	// init + testing
	CubeCase cbc(cubeIndex);
	cbc.calcFaceIndices();
	for (int i = 0; i < 6; i++) {
		int fi1 = cbc.faceIndex[i], fi2 = faceIndex[i];
		assert(fi1 == fi2 || (fi1 == 5 && fi2 == 16) || (fi1 == 10 && fi2 == 17));
	}
	bool is_nonempty_edge[12];
	std::vector<int> cube_edges;
	for (int i = 0; i < 12; i++) {
		is_nonempty_edge[i] = cbc.v[EDGE_LIST[i].x] * cbc.v[EDGE_LIST[i].y] < 0.0f;
		if (is_nonempty_edge[i]) cube_edges.push_back(i);
	}
	std::vector<ivec2> required_trig_edges;
	for (int i = 0; i < 6; i++) {
		int fi = faceIndex[i];
		for (int j = 0; j < 4 && LUT_FACE[fi][j] != -1; j += 2) {
			int e1 = FACE_EDGE_LIST[i][LUT_FACE[fi][j]], e2 = FACE_EDGE_LIST[i][LUT_FACE[fi][j + 1]];
			if (i >= 3) std::swap(e1, e2);
			assert(is_nonempty_edge[e1] && is_nonempty_edge[e2]);
			required_trig_edges.push_back(ivec2(e1, e2));
		}
	}
	// recursively iterate through all triangulations
	std::vector<std::vector<ivec3>> res;
	std::vector<ivec3> trigs;
	_completeTriangulation_recurse(
		res,
		cube_edges, required_trig_edges,
		trigs
	);
	return res;
}


/* visualization */

void visualizeTriangulations(std::vector<std::vector<ivec3>> triangulations, const char* filepath) {
	std::vector<vec3> vertices;
	std::vector<ivec3> faces;
	auto addTriangle = [&](vec3 a, vec3 b, vec3 c) {
		int n = (int)vertices.size();
		vertices.push_back(a);
		vertices.push_back(b);
		vertices.push_back(c);
		faces.push_back(ivec3(n, n + 1, n + 2));
	};
	vec3 translate = vec3(0.0);
	for (std::vector<ivec3> triangulation : triangulations) {
		for (ivec3 e : triangulation) {
			vec3 v0 = 0.5f * vec3(VERTICE_LIST[EDGE_LIST[e.x].x] + VERTICE_LIST[EDGE_LIST[e.x].y]);
			vec3 v1 = 0.5f * vec3(VERTICE_LIST[EDGE_LIST[e.y].x] + VERTICE_LIST[EDGE_LIST[e.y].y]);
			vec3 v2 = 0.5f * vec3(VERTICE_LIST[EDGE_LIST[e.z].x] + VERTICE_LIST[EDGE_LIST[e.z].y]);
			addTriangle(translate + v0, translate + v1, translate + v2);
		}
		translate += vec3(1, 0, 0);
	}
	WritePLY(filepath, &vertices[0], (int)vertices.size(), &faces[0], (int)faces.size());
}

int main(int argc, char* argv[]) {
	initCubeCases();

	int index = 15;
	std::vector<std::vector<ivec3>> triangulations = completeTriangulation(index, cubes[index].faceIndex);

	printf("%d\n", triangulations.size());
	for (std::vector<ivec3> triangulation : triangulations) {
		for (ivec3 e : triangulation) printf("%d %d %d  ", e.x, e.y, e.z);
		printf("\n");
	}

	visualizeTriangulations(triangulations, "D:\\.ply");

	return 0;
}
