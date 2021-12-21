// Attempt to find an adaptive marching cube look up table that
// generates a mesh to best represent the isosurface.

// References:
// http://paulbourke.net/geometry/polygonise/ (indices of vertices/edges)
// https://www.cs.upc.edu/~virtual/SGI/docs/1.%20Theory/Unit%2010.%20Volume%20models.%20Marching%20Cubes/Marching%20Cubes.pdf


#include <vector>
#include <initializer_list>
#include <map>
#include <assert.h>

#include "../octatree/trigs2mesh.h"
#include "../octatree/ply_writer.h"

#include <chrono>

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
int calcIndex(const std::initializer_list<float> c) {
	return calcIndex(&*c.begin());
}

void toIndex(const int index, float v[8]) {
	for (int i = 0; i < 8; i++)
		v[i] = (index >> i) & 1 ? -1.0f : 1.0f;
}

int calcFaceIndex(const float v[4]) {
	if (isnan(v[0] + v[1] + v[2] + v[3]))
		return 0;
	return int(v[0] <= 0) |
		(int(v[1] <= 0) << 1) |
		(int(v[2] <= 0) << 2) |
		(int(v[3] <= 0) << 3);
}


bool isIdenticalEdge(ivec2 a, ivec2 b) {
	return a == b || a == ivec2(b.y, b.x);
}
bool isIdenticalTriangle(ivec3 a, ivec3 b) {
	std::sort(&a.x, &a.x + 3);
	std::sort(&b.x, &b.x + 3);
	return a == b;
}
bool isIdenticalTriangulation(std::vector<ivec3> tg1, std::vector<ivec3> tg2) {
	if (tg1.size() != tg2.size()) return false;
	for (int i = 0; i < (int)tg1.size(); i++) {
		bool found = false;
		for (int j = 0; j < (int)tg2.size(); j++) {
			if (isIdenticalTriangle(tg1[i], tg2[j])) { found = true; break; }
		}
		if (!found) return false;
	}
	return true;
}
bool isIdenticalFace(ivec4 a, ivec4 b) {  // not oriented
	std::sort(&a.x, &a.x + 4);
	std::sort(&b.x, &b.x + 4);
	return a == b;
}

/* cube cases */

struct TriangulationCase {
	std::vector<int> faceIndex = std::vector<int>({ 0, 0, 0, 0, 0, 0 });  // 0-17
	std::vector<std::vector<ivec3>> triangulations;
};

struct CubeCase {
	float v[8];
	int index = -1;
	int case_id = -1;  // 0-14
	bool isAmbiguous = false;
	int faceIndex[6] = { -1, -1, -1, -1, -1, -1 };  // 0-15
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
	std::vector<TriangulationCase> triangulationCases;
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
		//printf("%d(%d) ", i, cubes[i].case_id);
		cube_ids.push_back(i);
		cubes[i].calcFaceIndices();
		if (cubes[i].hasAmbiguousFace) assert(cubes[i].isAmbiguous);
		else if (cubes[i].isAmbiguous) assert(cubes[i].case_id == 4);
	}
	//printf("\n");  // 0(0) 1(1) 9(2) 14(5) 15(8) 30(12) 39(9) 46(14) 65(4) 71(11) 73(6) 85(10) 88(7) 129(3) 165(13)

	// generated by: bruteforce enumation + manually remove bad cases
	{
		int index = 1;  // cube case 1
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,1,1,0,0,0 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,3,0), }),
			}) });
	}
	{
		int index = 9;  // cube case 2
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,3,9,1,0,0 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,2,0), ivec3(11,2,8), }),
				std::vector<ivec3>({ ivec3(8,11,0), ivec3(0,11,2), }),
			}) });
	}
	{
		int index = 14;  // cube case 5
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 8,2,14,9,3,0 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(0,3,9), ivec3(11,9,3), ivec3(10,9,11), }),
				std::vector<ivec3>({ ivec3(0,3,9), ivec3(9,3,10), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(0,10,9), ivec3(3,10,0), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(0,10,9), ivec3(10,0,11), ivec3(11,0,3), }),
				std::vector<ivec3>({ ivec3(0,11,9), ivec3(11,0,3), ivec3(10,9,11), }),
			}) });
	}
	{
		int index = 15;  // cube case 8
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 9,3,15,9,3,0 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,10,9), ivec3(11,10,8), }),
				std::vector<ivec3>({ ivec3(8,11,9), ivec3(10,9,11), }),
			}) });
	}
	{
		int index = 30;  // cube case 12
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 10,10,14,9,3,1 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(0,3,9), ivec3(11,9,3), ivec3(10,9,11), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(0,3,9), ivec3(9,3,10), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(0,10,9), ivec3(3,10,0), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(0,10,9), ivec3(10,0,11), ivec3(11,0,3), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(0,11,9), ivec3(11,0,3), ivec3(10,9,11), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 17,10,14,9,3,1 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,7,9), ivec3(9,7,10), ivec3(8,10,7), ivec3(0,10,8), ivec3(3,10,0), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(4,7,9), ivec3(9,7,10), ivec3(8,10,7), ivec3(0,10,8), ivec3(10,0,11), ivec3(11,0,3), }),
				std::vector<ivec3>({ ivec3(4,10,9), ivec3(7,10,4), ivec3(8,10,7), ivec3(0,10,8), ivec3(3,10,0), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(4,10,9), ivec3(7,10,4), ivec3(8,10,7), ivec3(0,10,8), ivec3(10,0,11), ivec3(11,0,3), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 10,17,14,9,3,1 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(8,10,3), ivec3(3,10,0), ivec3(0,10,9), ivec3(10,4,11), ivec3(11,4,7), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(8,10,3), ivec3(3,10,0), ivec3(0,10,9), ivec3(7,10,4), ivec3(11,10,7), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(8,10,3), ivec3(10,4,11), ivec3(11,4,7), ivec3(9,3,10), ivec3(0,3,9), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(8,10,3), ivec3(9,3,10), ivec3(0,3,9), ivec3(7,10,4), ivec3(11,10,7), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 17,17,14,9,3,1 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,7,9), ivec3(0,3,8), ivec3(11,9,7), ivec3(10,9,11), }),
				std::vector<ivec3>({ ivec3(4,7,9), ivec3(0,3,8), ivec3(9,7,10), ivec3(11,10,7), }),
				std::vector<ivec3>({ ivec3(4,10,9), ivec3(0,3,8), ivec3(10,4,11), ivec3(11,4,7), }),
				std::vector<ivec3>({ ivec3(4,10,9), ivec3(0,3,8), ivec3(7,10,4), ivec3(11,10,7), }),
				std::vector<ivec3>({ ivec3(4,11,9), ivec3(0,3,8), ivec3(11,4,7), ivec3(10,9,11), }),
			}) });
	}
	{
		int index = 39;  // cube case 9
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 13,1,7,8,11,2 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,2,4), ivec3(3,2,8), ivec3(10,4,2), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(8,2,4), ivec3(3,2,8), ivec3(4,2,5), ivec3(10,5,2), }),
				std::vector<ivec3>({ ivec3(8,3,4), ivec3(2,4,3), ivec3(10,4,2), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(8,3,4), ivec3(2,4,3), ivec3(4,2,5), ivec3(10,5,2), }),
				std::vector<ivec3>({ ivec3(8,3,4), ivec3(4,3,5), ivec3(2,5,3), ivec3(10,5,2), }),
				std::vector<ivec3>({ ivec3(8,3,4), ivec3(4,3,5), ivec3(5,3,10), ivec3(2,10,3), }),
				std::vector<ivec3>({ ivec3(8,5,4), ivec3(3,5,8), ivec3(2,5,3), ivec3(10,5,2), }),
				std::vector<ivec3>({ ivec3(8,5,4), ivec3(3,5,8), ivec3(5,3,10), ivec3(2,10,3), }),
				std::vector<ivec3>({ ivec3(8,5,4), ivec3(5,8,10), ivec3(3,10,8), ivec3(2,10,3), }),
				std::vector<ivec3>({ ivec3(8,5,4), ivec3(5,8,10), ivec3(10,8,2), ivec3(3,2,8), }),
				std::vector<ivec3>({ ivec3(8,10,4), ivec3(3,10,8), ivec3(2,10,3), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(8,10,4), ivec3(10,8,2), ivec3(3,2,8), ivec3(5,4,10), }),
			}) });
	}
	{
		int index = 46;  // cube case 14
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 12,2,14,9,11,2 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(0,3,4), ivec3(11,4,3), ivec3(10,4,11), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(0,3,4), ivec3(11,4,3), ivec3(4,11,5), ivec3(10,5,11), }),
				std::vector<ivec3>({ ivec3(0,3,4), ivec3(4,3,5), ivec3(11,5,3), ivec3(10,5,11), }),
				std::vector<ivec3>({ ivec3(0,3,4), ivec3(4,3,5), ivec3(5,3,10), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(0,5,4), ivec3(3,5,0), ivec3(11,5,3), ivec3(10,5,11), }),
				std::vector<ivec3>({ ivec3(0,5,4), ivec3(3,5,0), ivec3(5,3,10), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(0,5,4), ivec3(5,0,10), ivec3(3,10,0), ivec3(11,10,3), }),
				std::vector<ivec3>({ ivec3(0,5,4), ivec3(5,0,10), ivec3(10,0,11), ivec3(11,0,3), }),
				std::vector<ivec3>({ ivec3(0,10,4), ivec3(3,10,0), ivec3(11,10,3), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(0,10,4), ivec3(10,0,11), ivec3(11,0,3), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(0,11,4), ivec3(11,0,3), ivec3(10,4,11), ivec3(5,4,10), }),
				std::vector<ivec3>({ ivec3(0,11,4), ivec3(11,0,3), ivec3(4,11,5), ivec3(10,5,11), }),
			}) });
	}
	{
		int index = 65;  // cube case 4
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,1,1,4,4,4 }),
			std::vector<std::vector<ivec3>>({  // has a lot of valid triangulations, choose two practical
				std::vector<ivec3>({ ivec3(8,3,0), ivec3(6,5,10), }),
				std::vector<ivec3>({ ivec3(8,5,0), ivec3(10,0,5), ivec3(0,10,3), ivec3(6,3,10), ivec3(3,6,8), ivec3(5,8,6), }),
			}) });
	}
	{
		int index = 71;  // cube case 11
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 9,1,7,12,7,4 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,2,9), ivec3(3,2,8), ivec3(6,9,2), ivec3(9,6,5), }),
				std::vector<ivec3>({ ivec3(8,2,9), ivec3(3,2,8), ivec3(9,2,5), ivec3(6,5,2), }),
				std::vector<ivec3>({ ivec3(8,3,9), ivec3(2,9,3), ivec3(6,9,2), ivec3(9,6,5), }),
				std::vector<ivec3>({ ivec3(8,3,9), ivec3(2,9,3), ivec3(9,2,5), ivec3(6,5,2), }),
				std::vector<ivec3>({ ivec3(8,3,9), ivec3(9,3,5), ivec3(2,5,3), ivec3(6,5,2), }),
				std::vector<ivec3>({ ivec3(8,3,9), ivec3(9,3,5), ivec3(5,3,6), ivec3(2,6,3), }),
				std::vector<ivec3>({ ivec3(8,5,9), ivec3(3,5,8), ivec3(2,5,3), ivec3(6,5,2), }),
				std::vector<ivec3>({ ivec3(8,5,9), ivec3(3,5,8), ivec3(5,3,6), ivec3(2,6,3), }),
				std::vector<ivec3>({ ivec3(8,5,9), ivec3(5,8,6), ivec3(3,6,8), ivec3(2,6,3), }),
				std::vector<ivec3>({ ivec3(8,5,9), ivec3(5,8,6), ivec3(6,8,2), ivec3(3,2,8), }),
				std::vector<ivec3>({ ivec3(8,6,9), ivec3(3,6,8), ivec3(2,6,3), ivec3(9,6,5), }),
				std::vector<ivec3>({ ivec3(8,6,9), ivec3(6,8,2), ivec3(3,2,8), ivec3(9,6,5), }),
			}) });
	}
	{
		int index = 73;  // cube case 6
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,3,9,5,4,4 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,2,0), ivec3(11,2,8), ivec3(6,5,10), }),
				std::vector<ivec3>({ ivec3(8,11,0), ivec3(0,11,2), ivec3(6,5,10), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,3,9,16,4,4 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,2,0), ivec3(2,8,10), ivec3(10,8,5), ivec3(11,5,8), ivec3(6,5,11), }),
				std::vector<ivec3>({ ivec3(8,2,0), ivec3(2,8,10), ivec3(10,8,5), ivec3(5,8,6), ivec3(11,6,8), }),
				std::vector<ivec3>({ ivec3(8,5,0), ivec3(11,5,8), ivec3(0,5,2), ivec3(6,5,11), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(8,5,0), ivec3(11,5,8), ivec3(6,5,11), ivec3(10,0,5), ivec3(0,10,2), }),
				std::vector<ivec3>({ ivec3(8,5,0), ivec3(0,5,2), ivec3(2,5,10), ivec3(5,8,6), ivec3(11,6,8), }),
				std::vector<ivec3>({ ivec3(8,5,0), ivec3(10,0,5), ivec3(0,10,2), ivec3(5,8,6), ivec3(11,6,8), }),
				std::vector<ivec3>({ ivec3(8,6,0), ivec3(11,6,8), ivec3(5,0,6), ivec3(0,5,2), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(8,6,0), ivec3(11,6,8), ivec3(5,0,6), ivec3(10,0,5), ivec3(0,10,2), }),
				std::vector<ivec3>({ ivec3(8,10,0), ivec3(0,10,2), ivec3(10,8,5), ivec3(11,5,8), ivec3(6,5,11), }),
				std::vector<ivec3>({ ivec3(8,10,0), ivec3(0,10,2), ivec3(10,8,5), ivec3(5,8,6), ivec3(11,6,8), }),
				std::vector<ivec3>({ ivec3(8,11,0), ivec3(6,0,11), ivec3(5,0,6), ivec3(0,5,2), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(8,11,0), ivec3(6,0,11), ivec3(5,0,6), ivec3(10,0,5), ivec3(0,10,2), }),
			}) });
	}
	{
		int index = 85;  // cube case 10
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 3,9,5,12,6,5 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,3,0), ivec3(3,4,7), ivec3(2,5,1), ivec3(6,5,2), }),
				std::vector<ivec3>({ ivec3(4,3,0), ivec3(3,4,7), ivec3(2,6,1), ivec3(1,6,5), }),
				std::vector<ivec3>({ ivec3(4,7,0), ivec3(3,0,7), ivec3(2,5,1), ivec3(6,5,2), }),
				std::vector<ivec3>({ ivec3(4,7,0), ivec3(3,0,7), ivec3(2,6,1), ivec3(1,6,5), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 3,9,16,12,6,16 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,1,0), ivec3(3,2,7), ivec3(6,7,2), ivec3(1,4,5), }),
				std::vector<ivec3>({ ivec3(4,1,0), ivec3(3,6,7), ivec3(2,6,3), ivec3(1,4,5), }),
				std::vector<ivec3>({ ivec3(4,5,0), ivec3(3,2,7), ivec3(0,5,1), ivec3(6,7,2), }),
				std::vector<ivec3>({ ivec3(4,5,0), ivec3(3,6,7), ivec3(0,5,1), ivec3(2,6,3), }),
			}) });
		// troll cases: { 3,9,16,12,6,5 }, { 3,9,5,12,6,16 }
	}
	{
		int index = 88;  // cube case 7
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 2,10,8,5,4,5 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(11,2,3), ivec3(6,5,10), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 2,17,8,5,4,5 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,2,8), ivec3(8,2,3), ivec3(2,4,11), ivec3(11,4,7), ivec3(6,5,10), }),
				std::vector<ivec3>({ ivec3(4,2,8), ivec3(8,2,3), ivec3(6,5,10), ivec3(7,2,4), ivec3(11,2,7), }),
				std::vector<ivec3>({ ivec3(4,3,8), ivec3(3,4,2), ivec3(2,4,11), ivec3(11,4,7), ivec3(6,5,10), }),
				std::vector<ivec3>({ ivec3(4,3,8), ivec3(3,4,2), ivec3(6,5,10), ivec3(7,2,4), ivec3(11,2,7), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 2,10,8,16,4,5 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(11,5,3), ivec3(3,5,2), ivec3(6,5,11), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(11,5,3), ivec3(6,5,11), ivec3(10,3,5), ivec3(3,10,2), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(11,6,3), ivec3(5,3,6), ivec3(3,5,2), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(4,7,8), ivec3(11,6,3), ivec3(5,3,6), ivec3(10,3,5), ivec3(3,10,2), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 2,10,8,5,4,16 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,5,8), ivec3(11,2,3), ivec3(10,8,5), ivec3(8,10,7), ivec3(6,7,10), }),
				std::vector<ivec3>({ ivec3(4,5,8), ivec3(11,2,3), ivec3(10,8,5), ivec3(6,8,10), ivec3(8,6,7), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(11,2,3), ivec3(8,10,7), ivec3(6,7,10), ivec3(10,4,5), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(11,2,3), ivec3(6,8,10), ivec3(8,6,7), ivec3(10,4,5), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 2,17,8,16,4,16 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,2,8), ivec3(11,6,7), ivec3(8,2,3), ivec3(2,4,10), ivec3(10,4,5), }),
				std::vector<ivec3>({ ivec3(4,2,8), ivec3(11,6,7), ivec3(8,2,3), ivec3(5,2,4), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(4,3,8), ivec3(11,6,7), ivec3(3,4,2), ivec3(2,4,10), ivec3(10,4,5), }),
				std::vector<ivec3>({ ivec3(4,3,8), ivec3(11,6,7), ivec3(3,4,2), ivec3(5,2,4), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(4,3,8), ivec3(11,6,7), ivec3(5,3,4), ivec3(3,5,2), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(4,3,8), ivec3(11,6,7), ivec3(5,3,4), ivec3(10,3,5), ivec3(3,10,2), }),
				std::vector<ivec3>({ ivec3(4,5,8), ivec3(11,6,7), ivec3(8,5,3), ivec3(3,5,2), ivec3(2,5,10), }),
				std::vector<ivec3>({ ivec3(4,5,8), ivec3(11,6,7), ivec3(8,5,3), ivec3(10,3,5), ivec3(3,10,2), }),
				std::vector<ivec3>({ ivec3(4,5,8), ivec3(11,6,7), ivec3(10,8,5), ivec3(8,10,3), ivec3(3,10,2), }),
				std::vector<ivec3>({ ivec3(4,5,8), ivec3(11,6,7), ivec3(10,8,5), ivec3(2,8,10), ivec3(8,2,3), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(11,6,7), ivec3(8,10,3), ivec3(3,10,2), ivec3(10,4,5), }),
				std::vector<ivec3>({ ivec3(4,10,8), ivec3(11,6,7), ivec3(2,8,10), ivec3(8,2,3), ivec3(10,4,5), }),
			}) });
		// troll cases: { 2,17,8,16,4,5 }, { 2,17,8,5,4,16 }, { 2,10,8,16,4,16 }
	}
	{
		int index = 129;  // cube case 4
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,5,1,2,0,8 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,3,0), ivec3(7,6,11), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 1,16,1,2,0,8 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,6,0), ivec3(7,6,8), ivec3(0,6,3), ivec3(3,6,11), }),
				std::vector<ivec3>({ ivec3(8,6,0), ivec3(7,6,8), ivec3(11,0,6), ivec3(3,0,11), }),
				std::vector<ivec3>({ ivec3(8,7,0), ivec3(6,0,7), ivec3(0,6,3), ivec3(3,6,11), }),
				std::vector<ivec3>({ ivec3(8,7,0), ivec3(6,0,7), ivec3(11,0,6), ivec3(3,0,11), }),
			}) });
	}
	{
		int index = 165;  // cube case 13 (the troll case)
		cubes[index].triangulationCases.push_back(TriangulationCase{
			std::vector<int>({ 5,5,5,10,10,10 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(8,3,0), ivec3(9,5,4), ivec3(7,6,11), ivec3(2,10,1), }),
			}) });
		cubes[index].triangulationCases.push_back(TriangulationCase{  // manually added
			std::vector<int>({ 16,16,16,17,17,17 }),
			std::vector<std::vector<ivec3>>({
				std::vector<ivec3>({ ivec3(4,8,7), ivec3(2,11,3), ivec3(5,6,10), ivec3(0,9,1), }),
				}) });
		// only two solutions here, subjects to change
	}
	for (int ci : cube_ids) {
		if (ci) assert(!cubes[ci].triangulationCases.empty());
	}

	// spread to all cube cases (BFS)
	for (int i_ = 0; i_ < (int)cube_ids.size(); i_++) {
		int i = cube_ids[i_];
		int index;

		// bit reverse
		index = 255 - cubes[i].index;
		if (cubes[index].index == -1) {
			CubeCase cbc = CubeCase(index);
			cbc.calcFaceIndices();
			cbc.case_id = cubes[i].case_id;
			cbc.isAmbiguous = cubes[i].isAmbiguous;
			for (TriangulationCase tgl : cubes[i].triangulationCases) {
				TriangulationCase tgl1;
				for (int _ = 0; _ < 6; _++) {
					int fi = tgl.faceIndex[_];
					assert(fi >= 0 && fi < 18);
					fi = fi == 5 ? 16 : fi == 10 ? 17 : fi == 16 ? 5 : fi == 17 ? 10 : 15 - fi;
					tgl1.faceIndex[_] = fi;
					if (fi == 5 || fi == 10 || fi == 16 || fi == 17) {
						assert(cbc.isAmbiguousFace[_]);
						assert(cbc.faceIndex[_] == 5 || cbc.faceIndex[_] == 10 || cbc.faceIndex[_] == 16 || cbc.faceIndex[_] == 17);
					}
					else assert(cbc.faceIndex[_] == fi);
				}
				tgl1.triangulations = tgl.triangulations;
				for (int k = 0; k < (int)tgl1.triangulations.size(); k++) {
					for (int j = 0; j < (int)tgl1.triangulations[k].size(); j++)
						tgl1.triangulations[k][j] = ivec3(vec3(tgl1.triangulations[k][j]).zyx());
				}
				cbc.triangulationCases.push_back(tgl1);
			}
			cubes[index] = cbc;
			cube_ids.push_back(index);
		}

		// reflection / rotation
		auto applyTransform = [&](const int EDGE_MAP[12], const int VERTEX_MAP[8], bool hasReflection) {
			// calculate transformed index
			int index = 0;
			for (int j = 0; j < 8; j++) {
				int k = VERTEX_MAP[j];
				index |= ((cubes[i].index >> j) & 1) << k;
			}
			if (cubes[index].index == -1) {
				// create transformed CubeCase
				CubeCase cbc = CubeCase(index);
				cbc.calcFaceIndices();
				cbc.case_id = cubes[i].case_id;
				cbc.isAmbiguous = cubes[i].isAmbiguous;
				for (TriangulationCase tgl : cubes[i].triangulationCases) {
					// apply transform to triangles
					TriangulationCase tgl1;
					tgl1.triangulations = tgl.triangulations;
					for (int k = 0; k < (int)tgl.triangulations.size(); k++) {
						for (int j = 0; j < (int)tgl.triangulations[k].size(); j++) {
							ivec3 t;
							t.x = EDGE_MAP[tgl.triangulations[k][j].x];
							t.y = EDGE_MAP[tgl.triangulations[k][j].y];
							t.z = EDGE_MAP[tgl.triangulations[k][j].z];
							if (hasReflection) std::swap(t.x, t.z);
							tgl1.triangulations[k][j] = t;
						}
					}
					// apply transform to faceIndex
					for (int fi = 0; fi < 6; fi++) {
						int edges[4];
						for (int _ = 0; _ < 4; _++) edges[_] = EDGE_MAP[FACE_EDGE_LIST[fi][_]];
						int new_fi = -1;
						for (int _ = 0; _ < 16; _++) {
							if (isIdenticalFace(*(ivec4*)edges, *(ivec4*)FACE_EDGE_LIST[_])) {
								new_fi = _;
								break;
							}
						}
						assert(new_fi != -1);
						const int* fel = FACE_EDGE_LIST[new_fi];
						int fI = tgl.faceIndex[fi];
						std::vector<ivec2> edge_occur;
						for (int k = 0; k < (int)tgl1.triangulations.size(); k++) {
							for (int j = 0; j < (int)tgl1.triangulations[k].size(); j++) {
								ivec3 t = tgl1.triangulations[k][j];
								for (int _ = 0; _ < 3; _++) {
									ivec2 e = ivec2((&t.x)[_], (&t.x)[(_ + 1) % 3]);
									ivec2 found = ivec2(
										std::find(fel, fel + 4, e.x) - fel,
										std::find(fel, fel + 4, e.y) - fel);
									if (new_fi < 3) std::swap(found.x, found.y);
									if (found.x != 4 && found.y != 4 &&
										std::find(edge_occur.begin(), edge_occur.end(), found) == edge_occur.end())
										edge_occur.push_back(found);
								}
							}
						}
						assert((int)edge_occur.size() <= 2);
						if (edge_occur.empty()) {
							assert(fI == 0 || fI == 15);
						}
						else if (edge_occur.size() == 1) {
							assert(fI < 15 && fI != 5 && fI != 10);
							fI = -1;
							for (int _ : std::vector<int>({ 1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14 })) {
								if (edge_occur[0] == *(ivec2*)LUT_FACE[_]) {
									fI = _; break;
								}
							}
							assert(fI != -1);
						}
						else {
							assert(fI == 5 || fI == 10 || fI == 16 || fI == 17);
							int fI1 = -1;
							for (int _ : std::vector<int>({ 5, 10, 16, 17 })) {
								if ((edge_occur[0] == ((ivec2*)LUT_FACE[_])[0] && edge_occur[1] == ((ivec2*)LUT_FACE[_])[1]) ||
									(edge_occur[1] == ((ivec2*)LUT_FACE[_])[0] && edge_occur[0] == ((ivec2*)LUT_FACE[_])[1])) {
									fI1 = _; break;
								}
							}
							assert(fI1 != -1);
							int fI0 = cbc.faceIndex[new_fi];
							assert((fI0 == 5 || fI0 == 16) == (fI1 == 5 || fI1 == 16));
							assert((fI0 == 10 || fI0 == 17) == (fI1 == 10 || fI1 == 17));
							fI = fI1;
						}
						tgl1.faceIndex[new_fi] = fI;
					}
					// push transform
					cbc.triangulationCases.push_back(tgl1);
				}
				cubes[index] = cbc;
				cube_ids.push_back(index);
			}
		};

		// reflection about z-axis
		const int REFLECT_Z_EDGE_MAP[12] = {
			4, 5, 6, 7, 0, 1, 2, 3, 8, 9, 10, 11
		};
		applyTransform(REFLECT_Z_EDGE_MAP, REFLECT_Z_EDGE_MAP, true);

		// rotation about z-axis
		const int ROTATE_Z_EDGE_MAP[12] = {
			3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10
		};
		const int ROTATE_Z_VERTEX_MAP[8] = {
			3, 0, 1, 2, 7, 4, 5, 6
		};
		applyTransform(ROTATE_Z_EDGE_MAP, ROTATE_Z_VERTEX_MAP, false);

		// rotation about x-axis clockwise
		const int ROTATE_X_EDGE_MAP[12] = {
			8, 3, 11, 7, 9, 1, 10, 5, 4, 0, 2, 6
		};
		const int ROTATE_X_VERTEX_MAP[8] = {
			4, 0, 3, 7, 5, 1, 2, 6
		};
		applyTransform(ROTATE_X_EDGE_MAP, ROTATE_X_VERTEX_MAP, false);

	}

	int complete_count = 0;
	for (int i = 0; i < 256; i++) {
		//printf("%d%c", cubes[i].triangulationCases.size(), i % 16 == 15 ? '\n' : ' ');
		complete_count += cubes[i].index != -1;
	}
	//printf("%d cubes complete\n", complete_count);
}

// export lookup table to C++ code
void exportTable(const char* filepath) {
	FILE* fp = fopen(filepath, "w");
	fprintf(fp, "/*const*/ std::map<int, std::vector<std::vector<ivec3>>> DISAMBIGUATION_LUT[256] = {\n");

	for (int i = 0; i < 256; i++) {
		fprintf(fp, "/*%d*/ std::map<int, std::vector<std::vector<ivec3>>>({", i);
		if (!cubes[i].triangulationCases.empty()) {
			fprintf(fp, "\n");
			bool has_zero = false;
			for (TriangulationCase tc : cubes[i].triangulationCases) {
				int amb = 0;
				for (int b = 0; b < 6; b++) {
					int id = tc.faceIndex[b];
					amb = (amb << 1) | int(id >= 16);
				}
				has_zero |= amb == 0;
				fprintf(fp, "{ %d, std::vector<std::vector<ivec3>>({ ", amb);
				assert(!tc.triangulations.empty());
				for (std::vector<ivec3> triangulation : tc.triangulations) {
					fprintf(fp, "std::vector<ivec3>({ ");
					for (ivec3 t : triangulation) fprintf(fp, "ivec3(%d,%d,%d),", t.x, t.y, t.z);
					fprintf(fp, " }),");
				}
				fprintf(fp, "}) },\n");
			}
			assert(has_zero);
		}
		fprintf(fp, "}),\n");
	}

	fprintf(fp, "};");
	fclose(fp);
}


/* brute-force iterate through all possible triangulations */

// check if it's a valid triangulation
bool isValidTriangulation(
	const std::vector<int> &cube_edges,
	const std::vector<ivec2> &required_trig_edges,
	std::vector<ivec3> trigs  // ccw
) {
	std::map<uint64_t, int> edge_count;
	for (ivec3 t : trigs) {
		uint64_t h = ((uint64_t)min(t.x, t.y) << 32ull) | (uint64_t)max(t.x, t.y);
		if (edge_count.find(h) != edge_count.end()) edge_count[h] += 1;
		else edge_count[h] = 1;
		h = ((uint64_t)min(t.x, t.z) << 32ull) | (uint64_t)max(t.x, t.z);
		if (edge_count.find(h) != edge_count.end()) edge_count[h] += 1;
		else edge_count[h] = 1;
		h = ((uint64_t)min(t.y, t.z) << 32ull) | (uint64_t)max(t.y, t.z);
		if (edge_count.find(h) != edge_count.end()) edge_count[h] += 1;
		else edge_count[h] = 1;
	}
	// each surface edge must be connected to one triangle edge
	for (ivec2 edge : required_trig_edges) {
		uint64_t h = ((uint64_t)min(edge.x, edge.y) << 32) | (uint64_t)max(edge.x, edge.y);
		if (edge_count.find(h) == edge_count.end()) return false;
		if (edge_count[h] != 1) return false;
		edge_count[h] = 0;
	}
	// each non-surface edge must be connected to exactly 2 triangles
	for (std::pair<uint64_t, int> ec : edge_count) {
		if (ec.second != 0 && ec.second != 2) return false;
	}
	// triangles connected to a vertex must be connected by common edges
	int tn = (int)trigs.size();
	DisjointSet dsj(tn);
	for (int i = 0; i < tn; i++) {
		for (int j = 0; j < i; j++) {
			ivec3 t1 = trigs[i], t2 = trigs[j];
			int equ_count = int(t1.x == t2.x) + int(t1.x == t2.y) + int(t1.x == t2.z)
				+ int(t1.y == t2.x) + int(t1.y == t2.y) + int(t1.y == t2.z)
				+ int(t1.z == t2.x) + int(t1.z == t2.y) + int(t1.z == t2.z);
			if (equ_count == 2) dsj.unionSet(i, j);
		}
	}
	std::vector<int> vertice[12];
	for (int i = 0; i < tn; i++) {
		vertice[trigs[i].x].push_back(i);
		vertice[trigs[i].y].push_back(i);
		vertice[trigs[i].z].push_back(i);
	}
	for (int vi = 0; vi < 12; vi++) {
		int vn = (int)vertice[vi].size();
		if (vn == 0) continue;
		for (int i = 0; i < vn; i++) for (int j = 0; j < i; j++) {
			if (dsj.findRepresentative(vertice[vi][i]) !=
				dsj.findRepresentative(vertice[vi][j])) return false;
		}
	}
	return true;
}

// check if it's impossible to form a valid triangulation by adding more trigs
bool isImpossibleTriangulation(
	const std::vector<int> &cube_edges,
	const std::vector<ivec2> &required_trig_edges,
	std::vector<ivec3> trigs  // ccw
) {
	std::map<uint64_t, int> edge_count;
	for (ivec3 t : trigs) {
		uint64_t h = ((uint64_t)min(t.x, t.y) << 32ull) | (uint64_t)max(t.x, t.y);
		if (edge_count.find(h) != edge_count.end()) edge_count[h] += 1;
		else edge_count[h] = 1;
		h = ((uint64_t)min(t.x, t.z) << 32ull) | (uint64_t)max(t.x, t.z);
		if (edge_count.find(h) != edge_count.end()) edge_count[h] += 1;
		else edge_count[h] = 1;
		h = ((uint64_t)min(t.y, t.z) << 32ull) | (uint64_t)max(t.y, t.z);
		if (edge_count.find(h) != edge_count.end()) edge_count[h] += 1;
		else edge_count[h] = 1;
	}
	// if there are more than 1 triangles connected to a surface edge
	for (ivec2 edge : required_trig_edges) {
		uint64_t h = ((uint64_t)min(edge.x, edge.y) << 32) | (uint64_t)max(edge.x, edge.y);
		if (edge_count[h] > 1) return true;
		edge_count[h] = 0;
	}
	// if there are more than 2 triangles connected to an edge
	for (std::pair<uint64_t, int> ec : edge_count) {
		if (ec.second > 2) return true;
	}
	// if there is a surface edge that is not a required edge
	for (std::pair<uint64_t, int> ec : edge_count) {
		if (ec.second == 0) continue;
		ivec2 e = *(ivec2*)&ec.first;
		for (int i = 0; i < 6; i++) {
			auto fi = FACE_EDGE_LIST[i];
			ivec2 on_face = ivec2(0);
			for (int j = 0; j < 4; j++) {
				if (fi[j] == e.x) on_face.x = 1;
				if (fi[j] == e.y) on_face.y = 1;
			}
			if (on_face == ivec2(1)) return true;
		}
	}
	// triangles connected to a vertex must be connected by common edges
	// this should not affect the result because it bruteforces through all triangle-adding orders
	int tn = (int)trigs.size();
	DisjointSet dsj(tn);
	for (int i = 0; i < tn; i++) {
		for (int j = 0; j < i; j++) {
			ivec3 t1 = trigs[i], t2 = trigs[j];
			int equ_count = int(t1.x == t2.x) + int(t1.x == t2.y) + int(t1.x == t2.z)
				+ int(t1.y == t2.x) + int(t1.y == t2.y) + int(t1.y == t2.z)
				+ int(t1.z == t2.x) + int(t1.z == t2.y) + int(t1.z == t2.z);
			if (equ_count == 2) dsj.unionSet(i, j);
		}
	}
	std::vector<int> vertice[12];
	for (int i = 0; i < tn; i++) {
		vertice[trigs[i].x].push_back(i);
		vertice[trigs[i].y].push_back(i);
		vertice[trigs[i].z].push_back(i);
	}
	for (int vi = 0; vi < 12; vi++) {
		int vn = (int)vertice[vi].size();
		if (vn == 0) continue;
		for (int i = 0; i < vn; i++) for (int j = 0; j < i; j++) {
			if (dsj.findRepresentative(vertice[vi][i]) !=
				dsj.findRepresentative(vertice[vi][j])) return true;
		}
	}
	return false;
}

int _completeTriangulation_recurse_callCount = 0;
void _completeTriangulation_recurse(
	std::vector<std::vector<ivec3>> &valid_triangulations,
	const std::vector<int> &cube_edges,
	const std::vector<ivec2> &required_trig_edges,
	std::vector<ivec3> &trigs
) {
	// factorial time complexity
	const int MAX_CALL = 1 << 21;
	if (++_completeTriangulation_recurse_callCount >= MAX_CALL) {
		if (_completeTriangulation_recurse_callCount == MAX_CALL)
			fprintf(stderr, "/* Aborted */\n");
		return;
	}
	if (isImpossibleTriangulation(cube_edges, required_trig_edges, trigs))
		return;
	// attempt to construct a new triangle
	for (ivec2 ab : required_trig_edges) for (int c : cube_edges) {
		if (c == ab.x || c == ab.y) continue;
		ivec3 t = ivec3(ab.x, c, ab.y);
		bool is_dup = false;
		for (ivec3 t1 : trigs)
			if (isIdenticalTriangle(t, t1)) { is_dup = true; break; }
		if (is_dup) continue;
		bool is_on_face = false;
		for (int i = 0; i < 6; i++) {
			auto fi = FACE_EDGE_LIST[i];
			ivec3 on_face = ivec3(0);
			for (int j = 0; j < 4; j++) {
				if (fi[j] == t.x) on_face.x = 1;
				if (fi[j] == t.y) on_face.y = 1;
				if (fi[j] == t.z) on_face.z = 1;
			}
			if (on_face == ivec3(1)) { is_on_face = true; break; }
		}
		if (is_on_face) continue;
		// recursion
		trigs.push_back(t);
		//for (ivec3 e : trigs) printf("%d %d %d  ", e.x, e.y, e.z); printf("\n");
		if (isValidTriangulation(cube_edges, required_trig_edges, trigs)) {
			bool is_new = true;
			for (int i = 0; i < (int)valid_triangulations.size(); i++) {
				if (isIdenticalTriangulation(trigs, valid_triangulations[i])) {
					is_new = false; break;
				}
			}
			if (is_new) {
				valid_triangulations.push_back(trigs);
				//printf("%d.", valid_triangulations.size());
			}
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
	_completeTriangulation_recurse_callCount = 0;
	_completeTriangulation_recurse(
		res,
		cube_edges, required_trig_edges,
		trigs
	);
	return res;
}

// designed to automatically bruteforce the tetrahedron case
void bruteforceCase13() {
	const int index = 165;
	const int N = 64;  // 2**6
	for (int bi = 0; bi < N; bi++) {
		int faceIndex[6];
		for (int i = 0; i < 6; i++) {
			faceIndex[i] = cubes[index].faceIndex[i];
			assert(faceIndex[i] == 5 || faceIndex[i] == 10);
			if (((bi >> i) & 1) == 1) {
				if (faceIndex[i] == 5) faceIndex[i] = 16;
				if (faceIndex[i] == 10) faceIndex[i] = 17;
			}
		}
		std::vector<std::vector<ivec3>> triangulations = completeTriangulation(index, faceIndex);
		if (!triangulations.empty()) {
			printf("\
cubes[index].triangulationCases.push_back(TriangulationCase{\n\
    std::vector<int>({ %d,%d,%d,%d,%d,%d }),\n",
				faceIndex[0], faceIndex[1], faceIndex[2], faceIndex[3], faceIndex[4], faceIndex[5]);
			printf("    std::vector<std::vector<ivec3>>({\n");
			for (std::vector<ivec3> triangulation : triangulations) {
				printf("        std::vector<ivec3>({ ");
				for (ivec3 e : triangulation) printf("ivec3(%d,%d,%d), ", e.x, e.y, e.z);
				printf("}),\n");
			}
			printf("    }) });\n");
		}
		else {
			printf("// troll case { %d,%d,%d,%d,%d,%d }\n",
				faceIndex[0], faceIndex[1], faceIndex[2], faceIndex[3], faceIndex[4], faceIndex[5]);
		}
	}
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
	if (vertices.empty() || faces.empty()) fprintf(stderr, "Empty triangulation\n");
	else WritePLY(filepath, &vertices[0], (int)vertices.size(), &faces[0], (int)faces.size());
}

int main(int argc, char* argv[]) {
	initCubeCases();

	exportTable("D:\\disambiguation_lut.h"); return 0;

	//bruteforceCase13(); return 0;

	int index = 165;
	int faceIndex[6];
	for (int i = 0; i < 6; i++) {
		faceIndex[i] = cubes[index].faceIndex[i];
		if ((i >= 0) && faceIndex[i] == 5) faceIndex[i] = 16;
		if ((i >= 0) && faceIndex[i] == 10) faceIndex[i] = 17;
		printf("%d%s", faceIndex[i], i == 5 ? "\n" : ",");
	}

	auto t0 = std::chrono::high_resolution_clock::now();
	std::vector<std::vector<ivec3>> triangulations = completeTriangulation(index, faceIndex);
	double time_elapsed = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
	printf("\n%lf ms\n", 1000.0 * time_elapsed);

	printf("%d\n", triangulations.size());
	printf("std::vector<std::vector<ivec3>>({\n");
	for (std::vector<ivec3> triangulation : triangulations) {
		printf("    std::vector<ivec3>({ ");
		for (ivec3 e : triangulation) printf("ivec3(%d,%d,%d), ", e.x, e.y, e.z);
		printf("}),\n");
	}
	printf("})");

	visualizeTriangulations(triangulations, "D:\\.ply");

	return 0;
}
