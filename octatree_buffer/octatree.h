// Implicit surface triangulation based on marching cube.
// Support exporting byte buffer that can be placed inside an OpenGL texture.


#include <vector>
#include <functional>
#include <algorithm>
#include <map>

#include <chrono>
#include <string>

#include "../octatree/trigs2mesh.h"


#define ScalarFieldTriangulator_octatree_BEGIN_ namespace ScalarFieldTriangulator_octatree {
#define ScalarFieldTriangulator_octatree_END_ }
#define ScalarFieldTriangulator_octatree_PRIVATE_BEGIN_ namespace __private__ {
#define ScalarFieldTriangulator_octatree_PRIVATE_END_ }


ScalarFieldTriangulator_octatree_BEGIN_

bool verbose = true;

ScalarFieldTriangulator_octatree_PRIVATE_BEGIN_

// Timing

std::map<std::string, std::chrono::high_resolution_clock::time_point> time_events;
void timeEventStart(std::string name) {
	int stack_size = time_events.size();
	auto t = std::chrono::high_resolution_clock::now();
	time_events[name] = t;
	if (verbose) {
		while (stack_size--) printf("  ");
		//printf("S: %s\n", &name[0]);
		printf("%s", &name[0]);
	}
}
void timeEventEnd(std::string name) {
	if (time_events.find(name) == time_events.end()) {
		fprintf(stderr, "No event \"%s\" found.\n", &name[0]);
	}
	else {
		auto t0 = time_events[name];
		time_events.erase(name);
		double dt = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
		int stack_size = time_events.size();
		if (verbose) {
			while (stack_size--) printf("  ");
			//printf("E: %s, %.1lf ms\n", &name[0], 1000.0*dt);
			printf(" - %.1lf ms\n", 1000.0*dt);
		}
	}
}


/* LOOKUP TABLES */

// list of vertice on a unit cube
const ivec3 VERTICE_LIST[8] = {
	ivec3(0,0,0), ivec3(0,1,0), ivec3(1,1,0), ivec3(1,0,0),
	ivec3(0,0,1), ivec3(0,1,1), ivec3(1,1,1), ivec3(1,0,1)
};
const int VERTICE_LIST_INV[2][2][2] = {
	{{0, 4}, {1, 5}}, {{3, 7}, {2, 6}}
};

// list of edges connecting two vertices on a unit cube
const ivec2 EDGE_LIST[12] = {
	ivec2(0,1), ivec2(1,2), ivec2(2,3), ivec2(3,0),
	ivec2(4,5), ivec2(5,6), ivec2(6,7), ivec2(7,4),
	ivec2(0,4), ivec2(1,5), ivec2(2,6), ivec2(3,7)
};

// list of faces; opposite face: (+3)%6
const int FACE_LIST[6][4] = {
	{0, 1, 5, 4}, {0, 3, 7, 4}, {0, 1, 2, 3},
	{2, 3, 7, 6}, {2, 1, 5, 6}, {4, 5, 6, 7}
};
const ivec3 FACE_DIR[6] = {
	ivec3(-1,0,0), ivec3(0,-1,0), ivec3(0,0,-1),
	ivec3(1,0,0), ivec3(0,1,0), ivec3(0,0,1)
};

// lookup tables for reconstruction
// http://paulbourke.net/geometry/polygonise/
const int EDGE_TABLE[256] = {
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc , 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};
const int TRIG_TABLE[256][16] = {
	{-1}, {0,8,3,-1}, {0,1,9,-1}, {1,8,3,9,8,1,-1}, {1,2,10,-1}, {0,8,3,1,2,10,-1}, {9,2,10,0,2,9,-1}, {2,8,3,2,10,8,10,9,8,-1}, {3,11,2,-1}, {0,11,2,8,11,0,-1}, {1,9,0,2,3,11,-1}, {1,11,2,1,9,11,9,8,11,-1}, {3,10,1,11,10,3,-1}, {0,10,1,0,8,10,8,11,10,-1}, {3,9,0,3,11,9,11,10,9,-1}, {9,8,10,10,8,11,-1},
	{4,7,8,-1}, {4,3,0,7,3,4,-1}, {0,1,9,8,4,7,-1}, {4,1,9,4,7,1,7,3,1,-1}, {1,2,10,8,4,7,-1}, {3,4,7,3,0,4,1,2,10,-1}, {9,2,10,9,0,2,8,4,7,-1}, {2,10,9,2,9,7,2,7,3,7,9,4,-1}, {8,4,7,3,11,2,-1}, {11,4,7,11,2,4,2,0,4,-1}, {9,0,1,8,4,7,2,3,11,-1}, {4,7,11,9,4,11,9,11,2,9,2,1,-1}, {3,10,1,3,11,10,7,8,4,-1}, {1,11,10,1,4,11,1,0,4,7,11,4,-1}, {4,7,8,9,0,11,9,11,10,11,0,3,-1}, {4,7,11,4,11,9,9,11,10,-1},
	{9,5,4,-1}, {9,5,4,0,8,3,-1}, {0,5,4,1,5,0,-1}, {8,5,4,8,3,5,3,1,5,-1}, {1,2,10,9,5,4,-1}, {3,0,8,1,2,10,4,9,5,-1}, {5,2,10,5,4,2,4,0,2,-1}, {2,10,5,3,2,5,3,5,4,3,4,8,-1}, {9,5,4,2,3,11,-1}, {0,11,2,0,8,11,4,9,5,-1}, {0,5,4,0,1,5,2,3,11,-1}, {2,1,5,2,5,8,2,8,11,4,8,5,-1}, {10,3,11,10,1,3,9,5,4,-1}, {4,9,5,0,8,1,8,10,1,8,11,10,-1}, {5,4,0,5,0,11,5,11,10,11,0,3,-1}, {5,4,8,5,8,10,10,8,11,-1},
	{9,7,8,5,7,9,-1}, {9,3,0,9,5,3,5,7,3,-1}, {0,7,8,0,1,7,1,5,7,-1}, {1,5,3,3,5,7,-1}, {9,7,8,9,5,7,10,1,2,-1}, {10,1,2,9,5,0,5,3,0,5,7,3,-1}, {8,0,2,8,2,5,8,5,7,10,5,2,-1}, {2,10,5,2,5,3,3,5,7,-1}, {7,9,5,7,8,9,3,11,2,-1}, {9,5,7,9,7,2,9,2,0,2,7,11,-1}, {2,3,11,0,1,8,1,7,8,1,5,7,-1}, {11,2,1,11,1,7,7,1,5,-1}, {9,5,8,8,5,7,10,1,3,10,3,11,-1}, {5,7,0,5,0,9,7,11,0,1,0,10,11,10,0,-1}, {11,10,0,11,0,3,10,5,0,8,0,7,5,7,0,-1}, {11,10,5,7,11,5,-1},
	{10,6,5,-1}, {0,8,3,5,10,6,-1}, {9,0,1,5,10,6,-1}, {1,8,3,1,9,8,5,10,6,-1}, {1,6,5,2,6,1,-1}, {1,6,5,1,2,6,3,0,8,-1}, {9,6,5,9,0,6,0,2,6,-1}, {5,9,8,5,8,2,5,2,6,3,2,8,-1}, {2,3,11,10,6,5,-1}, {11,0,8,11,2,0,10,6,5,-1}, {0,1,9,2,3,11,5,10,6,-1}, {5,10,6,1,9,2,9,11,2,9,8,11,-1}, {6,3,11,6,5,3,5,1,3,-1}, {0,8,11,0,11,5,0,5,1,5,11,6,-1}, {3,11,6,0,3,6,0,6,5,0,5,9,-1}, {6,5,9,6,9,11,11,9,8,-1},
	{5,10,6,4,7,8,-1}, {4,3,0,4,7,3,6,5,10,-1}, {1,9,0,5,10,6,8,4,7,-1}, {10,6,5,1,9,7,1,7,3,7,9,4,-1}, {6,1,2,6,5,1,4,7,8,-1}, {1,2,5,5,2,6,3,0,4,3,4,7,-1}, {8,4,7,9,0,5,0,6,5,0,2,6,-1}, {7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,-1}, {3,11,2,7,8,4,10,6,5,-1}, {5,10,6,4,7,2,4,2,0,2,7,11,-1}, {0,1,9,4,7,8,2,3,11,5,10,6,-1}, {9,2,1,9,11,2,9,4,11,7,11,4,5,10,6,-1}, {8,4,7,3,11,5,3,5,1,5,11,6,-1}, {5,1,11,5,11,6,1,0,11,7,11,4,0,4,11,-1}, {0,5,9,0,6,5,0,3,6,11,6,3,8,4,7,-1}, {6,5,9,6,9,11,4,7,9,7,11,9,-1},
	{10,4,9,6,4,10,-1}, {4,10,6,4,9,10,0,8,3,-1}, {10,0,1,10,6,0,6,4,0,-1}, {8,3,1,8,1,6,8,6,4,6,1,10,-1}, {1,4,9,1,2,4,2,6,4,-1}, {3,0,8,1,2,9,2,4,9,2,6,4,-1}, {0,2,4,4,2,6,-1}, {8,3,2,8,2,4,4,2,6,-1}, {10,4,9,10,6,4,11,2,3,-1}, {0,8,2,2,8,11,4,9,10,4,10,6,-1}, {3,11,2,0,1,6,0,6,4,6,1,10,-1}, {6,4,1,6,1,10,4,8,1,2,1,11,8,11,1,-1}, {9,6,4,9,3,6,9,1,3,11,6,3,-1}, {8,11,1,8,1,0,11,6,1,9,1,4,6,4,1,-1}, {3,11,6,3,6,0,0,6,4,-1}, {6,4,8,11,6,8,-1},
	{7,10,6,7,8,10,8,9,10,-1}, {0,7,3,0,10,7,0,9,10,6,7,10,-1}, {10,6,7,1,10,7,1,7,8,1,8,0,-1}, {10,6,7,10,7,1,1,7,3,-1}, {1,2,6,1,6,8,1,8,9,8,6,7,-1}, {2,6,9,2,9,1,6,7,9,0,9,3,7,3,9,-1}, {7,8,0,7,0,6,6,0,2,-1}, {7,3,2,6,7,2,-1}, {2,3,11,10,6,8,10,8,9,8,6,7,-1}, {2,0,7,2,7,11,0,9,7,6,7,10,9,10,7,-1}, {1,8,0,1,7,8,1,10,7,6,7,10,2,3,11,-1}, {11,2,1,11,1,7,10,6,1,6,7,1,-1}, {8,9,6,8,6,7,9,1,6,11,6,3,1,3,6,-1}, {0,9,1,11,6,7,-1}, {7,8,0,7,0,6,3,11,0,11,6,0,-1}, {7,11,6,-1},
	{7,6,11,-1}, {3,0,8,11,7,6,-1}, {0,1,9,11,7,6,-1}, {8,1,9,8,3,1,11,7,6,-1}, {10,1,2,6,11,7,-1}, {1,2,10,3,0,8,6,11,7,-1}, {2,9,0,2,10,9,6,11,7,-1}, {6,11,7,2,10,3,10,8,3,10,9,8,-1}, {7,2,3,6,2,7,-1}, {7,0,8,7,6,0,6,2,0,-1}, {2,7,6,2,3,7,0,1,9,-1}, {1,6,2,1,8,6,1,9,8,8,7,6,-1}, {10,7,6,10,1,7,1,3,7,-1}, {10,7,6,1,7,10,1,8,7,1,0,8,-1}, {0,3,7,0,7,10,0,10,9,6,10,7,-1}, {7,6,10,7,10,8,8,10,9,-1},
	{6,8,4,11,8,6,-1}, {3,6,11,3,0,6,0,4,6,-1}, {8,6,11,8,4,6,9,0,1,-1}, {9,4,6,9,6,3,9,3,1,11,3,6,-1}, {6,8,4,6,11,8,2,10,1,-1}, {1,2,10,3,0,11,0,6,11,0,4,6,-1}, {4,11,8,4,6,11,0,2,9,2,10,9,-1}, {10,9,3,10,3,2,9,4,3,11,3,6,4,6,3,-1}, {8,2,3,8,4,2,4,6,2,-1}, {0,4,2,4,6,2,-1}, {1,9,0,2,3,4,2,4,6,4,3,8,-1}, {1,9,4,1,4,2,2,4,6,-1}, {8,1,3,8,6,1,8,4,6,6,10,1,-1}, {10,1,0,10,0,6,6,0,4,-1}, {4,6,3,4,3,8,6,10,3,0,3,9,10,9,3,-1}, {10,9,4,6,10,4,-1},
	{4,9,5,7,6,11,-1}, {0,8,3,4,9,5,11,7,6,-1}, {5,0,1,5,4,0,7,6,11,-1}, {11,7,6,8,3,4,3,5,4,3,1,5,-1}, {9,5,4,10,1,2,7,6,11,-1}, {6,11,7,1,2,10,0,8,3,4,9,5,-1}, {7,6,11,5,4,10,4,2,10,4,0,2,-1}, {3,4,8,3,5,4,3,2,5,10,5,2,11,7,6,-1}, {7,2,3,7,6,2,5,4,9,-1}, {9,5,4,0,8,6,0,6,2,6,8,7,-1}, {3,6,2,3,7,6,1,5,0,5,4,0,-1}, {6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,-1}, {9,5,4,10,1,6,1,7,6,1,3,7,-1}, {1,6,10,1,7,6,1,0,7,8,7,0,9,5,4,-1}, {4,0,10,4,10,5,0,3,10,6,10,7,3,7,10,-1}, {7,6,10,7,10,8,5,4,10,4,8,10,-1},
	{6,9,5,6,11,9,11,8,9,-1}, {3,6,11,0,6,3,0,5,6,0,9,5,-1}, {0,11,8,0,5,11,0,1,5,5,6,11,-1}, {6,11,3,6,3,5,5,3,1,-1}, {1,2,10,9,5,11,9,11,8,11,5,6,-1}, {0,11,3,0,6,11,0,9,6,5,6,9,1,2,10,-1}, {11,8,5,11,5,6,8,0,5,10,5,2,0,2,5,-1}, {6,11,3,6,3,5,2,10,3,10,5,3,-1}, {5,8,9,5,2,8,5,6,2,3,8,2,-1}, {9,5,6,9,6,0,0,6,2,-1}, {1,5,8,1,8,0,5,6,8,3,8,2,6,2,8,-1}, {1,5,6,2,1,6,-1}, {1,3,6,1,6,10,3,8,6,5,6,9,8,9,6,-1}, {10,1,0,10,0,6,9,5,0,5,6,0,-1}, {0,3,8,5,6,10,-1}, {10,5,6,-1},
	{11,5,10,7,5,11,-1}, {11,5,10,11,7,5,8,3,0,-1}, {5,11,7,5,10,11,1,9,0,-1}, {10,7,5,10,11,7,9,8,1,8,3,1,-1}, {11,1,2,11,7,1,7,5,1,-1}, {0,8,3,1,2,7,1,7,5,7,2,11,-1}, {9,7,5,9,2,7,9,0,2,2,11,7,-1}, {7,5,2,7,2,11,5,9,2,3,2,8,9,8,2,-1}, {2,5,10,2,3,5,3,7,5,-1}, {8,2,0,8,5,2,8,7,5,10,2,5,-1}, {9,0,1,5,10,3,5,3,7,3,10,2,-1}, {9,8,2,9,2,1,8,7,2,10,2,5,7,5,2,-1}, {1,3,5,3,7,5,-1}, {0,8,7,0,7,1,1,7,5,-1}, {9,0,3,9,3,5,5,3,7,-1}, {9,8,7,5,9,7,-1},
	{5,8,4,5,10,8,10,11,8,-1}, {5,0,4,5,11,0,5,10,11,11,3,0,-1}, {0,1,9,8,4,10,8,10,11,10,4,5,-1}, {10,11,4,10,4,5,11,3,4,9,4,1,3,1,4,-1}, {2,5,1,2,8,5,2,11,8,4,5,8,-1}, {0,4,11,0,11,3,4,5,11,2,11,1,5,1,11,-1}, {0,2,5,0,5,9,2,11,5,4,5,8,11,8,5,-1}, {9,4,5,2,11,3,-1}, {2,5,10,3,5,2,3,4,5,3,8,4,-1}, {5,10,2,5,2,4,4,2,0,-1}, {3,10,2,3,5,10,3,8,5,4,5,8,0,1,9,-1}, {5,10,2,5,2,4,1,9,2,9,4,2,-1}, {8,4,5,8,5,3,3,5,1,-1}, {0,4,5,1,0,5,-1}, {8,4,5,8,5,3,9,0,5,0,3,5,-1}, {9,4,5,-1},
	{4,11,7,4,9,11,9,10,11,-1}, {0,8,3,4,9,7,9,11,7,9,10,11,-1}, {1,10,11,1,11,4,1,4,0,7,4,11,-1}, {3,1,4,3,4,8,1,10,4,7,4,11,10,11,4,-1}, {4,11,7,9,11,4,9,2,11,9,1,2,-1}, {9,7,4,9,11,7,9,1,11,2,11,1,0,8,3,-1}, {11,7,4,11,4,2,2,4,0,-1}, {11,7,4,11,4,2,8,3,4,3,2,4,-1}, {2,9,10,2,7,9,2,3,7,7,4,9,-1}, {9,10,7,9,7,4,10,2,7,8,7,0,2,0,7,-1}, {3,7,10,3,10,2,7,4,10,1,10,0,4,0,10,-1}, {1,10,2,8,7,4,-1}, {4,9,1,4,1,7,7,1,3,-1}, {4,9,1,4,1,7,0,8,1,8,7,1,-1}, {4,0,3,7,4,3,-1}, {4,8,7,-1},
	{9,10,8,10,11,8,-1}, {3,0,9,3,9,11,11,9,10,-1}, {0,1,10,0,10,8,8,10,11,-1}, {3,1,10,11,3,10,-1}, {1,2,11,1,11,9,9,11,8,-1}, {3,0,9,3,9,11,1,2,9,2,11,9,-1}, {0,2,11,8,0,11,-1}, {3,2,11,-1}, {2,3,8,2,8,10,10,8,9,-1}, {9,10,2,0,9,2,-1}, {2,3,8,2,8,10,0,1,8,1,10,8,-1}, {1,10,2,-1}, {1,3,8,9,1,8,-1}, {0,9,1,-1}, {0,3,8,-1}, {-1}
};

// calculate index for table lookup
static int CalcIndex(const float v[8]) {
	if (isnan(v[0] + v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[7])) return 0;
	return int(v[0] < 0) |
		(int(v[1] < 0) << 1) |
		(int(v[2] < 0) << 2) |
		(int(v[3] < 0) << 3) |
		(int(v[4] < 0) << 4) |
		(int(v[5] < 0) << 5) |
		(int(v[6] < 0) << 6) |
		(int(v[7] < 0) << 7);
}

// linear interpolation on an edge
static vec3 getInterpolation(const vec3 pos[8], const float val[8], int i) {
	float v0 = val[EDGE_LIST[i].x];
	float v1 = val[EDGE_LIST[i].y];
	vec3 p0 = pos[EDGE_LIST[i].x];
	vec3 p1 = pos[EDGE_LIST[i].y];
	if (v0 > v1) std::swap(v0, v1), std::swap(p0, p1);  // make it consistant
	float t = v0 / (v0 - v1);
	return p0 * (1 - t) + p1 * t;
};
static ivec3 getInterpolation(const ivec3 pos[8], const float val[8], int i) {
	double v0 = (double)val[EDGE_LIST[i].x];
	double v1 = (double)val[EDGE_LIST[i].y];
	ivec3 p0 = pos[EDGE_LIST[i].x];
	ivec3 p1 = pos[EDGE_LIST[i].y];
	if (v0 > v1) std::swap(v0, v1), std::swap(p0, p1);  // make it consistant
	double t = v0 / (v0 - v1);
	ivec3 p;
	for (int i = 0; i < 3; i++) {
		((int*)&p)[i] = (int)round(((int*)&p0)[i] * (1.0 - t) + ((int*)&p1)[i] * t);
	}
	return p;
};



/* DATA */

std::function<float(vec3)> fun;  // float/float fun(vec3)
vec3 P0, P1;  // search interval
ivec3 SEARCH_DIF;  // initial diff
int PLOT_DEPTH;  // octatree depth
int PLOT_SIZE;  // 1 << PLOT_DEPTH
ivec3 GRID_SIZE;  // SEARCH_DIF * PLOT_SIZE
int EDGE_ROUNDING;  // subdivide an edge into this number of intervals and round coordinates to integer
ivec3 MESH_SIZE;  // precision of mesh: (P1-P0)/MESH_SIZE
void calcParams() {
	if (PLOT_DEPTH <= 0)
		throw "Plot depth should be a positive integer.";
	PLOT_SIZE = 1 << PLOT_DEPTH;
	GRID_SIZE = SEARCH_DIF * PLOT_SIZE;
	MESH_SIZE = GRID_SIZE * EDGE_ROUNDING;
	if (EDGE_ROUNDING >= 0x100)
		throw "Edge rounding beyond 8 bit integer limit.";
	/*if (max(max(MESH_SIZE.x, MESH_SIZE.y), MESH_SIZE.z) >= 0x10000)
		throw "Grid size beyond 16 bit integer limit.";*/
}

// position ID to position
vec3 i2f(ivec3 p) {
	vec3 d = vec3(p) / vec3(GRID_SIZE);
	return P0 * (vec3(1.) - d) + P1 * d;
}


/* OCTATREE */

// octatree
float getSample_global(ivec3 p);
class octatree_node {
public:
	float v[8]; // 32 bytes, only the first one really matters
	octatree_node *c[8];  // 32 or 64 bytes, child nodes
	octatree_node *parent;  // 4 or 8 bytes, parent node
	ivec3 _p; // 12 bytes, vertex ID
	ivec3 p(int i) { return _p + VERTICE_LIST[i] * size; }
	int size;  // 4 bytes, top right: p+ivec3(size)
	int index;  // 4 bytes, calculated according to signs of v for table lookup
	bool has_sign_change[6];  // 6 bytes, indicate whether there is a sign change in each face
	bool edge_checked[6];  // 6 bytes, used in looking for missed samples, indicate whether the face is already checked
	octatree_node(int size = 0, ivec3 p = ivec3(-1)) {
		for (int i = 0; i < 8; i++) this->_p = p;
		for (int i = 0; i < 8; i++) v[i] = NAN;
		this->size = size;
		this->index = -1;
		this->parent = nullptr;
		for (int i = 0; i < 8; i++) c[i] = nullptr;
		for (int i = 0; i < 6; i++) has_sign_change[i] = edge_checked[i] = false;
	}
	~octatree_node() {
		for (int i = 0; i < 8; i++) if (c[i]) {
			delete c[i]; c[i] = 0;
		}
	}
	float getSample(ivec3 q) {
		if (q == _p) {
			if (isnan(v[0])) {
				if (c[0]) {
					v[0] = c[0]->getSample(q);
				}
				else v[0] = (float)fun(i2f(_p));
			}
			return v[0];
		}
		ivec3 d = (q - _p) / (size >> 1);
		int i = VERTICE_LIST_INV[d.x][d.y][d.z];
		if (!c[i]) {
			c[i] = new octatree_node(size / 2, _p + d * (size / 2));
			c[i]->parent = this;
		}
		return c[i]->getSample(q);
	}
	octatree_node* getGrid(ivec3 q, int sz) {
		if (q == _p && sz == size) {
			if (isnan(v[0])) {
				v[0] = getSample_global(_p);
			}
			for (int i = 1; i < 8; i++) if (isnan(v[i])) {
				v[i] = getSample_global(p(i));
			}
			return this;
		}
		ivec3 d = (q - _p) / (size >> 1);
		int i = VERTICE_LIST_INV[d.x][d.y][d.z];
		if (!c[i]) {
			c[i] = new octatree_node(size / 2, _p + d * (size / 2));
			c[i]->parent = this;
		}
		return c[i]->getGrid(q, sz);
	}
	int calcIndex() {
		return (index = CalcIndex(v));
	}

	// grid subdivision
	void subdivide(int sz);

};

octatree_node*** octatree_grid = 0;  // a grid [x][y][z]
void create_octatree() {  // sample tree initialization
	octatree_grid = new octatree_node**[SEARCH_DIF.x + 1];
	for (int x = 0; x <= SEARCH_DIF.x; x++) {
		octatree_grid[x] = new octatree_node*[SEARCH_DIF.y + 1];
		for (int y = 0; y <= SEARCH_DIF.y; y++) {
			octatree_grid[x][y] = new octatree_node[SEARCH_DIF.z + 1];
			for (int z = 0; z <= SEARCH_DIF.z; z++) {
				octatree_grid[x][y][z] = octatree_node(PLOT_SIZE, ivec3(x, y, z)*PLOT_SIZE);
			}
		}
	}
}
void destroy_octatree() {  // sample tree destruction
	for (int x = 0; x <= SEARCH_DIF.x; x++) {
		for (int y = 0; y <= SEARCH_DIF.y; y++) {
			delete[] octatree_grid[x][y];
		}
		delete octatree_grid[x];
	}
	delete octatree_grid;
	octatree_grid = 0;
}
float getSample_global(ivec3 p) {  // access a sample on the sample tree
	ivec3 pi = p / PLOT_SIZE;
	return octatree_grid[pi.x][pi.y][pi.z].getSample(p);
}
octatree_node* getGrid_global(ivec3 p, int sz) {
	ivec3 pi = p / PLOT_SIZE;
	return octatree_grid[pi.x][pi.y][pi.z].getGrid(p, sz);
}


void octatree_node::subdivide(int sz) {

	for (int u = 0; u < 8; u++) {
		if (!c[u]) c[u] = new octatree_node(sz / 2, _p + VERTICE_LIST[u] * (sz / 2));
		c[u]->parent = this;
	}

	float samples[27];
	samples[0] = v[0];
	samples[2] = v[1];
	samples[4] = v[2];
	samples[6] = v[3];
	samples[18] = v[4];
	samples[20] = v[5];
	samples[22] = v[6];
	samples[24] = v[7];
	samples[1] = isnan(c[1]->v[0]) ? (float)fun(i2f(c[1]->_p)) : c[1]->v[0];
	samples[3] = getSample_global(c[1]->p(2));
	samples[5] = getSample_global(c[3]->p(2));
	samples[7] = isnan(c[3]->v[0]) ? (float)fun(i2f(c[3]->_p)) : c[3]->v[0];
	samples[8] = isnan(c[2]->v[0]) ? (float)fun(i2f(c[2]->_p)) : c[2]->v[0];
	samples[9] = isnan(c[4]->v[0]) ? (float)fun(i2f(c[4]->_p)) : c[4]->v[0];
	samples[10] = isnan(c[5]->v[0]) ? (float)fun(i2f(c[5]->_p)) : c[5]->v[0];
	samples[11] = getSample_global(c[5]->p(1));
	samples[12] = getSample_global(c[5]->p(2));
	samples[13] = getSample_global(c[6]->p(2));
	samples[14] = getSample_global(c[7]->p(2));
	samples[15] = getSample_global(c[7]->p(3));
	samples[16] = isnan(c[7]->v[0]) ? (float)fun(i2f(c[7]->_p)) : c[7]->v[0];
	samples[17] = isnan(c[6]->v[0]) ? (float)fun(i2f(c[6]->_p)) : c[6]->v[0];
	samples[19] = getSample_global(c[4]->p(5));
	samples[21] = getSample_global(c[5]->p(6));
	samples[23] = getSample_global(c[7]->p(6));
	samples[25] = getSample_global(c[4]->p(7));
	samples[26] = getSample_global(c[4]->p(6));

	const static int SUBDIV_LOOKUP[8][8] = {
		{0, 1, 8, 7, 9, 10, 17, 16},
		{1, 2, 3, 8, 10, 11, 12, 17},
		{8, 3, 4, 5, 17, 12, 13, 14},
		{7, 8, 5, 6, 16, 17, 14, 15},
		{9, 10, 17, 16, 18, 19, 26, 25},
		{10, 11, 12, 17, 19, 20, 21, 26},
		{17, 12, 13, 14, 26, 21, 22, 23},
		{16, 17, 14, 15, 25, 26, 23, 24}
	};
	for (int u = 0; u < 8; u++) for (int v = 0; v < 8; v++) {
		c[u]->v[v] = samples[SUBDIV_LOOKUP[u][v]];
	}

}


/* CALL FUNCTIONS */

std::vector<octatree_node*> cells;

// octatree main function, construct octatree and march cubes to @cells
void octatree_main() {

	// initialize octatree root
	timeEventStart("sample octatree root");
	create_octatree();
	for (int x = 0; x <= SEARCH_DIF.x; x++) {
		for (int y = 0; y <= SEARCH_DIF.y; y++) {
			for (int z = 0; z <= SEARCH_DIF.z; z++) {
				octatree_grid[x][y][z].v[0] = (float)fun(i2f(octatree_grid[x][y][z]._p = ivec3(x, y, z)*PLOT_SIZE));
			}
		}
	}
	for (int x = 0; x < SEARCH_DIF.x; x++) {
		for (int y = 0; y < SEARCH_DIF.y; y++) {
			for (int z = 0; z < SEARCH_DIF.z; z++) {
				for (int u = 1; u < 8; u++) {
					ivec3 p = ivec3(x, y, z) + VERTICE_LIST[u];
					octatree_grid[x][y][z].v[u] = octatree_grid[p.x][p.y][p.z].v[0];
					octatree_grid[x][y][z].calcIndex();
				}
			}
		}
	}

	// initial sample cells
	int CELL_PLOT_SIZE = PLOT_SIZE;
	cells.clear();
	for (int x = 0; x < SEARCH_DIF.x; x++) {
		for (int y = 0; y < SEARCH_DIF.y; y++) {
			for (int z = 0; z < SEARCH_DIF.z; z++) {
				octatree_node *n = &octatree_grid[x][y][z];
				if (TRIG_TABLE[n->index][0] != -1) {
					cells.push_back(n);
				}
			}
		}
	}
	timeEventEnd("sample octatree root");

	// subdivide grid cells
	for (int layer = 1, size = CELL_PLOT_SIZE; size > 1; layer++) {
		std::string event_name = "subdivide octatree layer " + std::to_string(layer) + "/" + std::to_string(PLOT_DEPTH);
		timeEventStart(event_name);

		std::vector<octatree_node*> new_cells;
		int s2 = size / 2;
		for (int i = 0; i < (int)cells.size(); i++) {
			octatree_node* ci = cells[i];
			ci->subdivide(size);
			for (int u = 0; u < 8; u++) {
				if (TRIG_TABLE[ci->c[u]->calcIndex()][0] != -1) {
					new_cells.push_back(ci->c[u]);
				}
			}
		}
		cells = new_cells;

		size >>= 1;

		// add missed samples (similar to BFS)
		for (int i = 0; i < (int)cells.size(); i++) {
			octatree_node* ci = cells[i];
			for (int u = 0; u < 6; u++) {
				ci->has_sign_change[u] = ((int)signbit(ci->v[FACE_LIST[u][0]]) + (int)signbit(ci->v[FACE_LIST[u][1]]) + (int)signbit(ci->v[FACE_LIST[u][2]]) + (int)signbit(ci->v[FACE_LIST[u][3]])) % 4 != 0;
			}
		}
		for (int i = 0; i < (int)cells.size(); i++) {
			octatree_node* ci = cells[i];
			for (int u = 0; u < 6; u++) if (ci->has_sign_change[u] && !ci->edge_checked[u]) {
				ivec3 nb_p = ci->p(0) + FACE_DIR[u] * ci->size;
				if (nb_p.x >= 0 && nb_p.y >= 0 && nb_p.z >= 0 && nb_p.x < GRID_SIZE.x && nb_p.y < GRID_SIZE.y && nb_p.z < GRID_SIZE.z) {
					octatree_node* nb = getGrid_global(nb_p, ci->size);
					if (!nb->has_sign_change[(u + 3) % 6]) {
						for (int u = 0; u < 6; u++)
							nb->has_sign_change[u] = ((int)signbit(nb->v[FACE_LIST[u][0]]) + (int)signbit(nb->v[FACE_LIST[u][1]]) + (int)signbit(nb->v[FACE_LIST[u][2]]) + (int)signbit(nb->v[FACE_LIST[u][3]])) % 4 != 0;
						cells.push_back(nb);
					}
					nb->edge_checked[(u + 3) % 6] = true;
				}
				ci->edge_checked[u] = true;
			}
		}

		timeEventEnd(event_name);
	}
}


// export a discrete list of triangles
void triangulate(std::vector<triangle_3d> &trigs) {
	octatree_main();
	for (int i = 0, cn = (int)cells.size(); i < cn; i++) {
		vec3 p[8];
		for (int j = 0; j < 8; j++) p[j] = i2f(cells[i]->p(j));
		auto v_table = TRIG_TABLE[CalcIndex(cells[i]->v)];
		for (int u = 0; ; u += 3) {
			if (v_table[u] == -1) break;
			vec3 a = getInterpolation(p, cells[i]->v, v_table[u]);
			vec3 b = getInterpolation(p, cells[i]->v, v_table[u + 1]);
			vec3 c = getInterpolation(p, cells[i]->v, v_table[u + 2]);
			trigs.push_back(triangle_3d(a, b, c));
		}
	}
	destroy_octatree();
}


// export the entire tree. see `README.md` for details
// return vector may not be consistent due to the use of memory address
std::vector<uint8_t> to_buffer(vec3(*colorf)(vec3), bool shrink_grid = true) {
	octatree_main();

	// triangles + bottom layer
	timeEventStart("reconstruct triangles");
	struct Triangle {
		ivec3 a, b, c;
	};
	struct TriangleNode {
		octatree_node* p;  // cell ID
		ivec3 po;  // origin of coordinates of triangles
		int n = 0;  // number of triangles
		int t[5] = { -1, -1, -1, -1, -1 };
	};
	std::vector<Triangle> triangles;
	std::vector<TriangleNode> bottom_layer;
	for (int i = 0; i < (int)cells.size(); i++) {
		ivec3 p[8];
		for (int j = 0; j < 8; j++) p[j] = (cells[i]->_p + VERTICE_LIST[j]) * EDGE_ROUNDING;
		auto v_table = TRIG_TABLE[CalcIndex(cells[i]->v)];
		if (v_table[0] != -1) {
			TriangleNode n;
			n.p = cells[i];
			n.po = cells[i]->_p;
			for (int u = 0; ; u += 3) {
				if (v_table[u] == -1) break;
				ivec3 a = getInterpolation(p, cells[i]->v, v_table[u]);
				ivec3 b = getInterpolation(p, cells[i]->v, v_table[u + 1]);
				ivec3 c = getInterpolation(p, cells[i]->v, v_table[u + 2]);
				if (a != b && a != c && b != c) {
					n.t[n.n++] = (int)triangles.size();
					triangles.push_back(Triangle{ a, b, c });
				}
			}
			if (n.n != 0)
				bottom_layer.push_back(n);
		}
	}
	//std::sort(bottom_layer.begin(), bottom_layer.end(), [](TriangleNode a, TriangleNode b) { return a.p < b.p; });
	timeEventEnd("reconstruct triangles");

	// middle layers
	timeEventStart("restore layers");
	struct Node {
		octatree_node *p;  // cell ID
		int c[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
	};
	struct Pi {
		octatree_node *p;
		int i;
	};
	std::vector<std::vector<Node>> middle_layers;
	std::vector<Node> layer;
	for (int i = 0; i < (int)bottom_layer.size(); i++) {  // fill using bottom layer
		Node n;
		n.p = bottom_layer[i].p;
		layer.push_back(n);
	}
	for (int h = PLOT_DEPTH; h > 0; h--) {
		std::vector<Pi> pc_map;  // parent, child
		for (int i = 0; i < (int)layer.size(); i++) {
			pc_map.push_back(Pi{ layer[i].p->parent, i });
		}
		std::sort(pc_map.begin(), pc_map.end(), [](Pi a, Pi b) { return a.p < b.p; });
		std::vector<Node> prev_layer;
		for (Pi pc : pc_map) {
			octatree_node *p = pc.p;
			int i = pc.i;
			if (prev_layer.empty() || prev_layer.back().p != p) {
				Node n;
				n.p = p;
				for (int _ = 0; _ < 8; _++) if (p->c[_] == layer[i].p) n.c[_] = i;
				prev_layer.push_back(n);
			}
			else {
				Node n = prev_layer.back();
				for (int _ = 0; _ < 8; _++) if (p->c[_] == layer[i].p) n.c[_] = i;
				prev_layer.back() = n;
			}
		}
		middle_layers.push_back(prev_layer);
		layer = prev_layer;
	}
	timeEventEnd("restore layers");

	// top layer (grid)
	timeEventStart("reconstruct grid");
	std::vector<int> top_layer;
	for (int i = SEARCH_DIF.x*SEARCH_DIF.y*SEARCH_DIF.z; i--;) top_layer.push_back(-1);
	std::map<octatree_node*, int> pimap;
	for (int i = 0; i < (int)middle_layers.back().size(); i++) pimap[middle_layers.back()[i].p] = i;
	for (int z = 0; z < SEARCH_DIF.z; z++) {
		for (int y = 0; y < SEARCH_DIF.y; y++) {
			for (int x = 0; x < SEARCH_DIF.x; x++) {
				int j = (z * SEARCH_DIF.y + y) * SEARCH_DIF.x + x;
				octatree_node *p = &octatree_grid[x][y][z];
				if (pimap.find(p) != pimap.end()) top_layer[j] = pimap[p];
			}
		}
	}
	timeEventEnd("reconstruct grid");

	// release memory
	timeEventStart("destroy octatree");
	destroy_octatree();
	timeEventEnd("destroy octatree");

	// shrink grid
	timeEventStart("shrink grid");
	ivec3 GRID_DIF = SEARCH_DIF;
	while (shrink_grid && (GRID_DIF.x % 2 == 0 && GRID_DIF.y % 2 == 0 && GRID_DIF.z % 2 == 0)) {

		std::vector<int> cell_ptr_map;
		layer.clear();
		for (int z = 0; z < GRID_DIF.z; z++) for (int y = 0; y < GRID_DIF.y; y++) for (int x = 0; x < GRID_DIF.x; x++) {
			int j = (z * GRID_DIF.y + y) * GRID_DIF.x + x;
			if (top_layer[j] != -1) {
				cell_ptr_map.push_back((int)layer.size());
				Node n = middle_layers.back()[top_layer[j]];
				layer.push_back(n);
			}
			else cell_ptr_map.push_back(-1);
		}
		middle_layers.back() = layer;

		GRID_DIF /= 2;
		std::vector<int> prev_top_layer;
		layer.clear();
		for (int z = 0; z < GRID_DIF.z; z++) for (int y = 0; y < GRID_DIF.y; y++) for (int x = 0; x < GRID_DIF.x; x++) {
			Node n;
			bool is_empty = true;
			for (int i = 0; i < 8; i++) {
				ivec3 p = ivec3(x, y, z) * 2 + VERTICE_LIST[i];
				int j = (p.z * 2 * GRID_DIF.y + p.y) * 2 * GRID_DIF.x + p.x;
				n.c[i] = cell_ptr_map[j];
				if (n.c[i] != -1) is_empty = false;
			}
			if (is_empty) prev_top_layer.push_back(-1);
			else prev_top_layer.push_back((int)layer.size()), layer.push_back(n);
		}
		top_layer = prev_top_layer;
		middle_layers.push_back(layer);
	}
	timeEventEnd("shrink grid");

	// triangles to mesh
	timeEventStart("triangles to mesh");
	struct PosI {
		ivec3 p;
		int i;
	};
	std::vector<PosI> pi_map;
	pi_map.reserve(triangles.size() * 3);
	for (int i = 0; i < (int)triangles.size(); i++) {
		pi_map.push_back(PosI{ triangles[i].a, i });
		pi_map.push_back(PosI{ triangles[i].b, i });
		pi_map.push_back(PosI{ triangles[i].c, i });
	}
	std::sort(pi_map.begin(), pi_map.end(), [](PosI a, PosI b) {
		return a.p.z == b.p.z ? a.p.y == b.p.y ? a.p.x < b.p.x : a.p.y < b.p.y : a.p.z < b.p.z;
	});
	std::vector<ivec3> vertices;
	std::vector<ivec3> triangles_v;
	triangles_v.resize(triangles.size());
	for (int j = 0; j < (int)pi_map.size(); j++) {
		if (vertices.empty() || vertices.back() != pi_map[j].p) {
			vertices.push_back(pi_map[j].p);
		}
		int i = pi_map[j].i;
		for (int _ = 0; _ < 3; _++) {
			if (((ivec3*)&triangles[i])[_] == pi_map[j].p)
				((int*)&triangles_v[i])[_] = (int)vertices.size() - 1;
		}
	}
	timeEventEnd("triangles to mesh");
	timeEventStart("coloring");
	std::vector<vec3> vert_cols;
	for (ivec3 p : vertices) {
		vec3 q = mix(P0, P1, vec3(p) / vec3(MESH_SIZE));
		vert_cols.push_back(colorf(q));
	}
	timeEventEnd("coloring");

	// put layers together
	timeEventStart("encode buffer");
	std::vector<uint8_t> result;
	auto pushUint16 = [](std::vector<uint8_t> &buf, uint16_t val) {
		buf.push_back((uint8_t)(val & 0xff));
		buf.push_back((uint8_t)(val >> 8));
	};
	auto pushUint32 = [](std::vector<uint8_t> &buf, uint32_t val) {
		buf.push_back((uint8_t)(val & 0xff));
		buf.push_back((uint8_t)((val >> 8) & 0xff));
		buf.push_back((uint8_t)((val >> 16) & 0xff));
		buf.push_back((uint8_t)(val >> 24));
	};
	int add = 4 * (int)top_layer.size();  // top layer (grid)
	for (int i = 0; i < (int)top_layer.size(); i++) {
		int data = top_layer[i];
		int p = data == -1 ? 0 : 4 * 8 * data + add;
		pushUint32(result, p);
	}
	std::vector<int> bottom_layer_add;  // pointers to bottom layer
	bottom_layer_add.push_back(0);
	for (int i = 0; i < (int)bottom_layer.size(); i++) {
		int sz = 1 + 12 * bottom_layer[i].n;
		bottom_layer_add.push_back(bottom_layer_add.back() + sz);
	}
	for (int j = (int)middle_layers.size(); j--;) {  // middle layers
		add += 4 * 8 * (int)middle_layers[j].size();
		for (int i = 0; i < (int)middle_layers[j].size(); i++) {
			for (int _ = 0; _ < 8; _++) {
				int data = middle_layers[j][i].c[_];
				int p = data == -1 ? 0 : (j > 0 ? 4 * 8 * data : bottom_layer_add[data]) + add;
				pushUint32(result, p);
			}
		}
	}
	add += bottom_layer_add.back();  // bottom layer
	for (int i = 0; i < (int)bottom_layer.size(); i++) {
		int data = bottom_layer[i].n;
		result.push_back((uint8_t)data);
		for (int j = 0; j < bottom_layer[i].n; j++) {
			int trig_i = bottom_layer[i].t[j];
			vec3 col = vec3(0.0);
			for (int _ = 0; _ < 3; _++) {
				int vi = ((int*)&triangles_v[trig_i])[_];
				ivec3 v = vertices[vi] - bottom_layer[i].po * EDGE_ROUNDING;
				result.push_back((uint8_t)(v.x));
				result.push_back((uint8_t)(v.y));
				result.push_back((uint8_t)(v.z));
				col += vert_cols[vi];
			}
			ivec3 rgb = ivec3(255.0*clamp(col / 3.0, 0.0, 1.0) + 0.5);
			result.push_back((uint8_t)rgb.x);
			result.push_back((uint8_t)rgb.y);
			result.push_back((uint8_t)rgb.z);
		}
	}
	timeEventEnd("encode buffer");
	return result;
}


ScalarFieldTriangulator_octatree_PRIVATE_END_


// octatree
template<typename Fun>
std::vector<triangle_3d> octatree(Fun fun, vec3 p0, vec3 p1, ivec3 dif, int plot_depth) {
	__private__::fun = fun;
	__private__::P0 = p0, __private__::P1 = p1;
	__private__::SEARCH_DIF = dif, __private__::PLOT_DEPTH = plot_depth;
	__private__::EDGE_ROUNDING = 0;
	__private__::calcParams();
	std::vector<triangle_3d> Trigs;
	__private__::triangulate(Trigs);
	return Trigs;
}

// octatree to buffer
template<typename Fun>
std::vector<uint8_t> octatree_buffer(Fun fun, vec3 p0, vec3 p1, ivec3 dif, int plot_depth, int edge_rounding,
	vec3(*colorf)(vec3)) {
	__private__::fun = fun;
	__private__::P0 = p0, __private__::P1 = p1;
	__private__::SEARCH_DIF = dif, __private__::PLOT_DEPTH = plot_depth;
	__private__::EDGE_ROUNDING = edge_rounding;
	__private__::calcParams();
	return __private__::to_buffer(colorf, true);
}


ScalarFieldTriangulator_octatree_END_

#undef ScalarFieldTriangulator_octatree_PRIVATE_
#undef ScalarFieldTriangulator_octatree__PRIVATE
