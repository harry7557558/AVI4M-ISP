// Implement the following functions or global variables in .glsl.cpp:

// vec4 map(vec3 p, bool col_required);  // return (r, g, b, sd)
// vec3 BOX_RADIUS;  // march between -BOX_RADIUS and BOX_RADIUS
// float STEP;  // help determining size of initial marching box
// float MAX_STEP;  // help determining the marching precision

#include <cstdio>
#include <chrono>
#include <string>

#include "octatree.h"
#include "trigs2mesh.h"
#include "ply_writer.h"


namespace GLSL {

	float iTime = 0.0;
	int iFrame = 0;
	vec3 iResolution = vec3(0, 0, 1);
	vec4 iMouse = vec4(0, 0, 0, 0);

	vec4 gl_FragCoord;

#include "../.glsl.cpp"

#ifndef BOX_RADIUS
#define BOX_RADIUS GLSL::BOX_RADIUS
#endif
#ifndef STEP
#define STEP GLSL::STEP
#endif
#ifndef MIN_STEP
#define MIN_STEP GLSL::MIN_STEP
#endif
}


int main(int argc, char* argv[]) {

	// compute parameters
	ivec3 box0 = ivec3(
		1 << (int)ceil(log2(2.0*BOX_RADIUS.x / STEP)),
		1 << (int)ceil(log2(2.0*BOX_RADIUS.y / STEP)),
		1 << (int)ceil(log2(2.0*BOX_RADIUS.z / STEP))
	);
	int boxe = (int)round(log2(STEP / MIN_STEP));
	boxe = std::stoi(argv[2]);
	vec3 cell_size = 2.0*BOX_RADIUS / (vec3(box0) * exp2(boxe));
	float epsilon = 0.001f * std::min({ cell_size.x, cell_size.y, cell_size.z });
	printf("%d %d %d  %d  epsilon=%.2g\n", box0.x, box0.y, box0.z, 1 << boxe, epsilon);

	// marching cube
	int eval_count = 0;
	auto triangulate_start = std::chrono::high_resolution_clock::now();
	std::vector<triangle_3d> trigs = ScalarFieldTriangulator_octatree::octatree(
		[&](vec3 p) { eval_count++; return GLSL::map(p, false).w; },
		-BOX_RADIUS, BOX_RADIUS,
		box0, boxe
	);
	float triangulate_time = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - triangulate_start).count();
	printf("%.1fms, %d evaluations\n", 1000.0f*triangulate_time, eval_count);
	printf("%d triangles => ", (int)trigs.size());

	// triangles to mesh
	std::vector<vec3> vertices;
	std::vector<ivec3> faces;
	TrigsToMesh(trigs, epsilon, vertices, faces);
	printf("%d vertices, %d faces\n", (int)vertices.size(), (int)faces.size());

	// color mesh
	auto coloring_start = std::chrono::high_resolution_clock::now();
	std::vector<vec3> colors;
	colors.resize((int)vertices.size());
	for (int i = 0; i < (int)vertices.size(); i++) {
		colors[i] = GLSL::map(vertices[i], true).xyz();
	}
	float coloring_time = std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - coloring_start).count();
	printf("%.1fms coloring\n", 1000.0f*coloring_time);

	// output
	WritePLY(argv[1],
		&vertices[0], (int)vertices.size(),
		&faces[0], (int)faces.size(),
		&colors[0]
	);

	return 0;
}
