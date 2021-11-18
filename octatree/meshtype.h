#pragma once

#include "../glsl2cpp/glslmath.h"

struct triangle_3d {
	vec3 v[3];
	triangle_3d() {}
	triangle_3d(vec3 a, vec3 b, vec3 c) { v[0] = a, v[1] = b, v[2] = c; }
	vec3 operator[](int i) const { return v[i]; }
	vec3& operator[](int i) { return v[i]; }
};
