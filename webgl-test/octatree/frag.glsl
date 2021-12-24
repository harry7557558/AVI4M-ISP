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
uniform highp sampler2D fTreeBuffer;

#define ZERO min(iFrame,0)

#define P0 vec3(-2.0, -2.0, -2.0) /* min coordinates of grid */
#define P1 vec3(2.0, 2.0, 2.0) /* max coordinates of grid */
#define GRID_DIF ivec3(1, 1, 1) /* initial grid size, at least one odd component */
#define PLOT_DEPTH 6 /* depth of the tree */
#define GRID_SIZE (GRID_DIF*(1<<PLOT_DEPTH))
#define EDGE_ROUNDING 255 /* divide edge into # intervals and round to integer coordinate */
#define MESH_SIZE (GRID_SIZE*EDGE_ROUNDING)


// sample integer
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

// test if float is faster than int
float getUint8f(float i) {
	float j = floor(i/4.);
    vec2 pos = vec2(mod(j,float(TREEBUFFER_SIZE)), floor(j/float(TREEBUFFER_SIZE)));
    vec4 sp = 255.0*texelFetch(fTreeBuffer, ivec2(pos), 0);
	float f = i - 4.0*j;
    return f==0. ? sp.x : f==1. ? sp.y : f==2. ? sp.z : sp.w;
}
float getUint32f(float i) {
	float j = floor(i/4.);
    vec2 pos = vec2(mod(j,float(TREEBUFFER_SIZE)), floor(j/float(TREEBUFFER_SIZE)));
    vec4 sp = 255.0*texelFetch(fTreeBuffer, ivec2(pos), 0);
	return sp.x + 256. * (sp.y + 256. * (sp.z + 256. * sp.w));
}
vec3 getUvec3f(float i) {
	float x = getUint8f(i);
	float y = getUint8f(i + 1.);
	float z = getUint8f(i + 2.);
	return vec3(x, y, z);
}


#include "intersector2.glsl"


#if 0
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
	//float i = gl_FragCoord.y*iResolution.x+gl_FragCoord.x;
	//fragColor = vec4(vec3(getUint8f(i)),1.0);
}
