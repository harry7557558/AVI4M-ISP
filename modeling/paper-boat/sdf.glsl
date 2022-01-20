#include "common.glsl"

#define NO_MAP
#include "sdf_boat.glsl"
#include "sdf_coraltrees.glsl"
#include "sdf_crab.glsl"

#ifndef __cplusplus
throw "Prevent this file from being compiled as GLSL because it freezes the IDE.";
#endif

vec4 map(vec3 p, bool col_required) {
    vec4 d = vec4(0,0,0, 1e3);

d = smin(d, mapBoat(((p-vec3(0.00,0.00,-0.97)))/vec3(1.00,1.00,1.00), col_required)*vec4(1,1,1,1.00), 0.01);
d = smin(d, mapCoralTree01((rotx(-1.00*PI)*roty(-0.02*PI)*rotz(-0.02*PI)*(p-vec3(1.30,-0.11,0.04)))/vec3(-0.75,-0.75,-0.81), col_required)*vec4(1,1,1,0.75), 0.01);
d = smin(d, mapCoralTree02((rotx(0.04*PI)*roty(0.02*PI)*rotz(-0.05*PI)*(p-vec3(1.19,0.23,-0.01)))/vec3(0.78,0.78,0.78), col_required)*vec4(1,1,1,0.78), 0.01);
d = smin(d, mapCoralTree03((roty(0.01*PI)*(p-vec3(1.12,-0.04,0.06)))/vec3(0.85,0.85,0.85), col_required)*vec4(1,1,1,0.85), 0.01);
d = smin(d, mapCoralTree04((rotx(-0.06*PI)*roty(0.02*PI)*rotz(0.07*PI)*(p-vec3(1.05,-0.26,-0.22)))/vec3(0.77,0.77,0.77), col_required)*vec4(1,1,1,0.77), 0.01);
d = smin(d, mapCrab((rotx(0.06*PI)*roty(0.02*PI)*rotz(0.07*PI)*(p-vec3(-0.83,0.07,-0.42)))/vec3(0.35,0.35,0.37), col_required)*vec4(1,1,1,0.35), 0.01);

    return d;
}
