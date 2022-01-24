#include "common.glsl"

#define NO_MAP
#include "sdf_boat.glsl"
#include "sdf_coraltrees.glsl"
#include "sdf_crab.glsl"
#include "sdf_chip.glsl"
#include "sdf_fish.glsl"
#include "sdf_seafloor.glsl"

#ifndef __cplusplus
throw "Prevent this file from being compiled as GLSL because it freezes the IDE.";
#endif

vec4 map(vec3 p, bool col_required) {
    vec4 d = vec4(0,0,0, 1e3);

d = smin(d, mapBoat(((p-vec3(0.00,0.00,-0.97)))/vec3(1.00,1.00,1.00), col_required)*vec4(1,1,1,1.00), 0.01);
d = smin(d, mapChip((rotx(-0.45*PI)*roty(0.06*PI)*(p-vec3(0.00,0.04,0.01)))/vec3(0.16,0.16,0.16), col_required)*vec4(1,1,1,0.16), 0.01);
d = smin(d, mapCoralTree01((rotx(-1.00*PI)*roty(-0.02*PI)*rotz(-0.02*PI)*(p-vec3(1.30,-0.11,0.07)))/vec3(-0.75,-0.75,-0.81), col_required)*vec4(1,1,1,0.75), 0.01);
d = smin(d, mapCoralTree02((rotx(0.04*PI)*roty(0.02*PI)*rotz(-0.05*PI)*(p-vec3(1.19,0.23,0.02)))/vec3(0.78,0.78,0.78), col_required)*vec4(1,1,1,0.78), 0.01);
d = smin(d, mapCoralTree03((roty(0.01*PI)*(p-vec3(1.12,-0.04,0.08)))/vec3(0.85,0.85,0.85), col_required)*vec4(1,1,1,0.85), 0.01);
d = smin(d, mapCoralTree04((rotx(-0.06*PI)*roty(0.02*PI)*rotz(0.07*PI)*(p-vec3(1.05,-0.26,-0.19)))/vec3(0.77,0.77,0.77), col_required)*vec4(1,1,1,0.77), 0.01);
d = smin(d, mapCrab((rotx(0.06*PI)*roty(0.02*PI)*rotz(0.07*PI)*(p-vec3(-0.87,0.06,-0.41)))/vec3(0.37,0.37,0.39), col_required)*vec4(1,1,1,0.37), 0.01);
d = smin(d, mapFish((rotx(0.01*PI)*roty(0.03*PI)*rotz(-0.05*PI)*(p-vec3(1.29,-0.77,-0.94)))/vec3(0.35,0.32,0.32), col_required)*vec4(1,1,1,0.32), 0.01);
d = smin(d, mapFish((rotx(0.01*PI)*roty(-0.05*PI)*rotz(0.14*PI)*(p-vec3(-0.99,0.75,0.55)))/vec3(0.31,0.28,0.28), col_required)*vec4(1,1,1,0.28), 0.01);
d = smin(d, mapFish((rotx(0.01*PI)*roty(0.06*PI)*rotz(0.09*PI)*(p-vec3(-1.62,0.79,0.23)))/vec3(0.35,0.32,0.32), col_required)*vec4(1,1,1,0.32), 0.01);
d = smin(d, mapRock((rotx(0.06*PI)*roty(0.09*PI)*rotz(-0.02*PI)*(p-vec3(1.11,-0.01,-0.80)))/vec3(0.43,0.43,0.43), col_required)*vec4(1,1,1,0.43), 0.01);
d = smin(d, mapRock((rotx(0.97*PI)*roty(-0.06*PI)*rotz(-0.08*PI)*(p-vec3(-0.85,0.22,-0.98)))/vec3(0.49,0.49,0.49), col_required)*vec4(1,1,1,0.49), 0.01);
d = smin(d, mapRock((rotx(0.96*PI)*roty(0.03*PI)*rotz(0.11*PI)*(p-vec3(-0.03,-0.65,-1.01)))/vec3(0.39,0.39,0.31), col_required)*vec4(1,1,1,0.31), 0.01);
d = smin(d, mapSeashell((rotx(-0.06*PI)*roty(0.05*PI)*rotz(0.01*PI)*(p-vec3(1.01,-0.33,-0.76)))/vec3(0.20,0.16,0.20), col_required)*vec4(1,1,1,0.16), 0.01);
d = smin(d, mapSeashell((rotx(-0.05*PI)*roty(0.03*PI)*rotz(0.35*PI)*(p-vec3(0.90,-0.04,-0.74)))/vec3(0.17,0.15,0.17), col_required)*vec4(1,1,1,0.15), 0.01);
d = smin(d, mapSeashell((rotx(-0.16*PI)*roty(-0.08*PI)*rotz(0.30*PI)*(p-vec3(0.86,-0.22,-0.77)))/vec3(0.17,0.15,0.19), col_required)*vec4(1,1,1,0.15), 0.01);
d = smin(d, mapSeashell((rotx(-0.01*PI)*roty(0.01*PI)*rotz(0.12*PI)*(p-vec3(-0.24,-0.33,-0.93)))/vec3(0.20,0.17,0.20), col_required)*vec4(1,1,1,0.17), 0.01);
d = smin(d, mapSeashell((rotx(0.00*PI)*roty(0.02*PI)*rotz(0.24*PI)*(p-vec3(0.18,-0.45,-0.91)))/vec3(0.18,0.15,0.18), col_required)*vec4(1,1,1,0.15), 0.01);
d = smin(d, mapSeashell((rotx(-0.03*PI)*roty(0.01*PI)*rotz(0.18*PI)*(p-vec3(0.14,0.40,-0.91)))/vec3(0.19,0.17,0.19), col_required)*vec4(1,1,1,0.17), 0.01);
d = smin(d, mapSeashell((rotx(0.02*PI)*roty(-0.04*PI)*rotz(-0.09*PI)*(p-vec3(-0.58,-0.62,-0.91)))/vec3(0.22,0.18,0.22), col_required)*vec4(1,1,1,0.18), 0.01);
d = smin(d, mapStarfish01((rotx(0.06*PI)*roty(0.02*PI)*rotz(-0.32*PI)*(p-vec3(0.61,-0.26,-0.85)))/vec3(0.19,0.19,0.19), col_required)*vec4(1,1,1,0.19), 0.01);
d = smin(d, mapStarfish02((rotx(0.02*PI)*roty(0.07*PI)*rotz(-0.06*PI)*(p-vec3(0.46,-0.59,-0.88)))/vec3(0.18,0.18,0.18), col_required)*vec4(1,1,1,0.18), 0.01);
d = smin(d, mapStarfish02((rotx(-0.02*PI)*roty(-0.02*PI)*rotz(-0.25*PI)*(p-vec3(-0.25,0.63,-0.94)))/vec3(0.19,0.19,0.19), col_required)*vec4(1,1,1,0.19), 0.01);
d = smin(d, mapStarfish03((rotx(0.06*PI)*roty(0.06*PI)*rotz(-0.26*PI)*(p-vec3(0.64,0.37,-0.84)))/vec3(0.19,0.19,0.19), col_required)*vec4(1,1,1,0.19), 0.01);

    return d;
}
