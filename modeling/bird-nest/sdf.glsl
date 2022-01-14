#include "common.glsl"

#define NO_MAP
#include "sdf_bird.glsl"
#include "sdf_nest.glsl"
#include "sdf_baby.glsl"
#include "sdf_berries.glsl"
#include "sdf_liriope.glsl"

#ifndef __cplusplus
throw "Prevent this file from being compiled as GLSL because it freezes the IDE.";
#endif

vec4 map(vec3 p, bool col_required) {
    vec4 d = vec4(0,0,0, 1e3);

d = smin(d, mapBaby((rotx(-0.01*PI)*roty(-0.10*PI)*rotz(-0.06*PI)*(p-vec3(0.11,-0.18,-0.29)))/vec3(0.35,0.35,0.37), col_required)*vec4(1,1,1,0.35), 0.01);
d = smin(d, mapBaby((rotx(0.02*PI)*roty(-0.14*PI)*rotz(0.07*PI)*(p-vec3(0.25,0.19,-0.34)))/vec3(0.32,0.32,0.34), col_required)*vec4(1,1,1,0.32), 0.01);
d = smin(d, mapBaby((rotx(0.01*PI)*roty(-0.17*PI)*rotz(0.04*PI)*(p-vec3(0.47,-0.04,-0.54)))/vec3(0.32,0.32,0.34), col_required)*vec4(1,1,1,0.32), 0.01);
d = smin(d, mapBerriesFlowerT((rotx(-0.23*PI)*roty(-0.10*PI)*rotz(-0.05*PI)*(p-vec3(0.32,-0.95,-0.80)))/vec3(0.19,0.19,0.20), col_required)*vec4(1,1,1,0.19), 0.01);
d = smin(d, mapBerriesFruitT((rotx(-0.11*PI)*roty(0.05*PI)*rotz(-0.03*PI)*(p-vec3(0.16,-0.79,-0.82)))/vec3(0.24,0.24,0.24), col_required)*vec4(1,1,1,0.24), 0.01);
d = smin(d, mapBerriesFruitT((rotx(-0.01*PI)*roty(-0.09*PI)*rotz(0.08*PI)*(p-vec3(0.38,-0.77,-0.82)))/vec3(0.19,0.19,0.19), col_required)*vec4(1,1,1,0.19), 0.01);
d = smin(d, mapBerriesFruitT((rotx(0.08*PI)*roty(-0.13*PI)*rotz(0.19*PI)*(p-vec3(0.96,-0.25,-0.95)))/vec3(0.20,0.20,0.20), col_required)*vec4(1,1,1,0.20), 0.01);
d = smin(d, mapBerriesFruitT((rotx(0.00*PI)*roty(-0.02*PI)*rotz(-0.01*PI)*(p-vec3(0.83,-0.44,-0.96)))/vec3(0.17,0.19,0.17), col_required)*vec4(1,1,1,0.17), 0.01);
d = smin(d, mapBerriesFruitT((rotx(0.14*PI)*roty(-0.18*PI)*rotz(0.38*PI)*(p-vec3(1.00,-0.43,-1.00)))/vec3(0.16,0.16,0.16), col_required)*vec4(1,1,1,0.16), 0.01);
d = smin(d, mapBerriesLeafT((rotx(-0.01*PI)*roty(-0.11*PI)*rotz(0.13*PI)*(p-vec3(0.25,-0.76,-0.82)))/vec3(0.31,0.31,0.31), col_required)*vec4(1,1,1,0.31), 0.01);
d = smin(d, mapBerriesLeafT((rotx(0.13*PI)*roty(-0.09*PI)*rotz(0.61*PI)*(p-vec3(0.15,-0.70,-0.84)))/vec3(0.34,0.27,0.27), col_required)*vec4(1,1,1,0.27), 0.01);
d = smin(d, mapBerriesLeafT((rotx(0.03*PI)*roty(-0.10*PI)*rotz(0.41*PI)*(p-vec3(0.17,-0.81,-0.83)))/vec3(0.25,0.25,0.25), col_required)*vec4(1,1,1,0.25), 0.01);
d = smin(d, mapBerriesLeafT((rotx(0.11*PI)*roty(-0.06*PI)*rotz(-0.08*PI)*(p-vec3(0.36,-0.72,-0.84)))/vec3(0.22,0.22,0.22), col_required)*vec4(1,1,1,0.22), 0.01);
d = smin(d, mapBerriesLeafT((rotx(-0.04*PI)*roty(-0.10*PI)*rotz(0.49*PI)*(p-vec3(0.86,-0.42,-0.96)))/vec3(0.26,0.26,0.26), col_required)*vec4(1,1,1,0.26), 0.01);
d = smin(d, mapBerriesLeafT((rotx(0.01*PI)*roty(-0.17*PI)*rotz(0.05*PI)*(p-vec3(0.93,-0.35,-0.99)))/vec3(0.32,0.32,0.32), col_required)*vec4(1,1,1,0.32), 0.01);
d = smin(d, mapBerriesLeafT((rotx(0.12*PI)*roty(-0.13*PI)*rotz(0.35*PI)*(p-vec3(0.91,-0.37,-1.01)))/vec3(0.29,0.29,0.29), col_required)*vec4(1,1,1,0.29), 0.01);
d = smin(d, mapBerriesLeafT((rotx(-0.12*PI)*roty(-0.10*PI)*rotz(-0.31*PI)*(p-vec3(0.96,-0.32,-0.98)))/vec3(0.26,0.26,0.26), col_required)*vec4(1,1,1,0.26), 0.01);
d = smin(d, mapBerriesLeafT((rotx(0.03*PI)*roty(-0.15*PI)*rotz(0.31*PI)*(p-vec3(0.28,-0.76,-0.85)))/vec3(0.30,0.26,0.26), col_required)*vec4(1,1,1,0.26), 0.01);
d = smin(d, mapBird((roty(0.04*PI)*(p-vec3(1.55,0.00,0.36)))/vec3(0.52,0.52,0.52), col_required)*vec4(1,1,1,0.52), 0.01);
d = smin(d, mapEggBrokenBE((rotx(-0.24*PI)*roty(0.11*PI)*rotz(-0.42*PI)*(p-vec3(-0.34,-0.10,-0.44)))/vec3(0.30,0.30,0.30), col_required)*vec4(1,1,1,0.30), 0.01);
d = smin(d, mapEggBrokenLE((rotx(0.85*PI)*roty(-0.81*PI)*rotz(0.53*PI)*(p-vec3(-0.32,0.10,-0.39)))/vec3(0.30,0.30,0.30), col_required)*vec4(1,1,1,0.30), 0.01);
d = smin(d, mapEggBrokenS1((rotx(-0.00*PI)*roty(-0.78*PI)*rotz(1.00*PI)*(p-vec3(-0.26,0.01,-0.78)))/vec3(0.30,0.30,0.30), col_required)*vec4(1,1,1,0.30), 0.01);
d = smin(d, mapEggBrokenS2((rotx(-0.41*PI)*roty(0.17*PI)*rotz(-0.43*PI)*(p-vec3(-0.04,0.29,-0.54)))/vec3(0.30,0.30,0.30), col_required)*vec4(1,1,1,0.30), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(-0.09*PI)*roty(0.05*PI)*rotz(-0.02*PI)*(p-vec3(-0.91,-0.54,-0.38)))/vec3(0.45,0.45,0.45), col_required)*vec4(1,1,1,0.45), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(-0.07*PI)*roty(0.00*PI)*rotz(-0.06*PI)*(p-vec3(-0.71,-0.63,-0.51)))/vec3(0.45,0.45,0.45), col_required)*vec4(1,1,1,0.45), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(0.07*PI)*roty(-0.07*PI)*rotz(0.01*PI)*(p-vec3(0.47,0.69,-1.00)))/vec3(0.50,0.50,0.50), col_required)*vec4(1,1,1,0.50), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(0.07*PI)*roty(0.11*PI)*rotz(0.17*PI)*(p-vec3(0.27,0.81,-0.70)))/vec3(0.50,0.50,0.50), col_required)*vec4(1,1,1,0.50), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(0.03*PI)*roty(0.20*PI)*rotz(0.40*PI)*(p-vec3(0.34,0.94,-0.79)))/vec3(0.40,0.40,0.40), col_required)*vec4(1,1,1,0.40), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(-0.07*PI)*roty(0.07*PI)*rotz(-0.13*PI)*(p-vec3(-0.12,-0.84,-0.72)))/vec3(0.38,0.38,0.38), col_required)*vec4(1,1,1,0.38), 0.01);
d = smin(d, mapLiriopeFlowersT((rotx(-0.19*PI)*roty(0.09*PI)*rotz(-0.13*PI)*(p-vec3(-0.45,-0.99,-0.79)))/vec3(0.38,0.38,0.38), col_required)*vec4(1,1,1,0.38), 0.01);
d = smin(d, mapLiriopeLeaves3((rotx(0.12*PI)*roty(0.09*PI)*rotz(0.18*PI)*(p-vec3(0.39,0.91,-0.85)))/vec3(0.36,0.36,0.67), col_required)*vec4(1,1,1,0.36), 0.01);
d = smin(d, mapLiriopeLeaves3((rotx(-0.15*PI)*roty(-0.08*PI)*rotz(0.10*PI)*(p-vec3(-0.53,-0.73,-0.63)))/vec3(0.37,0.37,0.51), col_required)*vec4(1,1,1,0.37), 0.01);
d = smin(d, mapLiriopeLeaves3((rotx(-0.07*PI)*roty(-0.12*PI)*rotz(0.50*PI)*(p-vec3(-0.42,-0.73,-0.72)))/vec3(0.28,0.28,0.40), col_required)*vec4(1,1,1,0.28), 0.01);
d = smin(d, mapLiriopeLeaves4((rotx(-0.13*PI)*roty(0.06*PI)*rotz(-0.03*PI)*(p-vec3(-0.67,-0.60,-0.59)))/vec3(0.53,0.53,0.53), col_required)*vec4(1,1,1,0.53), 0.01);
d = smin(d, mapLiriopeLeaves4((rotx(-0.03*PI)*roty(0.04*PI)*rotz(0.05*PI)*(p-vec3(-1.08,-0.41,-0.41)))/vec3(0.39,0.39,0.46), col_required)*vec4(1,1,1,0.39), 0.01);
d = smin(d, mapLiriopeLeaves4((rotx(0.08*PI)*roty(0.19*PI)*rotz(0.28*PI)*(p-vec3(0.29,0.79,-0.76)))/vec3(0.52,0.52,0.52), col_required)*vec4(1,1,1,0.52), 0.01);
d = smin(d, mapLiriopeLeaves4((rotx(-0.13*PI)*roty(0.13*PI)*rotz(-0.10*PI)*(p-vec3(-0.09,-0.80,-0.76)))/vec3(0.46,0.46,0.46), col_required)*vec4(1,1,1,0.46), 0.01);
d = smin(d, mapLiriopeLeaves5((rotx(-0.14*PI)*roty(0.11*PI)*rotz(-0.08*PI)*(p-vec3(-0.86,-0.43,-0.48)))/vec3(0.55,0.55,0.55), col_required)*vec4(1,1,1,0.55), 0.01);
d = smin(d, mapLiriopeLeaves5((rotx(-0.07*PI)*roty(-0.14*PI)*rotz(-0.51*PI)*(p-vec3(0.49,0.65,-0.98)))/vec3(0.40,0.40,0.55), col_required)*vec4(1,1,1,0.40), 0.01);
d = smin(d, mapLiriopeLeaves5((rotx(-0.17*PI)*roty(0.03*PI)*rotz(0.06*PI)*(p-vec3(-0.44,-0.91,-0.85)))/vec3(0.43,0.43,0.43), col_required)*vec4(1,1,1,0.43), 0.01);
d = smin(d, mapNest((roty(0.82*PI)*(p-vec3(0.04,0.00,-0.40)))/vec3(-0.54,-0.45,-0.54), col_required)*vec4(1,1,1,0.45), 0.01);

    return d;
}
