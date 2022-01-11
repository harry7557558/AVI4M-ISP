#include "common.glsl"


vec4 mapLiriopeFlowersLayer(float elev, vec3 p, bool col_required) {
    float r = length(p.xy), a = atan(p.y, p.x);
    a += 0.05*sin(4.0*a) + 0.05*cos(a);
    vec3 q = vec3(r*cossin(asin(0.97*sin(2.5*a))/2.5), p.z);
    q = roty(elev)*q;
    vec4 stem = vec4(0.9,0.85,0.95, sdSegment(q, vec3(0.0), vec3(1.0,0,0))-0.1);
    float r1 = length(q.yz), a1 = atan(q.z, q.y);
    vec3 q1 = vec3(q.x, r1*cossin(asin(0.93*cos(2.5*a1))/2.5));
    q1 = roty(0.3*PI)*(q1-vec3(1,0,0))+vec3(1,0,0);
    vec4 flower = vec4(mix(vec3(0.65,0.45,0.7), vec3(0.9,0.8,0.95), (3.0*(q1.x-1.0)+0.5)+0.15*sin(5.0*a1)),
        sdEllipsoid(q1-vec3(1.0,0.0,0.05), vec3(0.45,vec2(0.3+0.1*(q1.x-1.0)))));
    vec4 d = smin(stem, flower, 0.1);
    return d;
}

vec4 mapLiriopeFlowers(vec3 p, bool col_required) {
    p /= 0.8;
    p.x += 0.1*cos(p.z);
    float bound = sdSegment(p,vec3(0,0,-2.8),vec3(0,0,2.0))-1.0, boundw = 0.4;
    if (bound > 0.0) return vec4(1,0,0, (bound+boundw)*0.8);  // clipping
    vec4 stem = vec4(0,0,0, sdSegment(p, vec3(0,0,-3.0), vec3(0,0,2))-0.06);
    if (col_required) {
        stem.xyz = mix(vec3(0.6,0.35,0.65), vec3(0.9,0.8,0.95), smootherstep((p.z+1.0)/3.0));
        stem.xyz *= mix(vec3(0.35,0.2,0.25), vec3(1.0), 0.2+0.8*smootherstep((p.z+2.5)/2.5));
    }
    float seed = 0.0;
    for (float t=0.0; t<1.0; t+=1.0/11.0) {
        float h = mix(-1.5, 2.0, 1.0-pow(1.0-t,1.2));
        h += 0.5*t*(1.0-t)*(2.0*rand(seed)-1.0);
        vec3 q = p-vec3(0,0,h);
        q.xy = rot2(2.0*PI*rand(seed))*q.xy;
        float elev = 0.05*PI + 0.15*PI*t + mix(0.05*PI, 0.2*PI, rand(seed));
        float sc = 0.08*(1.0+0.3*t-1.3*t*t) + mix(0.25, 0.4, smoothstep(0.,1.,(p.x+1.0)/2.0)) * (1.0+0.1*rand(seed));
        vec4 d = mapLiriopeFlowersLayer(elev, q/sc, col_required)*vec4(1,1,1,sc);
        if (col_required) d.xyz = mix(d.xyz, mix(vec3(0.65,0.45,0.75), vec3(0.8,0.7,0.9), t), 0.2);
        stem = smin(stem, d, 0.05);
    }
    return stem * vec4(1,1,1, 0.8);
}
vec4 mapLiriopeFlowersT(vec3 p, bool col_required) {
    return mapLiriopeFlowers( p/0.5-vec3(0,0,2.5), col_required) * vec4(1,1,1,0.5);
}

vec4 mapLiriopeLeaf(vec3 p, bool col_required) {
    p.z -= 0.5*cos(0.25*PI*p.x);
    float bound = sdSegment(p,vec3(-2.0,0,0),vec3(2.0,0,0))-0.5, boundw = 0.3;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    float near_stem = 1.0-0.3/(1.0+pow(abs(0.6*(p.x+2.0)),4.0));
    float w = pow(max(1.0-sqr(p.x/2.0), 0.0), max(0.3+0.1*p.x,0.05)) * (0.25/(1.0+sqr(0.3*(p.x-0.5)))) * near_stem;
    float u = clamp(p.y/max(w,1e-4), -1.0, 1.0);
    float thickness = 0.2 * 0.2*pow(w/0.2,0.4) * pow(max(1.0-u*u,0.0),0.5) * (1.0+sqr(p.y/0.3)) * (exp(-(-0.08)*(p.x+2.0))) / pow(near_stem, 2.0);
    float veins = 0.03*cos(15.0*u)*(1.0-u*u);
    float zd = 0.05/sqrt(1.0+sqr(p.y/0.1));
    vec4 leaf = vec4(0,0,0, (w==0. ? length(p.yz) : sdSegment(p.yz+vec2(0,zd), vec2(-w,0), vec2(w,0))) - thickness * (1.0+veins));
    if (col_required) {
        leaf.xyz = pow(mix(vec3(0.35,0.55,0.25), vec3(0.65,0.8,0.5), (p.x+2.0)/4.0), vec3(1.8));
        leaf.xyz = mix(leaf.xyz, vec3(0.45,0.65,0.3), 1.0-20.0*zd) * vec3(1.0+1.0*veins);
    }
    leaf.w += 0.01;
    return leaf;
}
vec4 mapLiriopeLeafT(vec3 p, bool col_required) {
    return mapLiriopeLeaf( p/0.5-vec3(2.0,0,0), col_required) * vec4(1,1,1,0.5);
}
vec4 mapLiriopeLeaves3(vec3 p, bool col_required) {
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q = vec3(r*cossin(asin(0.95*sin(1.5*a))/1.5), p.z);
    return mapLiriopeLeafT(roty(0.2*PI+0.1*p.y)*q, col_required);
}
vec4 mapLiriopeLeaves4(vec3 p, bool col_required) {
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q = vec3(r*cossin(asin(0.95*sin(2.0*a))/2.0), p.z);
    return mapLiriopeLeafT(roty(0.2*PI+0.1*p.y)*q, col_required);
}
vec4 mapLiriopeLeaves5(vec3 p, bool col_required) {
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q = vec3(r*cossin(asin(0.95*sin(2.5*a))/2.5), p.z);
    return mapLiriopeLeafT(roty(0.2*PI+0.1*p.y)*q, col_required);
}


// test Blender export
vec4 mapLiriopeGroupTest01(vec3 p, bool col_required) {
    vec4 d = vec4(1,1,1, 100.);
d = smin(d, mapLiriopeFlowers(((p-vec3(-0.19,0.24,0.59)))/vec3(0.68,0.68,0.68), col_required)*vec4(1,1,1,0.68), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.06*PI)*roty(0.20*PI)*rotz(0.01*PI)*(p-vec3(0.74,0.34,-0.29)))/vec3(0.54,0.48,0.61), col_required)*vec4(1,1,1,0.48), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.06*PI)*roty(0.20*PI)*rotz(-0.18*PI)*(p-vec3(0.70,1.00,-0.14)))/vec3(0.79,0.68,0.79), col_required)*vec4(1,1,1,0.68), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.07*PI)*roty(0.27*PI)*rotz(0.95*PI)*(p-vec3(-1.18,0.15,-0.18)))/vec3(0.67,0.45,0.83), col_required)*vec4(1,1,1,0.45), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.11*PI)*roty(0.14*PI)*rotz(0.32*PI)*(p-vec3(0.54,-0.90,-0.50)))/vec3(0.78,0.69,0.88), col_required)*vec4(1,1,1,0.69), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.09*PI)*roty(0.24*PI)*rotz(-0.83*PI)*(p-vec3(-1.15,0.94,0.06)))/vec3(0.72,0.58,0.89), col_required)*vec4(1,1,1,0.58), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.06*PI)*roty(0.20*PI)*rotz(0.47*PI)*(p-vec3(-0.17,-0.65,-0.19)))/vec3(0.56,0.47,0.65), col_required)*vec4(1,1,1,0.47), 0.01);
d = smin(d, mapLiriopeLeaf((rotx(0.06*PI)*roty(0.20*PI)*rotz(0.66*PI)*(p-vec3(-0.80,-0.64,-0.20)))/vec3(0.79,0.45,0.95), col_required)*vec4(1,1,1,0.45), 0.01);
    return d;
}


// map
vec4 map(vec3 p, bool col_required) {
    //return mapLiriopeFlowersLayer(0.2*PI, p, col_required);
    //return mapLiriopeFlowers(p, col_required);
    //return mapLiriopeLeaf(p, col_required);
    return mapLiriopeLeaves5(p, col_required);
    return mapLiriopeGroupTest01(p, col_required);
}
