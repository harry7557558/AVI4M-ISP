#include "common.glsl"


vec4 mapStarfish01(vec3 p, bool col_required) {
    vec4 bound = vec4(1,0,0,sdEllipsoid(p-vec3(0,0,0.1),vec3(1.5,1.5,0.7))), boundw = vec4(0,0,0,0.4);
    if (bound.w > 0.0) return bound+boundw;  // clipping
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q = vec3(r*cossin(asin(0.7*sin(2.5*a))/2.5), p.z-0.1*exp(-sqr(2.0*r)));
    q.y = length(vec2(q.y,0.01));
    float noise = exp(2.0*SimplexNoise3D(10.0*vec3(r*cossin(asin(0.9*sin(2.5*a))/2.5),p.z)));
    vec4 mid = vec4(vec3(0.6,0.2,0.1)*exp(-0.5*r)+0.1*noise-0.1,
        sdSegment(q-vec3(0,0,0.1*exp(-sqr(10.0*q.y))),vec3(0.0),vec3(1,0,0))-0.2*exp(-0.8*r*r)-0.005*noise);
    vec4 disk = vec4(0.7,0.3,0.15, sdEllipsoid(q-vec3(0,0,0.0),vec3(0.4,0.4,0.2)));
    vec4 d = smin(mid, disk, 0.4) - vec4(0,0,0,0.01);
    if (col_required) d.xyz = pow(d.xyz, vec3(0.9));
    return d;
}

vec4 mapStarfish02(vec3 p, bool col_required) {
    vec4 bound = vec4(1,0,0,sdEllipsoid(p-vec3(0,0,0.1),vec3(1.5,1.5,0.7))), boundw = vec4(0,0,0,0.4);
    if (bound.w > 0.0) return bound+boundw;  // clipping
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q = vec3(r*cossin(asin(0.7*sin(2.5*a))/2.5), p.z-0.1*exp(-sqr(2.0*r)));
    q.y = length(vec2(q.y,0.01));
    float noise = exp(2.0*SimplexNoise3D(10.0*vec3(r*cossin(asin(0.98*sin(2.5*a))/2.5),p.z)));
    vec4 mid = vec4(vec3(0.75,0.45,0.2)*exp(-0.5*r)*exp(1.0*q.y)+0.1*noise-0.1,
        sdSegment(q-vec3(0,0,0.1*exp(-sqr(10.0*q.y))-0.05*r*r),vec3(0.0),vec3(1,0,0))-0.2*exp(-0.8*r*r)-0.008*noise);
    vec4 disk = vec4(0.85,0.45,0.2, sdEllipsoid(q-vec3(0,0,0.0),vec3(0.4,0.4,0.2)));
    vec4 d = smin(mid, disk, 0.4) - vec4(0,0,0,0.00);
    if (col_required) d.xyz = pow(d.xyz, vec3(1.1));
    return d;
}

vec4 mapStarfish03(vec3 p, bool col_required) {
    vec4 bound = vec4(1,0,0,sdEllipsoid(p-vec3(0,0,0.1),vec3(1.5,1.5,0.7))), boundw = vec4(0,0,0,0.4);
    if (bound.w > 0.0) return bound+boundw;  // clipping
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q = vec3(r*cossin(asin(0.7*sin(2.5*a))/2.5), p.z-0.1*exp(-sqr(2.0*r)));
    q.y = length(vec2(q.y,0.01));
    float noise = exp(2.5*SimplexNoise3D(8.0*vec3(r*cossin(asin(0.98*sin(2.5*a))/2.5),p.z)));
    vec4 mid = vec4(vec3(0.95,0.8,0.5)*exp(0.2*r)*exp(-1.0*q.y)-0.15*noise+0.1,
        sdSegment(q-vec3(0,0,0.1*exp(-sqr(10.0*q.y))-0.05*r*r),vec3(0.0),vec3(1,0,0))-0.2*exp(-0.8*r*r)-0.005*noise);
    vec4 disk = vec4(0.85,0.7,0.45, sdEllipsoid(q-vec3(0,0,0.0),vec3(0.4,0.4,0.2)));
    vec4 d = smin(mid, disk, 0.3) - vec4(0,0,0,0.00);
    if (col_required) d.xyz = 0.8*pow(d.xyz, vec3(0.9));
    return d;
}


vec4 mapSeashell(vec3 p, bool col_required) {
    vec4 bound = vec4(1,0,0,length(p)-0.6), boundw = vec4(0,0,0,0.25);
    if (bound.w > 0.0) return bound+boundw;  // clipping
    vec4 body = vec4(0.95,0.85,0.55, sdEllipsoid(p, vec3(0.4,0.3,0.25)));
    vec4 opening = vec4(0.95,0.9,0.6, sdTorus(p-vec3(0,0,-0.1),0.3+0.05*cos(2.0*atan(p.y,p.x)),0.05));
    vec4 shell = smin(body, opening, 0.2);
    shell.w += 0.05*exp(-sqr(10.0*p.y))*exp(-sqr(4.0*p.x));
    if (col_required) shell.xyz = pow(shell.xyz, vec3(1.6));
    return shell;
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    // return mapStarfish03(p, col_required);
    return mapSeashell(p, col_required);
}
#endif
