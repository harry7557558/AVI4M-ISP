#include "common.glsl"


vec4 mapBabyHead(vec3 p, bool col_required) {
    p.x = -p.x;
    vec3 q = roty(0.1)*p-vec3(0,0,-0.1);
    q.xy /= 0.9+0.05*q.z*q.z;
    q.z = q.z-0.1 - 0.47*clamp(q.z,0.0,2.0)*(1.0+exp(-1.2*length(vec2(p.y,0.2))))+0.1*exp(-sqr(4.0*p.x));
    q.x = 1.2*q.x + 0.1+0.1*q.x*cos(10.0*q.z) - 0.2*tanh(10.0*p.x)*abs(p.z);
    vec4 beak = vec4(mix(vec3(0.8,0.5,0.1),vec3(0.9,0.8,0.5),smoothstep(0.,1.,0.5+4.0*(length(q.xy)-0.5))), sdTorus(q, 0.5, 0.17));
    q = p + vec3(0.03*cos(8.0*p.y), 0,0);
    vec4 head = smax(
        vec4(0.85,0.55,0.4, sdEllipsoid(q, vec3(0.75*exp(0.05*q.x)*exp(-0.1*q.z),0.7,0.85))),
        vec4(0.55,0.05,0.05, -(length(q-vec3(-0.1,0,1))-0.5)),
        0.2);
    if (col_required) {
        float t = 4.0*q.x+sin(4.0*q.z);
        head.xyz = mix(head.xyz, vec3(0.5,0.45,0.4), smoothstep(0.,1.,t));
        head.xyz = mix(head.xyz, vec3(0.35,0.3,0.25), smoothstep(0.,1.,0.5*t-1.0));
    }
    head = smin(head, beak, 0.3);
    q = vec3(p.x, abs(p.y), p.z);
    vec4 eye = vec4(0.2,0,0.1, length(q-vec3(0.3,0.58,0.2))-0.1);
    //head = smax(head, vec4(0.8,0.9,0.98,-eye.w), 0.05);
    head = smin(head, eye, 0.05);
    return head;
}

vec4 mapBaby(vec3 p, bool col_required) {
    float bound = sdSegment(p, vec3(0,0,-1.3), vec3(0,0,1.3)) - 1.1, boundw = 0.4;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    p = roty(0.05*PI)*p / 0.7;
    p.z += 1.1;
    vec3 q;
    q = roty(-0.05*PI)*(p-vec3(-0.2,0,0));
    vec4 body = vec4(vec3(0.7,0.45,0.35)-vec3(0.15)*q.x, sdEllipsoid(q, vec3(vec2(0.8*(0.95+0.25*tanh(-p.z))),1.0)));
    q = roty(-0.3*PI)*(p-vec3(-0.7,0,-0.7));
    vec4 tail = vec4(0.6,0.4,0.25, sdEllipsoid(q, vec3(0.3,0.3,0.5)));
    body = smin(body, tail, 0.5);
    q = p;
    vec4 neck = vec4(0.75,0.55,0.2, sdSegment(q, vec3(0,0,0), vec3(0,0,1.5))-0.3);
    body = smin(body, neck, 0.5);
    q = roty(-0.2*PI)*(p-vec3(-0.1,0,2.1));
    vec4 head = mapBabyHead(q/0.9, col_required)*0.9;
    body = smin(body, head, 0.2*(1.0-tanh(3.0*p.x)));
    q = vec3(p.x, abs(p.y), p.z);
    vec4 thigh = vec4(0.7,0.45,0.25, sdSegment(q, vec3(-0.3,0.45,-0.7), vec3(0.1,0.7,-0.4))-0.15);
    vec4 shank = vec4(0.7,0.45,0.25, sdSegment(q, vec3(0.1,0.7,-0.4), vec3(0.0,0.8,-1.1))-(0.12+0.05*tanh(4.0*(q.z+0.8))));
    vec4 feet = vec4(0.55,0.35,0.3, sdSegment(q, vec3(0.0,0.8,-1.1), vec3(0.3,0.8,-1.3))-0.08);
    body = smin(body, smin(smin(thigh, shank, 0.2), feet, 0.2), 0.2);
    return body * vec4(1,1,1,0.7);
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    //return mapBabyHead(p, col_required);
    return mapBaby(p, col_required);
}
#endif
