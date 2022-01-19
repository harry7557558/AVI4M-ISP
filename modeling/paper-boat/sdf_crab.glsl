#include "common.glsl"


vec4 mapCrabShell(vec3 p, bool col_required) {
    const float sc = 0.5;
    p = p * vec3(1,-1,1) / sc - vec3(1.0, -1.0, 0.9);
    const float b = 0.15;
    float r = length(p.xy);
    float a = mix(0.0, -0.4, smoothstep(0.0, 1.0, 0.5*(r-0.6))) - 0.2*r;
    p.xy = mat2(cos(a),-sin(a),sin(a),cos(a))*p.xy;
    float t = atan(p.y, p.x);
    float n = min( (log(r)/b-t)/(2.0*PI), 0.0);
    float n0 = floor(n), n1 = ceil(n);
    float x0 = exp(b*(t+2.0*PI*n0));
    float x1 = exp(b*(t+2.0*PI*n1));
    float r0 = 1.0*x0;
    float r1 = 1.0*x1;
    float h0 = p.z + 1.5*(x0-1.0);
    float h1 = p.z + 1.5*(x1-1.0);
    float d0 = length(vec2(x0-r,h0)) - r0;
    float d1 = length(vec2(x1-r,h1)) - r1;
    float d, dx, dy;
    if (d0 < 0.0) d = d0, dx = x0-r, dy = h0;
    else if (d1 < 0.0 && d1<-d0) d = -d0, dx = x0-r, dy = h0;
    else if (d1 < 0.0) d = d1, dx = x1-r, dy = h1;
    else if (d0 < d1) d = d0, dx = x0-r, dy = h0;
    else d = d1, dx = x1-r, dy = h1;
    d += 0.002*r*sin(40.*t);
    d += 0.002*r*sin(40.*atan(dy,dx));
    d = abs(d)-0.0*r-0.25;
    vec4 shell = vec4(1, 0, 0, d);
    vec4 opening = vec4(1, 0, 0, sdTorus((p-vec3(-1.65,0,-0.91)).xzy*vec3(1,0.96,1), 1.55, 0.04*(0.0*r+8.0)));
    shell = smin(shell, opening, 0.1) + vec4(0,0,0,0.16);
    return shell * vec4(1,1,1,sc);
}


vec4 mapCrabLeg(vec3 p, vec3 joint1, vec3 joint2, vec3 joint3, vec3 tip, bool col_required) {
    vec3 p0, p1, dp; float h;
    vec4 seg0 = vec4(1,0,0, length(p-vec3(0,0.05,0))-(0.09+0.1*p.y));
    p0 = mix(vec3(0.), joint1, 0.1), p1 = mix(vec3(0.), joint1, 0.9);
    float seg1_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg1_d -= 0.12-0.05*cos(normalize(dp).z)+0.1*h*(1.0-h);
    vec4 seg1 = vec4(1,0,0, seg1_d);
    p0 = mix(joint1, joint2, 0.1), p1 = mix(joint1, joint2, 0.9);
    float seg2_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg2_d -= 0.12-0.05*cos(normalize(dp).x)-0.05*dp.x+0.1*h*(1.0-h);
    vec4 seg2 = vec4(1,0,0, seg2_d);
    p0 = mix(joint2, joint3, 0.1), p1 = mix(joint2, joint3, 0.9);
    float seg3_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg3_d -= 0.12-0.05*cos(normalize(dp).x)-0.05*dp.x+0.1*h*(1.0-h);
    vec4 seg3 = vec4(1,0,0, seg3_d);
    p0 = mix(joint3, tip, 0.1), p1 = mix(joint3, tip, 0.9);
    float seg4_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg4_d -= (0.12-0.05*cos(normalize(dp).x)-0.05*dp.x+0.1*h*(1.0-h))*exp(-0.8*h);
    vec4 seg4 = vec4(1,0,0, seg4_d);
    return smin(smin(smin(seg0, seg1, 0.05), smin(seg2, seg3, 0.01), 0.01), seg4, 0.01) - vec4(0,0,0,0.02);
}


vec4 mapCrabBody(vec3 p, bool col_required) {
    vec3 q = roty(0.1*PI)*p; q.y=0.2*log(2.0*cosh(q.y/0.2)), q.z-=0.1;
    float body_d = sdLnNormEllipsoid(q, vec3(0.8,0.5,0.38)+0.02*cos(4.0*q.x), min(2.0+exp(-q.x),8.0));
    body_d = smin(body_d, length(q-vec3(0,0.1,0.2))-0.2, 0.4);
    body_d = smin(body_d, sdSegment(q, vec3(0.0,0.1,0.1), vec3(0.7,0.18,0.1))-0.1, 0.3);
    vec4 body = vec4(1,0,0, body_d);
    q = p; q.y = -length(vec2(q.y,0.05));
    q = rotx(0.1*PI)*rotz(0.1*PI)*roty(0.12*PI)*(q-vec3(-0.2,-0.5,-0.2));
    vec4 legs = mapCrabLeg(q, vec3(0.8,0,-0.1), vec3(0.7,-0.1,-0.7), vec3(0.3,-0.2,-1.1), vec3(-0.3,-0.3,-1.3), col_required);
    q = p; q.y = -length(vec2(q.y,0.05));
    q = rotx(0.05*PI)*rotz(0.05*PI)*roty(0.2*PI)*(q-vec3(0.0,-0.3,-0.2));
    legs = cmin(legs, mapCrabLeg(q/0.95, vec3(0.8,0,-0.1), vec3(0.7,-0.1,-0.7), vec3(0.4,-0.2,-1.2), vec3(-0.3,-0.3,-1.4), col_required)*vec4(1,1,1,0.95));
    return smin(body, legs, 0.1);
}

vec4 mapCrab(vec3 p, bool col_required) {
    vec3 q = rotx(0.5*PI)*rotz(PI)*roty(-0.05*PI)*rotx(0.08*PI)*(p-vec3(-0.2,0,-0.2));
    vec4 shell = mapCrabShell(q, col_required);
    q = p;
    vec4 body = mapCrabBody(q/0.9, col_required)*vec4(1,1,1,0.9);
    return cmin(shell, body);
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    //return mapCrabShell(p, col_required);
    //return mapCrabLeg(p, vec3(0.8,0,-0.1), vec3(0.7,-0.1,-0.7), vec3(0.4,-0.2,-1.2), vec3(0.0,-0.3,-1.7), col_required);
    //return mapCrabBody(p, col_required);
    return mapCrab(p, col_required);
}
#endif
