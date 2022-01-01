#include "common.glsl"


vec4 mapCrabShell(vec3 p, bool col_required) {
    float bound = length(p-vec3(0.55,0.15,0.2))-1.6, boundw = 0.4;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);
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
    float texture = 0.001*r*sin(40.*t) + 0.003*r*sin(40.*atan(dy,dx));
    d = abs(d)-0.0*r-0.25 + texture;
    vec4 shell = vec4(0.85,0.8,0.75, d);
    if (col_required) {
        shell.xyz = mix(vec3(0.65,0.55,0.45), vec3(0.9,0.8,0.75), smoothstep(0.,1.,exp(2.0*PI*n+t-PI)));
        shell.xyz *= pow(vec3(0.57,0.35,0.13), vec3(-20.0*texture));
        shell.xyz *= pow(vec3(0.5,0.35,0.22), vec3(-0.1*sin(2.0*atan(dy,dx))));
    }
    vec4 opening = vec4(0.85,0.8,0.75, sdTorus((p-vec3(-1.65,0,-0.91)).xzy*vec3(1,0.96,1), 1.55, 0.04*(0.0*r+8.0)));
    shell = smin(shell, opening, 0.1) + vec4(0,0,0,0.16);
    return shell * vec4(1,1,1,sc);
}


vec4 mapCrabLeg(vec3 p, vec3 joint1, vec3 joint2, vec3 joint3, vec3 tip, bool col_required) {
    vec3 center = (joint1+joint2+joint3+tip)/5.0;
    float bound = sdEllipsoid(p-center, length(center)*vec3(1.5,1.0,2.0)), boundw = 0.6;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);
    vec3 p0, p1, dp; float h;
    vec4 seg0 = vec4(0.55,0.35,0.15, length(p-vec3(0,0.05,0))-(0.09+0.1*p.y));
    p0 = mix(vec3(0.), joint1, 0.1), p1 = mix(vec3(0.), joint1, 0.9);
    float seg1_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg1_d -= 0.12-0.05*cos(normalize(dp).z)+0.1*h*(1.0-h);
    vec4 seg1 = vec4(0.8,0.2,0.05, seg1_d);
    p0 = mix(joint1, joint2, 0.1), p1 = mix(joint1, joint2, 0.9);
    float seg2_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg2_d -= 0.12-0.05*cos(normalize(dp).x)-0.05*dp.x+0.1*h*(1.0-h);
    vec4 seg2 = vec4(0.9,0.45,0.15, seg2_d);
    p0 = mix(joint2, joint3, 0.1), p1 = mix(joint2, joint3, 0.9);
    float seg3_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg3_d -= 0.12-0.05*cos(normalize(dp).x)-0.05*dp.x+0.1*h*(1.0-h);
    vec4 seg3 = vec4(0.95,0.6,0.3, seg3_d);
    p0 = mix(joint3, tip, 0.1), p1 = mix(joint3, tip, 0.9);
    float seg4_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg4_d -= (0.12-0.05*cos(normalize(dp).x)-0.05*dp.x+0.1*h*(1.0-h))*exp(-0.8*h);
    vec4 seg4 = vec4(0.95,0.5,0.1, seg4_d);
    return smin(smin(smin(seg0, seg1, 0.05), smin(seg2, seg3, 0.01), 0.01), seg4, 0.01) - vec4(0,0,0,0.02);
}

vec4 mapCrabClaw(vec3 p, vec3 joint1, vec3 joint2, vec3 tip1, vec3 tip2, bool col_required) {
    vec3 center = (joint1+joint2+0.5*(tip1+tip2))/4.0;
    float bound = sdEllipsoid(p-center, length(center)*vec3(1.5,1.0,1.3)), boundw = 0.6;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);
    vec3 p0, p1, dp; float h, t;
    vec4 seg0 = vec4(0.55,0.35,0.15, length(p-vec3(0,0.05,0))-(0.09+0.1*p.y));
    p0 = mix(vec3(0.), joint1, 0.1), p1 = mix(vec3(0.), joint1, 0.9);
    float seg1_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg1_d -= (0.15+0.02*cos(3.0*atan(dp.y,dp.z))-0.02*cos(PI*normalize(dp).z)+0.15*h*(1.0-h))*exp(-(1.0-h));
    vec4 seg1 = vec4(0.85,0.1,0.05, seg1_d);
    p0 = mix(joint1, joint2, 0.2), p1 = mix(joint1, joint2, 0.9);
    float seg2_d = sdSegmentH(p, p0, p1, h); dp = p-mix(p0,p1,h);
    seg2_d -= 0.14+0.01*cos(3.0*atan(dp.x,dp.y))-0.02*normalize(dp).x*h+0.05*sin(1.4*h);
    vec4 seg2 = vec4(0.95,0.45,0.1, seg2_d);
    p0 = mix(joint2, tip1, 0.1), p1 = mix(joint2, tip1, 0.9);
    t = abs(p1-p0).x>abs(p1-p0).z ? (p.x-p0.x)/(p1.x-p0.x) : (p.z-p0.z)/(p1.z-p0.z);
    float seg3_d = sdSegmentH(p-vec3(-0.08*sin(PI*t)*vec2(1),0.), p0, p1, h); dp = p-mix(p0,p1,h);
    seg3_d -= 0.08-0.05*dp.x+0.07*sin(3.5*h);
    p0 = mix(joint2, mix(tip1,tip2,1.2), 0.4);
    seg3_d = smin(seg3_d, length(p-p0)-0.15, 0.2);
    p0 = mix(joint2, mix(tip1,tip2,-0.1), 0.3);
    seg3_d = smin(seg3_d, length(p-p0)-0.1, 0.2);
    vec4 seg3 = vec4(0.8,0.55,0.25, seg3_d);
    p0 = mix(joint2, tip2, 0.3), p1 = mix(joint2, tip2, 0.9);
    t = abs(p1-p0).x>abs(p1-p0).z ? (p.x-p0.x)/(p1.x-p0.x) : (p.z-p0.z)/(p1.z-p0.z);
    float seg4_d = sdSegmentH(p-vec3(0.08*sin(PI*t)*vec2(1),0.), p0, p1, h); dp = p-mix(p0,p1,h);
    seg4_d -= 0.07-0.03*normalize(dp).z+0.05*sin(3.5*h);
    vec4 seg4 = vec4(0.8,0.35,0.11, seg4_d);
    return smin(smin(smin(seg1, seg2+0.02, 0.02), smin(seg3, seg4, 0.03)+0.01, 0.02), seg0, 0.05) - vec4(0,0,0,0.0);
}


vec4 mapCrabEye(vec3 p, bool col_required) {
    p.z *= 0.9;
    vec3 q = p-vec3(0,0,0.2*p.x*(1.0-p.x));
    vec4 rod = vec4(1,0,0, sdSegment(q, vec3(0.), vec3(1.,0.,0.)) - (0.3+0.03*sin(1.5*p.x)+0.05*sin(3.0*p.x)));
    q = roty(-0.2*PI)*(p-vec3(1.0,0,0.1));
    vec4 eye = vec4(1,0,0, sdEllipsoid(q, vec3(0.4,0.35,0.45)));
    vec4 res = smin(rod, eye, 0.2) + vec4(0,0,0,0.1);
    if (col_required) res.xyz = mix(vec3(0.9,0.75,0.65), vec3(0.1,0.05,0.02), smoothstep(0.,1.,20.0*(p.x-1.0)));
    return res;
}

vec4 mapCrabBody(vec3 p, bool col_required) {
    p.z -= 0.1*p.y, p.x -= 0.15*p.y;
    vec3 q = roty(0.1*PI)*p; q.y=0.2*log(2.0*cosh(q.y/0.2)), q.z-=0.1;
    float body_d = sdLnNormEllipsoid(q, vec3(0.8,0.5,0.38)+0.02*cos(4.0*q.x), min(2.0+exp(-q.x),8.0));
    body_d = smin(body_d, length(q-vec3(0,0.1,0.2))-0.2, 0.4);
    body_d = smin(body_d, sdSegment(q, vec3(0.0,0.1,0.1), vec3(0.7,0.18,0.1))-0.1, 0.3);
    vec4 body = vec4(mix(vec3(0.95,0.65,0.25),vec3(0.85,0.25,0.1),smoothstep(0.,1.,5.0*q.z+0.5)), body_d);
    // legs
    vec4 legs = vec4(1,0,0, 100.);
    float ls = 0.5+0.5*tanh(4.0*p.y);
    q = p; q.y = -length(vec2(q.y,0.05));
    q = rotx(0.1*PI)*rotz(0.1*PI)*roty((0.1+0.1*ls)*PI)*(q-vec3(-0.15,-0.5,-0.2));
    legs = cmin(legs, mapCrabLeg(q, vec3(0.8,0,-0.1), vec3(0.7,-0.1,-0.7), vec3(0.3+0.2*ls,-0.2,-1.1-0.1*ls), vec3(-0.3+0.3*ls,-0.3,-1.3-0.3*ls), col_required));
    q = p; q.y = -length(vec2(q.y,0.05));
    q = rotx(0.05*PI)*rotz(0.05*PI)*roty((0.12+0.1*ls)*PI)*(q-vec3(0.05,-0.35,-0.2));
    legs = cmin(legs, mapCrabLeg(q/0.95, vec3(0.8,0,-0.1), vec3(0.7,-0.1,-0.7), vec3(0.4,-0.2,-1.2), vec3(-0.3,-0.3,-1.4-0.2*ls), col_required)*vec4(1,1,1,0.95));
    q = rotz(0.05*PI)*roty(0.0*PI)*(p-vec3(0.15,-0.15,-0.2));
    legs = cmin(legs, mapCrabClaw(q/1.05, vec3(1.0,0,0), vec3(1.1,0,-0.4), vec3(0.8,-0.1,-1.3), vec3(0.9,-0.05,-1.3), col_required)*vec4(1,1,1,1.05));
    q = rotz(0.03*PI)*roty(0.1*PI)*(vec3(1,-1,1)*p-vec3(0.15,-0.2,-0.2));
    legs = cmin(legs, mapCrabClaw(q/1.05, vec3(1.0,0,0), vec3(1.3,0,0.1), vec3(2.0,-0.1,0.3), vec3(2.05,0.1,0.3), col_required)*vec4(1,1,1,1.05));
    body = smin(body, legs, 0.1);
    if (col_required) {
        body.xyz += 0.4*SimplexNoise3D(2.0*p);
        body.xyz *= pow(vec3(0.55,0.2,0.05), vec3(1.5*smax(SimplexNoise3D(20.0*p)-0.1, 0., 0.05)));
        body.xyz = pow(clamp(body.xyz, 0.0, 1.0), vec3(1.0));
    }
    // eyes/tentacles
    q = roty(0.15*PI)*(p-vec3(0.6,0,0.3)); q.y=0.05*log(2.0*cosh(q.y/0.05));
    vec4 eyes = mapCrabEye((q-vec3(0,0.15,0))/0.35, col_required)*vec4(1,1,1,0.35);
    body = smin(body, eyes, 0.03);
    vec4 tentacles = vec4(0.9,0.3,0.05, sdSegment(q, vec3(0,0.05,0.1), vec3(0.2,0.12,-0.3))-0.03);
    body = smin(body, tentacles, 0.04);
    return body;
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
    //return mapCrabClaw(p, vec3(1.0,0,0), vec3(1.1,0,-0.4), vec3(0.8,-0.1,-1.3), vec3(0.9,-0.05,-1.3), col_required);
    //return mapCrabClaw(p, vec3(1.0,0,0), vec3(1.3,0,0.1), vec3(2.0,-0.1,0.3), vec3(2.05,0.1,0.3), col_required);
    //return mapCrabEye(p, col_required);
    //return mapCrabBody(p, col_required);
    return mapCrab(p, col_required);
}
#endif
