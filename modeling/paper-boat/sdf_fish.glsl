#include "common.glsl"


vec4 mapFish(vec3 p, bool col_required) {
    const float sc = 0.7; p /= sc;
    float bound = sdEllipsoid(p-vec3(-0.3,0,0),vec3(2.0,0.8,1.2)), boundw = 0.3;
    if (bound > 0.0) return vec4(1,0,0, (bound+boundw)*sc);
    float rad = 1.0+0.2*exp(-sqr(1.0*(p.x-0.8)))*(p.z+0.5)-0.05*exp(-sqr(1.0*(p.x+0.6)))-0.2*exp(-sqr(3.0*(p.x-1.3)));
    float lateral_line = p.z-0.15*sin(2.0*p.x)-0.15;
    float lateral_line_d = exp(-sqr(40.0*lateral_line))*exp(-sqr(1.5*(p.x+0.5)));
    vec4 body = vec4(0.85,0.65,0.7, sdEllipsoid(p, vec3(1.3, vec2(0.35*(0.8+0.15*tanh(p.x)),0.5*(0.9+0.1*tanh(p.x)))*rad))
        +0.005*lateral_line_d +0.02*cos(20.0*p.z)*exp(-sqr(20.0*(p.z+0.05)))*exp(-sqr(10.0*(p.x-1.3))) );
    if (col_required)
        body.xyz = mix(body.xyz,
            mix(vec3(0.95,0.6,0.4),
                vec3(0.25,0.4,0.7)*(1.0-max(SimplexNoise3D(8.0*p),0.0)),
                smoothstep(0.,1.,1.5*lateral_line+0.5)) + 0.08*lateral_line_d,
            smoothstep(0.,1.,-1.0*p.x+1.3));
    vec3 q = p; q.y=length(vec2(q.y,0.1));
    q = rotx(-0.04*PI)*rotz(0.12*PI)*(q-vec3(0.8,0.15,-0.05));
    vec4 gill = vec4(0.65,0.35,0.15, sdEllipsoid(q, vec3(0.3+0.04*cos(8.0*q.z),0.1,0.35)));
    body = smin(body, gill, 0.05);
    q = rotz(0.1*PI)*(vec3(p.x,abs(p.y),p.z)-vec3(0.9,0.13,0.1));
    vec4 eyes = vec4(1,0,0, length(q)-0.13);
    if (col_required) eyes.xyz = mix(vec3(0.1,0.05,0.1), vec3(0.55,0.3,0.45), 0.4+0.6*clamp(40.0*(length(q.xz)-0.05)+0.5,0.,1.));
    body = smin(body, eyes, 0.02);
    q = roty(0.05*PI)*(p-vec3(-0.2,0,0.4)); q.z += 0.3*q.x*q.x*q.x;
    float spines = exp(sin(80.0*atan(p.z+0.5,p.x-1.0)));
    vec4 fin_dorsal = vec4(1,0,0, sdEllipsoid(q,
        vec3(0.8,max(0.08*exp(-2.0*q.z),0.01),max(0.2+0.1*exp(-cos(4.0*q.x))+0.05*q.x,0.01)))-0.005*spines);
    if (col_required) fin_dorsal.xyz = pow(vec3(0.9,0.75,0.15),vec3(1.0+0.2*spines));
    body = smin(body, fin_dorsal, 0.05);
    q = roty(-0.05*PI)*(p-vec3(-1.4,0,0.0)); q.z = length(vec2(q.z,0.15));
    q = roty(-0.2*PI)*(q-vec3(-0.0,0,0.2));
    spines = exp(sin(35.0*atan(p.z,p.x+1.0)))*length(vec3(p.z,p.x+1.0,0.1));
    vec4 fin_tail = vec4(1,0,0, sdEllipsoid(q, vec3(0.4,max(0.08*exp(1.0*q.x),0.02),0.15))-0.01*spines);
    if (col_required) fin_tail.xyz = pow(vec3(1.0,0.7,0.2),vec3(1.0+0.4*spines));
    body = smin(body, fin_tail, 0.05);
    q = roty(0.1*PI)*(p-vec3(-0.8,0,-0.3));
    spines = exp(sin(28.0*atan(p.z+0.15,p.x+0.5)));
    vec4 fin_anal = vec4(1,0,0, sdEllipsoid(q, vec3(0.3,max(0.08*exp(2.0*q.x),0.015),0.2))-0.008*spines);
    if (col_required) fin_anal.xyz = pow(vec3(0.95,0.75,0.15),vec3(1.0+0.2*spines));
    body = smin(body, fin_anal, 0.05);
    q = p; q.y = length(vec2(q.y,0.01));
    q = roty(0.2*PI)*rotz(0.1*PI)*(q-vec3(0.1,0.1,-0.4));
    spines = exp(sin(28.0*atan(p.z+0.25,p.x-0.25)));
    vec4 fin_pelvic = vec4(1,0,0, sdEllipsoid(q, vec3(0.25,0.05,0.12))-0.004*spines);
    if (col_required) fin_pelvic.xyz = pow(vec3(0.9,0.6,0.25),vec3(1.0+0.2*spines));
    body = smin(body, fin_pelvic, 0.05);
    q = p; q.y = length(vec2(p.y,0.01));
    q = roty(0.05*PI)*rotz(0.1*PI)*(q-vec3(0.35,0.25,-0.1));
    spines = exp(sin(40.0*atan(p.x-0.6,p.z+0.1)));
    vec4 fin_pectoral = vec4(1,0,0, sdEllipsoid(q, vec3(0.25,0.05,max(0.1-0.2*q.x,0.01)))-0.003*spines);
    if (col_required) fin_pectoral.xyz = pow(vec3(0.95,0.75,0.1),vec3(1.0+0.2*spines));
    body = smin(body, fin_pectoral, 0.05);
    if (col_required) body.xyz = pow(body.xyz, vec3(0.9));
    return body * vec4(1,1,1, sc);
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    return mapFish(p, col_required);
}
#endif
