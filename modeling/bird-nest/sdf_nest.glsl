#include "common.glsl"


vec4 mapNest(vec3 p, bool col_required) {
    float bound = max(sdEllipsoid(p, vec3(2.2,2.2,1.8)), p.z-0.6), boundw = 0.2;
    if (bound > 0.0) return vec4(1,0,0, (bound+boundw)*0.8);  // clipping
    p.y = 1.05*p.y + 0.25*tanh(p.y)*cos(p.x)*(0.7+0.3*sin(2.0*p.z))*(1.0+0.2*sin(4.0*p.y)) - 0.03*cos(2.5*p.x);
    p.z -= -0.1*sin(p.y) + 0.05*cos(2.5*p.y) + 0.05*sin(p.x);
    vec4 res = vec4(0, 0, 0, 1e8);
    vec3 col1 = vec3(0.0);
    float seed = 0., seedc = 0.;
    for (float i=ZERO; i<24.; i++) {  // frame
        float t = 0.2 + 1.5 * (i+0.5) / 24.;
        float r = 1.6*sin(t)*(0.9+0.2*rand(seed));
        float z = pow(1.0-cos(t),1.4)-1.0;
        float rot_phi = 0.05*PI*(-1.+2.*rand(seed));
        vec3 q = rotx(rot_phi)*(p-vec3(0.,0.,z));
        float d1 = sdTorus(q, r, 0.04);
        if (col_required) col1 = mix(vec3(0.25,0.15,0.1), vec3(0.65,0.66,0.5), rand(seedc));
        res = cmin(res, vec4(col1, d1));
    }
    for (float i=ZERO; i<16.; i++) {  // lower wires
        float t = 1.4 * (i+0.5) / 16.;
        float r = 1.6*sin(t)*(0.9+0.1*rand(seed)) * (1.0+0.02*sin(7.0*atan(p.x,p.y)));
        float z = pow(1.0-cos(t),1.4)-1.0+0.01*r*sin(8.0*atan(p.y,p.x));
        float rot_theta = 2.0*PI*rand(seed);
        float rot_phi = 0.05*PI*(1.0+rand(seed));
        vec3 q = rotx(rot_phi)*rotz(rot_theta)*(p-vec3(0.,0.,z));
        float d1 = sdTorus(q, r, 0.03);
        if (col_required) col1 = mix(vec3(0.35,0.35,0.3), vec3(0.4,0.35,0.22), rand(seedc));
        res = cmin(res, vec4(col1, d1));
    }
    for (float i=ZERO; i<8.; i++) {  // upper wires
        float t = 1.2 + 0.5 * (i+0.5) / 8.;
        float r = 1.6*sin(t)*(0.9+0.2*rand(seed));
        float z = pow(1.0-cos(t),1.4)-1.0+0.06*sin(6.0*atan(p.y,p.x));
        float rot_phi = 0.02*PI*(-1.0+2.0*rand(seed));
        vec3 q = rotx(rot_phi)*(p-vec3(0.,0.,z));
        float d1 = sdTorus(q, r, 0.04);
        if (col_required) col1 = mix(vec3(0.55,0.35,0.25), vec3(0.95,0.85,0.7), rand(seedc));
        res = cmin(res, vec4(col1, d1));
    }
    for (float i=ZERO; i<18.; i++) {  // cross wires
        float t1 = 2.0*PI * i / 18.;
        float t2 = t1+PI + 0.4*PI*(-1.+2.*vanDerCorput(i+100.,2.));
        vec3 v1 = vec3((1.0+0.3*rand(seed))*vec2(cos(t1), sin(t1)), 0.1*randt(seed));
        vec3 v2 = vec3((1.0+0.3*rand(seed))*vec2(cos(t2), sin(t2)), 0.1*randt(seed));
        vec3 q = p + vec3(0,0, 0.7*sqrt(max(2.6-dot(p.xy,p.xy),0.))-0.12);
        float d1 = sdSegment(q, v1, v2) - 0.03;
        if (col_required) col1 = mix(vec3(0.6,0.6,0.55), vec3(0.4,0.4,0.25), rand(seedc));
        res = cmin(res, vec4(col1, d1));
    }
    return res;
}


vec4 map(vec3 p, bool col_required) {
    vec4 d = mapNest(p, col_required);
    return d;
}

