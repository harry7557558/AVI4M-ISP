#include "common.glsl"


// paper boat with sand inside it
vec4 mapBoat(vec3 p, bool col_required) {
    vec4 bound = vec4(1,1,1,sdEllipsoid(p-vec3(0.2,0,0.1),vec3(2.6,1.6,1.6))), boundw = vec4(0,0,0,0.3);
    vec3 q = vec3(abs(p.xy), p.z+0.8);
    if (bound.w > 0.0) return bound+boundw;  // clipping
    q *= 1.0-0.1*cos(length(p.xy));
    // modified from https://www.shadertoy.com/view/tlXyzr by blackle
	vec3 mast = vec3(0.02, 0.0, 1.4);
	vec3 mid = vec3(0.01, 0.55, 0.25);
	vec3 port_bow = vec3(1.0, 0.0, 0.0);
	vec3 keel = vec3(0.0, 0.3, 0.0);
	vec3 port = vec3(0.0, 1.1-0.0*sign(p.y), 0.75);
	vec3 bow = vec3(2.1, 0.0, 1.3) + sign(p.x)*vec3(0.1,0.0,0.25);
	float tri1 = sdTriangle(q, port_bow, mid, mast);
	float tri2 = sdTriangle(q, port_bow, bow, port);
	float tri3 = sdTriangle(q, port_bow, port, keel);
    vec4 boat = vec4(0.95,0.9,0.85, min(tri1,min(tri2,tri3))-0.025);
    // sand
    float sign1 = dot(cross(mid-port_bow, mast-port_bow), q-port_bow);
    float sign2 = dot(cross(bow-port_bow, port-port_bow), q-port_bow);
    float sign3 = dot(cross(port-port_bow, keel-port_bow), q-port_bow);
    float sand_d = min(tri1, min(tri2, tri3)) * -sign(min(sign1, min(sign2, sign3)));
    float noise_lf = 0.02*GradientNoise2D(4.0*p.xy);
    float noise_hf = 0.004*GradientNoise2D(18.0*p.xy);
    float sand_face = p.z + 0.1*cos(2.0*p.x)-0.05*p.x-0.12 + noise_lf+noise_hf;
    vec4 sand = vec4(1,0,0, max(sand_d,sand_face)+0.02);
    if (col_required) {
        sand.xyz = vec3(0.7,0.65,0.5)-20.0*vec3(0.95,0.85,0.55)*noise_hf;
        sand.xyz = vec3(1.1,1.0,0.9)*pow(sand.xyz,vec3(0.6));
    }
    // return
    vec4 res = smin(boat, sand, 0.04);
    res = mix(res, bound+boundw, smoothstep(0.,1.,(bound.w+boundw.w)/boundw.w));  // smooth clipping
    return res;
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    return mapBoat(p, col_required);
}
#endif
