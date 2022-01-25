#include "common.glsl"

#iChannel1 "keep-off-grass.png"


vec4 mapSign(vec3 p, bool col_required) {
    float bound = sdBox(p,vec3(1.6,0.6,1.6)), boundw = 0.4;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);
    p -= vec3(0,0,0.5);
    vec4 sign_body = vec4(1,0,0, sdLnNormEllipsoid(p,vec3(1.2,0.15,0.55), 6.0));
    if (col_required) {
        vec2 uv = p.xz;
        float swirl = sin(10.0*(4.0*uv.y+GradientNoise2D(2.0*uv)+0.3*GradientNoise2D(4.0*uv))) + 2.0*GradientNoise2D(4.0*uv);
        sign_body.xyz = pow(mix(vec3(0.2,0.4,0.35),vec3(0.2,0.6,0.45),smoothstep(0.,1.,0.5+0.3*swirl)),vec3(2.0));
    }
    vec2 uv = clamp(0.5+0.6*p.xz/vec2(1.2,0.6), 0., 1.);
    float word = 1.0-texture(iChannel1, uv).x;
    if (p.y < 0.) {
        sign_body.w -= 0.01*word;
        sign_body.xyz = mix(sign_body.xyz, vec3(0.88,0.9,0.92), smoothstep(0.,1.,2.0*(word)));
    }
    vec3 q = vec3(abs(p.x)-0.6, p.y-0.05, p.z+0.5);
    vec4 sign_stick = vec4(pow(vec3(0.2,0.4,0.35),vec3(1.2)), sdLnNormEllipsoid(q,vec3(0.2*exp(0.2*q.z),0.1,1.2), 6.0));
    vec4 sign_ = smin(sign_body, sign_stick, 0.02);
    if (col_required) sign_.xyz = pow(sign_.xyz, vec3(0.6));
    return sign_;
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    return mapSign(p, col_required);
}
#endif
