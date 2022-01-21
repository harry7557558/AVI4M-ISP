#include "common.glsl"

#iChannel0 "nvidia-logo.png"


vec4 mapChip(vec3 p, bool col_required) {
    vec4 board = vec4(0.2,0.4,0.3, sdBox(p,vec3(1.25,1.25,0.07))-0.05);
    vec4 logo = vec4(0.85,0.85,0.87, sdBox(p-vec3(0,0,0.05),vec3(1,1,0.05))-0.1);
    if (col_required && abs(logo.w)<0.2) {
        vec2 uv = 0.5+0.50*(p.xy-vec2(0.06,0.02))/vec2(1.3125,1.);
        logo.xyz *= texture(iChannel0, clamp(uv,0.,1.)).xyz;
    }
    return cmin(board, logo);
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    return mapChip(p, col_required);
}
#endif
