#include "common.glsl"


vec4 mapTreeBranch(vec3 p, bool col_required) {
    float th = 0.05 + 0.02*(4.*p.z*(1.0-p.z)) - 0.02*p.z;
    return vec4(mix(vec3(0.7,0.15,0.1), vec3(0.95,0.3,0.3), pow(max(p.z,0.),2.0)),
        sdSegment(p,vec3(0,0,0),vec3(0,0,1.0))-th);
}

vec4 mapCoralTree01(vec3 p, bool col_required) {
    vec3 q = p+vec3(0,0,1);
    const int nbranch = 3;
    const int depth = 4;
    int stk_b[depth];
    mat4 stk_t[depth];
    int stkptr = 0;
    mat4 transform = mat4(1.2); transform[3] = vec4(0,0,0,1);
    stk_b[stkptr] = -1, stk_t[stkptr] = transform, stkptr++;
    vec3 q1 = (stk_t[0]*vec4(q,1.)).xyz;
    float sc = pow(determinant(mat3(stk_t[0])),1./3.);
    vec4 tree = mapTreeBranch(q1, col_required) * vec4(0.55,0.3,0.5,1./sc);
    while (stkptr >= 0) {
        do { stkptr--; } while (stkptr>=0 && stk_b[stkptr]>=nbranch-1);
        if (stkptr<0) break;
        stk_b[stkptr] += 1;
        for (stkptr++; stkptr<depth; stkptr++) {
            float t = exp(-float(stkptr-1));
            if (stk_b[stkptr-1] == 0) {  // bottom branch
                transform = mat4(1.6*rotz(0.6*PI)*roty(mix(0.2,0.3,1.0-t)*PI)*rotx(0.1*PI));
                transform[3] = transform*vec4(-vec3(0,0,mix(0.2,0.8,1.0-t)), 1);
            }
            else if (stk_b[stkptr-1] == 1) {  // top branch
                transform = mat4(1.1*rotz(0.1*PI)*roty(0.1*PI)*rotx(0.0*PI));
                transform[3] = transform*vec4(-vec3(0,0,1.0), 1);
            }
            else {  // side branch
                transform = mat4(1.4*rotz(-0.5*PI+float(stkptr))*roty(-0.1*PI)*rotx(-0.1*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.8), 1);
            }
            stk_b[stkptr] = 0;
            stk_t[stkptr] = transform * stk_t[stkptr-1];
            q1 = (stk_t[stkptr]*vec4(q,1.)).xyz;
            sc = pow(determinant(mat3(stk_t[stkptr])),1./3.);
            if (sc > 3.0) break;
            vec4 branch = mapTreeBranch(q1, col_required) / vec4(1,1,1,sc);
            branch.xyz *= pow(mix(vec3(0.55,0.3,0.5), vec3(1.0), 1.0-t), vec3(0.6));
            tree = smin(tree, branch, 0.05);
        }
    }
    tree.w -= 0.02;
    tree.xyz = pow(vec3(1.1,1.0,1.1)*tree.xyz, vec3(0.7));
    return tree;
}

vec4 mapCoralTree02(vec3 p, bool col_required) {
    vec3 q = p+vec3(0,0,1);
    const int nbranch = 3;
    const int depth = 4;
    int stk_b[depth];
    mat4 stk_t[depth];
    int stkptr = 0;
    mat4 transform = mat4(1.2); transform[3] = vec4(0,0,0,1);
    stk_b[stkptr] = -1, stk_t[stkptr] = transform, stkptr++;
    vec3 q1 = (stk_t[0]*vec4(q,1.)).xyz;
    float sc = pow(determinant(mat3(stk_t[0])),1./3.);
    vec4 tree = mapTreeBranch(q1, col_required) * vec4(0.55,0.3,0.5, 1./sc);
    while (stkptr >= 0) {
        do { stkptr--; } while (stkptr>=0 && stk_b[stkptr]>=nbranch-1);
        if (stkptr<0) break;
        stk_b[stkptr] += 1;
        for (stkptr++; stkptr<depth; stkptr++) {
            float t = exp(-float(stkptr-1));
            if (stk_b[stkptr-1] == 0) {  // top branch
                transform = mat4(1.1*rotz((0.9+float(stkptr))*PI)*roty(cos(PI*float(stkptr))*0.05*PI));
                transform[3] = transform*vec4(-vec3(0,0,1.0), 1);
            }
            else if (stk_b[stkptr-1] == 1) {  // side branch
                transform = mat4(1.8*rotz(-0.9*PI)*roty(0.2*PI)*rotx(0.1*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.5), 1);
            }
            else {  // side branch 2
                transform = mat4(1.8*roty(0.2*PI)*rotx(0.1*PI)*rotz(0.9*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.9), 1);
            }
            stk_b[stkptr] = 0;
            stk_t[stkptr] = transform * stk_t[stkptr-1];
            q1 = (stk_t[stkptr]*vec4(q,1.)).xyz;
            sc = pow(determinant(mat3(stk_t[stkptr])),1./3.);
            if (length(stk_t[stkptr-1][3].xy)+0.3*float(stkptr-3)+0.3*sc > 2.0) break;
            vec4 branch = mapTreeBranch(q1, col_required) / vec4(1,1,1,sc);
            branch.xyz *= pow(mix(vec3(0.55,0.3,0.5), vec3(1.0), 1.0-t), vec3(0.6));
            tree = smin(tree, branch, 0.05);
        }
    }
    tree.w -= 0.02;
    tree.xyz = pow(vec3(1.1,1.0,1.1)*tree.xyz, vec3(0.7));
    return tree;
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    return mapCoralTree02(p, col_required);
}
#endif
