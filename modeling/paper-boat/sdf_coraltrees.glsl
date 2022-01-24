#include "common.glsl"


vec4 mapTreeBranch(vec3 p, bool col_required) {
    float th = 0.05 + 0.02*(4.*p.z*(1.0-p.z)) - 0.02*p.z;
    return vec4(mix(vec3(0.7,0.15,0.1), vec3(0.95,0.3,0.3), pow(max(p.z,0.),2.0)),
        sdSegment(p,vec3(0,0,0),vec3(0,0,1.0))-th);
}

vec4 mapCoralTree01(vec3 p, bool col_required) {
    float bound = max(max(abs(p.x+0.25)-0.88,abs(p.y-0.01)-0.72),abs(p.z-0.24)-1.30), boundw = 0.3;
    if (bound-boundw > 0.0) return vec4(1,0,0, bound);
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
    float bound = max(max(abs(p.x-0.01)-0.63,abs(p.y-0.07)-0.27),abs(p.z-0.43)-1.49), boundw = 0.3;
    if (bound-boundw > 0.0) return vec4(1,0,0, bound);
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
            else {  // side branch
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

vec4 mapCoralTree03(vec3 p, bool col_required) {
    float bound = max(max(abs(p.x+0.66)-1.21,abs(p.y-0.24)-0.26),abs(p.z-0.39)-1.45), boundw = 0.3;
    if (bound-boundw > 0.0) return vec4(1,0,0, bound);
    p.y -= 0.15*cos(1.5*p.x)-0.1;
    vec3 q = p+vec3(0,0,1);
    const int nbranch = 3;
    const int depth = 7;
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
            float t = float(stkptr);
            sc = pow(determinant(mat3(stk_t[stkptr-1])),1./3.);
            float tz = -stk_t[stkptr-1][3].z, tr = length(stk_t[stkptr-1][3].xy);
            if (stk_b[stkptr-1] == 0) {  // top branch
                transform = mat4((1.1+0.2*exp(-2.0*tz))*rotz(0.05*PI)*roty(clamp(0.0+0.04*t*t*exp(-0.5*sc)-0.01*tr*tr,0.05,0.4)*PI));
                transform[3] = transform*vec4(-vec3(0,0,1.0), 1);
            }
            else if (stk_b[stkptr-1] == 1) {  // side branch (-x)
                transform = mat4((2.0-0.1*t+0.5*exp(-4.0*tz))*roty((0.2+0.01*t)*PI)*rotz(0.1*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.5), 1);
            }
            else {  // side branch (+x)
                transform = mat4(max(2.2-0.3*t,1.0)*roty(-0.15*PI)*rotz(-0.1*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.9), 1);
            }
            stk_b[stkptr] = 0;
            stk_t[stkptr] = transform * stk_t[stkptr-1];
            q1 = (stk_t[stkptr]*vec4(q,1.)).xyz;
            sc = pow(determinant(mat3(stk_t[stkptr])),1./3.);
            if (1.15*sc+0.9*t + 0.8*tz > 10.0) break;
            vec4 branch = mapTreeBranch(q1, col_required) / vec4(1,1,1,sc);
            branch.xyz *= pow(mix(vec3(0.55,0.3,0.5), vec3(1.0), 1.0-exp(-t)), vec3(0.6));
            tree = smin(tree, branch, 0.05);
        }
    }
    tree.w -= 0.02;
    tree.xyz = pow(vec3(1.1,1.0,1.1)*tree.xyz, vec3(0.7));
    return tree;
}

vec4 mapCoralTree04(vec3 p, bool col_required) {
    float bound = max(max(abs(p.x+0.13)-0.47,abs(p.y-0.09)-0.16),abs(p.z+0.22)-0.83), boundw = 0.3;
    if (bound-boundw > 0.0) return vec4(1,0,0, bound);
    vec3 q = p+vec3(0,0,1);
    const int nbranch = 3;
    const int depth = 3;
    int stk_b[depth];
    mat4 stk_t[depth];
    int stkptr = 0;
    mat4 transform = mat4(1.6); transform[3] = vec4(0,0,0,1);
    stk_b[stkptr] = -1, stk_t[stkptr] = transform, stkptr++;
    vec3 q1 = (stk_t[0]*vec4(q,1.)).xyz;
    float sc = pow(determinant(mat3(stk_t[0])),1./3.);
    vec4 tree = mapTreeBranch(q1, col_required) * vec4(0.55,0.3,0.5, 1./sc);
    while (stkptr >= 0) {
        do { stkptr--; } while (stkptr>=0 && stk_b[stkptr]>=nbranch-1);
        if (stkptr<0) break;
        stk_b[stkptr] += 1;
        for (stkptr++; stkptr<depth; stkptr++) {
            float t = float(stkptr-1);
            if (stk_b[stkptr-1] == 0) {  // top branch
                transform = mat4(1.2*rotz((0.9+float(stkptr))*PI)*roty(cos(PI*t)*0.05*PI));
                transform[3] = transform*vec4(-vec3(0,0,1.0), 1);
            }
            else if (stk_b[stkptr-1] == 1) {  // side branch
                transform = mat4(1.3*rotz(-0.9*PI)*roty(0.2*PI)*rotx(0.1*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.7*(1.0-0.5*exp(-t))), 1);
            }
            else {  // side branch
                transform = mat4(2.2*roty(0.2*PI)*rotx(0.1*PI)*rotz(0.9*PI));
                transform[3] = transform*vec4(-vec3(0,0,0.8), 1);
            }
            stk_b[stkptr] = 0;
            stk_t[stkptr] = transform * stk_t[stkptr-1];
            q1 = (stk_t[stkptr]*vec4(q,1.)).xyz;
            sc = pow(determinant(mat3(stk_t[stkptr])),1./3.);
            vec4 branch = mapTreeBranch(q1, col_required) / vec4(1,1,1,sc);
            branch.xyz *= pow(mix(vec3(0.55,0.3,0.5), vec3(1.0), 1.0-exp(-t)), vec3(0.6));
            tree = smin(tree, branch, 0.05);
        }
    }
    tree.w -= 0.02;
    tree.xyz = pow(vec3(1.1,1.0,1.1)*tree.xyz, vec3(0.7));
    return tree;
}


vec4 mapRock(vec3 p, bool col_required) {
    float d = sdEllipsoid(p, vec3(0.8,0.6,0.3));
    if (abs(d)>0.5) return vec4(1,0,0,d);
    float dcol = 0.0;
    for (float k=0.0; k<4.0; k+=1.0) {
        float noise = SimplexNoise3D(1.0*exp2(k)*p);
        d += 0.5*exp2(-k)*noise;
        dcol += 1.0*pow(1.0,-k)*noise;
    }
    vec3 col = mix(vec3(0.4,0.4,0.3), vec3(0.4,0.35,0.15), smoothstep(0.,1.,0.5+dcol));
    return vec4(pow(col,vec3(0.7)), d*0.7);
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    //return mapRock(p, col_required);
    return mapCoralTree01(p, col_required);
}
#endif
