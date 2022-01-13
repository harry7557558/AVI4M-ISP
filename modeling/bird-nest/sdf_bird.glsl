#include "common.glsl"


vec4 mapBirdHead(vec3 p, bool col_required) {
    vec4 head = vec4(0.95,0.9,0.8,
        sdEllipsoid(p, vec3(0.5, 0.4/(1.0+sqr(0.5*p.x)), (0.32-0.001*cos(-12.0*p.x)+0.005*sin(10.0*p.z)))));
    if (col_required) {
        float black = smax(-p.x + 0.1*exp(sin(12.0*p.z-0.3))-0.27, -1e2 + 0.15-length(p.xz*vec2(0.8,1.1)-vec2(0.0,0.1)), 0.2);
        head.xyz = mix(head.xyz, vec3(0.1,0.05,0.03), clamp(20.0*black-0.5,0.,1.));
        float brown = smin(p.z+0.7*p.x*p.x-0.08-0.1*(1.0-exp(2.0*p.x)), p.z-p.x+0.3+0.15*exp(-sqr(4.0*p.y)), 0.1);
        head.xyz = mix(head.xyz, mix(vec3(0.8,0.4,0.3),vec3(0.4,0.15,0.1),0.5-0.5*cos(6.0*p.x+0.5)+0.5*p.x), clamp(20.0*brown-0.5,0.,1.));
    }
    float beak_z = mix(-0.1,0.0,smootherstep(p.x+1.25));
    vec3 beak_r = vec3((1.0-0.05*exp(-sqr(80.0*(p.z-0.02-beak_z))))*vec2(max(0.2-0.12*p.x,0.01),0.1), 0.1);
    vec4 beak = vec4(0,0,0, sdEllipsoid(p-vec3(-0.45,0,beak_z), beak_r));
    if (col_required) beak.xyz = mix(0.5*vec3(0.2,0.22,0.25), vec3(0.3), 0.5-0.5*cos(80.0*(p.z-0.2*(p.x+0.45))));
    head = smin(head, beak, 0.04);
    vec4 eye = vec4(0.05,0,0.05, sdEllipsoid(p-vec3(-0.32,0.25,0.08), vec3(0.06,0.04,0.04)));
    head = smax(head, vec4(0.8,0.9,0.98,-eye.w), 0.05);
    head = smin(head, eye, 0.01);
    return head;
}

vec4 mapBirdFeather(vec3 p, bool col_required) {
    p -= vec3(0,0,1);
    vec4 bound = vec4(0.55,0.35,0.2,sdEllipsoid(p-vec3(0.04,0,0.25), vec3(0.4,0.2,1.5))), boundw = vec4(0,0,0,0.1);
    if (bound.w > 0.0) return bound+boundw;  // clipping
    vec4 rachis = vec4(0.5,0.45,0.4, sdCapsule(p-vec3(0,0,-1), 2.2, 0.04*(0.5-0.5*tanh(1.5*(p.z-1.2)))));
    vec4 barbs = vec4(0,0,0, 1e8);
    //float barb = sin(120.0*(2.0*pow(abs(p.x),1.4)-p.z));
    float barb = 0.0;
    barbs.w = sdEllipsoid(p-vec3(0,0,0.5), vec3(
        clamp((0.2+0.2*p.x)/(1.0+sqr(0.7*(p.z-1.1))),0.01,0.5),
        0.038/(1.0-sqr(0.1*p.x))*(1.0+0.05*barb),
        1.0));
    if (col_required) barbs.xyz = mix(vec3(0.25,0.2,0.2), vec3(0.95,0.75,0.5),
        smootherstep(-0.4+8.0*(abs(p.x)+1.5*p.y)*(2.0-1.3/(1.0+sqr(0.9*(p.z-0.8))))));
    vec4 feather = smin(rachis, barbs, 0.01);
    feather = mix(feather, bound+boundw, smoothstep(0.,1.,(bound.w+boundw.w)/boundw.w));  // smooth clipping
    return feather;
}

vec4 mapBirdLeg(vec3 p, vec3 joint1, vec3 joint2, vec3 joint3, float digitl, float digitr, float digita0, float digita, bool col_required) {
    float bound = length(p-0.6*joint3)-2.0*length(joint3), boundw = 0.5;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    vec4 thigh = vec4(0.85,0.65,0.5, sdSegment(p, vec3(0.), joint1)-0.07);
    vec4 shank = vec4(0.95,0.7,0.5, sdSegment(p, joint1, joint2)-0.05);
    vec4 foot = vec4(0.95,0.6,0.5, sdSegment(p, joint2, joint3)-0.04);
    vec4 leg = smin(smin(thigh, shank, 0.05), foot, 0.03);
    vec3 w = roty(digitr)*normalize(joint3-joint2);
    vec3 v0 = vec3(0,1,0);
    vec3 u0 = cross(v0, w);
    for (float digit=ZERO; digit<4.; digit++) {
        float a = mix(0.1*PI, 1.1*PI, pow(digit/3.,1.0));
        vec3 u = cos(a)*v0+sin(a)*u0;
        vec3 q = vec3(0.0);
        vec4 toe = vec4(0,0,0, 1e8);
        for (float i=ZERO; i<digit+2.; i++) {
            float a = digita0 + i*digita;
            vec3 dq = digitl/pow(i+1.,1.0) * (u * cos(a) + w * sin(a));
            vec4 bone = vec4(0,0,0, sdSegment(p-joint3, q, q+dq)-0.04*exp(-0.2*i));
            if (col_required) bone.xyz = pow(mix(vec3(0.5,0.1,0.05), vec3(0.55,0.3,0.1), exp(-0.3*i)), vec3(0.4545));
            toe = smin(toe, bone, 0.01);
            q += dq;
        }
        leg = smin(leg, toe, 0.05);
    }
    leg.w = mix(leg.w, bound+boundw, smoothstep(0.,1.,(bound+boundw)/boundw));  // smooth clipping
    return leg;
}

vec4 mapBirdBody(vec3 p, bool col_required) {
    vec4 body = vec4(0.95,0.8,0.7,
        sdEllipsoid(p-vec3(0.05*cos(2.5*p.z)+0.05*sin(1.5*p.z)*sin(p.x)-0.02*sin(3.0*p.z),0,0),
            vec3(0.48+0.07*p.z, 0.5/(1.0+sqr(0.3*(p.z-0.5))), 1.2)) );
    if (col_required) {  // body texture
        float t = 0.5*(2.0*p.x-0.4*p.z+0.3*p.z*p.z+0.5);
        vec2 uv = vec2(cos(atan(p.y,p.x+0.2*p.z+0.2))+exp(0.4*p.z),p.z);
        t += 0.2*GradientNoise2D(vec2(14.0,8.0)*uv);
        body.xyz = mix(body.xyz, vec3(0.55,0.35,0.2), 0.85*smootherstep(t));
    }
    //return body;
    vec4 tail = vec4(0,0,0, 1e8);
    for (float i=ZERO, t; bool(t=(i+0.5)/4.) && i<4.; i++) {  // tail feathers
        vec3 q = mix(vec3(0,0.05,-0.9), vec3(0,0.2,-0.9), t);
        float a = mix(-0.05,0.1,t);
        float r = mix(0.05*PI,0.4*PI,t*t);
        vec3 w = vec3(0.,sin(a),-cos(a));
        vec3 v = vec3(cos(r),sin(r),0.); v = normalize(v-dot(v,w)*w);
        vec3 u = cross(v, w);
        float s = mix(0.45,0.55,5.7*t*t*(1.05-t));
        vec4 feather = mapBirdFeather(transpose(roty(0.1)*mat3(-u,-v,w))*(p-q)/s, col_required)*s;
        tail = cmin(tail, feather);
    }
    for (float i=ZERO, t; bool(t=i/4.+1e-3) && i<=3.; i++) {  // tail converts
        vec3 q = mix(vec3(0.03,0.0,-0.8), vec3(0.03,0.15,-0.8), t);
        float a = mix(0.0,0.1,t);
        float r = mix(0.05*PI,0.3*PI,t*t);
        vec3 w = vec3(0.,sin(a),-cos(a));
        vec3 v = vec3(cos(r),sin(r),0.); v = normalize(v-dot(v,w)*w);
        vec3 u = cross(v, w);
        float s = mix(0.35,0.45,4.0*t*(1.0-t));
        vec4 feather = mapBirdFeather(transpose(roty(0.0)*mat3(u,v,w))*(p-q)/s, col_required)*s;
        tail = cmin(tail, feather+vec4(0.1,0.08,0.05,0));
    }
    if (col_required) tail.xyz = tail.xyz * 1.4*vec3(0.9,0.85,0.8) + vec3(0.05);
    return smin(body, tail, 0.1);
}

vec4 mapBirdWing(vec3 p, vec3 joint1, vec3 joint2, vec3 tip, bool col_required) {
    float bound = sdEllipsoid(p-vec3(1.0,0.0,1.0), vec3(2.0,0.8,2.0)), boundw = 0.3;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    p.y += 0.08*pow(abs(p.z-joint2.z),2.0);
    // approximately on xOz plane
    // https://www.deviantart.com/key-feathers/art/Bird-Wings-Tutorial-201057544
    // https://monikazagrobelna.com/wp-content/uploads/2019/07/how-to-draw-wings-ventral-dorsal-bird-3-1.png?resize=700%2C810
    //float it = 0.75+0.25*sin(iTime); joint2 = it*joint2, tip = mix(joint1, tip, it);
    vec4 upperarm = vec4(0.45,0.3,0.2, sdSegment(p, vec3(0.), joint1)-0.1);
    vec4 forearm = vec4(0.35,0.3,0.3, sdSegment(p, joint1, joint2)-0.07);
    vec4 digits = vec4(0.5,0.45,0.45, smin(sdSegment(p, joint2, mix(joint2,tip,0.5))-0.05, sdSegment(p, mix(joint2,tip,0.5), tip)-0.04, 0.05));
    //vec4 wing = smin(smin(digits, forearm, 0.08), upperarm, 0.1);
    vec4 wing = smin(digits, upperarm, 0.1);
    vec4 skin = smin(
        vec4(0.55,0.5,0.55, sdTriangle(p, vec3(0), joint1, joint2) - 0.08),
        vec4(0.5,0.45,0.4, sdTriangle(p, joint1, joint2, tip) - 0.05),
        0.06);
    wing = smin(wing+vec4(0,0,0,0.03), skin, 0.1);
    vec4 feathers = vec4(0,0,0, 1e8);
    for (float i=ZERO+1., t; bool(t=i/7.) && i<=7.; i++) {  // primary feathers
        vec3 q = mix(joint2, tip, t);
        float a = mix(0.35*PI,0.0,t);
        float r = mix(0.05*PI,0.25*PI,t*t);
        vec3 w = normalize(cos(a)*(tip-joint2)+sin(a)*(tip-joint2).zyx*vec3(1,1,-1));
        vec3 v0 = vec3(-0.1,1.0,0.0); v0 = normalize(v0-dot(v0,w)*w);
        vec3 u0 = cross(v0, w);
        vec3 u = cos(r)*u0+sin(r)*v0, v=-sin(r)*u0+cos(r)*v0;
        float s = mix(0.6,0.7,5.7*t*t*(1.05-t));
        vec4 feather = mapBirdFeather(transpose(mat3(u,v,w))*(p-q)/s, col_required)*s;
        if (col_required) feather.xyz = pow(feather.xyz, vec3(0.85));
        feathers = cmin(feathers, feather);
    }
    for (float i=ZERO, t; bool(t=(i+0.5)/8.) && i<8.; i++) {  // secondary feathers
        vec3 q = mix(joint1, joint2, t);
        float a = mix(0.75*PI,0.55*PI,t);
        float r = mix(0.0*PI,0.1*PI,t);
        vec3 w = normalize(cos(a)*(joint2-joint1)+sin(a)*(joint2-joint1).zyx*vec3(1,1,-1));
        vec3 v0 = vec3(mix(-0.08,-0.02,t),1.0,0.0); //v0 = normalize(v0-dot(v0,w)*w);
        vec3 u0 = cross(v0, w);
        vec3 u = cos(r)*u0+sin(r)*v0, v=-sin(r)*u0+cos(r)*v0;
        float s = mix(0.4,0.6,t);
        vec4 feather = mapBirdFeather(transpose(mat3(u,v,w))*(p-q)/s, col_required)*s;
        if (col_required) feather.xyz = pow(feather.xyz, vec3(0.85));
        feathers = cmin(feathers, feather);
    }
    for (float i=ZERO, t; bool(t=(i+0.5)/7.) && i<7.; i++) {  // primary converts
        vec3 q = mix(joint2, tip, t);
        float a = mix(0.35*PI,0.0*PI,t);
        float r = mix(0.0*PI,0.1*PI,t);
        vec3 w = normalize(cos(a)*(tip-joint2)+sin(a)*(tip-joint2).zyx*vec3(1,1,-1));
        vec3 v0 = vec3(mix(-0.2,-0.3,t),1.0,0.0); //v0 = normalize(v0-dot(v0,w)*w);
        vec3 u0 = cross(v0, w);
        vec3 u = cos(r)*u0+sin(r)*v0, v=-sin(r)*u0+cos(r)*v0;
        float s = mix(0.3,0.35,5.7*t*t*(1.05-t));
        vec4 feather = mapBirdFeather(transpose(mat3(u,v,w))*(p-q)/s, col_required)*s;
        if (col_required) feather.xyz = pow(feather.xyz, vec3(0.8)) + vec3(0.05);
        feathers = cmin(feathers, feather);
    }
    for (float i=ZERO+1., t; bool(t=i/8.) && i<=8.; i++) {  // secondary converts
        vec3 q = mix(joint1, joint2, t);
        float a = mix(0.75*PI,0.55*PI,t);
        float r = mix(0.0*PI,0.1*PI,t);
        vec3 w = normalize(cos(a)*(joint2-joint1)+sin(a)*(joint2-joint1).zyx*vec3(1,1,-1));
        vec3 v0 = vec3(-0.1,1.0,0.0); //v0 = normalize(v0-dot(v0,w)*w);
        vec3 u0 = cross(v0, w);
        vec3 u = cos(r)*u0+sin(r)*v0, v=-sin(r)*u0+cos(r)*v0;
        float s = mix(0.2,0.3,3.0*t*(1.0-t)+0.5*t);
        vec4 feather = mapBirdFeather(transpose(mat3(u,v,w))*(p-q)/s, col_required)*s;
        if (col_required) feather.xyz = pow(feather.xyz, vec3(0.8)) + vec3(0.05);
        feathers = cmin(feathers, feather);
    }
    if (col_required) feathers.xyz = vec3(1.25,1.2,1.15)*pow(feathers.xyz,vec3(0.9)) + vec3(0.1,0.08,0.05);
    wing = smin(wing, feathers, 0.1);
    wing.w = mix(wing.w, bound+boundw, smoothstep(0.,1.,(bound+boundw)/boundw));  // smooth clipping
    return wing;
}

vec4 mapBird(vec3 p, bool col_required) {
    float bound = max(max(max(max(p.x-1.5,-p.x-1.4),abs(p.y)-1.2),max(p.z-1.7,-p.z-1.8)),-0.707*(1.1+p.x+p.z))-0.2, boundw = 0.2;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    float y0 = p.y;
    p.y = length(vec2(p.y,0.02));
    vec4 head = mapBirdHead(roty(0.2)*(p-vec3(-0.75,0,0.75))/0.9, col_required)*vec4(1,1,1,0.9);
    vec4 body = mapBirdBody(roty(0.7)*p, col_required);
    body = smin(head, body, 0.2);
    vec4 leg = mapBirdLeg(roty(-0.45)*(p-vec3(0.25,0.25-0.08*y0,-0.35))*vec3(1,-1,1)/0.7,
        vec3(-0.6,0.05,-0.15), vec3(-0.2+0.1*y0,-0.0,-0.3), vec3(-0.6+0.1*y0,-0.1,-0.45-0.25*y0),
        0.17-0.05*y0, -0.4*PI, 0.15*PI-0.5*y0, 0.1*PI-0.2*y0, col_required)*0.7;
    body = smin(body, leg, 0.15);
    //return body;
    vec3 rzyx = mix(vec3(-0.2,-0.3,0.3), vec3(-1.2,0.2,0.9), 0.2);
    vec3 wing_p = rotz(rzyx.x)*roty(rzyx.y)*rotx(rzyx.z)*(p-vec3(-0.1,0.25,0.55))*vec3(1,-1,1);
    vec4 wing = mapBirdWing(wing_p/0.7,
         vec3(0.5,0.0,0.05),vec3(0.4,0.05,0.7),vec3(0.7,0.1,1.1), col_required)*0.7;
    body = smin(body, wing, 0.15*exp(-8.0*dot(wing_p,wing_p)));
    return body;
}


#ifndef NO_MAP
vec4 map(vec3 p, bool col_required) {
    //return mapBirdFeather(p+vec3(0,0,1), col_required);
    //return mapBirdWing(p-vec3(-1,0,-1), vec3(0.5,0.0,0.15), vec3(0.4,0.05,0.8), vec3(0.7,0.1,1.2), col_required);
    return mapBird(p, col_required);
}
#endif
