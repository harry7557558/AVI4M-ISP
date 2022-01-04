// SDF visualizer v2
// Allows a usually larger average marching step without missing thin surface
// Expected to be faster than v1

// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926
#define ZERO min(iTime,0.)
//#define ZERO 0.0


float sqr(float x) { return x*x; }
float sqr(vec2 p) { return dot(p,p); }
float sqr(vec3 p) { return dot(p,p); }

// rotation matrices
mat2 rot2(float a) { return mat2(cos(a), sin(a), -sin(a), cos(a)); }
mat3 rotx(float a) { return mat3(1, 0, 0, 0, cos(a), sin(a), 0, -sin(a), cos(a)); }
mat3 rotz(float a) { return mat3(cos(a), sin(a), 0, -sin(a), cos(a), 0, 0, 0, 1); }
mat3 roty(float a) { return mat3(cos(a), 0, -sin(a), 0, 1, 0, sin(a), 0, cos(a)); }

// geometric primitives
float sdSegment(vec2 p, vec2 a, vec2 b) {
    vec2 q = p-a, d = b-a;
    float h = dot(q,d)/dot(d,d);
    return length(q-d*clamp(h,0.,1.));
}
float sdSegment(vec3 p, vec3 a, vec3 b) {
    vec3 q = p-a, d = b-a;
    float h = dot(q,d)/dot(d,d);
    return length(q-d*clamp(h,0.,1.));
}
float sdCapsule(vec3 p, float h, float r) {
    p.z = p.z-clamp(p.z, 0.0, h);
    return length(p) - r;
}
float sdEllipsoid(vec3 p, vec3 r) {
    float k1 = length(p/r);
    float k2 = length(p/(r*r));
    return k1*(k1-1.0)/k2;
}
float sdTorus(vec3 p, float R, float r) {
    vec2 q = vec2(length(p.xy)-R,p.z);
    return length(q)-r;
}
float sdTriangle(vec3 p, vec3 a, vec3 b, vec3 c) {
    // https://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
    vec3 ba = b - a; vec3 pa = p - a;
    vec3 cb = c - b; vec3 pb = p - b;
    vec3 ac = a - c; vec3 pc = p - c;
    vec3 n = cross(ba, ac);
    return sqrt(
      (sign(dot(cross(ba,n),pa)) +
       sign(dot(cross(cb,n),pb)) +
       sign(dot(cross(ac,n),pc))<2.0)
       ?
       min(min(
       sqr(ba*clamp(dot(ba,pa)/dot(ba,ba),0.0,1.0)-pa),
       sqr(cb*clamp(dot(cb,pb)/dot(cb,cb),0.0,1.0)-pb)),
       sqr(ac*clamp(dot(ac,pc)/dot(ac,ac),0.0,1.0)-pc))
       :
       dot(n,pa)*dot(n,pa)/dot(n,n) );
}


// smoothed blending functions
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return mix(b, a, h) - k * h * (1.0 - h);
}
float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
}
float smootherstep(float x) {
    x = clamp(x, 0., 1.);
    return x*x*x*(10.+x*(-15.+6.*x));
}

// color blending functions
vec4 cmin(vec4 c1, vec4 c2) {
    return c1.w<c2.w ? c1 : c2;
}
vec4 smin(vec4 a, vec4 b, float k) {
    float h = clamp(0.5 + 0.5 * (b.w - a.w) / k, 0., 1.);
    float d = mix(b.w, a.w, h) - k * h * (1.0 - h);
    return vec4(mix(b.xyz, a.xyz, h), d);
}
vec4 smax(vec4 a, vec4 b, float k) {
    return smin(a*vec4(1,1,1,-1), b*vec4(1,1,1,-1), k) * vec4(1,1,1,-1);
}

// random/hash functions from https://www.shadertoy.com/view/4djSRW
vec2 hash22(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

// noise function(s)
float GradientNoise2D(vec2 xy) {
    float i0 = floor(xy.x), i1 = i0 + 1.0;
    float j0 = floor(xy.y), j1 = j0 + 1.0;
    float v00 = dot(2.0 * hash22(vec2(i0, j0)) - 1.0, xy - vec2(i0, j0));
    float v01 = dot(2.0 * hash22(vec2(i0, j1)) - 1.0, xy - vec2(i0, j1));
    float v10 = dot(2.0 * hash22(vec2(i1, j0)) - 1.0, xy - vec2(i1, j0));
    float v11 = dot(2.0 * hash22(vec2(i1, j1)) - 1.0, xy - vec2(i1, j1));
    float xf = xy.x - i0; xf = xf * xf * xf * (10.0 + xf * (-15.0 + xf * 6.0));
    float yf = xy.y - j0; yf = yf * yf * yf * (10.0 + yf * (-15.0 + yf * 6.0));
    return v00 + (v10 - v00)*xf + (v01 - v00)*yf + (v00 + v11 - v01 - v10) * xf*yf;
}


// modeling

vec4 mapBirdHead(vec3 p, bool col_required) {
    vec4 head = vec4(0.95,0.9,0.8,
        sdEllipsoid(p, vec3(0.5, 0.4/(1.0+sqr(0.5*p.x)), (0.32-0.001*cos(-12.0*p.x)+0.005*sin(10.0*p.z)))));
    if (col_required) {
        float black = smax(-p.x + 0.1*exp(sin(12.0*p.z-0.3))-0.27, -1e2 + 0.15-length(p.xz*vec2(0.8,1.1)-vec2(0.0,0.1)), 0.2);
        head.xyz = mix(head.xyz, vec3(0.1,0.05,0.03), saturate(20.0*black-0.5));
        float brown = smin(p.z+0.7*p.x*p.x-0.08-0.1*(1.0-exp(2.0*p.x)), p.z-p.x+0.3+0.15*exp(-sqr(4.0*p.y)), 0.1);
        head.xyz = mix(head.xyz, mix(vec3(0.8,0.4,0.3),vec3(0.4,0.15,0.1),0.5-0.5*cos(6.0*p.x+0.5)+0.5*p.x), saturate(20.0*brown-0.5));
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
    float bound = sdEllipsoid(p-vec3(0.04,0,0.25), vec3(0.4,0.2,1.5)), boundw = 0.1;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
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
    feather.w = mix(feather.w, bound+boundw, smoothstep(0.,1.,(bound+boundw)/boundw));  // smooth clipping
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
        float a = mix(0.1*PI, 1.2*PI, pow(digit/3.,0.7));
        vec3 u = cos(a)*v0+sin(a)*u0;
        vec3 q = vec3(0.0);
        vec4 toe = vec4(0,0,0, 1e8);
        for (float i=ZERO; i<digit+2.; i++) {
            float a = digita0 + i*digita;
            vec3 dq = digitl/sqrt(i+1.) * (u * cos(a) + w * sin(a));
            vec4 bone = vec4(0,0,0, sdSegment(p-joint3, q, q+dq)-0.03*exp(-0.2*i));
            if (col_required) bone.xyz = pow(mix(vec3(0.5,0.1,0.05), vec3(0.55,0.3,0.1), exp(-0.3*i)), vec3(0.5));
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
        float t = 0.5*(2.0*p.x-0.2*p.z+0.3*p.z*p.z+0.5);
        vec2 uv = vec2(cos(atan(p.y,p.x+0.2*p.z))+exp(0.4*p.z),p.z);
        t += 0.2*GradientNoise2D(vec2(14.0,8.0)*uv);
        body.xyz = mix(body.xyz, vec3(0.55,0.35,0.2), 0.9*smootherstep(t));
    }
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
    p.y = length(vec2(p.y,0.02));
    vec4 head = mapBirdHead(roty(0.2)*(p-vec3(-0.75,0,0.75)), col_required);
    vec4 body = mapBirdBody(roty(0.7)*p, col_required);
    body = smin(head, body, 0.2);
    vec4 leg = mapBirdLeg(roty(-0.3)*(p-vec3(-0.1,0.2,-0.4))*vec3(1,-1,1)/0.7,
        vec3(-0.25,0,-0.15), vec3(-0.05,0,-0.4), vec3(-0.2,0,-0.55), 0.12, -0.4, 0.05*PI, 0.1*PI, col_required)*0.7;
    body = smin(body, leg, 0.2);
    vec3 rzyx = mix(vec3(-0.2,-0.3,0.3), vec3(-1.2,0.2,0.9), 0.2);
    vec4 wing = mapBirdWing(rotz(rzyx.x)*roty(rzyx.y)*rotx(rzyx.z)*(p-vec3(-0.1,0.25,0.55))*vec3(1,-1,1)/0.7,
         vec3(0.5,0.0,0.05),vec3(0.4,0.05,0.7),vec3(0.7,0.1,1.1), col_required)*0.7;
    body = smin(body, wing, 0.05);
    return body;
}


vec4 map(vec3 p, bool col_required) {
    //return mapBirdFeather(p+vec3(0,0,1), col_required);
    //return mapBirdWing(p-vec3(-1,0,-1), vec3(0.5,0.0,0.15), vec3(0.4,0.05,0.8), vec3(0.7,0.1,1.2), col_required);
    return mapBird(p, col_required);
}


float sdf(vec3 p) {
    const float sc = 1.0;
    return map(p/sc, false).w*sc;
}


vec3 sdfGrad(in vec3 p, in float e) {
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}



// color surface or not
#define COLOR 1


// raymarching parameters
#define BOX_RADIUS vec3(2.5, 2.5, 2.0)
#define STEP 0.1
#define MIN_STEP 0.005
#define MAX_STEP 100.

// rendering parameters
#define FIELD_EMISSION 0.1
#define ISOSURFACE_FREQUENCY 8.0
#define DISCONTINUITY_OPACITY 0.01
#define SURFACE_GRADIENT 10.0

// projection parameters
#define PERSPECTIVE 10.0  /* larger: less perspective effect */
#define SCALE 5.0  /* image appears smaller when this is set larger */


// light direction as global variable
vec3 light = normalize(vec3(0.5,0.5,1.0));

// colormaps - https://www.shadertoy.com/view/NsSSRK
vec3 colorSdf(float t) {
    float r = .385+.619*t+.238*cos(4.903*t-2.61);
    float g = -5.491+.959*t+6.089*cos(.968*t-.329);
    float b = 1.107-.734*t+.172*cos(6.07*t-2.741);
    return clamp(vec3(r, g, b), 0.0, 1.0);
}
vec3 colorNormal(float t) {
    float r = .529-.054*t+.55*cos(5.498*t+2.779);
    float g = .21+.512*t+.622*cos(4.817*t-1.552);
    float b = .602-.212*t+.569*cos(5.266*t+2.861);
    return clamp(vec3(r, g, b), 0.0, 1.0);
}


// modified from a volume rendering demo
// https://github.com/harry7557558/Graphics/blob/master/raytracing/webgl_volume/fs-source.glsl
vec3 render(in vec3 ro, in vec3 rd, float t0, float t1) {
    float t = t0;
    vec3 totcol = vec3(0.0);
    float totabs = 1.0;
    float v_old = sdf(ro+rd*t), v;
    float dt = min(STEP, abs(v_old));
    for (float i=ZERO; i<MAX_STEP; i++) {
        t += dt;
        if (t > t1) return totcol;
        v = sdf(ro+rd*t);
        if (v*v_old<0.) break;
        vec3 col = colorSdf(0.5+0.5*sin(ISOSURFACE_FREQUENCY*PI*0.5*(v_old+v)));
        float grad = abs(v-v_old)/dt;
        float grad_abs = (1.0-grad)/dt;
        col = mix(vec3(1,0,0), col, clamp(exp(grad_abs),0.0,1.0));
        float absorb = FIELD_EMISSION+DISCONTINUITY_OPACITY*max(-grad_abs,0.0);
        totabs *= exp(-absorb*dt);
        totcol += col*absorb*totabs*dt;
        dt = min(abs(v_old=v), STEP);
        //if (dt < 1e-3) break;
        dt = max(dt, MIN_STEP);
    }
    if (v*v_old<0.) {
        for (int s = int(ZERO); s < 6; s += 1) {
            v_old = v, dt *= -0.5;
            for (int i = int(ZERO); i < 2; i++) {
                t += dt, v = sdf(ro+rd*t);
                if (v*v_old<0.0) break;
            }
        }
    }
    vec3 grad = sdfGrad(ro+rd*t, 1e-3);
    float grad_col = SURFACE_GRADIENT*(0.5*length(grad)-0.5);
    vec3 col = colorNormal(1.0-1.0/(1.0+exp(2.0*grad_col)));  // 0.5+0.5*tanh(grad_col)
#if COLOR
    col = map(ro+rd*t, true).xyz;
    col *= 0.2+0.05*normalize(grad).y+max(dot(normalize(grad), light),0.0);
#else
    col = 0.2+0.05*normalize(grad).y+col*max(dot(normalize(grad), light),0.0);
#endif
    return totcol + col * totabs;
}


// ray intersection with a box
bool boxIntersection(vec3 ro, vec3 rd, out float tn, out float tf) {
    vec3 inv_rd = 1.0 / rd;
    vec3 n = inv_rd*(ro);
    vec3 k = abs(inv_rd)*BOX_RADIUS;
    vec3 t1 = -n - k, t2 = -n + k;
    tn = max(max(t1.x, t1.y), t1.z);
    tf = min(min(t2.x, t2.y), t2.z);
    if (tn > tf) return false;
    return true;
}

// main
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // set camera
    float rx = iMouse.z!=0. ? 2.0*(3.14*(iMouse.y/iResolution.y)-1.57) : 0.3;
    float rz = iMouse.z!=0. ? -iMouse.x/iResolution.x*4.0*3.14 : -0.6;
    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);
    vec3 ro = SCALE*PERSPECTIVE*w;
    vec2 uv = 2.0*fragCoord.xy/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, PERSPECTIVE*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    light = normalize(w+0.5*u+0.1*v);

    float t0, t1;
    if (!boxIntersection(ro, rd, t0, t1)) {
        fragColor = vec4(vec3(0.0), 1.0);
        return;
    }
    vec3 col = render(ro, rd, t0, t1);;
    fragColor = vec4(col, 1.0);
}
