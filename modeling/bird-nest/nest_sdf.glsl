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

// rotation matrices
mat3 rotx(float a) {
    return mat3(
        1, 0, 0,
        0, cos(a), sin(a),
        0, -sin(a), cos(a)
    );
}
mat3 rotz(float a) {
    return mat3(
        cos(a), sin(a), 0,
        -sin(a), cos(a), 0,
        0, 0, 1
    );
}

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


// random/hash functions from https://www.shadertoy.com/view/4djSRW
float hash11(float p) {
    p = fract(p * .1031);
    p *= p + 33.33;
    return fract(p * (p + p));
}
float hash12(vec2 p) {
    vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}
vec2 hash22(vec2 p) {
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

// hash-based random functions
float rand(inout float seed) {
    return hash11(float(seed=seed+0.1));
}
float randt(inout float seed) {  // pdf(x) = max(1-abs(x),0)
    float t1 = -1.0 + 2.0 * rand(seed);
    float t2 = -1.0 + 2.0 * rand(seed);
    return 0.5*(t1+t2);
}
float randt(inout float seed, float mu, float k) {
    return mu + k * randt(seed);
}

// quasi-random
float vanDerCorput(float n, float b) {
    float x = 0.0;
    float e = 1.0 / b;
    while (n > 0.0) {
        float d = mod(n, b);
        x += d * e;
        e /= b;
        n = floor(n / b);
    }
    return x;
}


// modeling
vec4 mapNest(vec3 p, bool col_required) {
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
