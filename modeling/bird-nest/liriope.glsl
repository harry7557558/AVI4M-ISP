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
vec2 cossin(float x) { return vec2(cos(x), sin(x)); }
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

vec4 mapLiriopeFlowersLayer(float elev, vec3 p, bool col_required) {
    float r = length(p.xy), a = atan(p.y, p.x);
    a += 0.05*sin(4.0*a) + 0.05*cos(a);
    vec3 q = vec3(r*cossin(asin(0.97*sin(2.5*a))/2.5), p.z);
    q = roty(elev)*q;
    vec4 stem = vec4(0.9,0.85,0.95, sdSegment(q, vec3(0.0), vec3(1.0,0,0))-0.1);
    float r1 = length(q.yz), a1 = atan(q.z, q.y);
    vec3 q1 = vec3(q.x, r1*cossin(asin(0.93*cos(2.5*a1))/2.5));
    q1 = roty(0.3*PI)*(q1-vec3(1,0,0))+vec3(1,0,0);
    vec4 flower = vec4(mix(vec3(0.65,0.45,0.7), vec3(0.9,0.8,0.95), (3.0*(q1.x-1.0)+0.5)+0.15*sin(5.0*a1)),
        sdEllipsoid(q1-vec3(1.0,0.0,0.05), vec3(0.45,vec2(0.3+0.1*(q1.x-1.0)))));
    vec4 d = smin(stem, flower, 0.1);
    return d;
}

vec4 mapLiriopeFlowers(vec3 p, bool col_required) {
    p.x += 0.1*cos(p.z);
    float bound = sdSegment(p,vec3(0,0,-2.8),vec3(0,0,2.0))-1.0, boundw = 0.4;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    vec4 stem = vec4(0,0,0, sdSegment(p, vec3(0,0,-3.0), vec3(0,0,2))-0.06);
    if (col_required) {
        stem.xyz = mix(vec3(0.6,0.35,0.65), vec3(0.9,0.8,0.95), smootherstep((p.z+1.0)/3.0));
        stem.xyz *= mix(vec3(0.35,0.2,0.25), vec3(1.0), 0.2+0.8*smootherstep((p.z+2.5)/2.5));
    }
    float seed = 0.0;
    for (float t=0.0; t<1.0; t+=1.0/11.0) {
        float h = mix(-1.5, 2.0, 1.0-pow(1.0-t,1.2));
        h += 0.5*t*(1.0-t)*(2.0*rand(seed)-1.0);
        vec3 q = p-vec3(0,0,h);
        q.xy = rot2(2.0*PI*rand(seed))*q.xy;
        float elev = 0.05*PI + 0.15*PI*t + mix(0.05*PI, 0.2*PI, rand(seed));
        float sc = 0.08*(1.0+0.3*t-1.3*t*t) + mix(0.25, 0.4, smoothstep(0.,1.,(p.x+1.0)/2.0)) * (1.0+0.1*rand(seed));
        vec4 d = mapLiriopeFlowersLayer(elev, q/sc, col_required)*vec4(1,1,1,sc);
        if (col_required) d.xyz = mix(d.xyz, mix(vec3(0.65,0.45,0.75), vec3(0.8,0.7,0.9), t), 0.2);
        stem = smin(stem, d, 0.05);
    }
    return stem;
}

vec4 mapLiriopeLeaf(vec3 p, bool col_required) {
    p.z -= 0.4*cos(0.25*PI*p.x);
    float bound = sdSegment(p,vec3(-2.0,0,0),vec3(2.0,0,0))-0.5, boundw = 0.3;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    float near_stem = 1.0-0.3/(1.0+pow(abs(0.6*(p.x+2.0)),4.0));
    float w = pow(max(1.0-sqr(p.x/2.0), 0.0), 0.5) * (0.2/(1.0+sqr(0.3*(p.x-0.5)))) * near_stem;
    float u = clamp(p.y/w, -1.0, 1.0);
    float thickness = 0.2*w * pow(max(1.0-u*u,0.0),0.5) * (1.0+sqr(p.y/0.3)) * (exp(-0.1*(p.x+2.0))) / pow(near_stem, 3.5);
    float veins = 0.03*cos(15.0*u)*(1.0-u*u);
    float zd = 0.05/sqrt(1.0+sqr(p.y/0.1));
    vec4 leaf = vec4(0,0,0, sdSegment(p.yz+vec2(0,zd), vec2(-w,0), vec2(w,0)) - thickness * (1.0+veins));
    if (col_required) {
        leaf.xyz = pow(mix(vec3(0.35,0.55,0.25), vec3(0.65,0.8,0.5), (p.x+2.0)/4.0), vec3(1.8));
        leaf.xyz = mix(leaf.xyz, vec3(0.45,0.65,0.3), 1.0-20.0*zd) * vec3(1.0+1.0*veins);
    }
    return leaf;
}


vec4 map(vec3 p, bool col_required) {
    //return mapLiriopeFlowersLayer(0.2*PI, p, col_required);
    //return mapLiriopeFlowers(p/0.8, col_required)*vec4(1,1,1,0.8);
    return mapLiriopeLeaf(p, col_required);
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
#define BOX_RADIUS vec3(2.5, 2.5, 2.5)
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
