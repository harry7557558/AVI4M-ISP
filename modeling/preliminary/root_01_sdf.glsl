// SDF visualizer v2
// Allows a usually larger average marching step without missing thin surface
// Expected to be faster than v1

// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926
//#define ZERO min(iTime,0.)
#define ZERO 0.0


float sqr(float x) { return x*x; }
vec2 cossin(float x) { return vec2(cos(x),sin(x)); }

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
    vec2 q = vec2(length(p.xz)-R,p.y);
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
vec2 smind(vec2 ad, vec2 bd, float k) {  // with autodiff
    float h = 0.5+0.5*(bd.x-ad.x)/k;
    if (h<0.0) return bd;
    if (h>1.0) return ad;
  #if 0
    return mix(bd, ad, h) - k*h*(1.0-h);  // approximation
  #else
    float dh = 0.5*(bd.y-ad.y)/k;
    return vec2(
        mix(bd.x,ad.x,h) - k*h*(1.0-h),
        mix(bd.y,ad.y,h)+(ad.x-bd.x)*dh - k*dh*(1.0-2.0*h)
    );
  #endif
    // might be interesting to check this:
    // http://web.archive.org/web/20190830224228/http://www.iquilezles.org/www/articles/smin/smin.htm
}
vec2 smaxd(vec2 ad, vec2 bd, float k) {
    return -smind(-ad, -bd, k);
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
float rand(inout int seed) {
    return hash11(float(seed=seed+1));
}
float randt(inout int seed) {  // pdf(x) = max(1-abs(x),0)
    float t1 = -1.0 + 2.0 * rand(seed);
    float t2 = -1.0 + 2.0 * rand(seed);
    return 0.5*(t1+t2);
}
float randt(inout int seed, float mu, float k) {
    return mu + k * randt(seed);
}


// Modeling

float mapRoot1(vec3 p, int seed) {
    p.x += randt(seed,0.0,0.2)*sin(randt(seed,4.0,4.0)*p.z+2.0*PI*rand(seed));
    p.y += randt(seed,0.0,0.2)*sin(randt(seed,4.0,4.0)*p.z+2.0*PI*rand(seed));
    p.z += 1.0+randt(seed,0.0,0.2)*p.x+randt(seed,0.0,0.2)*p.y;
    if (length(p.xy) > 0.7) return length(p.xy)-0.6;
    if (abs(p.z-1.2)-1.4 > 0.0) return abs(p.z-1.2)-1.3;
    float f = randt(seed,8.0,4.0);
    float a0 = randt(seed,0.0,4.0)*PI;
    float a1 = f*randt(seed,0.9,0.3)*PI;
    float as = sign(rand(seed)-0.5);
    vec2 a = as * smind(smaxd(f*vec2(p.z,1), vec2(a0,0), 5.0), vec2(a1,0), 5.0);
    p.xy = rot2(a.x)*p.xy;
    float w = randt(seed,0.2,0.1) + randt(seed,0.1,0.05)*sin(randt(seed,4.0,4.0)*p.z+2.0*PI*rand(seed));
    float r = randt(seed,0.1,0.03)*exp(-(randt(seed,0.5,0.2)*p.z+randt(seed,0.1,0.1)*sin(randt(seed,3.0,3.0)*p.z+2.0*PI*rand(seed))));
    r = smin(r, 0.5, 0.1);
    float k = length(vec2(length(p.xy)*a.y,1.0));
    float d = sdCapsule((p-vec3(w,0,0))/vec3(1.0,1.0+0.2*k,1.0), 2.5, r);
    d = d / mix(1.0, max(1.0*pow(k,0.7),0.9), smoothstep(1.0,0.0,10.0*(d-r))) - 0.02;
    d = smin(1.1*d, 0.5*(d+0.7), 0.1);
    return d;
}

vec4 mapRoots1(vec3 p, bool col_required) {
    p.y += 0.3*sin(1.8*p.z);
    vec4 c = vec4(1,0,0, 1e12), d;
    for (float i=ZERO-1.; i<=1.; i+=1.) {
        for (float j=ZERO-1.; j<=1.; j+=1.) {
            int seed = int(65536.0*hash12(vec2(i,j)+0.1));
            vec3 q = p;
            q.xy -= vec2(0.5,0.4) * (vec2(i,j) + 0.7 * (hash22(vec2(i,j))-0.5)) * smoothstep(1.0,0.5,-0.5*p.z);
            q.z = randt(seed,1.0,0.3)*q.z - randt(seed,0.0,0.5);
            d.w = mapRoot1(vec3(1,1,-1)*q, seed);
            if (col_required) {
                float t = 0.5+randt(seed,0.5,0.2)*sin(randt(seed,4.0,4.0)*p.z+2.0*PI*rand(seed));
                d.xyz = mix(vec3(0.3,0.2,0.05), vec3(0.75,0.55,0.25), smootherstep(t));
            }
            c = smin(c, d, 0.02);
        }
    }
    return c;
}


vec4 map(vec3 p, bool col_required) {
    vec4 c = mapRoots1((p-vec3(0,0,0.2))/0.8, col_required)*0.8;
    return c;
}


float sdf(vec3 p) {
    const float sc = 1.0;
    return map(p/sc, false).w*sc;
}

vec3 sdfGrad(in vec3 p, in float e) {
    // https://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
  #if 0
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
  #else
    vec3 n = vec3(0.0);
    for(int i=int(ZERO); i<4; i++) {
        vec3 r = -1.0 + 2.0 * vec3(((i+3)>>1)&1, (i>>1)&1, i&1);
        n += r*sdf(p+r*e);
    }
    return n * (.25/e);
  #endif
}



// color surface or not
#define COLOR 1


// raymarching parameters
#define BOX_RADIUS vec3(2.0, 2.0, 2.0)
#define STEP 0.1
#define MIN_STEP 0.005
#define MAX_STEP 200.

// rendering parameters
#define FIELD_EMISSION 0.1
#define ISOSURFACE_FREQUENCY 8.0
#define DISCONTINUITY_OPACITY 0.01
#define SURFACE_GRADIENT 10.0

// projection parameters
#define PERSPECTIVE 10.0  /* larger: less perspective effect */
#define SCALE 4.0  /* image appears smaller when this is set larger */


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
