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
vec3 hash33(vec3 p3) {
    p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz + 33.33);
    return fract((vec3(p3.x, p3.x, p3.y) + vec3(p3.y, p3.z, p3.z))*p3.zyx);
}

// noise
float SimplexNoise3D(vec3 xyz) {
    const float K1 = 0.3333333333;
    const float K2 = 0.1666666667;
    vec3 p = xyz + (xyz.x + xyz.y + xyz.z)*K1;
    vec3 i = floor(p);
    vec3 f0 = xyz - (i - (i.x + i.y + i.z)*K2);
    //vec3f e = step(f0.yzx(), f0);  // possibly result in degenerated simplex
    vec3 e = vec3(f0.y > f0.x ? 0.0 : 1.0, f0.z >= f0.y ? 0.0 : 1.0, f0.x > f0.z ? 0.0 : 1.0);
    vec3 i1 = e * (vec3(1.0) - e.zxy);
    vec3 i2 = vec3(1.0) - e.zxy * (vec3(1.0) - e);
    vec3 f1 = f0 - i1 + K2;
    vec3 f2 = f0 - i2 + 2.0*K2;
    vec3 f3 = f0 - 1.0 + 3.0*K2;
    vec3 n0 = 2.0 * hash33(i) - 1.0;
    vec3 n1 = 2.0 * hash33(i + i1) - 1.0;
    vec3 n2 = 2.0 * hash33(i + i2) - 1.0;
    vec3 n3 = 2.0 * hash33(i + 1.0) - 1.0;
    vec4 v = vec4(dot(f0, n0), dot(f1, n1), dot(f2, n2), dot(f3, n3));
    vec4 w = max(-vec4(dot(f0, f0), dot(f1, f1), dot(f2, f2), dot(f3, f3)) + 0.5, vec4(0.0));
    return dot((w*w*w*w) * v, vec4(32.0));
}


// modeling

vec4 mapBerriesFruit(vec3 p, bool col_required) {
    float bound = length(p)-2.2, boundw = 0.5;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    p += vec3(0,0,0.5);
    float r = length(p.xy), a = atan(p.y, p.x);
    float x, y; vec3 q;
    // fruit
    q = vec3(r*cossin(asin(0.97*sin(2.5*a))/2.5), p.z);
    vec4 fruit = vec4(0,0,0, sdEllipsoid(q-vec3(0.4,0,0.68), vec3(0.8+0.1*p.z,1.1,0.95)));
    if (col_required) {
        float noise = SimplexNoise3D(4.0*p);
        fruit.xyz = mix(vec3(0.5,0.05,0.1),vec3(0.8,0.1,0.05),
            smootherstep(0.6*(q.z+1.0*(r-0.9)-0.5+smax(q.z-1.2,0.0,0.1)))
        ) + vec3(0.15)*(-noise+0.5);
        fruit.xyz = mix(fruit.xyz, vec3(0.8,0.8,0.0), 0.25+0.2*p.z);
        fruit.xyz = mix(fruit.xyz, vec3(0.8,0.0,0.5), 0.2);
    }
    q = vec3(r*cossin(asin(0.95*sin(2.5*a+0.8))/2.5), p.z);
    fruit = smin(fruit, vec4(  // hair
        mix(mix(vec3(0.8,0.75,0.0),vec3(0.9,0.85,0.5),smootherstep(r/0.3)), vec3(0.8,0.4,0.0), 0.2),
        sdEllipsoid(roty(0.2-0.05*cos(3.0*a))*(q-vec3(0.08,0,1.61)), vec3(0.2+0.03*sin(4.0*a),0.12,0.05))), 0.05);
    // sepal/stem
    q = vec3(r*cossin(asin(0.98*sin(2.5*a-1.2))/2.5), p.z);
    q -= vec3(0.2,0,0);
    vec3 br = 1.2*vec3(
        1.0+0.1*sin(3.0*a)-0.1*cos(2.0*a),
        1.25*(0.3-0.2*exp(-sqr(1.8*(q.x-1.6)))-0.2*exp(-sqr(3.0*(q.x-0.0)))),
        0.05);
    q = roty(smin(0.2*r,0.5,0.1))*(q-1.4*vec3(0.2,0,-0.15));
    q.z -= clamp(
        (0.1*pow(abs(q.x),4.0)+0.02*r*r*sin(3.0*a))*sin(a) + 0.02*sin(12.0*r)-0.1*sin(4.0*pow(r,1.3))
        - 0.05*exp(-sqr(12.0*q.y))*exp(-r),
        -0.5, 0.5);
    vec4 disk = vec4(mix(0.8*vec3(0.35,0.5,0.0),vec3(0.4,0.55,0.0),smootherstep(2.0*r-1.0)), sdEllipsoid(q, br));  // sepals
    q = p;
    vec4 stem = vec4(0.4,0.4,0.15, sdCapsule(q-vec3(0.05,0,-0.5), 0.2, 0.1));  // stem
    disk = smin(disk, stem, 0.1);
    vec4 d = smin(disk, fruit, 0.03);
    return d;
}

vec4 mapBerriesLeaf(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(-0.1,0.0,0.1), vec3(2.2,1.5,0.6)), boundw = 0.3;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    vec3 br;
    br.x = 1.6;
    br.y = 0.85 -0.45*exp(-sqr(1.2*(p.x-1.5))) -0.3*exp(-sqr(2.5*(p.x+1.6)));
    br.z = 0.08;
    float dqz = 0.0;
    float d2e = abs( br.y*sqrt(smax(br.x*br.x-p.x*p.x,0.0,0.1))/br.x - abs(p.y) );
    float veins = asin(0.9*sin((16.0+2.0*p.x)*(p.x-0.5*abs(p.y)-0.5*pow(abs(p.y),1.3)+0.2*pow(d2e,0.4))+sign(p.y)*0.4*PI));
    float veins_fade = (1.0-exp(-(6.0/br.y)*abs(p.y))) * (1.0-exp(-(2.0/br.y)*d2e)) * (exp(-0.2*p.x));
    dqz += 0.01 * veins*veins_fade;
    dqz += 0.05 * smax(veins_fade*(veins-0.9),0.0,0.1);
    float midrib = abs(p.y)-1.0;
    float midrib_fade = (1.0-exp(-(4.0/br.y)*abs(p.y))) * (exp(-0.4*p.x));
    dqz += 0.05 * midrib*midrib_fade;
    vec3 q = p;
    q.z -= 0.5*(1.0-exp(-0.5*length(vec2(p.y,0.01))));
    q.z -= 0.1*cos(p.x);
    vec4 leaf = vec4(0,0,0, sdEllipsoid(q+vec3(0,0,dqz), br));
    if (col_required) leaf.xyz = mix(vec3(0.3,0.45,0.05), vec3(0.55,0.7,0.15), -0.5-20.0*dqz);
    vec4 stem = vec4(0.3,0.35,0.05, sdSegment(q-vec3(0,0,-0.00), vec3(-1.8,0,0), vec3(1.0,0,0)) - 0.06*exp(-0.2*abs(p.x+1.8)));
    return smin(leaf, stem, mix(0.05,0.001,clamp(p.x+1.8,0.,1.)));
}

vec4 mapBerries1(vec3 p, bool col_required) {
    vec4 res = vec4(1, 0, 0, 1e8);
    return mapBerriesFruit(p, col_required);
}


vec4 map(vec3 p, bool col_required) {
    //return mapBerriesLeaf(p, col_required);
    //return mapBerriesFruit(p, col_required);
    return mapBerries1(p, col_required);
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
