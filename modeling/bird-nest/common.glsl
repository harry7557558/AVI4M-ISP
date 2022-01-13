#ifndef COMMON_GLSL

#define COMMON_GLSL

// constants
#define PI 3.1415926
#define ZERO min(iTime,0.)

// square
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

// hash functions from https://www.shadertoy.com/view/4djSRW
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

// random functions
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

// noise functions
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

#endif // COMMON_GLSL
