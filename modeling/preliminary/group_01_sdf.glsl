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
    vec2 q = vec2(length(p.xy)-R,p.z);
    return length(q)-r;
}
float sdCapsuleY(vec3 p, float h, float r) {
    p.y = abs(p.y)-min(abs(p.y), h);
    return length(p) - r;
}
float sdTorusY(vec3 p, float R, float r) {
    vec2 q = vec2(length(p.xz)-R,p.y);
    return length(q)-r;
}
float sdBox(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    return length(max(q,vec3(0))) + min(max(q.x,max(q.y,q.z)),0.0);
}
float sdLnNormEllipsoid(vec3 p, vec3 r, float n) {
    p = abs(p);
    float k1 = pow(dot(pow(p/r,vec3(n)),vec3(1)),1.0/n);
    float k2 = pow(dot(pow(p/(r*r),vec3(n)),vec3(1)),1.0/n);
    return k1*(k1-1.0)/k2;
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
vec4 smax(vec4 a, vec4 b, float k) {
    return smin(a*vec4(1,1,1,-1), b*vec4(1,1,1,-1), k) * vec4(1,1,1,-1);
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



// Modeling

vec4 mapFlower01(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(0,0,0.25), vec3(2.5,2.5,1.5));
    if (bound > 0.0) return vec4(1,0,0, bound+0.1);
    p += 0.08*sin(4.0*p.yzx)*cos(4.0*p.xyz)*sin(4.0*p.zxy);
    float x, y; vec3 q;
    float r = length(p.xy), a = atan(p.y, p.x);
    float a_ = asin(0.99*sin(2.5*a))/2.5;
    vec3 p_ = vec3(r*vec2(cos(a_), sin(a_)), p.z);
    // petals
    q = p_;
    x = q.x-1.0, y = q.y;
    q.x -= 1.2;
    q.z += 0.1*exp(-sqr(3.0*y)) - 0.1*length(vec2(r,0.2)) - 0.4/(1.0+exp(-4.0*(r-1.6)));
    float w = 1.05*(0.2
        + 0.3*exp(-sqr(1.2*(x-0.2)))
        + 0.1*exp(-sqr(1.5*(x+0.8)))
        - 0.05*exp(-sqr(4.0*(x-1.0)))
        + (0.1+0.1*cos(4.0*a)*sin(5.0*a))*exp(-sqr(1.0*(x-1.0))));
    float h = 0.2*w+0.1*exp(-sqr(3.0*(x+0.8)));
    float petal_d = sdEllipsoid(q, vec3(1.0, w, h));
    vec3 petal_c = vec3(0.0);
    if (col_required) {
        float t = 0.5+0.5*cos(10.0*a_)-exp(-4.0*(r-0.5));
        petal_c = mix(vec3(0.8,0.55,0.65), vec3(0.85,0.7,0.8), smoothstep(0.,1.,t));
        t = sdSegment(q.xy, vec2(-1.0,0.0), vec2(0.5,0.0))-0.05*exp(-2.0*r);
        petal_c = mix(vec3(0.9,0.25,0.25), petal_c, 0.7+0.3*smootherstep(8.0*t));
    }
    vec4 petal = vec4(petal_c, petal_d);
    // style/sepal
    q = vec3(r, a_, p.z);
    float leaf_d = smin(  // green part
        sdSegment(q.xz, vec2(0.24,-0.05), vec2(0.32,0.0))-0.05,  // ring
        sdSegment(q.xz+0.05*sin(8.0*q.z), vec2(0.12,-0.4), vec2(0.22,-0.05))-0.05,  // wall
        0.05);
    leaf_d = smin(leaf_d, sdEllipsoid(p_-vec3(0.55,0,-0.1), vec3(0.4,0.1,0.05)), 0.05);  // sepal
    leaf_d = smin(leaf_d,
        sdSegment(p-vec3(0.01*sin(10.0*r),0,0), vec3(0,0,-0.6),vec3(0,0,1.0))
            -max(0.06/(1.0+exp(10.0*(q.z-0.5)))+0.02,0.),  // style
        0.1);
    vec4 leaf = vec4(0.5,0.7,0, leaf_d);
    // filament/anther
    a_ = asin(0.99*sin(7.5*atan(p_.y,p_.x)))/7.5;
    p_ = vec3(r*vec2(cos(a_), sin(a_)), p.z);
    q = p_ - vec3(0.1+(0.7+0.2*sin(7.0*a))/(1.0+exp(-8.0*(p.z-0.2))),0,-0.2);
    h = 0.8+0.2*sin(4.0*a);
    vec4 filament = vec4(0.85,0.85,0.,
        0.5*sdCapsule(q, h, 0.01+0.01*exp(-sqr(4.0*(p.z-0.2)))));
    vec4 anther = vec4(0.75,0.55,0.,
        sdEllipsoid(q-vec3(0.02,0,h), vec3(0.08,0.04,0.05)));
    // put all together
    vec4 d = cmin(smin(petal, leaf, 0.05-0.03*exp(-r)), smin(filament, anther, 0.01)-0.005);
    return d;
}

vec4 mapPetal02(vec3 p, float sc, float curve, bool col_required) {
    p /= sc;
    p.z -= min(curve*(sqr(p.x)+sqr(p.y)),1.0);
    vec3 pr = vec3(1.0, 0.5, 0.1);
    pr.y = 0.7 - 0.0*exp(-sqr(2.0*(p.x+0.8))) + 0.6*exp(-sqr(1.2*(p.x-0.6)));
    float d = sdEllipsoid(p, 0.9*pr);
    d *= exp(-sqr(curve));
    vec3 col = mix(vec3(0.9,0.6,0.55), vec3(0.98,0.95,0.97), 0.5+0.5*tanh(p.x));
    return vec4(col, d*sc);
}
vec4 mapFilament02(vec3 p, float h, bool col_required) {
    vec4 fil = vec4(1.0,0.8,0.25, sdCapsule(p, h, 0.015));
    vec4 ant = vec4(1.0,0.6,0.0, sdEllipsoid(p-vec3(0,0,h), vec3(0.04,0.03,0.02)));
    vec4 d = smin(fil, ant, 0.01);
    return d;
}
vec4 mapFlower02(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(0,0,0.3), vec3(2.0,2.0,1.5));
    if (bound > 0.0) return vec4(1,0,0, bound+0.1);
    p += 0.04*sin(4.0*p.yzx)*cos(4.0*p.xyz)*sin(4.0*p.zxy);
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q;
    // petals
    q = vec3(r*cossin(asin(0.999*sin(2.5*a-1.3))/2.5), p.z);
    vec4 p0 = mapPetal02(roty(0.2+0.15*sin(3.0*a))*(q-vec3(1.0,0,-0.2)), 0.5, 0.3, col_required);
    q = vec3(r*cossin(asin(0.999*sin(2.5*a))/2.5), p.z);
    vec4 p1 = mapPetal02(roty(0.75+0.2*cos(4.0*a))*(q-vec3(1.2,0,0)), 0.7, 0.5, col_required);
    q = vec3(r*cossin(asin(0.999*sin(2.5*a+0.8))/2.5), p.z);
    vec4 p2 = mapPetal02(roty(0.85-0.2*sin(5.0*a-1.0))*(q-vec3(1.15,0,0.25+0.05*sin(a+1.0))), 0.65, 0.45, col_required);
    q = vec3(r*cossin(asin(0.999*sin(2.5*a+1.5))/2.5), p.z);
    vec4 p3 = mapPetal02(roty(1.05+0.1*sin(4.0*a+1.0))*(q-vec3(1.0,0,0.45+0.1*cos(a))), 0.55, 0.3, col_required);
    q = vec3(r*cossin(asin(0.999*sin(2.5*a+2.7))/2.5), p.z);
    vec4 p4 = mapPetal02(roty(1.0+0.1*cos(3.0*a))*(q-vec3(0.85-0.05*sin(4.0*a),0,0.35-0.05*sin(a))), 0.45, 0.55, col_required);
    const float petal_k = 0.02;
    vec4 petal = smin(p0, smin(smin(p1, p2, petal_k), smin(p3, p4, petal_k), petal_k), petal_k);
    if (col_required) petal.xyz *= mix(vec3(0.75,0.25,0.2), vec3(1.0,0.9,0.85), smootherstep(1.3-0.8*exp(-0.5*p.z)));
    // filament
    q = vec3(r*cossin(asin(0.999*sin(11.0*a-1.3))/11.0), p.z);
    vec4 f1 = mapFilament02(roty(-0.4)*(q-vec3(0.25+0.15/(1.0+exp(-4.0*(p.z-0.8)))+0.05*cos(5.0*a),0,0)), 1.0+0.1*sin(7.0*a), col_required);
    q = vec3(r*cossin(asin(0.999*sin(13.0*a-0.8))/13.0), p.z);
    vec4 f2 = mapFilament02(roty(-0.2)*(q-vec3(0.25+0.15/(1.0+exp(-4.0*(p.z-0.8)))+0.05*cos(6.0*a+1.2),0,0)), 1.0+0.1*sin(4.0*a-0.3), col_required);
    vec4 filament = 0.9*smin(f1, f2, 0.01);
    // disk
    q = vec3(r*cossin(asin(0.98*sin(2.5*a))/2.5), p.z);
    vec3 br = vec3(0.6, 0.4-0.3*exp(-sqr(2.0*(q.x-1.0))), 0.1);
    vec4 d1 = vec4(0.75,0.55,0.0, sdEllipsoid(roty(0.1)*(q-vec3(0.5,0,-0.25)), br));
    q = vec3(r*cossin(asin(0.9*sin(2.5*a))/2.5), p.z);
    vec4 d2 = vec4(0.95,0.7,0.0, sdEllipsoid(q-vec3(0.1,0,-0.1), vec3(0.3,0.2,0.2)));
    vec4 d3 = vec4(0.65,0.6,0.05, sdSegment(q-vec3(0.05,0,0), vec3(0,0,-0.55),vec3(0,0,-0.2))-0.08);
    vec4 disk = smin(smin(d1, d2, 0.1), d3, 0.05);
    // put them together
    vec4 d = smin(cmin(petal, filament), disk, 0.05);
    return d;
}

vec4 mapPetal03(vec3 p, float sc, float curve, bool col_required) {
    p /= sc;
    p.z -= curve*min(dot(p.xy,p.xy), 1.0);
    p.z += 0.1/(1.0+sqr(2.0*p.y));
    vec3 pr = vec3(1.0, 0.5, 0.1);
    pr.y = 0.3 + 0.35*exp(-sqr(1.2*(p.x+0.4)));
    float d = sdEllipsoid(p, pr);
    if (dot(p.xy,p.xy)<1.0) d /= length(vec3(curve*2.0*p.xy,1.0));
    vec3 col = mix(vec3(0.6,0.55,0.05), vec3(0.85,0.9,0.3), 0.5+0.5*tanh(2.0*p.x));
    return vec4(col, d*sc);
}
vec4 mapFilament03(vec3 p, float h, bool col_required) {
    vec4 fil = vec4(0.95,0.95,0.9, sdCapsule(p, h, 0.024));
    vec4 ant = vec4(0.98,0.95,0.75, sdEllipsoid(p-vec3(0,0,h), vec3(0.08,0.07,0.06)));
    vec4 d = smin(fil, ant, 0.01);
    return d;
}
vec4 mapFlower03(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(0,0,0.1), vec3(2.0,2.0,1.5));
    if (bound > 0.0) return vec4(1,0,0, bound+0.1);
    p /= vec3(1.0, 1.0, 1.1);
    p.y += 0.02*sin(8.0*p.z);
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q;
    // petals
    q = vec3(r*cossin(asin(0.999*sin(3.5*a-0.2*cos(4.0*r)*sin(5.0*a)+0.9))/3.5), p.z);
    vec4 p0 = mapPetal03(roty(-0.4)*(q-vec3(0.8,0,-0.2)), 0.45, -0.1, col_required);
    q = vec3(r*cossin(asin(0.999*sin(3.5*a-0.2*sin(4.0*r)*cos(5.0*a)+0.0))/3.5), p.z);
    vec4 p1 = mapPetal03(roty(-0.1)*(q-vec3(1.0,0,-0.05)), 0.5, -0.2, col_required);
    q = vec3(r*cossin(asin(0.999*sin(3.0*a-0.2*cos(4.0*r)*cos(3.0*a)-0.9))/3.0), p.z);
    vec4 p2 = mapPetal03(rotx(0.1)*roty(0.1)*(q-vec3(1.0,0,0.1)), 0.55, -0.1, col_required);
    q = vec3(r*cossin(asin(0.999*sin(3.0*a-0.2*cos(4.0*r)*cos(4.0*a)-2.0))/3.0), p.z);
    vec4 p3 = mapPetal03(roty(0.4)*(q-vec3(0.7,0,0.16)), 0.55, -0.05, col_required);
    q = vec3(r*cossin(asin(0.999*sin(2.5*a-0.2*sin(3.0*r)*sin(3.0*a)-2.7))/2.5), p.z);
    vec4 p4 = mapPetal03(roty(0.6)*(q-vec3(0.6,0,0.2)), 0.45, 0.05, col_required);
    q = vec3(r*cossin(asin(0.999*sin(2.0*a-0.3*sin(3.0*r)*sin(3.0*a)-3.5))/2.0), p.z);
    vec4 p5 = mapPetal03(roty(0.7)*(q-vec3(0.4,0,0.15)), 0.4, 0.2, col_required);
    vec4 petal = smin(smin(smin(p0, p1, 0.01), smin(p2, p3, 0.01), 0.01), smin(p4, p5, 0.01), 0.01);
    // filaments
    q = vec3(r*cossin(asin(0.999*sin(3.5*a+0.9))/3.5), p.z);
    vec4 f1 = mapFilament03(roty(-0.7-0.1*cos(4.0*a))*(q-vec3(0.2,0,0.0)), 0.8, col_required);
    q = vec3(r*cossin(asin(0.999*sin(4.0*a+1.5))/4.0), p.z);
    vec4 f2 = mapFilament03(roty(-0.8+0.1*sin(3.0*a))*(q-vec3(0.18,0,0.0)), 0.8, col_required);
    vec4 filament = smin(f1, f2, 0.01);
    // disk
    q = vec3(r*cossin(asin(0.98*sin(3.5*a-0.5))/3.5), p.z);
    vec3 br = vec3(0.4, 0.3-0.2*exp(-sqr(2.0*(q.x-1.0))), 0.05);
    vec4 d1 = vec4(0.4,0.35,0.0, sdEllipsoid(roty(0.2)*(q-vec3(0.2,0,-0.25)), br));  // sepals
    q = vec3(r*cossin(asin(0.8*sin(3.5*a))/3.5), p.z);
    vec4 d2 = vec4(0.6,0.5,0.0, sdEllipsoid(q-vec3(0.05,0,-0.1), vec3(0.25,0.2,0.3)));  // ovary
    vec4 d3 = vec4(0.9,0.8,0.0, sdCapsule(q-vec3(0.05,0,0),0.56,0.08-0.03*sin(2.0*p.z)));  // style
    vec4 d4 = vec4(0.95,0.95,0.6, length(q-vec3(0.01,0,0.56))-0.06);  // stigma
    vec4 disk = smin(smin(d1, d2, 0.1), smin(d3, d4, 0.05), 0.1);
    // put them together
    vec4 d = smin(cmin(petal, filament), disk, 0.01);
    return d;
}

vec4 mapFruit01(vec3 p, bool col_required) {
    p *= vec3(1.0,1.05,1.0);
    float bound = sdEllipsoid(p-vec3(0,0,0.3), vec3(1.4,1.4,1.8));
    if (bound > 0.0) return vec4(1,0,0, bound+0.1);
    float r = length(p.xy), a = atan(p.y, p.x), b = atan(r, p.z);
    vec3 q;
    q = vec3(length(vec2(r,0.02))*cossin(asin(0.9*sin(2.5*a-0.5))/2.5), p.z);
    float bottom_hole = -smin(
            sdLnNormEllipsoid(q-vec3(0,0,-1.0), vec3(0.3,0.3,0.3), 5.0),
            sdTorus(q-vec3(0,0,-0.7), 0.04, 0.05),
            0.05);
    vec4 d = vec4(0.05,0.7,0.9, smax(
        sdEllipsoid(q, vec3(1.15,1.15,1.0-0.4/(1.0+10.0*length(vec3(q.xy,0.2)))))
            + 0.02*sin(5.0*p.x)*sin(5.0*p.y)*sin(5.0*p.z),
        bottom_hole, 0.02));
    if (col_required) {
        const vec3 c1 = vec3(0.15,0.25,0.4);
        const vec3 c2 = vec3(0.3,0.45,0.6);
        const vec3 c3 = vec3(0.8,0.75,0.85);
        float t = dot(p, vec3(0.5,0.5,1)) + 0.8*SimplexNoise3D(1.5*p);
        d.xyz = mix(mix(c1, c2, 0.5+0.5*tanh(1.0*t)), c3, 0.5+0.5*tanh(1.0*(t-1.0)));
        t = 0.4*GradientNoise2D(vec2(6.0*a,3.0*b))-4.0*(b/PI-0.5)+0.5;
        d.xyz *= smootherstep(0.7+0.3*tanh(1.2*t));
    }
    d = smin(d, vec4(0.5*vec3(0.1,0.12,0.15),
        sdTorus(q-vec3(0,0,-0.92), 0.33, 0.05)),
        0.05);
    q = vec3(length(vec2(r,0.02))*cossin(asin(0.999*sin(2.5*a-0.5))/2.5), p.z);
    d = smin(d, vec4(0.5*vec3(0.25,0.3,0.4),
        smax(
            sdLnNormEllipsoid(roty(-1.1)*(q-vec3(0.2+0.2/(1.0+10.0*length(vec2(q.y,0.05))),0,-0.8)), vec3(0.3,0.25,0.05), 2.5),
            bottom_hole, 0.02)),
        0.08);
    if (col_required) {
        d.xyz *= 0.2+0.8*smootherstep(10.0*sdTorus(p-vec3(0,0,-0.78),0.18,0.08));
        d.xyz = tanh(d.xyz * vec3(0.7,1.1,1.3));
    }
    d = smin(d, vec4(0.75,0.75,0.5,
        sdCapsule(q-vec3(0,0,0.8), 0.8, 0.2+0.05*sin(4.0*q.z))), 0.05);
    return d;
}

vec4 mapFruit02(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(0,0,0.1), vec3(1.8,1.8,1.5));
    if (bound > 0.0) return vec4(1,0,0, bound+0.1);
    p *= vec3(1.05, 1.0, 1.0);
    p += vec3(0,0,0.5);
    float r = length(p.xy), a = atan(p.y, p.x);
    vec3 q;
    // fruit
    q = vec3(r*cossin(asin(0.99*sin(2.5*a+0.3*sin(4.0*a)))/2.5), p.z);
    vec4 fruit = vec4(0,0,0, sdEllipsoid(q-vec3(0.4,0,0.68), vec3(1.0,0.8,0.85+0.05*r)));
    if (col_required) {
        fruit.xyz = mix(vec3(0.55,0.2,0.0),vec3(0.95,0.85,0.0),
            smootherstep(0.4*(q.z+1.0*(r-0.9)-0.0+smax(q.z-1.2,0.0,0.1))
                +0.012*sin(20.0*p.x)*sin(20.0*p.y)*sin(20.0*p.z))
        );
    }
    q = vec3(r*cossin(asin(0.99*sin(2.5*a+0.8))/2.5), p.z);
    fruit = smin(fruit, vec4(  // hair
        mix(vec3(0.1,0.02,0.06),vec3(0.25,0.02,0.0),smootherstep(r/0.3)),
        sdEllipsoid(roty(0.2-0.05*cos(3.0*a))*(q-vec3(0.08,0,1.54)), vec3(0.2+0.03*sin(4.0*a),0.12,0.05))), 0.05);
    // sepal/stem
    q = vec3(r*cossin(asin(0.98*sin(2.5*a-1.2))/2.5), p.z);
    q -= vec3(0.2,0,0);
    vec3 br;
    br.x = 1.0+0.1*sin(3.0*a);
    br.y = 1.2*(0.3-0.2*exp(-sqr(1.8*(q.x-1.0)))+0.2*exp(-sqr(4.0*(q.x-0.7))));
    br.z = 0.1*pow(smax(br.y,0.0,0.1),0.2);
    br *= 1.2;
    q = roty(0.3)*(q-vec3(0.2,0,-0.25));
    q.z += 0.05*sin(6.0*q.x);
    vec4 d1 = vec4(0.4,0.35,0.0, sdEllipsoid(q, br));  // sepals
    q = vec3(r*cossin(asin(0.95*sin(2.5*a-0.9))/2.5), p.z);
    vec4 d2 = vec4(0.25,0.15,0.05, sdCapsule(q-vec3(0.05,0,-0.6), 1.0, 0.1));  // stem
    vec4 disk = smin(d1, d2, 0.1);
    // put them together
    vec4 d = smin(disk, fruit, 0.03);
    return d;
}

vec4 mapSwirl01(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(0,0,0.0), vec3(0.8,0.8,2.0));
    if (bound > 0.0f) return vec4(1,0,0, bound+0.1);
    p.z += 1.0+0.1*p.y;
    vec2 a = -smind(smaxd(8.0*vec2(p.z,1), vec2(0.0,0), 3.0), vec2(4.0*PI,0), 5.0);
    p.xy = rot2(a.x)*p.xy;
    float w = 0.25 + 0.3*exp(-sqr(6.0*p.z));
    float r = 0.1*exp(-1.5*(0.2*p.z+0.1*sin(16.0*p.z)));
    float k = length(vec2(length(p.xy)*a.y,1.0));
    float d = sdCapsule((p-vec3(w,0,0))/vec3(1.0,1.0+0.3*k,1.0), 2.5, r);
    d = d / mix(1.0, max(1.0*pow(k,0.7),0.9), smoothstep(1.0,0.0,10.0*(d-r))) - 0.01;
    d = smin(1.1*d, 0.5*(d+0.7), 0.1);
    return vec4(
        mix(vec3(0.6,0.65,0.0), vec3(0.45,0.4,0.0), 0.5+0.5*sin(10.0*p.x)*sin(10.0*p.y)*sin(10.0*p.z)),
        d);
}
vec4 mapSwirl02(vec3 p, bool col_required) {
    p.z += 1.0-0.1*p.y;
    float bound = sdEllipsoid(p-vec3(0,0,1.2), vec3(0.8,0.8,2.0));
    if (bound > 0.0f) return vec4(1,0,0, bound+0.1);
    vec2 a = -smind(smaxd(6.0*vec2(p.z-0.5,1), vec2(0.0,0), 2.0), vec2(3.0*PI,0), 5.0);
    p.xy = rot2(a.x)*p.xy;
    float w = 0.25 + 0.15*exp(-sqr(6.0*p.z));
    float r = 0.1*exp(-1.5*(0.2*p.z+0.1*sin(16.0*p.z)));
    float k = length(vec2(length(p.xy)*a.y,1.0));
    float d = sdCapsule((p-vec3(w,0,0))/vec3(1.0,1.0+0.3*k,1.0), 2.5, r);
    d = d / mix(1.0, max(1.0*pow(k,0.7),0.9), smoothstep(1.0,0.0,10.0*(d-r))) - 0.005;
    d = smin(1.1*d, 0.5*(d+0.7), 0.1);
    return vec4(
        mix(vec3(0.6,0.65,0.0), vec3(0.45,0.4,0.0), 0.5+0.5*sin(10.0*p.x)*sin(10.0*p.y)*sin(10.0*p.z)),
        d);
}
vec4 mapSwirl04(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(0.3,-0.2,0.2), vec3(1.8,1.6,2.0));
    if (bound > 0.0f) return vec4(1,0,0, bound+0.1);
    p.xz = rot2(0.5)*p.xz;
    p.z -= 1.0;
    p.yz = mix(vec2(p.y,p.z), vec2(p.z,-1.6*p.y-0.4*p.z), smootherstep(0.5*(p.z-p.y)+0.5));
    p.z += 2.0;
    p.xz = rot2(-0.8*exp(-sqr(p.z-0.0)))*p.xz;
    vec2 a = -smind(smaxd(9.0*vec2(p.z,1.0), vec2(2.0*PI,0), 3.0), vec2(10.0*PI,0), 5.0);
    p.xy = rot2(a.x)*p.xy;
    float k = length(vec2(length(p.xy)*a.y,1.0));
    float w = 0.25 + 0.01*k + 0.15*exp(-sqr(1.0*(p.z-3.0)));
    float r = 0.1 - 0.03*exp(-sqr(1.5*(p.z+0.0)));
    float d = sdCapsule((p-vec3(w,0,-0.5))/vec3(1.1,1.2+0.3*k,1.0), 4.5, r);
    d = d / mix(1.0, max(1.0*pow(k,0.9),0.9), smoothstep(1.0,0.0,5.0*(d-r))) - 0.0;
    d = smin(1.1*d, 0.5*(d+0.7), 0.1);
    return vec4(
        mix(vec3(0.6,0.65,0.0), vec3(0.45,0.4,0.0), 0.5+0.5*sin(10.0*p.x)*sin(10.0*p.y)*sin(10.0*p.z)),
        d);
}

float mapRoot01(vec3 p, int seed) {
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
    float r = randt(seed,0.1,0.03)*exp(-(randt(seed,0.3,0.15)*p.z+randt(seed,0.1,0.1)*sin(randt(seed,3.0,3.0)*p.z+2.0*PI*rand(seed))));
    r = smin(r, 0.5, 0.1);
    float k = length(vec2(length(p.xy)*a.y,1.0));
    float d = sdCapsule((p-vec3(w,0,0))/vec3(1.0,1.0+0.2*k,1.0), 2.5, r);
    d = d / mix(1.0, max(1.0*pow(k,0.7),0.9), smoothstep(1.0,0.0,10.0*(d-r))) - 0.02;
    d = smin(1.1*d, 0.5*(d+0.7), 0.1);
    return d;
}
vec4 mapRoots01(vec3 p, bool col_required) {
    p.y += 0.3*sin(1.8*p.z);
    float bound = sdEllipsoid(p-vec3(0,0,-0.3), vec3(1.6,1.6,2.5));
    if (bound > 0.0f) return vec4(1,0,0, bound+0.1);
    vec4 c = vec4(1,0,0, 1e12), d;
    for (float i=ZERO-1.; i<=1.; i+=1.) {
        for (float j=ZERO-1.; j<=1.; j+=1.) {
            int seed = int(65536.0*hash12(vec2(i,j)+0.1));
            vec3 q = p;
            q.xy -= vec2(0.5,0.4) * (vec2(i,j) + 0.7 * (hash22(vec2(i,j))-0.5)) * smoothstep(1.0,0.5,-0.5*p.z);
            q.z = randt(seed,1.0,0.3)*q.z - randt(seed,0.0,0.5);
            d.w = mapRoot01(vec3(1,1,-1)*q, seed);
            if (col_required) {
                float t = 0.5+randt(seed,0.5,0.2)*sin(randt(seed,4.0,4.0)*p.z+2.0*PI*rand(seed))*sin(4.0*p.y)*sin(4.0*p.x);
                d.xyz = mix(vec3(0.3,0.2,0.05), vec3(0.75,0.55,0.25), smootherstep(t));
            }
            c = smin(c, d, 0.02);
        }
    }
    c.w += 0.02;
    return c;
}

vec4 mapBody(vec3 p, bool col_required) {
    p.x = length(vec2(p.x, 0.01));
    p.z -= 0.1/(pow(p.y-0.5,2.)+1.);
    p.z += 0.2/(pow(p.y-1.2,2.)+1.);
    vec3 q;
    // segments
    q = vec3(1,1,1.1)*(p-vec3(0,0.05,0));;
    vec4 s1 = vec4(mix(vec3(0.7,0.1,0.1), vec3(0.8,0.07,0.1), saturate(0.5+3.0*q.y)),
        sdCapsuleY(q-vec3(0,0,0.04*(cos(8.0*q.y)-1.)), 0.28, 0.15-0.05*q.y));
    q = vec3(1,1,1.1)*(p-vec3(0,0.8,0));
    vec4 s2 = vec4(mix(vec3(0.9,0.3,0.4), vec3(1.0,0.7,0.4), saturate(0.5+5.0*q.z)),
        sdCapsuleY(q-vec3(0,0,1.5*(cos(q.y+0.05)-1.)), 0.22, 0.13-0.05*q.y));
    q = vec3(1,1,1.1)*(p-vec3(0,1.55,0));
    vec4 s3 = vec4(mix(vec3(0.9,0.3,0.4), vec3(1.0,0.7,0.4), saturate(0.5+5.0*q.z)),
        sdCapsuleY(q-vec3(0,0,0.03*(-sin(8.0*q.y)-2.)), 0.3, 0.11-0.03*q.y));
    q = vec3(1,1,1.2)*(p-vec3(0,2.2,0));
    vec4 s4 = vec4(mix(vec3(0.9,0.3,0.4), vec3(1.0,0.7,0.4), saturate(0.5+5.0*q.z)),
        sdCapsuleY(q-vec3(0,0,0.02*(sin(8.0*q.y)-1.)), 0.12, 0.11+0.02*q.y));
    vec4 s = cmin(cmin(s1, s2), cmin(s3, s4));
    // rings
    q = rotx(0.1)*(vec3(1,1,1.1)*(p-vec3(0,-0.39,-0.08)));
    float r0 = sdTorusY(q, 0.12, 0.06);
    q = rotx(0.2)*(vec3(1,1,1.1)*(p-vec3(0,0.41,-0.05)));
    float r1 = sdTorusY(q, 0.08, 0.04);
    q = rotx(0.3)*(vec3(1,0.8,1.1)*(p-vec3(0,1.13,-0.04)));
    float r2 = sdTorusY(q, 0.09, 0.04);
    q = rotx(0.4)*(vec3(1,1,1.1)*(p-vec3(0,1.95,-0.04)));
    float r3 = sdTorusY(q, 0.08, 0.04);
    q = rotx(0.1)*(vec3(1.3,1,1.1)*(p-vec3(0.03,2.6,0)));
    float r4 = sdEllipsoid(q, vec2(0.08-0.2*q.y+0.004*sin(40.0*q.y), 0.2).xyx);
    vec4 r = vec4(vec3(0.05,0.00,0.02), min(min(r0, r1), min(min(r2, r3), r4)));
    // head
    q = vec3(1,1,1.1)*(p-vec3(0,-0.8,-0.1));
    vec4 h0 = vec4(vec3(0.7,0.2,0.05),
        pow(dot(pow(max(abs(q)-vec3(0.0,0.3,0.0),vec3(0.)),vec3(4.)),vec3(1.)),0.25)-max(0.2-0.5*q.y*q.y+0.1*q.y,0.));
    q = p - vec3(0,-0.95,-0.1);
    vec4 h1 = vec4(vec3(0.3,0.0,0.05), length(q)-0.23);
    q = p - vec3(0.1,-1.2,-0.1);
    vec4 eye = vec4(vec3(0.5,0.15,0.07), sdEllipsoid(q, vec3(0.12,0.1,0.15)));
    vec4 h = smin(smin(h0, h1, 0.1), eye, 0.05);
    // put them together
    vec4 d = smin(smin(s, r, 0.03), h, 0.04);
    return d;
}
vec4 mapWing1(vec3 p, bool col_required) {
    p.x = abs(p.x);
    vec3 q;
    q = rotx(0.4)*rotz(0.3)*(p-vec3(1.35,-1.1,-0.05));
    q.z += min(0.4*q.y*q.y, 1.0);
    vec3 r = vec3(
        1.3,
        max(0.16+0.3/(pow(q.x-0.3,2.)+1.5)+0.2*q.y+0.05*q.x,0.1),
        0.4*max(0.05+0.03*q.y,0.01));
    float d = sdEllipsoid(q, r);
    if (!col_required) return vec4(vec3(1.0), d);
    q = q / r + vec3(1, 0, 0);
    float u = 2.0*atan(q.y,q.x);
    float v = dot(q.xy,q.xy)/(2.0*q.x);
    float tu = 0.1*sin(47.0*u)+0.1*sin(31.0*u);
    float tv = -0.1*sin(16.0*v)+0.1*sin(137.0*v)+0.05*sin(73.0*v);
    float t = (0.5+tv+tu) * 1.0/(1.0+exp(-8.0*(u-v+1.4)));
    vec3 col = mix(vec3(0.06,0.03,0.01), vec3(0.6,0.5,0.4), t);
    return vec4(col, d);
}
vec4 mapWing2(vec3 p, bool col_required) {
    p.x = abs(p.x);
    vec3 q;
    q = rotx(0.2)*rotz(-0.35)*(p-vec3(1.13,-0.28,-0.05));
    q.z += min(0.5*q.y*q.y, 1.0);
    vec3 r = vec3(
        1.05,
        max(0.1+0.3/(pow(q.x+0.2,2.)+1.5)+0.2*q.y-0.05*q.x+0.02*exp(2.0*q.x),0.1),
        0.4*max(0.05+0.03*q.y,0.01));
    float d = sdEllipsoid(q, r);
    if (!col_required) return vec4(vec3(1.0), d);
    q = q / r + vec3(1, 0, 0);
    float u = 2.0*atan(q.y,q.x);
    float v = dot(q.xy,q.xy)/(2.0*q.x);
    float tu = 0.1*sin(47.0*u)+0.1*sin(31.0*u);
    float tv = -0.1*sin(16.0*v)+0.1*sin(137.0*v)+0.05*sin(73.0*v);
    float t = (0.5+tv+tu) * 1.0/(1.0+exp(-4.0*(u+0.9)));
    vec3 col = mix(vec3(0.06,0.03,0.01), vec3(0.6,0.5,0.4), t);
    return vec4(col, d);
}
vec4 mapDragonfly(vec3 p, bool col_required) {
    float bound = sdBox(p-vec3(0,0.3,-0.1), vec3(2.8,2.4,0.8));
    if (bound > 0.0f) return vec4(1,0,0, bound+0.1);
    p.z -= 0.2*length(vec2(p.x, 0.5))-0.2;
    vec4 body = mapBody(1.2*p, col_required)/1.2;
    vec4 wing1 = mapWing1(p, col_required);
    vec4 wing2 = mapWing2(p, col_required);
    vec4 d = smin(body, cmin(wing1, wing2), 0.02);
    return d;
}


vec4 mapFlowers(vec3 p, bool col_required) {
    vec3 q;
    q = roty(0.4)*rotz(-0.4)*(p-vec3(-0.75,-0.2,0.2));
    vec4 flower1 = mapFlower02(q/vec3(0.35,0.35,0.45), col_required)*vec4(1,1,1,0.35);
    q = roty(0.3)*rotz(-0.4)*(p-vec3(-0.55,-0.85,0.15));
    vec4 fruit1 = mapFruit02(q/vec3(0.2,0.2,0.2), col_required)*vec4(1,1,1,0.2);
    q = rotx(-0.6)*rotz(-0.5)*(p-vec3(0.05,-0.6,0.2));
    vec4 flower2 = mapFlower01(q/vec3(0.3,0.3,0.4), col_required)*vec4(1,1,1,0.3);
    q = roty(-0.5)*rotz(0.2)*(p-vec3(0.9,0.6,0.35));
    vec4 flower3 = mapFlower02(q/vec3(0.4,0.4,0.4), col_required)*vec4(1,1,1,0.4);
    q = roty(-0.4)*rotz(0.3)*(p-vec3(0.8,-0.2,0.1));
    vec4 flower4 = mapFlower03(q/vec3(0.35,0.35,0.35), col_required)*vec4(1,1,1,0.35);
    q = rotx(1.0)*rotz(-0.6)*(p-vec3(-0.55,0.6,0.25));
    vec4 flower5 = mapFlower01(q/vec3(0.3,0.3,0.4), col_required)*vec4(1,1,1,0.3);
    q = rotx(0.4)*rotz(-0.2)*(p-vec3(0.3,0.95,0.3));
    vec4 flower6 = mapFlower03(q/vec3(0.35,0.35,0.4), col_required)*vec4(1,1,1,0.35);
    q = roty(-0.3)*rotz(0.0)*(p-vec3(0.27,0.1,0.4));
    vec4 fruit2 = mapFruit02(q/vec3(0.2,0.2,0.2), col_required)*vec4(1,1,1,0.2);
    q = roty(-0.3)*rotz(0.5)*(p-vec3(-0.12,-0.15,0.4));
    vec4 fruit3 = mapFruit01(q/vec3(0.18,0.18,-0.18), col_required)*vec4(1,1,1,0.18);
    q = roty(0.3)*rotz(-0.5)*(p-vec3(-0.3,0.2,0.5));
    vec4 fruit4 = mapFruit01(q/vec3(0.18,0.18,-0.18), col_required)*vec4(1,1,1,0.18);
    q = rotx(0.3)*rotz(0.1)*(p-vec3(-0.14,0.53,0.45));
    vec4 fruit5 = mapFruit01(q/vec3(0.17,0.16,-0.18), col_required)*vec4(1,1,1,0.16);
    q = roty(-0.3)*(p-vec3(0.35,0.4,0.8));
    vec4 leaf1 = mapSwirl01(q/vec3(0.5,0.4,0.5), col_required)*vec4(1,1,1,0.4);
    q = roty(-0.1)*(p-vec3(0.35,0.4,0.8));
    vec4 leaf2 = mapSwirl02(q/vec3(0.5,0.4,0.5), col_required)*vec4(1,1,1,0.4);
    q = rotz(1.5)*(p-vec3(-0.08,0.2,0.7));
    vec4 leaf3 = mapSwirl04(q/vec3(0.4,0.4,0.4), col_required)*vec4(1,1,1,0.4);
    const float k = 0.01;
    vec4 c = smin(
        smin(
            smin(smin(smin(flower1, fruit1, k), flower5, k), smin(flower2, flower3, k), k),
            smin(smin(flower4, flower6, k), fruit2, k),
        k),
        smin(
            smin(smin(fruit3, fruit4, k), fruit5, k),
            smin(smin(leaf1, leaf2, k), leaf3, k),
        k),
    k);
    return c;
}
vec4 map(vec3 p, bool col_required) {
    vec4 axes = cmin(cmin(
        vec4(1,0,0, length(p-vec3(2,0,0))-0.1), vec4(0,0.5,0, length(p-vec3(0,1.5,0))-0.1)),
        vec4(0,0,1, length(p-vec3(0,0,3))-0.1));
    vec3 q;
    q = roty(-0.05)*rotx(-0.15)*(p-vec3(0,0,0.0));
    vec4 flowers = mapFlowers(q/vec3(0.8,0.8,0.95), col_required)*vec4(1,1,1,0.8);
    q = p.yxz-vec3(0.1,0.3,-0.7);
    vec4 roots = mapRoots01(q/vec3(0.55,0.55,0.65), col_required)*vec4(1,1,1,0.55);
    q = roty(-0.1)*rotx(-0.5)*(p-vec3(0,-0.1,2.0));
    vec4 dragonfly = mapDragonfly(q/vec3(0.35,-0.4,0.45), col_required)*vec4(1,1,1,0.35);
    const float k = 0.01;
    vec4 c = cmin(
        dragonfly,
        smin(flowers, roots, k)
    );
    return c;
}


float sdf(vec3 p) {
    const float sc = 1.0;
    return map(p/sc, false).w*sc;
}

vec3 sdfGrad(in vec3 p, in float e) {
    // https://iquilezles.org/www/articles/normalsSDF/normalsSDF.htm
    vec3 n = vec3(0.0);
    for(int i=int(ZERO); i<4; i++) {
        vec3 r = -1.0 + 2.0 * vec3(((i+3)>>1)&1, (i>>1)&1, i&1);
        n += r*sdf(p+r*e);
    }
    return n * (.25/e);
}



// color surface or not
#define COLOR 1


// raymarching parameters
#define BOX_RADIUS vec3(2.2, 1.7, 3.5)
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
