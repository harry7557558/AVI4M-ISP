#define PI 3.1415926f

float sqr(float x) { return x*x; }
vec2 cossin(float x) { return vec2(cos(x),sin(x)); }

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
mat3 roty(float a) {
    return mat3(
        cos(a), 0, -sin(a),
        0, 1, 0,
        sin(a), 0, cos(a)
    );
}

// geometric primitives
float sdSegment(vec2 p, vec2 a, vec2 b) {
    vec2 q = p-a, d = b-a;
    float h = dot(q,d)/dot(d,d);
    return length(q-d*clamp(h,0.f,1.f));
}
float sdSegment(vec3 p, vec3 a, vec3 b) {
    vec3 q = p-a, d = b-a;
    float h = dot(q,d)/dot(d,d);
    return length(q-d*clamp(h,0.f,1.f));
}
float sdCapsule(vec3 p, float h, float r) {
    p.z = p.z-clamp(p.z, 0.0f, h);
    return length(p) - r;
}
float sdEllipsoid(vec3 p, vec3 r) {
    float k1 = length(p/r);
    float k2 = length(p/(r*r));
    return k1*(k1-1.0f)/k2;
}
float sdTorus(vec3 p, float R, float r) {
    vec2 q = vec2(length(p.xz())-R,p.y);
    return length(q)-r;
}

// smoothed blending functions
float smin(float a, float b, float k) {
    float h = clamp(0.5f + 0.5f * (b - a) / k, 0.f, 1.f);
    return mix(b, a, h) - k * h * (1.0f - h);
}
float smax(float a, float b, float k) {
    return -smin(-a, -b, k);
}
float smootherstep(float x) {
    x = clamp(x, 0.f, 1.f);
    return x*x*x*(10.f+x*(-15.f+6.f*x));
}

// color blending functions
vec4 cmin(vec4 c1, vec4 c2) {
    return c1.w<c2.w ? c1 : c2;
}
vec4 smin(vec4 a, vec4 b, float k) {
    float h = clamp(0.5f + 0.5f * (b.w - a.w) / k, 0.f, 1.f);
    float d = mix(b.w, a.w, h) - k * h * (1.0f - h);
    return vec4(mix(b.xyz(), a.xyz(), h), d);
}

vec4 mapPetal(vec3 p, float sc, float curve, bool col_required) {
    p /= sc;
    p.z -= curve*min(dot(p.xy(),p.xy()), 1.0f);
    p.z += 0.1f/(1.0f+sqr(2.0f*p.y));
    vec3 pr = vec3(1.0f, 0.5f, 0.1f);
    pr.y = 0.3f + 0.35f*exp(-sqr(1.2f*(p.x+0.4f)));
    float d = sdEllipsoid(p, pr);
    if (dot(p.xy(),p.xy())<1.0f) d /= length(vec3(curve*2.0f*p.xy(),1.0f));
    vec3 col = mix(vec3(0.6f,0.55f,0.05f), vec3(0.85f,0.9f,0.3f), 0.5f+0.5f*tanh(2.0f*p.x));
    return vec4(col, d*sc);
}

vec4 mapFilament(vec3 p, float h, bool col_required) {
    vec4 fil = vec4(0.95f,0.95f,0.9f, sdCapsule(p, h, 0.024f));
    vec4 ant = vec4(0.98f,0.95f,0.75f, sdEllipsoid(p-vec3(0,0,h), vec3(0.08f,0.07f,0.06f)));
    vec4 d = smin(fil, ant, 0.01f);
    return d;
}

vec4 mapFlower(vec3 p, bool col_required) {
    p /= vec3(1.0f, 1.0f, 1.1f);
    //p += 0.04f*sin(4.0f*p.yzx())*cos(4.0f*p.xyz())*sin(4.0f*p.zxy());
    p.y += 0.02f*sin(8.0f*p.z);
    float r = length(p.xy()), a = atan(p.y, p.x);
    float x, y; vec3 q;
    // petals
    q = vec3(r*cossin(asin(0.999f*sin(3.5f*a-0.2f*cos(4.0f*r)*sin(5.0f*a)+0.9f))/3.5f), p.z);
    vec4 p0 = mapPetal(roty(-0.4f)*(q-vec3(0.8f,0,-0.2f)), 0.45f, -0.1f, col_required);
    q = vec3(r*cossin(asin(0.999f*sin(3.5f*a-0.2f*sin(4.0f*r)*cos(5.0f*a)+0.0f))/3.5f), p.z);
    vec4 p1 = mapPetal(roty(-0.1f)*(q-vec3(1.0f,0,-0.05f)), 0.5f, -0.2f, col_required);
    q = vec3(r*cossin(asin(0.999f*sin(3.0f*a-0.2f*cos(4.0f*r)*cos(3.0f*a)-0.9f))/3.0f), p.z);
    vec4 p2 = mapPetal(rotx(0.1f)*roty(0.1f)*(q-vec3(1.0f,0,0.1f)), 0.55f, -0.1f, col_required);
    q = vec3(r*cossin(asin(0.999f*sin(3.0f*a-0.2f*cos(4.0f*r)*cos(4.0f*a)-2.0f))/3.0f), p.z);
    vec4 p3 = mapPetal(roty(0.4f)*(q-vec3(0.7f,0,0.16f)), 0.55f, -0.05f, col_required);
    q = vec3(r*cossin(asin(0.999f*sin(2.5f*a-0.2f*sin(3.0f*r)*sin(3.0f*a)-2.7f))/2.5f), p.z);
    vec4 p4 = mapPetal(roty(0.6f)*(q-vec3(0.6f,0,0.2f)), 0.45f, 0.05f, col_required);
    q = vec3(r*cossin(asin(0.999f*sin(2.0f*a-0.3f*sin(3.0f*r)*sin(3.0f*a)-3.5f))/2.0f), p.z);
    vec4 p5 = mapPetal(roty(0.7f)*(q-vec3(0.4f,0,0.15f)), 0.4f, 0.2f, col_required);
    vec4 petal = smin(smin(smin(p0, p1, 0.01f), smin(p2, p3, 0.01f), 0.01f), smin(p4, p5, 0.01f), 0.01f);
    // filaments
    q = vec3(r*cossin(asin(0.999f*sin(3.5f*a+0.9f))/3.5f), p.z);
    vec4 f1 = mapFilament(roty(-0.7f-0.1f*cos(4.0f*a))*(q-vec3(0.2f,0,0.0f)), 0.8f, col_required);
    q = vec3(r*cossin(asin(0.999f*sin(4.0f*a+1.5f))/4.0f), p.z);
    vec4 f2 = mapFilament(roty(-0.8f+0.1f*sin(3.0f*a))*(q-vec3(0.18f,0,0.0f)), 0.8f, col_required);
    vec4 filament = smin(f1, f2, 0.01f);
    // disk
    q = vec3(r*cossin(asin(0.98f*sin(3.5f*a-0.5f))/3.5f), p.z);
    vec3 br = vec3(0.4f, 0.3f-0.2f*exp(-sqr(2.0f*(q.x-1.0f))), 0.05f);
    vec4 d1 = vec4(0.4f,0.35f,0.0f, sdEllipsoid(roty(0.2f)*(q-vec3(0.2f,0,-0.25f)), br));  // sepals
    q = vec3(r*cossin(asin(0.8f*sin(3.5f*a))/3.5f), p.z);
    vec4 d2 = vec4(0.6f,0.5f,0.0f, sdEllipsoid(q-vec3(0.05f,0,-0.1f), vec3(0.25f,0.2f,0.3f)));  // ovary
    vec4 d3 = vec4(0.9f,0.8f,0.0f, sdCapsule(q-vec3(0.05f,0,0),0.56f,0.08f-0.03f*sin(2.0f*p.z)));  // style
    vec4 d4 = vec4(0.95f,0.95f,0.6f, length(q-vec3(0.01f,0,0.56f))-0.06f);  // stigma
    vec4 disk = smin(smin(d1, d2, 0.1f), smin(d3, d4, 0.05f), 0.1f);
    // put them together
    vec4 d = smin(cmin(petal, filament), disk, 0.01f);
    return d;
}

vec4 map(vec3 p, bool col_required) {
    vec4 d = mapFlower(p, col_required);
    return d;
}
