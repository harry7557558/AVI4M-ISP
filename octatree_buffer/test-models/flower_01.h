#define PI 3.1415926f

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


vec4 mapFlower(vec3 p, bool col_required) {
    float x, y; vec3 q;
    float r = length(p.xy()), a = atan(p.y, p.x);
    float a_ = asin(0.99f*sin(2.5f*a))/2.5f;
    vec3 p_ = vec3(r*vec2(cos(a_), sin(a_)), p.z);
    // petals
    q = p_;
    x = q.x-1.0f, y = q.y;
    q.x -= 1.2f;
    q.z += 0.1f*exp(-sqr(3.0f*y)) - 0.1f*length(vec2(r,0.2f)) - 0.4f/(1.0f+exp(-4.0f*(r-1.6f)));
    float w = 1.05f*(0.2f
        + 0.3f*exp(-sqr(1.2f*(x-0.2f)))
        + 0.1f*exp(-sqr(1.5f*(x+0.8f)))
        - 0.05f*exp(-sqr(4.0f*(x-1.0f)))
        + (0.1f+0.1f*cos(4.0f*a)*sin(5.0f*a))*exp(-sqr(1.0f*(x-1.0f))));
    float h = 0.2f*w+0.1f*exp(-sqr(3.0f*(x+0.8f)));
    float petal_d = sdEllipsoid(q, vec3(1.0f, w, h));
    vec3 petal_c = vec3(0.0f);
    if (col_required) {
        float t = 0.5f+0.5f*cos(10.0f*a_)-exp(-4.0f*(r-0.5f));
        petal_c = mix(vec3(0.8f,0.55f,0.65f), vec3(0.85f,0.7f,0.8f), smoothstep(0.f,1.f,t));
    }
    vec4 petal = vec4(petal_c, petal_d);
    // style/sepal
    q = vec3(r, a_, p.z);
    float leaf_d = smin(  // green part
        sdSegment(q.xz(), vec2(0.24f,-0.05f), vec2(0.32f,0.0f))-0.05f,  // ring
        sdSegment(q.xz()+0.05f*sin(8.0f*q.z), vec2(0.12f,-0.4f), vec2(0.22f,-0.05f))-0.05f,  // wall
        0.05f);
    leaf_d = smin(leaf_d, sdEllipsoid(p_-vec3(0.55f,0,-0.1f), vec3(0.4f,0.1f,0.05f)), 0.05f);  // sepal
    leaf_d = smin(leaf_d,
        sdSegment(p-vec3(0.01f*sin(10.0f*r),0,0), vec3(0,0,-0.6f),vec3(0,0,1.0f))
            -max(0.06f/(1.0f+exp(10.0f*(q.z-0.5f)))+0.02f,0.f),  // style
        0.1f);
    vec4 leaf = vec4(0.5f,0.7f,0, leaf_d);
    // filament/anther
    a_ = asin(0.99f*sin(7.5f*atan(p_.y,p_.x)))/7.5f;
    p_ = vec3(r*vec2(cos(a_), sin(a_)), p.z);
    q = p_ - vec3(0.1f+(0.7f+0.2f*sin(7.0f*a))/(1.0f+exp(-8.0f*(p.z-0.2f))),0,-0.2f);
    h = 0.8f+0.2f*sin(4.0f*a);
    vec4 filament = vec4(0.85f,0.85f,0.f,
        0.5f*sdCapsule(q, h, 0.01f+0.01f*exp(-sqr(4.0f*(p.z-0.2f)))));
    vec4 anther = vec4(0.75f,0.55f,0.f,
        sdEllipsoid(q-vec3(0.02f,0,h), vec3(0.08f,0.04f,0.05f)));
    // put all together
    vec4 d = cmin(smin(petal, leaf, 0.05f-0.03f*exp(-r)), smin(filament, anther, 0.01f));
    //d = max(d, p.y);
    return d;
}

vec4 map(vec3 p, bool col_required) {
	p *= 1.2f;
    p += 0.08f*sin(4.0f*p.yzx())*cos(4.0f*p.xyz())*sin(4.0f*p.zxy());
    vec4 d = mapFlower(p, col_required);
    d.w *= 0.8f;
    return d;
}

