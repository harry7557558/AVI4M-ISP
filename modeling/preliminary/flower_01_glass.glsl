#iChannel0 "self"

#iChannel1 "../cubemaps/shadertoy_uffizi_gallery/{}.jpg"
#iChannel1::Type "CubeMap"



uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }

#define PI 3.1415926
#define ZERO min(iTime,0.)
#define EPSILON 1e-4


float sqr(float x) { return x*x; }


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
float sdEllipsoid(vec3 p, vec3 r) {
    float k1 = length(p/r);
    float k2 = length(p/(r*r));
    return k1*(k1-1.0)/k2;
}
float sdCapsule(vec3 p, float h, float r) {
    p.z = p.z-clamp(p.z, 0.0, h);
    return length(p) - r;
}
float sdTorus(vec3 p, float R, float r) {
    vec2 q = vec2(length(p.xz)-R,p.y);
    return length(q)-r;
}
float sdLnNormEllipsoid(vec3 p, vec3 r, vec3 n) {
    // not so good, results in too small gradient
    float d = pow(dot(pow(abs(p)/r,n),vec3(1)),1.0/max(max(n.x,n.y),n.z))-1.0;
    float m = min(min(r.x,r.y),r.z);
    return d * m;
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


// Modeling

vec4 mapFlower0(vec3 p, bool col_required) {
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
    vec4 leaf = vec4(0.6,0.65,0, leaf_d);
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
    vec4 d = cmin(smin(petal, leaf, 0.05-0.03*exp(-r)), smin(filament, anther, 0.01));
    //d = max(d, p.y);
    return d;
}

vec4 mapFlower(vec3 p, bool col_required) {
    p -= vec3(0, 0, 1.2);
    vec4 d = mapFlower0(p, col_required);
    d.w *= 0.8;
    return d;
}


float mapContent(vec3 p) {
    return mapFlower(p, false).w;
}

vec3 gradContent(in vec3 p) {
    const float e = 0.001;
	float a = mapContent(p+vec3(e,e,e));
	float b = mapContent(p+vec3(e,-e,-e));
	float c = mapContent(p+vec3(-e,e,-e));
	float d = mapContent(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}

float mapGlass(vec3 p) {
    const vec3 r = vec3(2.5,2.5,1.3);
    p.z -= r.z+0.01;
    return sdLnNormEllipsoid(p, vec3(r), vec3(6.0,6.0,12.0));
}

vec3 gradGlass(in vec3 p) {
    const float e = 0.001;
	float a = mapGlass(p+vec3(e,e,e));
	float b = mapGlass(p+vec3(e,-e,-e));
	float c = mapGlass(p+vec3(-e,e,-e));
	float d = mapGlass(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}

// Ray sampling and scattering

vec3 sampleCosWeighted(vec3 n) {
    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    float rn = rand01();
    float an = 2.0*PI*rand01();
    vec2 rh = sqrt(rn) * vec2(cos(an), sin(an));
    float rz = sqrt(1. - rn);
    return rh.x * u + rh.y * v + rz * n;
}

vec3 sampleFresnelDielectric(vec3 rd, vec3 n, float n1, float n2) {
    float eta = n1 / n2;
    float ci = -dot(n, rd);
    if (ci < 0.0) ci = -ci, n = -n;
    float ct = 1.0 - eta * eta * (1.0 - ci * ci);
    if (ct < 0.0) return rd + 2.0*ci*n;
    ct = sqrt(ct);
    float Rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct);
    float Rp = (n1 * ct - n2 * ci) / (n1 * ct + n2 * ci);
    float R = 0.5 * (Rs * Rs + Rp * Rp);
    return rand01() > R ?
        rd * eta + n * (eta * ci - ct)  // refraction
        : rd + 2.0*ci*n;  // reflection
}

vec3 sampleCookTorrance(
    vec3 wi, vec3 n,
    float alpha,  // roughness
    float f0,  // ratio of reflection along the normal
    float lambertian,  // ratio of lambertian coefficient
    vec3 lambert_col,  // lambertian color
    vec3 microfacet_col,  // microfacet color
    inout vec3 m_col
    ) {

    // importance sample Lambertian
    if (rand01() < lambertian) {
        vec3 wo = sampleCosWeighted(n);
        m_col *= lambert_col;
        return wo;
    }

    // transform
    vec3 u = normalize(cross(n, vec3(1.2345, 2.3456, -3.4561)));
    vec3 v = cross(u, n);
    wi = vec3(dot(wi, u), dot(wi, v), dot(wi, n));
    vec3 wo, m;  // out and half vector

    // GGX importance sampling
    float su = 2.0*PI*rand01();
    float sv = rand01();
    sv = atan(alpha*sqrt(sv/(1.0-sv)));
    m = vec3(sin(sv)*vec2(cos(su),sin(su)), cos(sv));
    wo = -(wi-2.0*dot(wi,m)*m);
    float D = wo.z<0. ? 0. : 4.0*dot(wi, m);

    // Geometry
    float tan2_theta_i = (1.0-wi.z*wi.z)/(wi.z*wi.z);
    float tan2_theta_o = (1.0-wo.z*wo.z)/(wo.z*wo.z);
    float lambda_i = 0.5*(sqrt(1.0+alpha*alpha*tan2_theta_i)-1.0);
    float lambda_o = 0.5*(sqrt(1.0+alpha*alpha*tan2_theta_o)-1.0);
    float G = 1.0/(1.0+lambda_i+lambda_o);

    // Fresnel
    float F = f0 + (1.0-f0)*pow(1.0-dot(wi, m), 5.0);

    // Put all together
    float Fr = D*G*F / (4.0*wi.z*wo.z+EPSILON);
    float Fr_cos = Fr * wo.z;  // wo is the direction of light in path tracing
    m_col *= Fr_cos * microfacet_col;
    return wo.x * u + wo.y * v + wo.z * n;
}


// Raymarching functions

const int MAT_BACKGROUND = -1;
const int MAT_NONE = 0;
const int MAT_PLANE = 1;
const int MAT_GLASS = 2;
const int MAT_CONTENT = 3;

bool intersectGlass(vec3 ro, vec3 rd, inout float t, in float t1) {
    const float STEP = 100.0, MIN_STEP = 0.01;
    float v_old = mapGlass(ro+rd*t), v;
    float dt = min(STEP, abs(v_old));
    for (int i=int(ZERO); i<128; i++) {
        t += dt;
        if (t > t1) return false;
        v = mapGlass(ro+rd*t);
        if (v*v_old<0.) break;
        dt = clamp(abs(v_old=v), MIN_STEP, STEP);
    }
    if (v*v_old<0.) {
        for (int s = int(ZERO); s < 8; s += 1) {
            v_old = v, dt *= -0.5;
            for (int i = int(ZERO); i < 2; i++) {
                t += dt, v = mapGlass(ro+rd*t);
                if (v*v_old<0.0) break;
            }
        }
        return true;
    }
    return false;
}

bool intersectContent(vec3 ro, vec3 rd, inout float t, in float t1) {
    const float STEP = 1.0, MIN_STEP = 0.005;
    float v_old = mapContent(ro+rd*t), v;
    float dt = min(STEP, abs(v_old));
    for (int i=int(ZERO); i<128; i++) {
        t += dt;
        if (t > t1) return false;
        v = mapContent(ro+rd*t);
        if (v*v_old<0.) break;
        dt = clamp(abs(v_old=v), MIN_STEP, STEP);
    }
    if (v*v_old<0.) {
        // binary search
        float t0 = t-dt, t1 = t;
        float v0 = v_old, v1 = v;
        for (int s = int(ZERO); s < 8; s++) {
            float t = 0.5*(t0+t1);
            float v = mapContent(ro+rd*t);
            if (v*v0 < 0.0) t1 = t, v1 = v;
            else t0 = t, v0 = v;
        }
        t = t0;
        return true;
    }
    return false;
}


// Rendering

vec3 light(vec3 rd) {
    const vec3 sunpos = normalize(vec3(-0.2, -0.5, 0.5));
    vec3 col = texture(iChannel1, rd.xyz).xyz;
    vec3 amb = vec3(1.0) + vec3(2.0) * pow(max(dot(rd, sunpos), 0.), 4.);
    vec3 sun = (dot(rd,sunpos)>0.9 ? 1.0 : 0.0) * vec3(10.0);
    return col * 0.1*amb + sun;
}

vec3 mainRender(vec3 ro, vec3 rd) {

    vec3 m_col = vec3(1.0), t_col = vec3(0.0), col;
    bool inside_glass = false, inside_object = false;

    for (int iter = 0; iter < 64; iter++) {
        ro += EPSILON*rd;
        float t, min_t = 1e12;
        vec3 n, min_n;
        vec3 min_ro = ro, min_rd = rd;
        vec3 min_emi = vec3(0.0);
        int material = MAT_BACKGROUND;

        // plane
        t = -ro.z / rd.z;
        if (t > 0.0) {
            min_t = t, min_n = vec3(0, 0, 1);
            min_ro = ro + rd * t, min_rd = rd;
            col = vec3(0.9, 0.95, 0.98);
            material = MAT_PLANE;
        }

        // glass
        t = 0.0;
        if (intersectGlass(ro, rd, t, min_t)) {
            min_t = t;
            min_ro = ro + rd * t, min_rd = rd;
            min_n = normalize(gradGlass(min_ro));
            col = vec3(1.0);
            material = MAT_GLASS;
        }

        // content
        t = 0.0;
        if (inside_glass) {
            if (intersectContent(ro, rd, t, min_t)) {
                min_t = t;
                min_ro = ro + rd * t, min_rd = rd;
                min_n = normalize(gradContent(min_ro));
                col = mapFlower(ro+rd*t, true).xyz;
                material = MAT_CONTENT;
            }
        }

        // update ray
        if (material == MAT_BACKGROUND) {
            col = light(rd);
            return m_col * col + t_col;
        }
        if (inside_object);
        else if (inside_glass) m_col *= exp(-0.1*vec3(0.0,0.2,0.4)*min_t);
        min_n = dot(rd, min_n) < 0. ? min_n : -min_n;  // ray hits into the surface
        ro = min_ro, rd = min_rd;
        if (material == MAT_PLANE) {
            rd = sampleCookTorrance(-rd, min_n, 0.05, 0.3, 0.2, col, col, m_col);
        }
        else if (material == MAT_GLASS) {
            vec2 eta = inside_glass ? vec2(1.5, 1.0) : vec2(1.0, 1.5);
            rd = sampleFresnelDielectric(rd, min_n, eta.x, eta.y);
            if (dot(rd, min_n) < 0.0) inside_glass = !inside_glass;
        }
        else if (material == MAT_CONTENT) {
            //t_col += 0.5 * m_col * col;
            rd = sampleCookTorrance(-rd, min_n, 0.5, 0.8, 0.4, vec3(1.0), vec3(1.0), m_col);
            m_col *= 1.2*col*(col+0.5);
            if (dot(rd, min_n) < 0.0) inside_object = !inside_object;
        }
        if (m_col == vec3(0.0)) return t_col;
        //if (mapContent(ro) < 0.0) return vec3(100.0, -100.0, -100.0);
        if (inside_object) return vec3(100.0, -100.0, -100.0);
    }
    return m_col + t_col;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // random number seed
    seed = uint(fragCoord.x)*uint(fragCoord.y)*uint(iFrame+1);
    seed = randu() + 161u*uint(fragCoord.y);
    seed = randu() + 239u*uint(fragCoord.x);
    seed = randu() + 197u*uint(iFrame+1);

    // camera
    float rx = 1.8*(iMouse.y/iResolution.y)-0.3;
    //float rx = 3.14*(iMouse.y/iResolution.y)-1.57;
    float rz = -iMouse.x/iResolution.x*4.0*3.14;
    vec3 w = vec3(cos(rx)*vec2(cos(rz),sin(rz)), sin(rx));
    vec3 u = vec3(-sin(rz),cos(rz),0);
    vec3 v = cross(w,u);
    vec3 ro = 20.0*w + vec3(0, 0, 1.2);
    vec2 uv = 2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 3.0*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    vec3 col = mainRender(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(fragCoord), 0);
    if (iMouse.z>0.) rgbn.w = 0.0;
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
