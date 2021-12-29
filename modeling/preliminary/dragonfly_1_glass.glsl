#iChannel0 "self"

#iChannel1 "../cubemaps/shadertoy_uffizi_gallery/{}.jpg"
#iChannel1::Type "CubeMap"



uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }

#define PI 3.1415926
#define ZERO min(iTime,0.)
#define EPSILON 1e-4


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
float sdEllipsoid(vec3 p, vec3 r) {
    float k1 = length(p/r);
    float k2 = length(p/(r*r));
    return k1*(k1-1.0)/k2;
}
float sdCapsule(vec3 p, float h, float r) {
    p.y = abs(p.y)-min(abs(p.y), h);
    return length(p) - r;
}
float sdTorus(vec3 p, float R, float r) {
    vec2 q = vec2(length(p.xz)-R,p.y);
    return length(q)-r;
}
float sdLnNormEllipsoid(vec3 p, vec3 r, float n) {
    // not so good, results in too small gradient
    float d = pow(dot(pow(abs(p)/r,vec3(n)),vec3(1)),1./n)-1.0;
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

vec4 mapBody(vec3 p, bool col_required) {
    p.x = length(vec2(p.x, 0.01));
    p.z -= 0.1/(pow(p.y-0.5,2.)+1.);
    p.z += 0.2/(pow(p.y-1.2,2.)+1.);
    vec3 q;
    // segments
    q = vec3(1,1,1.1)*(p-vec3(0,0.05,0));;
    vec4 s1 = vec4(mix(vec3(0.7,0.1,0.1), vec3(0.8,0.07,0.1), saturate(0.5+3.0*q.y)),
        sdCapsule(q-vec3(0,0,0.04*(cos(8.0*q.y)-1.)), 0.28, 0.15-0.05*q.y));
    q = vec3(1,1,1.1)*(p-vec3(0,0.8,0));
    vec4 s2 = vec4(mix(vec3(0.9,0.3,0.4), vec3(1.0,0.7,0.4), saturate(0.5+5.0*q.z)),
        sdCapsule(q-vec3(0,0,1.5*(cos(q.y+0.05)-1.)), 0.22, 0.13-0.05*q.y));
    q = vec3(1,1,1.1)*(p-vec3(0,1.55,0));
    vec4 s3 = vec4(mix(vec3(0.9,0.3,0.4), vec3(1.0,0.7,0.4), saturate(0.5+5.0*q.z)),
        sdCapsule(q-vec3(0,0,0.03*(-sin(8.0*q.y)-2.)), 0.3, 0.11-0.03*q.y));
    q = vec3(1,1,1.2)*(p-vec3(0,2.2,0));
    vec4 s4 = vec4(mix(vec3(0.9,0.3,0.4), vec3(1.0,0.7,0.4), saturate(0.5+5.0*q.z)),
        sdCapsule(q-vec3(0,0,0.02*(sin(8.0*q.y)-1.)), 0.12, 0.11+0.02*q.y));
    vec4 s = cmin(cmin(s1, s2), cmin(s3, s4));
    // rings
    q = rotx(0.1)*(vec3(1,1,1.1)*(p-vec3(0,-0.39,-0.08)));
    float r0 = sdTorus(q, 0.12, 0.06);
    //q = rotx(0.1)*rotz(-0.1)*(p-vec3(0.14,-0.65,0.02));
    //r0 = smin(r0, sdCapsule(q, 0.2, 0.08), 0.01);
    q = rotx(0.2)*(vec3(1,1,1.1)*(p-vec3(0,0.41,-0.05)));
    float r1 = sdTorus(q, 0.08, 0.04);
    q = rotx(0.3)*(vec3(1,0.8,1.1)*(p-vec3(0,1.13,-0.04)));
    float r2 = sdTorus(q, 0.09, 0.04);
    q = rotx(0.4)*(vec3(1,1,1.1)*(p-vec3(0,1.95,-0.04)));
    float r3 = sdTorus(q, 0.08, 0.04);
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
    p -= vec3(0, 0, 1);
    p = rotx(0.2)*p;
    p.z -= 0.2*length(vec2(p.x, 0.5))-0.2;
    vec4 body = mapBody(1.2*p, col_required)/1.2;
    vec4 wing1 = mapWing1(p, col_required);
    vec4 wing2 = mapWing2(p, col_required);
    //return cmin(wing1, wing2);
    vec4 d = smin(body, cmin(wing1, wing2), 0.02);
    //d.xyz = 0.05+0.95*pow(d.xyz, vec3(0.8));
    d.xyz = saturate(d.xyz);
    return d;
}

float mapContent(vec3 p) {
    return mapDragonfly(p, false).w;
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
    const vec3 r = vec3(3.0,3.0,1.0);
    p.z -= r.z+0.01;
    return sdLnNormEllipsoid(p, vec3(r), 12.0);
    return max(max(abs(p.x)-r.x,abs(p.y)-r.y),abs(p.z)-r.z);
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
    const float STEP = 0.5, MIN_STEP = 0.01;
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
    const float STEP = 0.1, MIN_STEP = 0.005;
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
        for (int s = int(ZERO); s < 8; s += 1) {
            v_old = v, dt *= -0.5;
            for (int i = int(ZERO); i < 2; i++) {
                t += dt, v = mapContent(ro+rd*t);
                if (v*v_old<0.0) break;
            }
        }
        return true;
    }
    return false;
}


// Rendering

vec3 light(vec3 rd) {
    vec3 col = texture(iChannel1, rd.xyz).xyz;
    vec3 bri = vec3(1.0) + vec3(2.0) * pow(max(dot(rd, normalize(vec3(-0.2, -0.5, 0.5))), 0.), 4.);
    return col * bri;
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
                col = mapDragonfly(ro+rd*t, true).xyz;
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
            t_col += m_col * col;
            rd = sampleCookTorrance(-rd, min_n, 0.05, 0.3, 0.2, col, col, m_col);
            if (dot(rd, min_n) < 0.0) inside_object = !inside_object;
        }
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
    vec3 ro = 20.0*w + vec3(0, 0, 0.7);
    vec2 uv = 2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy - vec2(1.0);
    vec3 rd = mat3(u,v,-w)*vec3(uv*iResolution.xy, 3.0*length(iResolution.xy));
    rd = normalize(rd);

    // calculate pixel color
    vec3 col = mainRender(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(fragCoord), 0);
    if (iMouse.z>0.) rgbn.w = 0.0;
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}
