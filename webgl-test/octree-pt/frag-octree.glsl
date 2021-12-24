#version 300 es
precision highp float;

out vec4 fragColor;

uniform float iRx;
uniform float iRz;
uniform float iSc;
uniform vec2 iResolution;
uniform int iFrame;
uniform vec4 iMouse;
uniform sampler2D iChannel0;

#define TREE_TEXTURE_SIZE 1024
uniform mediump usampler2D treeSamplerGlass;
uniform mediump usampler2D treeSamplerContent;

uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }

#define PI 3.1415926
#define ZERO min(float(iFrame),0.)
#define EPSILON 1e-4


// Intersection functions

int getUint8(usampler2D treeSampler, int i) {
    ivec2 pos = ivec2((i/4)%TREE_TEXTURE_SIZE, (i/4)/TREE_TEXTURE_SIZE);
    uvec4 sp = texelFetch(treeSampler, pos, 0);
    return int(i%4==0 ? sp.x : i%4==1 ? sp.y : i%4==2 ? sp.z : sp.w);
}
int getUint32(usampler2D treeSampler, int i) {
    ivec2 pos = ivec2((i/4)%TREE_TEXTURE_SIZE, (i/4)/TREE_TEXTURE_SIZE);
    ivec4 sp = ivec4(texelFetch(treeSampler, pos, 0));
	return sp.x + 256 * (sp.y + 256 * (sp.z + 256 * sp.w));
}
ivec3 getUvec3(usampler2D treeSampler, int i) {
	int x = getUint8(treeSampler, i);
	int y = getUint8(treeSampler, i + 1);
	int z = getUint8(treeSampler, i + 2);
	return ivec3(x, y, z);
}

#include "intersector2.glsl"

bool intersectGlass(vec3 ro, vec3 rd, inout float t, in float t1, out vec3 n) {
    vec3 col;
    return intersectOctree(treeSamplerGlass, vec3(-2,-2,-1), vec3(2,2,7), ivec3(1,1,2), 6, 255,
        ro, rd, t, t1, n, col);
}

bool intersectContent(vec3 ro, vec3 rd, inout float t, in float t1, out vec3 n, out vec3 col) {
    return intersectOctree(treeSamplerContent, vec3(-1.5,-1.5,-3.0), vec3(1.5,1.5,3.0), ivec3(1,1,2), 6, 255,
        ro, rd, t, t1, n, col);
}


// Ray sampling and scattering

const int MAT_BACKGROUND = -1;
const int MAT_NONE = 0;
const int MAT_PLANE = 1;
const int MAT_GLASS = 2;
const int MAT_CONTENT = 3;

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


// Rendering

mat2 rot2(float a) { return mat2(cos(a), sin(a), -sin(a), cos(a)); }

vec3 light(vec3 rd) {
    const vec3 sunpos = normalize(vec3(0.2, -0.5, 0.2));
    vec3 col = vec3(0.8);
    col *= 1.0 + 2.0 * pow(max(dot(rd, sunpos), 0.), 4.);
    col *= mix(1.0, 4.0, pow(1.0-rd.z, 2.0));
    vec3 sun = (dot(rd,sunpos)>0.9 ? 1.0 : 0.0) * vec3(8.0);
    return 0.5*col + sun;
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
            material = MAT_PLANE;
        }

        // glass
        t = 0.0;
        if (intersectGlass(ro, rd, t, min_t, min_n)) {
            min_t = t;
            min_ro = ro + rd * t, min_rd = rd;
            min_n = normalize(min_n);
            col = vec3(1.0);
            material = MAT_GLASS;
        }

        // content
        t = 0.0;
        if (inside_glass) {
            if (intersectContent(ro-vec3(0.0,0.0,3.0), rd, t, min_t, min_n, col)) {
                min_t = t;
                min_ro = ro + rd * t, min_rd = rd;
                min_n = normalize(min_n);
                material = MAT_CONTENT;
            }
        }

        // update ray
        if (material == MAT_BACKGROUND) {
            m_col *= light(rd);
            break;
        }
        if (inside_object);
        else if (inside_glass) m_col *= exp(-0.1*vec3(0.0,0.2,0.4)*min_t);
        min_n = dot(rd, min_n) < 0. ? min_n : -min_n;  // ray hits into the surface
        ro = min_ro, rd = min_rd;
        if (material == MAT_PLANE) {
            // faked light behind the glass
            vec2 xy = min_ro.xy;
            float c = length(rot2(-1.0)*(xy-vec2(0.0,15.0))/vec2(2.0,1.0))-10.0;
            col = vec3(0.5)-0.3*tanh(0.4*c);
            rd = sampleCookTorrance(-rd, min_n, 0.01, 0.1, 0.01, col, col, m_col);
        }
        else if (material == MAT_GLASS) {
            vec2 eta = inside_glass ? vec2(1.5, 1.0) : vec2(1.0, 1.5);
            rd = sampleFresnelDielectric(rd, min_n, eta.x, eta.y);
            if (dot(rd, min_n) < 0.0) inside_glass = !inside_glass;
        }
        else if (material == MAT_CONTENT) {
            rd = sampleCookTorrance(-rd, min_n, 0.1, 1.0, 0.1, vec3(1.0), vec3(1.0), m_col);
            m_col *= 1.0*pow(col, vec3(1.0));
            if (dot(rd, min_n) < 0.0) inside_object = !inside_object;
        }
        if (m_col == vec3(0.0)) break;
        if (inside_object) return 1e12f*vec3(1,-1,-1);  // red warning
    }
    return m_col + t_col;
}

vec2 randomUnitDisk() {
    float a = 2.0*PI*rand01();
    return sqrt(rand01())*vec2(cos(a), sin(a));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // random number seed
    seed = uint(fragCoord.x)*uint(fragCoord.y)*uint(iFrame+1);
    seed = randu() + 161u*uint(fragCoord.y);
    seed = randu() + 239u*uint(fragCoord.x);
    seed = randu() + 197u*uint(iFrame+1);

    // camera parameters
    const vec3 CENTER = vec3(0, 0, 3.0);
    const float DIST = 20.0;  // larger = smaller
    const float VIEW_FIELD = 0.4;  // larger = larger + more perspective
    const float APERTURE = 0.2;  // larger = blurred

    // camera
    vec3 w = vec3(cos(iRx)*vec2(cos(iRz),sin(iRz)), sin(iRx));
    vec3 u = vec3(-sin(iRz),cos(iRz),0);
    vec3 v = cross(w,u);
    vec3 ro = DIST*w + CENTER;
    vec2 uv = iSc*(2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy-1.0);
    vec2 sc = iResolution.xy/length(iResolution.xy);
    vec2 offset = APERTURE*randomUnitDisk();
    vec3 rd = mat3(u,v,-w)*vec3(VIEW_FIELD*uv*sc+offset/DIST, 1.0);
    rd = normalize(rd);
    ro -= offset.x*u+offset.y*v;

    // calculate pixel color
    vec3 col = mainRender(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(int(fragCoord.x), int(fragCoord.y)), 0);
    if (iMouse.z>0. || iFrame==0) rgbn = vec4(0.0);
    fragColor = rgbn + vec4(col, 1.0);
}


void main(void) {
    mainImage(fragColor, gl_FragCoord.xy);
}
