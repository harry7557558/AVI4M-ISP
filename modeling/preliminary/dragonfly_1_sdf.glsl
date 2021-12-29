// SDF visualizer v2
// Allows a usually larger average marching step without missing thin surface
// Expected to be faster than v1

// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926
//#define ZERO min(iTime,0.)
#define ZERO 0.0


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

vec4 map(vec3 p, bool col_required) {
    p = rotx(0.2)*p;
    p.z -= 0.2*length(vec2(p.x, 0.5))-0.2;
    vec4 body = mapBody(1.2*p, col_required)/1.2;
    vec4 wing1 = mapWing1(p, col_required);
    vec4 wing2 = mapWing2(p, col_required);
    vec4 d = smin(body, cmin(wing1, wing2), 0.02);
    //d.xyz = 0.05+0.95*pow(d.xyz, vec3(0.8));
    d.xyz = saturate(d.xyz);
    return d;
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
#define BOX_RADIUS vec3(3.0, 3.0, 1.0)
#define STEP 0.1
#define MIN_STEP 0.005
#define MAX_STEP 100.

// rendering parameters
#define FIELD_EMISSION 0.1
#define ISOSURFACE_FREQUENCY 8.0
#define DISCONTINUITY_OPACITY 0.01
#define SURFACE_GRADIENT 5.0

// projection parameters
#define PERSPECTIVE 10.0  /* larger: less perspective effect */
#define SCALE 6.0  /* image appears smaller when this is set larger */


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
    col *= 0.2+0.05*grad.y+max(dot(normalize(grad), light),0.0);
#else
    col = 0.2+0.05*grad.y+col*max(dot(normalize(grad), light),0.0);
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
