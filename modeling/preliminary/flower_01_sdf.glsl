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



vec4 mapFlower(vec3 p, bool col_required) {
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
    vec4 d = cmin(smin(petal, leaf, 0.05-0.03*exp(-r)), smin(filament, anther, 0.01));
    //d = max(d, p.y);
    return d;
}


vec4 map(vec3 p, bool col_required) {
    p += 0.08*sin(4.0*p.yzx)*cos(4.0*p.xyz)*sin(4.0*p.zxy);
    vec4 d = mapFlower(p, col_required);
    d.w *= 0.8;
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
