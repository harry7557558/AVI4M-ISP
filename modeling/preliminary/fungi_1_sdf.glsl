// SDF visualizer v2
// Allows a usually larger average marching step without missing thin surface
// Expected to be faster than v1

// orange-blue: SDF isosurfaces
// red-black: discontinuity (high numerical gradient)
// green-pink: surface gradient lower/higher than 1

#define PI 3.1415926
//#define ZERO min(iTime,0.)
#define ZERO 0.0


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

float sdEllipsoid(vec3 p, vec3 r) {
    float k1 = length(p/r);
    float k2 = length(p/(r*r));
    return k1*(k1-1.0)/k2;
}


float mapCapBack(vec3 p) {
    float r = length(p.xy);
    float a = atan(p.y, p.x);
    for (float i=0.;i<1.;i+=1.) {
        float f = 2.5*pow(2.,i);
        float dr = 0.3/f*sin(f*r)*cos(f*a);
        float da = 0.3/f*cos(f*a)*sin(f*r);
        r+=abs(dr), a+=da;
    }
    float f1 = cos(80.*a);
    float f2 = cos(160.*a);
    float f = 0.02*mix(f2,f1,clamp(2.0*r,0.,1.));
    float d = smax(f, p.z+0.03, 0.02);
    return d;
}

float mapCap1(vec3 p) {
    p.z -= 0.5/(1.0+dot(p.xy,p.xy))-0.4;
    p.z -= 0.05*sin(p.x)*cos(p.y);
    p.x += 0.02*cos(4.0*(p.y-0.1));
    float d = sdEllipsoid(p, vec3(1.0,1.0,max(0.2-0.1*length(p.xy),0.05)));
    d = smax(d, -mapCapBack(p), 0.01);
    return d;
}
float mapStem1(vec3 p) {
    p.xy += mix(vec2(0.6), vec2(0.0), 1.0-1.0/(1.0+exp(-4.0*(-p.z-2.4))));
    float r = mix(0.1, 0.14, smootherstep(0.2*p.z*p.z));
    float d = length(p.xy)-r;
    d = max(d, p.z);
    d = smax(d, -p.z-3.4, 0.2);
    d = smin(d, sdEllipsoid(p-vec3(0,0,-3.2),vec3(0.1,0.1,0.14)), 0.3);
    return d;
}
float mapFungi1(vec3 p) {
    p -= vec3(0.,0.,1.);
    float cap = mapCap1(p/1.1)*1.1;
    float stem = mapStem1(p);
    cap = smax(cap, -stem, 0.2);
    return smin(cap, stem, 0.1);
}

float mapCap2(vec3 p) {
    float r2 = dot(p.xy,p.xy);
    p.z -= 0.3/(1.0+pow(r2,4.))+0.2*exp(-10.*r2)-0.35;
    p.z += smin(0.1*r2*r2, 1.0, 0.5);
    p.z -= 0.05*cos(p.x)*sin(p.y);
    p.x += 0.02*cos(4.0*(p.y-0.1));
    float d = sdEllipsoid(p, vec3(1.0,1.0,max(0.2-0.1*length(p.xy),0.05)));
    d = smax(d, -mapCapBack(p), 0.01);
    return d;
}
float mapStem2(vec3 p) {
    p.xy += mix(vec2(-0.7), vec2(0.0), 1.0-1.0/(1.0+exp(-2.0*(-p.z-2.4))));
    float r = mix(0.1, 0.14, smootherstep(0.2*p.z*p.z)) - 0.02/(1.0+2.0*(p.z+1.6)*(p.z+1.6));
    float d = length(p.xy)-r;
    d = max(d, p.z);
    d = smax(d, -p.z-3.4, 0.2);
    d = smin(d, sdEllipsoid(p-vec3(0,0,-3.2),vec3(0.1,0.1,0.14)), 0.3);
    return d;
}
float mapFungi2(vec3 p) {
    p -= vec3(0.,0.,1.);
    float cap = mapCap2(p/1.0)*1.0;
    float stem = mapStem2(p);
    cap = smax(cap, -stem, 0.2);
    return smin(cap, stem, 0.1);
}

float sdf(vec3 p) {
    float d1 = mapFungi1((p-vec3(0.4,0.4,0.4))/1.1)*1.1;
    float d2 = mapFungi2((p-vec3(-0.9,-0.9,-0.2))/0.8)*0.8;
    float d = smin(d1, d2, 0.05);
    return d;
}


vec3 sdfGrad(in vec3 p, in float e) {
	float a = sdf(p+vec3(e,e,e));
	float b = sdf(p+vec3(e,-e,-e));
	float c = sdf(p+vec3(-e,e,-e));
	float d = sdf(p+vec3(-e,-e,e));
	return (.25/e)*vec3(a+b-c-d,a-b+c-d,a-b-c+d);
}




// raymarching parameters
#define BOX_RADIUS vec3(2.0, 2.0, 3.0)
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
    col = 0.2+0.05*grad.y+col*max(dot(normalize(grad), light),0.0);
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
