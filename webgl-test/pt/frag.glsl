#version 300 es
precision highp float;

out vec4 fragColor;

uniform vec2 iResolution;
uniform int iFrame;
uniform sampler2D iChannel0;

// baby path tracing

#define PI 3.1415926

uint seed = 0u;
uint randu() { return seed = seed * 1664525u + 1013904223u; }
float rand01() { return float(randu()) * (1./4294967296.); }


vec3 intersectRay(vec3 ro, vec3 rd) {
    // inspired by https://www.shadertoy.com/view/MdfBRX by BigWIngs
    float t = -(ro.z+2.0)/rd.z;
    vec3 p = ro+rd*t;
    vec3 col = rd.z > 0.0 ?
         mix(vec3(0.0), 0.3*vec3(0.2,0.5,0.8), pow(rd.z,1.0)) :  // sky
         abs(p.x)<12.0 ? 0.3*vec3(0.05,0.06,0.06) : 0.3*vec3(0.04,0.05,0.02);  // road
    // all distance are in meters
    for (float d=10.0; d<=120.0; d+=10.0) {  // street light
        t = -(ro.y+d)/rd.y;
        p = ro+rd*t;
        p.x = abs(p.x);
        if (length(p.xz-vec2(8.0,5.0))<0.15)
            col += 20.0*vec3(1.0,0.8,0.5);
    }
    for (float d=4.0; d<=54.0; d+=10.0) {  // head light
        t = -(ro.y+d)/rd.y;
        p = ro+rd*t;
        p.x = abs(p.x+4.0);
        if (length(p.xz-vec2(2.0,-0.5))<0.1)
            col += 15.0*vec3(0.8,0.8,1.0);
    }
    for (float d=6.0; d<=56.0; d+=10.0) {  // tail light
        t = -(ro.y+d)/rd.y;
        p = ro+rd*t;
        p.x = abs(p.x-4.0);
        if (length(p.xz-vec2(2.0,-0.5))<0.08)
            col += 10.0*vec3(1.0,0.1,0.0);
    }
    return col;
}

vec2 rndHeart() {  // x²-|x|y+y² < 1
    float u = 2.0*PI*rand01();
    float v = rand01();
    vec2 c = sqrt(v)*vec2(cos(u), sin(u));
    c = mat2(1.0,1.0,-0.577,0.577)*c;
    if (c.x<0.0) c.y=-c.y;
    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // random number seed
    vec3 p3 = fract(vec3(fragCoord/iResolution.xy, sin(0.001*float(iFrame))) * .1031);
    p3 += dot(p3, p3.zyx + 31.32);
    float h = fract((p3.x + p3.y) * p3.z);
    seed = uint(16777216.*h);

    // camera parameters
    const vec3 POS = vec3(3.0,3.0,0.0);
    const float SCALE = 1.5;  // larger = smaller (more view field)
    const float DIST = 0.2;  // larger = smaller
    const float VIEW_FIELD = 0.5;  // larger = larger + more perspective
    const float APERTURE = 0.01;  // larger = blurred

    // sample aperture shape
    vec2 rnd = rndHeart();
    float a = 0.2;
    rnd = mat2(cos(a),sin(a),-sin(a),cos(a))*rnd;

    // camera
    vec3 ro = POS+vec3(0,DIST,0);
    vec2 uv = SCALE*(2.0*(fragCoord.xy+vec2(rand01(),rand01())-0.5)/iResolution.xy-1.0);
    vec2 sc = iResolution.xy/length(iResolution.xy);
    vec2 offset = APERTURE*rnd;
    ro.xz += offset;
    vec3 rd = vec3(VIEW_FIELD*uv*sc+vec2(-0.2,0.0)-offset/DIST, -1.0).xzy;

    // calculate pixel color
    vec3 col = intersectRay(ro, rd);
    vec4 rgbn = texelFetch(iChannel0, ivec2(int(fragCoord.x), int(fragCoord.y)), 0);
    fragColor = vec4((rgbn.xyz*rgbn.w + col)/(rgbn.w+1.0), rgbn.w+1.0);
}

void main(void) {
    mainImage(fragColor, gl_FragCoord.xy);
}
