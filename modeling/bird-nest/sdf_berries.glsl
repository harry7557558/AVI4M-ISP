#include "common.glsl"


vec4 mapBerriesFruit(vec3 p, bool col_required) {
    float bound = length(p)-2.2, boundw = 0.5;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    p += vec3(0,0,0.5);
    float r = length(p.xy), a = atan(p.y, p.x);
    float x, y; vec3 q;
    // fruit
    q = vec3(r*cossin(asin(0.97*sin(2.5*a))/2.5), p.z);
    vec4 fruit = vec4(0,0,0, sdEllipsoid(q-vec3(0.4,0,0.68), vec3(0.8+0.1*p.z,1.1,0.95)));
    if (col_required) {
        float noise = SimplexNoise3D(4.0*p);
        fruit.xyz = mix(vec3(0.5,0.05,0.1),vec3(0.8,0.1,0.05),
            smootherstep(0.6*(q.z+1.0*(r-0.9)-0.5+smax(q.z-1.2,0.0,0.1)))
        ) + vec3(0.15)*(-noise+0.5);
        fruit.xyz = mix(fruit.xyz, vec3(0.8,0.8,0.0), 0.25+0.2*p.z);
        fruit.xyz = mix(fruit.xyz, vec3(0.8,0.0,0.5), 0.2);
    }
    q = vec3(r*cossin(asin(0.95*sin(2.5*a+0.8))/2.5), p.z);
    fruit = smin(fruit, vec4(  // hair
        mix(mix(vec3(0.8,0.75,0.0),vec3(0.9,0.85,0.5),smootherstep(r/0.3)), vec3(0.8,0.4,0.0), 0.2),
        sdEllipsoid(roty(0.2-0.05*cos(3.0*a))*(q-vec3(0.08,0,1.61)), vec3(0.2+0.03*sin(4.0*a),0.12,0.05))), 0.05);
    // sepal/stem
    q = vec3(r*cossin(asin(0.98*sin(2.5*a-1.2))/2.5), p.z);
    q -= vec3(0.2,0,0);
    vec3 br = 1.2*vec3(
        1.0+0.1*sin(3.0*a)-0.1*cos(2.0*a),
        1.25*(0.3-0.2*exp(-sqr(1.8*(q.x-1.6)))-0.2*exp(-sqr(3.0*(q.x-0.0)))),
        0.05);
    q = roty(smin(0.2*r,0.5,0.1))*(q-1.4*vec3(0.2,0,-0.15));
    q.z -= clamp(
        (0.1*pow(abs(q.x),4.0)+0.02*r*r*sin(3.0*a))*sin(a) + 0.02*sin(12.0*r)-0.1*sin(4.0*pow(r,1.3))
        - 0.05*exp(-sqr(12.0*q.y))*exp(-r),
        -0.5, 0.5);
    vec4 disk = vec4(mix(0.8*vec3(0.35,0.5,0.0),vec3(0.4,0.55,0.0),smootherstep(2.0*r-1.0)), sdEllipsoid(q, br));  // sepals
    q = p;
    vec4 stem = vec4(0.4,0.4,0.15, sdCapsule(q-vec3(0.05,0,-0.5), 0.2, 0.1));  // stem
    disk = smin(disk, stem, 0.1);
    vec4 d = smin(disk, fruit, 0.03);
    return d;
}
vec4 mapBerriesFruitT(vec3 p, bool col_required) {
    return mapBerriesFruit( p/0.5-vec3(0,0,1.0), col_required) * vec4(1,1,1,0.5);
}

vec4 mapBerriesLeaf(vec3 p, bool col_required) {
    float bound = sdEllipsoid(p-vec3(-0.1,0.0,0.1), vec3(2.2,1.5,0.6)), boundw = 0.3;
    if (bound > 0.0) return vec4(1,0,0, bound+boundw);  // clipping
    vec3 br;
    br.x = 1.6;
    br.y = 0.85 -0.45*exp(-sqr(1.2*(p.x-1.5))) -0.3*exp(-sqr(2.5*(p.x+1.6)));
    br.z = 0.08;
    float dqz = 0.0;
    float d2e = abs( br.y*sqrt(smax(br.x*br.x-p.x*p.x,0.0,0.1))/br.x - abs(p.y) );
    float veins = asin(0.9*sin((16.0+2.0*p.x)*(p.x-0.5*abs(p.y)-0.5*pow(abs(p.y),1.3)+0.2*pow(d2e,0.4))+sign(p.y)*0.4*PI));
    float veins_fade = (1.0-exp(-(6.0/br.y)*abs(p.y))) * (1.0-exp(-(2.0/br.y)*d2e)) * (exp(-0.2*p.x));
    dqz += 0.01 * veins*veins_fade;
    dqz += 0.05 * smax(veins_fade*(veins-0.9),0.0,0.1);
    float midrib = length(vec2(p.y,0.05))-1.0;
    float midrib_fade = (1.0-exp(-(4.0/br.y)*length(vec2(p.y,0.05)))) * (exp(-0.4*p.x));
    dqz += 0.05 * midrib*midrib_fade;
    vec3 q = p;
    q.z -= 0.5*(1.0-exp(-0.5*length(vec2(p.y,0.01))));
    q.z -= 0.1*cos(p.x);
    vec4 leaf = vec4(0,0,0, sdEllipsoid(q+vec3(0,0,dqz), br));
    if (col_required) leaf.xyz = mix(vec3(0.3,0.45,0.05), vec3(0.55,0.7,0.15), -0.5-20.0*dqz);
    vec4 stem = vec4(0.3,0.35,0.05, sdSegment(q-vec3(0,0,-0.00), vec3(-1.8,0,0), vec3(1.0,0,0)) - 0.06*exp(-0.2*abs(p.x+1.8)));
    return smin(leaf, stem, mix(0.05,0.001,clamp(p.x+1.8,0.,1.)));
}
vec4 mapBerriesLeafT(vec3 p, bool col_required) {
    return mapBerriesLeaf( p/0.7-vec3(1.8,0,0), col_required) * vec4(1,1,1,0.7);
}


vec4 map(vec3 p, bool col_required) {
    //return mapBerriesLeaf(p, col_required);
    return mapBerriesFruit(p, col_required);
}

