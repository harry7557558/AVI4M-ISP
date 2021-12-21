float fun(vec3 p) {
	vec3 u = p*p;
	float d = u.x+2.*u.y+u.z-1.;
	if (d>3.0) return d;
	return 4.*d*d-p.z*(5.*u.x*u.x-10.*u.x*u.z+u.z*u.z)-1.;
}
vec3 nGrad(vec3 p) {
	const float e = 1e-5;
	float a = fun(p+vec3(e,e,e));
	float b = fun(p+vec3(e,-e,-e));
	float c = fun(p+vec3(-e,e,-e));
	float d = fun(p+vec3(-e,-e,e));
	return vec3(a+b-c-d,a-b+c-d,a-b-c+d)*(.25/e);
}
vec4 map(vec3 p, bool col_required) {
	float sdf = fun(p) / length(vec4(nGrad(p), 0.1));
	sdf = abs(sdf) - 0.008;
	vec3 col = vec3(1.0,0.9,0.9)-sqrt(p.y*p.y+0.5)*vec3(0.1,0.4,0.9);
	return vec4(col, sdf);
}
