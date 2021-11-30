// https://www.shadertoy.com/view/sdVGWh

#define PI 3.1415926f

vec4 mapShell(vec3 p, bool req_color) {
	p -= vec3(0.7f, 0, 0);

	// r=exp(b*θ)
	const float b = 0.17f;

	// Catesian to cylindrical
	float r = length(p.xy());  // r
	float a = mix(0.0f, 0.45f, smoothstep(0.0f, 1.0f, 0.5f*(r - 0.6f)));  // rotate by this angle
	p.xy() = mat2(cos(a), -sin(a), sin(a), cos(a))*p.xy();  // rotation
	float t = atan(p.y, p.x);  // θ

	// shell opening, kill discontinuities of the spiral
	float ro = exp(b*PI);  // center of the "ring"
	float d = length(vec2(length(p.xz() - vec2(-ro, 0)) - ro, p.y));  // distance to the "ring"
	float u = t, dx = r - ro, dy = p.z;  // longitude and two numbers to determine latitude

	// spiral
	// r(n) = exp(b*(2.f*PI*n+t)), (x-r)^2+y^2=r^2, solve for n
	float n = (log((r*r + p.z*p.z) / (2.f*r)) / b - t) / (2.0f*PI);  // decimal n
	n = min(n, 0.0f);  // clamp to opening
	float n0 = floor(n), n1 = ceil(n);  // test two boundaries
	float r0 = exp(b*(2.f*PI*n0 + t)), r1 = exp(b*(2.f*PI*n1 + t));  // two r
	float d0 = abs(length(vec2(r - r0, p.z)) - r0);  // distance to inner
	float d1 = abs(length(vec2(r - r1, p.z)) - r1);  // distance to outer
	if (d0 < d) d = d0, u = 2.f*PI*n0 + t, dx = r - r0, dy = p.z;  // update distance
	if (d1 < d) d = d1, u = 2.f*PI*n1 + t, dx = r - r1, dy = p.z;  // update distance

	// septa/chambers
	const float f = 2.4f;  // "frequency" of chambers
	float s0 = t + 2.0f*PI*(n0 + 0.5f);  // longitude parameter
	float v = fract(n);  // 0-1, distance from inner circle
	float s = f * s0 + 1.0f*pow(0.25f - (v - 0.5f)*(v - 0.5f), 0.5f) + 0.5f*v;  // curve of septa
	s += pow(min(1.0f / (40.0f*length(vec2(v - 0.5f, p.z)) + 1.0f), 0.5f), 2.0f);  // hole on septa
	float sf = fract(s);  // periodic
	sf = s0 > -1.8f ? abs(s + 3.25f) :  // outer-most septa, possibly cause discontinuities
		min(sf, 1.0f - sf);  // inner septa
	float w = sf / f * exp(b*(s0 + PI));  // adjust distance field
	if (length(p*vec3(1, 1, 1.5f)) < 3.0f)  // prevent outer discontinuity
		d = min(d, 0.5f*w + 0.012f);  // union chambers

	d += 0.00012f*r*sin(200.f*u);  // geometric texture
	d = abs(d) - 0.8f*max(0.02f*pow(r, 0.4f), 0.02f);  // thickness of shell
	d = max(d, p.z);
	if (!req_color) return vec4(0, 0, 0, d);  // distance calculation finished

	// color
	vec3 col;
	v = atan(dy, dx);  // latitude parameter
	w = length(vec2(dx, dy)) / exp(b*u);  // section radius parameter
	for (float i = 0.f; i < 6.f; i += 1.f) {  // distort the parameters
		float f = pow(2.f, i);
		float du = 0.15f / f * sin(f*u)*cos(f*v);
		float dv = 0.15f / f * cos(f*u)*sin(f*v);
		u += du, v += dv;
	}
	float f1 = cos(50.f*u);  // middle stripes
	float f2 = cos(21.3f*u) + 0.1f;  // side stripes
	float tex = mix(f1, f2, 0.5f - 0.5f*tanh(1.0f - 3.0f*sin(v)*sin(v)))  // blend stripes
		+ 0.5f - 0.6f*cos(v);  // fading at sides
	tex += 0.5f + 0.5f*tanh(4.0f*(u - 2.0f));  // fading near opening
	col = n == 0.0f ? vec3(0.9f, 0.85f, 0.8f) : vec3(0.95f, 0.85f, 0.7f);  // base color, outer and inner
	if (w > 1.0f && w < 1.1f)  // on the surface of the shell
		col = (u - 0.3f*cos(v) < -2.6f ? 1.0f - 0.6f*min(exp(2.f + 0.5f*u), 1.0f) : 1.0f)  // black inside the opening
		* mix(vec3(0.6f, 0.3f, 0.2f), col, clamp(8.0f*tex + 0.5f, 0.f, 1.f));  // apply stripes

	return vec4(col, d);
}

vec4 map(vec3 p, bool col_required) {
	return mapShell(p / 0.7, col_required) * 0.7;
}
