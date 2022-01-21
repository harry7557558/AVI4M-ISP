#pragma once

#include <cmath>
#include <cstdlib>
#include <cstdint>

#define ENABLE_SWIZZLING 1

typedef uint32_t uint;

struct ivec2;
struct ivec3;
struct ivec4;
struct vec2;
struct vec3;
struct vec4;
struct mat2;
struct mat3;
struct mat4;

struct ivec2
{
	int x, y;
	explicit ivec2() {}
	explicit ivec2(const int &a) : x(a), y(a) {}
	explicit ivec2(const int &x, const int &y) : x((int)x), y((int)y) {}
	explicit ivec2(const vec2 &p);
	bool operator==(const ivec2 &v) const { return x == v.x && y == v.y; }
	bool operator!=(const ivec2 &v) const { return x != v.x || y != v.y; }
	ivec2 operator-() const { return ivec2(-x, -y); }
	ivec2 operator+(const ivec2 &v) const { return ivec2(x + v.x, y + v.y); }
	ivec2 operator-(const ivec2 &v) const { return ivec2(x - v.x, y - v.y); }
	ivec2 operator*(const ivec2 &v) const { return ivec2(x * v.x, y * v.y); }
	ivec2 operator/(const ivec2 &v) const { return ivec2(x / v.x, y / v.y); }
	ivec2 operator%(const ivec2 &v) const { return ivec2(x % v.x, y % v.y); }
	ivec2 operator+(const int &a) const { return ivec2(x + a, y + a); }
	ivec2 operator-(const int &a) const { return ivec2(x - a, y - a); }
	ivec2 operator*(const int &a) const { return ivec2(x * a, y * a); }
	ivec2 operator/(const int &a) const { return ivec2(x / a, y / a); }
	ivec2 operator%(const int &a) const { return ivec2(x % a, y % a); }
	ivec2 operator+=(const ivec2 &v) { x += v.x, y += v.y; return *this; }
	ivec2 operator-=(const ivec2 &v) { x -= v.x, y -= v.y; return *this; }
	ivec2 operator*=(const ivec2 &v) { x *= v.x, y *= v.y; return *this; }
	ivec2 operator/=(const ivec2 &v) { x /= v.x, y /= v.y; return *this; }
	ivec2 operator%=(const ivec2 &v) { x %= v.x, y %= v.y; return *this; }
	ivec2 operator+=(const int &a) { x += a, y += a; return *this; }
	ivec2 operator-=(const int &a) { x -= a, y -= a; return *this; }
	ivec2 operator*=(const int &a) { x *= a, y *= a; return *this; }
	ivec2 operator/=(const int &a) { x /= a, y /= a; return *this; }
	ivec2 operator%=(const int &a) { x %= a, y %= a; return *this; }
	friend ivec2 operator+(const int &a, const ivec2 &v) { return ivec2(a + v.x, a + v.y); }
	friend ivec2 operator-(const int &a, const ivec2 &v) { return ivec2(a - v.x, a - v.y); }
	friend ivec2 operator*(const int &a, const ivec2 &v) { return ivec2(a * v.x, a * v.y); }
	friend ivec2 operator/(const int &a, const ivec2 &v) { return ivec2(a / v.x, a / v.y); }
	friend ivec2 operator%(const int &a, const ivec2 &v) { return ivec2(a % v.x, a % v.y); }
};

struct ivec3
{
	int x, y, z;
	explicit ivec3() {}
	explicit ivec3(const int &a) : x(a), y(a), z(a) {}
	explicit ivec3(const int &x, const int &y, const int &z) : x((int)x), y((int)y), z((int)z) {}
	explicit ivec3(const ivec2 &xy, const int &z) : x((int)xy.x), y((int)xy.y), z((int)z) {}
	explicit ivec3(const int &x, const ivec2 &yz) : x((int)x), y((int)yz.x), z((int)yz.y) {}
	explicit ivec3(const vec3 &p);
	bool operator==(const ivec3 &v) const { return x == v.x && y == v.y && z == v.z; }
	bool operator!=(const ivec3 &v) const { return x != v.x || y != v.y || z != v.z; }
	ivec3 operator-() const { return ivec3(-x, -y, -z); }
	ivec3 operator+(const ivec3 &v) const { return ivec3(x + v.x, y + v.y, z + v.z); }
	ivec3 operator-(const ivec3 &v) const { return ivec3(x - v.x, y - v.y, z - v.z); }
	ivec3 operator*(const ivec3 &v) const { return ivec3(x * v.x, y * v.y, z * v.z); }
	ivec3 operator/(const ivec3 &v) const { return ivec3(x / v.x, y / v.y, z / v.z); }
	ivec3 operator%(const ivec3 &v) const { return ivec3(x % v.x, y % v.y, z % v.z); }
	ivec3 operator+(const int &a) const { return ivec3(x + a, y + a, z + a); }
	ivec3 operator-(const int &a) const { return ivec3(x - a, y - a, z - a); }
	ivec3 operator*(const int &a) const { return ivec3(x * a, y * a, z * a); }
	ivec3 operator/(const int &a) const { return ivec3(x / a, y / a, z / a); }
	ivec3 operator%(const int &a) const { return ivec3(x % a, y % a, z % a); }
	ivec3 operator+=(const ivec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
	ivec3 operator-=(const ivec3 &v) { x -= v.x, y -= v.y, z -= v.z; return *this; }
	ivec3 operator*=(const ivec3 &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	ivec3 operator/=(const ivec3 &v) { x /= v.x, y /= v.y, z /= v.z; return *this; }
	ivec3 operator%=(const ivec3 &v) { x %= v.x, y %= v.y, z %= v.z; return *this; }
	ivec3 operator+=(const int &a) { x += a, y += a, z += a; return *this; }
	ivec3 operator-=(const int &a) { x -= a, y -= a, z -= a; return *this; }
	ivec3 operator*=(const int &a) { x *= a, y *= a, z *= a; return *this; }
	ivec3 operator/=(const int &a) { x /= a, y /= a, z /= a; return *this; }
	ivec3 operator%=(const int &a) { x %= a, y %= a, z %= a; return *this; }
	friend ivec3 operator+(const int &a, const ivec3 &v) { return ivec3(a + v.x, a + v.y, a + v.z); }
	friend ivec3 operator-(const int &a, const ivec3 &v) { return ivec3(a - v.x, a - v.y, a - v.z); }
	friend ivec3 operator*(const int &a, const ivec3 &v) { return ivec3(a * v.x, a * v.y, a * v.z); }
	friend ivec3 operator/(const int &a, const ivec3 &v) { return ivec3(a / v.x, a / v.y, a / v.z); }
	friend ivec3 operator%(const int &a, const ivec3 &v) { return ivec3(a % v.x, a % v.y, a % v.z); }
};

struct ivec4
{
	int x, y, z, w;
	explicit ivec4() {}
	explicit ivec4(const int &a) : x(a), y(a), z(a), w(a) {}
	explicit ivec4(const int &x, const int &y, const int &z, const int &w) : x((int)x), y((int)y), z((int)z), w((int)w) {}
	explicit ivec4(const ivec2 &xy, const int &z, const int &w) : x((int)xy.x), y((int)xy.y), z((int)z), w((int)w) {}
	explicit ivec4(const int &x, const ivec2 &yz, const int &w) : x((int)x), y((int)yz.x), z((int)yz.y), w((int)w) {}
	explicit ivec4(const int &x, const int &y, const ivec2 &zw) : x((int)x), y((int)y), z((int)zw.x), w((int)zw.y) {}
	explicit ivec4(const ivec3 &xyz, const int &w) : x((int)xyz.x), y((int)xyz.y), z((int)xyz.z), w((int)w) {}
	explicit ivec4(const int &x, const ivec3 &yzw) : x((int)x), y((int)yzw.x), z((int)yzw.y), w((int)yzw.z) {}
	explicit ivec4(const vec4 &p);
	bool operator==(const ivec4 &v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
	bool operator!=(const ivec4 &v) const { return x != v.x || y != v.y || z != v.z || w != v.w; }
	ivec4 operator-() const { return ivec4(-x, -y, -z, -w); }
	ivec4 operator+(const ivec4 &v) const { return ivec4(x + v.x, y + v.y, z + v.z, w + v.w); }
	ivec4 operator-(const ivec4 &v) const { return ivec4(x - v.x, y - v.y, z - v.z, w - v.w); }
	ivec4 operator*(const ivec4 &v) const { return ivec4(x * v.x, y * v.y, z * v.z, w * v.w); }
	ivec4 operator/(const ivec4 &v) const { return ivec4(x / v.x, y / v.y, z / v.z, w / v.w); }
	ivec4 operator%(const ivec4 &v) const { return ivec4(x % v.x, y % v.y, z % v.z, w % v.w); }
	ivec4 operator+(const int &a) const { return ivec4(x + a, y + a, z + a, w + a); }
	ivec4 operator-(const int &a) const { return ivec4(x - a, y - a, z - a, w - a); }
	ivec4 operator*(const int &a) const { return ivec4(x * a, y * a, z * a, w * a); }
	ivec4 operator/(const int &a) const { return ivec4(x / a, y / a, z / a, w / a); }
	ivec4 operator%(const int &a) const { return ivec4(x % a, y % a, z % a, w % a); }
	ivec4 operator+=(const ivec4 &v) { x += v.x, y += v.y, z += v.z, w += v.w; return *this; }
	ivec4 operator-=(const ivec4 &v) { x -= v.x, y -= v.y, z -= v.z, w -= v.w; return *this; }
	ivec4 operator*=(const ivec4 &v) { x *= v.x, y *= v.y, z *= v.z, w *= v.w; return *this; }
	ivec4 operator/=(const ivec4 &v) { x /= v.x, y /= v.y, z /= v.z, w /= v.w; return *this; }
	ivec4 operator%=(const ivec4 &v) { x %= v.x, y %= v.y, z %= v.z, w %= v.w; return *this; }
	ivec4 operator+=(const int &a) { x += a, y += a, z += a, w += a; return *this; }
	ivec4 operator-=(const int &a) { x -= a, y -= a, z -= a, w -= a; return *this; }
	ivec4 operator*=(const int &a) { x *= a, y *= a, z *= a, w *= a; return *this; }
	ivec4 operator/=(const int &a) { x /= a, y /= a, z /= a, w /= a; return *this; }
	ivec4 operator%=(const int &a) { x %= a, y %= a, z %= a, w %= a; return *this; }
	friend ivec4 operator+(const int &a, const ivec4 &v) { return ivec4(a + v.x, a + v.y, a + v.z, a + v.w); }
	friend ivec4 operator-(const int &a, const ivec4 &v) { return ivec4(a - v.x, a - v.y, a - v.z, a - v.w); }
	friend ivec4 operator*(const int &a, const ivec4 &v) { return ivec4(a * v.x, a * v.y, a * v.z, a * v.w); }
	friend ivec4 operator/(const int &a, const ivec4 &v) { return ivec4(a / v.x, a / v.y, a / v.z, a / v.w); }
	friend ivec4 operator%(const int &a, const ivec4 &v) { return ivec4(a % v.x, a % v.y, a % v.z, a % v.w); }
};

struct vec2
{
	float x, y;
	explicit vec2() {}
	explicit vec2(const float &a) : x(a), y(a) {}
	explicit vec2(const float &x, const float &y) : x((float)x), y((float)y) {}
	explicit vec2(const ivec2 &v) : x((float)v.x), y((float)v.y) {}
	bool operator==(const vec2 &v) const { return x == v.x && y == v.y; }
	bool operator!=(const vec2 &v) const { return x != v.x || y != v.y; }
	vec2 operator-() const { return vec2(-x, -y); }
	vec2 operator+(const vec2 &v) const { return vec2(x + v.x, y + v.y); }
	vec2 operator-(const vec2 &v) const { return vec2(x - v.x, y - v.y); }
	vec2 operator*(const vec2 &v) const { return vec2(x * v.x, y * v.y); }
	vec2 operator/(const vec2 &v) const { return vec2(x / v.x, y / v.y); }
	vec2 operator+(const float &a) const { return vec2(x + a, y + a); }
	vec2 operator-(const float &a) const { return vec2(x - a, y - a); }
	vec2 operator*(const float &a) const { return vec2(x * a, y * a); }
	vec2 operator/(const float &a) const { return vec2(x / a, y / a); }
	vec2 operator+=(const vec2 &v) { x += v.x, y += v.y; return *this; }
	vec2 operator-=(const vec2 &v) { x -= v.x, y -= v.y; return *this; }
	vec2 operator*=(const vec2 &v) { x *= v.x, y *= v.y; return *this; }
	vec2 operator/=(const vec2 &v) { x /= v.x, y /= v.y; return *this; }
	vec2 operator+=(const float &a) { x += a, y += a; return *this; }
	vec2 operator-=(const float &a) { x -= a, y -= a; return *this; }
	vec2 operator*=(const float &a) { x *= a, y *= a; return *this; }
	vec2 operator/=(const float &a) { x /= a, y /= a; return *this; }
	friend vec2 operator+(const float &a, const vec2 &v) { return vec2(a + v.x, a + v.y); }
	friend vec2 operator-(const float &a, const vec2 &v) { return vec2(a - v.x, a - v.y); }
	friend vec2 operator*(const float &a, const vec2 &v) { return vec2(a * v.x, a * v.y); }
	friend vec2 operator/(const float &a, const vec2 &v) { return vec2(a / v.x, a / v.y); }
#if ENABLE_SWIZZLING
	vec2 xx() const { return vec2(x, x); }
	vec2 yx() const { return vec2(y, x); }
	vec2& xy() { return *(vec2*)&this->x; }
	vec2 yy() const { return vec2(y, y); }
	vec3 xxx() const;
	vec3 yxx() const;
	vec3 xyx() const;
	vec3 yyx() const;
	vec3 xxy() const;
	vec3 yxy() const;
	vec3 xyy() const;
	vec3 yyy() const;
	vec4 xxxx() const;
	vec4 yxxx() const;
	vec4 xyxx() const;
	vec4 yyxx() const;
	vec4 xxyx() const;
	vec4 yxyx() const;
	vec4 xyyx() const;
	vec4 yyyx() const;
	vec4 xxxy() const;
	vec4 yxxy() const;
	vec4 xyxy() const;
	vec4 yyxy() const;
	vec4 xxyy() const;
	vec4 yxyy() const;
	vec4 xyyy() const;
	vec4 yyyy() const;
#endif
};

struct vec3
{
	float x, y, z;
	explicit vec3() {}
	explicit vec3(const float &a) : x(a), y(a), z(a) {}
	explicit vec3(const float &x, const float &y, const float &z) : x((float)x), y((float)y), z((float)z) {}
	explicit vec3(const vec2 &xy, const float &z) : x((float)xy.x), y((float)xy.y), z((float)z) {}
	explicit vec3(const float &x, const vec2 &yz) : x((float)x), y((float)yz.x), z((float)yz.y) {}
	explicit vec3(const ivec3 &v) : x((float)v.x), y((float)v.y), z((float)v.z) {}
	bool operator==(const vec3 &v) const { return x == v.x && y == v.y && z == v.z; }
	bool operator!=(const vec3 &v) const { return x != v.x || y != v.y || z != v.z; }
	vec3 operator-() const { return vec3(-x, -y, -z); }
	vec3 operator+(const vec3 &v) const { return vec3(x + v.x, y + v.y, z + v.z); }
	vec3 operator-(const vec3 &v) const { return vec3(x - v.x, y - v.y, z - v.z); }
	vec3 operator*(const vec3 &v) const { return vec3(x * v.x, y * v.y, z * v.z); }
	vec3 operator/(const vec3 &v) const { return vec3(x / v.x, y / v.y, z / v.z); }
	vec3 operator+(const float &a) const { return vec3(x + a, y + a, z + a); }
	vec3 operator-(const float &a) const { return vec3(x - a, y - a, z - a); }
	vec3 operator*(const float &a) const { return vec3(x * a, y * a, z * a); }
	vec3 operator/(const float &a) const { return vec3(x / a, y / a, z / a); }
	vec3 operator+=(const vec3 &v) { x += v.x, y += v.y, z += v.z; return *this; }
	vec3 operator-=(const vec3 &v) { x -= v.x, y -= v.y, z -= v.z; return *this; }
	vec3 operator*=(const vec3 &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	vec3 operator/=(const vec3 &v) { x /= v.x, y /= v.y, z /= v.z; return *this; }
	vec3 operator+=(const float &a) { x += a, y += a, z += a; return *this; }
	vec3 operator-=(const float &a) { x -= a, y -= a, z -= a; return *this; }
	vec3 operator*=(const float &a) { x *= a, y *= a, z *= a; return *this; }
	vec3 operator/=(const float &a) { x /= a, y /= a, z /= a; return *this; }
	friend vec3 operator+(const float &a, const vec3 &v) { return vec3(a + v.x, a + v.y, a + v.z); }
	friend vec3 operator-(const float &a, const vec3 &v) { return vec3(a - v.x, a - v.y, a - v.z); }
	friend vec3 operator*(const float &a, const vec3 &v) { return vec3(a * v.x, a * v.y, a * v.z); }
	friend vec3 operator/(const float &a, const vec3 &v) { return vec3(a / v.x, a / v.y, a / v.z); }
#if ENABLE_SWIZZLING
	vec2 xx() const { return vec2(x, x); }
	vec2 yx() const { return vec2(y, x); }
	vec2 zx() const { return vec2(z, x); }
	vec2& xy() { return *(vec2*)&this->x; }
	vec2 yy() const { return vec2(y, y); }
	vec2 zy() const { return vec2(z, y); }
	vec2 xz() const { return vec2(x, z); }
	vec2& yz() const { return *(vec2*)&this->y; }
	vec2 zz() const { return vec2(z, z); }
	vec3 xxx() const { return vec3(x, x, x); }
	vec3 yxx() const { return vec3(y, x, x); }
	vec3 zxx() const { return vec3(z, x, x); }
	vec3 xyx() const { return vec3(x, y, x); }
	vec3 yyx() const { return vec3(y, y, x); }
	vec3 zyx() const { return vec3(z, y, x); }
	vec3 xzx() const { return vec3(x, z, x); }
	vec3 yzx() const { return vec3(y, z, x); }
	vec3 zzx() const { return vec3(z, z, x); }
	vec3 xxy() const { return vec3(x, x, y); }
	vec3 yxy() const { return vec3(y, x, y); }
	vec3 zxy() const { return vec3(z, x, y); }
	vec3 xyy() const { return vec3(x, y, y); }
	vec3 yyy() const { return vec3(y, y, y); }
	vec3 zyy() const { return vec3(z, y, y); }
	vec3 xzy() const { return vec3(x, z, y); }
	vec3 yzy() const { return vec3(y, z, y); }
	vec3 zzy() const { return vec3(z, z, y); }
	vec3 xxz() const { return vec3(x, x, z); }
	vec3 yxz() const { return vec3(y, x, z); }
	vec3 zxz() const { return vec3(z, x, z); }
	vec3& xyz() const { return *(vec3*)&this->x; }
	vec3 yyz() const { return vec3(y, y, z); }
	vec3 zyz() const { return vec3(z, y, z); }
	vec3 xzz() const { return vec3(x, z, z); }
	vec3 yzz() const { return vec3(y, z, z); }
	vec3 zzz() const { return vec3(z, z, z); }
	vec4 xxxx() const;
	vec4 yxxx() const;
	vec4 zxxx() const;
	vec4 xyxx() const;
	vec4 yyxx() const;
	vec4 zyxx() const;
	vec4 xzxx() const;
	vec4 yzxx() const;
	vec4 zzxx() const;
	vec4 xxyx() const;
	vec4 yxyx() const;
	vec4 zxyx() const;
	vec4 xyyx() const;
	vec4 yyyx() const;
	vec4 zyyx() const;
	vec4 xzyx() const;
	vec4 yzyx() const;
	vec4 zzyx() const;
	vec4 xxzx() const;
	vec4 yxzx() const;
	vec4 zxzx() const;
	vec4 xyzx() const;
	vec4 yyzx() const;
	vec4 zyzx() const;
	vec4 xzzx() const;
	vec4 yzzx() const;
	vec4 zzzx() const;
	vec4 xxxy() const;
	vec4 yxxy() const;
	vec4 zxxy() const;
	vec4 xyxy() const;
	vec4 yyxy() const;
	vec4 zyxy() const;
	vec4 xzxy() const;
	vec4 yzxy() const;
	vec4 zzxy() const;
	vec4 xxyy() const;
	vec4 yxyy() const;
	vec4 zxyy() const;
	vec4 xyyy() const;
	vec4 yyyy() const;
	vec4 zyyy() const;
	vec4 xzyy() const;
	vec4 yzyy() const;
	vec4 zzyy() const;
	vec4 xxzy() const;
	vec4 yxzy() const;
	vec4 zxzy() const;
	vec4 xyzy() const;
	vec4 yyzy() const;
	vec4 zyzy() const;
	vec4 xzzy() const;
	vec4 yzzy() const;
	vec4 zzzy() const;
	vec4 xxxz() const;
	vec4 yxxz() const;
	vec4 zxxz() const;
	vec4 xyxz() const;
	vec4 yyxz() const;
	vec4 zyxz() const;
	vec4 xzxz() const;
	vec4 yzxz() const;
	vec4 zzxz() const;
	vec4 xxyz() const;
	vec4 yxyz() const;
	vec4 zxyz() const;
	vec4 xyyz() const;
	vec4 yyyz() const;
	vec4 zyyz() const;
	vec4 xzyz() const;
	vec4 yzyz() const;
	vec4 zzyz() const;
	vec4 xxzz() const;
	vec4 yxzz() const;
	vec4 zxzz() const;
	vec4 xyzz() const;
	vec4 yyzz() const;
	vec4 zyzz() const;
	vec4 xzzz() const;
	vec4 yzzz() const;
	vec4 zzzz() const;
#endif
};

struct vec4
{
	float x, y, z, w;
	explicit vec4() {}
	explicit vec4(const float &a) : x(a), y(a), z(a), w(a) {}
	explicit vec4(const float &x, const float &y, const float &z, const float &w) : x((float)x), y((float)y), z((float)z), w((float)w) {}
	explicit vec4(const vec2 &xy, const float &z, const float &w) : x((float)xy.x), y((float)xy.y), z((float)z), w((float)w) {}
	explicit vec4(const float &x, const vec2 &yz, const float &w) : x((float)x), y((float)yz.x), z((float)yz.y), w((float)w) {}
	explicit vec4(const float &x, const float &y, const vec2 &zw) : x((float)x), y((float)y), z((float)zw.x), w((float)zw.y) {}
	explicit vec4(const vec3 &xyz, const float &w) : x((float)xyz.x), y((float)xyz.y), z((float)xyz.z), w((float)w) {}
	explicit vec4(const float &x, const vec3 &yzw) : x((float)x), y((float)yzw.x), z((float)yzw.y), w((float)yzw.z) {}
	explicit vec4(const ivec4 &v) : x((float)v.x), y((float)v.y), z((float)v.z), w((float)v.w) {}
	bool operator==(const vec4 &v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
	bool operator!=(const vec4 &v) const { return x != v.x || y != v.y || z != v.z || w != v.w; }
	vec4 operator-() const { return vec4(-x, -y, -z, -w); }
	vec4 operator+(const vec4 &v) const { return vec4(x + v.x, y + v.y, z + v.z, w + v.w); }
	vec4 operator-(const vec4 &v) const { return vec4(x - v.x, y - v.y, z - v.z, w - v.w); }
	vec4 operator*(const vec4 &v) const { return vec4(x * v.x, y * v.y, z * v.z, w * v.w); }
	vec4 operator/(const vec4 &v) const { return vec4(x / v.x, y / v.y, z / v.z, w / v.w); }
	vec4 operator+(const float &a) const { return vec4(x + a, y + a, z + a, w + a); }
	vec4 operator-(const float &a) const { return vec4(x - a, y - a, z - a, w - a); }
	vec4 operator*(const float &a) const { return vec4(x * a, y * a, z * a, w * a); }
	vec4 operator/(const float &a) const { return vec4(x / a, y / a, z / a, w / a); }
	vec4 operator+=(const vec4 &v) { x += v.x, y += v.y, z += v.z, w += v.w; return *this; }
	vec4 operator-=(const vec4 &v) { x -= v.x, y -= v.y, z -= v.z, w -= v.w; return *this; }
	vec4 operator*=(const vec4 &v) { x *= v.x, y *= v.y, z *= v.z, w *= v.w; return *this; }
	vec4 operator/=(const vec4 &v) { x /= v.x, y /= v.y, z /= v.z, w /= v.w; return *this; }
	vec4 operator+=(const float &a) { x += a, y += a, z += a, w += a; return *this; }
	vec4 operator-=(const float &a) { x -= a, y -= a, z -= a, w -= a; return *this; }
	vec4 operator*=(const float &a) { x *= a, y *= a, z *= a, w *= a; return *this; }
	vec4 operator/=(const float &a) { x /= a, y /= a, z /= a, w /= a; return *this; }
	friend vec4 operator+(const float &a, const vec4 &v) { return vec4(a + v.x, a + v.y, a + v.z, a + v.w); }
	friend vec4 operator-(const float &a, const vec4 &v) { return vec4(a - v.x, a - v.y, a - v.z, a - v.w); }
	friend vec4 operator*(const float &a, const vec4 &v) { return vec4(a * v.x, a * v.y, a * v.z, a * v.w); }
	friend vec4 operator/(const float &a, const vec4 &v) { return vec4(a / v.x, a / v.y, a / v.z, a / v.w); }
#if ENABLE_SWIZZLING
	vec2 xx() const { return vec2(x, x); }
	vec2 yx() const { return vec2(y, x); }
	vec2 zx() const { return vec2(z, x); }
	vec2 wx() const { return vec2(w, x); }
	vec2& xy() { return *(vec2*)&this->x; }
	vec2 yy() const { return vec2(y, y); }
	vec2 zy() const { return vec2(z, y); }
	vec2 wy() const { return vec2(w, y); }
	vec2 xz() const { return vec2(x, z); }
	vec2& yz() { return *(vec2*)&this->y; }
	vec2 zz() const { return vec2(z, z); }
	vec2 wz() const { return vec2(w, z); }
	vec2 xw() const { return vec2(x, w); }
	vec2 yw() const { return vec2(y, w); }
	vec2& zw() { return *(vec2*)&this->z; }
	vec2 ww() const { return vec2(w, w); }
	vec3 xxx() const { return vec3(x, x, x); }
	vec3 yxx() const { return vec3(y, x, x); }
	vec3 zxx() const { return vec3(z, x, x); }
	vec3 wxx() const { return vec3(w, x, x); }
	vec3 xyx() const { return vec3(x, y, x); }
	vec3 yyx() const { return vec3(y, y, x); }
	vec3 zyx() const { return vec3(z, y, x); }
	vec3 wyx() const { return vec3(w, y, x); }
	vec3 xzx() const { return vec3(x, z, x); }
	vec3 yzx() const { return vec3(y, z, x); }
	vec3 zzx() const { return vec3(z, z, x); }
	vec3 wzx() const { return vec3(w, z, x); }
	vec3 xwx() const { return vec3(x, w, x); }
	vec3 ywx() const { return vec3(y, w, x); }
	vec3 zwx() const { return vec3(z, w, x); }
	vec3 wwx() const { return vec3(w, w, x); }
	vec3 xxy() const { return vec3(x, x, y); }
	vec3 yxy() const { return vec3(y, x, y); }
	vec3 zxy() const { return vec3(z, x, y); }
	vec3 wxy() const { return vec3(w, x, y); }
	vec3 xyy() const { return vec3(x, y, y); }
	vec3 yyy() const { return vec3(y, y, y); }
	vec3 zyy() const { return vec3(z, y, y); }
	vec3 wyy() const { return vec3(w, y, y); }
	vec3 xzy() const { return vec3(x, z, y); }
	vec3 yzy() const { return vec3(y, z, y); }
	vec3 zzy() const { return vec3(z, z, y); }
	vec3 wzy() const { return vec3(w, z, y); }
	vec3 xwy() const { return vec3(x, w, y); }
	vec3 ywy() const { return vec3(y, w, y); }
	vec3 zwy() const { return vec3(z, w, y); }
	vec3 wwy() const { return vec3(w, w, y); }
	vec3 xxz() const { return vec3(x, x, z); }
	vec3 yxz() const { return vec3(y, x, z); }
	vec3 zxz() const { return vec3(z, x, z); }
	vec3 wxz() const { return vec3(w, x, z); }
	vec3& xyz() { return *(vec3*)&this->x; }
	vec3 yyz() const { return vec3(y, y, z); }
	vec3 zyz() const { return vec3(z, y, z); }
	vec3 wyz() const { return vec3(w, y, z); }
	vec3 xzz() const { return vec3(x, z, z); }
	vec3 yzz() const { return vec3(y, z, z); }
	vec3 zzz() const { return vec3(z, z, z); }
	vec3 wzz() const { return vec3(w, z, z); }
	vec3 xwz() const { return vec3(x, w, z); }
	vec3 ywz() const { return vec3(y, w, z); }
	vec3 zwz() const { return vec3(z, w, z); }
	vec3 wwz() const { return vec3(w, w, z); }
	vec3 xxw() const { return vec3(x, x, w); }
	vec3 yxw() const { return vec3(y, x, w); }
	vec3 zxw() const { return vec3(z, x, w); }
	vec3 wxw() const { return vec3(w, x, w); }
	vec3 xyw() const { return vec3(x, y, w); }
	vec3 yyw() const { return vec3(y, y, w); }
	vec3 zyw() const { return vec3(z, y, w); }
	vec3 wyw() const { return vec3(w, y, w); }
	vec3 xzw() const { return vec3(x, z, w); }
	vec3& yzw() { return *(vec3*)&this->y; }
	vec3 zzw() const { return vec3(z, z, w); }
	vec3 wzw() const { return vec3(w, z, w); }
	vec3 xww() const { return vec3(x, w, w); }
	vec3 yww() const { return vec3(y, w, w); }
	vec3 zww() const { return vec3(z, w, w); }
	vec3 www() const { return vec3(w, w, w); }
	vec4 xxxx() const { return vec4(x, x, x, x); }
	vec4 yxxx() const { return vec4(y, x, x, x); }
	vec4 zxxx() const { return vec4(z, x, x, x); }
	vec4 wxxx() const { return vec4(w, x, x, x); }
	vec4 xyxx() const { return vec4(x, y, x, x); }
	vec4 yyxx() const { return vec4(y, y, x, x); }
	vec4 zyxx() const { return vec4(z, y, x, x); }
	vec4 wyxx() const { return vec4(w, y, x, x); }
	vec4 xzxx() const { return vec4(x, z, x, x); }
	vec4 yzxx() const { return vec4(y, z, x, x); }
	vec4 zzxx() const { return vec4(z, z, x, x); }
	vec4 wzxx() const { return vec4(w, z, x, x); }
	vec4 xwxx() const { return vec4(x, w, x, x); }
	vec4 ywxx() const { return vec4(y, w, x, x); }
	vec4 zwxx() const { return vec4(z, w, x, x); }
	vec4 wwxx() const { return vec4(w, w, x, x); }
	vec4 xxyx() const { return vec4(x, x, y, x); }
	vec4 yxyx() const { return vec4(y, x, y, x); }
	vec4 zxyx() const { return vec4(z, x, y, x); }
	vec4 wxyx() const { return vec4(w, x, y, x); }
	vec4 xyyx() const { return vec4(x, y, y, x); }
	vec4 yyyx() const { return vec4(y, y, y, x); }
	vec4 zyyx() const { return vec4(z, y, y, x); }
	vec4 wyyx() const { return vec4(w, y, y, x); }
	vec4 xzyx() const { return vec4(x, z, y, x); }
	vec4 yzyx() const { return vec4(y, z, y, x); }
	vec4 zzyx() const { return vec4(z, z, y, x); }
	vec4 wzyx() const { return vec4(w, z, y, x); }
	vec4 xwyx() const { return vec4(x, w, y, x); }
	vec4 ywyx() const { return vec4(y, w, y, x); }
	vec4 zwyx() const { return vec4(z, w, y, x); }
	vec4 wwyx() const { return vec4(w, w, y, x); }
	vec4 xxzx() const { return vec4(x, x, z, x); }
	vec4 yxzx() const { return vec4(y, x, z, x); }
	vec4 zxzx() const { return vec4(z, x, z, x); }
	vec4 wxzx() const { return vec4(w, x, z, x); }
	vec4 xyzx() const { return vec4(x, y, z, x); }
	vec4 yyzx() const { return vec4(y, y, z, x); }
	vec4 zyzx() const { return vec4(z, y, z, x); }
	vec4 wyzx() const { return vec4(w, y, z, x); }
	vec4 xzzx() const { return vec4(x, z, z, x); }
	vec4 yzzx() const { return vec4(y, z, z, x); }
	vec4 zzzx() const { return vec4(z, z, z, x); }
	vec4 wzzx() const { return vec4(w, z, z, x); }
	vec4 xwzx() const { return vec4(x, w, z, x); }
	vec4 ywzx() const { return vec4(y, w, z, x); }
	vec4 zwzx() const { return vec4(z, w, z, x); }
	vec4 wwzx() const { return vec4(w, w, z, x); }
	vec4 xxwx() const { return vec4(x, x, w, x); }
	vec4 yxwx() const { return vec4(y, x, w, x); }
	vec4 zxwx() const { return vec4(z, x, w, x); }
	vec4 wxwx() const { return vec4(w, x, w, x); }
	vec4 xywx() const { return vec4(x, y, w, x); }
	vec4 yywx() const { return vec4(y, y, w, x); }
	vec4 zywx() const { return vec4(z, y, w, x); }
	vec4 wywx() const { return vec4(w, y, w, x); }
	vec4 xzwx() const { return vec4(x, z, w, x); }
	vec4 yzwx() const { return vec4(y, z, w, x); }
	vec4 zzwx() const { return vec4(z, z, w, x); }
	vec4 wzwx() const { return vec4(w, z, w, x); }
	vec4 xwwx() const { return vec4(x, w, w, x); }
	vec4 ywwx() const { return vec4(y, w, w, x); }
	vec4 zwwx() const { return vec4(z, w, w, x); }
	vec4 wwwx() const { return vec4(w, w, w, x); }
	vec4 xxxy() const { return vec4(x, x, x, y); }
	vec4 yxxy() const { return vec4(y, x, x, y); }
	vec4 zxxy() const { return vec4(z, x, x, y); }
	vec4 wxxy() const { return vec4(w, x, x, y); }
	vec4 xyxy() const { return vec4(x, y, x, y); }
	vec4 yyxy() const { return vec4(y, y, x, y); }
	vec4 zyxy() const { return vec4(z, y, x, y); }
	vec4 wyxy() const { return vec4(w, y, x, y); }
	vec4 xzxy() const { return vec4(x, z, x, y); }
	vec4 yzxy() const { return vec4(y, z, x, y); }
	vec4 zzxy() const { return vec4(z, z, x, y); }
	vec4 wzxy() const { return vec4(w, z, x, y); }
	vec4 xwxy() const { return vec4(x, w, x, y); }
	vec4 ywxy() const { return vec4(y, w, x, y); }
	vec4 zwxy() const { return vec4(z, w, x, y); }
	vec4 wwxy() const { return vec4(w, w, x, y); }
	vec4 xxyy() const { return vec4(x, x, y, y); }
	vec4 yxyy() const { return vec4(y, x, y, y); }
	vec4 zxyy() const { return vec4(z, x, y, y); }
	vec4 wxyy() const { return vec4(w, x, y, y); }
	vec4 xyyy() const { return vec4(x, y, y, y); }
	vec4 yyyy() const { return vec4(y, y, y, y); }
	vec4 zyyy() const { return vec4(z, y, y, y); }
	vec4 wyyy() const { return vec4(w, y, y, y); }
	vec4 xzyy() const { return vec4(x, z, y, y); }
	vec4 yzyy() const { return vec4(y, z, y, y); }
	vec4 zzyy() const { return vec4(z, z, y, y); }
	vec4 wzyy() const { return vec4(w, z, y, y); }
	vec4 xwyy() const { return vec4(x, w, y, y); }
	vec4 ywyy() const { return vec4(y, w, y, y); }
	vec4 zwyy() const { return vec4(z, w, y, y); }
	vec4 wwyy() const { return vec4(w, w, y, y); }
	vec4 xxzy() const { return vec4(x, x, z, y); }
	vec4 yxzy() const { return vec4(y, x, z, y); }
	vec4 zxzy() const { return vec4(z, x, z, y); }
	vec4 wxzy() const { return vec4(w, x, z, y); }
	vec4 xyzy() const { return vec4(x, y, z, y); }
	vec4 yyzy() const { return vec4(y, y, z, y); }
	vec4 zyzy() const { return vec4(z, y, z, y); }
	vec4 wyzy() const { return vec4(w, y, z, y); }
	vec4 xzzy() const { return vec4(x, z, z, y); }
	vec4 yzzy() const { return vec4(y, z, z, y); }
	vec4 zzzy() const { return vec4(z, z, z, y); }
	vec4 wzzy() const { return vec4(w, z, z, y); }
	vec4 xwzy() const { return vec4(x, w, z, y); }
	vec4 ywzy() const { return vec4(y, w, z, y); }
	vec4 zwzy() const { return vec4(z, w, z, y); }
	vec4 wwzy() const { return vec4(w, w, z, y); }
	vec4 xxwy() const { return vec4(x, x, w, y); }
	vec4 yxwy() const { return vec4(y, x, w, y); }
	vec4 zxwy() const { return vec4(z, x, w, y); }
	vec4 wxwy() const { return vec4(w, x, w, y); }
	vec4 xywy() const { return vec4(x, y, w, y); }
	vec4 yywy() const { return vec4(y, y, w, y); }
	vec4 zywy() const { return vec4(z, y, w, y); }
	vec4 wywy() const { return vec4(w, y, w, y); }
	vec4 xzwy() const { return vec4(x, z, w, y); }
	vec4 yzwy() const { return vec4(y, z, w, y); }
	vec4 zzwy() const { return vec4(z, z, w, y); }
	vec4 wzwy() const { return vec4(w, z, w, y); }
	vec4 xwwy() const { return vec4(x, w, w, y); }
	vec4 ywwy() const { return vec4(y, w, w, y); }
	vec4 zwwy() const { return vec4(z, w, w, y); }
	vec4 wwwy() const { return vec4(w, w, w, y); }
	vec4 xxxz() const { return vec4(x, x, x, z); }
	vec4 yxxz() const { return vec4(y, x, x, z); }
	vec4 zxxz() const { return vec4(z, x, x, z); }
	vec4 wxxz() const { return vec4(w, x, x, z); }
	vec4 xyxz() const { return vec4(x, y, x, z); }
	vec4 yyxz() const { return vec4(y, y, x, z); }
	vec4 zyxz() const { return vec4(z, y, x, z); }
	vec4 wyxz() const { return vec4(w, y, x, z); }
	vec4 xzxz() const { return vec4(x, z, x, z); }
	vec4 yzxz() const { return vec4(y, z, x, z); }
	vec4 zzxz() const { return vec4(z, z, x, z); }
	vec4 wzxz() const { return vec4(w, z, x, z); }
	vec4 xwxz() const { return vec4(x, w, x, z); }
	vec4 ywxz() const { return vec4(y, w, x, z); }
	vec4 zwxz() const { return vec4(z, w, x, z); }
	vec4 wwxz() const { return vec4(w, w, x, z); }
	vec4 xxyz() const { return vec4(x, x, y, z); }
	vec4 yxyz() const { return vec4(y, x, y, z); }
	vec4 zxyz() const { return vec4(z, x, y, z); }
	vec4 wxyz() const { return vec4(w, x, y, z); }
	vec4 xyyz() const { return vec4(x, y, y, z); }
	vec4 yyyz() const { return vec4(y, y, y, z); }
	vec4 zyyz() const { return vec4(z, y, y, z); }
	vec4 wyyz() const { return vec4(w, y, y, z); }
	vec4 xzyz() const { return vec4(x, z, y, z); }
	vec4 yzyz() const { return vec4(y, z, y, z); }
	vec4 zzyz() const { return vec4(z, z, y, z); }
	vec4 wzyz() const { return vec4(w, z, y, z); }
	vec4 xwyz() const { return vec4(x, w, y, z); }
	vec4 ywyz() const { return vec4(y, w, y, z); }
	vec4 zwyz() const { return vec4(z, w, y, z); }
	vec4 wwyz() const { return vec4(w, w, y, z); }
	vec4 xxzz() const { return vec4(x, x, z, z); }
	vec4 yxzz() const { return vec4(y, x, z, z); }
	vec4 zxzz() const { return vec4(z, x, z, z); }
	vec4 wxzz() const { return vec4(w, x, z, z); }
	vec4 xyzz() const { return vec4(x, y, z, z); }
	vec4 yyzz() const { return vec4(y, y, z, z); }
	vec4 zyzz() const { return vec4(z, y, z, z); }
	vec4 wyzz() const { return vec4(w, y, z, z); }
	vec4 xzzz() const { return vec4(x, z, z, z); }
	vec4 yzzz() const { return vec4(y, z, z, z); }
	vec4 zzzz() const { return vec4(z, z, z, z); }
	vec4 wzzz() const { return vec4(w, z, z, z); }
	vec4 xwzz() const { return vec4(x, w, z, z); }
	vec4 ywzz() const { return vec4(y, w, z, z); }
	vec4 zwzz() const { return vec4(z, w, z, z); }
	vec4 wwzz() const { return vec4(w, w, z, z); }
	vec4 xxwz() const { return vec4(x, x, w, z); }
	vec4 yxwz() const { return vec4(y, x, w, z); }
	vec4 zxwz() const { return vec4(z, x, w, z); }
	vec4 wxwz() const { return vec4(w, x, w, z); }
	vec4 xywz() const { return vec4(x, y, w, z); }
	vec4 yywz() const { return vec4(y, y, w, z); }
	vec4 zywz() const { return vec4(z, y, w, z); }
	vec4 wywz() const { return vec4(w, y, w, z); }
	vec4 xzwz() const { return vec4(x, z, w, z); }
	vec4 yzwz() const { return vec4(y, z, w, z); }
	vec4 zzwz() const { return vec4(z, z, w, z); }
	vec4 wzwz() const { return vec4(w, z, w, z); }
	vec4 xwwz() const { return vec4(x, w, w, z); }
	vec4 ywwz() const { return vec4(y, w, w, z); }
	vec4 zwwz() const { return vec4(z, w, w, z); }
	vec4 wwwz() const { return vec4(w, w, w, z); }
	vec4 xxxw() const { return vec4(x, x, x, w); }
	vec4 yxxw() const { return vec4(y, x, x, w); }
	vec4 zxxw() const { return vec4(z, x, x, w); }
	vec4 wxxw() const { return vec4(w, x, x, w); }
	vec4 xyxw() const { return vec4(x, y, x, w); }
	vec4 yyxw() const { return vec4(y, y, x, w); }
	vec4 zyxw() const { return vec4(z, y, x, w); }
	vec4 wyxw() const { return vec4(w, y, x, w); }
	vec4 xzxw() const { return vec4(x, z, x, w); }
	vec4 yzxw() const { return vec4(y, z, x, w); }
	vec4 zzxw() const { return vec4(z, z, x, w); }
	vec4 wzxw() const { return vec4(w, z, x, w); }
	vec4 xwxw() const { return vec4(x, w, x, w); }
	vec4 ywxw() const { return vec4(y, w, x, w); }
	vec4 zwxw() const { return vec4(z, w, x, w); }
	vec4 wwxw() const { return vec4(w, w, x, w); }
	vec4 xxyw() const { return vec4(x, x, y, w); }
	vec4 yxyw() const { return vec4(y, x, y, w); }
	vec4 zxyw() const { return vec4(z, x, y, w); }
	vec4 wxyw() const { return vec4(w, x, y, w); }
	vec4 xyyw() const { return vec4(x, y, y, w); }
	vec4 yyyw() const { return vec4(y, y, y, w); }
	vec4 zyyw() const { return vec4(z, y, y, w); }
	vec4 wyyw() const { return vec4(w, y, y, w); }
	vec4 xzyw() const { return vec4(x, z, y, w); }
	vec4 yzyw() const { return vec4(y, z, y, w); }
	vec4 zzyw() const { return vec4(z, z, y, w); }
	vec4 wzyw() const { return vec4(w, z, y, w); }
	vec4 xwyw() const { return vec4(x, w, y, w); }
	vec4 ywyw() const { return vec4(y, w, y, w); }
	vec4 zwyw() const { return vec4(z, w, y, w); }
	vec4 wwyw() const { return vec4(w, w, y, w); }
	vec4 xxzw() const { return vec4(x, x, z, w); }
	vec4 yxzw() const { return vec4(y, x, z, w); }
	vec4 zxzw() const { return vec4(z, x, z, w); }
	vec4 wxzw() const { return vec4(w, x, z, w); }
	vec4 xyzw() const { return vec4(x, y, z, w); }
	vec4 yyzw() const { return vec4(y, y, z, w); }
	vec4 zyzw() const { return vec4(z, y, z, w); }
	vec4 wyzw() const { return vec4(w, y, z, w); }
	vec4 xzzw() const { return vec4(x, z, z, w); }
	vec4 yzzw() const { return vec4(y, z, z, w); }
	vec4 zzzw() const { return vec4(z, z, z, w); }
	vec4 wzzw() const { return vec4(w, z, z, w); }
	vec4 xwzw() const { return vec4(x, w, z, w); }
	vec4 ywzw() const { return vec4(y, w, z, w); }
	vec4 zwzw() const { return vec4(z, w, z, w); }
	vec4 wwzw() const { return vec4(w, w, z, w); }
	vec4 xxww() const { return vec4(x, x, w, w); }
	vec4 yxww() const { return vec4(y, x, w, w); }
	vec4 zxww() const { return vec4(z, x, w, w); }
	vec4 wxww() const { return vec4(w, x, w, w); }
	vec4 xyww() const { return vec4(x, y, w, w); }
	vec4 yyww() const { return vec4(y, y, w, w); }
	vec4 zyww() const { return vec4(z, y, w, w); }
	vec4 wyww() const { return vec4(w, y, w, w); }
	vec4 xzww() const { return vec4(x, z, w, w); }
	vec4 yzww() const { return vec4(y, z, w, w); }
	vec4 zzww() const { return vec4(z, z, w, w); }
	vec4 wzww() const { return vec4(w, z, w, w); }
	vec4 xwww() const { return vec4(x, w, w, w); }
	vec4 ywww() const { return vec4(y, w, w, w); }
	vec4 zwww() const { return vec4(z, w, w, w); }
	vec4 wwww() const { return vec4(w, w, w, w); }
#endif
};

ivec2::ivec2(const vec2 &p) : x((int)p.x), y((int)p.y) {}
ivec3::ivec3(const vec3 &p) : x((int)p.x), y((int)p.y), z((int)p.z) {}
ivec4::ivec4(const vec4 &p) : x((int)p.x), y((int)p.y), z((int)p.z), w((int)p.w) {}

#if ENABLE_SWIZZLING
vec3 vec2::xxx() const { return vec3(x, x, x); }
vec3 vec2::yxx() const { return vec3(y, x, x); }
vec3 vec2::xyx() const { return vec3(x, y, x); }
vec3 vec2::yyx() const { return vec3(y, y, x); }
vec3 vec2::xxy() const { return vec3(x, x, y); }
vec3 vec2::yxy() const { return vec3(y, x, y); }
vec3 vec2::xyy() const { return vec3(x, y, y); }
vec3 vec2::yyy() const { return vec3(y, y, y); }
vec4 vec2::xxxx() const { return vec4(x, x, x, x); }
vec4 vec2::yxxx() const { return vec4(y, x, x, x); }
vec4 vec2::xyxx() const { return vec4(x, y, x, x); }
vec4 vec2::yyxx() const { return vec4(y, y, x, x); }
vec4 vec2::xxyx() const { return vec4(x, x, y, x); }
vec4 vec2::yxyx() const { return vec4(y, x, y, x); }
vec4 vec2::xyyx() const { return vec4(x, y, y, x); }
vec4 vec2::yyyx() const { return vec4(y, y, y, x); }
vec4 vec2::xxxy() const { return vec4(x, x, x, y); }
vec4 vec2::yxxy() const { return vec4(y, x, x, y); }
vec4 vec2::xyxy() const { return vec4(x, y, x, y); }
vec4 vec2::yyxy() const { return vec4(y, y, x, y); }
vec4 vec2::xxyy() const { return vec4(x, x, y, y); }
vec4 vec2::yxyy() const { return vec4(y, x, y, y); }
vec4 vec2::xyyy() const { return vec4(x, y, y, y); }
vec4 vec2::yyyy() const { return vec4(y, y, y, y); }
vec4 vec3::xxxx() const { return vec4(x, x, x, x); }
vec4 vec3::yxxx() const { return vec4(y, x, x, x); }
vec4 vec3::zxxx() const { return vec4(z, x, x, x); }
vec4 vec3::xyxx() const { return vec4(x, y, x, x); }
vec4 vec3::yyxx() const { return vec4(y, y, x, x); }
vec4 vec3::zyxx() const { return vec4(z, y, x, x); }
vec4 vec3::xzxx() const { return vec4(x, z, x, x); }
vec4 vec3::yzxx() const { return vec4(y, z, x, x); }
vec4 vec3::zzxx() const { return vec4(z, z, x, x); }
vec4 vec3::xxyx() const { return vec4(x, x, y, x); }
vec4 vec3::yxyx() const { return vec4(y, x, y, x); }
vec4 vec3::zxyx() const { return vec4(z, x, y, x); }
vec4 vec3::xyyx() const { return vec4(x, y, y, x); }
vec4 vec3::yyyx() const { return vec4(y, y, y, x); }
vec4 vec3::zyyx() const { return vec4(z, y, y, x); }
vec4 vec3::xzyx() const { return vec4(x, z, y, x); }
vec4 vec3::yzyx() const { return vec4(y, z, y, x); }
vec4 vec3::zzyx() const { return vec4(z, z, y, x); }
vec4 vec3::xxzx() const { return vec4(x, x, z, x); }
vec4 vec3::yxzx() const { return vec4(y, x, z, x); }
vec4 vec3::zxzx() const { return vec4(z, x, z, x); }
vec4 vec3::xyzx() const { return vec4(x, y, z, x); }
vec4 vec3::yyzx() const { return vec4(y, y, z, x); }
vec4 vec3::zyzx() const { return vec4(z, y, z, x); }
vec4 vec3::xzzx() const { return vec4(x, z, z, x); }
vec4 vec3::yzzx() const { return vec4(y, z, z, x); }
vec4 vec3::zzzx() const { return vec4(z, z, z, x); }
vec4 vec3::xxxy() const { return vec4(x, x, x, y); }
vec4 vec3::yxxy() const { return vec4(y, x, x, y); }
vec4 vec3::zxxy() const { return vec4(z, x, x, y); }
vec4 vec3::xyxy() const { return vec4(x, y, x, y); }
vec4 vec3::yyxy() const { return vec4(y, y, x, y); }
vec4 vec3::zyxy() const { return vec4(z, y, x, y); }
vec4 vec3::xzxy() const { return vec4(x, z, x, y); }
vec4 vec3::yzxy() const { return vec4(y, z, x, y); }
vec4 vec3::zzxy() const { return vec4(z, z, x, y); }
vec4 vec3::xxyy() const { return vec4(x, x, y, y); }
vec4 vec3::yxyy() const { return vec4(y, x, y, y); }
vec4 vec3::zxyy() const { return vec4(z, x, y, y); }
vec4 vec3::xyyy() const { return vec4(x, y, y, y); }
vec4 vec3::yyyy() const { return vec4(y, y, y, y); }
vec4 vec3::zyyy() const { return vec4(z, y, y, y); }
vec4 vec3::xzyy() const { return vec4(x, z, y, y); }
vec4 vec3::yzyy() const { return vec4(y, z, y, y); }
vec4 vec3::zzyy() const { return vec4(z, z, y, y); }
vec4 vec3::xxzy() const { return vec4(x, x, z, y); }
vec4 vec3::yxzy() const { return vec4(y, x, z, y); }
vec4 vec3::zxzy() const { return vec4(z, x, z, y); }
vec4 vec3::xyzy() const { return vec4(x, y, z, y); }
vec4 vec3::yyzy() const { return vec4(y, y, z, y); }
vec4 vec3::zyzy() const { return vec4(z, y, z, y); }
vec4 vec3::xzzy() const { return vec4(x, z, z, y); }
vec4 vec3::yzzy() const { return vec4(y, z, z, y); }
vec4 vec3::zzzy() const { return vec4(z, z, z, y); }
vec4 vec3::xxxz() const { return vec4(x, x, x, z); }
vec4 vec3::yxxz() const { return vec4(y, x, x, z); }
vec4 vec3::zxxz() const { return vec4(z, x, x, z); }
vec4 vec3::xyxz() const { return vec4(x, y, x, z); }
vec4 vec3::yyxz() const { return vec4(y, y, x, z); }
vec4 vec3::zyxz() const { return vec4(z, y, x, z); }
vec4 vec3::xzxz() const { return vec4(x, z, x, z); }
vec4 vec3::yzxz() const { return vec4(y, z, x, z); }
vec4 vec3::zzxz() const { return vec4(z, z, x, z); }
vec4 vec3::xxyz() const { return vec4(x, x, y, z); }
vec4 vec3::yxyz() const { return vec4(y, x, y, z); }
vec4 vec3::zxyz() const { return vec4(z, x, y, z); }
vec4 vec3::xyyz() const { return vec4(x, y, y, z); }
vec4 vec3::yyyz() const { return vec4(y, y, y, z); }
vec4 vec3::zyyz() const { return vec4(z, y, y, z); }
vec4 vec3::xzyz() const { return vec4(x, z, y, z); }
vec4 vec3::yzyz() const { return vec4(y, z, y, z); }
vec4 vec3::zzyz() const { return vec4(z, z, y, z); }
vec4 vec3::xxzz() const { return vec4(x, x, z, z); }
vec4 vec3::yxzz() const { return vec4(y, x, z, z); }
vec4 vec3::zxzz() const { return vec4(z, x, z, z); }
vec4 vec3::xyzz() const { return vec4(x, y, z, z); }
vec4 vec3::yyzz() const { return vec4(y, y, z, z); }
vec4 vec3::zyzz() const { return vec4(z, y, z, z); }
vec4 vec3::xzzz() const { return vec4(x, z, z, z); }
vec4 vec3::yzzz() const { return vec4(y, z, z, z); }
vec4 vec3::zzzz() const { return vec4(z, z, z, z); }
#endif

struct mat2 {
	float v[2][2];  // row, col
	mat2() {}
	explicit mat2(float k) {
		v[0][0] = v[1][1] = k, v[1][0] = v[0][1] = 0;
	}
	explicit mat2(vec2 i, vec2 j) {  // by column vectors
		*(vec2*)&v[0] = i, *(vec2*)&v[1] = j;
	}
	explicit mat2(float _00, float _10, float _01, float _11) {  // by columns
		v[0][0] = _00, v[0][1] = _10, v[1][0] = _01, v[1][1] = _11;
	}
	//float* operator[](int i) { return (float*)&v[i][0]; }
	vec2& operator[](int i) { return *(vec2*)&v[i][0]; }
	mat2 operator-() const { mat2 r; for (int i = 0; i < 4; i++) (&r.v[0][0])[i] = -(&v[0][0])[i]; return r; }
	void operator+=(const mat2 &m) { for (int i = 0; i < 4; i++) (&v[0][0])[i] += (&m.v[0][0])[i]; }
	void operator-=(const mat2 &m) { for (int i = 0; i < 4; i++) (&v[0][0])[i] -= (&m.v[0][0])[i]; }
	void operator*=(float k) { for (int i = 0; i < 4; i++) (&v[0][0])[i] *= k; }
	void operator/=(float k) { for (int i = 0; i < 4; i++) (&v[0][0])[i] /= k; }
	mat2 operator+(const mat2 &m) const { mat2 r; for (int i = 0; i < 4; i++) (&r.v[0][0])[i] = (&v[0][0])[i] + (&m.v[0][0])[i]; return r; }
	mat2 operator-(const mat2 &m) const { mat2 r; for (int i = 0; i < 4; i++) (&r.v[0][0])[i] = (&v[0][0])[i] - (&m.v[0][0])[i]; return r; }
	mat2 operator*(float k) const { mat2 r; for (int i = 0; i < 4; i++) (&r.v[0][0])[i] = (&v[0][0])[i] * k; return r; }
	mat2 operator/(float k) const { mat2 r; for (int i = 0; i < 4; i++) (&r.v[0][0])[i] = (&v[0][0])[i] / k; return r; }
	mat2 operator*(const mat2 &a) const {
		mat2 r;
		for (int m = 0; m < 2; m++) for (int n = 0; n < 2; n++) {
			r.v[n][m] = 0;
			for (int i = 0; i < 2; i++) r.v[n][m] += v[i][m] * a.v[n][i];
		}
		return r;
	}
	vec2 operator*(const vec2 &a) const {
		return vec2(
			v[0][0] * a.x + v[1][0] * a.y,
			v[0][1] * a.x + v[1][1] * a.y);
	}
};
mat2 operator*(float k, const mat2 &m) { mat2 r; for (int i = 0; i < 4; i++) (&r.v[0][0])[i] = k * (&m.v[0][0])[i]; return r; }
float determinant(const mat2 &m) { return m.v[0][0] * m.v[1][1] - m.v[0][1] * m.v[1][0]; }
mat2 transpose(const mat2 &m) { mat2 r; for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) r.v[i][j] = m.v[j][i]; return r; }
mat2 inverse(const mat2 &m) {
	float d = 1.0f / (m.v[0][0] * m.v[1][1] - m.v[0][1] * m.v[1][0]);
	return mat2(d * m.v[1][1], -d * m.v[1][0], -d * m.v[0][1], d * m.v[0][0]);
}

struct mat3 {
	float v[3][3];  // row, col
	mat3() {}
	explicit mat3(float k) {  // diagonal matrix
		for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) v[i][j] = k * (i == j);
	}
	explicit mat3(vec3 i, vec3 j, vec3 k) {  // matrix by column vectors
		*(vec3*)&v[0] = i, *(vec3*)&v[1] = j, *(vec3*)&v[2] = k;
	}
	explicit mat3(float _00, float _10, float _20, float _01, float _11, float _21, float _02, float _12, float _22) {  // ordered column-wise
		v[0][0] = _00, v[0][1] = _10, v[0][2] = _20, v[1][0] = _01, v[1][1] = _11, v[1][2] = _21, v[2][0] = _02, v[2][1] = _12, v[2][2] = _22;
	}
	explicit mat3(mat4 m);
	//float* operator[](int i) { return (float*)&v[i][0]; }
	vec3& operator[](int i) { return *(vec3*)&v[i][0]; }
	mat3 operator-() const { mat3 r; for (int i = 0; i < 9; i++) (&r.v[0][0])[i] = -(&v[0][0])[i]; return r; }
	void operator+=(const mat3 &m) { for (int i = 0; i < 9; i++) (&v[0][0])[i] += (&m.v[0][0])[i]; }
	void operator-=(const mat3 &m) { for (int i = 0; i < 9; i++) (&v[0][0])[i] -= (&m.v[0][0])[i]; }
	void operator*=(float k) { for (int i = 0; i < 9; i++) (&v[0][0])[i] *= k; }
	void operator/=(float k) { for (int i = 0; i < 9; i++) (&v[0][0])[i] /= k; }
	mat3 operator+(const mat3 &m) const { mat3 r; for (int i = 0; i < 9; i++) (&r.v[0][0])[i] = (&v[0][0])[i] + (&m.v[0][0])[i]; return r; }
	mat3 operator-(const mat3 &m) const { mat3 r; for (int i = 0; i < 9; i++) (&r.v[0][0])[i] = (&v[0][0])[i] - (&m.v[0][0])[i]; return r; }
	mat3 operator*(float k) const { mat3 r; for (int i = 0; i < 9; i++) (&r.v[0][0])[i] = (&v[0][0])[i] * k; return r; }
	mat3 operator/(float k) const { mat3 r; for (int i = 0; i < 9; i++) (&r.v[0][0])[i] = (&v[0][0])[i] / k; return r; }
	mat3 operator*(const mat3 &a) const {
		mat3 r;
		for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++) {
			r.v[n][m] = 0;
			for (int i = 0; i < 3; i++) r.v[n][m] += v[i][m] * a.v[n][i];
		}
		return r;
	}
	vec3 operator*(const vec3 &a) const {
		return vec3(
			v[0][0] * a.x + v[1][0] * a.y + v[2][0] * a.z,
			v[0][1] * a.x + v[1][1] * a.y + v[2][1] * a.z,
			v[0][2] * a.x + v[1][2] * a.y + v[2][2] * a.z);
	}
};
mat3 operator*(float k, const mat3 &m) { mat3 r; for (int i = 0; i < 9; i++) (&r.v[0][0])[i] = k * (&m.v[0][0])[i]; return r; }
float determinant(const mat3 &m) { return m.v[0][0] * (m.v[1][1] * m.v[2][2] - m.v[1][2] * m.v[2][1]) - m.v[0][1] * (m.v[1][0] * m.v[2][2] - m.v[1][2] * m.v[2][0]) + m.v[0][2] * (m.v[1][0] * m.v[2][1] - m.v[1][1] * m.v[2][0]); }
mat3 transpose(const mat3 &m) { mat3 r; for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) r.v[i][j] = m.v[j][i]; return r; }
mat3 inverse(const mat3 &m) {
	mat3 r;
	r.v[0][0] = m.v[1][1] * m.v[2][2] - m.v[1][2] * m.v[2][1], r.v[0][1] = m.v[0][2] * m.v[2][1] - m.v[0][1] * m.v[2][2], r.v[0][2] = m.v[0][1] * m.v[1][2] - m.v[0][2] * m.v[1][1];
	r.v[1][0] = m.v[1][2] * m.v[2][0] - m.v[1][0] * m.v[2][2], r.v[1][1] = m.v[0][0] * m.v[2][2] - m.v[0][2] * m.v[2][0], r.v[1][2] = m.v[0][2] * m.v[1][0] - m.v[0][0] * m.v[1][2];
	r.v[2][0] = m.v[1][0] * m.v[2][1] - m.v[1][1] * m.v[2][0], r.v[2][1] = m.v[0][1] * m.v[2][0] - m.v[0][0] * m.v[2][1], r.v[2][2] = m.v[0][0] * m.v[1][1] - m.v[0][1] * m.v[1][0];
	float k = 1.0f / (m.v[0][0] * r.v[0][0] + m.v[0][1] * r.v[1][0] + m.v[0][2] * r.v[2][0]);
	return r * k;
}

struct mat4 {
	float v[4][4];  // row, col
	mat4() {}
	explicit mat4(float k) {  // diagonal matrix
		for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) v[i][j] = k * (i == j);
	}
	explicit mat4(vec4 i, vec4 j, vec4 k, vec4 l) {  // matrix by column vectors
		*(vec4*)&v[0] = i, *(vec4*)&v[1] = j, *(vec4*)&v[2] = k, *(vec4*)&v[3] = k;
	}
	explicit mat4(mat3 m) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) v[i][j] = m.v[i][j];
			v[i][3] = v[3][i] = 0.0f;
		}
		v[3][3] = 1.0f;
	}
	explicit mat4(
		float _00, float _10, float _20, float _30,
		float _01, float _11, float _21, float _31,
		float _02, float _12, float _22, float _32,
		float _03, float _13, float _23, float _33) {  // ordered column-wise
		v[0][0] = _00, v[0][1] = _10, v[0][2] = _20, v[0][3] = _30;
		v[1][0] = _01, v[1][1] = _11, v[1][2] = _21, v[1][3] = _31;
		v[2][0] = _02, v[2][1] = _12, v[2][2] = _22, v[2][3] = _32;
		v[3][0] = _03, v[3][1] = _13, v[3][2] = _23, v[3][3] = _33;
	}
	//float* operator[] (int i) { return (float*)&v[i][0]; }
	vec4& operator[] (int i) { return *(vec4*)&v[i][0]; }
	mat4 operator-() const { mat4 r; for (int i = 0; i < 16; i++) (&r.v[0][0])[i] = -(&v[0][0])[i]; return r; }
	void operator+=(const mat4 &m) { for (int i = 0; i < 16; i++) (&v[0][0])[i] += (&m.v[0][0])[i]; }
	void operator-=(const mat4 &m) { for (int i = 0; i < 16; i++) (&v[0][0])[i] -= (&m.v[0][0])[i]; }
	void operator*=(float k) { for (int i = 0; i < 16; i++) (&v[0][0])[i] *= k; }
	void operator/=(float k) { for (int i = 0; i < 16; i++) (&v[0][0])[i] /= k; }
	mat4 operator+(const mat4 &m) const { mat4 r; for (int i = 0; i < 16; i++) (&r.v[0][0])[i] = (&v[0][0])[i] + (&m.v[0][0])[i]; return r; }
	mat4 operator-(const mat4 &m) const { mat4 r; for (int i = 0; i < 16; i++) (&r.v[0][0])[i] = (&v[0][0])[i] - (&m.v[0][0])[i]; return r; }
	mat4 operator*(float k) const { mat4 r; for (int i = 0; i < 16; i++) (&r.v[0][0])[i] = (&v[0][0])[i] * k; return r; }
	mat4 operator/(float k) const { mat4 r; for (int i = 0; i < 16; i++) (&r.v[0][0])[i] = (&v[0][0])[i] / k; return r; }
	mat4 operator*(const mat4 &a) const {
		mat4 r;
		for (int m = 0; m < 4; m++) for (int n = 0; n < 4; n++) {
			r.v[n][m] = 0;
			for (int i = 0; i < 4; i++) r.v[n][m] += v[i][m] * a.v[n][i];
		}
		return r;
	}
	vec4 operator*(const vec4 &a) const {
		return vec4(
			v[0][0] * a.x + v[1][0] * a.y + v[2][0] * a.z + v[3][0] * a.w,
			v[0][1] * a.x + v[1][1] * a.y + v[2][1] * a.z + v[3][1] * a.w,
			v[0][2] * a.x + v[1][2] * a.y + v[2][2] * a.z + v[3][2] * a.w,
			v[0][3] * a.x + v[1][3] * a.y + v[2][3] * a.z + v[3][3] * a.w);
	}
};
mat4 transpose(const mat4 &m) { mat4 r; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) r.v[i][j] = m.v[j][i]; return r; }
mat4 operator*(float k, const mat4 &m) { mat4 r; for (int i = 0; i < 16; i++) (&r.v[0][0])[i] = k * (&m.v[0][0])[i]; return r; }

mat3::mat3(mat4 m) {
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) v[i][j] = m.v[i][j];
}

#undef max
#undef min
#undef abs

#ifdef GLSLMATH_USE_STD_MATH
float sin(float x) { return std::sin(x); }
float cos(float x) { return std::cos(x); }
float tan(float x) { return std::tan(x); }
float asin(float x) { return std::asin(x); }
float acos(float x) { return std::acos(x); }
float atan(float x) { return std::atan(x); }
float sinh(float x) { return std::sinh(x); }
float cosh(float x) { return std::cosh(x); }
float tanh(float x) { return std::tanh(x); }
float asinh(float x) { return std::asinh(x); }
float acosh(float x) { return std::acosh(x); }
float atanh(float x) { return std::atanh(x); }
float exp(float x) { return std::exp(x); }
float log(float x) { return std::log(x); }
float exp2(float x) { return std::exp2(x); }
float log2(float x) { return std::log2(x); }
float sqrt(float x) { return std::sqrt(x); }
float floor(float x) { return std::floor(x); }
float ceil(float x) { return std::ceil(x); }
float round(float x) { return std::round(x); }
float trunc(float x) { return std::trunc(x); }
float pow(float x, float y) { return std::pow(x, y); }
float abs(float x) { return std::abs(x); }
#endif
float atan(float y, float x) { return std::atan2(y, x); }

float radians(float x) { return 0.017453292519943295f * x; }
float degrees(float x) { return 57.29577951308232f * x; }
float inversesqrt(float x) { return 1.0f / sqrt(x); }
float fract(float x) { return x - floor(x); }
float sign(float x) { return x > 0.0f ? 1.0f : x < 0.0f ? -1.0f : 0.0f; }

float length(vec2 v) { return sqrt(v.x * v.x + v.y * v.y); }
float length(vec3 v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }
float length(vec4 v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w); }
vec2 normalize(vec2 v) { return v * (1.0f / length(v)); }
vec3 normalize(vec3 v) { return v * (1.0f / length(v)); }
vec4 normalize(vec4 v) { return v * (1.0f / length(v)); }
float dot(vec2 a, vec2 b) { return a.x * b.x + a.y * b.y; }
float dot(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
float dot(vec4 a, vec4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
vec3 cross(vec3 u, vec3 v) { return vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x); }

float max(float x, float y) { return x > y ? x : y; }
float min(float x, float y) { return x < y ? x : y; }
int max(int x, int y) { return x > y ? x : y; }
int min(int x, int y) { return x < y ? x : y; }
float clamp(float x, float a, float b) { return x < a ? a : x > b ? b : x; }
float saturate(float x) { return clamp(x, 0.0f, 1.0f); }
float mod(float x, float y) { return x - y * floorf(x / y); }
float step(float edge, float x) { return x < edge ? 0.0f : 1.0f; }
float mix(float x, float y, float a) { return x * (1.0f - a) + y * a; }
template <typename vec>
vec mix(vec x, vec y, float a) { return x * (1.0f - a) + y * a; }
template <typename vec>
vec mix(vec x, vec y, vec a) { return x * (1.0f - a) + y * a; }
template <typename T>
T smoothstep(T edge0, T edge1, T x)
{
	T t = clamp((x - edge0) / (edge1 - edge0), T(0.0), T(1.0));
	return t * t * (T(3.0) - T(2.0) * t);
}
vec2 max(vec2 a, vec2 b) { return vec2(max(a.x, b.x), max(a.y, b.y)); }
vec2 min(vec2 a, vec2 b) { return vec2(min(a.x, b.x), min(a.y, b.y)); }
vec3 max(vec3 a, vec3 b) { return vec3(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)); }
vec3 min(vec3 a, vec3 b) { return vec3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)); }
vec4 max(vec4 a, vec4 b) { return vec4(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)); }
vec4 min(vec4 a, vec4 b) { return vec4(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)); }
vec2 clamp(vec2 v, float a, float b) { return vec2(clamp(v.x, a, b), clamp(v.y, a, b)); }
vec3 clamp(vec3 v, float a, float b) { return vec3(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b)); }
vec4 clamp(vec4 v, float a, float b) { return vec4(clamp(v.x, a, b), clamp(v.y, a, b), clamp(v.z, a, b), clamp(v.w, a, b)); }

vec2 saturate(const vec2 &v) { return vec2(saturate(v.x), saturate(v.y)); }
vec2 radians(const vec2 &v) { return vec2(radians(v.x), radians(v.y)); }
vec2 degrees(const vec2 &v) { return vec2(degrees(v.x), degrees(v.y)); }
vec2 sin(const vec2 &v) { return vec2(sin(v.x), sin(v.y)); }
vec2 cos(const vec2 &v) { return vec2(cos(v.x), cos(v.y)); }
vec2 tan(const vec2 &v) { return vec2(tan(v.x), tan(v.y)); }
vec2 asin(const vec2 &v) { return vec2(asin(v.x), asin(v.y)); }
vec2 acos(const vec2 &v) { return vec2(acos(v.x), acos(v.y)); }
vec2 atan(const vec2 &v) { return vec2(atan(v.x), atan(v.y)); }
vec2 sinh(const vec2 &v) { return vec2(sinh(v.x), sinh(v.y)); }
vec2 cosh(const vec2 &v) { return vec2(cosh(v.x), cosh(v.y)); }
vec2 tanh(const vec2 &v) { return vec2(tanh(v.x), tanh(v.y)); }
vec2 asinh(const vec2 &v) { return vec2(asinh(v.x), asinh(v.y)); }
vec2 acosh(const vec2 &v) { return vec2(acosh(v.x), acosh(v.y)); }
vec2 atanh(const vec2 &v) { return vec2(atanh(v.x), atanh(v.y)); }
vec2 exp(const vec2 &v) { return vec2(exp(v.x), exp(v.y)); }
vec2 log(const vec2 &v) { return vec2(log(v.x), log(v.y)); }
vec2 exp2(const vec2 &v) { return vec2(exp2(v.x), exp2(v.y)); }
vec2 log2(const vec2 &v) { return vec2(log2(v.x), log2(v.y)); }
vec2 sqrt(const vec2 &v) { return vec2(sqrt(v.x), sqrt(v.y)); }
vec2 inversesqrt(const vec2 &v) { return vec2(inversesqrt(v.x), inversesqrt(v.y)); }
vec2 abs(const vec2 &v) { return vec2(abs(v.x), abs(v.y)); }
vec2 sign(const vec2 &v) { return vec2(sign(v.x), sign(v.y)); }
vec2 floor(const vec2 &v) { return vec2(floor(v.x), floor(v.y)); }
vec2 ceil(const vec2 &v) { return vec2(ceil(v.x), ceil(v.y)); }
vec2 round(const vec2 &v) { return vec2(round(v.x), round(v.y)); }
vec2 trunc(const vec2 &v) { return vec2(trunc(v.x), trunc(v.y)); }
vec2 fract(const vec2 &v) { return vec2(fract(v.x), fract(v.y)); }

vec3 saturate(const vec3 &v) { return vec3(saturate(v.x), saturate(v.y), saturate(v.z)); }
vec3 radians(const vec3 &v) { return vec3(radians(v.x), radians(v.y), radians(v.z)); }
vec3 degrees(const vec3 &v) { return vec3(degrees(v.x), degrees(v.y), degrees(v.z)); }
vec3 sin(const vec3 &v) { return vec3(sin(v.x), sin(v.y), sin(v.z)); }
vec3 cos(const vec3 &v) { return vec3(cos(v.x), cos(v.y), cos(v.z)); }
vec3 tan(const vec3 &v) { return vec3(tan(v.x), tan(v.y), tan(v.z)); }
vec3 asin(const vec3 &v) { return vec3(asin(v.x), asin(v.y), asin(v.z)); }
vec3 acos(const vec3 &v) { return vec3(acos(v.x), acos(v.y), acos(v.z)); }
vec3 atan(const vec3 &v) { return vec3(atan(v.x), atan(v.y), atan(v.z)); }
vec3 sinh(const vec3 &v) { return vec3(sinh(v.x), sinh(v.y), sinh(v.z)); }
vec3 cosh(const vec3 &v) { return vec3(cosh(v.x), cosh(v.y), cosh(v.z)); }
vec3 tanh(const vec3 &v) { return vec3(tanh(v.x), tanh(v.y), tanh(v.z)); }
vec3 asinh(const vec3 &v) { return vec3(asinh(v.x), asinh(v.y), asinh(v.z)); }
vec3 acosh(const vec3 &v) { return vec3(acosh(v.x), acosh(v.y), acosh(v.z)); }
vec3 atanh(const vec3 &v) { return vec3(atanh(v.x), atanh(v.y), atanh(v.z)); }
vec3 exp(const vec3 &v) { return vec3(exp(v.x), exp(v.y), exp(v.z)); }
vec3 log(const vec3 &v) { return vec3(log(v.x), log(v.y), log(v.z)); }
vec3 exp2(const vec3 &v) { return vec3(exp2(v.x), exp2(v.y), exp2(v.z)); }
vec3 log2(const vec3 &v) { return vec3(log2(v.x), log2(v.y), log2(v.z)); }
vec3 sqrt(const vec3 &v) { return vec3(sqrt(v.x), sqrt(v.y), sqrt(v.z)); }
vec3 inversesqrt(const vec3 &v) { return vec3(inversesqrt(v.x), inversesqrt(v.y), inversesqrt(v.z)); }
vec3 abs(const vec3 &v) { return vec3(abs(v.x), abs(v.y), abs(v.z)); }
vec3 sign(const vec3 &v) { return vec3(sign(v.x), sign(v.y), sign(v.z)); }
vec3 floor(const vec3 &v) { return vec3(floor(v.x), floor(v.y), floor(v.z)); }
vec3 ceil(const vec3 &v) { return vec3(ceil(v.x), ceil(v.y), ceil(v.z)); }
vec3 round(const vec3 &v) { return vec3(round(v.x), round(v.y), round(v.z)); }
vec3 trunc(const vec3 &v) { return vec3(trunc(v.x), trunc(v.y), trunc(v.z)); }
vec3 fract(const vec3 &v) { return vec3(fract(v.x), fract(v.y), fract(v.z)); }

vec4 saturate(const vec4 &v) { return vec4(saturate(v.x), saturate(v.y), saturate(v.z), saturate(v.w)); }
vec4 radians(const vec4 &v) { return vec4(radians(v.x), radians(v.y), radians(v.z), radians(v.w)); }
vec4 degrees(const vec4 &v) { return vec4(degrees(v.x), degrees(v.y), degrees(v.z), degrees(v.w)); }
vec4 sin(const vec4 &v) { return vec4(sin(v.x), sin(v.y), sin(v.z), sin(v.w)); }
vec4 cos(const vec4 &v) { return vec4(cos(v.x), cos(v.y), cos(v.z), cos(v.w)); }
vec4 tan(const vec4 &v) { return vec4(tan(v.x), tan(v.y), tan(v.z), tan(v.w)); }
vec4 asin(const vec4 &v) { return vec4(asin(v.x), asin(v.y), asin(v.z), asin(v.w)); }
vec4 acos(const vec4 &v) { return vec4(acos(v.x), acos(v.y), acos(v.z), acos(v.w)); }
vec4 atan(const vec4 &v) { return vec4(atan(v.x), atan(v.y), atan(v.z), atan(v.w)); }
vec4 sinh(const vec4 &v) { return vec4(sinh(v.x), sinh(v.y), sinh(v.z), sinh(v.w)); }
vec4 cosh(const vec4 &v) { return vec4(cosh(v.x), cosh(v.y), cosh(v.z), cosh(v.w)); }
vec4 tanh(const vec4 &v) { return vec4(tanh(v.x), tanh(v.y), tanh(v.z), tanh(v.w)); }
vec4 asinh(const vec4 &v) { return vec4(asinh(v.x), asinh(v.y), asinh(v.z), asinh(v.w)); }
vec4 acosh(const vec4 &v) { return vec4(acosh(v.x), acosh(v.y), acosh(v.z), acosh(v.w)); }
vec4 atanh(const vec4 &v) { return vec4(atanh(v.x), atanh(v.y), atanh(v.z), atanh(v.w)); }
vec4 exp(const vec4 &v) { return vec4(exp(v.x), exp(v.y), exp(v.z), exp(v.w)); }
vec4 log(const vec4 &v) { return vec4(log(v.x), log(v.y), log(v.z), log(v.w)); }
vec4 exp2(const vec4 &v) { return vec4(exp2(v.x), exp2(v.y), exp2(v.z), exp2(v.w)); }
vec4 log2(const vec4 &v) { return vec4(log2(v.x), log2(v.y), log2(v.z), log2(v.w)); }
vec4 sqrt(const vec4 &v) { return vec4(sqrt(v.x), sqrt(v.y), sqrt(v.z), sqrt(v.w)); }
vec4 inversesqrt(const vec4 &v) { return vec4(inversesqrt(v.x), inversesqrt(v.y), inversesqrt(v.z), inversesqrt(v.w)); }
vec4 abs(const vec4 &v) { return vec4(abs(v.x), abs(v.y), abs(v.z), abs(v.w)); }
vec4 sign(const vec4 &v) { return vec4(sign(v.x), sign(v.y), sign(v.z), sign(v.w)); }
vec4 floor(const vec4 &v) { return vec4(floor(v.x), floor(v.y), floor(v.z), floor(v.w)); }
vec4 ceil(const vec4 &v) { return vec4(ceil(v.x), ceil(v.y), ceil(v.z), ceil(v.w)); }
vec4 round(const vec4 &v) { return vec4(round(v.x), round(v.y), round(v.z), round(v.w)); }
vec4 trunc(const vec4 &v) { return vec4(trunc(v.x), trunc(v.y), trunc(v.z), trunc(v.w)); }
vec4 fract(const vec4 &v) { return vec4(fract(v.x), fract(v.y), fract(v.z), fract(v.w)); }

vec2 pow(vec2 v, vec2 k) { return vec2(pow(v.x, k.x), pow(v.y, k.y)); }
vec3 pow(vec3 v, vec3 k) { return vec3(pow(v.x, k.x), pow(v.y, k.y), pow(v.z, k.z)); }
vec4 pow(vec4 v, vec4 k) { return vec4(pow(v.x, k.x), pow(v.y, k.y), pow(v.z, k.z), pow(v.w, k.w)); }

template<typename vec>
vec reflect(vec I, vec N) {
	return I - 2.0f * dot(N, I) * N;
}
template<typename vec>
vec refract(vec I, vec N, float eta) {
	float k = 1.0f - eta * eta * (1.0f - dot(N, I) * dot(N, I));
	if (k < 0.0f)
		return vec(0.0f);
	else
		return eta * I - (eta * dot(N, I) + sqrt(k)) * N;
}


// texture

#define STB_IMAGE_IMPLEMENTATION
#include <libraries/stb_image.h>

class sampler2D {
	int w, h;
	vec4 *data;
public:
	sampler2D(const char filename[]) {
		uint8_t *pixels = stbi_load(filename, &w, &h, nullptr, 4);
		if (pixels == nullptr) {
			w = h = 0, data = nullptr;
			return;
		}
		data = new vec4[w*h];
		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				uint8_t *pixel = &pixels[4 * ((h - 1 - j)*w + i)];
				data[j*w + i] = vec4(pixel[0], pixel[1], pixel[2], pixel[3]) / 255.0f;
			}
		}
		delete pixels;
	}
	~sampler2D() {
		delete data; data = nullptr;
		w = h = 0;
	}
	friend vec4 texelFetch(const sampler2D &sampler, ivec2 pos, int plane);
	friend vec4 texture(const sampler2D &sampler, vec2 uv);
};
vec4 texelFetch(const sampler2D &sampler, ivec2 pos, int plane) {
	ivec2 res = ivec2(sampler.w, sampler.h);
	pos = (pos % res + res) % res;
	return sampler.data[pos.y*sampler.w + pos.x];
}
vec4 texture(const sampler2D &sampler, vec2 uv) {
	uv = uv * vec2(sampler.w, sampler.h) - 0.5;
	vec2 p0 = floor(uv), pf = uv - p0;
	vec2 p1 = p0 + 1.0;
	return mix(  // bilinear interpolation
		mix(texelFetch(sampler, ivec2(p0) + ivec2(0, 0), 0),
			texelFetch(sampler, ivec2(p0) + ivec2(0, 1), 0), pf.x),
		mix(texelFetch(sampler, ivec2(p0) + ivec2(1, 0), 0),
			texelFetch(sampler, ivec2(p0) + ivec2(1, 1), 0), pf.x),
		pf.y);
}
