#ifndef MMATH_HEADER_FILE
#define MMATH_HEADER_FILE

/* MMath.h -- MMath library file
 * 
 * Copyright (C) 2017 Zachary Wells
 *
 * This software may be modified and distributed under the terms
 * of the MIT license.  See the LICENSE file for details.
 */

//Open Ended math library v0.1

/*
Library cheat sheet:

Parameters:
	a - left side operand
	b - right side operand

	r - angle in radians
	pitch \
	yaw    |- angles in radians
	roll  /

	s - scale
	t - translation

	aspect - aspect ratio (width / height)
	fovY - field of view of the Y axis ( fovX calculated with 2 * atan(aspect * tan(fovY / 2)) )
	near - near Z clipping plane
	far - far Z clipping plane

	axis - unit vector
	
WARNINGS (If you're having a problem, read this first):
	-Unit vector parameters are NOT checked, make sure they are normalized!
		-vecXNormalize(vecX a)

	-Quaternions are NOT guaranteed to be normalized after functions (most still are though).
		-quatNormalize(quat a)

	-All angles are ALWAYS in radians
		-radians(scalar degrees)
		-degrees(scalar radians)
*/

#if defined(__cplusplus)
extern "C" {
#endif

	#include <math.h>

	#define MMATH_INLINE inline

	//Types
	#pragma region Types
	#if defined(MMATH_DOUBLE)
	typedef double scalar;
	#else
	typedef float scalar;
	#endif
	typedef struct vec2_t {
		union {
			scalar data[2];
			struct { scalar x, y; };
			//struct { scalar g, a; };
		};
	} vec2;
	typedef struct vec3_t {
		union {
			scalar data[3];
			struct { scalar x, y, z; };
			struct { scalar r, g, b; };
		};
	} vec3;
	typedef struct vec4_t {
		union {
			scalar data[4];
			struct { scalar x, y, z, w; };
			struct { scalar r, g, b, a; };
		};
	} vec4;

	typedef struct quat_t {
		union {
			scalar data[4];
			vec4 vec;
			vec3 axis;
			struct { float x, y, z, w; };
		};
	} quat;

	typedef struct mat2_t {
		union {
			scalar data[2 * 2];
			vec2 row[2];
			struct {
				scalar x0, y0;
				scalar x1, y1;
			};
			struct {
				vec2 r0;
				vec2 r1;
			};
		};
	} mat2;
	typedef struct mat3_t {
		union {
			scalar data[3 * 3];
			vec3  row[3];
			struct {
				scalar x0, y0, z0;
				scalar x1, y1, z1;
				scalar x2, y2, z2;
			};
			struct {
				vec3 r0;
				vec3 r1;
				vec3 r2;
			};
		};
	} mat3;
	typedef struct mat4_t {
		union {
			scalar data[4 * 4];
			vec4 row[4];
			struct {
				scalar x0, y0, z0, w0;
				scalar x1, y1, z1, w1;
				scalar x2, y2, z2, w2;
				scalar x3, y3, z3, w3;
			};
			struct {
				vec4 r0;
				vec4 r1;
				vec4 r2;
				vec4 r3;
			};
		};
	} mat4;
	#pragma endregion Types

	#define mm_sqrt(var) ((scalar)sqrt(var))
	#define mm_cos(var)  ((scalar)cos(var))
	#define mm_sin(var)  ((scalar)sin(var))
	#define mm_tan(var)  ((scalar)tan(var))

	//constants
	#define mm_dpi ((scalar)6.283185307179586)
	#define mm_pi  ((scalar)3.141592653589793)
	#define mm_hpi ((scalar)1.570796326794896)

	MMATH_INLINE scalar radians(scalar degrees) {
		return degrees * (scalar)0.017453292519943295; //PI / 180
	}
	MMATH_INLINE scalar degrees(scalar radians) {
		return radians * (scalar)57.29577951308232; //180 / PI
	}

	//Vector Math
	#pragma region Vec_TempDef
	#define VEC_FOR(integer) for (int i = 0; i < integer; i++)
	#define GENFUNC_VECADDSCALAR(integer) \
	MMATH_INLINE vec##integer vec##integer##AddScalar(vec##integer a, scalar b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] + b; \
		} \
		return ret; \
	}
	#define GENFUNC_VECSUBSCALAR(integer) \
	MMATH_INLINE vec##integer vec##integer##SubScalar(vec##integer a, scalar b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] - b; \
		} \
		return ret; \
	}
	#define GENFUNC_VECMULSCALAR(integer) \
	MMATH_INLINE vec##integer vec##integer##MulScalar(vec##integer a, scalar b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] * b; \
		} \
		return ret; \
	}
	#define GENFUNC_VECDIVSCALAR(integer) \
	MMATH_INLINE vec##integer vec##integer##DivScalar(vec##integer a, scalar b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] / b; \
		} \
		return ret; \
	}
	#define GENFUNC_VECADD(integer) \
	MMATH_INLINE vec##integer vec##integer##Add(vec##integer a, vec##integer b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] + b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_VECSUB(integer) \
	MMATH_INLINE vec##integer vec##integer##Sub(vec##integer a, vec##integer b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] - b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_VECMUL(integer) \
	MMATH_INLINE vec##integer vec##integer##Mul(vec##integer a, vec##integer b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] * b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_VECDIV(integer) \
	MMATH_INLINE vec##integer vec##integer##Div(vec##integer a, vec##integer b) { \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] / b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_VECDOT(integer) \
	MMATH_INLINE scalar vec##integer##Dot(vec##integer a, vec##integer b) { \
		scalar ret = 0; \
		VEC_FOR(integer) { \
			ret += a.data[i] * b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_VECLEN(integer) \
	MMATH_INLINE scalar vec##integer##Length(vec##integer a) { \
		scalar sum = 0; \
		VEC_FOR(integer) { \
			sum += a.data[i] * a.data[i]; \
		} \
		return mm_sqrt(sum); \
	}
	#define GENFUNC_VECNORM(integer) \
	MMATH_INLINE vec##integer vec##integer##Normalize(vec##integer a) { \
		scalar len = vec##integer##Length(a); \
		vec##integer ret; \
		VEC_FOR(integer) { \
			ret.data[i] = a.data[i] / len; \
		} \
		return ret; \
	}
	#define GENFUNC_VECSTANDARD(integer) \
		GENFUNC_VECADDSCALAR(integer) \
		GENFUNC_VECSUBSCALAR(integer) \
		GENFUNC_VECMULSCALAR(integer) \
		GENFUNC_VECDIVSCALAR(integer) \
		GENFUNC_VECADD(integer); \
		GENFUNC_VECSUB(integer); \
		GENFUNC_VECMUL(integer); \
		GENFUNC_VECDIV(integer); \
		GENFUNC_VECDOT(integer); \
		GENFUNC_VECLEN(integer); \
		GENFUNC_VECNORM(integer);
	#pragma endregion Vec_TempDef

	#pragma region Vec2_Functions
	GENFUNC_VECSTANDARD(2);
	MMATH_INLINE vec3 vec2ToVec3(vec2 a, scalar z) {
		vec3 ret = {a.x, a.y, z};
		return ret;
	}
	MMATH_INLINE vec4 vec2ToVec4(vec2 a, scalar z, scalar w) {
		vec4 ret = {a.x, a.y, z, w};
		return ret;
	}
	#pragma endregion Vec2_Functions

	#pragma region Vec3_Functions
	GENFUNC_VECSTANDARD(3);
	MMATH_INLINE vec2 vec3ToVec2(vec3 a) {
		vec2 ret = {a.x, a.y};
		return ret;
	}
	MMATH_INLINE vec4 vec3ToVec4(vec3 a, scalar w) {
		vec4 ret = {a.x, a.y, a.z, w};
		return ret;
	}
	MMATH_INLINE quat vec3ToQuat(vec3 a) {
		quat ret;

		return ret;
	}
	MMATH_INLINE vec3 vec3Cross(vec3 a, vec3 b) {
		vec3 ret;
		ret.x = a.y * b.z - a.z * b.y;
		ret.y = a.z * b.x - a.x * b.z;
		ret.z = a.x * b.y - a.y * b.x;
		return ret;
	}
	#pragma endregion Vec3_Functions

	#pragma region Vec4_Functions
	GENFUNC_VECSTANDARD(4);
	MMATH_INLINE vec2 vec4ToVec2(vec4 a) {
		vec2 ret = {a.x, a.y};
		return ret;
	}
	MMATH_INLINE vec3 vec4ToVec3(vec4 a) {
		vec3 ret = {a.x, a.y, a.z};
		return ret;
	}
	MMATH_INLINE vec3 vec4DivW(vec4 a) {
		vec3 ret = {a.x / a.w, a.y / a.w, a.z / a.w};
		return ret;
	}
	#pragma endregion Vec4_Functions

	//Quaternion Math
	#pragma region Quaternion_Functions
	MMATH_INLINE scalar quatLength(quat a) {
		return vec4Length(a.vec);
	}
	MMATH_INLINE quat quatNormalize(quat a) {
		quat ret;
		ret.vec = vec4Normalize(a.vec);
		return ret;
	}
	MMATH_INLINE quat quatInverse(quat a) {
		quat ret = {
			.x = (-a.x) / (a.x * a.x),
			.y = (-a.y) / (a.y * a.y),
			.z = (-a.z) / (a.z * a.z),
			.w = (a.w)  / (a.w * a.w)
		};
		return ret;
	}
	MMATH_INLINE quat quatEurler(scalar pitch, scalar yaw, scalar roll) {
		quat ret;

		scalar cy = mm_cos(yaw   * 0.5f);
		scalar sy = mm_sin(yaw   * 0.5f);
		scalar cr = mm_cos(roll  * 0.5f);
		scalar sr = mm_sin(roll  * 0.5f);
		scalar cp = mm_cos(pitch * 0.5f);
		scalar sp = mm_sin(pitch * 0.5f);

		ret.x = cy * sr * cp - sy * cr * sp;
		ret.y = cy * cr * sp + sy * sr * cp;
		ret.z = sy * cr * cp - cy * sr * sp;
		ret.w = cy * cr * cp + sy * sr * sp;

		return ret;
	}
	MMATH_INLINE quat quatAxisAngle(vec3 axis, scalar r) {
		scalar a2 = r * (scalar)0.5;
		quat ret;
		ret.axis = vec3MulScalar(axis, mm_sin(a2));
		ret.w    = mm_cos(a2);
		return ret;
	}
	MMATH_INLINE quat quatMul(quat a, quat b) {
		quat ret;
		ret.w = a.w * b.w - vec3Dot(a.axis, b.axis);

		vec3 BwAv = vec3MulScalar(a.axis, b.w);
		vec3 AwBv = vec3MulScalar(b.axis, a.w);
		vec3 abv  = vec3Add(BwAv, AwBv);
		vec3 AxB  = vec3Cross(a.axis, b.axis);

		ret.axis = vec3Add(abv, AxB);

		return ret;
	}
	MMATH_INLINE vec3 quatMulVec3(quat a, vec3 b) {
		vec3 cross1  = vec3Cross(a.axis, b);
		vec3 t       = vec3MulScalar(cross1, (scalar)2.0);
		vec3 tw      = vec3MulScalar(t, a.w);
		vec3 cross2  = vec3Cross(a.axis, t);
		vec3 twcross = vec3Add(tw, cross2);

		vec3 ret = vec3Add(b, twcross);
		return ret;
	}
	MMATH_INLINE mat3 quatToMat3(quat a) {
		scalar x2 = 2 * a.x * a.x,
			   y2 = 2 * a.y * a.y,
			   z2 = 2 * a.z * a.z,
			   xy = 2 * a.x * a.y,
			   xz = 2 * a.x * a.z,
			   yz = 2 * a.y * a.z,
			   xw = 2 * a.x * a.w,
			   yw = 2 * a.y * a.w,
			   zw = 2 * a.z * a.w;
		mat3 ret = {
			1-y2-z2, xy+zw,   xz-yw,
			xy-zw,   1-x2-z2, yz+xw,
			xz+yz,   yz-xw,   1-x2-y2
		};
		return ret;
	}
	#pragma endregion Quaternion_Functions

	//Matrix Math
	#pragma region Mat_TempDef
	#define MAT_FOR_FLAT(integer) for (int i = 0; i < integer * integer; i++)
	#define MAT_FOR(integer) for (int x = 0; x < integer; x++) for (int y = 0; y < integer; y++)
	#define GENFUNC_MATTRPOSE(integer) \
	MMATH_INLINE mat##integer mat##integer##Transpose(mat##integer a) { \
		mat##integer ret; \
		for (int x = 0; x < integer; x++) { \
			for (int y = 0; y < integer; y++) { \
				ret.row[x].data[y] = a.row[y].data[x]; \
			} \
		} \
		return ret; \
	}
	#define GENFUNC_MATIDENT(integer) \
	MMATH_INLINE mat##integer mat##integer##Identity(void) { \
		mat##integer ret = {0}; \
		for (int i = 0; i < integer * integer; i += integer + 1) { \
			ret.data[i] = (scalar)1.0; \
		} \
		return ret; \
	}
	#define GENFUNC_MATDIAG(integer) \
	MMATH_INLINE mat##integer mat##integer##Diagonal(scalar f) { \
		mat##integer ret = {0}; \
		for (int i = 0; i < integer * integer; i += integer + 1) { \
			ret.data[i] = f; \
		} \
		return ret; \
	}
	#define GENFUNC_MATADD(integer) \
	MMATH_INLINE mat##integer mat##integer##Add(mat##integer a, mat##integer b) { \
		mat##integer ret; \
		MAT_FOR_FLAT(integer) { \
			ret.data[i] = a.data[i] + b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_MATSUB(integer) \
	MMATH_INLINE mat##integer mat##integer##Sub(mat##integer a, mat##integer b) { \
		mat##integer ret; \
		MAT_FOR_FLAT(integer) { \
			ret.data[i] = a.data[i] - b.data[i]; \
		} \
		return ret; \
	}
	#define GENFUNC_MATMUL(integer) \
	MMATH_INLINE mat##integer mat##integer##Mul(mat##integer a, mat##integer b) { \
		mat##integer ret = {0}; \
		MAT_FOR(integer) { \
			VEC_FOR(integer) { \
				ret.row[x].data[y] += a.row[x].data[i] * b.row[i].data[y]; \
			} \
		} \
		return ret; \
	}
	#define GENFUNC_MATMULSCALAR(integer) \
	MMATH_INLINE mat##integer mat##integer##MulScalar(mat##integer a, scalar b) { \
		mat##integer ret; \
		MAT_FOR_FLAT(integer) { \
			ret.data[i] = a.data[i] * b; \
		} \
		return ret; \
	}
	#define GENFUNC_MATMULVEC(integer) \
	MMATH_INLINE vec##integer mat##integer##MulVec##integer(mat##integer a, vec##integer b) { \
		vec##integer ret = {0}; \
		VEC_FOR(integer) { \
			for (int c = 0; c < integer; c++) { \
				ret.data[i] += a.row[c].data[i] * b.data[c]; \
			} \
		} \
		return ret; \
	}
	#define GENFUNC_MATSTANDARD(integer) \
		GENFUNC_MATTRPOSE(integer); \
		GENFUNC_MATIDENT(integer); \
		GENFUNC_MATDIAG(integer); \
		GENFUNC_MATADD(integer); \
		GENFUNC_MATSUB(integer); \
		GENFUNC_MATMUL(integer); \
		GENFUNC_MATMULSCALAR(integer); \
		GENFUNC_MATMULVEC(integer);
	#pragma endregion Mat_TempDef
	
	#pragma region Mat2_Functions
	GENFUNC_MATSTANDARD(2);
	MMATH_INLINE mat3 mat2ToMat3(mat2 a) {
		mat3 ret = {
			a.x0, a.y0, 0,
			a.x1, a.y1, 0,
			0,    0,    1
		};
		return ret;
	}
	MMATH_INLINE mat4 mat2ToMat4(mat2 a) {
		mat4 ret = {
			a.x0, a.y0, 0, 0,
			a.x1, a.y1, 0, 0,
			0,    0,    1, 0,
			0,    0,    0, 1
		};
		return ret;
	}
	#pragma endregion Mat2_Functions

	#pragma region Mat3_Functions
	GENFUNC_MATSTANDARD(3);
	MMATH_INLINE mat3 mat3RotateX(scalar r) {
		scalar c = mm_cos(r);
		scalar s = mm_sin(r);
		mat3 ret = {
			1, 0, 0,
			0, c, s,
			0,-s, c
		};
		return ret;
	}
	MMATH_INLINE mat3 mat3RotateY(scalar r) {
		scalar c = mm_cos(r);
		scalar s = mm_sin(r);
		mat3 ret = {
			c, 0,-s,
			0, 1, 0,
		    s, 0, c
		};
		return ret;
	}
	MMATH_INLINE mat3 mat3RotateZ(scalar r) {
		scalar c = mm_cos(r);
		scalar s = mm_sin(r);
		mat3 ret = {
			c, s, 0,
		   -s, c, 0,
			0, 0, 1
		};
		return ret;
	}
	MMATH_INLINE mat2 mat3ToMat2(mat3 a) {
		mat2 ret = {
			a.x0, a.y0,
			a.x1, a.y1
		};
		return ret;
	}
	MMATH_INLINE mat4 mat3ToMat4(mat3 a) {
		mat4 ret = {
			a.x0, a.y0, a.z0, 0,
			a.x1, a.y1, a.z1, 0,
			a.x2, a.y2, a.z2, 0,
			0,    0,    0,    1
		};
		return ret;
	}
	#pragma endregion Mat3_Functions

	#pragma region Mat4_Functions
	GENFUNC_MATSTANDARD(4);
	MMATH_INLINE mat4 mat4Scale(vec3 s) {
		mat4 ret = {
			s.x, 0,   0,   0,
			0,   s.y, 0,   0,
			0,   0,   s.z, 0,
			0,   0,   0,   1
		};
		return ret;
	}
	MMATH_INLINE mat4 mat4Translate(vec3 t) {
		mat4 ret = {
			1,   0,   0,   0,
			0,   1,   0,   0,
			0,   0,   1,   0,
			t.x, t.y, t.z, 1
		};
		return ret;
	}
	MMATH_INLINE mat4 mat4Perspective(scalar aspect, scalar fovY, scalar near, scalar far) {
		scalar f   = (scalar)1.0 / mm_tan(fovY * 0.5);
		scalar nf  = (scalar)1.0 / (near - far);
		mat4 ret = {
			f/aspect, 0,  0,             0,
			0,        f,  0,             0,
			0,        0, (far+near)*nf, -1,
			0,        0, (2*far*near)*nf,0
		};
		return ret;
	}
	MMATH_INLINE mat2 mat4ToMat2(mat4 a) {
		mat2 ret = {
			a.x0, a.y0,
			a.x1, a.y1
		};
		return ret;
	}
	MMATH_INLINE mat3 mat4ToMat3(mat4 a) {
		mat3 ret = {
			a.x0, a.y0, a.z0,
			a.x1, a.y1, a.z1,
			a.x2, a.y2, a.z2
		};
		return ret;
	}
	#pragma endregion Mat4_Functions

#if defined(__cplusplus)
}
#endif

#endif //MMATH_HEADER_FILE
