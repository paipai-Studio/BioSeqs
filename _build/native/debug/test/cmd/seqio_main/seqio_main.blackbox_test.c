#ifdef __cplusplus
extern "C" {
#endif

#include "moonbit.h"
#include "moonbit_simd.h"

#ifdef _MSC_VER
#define _Noreturn __declspec(noreturn)
#endif

#if defined(__clang__)
#pragma clang diagnostic ignored "-Wshift-op-parentheses"
#pragma clang diagnostic ignored "-Wtautological-compare"
#endif

MOONBIT_EXPORT _Noreturn void moonbit_panic(void);
MOONBIT_EXPORT void *moonbit_malloc_array(enum moonbit_block_kind kind,
                                          int elem_size_shift, int32_t len);
int memcmp(const void *s1, const void *s2, size_t n);
MOONBIT_EXPORT int moonbit_val_array_equal_sized(const void *lhs,
                                                 const void *rhs,
                                                 int32_t elem_size);
MOONBIT_EXPORT moonbit_string_t moonbit_add_string(moonbit_string_t s1,
                                                   moonbit_string_t s2);
MOONBIT_EXPORT void moonbit_unsafe_bytes_blit(moonbit_bytes_t dst,
                                              int32_t dst_start,
                                              moonbit_bytes_t src,
                                              int32_t src_offset, int32_t len);
MOONBIT_EXPORT moonbit_string_t moonbit_unsafe_bytes_sub_string(
    moonbit_bytes_t bytes, int32_t start, int32_t len);
MOONBIT_EXPORT int32_t moonbit_unsafe_val_array_blit(void *dst,
                                                     int32_t dst_offset,
                                                     void *src,
                                                     int32_t src_offset,
                                                     int32_t len,
                                                     int32_t elem_size);
MOONBIT_EXPORT int32_t moonbit_unsafe_ref_array_blit(void *dst,
                                                     int32_t dst_offset,
                                                     void *src,
                                                     int32_t src_offset,
                                                     int32_t len);
MOONBIT_EXPORT void moonbit_println(moonbit_string_t str);
MOONBIT_EXPORT moonbit_bytes_t *moonbit_get_cli_args(void);
MOONBIT_EXPORT void moonbit_runtime_init(int argc, char **argv);
MOONBIT_EXPORT void moonbit_drop_object(void *);
MOONBIT_EXPORT int32_t moonbit_utf16_len_from_utf8(moonbit_bytes_t src,
                                                   int32_t src_offset,
                                                   int32_t src_length);
MOONBIT_EXPORT int32_t moonbit_utf8_decode_into_utf16(
    moonbit_bytes_t src, int32_t src_offset, int32_t src_length,
    moonbit_string_t dst, int32_t dst_offset);
MOONBIT_EXPORT int32_t moonbit_utf8_decode_lossy_into_utf16(
    moonbit_bytes_t src, int32_t src_offset, int32_t src_length,
    moonbit_string_t dst, int32_t dst_offset);
MOONBIT_EXPORT int32_t moonbit_utf8_len_from_utf16(moonbit_string_t src,
                                                   int32_t src_offset,
                                                   int32_t src_length);
MOONBIT_EXPORT int32_t moonbit_utf8_encode_from_utf16(
    moonbit_string_t src, int32_t src_offset, int32_t src_length,
    moonbit_bytes_t dst, int32_t dst_offset);

#if !defined(_WIN64) && !defined(_WIN32)
void *malloc(size_t size);
void free(void *ptr);
#define libc_malloc malloc
#define libc_free free
#endif

// several important runtime functions are inlined
static void *moonbit_malloc_inlined(size_t size) {
  struct moonbit_object *ptr = (struct moonbit_object *)libc_malloc(
      sizeof(struct moonbit_object) + size);
  Moonbit_init_dynamic_rc(ptr, moonbit_BLOCK_KIND_REGULAR);
  return ptr + 1;
}

#define moonbit_malloc(obj) moonbit_malloc_inlined(obj)
#define moonbit_free(obj) libc_free(Moonbit_object_header(obj))

#define MOONBIT_RC_COUNT_UNIT ((int32_t)(1u << MOONBIT_RC_COUNT_SHIFT))
#define raw_rc_is_dynamic(rc) ((int32_t)(rc) >= MOONBIT_RC_COUNT_UNIT)
#define raw_rc_is_shared(rc) ((int32_t)(rc) >= (MOONBIT_RC_COUNT_UNIT * 2))

extern const uint32_t *moonbit_layout_table;

static void moonbit_incref_inlined(void *ptr) {
  struct moonbit_object *header = Moonbit_object_header(ptr);
  int32_t const rc = header->rc;
  if (raw_rc_is_dynamic(rc)) {
    Moonbit_increase_rc_count(header);
  }
}

#define moonbit_incref moonbit_incref_inlined

static void moonbit_decref_inlined(void *ptr) {
  struct moonbit_object *header = Moonbit_object_header(ptr);
  int32_t const rc = header->rc;
  if (raw_rc_is_shared(rc)) {
    header->rc = rc - MOONBIT_RC_COUNT_UNIT;
  } else if (raw_rc_is_dynamic(rc)) {
    moonbit_drop_object(ptr);
  }
}

#define moonbit_decref moonbit_decref_inlined

#define moonbit_unsafe_make_string moonbit_make_string

#if defined(MOONBIT_V128_NEON)
#define Moonbit_v128_make(lo, hi)                                             \
  vreinterpretq_u8_u64(                                                       \
      vcombine_u64(vcreate_u64((uint64_t)(lo)), vcreate_u64((uint64_t)(hi))))
#define Moonbit_v128_lo(v) vgetq_lane_u64(vreinterpretq_u64_u8(v), 0)
#define Moonbit_v128_hi(v) vgetq_lane_u64(vreinterpretq_u64_u8(v), 1)
#define Moonbit_v128_load_storage(p) vld1q_u8((const uint8_t *)(p))
#define Moonbit_v128_store_storage(p, v) vst1q_u8((uint8_t *)(p), (v))
#elif defined(MOONBIT_V128_SSE2)
#define Moonbit_v128_make(lo, hi) _mm_set_epi64x((int64_t)(hi), (int64_t)(lo))
#define Moonbit_v128_lo(v) ((uint64_t)_mm_cvtsi128_si64(v))
#define Moonbit_v128_hi(v) ((uint64_t)_mm_cvtsi128_si64(_mm_srli_si128((v), 8)))
#define Moonbit_v128_load_storage(p) _mm_loadu_si128((const __m128i *)(p))
#define Moonbit_v128_store_storage(p, v) _mm_storeu_si128((__m128i *)(p), (v))
#else
#define Moonbit_v128_make(lo, hi) ((moonbit_v128_t){(lo), (hi)})
#define Moonbit_v128_lo(v) ((v).lo)
#define Moonbit_v128_hi(v) ((v).hi)
#define Moonbit_v128_load_storage(p) (*(p))
#define Moonbit_v128_store_storage(p, v) (*(p) = (v))
#endif

// detect whether compiler builtins exist for advanced bitwise operations
#ifdef __has_builtin

#if __has_builtin(__builtin_clz)
#define HAS_BUILTIN_CLZ
#endif

#if __has_builtin(__builtin_ctz)
#define HAS_BUILTIN_CTZ
#endif

#if __has_builtin(__builtin_popcount)
#define HAS_BUILTIN_POPCNT
#endif

#if __has_builtin(__builtin_sqrt)
#define HAS_BUILTIN_SQRT
#endif

#if __has_builtin(__builtin_sqrtf)
#define HAS_BUILTIN_SQRTF
#endif

#if __has_builtin(__builtin_fabs)
#define HAS_BUILTIN_FABS
#endif

#if __has_builtin(__builtin_fabsf)
#define HAS_BUILTIN_FABSF
#endif

#endif

// if there is no builtin operators, use software implementation
#ifdef HAS_BUILTIN_CLZ
static inline int32_t moonbit_clz32(int32_t x) {
  return x == 0 ? 32 : __builtin_clz(x);
}

static inline int32_t moonbit_clz64(int64_t x) {
  return x == 0 ? 64 : __builtin_clzll(x);
}

#undef HAS_BUILTIN_CLZ
#else
// table for [clz] value of 4bit integer.
static const uint8_t moonbit_clz4[] = {4, 3, 2, 2, 1, 1, 1, 1,
                                       0, 0, 0, 0, 0, 0, 0, 0};

int32_t moonbit_clz32(uint32_t x) {
  /* The ideas is to:

     1. narrow down the 4bit block where the most signficant "1" bit lies,
        using binary search
     2. find the number of leading zeros in that 4bit block via table lookup

     Different time/space tradeoff can be made here by enlarging the table
     and do less binary search.
     One benefit of the 4bit lookup table is that it can fit into a single cache
     line.
  */
  int32_t result = 0;
  if (x > 0xffff) {
    x >>= 16;
  } else {
    result += 16;
  }
  if (x > 0xff) {
    x >>= 8;
  } else {
    result += 8;
  }
  if (x > 0xf) {
    x >>= 4;
  } else {
    result += 4;
  }
  return result + moonbit_clz4[x];
}

int32_t moonbit_clz64(uint64_t x) {
  int32_t result = 0;
  if (x > 0xffffffff) {
    x >>= 32;
  } else {
    result += 32;
  }
  return result + moonbit_clz32((uint32_t)x);
}
#endif

#ifdef HAS_BUILTIN_CTZ
static inline int32_t moonbit_ctz32(int32_t x) {
  return x == 0 ? 32 : __builtin_ctz(x);
}

static inline int32_t moonbit_ctz64(int64_t x) {
  return x == 0 ? 64 : __builtin_ctzll(x);
}

#undef HAS_BUILTIN_CTZ
#else
int32_t moonbit_ctz32(int32_t x) {
  /* The algorithm comes from:

       Leiserson, Charles E. et al. “Using de Bruijn Sequences to Index a 1 in a
     Computer Word.” (1998).

     The ideas is:

     1. leave only the least significant "1" bit in the input,
        set all other bits to "0". This is achieved via [x & -x]
     2. now we have [x * n == n << ctz(x)], if [n] is a de bruijn sequence
        (every 5bit pattern occurn exactly once when you cycle through the bit
     string), we can find [ctz(x)] from the most significant 5 bits of [x * n]
 */
  static const uint32_t de_bruijn_32 = 0x077CB531;
  static const uint8_t index32[] = {0,  1,  28, 2,  29, 14, 24, 3,  30, 22, 20,
                                    15, 25, 17, 4,  8,  31, 27, 13, 23, 21, 19,
                                    16, 7,  26, 12, 18, 6,  11, 5,  10, 9};
  return (x == 0) * 32 + index32[(de_bruijn_32 * (x & -x)) >> 27];
}

int32_t moonbit_ctz64(int64_t x) {
  static const uint64_t de_bruijn_64 = 0x0218A392CD3D5DBF;
  static const uint8_t index64[] = {
      0,  1,  2,  7,  3,  13, 8,  19, 4,  25, 14, 28, 9,  34, 20, 40,
      5,  17, 26, 38, 15, 46, 29, 48, 10, 31, 35, 54, 21, 50, 41, 57,
      63, 6,  12, 18, 24, 27, 33, 39, 16, 37, 45, 47, 30, 53, 49, 56,
      62, 11, 23, 32, 36, 44, 52, 55, 61, 22, 43, 51, 60, 42, 59, 58};
  return (x == 0) * 64 + index64[(de_bruijn_64 * (x & -x)) >> 58];
}
#endif

#ifdef HAS_BUILTIN_POPCNT

#define moonbit_popcnt32 __builtin_popcount
#define moonbit_popcnt64 __builtin_popcountll
#undef HAS_BUILTIN_POPCNT

#else
int32_t moonbit_popcnt32(uint32_t x) {
  /* The classic SIMD Within A Register algorithm.
     ref: [https://nimrod.blog/posts/algorithms-behind-popcount/]
 */
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x + (x >> 4)) & 0x0F0F0F0F;
  return (x * 0x01010101) >> 24;
}

int32_t moonbit_popcnt64(uint64_t x) {
  x = x - ((x >> 1) & 0x5555555555555555);
  x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
  x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F;
  return (x * 0x0101010101010101) >> 56;
}
#endif

/* The following sqrt implementation comes from
   [musl](https://git.musl-libc.org/cgit/musl),
   with some helpers inlined to make it zero dependency.
 */
#ifdef MOONBIT_NATIVE_NO_SYS_HEADER
const uint16_t __rsqrt_tab[128] = {
    0xb451, 0xb2f0, 0xb196, 0xb044, 0xaef9, 0xadb6, 0xac79, 0xab43, 0xaa14,
    0xa8eb, 0xa7c8, 0xa6aa, 0xa592, 0xa480, 0xa373, 0xa26b, 0xa168, 0xa06a,
    0x9f70, 0x9e7b, 0x9d8a, 0x9c9d, 0x9bb5, 0x9ad1, 0x99f0, 0x9913, 0x983a,
    0x9765, 0x9693, 0x95c4, 0x94f8, 0x9430, 0x936b, 0x92a9, 0x91ea, 0x912e,
    0x9075, 0x8fbe, 0x8f0a, 0x8e59, 0x8daa, 0x8cfe, 0x8c54, 0x8bac, 0x8b07,
    0x8a64, 0x89c4, 0x8925, 0x8889, 0x87ee, 0x8756, 0x86c0, 0x862b, 0x8599,
    0x8508, 0x8479, 0x83ec, 0x8361, 0x82d8, 0x8250, 0x81c9, 0x8145, 0x80c2,
    0x8040, 0xff02, 0xfd0e, 0xfb25, 0xf947, 0xf773, 0xf5aa, 0xf3ea, 0xf234,
    0xf087, 0xeee3, 0xed47, 0xebb3, 0xea27, 0xe8a3, 0xe727, 0xe5b2, 0xe443,
    0xe2dc, 0xe17a, 0xe020, 0xdecb, 0xdd7d, 0xdc34, 0xdaf1, 0xd9b3, 0xd87b,
    0xd748, 0xd61a, 0xd4f1, 0xd3cd, 0xd2ad, 0xd192, 0xd07b, 0xcf69, 0xce5b,
    0xcd51, 0xcc4a, 0xcb48, 0xca4a, 0xc94f, 0xc858, 0xc764, 0xc674, 0xc587,
    0xc49d, 0xc3b7, 0xc2d4, 0xc1f4, 0xc116, 0xc03c, 0xbf65, 0xbe90, 0xbdbe,
    0xbcef, 0xbc23, 0xbb59, 0xba91, 0xb9cc, 0xb90a, 0xb84a, 0xb78c, 0xb6d0,
    0xb617, 0xb560,
};

/* returns a*b*2^-32 - e, with error 0 <= e < 1.  */
static inline uint32_t mul32(uint32_t a, uint32_t b) {
  return (uint64_t)a * b >> 32;
}
#endif

#ifdef MOONBIT_NATIVE_NO_SYS_HEADER
float sqrtf(float x) {
  uint32_t ix, m, m1, m0, even, ey;

  ix = *(uint32_t *)&x;
  if (ix - 0x00800000 >= 0x7f800000 - 0x00800000) {
    /* x < 0x1p-126 or inf or nan.  */
    if (ix * 2 == 0)
      return x;
    if (ix == 0x7f800000)
      return x;
    if (ix > 0x7f800000)
      return (x - x) / (x - x);
    /* x is subnormal, normalize it.  */
    x *= 0x1p23f;
    ix = *(uint32_t *)&x;
    ix -= 23 << 23;
  }

  /* x = 4^e m; with int e and m in [1, 4).  */
  even = ix & 0x00800000;
  m1 = (ix << 8) | 0x80000000;
  m0 = (ix << 7) & 0x7fffffff;
  m = even ? m0 : m1;

  /* 2^e is the exponent part of the return value.  */
  ey = ix >> 1;
  ey += 0x3f800000 >> 1;
  ey &= 0x7f800000;

  /* compute r ~ 1/sqrt(m), s ~ sqrt(m) with 2 goldschmidt iterations.  */
  static const uint32_t three = 0xc0000000;
  uint32_t r, s, d, u, i;
  i = (ix >> 17) % 128;
  r = (uint32_t)__rsqrt_tab[i] << 16;
  /* |r*sqrt(m) - 1| < 0x1p-8 */
  s = mul32(m, r);
  /* |s/sqrt(m) - 1| < 0x1p-8 */
  d = mul32(s, r);
  u = three - d;
  r = mul32(r, u) << 1;
  /* |r*sqrt(m) - 1| < 0x1.7bp-16 */
  s = mul32(s, u) << 1;
  /* |s/sqrt(m) - 1| < 0x1.7bp-16 */
  d = mul32(s, r);
  u = three - d;
  s = mul32(s, u);
  /* -0x1.03p-28 < s/sqrt(m) - 1 < 0x1.fp-31 */
  s = (s - 1) >> 6;
  /* s < sqrt(m) < s + 0x1.08p-23 */

  /* compute nearest rounded result.  */
  uint32_t d0, d1, d2;
  float y, t;
  d0 = (m << 16) - s * s;
  d1 = s - d0;
  d2 = d1 + s + 1;
  s += d1 >> 31;
  s &= 0x007fffff;
  s |= ey;
  y = *(float *)&s;
  /* handle rounding and inexact exception. */
  uint32_t tiny = d2 == 0 ? 0 : 0x01000000;
  tiny |= (d1 ^ d2) & 0x80000000;
  t = *(float *)&tiny;
  y = y + t;
  return y;
}
#endif

#ifdef MOONBIT_NATIVE_NO_SYS_HEADER
/* returns a*b*2^-64 - e, with error 0 <= e < 3.  */
static inline uint64_t mul64(uint64_t a, uint64_t b) {
  uint64_t ahi = a >> 32;
  uint64_t alo = a & 0xffffffff;
  uint64_t bhi = b >> 32;
  uint64_t blo = b & 0xffffffff;
  return ahi * bhi + (ahi * blo >> 32) + (alo * bhi >> 32);
}

double sqrt(double x) {
  uint64_t ix, top, m;

  /* special case handling.  */
  ix = *(uint64_t *)&x;
  top = ix >> 52;
  if (top - 0x001 >= 0x7ff - 0x001) {
    /* x < 0x1p-1022 or inf or nan.  */
    if (ix * 2 == 0)
      return x;
    if (ix == 0x7ff0000000000000)
      return x;
    if (ix > 0x7ff0000000000000)
      return (x - x) / (x - x);
    /* x is subnormal, normalize it.  */
    x *= 0x1p52;
    ix = *(uint64_t *)&x;
    top = ix >> 52;
    top -= 52;
  }

  /* argument reduction:
     x = 4^e m; with integer e, and m in [1, 4)
     m: fixed point representation [2.62]
     2^e is the exponent part of the result.  */
  int even = top & 1;
  m = (ix << 11) | 0x8000000000000000;
  if (even)
    m >>= 1;
  top = (top + 0x3ff) >> 1;

  /* approximate r ~ 1/sqrt(m) and s ~ sqrt(m) when m in [1,4)

     initial estimate:
     7bit table lookup (1bit exponent and 6bit significand).

     iterative approximation:
     using 2 goldschmidt iterations with 32bit int arithmetics
     and a final iteration with 64bit int arithmetics.

     details:

     the relative error (e = r0 sqrt(m)-1) of a linear estimate
     (r0 = a m + b) is |e| < 0.085955 ~ 0x1.6p-4 at best,
     a table lookup is faster and needs one less iteration
     6 bit lookup table (128b) gives |e| < 0x1.f9p-8
     7 bit lookup table (256b) gives |e| < 0x1.fdp-9
     for single and double prec 6bit is enough but for quad
     prec 7bit is needed (or modified iterations). to avoid
     one more iteration >=13bit table would be needed (16k).

     a newton-raphson iteration for r is
       w = r*r
       u = 3 - m*w
       r = r*u/2
     can use a goldschmidt iteration for s at the end or
       s = m*r

     first goldschmidt iteration is
       s = m*r
       u = 3 - s*r
       r = r*u/2
       s = s*u/2
     next goldschmidt iteration is
       u = 3 - s*r
       r = r*u/2
       s = s*u/2
     and at the end r is not computed only s.

     they use the same amount of operations and converge at the
     same quadratic rate, i.e. if
       r1 sqrt(m) - 1 = e, then
       r2 sqrt(m) - 1 = -3/2 e^2 - 1/2 e^3
     the advantage of goldschmidt is that the mul for s and r
     are independent (computed in parallel), however it is not
     "self synchronizing": it only uses the input m in the
     first iteration so rounding errors accumulate. at the end
     or when switching to larger precision arithmetics rounding
     errors dominate so the first iteration should be used.

     the fixed point representations are
       m: 2.30 r: 0.32, s: 2.30, d: 2.30, u: 2.30, three: 2.30
     and after switching to 64 bit
       m: 2.62 r: 0.64, s: 2.62, d: 2.62, u: 2.62, three: 2.62  */

  static const uint64_t three = 0xc0000000;
  uint64_t r, s, d, u, i;

  i = (ix >> 46) % 128;
  r = (uint32_t)__rsqrt_tab[i] << 16;
  /* |r sqrt(m) - 1| < 0x1.fdp-9 */
  s = mul32(m >> 32, r);
  /* |s/sqrt(m) - 1| < 0x1.fdp-9 */
  d = mul32(s, r);
  u = three - d;
  r = mul32(r, u) << 1;
  /* |r sqrt(m) - 1| < 0x1.7bp-16 */
  s = mul32(s, u) << 1;
  /* |s/sqrt(m) - 1| < 0x1.7bp-16 */
  d = mul32(s, r);
  u = three - d;
  r = mul32(r, u) << 1;
  /* |r sqrt(m) - 1| < 0x1.3704p-29 (measured worst-case) */
  r = r << 32;
  s = mul64(m, r);
  d = mul64(s, r);
  u = (three << 32) - d;
  s = mul64(s, u); /* repr: 3.61 */
  /* -0x1p-57 < s - sqrt(m) < 0x1.8001p-61 */
  s = (s - 2) >> 9; /* repr: 12.52 */
  /* -0x1.09p-52 < s - sqrt(m) < -0x1.fffcp-63 */

  /* s < sqrt(m) < s + 0x1.09p-52,
     compute nearest rounded result:
     the nearest result to 52 bits is either s or s+0x1p-52,
     we can decide by comparing (2^52 s + 0.5)^2 to 2^104 m.  */
  uint64_t d0, d1, d2;
  double y, t;
  d0 = (m << 42) - s * s;
  d1 = s - d0;
  d2 = d1 + s + 1;
  s += d1 >> 63;
  s &= 0x000fffffffffffff;
  s |= top << 52;
  y = *(double *)&s;
  return y;
}
#endif

#ifdef MOONBIT_NATIVE_NO_SYS_HEADER
double fabs(double x) {
  union {
    double f;
    uint64_t i;
  } u = {x};
  u.i &= 0x7fffffffffffffffULL;
  return u.f;
}
#endif

#ifdef MOONBIT_NATIVE_NO_SYS_HEADER
float fabsf(float x) {
  union {
    float f;
    uint32_t i;
  } u = {x};
  u.i &= 0x7fffffff;
  return u.f;
}
#endif

#ifdef _MSC_VER
/* MSVC treats syntactic division by zero as fatal error,
   even for float point numbers,
   so we have to use a constant variable to work around this */
static const int MOONBIT_ZERO = 0;
#else
#define MOONBIT_ZERO 0
#endif

#ifdef __cplusplus
}
#endif
struct _M0TPB5ArrayGRPC16string10StringViewE;

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure;

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError;

struct _M0TPB8MutLocalGORPC16string10StringViewE;

struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__;

struct _M0TPB4IterGUsRPB5ArrayGiEEE;

struct _M0TWRPC15error5ErrorEs;

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE3Err;

struct _M0DTPC15error5Error123biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError;

struct _M0TWssbEu;

struct _M0TUsiE;

struct _M0TPB8MutLocalGsE;

struct _M0TPB13StringBuilder;

struct _M0TPB8MutLocalGRPB13StringBuilderE;

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some;

struct _M0TWEOUsRPB5ArrayGiEE;

struct _M0TPB3MapGsRPB5ArrayGiEE;

struct _M0TP26biolab8bio__seq3Seq;

struct _M0BTPB6Logger;

struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE;

struct _M0TPB4IterGcE;

struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError;

struct _M0TPB6Logger;

struct _M0TWcERPC16string10StringView;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd27seqio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE2Ok;

struct _M0TPB5ArrayGUsiEE;

struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE;

struct _M0TPB8MutLocalGOsE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok;

struct _M0TPB5ArrayGiE;

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq9SeqRecordRP26biolab8bio__seq10SeqIOErrorE2Ok;

struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__;

struct _M0TPB8MutLocalGAkE;

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq10SeqIOErrorE3Err;

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError;

struct _M0TWRPC15error5ErrorEu;

struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__;

struct _M0TURPC16string10StringViewRPB6LoggerE;

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq10SeqIOErrorE2Ok;

struct _M0TPB4IterGRPC16string10StringViewE;

struct _M0TPB5EntryGsRPB5ArrayGiEE;

struct _M0TPB8MutLocalGiE;

struct _M0R38String_3a_3aiter_2eanon__u1750__l256__;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq9SeqRecordRP26biolab8bio__seq10SeqIOErrorE3Err;

struct _M0TWERPC16option6OptionGRPC16string10StringViewE;

struct _M0TPB4Show;

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE;

struct _M0TWEOc;

struct _M0TUsRPB5ArrayGiEE;

struct _M0TP26biolab8bio__seq9SeqRecord;

struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd27seqio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE3Err;

struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190;

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error;

struct _M0TUsRP26biolab8bio__seq9SeqRecordE;

struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err;

struct _M0BTPB4Show;

struct _M0TWuEu;

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE2Ok;

struct _M0TPC16string10StringView;

struct _M0TPB8MutLocalGbE;

struct _M0KTPB6LoggerTPB13StringBuilder;

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE;

struct _M0TPB5ArrayGsE;

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE;

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__;

struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest;

struct _M0TPB5ArrayGRPC16string10StringViewE {
  struct _M0TPC16string10StringView* $0;
  int32_t $1;
  
};

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure {
  moonbit_string_t $0;
  
};

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError {
  moonbit_string_t $0;
  
};

struct _M0TPB8MutLocalGORPC16string10StringViewE {
  void* $0;
  
};

struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__ {
  struct _M0TUsRPB5ArrayGiEE*(* code)(struct _M0TWEOUsRPB5ArrayGiEE*);
  struct _M0TPB8MutLocalGiE* $0;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE* $1;
  
};

struct _M0TPB4IterGUsRPB5ArrayGiEEE {
  struct _M0TWEOUsRPB5ArrayGiEE* $0;
  int64_t $1;
  
};

struct _M0TWRPC15error5ErrorEs {
  moonbit_string_t(* code)(struct _M0TWRPC15error5ErrorEs*, void*);
  
};

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE3Err {
  void* $0;
  
};

struct _M0DTPC15error5Error123biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError {
  moonbit_string_t $0;
  
};

struct _M0TWssbEu {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  
};

struct _M0TUsiE {
  moonbit_string_t $0;
  int32_t $1;
  
};

struct _M0TPB8MutLocalGsE {
  moonbit_string_t $0;
  
};

struct _M0TPB13StringBuilder {
  uint16_t* $0;
  int32_t $1;
  
};

struct _M0TPB8MutLocalGRPB13StringBuilderE {
  struct _M0TPB13StringBuilder* $0;
  
};

struct _M0TPC16string10StringView {
  moonbit_string_t $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some {
  struct _M0TPC16string10StringView $0;
  
};

struct _M0TWEOUsRPB5ArrayGiEE {
  struct _M0TUsRPB5ArrayGiEE*(* code)(struct _M0TWEOUsRPB5ArrayGiEE*);
  
};

struct _M0TPB3MapGsRPB5ArrayGiEE {
  struct _M0TPB5EntryGsRPB5ArrayGiEE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* $5;
  int32_t $6;
  
};

struct _M0TP26biolab8bio__seq3Seq {
  moonbit_string_t $0;
  
};

struct _M0BTPB6Logger {
  int32_t(* $method_0)(void*, moonbit_string_t);
  int32_t(* $method_1)(void*, moonbit_string_t, int32_t, int32_t);
  int32_t(* $method_2)(void*, struct _M0TPC16string10StringView);
  int32_t(* $method_3)(void*, int32_t);
  int32_t(* $method_4)(void*, struct _M0TPB4Show);
  int32_t(* $method_5)(void*, struct _M0TPB4Show);
  
};

struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE {
  struct _M0TUsRP26biolab8bio__seq9SeqRecordE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0TPB4IterGcE {
  struct _M0TWEOc* $0;
  int64_t $1;
  
};

struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError {
  moonbit_string_t $0;
  
};

struct _M0TPB6Logger {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0TWcERPC16string10StringView {
  struct _M0TPC16string10StringView(* code)(
    struct _M0TWcERPC16string10StringView*,
    int32_t
  );
  
};

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd27seqio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE2Ok {
  int32_t $0;
  
};

struct _M0TPB5ArrayGUsiEE {
  struct _M0TUsiE** $0;
  int32_t $1;
  
};

struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE {
  struct _M0TPB5EntryGsRPB5ArrayGiEE* $0;
  
};

struct _M0TPB8MutLocalGOsE {
  moonbit_string_t $0;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok {
  int32_t $0;
  
};

struct _M0TPB5ArrayGiE {
  int32_t* $0;
  int32_t $1;
  
};

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** $0;
  int32_t $1;
  
};

struct _M0DTPC16result6ResultGRP26biolab8bio__seq9SeqRecordRP26biolab8bio__seq10SeqIOErrorE2Ok {
  struct _M0TP26biolab8bio__seq9SeqRecord* $0;
  
};

struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TPB8MutLocalGORPC16string10StringViewE* $0;
  struct _M0TPC16string10StringView $1;
  int32_t $2;
  
};

struct _M0TPB8MutLocalGAkE {
  uint16_t* $0;
  
};

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq10SeqIOErrorE3Err {
  void* $0;
  
};

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError {
  moonbit_string_t $0;
  
};

struct _M0TWRPC15error5ErrorEu {
  int32_t(* code)(struct _M0TWRPC15error5ErrorEu*, void*);
  
};

struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  int32_t $1;
  struct _M0TPC16string10StringView $2;
  
};

struct _M0TURPC16string10StringViewRPB6LoggerE {
  struct _M0TPC16string10StringView $0;
  struct _M0TPB6Logger $1;
  
};

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq10SeqIOErrorE2Ok {
  moonbit_string_t $0;
  
};

struct _M0TPB4IterGRPC16string10StringViewE {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* $0;
  int64_t $1;
  
};

struct _M0TPB5EntryGsRPB5ArrayGiEE {
  int32_t $0;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TPB5ArrayGiE* $5;
  
};

struct _M0TPB8MutLocalGiE {
  int32_t $0;
  
};

struct _M0R38String_3a_3aiter_2eanon__u1750__l256__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  moonbit_string_t $1;
  int32_t $2;
  
};

struct _M0DTPC16result6ResultGRP26biolab8bio__seq9SeqRecordRP26biolab8bio__seq10SeqIOErrorE3Err {
  void* $0;
  
};

struct _M0TWERPC16option6OptionGRPC16string10StringViewE {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  
};

struct _M0TPB4Show {
  struct _M0BTPB4Show* $0;
  void* $1;
  
};

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE {
  int32_t $0;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TP26biolab8bio__seq9SeqRecord* $5;
  
};

struct _M0TWEOc {
  int32_t(* code)(struct _M0TWEOc*);
  
};

struct _M0TUsRPB5ArrayGiEE {
  moonbit_string_t $0;
  struct _M0TPB5ArrayGiE* $1;
  
};

struct _M0TP26biolab8bio__seq9SeqRecord {
  struct _M0TP26biolab8bio__seq3Seq* $0;
  moonbit_string_t $1;
  moonbit_string_t $2;
  moonbit_string_t $3;
  struct _M0TPB3MapGsRPB5ArrayGiEE* $4;
  
};

struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195 {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd27seqio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE3Err {
  void* $0;
  
};

struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190 {
  int32_t(* code)(struct _M0TWEOc*);
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error {
  struct moonbit_result_0(* code)(
    struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error*,
    struct _M0TWuEu*,
    struct _M0TWRPC15error5ErrorEu*
  );
  
};

struct _M0TUsRP26biolab8bio__seq9SeqRecordE {
  moonbit_string_t $0;
  struct _M0TP26biolab8bio__seq9SeqRecord* $1;
  
};

struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE {
  struct _M0TUsRPB5ArrayGiEE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err {
  void* $0;
  
};

struct _M0BTPB4Show {
  int32_t(* $method_0)(void*, struct _M0TPB6Logger);
  moonbit_string_t(* $method_1)(void*);
  
};

struct _M0TWuEu {
  int32_t(* code)(struct _M0TWuEu*, int32_t);
  
};

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE2Ok {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* $0;
  
};

struct _M0TPB8MutLocalGbE {
  int32_t $0;
  
};

struct _M0KTPB6LoggerTPB13StringBuilder {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE {
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* $5;
  int32_t $6;
  
};

struct _M0TPB5ArrayGsE {
  moonbit_string_t* $0;
  int32_t $1;
  
};

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE {
  struct _M0TP26biolab8bio__seq9SeqRecord** $0;
  int32_t $1;
  
};

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TWcERPC16string10StringView* $0;
  struct _M0TPB4IterGcE* $1;
  
};

struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest {
  moonbit_string_t $0;
  
};

struct moonbit_result_2 {
  int tag;
  union { struct _M0TP26biolab8bio__seq9SeqRecord* ok; void* err;  } data;
  
};

struct moonbit_result_1 {
  int tag;
  union { struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* ok; void* err; 
  } data;
  
};

struct moonbit_result_3 {
  int tag;
  union { moonbit_string_t ok; void* err;  } data;
  
};

struct moonbit_result_0 {
  int tag;
  union { int32_t ok; void* err;  } data;
  
};

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1202(
  struct _M0TWRPC15error5ErrorEs*,
  void*
);

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1195(
  struct _M0TWssbEu*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1190(
  struct _M0TWEOc*
);

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__args(
  
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1168(
  int32_t,
  moonbit_string_t,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1161(
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154(
  int32_t,
  moonbit_bytes_t
);

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1147(
  int32_t,
  moonbit_string_t
);

#define _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__get__cli__args__ffi moonbit_get_cli_args

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*
);

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main11seqio__test();

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main4emit(
  moonbit_string_t,
  moonbit_string_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main25records__to__json_2einner(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main24record__to__json_2einner(
  struct _M0TP26biolab8bio__seq9SeqRecord*,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(
  moonbit_string_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main12json__escape(
  moonbit_string_t
);

struct _M0TPB5ArrayGiE* _M0MP26biolab8bio__seq9SeqRecord14phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord*,
  struct _M0TPB5ArrayGiE*
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord9from__seq(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq15seqio__to__dict(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

struct moonbit_result_3 _M0FP26biolab8bio__seq12seqio__write(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t
);

struct moonbit_result_3 _M0FP26biolab8bio__seq12write__fastq(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

moonbit_string_t _M0FP26biolab8bio__seq12write__fasta(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

moonbit_string_t _M0FP26biolab8bio__seq19build__fasta__title(
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_2 _M0FP26biolab8bio__seq11seqio__read(
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_1 _M0FP26biolab8bio__seq12seqio__parse(
  moonbit_string_t,
  moonbit_string_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq14parse__genbank(
  moonbit_string_t
);

struct moonbit_result_1 _M0FP26biolab8bio__seq12parse__fastq(
  moonbit_string_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq12parse__fasta(
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0FP26biolab8bio__seq24fasta__title__to__record(
  moonbit_string_t,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord11new_2einner(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t,
  moonbit_string_t,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq9SeqRecord6length(
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

int32_t _M0MP26biolab8bio__seq3Seq6length(struct _M0TP26biolab8bio__seq3Seq*);

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

moonbit_string_t _M0MPC15array5Array2atGsE(struct _M0TPB5ArrayGsE*, int32_t);

struct _M0TPC16string10StringView _M0MPC15array5Array2atGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  int32_t
);

int32_t _M0MPC15array5Array2atGiE(struct _M0TPB5ArrayGiE*, int32_t);

int32_t _M0FPB7printlnGsE(moonbit_string_t);

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t);

struct _M0TUsRPB5ArrayGiEE* _M0MPB5Iter24nextGsRPB5ArrayGiEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE*
);

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB3Map5iter2GsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*
);

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB3Map4iterGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*
);

struct _M0TUsRPB5ArrayGiEE* _M0MPB3Map4iterGsRPB5ArrayGiEEC1984l711(
  struct _M0TWEOUsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map6lengthGsRPB5ArrayGiEE(struct _M0TPB3MapGsRPB5ArrayGiEE*);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t
);

struct _M0TPB5ArrayGiE* _M0MPB3Map3getGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t
);

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0MPB3Map3MapGsRPB5ArrayGiEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE,
  int64_t
);

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0MPB3Map3MapGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE,
  int64_t
);

int32_t _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t,
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

int32_t _M0MPB3Map3setGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGiE*
);

int32_t _M0MPB3Map15set__with__hashGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t,
  struct _M0TP26biolab8bio__seq9SeqRecord*,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGiE*,
  int32_t
);

int32_t _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map4growGsRPB5ArrayGiEE(struct _M0TPB3MapGsRPB5ArrayGiEE*);

int32_t _M0MPB3Map20rehash__place__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  int32_t,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*,
  int32_t
);

int32_t _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  int32_t,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPC13int3Int3max(int32_t, int32_t);

int32_t _M0FPB21capacity__for__length(int32_t);

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0FPB8new__mapGsRPB5ArrayGiEE(int32_t);

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FPB8new__mapGsRP26biolab8bio__seq9SeqRecordE(
  int32_t
);

int32_t _M0MPC13int3Int20next__power__of__two(int32_t);

int32_t _M0FPB21calc__grow__threshold(int32_t);

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0MPC16option6Option6unwrapGRPB5EntryGsRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t);

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPB4Iter9to__arrayGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE*
);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string6String5split(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string10StringView5split(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

void* _M0MPC16string10StringView5splitC1779l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC1775l1145(
  struct _M0TWcERPC16string10StringView*,
  int32_t
);

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t);

moonbit_string_t _M0FPB16char__to__string(int32_t);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE*,
  struct _M0TWcERPC16string10StringView*
);

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC1768l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPB4IterGcE* _M0MPC16string6String4iter(moonbit_string_t);

int32_t _M0MPC16string6String4iterC1750l256(struct _M0TWEOc*);

int32_t _M0MPC15array5Array4pushGiE(struct _M0TPB5ArrayGiE*, int32_t);

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE*,
  moonbit_string_t
);

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE*,
  struct _M0TUsiE*
);

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

int32_t _M0MPC15array5Array4pushGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC15array5Array7reallocGiE(struct _M0TPB5ArrayGiE*);

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array7reallocGUsiEE(struct _M0TPB5ArrayGUsiEE*);

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPC15array5Array7reallocGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*
);

int32_t _M0MPC15array5Array14resize__bufferGiE(
  struct _M0TPB5ArrayGiE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  int32_t
);

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPC15array5Array6lengthGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*
);

int32_t _M0MPC15array5Array6lengthGiE(struct _M0TPB5ArrayGiE*);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(
  int32_t
);

struct _M0TPB5ArrayGiE* _M0MPC15array5Array11new_2einnerGiE(int32_t);

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(int32_t);

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(
  int32_t
);

int32_t _M0MPC16string6String11has__prefix(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView11has__prefix(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

moonbit_string_t _M0MPC16string6String6repeat(moonbit_string_t, int32_t);

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(moonbit_string_t);

int64_t _M0MPC16string6String4find(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

int64_t _M0MPC16string10StringView4find(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int64_t _M0FPB18brute__force__find(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int64_t _M0FPB28boyer__moore__horspool__find(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int32_t _M0IPB13StringBuilderPB6Logger11write__view(
  struct _M0TPB13StringBuilder*,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t,
  int32_t,
  int64_t
);

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView4iterC1621l219(struct _M0TWEOc*);

moonbit_string_t _M0MPC16string10StringView9to__owned(
  struct _M0TPC16string10StringView
);

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t,
  int32_t,
  int32_t
);

int32_t _M0IPC14byte4BytePB7Default7default();

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t,
  int32_t,
  int64_t
);

#define _M0FPB19unsafe__sub__string moonbit_unsafe_bytes_sub_string

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t,
  int32_t,
  moonbit_string_t,
  int32_t,
  int32_t
);

int32_t _M0MPC14uint4UInt8to__byte(uint32_t);

int32_t* _M0MPC15array5Array6bufferGiE(struct _M0TPB5ArrayGiE*);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

moonbit_string_t* _M0MPC15array5Array6bufferGsE(struct _M0TPB5ArrayGsE*);

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*
);

struct _M0TPC16string10StringView* _M0MPC15array5Array6bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*
);

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView,
  int32_t,
  int64_t
);

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(struct _M0TWEOc*, int64_t);

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB4Iter3newGUsRPB5ArrayGiEEE(
  struct _M0TWEOUsRPB5ArrayGiEE*,
  int64_t
);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*,
  int64_t
);

moonbit_string_t _M0MPC13int3Int18to__string_2einner(int32_t, int32_t);

int32_t _M0FPB14radix__count32(uint32_t, int32_t);

int32_t _M0FPB12hex__count32(uint32_t);

int32_t _M0FPB12dec__count32(uint32_t);

int32_t _M0FPB20int__to__string__dec(uint16_t*, uint32_t, int32_t, int32_t);

int32_t _M0FPB24int__to__string__generic(
  uint16_t*,
  uint32_t,
  int32_t,
  int32_t,
  int32_t
);

int32_t _M0FPB20int__to__string__hex(uint16_t*, uint32_t, int32_t, int32_t);

int32_t _M0MPB4Iter4nextGcE(struct _M0TPB4IterGcE*);

struct _M0TUsRPB5ArrayGiEE* _M0MPB4Iter4nextGUsRPB5ArrayGiEEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE*
);

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE*
);

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void*
);

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t,
  struct _M0TPB6Logger
);

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t,
  struct _M0TPB6Logger
);

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView
);

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView
);

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder*,
  moonbit_string_t,
  int32_t,
  int32_t
);

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t,
  int32_t,
  int64_t
);

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder*,
  struct _M0TPB4Show
);

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder*,
  struct _M0TPB4Show
);

int32_t _M0FPB13finalize__acc(uint32_t);

uint32_t _M0FPB14avalanche__acc(uint32_t);

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder*,
  moonbit_string_t
);

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t*,
  int32_t,
  moonbit_string_t,
  int32_t,
  int32_t
);

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t,
  int32_t
);

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView,
  struct _M0TPB6Logger,
  int32_t
);

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE*,
  int32_t,
  int32_t
);

int32_t _M0MPC16string10StringView11unsafe__get(
  struct _M0TPC16string10StringView,
  int32_t
);

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView,
  int32_t,
  int64_t
);

int32_t _M0MPC16string10StringView6length(struct _M0TPC16string10StringView);

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t);

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t);

int32_t _M0IPC14byte4BytePB3Sub3sub(int32_t, int32_t);

int32_t _M0IPC14byte4BytePB3Mod3mod(int32_t, int32_t);

int32_t _M0IPC14byte4BytePB3Div3div(int32_t, int32_t);

int32_t _M0IPC14byte4BytePB3Add3add(int32_t, int32_t);

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t);

int32_t _M0FPB32code__point__of__surrogate__pair(int32_t, int32_t);

int32_t _M0MPC16string6String20char__length_2einner(
  moonbit_string_t,
  int32_t,
  int64_t
);

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t);

int32_t _M0MPC16uint166UInt1622is__leading__surrogate(int32_t);

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder*,
  int32_t
);

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder*,
  int32_t
);

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t);

uint32_t _M0MPC14char4Char8to__uint(int32_t);

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder*
);

int32_t _M0IPC16uint166UInt16PB7Default7default();

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t*,
  int32_t,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t*,
  int32_t,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t
);

int32_t _M0MPC14byte4Byte8to__char(int32_t);

int32_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGiE(
  int32_t*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(
  struct _M0TPC16string10StringView*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder*,
  moonbit_string_t
);

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder*,
  int32_t
);

int32_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGiE(
  int32_t*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

moonbit_string_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGsE(
  moonbit_string_t*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TUsiE** _M0MPB18UninitializedArray23unsafe__make__and__blitGUsiEE(
  struct _M0TUsiE**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23unsafe__make__and__blitGRPC16string10StringViewE(
  struct _M0TPC16string10StringView*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGiE(
  int32_t*,
  int32_t,
  int32_t*,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGsE(
  moonbit_string_t*,
  int32_t,
  moonbit_string_t*,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGUsiEE(
  struct _M0TUsiE**,
  int32_t,
  struct _M0TUsiE**,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord**,
  int32_t,
  struct _M0TP26biolab8bio__seq9SeqRecord**,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGRPC16string10StringViewE(
  struct _M0TPC16string10StringView*,
  int32_t,
  struct _M0TPC16string10StringView*,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGkE(
  uint16_t*,
  int32_t,
  uint16_t*,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGiEE(
  int32_t*,
  int32_t,
  int32_t*,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGsEE(
  moonbit_string_t*,
  int32_t,
  moonbit_string_t*,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGUsiEEE(
  struct _M0TUsiE**,
  int32_t,
  struct _M0TUsiE**,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TP26biolab8bio__seq9SeqRecord**,
  int32_t,
  struct _M0TP26biolab8bio__seq9SeqRecord**,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPC16string10StringViewEE(
  struct _M0TPC16string10StringView*,
  int32_t,
  struct _M0TPC16string10StringView*,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray6lengthGiE(int32_t*);

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t*);

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(struct _M0TUsiE**);

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord**
);

int32_t _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(
  struct _M0TPC16string10StringView*
);

uint32_t _M0FPB13consume4__acc(uint32_t, uint32_t);

uint32_t _M0FPB4rotl(uint32_t, int32_t);

int32_t _M0IPB7FailurePB4Show6output(void*, struct _M0TPB6Logger);

int32_t _M0MPB6Logger13write__objectGsE(
  struct _M0TPB6Logger,
  moonbit_string_t
);

moonbit_string_t _M0FPC15abort5abortGsE(moonbit_string_t);

int32_t _M0FPC15abort5abortGuE(moonbit_string_t);

uint16_t* _M0FPC15abort5abortGAkE(moonbit_string_t);

int32_t _M0FPC15abort5abortGiE(moonbit_string_t);

int32_t* _M0FPC15abort5abortGRPB18UninitializedArrayGiEE(moonbit_string_t);

moonbit_string_t* _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(
  moonbit_string_t
);

struct _M0TUsiE** _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(
  moonbit_string_t
);

struct _M0TPC16string10StringView _M0FPC15abort5abortGRPC16string10StringViewE(
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(
  moonbit_string_t
);

struct _M0TPC16string10StringView* _M0FPC15abort5abortGRPB18UninitializedArrayGRPC16string10StringViewEE(
  moonbit_string_t
);

moonbit_string_t _M0FP15Error10to__string(void*);

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void*,
  struct _M0TPB4Show
);

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void*,
  struct _M0TPB4Show
);

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void*,
  int32_t
);

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void*,
  struct _M0TPC16string10StringView
);

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void*,
  moonbit_string_t,
  int32_t,
  int32_t
);

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void*,
  moonbit_string_t
);

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_44 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 103, 101, 
    110, 98, 97, 110, 107, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_1 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 115, 107, 
    105, 112, 112, 101, 100, 32, 116, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[1]; 
} const moonbit_string_literal_0 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 0, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_24 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 44, 32, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_25 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 58, 32, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[63]; 
} const moonbit_string_literal_12 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 62, 62, 111, 
    110, 108, 121, 32, 79, 110, 101, 32, 114, 101, 99, 111, 114, 100, 
    32, 111, 110, 108, 121, 10, 65, 84, 71, 71, 67, 67, 65, 84, 84, 71, 
    84, 65, 65, 84, 71, 71, 71, 67, 67, 71, 67, 84, 71, 65, 65, 65, 71, 
    71, 71, 84, 71, 67, 67, 67, 71, 65, 84, 65, 71, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_71 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 43, 10, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_66 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 102, 97, 
    115, 116, 113, 45, 115, 97, 110, 103, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_45 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 123, 125, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_21 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 123, 34, 
    101, 114, 114, 111, 114, 34, 58, 32, 110, 117, 108, 108, 125, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_90 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 117, 114, 114, 111, 103, 97, 116, 
    101, 32, 112, 97, 105, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[27]; 
} const moonbit_string_literal_73 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 78, 111, 
    32, 114, 101, 99, 111, 114, 100, 115, 32, 102, 111, 117, 110, 100, 
    32, 105, 110, 32, 104, 97, 110, 100, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_64 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 60, 117, 
    110, 107, 110, 111, 119, 110, 32, 100, 101, 115, 99, 114, 105, 112, 
    116, 105, 111, 110, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_63 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 60, 117, 
    110, 107, 110, 111, 119, 110, 32, 110, 97, 109, 101, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[53]; 
} const moonbit_string_literal_103 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 52, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_40 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 102, 113, 50, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_59 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 110, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[112]; 
} const moonbit_string_literal_101 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 111, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 115, 101, 113, 105, 111, 95, 109, 97, 105, 110, 95, 98, 
    108, 97, 99, 107, 98, 111, 120, 95, 116, 101, 115, 116, 46, 77, 111, 
    111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 
    114, 73, 110, 116, 101, 114, 110, 97, 108, 83, 107, 105, 112, 84, 
    101, 115, 116, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 
    116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 
    108, 83, 107, 105, 112, 84, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[31]; 
} const moonbit_string_literal_84 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 30, 114, 97, 
    100, 105, 120, 32, 109, 117, 115, 116, 32, 98, 101, 32, 98, 101, 
    116, 119, 101, 101, 110, 32, 50, 32, 97, 110, 100, 32, 51, 54, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_74 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 77, 111, 
    114, 101, 32, 116, 104, 97, 110, 32, 111, 110, 101, 32, 114, 101, 
    99, 111, 114, 100, 32, 102, 111, 117, 110, 100, 32, 105, 110, 32, 
    104, 97, 110, 100, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_28 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 119, 114, 
    105, 116, 116, 101, 110, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_13 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 73, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_3 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 114, 101, 115, 117, 108, 116, 34, 
    44, 34, 102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_85 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_61 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 112, 104, 
    114, 101, 100, 95, 113, 117, 97, 108, 105, 116, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_20 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 114, 101, 
    97, 100, 95, 115, 105, 110, 103, 108, 101, 95, 114, 101, 99, 111, 
    114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_18 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 108, 111, 119, 101, 
    114, 99, 97, 115, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_47 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 61, 61, 
    61, 32, 67, 65, 83, 69, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_96 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 44, 32, 
    108, 101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[561]; 
} const moonbit_string_literal_43 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 560, 76, 79, 
    67, 85, 83, 32, 32, 32, 32, 32, 32, 32, 78, 77, 95, 48, 48, 49, 32, 
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 52, 56, 32, 98, 
    112, 32, 32, 32, 32, 68, 78, 65, 32, 32, 32, 32, 32, 108, 105, 110, 
    101, 97, 114, 32, 32, 32, 80, 82, 73, 32, 48, 49, 45, 74, 65, 78, 
    45, 50, 48, 50, 48, 10, 68, 69, 70, 73, 78, 73, 84, 73, 79, 78, 32, 
    32, 84, 101, 115, 116, 32, 103, 101, 110, 101, 32, 102, 111, 114, 
    32, 99, 111, 109, 112, 97, 114, 105, 115, 111, 110, 46, 10, 65, 67, 
    67, 69, 83, 83, 73, 79, 78, 32, 32, 32, 78, 77, 95, 48, 48, 49, 10, 
    86, 69, 82, 83, 73, 79, 78, 32, 32, 32, 32, 32, 78, 77, 95, 48, 48, 
    49, 46, 49, 10, 83, 79, 85, 82, 67, 69, 32, 32, 32, 32, 32, 32, 84, 
    101, 115, 116, 32, 111, 114, 103, 97, 110, 105, 115, 109, 10, 32, 
    32, 79, 82, 71, 65, 78, 73, 83, 77, 32, 32, 84, 101, 115, 116, 32, 
    111, 114, 103, 97, 110, 105, 115, 109, 10, 32, 32, 32, 32, 32, 32, 
    32, 32, 32, 32, 32, 32, 69, 117, 107, 97, 114, 121, 111, 116, 97, 
    59, 32, 84, 101, 115, 116, 46, 10, 70, 69, 65, 84, 85, 82, 69, 83, 
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 76, 111, 99, 
    97, 116, 105, 111, 110, 47, 81, 117, 97, 108, 105, 102, 105, 101, 
    114, 115, 10, 32, 32, 32, 32, 32, 115, 111, 117, 114, 99, 101, 32, 
    32, 32, 32, 32, 32, 32, 32, 32, 32, 49, 46, 46, 52, 56, 10, 32, 32, 
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 
    32, 32, 47, 111, 114, 103, 97, 110, 105, 115, 109, 61, 34, 84, 101, 
    115, 116, 32, 111, 114, 103, 97, 110, 105, 115, 109, 34, 10, 32, 
    32, 32, 32, 32, 103, 101, 110, 101, 32, 32, 32, 32, 32, 32, 32, 32, 
    32, 32, 32, 32, 49, 46, 46, 52, 56, 10, 32, 32, 32, 32, 32, 32, 32, 
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 47, 103, 
    101, 110, 101, 61, 34, 116, 101, 115, 116, 95, 103, 101, 110, 101, 
    34, 10, 32, 32, 32, 32, 32, 67, 68, 83, 32, 32, 32, 32, 32, 32, 32, 
    32, 32, 32, 32, 32, 32, 49, 46, 46, 52, 56, 10, 32, 32, 32, 32, 32, 
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 47, 
    116, 114, 97, 110, 115, 108, 97, 116, 105, 111, 110, 61, 34, 77, 
    75, 82, 72, 80, 71, 83, 84, 72, 34, 10, 79, 82, 73, 71, 73, 78, 10, 
    32, 32, 32, 32, 32, 32, 32, 32, 49, 32, 97, 116, 103, 99, 97, 116, 
    103, 99, 97, 116, 32, 103, 99, 97, 116, 103, 99, 97, 116, 103, 99, 
    32, 97, 116, 103, 99, 97, 116, 103, 99, 97, 116, 32, 103, 99, 97, 
    116, 103, 99, 97, 116, 103, 99, 32, 97, 116, 103, 99, 97, 116, 103, 
    99, 10, 47, 47, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_16 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    115, 116, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_93 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 111, 
    117, 110, 100, 115, 32, 99, 104, 101, 99, 107, 32, 102, 97, 105, 
    108, 101, 100, 58, 32, 97, 108, 108, 111, 99, 97, 116, 101, 95, 108, 
    101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[29]; 
} const moonbit_string_literal_80 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 28, 32, 100, 
    111, 110, 39, 116, 32, 109, 97, 116, 99, 104, 32, 105, 110, 32, 70, 
    65, 83, 84, 81, 32, 114, 101, 99, 111, 114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[56]; 
} const moonbit_string_literal_91 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 55, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 116, 97, 114, 116, 32, 111, 114, 
    32, 101, 110, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 114, 
    32, 83, 116, 114, 105, 110, 103, 58, 58, 99, 111, 100, 101, 112, 
    111, 105, 110, 116, 95, 108, 101, 110, 103, 116, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_65 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 85, 110, 
    107, 110, 111, 119, 110, 32, 102, 111, 114, 109, 97, 116, 58, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_30 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 97, 99, 
    103, 116, 97, 99, 103, 116, 97, 99, 103, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_32 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 119, 114, 
    105, 116, 116, 101, 110, 50, 32, 108, 111, 119, 101, 114, 99, 97, 
    115, 101, 32, 105, 110, 112, 117, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[75]; 
} const moonbit_string_literal_22 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 74, 123, 34, 
    101, 114, 114, 111, 114, 34, 58, 32, 34, 86, 97, 108, 117, 101, 69, 
    114, 114, 111, 114, 34, 44, 32, 34, 109, 101, 115, 115, 97, 103, 
    101, 34, 58, 32, 34, 77, 111, 114, 101, 32, 116, 104, 97, 110, 32, 
    111, 110, 101, 32, 114, 101, 99, 111, 114, 100, 32, 102, 111, 117, 
    110, 100, 32, 105, 110, 32, 104, 97, 110, 100, 108, 101, 34, 125, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_99 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 41, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_36 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_26 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 116, 111, 
    95, 100, 105, 99, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_86 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 48, 49, 
    50, 51, 52, 53, 54, 55, 56, 57, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 119, 120, 121, 122, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_72 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 32, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[51]; 
} const moonbit_string_literal_100 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 50, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 73, 110, 115, 112, 101, 
    99, 116, 69, 114, 114, 111, 114, 46, 73, 110, 115, 112, 101, 99, 
    116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_97 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 46, 108, 101, 110, 103, 116, 104, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_76 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 10, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_34 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 119, 114, 
    105, 116, 101, 95, 102, 97, 115, 116, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_52 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 44, 32, 
    34, 100, 101, 115, 99, 114, 105, 112, 116, 105, 111, 110, 34, 58, 
    32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_7 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 45, 45, 
    45, 45, 45, 32, 69, 78, 68, 32, 77, 79, 79, 78, 32, 84, 69, 83, 84, 
    32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_94 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_8 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 115, 116, 97, 114, 116, 34, 44, 34, 
    102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_58 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 34, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_49 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 61, 61, 
    61, 32, 69, 78, 68, 32, 61, 61, 61, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_4 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 34, 
    105, 110, 100, 101, 120, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_33 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 119, 114, 
    105, 116, 116, 101, 110, 51, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_62 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 60, 117, 
    110, 107, 110, 111, 119, 110, 32, 105, 100, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_51 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 44, 32, 
    34, 110, 97, 109, 101, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_35 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    115, 116, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_98 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 70, 97, 
    105, 108, 117, 114, 101, 40, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_50 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 123, 34, 
    105, 100, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_82 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 114, 101, 
    112, 101, 97, 116, 32, 114, 101, 115, 117, 108, 116, 32, 116, 111, 
    111, 32, 108, 97, 114, 103, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_70 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 32, 100, 
    111, 110, 39, 116, 32, 109, 97, 116, 99, 104, 32, 105, 110, 32, 114, 
    101, 99, 111, 114, 100, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_41 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 102, 113, 
    50, 32, 108, 111, 119, 101, 114, 99, 97, 115, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[135]; 
} const moonbit_string_literal_10 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 134, 62, 108, 
    111, 119, 101, 114, 49, 32, 108, 111, 119, 101, 114, 99, 97, 115, 
    101, 32, 115, 101, 113, 117, 101, 110, 99, 101, 10, 97, 116, 103, 
    103, 99, 99, 97, 116, 116, 103, 116, 97, 97, 116, 103, 103, 103, 
    99, 99, 103, 99, 116, 103, 97, 97, 97, 103, 103, 103, 116, 103, 99, 
    99, 99, 103, 97, 116, 97, 103, 10, 62, 108, 111, 119, 101, 114, 50, 
    32, 77, 105, 120, 101, 100, 67, 97, 115, 101, 32, 83, 101, 113, 117, 
    101, 110, 99, 101, 10, 97, 116, 103, 71, 67, 67, 97, 116, 116, 71, 
    84, 65, 97, 116, 103, 71, 71, 67, 99, 103, 99, 116, 71, 65, 65, 97, 
    103, 103, 103, 116, 103, 99, 99, 67, 71, 65, 84, 65, 71, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_6 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 125, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_31 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 119, 114, 
    105, 116, 116, 101, 110, 50, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_68 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 83, 101, 
    113, 117, 101, 110, 99, 101, 32, 108, 101, 110, 103, 116, 104, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[35]; 
} const moonbit_string_literal_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 34, 45, 45, 
    45, 45, 45, 32, 66, 69, 71, 73, 78, 32, 77, 79, 79, 78, 32, 84, 69, 
    83, 84, 32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_87 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 114, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_60 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 116, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_46 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 112, 97, 
    114, 115, 101, 95, 103, 101, 110, 98, 97, 110, 107, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_95 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    100, 115, 116, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_92 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 105, 110, 
    118, 97, 108, 105, 100, 32, 99, 111, 100, 101, 32, 112, 111, 105, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_54 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 44, 32, 
    34, 108, 101, 110, 103, 116, 104, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_42 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 119, 114, 
    105, 116, 101, 95, 102, 97, 115, 116, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[20]; 
} const moonbit_string_literal_23 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 114, 101, 
    97, 100, 95, 109, 117, 108, 116, 105, 112, 108, 101, 95, 101, 114, 
    114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_5 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 44, 34, 
    109, 101, 115, 115, 97, 103, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_67 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 77, 105, 
    115, 115, 105, 110, 103, 32, 112, 104, 114, 101, 100, 95, 113, 117, 
    97, 108, 105, 116, 121, 32, 97, 110, 110, 111, 116, 97, 116, 105, 
    111, 110, 32, 102, 111, 114, 32, 114, 101, 99, 111, 114, 100, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_75 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 103, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[118]; 
} const moonbit_string_literal_11 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 117, 62, 109, 
    117, 108, 116, 105, 49, 32, 77, 117, 108, 116, 105, 45, 108, 105, 
    110, 101, 32, 115, 101, 113, 117, 101, 110, 99, 101, 10, 65, 84, 
    71, 71, 67, 67, 65, 84, 84, 71, 84, 65, 10, 65, 84, 71, 71, 71, 67, 
    67, 71, 67, 84, 71, 65, 65, 10, 65, 71, 71, 71, 84, 71, 67, 67, 67, 
    71, 65, 84, 65, 71, 10, 62, 109, 117, 108, 116, 105, 50, 32, 65, 
    110, 111, 116, 104, 101, 114, 32, 109, 117, 108, 116, 105, 45, 108, 
    105, 110, 101, 10, 65, 67, 71, 84, 10, 65, 67, 71, 84, 10, 65, 67, 
    71, 84, 10, 65, 67, 71, 84, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_83 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 73, 110, 
    118, 97, 108, 105, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 
    114, 32, 86, 105, 101, 119, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[70]; 
} const moonbit_string_literal_14 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 69, 64, 115, 
    101, 113, 49, 32, 68, 101, 115, 99, 114, 105, 112, 116, 105, 111, 
    110, 32, 102, 111, 114, 32, 115, 101, 113, 49, 10, 65, 84, 71, 71, 
    67, 67, 65, 84, 84, 71, 84, 65, 65, 84, 71, 71, 71, 67, 67, 71, 67, 
    84, 71, 65, 65, 65, 71, 71, 71, 84, 71, 67, 67, 67, 71, 65, 84, 65, 
    71, 10, 43, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_69 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 32, 97, 
    110, 100, 32, 113, 117, 97, 108, 105, 116, 121, 32, 108, 101, 110, 
    103, 116, 104, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[40]; 
} const moonbit_string_literal_27 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 39, 65, 84, 
    71, 71, 67, 67, 65, 84, 84, 71, 84, 65, 65, 84, 71, 71, 71, 67, 67, 
    71, 67, 84, 71, 65, 65, 65, 71, 71, 71, 84, 71, 67, 67, 67, 71, 65, 
    84, 65, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[157]; 
} const moonbit_string_literal_9 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 156, 62, 115, 
    101, 113, 49, 32, 68, 101, 115, 99, 114, 105, 112, 116, 105, 111, 
    110, 32, 102, 111, 114, 32, 115, 101, 113, 49, 10, 65, 84, 71, 71, 
    67, 67, 65, 84, 84, 71, 84, 65, 65, 84, 71, 71, 71, 67, 67, 71, 67, 
    84, 71, 65, 65, 65, 71, 71, 71, 84, 71, 67, 67, 67, 71, 65, 84, 65, 
    71, 10, 62, 115, 101, 113, 50, 32, 65, 110, 111, 116, 104, 101, 114, 
    32, 100, 101, 115, 99, 114, 105, 112, 116, 105, 111, 110, 10, 65, 
    84, 71, 71, 67, 67, 65, 84, 84, 71, 84, 65, 65, 84, 71, 71, 71, 67, 
    67, 71, 67, 84, 71, 65, 65, 65, 71, 71, 71, 84, 71, 67, 67, 67, 71, 
    65, 84, 65, 71, 10, 62, 115, 101, 113, 51, 10, 65, 67, 71, 84, 65, 
    67, 71, 84, 65, 67, 71, 84, 65, 67, 71, 84, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_78 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 69, 120, 
    112, 101, 99, 116, 101, 100, 32, 39, 64, 39, 32, 97, 116, 32, 115, 
    116, 97, 114, 116, 32, 111, 102, 32, 70, 65, 83, 84, 81, 32, 104, 
    101, 97, 100, 101, 114, 44, 32, 103, 111, 116, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_38 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 102, 113, 
    49, 32, 119, 114, 105, 116, 116, 101, 110, 32, 102, 97, 115, 116, 
    113, 32, 114, 101, 99, 111, 114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_37 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 102, 113, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[110]; 
} const moonbit_string_literal_104 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 109, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 115, 101, 113, 105, 111, 95, 109, 97, 105, 110, 95, 98, 
    108, 97, 99, 107, 98, 111, 120, 95, 116, 101, 115, 116, 46, 77, 111, 
    111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 
    114, 73, 110, 116, 101, 114, 110, 97, 108, 74, 115, 69, 114, 114, 
    111, 114, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 
    68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 
    74, 115, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_89 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 92, 117, 123, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_17 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 109, 117, 108, 116, 
    105, 112, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_55 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 44, 32, 
    34, 113, 117, 97, 108, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_19 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 109, 117, 108, 116, 
    105, 108, 105, 110, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_39 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 97, 99, 
    103, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[48]; 
} const moonbit_string_literal_79 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 47, 69, 120, 
    112, 101, 99, 116, 101, 100, 32, 39, 43, 39, 32, 97, 116, 32, 115, 
    116, 97, 114, 116, 32, 111, 102, 32, 70, 65, 83, 84, 81, 32, 115, 
    101, 112, 97, 114, 97, 116, 111, 114, 44, 32, 103, 111, 116, 58, 
    32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[98]; 
} const moonbit_string_literal_15 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 97, 10, 64, 
    115, 101, 113, 50, 32, 65, 110, 111, 116, 104, 101, 114, 32, 100, 
    101, 115, 99, 114, 105, 112, 116, 105, 111, 110, 10, 65, 67, 71, 
    84, 65, 67, 71, 84, 10, 43, 10, 73, 73, 73, 73, 73, 73, 73, 73, 10, 
    64, 115, 101, 113, 51, 32, 78, 111, 32, 113, 117, 97, 108, 105, 116, 
    121, 32, 119, 111, 114, 116, 104, 32, 109, 101, 110, 116, 105, 111, 
    110, 105, 110, 103, 10, 97, 116, 103, 103, 99, 99, 10, 43, 10, 59, 
    59, 59, 59, 59, 59, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_77 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 60, 117, 
    110, 107, 110, 111, 119, 110, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_88 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_56 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 34, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[30]; 
} const moonbit_string_literal_29 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 29, 119, 114, 
    105, 116, 116, 101, 110, 49, 32, 70, 105, 114, 115, 116, 32, 119, 
    114, 105, 116, 116, 101, 110, 32, 114, 101, 99, 111, 114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_53 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 32, 
    34, 115, 101, 113, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_81 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 110, 101, 
    103, 97, 116, 105, 118, 101, 32, 114, 101, 112, 101, 97, 116, 32, 
    99, 111, 117, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_57 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 92, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_48 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 32, 61, 
    61, 61, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_102 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 83, 101, 
    113, 73, 79, 69, 114, 114, 111, 114, 46, 83, 101, 113, 73, 79, 69, 
    114, 114, 111, 114, 0
  };

struct moonbit_object const moonbit_constant_constructor_0 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0)
  };

struct { int32_t rc; uint32_t meta; struct _M0TWRPC15error5ErrorEs data; 
} const _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1202$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1202
  };

struct {
  int32_t rc;
  uint32_t meta;
  struct _M0TWcERPC16string10StringView data;
  
} const _M0MPC16string10StringView5splitC1775l1145$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0MPC16string10StringView5splitC1775l1145
  };

uint32_t const moonbit_layout_table_data[130] =
  {
    sizeof(struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190)
    / 4, 1,
    offsetof(struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190, $1)
    / 4,
    sizeof(struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195)
    / 4, 1,
    offsetof(struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195, $1)
    / 4,
    sizeof(struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest, $0)
    / 4, sizeof(struct _M0TPB5ArrayGUsiEE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGUsiEE, $0) / 4, sizeof(struct _M0TUsiE) / 4,
    1, offsetof(struct _M0TUsiE, $0) / 4, sizeof(struct _M0TPB5ArrayGsE) / 4,
    1, offsetof(struct _M0TPB5ArrayGsE, $0) / 4,
    sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE, $0) / 4,
    sizeof(struct _M0TPB5ArrayGiE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGiE, $0) / 4,
    sizeof(struct _M0TP26biolab8bio__seq9SeqRecord) / 4, 5,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $0) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $1) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $2) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $3) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $4) / 4,
    sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError, $0)
    / 4, sizeof(struct _M0TPB8MutLocalGsE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGsE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGAkE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGAkE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGOsE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGOsE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGRPB13StringBuilderE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGRPB13StringBuilderE, $0) / 4,
    sizeof(struct _M0TP26biolab8bio__seq3Seq) / 4, 1,
    offsetof(struct _M0TP26biolab8bio__seq3Seq, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE, $0) / 4,
    sizeof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__)
    / 4, 2,
    offsetof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__, $0)
    / 4,
    offsetof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__, $1)
    / 4, sizeof(struct _M0TUsRPB5ArrayGiEE) / 4, 2,
    offsetof(struct _M0TUsRPB5ArrayGiEE, $0) / 4,
    offsetof(struct _M0TUsRPB5ArrayGiEE, $1) / 4,
    sizeof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE) / 4, 
    3, offsetof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE, $5) / 4,
    sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE) / 4, 3,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRPB5ArrayGiEE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGiEE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGiEE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE, $5) / 4,
    sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0TPC16string10StringView) / 4, 1,
    offsetof(struct _M0TPC16string10StringView, $0) / 4,
    sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some) / 4,
    1,
    (offsetof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some, $0)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE) / 4, 
    1, offsetof(struct _M0TPB8MutLocalGORPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__) / 4, 
    2,
    offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__, $0)
    / 4,
    (offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__, $1)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4,
    sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__)
    / 4, 2,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__, $0)
    / 4,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__, $1)
    / 4, sizeof(struct _M0TPB4IterGRPC16string10StringViewE) / 4, 1,
    offsetof(struct _M0TPB4IterGRPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R38String_3a_3aiter_2eanon__u1750__l256__) / 4, 
    2, offsetof(struct _M0R38String_3a_3aiter_2eanon__u1750__l256__, $0) / 4,
    offsetof(struct _M0R38String_3a_3aiter_2eanon__u1750__l256__, $1) / 4,
    sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__) / 4, 
    2,
    offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__, $0) / 4,
    (offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__, $2)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB4IterGcE) / 4, 1,
    offsetof(struct _M0TPB4IterGcE, $0) / 4,
    sizeof(struct _M0TPB4IterGUsRPB5ArrayGiEEE) / 4, 1,
    offsetof(struct _M0TPB4IterGUsRPB5ArrayGiEEE, $0) / 4,
    sizeof(struct _M0TPB6Logger) / 4, 2,
    offsetof(struct _M0TPB6Logger, $0) / 4,
    offsetof(struct _M0TPB6Logger, $1) / 4,
    sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE) / 4, 3,
    (offsetof(struct _M0TURPC16string10StringViewRPB6LoggerE, $0)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4,
    (offsetof(struct _M0TURPC16string10StringViewRPB6LoggerE, $1)
     + offsetof(struct _M0TPB6Logger, $0))
    / 4,
    (offsetof(struct _M0TURPC16string10StringViewRPB6LoggerE, $1)
     + offsetof(struct _M0TPB6Logger, $1))
    / 4, sizeof(struct _M0TPB13StringBuilder) / 4, 1,
    offsetof(struct _M0TPB13StringBuilder, $0) / 4,
    sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE, $0) / 4
  };

struct { int32_t rc; uint32_t meta; struct _M0BTPB6Logger data; 
} _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id$object =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    {.$method_0 = _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger,
       .$method_1 = _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE,
       .$method_2 = _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger,
       .$method_3 = _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger,
       .$method_4 = _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE,
       .$method_5 = _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE}
  };

struct _M0BTPB6Logger* _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id =
  &_M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id$object.data;

int64_t _M0MPB4Iter4nextN6constrS9980GcE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GcE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GUsRPB5ArrayGiEEE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GUsRPB5ArrayGiEEE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GcE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GUsRPB5ArrayGiEEE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE = 0ll;

int64_t _M0FPB28boyer__moore__horspool__findN6constrS9990 = 0ll;

int64_t _M0FPB18brute__force__findN6constrS9991 = 0ll;

uint16_t* _M0FP26biolab8bio__seq16uppercase__table;

int32_t _M0FP26biolab8bio__seq18fasta__line__width = 60;

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1223,
  moonbit_string_t _M0L8filenameS1192,
  int32_t _M0L5indexS1194
) {
  struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190* _closure_2932;
  struct _M0TWEOc* _M0L13handle__startS1190;
  struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195* _closure_2933;
  struct _M0TWssbEu* _M0L14handle__resultS1195;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1202;
  void* _M0L11_2atry__errS1217;
  struct moonbit_result_0 _tmp_2935;
  int32_t _handle__error__result_2936;
  int32_t _M0L6_2atmpS2665;
  void* _M0L3errS1218;
  moonbit_string_t _M0L4nameS1220;
  struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest* _M0L36_2aMoonBitTestDriverInternalSkipTestS1221;
  moonbit_string_t _M0L8_2afieldS2677;
  int32_t _M0L6_2acntS2893;
  moonbit_string_t _M0L7_2anameS1222;
  #line 515 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  moonbit_incref(_M0L8filenameS1192);
  _closure_2932
  = (struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190*)moonbit_malloc(sizeof(struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190));
  Moonbit_object_header(_closure_2932)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 0, 0);
  _closure_2932->code
  = &_M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1190;
  _closure_2932->$0 = _M0L5indexS1194;
  _closure_2932->$1 = _M0L8filenameS1192;
  _M0L13handle__startS1190 = (struct _M0TWEOc*)_closure_2932;
  moonbit_incref(_M0L8filenameS1192);
  _closure_2933
  = (struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195*)moonbit_malloc(sizeof(struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195));
  Moonbit_object_header(_closure_2933)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 3, 0);
  _closure_2933->code
  = &_M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1195;
  _closure_2933->$0 = _M0L5indexS1194;
  _closure_2933->$1 = _M0L8filenameS1192;
  _M0L14handle__resultS1195 = (struct _M0TWssbEu*)_closure_2933;
  _M0L17error__to__stringS1202
  = (struct _M0TWRPC15error5ErrorEs*)&_M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1202$closure.data;
  #line 557 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _tmp_2935
  = _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(_M0L12async__testsS1223, _M0L8filenameS1192, _M0L5indexS1194, _M0L13handle__startS1190, _M0L14handle__resultS1195, _M0L17error__to__stringS1202);
  if (_tmp_2935.tag) {
    int32_t const _M0L5_2aokS2674 = _tmp_2935.data.ok;
    _handle__error__result_2936 = _M0L5_2aokS2674;
  } else {
    void* const _M0L6_2aerrS2675 = _tmp_2935.data.err;
    moonbit_decref(_M0L17error__to__stringS1202);
    moonbit_decref(_M0L13handle__startS1190);
    _M0L11_2atry__errS1217 = _M0L6_2aerrS2675;
    goto join_1216;
  }
  if (_handle__error__result_2936) {
    moonbit_decref(_M0L17error__to__stringS1202);
    moonbit_decref(_M0L13handle__startS1190);
    _M0L6_2atmpS2665 = 1;
  } else {
    struct moonbit_result_0 _tmp_2937;
    int32_t _handle__error__result_2938;
    #line 560 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
    _tmp_2937
    = _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(_M0L12async__testsS1223, _M0L8filenameS1192, _M0L5indexS1194, _M0L13handle__startS1190, _M0L14handle__resultS1195, _M0L17error__to__stringS1202);
    if (_tmp_2937.tag) {
      int32_t const _M0L5_2aokS2672 = _tmp_2937.data.ok;
      _handle__error__result_2938 = _M0L5_2aokS2672;
    } else {
      void* const _M0L6_2aerrS2673 = _tmp_2937.data.err;
      moonbit_decref(_M0L17error__to__stringS1202);
      moonbit_decref(_M0L13handle__startS1190);
      _M0L11_2atry__errS1217 = _M0L6_2aerrS2673;
      goto join_1216;
    }
    if (_handle__error__result_2938) {
      moonbit_decref(_M0L17error__to__stringS1202);
      moonbit_decref(_M0L13handle__startS1190);
      _M0L6_2atmpS2665 = 1;
    } else {
      struct moonbit_result_0 _tmp_2939;
      int32_t _handle__error__result_2940;
      #line 563 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _tmp_2939
      = _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(_M0L12async__testsS1223, _M0L8filenameS1192, _M0L5indexS1194, _M0L13handle__startS1190, _M0L14handle__resultS1195, _M0L17error__to__stringS1202);
      if (_tmp_2939.tag) {
        int32_t const _M0L5_2aokS2670 = _tmp_2939.data.ok;
        _handle__error__result_2940 = _M0L5_2aokS2670;
      } else {
        void* const _M0L6_2aerrS2671 = _tmp_2939.data.err;
        moonbit_decref(_M0L17error__to__stringS1202);
        moonbit_decref(_M0L13handle__startS1190);
        _M0L11_2atry__errS1217 = _M0L6_2aerrS2671;
        goto join_1216;
      }
      if (_handle__error__result_2940) {
        moonbit_decref(_M0L17error__to__stringS1202);
        moonbit_decref(_M0L13handle__startS1190);
        _M0L6_2atmpS2665 = 1;
      } else {
        struct moonbit_result_0 _tmp_2941;
        int32_t _handle__error__result_2942;
        #line 566 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        _tmp_2941
        = _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(_M0L12async__testsS1223, _M0L8filenameS1192, _M0L5indexS1194, _M0L13handle__startS1190, _M0L14handle__resultS1195, _M0L17error__to__stringS1202);
        if (_tmp_2941.tag) {
          int32_t const _M0L5_2aokS2668 = _tmp_2941.data.ok;
          _handle__error__result_2942 = _M0L5_2aokS2668;
        } else {
          void* const _M0L6_2aerrS2669 = _tmp_2941.data.err;
          moonbit_decref(_M0L17error__to__stringS1202);
          moonbit_decref(_M0L13handle__startS1190);
          _M0L11_2atry__errS1217 = _M0L6_2aerrS2669;
          goto join_1216;
        }
        if (_handle__error__result_2942) {
          moonbit_decref(_M0L17error__to__stringS1202);
          moonbit_decref(_M0L13handle__startS1190);
          _M0L6_2atmpS2665 = 1;
        } else {
          struct moonbit_result_0 _tmp_2943;
          #line 569 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
          _tmp_2943
          = _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(_M0L12async__testsS1223, _M0L8filenameS1192, _M0L5indexS1194, _M0L13handle__startS1190, _M0L14handle__resultS1195, _M0L17error__to__stringS1202);
          moonbit_decref(_M0L13handle__startS1190);
          moonbit_decref(_M0L17error__to__stringS1202);
          if (_tmp_2943.tag) {
            int32_t const _M0L5_2aokS2666 = _tmp_2943.data.ok;
            _M0L6_2atmpS2665 = _M0L5_2aokS2666;
          } else {
            void* const _M0L6_2aerrS2667 = _tmp_2943.data.err;
            _M0L11_2atry__errS1217 = _M0L6_2aerrS2667;
            goto join_1216;
          }
        }
      }
    }
  }
  if (!_M0L6_2atmpS2665) {
    void* _M0L125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2676 =
      (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
    Moonbit_object_header(_M0L125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2676)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 2);
    ((struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2676)->$0
    = (moonbit_string_t)moonbit_string_literal_0.data;
    _M0L11_2atry__errS1217
    = _M0L125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2676;
    goto join_1216;
  } else {
    moonbit_decref(_M0L14handle__resultS1195);
  }
  goto joinlet_2934;
  join_1216:;
  _M0L3errS1218 = _M0L11_2atry__errS1217;
  _M0L36_2aMoonBitTestDriverInternalSkipTestS1221
  = (struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L3errS1218;
  _M0L8_2afieldS2677 = _M0L36_2aMoonBitTestDriverInternalSkipTestS1221->$0;
  _M0L6_2acntS2893
  = Moonbit_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1221));
  if (_M0L6_2acntS2893 > 1) {
    int32_t _M0L11_2anew__cntS2894 = _M0L6_2acntS2893 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1221), _M0L11_2anew__cntS2894);
    moonbit_incref(_M0L8_2afieldS2677);
  } else if (_M0L6_2acntS2893 == 1) {
    #line 576 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
    moonbit_free(_M0L36_2aMoonBitTestDriverInternalSkipTestS1221);
  }
  _M0L7_2anameS1222 = _M0L8_2afieldS2677;
  _M0L4nameS1220 = _M0L7_2anameS1222;
  goto join_1219;
  goto joinlet_2944;
  join_1219:;
  #line 577 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1195(_M0L14handle__resultS1195, _M0L4nameS1220, (moonbit_string_t)moonbit_string_literal_1.data, 1);
  moonbit_decref(_M0L14handle__resultS1195);
  moonbit_decref(_M0L4nameS1220);
  joinlet_2944:;
  joinlet_2934:;
  return 0;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1202(
  struct _M0TWRPC15error5ErrorEs* _M0L6_2aenvS2664,
  void* _M0L3errS1203
) {
  void* _M0L1eS1205;
  moonbit_string_t _M0L1eS1207;
  moonbit_string_t _result_2947;
  #line 546 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  switch (Moonbit_object_tag(_M0L3errS1203)) {
    case 0: {
      struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS1208 =
        (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L3errS1203;
      moonbit_string_t _M0L4_2aeS1209 = _M0L10_2aFailureS1208->$0;
      moonbit_incref(_M0L4_2aeS1209);
      _M0L1eS1207 = _M0L4_2aeS1209;
      goto join_1206;
      break;
    }
    
    case 3: {
      struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError* _M0L15_2aInspectErrorS1210 =
        (struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError*)_M0L3errS1203;
      moonbit_string_t _M0L4_2aeS1211 = _M0L15_2aInspectErrorS1210->$0;
      moonbit_incref(_M0L4_2aeS1211);
      _M0L1eS1207 = _M0L4_2aeS1211;
      goto join_1206;
      break;
    }
    
    case 4: {
      struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError* _M0L16_2aSnapshotErrorS1212 =
        (struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError*)_M0L3errS1203;
      moonbit_string_t _M0L4_2aeS1213 = _M0L16_2aSnapshotErrorS1212->$0;
      moonbit_incref(_M0L4_2aeS1213);
      _M0L1eS1207 = _M0L4_2aeS1213;
      goto join_1206;
      break;
    }
    
    case 5: {
      struct _M0DTPC15error5Error123biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError* _M0L35_2aMoonBitTestDriverInternalJsErrorS1214 =
        (struct _M0DTPC15error5Error123biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError*)_M0L3errS1203;
      moonbit_string_t _M0L4_2aeS1215 =
        _M0L35_2aMoonBitTestDriverInternalJsErrorS1214->$0;
      moonbit_incref(_M0L4_2aeS1215);
      _M0L1eS1207 = _M0L4_2aeS1215;
      goto join_1206;
      break;
    }
    default: {
      moonbit_incref(_M0L3errS1203);
      _M0L1eS1205 = _M0L3errS1203;
      goto join_1204;
      break;
    }
  }
  join_1206:;
  return _M0L1eS1207;
  join_1204:;
  #line 552 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2947 = _M0FP15Error10to__string(_M0L1eS1205);
  moonbit_decref(_M0L1eS1205);
  return _result_2947;
}

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1195(
  struct _M0TWssbEu* _M0L6_2aenvS2661,
  moonbit_string_t _M0L10__testnameS1196,
  moonbit_string_t _M0L7messageS1197,
  int32_t _M0L7skippedS1198
) {
  struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195* _M0L14_2acasted__envS2662;
  moonbit_string_t _M0L8filenameS1192;
  int32_t _M0L5indexS1194;
  int32_t _if__result_2948;
  moonbit_string_t _M0L10file__nameS1199;
  moonbit_string_t _M0L7messageS1200;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1201;
  moonbit_string_t _M0L6_2atmpS2663;
  #line 531 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L14_2acasted__envS2662
  = (struct _M0R127_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1195*)_M0L6_2aenvS2661;
  _M0L8filenameS1192 = _M0L14_2acasted__envS2662->$1;
  _M0L5indexS1194 = _M0L14_2acasted__envS2662->$0;
  if (!_M0L7skippedS1198) {
    _if__result_2948 = 1;
  } else {
    _if__result_2948 = 0;
  }
  if (_if__result_2948) {
    
  }
  #line 537 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L10file__nameS1199
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1192, 1);
  #line 538 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L7messageS1200
  = _M0MPC16string6String14escape_2einner(_M0L7messageS1197, 1);
  #line 539 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L18_2astring__builderS1201
  = _M0MPB13StringBuilder21StringBuilder_2einner(45);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1201, (moonbit_string_t)moonbit_string_literal_3.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1201, _M0L10file__nameS1199);
  moonbit_decref(_M0L10file__nameS1199);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1201, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1201, _M0L5indexS1194);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1201, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1201, _M0L7messageS1200);
  moonbit_decref(_M0L7messageS1200);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1201, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2663
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1201);
  moonbit_decref(_M0L18_2astring__builderS1201);
  #line 540 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS2663);
  moonbit_decref(_M0L6_2atmpS2663);
  #line 543 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1190(
  struct _M0TWEOc* _M0L6_2aenvS2658
) {
  struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190* _M0L14_2acasted__envS2659;
  moonbit_string_t _M0L8filenameS1192;
  int32_t _M0L5indexS1194;
  moonbit_string_t _M0L10file__nameS1191;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1193;
  moonbit_string_t _M0L6_2atmpS2660;
  #line 522 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L14_2acasted__envS2659
  = (struct _M0R126_24biolab_2fbio__seq_2fcmd_2fseqio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1190*)_M0L6_2aenvS2658;
  _M0L8filenameS1192 = _M0L14_2acasted__envS2659->$1;
  _M0L5indexS1194 = _M0L14_2acasted__envS2659->$0;
  #line 523 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L10file__nameS1191
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1192, 1);
  #line 524 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L18_2astring__builderS1193
  = _M0MPB13StringBuilder21StringBuilder_2einner(33);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1193, (moonbit_string_t)moonbit_string_literal_8.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1193, _M0L10file__nameS1191);
  moonbit_decref(_M0L10file__nameS1191);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1193, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1193, _M0L5indexS1194);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1193, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2660
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1193);
  moonbit_decref(_M0L18_2astring__builderS1193);
  #line 525 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS2660);
  moonbit_decref(_M0L6_2atmpS2660);
  #line 528 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__args(
  
) {
  int32_t _M0L45moonbit__test__driver__internal__parse__int__S1147;
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154;
  int32_t _M0L57moonbit__test__driver__internal__get__cli__args__internalS1161;
  int32_t _M0L51moonbit__test__driver__internal__split__mbt__stringS1168;
  struct _M0TUsiE** _M0L6_2atmpS2657;
  struct _M0TPB5ArrayGUsiEE* _M0L16file__and__indexS1175;
  struct _M0TPB5ArrayGsE* _M0L9cli__argsS1176;
  moonbit_string_t _M0L6_2atmpS2656;
  struct _M0TPB5ArrayGsE* _M0L10test__argsS1177;
  int32_t _M0L7_2abindS1178;
  int32_t _M0L2__S1179;
  #line 194 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L45moonbit__test__driver__internal__parse__int__S1147 = 0;
  _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154
  = 0;
  _M0L57moonbit__test__driver__internal__get__cli__args__internalS1161
  = _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154;
  _M0L51moonbit__test__driver__internal__split__mbt__stringS1168 = 0;
  _M0L6_2atmpS2657 = (struct _M0TUsiE**)moonbit_empty_ref_array;
  _M0L16file__and__indexS1175
  = (struct _M0TPB5ArrayGUsiEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGUsiEE));
  Moonbit_object_header(_M0L16file__and__indexS1175)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 9, 0);
  _M0L16file__and__indexS1175->$0 = _M0L6_2atmpS2657;
  _M0L16file__and__indexS1175->$1 = 0;
  #line 283 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L9cli__argsS1176
  = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1161(_M0L57moonbit__test__driver__internal__get__cli__args__internalS1161);
  #line 285 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2656 = _M0MPC15array5Array2atGsE(_M0L9cli__argsS1176, 1);
  moonbit_decref(_M0L9cli__argsS1176);
  #line 284 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L10test__argsS1177
  = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1168(_M0L51moonbit__test__driver__internal__split__mbt__stringS1168, _M0L6_2atmpS2656, 47);
  moonbit_decref(_M0L6_2atmpS2656);
  _M0L7_2abindS1178 = _M0L10test__argsS1177->$1;
  _M0L2__S1179 = 0;
  while (1) {
    if (_M0L2__S1179 < _M0L7_2abindS1178) {
      moonbit_string_t* _M0L3bufS2655 = _M0L10test__argsS1177->$0;
      moonbit_string_t _M0L3argS1180 =
        (moonbit_string_t)_M0L3bufS2655[_M0L2__S1179];
      struct _M0TPB5ArrayGsE* _M0L16file__and__rangeS1181;
      moonbit_string_t _M0L4fileS1182;
      moonbit_string_t _M0L5rangeS1183;
      struct _M0TPB5ArrayGsE* _M0L15start__and__endS1184;
      moonbit_string_t _M0L6_2atmpS2653;
      int32_t _M0L5startS1185;
      moonbit_string_t _M0L6_2atmpS2652;
      int32_t _M0L3endS1186;
      int32_t _M0L1iS1187;
      int32_t _M0L6_2atmpS2654;
      moonbit_incref(_M0L3argS1180);
      #line 289 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L16file__and__rangeS1181
      = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1168(_M0L51moonbit__test__driver__internal__split__mbt__stringS1168, _M0L3argS1180, 58);
      moonbit_decref(_M0L3argS1180);
      #line 290 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L4fileS1182
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1181, 0);
      #line 291 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L5rangeS1183
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1181, 1);
      moonbit_decref(_M0L16file__and__rangeS1181);
      #line 292 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L15start__and__endS1184
      = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1168(_M0L51moonbit__test__driver__internal__split__mbt__stringS1168, _M0L5rangeS1183, 45);
      moonbit_decref(_M0L5rangeS1183);
      #line 295 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS2653
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1184, 0);
      #line 295 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L5startS1185
      = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1147(_M0L45moonbit__test__driver__internal__parse__int__S1147, _M0L6_2atmpS2653);
      moonbit_decref(_M0L6_2atmpS2653);
      #line 296 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS2652
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1184, 1);
      moonbit_decref(_M0L15start__and__endS1184);
      #line 296 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L3endS1186
      = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1147(_M0L45moonbit__test__driver__internal__parse__int__S1147, _M0L6_2atmpS2652);
      moonbit_decref(_M0L6_2atmpS2652);
      _M0L1iS1187 = _M0L5startS1185;
      while (1) {
        if (_M0L1iS1187 < _M0L3endS1186) {
          struct _M0TUsiE* _M0L8_2atupleS2650;
          int32_t _M0L6_2atmpS2651;
          moonbit_incref(_M0L4fileS1182);
          _M0L8_2atupleS2650
          = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
          Moonbit_object_header(_M0L8_2atupleS2650)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 12, 0);
          _M0L8_2atupleS2650->$0 = _M0L4fileS1182;
          _M0L8_2atupleS2650->$1 = _M0L1iS1187;
          #line 298 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
          _M0MPC15array5Array4pushGUsiEE(_M0L16file__and__indexS1175, _M0L8_2atupleS2650);
          moonbit_decref(_M0L8_2atupleS2650);
          _M0L6_2atmpS2651 = _M0L1iS1187 + 1;
          _M0L1iS1187 = _M0L6_2atmpS2651;
          continue;
        } else {
          moonbit_decref(_M0L4fileS1182);
        }
        break;
      }
      _M0L6_2atmpS2654 = _M0L2__S1179 + 1;
      _M0L2__S1179 = _M0L6_2atmpS2654;
      continue;
    } else {
      moonbit_decref(_M0L10test__argsS1177);
    }
    break;
  }
  return _M0L16file__and__indexS1175;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1168(
  int32_t _M0L6_2aenvS2631,
  moonbit_string_t _M0L1sS1169,
  int32_t _M0L3sepS1170
) {
  moonbit_string_t* _M0L6_2atmpS2649;
  struct _M0TPB5ArrayGsE* _M0L3resS1171;
  struct _M0TPB8MutLocalGiE* _M0L1iS1172;
  struct _M0TPB8MutLocalGiE* _M0L5startS1173;
  int32_t _M0L3valS2644;
  int32_t _M0L6_2atmpS2645;
  #line 262 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2649 = (moonbit_string_t*)moonbit_empty_ref_array;
  _M0L3resS1171
  = (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
  Moonbit_object_header(_M0L3resS1171)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
  _M0L3resS1171->$0 = _M0L6_2atmpS2649;
  _M0L3resS1171->$1 = 0;
  _M0L1iS1172
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1172)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1172->$0 = 0;
  _M0L5startS1173
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS1173)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS1173->$0 = 0;
  while (1) {
    int32_t _M0L3valS2632 = _M0L1iS1172->$0;
    int32_t _M0L6_2atmpS2633 = Moonbit_array_length(_M0L1sS1169);
    if (_M0L3valS2632 < _M0L6_2atmpS2633) {
      int32_t _M0L3valS2636 = _M0L1iS1172->$0;
      int32_t _M0L6_2atmpS2635;
      int32_t _M0L6_2atmpS2634;
      int32_t _M0L3valS2643;
      int32_t _M0L6_2atmpS2642;
      if (
        _M0L3valS2636 < 0
        || _M0L3valS2636 >= Moonbit_array_length(_M0L1sS1169)
      ) {
        #line 270 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2635 = _M0L1sS1169[_M0L3valS2636];
      _M0L6_2atmpS2634 = _M0L6_2atmpS2635;
      if (_M0L6_2atmpS2634 == _M0L3sepS1170) {
        int32_t _M0L3valS2638 = _M0L5startS1173->$0;
        int32_t _M0L3valS2639 = _M0L1iS1172->$0;
        moonbit_string_t _M0L6_2atmpS2637;
        int32_t _M0L3valS2641;
        int32_t _M0L6_2atmpS2640;
        #line 271 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        _M0L6_2atmpS2637
        = _M0MPC16string6String17unsafe__substring(_M0L1sS1169, _M0L3valS2638, _M0L3valS2639);
        #line 271 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        _M0MPC15array5Array4pushGsE(_M0L3resS1171, _M0L6_2atmpS2637);
        moonbit_decref(_M0L6_2atmpS2637);
        _M0L3valS2641 = _M0L1iS1172->$0;
        _M0L6_2atmpS2640 = _M0L3valS2641 + 1;
        _M0L5startS1173->$0 = _M0L6_2atmpS2640;
      }
      _M0L3valS2643 = _M0L1iS1172->$0;
      _M0L6_2atmpS2642 = _M0L3valS2643 + 1;
      _M0L1iS1172->$0 = _M0L6_2atmpS2642;
      continue;
    } else {
      moonbit_decref(_M0L1iS1172);
    }
    break;
  }
  _M0L3valS2644 = _M0L5startS1173->$0;
  _M0L6_2atmpS2645 = Moonbit_array_length(_M0L1sS1169);
  if (_M0L3valS2644 < _M0L6_2atmpS2645) {
    int32_t _M0L3valS2647 = _M0L5startS1173->$0;
    int32_t _M0L6_2atmpS2648;
    moonbit_string_t _M0L6_2atmpS2646;
    moonbit_decref(_M0L5startS1173);
    _M0L6_2atmpS2648 = Moonbit_array_length(_M0L1sS1169);
    #line 277 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
    _M0L6_2atmpS2646
    = _M0MPC16string6String17unsafe__substring(_M0L1sS1169, _M0L3valS2647, _M0L6_2atmpS2648);
    #line 277 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
    _M0MPC15array5Array4pushGsE(_M0L3resS1171, _M0L6_2atmpS2646);
    moonbit_decref(_M0L6_2atmpS2646);
  } else {
    moonbit_decref(_M0L5startS1173);
  }
  return _M0L3resS1171;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1161(
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154
) {
  moonbit_bytes_t* _M0L3tmpS1162;
  int32_t _M0L6_2atmpS2630;
  struct _M0TPB5ArrayGsE* _M0L3resS1163;
  int32_t _M0L7_2abindS1164;
  int32_t _M0L7_2abindS1165;
  int32_t _M0L1iS1166;
  #line 251 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  #line 254 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3tmpS1162
  = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__get__cli__args__ffi();
  _M0L6_2atmpS2630 = Moonbit_array_length(_M0L3tmpS1162);
  #line 255 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1163 = _M0MPC15array5Array11new_2einnerGsE(_M0L6_2atmpS2630);
  _M0L7_2abindS1164 = 0;
  _M0L7_2abindS1165 = Moonbit_array_length(_M0L3tmpS1162);
  _M0L1iS1166 = _M0L7_2abindS1164;
  while (1) {
    if (_M0L1iS1166 < _M0L7_2abindS1165) {
      moonbit_bytes_t _M0L6_2atmpS2628;
      moonbit_string_t _M0L6_2atmpS2627;
      int32_t _M0L6_2atmpS2629;
      if (
        _M0L1iS1166 < 0 || _M0L1iS1166 >= Moonbit_array_length(_M0L3tmpS1162)
      ) {
        #line 257 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2628 = (moonbit_bytes_t)_M0L3tmpS1162[_M0L1iS1166];
      moonbit_incref(_M0L6_2atmpS2628);
      #line 257 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS2627
      = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154(_M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154, _M0L6_2atmpS2628);
      moonbit_decref(_M0L6_2atmpS2628);
      #line 257 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0MPC15array5Array4pushGsE(_M0L3resS1163, _M0L6_2atmpS2627);
      moonbit_decref(_M0L6_2atmpS2627);
      _M0L6_2atmpS2629 = _M0L1iS1166 + 1;
      _M0L1iS1166 = _M0L6_2atmpS2629;
      continue;
    } else {
      moonbit_decref(_M0L3tmpS1162);
    }
    break;
  }
  return _M0L3resS1163;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1154(
  int32_t _M0L6_2aenvS2541,
  moonbit_bytes_t _M0L5bytesS1155
) {
  struct _M0TPB13StringBuilder* _M0L3resS1156;
  int32_t _M0L3lenS1157;
  struct _M0TPB8MutLocalGiE* _M0L1iS1158;
  moonbit_string_t _result_2954;
  #line 207 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  #line 210 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1156 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3lenS1157 = Moonbit_array_length(_M0L5bytesS1155);
  _M0L1iS1158
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1158)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1158->$0 = 0;
  while (1) {
    int32_t _M0L3valS2542 = _M0L1iS1158->$0;
    if (_M0L3valS2542 < _M0L3lenS1157) {
      int32_t _M0L3valS2626 = _M0L1iS1158->$0;
      int32_t _M0L6_2atmpS2625;
      int32_t _M0L6_2atmpS2624;
      struct _M0TPB8MutLocalGiE* _M0L1cS1159;
      int32_t _M0L3valS2543;
      if (
        _M0L3valS2626 < 0
        || _M0L3valS2626 >= Moonbit_array_length(_M0L5bytesS1155)
      ) {
        #line 214 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2625 = _M0L5bytesS1155[_M0L3valS2626];
      _M0L6_2atmpS2624 = (int32_t)_M0L6_2atmpS2625;
      _M0L1cS1159
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1cS1159)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1cS1159->$0 = _M0L6_2atmpS2624;
      _M0L3valS2543 = _M0L1cS1159->$0;
      if (_M0L3valS2543 < 128) {
        int32_t _M0L3valS2545 = _M0L1cS1159->$0;
        int32_t _M0L6_2atmpS2544;
        int32_t _M0L3valS2547;
        int32_t _M0L6_2atmpS2546;
        moonbit_decref(_M0L1cS1159);
        _M0L6_2atmpS2544 = _M0L3valS2545;
        #line 216 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1156, _M0L6_2atmpS2544);
        _M0L3valS2547 = _M0L1iS1158->$0;
        _M0L6_2atmpS2546 = _M0L3valS2547 + 1;
        _M0L1iS1158->$0 = _M0L6_2atmpS2546;
      } else {
        int32_t _M0L3valS2548 = _M0L1cS1159->$0;
        if (_M0L3valS2548 < 224) {
          int32_t _M0L3valS2550 = _M0L1iS1158->$0;
          int32_t _M0L6_2atmpS2549 = _M0L3valS2550 + 1;
          int32_t _M0L3valS2559;
          int32_t _M0L6_2atmpS2558;
          int32_t _M0L6_2atmpS2552;
          int32_t _M0L3valS2557;
          int32_t _M0L6_2atmpS2556;
          int32_t _M0L6_2atmpS2555;
          int32_t _M0L6_2atmpS2554;
          int32_t _M0L6_2atmpS2553;
          int32_t _M0L6_2atmpS2551;
          int32_t _M0L3valS2561;
          int32_t _M0L6_2atmpS2560;
          int32_t _M0L3valS2563;
          int32_t _M0L6_2atmpS2562;
          if (_M0L6_2atmpS2549 >= _M0L3lenS1157) {
            moonbit_decref(_M0L1cS1159);
            moonbit_decref(_M0L1iS1158);
            break;
          }
          _M0L3valS2559 = _M0L1cS1159->$0;
          _M0L6_2atmpS2558 = _M0L3valS2559 & 31;
          _M0L6_2atmpS2552 = _M0L6_2atmpS2558 << 6;
          _M0L3valS2557 = _M0L1iS1158->$0;
          _M0L6_2atmpS2556 = _M0L3valS2557 + 1;
          if (
            _M0L6_2atmpS2556 < 0
            || _M0L6_2atmpS2556 >= Moonbit_array_length(_M0L5bytesS1155)
          ) {
            #line 222 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2555 = _M0L5bytesS1155[_M0L6_2atmpS2556];
          _M0L6_2atmpS2554 = (int32_t)_M0L6_2atmpS2555;
          _M0L6_2atmpS2553 = _M0L6_2atmpS2554 & 63;
          _M0L6_2atmpS2551 = _M0L6_2atmpS2552 | _M0L6_2atmpS2553;
          _M0L1cS1159->$0 = _M0L6_2atmpS2551;
          _M0L3valS2561 = _M0L1cS1159->$0;
          moonbit_decref(_M0L1cS1159);
          _M0L6_2atmpS2560 = _M0L3valS2561;
          #line 223 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1156, _M0L6_2atmpS2560);
          _M0L3valS2563 = _M0L1iS1158->$0;
          _M0L6_2atmpS2562 = _M0L3valS2563 + 2;
          _M0L1iS1158->$0 = _M0L6_2atmpS2562;
        } else {
          int32_t _M0L3valS2564 = _M0L1cS1159->$0;
          if (_M0L3valS2564 < 240) {
            int32_t _M0L3valS2566 = _M0L1iS1158->$0;
            int32_t _M0L6_2atmpS2565 = _M0L3valS2566 + 2;
            int32_t _M0L3valS2582;
            int32_t _M0L6_2atmpS2581;
            int32_t _M0L6_2atmpS2574;
            int32_t _M0L3valS2580;
            int32_t _M0L6_2atmpS2579;
            int32_t _M0L6_2atmpS2578;
            int32_t _M0L6_2atmpS2577;
            int32_t _M0L6_2atmpS2576;
            int32_t _M0L6_2atmpS2575;
            int32_t _M0L6_2atmpS2568;
            int32_t _M0L3valS2573;
            int32_t _M0L6_2atmpS2572;
            int32_t _M0L6_2atmpS2571;
            int32_t _M0L6_2atmpS2570;
            int32_t _M0L6_2atmpS2569;
            int32_t _M0L6_2atmpS2567;
            int32_t _M0L3valS2584;
            int32_t _M0L6_2atmpS2583;
            int32_t _M0L3valS2586;
            int32_t _M0L6_2atmpS2585;
            if (_M0L6_2atmpS2565 >= _M0L3lenS1157) {
              moonbit_decref(_M0L1cS1159);
              moonbit_decref(_M0L1iS1158);
              break;
            }
            _M0L3valS2582 = _M0L1cS1159->$0;
            _M0L6_2atmpS2581 = _M0L3valS2582 & 15;
            _M0L6_2atmpS2574 = _M0L6_2atmpS2581 << 12;
            _M0L3valS2580 = _M0L1iS1158->$0;
            _M0L6_2atmpS2579 = _M0L3valS2580 + 1;
            if (
              _M0L6_2atmpS2579 < 0
              || _M0L6_2atmpS2579 >= Moonbit_array_length(_M0L5bytesS1155)
            ) {
              #line 230 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2578 = _M0L5bytesS1155[_M0L6_2atmpS2579];
            _M0L6_2atmpS2577 = (int32_t)_M0L6_2atmpS2578;
            _M0L6_2atmpS2576 = _M0L6_2atmpS2577 & 63;
            _M0L6_2atmpS2575 = _M0L6_2atmpS2576 << 6;
            _M0L6_2atmpS2568 = _M0L6_2atmpS2574 | _M0L6_2atmpS2575;
            _M0L3valS2573 = _M0L1iS1158->$0;
            _M0L6_2atmpS2572 = _M0L3valS2573 + 2;
            if (
              _M0L6_2atmpS2572 < 0
              || _M0L6_2atmpS2572 >= Moonbit_array_length(_M0L5bytesS1155)
            ) {
              #line 231 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2571 = _M0L5bytesS1155[_M0L6_2atmpS2572];
            _M0L6_2atmpS2570 = (int32_t)_M0L6_2atmpS2571;
            _M0L6_2atmpS2569 = _M0L6_2atmpS2570 & 63;
            _M0L6_2atmpS2567 = _M0L6_2atmpS2568 | _M0L6_2atmpS2569;
            _M0L1cS1159->$0 = _M0L6_2atmpS2567;
            _M0L3valS2584 = _M0L1cS1159->$0;
            moonbit_decref(_M0L1cS1159);
            _M0L6_2atmpS2583 = _M0L3valS2584;
            #line 232 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1156, _M0L6_2atmpS2583);
            _M0L3valS2586 = _M0L1iS1158->$0;
            _M0L6_2atmpS2585 = _M0L3valS2586 + 3;
            _M0L1iS1158->$0 = _M0L6_2atmpS2585;
          } else {
            int32_t _M0L3valS2588 = _M0L1iS1158->$0;
            int32_t _M0L6_2atmpS2587 = _M0L3valS2588 + 3;
            int32_t _M0L3valS2611;
            int32_t _M0L6_2atmpS2610;
            int32_t _M0L6_2atmpS2603;
            int32_t _M0L3valS2609;
            int32_t _M0L6_2atmpS2608;
            int32_t _M0L6_2atmpS2607;
            int32_t _M0L6_2atmpS2606;
            int32_t _M0L6_2atmpS2605;
            int32_t _M0L6_2atmpS2604;
            int32_t _M0L6_2atmpS2596;
            int32_t _M0L3valS2602;
            int32_t _M0L6_2atmpS2601;
            int32_t _M0L6_2atmpS2600;
            int32_t _M0L6_2atmpS2599;
            int32_t _M0L6_2atmpS2598;
            int32_t _M0L6_2atmpS2597;
            int32_t _M0L6_2atmpS2590;
            int32_t _M0L3valS2595;
            int32_t _M0L6_2atmpS2594;
            int32_t _M0L6_2atmpS2593;
            int32_t _M0L6_2atmpS2592;
            int32_t _M0L6_2atmpS2591;
            int32_t _M0L6_2atmpS2589;
            int32_t _M0L3valS2613;
            int32_t _M0L6_2atmpS2612;
            int32_t _M0L3valS2617;
            int32_t _M0L6_2atmpS2616;
            int32_t _M0L6_2atmpS2615;
            int32_t _M0L6_2atmpS2614;
            int32_t _M0L3valS2621;
            int32_t _M0L6_2atmpS2620;
            int32_t _M0L6_2atmpS2619;
            int32_t _M0L6_2atmpS2618;
            int32_t _M0L3valS2623;
            int32_t _M0L6_2atmpS2622;
            if (_M0L6_2atmpS2587 >= _M0L3lenS1157) {
              moonbit_decref(_M0L1cS1159);
              moonbit_decref(_M0L1iS1158);
              break;
            }
            _M0L3valS2611 = _M0L1cS1159->$0;
            _M0L6_2atmpS2610 = _M0L3valS2611 & 7;
            _M0L6_2atmpS2603 = _M0L6_2atmpS2610 << 18;
            _M0L3valS2609 = _M0L1iS1158->$0;
            _M0L6_2atmpS2608 = _M0L3valS2609 + 1;
            if (
              _M0L6_2atmpS2608 < 0
              || _M0L6_2atmpS2608 >= Moonbit_array_length(_M0L5bytesS1155)
            ) {
              #line 239 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2607 = _M0L5bytesS1155[_M0L6_2atmpS2608];
            _M0L6_2atmpS2606 = (int32_t)_M0L6_2atmpS2607;
            _M0L6_2atmpS2605 = _M0L6_2atmpS2606 & 63;
            _M0L6_2atmpS2604 = _M0L6_2atmpS2605 << 12;
            _M0L6_2atmpS2596 = _M0L6_2atmpS2603 | _M0L6_2atmpS2604;
            _M0L3valS2602 = _M0L1iS1158->$0;
            _M0L6_2atmpS2601 = _M0L3valS2602 + 2;
            if (
              _M0L6_2atmpS2601 < 0
              || _M0L6_2atmpS2601 >= Moonbit_array_length(_M0L5bytesS1155)
            ) {
              #line 240 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2600 = _M0L5bytesS1155[_M0L6_2atmpS2601];
            _M0L6_2atmpS2599 = (int32_t)_M0L6_2atmpS2600;
            _M0L6_2atmpS2598 = _M0L6_2atmpS2599 & 63;
            _M0L6_2atmpS2597 = _M0L6_2atmpS2598 << 6;
            _M0L6_2atmpS2590 = _M0L6_2atmpS2596 | _M0L6_2atmpS2597;
            _M0L3valS2595 = _M0L1iS1158->$0;
            _M0L6_2atmpS2594 = _M0L3valS2595 + 3;
            if (
              _M0L6_2atmpS2594 < 0
              || _M0L6_2atmpS2594 >= Moonbit_array_length(_M0L5bytesS1155)
            ) {
              #line 241 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2593 = _M0L5bytesS1155[_M0L6_2atmpS2594];
            _M0L6_2atmpS2592 = (int32_t)_M0L6_2atmpS2593;
            _M0L6_2atmpS2591 = _M0L6_2atmpS2592 & 63;
            _M0L6_2atmpS2589 = _M0L6_2atmpS2590 | _M0L6_2atmpS2591;
            _M0L1cS1159->$0 = _M0L6_2atmpS2589;
            _M0L3valS2613 = _M0L1cS1159->$0;
            _M0L6_2atmpS2612 = _M0L3valS2613 - 65536;
            _M0L1cS1159->$0 = _M0L6_2atmpS2612;
            _M0L3valS2617 = _M0L1cS1159->$0;
            _M0L6_2atmpS2616 = _M0L3valS2617 >> 10;
            _M0L6_2atmpS2615 = _M0L6_2atmpS2616 + 55296;
            _M0L6_2atmpS2614 = _M0L6_2atmpS2615;
            #line 243 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1156, _M0L6_2atmpS2614);
            _M0L3valS2621 = _M0L1cS1159->$0;
            moonbit_decref(_M0L1cS1159);
            _M0L6_2atmpS2620 = _M0L3valS2621 & 1023;
            _M0L6_2atmpS2619 = _M0L6_2atmpS2620 + 56320;
            _M0L6_2atmpS2618 = _M0L6_2atmpS2619;
            #line 244 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1156, _M0L6_2atmpS2618);
            _M0L3valS2623 = _M0L1iS1158->$0;
            _M0L6_2atmpS2622 = _M0L3valS2623 + 4;
            _M0L1iS1158->$0 = _M0L6_2atmpS2622;
          }
        }
      }
      continue;
    } else {
      moonbit_decref(_M0L1iS1158);
    }
    break;
  }
  #line 248 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2954 = _M0MPB13StringBuilder10to__string(_M0L3resS1156);
  moonbit_decref(_M0L3resS1156);
  return _result_2954;
}

int32_t _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1147(
  int32_t _M0L6_2aenvS2534,
  moonbit_string_t _M0L1sS1148
) {
  struct _M0TPB8MutLocalGiE* _M0L3resS1149;
  int32_t _M0L3lenS1150;
  int32_t _M0L7_2abindS1151;
  int32_t _M0L1iS1152;
  int32_t _result_2956;
  #line 198 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1149
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3resS1149)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3resS1149->$0 = 0;
  _M0L3lenS1150 = Moonbit_array_length(_M0L1sS1148);
  _M0L7_2abindS1151 = 0;
  _M0L1iS1152 = _M0L7_2abindS1151;
  while (1) {
    if (_M0L1iS1152 < _M0L3lenS1150) {
      int32_t _M0L3valS2539 = _M0L3resS1149->$0;
      int32_t _M0L6_2atmpS2536 = _M0L3valS2539 * 10;
      int32_t _M0L6_2atmpS2538;
      int32_t _M0L6_2atmpS2537;
      int32_t _M0L6_2atmpS2535;
      int32_t _M0L6_2atmpS2540;
      if (
        _M0L1iS1152 < 0 || _M0L1iS1152 >= Moonbit_array_length(_M0L1sS1148)
      ) {
        #line 202 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2538 = _M0L1sS1148[_M0L1iS1152];
      _M0L6_2atmpS2537 = _M0L6_2atmpS2538 - 48;
      _M0L6_2atmpS2535 = _M0L6_2atmpS2536 + _M0L6_2atmpS2537;
      _M0L3resS1149->$0 = _M0L6_2atmpS2535;
      _M0L6_2atmpS2540 = _M0L1iS1152 + 1;
      _M0L1iS1152 = _M0L6_2atmpS2540;
      continue;
    }
    break;
  }
  _result_2956 = _M0L3resS1149->$0;
  moonbit_decref(_M0L3resS1149);
  return _result_2956;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1117,
  moonbit_string_t _M0L12_2adiscard__S1118,
  int32_t _M0L12_2adiscard__S1119,
  struct _M0TWEOc* _M0L12_2adiscard__S1120,
  struct _M0TWssbEu* _M0L12_2adiscard__S1121,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1122
) {
  struct moonbit_result_0 _result_2957;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2957.tag = 1;
  _result_2957.data.ok = 0;
  return _result_2957;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1123,
  moonbit_string_t _M0L12_2adiscard__S1124,
  int32_t _M0L12_2adiscard__S1125,
  struct _M0TWEOc* _M0L12_2adiscard__S1126,
  struct _M0TWssbEu* _M0L12_2adiscard__S1127,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1128
) {
  struct moonbit_result_0 _result_2958;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2958.tag = 1;
  _result_2958.data.ok = 0;
  return _result_2958;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1129,
  moonbit_string_t _M0L12_2adiscard__S1130,
  int32_t _M0L12_2adiscard__S1131,
  struct _M0TWEOc* _M0L12_2adiscard__S1132,
  struct _M0TWssbEu* _M0L12_2adiscard__S1133,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1134
) {
  struct moonbit_result_0 _result_2959;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2959.tag = 1;
  _result_2959.data.ok = 0;
  return _result_2959;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1135,
  moonbit_string_t _M0L12_2adiscard__S1136,
  int32_t _M0L12_2adiscard__S1137,
  struct _M0TWEOc* _M0L12_2adiscard__S1138,
  struct _M0TWssbEu* _M0L12_2adiscard__S1139,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1140
) {
  struct moonbit_result_0 _result_2960;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2960.tag = 1;
  _result_2960.data.ok = 0;
  return _result_2960;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1141,
  moonbit_string_t _M0L12_2adiscard__S1142,
  int32_t _M0L12_2adiscard__S1143,
  struct _M0TWEOc* _M0L12_2adiscard__S1144,
  struct _M0TWssbEu* _M0L12_2adiscard__S1145,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1146
) {
  struct moonbit_result_0 _result_2961;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _result_2961.tag = 1;
  _result_2961.data.ok = 0;
  return _result_2961;
}

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1116
) {
  #line 12 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main11seqio__test() {
  moonbit_string_t _M0L12fasta__multiS1064;
  moonbit_string_t _M0L16fasta__lowercaseS1065;
  moonbit_string_t _M0L16fasta__multilineS1066;
  moonbit_string_t _M0L13fasta__singleS1067;
  moonbit_string_t _M0L6_2atmpS2533;
  moonbit_string_t _M0L6_2atmpS2532;
  moonbit_string_t _M0L12fastq__multiS1068;
  void* _M0L11_2atry__errS1071;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4recsS1069;
  struct moonbit_result_1 _tmp_2963;
  moonbit_string_t _M0L6_2atmpS2480;
  void* _M0L11_2atry__errS1074;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recs__lS1072;
  struct moonbit_result_1 _tmp_2965;
  moonbit_string_t _M0L6_2atmpS2481;
  void* _M0L11_2atry__errS1077;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recs__mS1075;
  struct moonbit_result_1 _tmp_2967;
  moonbit_string_t _M0L6_2atmpS2482;
  void* _M0L11_2atry__errS1080;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6singleS1078;
  struct moonbit_result_2 _tmp_2969;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2523;
  moonbit_string_t _M0L6_2atmpS2483;
  void* _M0L11_2atry__errS1083;
  moonbit_string_t _M0L12err__payloadS1081;
  struct moonbit_result_2 _tmp_2971;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2520;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L1dS1084;
  struct _M0TPB13StringBuilder* _M0L4dbufS1085;
  int32_t _M0L7_2abindS1086;
  int32_t _M0L1iS1087;
  moonbit_string_t _M0L6_2atmpS2490;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2519;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2514;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2518;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2515;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2517;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2516;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2513;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L11write__recsS1096;
  void* _M0L11_2atry__errS1099;
  moonbit_string_t _M0L10fasta__outS1097;
  struct moonbit_result_3 _tmp_2975;
  moonbit_string_t _M0L6_2atmpS2491;
  void* _M0L11_2atry__errS1102;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L8fq__recsS1100;
  struct moonbit_result_1 _tmp_2977;
  moonbit_string_t _M0L6_2atmpS2492;
  struct _M0TPB5ArrayGiE* _M0L10fq1__qualsS1103;
  int32_t _M0L1iS1104;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2508;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2507;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2501;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2506;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2503;
  int32_t* _M0L6_2atmpS2505;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS2504;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2502;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2500;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L15write__fq__recsS1106;
  void* _M0L11_2atry__errS1109;
  moonbit_string_t _M0L7fq__outS1107;
  struct moonbit_result_3 _tmp_2980;
  moonbit_string_t _M0L6_2atmpS2494;
  moonbit_string_t _M0L15genbank__sampleS1110;
  void* _M0L11_2atry__errS1113;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L8gb__recsS1111;
  struct moonbit_result_1 _tmp_2982;
  int32_t _M0L7_2abindS1115;
  moonbit_string_t _M0L8gb__jsonS1114;
  #line 93 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L12fasta__multiS1064 = (moonbit_string_t)moonbit_string_literal_9.data;
  _M0L16fasta__lowercaseS1065
  = (moonbit_string_t)moonbit_string_literal_10.data;
  _M0L16fasta__multilineS1066
  = (moonbit_string_t)moonbit_string_literal_11.data;
  _M0L13fasta__singleS1067 = (moonbit_string_t)moonbit_string_literal_12.data;
  #line 102 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2533
  = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_13.data, 39);
  #line 101 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2532
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_14.data, _M0L6_2atmpS2533);
  moonbit_decref(_M0L6_2atmpS2533);
  #line 101 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L12fastq__multiS1068
  = moonbit_add_string(_M0L6_2atmpS2532, (moonbit_string_t)moonbit_string_literal_15.data);
  moonbit_decref(_M0L6_2atmpS2532);
  #line 106 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2963
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L12fasta__multiS1064, (moonbit_string_t)moonbit_string_literal_16.data);
  if (_tmp_2963.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2530 =
      _tmp_2963.data.ok;
    _M0L4recsS1069 = _M0L5_2aokS2530;
  } else {
    void* const _M0L6_2aerrS2531 = _tmp_2963.data.err;
    _M0L11_2atry__errS1071 = _M0L6_2aerrS2531;
    goto join_1070;
  }
  goto joinlet_2962;
  join_1070:;
  moonbit_decref(_M0L11_2atry__errS1071);
  #line 106 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L4recsS1069
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  joinlet_2962:;
  #line 107 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2480
  = _M0FP46biolab8bio__seq3cmd11seqio__main25records__to__json_2einner(_M0L4recsS1069, 0);
  #line 107 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_17.data, _M0L6_2atmpS2480);
  moonbit_decref(_M0L6_2atmpS2480);
  #line 110 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2965
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L16fasta__lowercaseS1065, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L16fasta__lowercaseS1065);
  if (_tmp_2965.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2528 =
      _tmp_2965.data.ok;
    _M0L7recs__lS1072 = _M0L5_2aokS2528;
  } else {
    void* const _M0L6_2aerrS2529 = _tmp_2965.data.err;
    _M0L11_2atry__errS1074 = _M0L6_2aerrS2529;
    goto join_1073;
  }
  goto joinlet_2964;
  join_1073:;
  moonbit_decref(_M0L11_2atry__errS1074);
  #line 110 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L7recs__lS1072
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  joinlet_2964:;
  #line 113 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2481
  = _M0FP46biolab8bio__seq3cmd11seqio__main25records__to__json_2einner(_M0L7recs__lS1072, 0);
  moonbit_decref(_M0L7recs__lS1072);
  #line 113 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_18.data, _M0L6_2atmpS2481);
  moonbit_decref(_M0L6_2atmpS2481);
  #line 116 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2967
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L16fasta__multilineS1066, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L16fasta__multilineS1066);
  if (_tmp_2967.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2526 =
      _tmp_2967.data.ok;
    _M0L7recs__mS1075 = _M0L5_2aokS2526;
  } else {
    void* const _M0L6_2aerrS2527 = _tmp_2967.data.err;
    _M0L11_2atry__errS1077 = _M0L6_2aerrS2527;
    goto join_1076;
  }
  goto joinlet_2966;
  join_1076:;
  moonbit_decref(_M0L11_2atry__errS1077);
  #line 116 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L7recs__mS1075
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  joinlet_2966:;
  #line 119 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2482
  = _M0FP46biolab8bio__seq3cmd11seqio__main25records__to__json_2einner(_M0L7recs__mS1075, 0);
  moonbit_decref(_M0L7recs__mS1075);
  #line 119 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_19.data, _M0L6_2atmpS2482);
  moonbit_decref(_M0L6_2atmpS2482);
  #line 122 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2969
  = _M0FP26biolab8bio__seq11seqio__read(_M0L13fasta__singleS1067, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L13fasta__singleS1067);
  if (_tmp_2969.tag) {
    struct _M0TP26biolab8bio__seq9SeqRecord* const _M0L5_2aokS2524 =
      _tmp_2969.data.ok;
    _M0L6singleS1078 = _M0L5_2aokS2524;
  } else {
    void* const _M0L6_2aerrS2525 = _tmp_2969.data.err;
    _M0L11_2atry__errS1080 = _M0L6_2aerrS2525;
    goto join_1079;
  }
  goto joinlet_2968;
  join_1079:;
  moonbit_decref(_M0L11_2atry__errS1080);
  #line 123 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2523
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_0.data);
  #line 122 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6singleS1078
  = _M0MP26biolab8bio__seq9SeqRecord9from__seq(_M0L6_2atmpS2523);
  moonbit_decref(_M0L6_2atmpS2523);
  joinlet_2968:;
  #line 125 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2483
  = _M0FP46biolab8bio__seq3cmd11seqio__main24record__to__json_2einner(_M0L6singleS1078, 0);
  moonbit_decref(_M0L6singleS1078);
  #line 125 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_20.data, _M0L6_2atmpS2483);
  moonbit_decref(_M0L6_2atmpS2483);
  #line 129 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2971
  = _M0FP26biolab8bio__seq11seqio__read(_M0L12fasta__multiS1064, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L12fasta__multiS1064);
  if (_tmp_2971.tag) {
    struct _M0TP26biolab8bio__seq9SeqRecord* const _M0L5_2aokS2521 =
      _tmp_2971.data.ok;
    _M0L6_2atmpS2520 = _M0L5_2aokS2521;
  } else {
    void* const _M0L6_2aerrS2522 = _tmp_2971.data.err;
    _M0L11_2atry__errS1083 = _M0L6_2aerrS2522;
    goto join_1082;
  }
  moonbit_decref(_M0L6_2atmpS2520);
  _M0L12err__payloadS1081 = (moonbit_string_t)moonbit_string_literal_21.data;
  goto joinlet_2970;
  join_1082:;
  moonbit_decref(_M0L11_2atry__errS1083);
  _M0L12err__payloadS1081 = (moonbit_string_t)moonbit_string_literal_22.data;
  joinlet_2970:;
  #line 135 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_23.data, _M0L12err__payloadS1081);
  moonbit_decref(_M0L12err__payloadS1081);
  #line 139 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L1dS1084 = _M0FP26biolab8bio__seq15seqio__to__dict(_M0L4recsS1069);
  #line 140 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L4dbufS1085 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 141 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L4dbufS1085, 123);
  _M0L7_2abindS1086 = _M0L4recsS1069->$1;
  _M0L1iS1087 = 0;
  while (1) {
    if (_M0L1iS1087 < _M0L7_2abindS1086) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2489 =
        _M0L4recsS1069->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1088 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2489[_M0L1iS1087];
      moonbit_string_t _M0L2idS2485;
      moonbit_string_t _M0L6_2atmpS2484;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L1rS1091;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L1vS1089;
      moonbit_string_t _M0L2idS2487;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L7_2abindS1092;
      moonbit_string_t _M0L6_2atmpS2486;
      int32_t _M0L6_2atmpS2488;
      if (_M0L1iS1087 > 0) {
        moonbit_incref(_M0L3recS1088);
        #line 144 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L4dbufS1085, (moonbit_string_t)moonbit_string_literal_24.data);
      } else {
        moonbit_incref(_M0L3recS1088);
      }
      _M0L2idS2485 = _M0L3recS1088->$1;
      moonbit_incref(_M0L2idS2485);
      #line 146 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L6_2atmpS2484
      = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L2idS2485);
      moonbit_decref(_M0L2idS2485);
      #line 146 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L4dbufS1085, _M0L6_2atmpS2484);
      moonbit_decref(_M0L6_2atmpS2484);
      #line 147 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L4dbufS1085, (moonbit_string_t)moonbit_string_literal_25.data);
      _M0L2idS2487 = _M0L3recS1088->$1;
      moonbit_incref(_M0L2idS2487);
      #line 148 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L7_2abindS1092
      = _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(_M0L1dS1084, _M0L2idS2487);
      moonbit_decref(_M0L2idS2487);
      if (_M0L7_2abindS1092 == 0) {
        if (_M0L7_2abindS1092) {
          moonbit_decref(_M0L7_2abindS1092);
        }
        _M0L1vS1089 = _M0L3recS1088;
      } else {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L7_2aSomeS1093;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4_2arS1094;
        moonbit_decref(_M0L3recS1088);
        _M0L7_2aSomeS1093 = _M0L7_2abindS1092;
        _M0L4_2arS1094 = _M0L7_2aSomeS1093;
        _M0L1rS1091 = _M0L4_2arS1094;
        goto join_1090;
      }
      goto joinlet_2973;
      join_1090:;
      _M0L1vS1089 = _M0L1rS1091;
      joinlet_2973:;
      #line 152 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L6_2atmpS2486
      = _M0FP46biolab8bio__seq3cmd11seqio__main24record__to__json_2einner(_M0L1vS1089, 0);
      moonbit_decref(_M0L1vS1089);
      #line 152 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L4dbufS1085, _M0L6_2atmpS2486);
      moonbit_decref(_M0L6_2atmpS2486);
      _M0L6_2atmpS2488 = _M0L1iS1087 + 1;
      _M0L1iS1087 = _M0L6_2atmpS2488;
      continue;
    } else {
      moonbit_decref(_M0L1dS1084);
      moonbit_decref(_M0L4recsS1069);
    }
    break;
  }
  #line 154 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L4dbufS1085, 125);
  #line 155 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2490 = _M0MPB13StringBuilder10to__string(_M0L4dbufS1085);
  moonbit_decref(_M0L4dbufS1085);
  #line 155 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_26.data, _M0L6_2atmpS2490);
  moonbit_decref(_M0L6_2atmpS2490);
  #line 160 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2519
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_27.data);
  #line 159 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2514
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2519, (moonbit_string_t)moonbit_string_literal_28.data, (moonbit_string_t)moonbit_string_literal_28.data, (moonbit_string_t)moonbit_string_literal_29.data);
  moonbit_decref(_M0L6_2atmpS2519);
  #line 166 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2518
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_30.data);
  #line 165 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2515
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2518, (moonbit_string_t)moonbit_string_literal_31.data, (moonbit_string_t)moonbit_string_literal_31.data, (moonbit_string_t)moonbit_string_literal_32.data);
  moonbit_decref(_M0L6_2atmpS2518);
  #line 172 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2517
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_27.data);
  #line 171 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2516
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2517, (moonbit_string_t)moonbit_string_literal_33.data, (moonbit_string_t)moonbit_string_literal_33.data, (moonbit_string_t)moonbit_string_literal_33.data);
  moonbit_decref(_M0L6_2atmpS2517);
  _M0L6_2atmpS2513
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(3);
  _M0L6_2atmpS2513[0] = _M0L6_2atmpS2514;
  _M0L6_2atmpS2513[1] = _M0L6_2atmpS2515;
  _M0L6_2atmpS2513[2] = _M0L6_2atmpS2516;
  _M0L11write__recsS1096
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L11write__recsS1096)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L11write__recsS1096->$0 = _M0L6_2atmpS2513;
  _M0L11write__recsS1096->$1 = 3;
  #line 178 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2975
  = _M0FP26biolab8bio__seq12seqio__write(_M0L11write__recsS1096, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L11write__recsS1096);
  if (_tmp_2975.tag) {
    moonbit_string_t const _M0L5_2aokS2511 = _tmp_2975.data.ok;
    _M0L10fasta__outS1097 = _M0L5_2aokS2511;
  } else {
    void* const _M0L6_2aerrS2512 = _tmp_2975.data.err;
    _M0L11_2atry__errS1099 = _M0L6_2aerrS2512;
    goto join_1098;
  }
  goto joinlet_2974;
  join_1098:;
  moonbit_decref(_M0L11_2atry__errS1099);
  _M0L10fasta__outS1097 = (moonbit_string_t)moonbit_string_literal_0.data;
  joinlet_2974:;
  #line 179 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2491
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L10fasta__outS1097);
  moonbit_decref(_M0L10fasta__outS1097);
  #line 179 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_34.data, _M0L6_2atmpS2491);
  moonbit_decref(_M0L6_2atmpS2491);
  #line 182 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2977
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L12fastq__multiS1068, (moonbit_string_t)moonbit_string_literal_35.data);
  moonbit_decref(_M0L12fastq__multiS1068);
  if (_tmp_2977.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2509 =
      _tmp_2977.data.ok;
    _M0L8fq__recsS1100 = _M0L5_2aokS2509;
  } else {
    void* const _M0L6_2aerrS2510 = _tmp_2977.data.err;
    _M0L11_2atry__errS1102 = _M0L6_2aerrS2510;
    goto join_1101;
  }
  goto joinlet_2976;
  join_1101:;
  moonbit_decref(_M0L11_2atry__errS1102);
  #line 182 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L8fq__recsS1100
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  joinlet_2976:;
  #line 185 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2492
  = _M0FP46biolab8bio__seq3cmd11seqio__main25records__to__json_2einner(_M0L8fq__recsS1100, 1);
  moonbit_decref(_M0L8fq__recsS1100);
  #line 185 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_36.data, _M0L6_2atmpS2492);
  moonbit_decref(_M0L6_2atmpS2492);
  #line 188 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L10fq1__qualsS1103 = _M0MPC15array5Array11new_2einnerGiE(0);
  _M0L1iS1104 = 0;
  while (1) {
    if (_M0L1iS1104 < 39) {
      int32_t _M0L6_2atmpS2493;
      #line 190 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0MPC15array5Array4pushGiE(_M0L10fq1__qualsS1103, 40);
      _M0L6_2atmpS2493 = _M0L1iS1104 + 1;
      _M0L1iS1104 = _M0L6_2atmpS2493;
      continue;
    }
    break;
  }
  #line 194 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2508
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_27.data);
  #line 193 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2507
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2508, (moonbit_string_t)moonbit_string_literal_37.data, (moonbit_string_t)moonbit_string_literal_37.data, (moonbit_string_t)moonbit_string_literal_38.data);
  moonbit_decref(_M0L6_2atmpS2508);
  #line 193 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2501
  = _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(_M0L6_2atmpS2507, _M0L10fq1__qualsS1103);
  moonbit_decref(_M0L6_2atmpS2507);
  moonbit_decref(_M0L10fq1__qualsS1103);
  #line 200 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2506
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_39.data);
  #line 199 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2503
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2506, (moonbit_string_t)moonbit_string_literal_40.data, (moonbit_string_t)moonbit_string_literal_40.data, (moonbit_string_t)moonbit_string_literal_41.data);
  moonbit_decref(_M0L6_2atmpS2506);
  _M0L6_2atmpS2505 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS2505[0] = 10;
  _M0L6_2atmpS2505[1] = 20;
  _M0L6_2atmpS2505[2] = 30;
  _M0L6_2atmpS2505[3] = 40;
  _M0L6_2atmpS2504
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS2504)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
  _M0L6_2atmpS2504->$0 = _M0L6_2atmpS2505;
  _M0L6_2atmpS2504->$1 = 4;
  #line 199 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2502
  = _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(_M0L6_2atmpS2503, _M0L6_2atmpS2504);
  moonbit_decref(_M0L6_2atmpS2503);
  moonbit_decref(_M0L6_2atmpS2504);
  _M0L6_2atmpS2500
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(2);
  _M0L6_2atmpS2500[0] = _M0L6_2atmpS2501;
  _M0L6_2atmpS2500[1] = _M0L6_2atmpS2502;
  _M0L15write__fq__recsS1106
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L15write__fq__recsS1106)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L15write__fq__recsS1106->$0 = _M0L6_2atmpS2500;
  _M0L15write__fq__recsS1106->$1 = 2;
  #line 206 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2980
  = _M0FP26biolab8bio__seq12seqio__write(_M0L15write__fq__recsS1106, (moonbit_string_t)moonbit_string_literal_35.data);
  moonbit_decref(_M0L15write__fq__recsS1106);
  if (_tmp_2980.tag) {
    moonbit_string_t const _M0L5_2aokS2498 = _tmp_2980.data.ok;
    _M0L7fq__outS1107 = _M0L5_2aokS2498;
  } else {
    void* const _M0L6_2aerrS2499 = _tmp_2980.data.err;
    _M0L11_2atry__errS1109 = _M0L6_2aerrS2499;
    goto join_1108;
  }
  goto joinlet_2979;
  join_1108:;
  moonbit_decref(_M0L11_2atry__errS1109);
  _M0L7fq__outS1107 = (moonbit_string_t)moonbit_string_literal_0.data;
  joinlet_2979:;
  #line 207 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2494
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L7fq__outS1107);
  moonbit_decref(_M0L7fq__outS1107);
  #line 207 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_42.data, _M0L6_2atmpS2494);
  moonbit_decref(_M0L6_2atmpS2494);
  _M0L15genbank__sampleS1110
  = (moonbit_string_t)moonbit_string_literal_43.data;
  #line 211 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _tmp_2982
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L15genbank__sampleS1110, (moonbit_string_t)moonbit_string_literal_44.data);
  moonbit_decref(_M0L15genbank__sampleS1110);
  if (_tmp_2982.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2496 =
      _tmp_2982.data.ok;
    _M0L8gb__recsS1111 = _M0L5_2aokS2496;
  } else {
    void* const _M0L6_2aerrS2497 = _tmp_2982.data.err;
    _M0L11_2atry__errS1113 = _M0L6_2aerrS2497;
    goto join_1112;
  }
  goto joinlet_2981;
  join_1112:;
  moonbit_decref(_M0L11_2atry__errS1113);
  #line 211 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L8gb__recsS1111
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  joinlet_2981:;
  #line 214 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L7_2abindS1115
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L8gb__recsS1111);
  switch (_M0L7_2abindS1115) {
    case 0: {
      moonbit_decref(_M0L8gb__recsS1111);
      _M0L8gb__jsonS1114 = (moonbit_string_t)moonbit_string_literal_45.data;
      break;
    }
    default: {
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2495;
      #line 216 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L6_2atmpS2495
      = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8gb__recsS1111, 0);
      moonbit_decref(_M0L8gb__recsS1111);
      #line 216 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L8gb__jsonS1114
      = _M0FP46biolab8bio__seq3cmd11seqio__main24record__to__json_2einner(_M0L6_2atmpS2495, 0);
      moonbit_decref(_M0L6_2atmpS2495);
      break;
    }
  }
  #line 218 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main4emit((moonbit_string_t)moonbit_string_literal_46.data, _M0L8gb__jsonS1114);
  moonbit_decref(_M0L8gb__jsonS1114);
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main4emit(
  moonbit_string_t _M0L4nameS1062,
  moonbit_string_t _M0L7payloadS1063
) {
  moonbit_string_t _M0L6_2atmpS2479;
  moonbit_string_t _M0L6_2atmpS2478;
  #line 85 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  #line 86 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2479
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_47.data, _M0L4nameS1062);
  #line 86 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2478
  = moonbit_add_string(_M0L6_2atmpS2479, (moonbit_string_t)moonbit_string_literal_48.data);
  moonbit_decref(_M0L6_2atmpS2479);
  #line 86 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS2478);
  moonbit_decref(_M0L6_2atmpS2478);
  #line 87 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FPB7printlnGsE(_M0L7payloadS1063);
  #line 88 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_49.data);
  return 0;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main25records__to__json_2einner(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1057,
  int32_t _M0L13include__qualS1060
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1055;
  int32_t _M0L7_2abindS1056;
  int32_t _M0L1iS1058;
  moonbit_string_t _result_2984;
  #line 67 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  #line 71 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L3bufS1055 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 72 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1055, 91);
  _M0L7_2abindS1056 = _M0L7recordsS1057->$1;
  _M0L1iS1058 = 0;
  while (1) {
    if (_M0L1iS1058 < _M0L7_2abindS1056) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2477 =
        _M0L7recordsS1057->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1059 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2477[_M0L1iS1058];
      moonbit_string_t _M0L6_2atmpS2475;
      int32_t _M0L6_2atmpS2476;
      if (_M0L1iS1058 > 0) {
        moonbit_incref(_M0L3recS1059);
        #line 75 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1055, (moonbit_string_t)moonbit_string_literal_24.data);
      } else {
        moonbit_incref(_M0L3recS1059);
      }
      #line 77 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L6_2atmpS2475
      = _M0FP46biolab8bio__seq3cmd11seqio__main24record__to__json_2einner(_M0L3recS1059, _M0L13include__qualS1060);
      moonbit_decref(_M0L3recS1059);
      #line 77 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1055, _M0L6_2atmpS2475);
      moonbit_decref(_M0L6_2atmpS2475);
      _M0L6_2atmpS2476 = _M0L1iS1058 + 1;
      _M0L1iS1058 = _M0L6_2atmpS2476;
      continue;
    }
    break;
  }
  #line 79 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1055, 93);
  #line 80 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _result_2984 = _M0MPB13StringBuilder10to__string(_M0L3bufS1055);
  moonbit_decref(_M0L3bufS1055);
  return _result_2984;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main24record__to__json_2einner(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1043,
  int32_t _M0L13include__qualS1044
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1042;
  moonbit_string_t _M0L2idS2462;
  moonbit_string_t _M0L6_2atmpS2461;
  moonbit_string_t _M0L4nameS2464;
  moonbit_string_t _M0L6_2atmpS2463;
  moonbit_string_t _M0L11descriptionS2466;
  moonbit_string_t _M0L6_2atmpS2465;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2469;
  moonbit_string_t _M0L6_2atmpS2468;
  moonbit_string_t _M0L6_2atmpS2467;
  int32_t _M0L6_2atmpS2471;
  moonbit_string_t _M0L6_2atmpS2470;
  moonbit_string_t _result_2987;
  #line 34 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L3bufS1042 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 36 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_50.data);
  _M0L2idS2462 = _M0L3recS1043->$1;
  moonbit_incref(_M0L2idS2462);
  #line 37 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2461
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L2idS2462);
  moonbit_decref(_M0L2idS2462);
  #line 37 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, _M0L6_2atmpS2461);
  moonbit_decref(_M0L6_2atmpS2461);
  #line 38 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_51.data);
  _M0L4nameS2464 = _M0L3recS1043->$2;
  moonbit_incref(_M0L4nameS2464);
  #line 39 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2463
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L4nameS2464);
  moonbit_decref(_M0L4nameS2464);
  #line 39 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, _M0L6_2atmpS2463);
  moonbit_decref(_M0L6_2atmpS2463);
  #line 40 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_52.data);
  _M0L11descriptionS2466 = _M0L3recS1043->$3;
  moonbit_incref(_M0L11descriptionS2466);
  #line 41 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2465
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L11descriptionS2466);
  moonbit_decref(_M0L11descriptionS2466);
  #line 41 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, _M0L6_2atmpS2465);
  moonbit_decref(_M0L6_2atmpS2465);
  #line 42 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_53.data);
  _M0L3seqS2469 = _M0L3recS1043->$0;
  moonbit_incref(_M0L3seqS2469);
  #line 43 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2468 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2469);
  moonbit_decref(_M0L3seqS2469);
  #line 43 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2467
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(_M0L6_2atmpS2468);
  moonbit_decref(_M0L6_2atmpS2468);
  #line 43 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, _M0L6_2atmpS2467);
  moonbit_decref(_M0L6_2atmpS2467);
  #line 44 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_54.data);
  #line 45 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2471 = _M0MP26biolab8bio__seq9SeqRecord6length(_M0L3recS1043);
  #line 45 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2470
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2471, 10);
  #line 45 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, _M0L6_2atmpS2470);
  moonbit_decref(_M0L6_2atmpS2470);
  if (_M0L13include__qualS1044) {
    struct _M0TPB5ArrayGiE* _M0L1qS1047;
    struct _M0TPB5ArrayGiE* _M0L5qualsS1045;
    struct _M0TPB5ArrayGiE* _M0L7_2abindS1048;
    int32_t _M0L7_2abindS1051;
    int32_t _M0L1iS1052;
    #line 47 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_55.data);
    #line 48 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
    _M0L7_2abindS1048
    = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L3recS1043);
    if (_M0L7_2abindS1048 == 0) {
      if (_M0L7_2abindS1048) {
        moonbit_decref(_M0L7_2abindS1048);
      }
      #line 50 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0L5qualsS1045 = _M0MPC15array5Array11new_2einnerGiE(0);
    } else {
      struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1049 = _M0L7_2abindS1048;
      struct _M0TPB5ArrayGiE* _M0L4_2aqS1050 = _M0L7_2aSomeS1049;
      _M0L1qS1047 = _M0L4_2aqS1050;
      goto join_1046;
    }
    goto joinlet_2985;
    join_1046:;
    _M0L5qualsS1045 = _M0L1qS1047;
    joinlet_2985:;
    #line 52 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
    _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1042, 91);
    _M0L7_2abindS1051 = _M0L5qualsS1045->$1;
    _M0L1iS1052 = 0;
    while (1) {
      if (_M0L1iS1052 < _M0L7_2abindS1051) {
        int32_t* _M0L3bufS2474 = _M0L5qualsS1045->$0;
        int32_t _M0L1qS1053 = (int32_t)_M0L3bufS2474[_M0L1iS1052];
        moonbit_string_t _M0L6_2atmpS2472;
        int32_t _M0L6_2atmpS2473;
        if (_M0L1iS1052 > 0) {
          #line 55 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, (moonbit_string_t)moonbit_string_literal_24.data);
        }
        #line 57 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
        _M0L6_2atmpS2472
        = _M0MPC13int3Int18to__string_2einner(_M0L1qS1053, 10);
        #line 57 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1042, _M0L6_2atmpS2472);
        moonbit_decref(_M0L6_2atmpS2472);
        _M0L6_2atmpS2473 = _M0L1iS1052 + 1;
        _M0L1iS1052 = _M0L6_2atmpS2473;
        continue;
      } else {
        moonbit_decref(_M0L5qualsS1045);
      }
      break;
    }
    #line 59 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
    _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1042, 93);
  }
  #line 61 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1042, 125);
  #line 62 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _result_2987 = _M0MPB13StringBuilder10to__string(_M0L3bufS1042);
  moonbit_decref(_M0L3bufS1042);
  return _result_2987;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main12json__string(
  moonbit_string_t _M0L1sS1041
) {
  moonbit_string_t _M0L6_2atmpS2460;
  moonbit_string_t _M0L6_2atmpS2459;
  moonbit_string_t _result_2988;
  #line 27 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  #line 28 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2460
  = _M0FP46biolab8bio__seq3cmd11seqio__main12json__escape(_M0L1sS1041);
  #line 28 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2459
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_56.data, _M0L6_2atmpS2460);
  moonbit_decref(_M0L6_2atmpS2460);
  #line 28 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _result_2988
  = moonbit_add_string(_M0L6_2atmpS2459, (moonbit_string_t)moonbit_string_literal_56.data);
  moonbit_decref(_M0L6_2atmpS2459);
  return _result_2988;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main12json__escape(
  moonbit_string_t _M0L1sS1033
) {
  int32_t _M0L6_2atmpS2458;
  int32_t _M0L6_2atmpS2457;
  struct _M0TPB13StringBuilder* _M0L3bufS1032;
  struct _M0TPB4IterGcE* _M0L5_2aitS1034;
  moonbit_string_t _result_2991;
  #line 11 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L6_2atmpS2458 = Moonbit_array_length(_M0L1sS1033);
  _M0L6_2atmpS2457 = _M0L6_2atmpS2458 + 2;
  #line 12 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L3bufS1032
  = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L6_2atmpS2457);
  #line 13 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _M0L5_2aitS1034 = _M0MPC16string6String4iter(_M0L1sS1033);
  while (1) {
    int32_t _M0L1cS1036;
    int32_t _M0L7_2abindS1038;
    #line 13 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
    _M0L7_2abindS1038 = _M0MPB4Iter4nextGcE(_M0L5_2aitS1034);
    if (_M0L7_2abindS1038 == -1) {
      moonbit_decref(_M0L5_2aitS1034);
    } else {
      int32_t _M0L7_2aSomeS1039 = _M0L7_2abindS1038;
      int32_t _M0L4_2acS1040 = _M0L7_2aSomeS1039;
      _M0L1cS1036 = _M0L4_2acS1040;
      goto join_1035;
    }
    goto joinlet_2990;
    join_1035:;
    if (_M0L1cS1036 == 92) {
      #line 15 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1032, (moonbit_string_t)moonbit_string_literal_57.data);
    } else if (_M0L1cS1036 == 34) {
      #line 16 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1032, (moonbit_string_t)moonbit_string_literal_58.data);
    } else if (_M0L1cS1036 == 10) {
      #line 17 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1032, (moonbit_string_t)moonbit_string_literal_59.data);
    } else if (_M0L1cS1036 == 9) {
      #line 18 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1032, (moonbit_string_t)moonbit_string_literal_60.data);
    } else {
      #line 19 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1032, _M0L1cS1036);
    }
    continue;
    joinlet_2990:;
    break;
  }
  #line 22 "/home/developer/Documents/1/cmd/seqio_main/seqio_tests.mbt"
  _result_2991 = _M0MPB13StringBuilder10to__string(_M0L3bufS1032);
  moonbit_decref(_M0L3bufS1032);
  return _result_2991;
}

struct _M0TPB5ArrayGiE* _M0MP26biolab8bio__seq9SeqRecord14phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS1031
) {
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2456;
  struct _M0TPB5ArrayGiE* _result_2992;
  #line 78 "/home/developer/Documents/1/seq_record.mbt"
  _M0L19letter__annotationsS2456 = _M0L4selfS1031->$4;
  moonbit_incref(_M0L19letter__annotationsS2456);
  #line 79 "/home/developer/Documents/1/seq_record.mbt"
  _result_2992
  = _M0MPB3Map3getGsRPB5ArrayGiEE(_M0L19letter__annotationsS2456, (moonbit_string_t)moonbit_string_literal_61.data);
  moonbit_decref(_M0L19letter__annotationsS2456);
  return _result_2992;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS1019,
  struct _M0TPB5ArrayGiE* _M0L5qualsS1030
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1018;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2455;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2450;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2454;
  int32_t _M0L6_2atmpS2453;
  int32_t _M0L6_2atmpS2452;
  int64_t _M0L6_2atmpS2451;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L16new__annotationsS1017;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2445;
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0L5_2aitS1020;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2446;
  moonbit_string_t _M0L2idS2447;
  moonbit_string_t _M0L4nameS2448;
  moonbit_string_t _M0L11descriptionS2449;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_2995;
  #line 61 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS1018 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2455 = _M0L7_2abindS1018;
  _M0L6_2atmpS2450
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2455, .$1 = 0, .$2 = 0
  };
  _M0L19letter__annotationsS2454 = _M0L4selfS1019->$4;
  moonbit_incref(_M0L19letter__annotationsS2454);
  #line 67 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS2453
  = _M0MPB3Map6lengthGsRPB5ArrayGiEE(_M0L19letter__annotationsS2454);
  moonbit_decref(_M0L19letter__annotationsS2454);
  _M0L6_2atmpS2452 = _M0L6_2atmpS2453 + 1;
  _M0L6_2atmpS2451 = (int64_t)_M0L6_2atmpS2452;
  #line 65 "/home/developer/Documents/1/seq_record.mbt"
  _M0L16new__annotationsS1017
  = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2450, _M0L6_2atmpS2451);
  moonbit_decref(_M0L6_2atmpS2450.$0);
  _M0L19letter__annotationsS2445 = _M0L4selfS1019->$4;
  moonbit_incref(_M0L19letter__annotationsS2445);
  #line 65 "/home/developer/Documents/1/seq_record.mbt"
  _M0L5_2aitS1020
  = _M0MPB3Map5iter2GsRPB5ArrayGiEE(_M0L19letter__annotationsS2445);
  moonbit_decref(_M0L19letter__annotationsS2445);
  while (1) {
    moonbit_string_t _M0L1kS1022;
    struct _M0TPB5ArrayGiE* _M0L1vS1023;
    struct _M0TUsRPB5ArrayGiEE* _M0L7_2abindS1025;
    #line 69 "/home/developer/Documents/1/seq_record.mbt"
    _M0L7_2abindS1025 = _M0MPB5Iter24nextGsRPB5ArrayGiEE(_M0L5_2aitS1020);
    if (_M0L7_2abindS1025 == 0) {
      if (_M0L7_2abindS1025) {
        moonbit_decref(_M0L7_2abindS1025);
      }
      moonbit_decref(_M0L5_2aitS1020);
    } else {
      struct _M0TUsRPB5ArrayGiEE* _M0L7_2aSomeS1026 = _M0L7_2abindS1025;
      struct _M0TUsRPB5ArrayGiEE* _M0L4_2axS1027 = _M0L7_2aSomeS1026;
      moonbit_string_t _M0L4_2akS1028 = _M0L4_2axS1027->$0;
      struct _M0TPB5ArrayGiE* _M0L8_2afieldS2693 = _M0L4_2axS1027->$1;
      int32_t _M0L6_2acntS2895 =
        Moonbit_rc_count(Moonbit_object_header(_M0L4_2axS1027));
      struct _M0TPB5ArrayGiE* _M0L4_2avS1029;
      if (_M0L6_2acntS2895 > 1) {
        int32_t _M0L11_2anew__cntS2896 = _M0L6_2acntS2895 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L4_2axS1027), _M0L11_2anew__cntS2896);
        moonbit_incref(_M0L8_2afieldS2693);
        moonbit_incref(_M0L4_2akS1028);
      } else if (_M0L6_2acntS2895 == 1) {
        #line 69 "/home/developer/Documents/1/seq_record.mbt"
        moonbit_free(_M0L4_2axS1027);
      }
      _M0L4_2avS1029 = _M0L8_2afieldS2693;
      _M0L1kS1022 = _M0L4_2akS1028;
      _M0L1vS1023 = _M0L4_2avS1029;
      goto join_1021;
    }
    goto joinlet_2994;
    join_1021:;
    #line 70 "/home/developer/Documents/1/seq_record.mbt"
    _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L16new__annotationsS1017, _M0L1kS1022, _M0L1vS1023);
    moonbit_decref(_M0L1kS1022);
    moonbit_decref(_M0L1vS1023);
    continue;
    joinlet_2994:;
    break;
  }
  #line 72 "/home/developer/Documents/1/seq_record.mbt"
  _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L16new__annotationsS1017, (moonbit_string_t)moonbit_string_literal_61.data, _M0L5qualsS1030);
  _M0L3seqS2446 = _M0L4selfS1019->$0;
  _M0L2idS2447 = _M0L4selfS1019->$1;
  _M0L4nameS2448 = _M0L4selfS1019->$2;
  _M0L11descriptionS2449 = _M0L4selfS1019->$3;
  moonbit_incref(_M0L3seqS2446);
  moonbit_incref(_M0L2idS2447);
  moonbit_incref(_M0L4nameS2448);
  moonbit_incref(_M0L11descriptionS2449);
  _block_2995
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_2995)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
  _block_2995->$0 = _M0L3seqS2446;
  _block_2995->$1 = _M0L2idS2447;
  _block_2995->$2 = _M0L4nameS2448;
  _block_2995->$3 = _M0L11descriptionS2449;
  _block_2995->$4 = _M0L16new__annotationsS1017;
  return _block_2995;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord9from__seq(
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS1015
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1016;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2444;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2443;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2442;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_2996;
  #line 48 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS1016 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2444 = _M0L7_2abindS1016;
  _M0L6_2atmpS2443
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2444, .$1 = 0, .$2 = 0
  };
  #line 54 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS2442 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2443, 0ll);
  moonbit_decref(_M0L6_2atmpS2443.$0);
  moonbit_incref(_M0L3seqS1015);
  _block_2996
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_2996)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
  _block_2996->$0 = _M0L3seqS1015;
  _block_2996->$1 = (moonbit_string_t)moonbit_string_literal_62.data;
  _block_2996->$2 = (moonbit_string_t)moonbit_string_literal_63.data;
  _block_2996->$3 = (moonbit_string_t)moonbit_string_literal_64.data;
  _block_2996->$4 = _M0L6_2atmpS2442;
  return _block_2996;
}

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq15seqio__to__dict(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1011
) {
  struct _M0TUsRP26biolab8bio__seq9SeqRecordE** _M0L7_2abindS1009;
  struct _M0TUsRP26biolab8bio__seq9SeqRecordE** _M0L6_2atmpS2441;
  struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE _M0L6_2atmpS2440;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L1dS1008;
  int32_t _M0L7_2abindS1010;
  int32_t _M0L2__S1012;
  #line 67 "/home/developer/Documents/1/seqio.mbt"
  _M0L7_2abindS1009
  = (struct _M0TUsRP26biolab8bio__seq9SeqRecordE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2441 = _M0L7_2abindS1009;
  _M0L6_2atmpS2440
  = (struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE){
    .$0 = _M0L6_2atmpS2441, .$1 = 0, .$2 = 0
  };
  #line 68 "/home/developer/Documents/1/seqio.mbt"
  _M0L1dS1008
  = _M0MPB3Map3MapGsRP26biolab8bio__seq9SeqRecordE(_M0L6_2atmpS2440, 16ll);
  moonbit_decref(_M0L6_2atmpS2440.$0);
  _M0L7_2abindS1010 = _M0L7recordsS1011->$1;
  _M0L2__S1012 = 0;
  while (1) {
    if (_M0L2__S1012 < _M0L7_2abindS1010) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2439 =
        _M0L7recordsS1011->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1013 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2439[_M0L2__S1012];
      moonbit_string_t _M0L2idS2437 = _M0L3recS1013->$1;
      int32_t _M0L6_2atmpS2438;
      moonbit_incref(_M0L2idS2437);
      moonbit_incref(_M0L3recS1013);
      #line 70 "/home/developer/Documents/1/seqio.mbt"
      _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(_M0L1dS1008, _M0L2idS2437, _M0L3recS1013);
      moonbit_decref(_M0L2idS2437);
      moonbit_decref(_M0L3recS1013);
      _M0L6_2atmpS2438 = _M0L2__S1012 + 1;
      _M0L2__S1012 = _M0L6_2atmpS2438;
      continue;
    }
    break;
  }
  return _M0L1dS1008;
}

struct moonbit_result_3 _M0FP26biolab8bio__seq12seqio__write(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1006,
  moonbit_string_t _M0L6formatS1007
) {
  #line 51 "/home/developer/Documents/1/seqio.mbt"
  if (
    _M0L6formatS1007 == (moonbit_string_t)moonbit_string_literal_16.data
    || Moonbit_array_length(_M0L6formatS1007) == 5
       && 0
          == memcmp(_M0L6formatS1007, (moonbit_string_t)moonbit_string_literal_16.data, 10)
  ) {
    moonbit_string_t _M0L6_2atmpS2434;
    struct moonbit_result_3 _result_3000;
    #line 56 "/home/developer/Documents/1/seqio.mbt"
    _M0L6_2atmpS2434
    = _M0FP26biolab8bio__seq12write__fasta(_M0L7recordsS1006);
    _result_3000.tag = 1;
    _result_3000.data.ok = _M0L6_2atmpS2434;
    return _result_3000;
  } else if (
           _M0L6formatS1007
           == (moonbit_string_t)moonbit_string_literal_35.data
           || Moonbit_array_length(_M0L6formatS1007) == 5
              && 0
                 == memcmp(_M0L6formatS1007, (moonbit_string_t)moonbit_string_literal_35.data, 10)
         ) {
    goto join_1005;
  } else if (
           _M0L6formatS1007
           == (moonbit_string_t)moonbit_string_literal_66.data
           || Moonbit_array_length(_M0L6formatS1007) == 12
              && 0
                 == memcmp(_M0L6formatS1007, (moonbit_string_t)moonbit_string_literal_66.data, 24)
         ) {
    goto join_1005;
  } else {
    moonbit_string_t _M0L6_2atmpS2436;
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2435;
    struct moonbit_result_3 _result_2999;
    #line 58 "/home/developer/Documents/1/seqio.mbt"
    _M0L6_2atmpS2436
    = moonbit_add_string((moonbit_string_t)moonbit_string_literal_65.data, _M0L6formatS1007);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2435
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2435)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2435)->$0
    = _M0L6_2atmpS2436;
    _result_2999.tag = 0;
    _result_2999.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2435;
    return _result_2999;
  }
  join_1005:;
  #line 57 "/home/developer/Documents/1/seqio.mbt"
  return _M0FP26biolab8bio__seq12write__fastq(_M0L7recordsS1006);
}

struct moonbit_result_3 _M0FP26biolab8bio__seq12write__fastq(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS991
) {
  struct _M0TPB13StringBuilder* _M0L3bufS989;
  int32_t _M0L7_2abindS990;
  int32_t _M0L2__S992;
  moonbit_string_t _M0L6_2atmpS2433;
  struct moonbit_result_3 _result_3006;
  #line 123 "/home/developer/Documents/1/fastq_io.mbt"
  #line 124 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L3bufS989 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS990 = _M0L7recordsS991->$1;
  _M0L2__S992 = 0;
  while (1) {
    if (_M0L2__S992 < _M0L7_2abindS990) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2432 =
        _M0L7recordsS991->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS993 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2432[_M0L2__S992];
      struct _M0TPB5ArrayGiE* _M0L1qS996;
      struct _M0TPB5ArrayGiE* _M0L5qualsS994;
      struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2430 =
        _M0L3recS993->$4;
      struct _M0TPB5ArrayGiE* _M0L7_2abindS997;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2426;
      moonbit_string_t _M0L3seqS1000;
      int32_t _M0L6_2atmpS2403;
      int32_t _M0L6_2atmpS2404;
      moonbit_string_t _M0L2idS2417;
      moonbit_string_t _M0L8_2afieldS2700;
      int32_t _M0L6_2acntS2909;
      moonbit_string_t _M0L11descriptionS2418;
      moonbit_string_t _M0L6_2atmpS2416;
      int32_t _M0L6_2atmpS2425;
      uint16_t* _M0L4qbufS1001;
      int32_t _M0L1jS1002;
      moonbit_string_t _M0L6_2atmpS2424;
      int32_t _M0L6_2atmpS2431;
      moonbit_incref(_M0L19letter__annotationsS2430);
      moonbit_incref(_M0L3recS993);
      #line 126 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L7_2abindS997
      = _M0MPB3Map3getGsRPB5ArrayGiEE(_M0L19letter__annotationsS2430, (moonbit_string_t)moonbit_string_literal_61.data);
      moonbit_decref(_M0L19letter__annotationsS2430);
      if (_M0L7_2abindS997 == 0) {
        moonbit_string_t _M0L8_2afieldS2704;
        int32_t _M0L6_2acntS2897;
        moonbit_string_t _M0L2idS2429;
        moonbit_string_t _M0L6_2atmpS2428;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2427;
        struct moonbit_result_3 _result_3003;
        if (_M0L7_2abindS997) {
          moonbit_decref(_M0L7_2abindS997);
        }
        moonbit_decref(_M0L3bufS989);
        _M0L8_2afieldS2704 = _M0L3recS993->$1;
        _M0L6_2acntS2897
        = Moonbit_rc_count(Moonbit_object_header(_M0L3recS993));
        if (_M0L6_2acntS2897 > 1) {
          int32_t _M0L11_2anew__cntS2902 = _M0L6_2acntS2897 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS993), _M0L11_2anew__cntS2902);
          moonbit_incref(_M0L8_2afieldS2704);
        } else if (_M0L6_2acntS2897 == 1) {
          struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS2901 =
            _M0L3recS993->$4;
          moonbit_string_t _M0L8_2afieldS2900;
          moonbit_string_t _M0L8_2afieldS2899;
          struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2898;
          moonbit_decref(_M0L8_2afieldS2901);
          _M0L8_2afieldS2900 = _M0L3recS993->$3;
          moonbit_decref(_M0L8_2afieldS2900);
          _M0L8_2afieldS2899 = _M0L3recS993->$2;
          moonbit_decref(_M0L8_2afieldS2899);
          _M0L8_2afieldS2898 = _M0L3recS993->$0;
          moonbit_decref(_M0L8_2afieldS2898);
          #line 130 "/home/developer/Documents/1/fastq_io.mbt"
          moonbit_free(_M0L3recS993);
        }
        _M0L2idS2429 = _M0L8_2afieldS2704;
        #line 130 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2428
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_67.data, _M0L2idS2429);
        moonbit_decref(_M0L2idS2429);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2427
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2427)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2427)->$0
        = _M0L6_2atmpS2428;
        _result_3003.tag = 0;
        _result_3003.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2427;
        return _result_3003;
      } else {
        struct _M0TPB5ArrayGiE* _M0L7_2aSomeS998 = _M0L7_2abindS997;
        struct _M0TPB5ArrayGiE* _M0L4_2aqS999 = _M0L7_2aSomeS998;
        _M0L1qS996 = _M0L4_2aqS999;
        goto join_995;
      }
      goto joinlet_3002;
      join_995:;
      _M0L5qualsS994 = _M0L1qS996;
      joinlet_3002:;
      _M0L3seqS2426 = _M0L3recS993->$0;
      moonbit_incref(_M0L3seqS2426);
      #line 133 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L3seqS1000 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2426);
      moonbit_decref(_M0L3seqS2426);
      _M0L6_2atmpS2403 = Moonbit_array_length(_M0L3seqS1000);
      #line 134 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2404 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS994);
      if (_M0L6_2atmpS2403 != _M0L6_2atmpS2404) {
        int32_t _M0L6_2atmpS2415;
        moonbit_string_t _M0L6_2atmpS2414;
        moonbit_string_t _M0L6_2atmpS2413;
        moonbit_string_t _M0L6_2atmpS2410;
        int32_t _M0L6_2atmpS2412;
        moonbit_string_t _M0L6_2atmpS2411;
        moonbit_string_t _M0L6_2atmpS2409;
        moonbit_string_t _M0L6_2atmpS2407;
        moonbit_string_t _M0L8_2afieldS2702;
        int32_t _M0L6_2acntS2903;
        moonbit_string_t _M0L2idS2408;
        moonbit_string_t _M0L6_2atmpS2406;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2405;
        struct moonbit_result_3 _result_3004;
        moonbit_decref(_M0L3bufS989);
        _M0L6_2atmpS2415 = Moonbit_array_length(_M0L3seqS1000);
        moonbit_decref(_M0L3seqS1000);
        #line 137 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2414
        = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2415, 10);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2413
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_68.data, _M0L6_2atmpS2414);
        moonbit_decref(_M0L6_2atmpS2414);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2410
        = moonbit_add_string(_M0L6_2atmpS2413, (moonbit_string_t)moonbit_string_literal_69.data);
        moonbit_decref(_M0L6_2atmpS2413);
        #line 139 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2412 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS994);
        moonbit_decref(_M0L5qualsS994);
        #line 139 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2411
        = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2412, 10);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2409
        = moonbit_add_string(_M0L6_2atmpS2410, _M0L6_2atmpS2411);
        moonbit_decref(_M0L6_2atmpS2411);
        moonbit_decref(_M0L6_2atmpS2410);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2407
        = moonbit_add_string(_M0L6_2atmpS2409, (moonbit_string_t)moonbit_string_literal_70.data);
        moonbit_decref(_M0L6_2atmpS2409);
        _M0L8_2afieldS2702 = _M0L3recS993->$1;
        _M0L6_2acntS2903
        = Moonbit_rc_count(Moonbit_object_header(_M0L3recS993));
        if (_M0L6_2acntS2903 > 1) {
          int32_t _M0L11_2anew__cntS2908 = _M0L6_2acntS2903 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS993), _M0L11_2anew__cntS2908);
          moonbit_incref(_M0L8_2afieldS2702);
        } else if (_M0L6_2acntS2903 == 1) {
          struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS2907 =
            _M0L3recS993->$4;
          moonbit_string_t _M0L8_2afieldS2906;
          moonbit_string_t _M0L8_2afieldS2905;
          struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2904;
          moonbit_decref(_M0L8_2afieldS2907);
          _M0L8_2afieldS2906 = _M0L3recS993->$3;
          moonbit_decref(_M0L8_2afieldS2906);
          _M0L8_2afieldS2905 = _M0L3recS993->$2;
          moonbit_decref(_M0L8_2afieldS2905);
          _M0L8_2afieldS2904 = _M0L3recS993->$0;
          moonbit_decref(_M0L8_2afieldS2904);
          #line 141 "/home/developer/Documents/1/fastq_io.mbt"
          moonbit_free(_M0L3recS993);
        }
        _M0L2idS2408 = _M0L8_2afieldS2702;
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2406 = moonbit_add_string(_M0L6_2atmpS2407, _M0L2idS2408);
        moonbit_decref(_M0L2idS2408);
        moonbit_decref(_M0L6_2atmpS2407);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2405
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2405)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2405)->$0
        = _M0L6_2atmpS2406;
        _result_3004.tag = 0;
        _result_3004.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2405;
        return _result_3004;
      }
      #line 145 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS989, 64);
      _M0L2idS2417 = _M0L3recS993->$1;
      _M0L8_2afieldS2700 = _M0L3recS993->$3;
      _M0L6_2acntS2909
      = Moonbit_rc_count(Moonbit_object_header(_M0L3recS993));
      if (_M0L6_2acntS2909 > 1) {
        int32_t _M0L11_2anew__cntS2913 = _M0L6_2acntS2909 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS993), _M0L11_2anew__cntS2913);
        moonbit_incref(_M0L8_2afieldS2700);
        moonbit_incref(_M0L2idS2417);
      } else if (_M0L6_2acntS2909 == 1) {
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS2912 =
          _M0L3recS993->$4;
        moonbit_string_t _M0L8_2afieldS2911;
        struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2910;
        moonbit_decref(_M0L8_2afieldS2912);
        _M0L8_2afieldS2911 = _M0L3recS993->$2;
        moonbit_decref(_M0L8_2afieldS2911);
        _M0L8_2afieldS2910 = _M0L3recS993->$0;
        moonbit_decref(_M0L8_2afieldS2910);
        #line 146 "/home/developer/Documents/1/fastq_io.mbt"
        moonbit_free(_M0L3recS993);
      }
      _M0L11descriptionS2418 = _M0L8_2afieldS2700;
      #line 146 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2416
      = _M0FP26biolab8bio__seq19build__fasta__title(_M0L2idS2417, _M0L11descriptionS2418);
      moonbit_decref(_M0L2idS2417);
      moonbit_decref(_M0L11descriptionS2418);
      #line 146 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS989, _M0L6_2atmpS2416);
      moonbit_decref(_M0L6_2atmpS2416);
      #line 147 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS989, 10);
      #line 149 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS989, _M0L3seqS1000);
      moonbit_decref(_M0L3seqS1000);
      #line 150 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS989, 10);
      #line 152 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS989, (moonbit_string_t)moonbit_string_literal_71.data);
      #line 154 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2425 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS994);
      _M0L4qbufS1001 = (uint16_t*)moonbit_make_string(_M0L6_2atmpS2425, 0);
      _M0L1jS1002 = 0;
      while (1) {
        int32_t _M0L6_2atmpS2419;
        #line 155 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2419 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS994);
        if (_M0L1jS1002 < _M0L6_2atmpS2419) {
          int32_t _M0L6_2atmpS2422;
          int32_t _M0L6_2atmpS2421;
          int32_t _M0L6_2atmpS2420;
          int32_t _M0L6_2atmpS2423;
          #line 156 "/home/developer/Documents/1/fastq_io.mbt"
          _M0L6_2atmpS2422
          = _M0MPC15array5Array2atGiE(_M0L5qualsS994, _M0L1jS1002);
          _M0L6_2atmpS2421 = _M0L6_2atmpS2422 + 33;
          _M0L6_2atmpS2420 = (uint16_t)_M0L6_2atmpS2421;
          if (
            _M0L1jS1002 < 0
            || _M0L1jS1002 >= Moonbit_array_length(_M0L4qbufS1001)
          ) {
            #line 156 "/home/developer/Documents/1/fastq_io.mbt"
            moonbit_panic();
          }
          _M0L4qbufS1001[_M0L1jS1002] = _M0L6_2atmpS2420;
          _M0L6_2atmpS2423 = _M0L1jS1002 + 1;
          _M0L1jS1002 = _M0L6_2atmpS2423;
          continue;
        } else {
          moonbit_decref(_M0L5qualsS994);
        }
        break;
      }
      _M0L6_2atmpS2424 = _M0L4qbufS1001;
      #line 158 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS989, _M0L6_2atmpS2424);
      moonbit_decref(_M0L6_2atmpS2424);
      #line 159 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS989, 10);
      _M0L6_2atmpS2431 = _M0L2__S992 + 1;
      _M0L2__S992 = _M0L6_2atmpS2431;
      continue;
    }
    break;
  }
  #line 161 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L6_2atmpS2433 = _M0MPB13StringBuilder10to__string(_M0L3bufS989);
  moonbit_decref(_M0L3bufS989);
  _result_3006.tag = 1;
  _result_3006.data.ok = _M0L6_2atmpS2433;
  return _result_3006;
}

moonbit_string_t _M0FP26biolab8bio__seq12write__fasta(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS979
) {
  struct _M0TPB13StringBuilder* _M0L3bufS977;
  int32_t _M0L7_2abindS978;
  int32_t _M0L2__S980;
  moonbit_string_t _result_3009;
  #line 99 "/home/developer/Documents/1/fasta_io.mbt"
  #line 100 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L3bufS977 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS978 = _M0L7recordsS979->$1;
  _M0L2__S980 = 0;
  while (1) {
    if (_M0L2__S980 < _M0L7_2abindS978) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2402 =
        _M0L7recordsS979->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS981 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2402[_M0L2__S980];
      moonbit_string_t _M0L2idS2399 = _M0L3recS981->$1;
      moonbit_string_t _M0L11descriptionS2400 = _M0L3recS981->$3;
      moonbit_string_t _M0L5titleS982;
      struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2708;
      int32_t _M0L6_2acntS2914;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2398;
      moonbit_string_t _M0L3seqS983;
      int32_t _M0L1nS984;
      struct _M0TPB8MutLocalGiE* _M0L1iS985;
      int32_t _M0L6_2atmpS2401;
      moonbit_incref(_M0L11descriptionS2400);
      moonbit_incref(_M0L2idS2399);
      moonbit_incref(_M0L3recS981);
      #line 107 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L5titleS982
      = _M0FP26biolab8bio__seq19build__fasta__title(_M0L2idS2399, _M0L11descriptionS2400);
      moonbit_decref(_M0L2idS2399);
      moonbit_decref(_M0L11descriptionS2400);
      #line 108 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS977, 62);
      #line 109 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS977, _M0L5titleS982);
      moonbit_decref(_M0L5titleS982);
      #line 110 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS977, 10);
      _M0L8_2afieldS2708 = _M0L3recS981->$0;
      _M0L6_2acntS2914
      = Moonbit_rc_count(Moonbit_object_header(_M0L3recS981));
      if (_M0L6_2acntS2914 > 1) {
        int32_t _M0L11_2anew__cntS2919 = _M0L6_2acntS2914 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS981), _M0L11_2anew__cntS2919);
        moonbit_incref(_M0L8_2afieldS2708);
      } else if (_M0L6_2acntS2914 == 1) {
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS2918 =
          _M0L3recS981->$4;
        moonbit_string_t _M0L8_2afieldS2917;
        moonbit_string_t _M0L8_2afieldS2916;
        moonbit_string_t _M0L8_2afieldS2915;
        moonbit_decref(_M0L8_2afieldS2918);
        _M0L8_2afieldS2917 = _M0L3recS981->$3;
        moonbit_decref(_M0L8_2afieldS2917);
        _M0L8_2afieldS2916 = _M0L3recS981->$2;
        moonbit_decref(_M0L8_2afieldS2916);
        _M0L8_2afieldS2915 = _M0L3recS981->$1;
        moonbit_decref(_M0L8_2afieldS2915);
        #line 112 "/home/developer/Documents/1/fasta_io.mbt"
        moonbit_free(_M0L3recS981);
      }
      _M0L3seqS2398 = _M0L8_2afieldS2708;
      #line 112 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L3seqS983 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2398);
      moonbit_decref(_M0L3seqS2398);
      _M0L1nS984 = Moonbit_array_length(_M0L3seqS983);
      _M0L1iS985
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1iS985)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1iS985->$0 = 0;
      while (1) {
        int32_t _M0L3valS2390 = _M0L1iS985->$0;
        if (_M0L3valS2390 < _M0L1nS984) {
          int32_t _M0L3valS2396 = _M0L1iS985->$0;
          int32_t _M0L6_2atmpS2395 =
            _M0L3valS2396 + _M0FP26biolab8bio__seq18fasta__line__width;
          int32_t _M0L3endS986;
          int32_t _M0L3valS2393;
          int64_t _M0L6_2atmpS2394;
          struct _M0TPC16string10StringView _M0L6_2atmpS2392;
          moonbit_string_t _M0L6_2atmpS2391;
          if (_M0L6_2atmpS2395 < _M0L1nS984) {
            int32_t _M0L3valS2397 = _M0L1iS985->$0;
            _M0L3endS986
            = _M0L3valS2397 + _M0FP26biolab8bio__seq18fasta__line__width;
          } else {
            _M0L3endS986 = _M0L1nS984;
          }
          _M0L3valS2393 = _M0L1iS985->$0;
          _M0L6_2atmpS2394 = (int64_t)_M0L3endS986;
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2392
          = _M0MPC16string6String11sub_2einner(_M0L3seqS983, _M0L3valS2393, _M0L6_2atmpS2394);
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2391
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2392);
          moonbit_decref(_M0L6_2atmpS2392.$0);
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS977, _M0L6_2atmpS2391);
          moonbit_decref(_M0L6_2atmpS2391);
          #line 118 "/home/developer/Documents/1/fasta_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS977, 10);
          _M0L1iS985->$0 = _M0L3endS986;
          continue;
        } else {
          moonbit_decref(_M0L1iS985);
          moonbit_decref(_M0L3seqS983);
        }
        break;
      }
      _M0L6_2atmpS2401 = _M0L2__S980 + 1;
      _M0L2__S980 = _M0L6_2atmpS2401;
      continue;
    }
    break;
  }
  #line 122 "/home/developer/Documents/1/fasta_io.mbt"
  _result_3009 = _M0MPB13StringBuilder10to__string(_M0L3bufS977);
  moonbit_decref(_M0L3bufS977);
  return _result_3009;
}

moonbit_string_t _M0FP26biolab8bio__seq19build__fasta__title(
  moonbit_string_t _M0L2idS976,
  moonbit_string_t _M0L11descriptionS975
) {
  int32_t _M0L6_2atmpS2386;
  #line 129 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2386 = Moonbit_array_length(_M0L11descriptionS975);
  if (_M0L6_2atmpS2386 == 0) {
    moonbit_incref(_M0L2idS976);
    return _M0L2idS976;
  } else {
    int32_t _M0L6_2atmpS2388 = Moonbit_array_length(_M0L2idS976);
    struct _M0TPC16string10StringView _M0L6_2atmpS2387;
    int32_t _result_3010;
    moonbit_incref(_M0L2idS976);
    _M0L6_2atmpS2387
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L2idS976, .$1 = 0, .$2 = _M0L6_2atmpS2388
    };
    #line 132 "/home/developer/Documents/1/fasta_io.mbt"
    _result_3010
    = _M0MPC16string6String11has__prefix(_M0L11descriptionS975, _M0L6_2atmpS2387);
    moonbit_decref(_M0L6_2atmpS2387.$0);
    if (_result_3010) {
      moonbit_incref(_M0L11descriptionS975);
      return _M0L11descriptionS975;
    } else {
      moonbit_string_t _M0L6_2atmpS2389;
      moonbit_string_t _result_3011;
      #line 135 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L6_2atmpS2389
      = moonbit_add_string(_M0L2idS976, (moonbit_string_t)moonbit_string_literal_72.data);
      #line 135 "/home/developer/Documents/1/fasta_io.mbt"
      _result_3011
      = moonbit_add_string(_M0L6_2atmpS2389, _M0L11descriptionS975);
      moonbit_decref(_M0L6_2atmpS2389);
      return _result_3011;
    }
  }
}

struct moonbit_result_2 _M0FP26biolab8bio__seq11seqio__read(
  moonbit_string_t _M0L7contentS973,
  moonbit_string_t _M0L6formatS974
) {
  struct moonbit_result_1 _tmp_3012;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS972;
  int32_t _M0L6_2atmpS2379;
  int32_t _M0L6_2atmpS2381;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2383;
  struct moonbit_result_2 _result_3016;
  #line 32 "/home/developer/Documents/1/seqio.mbt"
  #line 36 "/home/developer/Documents/1/seqio.mbt"
  _tmp_3012
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L7contentS973, _M0L6formatS974);
  if (_tmp_3012.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2384 =
      _tmp_3012.data.ok;
    _M0L7recordsS972 = _M0L5_2aokS2384;
  } else {
    void* const _M0L6_2aerrS2385 = _tmp_3012.data.err;
    struct moonbit_result_2 _result_3013;
    _result_3013.tag = 0;
    _result_3013.data.err = _M0L6_2aerrS2385;
    return _result_3013;
  }
  #line 37 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS2379
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS972);
  if (_M0L6_2atmpS2379 == 0) {
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2380;
    struct moonbit_result_2 _result_3014;
    moonbit_decref(_M0L7recordsS972);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2380
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2380)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2380)->$0
    = (moonbit_string_t)moonbit_string_literal_73.data;
    _result_3014.tag = 0;
    _result_3014.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2380;
    return _result_3014;
  }
  #line 40 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS2381
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS972);
  if (_M0L6_2atmpS2381 > 1) {
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2382;
    struct moonbit_result_2 _result_3015;
    moonbit_decref(_M0L7recordsS972);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2382
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2382)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2382)->$0
    = (moonbit_string_t)moonbit_string_literal_74.data;
    _result_3015.tag = 0;
    _result_3015.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2382;
    return _result_3015;
  }
  #line 43 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS2383
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS972, 0);
  moonbit_decref(_M0L7recordsS972);
  _result_3016.tag = 1;
  _result_3016.data.ok = _M0L6_2atmpS2383;
  return _result_3016;
}

struct moonbit_result_1 _M0FP26biolab8bio__seq12seqio__parse(
  moonbit_string_t _M0L7contentS969,
  moonbit_string_t _M0L6formatS971
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2375;
  struct moonbit_result_1 _result_3021;
  #line 15 "/home/developer/Documents/1/seqio.mbt"
  if (
    _M0L6formatS971 == (moonbit_string_t)moonbit_string_literal_16.data
    || Moonbit_array_length(_M0L6formatS971) == 5
       && 0
          == memcmp(_M0L6formatS971, (moonbit_string_t)moonbit_string_literal_16.data, 10)
  ) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2376;
    struct moonbit_result_1 _result_3020;
    #line 20 "/home/developer/Documents/1/seqio.mbt"
    _M0L6_2atmpS2376 = _M0FP26biolab8bio__seq12parse__fasta(_M0L7contentS969);
    _result_3020.tag = 1;
    _result_3020.data.ok = _M0L6_2atmpS2376;
    return _result_3020;
  } else if (
           _M0L6formatS971
           == (moonbit_string_t)moonbit_string_literal_35.data
           || Moonbit_array_length(_M0L6formatS971) == 5
              && 0
                 == memcmp(_M0L6formatS971, (moonbit_string_t)moonbit_string_literal_35.data, 10)
         ) {
    goto join_970;
  } else if (
           _M0L6formatS971
           == (moonbit_string_t)moonbit_string_literal_66.data
           || Moonbit_array_length(_M0L6formatS971) == 12
              && 0
                 == memcmp(_M0L6formatS971, (moonbit_string_t)moonbit_string_literal_66.data, 24)
         ) {
    goto join_970;
  } else if (
           _M0L6formatS971
           == (moonbit_string_t)moonbit_string_literal_44.data
           || Moonbit_array_length(_M0L6formatS971) == 7
              && 0
                 == memcmp(_M0L6formatS971, (moonbit_string_t)moonbit_string_literal_44.data, 14)
         ) {
    goto join_968;
  } else if (
           _M0L6formatS971
           == (moonbit_string_t)moonbit_string_literal_75.data
           || Moonbit_array_length(_M0L6formatS971) == 2
              && 0
                 == memcmp(_M0L6formatS971, (moonbit_string_t)moonbit_string_literal_75.data, 4)
         ) {
    goto join_968;
  } else {
    moonbit_string_t _M0L6_2atmpS2378;
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2377;
    struct moonbit_result_1 _result_3019;
    #line 23 "/home/developer/Documents/1/seqio.mbt"
    _M0L6_2atmpS2378
    = moonbit_add_string((moonbit_string_t)moonbit_string_literal_65.data, _M0L6formatS971);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2377
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2377)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2377)->$0
    = _M0L6_2atmpS2378;
    _result_3019.tag = 0;
    _result_3019.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2377;
    return _result_3019;
  }
  join_970:;
  #line 21 "/home/developer/Documents/1/seqio.mbt"
  return _M0FP26biolab8bio__seq12parse__fastq(_M0L7contentS969);
  join_968:;
  #line 22 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS2375 = _M0FP26biolab8bio__seq14parse__genbank(_M0L7contentS969);
  _result_3021.tag = 1;
  _result_3021.data.ok = _M0L6_2atmpS2375;
  return _result_3021;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq14parse__genbank(
  moonbit_string_t _M0L7contentS913
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS911;
  moonbit_string_t _M0L7_2abindS914;
  int32_t _M0L6_2atmpS2374;
  struct _M0TPC16string10StringView _M0L6_2atmpS2373;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2372;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS912;
  int32_t _M0L1nS915;
  struct _M0TPB8MutLocalGiE* _M0L1iS916;
  int32_t _M0L6_2atmpS2371;
  int32_t _M0L2crS917;
  int32_t _M0L6_2atmpS2370;
  int32_t _M0L2spS918;
  int32_t _M0L6_2atmpS2369;
  int32_t _M0L3dotS919;
  #line 41 "/home/developer/Documents/1/genbank_io.mbt"
  #line 42 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L7recordsS911
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS914 = (moonbit_string_t)moonbit_string_literal_76.data;
  _M0L6_2atmpS2374 = Moonbit_array_length(_M0L7_2abindS914);
  _M0L6_2atmpS2373
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS914, .$1 = 0, .$2 = _M0L6_2atmpS2374
  };
  #line 43 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L6_2atmpS2372
  = _M0MPC16string6String5split(_M0L7contentS913, _M0L6_2atmpS2373);
  moonbit_decref(_M0L6_2atmpS2373.$0);
  #line 43 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L5linesS912
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2372);
  moonbit_decref(_M0L6_2atmpS2372);
  #line 44 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L1nS915
  = _M0MPC15array5Array6lengthGRPC16string10StringViewE(_M0L5linesS912);
  _M0L1iS916
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS916)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS916->$0 = 0;
  _M0L6_2atmpS2371 = 13;
  _M0L2crS917 = (uint16_t)_M0L6_2atmpS2371;
  _M0L6_2atmpS2370 = 32;
  _M0L2spS918 = (uint16_t)_M0L6_2atmpS2370;
  _M0L6_2atmpS2369 = 46;
  _M0L3dotS919 = (uint16_t)_M0L6_2atmpS2369;
  while (1) {
    int32_t _M0L3valS2128 = _M0L1iS916->$0;
    if (_M0L3valS2128 < _M0L1nS915) {
      int32_t _M0L3valS2368 = _M0L1iS916->$0;
      struct _M0TPC16string10StringView _M0L4lineS920;
      int32_t _M0L3lenS921;
      struct _M0TPB8MutLocalGiE* _M0L7trimmedS922;
      int32_t _if__result_3023;
      int32_t _M0L3valS2147;
      int32_t _if__result_3024;
      int32_t _M0L3valS2367;
      int32_t _M0L6_2atmpS2366;
      #line 51 "/home/developer/Documents/1/genbank_io.mbt"
      _M0L4lineS920
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS912, _M0L3valS2368);
      #line 52 "/home/developer/Documents/1/genbank_io.mbt"
      _M0L3lenS921 = _M0MPC16string10StringView6length(_M0L4lineS920);
      _M0L7trimmedS922
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L7trimmedS922)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L7trimmedS922->$0 = _M0L3lenS921;
      if (_M0L3lenS921 > 0) {
        int32_t _M0L6_2atmpS2130 = _M0L3lenS921 - 1;
        int32_t _M0L6_2atmpS2129;
        #line 54 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2129
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS920, _M0L6_2atmpS2130);
        _if__result_3023 = _M0L6_2atmpS2129 == _M0L2crS917;
      } else {
        _if__result_3023 = 0;
      }
      if (_if__result_3023) {
        int32_t _M0L6_2atmpS2131 = _M0L3lenS921 - 1;
        _M0L7trimmedS922->$0 = _M0L6_2atmpS2131;
      }
      _M0L3valS2147 = _M0L7trimmedS922->$0;
      if (_M0L3valS2147 >= 5) {
        int32_t _M0L6_2atmpS2144;
        int32_t _M0L6_2atmpS2146;
        int32_t _M0L6_2atmpS2145;
        #line 58 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2144
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS920, 0);
        _M0L6_2atmpS2146 = 76;
        _M0L6_2atmpS2145 = (uint16_t)_M0L6_2atmpS2146;
        if (_M0L6_2atmpS2144 == _M0L6_2atmpS2145) {
          int32_t _M0L6_2atmpS2141;
          int32_t _M0L6_2atmpS2143;
          int32_t _M0L6_2atmpS2142;
          #line 59 "/home/developer/Documents/1/genbank_io.mbt"
          _M0L6_2atmpS2141
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS920, 1);
          _M0L6_2atmpS2143 = 79;
          _M0L6_2atmpS2142 = (uint16_t)_M0L6_2atmpS2143;
          if (_M0L6_2atmpS2141 == _M0L6_2atmpS2142) {
            int32_t _M0L6_2atmpS2138;
            int32_t _M0L6_2atmpS2140;
            int32_t _M0L6_2atmpS2139;
            #line 60 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L6_2atmpS2138
            = _M0MPC16string10StringView11unsafe__get(_M0L4lineS920, 2);
            _M0L6_2atmpS2140 = 67;
            _M0L6_2atmpS2139 = (uint16_t)_M0L6_2atmpS2140;
            if (_M0L6_2atmpS2138 == _M0L6_2atmpS2139) {
              int32_t _M0L6_2atmpS2135;
              int32_t _M0L6_2atmpS2137;
              int32_t _M0L6_2atmpS2136;
              #line 61 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2135
              = _M0MPC16string10StringView11unsafe__get(_M0L4lineS920, 3);
              _M0L6_2atmpS2137 = 85;
              _M0L6_2atmpS2136 = (uint16_t)_M0L6_2atmpS2137;
              if (_M0L6_2atmpS2135 == _M0L6_2atmpS2136) {
                int32_t _M0L6_2atmpS2132;
                int32_t _M0L6_2atmpS2134;
                int32_t _M0L6_2atmpS2133;
                #line 62 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2132
                = _M0MPC16string10StringView11unsafe__get(_M0L4lineS920, 4);
                _M0L6_2atmpS2134 = 83;
                _M0L6_2atmpS2133 = (uint16_t)_M0L6_2atmpS2134;
                _if__result_3024 = _M0L6_2atmpS2132 == _M0L6_2atmpS2133;
              } else {
                _if__result_3024 = 0;
              }
            } else {
              _if__result_3024 = 0;
            }
          } else {
            _if__result_3024 = 0;
          }
        } else {
          _if__result_3024 = 0;
        }
      } else {
        _if__result_3024 = 0;
      }
      if (_if__result_3024) {
        struct _M0TPB8MutLocalGsE* _M0L4nameS923 =
          (struct _M0TPB8MutLocalGsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGsE));
        struct _M0TPB8MutLocalGsE* _M0L2idS924;
        struct _M0TPB8MutLocalGsE* _M0L11descriptionS925;
        uint16_t* _M0L6_2atmpS2365;
        struct _M0TPB8MutLocalGAkE* _M0L8seq__bufS926;
        struct _M0TPB8MutLocalGiE* _M0L8seq__lenS927;
        struct _M0TPB8MutLocalGbE* _M0L10in__originS928;
        int32_t _M0L3valS2364;
        int64_t _M0L6_2atmpS2363;
        struct _M0TPC16string10StringView _M0L11locus__restS929;
        struct _M0TPB8MutLocalGiE* _M0L1jS930;
        int32_t _M0L9rest__lenS931;
        int32_t _M0L3valS2362;
        struct _M0TPB8MutLocalGiE* _M0L9name__endS933;
        int32_t _M0L3valS2160;
        int32_t _M0L3valS2162;
        int64_t _M0L6_2atmpS2161;
        struct _M0TPC16string10StringView _M0L6_2atmpS2159;
        moonbit_string_t _M0L6_2atmpS2158;
        moonbit_string_t _M0L6_2aoldS2730;
        int32_t _M0L3valS2164;
        int32_t _M0L6_2atmpS2163;
        int32_t _M0L3valS2361;
        uint16_t* _M0L15final__seq__bufS963;
        int32_t _M0L1kS964;
        moonbit_string_t _M0L6_2atmpS2360;
        struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2353;
        moonbit_string_t _M0L8_2afieldS2715;
        int32_t _M0L6_2acntS2920;
        moonbit_string_t _M0L3valS2354;
        moonbit_string_t _M0L8_2afieldS2714;
        int32_t _M0L6_2acntS2922;
        moonbit_string_t _M0L3valS2355;
        moonbit_string_t _M0L8_2afieldS2713;
        int32_t _M0L6_2acntS2924;
        moonbit_string_t _M0L3valS2356;
        struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS966;
        struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2359;
        struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2358;
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2357;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2352;
        Moonbit_object_header(_M0L4nameS923)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 34, 0);
        _M0L4nameS923->$0 = (moonbit_string_t)moonbit_string_literal_77.data;
        _M0L2idS924
        = (struct _M0TPB8MutLocalGsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGsE));
        Moonbit_object_header(_M0L2idS924)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 34, 0);
        _M0L2idS924->$0 = (moonbit_string_t)moonbit_string_literal_62.data;
        _M0L11descriptionS925
        = (struct _M0TPB8MutLocalGsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGsE));
        Moonbit_object_header(_M0L11descriptionS925)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 34, 0);
        _M0L11descriptionS925->$0
        = (moonbit_string_t)moonbit_string_literal_64.data;
        _M0L6_2atmpS2365 = (uint16_t*)moonbit_make_string(0, 0);
        _M0L8seq__bufS926
        = (struct _M0TPB8MutLocalGAkE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGAkE));
        Moonbit_object_header(_M0L8seq__bufS926)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 0);
        _M0L8seq__bufS926->$0 = _M0L6_2atmpS2365;
        _M0L8seq__lenS927
        = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
        Moonbit_object_header(_M0L8seq__lenS927)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L8seq__lenS927->$0 = 0;
        _M0L10in__originS928
        = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
        Moonbit_object_header(_M0L10in__originS928)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L10in__originS928->$0 = 0;
        _M0L3valS2364 = _M0L7trimmedS922->$0;
        moonbit_decref(_M0L7trimmedS922);
        _M0L6_2atmpS2363 = (int64_t)_M0L3valS2364;
        #line 70 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L11locus__restS929
        = _M0MPC16string10StringView11sub_2einner(_M0L4lineS920, 5, _M0L6_2atmpS2363);
        moonbit_decref(_M0L4lineS920.$0);
        _M0L1jS930
        = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
        Moonbit_object_header(_M0L1jS930)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L1jS930->$0 = 0;
        #line 72 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L9rest__lenS931
        = _M0MPC16string10StringView6length(_M0L11locus__restS929);
        while (1) {
          int32_t _M0L3valS2150 = _M0L1jS930->$0;
          int32_t _if__result_3026;
          if (_M0L3valS2150 < _M0L9rest__lenS931) {
            int32_t _M0L3valS2149 = _M0L1jS930->$0;
            int32_t _M0L6_2atmpS2148;
            #line 73 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L6_2atmpS2148
            = _M0MPC16string10StringView11unsafe__get(_M0L11locus__restS929, _M0L3valS2149);
            _if__result_3026 = _M0L6_2atmpS2148 == _M0L2spS918;
          } else {
            _if__result_3026 = 0;
          }
          if (_if__result_3026) {
            int32_t _M0L3valS2152 = _M0L1jS930->$0;
            int32_t _M0L6_2atmpS2151 = _M0L3valS2152 + 1;
            _M0L1jS930->$0 = _M0L6_2atmpS2151;
            continue;
          }
          break;
        }
        _M0L3valS2362 = _M0L1jS930->$0;
        _M0L9name__endS933
        = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
        Moonbit_object_header(_M0L9name__endS933)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L9name__endS933->$0 = _M0L3valS2362;
        while (1) {
          int32_t _M0L3valS2155 = _M0L9name__endS933->$0;
          int32_t _if__result_3028;
          if (_M0L3valS2155 < _M0L9rest__lenS931) {
            int32_t _M0L3valS2154 = _M0L9name__endS933->$0;
            int32_t _M0L6_2atmpS2153;
            #line 77 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L6_2atmpS2153
            = _M0MPC16string10StringView11unsafe__get(_M0L11locus__restS929, _M0L3valS2154);
            _if__result_3028 = _M0L6_2atmpS2153 != _M0L2spS918;
          } else {
            _if__result_3028 = 0;
          }
          if (_if__result_3028) {
            int32_t _M0L3valS2157 = _M0L9name__endS933->$0;
            int32_t _M0L6_2atmpS2156 = _M0L3valS2157 + 1;
            _M0L9name__endS933->$0 = _M0L6_2atmpS2156;
            continue;
          }
          break;
        }
        _M0L3valS2160 = _M0L1jS930->$0;
        moonbit_decref(_M0L1jS930);
        _M0L3valS2162 = _M0L9name__endS933->$0;
        moonbit_decref(_M0L9name__endS933);
        _M0L6_2atmpS2161 = (int64_t)_M0L3valS2162;
        #line 80 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2159
        = _M0MPC16string10StringView11sub_2einner(_M0L11locus__restS929, _M0L3valS2160, _M0L6_2atmpS2161);
        moonbit_decref(_M0L11locus__restS929.$0);
        #line 80 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2158
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2159);
        moonbit_decref(_M0L6_2atmpS2159.$0);
        _M0L6_2aoldS2730 = _M0L4nameS923->$0;
        moonbit_decref(_M0L6_2aoldS2730);
        _M0L4nameS923->$0 = _M0L6_2atmpS2158;
        _M0L3valS2164 = _M0L1iS916->$0;
        _M0L6_2atmpS2163 = _M0L3valS2164 + 1;
        _M0L1iS916->$0 = _M0L6_2atmpS2163;
        while (1) {
          int32_t _M0L3valS2165 = _M0L1iS916->$0;
          if (_M0L3valS2165 < _M0L1nS915) {
            int32_t _M0L3valS2347 = _M0L1iS916->$0;
            struct _M0TPC16string10StringView _M0L1lS935;
            int32_t _M0L2llS936;
            struct _M0TPB8MutLocalGiE* _M0L2ltS937;
            int32_t _if__result_3030;
            int32_t _M0L3valS2175;
            int32_t _if__result_3031;
            int32_t _M0L3valS2206;
            int32_t _if__result_3032;
            int32_t _M0L3valS2276;
            int32_t _if__result_3041;
            int32_t _M0L3valS2313;
            int32_t _if__result_3046;
            int32_t _M0L3valS2346;
            int32_t _M0L6_2atmpS2345;
            #line 84 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L1lS935
            = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS912, _M0L3valS2347);
            #line 85 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L2llS936 = _M0MPC16string10StringView6length(_M0L1lS935);
            _M0L2ltS937
            = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
            Moonbit_object_header(_M0L2ltS937)->meta
            = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
            _M0L2ltS937->$0 = _M0L2llS936;
            if (_M0L2llS936 > 0) {
              int32_t _M0L6_2atmpS2167 = _M0L2llS936 - 1;
              int32_t _M0L6_2atmpS2166;
              #line 87 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2166
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, _M0L6_2atmpS2167);
              _if__result_3030 = _M0L6_2atmpS2166 == _M0L2crS917;
            } else {
              _if__result_3030 = 0;
            }
            if (_if__result_3030) {
              int32_t _M0L6_2atmpS2168 = _M0L2llS936 - 1;
              _M0L2ltS937->$0 = _M0L6_2atmpS2168;
            }
            _M0L3valS2175 = _M0L2ltS937->$0;
            if (_M0L3valS2175 >= 2) {
              int32_t _M0L6_2atmpS2172;
              int32_t _M0L6_2atmpS2174;
              int32_t _M0L6_2atmpS2173;
              #line 91 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2172
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 0);
              _M0L6_2atmpS2174 = 47;
              _M0L6_2atmpS2173 = (uint16_t)_M0L6_2atmpS2174;
              if (_M0L6_2atmpS2172 == _M0L6_2atmpS2173) {
                int32_t _M0L6_2atmpS2169;
                int32_t _M0L6_2atmpS2171;
                int32_t _M0L6_2atmpS2170;
                #line 92 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2169
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 1);
                _M0L6_2atmpS2171 = 47;
                _M0L6_2atmpS2170 = (uint16_t)_M0L6_2atmpS2171;
                _if__result_3031 = _M0L6_2atmpS2169 == _M0L6_2atmpS2170;
              } else {
                _if__result_3031 = 0;
              }
            } else {
              _if__result_3031 = 0;
            }
            if (_if__result_3031) {
              moonbit_decref(_M0L2ltS937);
              moonbit_decref(_M0L1lS935.$0);
              moonbit_decref(_M0L10in__originS928);
              break;
            }
            _M0L3valS2206 = _M0L2ltS937->$0;
            if (_M0L3valS2206 >= 10) {
              int32_t _M0L6_2atmpS2203;
              int32_t _M0L6_2atmpS2205;
              int32_t _M0L6_2atmpS2204;
              #line 96 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2203
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 0);
              _M0L6_2atmpS2205 = 68;
              _M0L6_2atmpS2204 = (uint16_t)_M0L6_2atmpS2205;
              if (_M0L6_2atmpS2203 == _M0L6_2atmpS2204) {
                int32_t _M0L6_2atmpS2200;
                int32_t _M0L6_2atmpS2202;
                int32_t _M0L6_2atmpS2201;
                #line 97 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2200
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 1);
                _M0L6_2atmpS2202 = 69;
                _M0L6_2atmpS2201 = (uint16_t)_M0L6_2atmpS2202;
                if (_M0L6_2atmpS2200 == _M0L6_2atmpS2201) {
                  int32_t _M0L6_2atmpS2197;
                  int32_t _M0L6_2atmpS2199;
                  int32_t _M0L6_2atmpS2198;
                  #line 98 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2197
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 2);
                  _M0L6_2atmpS2199 = 70;
                  _M0L6_2atmpS2198 = (uint16_t)_M0L6_2atmpS2199;
                  if (_M0L6_2atmpS2197 == _M0L6_2atmpS2198) {
                    int32_t _M0L6_2atmpS2194;
                    int32_t _M0L6_2atmpS2196;
                    int32_t _M0L6_2atmpS2195;
                    #line 99 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS2194
                    = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 3);
                    _M0L6_2atmpS2196 = 73;
                    _M0L6_2atmpS2195 = (uint16_t)_M0L6_2atmpS2196;
                    if (_M0L6_2atmpS2194 == _M0L6_2atmpS2195) {
                      int32_t _M0L6_2atmpS2191;
                      int32_t _M0L6_2atmpS2193;
                      int32_t _M0L6_2atmpS2192;
                      #line 100 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS2191
                      = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 4);
                      _M0L6_2atmpS2193 = 78;
                      _M0L6_2atmpS2192 = (uint16_t)_M0L6_2atmpS2193;
                      if (_M0L6_2atmpS2191 == _M0L6_2atmpS2192) {
                        int32_t _M0L6_2atmpS2188;
                        int32_t _M0L6_2atmpS2190;
                        int32_t _M0L6_2atmpS2189;
                        #line 101 "/home/developer/Documents/1/genbank_io.mbt"
                        _M0L6_2atmpS2188
                        = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 5);
                        _M0L6_2atmpS2190 = 73;
                        _M0L6_2atmpS2189 = (uint16_t)_M0L6_2atmpS2190;
                        if (_M0L6_2atmpS2188 == _M0L6_2atmpS2189) {
                          int32_t _M0L6_2atmpS2185;
                          int32_t _M0L6_2atmpS2187;
                          int32_t _M0L6_2atmpS2186;
                          #line 102 "/home/developer/Documents/1/genbank_io.mbt"
                          _M0L6_2atmpS2185
                          = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 6);
                          _M0L6_2atmpS2187 = 84;
                          _M0L6_2atmpS2186 = (uint16_t)_M0L6_2atmpS2187;
                          if (_M0L6_2atmpS2185 == _M0L6_2atmpS2186) {
                            int32_t _M0L6_2atmpS2182;
                            int32_t _M0L6_2atmpS2184;
                            int32_t _M0L6_2atmpS2183;
                            #line 103 "/home/developer/Documents/1/genbank_io.mbt"
                            _M0L6_2atmpS2182
                            = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 7);
                            _M0L6_2atmpS2184 = 73;
                            _M0L6_2atmpS2183 = (uint16_t)_M0L6_2atmpS2184;
                            if (_M0L6_2atmpS2182 == _M0L6_2atmpS2183) {
                              int32_t _M0L6_2atmpS2179;
                              int32_t _M0L6_2atmpS2181;
                              int32_t _M0L6_2atmpS2180;
                              #line 104 "/home/developer/Documents/1/genbank_io.mbt"
                              _M0L6_2atmpS2179
                              = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 8);
                              _M0L6_2atmpS2181 = 79;
                              _M0L6_2atmpS2180 = (uint16_t)_M0L6_2atmpS2181;
                              if (_M0L6_2atmpS2179 == _M0L6_2atmpS2180) {
                                int32_t _M0L6_2atmpS2176;
                                int32_t _M0L6_2atmpS2178;
                                int32_t _M0L6_2atmpS2177;
                                #line 105 "/home/developer/Documents/1/genbank_io.mbt"
                                _M0L6_2atmpS2176
                                = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 9);
                                _M0L6_2atmpS2178 = 78;
                                _M0L6_2atmpS2177 = (uint16_t)_M0L6_2atmpS2178;
                                _if__result_3032
                                = _M0L6_2atmpS2176 == _M0L6_2atmpS2177;
                              } else {
                                _if__result_3032 = 0;
                              }
                            } else {
                              _if__result_3032 = 0;
                            }
                          } else {
                            _if__result_3032 = 0;
                          }
                        } else {
                          _if__result_3032 = 0;
                        }
                      } else {
                        _if__result_3032 = 0;
                      }
                    } else {
                      _if__result_3032 = 0;
                    }
                  } else {
                    _if__result_3032 = 0;
                  }
                } else {
                  _if__result_3032 = 0;
                }
              } else {
                _if__result_3032 = 0;
              }
            } else {
              _if__result_3032 = 0;
            }
            if (_if__result_3032) {
              int32_t _M0L3valS2254 = _M0L2ltS937->$0;
              int64_t _M0L6_2atmpS2253;
              struct _M0TPC16string10StringView _M0L9def__partS939;
              struct _M0TPB8MutLocalGiE* _M0L2djS940;
              int32_t _M0L3dplS941;
              struct _M0TPB13StringBuilder* _M0L8def__bufS943;
              int32_t _M0L3valS2214;
              int64_t _M0L6_2atmpS2215;
              struct _M0TPC16string10StringView _M0L6_2atmpS2213;
              moonbit_string_t _M0L6_2atmpS2212;
              int32_t _M0L3valS2217;
              int32_t _M0L6_2atmpS2216;
              moonbit_string_t _M0L6_2atmpS2238;
              moonbit_string_t _M0L6_2aoldS2729;
              moonbit_string_t _M0L3valS2245;
              int32_t _M0L6_2atmpS2244;
              int32_t _if__result_3040;
              moonbit_decref(_M0L2ltS937);
              _M0L6_2atmpS2253 = (int64_t)_M0L3valS2254;
              #line 106 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L9def__partS939
              = _M0MPC16string10StringView11sub_2einner(_M0L1lS935, 10, _M0L6_2atmpS2253);
              moonbit_decref(_M0L1lS935.$0);
              _M0L2djS940
              = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
              Moonbit_object_header(_M0L2djS940)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
              _M0L2djS940->$0 = 0;
              #line 108 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L3dplS941
              = _M0MPC16string10StringView6length(_M0L9def__partS939);
              while (1) {
                int32_t _M0L3valS2209 = _M0L2djS940->$0;
                int32_t _if__result_3034;
                if (_M0L3valS2209 < _M0L3dplS941) {
                  int32_t _M0L3valS2208 = _M0L2djS940->$0;
                  int32_t _M0L6_2atmpS2207;
                  #line 109 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2207
                  = _M0MPC16string10StringView11unsafe__get(_M0L9def__partS939, _M0L3valS2208);
                  _if__result_3034 = _M0L6_2atmpS2207 == _M0L2spS918;
                } else {
                  _if__result_3034 = 0;
                }
                if (_if__result_3034) {
                  int32_t _M0L3valS2211 = _M0L2djS940->$0;
                  int32_t _M0L6_2atmpS2210 = _M0L3valS2211 + 1;
                  _M0L2djS940->$0 = _M0L6_2atmpS2210;
                  continue;
                }
                break;
              }
              #line 112 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L8def__bufS943
              = _M0MPB13StringBuilder21StringBuilder_2einner(0);
              _M0L3valS2214 = _M0L2djS940->$0;
              moonbit_decref(_M0L2djS940);
              _M0L6_2atmpS2215 = (int64_t)_M0L3dplS941;
              #line 113 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2213
              = _M0MPC16string10StringView11sub_2einner(_M0L9def__partS939, _M0L3valS2214, _M0L6_2atmpS2215);
              moonbit_decref(_M0L9def__partS939.$0);
              #line 113 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2212
              = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2213);
              moonbit_decref(_M0L6_2atmpS2213.$0);
              #line 113 "/home/developer/Documents/1/genbank_io.mbt"
              _M0IPB13StringBuilderPB6Logger13write__string(_M0L8def__bufS943, _M0L6_2atmpS2212);
              moonbit_decref(_M0L6_2atmpS2212);
              _M0L3valS2217 = _M0L1iS916->$0;
              _M0L6_2atmpS2216 = _M0L3valS2217 + 1;
              _M0L1iS916->$0 = _M0L6_2atmpS2216;
              while (1) {
                int32_t _M0L3valS2218 = _M0L1iS916->$0;
                if (_M0L3valS2218 < _M0L1nS915) {
                  int32_t _M0L3valS2237 = _M0L1iS916->$0;
                  struct _M0TPC16string10StringView _M0L4contS944;
                  int32_t _M0L2clS945;
                  struct _M0TPB8MutLocalGiE* _M0L2ctS946;
                  int32_t _if__result_3036;
                  int32_t _M0L3valS2223;
                  int32_t _if__result_3037;
                  struct _M0TPB8MutLocalGiE* _M0L2cjS948;
                  int32_t _M0L3valS2232;
                  int32_t _M0L3valS2234;
                  int64_t _M0L6_2atmpS2233;
                  struct _M0TPC16string10StringView _M0L6_2atmpS2231;
                  moonbit_string_t _M0L6_2atmpS2230;
                  int32_t _M0L3valS2236;
                  int32_t _M0L6_2atmpS2235;
                  #line 116 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L4contS944
                  = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS912, _M0L3valS2237);
                  #line 117 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L2clS945
                  = _M0MPC16string10StringView6length(_M0L4contS944);
                  _M0L2ctS946
                  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
                  Moonbit_object_header(_M0L2ctS946)->meta
                  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
                  _M0L2ctS946->$0 = _M0L2clS945;
                  if (_M0L2clS945 > 0) {
                    int32_t _M0L6_2atmpS2220 = _M0L2clS945 - 1;
                    int32_t _M0L6_2atmpS2219;
                    #line 119 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS2219
                    = _M0MPC16string10StringView11unsafe__get(_M0L4contS944, _M0L6_2atmpS2220);
                    _if__result_3036 = _M0L6_2atmpS2219 == _M0L2crS917;
                  } else {
                    _if__result_3036 = 0;
                  }
                  if (_if__result_3036) {
                    int32_t _M0L6_2atmpS2221 = _M0L2clS945 - 1;
                    _M0L2ctS946->$0 = _M0L6_2atmpS2221;
                  }
                  _M0L3valS2223 = _M0L2ctS946->$0;
                  if (_M0L3valS2223 == 0) {
                    _if__result_3037 = 1;
                  } else {
                    int32_t _M0L6_2atmpS2222;
                    #line 122 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS2222
                    = _M0MPC16string10StringView11unsafe__get(_M0L4contS944, 0);
                    _if__result_3037 = _M0L6_2atmpS2222 != _M0L2spS918;
                  }
                  if (_if__result_3037) {
                    moonbit_decref(_M0L2ctS946);
                    moonbit_decref(_M0L4contS944.$0);
                    break;
                  }
                  #line 125 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0IPB13StringBuilderPB6Logger11write__char(_M0L8def__bufS943, 32);
                  _M0L2cjS948
                  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
                  Moonbit_object_header(_M0L2cjS948)->meta
                  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
                  _M0L2cjS948->$0 = 0;
                  while (1) {
                    int32_t _M0L3valS2226 = _M0L2cjS948->$0;
                    int32_t _M0L3valS2227 = _M0L2ctS946->$0;
                    int32_t _if__result_3039;
                    if (_M0L3valS2226 < _M0L3valS2227) {
                      int32_t _M0L3valS2225 = _M0L2cjS948->$0;
                      int32_t _M0L6_2atmpS2224;
                      #line 127 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS2224
                      = _M0MPC16string10StringView11unsafe__get(_M0L4contS944, _M0L3valS2225);
                      _if__result_3039 = _M0L6_2atmpS2224 == _M0L2spS918;
                    } else {
                      _if__result_3039 = 0;
                    }
                    if (_if__result_3039) {
                      int32_t _M0L3valS2229 = _M0L2cjS948->$0;
                      int32_t _M0L6_2atmpS2228 = _M0L3valS2229 + 1;
                      _M0L2cjS948->$0 = _M0L6_2atmpS2228;
                      continue;
                    }
                    break;
                  }
                  _M0L3valS2232 = _M0L2cjS948->$0;
                  moonbit_decref(_M0L2cjS948);
                  _M0L3valS2234 = _M0L2ctS946->$0;
                  moonbit_decref(_M0L2ctS946);
                  _M0L6_2atmpS2233 = (int64_t)_M0L3valS2234;
                  #line 130 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2231
                  = _M0MPC16string10StringView11sub_2einner(_M0L4contS944, _M0L3valS2232, _M0L6_2atmpS2233);
                  moonbit_decref(_M0L4contS944.$0);
                  #line 130 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2230
                  = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2231);
                  moonbit_decref(_M0L6_2atmpS2231.$0);
                  #line 130 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0IPB13StringBuilderPB6Logger13write__string(_M0L8def__bufS943, _M0L6_2atmpS2230);
                  moonbit_decref(_M0L6_2atmpS2230);
                  _M0L3valS2236 = _M0L1iS916->$0;
                  _M0L6_2atmpS2235 = _M0L3valS2236 + 1;
                  _M0L1iS916->$0 = _M0L6_2atmpS2235;
                  continue;
                }
                break;
              }
              #line 133 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2238
              = _M0MPB13StringBuilder10to__string(_M0L8def__bufS943);
              moonbit_decref(_M0L8def__bufS943);
              _M0L6_2aoldS2729 = _M0L11descriptionS925->$0;
              moonbit_decref(_M0L6_2aoldS2729);
              _M0L11descriptionS925->$0 = _M0L6_2atmpS2238;
              _M0L3valS2245 = _M0L11descriptionS925->$0;
              _M0L6_2atmpS2244 = Moonbit_array_length(_M0L3valS2245);
              if (_M0L6_2atmpS2244 > 0) {
                moonbit_string_t _M0L3valS2240 = _M0L11descriptionS925->$0;
                moonbit_string_t _M0L3valS2243 = _M0L11descriptionS925->$0;
                int32_t _M0L6_2atmpS2242 =
                  Moonbit_array_length(_M0L3valS2243);
                int32_t _M0L6_2atmpS2241 = _M0L6_2atmpS2242 - 1;
                int32_t _M0L6_2atmpS2239 = _M0L3valS2240[_M0L6_2atmpS2241];
                _if__result_3040 = _M0L6_2atmpS2239 == _M0L3dotS919;
              } else {
                _if__result_3040 = 0;
              }
              if (_if__result_3040) {
                moonbit_string_t _M0L3valS2248 = _M0L11descriptionS925->$0;
                moonbit_string_t _M0L3valS2252 = _M0L11descriptionS925->$0;
                int32_t _M0L6_2atmpS2251 =
                  Moonbit_array_length(_M0L3valS2252);
                int32_t _M0L6_2atmpS2250 = _M0L6_2atmpS2251 - 1;
                int64_t _M0L6_2atmpS2249 = (int64_t)_M0L6_2atmpS2250;
                struct _M0TPC16string10StringView _M0L6_2atmpS2247;
                moonbit_string_t _M0L6_2atmpS2246;
                moonbit_string_t _M0L6_2aoldS2723;
                moonbit_incref(_M0L3valS2248);
                #line 136 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2247
                = _M0MPC16string6String11sub_2einner(_M0L3valS2248, 0, _M0L6_2atmpS2249);
                moonbit_decref(_M0L3valS2248);
                #line 136 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2246
                = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2247);
                moonbit_decref(_M0L6_2atmpS2247.$0);
                _M0L6_2aoldS2723 = _M0L11descriptionS925->$0;
                moonbit_decref(_M0L6_2aoldS2723);
                _M0L11descriptionS925->$0 = _M0L6_2atmpS2246;
              }
              continue;
            }
            _M0L3valS2276 = _M0L2ltS937->$0;
            if (_M0L3valS2276 >= 7) {
              int32_t _M0L6_2atmpS2273;
              int32_t _M0L6_2atmpS2275;
              int32_t _M0L6_2atmpS2274;
              #line 141 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2273
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 0);
              _M0L6_2atmpS2275 = 86;
              _M0L6_2atmpS2274 = (uint16_t)_M0L6_2atmpS2275;
              if (_M0L6_2atmpS2273 == _M0L6_2atmpS2274) {
                int32_t _M0L6_2atmpS2270;
                int32_t _M0L6_2atmpS2272;
                int32_t _M0L6_2atmpS2271;
                #line 142 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2270
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 1);
                _M0L6_2atmpS2272 = 69;
                _M0L6_2atmpS2271 = (uint16_t)_M0L6_2atmpS2272;
                if (_M0L6_2atmpS2270 == _M0L6_2atmpS2271) {
                  int32_t _M0L6_2atmpS2267;
                  int32_t _M0L6_2atmpS2269;
                  int32_t _M0L6_2atmpS2268;
                  #line 143 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2267
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 2);
                  _M0L6_2atmpS2269 = 82;
                  _M0L6_2atmpS2268 = (uint16_t)_M0L6_2atmpS2269;
                  if (_M0L6_2atmpS2267 == _M0L6_2atmpS2268) {
                    int32_t _M0L6_2atmpS2264;
                    int32_t _M0L6_2atmpS2266;
                    int32_t _M0L6_2atmpS2265;
                    #line 144 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS2264
                    = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 3);
                    _M0L6_2atmpS2266 = 83;
                    _M0L6_2atmpS2265 = (uint16_t)_M0L6_2atmpS2266;
                    if (_M0L6_2atmpS2264 == _M0L6_2atmpS2265) {
                      int32_t _M0L6_2atmpS2261;
                      int32_t _M0L6_2atmpS2263;
                      int32_t _M0L6_2atmpS2262;
                      #line 145 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS2261
                      = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 4);
                      _M0L6_2atmpS2263 = 73;
                      _M0L6_2atmpS2262 = (uint16_t)_M0L6_2atmpS2263;
                      if (_M0L6_2atmpS2261 == _M0L6_2atmpS2262) {
                        int32_t _M0L6_2atmpS2258;
                        int32_t _M0L6_2atmpS2260;
                        int32_t _M0L6_2atmpS2259;
                        #line 146 "/home/developer/Documents/1/genbank_io.mbt"
                        _M0L6_2atmpS2258
                        = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 5);
                        _M0L6_2atmpS2260 = 79;
                        _M0L6_2atmpS2259 = (uint16_t)_M0L6_2atmpS2260;
                        if (_M0L6_2atmpS2258 == _M0L6_2atmpS2259) {
                          int32_t _M0L6_2atmpS2255;
                          int32_t _M0L6_2atmpS2257;
                          int32_t _M0L6_2atmpS2256;
                          #line 147 "/home/developer/Documents/1/genbank_io.mbt"
                          _M0L6_2atmpS2255
                          = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 6);
                          _M0L6_2atmpS2257 = 78;
                          _M0L6_2atmpS2256 = (uint16_t)_M0L6_2atmpS2257;
                          _if__result_3041
                          = _M0L6_2atmpS2255 == _M0L6_2atmpS2256;
                        } else {
                          _if__result_3041 = 0;
                        }
                      } else {
                        _if__result_3041 = 0;
                      }
                    } else {
                      _if__result_3041 = 0;
                    }
                  } else {
                    _if__result_3041 = 0;
                  }
                } else {
                  _if__result_3041 = 0;
                }
              } else {
                _if__result_3041 = 0;
              }
            } else {
              _if__result_3041 = 0;
            }
            if (_if__result_3041) {
              int32_t _M0L3valS2294 = _M0L2ltS937->$0;
              int64_t _M0L6_2atmpS2293 = (int64_t)_M0L3valS2294;
              struct _M0TPC16string10StringView _M0L9ver__partS950;
              struct _M0TPB8MutLocalGiE* _M0L2vjS951;
              int32_t _M0L3vplS952;
              int32_t _M0L3valS2292;
              struct _M0TPB8MutLocalGiE* _M0L8vid__endS954;
              int32_t _M0L3valS2289;
              int32_t _M0L3valS2291;
              int64_t _M0L6_2atmpS2290;
              struct _M0TPC16string10StringView _M0L6_2atmpS2288;
              moonbit_string_t _M0L6_2atmpS2287;
              moonbit_string_t _M0L6_2aoldS2722;
              #line 148 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L9ver__partS950
              = _M0MPC16string10StringView11sub_2einner(_M0L1lS935, 7, _M0L6_2atmpS2293);
              _M0L2vjS951
              = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
              Moonbit_object_header(_M0L2vjS951)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
              _M0L2vjS951->$0 = 0;
              #line 150 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L3vplS952
              = _M0MPC16string10StringView6length(_M0L9ver__partS950);
              while (1) {
                int32_t _M0L3valS2279 = _M0L2vjS951->$0;
                int32_t _if__result_3043;
                if (_M0L3valS2279 < _M0L3vplS952) {
                  int32_t _M0L3valS2278 = _M0L2vjS951->$0;
                  int32_t _M0L6_2atmpS2277;
                  #line 151 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2277
                  = _M0MPC16string10StringView11unsafe__get(_M0L9ver__partS950, _M0L3valS2278);
                  _if__result_3043 = _M0L6_2atmpS2277 == _M0L2spS918;
                } else {
                  _if__result_3043 = 0;
                }
                if (_if__result_3043) {
                  int32_t _M0L3valS2281 = _M0L2vjS951->$0;
                  int32_t _M0L6_2atmpS2280 = _M0L3valS2281 + 1;
                  _M0L2vjS951->$0 = _M0L6_2atmpS2280;
                  continue;
                }
                break;
              }
              _M0L3valS2292 = _M0L2vjS951->$0;
              _M0L8vid__endS954
              = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
              Moonbit_object_header(_M0L8vid__endS954)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
              _M0L8vid__endS954->$0 = _M0L3valS2292;
              while (1) {
                int32_t _M0L3valS2284 = _M0L8vid__endS954->$0;
                int32_t _if__result_3045;
                if (_M0L3valS2284 < _M0L3vplS952) {
                  int32_t _M0L3valS2283 = _M0L8vid__endS954->$0;
                  int32_t _M0L6_2atmpS2282;
                  #line 155 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2282
                  = _M0MPC16string10StringView11unsafe__get(_M0L9ver__partS950, _M0L3valS2283);
                  _if__result_3045 = _M0L6_2atmpS2282 != _M0L2spS918;
                } else {
                  _if__result_3045 = 0;
                }
                if (_if__result_3045) {
                  int32_t _M0L3valS2286 = _M0L8vid__endS954->$0;
                  int32_t _M0L6_2atmpS2285 = _M0L3valS2286 + 1;
                  _M0L8vid__endS954->$0 = _M0L6_2atmpS2285;
                  continue;
                }
                break;
              }
              _M0L3valS2289 = _M0L2vjS951->$0;
              moonbit_decref(_M0L2vjS951);
              _M0L3valS2291 = _M0L8vid__endS954->$0;
              moonbit_decref(_M0L8vid__endS954);
              _M0L6_2atmpS2290 = (int64_t)_M0L3valS2291;
              #line 158 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2288
              = _M0MPC16string10StringView11sub_2einner(_M0L9ver__partS950, _M0L3valS2289, _M0L6_2atmpS2290);
              moonbit_decref(_M0L9ver__partS950.$0);
              #line 158 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2287
              = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2288);
              moonbit_decref(_M0L6_2atmpS2288.$0);
              _M0L6_2aoldS2722 = _M0L2idS924->$0;
              moonbit_decref(_M0L6_2aoldS2722);
              _M0L2idS924->$0 = _M0L6_2atmpS2287;
            }
            _M0L3valS2313 = _M0L2ltS937->$0;
            if (_M0L3valS2313 >= 6) {
              int32_t _M0L6_2atmpS2310;
              int32_t _M0L6_2atmpS2312;
              int32_t _M0L6_2atmpS2311;
              #line 161 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2310
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 0);
              _M0L6_2atmpS2312 = 79;
              _M0L6_2atmpS2311 = (uint16_t)_M0L6_2atmpS2312;
              if (_M0L6_2atmpS2310 == _M0L6_2atmpS2311) {
                int32_t _M0L6_2atmpS2307;
                int32_t _M0L6_2atmpS2309;
                int32_t _M0L6_2atmpS2308;
                #line 162 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2307
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 1);
                _M0L6_2atmpS2309 = 82;
                _M0L6_2atmpS2308 = (uint16_t)_M0L6_2atmpS2309;
                if (_M0L6_2atmpS2307 == _M0L6_2atmpS2308) {
                  int32_t _M0L6_2atmpS2304;
                  int32_t _M0L6_2atmpS2306;
                  int32_t _M0L6_2atmpS2305;
                  #line 163 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS2304
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 2);
                  _M0L6_2atmpS2306 = 73;
                  _M0L6_2atmpS2305 = (uint16_t)_M0L6_2atmpS2306;
                  if (_M0L6_2atmpS2304 == _M0L6_2atmpS2305) {
                    int32_t _M0L6_2atmpS2301;
                    int32_t _M0L6_2atmpS2303;
                    int32_t _M0L6_2atmpS2302;
                    #line 164 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS2301
                    = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 3);
                    _M0L6_2atmpS2303 = 71;
                    _M0L6_2atmpS2302 = (uint16_t)_M0L6_2atmpS2303;
                    if (_M0L6_2atmpS2301 == _M0L6_2atmpS2302) {
                      int32_t _M0L6_2atmpS2298;
                      int32_t _M0L6_2atmpS2300;
                      int32_t _M0L6_2atmpS2299;
                      #line 165 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS2298
                      = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 4);
                      _M0L6_2atmpS2300 = 73;
                      _M0L6_2atmpS2299 = (uint16_t)_M0L6_2atmpS2300;
                      if (_M0L6_2atmpS2298 == _M0L6_2atmpS2299) {
                        int32_t _M0L6_2atmpS2295;
                        int32_t _M0L6_2atmpS2297;
                        int32_t _M0L6_2atmpS2296;
                        #line 166 "/home/developer/Documents/1/genbank_io.mbt"
                        _M0L6_2atmpS2295
                        = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, 5);
                        _M0L6_2atmpS2297 = 78;
                        _M0L6_2atmpS2296 = (uint16_t)_M0L6_2atmpS2297;
                        _if__result_3046
                        = _M0L6_2atmpS2295 == _M0L6_2atmpS2296;
                      } else {
                        _if__result_3046 = 0;
                      }
                    } else {
                      _if__result_3046 = 0;
                    }
                  } else {
                    _if__result_3046 = 0;
                  }
                } else {
                  _if__result_3046 = 0;
                }
              } else {
                _if__result_3046 = 0;
              }
            } else {
              _if__result_3046 = 0;
            }
            if (_if__result_3046) {
              int32_t _M0L3valS2315;
              int32_t _M0L6_2atmpS2314;
              moonbit_decref(_M0L2ltS937);
              moonbit_decref(_M0L1lS935.$0);
              _M0L10in__originS928->$0 = 1;
              _M0L3valS2315 = _M0L1iS916->$0;
              _M0L6_2atmpS2314 = _M0L3valS2315 + 1;
              _M0L1iS916->$0 = _M0L6_2atmpS2314;
              continue;
            }
            if (_M0L10in__originS928->$0) {
              int32_t _M0L3valS2343 = _M0L8seq__lenS927->$0;
              int32_t _M0L3valS2344 = _M0L2ltS937->$0;
              int32_t _M0L8est__capS956 = _M0L3valS2343 + _M0L3valS2344;
              uint16_t* _M0L3valS2317 = _M0L8seq__bufS926->$0;
              int32_t _M0L6_2atmpS2316 = Moonbit_array_length(_M0L3valS2317);
              int32_t _M0L1jS960;
              if (_M0L6_2atmpS2316 < _M0L8est__capS956) {
                int32_t _M0L6_2atmpS2322 = _M0L8est__capS956 * 2;
                uint16_t* _M0L8new__bufS957 =
                  (uint16_t*)moonbit_make_string(_M0L6_2atmpS2322, 0);
                int32_t _M0L1kS958 = 0;
                uint16_t* _M0L6_2aoldS2719;
                while (1) {
                  int32_t _M0L3valS2318 = _M0L8seq__lenS927->$0;
                  if (_M0L1kS958 < _M0L3valS2318) {
                    uint16_t* _M0L3valS2320 = _M0L8seq__bufS926->$0;
                    int32_t _M0L6_2atmpS2319;
                    int32_t _M0L6_2atmpS2321;
                    if (
                      _M0L1kS958 < 0
                      || _M0L1kS958 >= Moonbit_array_length(_M0L3valS2320)
                    ) {
                      #line 176 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L6_2atmpS2319 = (int32_t)_M0L3valS2320[_M0L1kS958];
                    if (
                      _M0L1kS958 < 0
                      || _M0L1kS958
                         >= Moonbit_array_length(_M0L8new__bufS957)
                    ) {
                      #line 176 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L8new__bufS957[_M0L1kS958] = _M0L6_2atmpS2319;
                    _M0L6_2atmpS2321 = _M0L1kS958 + 1;
                    _M0L1kS958 = _M0L6_2atmpS2321;
                    continue;
                  }
                  break;
                }
                _M0L6_2aoldS2719 = _M0L8seq__bufS926->$0;
                moonbit_decref(_M0L6_2aoldS2719);
                _M0L8seq__bufS926->$0 = _M0L8new__bufS957;
              }
              _M0L1jS960 = 0;
              while (1) {
                int32_t _M0L3valS2323 = _M0L2ltS937->$0;
                if (_M0L1jS960 < _M0L3valS2323) {
                  int32_t _M0L1cS961;
                  int32_t _M0L6_2atmpS2327;
                  int32_t _M0L6_2atmpS2326;
                  int32_t _if__result_3049;
                  int32_t _M0L6_2atmpS2342;
                  #line 181 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L1cS961
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS935, _M0L1jS960);
                  _M0L6_2atmpS2327 = 97;
                  _M0L6_2atmpS2326 = (uint16_t)_M0L6_2atmpS2327;
                  if (_M0L1cS961 >= _M0L6_2atmpS2326) {
                    int32_t _M0L6_2atmpS2325 = 122;
                    int32_t _M0L6_2atmpS2324 = (uint16_t)_M0L6_2atmpS2325;
                    _if__result_3049 = _M0L1cS961 <= _M0L6_2atmpS2324;
                  } else {
                    _if__result_3049 = 0;
                  }
                  if (_if__result_3049) {
                    uint16_t* _M0L3valS2328 = _M0L8seq__bufS926->$0;
                    int32_t _M0L3valS2329 = _M0L8seq__lenS927->$0;
                    int32_t _M0L6_2atmpS2331 = (int32_t)_M0L1cS961;
                    int32_t _M0L6_2atmpS2330;
                    int32_t _M0L3valS2333;
                    int32_t _M0L6_2atmpS2332;
                    if (
                      _M0L6_2atmpS2331 < 0
                      || _M0L6_2atmpS2331
                         >= Moonbit_array_length(_M0FP26biolab8bio__seq16uppercase__table)
                    ) {
                      #line 183 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L6_2atmpS2330
                    = (int32_t)_M0FP26biolab8bio__seq16uppercase__table[
                        _M0L6_2atmpS2331
                      ];
                    if (
                      _M0L3valS2329 < 0
                      || _M0L3valS2329 >= Moonbit_array_length(_M0L3valS2328)
                    ) {
                      #line 183 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L3valS2328[_M0L3valS2329] = _M0L6_2atmpS2330;
                    _M0L3valS2333 = _M0L8seq__lenS927->$0;
                    _M0L6_2atmpS2332 = _M0L3valS2333 + 1;
                    _M0L8seq__lenS927->$0 = _M0L6_2atmpS2332;
                  } else {
                    int32_t _M0L6_2atmpS2337 = 65;
                    int32_t _M0L6_2atmpS2336 = (uint16_t)_M0L6_2atmpS2337;
                    int32_t _if__result_3050;
                    if (_M0L1cS961 >= _M0L6_2atmpS2336) {
                      int32_t _M0L6_2atmpS2335 = 90;
                      int32_t _M0L6_2atmpS2334 = (uint16_t)_M0L6_2atmpS2335;
                      _if__result_3050 = _M0L1cS961 <= _M0L6_2atmpS2334;
                    } else {
                      _if__result_3050 = 0;
                    }
                    if (_if__result_3050) {
                      uint16_t* _M0L3valS2338 = _M0L8seq__bufS926->$0;
                      int32_t _M0L3valS2339 = _M0L8seq__lenS927->$0;
                      int32_t _M0L3valS2341;
                      int32_t _M0L6_2atmpS2340;
                      if (
                        _M0L3valS2339 < 0
                        || _M0L3valS2339
                           >= Moonbit_array_length(_M0L3valS2338)
                      ) {
                        #line 187 "/home/developer/Documents/1/genbank_io.mbt"
                        moonbit_panic();
                      }
                      _M0L3valS2338[_M0L3valS2339] = _M0L1cS961;
                      _M0L3valS2341 = _M0L8seq__lenS927->$0;
                      _M0L6_2atmpS2340 = _M0L3valS2341 + 1;
                      _M0L8seq__lenS927->$0 = _M0L6_2atmpS2340;
                    }
                  }
                  _M0L6_2atmpS2342 = _M0L1jS960 + 1;
                  _M0L1jS960 = _M0L6_2atmpS2342;
                  continue;
                } else {
                  moonbit_decref(_M0L2ltS937);
                  moonbit_decref(_M0L1lS935.$0);
                }
                break;
              }
            } else {
              moonbit_decref(_M0L2ltS937);
              moonbit_decref(_M0L1lS935.$0);
            }
            _M0L3valS2346 = _M0L1iS916->$0;
            _M0L6_2atmpS2345 = _M0L3valS2346 + 1;
            _M0L1iS916->$0 = _M0L6_2atmpS2345;
            continue;
          } else {
            moonbit_decref(_M0L10in__originS928);
          }
          break;
        }
        _M0L3valS2361 = _M0L8seq__lenS927->$0;
        _M0L15final__seq__bufS963
        = (uint16_t*)moonbit_make_string(_M0L3valS2361, 0);
        _M0L1kS964 = 0;
        while (1) {
          int32_t _M0L3valS2348 = _M0L8seq__lenS927->$0;
          if (_M0L1kS964 < _M0L3valS2348) {
            uint16_t* _M0L3valS2350 = _M0L8seq__bufS926->$0;
            int32_t _M0L6_2atmpS2349;
            int32_t _M0L6_2atmpS2351;
            if (
              _M0L1kS964 < 0
              || _M0L1kS964 >= Moonbit_array_length(_M0L3valS2350)
            ) {
              #line 196 "/home/developer/Documents/1/genbank_io.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2349 = (int32_t)_M0L3valS2350[_M0L1kS964];
            if (
              _M0L1kS964 < 0
              || _M0L1kS964
                 >= Moonbit_array_length(_M0L15final__seq__bufS963)
            ) {
              #line 196 "/home/developer/Documents/1/genbank_io.mbt"
              moonbit_panic();
            }
            _M0L15final__seq__bufS963[_M0L1kS964] = _M0L6_2atmpS2349;
            _M0L6_2atmpS2351 = _M0L1kS964 + 1;
            _M0L1kS964 = _M0L6_2atmpS2351;
            continue;
          } else {
            moonbit_decref(_M0L8seq__lenS927);
            moonbit_decref(_M0L8seq__bufS926);
          }
          break;
        }
        _M0L6_2atmpS2360 = _M0L15final__seq__bufS963;
        #line 199 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2353 = _M0MP26biolab8bio__seq3Seq3new(_M0L6_2atmpS2360);
        moonbit_decref(_M0L6_2atmpS2360);
        _M0L8_2afieldS2715 = _M0L2idS924->$0;
        _M0L6_2acntS2920
        = Moonbit_rc_count(Moonbit_object_header(_M0L2idS924));
        if (_M0L6_2acntS2920 > 1) {
          int32_t _M0L11_2anew__cntS2921 = _M0L6_2acntS2920 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L2idS924), _M0L11_2anew__cntS2921);
          moonbit_incref(_M0L8_2afieldS2715);
        } else if (_M0L6_2acntS2920 == 1) {
          #line 198 "/home/developer/Documents/1/genbank_io.mbt"
          moonbit_free(_M0L2idS924);
        }
        _M0L3valS2354 = _M0L8_2afieldS2715;
        _M0L8_2afieldS2714 = _M0L4nameS923->$0;
        _M0L6_2acntS2922
        = Moonbit_rc_count(Moonbit_object_header(_M0L4nameS923));
        if (_M0L6_2acntS2922 > 1) {
          int32_t _M0L11_2anew__cntS2923 = _M0L6_2acntS2922 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L4nameS923), _M0L11_2anew__cntS2923);
          moonbit_incref(_M0L8_2afieldS2714);
        } else if (_M0L6_2acntS2922 == 1) {
          #line 198 "/home/developer/Documents/1/genbank_io.mbt"
          moonbit_free(_M0L4nameS923);
        }
        _M0L3valS2355 = _M0L8_2afieldS2714;
        _M0L8_2afieldS2713 = _M0L11descriptionS925->$0;
        _M0L6_2acntS2924
        = Moonbit_rc_count(Moonbit_object_header(_M0L11descriptionS925));
        if (_M0L6_2acntS2924 > 1) {
          int32_t _M0L11_2anew__cntS2925 = _M0L6_2acntS2924 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L11descriptionS925), _M0L11_2anew__cntS2925);
          moonbit_incref(_M0L8_2afieldS2713);
        } else if (_M0L6_2acntS2924 == 1) {
          #line 198 "/home/developer/Documents/1/genbank_io.mbt"
          moonbit_free(_M0L11descriptionS925);
        }
        _M0L3valS2356 = _M0L8_2afieldS2713;
        _M0L7_2abindS966
        = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
        _M0L6_2atmpS2359 = _M0L7_2abindS966;
        _M0L6_2atmpS2358
        = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
          .$0 = _M0L6_2atmpS2359, .$1 = 0, .$2 = 0
        };
        #line 203 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2357
        = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2358, 0ll);
        moonbit_decref(_M0L6_2atmpS2358.$0);
        _M0L6_2atmpS2352
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
        Moonbit_object_header(_M0L6_2atmpS2352)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
        _M0L6_2atmpS2352->$0 = _M0L6_2atmpS2353;
        _M0L6_2atmpS2352->$1 = _M0L3valS2354;
        _M0L6_2atmpS2352->$2 = _M0L3valS2355;
        _M0L6_2atmpS2352->$3 = _M0L3valS2356;
        _M0L6_2atmpS2352->$4 = _M0L6_2atmpS2357;
        #line 198 "/home/developer/Documents/1/genbank_io.mbt"
        _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS911, _M0L6_2atmpS2352);
        moonbit_decref(_M0L6_2atmpS2352);
      } else {
        moonbit_decref(_M0L7trimmedS922);
        moonbit_decref(_M0L4lineS920.$0);
      }
      _M0L3valS2367 = _M0L1iS916->$0;
      _M0L6_2atmpS2366 = _M0L3valS2367 + 1;
      _M0L1iS916->$0 = _M0L6_2atmpS2366;
      continue;
    } else {
      moonbit_decref(_M0L1iS916);
      moonbit_decref(_M0L5linesS912);
    }
    break;
  }
  return _M0L7recordsS911;
}

struct moonbit_result_1 _M0FP26biolab8bio__seq12parse__fastq(
  moonbit_string_t _M0L7contentS867
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS865;
  moonbit_string_t _M0L7_2abindS868;
  int32_t _M0L6_2atmpS2127;
  struct _M0TPC16string10StringView _M0L6_2atmpS2126;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2125;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS866;
  int32_t _M0L1nS869;
  struct _M0TPB8MutLocalGiE* _M0L1iS870;
  int32_t _M0L6_2atmpS2124;
  int32_t _M0L2atS871;
  int32_t _M0L6_2atmpS2123;
  int32_t _M0L4plusS872;
  int32_t _M0L6_2atmpS2122;
  int32_t _M0L2crS873;
  struct moonbit_result_1 _result_3067;
  #line 22 "/home/developer/Documents/1/fastq_io.mbt"
  #line 23 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L7recordsS865
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS868 = (moonbit_string_t)moonbit_string_literal_76.data;
  _M0L6_2atmpS2127 = Moonbit_array_length(_M0L7_2abindS868);
  _M0L6_2atmpS2126
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS868, .$1 = 0, .$2 = _M0L6_2atmpS2127
  };
  #line 24 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L6_2atmpS2125
  = _M0MPC16string6String5split(_M0L7contentS867, _M0L6_2atmpS2126);
  moonbit_decref(_M0L6_2atmpS2126.$0);
  #line 24 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L5linesS866
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2125);
  moonbit_decref(_M0L6_2atmpS2125);
  #line 25 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L1nS869
  = _M0MPC15array5Array6lengthGRPC16string10StringViewE(_M0L5linesS866);
  _M0L1iS870
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS870)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS870->$0 = 0;
  _M0L6_2atmpS2124 = 64;
  _M0L2atS871 = (uint16_t)_M0L6_2atmpS2124;
  _M0L6_2atmpS2123 = 43;
  _M0L4plusS872 = (uint16_t)_M0L6_2atmpS2123;
  _M0L6_2atmpS2122 = 13;
  _M0L2crS873 = (uint16_t)_M0L6_2atmpS2122;
  while (1) {
    int32_t _M0L3valS2052 = _M0L1iS870->$0;
    int32_t _M0L6_2atmpS2051 = _M0L3valS2052 + 3;
    if (_M0L6_2atmpS2051 < _M0L1nS869) {
      int32_t _M0L3valS2121 = _M0L1iS870->$0;
      struct _M0TPC16string10StringView _M0L5line0S874;
      int32_t _M0L3valS2120;
      int32_t _M0L6_2atmpS2119;
      struct _M0TPC16string10StringView _M0L5line1S875;
      int32_t _M0L3valS2118;
      int32_t _M0L6_2atmpS2117;
      struct _M0TPC16string10StringView _M0L5line2S876;
      int32_t _M0L3valS2116;
      int32_t _M0L6_2atmpS2115;
      struct _M0TPC16string10StringView _M0L5line3S877;
      int32_t _M0L4len0S878;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed0S879;
      int32_t _if__result_3053;
      int32_t _M0L3valS2058;
      int32_t _if__result_3054;
      int32_t _M0L3valS2060;
      int32_t _if__result_3055;
      int32_t _M0L4len2S882;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed2S883;
      int32_t _if__result_3057;
      int32_t _M0L3valS2067;
      int32_t _if__result_3058;
      int32_t _M0L4len1S885;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed1S886;
      int32_t _if__result_3060;
      int32_t _M0L4len3S887;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed3S888;
      int32_t _if__result_3061;
      int32_t _M0L3valS2076;
      int32_t _M0L3valS2077;
      int32_t _M0L3valS2111;
      moonbit_string_t _M0L5titleS889;
      moonbit_string_t _M0L2idS891;
      moonbit_string_t _M0L11descriptionS892;
      int32_t _M0L3posS903;
      moonbit_string_t _M0L7_2abindS905;
      int32_t _M0L6_2atmpS2110;
      struct _M0TPC16string10StringView _M0L6_2atmpS2109;
      int64_t _M0L7_2abindS904;
      int64_t _M0L6_2atmpS2108;
      struct _M0TPC16string10StringView _M0L6_2atmpS2107;
      moonbit_string_t _M0L6_2atmpS2106;
      int32_t _M0L3valS2102;
      moonbit_string_t _M0L8seq__strS893;
      int32_t _M0L3valS2101;
      int32_t* _M0L5qualsS894;
      int32_t _M0L1jS895;
      struct _M0TPB5ArrayGiE* _M0L10quals__arrS897;
      int32_t _M0L1jS898;
      struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS901;
      struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2100;
      struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2099;
      struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L11letter__annS900;
      struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2096;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2095;
      int32_t _M0L3valS2098;
      int32_t _M0L6_2atmpS2097;
      #line 32 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line0S874
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS866, _M0L3valS2121);
      _M0L3valS2120 = _M0L1iS870->$0;
      _M0L6_2atmpS2119 = _M0L3valS2120 + 1;
      #line 33 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line1S875
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS866, _M0L6_2atmpS2119);
      _M0L3valS2118 = _M0L1iS870->$0;
      _M0L6_2atmpS2117 = _M0L3valS2118 + 2;
      #line 34 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line2S876
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS866, _M0L6_2atmpS2117);
      _M0L3valS2116 = _M0L1iS870->$0;
      _M0L6_2atmpS2115 = _M0L3valS2116 + 3;
      #line 35 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line3S877
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS866, _M0L6_2atmpS2115);
      #line 37 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len0S878 = _M0MPC16string10StringView6length(_M0L5line0S874);
      _M0L8trimmed0S879
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed0S879)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed0S879->$0 = _M0L4len0S878;
      if (_M0L4len0S878 > 0) {
        int32_t _M0L6_2atmpS2054 = _M0L4len0S878 - 1;
        int32_t _M0L6_2atmpS2053;
        #line 39 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2053
        = _M0MPC16string10StringView11unsafe__get(_M0L5line0S874, _M0L6_2atmpS2054);
        _if__result_3053 = _M0L6_2atmpS2053 == _M0L2crS873;
      } else {
        _if__result_3053 = 0;
      }
      if (_if__result_3053) {
        int32_t _M0L6_2atmpS2055 = _M0L4len0S878 - 1;
        _M0L8trimmed0S879->$0 = _M0L6_2atmpS2055;
      }
      _M0L3valS2058 = _M0L8trimmed0S879->$0;
      if (_M0L3valS2058 == 0) {
        int32_t _M0L3valS2057 = _M0L1iS870->$0;
        int32_t _M0L6_2atmpS2056 = _M0L3valS2057 + 4;
        _if__result_3054 = _M0L6_2atmpS2056 >= _M0L1nS869;
      } else {
        _if__result_3054 = 0;
      }
      if (_if__result_3054) {
        moonbit_decref(_M0L8trimmed0S879);
        moonbit_decref(_M0L5line3S877.$0);
        moonbit_decref(_M0L5line2S876.$0);
        moonbit_decref(_M0L5line1S875.$0);
        moonbit_decref(_M0L5line0S874.$0);
        moonbit_decref(_M0L1iS870);
        moonbit_decref(_M0L5linesS866);
        break;
      }
      _M0L3valS2060 = _M0L8trimmed0S879->$0;
      if (_M0L3valS2060 == 0) {
        _if__result_3055 = 1;
      } else {
        int32_t _M0L6_2atmpS2059;
        #line 45 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2059
        = _M0MPC16string10StringView11unsafe__get(_M0L5line0S874, 0);
        _if__result_3055 = _M0L6_2atmpS2059 != _M0L2atS871;
      }
      if (_if__result_3055) {
        moonbit_string_t _M0L10line0__strS881;
        moonbit_string_t _M0L6_2atmpS2062;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2061;
        struct moonbit_result_1 _result_3056;
        moonbit_decref(_M0L8trimmed0S879);
        moonbit_decref(_M0L5line3S877.$0);
        moonbit_decref(_M0L5line2S876.$0);
        moonbit_decref(_M0L5line1S875.$0);
        moonbit_decref(_M0L1iS870);
        moonbit_decref(_M0L5linesS866);
        moonbit_decref(_M0L7recordsS865);
        #line 46 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L10line0__strS881
        = _M0MPC16string10StringView9to__owned(_M0L5line0S874);
        moonbit_decref(_M0L5line0S874.$0);
        #line 48 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2062
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_78.data, _M0L10line0__strS881);
        moonbit_decref(_M0L10line0__strS881);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2061
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2061)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2061)->$0
        = _M0L6_2atmpS2062;
        _result_3056.tag = 0;
        _result_3056.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2061;
        return _result_3056;
      }
      #line 52 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len2S882 = _M0MPC16string10StringView6length(_M0L5line2S876);
      _M0L8trimmed2S883
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed2S883)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed2S883->$0 = _M0L4len2S882;
      if (_M0L4len2S882 > 0) {
        int32_t _M0L6_2atmpS2064 = _M0L4len2S882 - 1;
        int32_t _M0L6_2atmpS2063;
        #line 54 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2063
        = _M0MPC16string10StringView11unsafe__get(_M0L5line2S876, _M0L6_2atmpS2064);
        _if__result_3057 = _M0L6_2atmpS2063 == _M0L2crS873;
      } else {
        _if__result_3057 = 0;
      }
      if (_if__result_3057) {
        int32_t _M0L6_2atmpS2065 = _M0L4len2S882 - 1;
        _M0L8trimmed2S883->$0 = _M0L6_2atmpS2065;
      }
      _M0L3valS2067 = _M0L8trimmed2S883->$0;
      moonbit_decref(_M0L8trimmed2S883);
      if (_M0L3valS2067 == 0) {
        _if__result_3058 = 1;
      } else {
        int32_t _M0L6_2atmpS2066;
        #line 57 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2066
        = _M0MPC16string10StringView11unsafe__get(_M0L5line2S876, 0);
        _if__result_3058 = _M0L6_2atmpS2066 != _M0L4plusS872;
      }
      if (_if__result_3058) {
        moonbit_string_t _M0L10line2__strS884;
        moonbit_string_t _M0L6_2atmpS2069;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2068;
        struct moonbit_result_1 _result_3059;
        moonbit_decref(_M0L8trimmed0S879);
        moonbit_decref(_M0L5line3S877.$0);
        moonbit_decref(_M0L5line1S875.$0);
        moonbit_decref(_M0L5line0S874.$0);
        moonbit_decref(_M0L1iS870);
        moonbit_decref(_M0L5linesS866);
        moonbit_decref(_M0L7recordsS865);
        #line 58 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L10line2__strS884
        = _M0MPC16string10StringView9to__owned(_M0L5line2S876);
        moonbit_decref(_M0L5line2S876.$0);
        #line 60 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2069
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_79.data, _M0L10line2__strS884);
        moonbit_decref(_M0L10line2__strS884);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2068
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2068)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2068)->$0
        = _M0L6_2atmpS2069;
        _result_3059.tag = 0;
        _result_3059.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2068;
        return _result_3059;
      } else {
        moonbit_decref(_M0L5line2S876.$0);
      }
      #line 64 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len1S885 = _M0MPC16string10StringView6length(_M0L5line1S875);
      _M0L8trimmed1S886
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed1S886)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed1S886->$0 = _M0L4len1S885;
      if (_M0L4len1S885 > 0) {
        int32_t _M0L6_2atmpS2071 = _M0L4len1S885 - 1;
        int32_t _M0L6_2atmpS2070;
        #line 66 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2070
        = _M0MPC16string10StringView11unsafe__get(_M0L5line1S875, _M0L6_2atmpS2071);
        _if__result_3060 = _M0L6_2atmpS2070 == _M0L2crS873;
      } else {
        _if__result_3060 = 0;
      }
      if (_if__result_3060) {
        int32_t _M0L6_2atmpS2072 = _M0L4len1S885 - 1;
        _M0L8trimmed1S886->$0 = _M0L6_2atmpS2072;
      }
      #line 70 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len3S887 = _M0MPC16string10StringView6length(_M0L5line3S877);
      _M0L8trimmed3S888
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed3S888)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed3S888->$0 = _M0L4len3S887;
      if (_M0L4len3S887 > 0) {
        int32_t _M0L6_2atmpS2074 = _M0L4len3S887 - 1;
        int32_t _M0L6_2atmpS2073;
        #line 72 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2073
        = _M0MPC16string10StringView11unsafe__get(_M0L5line3S877, _M0L6_2atmpS2074);
        _if__result_3061 = _M0L6_2atmpS2073 == _M0L2crS873;
      } else {
        _if__result_3061 = 0;
      }
      if (_if__result_3061) {
        int32_t _M0L6_2atmpS2075 = _M0L4len3S887 - 1;
        _M0L8trimmed3S888->$0 = _M0L6_2atmpS2075;
      }
      _M0L3valS2076 = _M0L8trimmed1S886->$0;
      _M0L3valS2077 = _M0L8trimmed3S888->$0;
      if (_M0L3valS2076 != _M0L3valS2077) {
        int32_t _M0L3valS2086;
        moonbit_string_t _M0L6_2atmpS2085;
        moonbit_string_t _M0L6_2atmpS2084;
        moonbit_string_t _M0L6_2atmpS2081;
        int32_t _M0L3valS2083;
        moonbit_string_t _M0L6_2atmpS2082;
        moonbit_string_t _M0L6_2atmpS2080;
        moonbit_string_t _M0L6_2atmpS2079;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2078;
        struct moonbit_result_1 _result_3062;
        moonbit_decref(_M0L8trimmed0S879);
        moonbit_decref(_M0L5line3S877.$0);
        moonbit_decref(_M0L5line1S875.$0);
        moonbit_decref(_M0L5line0S874.$0);
        moonbit_decref(_M0L1iS870);
        moonbit_decref(_M0L5linesS866);
        moonbit_decref(_M0L7recordsS865);
        _M0L3valS2086 = _M0L8trimmed1S886->$0;
        moonbit_decref(_M0L8trimmed1S886);
        #line 79 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2085
        = _M0MPC13int3Int18to__string_2einner(_M0L3valS2086, 10);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2084
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_68.data, _M0L6_2atmpS2085);
        moonbit_decref(_M0L6_2atmpS2085);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2081
        = moonbit_add_string(_M0L6_2atmpS2084, (moonbit_string_t)moonbit_string_literal_69.data);
        moonbit_decref(_M0L6_2atmpS2084);
        _M0L3valS2083 = _M0L8trimmed3S888->$0;
        moonbit_decref(_M0L8trimmed3S888);
        #line 81 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2082
        = _M0MPC13int3Int18to__string_2einner(_M0L3valS2083, 10);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2080
        = moonbit_add_string(_M0L6_2atmpS2081, _M0L6_2atmpS2082);
        moonbit_decref(_M0L6_2atmpS2082);
        moonbit_decref(_M0L6_2atmpS2081);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2079
        = moonbit_add_string(_M0L6_2atmpS2080, (moonbit_string_t)moonbit_string_literal_80.data);
        moonbit_decref(_M0L6_2atmpS2080);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2078
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2078)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2078)->$0
        = _M0L6_2atmpS2079;
        _result_3062.tag = 0;
        _result_3062.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2078;
        return _result_3062;
      }
      _M0L3valS2111 = _M0L8trimmed0S879->$0;
      if (_M0L3valS2111 > 1) {
        int32_t _M0L3valS2114 = _M0L8trimmed0S879->$0;
        int64_t _M0L6_2atmpS2113;
        struct _M0TPC16string10StringView _M0L6_2atmpS2112;
        moonbit_decref(_M0L8trimmed0S879);
        _M0L6_2atmpS2113 = (int64_t)_M0L3valS2114;
        #line 86 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2112
        = _M0MPC16string10StringView11sub_2einner(_M0L5line0S874, 1, _M0L6_2atmpS2113);
        moonbit_decref(_M0L5line0S874.$0);
        #line 86 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L5titleS889
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2112);
        moonbit_decref(_M0L6_2atmpS2112.$0);
      } else {
        moonbit_decref(_M0L8trimmed0S879);
        moonbit_decref(_M0L5line0S874.$0);
        _M0L5titleS889 = (moonbit_string_t)moonbit_string_literal_0.data;
      }
      _M0L7_2abindS905 = (moonbit_string_t)moonbit_string_literal_72.data;
      _M0L6_2atmpS2110 = Moonbit_array_length(_M0L7_2abindS905);
      _M0L6_2atmpS2109
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L7_2abindS905, .$1 = 0, .$2 = _M0L6_2atmpS2110
      };
      #line 87 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L7_2abindS904
      = _M0MPC16string6String4find(_M0L5titleS889, _M0L6_2atmpS2109);
      moonbit_decref(_M0L6_2atmpS2109.$0);
      if (_M0L7_2abindS904 == 4294967296ll) {
        moonbit_incref(_M0L5titleS889);
        _M0L2idS891 = _M0L5titleS889;
        _M0L11descriptionS892 = _M0L5titleS889;
        goto join_890;
      } else {
        int64_t _M0L7_2aSomeS906 = _M0L7_2abindS904;
        int32_t _M0L6_2aposS907 = (int32_t)_M0L7_2aSomeS906;
        _M0L3posS903 = _M0L6_2aposS907;
        goto join_902;
      }
      goto joinlet_3064;
      join_902:;
      _M0L6_2atmpS2108 = (int64_t)_M0L3posS903;
      #line 88 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2107
      = _M0MPC16string6String11sub_2einner(_M0L5titleS889, 0, _M0L6_2atmpS2108);
      #line 88 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2106
      = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2107);
      moonbit_decref(_M0L6_2atmpS2107.$0);
      _M0L2idS891 = _M0L6_2atmpS2106;
      _M0L11descriptionS892 = _M0L5titleS889;
      goto join_890;
      joinlet_3064:;
      goto joinlet_3063;
      join_890:;
      _M0L3valS2102 = _M0L8trimmed1S886->$0;
      if (_M0L3valS2102 > 0) {
        int32_t _M0L3valS2105 = _M0L8trimmed1S886->$0;
        int64_t _M0L6_2atmpS2104;
        struct _M0TPC16string10StringView _M0L6_2atmpS2103;
        moonbit_decref(_M0L8trimmed1S886);
        _M0L6_2atmpS2104 = (int64_t)_M0L3valS2105;
        #line 92 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2103
        = _M0MPC16string10StringView11sub_2einner(_M0L5line1S875, 0, _M0L6_2atmpS2104);
        moonbit_decref(_M0L5line1S875.$0);
        #line 92 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L8seq__strS893
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2103);
        moonbit_decref(_M0L6_2atmpS2103.$0);
      } else {
        moonbit_decref(_M0L8trimmed1S886);
        moonbit_decref(_M0L5line1S875.$0);
        _M0L8seq__strS893 = (moonbit_string_t)moonbit_string_literal_0.data;
      }
      _M0L3valS2101 = _M0L8trimmed3S888->$0;
      _M0L5qualsS894 = (int32_t*)moonbit_make_int32_array(_M0L3valS2101, 0);
      _M0L1jS895 = 0;
      while (1) {
        int32_t _M0L3valS2087 = _M0L8trimmed3S888->$0;
        if (_M0L1jS895 < _M0L3valS2087) {
          int32_t _M0L6_2atmpS2090;
          int32_t _M0L6_2atmpS2089;
          int32_t _M0L6_2atmpS2088;
          int32_t _M0L6_2atmpS2091;
          #line 96 "/home/developer/Documents/1/fastq_io.mbt"
          _M0L6_2atmpS2090
          = _M0MPC16string10StringView11unsafe__get(_M0L5line3S877, _M0L1jS895);
          _M0L6_2atmpS2089 = (int32_t)_M0L6_2atmpS2090;
          _M0L6_2atmpS2088 = _M0L6_2atmpS2089 - 33;
          if (
            _M0L1jS895 < 0
            || _M0L1jS895 >= Moonbit_array_length(_M0L5qualsS894)
          ) {
            #line 96 "/home/developer/Documents/1/fastq_io.mbt"
            moonbit_panic();
          }
          _M0L5qualsS894[_M0L1jS895] = _M0L6_2atmpS2088;
          _M0L6_2atmpS2091 = _M0L1jS895 + 1;
          _M0L1jS895 = _M0L6_2atmpS2091;
          continue;
        } else {
          moonbit_decref(_M0L5line3S877.$0);
        }
        break;
      }
      #line 98 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L10quals__arrS897 = _M0MPC15array5Array11new_2einnerGiE(0);
      _M0L1jS898 = 0;
      while (1) {
        int32_t _M0L3valS2092 = _M0L8trimmed3S888->$0;
        if (_M0L1jS898 < _M0L3valS2092) {
          int32_t _M0L6_2atmpS2093;
          int32_t _M0L6_2atmpS2094;
          if (
            _M0L1jS898 < 0
            || _M0L1jS898 >= Moonbit_array_length(_M0L5qualsS894)
          ) {
            #line 100 "/home/developer/Documents/1/fastq_io.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2093 = (int32_t)_M0L5qualsS894[_M0L1jS898];
          #line 100 "/home/developer/Documents/1/fastq_io.mbt"
          _M0MPC15array5Array4pushGiE(_M0L10quals__arrS897, _M0L6_2atmpS2093);
          _M0L6_2atmpS2094 = _M0L1jS898 + 1;
          _M0L1jS898 = _M0L6_2atmpS2094;
          continue;
        } else {
          moonbit_decref(_M0L5qualsS894);
          moonbit_decref(_M0L8trimmed3S888);
        }
        break;
      }
      _M0L7_2abindS901
      = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
      _M0L6_2atmpS2100 = _M0L7_2abindS901;
      _M0L6_2atmpS2099
      = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
        .$0 = _M0L6_2atmpS2100, .$1 = 0, .$2 = 0
      };
      #line 103 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L11letter__annS900
      = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2099, 1ll);
      moonbit_decref(_M0L6_2atmpS2099.$0);
      #line 104 "/home/developer/Documents/1/fastq_io.mbt"
      _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L11letter__annS900, (moonbit_string_t)moonbit_string_literal_61.data, _M0L10quals__arrS897);
      moonbit_decref(_M0L10quals__arrS897);
      #line 106 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2096 = _M0MP26biolab8bio__seq3Seq3new(_M0L8seq__strS893);
      moonbit_decref(_M0L8seq__strS893);
      moonbit_incref(_M0L2idS891);
      _M0L6_2atmpS2095
      = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
      Moonbit_object_header(_M0L6_2atmpS2095)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
      _M0L6_2atmpS2095->$0 = _M0L6_2atmpS2096;
      _M0L6_2atmpS2095->$1 = _M0L2idS891;
      _M0L6_2atmpS2095->$2 = _M0L2idS891;
      _M0L6_2atmpS2095->$3 = _M0L11descriptionS892;
      _M0L6_2atmpS2095->$4 = _M0L11letter__annS900;
      #line 105 "/home/developer/Documents/1/fastq_io.mbt"
      _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS865, _M0L6_2atmpS2095);
      moonbit_decref(_M0L6_2atmpS2095);
      _M0L3valS2098 = _M0L1iS870->$0;
      _M0L6_2atmpS2097 = _M0L3valS2098 + 4;
      _M0L1iS870->$0 = _M0L6_2atmpS2097;
      joinlet_3063:;
      continue;
    } else {
      moonbit_decref(_M0L1iS870);
      moonbit_decref(_M0L5linesS866);
    }
    break;
  }
  _result_3067.tag = 1;
  _result_3067.data.ok = _M0L7recordsS865;
  return _result_3067;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq12parse__fasta(
  moonbit_string_t _M0L7contentS837
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS835;
  moonbit_string_t _M0L7_2abindS838;
  int32_t _M0L6_2atmpS2050;
  struct _M0TPC16string10StringView _M0L6_2atmpS2049;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2048;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS836;
  moonbit_string_t _M0L6_2atmpS2047;
  struct _M0TPB8MutLocalGOsE* _M0L14current__titleS839;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS2046;
  struct _M0TPB8MutLocalGRPB13StringBuilderE* _M0L12current__seqS840;
  int32_t _M0L6_2atmpS2045;
  int32_t _M0L2gtS841;
  int32_t _M0L6_2atmpS2044;
  int32_t _M0L2crS842;
  int32_t _M0L6_2atmpS2043;
  int32_t _M0L2spS843;
  int32_t _M0L6_2atmpS2042;
  int32_t _M0L3tabS844;
  int32_t _M0L7_2abindS845;
  int32_t _M0L2__S846;
  moonbit_string_t _M0L5titleS861;
  moonbit_string_t _M0L8_2afieldS2732;
  int32_t _M0L6_2acntS2928;
  moonbit_string_t _M0L7_2abindS862;
  struct _M0TPB13StringBuilder* _M0L8_2afieldS2731;
  int32_t _M0L6_2acntS2926;
  struct _M0TPB13StringBuilder* _M0L3valS2041;
  moonbit_string_t _M0L6_2atmpS2040;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2039;
  #line 24 "/home/developer/Documents/1/fasta_io.mbt"
  #line 25 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7recordsS835
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS838 = (moonbit_string_t)moonbit_string_literal_76.data;
  _M0L6_2atmpS2050 = Moonbit_array_length(_M0L7_2abindS838);
  _M0L6_2atmpS2049
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS838, .$1 = 0, .$2 = _M0L6_2atmpS2050
  };
  #line 26 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2048
  = _M0MPC16string6String5split(_M0L7contentS837, _M0L6_2atmpS2049);
  moonbit_decref(_M0L6_2atmpS2049.$0);
  #line 26 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L5linesS836
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2048);
  moonbit_decref(_M0L6_2atmpS2048);
  _M0L6_2atmpS2047 = 0;
  _M0L14current__titleS839
  = (struct _M0TPB8MutLocalGOsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGOsE));
  Moonbit_object_header(_M0L14current__titleS839)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 40, 0);
  _M0L14current__titleS839->$0 = _M0L6_2atmpS2047;
  #line 28 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2046 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L12current__seqS840
  = (struct _M0TPB8MutLocalGRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGRPB13StringBuilderE));
  Moonbit_object_header(_M0L12current__seqS840)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 43, 0);
  _M0L12current__seqS840->$0 = _M0L6_2atmpS2046;
  _M0L6_2atmpS2045 = 62;
  _M0L2gtS841 = (uint16_t)_M0L6_2atmpS2045;
  _M0L6_2atmpS2044 = 13;
  _M0L2crS842 = (uint16_t)_M0L6_2atmpS2044;
  _M0L6_2atmpS2043 = 32;
  _M0L2spS843 = (uint16_t)_M0L6_2atmpS2043;
  _M0L6_2atmpS2042 = 9;
  _M0L3tabS844 = (uint16_t)_M0L6_2atmpS2042;
  _M0L7_2abindS845 = _M0L5linesS836->$1;
  _M0L2__S846 = 0;
  while (1) {
    if (_M0L2__S846 < _M0L7_2abindS845) {
      struct _M0TPC16string10StringView* _M0L3bufS2038 = _M0L5linesS836->$0;
      struct _M0TPC16string10StringView _M0L4lineS847 =
        _M0L3bufS2038[_M0L2__S846];
      int32_t _M0L1nS850;
      struct _M0TPB8MutLocalGiE* _M0L12trimmed__lenS851;
      int32_t _if__result_3070;
      int32_t _M0L3valS2022;
      int32_t _M0L6_2atmpS2023;
      int32_t _M0L6_2atmpS2018;
      moonbit_incref(_M0L4lineS847.$0);
      #line 35 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L1nS850 = _M0MPC16string10StringView6length(_M0L4lineS847);
      _M0L12trimmed__lenS851
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L12trimmed__lenS851)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L12trimmed__lenS851->$0 = _M0L1nS850;
      if (_M0L1nS850 > 0) {
        int32_t _M0L6_2atmpS2020 = _M0L1nS850 - 1;
        int32_t _M0L6_2atmpS2019;
        #line 37 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2019
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS847, _M0L6_2atmpS2020);
        _if__result_3070 = _M0L6_2atmpS2019 == _M0L2crS842;
      } else {
        _if__result_3070 = 0;
      }
      if (_if__result_3070) {
        int32_t _M0L6_2atmpS2021 = _M0L1nS850 - 1;
        _M0L12trimmed__lenS851->$0 = _M0L6_2atmpS2021;
      }
      _M0L3valS2022 = _M0L12trimmed__lenS851->$0;
      if (_M0L3valS2022 == 0) {
        moonbit_decref(_M0L12trimmed__lenS851);
        moonbit_decref(_M0L4lineS847.$0);
        goto join_848;
      }
      #line 43 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L6_2atmpS2023
      = _M0MPC16string10StringView11unsafe__get(_M0L4lineS847, 0);
      if (_M0L6_2atmpS2023 == _M0L2gtS841) {
        moonbit_string_t _M0L5titleS853;
        moonbit_string_t _M0L7_2abindS854 = _M0L14current__titleS839->$0;
        struct _M0TPB13StringBuilder* _M0L3valS2026;
        moonbit_string_t _M0L6_2atmpS2025;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2024;
        int32_t _M0L3valS2029;
        moonbit_string_t _M0L6_2atmpS2028;
        moonbit_string_t _M0L6_2atmpS2027;
        moonbit_string_t _M0L6_2aoldS2734;
        struct _M0TPB13StringBuilder* _M0L6_2atmpS2033;
        struct _M0TPB13StringBuilder* _M0L6_2aoldS2733;
        if (_M0L7_2abindS854 == 0) {
          
        } else {
          moonbit_string_t _M0L7_2aSomeS855 = _M0L7_2abindS854;
          moonbit_string_t _M0L8_2atitleS856 = _M0L7_2aSomeS855;
          moonbit_incref(_M0L8_2atitleS856);
          _M0L5titleS853 = _M0L8_2atitleS856;
          goto join_852;
        }
        goto joinlet_3071;
        join_852:;
        _M0L3valS2026 = _M0L12current__seqS840->$0;
        moonbit_incref(_M0L3valS2026);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2025 = _M0MPB13StringBuilder10to__string(_M0L3valS2026);
        moonbit_decref(_M0L3valS2026);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2024
        = _M0FP26biolab8bio__seq24fasta__title__to__record(_M0L5titleS853, _M0L6_2atmpS2025);
        moonbit_decref(_M0L5titleS853);
        moonbit_decref(_M0L6_2atmpS2025);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS835, _M0L6_2atmpS2024);
        moonbit_decref(_M0L6_2atmpS2024);
        joinlet_3071:;
        _M0L3valS2029 = _M0L12trimmed__lenS851->$0;
        if (_M0L3valS2029 > 1) {
          int32_t _M0L3valS2032 = _M0L12trimmed__lenS851->$0;
          int64_t _M0L6_2atmpS2031;
          struct _M0TPC16string10StringView _M0L6_2atmpS2030;
          moonbit_decref(_M0L12trimmed__lenS851);
          _M0L6_2atmpS2031 = (int64_t)_M0L3valS2032;
          #line 51 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2030
          = _M0MPC16string10StringView11sub_2einner(_M0L4lineS847, 1, _M0L6_2atmpS2031);
          moonbit_decref(_M0L4lineS847.$0);
          #line 51 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2028
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2030);
          moonbit_decref(_M0L6_2atmpS2030.$0);
        } else {
          moonbit_decref(_M0L12trimmed__lenS851);
          moonbit_decref(_M0L4lineS847.$0);
          _M0L6_2atmpS2028 = (moonbit_string_t)moonbit_string_literal_0.data;
        }
        _M0L6_2atmpS2027 = _M0L6_2atmpS2028;
        _M0L6_2aoldS2734 = _M0L14current__titleS839->$0;
        if (_M0L6_2aoldS2734) {
          moonbit_decref(_M0L6_2aoldS2734);
        }
        _M0L14current__titleS839->$0 = _M0L6_2atmpS2027;
        #line 56 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2033 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
        _M0L6_2aoldS2733 = _M0L12current__seqS840->$0;
        moonbit_decref(_M0L6_2aoldS2733);
        _M0L12current__seqS840->$0 = _M0L6_2atmpS2033;
      } else {
        int32_t _M0L1iS857 = 0;
        while (1) {
          int32_t _M0L3valS2034 = _M0L12trimmed__lenS851->$0;
          if (_M0L1iS857 < _M0L3valS2034) {
            int32_t _M0L1cS858;
            int32_t _if__result_3073;
            int32_t _M0L6_2atmpS2037;
            #line 59 "/home/developer/Documents/1/fasta_io.mbt"
            _M0L1cS858
            = _M0MPC16string10StringView11unsafe__get(_M0L4lineS847, _M0L1iS857);
            if (_M0L1cS858 != _M0L2spS843) {
              _if__result_3073 = _M0L1cS858 != _M0L3tabS844;
            } else {
              _if__result_3073 = 0;
            }
            if (_if__result_3073) {
              struct _M0TPB13StringBuilder* _M0L3valS2035 =
                _M0L12current__seqS840->$0;
              int32_t _M0L6_2atmpS2036;
              moonbit_incref(_M0L3valS2035);
              #line 61 "/home/developer/Documents/1/fasta_io.mbt"
              _M0L6_2atmpS2036
              = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS858);
              #line 61 "/home/developer/Documents/1/fasta_io.mbt"
              _M0IPB13StringBuilderPB6Logger11write__char(_M0L3valS2035, _M0L6_2atmpS2036);
              moonbit_decref(_M0L3valS2035);
            }
            _M0L6_2atmpS2037 = _M0L1iS857 + 1;
            _M0L1iS857 = _M0L6_2atmpS2037;
            continue;
          } else {
            moonbit_decref(_M0L12trimmed__lenS851);
            moonbit_decref(_M0L4lineS847.$0);
          }
          break;
        }
      }
      goto join_848;
      goto joinlet_3069;
      join_848:;
      _M0L6_2atmpS2018 = _M0L2__S846 + 1;
      _M0L2__S846 = _M0L6_2atmpS2018;
      continue;
      joinlet_3069:;
    } else {
      moonbit_decref(_M0L5linesS836);
    }
    break;
  }
  _M0L8_2afieldS2732 = _M0L14current__titleS839->$0;
  _M0L6_2acntS2928
  = Moonbit_rc_count(Moonbit_object_header(_M0L14current__titleS839));
  if (_M0L6_2acntS2928 > 1) {
    int32_t _M0L11_2anew__cntS2929 = _M0L6_2acntS2928 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L14current__titleS839), _M0L11_2anew__cntS2929);
    if (_M0L8_2afieldS2732) {
      moonbit_incref(_M0L8_2afieldS2732);
    }
  } else if (_M0L6_2acntS2928 == 1) {
    #line 66 "/home/developer/Documents/1/fasta_io.mbt"
    moonbit_free(_M0L14current__titleS839);
  }
  _M0L7_2abindS862 = _M0L8_2afieldS2732;
  if (_M0L7_2abindS862 == 0) {
    if (_M0L7_2abindS862) {
      moonbit_decref(_M0L7_2abindS862);
    }
    moonbit_decref(_M0L12current__seqS840);
  } else {
    moonbit_string_t _M0L7_2aSomeS863 = _M0L7_2abindS862;
    moonbit_string_t _M0L8_2atitleS864 = _M0L7_2aSomeS863;
    _M0L5titleS861 = _M0L8_2atitleS864;
    goto join_860;
  }
  goto joinlet_3074;
  join_860:;
  _M0L8_2afieldS2731 = _M0L12current__seqS840->$0;
  _M0L6_2acntS2926
  = Moonbit_rc_count(Moonbit_object_header(_M0L12current__seqS840));
  if (_M0L6_2acntS2926 > 1) {
    int32_t _M0L11_2anew__cntS2927 = _M0L6_2acntS2926 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L12current__seqS840), _M0L11_2anew__cntS2927);
    moonbit_incref(_M0L8_2afieldS2731);
  } else if (_M0L6_2acntS2926 == 1) {
    #line 68 "/home/developer/Documents/1/fasta_io.mbt"
    moonbit_free(_M0L12current__seqS840);
  }
  _M0L3valS2041 = _M0L8_2afieldS2731;
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2040 = _M0MPB13StringBuilder10to__string(_M0L3valS2041);
  moonbit_decref(_M0L3valS2041);
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2039
  = _M0FP26biolab8bio__seq24fasta__title__to__record(_M0L5titleS861, _M0L6_2atmpS2040);
  moonbit_decref(_M0L5titleS861);
  moonbit_decref(_M0L6_2atmpS2040);
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS835, _M0L6_2atmpS2039);
  moonbit_decref(_M0L6_2atmpS2039);
  joinlet_3074:;
  return _M0L7recordsS835;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0FP26biolab8bio__seq24fasta__title__to__record(
  moonbit_string_t _M0L5titleS830,
  moonbit_string_t _M0L8seq__strS826
) {
  moonbit_string_t _M0L2idS824;
  moonbit_string_t _M0L11descriptionS825;
  int32_t _M0L3posS829;
  moonbit_string_t _M0L7_2abindS832;
  int32_t _M0L6_2atmpS2017;
  struct _M0TPC16string10StringView _M0L6_2atmpS2016;
  int64_t _M0L7_2abindS831;
  int64_t _M0L6_2atmpS2015;
  struct _M0TPC16string10StringView _M0L6_2atmpS2014;
  moonbit_string_t _M0L6_2atmpS2013;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2009;
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS827;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2012;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2011;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2010;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_3077;
  #line 78 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7_2abindS832 = (moonbit_string_t)moonbit_string_literal_72.data;
  _M0L6_2atmpS2017 = Moonbit_array_length(_M0L7_2abindS832);
  _M0L6_2atmpS2016
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS832, .$1 = 0, .$2 = _M0L6_2atmpS2017
  };
  #line 79 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7_2abindS831
  = _M0MPC16string6String4find(_M0L5titleS830, _M0L6_2atmpS2016);
  moonbit_decref(_M0L6_2atmpS2016.$0);
  if (_M0L7_2abindS831 == 4294967296ll) {
    moonbit_incref(_M0L5titleS830);
    moonbit_incref(_M0L5titleS830);
    _M0L2idS824 = _M0L5titleS830;
    _M0L11descriptionS825 = _M0L5titleS830;
    goto join_823;
  } else {
    int64_t _M0L7_2aSomeS833 = _M0L7_2abindS831;
    int32_t _M0L6_2aposS834 = (int32_t)_M0L7_2aSomeS833;
    _M0L3posS829 = _M0L6_2aposS834;
    goto join_828;
  }
  join_828:;
  _M0L6_2atmpS2015 = (int64_t)_M0L3posS829;
  #line 80 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2014
  = _M0MPC16string6String11sub_2einner(_M0L5titleS830, 0, _M0L6_2atmpS2015);
  #line 80 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2013 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2014);
  moonbit_decref(_M0L6_2atmpS2014.$0);
  moonbit_incref(_M0L5titleS830);
  _M0L2idS824 = _M0L6_2atmpS2013;
  _M0L11descriptionS825 = _M0L5titleS830;
  goto join_823;
  join_823:;
  #line 84 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2009 = _M0MP26biolab8bio__seq3Seq3new(_M0L8seq__strS826);
  _M0L7_2abindS827 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2012 = _M0L7_2abindS827;
  _M0L6_2atmpS2011
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2012, .$1 = 0, .$2 = 0
  };
  #line 88 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2010 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2011, 0ll);
  moonbit_decref(_M0L6_2atmpS2011.$0);
  moonbit_incref(_M0L2idS824);
  _block_3077
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_3077)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
  _block_3077->$0 = _M0L6_2atmpS2009;
  _block_3077->$1 = _M0L2idS824;
  _block_3077->$2 = _M0L2idS824;
  _block_3077->$3 = _M0L11descriptionS825;
  _block_3077->$4 = _M0L6_2atmpS2010;
  return _block_3077;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord11new_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS818,
  moonbit_string_t _M0L2idS819,
  moonbit_string_t _M0L4nameS820,
  moonbit_string_t _M0L11descriptionS821
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS822;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2008;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2007;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2006;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_3078;
  #line 31 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS822 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2008 = _M0L7_2abindS822;
  _M0L6_2atmpS2007
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2008, .$1 = 0, .$2 = 0
  };
  #line 42 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS2006 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2007, 0ll);
  moonbit_decref(_M0L6_2atmpS2007.$0);
  moonbit_incref(_M0L3seqS818);
  moonbit_incref(_M0L2idS819);
  moonbit_incref(_M0L4nameS820);
  moonbit_incref(_M0L11descriptionS821);
  _block_3078
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_3078)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
  _block_3078->$0 = _M0L3seqS818;
  _block_3078->$1 = _M0L2idS819;
  _block_3078->$2 = _M0L4nameS820;
  _block_3078->$3 = _M0L11descriptionS821;
  _block_3078->$4 = _M0L6_2atmpS2006;
  return _block_3078;
}

int32_t _M0MP26biolab8bio__seq9SeqRecord6length(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS817
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2005;
  int32_t _result_3079;
  #line 84 "/home/developer/Documents/1/seq_record.mbt"
  _M0L3seqS2005 = _M0L4selfS817->$0;
  moonbit_incref(_M0L3seqS2005);
  #line 85 "/home/developer/Documents/1/seq_record.mbt"
  _result_3079 = _M0MP26biolab8bio__seq3Seq6length(_M0L3seqS2005);
  moonbit_decref(_M0L3seqS2005);
  return _result_3079;
}

int32_t _M0MP26biolab8bio__seq3Seq6length(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS816
) {
  moonbit_string_t _M0L4dataS2004;
  int32_t _result_3080;
  #line 86 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2004 = _M0L4selfS816->$0;
  moonbit_incref(_M0L4dataS2004);
  #line 87 "/home/developer/Documents/1/seq.mbt"
  _result_3080
  = _M0MPC16string6String20char__length_2einner(_M0L4dataS2004, 0, 4294967296ll);
  moonbit_decref(_M0L4dataS2004);
  return _result_3080;
}

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS815
) {
  moonbit_string_t _M0L8_2afieldS2742;
  #line 74 "/home/developer/Documents/1/seq.mbt"
  _M0L8_2afieldS2742 = _M0L4selfS815->$0;
  moonbit_incref(_M0L8_2afieldS2742);
  return _M0L8_2afieldS2742;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t _M0L1sS814
) {
  struct _M0TP26biolab8bio__seq3Seq* _block_3081;
  #line 62 "/home/developer/Documents/1/seq.mbt"
  moonbit_incref(_M0L1sS814);
  _block_3081
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3081)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 46, 0);
  _block_3081->$0 = _M0L1sS814;
  return _block_3081;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS803,
  int32_t _M0L5indexS804
) {
  int32_t _M0L3lenS802;
  int32_t _if__result_3082;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS802 = _M0L4selfS803->$1;
  if (_M0L5indexS804 >= 0) {
    _if__result_3082 = _M0L5indexS804 < _M0L3lenS802;
  } else {
    _if__result_3082 = 0;
  }
  if (_if__result_3082) {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2000;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2743;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2000
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS803);
    _M0L6_2atmpS2743
    = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L6_2atmpS2000[
        _M0L5indexS804
      ];
    if (_M0L6_2atmpS2743) {
      moonbit_incref(_M0L6_2atmpS2743);
    }
    moonbit_decref(_M0L6_2atmpS2000);
    return _M0L6_2atmpS2743;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

moonbit_string_t _M0MPC15array5Array2atGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS806,
  int32_t _M0L5indexS807
) {
  int32_t _M0L3lenS805;
  int32_t _if__result_3083;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS805 = _M0L4selfS806->$1;
  if (_M0L5indexS807 >= 0) {
    _if__result_3083 = _M0L5indexS807 < _M0L3lenS805;
  } else {
    _if__result_3083 = 0;
  }
  if (_if__result_3083) {
    moonbit_string_t* _M0L6_2atmpS2001;
    moonbit_string_t _M0L6_2atmpS2744;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2001 = _M0MPC15array5Array6bufferGsE(_M0L4selfS806);
    _M0L6_2atmpS2744 = (moonbit_string_t)_M0L6_2atmpS2001[_M0L5indexS807];
    moonbit_incref(_M0L6_2atmpS2744);
    moonbit_decref(_M0L6_2atmpS2001);
    return _M0L6_2atmpS2744;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TPC16string10StringView _M0MPC15array5Array2atGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS809,
  int32_t _M0L5indexS810
) {
  int32_t _M0L3lenS808;
  int32_t _if__result_3084;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS808 = _M0L4selfS809->$1;
  if (_M0L5indexS810 >= 0) {
    _if__result_3084 = _M0L5indexS810 < _M0L3lenS808;
  } else {
    _if__result_3084 = 0;
  }
  if (_if__result_3084) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS2002;
    struct _M0TPC16string10StringView _M0L6_2atmpS2745;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2002
    = _M0MPC15array5Array6bufferGRPC16string10StringViewE(_M0L4selfS809);
    _M0L6_2atmpS2745 = _M0L6_2atmpS2002[_M0L5indexS810];
    moonbit_incref(_M0L6_2atmpS2745.$0);
    moonbit_decref(_M0L6_2atmpS2002);
    return _M0L6_2atmpS2745;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array5Array2atGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS812,
  int32_t _M0L5indexS813
) {
  int32_t _M0L3lenS811;
  int32_t _if__result_3085;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS811 = _M0L4selfS812->$1;
  if (_M0L5indexS813 >= 0) {
    _if__result_3085 = _M0L5indexS813 < _M0L3lenS811;
  } else {
    _if__result_3085 = 0;
  }
  if (_if__result_3085) {
    int32_t* _M0L6_2atmpS2003;
    int32_t _result_3086;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2003 = _M0MPC15array5Array6bufferGiE(_M0L4selfS812);
    _result_3086 = (int32_t)_M0L6_2atmpS2003[_M0L5indexS813];
    moonbit_decref(_M0L6_2atmpS2003);
    return _result_3086;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB7printlnGsE(moonbit_string_t _M0L5inputS801) {
  moonbit_string_t _M0L6_2atmpS1999;
  #line 36 "/home/developer/.moon/lib/core/builtin/console.mbt"
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  _M0L6_2atmpS1999 = _M0IPC16string6StringPB4Show10to__string(_M0L5inputS801);
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  moonbit_println(_M0L6_2atmpS1999);
  moonbit_decref(_M0L6_2atmpS1999);
  return 0;
}

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t _M0L4selfS797) {
  int32_t _tmp_3087;
  uint32_t _M0L6_2atmpS1998;
  uint32_t _M0Lm3accS795;
  int32_t _M0L7_2abindS796;
  int32_t _M0L1iS798;
  uint32_t _M0L6_2atmpS1997;
  #line 522 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_3087 = 0;
  _M0L6_2atmpS1998 = *(uint32_t*)&_tmp_3087;
  _M0Lm3accS795 = _M0L6_2atmpS1998 + 374761393u;
  _M0L7_2abindS796 = Moonbit_array_length(_M0L4selfS797);
  _M0L1iS798 = 0;
  while (1) {
    if (_M0L1iS798 < _M0L7_2abindS796) {
      uint32_t _M0L6_2atmpS1992 = _M0Lm3accS795;
      int32_t _M0L6_2atmpS1995;
      int32_t _M0L6_2atmpS1994;
      uint32_t _M0L1vS799;
      uint32_t _M0L6_2atmpS1993;
      int32_t _M0L6_2atmpS1996;
      _M0Lm3accS795 = _M0L6_2atmpS1992 + 4u;
      _M0L6_2atmpS1995 = _M0L4selfS797[_M0L1iS798];
      _M0L6_2atmpS1994 = (int32_t)_M0L6_2atmpS1995;
      _M0L1vS799 = *(uint32_t*)&_M0L6_2atmpS1994;
      _M0L6_2atmpS1993 = _M0Lm3accS795;
      #line 527 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
      _M0Lm3accS795 = _M0FPB13consume4__acc(_M0L6_2atmpS1993, _M0L1vS799);
      _M0L6_2atmpS1996 = _M0L1iS798 + 1;
      _M0L1iS798 = _M0L6_2atmpS1996;
      continue;
    }
    break;
  }
  _M0L6_2atmpS1997 = _M0Lm3accS795;
  #line 529 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS1997);
}

struct _M0TUsRPB5ArrayGiEE* _M0MPB5Iter24nextGsRPB5ArrayGiEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0L4selfS794
) {
  #line 1099 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  #line 1100 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  return _M0MPB4Iter4nextGUsRPB5ArrayGiEEE(_M0L4selfS794);
}

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB3Map5iter2GsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS793
) {
  #line 725 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 727 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  return _M0MPB3Map4iterGsRPB5ArrayGiEE(_M0L4selfS793);
}

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB3Map4iterGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS783
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4headS1991;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE* _M0L11curr__entryS782;
  int32_t _M0L3lenS784;
  struct _M0TPB8MutLocalGiE* _M0L9remainingS785;
  struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__* _closure_3089;
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L6_2atmpS1982;
  int64_t _M0L6_2atmpS1983;
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _result_3090;
  #line 705 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4headS1991 = _M0L4selfS783->$5;
  if (_M0L4headS1991) {
    moonbit_incref(_M0L4headS1991);
  }
  _M0L11curr__entryS782
  = (struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE));
  Moonbit_object_header(_M0L11curr__entryS782)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 49, 0);
  _M0L11curr__entryS782->$0 = _M0L4headS1991;
  _M0L3lenS784 = _M0L4selfS783->$1;
  _M0L9remainingS785
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L9remainingS785)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L9remainingS785->$0 = _M0L3lenS784;
  _closure_3089
  = (struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__*)moonbit_malloc(sizeof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__));
  Moonbit_object_header(_closure_3089)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 52, 0);
  _closure_3089->code = &_M0MPB3Map4iterGsRPB5ArrayGiEEC1984l711;
  _closure_3089->$0 = _M0L9remainingS785;
  _closure_3089->$1 = _M0L11curr__entryS782;
  _M0L6_2atmpS1982 = (struct _M0TWEOUsRPB5ArrayGiEE*)_closure_3089;
  _M0L6_2atmpS1983 = (int64_t)_M0L3lenS784;
  #line 710 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _result_3090
  = _M0MPB4Iter3newGUsRPB5ArrayGiEEE(_M0L6_2atmpS1982, _M0L6_2atmpS1983);
  moonbit_decref(_M0L6_2atmpS1982);
  return _result_3090;
}

struct _M0TUsRPB5ArrayGiEE* _M0MPB3Map4iterGsRPB5ArrayGiEEC1984l711(
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L6_2aenvS1985
) {
  struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__* _M0L14_2acasted__envS1986;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE* _M0L11curr__entryS782;
  struct _M0TPB8MutLocalGiE* _M0L9remainingS785;
  int32_t _M0L3valS1987;
  #line 711 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14_2acasted__envS1986
  = (struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u1984__l711__*)_M0L6_2aenvS1985;
  _M0L11curr__entryS782 = _M0L14_2acasted__envS1986->$1;
  _M0L9remainingS785 = _M0L14_2acasted__envS1986->$0;
  _M0L3valS1987 = _M0L9remainingS785->$0;
  if (_M0L3valS1987 > 0) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS787 =
      _M0L11curr__entryS782->$0;
    if (_M0L7_2abindS787 == 0) {
      goto join_786;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS788 = _M0L7_2abindS787;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4_2axS789 = _M0L7_2aSomeS788;
      moonbit_string_t _M0L6_2akeyS790 = _M0L4_2axS789->$4;
      struct _M0TPB5ArrayGiE* _M0L8_2avalueS791 = _M0L4_2axS789->$5;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2anextS792 =
        _M0L4_2axS789->$1;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2747 =
        _M0L11curr__entryS782->$0;
      int32_t _M0L3valS1989;
      int32_t _M0L6_2atmpS1988;
      struct _M0TUsRPB5ArrayGiEE* _M0L8_2atupleS1990;
      if (_M0L7_2anextS792) {
        moonbit_incref(_M0L7_2anextS792);
      }
      moonbit_incref(_M0L8_2avalueS791);
      moonbit_incref(_M0L6_2akeyS790);
      if (_M0L6_2aoldS2747) {
        moonbit_decref(_M0L6_2aoldS2747);
      }
      _M0L11curr__entryS782->$0 = _M0L7_2anextS792;
      _M0L3valS1989 = _M0L9remainingS785->$0;
      _M0L6_2atmpS1988 = _M0L3valS1989 - 1;
      _M0L9remainingS785->$0 = _M0L6_2atmpS1988;
      _M0L8_2atupleS1990
      = (struct _M0TUsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TUsRPB5ArrayGiEE));
      Moonbit_object_header(_M0L8_2atupleS1990)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 56, 0);
      _M0L8_2atupleS1990->$0 = _M0L6_2akeyS790;
      _M0L8_2atupleS1990->$1 = _M0L8_2avalueS791;
      return _M0L8_2atupleS1990;
    }
  } else {
    goto join_786;
  }
  join_786:;
  return 0;
}

int32_t _M0MPB3Map6lengthGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS781
) {
  #line 641 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  return _M0L4selfS781->$1;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS768,
  moonbit_string_t _M0L3keyS764
) {
  int32_t _M0L4hashS763;
  int32_t _M0L14capacity__maskS1967;
  int32_t _M0L6_2atmpS1966;
  int32_t _M0L1iS765;
  int32_t _M0L3idxS766;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS763 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS764);
  _M0L14capacity__maskS1967 = _M0L4selfS768->$3;
  _M0L6_2atmpS1966 = _M0L4hashS763 & _M0L14capacity__maskS1967;
  _M0L1iS765 = 0;
  _M0L3idxS766 = _M0L6_2atmpS1966;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1965 =
      _M0L4selfS768->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS767;
    if (
      _M0L3idxS766 < 0
      || _M0L3idxS766 >= Moonbit_array_length(_M0L7entriesS1965)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS767
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1965[
        _M0L3idxS766
      ];
    if (_M0L7_2abindS767 == 0) {
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1954 = 0;
      return _M0L6_2atmpS1954;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS769 =
        _M0L7_2abindS767;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L8_2aentryS770 =
        _M0L7_2aSomeS769;
      int32_t _M0L4hashS1956 = _M0L8_2aentryS770->$3;
      int32_t _if__result_3093;
      int32_t _M0L3pslS1959;
      int32_t _M0L6_2atmpS1961;
      int32_t _M0L6_2atmpS1963;
      int32_t _M0L14capacity__maskS1964;
      int32_t _M0L6_2atmpS1962;
      if (_M0L4hashS1956 == _M0L4hashS763) {
        moonbit_string_t _M0L3keyS1955 = _M0L8_2aentryS770->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3093
        = _M0L3keyS1955 == _M0L3keyS764
          || Moonbit_array_length(_M0L3keyS1955)
             == Moonbit_array_length(_M0L3keyS764)
             && 0
                == memcmp(_M0L3keyS1955, _M0L3keyS764, Moonbit_array_length(_M0L3keyS1955) * 2);
      } else {
        _if__result_3093 = 0;
      }
      if (_if__result_3093) {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS1958 =
          _M0L8_2aentryS770->$5;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1957;
        moonbit_incref(_M0L5valueS1958);
        _M0L6_2atmpS1957 = _M0L5valueS1958;
        return _M0L6_2atmpS1957;
      } else {
        moonbit_incref(_M0L8_2aentryS770);
      }
      _M0L3pslS1959 = _M0L8_2aentryS770->$2;
      moonbit_decref(_M0L8_2aentryS770);
      if (_M0L1iS765 > _M0L3pslS1959) {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1960 = 0;
        return _M0L6_2atmpS1960;
      }
      _M0L6_2atmpS1961 = _M0L1iS765 + 1;
      _M0L6_2atmpS1963 = _M0L3idxS766 + 1;
      _M0L14capacity__maskS1964 = _M0L4selfS768->$3;
      _M0L6_2atmpS1962 = _M0L6_2atmpS1963 & _M0L14capacity__maskS1964;
      _M0L1iS765 = _M0L6_2atmpS1961;
      _M0L3idxS766 = _M0L6_2atmpS1962;
      continue;
    }
    break;
  }
}

struct _M0TPB5ArrayGiE* _M0MPB3Map3getGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS777,
  moonbit_string_t _M0L3keyS773
) {
  int32_t _M0L4hashS772;
  int32_t _M0L14capacity__maskS1981;
  int32_t _M0L6_2atmpS1980;
  int32_t _M0L1iS774;
  int32_t _M0L3idxS775;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS772 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS773);
  _M0L14capacity__maskS1981 = _M0L4selfS777->$3;
  _M0L6_2atmpS1980 = _M0L4hashS772 & _M0L14capacity__maskS1981;
  _M0L1iS774 = 0;
  _M0L3idxS775 = _M0L6_2atmpS1980;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1979 =
      _M0L4selfS777->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS776;
    if (
      _M0L3idxS775 < 0
      || _M0L3idxS775 >= Moonbit_array_length(_M0L7entriesS1979)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS776
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1979[_M0L3idxS775];
    if (_M0L7_2abindS776 == 0) {
      struct _M0TPB5ArrayGiE* _M0L6_2atmpS1968 = 0;
      return _M0L6_2atmpS1968;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS778 = _M0L7_2abindS776;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L8_2aentryS779 =
        _M0L7_2aSomeS778;
      int32_t _M0L4hashS1970 = _M0L8_2aentryS779->$3;
      int32_t _if__result_3095;
      int32_t _M0L3pslS1973;
      int32_t _M0L6_2atmpS1975;
      int32_t _M0L6_2atmpS1977;
      int32_t _M0L14capacity__maskS1978;
      int32_t _M0L6_2atmpS1976;
      if (_M0L4hashS1970 == _M0L4hashS772) {
        moonbit_string_t _M0L3keyS1969 = _M0L8_2aentryS779->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3095
        = _M0L3keyS1969 == _M0L3keyS773
          || Moonbit_array_length(_M0L3keyS1969)
             == Moonbit_array_length(_M0L3keyS773)
             && 0
                == memcmp(_M0L3keyS1969, _M0L3keyS773, Moonbit_array_length(_M0L3keyS1969) * 2);
      } else {
        _if__result_3095 = 0;
      }
      if (_if__result_3095) {
        struct _M0TPB5ArrayGiE* _M0L5valueS1972 = _M0L8_2aentryS779->$5;
        struct _M0TPB5ArrayGiE* _M0L6_2atmpS1971;
        moonbit_incref(_M0L5valueS1972);
        _M0L6_2atmpS1971 = _M0L5valueS1972;
        return _M0L6_2atmpS1971;
      } else {
        moonbit_incref(_M0L8_2aentryS779);
      }
      _M0L3pslS1973 = _M0L8_2aentryS779->$2;
      moonbit_decref(_M0L8_2aentryS779);
      if (_M0L1iS774 > _M0L3pslS1973) {
        struct _M0TPB5ArrayGiE* _M0L6_2atmpS1974 = 0;
        return _M0L6_2atmpS1974;
      }
      _M0L6_2atmpS1975 = _M0L1iS774 + 1;
      _M0L6_2atmpS1977 = _M0L3idxS775 + 1;
      _M0L14capacity__maskS1978 = _M0L4selfS777->$3;
      _M0L6_2atmpS1976 = _M0L6_2atmpS1977 & _M0L14capacity__maskS1978;
      _M0L1iS774 = _M0L6_2atmpS1975;
      _M0L3idxS775 = _M0L6_2atmpS1976;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0MPB3Map3MapGsRPB5ArrayGiEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L3arrS742,
  int64_t _M0L8capacityS744
) {
  int32_t _M0L3endS1941;
  int32_t _M0L5startS1942;
  int32_t _M0L6lengthS741;
  int32_t _M0L8capacityS743;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L1mS747;
  int32_t _M0L3endS1938;
  int32_t _M0L5startS1939;
  int32_t _M0L7_2abindS748;
  int32_t _M0L2__S749;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS1941 = _M0L3arrS742.$2;
  _M0L5startS1942 = _M0L3arrS742.$1;
  _M0L6lengthS741 = _M0L3endS1941 - _M0L5startS1942;
  if (_M0L8capacityS744 == 4294967296ll) {
    if (_M0L6lengthS741 == 0) {
      _M0L8capacityS743 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS743 = _M0FPB21capacity__for__length(_M0L6lengthS741);
    }
  } else {
    int64_t _M0L7_2aSomeS745 = _M0L8capacityS744;
    int32_t _M0L11_2acapacityS746 = (int32_t)_M0L7_2aSomeS745;
    int32_t _M0L6_2atmpS1940;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS1940 = _M0FPB21capacity__for__length(_M0L6lengthS741);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS743
    = _M0MPC13int3Int3max(_M0L11_2acapacityS746, _M0L6_2atmpS1940);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS747 = _M0FPB8new__mapGsRPB5ArrayGiEE(_M0L8capacityS743);
  _M0L3endS1938 = _M0L3arrS742.$2;
  _M0L5startS1939 = _M0L3arrS742.$1;
  _M0L7_2abindS748 = _M0L3endS1938 - _M0L5startS1939;
  _M0L2__S749 = 0;
  while (1) {
    if (_M0L2__S749 < _M0L7_2abindS748) {
      struct _M0TUsRPB5ArrayGiEE** _M0L3bufS1935 = _M0L3arrS742.$0;
      int32_t _M0L5startS1937 = _M0L3arrS742.$1;
      int32_t _M0L6_2atmpS1936 = _M0L5startS1937 + _M0L2__S749;
      struct _M0TUsRPB5ArrayGiEE* _M0L1eS750 =
        (struct _M0TUsRPB5ArrayGiEE*)_M0L3bufS1935[_M0L6_2atmpS1936];
      moonbit_string_t _M0L6_2atmpS1932 = _M0L1eS750->$0;
      struct _M0TPB5ArrayGiE* _M0L6_2atmpS1933 = _M0L1eS750->$1;
      int32_t _M0L6_2atmpS1934;
      moonbit_incref(_M0L6_2atmpS1933);
      moonbit_incref(_M0L6_2atmpS1932);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L1mS747, _M0L6_2atmpS1932, _M0L6_2atmpS1933);
      moonbit_decref(_M0L6_2atmpS1932);
      moonbit_decref(_M0L6_2atmpS1933);
      _M0L6_2atmpS1934 = _M0L2__S749 + 1;
      _M0L2__S749 = _M0L6_2atmpS1934;
      continue;
    }
    break;
  }
  return _M0L1mS747;
}

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0MPB3Map3MapGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE _M0L3arrS753,
  int64_t _M0L8capacityS755
) {
  int32_t _M0L3endS1952;
  int32_t _M0L5startS1953;
  int32_t _M0L6lengthS752;
  int32_t _M0L8capacityS754;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L1mS758;
  int32_t _M0L3endS1949;
  int32_t _M0L5startS1950;
  int32_t _M0L7_2abindS759;
  int32_t _M0L2__S760;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS1952 = _M0L3arrS753.$2;
  _M0L5startS1953 = _M0L3arrS753.$1;
  _M0L6lengthS752 = _M0L3endS1952 - _M0L5startS1953;
  if (_M0L8capacityS755 == 4294967296ll) {
    if (_M0L6lengthS752 == 0) {
      _M0L8capacityS754 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS754 = _M0FPB21capacity__for__length(_M0L6lengthS752);
    }
  } else {
    int64_t _M0L7_2aSomeS756 = _M0L8capacityS755;
    int32_t _M0L11_2acapacityS757 = (int32_t)_M0L7_2aSomeS756;
    int32_t _M0L6_2atmpS1951;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS1951 = _M0FPB21capacity__for__length(_M0L6lengthS752);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS754
    = _M0MPC13int3Int3max(_M0L11_2acapacityS757, _M0L6_2atmpS1951);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS758
  = _M0FPB8new__mapGsRP26biolab8bio__seq9SeqRecordE(_M0L8capacityS754);
  _M0L3endS1949 = _M0L3arrS753.$2;
  _M0L5startS1950 = _M0L3arrS753.$1;
  _M0L7_2abindS759 = _M0L3endS1949 - _M0L5startS1950;
  _M0L2__S760 = 0;
  while (1) {
    if (_M0L2__S760 < _M0L7_2abindS759) {
      struct _M0TUsRP26biolab8bio__seq9SeqRecordE** _M0L3bufS1946 =
        _M0L3arrS753.$0;
      int32_t _M0L5startS1948 = _M0L3arrS753.$1;
      int32_t _M0L6_2atmpS1947 = _M0L5startS1948 + _M0L2__S760;
      struct _M0TUsRP26biolab8bio__seq9SeqRecordE* _M0L1eS761 =
        (struct _M0TUsRP26biolab8bio__seq9SeqRecordE*)_M0L3bufS1946[
          _M0L6_2atmpS1947
        ];
      moonbit_string_t _M0L6_2atmpS1943 = _M0L1eS761->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1944 =
        _M0L1eS761->$1;
      int32_t _M0L6_2atmpS1945;
      moonbit_incref(_M0L6_2atmpS1944);
      moonbit_incref(_M0L6_2atmpS1943);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(_M0L1mS758, _M0L6_2atmpS1943, _M0L6_2atmpS1944);
      moonbit_decref(_M0L6_2atmpS1943);
      moonbit_decref(_M0L6_2atmpS1944);
      _M0L6_2atmpS1945 = _M0L2__S760 + 1;
      _M0L2__S760 = _M0L6_2atmpS1945;
      continue;
    }
    break;
  }
  return _M0L1mS758;
}

int32_t _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS735,
  moonbit_string_t _M0L3keyS736,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS737
) {
  int32_t _M0L6_2atmpS1930;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1930 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS736);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS735, _M0L3keyS736, _M0L5valueS737, _M0L6_2atmpS1930);
  return 0;
}

int32_t _M0MPB3Map3setGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS738,
  moonbit_string_t _M0L3keyS739,
  struct _M0TPB5ArrayGiE* _M0L5valueS740
) {
  int32_t _M0L6_2atmpS1931;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1931 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS739);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(_M0L4selfS738, _M0L3keyS739, _M0L5valueS740, _M0L6_2atmpS1931);
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS706,
  moonbit_string_t _M0L3keyS712,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS713,
  int32_t _M0L4hashS708
) {
  int32_t _M0L14capacity__maskS1911;
  int32_t _M0L6_2atmpS1910;
  int32_t _M0L3pslS703;
  int32_t _M0L3idxS704;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS1911 = _M0L4selfS706->$3;
  _M0L6_2atmpS1910 = _M0L4hashS708 & _M0L14capacity__maskS1911;
  _M0L3pslS703 = 0;
  _M0L3idxS704 = _M0L6_2atmpS1910;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1909 =
      _M0L4selfS706->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS705;
    if (
      _M0L3idxS704 < 0
      || _M0L3idxS704 >= Moonbit_array_length(_M0L7entriesS1909)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS705
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1909[
        _M0L3idxS704
      ];
    if (_M0L7_2abindS705 == 0) {
      int32_t _M0L4sizeS1894 = _M0L4selfS706->$1;
      int32_t _M0L8grow__atS1895 = _M0L4selfS706->$4;
      int32_t _M0L7_2abindS709;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS710;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS711;
      if (_M0L4sizeS1894 >= _M0L8grow__atS1895) {
        int32_t _M0L14capacity__maskS1897;
        int32_t _M0L6_2atmpS1896;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS706);
        _M0L14capacity__maskS1897 = _M0L4selfS706->$3;
        _M0L6_2atmpS1896 = _M0L4hashS708 & _M0L14capacity__maskS1897;
        _M0L3pslS703 = 0;
        _M0L3idxS704 = _M0L6_2atmpS1896;
        continue;
      }
      _M0L7_2abindS709 = _M0L4selfS706->$6;
      _M0L7_2abindS710 = 0;
      moonbit_incref(_M0L3keyS712);
      moonbit_incref(_M0L5valueS713);
      _M0L5entryS711
      = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE));
      Moonbit_object_header(_M0L5entryS711)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 60, 0);
      _M0L5entryS711->$0 = _M0L7_2abindS709;
      _M0L5entryS711->$1 = _M0L7_2abindS710;
      _M0L5entryS711->$2 = _M0L3pslS703;
      _M0L5entryS711->$3 = _M0L4hashS708;
      _M0L5entryS711->$4 = _M0L3keyS712;
      _M0L5entryS711->$5 = _M0L5valueS713;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS706, _M0L3idxS704, _M0L5entryS711);
      moonbit_decref(_M0L5entryS711);
      return 0;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS714 =
        _M0L7_2abindS705;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L14_2acurr__entryS715 =
        _M0L7_2aSomeS714;
      int32_t _M0L4hashS1899 = _M0L14_2acurr__entryS715->$3;
      int32_t _if__result_3099;
      int32_t _M0L3pslS1900;
      int32_t _M0L6_2atmpS1905;
      int32_t _M0L6_2atmpS1907;
      int32_t _M0L14capacity__maskS1908;
      int32_t _M0L6_2atmpS1906;
      if (_M0L4hashS1899 == _M0L4hashS708) {
        moonbit_string_t _M0L3keyS1898 = _M0L14_2acurr__entryS715->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3099
        = _M0L3keyS1898 == _M0L3keyS712
          || Moonbit_array_length(_M0L3keyS1898)
             == Moonbit_array_length(_M0L3keyS712)
             && 0
                == memcmp(_M0L3keyS1898, _M0L3keyS712, Moonbit_array_length(_M0L3keyS1898) * 2);
      } else {
        _if__result_3099 = 0;
      }
      if (_if__result_3099) {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS2768 =
          _M0L14_2acurr__entryS715->$5;
        moonbit_incref(_M0L5valueS713);
        moonbit_decref(_M0L6_2aoldS2768);
        _M0L14_2acurr__entryS715->$5 = _M0L5valueS713;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS715);
      }
      _M0L3pslS1900 = _M0L14_2acurr__entryS715->$2;
      if (_M0L3pslS703 > _M0L3pslS1900) {
        int32_t _M0L4sizeS1901 = _M0L4selfS706->$1;
        int32_t _M0L8grow__atS1902 = _M0L4selfS706->$4;
        int32_t _M0L7_2abindS716;
        struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS717;
        struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS718;
        if (_M0L4sizeS1901 >= _M0L8grow__atS1902) {
          int32_t _M0L14capacity__maskS1904;
          int32_t _M0L6_2atmpS1903;
          moonbit_decref(_M0L14_2acurr__entryS715);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS706);
          _M0L14capacity__maskS1904 = _M0L4selfS706->$3;
          _M0L6_2atmpS1903 = _M0L4hashS708 & _M0L14capacity__maskS1904;
          _M0L3pslS703 = 0;
          _M0L3idxS704 = _M0L6_2atmpS1903;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS706, _M0L3idxS704, _M0L14_2acurr__entryS715);
        moonbit_decref(_M0L14_2acurr__entryS715);
        _M0L7_2abindS716 = _M0L4selfS706->$6;
        _M0L7_2abindS717 = 0;
        moonbit_incref(_M0L3keyS712);
        moonbit_incref(_M0L5valueS713);
        _M0L5entryS718
        = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE));
        Moonbit_object_header(_M0L5entryS718)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 60, 0);
        _M0L5entryS718->$0 = _M0L7_2abindS716;
        _M0L5entryS718->$1 = _M0L7_2abindS717;
        _M0L5entryS718->$2 = _M0L3pslS703;
        _M0L5entryS718->$3 = _M0L4hashS708;
        _M0L5entryS718->$4 = _M0L3keyS712;
        _M0L5entryS718->$5 = _M0L5valueS713;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS706, _M0L3idxS704, _M0L5entryS718);
        moonbit_decref(_M0L5entryS718);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS715);
      }
      _M0L6_2atmpS1905 = _M0L3pslS703 + 1;
      _M0L6_2atmpS1907 = _M0L3idxS704 + 1;
      _M0L14capacity__maskS1908 = _M0L4selfS706->$3;
      _M0L6_2atmpS1906 = _M0L6_2atmpS1907 & _M0L14capacity__maskS1908;
      _M0L3pslS703 = _M0L6_2atmpS1905;
      _M0L3idxS704 = _M0L6_2atmpS1906;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS722,
  moonbit_string_t _M0L3keyS728,
  struct _M0TPB5ArrayGiE* _M0L5valueS729,
  int32_t _M0L4hashS724
) {
  int32_t _M0L14capacity__maskS1929;
  int32_t _M0L6_2atmpS1928;
  int32_t _M0L3pslS719;
  int32_t _M0L3idxS720;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS1929 = _M0L4selfS722->$3;
  _M0L6_2atmpS1928 = _M0L4hashS724 & _M0L14capacity__maskS1929;
  _M0L3pslS719 = 0;
  _M0L3idxS720 = _M0L6_2atmpS1928;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1927 =
      _M0L4selfS722->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS721;
    if (
      _M0L3idxS720 < 0
      || _M0L3idxS720 >= Moonbit_array_length(_M0L7entriesS1927)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS721
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1927[_M0L3idxS720];
    if (_M0L7_2abindS721 == 0) {
      int32_t _M0L4sizeS1912 = _M0L4selfS722->$1;
      int32_t _M0L8grow__atS1913 = _M0L4selfS722->$4;
      int32_t _M0L7_2abindS725;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS726;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS727;
      if (_M0L4sizeS1912 >= _M0L8grow__atS1913) {
        int32_t _M0L14capacity__maskS1915;
        int32_t _M0L6_2atmpS1914;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB5ArrayGiEE(_M0L4selfS722);
        _M0L14capacity__maskS1915 = _M0L4selfS722->$3;
        _M0L6_2atmpS1914 = _M0L4hashS724 & _M0L14capacity__maskS1915;
        _M0L3pslS719 = 0;
        _M0L3idxS720 = _M0L6_2atmpS1914;
        continue;
      }
      _M0L7_2abindS725 = _M0L4selfS722->$6;
      _M0L7_2abindS726 = 0;
      moonbit_incref(_M0L3keyS728);
      moonbit_incref(_M0L5valueS729);
      _M0L5entryS727
      = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE));
      Moonbit_object_header(_M0L5entryS727)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 65, 0);
      _M0L5entryS727->$0 = _M0L7_2abindS725;
      _M0L5entryS727->$1 = _M0L7_2abindS726;
      _M0L5entryS727->$2 = _M0L3pslS719;
      _M0L5entryS727->$3 = _M0L4hashS724;
      _M0L5entryS727->$4 = _M0L3keyS728;
      _M0L5entryS727->$5 = _M0L5valueS729;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS722, _M0L3idxS720, _M0L5entryS727);
      moonbit_decref(_M0L5entryS727);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS730 = _M0L7_2abindS721;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L14_2acurr__entryS731 =
        _M0L7_2aSomeS730;
      int32_t _M0L4hashS1917 = _M0L14_2acurr__entryS731->$3;
      int32_t _if__result_3101;
      int32_t _M0L3pslS1918;
      int32_t _M0L6_2atmpS1923;
      int32_t _M0L6_2atmpS1925;
      int32_t _M0L14capacity__maskS1926;
      int32_t _M0L6_2atmpS1924;
      if (_M0L4hashS1917 == _M0L4hashS724) {
        moonbit_string_t _M0L3keyS1916 = _M0L14_2acurr__entryS731->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3101
        = _M0L3keyS1916 == _M0L3keyS728
          || Moonbit_array_length(_M0L3keyS1916)
             == Moonbit_array_length(_M0L3keyS728)
             && 0
                == memcmp(_M0L3keyS1916, _M0L3keyS728, Moonbit_array_length(_M0L3keyS1916) * 2);
      } else {
        _if__result_3101 = 0;
      }
      if (_if__result_3101) {
        struct _M0TPB5ArrayGiE* _M0L6_2aoldS2772 =
          _M0L14_2acurr__entryS731->$5;
        moonbit_incref(_M0L5valueS729);
        moonbit_decref(_M0L6_2aoldS2772);
        _M0L14_2acurr__entryS731->$5 = _M0L5valueS729;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS731);
      }
      _M0L3pslS1918 = _M0L14_2acurr__entryS731->$2;
      if (_M0L3pslS719 > _M0L3pslS1918) {
        int32_t _M0L4sizeS1919 = _M0L4selfS722->$1;
        int32_t _M0L8grow__atS1920 = _M0L4selfS722->$4;
        int32_t _M0L7_2abindS732;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS733;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS734;
        if (_M0L4sizeS1919 >= _M0L8grow__atS1920) {
          int32_t _M0L14capacity__maskS1922;
          int32_t _M0L6_2atmpS1921;
          moonbit_decref(_M0L14_2acurr__entryS731);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB5ArrayGiEE(_M0L4selfS722);
          _M0L14capacity__maskS1922 = _M0L4selfS722->$3;
          _M0L6_2atmpS1921 = _M0L4hashS724 & _M0L14capacity__maskS1922;
          _M0L3pslS719 = 0;
          _M0L3idxS720 = _M0L6_2atmpS1921;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGiEE(_M0L4selfS722, _M0L3idxS720, _M0L14_2acurr__entryS731);
        moonbit_decref(_M0L14_2acurr__entryS731);
        _M0L7_2abindS732 = _M0L4selfS722->$6;
        _M0L7_2abindS733 = 0;
        moonbit_incref(_M0L3keyS728);
        moonbit_incref(_M0L5valueS729);
        _M0L5entryS734
        = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE));
        Moonbit_object_header(_M0L5entryS734)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 65, 0);
        _M0L5entryS734->$0 = _M0L7_2abindS732;
        _M0L5entryS734->$1 = _M0L7_2abindS733;
        _M0L5entryS734->$2 = _M0L3pslS719;
        _M0L5entryS734->$3 = _M0L4hashS724;
        _M0L5entryS734->$4 = _M0L3keyS728;
        _M0L5entryS734->$5 = _M0L5valueS729;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS722, _M0L3idxS720, _M0L5entryS734);
        moonbit_decref(_M0L5entryS734);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS731);
      }
      _M0L6_2atmpS1923 = _M0L3pslS719 + 1;
      _M0L6_2atmpS1925 = _M0L3idxS720 + 1;
      _M0L14capacity__maskS1926 = _M0L4selfS722->$3;
      _M0L6_2atmpS1924 = _M0L6_2atmpS1925 & _M0L14capacity__maskS1926;
      _M0L3pslS719 = _M0L6_2atmpS1923;
      _M0L3idxS720 = _M0L6_2atmpS1924;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS688
) {
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L9old__headS687;
  int32_t _M0L8capacityS1885;
  int32_t _M0L13new__capacityS689;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1879;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L6_2atmpS1878;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L6_2aoldS2779;
  int32_t _M0L6_2atmpS1880;
  int32_t _M0L8capacityS1882;
  int32_t _M0L6_2atmpS1881;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1883;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS2778;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L1xS690;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS687 = _M0L4selfS688->$5;
  _M0L8capacityS1885 = _M0L4selfS688->$2;
  _M0L13new__capacityS689 = _M0L8capacityS1885 << 1;
  _M0L6_2atmpS1879 = 0;
  _M0L6_2atmpS1878
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE**)moonbit_make_ref_array(_M0L13new__capacityS689, _M0L6_2atmpS1879);
  _M0L6_2aoldS2779 = _M0L4selfS688->$0;
  if (_M0L9old__headS687) {
    moonbit_incref(_M0L9old__headS687);
  }
  moonbit_decref(_M0L6_2aoldS2779);
  _M0L4selfS688->$0 = _M0L6_2atmpS1878;
  _M0L4selfS688->$2 = _M0L13new__capacityS689;
  _M0L6_2atmpS1880 = _M0L13new__capacityS689 - 1;
  _M0L4selfS688->$3 = _M0L6_2atmpS1880;
  _M0L8capacityS1882 = _M0L4selfS688->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1881 = _M0FPB21calc__grow__threshold(_M0L8capacityS1882);
  _M0L4selfS688->$4 = _M0L6_2atmpS1881;
  _M0L4selfS688->$1 = 0;
  _M0L6_2atmpS1883 = 0;
  _M0L6_2aoldS2778 = _M0L4selfS688->$5;
  if (_M0L6_2aoldS2778) {
    moonbit_decref(_M0L6_2aoldS2778);
  }
  _M0L4selfS688->$5 = _M0L6_2atmpS1883;
  _M0L4selfS688->$6 = -1;
  _M0L1xS690 = _M0L9old__headS687;
  while (1) {
    if (_M0L1xS690 == 0) {
      if (_M0L1xS690) {
        moonbit_decref(_M0L1xS690);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS692 =
        _M0L1xS690;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L4_2aeS693 =
        _M0L7_2aSomeS692;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L15next__in__chainS694 =
        _M0L4_2aeS693->$1;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1884 =
        0;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS2776 =
        _M0L4_2aeS693->$1;
      if (_M0L15next__in__chainS694) {
        moonbit_incref(_M0L15next__in__chainS694);
      }
      if (_M0L6_2aoldS2776) {
        moonbit_decref(_M0L6_2aoldS2776);
      }
      _M0L4_2aeS693->$1 = _M0L6_2atmpS1884;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS688, _M0L4_2aeS693);
      moonbit_decref(_M0L4_2aeS693);
      _M0L1xS690 = _M0L15next__in__chainS694;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS696
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L9old__headS695;
  int32_t _M0L8capacityS1893;
  int32_t _M0L13new__capacityS697;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1887;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L6_2atmpS1886;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L6_2aoldS2784;
  int32_t _M0L6_2atmpS1888;
  int32_t _M0L8capacityS1890;
  int32_t _M0L6_2atmpS1889;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1891;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2783;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L1xS698;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS695 = _M0L4selfS696->$5;
  _M0L8capacityS1893 = _M0L4selfS696->$2;
  _M0L13new__capacityS697 = _M0L8capacityS1893 << 1;
  _M0L6_2atmpS1887 = 0;
  _M0L6_2atmpS1886
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE**)moonbit_make_ref_array(_M0L13new__capacityS697, _M0L6_2atmpS1887);
  _M0L6_2aoldS2784 = _M0L4selfS696->$0;
  if (_M0L9old__headS695) {
    moonbit_incref(_M0L9old__headS695);
  }
  moonbit_decref(_M0L6_2aoldS2784);
  _M0L4selfS696->$0 = _M0L6_2atmpS1886;
  _M0L4selfS696->$2 = _M0L13new__capacityS697;
  _M0L6_2atmpS1888 = _M0L13new__capacityS697 - 1;
  _M0L4selfS696->$3 = _M0L6_2atmpS1888;
  _M0L8capacityS1890 = _M0L4selfS696->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1889 = _M0FPB21calc__grow__threshold(_M0L8capacityS1890);
  _M0L4selfS696->$4 = _M0L6_2atmpS1889;
  _M0L4selfS696->$1 = 0;
  _M0L6_2atmpS1891 = 0;
  _M0L6_2aoldS2783 = _M0L4selfS696->$5;
  if (_M0L6_2aoldS2783) {
    moonbit_decref(_M0L6_2aoldS2783);
  }
  _M0L4selfS696->$5 = _M0L6_2atmpS1891;
  _M0L4selfS696->$6 = -1;
  _M0L1xS698 = _M0L9old__headS695;
  while (1) {
    if (_M0L1xS698 == 0) {
      if (_M0L1xS698) {
        moonbit_decref(_M0L1xS698);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS700 = _M0L1xS698;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4_2aeS701 = _M0L7_2aSomeS700;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L15next__in__chainS702 =
        _M0L4_2aeS701->$1;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1892 = 0;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2781 =
        _M0L4_2aeS701->$1;
      if (_M0L15next__in__chainS702) {
        moonbit_incref(_M0L15next__in__chainS702);
      }
      if (_M0L6_2aoldS2781) {
        moonbit_decref(_M0L6_2aoldS2781);
      }
      _M0L4_2aeS701->$1 = _M0L6_2atmpS1892;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(_M0L4selfS696, _M0L4_2aeS701);
      moonbit_decref(_M0L4_2aeS701);
      _M0L1xS698 = _M0L15next__in__chainS702;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS674,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5outerS670
) {
  int32_t _M0L4hashS669;
  int32_t _M0L14capacity__maskS1867;
  int32_t _M0L6_2atmpS1866;
  int32_t _M0L3pslS671;
  int32_t _M0L3idxS672;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS669 = _M0L5outerS670->$3;
  _M0L14capacity__maskS1867 = _M0L4selfS674->$3;
  _M0L6_2atmpS1866 = _M0L4hashS669 & _M0L14capacity__maskS1867;
  _M0L3pslS671 = 0;
  _M0L3idxS672 = _M0L6_2atmpS1866;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1865 =
      _M0L4selfS674->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS673;
    if (
      _M0L3idxS672 < 0
      || _M0L3idxS672 >= Moonbit_array_length(_M0L7entriesS1865)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS673
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1865[
        _M0L3idxS672
      ];
    if (_M0L7_2abindS673 == 0) {
      int32_t _M0L4tailS1858;
      _M0L5outerS670->$2 = _M0L3pslS671;
      _M0L4tailS1858 = _M0L4selfS674->$6;
      _M0L5outerS670->$0 = _M0L4tailS1858;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS674, _M0L3idxS672, _M0L5outerS670);
      return 0;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS675 =
        _M0L7_2abindS673;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2acurrS676 =
        _M0L7_2aSomeS675;
      int32_t _M0L3pslS1859 = _M0L7_2acurrS676->$2;
      if (_M0L3pslS671 > _M0L3pslS1859) {
        int32_t _M0L4tailS1860;
        moonbit_incref(_M0L7_2acurrS676);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS674, _M0L3idxS672, _M0L7_2acurrS676);
        moonbit_decref(_M0L7_2acurrS676);
        _M0L5outerS670->$2 = _M0L3pslS671;
        _M0L4tailS1860 = _M0L4selfS674->$6;
        _M0L5outerS670->$0 = _M0L4tailS1860;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS674, _M0L3idxS672, _M0L5outerS670);
        return 0;
      } else {
        int32_t _M0L6_2atmpS1861 = _M0L3pslS671 + 1;
        int32_t _M0L6_2atmpS1863 = _M0L3idxS672 + 1;
        int32_t _M0L14capacity__maskS1864 = _M0L4selfS674->$3;
        int32_t _M0L6_2atmpS1862 =
          _M0L6_2atmpS1863 & _M0L14capacity__maskS1864;
        _M0L3pslS671 = _M0L6_2atmpS1861;
        _M0L3idxS672 = _M0L6_2atmpS1862;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS683,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5outerS679
) {
  int32_t _M0L4hashS678;
  int32_t _M0L14capacity__maskS1877;
  int32_t _M0L6_2atmpS1876;
  int32_t _M0L3pslS680;
  int32_t _M0L3idxS681;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS678 = _M0L5outerS679->$3;
  _M0L14capacity__maskS1877 = _M0L4selfS683->$3;
  _M0L6_2atmpS1876 = _M0L4hashS678 & _M0L14capacity__maskS1877;
  _M0L3pslS680 = 0;
  _M0L3idxS681 = _M0L6_2atmpS1876;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1875 =
      _M0L4selfS683->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS682;
    if (
      _M0L3idxS681 < 0
      || _M0L3idxS681 >= Moonbit_array_length(_M0L7entriesS1875)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS682
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1875[_M0L3idxS681];
    if (_M0L7_2abindS682 == 0) {
      int32_t _M0L4tailS1868;
      _M0L5outerS679->$2 = _M0L3pslS680;
      _M0L4tailS1868 = _M0L4selfS683->$6;
      _M0L5outerS679->$0 = _M0L4tailS1868;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS683, _M0L3idxS681, _M0L5outerS679);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS684 = _M0L7_2abindS682;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2acurrS685 = _M0L7_2aSomeS684;
      int32_t _M0L3pslS1869 = _M0L7_2acurrS685->$2;
      if (_M0L3pslS680 > _M0L3pslS1869) {
        int32_t _M0L4tailS1870;
        moonbit_incref(_M0L7_2acurrS685);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGiEE(_M0L4selfS683, _M0L3idxS681, _M0L7_2acurrS685);
        moonbit_decref(_M0L7_2acurrS685);
        _M0L5outerS679->$2 = _M0L3pslS680;
        _M0L4tailS1870 = _M0L4selfS683->$6;
        _M0L5outerS679->$0 = _M0L4tailS1870;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS683, _M0L3idxS681, _M0L5outerS679);
        return 0;
      } else {
        int32_t _M0L6_2atmpS1871 = _M0L3pslS680 + 1;
        int32_t _M0L6_2atmpS1873 = _M0L3idxS681 + 1;
        int32_t _M0L14capacity__maskS1874 = _M0L4selfS683->$3;
        int32_t _M0L6_2atmpS1872 =
          _M0L6_2atmpS1873 & _M0L14capacity__maskS1874;
        _M0L3pslS680 = _M0L6_2atmpS1871;
        _M0L3idxS681 = _M0L6_2atmpS1872;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS653,
  int32_t _M0L3idxS658,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS657
) {
  int32_t _M0L3pslS1841;
  int32_t _M0L6_2atmpS1837;
  int32_t _M0L6_2atmpS1839;
  int32_t _M0L14capacity__maskS1840;
  int32_t _M0L6_2atmpS1838;
  int32_t _M0L3pslS649;
  int32_t _M0L3idxS650;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS651;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS1841 = _M0L5entryS657->$2;
  _M0L6_2atmpS1837 = _M0L3pslS1841 + 1;
  _M0L6_2atmpS1839 = _M0L3idxS658 + 1;
  _M0L14capacity__maskS1840 = _M0L4selfS653->$3;
  _M0L6_2atmpS1838 = _M0L6_2atmpS1839 & _M0L14capacity__maskS1840;
  moonbit_incref(_M0L5entryS657);
  _M0L3pslS649 = _M0L6_2atmpS1837;
  _M0L3idxS650 = _M0L6_2atmpS1838;
  _M0L5entryS651 = _M0L5entryS657;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1836 =
      _M0L4selfS653->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS652;
    if (
      _M0L3idxS650 < 0
      || _M0L3idxS650 >= Moonbit_array_length(_M0L7entriesS1836)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS652
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1836[
        _M0L3idxS650
      ];
    if (_M0L7_2abindS652 == 0) {
      _M0L5entryS651->$2 = _M0L3pslS649;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS653, _M0L5entryS651, _M0L3idxS650);
      moonbit_decref(_M0L5entryS651);
      break;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS655 =
        _M0L7_2abindS652;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L14_2acurr__entryS656 =
        _M0L7_2aSomeS655;
      int32_t _M0L3pslS1826 = _M0L14_2acurr__entryS656->$2;
      if (_M0L3pslS649 > _M0L3pslS1826) {
        int32_t _M0L3pslS1831;
        int32_t _M0L6_2atmpS1827;
        int32_t _M0L6_2atmpS1829;
        int32_t _M0L14capacity__maskS1830;
        int32_t _M0L6_2atmpS1828;
        _M0L5entryS651->$2 = _M0L3pslS649;
        moonbit_incref(_M0L14_2acurr__entryS656);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS653, _M0L5entryS651, _M0L3idxS650);
        moonbit_decref(_M0L5entryS651);
        _M0L3pslS1831 = _M0L14_2acurr__entryS656->$2;
        _M0L6_2atmpS1827 = _M0L3pslS1831 + 1;
        _M0L6_2atmpS1829 = _M0L3idxS650 + 1;
        _M0L14capacity__maskS1830 = _M0L4selfS653->$3;
        _M0L6_2atmpS1828 = _M0L6_2atmpS1829 & _M0L14capacity__maskS1830;
        _M0L3pslS649 = _M0L6_2atmpS1827;
        _M0L3idxS650 = _M0L6_2atmpS1828;
        _M0L5entryS651 = _M0L14_2acurr__entryS656;
        continue;
      } else {
        int32_t _M0L6_2atmpS1832 = _M0L3pslS649 + 1;
        int32_t _M0L6_2atmpS1834 = _M0L3idxS650 + 1;
        int32_t _M0L14capacity__maskS1835 = _M0L4selfS653->$3;
        int32_t _M0L6_2atmpS1833 =
          _M0L6_2atmpS1834 & _M0L14capacity__maskS1835;
        struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _tmp_3107 =
          _M0L5entryS651;
        _M0L3pslS649 = _M0L6_2atmpS1832;
        _M0L3idxS650 = _M0L6_2atmpS1833;
        _M0L5entryS651 = _tmp_3107;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS663,
  int32_t _M0L3idxS668,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS667
) {
  int32_t _M0L3pslS1857;
  int32_t _M0L6_2atmpS1853;
  int32_t _M0L6_2atmpS1855;
  int32_t _M0L14capacity__maskS1856;
  int32_t _M0L6_2atmpS1854;
  int32_t _M0L3pslS659;
  int32_t _M0L3idxS660;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS661;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS1857 = _M0L5entryS667->$2;
  _M0L6_2atmpS1853 = _M0L3pslS1857 + 1;
  _M0L6_2atmpS1855 = _M0L3idxS668 + 1;
  _M0L14capacity__maskS1856 = _M0L4selfS663->$3;
  _M0L6_2atmpS1854 = _M0L6_2atmpS1855 & _M0L14capacity__maskS1856;
  moonbit_incref(_M0L5entryS667);
  _M0L3pslS659 = _M0L6_2atmpS1853;
  _M0L3idxS660 = _M0L6_2atmpS1854;
  _M0L5entryS661 = _M0L5entryS667;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1852 =
      _M0L4selfS663->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS662;
    if (
      _M0L3idxS660 < 0
      || _M0L3idxS660 >= Moonbit_array_length(_M0L7entriesS1852)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS662
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1852[_M0L3idxS660];
    if (_M0L7_2abindS662 == 0) {
      _M0L5entryS661->$2 = _M0L3pslS659;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB5ArrayGiEE(_M0L4selfS663, _M0L5entryS661, _M0L3idxS660);
      moonbit_decref(_M0L5entryS661);
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS665 = _M0L7_2abindS662;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L14_2acurr__entryS666 =
        _M0L7_2aSomeS665;
      int32_t _M0L3pslS1842 = _M0L14_2acurr__entryS666->$2;
      if (_M0L3pslS659 > _M0L3pslS1842) {
        int32_t _M0L3pslS1847;
        int32_t _M0L6_2atmpS1843;
        int32_t _M0L6_2atmpS1845;
        int32_t _M0L14capacity__maskS1846;
        int32_t _M0L6_2atmpS1844;
        _M0L5entryS661->$2 = _M0L3pslS659;
        moonbit_incref(_M0L14_2acurr__entryS666);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB5ArrayGiEE(_M0L4selfS663, _M0L5entryS661, _M0L3idxS660);
        moonbit_decref(_M0L5entryS661);
        _M0L3pslS1847 = _M0L14_2acurr__entryS666->$2;
        _M0L6_2atmpS1843 = _M0L3pslS1847 + 1;
        _M0L6_2atmpS1845 = _M0L3idxS660 + 1;
        _M0L14capacity__maskS1846 = _M0L4selfS663->$3;
        _M0L6_2atmpS1844 = _M0L6_2atmpS1845 & _M0L14capacity__maskS1846;
        _M0L3pslS659 = _M0L6_2atmpS1843;
        _M0L3idxS660 = _M0L6_2atmpS1844;
        _M0L5entryS661 = _M0L14_2acurr__entryS666;
        continue;
      } else {
        int32_t _M0L6_2atmpS1848 = _M0L3pslS659 + 1;
        int32_t _M0L6_2atmpS1850 = _M0L3idxS660 + 1;
        int32_t _M0L14capacity__maskS1851 = _M0L4selfS663->$3;
        int32_t _M0L6_2atmpS1849 =
          _M0L6_2atmpS1850 & _M0L14capacity__maskS1851;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _tmp_3109 = _M0L5entryS661;
        _M0L3pslS659 = _M0L6_2atmpS1848;
        _M0L3idxS660 = _M0L6_2atmpS1849;
        _M0L5entryS661 = _tmp_3109;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS637,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS639,
  int32_t _M0L8new__idxS638
) {
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1822;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1823;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS2795;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS640;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS1822 = _M0L4selfS637->$0;
  _M0L6_2atmpS1823 = _M0L5entryS639;
  if (
    _M0L8new__idxS638 < 0
    || _M0L8new__idxS638 >= Moonbit_array_length(_M0L7entriesS1822)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS2795
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1822[
      _M0L8new__idxS638
    ];
  if (_M0L6_2atmpS1823) {
    moonbit_incref(_M0L6_2atmpS1823);
  }
  if (_M0L6_2aoldS2795) {
    moonbit_decref(_M0L6_2aoldS2795);
  }
  _M0L7entriesS1822[_M0L8new__idxS638] = _M0L6_2atmpS1823;
  _M0L7_2abindS640 = _M0L5entryS639->$1;
  if (_M0L7_2abindS640 == 0) {
    _M0L4selfS637->$6 = _M0L8new__idxS638;
  } else {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS641 =
      _M0L7_2abindS640;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2anextS642 =
      _M0L7_2aSomeS641;
    _M0L7_2anextS642->$0 = _M0L8new__idxS638;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS643,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS645,
  int32_t _M0L8new__idxS644
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1824;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1825;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2798;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS646;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS1824 = _M0L4selfS643->$0;
  _M0L6_2atmpS1825 = _M0L5entryS645;
  if (
    _M0L8new__idxS644 < 0
    || _M0L8new__idxS644 >= Moonbit_array_length(_M0L7entriesS1824)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS2798
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1824[_M0L8new__idxS644];
  if (_M0L6_2atmpS1825) {
    moonbit_incref(_M0L6_2atmpS1825);
  }
  if (_M0L6_2aoldS2798) {
    moonbit_decref(_M0L6_2aoldS2798);
  }
  _M0L7entriesS1824[_M0L8new__idxS644] = _M0L6_2atmpS1825;
  _M0L7_2abindS646 = _M0L5entryS645->$1;
  if (_M0L7_2abindS646 == 0) {
    _M0L4selfS643->$6 = _M0L8new__idxS644;
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS647 = _M0L7_2abindS646;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2anextS648 = _M0L7_2aSomeS647;
    _M0L7_2anextS648->$0 = _M0L8new__idxS644;
  }
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS630,
  int32_t _M0L3idxS632,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS631
) {
  int32_t _M0L7_2abindS629;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1809;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1810;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS2800;
  int32_t _M0L4sizeS1812;
  int32_t _M0L6_2atmpS1811;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS629 = _M0L4selfS630->$6;
  switch (_M0L7_2abindS629) {
    case -1: {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1804 =
        _M0L5entryS631;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS2802 =
        _M0L4selfS630->$5;
      if (_M0L6_2atmpS1804) {
        moonbit_incref(_M0L6_2atmpS1804);
      }
      if (_M0L6_2aoldS2802) {
        moonbit_decref(_M0L6_2aoldS2802);
      }
      _M0L4selfS630->$5 = _M0L6_2atmpS1804;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS1808 =
        _M0L4selfS630->$0;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1807;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1805;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1806;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS2803;
      if (
        _M0L7_2abindS629 < 0
        || _M0L7_2abindS629 >= Moonbit_array_length(_M0L7entriesS1808)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS1807
      = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1808[
          _M0L7_2abindS629
        ];
      if (_M0L6_2atmpS1807) {
        moonbit_incref(_M0L6_2atmpS1807);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS1805
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRP26biolab8bio__seq9SeqRecordEE(_M0L6_2atmpS1807);
      if (_M0L6_2atmpS1807) {
        moonbit_decref(_M0L6_2atmpS1807);
      }
      _M0L6_2atmpS1806 = _M0L5entryS631;
      _M0L6_2aoldS2803 = _M0L6_2atmpS1805->$1;
      if (_M0L6_2atmpS1806) {
        moonbit_incref(_M0L6_2atmpS1806);
      }
      if (_M0L6_2aoldS2803) {
        moonbit_decref(_M0L6_2aoldS2803);
      }
      _M0L6_2atmpS1805->$1 = _M0L6_2atmpS1806;
      moonbit_decref(_M0L6_2atmpS1805);
      break;
    }
  }
  _M0L4selfS630->$6 = _M0L3idxS632;
  _M0L7entriesS1809 = _M0L4selfS630->$0;
  _M0L6_2atmpS1810 = _M0L5entryS631;
  if (
    _M0L3idxS632 < 0
    || _M0L3idxS632 >= Moonbit_array_length(_M0L7entriesS1809)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS2800
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS1809[
      _M0L3idxS632
    ];
  if (_M0L6_2atmpS1810) {
    moonbit_incref(_M0L6_2atmpS1810);
  }
  if (_M0L6_2aoldS2800) {
    moonbit_decref(_M0L6_2aoldS2800);
  }
  _M0L7entriesS1809[_M0L3idxS632] = _M0L6_2atmpS1810;
  _M0L4sizeS1812 = _M0L4selfS630->$1;
  _M0L6_2atmpS1811 = _M0L4sizeS1812 + 1;
  _M0L4selfS630->$1 = _M0L6_2atmpS1811;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS634,
  int32_t _M0L3idxS636,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS635
) {
  int32_t _M0L7_2abindS633;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1818;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1819;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2806;
  int32_t _M0L4sizeS1821;
  int32_t _M0L6_2atmpS1820;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS633 = _M0L4selfS634->$6;
  switch (_M0L7_2abindS633) {
    case -1: {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1813 = _M0L5entryS635;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2808 =
        _M0L4selfS634->$5;
      if (_M0L6_2atmpS1813) {
        moonbit_incref(_M0L6_2atmpS1813);
      }
      if (_M0L6_2aoldS2808) {
        moonbit_decref(_M0L6_2aoldS2808);
      }
      _M0L4selfS634->$5 = _M0L6_2atmpS1813;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1817 =
        _M0L4selfS634->$0;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1816;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1814;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1815;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS2809;
      if (
        _M0L7_2abindS633 < 0
        || _M0L7_2abindS633 >= Moonbit_array_length(_M0L7entriesS1817)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS1816
      = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1817[
          _M0L7_2abindS633
        ];
      if (_M0L6_2atmpS1816) {
        moonbit_incref(_M0L6_2atmpS1816);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS1814
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(_M0L6_2atmpS1816);
      if (_M0L6_2atmpS1816) {
        moonbit_decref(_M0L6_2atmpS1816);
      }
      _M0L6_2atmpS1815 = _M0L5entryS635;
      _M0L6_2aoldS2809 = _M0L6_2atmpS1814->$1;
      if (_M0L6_2atmpS1815) {
        moonbit_incref(_M0L6_2atmpS1815);
      }
      if (_M0L6_2aoldS2809) {
        moonbit_decref(_M0L6_2aoldS2809);
      }
      _M0L6_2atmpS1814->$1 = _M0L6_2atmpS1815;
      moonbit_decref(_M0L6_2atmpS1814);
      break;
    }
  }
  _M0L4selfS634->$6 = _M0L3idxS636;
  _M0L7entriesS1818 = _M0L4selfS634->$0;
  _M0L6_2atmpS1819 = _M0L5entryS635;
  if (
    _M0L3idxS636 < 0
    || _M0L3idxS636 >= Moonbit_array_length(_M0L7entriesS1818)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS2806
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1818[_M0L3idxS636];
  if (_M0L6_2atmpS1819) {
    moonbit_incref(_M0L6_2atmpS1819);
  }
  if (_M0L6_2aoldS2806) {
    moonbit_decref(_M0L6_2aoldS2806);
  }
  _M0L7entriesS1818[_M0L3idxS636] = _M0L6_2atmpS1819;
  _M0L4sizeS1821 = _M0L4selfS634->$1;
  _M0L6_2atmpS1820 = _M0L4sizeS1821 + 1;
  _M0L4selfS634->$1 = _M0L6_2atmpS1820;
  return 0;
}

int32_t _M0MPC13int3Int3max(int32_t _M0L4selfS627, int32_t _M0L5otherS628) {
  #line 75 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS627 > _M0L5otherS628) {
    return _M0L4selfS627;
  } else {
    return _M0L5otherS628;
  }
}

int32_t _M0FPB21capacity__for__length(int32_t _M0L6lengthS626) {
  int32_t _M0Lm8capacityS625;
  int32_t _M0L6_2atmpS1802;
  int32_t _M0L6_2atmpS1801;
  #line 71 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 72 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0Lm8capacityS625 = _M0MPC13int3Int20next__power__of__two(_M0L6lengthS626);
  _M0L6_2atmpS1802 = _M0Lm8capacityS625;
  #line 73 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1801 = _M0FPB21calc__grow__threshold(_M0L6_2atmpS1802);
  if (_M0L6lengthS626 > _M0L6_2atmpS1801) {
    int32_t _M0L6_2atmpS1803 = _M0Lm8capacityS625;
    _M0Lm8capacityS625 = _M0L6_2atmpS1803 * 2;
  }
  return _M0Lm8capacityS625;
}

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0FPB8new__mapGsRPB5ArrayGiEE(
  int32_t _M0L8capacityS614
) {
  int32_t _M0L8capacityS613;
  int32_t _M0L7_2abindS615;
  int32_t _M0L7_2abindS616;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1799;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7_2abindS617;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS618;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _block_3110;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS613
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS614);
  _M0L7_2abindS615 = _M0L8capacityS613 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS616 = _M0FPB21calc__grow__threshold(_M0L8capacityS613);
  _M0L6_2atmpS1799 = 0;
  _M0L7_2abindS617
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE**)moonbit_make_ref_array(_M0L8capacityS613, _M0L6_2atmpS1799);
  _M0L7_2abindS618 = 0;
  _block_3110
  = (struct _M0TPB3MapGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB5ArrayGiEE));
  Moonbit_object_header(_block_3110)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 70, 0);
  _block_3110->$0 = _M0L7_2abindS617;
  _block_3110->$1 = 0;
  _block_3110->$2 = _M0L8capacityS613;
  _block_3110->$3 = _M0L7_2abindS615;
  _block_3110->$4 = _M0L7_2abindS616;
  _block_3110->$5 = _M0L7_2abindS618;
  _block_3110->$6 = -1;
  return _block_3110;
}

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FPB8new__mapGsRP26biolab8bio__seq9SeqRecordE(
  int32_t _M0L8capacityS620
) {
  int32_t _M0L8capacityS619;
  int32_t _M0L7_2abindS621;
  int32_t _M0L7_2abindS622;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS1800;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7_2abindS623;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS624;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _block_3111;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS619
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS620);
  _M0L7_2abindS621 = _M0L8capacityS619 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS622 = _M0FPB21calc__grow__threshold(_M0L8capacityS619);
  _M0L6_2atmpS1800 = 0;
  _M0L7_2abindS623
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE**)moonbit_make_ref_array(_M0L8capacityS619, _M0L6_2atmpS1800);
  _M0L7_2abindS624 = 0;
  _block_3111
  = (struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_block_3111)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 74, 0);
  _block_3111->$0 = _M0L7_2abindS623;
  _block_3111->$1 = 0;
  _block_3111->$2 = _M0L8capacityS619;
  _block_3111->$3 = _M0L7_2abindS621;
  _block_3111->$4 = _M0L7_2abindS622;
  _block_3111->$5 = _M0L7_2abindS624;
  _block_3111->$6 = -1;
  return _block_3111;
}

int32_t _M0MPC13int3Int20next__power__of__two(int32_t _M0L4selfS612) {
  #line 33 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS612 >= 0) {
    int32_t _M0L6_2atmpS1798;
    int32_t _M0L6_2atmpS1797;
    int32_t _M0L6_2atmpS1796;
    int32_t _M0L6_2atmpS1795;
    if (_M0L4selfS612 <= 1) {
      return 1;
    }
    if (_M0L4selfS612 > 1073741824) {
      return 1073741824;
    }
    _M0L6_2atmpS1798 = _M0L4selfS612 - 1;
    #line 44 "/home/developer/.moon/lib/core/builtin/int.mbt"
    _M0L6_2atmpS1797 = moonbit_clz32(_M0L6_2atmpS1798);
    _M0L6_2atmpS1796 = _M0L6_2atmpS1797 - 1;
    _M0L6_2atmpS1795 = 2147483647 >> (_M0L6_2atmpS1796 & 31);
    return _M0L6_2atmpS1795 + 1;
  } else {
    #line 34 "/home/developer/.moon/lib/core/builtin/int.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB21calc__grow__threshold(int32_t _M0L8capacityS611) {
  int32_t _M0L6_2atmpS1794;
  #line 610 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1794 = _M0L8capacityS611 * 13;
  return _M0L6_2atmpS1794 / 16;
}

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0MPC16option6Option6unwrapGRPB5EntryGsRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS607
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS607 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS608 =
      _M0L4selfS607;
    if (_M0L7_2aSomeS608) {
      moonbit_incref(_M0L7_2aSomeS608);
    }
    return _M0L7_2aSomeS608;
  }
}

struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4selfS609
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS609 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS610 = _M0L4selfS609;
    if (_M0L7_2aSomeS610) {
      moonbit_incref(_M0L7_2aSomeS610);
    }
    return _M0L7_2aSomeS610;
  }
}

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t _M0L4selfS606) {
  #line 35 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 36 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC13int3Int18to__string_2einner(_M0L4selfS606, 10);
}

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPB4Iter9to__arrayGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS598
) {
  int64_t _M0L7_2abindS597;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L6resultS599;
  #line 841 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS597 = _M0L4selfS598->$1;
  if (_M0L7_2abindS597 == 4294967296ll) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS1793 =
      (struct _M0TPC16string10StringView*)moonbit_empty_ref_valtype_array;
    _M0L6resultS599
    = (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_M0L6resultS599)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 78, 0);
    _M0L6resultS599->$0 = _M0L6_2atmpS1793;
    _M0L6resultS599->$1 = 0;
  } else {
    int64_t _M0L7_2aSomeS600 = _M0L7_2abindS597;
    int32_t _M0L4_2anS601 = (int32_t)_M0L7_2aSomeS600;
    #line 844 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6resultS599
    = _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(_M0L4_2anS601);
  }
  _2awhile_605:;
  while (1) {
    void* _M0L7_2abindS602;
    #line 847 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L7_2abindS602
    = _M0MPB4Iter4nextGRPC16string10StringViewE(_M0L4selfS598);
    switch (Moonbit_object_tag(_M0L7_2abindS602)) {
      case 1: {
        struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS603 =
          (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS602;
        struct _M0TPC16string10StringView _M0L8_2afieldS2812 =
          _M0L7_2aSomeS603->$0;
        int32_t _M0L6_2acntS2930 =
          Moonbit_rc_count(Moonbit_object_header(_M0L7_2aSomeS603));
        struct _M0TPC16string10StringView _M0L4_2axS604;
        if (_M0L6_2acntS2930 > 1) {
          int32_t _M0L11_2anew__cntS2931 = _M0L6_2acntS2930 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2aSomeS603), _M0L11_2anew__cntS2931);
          moonbit_incref(_M0L8_2afieldS2812.$0);
        } else if (_M0L6_2acntS2930 == 1) {
          #line 847 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
          moonbit_free(_M0L7_2aSomeS603);
        }
        _M0L4_2axS604 = _M0L8_2afieldS2812;
        #line 848 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
        _M0MPC15array5Array4pushGRPC16string10StringViewE(_M0L6resultS599, _M0L4_2axS604);
        moonbit_decref(_M0L4_2axS604.$0);
        goto _2awhile_605;
        break;
      }
      default: {
        moonbit_decref(_M0L7_2abindS602);
        break;
      }
    }
    break;
  }
  return _M0L6resultS599;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string6String5split(
  moonbit_string_t _M0L4selfS595,
  struct _M0TPC16string10StringView _M0L3sepS596
) {
  int32_t _M0L6_2atmpS1792;
  struct _M0TPC16string10StringView _M0L6_2atmpS1791;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_3113;
  #line 1168 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1792 = Moonbit_array_length(_M0L4selfS595);
  moonbit_incref(_M0L4selfS595);
  _M0L6_2atmpS1791
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS595, .$1 = 0, .$2 = _M0L6_2atmpS1792
  };
  #line 1169 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3113
  = _M0MPC16string10StringView5split(_M0L6_2atmpS1791, _M0L3sepS596);
  moonbit_decref(_M0L6_2atmpS1791.$0);
  return _result_3113;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string10StringView5split(
  struct _M0TPC16string10StringView _M0L4selfS586,
  struct _M0TPC16string10StringView _M0L3sepS585
) {
  int32_t _M0L3endS1789;
  int32_t _M0L5startS1790;
  int32_t _M0L8sep__lenS584;
  void* _M0L4SomeS1788;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS588;
  struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__* _closure_3115;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS1778;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_3116;
  #line 1139 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1789 = _M0L3sepS585.$2;
  _M0L5startS1790 = _M0L3sepS585.$1;
  _M0L8sep__lenS584 = _M0L3endS1789 - _M0L5startS1790;
  if (_M0L8sep__lenS584 == 0) {
    struct _M0TPB4IterGcE* _M0L6_2atmpS1773;
    struct _M0TWcERPC16string10StringView* _M0L6_2atmpS1774;
    struct _M0TPB4IterGRPC16string10StringViewE* _result_3114;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1773 = _M0MPC16string10StringView4iter(_M0L4selfS586);
    _M0L6_2atmpS1774
    = (struct _M0TWcERPC16string10StringView*)&_M0MPC16string10StringView5splitC1775l1145$closure.data;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_3114
    = _M0MPB4Iter3mapGcRPC16string10StringViewE(_M0L6_2atmpS1773, _M0L6_2atmpS1774);
    moonbit_decref(_M0L6_2atmpS1773);
    moonbit_decref(_M0L6_2atmpS1774);
    return _result_3114;
  }
  moonbit_incref(_M0L4selfS586.$0);
  _M0L4SomeS1788
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
  Moonbit_object_header(_M0L4SomeS1788)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 1);
  ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS1788)->$0
  = _M0L4selfS586;
  _M0L9remainingS588
  = (struct _M0TPB8MutLocalGORPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE));
  Moonbit_object_header(_M0L9remainingS588)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 87, 0);
  _M0L9remainingS588->$0 = _M0L4SomeS1788;
  moonbit_incref(_M0L3sepS585.$0);
  _closure_3115
  = (struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__*)moonbit_malloc(sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__));
  Moonbit_object_header(_closure_3115)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 90, 0);
  _closure_3115->code = &_M0MPC16string10StringView5splitC1779l1148;
  _closure_3115->$0 = _M0L9remainingS588;
  _closure_3115->$1 = _M0L3sepS585;
  _closure_3115->$2 = _M0L8sep__lenS584;
  _M0L6_2atmpS1778
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_3115;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3116
  = _M0MPB4Iter3newGRPC16string10StringViewE(_M0L6_2atmpS1778, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1778);
  return _result_3116;
}

void* _M0MPC16string10StringView5splitC1779l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS1780
) {
  struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__* _M0L14_2acasted__envS1781;
  int32_t _M0L8sep__lenS584;
  struct _M0TPC16string10StringView _M0L3sepS585;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS588;
  void* _M0L7_2abindS589;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L14_2acasted__envS1781
  = (struct _M0R44StringView_3a_3asplit_2eanon__u1779__l1148__*)_M0L6_2aenvS1780;
  _M0L8sep__lenS584 = _M0L14_2acasted__envS1781->$2;
  _M0L3sepS585 = _M0L14_2acasted__envS1781->$1;
  _M0L9remainingS588 = _M0L14_2acasted__envS1781->$0;
  _M0L7_2abindS589 = _M0L9remainingS588->$0;
  switch (Moonbit_object_tag(_M0L7_2abindS589)) {
    case 1: {
      struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS590 =
        (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS589;
      struct _M0TPC16string10StringView _M0L7_2aviewS591 =
        _M0L7_2aSomeS590->$0;
      int64_t _M0L7_2abindS592;
      moonbit_incref(_M0L7_2aviewS591.$0);
      #line 1150 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS592
      = _M0MPC16string10StringView4find(_M0L7_2aviewS591, _M0L3sepS585);
      if (_M0L7_2abindS592 == 4294967296ll) {
        void* _M0L4NoneS1782 =
          (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
        void* _M0L6_2aoldS2813 = _M0L9remainingS588->$0;
        void* _block_3117;
        moonbit_decref(_M0L6_2aoldS2813);
        _M0L9remainingS588->$0 = _M0L4NoneS1782;
        _block_3117
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_3117)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3117)->$0
        = _M0L7_2aviewS591;
        return _block_3117;
      } else {
        int64_t _M0L7_2aSomeS593 = _M0L7_2abindS592;
        int32_t _M0L6_2aendS594 = (int32_t)_M0L7_2aSomeS593;
        int32_t _M0L6_2atmpS1785 = _M0L6_2aendS594 + _M0L8sep__lenS584;
        struct _M0TPC16string10StringView _M0L6_2atmpS1784;
        void* _M0L4SomeS1783;
        void* _M0L6_2aoldS2814;
        int64_t _M0L6_2atmpS1787;
        struct _M0TPC16string10StringView _M0L6_2atmpS1786;
        void* _block_3118;
        #line 1154 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS1784
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS591, _M0L6_2atmpS1785, 4294967296ll);
        _M0L4SomeS1783
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_M0L4SomeS1783)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS1783)->$0
        = _M0L6_2atmpS1784;
        _M0L6_2aoldS2814 = _M0L9remainingS588->$0;
        moonbit_decref(_M0L6_2aoldS2814);
        _M0L9remainingS588->$0 = _M0L4SomeS1783;
        _M0L6_2atmpS1787 = (int64_t)_M0L6_2aendS594;
        #line 1155 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS1786
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS591, 0, _M0L6_2atmpS1787);
        moonbit_decref(_M0L7_2aviewS591.$0);
        _block_3118
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_3118)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3118)->$0
        = _M0L6_2atmpS1786;
        return _block_3118;
      }
      break;
    }
    default: {
      return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
      break;
    }
  }
}

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC1775l1145(
  struct _M0TWcERPC16string10StringView* _M0L6_2aenvS1776,
  int32_t _M0L1cS587
) {
  moonbit_string_t _M0L6_2atmpS1777;
  struct _M0TPC16string10StringView _result_3119;
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1777 = _M0IPC14char4CharPB4Show10to__string(_M0L1cS587);
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3119
  = _M0MPC16string6String12view_2einner(_M0L6_2atmpS1777, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1777);
  return _result_3119;
}

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t _M0L4selfS583) {
  #line 446 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 447 "/home/developer/.moon/lib/core/builtin/char.mbt"
  return _M0FPB16char__to__string(_M0L4selfS583);
}

moonbit_string_t _M0FPB16char__to__string(int32_t _M0L4charS582) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS581;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1772;
  moonbit_string_t _result_3120;
  #line 452 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0L7_2aselfS581 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS581, _M0L4charS582);
  _M0L6_2atmpS1772 = _M0L7_2aselfS581;
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _result_3120 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1772);
  moonbit_decref(_M0L6_2atmpS1772);
  return _result_3120;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE* _M0L4selfS577,
  struct _M0TWcERPC16string10StringView* _M0L1fS580
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__* _closure_3121;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS1766;
  int64_t _M0L10size__hintS1767;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_3122;
  #line 389 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  moonbit_incref(_M0L1fS580);
  moonbit_incref(_M0L4selfS577);
  _closure_3121
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__*)moonbit_malloc(sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__));
  Moonbit_object_header(_closure_3121)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 94, 0);
  _closure_3121->code = &_M0MPB4Iter3mapGcRPC16string10StringViewEC1768l391;
  _closure_3121->$0 = _M0L1fS580;
  _closure_3121->$1 = _M0L4selfS577;
  _M0L6_2atmpS1766
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_3121;
  _M0L10size__hintS1767 = _M0L4selfS577->$1;
  _block_3122
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_3122)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 98, 0);
  _block_3122->$0 = _M0L6_2atmpS1766;
  _block_3122->$1 = _M0L10size__hintS1767;
  return _block_3122;
}

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC1768l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS1769
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__* _M0L14_2acasted__envS1770;
  struct _M0TPB4IterGcE* _M0L4selfS577;
  struct _M0TWcERPC16string10StringView* _M0L1fS580;
  int32_t _M0L7_2abindS576;
  #line 391 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L14_2acasted__envS1770
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1768__l391__*)_M0L6_2aenvS1769;
  _M0L4selfS577 = _M0L14_2acasted__envS1770->$1;
  _M0L1fS580 = _M0L14_2acasted__envS1770->$0;
  #line 392 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS576 = _M0MPB4Iter4nextGcE(_M0L4selfS577);
  if (_M0L7_2abindS576 == -1) {
    return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  } else {
    int32_t _M0L7_2aSomeS578 = _M0L7_2abindS576;
    int32_t _M0L4_2axS579 = _M0L7_2aSomeS578;
    struct _M0TPC16string10StringView _M0L6_2atmpS1771;
    void* _block_3123;
    #line 393 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6_2atmpS1771 = _M0L1fS580->code(_M0L1fS580, _M0L4_2axS579);
    _block_3123
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
    Moonbit_object_header(_block_3123)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 1);
    ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3123)->$0
    = _M0L6_2atmpS1771;
    return _block_3123;
  }
}

struct _M0TPB4IterGcE* _M0MPC16string6String4iter(
  moonbit_string_t _M0L4selfS571
) {
  int32_t _M0L3lenS570;
  struct _M0TPB8MutLocalGiE* _M0L5indexS572;
  struct _M0R38String_3a_3aiter_2eanon__u1750__l256__* _closure_3124;
  struct _M0TWEOc* _M0L6_2atmpS1749;
  struct _M0TPB4IterGcE* _result_3125;
  #line 252 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L3lenS570 = Moonbit_array_length(_M0L4selfS571);
  _M0L5indexS572
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS572)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS572->$0 = 0;
  moonbit_incref(_M0L4selfS571);
  _closure_3124
  = (struct _M0R38String_3a_3aiter_2eanon__u1750__l256__*)moonbit_malloc(sizeof(struct _M0R38String_3a_3aiter_2eanon__u1750__l256__));
  Moonbit_object_header(_closure_3124)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 101, 0);
  _closure_3124->code = &_M0MPC16string6String4iterC1750l256;
  _closure_3124->$0 = _M0L5indexS572;
  _closure_3124->$1 = _M0L4selfS571;
  _closure_3124->$2 = _M0L3lenS570;
  _M0L6_2atmpS1749 = (struct _M0TWEOc*)_closure_3124;
  #line 256 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_3125 = _M0MPB4Iter3newGcE(_M0L6_2atmpS1749, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1749);
  return _result_3125;
}

int32_t _M0MPC16string6String4iterC1750l256(
  struct _M0TWEOc* _M0L6_2aenvS1751
) {
  struct _M0R38String_3a_3aiter_2eanon__u1750__l256__* _M0L14_2acasted__envS1752;
  int32_t _M0L3lenS570;
  moonbit_string_t _M0L4selfS571;
  struct _M0TPB8MutLocalGiE* _M0L5indexS572;
  int32_t _M0L3valS1753;
  #line 256 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L14_2acasted__envS1752
  = (struct _M0R38String_3a_3aiter_2eanon__u1750__l256__*)_M0L6_2aenvS1751;
  _M0L3lenS570 = _M0L14_2acasted__envS1752->$2;
  _M0L4selfS571 = _M0L14_2acasted__envS1752->$1;
  _M0L5indexS572 = _M0L14_2acasted__envS1752->$0;
  _M0L3valS1753 = _M0L5indexS572->$0;
  if (_M0L3valS1753 < _M0L3lenS570) {
    int32_t _M0L3valS1765 = _M0L5indexS572->$0;
    int32_t _M0L2c1S573 = _M0L4selfS571[_M0L3valS1765];
    int32_t _if__result_3126;
    int32_t _M0L3valS1763;
    int32_t _M0L6_2atmpS1762;
    int32_t _M0L6_2atmpS1764;
    #line 259 "/home/developer/.moon/lib/core/builtin/string.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S573)) {
      int32_t _M0L3valS1755 = _M0L5indexS572->$0;
      int32_t _M0L6_2atmpS1754 = _M0L3valS1755 + 1;
      _if__result_3126 = _M0L6_2atmpS1754 < _M0L3lenS570;
    } else {
      _if__result_3126 = 0;
    }
    if (_if__result_3126) {
      int32_t _M0L3valS1761 = _M0L5indexS572->$0;
      int32_t _M0L6_2atmpS1760 = _M0L3valS1761 + 1;
      int32_t _M0L2c2S574 = _M0L4selfS571[_M0L6_2atmpS1760];
      #line 261 "/home/developer/.moon/lib/core/builtin/string.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S574)) {
        int32_t _M0L6_2atmpS1758 = (int32_t)_M0L2c1S573;
        int32_t _M0L6_2atmpS1759 = (int32_t)_M0L2c2S574;
        int32_t _M0L1cS575;
        int32_t _M0L3valS1757;
        int32_t _M0L6_2atmpS1756;
        #line 262 "/home/developer/.moon/lib/core/builtin/string.mbt"
        _M0L1cS575
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1758, _M0L6_2atmpS1759);
        _M0L3valS1757 = _M0L5indexS572->$0;
        _M0L6_2atmpS1756 = _M0L3valS1757 + 2;
        _M0L5indexS572->$0 = _M0L6_2atmpS1756;
        return _M0L1cS575;
      }
    }
    _M0L3valS1763 = _M0L5indexS572->$0;
    _M0L6_2atmpS1762 = _M0L3valS1763 + 1;
    _M0L5indexS572->$0 = _M0L6_2atmpS1762;
    #line 269 "/home/developer/.moon/lib/core/builtin/string.mbt"
    _M0L6_2atmpS1764 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S573);
    return _M0L6_2atmpS1764;
  } else {
    return -1;
  }
}

int32_t _M0MPC15array5Array4pushGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS555,
  int32_t _M0L5valueS557
) {
  int32_t _M0L3lenS1724;
  int32_t* _M0L6_2atmpS1726;
  int32_t _M0L6_2atmpS1725;
  int32_t _M0L6lengthS556;
  int32_t* _M0L3bufS1727;
  int32_t _M0L6_2atmpS1728;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1724 = _M0L4selfS555->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1726 = _M0MPC15array5Array6bufferGiE(_M0L4selfS555);
  _M0L6_2atmpS1725 = Moonbit_array_length(_M0L6_2atmpS1726);
  moonbit_decref(_M0L6_2atmpS1726);
  if (_M0L3lenS1724 == _M0L6_2atmpS1725) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGiE(_M0L4selfS555);
  }
  _M0L6lengthS556 = _M0L4selfS555->$1;
  _M0L3bufS1727 = _M0L4selfS555->$0;
  _M0L3bufS1727[_M0L6lengthS556] = _M0L5valueS557;
  _M0L6_2atmpS1728 = _M0L6lengthS556 + 1;
  _M0L4selfS555->$1 = _M0L6_2atmpS1728;
  return 0;
}

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS558,
  moonbit_string_t _M0L5valueS560
) {
  int32_t _M0L3lenS1729;
  moonbit_string_t* _M0L6_2atmpS1731;
  int32_t _M0L6_2atmpS1730;
  int32_t _M0L6lengthS559;
  moonbit_string_t* _M0L3bufS1732;
  moonbit_string_t _M0L6_2aoldS2818;
  int32_t _M0L6_2atmpS1733;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1729 = _M0L4selfS558->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1731 = _M0MPC15array5Array6bufferGsE(_M0L4selfS558);
  _M0L6_2atmpS1730 = Moonbit_array_length(_M0L6_2atmpS1731);
  moonbit_decref(_M0L6_2atmpS1731);
  if (_M0L3lenS1729 == _M0L6_2atmpS1730) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGsE(_M0L4selfS558);
  }
  _M0L6lengthS559 = _M0L4selfS558->$1;
  _M0L3bufS1732 = _M0L4selfS558->$0;
  _M0L6_2aoldS2818 = (moonbit_string_t)_M0L3bufS1732[_M0L6lengthS559];
  moonbit_incref(_M0L5valueS560);
  moonbit_decref(_M0L6_2aoldS2818);
  _M0L3bufS1732[_M0L6lengthS559] = _M0L5valueS560;
  _M0L6_2atmpS1733 = _M0L6lengthS559 + 1;
  _M0L4selfS558->$1 = _M0L6_2atmpS1733;
  return 0;
}

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS561,
  struct _M0TUsiE* _M0L5valueS563
) {
  int32_t _M0L3lenS1734;
  struct _M0TUsiE** _M0L6_2atmpS1736;
  int32_t _M0L6_2atmpS1735;
  int32_t _M0L6lengthS562;
  struct _M0TUsiE** _M0L3bufS1737;
  struct _M0TUsiE* _M0L6_2aoldS2820;
  int32_t _M0L6_2atmpS1738;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1734 = _M0L4selfS561->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1736 = _M0MPC15array5Array6bufferGUsiEE(_M0L4selfS561);
  _M0L6_2atmpS1735 = Moonbit_array_length(_M0L6_2atmpS1736);
  moonbit_decref(_M0L6_2atmpS1736);
  if (_M0L3lenS1734 == _M0L6_2atmpS1735) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGUsiEE(_M0L4selfS561);
  }
  _M0L6lengthS562 = _M0L4selfS561->$1;
  _M0L3bufS1737 = _M0L4selfS561->$0;
  _M0L6_2aoldS2820 = (struct _M0TUsiE*)_M0L3bufS1737[_M0L6lengthS562];
  moonbit_incref(_M0L5valueS563);
  if (_M0L6_2aoldS2820) {
    moonbit_decref(_M0L6_2aoldS2820);
  }
  _M0L3bufS1737[_M0L6lengthS562] = _M0L5valueS563;
  _M0L6_2atmpS1738 = _M0L6lengthS562 + 1;
  _M0L4selfS561->$1 = _M0L6_2atmpS1738;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS564,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS566
) {
  int32_t _M0L3lenS1739;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS1741;
  int32_t _M0L6_2atmpS1740;
  int32_t _M0L6lengthS565;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS1742;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS2822;
  int32_t _M0L6_2atmpS1743;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1739 = _M0L4selfS564->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1741
  = _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS564);
  _M0L6_2atmpS1740 = Moonbit_array_length(_M0L6_2atmpS1741);
  moonbit_decref(_M0L6_2atmpS1741);
  if (_M0L3lenS1739 == _M0L6_2atmpS1740) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS564);
  }
  _M0L6lengthS565 = _M0L4selfS564->$1;
  _M0L3bufS1742 = _M0L4selfS564->$0;
  _M0L6_2aoldS2822
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS1742[_M0L6lengthS565];
  moonbit_incref(_M0L5valueS566);
  if (_M0L6_2aoldS2822) {
    moonbit_decref(_M0L6_2aoldS2822);
  }
  _M0L3bufS1742[_M0L6lengthS565] = _M0L5valueS566;
  _M0L6_2atmpS1743 = _M0L6lengthS565 + 1;
  _M0L4selfS564->$1 = _M0L6_2atmpS1743;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS567,
  struct _M0TPC16string10StringView _M0L5valueS569
) {
  int32_t _M0L3lenS1744;
  struct _M0TPC16string10StringView* _M0L6_2atmpS1746;
  int32_t _M0L6_2atmpS1745;
  int32_t _M0L6lengthS568;
  struct _M0TPC16string10StringView* _M0L3bufS1747;
  struct _M0TPC16string10StringView _M0L6_2aoldS2824;
  int32_t _M0L6_2atmpS1748;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1744 = _M0L4selfS567->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1746
  = _M0MPC15array5Array6bufferGRPC16string10StringViewE(_M0L4selfS567);
  _M0L6_2atmpS1745 = Moonbit_array_length(_M0L6_2atmpS1746);
  moonbit_decref(_M0L6_2atmpS1746);
  if (_M0L3lenS1744 == _M0L6_2atmpS1745) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRPC16string10StringViewE(_M0L4selfS567);
  }
  _M0L6lengthS568 = _M0L4selfS567->$1;
  _M0L3bufS1747 = _M0L4selfS567->$0;
  _M0L6_2aoldS2824 = _M0L3bufS1747[_M0L6lengthS568];
  moonbit_incref(_M0L5valueS569.$0);
  moonbit_decref(_M0L6_2aoldS2824.$0);
  _M0L3bufS1747[_M0L6lengthS568] = _M0L5valueS569;
  _M0L6_2atmpS1748 = _M0L6lengthS568 + 1;
  _M0L4selfS567->$1 = _M0L6_2atmpS1748;
  return 0;
}

int32_t _M0MPC15array5Array7reallocGiE(struct _M0TPB5ArrayGiE* _M0L4selfS541) {
  int32_t _M0L8old__capS540;
  int32_t _M0L8new__capS542;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS540 = _M0L4selfS541->$1;
  if (_M0L8old__capS540 == 0) {
    _M0L8new__capS542 = 8;
  } else {
    _M0L8new__capS542 = _M0L8old__capS540 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGiE(_M0L4selfS541, _M0L8new__capS542);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE* _M0L4selfS544) {
  int32_t _M0L8old__capS543;
  int32_t _M0L8new__capS545;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS543 = _M0L4selfS544->$1;
  if (_M0L8old__capS543 == 0) {
    _M0L8new__capS545 = 8;
  } else {
    _M0L8new__capS545 = _M0L8old__capS543 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGsE(_M0L4selfS544, _M0L8new__capS545);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS547
) {
  int32_t _M0L8old__capS546;
  int32_t _M0L8new__capS548;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS546 = _M0L4selfS547->$1;
  if (_M0L8old__capS546 == 0) {
    _M0L8new__capS548 = 8;
  } else {
    _M0L8new__capS548 = _M0L8old__capS546 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGUsiEE(_M0L4selfS547, _M0L8new__capS548);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS550
) {
  int32_t _M0L8old__capS549;
  int32_t _M0L8new__capS551;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS549 = _M0L4selfS550->$1;
  if (_M0L8old__capS549 == 0) {
    _M0L8new__capS551 = 8;
  } else {
    _M0L8new__capS551 = _M0L8old__capS549 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS550, _M0L8new__capS551);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS553
) {
  int32_t _M0L8old__capS552;
  int32_t _M0L8new__capS554;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS552 = _M0L4selfS553->$1;
  if (_M0L8old__capS552 == 0) {
    _M0L8new__capS554 = 8;
  } else {
    _M0L8new__capS554 = _M0L8old__capS552 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(_M0L4selfS553, _M0L8new__capS554);
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS511,
  int32_t _M0L13new__capacityS514
) {
  int32_t* _M0L8old__bufS510;
  int32_t _M0L8old__capS512;
  int32_t _M0L9copy__lenS513;
  int32_t* _M0L8new__bufS515;
  int32_t* _M0L6_2aoldS2826;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS510 = _M0L4selfS511->$0;
  _M0L8old__capS512 = Moonbit_array_length(_M0L8old__bufS510);
  if (_M0L8old__capS512 < _M0L13new__capacityS514) {
    _M0L9copy__lenS513 = _M0L8old__capS512;
  } else {
    _M0L9copy__lenS513 = _M0L13new__capacityS514;
  }
  moonbit_incref(_M0L8old__bufS510);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS515
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGiE(_M0L8old__bufS510, _M0L13new__capacityS514, _M0L9copy__lenS513, 0, 0);
  moonbit_decref(_M0L8old__bufS510);
  _M0L6_2aoldS2826 = _M0L4selfS511->$0;
  moonbit_decref(_M0L6_2aoldS2826);
  _M0L4selfS511->$0 = _M0L8new__bufS515;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS517,
  int32_t _M0L13new__capacityS520
) {
  moonbit_string_t* _M0L8old__bufS516;
  int32_t _M0L8old__capS518;
  int32_t _M0L9copy__lenS519;
  moonbit_string_t* _M0L8new__bufS521;
  moonbit_string_t* _M0L6_2aoldS2828;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS516 = _M0L4selfS517->$0;
  _M0L8old__capS518 = Moonbit_array_length(_M0L8old__bufS516);
  if (_M0L8old__capS518 < _M0L13new__capacityS520) {
    _M0L9copy__lenS519 = _M0L8old__capS518;
  } else {
    _M0L9copy__lenS519 = _M0L13new__capacityS520;
  }
  moonbit_incref(_M0L8old__bufS516);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS521
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(_M0L8old__bufS516, _M0L13new__capacityS520, _M0L9copy__lenS519, 0, 0);
  moonbit_decref(_M0L8old__bufS516);
  _M0L6_2aoldS2828 = _M0L4selfS517->$0;
  moonbit_decref(_M0L6_2aoldS2828);
  _M0L4selfS517->$0 = _M0L8new__bufS521;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS523,
  int32_t _M0L13new__capacityS526
) {
  struct _M0TUsiE** _M0L8old__bufS522;
  int32_t _M0L8old__capS524;
  int32_t _M0L9copy__lenS525;
  struct _M0TUsiE** _M0L8new__bufS527;
  struct _M0TUsiE** _M0L6_2aoldS2830;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS522 = _M0L4selfS523->$0;
  _M0L8old__capS524 = Moonbit_array_length(_M0L8old__bufS522);
  if (_M0L8old__capS524 < _M0L13new__capacityS526) {
    _M0L9copy__lenS525 = _M0L8old__capS524;
  } else {
    _M0L9copy__lenS525 = _M0L13new__capacityS526;
  }
  moonbit_incref(_M0L8old__bufS522);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS527
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(_M0L8old__bufS522, _M0L13new__capacityS526, _M0L9copy__lenS525, 0, 0);
  moonbit_decref(_M0L8old__bufS522);
  _M0L6_2aoldS2830 = _M0L4selfS523->$0;
  moonbit_decref(_M0L6_2aoldS2830);
  _M0L4selfS523->$0 = _M0L8new__bufS527;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS529,
  int32_t _M0L13new__capacityS532
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8old__bufS528;
  int32_t _M0L8old__capS530;
  int32_t _M0L9copy__lenS531;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8new__bufS533;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2aoldS2832;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS528 = _M0L4selfS529->$0;
  _M0L8old__capS530 = Moonbit_array_length(_M0L8old__bufS528);
  if (_M0L8old__capS530 < _M0L13new__capacityS532) {
    _M0L9copy__lenS531 = _M0L8old__capS530;
  } else {
    _M0L9copy__lenS531 = _M0L13new__capacityS532;
  }
  moonbit_incref(_M0L8old__bufS528);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS533
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(_M0L8old__bufS528, _M0L13new__capacityS532, _M0L9copy__lenS531, 0, 0);
  moonbit_decref(_M0L8old__bufS528);
  _M0L6_2aoldS2832 = _M0L4selfS529->$0;
  moonbit_decref(_M0L6_2aoldS2832);
  _M0L4selfS529->$0 = _M0L8new__bufS533;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS535,
  int32_t _M0L13new__capacityS538
) {
  struct _M0TPC16string10StringView* _M0L8old__bufS534;
  int32_t _M0L8old__capS536;
  int32_t _M0L9copy__lenS537;
  struct _M0TPC16string10StringView* _M0L8new__bufS539;
  struct _M0TPC16string10StringView* _M0L6_2aoldS2834;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS534 = _M0L4selfS535->$0;
  _M0L8old__capS536 = Moonbit_array_length(_M0L8old__bufS534);
  if (_M0L8old__capS536 < _M0L13new__capacityS538) {
    _M0L9copy__lenS537 = _M0L8old__capS536;
  } else {
    _M0L9copy__lenS537 = _M0L13new__capacityS538;
  }
  moonbit_incref(_M0L8old__bufS534);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS539
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(_M0L8old__bufS534, _M0L13new__capacityS538, _M0L9copy__lenS537, 0, 0);
  moonbit_decref(_M0L8old__bufS534);
  _M0L6_2aoldS2834 = _M0L4selfS535->$0;
  moonbit_decref(_M0L6_2aoldS2834);
  _M0L4selfS535->$0 = _M0L8new__bufS539;
  return 0;
}

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS507
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS507->$1;
}

int32_t _M0MPC15array5Array6lengthGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS508
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS508->$1;
}

int32_t _M0MPC15array5Array6lengthGiE(struct _M0TPB5ArrayGiE* _M0L4selfS509) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS509->$1;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(
  int32_t _M0L8capacityS503
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS503 == 0) {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS1716 =
      (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _block_3127 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
    Moonbit_object_header(_block_3127)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
    _block_3127->$0 = _M0L6_2atmpS1716;
    _block_3127->$1 = 0;
    return _block_3127;
  } else {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS1717 =
      (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array(_M0L8capacityS503, 0);
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _block_3128 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
    Moonbit_object_header(_block_3128)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
    _block_3128->$0 = _M0L6_2atmpS1717;
    _block_3128->$1 = 0;
    return _block_3128;
  }
}

struct _M0TPB5ArrayGiE* _M0MPC15array5Array11new_2einnerGiE(
  int32_t _M0L8capacityS504
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS504 == 0) {
    int32_t* _M0L6_2atmpS1718 = (int32_t*)moonbit_empty_int32_array;
    struct _M0TPB5ArrayGiE* _block_3129 =
      (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
    Moonbit_object_header(_block_3129)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
    _block_3129->$0 = _M0L6_2atmpS1718;
    _block_3129->$1 = 0;
    return _block_3129;
  } else {
    int32_t* _M0L6_2atmpS1719 =
      (int32_t*)moonbit_make_int32_array_raw(_M0L8capacityS504);
    struct _M0TPB5ArrayGiE* _block_3130 =
      (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
    Moonbit_object_header(_block_3130)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
    _block_3130->$0 = _M0L6_2atmpS1719;
    _block_3130->$1 = 0;
    return _block_3130;
  }
}

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(
  int32_t _M0L8capacityS505
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS505 == 0) {
    moonbit_string_t* _M0L6_2atmpS1720 =
      (moonbit_string_t*)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGsE* _block_3131 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_3131)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
    _block_3131->$0 = _M0L6_2atmpS1720;
    _block_3131->$1 = 0;
    return _block_3131;
  } else {
    moonbit_string_t* _M0L6_2atmpS1721 =
      (moonbit_string_t*)moonbit_make_ref_array(_M0L8capacityS505, (moonbit_string_t)moonbit_string_literal_0.data);
    struct _M0TPB5ArrayGsE* _block_3132 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_3132)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
    _block_3132->$0 = _M0L6_2atmpS1721;
    _block_3132->$1 = 0;
    return _block_3132;
  }
}

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(
  int32_t _M0L8capacityS506
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS506 == 0) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS1722 =
      (struct _M0TPC16string10StringView*)moonbit_empty_ref_valtype_array;
    struct _M0TPB5ArrayGRPC16string10StringViewE* _block_3133 =
      (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_block_3133)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 78, 0);
    _block_3133->$0 = _M0L6_2atmpS1722;
    _block_3133->$1 = 0;
    return _block_3133;
  } else {
    struct _M0TPC16string10StringView* _M0L6_2atmpS1723 =
      (struct _M0TPC16string10StringView*)moonbit_make_ref_valtype_array(_M0L8capacityS506, sizeof(struct _M0TPC16string10StringView), Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 81, 0), &(struct _M0TPC16string10StringView){.$0 = (moonbit_string_t)moonbit_string_literal_0.data, .$1 = 0, .$2 = 0});
    struct _M0TPB5ArrayGRPC16string10StringViewE* _block_3134 =
      (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_block_3134)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 78, 0);
    _block_3134->$0 = _M0L6_2atmpS1723;
    _block_3134->$1 = 0;
    return _block_3134;
  }
}

int32_t _M0MPC16string6String11has__prefix(
  moonbit_string_t _M0L4selfS501,
  struct _M0TPC16string10StringView _M0L3strS502
) {
  int32_t _M0L6_2atmpS1715;
  struct _M0TPC16string10StringView _M0L6_2atmpS1714;
  int32_t _result_3135;
  #line 297 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1715 = Moonbit_array_length(_M0L4selfS501);
  moonbit_incref(_M0L4selfS501);
  _M0L6_2atmpS1714
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS501, .$1 = 0, .$2 = _M0L6_2atmpS1715
  };
  #line 299 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3135
  = _M0MPC16string10StringView11has__prefix(_M0L6_2atmpS1714, _M0L3strS502);
  moonbit_decref(_M0L6_2atmpS1714.$0);
  return _result_3135;
}

int32_t _M0MPC16string10StringView11has__prefix(
  struct _M0TPC16string10StringView _M0L4selfS497,
  struct _M0TPC16string10StringView _M0L3strS498
) {
  int64_t _M0L7_2abindS496;
  #line 290 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 292 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS496
  = _M0MPC16string10StringView4find(_M0L4selfS497, _M0L3strS498);
  if (_M0L7_2abindS496 == 4294967296ll) {
    return 0;
  } else {
    int64_t _M0L7_2aSomeS499 = _M0L7_2abindS496;
    int32_t _M0L4_2aiS500 = (int32_t)_M0L7_2aSomeS499;
    return _M0L4_2aiS500 == 0;
  }
}

moonbit_string_t _M0MPC16string6String6repeat(
  moonbit_string_t _M0L4selfS489,
  int32_t _M0L1nS488
) {
  #line 1049 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  if (_M0L1nS488 < 0) {
    #line 1051 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_81.data);
  } else if (_M0L1nS488 == 0) {
    return (moonbit_string_t)moonbit_string_literal_0.data;
  } else if (_M0L1nS488 == 1) {
    moonbit_incref(_M0L4selfS489);
    return _M0L4selfS489;
  } else {
    int32_t _M0L3lenS490 = Moonbit_array_length(_M0L4selfS489);
    int32_t _M0L5totalS491 = _M0L3lenS490 * _M0L1nS488;
    int32_t _if__result_3136;
    if (_M0L3lenS490 == 0) {
      _if__result_3136 = 1;
    } else {
      int32_t _M0L6_2atmpS1712 = _M0L5totalS491 / _M0L1nS488;
      _if__result_3136 = _M0L6_2atmpS1712 == _M0L3lenS490;
    }
    if (_if__result_3136) {
      struct _M0TPB13StringBuilder* _M0L3bufS492;
      moonbit_string_t _M0L3strS493;
      int32_t _M0L2__S494;
      moonbit_string_t _result_3138;
      #line 1060 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3bufS492
      = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L5totalS491);
      #line 1061 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3strS493 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS489);
      _M0L2__S494 = 0;
      while (1) {
        if (_M0L2__S494 < _M0L1nS488) {
          int32_t _M0L6_2atmpS1713;
          #line 1063 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS492, _M0L3strS493);
          _M0L6_2atmpS1713 = _M0L2__S494 + 1;
          _M0L2__S494 = _M0L6_2atmpS1713;
          continue;
        } else {
          moonbit_decref(_M0L3strS493);
        }
        break;
      }
      #line 1065 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _result_3138 = _M0MPB13StringBuilder10to__string(_M0L3bufS492);
      moonbit_decref(_M0L3bufS492);
      return _result_3138;
    } else {
      #line 1058 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_82.data);
    }
  }
}

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(
  moonbit_string_t _M0L4selfS487
) {
  #line 222 "/home/developer/.moon/lib/core/builtin/show.mbt"
  moonbit_incref(_M0L4selfS487);
  return _M0L4selfS487;
}

int64_t _M0MPC16string6String4find(
  moonbit_string_t _M0L4selfS485,
  struct _M0TPC16string10StringView _M0L3strS486
) {
  int32_t _M0L6_2atmpS1711;
  struct _M0TPC16string10StringView _M0L6_2atmpS1710;
  int64_t _result_3139;
  #line 101 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1711 = Moonbit_array_length(_M0L4selfS485);
  moonbit_incref(_M0L4selfS485);
  _M0L6_2atmpS1710
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS485, .$1 = 0, .$2 = _M0L6_2atmpS1711
  };
  #line 102 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3139
  = _M0MPC16string10StringView4find(_M0L6_2atmpS1710, _M0L3strS486);
  moonbit_decref(_M0L6_2atmpS1710.$0);
  return _result_3139;
}

int64_t _M0MPC16string10StringView4find(
  struct _M0TPC16string10StringView _M0L4selfS484,
  struct _M0TPC16string10StringView _M0L3strS483
) {
  int32_t _M0L3endS1708;
  int32_t _M0L5startS1709;
  int32_t _M0L6_2atmpS1707;
  #line 18 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1708 = _M0L3strS483.$2;
  _M0L5startS1709 = _M0L3strS483.$1;
  _M0L6_2atmpS1707 = _M0L3endS1708 - _M0L5startS1709;
  if (_M0L6_2atmpS1707 <= 4) {
    #line 20 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB18brute__force__find(_M0L4selfS484, _M0L3strS483);
  } else {
    #line 22 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB28boyer__moore__horspool__find(_M0L4selfS484, _M0L3strS483);
  }
}

int64_t _M0FPB18brute__force__find(
  struct _M0TPC16string10StringView _M0L8haystackS473,
  struct _M0TPC16string10StringView _M0L6needleS475
) {
  int32_t _M0L3endS1705;
  int32_t _M0L5startS1706;
  int32_t _M0L13haystack__lenS472;
  int32_t _M0L3endS1703;
  int32_t _M0L5startS1704;
  int32_t _M0L11needle__lenS474;
  #line 31 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1705 = _M0L8haystackS473.$2;
  _M0L5startS1706 = _M0L8haystackS473.$1;
  _M0L13haystack__lenS472 = _M0L3endS1705 - _M0L5startS1706;
  _M0L3endS1703 = _M0L6needleS475.$2;
  _M0L5startS1704 = _M0L6needleS475.$1;
  _M0L11needle__lenS474 = _M0L3endS1703 - _M0L5startS1704;
  if (_M0L11needle__lenS474 > 0) {
    if (_M0L13haystack__lenS472 >= _M0L11needle__lenS474) {
      moonbit_string_t _M0L3strS1701 = _M0L6needleS475.$0;
      int32_t _M0L5startS1702 = _M0L6needleS475.$1;
      int32_t _M0L13needle__firstS476 = _M0L3strS1701[_M0L5startS1702];
      int32_t _M0L12forward__lenS477 =
        _M0L13haystack__lenS472 - _M0L11needle__lenS474;
      int32_t _M0L1iS478 = 0;
      while (1) {
        if (_M0L1iS478 <= _M0L12forward__lenS477) {
          moonbit_string_t _M0L3strS1688 = _M0L8haystackS473.$0;
          int32_t _M0L5startS1690 = _M0L8haystackS473.$1;
          int32_t _M0L6_2atmpS1689 = _M0L5startS1690 + _M0L1iS478;
          int32_t _M0L6_2atmpS1687 = _M0L3strS1688[_M0L6_2atmpS1689];
          int32_t _M0L1jS481;
          int32_t _M0L6_2atmpS1686;
          if (_M0L6_2atmpS1687 != _M0L13needle__firstS476) {
            goto join_479;
          }
          _M0L1jS481 = 1;
          while (1) {
            if (_M0L1jS481 < _M0L11needle__lenS474) {
              moonbit_string_t _M0L3strS1696 = _M0L8haystackS473.$0;
              int32_t _M0L5startS1698 = _M0L8haystackS473.$1;
              int32_t _M0L6_2atmpS1699 = _M0L1iS478 + _M0L1jS481;
              int32_t _M0L6_2atmpS1697 = _M0L5startS1698 + _M0L6_2atmpS1699;
              int32_t _M0L6_2atmpS1691 = _M0L3strS1696[_M0L6_2atmpS1697];
              moonbit_string_t _M0L3strS1693 = _M0L6needleS475.$0;
              int32_t _M0L5startS1695 = _M0L6needleS475.$1;
              int32_t _M0L6_2atmpS1694 = _M0L5startS1695 + _M0L1jS481;
              int32_t _M0L6_2atmpS1692 = _M0L3strS1693[_M0L6_2atmpS1694];
              int32_t _M0L6_2atmpS1700;
              if (_M0L6_2atmpS1691 != _M0L6_2atmpS1692) {
                break;
              }
              _M0L6_2atmpS1700 = _M0L1jS481 + 1;
              _M0L1jS481 = _M0L6_2atmpS1700;
              continue;
            } else {
              return (int64_t)_M0L1iS478;
            }
            break;
          }
          goto join_479;
          goto joinlet_3141;
          join_479:;
          _M0L6_2atmpS1686 = _M0L1iS478 + 1;
          _M0L1iS478 = _M0L6_2atmpS1686;
          continue;
          joinlet_3141:;
        }
        break;
      }
      return 4294967296ll;
    } else {
      return 4294967296ll;
    }
  } else {
    return _M0FPB18brute__force__findN6constrS9991;
  }
}

int64_t _M0FPB28boyer__moore__horspool__find(
  struct _M0TPC16string10StringView _M0L8haystackS460,
  struct _M0TPC16string10StringView _M0L6needleS462
) {
  int32_t _M0L3endS1684;
  int32_t _M0L5startS1685;
  int32_t _M0L13haystack__lenS459;
  int32_t _M0L3endS1682;
  int32_t _M0L5startS1683;
  int32_t _M0L11needle__lenS461;
  #line 57 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1684 = _M0L8haystackS460.$2;
  _M0L5startS1685 = _M0L8haystackS460.$1;
  _M0L13haystack__lenS459 = _M0L3endS1684 - _M0L5startS1685;
  _M0L3endS1682 = _M0L6needleS462.$2;
  _M0L5startS1683 = _M0L6needleS462.$1;
  _M0L11needle__lenS461 = _M0L3endS1682 - _M0L5startS1683;
  if (_M0L11needle__lenS461 > 0) {
    if (_M0L13haystack__lenS459 >= _M0L11needle__lenS461) {
      int32_t* _M0L11skip__tableS463 =
        (int32_t*)moonbit_make_int32_array(256, _M0L11needle__lenS461);
      int32_t _M0L7_2abindS464 = _M0L11needle__lenS461 - 1;
      int32_t _M0L1iS465 = 0;
      int32_t _M0L1iS467;
      while (1) {
        if (_M0L1iS465 < _M0L7_2abindS464) {
          moonbit_string_t _M0L3strS1657 = _M0L6needleS462.$0;
          int32_t _M0L5startS1659 = _M0L6needleS462.$1;
          int32_t _M0L6_2atmpS1658 = _M0L5startS1659 + _M0L1iS465;
          int32_t _M0L6_2atmpS1656 = _M0L3strS1657[_M0L6_2atmpS1658];
          int32_t _M0L6_2atmpS1655 = (int32_t)_M0L6_2atmpS1656;
          int32_t _M0L6_2atmpS1652 = _M0L6_2atmpS1655 & 255;
          int32_t _M0L6_2atmpS1654 = _M0L11needle__lenS461 - 1;
          int32_t _M0L6_2atmpS1653 = _M0L6_2atmpS1654 - _M0L1iS465;
          int32_t _M0L6_2atmpS1660;
          if (
            _M0L6_2atmpS1652 < 0
            || _M0L6_2atmpS1652
               >= Moonbit_array_length(_M0L11skip__tableS463)
          ) {
            #line 68 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L11skip__tableS463[_M0L6_2atmpS1652] = _M0L6_2atmpS1653;
          _M0L6_2atmpS1660 = _M0L1iS465 + 1;
          _M0L1iS465 = _M0L6_2atmpS1660;
          continue;
        }
        break;
      }
      _M0L1iS467 = 0;
      while (1) {
        int32_t _M0L6_2atmpS1661 =
          _M0L13haystack__lenS459 - _M0L11needle__lenS461;
        if (_M0L1iS467 <= _M0L6_2atmpS1661) {
          int32_t _M0L7_2abindS468 = _M0L11needle__lenS461 - 1;
          int32_t _M0L1jS469 = 0;
          moonbit_string_t _M0L3strS1677;
          int32_t _M0L5startS1679;
          int32_t _M0L6_2atmpS1681;
          int32_t _M0L6_2atmpS1680;
          int32_t _M0L6_2atmpS1678;
          int32_t _M0L6_2atmpS1676;
          int32_t _M0L6_2atmpS1675;
          int32_t _M0L6_2atmpS1674;
          int32_t _M0L6_2atmpS1673;
          int32_t _M0L6_2atmpS1672;
          while (1) {
            if (_M0L1jS469 <= _M0L7_2abindS468) {
              moonbit_string_t _M0L3strS1667 = _M0L8haystackS460.$0;
              int32_t _M0L5startS1669 = _M0L8haystackS460.$1;
              int32_t _M0L6_2atmpS1670 = _M0L1iS467 + _M0L1jS469;
              int32_t _M0L6_2atmpS1668 = _M0L5startS1669 + _M0L6_2atmpS1670;
              int32_t _M0L6_2atmpS1662 = _M0L3strS1667[_M0L6_2atmpS1668];
              moonbit_string_t _M0L3strS1664 = _M0L6needleS462.$0;
              int32_t _M0L5startS1666 = _M0L6needleS462.$1;
              int32_t _M0L6_2atmpS1665 = _M0L5startS1666 + _M0L1jS469;
              int32_t _M0L6_2atmpS1663 = _M0L3strS1664[_M0L6_2atmpS1665];
              int32_t _M0L6_2atmpS1671;
              if (_M0L6_2atmpS1662 != _M0L6_2atmpS1663) {
                break;
              }
              _M0L6_2atmpS1671 = _M0L1jS469 + 1;
              _M0L1jS469 = _M0L6_2atmpS1671;
              continue;
            } else {
              moonbit_decref(_M0L11skip__tableS463);
              return (int64_t)_M0L1iS467;
            }
            break;
          }
          _M0L3strS1677 = _M0L8haystackS460.$0;
          _M0L5startS1679 = _M0L8haystackS460.$1;
          _M0L6_2atmpS1681 = _M0L1iS467 + _M0L11needle__lenS461;
          _M0L6_2atmpS1680 = _M0L6_2atmpS1681 - 1;
          _M0L6_2atmpS1678 = _M0L5startS1679 + _M0L6_2atmpS1680;
          _M0L6_2atmpS1676 = _M0L3strS1677[_M0L6_2atmpS1678];
          _M0L6_2atmpS1675 = (int32_t)_M0L6_2atmpS1676;
          _M0L6_2atmpS1674 = _M0L6_2atmpS1675 & 255;
          if (
            _M0L6_2atmpS1674 < 0
            || _M0L6_2atmpS1674
               >= Moonbit_array_length(_M0L11skip__tableS463)
          ) {
            #line 73 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS1673 = (int32_t)_M0L11skip__tableS463[_M0L6_2atmpS1674];
          _M0L6_2atmpS1672 = _M0L1iS467 + _M0L6_2atmpS1673;
          _M0L1iS467 = _M0L6_2atmpS1672;
          continue;
        } else {
          moonbit_decref(_M0L11skip__tableS463);
        }
        break;
      }
      return 4294967296ll;
    } else {
      return 4294967296ll;
    }
  } else {
    return _M0FPB28boyer__moore__horspool__findN6constrS9990;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__view(
  struct _M0TPB13StringBuilder* _M0L4selfS458,
  struct _M0TPC16string10StringView _M0L3strS457
) {
  int32_t _M0L3endS1650;
  int32_t _M0L5startS1651;
  int32_t _M0L8str__lenS456;
  int32_t _M0L3lenS1643;
  int32_t _M0L6_2atmpS1642;
  uint16_t* _M0L4dataS1644;
  int32_t _M0L3lenS1645;
  moonbit_string_t _M0L6_2atmpS1646;
  int32_t _M0L6_2atmpS1647;
  int32_t _M0L3lenS1649;
  int32_t _M0L6_2atmpS1648;
  #line 131 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3endS1650 = _M0L3strS457.$2;
  _M0L5startS1651 = _M0L3strS457.$1;
  _M0L8str__lenS456 = _M0L3endS1650 - _M0L5startS1651;
  _M0L3lenS1643 = _M0L4selfS458->$1;
  _M0L6_2atmpS1642 = _M0L3lenS1643 + _M0L8str__lenS456;
  #line 136 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS458, _M0L6_2atmpS1642);
  _M0L4dataS1644 = _M0L4selfS458->$0;
  _M0L3lenS1645 = _M0L4selfS458->$1;
  moonbit_incref(_M0L4dataS1644);
  #line 139 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1646 = _M0MPC16string10StringView4data(_M0L3strS457);
  #line 140 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1647 = _M0MPC16string10StringView13start__offset(_M0L3strS457);
  #line 137 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1644, _M0L3lenS1645, _M0L6_2atmpS1646, _M0L6_2atmpS1647, _M0L8str__lenS456);
  moonbit_decref(_M0L4dataS1644);
  moonbit_decref(_M0L6_2atmpS1646);
  _M0L3lenS1649 = _M0L4selfS458->$1;
  _M0L6_2atmpS1648 = _M0L3lenS1649 + _M0L8str__lenS456;
  _M0L4selfS458->$1 = _M0L6_2atmpS1648;
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t _M0L4selfS454,
  int32_t _M0L13start__offsetS455,
  int64_t _M0L11end__offsetS452
) {
  int32_t _M0L11end__offsetS451;
  int32_t _if__result_3146;
  #line 614 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS452 == 4294967296ll) {
    _M0L11end__offsetS451 = Moonbit_array_length(_M0L4selfS454);
  } else {
    int64_t _M0L7_2aSomeS453 = _M0L11end__offsetS452;
    _M0L11end__offsetS451 = (int32_t)_M0L7_2aSomeS453;
  }
  if (_M0L13start__offsetS455 >= 0) {
    if (_M0L13start__offsetS455 <= _M0L11end__offsetS451) {
      int32_t _M0L6_2atmpS1641 = Moonbit_array_length(_M0L4selfS454);
      _if__result_3146 = _M0L11end__offsetS451 <= _M0L6_2atmpS1641;
    } else {
      _if__result_3146 = 0;
    }
  } else {
    _if__result_3146 = 0;
  }
  if (_if__result_3146) {
    moonbit_incref(_M0L4selfS454);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS454,
                                                 .$1 = _M0L13start__offsetS455,
                                                 .$2 = _M0L11end__offsetS451};
  } else {
    #line 623 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_83.data);
  }
}

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView _M0L4selfS446
) {
  int32_t _M0L5startS445;
  int32_t _M0L3endS447;
  struct _M0TPB8MutLocalGiE* _M0L5indexS448;
  struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__* _closure_3147;
  struct _M0TWEOc* _M0L6_2atmpS1620;
  struct _M0TPB4IterGcE* _result_3148;
  #line 214 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L5startS445 = _M0L4selfS446.$1;
  _M0L3endS447 = _M0L4selfS446.$2;
  _M0L5indexS448
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS448)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS448->$0 = _M0L5startS445;
  moonbit_incref(_M0L4selfS446.$0);
  _closure_3147
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__*)moonbit_malloc(sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__));
  Moonbit_object_header(_closure_3147)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _closure_3147->code = &_M0MPC16string10StringView4iterC1621l219;
  _closure_3147->$0 = _M0L5indexS448;
  _closure_3147->$1 = _M0L3endS447;
  _closure_3147->$2 = _M0L4selfS446;
  _M0L6_2atmpS1620 = (struct _M0TWEOc*)_closure_3147;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_3148 = _M0MPB4Iter3newGcE(_M0L6_2atmpS1620, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1620);
  return _result_3148;
}

int32_t _M0MPC16string10StringView4iterC1621l219(
  struct _M0TWEOc* _M0L6_2aenvS1622
) {
  struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__* _M0L14_2acasted__envS1623;
  struct _M0TPC16string10StringView _M0L4selfS446;
  int32_t _M0L3endS447;
  struct _M0TPB8MutLocalGiE* _M0L5indexS448;
  int32_t _M0L3valS1624;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L14_2acasted__envS1623
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1621__l219__*)_M0L6_2aenvS1622;
  _M0L4selfS446 = _M0L14_2acasted__envS1623->$2;
  _M0L3endS447 = _M0L14_2acasted__envS1623->$1;
  _M0L5indexS448 = _M0L14_2acasted__envS1623->$0;
  _M0L3valS1624 = _M0L5indexS448->$0;
  if (_M0L3valS1624 < _M0L3endS447) {
    moonbit_string_t _M0L3strS1639 = _M0L4selfS446.$0;
    int32_t _M0L3valS1640 = _M0L5indexS448->$0;
    int32_t _M0L2c1S449 = _M0L3strS1639[_M0L3valS1640];
    int32_t _if__result_3149;
    int32_t _M0L3valS1637;
    int32_t _M0L6_2atmpS1636;
    int32_t _M0L6_2atmpS1638;
    #line 222 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S449)) {
      int32_t _M0L3valS1627 = _M0L5indexS448->$0;
      int32_t _M0L6_2atmpS1625 = _M0L3valS1627 + 1;
      int32_t _M0L3endS1626 = _M0L4selfS446.$2;
      _if__result_3149 = _M0L6_2atmpS1625 < _M0L3endS1626;
    } else {
      _if__result_3149 = 0;
    }
    if (_if__result_3149) {
      moonbit_string_t _M0L3strS1633 = _M0L4selfS446.$0;
      int32_t _M0L3valS1635 = _M0L5indexS448->$0;
      int32_t _M0L6_2atmpS1634 = _M0L3valS1635 + 1;
      int32_t _M0L2c2S450 = _M0L3strS1633[_M0L6_2atmpS1634];
      #line 224 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S450)) {
        int32_t _M0L3valS1629 = _M0L5indexS448->$0;
        int32_t _M0L6_2atmpS1628 = _M0L3valS1629 + 2;
        int32_t _M0L6_2atmpS1631;
        int32_t _M0L6_2atmpS1632;
        int32_t _M0L6_2atmpS1630;
        _M0L5indexS448->$0 = _M0L6_2atmpS1628;
        _M0L6_2atmpS1631 = (int32_t)_M0L2c1S449;
        _M0L6_2atmpS1632 = (int32_t)_M0L2c2S450;
        #line 226 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        _M0L6_2atmpS1630
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1631, _M0L6_2atmpS1632);
        return _M0L6_2atmpS1630;
      }
    }
    _M0L3valS1637 = _M0L5indexS448->$0;
    _M0L6_2atmpS1636 = _M0L3valS1637 + 1;
    _M0L5indexS448->$0 = _M0L6_2atmpS1636;
    #line 230 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    _M0L6_2atmpS1638 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S449);
    return _M0L6_2atmpS1638;
  } else {
    return -1;
  }
}

moonbit_string_t _M0MPC16string10StringView9to__owned(
  struct _M0TPC16string10StringView _M0L4selfS444
) {
  moonbit_string_t _M0L3strS1617;
  int32_t _M0L5startS1618;
  int32_t _M0L3endS1619;
  moonbit_string_t _result_3150;
  #line 196 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1617 = _M0L4selfS444.$0;
  _M0L5startS1618 = _M0L4selfS444.$1;
  _M0L3endS1619 = _M0L4selfS444.$2;
  moonbit_incref(_M0L3strS1617);
  #line 199 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_3150
  = _M0MPC16string6String17unsafe__substring(_M0L3strS1617, _M0L5startS1618, _M0L3endS1619);
  moonbit_decref(_M0L3strS1617);
  return _result_3150;
}

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t _M0L3strS441,
  int32_t _M0L5startS439,
  int32_t _M0L3endS440
) {
  int32_t _if__result_3151;
  int32_t _M0L3lenS442;
  int32_t _M0L6_2atmpS1615;
  int32_t _M0L6_2atmpS1616;
  moonbit_bytes_t _M0L5bytesS443;
  moonbit_bytes_t _M0L6_2atmpS1614;
  moonbit_string_t _result_3152;
  #line 91 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L5startS439 == 0) {
    int32_t _M0L6_2atmpS1613 = Moonbit_array_length(_M0L3strS441);
    _if__result_3151 = _M0L3endS440 == _M0L6_2atmpS1613;
  } else {
    _if__result_3151 = 0;
  }
  if (_if__result_3151) {
    moonbit_incref(_M0L3strS441);
    return _M0L3strS441;
  }
  _M0L3lenS442 = _M0L3endS440 - _M0L5startS439;
  _M0L6_2atmpS1615 = _M0L3lenS442 * 2;
  #line 101 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1616 = _M0IPC14byte4BytePB7Default7default();
  _M0L5bytesS443
  = (moonbit_bytes_t)moonbit_make_bytes(_M0L6_2atmpS1615, _M0L6_2atmpS1616);
  #line 102 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0MPC15array10FixedArray18blit__from__string(_M0L5bytesS443, 0, _M0L3strS441, _M0L5startS439, _M0L3lenS442);
  _M0L6_2atmpS1614 = _M0L5bytesS443;
  #line 103 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_3152
  = _M0MPC15bytes5Bytes29to__unchecked__string_2einner(_M0L6_2atmpS1614, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1614);
  return _result_3152;
}

int32_t _M0IPC14byte4BytePB7Default7default() {
  #line 231 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  return 0;
}

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t _M0L4selfS434,
  int32_t _M0L6offsetS438,
  int64_t _M0L6lengthS436
) {
  int32_t _M0L3lenS433;
  int32_t _M0L6lengthS435;
  int32_t _if__result_3153;
  #line 77 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L3lenS433 = Moonbit_array_length(_M0L4selfS434);
  if (_M0L6lengthS436 == 4294967296ll) {
    _M0L6lengthS435 = _M0L3lenS433 - _M0L6offsetS438;
  } else {
    int64_t _M0L7_2aSomeS437 = _M0L6lengthS436;
    _M0L6lengthS435 = (int32_t)_M0L7_2aSomeS437;
  }
  if (_M0L6offsetS438 >= 0) {
    if (_M0L6lengthS435 >= 0) {
      int32_t _M0L6_2atmpS1612 = _M0L6offsetS438 + _M0L6lengthS435;
      _if__result_3153 = _M0L6_2atmpS1612 <= _M0L3lenS433;
    } else {
      _if__result_3153 = 0;
    }
  } else {
    _if__result_3153 = 0;
  }
  if (_if__result_3153) {
    moonbit_incref(_M0L4selfS434);
    #line 85 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    return _M0FPB19unsafe__sub__string(_M0L4selfS434, _M0L6offsetS438, _M0L6lengthS435);
  } else {
    #line 84 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t _M0L4selfS425,
  int32_t _M0L13bytes__offsetS420,
  moonbit_string_t _M0L3strS427,
  int32_t _M0L11str__offsetS423,
  int32_t _M0L6lengthS421
) {
  int32_t _M0L6_2atmpS1611;
  int32_t _M0L6_2atmpS1610;
  int32_t _M0L2e1S419;
  int32_t _M0L6_2atmpS1609;
  int32_t _M0L2e2S422;
  int32_t _M0L4len1S424;
  int32_t _M0L4len2S426;
  int32_t _if__result_3154;
  #line 125 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L6_2atmpS1611 = _M0L6lengthS421 * 2;
  _M0L6_2atmpS1610 = _M0L13bytes__offsetS420 + _M0L6_2atmpS1611;
  _M0L2e1S419 = _M0L6_2atmpS1610 - 1;
  _M0L6_2atmpS1609 = _M0L11str__offsetS423 + _M0L6lengthS421;
  _M0L2e2S422 = _M0L6_2atmpS1609 - 1;
  _M0L4len1S424 = Moonbit_array_length(_M0L4selfS425);
  _M0L4len2S426 = Moonbit_array_length(_M0L3strS427);
  if (_M0L6lengthS421 >= 0) {
    if (_M0L13bytes__offsetS420 >= 0) {
      if (_M0L2e1S419 < _M0L4len1S424) {
        if (_M0L11str__offsetS423 >= 0) {
          _if__result_3154 = _M0L2e2S422 < _M0L4len2S426;
        } else {
          _if__result_3154 = 0;
        }
      } else {
        _if__result_3154 = 0;
      }
    } else {
      _if__result_3154 = 0;
    }
  } else {
    _if__result_3154 = 0;
  }
  if (_if__result_3154) {
    int32_t _M0L16end__str__offsetS428 =
      _M0L11str__offsetS423 + _M0L6lengthS421;
    int32_t _M0L1iS429 = _M0L11str__offsetS423;
    int32_t _M0L1jS430 = _M0L13bytes__offsetS420;
    while (1) {
      if (_M0L1iS429 < _M0L16end__str__offsetS428) {
        int32_t _M0L6_2atmpS1606 = _M0L3strS427[_M0L1iS429];
        int32_t _M0L6_2atmpS1605 = (int32_t)_M0L6_2atmpS1606;
        uint32_t _M0L1cS431 = *(uint32_t*)&_M0L6_2atmpS1605;
        uint32_t _M0L6_2atmpS1601 = _M0L1cS431 & 255u;
        int32_t _M0L6_2atmpS1600;
        int32_t _M0L6_2atmpS1602;
        uint32_t _M0L6_2atmpS1604;
        int32_t _M0L6_2atmpS1603;
        int32_t _M0L6_2atmpS1607;
        int32_t _M0L6_2atmpS1608;
        #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1600 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1601);
        if (
          _M0L1jS430 < 0 || _M0L1jS430 >= Moonbit_array_length(_M0L4selfS425)
        ) {
          #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS425[_M0L1jS430] = _M0L6_2atmpS1600;
        _M0L6_2atmpS1602 = _M0L1jS430 + 1;
        _M0L6_2atmpS1604 = _M0L1cS431 >> 8;
        #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1603 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1604);
        if (
          _M0L6_2atmpS1602 < 0
          || _M0L6_2atmpS1602 >= Moonbit_array_length(_M0L4selfS425)
        ) {
          #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS425[_M0L6_2atmpS1602] = _M0L6_2atmpS1603;
        _M0L6_2atmpS1607 = _M0L1iS429 + 1;
        _M0L6_2atmpS1608 = _M0L1jS430 + 2;
        _M0L1iS429 = _M0L6_2atmpS1607;
        _M0L1jS430 = _M0L6_2atmpS1608;
        continue;
      }
      break;
    }
  } else {
    #line 138 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
  return 0;
}

int32_t _M0MPC14uint4UInt8to__byte(uint32_t _M0L4selfS418) {
  int32_t _M0L6_2atmpS1599;
  #line 2519 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1599 = *(int32_t*)&_M0L4selfS418;
  return _M0L6_2atmpS1599 & 0xff;
}

int32_t* _M0MPC15array5Array6bufferGiE(struct _M0TPB5ArrayGiE* _M0L4selfS413) {
  int32_t* _M0L8_2afieldS2848;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS2848 = _M0L4selfS413->$0;
  moonbit_incref(_M0L8_2afieldS2848);
  return _M0L8_2afieldS2848;
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS414
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8_2afieldS2849;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS2849 = _M0L4selfS414->$0;
  moonbit_incref(_M0L8_2afieldS2849);
  return _M0L8_2afieldS2849;
}

moonbit_string_t* _M0MPC15array5Array6bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS415
) {
  moonbit_string_t* _M0L8_2afieldS2850;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS2850 = _M0L4selfS415->$0;
  moonbit_incref(_M0L8_2afieldS2850);
  return _M0L8_2afieldS2850;
}

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS416
) {
  struct _M0TUsiE** _M0L8_2afieldS2851;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS2851 = _M0L4selfS416->$0;
  moonbit_incref(_M0L8_2afieldS2851);
  return _M0L8_2afieldS2851;
}

struct _M0TPC16string10StringView* _M0MPC15array5Array6bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS417
) {
  struct _M0TPC16string10StringView* _M0L8_2afieldS2852;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS2852 = _M0L4selfS417->$0;
  moonbit_incref(_M0L8_2afieldS2852);
  return _M0L8_2afieldS2852;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView _M0L4selfS411,
  int32_t _M0L13start__offsetS412,
  int64_t _M0L11end__offsetS409
) {
  int32_t _M0L11end__offsetS408;
  int32_t _if__result_3156;
  #line 105 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS409 == 4294967296ll) {
    int32_t _M0L3endS1597 = _M0L4selfS411.$2;
    int32_t _M0L5startS1598 = _M0L4selfS411.$1;
    _M0L11end__offsetS408 = _M0L3endS1597 - _M0L5startS1598;
  } else {
    int64_t _M0L7_2aSomeS410 = _M0L11end__offsetS409;
    _M0L11end__offsetS408 = (int32_t)_M0L7_2aSomeS410;
  }
  if (_M0L13start__offsetS412 >= 0) {
    if (_M0L13start__offsetS412 <= _M0L11end__offsetS408) {
      int32_t _M0L3endS1590 = _M0L4selfS411.$2;
      int32_t _M0L5startS1591 = _M0L4selfS411.$1;
      int32_t _M0L6_2atmpS1589 = _M0L3endS1590 - _M0L5startS1591;
      _if__result_3156 = _M0L11end__offsetS408 <= _M0L6_2atmpS1589;
    } else {
      _if__result_3156 = 0;
    }
  } else {
    _if__result_3156 = 0;
  }
  if (_if__result_3156) {
    moonbit_string_t _M0L3strS1592 = _M0L4selfS411.$0;
    int32_t _M0L5startS1596 = _M0L4selfS411.$1;
    int32_t _M0L6_2atmpS1593 = _M0L5startS1596 + _M0L13start__offsetS412;
    int32_t _M0L5startS1595 = _M0L4selfS411.$1;
    int32_t _M0L6_2atmpS1594 = _M0L5startS1595 + _M0L11end__offsetS408;
    moonbit_incref(_M0L3strS1592);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1592,
                                                 .$1 = _M0L6_2atmpS1593,
                                                 .$2 = _M0L6_2atmpS1594};
  } else {
    #line 114 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_83.data);
  }
}

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(
  struct _M0TWEOc* _M0L1fS397,
  int64_t _M0L10size__hintS394
) {
  int64_t _M0L10size__hintS393;
  struct _M0TPB4IterGcE* _block_3157;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS394 == 4294967296ll) {
    _M0L10size__hintS393 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS395 = _M0L10size__hintS394;
    int32_t _M0L4_2anS396 = (int32_t)_M0L7_2aSomeS395;
    if (_M0L4_2anS396 > 0) {
      _M0L10size__hintS393 = (int64_t)_M0L4_2anS396;
    } else {
      _M0L10size__hintS393 = _M0MPB4Iter3newN6constrS9988GcE;
    }
  }
  moonbit_incref(_M0L1fS397);
  _block_3157
  = (struct _M0TPB4IterGcE*)moonbit_malloc(sizeof(struct _M0TPB4IterGcE));
  Moonbit_object_header(_block_3157)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _block_3157->$0 = _M0L1fS397;
  _block_3157->$1 = _M0L10size__hintS393;
  return _block_3157;
}

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB4Iter3newGUsRPB5ArrayGiEEE(
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L1fS402,
  int64_t _M0L10size__hintS399
) {
  int64_t _M0L10size__hintS398;
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _block_3158;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS399 == 4294967296ll) {
    _M0L10size__hintS398 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS400 = _M0L10size__hintS399;
    int32_t _M0L4_2anS401 = (int32_t)_M0L7_2aSomeS400;
    if (_M0L4_2anS401 > 0) {
      _M0L10size__hintS398 = (int64_t)_M0L4_2anS401;
    } else {
      _M0L10size__hintS398 = _M0MPB4Iter3newN6constrS9988GUsRPB5ArrayGiEEE;
    }
  }
  moonbit_incref(_M0L1fS402);
  _block_3158
  = (struct _M0TPB4IterGUsRPB5ArrayGiEEE*)moonbit_malloc(sizeof(struct _M0TPB4IterGUsRPB5ArrayGiEEE));
  Moonbit_object_header(_block_3158)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _block_3158->$0 = _M0L1fS402;
  _block_3158->$1 = _M0L10size__hintS398;
  return _block_3158;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L1fS407,
  int64_t _M0L10size__hintS404
) {
  int64_t _M0L10size__hintS403;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_3159;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS404 == 4294967296ll) {
    _M0L10size__hintS403 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS405 = _M0L10size__hintS404;
    int32_t _M0L4_2anS406 = (int32_t)_M0L7_2aSomeS405;
    if (_M0L4_2anS406 > 0) {
      _M0L10size__hintS403 = (int64_t)_M0L4_2anS406;
    } else {
      _M0L10size__hintS403
      = _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE;
    }
  }
  moonbit_incref(_M0L1fS407);
  _block_3159
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_3159)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 98, 0);
  _block_3159->$0 = _M0L1fS407;
  _block_3159->$1 = _M0L10size__hintS403;
  return _block_3159;
}

moonbit_string_t _M0MPC13int3Int18to__string_2einner(
  int32_t _M0L4selfS377,
  int32_t _M0L5radixS376
) {
  int32_t _if__result_3160;
  int32_t _M0L12is__negativeS378;
  uint32_t _M0L3numS379;
  uint16_t* _M0L6bufferS380;
  #line 209 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS376 < 2) {
    _if__result_3160 = 1;
  } else {
    _if__result_3160 = _M0L5radixS376 > 36;
  }
  if (_if__result_3160) {
    #line 213 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_84.data);
  }
  if (_M0L4selfS377 == 0) {
    return (moonbit_string_t)moonbit_string_literal_85.data;
  }
  _M0L12is__negativeS378 = _M0L4selfS377 < 0;
  if (_M0L12is__negativeS378) {
    int32_t _M0L6_2atmpS1588 = -_M0L4selfS377;
    _M0L3numS379 = *(uint32_t*)&_M0L6_2atmpS1588;
  } else {
    _M0L3numS379 = *(uint32_t*)&_M0L4selfS377;
  }
  switch (_M0L5radixS376) {
    case 10: {
      int32_t _M0L10digit__lenS381;
      int32_t _M0L6_2atmpS1585;
      int32_t _M0L10total__lenS382;
      uint16_t* _M0L6bufferS383;
      int32_t _M0L12digit__startS384;
      #line 235 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS381 = _M0FPB12dec__count32(_M0L3numS379);
      if (_M0L12is__negativeS378) {
        _M0L6_2atmpS1585 = 1;
      } else {
        _M0L6_2atmpS1585 = 0;
      }
      _M0L10total__lenS382 = _M0L10digit__lenS381 + _M0L6_2atmpS1585;
      _M0L6bufferS383
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS382, 0);
      if (_M0L12is__negativeS378) {
        _M0L12digit__startS384 = 1;
      } else {
        _M0L12digit__startS384 = 0;
      }
      #line 239 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__dec(_M0L6bufferS383, _M0L3numS379, _M0L12digit__startS384, _M0L10total__lenS382);
      _M0L6bufferS380 = _M0L6bufferS383;
      break;
    }
    
    case 16: {
      int32_t _M0L10digit__lenS385;
      int32_t _M0L6_2atmpS1586;
      int32_t _M0L10total__lenS386;
      uint16_t* _M0L6bufferS387;
      int32_t _M0L12digit__startS388;
      #line 243 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS385 = _M0FPB12hex__count32(_M0L3numS379);
      if (_M0L12is__negativeS378) {
        _M0L6_2atmpS1586 = 1;
      } else {
        _M0L6_2atmpS1586 = 0;
      }
      _M0L10total__lenS386 = _M0L10digit__lenS385 + _M0L6_2atmpS1586;
      _M0L6bufferS387
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS386, 0);
      if (_M0L12is__negativeS378) {
        _M0L12digit__startS388 = 1;
      } else {
        _M0L12digit__startS388 = 0;
      }
      #line 247 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__hex(_M0L6bufferS387, _M0L3numS379, _M0L12digit__startS388, _M0L10total__lenS386);
      _M0L6bufferS380 = _M0L6bufferS387;
      break;
    }
    default: {
      int32_t _M0L10digit__lenS389;
      int32_t _M0L6_2atmpS1587;
      int32_t _M0L10total__lenS390;
      uint16_t* _M0L6bufferS391;
      int32_t _M0L12digit__startS392;
      #line 251 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS389
      = _M0FPB14radix__count32(_M0L3numS379, _M0L5radixS376);
      if (_M0L12is__negativeS378) {
        _M0L6_2atmpS1587 = 1;
      } else {
        _M0L6_2atmpS1587 = 0;
      }
      _M0L10total__lenS390 = _M0L10digit__lenS389 + _M0L6_2atmpS1587;
      _M0L6bufferS391
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS390, 0);
      if (_M0L12is__negativeS378) {
        _M0L12digit__startS392 = 1;
      } else {
        _M0L12digit__startS392 = 0;
      }
      #line 255 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB24int__to__string__generic(_M0L6bufferS391, _M0L3numS379, _M0L12digit__startS392, _M0L10total__lenS390, _M0L5radixS376);
      _M0L6bufferS380 = _M0L6bufferS391;
      break;
    }
  }
  if (_M0L12is__negativeS378) {
    _M0L6bufferS380[0] = 45;
  }
  return _M0L6bufferS380;
}

int32_t _M0FPB14radix__count32(
  uint32_t _M0L5valueS370,
  int32_t _M0L5radixS372
) {
  uint32_t _M0L4baseS371;
  uint32_t _M0L3numS373;
  int32_t _M0L5countS374;
  #line 189 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS370 == 0u) {
    return 1;
  }
  _M0L4baseS371 = *(uint32_t*)&_M0L5radixS372;
  _M0L3numS373 = _M0L5valueS370;
  _M0L5countS374 = 0;
  while (1) {
    if (_M0L3numS373 > 0u) {
      uint32_t _M0L6_2atmpS1583 = _M0L3numS373 / _M0L4baseS371;
      int32_t _M0L6_2atmpS1584 = _M0L5countS374 + 1;
      _M0L3numS373 = _M0L6_2atmpS1583;
      _M0L5countS374 = _M0L6_2atmpS1584;
      continue;
    } else {
      return _M0L5countS374;
    }
    break;
  }
}

int32_t _M0FPB12hex__count32(uint32_t _M0L5valueS368) {
  #line 177 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS368 == 0u) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS369;
    int32_t _M0L6_2atmpS1582;
    int32_t _M0L6_2atmpS1581;
    #line 182 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS369 = moonbit_clz32(_M0L5valueS368);
    _M0L6_2atmpS1582 = 31 - _M0L14leading__zerosS369;
    _M0L6_2atmpS1581 = _M0L6_2atmpS1582 / 4;
    return _M0L6_2atmpS1581 + 1;
  }
}

int32_t _M0FPB12dec__count32(uint32_t _M0L5valueS367) {
  #line 143 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS367 >= 100000u) {
    if (_M0L5valueS367 >= 10000000u) {
      if (_M0L5valueS367 >= 1000000000u) {
        return 10;
      } else if (_M0L5valueS367 >= 100000000u) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS367 >= 1000000u) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS367 >= 1000u) {
    if (_M0L5valueS367 >= 10000u) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS367 >= 100u) {
    return 3;
  } else if (_M0L5valueS367 >= 10u) {
    return 2;
  } else {
    return 1;
  }
}

int32_t _M0FPB20int__to__string__dec(
  uint16_t* _M0L6bufferS353,
  uint32_t _M0L3numS365,
  int32_t _M0L12digit__startS354,
  int32_t _M0L10total__lenS366
) {
  int32_t _M0L6_2atmpS1580;
  uint32_t _M0L3numS343;
  int32_t _M0L6offsetS344;
  #line 88 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1580 = _M0L10total__lenS366 - _M0L12digit__startS354;
  _M0L3numS343 = _M0L3numS365;
  _M0L6offsetS344 = _M0L6_2atmpS1580;
  while (1) {
    if (_M0L3numS343 >= 10000u) {
      uint32_t _M0L1tS345 = _M0L3numS343 / 10000u;
      uint32_t _M0L6_2atmpS1557 = _M0L3numS343 % 10000u;
      int32_t _M0L1rS346 = *(int32_t*)&_M0L6_2atmpS1557;
      int32_t _M0L2d1S347 = _M0L1rS346 / 100;
      int32_t _M0L2d2S348 = _M0L1rS346 % 100;
      int32_t _M0L6_2atmpS1556 = _M0L2d1S347 / 10;
      int32_t _M0L6_2atmpS1555 = 48 + _M0L6_2atmpS1556;
      int32_t _M0L6d1__hiS349 = (uint16_t)_M0L6_2atmpS1555;
      int32_t _M0L6_2atmpS1554 = _M0L2d1S347 % 10;
      int32_t _M0L6_2atmpS1553 = 48 + _M0L6_2atmpS1554;
      int32_t _M0L6d1__loS350 = (uint16_t)_M0L6_2atmpS1553;
      int32_t _M0L6_2atmpS1552 = _M0L2d2S348 / 10;
      int32_t _M0L6_2atmpS1551 = 48 + _M0L6_2atmpS1552;
      int32_t _M0L6d2__hiS351 = (uint16_t)_M0L6_2atmpS1551;
      int32_t _M0L6_2atmpS1550 = _M0L2d2S348 % 10;
      int32_t _M0L6_2atmpS1549 = 48 + _M0L6_2atmpS1550;
      int32_t _M0L6d2__loS352 = (uint16_t)_M0L6_2atmpS1549;
      int32_t _M0L6_2atmpS1541 = _M0L12digit__startS354 + _M0L6offsetS344;
      int32_t _M0L6_2atmpS1540 = _M0L6_2atmpS1541 - 4;
      int32_t _M0L6_2atmpS1543;
      int32_t _M0L6_2atmpS1542;
      int32_t _M0L6_2atmpS1545;
      int32_t _M0L6_2atmpS1544;
      int32_t _M0L6_2atmpS1547;
      int32_t _M0L6_2atmpS1546;
      int32_t _M0L6_2atmpS1548;
      _M0L6bufferS353[_M0L6_2atmpS1540] = _M0L6d1__hiS349;
      _M0L6_2atmpS1543 = _M0L12digit__startS354 + _M0L6offsetS344;
      _M0L6_2atmpS1542 = _M0L6_2atmpS1543 - 3;
      _M0L6bufferS353[_M0L6_2atmpS1542] = _M0L6d1__loS350;
      _M0L6_2atmpS1545 = _M0L12digit__startS354 + _M0L6offsetS344;
      _M0L6_2atmpS1544 = _M0L6_2atmpS1545 - 2;
      _M0L6bufferS353[_M0L6_2atmpS1544] = _M0L6d2__hiS351;
      _M0L6_2atmpS1547 = _M0L12digit__startS354 + _M0L6offsetS344;
      _M0L6_2atmpS1546 = _M0L6_2atmpS1547 - 1;
      _M0L6bufferS353[_M0L6_2atmpS1546] = _M0L6d2__loS352;
      _M0L6_2atmpS1548 = _M0L6offsetS344 - 4;
      _M0L3numS343 = _M0L1tS345;
      _M0L6offsetS344 = _M0L6_2atmpS1548;
      continue;
    } else {
      int32_t _M0L6_2atmpS1579 = *(int32_t*)&_M0L3numS343;
      int32_t _M0L9remainingS356 = _M0L6_2atmpS1579;
      int32_t _M0L6offsetS357 = _M0L6offsetS344;
      while (1) {
        if (_M0L9remainingS356 >= 100) {
          int32_t _M0L1tS358 = _M0L9remainingS356 / 100;
          int32_t _M0L1dS359 = _M0L9remainingS356 % 100;
          int32_t _M0L6_2atmpS1566 = _M0L1dS359 / 10;
          int32_t _M0L6_2atmpS1565 = 48 + _M0L6_2atmpS1566;
          int32_t _M0L5d__hiS360 = (uint16_t)_M0L6_2atmpS1565;
          int32_t _M0L6_2atmpS1564 = _M0L1dS359 % 10;
          int32_t _M0L6_2atmpS1563 = 48 + _M0L6_2atmpS1564;
          int32_t _M0L5d__loS361 = (uint16_t)_M0L6_2atmpS1563;
          int32_t _M0L6_2atmpS1559 = _M0L12digit__startS354 + _M0L6offsetS357;
          int32_t _M0L6_2atmpS1558 = _M0L6_2atmpS1559 - 2;
          int32_t _M0L6_2atmpS1561;
          int32_t _M0L6_2atmpS1560;
          int32_t _M0L6_2atmpS1562;
          _M0L6bufferS353[_M0L6_2atmpS1558] = _M0L5d__hiS360;
          _M0L6_2atmpS1561 = _M0L12digit__startS354 + _M0L6offsetS357;
          _M0L6_2atmpS1560 = _M0L6_2atmpS1561 - 1;
          _M0L6bufferS353[_M0L6_2atmpS1560] = _M0L5d__loS361;
          _M0L6_2atmpS1562 = _M0L6offsetS357 - 2;
          _M0L9remainingS356 = _M0L1tS358;
          _M0L6offsetS357 = _M0L6_2atmpS1562;
          continue;
        } else if (_M0L9remainingS356 >= 10) {
          int32_t _M0L6_2atmpS1574 = _M0L9remainingS356 / 10;
          int32_t _M0L6_2atmpS1573 = 48 + _M0L6_2atmpS1574;
          int32_t _M0L5d__hiS363 = (uint16_t)_M0L6_2atmpS1573;
          int32_t _M0L6_2atmpS1572 = _M0L9remainingS356 % 10;
          int32_t _M0L6_2atmpS1571 = 48 + _M0L6_2atmpS1572;
          int32_t _M0L5d__loS364 = (uint16_t)_M0L6_2atmpS1571;
          int32_t _M0L6_2atmpS1568 = _M0L12digit__startS354 + _M0L6offsetS357;
          int32_t _M0L6_2atmpS1567 = _M0L6_2atmpS1568 - 2;
          int32_t _M0L6_2atmpS1570;
          int32_t _M0L6_2atmpS1569;
          _M0L6bufferS353[_M0L6_2atmpS1567] = _M0L5d__hiS363;
          _M0L6_2atmpS1570 = _M0L12digit__startS354 + _M0L6offsetS357;
          _M0L6_2atmpS1569 = _M0L6_2atmpS1570 - 1;
          _M0L6bufferS353[_M0L6_2atmpS1569] = _M0L5d__loS364;
        } else {
          int32_t _M0L6_2atmpS1578 = _M0L12digit__startS354 + _M0L6offsetS357;
          int32_t _M0L6_2atmpS1575 = _M0L6_2atmpS1578 - 1;
          int32_t _M0L6_2atmpS1577 = 48 + _M0L9remainingS356;
          int32_t _M0L6_2atmpS1576 = (uint16_t)_M0L6_2atmpS1577;
          _M0L6bufferS353[_M0L6_2atmpS1575] = _M0L6_2atmpS1576;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB24int__to__string__generic(
  uint16_t* _M0L6bufferS333,
  uint32_t _M0L3numS337,
  int32_t _M0L12digit__startS334,
  int32_t _M0L10total__lenS336,
  int32_t _M0L5radixS327
) {
  uint32_t _M0L4baseS326;
  int32_t _M0L6_2atmpS1525;
  int32_t _M0L6_2atmpS1524;
  #line 57 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS326 = *(uint32_t*)&_M0L5radixS327;
  _M0L6_2atmpS1525 = _M0L5radixS327 - 1;
  _M0L6_2atmpS1524 = _M0L5radixS327 & _M0L6_2atmpS1525;
  if (_M0L6_2atmpS1524 == 0) {
    int32_t _M0L5shiftS328;
    uint32_t _M0L4maskS329;
    int32_t _M0L6_2atmpS1532;
    int32_t _M0L6offsetS330;
    uint32_t _M0L1nS331;
    #line 68 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS328 = moonbit_ctz32(_M0L5radixS327);
    _M0L4maskS329 = _M0L4baseS326 - 1u;
    _M0L6_2atmpS1532 = _M0L10total__lenS336 - _M0L12digit__startS334;
    _M0L6offsetS330 = _M0L6_2atmpS1532;
    _M0L1nS331 = _M0L3numS337;
    while (1) {
      if (_M0L1nS331 > 0u) {
        uint32_t _M0L6_2atmpS1531 = _M0L1nS331 & _M0L4maskS329;
        int32_t _M0L5digitS332 = *(int32_t*)&_M0L6_2atmpS1531;
        int32_t _M0L6_2atmpS1528 = _M0L12digit__startS334 + _M0L6offsetS330;
        int32_t _M0L6_2atmpS1526 = _M0L6_2atmpS1528 - 1;
        int32_t _M0L6_2atmpS1527 =
          ((moonbit_string_t)moonbit_string_literal_86.data)[_M0L5digitS332];
        int32_t _M0L6_2atmpS1529;
        uint32_t _M0L6_2atmpS1530;
        _M0L6bufferS333[_M0L6_2atmpS1526] = _M0L6_2atmpS1527;
        _M0L6_2atmpS1529 = _M0L6offsetS330 - 1;
        _M0L6_2atmpS1530 = _M0L1nS331 >> (_M0L5shiftS328 & 31);
        _M0L6offsetS330 = _M0L6_2atmpS1529;
        _M0L1nS331 = _M0L6_2atmpS1530;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1539 = _M0L10total__lenS336 - _M0L12digit__startS334;
    int32_t _M0L6offsetS338 = _M0L6_2atmpS1539;
    uint32_t _M0L1nS339 = _M0L3numS337;
    while (1) {
      if (_M0L1nS339 > 0u) {
        uint32_t _M0L1qS340 = _M0L1nS339 / _M0L4baseS326;
        uint32_t _M0L6_2atmpS1538 = _M0L1qS340 * _M0L4baseS326;
        uint32_t _M0L6_2atmpS1537 = _M0L1nS339 - _M0L6_2atmpS1538;
        int32_t _M0L5digitS341 = *(int32_t*)&_M0L6_2atmpS1537;
        int32_t _M0L6_2atmpS1535 = _M0L12digit__startS334 + _M0L6offsetS338;
        int32_t _M0L6_2atmpS1533 = _M0L6_2atmpS1535 - 1;
        int32_t _M0L6_2atmpS1534 =
          ((moonbit_string_t)moonbit_string_literal_86.data)[_M0L5digitS341];
        int32_t _M0L6_2atmpS1536;
        _M0L6bufferS333[_M0L6_2atmpS1533] = _M0L6_2atmpS1534;
        _M0L6_2atmpS1536 = _M0L6offsetS338 - 1;
        _M0L6offsetS338 = _M0L6_2atmpS1536;
        _M0L1nS339 = _M0L1qS340;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB20int__to__string__hex(
  uint16_t* _M0L6bufferS320,
  uint32_t _M0L3numS325,
  int32_t _M0L12digit__startS321,
  int32_t _M0L10total__lenS324
) {
  int32_t _M0L6_2atmpS1523;
  int32_t _M0L6offsetS315;
  uint32_t _M0L1nS316;
  #line 29 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1523 = _M0L10total__lenS324 - _M0L12digit__startS321;
  _M0L6offsetS315 = _M0L6_2atmpS1523;
  _M0L1nS316 = _M0L3numS325;
  while (1) {
    if (_M0L6offsetS315 >= 2) {
      uint32_t _M0L6_2atmpS1520 = _M0L1nS316 & 255u;
      int32_t _M0L9byte__valS317 = *(int32_t*)&_M0L6_2atmpS1520;
      int32_t _M0L2hiS318 = _M0L9byte__valS317 / 16;
      int32_t _M0L2loS319 = _M0L9byte__valS317 % 16;
      int32_t _M0L6_2atmpS1514 = _M0L12digit__startS321 + _M0L6offsetS315;
      int32_t _M0L6_2atmpS1512 = _M0L6_2atmpS1514 - 2;
      int32_t _M0L6_2atmpS1513 =
        ((moonbit_string_t)moonbit_string_literal_86.data)[_M0L2hiS318];
      int32_t _M0L6_2atmpS1517;
      int32_t _M0L6_2atmpS1515;
      int32_t _M0L6_2atmpS1516;
      int32_t _M0L6_2atmpS1518;
      uint32_t _M0L6_2atmpS1519;
      _M0L6bufferS320[_M0L6_2atmpS1512] = _M0L6_2atmpS1513;
      _M0L6_2atmpS1517 = _M0L12digit__startS321 + _M0L6offsetS315;
      _M0L6_2atmpS1515 = _M0L6_2atmpS1517 - 1;
      _M0L6_2atmpS1516
      = ((moonbit_string_t)moonbit_string_literal_86.data)[
        _M0L2loS319
      ];
      _M0L6bufferS320[_M0L6_2atmpS1515] = _M0L6_2atmpS1516;
      _M0L6_2atmpS1518 = _M0L6offsetS315 - 2;
      _M0L6_2atmpS1519 = _M0L1nS316 >> 8;
      _M0L6offsetS315 = _M0L6_2atmpS1518;
      _M0L1nS316 = _M0L6_2atmpS1519;
      continue;
    } else if (_M0L6offsetS315 == 1) {
      uint32_t _M0L6_2atmpS1522 = _M0L1nS316 & 15u;
      int32_t _M0L6nibbleS323 = *(int32_t*)&_M0L6_2atmpS1522;
      int32_t _M0L6_2atmpS1521 =
        ((moonbit_string_t)moonbit_string_literal_86.data)[_M0L6nibbleS323];
      _M0L6bufferS320[_M0L12digit__startS321] = _M0L6_2atmpS1521;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB4Iter4nextGcE(struct _M0TPB4IterGcE* _M0L4selfS298) {
  struct _M0TWEOc* _M0L7_2afuncS297;
  int32_t _M0L6resultS299;
  int64_t _M0L7_2abindS300;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS297 = _M0L4selfS298->$0;
  moonbit_incref(_M0L7_2afuncS297);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS299 = _M0L7_2afuncS297->code(_M0L7_2afuncS297);
  moonbit_decref(_M0L7_2afuncS297);
  _M0L7_2abindS300 = _M0L4selfS298->$1;
  if (_M0L6resultS299 == -1) {
    _M0L4selfS298->$1 = _M0MPB4Iter4nextN6constrS9981GcE;
  } else if (_M0L7_2abindS300 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS301 = _M0L7_2abindS300;
    int32_t _M0L4_2anS302 = (int32_t)_M0L7_2aSomeS301;
    int64_t _M0L6_2atmpS1506;
    if (_M0L4_2anS302 > 0) {
      int32_t _M0L6_2atmpS1507 = _M0L4_2anS302 - 1;
      _M0L6_2atmpS1506 = (int64_t)_M0L6_2atmpS1507;
    } else {
      _M0L6_2atmpS1506 = _M0MPB4Iter4nextN6constrS9980GcE;
    }
    _M0L4selfS298->$1 = _M0L6_2atmpS1506;
  }
  return _M0L6resultS299;
}

struct _M0TUsRPB5ArrayGiEE* _M0MPB4Iter4nextGUsRPB5ArrayGiEEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0L4selfS304
) {
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L7_2afuncS303;
  struct _M0TUsRPB5ArrayGiEE* _M0L6resultS305;
  int64_t _M0L7_2abindS306;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS303 = _M0L4selfS304->$0;
  moonbit_incref(_M0L7_2afuncS303);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS305 = _M0L7_2afuncS303->code(_M0L7_2afuncS303);
  moonbit_decref(_M0L7_2afuncS303);
  _M0L7_2abindS306 = _M0L4selfS304->$1;
  if (_M0L6resultS305 == 0) {
    _M0L4selfS304->$1 = _M0MPB4Iter4nextN6constrS9981GUsRPB5ArrayGiEEE;
  } else if (_M0L7_2abindS306 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS307 = _M0L7_2abindS306;
    int32_t _M0L4_2anS308 = (int32_t)_M0L7_2aSomeS307;
    int64_t _M0L6_2atmpS1508;
    if (_M0L4_2anS308 > 0) {
      int32_t _M0L6_2atmpS1509 = _M0L4_2anS308 - 1;
      _M0L6_2atmpS1508 = (int64_t)_M0L6_2atmpS1509;
    } else {
      _M0L6_2atmpS1508 = _M0MPB4Iter4nextN6constrS9980GUsRPB5ArrayGiEEE;
    }
    _M0L4selfS304->$1 = _M0L6_2atmpS1508;
  }
  return _M0L6resultS305;
}

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS310
) {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L7_2afuncS309;
  void* _M0L6resultS311;
  int64_t _M0L7_2abindS312;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS309 = _M0L4selfS310->$0;
  moonbit_incref(_M0L7_2afuncS309);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS311 = _M0L7_2afuncS309->code(_M0L7_2afuncS309);
  moonbit_decref(_M0L7_2afuncS309);
  _M0L7_2abindS312 = _M0L4selfS310->$1;
  switch (Moonbit_object_tag(_M0L6resultS311)) {
    case 1: {
      if (_M0L7_2abindS312 == 4294967296ll) {
        
      } else {
        int64_t _M0L7_2aSomeS313 = _M0L7_2abindS312;
        int32_t _M0L4_2anS314 = (int32_t)_M0L7_2aSomeS313;
        int64_t _M0L6_2atmpS1510;
        if (_M0L4_2anS314 > 0) {
          int32_t _M0L6_2atmpS1511 = _M0L4_2anS314 - 1;
          _M0L6_2atmpS1510 = (int64_t)_M0L6_2atmpS1511;
        } else {
          _M0L6_2atmpS1510
          = _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE;
        }
        _M0L4selfS310->$1 = _M0L6_2atmpS1510;
      }
      break;
    }
    default: {
      _M0L4selfS310->$1
      = _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE;
      break;
    }
  }
  return _M0L6resultS311;
}

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void* _M0L4selfS296
) {
  struct _M0TPB13StringBuilder* _M0L6loggerS295;
  struct _M0TPB6Logger _M0L6_2atmpS1505;
  moonbit_string_t _result_3167;
  #line 165 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 166 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS295 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  moonbit_incref(_M0L6loggerS295);
  _M0L6_2atmpS1505
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L6loggerS295
  };
  #line 167 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB7FailurePB4Show6output(_M0L4selfS296, _M0L6_2atmpS1505);
  if (_M0L6_2atmpS1505.$1) {
    moonbit_decref(_M0L6_2atmpS1505.$1);
  }
  #line 168 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _result_3167 = _M0MPB13StringBuilder10to__string(_M0L6loggerS295);
  moonbit_decref(_M0L6loggerS295);
  return _result_3167;
}

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t _M0L4selfS292,
  struct _M0TPB6Logger _M0L6loggerS291
) {
  moonbit_string_t _M0L6_2atmpS1503;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1503 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS292);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS291.$0->$method_0(_M0L6loggerS291.$1, _M0L6_2atmpS1503);
  moonbit_decref(_M0L6_2atmpS1503);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t _M0L4selfS294,
  struct _M0TPB6Logger _M0L6loggerS293
) {
  moonbit_string_t _M0L6_2atmpS1504;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1504 = _M0IPC13int3IntPB4Show10to__string(_M0L4selfS294);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS293.$0->$method_0(_M0L6loggerS293.$1, _M0L6_2atmpS1504);
  moonbit_decref(_M0L6_2atmpS1504);
  return 0;
}

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView _M0L4selfS290
) {
  #line 99 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  return _M0L4selfS290.$1;
}

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView _M0L4selfS289
) {
  moonbit_string_t _M0L8_2afieldS2857;
  #line 92 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L8_2afieldS2857 = _M0L4selfS289.$0;
  moonbit_incref(_M0L8_2afieldS2857);
  return _M0L8_2afieldS2857;
}

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS285,
  moonbit_string_t _M0L5valueS286,
  int32_t _M0L5startS287,
  int32_t _M0L3lenS288
) {
  int32_t _M0L6_2atmpS1502;
  int64_t _M0L6_2atmpS1501;
  struct _M0TPC16string10StringView _M0L6_2atmpS1500;
  #line 122 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1502 = _M0L5startS287 + _M0L3lenS288;
  _M0L6_2atmpS1501 = (int64_t)_M0L6_2atmpS1502;
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1500
  = _M0MPC16string6String11sub_2einner(_M0L5valueS286, _M0L5startS287, _M0L6_2atmpS1501);
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L4selfS285, _M0L6_2atmpS1500);
  moonbit_decref(_M0L6_2atmpS1500.$0);
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t _M0L4selfS278,
  int32_t _M0L5startS284,
  int64_t _M0L3endS280
) {
  int32_t _M0L3lenS277;
  int32_t _M0L3endS279;
  int32_t _M0L5startS283;
  int32_t _if__result_3168;
  #line 755 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3lenS277 = Moonbit_array_length(_M0L4selfS278);
  if (_M0L3endS280 == 4294967296ll) {
    _M0L3endS279 = _M0L3lenS277;
  } else {
    int64_t _M0L7_2aSomeS281 = _M0L3endS280;
    int32_t _M0L6_2aendS282 = (int32_t)_M0L7_2aSomeS281;
    if (_M0L6_2aendS282 < 0) {
      _M0L3endS279 = _M0L3lenS277 + _M0L6_2aendS282;
    } else {
      _M0L3endS279 = _M0L6_2aendS282;
    }
  }
  if (_M0L5startS284 < 0) {
    _M0L5startS283 = _M0L3lenS277 + _M0L5startS284;
  } else {
    _M0L5startS283 = _M0L5startS284;
  }
  if (_M0L5startS283 >= 0) {
    if (_M0L5startS283 <= _M0L3endS279) {
      _if__result_3168 = _M0L3endS279 <= _M0L3lenS277;
    } else {
      _if__result_3168 = 0;
    }
  } else {
    _if__result_3168 = 0;
  }
  if (_if__result_3168) {
    if (_M0L5startS283 < _M0L3lenS277) {
      int32_t _M0L6_2atmpS1497 = _M0L4selfS278[_M0L5startS283];
      int32_t _M0L6_2atmpS1496;
      #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1496
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1497);
      if (!_M0L6_2atmpS1496) {
        
      } else {
        #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L3endS279 < _M0L3lenS277) {
      int32_t _M0L6_2atmpS1499 = _M0L4selfS278[_M0L3endS279];
      int32_t _M0L6_2atmpS1498;
      #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1498
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1499);
      if (!_M0L6_2atmpS1498) {
        
      } else {
        #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    moonbit_incref(_M0L4selfS278);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS278,
                                                 .$1 = _M0L5startS283,
                                                 .$2 = _M0L3endS279};
  } else {
    #line 763 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS276,
  struct _M0TPB4Show _M0L4showS275
) {
  struct _M0TPB6Logger _M0L6_2atmpS1495;
  #line 116 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS276);
  _M0L6_2atmpS1495
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS276
  };
  #line 117 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS275.$0->$method_0(_M0L4showS275.$1, _M0L6_2atmpS1495);
  if (_M0L6_2atmpS1495.$1) {
    moonbit_decref(_M0L6_2atmpS1495.$1);
  }
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS274,
  struct _M0TPB4Show _M0L4showS273
) {
  struct _M0TPB6Logger _M0L6_2atmpS1494;
  #line 111 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS274);
  _M0L6_2atmpS1494
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS274
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS273.$0->$method_0(_M0L4showS273.$1, _M0L6_2atmpS1494);
  if (_M0L6_2atmpS1494.$1) {
    moonbit_decref(_M0L6_2atmpS1494.$1);
  }
  return 0;
}

int32_t _M0FPB13finalize__acc(uint32_t _M0L3accS272) {
  uint32_t _M0L6_2atmpS1493;
  #line 444 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  #line 445 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1493 = _M0FPB14avalanche__acc(_M0L3accS272);
  return *(int32_t*)&_M0L6_2atmpS1493;
}

uint32_t _M0FPB14avalanche__acc(uint32_t _M0L3accS271) {
  uint32_t _M0Lm3accS270;
  uint32_t _M0L6_2atmpS1482;
  uint32_t _M0L6_2atmpS1484;
  uint32_t _M0L6_2atmpS1483;
  uint32_t _M0L6_2atmpS1485;
  uint32_t _M0L6_2atmpS1486;
  uint32_t _M0L6_2atmpS1488;
  uint32_t _M0L6_2atmpS1487;
  uint32_t _M0L6_2atmpS1489;
  uint32_t _M0L6_2atmpS1490;
  uint32_t _M0L6_2atmpS1492;
  uint32_t _M0L6_2atmpS1491;
  #line 449 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0Lm3accS270 = _M0L3accS271;
  _M0L6_2atmpS1482 = _M0Lm3accS270;
  _M0L6_2atmpS1484 = _M0Lm3accS270;
  _M0L6_2atmpS1483 = _M0L6_2atmpS1484 >> 15;
  _M0Lm3accS270 = _M0L6_2atmpS1482 ^ _M0L6_2atmpS1483;
  _M0L6_2atmpS1485 = _M0Lm3accS270;
  _M0Lm3accS270 = _M0L6_2atmpS1485 * 2246822519u;
  _M0L6_2atmpS1486 = _M0Lm3accS270;
  _M0L6_2atmpS1488 = _M0Lm3accS270;
  _M0L6_2atmpS1487 = _M0L6_2atmpS1488 >> 13;
  _M0Lm3accS270 = _M0L6_2atmpS1486 ^ _M0L6_2atmpS1487;
  _M0L6_2atmpS1489 = _M0Lm3accS270;
  _M0Lm3accS270 = _M0L6_2atmpS1489 * 3266489917u;
  _M0L6_2atmpS1490 = _M0Lm3accS270;
  _M0L6_2atmpS1492 = _M0Lm3accS270;
  _M0L6_2atmpS1491 = _M0L6_2atmpS1492 >> 16;
  _M0Lm3accS270 = _M0L6_2atmpS1490 ^ _M0L6_2atmpS1491;
  return _M0Lm3accS270;
}

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder* _M0L4selfS269,
  moonbit_string_t _M0L3strS268
) {
  int32_t _M0L8str__lenS267;
  int32_t _M0L3lenS1477;
  int32_t _M0L6_2atmpS1476;
  uint16_t* _M0L4dataS1478;
  int32_t _M0L3lenS1479;
  int32_t _M0L3lenS1481;
  int32_t _M0L6_2atmpS1480;
  #line 86 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L8str__lenS267 = Moonbit_array_length(_M0L3strS268);
  _M0L3lenS1477 = _M0L4selfS269->$1;
  _M0L6_2atmpS1476 = _M0L3lenS1477 + _M0L8str__lenS267;
  #line 88 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS269, _M0L6_2atmpS1476);
  _M0L4dataS1478 = _M0L4selfS269->$0;
  _M0L3lenS1479 = _M0L4selfS269->$1;
  moonbit_incref(_M0L4dataS1478);
  #line 89 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1478, _M0L3lenS1479, _M0L3strS268, 0, _M0L8str__lenS267);
  moonbit_decref(_M0L4dataS1478);
  _M0L3lenS1481 = _M0L4selfS269->$1;
  _M0L6_2atmpS1480 = _M0L3lenS1481 + _M0L8str__lenS267;
  _M0L4selfS269->$1 = _M0L6_2atmpS1480;
  return 0;
}

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t* _M0L4selfS263,
  int32_t _M0L11dst__offsetS266,
  moonbit_string_t _M0L3strS264,
  int32_t _M0L11str__offsetS259,
  int32_t _M0L3lenS260
) {
  int32_t _M0L16end__str__offsetS258;
  int32_t _M0L1iS261;
  int32_t _M0L1jS262;
  #line 71 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L16end__str__offsetS258 = _M0L11str__offsetS259 + _M0L3lenS260;
  _M0L1iS261 = _M0L11str__offsetS259;
  _M0L1jS262 = _M0L11dst__offsetS266;
  while (1) {
    if (_M0L1iS261 < _M0L16end__str__offsetS258) {
      int32_t _M0L6_2atmpS1473 = _M0L3strS264[_M0L1iS261];
      int32_t _M0L6_2atmpS1474;
      int32_t _M0L6_2atmpS1475;
      if (
        _M0L1jS262 < 0 || _M0L1jS262 >= Moonbit_array_length(_M0L4selfS263)
      ) {
        #line 80 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
        moonbit_panic();
      }
      _M0L4selfS263[_M0L1jS262] = _M0L6_2atmpS1473;
      _M0L6_2atmpS1474 = _M0L1iS261 + 1;
      _M0L6_2atmpS1475 = _M0L1jS262 + 1;
      _M0L1iS261 = _M0L6_2atmpS1474;
      _M0L1jS262 = _M0L6_2atmpS1475;
      continue;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t _M0L4selfS256,
  int32_t _M0L5quoteS257
) {
  struct _M0TPB13StringBuilder* _M0L3bufS255;
  int32_t _M0L6_2atmpS1472;
  struct _M0TPC16string10StringView _M0L6_2atmpS1470;
  struct _M0TPB6Logger _M0L6_2atmpS1471;
  moonbit_string_t _result_3170;
  #line 110 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 111 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L3bufS255 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS1472 = Moonbit_array_length(_M0L4selfS256);
  moonbit_incref(_M0L4selfS256);
  _M0L6_2atmpS1470
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS256, .$1 = 0, .$2 = _M0L6_2atmpS1472
  };
  moonbit_incref(_M0L3bufS255);
  _M0L6_2atmpS1471
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS255
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0MPC16string10StringView18escape__to_2einner(_M0L6_2atmpS1470, _M0L6_2atmpS1471, _M0L5quoteS257);
  moonbit_decref(_M0L6_2atmpS1470.$0);
  if (_M0L6_2atmpS1471.$1) {
    moonbit_decref(_M0L6_2atmpS1471.$1);
  }
  #line 113 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_3170 = _M0MPB13StringBuilder10to__string(_M0L3bufS255);
  moonbit_decref(_M0L3bufS255);
  return _result_3170;
}

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView _M0L4selfS247,
  struct _M0TPB6Logger _M0L6loggerS245,
  int32_t _M0L5quoteS244
) {
  int32_t _M0L3endS1468;
  int32_t _M0L5startS1469;
  int32_t _M0L3lenS246;
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS248;
  int32_t _M0L1iS249;
  int32_t _M0L3segS250;
  #line 144 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L5quoteS244) {
    #line 150 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS245.$0->$method_3(_M0L6loggerS245.$1, 34);
  }
  _M0L3endS1468 = _M0L4selfS247.$2;
  _M0L5startS1469 = _M0L4selfS247.$1;
  _M0L3lenS246 = _M0L3endS1468 - _M0L5startS1469;
  moonbit_incref(_M0L4selfS247.$0);
  if (_M0L6loggerS245.$1) {
    moonbit_incref(_M0L6loggerS245.$1);
  }
  _M0L6_2aenvS248
  = (struct _M0TURPC16string10StringViewRPB6LoggerE*)moonbit_malloc(sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE));
  Moonbit_object_header(_M0L6_2aenvS248)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 119, 0);
  _M0L6_2aenvS248->$0 = _M0L4selfS247;
  _M0L6_2aenvS248->$1 = _M0L6loggerS245;
  _M0L1iS249 = 0;
  _M0L3segS250 = 0;
  _2afor_251:;
  while (1) {
    moonbit_string_t _M0L3strS1465;
    int32_t _M0L5startS1467;
    int32_t _M0L6_2atmpS1466;
    int32_t _M0L4codeS252;
    int32_t _M0L1cS254;
    int32_t _M0L6_2atmpS1449;
    int32_t _M0L6_2atmpS1450;
    int32_t _M0L6_2atmpS1451;
    int32_t _tmp_3174;
    int32_t _tmp_3175;
    if (_M0L1iS249 >= _M0L3lenS246) {
      #line 160 "/home/developer/.moon/lib/core/builtin/show.mbt"
      _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
      moonbit_decref(_M0L6_2aenvS248);
      break;
    }
    _M0L3strS1465 = _M0L4selfS247.$0;
    _M0L5startS1467 = _M0L4selfS247.$1;
    _M0L6_2atmpS1466 = _M0L5startS1467 + _M0L1iS249;
    _M0L4codeS252 = _M0L3strS1465[_M0L6_2atmpS1466];
    switch (_M0L4codeS252) {
      case 34: {
        _M0L1cS254 = _M0L4codeS252;
        goto join_253;
        break;
      }
      
      case 92: {
        _M0L1cS254 = _M0L4codeS252;
        goto join_253;
        break;
      }
      
      case 10: {
        int32_t _M0L6_2atmpS1452;
        int32_t _M0L6_2atmpS1453;
        #line 172 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
        #line 173 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS245.$0->$method_0(_M0L6loggerS245.$1, (moonbit_string_t)moonbit_string_literal_59.data);
        _M0L6_2atmpS1452 = _M0L1iS249 + 1;
        _M0L6_2atmpS1453 = _M0L1iS249 + 1;
        _M0L1iS249 = _M0L6_2atmpS1452;
        _M0L3segS250 = _M0L6_2atmpS1453;
        goto _2afor_251;
        break;
      }
      
      case 13: {
        int32_t _M0L6_2atmpS1454;
        int32_t _M0L6_2atmpS1455;
        #line 177 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
        #line 178 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS245.$0->$method_0(_M0L6loggerS245.$1, (moonbit_string_t)moonbit_string_literal_87.data);
        _M0L6_2atmpS1454 = _M0L1iS249 + 1;
        _M0L6_2atmpS1455 = _M0L1iS249 + 1;
        _M0L1iS249 = _M0L6_2atmpS1454;
        _M0L3segS250 = _M0L6_2atmpS1455;
        goto _2afor_251;
        break;
      }
      
      case 8: {
        int32_t _M0L6_2atmpS1456;
        int32_t _M0L6_2atmpS1457;
        #line 182 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
        #line 183 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS245.$0->$method_0(_M0L6loggerS245.$1, (moonbit_string_t)moonbit_string_literal_88.data);
        _M0L6_2atmpS1456 = _M0L1iS249 + 1;
        _M0L6_2atmpS1457 = _M0L1iS249 + 1;
        _M0L1iS249 = _M0L6_2atmpS1456;
        _M0L3segS250 = _M0L6_2atmpS1457;
        goto _2afor_251;
        break;
      }
      
      case 9: {
        int32_t _M0L6_2atmpS1458;
        int32_t _M0L6_2atmpS1459;
        #line 187 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
        #line 188 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS245.$0->$method_0(_M0L6loggerS245.$1, (moonbit_string_t)moonbit_string_literal_60.data);
        _M0L6_2atmpS1458 = _M0L1iS249 + 1;
        _M0L6_2atmpS1459 = _M0L1iS249 + 1;
        _M0L1iS249 = _M0L6_2atmpS1458;
        _M0L3segS250 = _M0L6_2atmpS1459;
        goto _2afor_251;
        break;
      }
      default: {
        if (_M0L4codeS252 < 32) {
          int32_t _M0L6_2atmpS1461;
          moonbit_string_t _M0L6_2atmpS1460;
          int32_t _M0L6_2atmpS1462;
          int32_t _M0L6_2atmpS1463;
          #line 193 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
          #line 194 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS245.$0->$method_0(_M0L6loggerS245.$1, (moonbit_string_t)moonbit_string_literal_89.data);
          _M0L6_2atmpS1461 = _M0L4codeS252 & 0xff;
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6_2atmpS1460 = _M0MPC14byte4Byte7to__hex(_M0L6_2atmpS1461);
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS245.$0->$method_0(_M0L6loggerS245.$1, _M0L6_2atmpS1460);
          moonbit_decref(_M0L6_2atmpS1460);
          #line 196 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS245.$0->$method_3(_M0L6loggerS245.$1, 125);
          _M0L6_2atmpS1462 = _M0L1iS249 + 1;
          _M0L6_2atmpS1463 = _M0L1iS249 + 1;
          _M0L1iS249 = _M0L6_2atmpS1462;
          _M0L3segS250 = _M0L6_2atmpS1463;
          goto _2afor_251;
        } else {
          int32_t _M0L6_2atmpS1464 = _M0L1iS249 + 1;
          int32_t _tmp_3173 = _M0L3segS250;
          _M0L1iS249 = _M0L6_2atmpS1464;
          _M0L3segS250 = _tmp_3173;
          goto _2afor_251;
        }
        break;
      }
    }
    goto joinlet_3172;
    join_253:;
    #line 166 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS248, _M0L3segS250, _M0L1iS249);
    #line 167 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS245.$0->$method_3(_M0L6loggerS245.$1, 92);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1449 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS254);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS245.$0->$method_3(_M0L6loggerS245.$1, _M0L6_2atmpS1449);
    _M0L6_2atmpS1450 = _M0L1iS249 + 1;
    _M0L6_2atmpS1451 = _M0L1iS249 + 1;
    _M0L1iS249 = _M0L6_2atmpS1450;
    _M0L3segS250 = _M0L6_2atmpS1451;
    continue;
    joinlet_3172:;
    _tmp_3174 = _M0L1iS249;
    _tmp_3175 = _M0L3segS250;
    _M0L1iS249 = _tmp_3174;
    _M0L3segS250 = _tmp_3175;
    continue;
    break;
  }
  if (_M0L5quoteS244) {
    #line 204 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS245.$0->$method_3(_M0L6loggerS245.$1, 34);
  }
  return 0;
}

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS240,
  int32_t _M0L3segS243,
  int32_t _M0L1iS242
) {
  struct _M0TPB6Logger _M0L6loggerS239;
  struct _M0TPC16string10StringView _M0L4selfS241;
  #line 153 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6loggerS239 = _M0L6_2aenvS240->$1;
  _M0L4selfS241 = _M0L6_2aenvS240->$0;
  if (_M0L1iS242 > _M0L3segS243) {
    int64_t _M0L6_2atmpS1448 = (int64_t)_M0L1iS242;
    struct _M0TPC16string10StringView _M0L6_2atmpS1447;
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1447
    = _M0MPC16string10StringView11sub_2einner(_M0L4selfS241, _M0L3segS243, _M0L6_2atmpS1448);
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS239.$0->$method_2(_M0L6loggerS239.$1, _M0L6_2atmpS1447);
    moonbit_decref(_M0L6_2atmpS1447.$0);
  }
  return 0;
}

int32_t _M0MPC16string10StringView11unsafe__get(
  struct _M0TPC16string10StringView _M0L4selfS237,
  int32_t _M0L5indexS238
) {
  moonbit_string_t _M0L3strS1444;
  int32_t _M0L5startS1446;
  int32_t _M0L6_2atmpS1445;
  #line 128 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1444 = _M0L4selfS237.$0;
  _M0L5startS1446 = _M0L4selfS237.$1;
  _M0L6_2atmpS1445 = _M0L5startS1446 + _M0L5indexS238;
  return _M0L3strS1444[_M0L6_2atmpS1445];
}

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView _M0L4selfS230,
  int32_t _M0L5startS236,
  int64_t _M0L3endS232
) {
  moonbit_string_t _M0L3strS1443;
  int32_t _M0L8str__lenS229;
  int32_t _M0L8abs__endS231;
  int32_t _M0L10abs__startS235;
  int32_t _M0L5startS1431;
  int32_t _if__result_3176;
  #line 814 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1443 = _M0L4selfS230.$0;
  _M0L8str__lenS229 = Moonbit_array_length(_M0L3strS1443);
  if (_M0L3endS232 == 4294967296ll) {
    _M0L8abs__endS231 = _M0L4selfS230.$2;
  } else {
    int64_t _M0L7_2aSomeS233 = _M0L3endS232;
    int32_t _M0L6_2aendS234 = (int32_t)_M0L7_2aSomeS233;
    if (_M0L6_2aendS234 < 0) {
      int32_t _M0L3endS1441 = _M0L4selfS230.$2;
      _M0L8abs__endS231 = _M0L3endS1441 + _M0L6_2aendS234;
    } else {
      int32_t _M0L5startS1442 = _M0L4selfS230.$1;
      _M0L8abs__endS231 = _M0L5startS1442 + _M0L6_2aendS234;
    }
  }
  if (_M0L5startS236 < 0) {
    int32_t _M0L3endS1439 = _M0L4selfS230.$2;
    _M0L10abs__startS235 = _M0L3endS1439 + _M0L5startS236;
  } else {
    int32_t _M0L5startS1440 = _M0L4selfS230.$1;
    _M0L10abs__startS235 = _M0L5startS1440 + _M0L5startS236;
  }
  _M0L5startS1431 = _M0L4selfS230.$1;
  if (_M0L10abs__startS235 >= _M0L5startS1431) {
    if (_M0L10abs__startS235 <= _M0L8abs__endS231) {
      int32_t _M0L3endS1430 = _M0L4selfS230.$2;
      _if__result_3176 = _M0L8abs__endS231 <= _M0L3endS1430;
    } else {
      _if__result_3176 = 0;
    }
  } else {
    _if__result_3176 = 0;
  }
  if (_if__result_3176) {
    moonbit_string_t _M0L3strS1438;
    if (_M0L10abs__startS235 < _M0L8str__lenS229) {
      moonbit_string_t _M0L3strS1434 = _M0L4selfS230.$0;
      int32_t _M0L6_2atmpS1433 = _M0L3strS1434[_M0L10abs__startS235];
      int32_t _M0L6_2atmpS1432;
      #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1432
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1433);
      if (!_M0L6_2atmpS1432) {
        
      } else {
        #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L8abs__endS231 < _M0L8str__lenS229) {
      moonbit_string_t _M0L3strS1437 = _M0L4selfS230.$0;
      int32_t _M0L6_2atmpS1436 = _M0L3strS1437[_M0L8abs__endS231];
      int32_t _M0L6_2atmpS1435;
      #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1435
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1436);
      if (!_M0L6_2atmpS1435) {
        
      } else {
        #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    _M0L3strS1438 = _M0L4selfS230.$0;
    moonbit_incref(_M0L3strS1438);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1438,
                                                 .$1 = _M0L10abs__startS235,
                                                 .$2 = _M0L8abs__endS231};
  } else {
    #line 834 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC16string10StringView6length(
  struct _M0TPC16string10StringView _M0L4selfS228
) {
  int32_t _M0L3endS1428;
  int32_t _M0L5startS1429;
  #line 49 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3endS1428 = _M0L4selfS228.$2;
  _M0L5startS1429 = _M0L4selfS228.$1;
  return _M0L3endS1428 - _M0L5startS1429;
}

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t _M0L1bS227) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS226;
  int32_t _M0L6_2atmpS1425;
  int32_t _M0L6_2atmpS1424;
  int32_t _M0L6_2atmpS1427;
  int32_t _M0L6_2atmpS1426;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1423;
  moonbit_string_t _result_3177;
  #line 74 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L7_2aselfS226 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1425 = _M0IPC14byte4BytePB3Div3div(_M0L1bS227, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1424
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1425);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS226, _M0L6_2atmpS1424);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1427 = _M0IPC14byte4BytePB3Mod3mod(_M0L1bS227, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1426
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1427);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS226, _M0L6_2atmpS1426);
  _M0L6_2atmpS1423 = _M0L7_2aselfS226;
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_3177 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1423);
  moonbit_decref(_M0L6_2atmpS1423);
  return _result_3177;
}

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t _M0L1iS225) {
  #line 75 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L1iS225 < 10) {
    int32_t _M0L6_2atmpS1420;
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1420 = _M0IPC14byte4BytePB3Add3add(_M0L1iS225, 48);
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1420);
  } else {
    int32_t _M0L6_2atmpS1422;
    int32_t _M0L6_2atmpS1421;
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1422 = _M0IPC14byte4BytePB3Add3add(_M0L1iS225, 97);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1421 = _M0IPC14byte4BytePB3Sub3sub(_M0L6_2atmpS1422, 10);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1421);
  }
}

int32_t _M0IPC14byte4BytePB3Sub3sub(
  int32_t _M0L4selfS223,
  int32_t _M0L4thatS224
) {
  int32_t _M0L6_2atmpS1418;
  int32_t _M0L6_2atmpS1419;
  int32_t _M0L6_2atmpS1417;
  #line 120 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1418 = (int32_t)_M0L4selfS223;
  _M0L6_2atmpS1419 = (int32_t)_M0L4thatS224;
  _M0L6_2atmpS1417 = _M0L6_2atmpS1418 - _M0L6_2atmpS1419;
  return _M0L6_2atmpS1417 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Mod3mod(
  int32_t _M0L4selfS221,
  int32_t _M0L4thatS222
) {
  int32_t _M0L6_2atmpS1415;
  int32_t _M0L6_2atmpS1416;
  int32_t _M0L6_2atmpS1414;
  #line 67 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1415 = (int32_t)_M0L4selfS221;
  _M0L6_2atmpS1416 = (int32_t)_M0L4thatS222;
  _M0L6_2atmpS1414 = _M0L6_2atmpS1415 % _M0L6_2atmpS1416;
  return _M0L6_2atmpS1414 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Div3div(
  int32_t _M0L4selfS219,
  int32_t _M0L4thatS220
) {
  int32_t _M0L6_2atmpS1412;
  int32_t _M0L6_2atmpS1413;
  int32_t _M0L6_2atmpS1411;
  #line 62 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1412 = (int32_t)_M0L4selfS219;
  _M0L6_2atmpS1413 = (int32_t)_M0L4thatS220;
  _M0L6_2atmpS1411 = _M0L6_2atmpS1412 / _M0L6_2atmpS1413;
  return _M0L6_2atmpS1411 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Add3add(
  int32_t _M0L4selfS217,
  int32_t _M0L4thatS218
) {
  int32_t _M0L6_2atmpS1409;
  int32_t _M0L6_2atmpS1410;
  int32_t _M0L6_2atmpS1408;
  #line 106 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1409 = (int32_t)_M0L4selfS217;
  _M0L6_2atmpS1410 = (int32_t)_M0L4thatS218;
  _M0L6_2atmpS1408 = _M0L6_2atmpS1409 + _M0L6_2atmpS1410;
  return _M0L6_2atmpS1408 & 0xff;
}

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t _M0L4selfS216) {
  int32_t _M0L6_2atmpS1407;
  #line 68 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  _M0L6_2atmpS1407 = (int32_t)_M0L4selfS216;
  return _M0L6_2atmpS1407;
}

int32_t _M0FPB32code__point__of__surrogate__pair(
  int32_t _M0L7leadingS214,
  int32_t _M0L8trailingS215
) {
  int32_t _M0L6_2atmpS1406;
  int32_t _M0L6_2atmpS1405;
  int32_t _M0L6_2atmpS1404;
  int32_t _M0L6_2atmpS1403;
  int32_t _M0L6_2atmpS1402;
  #line 40 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1406 = _M0L7leadingS214 - 55296;
  _M0L6_2atmpS1405 = _M0L6_2atmpS1406 * 1024;
  _M0L6_2atmpS1404 = _M0L6_2atmpS1405 + _M0L8trailingS215;
  _M0L6_2atmpS1403 = _M0L6_2atmpS1404 - 56320;
  _M0L6_2atmpS1402 = _M0L6_2atmpS1403 + 65536;
  return _M0L6_2atmpS1402;
}

int32_t _M0MPC16string6String20char__length_2einner(
  moonbit_string_t _M0L4selfS207,
  int32_t _M0L13start__offsetS208,
  int64_t _M0L11end__offsetS205
) {
  int32_t _M0L11end__offsetS204;
  int32_t _if__result_3178;
  #line 60 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L11end__offsetS205 == 4294967296ll) {
    _M0L11end__offsetS204 = Moonbit_array_length(_M0L4selfS207);
  } else {
    int64_t _M0L7_2aSomeS206 = _M0L11end__offsetS205;
    _M0L11end__offsetS204 = (int32_t)_M0L7_2aSomeS206;
  }
  if (_M0L13start__offsetS208 >= 0) {
    if (_M0L13start__offsetS208 <= _M0L11end__offsetS204) {
      int32_t _M0L6_2atmpS1395 = Moonbit_array_length(_M0L4selfS207);
      _if__result_3178 = _M0L11end__offsetS204 <= _M0L6_2atmpS1395;
    } else {
      _if__result_3178 = 0;
    }
  } else {
    _if__result_3178 = 0;
  }
  if (_if__result_3178) {
    int32_t _M0L12utf16__indexS209 = _M0L13start__offsetS208;
    int32_t _M0L11char__countS210 = 0;
    while (1) {
      if (_M0L12utf16__indexS209 < _M0L11end__offsetS204) {
        int32_t _M0L2c1S211 = _M0L4selfS207[_M0L12utf16__indexS209];
        int32_t _if__result_3180;
        int32_t _M0L6_2atmpS1400;
        int32_t _M0L6_2atmpS1401;
        #line 76 "/home/developer/.moon/lib/core/builtin/string.mbt"
        if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S211)) {
          int32_t _M0L6_2atmpS1396 = _M0L12utf16__indexS209 + 1;
          _if__result_3180 = _M0L6_2atmpS1396 < _M0L11end__offsetS204;
        } else {
          _if__result_3180 = 0;
        }
        if (_if__result_3180) {
          int32_t _M0L6_2atmpS1399 = _M0L12utf16__indexS209 + 1;
          int32_t _M0L2c2S212 = _M0L4selfS207[_M0L6_2atmpS1399];
          #line 78 "/home/developer/.moon/lib/core/builtin/string.mbt"
          if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S212)) {
            int32_t _M0L6_2atmpS1397 = _M0L12utf16__indexS209 + 2;
            int32_t _M0L6_2atmpS1398 = _M0L11char__countS210 + 1;
            _M0L12utf16__indexS209 = _M0L6_2atmpS1397;
            _M0L11char__countS210 = _M0L6_2atmpS1398;
            continue;
          } else {
            #line 81 "/home/developer/.moon/lib/core/builtin/string.mbt"
            _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_90.data);
          }
        }
        _M0L6_2atmpS1400 = _M0L12utf16__indexS209 + 1;
        _M0L6_2atmpS1401 = _M0L11char__countS210 + 1;
        _M0L12utf16__indexS209 = _M0L6_2atmpS1400;
        _M0L11char__countS210 = _M0L6_2atmpS1401;
        continue;
      } else {
        return _M0L11char__countS210;
      }
      break;
    }
  } else {
    #line 70 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0FPC15abort5abortGiE((moonbit_string_t)moonbit_string_literal_91.data);
  }
}

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t _M0L4selfS203) {
  #line 45 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS203 >= 56320) {
    return _M0L4selfS203 <= 57343;
  } else {
    return 0;
  }
}

int32_t _M0MPC16uint166UInt1622is__leading__surrogate(int32_t _M0L4selfS202) {
  #line 28 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS202 >= 55296) {
    return _M0L4selfS202 <= 56319;
  } else {
    return 0;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder* _M0L4selfS200,
  int32_t _M0L2chS199
) {
  uint32_t _M0L4codeS198;
  #line 95 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  #line 96 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4codeS198 = _M0MPC14char4Char8to__uint(_M0L2chS199);
  if (_M0L4codeS198 <= 65535u) {
    int32_t _M0L3lenS1374 = _M0L4selfS200->$1;
    int32_t _M0L6_2atmpS1373 = _M0L3lenS1374 + 1;
    uint16_t* _M0L4dataS1375;
    int32_t _M0L3lenS1376;
    int32_t _M0L6_2atmpS1377;
    int32_t _M0L3lenS1379;
    int32_t _M0L6_2atmpS1378;
    #line 98 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS200, _M0L6_2atmpS1373);
    _M0L4dataS1375 = _M0L4selfS200->$0;
    _M0L3lenS1376 = _M0L4selfS200->$1;
    moonbit_incref(_M0L4dataS1375);
    #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1377 = _M0MPC14uint4UInt10to__uint16(_M0L4codeS198);
    if (
      _M0L3lenS1376 < 0
      || _M0L3lenS1376 >= Moonbit_array_length(_M0L4dataS1375)
    ) {
      #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1375[_M0L3lenS1376] = _M0L6_2atmpS1377;
    moonbit_decref(_M0L4dataS1375);
    _M0L3lenS1379 = _M0L4selfS200->$1;
    _M0L6_2atmpS1378 = _M0L3lenS1379 + 1;
    _M0L4selfS200->$1 = _M0L6_2atmpS1378;
  } else if (_M0L4codeS198 <= 1114111u) {
    int32_t _M0L3lenS1381 = _M0L4selfS200->$1;
    int32_t _M0L6_2atmpS1380 = _M0L3lenS1381 + 2;
    uint32_t _M0L4codeS201;
    uint16_t* _M0L4dataS1382;
    int32_t _M0L3lenS1383;
    uint32_t _M0L6_2atmpS1386;
    uint32_t _M0L6_2atmpS1385;
    int32_t _M0L6_2atmpS1384;
    uint16_t* _M0L4dataS1387;
    int32_t _M0L3lenS1392;
    int32_t _M0L6_2atmpS1388;
    uint32_t _M0L6_2atmpS1391;
    uint32_t _M0L6_2atmpS1390;
    int32_t _M0L6_2atmpS1389;
    int32_t _M0L3lenS1394;
    int32_t _M0L6_2atmpS1393;
    #line 102 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS200, _M0L6_2atmpS1380);
    _M0L4codeS201 = _M0L4codeS198 - 65536u;
    _M0L4dataS1382 = _M0L4selfS200->$0;
    _M0L3lenS1383 = _M0L4selfS200->$1;
    _M0L6_2atmpS1386 = _M0L4codeS201 >> 10;
    _M0L6_2atmpS1385 = 55296u + _M0L6_2atmpS1386;
    moonbit_incref(_M0L4dataS1382);
    #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1384 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1385);
    if (
      _M0L3lenS1383 < 0
      || _M0L3lenS1383 >= Moonbit_array_length(_M0L4dataS1382)
    ) {
      #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1382[_M0L3lenS1383] = _M0L6_2atmpS1384;
    moonbit_decref(_M0L4dataS1382);
    _M0L4dataS1387 = _M0L4selfS200->$0;
    _M0L3lenS1392 = _M0L4selfS200->$1;
    _M0L6_2atmpS1388 = _M0L3lenS1392 + 1;
    _M0L6_2atmpS1391 = _M0L4codeS201 & 1023u;
    _M0L6_2atmpS1390 = 56320u + _M0L6_2atmpS1391;
    moonbit_incref(_M0L4dataS1387);
    #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1389 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1390);
    if (
      _M0L6_2atmpS1388 < 0
      || _M0L6_2atmpS1388 >= Moonbit_array_length(_M0L4dataS1387)
    ) {
      #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1387[_M0L6_2atmpS1388] = _M0L6_2atmpS1389;
    moonbit_decref(_M0L4dataS1387);
    _M0L3lenS1394 = _M0L4selfS200->$1;
    _M0L6_2atmpS1393 = _M0L3lenS1394 + 2;
    _M0L4selfS200->$1 = _M0L6_2atmpS1393;
  } else {
    #line 108 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_92.data);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder* _M0L4selfS192,
  int32_t _M0L8requiredS193
) {
  uint16_t* _M0L4dataS1372;
  int32_t _M0L12current__lenS191;
  int32_t _M0L13enough__spaceS194;
  int32_t _M0L13enough__spaceS195;
  uint16_t* _M0L4dataS1368;
  int32_t _M0L6_2atmpS1369;
  int32_t _M0L3lenS1370;
  uint16_t* _M0L9new__dataS197;
  uint16_t* _M0L6_2aoldS2868;
  #line 46 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4dataS1372 = _M0L4selfS192->$0;
  _M0L12current__lenS191 = Moonbit_array_length(_M0L4dataS1372);
  if (_M0L8requiredS193 <= _M0L12current__lenS191) {
    return 0;
  }
  _M0L13enough__spaceS195 = _M0L12current__lenS191;
  while (1) {
    if (_M0L13enough__spaceS195 < _M0L8requiredS193) {
      int32_t _M0L6_2atmpS1371 = _M0L13enough__spaceS195 * 2;
      _M0L13enough__spaceS195 = _M0L6_2atmpS1371;
      continue;
    } else {
      _M0L13enough__spaceS194 = _M0L13enough__spaceS195;
    }
    break;
  }
  _M0L4dataS1368 = _M0L4selfS192->$0;
  moonbit_incref(_M0L4dataS1368);
  #line 64 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1369 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L3lenS1370 = _M0L4selfS192->$1;
  #line 61 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L9new__dataS197
  = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1368, _M0L13enough__spaceS194, _M0L6_2atmpS1369, _M0L3lenS1370, 0, 0);
  moonbit_decref(_M0L4dataS1368);
  _M0L6_2aoldS2868 = _M0L4selfS192->$0;
  moonbit_decref(_M0L6_2aoldS2868);
  _M0L4selfS192->$0 = _M0L9new__dataS197;
  return 0;
}

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t _M0L4selfS190) {
  int32_t _M0L6_2atmpS1367;
  #line 2676 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1367 = *(int32_t*)&_M0L4selfS190;
  return (uint16_t)_M0L6_2atmpS1367;
}

uint32_t _M0MPC14char4Char8to__uint(int32_t _M0L4selfS189) {
  int32_t _M0L6_2atmpS1366;
  #line 1254 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1366 = _M0L4selfS189;
  return *(uint32_t*)&_M0L6_2atmpS1366;
}

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder* _M0L4selfS187
) {
  int32_t _M0L3lenS1358;
  uint16_t* _M0L4dataS1360;
  int32_t _M0L6_2atmpS1359;
  #line 148 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3lenS1358 = _M0L4selfS187->$1;
  _M0L4dataS1360 = _M0L4selfS187->$0;
  _M0L6_2atmpS1359 = Moonbit_array_length(_M0L4dataS1360);
  if (_M0L3lenS1358 == _M0L6_2atmpS1359) {
    uint16_t* _M0L4dataS1361 = _M0L4selfS187->$0;
    moonbit_incref(_M0L4dataS1361);
    return _M0L4dataS1361;
  } else {
    uint16_t* _M0L4dataS1362 = _M0L4selfS187->$0;
    int32_t _M0L3lenS1363 = _M0L4selfS187->$1;
    int32_t _M0L6_2atmpS1364;
    int32_t _M0L3lenS1365;
    uint16_t* _M0L4dataS188;
    moonbit_incref(_M0L4dataS1362);
    #line 155 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1364 = _M0IPC16uint166UInt16PB7Default7default();
    _M0L3lenS1365 = _M0L4selfS187->$1;
    #line 152 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L4dataS188
    = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1362, _M0L3lenS1363, _M0L6_2atmpS1364, _M0L3lenS1365, 0, 0);
    moonbit_decref(_M0L4dataS1362);
    return _M0L4dataS188;
  }
}

int32_t _M0IPC16uint166UInt16PB7Default7default() {
  #line 176 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  return 0;
}

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t* _M0L3srcS184,
  int32_t _M0L13allocate__lenS180,
  int32_t _M0L4initS185,
  int32_t _M0L3lenS181,
  int32_t _M0L11src__offsetS182,
  int32_t _M0L11dst__offsetS183
) {
  int32_t _if__result_3182;
  #line 97 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L13allocate__lenS180 >= 0) {
    if (_M0L3lenS181 >= 0) {
      if (_M0L11src__offsetS182 >= 0) {
        if (_M0L11dst__offsetS183 >= 0) {
          int32_t _M0L6_2atmpS1354 = _M0L11src__offsetS182 + _M0L3lenS181;
          int32_t _M0L6_2atmpS1355 = Moonbit_array_length(_M0L3srcS184);
          if (_M0L6_2atmpS1354 <= _M0L6_2atmpS1355) {
            int32_t _M0L6_2atmpS1353 = _M0L11dst__offsetS183 + _M0L3lenS181;
            _if__result_3182 = _M0L6_2atmpS1353 <= _M0L13allocate__lenS180;
          } else {
            _if__result_3182 = 0;
          }
        } else {
          _if__result_3182 = 0;
        }
      } else {
        _if__result_3182 = 0;
      }
    } else {
      _if__result_3182 = 0;
    }
  } else {
    _if__result_3182 = 0;
  }
  if (_if__result_3182) {
    moonbit_incref(_M0L3srcS184);
    #line 115 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    return _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(_M0L3srcS184, _M0L13allocate__lenS180, _M0L4initS185, _M0L11src__offsetS182, _M0L11dst__offsetS183, _M0L3lenS181);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS186;
    int32_t _M0L6_2atmpS1357;
    moonbit_string_t _M0L6_2atmpS1356;
    uint16_t* _result_3183;
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L18_2astring__builderS186
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS186, (moonbit_string_t)moonbit_string_literal_93.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS186, _M0L13allocate__lenS180);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS186, (moonbit_string_t)moonbit_string_literal_94.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS186, _M0L11src__offsetS182);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS186, (moonbit_string_t)moonbit_string_literal_95.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS186, _M0L11dst__offsetS183);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS186, (moonbit_string_t)moonbit_string_literal_96.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS186, _M0L3lenS181);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS186, (moonbit_string_t)moonbit_string_literal_97.data);
    _M0L6_2atmpS1357 = Moonbit_array_length(_M0L3srcS184);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS186, _M0L6_2atmpS1357);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L6_2atmpS1356
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS186);
    moonbit_decref(_M0L18_2astring__builderS186);
    #line 111 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _result_3183 = _M0FPC15abort5abortGAkE(_M0L6_2atmpS1356);
    moonbit_decref(_M0L6_2atmpS1356);
    return _result_3183;
  }
}

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t* _M0L3srcS177,
  int32_t _M0L13allocate__lenS174,
  int32_t _M0L4initS175,
  int32_t _M0L11src__offsetS178,
  int32_t _M0L11dst__offsetS176,
  int32_t _M0L9blit__lenS179
) {
  uint16_t* _M0L3dstS173;
  #line 79 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  _M0L3dstS173
  = (uint16_t*)moonbit_make_string(_M0L13allocate__lenS174, _M0L4initS175);
  moonbit_incref(_M0L3dstS173);
  #line 90 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS173, _M0L11dst__offsetS176, _M0L3srcS177, _M0L11src__offsetS178, _M0L9blit__lenS179, sizeof(uint16_t));
  return _M0L3dstS173;
}

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t _M0L10size__hintS171
) {
  int32_t _M0L7initialS170;
  uint16_t* _M0L4dataS172;
  struct _M0TPB13StringBuilder* _block_3184;
  #line 32 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  if (_M0L10size__hintS171 < 1) {
    _M0L7initialS170 = 1;
  } else {
    int32_t _M0L6_2atmpS1352 = _M0L10size__hintS171 + 1;
    _M0L7initialS170 = _M0L6_2atmpS1352 / 2;
  }
  _M0L4dataS172 = (uint16_t*)moonbit_make_string(_M0L7initialS170, 0);
  _block_3184
  = (struct _M0TPB13StringBuilder*)moonbit_malloc(sizeof(struct _M0TPB13StringBuilder));
  Moonbit_object_header(_block_3184)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 124, 0);
  _block_3184->$0 = _M0L4dataS172;
  _block_3184->$1 = 0;
  return _block_3184;
}

int32_t _M0MPC14byte4Byte8to__char(int32_t _M0L4selfS169) {
  int32_t _M0L6_2atmpS1351;
  #line 1868 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1351 = (int32_t)_M0L4selfS169;
  return _M0L6_2atmpS1351;
}

int32_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGiE(
  int32_t* _M0L3srcS143,
  int32_t _M0L13allocate__lenS139,
  int32_t _M0L3lenS140,
  int32_t _M0L11src__offsetS141,
  int32_t _M0L11dst__offsetS142
) {
  int32_t _if__result_3185;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS139 >= 0) {
    if (_M0L3lenS140 >= 0) {
      if (_M0L11src__offsetS141 >= 0) {
        if (_M0L11dst__offsetS142 >= 0) {
          int32_t _M0L6_2atmpS1327 = _M0L11src__offsetS141 + _M0L3lenS140;
          int32_t _M0L6_2atmpS1328;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1328
          = _M0MPB18UninitializedArray6lengthGiE(_M0L3srcS143);
          if (_M0L6_2atmpS1327 <= _M0L6_2atmpS1328) {
            int32_t _M0L6_2atmpS1326 = _M0L11dst__offsetS142 + _M0L3lenS140;
            _if__result_3185 = _M0L6_2atmpS1326 <= _M0L13allocate__lenS139;
          } else {
            _if__result_3185 = 0;
          }
        } else {
          _if__result_3185 = 0;
        }
      } else {
        _if__result_3185 = 0;
      }
    } else {
      _if__result_3185 = 0;
    }
  } else {
    _if__result_3185 = 0;
  }
  if (_if__result_3185) {
    moonbit_incref(_M0L3srcS143);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return _M0MPB18UninitializedArray23unsafe__make__and__blitGiE(_M0L3srcS143, _M0L13allocate__lenS139, _M0L11src__offsetS141, _M0L11dst__offsetS142, _M0L3lenS140);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS144;
    int32_t _M0L6_2atmpS1330;
    moonbit_string_t _M0L6_2atmpS1329;
    int32_t* _result_3186;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS144
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS144, (moonbit_string_t)moonbit_string_literal_93.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS144, _M0L13allocate__lenS139);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS144, (moonbit_string_t)moonbit_string_literal_94.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS144, _M0L11src__offsetS141);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS144, (moonbit_string_t)moonbit_string_literal_95.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS144, _M0L11dst__offsetS142);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS144, (moonbit_string_t)moonbit_string_literal_96.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS144, _M0L3lenS140);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS144, (moonbit_string_t)moonbit_string_literal_97.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1330 = _M0MPB18UninitializedArray6lengthGiE(_M0L3srcS143);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS144, _M0L6_2atmpS1330);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1329
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS144);
    moonbit_decref(_M0L18_2astring__builderS144);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3186
    = _M0FPC15abort5abortGRPB18UninitializedArrayGiEE(_M0L6_2atmpS1329);
    moonbit_decref(_M0L6_2atmpS1329);
    return _result_3186;
  }
}

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t* _M0L3srcS149,
  int32_t _M0L13allocate__lenS145,
  int32_t _M0L3lenS146,
  int32_t _M0L11src__offsetS147,
  int32_t _M0L11dst__offsetS148
) {
  int32_t _if__result_3187;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS145 >= 0) {
    if (_M0L3lenS146 >= 0) {
      if (_M0L11src__offsetS147 >= 0) {
        if (_M0L11dst__offsetS148 >= 0) {
          int32_t _M0L6_2atmpS1332 = _M0L11src__offsetS147 + _M0L3lenS146;
          int32_t _M0L6_2atmpS1333;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1333
          = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS149);
          if (_M0L6_2atmpS1332 <= _M0L6_2atmpS1333) {
            int32_t _M0L6_2atmpS1331 = _M0L11dst__offsetS148 + _M0L3lenS146;
            _if__result_3187 = _M0L6_2atmpS1331 <= _M0L13allocate__lenS145;
          } else {
            _if__result_3187 = 0;
          }
        } else {
          _if__result_3187 = 0;
        }
      } else {
        _if__result_3187 = 0;
      }
    } else {
      _if__result_3187 = 0;
    }
  } else {
    _if__result_3187 = 0;
  }
  if (_if__result_3187) {
    moonbit_incref(_M0L3srcS149);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (moonbit_string_t*)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS145, (moonbit_string_t)moonbit_string_literal_0.data, _M0L3srcS149, _M0L11src__offsetS147, _M0L11dst__offsetS148, _M0L3lenS146);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS150;
    int32_t _M0L6_2atmpS1335;
    moonbit_string_t _M0L6_2atmpS1334;
    moonbit_string_t* _result_3188;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS150
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_93.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L13allocate__lenS145);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_94.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L11src__offsetS147);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_95.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L11dst__offsetS148);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_96.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L3lenS146);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_97.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1335 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS149);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L6_2atmpS1335);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1334
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS150);
    moonbit_decref(_M0L18_2astring__builderS150);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3188
    = _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(_M0L6_2atmpS1334);
    moonbit_decref(_M0L6_2atmpS1334);
    return _result_3188;
  }
}

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE** _M0L3srcS155,
  int32_t _M0L13allocate__lenS151,
  int32_t _M0L3lenS152,
  int32_t _M0L11src__offsetS153,
  int32_t _M0L11dst__offsetS154
) {
  int32_t _if__result_3189;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS151 >= 0) {
    if (_M0L3lenS152 >= 0) {
      if (_M0L11src__offsetS153 >= 0) {
        if (_M0L11dst__offsetS154 >= 0) {
          int32_t _M0L6_2atmpS1337 = _M0L11src__offsetS153 + _M0L3lenS152;
          int32_t _M0L6_2atmpS1338;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1338
          = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS155);
          if (_M0L6_2atmpS1337 <= _M0L6_2atmpS1338) {
            int32_t _M0L6_2atmpS1336 = _M0L11dst__offsetS154 + _M0L3lenS152;
            _if__result_3189 = _M0L6_2atmpS1336 <= _M0L13allocate__lenS151;
          } else {
            _if__result_3189 = 0;
          }
        } else {
          _if__result_3189 = 0;
        }
      } else {
        _if__result_3189 = 0;
      }
    } else {
      _if__result_3189 = 0;
    }
  } else {
    _if__result_3189 = 0;
  }
  if (_if__result_3189) {
    moonbit_incref(_M0L3srcS155);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TUsiE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS151, 0, _M0L3srcS155, _M0L11src__offsetS153, _M0L11dst__offsetS154, _M0L3lenS152);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS156;
    int32_t _M0L6_2atmpS1340;
    moonbit_string_t _M0L6_2atmpS1339;
    struct _M0TUsiE** _result_3190;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS156
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_93.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L13allocate__lenS151);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_94.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L11src__offsetS153);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_95.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L11dst__offsetS154);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_96.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L3lenS152);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_97.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1340 = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS155);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L6_2atmpS1340);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1339
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS156);
    moonbit_decref(_M0L18_2astring__builderS156);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3190
    = _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(_M0L6_2atmpS1339);
    moonbit_decref(_M0L6_2atmpS1339);
    return _result_3190;
  }
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS161,
  int32_t _M0L13allocate__lenS157,
  int32_t _M0L3lenS158,
  int32_t _M0L11src__offsetS159,
  int32_t _M0L11dst__offsetS160
) {
  int32_t _if__result_3191;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS157 >= 0) {
    if (_M0L3lenS158 >= 0) {
      if (_M0L11src__offsetS159 >= 0) {
        if (_M0L11dst__offsetS160 >= 0) {
          int32_t _M0L6_2atmpS1342 = _M0L11src__offsetS159 + _M0L3lenS158;
          int32_t _M0L6_2atmpS1343;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1343
          = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L3srcS161);
          if (_M0L6_2atmpS1342 <= _M0L6_2atmpS1343) {
            int32_t _M0L6_2atmpS1341 = _M0L11dst__offsetS160 + _M0L3lenS158;
            _if__result_3191 = _M0L6_2atmpS1341 <= _M0L13allocate__lenS157;
          } else {
            _if__result_3191 = 0;
          }
        } else {
          _if__result_3191 = 0;
        }
      } else {
        _if__result_3191 = 0;
      }
    } else {
      _if__result_3191 = 0;
    }
  } else {
    _if__result_3191 = 0;
  }
  if (_if__result_3191) {
    moonbit_incref(_M0L3srcS161);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS157, 0, _M0L3srcS161, _M0L11src__offsetS159, _M0L11dst__offsetS160, _M0L3lenS158);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS162;
    int32_t _M0L6_2atmpS1345;
    moonbit_string_t _M0L6_2atmpS1344;
    struct _M0TP26biolab8bio__seq9SeqRecord** _result_3192;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS162
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_93.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L13allocate__lenS157);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_94.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L11src__offsetS159);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_95.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L11dst__offsetS160);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_96.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L3lenS158);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_97.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1345
    = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L3srcS161);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L6_2atmpS1345);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1344
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS162);
    moonbit_decref(_M0L18_2astring__builderS162);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3192
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(_M0L6_2atmpS1344);
    moonbit_decref(_M0L6_2atmpS1344);
    return _result_3192;
  }
}

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3srcS167,
  int32_t _M0L13allocate__lenS163,
  int32_t _M0L3lenS164,
  int32_t _M0L11src__offsetS165,
  int32_t _M0L11dst__offsetS166
) {
  int32_t _if__result_3193;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS163 >= 0) {
    if (_M0L3lenS164 >= 0) {
      if (_M0L11src__offsetS165 >= 0) {
        if (_M0L11dst__offsetS166 >= 0) {
          int32_t _M0L6_2atmpS1347 = _M0L11src__offsetS165 + _M0L3lenS164;
          int32_t _M0L6_2atmpS1348;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1348
          = _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(_M0L3srcS167);
          if (_M0L6_2atmpS1347 <= _M0L6_2atmpS1348) {
            int32_t _M0L6_2atmpS1346 = _M0L11dst__offsetS166 + _M0L3lenS164;
            _if__result_3193 = _M0L6_2atmpS1346 <= _M0L13allocate__lenS163;
          } else {
            _if__result_3193 = 0;
          }
        } else {
          _if__result_3193 = 0;
        }
      } else {
        _if__result_3193 = 0;
      }
    } else {
      _if__result_3193 = 0;
    }
  } else {
    _if__result_3193 = 0;
  }
  if (_if__result_3193) {
    moonbit_incref(_M0L3srcS167);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return _M0MPB18UninitializedArray23unsafe__make__and__blitGRPC16string10StringViewE(_M0L3srcS167, _M0L13allocate__lenS163, _M0L11src__offsetS165, _M0L11dst__offsetS166, _M0L3lenS164);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS168;
    int32_t _M0L6_2atmpS1350;
    moonbit_string_t _M0L6_2atmpS1349;
    struct _M0TPC16string10StringView* _result_3194;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS168
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_93.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L13allocate__lenS163);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_94.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L11src__offsetS165);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_95.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L11dst__offsetS166);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_96.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L3lenS164);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_97.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1350
    = _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(_M0L3srcS167);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L6_2atmpS1350);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1349
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS168);
    moonbit_decref(_M0L18_2astring__builderS168);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3194
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRPC16string10StringViewEE(_M0L6_2atmpS1349);
    moonbit_decref(_M0L6_2atmpS1349);
    return _result_3194;
  }
}

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder* _M0L4selfS136,
  moonbit_string_t _M0L3objS135
) {
  struct _M0TPB6Logger _M0L6_2atmpS1324;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS136);
  _M0L6_2atmpS1324
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS136
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS135, _M0L6_2atmpS1324);
  if (_M0L6_2atmpS1324.$1) {
    moonbit_decref(_M0L6_2atmpS1324.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder* _M0L4selfS138,
  int32_t _M0L3objS137
) {
  struct _M0TPB6Logger _M0L6_2atmpS1325;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS138);
  _M0L6_2atmpS1325
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS138
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L3objS137, _M0L6_2atmpS1325);
  if (_M0L6_2atmpS1325.$1) {
    moonbit_decref(_M0L6_2atmpS1325.$1);
  }
  return 0;
}

int32_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGiE(
  int32_t* _M0L3srcS108,
  int32_t _M0L13allocate__lenS106,
  int32_t _M0L11src__offsetS109,
  int32_t _M0L11dst__offsetS107,
  int32_t _M0L9blit__lenS110
) {
  int32_t* _M0L3dstS105;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS105
  = (int32_t*)moonbit_make_int32_array_raw(_M0L13allocate__lenS106);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGiE(_M0L3dstS105, _M0L11dst__offsetS107, _M0L3srcS108, _M0L11src__offsetS109, _M0L9blit__lenS110);
  moonbit_decref(_M0L3srcS108);
  return _M0L3dstS105;
}

moonbit_string_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGsE(
  moonbit_string_t* _M0L3srcS114,
  int32_t _M0L13allocate__lenS112,
  int32_t _M0L11src__offsetS115,
  int32_t _M0L11dst__offsetS113,
  int32_t _M0L9blit__lenS116
) {
  moonbit_string_t* _M0L3dstS111;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS111
  = (moonbit_string_t*)moonbit_make_ref_array(_M0L13allocate__lenS112, (moonbit_string_t)moonbit_string_literal_0.data);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGsE(_M0L3dstS111, _M0L11dst__offsetS113, _M0L3srcS114, _M0L11src__offsetS115, _M0L9blit__lenS116);
  moonbit_decref(_M0L3srcS114);
  return _M0L3dstS111;
}

struct _M0TUsiE** _M0MPB18UninitializedArray23unsafe__make__and__blitGUsiEE(
  struct _M0TUsiE** _M0L3srcS120,
  int32_t _M0L13allocate__lenS118,
  int32_t _M0L11src__offsetS121,
  int32_t _M0L11dst__offsetS119,
  int32_t _M0L9blit__lenS122
) {
  struct _M0TUsiE** _M0L3dstS117;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS117
  = (struct _M0TUsiE**)moonbit_make_ref_array(_M0L13allocate__lenS118, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGUsiEE(_M0L3dstS117, _M0L11dst__offsetS119, _M0L3srcS120, _M0L11src__offsetS121, _M0L9blit__lenS122);
  moonbit_decref(_M0L3srcS120);
  return _M0L3dstS117;
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS126,
  int32_t _M0L13allocate__lenS124,
  int32_t _M0L11src__offsetS127,
  int32_t _M0L11dst__offsetS125,
  int32_t _M0L9blit__lenS128
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS123;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS123
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array(_M0L13allocate__lenS124, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(_M0L3dstS123, _M0L11dst__offsetS125, _M0L3srcS126, _M0L11src__offsetS127, _M0L9blit__lenS128);
  moonbit_decref(_M0L3srcS126);
  return _M0L3dstS123;
}

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23unsafe__make__and__blitGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3srcS132,
  int32_t _M0L13allocate__lenS130,
  int32_t _M0L11src__offsetS133,
  int32_t _M0L11dst__offsetS131,
  int32_t _M0L9blit__lenS134
) {
  struct _M0TPC16string10StringView* _M0L3dstS129;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS129
  = (struct _M0TPC16string10StringView*)moonbit_make_ref_valtype_array(_M0L13allocate__lenS130, sizeof(struct _M0TPC16string10StringView), Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 81, 0), &(struct _M0TPC16string10StringView){.$0 = (moonbit_string_t)moonbit_string_literal_0.data, .$1 = 0, .$2 = 0});
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRPC16string10StringViewE(_M0L3dstS129, _M0L11dst__offsetS131, _M0L3srcS132, _M0L11src__offsetS133, _M0L9blit__lenS134);
  moonbit_decref(_M0L3srcS132);
  return _M0L3dstS129;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGiE(
  int32_t* _M0L3dstS80,
  int32_t _M0L11dst__offsetS81,
  int32_t* _M0L3srcS82,
  int32_t _M0L11src__offsetS83,
  int32_t _M0L3lenS84
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS82);
  moonbit_incref(_M0L3dstS80);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS80, _M0L11dst__offsetS81, _M0L3srcS82, _M0L11src__offsetS83, _M0L3lenS84, sizeof(int32_t));
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGsE(
  moonbit_string_t* _M0L3dstS85,
  int32_t _M0L11dst__offsetS86,
  moonbit_string_t* _M0L3srcS87,
  int32_t _M0L11src__offsetS88,
  int32_t _M0L3lenS89
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS87);
  moonbit_incref(_M0L3dstS85);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS85, _M0L11dst__offsetS86, _M0L3srcS87, _M0L11src__offsetS88, _M0L3lenS89);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGUsiEE(
  struct _M0TUsiE** _M0L3dstS90,
  int32_t _M0L11dst__offsetS91,
  struct _M0TUsiE** _M0L3srcS92,
  int32_t _M0L11src__offsetS93,
  int32_t _M0L3lenS94
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS92);
  moonbit_incref(_M0L3dstS90);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS90, _M0L11dst__offsetS91, _M0L3srcS92, _M0L11src__offsetS93, _M0L3lenS94);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS95,
  int32_t _M0L11dst__offsetS96,
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS97,
  int32_t _M0L11src__offsetS98,
  int32_t _M0L3lenS99
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS97);
  moonbit_incref(_M0L3dstS95);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS95, _M0L11dst__offsetS96, _M0L3srcS97, _M0L11src__offsetS98, _M0L3lenS99);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3dstS100,
  int32_t _M0L11dst__offsetS101,
  struct _M0TPC16string10StringView* _M0L3srcS102,
  int32_t _M0L11src__offsetS103,
  int32_t _M0L3lenS104
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS102);
  moonbit_incref(_M0L3dstS100);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPC16string10StringViewEE(_M0L3dstS100, _M0L11dst__offsetS101, _M0L3srcS102, _M0L11src__offsetS103, _M0L3lenS104);
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGkE(
  uint16_t* _M0L3dstS26,
  int32_t _M0L11dst__offsetS28,
  uint16_t* _M0L3srcS27,
  int32_t _M0L11src__offsetS29,
  int32_t _M0L3lenS31
) {
  int32_t _if__result_3195;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS26 == _M0L3srcS27) {
    _if__result_3195 = _M0L11dst__offsetS28 < _M0L11src__offsetS29;
  } else {
    _if__result_3195 = 0;
  }
  if (_if__result_3195) {
    int32_t _M0L1iS30 = 0;
    while (1) {
      if (_M0L1iS30 < _M0L3lenS31) {
        int32_t _M0L6_2atmpS1270 = _M0L11dst__offsetS28 + _M0L1iS30;
        int32_t _M0L6_2atmpS1272 = _M0L11src__offsetS29 + _M0L1iS30;
        int32_t _M0L6_2atmpS1271;
        int32_t _M0L6_2atmpS1273;
        if (
          _M0L6_2atmpS1272 < 0
          || _M0L6_2atmpS1272 >= Moonbit_array_length(_M0L3srcS27)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1271 = (int32_t)_M0L3srcS27[_M0L6_2atmpS1272];
        if (
          _M0L6_2atmpS1270 < 0
          || _M0L6_2atmpS1270 >= Moonbit_array_length(_M0L3dstS26)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS26[_M0L6_2atmpS1270] = _M0L6_2atmpS1271;
        _M0L6_2atmpS1273 = _M0L1iS30 + 1;
        _M0L1iS30 = _M0L6_2atmpS1273;
        continue;
      } else {
        moonbit_decref(_M0L3srcS27);
        moonbit_decref(_M0L3dstS26);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1278 = _M0L3lenS31 - 1;
    int32_t _M0L1iS33 = _M0L6_2atmpS1278;
    while (1) {
      if (_M0L1iS33 >= 0) {
        int32_t _M0L6_2atmpS1274 = _M0L11dst__offsetS28 + _M0L1iS33;
        int32_t _M0L6_2atmpS1276 = _M0L11src__offsetS29 + _M0L1iS33;
        int32_t _M0L6_2atmpS1275;
        int32_t _M0L6_2atmpS1277;
        if (
          _M0L6_2atmpS1276 < 0
          || _M0L6_2atmpS1276 >= Moonbit_array_length(_M0L3srcS27)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1275 = (int32_t)_M0L3srcS27[_M0L6_2atmpS1276];
        if (
          _M0L6_2atmpS1274 < 0
          || _M0L6_2atmpS1274 >= Moonbit_array_length(_M0L3dstS26)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS26[_M0L6_2atmpS1274] = _M0L6_2atmpS1275;
        _M0L6_2atmpS1277 = _M0L1iS33 - 1;
        _M0L1iS33 = _M0L6_2atmpS1277;
        continue;
      } else {
        moonbit_decref(_M0L3srcS27);
        moonbit_decref(_M0L3dstS26);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGiEE(
  int32_t* _M0L3dstS35,
  int32_t _M0L11dst__offsetS37,
  int32_t* _M0L3srcS36,
  int32_t _M0L11src__offsetS38,
  int32_t _M0L3lenS40
) {
  int32_t _if__result_3198;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS35 == _M0L3srcS36) {
    _if__result_3198 = _M0L11dst__offsetS37 < _M0L11src__offsetS38;
  } else {
    _if__result_3198 = 0;
  }
  if (_if__result_3198) {
    int32_t _M0L1iS39 = 0;
    while (1) {
      if (_M0L1iS39 < _M0L3lenS40) {
        int32_t _M0L6_2atmpS1279 = _M0L11dst__offsetS37 + _M0L1iS39;
        int32_t _M0L6_2atmpS1281 = _M0L11src__offsetS38 + _M0L1iS39;
        int32_t _M0L6_2atmpS1280;
        int32_t _M0L6_2atmpS1282;
        if (
          _M0L6_2atmpS1281 < 0
          || _M0L6_2atmpS1281 >= Moonbit_array_length(_M0L3srcS36)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1280 = (int32_t)_M0L3srcS36[_M0L6_2atmpS1281];
        if (
          _M0L6_2atmpS1279 < 0
          || _M0L6_2atmpS1279 >= Moonbit_array_length(_M0L3dstS35)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS35[_M0L6_2atmpS1279] = _M0L6_2atmpS1280;
        _M0L6_2atmpS1282 = _M0L1iS39 + 1;
        _M0L1iS39 = _M0L6_2atmpS1282;
        continue;
      } else {
        moonbit_decref(_M0L3srcS36);
        moonbit_decref(_M0L3dstS35);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1287 = _M0L3lenS40 - 1;
    int32_t _M0L1iS42 = _M0L6_2atmpS1287;
    while (1) {
      if (_M0L1iS42 >= 0) {
        int32_t _M0L6_2atmpS1283 = _M0L11dst__offsetS37 + _M0L1iS42;
        int32_t _M0L6_2atmpS1285 = _M0L11src__offsetS38 + _M0L1iS42;
        int32_t _M0L6_2atmpS1284;
        int32_t _M0L6_2atmpS1286;
        if (
          _M0L6_2atmpS1285 < 0
          || _M0L6_2atmpS1285 >= Moonbit_array_length(_M0L3srcS36)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1284 = (int32_t)_M0L3srcS36[_M0L6_2atmpS1285];
        if (
          _M0L6_2atmpS1283 < 0
          || _M0L6_2atmpS1283 >= Moonbit_array_length(_M0L3dstS35)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS35[_M0L6_2atmpS1283] = _M0L6_2atmpS1284;
        _M0L6_2atmpS1286 = _M0L1iS42 - 1;
        _M0L1iS42 = _M0L6_2atmpS1286;
        continue;
      } else {
        moonbit_decref(_M0L3srcS36);
        moonbit_decref(_M0L3dstS35);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGsEE(
  moonbit_string_t* _M0L3dstS44,
  int32_t _M0L11dst__offsetS46,
  moonbit_string_t* _M0L3srcS45,
  int32_t _M0L11src__offsetS47,
  int32_t _M0L3lenS49
) {
  int32_t _if__result_3201;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS44 == _M0L3srcS45) {
    _if__result_3201 = _M0L11dst__offsetS46 < _M0L11src__offsetS47;
  } else {
    _if__result_3201 = 0;
  }
  if (_if__result_3201) {
    int32_t _M0L1iS48 = 0;
    while (1) {
      if (_M0L1iS48 < _M0L3lenS49) {
        int32_t _M0L6_2atmpS1288 = _M0L11dst__offsetS46 + _M0L1iS48;
        int32_t _M0L6_2atmpS1290 = _M0L11src__offsetS47 + _M0L1iS48;
        moonbit_string_t _M0L6_2atmpS1289;
        moonbit_string_t _M0L6_2aoldS2874;
        int32_t _M0L6_2atmpS1291;
        if (
          _M0L6_2atmpS1290 < 0
          || _M0L6_2atmpS1290 >= Moonbit_array_length(_M0L3srcS45)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1289 = (moonbit_string_t)_M0L3srcS45[_M0L6_2atmpS1290];
        if (
          _M0L6_2atmpS1288 < 0
          || _M0L6_2atmpS1288 >= Moonbit_array_length(_M0L3dstS44)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2874 = (moonbit_string_t)_M0L3dstS44[_M0L6_2atmpS1288];
        moonbit_incref(_M0L6_2atmpS1289);
        moonbit_decref(_M0L6_2aoldS2874);
        _M0L3dstS44[_M0L6_2atmpS1288] = _M0L6_2atmpS1289;
        _M0L6_2atmpS1291 = _M0L1iS48 + 1;
        _M0L1iS48 = _M0L6_2atmpS1291;
        continue;
      } else {
        moonbit_decref(_M0L3srcS45);
        moonbit_decref(_M0L3dstS44);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1296 = _M0L3lenS49 - 1;
    int32_t _M0L1iS51 = _M0L6_2atmpS1296;
    while (1) {
      if (_M0L1iS51 >= 0) {
        int32_t _M0L6_2atmpS1292 = _M0L11dst__offsetS46 + _M0L1iS51;
        int32_t _M0L6_2atmpS1294 = _M0L11src__offsetS47 + _M0L1iS51;
        moonbit_string_t _M0L6_2atmpS1293;
        moonbit_string_t _M0L6_2aoldS2876;
        int32_t _M0L6_2atmpS1295;
        if (
          _M0L6_2atmpS1294 < 0
          || _M0L6_2atmpS1294 >= Moonbit_array_length(_M0L3srcS45)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1293 = (moonbit_string_t)_M0L3srcS45[_M0L6_2atmpS1294];
        if (
          _M0L6_2atmpS1292 < 0
          || _M0L6_2atmpS1292 >= Moonbit_array_length(_M0L3dstS44)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2876 = (moonbit_string_t)_M0L3dstS44[_M0L6_2atmpS1292];
        moonbit_incref(_M0L6_2atmpS1293);
        moonbit_decref(_M0L6_2aoldS2876);
        _M0L3dstS44[_M0L6_2atmpS1292] = _M0L6_2atmpS1293;
        _M0L6_2atmpS1295 = _M0L1iS51 - 1;
        _M0L1iS51 = _M0L6_2atmpS1295;
        continue;
      } else {
        moonbit_decref(_M0L3srcS45);
        moonbit_decref(_M0L3dstS44);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGUsiEEE(
  struct _M0TUsiE** _M0L3dstS53,
  int32_t _M0L11dst__offsetS55,
  struct _M0TUsiE** _M0L3srcS54,
  int32_t _M0L11src__offsetS56,
  int32_t _M0L3lenS58
) {
  int32_t _if__result_3204;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS53 == _M0L3srcS54) {
    _if__result_3204 = _M0L11dst__offsetS55 < _M0L11src__offsetS56;
  } else {
    _if__result_3204 = 0;
  }
  if (_if__result_3204) {
    int32_t _M0L1iS57 = 0;
    while (1) {
      if (_M0L1iS57 < _M0L3lenS58) {
        int32_t _M0L6_2atmpS1297 = _M0L11dst__offsetS55 + _M0L1iS57;
        int32_t _M0L6_2atmpS1299 = _M0L11src__offsetS56 + _M0L1iS57;
        struct _M0TUsiE* _M0L6_2atmpS1298;
        struct _M0TUsiE* _M0L6_2aoldS2878;
        int32_t _M0L6_2atmpS1300;
        if (
          _M0L6_2atmpS1299 < 0
          || _M0L6_2atmpS1299 >= Moonbit_array_length(_M0L3srcS54)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1298 = (struct _M0TUsiE*)_M0L3srcS54[_M0L6_2atmpS1299];
        if (
          _M0L6_2atmpS1297 < 0
          || _M0L6_2atmpS1297 >= Moonbit_array_length(_M0L3dstS53)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2878 = (struct _M0TUsiE*)_M0L3dstS53[_M0L6_2atmpS1297];
        if (_M0L6_2atmpS1298) {
          moonbit_incref(_M0L6_2atmpS1298);
        }
        if (_M0L6_2aoldS2878) {
          moonbit_decref(_M0L6_2aoldS2878);
        }
        _M0L3dstS53[_M0L6_2atmpS1297] = _M0L6_2atmpS1298;
        _M0L6_2atmpS1300 = _M0L1iS57 + 1;
        _M0L1iS57 = _M0L6_2atmpS1300;
        continue;
      } else {
        moonbit_decref(_M0L3srcS54);
        moonbit_decref(_M0L3dstS53);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1305 = _M0L3lenS58 - 1;
    int32_t _M0L1iS60 = _M0L6_2atmpS1305;
    while (1) {
      if (_M0L1iS60 >= 0) {
        int32_t _M0L6_2atmpS1301 = _M0L11dst__offsetS55 + _M0L1iS60;
        int32_t _M0L6_2atmpS1303 = _M0L11src__offsetS56 + _M0L1iS60;
        struct _M0TUsiE* _M0L6_2atmpS1302;
        struct _M0TUsiE* _M0L6_2aoldS2880;
        int32_t _M0L6_2atmpS1304;
        if (
          _M0L6_2atmpS1303 < 0
          || _M0L6_2atmpS1303 >= Moonbit_array_length(_M0L3srcS54)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1302 = (struct _M0TUsiE*)_M0L3srcS54[_M0L6_2atmpS1303];
        if (
          _M0L6_2atmpS1301 < 0
          || _M0L6_2atmpS1301 >= Moonbit_array_length(_M0L3dstS53)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2880 = (struct _M0TUsiE*)_M0L3dstS53[_M0L6_2atmpS1301];
        if (_M0L6_2atmpS1302) {
          moonbit_incref(_M0L6_2atmpS1302);
        }
        if (_M0L6_2aoldS2880) {
          moonbit_decref(_M0L6_2aoldS2880);
        }
        _M0L3dstS53[_M0L6_2atmpS1301] = _M0L6_2atmpS1302;
        _M0L6_2atmpS1304 = _M0L1iS60 - 1;
        _M0L1iS60 = _M0L6_2atmpS1304;
        continue;
      } else {
        moonbit_decref(_M0L3srcS54);
        moonbit_decref(_M0L3dstS53);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS62,
  int32_t _M0L11dst__offsetS64,
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS63,
  int32_t _M0L11src__offsetS65,
  int32_t _M0L3lenS67
) {
  int32_t _if__result_3207;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS62 == _M0L3srcS63) {
    _if__result_3207 = _M0L11dst__offsetS64 < _M0L11src__offsetS65;
  } else {
    _if__result_3207 = 0;
  }
  if (_if__result_3207) {
    int32_t _M0L1iS66 = 0;
    while (1) {
      if (_M0L1iS66 < _M0L3lenS67) {
        int32_t _M0L6_2atmpS1306 = _M0L11dst__offsetS64 + _M0L1iS66;
        int32_t _M0L6_2atmpS1308 = _M0L11src__offsetS65 + _M0L1iS66;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1307;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS2882;
        int32_t _M0L6_2atmpS1309;
        if (
          _M0L6_2atmpS1308 < 0
          || _M0L6_2atmpS1308 >= Moonbit_array_length(_M0L3srcS63)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1307
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3srcS63[
            _M0L6_2atmpS1308
          ];
        if (
          _M0L6_2atmpS1306 < 0
          || _M0L6_2atmpS1306 >= Moonbit_array_length(_M0L3dstS62)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2882
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3dstS62[
            _M0L6_2atmpS1306
          ];
        if (_M0L6_2atmpS1307) {
          moonbit_incref(_M0L6_2atmpS1307);
        }
        if (_M0L6_2aoldS2882) {
          moonbit_decref(_M0L6_2aoldS2882);
        }
        _M0L3dstS62[_M0L6_2atmpS1306] = _M0L6_2atmpS1307;
        _M0L6_2atmpS1309 = _M0L1iS66 + 1;
        _M0L1iS66 = _M0L6_2atmpS1309;
        continue;
      } else {
        moonbit_decref(_M0L3srcS63);
        moonbit_decref(_M0L3dstS62);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1314 = _M0L3lenS67 - 1;
    int32_t _M0L1iS69 = _M0L6_2atmpS1314;
    while (1) {
      if (_M0L1iS69 >= 0) {
        int32_t _M0L6_2atmpS1310 = _M0L11dst__offsetS64 + _M0L1iS69;
        int32_t _M0L6_2atmpS1312 = _M0L11src__offsetS65 + _M0L1iS69;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1311;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS2884;
        int32_t _M0L6_2atmpS1313;
        if (
          _M0L6_2atmpS1312 < 0
          || _M0L6_2atmpS1312 >= Moonbit_array_length(_M0L3srcS63)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1311
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3srcS63[
            _M0L6_2atmpS1312
          ];
        if (
          _M0L6_2atmpS1310 < 0
          || _M0L6_2atmpS1310 >= Moonbit_array_length(_M0L3dstS62)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2884
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3dstS62[
            _M0L6_2atmpS1310
          ];
        if (_M0L6_2atmpS1311) {
          moonbit_incref(_M0L6_2atmpS1311);
        }
        if (_M0L6_2aoldS2884) {
          moonbit_decref(_M0L6_2aoldS2884);
        }
        _M0L3dstS62[_M0L6_2atmpS1310] = _M0L6_2atmpS1311;
        _M0L6_2atmpS1313 = _M0L1iS69 - 1;
        _M0L1iS69 = _M0L6_2atmpS1313;
        continue;
      } else {
        moonbit_decref(_M0L3srcS63);
        moonbit_decref(_M0L3dstS62);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPC16string10StringViewEE(
  struct _M0TPC16string10StringView* _M0L3dstS71,
  int32_t _M0L11dst__offsetS73,
  struct _M0TPC16string10StringView* _M0L3srcS72,
  int32_t _M0L11src__offsetS74,
  int32_t _M0L3lenS76
) {
  int32_t _if__result_3210;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS71 == _M0L3srcS72) {
    _if__result_3210 = _M0L11dst__offsetS73 < _M0L11src__offsetS74;
  } else {
    _if__result_3210 = 0;
  }
  if (_if__result_3210) {
    int32_t _M0L1iS75 = 0;
    while (1) {
      if (_M0L1iS75 < _M0L3lenS76) {
        int32_t _M0L6_2atmpS1315 = _M0L11dst__offsetS73 + _M0L1iS75;
        int32_t _M0L6_2atmpS1317 = _M0L11src__offsetS74 + _M0L1iS75;
        struct _M0TPC16string10StringView _M0L6_2atmpS1316;
        struct _M0TPC16string10StringView _M0L6_2aoldS2886;
        int32_t _M0L6_2atmpS1318;
        if (
          _M0L6_2atmpS1317 < 0
          || _M0L6_2atmpS1317 >= Moonbit_array_length(_M0L3srcS72)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1316 = _M0L3srcS72[_M0L6_2atmpS1317];
        if (
          _M0L6_2atmpS1315 < 0
          || _M0L6_2atmpS1315 >= Moonbit_array_length(_M0L3dstS71)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2886 = _M0L3dstS71[_M0L6_2atmpS1315];
        moonbit_incref(_M0L6_2atmpS1316.$0);
        moonbit_decref(_M0L6_2aoldS2886.$0);
        _M0L3dstS71[_M0L6_2atmpS1315] = _M0L6_2atmpS1316;
        _M0L6_2atmpS1318 = _M0L1iS75 + 1;
        _M0L1iS75 = _M0L6_2atmpS1318;
        continue;
      } else {
        moonbit_decref(_M0L3srcS72);
        moonbit_decref(_M0L3dstS71);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1323 = _M0L3lenS76 - 1;
    int32_t _M0L1iS78 = _M0L6_2atmpS1323;
    while (1) {
      if (_M0L1iS78 >= 0) {
        int32_t _M0L6_2atmpS1319 = _M0L11dst__offsetS73 + _M0L1iS78;
        int32_t _M0L6_2atmpS1321 = _M0L11src__offsetS74 + _M0L1iS78;
        struct _M0TPC16string10StringView _M0L6_2atmpS1320;
        struct _M0TPC16string10StringView _M0L6_2aoldS2888;
        int32_t _M0L6_2atmpS1322;
        if (
          _M0L6_2atmpS1321 < 0
          || _M0L6_2atmpS1321 >= Moonbit_array_length(_M0L3srcS72)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1320 = _M0L3srcS72[_M0L6_2atmpS1321];
        if (
          _M0L6_2atmpS1319 < 0
          || _M0L6_2atmpS1319 >= Moonbit_array_length(_M0L3dstS71)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS2888 = _M0L3dstS71[_M0L6_2atmpS1319];
        moonbit_incref(_M0L6_2atmpS1320.$0);
        moonbit_decref(_M0L6_2aoldS2888.$0);
        _M0L3dstS71[_M0L6_2atmpS1319] = _M0L6_2atmpS1320;
        _M0L6_2atmpS1322 = _M0L1iS78 - 1;
        _M0L1iS78 = _M0L6_2atmpS1322;
        continue;
      } else {
        moonbit_decref(_M0L3srcS72);
        moonbit_decref(_M0L3dstS71);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPB18UninitializedArray6lengthGiE(int32_t* _M0L4selfS21) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS21);
}

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t* _M0L4selfS22) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS22);
}

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(
  struct _M0TUsiE** _M0L4selfS23
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS23);
}

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L4selfS24
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS24);
}

int32_t _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L4selfS25
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS25);
}

uint32_t _M0FPB13consume4__acc(uint32_t _M0L3accS19, uint32_t _M0L5inputS20) {
  uint32_t _M0L6_2atmpS1269;
  uint32_t _M0L6_2atmpS1268;
  uint32_t _M0L6_2atmpS1267;
  #line 465 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1269 = _M0L5inputS20 * 3266489917u;
  _M0L6_2atmpS1268 = _M0L3accS19 + _M0L6_2atmpS1269;
  #line 466 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1267 = _M0FPB4rotl(_M0L6_2atmpS1268, 17);
  return _M0L6_2atmpS1267 * 668265263u;
}

uint32_t _M0FPB4rotl(uint32_t _M0L1xS17, int32_t _M0L1rS18) {
  uint32_t _M0L6_2atmpS1264;
  int32_t _M0L6_2atmpS1266;
  uint32_t _M0L6_2atmpS1265;
  #line 475 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1264 = _M0L1xS17 << (_M0L1rS18 & 31);
  _M0L6_2atmpS1266 = 32 - _M0L1rS18;
  _M0L6_2atmpS1265 = _M0L1xS17 >> (_M0L6_2atmpS1266 & 31);
  return _M0L6_2atmpS1264 | _M0L6_2atmpS1265;
}

int32_t _M0IPB7FailurePB4Show6output(
  void* _M0L10_2ax__5721S13,
  struct _M0TPB6Logger _M0L10_2ax__5722S16
) {
  struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS14;
  moonbit_string_t _M0L15_2a_2aarg__5723S15;
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2aFailureS14
  = (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L10_2ax__5721S13;
  _M0L15_2a_2aarg__5723S15 = _M0L10_2aFailureS14->$0;
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S16.$0->$method_0(_M0L10_2ax__5722S16.$1, (moonbit_string_t)moonbit_string_literal_98.data);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB6Logger13write__objectGsE(_M0L10_2ax__5722S16, _M0L15_2a_2aarg__5723S15);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S16.$0->$method_0(_M0L10_2ax__5722S16.$1, (moonbit_string_t)moonbit_string_literal_99.data);
  return 0;
}

int32_t _M0MPB6Logger13write__objectGsE(
  struct _M0TPB6Logger _M0L4selfS12,
  moonbit_string_t _M0L3objS11
) {
  #line 173 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 174 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS11, _M0L4selfS12);
  return 0;
}

moonbit_string_t _M0FPC15abort5abortGsE(moonbit_string_t _M0L3msgS1) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS1);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t _M0FPC15abort5abortGuE(moonbit_string_t _M0L3msgS2) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS2);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
  return 0;
}

uint16_t* _M0FPC15abort5abortGAkE(moonbit_string_t _M0L3msgS3) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS3);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t _M0FPC15abort5abortGiE(moonbit_string_t _M0L3msgS4) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS4);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t* _M0FPC15abort5abortGRPB18UninitializedArrayGiEE(
  moonbit_string_t _M0L3msgS5
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS5);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

moonbit_string_t* _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(
  moonbit_string_t _M0L3msgS6
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS6);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TUsiE** _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(
  moonbit_string_t _M0L3msgS7
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS7);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TPC16string10StringView _M0FPC15abort5abortGRPC16string10StringViewE(
  moonbit_string_t _M0L3msgS8
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS8);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(
  moonbit_string_t _M0L3msgS9
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS9);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TPC16string10StringView* _M0FPC15abort5abortGRPB18UninitializedArrayGRPC16string10StringViewEE(
  moonbit_string_t _M0L3msgS10
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS10);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

moonbit_string_t _M0FP15Error10to__string(void* _M0L4_2aeS1230) {
  switch (Moonbit_object_tag(_M0L4_2aeS1230)) {
    case 3: {
      return (moonbit_string_t)moonbit_string_literal_100.data;
      break;
    }
    
    case 2: {
      return (moonbit_string_t)moonbit_string_literal_101.data;
      break;
    }
    
    case 1: {
      return (moonbit_string_t)moonbit_string_literal_102.data;
      break;
    }
    
    case 0: {
      return _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(_M0L4_2aeS1230);
      break;
    }
    
    case 4: {
      return (moonbit_string_t)moonbit_string_literal_103.data;
      break;
    }
    default: {
      return (moonbit_string_t)moonbit_string_literal_104.data;
      break;
    }
  }
}

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1253,
  struct _M0TPB4Show _M0L8_2aparamS1252
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1251 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1253;
  _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(_M0L7_2aselfS1251, _M0L8_2aparamS1252);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1250,
  struct _M0TPB4Show _M0L8_2aparamS1249
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1248 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1250;
  _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(_M0L7_2aselfS1248, _M0L8_2aparamS1249);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1247,
  int32_t _M0L8_2aparamS1246
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1245 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1247;
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS1245, _M0L8_2aparamS1246);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1244,
  struct _M0TPC16string10StringView _M0L8_2aparamS1243
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1242 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1244;
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L7_2aselfS1242, _M0L8_2aparamS1243);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1241,
  moonbit_string_t _M0L8_2aparamS1238,
  int32_t _M0L8_2aparamS1239,
  int32_t _M0L8_2aparamS1240
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1237 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1241;
  _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L7_2aselfS1237, _M0L8_2aparamS1238, _M0L8_2aparamS1239, _M0L8_2aparamS1240);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1236,
  moonbit_string_t _M0L8_2aparamS1235
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1234 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1236;
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L7_2aselfS1234, _M0L8_2aparamS1235);
  return 0;
}

void moonbit_init() {
  moonbit_layout_table = moonbit_layout_table_data;
  uint16_t* _M0L1tS908 = (uint16_t*)moonbit_make_string(128, 0);
  int32_t _M0L6_2atmpS1263 = 97;
  int32_t _M0L1cS909 = _M0L6_2atmpS1263;
  while (1) {
    int32_t _M0L6_2atmpS1259 = 122;
    if (_M0L1cS909 <= _M0L6_2atmpS1259) {
      int32_t _M0L6_2atmpS1261 = _M0L1cS909 - 32;
      int32_t _M0L6_2atmpS1260 = (uint16_t)_M0L6_2atmpS1261;
      int32_t _M0L6_2atmpS1262;
      if (_M0L1cS909 < 0 || _M0L1cS909 >= Moonbit_array_length(_M0L1tS908)) {
        #line 35 "/home/developer/Documents/1/genbank_io.mbt"
        moonbit_panic();
      }
      _M0L1tS908[_M0L1cS909] = _M0L6_2atmpS1260;
      _M0L6_2atmpS1262 = _M0L1cS909 + 1;
      _M0L1cS909 = _M0L6_2atmpS1262;
      continue;
    }
    break;
  }
  _M0FP26biolab8bio__seq16uppercase__table = _M0L1tS908;
}

int main(int argc, char** argv) {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** _M0L6_2atmpS1258;
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1224;
  struct _M0TPB5ArrayGUsiEE* _M0L7_2abindS1225;
  int32_t _M0L7_2abindS1226;
  int32_t _M0L2__S1227;
  moonbit_runtime_init(argc, argv);
  moonbit_init();
  _M0L6_2atmpS1258
  = (struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error**)moonbit_empty_ref_array;
  _M0L12async__testsS1224
  = (struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE));
  Moonbit_object_header(_M0L12async__testsS1224)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 127, 0);
  _M0L12async__testsS1224->$0 = _M0L6_2atmpS1258;
  _M0L12async__testsS1224->$1 = 0;
  #line 399 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L7_2abindS1225
  = _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test52moonbit__test__driver__internal__native__parse__args();
  _M0L7_2abindS1226 = _M0L7_2abindS1225->$1;
  _M0L2__S1227 = 0;
  while (1) {
    if (_M0L2__S1227 < _M0L7_2abindS1226) {
      struct _M0TUsiE** _M0L3bufS1257 = _M0L7_2abindS1225->$0;
      struct _M0TUsiE* _M0L3argS1228 =
        (struct _M0TUsiE*)_M0L3bufS1257[_M0L2__S1227];
      moonbit_string_t _M0L6_2atmpS1254 = _M0L3argS1228->$0;
      int32_t _M0L6_2atmpS1255 = _M0L3argS1228->$1;
      int32_t _M0L6_2atmpS1256;
      moonbit_incref(_M0L6_2atmpS1254);
      #line 400 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
      _M0FP46biolab8bio__seq3cmd27seqio__main__blackbox__test44moonbit__test__driver__internal__do__execute(_M0L12async__testsS1224, _M0L6_2atmpS1254, _M0L6_2atmpS1255);
      moonbit_decref(_M0L6_2atmpS1254);
      _M0L6_2atmpS1256 = _M0L2__S1227 + 1;
      _M0L2__S1227 = _M0L6_2atmpS1256;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1225);
    }
    break;
  }
  #line 402 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IP016_24default__implP46biolab8bio__seq3cmd27seqio__main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd27seqio__main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(_M0L12async__testsS1224);
  moonbit_decref(_M0L12async__testsS1224);
  return 0;
}