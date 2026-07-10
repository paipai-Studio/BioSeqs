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
struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__;

struct _M0TPB5ArrayGRPC16string10StringViewE;

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure;

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError;

struct _M0TPB4IterGsE;

struct _M0TPB8MutLocalGORPC16string10StringViewE;

struct _M0TPB4IterGUsRPB5ArrayGiEEE;

struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__;

struct _M0TWRPC15error5ErrorEs;

struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest;

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE3Err;

struct _M0TWssbEu;

struct _M0TUsiE;

struct _M0TPB8MutLocalGsE;

struct _M0TPB13StringBuilder;

struct _M0TPB8MutLocalGRPB13StringBuilderE;

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some;

struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477;

struct _M0TWEOUsRPB5ArrayGiEE;

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__;

struct _M0TWEuQRPC15error5Error;

struct _M0TPB3MapGsRPB5ArrayGiEE;

struct _M0TP26biolab8bio__seq3Seq;

struct _M0BTPB6Logger;

struct _M0DTPC16result6ResultGbRP26biolab8bio__seq33MoonBitTestDriverInternalSkipTestE3Err;

struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE;

struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__;

struct _M0TPB4IterGcE;

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE;

struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError;

struct _M0TPB6Logger;

struct _M0TWcERPC16string10StringView;

struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE;

struct _M0TPB5ArrayGUsiEE;

struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE;

struct _M0TPB8MutLocalGOsE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok;

struct _M0TPB5ArrayGiE;

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq9SeqRecordRP26biolab8bio__seq10SeqIOErrorE2Ok;

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__;

struct _M0TPB8MutLocalGAkE;

struct _M0TUiUWEuQRPC15error5ErrorNsEE;

struct _M0DTPC16result6ResultGuRPB7FailureE3Err;

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq10SeqIOErrorE3Err;

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError;

struct _M0TWRPC15error5ErrorEu;

struct _M0TURPC16string10StringViewRPB6LoggerE;

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq10SeqIOErrorE2Ok;

struct _M0TPB4IterGRPC16string10StringViewE;

struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError;

struct _M0TPB5EntryGsRPB5ArrayGiEE;

struct _M0TPB8MutLocalGiE;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq9SeqRecordRP26biolab8bio__seq10SeqIOErrorE3Err;

struct _M0TWERPC16option6OptionGRPC16string10StringViewE;

struct _M0TPB4Show;

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE;

struct _M0TWEOc;

struct _M0TUsRPB5ArrayGiEE;

struct _M0TP26biolab8bio__seq9SeqRecord;

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE;

struct _M0DTPC16result6ResultGuRPC15error5ErrorE2Ok;

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error;

struct _M0TUsRP26biolab8bio__seq9SeqRecordE;

struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err;

struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE;

struct _M0BTPB4Show;

struct _M0TWuEu;

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE2Ok;

struct _M0TPC16string10StringView;

struct _M0TPB8MutLocalGbE;

struct _M0DTPC16result6ResultGuRPB7FailureE2Ok;

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__;

struct _M0KTPB6LoggerTPB13StringBuilder;

struct _M0DTPC16result6ResultGbRP26biolab8bio__seq33MoonBitTestDriverInternalSkipTestE2Ok;

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE;

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE;

struct _M0TPB5ArrayGsE;

struct _M0TUWEuQRPC15error5ErrorNsE;

struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__;

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE;

struct _M0TPB9ArrayViewGsE;

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE;

struct _M0TWEOs;

struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482;

struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE;

struct _M0DTPC16result6ResultGuRPC15error5ErrorE3Err;

struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__ {
  struct _M0TUsRPB5ArrayGiEE*(* code)(struct _M0TWEOUsRPB5ArrayGiEE*);
  struct _M0TPB8MutLocalGiE* $0;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE* $1;
  
};

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

struct _M0TPB4IterGsE {
  struct _M0TWEOs* $0;
  int64_t $1;
  
};

struct _M0TPB8MutLocalGORPC16string10StringViewE {
  void* $0;
  
};

struct _M0TPB4IterGUsRPB5ArrayGiEEE {
  struct _M0TWEOUsRPB5ArrayGiEE* $0;
  int64_t $1;
  
};

struct _M0TPC16string10StringView {
  moonbit_string_t $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  int32_t $1;
  struct _M0TPC16string10StringView $2;
  
};

struct _M0TWRPC15error5ErrorEs {
  moonbit_string_t(* code)(struct _M0TWRPC15error5ErrorEs*, void*);
  
};

struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest {
  moonbit_string_t $0;
  
};

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq9SeqRecordERP26biolab8bio__seq10SeqIOErrorE3Err {
  void* $0;
  
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

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some {
  struct _M0TPC16string10StringView $0;
  
};

struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477 {
  int32_t(* code)(struct _M0TWEOc*);
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TWEOUsRPB5ArrayGiEE {
  struct _M0TUsRPB5ArrayGiEE*(* code)(struct _M0TWEOUsRPB5ArrayGiEE*);
  
};

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TWssbEu* $0;
  moonbit_string_t $1;
  
};

struct _M0TWEuQRPC15error5Error {
  struct moonbit_result_0(* code)(struct _M0TWEuQRPC15error5Error*);
  
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

struct _M0DTPC16result6ResultGbRP26biolab8bio__seq33MoonBitTestDriverInternalSkipTestE3Err {
  void* $0;
  
};

struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE {
  struct _M0TUsRP26biolab8bio__seq9SeqRecordE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TPB8MutLocalGORPC16string10StringViewE* $0;
  struct _M0TPC16string10StringView $1;
  int32_t $2;
  
};

struct _M0TPB4IterGcE {
  struct _M0TWEOc* $0;
  int64_t $1;
  
};

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE {
  int32_t $0;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* $5;
  
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

struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE {
  moonbit_string_t $0;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* $1;
  
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

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TWcERPC16string10StringView* $0;
  struct _M0TPB4IterGcE* $1;
  
};

struct _M0TPB8MutLocalGAkE {
  uint16_t* $0;
  
};

struct _M0TUiUWEuQRPC15error5ErrorNsEE {
  int32_t $0;
  struct _M0TUWEuQRPC15error5ErrorNsE* $1;
  
};

struct _M0DTPC16result6ResultGuRPB7FailureE3Err {
  void* $0;
  
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

struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError {
  moonbit_string_t $0;
  
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

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE {
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* $5;
  int32_t $6;
  
};

struct _M0DTPC16result6ResultGuRPC15error5ErrorE2Ok {
  int32_t $0;
  
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

struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE {
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** $0;
  int32_t $1;
  int32_t $2;
  
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

struct _M0DTPC16result6ResultGuRPB7FailureE2Ok {
  int32_t $0;
  
};

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__ {
  int32_t(* code)(struct _M0TWRPC15error5ErrorEu*, void*);
  struct _M0TWRPC15error5ErrorEs* $0;
  struct _M0TWssbEu* $1;
  moonbit_string_t $2;
  
};

struct _M0KTPB6LoggerTPB13StringBuilder {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0DTPC16result6ResultGbRP26biolab8bio__seq33MoonBitTestDriverInternalSkipTestE2Ok {
  int32_t $0;
  
};

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE {
  int32_t $0;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TUWEuQRPC15error5ErrorNsE* $5;
  
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

struct _M0TUWEuQRPC15error5ErrorNsE {
  struct _M0TWEuQRPC15error5Error* $0;
  moonbit_string_t* $1;
  
};

struct _M0TPB9ArrayViewGsE {
  moonbit_string_t* $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__ {
  moonbit_string_t(* code)(struct _M0TWEOs*);
  struct _M0TPB9ArrayViewGsE $0;
  int32_t $1;
  struct _M0TPB8MutLocalGiE* $2;
  
};

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE {
  struct _M0TP26biolab8bio__seq9SeqRecord** $0;
  int32_t $1;
  
};

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE {
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* $5;
  int32_t $6;
  
};

struct _M0TWEOs {
  moonbit_string_t(* code)(struct _M0TWEOs*);
  
};

struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482 {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE {
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16result6ResultGuRPC15error5ErrorE3Err {
  void* $0;
  
};

struct moonbit_result_3 {
  int tag;
  union { struct _M0TP26biolab8bio__seq9SeqRecord* ok; void* err;  } data;
  
};

struct moonbit_result_1 {
  int tag;
  union { struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* ok; void* err; 
  } data;
  
};

struct moonbit_result_2 {
  int tag;
  union { moonbit_string_t ok; void* err;  } data;
  
};

struct moonbit_result_0 {
  int tag;
  union { int32_t ok; void* err;  } data;
  
};

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__16_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__10_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__7_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__15_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__8_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__3_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__2_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__1_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__13_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__5_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__4_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__9_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__17_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__19_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__12_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__6_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__14_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__11_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__18_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__0_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t
);

moonbit_string_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1489(
  struct _M0TWRPC15error5ErrorEs*,
  void*
);

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1482(
  struct _M0TWssbEu*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN13handle__startS1477(
  struct _M0TWEOc*
);

struct moonbit_result_0 _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__test(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3562l450(
  struct _M0TWEOc*
);

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3558l451(
  struct _M0TWRPC15error5ErrorEu*,
  void*
);

int32_t _M0FP26biolab8bio__seq45moonbit__test__driver__internal__catch__error(
  struct _M0TWEuQRPC15error5Error*,
  struct _M0TWEOc*,
  struct _M0TWRPC15error5ErrorEu*
);

struct _M0TPB5ArrayGUsiEE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__args(
  
);

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1409(
  int32_t,
  moonbit_string_t,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1402(
  int32_t
);

moonbit_string_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395(
  int32_t,
  moonbit_bytes_t
);

int32_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1388(
  int32_t,
  moonbit_string_t
);

#define _M0FP26biolab8bio__seq52moonbit__test__driver__internal__get__cli__args__ffi moonbit_get_cli_args

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP016_24default__implP26biolab8bio__seq28MoonBit__Async__Test__Driver17run__async__testsGRP26biolab8bio__seq34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__19(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__18(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__17(
  
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord9from__seq(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__16(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__15(
  
);

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq15seqio__to__dict(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__14(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__13(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__12(
  
);

struct moonbit_result_3 _M0FP26biolab8bio__seq11seqio__read(
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__11(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__10(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__9(
  
);

struct moonbit_result_1 _M0FP26biolab8bio__seq12seqio__parse(
  moonbit_string_t,
  moonbit_string_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq14parse__genbank(
  moonbit_string_t
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__8(
  
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord*,
  struct _M0TPB5ArrayGiE*
);

struct moonbit_result_2 _M0FP26biolab8bio__seq12write__fastq(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__7(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__6(
  
);

moonbit_string_t _M0FP26biolab8bio__seq12write__fasta(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

moonbit_string_t _M0FP26biolab8bio__seq19build__fasta__title(
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__5(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__4(
  
);

struct _M0TPB5ArrayGiE* _M0MP26biolab8bio__seq9SeqRecord14phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

struct moonbit_result_1 _M0FP26biolab8bio__seq12parse__fastq(
  moonbit_string_t
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__3(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__2(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__1(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__0(
  
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

int32_t _M0IPC15array5ArrayPB2Eq5equalGiE(
  struct _M0TPB5ArrayGiE*,
  struct _M0TPB5ArrayGiE*
);

moonbit_string_t _M0MPC15array5Array2atGsE(struct _M0TPB5ArrayGsE*, int32_t);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

struct _M0TPC16string10StringView _M0MPC15array5Array2atGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  int32_t
);

int32_t _M0MPC15array5Array2atGiE(struct _M0TPB5ArrayGiE*, int32_t);

int32_t _M0IPB9SourceLocPB4Show6output(
  moonbit_string_t,
  struct _M0TPB6Logger
);

int32_t _M0FPB7printlnGsE(moonbit_string_t);

int32_t _M0IPC13int3IntPB4Hash4hash(int32_t);

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

struct _M0TUsRPB5ArrayGiEE* _M0MPB3Map4iterGsRPB5ArrayGiEEC2560l711(
  struct _M0TWEOUsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map6lengthGsRPB5ArrayGiEE(struct _M0TPB3MapGsRPB5ArrayGiEE*);

struct _M0TUWEuQRPC15error5ErrorNsE* _M0MPB3Map3getGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t
);

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3getGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t
);

struct _M0TPB5ArrayGiE* _M0MPB3Map3getGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t
);

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPB3Map3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE,
  int64_t
);

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE,
  int64_t
);

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0MPB3Map3MapGsRPB5ArrayGiEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE,
  int64_t
);

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0MPB3Map3MapGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE,
  int64_t
);

int32_t _M0MPB3Map3setGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  moonbit_string_t,
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0MPB3Map3setGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t,
  struct _M0TUWEuQRPC15error5ErrorNsE*
);

int32_t _M0MPB3Map3setGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGiE*
);

int32_t _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t,
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

int32_t _M0MPB3Map15set__with__hashGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  moonbit_string_t,
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t,
  struct _M0TUWEuQRPC15error5ErrorNsE*,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGiE*,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  moonbit_string_t,
  struct _M0TP26biolab8bio__seq9SeqRecord*,
  int32_t
);

int32_t _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

int32_t _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0MPB3Map4growGsRPB5ArrayGiEE(struct _M0TPB3MapGsRPB5ArrayGiEE*);

int32_t _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

int32_t _M0MPB3Map20rehash__place__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

int32_t _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  int32_t,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

int32_t _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

int32_t _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*,
  int32_t,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPC13int3Int3max(int32_t, int32_t);

int32_t _M0FPB21capacity__for__length(int32_t);

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0FPB8new__mapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  int32_t
);

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0FPB8new__mapGiUWEuQRPC15error5ErrorNsEE(
  int32_t
);

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0FPB8new__mapGsRPB5ArrayGiEE(int32_t);

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FPB8new__mapGsRP26biolab8bio__seq9SeqRecordE(
  int32_t
);

int32_t _M0MPC13int3Int20next__power__of__two(int32_t);

int32_t _M0FPB21calc__grow__threshold(int32_t);

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE(
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0MPC16option6Option6unwrapGRPB5EntryGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*
);

struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0MPC16option6Option6unwrapGRPB5EntryGsRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*
);

struct _M0TPB4IterGsE* _M0MPC15array13ReadOnlyArray4iterGsE(
  moonbit_string_t*
);

struct _M0TPB4IterGsE* _M0MPC15array10FixedArray4iterGsE(moonbit_string_t*);

struct _M0TPB4IterGsE* _M0MPC15array9ArrayView4iterGsE(
  struct _M0TPB9ArrayViewGsE
);

moonbit_string_t _M0MPC15array9ArrayView4iterGsEC2174l729(struct _M0TWEOs*);

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t);

moonbit_string_t _M0IPC14bool4BoolPB4Show10to__string(int32_t);

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

void* _M0MPC16string10StringView5splitC2157l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC2153l1145(
  struct _M0TWcERPC16string10StringView*,
  int32_t
);

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t);

moonbit_string_t _M0FPB16char__to__string(int32_t);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE*,
  struct _M0TWcERPC16string10StringView*
);

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC2146l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

int32_t _M0MPC16string6String8contains(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView8contains(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView20contains__code__unit(
  struct _M0TPC16string10StringView,
  int32_t
);

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

int32_t _M0MPC15array5Array4pushGiE(struct _M0TPB5ArrayGiE*, int32_t);

int32_t _M0MPC15array5Array4pushGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array7reallocGUsiEE(struct _M0TPB5ArrayGUsiEE*);

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPC15array5Array7reallocGiE(struct _M0TPB5ArrayGiE*);

int32_t _M0MPC15array5Array7reallocGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*
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

int32_t _M0MPC15array5Array14resize__bufferGiE(
  struct _M0TPB5ArrayGiE*,
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

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(int32_t);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(
  int32_t
);

struct _M0TPB5ArrayGiE* _M0MPC15array5Array11new_2einnerGiE(int32_t);

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

struct moonbit_result_0 _M0FPB12assert__true(
  int32_t,
  void*,
  moonbit_string_t
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

int32_t _M0MPC16string10StringView4iterC1999l219(struct _M0TWEOc*);

int32_t _M0IPC16string10StringViewPB4Show6output(
  struct _M0TPC16string10StringView,
  struct _M0TPB6Logger
);

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

moonbit_string_t* _M0MPC15array5Array6bufferGsE(struct _M0TPB5ArrayGsE*);

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*
);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

struct _M0TPC16string10StringView* _M0MPC15array5Array6bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*
);

int32_t* _M0MPC15array5Array6bufferGiE(struct _M0TPB5ArrayGiE*);

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView,
  int32_t,
  int64_t
);

struct _M0TPB4IterGsE* _M0MPB4Iter3newGsE(struct _M0TWEOs*, int64_t);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*,
  int64_t
);

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB4Iter3newGUsRPB5ArrayGiEEE(
  struct _M0TWEOUsRPB5ArrayGiEE*,
  int64_t
);

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(struct _M0TWEOc*, int64_t);

struct moonbit_result_0 _M0FPB10assert__eqGiE(
  int32_t,
  int32_t,
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_0 _M0FPB10assert__eqGsE(
  moonbit_string_t,
  moonbit_string_t,
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_0 _M0FPB4failGuE(
  struct _M0TPC16string10StringView,
  moonbit_string_t
);

moonbit_string_t _M0FPB13debug__stringGiE(int32_t);

moonbit_string_t _M0FPB13debug__stringGsE(moonbit_string_t);

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

moonbit_string_t _M0MPB4Iter4nextGsE(struct _M0TPB4IterGsE*);

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE*
);

struct _M0TUsRPB5ArrayGiEE* _M0MPB4Iter4nextGUsRPB5ArrayGiEEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE*
);

int32_t _M0MPB4Iter4nextGcE(struct _M0TPB4IterGcE*);

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

int32_t _M0IP016_24default__implPB4Show6outputGbE(
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

int32_t _M0IP016_24default__implPB2Eq10not__equalGsE(
  moonbit_string_t,
  moonbit_string_t
);

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

int32_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGiE(
  int32_t*,
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

int32_t _M0MPB13StringBuilder13write__objectGbE(
  struct _M0TPB13StringBuilder*,
  int32_t
);

int32_t _M0MPB13StringBuilder13write__objectGRPB9SourceLocE(
  struct _M0TPB13StringBuilder*,
  moonbit_string_t
);

int32_t _M0MPB13StringBuilder13write__objectGRPC16string10StringViewE(
  struct _M0TPB13StringBuilder*,
  struct _M0TPC16string10StringView
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

int32_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGiE(
  int32_t*,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGiE(
  int32_t*,
  int32_t,
  int32_t*,
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGiEE(
  int32_t*,
  int32_t,
  int32_t*,
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

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t*);

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(struct _M0TUsiE**);

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord**
);

int32_t _M0MPB18UninitializedArray6lengthGiE(int32_t*);

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

int32_t _M0FPC15abort5abortGuE(moonbit_string_t);

uint16_t* _M0FPC15abort5abortGAkE(moonbit_string_t);

moonbit_string_t* _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(
  moonbit_string_t
);

struct _M0TUsiE** _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(
  moonbit_string_t
);

moonbit_string_t _M0FPC15abort5abortGsE(moonbit_string_t);

int32_t _M0FPC15abort5abortGiE(moonbit_string_t);

struct _M0TPC16string10StringView _M0FPC15abort5abortGRPC16string10StringViewE(
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(
  moonbit_string_t
);

int32_t* _M0FPC15abort5abortGRPB18UninitializedArrayGiEE(moonbit_string_t);

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

struct { int32_t rc; uint32_t meta; uint16_t const data[1]; 
} const moonbit_string_literal_0 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 0, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_71 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 52, 48, 58, 51, 45, 49, 52, 48, 58, 50, 50, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_51 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 55, 57, 58, 55, 45, 49, 55, 57, 58, 50, 55, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_48 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 97, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_213 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 115, 101, 
    113, 114, 101, 99, 111, 114, 100, 95, 110, 101, 119, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_154 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 50, 58, 51, 45, 50, 50, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_42 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 57, 56, 58, 51, 45, 49, 57, 56, 58, 50, 56, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_212 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 115, 101, 
    113, 105, 111, 95, 116, 111, 95, 100, 105, 99, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_100 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 43, 10, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_83 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 49, 54, 58, 51, 45, 49, 49, 54, 58, 52, 56, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_61 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 62, 111, 
    110, 108, 121, 32, 79, 110, 101, 10, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_37 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 60, 117, 
    110, 107, 110, 111, 119, 110, 32, 110, 97, 109, 101, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[53]; 
} const moonbit_string_literal_196 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 52, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_177 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 110, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_204 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 119, 114, 
    105, 116, 101, 95, 102, 97, 115, 116, 97, 95, 119, 114, 97, 112, 
    112, 105, 110, 103, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_117 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 48, 58, 51, 45, 54, 48, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_68 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 77, 111, 
    114, 101, 32, 116, 104, 97, 110, 32, 111, 110, 101, 32, 114, 101, 
    99, 111, 114, 100, 32, 102, 111, 117, 110, 100, 32, 105, 110, 32, 
    104, 97, 110, 100, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_202 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 113, 95, 109, 117, 108, 116, 
    105, 112, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_153 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 109, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[44]; 
} const moonbit_string_literal_121 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 43, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 53, 58, 49, 51, 45, 54, 53, 58, 53, 48, 64, 98, 105, 111, 
    108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_175 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_95 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 112, 104, 
    114, 101, 100, 95, 113, 117, 97, 108, 105, 116, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_9 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 64, 115, 
    49, 32, 68, 101, 115, 99, 10, 65, 67, 71, 84, 10, 43, 10, 73, 73, 
    73, 73, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_21 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 50, 55, 58, 51, 45, 50, 50, 55, 58, 51, 52, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_57 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    115, 116, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[29]; 
} const moonbit_string_literal_138 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 28, 32, 100, 
    111, 110, 39, 116, 32, 109, 97, 116, 99, 104, 32, 105, 110, 32, 70, 
    65, 83, 84, 81, 32, 114, 101, 99, 111, 114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_75 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 50, 53, 58, 51, 45, 49, 50, 53, 58, 52, 56, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_105 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 57, 55, 58, 51, 45, 57, 55, 58, 52, 50, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_85 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 85, 110, 
    107, 110, 111, 119, 110, 32, 102, 111, 114, 109, 97, 116, 58, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_139 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 62, 108, 
    49, 32, 76, 111, 119, 101, 114, 10, 97, 99, 103, 116, 97, 99, 103, 
    116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_65 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 52, 55, 58, 51, 45, 49, 52, 55, 58, 52, 49, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_191 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 41, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_164 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 53, 58, 51, 45, 49, 53, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_30 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 51, 50, 58, 51, 45, 50, 51, 50, 58, 52, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_115 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 32, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_189 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 46, 108, 101, 110, 103, 116, 104, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_101 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 65, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_10 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 49, 58, 51, 45, 50, 52, 49, 58, 51, 52, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_186 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_147 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 115, 101, 
    113, 51, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[30]; 
} const moonbit_string_literal_126 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 29, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 112, 104, 114, 101, 100, 95, 113, 
    117, 97, 108, 105, 116, 121, 32, 102, 111, 114, 32, 115, 50, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_18 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 55, 58, 49, 51, 45, 50, 52, 55, 58, 53, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_142 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 52, 48, 58, 51, 45, 52, 48, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_25 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 50, 57, 58, 51, 45, 50, 50, 57, 58, 53, 49, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_72 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    115, 116, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_36 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 60, 117, 
    110, 107, 110, 111, 119, 110, 32, 105, 100, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[36]; 
} const moonbit_string_literal_31 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 35, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 110, 111, 32, 112, 104, 114, 101, 
    100, 95, 113, 117, 97, 108, 105, 116, 121, 32, 105, 110, 105, 116, 
    105, 97, 108, 108, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_56 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 56, 53, 58, 49, 54, 45, 49, 56, 53, 58, 53, 50, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[30]; 
} const moonbit_string_literal_120 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 29, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 112, 104, 114, 101, 100, 95, 113, 
    117, 97, 108, 105, 116, 121, 32, 102, 111, 114, 32, 115, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_180 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 116, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_178 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 114, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_170 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 96, 32, 
    105, 115, 32, 110, 111, 116, 32, 116, 114, 117, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_124 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 55, 58, 51, 45, 54, 55, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_66 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 52, 56, 58, 51, 45, 49, 52, 56, 58, 52, 49, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[44]; 
} const moonbit_string_literal_135 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 43, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 53, 50, 58, 49, 54, 45, 53, 50, 58, 53, 48, 64, 98, 105, 111, 
    108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_15 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 65, 67, 
    71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_114 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 56, 53, 58, 51, 45, 56, 53, 58, 53, 54, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_103 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 62, 108, 
    111, 110, 103, 49, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_96 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 77, 105, 
    115, 115, 105, 110, 103, 32, 112, 104, 114, 101, 100, 95, 113, 117, 
    97, 108, 105, 116, 121, 32, 97, 110, 110, 111, 116, 97, 116, 105, 
    111, 110, 32, 102, 111, 114, 32, 114, 101, 99, 111, 114, 100, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_155 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 109, 49, 
    32, 77, 117, 108, 116, 105, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_86 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 103, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_78 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 50, 55, 58, 49, 54, 45, 49, 50, 55, 58, 53, 48, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_53 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_92 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 102, 113, 
    49, 32, 114, 101, 99, 111, 114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_33 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 112, 104, 114, 101, 100, 95, 113, 
    117, 97, 108, 105, 116, 121, 32, 97, 102, 116, 101, 114, 32, 115, 
    101, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_22 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 115, 101, 
    113, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_161 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 50, 58, 51, 45, 49, 50, 58, 53, 48, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_52 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 56, 48, 58, 55, 45, 49, 56, 48, 58, 52, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_98 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 32, 97, 
    110, 100, 32, 113, 117, 97, 108, 105, 116, 121, 32, 108, 101, 110, 
    103, 116, 104, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_205 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 119, 114, 
    105, 116, 101, 95, 102, 97, 115, 116, 113, 95, 98, 97, 115, 105, 
    99, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_160 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 49, 58, 51, 45, 49, 49, 58, 51, 53, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_102 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 108, 111, 
    110, 103, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_59 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 62, 115, 
    49, 10, 65, 67, 71, 84, 10, 62, 115, 50, 10, 84, 84, 84, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_91 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 102, 113, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[26]; 
} const moonbit_string_literal_210 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 25, 115, 101, 
    113, 105, 111, 95, 114, 101, 97, 100, 95, 109, 117, 108, 116, 105, 
    112, 108, 101, 95, 101, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_44 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 48, 48, 58, 51, 45, 50, 48, 48, 58, 51, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_39 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 109, 121, 
    105, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[48]; 
} const moonbit_string_literal_137 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 47, 69, 120, 
    112, 101, 99, 116, 101, 100, 32, 39, 43, 39, 32, 97, 116, 32, 115, 
    116, 97, 114, 116, 32, 111, 102, 32, 70, 65, 83, 84, 81, 32, 115, 
    101, 112, 97, 114, 97, 116, 111, 114, 44, 32, 103, 111, 116, 58, 
    32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_79 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 62, 115, 
    49, 32, 72, 101, 108, 108, 111, 10, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_58 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 55, 48, 58, 51, 45, 49, 55, 48, 58, 50, 50, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_162 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 51, 58, 51, 45, 49, 51, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_90 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 60, 117, 
    110, 107, 110, 111, 119, 110, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_69 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 62, 115, 
    49, 10, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_158 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 58, 51, 45, 50, 52, 58, 53, 54, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_35 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 49, 54, 58, 49, 54, 45, 50, 49, 54, 58, 53, 48, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_28 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 51, 49, 58, 51, 45, 50, 51, 49, 58, 51, 54, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_14 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 51, 58, 51, 45, 50, 52, 51, 58, 52, 56, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_133 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 112, 104, 114, 101, 100, 95, 113, 
    117, 97, 108, 105, 116, 121, 32, 116, 111, 32, 98, 101, 32, 112, 
    114, 101, 115, 101, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_84 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 49, 55, 58, 51, 45, 49, 49, 55, 58, 52, 56, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_12 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 50, 58, 51, 45, 50, 52, 50, 58, 51, 52, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_150 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 51, 51, 58, 51, 45, 51, 51, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_27 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 115, 101, 
    113, 50, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[83]; 
} const moonbit_string_literal_194 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 82, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 77, 111, 
    111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 
    114, 73, 110, 116, 101, 114, 110, 97, 108, 83, 107, 105, 112, 84, 
    101, 115, 116, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 
    116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 
    108, 83, 107, 105, 112, 84, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_172 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 32, 33, 
    61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_163 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 52, 58, 51, 45, 49, 52, 58, 51, 53, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_64 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 111, 110, 
    108, 121, 32, 79, 110, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_216 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 114, 111, 
    117, 110, 100, 116, 114, 105, 112, 95, 102, 97, 115, 116, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_201 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 113, 95, 98, 97, 115, 105, 
    99, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_107 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 114, 49, 
    32, 70, 105, 114, 115, 116, 32, 114, 101, 99, 111, 114, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_217 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_125 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 56, 58, 51, 45, 54, 56, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_87 =
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

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_132 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 53, 48, 58, 51, 45, 53, 48, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_88 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 102, 97, 
    115, 116, 113, 45, 115, 97, 110, 103, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[29]; 
} const moonbit_string_literal_47 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 28, 62, 97, 
    32, 70, 105, 114, 115, 116, 10, 65, 67, 71, 84, 10, 62, 98, 32, 83, 
    101, 99, 111, 110, 100, 10, 84, 84, 84, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_206 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 115, 101, 
    113, 105, 111, 95, 112, 97, 114, 115, 101, 95, 102, 97, 115, 116, 
    97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_182 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 117, 114, 114, 111, 103, 97, 116, 
    101, 32, 112, 97, 105, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_106 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 114, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[27]; 
} const moonbit_string_literal_67 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 78, 111, 
    32, 114, 101, 99, 111, 114, 100, 115, 32, 102, 111, 117, 110, 100, 
    32, 105, 110, 32, 104, 97, 110, 100, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_62 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 111, 110, 
    108, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_38 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 60, 117, 
    110, 107, 110, 111, 119, 110, 32, 100, 101, 115, 99, 114, 105, 112, 
    116, 105, 111, 110, 62, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_76 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 112, 104, 114, 101, 100, 95, 113, 
    117, 97, 108, 105, 116, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_23 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 50, 56, 58, 51, 45, 50, 50, 56, 58, 51, 54, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_32 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 49, 48, 58, 49, 54, 45, 50, 49, 48, 58, 53, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_218 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 98, 105, 
    111, 95, 115, 101, 113, 95, 119, 98, 116, 101, 115, 116, 46, 109, 
    98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[31]; 
} const moonbit_string_literal_174 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 30, 114, 97, 
    100, 105, 120, 32, 109, 117, 115, 116, 32, 98, 101, 32, 98, 101, 
    116, 119, 101, 101, 110, 32, 50, 32, 97, 110, 100, 32, 51, 54, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_3 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 114, 101, 115, 117, 108, 116, 34, 
    44, 34, 102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_173 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 32, 70, 
    65, 73, 76, 69, 68, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_207 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 115, 101, 
    113, 105, 111, 95, 112, 97, 114, 115, 101, 95, 102, 97, 115, 116, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_81 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 49, 53, 58, 51, 45, 49, 49, 53, 58, 51, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[28]; 
} const moonbit_string_literal_20 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 27, 62, 115, 
    101, 113, 49, 32, 72, 101, 108, 108, 111, 10, 65, 67, 71, 84, 10, 
    62, 115, 101, 113, 50, 10, 84, 84, 84, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_200 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 108, 111, 119, 101, 
    114, 99, 97, 115, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_41 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 109, 121, 
    100, 101, 115, 99, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_26 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 51, 48, 58, 51, 45, 50, 51, 48, 58, 52, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_151 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 62, 109, 
    49, 32, 77, 117, 108, 116, 105, 10, 65, 67, 71, 84, 10, 84, 84, 84, 
    84, 10, 71, 71, 71, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_93 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 64, 102, 
    113, 49, 32, 114, 101, 99, 111, 114, 100, 10, 65, 67, 71, 84, 10, 
    43, 10, 73, 63, 53, 43, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_215 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 114, 111, 
    117, 110, 100, 116, 114, 105, 112, 95, 102, 97, 115, 116, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_188 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 44, 32, 
    108, 101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_166 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    108, 115, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[32]; 
} const moonbit_string_literal_113 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 31, 62, 114, 
    49, 32, 70, 105, 114, 115, 116, 32, 114, 101, 99, 111, 114, 100, 
    10, 65, 67, 71, 84, 10, 62, 114, 50, 10, 84, 84, 84, 84, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_185 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 111, 
    117, 110, 100, 115, 32, 99, 104, 101, 99, 107, 32, 102, 97, 105, 
    108, 101, 100, 58, 32, 97, 108, 108, 111, 99, 97, 116, 101, 95, 108, 
    101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[56]; 
} const moonbit_string_literal_183 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 55, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 116, 97, 114, 116, 32, 111, 114, 
    32, 101, 110, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 114, 
    32, 83, 116, 114, 105, 110, 103, 58, 58, 99, 111, 100, 101, 112, 
    111, 105, 110, 116, 95, 108, 101, 110, 103, 116, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_140 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 51, 57, 58, 51, 45, 51, 57, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_45 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 48, 49, 58, 51, 45, 50, 48, 49, 58, 52, 49, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_104 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 57, 54, 58, 51, 45, 57, 54, 58, 50, 55, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_16 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 52, 58, 51, 45, 50, 52, 52, 58, 52, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_144 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 52, 49, 58, 51, 45, 52, 49, 58, 53, 50, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_50 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 56, 50, 58, 49, 51, 45, 49, 56, 50, 58, 52, 53, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_19 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 52, 54, 58, 49, 54, 45, 50, 52, 54, 58, 53, 48, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[27]; 
} const moonbit_string_literal_208 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 115, 101, 
    113, 105, 111, 95, 112, 97, 114, 115, 101, 95, 117, 110, 107, 110, 
    111, 119, 110, 95, 102, 111, 114, 109, 97, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_176 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 48, 49, 
    50, 51, 52, 53, 54, 55, 56, 57, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 119, 120, 121, 122, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_148 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 51, 49, 58, 51, 45, 51, 49, 58, 51, 53, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[51]; 
} const moonbit_string_literal_192 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 50, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 73, 110, 115, 112, 101, 
    99, 116, 69, 114, 114, 111, 114, 46, 73, 110, 115, 112, 101, 99, 
    116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_89 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 10, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[44]; 
} const moonbit_string_literal_134 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 43, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 53, 51, 58, 49, 51, 45, 53, 51, 58, 53, 55, 64, 98, 105, 111, 
    108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_118 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 49, 58, 51, 45, 54, 49, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_108 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 114, 50, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_7 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 45, 45, 
    45, 45, 45, 32, 69, 78, 68, 32, 77, 79, 79, 78, 32, 84, 69, 83, 84, 
    32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_159 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 48, 58, 51, 45, 49, 48, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_129 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 52, 55, 58, 51, 45, 52, 55, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_8 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 115, 116, 97, 114, 116, 34, 44, 34, 
    102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_4 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 34, 
    105, 110, 100, 101, 120, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[81]; 
} const moonbit_string_literal_193 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 80, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 77, 111, 
    111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 
    114, 73, 110, 116, 101, 114, 110, 97, 108, 74, 115, 69, 114, 114, 
    111, 114, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 
    68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 
    74, 115, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_123 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 115, 50, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_55 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 56, 54, 58, 49, 51, 45, 49, 56, 54, 58, 52, 53, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_209 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 115, 101, 
    113, 105, 111, 95, 114, 101, 97, 100, 95, 115, 105, 110, 103, 108, 
    101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_190 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 70, 97, 
    105, 108, 117, 114, 101, 40, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[44]; 
} const moonbit_string_literal_127 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 43, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 55, 49, 58, 49, 51, 45, 55, 49, 58, 53, 48, 64, 98, 105, 111, 
    108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_109 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 62, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_80 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 49, 52, 58, 51, 45, 49, 49, 52, 58, 51, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[27]; 
} const moonbit_string_literal_199 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 110, 111, 95, 100, 
    101, 115, 99, 114, 105, 112, 116, 105, 111, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_168 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 114, 101, 
    112, 101, 97, 116, 32, 114, 101, 115, 117, 108, 116, 32, 116, 111, 
    111, 32, 108, 97, 114, 103, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_99 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 32, 100, 
    111, 110, 39, 116, 32, 109, 97, 116, 99, 104, 32, 105, 110, 32, 114, 
    101, 99, 111, 114, 100, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_49 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 107, 101, 121, 32, 39, 97, 39, 32, 
    105, 110, 32, 100, 105, 99, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_6 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 125, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_97 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 83, 101, 
    113, 117, 101, 110, 99, 101, 32, 108, 101, 110, 103, 116, 104, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_29 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 84, 84, 
    84, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_24 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 115, 101, 
    113, 49, 32, 72, 101, 108, 108, 111, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[35]; 
} const moonbit_string_literal_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 34, 45, 45, 
    45, 45, 45, 32, 66, 69, 71, 73, 78, 32, 77, 79, 79, 78, 32, 84, 69, 
    83, 84, 32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_203 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 119, 114, 
    105, 116, 101, 95, 102, 97, 115, 116, 97, 95, 98, 97, 115, 105, 99, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_82 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 115, 49, 
    32, 72, 101, 108, 108, 111, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_187 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    100, 115, 116, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_184 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 105, 110, 
    118, 97, 108, 105, 100, 32, 99, 111, 100, 101, 32, 112, 111, 105, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_63 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 52, 54, 58, 51, 45, 49, 52, 54, 58, 50, 56, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_46 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 48, 50, 58, 51, 45, 50, 48, 50, 58, 50, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[28]; 
} const moonbit_string_literal_214 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 27, 115, 101, 
    113, 114, 101, 99, 111, 114, 100, 95, 115, 101, 116, 95, 112, 104, 
    114, 101, 100, 95, 113, 117, 97, 108, 105, 116, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_112 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 56, 52, 58, 51, 45, 56, 52, 58, 51, 54, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_5 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 44, 34, 
    109, 101, 115, 115, 97, 103, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_143 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 97, 99, 
    103, 116, 97, 99, 103, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_149 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 51, 50, 58, 51, 45, 51, 50, 58, 52, 52, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_111 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 56, 51, 58, 51, 45, 56, 51, 58, 51, 54, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_171 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 73, 110, 
    118, 97, 108, 105, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 
    114, 32, 86, 105, 101, 119, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[39]; 
} const moonbit_string_literal_17 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 38, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 112, 104, 114, 101, 100, 95, 113, 
    117, 97, 108, 105, 116, 121, 32, 97, 102, 116, 101, 114, 32, 114, 
    111, 117, 110, 100, 116, 114, 105, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_197 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 98, 97, 115, 105, 99, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_169 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 96, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_136 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 69, 120, 
    112, 101, 99, 116, 101, 100, 32, 39, 64, 39, 32, 97, 116, 32, 115, 
    116, 97, 114, 116, 32, 111, 102, 32, 70, 65, 83, 84, 81, 32, 104, 
    101, 97, 100, 101, 114, 44, 32, 103, 111, 116, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_94 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 48, 56, 58, 51, 45, 49, 48, 56, 58, 52, 57, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_13 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 115, 49, 
    32, 68, 101, 115, 99, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_141 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 108, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_130 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 52, 56, 58, 51, 45, 52, 56, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_181 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 92, 117, 123, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_157 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 65, 67, 
    71, 84, 84, 84, 84, 84, 71, 71, 71, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_74 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 50, 52, 58, 51, 45, 49, 50, 52, 58, 51, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_198 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 109, 117, 108, 116, 
    105, 108, 105, 110, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_11 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 115, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_165 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 116, 114, 
    117, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_131 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 52, 57, 58, 51, 45, 52, 57, 58, 52, 55, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[44]; 
} const moonbit_string_literal_128 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 43, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 55, 48, 58, 49, 54, 45, 55, 48, 58, 53, 48, 64, 98, 105, 111, 
    108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_179 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_211 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 115, 101, 
    113, 105, 111, 95, 114, 101, 97, 100, 95, 101, 109, 112, 116, 121, 
    95, 101, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_43 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 57, 57, 58, 51, 45, 49, 57, 57, 58, 51, 50, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_73 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 50, 51, 58, 51, 45, 49, 50, 51, 58, 51, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_167 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 110, 101, 
    103, 97, 116, 105, 118, 101, 32, 114, 101, 112, 101, 97, 116, 32, 
    99, 111, 117, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_116 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 64, 115, 
    49, 32, 68, 101, 115, 99, 10, 65, 67, 71, 84, 10, 43, 10, 73, 73, 
    73, 73, 10, 64, 115, 50, 32, 79, 116, 104, 101, 114, 10, 84, 84, 
    84, 84, 10, 43, 10, 74, 74, 74, 74, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_145 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 62, 115, 
    101, 113, 51, 10, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_119 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 50, 58, 51, 45, 54, 50, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_77 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 50, 56, 58, 49, 51, 45, 49, 50, 56, 58, 52, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_146 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 51, 48, 58, 51, 45, 51, 48, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_54 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 101, 120, 
    112, 101, 99, 116, 101, 100, 32, 107, 101, 121, 32, 39, 98, 39, 32, 
    105, 110, 32, 100, 105, 99, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_156 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 51, 58, 51, 45, 50, 51, 58, 52, 56, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_110 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 56, 50, 58, 51, 45, 56, 50, 58, 51, 53, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[45]; 
} const moonbit_string_literal_60 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 44, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 49, 53, 57, 58, 51, 45, 49, 53, 57, 58, 50, 50, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_152 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 49, 58, 51, 45, 50, 49, 58, 51, 51, 64, 98, 105, 111, 108, 
    97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[44]; 
} const moonbit_string_literal_122 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 43, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 54, 52, 58, 49, 54, 45, 54, 52, 58, 53, 48, 64, 98, 105, 111, 
    108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[46]; 
} const moonbit_string_literal_34 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 45, 115, 101, 
    113, 105, 111, 95, 119, 98, 116, 101, 115, 116, 46, 109, 98, 116, 
    58, 50, 49, 55, 58, 49, 51, 45, 50, 49, 55, 58, 53, 51, 64, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_70 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 120, 121, 
    122, 95, 117, 110, 107, 110, 111, 119, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_195 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 83, 101, 
    113, 73, 79, 69, 114, 114, 111, 114, 46, 83, 101, 113, 73, 79, 69, 
    114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_40 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 109, 121, 
    110, 97, 109, 101, 0
  };

struct moonbit_object const moonbit_constant_constructor_0 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0)
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__0_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__0_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__3_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__3_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__17_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__17_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__19_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__19_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__16_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__16_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__1_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__1_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__10_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__10_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__11_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__11_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__18_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__18_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__2_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__2_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__13_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__13_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__5_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__5_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__7_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__7_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__4_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__4_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__15_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__15_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__14_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__14_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__6_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__6_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__8_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__8_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__12_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__12_2edyncall
  };

struct {
  int32_t rc;
  uint32_t meta;
  struct _M0TWcERPC16string10StringView data;
  
} const _M0MPC16string10StringView5splitC2153l1145$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0MPC16string10StringView5splitC2153l1145
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__9_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__9_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWRPC15error5ErrorEs data; 
} const _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1489$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1489
  };

uint32_t const moonbit_layout_table_data[176] =
  {
    sizeof(struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477)
    / 4, 1,
    offsetof(struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477, $1)
    / 4,
    sizeof(struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482)
    / 4, 1,
    offsetof(struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482, $1)
    / 4,
    sizeof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest, $0)
    / 4,
    sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__)
    / 4, 2,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__, $0)
    / 4,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__, $1)
    / 4,
    sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__)
    / 4, 3,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__, $0)
    / 4,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__, $1)
    / 4,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__, $2)
    / 4, sizeof(struct _M0TPB5ArrayGUsiEE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGUsiEE, $0) / 4, sizeof(struct _M0TUsiE) / 4,
    1, offsetof(struct _M0TUsiE, $0) / 4, sizeof(struct _M0TPB5ArrayGsE) / 4,
    1, offsetof(struct _M0TPB5ArrayGsE, $0) / 4,
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
    sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGOsE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGOsE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGRPB13StringBuilderE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGRPB13StringBuilderE, $0) / 4,
    sizeof(struct _M0TP26biolab8bio__seq3Seq) / 4, 1,
    offsetof(struct _M0TP26biolab8bio__seq3Seq, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE, $0) / 4,
    sizeof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__)
    / 4, 2,
    offsetof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__, $0)
    / 4,
    offsetof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__, $1)
    / 4, sizeof(struct _M0TUsRPB5ArrayGiEE) / 4, 2,
    offsetof(struct _M0TUsRPB5ArrayGiEE, $0) / 4,
    offsetof(struct _M0TUsRPB5ArrayGiEE, $1) / 4,
    sizeof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE) / 4, 
    3,
    offsetof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $1)
    / 4,
    offsetof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $4)
    / 4,
    offsetof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $5)
    / 4, sizeof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE) / 4, 
    2, offsetof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE, $1) / 4,
    offsetof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE, $5) / 4,
    sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE) / 4, 3,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $5) / 4,
    sizeof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE) / 4, 
    3, offsetof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE) / 4, 
    2,
    offsetof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $5) / 4,
    sizeof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE) / 4, 2,
    offsetof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE, $0) / 4,
    offsetof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRPB5ArrayGiEE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGiEE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGiEE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE, $5) / 4,
    sizeof(struct _M0TPB9ArrayViewGsE) / 4, 1,
    offsetof(struct _M0TPB9ArrayViewGsE, $0) / 4,
    sizeof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__)
    / 4, 2,
    (offsetof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__, $0)
     + offsetof(struct _M0TPB9ArrayViewGsE, $0))
    / 4,
    offsetof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__, $2)
    / 4, sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0TPC16string10StringView) / 4, 1,
    offsetof(struct _M0TPC16string10StringView, $0) / 4,
    sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some) / 4,
    1,
    (offsetof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some, $0)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE) / 4, 
    1, offsetof(struct _M0TPB8MutLocalGORPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__) / 4, 
    2,
    offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__, $0)
    / 4,
    (offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__, $1)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4,
    sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__)
    / 4, 2,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__, $0)
    / 4,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__, $1)
    / 4, sizeof(struct _M0TPB4IterGRPC16string10StringViewE) / 4, 1,
    offsetof(struct _M0TPB4IterGRPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__) / 4, 
    2,
    offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__, $0) / 4,
    (offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__, $2)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB4IterGsE) / 4, 1,
    offsetof(struct _M0TPB4IterGsE, $0) / 4,
    sizeof(struct _M0TPB4IterGUsRPB5ArrayGiEEE) / 4, 1,
    offsetof(struct _M0TPB4IterGUsRPB5ArrayGiEEE, $0) / 4,
    sizeof(struct _M0TPB4IterGcE) / 4, 1,
    offsetof(struct _M0TPB4IterGcE, $0) / 4,
    sizeof(struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure, $0)
    / 4, sizeof(struct _M0TPB6Logger) / 4, 2,
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
    sizeof(struct _M0TUWEuQRPC15error5ErrorNsE) / 4, 2,
    offsetof(struct _M0TUWEuQRPC15error5ErrorNsE, $0) / 4,
    offsetof(struct _M0TUWEuQRPC15error5ErrorNsE, $1) / 4,
    sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE) / 4, 1,
    offsetof(struct _M0TUiUWEuQRPC15error5ErrorNsEE, $1) / 4,
    sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE) / 4, 2,
    offsetof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $0) / 4,
    offsetof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $1) / 4,
    sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE, $0) / 4
  };

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__0_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__0_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__18_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__18_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__11_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__11_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__14_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__14_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__6_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__6_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__12_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__12_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__19_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__19_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__17_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__17_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__9_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__9_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__4_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__4_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__5_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__5_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__13_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__13_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__1_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__1_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__2_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__2_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__3_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__3_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__8_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__8_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__15_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__15_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__7_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__7_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__10_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__10_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__16_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__16_2edyncall$closure.data;

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

int64_t _M0MPB4Iter4nextN6constrS9980GsE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GsE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GUsRPB5ArrayGiEEE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GUsRPB5ArrayGiEEE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GcE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GcE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GsE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GUsRPB5ArrayGiEEE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GcE = 0ll;

int64_t _M0FPB28boyer__moore__horspool__findN6constrS9990 = 0ll;

int64_t _M0FPB18brute__force__findN6constrS9991 = 0ll;

int32_t _M0FP26biolab8bio__seq18fasta__line__width = 60;

uint16_t* _M0FP26biolab8bio__seq16uppercase__table;

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0FP26biolab8bio__seq48moonbit__test__driver__internal__no__args__tests;

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__16_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3604
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__16();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__10_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3603
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__10();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__7_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3602
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__7();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__15_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3601
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__15();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__8_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3600
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__8();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__3_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3599
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__3();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__2_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3598
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__2();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__1_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3597
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__1();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__13_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3596
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__13();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__5_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3595
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__5();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__4_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3594
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__4();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__9_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3593
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__9();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__17_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3592
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__17();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__19_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3591
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__19();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__12_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3590
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__12();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__6_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3589
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__6();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__14_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3588
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__14();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__11_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3587
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__11();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq56____test__736571696f5f7762746573742e6d6274__18_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3586
) {
  return _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__18();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq55____test__736571696f5f7762746573742e6d6274__0_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3585
) {
  return _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__0();
}

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1510,
  moonbit_string_t _M0L8filenameS1479,
  int32_t _M0L5indexS1481
) {
  struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477* _closure_4178;
  struct _M0TWEOc* _M0L13handle__startS1477;
  struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482* _closure_4179;
  struct _M0TWssbEu* _M0L14handle__resultS1482;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1489;
  void* _M0L11_2atry__errS1504;
  struct moonbit_result_0 _tmp_4181;
  int32_t _handle__error__result_4182;
  int32_t _M0L6_2atmpS3573;
  void* _M0L3errS1505;
  moonbit_string_t _M0L4nameS1507;
  struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest* _M0L36_2aMoonBitTestDriverInternalSkipTestS1508;
  moonbit_string_t _M0L8_2afieldS3605;
  int32_t _M0L6_2acntS3921;
  moonbit_string_t _M0L7_2anameS1509;
  #line 579 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  moonbit_incref(_M0L8filenameS1479);
  _closure_4178
  = (struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477*)moonbit_malloc(sizeof(struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477));
  Moonbit_object_header(_closure_4178)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 0, 0);
  _closure_4178->code
  = &_M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN13handle__startS1477;
  _closure_4178->$0 = _M0L5indexS1481;
  _closure_4178->$1 = _M0L8filenameS1479;
  _M0L13handle__startS1477 = (struct _M0TWEOc*)_closure_4178;
  moonbit_incref(_M0L8filenameS1479);
  _closure_4179
  = (struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482*)moonbit_malloc(sizeof(struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482));
  Moonbit_object_header(_closure_4179)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 3, 0);
  _closure_4179->code
  = &_M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1482;
  _closure_4179->$0 = _M0L5indexS1481;
  _closure_4179->$1 = _M0L8filenameS1479;
  _M0L14handle__resultS1482 = (struct _M0TWssbEu*)_closure_4179;
  _M0L17error__to__stringS1489
  = (struct _M0TWRPC15error5ErrorEs*)&_M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1489$closure.data;
  #line 621 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _tmp_4181
  = _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__test(_M0L12async__testsS1510, _M0L8filenameS1479, _M0L5indexS1481, _M0L13handle__startS1477, _M0L14handle__resultS1482, _M0L17error__to__stringS1489);
  if (_tmp_4181.tag) {
    int32_t const _M0L5_2aokS3582 = _tmp_4181.data.ok;
    _handle__error__result_4182 = _M0L5_2aokS3582;
  } else {
    void* const _M0L6_2aerrS3583 = _tmp_4181.data.err;
    moonbit_decref(_M0L17error__to__stringS1489);
    moonbit_decref(_M0L13handle__startS1477);
    _M0L11_2atry__errS1504 = _M0L6_2aerrS3583;
    goto join_1503;
  }
  if (_handle__error__result_4182) {
    moonbit_decref(_M0L17error__to__stringS1489);
    moonbit_decref(_M0L13handle__startS1477);
    _M0L6_2atmpS3573 = 1;
  } else {
    struct moonbit_result_0 _tmp_4183;
    int32_t _handle__error__result_4184;
    #line 624 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
    _tmp_4183
    = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq43MoonBit__Test__Driver__Internal__With__ArgsE(_M0L12async__testsS1510, _M0L8filenameS1479, _M0L5indexS1481, _M0L13handle__startS1477, _M0L14handle__resultS1482, _M0L17error__to__stringS1489);
    if (_tmp_4183.tag) {
      int32_t const _M0L5_2aokS3580 = _tmp_4183.data.ok;
      _handle__error__result_4184 = _M0L5_2aokS3580;
    } else {
      void* const _M0L6_2aerrS3581 = _tmp_4183.data.err;
      moonbit_decref(_M0L17error__to__stringS1489);
      moonbit_decref(_M0L13handle__startS1477);
      _M0L11_2atry__errS1504 = _M0L6_2aerrS3581;
      goto join_1503;
    }
    if (_handle__error__result_4184) {
      moonbit_decref(_M0L17error__to__stringS1489);
      moonbit_decref(_M0L13handle__startS1477);
      _M0L6_2atmpS3573 = 1;
    } else {
      struct moonbit_result_0 _tmp_4185;
      int32_t _handle__error__result_4186;
      #line 627 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _tmp_4185
      = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq48MoonBit__Test__Driver__Internal__Async__No__ArgsE(_M0L12async__testsS1510, _M0L8filenameS1479, _M0L5indexS1481, _M0L13handle__startS1477, _M0L14handle__resultS1482, _M0L17error__to__stringS1489);
      if (_tmp_4185.tag) {
        int32_t const _M0L5_2aokS3578 = _tmp_4185.data.ok;
        _handle__error__result_4186 = _M0L5_2aokS3578;
      } else {
        void* const _M0L6_2aerrS3579 = _tmp_4185.data.err;
        moonbit_decref(_M0L17error__to__stringS1489);
        moonbit_decref(_M0L13handle__startS1477);
        _M0L11_2atry__errS1504 = _M0L6_2aerrS3579;
        goto join_1503;
      }
      if (_handle__error__result_4186) {
        moonbit_decref(_M0L17error__to__stringS1489);
        moonbit_decref(_M0L13handle__startS1477);
        _M0L6_2atmpS3573 = 1;
      } else {
        struct moonbit_result_0 _tmp_4187;
        int32_t _handle__error__result_4188;
        #line 630 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        _tmp_4187
        = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__Async__With__ArgsE(_M0L12async__testsS1510, _M0L8filenameS1479, _M0L5indexS1481, _M0L13handle__startS1477, _M0L14handle__resultS1482, _M0L17error__to__stringS1489);
        if (_tmp_4187.tag) {
          int32_t const _M0L5_2aokS3576 = _tmp_4187.data.ok;
          _handle__error__result_4188 = _M0L5_2aokS3576;
        } else {
          void* const _M0L6_2aerrS3577 = _tmp_4187.data.err;
          moonbit_decref(_M0L17error__to__stringS1489);
          moonbit_decref(_M0L13handle__startS1477);
          _M0L11_2atry__errS1504 = _M0L6_2aerrS3577;
          goto join_1503;
        }
        if (_handle__error__result_4188) {
          moonbit_decref(_M0L17error__to__stringS1489);
          moonbit_decref(_M0L13handle__startS1477);
          _M0L6_2atmpS3573 = 1;
        } else {
          struct moonbit_result_0 _tmp_4189;
          #line 633 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
          _tmp_4189
          = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(_M0L12async__testsS1510, _M0L8filenameS1479, _M0L5indexS1481, _M0L13handle__startS1477, _M0L14handle__resultS1482, _M0L17error__to__stringS1489);
          moonbit_decref(_M0L13handle__startS1477);
          moonbit_decref(_M0L17error__to__stringS1489);
          if (_tmp_4189.tag) {
            int32_t const _M0L5_2aokS3574 = _tmp_4189.data.ok;
            _M0L6_2atmpS3573 = _M0L5_2aokS3574;
          } else {
            void* const _M0L6_2aerrS3575 = _tmp_4189.data.err;
            _M0L11_2atry__errS1504 = _M0L6_2aerrS3575;
            goto join_1503;
          }
        }
      }
    }
  }
  if (!_M0L6_2atmpS3573) {
    void* _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3584 =
      (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
    Moonbit_object_header(_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3584)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 2);
    ((struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3584)->$0
    = (moonbit_string_t)moonbit_string_literal_0.data;
    _M0L11_2atry__errS1504
    = _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3584;
    goto join_1503;
  } else {
    moonbit_decref(_M0L14handle__resultS1482);
  }
  goto joinlet_4180;
  join_1503:;
  _M0L3errS1505 = _M0L11_2atry__errS1504;
  _M0L36_2aMoonBitTestDriverInternalSkipTestS1508
  = (struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L3errS1505;
  _M0L8_2afieldS3605 = _M0L36_2aMoonBitTestDriverInternalSkipTestS1508->$0;
  _M0L6_2acntS3921
  = Moonbit_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1508));
  if (_M0L6_2acntS3921 > 1) {
    int32_t _M0L11_2anew__cntS3922 = _M0L6_2acntS3921 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1508), _M0L11_2anew__cntS3922);
    moonbit_incref(_M0L8_2afieldS3605);
  } else if (_M0L6_2acntS3921 == 1) {
    #line 640 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
    moonbit_free(_M0L36_2aMoonBitTestDriverInternalSkipTestS1508);
  }
  _M0L7_2anameS1509 = _M0L8_2afieldS3605;
  _M0L4nameS1507 = _M0L7_2anameS1509;
  goto join_1506;
  goto joinlet_4190;
  join_1506:;
  #line 641 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1482(_M0L14handle__resultS1482, _M0L4nameS1507, (moonbit_string_t)moonbit_string_literal_1.data, 1);
  moonbit_decref(_M0L14handle__resultS1482);
  moonbit_decref(_M0L4nameS1507);
  joinlet_4190:;
  joinlet_4180:;
  return 0;
}

moonbit_string_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1489(
  struct _M0TWRPC15error5ErrorEs* _M0L6_2aenvS3572,
  void* _M0L3errS1490
) {
  void* _M0L1eS1492;
  moonbit_string_t _M0L1eS1494;
  moonbit_string_t _result_4193;
  #line 610 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  switch (Moonbit_object_tag(_M0L3errS1490)) {
    case 0: {
      struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS1495 =
        (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L3errS1490;
      moonbit_string_t _M0L4_2aeS1496 = _M0L10_2aFailureS1495->$0;
      moonbit_incref(_M0L4_2aeS1496);
      _M0L1eS1494 = _M0L4_2aeS1496;
      goto join_1493;
      break;
    }
    
    case 3: {
      struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError* _M0L15_2aInspectErrorS1497 =
        (struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError*)_M0L3errS1490;
      moonbit_string_t _M0L4_2aeS1498 = _M0L15_2aInspectErrorS1497->$0;
      moonbit_incref(_M0L4_2aeS1498);
      _M0L1eS1494 = _M0L4_2aeS1498;
      goto join_1493;
      break;
    }
    
    case 4: {
      struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError* _M0L16_2aSnapshotErrorS1499 =
        (struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError*)_M0L3errS1490;
      moonbit_string_t _M0L4_2aeS1500 = _M0L16_2aSnapshotErrorS1499->$0;
      moonbit_incref(_M0L4_2aeS1500);
      _M0L1eS1494 = _M0L4_2aeS1500;
      goto join_1493;
      break;
    }
    
    case 5: {
      struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError* _M0L35_2aMoonBitTestDriverInternalJsErrorS1501 =
        (struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError*)_M0L3errS1490;
      moonbit_string_t _M0L4_2aeS1502 =
        _M0L35_2aMoonBitTestDriverInternalJsErrorS1501->$0;
      moonbit_incref(_M0L4_2aeS1502);
      _M0L1eS1494 = _M0L4_2aeS1502;
      goto join_1493;
      break;
    }
    default: {
      moonbit_incref(_M0L3errS1490);
      _M0L1eS1492 = _M0L3errS1490;
      goto join_1491;
      break;
    }
  }
  join_1493:;
  return _M0L1eS1494;
  join_1491:;
  #line 616 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _result_4193 = _M0FP15Error10to__string(_M0L1eS1492);
  moonbit_decref(_M0L1eS1492);
  return _result_4193;
}

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1482(
  struct _M0TWssbEu* _M0L6_2aenvS3569,
  moonbit_string_t _M0L10__testnameS1483,
  moonbit_string_t _M0L7messageS1484,
  int32_t _M0L7skippedS1485
) {
  struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482* _M0L14_2acasted__envS3570;
  moonbit_string_t _M0L8filenameS1479;
  int32_t _M0L5indexS1481;
  int32_t _if__result_4194;
  moonbit_string_t _M0L10file__nameS1486;
  moonbit_string_t _M0L7messageS1487;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1488;
  moonbit_string_t _M0L6_2atmpS3571;
  #line 595 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L14_2acasted__envS3570
  = (struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1482*)_M0L6_2aenvS3569;
  _M0L8filenameS1479 = _M0L14_2acasted__envS3570->$1;
  _M0L5indexS1481 = _M0L14_2acasted__envS3570->$0;
  if (!_M0L7skippedS1485) {
    _if__result_4194 = 1;
  } else {
    _if__result_4194 = 0;
  }
  if (_if__result_4194) {
    
  }
  #line 601 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L10file__nameS1486
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1479, 1);
  #line 602 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L7messageS1487
  = _M0MPC16string6String14escape_2einner(_M0L7messageS1484, 1);
  #line 603 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L18_2astring__builderS1488
  = _M0MPB13StringBuilder21StringBuilder_2einner(45);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1488, (moonbit_string_t)moonbit_string_literal_3.data);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1488, _M0L10file__nameS1486);
  moonbit_decref(_M0L10file__nameS1486);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1488, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1488, _M0L5indexS1481);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1488, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1488, _M0L7messageS1487);
  moonbit_decref(_M0L7messageS1487);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1488, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 605 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS3571
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1488);
  moonbit_decref(_M0L18_2astring__builderS1488);
  #line 604 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS3571);
  moonbit_decref(_M0L6_2atmpS3571);
  #line 607 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN13handle__startS1477(
  struct _M0TWEOc* _M0L6_2aenvS3566
) {
  struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477* _M0L14_2acasted__envS3567;
  moonbit_string_t _M0L8filenameS1479;
  int32_t _M0L5indexS1481;
  moonbit_string_t _M0L10file__nameS1478;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1480;
  moonbit_string_t _M0L6_2atmpS3568;
  #line 586 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L14_2acasted__envS3567
  = (struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1477*)_M0L6_2aenvS3566;
  _M0L8filenameS1479 = _M0L14_2acasted__envS3567->$1;
  _M0L5indexS1481 = _M0L14_2acasted__envS3567->$0;
  #line 587 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L10file__nameS1478
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1479, 1);
  #line 588 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L18_2astring__builderS1480
  = _M0MPB13StringBuilder21StringBuilder_2einner(33);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1480, (moonbit_string_t)moonbit_string_literal_8.data);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1480, _M0L10file__nameS1478);
  moonbit_decref(_M0L10file__nameS1478);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1480, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1480, _M0L5indexS1481);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1480, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 590 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS3568
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1480);
  moonbit_decref(_M0L18_2astring__builderS1480);
  #line 589 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS3568);
  moonbit_decref(_M0L6_2atmpS3568);
  #line 592 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

struct moonbit_result_0 _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__test(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1476,
  moonbit_string_t _M0L8filenameS1473,
  int32_t _M0L5indexS1467,
  struct _M0TWEOc* _M0L13handle__startS1462,
  struct _M0TWssbEu* _M0L14handle__resultS1463,
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1465
) {
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L10index__mapS1442;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS1472;
  struct _M0TWEuQRPC15error5Error* _M0L1fS1444;
  moonbit_string_t* _M0L5attrsS1445;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L7_2abindS1466;
  moonbit_string_t _M0L4nameS1448;
  moonbit_string_t _M0L4nameS1446;
  int32_t _M0L6_2atmpS3565;
  struct _M0TPB4IterGsE* _M0L5_2aitS1450;
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__* _closure_4203;
  struct _M0TWEOc* _M0L6_2atmpS3556;
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__* _closure_4204;
  struct _M0TWRPC15error5ErrorEu* _M0L6_2atmpS3557;
  struct moonbit_result_0 _result_4205;
  #line 428 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  #line 436 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L7_2abindS1472
  = _M0MPB3Map3getGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0FP26biolab8bio__seq48moonbit__test__driver__internal__no__args__tests, _M0L8filenameS1473);
  if (_M0L7_2abindS1472 == 0) {
    struct moonbit_result_0 _result_4196;
    if (_M0L7_2abindS1472) {
      moonbit_decref(_M0L7_2abindS1472);
    }
    _result_4196.tag = 1;
    _result_4196.data.ok = 0;
    return _result_4196;
  } else {
    struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS1474 =
      _M0L7_2abindS1472;
    struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L13_2aindex__mapS1475 =
      _M0L7_2aSomeS1474;
    _M0L10index__mapS1442 = _M0L13_2aindex__mapS1475;
    goto join_1441;
  }
  join_1441:;
  #line 438 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L7_2abindS1466
  = _M0MPB3Map3getGiUWEuQRPC15error5ErrorNsEE(_M0L10index__mapS1442, _M0L5indexS1467);
  moonbit_decref(_M0L10index__mapS1442);
  if (_M0L7_2abindS1466 == 0) {
    struct moonbit_result_0 _result_4198;
    if (_M0L7_2abindS1466) {
      moonbit_decref(_M0L7_2abindS1466);
    }
    _result_4198.tag = 1;
    _result_4198.data.ok = 0;
    return _result_4198;
  } else {
    struct _M0TUWEuQRPC15error5ErrorNsE* _M0L7_2aSomeS1468 =
      _M0L7_2abindS1466;
    struct _M0TUWEuQRPC15error5ErrorNsE* _M0L4_2axS1469 = _M0L7_2aSomeS1468;
    struct _M0TWEuQRPC15error5Error* _M0L4_2afS1470 = _M0L4_2axS1469->$0;
    moonbit_string_t* _M0L8_2afieldS3607 = _M0L4_2axS1469->$1;
    int32_t _M0L6_2acntS3923 =
      Moonbit_rc_count(Moonbit_object_header(_M0L4_2axS1469));
    moonbit_string_t* _M0L8_2aattrsS1471;
    if (_M0L6_2acntS3923 > 1) {
      int32_t _M0L11_2anew__cntS3924 = _M0L6_2acntS3923 - 1;
      Moonbit_set_rc_count(Moonbit_object_header(_M0L4_2axS1469), _M0L11_2anew__cntS3924);
      moonbit_incref(_M0L8_2afieldS3607);
      moonbit_incref(_M0L4_2afS1470);
    } else if (_M0L6_2acntS3923 == 1) {
      #line 436 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      moonbit_free(_M0L4_2axS1469);
    }
    _M0L8_2aattrsS1471 = _M0L8_2afieldS3607;
    _M0L1fS1444 = _M0L4_2afS1470;
    _M0L5attrsS1445 = _M0L8_2aattrsS1471;
    goto join_1443;
  }
  join_1443:;
  _M0L6_2atmpS3565 = Moonbit_array_length(_M0L5attrsS1445);
  if (_M0L6_2atmpS3565 >= 1) {
    moonbit_string_t _M0L7_2anameS1449 = (moonbit_string_t)_M0L5attrsS1445[0];
    moonbit_incref(_M0L7_2anameS1449);
    _M0L4nameS1448 = _M0L7_2anameS1449;
    goto join_1447;
  } else {
    _M0L4nameS1446 = (moonbit_string_t)moonbit_string_literal_0.data;
  }
  goto joinlet_4199;
  join_1447:;
  _M0L4nameS1446 = _M0L4nameS1448;
  joinlet_4199:;
  #line 439 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L5_2aitS1450 = _M0MPC15array13ReadOnlyArray4iterGsE(_M0L5attrsS1445);
  moonbit_decref(_M0L5attrsS1445);
  while (1) {
    moonbit_string_t _M0L4attrS1452;
    moonbit_string_t _M0L7_2abindS1459;
    int32_t _M0L6_2atmpS3549;
    #line 441 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
    _M0L7_2abindS1459 = _M0MPB4Iter4nextGsE(_M0L5_2aitS1450);
    if (_M0L7_2abindS1459 == 0) {
      if (_M0L7_2abindS1459) {
        moonbit_decref(_M0L7_2abindS1459);
      }
      moonbit_decref(_M0L5_2aitS1450);
    } else {
      moonbit_string_t _M0L7_2aSomeS1460 = _M0L7_2abindS1459;
      moonbit_string_t _M0L7_2aattrS1461 = _M0L7_2aSomeS1460;
      _M0L4attrS1452 = _M0L7_2aattrS1461;
      goto join_1451;
    }
    goto joinlet_4201;
    join_1451:;
    _M0L6_2atmpS3549 = Moonbit_array_length(_M0L4attrS1452);
    if (_M0L6_2atmpS3549 >= 5) {
      int32_t _M0L6_2atmpS3555 = _M0L4attrS1452[0];
      int32_t _M0L4_2axS1453 = _M0L6_2atmpS3555;
      if (_M0L4_2axS1453 == 112) {
        int32_t _M0L6_2atmpS3554 = _M0L4attrS1452[1];
        int32_t _M0L4_2axS1454 = _M0L6_2atmpS3554;
        if (_M0L4_2axS1454 == 97) {
          int32_t _M0L6_2atmpS3553 = _M0L4attrS1452[2];
          int32_t _M0L4_2axS1455 = _M0L6_2atmpS3553;
          if (_M0L4_2axS1455 == 110) {
            int32_t _M0L6_2atmpS3552 = _M0L4attrS1452[3];
            int32_t _M0L4_2axS1456 = _M0L6_2atmpS3552;
            if (_M0L4_2axS1456 == 105) {
              int32_t _M0L6_2atmpS3551 = _M0L4attrS1452[4];
              int32_t _M0L4_2axS1457;
              moonbit_decref(_M0L4attrS1452);
              _M0L4_2axS1457 = _M0L6_2atmpS3551;
              if (_M0L4_2axS1457 == 99) {
                void* _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3550;
                struct moonbit_result_0 _result_4202;
                moonbit_decref(_M0L5_2aitS1450);
                moonbit_decref(_M0L1fS1444);
                _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3550
                = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
                Moonbit_object_header(_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3550)->meta
                = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 2);
                ((struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3550)->$0
                = _M0L4nameS1446;
                _result_4202.tag = 0;
                _result_4202.data.err
                = _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3550;
                return _result_4202;
              }
            } else {
              moonbit_decref(_M0L4attrS1452);
            }
          } else {
            moonbit_decref(_M0L4attrS1452);
          }
        } else {
          moonbit_decref(_M0L4attrS1452);
        }
      } else {
        moonbit_decref(_M0L4attrS1452);
      }
    } else {
      moonbit_decref(_M0L4attrS1452);
    }
    continue;
    joinlet_4201:;
    break;
  }
  #line 447 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L13handle__startS1462->code(_M0L13handle__startS1462);
  moonbit_incref(_M0L14handle__resultS1463);
  moonbit_incref(_M0L4nameS1446);
  _closure_4203
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__*)moonbit_malloc(sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__));
  Moonbit_object_header(_closure_4203)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 9, 0);
  _closure_4203->code
  = &_M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3562l450;
  _closure_4203->$0 = _M0L14handle__resultS1463;
  _closure_4203->$1 = _M0L4nameS1446;
  _M0L6_2atmpS3556 = (struct _M0TWEOc*)_closure_4203;
  moonbit_incref(_M0L17error__to__stringS1465);
  moonbit_incref(_M0L14handle__resultS1463);
  _closure_4204
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__*)moonbit_malloc(sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__));
  Moonbit_object_header(_closure_4204)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 13, 0);
  _closure_4204->code
  = &_M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3558l451;
  _closure_4204->$0 = _M0L17error__to__stringS1465;
  _closure_4204->$1 = _M0L14handle__resultS1463;
  _closure_4204->$2 = _M0L4nameS1446;
  _M0L6_2atmpS3557 = (struct _M0TWRPC15error5ErrorEu*)_closure_4204;
  #line 448 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FP26biolab8bio__seq45moonbit__test__driver__internal__catch__error(_M0L1fS1444, _M0L6_2atmpS3556, _M0L6_2atmpS3557);
  moonbit_decref(_M0L1fS1444);
  moonbit_decref(_M0L6_2atmpS3556);
  moonbit_decref(_M0L6_2atmpS3557);
  _result_4205.tag = 1;
  _result_4205.data.ok = 1;
  return _result_4205;
}

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3562l450(
  struct _M0TWEOc* _M0L6_2aenvS3563
) {
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__* _M0L14_2acasted__envS3564;
  moonbit_string_t _M0L4nameS1446;
  struct _M0TWssbEu* _M0L14handle__resultS1463;
  #line 450 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L14_2acasted__envS3564
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3562__l450__*)_M0L6_2aenvS3563;
  _M0L4nameS1446 = _M0L14_2acasted__envS3564->$1;
  _M0L14handle__resultS1463 = _M0L14_2acasted__envS3564->$0;
  #line 450 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L14handle__resultS1463->code(_M0L14handle__resultS1463, _M0L4nameS1446, (moonbit_string_t)moonbit_string_literal_0.data, 0);
  return 0;
}

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3558l451(
  struct _M0TWRPC15error5ErrorEu* _M0L6_2aenvS3559,
  void* _M0L3errS1464
) {
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__* _M0L14_2acasted__envS3560;
  moonbit_string_t _M0L4nameS1446;
  struct _M0TWssbEu* _M0L14handle__resultS1463;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1465;
  moonbit_string_t _M0L6_2atmpS3561;
  #line 451 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L14_2acasted__envS3560
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3558__l451__*)_M0L6_2aenvS3559;
  _M0L4nameS1446 = _M0L14_2acasted__envS3560->$2;
  _M0L14handle__resultS1463 = _M0L14_2acasted__envS3560->$1;
  _M0L17error__to__stringS1465 = _M0L14_2acasted__envS3560->$0;
  #line 451 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS3561
  = _M0L17error__to__stringS1465->code(_M0L17error__to__stringS1465, _M0L3errS1464);
  #line 451 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L14handle__resultS1463->code(_M0L14handle__resultS1463, _M0L4nameS1446, _M0L6_2atmpS3561, 0);
  moonbit_decref(_M0L6_2atmpS3561);
  return 0;
}

int32_t _M0FP26biolab8bio__seq45moonbit__test__driver__internal__catch__error(
  struct _M0TWEuQRPC15error5Error* _M0L1fS1436,
  struct _M0TWEOc* _M0L6on__okS1437,
  struct _M0TWRPC15error5ErrorEu* _M0L7on__errS1434
) {
  void* _M0L11_2atry__errS1432;
  struct moonbit_result_0 _tmp_4207;
  void* _M0L3errS1433;
  #line 377 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  #line 384 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _tmp_4207 = _M0L1fS1436->code(_M0L1fS1436);
  if (_tmp_4207.tag) {
    int32_t const _M0L5_2aokS3547 = _tmp_4207.data.ok;
  } else {
    void* const _M0L6_2aerrS3548 = _tmp_4207.data.err;
    _M0L11_2atry__errS1432 = _M0L6_2aerrS3548;
    goto join_1431;
  }
  #line 384 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6on__okS1437->code(_M0L6on__okS1437);
  goto joinlet_4206;
  join_1431:;
  _M0L3errS1433 = _M0L11_2atry__errS1432;
  #line 385 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L7on__errS1434->code(_M0L7on__errS1434, _M0L3errS1433);
  moonbit_decref(_M0L3errS1433);
  joinlet_4206:;
  return 0;
}

struct _M0TPB5ArrayGUsiEE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__args(
  
) {
  int32_t _M0L45moonbit__test__driver__internal__parse__int__S1388;
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395;
  int32_t _M0L57moonbit__test__driver__internal__get__cli__args__internalS1402;
  int32_t _M0L51moonbit__test__driver__internal__split__mbt__stringS1409;
  struct _M0TUsiE** _M0L6_2atmpS3546;
  struct _M0TPB5ArrayGUsiEE* _M0L16file__and__indexS1416;
  struct _M0TPB5ArrayGsE* _M0L9cli__argsS1417;
  moonbit_string_t _M0L6_2atmpS3545;
  struct _M0TPB5ArrayGsE* _M0L10test__argsS1418;
  int32_t _M0L7_2abindS1419;
  int32_t _M0L2__S1420;
  #line 194 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L45moonbit__test__driver__internal__parse__int__S1388 = 0;
  _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395
  = 0;
  _M0L57moonbit__test__driver__internal__get__cli__args__internalS1402
  = _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395;
  _M0L51moonbit__test__driver__internal__split__mbt__stringS1409 = 0;
  _M0L6_2atmpS3546 = (struct _M0TUsiE**)moonbit_empty_ref_array;
  _M0L16file__and__indexS1416
  = (struct _M0TPB5ArrayGUsiEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGUsiEE));
  Moonbit_object_header(_M0L16file__and__indexS1416)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L16file__and__indexS1416->$0 = _M0L6_2atmpS3546;
  _M0L16file__and__indexS1416->$1 = 0;
  #line 283 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L9cli__argsS1417
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1402(_M0L57moonbit__test__driver__internal__get__cli__args__internalS1402);
  #line 285 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS3545 = _M0MPC15array5Array2atGsE(_M0L9cli__argsS1417, 1);
  moonbit_decref(_M0L9cli__argsS1417);
  #line 284 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L10test__argsS1418
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1409(_M0L51moonbit__test__driver__internal__split__mbt__stringS1409, _M0L6_2atmpS3545, 47);
  moonbit_decref(_M0L6_2atmpS3545);
  _M0L7_2abindS1419 = _M0L10test__argsS1418->$1;
  _M0L2__S1420 = 0;
  while (1) {
    if (_M0L2__S1420 < _M0L7_2abindS1419) {
      moonbit_string_t* _M0L3bufS3544 = _M0L10test__argsS1418->$0;
      moonbit_string_t _M0L3argS1421 =
        (moonbit_string_t)_M0L3bufS3544[_M0L2__S1420];
      struct _M0TPB5ArrayGsE* _M0L16file__and__rangeS1422;
      moonbit_string_t _M0L4fileS1423;
      moonbit_string_t _M0L5rangeS1424;
      struct _M0TPB5ArrayGsE* _M0L15start__and__endS1425;
      moonbit_string_t _M0L6_2atmpS3542;
      int32_t _M0L5startS1426;
      moonbit_string_t _M0L6_2atmpS3541;
      int32_t _M0L3endS1427;
      int32_t _M0L1iS1428;
      int32_t _M0L6_2atmpS3543;
      moonbit_incref(_M0L3argS1421);
      #line 289 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L16file__and__rangeS1422
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1409(_M0L51moonbit__test__driver__internal__split__mbt__stringS1409, _M0L3argS1421, 58);
      moonbit_decref(_M0L3argS1421);
      #line 290 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L4fileS1423
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1422, 0);
      #line 291 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L5rangeS1424
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1422, 1);
      moonbit_decref(_M0L16file__and__rangeS1422);
      #line 292 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L15start__and__endS1425
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1409(_M0L51moonbit__test__driver__internal__split__mbt__stringS1409, _M0L5rangeS1424, 45);
      moonbit_decref(_M0L5rangeS1424);
      #line 295 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L6_2atmpS3542
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1425, 0);
      #line 295 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L5startS1426
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1388(_M0L45moonbit__test__driver__internal__parse__int__S1388, _M0L6_2atmpS3542);
      moonbit_decref(_M0L6_2atmpS3542);
      #line 296 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L6_2atmpS3541
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1425, 1);
      moonbit_decref(_M0L15start__and__endS1425);
      #line 296 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L3endS1427
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1388(_M0L45moonbit__test__driver__internal__parse__int__S1388, _M0L6_2atmpS3541);
      moonbit_decref(_M0L6_2atmpS3541);
      _M0L1iS1428 = _M0L5startS1426;
      while (1) {
        if (_M0L1iS1428 < _M0L3endS1427) {
          struct _M0TUsiE* _M0L8_2atupleS3539;
          int32_t _M0L6_2atmpS3540;
          moonbit_incref(_M0L4fileS1423);
          _M0L8_2atupleS3539
          = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
          Moonbit_object_header(_M0L8_2atupleS3539)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
          _M0L8_2atupleS3539->$0 = _M0L4fileS1423;
          _M0L8_2atupleS3539->$1 = _M0L1iS1428;
          #line 298 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
          _M0MPC15array5Array4pushGUsiEE(_M0L16file__and__indexS1416, _M0L8_2atupleS3539);
          moonbit_decref(_M0L8_2atupleS3539);
          _M0L6_2atmpS3540 = _M0L1iS1428 + 1;
          _M0L1iS1428 = _M0L6_2atmpS3540;
          continue;
        } else {
          moonbit_decref(_M0L4fileS1423);
        }
        break;
      }
      _M0L6_2atmpS3543 = _M0L2__S1420 + 1;
      _M0L2__S1420 = _M0L6_2atmpS3543;
      continue;
    } else {
      moonbit_decref(_M0L10test__argsS1418);
    }
    break;
  }
  return _M0L16file__and__indexS1416;
}

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1409(
  int32_t _M0L6_2aenvS3520,
  moonbit_string_t _M0L1sS1410,
  int32_t _M0L3sepS1411
) {
  moonbit_string_t* _M0L6_2atmpS3538;
  struct _M0TPB5ArrayGsE* _M0L3resS1412;
  struct _M0TPB8MutLocalGiE* _M0L1iS1413;
  struct _M0TPB8MutLocalGiE* _M0L5startS1414;
  int32_t _M0L3valS3533;
  int32_t _M0L6_2atmpS3534;
  #line 262 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS3538 = (moonbit_string_t*)moonbit_empty_ref_array;
  _M0L3resS1412
  = (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
  Moonbit_object_header(_M0L3resS1412)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
  _M0L3resS1412->$0 = _M0L6_2atmpS3538;
  _M0L3resS1412->$1 = 0;
  _M0L1iS1413
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1413)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1413->$0 = 0;
  _M0L5startS1414
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS1414)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS1414->$0 = 0;
  while (1) {
    int32_t _M0L3valS3521 = _M0L1iS1413->$0;
    int32_t _M0L6_2atmpS3522 = Moonbit_array_length(_M0L1sS1410);
    if (_M0L3valS3521 < _M0L6_2atmpS3522) {
      int32_t _M0L3valS3525 = _M0L1iS1413->$0;
      int32_t _M0L6_2atmpS3524;
      int32_t _M0L6_2atmpS3523;
      int32_t _M0L3valS3532;
      int32_t _M0L6_2atmpS3531;
      if (
        _M0L3valS3525 < 0
        || _M0L3valS3525 >= Moonbit_array_length(_M0L1sS1410)
      ) {
        #line 270 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3524 = _M0L1sS1410[_M0L3valS3525];
      _M0L6_2atmpS3523 = _M0L6_2atmpS3524;
      if (_M0L6_2atmpS3523 == _M0L3sepS1411) {
        int32_t _M0L3valS3527 = _M0L5startS1414->$0;
        int32_t _M0L3valS3528 = _M0L1iS1413->$0;
        moonbit_string_t _M0L6_2atmpS3526;
        int32_t _M0L3valS3530;
        int32_t _M0L6_2atmpS3529;
        #line 271 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        _M0L6_2atmpS3526
        = _M0MPC16string6String17unsafe__substring(_M0L1sS1410, _M0L3valS3527, _M0L3valS3528);
        #line 271 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        _M0MPC15array5Array4pushGsE(_M0L3resS1412, _M0L6_2atmpS3526);
        moonbit_decref(_M0L6_2atmpS3526);
        _M0L3valS3530 = _M0L1iS1413->$0;
        _M0L6_2atmpS3529 = _M0L3valS3530 + 1;
        _M0L5startS1414->$0 = _M0L6_2atmpS3529;
      }
      _M0L3valS3532 = _M0L1iS1413->$0;
      _M0L6_2atmpS3531 = _M0L3valS3532 + 1;
      _M0L1iS1413->$0 = _M0L6_2atmpS3531;
      continue;
    } else {
      moonbit_decref(_M0L1iS1413);
    }
    break;
  }
  _M0L3valS3533 = _M0L5startS1414->$0;
  _M0L6_2atmpS3534 = Moonbit_array_length(_M0L1sS1410);
  if (_M0L3valS3533 < _M0L6_2atmpS3534) {
    int32_t _M0L3valS3536 = _M0L5startS1414->$0;
    int32_t _M0L6_2atmpS3537;
    moonbit_string_t _M0L6_2atmpS3535;
    moonbit_decref(_M0L5startS1414);
    _M0L6_2atmpS3537 = Moonbit_array_length(_M0L1sS1410);
    #line 277 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
    _M0L6_2atmpS3535
    = _M0MPC16string6String17unsafe__substring(_M0L1sS1410, _M0L3valS3536, _M0L6_2atmpS3537);
    #line 277 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
    _M0MPC15array5Array4pushGsE(_M0L3resS1412, _M0L6_2atmpS3535);
    moonbit_decref(_M0L6_2atmpS3535);
  } else {
    moonbit_decref(_M0L5startS1414);
  }
  return _M0L3resS1412;
}

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1402(
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395
) {
  moonbit_bytes_t* _M0L3tmpS1403;
  int32_t _M0L6_2atmpS3519;
  struct _M0TPB5ArrayGsE* _M0L3resS1404;
  int32_t _M0L7_2abindS1405;
  int32_t _M0L7_2abindS1406;
  int32_t _M0L1iS1407;
  #line 251 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  #line 254 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L3tmpS1403
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__get__cli__args__ffi();
  _M0L6_2atmpS3519 = Moonbit_array_length(_M0L3tmpS1403);
  #line 255 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L3resS1404 = _M0MPC15array5Array11new_2einnerGsE(_M0L6_2atmpS3519);
  _M0L7_2abindS1405 = 0;
  _M0L7_2abindS1406 = Moonbit_array_length(_M0L3tmpS1403);
  _M0L1iS1407 = _M0L7_2abindS1405;
  while (1) {
    if (_M0L1iS1407 < _M0L7_2abindS1406) {
      moonbit_bytes_t _M0L6_2atmpS3517;
      moonbit_string_t _M0L6_2atmpS3516;
      int32_t _M0L6_2atmpS3518;
      if (
        _M0L1iS1407 < 0 || _M0L1iS1407 >= Moonbit_array_length(_M0L3tmpS1403)
      ) {
        #line 257 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3517 = (moonbit_bytes_t)_M0L3tmpS1403[_M0L1iS1407];
      moonbit_incref(_M0L6_2atmpS3517);
      #line 257 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0L6_2atmpS3516
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395(_M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395, _M0L6_2atmpS3517);
      moonbit_decref(_M0L6_2atmpS3517);
      #line 257 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0MPC15array5Array4pushGsE(_M0L3resS1404, _M0L6_2atmpS3516);
      moonbit_decref(_M0L6_2atmpS3516);
      _M0L6_2atmpS3518 = _M0L1iS1407 + 1;
      _M0L1iS1407 = _M0L6_2atmpS3518;
      continue;
    } else {
      moonbit_decref(_M0L3tmpS1403);
    }
    break;
  }
  return _M0L3resS1404;
}

moonbit_string_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1395(
  int32_t _M0L6_2aenvS3430,
  moonbit_bytes_t _M0L5bytesS1396
) {
  struct _M0TPB13StringBuilder* _M0L3resS1397;
  int32_t _M0L3lenS1398;
  struct _M0TPB8MutLocalGiE* _M0L1iS1399;
  moonbit_string_t _result_4213;
  #line 207 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  #line 210 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L3resS1397 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3lenS1398 = Moonbit_array_length(_M0L5bytesS1396);
  _M0L1iS1399
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1399)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1399->$0 = 0;
  while (1) {
    int32_t _M0L3valS3431 = _M0L1iS1399->$0;
    if (_M0L3valS3431 < _M0L3lenS1398) {
      int32_t _M0L3valS3515 = _M0L1iS1399->$0;
      int32_t _M0L6_2atmpS3514;
      int32_t _M0L6_2atmpS3513;
      struct _M0TPB8MutLocalGiE* _M0L1cS1400;
      int32_t _M0L3valS3432;
      if (
        _M0L3valS3515 < 0
        || _M0L3valS3515 >= Moonbit_array_length(_M0L5bytesS1396)
      ) {
        #line 214 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3514 = _M0L5bytesS1396[_M0L3valS3515];
      _M0L6_2atmpS3513 = (int32_t)_M0L6_2atmpS3514;
      _M0L1cS1400
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1cS1400)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1cS1400->$0 = _M0L6_2atmpS3513;
      _M0L3valS3432 = _M0L1cS1400->$0;
      if (_M0L3valS3432 < 128) {
        int32_t _M0L3valS3434 = _M0L1cS1400->$0;
        int32_t _M0L6_2atmpS3433;
        int32_t _M0L3valS3436;
        int32_t _M0L6_2atmpS3435;
        moonbit_decref(_M0L1cS1400);
        _M0L6_2atmpS3433 = _M0L3valS3434;
        #line 216 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1397, _M0L6_2atmpS3433);
        _M0L3valS3436 = _M0L1iS1399->$0;
        _M0L6_2atmpS3435 = _M0L3valS3436 + 1;
        _M0L1iS1399->$0 = _M0L6_2atmpS3435;
      } else {
        int32_t _M0L3valS3437 = _M0L1cS1400->$0;
        if (_M0L3valS3437 < 224) {
          int32_t _M0L3valS3439 = _M0L1iS1399->$0;
          int32_t _M0L6_2atmpS3438 = _M0L3valS3439 + 1;
          int32_t _M0L3valS3448;
          int32_t _M0L6_2atmpS3447;
          int32_t _M0L6_2atmpS3441;
          int32_t _M0L3valS3446;
          int32_t _M0L6_2atmpS3445;
          int32_t _M0L6_2atmpS3444;
          int32_t _M0L6_2atmpS3443;
          int32_t _M0L6_2atmpS3442;
          int32_t _M0L6_2atmpS3440;
          int32_t _M0L3valS3450;
          int32_t _M0L6_2atmpS3449;
          int32_t _M0L3valS3452;
          int32_t _M0L6_2atmpS3451;
          if (_M0L6_2atmpS3438 >= _M0L3lenS1398) {
            moonbit_decref(_M0L1cS1400);
            moonbit_decref(_M0L1iS1399);
            break;
          }
          _M0L3valS3448 = _M0L1cS1400->$0;
          _M0L6_2atmpS3447 = _M0L3valS3448 & 31;
          _M0L6_2atmpS3441 = _M0L6_2atmpS3447 << 6;
          _M0L3valS3446 = _M0L1iS1399->$0;
          _M0L6_2atmpS3445 = _M0L3valS3446 + 1;
          if (
            _M0L6_2atmpS3445 < 0
            || _M0L6_2atmpS3445 >= Moonbit_array_length(_M0L5bytesS1396)
          ) {
            #line 222 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS3444 = _M0L5bytesS1396[_M0L6_2atmpS3445];
          _M0L6_2atmpS3443 = (int32_t)_M0L6_2atmpS3444;
          _M0L6_2atmpS3442 = _M0L6_2atmpS3443 & 63;
          _M0L6_2atmpS3440 = _M0L6_2atmpS3441 | _M0L6_2atmpS3442;
          _M0L1cS1400->$0 = _M0L6_2atmpS3440;
          _M0L3valS3450 = _M0L1cS1400->$0;
          moonbit_decref(_M0L1cS1400);
          _M0L6_2atmpS3449 = _M0L3valS3450;
          #line 223 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1397, _M0L6_2atmpS3449);
          _M0L3valS3452 = _M0L1iS1399->$0;
          _M0L6_2atmpS3451 = _M0L3valS3452 + 2;
          _M0L1iS1399->$0 = _M0L6_2atmpS3451;
        } else {
          int32_t _M0L3valS3453 = _M0L1cS1400->$0;
          if (_M0L3valS3453 < 240) {
            int32_t _M0L3valS3455 = _M0L1iS1399->$0;
            int32_t _M0L6_2atmpS3454 = _M0L3valS3455 + 2;
            int32_t _M0L3valS3471;
            int32_t _M0L6_2atmpS3470;
            int32_t _M0L6_2atmpS3463;
            int32_t _M0L3valS3469;
            int32_t _M0L6_2atmpS3468;
            int32_t _M0L6_2atmpS3467;
            int32_t _M0L6_2atmpS3466;
            int32_t _M0L6_2atmpS3465;
            int32_t _M0L6_2atmpS3464;
            int32_t _M0L6_2atmpS3457;
            int32_t _M0L3valS3462;
            int32_t _M0L6_2atmpS3461;
            int32_t _M0L6_2atmpS3460;
            int32_t _M0L6_2atmpS3459;
            int32_t _M0L6_2atmpS3458;
            int32_t _M0L6_2atmpS3456;
            int32_t _M0L3valS3473;
            int32_t _M0L6_2atmpS3472;
            int32_t _M0L3valS3475;
            int32_t _M0L6_2atmpS3474;
            if (_M0L6_2atmpS3454 >= _M0L3lenS1398) {
              moonbit_decref(_M0L1cS1400);
              moonbit_decref(_M0L1iS1399);
              break;
            }
            _M0L3valS3471 = _M0L1cS1400->$0;
            _M0L6_2atmpS3470 = _M0L3valS3471 & 15;
            _M0L6_2atmpS3463 = _M0L6_2atmpS3470 << 12;
            _M0L3valS3469 = _M0L1iS1399->$0;
            _M0L6_2atmpS3468 = _M0L3valS3469 + 1;
            if (
              _M0L6_2atmpS3468 < 0
              || _M0L6_2atmpS3468 >= Moonbit_array_length(_M0L5bytesS1396)
            ) {
              #line 230 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3467 = _M0L5bytesS1396[_M0L6_2atmpS3468];
            _M0L6_2atmpS3466 = (int32_t)_M0L6_2atmpS3467;
            _M0L6_2atmpS3465 = _M0L6_2atmpS3466 & 63;
            _M0L6_2atmpS3464 = _M0L6_2atmpS3465 << 6;
            _M0L6_2atmpS3457 = _M0L6_2atmpS3463 | _M0L6_2atmpS3464;
            _M0L3valS3462 = _M0L1iS1399->$0;
            _M0L6_2atmpS3461 = _M0L3valS3462 + 2;
            if (
              _M0L6_2atmpS3461 < 0
              || _M0L6_2atmpS3461 >= Moonbit_array_length(_M0L5bytesS1396)
            ) {
              #line 231 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3460 = _M0L5bytesS1396[_M0L6_2atmpS3461];
            _M0L6_2atmpS3459 = (int32_t)_M0L6_2atmpS3460;
            _M0L6_2atmpS3458 = _M0L6_2atmpS3459 & 63;
            _M0L6_2atmpS3456 = _M0L6_2atmpS3457 | _M0L6_2atmpS3458;
            _M0L1cS1400->$0 = _M0L6_2atmpS3456;
            _M0L3valS3473 = _M0L1cS1400->$0;
            moonbit_decref(_M0L1cS1400);
            _M0L6_2atmpS3472 = _M0L3valS3473;
            #line 232 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1397, _M0L6_2atmpS3472);
            _M0L3valS3475 = _M0L1iS1399->$0;
            _M0L6_2atmpS3474 = _M0L3valS3475 + 3;
            _M0L1iS1399->$0 = _M0L6_2atmpS3474;
          } else {
            int32_t _M0L3valS3477 = _M0L1iS1399->$0;
            int32_t _M0L6_2atmpS3476 = _M0L3valS3477 + 3;
            int32_t _M0L3valS3500;
            int32_t _M0L6_2atmpS3499;
            int32_t _M0L6_2atmpS3492;
            int32_t _M0L3valS3498;
            int32_t _M0L6_2atmpS3497;
            int32_t _M0L6_2atmpS3496;
            int32_t _M0L6_2atmpS3495;
            int32_t _M0L6_2atmpS3494;
            int32_t _M0L6_2atmpS3493;
            int32_t _M0L6_2atmpS3485;
            int32_t _M0L3valS3491;
            int32_t _M0L6_2atmpS3490;
            int32_t _M0L6_2atmpS3489;
            int32_t _M0L6_2atmpS3488;
            int32_t _M0L6_2atmpS3487;
            int32_t _M0L6_2atmpS3486;
            int32_t _M0L6_2atmpS3479;
            int32_t _M0L3valS3484;
            int32_t _M0L6_2atmpS3483;
            int32_t _M0L6_2atmpS3482;
            int32_t _M0L6_2atmpS3481;
            int32_t _M0L6_2atmpS3480;
            int32_t _M0L6_2atmpS3478;
            int32_t _M0L3valS3502;
            int32_t _M0L6_2atmpS3501;
            int32_t _M0L3valS3506;
            int32_t _M0L6_2atmpS3505;
            int32_t _M0L6_2atmpS3504;
            int32_t _M0L6_2atmpS3503;
            int32_t _M0L3valS3510;
            int32_t _M0L6_2atmpS3509;
            int32_t _M0L6_2atmpS3508;
            int32_t _M0L6_2atmpS3507;
            int32_t _M0L3valS3512;
            int32_t _M0L6_2atmpS3511;
            if (_M0L6_2atmpS3476 >= _M0L3lenS1398) {
              moonbit_decref(_M0L1cS1400);
              moonbit_decref(_M0L1iS1399);
              break;
            }
            _M0L3valS3500 = _M0L1cS1400->$0;
            _M0L6_2atmpS3499 = _M0L3valS3500 & 7;
            _M0L6_2atmpS3492 = _M0L6_2atmpS3499 << 18;
            _M0L3valS3498 = _M0L1iS1399->$0;
            _M0L6_2atmpS3497 = _M0L3valS3498 + 1;
            if (
              _M0L6_2atmpS3497 < 0
              || _M0L6_2atmpS3497 >= Moonbit_array_length(_M0L5bytesS1396)
            ) {
              #line 239 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3496 = _M0L5bytesS1396[_M0L6_2atmpS3497];
            _M0L6_2atmpS3495 = (int32_t)_M0L6_2atmpS3496;
            _M0L6_2atmpS3494 = _M0L6_2atmpS3495 & 63;
            _M0L6_2atmpS3493 = _M0L6_2atmpS3494 << 12;
            _M0L6_2atmpS3485 = _M0L6_2atmpS3492 | _M0L6_2atmpS3493;
            _M0L3valS3491 = _M0L1iS1399->$0;
            _M0L6_2atmpS3490 = _M0L3valS3491 + 2;
            if (
              _M0L6_2atmpS3490 < 0
              || _M0L6_2atmpS3490 >= Moonbit_array_length(_M0L5bytesS1396)
            ) {
              #line 240 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3489 = _M0L5bytesS1396[_M0L6_2atmpS3490];
            _M0L6_2atmpS3488 = (int32_t)_M0L6_2atmpS3489;
            _M0L6_2atmpS3487 = _M0L6_2atmpS3488 & 63;
            _M0L6_2atmpS3486 = _M0L6_2atmpS3487 << 6;
            _M0L6_2atmpS3479 = _M0L6_2atmpS3485 | _M0L6_2atmpS3486;
            _M0L3valS3484 = _M0L1iS1399->$0;
            _M0L6_2atmpS3483 = _M0L3valS3484 + 3;
            if (
              _M0L6_2atmpS3483 < 0
              || _M0L6_2atmpS3483 >= Moonbit_array_length(_M0L5bytesS1396)
            ) {
              #line 241 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3482 = _M0L5bytesS1396[_M0L6_2atmpS3483];
            _M0L6_2atmpS3481 = (int32_t)_M0L6_2atmpS3482;
            _M0L6_2atmpS3480 = _M0L6_2atmpS3481 & 63;
            _M0L6_2atmpS3478 = _M0L6_2atmpS3479 | _M0L6_2atmpS3480;
            _M0L1cS1400->$0 = _M0L6_2atmpS3478;
            _M0L3valS3502 = _M0L1cS1400->$0;
            _M0L6_2atmpS3501 = _M0L3valS3502 - 65536;
            _M0L1cS1400->$0 = _M0L6_2atmpS3501;
            _M0L3valS3506 = _M0L1cS1400->$0;
            _M0L6_2atmpS3505 = _M0L3valS3506 >> 10;
            _M0L6_2atmpS3504 = _M0L6_2atmpS3505 + 55296;
            _M0L6_2atmpS3503 = _M0L6_2atmpS3504;
            #line 243 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1397, _M0L6_2atmpS3503);
            _M0L3valS3510 = _M0L1cS1400->$0;
            moonbit_decref(_M0L1cS1400);
            _M0L6_2atmpS3509 = _M0L3valS3510 & 1023;
            _M0L6_2atmpS3508 = _M0L6_2atmpS3509 + 56320;
            _M0L6_2atmpS3507 = _M0L6_2atmpS3508;
            #line 244 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1397, _M0L6_2atmpS3507);
            _M0L3valS3512 = _M0L1iS1399->$0;
            _M0L6_2atmpS3511 = _M0L3valS3512 + 4;
            _M0L1iS1399->$0 = _M0L6_2atmpS3511;
          }
        }
      }
      continue;
    } else {
      moonbit_decref(_M0L1iS1399);
    }
    break;
  }
  #line 248 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _result_4213 = _M0MPB13StringBuilder10to__string(_M0L3resS1397);
  moonbit_decref(_M0L3resS1397);
  return _result_4213;
}

int32_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1388(
  int32_t _M0L6_2aenvS3423,
  moonbit_string_t _M0L1sS1389
) {
  struct _M0TPB8MutLocalGiE* _M0L3resS1390;
  int32_t _M0L3lenS1391;
  int32_t _M0L7_2abindS1392;
  int32_t _M0L1iS1393;
  int32_t _result_4215;
  #line 198 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L3resS1390
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3resS1390)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3resS1390->$0 = 0;
  _M0L3lenS1391 = Moonbit_array_length(_M0L1sS1389);
  _M0L7_2abindS1392 = 0;
  _M0L1iS1393 = _M0L7_2abindS1392;
  while (1) {
    if (_M0L1iS1393 < _M0L3lenS1391) {
      int32_t _M0L3valS3428 = _M0L3resS1390->$0;
      int32_t _M0L6_2atmpS3425 = _M0L3valS3428 * 10;
      int32_t _M0L6_2atmpS3427;
      int32_t _M0L6_2atmpS3426;
      int32_t _M0L6_2atmpS3424;
      int32_t _M0L6_2atmpS3429;
      if (
        _M0L1iS1393 < 0 || _M0L1iS1393 >= Moonbit_array_length(_M0L1sS1389)
      ) {
        #line 202 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3427 = _M0L1sS1389[_M0L1iS1393];
      _M0L6_2atmpS3426 = _M0L6_2atmpS3427 - 48;
      _M0L6_2atmpS3424 = _M0L6_2atmpS3425 + _M0L6_2atmpS3426;
      _M0L3resS1390->$0 = _M0L6_2atmpS3424;
      _M0L6_2atmpS3429 = _M0L1iS1393 + 1;
      _M0L1iS1393 = _M0L6_2atmpS3429;
      continue;
    }
    break;
  }
  _result_4215 = _M0L3resS1390->$0;
  moonbit_decref(_M0L3resS1390);
  return _result_4215;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1364,
  moonbit_string_t _M0L12_2adiscard__S1365,
  int32_t _M0L12_2adiscard__S1366,
  struct _M0TWEOc* _M0L12_2adiscard__S1367,
  struct _M0TWssbEu* _M0L12_2adiscard__S1368,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1369
) {
  struct moonbit_result_0 _result_4216;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _result_4216.tag = 1;
  _result_4216.data.ok = 0;
  return _result_4216;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1370,
  moonbit_string_t _M0L12_2adiscard__S1371,
  int32_t _M0L12_2adiscard__S1372,
  struct _M0TWEOc* _M0L12_2adiscard__S1373,
  struct _M0TWssbEu* _M0L12_2adiscard__S1374,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1375
) {
  struct moonbit_result_0 _result_4217;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _result_4217.tag = 1;
  _result_4217.data.ok = 0;
  return _result_4217;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1376,
  moonbit_string_t _M0L12_2adiscard__S1377,
  int32_t _M0L12_2adiscard__S1378,
  struct _M0TWEOc* _M0L12_2adiscard__S1379,
  struct _M0TWssbEu* _M0L12_2adiscard__S1380,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1381
) {
  struct moonbit_result_0 _result_4218;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _result_4218.tag = 1;
  _result_4218.data.ok = 0;
  return _result_4218;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1382,
  moonbit_string_t _M0L12_2adiscard__S1383,
  int32_t _M0L12_2adiscard__S1384,
  struct _M0TWEOc* _M0L12_2adiscard__S1385,
  struct _M0TWssbEu* _M0L12_2adiscard__S1386,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1387
) {
  struct moonbit_result_0 _result_4219;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _result_4219.tag = 1;
  _result_4219.data.ok = 0;
  return _result_4219;
}

int32_t _M0IP016_24default__implP26biolab8bio__seq28MoonBit__Async__Test__Driver17run__async__testsGRP26biolab8bio__seq34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1363
) {
  #line 12 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  return 0;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__19(
  
) {
  moonbit_string_t _M0L8originalS1353;
  struct moonbit_result_1 _tmp_4220;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1354;
  struct moonbit_result_2 _tmp_4222;
  moonbit_string_t _M0L9rewrittenS1355;
  struct moonbit_result_1 _tmp_4224;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L8records2S1356;
  int32_t _M0L6_2atmpS3390;
  moonbit_string_t _M0L6_2atmpS3391;
  struct moonbit_result_0 _tmp_4226;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3396;
  moonbit_string_t _M0L8_2afieldS3614;
  int32_t _M0L6_2acntS3925;
  moonbit_string_t _M0L2idS3394;
  moonbit_string_t _M0L6_2atmpS3395;
  struct moonbit_result_0 _tmp_4228;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3401;
  moonbit_string_t _M0L8_2afieldS3613;
  int32_t _M0L6_2acntS3931;
  moonbit_string_t _M0L11descriptionS3399;
  moonbit_string_t _M0L6_2atmpS3400;
  struct moonbit_result_0 _tmp_4230;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3407;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3612;
  int32_t _M0L6_2acntS3937;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3406;
  moonbit_string_t _M0L6_2atmpS3404;
  moonbit_string_t _M0L6_2atmpS3405;
  struct moonbit_result_0 _tmp_4232;
  struct _M0TPB5ArrayGiE* _M0L1qS1358;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3416;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1359;
  int32_t* _M0L6_2atmpS3413;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS3412;
  int32_t _M0L6_2atmpS3410;
  void* _M0L4NoneS3411;
  struct moonbit_result_0 _result_4236;
  #line 236 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L8originalS1353 = (moonbit_string_t)moonbit_string_literal_9.data;
  #line 238 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4220 = _M0FP26biolab8bio__seq12parse__fastq(_M0L8originalS1353);
  moonbit_decref(_M0L8originalS1353);
  if (_tmp_4220.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS3421 =
      _tmp_4220.data.ok;
    _M0L7recordsS1354 = _M0L5_2aokS3421;
  } else {
    void* const _M0L6_2aerrS3422 = _tmp_4220.data.err;
    struct moonbit_result_0 _result_4221;
    _result_4221.tag = 0;
    _result_4221.data.err = _M0L6_2aerrS3422;
    return _result_4221;
  }
  #line 239 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4222 = _M0FP26biolab8bio__seq12write__fastq(_M0L7recordsS1354);
  moonbit_decref(_M0L7recordsS1354);
  if (_tmp_4222.tag) {
    moonbit_string_t const _M0L5_2aokS3419 = _tmp_4222.data.ok;
    _M0L9rewrittenS1355 = _M0L5_2aokS3419;
  } else {
    void* const _M0L6_2aerrS3420 = _tmp_4222.data.err;
    struct moonbit_result_0 _result_4223;
    _result_4223.tag = 0;
    _result_4223.data.err = _M0L6_2aerrS3420;
    return _result_4223;
  }
  #line 240 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4224 = _M0FP26biolab8bio__seq12parse__fastq(_M0L9rewrittenS1355);
  moonbit_decref(_M0L9rewrittenS1355);
  if (_tmp_4224.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS3417 =
      _tmp_4224.data.ok;
    _M0L8records2S1356 = _M0L5_2aokS3417;
  } else {
    void* const _M0L6_2aerrS3418 = _tmp_4224.data.err;
    struct moonbit_result_0 _result_4225;
    _result_4225.tag = 0;
    _result_4225.data.err = _M0L6_2aerrS3418;
    return _result_4225;
  }
  #line 241 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3390
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1356);
  _M0L6_2atmpS3391 = 0;
  #line 241 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4226
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3390, 1, _M0L6_2atmpS3391, (moonbit_string_t)moonbit_string_literal_10.data);
  if (_M0L6_2atmpS3391) {
    moonbit_decref(_M0L6_2atmpS3391);
  }
  if (_tmp_4226.tag) {
    int32_t const _M0L5_2aokS3392 = _tmp_4226.data.ok;
  } else {
    void* const _M0L6_2aerrS3393 = _tmp_4226.data.err;
    struct moonbit_result_0 _result_4227;
    moonbit_decref(_M0L8records2S1356);
    _result_4227.tag = 0;
    _result_4227.data.err = _M0L6_2aerrS3393;
    return _result_4227;
  }
  #line 242 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3396
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1356, 0);
  _M0L8_2afieldS3614 = _M0L6_2atmpS3396->$1;
  _M0L6_2acntS3925
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3396));
  if (_M0L6_2acntS3925 > 1) {
    int32_t _M0L11_2anew__cntS3930 = _M0L6_2acntS3925 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3396), _M0L11_2anew__cntS3930);
    moonbit_incref(_M0L8_2afieldS3614);
  } else if (_M0L6_2acntS3925 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3929 =
      _M0L6_2atmpS3396->$4;
    moonbit_string_t _M0L8_2afieldS3928;
    moonbit_string_t _M0L8_2afieldS3927;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3926;
    moonbit_decref(_M0L8_2afieldS3929);
    _M0L8_2afieldS3928 = _M0L6_2atmpS3396->$3;
    moonbit_decref(_M0L8_2afieldS3928);
    _M0L8_2afieldS3927 = _M0L6_2atmpS3396->$2;
    moonbit_decref(_M0L8_2afieldS3927);
    _M0L8_2afieldS3926 = _M0L6_2atmpS3396->$0;
    moonbit_decref(_M0L8_2afieldS3926);
    #line 242 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3396);
  }
  _M0L2idS3394 = _M0L8_2afieldS3614;
  _M0L6_2atmpS3395 = 0;
  #line 242 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4228
  = _M0FPB10assert__eqGsE(_M0L2idS3394, (moonbit_string_t)moonbit_string_literal_11.data, _M0L6_2atmpS3395, (moonbit_string_t)moonbit_string_literal_12.data);
  moonbit_decref(_M0L2idS3394);
  if (_M0L6_2atmpS3395) {
    moonbit_decref(_M0L6_2atmpS3395);
  }
  if (_tmp_4228.tag) {
    int32_t const _M0L5_2aokS3397 = _tmp_4228.data.ok;
  } else {
    void* const _M0L6_2aerrS3398 = _tmp_4228.data.err;
    struct moonbit_result_0 _result_4229;
    moonbit_decref(_M0L8records2S1356);
    _result_4229.tag = 0;
    _result_4229.data.err = _M0L6_2aerrS3398;
    return _result_4229;
  }
  #line 243 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3401
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1356, 0);
  _M0L8_2afieldS3613 = _M0L6_2atmpS3401->$3;
  _M0L6_2acntS3931
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3401));
  if (_M0L6_2acntS3931 > 1) {
    int32_t _M0L11_2anew__cntS3936 = _M0L6_2acntS3931 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3401), _M0L11_2anew__cntS3936);
    moonbit_incref(_M0L8_2afieldS3613);
  } else if (_M0L6_2acntS3931 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3935 =
      _M0L6_2atmpS3401->$4;
    moonbit_string_t _M0L8_2afieldS3934;
    moonbit_string_t _M0L8_2afieldS3933;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3932;
    moonbit_decref(_M0L8_2afieldS3935);
    _M0L8_2afieldS3934 = _M0L6_2atmpS3401->$2;
    moonbit_decref(_M0L8_2afieldS3934);
    _M0L8_2afieldS3933 = _M0L6_2atmpS3401->$1;
    moonbit_decref(_M0L8_2afieldS3933);
    _M0L8_2afieldS3932 = _M0L6_2atmpS3401->$0;
    moonbit_decref(_M0L8_2afieldS3932);
    #line 243 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3401);
  }
  _M0L11descriptionS3399 = _M0L8_2afieldS3613;
  _M0L6_2atmpS3400 = 0;
  #line 243 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4230
  = _M0FPB10assert__eqGsE(_M0L11descriptionS3399, (moonbit_string_t)moonbit_string_literal_13.data, _M0L6_2atmpS3400, (moonbit_string_t)moonbit_string_literal_14.data);
  moonbit_decref(_M0L11descriptionS3399);
  if (_M0L6_2atmpS3400) {
    moonbit_decref(_M0L6_2atmpS3400);
  }
  if (_tmp_4230.tag) {
    int32_t const _M0L5_2aokS3402 = _tmp_4230.data.ok;
  } else {
    void* const _M0L6_2aerrS3403 = _tmp_4230.data.err;
    struct moonbit_result_0 _result_4231;
    moonbit_decref(_M0L8records2S1356);
    _result_4231.tag = 0;
    _result_4231.data.err = _M0L6_2aerrS3403;
    return _result_4231;
  }
  #line 244 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3407
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1356, 0);
  _M0L8_2afieldS3612 = _M0L6_2atmpS3407->$0;
  _M0L6_2acntS3937
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3407));
  if (_M0L6_2acntS3937 > 1) {
    int32_t _M0L11_2anew__cntS3942 = _M0L6_2acntS3937 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3407), _M0L11_2anew__cntS3942);
    moonbit_incref(_M0L8_2afieldS3612);
  } else if (_M0L6_2acntS3937 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3941 =
      _M0L6_2atmpS3407->$4;
    moonbit_string_t _M0L8_2afieldS3940;
    moonbit_string_t _M0L8_2afieldS3939;
    moonbit_string_t _M0L8_2afieldS3938;
    moonbit_decref(_M0L8_2afieldS3941);
    _M0L8_2afieldS3940 = _M0L6_2atmpS3407->$3;
    moonbit_decref(_M0L8_2afieldS3940);
    _M0L8_2afieldS3939 = _M0L6_2atmpS3407->$2;
    moonbit_decref(_M0L8_2afieldS3939);
    _M0L8_2afieldS3938 = _M0L6_2atmpS3407->$1;
    moonbit_decref(_M0L8_2afieldS3938);
    #line 244 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3407);
  }
  _M0L3seqS3406 = _M0L8_2afieldS3612;
  #line 244 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3404 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3406);
  moonbit_decref(_M0L3seqS3406);
  _M0L6_2atmpS3405 = 0;
  #line 244 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4232
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3404, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3405, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L6_2atmpS3404);
  if (_M0L6_2atmpS3405) {
    moonbit_decref(_M0L6_2atmpS3405);
  }
  if (_tmp_4232.tag) {
    int32_t const _M0L5_2aokS3408 = _tmp_4232.data.ok;
  } else {
    void* const _M0L6_2aerrS3409 = _tmp_4232.data.err;
    struct moonbit_result_0 _result_4233;
    moonbit_decref(_M0L8records2S1356);
    _result_4233.tag = 0;
    _result_4233.data.err = _M0L6_2aerrS3409;
    return _result_4233;
  }
  #line 245 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3416
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1356, 0);
  moonbit_decref(_M0L8records2S1356);
  #line 245 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1359
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L6_2atmpS3416);
  moonbit_decref(_M0L6_2atmpS3416);
  if (_M0L7_2abindS1359 == 0) {
    moonbit_string_t _M0L7_2abindS1362;
    int32_t _M0L6_2atmpS3415;
    struct _M0TPC16string10StringView _M0L6_2atmpS3414;
    struct moonbit_result_0 _result_4235;
    if (_M0L7_2abindS1359) {
      moonbit_decref(_M0L7_2abindS1359);
    }
    _M0L7_2abindS1362 = (moonbit_string_t)moonbit_string_literal_17.data;
    _M0L6_2atmpS3415 = Moonbit_array_length(_M0L7_2abindS1362);
    _M0L6_2atmpS3414
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1362, .$1 = 0, .$2 = _M0L6_2atmpS3415
    };
    #line 247 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _result_4235
    = _M0FPB4failGuE(_M0L6_2atmpS3414, (moonbit_string_t)moonbit_string_literal_18.data);
    moonbit_decref(_M0L6_2atmpS3414.$0);
    return _result_4235;
  } else {
    struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1360 = _M0L7_2abindS1359;
    struct _M0TPB5ArrayGiE* _M0L4_2aqS1361 = _M0L7_2aSomeS1360;
    _M0L1qS1358 = _M0L4_2aqS1361;
    goto join_1357;
  }
  join_1357:;
  _M0L6_2atmpS3413 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS3413[0] = 40;
  _M0L6_2atmpS3413[1] = 40;
  _M0L6_2atmpS3413[2] = 40;
  _M0L6_2atmpS3413[3] = 40;
  _M0L6_2atmpS3412
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS3412)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS3412->$0 = _M0L6_2atmpS3413;
  _M0L6_2atmpS3412->$1 = 4;
  #line 246 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3410
  = _M0IPC15array5ArrayPB2Eq5equalGiE(_M0L1qS1358, _M0L6_2atmpS3412);
  moonbit_decref(_M0L1qS1358);
  moonbit_decref(_M0L6_2atmpS3412);
  _M0L4NoneS3411
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 246 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4236
  = _M0FPB12assert__true(_M0L6_2atmpS3410, _M0L4NoneS3411, (moonbit_string_t)moonbit_string_literal_19.data);
  moonbit_decref(_M0L4NoneS3411);
  return _result_4236;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__18(
  
) {
  moonbit_string_t _M0L8originalS1349;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1350;
  moonbit_string_t _M0L9rewrittenS1351;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L8records2S1352;
  int32_t _M0L6_2atmpS3361;
  moonbit_string_t _M0L6_2atmpS3362;
  struct moonbit_result_0 _tmp_4237;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3367;
  moonbit_string_t _M0L8_2afieldS3619;
  int32_t _M0L6_2acntS3943;
  moonbit_string_t _M0L2idS3365;
  moonbit_string_t _M0L6_2atmpS3366;
  struct moonbit_result_0 _tmp_4239;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3372;
  moonbit_string_t _M0L8_2afieldS3618;
  int32_t _M0L6_2acntS3949;
  moonbit_string_t _M0L11descriptionS3370;
  moonbit_string_t _M0L6_2atmpS3371;
  struct moonbit_result_0 _tmp_4241;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3378;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3617;
  int32_t _M0L6_2acntS3955;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3377;
  moonbit_string_t _M0L6_2atmpS3375;
  moonbit_string_t _M0L6_2atmpS3376;
  struct moonbit_result_0 _tmp_4243;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3383;
  moonbit_string_t _M0L8_2afieldS3616;
  int32_t _M0L6_2acntS3961;
  moonbit_string_t _M0L2idS3381;
  moonbit_string_t _M0L6_2atmpS3382;
  struct moonbit_result_0 _tmp_4245;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3389;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3615;
  int32_t _M0L6_2acntS3967;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3388;
  moonbit_string_t _M0L6_2atmpS3386;
  moonbit_string_t _M0L6_2atmpS3387;
  struct moonbit_result_0 _result_4247;
  #line 222 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L8originalS1349 = (moonbit_string_t)moonbit_string_literal_20.data;
  #line 224 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7recordsS1350
  = _M0FP26biolab8bio__seq12parse__fasta(_M0L8originalS1349);
  moonbit_decref(_M0L8originalS1349);
  #line 225 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L9rewrittenS1351
  = _M0FP26biolab8bio__seq12write__fasta(_M0L7recordsS1350);
  moonbit_decref(_M0L7recordsS1350);
  #line 226 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L8records2S1352
  = _M0FP26biolab8bio__seq12parse__fasta(_M0L9rewrittenS1351);
  moonbit_decref(_M0L9rewrittenS1351);
  #line 227 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3361
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1352);
  _M0L6_2atmpS3362 = 0;
  #line 227 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4237
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3361, 2, _M0L6_2atmpS3362, (moonbit_string_t)moonbit_string_literal_21.data);
  if (_M0L6_2atmpS3362) {
    moonbit_decref(_M0L6_2atmpS3362);
  }
  if (_tmp_4237.tag) {
    int32_t const _M0L5_2aokS3363 = _tmp_4237.data.ok;
  } else {
    void* const _M0L6_2aerrS3364 = _tmp_4237.data.err;
    struct moonbit_result_0 _result_4238;
    moonbit_decref(_M0L8records2S1352);
    _result_4238.tag = 0;
    _result_4238.data.err = _M0L6_2aerrS3364;
    return _result_4238;
  }
  #line 228 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3367
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1352, 0);
  _M0L8_2afieldS3619 = _M0L6_2atmpS3367->$1;
  _M0L6_2acntS3943
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3367));
  if (_M0L6_2acntS3943 > 1) {
    int32_t _M0L11_2anew__cntS3948 = _M0L6_2acntS3943 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3367), _M0L11_2anew__cntS3948);
    moonbit_incref(_M0L8_2afieldS3619);
  } else if (_M0L6_2acntS3943 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3947 =
      _M0L6_2atmpS3367->$4;
    moonbit_string_t _M0L8_2afieldS3946;
    moonbit_string_t _M0L8_2afieldS3945;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3944;
    moonbit_decref(_M0L8_2afieldS3947);
    _M0L8_2afieldS3946 = _M0L6_2atmpS3367->$3;
    moonbit_decref(_M0L8_2afieldS3946);
    _M0L8_2afieldS3945 = _M0L6_2atmpS3367->$2;
    moonbit_decref(_M0L8_2afieldS3945);
    _M0L8_2afieldS3944 = _M0L6_2atmpS3367->$0;
    moonbit_decref(_M0L8_2afieldS3944);
    #line 228 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3367);
  }
  _M0L2idS3365 = _M0L8_2afieldS3619;
  _M0L6_2atmpS3366 = 0;
  #line 228 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4239
  = _M0FPB10assert__eqGsE(_M0L2idS3365, (moonbit_string_t)moonbit_string_literal_22.data, _M0L6_2atmpS3366, (moonbit_string_t)moonbit_string_literal_23.data);
  moonbit_decref(_M0L2idS3365);
  if (_M0L6_2atmpS3366) {
    moonbit_decref(_M0L6_2atmpS3366);
  }
  if (_tmp_4239.tag) {
    int32_t const _M0L5_2aokS3368 = _tmp_4239.data.ok;
  } else {
    void* const _M0L6_2aerrS3369 = _tmp_4239.data.err;
    struct moonbit_result_0 _result_4240;
    moonbit_decref(_M0L8records2S1352);
    _result_4240.tag = 0;
    _result_4240.data.err = _M0L6_2aerrS3369;
    return _result_4240;
  }
  #line 229 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3372
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1352, 0);
  _M0L8_2afieldS3618 = _M0L6_2atmpS3372->$3;
  _M0L6_2acntS3949
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3372));
  if (_M0L6_2acntS3949 > 1) {
    int32_t _M0L11_2anew__cntS3954 = _M0L6_2acntS3949 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3372), _M0L11_2anew__cntS3954);
    moonbit_incref(_M0L8_2afieldS3618);
  } else if (_M0L6_2acntS3949 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3953 =
      _M0L6_2atmpS3372->$4;
    moonbit_string_t _M0L8_2afieldS3952;
    moonbit_string_t _M0L8_2afieldS3951;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3950;
    moonbit_decref(_M0L8_2afieldS3953);
    _M0L8_2afieldS3952 = _M0L6_2atmpS3372->$2;
    moonbit_decref(_M0L8_2afieldS3952);
    _M0L8_2afieldS3951 = _M0L6_2atmpS3372->$1;
    moonbit_decref(_M0L8_2afieldS3951);
    _M0L8_2afieldS3950 = _M0L6_2atmpS3372->$0;
    moonbit_decref(_M0L8_2afieldS3950);
    #line 229 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3372);
  }
  _M0L11descriptionS3370 = _M0L8_2afieldS3618;
  _M0L6_2atmpS3371 = 0;
  #line 229 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4241
  = _M0FPB10assert__eqGsE(_M0L11descriptionS3370, (moonbit_string_t)moonbit_string_literal_24.data, _M0L6_2atmpS3371, (moonbit_string_t)moonbit_string_literal_25.data);
  moonbit_decref(_M0L11descriptionS3370);
  if (_M0L6_2atmpS3371) {
    moonbit_decref(_M0L6_2atmpS3371);
  }
  if (_tmp_4241.tag) {
    int32_t const _M0L5_2aokS3373 = _tmp_4241.data.ok;
  } else {
    void* const _M0L6_2aerrS3374 = _tmp_4241.data.err;
    struct moonbit_result_0 _result_4242;
    moonbit_decref(_M0L8records2S1352);
    _result_4242.tag = 0;
    _result_4242.data.err = _M0L6_2aerrS3374;
    return _result_4242;
  }
  #line 230 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3378
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1352, 0);
  _M0L8_2afieldS3617 = _M0L6_2atmpS3378->$0;
  _M0L6_2acntS3955
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3378));
  if (_M0L6_2acntS3955 > 1) {
    int32_t _M0L11_2anew__cntS3960 = _M0L6_2acntS3955 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3378), _M0L11_2anew__cntS3960);
    moonbit_incref(_M0L8_2afieldS3617);
  } else if (_M0L6_2acntS3955 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3959 =
      _M0L6_2atmpS3378->$4;
    moonbit_string_t _M0L8_2afieldS3958;
    moonbit_string_t _M0L8_2afieldS3957;
    moonbit_string_t _M0L8_2afieldS3956;
    moonbit_decref(_M0L8_2afieldS3959);
    _M0L8_2afieldS3958 = _M0L6_2atmpS3378->$3;
    moonbit_decref(_M0L8_2afieldS3958);
    _M0L8_2afieldS3957 = _M0L6_2atmpS3378->$2;
    moonbit_decref(_M0L8_2afieldS3957);
    _M0L8_2afieldS3956 = _M0L6_2atmpS3378->$1;
    moonbit_decref(_M0L8_2afieldS3956);
    #line 230 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3378);
  }
  _M0L3seqS3377 = _M0L8_2afieldS3617;
  #line 230 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3375 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3377);
  moonbit_decref(_M0L3seqS3377);
  _M0L6_2atmpS3376 = 0;
  #line 230 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4243
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3375, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3376, (moonbit_string_t)moonbit_string_literal_26.data);
  moonbit_decref(_M0L6_2atmpS3375);
  if (_M0L6_2atmpS3376) {
    moonbit_decref(_M0L6_2atmpS3376);
  }
  if (_tmp_4243.tag) {
    int32_t const _M0L5_2aokS3379 = _tmp_4243.data.ok;
  } else {
    void* const _M0L6_2aerrS3380 = _tmp_4243.data.err;
    struct moonbit_result_0 _result_4244;
    moonbit_decref(_M0L8records2S1352);
    _result_4244.tag = 0;
    _result_4244.data.err = _M0L6_2aerrS3380;
    return _result_4244;
  }
  #line 231 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3383
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1352, 1);
  _M0L8_2afieldS3616 = _M0L6_2atmpS3383->$1;
  _M0L6_2acntS3961
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3383));
  if (_M0L6_2acntS3961 > 1) {
    int32_t _M0L11_2anew__cntS3966 = _M0L6_2acntS3961 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3383), _M0L11_2anew__cntS3966);
    moonbit_incref(_M0L8_2afieldS3616);
  } else if (_M0L6_2acntS3961 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3965 =
      _M0L6_2atmpS3383->$4;
    moonbit_string_t _M0L8_2afieldS3964;
    moonbit_string_t _M0L8_2afieldS3963;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3962;
    moonbit_decref(_M0L8_2afieldS3965);
    _M0L8_2afieldS3964 = _M0L6_2atmpS3383->$3;
    moonbit_decref(_M0L8_2afieldS3964);
    _M0L8_2afieldS3963 = _M0L6_2atmpS3383->$2;
    moonbit_decref(_M0L8_2afieldS3963);
    _M0L8_2afieldS3962 = _M0L6_2atmpS3383->$0;
    moonbit_decref(_M0L8_2afieldS3962);
    #line 231 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3383);
  }
  _M0L2idS3381 = _M0L8_2afieldS3616;
  _M0L6_2atmpS3382 = 0;
  #line 231 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4245
  = _M0FPB10assert__eqGsE(_M0L2idS3381, (moonbit_string_t)moonbit_string_literal_27.data, _M0L6_2atmpS3382, (moonbit_string_t)moonbit_string_literal_28.data);
  moonbit_decref(_M0L2idS3381);
  if (_M0L6_2atmpS3382) {
    moonbit_decref(_M0L6_2atmpS3382);
  }
  if (_tmp_4245.tag) {
    int32_t const _M0L5_2aokS3384 = _tmp_4245.data.ok;
  } else {
    void* const _M0L6_2aerrS3385 = _tmp_4245.data.err;
    struct moonbit_result_0 _result_4246;
    moonbit_decref(_M0L8records2S1352);
    _result_4246.tag = 0;
    _result_4246.data.err = _M0L6_2aerrS3385;
    return _result_4246;
  }
  #line 232 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3389
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L8records2S1352, 1);
  moonbit_decref(_M0L8records2S1352);
  _M0L8_2afieldS3615 = _M0L6_2atmpS3389->$0;
  _M0L6_2acntS3967
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3389));
  if (_M0L6_2acntS3967 > 1) {
    int32_t _M0L11_2anew__cntS3972 = _M0L6_2acntS3967 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3389), _M0L11_2anew__cntS3972);
    moonbit_incref(_M0L8_2afieldS3615);
  } else if (_M0L6_2acntS3967 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3971 =
      _M0L6_2atmpS3389->$4;
    moonbit_string_t _M0L8_2afieldS3970;
    moonbit_string_t _M0L8_2afieldS3969;
    moonbit_string_t _M0L8_2afieldS3968;
    moonbit_decref(_M0L8_2afieldS3971);
    _M0L8_2afieldS3970 = _M0L6_2atmpS3389->$3;
    moonbit_decref(_M0L8_2afieldS3970);
    _M0L8_2afieldS3969 = _M0L6_2atmpS3389->$2;
    moonbit_decref(_M0L8_2afieldS3969);
    _M0L8_2afieldS3968 = _M0L6_2atmpS3389->$1;
    moonbit_decref(_M0L8_2afieldS3968);
    #line 232 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3389);
  }
  _M0L3seqS3388 = _M0L8_2afieldS3615;
  #line 232 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3386 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3388);
  moonbit_decref(_M0L3seqS3388);
  _M0L6_2atmpS3387 = 0;
  #line 232 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4247
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3386, (moonbit_string_t)moonbit_string_literal_29.data, _M0L6_2atmpS3387, (moonbit_string_t)moonbit_string_literal_30.data);
  moonbit_decref(_M0L6_2atmpS3386);
  if (_M0L6_2atmpS3387) {
    moonbit_decref(_M0L6_2atmpS3387);
  }
  return _result_4247;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__17(
  
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS3360;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1339;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1340;
  int32_t _result_4248;
  int32_t* _M0L6_2atmpS3359;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS3358;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4rec2S1342;
  struct _M0TPB5ArrayGiE* _M0L1qS1344;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1345;
  int32_t* _M0L6_2atmpS3355;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS3354;
  int32_t _M0L6_2atmpS3352;
  void* _M0L4NoneS3353;
  struct moonbit_result_0 _result_4253;
  #line 206 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 207 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3360
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_15.data);
  #line 207 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L3recS1339
  = _M0MP26biolab8bio__seq9SeqRecord9from__seq(_M0L6_2atmpS3360);
  moonbit_decref(_M0L6_2atmpS3360);
  #line 209 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1340
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L3recS1339);
  _result_4248 = _M0L7_2abindS1340 == 0;
  if (_M0L7_2abindS1340) {
    moonbit_decref(_M0L7_2abindS1340);
  }
  if (_result_4248) {
    
  } else {
    moonbit_string_t _M0L7_2abindS1341 =
      (moonbit_string_t)moonbit_string_literal_31.data;
    int32_t _M0L6_2atmpS3349 = Moonbit_array_length(_M0L7_2abindS1341);
    struct _M0TPC16string10StringView _M0L6_2atmpS3348 =
      (struct _M0TPC16string10StringView){.$0 = _M0L7_2abindS1341,
                                            .$1 = 0,
                                            .$2 = _M0L6_2atmpS3349};
    struct moonbit_result_0 _tmp_4249;
    #line 210 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _tmp_4249
    = _M0FPB4failGuE(_M0L6_2atmpS3348, (moonbit_string_t)moonbit_string_literal_32.data);
    moonbit_decref(_M0L6_2atmpS3348.$0);
    if (_tmp_4249.tag) {
      int32_t const _M0L5_2aokS3350 = _tmp_4249.data.ok;
    } else {
      void* const _M0L6_2aerrS3351 = _tmp_4249.data.err;
      struct moonbit_result_0 _result_4250;
      moonbit_decref(_M0L3recS1339);
      _result_4250.tag = 0;
      _result_4250.data.err = _M0L6_2aerrS3351;
      return _result_4250;
    }
  }
  _M0L6_2atmpS3359 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS3359[0] = 40;
  _M0L6_2atmpS3359[1] = 30;
  _M0L6_2atmpS3359[2] = 20;
  _M0L6_2atmpS3359[3] = 10;
  _M0L6_2atmpS3358
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS3358)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS3358->$0 = _M0L6_2atmpS3359;
  _M0L6_2atmpS3358->$1 = 4;
  #line 214 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L4rec2S1342
  = _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(_M0L3recS1339, _M0L6_2atmpS3358);
  moonbit_decref(_M0L3recS1339);
  moonbit_decref(_M0L6_2atmpS3358);
  #line 215 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1345
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L4rec2S1342);
  moonbit_decref(_M0L4rec2S1342);
  if (_M0L7_2abindS1345 == 0) {
    moonbit_string_t _M0L7_2abindS1348;
    int32_t _M0L6_2atmpS3357;
    struct _M0TPC16string10StringView _M0L6_2atmpS3356;
    struct moonbit_result_0 _result_4252;
    if (_M0L7_2abindS1345) {
      moonbit_decref(_M0L7_2abindS1345);
    }
    _M0L7_2abindS1348 = (moonbit_string_t)moonbit_string_literal_33.data;
    _M0L6_2atmpS3357 = Moonbit_array_length(_M0L7_2abindS1348);
    _M0L6_2atmpS3356
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1348, .$1 = 0, .$2 = _M0L6_2atmpS3357
    };
    #line 217 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _result_4252
    = _M0FPB4failGuE(_M0L6_2atmpS3356, (moonbit_string_t)moonbit_string_literal_34.data);
    moonbit_decref(_M0L6_2atmpS3356.$0);
    return _result_4252;
  } else {
    struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1346 = _M0L7_2abindS1345;
    struct _M0TPB5ArrayGiE* _M0L4_2aqS1347 = _M0L7_2aSomeS1346;
    _M0L1qS1344 = _M0L4_2aqS1347;
    goto join_1343;
  }
  join_1343:;
  _M0L6_2atmpS3355 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS3355[0] = 40;
  _M0L6_2atmpS3355[1] = 30;
  _M0L6_2atmpS3355[2] = 20;
  _M0L6_2atmpS3355[3] = 10;
  _M0L6_2atmpS3354
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS3354)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS3354->$0 = _M0L6_2atmpS3355;
  _M0L6_2atmpS3354->$1 = 4;
  #line 216 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3352
  = _M0IPC15array5ArrayPB2Eq5equalGiE(_M0L1qS1344, _M0L6_2atmpS3354);
  moonbit_decref(_M0L1qS1344);
  moonbit_decref(_M0L6_2atmpS3354);
  _M0L4NoneS3353
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 216 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4253
  = _M0FPB12assert__true(_M0L6_2atmpS3352, _M0L4NoneS3353, (moonbit_string_t)moonbit_string_literal_35.data);
  moonbit_decref(_M0L4NoneS3353);
  return _result_4253;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord9from__seq(
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS1337
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1338;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS3347;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS3346;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS3345;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_4254;
  #line 48 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS1338 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS3347 = _M0L7_2abindS1338;
  _M0L6_2atmpS3346
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS3347, .$1 = 0, .$2 = 0
  };
  #line 54 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS3345 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS3346, 0ll);
  moonbit_decref(_M0L6_2atmpS3346.$0);
  moonbit_incref(_M0L3seqS1337);
  _block_4254
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_4254)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
  _block_4254->$0 = _M0L3seqS1337;
  _block_4254->$1 = (moonbit_string_t)moonbit_string_literal_36.data;
  _block_4254->$2 = (moonbit_string_t)moonbit_string_literal_37.data;
  _block_4254->$3 = (moonbit_string_t)moonbit_string_literal_38.data;
  _block_4254->$4 = _M0L6_2atmpS3345;
  return _block_4254;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__16(
  
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS3344;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1336;
  moonbit_string_t _M0L2idS3325;
  moonbit_string_t _M0L6_2atmpS3326;
  struct moonbit_result_0 _tmp_4255;
  moonbit_string_t _M0L4nameS3329;
  moonbit_string_t _M0L6_2atmpS3330;
  struct moonbit_result_0 _tmp_4257;
  moonbit_string_t _M0L11descriptionS3333;
  moonbit_string_t _M0L6_2atmpS3334;
  struct moonbit_result_0 _tmp_4259;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3339;
  moonbit_string_t _M0L6_2atmpS3337;
  moonbit_string_t _M0L6_2atmpS3338;
  struct moonbit_result_0 _tmp_4261;
  int32_t _M0L6_2atmpS3342;
  moonbit_string_t _M0L6_2atmpS3343;
  struct moonbit_result_0 _result_4263;
  #line 191 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 193 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3344
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_15.data);
  #line 192 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L3recS1336
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS3344, (moonbit_string_t)moonbit_string_literal_39.data, (moonbit_string_t)moonbit_string_literal_40.data, (moonbit_string_t)moonbit_string_literal_41.data);
  moonbit_decref(_M0L6_2atmpS3344);
  _M0L2idS3325 = _M0L3recS1336->$1;
  _M0L6_2atmpS3326 = 0;
  moonbit_incref(_M0L2idS3325);
  #line 198 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4255
  = _M0FPB10assert__eqGsE(_M0L2idS3325, (moonbit_string_t)moonbit_string_literal_39.data, _M0L6_2atmpS3326, (moonbit_string_t)moonbit_string_literal_42.data);
  moonbit_decref(_M0L2idS3325);
  if (_M0L6_2atmpS3326) {
    moonbit_decref(_M0L6_2atmpS3326);
  }
  if (_tmp_4255.tag) {
    int32_t const _M0L5_2aokS3327 = _tmp_4255.data.ok;
  } else {
    void* const _M0L6_2aerrS3328 = _tmp_4255.data.err;
    struct moonbit_result_0 _result_4256;
    moonbit_decref(_M0L3recS1336);
    _result_4256.tag = 0;
    _result_4256.data.err = _M0L6_2aerrS3328;
    return _result_4256;
  }
  _M0L4nameS3329 = _M0L3recS1336->$2;
  _M0L6_2atmpS3330 = 0;
  moonbit_incref(_M0L4nameS3329);
  #line 199 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4257
  = _M0FPB10assert__eqGsE(_M0L4nameS3329, (moonbit_string_t)moonbit_string_literal_40.data, _M0L6_2atmpS3330, (moonbit_string_t)moonbit_string_literal_43.data);
  moonbit_decref(_M0L4nameS3329);
  if (_M0L6_2atmpS3330) {
    moonbit_decref(_M0L6_2atmpS3330);
  }
  if (_tmp_4257.tag) {
    int32_t const _M0L5_2aokS3331 = _tmp_4257.data.ok;
  } else {
    void* const _M0L6_2aerrS3332 = _tmp_4257.data.err;
    struct moonbit_result_0 _result_4258;
    moonbit_decref(_M0L3recS1336);
    _result_4258.tag = 0;
    _result_4258.data.err = _M0L6_2aerrS3332;
    return _result_4258;
  }
  _M0L11descriptionS3333 = _M0L3recS1336->$3;
  _M0L6_2atmpS3334 = 0;
  moonbit_incref(_M0L11descriptionS3333);
  #line 200 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4259
  = _M0FPB10assert__eqGsE(_M0L11descriptionS3333, (moonbit_string_t)moonbit_string_literal_41.data, _M0L6_2atmpS3334, (moonbit_string_t)moonbit_string_literal_44.data);
  moonbit_decref(_M0L11descriptionS3333);
  if (_M0L6_2atmpS3334) {
    moonbit_decref(_M0L6_2atmpS3334);
  }
  if (_tmp_4259.tag) {
    int32_t const _M0L5_2aokS3335 = _tmp_4259.data.ok;
  } else {
    void* const _M0L6_2aerrS3336 = _tmp_4259.data.err;
    struct moonbit_result_0 _result_4260;
    moonbit_decref(_M0L3recS1336);
    _result_4260.tag = 0;
    _result_4260.data.err = _M0L6_2aerrS3336;
    return _result_4260;
  }
  _M0L3seqS3339 = _M0L3recS1336->$0;
  moonbit_incref(_M0L3seqS3339);
  #line 201 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3337 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3339);
  moonbit_decref(_M0L3seqS3339);
  _M0L6_2atmpS3338 = 0;
  #line 201 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4261
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3337, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3338, (moonbit_string_t)moonbit_string_literal_45.data);
  moonbit_decref(_M0L6_2atmpS3337);
  if (_M0L6_2atmpS3338) {
    moonbit_decref(_M0L6_2atmpS3338);
  }
  if (_tmp_4261.tag) {
    int32_t const _M0L5_2aokS3340 = _tmp_4261.data.ok;
  } else {
    void* const _M0L6_2aerrS3341 = _tmp_4261.data.err;
    struct moonbit_result_0 _result_4262;
    moonbit_decref(_M0L3recS1336);
    _result_4262.tag = 0;
    _result_4262.data.err = _M0L6_2aerrS3341;
    return _result_4262;
  }
  #line 202 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3342 = _M0MP26biolab8bio__seq9SeqRecord6length(_M0L3recS1336);
  moonbit_decref(_M0L3recS1336);
  _M0L6_2atmpS3343 = 0;
  #line 202 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4263
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3342, 4, _M0L6_2atmpS3343, (moonbit_string_t)moonbit_string_literal_46.data);
  if (_M0L6_2atmpS3343) {
    moonbit_decref(_M0L6_2atmpS3343);
  }
  return _result_4263;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__15(
  
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1322;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L1dS1323;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L1rS1325;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L7_2abindS1326;
  moonbit_string_t _M0L2idS3307;
  moonbit_string_t _M0L6_2atmpS3308;
  struct moonbit_result_0 _tmp_4267;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3625;
  int32_t _M0L6_2acntS3973;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3313;
  moonbit_string_t _M0L6_2atmpS3311;
  moonbit_string_t _M0L6_2atmpS3312;
  struct moonbit_result_0 _tmp_4269;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L1rS1331;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L7_2abindS1332;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3624;
  int32_t _M0L6_2acntS3979;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3322;
  moonbit_string_t _M0L6_2atmpS3320;
  moonbit_string_t _M0L6_2atmpS3321;
  struct moonbit_result_0 _result_4273;
  #line 174 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 175 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7recordsS1322
  = _M0FP26biolab8bio__seq12parse__fasta((moonbit_string_t)moonbit_string_literal_47.data);
  #line 176 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L1dS1323 = _M0FP26biolab8bio__seq15seqio__to__dict(_M0L7recordsS1322);
  moonbit_decref(_M0L7recordsS1322);
  #line 177 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1326
  = _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(_M0L1dS1323, (moonbit_string_t)moonbit_string_literal_48.data);
  if (_M0L7_2abindS1326 == 0) {
    moonbit_string_t _M0L7_2abindS1329;
    int32_t _M0L6_2atmpS3317;
    struct _M0TPC16string10StringView _M0L6_2atmpS3316;
    struct moonbit_result_0 _tmp_4265;
    if (_M0L7_2abindS1326) {
      moonbit_decref(_M0L7_2abindS1326);
    }
    _M0L7_2abindS1329 = (moonbit_string_t)moonbit_string_literal_49.data;
    _M0L6_2atmpS3317 = Moonbit_array_length(_M0L7_2abindS1329);
    _M0L6_2atmpS3316
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1329, .$1 = 0, .$2 = _M0L6_2atmpS3317
    };
    #line 182 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _tmp_4265
    = _M0FPB4failGuE(_M0L6_2atmpS3316, (moonbit_string_t)moonbit_string_literal_50.data);
    moonbit_decref(_M0L6_2atmpS3316.$0);
    if (_tmp_4265.tag) {
      int32_t const _M0L5_2aokS3318 = _tmp_4265.data.ok;
    } else {
      void* const _M0L6_2aerrS3319 = _tmp_4265.data.err;
      struct moonbit_result_0 _result_4266;
      moonbit_decref(_M0L1dS1323);
      _result_4266.tag = 0;
      _result_4266.data.err = _M0L6_2aerrS3319;
      return _result_4266;
    }
  } else {
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L7_2aSomeS1327 =
      _M0L7_2abindS1326;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4_2arS1328 =
      _M0L7_2aSomeS1327;
    _M0L1rS1325 = _M0L4_2arS1328;
    goto join_1324;
  }
  goto joinlet_4264;
  join_1324:;
  _M0L2idS3307 = _M0L1rS1325->$1;
  _M0L6_2atmpS3308 = 0;
  moonbit_incref(_M0L2idS3307);
  #line 179 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4267
  = _M0FPB10assert__eqGsE(_M0L2idS3307, (moonbit_string_t)moonbit_string_literal_48.data, _M0L6_2atmpS3308, (moonbit_string_t)moonbit_string_literal_51.data);
  moonbit_decref(_M0L2idS3307);
  if (_M0L6_2atmpS3308) {
    moonbit_decref(_M0L6_2atmpS3308);
  }
  if (_tmp_4267.tag) {
    int32_t const _M0L5_2aokS3309 = _tmp_4267.data.ok;
  } else {
    void* const _M0L6_2aerrS3310 = _tmp_4267.data.err;
    struct moonbit_result_0 _result_4268;
    moonbit_decref(_M0L1rS1325);
    moonbit_decref(_M0L1dS1323);
    _result_4268.tag = 0;
    _result_4268.data.err = _M0L6_2aerrS3310;
    return _result_4268;
  }
  _M0L8_2afieldS3625 = _M0L1rS1325->$0;
  _M0L6_2acntS3973 = Moonbit_rc_count(Moonbit_object_header(_M0L1rS1325));
  if (_M0L6_2acntS3973 > 1) {
    int32_t _M0L11_2anew__cntS3978 = _M0L6_2acntS3973 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L1rS1325), _M0L11_2anew__cntS3978);
    moonbit_incref(_M0L8_2afieldS3625);
  } else if (_M0L6_2acntS3973 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3977 = _M0L1rS1325->$4;
    moonbit_string_t _M0L8_2afieldS3976;
    moonbit_string_t _M0L8_2afieldS3975;
    moonbit_string_t _M0L8_2afieldS3974;
    moonbit_decref(_M0L8_2afieldS3977);
    _M0L8_2afieldS3976 = _M0L1rS1325->$3;
    moonbit_decref(_M0L8_2afieldS3976);
    _M0L8_2afieldS3975 = _M0L1rS1325->$2;
    moonbit_decref(_M0L8_2afieldS3975);
    _M0L8_2afieldS3974 = _M0L1rS1325->$1;
    moonbit_decref(_M0L8_2afieldS3974);
    #line 180 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L1rS1325);
  }
  _M0L3seqS3313 = _M0L8_2afieldS3625;
  #line 180 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3311 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3313);
  moonbit_decref(_M0L3seqS3313);
  _M0L6_2atmpS3312 = 0;
  #line 180 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4269
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3311, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3312, (moonbit_string_t)moonbit_string_literal_52.data);
  moonbit_decref(_M0L6_2atmpS3311);
  if (_M0L6_2atmpS3312) {
    moonbit_decref(_M0L6_2atmpS3312);
  }
  if (_tmp_4269.tag) {
    int32_t const _M0L5_2aokS3314 = _tmp_4269.data.ok;
  } else {
    void* const _M0L6_2aerrS3315 = _tmp_4269.data.err;
    struct moonbit_result_0 _result_4270;
    moonbit_decref(_M0L1dS1323);
    _result_4270.tag = 0;
    _result_4270.data.err = _M0L6_2aerrS3315;
    return _result_4270;
  }
  joinlet_4264:;
  #line 184 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1332
  = _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(_M0L1dS1323, (moonbit_string_t)moonbit_string_literal_53.data);
  moonbit_decref(_M0L1dS1323);
  if (_M0L7_2abindS1332 == 0) {
    moonbit_string_t _M0L7_2abindS1335;
    int32_t _M0L6_2atmpS3324;
    struct _M0TPC16string10StringView _M0L6_2atmpS3323;
    struct moonbit_result_0 _result_4272;
    if (_M0L7_2abindS1332) {
      moonbit_decref(_M0L7_2abindS1332);
    }
    _M0L7_2abindS1335 = (moonbit_string_t)moonbit_string_literal_54.data;
    _M0L6_2atmpS3324 = Moonbit_array_length(_M0L7_2abindS1335);
    _M0L6_2atmpS3323
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1335, .$1 = 0, .$2 = _M0L6_2atmpS3324
    };
    #line 186 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _result_4272
    = _M0FPB4failGuE(_M0L6_2atmpS3323, (moonbit_string_t)moonbit_string_literal_55.data);
    moonbit_decref(_M0L6_2atmpS3323.$0);
    return _result_4272;
  } else {
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L7_2aSomeS1333 =
      _M0L7_2abindS1332;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4_2arS1334 =
      _M0L7_2aSomeS1333;
    _M0L1rS1331 = _M0L4_2arS1334;
    goto join_1330;
  }
  join_1330:;
  _M0L8_2afieldS3624 = _M0L1rS1331->$0;
  _M0L6_2acntS3979 = Moonbit_rc_count(Moonbit_object_header(_M0L1rS1331));
  if (_M0L6_2acntS3979 > 1) {
    int32_t _M0L11_2anew__cntS3984 = _M0L6_2acntS3979 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L1rS1331), _M0L11_2anew__cntS3984);
    moonbit_incref(_M0L8_2afieldS3624);
  } else if (_M0L6_2acntS3979 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3983 = _M0L1rS1331->$4;
    moonbit_string_t _M0L8_2afieldS3982;
    moonbit_string_t _M0L8_2afieldS3981;
    moonbit_string_t _M0L8_2afieldS3980;
    moonbit_decref(_M0L8_2afieldS3983);
    _M0L8_2afieldS3982 = _M0L1rS1331->$3;
    moonbit_decref(_M0L8_2afieldS3982);
    _M0L8_2afieldS3981 = _M0L1rS1331->$2;
    moonbit_decref(_M0L8_2afieldS3981);
    _M0L8_2afieldS3980 = _M0L1rS1331->$1;
    moonbit_decref(_M0L8_2afieldS3980);
    #line 185 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L1rS1331);
  }
  _M0L3seqS3322 = _M0L8_2afieldS3624;
  #line 185 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3320 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3322);
  moonbit_decref(_M0L3seqS3322);
  _M0L6_2atmpS3321 = 0;
  #line 185 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4273
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3320, (moonbit_string_t)moonbit_string_literal_29.data, _M0L6_2atmpS3321, (moonbit_string_t)moonbit_string_literal_56.data);
  moonbit_decref(_M0L6_2atmpS3320);
  if (_M0L6_2atmpS3321) {
    moonbit_decref(_M0L6_2atmpS3321);
  }
  return _result_4273;
}

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq15seqio__to__dict(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1318
) {
  struct _M0TUsRP26biolab8bio__seq9SeqRecordE** _M0L7_2abindS1316;
  struct _M0TUsRP26biolab8bio__seq9SeqRecordE** _M0L6_2atmpS3306;
  struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE _M0L6_2atmpS3305;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L1dS1315;
  int32_t _M0L7_2abindS1317;
  int32_t _M0L2__S1319;
  #line 67 "/home/developer/Documents/1/seqio.mbt"
  _M0L7_2abindS1316
  = (struct _M0TUsRP26biolab8bio__seq9SeqRecordE**)moonbit_empty_ref_array;
  _M0L6_2atmpS3306 = _M0L7_2abindS1316;
  _M0L6_2atmpS3305
  = (struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE){
    .$0 = _M0L6_2atmpS3306, .$1 = 0, .$2 = 0
  };
  #line 68 "/home/developer/Documents/1/seqio.mbt"
  _M0L1dS1315
  = _M0MPB3Map3MapGsRP26biolab8bio__seq9SeqRecordE(_M0L6_2atmpS3305, 16ll);
  moonbit_decref(_M0L6_2atmpS3305.$0);
  _M0L7_2abindS1317 = _M0L7recordsS1318->$1;
  _M0L2__S1319 = 0;
  while (1) {
    if (_M0L2__S1319 < _M0L7_2abindS1317) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS3304 =
        _M0L7recordsS1318->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1320 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS3304[_M0L2__S1319];
      moonbit_string_t _M0L2idS3302 = _M0L3recS1320->$1;
      int32_t _M0L6_2atmpS3303;
      moonbit_incref(_M0L2idS3302);
      moonbit_incref(_M0L3recS1320);
      #line 70 "/home/developer/Documents/1/seqio.mbt"
      _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(_M0L1dS1315, _M0L2idS3302, _M0L3recS1320);
      moonbit_decref(_M0L2idS3302);
      moonbit_decref(_M0L3recS1320);
      _M0L6_2atmpS3303 = _M0L2__S1319 + 1;
      _M0L2__S1319 = _M0L6_2atmpS3303;
      continue;
    }
    break;
  }
  return _M0L1dS1315;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__14(
  
) {
  void* _M0L11_2atry__errS1314;
  int32_t _M0L6raisedS1312;
  struct moonbit_result_3 _tmp_4276;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3299;
  void* _M0L4NoneS3298;
  struct moonbit_result_0 _result_4277;
  #line 163 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 165 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4276
  = _M0FP26biolab8bio__seq11seqio__read((moonbit_string_t)moonbit_string_literal_0.data, (moonbit_string_t)moonbit_string_literal_57.data);
  if (_tmp_4276.tag) {
    struct _M0TP26biolab8bio__seq9SeqRecord* const _M0L5_2aokS3300 =
      _tmp_4276.data.ok;
    _M0L6_2atmpS3299 = _M0L5_2aokS3300;
  } else {
    void* const _M0L6_2aerrS3301 = _tmp_4276.data.err;
    _M0L11_2atry__errS1314 = _M0L6_2aerrS3301;
    goto join_1313;
  }
  moonbit_decref(_M0L6_2atmpS3299);
  _M0L6raisedS1312 = 0;
  goto joinlet_4275;
  join_1313:;
  moonbit_decref(_M0L11_2atry__errS1314);
  _M0L6raisedS1312 = 1;
  joinlet_4275:;
  _M0L4NoneS3298
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 170 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4277
  = _M0FPB12assert__true(_M0L6raisedS1312, _M0L4NoneS3298, (moonbit_string_t)moonbit_string_literal_58.data);
  moonbit_decref(_M0L4NoneS3298);
  return _result_4277;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__13(
  
) {
  void* _M0L11_2atry__errS1311;
  int32_t _M0L6raisedS1309;
  struct moonbit_result_3 _tmp_4279;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3295;
  void* _M0L4NoneS3294;
  struct moonbit_result_0 _result_4280;
  #line 152 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 154 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4279
  = _M0FP26biolab8bio__seq11seqio__read((moonbit_string_t)moonbit_string_literal_59.data, (moonbit_string_t)moonbit_string_literal_57.data);
  if (_tmp_4279.tag) {
    struct _M0TP26biolab8bio__seq9SeqRecord* const _M0L5_2aokS3296 =
      _tmp_4279.data.ok;
    _M0L6_2atmpS3295 = _M0L5_2aokS3296;
  } else {
    void* const _M0L6_2aerrS3297 = _tmp_4279.data.err;
    _M0L11_2atry__errS1311 = _M0L6_2aerrS3297;
    goto join_1310;
  }
  moonbit_decref(_M0L6_2atmpS3295);
  _M0L6raisedS1309 = 0;
  goto joinlet_4278;
  join_1310:;
  moonbit_decref(_M0L11_2atry__errS1311);
  _M0L6raisedS1309 = 1;
  joinlet_4278:;
  _M0L4NoneS3294
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 159 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4280
  = _M0FPB12assert__true(_M0L6raisedS1309, _M0L4NoneS3294, (moonbit_string_t)moonbit_string_literal_60.data);
  moonbit_decref(_M0L4NoneS3294);
  return _result_4280;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__12(
  
) {
  struct moonbit_result_3 _tmp_4281;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1308;
  moonbit_string_t _M0L2idS3281;
  moonbit_string_t _M0L6_2atmpS3282;
  struct moonbit_result_0 _tmp_4283;
  moonbit_string_t _M0L11descriptionS3285;
  moonbit_string_t _M0L6_2atmpS3286;
  struct moonbit_result_0 _tmp_4285;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3630;
  int32_t _M0L6_2acntS3985;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3291;
  moonbit_string_t _M0L6_2atmpS3289;
  moonbit_string_t _M0L6_2atmpS3290;
  struct moonbit_result_0 _result_4287;
  #line 144 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 145 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4281
  = _M0FP26biolab8bio__seq11seqio__read((moonbit_string_t)moonbit_string_literal_61.data, (moonbit_string_t)moonbit_string_literal_57.data);
  if (_tmp_4281.tag) {
    struct _M0TP26biolab8bio__seq9SeqRecord* const _M0L5_2aokS3292 =
      _tmp_4281.data.ok;
    _M0L3recS1308 = _M0L5_2aokS3292;
  } else {
    void* const _M0L6_2aerrS3293 = _tmp_4281.data.err;
    struct moonbit_result_0 _result_4282;
    _result_4282.tag = 0;
    _result_4282.data.err = _M0L6_2aerrS3293;
    return _result_4282;
  }
  _M0L2idS3281 = _M0L3recS1308->$1;
  _M0L6_2atmpS3282 = 0;
  moonbit_incref(_M0L2idS3281);
  #line 146 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4283
  = _M0FPB10assert__eqGsE(_M0L2idS3281, (moonbit_string_t)moonbit_string_literal_62.data, _M0L6_2atmpS3282, (moonbit_string_t)moonbit_string_literal_63.data);
  moonbit_decref(_M0L2idS3281);
  if (_M0L6_2atmpS3282) {
    moonbit_decref(_M0L6_2atmpS3282);
  }
  if (_tmp_4283.tag) {
    int32_t const _M0L5_2aokS3283 = _tmp_4283.data.ok;
  } else {
    void* const _M0L6_2aerrS3284 = _tmp_4283.data.err;
    struct moonbit_result_0 _result_4284;
    moonbit_decref(_M0L3recS1308);
    _result_4284.tag = 0;
    _result_4284.data.err = _M0L6_2aerrS3284;
    return _result_4284;
  }
  _M0L11descriptionS3285 = _M0L3recS1308->$3;
  _M0L6_2atmpS3286 = 0;
  moonbit_incref(_M0L11descriptionS3285);
  #line 147 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4285
  = _M0FPB10assert__eqGsE(_M0L11descriptionS3285, (moonbit_string_t)moonbit_string_literal_64.data, _M0L6_2atmpS3286, (moonbit_string_t)moonbit_string_literal_65.data);
  moonbit_decref(_M0L11descriptionS3285);
  if (_M0L6_2atmpS3286) {
    moonbit_decref(_M0L6_2atmpS3286);
  }
  if (_tmp_4285.tag) {
    int32_t const _M0L5_2aokS3287 = _tmp_4285.data.ok;
  } else {
    void* const _M0L6_2aerrS3288 = _tmp_4285.data.err;
    struct moonbit_result_0 _result_4286;
    moonbit_decref(_M0L3recS1308);
    _result_4286.tag = 0;
    _result_4286.data.err = _M0L6_2aerrS3288;
    return _result_4286;
  }
  _M0L8_2afieldS3630 = _M0L3recS1308->$0;
  _M0L6_2acntS3985 = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1308));
  if (_M0L6_2acntS3985 > 1) {
    int32_t _M0L11_2anew__cntS3990 = _M0L6_2acntS3985 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1308), _M0L11_2anew__cntS3990);
    moonbit_incref(_M0L8_2afieldS3630);
  } else if (_M0L6_2acntS3985 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3989 = _M0L3recS1308->$4;
    moonbit_string_t _M0L8_2afieldS3988;
    moonbit_string_t _M0L8_2afieldS3987;
    moonbit_string_t _M0L8_2afieldS3986;
    moonbit_decref(_M0L8_2afieldS3989);
    _M0L8_2afieldS3988 = _M0L3recS1308->$3;
    moonbit_decref(_M0L8_2afieldS3988);
    _M0L8_2afieldS3987 = _M0L3recS1308->$2;
    moonbit_decref(_M0L8_2afieldS3987);
    _M0L8_2afieldS3986 = _M0L3recS1308->$1;
    moonbit_decref(_M0L8_2afieldS3986);
    #line 148 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L3recS1308);
  }
  _M0L3seqS3291 = _M0L8_2afieldS3630;
  #line 148 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3289 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3291);
  moonbit_decref(_M0L3seqS3291);
  _M0L6_2atmpS3290 = 0;
  #line 148 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4287
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3289, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3290, (moonbit_string_t)moonbit_string_literal_66.data);
  moonbit_decref(_M0L6_2atmpS3289);
  if (_M0L6_2atmpS3290) {
    moonbit_decref(_M0L6_2atmpS3290);
  }
  return _result_4287;
}

struct moonbit_result_3 _M0FP26biolab8bio__seq11seqio__read(
  moonbit_string_t _M0L7contentS1306,
  moonbit_string_t _M0L6formatS1307
) {
  struct moonbit_result_1 _tmp_4288;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1305;
  int32_t _M0L6_2atmpS3274;
  int32_t _M0L6_2atmpS3276;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3278;
  struct moonbit_result_3 _result_4292;
  #line 32 "/home/developer/Documents/1/seqio.mbt"
  #line 36 "/home/developer/Documents/1/seqio.mbt"
  _tmp_4288
  = _M0FP26biolab8bio__seq12seqio__parse(_M0L7contentS1306, _M0L6formatS1307);
  if (_tmp_4288.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS3279 =
      _tmp_4288.data.ok;
    _M0L7recordsS1305 = _M0L5_2aokS3279;
  } else {
    void* const _M0L6_2aerrS3280 = _tmp_4288.data.err;
    struct moonbit_result_3 _result_4289;
    _result_4289.tag = 0;
    _result_4289.data.err = _M0L6_2aerrS3280;
    return _result_4289;
  }
  #line 37 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS3274
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1305);
  if (_M0L6_2atmpS3274 == 0) {
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3275;
    struct moonbit_result_3 _result_4290;
    moonbit_decref(_M0L7recordsS1305);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3275
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3275)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3275)->$0
    = (moonbit_string_t)moonbit_string_literal_67.data;
    _result_4290.tag = 0;
    _result_4290.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3275;
    return _result_4290;
  }
  #line 40 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS3276
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1305);
  if (_M0L6_2atmpS3276 > 1) {
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3277;
    struct moonbit_result_3 _result_4291;
    moonbit_decref(_M0L7recordsS1305);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3277
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3277)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3277)->$0
    = (moonbit_string_t)moonbit_string_literal_68.data;
    _result_4291.tag = 0;
    _result_4291.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3277;
    return _result_4291;
  }
  #line 43 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS3278
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1305, 0);
  moonbit_decref(_M0L7recordsS1305);
  _result_4292.tag = 1;
  _result_4292.data.ok = _M0L6_2atmpS3278;
  return _result_4292;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__11(
  
) {
  void* _M0L11_2atry__errS1304;
  int32_t _M0L6raisedS1302;
  struct moonbit_result_1 _tmp_4294;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS3271;
  void* _M0L4NoneS3270;
  struct moonbit_result_0 _result_4295;
  #line 133 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 135 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4294
  = _M0FP26biolab8bio__seq12seqio__parse((moonbit_string_t)moonbit_string_literal_69.data, (moonbit_string_t)moonbit_string_literal_70.data);
  if (_tmp_4294.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS3272 =
      _tmp_4294.data.ok;
    _M0L6_2atmpS3271 = _M0L5_2aokS3272;
  } else {
    void* const _M0L6_2aerrS3273 = _tmp_4294.data.err;
    _M0L11_2atry__errS1304 = _M0L6_2aerrS3273;
    goto join_1303;
  }
  moonbit_decref(_M0L6_2atmpS3271);
  _M0L6raisedS1302 = 0;
  goto joinlet_4293;
  join_1303:;
  moonbit_decref(_M0L11_2atry__errS1304);
  _M0L6raisedS1302 = 1;
  joinlet_4293:;
  _M0L4NoneS3270
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 140 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4295
  = _M0FPB12assert__true(_M0L6raisedS1302, _M0L4NoneS3270, (moonbit_string_t)moonbit_string_literal_71.data);
  moonbit_decref(_M0L4NoneS3270);
  return _result_4295;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq46____test__736571696f5f7762746573742e6d6274__10(
  
) {
  struct moonbit_result_1 _tmp_4296;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1295;
  int32_t _M0L6_2atmpS3246;
  moonbit_string_t _M0L6_2atmpS3247;
  struct moonbit_result_0 _tmp_4298;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3252;
  moonbit_string_t _M0L8_2afieldS3634;
  int32_t _M0L6_2acntS3991;
  moonbit_string_t _M0L2idS3250;
  moonbit_string_t _M0L6_2atmpS3251;
  struct moonbit_result_0 _tmp_4300;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3258;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3633;
  int32_t _M0L6_2acntS3997;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3257;
  moonbit_string_t _M0L6_2atmpS3255;
  moonbit_string_t _M0L6_2atmpS3256;
  struct moonbit_result_0 _tmp_4302;
  struct _M0TPB5ArrayGiE* _M0L1qS1297;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3267;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1298;
  int32_t* _M0L6_2atmpS3264;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS3263;
  int32_t _M0L6_2atmpS3261;
  void* _M0L4NoneS3262;
  struct moonbit_result_0 _result_4306;
  #line 121 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 122 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4296
  = _M0FP26biolab8bio__seq12seqio__parse((moonbit_string_t)moonbit_string_literal_9.data, (moonbit_string_t)moonbit_string_literal_72.data);
  if (_tmp_4296.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS3268 =
      _tmp_4296.data.ok;
    _M0L7recordsS1295 = _M0L5_2aokS3268;
  } else {
    void* const _M0L6_2aerrS3269 = _tmp_4296.data.err;
    struct moonbit_result_0 _result_4297;
    _result_4297.tag = 0;
    _result_4297.data.err = _M0L6_2aerrS3269;
    return _result_4297;
  }
  #line 123 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3246
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1295);
  _M0L6_2atmpS3247 = 0;
  #line 123 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4298
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3246, 1, _M0L6_2atmpS3247, (moonbit_string_t)moonbit_string_literal_73.data);
  if (_M0L6_2atmpS3247) {
    moonbit_decref(_M0L6_2atmpS3247);
  }
  if (_tmp_4298.tag) {
    int32_t const _M0L5_2aokS3248 = _tmp_4298.data.ok;
  } else {
    void* const _M0L6_2aerrS3249 = _tmp_4298.data.err;
    struct moonbit_result_0 _result_4299;
    moonbit_decref(_M0L7recordsS1295);
    _result_4299.tag = 0;
    _result_4299.data.err = _M0L6_2aerrS3249;
    return _result_4299;
  }
  #line 124 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3252
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1295, 0);
  _M0L8_2afieldS3634 = _M0L6_2atmpS3252->$1;
  _M0L6_2acntS3991
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3252));
  if (_M0L6_2acntS3991 > 1) {
    int32_t _M0L11_2anew__cntS3996 = _M0L6_2acntS3991 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3252), _M0L11_2anew__cntS3996);
    moonbit_incref(_M0L8_2afieldS3634);
  } else if (_M0L6_2acntS3991 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3995 =
      _M0L6_2atmpS3252->$4;
    moonbit_string_t _M0L8_2afieldS3994;
    moonbit_string_t _M0L8_2afieldS3993;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3992;
    moonbit_decref(_M0L8_2afieldS3995);
    _M0L8_2afieldS3994 = _M0L6_2atmpS3252->$3;
    moonbit_decref(_M0L8_2afieldS3994);
    _M0L8_2afieldS3993 = _M0L6_2atmpS3252->$2;
    moonbit_decref(_M0L8_2afieldS3993);
    _M0L8_2afieldS3992 = _M0L6_2atmpS3252->$0;
    moonbit_decref(_M0L8_2afieldS3992);
    #line 124 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3252);
  }
  _M0L2idS3250 = _M0L8_2afieldS3634;
  _M0L6_2atmpS3251 = 0;
  #line 124 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4300
  = _M0FPB10assert__eqGsE(_M0L2idS3250, (moonbit_string_t)moonbit_string_literal_11.data, _M0L6_2atmpS3251, (moonbit_string_t)moonbit_string_literal_74.data);
  moonbit_decref(_M0L2idS3250);
  if (_M0L6_2atmpS3251) {
    moonbit_decref(_M0L6_2atmpS3251);
  }
  if (_tmp_4300.tag) {
    int32_t const _M0L5_2aokS3253 = _tmp_4300.data.ok;
  } else {
    void* const _M0L6_2aerrS3254 = _tmp_4300.data.err;
    struct moonbit_result_0 _result_4301;
    moonbit_decref(_M0L7recordsS1295);
    _result_4301.tag = 0;
    _result_4301.data.err = _M0L6_2aerrS3254;
    return _result_4301;
  }
  #line 125 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3258
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1295, 0);
  _M0L8_2afieldS3633 = _M0L6_2atmpS3258->$0;
  _M0L6_2acntS3997
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3258));
  if (_M0L6_2acntS3997 > 1) {
    int32_t _M0L11_2anew__cntS4002 = _M0L6_2acntS3997 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3258), _M0L11_2anew__cntS4002);
    moonbit_incref(_M0L8_2afieldS3633);
  } else if (_M0L6_2acntS3997 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4001 =
      _M0L6_2atmpS3258->$4;
    moonbit_string_t _M0L8_2afieldS4000;
    moonbit_string_t _M0L8_2afieldS3999;
    moonbit_string_t _M0L8_2afieldS3998;
    moonbit_decref(_M0L8_2afieldS4001);
    _M0L8_2afieldS4000 = _M0L6_2atmpS3258->$3;
    moonbit_decref(_M0L8_2afieldS4000);
    _M0L8_2afieldS3999 = _M0L6_2atmpS3258->$2;
    moonbit_decref(_M0L8_2afieldS3999);
    _M0L8_2afieldS3998 = _M0L6_2atmpS3258->$1;
    moonbit_decref(_M0L8_2afieldS3998);
    #line 125 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3258);
  }
  _M0L3seqS3257 = _M0L8_2afieldS3633;
  #line 125 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3255 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3257);
  moonbit_decref(_M0L3seqS3257);
  _M0L6_2atmpS3256 = 0;
  #line 125 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4302
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3255, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3256, (moonbit_string_t)moonbit_string_literal_75.data);
  moonbit_decref(_M0L6_2atmpS3255);
  if (_M0L6_2atmpS3256) {
    moonbit_decref(_M0L6_2atmpS3256);
  }
  if (_tmp_4302.tag) {
    int32_t const _M0L5_2aokS3259 = _tmp_4302.data.ok;
  } else {
    void* const _M0L6_2aerrS3260 = _tmp_4302.data.err;
    struct moonbit_result_0 _result_4303;
    moonbit_decref(_M0L7recordsS1295);
    _result_4303.tag = 0;
    _result_4303.data.err = _M0L6_2aerrS3260;
    return _result_4303;
  }
  #line 126 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3267
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1295, 0);
  moonbit_decref(_M0L7recordsS1295);
  #line 126 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1298
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L6_2atmpS3267);
  moonbit_decref(_M0L6_2atmpS3267);
  if (_M0L7_2abindS1298 == 0) {
    moonbit_string_t _M0L7_2abindS1301;
    int32_t _M0L6_2atmpS3266;
    struct _M0TPC16string10StringView _M0L6_2atmpS3265;
    struct moonbit_result_0 _result_4305;
    if (_M0L7_2abindS1298) {
      moonbit_decref(_M0L7_2abindS1298);
    }
    _M0L7_2abindS1301 = (moonbit_string_t)moonbit_string_literal_76.data;
    _M0L6_2atmpS3266 = Moonbit_array_length(_M0L7_2abindS1301);
    _M0L6_2atmpS3265
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1301, .$1 = 0, .$2 = _M0L6_2atmpS3266
    };
    #line 128 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _result_4305
    = _M0FPB4failGuE(_M0L6_2atmpS3265, (moonbit_string_t)moonbit_string_literal_77.data);
    moonbit_decref(_M0L6_2atmpS3265.$0);
    return _result_4305;
  } else {
    struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1299 = _M0L7_2abindS1298;
    struct _M0TPB5ArrayGiE* _M0L4_2aqS1300 = _M0L7_2aSomeS1299;
    _M0L1qS1297 = _M0L4_2aqS1300;
    goto join_1296;
  }
  join_1296:;
  _M0L6_2atmpS3264 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS3264[0] = 40;
  _M0L6_2atmpS3264[1] = 40;
  _M0L6_2atmpS3264[2] = 40;
  _M0L6_2atmpS3264[3] = 40;
  _M0L6_2atmpS3263
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS3263)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS3263->$0 = _M0L6_2atmpS3264;
  _M0L6_2atmpS3263->$1 = 4;
  #line 127 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3261
  = _M0IPC15array5ArrayPB2Eq5equalGiE(_M0L1qS1297, _M0L6_2atmpS3263);
  moonbit_decref(_M0L1qS1297);
  moonbit_decref(_M0L6_2atmpS3263);
  _M0L4NoneS3262
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 127 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4306
  = _M0FPB12assert__true(_M0L6_2atmpS3261, _M0L4NoneS3262, (moonbit_string_t)moonbit_string_literal_78.data);
  moonbit_decref(_M0L4NoneS3262);
  return _result_4306;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__9(
  
) {
  struct moonbit_result_1 _tmp_4307;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1294;
  int32_t _M0L6_2atmpS3226;
  moonbit_string_t _M0L6_2atmpS3227;
  struct moonbit_result_0 _tmp_4309;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3232;
  moonbit_string_t _M0L8_2afieldS3637;
  int32_t _M0L6_2acntS4003;
  moonbit_string_t _M0L2idS3230;
  moonbit_string_t _M0L6_2atmpS3231;
  struct moonbit_result_0 _tmp_4311;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3237;
  moonbit_string_t _M0L8_2afieldS3636;
  int32_t _M0L6_2acntS4009;
  moonbit_string_t _M0L11descriptionS3235;
  moonbit_string_t _M0L6_2atmpS3236;
  struct moonbit_result_0 _tmp_4313;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3243;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3635;
  int32_t _M0L6_2acntS4015;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS3242;
  moonbit_string_t _M0L6_2atmpS3240;
  moonbit_string_t _M0L6_2atmpS3241;
  struct moonbit_result_0 _result_4315;
  #line 112 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 113 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4307
  = _M0FP26biolab8bio__seq12seqio__parse((moonbit_string_t)moonbit_string_literal_79.data, (moonbit_string_t)moonbit_string_literal_57.data);
  if (_tmp_4307.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS3244 =
      _tmp_4307.data.ok;
    _M0L7recordsS1294 = _M0L5_2aokS3244;
  } else {
    void* const _M0L6_2aerrS3245 = _tmp_4307.data.err;
    struct moonbit_result_0 _result_4308;
    _result_4308.tag = 0;
    _result_4308.data.err = _M0L6_2aerrS3245;
    return _result_4308;
  }
  #line 114 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3226
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1294);
  _M0L6_2atmpS3227 = 0;
  #line 114 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4309
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3226, 1, _M0L6_2atmpS3227, (moonbit_string_t)moonbit_string_literal_80.data);
  if (_M0L6_2atmpS3227) {
    moonbit_decref(_M0L6_2atmpS3227);
  }
  if (_tmp_4309.tag) {
    int32_t const _M0L5_2aokS3228 = _tmp_4309.data.ok;
  } else {
    void* const _M0L6_2aerrS3229 = _tmp_4309.data.err;
    struct moonbit_result_0 _result_4310;
    moonbit_decref(_M0L7recordsS1294);
    _result_4310.tag = 0;
    _result_4310.data.err = _M0L6_2aerrS3229;
    return _result_4310;
  }
  #line 115 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3232
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1294, 0);
  _M0L8_2afieldS3637 = _M0L6_2atmpS3232->$1;
  _M0L6_2acntS4003
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3232));
  if (_M0L6_2acntS4003 > 1) {
    int32_t _M0L11_2anew__cntS4008 = _M0L6_2acntS4003 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3232), _M0L11_2anew__cntS4008);
    moonbit_incref(_M0L8_2afieldS3637);
  } else if (_M0L6_2acntS4003 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4007 =
      _M0L6_2atmpS3232->$4;
    moonbit_string_t _M0L8_2afieldS4006;
    moonbit_string_t _M0L8_2afieldS4005;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4004;
    moonbit_decref(_M0L8_2afieldS4007);
    _M0L8_2afieldS4006 = _M0L6_2atmpS3232->$3;
    moonbit_decref(_M0L8_2afieldS4006);
    _M0L8_2afieldS4005 = _M0L6_2atmpS3232->$2;
    moonbit_decref(_M0L8_2afieldS4005);
    _M0L8_2afieldS4004 = _M0L6_2atmpS3232->$0;
    moonbit_decref(_M0L8_2afieldS4004);
    #line 115 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3232);
  }
  _M0L2idS3230 = _M0L8_2afieldS3637;
  _M0L6_2atmpS3231 = 0;
  #line 115 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4311
  = _M0FPB10assert__eqGsE(_M0L2idS3230, (moonbit_string_t)moonbit_string_literal_11.data, _M0L6_2atmpS3231, (moonbit_string_t)moonbit_string_literal_81.data);
  moonbit_decref(_M0L2idS3230);
  if (_M0L6_2atmpS3231) {
    moonbit_decref(_M0L6_2atmpS3231);
  }
  if (_tmp_4311.tag) {
    int32_t const _M0L5_2aokS3233 = _tmp_4311.data.ok;
  } else {
    void* const _M0L6_2aerrS3234 = _tmp_4311.data.err;
    struct moonbit_result_0 _result_4312;
    moonbit_decref(_M0L7recordsS1294);
    _result_4312.tag = 0;
    _result_4312.data.err = _M0L6_2aerrS3234;
    return _result_4312;
  }
  #line 116 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3237
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1294, 0);
  _M0L8_2afieldS3636 = _M0L6_2atmpS3237->$3;
  _M0L6_2acntS4009
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3237));
  if (_M0L6_2acntS4009 > 1) {
    int32_t _M0L11_2anew__cntS4014 = _M0L6_2acntS4009 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3237), _M0L11_2anew__cntS4014);
    moonbit_incref(_M0L8_2afieldS3636);
  } else if (_M0L6_2acntS4009 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4013 =
      _M0L6_2atmpS3237->$4;
    moonbit_string_t _M0L8_2afieldS4012;
    moonbit_string_t _M0L8_2afieldS4011;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4010;
    moonbit_decref(_M0L8_2afieldS4013);
    _M0L8_2afieldS4012 = _M0L6_2atmpS3237->$2;
    moonbit_decref(_M0L8_2afieldS4012);
    _M0L8_2afieldS4011 = _M0L6_2atmpS3237->$1;
    moonbit_decref(_M0L8_2afieldS4011);
    _M0L8_2afieldS4010 = _M0L6_2atmpS3237->$0;
    moonbit_decref(_M0L8_2afieldS4010);
    #line 116 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3237);
  }
  _M0L11descriptionS3235 = _M0L8_2afieldS3636;
  _M0L6_2atmpS3236 = 0;
  #line 116 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4313
  = _M0FPB10assert__eqGsE(_M0L11descriptionS3235, (moonbit_string_t)moonbit_string_literal_82.data, _M0L6_2atmpS3236, (moonbit_string_t)moonbit_string_literal_83.data);
  moonbit_decref(_M0L11descriptionS3235);
  if (_M0L6_2atmpS3236) {
    moonbit_decref(_M0L6_2atmpS3236);
  }
  if (_tmp_4313.tag) {
    int32_t const _M0L5_2aokS3238 = _tmp_4313.data.ok;
  } else {
    void* const _M0L6_2aerrS3239 = _tmp_4313.data.err;
    struct moonbit_result_0 _result_4314;
    moonbit_decref(_M0L7recordsS1294);
    _result_4314.tag = 0;
    _result_4314.data.err = _M0L6_2aerrS3239;
    return _result_4314;
  }
  #line 117 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3243
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1294, 0);
  moonbit_decref(_M0L7recordsS1294);
  _M0L8_2afieldS3635 = _M0L6_2atmpS3243->$0;
  _M0L6_2acntS4015
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3243));
  if (_M0L6_2acntS4015 > 1) {
    int32_t _M0L11_2anew__cntS4020 = _M0L6_2acntS4015 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3243), _M0L11_2anew__cntS4020);
    moonbit_incref(_M0L8_2afieldS3635);
  } else if (_M0L6_2acntS4015 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4019 =
      _M0L6_2atmpS3243->$4;
    moonbit_string_t _M0L8_2afieldS4018;
    moonbit_string_t _M0L8_2afieldS4017;
    moonbit_string_t _M0L8_2afieldS4016;
    moonbit_decref(_M0L8_2afieldS4019);
    _M0L8_2afieldS4018 = _M0L6_2atmpS3243->$3;
    moonbit_decref(_M0L8_2afieldS4018);
    _M0L8_2afieldS4017 = _M0L6_2atmpS3243->$2;
    moonbit_decref(_M0L8_2afieldS4017);
    _M0L8_2afieldS4016 = _M0L6_2atmpS3243->$1;
    moonbit_decref(_M0L8_2afieldS4016);
    #line 117 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS3243);
  }
  _M0L3seqS3242 = _M0L8_2afieldS3635;
  #line 117 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS3240 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS3242);
  moonbit_decref(_M0L3seqS3242);
  _M0L6_2atmpS3241 = 0;
  #line 117 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4315
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS3240, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS3241, (moonbit_string_t)moonbit_string_literal_84.data);
  moonbit_decref(_M0L6_2atmpS3240);
  if (_M0L6_2atmpS3241) {
    moonbit_decref(_M0L6_2atmpS3241);
  }
  return _result_4315;
}

struct moonbit_result_1 _M0FP26biolab8bio__seq12seqio__parse(
  moonbit_string_t _M0L7contentS1291,
  moonbit_string_t _M0L6formatS1293
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS3222;
  struct moonbit_result_1 _result_4320;
  #line 15 "/home/developer/Documents/1/seqio.mbt"
  if (
    _M0L6formatS1293 == (moonbit_string_t)moonbit_string_literal_57.data
    || Moonbit_array_length(_M0L6formatS1293) == 5
       && 0
          == memcmp(_M0L6formatS1293, (moonbit_string_t)moonbit_string_literal_57.data, 10)
  ) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS3223;
    struct moonbit_result_1 _result_4319;
    #line 20 "/home/developer/Documents/1/seqio.mbt"
    _M0L6_2atmpS3223
    = _M0FP26biolab8bio__seq12parse__fasta(_M0L7contentS1291);
    _result_4319.tag = 1;
    _result_4319.data.ok = _M0L6_2atmpS3223;
    return _result_4319;
  } else if (
           _M0L6formatS1293
           == (moonbit_string_t)moonbit_string_literal_72.data
           || Moonbit_array_length(_M0L6formatS1293) == 5
              && 0
                 == memcmp(_M0L6formatS1293, (moonbit_string_t)moonbit_string_literal_72.data, 10)
         ) {
    goto join_1292;
  } else if (
           _M0L6formatS1293
           == (moonbit_string_t)moonbit_string_literal_88.data
           || Moonbit_array_length(_M0L6formatS1293) == 12
              && 0
                 == memcmp(_M0L6formatS1293, (moonbit_string_t)moonbit_string_literal_88.data, 24)
         ) {
    goto join_1292;
  } else if (
           _M0L6formatS1293
           == (moonbit_string_t)moonbit_string_literal_87.data
           || Moonbit_array_length(_M0L6formatS1293) == 7
              && 0
                 == memcmp(_M0L6formatS1293, (moonbit_string_t)moonbit_string_literal_87.data, 14)
         ) {
    goto join_1290;
  } else if (
           _M0L6formatS1293
           == (moonbit_string_t)moonbit_string_literal_86.data
           || Moonbit_array_length(_M0L6formatS1293) == 2
              && 0
                 == memcmp(_M0L6formatS1293, (moonbit_string_t)moonbit_string_literal_86.data, 4)
         ) {
    goto join_1290;
  } else {
    moonbit_string_t _M0L6_2atmpS3225;
    void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3224;
    struct moonbit_result_1 _result_4318;
    #line 23 "/home/developer/Documents/1/seqio.mbt"
    _M0L6_2atmpS3225
    = moonbit_add_string((moonbit_string_t)moonbit_string_literal_85.data, _M0L6formatS1293);
    _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3224
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
    Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3224)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
    ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3224)->$0
    = _M0L6_2atmpS3225;
    _result_4318.tag = 0;
    _result_4318.data.err
    = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS3224;
    return _result_4318;
  }
  join_1292:;
  #line 21 "/home/developer/Documents/1/seqio.mbt"
  return _M0FP26biolab8bio__seq12parse__fastq(_M0L7contentS1291);
  join_1290:;
  #line 22 "/home/developer/Documents/1/seqio.mbt"
  _M0L6_2atmpS3222
  = _M0FP26biolab8bio__seq14parse__genbank(_M0L7contentS1291);
  _result_4320.tag = 1;
  _result_4320.data.ok = _M0L6_2atmpS3222;
  return _result_4320;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq14parse__genbank(
  moonbit_string_t _M0L7contentS1235
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1233;
  moonbit_string_t _M0L7_2abindS1236;
  int32_t _M0L6_2atmpS3221;
  struct _M0TPC16string10StringView _M0L6_2atmpS3220;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS3219;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS1234;
  int32_t _M0L1nS1237;
  struct _M0TPB8MutLocalGiE* _M0L1iS1238;
  int32_t _M0L6_2atmpS3218;
  int32_t _M0L2crS1239;
  int32_t _M0L6_2atmpS3217;
  int32_t _M0L2spS1240;
  int32_t _M0L6_2atmpS3216;
  int32_t _M0L3dotS1241;
  #line 41 "/home/developer/Documents/1/genbank_io.mbt"
  #line 42 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L7recordsS1233
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS1236 = (moonbit_string_t)moonbit_string_literal_89.data;
  _M0L6_2atmpS3221 = Moonbit_array_length(_M0L7_2abindS1236);
  _M0L6_2atmpS3220
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1236, .$1 = 0, .$2 = _M0L6_2atmpS3221
  };
  #line 43 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L6_2atmpS3219
  = _M0MPC16string6String5split(_M0L7contentS1235, _M0L6_2atmpS3220);
  moonbit_decref(_M0L6_2atmpS3220.$0);
  #line 43 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L5linesS1234
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS3219);
  moonbit_decref(_M0L6_2atmpS3219);
  #line 44 "/home/developer/Documents/1/genbank_io.mbt"
  _M0L1nS1237
  = _M0MPC15array5Array6lengthGRPC16string10StringViewE(_M0L5linesS1234);
  _M0L1iS1238
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1238)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1238->$0 = 0;
  _M0L6_2atmpS3218 = 13;
  _M0L2crS1239 = (uint16_t)_M0L6_2atmpS3218;
  _M0L6_2atmpS3217 = 32;
  _M0L2spS1240 = (uint16_t)_M0L6_2atmpS3217;
  _M0L6_2atmpS3216 = 46;
  _M0L3dotS1241 = (uint16_t)_M0L6_2atmpS3216;
  while (1) {
    int32_t _M0L3valS2975 = _M0L1iS1238->$0;
    if (_M0L3valS2975 < _M0L1nS1237) {
      int32_t _M0L3valS3215 = _M0L1iS1238->$0;
      struct _M0TPC16string10StringView _M0L4lineS1242;
      int32_t _M0L3lenS1243;
      struct _M0TPB8MutLocalGiE* _M0L7trimmedS1244;
      int32_t _if__result_4322;
      int32_t _M0L3valS2994;
      int32_t _if__result_4323;
      int32_t _M0L3valS3214;
      int32_t _M0L6_2atmpS3213;
      #line 51 "/home/developer/Documents/1/genbank_io.mbt"
      _M0L4lineS1242
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1234, _M0L3valS3215);
      #line 52 "/home/developer/Documents/1/genbank_io.mbt"
      _M0L3lenS1243 = _M0MPC16string10StringView6length(_M0L4lineS1242);
      _M0L7trimmedS1244
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L7trimmedS1244)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L7trimmedS1244->$0 = _M0L3lenS1243;
      if (_M0L3lenS1243 > 0) {
        int32_t _M0L6_2atmpS2977 = _M0L3lenS1243 - 1;
        int32_t _M0L6_2atmpS2976;
        #line 54 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2976
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1242, _M0L6_2atmpS2977);
        _if__result_4322 = _M0L6_2atmpS2976 == _M0L2crS1239;
      } else {
        _if__result_4322 = 0;
      }
      if (_if__result_4322) {
        int32_t _M0L6_2atmpS2978 = _M0L3lenS1243 - 1;
        _M0L7trimmedS1244->$0 = _M0L6_2atmpS2978;
      }
      _M0L3valS2994 = _M0L7trimmedS1244->$0;
      if (_M0L3valS2994 >= 5) {
        int32_t _M0L6_2atmpS2991;
        int32_t _M0L6_2atmpS2993;
        int32_t _M0L6_2atmpS2992;
        #line 58 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS2991
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1242, 0);
        _M0L6_2atmpS2993 = 76;
        _M0L6_2atmpS2992 = (uint16_t)_M0L6_2atmpS2993;
        if (_M0L6_2atmpS2991 == _M0L6_2atmpS2992) {
          int32_t _M0L6_2atmpS2988;
          int32_t _M0L6_2atmpS2990;
          int32_t _M0L6_2atmpS2989;
          #line 59 "/home/developer/Documents/1/genbank_io.mbt"
          _M0L6_2atmpS2988
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1242, 1);
          _M0L6_2atmpS2990 = 79;
          _M0L6_2atmpS2989 = (uint16_t)_M0L6_2atmpS2990;
          if (_M0L6_2atmpS2988 == _M0L6_2atmpS2989) {
            int32_t _M0L6_2atmpS2985;
            int32_t _M0L6_2atmpS2987;
            int32_t _M0L6_2atmpS2986;
            #line 60 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L6_2atmpS2985
            = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1242, 2);
            _M0L6_2atmpS2987 = 67;
            _M0L6_2atmpS2986 = (uint16_t)_M0L6_2atmpS2987;
            if (_M0L6_2atmpS2985 == _M0L6_2atmpS2986) {
              int32_t _M0L6_2atmpS2982;
              int32_t _M0L6_2atmpS2984;
              int32_t _M0L6_2atmpS2983;
              #line 61 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS2982
              = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1242, 3);
              _M0L6_2atmpS2984 = 85;
              _M0L6_2atmpS2983 = (uint16_t)_M0L6_2atmpS2984;
              if (_M0L6_2atmpS2982 == _M0L6_2atmpS2983) {
                int32_t _M0L6_2atmpS2979;
                int32_t _M0L6_2atmpS2981;
                int32_t _M0L6_2atmpS2980;
                #line 62 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS2979
                = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1242, 4);
                _M0L6_2atmpS2981 = 83;
                _M0L6_2atmpS2980 = (uint16_t)_M0L6_2atmpS2981;
                _if__result_4323 = _M0L6_2atmpS2979 == _M0L6_2atmpS2980;
              } else {
                _if__result_4323 = 0;
              }
            } else {
              _if__result_4323 = 0;
            }
          } else {
            _if__result_4323 = 0;
          }
        } else {
          _if__result_4323 = 0;
        }
      } else {
        _if__result_4323 = 0;
      }
      if (_if__result_4323) {
        struct _M0TPB8MutLocalGsE* _M0L4nameS1245 =
          (struct _M0TPB8MutLocalGsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGsE));
        struct _M0TPB8MutLocalGsE* _M0L2idS1246;
        struct _M0TPB8MutLocalGsE* _M0L11descriptionS1247;
        uint16_t* _M0L6_2atmpS3212;
        struct _M0TPB8MutLocalGAkE* _M0L8seq__bufS1248;
        struct _M0TPB8MutLocalGiE* _M0L8seq__lenS1249;
        struct _M0TPB8MutLocalGbE* _M0L10in__originS1250;
        int32_t _M0L3valS3211;
        int64_t _M0L6_2atmpS3210;
        struct _M0TPC16string10StringView _M0L11locus__restS1251;
        struct _M0TPB8MutLocalGiE* _M0L1jS1252;
        int32_t _M0L9rest__lenS1253;
        int32_t _M0L3valS3209;
        struct _M0TPB8MutLocalGiE* _M0L9name__endS1255;
        int32_t _M0L3valS3007;
        int32_t _M0L3valS3009;
        int64_t _M0L6_2atmpS3008;
        struct _M0TPC16string10StringView _M0L6_2atmpS3006;
        moonbit_string_t _M0L6_2atmpS3005;
        moonbit_string_t _M0L6_2aoldS3655;
        int32_t _M0L3valS3011;
        int32_t _M0L6_2atmpS3010;
        int32_t _M0L3valS3208;
        uint16_t* _M0L15final__seq__bufS1285;
        int32_t _M0L1kS1286;
        moonbit_string_t _M0L6_2atmpS3207;
        struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS3200;
        moonbit_string_t _M0L8_2afieldS3640;
        int32_t _M0L6_2acntS4021;
        moonbit_string_t _M0L3valS3201;
        moonbit_string_t _M0L8_2afieldS3639;
        int32_t _M0L6_2acntS4023;
        moonbit_string_t _M0L3valS3202;
        moonbit_string_t _M0L8_2afieldS3638;
        int32_t _M0L6_2acntS4025;
        moonbit_string_t _M0L3valS3203;
        struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1288;
        struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS3206;
        struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS3205;
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS3204;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3199;
        Moonbit_object_header(_M0L4nameS1245)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 40, 0);
        _M0L4nameS1245->$0 = (moonbit_string_t)moonbit_string_literal_90.data;
        _M0L2idS1246
        = (struct _M0TPB8MutLocalGsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGsE));
        Moonbit_object_header(_M0L2idS1246)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 40, 0);
        _M0L2idS1246->$0 = (moonbit_string_t)moonbit_string_literal_36.data;
        _M0L11descriptionS1247
        = (struct _M0TPB8MutLocalGsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGsE));
        Moonbit_object_header(_M0L11descriptionS1247)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 40, 0);
        _M0L11descriptionS1247->$0
        = (moonbit_string_t)moonbit_string_literal_38.data;
        _M0L6_2atmpS3212 = (uint16_t*)moonbit_make_string(0, 0);
        _M0L8seq__bufS1248
        = (struct _M0TPB8MutLocalGAkE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGAkE));
        Moonbit_object_header(_M0L8seq__bufS1248)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 43, 0);
        _M0L8seq__bufS1248->$0 = _M0L6_2atmpS3212;
        _M0L8seq__lenS1249
        = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
        Moonbit_object_header(_M0L8seq__lenS1249)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L8seq__lenS1249->$0 = 0;
        _M0L10in__originS1250
        = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
        Moonbit_object_header(_M0L10in__originS1250)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L10in__originS1250->$0 = 0;
        _M0L3valS3211 = _M0L7trimmedS1244->$0;
        moonbit_decref(_M0L7trimmedS1244);
        _M0L6_2atmpS3210 = (int64_t)_M0L3valS3211;
        #line 70 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L11locus__restS1251
        = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1242, 5, _M0L6_2atmpS3210);
        moonbit_decref(_M0L4lineS1242.$0);
        _M0L1jS1252
        = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
        Moonbit_object_header(_M0L1jS1252)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L1jS1252->$0 = 0;
        #line 72 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L9rest__lenS1253
        = _M0MPC16string10StringView6length(_M0L11locus__restS1251);
        while (1) {
          int32_t _M0L3valS2997 = _M0L1jS1252->$0;
          int32_t _if__result_4325;
          if (_M0L3valS2997 < _M0L9rest__lenS1253) {
            int32_t _M0L3valS2996 = _M0L1jS1252->$0;
            int32_t _M0L6_2atmpS2995;
            #line 73 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L6_2atmpS2995
            = _M0MPC16string10StringView11unsafe__get(_M0L11locus__restS1251, _M0L3valS2996);
            _if__result_4325 = _M0L6_2atmpS2995 == _M0L2spS1240;
          } else {
            _if__result_4325 = 0;
          }
          if (_if__result_4325) {
            int32_t _M0L3valS2999 = _M0L1jS1252->$0;
            int32_t _M0L6_2atmpS2998 = _M0L3valS2999 + 1;
            _M0L1jS1252->$0 = _M0L6_2atmpS2998;
            continue;
          }
          break;
        }
        _M0L3valS3209 = _M0L1jS1252->$0;
        _M0L9name__endS1255
        = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
        Moonbit_object_header(_M0L9name__endS1255)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
        _M0L9name__endS1255->$0 = _M0L3valS3209;
        while (1) {
          int32_t _M0L3valS3002 = _M0L9name__endS1255->$0;
          int32_t _if__result_4327;
          if (_M0L3valS3002 < _M0L9rest__lenS1253) {
            int32_t _M0L3valS3001 = _M0L9name__endS1255->$0;
            int32_t _M0L6_2atmpS3000;
            #line 77 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L6_2atmpS3000
            = _M0MPC16string10StringView11unsafe__get(_M0L11locus__restS1251, _M0L3valS3001);
            _if__result_4327 = _M0L6_2atmpS3000 != _M0L2spS1240;
          } else {
            _if__result_4327 = 0;
          }
          if (_if__result_4327) {
            int32_t _M0L3valS3004 = _M0L9name__endS1255->$0;
            int32_t _M0L6_2atmpS3003 = _M0L3valS3004 + 1;
            _M0L9name__endS1255->$0 = _M0L6_2atmpS3003;
            continue;
          }
          break;
        }
        _M0L3valS3007 = _M0L1jS1252->$0;
        moonbit_decref(_M0L1jS1252);
        _M0L3valS3009 = _M0L9name__endS1255->$0;
        moonbit_decref(_M0L9name__endS1255);
        _M0L6_2atmpS3008 = (int64_t)_M0L3valS3009;
        #line 80 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS3006
        = _M0MPC16string10StringView11sub_2einner(_M0L11locus__restS1251, _M0L3valS3007, _M0L6_2atmpS3008);
        moonbit_decref(_M0L11locus__restS1251.$0);
        #line 80 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS3005
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3006);
        moonbit_decref(_M0L6_2atmpS3006.$0);
        _M0L6_2aoldS3655 = _M0L4nameS1245->$0;
        moonbit_decref(_M0L6_2aoldS3655);
        _M0L4nameS1245->$0 = _M0L6_2atmpS3005;
        _M0L3valS3011 = _M0L1iS1238->$0;
        _M0L6_2atmpS3010 = _M0L3valS3011 + 1;
        _M0L1iS1238->$0 = _M0L6_2atmpS3010;
        while (1) {
          int32_t _M0L3valS3012 = _M0L1iS1238->$0;
          if (_M0L3valS3012 < _M0L1nS1237) {
            int32_t _M0L3valS3194 = _M0L1iS1238->$0;
            struct _M0TPC16string10StringView _M0L1lS1257;
            int32_t _M0L2llS1258;
            struct _M0TPB8MutLocalGiE* _M0L2ltS1259;
            int32_t _if__result_4329;
            int32_t _M0L3valS3022;
            int32_t _if__result_4330;
            int32_t _M0L3valS3053;
            int32_t _if__result_4331;
            int32_t _M0L3valS3123;
            int32_t _if__result_4340;
            int32_t _M0L3valS3160;
            int32_t _if__result_4345;
            int32_t _M0L3valS3193;
            int32_t _M0L6_2atmpS3192;
            #line 84 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L1lS1257
            = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1234, _M0L3valS3194);
            #line 85 "/home/developer/Documents/1/genbank_io.mbt"
            _M0L2llS1258 = _M0MPC16string10StringView6length(_M0L1lS1257);
            _M0L2ltS1259
            = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
            Moonbit_object_header(_M0L2ltS1259)->meta
            = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
            _M0L2ltS1259->$0 = _M0L2llS1258;
            if (_M0L2llS1258 > 0) {
              int32_t _M0L6_2atmpS3014 = _M0L2llS1258 - 1;
              int32_t _M0L6_2atmpS3013;
              #line 87 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3013
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, _M0L6_2atmpS3014);
              _if__result_4329 = _M0L6_2atmpS3013 == _M0L2crS1239;
            } else {
              _if__result_4329 = 0;
            }
            if (_if__result_4329) {
              int32_t _M0L6_2atmpS3015 = _M0L2llS1258 - 1;
              _M0L2ltS1259->$0 = _M0L6_2atmpS3015;
            }
            _M0L3valS3022 = _M0L2ltS1259->$0;
            if (_M0L3valS3022 >= 2) {
              int32_t _M0L6_2atmpS3019;
              int32_t _M0L6_2atmpS3021;
              int32_t _M0L6_2atmpS3020;
              #line 91 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3019
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 0);
              _M0L6_2atmpS3021 = 47;
              _M0L6_2atmpS3020 = (uint16_t)_M0L6_2atmpS3021;
              if (_M0L6_2atmpS3019 == _M0L6_2atmpS3020) {
                int32_t _M0L6_2atmpS3016;
                int32_t _M0L6_2atmpS3018;
                int32_t _M0L6_2atmpS3017;
                #line 92 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS3016
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 1);
                _M0L6_2atmpS3018 = 47;
                _M0L6_2atmpS3017 = (uint16_t)_M0L6_2atmpS3018;
                _if__result_4330 = _M0L6_2atmpS3016 == _M0L6_2atmpS3017;
              } else {
                _if__result_4330 = 0;
              }
            } else {
              _if__result_4330 = 0;
            }
            if (_if__result_4330) {
              moonbit_decref(_M0L2ltS1259);
              moonbit_decref(_M0L1lS1257.$0);
              moonbit_decref(_M0L10in__originS1250);
              break;
            }
            _M0L3valS3053 = _M0L2ltS1259->$0;
            if (_M0L3valS3053 >= 10) {
              int32_t _M0L6_2atmpS3050;
              int32_t _M0L6_2atmpS3052;
              int32_t _M0L6_2atmpS3051;
              #line 96 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3050
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 0);
              _M0L6_2atmpS3052 = 68;
              _M0L6_2atmpS3051 = (uint16_t)_M0L6_2atmpS3052;
              if (_M0L6_2atmpS3050 == _M0L6_2atmpS3051) {
                int32_t _M0L6_2atmpS3047;
                int32_t _M0L6_2atmpS3049;
                int32_t _M0L6_2atmpS3048;
                #line 97 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS3047
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 1);
                _M0L6_2atmpS3049 = 69;
                _M0L6_2atmpS3048 = (uint16_t)_M0L6_2atmpS3049;
                if (_M0L6_2atmpS3047 == _M0L6_2atmpS3048) {
                  int32_t _M0L6_2atmpS3044;
                  int32_t _M0L6_2atmpS3046;
                  int32_t _M0L6_2atmpS3045;
                  #line 98 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3044
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 2);
                  _M0L6_2atmpS3046 = 70;
                  _M0L6_2atmpS3045 = (uint16_t)_M0L6_2atmpS3046;
                  if (_M0L6_2atmpS3044 == _M0L6_2atmpS3045) {
                    int32_t _M0L6_2atmpS3041;
                    int32_t _M0L6_2atmpS3043;
                    int32_t _M0L6_2atmpS3042;
                    #line 99 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS3041
                    = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 3);
                    _M0L6_2atmpS3043 = 73;
                    _M0L6_2atmpS3042 = (uint16_t)_M0L6_2atmpS3043;
                    if (_M0L6_2atmpS3041 == _M0L6_2atmpS3042) {
                      int32_t _M0L6_2atmpS3038;
                      int32_t _M0L6_2atmpS3040;
                      int32_t _M0L6_2atmpS3039;
                      #line 100 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS3038
                      = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 4);
                      _M0L6_2atmpS3040 = 78;
                      _M0L6_2atmpS3039 = (uint16_t)_M0L6_2atmpS3040;
                      if (_M0L6_2atmpS3038 == _M0L6_2atmpS3039) {
                        int32_t _M0L6_2atmpS3035;
                        int32_t _M0L6_2atmpS3037;
                        int32_t _M0L6_2atmpS3036;
                        #line 101 "/home/developer/Documents/1/genbank_io.mbt"
                        _M0L6_2atmpS3035
                        = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 5);
                        _M0L6_2atmpS3037 = 73;
                        _M0L6_2atmpS3036 = (uint16_t)_M0L6_2atmpS3037;
                        if (_M0L6_2atmpS3035 == _M0L6_2atmpS3036) {
                          int32_t _M0L6_2atmpS3032;
                          int32_t _M0L6_2atmpS3034;
                          int32_t _M0L6_2atmpS3033;
                          #line 102 "/home/developer/Documents/1/genbank_io.mbt"
                          _M0L6_2atmpS3032
                          = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 6);
                          _M0L6_2atmpS3034 = 84;
                          _M0L6_2atmpS3033 = (uint16_t)_M0L6_2atmpS3034;
                          if (_M0L6_2atmpS3032 == _M0L6_2atmpS3033) {
                            int32_t _M0L6_2atmpS3029;
                            int32_t _M0L6_2atmpS3031;
                            int32_t _M0L6_2atmpS3030;
                            #line 103 "/home/developer/Documents/1/genbank_io.mbt"
                            _M0L6_2atmpS3029
                            = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 7);
                            _M0L6_2atmpS3031 = 73;
                            _M0L6_2atmpS3030 = (uint16_t)_M0L6_2atmpS3031;
                            if (_M0L6_2atmpS3029 == _M0L6_2atmpS3030) {
                              int32_t _M0L6_2atmpS3026;
                              int32_t _M0L6_2atmpS3028;
                              int32_t _M0L6_2atmpS3027;
                              #line 104 "/home/developer/Documents/1/genbank_io.mbt"
                              _M0L6_2atmpS3026
                              = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 8);
                              _M0L6_2atmpS3028 = 79;
                              _M0L6_2atmpS3027 = (uint16_t)_M0L6_2atmpS3028;
                              if (_M0L6_2atmpS3026 == _M0L6_2atmpS3027) {
                                int32_t _M0L6_2atmpS3023;
                                int32_t _M0L6_2atmpS3025;
                                int32_t _M0L6_2atmpS3024;
                                #line 105 "/home/developer/Documents/1/genbank_io.mbt"
                                _M0L6_2atmpS3023
                                = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 9);
                                _M0L6_2atmpS3025 = 78;
                                _M0L6_2atmpS3024 = (uint16_t)_M0L6_2atmpS3025;
                                _if__result_4331
                                = _M0L6_2atmpS3023 == _M0L6_2atmpS3024;
                              } else {
                                _if__result_4331 = 0;
                              }
                            } else {
                              _if__result_4331 = 0;
                            }
                          } else {
                            _if__result_4331 = 0;
                          }
                        } else {
                          _if__result_4331 = 0;
                        }
                      } else {
                        _if__result_4331 = 0;
                      }
                    } else {
                      _if__result_4331 = 0;
                    }
                  } else {
                    _if__result_4331 = 0;
                  }
                } else {
                  _if__result_4331 = 0;
                }
              } else {
                _if__result_4331 = 0;
              }
            } else {
              _if__result_4331 = 0;
            }
            if (_if__result_4331) {
              int32_t _M0L3valS3101 = _M0L2ltS1259->$0;
              int64_t _M0L6_2atmpS3100;
              struct _M0TPC16string10StringView _M0L9def__partS1261;
              struct _M0TPB8MutLocalGiE* _M0L2djS1262;
              int32_t _M0L3dplS1263;
              struct _M0TPB13StringBuilder* _M0L8def__bufS1265;
              int32_t _M0L3valS3061;
              int64_t _M0L6_2atmpS3062;
              struct _M0TPC16string10StringView _M0L6_2atmpS3060;
              moonbit_string_t _M0L6_2atmpS3059;
              int32_t _M0L3valS3064;
              int32_t _M0L6_2atmpS3063;
              moonbit_string_t _M0L6_2atmpS3085;
              moonbit_string_t _M0L6_2aoldS3654;
              moonbit_string_t _M0L3valS3092;
              int32_t _M0L6_2atmpS3091;
              int32_t _if__result_4339;
              moonbit_decref(_M0L2ltS1259);
              _M0L6_2atmpS3100 = (int64_t)_M0L3valS3101;
              #line 106 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L9def__partS1261
              = _M0MPC16string10StringView11sub_2einner(_M0L1lS1257, 10, _M0L6_2atmpS3100);
              moonbit_decref(_M0L1lS1257.$0);
              _M0L2djS1262
              = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
              Moonbit_object_header(_M0L2djS1262)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
              _M0L2djS1262->$0 = 0;
              #line 108 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L3dplS1263
              = _M0MPC16string10StringView6length(_M0L9def__partS1261);
              while (1) {
                int32_t _M0L3valS3056 = _M0L2djS1262->$0;
                int32_t _if__result_4333;
                if (_M0L3valS3056 < _M0L3dplS1263) {
                  int32_t _M0L3valS3055 = _M0L2djS1262->$0;
                  int32_t _M0L6_2atmpS3054;
                  #line 109 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3054
                  = _M0MPC16string10StringView11unsafe__get(_M0L9def__partS1261, _M0L3valS3055);
                  _if__result_4333 = _M0L6_2atmpS3054 == _M0L2spS1240;
                } else {
                  _if__result_4333 = 0;
                }
                if (_if__result_4333) {
                  int32_t _M0L3valS3058 = _M0L2djS1262->$0;
                  int32_t _M0L6_2atmpS3057 = _M0L3valS3058 + 1;
                  _M0L2djS1262->$0 = _M0L6_2atmpS3057;
                  continue;
                }
                break;
              }
              #line 112 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L8def__bufS1265
              = _M0MPB13StringBuilder21StringBuilder_2einner(0);
              _M0L3valS3061 = _M0L2djS1262->$0;
              moonbit_decref(_M0L2djS1262);
              _M0L6_2atmpS3062 = (int64_t)_M0L3dplS1263;
              #line 113 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3060
              = _M0MPC16string10StringView11sub_2einner(_M0L9def__partS1261, _M0L3valS3061, _M0L6_2atmpS3062);
              moonbit_decref(_M0L9def__partS1261.$0);
              #line 113 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3059
              = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3060);
              moonbit_decref(_M0L6_2atmpS3060.$0);
              #line 113 "/home/developer/Documents/1/genbank_io.mbt"
              _M0IPB13StringBuilderPB6Logger13write__string(_M0L8def__bufS1265, _M0L6_2atmpS3059);
              moonbit_decref(_M0L6_2atmpS3059);
              _M0L3valS3064 = _M0L1iS1238->$0;
              _M0L6_2atmpS3063 = _M0L3valS3064 + 1;
              _M0L1iS1238->$0 = _M0L6_2atmpS3063;
              while (1) {
                int32_t _M0L3valS3065 = _M0L1iS1238->$0;
                if (_M0L3valS3065 < _M0L1nS1237) {
                  int32_t _M0L3valS3084 = _M0L1iS1238->$0;
                  struct _M0TPC16string10StringView _M0L4contS1266;
                  int32_t _M0L2clS1267;
                  struct _M0TPB8MutLocalGiE* _M0L2ctS1268;
                  int32_t _if__result_4335;
                  int32_t _M0L3valS3070;
                  int32_t _if__result_4336;
                  struct _M0TPB8MutLocalGiE* _M0L2cjS1270;
                  int32_t _M0L3valS3079;
                  int32_t _M0L3valS3081;
                  int64_t _M0L6_2atmpS3080;
                  struct _M0TPC16string10StringView _M0L6_2atmpS3078;
                  moonbit_string_t _M0L6_2atmpS3077;
                  int32_t _M0L3valS3083;
                  int32_t _M0L6_2atmpS3082;
                  #line 116 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L4contS1266
                  = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1234, _M0L3valS3084);
                  #line 117 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L2clS1267
                  = _M0MPC16string10StringView6length(_M0L4contS1266);
                  _M0L2ctS1268
                  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
                  Moonbit_object_header(_M0L2ctS1268)->meta
                  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
                  _M0L2ctS1268->$0 = _M0L2clS1267;
                  if (_M0L2clS1267 > 0) {
                    int32_t _M0L6_2atmpS3067 = _M0L2clS1267 - 1;
                    int32_t _M0L6_2atmpS3066;
                    #line 119 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS3066
                    = _M0MPC16string10StringView11unsafe__get(_M0L4contS1266, _M0L6_2atmpS3067);
                    _if__result_4335 = _M0L6_2atmpS3066 == _M0L2crS1239;
                  } else {
                    _if__result_4335 = 0;
                  }
                  if (_if__result_4335) {
                    int32_t _M0L6_2atmpS3068 = _M0L2clS1267 - 1;
                    _M0L2ctS1268->$0 = _M0L6_2atmpS3068;
                  }
                  _M0L3valS3070 = _M0L2ctS1268->$0;
                  if (_M0L3valS3070 == 0) {
                    _if__result_4336 = 1;
                  } else {
                    int32_t _M0L6_2atmpS3069;
                    #line 122 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS3069
                    = _M0MPC16string10StringView11unsafe__get(_M0L4contS1266, 0);
                    _if__result_4336 = _M0L6_2atmpS3069 != _M0L2spS1240;
                  }
                  if (_if__result_4336) {
                    moonbit_decref(_M0L2ctS1268);
                    moonbit_decref(_M0L4contS1266.$0);
                    break;
                  }
                  #line 125 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0IPB13StringBuilderPB6Logger11write__char(_M0L8def__bufS1265, 32);
                  _M0L2cjS1270
                  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
                  Moonbit_object_header(_M0L2cjS1270)->meta
                  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
                  _M0L2cjS1270->$0 = 0;
                  while (1) {
                    int32_t _M0L3valS3073 = _M0L2cjS1270->$0;
                    int32_t _M0L3valS3074 = _M0L2ctS1268->$0;
                    int32_t _if__result_4338;
                    if (_M0L3valS3073 < _M0L3valS3074) {
                      int32_t _M0L3valS3072 = _M0L2cjS1270->$0;
                      int32_t _M0L6_2atmpS3071;
                      #line 127 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS3071
                      = _M0MPC16string10StringView11unsafe__get(_M0L4contS1266, _M0L3valS3072);
                      _if__result_4338 = _M0L6_2atmpS3071 == _M0L2spS1240;
                    } else {
                      _if__result_4338 = 0;
                    }
                    if (_if__result_4338) {
                      int32_t _M0L3valS3076 = _M0L2cjS1270->$0;
                      int32_t _M0L6_2atmpS3075 = _M0L3valS3076 + 1;
                      _M0L2cjS1270->$0 = _M0L6_2atmpS3075;
                      continue;
                    }
                    break;
                  }
                  _M0L3valS3079 = _M0L2cjS1270->$0;
                  moonbit_decref(_M0L2cjS1270);
                  _M0L3valS3081 = _M0L2ctS1268->$0;
                  moonbit_decref(_M0L2ctS1268);
                  _M0L6_2atmpS3080 = (int64_t)_M0L3valS3081;
                  #line 130 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3078
                  = _M0MPC16string10StringView11sub_2einner(_M0L4contS1266, _M0L3valS3079, _M0L6_2atmpS3080);
                  moonbit_decref(_M0L4contS1266.$0);
                  #line 130 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3077
                  = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3078);
                  moonbit_decref(_M0L6_2atmpS3078.$0);
                  #line 130 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0IPB13StringBuilderPB6Logger13write__string(_M0L8def__bufS1265, _M0L6_2atmpS3077);
                  moonbit_decref(_M0L6_2atmpS3077);
                  _M0L3valS3083 = _M0L1iS1238->$0;
                  _M0L6_2atmpS3082 = _M0L3valS3083 + 1;
                  _M0L1iS1238->$0 = _M0L6_2atmpS3082;
                  continue;
                }
                break;
              }
              #line 133 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3085
              = _M0MPB13StringBuilder10to__string(_M0L8def__bufS1265);
              moonbit_decref(_M0L8def__bufS1265);
              _M0L6_2aoldS3654 = _M0L11descriptionS1247->$0;
              moonbit_decref(_M0L6_2aoldS3654);
              _M0L11descriptionS1247->$0 = _M0L6_2atmpS3085;
              _M0L3valS3092 = _M0L11descriptionS1247->$0;
              _M0L6_2atmpS3091 = Moonbit_array_length(_M0L3valS3092);
              if (_M0L6_2atmpS3091 > 0) {
                moonbit_string_t _M0L3valS3087 = _M0L11descriptionS1247->$0;
                moonbit_string_t _M0L3valS3090 = _M0L11descriptionS1247->$0;
                int32_t _M0L6_2atmpS3089 =
                  Moonbit_array_length(_M0L3valS3090);
                int32_t _M0L6_2atmpS3088 = _M0L6_2atmpS3089 - 1;
                int32_t _M0L6_2atmpS3086 = _M0L3valS3087[_M0L6_2atmpS3088];
                _if__result_4339 = _M0L6_2atmpS3086 == _M0L3dotS1241;
              } else {
                _if__result_4339 = 0;
              }
              if (_if__result_4339) {
                moonbit_string_t _M0L3valS3095 = _M0L11descriptionS1247->$0;
                moonbit_string_t _M0L3valS3099 = _M0L11descriptionS1247->$0;
                int32_t _M0L6_2atmpS3098 =
                  Moonbit_array_length(_M0L3valS3099);
                int32_t _M0L6_2atmpS3097 = _M0L6_2atmpS3098 - 1;
                int64_t _M0L6_2atmpS3096 = (int64_t)_M0L6_2atmpS3097;
                struct _M0TPC16string10StringView _M0L6_2atmpS3094;
                moonbit_string_t _M0L6_2atmpS3093;
                moonbit_string_t _M0L6_2aoldS3648;
                moonbit_incref(_M0L3valS3095);
                #line 136 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS3094
                = _M0MPC16string6String11sub_2einner(_M0L3valS3095, 0, _M0L6_2atmpS3096);
                moonbit_decref(_M0L3valS3095);
                #line 136 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS3093
                = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3094);
                moonbit_decref(_M0L6_2atmpS3094.$0);
                _M0L6_2aoldS3648 = _M0L11descriptionS1247->$0;
                moonbit_decref(_M0L6_2aoldS3648);
                _M0L11descriptionS1247->$0 = _M0L6_2atmpS3093;
              }
              continue;
            }
            _M0L3valS3123 = _M0L2ltS1259->$0;
            if (_M0L3valS3123 >= 7) {
              int32_t _M0L6_2atmpS3120;
              int32_t _M0L6_2atmpS3122;
              int32_t _M0L6_2atmpS3121;
              #line 141 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3120
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 0);
              _M0L6_2atmpS3122 = 86;
              _M0L6_2atmpS3121 = (uint16_t)_M0L6_2atmpS3122;
              if (_M0L6_2atmpS3120 == _M0L6_2atmpS3121) {
                int32_t _M0L6_2atmpS3117;
                int32_t _M0L6_2atmpS3119;
                int32_t _M0L6_2atmpS3118;
                #line 142 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS3117
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 1);
                _M0L6_2atmpS3119 = 69;
                _M0L6_2atmpS3118 = (uint16_t)_M0L6_2atmpS3119;
                if (_M0L6_2atmpS3117 == _M0L6_2atmpS3118) {
                  int32_t _M0L6_2atmpS3114;
                  int32_t _M0L6_2atmpS3116;
                  int32_t _M0L6_2atmpS3115;
                  #line 143 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3114
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 2);
                  _M0L6_2atmpS3116 = 82;
                  _M0L6_2atmpS3115 = (uint16_t)_M0L6_2atmpS3116;
                  if (_M0L6_2atmpS3114 == _M0L6_2atmpS3115) {
                    int32_t _M0L6_2atmpS3111;
                    int32_t _M0L6_2atmpS3113;
                    int32_t _M0L6_2atmpS3112;
                    #line 144 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS3111
                    = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 3);
                    _M0L6_2atmpS3113 = 83;
                    _M0L6_2atmpS3112 = (uint16_t)_M0L6_2atmpS3113;
                    if (_M0L6_2atmpS3111 == _M0L6_2atmpS3112) {
                      int32_t _M0L6_2atmpS3108;
                      int32_t _M0L6_2atmpS3110;
                      int32_t _M0L6_2atmpS3109;
                      #line 145 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS3108
                      = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 4);
                      _M0L6_2atmpS3110 = 73;
                      _M0L6_2atmpS3109 = (uint16_t)_M0L6_2atmpS3110;
                      if (_M0L6_2atmpS3108 == _M0L6_2atmpS3109) {
                        int32_t _M0L6_2atmpS3105;
                        int32_t _M0L6_2atmpS3107;
                        int32_t _M0L6_2atmpS3106;
                        #line 146 "/home/developer/Documents/1/genbank_io.mbt"
                        _M0L6_2atmpS3105
                        = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 5);
                        _M0L6_2atmpS3107 = 79;
                        _M0L6_2atmpS3106 = (uint16_t)_M0L6_2atmpS3107;
                        if (_M0L6_2atmpS3105 == _M0L6_2atmpS3106) {
                          int32_t _M0L6_2atmpS3102;
                          int32_t _M0L6_2atmpS3104;
                          int32_t _M0L6_2atmpS3103;
                          #line 147 "/home/developer/Documents/1/genbank_io.mbt"
                          _M0L6_2atmpS3102
                          = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 6);
                          _M0L6_2atmpS3104 = 78;
                          _M0L6_2atmpS3103 = (uint16_t)_M0L6_2atmpS3104;
                          _if__result_4340
                          = _M0L6_2atmpS3102 == _M0L6_2atmpS3103;
                        } else {
                          _if__result_4340 = 0;
                        }
                      } else {
                        _if__result_4340 = 0;
                      }
                    } else {
                      _if__result_4340 = 0;
                    }
                  } else {
                    _if__result_4340 = 0;
                  }
                } else {
                  _if__result_4340 = 0;
                }
              } else {
                _if__result_4340 = 0;
              }
            } else {
              _if__result_4340 = 0;
            }
            if (_if__result_4340) {
              int32_t _M0L3valS3141 = _M0L2ltS1259->$0;
              int64_t _M0L6_2atmpS3140 = (int64_t)_M0L3valS3141;
              struct _M0TPC16string10StringView _M0L9ver__partS1272;
              struct _M0TPB8MutLocalGiE* _M0L2vjS1273;
              int32_t _M0L3vplS1274;
              int32_t _M0L3valS3139;
              struct _M0TPB8MutLocalGiE* _M0L8vid__endS1276;
              int32_t _M0L3valS3136;
              int32_t _M0L3valS3138;
              int64_t _M0L6_2atmpS3137;
              struct _M0TPC16string10StringView _M0L6_2atmpS3135;
              moonbit_string_t _M0L6_2atmpS3134;
              moonbit_string_t _M0L6_2aoldS3647;
              #line 148 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L9ver__partS1272
              = _M0MPC16string10StringView11sub_2einner(_M0L1lS1257, 7, _M0L6_2atmpS3140);
              _M0L2vjS1273
              = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
              Moonbit_object_header(_M0L2vjS1273)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
              _M0L2vjS1273->$0 = 0;
              #line 150 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L3vplS1274
              = _M0MPC16string10StringView6length(_M0L9ver__partS1272);
              while (1) {
                int32_t _M0L3valS3126 = _M0L2vjS1273->$0;
                int32_t _if__result_4342;
                if (_M0L3valS3126 < _M0L3vplS1274) {
                  int32_t _M0L3valS3125 = _M0L2vjS1273->$0;
                  int32_t _M0L6_2atmpS3124;
                  #line 151 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3124
                  = _M0MPC16string10StringView11unsafe__get(_M0L9ver__partS1272, _M0L3valS3125);
                  _if__result_4342 = _M0L6_2atmpS3124 == _M0L2spS1240;
                } else {
                  _if__result_4342 = 0;
                }
                if (_if__result_4342) {
                  int32_t _M0L3valS3128 = _M0L2vjS1273->$0;
                  int32_t _M0L6_2atmpS3127 = _M0L3valS3128 + 1;
                  _M0L2vjS1273->$0 = _M0L6_2atmpS3127;
                  continue;
                }
                break;
              }
              _M0L3valS3139 = _M0L2vjS1273->$0;
              _M0L8vid__endS1276
              = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
              Moonbit_object_header(_M0L8vid__endS1276)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
              _M0L8vid__endS1276->$0 = _M0L3valS3139;
              while (1) {
                int32_t _M0L3valS3131 = _M0L8vid__endS1276->$0;
                int32_t _if__result_4344;
                if (_M0L3valS3131 < _M0L3vplS1274) {
                  int32_t _M0L3valS3130 = _M0L8vid__endS1276->$0;
                  int32_t _M0L6_2atmpS3129;
                  #line 155 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3129
                  = _M0MPC16string10StringView11unsafe__get(_M0L9ver__partS1272, _M0L3valS3130);
                  _if__result_4344 = _M0L6_2atmpS3129 != _M0L2spS1240;
                } else {
                  _if__result_4344 = 0;
                }
                if (_if__result_4344) {
                  int32_t _M0L3valS3133 = _M0L8vid__endS1276->$0;
                  int32_t _M0L6_2atmpS3132 = _M0L3valS3133 + 1;
                  _M0L8vid__endS1276->$0 = _M0L6_2atmpS3132;
                  continue;
                }
                break;
              }
              _M0L3valS3136 = _M0L2vjS1273->$0;
              moonbit_decref(_M0L2vjS1273);
              _M0L3valS3138 = _M0L8vid__endS1276->$0;
              moonbit_decref(_M0L8vid__endS1276);
              _M0L6_2atmpS3137 = (int64_t)_M0L3valS3138;
              #line 158 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3135
              = _M0MPC16string10StringView11sub_2einner(_M0L9ver__partS1272, _M0L3valS3136, _M0L6_2atmpS3137);
              moonbit_decref(_M0L9ver__partS1272.$0);
              #line 158 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3134
              = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3135);
              moonbit_decref(_M0L6_2atmpS3135.$0);
              _M0L6_2aoldS3647 = _M0L2idS1246->$0;
              moonbit_decref(_M0L6_2aoldS3647);
              _M0L2idS1246->$0 = _M0L6_2atmpS3134;
            }
            _M0L3valS3160 = _M0L2ltS1259->$0;
            if (_M0L3valS3160 >= 6) {
              int32_t _M0L6_2atmpS3157;
              int32_t _M0L6_2atmpS3159;
              int32_t _M0L6_2atmpS3158;
              #line 161 "/home/developer/Documents/1/genbank_io.mbt"
              _M0L6_2atmpS3157
              = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 0);
              _M0L6_2atmpS3159 = 79;
              _M0L6_2atmpS3158 = (uint16_t)_M0L6_2atmpS3159;
              if (_M0L6_2atmpS3157 == _M0L6_2atmpS3158) {
                int32_t _M0L6_2atmpS3154;
                int32_t _M0L6_2atmpS3156;
                int32_t _M0L6_2atmpS3155;
                #line 162 "/home/developer/Documents/1/genbank_io.mbt"
                _M0L6_2atmpS3154
                = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 1);
                _M0L6_2atmpS3156 = 82;
                _M0L6_2atmpS3155 = (uint16_t)_M0L6_2atmpS3156;
                if (_M0L6_2atmpS3154 == _M0L6_2atmpS3155) {
                  int32_t _M0L6_2atmpS3151;
                  int32_t _M0L6_2atmpS3153;
                  int32_t _M0L6_2atmpS3152;
                  #line 163 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L6_2atmpS3151
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 2);
                  _M0L6_2atmpS3153 = 73;
                  _M0L6_2atmpS3152 = (uint16_t)_M0L6_2atmpS3153;
                  if (_M0L6_2atmpS3151 == _M0L6_2atmpS3152) {
                    int32_t _M0L6_2atmpS3148;
                    int32_t _M0L6_2atmpS3150;
                    int32_t _M0L6_2atmpS3149;
                    #line 164 "/home/developer/Documents/1/genbank_io.mbt"
                    _M0L6_2atmpS3148
                    = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 3);
                    _M0L6_2atmpS3150 = 71;
                    _M0L6_2atmpS3149 = (uint16_t)_M0L6_2atmpS3150;
                    if (_M0L6_2atmpS3148 == _M0L6_2atmpS3149) {
                      int32_t _M0L6_2atmpS3145;
                      int32_t _M0L6_2atmpS3147;
                      int32_t _M0L6_2atmpS3146;
                      #line 165 "/home/developer/Documents/1/genbank_io.mbt"
                      _M0L6_2atmpS3145
                      = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 4);
                      _M0L6_2atmpS3147 = 73;
                      _M0L6_2atmpS3146 = (uint16_t)_M0L6_2atmpS3147;
                      if (_M0L6_2atmpS3145 == _M0L6_2atmpS3146) {
                        int32_t _M0L6_2atmpS3142;
                        int32_t _M0L6_2atmpS3144;
                        int32_t _M0L6_2atmpS3143;
                        #line 166 "/home/developer/Documents/1/genbank_io.mbt"
                        _M0L6_2atmpS3142
                        = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, 5);
                        _M0L6_2atmpS3144 = 78;
                        _M0L6_2atmpS3143 = (uint16_t)_M0L6_2atmpS3144;
                        _if__result_4345
                        = _M0L6_2atmpS3142 == _M0L6_2atmpS3143;
                      } else {
                        _if__result_4345 = 0;
                      }
                    } else {
                      _if__result_4345 = 0;
                    }
                  } else {
                    _if__result_4345 = 0;
                  }
                } else {
                  _if__result_4345 = 0;
                }
              } else {
                _if__result_4345 = 0;
              }
            } else {
              _if__result_4345 = 0;
            }
            if (_if__result_4345) {
              int32_t _M0L3valS3162;
              int32_t _M0L6_2atmpS3161;
              moonbit_decref(_M0L2ltS1259);
              moonbit_decref(_M0L1lS1257.$0);
              _M0L10in__originS1250->$0 = 1;
              _M0L3valS3162 = _M0L1iS1238->$0;
              _M0L6_2atmpS3161 = _M0L3valS3162 + 1;
              _M0L1iS1238->$0 = _M0L6_2atmpS3161;
              continue;
            }
            if (_M0L10in__originS1250->$0) {
              int32_t _M0L3valS3190 = _M0L8seq__lenS1249->$0;
              int32_t _M0L3valS3191 = _M0L2ltS1259->$0;
              int32_t _M0L8est__capS1278 = _M0L3valS3190 + _M0L3valS3191;
              uint16_t* _M0L3valS3164 = _M0L8seq__bufS1248->$0;
              int32_t _M0L6_2atmpS3163 = Moonbit_array_length(_M0L3valS3164);
              int32_t _M0L1jS1282;
              if (_M0L6_2atmpS3163 < _M0L8est__capS1278) {
                int32_t _M0L6_2atmpS3169 = _M0L8est__capS1278 * 2;
                uint16_t* _M0L8new__bufS1279 =
                  (uint16_t*)moonbit_make_string(_M0L6_2atmpS3169, 0);
                int32_t _M0L1kS1280 = 0;
                uint16_t* _M0L6_2aoldS3644;
                while (1) {
                  int32_t _M0L3valS3165 = _M0L8seq__lenS1249->$0;
                  if (_M0L1kS1280 < _M0L3valS3165) {
                    uint16_t* _M0L3valS3167 = _M0L8seq__bufS1248->$0;
                    int32_t _M0L6_2atmpS3166;
                    int32_t _M0L6_2atmpS3168;
                    if (
                      _M0L1kS1280 < 0
                      || _M0L1kS1280 >= Moonbit_array_length(_M0L3valS3167)
                    ) {
                      #line 176 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L6_2atmpS3166 = (int32_t)_M0L3valS3167[_M0L1kS1280];
                    if (
                      _M0L1kS1280 < 0
                      || _M0L1kS1280
                         >= Moonbit_array_length(_M0L8new__bufS1279)
                    ) {
                      #line 176 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L8new__bufS1279[_M0L1kS1280] = _M0L6_2atmpS3166;
                    _M0L6_2atmpS3168 = _M0L1kS1280 + 1;
                    _M0L1kS1280 = _M0L6_2atmpS3168;
                    continue;
                  }
                  break;
                }
                _M0L6_2aoldS3644 = _M0L8seq__bufS1248->$0;
                moonbit_decref(_M0L6_2aoldS3644);
                _M0L8seq__bufS1248->$0 = _M0L8new__bufS1279;
              }
              _M0L1jS1282 = 0;
              while (1) {
                int32_t _M0L3valS3170 = _M0L2ltS1259->$0;
                if (_M0L1jS1282 < _M0L3valS3170) {
                  int32_t _M0L1cS1283;
                  int32_t _M0L6_2atmpS3174;
                  int32_t _M0L6_2atmpS3173;
                  int32_t _if__result_4348;
                  int32_t _M0L6_2atmpS3189;
                  #line 181 "/home/developer/Documents/1/genbank_io.mbt"
                  _M0L1cS1283
                  = _M0MPC16string10StringView11unsafe__get(_M0L1lS1257, _M0L1jS1282);
                  _M0L6_2atmpS3174 = 97;
                  _M0L6_2atmpS3173 = (uint16_t)_M0L6_2atmpS3174;
                  if (_M0L1cS1283 >= _M0L6_2atmpS3173) {
                    int32_t _M0L6_2atmpS3172 = 122;
                    int32_t _M0L6_2atmpS3171 = (uint16_t)_M0L6_2atmpS3172;
                    _if__result_4348 = _M0L1cS1283 <= _M0L6_2atmpS3171;
                  } else {
                    _if__result_4348 = 0;
                  }
                  if (_if__result_4348) {
                    uint16_t* _M0L3valS3175 = _M0L8seq__bufS1248->$0;
                    int32_t _M0L3valS3176 = _M0L8seq__lenS1249->$0;
                    int32_t _M0L6_2atmpS3178 = (int32_t)_M0L1cS1283;
                    int32_t _M0L6_2atmpS3177;
                    int32_t _M0L3valS3180;
                    int32_t _M0L6_2atmpS3179;
                    if (
                      _M0L6_2atmpS3178 < 0
                      || _M0L6_2atmpS3178
                         >= Moonbit_array_length(_M0FP26biolab8bio__seq16uppercase__table)
                    ) {
                      #line 183 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L6_2atmpS3177
                    = (int32_t)_M0FP26biolab8bio__seq16uppercase__table[
                        _M0L6_2atmpS3178
                      ];
                    if (
                      _M0L3valS3176 < 0
                      || _M0L3valS3176 >= Moonbit_array_length(_M0L3valS3175)
                    ) {
                      #line 183 "/home/developer/Documents/1/genbank_io.mbt"
                      moonbit_panic();
                    }
                    _M0L3valS3175[_M0L3valS3176] = _M0L6_2atmpS3177;
                    _M0L3valS3180 = _M0L8seq__lenS1249->$0;
                    _M0L6_2atmpS3179 = _M0L3valS3180 + 1;
                    _M0L8seq__lenS1249->$0 = _M0L6_2atmpS3179;
                  } else {
                    int32_t _M0L6_2atmpS3184 = 65;
                    int32_t _M0L6_2atmpS3183 = (uint16_t)_M0L6_2atmpS3184;
                    int32_t _if__result_4349;
                    if (_M0L1cS1283 >= _M0L6_2atmpS3183) {
                      int32_t _M0L6_2atmpS3182 = 90;
                      int32_t _M0L6_2atmpS3181 = (uint16_t)_M0L6_2atmpS3182;
                      _if__result_4349 = _M0L1cS1283 <= _M0L6_2atmpS3181;
                    } else {
                      _if__result_4349 = 0;
                    }
                    if (_if__result_4349) {
                      uint16_t* _M0L3valS3185 = _M0L8seq__bufS1248->$0;
                      int32_t _M0L3valS3186 = _M0L8seq__lenS1249->$0;
                      int32_t _M0L3valS3188;
                      int32_t _M0L6_2atmpS3187;
                      if (
                        _M0L3valS3186 < 0
                        || _M0L3valS3186
                           >= Moonbit_array_length(_M0L3valS3185)
                      ) {
                        #line 187 "/home/developer/Documents/1/genbank_io.mbt"
                        moonbit_panic();
                      }
                      _M0L3valS3185[_M0L3valS3186] = _M0L1cS1283;
                      _M0L3valS3188 = _M0L8seq__lenS1249->$0;
                      _M0L6_2atmpS3187 = _M0L3valS3188 + 1;
                      _M0L8seq__lenS1249->$0 = _M0L6_2atmpS3187;
                    }
                  }
                  _M0L6_2atmpS3189 = _M0L1jS1282 + 1;
                  _M0L1jS1282 = _M0L6_2atmpS3189;
                  continue;
                } else {
                  moonbit_decref(_M0L2ltS1259);
                  moonbit_decref(_M0L1lS1257.$0);
                }
                break;
              }
            } else {
              moonbit_decref(_M0L2ltS1259);
              moonbit_decref(_M0L1lS1257.$0);
            }
            _M0L3valS3193 = _M0L1iS1238->$0;
            _M0L6_2atmpS3192 = _M0L3valS3193 + 1;
            _M0L1iS1238->$0 = _M0L6_2atmpS3192;
            continue;
          } else {
            moonbit_decref(_M0L10in__originS1250);
          }
          break;
        }
        _M0L3valS3208 = _M0L8seq__lenS1249->$0;
        _M0L15final__seq__bufS1285
        = (uint16_t*)moonbit_make_string(_M0L3valS3208, 0);
        _M0L1kS1286 = 0;
        while (1) {
          int32_t _M0L3valS3195 = _M0L8seq__lenS1249->$0;
          if (_M0L1kS1286 < _M0L3valS3195) {
            uint16_t* _M0L3valS3197 = _M0L8seq__bufS1248->$0;
            int32_t _M0L6_2atmpS3196;
            int32_t _M0L6_2atmpS3198;
            if (
              _M0L1kS1286 < 0
              || _M0L1kS1286 >= Moonbit_array_length(_M0L3valS3197)
            ) {
              #line 196 "/home/developer/Documents/1/genbank_io.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3196 = (int32_t)_M0L3valS3197[_M0L1kS1286];
            if (
              _M0L1kS1286 < 0
              || _M0L1kS1286
                 >= Moonbit_array_length(_M0L15final__seq__bufS1285)
            ) {
              #line 196 "/home/developer/Documents/1/genbank_io.mbt"
              moonbit_panic();
            }
            _M0L15final__seq__bufS1285[_M0L1kS1286] = _M0L6_2atmpS3196;
            _M0L6_2atmpS3198 = _M0L1kS1286 + 1;
            _M0L1kS1286 = _M0L6_2atmpS3198;
            continue;
          } else {
            moonbit_decref(_M0L8seq__lenS1249);
            moonbit_decref(_M0L8seq__bufS1248);
          }
          break;
        }
        _M0L6_2atmpS3207 = _M0L15final__seq__bufS1285;
        #line 199 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS3200 = _M0MP26biolab8bio__seq3Seq3new(_M0L6_2atmpS3207);
        moonbit_decref(_M0L6_2atmpS3207);
        _M0L8_2afieldS3640 = _M0L2idS1246->$0;
        _M0L6_2acntS4021
        = Moonbit_rc_count(Moonbit_object_header(_M0L2idS1246));
        if (_M0L6_2acntS4021 > 1) {
          int32_t _M0L11_2anew__cntS4022 = _M0L6_2acntS4021 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L2idS1246), _M0L11_2anew__cntS4022);
          moonbit_incref(_M0L8_2afieldS3640);
        } else if (_M0L6_2acntS4021 == 1) {
          #line 198 "/home/developer/Documents/1/genbank_io.mbt"
          moonbit_free(_M0L2idS1246);
        }
        _M0L3valS3201 = _M0L8_2afieldS3640;
        _M0L8_2afieldS3639 = _M0L4nameS1245->$0;
        _M0L6_2acntS4023
        = Moonbit_rc_count(Moonbit_object_header(_M0L4nameS1245));
        if (_M0L6_2acntS4023 > 1) {
          int32_t _M0L11_2anew__cntS4024 = _M0L6_2acntS4023 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L4nameS1245), _M0L11_2anew__cntS4024);
          moonbit_incref(_M0L8_2afieldS3639);
        } else if (_M0L6_2acntS4023 == 1) {
          #line 198 "/home/developer/Documents/1/genbank_io.mbt"
          moonbit_free(_M0L4nameS1245);
        }
        _M0L3valS3202 = _M0L8_2afieldS3639;
        _M0L8_2afieldS3638 = _M0L11descriptionS1247->$0;
        _M0L6_2acntS4025
        = Moonbit_rc_count(Moonbit_object_header(_M0L11descriptionS1247));
        if (_M0L6_2acntS4025 > 1) {
          int32_t _M0L11_2anew__cntS4026 = _M0L6_2acntS4025 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L11descriptionS1247), _M0L11_2anew__cntS4026);
          moonbit_incref(_M0L8_2afieldS3638);
        } else if (_M0L6_2acntS4025 == 1) {
          #line 198 "/home/developer/Documents/1/genbank_io.mbt"
          moonbit_free(_M0L11descriptionS1247);
        }
        _M0L3valS3203 = _M0L8_2afieldS3638;
        _M0L7_2abindS1288
        = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
        _M0L6_2atmpS3206 = _M0L7_2abindS1288;
        _M0L6_2atmpS3205
        = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
          .$0 = _M0L6_2atmpS3206, .$1 = 0, .$2 = 0
        };
        #line 203 "/home/developer/Documents/1/genbank_io.mbt"
        _M0L6_2atmpS3204
        = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS3205, 0ll);
        moonbit_decref(_M0L6_2atmpS3205.$0);
        _M0L6_2atmpS3199
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
        Moonbit_object_header(_M0L6_2atmpS3199)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
        _M0L6_2atmpS3199->$0 = _M0L6_2atmpS3200;
        _M0L6_2atmpS3199->$1 = _M0L3valS3201;
        _M0L6_2atmpS3199->$2 = _M0L3valS3202;
        _M0L6_2atmpS3199->$3 = _M0L3valS3203;
        _M0L6_2atmpS3199->$4 = _M0L6_2atmpS3204;
        #line 198 "/home/developer/Documents/1/genbank_io.mbt"
        _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1233, _M0L6_2atmpS3199);
        moonbit_decref(_M0L6_2atmpS3199);
      } else {
        moonbit_decref(_M0L7trimmedS1244);
        moonbit_decref(_M0L4lineS1242.$0);
      }
      _M0L3valS3214 = _M0L1iS1238->$0;
      _M0L6_2atmpS3213 = _M0L3valS3214 + 1;
      _M0L1iS1238->$0 = _M0L6_2atmpS3213;
      continue;
    } else {
      moonbit_decref(_M0L1iS1238);
      moonbit_decref(_M0L5linesS1234);
    }
    break;
  }
  return _M0L7recordsS1233;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__8(
  
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2974;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2971;
  int32_t* _M0L6_2atmpS2973;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS2972;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1228;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2968;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2967;
  struct moonbit_result_2 _tmp_4351;
  moonbit_string_t _M0L3outS1229;
  moonbit_string_t _M0L6_2atmpS2966;
  struct moonbit_result_0 _result_4353;
  #line 101 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 102 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2974
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_15.data);
  #line 102 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2971
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2974, (moonbit_string_t)moonbit_string_literal_91.data, (moonbit_string_t)moonbit_string_literal_37.data, (moonbit_string_t)moonbit_string_literal_92.data);
  moonbit_decref(_M0L6_2atmpS2974);
  _M0L6_2atmpS2973 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS2973[0] = 40;
  _M0L6_2atmpS2973[1] = 30;
  _M0L6_2atmpS2973[2] = 20;
  _M0L6_2atmpS2973[3] = 10;
  _M0L6_2atmpS2972
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS2972)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS2972->$0 = _M0L6_2atmpS2973;
  _M0L6_2atmpS2972->$1 = 4;
  #line 102 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L3recS1228
  = _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(_M0L6_2atmpS2971, _M0L6_2atmpS2972);
  moonbit_decref(_M0L6_2atmpS2971);
  moonbit_decref(_M0L6_2atmpS2972);
  _M0L6_2atmpS2968
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS2968[0] = _M0L3recS1228;
  _M0L6_2atmpS2967
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L6_2atmpS2967)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 46, 0);
  _M0L6_2atmpS2967->$0 = _M0L6_2atmpS2968;
  _M0L6_2atmpS2967->$1 = 1;
  #line 106 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4351 = _M0FP26biolab8bio__seq12write__fastq(_M0L6_2atmpS2967);
  moonbit_decref(_M0L6_2atmpS2967);
  if (_tmp_4351.tag) {
    moonbit_string_t const _M0L5_2aokS2969 = _tmp_4351.data.ok;
    _M0L3outS1229 = _M0L5_2aokS2969;
  } else {
    void* const _M0L6_2aerrS2970 = _tmp_4351.data.err;
    struct moonbit_result_0 _result_4352;
    _result_4352.tag = 0;
    _result_4352.data.err = _M0L6_2aerrS2970;
    return _result_4352;
  }
  _M0L6_2atmpS2966 = 0;
  #line 108 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4353
  = _M0FPB10assert__eqGsE(_M0L3outS1229, (moonbit_string_t)moonbit_string_literal_93.data, _M0L6_2atmpS2966, (moonbit_string_t)moonbit_string_literal_94.data);
  moonbit_decref(_M0L3outS1229);
  if (_M0L6_2atmpS2966) {
    moonbit_decref(_M0L6_2atmpS2966);
  }
  return _result_4353;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord19set__phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS1216,
  struct _M0TPB5ArrayGiE* _M0L5qualsS1227
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1215;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2965;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2960;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2964;
  int32_t _M0L6_2atmpS2963;
  int32_t _M0L6_2atmpS2962;
  int64_t _M0L6_2atmpS2961;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L16new__annotationsS1214;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2955;
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0L5_2aitS1217;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2956;
  moonbit_string_t _M0L2idS2957;
  moonbit_string_t _M0L4nameS2958;
  moonbit_string_t _M0L11descriptionS2959;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_4356;
  #line 61 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS1215 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2965 = _M0L7_2abindS1215;
  _M0L6_2atmpS2960
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2965, .$1 = 0, .$2 = 0
  };
  _M0L19letter__annotationsS2964 = _M0L4selfS1216->$4;
  moonbit_incref(_M0L19letter__annotationsS2964);
  #line 67 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS2963
  = _M0MPB3Map6lengthGsRPB5ArrayGiEE(_M0L19letter__annotationsS2964);
  moonbit_decref(_M0L19letter__annotationsS2964);
  _M0L6_2atmpS2962 = _M0L6_2atmpS2963 + 1;
  _M0L6_2atmpS2961 = (int64_t)_M0L6_2atmpS2962;
  #line 65 "/home/developer/Documents/1/seq_record.mbt"
  _M0L16new__annotationsS1214
  = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2960, _M0L6_2atmpS2961);
  moonbit_decref(_M0L6_2atmpS2960.$0);
  _M0L19letter__annotationsS2955 = _M0L4selfS1216->$4;
  moonbit_incref(_M0L19letter__annotationsS2955);
  #line 65 "/home/developer/Documents/1/seq_record.mbt"
  _M0L5_2aitS1217
  = _M0MPB3Map5iter2GsRPB5ArrayGiEE(_M0L19letter__annotationsS2955);
  moonbit_decref(_M0L19letter__annotationsS2955);
  while (1) {
    moonbit_string_t _M0L1kS1219;
    struct _M0TPB5ArrayGiE* _M0L1vS1220;
    struct _M0TUsRPB5ArrayGiEE* _M0L7_2abindS1222;
    #line 69 "/home/developer/Documents/1/seq_record.mbt"
    _M0L7_2abindS1222 = _M0MPB5Iter24nextGsRPB5ArrayGiEE(_M0L5_2aitS1217);
    if (_M0L7_2abindS1222 == 0) {
      if (_M0L7_2abindS1222) {
        moonbit_decref(_M0L7_2abindS1222);
      }
      moonbit_decref(_M0L5_2aitS1217);
    } else {
      struct _M0TUsRPB5ArrayGiEE* _M0L7_2aSomeS1223 = _M0L7_2abindS1222;
      struct _M0TUsRPB5ArrayGiEE* _M0L4_2axS1224 = _M0L7_2aSomeS1223;
      moonbit_string_t _M0L4_2akS1225 = _M0L4_2axS1224->$0;
      struct _M0TPB5ArrayGiE* _M0L8_2afieldS3656 = _M0L4_2axS1224->$1;
      int32_t _M0L6_2acntS4027 =
        Moonbit_rc_count(Moonbit_object_header(_M0L4_2axS1224));
      struct _M0TPB5ArrayGiE* _M0L4_2avS1226;
      if (_M0L6_2acntS4027 > 1) {
        int32_t _M0L11_2anew__cntS4028 = _M0L6_2acntS4027 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L4_2axS1224), _M0L11_2anew__cntS4028);
        moonbit_incref(_M0L8_2afieldS3656);
        moonbit_incref(_M0L4_2akS1225);
      } else if (_M0L6_2acntS4027 == 1) {
        #line 69 "/home/developer/Documents/1/seq_record.mbt"
        moonbit_free(_M0L4_2axS1224);
      }
      _M0L4_2avS1226 = _M0L8_2afieldS3656;
      _M0L1kS1219 = _M0L4_2akS1225;
      _M0L1vS1220 = _M0L4_2avS1226;
      goto join_1218;
    }
    goto joinlet_4355;
    join_1218:;
    #line 70 "/home/developer/Documents/1/seq_record.mbt"
    _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L16new__annotationsS1214, _M0L1kS1219, _M0L1vS1220);
    moonbit_decref(_M0L1kS1219);
    moonbit_decref(_M0L1vS1220);
    continue;
    joinlet_4355:;
    break;
  }
  #line 72 "/home/developer/Documents/1/seq_record.mbt"
  _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L16new__annotationsS1214, (moonbit_string_t)moonbit_string_literal_95.data, _M0L5qualsS1227);
  _M0L3seqS2956 = _M0L4selfS1216->$0;
  _M0L2idS2957 = _M0L4selfS1216->$1;
  _M0L4nameS2958 = _M0L4selfS1216->$2;
  _M0L11descriptionS2959 = _M0L4selfS1216->$3;
  moonbit_incref(_M0L3seqS2956);
  moonbit_incref(_M0L2idS2957);
  moonbit_incref(_M0L4nameS2958);
  moonbit_incref(_M0L11descriptionS2959);
  _block_4356
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_4356)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
  _block_4356->$0 = _M0L3seqS2956;
  _block_4356->$1 = _M0L2idS2957;
  _block_4356->$2 = _M0L4nameS2958;
  _block_4356->$3 = _M0L11descriptionS2959;
  _block_4356->$4 = _M0L16new__annotationsS1214;
  return _block_4356;
}

struct moonbit_result_2 _M0FP26biolab8bio__seq12write__fastq(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1200
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1198;
  int32_t _M0L7_2abindS1199;
  int32_t _M0L2__S1201;
  moonbit_string_t _M0L6_2atmpS2954;
  struct moonbit_result_2 _result_4362;
  #line 123 "/home/developer/Documents/1/fastq_io.mbt"
  #line 124 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L3bufS1198 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS1199 = _M0L7recordsS1200->$1;
  _M0L2__S1201 = 0;
  while (1) {
    if (_M0L2__S1201 < _M0L7_2abindS1199) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2953 =
        _M0L7recordsS1200->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1202 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2953[_M0L2__S1201];
      struct _M0TPB5ArrayGiE* _M0L1qS1205;
      struct _M0TPB5ArrayGiE* _M0L5qualsS1203;
      struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2951 =
        _M0L3recS1202->$4;
      struct _M0TPB5ArrayGiE* _M0L7_2abindS1206;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2947;
      moonbit_string_t _M0L3seqS1209;
      int32_t _M0L6_2atmpS2924;
      int32_t _M0L6_2atmpS2925;
      moonbit_string_t _M0L2idS2938;
      moonbit_string_t _M0L8_2afieldS3660;
      int32_t _M0L6_2acntS4041;
      moonbit_string_t _M0L11descriptionS2939;
      moonbit_string_t _M0L6_2atmpS2937;
      int32_t _M0L6_2atmpS2946;
      uint16_t* _M0L4qbufS1210;
      int32_t _M0L1jS1211;
      moonbit_string_t _M0L6_2atmpS2945;
      int32_t _M0L6_2atmpS2952;
      moonbit_incref(_M0L19letter__annotationsS2951);
      moonbit_incref(_M0L3recS1202);
      #line 126 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L7_2abindS1206
      = _M0MPB3Map3getGsRPB5ArrayGiEE(_M0L19letter__annotationsS2951, (moonbit_string_t)moonbit_string_literal_95.data);
      moonbit_decref(_M0L19letter__annotationsS2951);
      if (_M0L7_2abindS1206 == 0) {
        moonbit_string_t _M0L8_2afieldS3664;
        int32_t _M0L6_2acntS4029;
        moonbit_string_t _M0L2idS2950;
        moonbit_string_t _M0L6_2atmpS2949;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2948;
        struct moonbit_result_2 _result_4359;
        if (_M0L7_2abindS1206) {
          moonbit_decref(_M0L7_2abindS1206);
        }
        moonbit_decref(_M0L3bufS1198);
        _M0L8_2afieldS3664 = _M0L3recS1202->$1;
        _M0L6_2acntS4029
        = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1202));
        if (_M0L6_2acntS4029 > 1) {
          int32_t _M0L11_2anew__cntS4034 = _M0L6_2acntS4029 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1202), _M0L11_2anew__cntS4034);
          moonbit_incref(_M0L8_2afieldS3664);
        } else if (_M0L6_2acntS4029 == 1) {
          struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4033 =
            _M0L3recS1202->$4;
          moonbit_string_t _M0L8_2afieldS4032;
          moonbit_string_t _M0L8_2afieldS4031;
          struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4030;
          moonbit_decref(_M0L8_2afieldS4033);
          _M0L8_2afieldS4032 = _M0L3recS1202->$3;
          moonbit_decref(_M0L8_2afieldS4032);
          _M0L8_2afieldS4031 = _M0L3recS1202->$2;
          moonbit_decref(_M0L8_2afieldS4031);
          _M0L8_2afieldS4030 = _M0L3recS1202->$0;
          moonbit_decref(_M0L8_2afieldS4030);
          #line 130 "/home/developer/Documents/1/fastq_io.mbt"
          moonbit_free(_M0L3recS1202);
        }
        _M0L2idS2950 = _M0L8_2afieldS3664;
        #line 130 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2949
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_96.data, _M0L2idS2950);
        moonbit_decref(_M0L2idS2950);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2948
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2948)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2948)->$0
        = _M0L6_2atmpS2949;
        _result_4359.tag = 0;
        _result_4359.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2948;
        return _result_4359;
      } else {
        struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1207 = _M0L7_2abindS1206;
        struct _M0TPB5ArrayGiE* _M0L4_2aqS1208 = _M0L7_2aSomeS1207;
        _M0L1qS1205 = _M0L4_2aqS1208;
        goto join_1204;
      }
      goto joinlet_4358;
      join_1204:;
      _M0L5qualsS1203 = _M0L1qS1205;
      joinlet_4358:;
      _M0L3seqS2947 = _M0L3recS1202->$0;
      moonbit_incref(_M0L3seqS2947);
      #line 133 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L3seqS1209 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2947);
      moonbit_decref(_M0L3seqS2947);
      _M0L6_2atmpS2924 = Moonbit_array_length(_M0L3seqS1209);
      #line 134 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2925 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS1203);
      if (_M0L6_2atmpS2924 != _M0L6_2atmpS2925) {
        int32_t _M0L6_2atmpS2936;
        moonbit_string_t _M0L6_2atmpS2935;
        moonbit_string_t _M0L6_2atmpS2934;
        moonbit_string_t _M0L6_2atmpS2931;
        int32_t _M0L6_2atmpS2933;
        moonbit_string_t _M0L6_2atmpS2932;
        moonbit_string_t _M0L6_2atmpS2930;
        moonbit_string_t _M0L6_2atmpS2928;
        moonbit_string_t _M0L8_2afieldS3662;
        int32_t _M0L6_2acntS4035;
        moonbit_string_t _M0L2idS2929;
        moonbit_string_t _M0L6_2atmpS2927;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2926;
        struct moonbit_result_2 _result_4360;
        moonbit_decref(_M0L3bufS1198);
        _M0L6_2atmpS2936 = Moonbit_array_length(_M0L3seqS1209);
        moonbit_decref(_M0L3seqS1209);
        #line 137 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2935
        = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2936, 10);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2934
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_97.data, _M0L6_2atmpS2935);
        moonbit_decref(_M0L6_2atmpS2935);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2931
        = moonbit_add_string(_M0L6_2atmpS2934, (moonbit_string_t)moonbit_string_literal_98.data);
        moonbit_decref(_M0L6_2atmpS2934);
        #line 139 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2933 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS1203);
        moonbit_decref(_M0L5qualsS1203);
        #line 139 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2932
        = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2933, 10);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2930
        = moonbit_add_string(_M0L6_2atmpS2931, _M0L6_2atmpS2932);
        moonbit_decref(_M0L6_2atmpS2932);
        moonbit_decref(_M0L6_2atmpS2931);
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2928
        = moonbit_add_string(_M0L6_2atmpS2930, (moonbit_string_t)moonbit_string_literal_99.data);
        moonbit_decref(_M0L6_2atmpS2930);
        _M0L8_2afieldS3662 = _M0L3recS1202->$1;
        _M0L6_2acntS4035
        = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1202));
        if (_M0L6_2acntS4035 > 1) {
          int32_t _M0L11_2anew__cntS4040 = _M0L6_2acntS4035 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1202), _M0L11_2anew__cntS4040);
          moonbit_incref(_M0L8_2afieldS3662);
        } else if (_M0L6_2acntS4035 == 1) {
          struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4039 =
            _M0L3recS1202->$4;
          moonbit_string_t _M0L8_2afieldS4038;
          moonbit_string_t _M0L8_2afieldS4037;
          struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4036;
          moonbit_decref(_M0L8_2afieldS4039);
          _M0L8_2afieldS4038 = _M0L3recS1202->$3;
          moonbit_decref(_M0L8_2afieldS4038);
          _M0L8_2afieldS4037 = _M0L3recS1202->$2;
          moonbit_decref(_M0L8_2afieldS4037);
          _M0L8_2afieldS4036 = _M0L3recS1202->$0;
          moonbit_decref(_M0L8_2afieldS4036);
          #line 141 "/home/developer/Documents/1/fastq_io.mbt"
          moonbit_free(_M0L3recS1202);
        }
        _M0L2idS2929 = _M0L8_2afieldS3662;
        #line 136 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2927 = moonbit_add_string(_M0L6_2atmpS2928, _M0L2idS2929);
        moonbit_decref(_M0L2idS2929);
        moonbit_decref(_M0L6_2atmpS2928);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2926
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2926)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2926)->$0
        = _M0L6_2atmpS2927;
        _result_4360.tag = 0;
        _result_4360.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2926;
        return _result_4360;
      }
      #line 145 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1198, 64);
      _M0L2idS2938 = _M0L3recS1202->$1;
      _M0L8_2afieldS3660 = _M0L3recS1202->$3;
      _M0L6_2acntS4041
      = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1202));
      if (_M0L6_2acntS4041 > 1) {
        int32_t _M0L11_2anew__cntS4045 = _M0L6_2acntS4041 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1202), _M0L11_2anew__cntS4045);
        moonbit_incref(_M0L8_2afieldS3660);
        moonbit_incref(_M0L2idS2938);
      } else if (_M0L6_2acntS4041 == 1) {
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4044 =
          _M0L3recS1202->$4;
        moonbit_string_t _M0L8_2afieldS4043;
        struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4042;
        moonbit_decref(_M0L8_2afieldS4044);
        _M0L8_2afieldS4043 = _M0L3recS1202->$2;
        moonbit_decref(_M0L8_2afieldS4043);
        _M0L8_2afieldS4042 = _M0L3recS1202->$0;
        moonbit_decref(_M0L8_2afieldS4042);
        #line 146 "/home/developer/Documents/1/fastq_io.mbt"
        moonbit_free(_M0L3recS1202);
      }
      _M0L11descriptionS2939 = _M0L8_2afieldS3660;
      #line 146 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2937
      = _M0FP26biolab8bio__seq19build__fasta__title(_M0L2idS2938, _M0L11descriptionS2939);
      moonbit_decref(_M0L2idS2938);
      moonbit_decref(_M0L11descriptionS2939);
      #line 146 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1198, _M0L6_2atmpS2937);
      moonbit_decref(_M0L6_2atmpS2937);
      #line 147 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1198, 10);
      #line 149 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1198, _M0L3seqS1209);
      moonbit_decref(_M0L3seqS1209);
      #line 150 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1198, 10);
      #line 152 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1198, (moonbit_string_t)moonbit_string_literal_100.data);
      #line 154 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2946 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS1203);
      _M0L4qbufS1210 = (uint16_t*)moonbit_make_string(_M0L6_2atmpS2946, 0);
      _M0L1jS1211 = 0;
      while (1) {
        int32_t _M0L6_2atmpS2940;
        #line 155 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2940 = _M0MPC15array5Array6lengthGiE(_M0L5qualsS1203);
        if (_M0L1jS1211 < _M0L6_2atmpS2940) {
          int32_t _M0L6_2atmpS2943;
          int32_t _M0L6_2atmpS2942;
          int32_t _M0L6_2atmpS2941;
          int32_t _M0L6_2atmpS2944;
          #line 156 "/home/developer/Documents/1/fastq_io.mbt"
          _M0L6_2atmpS2943
          = _M0MPC15array5Array2atGiE(_M0L5qualsS1203, _M0L1jS1211);
          _M0L6_2atmpS2942 = _M0L6_2atmpS2943 + 33;
          _M0L6_2atmpS2941 = (uint16_t)_M0L6_2atmpS2942;
          if (
            _M0L1jS1211 < 0
            || _M0L1jS1211 >= Moonbit_array_length(_M0L4qbufS1210)
          ) {
            #line 156 "/home/developer/Documents/1/fastq_io.mbt"
            moonbit_panic();
          }
          _M0L4qbufS1210[_M0L1jS1211] = _M0L6_2atmpS2941;
          _M0L6_2atmpS2944 = _M0L1jS1211 + 1;
          _M0L1jS1211 = _M0L6_2atmpS2944;
          continue;
        } else {
          moonbit_decref(_M0L5qualsS1203);
        }
        break;
      }
      _M0L6_2atmpS2945 = _M0L4qbufS1210;
      #line 158 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1198, _M0L6_2atmpS2945);
      moonbit_decref(_M0L6_2atmpS2945);
      #line 159 "/home/developer/Documents/1/fastq_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1198, 10);
      _M0L6_2atmpS2952 = _M0L2__S1201 + 1;
      _M0L2__S1201 = _M0L6_2atmpS2952;
      continue;
    }
    break;
  }
  #line 161 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L6_2atmpS2954 = _M0MPB13StringBuilder10to__string(_M0L3bufS1198);
  moonbit_decref(_M0L3bufS1198);
  _result_4362.tag = 1;
  _result_4362.data.ok = _M0L6_2atmpS2954;
  return _result_4362;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__7(
  
) {
  moonbit_string_t _M0L9long__seqS1193;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2923;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1194;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2922;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2921;
  moonbit_string_t _M0L3outS1195;
  moonbit_string_t _M0L6_2atmpS2920;
  moonbit_string_t _M0L6_2atmpS2919;
  moonbit_string_t _M0L6_2atmpS2917;
  moonbit_string_t _M0L6_2atmpS2918;
  moonbit_string_t _M0L6_2atmpS2916;
  moonbit_string_t _M0L8expectedS1196;
  moonbit_string_t _M0L6_2atmpS2909;
  struct moonbit_result_0 _tmp_4363;
  moonbit_string_t _M0L7_2abindS1197;
  int32_t _M0L6_2atmpS2915;
  struct _M0TPC16string10StringView _M0L6_2atmpS2914;
  int32_t _M0L6_2atmpS2912;
  void* _M0L4NoneS2913;
  struct moonbit_result_0 _result_4365;
  #line 89 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 92 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L9long__seqS1193
  = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_101.data, 70);
  #line 93 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2923 = _M0MP26biolab8bio__seq3Seq3new(_M0L9long__seqS1193);
  moonbit_decref(_M0L9long__seqS1193);
  #line 93 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L3recS1194
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2923, (moonbit_string_t)moonbit_string_literal_102.data, (moonbit_string_t)moonbit_string_literal_37.data, (moonbit_string_t)moonbit_string_literal_102.data);
  moonbit_decref(_M0L6_2atmpS2923);
  _M0L6_2atmpS2922
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS2922[0] = _M0L3recS1194;
  _M0L6_2atmpS2921
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L6_2atmpS2921)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 46, 0);
  _M0L6_2atmpS2921->$0 = _M0L6_2atmpS2922;
  _M0L6_2atmpS2921->$1 = 1;
  #line 94 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L3outS1195 = _M0FP26biolab8bio__seq12write__fasta(_M0L6_2atmpS2921);
  moonbit_decref(_M0L6_2atmpS2921);
  #line 95 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2920
  = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_101.data, 60);
  #line 95 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2919
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_103.data, _M0L6_2atmpS2920);
  moonbit_decref(_M0L6_2atmpS2920);
  #line 95 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2917
  = moonbit_add_string(_M0L6_2atmpS2919, (moonbit_string_t)moonbit_string_literal_89.data);
  moonbit_decref(_M0L6_2atmpS2919);
  #line 95 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2918
  = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_101.data, 10);
  #line 95 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2916 = moonbit_add_string(_M0L6_2atmpS2917, _M0L6_2atmpS2918);
  moonbit_decref(_M0L6_2atmpS2918);
  moonbit_decref(_M0L6_2atmpS2917);
  #line 95 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L8expectedS1196
  = moonbit_add_string(_M0L6_2atmpS2916, (moonbit_string_t)moonbit_string_literal_89.data);
  moonbit_decref(_M0L6_2atmpS2916);
  _M0L6_2atmpS2909 = 0;
  #line 96 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4363
  = _M0FPB10assert__eqGsE(_M0L3outS1195, _M0L8expectedS1196, _M0L6_2atmpS2909, (moonbit_string_t)moonbit_string_literal_104.data);
  moonbit_decref(_M0L8expectedS1196);
  if (_M0L6_2atmpS2909) {
    moonbit_decref(_M0L6_2atmpS2909);
  }
  if (_tmp_4363.tag) {
    int32_t const _M0L5_2aokS2910 = _tmp_4363.data.ok;
  } else {
    void* const _M0L6_2aerrS2911 = _tmp_4363.data.err;
    struct moonbit_result_0 _result_4364;
    moonbit_decref(_M0L3outS1195);
    _result_4364.tag = 0;
    _result_4364.data.err = _M0L6_2aerrS2911;
    return _result_4364;
  }
  _M0L7_2abindS1197 = (moonbit_string_t)moonbit_string_literal_103.data;
  _M0L6_2atmpS2915 = Moonbit_array_length(_M0L7_2abindS1197);
  _M0L6_2atmpS2914
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1197, .$1 = 0, .$2 = _M0L6_2atmpS2915
  };
  #line 97 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2912
  = _M0MPC16string6String11has__prefix(_M0L3outS1195, _M0L6_2atmpS2914);
  moonbit_decref(_M0L3outS1195);
  moonbit_decref(_M0L6_2atmpS2914.$0);
  _M0L4NoneS2913
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 97 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4365
  = _M0FPB12assert__true(_M0L6_2atmpS2912, _M0L4NoneS2913, (moonbit_string_t)moonbit_string_literal_105.data);
  moonbit_decref(_M0L4NoneS2913);
  return _result_4365;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__6(
  
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2908;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2905;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2907;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2906;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2904;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1188;
  moonbit_string_t _M0L3outS1189;
  moonbit_string_t _M0L7_2abindS1190;
  int32_t _M0L6_2atmpS2888;
  struct _M0TPC16string10StringView _M0L6_2atmpS2887;
  int32_t _M0L6_2atmpS2885;
  void* _M0L4NoneS2886;
  struct moonbit_result_0 _tmp_4366;
  moonbit_string_t _M0L7_2abindS1191;
  int32_t _M0L6_2atmpS2894;
  struct _M0TPC16string10StringView _M0L6_2atmpS2893;
  int32_t _M0L6_2atmpS2891;
  void* _M0L4NoneS2892;
  struct moonbit_result_0 _tmp_4368;
  moonbit_string_t _M0L7_2abindS1192;
  int32_t _M0L6_2atmpS2900;
  struct _M0TPC16string10StringView _M0L6_2atmpS2899;
  int32_t _M0L6_2atmpS2897;
  void* _M0L4NoneS2898;
  struct moonbit_result_0 _tmp_4370;
  moonbit_string_t _M0L6_2atmpS2903;
  struct moonbit_result_0 _result_4372;
  #line 76 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 78 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2908
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_15.data);
  #line 78 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2905
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2908, (moonbit_string_t)moonbit_string_literal_106.data, (moonbit_string_t)moonbit_string_literal_37.data, (moonbit_string_t)moonbit_string_literal_107.data);
  moonbit_decref(_M0L6_2atmpS2908);
  #line 79 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2907
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_29.data);
  #line 79 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2906
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2907, (moonbit_string_t)moonbit_string_literal_108.data, (moonbit_string_t)moonbit_string_literal_37.data, (moonbit_string_t)moonbit_string_literal_108.data);
  moonbit_decref(_M0L6_2atmpS2907);
  _M0L6_2atmpS2904
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(2);
  _M0L6_2atmpS2904[0] = _M0L6_2atmpS2905;
  _M0L6_2atmpS2904[1] = _M0L6_2atmpS2906;
  _M0L7recordsS1188
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L7recordsS1188)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 46, 0);
  _M0L7recordsS1188->$0 = _M0L6_2atmpS2904;
  _M0L7recordsS1188->$1 = 2;
  #line 81 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L3outS1189 = _M0FP26biolab8bio__seq12write__fasta(_M0L7recordsS1188);
  moonbit_decref(_M0L7recordsS1188);
  _M0L7_2abindS1190 = (moonbit_string_t)moonbit_string_literal_109.data;
  _M0L6_2atmpS2888 = Moonbit_array_length(_M0L7_2abindS1190);
  _M0L6_2atmpS2887
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1190, .$1 = 0, .$2 = _M0L6_2atmpS2888
  };
  #line 82 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2885
  = _M0MPC16string6String11has__prefix(_M0L3outS1189, _M0L6_2atmpS2887);
  moonbit_decref(_M0L6_2atmpS2887.$0);
  _M0L4NoneS2886
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 82 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4366
  = _M0FPB12assert__true(_M0L6_2atmpS2885, _M0L4NoneS2886, (moonbit_string_t)moonbit_string_literal_110.data);
  moonbit_decref(_M0L4NoneS2886);
  if (_tmp_4366.tag) {
    int32_t const _M0L5_2aokS2889 = _tmp_4366.data.ok;
  } else {
    void* const _M0L6_2aerrS2890 = _tmp_4366.data.err;
    struct moonbit_result_0 _result_4367;
    moonbit_decref(_M0L3outS1189);
    _result_4367.tag = 0;
    _result_4367.data.err = _M0L6_2aerrS2890;
    return _result_4367;
  }
  _M0L7_2abindS1191 = (moonbit_string_t)moonbit_string_literal_15.data;
  _M0L6_2atmpS2894 = Moonbit_array_length(_M0L7_2abindS1191);
  _M0L6_2atmpS2893
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1191, .$1 = 0, .$2 = _M0L6_2atmpS2894
  };
  #line 83 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2891
  = _M0MPC16string6String8contains(_M0L3outS1189, _M0L6_2atmpS2893);
  moonbit_decref(_M0L6_2atmpS2893.$0);
  _M0L4NoneS2892
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 83 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4368
  = _M0FPB12assert__true(_M0L6_2atmpS2891, _M0L4NoneS2892, (moonbit_string_t)moonbit_string_literal_111.data);
  moonbit_decref(_M0L4NoneS2892);
  if (_tmp_4368.tag) {
    int32_t const _M0L5_2aokS2895 = _tmp_4368.data.ok;
  } else {
    void* const _M0L6_2aerrS2896 = _tmp_4368.data.err;
    struct moonbit_result_0 _result_4369;
    moonbit_decref(_M0L3outS1189);
    _result_4369.tag = 0;
    _result_4369.data.err = _M0L6_2aerrS2896;
    return _result_4369;
  }
  _M0L7_2abindS1192 = (moonbit_string_t)moonbit_string_literal_29.data;
  _M0L6_2atmpS2900 = Moonbit_array_length(_M0L7_2abindS1192);
  _M0L6_2atmpS2899
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1192, .$1 = 0, .$2 = _M0L6_2atmpS2900
  };
  #line 84 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2897
  = _M0MPC16string6String8contains(_M0L3outS1189, _M0L6_2atmpS2899);
  moonbit_decref(_M0L6_2atmpS2899.$0);
  _M0L4NoneS2898
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 84 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4370
  = _M0FPB12assert__true(_M0L6_2atmpS2897, _M0L4NoneS2898, (moonbit_string_t)moonbit_string_literal_112.data);
  moonbit_decref(_M0L4NoneS2898);
  if (_tmp_4370.tag) {
    int32_t const _M0L5_2aokS2901 = _tmp_4370.data.ok;
  } else {
    void* const _M0L6_2aerrS2902 = _tmp_4370.data.err;
    struct moonbit_result_0 _result_4371;
    moonbit_decref(_M0L3outS1189);
    _result_4371.tag = 0;
    _result_4371.data.err = _M0L6_2aerrS2902;
    return _result_4371;
  }
  _M0L6_2atmpS2903 = 0;
  #line 85 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4372
  = _M0FPB10assert__eqGsE(_M0L3outS1189, (moonbit_string_t)moonbit_string_literal_113.data, _M0L6_2atmpS2903, (moonbit_string_t)moonbit_string_literal_114.data);
  moonbit_decref(_M0L3outS1189);
  if (_M0L6_2atmpS2903) {
    moonbit_decref(_M0L6_2atmpS2903);
  }
  return _result_4372;
}

moonbit_string_t _M0FP26biolab8bio__seq12write__fasta(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1178
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1176;
  int32_t _M0L7_2abindS1177;
  int32_t _M0L2__S1179;
  moonbit_string_t _result_4375;
  #line 99 "/home/developer/Documents/1/fasta_io.mbt"
  #line 100 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L3bufS1176 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS1177 = _M0L7recordsS1178->$1;
  _M0L2__S1179 = 0;
  while (1) {
    if (_M0L2__S1179 < _M0L7_2abindS1177) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2884 =
        _M0L7recordsS1178->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1180 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2884[_M0L2__S1179];
      moonbit_string_t _M0L2idS2881 = _M0L3recS1180->$1;
      moonbit_string_t _M0L11descriptionS2882 = _M0L3recS1180->$3;
      moonbit_string_t _M0L5titleS1181;
      struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3668;
      int32_t _M0L6_2acntS4046;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2880;
      moonbit_string_t _M0L3seqS1182;
      int32_t _M0L1nS1183;
      struct _M0TPB8MutLocalGiE* _M0L1iS1184;
      int32_t _M0L6_2atmpS2883;
      moonbit_incref(_M0L11descriptionS2882);
      moonbit_incref(_M0L2idS2881);
      moonbit_incref(_M0L3recS1180);
      #line 107 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L5titleS1181
      = _M0FP26biolab8bio__seq19build__fasta__title(_M0L2idS2881, _M0L11descriptionS2882);
      moonbit_decref(_M0L2idS2881);
      moonbit_decref(_M0L11descriptionS2882);
      #line 108 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1176, 62);
      #line 109 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1176, _M0L5titleS1181);
      moonbit_decref(_M0L5titleS1181);
      #line 110 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1176, 10);
      _M0L8_2afieldS3668 = _M0L3recS1180->$0;
      _M0L6_2acntS4046
      = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1180));
      if (_M0L6_2acntS4046 > 1) {
        int32_t _M0L11_2anew__cntS4051 = _M0L6_2acntS4046 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1180), _M0L11_2anew__cntS4051);
        moonbit_incref(_M0L8_2afieldS3668);
      } else if (_M0L6_2acntS4046 == 1) {
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4050 =
          _M0L3recS1180->$4;
        moonbit_string_t _M0L8_2afieldS4049;
        moonbit_string_t _M0L8_2afieldS4048;
        moonbit_string_t _M0L8_2afieldS4047;
        moonbit_decref(_M0L8_2afieldS4050);
        _M0L8_2afieldS4049 = _M0L3recS1180->$3;
        moonbit_decref(_M0L8_2afieldS4049);
        _M0L8_2afieldS4048 = _M0L3recS1180->$2;
        moonbit_decref(_M0L8_2afieldS4048);
        _M0L8_2afieldS4047 = _M0L3recS1180->$1;
        moonbit_decref(_M0L8_2afieldS4047);
        #line 112 "/home/developer/Documents/1/fasta_io.mbt"
        moonbit_free(_M0L3recS1180);
      }
      _M0L3seqS2880 = _M0L8_2afieldS3668;
      #line 112 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L3seqS1182 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2880);
      moonbit_decref(_M0L3seqS2880);
      _M0L1nS1183 = Moonbit_array_length(_M0L3seqS1182);
      _M0L1iS1184
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1iS1184)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1iS1184->$0 = 0;
      while (1) {
        int32_t _M0L3valS2872 = _M0L1iS1184->$0;
        if (_M0L3valS2872 < _M0L1nS1183) {
          int32_t _M0L3valS2878 = _M0L1iS1184->$0;
          int32_t _M0L6_2atmpS2877 =
            _M0L3valS2878 + _M0FP26biolab8bio__seq18fasta__line__width;
          int32_t _M0L3endS1185;
          int32_t _M0L3valS2875;
          int64_t _M0L6_2atmpS2876;
          struct _M0TPC16string10StringView _M0L6_2atmpS2874;
          moonbit_string_t _M0L6_2atmpS2873;
          if (_M0L6_2atmpS2877 < _M0L1nS1183) {
            int32_t _M0L3valS2879 = _M0L1iS1184->$0;
            _M0L3endS1185
            = _M0L3valS2879 + _M0FP26biolab8bio__seq18fasta__line__width;
          } else {
            _M0L3endS1185 = _M0L1nS1183;
          }
          _M0L3valS2875 = _M0L1iS1184->$0;
          _M0L6_2atmpS2876 = (int64_t)_M0L3endS1185;
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2874
          = _M0MPC16string6String11sub_2einner(_M0L3seqS1182, _M0L3valS2875, _M0L6_2atmpS2876);
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2873
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2874);
          moonbit_decref(_M0L6_2atmpS2874.$0);
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1176, _M0L6_2atmpS2873);
          moonbit_decref(_M0L6_2atmpS2873);
          #line 118 "/home/developer/Documents/1/fasta_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1176, 10);
          _M0L1iS1184->$0 = _M0L3endS1185;
          continue;
        } else {
          moonbit_decref(_M0L1iS1184);
          moonbit_decref(_M0L3seqS1182);
        }
        break;
      }
      _M0L6_2atmpS2883 = _M0L2__S1179 + 1;
      _M0L2__S1179 = _M0L6_2atmpS2883;
      continue;
    }
    break;
  }
  #line 122 "/home/developer/Documents/1/fasta_io.mbt"
  _result_4375 = _M0MPB13StringBuilder10to__string(_M0L3bufS1176);
  moonbit_decref(_M0L3bufS1176);
  return _result_4375;
}

moonbit_string_t _M0FP26biolab8bio__seq19build__fasta__title(
  moonbit_string_t _M0L2idS1175,
  moonbit_string_t _M0L11descriptionS1174
) {
  int32_t _M0L6_2atmpS2868;
  #line 129 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2868 = Moonbit_array_length(_M0L11descriptionS1174);
  if (_M0L6_2atmpS2868 == 0) {
    moonbit_incref(_M0L2idS1175);
    return _M0L2idS1175;
  } else {
    int32_t _M0L6_2atmpS2870 = Moonbit_array_length(_M0L2idS1175);
    struct _M0TPC16string10StringView _M0L6_2atmpS2869;
    int32_t _result_4376;
    moonbit_incref(_M0L2idS1175);
    _M0L6_2atmpS2869
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L2idS1175, .$1 = 0, .$2 = _M0L6_2atmpS2870
    };
    #line 132 "/home/developer/Documents/1/fasta_io.mbt"
    _result_4376
    = _M0MPC16string6String11has__prefix(_M0L11descriptionS1174, _M0L6_2atmpS2869);
    moonbit_decref(_M0L6_2atmpS2869.$0);
    if (_result_4376) {
      moonbit_incref(_M0L11descriptionS1174);
      return _M0L11descriptionS1174;
    } else {
      moonbit_string_t _M0L6_2atmpS2871;
      moonbit_string_t _result_4377;
      #line 135 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L6_2atmpS2871
      = moonbit_add_string(_M0L2idS1175, (moonbit_string_t)moonbit_string_literal_115.data);
      #line 135 "/home/developer/Documents/1/fasta_io.mbt"
      _result_4377
      = moonbit_add_string(_M0L6_2atmpS2871, _M0L11descriptionS1174);
      moonbit_decref(_M0L6_2atmpS2871);
      return _result_4377;
    }
  }
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__5(
  
) {
  struct moonbit_result_1 _tmp_4378;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1161;
  int32_t _M0L6_2atmpS2822;
  moonbit_string_t _M0L6_2atmpS2823;
  struct moonbit_result_0 _tmp_4380;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2828;
  moonbit_string_t _M0L8_2afieldS3676;
  int32_t _M0L6_2acntS4052;
  moonbit_string_t _M0L2idS2826;
  moonbit_string_t _M0L6_2atmpS2827;
  struct moonbit_result_0 _tmp_4382;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2834;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3675;
  int32_t _M0L6_2acntS4058;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2833;
  moonbit_string_t _M0L6_2atmpS2831;
  moonbit_string_t _M0L6_2atmpS2832;
  struct moonbit_result_0 _tmp_4384;
  struct _M0TPB5ArrayGiE* _M0L1qS1163;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2847;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1164;
  int32_t* _M0L6_2atmpS2840;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS2839;
  int32_t _M0L6_2atmpS2837;
  void* _M0L4NoneS2838;
  struct moonbit_result_0 _tmp_4389;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2850;
  moonbit_string_t _M0L8_2afieldS3674;
  int32_t _M0L6_2acntS4064;
  moonbit_string_t _M0L2idS2848;
  moonbit_string_t _M0L6_2atmpS2849;
  struct moonbit_result_0 _tmp_4391;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2856;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3673;
  int32_t _M0L6_2acntS4070;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2855;
  moonbit_string_t _M0L6_2atmpS2853;
  moonbit_string_t _M0L6_2atmpS2854;
  struct moonbit_result_0 _tmp_4393;
  struct _M0TPB5ArrayGiE* _M0L1qS1169;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2865;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1170;
  int32_t* _M0L6_2atmpS2862;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS2861;
  int32_t _M0L6_2atmpS2859;
  void* _M0L4NoneS2860;
  struct moonbit_result_0 _result_4397;
  #line 58 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 59 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4378
  = _M0FP26biolab8bio__seq12parse__fastq((moonbit_string_t)moonbit_string_literal_116.data);
  if (_tmp_4378.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2866 =
      _tmp_4378.data.ok;
    _M0L7recordsS1161 = _M0L5_2aokS2866;
  } else {
    void* const _M0L6_2aerrS2867 = _tmp_4378.data.err;
    struct moonbit_result_0 _result_4379;
    _result_4379.tag = 0;
    _result_4379.data.err = _M0L6_2aerrS2867;
    return _result_4379;
  }
  #line 60 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2822
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161);
  _M0L6_2atmpS2823 = 0;
  #line 60 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4380
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS2822, 2, _M0L6_2atmpS2823, (moonbit_string_t)moonbit_string_literal_117.data);
  if (_M0L6_2atmpS2823) {
    moonbit_decref(_M0L6_2atmpS2823);
  }
  if (_tmp_4380.tag) {
    int32_t const _M0L5_2aokS2824 = _tmp_4380.data.ok;
  } else {
    void* const _M0L6_2aerrS2825 = _tmp_4380.data.err;
    struct moonbit_result_0 _result_4381;
    moonbit_decref(_M0L7recordsS1161);
    _result_4381.tag = 0;
    _result_4381.data.err = _M0L6_2aerrS2825;
    return _result_4381;
  }
  #line 61 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2828
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161, 0);
  _M0L8_2afieldS3676 = _M0L6_2atmpS2828->$1;
  _M0L6_2acntS4052
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2828));
  if (_M0L6_2acntS4052 > 1) {
    int32_t _M0L11_2anew__cntS4057 = _M0L6_2acntS4052 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2828), _M0L11_2anew__cntS4057);
    moonbit_incref(_M0L8_2afieldS3676);
  } else if (_M0L6_2acntS4052 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4056 =
      _M0L6_2atmpS2828->$4;
    moonbit_string_t _M0L8_2afieldS4055;
    moonbit_string_t _M0L8_2afieldS4054;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4053;
    moonbit_decref(_M0L8_2afieldS4056);
    _M0L8_2afieldS4055 = _M0L6_2atmpS2828->$3;
    moonbit_decref(_M0L8_2afieldS4055);
    _M0L8_2afieldS4054 = _M0L6_2atmpS2828->$2;
    moonbit_decref(_M0L8_2afieldS4054);
    _M0L8_2afieldS4053 = _M0L6_2atmpS2828->$0;
    moonbit_decref(_M0L8_2afieldS4053);
    #line 61 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2828);
  }
  _M0L2idS2826 = _M0L8_2afieldS3676;
  _M0L6_2atmpS2827 = 0;
  #line 61 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4382
  = _M0FPB10assert__eqGsE(_M0L2idS2826, (moonbit_string_t)moonbit_string_literal_11.data, _M0L6_2atmpS2827, (moonbit_string_t)moonbit_string_literal_118.data);
  moonbit_decref(_M0L2idS2826);
  if (_M0L6_2atmpS2827) {
    moonbit_decref(_M0L6_2atmpS2827);
  }
  if (_tmp_4382.tag) {
    int32_t const _M0L5_2aokS2829 = _tmp_4382.data.ok;
  } else {
    void* const _M0L6_2aerrS2830 = _tmp_4382.data.err;
    struct moonbit_result_0 _result_4383;
    moonbit_decref(_M0L7recordsS1161);
    _result_4383.tag = 0;
    _result_4383.data.err = _M0L6_2aerrS2830;
    return _result_4383;
  }
  #line 62 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2834
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161, 0);
  _M0L8_2afieldS3675 = _M0L6_2atmpS2834->$0;
  _M0L6_2acntS4058
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2834));
  if (_M0L6_2acntS4058 > 1) {
    int32_t _M0L11_2anew__cntS4063 = _M0L6_2acntS4058 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2834), _M0L11_2anew__cntS4063);
    moonbit_incref(_M0L8_2afieldS3675);
  } else if (_M0L6_2acntS4058 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4062 =
      _M0L6_2atmpS2834->$4;
    moonbit_string_t _M0L8_2afieldS4061;
    moonbit_string_t _M0L8_2afieldS4060;
    moonbit_string_t _M0L8_2afieldS4059;
    moonbit_decref(_M0L8_2afieldS4062);
    _M0L8_2afieldS4061 = _M0L6_2atmpS2834->$3;
    moonbit_decref(_M0L8_2afieldS4061);
    _M0L8_2afieldS4060 = _M0L6_2atmpS2834->$2;
    moonbit_decref(_M0L8_2afieldS4060);
    _M0L8_2afieldS4059 = _M0L6_2atmpS2834->$1;
    moonbit_decref(_M0L8_2afieldS4059);
    #line 62 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2834);
  }
  _M0L3seqS2833 = _M0L8_2afieldS3675;
  #line 62 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2831 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2833);
  moonbit_decref(_M0L3seqS2833);
  _M0L6_2atmpS2832 = 0;
  #line 62 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4384
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2831, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS2832, (moonbit_string_t)moonbit_string_literal_119.data);
  moonbit_decref(_M0L6_2atmpS2831);
  if (_M0L6_2atmpS2832) {
    moonbit_decref(_M0L6_2atmpS2832);
  }
  if (_tmp_4384.tag) {
    int32_t const _M0L5_2aokS2835 = _tmp_4384.data.ok;
  } else {
    void* const _M0L6_2aerrS2836 = _tmp_4384.data.err;
    struct moonbit_result_0 _result_4385;
    moonbit_decref(_M0L7recordsS1161);
    _result_4385.tag = 0;
    _result_4385.data.err = _M0L6_2aerrS2836;
    return _result_4385;
  }
  #line 63 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2847
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161, 0);
  #line 63 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1164
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L6_2atmpS2847);
  moonbit_decref(_M0L6_2atmpS2847);
  if (_M0L7_2abindS1164 == 0) {
    moonbit_string_t _M0L7_2abindS1167;
    int32_t _M0L6_2atmpS2844;
    struct _M0TPC16string10StringView _M0L6_2atmpS2843;
    struct moonbit_result_0 _tmp_4387;
    if (_M0L7_2abindS1164) {
      moonbit_decref(_M0L7_2abindS1164);
    }
    _M0L7_2abindS1167 = (moonbit_string_t)moonbit_string_literal_120.data;
    _M0L6_2atmpS2844 = Moonbit_array_length(_M0L7_2abindS1167);
    _M0L6_2atmpS2843
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1167, .$1 = 0, .$2 = _M0L6_2atmpS2844
    };
    #line 65 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _tmp_4387
    = _M0FPB4failGuE(_M0L6_2atmpS2843, (moonbit_string_t)moonbit_string_literal_121.data);
    moonbit_decref(_M0L6_2atmpS2843.$0);
    if (_tmp_4387.tag) {
      int32_t const _M0L5_2aokS2845 = _tmp_4387.data.ok;
    } else {
      void* const _M0L6_2aerrS2846 = _tmp_4387.data.err;
      struct moonbit_result_0 _result_4388;
      moonbit_decref(_M0L7recordsS1161);
      _result_4388.tag = 0;
      _result_4388.data.err = _M0L6_2aerrS2846;
      return _result_4388;
    }
  } else {
    struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1165 = _M0L7_2abindS1164;
    struct _M0TPB5ArrayGiE* _M0L4_2aqS1166 = _M0L7_2aSomeS1165;
    _M0L1qS1163 = _M0L4_2aqS1166;
    goto join_1162;
  }
  goto joinlet_4386;
  join_1162:;
  _M0L6_2atmpS2840 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS2840[0] = 40;
  _M0L6_2atmpS2840[1] = 40;
  _M0L6_2atmpS2840[2] = 40;
  _M0L6_2atmpS2840[3] = 40;
  _M0L6_2atmpS2839
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS2839)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS2839->$0 = _M0L6_2atmpS2840;
  _M0L6_2atmpS2839->$1 = 4;
  #line 64 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2837
  = _M0IPC15array5ArrayPB2Eq5equalGiE(_M0L1qS1163, _M0L6_2atmpS2839);
  moonbit_decref(_M0L1qS1163);
  moonbit_decref(_M0L6_2atmpS2839);
  _M0L4NoneS2838
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 64 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4389
  = _M0FPB12assert__true(_M0L6_2atmpS2837, _M0L4NoneS2838, (moonbit_string_t)moonbit_string_literal_122.data);
  moonbit_decref(_M0L4NoneS2838);
  if (_tmp_4389.tag) {
    int32_t const _M0L5_2aokS2841 = _tmp_4389.data.ok;
  } else {
    void* const _M0L6_2aerrS2842 = _tmp_4389.data.err;
    struct moonbit_result_0 _result_4390;
    moonbit_decref(_M0L7recordsS1161);
    _result_4390.tag = 0;
    _result_4390.data.err = _M0L6_2aerrS2842;
    return _result_4390;
  }
  joinlet_4386:;
  #line 67 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2850
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161, 1);
  _M0L8_2afieldS3674 = _M0L6_2atmpS2850->$1;
  _M0L6_2acntS4064
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2850));
  if (_M0L6_2acntS4064 > 1) {
    int32_t _M0L11_2anew__cntS4069 = _M0L6_2acntS4064 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2850), _M0L11_2anew__cntS4069);
    moonbit_incref(_M0L8_2afieldS3674);
  } else if (_M0L6_2acntS4064 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4068 =
      _M0L6_2atmpS2850->$4;
    moonbit_string_t _M0L8_2afieldS4067;
    moonbit_string_t _M0L8_2afieldS4066;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4065;
    moonbit_decref(_M0L8_2afieldS4068);
    _M0L8_2afieldS4067 = _M0L6_2atmpS2850->$3;
    moonbit_decref(_M0L8_2afieldS4067);
    _M0L8_2afieldS4066 = _M0L6_2atmpS2850->$2;
    moonbit_decref(_M0L8_2afieldS4066);
    _M0L8_2afieldS4065 = _M0L6_2atmpS2850->$0;
    moonbit_decref(_M0L8_2afieldS4065);
    #line 67 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2850);
  }
  _M0L2idS2848 = _M0L8_2afieldS3674;
  _M0L6_2atmpS2849 = 0;
  #line 67 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4391
  = _M0FPB10assert__eqGsE(_M0L2idS2848, (moonbit_string_t)moonbit_string_literal_123.data, _M0L6_2atmpS2849, (moonbit_string_t)moonbit_string_literal_124.data);
  moonbit_decref(_M0L2idS2848);
  if (_M0L6_2atmpS2849) {
    moonbit_decref(_M0L6_2atmpS2849);
  }
  if (_tmp_4391.tag) {
    int32_t const _M0L5_2aokS2851 = _tmp_4391.data.ok;
  } else {
    void* const _M0L6_2aerrS2852 = _tmp_4391.data.err;
    struct moonbit_result_0 _result_4392;
    moonbit_decref(_M0L7recordsS1161);
    _result_4392.tag = 0;
    _result_4392.data.err = _M0L6_2aerrS2852;
    return _result_4392;
  }
  #line 68 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2856
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161, 1);
  _M0L8_2afieldS3673 = _M0L6_2atmpS2856->$0;
  _M0L6_2acntS4070
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2856));
  if (_M0L6_2acntS4070 > 1) {
    int32_t _M0L11_2anew__cntS4075 = _M0L6_2acntS4070 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2856), _M0L11_2anew__cntS4075);
    moonbit_incref(_M0L8_2afieldS3673);
  } else if (_M0L6_2acntS4070 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4074 =
      _M0L6_2atmpS2856->$4;
    moonbit_string_t _M0L8_2afieldS4073;
    moonbit_string_t _M0L8_2afieldS4072;
    moonbit_string_t _M0L8_2afieldS4071;
    moonbit_decref(_M0L8_2afieldS4074);
    _M0L8_2afieldS4073 = _M0L6_2atmpS2856->$3;
    moonbit_decref(_M0L8_2afieldS4073);
    _M0L8_2afieldS4072 = _M0L6_2atmpS2856->$2;
    moonbit_decref(_M0L8_2afieldS4072);
    _M0L8_2afieldS4071 = _M0L6_2atmpS2856->$1;
    moonbit_decref(_M0L8_2afieldS4071);
    #line 68 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2856);
  }
  _M0L3seqS2855 = _M0L8_2afieldS3673;
  #line 68 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2853 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2855);
  moonbit_decref(_M0L3seqS2855);
  _M0L6_2atmpS2854 = 0;
  #line 68 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4393
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2853, (moonbit_string_t)moonbit_string_literal_29.data, _M0L6_2atmpS2854, (moonbit_string_t)moonbit_string_literal_125.data);
  moonbit_decref(_M0L6_2atmpS2853);
  if (_M0L6_2atmpS2854) {
    moonbit_decref(_M0L6_2atmpS2854);
  }
  if (_tmp_4393.tag) {
    int32_t const _M0L5_2aokS2857 = _tmp_4393.data.ok;
  } else {
    void* const _M0L6_2aerrS2858 = _tmp_4393.data.err;
    struct moonbit_result_0 _result_4394;
    moonbit_decref(_M0L7recordsS1161);
    _result_4394.tag = 0;
    _result_4394.data.err = _M0L6_2aerrS2858;
    return _result_4394;
  }
  #line 69 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2865
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1161, 1);
  moonbit_decref(_M0L7recordsS1161);
  #line 69 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1170
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L6_2atmpS2865);
  moonbit_decref(_M0L6_2atmpS2865);
  if (_M0L7_2abindS1170 == 0) {
    moonbit_string_t _M0L7_2abindS1173;
    int32_t _M0L6_2atmpS2864;
    struct _M0TPC16string10StringView _M0L6_2atmpS2863;
    struct moonbit_result_0 _result_4396;
    if (_M0L7_2abindS1170) {
      moonbit_decref(_M0L7_2abindS1170);
    }
    _M0L7_2abindS1173 = (moonbit_string_t)moonbit_string_literal_126.data;
    _M0L6_2atmpS2864 = Moonbit_array_length(_M0L7_2abindS1173);
    _M0L6_2atmpS2863
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1173, .$1 = 0, .$2 = _M0L6_2atmpS2864
    };
    #line 71 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _result_4396
    = _M0FPB4failGuE(_M0L6_2atmpS2863, (moonbit_string_t)moonbit_string_literal_127.data);
    moonbit_decref(_M0L6_2atmpS2863.$0);
    return _result_4396;
  } else {
    struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1171 = _M0L7_2abindS1170;
    struct _M0TPB5ArrayGiE* _M0L4_2aqS1172 = _M0L7_2aSomeS1171;
    _M0L1qS1169 = _M0L4_2aqS1172;
    goto join_1168;
  }
  join_1168:;
  _M0L6_2atmpS2862 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS2862[0] = 41;
  _M0L6_2atmpS2862[1] = 41;
  _M0L6_2atmpS2862[2] = 41;
  _M0L6_2atmpS2862[3] = 41;
  _M0L6_2atmpS2861
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS2861)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS2861->$0 = _M0L6_2atmpS2862;
  _M0L6_2atmpS2861->$1 = 4;
  #line 70 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2859
  = _M0IPC15array5ArrayPB2Eq5equalGiE(_M0L1qS1169, _M0L6_2atmpS2861);
  moonbit_decref(_M0L1qS1169);
  moonbit_decref(_M0L6_2atmpS2861);
  _M0L4NoneS2860
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 70 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4397
  = _M0FPB12assert__true(_M0L6_2atmpS2859, _M0L4NoneS2860, (moonbit_string_t)moonbit_string_literal_128.data);
  moonbit_decref(_M0L4NoneS2860);
  return _result_4397;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__4(
  
) {
  struct moonbit_result_1 _tmp_4398;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1154;
  int32_t _M0L6_2atmpS2793;
  moonbit_string_t _M0L6_2atmpS2794;
  struct moonbit_result_0 _tmp_4400;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2799;
  moonbit_string_t _M0L8_2afieldS3679;
  int32_t _M0L6_2acntS4076;
  moonbit_string_t _M0L2idS2797;
  moonbit_string_t _M0L6_2atmpS2798;
  struct moonbit_result_0 _tmp_4402;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2804;
  moonbit_string_t _M0L8_2afieldS3678;
  int32_t _M0L6_2acntS4082;
  moonbit_string_t _M0L11descriptionS2802;
  moonbit_string_t _M0L6_2atmpS2803;
  struct moonbit_result_0 _tmp_4404;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2810;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3677;
  int32_t _M0L6_2acntS4088;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2809;
  moonbit_string_t _M0L6_2atmpS2807;
  moonbit_string_t _M0L6_2atmpS2808;
  struct moonbit_result_0 _tmp_4406;
  struct _M0TPB5ArrayGiE* _M0L1qS1156;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2819;
  struct _M0TPB5ArrayGiE* _M0L7_2abindS1157;
  int32_t* _M0L6_2atmpS2816;
  struct _M0TPB5ArrayGiE* _M0L6_2atmpS2815;
  int32_t _M0L6_2atmpS2813;
  void* _M0L4NoneS2814;
  struct moonbit_result_0 _result_4410;
  #line 45 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 46 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4398
  = _M0FP26biolab8bio__seq12parse__fastq((moonbit_string_t)moonbit_string_literal_9.data);
  if (_tmp_4398.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* const _M0L5_2aokS2820 =
      _tmp_4398.data.ok;
    _M0L7recordsS1154 = _M0L5_2aokS2820;
  } else {
    void* const _M0L6_2aerrS2821 = _tmp_4398.data.err;
    struct moonbit_result_0 _result_4399;
    _result_4399.tag = 0;
    _result_4399.data.err = _M0L6_2aerrS2821;
    return _result_4399;
  }
  #line 47 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2793
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1154);
  _M0L6_2atmpS2794 = 0;
  #line 47 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4400
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS2793, 1, _M0L6_2atmpS2794, (moonbit_string_t)moonbit_string_literal_129.data);
  if (_M0L6_2atmpS2794) {
    moonbit_decref(_M0L6_2atmpS2794);
  }
  if (_tmp_4400.tag) {
    int32_t const _M0L5_2aokS2795 = _tmp_4400.data.ok;
  } else {
    void* const _M0L6_2aerrS2796 = _tmp_4400.data.err;
    struct moonbit_result_0 _result_4401;
    moonbit_decref(_M0L7recordsS1154);
    _result_4401.tag = 0;
    _result_4401.data.err = _M0L6_2aerrS2796;
    return _result_4401;
  }
  #line 48 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2799
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1154, 0);
  _M0L8_2afieldS3679 = _M0L6_2atmpS2799->$1;
  _M0L6_2acntS4076
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2799));
  if (_M0L6_2acntS4076 > 1) {
    int32_t _M0L11_2anew__cntS4081 = _M0L6_2acntS4076 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2799), _M0L11_2anew__cntS4081);
    moonbit_incref(_M0L8_2afieldS3679);
  } else if (_M0L6_2acntS4076 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4080 =
      _M0L6_2atmpS2799->$4;
    moonbit_string_t _M0L8_2afieldS4079;
    moonbit_string_t _M0L8_2afieldS4078;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4077;
    moonbit_decref(_M0L8_2afieldS4080);
    _M0L8_2afieldS4079 = _M0L6_2atmpS2799->$3;
    moonbit_decref(_M0L8_2afieldS4079);
    _M0L8_2afieldS4078 = _M0L6_2atmpS2799->$2;
    moonbit_decref(_M0L8_2afieldS4078);
    _M0L8_2afieldS4077 = _M0L6_2atmpS2799->$0;
    moonbit_decref(_M0L8_2afieldS4077);
    #line 48 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2799);
  }
  _M0L2idS2797 = _M0L8_2afieldS3679;
  _M0L6_2atmpS2798 = 0;
  #line 48 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4402
  = _M0FPB10assert__eqGsE(_M0L2idS2797, (moonbit_string_t)moonbit_string_literal_11.data, _M0L6_2atmpS2798, (moonbit_string_t)moonbit_string_literal_130.data);
  moonbit_decref(_M0L2idS2797);
  if (_M0L6_2atmpS2798) {
    moonbit_decref(_M0L6_2atmpS2798);
  }
  if (_tmp_4402.tag) {
    int32_t const _M0L5_2aokS2800 = _tmp_4402.data.ok;
  } else {
    void* const _M0L6_2aerrS2801 = _tmp_4402.data.err;
    struct moonbit_result_0 _result_4403;
    moonbit_decref(_M0L7recordsS1154);
    _result_4403.tag = 0;
    _result_4403.data.err = _M0L6_2aerrS2801;
    return _result_4403;
  }
  #line 49 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2804
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1154, 0);
  _M0L8_2afieldS3678 = _M0L6_2atmpS2804->$3;
  _M0L6_2acntS4082
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2804));
  if (_M0L6_2acntS4082 > 1) {
    int32_t _M0L11_2anew__cntS4087 = _M0L6_2acntS4082 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2804), _M0L11_2anew__cntS4087);
    moonbit_incref(_M0L8_2afieldS3678);
  } else if (_M0L6_2acntS4082 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4086 =
      _M0L6_2atmpS2804->$4;
    moonbit_string_t _M0L8_2afieldS4085;
    moonbit_string_t _M0L8_2afieldS4084;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4083;
    moonbit_decref(_M0L8_2afieldS4086);
    _M0L8_2afieldS4085 = _M0L6_2atmpS2804->$2;
    moonbit_decref(_M0L8_2afieldS4085);
    _M0L8_2afieldS4084 = _M0L6_2atmpS2804->$1;
    moonbit_decref(_M0L8_2afieldS4084);
    _M0L8_2afieldS4083 = _M0L6_2atmpS2804->$0;
    moonbit_decref(_M0L8_2afieldS4083);
    #line 49 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2804);
  }
  _M0L11descriptionS2802 = _M0L8_2afieldS3678;
  _M0L6_2atmpS2803 = 0;
  #line 49 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4404
  = _M0FPB10assert__eqGsE(_M0L11descriptionS2802, (moonbit_string_t)moonbit_string_literal_13.data, _M0L6_2atmpS2803, (moonbit_string_t)moonbit_string_literal_131.data);
  moonbit_decref(_M0L11descriptionS2802);
  if (_M0L6_2atmpS2803) {
    moonbit_decref(_M0L6_2atmpS2803);
  }
  if (_tmp_4404.tag) {
    int32_t const _M0L5_2aokS2805 = _tmp_4404.data.ok;
  } else {
    void* const _M0L6_2aerrS2806 = _tmp_4404.data.err;
    struct moonbit_result_0 _result_4405;
    moonbit_decref(_M0L7recordsS1154);
    _result_4405.tag = 0;
    _result_4405.data.err = _M0L6_2aerrS2806;
    return _result_4405;
  }
  #line 50 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2810
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1154, 0);
  _M0L8_2afieldS3677 = _M0L6_2atmpS2810->$0;
  _M0L6_2acntS4088
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2810));
  if (_M0L6_2acntS4088 > 1) {
    int32_t _M0L11_2anew__cntS4093 = _M0L6_2acntS4088 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2810), _M0L11_2anew__cntS4093);
    moonbit_incref(_M0L8_2afieldS3677);
  } else if (_M0L6_2acntS4088 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4092 =
      _M0L6_2atmpS2810->$4;
    moonbit_string_t _M0L8_2afieldS4091;
    moonbit_string_t _M0L8_2afieldS4090;
    moonbit_string_t _M0L8_2afieldS4089;
    moonbit_decref(_M0L8_2afieldS4092);
    _M0L8_2afieldS4091 = _M0L6_2atmpS2810->$3;
    moonbit_decref(_M0L8_2afieldS4091);
    _M0L8_2afieldS4090 = _M0L6_2atmpS2810->$2;
    moonbit_decref(_M0L8_2afieldS4090);
    _M0L8_2afieldS4089 = _M0L6_2atmpS2810->$1;
    moonbit_decref(_M0L8_2afieldS4089);
    #line 50 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2810);
  }
  _M0L3seqS2809 = _M0L8_2afieldS3677;
  #line 50 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2807 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2809);
  moonbit_decref(_M0L3seqS2809);
  _M0L6_2atmpS2808 = 0;
  #line 50 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4406
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2807, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS2808, (moonbit_string_t)moonbit_string_literal_132.data);
  moonbit_decref(_M0L6_2atmpS2807);
  if (_M0L6_2atmpS2808) {
    moonbit_decref(_M0L6_2atmpS2808);
  }
  if (_tmp_4406.tag) {
    int32_t const _M0L5_2aokS2811 = _tmp_4406.data.ok;
  } else {
    void* const _M0L6_2aerrS2812 = _tmp_4406.data.err;
    struct moonbit_result_0 _result_4407;
    moonbit_decref(_M0L7recordsS1154);
    _result_4407.tag = 0;
    _result_4407.data.err = _M0L6_2aerrS2812;
    return _result_4407;
  }
  #line 51 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2819
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1154, 0);
  moonbit_decref(_M0L7recordsS1154);
  #line 51 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7_2abindS1157
  = _M0MP26biolab8bio__seq9SeqRecord14phred__quality(_M0L6_2atmpS2819);
  moonbit_decref(_M0L6_2atmpS2819);
  if (_M0L7_2abindS1157 == 0) {
    moonbit_string_t _M0L7_2abindS1160;
    int32_t _M0L6_2atmpS2818;
    struct _M0TPC16string10StringView _M0L6_2atmpS2817;
    struct moonbit_result_0 _result_4409;
    if (_M0L7_2abindS1157) {
      moonbit_decref(_M0L7_2abindS1157);
    }
    _M0L7_2abindS1160 = (moonbit_string_t)moonbit_string_literal_133.data;
    _M0L6_2atmpS2818 = Moonbit_array_length(_M0L7_2abindS1160);
    _M0L6_2atmpS2817
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L7_2abindS1160, .$1 = 0, .$2 = _M0L6_2atmpS2818
    };
    #line 53 "/home/developer/Documents/1/seqio_wbtest.mbt"
    _result_4409
    = _M0FPB4failGuE(_M0L6_2atmpS2817, (moonbit_string_t)moonbit_string_literal_134.data);
    moonbit_decref(_M0L6_2atmpS2817.$0);
    return _result_4409;
  } else {
    struct _M0TPB5ArrayGiE* _M0L7_2aSomeS1158 = _M0L7_2abindS1157;
    struct _M0TPB5ArrayGiE* _M0L4_2aqS1159 = _M0L7_2aSomeS1158;
    _M0L1qS1156 = _M0L4_2aqS1159;
    goto join_1155;
  }
  join_1155:;
  _M0L6_2atmpS2816 = (int32_t*)moonbit_make_int32_array_raw(4);
  _M0L6_2atmpS2816[0] = 40;
  _M0L6_2atmpS2816[1] = 40;
  _M0L6_2atmpS2816[2] = 40;
  _M0L6_2atmpS2816[3] = 40;
  _M0L6_2atmpS2815
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L6_2atmpS2815)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _M0L6_2atmpS2815->$0 = _M0L6_2atmpS2816;
  _M0L6_2atmpS2815->$1 = 4;
  #line 52 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2813
  = _M0IPC15array5ArrayPB2Eq5equalGiE(_M0L1qS1156, _M0L6_2atmpS2815);
  moonbit_decref(_M0L1qS1156);
  moonbit_decref(_M0L6_2atmpS2815);
  _M0L4NoneS2814
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 52 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4410
  = _M0FPB12assert__true(_M0L6_2atmpS2813, _M0L4NoneS2814, (moonbit_string_t)moonbit_string_literal_135.data);
  moonbit_decref(_M0L4NoneS2814);
  return _result_4410;
}

struct _M0TPB5ArrayGiE* _M0MP26biolab8bio__seq9SeqRecord14phred__quality(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS1153
) {
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L19letter__annotationsS2792;
  struct _M0TPB5ArrayGiE* _result_4411;
  #line 78 "/home/developer/Documents/1/seq_record.mbt"
  _M0L19letter__annotationsS2792 = _M0L4selfS1153->$4;
  moonbit_incref(_M0L19letter__annotationsS2792);
  #line 79 "/home/developer/Documents/1/seq_record.mbt"
  _result_4411
  = _M0MPB3Map3getGsRPB5ArrayGiEE(_M0L19letter__annotationsS2792, (moonbit_string_t)moonbit_string_literal_95.data);
  moonbit_decref(_M0L19letter__annotationsS2792);
  return _result_4411;
}

struct moonbit_result_1 _M0FP26biolab8bio__seq12parse__fastq(
  moonbit_string_t _M0L7contentS1112
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1110;
  moonbit_string_t _M0L7_2abindS1113;
  int32_t _M0L6_2atmpS2791;
  struct _M0TPC16string10StringView _M0L6_2atmpS2790;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2789;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS1111;
  int32_t _M0L1nS1114;
  struct _M0TPB8MutLocalGiE* _M0L1iS1115;
  int32_t _M0L6_2atmpS2788;
  int32_t _M0L2atS1116;
  int32_t _M0L6_2atmpS2787;
  int32_t _M0L4plusS1117;
  int32_t _M0L6_2atmpS2786;
  int32_t _M0L2crS1118;
  struct moonbit_result_1 _result_4427;
  #line 22 "/home/developer/Documents/1/fastq_io.mbt"
  #line 23 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L7recordsS1110
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS1113 = (moonbit_string_t)moonbit_string_literal_89.data;
  _M0L6_2atmpS2791 = Moonbit_array_length(_M0L7_2abindS1113);
  _M0L6_2atmpS2790
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1113, .$1 = 0, .$2 = _M0L6_2atmpS2791
  };
  #line 24 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L6_2atmpS2789
  = _M0MPC16string6String5split(_M0L7contentS1112, _M0L6_2atmpS2790);
  moonbit_decref(_M0L6_2atmpS2790.$0);
  #line 24 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L5linesS1111
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2789);
  moonbit_decref(_M0L6_2atmpS2789);
  #line 25 "/home/developer/Documents/1/fastq_io.mbt"
  _M0L1nS1114
  = _M0MPC15array5Array6lengthGRPC16string10StringViewE(_M0L5linesS1111);
  _M0L1iS1115
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1115)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1115->$0 = 0;
  _M0L6_2atmpS2788 = 64;
  _M0L2atS1116 = (uint16_t)_M0L6_2atmpS2788;
  _M0L6_2atmpS2787 = 43;
  _M0L4plusS1117 = (uint16_t)_M0L6_2atmpS2787;
  _M0L6_2atmpS2786 = 13;
  _M0L2crS1118 = (uint16_t)_M0L6_2atmpS2786;
  while (1) {
    int32_t _M0L3valS2716 = _M0L1iS1115->$0;
    int32_t _M0L6_2atmpS2715 = _M0L3valS2716 + 3;
    if (_M0L6_2atmpS2715 < _M0L1nS1114) {
      int32_t _M0L3valS2785 = _M0L1iS1115->$0;
      struct _M0TPC16string10StringView _M0L5line0S1119;
      int32_t _M0L3valS2784;
      int32_t _M0L6_2atmpS2783;
      struct _M0TPC16string10StringView _M0L5line1S1120;
      int32_t _M0L3valS2782;
      int32_t _M0L6_2atmpS2781;
      struct _M0TPC16string10StringView _M0L5line2S1121;
      int32_t _M0L3valS2780;
      int32_t _M0L6_2atmpS2779;
      struct _M0TPC16string10StringView _M0L5line3S1122;
      int32_t _M0L4len0S1123;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed0S1124;
      int32_t _if__result_4413;
      int32_t _M0L3valS2722;
      int32_t _if__result_4414;
      int32_t _M0L3valS2724;
      int32_t _if__result_4415;
      int32_t _M0L4len2S1127;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed2S1128;
      int32_t _if__result_4417;
      int32_t _M0L3valS2731;
      int32_t _if__result_4418;
      int32_t _M0L4len1S1130;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed1S1131;
      int32_t _if__result_4420;
      int32_t _M0L4len3S1132;
      struct _M0TPB8MutLocalGiE* _M0L8trimmed3S1133;
      int32_t _if__result_4421;
      int32_t _M0L3valS2740;
      int32_t _M0L3valS2741;
      int32_t _M0L3valS2775;
      moonbit_string_t _M0L5titleS1134;
      moonbit_string_t _M0L2idS1136;
      moonbit_string_t _M0L11descriptionS1137;
      int32_t _M0L3posS1148;
      moonbit_string_t _M0L7_2abindS1150;
      int32_t _M0L6_2atmpS2774;
      struct _M0TPC16string10StringView _M0L6_2atmpS2773;
      int64_t _M0L7_2abindS1149;
      int64_t _M0L6_2atmpS2772;
      struct _M0TPC16string10StringView _M0L6_2atmpS2771;
      moonbit_string_t _M0L6_2atmpS2770;
      int32_t _M0L3valS2766;
      moonbit_string_t _M0L8seq__strS1138;
      int32_t _M0L3valS2765;
      int32_t* _M0L5qualsS1139;
      int32_t _M0L1jS1140;
      struct _M0TPB5ArrayGiE* _M0L10quals__arrS1142;
      int32_t _M0L1jS1143;
      struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1146;
      struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2764;
      struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2763;
      struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L11letter__annS1145;
      struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2760;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2759;
      int32_t _M0L3valS2762;
      int32_t _M0L6_2atmpS2761;
      #line 32 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line0S1119
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1111, _M0L3valS2785);
      _M0L3valS2784 = _M0L1iS1115->$0;
      _M0L6_2atmpS2783 = _M0L3valS2784 + 1;
      #line 33 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line1S1120
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1111, _M0L6_2atmpS2783);
      _M0L3valS2782 = _M0L1iS1115->$0;
      _M0L6_2atmpS2781 = _M0L3valS2782 + 2;
      #line 34 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line2S1121
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1111, _M0L6_2atmpS2781);
      _M0L3valS2780 = _M0L1iS1115->$0;
      _M0L6_2atmpS2779 = _M0L3valS2780 + 3;
      #line 35 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L5line3S1122
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1111, _M0L6_2atmpS2779);
      #line 37 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len0S1123 = _M0MPC16string10StringView6length(_M0L5line0S1119);
      _M0L8trimmed0S1124
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed0S1124)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed0S1124->$0 = _M0L4len0S1123;
      if (_M0L4len0S1123 > 0) {
        int32_t _M0L6_2atmpS2718 = _M0L4len0S1123 - 1;
        int32_t _M0L6_2atmpS2717;
        #line 39 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2717
        = _M0MPC16string10StringView11unsafe__get(_M0L5line0S1119, _M0L6_2atmpS2718);
        _if__result_4413 = _M0L6_2atmpS2717 == _M0L2crS1118;
      } else {
        _if__result_4413 = 0;
      }
      if (_if__result_4413) {
        int32_t _M0L6_2atmpS2719 = _M0L4len0S1123 - 1;
        _M0L8trimmed0S1124->$0 = _M0L6_2atmpS2719;
      }
      _M0L3valS2722 = _M0L8trimmed0S1124->$0;
      if (_M0L3valS2722 == 0) {
        int32_t _M0L3valS2721 = _M0L1iS1115->$0;
        int32_t _M0L6_2atmpS2720 = _M0L3valS2721 + 4;
        _if__result_4414 = _M0L6_2atmpS2720 >= _M0L1nS1114;
      } else {
        _if__result_4414 = 0;
      }
      if (_if__result_4414) {
        moonbit_decref(_M0L8trimmed0S1124);
        moonbit_decref(_M0L5line3S1122.$0);
        moonbit_decref(_M0L5line2S1121.$0);
        moonbit_decref(_M0L5line1S1120.$0);
        moonbit_decref(_M0L5line0S1119.$0);
        moonbit_decref(_M0L1iS1115);
        moonbit_decref(_M0L5linesS1111);
        break;
      }
      _M0L3valS2724 = _M0L8trimmed0S1124->$0;
      if (_M0L3valS2724 == 0) {
        _if__result_4415 = 1;
      } else {
        int32_t _M0L6_2atmpS2723;
        #line 45 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2723
        = _M0MPC16string10StringView11unsafe__get(_M0L5line0S1119, 0);
        _if__result_4415 = _M0L6_2atmpS2723 != _M0L2atS1116;
      }
      if (_if__result_4415) {
        moonbit_string_t _M0L10line0__strS1126;
        moonbit_string_t _M0L6_2atmpS2726;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2725;
        struct moonbit_result_1 _result_4416;
        moonbit_decref(_M0L8trimmed0S1124);
        moonbit_decref(_M0L5line3S1122.$0);
        moonbit_decref(_M0L5line2S1121.$0);
        moonbit_decref(_M0L5line1S1120.$0);
        moonbit_decref(_M0L1iS1115);
        moonbit_decref(_M0L5linesS1111);
        moonbit_decref(_M0L7recordsS1110);
        #line 46 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L10line0__strS1126
        = _M0MPC16string10StringView9to__owned(_M0L5line0S1119);
        moonbit_decref(_M0L5line0S1119.$0);
        #line 48 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2726
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_136.data, _M0L10line0__strS1126);
        moonbit_decref(_M0L10line0__strS1126);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2725
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2725)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2725)->$0
        = _M0L6_2atmpS2726;
        _result_4416.tag = 0;
        _result_4416.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2725;
        return _result_4416;
      }
      #line 52 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len2S1127 = _M0MPC16string10StringView6length(_M0L5line2S1121);
      _M0L8trimmed2S1128
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed2S1128)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed2S1128->$0 = _M0L4len2S1127;
      if (_M0L4len2S1127 > 0) {
        int32_t _M0L6_2atmpS2728 = _M0L4len2S1127 - 1;
        int32_t _M0L6_2atmpS2727;
        #line 54 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2727
        = _M0MPC16string10StringView11unsafe__get(_M0L5line2S1121, _M0L6_2atmpS2728);
        _if__result_4417 = _M0L6_2atmpS2727 == _M0L2crS1118;
      } else {
        _if__result_4417 = 0;
      }
      if (_if__result_4417) {
        int32_t _M0L6_2atmpS2729 = _M0L4len2S1127 - 1;
        _M0L8trimmed2S1128->$0 = _M0L6_2atmpS2729;
      }
      _M0L3valS2731 = _M0L8trimmed2S1128->$0;
      moonbit_decref(_M0L8trimmed2S1128);
      if (_M0L3valS2731 == 0) {
        _if__result_4418 = 1;
      } else {
        int32_t _M0L6_2atmpS2730;
        #line 57 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2730
        = _M0MPC16string10StringView11unsafe__get(_M0L5line2S1121, 0);
        _if__result_4418 = _M0L6_2atmpS2730 != _M0L4plusS1117;
      }
      if (_if__result_4418) {
        moonbit_string_t _M0L10line2__strS1129;
        moonbit_string_t _M0L6_2atmpS2733;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2732;
        struct moonbit_result_1 _result_4419;
        moonbit_decref(_M0L8trimmed0S1124);
        moonbit_decref(_M0L5line3S1122.$0);
        moonbit_decref(_M0L5line1S1120.$0);
        moonbit_decref(_M0L5line0S1119.$0);
        moonbit_decref(_M0L1iS1115);
        moonbit_decref(_M0L5linesS1111);
        moonbit_decref(_M0L7recordsS1110);
        #line 58 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L10line2__strS1129
        = _M0MPC16string10StringView9to__owned(_M0L5line2S1121);
        moonbit_decref(_M0L5line2S1121.$0);
        #line 60 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2733
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_137.data, _M0L10line2__strS1129);
        moonbit_decref(_M0L10line2__strS1129);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2732
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2732)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2732)->$0
        = _M0L6_2atmpS2733;
        _result_4419.tag = 0;
        _result_4419.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2732;
        return _result_4419;
      } else {
        moonbit_decref(_M0L5line2S1121.$0);
      }
      #line 64 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len1S1130 = _M0MPC16string10StringView6length(_M0L5line1S1120);
      _M0L8trimmed1S1131
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed1S1131)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed1S1131->$0 = _M0L4len1S1130;
      if (_M0L4len1S1130 > 0) {
        int32_t _M0L6_2atmpS2735 = _M0L4len1S1130 - 1;
        int32_t _M0L6_2atmpS2734;
        #line 66 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2734
        = _M0MPC16string10StringView11unsafe__get(_M0L5line1S1120, _M0L6_2atmpS2735);
        _if__result_4420 = _M0L6_2atmpS2734 == _M0L2crS1118;
      } else {
        _if__result_4420 = 0;
      }
      if (_if__result_4420) {
        int32_t _M0L6_2atmpS2736 = _M0L4len1S1130 - 1;
        _M0L8trimmed1S1131->$0 = _M0L6_2atmpS2736;
      }
      #line 70 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L4len3S1132 = _M0MPC16string10StringView6length(_M0L5line3S1122);
      _M0L8trimmed3S1133
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L8trimmed3S1133)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L8trimmed3S1133->$0 = _M0L4len3S1132;
      if (_M0L4len3S1132 > 0) {
        int32_t _M0L6_2atmpS2738 = _M0L4len3S1132 - 1;
        int32_t _M0L6_2atmpS2737;
        #line 72 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2737
        = _M0MPC16string10StringView11unsafe__get(_M0L5line3S1122, _M0L6_2atmpS2738);
        _if__result_4421 = _M0L6_2atmpS2737 == _M0L2crS1118;
      } else {
        _if__result_4421 = 0;
      }
      if (_if__result_4421) {
        int32_t _M0L6_2atmpS2739 = _M0L4len3S1132 - 1;
        _M0L8trimmed3S1133->$0 = _M0L6_2atmpS2739;
      }
      _M0L3valS2740 = _M0L8trimmed1S1131->$0;
      _M0L3valS2741 = _M0L8trimmed3S1133->$0;
      if (_M0L3valS2740 != _M0L3valS2741) {
        int32_t _M0L3valS2750;
        moonbit_string_t _M0L6_2atmpS2749;
        moonbit_string_t _M0L6_2atmpS2748;
        moonbit_string_t _M0L6_2atmpS2745;
        int32_t _M0L3valS2747;
        moonbit_string_t _M0L6_2atmpS2746;
        moonbit_string_t _M0L6_2atmpS2744;
        moonbit_string_t _M0L6_2atmpS2743;
        void* _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2742;
        struct moonbit_result_1 _result_4422;
        moonbit_decref(_M0L8trimmed0S1124);
        moonbit_decref(_M0L5line3S1122.$0);
        moonbit_decref(_M0L5line1S1120.$0);
        moonbit_decref(_M0L5line0S1119.$0);
        moonbit_decref(_M0L1iS1115);
        moonbit_decref(_M0L5linesS1111);
        moonbit_decref(_M0L7recordsS1110);
        _M0L3valS2750 = _M0L8trimmed1S1131->$0;
        moonbit_decref(_M0L8trimmed1S1131);
        #line 79 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2749
        = _M0MPC13int3Int18to__string_2einner(_M0L3valS2750, 10);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2748
        = moonbit_add_string((moonbit_string_t)moonbit_string_literal_97.data, _M0L6_2atmpS2749);
        moonbit_decref(_M0L6_2atmpS2749);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2745
        = moonbit_add_string(_M0L6_2atmpS2748, (moonbit_string_t)moonbit_string_literal_98.data);
        moonbit_decref(_M0L6_2atmpS2748);
        _M0L3valS2747 = _M0L8trimmed3S1133->$0;
        moonbit_decref(_M0L8trimmed3S1133);
        #line 81 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2746
        = _M0MPC13int3Int18to__string_2einner(_M0L3valS2747, 10);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2744
        = moonbit_add_string(_M0L6_2atmpS2745, _M0L6_2atmpS2746);
        moonbit_decref(_M0L6_2atmpS2746);
        moonbit_decref(_M0L6_2atmpS2745);
        #line 78 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2743
        = moonbit_add_string(_M0L6_2atmpS2744, (moonbit_string_t)moonbit_string_literal_138.data);
        moonbit_decref(_M0L6_2atmpS2744);
        _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2742
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError));
        Moonbit_object_header(_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2742)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 37, 1);
        ((struct _M0DTPC15error5Error43biolab_2fbio__seq_2eSeqIOError_2eSeqIOError*)_M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2742)->$0
        = _M0L6_2atmpS2743;
        _result_4422.tag = 0;
        _result_4422.data.err
        = _M0L43biolab_2fbio__seq_2eSeqIOError_2eSeqIOErrorS2742;
        return _result_4422;
      }
      _M0L3valS2775 = _M0L8trimmed0S1124->$0;
      if (_M0L3valS2775 > 1) {
        int32_t _M0L3valS2778 = _M0L8trimmed0S1124->$0;
        int64_t _M0L6_2atmpS2777;
        struct _M0TPC16string10StringView _M0L6_2atmpS2776;
        moonbit_decref(_M0L8trimmed0S1124);
        _M0L6_2atmpS2777 = (int64_t)_M0L3valS2778;
        #line 86 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2776
        = _M0MPC16string10StringView11sub_2einner(_M0L5line0S1119, 1, _M0L6_2atmpS2777);
        moonbit_decref(_M0L5line0S1119.$0);
        #line 86 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L5titleS1134
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2776);
        moonbit_decref(_M0L6_2atmpS2776.$0);
      } else {
        moonbit_decref(_M0L8trimmed0S1124);
        moonbit_decref(_M0L5line0S1119.$0);
        _M0L5titleS1134 = (moonbit_string_t)moonbit_string_literal_0.data;
      }
      _M0L7_2abindS1150 = (moonbit_string_t)moonbit_string_literal_115.data;
      _M0L6_2atmpS2774 = Moonbit_array_length(_M0L7_2abindS1150);
      _M0L6_2atmpS2773
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L7_2abindS1150, .$1 = 0, .$2 = _M0L6_2atmpS2774
      };
      #line 87 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L7_2abindS1149
      = _M0MPC16string6String4find(_M0L5titleS1134, _M0L6_2atmpS2773);
      moonbit_decref(_M0L6_2atmpS2773.$0);
      if (_M0L7_2abindS1149 == 4294967296ll) {
        moonbit_incref(_M0L5titleS1134);
        _M0L2idS1136 = _M0L5titleS1134;
        _M0L11descriptionS1137 = _M0L5titleS1134;
        goto join_1135;
      } else {
        int64_t _M0L7_2aSomeS1151 = _M0L7_2abindS1149;
        int32_t _M0L6_2aposS1152 = (int32_t)_M0L7_2aSomeS1151;
        _M0L3posS1148 = _M0L6_2aposS1152;
        goto join_1147;
      }
      goto joinlet_4424;
      join_1147:;
      _M0L6_2atmpS2772 = (int64_t)_M0L3posS1148;
      #line 88 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2771
      = _M0MPC16string6String11sub_2einner(_M0L5titleS1134, 0, _M0L6_2atmpS2772);
      #line 88 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2770
      = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2771);
      moonbit_decref(_M0L6_2atmpS2771.$0);
      _M0L2idS1136 = _M0L6_2atmpS2770;
      _M0L11descriptionS1137 = _M0L5titleS1134;
      goto join_1135;
      joinlet_4424:;
      goto joinlet_4423;
      join_1135:;
      _M0L3valS2766 = _M0L8trimmed1S1131->$0;
      if (_M0L3valS2766 > 0) {
        int32_t _M0L3valS2769 = _M0L8trimmed1S1131->$0;
        int64_t _M0L6_2atmpS2768;
        struct _M0TPC16string10StringView _M0L6_2atmpS2767;
        moonbit_decref(_M0L8trimmed1S1131);
        _M0L6_2atmpS2768 = (int64_t)_M0L3valS2769;
        #line 92 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L6_2atmpS2767
        = _M0MPC16string10StringView11sub_2einner(_M0L5line1S1120, 0, _M0L6_2atmpS2768);
        moonbit_decref(_M0L5line1S1120.$0);
        #line 92 "/home/developer/Documents/1/fastq_io.mbt"
        _M0L8seq__strS1138
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2767);
        moonbit_decref(_M0L6_2atmpS2767.$0);
      } else {
        moonbit_decref(_M0L8trimmed1S1131);
        moonbit_decref(_M0L5line1S1120.$0);
        _M0L8seq__strS1138 = (moonbit_string_t)moonbit_string_literal_0.data;
      }
      _M0L3valS2765 = _M0L8trimmed3S1133->$0;
      _M0L5qualsS1139 = (int32_t*)moonbit_make_int32_array(_M0L3valS2765, 0);
      _M0L1jS1140 = 0;
      while (1) {
        int32_t _M0L3valS2751 = _M0L8trimmed3S1133->$0;
        if (_M0L1jS1140 < _M0L3valS2751) {
          int32_t _M0L6_2atmpS2754;
          int32_t _M0L6_2atmpS2753;
          int32_t _M0L6_2atmpS2752;
          int32_t _M0L6_2atmpS2755;
          #line 96 "/home/developer/Documents/1/fastq_io.mbt"
          _M0L6_2atmpS2754
          = _M0MPC16string10StringView11unsafe__get(_M0L5line3S1122, _M0L1jS1140);
          _M0L6_2atmpS2753 = (int32_t)_M0L6_2atmpS2754;
          _M0L6_2atmpS2752 = _M0L6_2atmpS2753 - 33;
          if (
            _M0L1jS1140 < 0
            || _M0L1jS1140 >= Moonbit_array_length(_M0L5qualsS1139)
          ) {
            #line 96 "/home/developer/Documents/1/fastq_io.mbt"
            moonbit_panic();
          }
          _M0L5qualsS1139[_M0L1jS1140] = _M0L6_2atmpS2752;
          _M0L6_2atmpS2755 = _M0L1jS1140 + 1;
          _M0L1jS1140 = _M0L6_2atmpS2755;
          continue;
        } else {
          moonbit_decref(_M0L5line3S1122.$0);
        }
        break;
      }
      #line 98 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L10quals__arrS1142 = _M0MPC15array5Array11new_2einnerGiE(0);
      _M0L1jS1143 = 0;
      while (1) {
        int32_t _M0L3valS2756 = _M0L8trimmed3S1133->$0;
        if (_M0L1jS1143 < _M0L3valS2756) {
          int32_t _M0L6_2atmpS2757;
          int32_t _M0L6_2atmpS2758;
          if (
            _M0L1jS1143 < 0
            || _M0L1jS1143 >= Moonbit_array_length(_M0L5qualsS1139)
          ) {
            #line 100 "/home/developer/Documents/1/fastq_io.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2757 = (int32_t)_M0L5qualsS1139[_M0L1jS1143];
          #line 100 "/home/developer/Documents/1/fastq_io.mbt"
          _M0MPC15array5Array4pushGiE(_M0L10quals__arrS1142, _M0L6_2atmpS2757);
          _M0L6_2atmpS2758 = _M0L1jS1143 + 1;
          _M0L1jS1143 = _M0L6_2atmpS2758;
          continue;
        } else {
          moonbit_decref(_M0L5qualsS1139);
          moonbit_decref(_M0L8trimmed3S1133);
        }
        break;
      }
      _M0L7_2abindS1146
      = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
      _M0L6_2atmpS2764 = _M0L7_2abindS1146;
      _M0L6_2atmpS2763
      = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
        .$0 = _M0L6_2atmpS2764, .$1 = 0, .$2 = 0
      };
      #line 103 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L11letter__annS1145
      = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2763, 1ll);
      moonbit_decref(_M0L6_2atmpS2763.$0);
      #line 104 "/home/developer/Documents/1/fastq_io.mbt"
      _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L11letter__annS1145, (moonbit_string_t)moonbit_string_literal_95.data, _M0L10quals__arrS1142);
      moonbit_decref(_M0L10quals__arrS1142);
      #line 106 "/home/developer/Documents/1/fastq_io.mbt"
      _M0L6_2atmpS2760 = _M0MP26biolab8bio__seq3Seq3new(_M0L8seq__strS1138);
      moonbit_decref(_M0L8seq__strS1138);
      moonbit_incref(_M0L2idS1136);
      _M0L6_2atmpS2759
      = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
      Moonbit_object_header(_M0L6_2atmpS2759)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
      _M0L6_2atmpS2759->$0 = _M0L6_2atmpS2760;
      _M0L6_2atmpS2759->$1 = _M0L2idS1136;
      _M0L6_2atmpS2759->$2 = _M0L2idS1136;
      _M0L6_2atmpS2759->$3 = _M0L11descriptionS1137;
      _M0L6_2atmpS2759->$4 = _M0L11letter__annS1145;
      #line 105 "/home/developer/Documents/1/fastq_io.mbt"
      _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1110, _M0L6_2atmpS2759);
      moonbit_decref(_M0L6_2atmpS2759);
      _M0L3valS2762 = _M0L1iS1115->$0;
      _M0L6_2atmpS2761 = _M0L3valS2762 + 4;
      _M0L1iS1115->$0 = _M0L6_2atmpS2761;
      joinlet_4423:;
      continue;
    } else {
      moonbit_decref(_M0L1iS1115);
      moonbit_decref(_M0L5linesS1111);
    }
    break;
  }
  _result_4427.tag = 1;
  _result_4427.data.ok = _M0L7recordsS1110;
  return _result_4427;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__3(
  
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1109;
  int32_t _M0L6_2atmpS2702;
  moonbit_string_t _M0L6_2atmpS2703;
  struct moonbit_result_0 _tmp_4428;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2708;
  moonbit_string_t _M0L8_2afieldS3682;
  int32_t _M0L6_2acntS4094;
  moonbit_string_t _M0L2idS2706;
  moonbit_string_t _M0L6_2atmpS2707;
  struct moonbit_result_0 _tmp_4430;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2714;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3681;
  int32_t _M0L6_2acntS4100;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2713;
  moonbit_string_t _M0L6_2atmpS2711;
  moonbit_string_t _M0L6_2atmpS2712;
  struct moonbit_result_0 _result_4432;
  #line 37 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 38 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7recordsS1109
  = _M0FP26biolab8bio__seq12parse__fasta((moonbit_string_t)moonbit_string_literal_139.data);
  #line 39 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2702
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1109);
  _M0L6_2atmpS2703 = 0;
  #line 39 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4428
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS2702, 1, _M0L6_2atmpS2703, (moonbit_string_t)moonbit_string_literal_140.data);
  if (_M0L6_2atmpS2703) {
    moonbit_decref(_M0L6_2atmpS2703);
  }
  if (_tmp_4428.tag) {
    int32_t const _M0L5_2aokS2704 = _tmp_4428.data.ok;
  } else {
    void* const _M0L6_2aerrS2705 = _tmp_4428.data.err;
    struct moonbit_result_0 _result_4429;
    moonbit_decref(_M0L7recordsS1109);
    _result_4429.tag = 0;
    _result_4429.data.err = _M0L6_2aerrS2705;
    return _result_4429;
  }
  #line 40 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2708
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1109, 0);
  _M0L8_2afieldS3682 = _M0L6_2atmpS2708->$1;
  _M0L6_2acntS4094
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2708));
  if (_M0L6_2acntS4094 > 1) {
    int32_t _M0L11_2anew__cntS4099 = _M0L6_2acntS4094 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2708), _M0L11_2anew__cntS4099);
    moonbit_incref(_M0L8_2afieldS3682);
  } else if (_M0L6_2acntS4094 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4098 =
      _M0L6_2atmpS2708->$4;
    moonbit_string_t _M0L8_2afieldS4097;
    moonbit_string_t _M0L8_2afieldS4096;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4095;
    moonbit_decref(_M0L8_2afieldS4098);
    _M0L8_2afieldS4097 = _M0L6_2atmpS2708->$3;
    moonbit_decref(_M0L8_2afieldS4097);
    _M0L8_2afieldS4096 = _M0L6_2atmpS2708->$2;
    moonbit_decref(_M0L8_2afieldS4096);
    _M0L8_2afieldS4095 = _M0L6_2atmpS2708->$0;
    moonbit_decref(_M0L8_2afieldS4095);
    #line 40 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2708);
  }
  _M0L2idS2706 = _M0L8_2afieldS3682;
  _M0L6_2atmpS2707 = 0;
  #line 40 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4430
  = _M0FPB10assert__eqGsE(_M0L2idS2706, (moonbit_string_t)moonbit_string_literal_141.data, _M0L6_2atmpS2707, (moonbit_string_t)moonbit_string_literal_142.data);
  moonbit_decref(_M0L2idS2706);
  if (_M0L6_2atmpS2707) {
    moonbit_decref(_M0L6_2atmpS2707);
  }
  if (_tmp_4430.tag) {
    int32_t const _M0L5_2aokS2709 = _tmp_4430.data.ok;
  } else {
    void* const _M0L6_2aerrS2710 = _tmp_4430.data.err;
    struct moonbit_result_0 _result_4431;
    moonbit_decref(_M0L7recordsS1109);
    _result_4431.tag = 0;
    _result_4431.data.err = _M0L6_2aerrS2710;
    return _result_4431;
  }
  #line 41 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2714
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1109, 0);
  moonbit_decref(_M0L7recordsS1109);
  _M0L8_2afieldS3681 = _M0L6_2atmpS2714->$0;
  _M0L6_2acntS4100
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2714));
  if (_M0L6_2acntS4100 > 1) {
    int32_t _M0L11_2anew__cntS4105 = _M0L6_2acntS4100 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2714), _M0L11_2anew__cntS4105);
    moonbit_incref(_M0L8_2afieldS3681);
  } else if (_M0L6_2acntS4100 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4104 =
      _M0L6_2atmpS2714->$4;
    moonbit_string_t _M0L8_2afieldS4103;
    moonbit_string_t _M0L8_2afieldS4102;
    moonbit_string_t _M0L8_2afieldS4101;
    moonbit_decref(_M0L8_2afieldS4104);
    _M0L8_2afieldS4103 = _M0L6_2atmpS2714->$3;
    moonbit_decref(_M0L8_2afieldS4103);
    _M0L8_2afieldS4102 = _M0L6_2atmpS2714->$2;
    moonbit_decref(_M0L8_2afieldS4102);
    _M0L8_2afieldS4101 = _M0L6_2atmpS2714->$1;
    moonbit_decref(_M0L8_2afieldS4101);
    #line 41 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2714);
  }
  _M0L3seqS2713 = _M0L8_2afieldS3681;
  #line 41 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2711 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2713);
  moonbit_decref(_M0L3seqS2713);
  _M0L6_2atmpS2712 = 0;
  #line 41 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4432
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2711, (moonbit_string_t)moonbit_string_literal_143.data, _M0L6_2atmpS2712, (moonbit_string_t)moonbit_string_literal_144.data);
  moonbit_decref(_M0L6_2atmpS2711);
  if (_M0L6_2atmpS2712) {
    moonbit_decref(_M0L6_2atmpS2712);
  }
  return _result_4432;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__2(
  
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1108;
  int32_t _M0L6_2atmpS2684;
  moonbit_string_t _M0L6_2atmpS2685;
  struct moonbit_result_0 _tmp_4433;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2690;
  moonbit_string_t _M0L8_2afieldS3685;
  int32_t _M0L6_2acntS4106;
  moonbit_string_t _M0L2idS2688;
  moonbit_string_t _M0L6_2atmpS2689;
  struct moonbit_result_0 _tmp_4435;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2695;
  moonbit_string_t _M0L8_2afieldS3684;
  int32_t _M0L6_2acntS4112;
  moonbit_string_t _M0L11descriptionS2693;
  moonbit_string_t _M0L6_2atmpS2694;
  struct moonbit_result_0 _tmp_4437;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2701;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3683;
  int32_t _M0L6_2acntS4118;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2700;
  moonbit_string_t _M0L6_2atmpS2698;
  moonbit_string_t _M0L6_2atmpS2699;
  struct moonbit_result_0 _result_4439;
  #line 28 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 29 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7recordsS1108
  = _M0FP26biolab8bio__seq12parse__fasta((moonbit_string_t)moonbit_string_literal_145.data);
  #line 30 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2684
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1108);
  _M0L6_2atmpS2685 = 0;
  #line 30 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4433
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS2684, 1, _M0L6_2atmpS2685, (moonbit_string_t)moonbit_string_literal_146.data);
  if (_M0L6_2atmpS2685) {
    moonbit_decref(_M0L6_2atmpS2685);
  }
  if (_tmp_4433.tag) {
    int32_t const _M0L5_2aokS2686 = _tmp_4433.data.ok;
  } else {
    void* const _M0L6_2aerrS2687 = _tmp_4433.data.err;
    struct moonbit_result_0 _result_4434;
    moonbit_decref(_M0L7recordsS1108);
    _result_4434.tag = 0;
    _result_4434.data.err = _M0L6_2aerrS2687;
    return _result_4434;
  }
  #line 31 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2690
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1108, 0);
  _M0L8_2afieldS3685 = _M0L6_2atmpS2690->$1;
  _M0L6_2acntS4106
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2690));
  if (_M0L6_2acntS4106 > 1) {
    int32_t _M0L11_2anew__cntS4111 = _M0L6_2acntS4106 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2690), _M0L11_2anew__cntS4111);
    moonbit_incref(_M0L8_2afieldS3685);
  } else if (_M0L6_2acntS4106 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4110 =
      _M0L6_2atmpS2690->$4;
    moonbit_string_t _M0L8_2afieldS4109;
    moonbit_string_t _M0L8_2afieldS4108;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4107;
    moonbit_decref(_M0L8_2afieldS4110);
    _M0L8_2afieldS4109 = _M0L6_2atmpS2690->$3;
    moonbit_decref(_M0L8_2afieldS4109);
    _M0L8_2afieldS4108 = _M0L6_2atmpS2690->$2;
    moonbit_decref(_M0L8_2afieldS4108);
    _M0L8_2afieldS4107 = _M0L6_2atmpS2690->$0;
    moonbit_decref(_M0L8_2afieldS4107);
    #line 31 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2690);
  }
  _M0L2idS2688 = _M0L8_2afieldS3685;
  _M0L6_2atmpS2689 = 0;
  #line 31 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4435
  = _M0FPB10assert__eqGsE(_M0L2idS2688, (moonbit_string_t)moonbit_string_literal_147.data, _M0L6_2atmpS2689, (moonbit_string_t)moonbit_string_literal_148.data);
  moonbit_decref(_M0L2idS2688);
  if (_M0L6_2atmpS2689) {
    moonbit_decref(_M0L6_2atmpS2689);
  }
  if (_tmp_4435.tag) {
    int32_t const _M0L5_2aokS2691 = _tmp_4435.data.ok;
  } else {
    void* const _M0L6_2aerrS2692 = _tmp_4435.data.err;
    struct moonbit_result_0 _result_4436;
    moonbit_decref(_M0L7recordsS1108);
    _result_4436.tag = 0;
    _result_4436.data.err = _M0L6_2aerrS2692;
    return _result_4436;
  }
  #line 32 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2695
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1108, 0);
  _M0L8_2afieldS3684 = _M0L6_2atmpS2695->$3;
  _M0L6_2acntS4112
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2695));
  if (_M0L6_2acntS4112 > 1) {
    int32_t _M0L11_2anew__cntS4117 = _M0L6_2acntS4112 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2695), _M0L11_2anew__cntS4117);
    moonbit_incref(_M0L8_2afieldS3684);
  } else if (_M0L6_2acntS4112 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4116 =
      _M0L6_2atmpS2695->$4;
    moonbit_string_t _M0L8_2afieldS4115;
    moonbit_string_t _M0L8_2afieldS4114;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4113;
    moonbit_decref(_M0L8_2afieldS4116);
    _M0L8_2afieldS4115 = _M0L6_2atmpS2695->$2;
    moonbit_decref(_M0L8_2afieldS4115);
    _M0L8_2afieldS4114 = _M0L6_2atmpS2695->$1;
    moonbit_decref(_M0L8_2afieldS4114);
    _M0L8_2afieldS4113 = _M0L6_2atmpS2695->$0;
    moonbit_decref(_M0L8_2afieldS4113);
    #line 32 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2695);
  }
  _M0L11descriptionS2693 = _M0L8_2afieldS3684;
  _M0L6_2atmpS2694 = 0;
  #line 32 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4437
  = _M0FPB10assert__eqGsE(_M0L11descriptionS2693, (moonbit_string_t)moonbit_string_literal_147.data, _M0L6_2atmpS2694, (moonbit_string_t)moonbit_string_literal_149.data);
  moonbit_decref(_M0L11descriptionS2693);
  if (_M0L6_2atmpS2694) {
    moonbit_decref(_M0L6_2atmpS2694);
  }
  if (_tmp_4437.tag) {
    int32_t const _M0L5_2aokS2696 = _tmp_4437.data.ok;
  } else {
    void* const _M0L6_2aerrS2697 = _tmp_4437.data.err;
    struct moonbit_result_0 _result_4438;
    moonbit_decref(_M0L7recordsS1108);
    _result_4438.tag = 0;
    _result_4438.data.err = _M0L6_2aerrS2697;
    return _result_4438;
  }
  #line 33 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2701
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1108, 0);
  moonbit_decref(_M0L7recordsS1108);
  _M0L8_2afieldS3683 = _M0L6_2atmpS2701->$0;
  _M0L6_2acntS4118
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2701));
  if (_M0L6_2acntS4118 > 1) {
    int32_t _M0L11_2anew__cntS4123 = _M0L6_2acntS4118 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2701), _M0L11_2anew__cntS4123);
    moonbit_incref(_M0L8_2afieldS3683);
  } else if (_M0L6_2acntS4118 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4122 =
      _M0L6_2atmpS2701->$4;
    moonbit_string_t _M0L8_2afieldS4121;
    moonbit_string_t _M0L8_2afieldS4120;
    moonbit_string_t _M0L8_2afieldS4119;
    moonbit_decref(_M0L8_2afieldS4122);
    _M0L8_2afieldS4121 = _M0L6_2atmpS2701->$3;
    moonbit_decref(_M0L8_2afieldS4121);
    _M0L8_2afieldS4120 = _M0L6_2atmpS2701->$2;
    moonbit_decref(_M0L8_2afieldS4120);
    _M0L8_2afieldS4119 = _M0L6_2atmpS2701->$1;
    moonbit_decref(_M0L8_2afieldS4119);
    #line 33 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2701);
  }
  _M0L3seqS2700 = _M0L8_2afieldS3683;
  #line 33 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2698 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2700);
  moonbit_decref(_M0L3seqS2700);
  _M0L6_2atmpS2699 = 0;
  #line 33 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4439
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2698, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS2699, (moonbit_string_t)moonbit_string_literal_150.data);
  moonbit_decref(_M0L6_2atmpS2698);
  if (_M0L6_2atmpS2699) {
    moonbit_decref(_M0L6_2atmpS2699);
  }
  return _result_4439;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__1(
  
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1107;
  int32_t _M0L6_2atmpS2666;
  moonbit_string_t _M0L6_2atmpS2667;
  struct moonbit_result_0 _tmp_4440;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2672;
  moonbit_string_t _M0L8_2afieldS3688;
  int32_t _M0L6_2acntS4124;
  moonbit_string_t _M0L2idS2670;
  moonbit_string_t _M0L6_2atmpS2671;
  struct moonbit_result_0 _tmp_4442;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2677;
  moonbit_string_t _M0L8_2afieldS3687;
  int32_t _M0L6_2acntS4130;
  moonbit_string_t _M0L11descriptionS2675;
  moonbit_string_t _M0L6_2atmpS2676;
  struct moonbit_result_0 _tmp_4444;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2683;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3686;
  int32_t _M0L6_2acntS4136;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2682;
  moonbit_string_t _M0L6_2atmpS2680;
  moonbit_string_t _M0L6_2atmpS2681;
  struct moonbit_result_0 _result_4446;
  #line 19 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 20 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7recordsS1107
  = _M0FP26biolab8bio__seq12parse__fasta((moonbit_string_t)moonbit_string_literal_151.data);
  #line 21 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2666
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1107);
  _M0L6_2atmpS2667 = 0;
  #line 21 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4440
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS2666, 1, _M0L6_2atmpS2667, (moonbit_string_t)moonbit_string_literal_152.data);
  if (_M0L6_2atmpS2667) {
    moonbit_decref(_M0L6_2atmpS2667);
  }
  if (_tmp_4440.tag) {
    int32_t const _M0L5_2aokS2668 = _tmp_4440.data.ok;
  } else {
    void* const _M0L6_2aerrS2669 = _tmp_4440.data.err;
    struct moonbit_result_0 _result_4441;
    moonbit_decref(_M0L7recordsS1107);
    _result_4441.tag = 0;
    _result_4441.data.err = _M0L6_2aerrS2669;
    return _result_4441;
  }
  #line 22 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2672
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1107, 0);
  _M0L8_2afieldS3688 = _M0L6_2atmpS2672->$1;
  _M0L6_2acntS4124
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2672));
  if (_M0L6_2acntS4124 > 1) {
    int32_t _M0L11_2anew__cntS4129 = _M0L6_2acntS4124 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2672), _M0L11_2anew__cntS4129);
    moonbit_incref(_M0L8_2afieldS3688);
  } else if (_M0L6_2acntS4124 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4128 =
      _M0L6_2atmpS2672->$4;
    moonbit_string_t _M0L8_2afieldS4127;
    moonbit_string_t _M0L8_2afieldS4126;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4125;
    moonbit_decref(_M0L8_2afieldS4128);
    _M0L8_2afieldS4127 = _M0L6_2atmpS2672->$3;
    moonbit_decref(_M0L8_2afieldS4127);
    _M0L8_2afieldS4126 = _M0L6_2atmpS2672->$2;
    moonbit_decref(_M0L8_2afieldS4126);
    _M0L8_2afieldS4125 = _M0L6_2atmpS2672->$0;
    moonbit_decref(_M0L8_2afieldS4125);
    #line 22 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2672);
  }
  _M0L2idS2670 = _M0L8_2afieldS3688;
  _M0L6_2atmpS2671 = 0;
  #line 22 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4442
  = _M0FPB10assert__eqGsE(_M0L2idS2670, (moonbit_string_t)moonbit_string_literal_153.data, _M0L6_2atmpS2671, (moonbit_string_t)moonbit_string_literal_154.data);
  moonbit_decref(_M0L2idS2670);
  if (_M0L6_2atmpS2671) {
    moonbit_decref(_M0L6_2atmpS2671);
  }
  if (_tmp_4442.tag) {
    int32_t const _M0L5_2aokS2673 = _tmp_4442.data.ok;
  } else {
    void* const _M0L6_2aerrS2674 = _tmp_4442.data.err;
    struct moonbit_result_0 _result_4443;
    moonbit_decref(_M0L7recordsS1107);
    _result_4443.tag = 0;
    _result_4443.data.err = _M0L6_2aerrS2674;
    return _result_4443;
  }
  #line 23 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2677
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1107, 0);
  _M0L8_2afieldS3687 = _M0L6_2atmpS2677->$3;
  _M0L6_2acntS4130
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2677));
  if (_M0L6_2acntS4130 > 1) {
    int32_t _M0L11_2anew__cntS4135 = _M0L6_2acntS4130 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2677), _M0L11_2anew__cntS4135);
    moonbit_incref(_M0L8_2afieldS3687);
  } else if (_M0L6_2acntS4130 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4134 =
      _M0L6_2atmpS2677->$4;
    moonbit_string_t _M0L8_2afieldS4133;
    moonbit_string_t _M0L8_2afieldS4132;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4131;
    moonbit_decref(_M0L8_2afieldS4134);
    _M0L8_2afieldS4133 = _M0L6_2atmpS2677->$2;
    moonbit_decref(_M0L8_2afieldS4133);
    _M0L8_2afieldS4132 = _M0L6_2atmpS2677->$1;
    moonbit_decref(_M0L8_2afieldS4132);
    _M0L8_2afieldS4131 = _M0L6_2atmpS2677->$0;
    moonbit_decref(_M0L8_2afieldS4131);
    #line 23 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2677);
  }
  _M0L11descriptionS2675 = _M0L8_2afieldS3687;
  _M0L6_2atmpS2676 = 0;
  #line 23 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4444
  = _M0FPB10assert__eqGsE(_M0L11descriptionS2675, (moonbit_string_t)moonbit_string_literal_155.data, _M0L6_2atmpS2676, (moonbit_string_t)moonbit_string_literal_156.data);
  moonbit_decref(_M0L11descriptionS2675);
  if (_M0L6_2atmpS2676) {
    moonbit_decref(_M0L6_2atmpS2676);
  }
  if (_tmp_4444.tag) {
    int32_t const _M0L5_2aokS2678 = _tmp_4444.data.ok;
  } else {
    void* const _M0L6_2aerrS2679 = _tmp_4444.data.err;
    struct moonbit_result_0 _result_4445;
    moonbit_decref(_M0L7recordsS1107);
    _result_4445.tag = 0;
    _result_4445.data.err = _M0L6_2aerrS2679;
    return _result_4445;
  }
  #line 24 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2683
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1107, 0);
  moonbit_decref(_M0L7recordsS1107);
  _M0L8_2afieldS3686 = _M0L6_2atmpS2683->$0;
  _M0L6_2acntS4136
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2683));
  if (_M0L6_2acntS4136 > 1) {
    int32_t _M0L11_2anew__cntS4141 = _M0L6_2acntS4136 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2683), _M0L11_2anew__cntS4141);
    moonbit_incref(_M0L8_2afieldS3686);
  } else if (_M0L6_2acntS4136 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4140 =
      _M0L6_2atmpS2683->$4;
    moonbit_string_t _M0L8_2afieldS4139;
    moonbit_string_t _M0L8_2afieldS4138;
    moonbit_string_t _M0L8_2afieldS4137;
    moonbit_decref(_M0L8_2afieldS4140);
    _M0L8_2afieldS4139 = _M0L6_2atmpS2683->$3;
    moonbit_decref(_M0L8_2afieldS4139);
    _M0L8_2afieldS4138 = _M0L6_2atmpS2683->$2;
    moonbit_decref(_M0L8_2afieldS4138);
    _M0L8_2afieldS4137 = _M0L6_2atmpS2683->$1;
    moonbit_decref(_M0L8_2afieldS4137);
    #line 24 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2683);
  }
  _M0L3seqS2682 = _M0L8_2afieldS3686;
  #line 24 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2680 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2682);
  moonbit_decref(_M0L3seqS2682);
  _M0L6_2atmpS2681 = 0;
  #line 24 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4446
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2680, (moonbit_string_t)moonbit_string_literal_157.data, _M0L6_2atmpS2681, (moonbit_string_t)moonbit_string_literal_158.data);
  moonbit_decref(_M0L6_2atmpS2680);
  if (_M0L6_2atmpS2681) {
    moonbit_decref(_M0L6_2atmpS2681);
  }
  return _result_4446;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq45____test__736571696f5f7762746573742e6d6274__0(
  
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1106;
  int32_t _M0L6_2atmpS2637;
  moonbit_string_t _M0L6_2atmpS2638;
  struct moonbit_result_0 _tmp_4447;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2643;
  moonbit_string_t _M0L8_2afieldS3693;
  int32_t _M0L6_2acntS4142;
  moonbit_string_t _M0L2idS2641;
  moonbit_string_t _M0L6_2atmpS2642;
  struct moonbit_result_0 _tmp_4449;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2648;
  moonbit_string_t _M0L8_2afieldS3692;
  int32_t _M0L6_2acntS4148;
  moonbit_string_t _M0L11descriptionS2646;
  moonbit_string_t _M0L6_2atmpS2647;
  struct moonbit_result_0 _tmp_4451;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2654;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3691;
  int32_t _M0L6_2acntS4154;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2653;
  moonbit_string_t _M0L6_2atmpS2651;
  moonbit_string_t _M0L6_2atmpS2652;
  struct moonbit_result_0 _tmp_4453;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2659;
  moonbit_string_t _M0L8_2afieldS3690;
  int32_t _M0L6_2acntS4160;
  moonbit_string_t _M0L2idS2657;
  moonbit_string_t _M0L6_2atmpS2658;
  struct moonbit_result_0 _tmp_4455;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2665;
  struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS3689;
  int32_t _M0L6_2acntS4166;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2664;
  moonbit_string_t _M0L6_2atmpS2662;
  moonbit_string_t _M0L6_2atmpS2663;
  struct moonbit_result_0 _result_4457;
  #line 8 "/home/developer/Documents/1/seqio_wbtest.mbt"
  #line 9 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L7recordsS1106
  = _M0FP26biolab8bio__seq12parse__fasta((moonbit_string_t)moonbit_string_literal_20.data);
  #line 10 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2637
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1106);
  _M0L6_2atmpS2638 = 0;
  #line 10 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4447
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS2637, 2, _M0L6_2atmpS2638, (moonbit_string_t)moonbit_string_literal_159.data);
  if (_M0L6_2atmpS2638) {
    moonbit_decref(_M0L6_2atmpS2638);
  }
  if (_tmp_4447.tag) {
    int32_t const _M0L5_2aokS2639 = _tmp_4447.data.ok;
  } else {
    void* const _M0L6_2aerrS2640 = _tmp_4447.data.err;
    struct moonbit_result_0 _result_4448;
    moonbit_decref(_M0L7recordsS1106);
    _result_4448.tag = 0;
    _result_4448.data.err = _M0L6_2aerrS2640;
    return _result_4448;
  }
  #line 11 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2643
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1106, 0);
  _M0L8_2afieldS3693 = _M0L6_2atmpS2643->$1;
  _M0L6_2acntS4142
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2643));
  if (_M0L6_2acntS4142 > 1) {
    int32_t _M0L11_2anew__cntS4147 = _M0L6_2acntS4142 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2643), _M0L11_2anew__cntS4147);
    moonbit_incref(_M0L8_2afieldS3693);
  } else if (_M0L6_2acntS4142 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4146 =
      _M0L6_2atmpS2643->$4;
    moonbit_string_t _M0L8_2afieldS4145;
    moonbit_string_t _M0L8_2afieldS4144;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4143;
    moonbit_decref(_M0L8_2afieldS4146);
    _M0L8_2afieldS4145 = _M0L6_2atmpS2643->$3;
    moonbit_decref(_M0L8_2afieldS4145);
    _M0L8_2afieldS4144 = _M0L6_2atmpS2643->$2;
    moonbit_decref(_M0L8_2afieldS4144);
    _M0L8_2afieldS4143 = _M0L6_2atmpS2643->$0;
    moonbit_decref(_M0L8_2afieldS4143);
    #line 11 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2643);
  }
  _M0L2idS2641 = _M0L8_2afieldS3693;
  _M0L6_2atmpS2642 = 0;
  #line 11 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4449
  = _M0FPB10assert__eqGsE(_M0L2idS2641, (moonbit_string_t)moonbit_string_literal_22.data, _M0L6_2atmpS2642, (moonbit_string_t)moonbit_string_literal_160.data);
  moonbit_decref(_M0L2idS2641);
  if (_M0L6_2atmpS2642) {
    moonbit_decref(_M0L6_2atmpS2642);
  }
  if (_tmp_4449.tag) {
    int32_t const _M0L5_2aokS2644 = _tmp_4449.data.ok;
  } else {
    void* const _M0L6_2aerrS2645 = _tmp_4449.data.err;
    struct moonbit_result_0 _result_4450;
    moonbit_decref(_M0L7recordsS1106);
    _result_4450.tag = 0;
    _result_4450.data.err = _M0L6_2aerrS2645;
    return _result_4450;
  }
  #line 12 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2648
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1106, 0);
  _M0L8_2afieldS3692 = _M0L6_2atmpS2648->$3;
  _M0L6_2acntS4148
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2648));
  if (_M0L6_2acntS4148 > 1) {
    int32_t _M0L11_2anew__cntS4153 = _M0L6_2acntS4148 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2648), _M0L11_2anew__cntS4153);
    moonbit_incref(_M0L8_2afieldS3692);
  } else if (_M0L6_2acntS4148 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4152 =
      _M0L6_2atmpS2648->$4;
    moonbit_string_t _M0L8_2afieldS4151;
    moonbit_string_t _M0L8_2afieldS4150;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4149;
    moonbit_decref(_M0L8_2afieldS4152);
    _M0L8_2afieldS4151 = _M0L6_2atmpS2648->$2;
    moonbit_decref(_M0L8_2afieldS4151);
    _M0L8_2afieldS4150 = _M0L6_2atmpS2648->$1;
    moonbit_decref(_M0L8_2afieldS4150);
    _M0L8_2afieldS4149 = _M0L6_2atmpS2648->$0;
    moonbit_decref(_M0L8_2afieldS4149);
    #line 12 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2648);
  }
  _M0L11descriptionS2646 = _M0L8_2afieldS3692;
  _M0L6_2atmpS2647 = 0;
  #line 12 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4451
  = _M0FPB10assert__eqGsE(_M0L11descriptionS2646, (moonbit_string_t)moonbit_string_literal_24.data, _M0L6_2atmpS2647, (moonbit_string_t)moonbit_string_literal_161.data);
  moonbit_decref(_M0L11descriptionS2646);
  if (_M0L6_2atmpS2647) {
    moonbit_decref(_M0L6_2atmpS2647);
  }
  if (_tmp_4451.tag) {
    int32_t const _M0L5_2aokS2649 = _tmp_4451.data.ok;
  } else {
    void* const _M0L6_2aerrS2650 = _tmp_4451.data.err;
    struct moonbit_result_0 _result_4452;
    moonbit_decref(_M0L7recordsS1106);
    _result_4452.tag = 0;
    _result_4452.data.err = _M0L6_2aerrS2650;
    return _result_4452;
  }
  #line 13 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2654
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1106, 0);
  _M0L8_2afieldS3691 = _M0L6_2atmpS2654->$0;
  _M0L6_2acntS4154
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2654));
  if (_M0L6_2acntS4154 > 1) {
    int32_t _M0L11_2anew__cntS4159 = _M0L6_2acntS4154 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2654), _M0L11_2anew__cntS4159);
    moonbit_incref(_M0L8_2afieldS3691);
  } else if (_M0L6_2acntS4154 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4158 =
      _M0L6_2atmpS2654->$4;
    moonbit_string_t _M0L8_2afieldS4157;
    moonbit_string_t _M0L8_2afieldS4156;
    moonbit_string_t _M0L8_2afieldS4155;
    moonbit_decref(_M0L8_2afieldS4158);
    _M0L8_2afieldS4157 = _M0L6_2atmpS2654->$3;
    moonbit_decref(_M0L8_2afieldS4157);
    _M0L8_2afieldS4156 = _M0L6_2atmpS2654->$2;
    moonbit_decref(_M0L8_2afieldS4156);
    _M0L8_2afieldS4155 = _M0L6_2atmpS2654->$1;
    moonbit_decref(_M0L8_2afieldS4155);
    #line 13 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2654);
  }
  _M0L3seqS2653 = _M0L8_2afieldS3691;
  #line 13 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2651 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2653);
  moonbit_decref(_M0L3seqS2653);
  _M0L6_2atmpS2652 = 0;
  #line 13 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4453
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2651, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS2652, (moonbit_string_t)moonbit_string_literal_162.data);
  moonbit_decref(_M0L6_2atmpS2651);
  if (_M0L6_2atmpS2652) {
    moonbit_decref(_M0L6_2atmpS2652);
  }
  if (_tmp_4453.tag) {
    int32_t const _M0L5_2aokS2655 = _tmp_4453.data.ok;
  } else {
    void* const _M0L6_2aerrS2656 = _tmp_4453.data.err;
    struct moonbit_result_0 _result_4454;
    moonbit_decref(_M0L7recordsS1106);
    _result_4454.tag = 0;
    _result_4454.data.err = _M0L6_2aerrS2656;
    return _result_4454;
  }
  #line 14 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2659
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1106, 1);
  _M0L8_2afieldS3690 = _M0L6_2atmpS2659->$1;
  _M0L6_2acntS4160
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2659));
  if (_M0L6_2acntS4160 > 1) {
    int32_t _M0L11_2anew__cntS4165 = _M0L6_2acntS4160 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2659), _M0L11_2anew__cntS4165);
    moonbit_incref(_M0L8_2afieldS3690);
  } else if (_M0L6_2acntS4160 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4164 =
      _M0L6_2atmpS2659->$4;
    moonbit_string_t _M0L8_2afieldS4163;
    moonbit_string_t _M0L8_2afieldS4162;
    struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS4161;
    moonbit_decref(_M0L8_2afieldS4164);
    _M0L8_2afieldS4163 = _M0L6_2atmpS2659->$3;
    moonbit_decref(_M0L8_2afieldS4163);
    _M0L8_2afieldS4162 = _M0L6_2atmpS2659->$2;
    moonbit_decref(_M0L8_2afieldS4162);
    _M0L8_2afieldS4161 = _M0L6_2atmpS2659->$0;
    moonbit_decref(_M0L8_2afieldS4161);
    #line 14 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2659);
  }
  _M0L2idS2657 = _M0L8_2afieldS3690;
  _M0L6_2atmpS2658 = 0;
  #line 14 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _tmp_4455
  = _M0FPB10assert__eqGsE(_M0L2idS2657, (moonbit_string_t)moonbit_string_literal_27.data, _M0L6_2atmpS2658, (moonbit_string_t)moonbit_string_literal_163.data);
  moonbit_decref(_M0L2idS2657);
  if (_M0L6_2atmpS2658) {
    moonbit_decref(_M0L6_2atmpS2658);
  }
  if (_tmp_4455.tag) {
    int32_t const _M0L5_2aokS2660 = _tmp_4455.data.ok;
  } else {
    void* const _M0L6_2aerrS2661 = _tmp_4455.data.err;
    struct moonbit_result_0 _result_4456;
    moonbit_decref(_M0L7recordsS1106);
    _result_4456.tag = 0;
    _result_4456.data.err = _M0L6_2aerrS2661;
    return _result_4456;
  }
  #line 15 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2665
  = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1106, 1);
  moonbit_decref(_M0L7recordsS1106);
  _M0L8_2afieldS3689 = _M0L6_2atmpS2665->$0;
  _M0L6_2acntS4166
  = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS2665));
  if (_M0L6_2acntS4166 > 1) {
    int32_t _M0L11_2anew__cntS4171 = _M0L6_2acntS4166 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS2665), _M0L11_2anew__cntS4171);
    moonbit_incref(_M0L8_2afieldS3689);
  } else if (_M0L6_2acntS4166 == 1) {
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS4170 =
      _M0L6_2atmpS2665->$4;
    moonbit_string_t _M0L8_2afieldS4169;
    moonbit_string_t _M0L8_2afieldS4168;
    moonbit_string_t _M0L8_2afieldS4167;
    moonbit_decref(_M0L8_2afieldS4170);
    _M0L8_2afieldS4169 = _M0L6_2atmpS2665->$3;
    moonbit_decref(_M0L8_2afieldS4169);
    _M0L8_2afieldS4168 = _M0L6_2atmpS2665->$2;
    moonbit_decref(_M0L8_2afieldS4168);
    _M0L8_2afieldS4167 = _M0L6_2atmpS2665->$1;
    moonbit_decref(_M0L8_2afieldS4167);
    #line 15 "/home/developer/Documents/1/seqio_wbtest.mbt"
    moonbit_free(_M0L6_2atmpS2665);
  }
  _M0L3seqS2664 = _M0L8_2afieldS3689;
  #line 15 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _M0L6_2atmpS2662 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2664);
  moonbit_decref(_M0L3seqS2664);
  _M0L6_2atmpS2663 = 0;
  #line 15 "/home/developer/Documents/1/seqio_wbtest.mbt"
  _result_4457
  = _M0FPB10assert__eqGsE(_M0L6_2atmpS2662, (moonbit_string_t)moonbit_string_literal_29.data, _M0L6_2atmpS2663, (moonbit_string_t)moonbit_string_literal_164.data);
  moonbit_decref(_M0L6_2atmpS2662);
  if (_M0L6_2atmpS2663) {
    moonbit_decref(_M0L6_2atmpS2663);
  }
  return _result_4457;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq12parse__fasta(
  moonbit_string_t _M0L7contentS1078
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1076;
  moonbit_string_t _M0L7_2abindS1079;
  int32_t _M0L6_2atmpS2636;
  struct _M0TPC16string10StringView _M0L6_2atmpS2635;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2634;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS1077;
  moonbit_string_t _M0L6_2atmpS2633;
  struct _M0TPB8MutLocalGOsE* _M0L14current__titleS1080;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS2632;
  struct _M0TPB8MutLocalGRPB13StringBuilderE* _M0L12current__seqS1081;
  int32_t _M0L6_2atmpS2631;
  int32_t _M0L2gtS1082;
  int32_t _M0L6_2atmpS2630;
  int32_t _M0L2crS1083;
  int32_t _M0L6_2atmpS2629;
  int32_t _M0L2spS1084;
  int32_t _M0L6_2atmpS2628;
  int32_t _M0L3tabS1085;
  int32_t _M0L7_2abindS1086;
  int32_t _M0L2__S1087;
  moonbit_string_t _M0L5titleS1102;
  moonbit_string_t _M0L8_2afieldS3695;
  int32_t _M0L6_2acntS4174;
  moonbit_string_t _M0L7_2abindS1103;
  struct _M0TPB13StringBuilder* _M0L8_2afieldS3694;
  int32_t _M0L6_2acntS4172;
  struct _M0TPB13StringBuilder* _M0L3valS2627;
  moonbit_string_t _M0L6_2atmpS2626;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2625;
  #line 24 "/home/developer/Documents/1/fasta_io.mbt"
  #line 25 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7recordsS1076
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS1079 = (moonbit_string_t)moonbit_string_literal_89.data;
  _M0L6_2atmpS2636 = Moonbit_array_length(_M0L7_2abindS1079);
  _M0L6_2atmpS2635
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1079, .$1 = 0, .$2 = _M0L6_2atmpS2636
  };
  #line 26 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2634
  = _M0MPC16string6String5split(_M0L7contentS1078, _M0L6_2atmpS2635);
  moonbit_decref(_M0L6_2atmpS2635.$0);
  #line 26 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L5linesS1077
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2634);
  moonbit_decref(_M0L6_2atmpS2634);
  _M0L6_2atmpS2633 = 0;
  _M0L14current__titleS1080
  = (struct _M0TPB8MutLocalGOsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGOsE));
  Moonbit_object_header(_M0L14current__titleS1080)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 49, 0);
  _M0L14current__titleS1080->$0 = _M0L6_2atmpS2633;
  #line 28 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2632 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L12current__seqS1081
  = (struct _M0TPB8MutLocalGRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGRPB13StringBuilderE));
  Moonbit_object_header(_M0L12current__seqS1081)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 52, 0);
  _M0L12current__seqS1081->$0 = _M0L6_2atmpS2632;
  _M0L6_2atmpS2631 = 62;
  _M0L2gtS1082 = (uint16_t)_M0L6_2atmpS2631;
  _M0L6_2atmpS2630 = 13;
  _M0L2crS1083 = (uint16_t)_M0L6_2atmpS2630;
  _M0L6_2atmpS2629 = 32;
  _M0L2spS1084 = (uint16_t)_M0L6_2atmpS2629;
  _M0L6_2atmpS2628 = 9;
  _M0L3tabS1085 = (uint16_t)_M0L6_2atmpS2628;
  _M0L7_2abindS1086 = _M0L5linesS1077->$1;
  _M0L2__S1087 = 0;
  while (1) {
    if (_M0L2__S1087 < _M0L7_2abindS1086) {
      struct _M0TPC16string10StringView* _M0L3bufS2624 = _M0L5linesS1077->$0;
      struct _M0TPC16string10StringView _M0L4lineS1088 =
        _M0L3bufS2624[_M0L2__S1087];
      int32_t _M0L1nS1091;
      struct _M0TPB8MutLocalGiE* _M0L12trimmed__lenS1092;
      int32_t _if__result_4460;
      int32_t _M0L3valS2608;
      int32_t _M0L6_2atmpS2609;
      int32_t _M0L6_2atmpS2604;
      moonbit_incref(_M0L4lineS1088.$0);
      #line 35 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L1nS1091 = _M0MPC16string10StringView6length(_M0L4lineS1088);
      _M0L12trimmed__lenS1092
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L12trimmed__lenS1092)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L12trimmed__lenS1092->$0 = _M0L1nS1091;
      if (_M0L1nS1091 > 0) {
        int32_t _M0L6_2atmpS2606 = _M0L1nS1091 - 1;
        int32_t _M0L6_2atmpS2605;
        #line 37 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2605
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1088, _M0L6_2atmpS2606);
        _if__result_4460 = _M0L6_2atmpS2605 == _M0L2crS1083;
      } else {
        _if__result_4460 = 0;
      }
      if (_if__result_4460) {
        int32_t _M0L6_2atmpS2607 = _M0L1nS1091 - 1;
        _M0L12trimmed__lenS1092->$0 = _M0L6_2atmpS2607;
      }
      _M0L3valS2608 = _M0L12trimmed__lenS1092->$0;
      if (_M0L3valS2608 == 0) {
        moonbit_decref(_M0L12trimmed__lenS1092);
        moonbit_decref(_M0L4lineS1088.$0);
        goto join_1089;
      }
      #line 43 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L6_2atmpS2609
      = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1088, 0);
      if (_M0L6_2atmpS2609 == _M0L2gtS1082) {
        moonbit_string_t _M0L5titleS1094;
        moonbit_string_t _M0L7_2abindS1095 = _M0L14current__titleS1080->$0;
        struct _M0TPB13StringBuilder* _M0L3valS2612;
        moonbit_string_t _M0L6_2atmpS2611;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2610;
        int32_t _M0L3valS2615;
        moonbit_string_t _M0L6_2atmpS2614;
        moonbit_string_t _M0L6_2atmpS2613;
        moonbit_string_t _M0L6_2aoldS3697;
        struct _M0TPB13StringBuilder* _M0L6_2atmpS2619;
        struct _M0TPB13StringBuilder* _M0L6_2aoldS3696;
        if (_M0L7_2abindS1095 == 0) {
          
        } else {
          moonbit_string_t _M0L7_2aSomeS1096 = _M0L7_2abindS1095;
          moonbit_string_t _M0L8_2atitleS1097 = _M0L7_2aSomeS1096;
          moonbit_incref(_M0L8_2atitleS1097);
          _M0L5titleS1094 = _M0L8_2atitleS1097;
          goto join_1093;
        }
        goto joinlet_4461;
        join_1093:;
        _M0L3valS2612 = _M0L12current__seqS1081->$0;
        moonbit_incref(_M0L3valS2612);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2611 = _M0MPB13StringBuilder10to__string(_M0L3valS2612);
        moonbit_decref(_M0L3valS2612);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2610
        = _M0FP26biolab8bio__seq24fasta__title__to__record(_M0L5titleS1094, _M0L6_2atmpS2611);
        moonbit_decref(_M0L5titleS1094);
        moonbit_decref(_M0L6_2atmpS2611);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1076, _M0L6_2atmpS2610);
        moonbit_decref(_M0L6_2atmpS2610);
        joinlet_4461:;
        _M0L3valS2615 = _M0L12trimmed__lenS1092->$0;
        if (_M0L3valS2615 > 1) {
          int32_t _M0L3valS2618 = _M0L12trimmed__lenS1092->$0;
          int64_t _M0L6_2atmpS2617;
          struct _M0TPC16string10StringView _M0L6_2atmpS2616;
          moonbit_decref(_M0L12trimmed__lenS1092);
          _M0L6_2atmpS2617 = (int64_t)_M0L3valS2618;
          #line 51 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2616
          = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1088, 1, _M0L6_2atmpS2617);
          moonbit_decref(_M0L4lineS1088.$0);
          #line 51 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2614
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2616);
          moonbit_decref(_M0L6_2atmpS2616.$0);
        } else {
          moonbit_decref(_M0L12trimmed__lenS1092);
          moonbit_decref(_M0L4lineS1088.$0);
          _M0L6_2atmpS2614 = (moonbit_string_t)moonbit_string_literal_0.data;
        }
        _M0L6_2atmpS2613 = _M0L6_2atmpS2614;
        _M0L6_2aoldS3697 = _M0L14current__titleS1080->$0;
        if (_M0L6_2aoldS3697) {
          moonbit_decref(_M0L6_2aoldS3697);
        }
        _M0L14current__titleS1080->$0 = _M0L6_2atmpS2613;
        #line 56 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2619 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
        _M0L6_2aoldS3696 = _M0L12current__seqS1081->$0;
        moonbit_decref(_M0L6_2aoldS3696);
        _M0L12current__seqS1081->$0 = _M0L6_2atmpS2619;
      } else {
        int32_t _M0L1iS1098 = 0;
        while (1) {
          int32_t _M0L3valS2620 = _M0L12trimmed__lenS1092->$0;
          if (_M0L1iS1098 < _M0L3valS2620) {
            int32_t _M0L1cS1099;
            int32_t _if__result_4463;
            int32_t _M0L6_2atmpS2623;
            #line 59 "/home/developer/Documents/1/fasta_io.mbt"
            _M0L1cS1099
            = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1088, _M0L1iS1098);
            if (_M0L1cS1099 != _M0L2spS1084) {
              _if__result_4463 = _M0L1cS1099 != _M0L3tabS1085;
            } else {
              _if__result_4463 = 0;
            }
            if (_if__result_4463) {
              struct _M0TPB13StringBuilder* _M0L3valS2621 =
                _M0L12current__seqS1081->$0;
              int32_t _M0L6_2atmpS2622;
              moonbit_incref(_M0L3valS2621);
              #line 61 "/home/developer/Documents/1/fasta_io.mbt"
              _M0L6_2atmpS2622
              = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS1099);
              #line 61 "/home/developer/Documents/1/fasta_io.mbt"
              _M0IPB13StringBuilderPB6Logger11write__char(_M0L3valS2621, _M0L6_2atmpS2622);
              moonbit_decref(_M0L3valS2621);
            }
            _M0L6_2atmpS2623 = _M0L1iS1098 + 1;
            _M0L1iS1098 = _M0L6_2atmpS2623;
            continue;
          } else {
            moonbit_decref(_M0L12trimmed__lenS1092);
            moonbit_decref(_M0L4lineS1088.$0);
          }
          break;
        }
      }
      goto join_1089;
      goto joinlet_4459;
      join_1089:;
      _M0L6_2atmpS2604 = _M0L2__S1087 + 1;
      _M0L2__S1087 = _M0L6_2atmpS2604;
      continue;
      joinlet_4459:;
    } else {
      moonbit_decref(_M0L5linesS1077);
    }
    break;
  }
  _M0L8_2afieldS3695 = _M0L14current__titleS1080->$0;
  _M0L6_2acntS4174
  = Moonbit_rc_count(Moonbit_object_header(_M0L14current__titleS1080));
  if (_M0L6_2acntS4174 > 1) {
    int32_t _M0L11_2anew__cntS4175 = _M0L6_2acntS4174 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L14current__titleS1080), _M0L11_2anew__cntS4175);
    if (_M0L8_2afieldS3695) {
      moonbit_incref(_M0L8_2afieldS3695);
    }
  } else if (_M0L6_2acntS4174 == 1) {
    #line 66 "/home/developer/Documents/1/fasta_io.mbt"
    moonbit_free(_M0L14current__titleS1080);
  }
  _M0L7_2abindS1103 = _M0L8_2afieldS3695;
  if (_M0L7_2abindS1103 == 0) {
    if (_M0L7_2abindS1103) {
      moonbit_decref(_M0L7_2abindS1103);
    }
    moonbit_decref(_M0L12current__seqS1081);
  } else {
    moonbit_string_t _M0L7_2aSomeS1104 = _M0L7_2abindS1103;
    moonbit_string_t _M0L8_2atitleS1105 = _M0L7_2aSomeS1104;
    _M0L5titleS1102 = _M0L8_2atitleS1105;
    goto join_1101;
  }
  goto joinlet_4464;
  join_1101:;
  _M0L8_2afieldS3694 = _M0L12current__seqS1081->$0;
  _M0L6_2acntS4172
  = Moonbit_rc_count(Moonbit_object_header(_M0L12current__seqS1081));
  if (_M0L6_2acntS4172 > 1) {
    int32_t _M0L11_2anew__cntS4173 = _M0L6_2acntS4172 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L12current__seqS1081), _M0L11_2anew__cntS4173);
    moonbit_incref(_M0L8_2afieldS3694);
  } else if (_M0L6_2acntS4172 == 1) {
    #line 68 "/home/developer/Documents/1/fasta_io.mbt"
    moonbit_free(_M0L12current__seqS1081);
  }
  _M0L3valS2627 = _M0L8_2afieldS3694;
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2626 = _M0MPB13StringBuilder10to__string(_M0L3valS2627);
  moonbit_decref(_M0L3valS2627);
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2625
  = _M0FP26biolab8bio__seq24fasta__title__to__record(_M0L5titleS1102, _M0L6_2atmpS2626);
  moonbit_decref(_M0L5titleS1102);
  moonbit_decref(_M0L6_2atmpS2626);
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1076, _M0L6_2atmpS2625);
  moonbit_decref(_M0L6_2atmpS2625);
  joinlet_4464:;
  return _M0L7recordsS1076;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0FP26biolab8bio__seq24fasta__title__to__record(
  moonbit_string_t _M0L5titleS1071,
  moonbit_string_t _M0L8seq__strS1067
) {
  moonbit_string_t _M0L2idS1065;
  moonbit_string_t _M0L11descriptionS1066;
  int32_t _M0L3posS1070;
  moonbit_string_t _M0L7_2abindS1073;
  int32_t _M0L6_2atmpS2603;
  struct _M0TPC16string10StringView _M0L6_2atmpS2602;
  int64_t _M0L7_2abindS1072;
  int64_t _M0L6_2atmpS2601;
  struct _M0TPC16string10StringView _M0L6_2atmpS2600;
  moonbit_string_t _M0L6_2atmpS2599;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2595;
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1068;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2598;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2597;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2596;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_4467;
  #line 78 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7_2abindS1073 = (moonbit_string_t)moonbit_string_literal_115.data;
  _M0L6_2atmpS2603 = Moonbit_array_length(_M0L7_2abindS1073);
  _M0L6_2atmpS2602
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1073, .$1 = 0, .$2 = _M0L6_2atmpS2603
  };
  #line 79 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7_2abindS1072
  = _M0MPC16string6String4find(_M0L5titleS1071, _M0L6_2atmpS2602);
  moonbit_decref(_M0L6_2atmpS2602.$0);
  if (_M0L7_2abindS1072 == 4294967296ll) {
    moonbit_incref(_M0L5titleS1071);
    moonbit_incref(_M0L5titleS1071);
    _M0L2idS1065 = _M0L5titleS1071;
    _M0L11descriptionS1066 = _M0L5titleS1071;
    goto join_1064;
  } else {
    int64_t _M0L7_2aSomeS1074 = _M0L7_2abindS1072;
    int32_t _M0L6_2aposS1075 = (int32_t)_M0L7_2aSomeS1074;
    _M0L3posS1070 = _M0L6_2aposS1075;
    goto join_1069;
  }
  join_1069:;
  _M0L6_2atmpS2601 = (int64_t)_M0L3posS1070;
  #line 80 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2600
  = _M0MPC16string6String11sub_2einner(_M0L5titleS1071, 0, _M0L6_2atmpS2601);
  #line 80 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2599 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2600);
  moonbit_decref(_M0L6_2atmpS2600.$0);
  moonbit_incref(_M0L5titleS1071);
  _M0L2idS1065 = _M0L6_2atmpS2599;
  _M0L11descriptionS1066 = _M0L5titleS1071;
  goto join_1064;
  join_1064:;
  #line 84 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2595 = _M0MP26biolab8bio__seq3Seq3new(_M0L8seq__strS1067);
  _M0L7_2abindS1068 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2598 = _M0L7_2abindS1068;
  _M0L6_2atmpS2597
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2598, .$1 = 0, .$2 = 0
  };
  #line 88 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2596 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2597, 0ll);
  moonbit_decref(_M0L6_2atmpS2597.$0);
  moonbit_incref(_M0L2idS1065);
  _block_4467
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_4467)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
  _block_4467->$0 = _M0L6_2atmpS2595;
  _block_4467->$1 = _M0L2idS1065;
  _block_4467->$2 = _M0L2idS1065;
  _block_4467->$3 = _M0L11descriptionS1066;
  _block_4467->$4 = _M0L6_2atmpS2596;
  return _block_4467;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord11new_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS1059,
  moonbit_string_t _M0L2idS1060,
  moonbit_string_t _M0L4nameS1061,
  moonbit_string_t _M0L11descriptionS1062
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1063;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2594;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2593;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2592;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_4468;
  #line 31 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS1063 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2594 = _M0L7_2abindS1063;
  _M0L6_2atmpS2593
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2594, .$1 = 0, .$2 = 0
  };
  #line 42 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS2592 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2593, 0ll);
  moonbit_decref(_M0L6_2atmpS2593.$0);
  moonbit_incref(_M0L3seqS1059);
  moonbit_incref(_M0L2idS1060);
  moonbit_incref(_M0L4nameS1061);
  moonbit_incref(_M0L11descriptionS1062);
  _block_4468
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_4468)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
  _block_4468->$0 = _M0L3seqS1059;
  _block_4468->$1 = _M0L2idS1060;
  _block_4468->$2 = _M0L4nameS1061;
  _block_4468->$3 = _M0L11descriptionS1062;
  _block_4468->$4 = _M0L6_2atmpS2592;
  return _block_4468;
}

int32_t _M0MP26biolab8bio__seq9SeqRecord6length(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS1058
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2591;
  int32_t _result_4469;
  #line 84 "/home/developer/Documents/1/seq_record.mbt"
  _M0L3seqS2591 = _M0L4selfS1058->$0;
  moonbit_incref(_M0L3seqS2591);
  #line 85 "/home/developer/Documents/1/seq_record.mbt"
  _result_4469 = _M0MP26biolab8bio__seq3Seq6length(_M0L3seqS2591);
  moonbit_decref(_M0L3seqS2591);
  return _result_4469;
}

int32_t _M0MP26biolab8bio__seq3Seq6length(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS1057
) {
  moonbit_string_t _M0L4dataS2590;
  int32_t _result_4470;
  #line 86 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2590 = _M0L4selfS1057->$0;
  moonbit_incref(_M0L4dataS2590);
  #line 87 "/home/developer/Documents/1/seq.mbt"
  _result_4470
  = _M0MPC16string6String20char__length_2einner(_M0L4dataS2590, 0, 4294967296ll);
  moonbit_decref(_M0L4dataS2590);
  return _result_4470;
}

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS1056
) {
  moonbit_string_t _M0L8_2afieldS3705;
  #line 74 "/home/developer/Documents/1/seq.mbt"
  _M0L8_2afieldS3705 = _M0L4selfS1056->$0;
  moonbit_incref(_M0L8_2afieldS3705);
  return _M0L8_2afieldS3705;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t _M0L1sS1055
) {
  struct _M0TP26biolab8bio__seq3Seq* _block_4471;
  #line 62 "/home/developer/Documents/1/seq.mbt"
  moonbit_incref(_M0L1sS1055);
  _block_4471
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_4471)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 55, 0);
  _block_4471->$0 = _M0L1sS1055;
  return _block_4471;
}

int32_t _M0IPC15array5ArrayPB2Eq5equalGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS1050,
  struct _M0TPB5ArrayGiE* _M0L5otherS1052
) {
  int32_t _M0L9self__lenS1049;
  int32_t _M0L10other__lenS1051;
  #line 288 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L9self__lenS1049 = _M0L4selfS1050->$1;
  _M0L10other__lenS1051 = _M0L5otherS1052->$1;
  if (_M0L9self__lenS1049 == _M0L10other__lenS1051) {
    int32_t _M0L1iS1053 = 0;
    while (1) {
      if (_M0L1iS1053 < _M0L9self__lenS1049) {
        int32_t* _M0L3bufS2588 = _M0L4selfS1050->$0;
        int32_t _M0L6_2atmpS2585 = (int32_t)_M0L3bufS2588[_M0L1iS1053];
        int32_t* _M0L3bufS2587 = _M0L5otherS1052->$0;
        int32_t _M0L6_2atmpS2586 = (int32_t)_M0L3bufS2587[_M0L1iS1053];
        int32_t _M0L6_2atmpS2589;
        if (_M0L6_2atmpS2585 == _M0L6_2atmpS2586) {
          
        } else {
          return 0;
        }
        _M0L6_2atmpS2589 = _M0L1iS1053 + 1;
        _M0L1iS1053 = _M0L6_2atmpS2589;
        continue;
      } else {
        return 1;
      }
      break;
    }
  } else {
    return 0;
  }
}

moonbit_string_t _M0MPC15array5Array2atGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS1038,
  int32_t _M0L5indexS1039
) {
  int32_t _M0L3lenS1037;
  int32_t _if__result_4473;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1037 = _M0L4selfS1038->$1;
  if (_M0L5indexS1039 >= 0) {
    _if__result_4473 = _M0L5indexS1039 < _M0L3lenS1037;
  } else {
    _if__result_4473 = 0;
  }
  if (_if__result_4473) {
    moonbit_string_t* _M0L6_2atmpS2581;
    moonbit_string_t _M0L6_2atmpS3708;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2581 = _M0MPC15array5Array6bufferGsE(_M0L4selfS1038);
    _M0L6_2atmpS3708 = (moonbit_string_t)_M0L6_2atmpS2581[_M0L5indexS1039];
    moonbit_incref(_M0L6_2atmpS3708);
    moonbit_decref(_M0L6_2atmpS2581);
    return _M0L6_2atmpS3708;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS1041,
  int32_t _M0L5indexS1042
) {
  int32_t _M0L3lenS1040;
  int32_t _if__result_4474;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1040 = _M0L4selfS1041->$1;
  if (_M0L5indexS1042 >= 0) {
    _if__result_4474 = _M0L5indexS1042 < _M0L3lenS1040;
  } else {
    _if__result_4474 = 0;
  }
  if (_if__result_4474) {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2582;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3709;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2582
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS1041);
    _M0L6_2atmpS3709
    = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L6_2atmpS2582[
        _M0L5indexS1042
      ];
    if (_M0L6_2atmpS3709) {
      moonbit_incref(_M0L6_2atmpS3709);
    }
    moonbit_decref(_M0L6_2atmpS2582);
    return _M0L6_2atmpS3709;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TPC16string10StringView _M0MPC15array5Array2atGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS1044,
  int32_t _M0L5indexS1045
) {
  int32_t _M0L3lenS1043;
  int32_t _if__result_4475;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1043 = _M0L4selfS1044->$1;
  if (_M0L5indexS1045 >= 0) {
    _if__result_4475 = _M0L5indexS1045 < _M0L3lenS1043;
  } else {
    _if__result_4475 = 0;
  }
  if (_if__result_4475) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS2583;
    struct _M0TPC16string10StringView _M0L6_2atmpS3710;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2583
    = _M0MPC15array5Array6bufferGRPC16string10StringViewE(_M0L4selfS1044);
    _M0L6_2atmpS3710 = _M0L6_2atmpS2583[_M0L5indexS1045];
    moonbit_incref(_M0L6_2atmpS3710.$0);
    moonbit_decref(_M0L6_2atmpS2583);
    return _M0L6_2atmpS3710;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array5Array2atGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS1047,
  int32_t _M0L5indexS1048
) {
  int32_t _M0L3lenS1046;
  int32_t _if__result_4476;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1046 = _M0L4selfS1047->$1;
  if (_M0L5indexS1048 >= 0) {
    _if__result_4476 = _M0L5indexS1048 < _M0L3lenS1046;
  } else {
    _if__result_4476 = 0;
  }
  if (_if__result_4476) {
    int32_t* _M0L6_2atmpS2584;
    int32_t _result_4477;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2584 = _M0MPC15array5Array6bufferGiE(_M0L4selfS1047);
    _result_4477 = (int32_t)_M0L6_2atmpS2584[_M0L5indexS1048];
    moonbit_decref(_M0L6_2atmpS2584);
    return _result_4477;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0IPB9SourceLocPB4Show6output(
  moonbit_string_t _M0L4selfS1036,
  struct _M0TPB6Logger _M0L6loggerS1035
) {
  moonbit_string_t _M0L6_2atmpS2580;
  #line 43 "/home/developer/.moon/lib/core/builtin/autoloc.mbt"
  moonbit_incref(_M0L4selfS1036);
  _M0L6_2atmpS2580 = _M0L4selfS1036;
  #line 44 "/home/developer/.moon/lib/core/builtin/autoloc.mbt"
  _M0L6loggerS1035.$0->$method_0(_M0L6loggerS1035.$1, _M0L6_2atmpS2580);
  moonbit_decref(_M0L6_2atmpS2580);
  return 0;
}

int32_t _M0FPB7printlnGsE(moonbit_string_t _M0L5inputS1034) {
  moonbit_string_t _M0L6_2atmpS2579;
  #line 36 "/home/developer/.moon/lib/core/builtin/console.mbt"
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  _M0L6_2atmpS2579
  = _M0IPC16string6StringPB4Show10to__string(_M0L5inputS1034);
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  moonbit_println(_M0L6_2atmpS2579);
  moonbit_decref(_M0L6_2atmpS2579);
  return 0;
}

int32_t _M0IPC13int3IntPB4Hash4hash(int32_t _M0L4selfS1033) {
  int32_t _tmp_4478;
  uint32_t _M0L6_2atmpS2578;
  uint32_t _M0L6_2atmpS2577;
  uint32_t _M0L3accS1032;
  uint32_t _M0L6_2atmpS2576;
  uint32_t _M0L6_2atmpS2575;
  #line 588 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_4478 = 0;
  _M0L6_2atmpS2578 = *(uint32_t*)&_tmp_4478;
  _M0L6_2atmpS2577 = _M0L6_2atmpS2578 + 374761393u;
  _M0L3accS1032 = _M0L6_2atmpS2577 + 4u;
  _M0L6_2atmpS2576 = *(uint32_t*)&_M0L4selfS1033;
  #line 590 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS2575 = _M0FPB13consume4__acc(_M0L3accS1032, _M0L6_2atmpS2576);
  #line 590 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS2575);
}

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t _M0L4selfS1028) {
  int32_t _tmp_4479;
  uint32_t _M0L6_2atmpS2574;
  uint32_t _M0Lm3accS1026;
  int32_t _M0L7_2abindS1027;
  int32_t _M0L1iS1029;
  uint32_t _M0L6_2atmpS2573;
  #line 522 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_4479 = 0;
  _M0L6_2atmpS2574 = *(uint32_t*)&_tmp_4479;
  _M0Lm3accS1026 = _M0L6_2atmpS2574 + 374761393u;
  _M0L7_2abindS1027 = Moonbit_array_length(_M0L4selfS1028);
  _M0L1iS1029 = 0;
  while (1) {
    if (_M0L1iS1029 < _M0L7_2abindS1027) {
      uint32_t _M0L6_2atmpS2568 = _M0Lm3accS1026;
      int32_t _M0L6_2atmpS2571;
      int32_t _M0L6_2atmpS2570;
      uint32_t _M0L1vS1030;
      uint32_t _M0L6_2atmpS2569;
      int32_t _M0L6_2atmpS2572;
      _M0Lm3accS1026 = _M0L6_2atmpS2568 + 4u;
      _M0L6_2atmpS2571 = _M0L4selfS1028[_M0L1iS1029];
      _M0L6_2atmpS2570 = (int32_t)_M0L6_2atmpS2571;
      _M0L1vS1030 = *(uint32_t*)&_M0L6_2atmpS2570;
      _M0L6_2atmpS2569 = _M0Lm3accS1026;
      #line 527 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
      _M0Lm3accS1026 = _M0FPB13consume4__acc(_M0L6_2atmpS2569, _M0L1vS1030);
      _M0L6_2atmpS2572 = _M0L1iS1029 + 1;
      _M0L1iS1029 = _M0L6_2atmpS2572;
      continue;
    }
    break;
  }
  _M0L6_2atmpS2573 = _M0Lm3accS1026;
  #line 529 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS2573);
}

struct _M0TUsRPB5ArrayGiEE* _M0MPB5Iter24nextGsRPB5ArrayGiEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0L4selfS1025
) {
  #line 1099 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  #line 1100 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  return _M0MPB4Iter4nextGUsRPB5ArrayGiEEE(_M0L4selfS1025);
}

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB3Map5iter2GsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS1024
) {
  #line 725 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 727 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  return _M0MPB3Map4iterGsRPB5ArrayGiEE(_M0L4selfS1024);
}

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB3Map4iterGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS1014
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4headS2567;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE* _M0L11curr__entryS1013;
  int32_t _M0L3lenS1015;
  struct _M0TPB8MutLocalGiE* _M0L9remainingS1016;
  struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__* _closure_4481;
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L6_2atmpS2558;
  int64_t _M0L6_2atmpS2559;
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _result_4482;
  #line 705 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4headS2567 = _M0L4selfS1014->$5;
  if (_M0L4headS2567) {
    moonbit_incref(_M0L4headS2567);
  }
  _M0L11curr__entryS1013
  = (struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE));
  Moonbit_object_header(_M0L11curr__entryS1013)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 58, 0);
  _M0L11curr__entryS1013->$0 = _M0L4headS2567;
  _M0L3lenS1015 = _M0L4selfS1014->$1;
  _M0L9remainingS1016
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L9remainingS1016)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L9remainingS1016->$0 = _M0L3lenS1015;
  _closure_4481
  = (struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__*)moonbit_malloc(sizeof(struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__));
  Moonbit_object_header(_closure_4481)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 61, 0);
  _closure_4481->code = &_M0MPB3Map4iterGsRPB5ArrayGiEEC2560l711;
  _closure_4481->$0 = _M0L9remainingS1016;
  _closure_4481->$1 = _M0L11curr__entryS1013;
  _M0L6_2atmpS2558 = (struct _M0TWEOUsRPB5ArrayGiEE*)_closure_4481;
  _M0L6_2atmpS2559 = (int64_t)_M0L3lenS1015;
  #line 710 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _result_4482
  = _M0MPB4Iter3newGUsRPB5ArrayGiEEE(_M0L6_2atmpS2558, _M0L6_2atmpS2559);
  moonbit_decref(_M0L6_2atmpS2558);
  return _result_4482;
}

struct _M0TUsRPB5ArrayGiEE* _M0MPB3Map4iterGsRPB5ArrayGiEEC2560l711(
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L6_2aenvS2561
) {
  struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__* _M0L14_2acasted__envS2562;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB5ArrayGiEEE* _M0L11curr__entryS1013;
  struct _M0TPB8MutLocalGiE* _M0L9remainingS1016;
  int32_t _M0L3valS2563;
  #line 711 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14_2acasted__envS2562
  = (struct _M0R104Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fArray_5bInt_5d_5d_7c_2eanon__u2560__l711__*)_M0L6_2aenvS2561;
  _M0L11curr__entryS1013 = _M0L14_2acasted__envS2562->$1;
  _M0L9remainingS1016 = _M0L14_2acasted__envS2562->$0;
  _M0L3valS2563 = _M0L9remainingS1016->$0;
  if (_M0L3valS2563 > 0) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS1018 =
      _M0L11curr__entryS1013->$0;
    if (_M0L7_2abindS1018 == 0) {
      goto join_1017;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS1019 =
        _M0L7_2abindS1018;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4_2axS1020 = _M0L7_2aSomeS1019;
      moonbit_string_t _M0L6_2akeyS1021 = _M0L4_2axS1020->$4;
      struct _M0TPB5ArrayGiE* _M0L8_2avalueS1022 = _M0L4_2axS1020->$5;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2anextS1023 =
        _M0L4_2axS1020->$1;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3712 =
        _M0L11curr__entryS1013->$0;
      int32_t _M0L3valS2565;
      int32_t _M0L6_2atmpS2564;
      struct _M0TUsRPB5ArrayGiEE* _M0L8_2atupleS2566;
      if (_M0L7_2anextS1023) {
        moonbit_incref(_M0L7_2anextS1023);
      }
      moonbit_incref(_M0L8_2avalueS1022);
      moonbit_incref(_M0L6_2akeyS1021);
      if (_M0L6_2aoldS3712) {
        moonbit_decref(_M0L6_2aoldS3712);
      }
      _M0L11curr__entryS1013->$0 = _M0L7_2anextS1023;
      _M0L3valS2565 = _M0L9remainingS1016->$0;
      _M0L6_2atmpS2564 = _M0L3valS2565 - 1;
      _M0L9remainingS1016->$0 = _M0L6_2atmpS2564;
      _M0L8_2atupleS2566
      = (struct _M0TUsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TUsRPB5ArrayGiEE));
      Moonbit_object_header(_M0L8_2atupleS2566)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 65, 0);
      _M0L8_2atupleS2566->$0 = _M0L6_2akeyS1021;
      _M0L8_2atupleS2566->$1 = _M0L8_2avalueS1022;
      return _M0L8_2atupleS2566;
    }
  } else {
    goto join_1017;
  }
  join_1017:;
  return 0;
}

int32_t _M0MPB3Map6lengthGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS1012
) {
  #line 641 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  return _M0L4selfS1012->$1;
}

struct _M0TUWEuQRPC15error5ErrorNsE* _M0MPB3Map3getGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS981,
  int32_t _M0L3keyS977
) {
  int32_t _M0L4hashS976;
  int32_t _M0L14capacity__maskS2515;
  int32_t _M0L6_2atmpS2514;
  int32_t _M0L1iS978;
  int32_t _M0L3idxS979;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS976 = _M0IPC13int3IntPB4Hash4hash(_M0L3keyS977);
  _M0L14capacity__maskS2515 = _M0L4selfS981->$3;
  _M0L6_2atmpS2514 = _M0L4hashS976 & _M0L14capacity__maskS2515;
  _M0L1iS978 = 0;
  _M0L3idxS979 = _M0L6_2atmpS2514;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2513 =
      _M0L4selfS981->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS980;
    if (
      _M0L3idxS979 < 0
      || _M0L3idxS979 >= Moonbit_array_length(_M0L7entriesS2513)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS980
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2513[
        _M0L3idxS979
      ];
    if (_M0L7_2abindS980 == 0) {
      struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2502 = 0;
      return _M0L6_2atmpS2502;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS982 =
        _M0L7_2abindS980;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L8_2aentryS983 =
        _M0L7_2aSomeS982;
      int32_t _M0L4hashS2504 = _M0L8_2aentryS983->$3;
      int32_t _if__result_4485;
      int32_t _M0L3pslS2507;
      int32_t _M0L6_2atmpS2509;
      int32_t _M0L6_2atmpS2511;
      int32_t _M0L14capacity__maskS2512;
      int32_t _M0L6_2atmpS2510;
      if (_M0L4hashS2504 == _M0L4hashS976) {
        int32_t _M0L3keyS2503 = _M0L8_2aentryS983->$4;
        _if__result_4485 = _M0L3keyS2503 == _M0L3keyS977;
      } else {
        _if__result_4485 = 0;
      }
      if (_if__result_4485) {
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L5valueS2506 =
          _M0L8_2aentryS983->$5;
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2505;
        moonbit_incref(_M0L5valueS2506);
        _M0L6_2atmpS2505 = _M0L5valueS2506;
        return _M0L6_2atmpS2505;
      } else {
        moonbit_incref(_M0L8_2aentryS983);
      }
      _M0L3pslS2507 = _M0L8_2aentryS983->$2;
      moonbit_decref(_M0L8_2aentryS983);
      if (_M0L1iS978 > _M0L3pslS2507) {
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2508 = 0;
        return _M0L6_2atmpS2508;
      }
      _M0L6_2atmpS2509 = _M0L1iS978 + 1;
      _M0L6_2atmpS2511 = _M0L3idxS979 + 1;
      _M0L14capacity__maskS2512 = _M0L4selfS981->$3;
      _M0L6_2atmpS2510 = _M0L6_2atmpS2511 & _M0L14capacity__maskS2512;
      _M0L1iS978 = _M0L6_2atmpS2509;
      _M0L3idxS979 = _M0L6_2atmpS2510;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3getGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS990,
  moonbit_string_t _M0L3keyS986
) {
  int32_t _M0L4hashS985;
  int32_t _M0L14capacity__maskS2529;
  int32_t _M0L6_2atmpS2528;
  int32_t _M0L1iS987;
  int32_t _M0L3idxS988;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS985 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS986);
  _M0L14capacity__maskS2529 = _M0L4selfS990->$3;
  _M0L6_2atmpS2528 = _M0L4hashS985 & _M0L14capacity__maskS2529;
  _M0L1iS987 = 0;
  _M0L3idxS988 = _M0L6_2atmpS2528;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2527 =
      _M0L4selfS990->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS989;
    if (
      _M0L3idxS988 < 0
      || _M0L3idxS988 >= Moonbit_array_length(_M0L7entriesS2527)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS989
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2527[
        _M0L3idxS988
      ];
    if (_M0L7_2abindS989 == 0) {
      struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2516 = 0;
      return _M0L6_2atmpS2516;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS991 =
        _M0L7_2abindS989;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2aentryS992 =
        _M0L7_2aSomeS991;
      int32_t _M0L4hashS2518 = _M0L8_2aentryS992->$3;
      int32_t _if__result_4487;
      int32_t _M0L3pslS2521;
      int32_t _M0L6_2atmpS2523;
      int32_t _M0L6_2atmpS2525;
      int32_t _M0L14capacity__maskS2526;
      int32_t _M0L6_2atmpS2524;
      if (_M0L4hashS2518 == _M0L4hashS985) {
        moonbit_string_t _M0L3keyS2517 = _M0L8_2aentryS992->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_4487
        = _M0L3keyS2517 == _M0L3keyS986
          || Moonbit_array_length(_M0L3keyS2517)
             == Moonbit_array_length(_M0L3keyS986)
             && 0
                == memcmp(_M0L3keyS2517, _M0L3keyS986, Moonbit_array_length(_M0L3keyS2517) * 2);
      } else {
        _if__result_4487 = 0;
      }
      if (_if__result_4487) {
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L5valueS2520 =
          _M0L8_2aentryS992->$5;
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2519;
        moonbit_incref(_M0L5valueS2520);
        _M0L6_2atmpS2519 = _M0L5valueS2520;
        return _M0L6_2atmpS2519;
      } else {
        moonbit_incref(_M0L8_2aentryS992);
      }
      _M0L3pslS2521 = _M0L8_2aentryS992->$2;
      moonbit_decref(_M0L8_2aentryS992);
      if (_M0L1iS987 > _M0L3pslS2521) {
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2522 = 0;
        return _M0L6_2atmpS2522;
      }
      _M0L6_2atmpS2523 = _M0L1iS987 + 1;
      _M0L6_2atmpS2525 = _M0L3idxS988 + 1;
      _M0L14capacity__maskS2526 = _M0L4selfS990->$3;
      _M0L6_2atmpS2524 = _M0L6_2atmpS2525 & _M0L14capacity__maskS2526;
      _M0L1iS987 = _M0L6_2atmpS2523;
      _M0L3idxS988 = _M0L6_2atmpS2524;
      continue;
    }
    break;
  }
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPB3Map3getGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS999,
  moonbit_string_t _M0L3keyS995
) {
  int32_t _M0L4hashS994;
  int32_t _M0L14capacity__maskS2543;
  int32_t _M0L6_2atmpS2542;
  int32_t _M0L1iS996;
  int32_t _M0L3idxS997;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS994 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS995);
  _M0L14capacity__maskS2543 = _M0L4selfS999->$3;
  _M0L6_2atmpS2542 = _M0L4hashS994 & _M0L14capacity__maskS2543;
  _M0L1iS996 = 0;
  _M0L3idxS997 = _M0L6_2atmpS2542;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2541 =
      _M0L4selfS999->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS998;
    if (
      _M0L3idxS997 < 0
      || _M0L3idxS997 >= Moonbit_array_length(_M0L7entriesS2541)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS998
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2541[
        _M0L3idxS997
      ];
    if (_M0L7_2abindS998 == 0) {
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2530 = 0;
      return _M0L6_2atmpS2530;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS1000 =
        _M0L7_2abindS998;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L8_2aentryS1001 =
        _M0L7_2aSomeS1000;
      int32_t _M0L4hashS2532 = _M0L8_2aentryS1001->$3;
      int32_t _if__result_4489;
      int32_t _M0L3pslS2535;
      int32_t _M0L6_2atmpS2537;
      int32_t _M0L6_2atmpS2539;
      int32_t _M0L14capacity__maskS2540;
      int32_t _M0L6_2atmpS2538;
      if (_M0L4hashS2532 == _M0L4hashS994) {
        moonbit_string_t _M0L3keyS2531 = _M0L8_2aentryS1001->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_4489
        = _M0L3keyS2531 == _M0L3keyS995
          || Moonbit_array_length(_M0L3keyS2531)
             == Moonbit_array_length(_M0L3keyS995)
             && 0
                == memcmp(_M0L3keyS2531, _M0L3keyS995, Moonbit_array_length(_M0L3keyS2531) * 2);
      } else {
        _if__result_4489 = 0;
      }
      if (_if__result_4489) {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS2534 =
          _M0L8_2aentryS1001->$5;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2533;
        moonbit_incref(_M0L5valueS2534);
        _M0L6_2atmpS2533 = _M0L5valueS2534;
        return _M0L6_2atmpS2533;
      } else {
        moonbit_incref(_M0L8_2aentryS1001);
      }
      _M0L3pslS2535 = _M0L8_2aentryS1001->$2;
      moonbit_decref(_M0L8_2aentryS1001);
      if (_M0L1iS996 > _M0L3pslS2535) {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2536 = 0;
        return _M0L6_2atmpS2536;
      }
      _M0L6_2atmpS2537 = _M0L1iS996 + 1;
      _M0L6_2atmpS2539 = _M0L3idxS997 + 1;
      _M0L14capacity__maskS2540 = _M0L4selfS999->$3;
      _M0L6_2atmpS2538 = _M0L6_2atmpS2539 & _M0L14capacity__maskS2540;
      _M0L1iS996 = _M0L6_2atmpS2537;
      _M0L3idxS997 = _M0L6_2atmpS2538;
      continue;
    }
    break;
  }
}

struct _M0TPB5ArrayGiE* _M0MPB3Map3getGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS1008,
  moonbit_string_t _M0L3keyS1004
) {
  int32_t _M0L4hashS1003;
  int32_t _M0L14capacity__maskS2557;
  int32_t _M0L6_2atmpS2556;
  int32_t _M0L1iS1005;
  int32_t _M0L3idxS1006;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS1003 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS1004);
  _M0L14capacity__maskS2557 = _M0L4selfS1008->$3;
  _M0L6_2atmpS2556 = _M0L4hashS1003 & _M0L14capacity__maskS2557;
  _M0L1iS1005 = 0;
  _M0L3idxS1006 = _M0L6_2atmpS2556;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2555 =
      _M0L4selfS1008->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS1007;
    if (
      _M0L3idxS1006 < 0
      || _M0L3idxS1006 >= Moonbit_array_length(_M0L7entriesS2555)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS1007
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2555[_M0L3idxS1006];
    if (_M0L7_2abindS1007 == 0) {
      struct _M0TPB5ArrayGiE* _M0L6_2atmpS2544 = 0;
      return _M0L6_2atmpS2544;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS1009 =
        _M0L7_2abindS1007;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L8_2aentryS1010 =
        _M0L7_2aSomeS1009;
      int32_t _M0L4hashS2546 = _M0L8_2aentryS1010->$3;
      int32_t _if__result_4491;
      int32_t _M0L3pslS2549;
      int32_t _M0L6_2atmpS2551;
      int32_t _M0L6_2atmpS2553;
      int32_t _M0L14capacity__maskS2554;
      int32_t _M0L6_2atmpS2552;
      if (_M0L4hashS2546 == _M0L4hashS1003) {
        moonbit_string_t _M0L3keyS2545 = _M0L8_2aentryS1010->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_4491
        = _M0L3keyS2545 == _M0L3keyS1004
          || Moonbit_array_length(_M0L3keyS2545)
             == Moonbit_array_length(_M0L3keyS1004)
             && 0
                == memcmp(_M0L3keyS2545, _M0L3keyS1004, Moonbit_array_length(_M0L3keyS2545) * 2);
      } else {
        _if__result_4491 = 0;
      }
      if (_if__result_4491) {
        struct _M0TPB5ArrayGiE* _M0L5valueS2548 = _M0L8_2aentryS1010->$5;
        struct _M0TPB5ArrayGiE* _M0L6_2atmpS2547;
        moonbit_incref(_M0L5valueS2548);
        _M0L6_2atmpS2547 = _M0L5valueS2548;
        return _M0L6_2atmpS2547;
      } else {
        moonbit_incref(_M0L8_2aentryS1010);
      }
      _M0L3pslS2549 = _M0L8_2aentryS1010->$2;
      moonbit_decref(_M0L8_2aentryS1010);
      if (_M0L1iS1005 > _M0L3pslS2549) {
        struct _M0TPB5ArrayGiE* _M0L6_2atmpS2550 = 0;
        return _M0L6_2atmpS2550;
      }
      _M0L6_2atmpS2551 = _M0L1iS1005 + 1;
      _M0L6_2atmpS2553 = _M0L3idxS1006 + 1;
      _M0L14capacity__maskS2554 = _M0L4selfS1008->$3;
      _M0L6_2atmpS2552 = _M0L6_2atmpS2553 & _M0L14capacity__maskS2554;
      _M0L1iS1005 = _M0L6_2atmpS2551;
      _M0L3idxS1006 = _M0L6_2atmpS2552;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPB3Map3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE _M0L3arrS933,
  int64_t _M0L8capacityS935
) {
  int32_t _M0L3endS2467;
  int32_t _M0L5startS2468;
  int32_t _M0L6lengthS932;
  int32_t _M0L8capacityS934;
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L1mS938;
  int32_t _M0L3endS2464;
  int32_t _M0L5startS2465;
  int32_t _M0L7_2abindS939;
  int32_t _M0L2__S940;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2467 = _M0L3arrS933.$2;
  _M0L5startS2468 = _M0L3arrS933.$1;
  _M0L6lengthS932 = _M0L3endS2467 - _M0L5startS2468;
  if (_M0L8capacityS935 == 4294967296ll) {
    if (_M0L6lengthS932 == 0) {
      _M0L8capacityS934 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS934 = _M0FPB21capacity__for__length(_M0L6lengthS932);
    }
  } else {
    int64_t _M0L7_2aSomeS936 = _M0L8capacityS935;
    int32_t _M0L11_2acapacityS937 = (int32_t)_M0L7_2aSomeS936;
    int32_t _M0L6_2atmpS2466;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2466 = _M0FPB21capacity__for__length(_M0L6lengthS932);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS934
    = _M0MPC13int3Int3max(_M0L11_2acapacityS937, _M0L6_2atmpS2466);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS938
  = _M0FPB8new__mapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L8capacityS934);
  _M0L3endS2464 = _M0L3arrS933.$2;
  _M0L5startS2465 = _M0L3arrS933.$1;
  _M0L7_2abindS939 = _M0L3endS2464 - _M0L5startS2465;
  _M0L2__S940 = 0;
  while (1) {
    if (_M0L2__S940 < _M0L7_2abindS939) {
      struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L3bufS2461 =
        _M0L3arrS933.$0;
      int32_t _M0L5startS2463 = _M0L3arrS933.$1;
      int32_t _M0L6_2atmpS2462 = _M0L5startS2463 + _M0L2__S940;
      struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L1eS941 =
        (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L3bufS2461[
          _M0L6_2atmpS2462
        ];
      moonbit_string_t _M0L6_2atmpS2458 = _M0L1eS941->$0;
      struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2459 =
        _M0L1eS941->$1;
      int32_t _M0L6_2atmpS2460;
      moonbit_incref(_M0L6_2atmpS2459);
      moonbit_incref(_M0L6_2atmpS2458);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L1mS938, _M0L6_2atmpS2458, _M0L6_2atmpS2459);
      moonbit_decref(_M0L6_2atmpS2458);
      moonbit_decref(_M0L6_2atmpS2459);
      _M0L6_2atmpS2460 = _M0L2__S940 + 1;
      _M0L2__S940 = _M0L6_2atmpS2460;
      continue;
    }
    break;
  }
  return _M0L1mS938;
}

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L3arrS944,
  int64_t _M0L8capacityS946
) {
  int32_t _M0L3endS2478;
  int32_t _M0L5startS2479;
  int32_t _M0L6lengthS943;
  int32_t _M0L8capacityS945;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L1mS949;
  int32_t _M0L3endS2475;
  int32_t _M0L5startS2476;
  int32_t _M0L7_2abindS950;
  int32_t _M0L2__S951;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2478 = _M0L3arrS944.$2;
  _M0L5startS2479 = _M0L3arrS944.$1;
  _M0L6lengthS943 = _M0L3endS2478 - _M0L5startS2479;
  if (_M0L8capacityS946 == 4294967296ll) {
    if (_M0L6lengthS943 == 0) {
      _M0L8capacityS945 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS945 = _M0FPB21capacity__for__length(_M0L6lengthS943);
    }
  } else {
    int64_t _M0L7_2aSomeS947 = _M0L8capacityS946;
    int32_t _M0L11_2acapacityS948 = (int32_t)_M0L7_2aSomeS947;
    int32_t _M0L6_2atmpS2477;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2477 = _M0FPB21capacity__for__length(_M0L6lengthS943);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS945
    = _M0MPC13int3Int3max(_M0L11_2acapacityS948, _M0L6_2atmpS2477);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS949 = _M0FPB8new__mapGiUWEuQRPC15error5ErrorNsEE(_M0L8capacityS945);
  _M0L3endS2475 = _M0L3arrS944.$2;
  _M0L5startS2476 = _M0L3arrS944.$1;
  _M0L7_2abindS950 = _M0L3endS2475 - _M0L5startS2476;
  _M0L2__S951 = 0;
  while (1) {
    if (_M0L2__S951 < _M0L7_2abindS950) {
      struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L3bufS2472 =
        _M0L3arrS944.$0;
      int32_t _M0L5startS2474 = _M0L3arrS944.$1;
      int32_t _M0L6_2atmpS2473 = _M0L5startS2474 + _M0L2__S951;
      struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L1eS952 =
        (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)_M0L3bufS2472[
          _M0L6_2atmpS2473
        ];
      int32_t _M0L6_2atmpS2469 = _M0L1eS952->$0;
      struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2470 = _M0L1eS952->$1;
      int32_t _M0L6_2atmpS2471;
      moonbit_incref(_M0L6_2atmpS2470);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGiUWEuQRPC15error5ErrorNsEE(_M0L1mS949, _M0L6_2atmpS2469, _M0L6_2atmpS2470);
      moonbit_decref(_M0L6_2atmpS2470);
      _M0L6_2atmpS2471 = _M0L2__S951 + 1;
      _M0L2__S951 = _M0L6_2atmpS2471;
      continue;
    }
    break;
  }
  return _M0L1mS949;
}

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0MPB3Map3MapGsRPB5ArrayGiEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L3arrS955,
  int64_t _M0L8capacityS957
) {
  int32_t _M0L3endS2489;
  int32_t _M0L5startS2490;
  int32_t _M0L6lengthS954;
  int32_t _M0L8capacityS956;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L1mS960;
  int32_t _M0L3endS2486;
  int32_t _M0L5startS2487;
  int32_t _M0L7_2abindS961;
  int32_t _M0L2__S962;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2489 = _M0L3arrS955.$2;
  _M0L5startS2490 = _M0L3arrS955.$1;
  _M0L6lengthS954 = _M0L3endS2489 - _M0L5startS2490;
  if (_M0L8capacityS957 == 4294967296ll) {
    if (_M0L6lengthS954 == 0) {
      _M0L8capacityS956 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS956 = _M0FPB21capacity__for__length(_M0L6lengthS954);
    }
  } else {
    int64_t _M0L7_2aSomeS958 = _M0L8capacityS957;
    int32_t _M0L11_2acapacityS959 = (int32_t)_M0L7_2aSomeS958;
    int32_t _M0L6_2atmpS2488;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2488 = _M0FPB21capacity__for__length(_M0L6lengthS954);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS956
    = _M0MPC13int3Int3max(_M0L11_2acapacityS959, _M0L6_2atmpS2488);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS960 = _M0FPB8new__mapGsRPB5ArrayGiEE(_M0L8capacityS956);
  _M0L3endS2486 = _M0L3arrS955.$2;
  _M0L5startS2487 = _M0L3arrS955.$1;
  _M0L7_2abindS961 = _M0L3endS2486 - _M0L5startS2487;
  _M0L2__S962 = 0;
  while (1) {
    if (_M0L2__S962 < _M0L7_2abindS961) {
      struct _M0TUsRPB5ArrayGiEE** _M0L3bufS2483 = _M0L3arrS955.$0;
      int32_t _M0L5startS2485 = _M0L3arrS955.$1;
      int32_t _M0L6_2atmpS2484 = _M0L5startS2485 + _M0L2__S962;
      struct _M0TUsRPB5ArrayGiEE* _M0L1eS963 =
        (struct _M0TUsRPB5ArrayGiEE*)_M0L3bufS2483[_M0L6_2atmpS2484];
      moonbit_string_t _M0L6_2atmpS2480 = _M0L1eS963->$0;
      struct _M0TPB5ArrayGiE* _M0L6_2atmpS2481 = _M0L1eS963->$1;
      int32_t _M0L6_2atmpS2482;
      moonbit_incref(_M0L6_2atmpS2481);
      moonbit_incref(_M0L6_2atmpS2480);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L1mS960, _M0L6_2atmpS2480, _M0L6_2atmpS2481);
      moonbit_decref(_M0L6_2atmpS2480);
      moonbit_decref(_M0L6_2atmpS2481);
      _M0L6_2atmpS2482 = _M0L2__S962 + 1;
      _M0L2__S962 = _M0L6_2atmpS2482;
      continue;
    }
    break;
  }
  return _M0L1mS960;
}

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0MPB3Map3MapGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB9ArrayViewGUsRP26biolab8bio__seq9SeqRecordEE _M0L3arrS966,
  int64_t _M0L8capacityS968
) {
  int32_t _M0L3endS2500;
  int32_t _M0L5startS2501;
  int32_t _M0L6lengthS965;
  int32_t _M0L8capacityS967;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L1mS971;
  int32_t _M0L3endS2497;
  int32_t _M0L5startS2498;
  int32_t _M0L7_2abindS972;
  int32_t _M0L2__S973;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2500 = _M0L3arrS966.$2;
  _M0L5startS2501 = _M0L3arrS966.$1;
  _M0L6lengthS965 = _M0L3endS2500 - _M0L5startS2501;
  if (_M0L8capacityS968 == 4294967296ll) {
    if (_M0L6lengthS965 == 0) {
      _M0L8capacityS967 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS967 = _M0FPB21capacity__for__length(_M0L6lengthS965);
    }
  } else {
    int64_t _M0L7_2aSomeS969 = _M0L8capacityS968;
    int32_t _M0L11_2acapacityS970 = (int32_t)_M0L7_2aSomeS969;
    int32_t _M0L6_2atmpS2499;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2499 = _M0FPB21capacity__for__length(_M0L6lengthS965);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS967
    = _M0MPC13int3Int3max(_M0L11_2acapacityS970, _M0L6_2atmpS2499);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS971
  = _M0FPB8new__mapGsRP26biolab8bio__seq9SeqRecordE(_M0L8capacityS967);
  _M0L3endS2497 = _M0L3arrS966.$2;
  _M0L5startS2498 = _M0L3arrS966.$1;
  _M0L7_2abindS972 = _M0L3endS2497 - _M0L5startS2498;
  _M0L2__S973 = 0;
  while (1) {
    if (_M0L2__S973 < _M0L7_2abindS972) {
      struct _M0TUsRP26biolab8bio__seq9SeqRecordE** _M0L3bufS2494 =
        _M0L3arrS966.$0;
      int32_t _M0L5startS2496 = _M0L3arrS966.$1;
      int32_t _M0L6_2atmpS2495 = _M0L5startS2496 + _M0L2__S973;
      struct _M0TUsRP26biolab8bio__seq9SeqRecordE* _M0L1eS974 =
        (struct _M0TUsRP26biolab8bio__seq9SeqRecordE*)_M0L3bufS2494[
          _M0L6_2atmpS2495
        ];
      moonbit_string_t _M0L6_2atmpS2491 = _M0L1eS974->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2492 =
        _M0L1eS974->$1;
      int32_t _M0L6_2atmpS2493;
      moonbit_incref(_M0L6_2atmpS2492);
      moonbit_incref(_M0L6_2atmpS2491);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(_M0L1mS971, _M0L6_2atmpS2491, _M0L6_2atmpS2492);
      moonbit_decref(_M0L6_2atmpS2491);
      moonbit_decref(_M0L6_2atmpS2492);
      _M0L6_2atmpS2493 = _M0L2__S973 + 1;
      _M0L2__S973 = _M0L6_2atmpS2493;
      continue;
    }
    break;
  }
  return _M0L1mS971;
}

int32_t _M0MPB3Map3setGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS920,
  moonbit_string_t _M0L3keyS921,
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L5valueS922
) {
  int32_t _M0L6_2atmpS2454;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2454 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS921);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS920, _M0L3keyS921, _M0L5valueS922, _M0L6_2atmpS2454);
  return 0;
}

int32_t _M0MPB3Map3setGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS923,
  int32_t _M0L3keyS924,
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L5valueS925
) {
  int32_t _M0L6_2atmpS2455;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2455 = _M0IPC13int3IntPB4Hash4hash(_M0L3keyS924);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS923, _M0L3keyS924, _M0L5valueS925, _M0L6_2atmpS2455);
  return 0;
}

int32_t _M0MPB3Map3setGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS926,
  moonbit_string_t _M0L3keyS927,
  struct _M0TPB5ArrayGiE* _M0L5valueS928
) {
  int32_t _M0L6_2atmpS2456;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2456 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS927);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(_M0L4selfS926, _M0L3keyS927, _M0L5valueS928, _M0L6_2atmpS2456);
  return 0;
}

int32_t _M0MPB3Map3setGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS929,
  moonbit_string_t _M0L3keyS930,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS931
) {
  int32_t _M0L6_2atmpS2457;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2457 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS930);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS929, _M0L3keyS930, _M0L5valueS931, _M0L6_2atmpS2457);
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS859,
  moonbit_string_t _M0L3keyS865,
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L5valueS866,
  int32_t _M0L4hashS861
) {
  int32_t _M0L14capacity__maskS2399;
  int32_t _M0L6_2atmpS2398;
  int32_t _M0L3pslS856;
  int32_t _M0L3idxS857;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2399 = _M0L4selfS859->$3;
  _M0L6_2atmpS2398 = _M0L4hashS861 & _M0L14capacity__maskS2399;
  _M0L3pslS856 = 0;
  _M0L3idxS857 = _M0L6_2atmpS2398;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2397 =
      _M0L4selfS859->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS858;
    if (
      _M0L3idxS857 < 0
      || _M0L3idxS857 >= Moonbit_array_length(_M0L7entriesS2397)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS858
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2397[
        _M0L3idxS857
      ];
    if (_M0L7_2abindS858 == 0) {
      int32_t _M0L4sizeS2382 = _M0L4selfS859->$1;
      int32_t _M0L8grow__atS2383 = _M0L4selfS859->$4;
      int32_t _M0L7_2abindS862;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS863;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS864;
      if (_M0L4sizeS2382 >= _M0L8grow__atS2383) {
        int32_t _M0L14capacity__maskS2385;
        int32_t _M0L6_2atmpS2384;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS859);
        _M0L14capacity__maskS2385 = _M0L4selfS859->$3;
        _M0L6_2atmpS2384 = _M0L4hashS861 & _M0L14capacity__maskS2385;
        _M0L3pslS856 = 0;
        _M0L3idxS857 = _M0L6_2atmpS2384;
        continue;
      }
      _M0L7_2abindS862 = _M0L4selfS859->$6;
      _M0L7_2abindS863 = 0;
      moonbit_incref(_M0L3keyS865);
      moonbit_incref(_M0L5valueS866);
      _M0L5entryS864
      = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
      Moonbit_object_header(_M0L5entryS864)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 69, 0);
      _M0L5entryS864->$0 = _M0L7_2abindS862;
      _M0L5entryS864->$1 = _M0L7_2abindS863;
      _M0L5entryS864->$2 = _M0L3pslS856;
      _M0L5entryS864->$3 = _M0L4hashS861;
      _M0L5entryS864->$4 = _M0L3keyS865;
      _M0L5entryS864->$5 = _M0L5valueS866;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS859, _M0L3idxS857, _M0L5entryS864);
      moonbit_decref(_M0L5entryS864);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS867 =
        _M0L7_2abindS858;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L14_2acurr__entryS868 =
        _M0L7_2aSomeS867;
      int32_t _M0L4hashS2387 = _M0L14_2acurr__entryS868->$3;
      int32_t _if__result_4497;
      int32_t _M0L3pslS2388;
      int32_t _M0L6_2atmpS2393;
      int32_t _M0L6_2atmpS2395;
      int32_t _M0L14capacity__maskS2396;
      int32_t _M0L6_2atmpS2394;
      if (_M0L4hashS2387 == _M0L4hashS861) {
        moonbit_string_t _M0L3keyS2386 = _M0L14_2acurr__entryS868->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_4497
        = _M0L3keyS2386 == _M0L3keyS865
          || Moonbit_array_length(_M0L3keyS2386)
             == Moonbit_array_length(_M0L3keyS865)
             && 0
                == memcmp(_M0L3keyS2386, _M0L3keyS865, Moonbit_array_length(_M0L3keyS2386) * 2);
      } else {
        _if__result_4497 = 0;
      }
      if (_if__result_4497) {
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3747 =
          _M0L14_2acurr__entryS868->$5;
        moonbit_incref(_M0L5valueS866);
        moonbit_decref(_M0L6_2aoldS3747);
        _M0L14_2acurr__entryS868->$5 = _M0L5valueS866;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS868);
      }
      _M0L3pslS2388 = _M0L14_2acurr__entryS868->$2;
      if (_M0L3pslS856 > _M0L3pslS2388) {
        int32_t _M0L4sizeS2389 = _M0L4selfS859->$1;
        int32_t _M0L8grow__atS2390 = _M0L4selfS859->$4;
        int32_t _M0L7_2abindS869;
        struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS870;
        struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS871;
        if (_M0L4sizeS2389 >= _M0L8grow__atS2390) {
          int32_t _M0L14capacity__maskS2392;
          int32_t _M0L6_2atmpS2391;
          moonbit_decref(_M0L14_2acurr__entryS868);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS859);
          _M0L14capacity__maskS2392 = _M0L4selfS859->$3;
          _M0L6_2atmpS2391 = _M0L4hashS861 & _M0L14capacity__maskS2392;
          _M0L3pslS856 = 0;
          _M0L3idxS857 = _M0L6_2atmpS2391;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS859, _M0L3idxS857, _M0L14_2acurr__entryS868);
        moonbit_decref(_M0L14_2acurr__entryS868);
        _M0L7_2abindS869 = _M0L4selfS859->$6;
        _M0L7_2abindS870 = 0;
        moonbit_incref(_M0L3keyS865);
        moonbit_incref(_M0L5valueS866);
        _M0L5entryS871
        = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
        Moonbit_object_header(_M0L5entryS871)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 69, 0);
        _M0L5entryS871->$0 = _M0L7_2abindS869;
        _M0L5entryS871->$1 = _M0L7_2abindS870;
        _M0L5entryS871->$2 = _M0L3pslS856;
        _M0L5entryS871->$3 = _M0L4hashS861;
        _M0L5entryS871->$4 = _M0L3keyS865;
        _M0L5entryS871->$5 = _M0L5valueS866;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS859, _M0L3idxS857, _M0L5entryS871);
        moonbit_decref(_M0L5entryS871);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS868);
      }
      _M0L6_2atmpS2393 = _M0L3pslS856 + 1;
      _M0L6_2atmpS2395 = _M0L3idxS857 + 1;
      _M0L14capacity__maskS2396 = _M0L4selfS859->$3;
      _M0L6_2atmpS2394 = _M0L6_2atmpS2395 & _M0L14capacity__maskS2396;
      _M0L3pslS856 = _M0L6_2atmpS2393;
      _M0L3idxS857 = _M0L6_2atmpS2394;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS875,
  int32_t _M0L3keyS881,
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L5valueS882,
  int32_t _M0L4hashS877
) {
  int32_t _M0L14capacity__maskS2417;
  int32_t _M0L6_2atmpS2416;
  int32_t _M0L3pslS872;
  int32_t _M0L3idxS873;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2417 = _M0L4selfS875->$3;
  _M0L6_2atmpS2416 = _M0L4hashS877 & _M0L14capacity__maskS2417;
  _M0L3pslS872 = 0;
  _M0L3idxS873 = _M0L6_2atmpS2416;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2415 =
      _M0L4selfS875->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS874;
    if (
      _M0L3idxS873 < 0
      || _M0L3idxS873 >= Moonbit_array_length(_M0L7entriesS2415)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS874
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2415[
        _M0L3idxS873
      ];
    if (_M0L7_2abindS874 == 0) {
      int32_t _M0L4sizeS2400 = _M0L4selfS875->$1;
      int32_t _M0L8grow__atS2401 = _M0L4selfS875->$4;
      int32_t _M0L7_2abindS878;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS879;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS880;
      if (_M0L4sizeS2400 >= _M0L8grow__atS2401) {
        int32_t _M0L14capacity__maskS2403;
        int32_t _M0L6_2atmpS2402;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS875);
        _M0L14capacity__maskS2403 = _M0L4selfS875->$3;
        _M0L6_2atmpS2402 = _M0L4hashS877 & _M0L14capacity__maskS2403;
        _M0L3pslS872 = 0;
        _M0L3idxS873 = _M0L6_2atmpS2402;
        continue;
      }
      _M0L7_2abindS878 = _M0L4selfS875->$6;
      _M0L7_2abindS879 = 0;
      moonbit_incref(_M0L5valueS882);
      _M0L5entryS880
      = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE));
      Moonbit_object_header(_M0L5entryS880)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 74, 0);
      _M0L5entryS880->$0 = _M0L7_2abindS878;
      _M0L5entryS880->$1 = _M0L7_2abindS879;
      _M0L5entryS880->$2 = _M0L3pslS872;
      _M0L5entryS880->$3 = _M0L4hashS877;
      _M0L5entryS880->$4 = _M0L3keyS881;
      _M0L5entryS880->$5 = _M0L5valueS882;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS875, _M0L3idxS873, _M0L5entryS880);
      moonbit_decref(_M0L5entryS880);
      return 0;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS883 =
        _M0L7_2abindS874;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L14_2acurr__entryS884 =
        _M0L7_2aSomeS883;
      int32_t _M0L4hashS2405 = _M0L14_2acurr__entryS884->$3;
      int32_t _if__result_4499;
      int32_t _M0L3pslS2406;
      int32_t _M0L6_2atmpS2411;
      int32_t _M0L6_2atmpS2413;
      int32_t _M0L14capacity__maskS2414;
      int32_t _M0L6_2atmpS2412;
      if (_M0L4hashS2405 == _M0L4hashS877) {
        int32_t _M0L3keyS2404 = _M0L14_2acurr__entryS884->$4;
        _if__result_4499 = _M0L3keyS2404 == _M0L3keyS881;
      } else {
        _if__result_4499 = 0;
      }
      if (_if__result_4499) {
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2aoldS3751 =
          _M0L14_2acurr__entryS884->$5;
        moonbit_incref(_M0L5valueS882);
        moonbit_decref(_M0L6_2aoldS3751);
        _M0L14_2acurr__entryS884->$5 = _M0L5valueS882;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS884);
      }
      _M0L3pslS2406 = _M0L14_2acurr__entryS884->$2;
      if (_M0L3pslS872 > _M0L3pslS2406) {
        int32_t _M0L4sizeS2407 = _M0L4selfS875->$1;
        int32_t _M0L8grow__atS2408 = _M0L4selfS875->$4;
        int32_t _M0L7_2abindS885;
        struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS886;
        struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS887;
        if (_M0L4sizeS2407 >= _M0L8grow__atS2408) {
          int32_t _M0L14capacity__maskS2410;
          int32_t _M0L6_2atmpS2409;
          moonbit_decref(_M0L14_2acurr__entryS884);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS875);
          _M0L14capacity__maskS2410 = _M0L4selfS875->$3;
          _M0L6_2atmpS2409 = _M0L4hashS877 & _M0L14capacity__maskS2410;
          _M0L3pslS872 = 0;
          _M0L3idxS873 = _M0L6_2atmpS2409;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS875, _M0L3idxS873, _M0L14_2acurr__entryS884);
        moonbit_decref(_M0L14_2acurr__entryS884);
        _M0L7_2abindS885 = _M0L4selfS875->$6;
        _M0L7_2abindS886 = 0;
        moonbit_incref(_M0L5valueS882);
        _M0L5entryS887
        = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE));
        Moonbit_object_header(_M0L5entryS887)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 74, 0);
        _M0L5entryS887->$0 = _M0L7_2abindS885;
        _M0L5entryS887->$1 = _M0L7_2abindS886;
        _M0L5entryS887->$2 = _M0L3pslS872;
        _M0L5entryS887->$3 = _M0L4hashS877;
        _M0L5entryS887->$4 = _M0L3keyS881;
        _M0L5entryS887->$5 = _M0L5valueS882;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS875, _M0L3idxS873, _M0L5entryS887);
        moonbit_decref(_M0L5entryS887);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS884);
      }
      _M0L6_2atmpS2411 = _M0L3pslS872 + 1;
      _M0L6_2atmpS2413 = _M0L3idxS873 + 1;
      _M0L14capacity__maskS2414 = _M0L4selfS875->$3;
      _M0L6_2atmpS2412 = _M0L6_2atmpS2413 & _M0L14capacity__maskS2414;
      _M0L3pslS872 = _M0L6_2atmpS2411;
      _M0L3idxS873 = _M0L6_2atmpS2412;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS891,
  moonbit_string_t _M0L3keyS897,
  struct _M0TPB5ArrayGiE* _M0L5valueS898,
  int32_t _M0L4hashS893
) {
  int32_t _M0L14capacity__maskS2435;
  int32_t _M0L6_2atmpS2434;
  int32_t _M0L3pslS888;
  int32_t _M0L3idxS889;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2435 = _M0L4selfS891->$3;
  _M0L6_2atmpS2434 = _M0L4hashS893 & _M0L14capacity__maskS2435;
  _M0L3pslS888 = 0;
  _M0L3idxS889 = _M0L6_2atmpS2434;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2433 =
      _M0L4selfS891->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS890;
    if (
      _M0L3idxS889 < 0
      || _M0L3idxS889 >= Moonbit_array_length(_M0L7entriesS2433)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS890
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2433[_M0L3idxS889];
    if (_M0L7_2abindS890 == 0) {
      int32_t _M0L4sizeS2418 = _M0L4selfS891->$1;
      int32_t _M0L8grow__atS2419 = _M0L4selfS891->$4;
      int32_t _M0L7_2abindS894;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS895;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS896;
      if (_M0L4sizeS2418 >= _M0L8grow__atS2419) {
        int32_t _M0L14capacity__maskS2421;
        int32_t _M0L6_2atmpS2420;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB5ArrayGiEE(_M0L4selfS891);
        _M0L14capacity__maskS2421 = _M0L4selfS891->$3;
        _M0L6_2atmpS2420 = _M0L4hashS893 & _M0L14capacity__maskS2421;
        _M0L3pslS888 = 0;
        _M0L3idxS889 = _M0L6_2atmpS2420;
        continue;
      }
      _M0L7_2abindS894 = _M0L4selfS891->$6;
      _M0L7_2abindS895 = 0;
      moonbit_incref(_M0L3keyS897);
      moonbit_incref(_M0L5valueS898);
      _M0L5entryS896
      = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE));
      Moonbit_object_header(_M0L5entryS896)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 78, 0);
      _M0L5entryS896->$0 = _M0L7_2abindS894;
      _M0L5entryS896->$1 = _M0L7_2abindS895;
      _M0L5entryS896->$2 = _M0L3pslS888;
      _M0L5entryS896->$3 = _M0L4hashS893;
      _M0L5entryS896->$4 = _M0L3keyS897;
      _M0L5entryS896->$5 = _M0L5valueS898;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS891, _M0L3idxS889, _M0L5entryS896);
      moonbit_decref(_M0L5entryS896);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS899 = _M0L7_2abindS890;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L14_2acurr__entryS900 =
        _M0L7_2aSomeS899;
      int32_t _M0L4hashS2423 = _M0L14_2acurr__entryS900->$3;
      int32_t _if__result_4501;
      int32_t _M0L3pslS2424;
      int32_t _M0L6_2atmpS2429;
      int32_t _M0L6_2atmpS2431;
      int32_t _M0L14capacity__maskS2432;
      int32_t _M0L6_2atmpS2430;
      if (_M0L4hashS2423 == _M0L4hashS893) {
        moonbit_string_t _M0L3keyS2422 = _M0L14_2acurr__entryS900->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_4501
        = _M0L3keyS2422 == _M0L3keyS897
          || Moonbit_array_length(_M0L3keyS2422)
             == Moonbit_array_length(_M0L3keyS897)
             && 0
                == memcmp(_M0L3keyS2422, _M0L3keyS897, Moonbit_array_length(_M0L3keyS2422) * 2);
      } else {
        _if__result_4501 = 0;
      }
      if (_if__result_4501) {
        struct _M0TPB5ArrayGiE* _M0L6_2aoldS3754 =
          _M0L14_2acurr__entryS900->$5;
        moonbit_incref(_M0L5valueS898);
        moonbit_decref(_M0L6_2aoldS3754);
        _M0L14_2acurr__entryS900->$5 = _M0L5valueS898;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS900);
      }
      _M0L3pslS2424 = _M0L14_2acurr__entryS900->$2;
      if (_M0L3pslS888 > _M0L3pslS2424) {
        int32_t _M0L4sizeS2425 = _M0L4selfS891->$1;
        int32_t _M0L8grow__atS2426 = _M0L4selfS891->$4;
        int32_t _M0L7_2abindS901;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS902;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS903;
        if (_M0L4sizeS2425 >= _M0L8grow__atS2426) {
          int32_t _M0L14capacity__maskS2428;
          int32_t _M0L6_2atmpS2427;
          moonbit_decref(_M0L14_2acurr__entryS900);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB5ArrayGiEE(_M0L4selfS891);
          _M0L14capacity__maskS2428 = _M0L4selfS891->$3;
          _M0L6_2atmpS2427 = _M0L4hashS893 & _M0L14capacity__maskS2428;
          _M0L3pslS888 = 0;
          _M0L3idxS889 = _M0L6_2atmpS2427;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGiEE(_M0L4selfS891, _M0L3idxS889, _M0L14_2acurr__entryS900);
        moonbit_decref(_M0L14_2acurr__entryS900);
        _M0L7_2abindS901 = _M0L4selfS891->$6;
        _M0L7_2abindS902 = 0;
        moonbit_incref(_M0L3keyS897);
        moonbit_incref(_M0L5valueS898);
        _M0L5entryS903
        = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE));
        Moonbit_object_header(_M0L5entryS903)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 78, 0);
        _M0L5entryS903->$0 = _M0L7_2abindS901;
        _M0L5entryS903->$1 = _M0L7_2abindS902;
        _M0L5entryS903->$2 = _M0L3pslS888;
        _M0L5entryS903->$3 = _M0L4hashS893;
        _M0L5entryS903->$4 = _M0L3keyS897;
        _M0L5entryS903->$5 = _M0L5valueS898;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS891, _M0L3idxS889, _M0L5entryS903);
        moonbit_decref(_M0L5entryS903);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS900);
      }
      _M0L6_2atmpS2429 = _M0L3pslS888 + 1;
      _M0L6_2atmpS2431 = _M0L3idxS889 + 1;
      _M0L14capacity__maskS2432 = _M0L4selfS891->$3;
      _M0L6_2atmpS2430 = _M0L6_2atmpS2431 & _M0L14capacity__maskS2432;
      _M0L3pslS888 = _M0L6_2atmpS2429;
      _M0L3idxS889 = _M0L6_2atmpS2430;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS907,
  moonbit_string_t _M0L3keyS913,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS914,
  int32_t _M0L4hashS909
) {
  int32_t _M0L14capacity__maskS2453;
  int32_t _M0L6_2atmpS2452;
  int32_t _M0L3pslS904;
  int32_t _M0L3idxS905;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2453 = _M0L4selfS907->$3;
  _M0L6_2atmpS2452 = _M0L4hashS909 & _M0L14capacity__maskS2453;
  _M0L3pslS904 = 0;
  _M0L3idxS905 = _M0L6_2atmpS2452;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2451 =
      _M0L4selfS907->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS906;
    if (
      _M0L3idxS905 < 0
      || _M0L3idxS905 >= Moonbit_array_length(_M0L7entriesS2451)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS906
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2451[
        _M0L3idxS905
      ];
    if (_M0L7_2abindS906 == 0) {
      int32_t _M0L4sizeS2436 = _M0L4selfS907->$1;
      int32_t _M0L8grow__atS2437 = _M0L4selfS907->$4;
      int32_t _M0L7_2abindS910;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS911;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS912;
      if (_M0L4sizeS2436 >= _M0L8grow__atS2437) {
        int32_t _M0L14capacity__maskS2439;
        int32_t _M0L6_2atmpS2438;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS907);
        _M0L14capacity__maskS2439 = _M0L4selfS907->$3;
        _M0L6_2atmpS2438 = _M0L4hashS909 & _M0L14capacity__maskS2439;
        _M0L3pslS904 = 0;
        _M0L3idxS905 = _M0L6_2atmpS2438;
        continue;
      }
      _M0L7_2abindS910 = _M0L4selfS907->$6;
      _M0L7_2abindS911 = 0;
      moonbit_incref(_M0L3keyS913);
      moonbit_incref(_M0L5valueS914);
      _M0L5entryS912
      = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE));
      Moonbit_object_header(_M0L5entryS912)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 83, 0);
      _M0L5entryS912->$0 = _M0L7_2abindS910;
      _M0L5entryS912->$1 = _M0L7_2abindS911;
      _M0L5entryS912->$2 = _M0L3pslS904;
      _M0L5entryS912->$3 = _M0L4hashS909;
      _M0L5entryS912->$4 = _M0L3keyS913;
      _M0L5entryS912->$5 = _M0L5valueS914;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS907, _M0L3idxS905, _M0L5entryS912);
      moonbit_decref(_M0L5entryS912);
      return 0;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS915 =
        _M0L7_2abindS906;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L14_2acurr__entryS916 =
        _M0L7_2aSomeS915;
      int32_t _M0L4hashS2441 = _M0L14_2acurr__entryS916->$3;
      int32_t _if__result_4503;
      int32_t _M0L3pslS2442;
      int32_t _M0L6_2atmpS2447;
      int32_t _M0L6_2atmpS2449;
      int32_t _M0L14capacity__maskS2450;
      int32_t _M0L6_2atmpS2448;
      if (_M0L4hashS2441 == _M0L4hashS909) {
        moonbit_string_t _M0L3keyS2440 = _M0L14_2acurr__entryS916->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_4503
        = _M0L3keyS2440 == _M0L3keyS913
          || Moonbit_array_length(_M0L3keyS2440)
             == Moonbit_array_length(_M0L3keyS913)
             && 0
                == memcmp(_M0L3keyS2440, _M0L3keyS913, Moonbit_array_length(_M0L3keyS2440) * 2);
      } else {
        _if__result_4503 = 0;
      }
      if (_if__result_4503) {
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3758 =
          _M0L14_2acurr__entryS916->$5;
        moonbit_incref(_M0L5valueS914);
        moonbit_decref(_M0L6_2aoldS3758);
        _M0L14_2acurr__entryS916->$5 = _M0L5valueS914;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS916);
      }
      _M0L3pslS2442 = _M0L14_2acurr__entryS916->$2;
      if (_M0L3pslS904 > _M0L3pslS2442) {
        int32_t _M0L4sizeS2443 = _M0L4selfS907->$1;
        int32_t _M0L8grow__atS2444 = _M0L4selfS907->$4;
        int32_t _M0L7_2abindS917;
        struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS918;
        struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS919;
        if (_M0L4sizeS2443 >= _M0L8grow__atS2444) {
          int32_t _M0L14capacity__maskS2446;
          int32_t _M0L6_2atmpS2445;
          moonbit_decref(_M0L14_2acurr__entryS916);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS907);
          _M0L14capacity__maskS2446 = _M0L4selfS907->$3;
          _M0L6_2atmpS2445 = _M0L4hashS909 & _M0L14capacity__maskS2446;
          _M0L3pslS904 = 0;
          _M0L3idxS905 = _M0L6_2atmpS2445;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS907, _M0L3idxS905, _M0L14_2acurr__entryS916);
        moonbit_decref(_M0L14_2acurr__entryS916);
        _M0L7_2abindS917 = _M0L4selfS907->$6;
        _M0L7_2abindS918 = 0;
        moonbit_incref(_M0L3keyS913);
        moonbit_incref(_M0L5valueS914);
        _M0L5entryS919
        = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE));
        Moonbit_object_header(_M0L5entryS919)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 83, 0);
        _M0L5entryS919->$0 = _M0L7_2abindS917;
        _M0L5entryS919->$1 = _M0L7_2abindS918;
        _M0L5entryS919->$2 = _M0L3pslS904;
        _M0L5entryS919->$3 = _M0L4hashS909;
        _M0L5entryS919->$4 = _M0L3keyS913;
        _M0L5entryS919->$5 = _M0L5valueS914;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS907, _M0L3idxS905, _M0L5entryS919);
        moonbit_decref(_M0L5entryS919);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS916);
      }
      _M0L6_2atmpS2447 = _M0L3pslS904 + 1;
      _M0L6_2atmpS2449 = _M0L3idxS905 + 1;
      _M0L14capacity__maskS2450 = _M0L4selfS907->$3;
      _M0L6_2atmpS2448 = _M0L6_2atmpS2449 & _M0L14capacity__maskS2450;
      _M0L3pslS904 = _M0L6_2atmpS2447;
      _M0L3idxS905 = _M0L6_2atmpS2448;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS825
) {
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L9old__headS824;
  int32_t _M0L8capacityS2357;
  int32_t _M0L13new__capacityS826;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2351;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L6_2atmpS2350;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L6_2aoldS3765;
  int32_t _M0L6_2atmpS2352;
  int32_t _M0L8capacityS2354;
  int32_t _M0L6_2atmpS2353;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2355;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3764;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L1xS827;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS824 = _M0L4selfS825->$5;
  _M0L8capacityS2357 = _M0L4selfS825->$2;
  _M0L13new__capacityS826 = _M0L8capacityS2357 << 1;
  _M0L6_2atmpS2351 = 0;
  _M0L6_2atmpS2350
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE**)moonbit_make_ref_array(_M0L13new__capacityS826, _M0L6_2atmpS2351);
  _M0L6_2aoldS3765 = _M0L4selfS825->$0;
  if (_M0L9old__headS824) {
    moonbit_incref(_M0L9old__headS824);
  }
  moonbit_decref(_M0L6_2aoldS3765);
  _M0L4selfS825->$0 = _M0L6_2atmpS2350;
  _M0L4selfS825->$2 = _M0L13new__capacityS826;
  _M0L6_2atmpS2352 = _M0L13new__capacityS826 - 1;
  _M0L4selfS825->$3 = _M0L6_2atmpS2352;
  _M0L8capacityS2354 = _M0L4selfS825->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2353 = _M0FPB21calc__grow__threshold(_M0L8capacityS2354);
  _M0L4selfS825->$4 = _M0L6_2atmpS2353;
  _M0L4selfS825->$1 = 0;
  _M0L6_2atmpS2355 = 0;
  _M0L6_2aoldS3764 = _M0L4selfS825->$5;
  if (_M0L6_2aoldS3764) {
    moonbit_decref(_M0L6_2aoldS3764);
  }
  _M0L4selfS825->$5 = _M0L6_2atmpS2355;
  _M0L4selfS825->$6 = -1;
  _M0L1xS827 = _M0L9old__headS824;
  while (1) {
    if (_M0L1xS827 == 0) {
      if (_M0L1xS827) {
        moonbit_decref(_M0L1xS827);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS829 =
        _M0L1xS827;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4_2aeS830 =
        _M0L7_2aSomeS829;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L15next__in__chainS831 =
        _M0L4_2aeS830->$1;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2356 =
        0;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3762 =
        _M0L4_2aeS830->$1;
      if (_M0L15next__in__chainS831) {
        moonbit_incref(_M0L15next__in__chainS831);
      }
      if (_M0L6_2aoldS3762) {
        moonbit_decref(_M0L6_2aoldS3762);
      }
      _M0L4_2aeS830->$1 = _M0L6_2atmpS2356;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS825, _M0L4_2aeS830);
      moonbit_decref(_M0L4_2aeS830);
      _M0L1xS827 = _M0L15next__in__chainS831;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS833
) {
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L9old__headS832;
  int32_t _M0L8capacityS2365;
  int32_t _M0L13new__capacityS834;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2359;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS2358;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L6_2aoldS3770;
  int32_t _M0L6_2atmpS2360;
  int32_t _M0L8capacityS2362;
  int32_t _M0L6_2atmpS2361;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2363;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3769;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L1xS835;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS832 = _M0L4selfS833->$5;
  _M0L8capacityS2365 = _M0L4selfS833->$2;
  _M0L13new__capacityS834 = _M0L8capacityS2365 << 1;
  _M0L6_2atmpS2359 = 0;
  _M0L6_2atmpS2358
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE**)moonbit_make_ref_array(_M0L13new__capacityS834, _M0L6_2atmpS2359);
  _M0L6_2aoldS3770 = _M0L4selfS833->$0;
  if (_M0L9old__headS832) {
    moonbit_incref(_M0L9old__headS832);
  }
  moonbit_decref(_M0L6_2aoldS3770);
  _M0L4selfS833->$0 = _M0L6_2atmpS2358;
  _M0L4selfS833->$2 = _M0L13new__capacityS834;
  _M0L6_2atmpS2360 = _M0L13new__capacityS834 - 1;
  _M0L4selfS833->$3 = _M0L6_2atmpS2360;
  _M0L8capacityS2362 = _M0L4selfS833->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2361 = _M0FPB21calc__grow__threshold(_M0L8capacityS2362);
  _M0L4selfS833->$4 = _M0L6_2atmpS2361;
  _M0L4selfS833->$1 = 0;
  _M0L6_2atmpS2363 = 0;
  _M0L6_2aoldS3769 = _M0L4selfS833->$5;
  if (_M0L6_2aoldS3769) {
    moonbit_decref(_M0L6_2aoldS3769);
  }
  _M0L4selfS833->$5 = _M0L6_2atmpS2363;
  _M0L4selfS833->$6 = -1;
  _M0L1xS835 = _M0L9old__headS832;
  while (1) {
    if (_M0L1xS835 == 0) {
      if (_M0L1xS835) {
        moonbit_decref(_M0L1xS835);
      }
      break;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS837 =
        _M0L1xS835;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L4_2aeS838 =
        _M0L7_2aSomeS837;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L15next__in__chainS839 =
        _M0L4_2aeS838->$1;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2364 = 0;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3767 =
        _M0L4_2aeS838->$1;
      if (_M0L15next__in__chainS839) {
        moonbit_incref(_M0L15next__in__chainS839);
      }
      if (_M0L6_2aoldS3767) {
        moonbit_decref(_M0L6_2aoldS3767);
      }
      _M0L4_2aeS838->$1 = _M0L6_2atmpS2364;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS833, _M0L4_2aeS838);
      moonbit_decref(_M0L4_2aeS838);
      _M0L1xS835 = _M0L15next__in__chainS839;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS841
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L9old__headS840;
  int32_t _M0L8capacityS2373;
  int32_t _M0L13new__capacityS842;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2367;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L6_2atmpS2366;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L6_2aoldS3775;
  int32_t _M0L6_2atmpS2368;
  int32_t _M0L8capacityS2370;
  int32_t _M0L6_2atmpS2369;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2371;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3774;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L1xS843;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS840 = _M0L4selfS841->$5;
  _M0L8capacityS2373 = _M0L4selfS841->$2;
  _M0L13new__capacityS842 = _M0L8capacityS2373 << 1;
  _M0L6_2atmpS2367 = 0;
  _M0L6_2atmpS2366
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE**)moonbit_make_ref_array(_M0L13new__capacityS842, _M0L6_2atmpS2367);
  _M0L6_2aoldS3775 = _M0L4selfS841->$0;
  if (_M0L9old__headS840) {
    moonbit_incref(_M0L9old__headS840);
  }
  moonbit_decref(_M0L6_2aoldS3775);
  _M0L4selfS841->$0 = _M0L6_2atmpS2366;
  _M0L4selfS841->$2 = _M0L13new__capacityS842;
  _M0L6_2atmpS2368 = _M0L13new__capacityS842 - 1;
  _M0L4selfS841->$3 = _M0L6_2atmpS2368;
  _M0L8capacityS2370 = _M0L4selfS841->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2369 = _M0FPB21calc__grow__threshold(_M0L8capacityS2370);
  _M0L4selfS841->$4 = _M0L6_2atmpS2369;
  _M0L4selfS841->$1 = 0;
  _M0L6_2atmpS2371 = 0;
  _M0L6_2aoldS3774 = _M0L4selfS841->$5;
  if (_M0L6_2aoldS3774) {
    moonbit_decref(_M0L6_2aoldS3774);
  }
  _M0L4selfS841->$5 = _M0L6_2atmpS2371;
  _M0L4selfS841->$6 = -1;
  _M0L1xS843 = _M0L9old__headS840;
  while (1) {
    if (_M0L1xS843 == 0) {
      if (_M0L1xS843) {
        moonbit_decref(_M0L1xS843);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS845 = _M0L1xS843;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4_2aeS846 = _M0L7_2aSomeS845;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L15next__in__chainS847 =
        _M0L4_2aeS846->$1;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2372 = 0;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3772 =
        _M0L4_2aeS846->$1;
      if (_M0L15next__in__chainS847) {
        moonbit_incref(_M0L15next__in__chainS847);
      }
      if (_M0L6_2aoldS3772) {
        moonbit_decref(_M0L6_2aoldS3772);
      }
      _M0L4_2aeS846->$1 = _M0L6_2atmpS2372;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(_M0L4selfS841, _M0L4_2aeS846);
      moonbit_decref(_M0L4_2aeS846);
      _M0L1xS843 = _M0L15next__in__chainS847;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS849
) {
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L9old__headS848;
  int32_t _M0L8capacityS2381;
  int32_t _M0L13new__capacityS850;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2375;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L6_2atmpS2374;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L6_2aoldS3780;
  int32_t _M0L6_2atmpS2376;
  int32_t _M0L8capacityS2378;
  int32_t _M0L6_2atmpS2377;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2379;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS3779;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L1xS851;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS848 = _M0L4selfS849->$5;
  _M0L8capacityS2381 = _M0L4selfS849->$2;
  _M0L13new__capacityS850 = _M0L8capacityS2381 << 1;
  _M0L6_2atmpS2375 = 0;
  _M0L6_2atmpS2374
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE**)moonbit_make_ref_array(_M0L13new__capacityS850, _M0L6_2atmpS2375);
  _M0L6_2aoldS3780 = _M0L4selfS849->$0;
  if (_M0L9old__headS848) {
    moonbit_incref(_M0L9old__headS848);
  }
  moonbit_decref(_M0L6_2aoldS3780);
  _M0L4selfS849->$0 = _M0L6_2atmpS2374;
  _M0L4selfS849->$2 = _M0L13new__capacityS850;
  _M0L6_2atmpS2376 = _M0L13new__capacityS850 - 1;
  _M0L4selfS849->$3 = _M0L6_2atmpS2376;
  _M0L8capacityS2378 = _M0L4selfS849->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2377 = _M0FPB21calc__grow__threshold(_M0L8capacityS2378);
  _M0L4selfS849->$4 = _M0L6_2atmpS2377;
  _M0L4selfS849->$1 = 0;
  _M0L6_2atmpS2379 = 0;
  _M0L6_2aoldS3779 = _M0L4selfS849->$5;
  if (_M0L6_2aoldS3779) {
    moonbit_decref(_M0L6_2aoldS3779);
  }
  _M0L4selfS849->$5 = _M0L6_2atmpS2379;
  _M0L4selfS849->$6 = -1;
  _M0L1xS851 = _M0L9old__headS848;
  while (1) {
    if (_M0L1xS851 == 0) {
      if (_M0L1xS851) {
        moonbit_decref(_M0L1xS851);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS853 =
        _M0L1xS851;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L4_2aeS854 =
        _M0L7_2aSomeS853;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L15next__in__chainS855 =
        _M0L4_2aeS854->$1;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2380 =
        0;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS3777 =
        _M0L4_2aeS854->$1;
      if (_M0L15next__in__chainS855) {
        moonbit_incref(_M0L15next__in__chainS855);
      }
      if (_M0L6_2aoldS3777) {
        moonbit_decref(_M0L6_2aoldS3777);
      }
      _M0L4_2aeS854->$1 = _M0L6_2atmpS2380;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS849, _M0L4_2aeS854);
      moonbit_decref(_M0L4_2aeS854);
      _M0L1xS851 = _M0L15next__in__chainS855;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS793,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5outerS789
) {
  int32_t _M0L4hashS788;
  int32_t _M0L14capacity__maskS2319;
  int32_t _M0L6_2atmpS2318;
  int32_t _M0L3pslS790;
  int32_t _M0L3idxS791;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS788 = _M0L5outerS789->$3;
  _M0L14capacity__maskS2319 = _M0L4selfS793->$3;
  _M0L6_2atmpS2318 = _M0L4hashS788 & _M0L14capacity__maskS2319;
  _M0L3pslS790 = 0;
  _M0L3idxS791 = _M0L6_2atmpS2318;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2317 =
      _M0L4selfS793->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS792;
    if (
      _M0L3idxS791 < 0
      || _M0L3idxS791 >= Moonbit_array_length(_M0L7entriesS2317)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS792
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2317[
        _M0L3idxS791
      ];
    if (_M0L7_2abindS792 == 0) {
      int32_t _M0L4tailS2310;
      _M0L5outerS789->$2 = _M0L3pslS790;
      _M0L4tailS2310 = _M0L4selfS793->$6;
      _M0L5outerS789->$0 = _M0L4tailS2310;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS793, _M0L3idxS791, _M0L5outerS789);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS794 =
        _M0L7_2abindS792;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2acurrS795 =
        _M0L7_2aSomeS794;
      int32_t _M0L3pslS2311 = _M0L7_2acurrS795->$2;
      if (_M0L3pslS790 > _M0L3pslS2311) {
        int32_t _M0L4tailS2312;
        moonbit_incref(_M0L7_2acurrS795);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS793, _M0L3idxS791, _M0L7_2acurrS795);
        moonbit_decref(_M0L7_2acurrS795);
        _M0L5outerS789->$2 = _M0L3pslS790;
        _M0L4tailS2312 = _M0L4selfS793->$6;
        _M0L5outerS789->$0 = _M0L4tailS2312;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS793, _M0L3idxS791, _M0L5outerS789);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2313 = _M0L3pslS790 + 1;
        int32_t _M0L6_2atmpS2315 = _M0L3idxS791 + 1;
        int32_t _M0L14capacity__maskS2316 = _M0L4selfS793->$3;
        int32_t _M0L6_2atmpS2314 =
          _M0L6_2atmpS2315 & _M0L14capacity__maskS2316;
        _M0L3pslS790 = _M0L6_2atmpS2313;
        _M0L3idxS791 = _M0L6_2atmpS2314;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS802,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5outerS798
) {
  int32_t _M0L4hashS797;
  int32_t _M0L14capacity__maskS2329;
  int32_t _M0L6_2atmpS2328;
  int32_t _M0L3pslS799;
  int32_t _M0L3idxS800;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS797 = _M0L5outerS798->$3;
  _M0L14capacity__maskS2329 = _M0L4selfS802->$3;
  _M0L6_2atmpS2328 = _M0L4hashS797 & _M0L14capacity__maskS2329;
  _M0L3pslS799 = 0;
  _M0L3idxS800 = _M0L6_2atmpS2328;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2327 =
      _M0L4selfS802->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS801;
    if (
      _M0L3idxS800 < 0
      || _M0L3idxS800 >= Moonbit_array_length(_M0L7entriesS2327)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS801
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2327[
        _M0L3idxS800
      ];
    if (_M0L7_2abindS801 == 0) {
      int32_t _M0L4tailS2320;
      _M0L5outerS798->$2 = _M0L3pslS799;
      _M0L4tailS2320 = _M0L4selfS802->$6;
      _M0L5outerS798->$0 = _M0L4tailS2320;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS802, _M0L3idxS800, _M0L5outerS798);
      return 0;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS803 =
        _M0L7_2abindS801;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2acurrS804 =
        _M0L7_2aSomeS803;
      int32_t _M0L3pslS2321 = _M0L7_2acurrS804->$2;
      if (_M0L3pslS799 > _M0L3pslS2321) {
        int32_t _M0L4tailS2322;
        moonbit_incref(_M0L7_2acurrS804);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS802, _M0L3idxS800, _M0L7_2acurrS804);
        moonbit_decref(_M0L7_2acurrS804);
        _M0L5outerS798->$2 = _M0L3pslS799;
        _M0L4tailS2322 = _M0L4selfS802->$6;
        _M0L5outerS798->$0 = _M0L4tailS2322;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS802, _M0L3idxS800, _M0L5outerS798);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2323 = _M0L3pslS799 + 1;
        int32_t _M0L6_2atmpS2325 = _M0L3idxS800 + 1;
        int32_t _M0L14capacity__maskS2326 = _M0L4selfS802->$3;
        int32_t _M0L6_2atmpS2324 =
          _M0L6_2atmpS2325 & _M0L14capacity__maskS2326;
        _M0L3pslS799 = _M0L6_2atmpS2323;
        _M0L3idxS800 = _M0L6_2atmpS2324;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS811,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5outerS807
) {
  int32_t _M0L4hashS806;
  int32_t _M0L14capacity__maskS2339;
  int32_t _M0L6_2atmpS2338;
  int32_t _M0L3pslS808;
  int32_t _M0L3idxS809;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS806 = _M0L5outerS807->$3;
  _M0L14capacity__maskS2339 = _M0L4selfS811->$3;
  _M0L6_2atmpS2338 = _M0L4hashS806 & _M0L14capacity__maskS2339;
  _M0L3pslS808 = 0;
  _M0L3idxS809 = _M0L6_2atmpS2338;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2337 =
      _M0L4selfS811->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS810;
    if (
      _M0L3idxS809 < 0
      || _M0L3idxS809 >= Moonbit_array_length(_M0L7entriesS2337)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS810
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2337[_M0L3idxS809];
    if (_M0L7_2abindS810 == 0) {
      int32_t _M0L4tailS2330;
      _M0L5outerS807->$2 = _M0L3pslS808;
      _M0L4tailS2330 = _M0L4selfS811->$6;
      _M0L5outerS807->$0 = _M0L4tailS2330;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS811, _M0L3idxS809, _M0L5outerS807);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS812 = _M0L7_2abindS810;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2acurrS813 = _M0L7_2aSomeS812;
      int32_t _M0L3pslS2331 = _M0L7_2acurrS813->$2;
      if (_M0L3pslS808 > _M0L3pslS2331) {
        int32_t _M0L4tailS2332;
        moonbit_incref(_M0L7_2acurrS813);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGiEE(_M0L4selfS811, _M0L3idxS809, _M0L7_2acurrS813);
        moonbit_decref(_M0L7_2acurrS813);
        _M0L5outerS807->$2 = _M0L3pslS808;
        _M0L4tailS2332 = _M0L4selfS811->$6;
        _M0L5outerS807->$0 = _M0L4tailS2332;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS811, _M0L3idxS809, _M0L5outerS807);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2333 = _M0L3pslS808 + 1;
        int32_t _M0L6_2atmpS2335 = _M0L3idxS809 + 1;
        int32_t _M0L14capacity__maskS2336 = _M0L4selfS811->$3;
        int32_t _M0L6_2atmpS2334 =
          _M0L6_2atmpS2335 & _M0L14capacity__maskS2336;
        _M0L3pslS808 = _M0L6_2atmpS2333;
        _M0L3idxS809 = _M0L6_2atmpS2334;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS820,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5outerS816
) {
  int32_t _M0L4hashS815;
  int32_t _M0L14capacity__maskS2349;
  int32_t _M0L6_2atmpS2348;
  int32_t _M0L3pslS817;
  int32_t _M0L3idxS818;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS815 = _M0L5outerS816->$3;
  _M0L14capacity__maskS2349 = _M0L4selfS820->$3;
  _M0L6_2atmpS2348 = _M0L4hashS815 & _M0L14capacity__maskS2349;
  _M0L3pslS817 = 0;
  _M0L3idxS818 = _M0L6_2atmpS2348;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2347 =
      _M0L4selfS820->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS819;
    if (
      _M0L3idxS818 < 0
      || _M0L3idxS818 >= Moonbit_array_length(_M0L7entriesS2347)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS819
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2347[
        _M0L3idxS818
      ];
    if (_M0L7_2abindS819 == 0) {
      int32_t _M0L4tailS2340;
      _M0L5outerS816->$2 = _M0L3pslS817;
      _M0L4tailS2340 = _M0L4selfS820->$6;
      _M0L5outerS816->$0 = _M0L4tailS2340;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS820, _M0L3idxS818, _M0L5outerS816);
      return 0;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS821 =
        _M0L7_2abindS819;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2acurrS822 =
        _M0L7_2aSomeS821;
      int32_t _M0L3pslS2341 = _M0L7_2acurrS822->$2;
      if (_M0L3pslS817 > _M0L3pslS2341) {
        int32_t _M0L4tailS2342;
        moonbit_incref(_M0L7_2acurrS822);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS820, _M0L3idxS818, _M0L7_2acurrS822);
        moonbit_decref(_M0L7_2acurrS822);
        _M0L5outerS816->$2 = _M0L3pslS817;
        _M0L4tailS2342 = _M0L4selfS820->$6;
        _M0L5outerS816->$0 = _M0L4tailS2342;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS820, _M0L3idxS818, _M0L5outerS816);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2343 = _M0L3pslS817 + 1;
        int32_t _M0L6_2atmpS2345 = _M0L3idxS818 + 1;
        int32_t _M0L14capacity__maskS2346 = _M0L4selfS820->$3;
        int32_t _M0L6_2atmpS2344 =
          _M0L6_2atmpS2345 & _M0L14capacity__maskS2346;
        _M0L3pslS817 = _M0L6_2atmpS2343;
        _M0L3idxS818 = _M0L6_2atmpS2344;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS752,
  int32_t _M0L3idxS757,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS756
) {
  int32_t _M0L3pslS2261;
  int32_t _M0L6_2atmpS2257;
  int32_t _M0L6_2atmpS2259;
  int32_t _M0L14capacity__maskS2260;
  int32_t _M0L6_2atmpS2258;
  int32_t _M0L3pslS748;
  int32_t _M0L3idxS749;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS750;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2261 = _M0L5entryS756->$2;
  _M0L6_2atmpS2257 = _M0L3pslS2261 + 1;
  _M0L6_2atmpS2259 = _M0L3idxS757 + 1;
  _M0L14capacity__maskS2260 = _M0L4selfS752->$3;
  _M0L6_2atmpS2258 = _M0L6_2atmpS2259 & _M0L14capacity__maskS2260;
  moonbit_incref(_M0L5entryS756);
  _M0L3pslS748 = _M0L6_2atmpS2257;
  _M0L3idxS749 = _M0L6_2atmpS2258;
  _M0L5entryS750 = _M0L5entryS756;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2256 =
      _M0L4selfS752->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS751;
    if (
      _M0L3idxS749 < 0
      || _M0L3idxS749 >= Moonbit_array_length(_M0L7entriesS2256)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS751
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2256[
        _M0L3idxS749
      ];
    if (_M0L7_2abindS751 == 0) {
      _M0L5entryS750->$2 = _M0L3pslS748;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS752, _M0L5entryS750, _M0L3idxS749);
      moonbit_decref(_M0L5entryS750);
      break;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS754 =
        _M0L7_2abindS751;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L14_2acurr__entryS755 =
        _M0L7_2aSomeS754;
      int32_t _M0L3pslS2246 = _M0L14_2acurr__entryS755->$2;
      if (_M0L3pslS748 > _M0L3pslS2246) {
        int32_t _M0L3pslS2251;
        int32_t _M0L6_2atmpS2247;
        int32_t _M0L6_2atmpS2249;
        int32_t _M0L14capacity__maskS2250;
        int32_t _M0L6_2atmpS2248;
        _M0L5entryS750->$2 = _M0L3pslS748;
        moonbit_incref(_M0L14_2acurr__entryS755);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS752, _M0L5entryS750, _M0L3idxS749);
        moonbit_decref(_M0L5entryS750);
        _M0L3pslS2251 = _M0L14_2acurr__entryS755->$2;
        _M0L6_2atmpS2247 = _M0L3pslS2251 + 1;
        _M0L6_2atmpS2249 = _M0L3idxS749 + 1;
        _M0L14capacity__maskS2250 = _M0L4selfS752->$3;
        _M0L6_2atmpS2248 = _M0L6_2atmpS2249 & _M0L14capacity__maskS2250;
        _M0L3pslS748 = _M0L6_2atmpS2247;
        _M0L3idxS749 = _M0L6_2atmpS2248;
        _M0L5entryS750 = _M0L14_2acurr__entryS755;
        continue;
      } else {
        int32_t _M0L6_2atmpS2252 = _M0L3pslS748 + 1;
        int32_t _M0L6_2atmpS2254 = _M0L3idxS749 + 1;
        int32_t _M0L14capacity__maskS2255 = _M0L4selfS752->$3;
        int32_t _M0L6_2atmpS2253 =
          _M0L6_2atmpS2254 & _M0L14capacity__maskS2255;
        struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _tmp_4513 =
          _M0L5entryS750;
        _M0L3pslS748 = _M0L6_2atmpS2252;
        _M0L3idxS749 = _M0L6_2atmpS2253;
        _M0L5entryS750 = _tmp_4513;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS762,
  int32_t _M0L3idxS767,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS766
) {
  int32_t _M0L3pslS2277;
  int32_t _M0L6_2atmpS2273;
  int32_t _M0L6_2atmpS2275;
  int32_t _M0L14capacity__maskS2276;
  int32_t _M0L6_2atmpS2274;
  int32_t _M0L3pslS758;
  int32_t _M0L3idxS759;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS760;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2277 = _M0L5entryS766->$2;
  _M0L6_2atmpS2273 = _M0L3pslS2277 + 1;
  _M0L6_2atmpS2275 = _M0L3idxS767 + 1;
  _M0L14capacity__maskS2276 = _M0L4selfS762->$3;
  _M0L6_2atmpS2274 = _M0L6_2atmpS2275 & _M0L14capacity__maskS2276;
  moonbit_incref(_M0L5entryS766);
  _M0L3pslS758 = _M0L6_2atmpS2273;
  _M0L3idxS759 = _M0L6_2atmpS2274;
  _M0L5entryS760 = _M0L5entryS766;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2272 =
      _M0L4selfS762->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS761;
    if (
      _M0L3idxS759 < 0
      || _M0L3idxS759 >= Moonbit_array_length(_M0L7entriesS2272)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS761
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2272[
        _M0L3idxS759
      ];
    if (_M0L7_2abindS761 == 0) {
      _M0L5entryS760->$2 = _M0L3pslS758;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS762, _M0L5entryS760, _M0L3idxS759);
      moonbit_decref(_M0L5entryS760);
      break;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS764 =
        _M0L7_2abindS761;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L14_2acurr__entryS765 =
        _M0L7_2aSomeS764;
      int32_t _M0L3pslS2262 = _M0L14_2acurr__entryS765->$2;
      if (_M0L3pslS758 > _M0L3pslS2262) {
        int32_t _M0L3pslS2267;
        int32_t _M0L6_2atmpS2263;
        int32_t _M0L6_2atmpS2265;
        int32_t _M0L14capacity__maskS2266;
        int32_t _M0L6_2atmpS2264;
        _M0L5entryS760->$2 = _M0L3pslS758;
        moonbit_incref(_M0L14_2acurr__entryS765);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS762, _M0L5entryS760, _M0L3idxS759);
        moonbit_decref(_M0L5entryS760);
        _M0L3pslS2267 = _M0L14_2acurr__entryS765->$2;
        _M0L6_2atmpS2263 = _M0L3pslS2267 + 1;
        _M0L6_2atmpS2265 = _M0L3idxS759 + 1;
        _M0L14capacity__maskS2266 = _M0L4selfS762->$3;
        _M0L6_2atmpS2264 = _M0L6_2atmpS2265 & _M0L14capacity__maskS2266;
        _M0L3pslS758 = _M0L6_2atmpS2263;
        _M0L3idxS759 = _M0L6_2atmpS2264;
        _M0L5entryS760 = _M0L14_2acurr__entryS765;
        continue;
      } else {
        int32_t _M0L6_2atmpS2268 = _M0L3pslS758 + 1;
        int32_t _M0L6_2atmpS2270 = _M0L3idxS759 + 1;
        int32_t _M0L14capacity__maskS2271 = _M0L4selfS762->$3;
        int32_t _M0L6_2atmpS2269 =
          _M0L6_2atmpS2270 & _M0L14capacity__maskS2271;
        struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _tmp_4515 =
          _M0L5entryS760;
        _M0L3pslS758 = _M0L6_2atmpS2268;
        _M0L3idxS759 = _M0L6_2atmpS2269;
        _M0L5entryS760 = _tmp_4515;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS772,
  int32_t _M0L3idxS777,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS776
) {
  int32_t _M0L3pslS2293;
  int32_t _M0L6_2atmpS2289;
  int32_t _M0L6_2atmpS2291;
  int32_t _M0L14capacity__maskS2292;
  int32_t _M0L6_2atmpS2290;
  int32_t _M0L3pslS768;
  int32_t _M0L3idxS769;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS770;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2293 = _M0L5entryS776->$2;
  _M0L6_2atmpS2289 = _M0L3pslS2293 + 1;
  _M0L6_2atmpS2291 = _M0L3idxS777 + 1;
  _M0L14capacity__maskS2292 = _M0L4selfS772->$3;
  _M0L6_2atmpS2290 = _M0L6_2atmpS2291 & _M0L14capacity__maskS2292;
  moonbit_incref(_M0L5entryS776);
  _M0L3pslS768 = _M0L6_2atmpS2289;
  _M0L3idxS769 = _M0L6_2atmpS2290;
  _M0L5entryS770 = _M0L5entryS776;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2288 =
      _M0L4selfS772->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS771;
    if (
      _M0L3idxS769 < 0
      || _M0L3idxS769 >= Moonbit_array_length(_M0L7entriesS2288)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS771
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2288[_M0L3idxS769];
    if (_M0L7_2abindS771 == 0) {
      _M0L5entryS770->$2 = _M0L3pslS768;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB5ArrayGiEE(_M0L4selfS772, _M0L5entryS770, _M0L3idxS769);
      moonbit_decref(_M0L5entryS770);
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS774 = _M0L7_2abindS771;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L14_2acurr__entryS775 =
        _M0L7_2aSomeS774;
      int32_t _M0L3pslS2278 = _M0L14_2acurr__entryS775->$2;
      if (_M0L3pslS768 > _M0L3pslS2278) {
        int32_t _M0L3pslS2283;
        int32_t _M0L6_2atmpS2279;
        int32_t _M0L6_2atmpS2281;
        int32_t _M0L14capacity__maskS2282;
        int32_t _M0L6_2atmpS2280;
        _M0L5entryS770->$2 = _M0L3pslS768;
        moonbit_incref(_M0L14_2acurr__entryS775);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB5ArrayGiEE(_M0L4selfS772, _M0L5entryS770, _M0L3idxS769);
        moonbit_decref(_M0L5entryS770);
        _M0L3pslS2283 = _M0L14_2acurr__entryS775->$2;
        _M0L6_2atmpS2279 = _M0L3pslS2283 + 1;
        _M0L6_2atmpS2281 = _M0L3idxS769 + 1;
        _M0L14capacity__maskS2282 = _M0L4selfS772->$3;
        _M0L6_2atmpS2280 = _M0L6_2atmpS2281 & _M0L14capacity__maskS2282;
        _M0L3pslS768 = _M0L6_2atmpS2279;
        _M0L3idxS769 = _M0L6_2atmpS2280;
        _M0L5entryS770 = _M0L14_2acurr__entryS775;
        continue;
      } else {
        int32_t _M0L6_2atmpS2284 = _M0L3pslS768 + 1;
        int32_t _M0L6_2atmpS2286 = _M0L3idxS769 + 1;
        int32_t _M0L14capacity__maskS2287 = _M0L4selfS772->$3;
        int32_t _M0L6_2atmpS2285 =
          _M0L6_2atmpS2286 & _M0L14capacity__maskS2287;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _tmp_4517 = _M0L5entryS770;
        _M0L3pslS768 = _M0L6_2atmpS2284;
        _M0L3idxS769 = _M0L6_2atmpS2285;
        _M0L5entryS770 = _tmp_4517;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS782,
  int32_t _M0L3idxS787,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS786
) {
  int32_t _M0L3pslS2309;
  int32_t _M0L6_2atmpS2305;
  int32_t _M0L6_2atmpS2307;
  int32_t _M0L14capacity__maskS2308;
  int32_t _M0L6_2atmpS2306;
  int32_t _M0L3pslS778;
  int32_t _M0L3idxS779;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS780;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2309 = _M0L5entryS786->$2;
  _M0L6_2atmpS2305 = _M0L3pslS2309 + 1;
  _M0L6_2atmpS2307 = _M0L3idxS787 + 1;
  _M0L14capacity__maskS2308 = _M0L4selfS782->$3;
  _M0L6_2atmpS2306 = _M0L6_2atmpS2307 & _M0L14capacity__maskS2308;
  moonbit_incref(_M0L5entryS786);
  _M0L3pslS778 = _M0L6_2atmpS2305;
  _M0L3idxS779 = _M0L6_2atmpS2306;
  _M0L5entryS780 = _M0L5entryS786;
  while (1) {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2304 =
      _M0L4selfS782->$0;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS781;
    if (
      _M0L3idxS779 < 0
      || _M0L3idxS779 >= Moonbit_array_length(_M0L7entriesS2304)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS781
    = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2304[
        _M0L3idxS779
      ];
    if (_M0L7_2abindS781 == 0) {
      _M0L5entryS780->$2 = _M0L3pslS778;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS782, _M0L5entryS780, _M0L3idxS779);
      moonbit_decref(_M0L5entryS780);
      break;
    } else {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS784 =
        _M0L7_2abindS781;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L14_2acurr__entryS785 =
        _M0L7_2aSomeS784;
      int32_t _M0L3pslS2294 = _M0L14_2acurr__entryS785->$2;
      if (_M0L3pslS778 > _M0L3pslS2294) {
        int32_t _M0L3pslS2299;
        int32_t _M0L6_2atmpS2295;
        int32_t _M0L6_2atmpS2297;
        int32_t _M0L14capacity__maskS2298;
        int32_t _M0L6_2atmpS2296;
        _M0L5entryS780->$2 = _M0L3pslS778;
        moonbit_incref(_M0L14_2acurr__entryS785);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(_M0L4selfS782, _M0L5entryS780, _M0L3idxS779);
        moonbit_decref(_M0L5entryS780);
        _M0L3pslS2299 = _M0L14_2acurr__entryS785->$2;
        _M0L6_2atmpS2295 = _M0L3pslS2299 + 1;
        _M0L6_2atmpS2297 = _M0L3idxS779 + 1;
        _M0L14capacity__maskS2298 = _M0L4selfS782->$3;
        _M0L6_2atmpS2296 = _M0L6_2atmpS2297 & _M0L14capacity__maskS2298;
        _M0L3pslS778 = _M0L6_2atmpS2295;
        _M0L3idxS779 = _M0L6_2atmpS2296;
        _M0L5entryS780 = _M0L14_2acurr__entryS785;
        continue;
      } else {
        int32_t _M0L6_2atmpS2300 = _M0L3pslS778 + 1;
        int32_t _M0L6_2atmpS2302 = _M0L3idxS779 + 1;
        int32_t _M0L14capacity__maskS2303 = _M0L4selfS782->$3;
        int32_t _M0L6_2atmpS2301 =
          _M0L6_2atmpS2302 & _M0L14capacity__maskS2303;
        struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _tmp_4519 =
          _M0L5entryS780;
        _M0L3pslS778 = _M0L6_2atmpS2300;
        _M0L3idxS779 = _M0L6_2atmpS2301;
        _M0L5entryS780 = _tmp_4519;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS724,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS726,
  int32_t _M0L8new__idxS725
) {
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2238;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2239;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3799;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS727;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2238 = _M0L4selfS724->$0;
  _M0L6_2atmpS2239 = _M0L5entryS726;
  if (
    _M0L8new__idxS725 < 0
    || _M0L8new__idxS725 >= Moonbit_array_length(_M0L7entriesS2238)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3799
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2238[
      _M0L8new__idxS725
    ];
  if (_M0L6_2atmpS2239) {
    moonbit_incref(_M0L6_2atmpS2239);
  }
  if (_M0L6_2aoldS3799) {
    moonbit_decref(_M0L6_2aoldS3799);
  }
  _M0L7entriesS2238[_M0L8new__idxS725] = _M0L6_2atmpS2239;
  _M0L7_2abindS727 = _M0L5entryS726->$1;
  if (_M0L7_2abindS727 == 0) {
    _M0L4selfS724->$6 = _M0L8new__idxS725;
  } else {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS728 =
      _M0L7_2abindS727;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2anextS729 =
      _M0L7_2aSomeS728;
    _M0L7_2anextS729->$0 = _M0L8new__idxS725;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS730,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS732,
  int32_t _M0L8new__idxS731
) {
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2240;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2241;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3802;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS733;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2240 = _M0L4selfS730->$0;
  _M0L6_2atmpS2241 = _M0L5entryS732;
  if (
    _M0L8new__idxS731 < 0
    || _M0L8new__idxS731 >= Moonbit_array_length(_M0L7entriesS2240)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3802
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2240[
      _M0L8new__idxS731
    ];
  if (_M0L6_2atmpS2241) {
    moonbit_incref(_M0L6_2atmpS2241);
  }
  if (_M0L6_2aoldS3802) {
    moonbit_decref(_M0L6_2aoldS3802);
  }
  _M0L7entriesS2240[_M0L8new__idxS731] = _M0L6_2atmpS2241;
  _M0L7_2abindS733 = _M0L5entryS732->$1;
  if (_M0L7_2abindS733 == 0) {
    _M0L4selfS730->$6 = _M0L8new__idxS731;
  } else {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS734 =
      _M0L7_2abindS733;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2anextS735 =
      _M0L7_2aSomeS734;
    _M0L7_2anextS735->$0 = _M0L8new__idxS731;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS736,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS738,
  int32_t _M0L8new__idxS737
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2242;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2243;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3805;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS739;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2242 = _M0L4selfS736->$0;
  _M0L6_2atmpS2243 = _M0L5entryS738;
  if (
    _M0L8new__idxS737 < 0
    || _M0L8new__idxS737 >= Moonbit_array_length(_M0L7entriesS2242)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3805
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2242[_M0L8new__idxS737];
  if (_M0L6_2atmpS2243) {
    moonbit_incref(_M0L6_2atmpS2243);
  }
  if (_M0L6_2aoldS3805) {
    moonbit_decref(_M0L6_2aoldS3805);
  }
  _M0L7entriesS2242[_M0L8new__idxS737] = _M0L6_2atmpS2243;
  _M0L7_2abindS739 = _M0L5entryS738->$1;
  if (_M0L7_2abindS739 == 0) {
    _M0L4selfS736->$6 = _M0L8new__idxS737;
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS740 = _M0L7_2abindS739;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2anextS741 = _M0L7_2aSomeS740;
    _M0L7_2anextS741->$0 = _M0L8new__idxS737;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS742,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS744,
  int32_t _M0L8new__idxS743
) {
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2244;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2245;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS3808;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS745;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2244 = _M0L4selfS742->$0;
  _M0L6_2atmpS2245 = _M0L5entryS744;
  if (
    _M0L8new__idxS743 < 0
    || _M0L8new__idxS743 >= Moonbit_array_length(_M0L7entriesS2244)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3808
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2244[
      _M0L8new__idxS743
    ];
  if (_M0L6_2atmpS2245) {
    moonbit_incref(_M0L6_2atmpS2245);
  }
  if (_M0L6_2aoldS3808) {
    moonbit_decref(_M0L6_2aoldS3808);
  }
  _M0L7entriesS2244[_M0L8new__idxS743] = _M0L6_2atmpS2245;
  _M0L7_2abindS745 = _M0L5entryS744->$1;
  if (_M0L7_2abindS745 == 0) {
    _M0L4selfS742->$6 = _M0L8new__idxS743;
  } else {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS746 =
      _M0L7_2abindS745;
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2anextS747 =
      _M0L7_2aSomeS746;
    _M0L7_2anextS747->$0 = _M0L8new__idxS743;
  }
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS709,
  int32_t _M0L3idxS711,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS710
) {
  int32_t _M0L7_2abindS708;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2207;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2208;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3810;
  int32_t _M0L4sizeS2210;
  int32_t _M0L6_2atmpS2209;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS708 = _M0L4selfS709->$6;
  switch (_M0L7_2abindS708) {
    case -1: {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2202 =
        _M0L5entryS710;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3812 =
        _M0L4selfS709->$5;
      if (_M0L6_2atmpS2202) {
        moonbit_incref(_M0L6_2atmpS2202);
      }
      if (_M0L6_2aoldS3812) {
        moonbit_decref(_M0L6_2aoldS3812);
      }
      _M0L4selfS709->$5 = _M0L6_2atmpS2202;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2206 =
        _M0L4selfS709->$0;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2205;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2203;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2204;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3813;
      if (
        _M0L7_2abindS708 < 0
        || _M0L7_2abindS708 >= Moonbit_array_length(_M0L7entriesS2206)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2205
      = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2206[
          _M0L7_2abindS708
        ];
      if (_M0L6_2atmpS2205) {
        moonbit_incref(_M0L6_2atmpS2205);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2203
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE(_M0L6_2atmpS2205);
      if (_M0L6_2atmpS2205) {
        moonbit_decref(_M0L6_2atmpS2205);
      }
      _M0L6_2atmpS2204 = _M0L5entryS710;
      _M0L6_2aoldS3813 = _M0L6_2atmpS2203->$1;
      if (_M0L6_2atmpS2204) {
        moonbit_incref(_M0L6_2atmpS2204);
      }
      if (_M0L6_2aoldS3813) {
        moonbit_decref(_M0L6_2aoldS3813);
      }
      _M0L6_2atmpS2203->$1 = _M0L6_2atmpS2204;
      moonbit_decref(_M0L6_2atmpS2203);
      break;
    }
  }
  _M0L4selfS709->$6 = _M0L3idxS711;
  _M0L7entriesS2207 = _M0L4selfS709->$0;
  _M0L6_2atmpS2208 = _M0L5entryS710;
  if (
    _M0L3idxS711 < 0
    || _M0L3idxS711 >= Moonbit_array_length(_M0L7entriesS2207)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3810
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2207[
      _M0L3idxS711
    ];
  if (_M0L6_2atmpS2208) {
    moonbit_incref(_M0L6_2atmpS2208);
  }
  if (_M0L6_2aoldS3810) {
    moonbit_decref(_M0L6_2aoldS3810);
  }
  _M0L7entriesS2207[_M0L3idxS711] = _M0L6_2atmpS2208;
  _M0L4sizeS2210 = _M0L4selfS709->$1;
  _M0L6_2atmpS2209 = _M0L4sizeS2210 + 1;
  _M0L4selfS709->$1 = _M0L6_2atmpS2209;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS713,
  int32_t _M0L3idxS715,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS714
) {
  int32_t _M0L7_2abindS712;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2216;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2217;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3816;
  int32_t _M0L4sizeS2219;
  int32_t _M0L6_2atmpS2218;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS712 = _M0L4selfS713->$6;
  switch (_M0L7_2abindS712) {
    case -1: {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2211 =
        _M0L5entryS714;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3818 =
        _M0L4selfS713->$5;
      if (_M0L6_2atmpS2211) {
        moonbit_incref(_M0L6_2atmpS2211);
      }
      if (_M0L6_2aoldS3818) {
        moonbit_decref(_M0L6_2aoldS3818);
      }
      _M0L4selfS713->$5 = _M0L6_2atmpS2211;
      break;
    }
    default: {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2215 =
        _M0L4selfS713->$0;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2214;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2212;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2213;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3819;
      if (
        _M0L7_2abindS712 < 0
        || _M0L7_2abindS712 >= Moonbit_array_length(_M0L7entriesS2215)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2214
      = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2215[
          _M0L7_2abindS712
        ];
      if (_M0L6_2atmpS2214) {
        moonbit_incref(_M0L6_2atmpS2214);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2212
      = _M0MPC16option6Option6unwrapGRPB5EntryGiUWEuQRPC15error5ErrorNsEEE(_M0L6_2atmpS2214);
      if (_M0L6_2atmpS2214) {
        moonbit_decref(_M0L6_2atmpS2214);
      }
      _M0L6_2atmpS2213 = _M0L5entryS714;
      _M0L6_2aoldS3819 = _M0L6_2atmpS2212->$1;
      if (_M0L6_2atmpS2213) {
        moonbit_incref(_M0L6_2atmpS2213);
      }
      if (_M0L6_2aoldS3819) {
        moonbit_decref(_M0L6_2aoldS3819);
      }
      _M0L6_2atmpS2212->$1 = _M0L6_2atmpS2213;
      moonbit_decref(_M0L6_2atmpS2212);
      break;
    }
  }
  _M0L4selfS713->$6 = _M0L3idxS715;
  _M0L7entriesS2216 = _M0L4selfS713->$0;
  _M0L6_2atmpS2217 = _M0L5entryS714;
  if (
    _M0L3idxS715 < 0
    || _M0L3idxS715 >= Moonbit_array_length(_M0L7entriesS2216)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3816
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2216[
      _M0L3idxS715
    ];
  if (_M0L6_2atmpS2217) {
    moonbit_incref(_M0L6_2atmpS2217);
  }
  if (_M0L6_2aoldS3816) {
    moonbit_decref(_M0L6_2aoldS3816);
  }
  _M0L7entriesS2216[_M0L3idxS715] = _M0L6_2atmpS2217;
  _M0L4sizeS2219 = _M0L4selfS713->$1;
  _M0L6_2atmpS2218 = _M0L4sizeS2219 + 1;
  _M0L4selfS713->$1 = _M0L6_2atmpS2218;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS717,
  int32_t _M0L3idxS719,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS718
) {
  int32_t _M0L7_2abindS716;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2225;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2226;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3822;
  int32_t _M0L4sizeS2228;
  int32_t _M0L6_2atmpS2227;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS716 = _M0L4selfS717->$6;
  switch (_M0L7_2abindS716) {
    case -1: {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2220 = _M0L5entryS718;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3824 =
        _M0L4selfS717->$5;
      if (_M0L6_2atmpS2220) {
        moonbit_incref(_M0L6_2atmpS2220);
      }
      if (_M0L6_2aoldS3824) {
        moonbit_decref(_M0L6_2aoldS3824);
      }
      _M0L4selfS717->$5 = _M0L6_2atmpS2220;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2224 =
        _M0L4selfS717->$0;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2223;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2221;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2222;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3825;
      if (
        _M0L7_2abindS716 < 0
        || _M0L7_2abindS716 >= Moonbit_array_length(_M0L7entriesS2224)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2223
      = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2224[
          _M0L7_2abindS716
        ];
      if (_M0L6_2atmpS2223) {
        moonbit_incref(_M0L6_2atmpS2223);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2221
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(_M0L6_2atmpS2223);
      if (_M0L6_2atmpS2223) {
        moonbit_decref(_M0L6_2atmpS2223);
      }
      _M0L6_2atmpS2222 = _M0L5entryS718;
      _M0L6_2aoldS3825 = _M0L6_2atmpS2221->$1;
      if (_M0L6_2atmpS2222) {
        moonbit_incref(_M0L6_2atmpS2222);
      }
      if (_M0L6_2aoldS3825) {
        moonbit_decref(_M0L6_2aoldS3825);
      }
      _M0L6_2atmpS2221->$1 = _M0L6_2atmpS2222;
      moonbit_decref(_M0L6_2atmpS2221);
      break;
    }
  }
  _M0L4selfS717->$6 = _M0L3idxS719;
  _M0L7entriesS2225 = _M0L4selfS717->$0;
  _M0L6_2atmpS2226 = _M0L5entryS718;
  if (
    _M0L3idxS719 < 0
    || _M0L3idxS719 >= Moonbit_array_length(_M0L7entriesS2225)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3822
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2225[_M0L3idxS719];
  if (_M0L6_2atmpS2226) {
    moonbit_incref(_M0L6_2atmpS2226);
  }
  if (_M0L6_2aoldS3822) {
    moonbit_decref(_M0L6_2aoldS3822);
  }
  _M0L7entriesS2225[_M0L3idxS719] = _M0L6_2atmpS2226;
  _M0L4sizeS2228 = _M0L4selfS717->$1;
  _M0L6_2atmpS2227 = _M0L4sizeS2228 + 1;
  _M0L4selfS717->$1 = _M0L6_2atmpS2227;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS721,
  int32_t _M0L3idxS723,
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L5entryS722
) {
  int32_t _M0L7_2abindS720;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2234;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2235;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS3828;
  int32_t _M0L4sizeS2237;
  int32_t _M0L6_2atmpS2236;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS720 = _M0L4selfS721->$6;
  switch (_M0L7_2abindS720) {
    case -1: {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2229 =
        _M0L5entryS722;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS3830 =
        _M0L4selfS721->$5;
      if (_M0L6_2atmpS2229) {
        moonbit_incref(_M0L6_2atmpS2229);
      }
      if (_M0L6_2aoldS3830) {
        moonbit_decref(_M0L6_2aoldS3830);
      }
      _M0L4selfS721->$5 = _M0L6_2atmpS2229;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7entriesS2233 =
        _M0L4selfS721->$0;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2232;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2230;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2231;
      struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2aoldS3831;
      if (
        _M0L7_2abindS720 < 0
        || _M0L7_2abindS720 >= Moonbit_array_length(_M0L7entriesS2233)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2232
      = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2233[
          _M0L7_2abindS720
        ];
      if (_M0L6_2atmpS2232) {
        moonbit_incref(_M0L6_2atmpS2232);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2230
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRP26biolab8bio__seq9SeqRecordEE(_M0L6_2atmpS2232);
      if (_M0L6_2atmpS2232) {
        moonbit_decref(_M0L6_2atmpS2232);
      }
      _M0L6_2atmpS2231 = _M0L5entryS722;
      _M0L6_2aoldS3831 = _M0L6_2atmpS2230->$1;
      if (_M0L6_2atmpS2231) {
        moonbit_incref(_M0L6_2atmpS2231);
      }
      if (_M0L6_2aoldS3831) {
        moonbit_decref(_M0L6_2aoldS3831);
      }
      _M0L6_2atmpS2230->$1 = _M0L6_2atmpS2231;
      moonbit_decref(_M0L6_2atmpS2230);
      break;
    }
  }
  _M0L4selfS721->$6 = _M0L3idxS723;
  _M0L7entriesS2234 = _M0L4selfS721->$0;
  _M0L6_2atmpS2235 = _M0L5entryS722;
  if (
    _M0L3idxS723 < 0
    || _M0L3idxS723 >= Moonbit_array_length(_M0L7entriesS2234)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3828
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE*)_M0L7entriesS2234[
      _M0L3idxS723
    ];
  if (_M0L6_2atmpS2235) {
    moonbit_incref(_M0L6_2atmpS2235);
  }
  if (_M0L6_2aoldS3828) {
    moonbit_decref(_M0L6_2aoldS3828);
  }
  _M0L7entriesS2234[_M0L3idxS723] = _M0L6_2atmpS2235;
  _M0L4sizeS2237 = _M0L4selfS721->$1;
  _M0L6_2atmpS2236 = _M0L4sizeS2237 + 1;
  _M0L4selfS721->$1 = _M0L6_2atmpS2236;
  return 0;
}

int32_t _M0MPC13int3Int3max(int32_t _M0L4selfS706, int32_t _M0L5otherS707) {
  #line 75 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS706 > _M0L5otherS707) {
    return _M0L4selfS706;
  } else {
    return _M0L5otherS707;
  }
}

int32_t _M0FPB21capacity__for__length(int32_t _M0L6lengthS705) {
  int32_t _M0Lm8capacityS704;
  int32_t _M0L6_2atmpS2200;
  int32_t _M0L6_2atmpS2199;
  #line 71 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 72 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0Lm8capacityS704 = _M0MPC13int3Int20next__power__of__two(_M0L6lengthS705);
  _M0L6_2atmpS2200 = _M0Lm8capacityS704;
  #line 73 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2199 = _M0FPB21calc__grow__threshold(_M0L6_2atmpS2200);
  if (_M0L6lengthS705 > _M0L6_2atmpS2199) {
    int32_t _M0L6_2atmpS2201 = _M0Lm8capacityS704;
    _M0Lm8capacityS704 = _M0L6_2atmpS2201 * 2;
  }
  return _M0Lm8capacityS704;
}

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0FPB8new__mapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  int32_t _M0L8capacityS681
) {
  int32_t _M0L8capacityS680;
  int32_t _M0L7_2abindS682;
  int32_t _M0L7_2abindS683;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2195;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7_2abindS684;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS685;
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _block_4520;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS680
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS681);
  _M0L7_2abindS682 = _M0L8capacityS680 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS683 = _M0FPB21calc__grow__threshold(_M0L8capacityS680);
  _M0L6_2atmpS2195 = 0;
  _M0L7_2abindS684
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE**)moonbit_make_ref_array(_M0L8capacityS680, _M0L6_2atmpS2195);
  _M0L7_2abindS685 = 0;
  _block_4520
  = (struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_block_4520)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 88, 0);
  _block_4520->$0 = _M0L7_2abindS684;
  _block_4520->$1 = 0;
  _block_4520->$2 = _M0L8capacityS680;
  _block_4520->$3 = _M0L7_2abindS682;
  _block_4520->$4 = _M0L7_2abindS683;
  _block_4520->$5 = _M0L7_2abindS685;
  _block_4520->$6 = -1;
  return _block_4520;
}

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0FPB8new__mapGiUWEuQRPC15error5ErrorNsEE(
  int32_t _M0L8capacityS687
) {
  int32_t _M0L8capacityS686;
  int32_t _M0L7_2abindS688;
  int32_t _M0L7_2abindS689;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2196;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS690;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS691;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _block_4521;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS686
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS687);
  _M0L7_2abindS688 = _M0L8capacityS686 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS689 = _M0FPB21calc__grow__threshold(_M0L8capacityS686);
  _M0L6_2atmpS2196 = 0;
  _M0L7_2abindS690
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE**)moonbit_make_ref_array(_M0L8capacityS686, _M0L6_2atmpS2196);
  _M0L7_2abindS691 = 0;
  _block_4521
  = (struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_block_4521)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 92, 0);
  _block_4521->$0 = _M0L7_2abindS690;
  _block_4521->$1 = 0;
  _block_4521->$2 = _M0L8capacityS686;
  _block_4521->$3 = _M0L7_2abindS688;
  _block_4521->$4 = _M0L7_2abindS689;
  _block_4521->$5 = _M0L7_2abindS691;
  _block_4521->$6 = -1;
  return _block_4521;
}

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0FPB8new__mapGsRPB5ArrayGiEE(
  int32_t _M0L8capacityS693
) {
  int32_t _M0L8capacityS692;
  int32_t _M0L7_2abindS694;
  int32_t _M0L7_2abindS695;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2197;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7_2abindS696;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS697;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _block_4522;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS692
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS693);
  _M0L7_2abindS694 = _M0L8capacityS692 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS695 = _M0FPB21calc__grow__threshold(_M0L8capacityS692);
  _M0L6_2atmpS2197 = 0;
  _M0L7_2abindS696
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE**)moonbit_make_ref_array(_M0L8capacityS692, _M0L6_2atmpS2197);
  _M0L7_2abindS697 = 0;
  _block_4522
  = (struct _M0TPB3MapGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB5ArrayGiEE));
  Moonbit_object_header(_block_4522)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 96, 0);
  _block_4522->$0 = _M0L7_2abindS696;
  _block_4522->$1 = 0;
  _block_4522->$2 = _M0L8capacityS692;
  _block_4522->$3 = _M0L7_2abindS694;
  _block_4522->$4 = _M0L7_2abindS695;
  _block_4522->$5 = _M0L7_2abindS697;
  _block_4522->$6 = -1;
  return _block_4522;
}

struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _M0FPB8new__mapGsRP26biolab8bio__seq9SeqRecordE(
  int32_t _M0L8capacityS699
) {
  int32_t _M0L8capacityS698;
  int32_t _M0L7_2abindS700;
  int32_t _M0L7_2abindS701;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2198;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE** _M0L7_2abindS702;
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS703;
  struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE* _block_4523;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS698
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS699);
  _M0L7_2abindS700 = _M0L8capacityS698 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS701 = _M0FPB21calc__grow__threshold(_M0L8capacityS698);
  _M0L6_2atmpS2198 = 0;
  _M0L7_2abindS702
  = (struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE**)moonbit_make_ref_array(_M0L8capacityS698, _M0L6_2atmpS2198);
  _M0L7_2abindS703 = 0;
  _block_4523
  = (struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_block_4523)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 100, 0);
  _block_4523->$0 = _M0L7_2abindS702;
  _block_4523->$1 = 0;
  _block_4523->$2 = _M0L8capacityS698;
  _block_4523->$3 = _M0L7_2abindS700;
  _block_4523->$4 = _M0L7_2abindS701;
  _block_4523->$5 = _M0L7_2abindS703;
  _block_4523->$6 = -1;
  return _block_4523;
}

int32_t _M0MPC13int3Int20next__power__of__two(int32_t _M0L4selfS679) {
  #line 33 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS679 >= 0) {
    int32_t _M0L6_2atmpS2194;
    int32_t _M0L6_2atmpS2193;
    int32_t _M0L6_2atmpS2192;
    int32_t _M0L6_2atmpS2191;
    if (_M0L4selfS679 <= 1) {
      return 1;
    }
    if (_M0L4selfS679 > 1073741824) {
      return 1073741824;
    }
    _M0L6_2atmpS2194 = _M0L4selfS679 - 1;
    #line 44 "/home/developer/.moon/lib/core/builtin/int.mbt"
    _M0L6_2atmpS2193 = moonbit_clz32(_M0L6_2atmpS2194);
    _M0L6_2atmpS2192 = _M0L6_2atmpS2193 - 1;
    _M0L6_2atmpS2191 = 2147483647 >> (_M0L6_2atmpS2192 & 31);
    return _M0L6_2atmpS2191 + 1;
  } else {
    #line 34 "/home/developer/.moon/lib/core/builtin/int.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB21calc__grow__threshold(int32_t _M0L8capacityS678) {
  int32_t _M0L6_2atmpS2190;
  #line 610 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2190 = _M0L8capacityS678 * 13;
  return _M0L6_2atmpS2190 / 16;
}

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE(
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS670
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS670 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS671 =
      _M0L4selfS670;
    if (_M0L7_2aSomeS671) {
      moonbit_incref(_M0L7_2aSomeS671);
    }
    return _M0L7_2aSomeS671;
  }
}

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0MPC16option6Option6unwrapGRPB5EntryGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS672
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS672 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS673 =
      _M0L4selfS672;
    if (_M0L7_2aSomeS673) {
      moonbit_incref(_M0L7_2aSomeS673);
    }
    return _M0L7_2aSomeS673;
  }
}

struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4selfS674
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS674 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS675 = _M0L4selfS674;
    if (_M0L7_2aSomeS675) {
      moonbit_incref(_M0L7_2aSomeS675);
    }
    return _M0L7_2aSomeS675;
  }
}

struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0MPC16option6Option6unwrapGRPB5EntryGsRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L4selfS676
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS676 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRP26biolab8bio__seq9SeqRecordE* _M0L7_2aSomeS677 =
      _M0L4selfS676;
    if (_M0L7_2aSomeS677) {
      moonbit_incref(_M0L7_2aSomeS677);
    }
    return _M0L7_2aSomeS677;
  }
}

struct _M0TPB4IterGsE* _M0MPC15array13ReadOnlyArray4iterGsE(
  moonbit_string_t* _M0L4selfS669
) {
  moonbit_string_t* _M0L6_2atmpS2189;
  struct _M0TPB4IterGsE* _result_4524;
  #line 194 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
  moonbit_incref(_M0L4selfS669);
  _M0L6_2atmpS2189 = _M0L4selfS669;
  #line 196 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
  _result_4524 = _M0MPC15array10FixedArray4iterGsE(_M0L6_2atmpS2189);
  moonbit_decref(_M0L6_2atmpS2189);
  return _result_4524;
}

struct _M0TPB4IterGsE* _M0MPC15array10FixedArray4iterGsE(
  moonbit_string_t* _M0L4selfS668
) {
  moonbit_string_t* _M0L6_2atmpS2187;
  int32_t _M0L6_2atmpS2188;
  struct _M0TPB9ArrayViewGsE _M0L6_2atmpS2186;
  struct _M0TPB4IterGsE* _result_4525;
  #line 1566 "/home/developer/.moon/lib/core/builtin/fixedarray.mbt"
  _M0L6_2atmpS2187 = _M0L4selfS668;
  _M0L6_2atmpS2188 = Moonbit_array_length(_M0L4selfS668);
  moonbit_incref(_M0L6_2atmpS2187);
  _M0L6_2atmpS2186
  = (struct _M0TPB9ArrayViewGsE){
    .$0 = _M0L6_2atmpS2187, .$1 = 0, .$2 = _M0L6_2atmpS2188
  };
  #line 1568 "/home/developer/.moon/lib/core/builtin/fixedarray.mbt"
  _result_4525 = _M0MPC15array9ArrayView4iterGsE(_M0L6_2atmpS2186);
  moonbit_decref(_M0L6_2atmpS2186.$0);
  return _result_4525;
}

struct _M0TPB4IterGsE* _M0MPC15array9ArrayView4iterGsE(
  struct _M0TPB9ArrayViewGsE _M0L4selfS666
) {
  struct _M0TPB8MutLocalGiE* _M0L1iS664;
  int32_t _M0L3endS2184;
  int32_t _M0L5startS2185;
  int32_t _M0L3lenS665;
  struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__* _closure_4526;
  struct _M0TWEOs* _M0L6_2atmpS2172;
  int64_t _M0L6_2atmpS2173;
  struct _M0TPB4IterGsE* _result_4527;
  #line 724 "/home/developer/.moon/lib/core/builtin/arrayview.mbt"
  _M0L1iS664
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS664)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS664->$0 = 0;
  _M0L3endS2184 = _M0L4selfS666.$2;
  _M0L5startS2185 = _M0L4selfS666.$1;
  _M0L3lenS665 = _M0L3endS2184 - _M0L5startS2185;
  moonbit_incref(_M0L4selfS666.$0);
  _closure_4526
  = (struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__*)moonbit_malloc(sizeof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__));
  Moonbit_object_header(_closure_4526)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 107, 0);
  _closure_4526->code = &_M0MPC15array9ArrayView4iterGsEC2174l729;
  _closure_4526->$0 = _M0L4selfS666;
  _closure_4526->$1 = _M0L3lenS665;
  _closure_4526->$2 = _M0L1iS664;
  _M0L6_2atmpS2172 = (struct _M0TWEOs*)_closure_4526;
  _M0L6_2atmpS2173 = (int64_t)_M0L3lenS665;
  #line 728 "/home/developer/.moon/lib/core/builtin/arrayview.mbt"
  _result_4527 = _M0MPB4Iter3newGsE(_M0L6_2atmpS2172, _M0L6_2atmpS2173);
  moonbit_decref(_M0L6_2atmpS2172);
  return _result_4527;
}

moonbit_string_t _M0MPC15array9ArrayView4iterGsEC2174l729(
  struct _M0TWEOs* _M0L6_2aenvS2175
) {
  struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__* _M0L14_2acasted__envS2176;
  struct _M0TPB8MutLocalGiE* _M0L1iS664;
  int32_t _M0L3lenS665;
  struct _M0TPB9ArrayViewGsE _M0L4selfS666;
  int32_t _M0L3valS2177;
  #line 729 "/home/developer/.moon/lib/core/builtin/arrayview.mbt"
  _M0L14_2acasted__envS2176
  = (struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2174__l729__*)_M0L6_2aenvS2175;
  _M0L1iS664 = _M0L14_2acasted__envS2176->$2;
  _M0L3lenS665 = _M0L14_2acasted__envS2176->$1;
  _M0L4selfS666 = _M0L14_2acasted__envS2176->$0;
  _M0L3valS2177 = _M0L1iS664->$0;
  if (_M0L3valS2177 < _M0L3lenS665) {
    moonbit_string_t* _M0L3bufS2180 = _M0L4selfS666.$0;
    int32_t _M0L5startS2182 = _M0L4selfS666.$1;
    int32_t _M0L3valS2183 = _M0L1iS664->$0;
    int32_t _M0L6_2atmpS2181 = _M0L5startS2182 + _M0L3valS2183;
    moonbit_string_t _M0L4elemS667 =
      (moonbit_string_t)_M0L3bufS2180[_M0L6_2atmpS2181];
    int32_t _M0L3valS2179 = _M0L1iS664->$0;
    int32_t _M0L6_2atmpS2178 = _M0L3valS2179 + 1;
    _M0L1iS664->$0 = _M0L6_2atmpS2178;
    moonbit_incref(_M0L4elemS667);
    return _M0L4elemS667;
  } else {
    return 0;
  }
}

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t _M0L4selfS663) {
  #line 35 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 36 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC13int3Int18to__string_2einner(_M0L4selfS663, 10);
}

moonbit_string_t _M0IPC14bool4BoolPB4Show10to__string(int32_t _M0L4selfS662) {
  #line 26 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L4selfS662) {
    return (moonbit_string_t)moonbit_string_literal_165.data;
  } else {
    return (moonbit_string_t)moonbit_string_literal_166.data;
  }
}

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPB4Iter9to__arrayGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS654
) {
  int64_t _M0L7_2abindS653;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L6resultS655;
  #line 841 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS653 = _M0L4selfS654->$1;
  if (_M0L7_2abindS653 == 4294967296ll) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS2171 =
      (struct _M0TPC16string10StringView*)moonbit_empty_ref_valtype_array;
    _M0L6resultS655
    = (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_M0L6resultS655)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 111, 0);
    _M0L6resultS655->$0 = _M0L6_2atmpS2171;
    _M0L6resultS655->$1 = 0;
  } else {
    int64_t _M0L7_2aSomeS656 = _M0L7_2abindS653;
    int32_t _M0L4_2anS657 = (int32_t)_M0L7_2aSomeS656;
    #line 844 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6resultS655
    = _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(_M0L4_2anS657);
  }
  _2awhile_661:;
  while (1) {
    void* _M0L7_2abindS658;
    #line 847 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L7_2abindS658
    = _M0MPB4Iter4nextGRPC16string10StringViewE(_M0L4selfS654);
    switch (Moonbit_object_tag(_M0L7_2abindS658)) {
      case 1: {
        struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS659 =
          (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS658;
        struct _M0TPC16string10StringView _M0L8_2afieldS3836 =
          _M0L7_2aSomeS659->$0;
        int32_t _M0L6_2acntS4176 =
          Moonbit_rc_count(Moonbit_object_header(_M0L7_2aSomeS659));
        struct _M0TPC16string10StringView _M0L4_2axS660;
        if (_M0L6_2acntS4176 > 1) {
          int32_t _M0L11_2anew__cntS4177 = _M0L6_2acntS4176 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2aSomeS659), _M0L11_2anew__cntS4177);
          moonbit_incref(_M0L8_2afieldS3836.$0);
        } else if (_M0L6_2acntS4176 == 1) {
          #line 847 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
          moonbit_free(_M0L7_2aSomeS659);
        }
        _M0L4_2axS660 = _M0L8_2afieldS3836;
        #line 848 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
        _M0MPC15array5Array4pushGRPC16string10StringViewE(_M0L6resultS655, _M0L4_2axS660);
        moonbit_decref(_M0L4_2axS660.$0);
        goto _2awhile_661;
        break;
      }
      default: {
        moonbit_decref(_M0L7_2abindS658);
        break;
      }
    }
    break;
  }
  return _M0L6resultS655;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string6String5split(
  moonbit_string_t _M0L4selfS651,
  struct _M0TPC16string10StringView _M0L3sepS652
) {
  int32_t _M0L6_2atmpS2170;
  struct _M0TPC16string10StringView _M0L6_2atmpS2169;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_4529;
  #line 1168 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2170 = Moonbit_array_length(_M0L4selfS651);
  moonbit_incref(_M0L4selfS651);
  _M0L6_2atmpS2169
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS651, .$1 = 0, .$2 = _M0L6_2atmpS2170
  };
  #line 1169 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4529
  = _M0MPC16string10StringView5split(_M0L6_2atmpS2169, _M0L3sepS652);
  moonbit_decref(_M0L6_2atmpS2169.$0);
  return _result_4529;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string10StringView5split(
  struct _M0TPC16string10StringView _M0L4selfS642,
  struct _M0TPC16string10StringView _M0L3sepS641
) {
  int32_t _M0L3endS2167;
  int32_t _M0L5startS2168;
  int32_t _M0L8sep__lenS640;
  void* _M0L4SomeS2166;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS644;
  struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__* _closure_4531;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS2156;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_4532;
  #line 1139 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2167 = _M0L3sepS641.$2;
  _M0L5startS2168 = _M0L3sepS641.$1;
  _M0L8sep__lenS640 = _M0L3endS2167 - _M0L5startS2168;
  if (_M0L8sep__lenS640 == 0) {
    struct _M0TPB4IterGcE* _M0L6_2atmpS2151;
    struct _M0TWcERPC16string10StringView* _M0L6_2atmpS2152;
    struct _M0TPB4IterGRPC16string10StringViewE* _result_4530;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS2151 = _M0MPC16string10StringView4iter(_M0L4selfS642);
    _M0L6_2atmpS2152
    = (struct _M0TWcERPC16string10StringView*)&_M0MPC16string10StringView5splitC2153l1145$closure.data;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_4530
    = _M0MPB4Iter3mapGcRPC16string10StringViewE(_M0L6_2atmpS2151, _M0L6_2atmpS2152);
    moonbit_decref(_M0L6_2atmpS2151);
    moonbit_decref(_M0L6_2atmpS2152);
    return _result_4530;
  }
  moonbit_incref(_M0L4selfS642.$0);
  _M0L4SomeS2166
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
  Moonbit_object_header(_M0L4SomeS2166)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 117, 1);
  ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS2166)->$0
  = _M0L4selfS642;
  _M0L9remainingS644
  = (struct _M0TPB8MutLocalGORPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE));
  Moonbit_object_header(_M0L9remainingS644)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 120, 0);
  _M0L9remainingS644->$0 = _M0L4SomeS2166;
  moonbit_incref(_M0L3sepS641.$0);
  _closure_4531
  = (struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__*)moonbit_malloc(sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__));
  Moonbit_object_header(_closure_4531)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 123, 0);
  _closure_4531->code = &_M0MPC16string10StringView5splitC2157l1148;
  _closure_4531->$0 = _M0L9remainingS644;
  _closure_4531->$1 = _M0L3sepS641;
  _closure_4531->$2 = _M0L8sep__lenS640;
  _M0L6_2atmpS2156
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_4531;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4532
  = _M0MPB4Iter3newGRPC16string10StringViewE(_M0L6_2atmpS2156, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS2156);
  return _result_4532;
}

void* _M0MPC16string10StringView5splitC2157l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS2158
) {
  struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__* _M0L14_2acasted__envS2159;
  int32_t _M0L8sep__lenS640;
  struct _M0TPC16string10StringView _M0L3sepS641;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS644;
  void* _M0L7_2abindS645;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L14_2acasted__envS2159
  = (struct _M0R44StringView_3a_3asplit_2eanon__u2157__l1148__*)_M0L6_2aenvS2158;
  _M0L8sep__lenS640 = _M0L14_2acasted__envS2159->$2;
  _M0L3sepS641 = _M0L14_2acasted__envS2159->$1;
  _M0L9remainingS644 = _M0L14_2acasted__envS2159->$0;
  _M0L7_2abindS645 = _M0L9remainingS644->$0;
  switch (Moonbit_object_tag(_M0L7_2abindS645)) {
    case 1: {
      struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS646 =
        (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS645;
      struct _M0TPC16string10StringView _M0L7_2aviewS647 =
        _M0L7_2aSomeS646->$0;
      int64_t _M0L7_2abindS648;
      moonbit_incref(_M0L7_2aviewS647.$0);
      #line 1150 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS648
      = _M0MPC16string10StringView4find(_M0L7_2aviewS647, _M0L3sepS641);
      if (_M0L7_2abindS648 == 4294967296ll) {
        void* _M0L4NoneS2160 =
          (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
        void* _M0L6_2aoldS3837 = _M0L9remainingS644->$0;
        void* _block_4533;
        moonbit_decref(_M0L6_2aoldS3837);
        _M0L9remainingS644->$0 = _M0L4NoneS2160;
        _block_4533
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_4533)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 117, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_4533)->$0
        = _M0L7_2aviewS647;
        return _block_4533;
      } else {
        int64_t _M0L7_2aSomeS649 = _M0L7_2abindS648;
        int32_t _M0L6_2aendS650 = (int32_t)_M0L7_2aSomeS649;
        int32_t _M0L6_2atmpS2163 = _M0L6_2aendS650 + _M0L8sep__lenS640;
        struct _M0TPC16string10StringView _M0L6_2atmpS2162;
        void* _M0L4SomeS2161;
        void* _M0L6_2aoldS3838;
        int64_t _M0L6_2atmpS2165;
        struct _M0TPC16string10StringView _M0L6_2atmpS2164;
        void* _block_4534;
        #line 1154 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2162
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS647, _M0L6_2atmpS2163, 4294967296ll);
        _M0L4SomeS2161
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_M0L4SomeS2161)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 117, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS2161)->$0
        = _M0L6_2atmpS2162;
        _M0L6_2aoldS3838 = _M0L9remainingS644->$0;
        moonbit_decref(_M0L6_2aoldS3838);
        _M0L9remainingS644->$0 = _M0L4SomeS2161;
        _M0L6_2atmpS2165 = (int64_t)_M0L6_2aendS650;
        #line 1155 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2164
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS647, 0, _M0L6_2atmpS2165);
        moonbit_decref(_M0L7_2aviewS647.$0);
        _block_4534
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_4534)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 117, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_4534)->$0
        = _M0L6_2atmpS2164;
        return _block_4534;
      }
      break;
    }
    default: {
      return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
      break;
    }
  }
}

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC2153l1145(
  struct _M0TWcERPC16string10StringView* _M0L6_2aenvS2154,
  int32_t _M0L1cS643
) {
  moonbit_string_t _M0L6_2atmpS2155;
  struct _M0TPC16string10StringView _result_4535;
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2155 = _M0IPC14char4CharPB4Show10to__string(_M0L1cS643);
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4535
  = _M0MPC16string6String12view_2einner(_M0L6_2atmpS2155, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS2155);
  return _result_4535;
}

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t _M0L4selfS639) {
  #line 446 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 447 "/home/developer/.moon/lib/core/builtin/char.mbt"
  return _M0FPB16char__to__string(_M0L4selfS639);
}

moonbit_string_t _M0FPB16char__to__string(int32_t _M0L4charS638) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS637;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS2150;
  moonbit_string_t _result_4536;
  #line 452 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0L7_2aselfS637 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS637, _M0L4charS638);
  _M0L6_2atmpS2150 = _M0L7_2aselfS637;
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _result_4536 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS2150);
  moonbit_decref(_M0L6_2atmpS2150);
  return _result_4536;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE* _M0L4selfS633,
  struct _M0TWcERPC16string10StringView* _M0L1fS636
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__* _closure_4537;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS2144;
  int64_t _M0L10size__hintS2145;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_4538;
  #line 389 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  moonbit_incref(_M0L1fS636);
  moonbit_incref(_M0L4selfS633);
  _closure_4537
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__*)moonbit_malloc(sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__));
  Moonbit_object_header(_closure_4537)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 127, 0);
  _closure_4537->code = &_M0MPB4Iter3mapGcRPC16string10StringViewEC2146l391;
  _closure_4537->$0 = _M0L1fS636;
  _closure_4537->$1 = _M0L4selfS633;
  _M0L6_2atmpS2144
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_4537;
  _M0L10size__hintS2145 = _M0L4selfS633->$1;
  _block_4538
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_4538)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 131, 0);
  _block_4538->$0 = _M0L6_2atmpS2144;
  _block_4538->$1 = _M0L10size__hintS2145;
  return _block_4538;
}

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC2146l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS2147
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__* _M0L14_2acasted__envS2148;
  struct _M0TPB4IterGcE* _M0L4selfS633;
  struct _M0TWcERPC16string10StringView* _M0L1fS636;
  int32_t _M0L7_2abindS632;
  #line 391 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L14_2acasted__envS2148
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2146__l391__*)_M0L6_2aenvS2147;
  _M0L4selfS633 = _M0L14_2acasted__envS2148->$1;
  _M0L1fS636 = _M0L14_2acasted__envS2148->$0;
  #line 392 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS632 = _M0MPB4Iter4nextGcE(_M0L4selfS633);
  if (_M0L7_2abindS632 == -1) {
    return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  } else {
    int32_t _M0L7_2aSomeS634 = _M0L7_2abindS632;
    int32_t _M0L4_2axS635 = _M0L7_2aSomeS634;
    struct _M0TPC16string10StringView _M0L6_2atmpS2149;
    void* _block_4539;
    #line 393 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6_2atmpS2149 = _M0L1fS636->code(_M0L1fS636, _M0L4_2axS635);
    _block_4539
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
    Moonbit_object_header(_block_4539)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 117, 1);
    ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_4539)->$0
    = _M0L6_2atmpS2149;
    return _block_4539;
  }
}

int32_t _M0MPC16string6String8contains(
  moonbit_string_t _M0L4selfS630,
  struct _M0TPC16string10StringView _M0L3strS631
) {
  int32_t _M0L6_2atmpS2143;
  struct _M0TPC16string10StringView _M0L6_2atmpS2142;
  int32_t _result_4540;
  #line 545 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2143 = Moonbit_array_length(_M0L4selfS630);
  moonbit_incref(_M0L4selfS630);
  _M0L6_2atmpS2142
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS630, .$1 = 0, .$2 = _M0L6_2atmpS2143
  };
  #line 546 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4540
  = _M0MPC16string10StringView8contains(_M0L6_2atmpS2142, _M0L3strS631);
  moonbit_decref(_M0L6_2atmpS2142.$0);
  return _result_4540;
}

int32_t _M0MPC16string10StringView8contains(
  struct _M0TPC16string10StringView _M0L4selfS628,
  struct _M0TPC16string10StringView _M0L3strS627
) {
  int32_t _M0L3endS2140;
  int32_t _M0L5startS2141;
  int32_t _M0L7_2abindS626;
  #line 535 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2140 = _M0L3strS627.$2;
  _M0L5startS2141 = _M0L3strS627.$1;
  _M0L7_2abindS626 = _M0L3endS2140 - _M0L5startS2141;
  switch (_M0L7_2abindS626) {
    case 0: {
      return 1;
      break;
    }
    
    case 1: {
      moonbit_string_t _M0L3strS2137 = _M0L3strS627.$0;
      int32_t _M0L5startS2138 = _M0L3strS627.$1;
      int32_t _M0L6_2atmpS2136 = _M0L3strS2137[_M0L5startS2138];
      #line 538 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0MPC16string10StringView20contains__code__unit(_M0L4selfS628, _M0L6_2atmpS2136);
      break;
    }
    default: {
      int64_t _M0L7_2abindS629;
      int32_t _M0L6_2atmpS2139;
      #line 539 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS629
      = _M0MPC16string10StringView4find(_M0L4selfS628, _M0L3strS627);
      _M0L6_2atmpS2139 = _M0L7_2abindS629 == 4294967296ll;
      return !_M0L6_2atmpS2139;
      break;
    }
  }
}

int32_t _M0MPC16string10StringView20contains__code__unit(
  struct _M0TPC16string10StringView _M0L4selfS622,
  int32_t _M0L4codeS624
) {
  int32_t _M0L3endS2134;
  int32_t _M0L5startS2135;
  int32_t _M0L7_2abindS621;
  int32_t _M0L1iS623;
  #line 524 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2134 = _M0L4selfS622.$2;
  _M0L5startS2135 = _M0L4selfS622.$1;
  _M0L7_2abindS621 = _M0L3endS2134 - _M0L5startS2135;
  _M0L1iS623 = 0;
  while (1) {
    if (_M0L1iS623 < _M0L7_2abindS621) {
      moonbit_string_t _M0L3strS2130 = _M0L4selfS622.$0;
      int32_t _M0L5startS2132 = _M0L4selfS622.$1;
      int32_t _M0L6_2atmpS2131 = _M0L5startS2132 + _M0L1iS623;
      int32_t _M0L6_2atmpS2129 = _M0L3strS2130[_M0L6_2atmpS2131];
      int32_t _M0L6_2atmpS2133;
      if (_M0L6_2atmpS2129 == _M0L4codeS624) {
        return 1;
      }
      _M0L6_2atmpS2133 = _M0L1iS623 + 1;
      _M0L1iS623 = _M0L6_2atmpS2133;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS606,
  moonbit_string_t _M0L5valueS608
) {
  int32_t _M0L3lenS2104;
  moonbit_string_t* _M0L6_2atmpS2106;
  int32_t _M0L6_2atmpS2105;
  int32_t _M0L6lengthS607;
  moonbit_string_t* _M0L3bufS2107;
  moonbit_string_t _M0L6_2aoldS3843;
  int32_t _M0L6_2atmpS2108;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2104 = _M0L4selfS606->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2106 = _M0MPC15array5Array6bufferGsE(_M0L4selfS606);
  _M0L6_2atmpS2105 = Moonbit_array_length(_M0L6_2atmpS2106);
  moonbit_decref(_M0L6_2atmpS2106);
  if (_M0L3lenS2104 == _M0L6_2atmpS2105) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGsE(_M0L4selfS606);
  }
  _M0L6lengthS607 = _M0L4selfS606->$1;
  _M0L3bufS2107 = _M0L4selfS606->$0;
  _M0L6_2aoldS3843 = (moonbit_string_t)_M0L3bufS2107[_M0L6lengthS607];
  moonbit_incref(_M0L5valueS608);
  moonbit_decref(_M0L6_2aoldS3843);
  _M0L3bufS2107[_M0L6lengthS607] = _M0L5valueS608;
  _M0L6_2atmpS2108 = _M0L6lengthS607 + 1;
  _M0L4selfS606->$1 = _M0L6_2atmpS2108;
  return 0;
}

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS609,
  struct _M0TUsiE* _M0L5valueS611
) {
  int32_t _M0L3lenS2109;
  struct _M0TUsiE** _M0L6_2atmpS2111;
  int32_t _M0L6_2atmpS2110;
  int32_t _M0L6lengthS610;
  struct _M0TUsiE** _M0L3bufS2112;
  struct _M0TUsiE* _M0L6_2aoldS3845;
  int32_t _M0L6_2atmpS2113;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2109 = _M0L4selfS609->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2111 = _M0MPC15array5Array6bufferGUsiEE(_M0L4selfS609);
  _M0L6_2atmpS2110 = Moonbit_array_length(_M0L6_2atmpS2111);
  moonbit_decref(_M0L6_2atmpS2111);
  if (_M0L3lenS2109 == _M0L6_2atmpS2110) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGUsiEE(_M0L4selfS609);
  }
  _M0L6lengthS610 = _M0L4selfS609->$1;
  _M0L3bufS2112 = _M0L4selfS609->$0;
  _M0L6_2aoldS3845 = (struct _M0TUsiE*)_M0L3bufS2112[_M0L6lengthS610];
  moonbit_incref(_M0L5valueS611);
  if (_M0L6_2aoldS3845) {
    moonbit_decref(_M0L6_2aoldS3845);
  }
  _M0L3bufS2112[_M0L6lengthS610] = _M0L5valueS611;
  _M0L6_2atmpS2113 = _M0L6lengthS610 + 1;
  _M0L4selfS609->$1 = _M0L6_2atmpS2113;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS612,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS614
) {
  int32_t _M0L3lenS2114;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2116;
  int32_t _M0L6_2atmpS2115;
  int32_t _M0L6lengthS613;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2117;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3847;
  int32_t _M0L6_2atmpS2118;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2114 = _M0L4selfS612->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2116
  = _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS612);
  _M0L6_2atmpS2115 = Moonbit_array_length(_M0L6_2atmpS2116);
  moonbit_decref(_M0L6_2atmpS2116);
  if (_M0L3lenS2114 == _M0L6_2atmpS2115) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS612);
  }
  _M0L6lengthS613 = _M0L4selfS612->$1;
  _M0L3bufS2117 = _M0L4selfS612->$0;
  _M0L6_2aoldS3847
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2117[_M0L6lengthS613];
  moonbit_incref(_M0L5valueS614);
  if (_M0L6_2aoldS3847) {
    moonbit_decref(_M0L6_2aoldS3847);
  }
  _M0L3bufS2117[_M0L6lengthS613] = _M0L5valueS614;
  _M0L6_2atmpS2118 = _M0L6lengthS613 + 1;
  _M0L4selfS612->$1 = _M0L6_2atmpS2118;
  return 0;
}

int32_t _M0MPC15array5Array4pushGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS615,
  int32_t _M0L5valueS617
) {
  int32_t _M0L3lenS2119;
  int32_t* _M0L6_2atmpS2121;
  int32_t _M0L6_2atmpS2120;
  int32_t _M0L6lengthS616;
  int32_t* _M0L3bufS2122;
  int32_t _M0L6_2atmpS2123;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2119 = _M0L4selfS615->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2121 = _M0MPC15array5Array6bufferGiE(_M0L4selfS615);
  _M0L6_2atmpS2120 = Moonbit_array_length(_M0L6_2atmpS2121);
  moonbit_decref(_M0L6_2atmpS2121);
  if (_M0L3lenS2119 == _M0L6_2atmpS2120) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGiE(_M0L4selfS615);
  }
  _M0L6lengthS616 = _M0L4selfS615->$1;
  _M0L3bufS2122 = _M0L4selfS615->$0;
  _M0L3bufS2122[_M0L6lengthS616] = _M0L5valueS617;
  _M0L6_2atmpS2123 = _M0L6lengthS616 + 1;
  _M0L4selfS615->$1 = _M0L6_2atmpS2123;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS618,
  struct _M0TPC16string10StringView _M0L5valueS620
) {
  int32_t _M0L3lenS2124;
  struct _M0TPC16string10StringView* _M0L6_2atmpS2126;
  int32_t _M0L6_2atmpS2125;
  int32_t _M0L6lengthS619;
  struct _M0TPC16string10StringView* _M0L3bufS2127;
  struct _M0TPC16string10StringView _M0L6_2aoldS3850;
  int32_t _M0L6_2atmpS2128;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2124 = _M0L4selfS618->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2126
  = _M0MPC15array5Array6bufferGRPC16string10StringViewE(_M0L4selfS618);
  _M0L6_2atmpS2125 = Moonbit_array_length(_M0L6_2atmpS2126);
  moonbit_decref(_M0L6_2atmpS2126);
  if (_M0L3lenS2124 == _M0L6_2atmpS2125) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRPC16string10StringViewE(_M0L4selfS618);
  }
  _M0L6lengthS619 = _M0L4selfS618->$1;
  _M0L3bufS2127 = _M0L4selfS618->$0;
  _M0L6_2aoldS3850 = _M0L3bufS2127[_M0L6lengthS619];
  moonbit_incref(_M0L5valueS620.$0);
  moonbit_decref(_M0L6_2aoldS3850.$0);
  _M0L3bufS2127[_M0L6lengthS619] = _M0L5valueS620;
  _M0L6_2atmpS2128 = _M0L6lengthS619 + 1;
  _M0L4selfS618->$1 = _M0L6_2atmpS2128;
  return 0;
}

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE* _M0L4selfS592) {
  int32_t _M0L8old__capS591;
  int32_t _M0L8new__capS593;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS591 = _M0L4selfS592->$1;
  if (_M0L8old__capS591 == 0) {
    _M0L8new__capS593 = 8;
  } else {
    _M0L8new__capS593 = _M0L8old__capS591 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGsE(_M0L4selfS592, _M0L8new__capS593);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS595
) {
  int32_t _M0L8old__capS594;
  int32_t _M0L8new__capS596;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS594 = _M0L4selfS595->$1;
  if (_M0L8old__capS594 == 0) {
    _M0L8new__capS596 = 8;
  } else {
    _M0L8new__capS596 = _M0L8old__capS594 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGUsiEE(_M0L4selfS595, _M0L8new__capS596);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS598
) {
  int32_t _M0L8old__capS597;
  int32_t _M0L8new__capS599;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS597 = _M0L4selfS598->$1;
  if (_M0L8old__capS597 == 0) {
    _M0L8new__capS599 = 8;
  } else {
    _M0L8new__capS599 = _M0L8old__capS597 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS598, _M0L8new__capS599);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGiE(struct _M0TPB5ArrayGiE* _M0L4selfS601) {
  int32_t _M0L8old__capS600;
  int32_t _M0L8new__capS602;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS600 = _M0L4selfS601->$1;
  if (_M0L8old__capS600 == 0) {
    _M0L8new__capS602 = 8;
  } else {
    _M0L8new__capS602 = _M0L8old__capS600 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGiE(_M0L4selfS601, _M0L8new__capS602);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS604
) {
  int32_t _M0L8old__capS603;
  int32_t _M0L8new__capS605;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS603 = _M0L4selfS604->$1;
  if (_M0L8old__capS603 == 0) {
    _M0L8new__capS605 = 8;
  } else {
    _M0L8new__capS605 = _M0L8old__capS603 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(_M0L4selfS604, _M0L8new__capS605);
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS562,
  int32_t _M0L13new__capacityS565
) {
  moonbit_string_t* _M0L8old__bufS561;
  int32_t _M0L8old__capS563;
  int32_t _M0L9copy__lenS564;
  moonbit_string_t* _M0L8new__bufS566;
  moonbit_string_t* _M0L6_2aoldS3852;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS561 = _M0L4selfS562->$0;
  _M0L8old__capS563 = Moonbit_array_length(_M0L8old__bufS561);
  if (_M0L8old__capS563 < _M0L13new__capacityS565) {
    _M0L9copy__lenS564 = _M0L8old__capS563;
  } else {
    _M0L9copy__lenS564 = _M0L13new__capacityS565;
  }
  moonbit_incref(_M0L8old__bufS561);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS566
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(_M0L8old__bufS561, _M0L13new__capacityS565, _M0L9copy__lenS564, 0, 0);
  moonbit_decref(_M0L8old__bufS561);
  _M0L6_2aoldS3852 = _M0L4selfS562->$0;
  moonbit_decref(_M0L6_2aoldS3852);
  _M0L4selfS562->$0 = _M0L8new__bufS566;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS568,
  int32_t _M0L13new__capacityS571
) {
  struct _M0TUsiE** _M0L8old__bufS567;
  int32_t _M0L8old__capS569;
  int32_t _M0L9copy__lenS570;
  struct _M0TUsiE** _M0L8new__bufS572;
  struct _M0TUsiE** _M0L6_2aoldS3854;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS567 = _M0L4selfS568->$0;
  _M0L8old__capS569 = Moonbit_array_length(_M0L8old__bufS567);
  if (_M0L8old__capS569 < _M0L13new__capacityS571) {
    _M0L9copy__lenS570 = _M0L8old__capS569;
  } else {
    _M0L9copy__lenS570 = _M0L13new__capacityS571;
  }
  moonbit_incref(_M0L8old__bufS567);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS572
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(_M0L8old__bufS567, _M0L13new__capacityS571, _M0L9copy__lenS570, 0, 0);
  moonbit_decref(_M0L8old__bufS567);
  _M0L6_2aoldS3854 = _M0L4selfS568->$0;
  moonbit_decref(_M0L6_2aoldS3854);
  _M0L4selfS568->$0 = _M0L8new__bufS572;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS574,
  int32_t _M0L13new__capacityS577
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8old__bufS573;
  int32_t _M0L8old__capS575;
  int32_t _M0L9copy__lenS576;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8new__bufS578;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2aoldS3856;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS573 = _M0L4selfS574->$0;
  _M0L8old__capS575 = Moonbit_array_length(_M0L8old__bufS573);
  if (_M0L8old__capS575 < _M0L13new__capacityS577) {
    _M0L9copy__lenS576 = _M0L8old__capS575;
  } else {
    _M0L9copy__lenS576 = _M0L13new__capacityS577;
  }
  moonbit_incref(_M0L8old__bufS573);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS578
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(_M0L8old__bufS573, _M0L13new__capacityS577, _M0L9copy__lenS576, 0, 0);
  moonbit_decref(_M0L8old__bufS573);
  _M0L6_2aoldS3856 = _M0L4selfS574->$0;
  moonbit_decref(_M0L6_2aoldS3856);
  _M0L4selfS574->$0 = _M0L8new__bufS578;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS580,
  int32_t _M0L13new__capacityS583
) {
  int32_t* _M0L8old__bufS579;
  int32_t _M0L8old__capS581;
  int32_t _M0L9copy__lenS582;
  int32_t* _M0L8new__bufS584;
  int32_t* _M0L6_2aoldS3858;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS579 = _M0L4selfS580->$0;
  _M0L8old__capS581 = Moonbit_array_length(_M0L8old__bufS579);
  if (_M0L8old__capS581 < _M0L13new__capacityS583) {
    _M0L9copy__lenS582 = _M0L8old__capS581;
  } else {
    _M0L9copy__lenS582 = _M0L13new__capacityS583;
  }
  moonbit_incref(_M0L8old__bufS579);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS584
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGiE(_M0L8old__bufS579, _M0L13new__capacityS583, _M0L9copy__lenS582, 0, 0);
  moonbit_decref(_M0L8old__bufS579);
  _M0L6_2aoldS3858 = _M0L4selfS580->$0;
  moonbit_decref(_M0L6_2aoldS3858);
  _M0L4selfS580->$0 = _M0L8new__bufS584;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS586,
  int32_t _M0L13new__capacityS589
) {
  struct _M0TPC16string10StringView* _M0L8old__bufS585;
  int32_t _M0L8old__capS587;
  int32_t _M0L9copy__lenS588;
  struct _M0TPC16string10StringView* _M0L8new__bufS590;
  struct _M0TPC16string10StringView* _M0L6_2aoldS3860;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS585 = _M0L4selfS586->$0;
  _M0L8old__capS587 = Moonbit_array_length(_M0L8old__bufS585);
  if (_M0L8old__capS587 < _M0L13new__capacityS589) {
    _M0L9copy__lenS588 = _M0L8old__capS587;
  } else {
    _M0L9copy__lenS588 = _M0L13new__capacityS589;
  }
  moonbit_incref(_M0L8old__bufS585);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS590
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(_M0L8old__bufS585, _M0L13new__capacityS589, _M0L9copy__lenS588, 0, 0);
  moonbit_decref(_M0L8old__bufS585);
  _M0L6_2aoldS3860 = _M0L4selfS586->$0;
  moonbit_decref(_M0L6_2aoldS3860);
  _M0L4selfS586->$0 = _M0L8new__bufS590;
  return 0;
}

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS558
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS558->$1;
}

int32_t _M0MPC15array5Array6lengthGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS559
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS559->$1;
}

int32_t _M0MPC15array5Array6lengthGiE(struct _M0TPB5ArrayGiE* _M0L4selfS560) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS560->$1;
}

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(
  int32_t _M0L8capacityS554
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS554 == 0) {
    moonbit_string_t* _M0L6_2atmpS2096 =
      (moonbit_string_t*)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGsE* _block_4542 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_4542)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
    _block_4542->$0 = _M0L6_2atmpS2096;
    _block_4542->$1 = 0;
    return _block_4542;
  } else {
    moonbit_string_t* _M0L6_2atmpS2097 =
      (moonbit_string_t*)moonbit_make_ref_array(_M0L8capacityS554, (moonbit_string_t)moonbit_string_literal_0.data);
    struct _M0TPB5ArrayGsE* _block_4543 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_4543)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
    _block_4543->$0 = _M0L6_2atmpS2097;
    _block_4543->$1 = 0;
    return _block_4543;
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(
  int32_t _M0L8capacityS555
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS555 == 0) {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2098 =
      (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _block_4544 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
    Moonbit_object_header(_block_4544)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 46, 0);
    _block_4544->$0 = _M0L6_2atmpS2098;
    _block_4544->$1 = 0;
    return _block_4544;
  } else {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2099 =
      (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array(_M0L8capacityS555, 0);
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _block_4545 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
    Moonbit_object_header(_block_4545)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 46, 0);
    _block_4545->$0 = _M0L6_2atmpS2099;
    _block_4545->$1 = 0;
    return _block_4545;
  }
}

struct _M0TPB5ArrayGiE* _M0MPC15array5Array11new_2einnerGiE(
  int32_t _M0L8capacityS556
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS556 == 0) {
    int32_t* _M0L6_2atmpS2100 = (int32_t*)moonbit_empty_int32_array;
    struct _M0TPB5ArrayGiE* _block_4546 =
      (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
    Moonbit_object_header(_block_4546)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
    _block_4546->$0 = _M0L6_2atmpS2100;
    _block_4546->$1 = 0;
    return _block_4546;
  } else {
    int32_t* _M0L6_2atmpS2101 =
      (int32_t*)moonbit_make_int32_array_raw(_M0L8capacityS556);
    struct _M0TPB5ArrayGiE* _block_4547 =
      (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
    Moonbit_object_header(_block_4547)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
    _block_4547->$0 = _M0L6_2atmpS2101;
    _block_4547->$1 = 0;
    return _block_4547;
  }
}

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(
  int32_t _M0L8capacityS557
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS557 == 0) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS2102 =
      (struct _M0TPC16string10StringView*)moonbit_empty_ref_valtype_array;
    struct _M0TPB5ArrayGRPC16string10StringViewE* _block_4548 =
      (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_block_4548)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 111, 0);
    _block_4548->$0 = _M0L6_2atmpS2102;
    _block_4548->$1 = 0;
    return _block_4548;
  } else {
    struct _M0TPC16string10StringView* _M0L6_2atmpS2103 =
      (struct _M0TPC16string10StringView*)moonbit_make_ref_valtype_array(_M0L8capacityS557, sizeof(struct _M0TPC16string10StringView), Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 114, 0), &(struct _M0TPC16string10StringView){.$0 = (moonbit_string_t)moonbit_string_literal_0.data, .$1 = 0, .$2 = 0});
    struct _M0TPB5ArrayGRPC16string10StringViewE* _block_4549 =
      (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_block_4549)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 111, 0);
    _block_4549->$0 = _M0L6_2atmpS2103;
    _block_4549->$1 = 0;
    return _block_4549;
  }
}

int32_t _M0MPC16string6String11has__prefix(
  moonbit_string_t _M0L4selfS552,
  struct _M0TPC16string10StringView _M0L3strS553
) {
  int32_t _M0L6_2atmpS2095;
  struct _M0TPC16string10StringView _M0L6_2atmpS2094;
  int32_t _result_4550;
  #line 297 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2095 = Moonbit_array_length(_M0L4selfS552);
  moonbit_incref(_M0L4selfS552);
  _M0L6_2atmpS2094
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS552, .$1 = 0, .$2 = _M0L6_2atmpS2095
  };
  #line 299 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4550
  = _M0MPC16string10StringView11has__prefix(_M0L6_2atmpS2094, _M0L3strS553);
  moonbit_decref(_M0L6_2atmpS2094.$0);
  return _result_4550;
}

int32_t _M0MPC16string10StringView11has__prefix(
  struct _M0TPC16string10StringView _M0L4selfS548,
  struct _M0TPC16string10StringView _M0L3strS549
) {
  int64_t _M0L7_2abindS547;
  #line 290 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 292 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS547
  = _M0MPC16string10StringView4find(_M0L4selfS548, _M0L3strS549);
  if (_M0L7_2abindS547 == 4294967296ll) {
    return 0;
  } else {
    int64_t _M0L7_2aSomeS550 = _M0L7_2abindS547;
    int32_t _M0L4_2aiS551 = (int32_t)_M0L7_2aSomeS550;
    return _M0L4_2aiS551 == 0;
  }
}

moonbit_string_t _M0MPC16string6String6repeat(
  moonbit_string_t _M0L4selfS540,
  int32_t _M0L1nS539
) {
  #line 1049 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  if (_M0L1nS539 < 0) {
    #line 1051 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_167.data);
  } else if (_M0L1nS539 == 0) {
    return (moonbit_string_t)moonbit_string_literal_0.data;
  } else if (_M0L1nS539 == 1) {
    moonbit_incref(_M0L4selfS540);
    return _M0L4selfS540;
  } else {
    int32_t _M0L3lenS541 = Moonbit_array_length(_M0L4selfS540);
    int32_t _M0L5totalS542 = _M0L3lenS541 * _M0L1nS539;
    int32_t _if__result_4551;
    if (_M0L3lenS541 == 0) {
      _if__result_4551 = 1;
    } else {
      int32_t _M0L6_2atmpS2092 = _M0L5totalS542 / _M0L1nS539;
      _if__result_4551 = _M0L6_2atmpS2092 == _M0L3lenS541;
    }
    if (_if__result_4551) {
      struct _M0TPB13StringBuilder* _M0L3bufS543;
      moonbit_string_t _M0L3strS544;
      int32_t _M0L2__S545;
      moonbit_string_t _result_4553;
      #line 1060 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3bufS543
      = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L5totalS542);
      #line 1061 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3strS544 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS540);
      _M0L2__S545 = 0;
      while (1) {
        if (_M0L2__S545 < _M0L1nS539) {
          int32_t _M0L6_2atmpS2093;
          #line 1063 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS543, _M0L3strS544);
          _M0L6_2atmpS2093 = _M0L2__S545 + 1;
          _M0L2__S545 = _M0L6_2atmpS2093;
          continue;
        } else {
          moonbit_decref(_M0L3strS544);
        }
        break;
      }
      #line 1065 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _result_4553 = _M0MPB13StringBuilder10to__string(_M0L3bufS543);
      moonbit_decref(_M0L3bufS543);
      return _result_4553;
    } else {
      #line 1058 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_168.data);
    }
  }
}

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(
  moonbit_string_t _M0L4selfS538
) {
  #line 222 "/home/developer/.moon/lib/core/builtin/show.mbt"
  moonbit_incref(_M0L4selfS538);
  return _M0L4selfS538;
}

int64_t _M0MPC16string6String4find(
  moonbit_string_t _M0L4selfS536,
  struct _M0TPC16string10StringView _M0L3strS537
) {
  int32_t _M0L6_2atmpS2091;
  struct _M0TPC16string10StringView _M0L6_2atmpS2090;
  int64_t _result_4554;
  #line 101 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2091 = Moonbit_array_length(_M0L4selfS536);
  moonbit_incref(_M0L4selfS536);
  _M0L6_2atmpS2090
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS536, .$1 = 0, .$2 = _M0L6_2atmpS2091
  };
  #line 102 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4554
  = _M0MPC16string10StringView4find(_M0L6_2atmpS2090, _M0L3strS537);
  moonbit_decref(_M0L6_2atmpS2090.$0);
  return _result_4554;
}

struct moonbit_result_0 _M0FPB12assert__true(
  int32_t _M0L1xS529,
  void* _M0L3msgS531,
  moonbit_string_t _M0L3locS535
) {
  #line 123 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  if (!_M0L1xS529) {
    struct _M0TPC16string10StringView _M0L9fail__msgS530;
    struct moonbit_result_0 _result_4555;
    switch (Moonbit_object_tag(_M0L3msgS531)) {
      case 1: {
        struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS532 =
          (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L3msgS531;
        struct _M0TPC16string10StringView _M0L8_2afieldS3862 =
          _M0L7_2aSomeS532->$0;
        moonbit_incref(_M0L8_2afieldS3862.$0);
        _M0L9fail__msgS530 = _M0L8_2afieldS3862;
        break;
      }
      default: {
        struct _M0TPB13StringBuilder* _M0L18_2astring__builderS533;
        moonbit_string_t _M0L7_2abindS534;
        int32_t _M0L6_2atmpS2088;
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0L18_2astring__builderS533
        = _M0MPB13StringBuilder21StringBuilder_2einner(14);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS533, (moonbit_string_t)moonbit_string_literal_169.data);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0MPB13StringBuilder13write__objectGbE(_M0L18_2astring__builderS533, _M0L1xS529);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS533, (moonbit_string_t)moonbit_string_literal_170.data);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0L7_2abindS534
        = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS533);
        moonbit_decref(_M0L18_2astring__builderS533);
        _M0L6_2atmpS2088 = Moonbit_array_length(_M0L7_2abindS534);
        _M0L9fail__msgS530
        = (struct _M0TPC16string10StringView){
          .$0 = _M0L7_2abindS534, .$1 = 0, .$2 = _M0L6_2atmpS2088
        };
        break;
      }
    }
    #line 131 "/home/developer/.moon/lib/core/builtin/assert.mbt"
    _result_4555 = _M0FPB4failGuE(_M0L9fail__msgS530, _M0L3locS535);
    moonbit_decref(_M0L9fail__msgS530.$0);
    return _result_4555;
  } else {
    int32_t _M0L6_2atmpS2089 = 0;
    struct moonbit_result_0 _result_4556;
    _result_4556.tag = 1;
    _result_4556.data.ok = _M0L6_2atmpS2089;
    return _result_4556;
  }
}

int64_t _M0MPC16string10StringView4find(
  struct _M0TPC16string10StringView _M0L4selfS528,
  struct _M0TPC16string10StringView _M0L3strS527
) {
  int32_t _M0L3endS2086;
  int32_t _M0L5startS2087;
  int32_t _M0L6_2atmpS2085;
  #line 18 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2086 = _M0L3strS527.$2;
  _M0L5startS2087 = _M0L3strS527.$1;
  _M0L6_2atmpS2085 = _M0L3endS2086 - _M0L5startS2087;
  if (_M0L6_2atmpS2085 <= 4) {
    #line 20 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB18brute__force__find(_M0L4selfS528, _M0L3strS527);
  } else {
    #line 22 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB28boyer__moore__horspool__find(_M0L4selfS528, _M0L3strS527);
  }
}

int64_t _M0FPB18brute__force__find(
  struct _M0TPC16string10StringView _M0L8haystackS517,
  struct _M0TPC16string10StringView _M0L6needleS519
) {
  int32_t _M0L3endS2083;
  int32_t _M0L5startS2084;
  int32_t _M0L13haystack__lenS516;
  int32_t _M0L3endS2081;
  int32_t _M0L5startS2082;
  int32_t _M0L11needle__lenS518;
  #line 31 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2083 = _M0L8haystackS517.$2;
  _M0L5startS2084 = _M0L8haystackS517.$1;
  _M0L13haystack__lenS516 = _M0L3endS2083 - _M0L5startS2084;
  _M0L3endS2081 = _M0L6needleS519.$2;
  _M0L5startS2082 = _M0L6needleS519.$1;
  _M0L11needle__lenS518 = _M0L3endS2081 - _M0L5startS2082;
  if (_M0L11needle__lenS518 > 0) {
    if (_M0L13haystack__lenS516 >= _M0L11needle__lenS518) {
      moonbit_string_t _M0L3strS2079 = _M0L6needleS519.$0;
      int32_t _M0L5startS2080 = _M0L6needleS519.$1;
      int32_t _M0L13needle__firstS520 = _M0L3strS2079[_M0L5startS2080];
      int32_t _M0L12forward__lenS521 =
        _M0L13haystack__lenS516 - _M0L11needle__lenS518;
      int32_t _M0L1iS522 = 0;
      while (1) {
        if (_M0L1iS522 <= _M0L12forward__lenS521) {
          moonbit_string_t _M0L3strS2066 = _M0L8haystackS517.$0;
          int32_t _M0L5startS2068 = _M0L8haystackS517.$1;
          int32_t _M0L6_2atmpS2067 = _M0L5startS2068 + _M0L1iS522;
          int32_t _M0L6_2atmpS2065 = _M0L3strS2066[_M0L6_2atmpS2067];
          int32_t _M0L1jS525;
          int32_t _M0L6_2atmpS2064;
          if (_M0L6_2atmpS2065 != _M0L13needle__firstS520) {
            goto join_523;
          }
          _M0L1jS525 = 1;
          while (1) {
            if (_M0L1jS525 < _M0L11needle__lenS518) {
              moonbit_string_t _M0L3strS2074 = _M0L8haystackS517.$0;
              int32_t _M0L5startS2076 = _M0L8haystackS517.$1;
              int32_t _M0L6_2atmpS2077 = _M0L1iS522 + _M0L1jS525;
              int32_t _M0L6_2atmpS2075 = _M0L5startS2076 + _M0L6_2atmpS2077;
              int32_t _M0L6_2atmpS2069 = _M0L3strS2074[_M0L6_2atmpS2075];
              moonbit_string_t _M0L3strS2071 = _M0L6needleS519.$0;
              int32_t _M0L5startS2073 = _M0L6needleS519.$1;
              int32_t _M0L6_2atmpS2072 = _M0L5startS2073 + _M0L1jS525;
              int32_t _M0L6_2atmpS2070 = _M0L3strS2071[_M0L6_2atmpS2072];
              int32_t _M0L6_2atmpS2078;
              if (_M0L6_2atmpS2069 != _M0L6_2atmpS2070) {
                break;
              }
              _M0L6_2atmpS2078 = _M0L1jS525 + 1;
              _M0L1jS525 = _M0L6_2atmpS2078;
              continue;
            } else {
              return (int64_t)_M0L1iS522;
            }
            break;
          }
          goto join_523;
          goto joinlet_4558;
          join_523:;
          _M0L6_2atmpS2064 = _M0L1iS522 + 1;
          _M0L1iS522 = _M0L6_2atmpS2064;
          continue;
          joinlet_4558:;
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
  struct _M0TPC16string10StringView _M0L8haystackS504,
  struct _M0TPC16string10StringView _M0L6needleS506
) {
  int32_t _M0L3endS2062;
  int32_t _M0L5startS2063;
  int32_t _M0L13haystack__lenS503;
  int32_t _M0L3endS2060;
  int32_t _M0L5startS2061;
  int32_t _M0L11needle__lenS505;
  #line 57 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2062 = _M0L8haystackS504.$2;
  _M0L5startS2063 = _M0L8haystackS504.$1;
  _M0L13haystack__lenS503 = _M0L3endS2062 - _M0L5startS2063;
  _M0L3endS2060 = _M0L6needleS506.$2;
  _M0L5startS2061 = _M0L6needleS506.$1;
  _M0L11needle__lenS505 = _M0L3endS2060 - _M0L5startS2061;
  if (_M0L11needle__lenS505 > 0) {
    if (_M0L13haystack__lenS503 >= _M0L11needle__lenS505) {
      int32_t* _M0L11skip__tableS507 =
        (int32_t*)moonbit_make_int32_array(256, _M0L11needle__lenS505);
      int32_t _M0L7_2abindS508 = _M0L11needle__lenS505 - 1;
      int32_t _M0L1iS509 = 0;
      int32_t _M0L1iS511;
      while (1) {
        if (_M0L1iS509 < _M0L7_2abindS508) {
          moonbit_string_t _M0L3strS2035 = _M0L6needleS506.$0;
          int32_t _M0L5startS2037 = _M0L6needleS506.$1;
          int32_t _M0L6_2atmpS2036 = _M0L5startS2037 + _M0L1iS509;
          int32_t _M0L6_2atmpS2034 = _M0L3strS2035[_M0L6_2atmpS2036];
          int32_t _M0L6_2atmpS2033 = (int32_t)_M0L6_2atmpS2034;
          int32_t _M0L6_2atmpS2030 = _M0L6_2atmpS2033 & 255;
          int32_t _M0L6_2atmpS2032 = _M0L11needle__lenS505 - 1;
          int32_t _M0L6_2atmpS2031 = _M0L6_2atmpS2032 - _M0L1iS509;
          int32_t _M0L6_2atmpS2038;
          if (
            _M0L6_2atmpS2030 < 0
            || _M0L6_2atmpS2030
               >= Moonbit_array_length(_M0L11skip__tableS507)
          ) {
            #line 68 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L11skip__tableS507[_M0L6_2atmpS2030] = _M0L6_2atmpS2031;
          _M0L6_2atmpS2038 = _M0L1iS509 + 1;
          _M0L1iS509 = _M0L6_2atmpS2038;
          continue;
        }
        break;
      }
      _M0L1iS511 = 0;
      while (1) {
        int32_t _M0L6_2atmpS2039 =
          _M0L13haystack__lenS503 - _M0L11needle__lenS505;
        if (_M0L1iS511 <= _M0L6_2atmpS2039) {
          int32_t _M0L7_2abindS512 = _M0L11needle__lenS505 - 1;
          int32_t _M0L1jS513 = 0;
          moonbit_string_t _M0L3strS2055;
          int32_t _M0L5startS2057;
          int32_t _M0L6_2atmpS2059;
          int32_t _M0L6_2atmpS2058;
          int32_t _M0L6_2atmpS2056;
          int32_t _M0L6_2atmpS2054;
          int32_t _M0L6_2atmpS2053;
          int32_t _M0L6_2atmpS2052;
          int32_t _M0L6_2atmpS2051;
          int32_t _M0L6_2atmpS2050;
          while (1) {
            if (_M0L1jS513 <= _M0L7_2abindS512) {
              moonbit_string_t _M0L3strS2045 = _M0L8haystackS504.$0;
              int32_t _M0L5startS2047 = _M0L8haystackS504.$1;
              int32_t _M0L6_2atmpS2048 = _M0L1iS511 + _M0L1jS513;
              int32_t _M0L6_2atmpS2046 = _M0L5startS2047 + _M0L6_2atmpS2048;
              int32_t _M0L6_2atmpS2040 = _M0L3strS2045[_M0L6_2atmpS2046];
              moonbit_string_t _M0L3strS2042 = _M0L6needleS506.$0;
              int32_t _M0L5startS2044 = _M0L6needleS506.$1;
              int32_t _M0L6_2atmpS2043 = _M0L5startS2044 + _M0L1jS513;
              int32_t _M0L6_2atmpS2041 = _M0L3strS2042[_M0L6_2atmpS2043];
              int32_t _M0L6_2atmpS2049;
              if (_M0L6_2atmpS2040 != _M0L6_2atmpS2041) {
                break;
              }
              _M0L6_2atmpS2049 = _M0L1jS513 + 1;
              _M0L1jS513 = _M0L6_2atmpS2049;
              continue;
            } else {
              moonbit_decref(_M0L11skip__tableS507);
              return (int64_t)_M0L1iS511;
            }
            break;
          }
          _M0L3strS2055 = _M0L8haystackS504.$0;
          _M0L5startS2057 = _M0L8haystackS504.$1;
          _M0L6_2atmpS2059 = _M0L1iS511 + _M0L11needle__lenS505;
          _M0L6_2atmpS2058 = _M0L6_2atmpS2059 - 1;
          _M0L6_2atmpS2056 = _M0L5startS2057 + _M0L6_2atmpS2058;
          _M0L6_2atmpS2054 = _M0L3strS2055[_M0L6_2atmpS2056];
          _M0L6_2atmpS2053 = (int32_t)_M0L6_2atmpS2054;
          _M0L6_2atmpS2052 = _M0L6_2atmpS2053 & 255;
          if (
            _M0L6_2atmpS2052 < 0
            || _M0L6_2atmpS2052
               >= Moonbit_array_length(_M0L11skip__tableS507)
          ) {
            #line 73 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2051 = (int32_t)_M0L11skip__tableS507[_M0L6_2atmpS2052];
          _M0L6_2atmpS2050 = _M0L1iS511 + _M0L6_2atmpS2051;
          _M0L1iS511 = _M0L6_2atmpS2050;
          continue;
        } else {
          moonbit_decref(_M0L11skip__tableS507);
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
  struct _M0TPB13StringBuilder* _M0L4selfS502,
  struct _M0TPC16string10StringView _M0L3strS501
) {
  int32_t _M0L3endS2028;
  int32_t _M0L5startS2029;
  int32_t _M0L8str__lenS500;
  int32_t _M0L3lenS2021;
  int32_t _M0L6_2atmpS2020;
  uint16_t* _M0L4dataS2022;
  int32_t _M0L3lenS2023;
  moonbit_string_t _M0L6_2atmpS2024;
  int32_t _M0L6_2atmpS2025;
  int32_t _M0L3lenS2027;
  int32_t _M0L6_2atmpS2026;
  #line 131 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3endS2028 = _M0L3strS501.$2;
  _M0L5startS2029 = _M0L3strS501.$1;
  _M0L8str__lenS500 = _M0L3endS2028 - _M0L5startS2029;
  _M0L3lenS2021 = _M0L4selfS502->$1;
  _M0L6_2atmpS2020 = _M0L3lenS2021 + _M0L8str__lenS500;
  #line 136 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS502, _M0L6_2atmpS2020);
  _M0L4dataS2022 = _M0L4selfS502->$0;
  _M0L3lenS2023 = _M0L4selfS502->$1;
  moonbit_incref(_M0L4dataS2022);
  #line 139 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS2024 = _M0MPC16string10StringView4data(_M0L3strS501);
  #line 140 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS2025 = _M0MPC16string10StringView13start__offset(_M0L3strS501);
  #line 137 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS2022, _M0L3lenS2023, _M0L6_2atmpS2024, _M0L6_2atmpS2025, _M0L8str__lenS500);
  moonbit_decref(_M0L4dataS2022);
  moonbit_decref(_M0L6_2atmpS2024);
  _M0L3lenS2027 = _M0L4selfS502->$1;
  _M0L6_2atmpS2026 = _M0L3lenS2027 + _M0L8str__lenS500;
  _M0L4selfS502->$1 = _M0L6_2atmpS2026;
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t _M0L4selfS498,
  int32_t _M0L13start__offsetS499,
  int64_t _M0L11end__offsetS496
) {
  int32_t _M0L11end__offsetS495;
  int32_t _if__result_4563;
  #line 614 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS496 == 4294967296ll) {
    _M0L11end__offsetS495 = Moonbit_array_length(_M0L4selfS498);
  } else {
    int64_t _M0L7_2aSomeS497 = _M0L11end__offsetS496;
    _M0L11end__offsetS495 = (int32_t)_M0L7_2aSomeS497;
  }
  if (_M0L13start__offsetS499 >= 0) {
    if (_M0L13start__offsetS499 <= _M0L11end__offsetS495) {
      int32_t _M0L6_2atmpS2019 = Moonbit_array_length(_M0L4selfS498);
      _if__result_4563 = _M0L11end__offsetS495 <= _M0L6_2atmpS2019;
    } else {
      _if__result_4563 = 0;
    }
  } else {
    _if__result_4563 = 0;
  }
  if (_if__result_4563) {
    moonbit_incref(_M0L4selfS498);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS498,
                                                 .$1 = _M0L13start__offsetS499,
                                                 .$2 = _M0L11end__offsetS495};
  } else {
    #line 623 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_171.data);
  }
}

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView _M0L4selfS490
) {
  int32_t _M0L5startS489;
  int32_t _M0L3endS491;
  struct _M0TPB8MutLocalGiE* _M0L5indexS492;
  struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__* _closure_4564;
  struct _M0TWEOc* _M0L6_2atmpS1998;
  struct _M0TPB4IterGcE* _result_4565;
  #line 214 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L5startS489 = _M0L4selfS490.$1;
  _M0L3endS491 = _M0L4selfS490.$2;
  _M0L5indexS492
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS492)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS492->$0 = _M0L5startS489;
  moonbit_incref(_M0L4selfS490.$0);
  _closure_4564
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__*)moonbit_malloc(sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__));
  Moonbit_object_header(_closure_4564)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 134, 0);
  _closure_4564->code = &_M0MPC16string10StringView4iterC1999l219;
  _closure_4564->$0 = _M0L5indexS492;
  _closure_4564->$1 = _M0L3endS491;
  _closure_4564->$2 = _M0L4selfS490;
  _M0L6_2atmpS1998 = (struct _M0TWEOc*)_closure_4564;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_4565 = _M0MPB4Iter3newGcE(_M0L6_2atmpS1998, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1998);
  return _result_4565;
}

int32_t _M0MPC16string10StringView4iterC1999l219(
  struct _M0TWEOc* _M0L6_2aenvS2000
) {
  struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__* _M0L14_2acasted__envS2001;
  struct _M0TPC16string10StringView _M0L4selfS490;
  int32_t _M0L3endS491;
  struct _M0TPB8MutLocalGiE* _M0L5indexS492;
  int32_t _M0L3valS2002;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L14_2acasted__envS2001
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1999__l219__*)_M0L6_2aenvS2000;
  _M0L4selfS490 = _M0L14_2acasted__envS2001->$2;
  _M0L3endS491 = _M0L14_2acasted__envS2001->$1;
  _M0L5indexS492 = _M0L14_2acasted__envS2001->$0;
  _M0L3valS2002 = _M0L5indexS492->$0;
  if (_M0L3valS2002 < _M0L3endS491) {
    moonbit_string_t _M0L3strS2017 = _M0L4selfS490.$0;
    int32_t _M0L3valS2018 = _M0L5indexS492->$0;
    int32_t _M0L2c1S493 = _M0L3strS2017[_M0L3valS2018];
    int32_t _if__result_4566;
    int32_t _M0L3valS2015;
    int32_t _M0L6_2atmpS2014;
    int32_t _M0L6_2atmpS2016;
    #line 222 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S493)) {
      int32_t _M0L3valS2005 = _M0L5indexS492->$0;
      int32_t _M0L6_2atmpS2003 = _M0L3valS2005 + 1;
      int32_t _M0L3endS2004 = _M0L4selfS490.$2;
      _if__result_4566 = _M0L6_2atmpS2003 < _M0L3endS2004;
    } else {
      _if__result_4566 = 0;
    }
    if (_if__result_4566) {
      moonbit_string_t _M0L3strS2011 = _M0L4selfS490.$0;
      int32_t _M0L3valS2013 = _M0L5indexS492->$0;
      int32_t _M0L6_2atmpS2012 = _M0L3valS2013 + 1;
      int32_t _M0L2c2S494 = _M0L3strS2011[_M0L6_2atmpS2012];
      #line 224 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S494)) {
        int32_t _M0L3valS2007 = _M0L5indexS492->$0;
        int32_t _M0L6_2atmpS2006 = _M0L3valS2007 + 2;
        int32_t _M0L6_2atmpS2009;
        int32_t _M0L6_2atmpS2010;
        int32_t _M0L6_2atmpS2008;
        _M0L5indexS492->$0 = _M0L6_2atmpS2006;
        _M0L6_2atmpS2009 = (int32_t)_M0L2c1S493;
        _M0L6_2atmpS2010 = (int32_t)_M0L2c2S494;
        #line 226 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        _M0L6_2atmpS2008
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS2009, _M0L6_2atmpS2010);
        return _M0L6_2atmpS2008;
      }
    }
    _M0L3valS2015 = _M0L5indexS492->$0;
    _M0L6_2atmpS2014 = _M0L3valS2015 + 1;
    _M0L5indexS492->$0 = _M0L6_2atmpS2014;
    #line 230 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    _M0L6_2atmpS2016 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S493);
    return _M0L6_2atmpS2016;
  } else {
    return -1;
  }
}

int32_t _M0IPC16string10StringViewPB4Show6output(
  struct _M0TPC16string10StringView _M0L4selfS488,
  struct _M0TPB6Logger _M0L6loggerS487
) {
  #line 203 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  #line 204 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L6loggerS487.$0->$method_2(_M0L6loggerS487.$1, _M0L4selfS488);
  return 0;
}

moonbit_string_t _M0MPC16string10StringView9to__owned(
  struct _M0TPC16string10StringView _M0L4selfS486
) {
  moonbit_string_t _M0L3strS1995;
  int32_t _M0L5startS1996;
  int32_t _M0L3endS1997;
  moonbit_string_t _result_4567;
  #line 196 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1995 = _M0L4selfS486.$0;
  _M0L5startS1996 = _M0L4selfS486.$1;
  _M0L3endS1997 = _M0L4selfS486.$2;
  moonbit_incref(_M0L3strS1995);
  #line 199 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_4567
  = _M0MPC16string6String17unsafe__substring(_M0L3strS1995, _M0L5startS1996, _M0L3endS1997);
  moonbit_decref(_M0L3strS1995);
  return _result_4567;
}

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t _M0L3strS483,
  int32_t _M0L5startS481,
  int32_t _M0L3endS482
) {
  int32_t _if__result_4568;
  int32_t _M0L3lenS484;
  int32_t _M0L6_2atmpS1993;
  int32_t _M0L6_2atmpS1994;
  moonbit_bytes_t _M0L5bytesS485;
  moonbit_bytes_t _M0L6_2atmpS1992;
  moonbit_string_t _result_4569;
  #line 91 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L5startS481 == 0) {
    int32_t _M0L6_2atmpS1991 = Moonbit_array_length(_M0L3strS483);
    _if__result_4568 = _M0L3endS482 == _M0L6_2atmpS1991;
  } else {
    _if__result_4568 = 0;
  }
  if (_if__result_4568) {
    moonbit_incref(_M0L3strS483);
    return _M0L3strS483;
  }
  _M0L3lenS484 = _M0L3endS482 - _M0L5startS481;
  _M0L6_2atmpS1993 = _M0L3lenS484 * 2;
  #line 101 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1994 = _M0IPC14byte4BytePB7Default7default();
  _M0L5bytesS485
  = (moonbit_bytes_t)moonbit_make_bytes(_M0L6_2atmpS1993, _M0L6_2atmpS1994);
  #line 102 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0MPC15array10FixedArray18blit__from__string(_M0L5bytesS485, 0, _M0L3strS483, _M0L5startS481, _M0L3lenS484);
  _M0L6_2atmpS1992 = _M0L5bytesS485;
  #line 103 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_4569
  = _M0MPC15bytes5Bytes29to__unchecked__string_2einner(_M0L6_2atmpS1992, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1992);
  return _result_4569;
}

int32_t _M0IPC14byte4BytePB7Default7default() {
  #line 231 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  return 0;
}

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t _M0L4selfS476,
  int32_t _M0L6offsetS480,
  int64_t _M0L6lengthS478
) {
  int32_t _M0L3lenS475;
  int32_t _M0L6lengthS477;
  int32_t _if__result_4570;
  #line 77 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L3lenS475 = Moonbit_array_length(_M0L4selfS476);
  if (_M0L6lengthS478 == 4294967296ll) {
    _M0L6lengthS477 = _M0L3lenS475 - _M0L6offsetS480;
  } else {
    int64_t _M0L7_2aSomeS479 = _M0L6lengthS478;
    _M0L6lengthS477 = (int32_t)_M0L7_2aSomeS479;
  }
  if (_M0L6offsetS480 >= 0) {
    if (_M0L6lengthS477 >= 0) {
      int32_t _M0L6_2atmpS1990 = _M0L6offsetS480 + _M0L6lengthS477;
      _if__result_4570 = _M0L6_2atmpS1990 <= _M0L3lenS475;
    } else {
      _if__result_4570 = 0;
    }
  } else {
    _if__result_4570 = 0;
  }
  if (_if__result_4570) {
    moonbit_incref(_M0L4selfS476);
    #line 85 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    return _M0FPB19unsafe__sub__string(_M0L4selfS476, _M0L6offsetS480, _M0L6lengthS477);
  } else {
    #line 84 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t _M0L4selfS467,
  int32_t _M0L13bytes__offsetS462,
  moonbit_string_t _M0L3strS469,
  int32_t _M0L11str__offsetS465,
  int32_t _M0L6lengthS463
) {
  int32_t _M0L6_2atmpS1989;
  int32_t _M0L6_2atmpS1988;
  int32_t _M0L2e1S461;
  int32_t _M0L6_2atmpS1987;
  int32_t _M0L2e2S464;
  int32_t _M0L4len1S466;
  int32_t _M0L4len2S468;
  int32_t _if__result_4571;
  #line 125 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L6_2atmpS1989 = _M0L6lengthS463 * 2;
  _M0L6_2atmpS1988 = _M0L13bytes__offsetS462 + _M0L6_2atmpS1989;
  _M0L2e1S461 = _M0L6_2atmpS1988 - 1;
  _M0L6_2atmpS1987 = _M0L11str__offsetS465 + _M0L6lengthS463;
  _M0L2e2S464 = _M0L6_2atmpS1987 - 1;
  _M0L4len1S466 = Moonbit_array_length(_M0L4selfS467);
  _M0L4len2S468 = Moonbit_array_length(_M0L3strS469);
  if (_M0L6lengthS463 >= 0) {
    if (_M0L13bytes__offsetS462 >= 0) {
      if (_M0L2e1S461 < _M0L4len1S466) {
        if (_M0L11str__offsetS465 >= 0) {
          _if__result_4571 = _M0L2e2S464 < _M0L4len2S468;
        } else {
          _if__result_4571 = 0;
        }
      } else {
        _if__result_4571 = 0;
      }
    } else {
      _if__result_4571 = 0;
    }
  } else {
    _if__result_4571 = 0;
  }
  if (_if__result_4571) {
    int32_t _M0L16end__str__offsetS470 =
      _M0L11str__offsetS465 + _M0L6lengthS463;
    int32_t _M0L1iS471 = _M0L11str__offsetS465;
    int32_t _M0L1jS472 = _M0L13bytes__offsetS462;
    while (1) {
      if (_M0L1iS471 < _M0L16end__str__offsetS470) {
        int32_t _M0L6_2atmpS1984 = _M0L3strS469[_M0L1iS471];
        int32_t _M0L6_2atmpS1983 = (int32_t)_M0L6_2atmpS1984;
        uint32_t _M0L1cS473 = *(uint32_t*)&_M0L6_2atmpS1983;
        uint32_t _M0L6_2atmpS1979 = _M0L1cS473 & 255u;
        int32_t _M0L6_2atmpS1978;
        int32_t _M0L6_2atmpS1980;
        uint32_t _M0L6_2atmpS1982;
        int32_t _M0L6_2atmpS1981;
        int32_t _M0L6_2atmpS1985;
        int32_t _M0L6_2atmpS1986;
        #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1978 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1979);
        if (
          _M0L1jS472 < 0 || _M0L1jS472 >= Moonbit_array_length(_M0L4selfS467)
        ) {
          #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS467[_M0L1jS472] = _M0L6_2atmpS1978;
        _M0L6_2atmpS1980 = _M0L1jS472 + 1;
        _M0L6_2atmpS1982 = _M0L1cS473 >> 8;
        #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1981 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1982);
        if (
          _M0L6_2atmpS1980 < 0
          || _M0L6_2atmpS1980 >= Moonbit_array_length(_M0L4selfS467)
        ) {
          #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS467[_M0L6_2atmpS1980] = _M0L6_2atmpS1981;
        _M0L6_2atmpS1985 = _M0L1iS471 + 1;
        _M0L6_2atmpS1986 = _M0L1jS472 + 2;
        _M0L1iS471 = _M0L6_2atmpS1985;
        _M0L1jS472 = _M0L6_2atmpS1986;
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

int32_t _M0MPC14uint4UInt8to__byte(uint32_t _M0L4selfS460) {
  int32_t _M0L6_2atmpS1977;
  #line 2519 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1977 = *(int32_t*)&_M0L4selfS460;
  return _M0L6_2atmpS1977 & 0xff;
}

moonbit_string_t* _M0MPC15array5Array6bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS455
) {
  moonbit_string_t* _M0L8_2afieldS3875;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3875 = _M0L4selfS455->$0;
  moonbit_incref(_M0L8_2afieldS3875);
  return _M0L8_2afieldS3875;
}

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS456
) {
  struct _M0TUsiE** _M0L8_2afieldS3876;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3876 = _M0L4selfS456->$0;
  moonbit_incref(_M0L8_2afieldS3876);
  return _M0L8_2afieldS3876;
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS457
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8_2afieldS3877;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3877 = _M0L4selfS457->$0;
  moonbit_incref(_M0L8_2afieldS3877);
  return _M0L8_2afieldS3877;
}

struct _M0TPC16string10StringView* _M0MPC15array5Array6bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS458
) {
  struct _M0TPC16string10StringView* _M0L8_2afieldS3878;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3878 = _M0L4selfS458->$0;
  moonbit_incref(_M0L8_2afieldS3878);
  return _M0L8_2afieldS3878;
}

int32_t* _M0MPC15array5Array6bufferGiE(struct _M0TPB5ArrayGiE* _M0L4selfS459) {
  int32_t* _M0L8_2afieldS3879;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3879 = _M0L4selfS459->$0;
  moonbit_incref(_M0L8_2afieldS3879);
  return _M0L8_2afieldS3879;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView _M0L4selfS453,
  int32_t _M0L13start__offsetS454,
  int64_t _M0L11end__offsetS451
) {
  int32_t _M0L11end__offsetS450;
  int32_t _if__result_4573;
  #line 105 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS451 == 4294967296ll) {
    int32_t _M0L3endS1975 = _M0L4selfS453.$2;
    int32_t _M0L5startS1976 = _M0L4selfS453.$1;
    _M0L11end__offsetS450 = _M0L3endS1975 - _M0L5startS1976;
  } else {
    int64_t _M0L7_2aSomeS452 = _M0L11end__offsetS451;
    _M0L11end__offsetS450 = (int32_t)_M0L7_2aSomeS452;
  }
  if (_M0L13start__offsetS454 >= 0) {
    if (_M0L13start__offsetS454 <= _M0L11end__offsetS450) {
      int32_t _M0L3endS1968 = _M0L4selfS453.$2;
      int32_t _M0L5startS1969 = _M0L4selfS453.$1;
      int32_t _M0L6_2atmpS1967 = _M0L3endS1968 - _M0L5startS1969;
      _if__result_4573 = _M0L11end__offsetS450 <= _M0L6_2atmpS1967;
    } else {
      _if__result_4573 = 0;
    }
  } else {
    _if__result_4573 = 0;
  }
  if (_if__result_4573) {
    moonbit_string_t _M0L3strS1970 = _M0L4selfS453.$0;
    int32_t _M0L5startS1974 = _M0L4selfS453.$1;
    int32_t _M0L6_2atmpS1971 = _M0L5startS1974 + _M0L13start__offsetS454;
    int32_t _M0L5startS1973 = _M0L4selfS453.$1;
    int32_t _M0L6_2atmpS1972 = _M0L5startS1973 + _M0L11end__offsetS450;
    moonbit_incref(_M0L3strS1970);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1970,
                                                 .$1 = _M0L6_2atmpS1971,
                                                 .$2 = _M0L6_2atmpS1972};
  } else {
    #line 114 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_171.data);
  }
}

struct _M0TPB4IterGsE* _M0MPB4Iter3newGsE(
  struct _M0TWEOs* _M0L1fS434,
  int64_t _M0L10size__hintS431
) {
  int64_t _M0L10size__hintS430;
  struct _M0TPB4IterGsE* _block_4574;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS431 == 4294967296ll) {
    _M0L10size__hintS430 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS432 = _M0L10size__hintS431;
    int32_t _M0L4_2anS433 = (int32_t)_M0L7_2aSomeS432;
    if (_M0L4_2anS433 > 0) {
      _M0L10size__hintS430 = (int64_t)_M0L4_2anS433;
    } else {
      _M0L10size__hintS430 = _M0MPB4Iter3newN6constrS9988GsE;
    }
  }
  moonbit_incref(_M0L1fS434);
  _block_4574
  = (struct _M0TPB4IterGsE*)moonbit_malloc(sizeof(struct _M0TPB4IterGsE));
  Moonbit_object_header(_block_4574)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 138, 0);
  _block_4574->$0 = _M0L1fS434;
  _block_4574->$1 = _M0L10size__hintS430;
  return _block_4574;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L1fS439,
  int64_t _M0L10size__hintS436
) {
  int64_t _M0L10size__hintS435;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_4575;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS436 == 4294967296ll) {
    _M0L10size__hintS435 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS437 = _M0L10size__hintS436;
    int32_t _M0L4_2anS438 = (int32_t)_M0L7_2aSomeS437;
    if (_M0L4_2anS438 > 0) {
      _M0L10size__hintS435 = (int64_t)_M0L4_2anS438;
    } else {
      _M0L10size__hintS435
      = _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE;
    }
  }
  moonbit_incref(_M0L1fS439);
  _block_4575
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_4575)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 131, 0);
  _block_4575->$0 = _M0L1fS439;
  _block_4575->$1 = _M0L10size__hintS435;
  return _block_4575;
}

struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0MPB4Iter3newGUsRPB5ArrayGiEEE(
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L1fS444,
  int64_t _M0L10size__hintS441
) {
  int64_t _M0L10size__hintS440;
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _block_4576;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS441 == 4294967296ll) {
    _M0L10size__hintS440 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS442 = _M0L10size__hintS441;
    int32_t _M0L4_2anS443 = (int32_t)_M0L7_2aSomeS442;
    if (_M0L4_2anS443 > 0) {
      _M0L10size__hintS440 = (int64_t)_M0L4_2anS443;
    } else {
      _M0L10size__hintS440 = _M0MPB4Iter3newN6constrS9988GUsRPB5ArrayGiEEE;
    }
  }
  moonbit_incref(_M0L1fS444);
  _block_4576
  = (struct _M0TPB4IterGUsRPB5ArrayGiEEE*)moonbit_malloc(sizeof(struct _M0TPB4IterGUsRPB5ArrayGiEEE));
  Moonbit_object_header(_block_4576)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 141, 0);
  _block_4576->$0 = _M0L1fS444;
  _block_4576->$1 = _M0L10size__hintS440;
  return _block_4576;
}

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(
  struct _M0TWEOc* _M0L1fS449,
  int64_t _M0L10size__hintS446
) {
  int64_t _M0L10size__hintS445;
  struct _M0TPB4IterGcE* _block_4577;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS446 == 4294967296ll) {
    _M0L10size__hintS445 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS447 = _M0L10size__hintS446;
    int32_t _M0L4_2anS448 = (int32_t)_M0L7_2aSomeS447;
    if (_M0L4_2anS448 > 0) {
      _M0L10size__hintS445 = (int64_t)_M0L4_2anS448;
    } else {
      _M0L10size__hintS445 = _M0MPB4Iter3newN6constrS9988GcE;
    }
  }
  moonbit_incref(_M0L1fS449);
  _block_4577
  = (struct _M0TPB4IterGcE*)moonbit_malloc(sizeof(struct _M0TPB4IterGcE));
  Moonbit_object_header(_block_4577)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 144, 0);
  _block_4577->$0 = _M0L1fS449;
  _block_4577->$1 = _M0L10size__hintS445;
  return _block_4577;
}

struct moonbit_result_0 _M0FPB10assert__eqGiE(
  int32_t _M0L1aS416,
  int32_t _M0L1bS417,
  moonbit_string_t _M0L3msgS419,
  moonbit_string_t _M0L3locS422
) {
  #line 45 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  if (_M0L1aS416 != _M0L1bS417) {
    moonbit_string_t _M0L9fail__msgS418;
    int32_t _M0L6_2atmpS1958;
    struct _M0TPC16string10StringView _M0L6_2atmpS1957;
    struct moonbit_result_0 _result_4578;
    if (_M0L3msgS419 == 0) {
      struct _M0TPB13StringBuilder* _M0L18_2astring__builderS421;
      moonbit_string_t _M0L6_2atmpS1959;
      moonbit_string_t _M0L6_2atmpS1960;
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L18_2astring__builderS421
      = _M0MPB13StringBuilder21StringBuilder_2einner(6);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS421, (moonbit_string_t)moonbit_string_literal_169.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L6_2atmpS1959 = _M0FPB13debug__stringGiE(_M0L1aS416);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS421, _M0L6_2atmpS1959);
      moonbit_decref(_M0L6_2atmpS1959);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS421, (moonbit_string_t)moonbit_string_literal_172.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L6_2atmpS1960 = _M0FPB13debug__stringGiE(_M0L1bS417);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS421, _M0L6_2atmpS1960);
      moonbit_decref(_M0L6_2atmpS1960);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS421, (moonbit_string_t)moonbit_string_literal_169.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L9fail__msgS418
      = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS421);
      moonbit_decref(_M0L18_2astring__builderS421);
    } else {
      moonbit_string_t _M0L7_2aSomeS420 = _M0L3msgS419;
      if (_M0L7_2aSomeS420) {
        moonbit_incref(_M0L7_2aSomeS420);
      }
      _M0L9fail__msgS418 = _M0L7_2aSomeS420;
    }
    _M0L6_2atmpS1958 = Moonbit_array_length(_M0L9fail__msgS418);
    _M0L6_2atmpS1957
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L9fail__msgS418, .$1 = 0, .$2 = _M0L6_2atmpS1958
    };
    #line 59 "/home/developer/.moon/lib/core/builtin/assert.mbt"
    _result_4578 = _M0FPB4failGuE(_M0L6_2atmpS1957, _M0L3locS422);
    moonbit_decref(_M0L6_2atmpS1957.$0);
    return _result_4578;
  } else {
    int32_t _M0L6_2atmpS1961 = 0;
    struct moonbit_result_0 _result_4579;
    _result_4579.tag = 1;
    _result_4579.data.ok = _M0L6_2atmpS1961;
    return _result_4579;
  }
}

struct moonbit_result_0 _M0FPB10assert__eqGsE(
  moonbit_string_t _M0L1aS423,
  moonbit_string_t _M0L1bS424,
  moonbit_string_t _M0L3msgS426,
  moonbit_string_t _M0L3locS429
) {
  #line 45 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  #line 54 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  if (_M0IP016_24default__implPB2Eq10not__equalGsE(_M0L1aS423, _M0L1bS424)) {
    moonbit_string_t _M0L9fail__msgS425;
    int32_t _M0L6_2atmpS1963;
    struct _M0TPC16string10StringView _M0L6_2atmpS1962;
    struct moonbit_result_0 _result_4580;
    if (_M0L3msgS426 == 0) {
      struct _M0TPB13StringBuilder* _M0L18_2astring__builderS428;
      moonbit_string_t _M0L6_2atmpS1964;
      moonbit_string_t _M0L6_2atmpS1965;
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L18_2astring__builderS428
      = _M0MPB13StringBuilder21StringBuilder_2einner(6);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS428, (moonbit_string_t)moonbit_string_literal_169.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L6_2atmpS1964 = _M0FPB13debug__stringGsE(_M0L1aS423);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS428, _M0L6_2atmpS1964);
      moonbit_decref(_M0L6_2atmpS1964);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS428, (moonbit_string_t)moonbit_string_literal_172.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L6_2atmpS1965 = _M0FPB13debug__stringGsE(_M0L1bS424);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS428, _M0L6_2atmpS1965);
      moonbit_decref(_M0L6_2atmpS1965);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS428, (moonbit_string_t)moonbit_string_literal_169.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L9fail__msgS425
      = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS428);
      moonbit_decref(_M0L18_2astring__builderS428);
    } else {
      moonbit_string_t _M0L7_2aSomeS427 = _M0L3msgS426;
      if (_M0L7_2aSomeS427) {
        moonbit_incref(_M0L7_2aSomeS427);
      }
      _M0L9fail__msgS425 = _M0L7_2aSomeS427;
    }
    _M0L6_2atmpS1963 = Moonbit_array_length(_M0L9fail__msgS425);
    _M0L6_2atmpS1962
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L9fail__msgS425, .$1 = 0, .$2 = _M0L6_2atmpS1963
    };
    #line 59 "/home/developer/.moon/lib/core/builtin/assert.mbt"
    _result_4580 = _M0FPB4failGuE(_M0L6_2atmpS1962, _M0L3locS429);
    moonbit_decref(_M0L6_2atmpS1962.$0);
    return _result_4580;
  } else {
    int32_t _M0L6_2atmpS1966 = 0;
    struct moonbit_result_0 _result_4581;
    _result_4581.tag = 1;
    _result_4581.data.ok = _M0L6_2atmpS1966;
    return _result_4581;
  }
}

struct moonbit_result_0 _M0FPB4failGuE(
  struct _M0TPC16string10StringView _M0L3msgS415,
  moonbit_string_t _M0L3locS414
) {
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS413;
  moonbit_string_t _M0L6_2atmpS1956;
  void* _M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS1955;
  struct moonbit_result_0 _result_4582;
  #line 56 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L18_2astring__builderS413
  = _M0MPB13StringBuilder21StringBuilder_2einner(9);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB13StringBuilder13write__objectGRPB9SourceLocE(_M0L18_2astring__builderS413, _M0L3locS414);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS413, (moonbit_string_t)moonbit_string_literal_173.data);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB13StringBuilder13write__objectGRPC16string10StringViewE(_M0L18_2astring__builderS413, _M0L3msgS415);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L6_2atmpS1956
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS413);
  moonbit_decref(_M0L18_2astring__builderS413);
  _M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS1955
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure));
  Moonbit_object_header(_M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS1955)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 147, 0);
  ((struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS1955)->$0
  = _M0L6_2atmpS1956;
  _result_4582.tag = 0;
  _result_4582.data.err
  = _M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS1955;
  return _result_4582;
}

moonbit_string_t _M0FPB13debug__stringGiE(int32_t _M0L1tS410) {
  struct _M0TPB13StringBuilder* _M0L3bufS409;
  struct _M0TPB6Logger _M0L6_2atmpS1953;
  moonbit_string_t _result_4583;
  #line 16 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  #line 17 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _M0L3bufS409 = _M0MPB13StringBuilder21StringBuilder_2einner(50);
  moonbit_incref(_M0L3bufS409);
  _M0L6_2atmpS1953
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS409
  };
  #line 18 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L1tS410, _M0L6_2atmpS1953);
  if (_M0L6_2atmpS1953.$1) {
    moonbit_decref(_M0L6_2atmpS1953.$1);
  }
  #line 19 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _result_4583 = _M0MPB13StringBuilder10to__string(_M0L3bufS409);
  moonbit_decref(_M0L3bufS409);
  return _result_4583;
}

moonbit_string_t _M0FPB13debug__stringGsE(moonbit_string_t _M0L1tS412) {
  struct _M0TPB13StringBuilder* _M0L3bufS411;
  struct _M0TPB6Logger _M0L6_2atmpS1954;
  moonbit_string_t _result_4584;
  #line 16 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  #line 17 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _M0L3bufS411 = _M0MPB13StringBuilder21StringBuilder_2einner(50);
  moonbit_incref(_M0L3bufS411);
  _M0L6_2atmpS1954
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS411
  };
  #line 18 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L1tS412, _M0L6_2atmpS1954);
  if (_M0L6_2atmpS1954.$1) {
    moonbit_decref(_M0L6_2atmpS1954.$1);
  }
  #line 19 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _result_4584 = _M0MPB13StringBuilder10to__string(_M0L3bufS411);
  moonbit_decref(_M0L3bufS411);
  return _result_4584;
}

moonbit_string_t _M0MPC13int3Int18to__string_2einner(
  int32_t _M0L4selfS393,
  int32_t _M0L5radixS392
) {
  int32_t _if__result_4585;
  int32_t _M0L12is__negativeS394;
  uint32_t _M0L3numS395;
  uint16_t* _M0L6bufferS396;
  #line 209 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS392 < 2) {
    _if__result_4585 = 1;
  } else {
    _if__result_4585 = _M0L5radixS392 > 36;
  }
  if (_if__result_4585) {
    #line 213 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_174.data);
  }
  if (_M0L4selfS393 == 0) {
    return (moonbit_string_t)moonbit_string_literal_175.data;
  }
  _M0L12is__negativeS394 = _M0L4selfS393 < 0;
  if (_M0L12is__negativeS394) {
    int32_t _M0L6_2atmpS1952 = -_M0L4selfS393;
    _M0L3numS395 = *(uint32_t*)&_M0L6_2atmpS1952;
  } else {
    _M0L3numS395 = *(uint32_t*)&_M0L4selfS393;
  }
  switch (_M0L5radixS392) {
    case 10: {
      int32_t _M0L10digit__lenS397;
      int32_t _M0L6_2atmpS1949;
      int32_t _M0L10total__lenS398;
      uint16_t* _M0L6bufferS399;
      int32_t _M0L12digit__startS400;
      #line 235 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS397 = _M0FPB12dec__count32(_M0L3numS395);
      if (_M0L12is__negativeS394) {
        _M0L6_2atmpS1949 = 1;
      } else {
        _M0L6_2atmpS1949 = 0;
      }
      _M0L10total__lenS398 = _M0L10digit__lenS397 + _M0L6_2atmpS1949;
      _M0L6bufferS399
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS398, 0);
      if (_M0L12is__negativeS394) {
        _M0L12digit__startS400 = 1;
      } else {
        _M0L12digit__startS400 = 0;
      }
      #line 239 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__dec(_M0L6bufferS399, _M0L3numS395, _M0L12digit__startS400, _M0L10total__lenS398);
      _M0L6bufferS396 = _M0L6bufferS399;
      break;
    }
    
    case 16: {
      int32_t _M0L10digit__lenS401;
      int32_t _M0L6_2atmpS1950;
      int32_t _M0L10total__lenS402;
      uint16_t* _M0L6bufferS403;
      int32_t _M0L12digit__startS404;
      #line 243 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS401 = _M0FPB12hex__count32(_M0L3numS395);
      if (_M0L12is__negativeS394) {
        _M0L6_2atmpS1950 = 1;
      } else {
        _M0L6_2atmpS1950 = 0;
      }
      _M0L10total__lenS402 = _M0L10digit__lenS401 + _M0L6_2atmpS1950;
      _M0L6bufferS403
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS402, 0);
      if (_M0L12is__negativeS394) {
        _M0L12digit__startS404 = 1;
      } else {
        _M0L12digit__startS404 = 0;
      }
      #line 247 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__hex(_M0L6bufferS403, _M0L3numS395, _M0L12digit__startS404, _M0L10total__lenS402);
      _M0L6bufferS396 = _M0L6bufferS403;
      break;
    }
    default: {
      int32_t _M0L10digit__lenS405;
      int32_t _M0L6_2atmpS1951;
      int32_t _M0L10total__lenS406;
      uint16_t* _M0L6bufferS407;
      int32_t _M0L12digit__startS408;
      #line 251 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS405
      = _M0FPB14radix__count32(_M0L3numS395, _M0L5radixS392);
      if (_M0L12is__negativeS394) {
        _M0L6_2atmpS1951 = 1;
      } else {
        _M0L6_2atmpS1951 = 0;
      }
      _M0L10total__lenS406 = _M0L10digit__lenS405 + _M0L6_2atmpS1951;
      _M0L6bufferS407
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS406, 0);
      if (_M0L12is__negativeS394) {
        _M0L12digit__startS408 = 1;
      } else {
        _M0L12digit__startS408 = 0;
      }
      #line 255 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB24int__to__string__generic(_M0L6bufferS407, _M0L3numS395, _M0L12digit__startS408, _M0L10total__lenS406, _M0L5radixS392);
      _M0L6bufferS396 = _M0L6bufferS407;
      break;
    }
  }
  if (_M0L12is__negativeS394) {
    _M0L6bufferS396[0] = 45;
  }
  return _M0L6bufferS396;
}

int32_t _M0FPB14radix__count32(
  uint32_t _M0L5valueS386,
  int32_t _M0L5radixS388
) {
  uint32_t _M0L4baseS387;
  uint32_t _M0L3numS389;
  int32_t _M0L5countS390;
  #line 189 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS386 == 0u) {
    return 1;
  }
  _M0L4baseS387 = *(uint32_t*)&_M0L5radixS388;
  _M0L3numS389 = _M0L5valueS386;
  _M0L5countS390 = 0;
  while (1) {
    if (_M0L3numS389 > 0u) {
      uint32_t _M0L6_2atmpS1947 = _M0L3numS389 / _M0L4baseS387;
      int32_t _M0L6_2atmpS1948 = _M0L5countS390 + 1;
      _M0L3numS389 = _M0L6_2atmpS1947;
      _M0L5countS390 = _M0L6_2atmpS1948;
      continue;
    } else {
      return _M0L5countS390;
    }
    break;
  }
}

int32_t _M0FPB12hex__count32(uint32_t _M0L5valueS384) {
  #line 177 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS384 == 0u) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS385;
    int32_t _M0L6_2atmpS1946;
    int32_t _M0L6_2atmpS1945;
    #line 182 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS385 = moonbit_clz32(_M0L5valueS384);
    _M0L6_2atmpS1946 = 31 - _M0L14leading__zerosS385;
    _M0L6_2atmpS1945 = _M0L6_2atmpS1946 / 4;
    return _M0L6_2atmpS1945 + 1;
  }
}

int32_t _M0FPB12dec__count32(uint32_t _M0L5valueS383) {
  #line 143 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS383 >= 100000u) {
    if (_M0L5valueS383 >= 10000000u) {
      if (_M0L5valueS383 >= 1000000000u) {
        return 10;
      } else if (_M0L5valueS383 >= 100000000u) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS383 >= 1000000u) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS383 >= 1000u) {
    if (_M0L5valueS383 >= 10000u) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS383 >= 100u) {
    return 3;
  } else if (_M0L5valueS383 >= 10u) {
    return 2;
  } else {
    return 1;
  }
}

int32_t _M0FPB20int__to__string__dec(
  uint16_t* _M0L6bufferS369,
  uint32_t _M0L3numS381,
  int32_t _M0L12digit__startS370,
  int32_t _M0L10total__lenS382
) {
  int32_t _M0L6_2atmpS1944;
  uint32_t _M0L3numS359;
  int32_t _M0L6offsetS360;
  #line 88 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1944 = _M0L10total__lenS382 - _M0L12digit__startS370;
  _M0L3numS359 = _M0L3numS381;
  _M0L6offsetS360 = _M0L6_2atmpS1944;
  while (1) {
    if (_M0L3numS359 >= 10000u) {
      uint32_t _M0L1tS361 = _M0L3numS359 / 10000u;
      uint32_t _M0L6_2atmpS1921 = _M0L3numS359 % 10000u;
      int32_t _M0L1rS362 = *(int32_t*)&_M0L6_2atmpS1921;
      int32_t _M0L2d1S363 = _M0L1rS362 / 100;
      int32_t _M0L2d2S364 = _M0L1rS362 % 100;
      int32_t _M0L6_2atmpS1920 = _M0L2d1S363 / 10;
      int32_t _M0L6_2atmpS1919 = 48 + _M0L6_2atmpS1920;
      int32_t _M0L6d1__hiS365 = (uint16_t)_M0L6_2atmpS1919;
      int32_t _M0L6_2atmpS1918 = _M0L2d1S363 % 10;
      int32_t _M0L6_2atmpS1917 = 48 + _M0L6_2atmpS1918;
      int32_t _M0L6d1__loS366 = (uint16_t)_M0L6_2atmpS1917;
      int32_t _M0L6_2atmpS1916 = _M0L2d2S364 / 10;
      int32_t _M0L6_2atmpS1915 = 48 + _M0L6_2atmpS1916;
      int32_t _M0L6d2__hiS367 = (uint16_t)_M0L6_2atmpS1915;
      int32_t _M0L6_2atmpS1914 = _M0L2d2S364 % 10;
      int32_t _M0L6_2atmpS1913 = 48 + _M0L6_2atmpS1914;
      int32_t _M0L6d2__loS368 = (uint16_t)_M0L6_2atmpS1913;
      int32_t _M0L6_2atmpS1905 = _M0L12digit__startS370 + _M0L6offsetS360;
      int32_t _M0L6_2atmpS1904 = _M0L6_2atmpS1905 - 4;
      int32_t _M0L6_2atmpS1907;
      int32_t _M0L6_2atmpS1906;
      int32_t _M0L6_2atmpS1909;
      int32_t _M0L6_2atmpS1908;
      int32_t _M0L6_2atmpS1911;
      int32_t _M0L6_2atmpS1910;
      int32_t _M0L6_2atmpS1912;
      _M0L6bufferS369[_M0L6_2atmpS1904] = _M0L6d1__hiS365;
      _M0L6_2atmpS1907 = _M0L12digit__startS370 + _M0L6offsetS360;
      _M0L6_2atmpS1906 = _M0L6_2atmpS1907 - 3;
      _M0L6bufferS369[_M0L6_2atmpS1906] = _M0L6d1__loS366;
      _M0L6_2atmpS1909 = _M0L12digit__startS370 + _M0L6offsetS360;
      _M0L6_2atmpS1908 = _M0L6_2atmpS1909 - 2;
      _M0L6bufferS369[_M0L6_2atmpS1908] = _M0L6d2__hiS367;
      _M0L6_2atmpS1911 = _M0L12digit__startS370 + _M0L6offsetS360;
      _M0L6_2atmpS1910 = _M0L6_2atmpS1911 - 1;
      _M0L6bufferS369[_M0L6_2atmpS1910] = _M0L6d2__loS368;
      _M0L6_2atmpS1912 = _M0L6offsetS360 - 4;
      _M0L3numS359 = _M0L1tS361;
      _M0L6offsetS360 = _M0L6_2atmpS1912;
      continue;
    } else {
      int32_t _M0L6_2atmpS1943 = *(int32_t*)&_M0L3numS359;
      int32_t _M0L9remainingS372 = _M0L6_2atmpS1943;
      int32_t _M0L6offsetS373 = _M0L6offsetS360;
      while (1) {
        if (_M0L9remainingS372 >= 100) {
          int32_t _M0L1tS374 = _M0L9remainingS372 / 100;
          int32_t _M0L1dS375 = _M0L9remainingS372 % 100;
          int32_t _M0L6_2atmpS1930 = _M0L1dS375 / 10;
          int32_t _M0L6_2atmpS1929 = 48 + _M0L6_2atmpS1930;
          int32_t _M0L5d__hiS376 = (uint16_t)_M0L6_2atmpS1929;
          int32_t _M0L6_2atmpS1928 = _M0L1dS375 % 10;
          int32_t _M0L6_2atmpS1927 = 48 + _M0L6_2atmpS1928;
          int32_t _M0L5d__loS377 = (uint16_t)_M0L6_2atmpS1927;
          int32_t _M0L6_2atmpS1923 = _M0L12digit__startS370 + _M0L6offsetS373;
          int32_t _M0L6_2atmpS1922 = _M0L6_2atmpS1923 - 2;
          int32_t _M0L6_2atmpS1925;
          int32_t _M0L6_2atmpS1924;
          int32_t _M0L6_2atmpS1926;
          _M0L6bufferS369[_M0L6_2atmpS1922] = _M0L5d__hiS376;
          _M0L6_2atmpS1925 = _M0L12digit__startS370 + _M0L6offsetS373;
          _M0L6_2atmpS1924 = _M0L6_2atmpS1925 - 1;
          _M0L6bufferS369[_M0L6_2atmpS1924] = _M0L5d__loS377;
          _M0L6_2atmpS1926 = _M0L6offsetS373 - 2;
          _M0L9remainingS372 = _M0L1tS374;
          _M0L6offsetS373 = _M0L6_2atmpS1926;
          continue;
        } else if (_M0L9remainingS372 >= 10) {
          int32_t _M0L6_2atmpS1938 = _M0L9remainingS372 / 10;
          int32_t _M0L6_2atmpS1937 = 48 + _M0L6_2atmpS1938;
          int32_t _M0L5d__hiS379 = (uint16_t)_M0L6_2atmpS1937;
          int32_t _M0L6_2atmpS1936 = _M0L9remainingS372 % 10;
          int32_t _M0L6_2atmpS1935 = 48 + _M0L6_2atmpS1936;
          int32_t _M0L5d__loS380 = (uint16_t)_M0L6_2atmpS1935;
          int32_t _M0L6_2atmpS1932 = _M0L12digit__startS370 + _M0L6offsetS373;
          int32_t _M0L6_2atmpS1931 = _M0L6_2atmpS1932 - 2;
          int32_t _M0L6_2atmpS1934;
          int32_t _M0L6_2atmpS1933;
          _M0L6bufferS369[_M0L6_2atmpS1931] = _M0L5d__hiS379;
          _M0L6_2atmpS1934 = _M0L12digit__startS370 + _M0L6offsetS373;
          _M0L6_2atmpS1933 = _M0L6_2atmpS1934 - 1;
          _M0L6bufferS369[_M0L6_2atmpS1933] = _M0L5d__loS380;
        } else {
          int32_t _M0L6_2atmpS1942 = _M0L12digit__startS370 + _M0L6offsetS373;
          int32_t _M0L6_2atmpS1939 = _M0L6_2atmpS1942 - 1;
          int32_t _M0L6_2atmpS1941 = 48 + _M0L9remainingS372;
          int32_t _M0L6_2atmpS1940 = (uint16_t)_M0L6_2atmpS1941;
          _M0L6bufferS369[_M0L6_2atmpS1939] = _M0L6_2atmpS1940;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB24int__to__string__generic(
  uint16_t* _M0L6bufferS349,
  uint32_t _M0L3numS353,
  int32_t _M0L12digit__startS350,
  int32_t _M0L10total__lenS352,
  int32_t _M0L5radixS343
) {
  uint32_t _M0L4baseS342;
  int32_t _M0L6_2atmpS1889;
  int32_t _M0L6_2atmpS1888;
  #line 57 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS342 = *(uint32_t*)&_M0L5radixS343;
  _M0L6_2atmpS1889 = _M0L5radixS343 - 1;
  _M0L6_2atmpS1888 = _M0L5radixS343 & _M0L6_2atmpS1889;
  if (_M0L6_2atmpS1888 == 0) {
    int32_t _M0L5shiftS344;
    uint32_t _M0L4maskS345;
    int32_t _M0L6_2atmpS1896;
    int32_t _M0L6offsetS346;
    uint32_t _M0L1nS347;
    #line 68 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS344 = moonbit_ctz32(_M0L5radixS343);
    _M0L4maskS345 = _M0L4baseS342 - 1u;
    _M0L6_2atmpS1896 = _M0L10total__lenS352 - _M0L12digit__startS350;
    _M0L6offsetS346 = _M0L6_2atmpS1896;
    _M0L1nS347 = _M0L3numS353;
    while (1) {
      if (_M0L1nS347 > 0u) {
        uint32_t _M0L6_2atmpS1895 = _M0L1nS347 & _M0L4maskS345;
        int32_t _M0L5digitS348 = *(int32_t*)&_M0L6_2atmpS1895;
        int32_t _M0L6_2atmpS1892 = _M0L12digit__startS350 + _M0L6offsetS346;
        int32_t _M0L6_2atmpS1890 = _M0L6_2atmpS1892 - 1;
        int32_t _M0L6_2atmpS1891 =
          ((moonbit_string_t)moonbit_string_literal_176.data)[_M0L5digitS348];
        int32_t _M0L6_2atmpS1893;
        uint32_t _M0L6_2atmpS1894;
        _M0L6bufferS349[_M0L6_2atmpS1890] = _M0L6_2atmpS1891;
        _M0L6_2atmpS1893 = _M0L6offsetS346 - 1;
        _M0L6_2atmpS1894 = _M0L1nS347 >> (_M0L5shiftS344 & 31);
        _M0L6offsetS346 = _M0L6_2atmpS1893;
        _M0L1nS347 = _M0L6_2atmpS1894;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1903 = _M0L10total__lenS352 - _M0L12digit__startS350;
    int32_t _M0L6offsetS354 = _M0L6_2atmpS1903;
    uint32_t _M0L1nS355 = _M0L3numS353;
    while (1) {
      if (_M0L1nS355 > 0u) {
        uint32_t _M0L1qS356 = _M0L1nS355 / _M0L4baseS342;
        uint32_t _M0L6_2atmpS1902 = _M0L1qS356 * _M0L4baseS342;
        uint32_t _M0L6_2atmpS1901 = _M0L1nS355 - _M0L6_2atmpS1902;
        int32_t _M0L5digitS357 = *(int32_t*)&_M0L6_2atmpS1901;
        int32_t _M0L6_2atmpS1899 = _M0L12digit__startS350 + _M0L6offsetS354;
        int32_t _M0L6_2atmpS1897 = _M0L6_2atmpS1899 - 1;
        int32_t _M0L6_2atmpS1898 =
          ((moonbit_string_t)moonbit_string_literal_176.data)[_M0L5digitS357];
        int32_t _M0L6_2atmpS1900;
        _M0L6bufferS349[_M0L6_2atmpS1897] = _M0L6_2atmpS1898;
        _M0L6_2atmpS1900 = _M0L6offsetS354 - 1;
        _M0L6offsetS354 = _M0L6_2atmpS1900;
        _M0L1nS355 = _M0L1qS356;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB20int__to__string__hex(
  uint16_t* _M0L6bufferS336,
  uint32_t _M0L3numS341,
  int32_t _M0L12digit__startS337,
  int32_t _M0L10total__lenS340
) {
  int32_t _M0L6_2atmpS1887;
  int32_t _M0L6offsetS331;
  uint32_t _M0L1nS332;
  #line 29 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1887 = _M0L10total__lenS340 - _M0L12digit__startS337;
  _M0L6offsetS331 = _M0L6_2atmpS1887;
  _M0L1nS332 = _M0L3numS341;
  while (1) {
    if (_M0L6offsetS331 >= 2) {
      uint32_t _M0L6_2atmpS1884 = _M0L1nS332 & 255u;
      int32_t _M0L9byte__valS333 = *(int32_t*)&_M0L6_2atmpS1884;
      int32_t _M0L2hiS334 = _M0L9byte__valS333 / 16;
      int32_t _M0L2loS335 = _M0L9byte__valS333 % 16;
      int32_t _M0L6_2atmpS1878 = _M0L12digit__startS337 + _M0L6offsetS331;
      int32_t _M0L6_2atmpS1876 = _M0L6_2atmpS1878 - 2;
      int32_t _M0L6_2atmpS1877 =
        ((moonbit_string_t)moonbit_string_literal_176.data)[_M0L2hiS334];
      int32_t _M0L6_2atmpS1881;
      int32_t _M0L6_2atmpS1879;
      int32_t _M0L6_2atmpS1880;
      int32_t _M0L6_2atmpS1882;
      uint32_t _M0L6_2atmpS1883;
      _M0L6bufferS336[_M0L6_2atmpS1876] = _M0L6_2atmpS1877;
      _M0L6_2atmpS1881 = _M0L12digit__startS337 + _M0L6offsetS331;
      _M0L6_2atmpS1879 = _M0L6_2atmpS1881 - 1;
      _M0L6_2atmpS1880
      = ((moonbit_string_t)moonbit_string_literal_176.data)[
        _M0L2loS335
      ];
      _M0L6bufferS336[_M0L6_2atmpS1879] = _M0L6_2atmpS1880;
      _M0L6_2atmpS1882 = _M0L6offsetS331 - 2;
      _M0L6_2atmpS1883 = _M0L1nS332 >> 8;
      _M0L6offsetS331 = _M0L6_2atmpS1882;
      _M0L1nS332 = _M0L6_2atmpS1883;
      continue;
    } else if (_M0L6offsetS331 == 1) {
      uint32_t _M0L6_2atmpS1886 = _M0L1nS332 & 15u;
      int32_t _M0L6nibbleS339 = *(int32_t*)&_M0L6_2atmpS1886;
      int32_t _M0L6_2atmpS1885 =
        ((moonbit_string_t)moonbit_string_literal_176.data)[_M0L6nibbleS339];
      _M0L6bufferS336[_M0L12digit__startS337] = _M0L6_2atmpS1885;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPB4Iter4nextGsE(struct _M0TPB4IterGsE* _M0L4selfS308) {
  struct _M0TWEOs* _M0L7_2afuncS307;
  moonbit_string_t _M0L6resultS309;
  int64_t _M0L7_2abindS310;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS307 = _M0L4selfS308->$0;
  moonbit_incref(_M0L7_2afuncS307);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS309 = _M0L7_2afuncS307->code(_M0L7_2afuncS307);
  moonbit_decref(_M0L7_2afuncS307);
  _M0L7_2abindS310 = _M0L4selfS308->$1;
  if (_M0L6resultS309 == 0) {
    _M0L4selfS308->$1 = _M0MPB4Iter4nextN6constrS9981GsE;
  } else if (_M0L7_2abindS310 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS311 = _M0L7_2abindS310;
    int32_t _M0L4_2anS312 = (int32_t)_M0L7_2aSomeS311;
    int64_t _M0L6_2atmpS1868;
    if (_M0L4_2anS312 > 0) {
      int32_t _M0L6_2atmpS1869 = _M0L4_2anS312 - 1;
      _M0L6_2atmpS1868 = (int64_t)_M0L6_2atmpS1869;
    } else {
      _M0L6_2atmpS1868 = _M0MPB4Iter4nextN6constrS9980GsE;
    }
    _M0L4selfS308->$1 = _M0L6_2atmpS1868;
  }
  return _M0L6resultS309;
}

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS314
) {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L7_2afuncS313;
  void* _M0L6resultS315;
  int64_t _M0L7_2abindS316;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS313 = _M0L4selfS314->$0;
  moonbit_incref(_M0L7_2afuncS313);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS315 = _M0L7_2afuncS313->code(_M0L7_2afuncS313);
  moonbit_decref(_M0L7_2afuncS313);
  _M0L7_2abindS316 = _M0L4selfS314->$1;
  switch (Moonbit_object_tag(_M0L6resultS315)) {
    case 1: {
      if (_M0L7_2abindS316 == 4294967296ll) {
        
      } else {
        int64_t _M0L7_2aSomeS317 = _M0L7_2abindS316;
        int32_t _M0L4_2anS318 = (int32_t)_M0L7_2aSomeS317;
        int64_t _M0L6_2atmpS1870;
        if (_M0L4_2anS318 > 0) {
          int32_t _M0L6_2atmpS1871 = _M0L4_2anS318 - 1;
          _M0L6_2atmpS1870 = (int64_t)_M0L6_2atmpS1871;
        } else {
          _M0L6_2atmpS1870
          = _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE;
        }
        _M0L4selfS314->$1 = _M0L6_2atmpS1870;
      }
      break;
    }
    default: {
      _M0L4selfS314->$1
      = _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE;
      break;
    }
  }
  return _M0L6resultS315;
}

struct _M0TUsRPB5ArrayGiEE* _M0MPB4Iter4nextGUsRPB5ArrayGiEEE(
  struct _M0TPB4IterGUsRPB5ArrayGiEEE* _M0L4selfS320
) {
  struct _M0TWEOUsRPB5ArrayGiEE* _M0L7_2afuncS319;
  struct _M0TUsRPB5ArrayGiEE* _M0L6resultS321;
  int64_t _M0L7_2abindS322;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS319 = _M0L4selfS320->$0;
  moonbit_incref(_M0L7_2afuncS319);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS321 = _M0L7_2afuncS319->code(_M0L7_2afuncS319);
  moonbit_decref(_M0L7_2afuncS319);
  _M0L7_2abindS322 = _M0L4selfS320->$1;
  if (_M0L6resultS321 == 0) {
    _M0L4selfS320->$1 = _M0MPB4Iter4nextN6constrS9981GUsRPB5ArrayGiEEE;
  } else if (_M0L7_2abindS322 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS323 = _M0L7_2abindS322;
    int32_t _M0L4_2anS324 = (int32_t)_M0L7_2aSomeS323;
    int64_t _M0L6_2atmpS1872;
    if (_M0L4_2anS324 > 0) {
      int32_t _M0L6_2atmpS1873 = _M0L4_2anS324 - 1;
      _M0L6_2atmpS1872 = (int64_t)_M0L6_2atmpS1873;
    } else {
      _M0L6_2atmpS1872 = _M0MPB4Iter4nextN6constrS9980GUsRPB5ArrayGiEEE;
    }
    _M0L4selfS320->$1 = _M0L6_2atmpS1872;
  }
  return _M0L6resultS321;
}

int32_t _M0MPB4Iter4nextGcE(struct _M0TPB4IterGcE* _M0L4selfS326) {
  struct _M0TWEOc* _M0L7_2afuncS325;
  int32_t _M0L6resultS327;
  int64_t _M0L7_2abindS328;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS325 = _M0L4selfS326->$0;
  moonbit_incref(_M0L7_2afuncS325);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS327 = _M0L7_2afuncS325->code(_M0L7_2afuncS325);
  moonbit_decref(_M0L7_2afuncS325);
  _M0L7_2abindS328 = _M0L4selfS326->$1;
  if (_M0L6resultS327 == -1) {
    _M0L4selfS326->$1 = _M0MPB4Iter4nextN6constrS9981GcE;
  } else if (_M0L7_2abindS328 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS329 = _M0L7_2abindS328;
    int32_t _M0L4_2anS330 = (int32_t)_M0L7_2aSomeS329;
    int64_t _M0L6_2atmpS1874;
    if (_M0L4_2anS330 > 0) {
      int32_t _M0L6_2atmpS1875 = _M0L4_2anS330 - 1;
      _M0L6_2atmpS1874 = (int64_t)_M0L6_2atmpS1875;
    } else {
      _M0L6_2atmpS1874 = _M0MPB4Iter4nextN6constrS9980GcE;
    }
    _M0L4selfS326->$1 = _M0L6_2atmpS1874;
  }
  return _M0L6resultS327;
}

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void* _M0L4selfS306
) {
  struct _M0TPB13StringBuilder* _M0L6loggerS305;
  struct _M0TPB6Logger _M0L6_2atmpS1867;
  moonbit_string_t _result_4592;
  #line 165 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 166 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS305 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  moonbit_incref(_M0L6loggerS305);
  _M0L6_2atmpS1867
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L6loggerS305
  };
  #line 167 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB7FailurePB4Show6output(_M0L4selfS306, _M0L6_2atmpS1867);
  if (_M0L6_2atmpS1867.$1) {
    moonbit_decref(_M0L6_2atmpS1867.$1);
  }
  #line 168 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _result_4592 = _M0MPB13StringBuilder10to__string(_M0L6loggerS305);
  moonbit_decref(_M0L6loggerS305);
  return _result_4592;
}

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t _M0L4selfS300,
  struct _M0TPB6Logger _M0L6loggerS299
) {
  moonbit_string_t _M0L6_2atmpS1864;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1864 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS300);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS299.$0->$method_0(_M0L6loggerS299.$1, _M0L6_2atmpS1864);
  moonbit_decref(_M0L6_2atmpS1864);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t _M0L4selfS302,
  struct _M0TPB6Logger _M0L6loggerS301
) {
  moonbit_string_t _M0L6_2atmpS1865;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1865 = _M0IPC13int3IntPB4Show10to__string(_M0L4selfS302);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS301.$0->$method_0(_M0L6loggerS301.$1, _M0L6_2atmpS1865);
  moonbit_decref(_M0L6_2atmpS1865);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGbE(
  int32_t _M0L4selfS304,
  struct _M0TPB6Logger _M0L6loggerS303
) {
  moonbit_string_t _M0L6_2atmpS1866;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1866 = _M0IPC14bool4BoolPB4Show10to__string(_M0L4selfS304);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS303.$0->$method_0(_M0L6loggerS303.$1, _M0L6_2atmpS1866);
  moonbit_decref(_M0L6_2atmpS1866);
  return 0;
}

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView _M0L4selfS298
) {
  #line 99 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  return _M0L4selfS298.$1;
}

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView _M0L4selfS297
) {
  moonbit_string_t _M0L8_2afieldS3885;
  #line 92 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L8_2afieldS3885 = _M0L4selfS297.$0;
  moonbit_incref(_M0L8_2afieldS3885);
  return _M0L8_2afieldS3885;
}

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS293,
  moonbit_string_t _M0L5valueS294,
  int32_t _M0L5startS295,
  int32_t _M0L3lenS296
) {
  int32_t _M0L6_2atmpS1863;
  int64_t _M0L6_2atmpS1862;
  struct _M0TPC16string10StringView _M0L6_2atmpS1861;
  #line 122 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1863 = _M0L5startS295 + _M0L3lenS296;
  _M0L6_2atmpS1862 = (int64_t)_M0L6_2atmpS1863;
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1861
  = _M0MPC16string6String11sub_2einner(_M0L5valueS294, _M0L5startS295, _M0L6_2atmpS1862);
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L4selfS293, _M0L6_2atmpS1861);
  moonbit_decref(_M0L6_2atmpS1861.$0);
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t _M0L4selfS286,
  int32_t _M0L5startS292,
  int64_t _M0L3endS288
) {
  int32_t _M0L3lenS285;
  int32_t _M0L3endS287;
  int32_t _M0L5startS291;
  int32_t _if__result_4593;
  #line 755 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3lenS285 = Moonbit_array_length(_M0L4selfS286);
  if (_M0L3endS288 == 4294967296ll) {
    _M0L3endS287 = _M0L3lenS285;
  } else {
    int64_t _M0L7_2aSomeS289 = _M0L3endS288;
    int32_t _M0L6_2aendS290 = (int32_t)_M0L7_2aSomeS289;
    if (_M0L6_2aendS290 < 0) {
      _M0L3endS287 = _M0L3lenS285 + _M0L6_2aendS290;
    } else {
      _M0L3endS287 = _M0L6_2aendS290;
    }
  }
  if (_M0L5startS292 < 0) {
    _M0L5startS291 = _M0L3lenS285 + _M0L5startS292;
  } else {
    _M0L5startS291 = _M0L5startS292;
  }
  if (_M0L5startS291 >= 0) {
    if (_M0L5startS291 <= _M0L3endS287) {
      _if__result_4593 = _M0L3endS287 <= _M0L3lenS285;
    } else {
      _if__result_4593 = 0;
    }
  } else {
    _if__result_4593 = 0;
  }
  if (_if__result_4593) {
    if (_M0L5startS291 < _M0L3lenS285) {
      int32_t _M0L6_2atmpS1858 = _M0L4selfS286[_M0L5startS291];
      int32_t _M0L6_2atmpS1857;
      #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1857
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1858);
      if (!_M0L6_2atmpS1857) {
        
      } else {
        #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L3endS287 < _M0L3lenS285) {
      int32_t _M0L6_2atmpS1860 = _M0L4selfS286[_M0L3endS287];
      int32_t _M0L6_2atmpS1859;
      #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1859
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1860);
      if (!_M0L6_2atmpS1859) {
        
      } else {
        #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    moonbit_incref(_M0L4selfS286);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS286,
                                                 .$1 = _M0L5startS291,
                                                 .$2 = _M0L3endS287};
  } else {
    #line 763 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS284,
  struct _M0TPB4Show _M0L4showS283
) {
  struct _M0TPB6Logger _M0L6_2atmpS1856;
  #line 116 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS284);
  _M0L6_2atmpS1856
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS284
  };
  #line 117 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS283.$0->$method_0(_M0L4showS283.$1, _M0L6_2atmpS1856);
  if (_M0L6_2atmpS1856.$1) {
    moonbit_decref(_M0L6_2atmpS1856.$1);
  }
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS282,
  struct _M0TPB4Show _M0L4showS281
) {
  struct _M0TPB6Logger _M0L6_2atmpS1855;
  #line 111 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS282);
  _M0L6_2atmpS1855
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS282
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS281.$0->$method_0(_M0L4showS281.$1, _M0L6_2atmpS1855);
  if (_M0L6_2atmpS1855.$1) {
    moonbit_decref(_M0L6_2atmpS1855.$1);
  }
  return 0;
}

int32_t _M0FPB13finalize__acc(uint32_t _M0L3accS280) {
  uint32_t _M0L6_2atmpS1854;
  #line 444 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  #line 445 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1854 = _M0FPB14avalanche__acc(_M0L3accS280);
  return *(int32_t*)&_M0L6_2atmpS1854;
}

uint32_t _M0FPB14avalanche__acc(uint32_t _M0L3accS279) {
  uint32_t _M0Lm3accS278;
  uint32_t _M0L6_2atmpS1843;
  uint32_t _M0L6_2atmpS1845;
  uint32_t _M0L6_2atmpS1844;
  uint32_t _M0L6_2atmpS1846;
  uint32_t _M0L6_2atmpS1847;
  uint32_t _M0L6_2atmpS1849;
  uint32_t _M0L6_2atmpS1848;
  uint32_t _M0L6_2atmpS1850;
  uint32_t _M0L6_2atmpS1851;
  uint32_t _M0L6_2atmpS1853;
  uint32_t _M0L6_2atmpS1852;
  #line 449 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0Lm3accS278 = _M0L3accS279;
  _M0L6_2atmpS1843 = _M0Lm3accS278;
  _M0L6_2atmpS1845 = _M0Lm3accS278;
  _M0L6_2atmpS1844 = _M0L6_2atmpS1845 >> 15;
  _M0Lm3accS278 = _M0L6_2atmpS1843 ^ _M0L6_2atmpS1844;
  _M0L6_2atmpS1846 = _M0Lm3accS278;
  _M0Lm3accS278 = _M0L6_2atmpS1846 * 2246822519u;
  _M0L6_2atmpS1847 = _M0Lm3accS278;
  _M0L6_2atmpS1849 = _M0Lm3accS278;
  _M0L6_2atmpS1848 = _M0L6_2atmpS1849 >> 13;
  _M0Lm3accS278 = _M0L6_2atmpS1847 ^ _M0L6_2atmpS1848;
  _M0L6_2atmpS1850 = _M0Lm3accS278;
  _M0Lm3accS278 = _M0L6_2atmpS1850 * 3266489917u;
  _M0L6_2atmpS1851 = _M0Lm3accS278;
  _M0L6_2atmpS1853 = _M0Lm3accS278;
  _M0L6_2atmpS1852 = _M0L6_2atmpS1853 >> 16;
  _M0Lm3accS278 = _M0L6_2atmpS1851 ^ _M0L6_2atmpS1852;
  return _M0Lm3accS278;
}

int32_t _M0IP016_24default__implPB2Eq10not__equalGsE(
  moonbit_string_t _M0L1xS276,
  moonbit_string_t _M0L1yS277
) {
  int32_t _M0L6_2atmpS1842;
  #line 23 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 24 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1842
  = _M0L1xS276 == _M0L1yS277
    || Moonbit_array_length(_M0L1xS276) == Moonbit_array_length(_M0L1yS277)
       && 0
          == memcmp(_M0L1xS276, _M0L1yS277, Moonbit_array_length(_M0L1xS276) * 2);
  return !_M0L6_2atmpS1842;
}

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder* _M0L4selfS275,
  moonbit_string_t _M0L3strS274
) {
  int32_t _M0L8str__lenS273;
  int32_t _M0L3lenS1837;
  int32_t _M0L6_2atmpS1836;
  uint16_t* _M0L4dataS1838;
  int32_t _M0L3lenS1839;
  int32_t _M0L3lenS1841;
  int32_t _M0L6_2atmpS1840;
  #line 86 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L8str__lenS273 = Moonbit_array_length(_M0L3strS274);
  _M0L3lenS1837 = _M0L4selfS275->$1;
  _M0L6_2atmpS1836 = _M0L3lenS1837 + _M0L8str__lenS273;
  #line 88 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS275, _M0L6_2atmpS1836);
  _M0L4dataS1838 = _M0L4selfS275->$0;
  _M0L3lenS1839 = _M0L4selfS275->$1;
  moonbit_incref(_M0L4dataS1838);
  #line 89 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1838, _M0L3lenS1839, _M0L3strS274, 0, _M0L8str__lenS273);
  moonbit_decref(_M0L4dataS1838);
  _M0L3lenS1841 = _M0L4selfS275->$1;
  _M0L6_2atmpS1840 = _M0L3lenS1841 + _M0L8str__lenS273;
  _M0L4selfS275->$1 = _M0L6_2atmpS1840;
  return 0;
}

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t* _M0L4selfS269,
  int32_t _M0L11dst__offsetS272,
  moonbit_string_t _M0L3strS270,
  int32_t _M0L11str__offsetS265,
  int32_t _M0L3lenS266
) {
  int32_t _M0L16end__str__offsetS264;
  int32_t _M0L1iS267;
  int32_t _M0L1jS268;
  #line 71 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L16end__str__offsetS264 = _M0L11str__offsetS265 + _M0L3lenS266;
  _M0L1iS267 = _M0L11str__offsetS265;
  _M0L1jS268 = _M0L11dst__offsetS272;
  while (1) {
    if (_M0L1iS267 < _M0L16end__str__offsetS264) {
      int32_t _M0L6_2atmpS1833 = _M0L3strS270[_M0L1iS267];
      int32_t _M0L6_2atmpS1834;
      int32_t _M0L6_2atmpS1835;
      if (
        _M0L1jS268 < 0 || _M0L1jS268 >= Moonbit_array_length(_M0L4selfS269)
      ) {
        #line 80 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
        moonbit_panic();
      }
      _M0L4selfS269[_M0L1jS268] = _M0L6_2atmpS1833;
      _M0L6_2atmpS1834 = _M0L1iS267 + 1;
      _M0L6_2atmpS1835 = _M0L1jS268 + 1;
      _M0L1iS267 = _M0L6_2atmpS1834;
      _M0L1jS268 = _M0L6_2atmpS1835;
      continue;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t _M0L4selfS262,
  int32_t _M0L5quoteS263
) {
  struct _M0TPB13StringBuilder* _M0L3bufS261;
  int32_t _M0L6_2atmpS1832;
  struct _M0TPC16string10StringView _M0L6_2atmpS1830;
  struct _M0TPB6Logger _M0L6_2atmpS1831;
  moonbit_string_t _result_4595;
  #line 110 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 111 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L3bufS261 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS1832 = Moonbit_array_length(_M0L4selfS262);
  moonbit_incref(_M0L4selfS262);
  _M0L6_2atmpS1830
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS262, .$1 = 0, .$2 = _M0L6_2atmpS1832
  };
  moonbit_incref(_M0L3bufS261);
  _M0L6_2atmpS1831
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS261
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0MPC16string10StringView18escape__to_2einner(_M0L6_2atmpS1830, _M0L6_2atmpS1831, _M0L5quoteS263);
  moonbit_decref(_M0L6_2atmpS1830.$0);
  if (_M0L6_2atmpS1831.$1) {
    moonbit_decref(_M0L6_2atmpS1831.$1);
  }
  #line 113 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_4595 = _M0MPB13StringBuilder10to__string(_M0L3bufS261);
  moonbit_decref(_M0L3bufS261);
  return _result_4595;
}

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView _M0L4selfS253,
  struct _M0TPB6Logger _M0L6loggerS251,
  int32_t _M0L5quoteS250
) {
  int32_t _M0L3endS1828;
  int32_t _M0L5startS1829;
  int32_t _M0L3lenS252;
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS254;
  int32_t _M0L1iS255;
  int32_t _M0L3segS256;
  #line 144 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L5quoteS250) {
    #line 150 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS251.$0->$method_3(_M0L6loggerS251.$1, 34);
  }
  _M0L3endS1828 = _M0L4selfS253.$2;
  _M0L5startS1829 = _M0L4selfS253.$1;
  _M0L3lenS252 = _M0L3endS1828 - _M0L5startS1829;
  moonbit_incref(_M0L4selfS253.$0);
  if (_M0L6loggerS251.$1) {
    moonbit_incref(_M0L6loggerS251.$1);
  }
  _M0L6_2aenvS254
  = (struct _M0TURPC16string10StringViewRPB6LoggerE*)moonbit_malloc(sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE));
  Moonbit_object_header(_M0L6_2aenvS254)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 154, 0);
  _M0L6_2aenvS254->$0 = _M0L4selfS253;
  _M0L6_2aenvS254->$1 = _M0L6loggerS251;
  _M0L1iS255 = 0;
  _M0L3segS256 = 0;
  _2afor_257:;
  while (1) {
    moonbit_string_t _M0L3strS1825;
    int32_t _M0L5startS1827;
    int32_t _M0L6_2atmpS1826;
    int32_t _M0L4codeS258;
    int32_t _M0L1cS260;
    int32_t _M0L6_2atmpS1809;
    int32_t _M0L6_2atmpS1810;
    int32_t _M0L6_2atmpS1811;
    int32_t _tmp_4599;
    int32_t _tmp_4600;
    if (_M0L1iS255 >= _M0L3lenS252) {
      #line 160 "/home/developer/.moon/lib/core/builtin/show.mbt"
      _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
      moonbit_decref(_M0L6_2aenvS254);
      break;
    }
    _M0L3strS1825 = _M0L4selfS253.$0;
    _M0L5startS1827 = _M0L4selfS253.$1;
    _M0L6_2atmpS1826 = _M0L5startS1827 + _M0L1iS255;
    _M0L4codeS258 = _M0L3strS1825[_M0L6_2atmpS1826];
    switch (_M0L4codeS258) {
      case 34: {
        _M0L1cS260 = _M0L4codeS258;
        goto join_259;
        break;
      }
      
      case 92: {
        _M0L1cS260 = _M0L4codeS258;
        goto join_259;
        break;
      }
      
      case 10: {
        int32_t _M0L6_2atmpS1812;
        int32_t _M0L6_2atmpS1813;
        #line 172 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
        #line 173 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS251.$0->$method_0(_M0L6loggerS251.$1, (moonbit_string_t)moonbit_string_literal_177.data);
        _M0L6_2atmpS1812 = _M0L1iS255 + 1;
        _M0L6_2atmpS1813 = _M0L1iS255 + 1;
        _M0L1iS255 = _M0L6_2atmpS1812;
        _M0L3segS256 = _M0L6_2atmpS1813;
        goto _2afor_257;
        break;
      }
      
      case 13: {
        int32_t _M0L6_2atmpS1814;
        int32_t _M0L6_2atmpS1815;
        #line 177 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
        #line 178 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS251.$0->$method_0(_M0L6loggerS251.$1, (moonbit_string_t)moonbit_string_literal_178.data);
        _M0L6_2atmpS1814 = _M0L1iS255 + 1;
        _M0L6_2atmpS1815 = _M0L1iS255 + 1;
        _M0L1iS255 = _M0L6_2atmpS1814;
        _M0L3segS256 = _M0L6_2atmpS1815;
        goto _2afor_257;
        break;
      }
      
      case 8: {
        int32_t _M0L6_2atmpS1816;
        int32_t _M0L6_2atmpS1817;
        #line 182 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
        #line 183 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS251.$0->$method_0(_M0L6loggerS251.$1, (moonbit_string_t)moonbit_string_literal_179.data);
        _M0L6_2atmpS1816 = _M0L1iS255 + 1;
        _M0L6_2atmpS1817 = _M0L1iS255 + 1;
        _M0L1iS255 = _M0L6_2atmpS1816;
        _M0L3segS256 = _M0L6_2atmpS1817;
        goto _2afor_257;
        break;
      }
      
      case 9: {
        int32_t _M0L6_2atmpS1818;
        int32_t _M0L6_2atmpS1819;
        #line 187 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
        #line 188 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS251.$0->$method_0(_M0L6loggerS251.$1, (moonbit_string_t)moonbit_string_literal_180.data);
        _M0L6_2atmpS1818 = _M0L1iS255 + 1;
        _M0L6_2atmpS1819 = _M0L1iS255 + 1;
        _M0L1iS255 = _M0L6_2atmpS1818;
        _M0L3segS256 = _M0L6_2atmpS1819;
        goto _2afor_257;
        break;
      }
      default: {
        if (_M0L4codeS258 < 32) {
          int32_t _M0L6_2atmpS1821;
          moonbit_string_t _M0L6_2atmpS1820;
          int32_t _M0L6_2atmpS1822;
          int32_t _M0L6_2atmpS1823;
          #line 193 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
          #line 194 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS251.$0->$method_0(_M0L6loggerS251.$1, (moonbit_string_t)moonbit_string_literal_181.data);
          _M0L6_2atmpS1821 = _M0L4codeS258 & 0xff;
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6_2atmpS1820 = _M0MPC14byte4Byte7to__hex(_M0L6_2atmpS1821);
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS251.$0->$method_0(_M0L6loggerS251.$1, _M0L6_2atmpS1820);
          moonbit_decref(_M0L6_2atmpS1820);
          #line 196 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS251.$0->$method_3(_M0L6loggerS251.$1, 125);
          _M0L6_2atmpS1822 = _M0L1iS255 + 1;
          _M0L6_2atmpS1823 = _M0L1iS255 + 1;
          _M0L1iS255 = _M0L6_2atmpS1822;
          _M0L3segS256 = _M0L6_2atmpS1823;
          goto _2afor_257;
        } else {
          int32_t _M0L6_2atmpS1824 = _M0L1iS255 + 1;
          int32_t _tmp_4598 = _M0L3segS256;
          _M0L1iS255 = _M0L6_2atmpS1824;
          _M0L3segS256 = _tmp_4598;
          goto _2afor_257;
        }
        break;
      }
    }
    goto joinlet_4597;
    join_259:;
    #line 166 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS254, _M0L3segS256, _M0L1iS255);
    #line 167 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS251.$0->$method_3(_M0L6loggerS251.$1, 92);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1809 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS260);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS251.$0->$method_3(_M0L6loggerS251.$1, _M0L6_2atmpS1809);
    _M0L6_2atmpS1810 = _M0L1iS255 + 1;
    _M0L6_2atmpS1811 = _M0L1iS255 + 1;
    _M0L1iS255 = _M0L6_2atmpS1810;
    _M0L3segS256 = _M0L6_2atmpS1811;
    continue;
    joinlet_4597:;
    _tmp_4599 = _M0L1iS255;
    _tmp_4600 = _M0L3segS256;
    _M0L1iS255 = _tmp_4599;
    _M0L3segS256 = _tmp_4600;
    continue;
    break;
  }
  if (_M0L5quoteS250) {
    #line 204 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS251.$0->$method_3(_M0L6loggerS251.$1, 34);
  }
  return 0;
}

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS246,
  int32_t _M0L3segS249,
  int32_t _M0L1iS248
) {
  struct _M0TPB6Logger _M0L6loggerS245;
  struct _M0TPC16string10StringView _M0L4selfS247;
  #line 153 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6loggerS245 = _M0L6_2aenvS246->$1;
  _M0L4selfS247 = _M0L6_2aenvS246->$0;
  if (_M0L1iS248 > _M0L3segS249) {
    int64_t _M0L6_2atmpS1808 = (int64_t)_M0L1iS248;
    struct _M0TPC16string10StringView _M0L6_2atmpS1807;
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1807
    = _M0MPC16string10StringView11sub_2einner(_M0L4selfS247, _M0L3segS249, _M0L6_2atmpS1808);
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS245.$0->$method_2(_M0L6loggerS245.$1, _M0L6_2atmpS1807);
    moonbit_decref(_M0L6_2atmpS1807.$0);
  }
  return 0;
}

int32_t _M0MPC16string10StringView11unsafe__get(
  struct _M0TPC16string10StringView _M0L4selfS243,
  int32_t _M0L5indexS244
) {
  moonbit_string_t _M0L3strS1804;
  int32_t _M0L5startS1806;
  int32_t _M0L6_2atmpS1805;
  #line 128 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1804 = _M0L4selfS243.$0;
  _M0L5startS1806 = _M0L4selfS243.$1;
  _M0L6_2atmpS1805 = _M0L5startS1806 + _M0L5indexS244;
  return _M0L3strS1804[_M0L6_2atmpS1805];
}

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView _M0L4selfS236,
  int32_t _M0L5startS242,
  int64_t _M0L3endS238
) {
  moonbit_string_t _M0L3strS1803;
  int32_t _M0L8str__lenS235;
  int32_t _M0L8abs__endS237;
  int32_t _M0L10abs__startS241;
  int32_t _M0L5startS1791;
  int32_t _if__result_4601;
  #line 814 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1803 = _M0L4selfS236.$0;
  _M0L8str__lenS235 = Moonbit_array_length(_M0L3strS1803);
  if (_M0L3endS238 == 4294967296ll) {
    _M0L8abs__endS237 = _M0L4selfS236.$2;
  } else {
    int64_t _M0L7_2aSomeS239 = _M0L3endS238;
    int32_t _M0L6_2aendS240 = (int32_t)_M0L7_2aSomeS239;
    if (_M0L6_2aendS240 < 0) {
      int32_t _M0L3endS1801 = _M0L4selfS236.$2;
      _M0L8abs__endS237 = _M0L3endS1801 + _M0L6_2aendS240;
    } else {
      int32_t _M0L5startS1802 = _M0L4selfS236.$1;
      _M0L8abs__endS237 = _M0L5startS1802 + _M0L6_2aendS240;
    }
  }
  if (_M0L5startS242 < 0) {
    int32_t _M0L3endS1799 = _M0L4selfS236.$2;
    _M0L10abs__startS241 = _M0L3endS1799 + _M0L5startS242;
  } else {
    int32_t _M0L5startS1800 = _M0L4selfS236.$1;
    _M0L10abs__startS241 = _M0L5startS1800 + _M0L5startS242;
  }
  _M0L5startS1791 = _M0L4selfS236.$1;
  if (_M0L10abs__startS241 >= _M0L5startS1791) {
    if (_M0L10abs__startS241 <= _M0L8abs__endS237) {
      int32_t _M0L3endS1790 = _M0L4selfS236.$2;
      _if__result_4601 = _M0L8abs__endS237 <= _M0L3endS1790;
    } else {
      _if__result_4601 = 0;
    }
  } else {
    _if__result_4601 = 0;
  }
  if (_if__result_4601) {
    moonbit_string_t _M0L3strS1798;
    if (_M0L10abs__startS241 < _M0L8str__lenS235) {
      moonbit_string_t _M0L3strS1794 = _M0L4selfS236.$0;
      int32_t _M0L6_2atmpS1793 = _M0L3strS1794[_M0L10abs__startS241];
      int32_t _M0L6_2atmpS1792;
      #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1792
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1793);
      if (!_M0L6_2atmpS1792) {
        
      } else {
        #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L8abs__endS237 < _M0L8str__lenS235) {
      moonbit_string_t _M0L3strS1797 = _M0L4selfS236.$0;
      int32_t _M0L6_2atmpS1796 = _M0L3strS1797[_M0L8abs__endS237];
      int32_t _M0L6_2atmpS1795;
      #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1795
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1796);
      if (!_M0L6_2atmpS1795) {
        
      } else {
        #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    _M0L3strS1798 = _M0L4selfS236.$0;
    moonbit_incref(_M0L3strS1798);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1798,
                                                 .$1 = _M0L10abs__startS241,
                                                 .$2 = _M0L8abs__endS237};
  } else {
    #line 834 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC16string10StringView6length(
  struct _M0TPC16string10StringView _M0L4selfS234
) {
  int32_t _M0L3endS1788;
  int32_t _M0L5startS1789;
  #line 49 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3endS1788 = _M0L4selfS234.$2;
  _M0L5startS1789 = _M0L4selfS234.$1;
  return _M0L3endS1788 - _M0L5startS1789;
}

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t _M0L1bS233) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS232;
  int32_t _M0L6_2atmpS1785;
  int32_t _M0L6_2atmpS1784;
  int32_t _M0L6_2atmpS1787;
  int32_t _M0L6_2atmpS1786;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1783;
  moonbit_string_t _result_4602;
  #line 74 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L7_2aselfS232 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1785 = _M0IPC14byte4BytePB3Div3div(_M0L1bS233, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1784
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1785);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS232, _M0L6_2atmpS1784);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1787 = _M0IPC14byte4BytePB3Mod3mod(_M0L1bS233, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1786
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1787);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS232, _M0L6_2atmpS1786);
  _M0L6_2atmpS1783 = _M0L7_2aselfS232;
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_4602 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1783);
  moonbit_decref(_M0L6_2atmpS1783);
  return _result_4602;
}

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t _M0L1iS231) {
  #line 75 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L1iS231 < 10) {
    int32_t _M0L6_2atmpS1780;
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1780 = _M0IPC14byte4BytePB3Add3add(_M0L1iS231, 48);
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1780);
  } else {
    int32_t _M0L6_2atmpS1782;
    int32_t _M0L6_2atmpS1781;
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1782 = _M0IPC14byte4BytePB3Add3add(_M0L1iS231, 97);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1781 = _M0IPC14byte4BytePB3Sub3sub(_M0L6_2atmpS1782, 10);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1781);
  }
}

int32_t _M0IPC14byte4BytePB3Sub3sub(
  int32_t _M0L4selfS229,
  int32_t _M0L4thatS230
) {
  int32_t _M0L6_2atmpS1778;
  int32_t _M0L6_2atmpS1779;
  int32_t _M0L6_2atmpS1777;
  #line 120 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1778 = (int32_t)_M0L4selfS229;
  _M0L6_2atmpS1779 = (int32_t)_M0L4thatS230;
  _M0L6_2atmpS1777 = _M0L6_2atmpS1778 - _M0L6_2atmpS1779;
  return _M0L6_2atmpS1777 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Mod3mod(
  int32_t _M0L4selfS227,
  int32_t _M0L4thatS228
) {
  int32_t _M0L6_2atmpS1775;
  int32_t _M0L6_2atmpS1776;
  int32_t _M0L6_2atmpS1774;
  #line 67 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1775 = (int32_t)_M0L4selfS227;
  _M0L6_2atmpS1776 = (int32_t)_M0L4thatS228;
  _M0L6_2atmpS1774 = _M0L6_2atmpS1775 % _M0L6_2atmpS1776;
  return _M0L6_2atmpS1774 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Div3div(
  int32_t _M0L4selfS225,
  int32_t _M0L4thatS226
) {
  int32_t _M0L6_2atmpS1772;
  int32_t _M0L6_2atmpS1773;
  int32_t _M0L6_2atmpS1771;
  #line 62 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1772 = (int32_t)_M0L4selfS225;
  _M0L6_2atmpS1773 = (int32_t)_M0L4thatS226;
  _M0L6_2atmpS1771 = _M0L6_2atmpS1772 / _M0L6_2atmpS1773;
  return _M0L6_2atmpS1771 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Add3add(
  int32_t _M0L4selfS223,
  int32_t _M0L4thatS224
) {
  int32_t _M0L6_2atmpS1769;
  int32_t _M0L6_2atmpS1770;
  int32_t _M0L6_2atmpS1768;
  #line 106 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1769 = (int32_t)_M0L4selfS223;
  _M0L6_2atmpS1770 = (int32_t)_M0L4thatS224;
  _M0L6_2atmpS1768 = _M0L6_2atmpS1769 + _M0L6_2atmpS1770;
  return _M0L6_2atmpS1768 & 0xff;
}

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t _M0L4selfS222) {
  int32_t _M0L6_2atmpS1767;
  #line 68 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  _M0L6_2atmpS1767 = (int32_t)_M0L4selfS222;
  return _M0L6_2atmpS1767;
}

int32_t _M0FPB32code__point__of__surrogate__pair(
  int32_t _M0L7leadingS220,
  int32_t _M0L8trailingS221
) {
  int32_t _M0L6_2atmpS1766;
  int32_t _M0L6_2atmpS1765;
  int32_t _M0L6_2atmpS1764;
  int32_t _M0L6_2atmpS1763;
  int32_t _M0L6_2atmpS1762;
  #line 40 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1766 = _M0L7leadingS220 - 55296;
  _M0L6_2atmpS1765 = _M0L6_2atmpS1766 * 1024;
  _M0L6_2atmpS1764 = _M0L6_2atmpS1765 + _M0L8trailingS221;
  _M0L6_2atmpS1763 = _M0L6_2atmpS1764 - 56320;
  _M0L6_2atmpS1762 = _M0L6_2atmpS1763 + 65536;
  return _M0L6_2atmpS1762;
}

int32_t _M0MPC16string6String20char__length_2einner(
  moonbit_string_t _M0L4selfS213,
  int32_t _M0L13start__offsetS214,
  int64_t _M0L11end__offsetS211
) {
  int32_t _M0L11end__offsetS210;
  int32_t _if__result_4603;
  #line 60 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L11end__offsetS211 == 4294967296ll) {
    _M0L11end__offsetS210 = Moonbit_array_length(_M0L4selfS213);
  } else {
    int64_t _M0L7_2aSomeS212 = _M0L11end__offsetS211;
    _M0L11end__offsetS210 = (int32_t)_M0L7_2aSomeS212;
  }
  if (_M0L13start__offsetS214 >= 0) {
    if (_M0L13start__offsetS214 <= _M0L11end__offsetS210) {
      int32_t _M0L6_2atmpS1755 = Moonbit_array_length(_M0L4selfS213);
      _if__result_4603 = _M0L11end__offsetS210 <= _M0L6_2atmpS1755;
    } else {
      _if__result_4603 = 0;
    }
  } else {
    _if__result_4603 = 0;
  }
  if (_if__result_4603) {
    int32_t _M0L12utf16__indexS215 = _M0L13start__offsetS214;
    int32_t _M0L11char__countS216 = 0;
    while (1) {
      if (_M0L12utf16__indexS215 < _M0L11end__offsetS210) {
        int32_t _M0L2c1S217 = _M0L4selfS213[_M0L12utf16__indexS215];
        int32_t _if__result_4605;
        int32_t _M0L6_2atmpS1760;
        int32_t _M0L6_2atmpS1761;
        #line 76 "/home/developer/.moon/lib/core/builtin/string.mbt"
        if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S217)) {
          int32_t _M0L6_2atmpS1756 = _M0L12utf16__indexS215 + 1;
          _if__result_4605 = _M0L6_2atmpS1756 < _M0L11end__offsetS210;
        } else {
          _if__result_4605 = 0;
        }
        if (_if__result_4605) {
          int32_t _M0L6_2atmpS1759 = _M0L12utf16__indexS215 + 1;
          int32_t _M0L2c2S218 = _M0L4selfS213[_M0L6_2atmpS1759];
          #line 78 "/home/developer/.moon/lib/core/builtin/string.mbt"
          if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S218)) {
            int32_t _M0L6_2atmpS1757 = _M0L12utf16__indexS215 + 2;
            int32_t _M0L6_2atmpS1758 = _M0L11char__countS216 + 1;
            _M0L12utf16__indexS215 = _M0L6_2atmpS1757;
            _M0L11char__countS216 = _M0L6_2atmpS1758;
            continue;
          } else {
            #line 81 "/home/developer/.moon/lib/core/builtin/string.mbt"
            _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_182.data);
          }
        }
        _M0L6_2atmpS1760 = _M0L12utf16__indexS215 + 1;
        _M0L6_2atmpS1761 = _M0L11char__countS216 + 1;
        _M0L12utf16__indexS215 = _M0L6_2atmpS1760;
        _M0L11char__countS216 = _M0L6_2atmpS1761;
        continue;
      } else {
        return _M0L11char__countS216;
      }
      break;
    }
  } else {
    #line 70 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0FPC15abort5abortGiE((moonbit_string_t)moonbit_string_literal_183.data);
  }
}

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t _M0L4selfS209) {
  #line 45 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS209 >= 56320) {
    return _M0L4selfS209 <= 57343;
  } else {
    return 0;
  }
}

int32_t _M0MPC16uint166UInt1622is__leading__surrogate(int32_t _M0L4selfS208) {
  #line 28 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS208 >= 55296) {
    return _M0L4selfS208 <= 56319;
  } else {
    return 0;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder* _M0L4selfS206,
  int32_t _M0L2chS205
) {
  uint32_t _M0L4codeS204;
  #line 95 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  #line 96 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4codeS204 = _M0MPC14char4Char8to__uint(_M0L2chS205);
  if (_M0L4codeS204 <= 65535u) {
    int32_t _M0L3lenS1734 = _M0L4selfS206->$1;
    int32_t _M0L6_2atmpS1733 = _M0L3lenS1734 + 1;
    uint16_t* _M0L4dataS1735;
    int32_t _M0L3lenS1736;
    int32_t _M0L6_2atmpS1737;
    int32_t _M0L3lenS1739;
    int32_t _M0L6_2atmpS1738;
    #line 98 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS206, _M0L6_2atmpS1733);
    _M0L4dataS1735 = _M0L4selfS206->$0;
    _M0L3lenS1736 = _M0L4selfS206->$1;
    moonbit_incref(_M0L4dataS1735);
    #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1737 = _M0MPC14uint4UInt10to__uint16(_M0L4codeS204);
    if (
      _M0L3lenS1736 < 0
      || _M0L3lenS1736 >= Moonbit_array_length(_M0L4dataS1735)
    ) {
      #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1735[_M0L3lenS1736] = _M0L6_2atmpS1737;
    moonbit_decref(_M0L4dataS1735);
    _M0L3lenS1739 = _M0L4selfS206->$1;
    _M0L6_2atmpS1738 = _M0L3lenS1739 + 1;
    _M0L4selfS206->$1 = _M0L6_2atmpS1738;
  } else if (_M0L4codeS204 <= 1114111u) {
    int32_t _M0L3lenS1741 = _M0L4selfS206->$1;
    int32_t _M0L6_2atmpS1740 = _M0L3lenS1741 + 2;
    uint32_t _M0L4codeS207;
    uint16_t* _M0L4dataS1742;
    int32_t _M0L3lenS1743;
    uint32_t _M0L6_2atmpS1746;
    uint32_t _M0L6_2atmpS1745;
    int32_t _M0L6_2atmpS1744;
    uint16_t* _M0L4dataS1747;
    int32_t _M0L3lenS1752;
    int32_t _M0L6_2atmpS1748;
    uint32_t _M0L6_2atmpS1751;
    uint32_t _M0L6_2atmpS1750;
    int32_t _M0L6_2atmpS1749;
    int32_t _M0L3lenS1754;
    int32_t _M0L6_2atmpS1753;
    #line 102 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS206, _M0L6_2atmpS1740);
    _M0L4codeS207 = _M0L4codeS204 - 65536u;
    _M0L4dataS1742 = _M0L4selfS206->$0;
    _M0L3lenS1743 = _M0L4selfS206->$1;
    _M0L6_2atmpS1746 = _M0L4codeS207 >> 10;
    _M0L6_2atmpS1745 = 55296u + _M0L6_2atmpS1746;
    moonbit_incref(_M0L4dataS1742);
    #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1744 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1745);
    if (
      _M0L3lenS1743 < 0
      || _M0L3lenS1743 >= Moonbit_array_length(_M0L4dataS1742)
    ) {
      #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1742[_M0L3lenS1743] = _M0L6_2atmpS1744;
    moonbit_decref(_M0L4dataS1742);
    _M0L4dataS1747 = _M0L4selfS206->$0;
    _M0L3lenS1752 = _M0L4selfS206->$1;
    _M0L6_2atmpS1748 = _M0L3lenS1752 + 1;
    _M0L6_2atmpS1751 = _M0L4codeS207 & 1023u;
    _M0L6_2atmpS1750 = 56320u + _M0L6_2atmpS1751;
    moonbit_incref(_M0L4dataS1747);
    #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1749 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1750);
    if (
      _M0L6_2atmpS1748 < 0
      || _M0L6_2atmpS1748 >= Moonbit_array_length(_M0L4dataS1747)
    ) {
      #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1747[_M0L6_2atmpS1748] = _M0L6_2atmpS1749;
    moonbit_decref(_M0L4dataS1747);
    _M0L3lenS1754 = _M0L4selfS206->$1;
    _M0L6_2atmpS1753 = _M0L3lenS1754 + 2;
    _M0L4selfS206->$1 = _M0L6_2atmpS1753;
  } else {
    #line 108 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_184.data);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder* _M0L4selfS198,
  int32_t _M0L8requiredS199
) {
  uint16_t* _M0L4dataS1732;
  int32_t _M0L12current__lenS197;
  int32_t _M0L13enough__spaceS200;
  int32_t _M0L13enough__spaceS201;
  uint16_t* _M0L4dataS1728;
  int32_t _M0L6_2atmpS1729;
  int32_t _M0L3lenS1730;
  uint16_t* _M0L9new__dataS203;
  uint16_t* _M0L6_2aoldS3896;
  #line 46 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4dataS1732 = _M0L4selfS198->$0;
  _M0L12current__lenS197 = Moonbit_array_length(_M0L4dataS1732);
  if (_M0L8requiredS199 <= _M0L12current__lenS197) {
    return 0;
  }
  _M0L13enough__spaceS201 = _M0L12current__lenS197;
  while (1) {
    if (_M0L13enough__spaceS201 < _M0L8requiredS199) {
      int32_t _M0L6_2atmpS1731 = _M0L13enough__spaceS201 * 2;
      _M0L13enough__spaceS201 = _M0L6_2atmpS1731;
      continue;
    } else {
      _M0L13enough__spaceS200 = _M0L13enough__spaceS201;
    }
    break;
  }
  _M0L4dataS1728 = _M0L4selfS198->$0;
  moonbit_incref(_M0L4dataS1728);
  #line 64 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1729 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L3lenS1730 = _M0L4selfS198->$1;
  #line 61 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L9new__dataS203
  = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1728, _M0L13enough__spaceS200, _M0L6_2atmpS1729, _M0L3lenS1730, 0, 0);
  moonbit_decref(_M0L4dataS1728);
  _M0L6_2aoldS3896 = _M0L4selfS198->$0;
  moonbit_decref(_M0L6_2aoldS3896);
  _M0L4selfS198->$0 = _M0L9new__dataS203;
  return 0;
}

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t _M0L4selfS196) {
  int32_t _M0L6_2atmpS1727;
  #line 2676 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1727 = *(int32_t*)&_M0L4selfS196;
  return (uint16_t)_M0L6_2atmpS1727;
}

uint32_t _M0MPC14char4Char8to__uint(int32_t _M0L4selfS195) {
  int32_t _M0L6_2atmpS1726;
  #line 1254 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1726 = _M0L4selfS195;
  return *(uint32_t*)&_M0L6_2atmpS1726;
}

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder* _M0L4selfS193
) {
  int32_t _M0L3lenS1718;
  uint16_t* _M0L4dataS1720;
  int32_t _M0L6_2atmpS1719;
  #line 148 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3lenS1718 = _M0L4selfS193->$1;
  _M0L4dataS1720 = _M0L4selfS193->$0;
  _M0L6_2atmpS1719 = Moonbit_array_length(_M0L4dataS1720);
  if (_M0L3lenS1718 == _M0L6_2atmpS1719) {
    uint16_t* _M0L4dataS1721 = _M0L4selfS193->$0;
    moonbit_incref(_M0L4dataS1721);
    return _M0L4dataS1721;
  } else {
    uint16_t* _M0L4dataS1722 = _M0L4selfS193->$0;
    int32_t _M0L3lenS1723 = _M0L4selfS193->$1;
    int32_t _M0L6_2atmpS1724;
    int32_t _M0L3lenS1725;
    uint16_t* _M0L4dataS194;
    moonbit_incref(_M0L4dataS1722);
    #line 155 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1724 = _M0IPC16uint166UInt16PB7Default7default();
    _M0L3lenS1725 = _M0L4selfS193->$1;
    #line 152 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L4dataS194
    = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1722, _M0L3lenS1723, _M0L6_2atmpS1724, _M0L3lenS1725, 0, 0);
    moonbit_decref(_M0L4dataS1722);
    return _M0L4dataS194;
  }
}

int32_t _M0IPC16uint166UInt16PB7Default7default() {
  #line 176 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  return 0;
}

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t* _M0L3srcS190,
  int32_t _M0L13allocate__lenS186,
  int32_t _M0L4initS191,
  int32_t _M0L3lenS187,
  int32_t _M0L11src__offsetS188,
  int32_t _M0L11dst__offsetS189
) {
  int32_t _if__result_4607;
  #line 97 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L13allocate__lenS186 >= 0) {
    if (_M0L3lenS187 >= 0) {
      if (_M0L11src__offsetS188 >= 0) {
        if (_M0L11dst__offsetS189 >= 0) {
          int32_t _M0L6_2atmpS1714 = _M0L11src__offsetS188 + _M0L3lenS187;
          int32_t _M0L6_2atmpS1715 = Moonbit_array_length(_M0L3srcS190);
          if (_M0L6_2atmpS1714 <= _M0L6_2atmpS1715) {
            int32_t _M0L6_2atmpS1713 = _M0L11dst__offsetS189 + _M0L3lenS187;
            _if__result_4607 = _M0L6_2atmpS1713 <= _M0L13allocate__lenS186;
          } else {
            _if__result_4607 = 0;
          }
        } else {
          _if__result_4607 = 0;
        }
      } else {
        _if__result_4607 = 0;
      }
    } else {
      _if__result_4607 = 0;
    }
  } else {
    _if__result_4607 = 0;
  }
  if (_if__result_4607) {
    moonbit_incref(_M0L3srcS190);
    #line 115 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    return _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(_M0L3srcS190, _M0L13allocate__lenS186, _M0L4initS191, _M0L11src__offsetS188, _M0L11dst__offsetS189, _M0L3lenS187);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS192;
    int32_t _M0L6_2atmpS1717;
    moonbit_string_t _M0L6_2atmpS1716;
    uint16_t* _result_4608;
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L18_2astring__builderS192
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS192, (moonbit_string_t)moonbit_string_literal_185.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS192, _M0L13allocate__lenS186);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS192, (moonbit_string_t)moonbit_string_literal_186.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS192, _M0L11src__offsetS188);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS192, (moonbit_string_t)moonbit_string_literal_187.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS192, _M0L11dst__offsetS189);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS192, (moonbit_string_t)moonbit_string_literal_188.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS192, _M0L3lenS187);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS192, (moonbit_string_t)moonbit_string_literal_189.data);
    _M0L6_2atmpS1717 = Moonbit_array_length(_M0L3srcS190);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS192, _M0L6_2atmpS1717);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L6_2atmpS1716
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS192);
    moonbit_decref(_M0L18_2astring__builderS192);
    #line 111 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _result_4608 = _M0FPC15abort5abortGAkE(_M0L6_2atmpS1716);
    moonbit_decref(_M0L6_2atmpS1716);
    return _result_4608;
  }
}

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t* _M0L3srcS183,
  int32_t _M0L13allocate__lenS180,
  int32_t _M0L4initS181,
  int32_t _M0L11src__offsetS184,
  int32_t _M0L11dst__offsetS182,
  int32_t _M0L9blit__lenS185
) {
  uint16_t* _M0L3dstS179;
  #line 79 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  _M0L3dstS179
  = (uint16_t*)moonbit_make_string(_M0L13allocate__lenS180, _M0L4initS181);
  moonbit_incref(_M0L3dstS179);
  #line 90 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS179, _M0L11dst__offsetS182, _M0L3srcS183, _M0L11src__offsetS184, _M0L9blit__lenS185, sizeof(uint16_t));
  return _M0L3dstS179;
}

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t _M0L10size__hintS177
) {
  int32_t _M0L7initialS176;
  uint16_t* _M0L4dataS178;
  struct _M0TPB13StringBuilder* _block_4609;
  #line 32 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  if (_M0L10size__hintS177 < 1) {
    _M0L7initialS176 = 1;
  } else {
    int32_t _M0L6_2atmpS1712 = _M0L10size__hintS177 + 1;
    _M0L7initialS176 = _M0L6_2atmpS1712 / 2;
  }
  _M0L4dataS178 = (uint16_t*)moonbit_make_string(_M0L7initialS176, 0);
  _block_4609
  = (struct _M0TPB13StringBuilder*)moonbit_malloc(sizeof(struct _M0TPB13StringBuilder));
  Moonbit_object_header(_block_4609)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 159, 0);
  _block_4609->$0 = _M0L4dataS178;
  _block_4609->$1 = 0;
  return _block_4609;
}

int32_t _M0MPC14byte4Byte8to__char(int32_t _M0L4selfS175) {
  int32_t _M0L6_2atmpS1711;
  #line 1868 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1711 = (int32_t)_M0L4selfS175;
  return _M0L6_2atmpS1711;
}

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t* _M0L3srcS149,
  int32_t _M0L13allocate__lenS145,
  int32_t _M0L3lenS146,
  int32_t _M0L11src__offsetS147,
  int32_t _M0L11dst__offsetS148
) {
  int32_t _if__result_4610;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS145 >= 0) {
    if (_M0L3lenS146 >= 0) {
      if (_M0L11src__offsetS147 >= 0) {
        if (_M0L11dst__offsetS148 >= 0) {
          int32_t _M0L6_2atmpS1687 = _M0L11src__offsetS147 + _M0L3lenS146;
          int32_t _M0L6_2atmpS1688;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1688
          = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS149);
          if (_M0L6_2atmpS1687 <= _M0L6_2atmpS1688) {
            int32_t _M0L6_2atmpS1686 = _M0L11dst__offsetS148 + _M0L3lenS146;
            _if__result_4610 = _M0L6_2atmpS1686 <= _M0L13allocate__lenS145;
          } else {
            _if__result_4610 = 0;
          }
        } else {
          _if__result_4610 = 0;
        }
      } else {
        _if__result_4610 = 0;
      }
    } else {
      _if__result_4610 = 0;
    }
  } else {
    _if__result_4610 = 0;
  }
  if (_if__result_4610) {
    moonbit_incref(_M0L3srcS149);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (moonbit_string_t*)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS145, (moonbit_string_t)moonbit_string_literal_0.data, _M0L3srcS149, _M0L11src__offsetS147, _M0L11dst__offsetS148, _M0L3lenS146);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS150;
    int32_t _M0L6_2atmpS1690;
    moonbit_string_t _M0L6_2atmpS1689;
    moonbit_string_t* _result_4611;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS150
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_185.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L13allocate__lenS145);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_186.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L11src__offsetS147);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_187.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L11dst__offsetS148);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_188.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L3lenS146);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS150, (moonbit_string_t)moonbit_string_literal_189.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1690 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS149);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS150, _M0L6_2atmpS1690);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1689
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS150);
    moonbit_decref(_M0L18_2astring__builderS150);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4611
    = _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(_M0L6_2atmpS1689);
    moonbit_decref(_M0L6_2atmpS1689);
    return _result_4611;
  }
}

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE** _M0L3srcS155,
  int32_t _M0L13allocate__lenS151,
  int32_t _M0L3lenS152,
  int32_t _M0L11src__offsetS153,
  int32_t _M0L11dst__offsetS154
) {
  int32_t _if__result_4612;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS151 >= 0) {
    if (_M0L3lenS152 >= 0) {
      if (_M0L11src__offsetS153 >= 0) {
        if (_M0L11dst__offsetS154 >= 0) {
          int32_t _M0L6_2atmpS1692 = _M0L11src__offsetS153 + _M0L3lenS152;
          int32_t _M0L6_2atmpS1693;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1693
          = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS155);
          if (_M0L6_2atmpS1692 <= _M0L6_2atmpS1693) {
            int32_t _M0L6_2atmpS1691 = _M0L11dst__offsetS154 + _M0L3lenS152;
            _if__result_4612 = _M0L6_2atmpS1691 <= _M0L13allocate__lenS151;
          } else {
            _if__result_4612 = 0;
          }
        } else {
          _if__result_4612 = 0;
        }
      } else {
        _if__result_4612 = 0;
      }
    } else {
      _if__result_4612 = 0;
    }
  } else {
    _if__result_4612 = 0;
  }
  if (_if__result_4612) {
    moonbit_incref(_M0L3srcS155);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TUsiE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS151, 0, _M0L3srcS155, _M0L11src__offsetS153, _M0L11dst__offsetS154, _M0L3lenS152);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS156;
    int32_t _M0L6_2atmpS1695;
    moonbit_string_t _M0L6_2atmpS1694;
    struct _M0TUsiE** _result_4613;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS156
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_185.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L13allocate__lenS151);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_186.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L11src__offsetS153);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_187.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L11dst__offsetS154);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_188.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L3lenS152);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS156, (moonbit_string_t)moonbit_string_literal_189.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1695 = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS155);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS156, _M0L6_2atmpS1695);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1694
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS156);
    moonbit_decref(_M0L18_2astring__builderS156);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4613
    = _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(_M0L6_2atmpS1694);
    moonbit_decref(_M0L6_2atmpS1694);
    return _result_4613;
  }
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS161,
  int32_t _M0L13allocate__lenS157,
  int32_t _M0L3lenS158,
  int32_t _M0L11src__offsetS159,
  int32_t _M0L11dst__offsetS160
) {
  int32_t _if__result_4614;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS157 >= 0) {
    if (_M0L3lenS158 >= 0) {
      if (_M0L11src__offsetS159 >= 0) {
        if (_M0L11dst__offsetS160 >= 0) {
          int32_t _M0L6_2atmpS1697 = _M0L11src__offsetS159 + _M0L3lenS158;
          int32_t _M0L6_2atmpS1698;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1698
          = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L3srcS161);
          if (_M0L6_2atmpS1697 <= _M0L6_2atmpS1698) {
            int32_t _M0L6_2atmpS1696 = _M0L11dst__offsetS160 + _M0L3lenS158;
            _if__result_4614 = _M0L6_2atmpS1696 <= _M0L13allocate__lenS157;
          } else {
            _if__result_4614 = 0;
          }
        } else {
          _if__result_4614 = 0;
        }
      } else {
        _if__result_4614 = 0;
      }
    } else {
      _if__result_4614 = 0;
    }
  } else {
    _if__result_4614 = 0;
  }
  if (_if__result_4614) {
    moonbit_incref(_M0L3srcS161);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS157, 0, _M0L3srcS161, _M0L11src__offsetS159, _M0L11dst__offsetS160, _M0L3lenS158);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS162;
    int32_t _M0L6_2atmpS1700;
    moonbit_string_t _M0L6_2atmpS1699;
    struct _M0TP26biolab8bio__seq9SeqRecord** _result_4615;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS162
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_185.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L13allocate__lenS157);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_186.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L11src__offsetS159);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_187.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L11dst__offsetS160);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_188.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L3lenS158);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS162, (moonbit_string_t)moonbit_string_literal_189.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1700
    = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L3srcS161);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS162, _M0L6_2atmpS1700);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1699
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS162);
    moonbit_decref(_M0L18_2astring__builderS162);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4615
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(_M0L6_2atmpS1699);
    moonbit_decref(_M0L6_2atmpS1699);
    return _result_4615;
  }
}

int32_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGiE(
  int32_t* _M0L3srcS167,
  int32_t _M0L13allocate__lenS163,
  int32_t _M0L3lenS164,
  int32_t _M0L11src__offsetS165,
  int32_t _M0L11dst__offsetS166
) {
  int32_t _if__result_4616;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS163 >= 0) {
    if (_M0L3lenS164 >= 0) {
      if (_M0L11src__offsetS165 >= 0) {
        if (_M0L11dst__offsetS166 >= 0) {
          int32_t _M0L6_2atmpS1702 = _M0L11src__offsetS165 + _M0L3lenS164;
          int32_t _M0L6_2atmpS1703;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1703
          = _M0MPB18UninitializedArray6lengthGiE(_M0L3srcS167);
          if (_M0L6_2atmpS1702 <= _M0L6_2atmpS1703) {
            int32_t _M0L6_2atmpS1701 = _M0L11dst__offsetS166 + _M0L3lenS164;
            _if__result_4616 = _M0L6_2atmpS1701 <= _M0L13allocate__lenS163;
          } else {
            _if__result_4616 = 0;
          }
        } else {
          _if__result_4616 = 0;
        }
      } else {
        _if__result_4616 = 0;
      }
    } else {
      _if__result_4616 = 0;
    }
  } else {
    _if__result_4616 = 0;
  }
  if (_if__result_4616) {
    moonbit_incref(_M0L3srcS167);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return _M0MPB18UninitializedArray23unsafe__make__and__blitGiE(_M0L3srcS167, _M0L13allocate__lenS163, _M0L11src__offsetS165, _M0L11dst__offsetS166, _M0L3lenS164);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS168;
    int32_t _M0L6_2atmpS1705;
    moonbit_string_t _M0L6_2atmpS1704;
    int32_t* _result_4617;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS168
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_185.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L13allocate__lenS163);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_186.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L11src__offsetS165);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_187.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L11dst__offsetS166);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_188.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L3lenS164);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS168, (moonbit_string_t)moonbit_string_literal_189.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1705 = _M0MPB18UninitializedArray6lengthGiE(_M0L3srcS167);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS168, _M0L6_2atmpS1705);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1704
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS168);
    moonbit_decref(_M0L18_2astring__builderS168);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4617
    = _M0FPC15abort5abortGRPB18UninitializedArrayGiEE(_M0L6_2atmpS1704);
    moonbit_decref(_M0L6_2atmpS1704);
    return _result_4617;
  }
}

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3srcS173,
  int32_t _M0L13allocate__lenS169,
  int32_t _M0L3lenS170,
  int32_t _M0L11src__offsetS171,
  int32_t _M0L11dst__offsetS172
) {
  int32_t _if__result_4618;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS169 >= 0) {
    if (_M0L3lenS170 >= 0) {
      if (_M0L11src__offsetS171 >= 0) {
        if (_M0L11dst__offsetS172 >= 0) {
          int32_t _M0L6_2atmpS1707 = _M0L11src__offsetS171 + _M0L3lenS170;
          int32_t _M0L6_2atmpS1708;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1708
          = _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(_M0L3srcS173);
          if (_M0L6_2atmpS1707 <= _M0L6_2atmpS1708) {
            int32_t _M0L6_2atmpS1706 = _M0L11dst__offsetS172 + _M0L3lenS170;
            _if__result_4618 = _M0L6_2atmpS1706 <= _M0L13allocate__lenS169;
          } else {
            _if__result_4618 = 0;
          }
        } else {
          _if__result_4618 = 0;
        }
      } else {
        _if__result_4618 = 0;
      }
    } else {
      _if__result_4618 = 0;
    }
  } else {
    _if__result_4618 = 0;
  }
  if (_if__result_4618) {
    moonbit_incref(_M0L3srcS173);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return _M0MPB18UninitializedArray23unsafe__make__and__blitGRPC16string10StringViewE(_M0L3srcS173, _M0L13allocate__lenS169, _M0L11src__offsetS171, _M0L11dst__offsetS172, _M0L3lenS170);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS174;
    int32_t _M0L6_2atmpS1710;
    moonbit_string_t _M0L6_2atmpS1709;
    struct _M0TPC16string10StringView* _result_4619;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS174
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS174, (moonbit_string_t)moonbit_string_literal_185.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS174, _M0L13allocate__lenS169);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS174, (moonbit_string_t)moonbit_string_literal_186.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS174, _M0L11src__offsetS171);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS174, (moonbit_string_t)moonbit_string_literal_187.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS174, _M0L11dst__offsetS172);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS174, (moonbit_string_t)moonbit_string_literal_188.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS174, _M0L3lenS170);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS174, (moonbit_string_t)moonbit_string_literal_189.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1710
    = _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(_M0L3srcS173);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS174, _M0L6_2atmpS1710);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1709
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS174);
    moonbit_decref(_M0L18_2astring__builderS174);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4619
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRPC16string10StringViewEE(_M0L6_2atmpS1709);
    moonbit_decref(_M0L6_2atmpS1709);
    return _result_4619;
  }
}

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder* _M0L4selfS136,
  moonbit_string_t _M0L3objS135
) {
  struct _M0TPB6Logger _M0L6_2atmpS1681;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS136);
  _M0L6_2atmpS1681
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS136
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS135, _M0L6_2atmpS1681);
  if (_M0L6_2atmpS1681.$1) {
    moonbit_decref(_M0L6_2atmpS1681.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder* _M0L4selfS138,
  int32_t _M0L3objS137
) {
  struct _M0TPB6Logger _M0L6_2atmpS1682;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS138);
  _M0L6_2atmpS1682
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS138
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L3objS137, _M0L6_2atmpS1682);
  if (_M0L6_2atmpS1682.$1) {
    moonbit_decref(_M0L6_2atmpS1682.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGbE(
  struct _M0TPB13StringBuilder* _M0L4selfS140,
  int32_t _M0L3objS139
) {
  struct _M0TPB6Logger _M0L6_2atmpS1683;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS140);
  _M0L6_2atmpS1683
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS140
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGbE(_M0L3objS139, _M0L6_2atmpS1683);
  if (_M0L6_2atmpS1683.$1) {
    moonbit_decref(_M0L6_2atmpS1683.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGRPB9SourceLocE(
  struct _M0TPB13StringBuilder* _M0L4selfS142,
  moonbit_string_t _M0L3objS141
) {
  struct _M0TPB6Logger _M0L6_2atmpS1684;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS142);
  _M0L6_2atmpS1684
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS142
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IPB9SourceLocPB4Show6output(_M0L3objS141, _M0L6_2atmpS1684);
  if (_M0L6_2atmpS1684.$1) {
    moonbit_decref(_M0L6_2atmpS1684.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGRPC16string10StringViewE(
  struct _M0TPB13StringBuilder* _M0L4selfS144,
  struct _M0TPC16string10StringView _M0L3objS143
) {
  struct _M0TPB6Logger _M0L6_2atmpS1685;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS144);
  _M0L6_2atmpS1685
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS144
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IPC16string10StringViewPB4Show6output(_M0L3objS143, _M0L6_2atmpS1685);
  if (_M0L6_2atmpS1685.$1) {
    moonbit_decref(_M0L6_2atmpS1685.$1);
  }
  return 0;
}

moonbit_string_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGsE(
  moonbit_string_t* _M0L3srcS108,
  int32_t _M0L13allocate__lenS106,
  int32_t _M0L11src__offsetS109,
  int32_t _M0L11dst__offsetS107,
  int32_t _M0L9blit__lenS110
) {
  moonbit_string_t* _M0L3dstS105;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS105
  = (moonbit_string_t*)moonbit_make_ref_array(_M0L13allocate__lenS106, (moonbit_string_t)moonbit_string_literal_0.data);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGsE(_M0L3dstS105, _M0L11dst__offsetS107, _M0L3srcS108, _M0L11src__offsetS109, _M0L9blit__lenS110);
  moonbit_decref(_M0L3srcS108);
  return _M0L3dstS105;
}

struct _M0TUsiE** _M0MPB18UninitializedArray23unsafe__make__and__blitGUsiEE(
  struct _M0TUsiE** _M0L3srcS114,
  int32_t _M0L13allocate__lenS112,
  int32_t _M0L11src__offsetS115,
  int32_t _M0L11dst__offsetS113,
  int32_t _M0L9blit__lenS116
) {
  struct _M0TUsiE** _M0L3dstS111;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS111
  = (struct _M0TUsiE**)moonbit_make_ref_array(_M0L13allocate__lenS112, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGUsiEE(_M0L3dstS111, _M0L11dst__offsetS113, _M0L3srcS114, _M0L11src__offsetS115, _M0L9blit__lenS116);
  moonbit_decref(_M0L3srcS114);
  return _M0L3dstS111;
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS120,
  int32_t _M0L13allocate__lenS118,
  int32_t _M0L11src__offsetS121,
  int32_t _M0L11dst__offsetS119,
  int32_t _M0L9blit__lenS122
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS117;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS117
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array(_M0L13allocate__lenS118, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(_M0L3dstS117, _M0L11dst__offsetS119, _M0L3srcS120, _M0L11src__offsetS121, _M0L9blit__lenS122);
  moonbit_decref(_M0L3srcS120);
  return _M0L3dstS117;
}

int32_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGiE(
  int32_t* _M0L3srcS126,
  int32_t _M0L13allocate__lenS124,
  int32_t _M0L11src__offsetS127,
  int32_t _M0L11dst__offsetS125,
  int32_t _M0L9blit__lenS128
) {
  int32_t* _M0L3dstS123;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS123
  = (int32_t*)moonbit_make_int32_array_raw(_M0L13allocate__lenS124);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGiE(_M0L3dstS123, _M0L11dst__offsetS125, _M0L3srcS126, _M0L11src__offsetS127, _M0L9blit__lenS128);
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
  = (struct _M0TPC16string10StringView*)moonbit_make_ref_valtype_array(_M0L13allocate__lenS130, sizeof(struct _M0TPC16string10StringView), Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 114, 0), &(struct _M0TPC16string10StringView){.$0 = (moonbit_string_t)moonbit_string_literal_0.data, .$1 = 0, .$2 = 0});
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRPC16string10StringViewE(_M0L3dstS129, _M0L11dst__offsetS131, _M0L3srcS132, _M0L11src__offsetS133, _M0L9blit__lenS134);
  moonbit_decref(_M0L3srcS132);
  return _M0L3dstS129;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGsE(
  moonbit_string_t* _M0L3dstS80,
  int32_t _M0L11dst__offsetS81,
  moonbit_string_t* _M0L3srcS82,
  int32_t _M0L11src__offsetS83,
  int32_t _M0L3lenS84
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS82);
  moonbit_incref(_M0L3dstS80);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS80, _M0L11dst__offsetS81, _M0L3srcS82, _M0L11src__offsetS83, _M0L3lenS84);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGUsiEE(
  struct _M0TUsiE** _M0L3dstS85,
  int32_t _M0L11dst__offsetS86,
  struct _M0TUsiE** _M0L3srcS87,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS90,
  int32_t _M0L11dst__offsetS91,
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS92,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGiE(
  int32_t* _M0L3dstS95,
  int32_t _M0L11dst__offsetS96,
  int32_t* _M0L3srcS97,
  int32_t _M0L11src__offsetS98,
  int32_t _M0L3lenS99
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS97);
  moonbit_incref(_M0L3dstS95);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS95, _M0L11dst__offsetS96, _M0L3srcS97, _M0L11src__offsetS98, _M0L3lenS99, sizeof(int32_t));
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
  int32_t _if__result_4620;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS26 == _M0L3srcS27) {
    _if__result_4620 = _M0L11dst__offsetS28 < _M0L11src__offsetS29;
  } else {
    _if__result_4620 = 0;
  }
  if (_if__result_4620) {
    int32_t _M0L1iS30 = 0;
    while (1) {
      if (_M0L1iS30 < _M0L3lenS31) {
        int32_t _M0L6_2atmpS1627 = _M0L11dst__offsetS28 + _M0L1iS30;
        int32_t _M0L6_2atmpS1629 = _M0L11src__offsetS29 + _M0L1iS30;
        int32_t _M0L6_2atmpS1628;
        int32_t _M0L6_2atmpS1630;
        if (
          _M0L6_2atmpS1629 < 0
          || _M0L6_2atmpS1629 >= Moonbit_array_length(_M0L3srcS27)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1628 = (int32_t)_M0L3srcS27[_M0L6_2atmpS1629];
        if (
          _M0L6_2atmpS1627 < 0
          || _M0L6_2atmpS1627 >= Moonbit_array_length(_M0L3dstS26)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS26[_M0L6_2atmpS1627] = _M0L6_2atmpS1628;
        _M0L6_2atmpS1630 = _M0L1iS30 + 1;
        _M0L1iS30 = _M0L6_2atmpS1630;
        continue;
      } else {
        moonbit_decref(_M0L3srcS27);
        moonbit_decref(_M0L3dstS26);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1635 = _M0L3lenS31 - 1;
    int32_t _M0L1iS33 = _M0L6_2atmpS1635;
    while (1) {
      if (_M0L1iS33 >= 0) {
        int32_t _M0L6_2atmpS1631 = _M0L11dst__offsetS28 + _M0L1iS33;
        int32_t _M0L6_2atmpS1633 = _M0L11src__offsetS29 + _M0L1iS33;
        int32_t _M0L6_2atmpS1632;
        int32_t _M0L6_2atmpS1634;
        if (
          _M0L6_2atmpS1633 < 0
          || _M0L6_2atmpS1633 >= Moonbit_array_length(_M0L3srcS27)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1632 = (int32_t)_M0L3srcS27[_M0L6_2atmpS1633];
        if (
          _M0L6_2atmpS1631 < 0
          || _M0L6_2atmpS1631 >= Moonbit_array_length(_M0L3dstS26)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS26[_M0L6_2atmpS1631] = _M0L6_2atmpS1632;
        _M0L6_2atmpS1634 = _M0L1iS33 - 1;
        _M0L1iS33 = _M0L6_2atmpS1634;
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGsEE(
  moonbit_string_t* _M0L3dstS35,
  int32_t _M0L11dst__offsetS37,
  moonbit_string_t* _M0L3srcS36,
  int32_t _M0L11src__offsetS38,
  int32_t _M0L3lenS40
) {
  int32_t _if__result_4623;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS35 == _M0L3srcS36) {
    _if__result_4623 = _M0L11dst__offsetS37 < _M0L11src__offsetS38;
  } else {
    _if__result_4623 = 0;
  }
  if (_if__result_4623) {
    int32_t _M0L1iS39 = 0;
    while (1) {
      if (_M0L1iS39 < _M0L3lenS40) {
        int32_t _M0L6_2atmpS1636 = _M0L11dst__offsetS37 + _M0L1iS39;
        int32_t _M0L6_2atmpS1638 = _M0L11src__offsetS38 + _M0L1iS39;
        moonbit_string_t _M0L6_2atmpS1637;
        moonbit_string_t _M0L6_2aoldS3902;
        int32_t _M0L6_2atmpS1639;
        if (
          _M0L6_2atmpS1638 < 0
          || _M0L6_2atmpS1638 >= Moonbit_array_length(_M0L3srcS36)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1637 = (moonbit_string_t)_M0L3srcS36[_M0L6_2atmpS1638];
        if (
          _M0L6_2atmpS1636 < 0
          || _M0L6_2atmpS1636 >= Moonbit_array_length(_M0L3dstS35)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3902 = (moonbit_string_t)_M0L3dstS35[_M0L6_2atmpS1636];
        moonbit_incref(_M0L6_2atmpS1637);
        moonbit_decref(_M0L6_2aoldS3902);
        _M0L3dstS35[_M0L6_2atmpS1636] = _M0L6_2atmpS1637;
        _M0L6_2atmpS1639 = _M0L1iS39 + 1;
        _M0L1iS39 = _M0L6_2atmpS1639;
        continue;
      } else {
        moonbit_decref(_M0L3srcS36);
        moonbit_decref(_M0L3dstS35);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1644 = _M0L3lenS40 - 1;
    int32_t _M0L1iS42 = _M0L6_2atmpS1644;
    while (1) {
      if (_M0L1iS42 >= 0) {
        int32_t _M0L6_2atmpS1640 = _M0L11dst__offsetS37 + _M0L1iS42;
        int32_t _M0L6_2atmpS1642 = _M0L11src__offsetS38 + _M0L1iS42;
        moonbit_string_t _M0L6_2atmpS1641;
        moonbit_string_t _M0L6_2aoldS3904;
        int32_t _M0L6_2atmpS1643;
        if (
          _M0L6_2atmpS1642 < 0
          || _M0L6_2atmpS1642 >= Moonbit_array_length(_M0L3srcS36)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1641 = (moonbit_string_t)_M0L3srcS36[_M0L6_2atmpS1642];
        if (
          _M0L6_2atmpS1640 < 0
          || _M0L6_2atmpS1640 >= Moonbit_array_length(_M0L3dstS35)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3904 = (moonbit_string_t)_M0L3dstS35[_M0L6_2atmpS1640];
        moonbit_incref(_M0L6_2atmpS1641);
        moonbit_decref(_M0L6_2aoldS3904);
        _M0L3dstS35[_M0L6_2atmpS1640] = _M0L6_2atmpS1641;
        _M0L6_2atmpS1643 = _M0L1iS42 - 1;
        _M0L1iS42 = _M0L6_2atmpS1643;
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGUsiEEE(
  struct _M0TUsiE** _M0L3dstS44,
  int32_t _M0L11dst__offsetS46,
  struct _M0TUsiE** _M0L3srcS45,
  int32_t _M0L11src__offsetS47,
  int32_t _M0L3lenS49
) {
  int32_t _if__result_4626;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS44 == _M0L3srcS45) {
    _if__result_4626 = _M0L11dst__offsetS46 < _M0L11src__offsetS47;
  } else {
    _if__result_4626 = 0;
  }
  if (_if__result_4626) {
    int32_t _M0L1iS48 = 0;
    while (1) {
      if (_M0L1iS48 < _M0L3lenS49) {
        int32_t _M0L6_2atmpS1645 = _M0L11dst__offsetS46 + _M0L1iS48;
        int32_t _M0L6_2atmpS1647 = _M0L11src__offsetS47 + _M0L1iS48;
        struct _M0TUsiE* _M0L6_2atmpS1646;
        struct _M0TUsiE* _M0L6_2aoldS3906;
        int32_t _M0L6_2atmpS1648;
        if (
          _M0L6_2atmpS1647 < 0
          || _M0L6_2atmpS1647 >= Moonbit_array_length(_M0L3srcS45)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1646 = (struct _M0TUsiE*)_M0L3srcS45[_M0L6_2atmpS1647];
        if (
          _M0L6_2atmpS1645 < 0
          || _M0L6_2atmpS1645 >= Moonbit_array_length(_M0L3dstS44)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3906 = (struct _M0TUsiE*)_M0L3dstS44[_M0L6_2atmpS1645];
        if (_M0L6_2atmpS1646) {
          moonbit_incref(_M0L6_2atmpS1646);
        }
        if (_M0L6_2aoldS3906) {
          moonbit_decref(_M0L6_2aoldS3906);
        }
        _M0L3dstS44[_M0L6_2atmpS1645] = _M0L6_2atmpS1646;
        _M0L6_2atmpS1648 = _M0L1iS48 + 1;
        _M0L1iS48 = _M0L6_2atmpS1648;
        continue;
      } else {
        moonbit_decref(_M0L3srcS45);
        moonbit_decref(_M0L3dstS44);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1653 = _M0L3lenS49 - 1;
    int32_t _M0L1iS51 = _M0L6_2atmpS1653;
    while (1) {
      if (_M0L1iS51 >= 0) {
        int32_t _M0L6_2atmpS1649 = _M0L11dst__offsetS46 + _M0L1iS51;
        int32_t _M0L6_2atmpS1651 = _M0L11src__offsetS47 + _M0L1iS51;
        struct _M0TUsiE* _M0L6_2atmpS1650;
        struct _M0TUsiE* _M0L6_2aoldS3908;
        int32_t _M0L6_2atmpS1652;
        if (
          _M0L6_2atmpS1651 < 0
          || _M0L6_2atmpS1651 >= Moonbit_array_length(_M0L3srcS45)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1650 = (struct _M0TUsiE*)_M0L3srcS45[_M0L6_2atmpS1651];
        if (
          _M0L6_2atmpS1649 < 0
          || _M0L6_2atmpS1649 >= Moonbit_array_length(_M0L3dstS44)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3908 = (struct _M0TUsiE*)_M0L3dstS44[_M0L6_2atmpS1649];
        if (_M0L6_2atmpS1650) {
          moonbit_incref(_M0L6_2atmpS1650);
        }
        if (_M0L6_2aoldS3908) {
          moonbit_decref(_M0L6_2aoldS3908);
        }
        _M0L3dstS44[_M0L6_2atmpS1649] = _M0L6_2atmpS1650;
        _M0L6_2atmpS1652 = _M0L1iS51 - 1;
        _M0L1iS51 = _M0L6_2atmpS1652;
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS53,
  int32_t _M0L11dst__offsetS55,
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS54,
  int32_t _M0L11src__offsetS56,
  int32_t _M0L3lenS58
) {
  int32_t _if__result_4629;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS53 == _M0L3srcS54) {
    _if__result_4629 = _M0L11dst__offsetS55 < _M0L11src__offsetS56;
  } else {
    _if__result_4629 = 0;
  }
  if (_if__result_4629) {
    int32_t _M0L1iS57 = 0;
    while (1) {
      if (_M0L1iS57 < _M0L3lenS58) {
        int32_t _M0L6_2atmpS1654 = _M0L11dst__offsetS55 + _M0L1iS57;
        int32_t _M0L6_2atmpS1656 = _M0L11src__offsetS56 + _M0L1iS57;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1655;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3910;
        int32_t _M0L6_2atmpS1657;
        if (
          _M0L6_2atmpS1656 < 0
          || _M0L6_2atmpS1656 >= Moonbit_array_length(_M0L3srcS54)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1655
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3srcS54[
            _M0L6_2atmpS1656
          ];
        if (
          _M0L6_2atmpS1654 < 0
          || _M0L6_2atmpS1654 >= Moonbit_array_length(_M0L3dstS53)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3910
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3dstS53[
            _M0L6_2atmpS1654
          ];
        if (_M0L6_2atmpS1655) {
          moonbit_incref(_M0L6_2atmpS1655);
        }
        if (_M0L6_2aoldS3910) {
          moonbit_decref(_M0L6_2aoldS3910);
        }
        _M0L3dstS53[_M0L6_2atmpS1654] = _M0L6_2atmpS1655;
        _M0L6_2atmpS1657 = _M0L1iS57 + 1;
        _M0L1iS57 = _M0L6_2atmpS1657;
        continue;
      } else {
        moonbit_decref(_M0L3srcS54);
        moonbit_decref(_M0L3dstS53);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1662 = _M0L3lenS58 - 1;
    int32_t _M0L1iS60 = _M0L6_2atmpS1662;
    while (1) {
      if (_M0L1iS60 >= 0) {
        int32_t _M0L6_2atmpS1658 = _M0L11dst__offsetS55 + _M0L1iS60;
        int32_t _M0L6_2atmpS1660 = _M0L11src__offsetS56 + _M0L1iS60;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1659;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3912;
        int32_t _M0L6_2atmpS1661;
        if (
          _M0L6_2atmpS1660 < 0
          || _M0L6_2atmpS1660 >= Moonbit_array_length(_M0L3srcS54)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1659
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3srcS54[
            _M0L6_2atmpS1660
          ];
        if (
          _M0L6_2atmpS1658 < 0
          || _M0L6_2atmpS1658 >= Moonbit_array_length(_M0L3dstS53)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3912
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3dstS53[
            _M0L6_2atmpS1658
          ];
        if (_M0L6_2atmpS1659) {
          moonbit_incref(_M0L6_2atmpS1659);
        }
        if (_M0L6_2aoldS3912) {
          moonbit_decref(_M0L6_2aoldS3912);
        }
        _M0L3dstS53[_M0L6_2atmpS1658] = _M0L6_2atmpS1659;
        _M0L6_2atmpS1661 = _M0L1iS60 - 1;
        _M0L1iS60 = _M0L6_2atmpS1661;
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGiEE(
  int32_t* _M0L3dstS62,
  int32_t _M0L11dst__offsetS64,
  int32_t* _M0L3srcS63,
  int32_t _M0L11src__offsetS65,
  int32_t _M0L3lenS67
) {
  int32_t _if__result_4632;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS62 == _M0L3srcS63) {
    _if__result_4632 = _M0L11dst__offsetS64 < _M0L11src__offsetS65;
  } else {
    _if__result_4632 = 0;
  }
  if (_if__result_4632) {
    int32_t _M0L1iS66 = 0;
    while (1) {
      if (_M0L1iS66 < _M0L3lenS67) {
        int32_t _M0L6_2atmpS1663 = _M0L11dst__offsetS64 + _M0L1iS66;
        int32_t _M0L6_2atmpS1665 = _M0L11src__offsetS65 + _M0L1iS66;
        int32_t _M0L6_2atmpS1664;
        int32_t _M0L6_2atmpS1666;
        if (
          _M0L6_2atmpS1665 < 0
          || _M0L6_2atmpS1665 >= Moonbit_array_length(_M0L3srcS63)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1664 = (int32_t)_M0L3srcS63[_M0L6_2atmpS1665];
        if (
          _M0L6_2atmpS1663 < 0
          || _M0L6_2atmpS1663 >= Moonbit_array_length(_M0L3dstS62)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS62[_M0L6_2atmpS1663] = _M0L6_2atmpS1664;
        _M0L6_2atmpS1666 = _M0L1iS66 + 1;
        _M0L1iS66 = _M0L6_2atmpS1666;
        continue;
      } else {
        moonbit_decref(_M0L3srcS63);
        moonbit_decref(_M0L3dstS62);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1671 = _M0L3lenS67 - 1;
    int32_t _M0L1iS69 = _M0L6_2atmpS1671;
    while (1) {
      if (_M0L1iS69 >= 0) {
        int32_t _M0L6_2atmpS1667 = _M0L11dst__offsetS64 + _M0L1iS69;
        int32_t _M0L6_2atmpS1669 = _M0L11src__offsetS65 + _M0L1iS69;
        int32_t _M0L6_2atmpS1668;
        int32_t _M0L6_2atmpS1670;
        if (
          _M0L6_2atmpS1669 < 0
          || _M0L6_2atmpS1669 >= Moonbit_array_length(_M0L3srcS63)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1668 = (int32_t)_M0L3srcS63[_M0L6_2atmpS1669];
        if (
          _M0L6_2atmpS1667 < 0
          || _M0L6_2atmpS1667 >= Moonbit_array_length(_M0L3dstS62)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS62[_M0L6_2atmpS1667] = _M0L6_2atmpS1668;
        _M0L6_2atmpS1670 = _M0L1iS69 - 1;
        _M0L1iS69 = _M0L6_2atmpS1670;
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
  int32_t _if__result_4635;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS71 == _M0L3srcS72) {
    _if__result_4635 = _M0L11dst__offsetS73 < _M0L11src__offsetS74;
  } else {
    _if__result_4635 = 0;
  }
  if (_if__result_4635) {
    int32_t _M0L1iS75 = 0;
    while (1) {
      if (_M0L1iS75 < _M0L3lenS76) {
        int32_t _M0L6_2atmpS1672 = _M0L11dst__offsetS73 + _M0L1iS75;
        int32_t _M0L6_2atmpS1674 = _M0L11src__offsetS74 + _M0L1iS75;
        struct _M0TPC16string10StringView _M0L6_2atmpS1673;
        struct _M0TPC16string10StringView _M0L6_2aoldS3914;
        int32_t _M0L6_2atmpS1675;
        if (
          _M0L6_2atmpS1674 < 0
          || _M0L6_2atmpS1674 >= Moonbit_array_length(_M0L3srcS72)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1673 = _M0L3srcS72[_M0L6_2atmpS1674];
        if (
          _M0L6_2atmpS1672 < 0
          || _M0L6_2atmpS1672 >= Moonbit_array_length(_M0L3dstS71)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3914 = _M0L3dstS71[_M0L6_2atmpS1672];
        moonbit_incref(_M0L6_2atmpS1673.$0);
        moonbit_decref(_M0L6_2aoldS3914.$0);
        _M0L3dstS71[_M0L6_2atmpS1672] = _M0L6_2atmpS1673;
        _M0L6_2atmpS1675 = _M0L1iS75 + 1;
        _M0L1iS75 = _M0L6_2atmpS1675;
        continue;
      } else {
        moonbit_decref(_M0L3srcS72);
        moonbit_decref(_M0L3dstS71);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1680 = _M0L3lenS76 - 1;
    int32_t _M0L1iS78 = _M0L6_2atmpS1680;
    while (1) {
      if (_M0L1iS78 >= 0) {
        int32_t _M0L6_2atmpS1676 = _M0L11dst__offsetS73 + _M0L1iS78;
        int32_t _M0L6_2atmpS1678 = _M0L11src__offsetS74 + _M0L1iS78;
        struct _M0TPC16string10StringView _M0L6_2atmpS1677;
        struct _M0TPC16string10StringView _M0L6_2aoldS3916;
        int32_t _M0L6_2atmpS1679;
        if (
          _M0L6_2atmpS1678 < 0
          || _M0L6_2atmpS1678 >= Moonbit_array_length(_M0L3srcS72)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1677 = _M0L3srcS72[_M0L6_2atmpS1678];
        if (
          _M0L6_2atmpS1676 < 0
          || _M0L6_2atmpS1676 >= Moonbit_array_length(_M0L3dstS71)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3916 = _M0L3dstS71[_M0L6_2atmpS1676];
        moonbit_incref(_M0L6_2atmpS1677.$0);
        moonbit_decref(_M0L6_2aoldS3916.$0);
        _M0L3dstS71[_M0L6_2atmpS1676] = _M0L6_2atmpS1677;
        _M0L6_2atmpS1679 = _M0L1iS78 - 1;
        _M0L1iS78 = _M0L6_2atmpS1679;
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

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t* _M0L4selfS21) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS21);
}

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(
  struct _M0TUsiE** _M0L4selfS22
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS22);
}

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L4selfS23
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS23);
}

int32_t _M0MPB18UninitializedArray6lengthGiE(int32_t* _M0L4selfS24) {
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
  uint32_t _M0L6_2atmpS1626;
  uint32_t _M0L6_2atmpS1625;
  uint32_t _M0L6_2atmpS1624;
  #line 465 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1626 = _M0L5inputS20 * 3266489917u;
  _M0L6_2atmpS1625 = _M0L3accS19 + _M0L6_2atmpS1626;
  #line 466 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1624 = _M0FPB4rotl(_M0L6_2atmpS1625, 17);
  return _M0L6_2atmpS1624 * 668265263u;
}

uint32_t _M0FPB4rotl(uint32_t _M0L1xS17, int32_t _M0L1rS18) {
  uint32_t _M0L6_2atmpS1621;
  int32_t _M0L6_2atmpS1623;
  uint32_t _M0L6_2atmpS1622;
  #line 475 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1621 = _M0L1xS17 << (_M0L1rS18 & 31);
  _M0L6_2atmpS1623 = 32 - _M0L1rS18;
  _M0L6_2atmpS1622 = _M0L1xS17 >> (_M0L6_2atmpS1623 & 31);
  return _M0L6_2atmpS1621 | _M0L6_2atmpS1622;
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
  _M0L10_2ax__5722S16.$0->$method_0(_M0L10_2ax__5722S16.$1, (moonbit_string_t)moonbit_string_literal_190.data);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB6Logger13write__objectGsE(_M0L10_2ax__5722S16, _M0L15_2a_2aarg__5723S15);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S16.$0->$method_0(_M0L10_2ax__5722S16.$1, (moonbit_string_t)moonbit_string_literal_191.data);
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

int32_t _M0FPC15abort5abortGuE(moonbit_string_t _M0L3msgS1) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS1);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
  return 0;
}

uint16_t* _M0FPC15abort5abortGAkE(moonbit_string_t _M0L3msgS2) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS2);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

moonbit_string_t* _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(
  moonbit_string_t _M0L3msgS3
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS3);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TUsiE** _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(
  moonbit_string_t _M0L3msgS4
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS4);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

moonbit_string_t _M0FPC15abort5abortGsE(moonbit_string_t _M0L3msgS5) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS5);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t _M0FPC15abort5abortGiE(moonbit_string_t _M0L3msgS6) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS6);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TPC16string10StringView _M0FPC15abort5abortGRPC16string10StringViewE(
  moonbit_string_t _M0L3msgS7
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS7);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(
  moonbit_string_t _M0L3msgS8
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS8);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t* _M0FPC15abort5abortGRPB18UninitializedArrayGiEE(
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

moonbit_string_t _M0FP15Error10to__string(void* _M0L4_2aeS1517) {
  switch (Moonbit_object_tag(_M0L4_2aeS1517)) {
    case 3: {
      return (moonbit_string_t)moonbit_string_literal_192.data;
      break;
    }
    
    case 5: {
      return (moonbit_string_t)moonbit_string_literal_193.data;
      break;
    }
    
    case 2: {
      return (moonbit_string_t)moonbit_string_literal_194.data;
      break;
    }
    
    case 1: {
      return (moonbit_string_t)moonbit_string_literal_195.data;
      break;
    }
    
    case 0: {
      return _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(_M0L4_2aeS1517);
      break;
    }
    default: {
      return (moonbit_string_t)moonbit_string_literal_196.data;
      break;
    }
  }
}

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1540,
  struct _M0TPB4Show _M0L8_2aparamS1539
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1538 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1540;
  _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(_M0L7_2aselfS1538, _M0L8_2aparamS1539);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1537,
  struct _M0TPB4Show _M0L8_2aparamS1536
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1535 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1537;
  _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(_M0L7_2aselfS1535, _M0L8_2aparamS1536);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1534,
  int32_t _M0L8_2aparamS1533
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1532 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1534;
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS1532, _M0L8_2aparamS1533);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1531,
  struct _M0TPC16string10StringView _M0L8_2aparamS1530
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1529 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1531;
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L7_2aselfS1529, _M0L8_2aparamS1530);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1528,
  moonbit_string_t _M0L8_2aparamS1525,
  int32_t _M0L8_2aparamS1526,
  int32_t _M0L8_2aparamS1527
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1524 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1528;
  _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L7_2aselfS1524, _M0L8_2aparamS1525, _M0L8_2aparamS1526, _M0L8_2aparamS1527);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1523,
  moonbit_string_t _M0L8_2aparamS1522
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1521 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1523;
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L7_2aselfS1521, _M0L8_2aparamS1522);
  return 0;
}

void moonbit_init() {
  moonbit_layout_table = moonbit_layout_table_data;
  uint16_t* _M0L1tS1230 = (uint16_t*)moonbit_make_string(128, 0);
  int32_t _M0L6_2atmpS1550 = 97;
  int32_t _M0L1cS1231 = _M0L6_2atmpS1550;
  moonbit_string_t* _M0L6_2atmpS1620;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1619;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1561;
  moonbit_string_t* _M0L6_2atmpS1618;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1617;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1562;
  moonbit_string_t* _M0L6_2atmpS1616;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1615;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1563;
  moonbit_string_t* _M0L6_2atmpS1614;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1613;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1564;
  moonbit_string_t* _M0L6_2atmpS1612;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1611;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1565;
  moonbit_string_t* _M0L6_2atmpS1610;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1609;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1566;
  moonbit_string_t* _M0L6_2atmpS1608;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1607;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1567;
  moonbit_string_t* _M0L6_2atmpS1606;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1605;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1568;
  moonbit_string_t* _M0L6_2atmpS1604;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1603;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1569;
  moonbit_string_t* _M0L6_2atmpS1602;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1601;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1570;
  moonbit_string_t* _M0L6_2atmpS1600;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1599;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1571;
  moonbit_string_t* _M0L6_2atmpS1598;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1597;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1572;
  moonbit_string_t* _M0L6_2atmpS1596;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1595;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1573;
  moonbit_string_t* _M0L6_2atmpS1594;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1593;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1574;
  moonbit_string_t* _M0L6_2atmpS1592;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1591;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1575;
  moonbit_string_t* _M0L6_2atmpS1590;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1589;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1576;
  moonbit_string_t* _M0L6_2atmpS1588;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1587;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1577;
  moonbit_string_t* _M0L6_2atmpS1586;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1585;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1578;
  moonbit_string_t* _M0L6_2atmpS1584;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1583;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1579;
  moonbit_string_t* _M0L6_2atmpS1582;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1581;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1580;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1439;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1560;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1559;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1558;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1553;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1440;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1557;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1556;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1555;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1554;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7_2abindS1438;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L6_2atmpS1552;
  struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE _M0L6_2atmpS1551;
  while (1) {
    int32_t _M0L6_2atmpS1546 = 122;
    if (_M0L1cS1231 <= _M0L6_2atmpS1546) {
      int32_t _M0L6_2atmpS1548 = _M0L1cS1231 - 32;
      int32_t _M0L6_2atmpS1547 = (uint16_t)_M0L6_2atmpS1548;
      int32_t _M0L6_2atmpS1549;
      if (
        _M0L1cS1231 < 0 || _M0L1cS1231 >= Moonbit_array_length(_M0L1tS1230)
      ) {
        #line 35 "/home/developer/Documents/1/genbank_io.mbt"
        moonbit_panic();
      }
      _M0L1tS1230[_M0L1cS1231] = _M0L6_2atmpS1547;
      _M0L6_2atmpS1549 = _M0L1cS1231 + 1;
      _M0L1cS1231 = _M0L6_2atmpS1549;
      continue;
    }
    break;
  }
  _M0FP26biolab8bio__seq16uppercase__table = _M0L1tS1230;
  _M0L6_2atmpS1620 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1620[0] = (moonbit_string_t)moonbit_string_literal_197.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__0_2eclo);
  _M0L8_2atupleS1619
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1619)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1619->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__0_2eclo;
  _M0L8_2atupleS1619->$1 = _M0L6_2atmpS1620;
  _M0L8_2atupleS1561
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1561)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1561->$0 = 0;
  _M0L8_2atupleS1561->$1 = _M0L8_2atupleS1619;
  _M0L6_2atmpS1618 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1618[0] = (moonbit_string_t)moonbit_string_literal_198.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__1_2eclo);
  _M0L8_2atupleS1617
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1617)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1617->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__1_2eclo;
  _M0L8_2atupleS1617->$1 = _M0L6_2atmpS1618;
  _M0L8_2atupleS1562
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1562)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1562->$0 = 1;
  _M0L8_2atupleS1562->$1 = _M0L8_2atupleS1617;
  _M0L6_2atmpS1616 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1616[0] = (moonbit_string_t)moonbit_string_literal_199.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__2_2eclo);
  _M0L8_2atupleS1615
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1615)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1615->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__2_2eclo;
  _M0L8_2atupleS1615->$1 = _M0L6_2atmpS1616;
  _M0L8_2atupleS1563
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1563)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1563->$0 = 2;
  _M0L8_2atupleS1563->$1 = _M0L8_2atupleS1615;
  _M0L6_2atmpS1614 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1614[0] = (moonbit_string_t)moonbit_string_literal_200.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__3_2eclo);
  _M0L8_2atupleS1613
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1613)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1613->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__3_2eclo;
  _M0L8_2atupleS1613->$1 = _M0L6_2atmpS1614;
  _M0L8_2atupleS1564
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1564)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1564->$0 = 3;
  _M0L8_2atupleS1564->$1 = _M0L8_2atupleS1613;
  _M0L6_2atmpS1612 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1612[0] = (moonbit_string_t)moonbit_string_literal_201.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__4_2eclo);
  _M0L8_2atupleS1611
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1611)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1611->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__4_2eclo;
  _M0L8_2atupleS1611->$1 = _M0L6_2atmpS1612;
  _M0L8_2atupleS1565
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1565)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1565->$0 = 4;
  _M0L8_2atupleS1565->$1 = _M0L8_2atupleS1611;
  _M0L6_2atmpS1610 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1610[0] = (moonbit_string_t)moonbit_string_literal_202.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__5_2eclo);
  _M0L8_2atupleS1609
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1609)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1609->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__5_2eclo;
  _M0L8_2atupleS1609->$1 = _M0L6_2atmpS1610;
  _M0L8_2atupleS1566
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1566)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1566->$0 = 5;
  _M0L8_2atupleS1566->$1 = _M0L8_2atupleS1609;
  _M0L6_2atmpS1608 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1608[0] = (moonbit_string_t)moonbit_string_literal_203.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__6_2eclo);
  _M0L8_2atupleS1607
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1607)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1607->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__6_2eclo;
  _M0L8_2atupleS1607->$1 = _M0L6_2atmpS1608;
  _M0L8_2atupleS1567
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1567)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1567->$0 = 6;
  _M0L8_2atupleS1567->$1 = _M0L8_2atupleS1607;
  _M0L6_2atmpS1606 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1606[0] = (moonbit_string_t)moonbit_string_literal_204.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__7_2eclo);
  _M0L8_2atupleS1605
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1605)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1605->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__7_2eclo;
  _M0L8_2atupleS1605->$1 = _M0L6_2atmpS1606;
  _M0L8_2atupleS1568
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1568)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1568->$0 = 7;
  _M0L8_2atupleS1568->$1 = _M0L8_2atupleS1605;
  _M0L6_2atmpS1604 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1604[0] = (moonbit_string_t)moonbit_string_literal_205.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__8_2eclo);
  _M0L8_2atupleS1603
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1603)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1603->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__8_2eclo;
  _M0L8_2atupleS1603->$1 = _M0L6_2atmpS1604;
  _M0L8_2atupleS1569
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1569)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1569->$0 = 8;
  _M0L8_2atupleS1569->$1 = _M0L8_2atupleS1603;
  _M0L6_2atmpS1602 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1602[0] = (moonbit_string_t)moonbit_string_literal_206.data;
  moonbit_incref(_M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__9_2eclo);
  _M0L8_2atupleS1601
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1601)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1601->$0
  = _M0FP26biolab8bio__seq51____test__736571696f5f7762746573742e6d6274__9_2eclo;
  _M0L8_2atupleS1601->$1 = _M0L6_2atmpS1602;
  _M0L8_2atupleS1570
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1570)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1570->$0 = 9;
  _M0L8_2atupleS1570->$1 = _M0L8_2atupleS1601;
  _M0L6_2atmpS1600 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1600[0] = (moonbit_string_t)moonbit_string_literal_207.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__10_2eclo);
  _M0L8_2atupleS1599
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1599)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1599->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__10_2eclo;
  _M0L8_2atupleS1599->$1 = _M0L6_2atmpS1600;
  _M0L8_2atupleS1571
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1571)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1571->$0 = 10;
  _M0L8_2atupleS1571->$1 = _M0L8_2atupleS1599;
  _M0L6_2atmpS1598 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1598[0] = (moonbit_string_t)moonbit_string_literal_208.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__11_2eclo);
  _M0L8_2atupleS1597
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1597)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1597->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__11_2eclo;
  _M0L8_2atupleS1597->$1 = _M0L6_2atmpS1598;
  _M0L8_2atupleS1572
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1572)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1572->$0 = 11;
  _M0L8_2atupleS1572->$1 = _M0L8_2atupleS1597;
  _M0L6_2atmpS1596 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1596[0] = (moonbit_string_t)moonbit_string_literal_209.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__12_2eclo);
  _M0L8_2atupleS1595
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1595)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1595->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__12_2eclo;
  _M0L8_2atupleS1595->$1 = _M0L6_2atmpS1596;
  _M0L8_2atupleS1573
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1573)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1573->$0 = 12;
  _M0L8_2atupleS1573->$1 = _M0L8_2atupleS1595;
  _M0L6_2atmpS1594 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1594[0] = (moonbit_string_t)moonbit_string_literal_210.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__13_2eclo);
  _M0L8_2atupleS1593
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1593)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1593->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__13_2eclo;
  _M0L8_2atupleS1593->$1 = _M0L6_2atmpS1594;
  _M0L8_2atupleS1574
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1574)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1574->$0 = 13;
  _M0L8_2atupleS1574->$1 = _M0L8_2atupleS1593;
  _M0L6_2atmpS1592 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1592[0] = (moonbit_string_t)moonbit_string_literal_211.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__14_2eclo);
  _M0L8_2atupleS1591
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1591)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1591->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__14_2eclo;
  _M0L8_2atupleS1591->$1 = _M0L6_2atmpS1592;
  _M0L8_2atupleS1575
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1575)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1575->$0 = 14;
  _M0L8_2atupleS1575->$1 = _M0L8_2atupleS1591;
  _M0L6_2atmpS1590 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1590[0] = (moonbit_string_t)moonbit_string_literal_212.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__15_2eclo);
  _M0L8_2atupleS1589
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1589)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1589->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__15_2eclo;
  _M0L8_2atupleS1589->$1 = _M0L6_2atmpS1590;
  _M0L8_2atupleS1576
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1576)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1576->$0 = 15;
  _M0L8_2atupleS1576->$1 = _M0L8_2atupleS1589;
  _M0L6_2atmpS1588 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1588[0] = (moonbit_string_t)moonbit_string_literal_213.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__16_2eclo);
  _M0L8_2atupleS1587
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1587)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1587->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__16_2eclo;
  _M0L8_2atupleS1587->$1 = _M0L6_2atmpS1588;
  _M0L8_2atupleS1577
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1577)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1577->$0 = 16;
  _M0L8_2atupleS1577->$1 = _M0L8_2atupleS1587;
  _M0L6_2atmpS1586 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1586[0] = (moonbit_string_t)moonbit_string_literal_214.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__17_2eclo);
  _M0L8_2atupleS1585
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1585)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1585->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__17_2eclo;
  _M0L8_2atupleS1585->$1 = _M0L6_2atmpS1586;
  _M0L8_2atupleS1578
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1578)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1578->$0 = 17;
  _M0L8_2atupleS1578->$1 = _M0L8_2atupleS1585;
  _M0L6_2atmpS1584 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1584[0] = (moonbit_string_t)moonbit_string_literal_215.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__18_2eclo);
  _M0L8_2atupleS1583
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1583)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1583->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__18_2eclo;
  _M0L8_2atupleS1583->$1 = _M0L6_2atmpS1584;
  _M0L8_2atupleS1579
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1579)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1579->$0 = 18;
  _M0L8_2atupleS1579->$1 = _M0L8_2atupleS1583;
  _M0L6_2atmpS1582 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1582[0] = (moonbit_string_t)moonbit_string_literal_216.data;
  moonbit_incref(_M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__19_2eclo);
  _M0L8_2atupleS1581
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1581)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 162, 0);
  _M0L8_2atupleS1581->$0
  = _M0FP26biolab8bio__seq52____test__736571696f5f7762746573742e6d6274__19_2eclo;
  _M0L8_2atupleS1581->$1 = _M0L6_2atmpS1582;
  _M0L8_2atupleS1580
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1580)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 166, 0);
  _M0L8_2atupleS1580->$0 = 19;
  _M0L8_2atupleS1580->$1 = _M0L8_2atupleS1581;
  _M0L7_2abindS1439
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_make_ref_array_raw(20);
  _M0L7_2abindS1439[0] = _M0L8_2atupleS1561;
  _M0L7_2abindS1439[1] = _M0L8_2atupleS1562;
  _M0L7_2abindS1439[2] = _M0L8_2atupleS1563;
  _M0L7_2abindS1439[3] = _M0L8_2atupleS1564;
  _M0L7_2abindS1439[4] = _M0L8_2atupleS1565;
  _M0L7_2abindS1439[5] = _M0L8_2atupleS1566;
  _M0L7_2abindS1439[6] = _M0L8_2atupleS1567;
  _M0L7_2abindS1439[7] = _M0L8_2atupleS1568;
  _M0L7_2abindS1439[8] = _M0L8_2atupleS1569;
  _M0L7_2abindS1439[9] = _M0L8_2atupleS1570;
  _M0L7_2abindS1439[10] = _M0L8_2atupleS1571;
  _M0L7_2abindS1439[11] = _M0L8_2atupleS1572;
  _M0L7_2abindS1439[12] = _M0L8_2atupleS1573;
  _M0L7_2abindS1439[13] = _M0L8_2atupleS1574;
  _M0L7_2abindS1439[14] = _M0L8_2atupleS1575;
  _M0L7_2abindS1439[15] = _M0L8_2atupleS1576;
  _M0L7_2abindS1439[16] = _M0L8_2atupleS1577;
  _M0L7_2abindS1439[17] = _M0L8_2atupleS1578;
  _M0L7_2abindS1439[18] = _M0L8_2atupleS1579;
  _M0L7_2abindS1439[19] = _M0L8_2atupleS1580;
  _M0L6_2atmpS1560 = _M0L7_2abindS1439;
  _M0L6_2atmpS1559
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1560, .$1 = 0, .$2 = 20
  };
  #line 400 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS1558
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1559, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1559.$0);
  _M0L8_2atupleS1553
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1553)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 169, 0);
  _M0L8_2atupleS1553->$0 = (moonbit_string_t)moonbit_string_literal_217.data;
  _M0L8_2atupleS1553->$1 = _M0L6_2atmpS1558;
  _M0L7_2abindS1440
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1557 = _M0L7_2abindS1440;
  _M0L6_2atmpS1556
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1557, .$1 = 0, .$2 = 0
  };
  #line 422 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L6_2atmpS1555
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1556, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1556.$0);
  _M0L8_2atupleS1554
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1554)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 169, 0);
  _M0L8_2atupleS1554->$0 = (moonbit_string_t)moonbit_string_literal_218.data;
  _M0L8_2atupleS1554->$1 = _M0L6_2atmpS1555;
  _M0L7_2abindS1438
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE**)moonbit_make_ref_array_raw(2);
  _M0L7_2abindS1438[0] = _M0L8_2atupleS1553;
  _M0L7_2abindS1438[1] = _M0L8_2atupleS1554;
  _M0L6_2atmpS1552 = _M0L7_2abindS1438;
  _M0L6_2atmpS1551
  = (struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE){
    .$0 = _M0L6_2atmpS1552, .$1 = 0, .$2 = 2
  };
  #line 399 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0FP26biolab8bio__seq48moonbit__test__driver__internal__no__args__tests
  = _M0MPB3Map3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L6_2atmpS1551, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1551.$0);
}

int main(int argc, char** argv) {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** _M0L6_2atmpS1545;
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1511;
  struct _M0TPB5ArrayGUsiEE* _M0L7_2abindS1512;
  int32_t _M0L7_2abindS1513;
  int32_t _M0L2__S1514;
  moonbit_runtime_init(argc, argv);
  moonbit_init();
  _M0L6_2atmpS1545
  = (struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error**)moonbit_empty_ref_array;
  _M0L12async__testsS1511
  = (struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE));
  Moonbit_object_header(_M0L12async__testsS1511)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 173, 0);
  _M0L12async__testsS1511->$0 = _M0L6_2atmpS1545;
  _M0L12async__testsS1511->$1 = 0;
  #line 463 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0L7_2abindS1512
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__args();
  _M0L7_2abindS1513 = _M0L7_2abindS1512->$1;
  _M0L2__S1514 = 0;
  while (1) {
    if (_M0L2__S1514 < _M0L7_2abindS1513) {
      struct _M0TUsiE** _M0L3bufS1544 = _M0L7_2abindS1512->$0;
      struct _M0TUsiE* _M0L3argS1515 =
        (struct _M0TUsiE*)_M0L3bufS1544[_M0L2__S1514];
      moonbit_string_t _M0L6_2atmpS1541 = _M0L3argS1515->$0;
      int32_t _M0L6_2atmpS1542 = _M0L3argS1515->$1;
      int32_t _M0L6_2atmpS1543;
      moonbit_incref(_M0L6_2atmpS1541);
      #line 464 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
      _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__execute(_M0L12async__testsS1511, _M0L6_2atmpS1541, _M0L6_2atmpS1542);
      moonbit_decref(_M0L6_2atmpS1541);
      _M0L6_2atmpS1543 = _M0L2__S1514 + 1;
      _M0L2__S1514 = _M0L6_2atmpS1543;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1512);
    }
    break;
  }
  #line 466 "/home/developer/Documents/1/__generated_driver_for_whitebox_test.mbt"
  _M0IP016_24default__implP26biolab8bio__seq28MoonBit__Async__Test__Driver17run__async__testsGRP26biolab8bio__seq34MoonBit__Async__Test__Driver__ImplE(_M0L12async__testsS1511);
  moonbit_decref(_M0L12async__testsS1511);
  return 0;
}