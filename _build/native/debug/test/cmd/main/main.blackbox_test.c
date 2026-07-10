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
struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd20main__blackbox__test33MoonBitTestDriverInternalSkipTestE2Ok;

struct _M0R38String_3a_3aiter_2eanon__u2264__l256__;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd20main__blackbox__test33MoonBitTestDriverInternalSkipTestE3Err;

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure;

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError;

struct _M0TPB8MutLocalGORPC16string10StringViewE;

struct _M0TWRPC15error5ErrorEs;

struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__;

struct _M0TWssbEu;

struct _M0TUsiE;

struct _M0TPB13StringBuilder;

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq3SeqRP26biolab8bio__seq8SeqErrorE3Err;

struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE;

struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError;

struct _M0TPB5EntryGicE;

struct _M0TP26biolab8bio__seq3Seq;

struct _M0BTPB6Logger;

struct _M0TPB4IterGcE;

struct _M0TPB6Logger;

struct _M0TWcERPC16string10StringView;

struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__;

struct _M0TPB5ArrayGUsiEE;

struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq3SeqRP26biolab8bio__seq8SeqErrorE2Ok;

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE;

struct _M0DTPC15error5Error116biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError;

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError;

struct _M0TWRPC15error5ErrorEu;

struct _M0TURPC16string10StringViewRPB6LoggerE;

struct _M0TUicE;

struct _M0TPB4IterGRPC16string10StringViewE;

struct _M0TPB9ArrayViewGUicEE;

struct _M0TPB8MutLocalGiE;

struct _M0TWERPC16option6OptionGRPC16string10StringViewE;

struct _M0TPB3MapGicE;

struct _M0TPB4Show;

struct _M0TWEOc;

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__;

struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE;

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err;

struct _M0BTPB4Show;

struct _M0TP46biolab8bio__seq3cmd4main4Case;

struct _M0TWuEu;

struct _M0TPC16string10StringView;

struct _M0TPB8MutLocalGbE;

struct _M0KTPB6LoggerTPB13StringBuilder;

struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest;

struct _M0TPB5ArrayGsE;

struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd20main__blackbox__test33MoonBitTestDriverInternalSkipTestE2Ok {
  int32_t $0;
  
};

struct _M0R38String_3a_3aiter_2eanon__u2264__l256__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  moonbit_string_t $1;
  int32_t $2;
  
};

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd20main__blackbox__test33MoonBitTestDriverInternalSkipTestE3Err {
  void* $0;
  
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

struct _M0TWRPC15error5ErrorEs {
  moonbit_string_t(* code)(struct _M0TWRPC15error5ErrorEs*, void*);
  
};

struct _M0TPC16string10StringView {
  moonbit_string_t $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  int32_t $1;
  struct _M0TPC16string10StringView $2;
  
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

struct _M0TPB13StringBuilder {
  uint16_t* $0;
  int32_t $1;
  
};

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some {
  struct _M0TPC16string10StringView $0;
  
};

struct _M0DTPC16result6ResultGRP26biolab8bio__seq3SeqRP26biolab8bio__seq8SeqErrorE3Err {
  void* $0;
  
};

struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE {
  struct _M0TP26biolab8bio__seq3Seq** $0;
  int32_t $1;
  
};

struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError {
  moonbit_string_t $0;
  
};

struct _M0TPB5EntryGicE {
  int32_t $0;
  struct _M0TPB5EntryGicE* $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  int32_t $5;
  
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

struct _M0TPB4IterGcE {
  struct _M0TWEOc* $0;
  int64_t $1;
  
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

struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TPB8MutLocalGORPC16string10StringViewE* $0;
  struct _M0TPC16string10StringView $1;
  int32_t $2;
  
};

struct _M0TPB5ArrayGUsiEE {
  struct _M0TUsiE** $0;
  int32_t $1;
  
};

struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175 {
  int32_t(* code)(struct _M0TWEOc*);
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok {
  int32_t $0;
  
};

struct _M0DTPC16result6ResultGRP26biolab8bio__seq3SeqRP26biolab8bio__seq8SeqErrorE2Ok {
  struct _M0TP26biolab8bio__seq3Seq* $0;
  
};

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** $0;
  int32_t $1;
  
};

struct _M0DTPC15error5Error116biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError {
  moonbit_string_t $0;
  
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

struct _M0TUicE {
  int32_t $0;
  int32_t $1;
  
};

struct _M0TPB4IterGRPC16string10StringViewE {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* $0;
  int64_t $1;
  
};

struct _M0TPB9ArrayViewGUicEE {
  struct _M0TUicE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0TPB8MutLocalGiE {
  int32_t $0;
  
};

struct _M0TWERPC16option6OptionGRPC16string10StringViewE {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  
};

struct _M0TPB3MapGicE {
  struct _M0TPB5EntryGicE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGicE* $5;
  int32_t $6;
  
};

struct _M0TPB4Show {
  struct _M0BTPB4Show* $0;
  void* $1;
  
};

struct _M0TWEOc {
  int32_t(* code)(struct _M0TWEOc*);
  
};

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TWcERPC16string10StringView* $0;
  struct _M0TPB4IterGcE* $1;
  
};

struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE {
  struct _M0TP46biolab8bio__seq3cmd4main4Case** $0;
  int32_t $1;
  
};

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error {
  struct moonbit_result_0(* code)(
    struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error*,
    struct _M0TWuEu*,
    struct _M0TWRPC15error5ErrorEu*
  );
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err {
  void* $0;
  
};

struct _M0BTPB4Show {
  int32_t(* $method_0)(void*, struct _M0TPB6Logger);
  moonbit_string_t(* $method_1)(void*);
  
};

struct _M0TP46biolab8bio__seq3cmd4main4Case {
  moonbit_string_t $0;
  moonbit_string_t $1;
  
};

struct _M0TWuEu {
  int32_t(* code)(struct _M0TWuEu*, int32_t);
  
};

struct _M0TPB8MutLocalGbE {
  int32_t $0;
  
};

struct _M0KTPB6LoggerTPB13StringBuilder {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest {
  moonbit_string_t $0;
  
};

struct _M0TPB5ArrayGsE {
  moonbit_string_t* $0;
  int32_t $1;
  
};

struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180 {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  int32_t $0;
  moonbit_string_t $1;
  
};

struct moonbit_result_0 {
  int tag;
  union { int32_t ok; void* err;  } data;
  
};

struct moonbit_result_1 {
  int tag;
  union { struct _M0TP26biolab8bio__seq3Seq* ok; void* err;  } data;
  
};

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1187(
  struct _M0TWRPC15error5ErrorEs*,
  void*
);

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1180(
  struct _M0TWssbEu*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1175(
  struct _M0TWEOc*
);

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__args(
  
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1153(
  int32_t,
  moonbit_string_t,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1146(
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139(
  int32_t,
  moonbit_bytes_t
);

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1132(
  int32_t,
  moonbit_string_t
);

#define _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__get__cli__args__ffi moonbit_get_cli_args

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd20main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*
);

struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0FP46biolab8bio__seq3cmd4main10run__cases(
  
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main19translate__to__stop(
  struct _M0TP26biolab8bio__seq3Seq*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main19translate__stop__at(
  struct _M0TP26biolab8bio__seq3Seq*,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main18translate__default(
  struct _M0TP26biolab8bio__seq3Seq*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main14translate__cds(
  struct _M0TP26biolab8bio__seq3Seq*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main14invalid__codon();

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main19cds__internal__stop(
  struct _M0TP26biolab8bio__seq3Seq*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main15cds__bad__start();

int32_t _M0FP46biolab8bio__seq3cmd4main12record__bool(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP46biolab8bio__seq3cmd4main6record(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*,
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_1 _M0MP26biolab8bio__seq3Seq17translate_2einner(
  struct _M0TP26biolab8bio__seq3Seq*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq13lookup__codon(moonbit_string_t, int32_t);

int32_t _M0FP26biolab8bio__seq10codon__key(moonbit_string_t, int32_t);

int32_t _M0FP26biolab8bio__seq19is__stop__codon__at(
  moonbit_string_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq20is__start__codon__at(
  moonbit_string_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq18is__gap__codon__at(
  moonbit_string_t,
  int32_t,
  int32_t
);

moonbit_string_t _M0FP26biolab8bio__seq17codon__to__string(
  moonbit_string_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq23all__valid__letters__at(
  moonbit_string_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq13valid__letter(int32_t);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq16back__transcribe(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq10transcribe(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq24reverse__complement__rna(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq7reverse(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq15complement__rna(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq10complement(
  struct _M0TP26biolab8bio__seq3Seq*
);

moonbit_string_t _M0FP26biolab8bio__seq17complement__chars(
  moonbit_string_t,
  uint16_t*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq4join(
  struct _M0TP26biolab8bio__seq3Seq*,
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*
);

struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0MP26biolab8bio__seq3Seq5split(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq7replace(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq14rstrip_2einner(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq14lstrip_2einner(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq13strip_2einner(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq8contains(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq8endswith(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq10startswith(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq5rfind(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq4find(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq14count__overlap(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

int32_t _M0MP26biolab8bio__seq3Seq5count(
  struct _M0TP26biolab8bio__seq3Seq*,
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq5lower(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq5upper(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq5slice(
  struct _M0TP26biolab8bio__seq3Seq*,
  int32_t,
  int32_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3mul(
  struct _M0TP26biolab8bio__seq3Seq*,
  int32_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0IP26biolab8bio__seq3SeqPB3Add3add(
  struct _M0TP26biolab8bio__seq3Seq*,
  struct _M0TP26biolab8bio__seq3Seq*
);

int64_t _M0FP26biolab8bio__seq10find__from(
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq19reverse__complement(
  struct _M0TP26biolab8bio__seq3Seq*
);

moonbit_string_t _M0FP26biolab8bio__seq26reverse__complement__chars(
  moonbit_string_t,
  uint16_t*
);

int32_t _M0MP26biolab8bio__seq3Seq6length(struct _M0TP26biolab8bio__seq3Seq*);

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t
);

int32_t _M0IP26biolab8bio__seq3SeqPB2Eq5equal(
  struct _M0TP26biolab8bio__seq3Seq*,
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MPC15array5Array2atGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*,
  int32_t
);

moonbit_string_t _M0MPC15array5Array2atGsE(struct _M0TPB5ArrayGsE*, int32_t);

int32_t _M0FPB7printlnGsE(moonbit_string_t);

int32_t _M0IPC13int3IntPB4Hash4hash(int32_t);

int32_t _M0MPB3Map3getGicE(struct _M0TPB3MapGicE*, int32_t);

struct _M0TPB3MapGicE* _M0MPB3Map3MapGicE(
  struct _M0TPB9ArrayViewGUicEE,
  int64_t
);

int32_t _M0MPB3Map3setGicE(struct _M0TPB3MapGicE*, int32_t, int32_t);

int32_t _M0MPB3Map15set__with__hashGicE(
  struct _M0TPB3MapGicE*,
  int32_t,
  int32_t,
  int32_t
);

int32_t _M0MPB3Map4growGicE(struct _M0TPB3MapGicE*);

int32_t _M0MPB3Map20rehash__place__entryGicE(
  struct _M0TPB3MapGicE*,
  struct _M0TPB5EntryGicE*
);

int32_t _M0MPB3Map10push__awayGicE(
  struct _M0TPB3MapGicE*,
  int32_t,
  struct _M0TPB5EntryGicE*
);

int32_t _M0MPB3Map10set__entryGicE(
  struct _M0TPB3MapGicE*,
  struct _M0TPB5EntryGicE*,
  int32_t
);

int32_t _M0MPB3Map20add__entry__to__tailGicE(
  struct _M0TPB3MapGicE*,
  int32_t,
  struct _M0TPB5EntryGicE*
);

int32_t _M0MPC13int3Int3max(int32_t, int32_t);

int32_t _M0FPB21capacity__for__length(int32_t);

struct _M0TPB3MapGicE* _M0FPB8new__mapGicE(int32_t);

int32_t _M0MPC13int3Int20next__power__of__two(int32_t);

int32_t _M0FPB21calc__grow__threshold(int32_t);

struct _M0TPB5EntryGicE* _M0MPC16option6Option6unwrapGRPB5EntryGicEE(
  struct _M0TPB5EntryGicE*
);

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t);

moonbit_string_t _M0MPC16string6String9to__upper(moonbit_string_t);

int32_t _M0MPC16string6String9to__upperC2411l1791(struct _M0TWuEu*, int32_t);

int32_t _M0MPC14char4Char20is__ascii__lowercase(int32_t);

moonbit_string_t _M0MPC16string6String9to__lower(moonbit_string_t);

int32_t _M0MPC16string6String9to__lowerC2377l1736(struct _M0TWuEu*, int32_t);

int32_t _M0MPC14char4Char20is__ascii__uppercase(int32_t);

moonbit_string_t _M0MPC16string6String12replace__all(
  moonbit_string_t,
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string6String5split(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string10StringView5split(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

void* _M0MPC16string10StringView5splitC2293l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC2289l1145(
  struct _M0TWcERPC16string10StringView*,
  int32_t
);

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t);

moonbit_string_t _M0FPB16char__to__string(int32_t);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE*,
  struct _M0TWcERPC16string10StringView*
);

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC2282l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPB4IterGcE* _M0MPC16string6String4iter(moonbit_string_t);

int32_t _M0MPC16string6String4iterC2264l256(struct _M0TWEOc*);

struct _M0TPC16string10StringView _M0MPC16string6String12trim_2einner(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string10StringView12trim_2einner(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string6String17trim__end_2einner(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string10StringView17trim__end_2einner(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string6String19trim__start_2einner(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string10StringView19trim__start_2einner(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView14contains__char(
  struct _M0TPC16string10StringView,
  int32_t
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

int32_t _M0MPC15array5Array4pushGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*,
  struct _M0TP46biolab8bio__seq3cmd4main4Case*
);

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*,
  struct _M0TP26biolab8bio__seq3Seq*
);

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array7reallocGUsiEE(struct _M0TPB5ArrayGUsiEE*);

int32_t _M0MPC15array5Array7reallocGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*
);

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*
);

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*,
  int32_t
);

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*
);

struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0MPC15array5Array11new_2einnerGRP46biolab8bio__seq3cmd4main4CaseE(
  int32_t
);

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(int32_t);

struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq3SeqE(
  int32_t
);

int64_t _M0MPC16string6String9rev__find(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

int64_t _M0MPC16string10StringView9rev__find(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int64_t _M0FPB23brute__force__rev__find(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

int64_t _M0MPC16string6String8find__by(moonbit_string_t, struct _M0TWuEu*);

int64_t _M0MPC16string10StringView8find__by(
  struct _M0TPC16string10StringView,
  struct _M0TWuEu*
);

moonbit_string_t _M0MPC16string6String6repeat(moonbit_string_t, int32_t);

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(moonbit_string_t);

int64_t _M0MPC16string6String4find(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

int64_t _M0FPB33boyer__moore__horspool__rev__find(
  struct _M0TPC16string10StringView,
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

int64_t _M0MPC16string6String29offset__of__nth__char_2einner(
  moonbit_string_t,
  int32_t,
  int32_t,
  int64_t
);

int64_t _M0MPC16string6String30offset__of__nth__char__forward(
  moonbit_string_t,
  int32_t,
  int32_t,
  int32_t
);

int64_t _M0MPC16string6String31offset__of__nth__char__backward(
  moonbit_string_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t,
  int32_t,
  int64_t
);

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView4iterC1978l219(struct _M0TWEOc*);

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

struct _M0TP26biolab8bio__seq3Seq** _M0MPC15array5Array6bufferGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*
);

moonbit_string_t* _M0MPC15array5Array6bufferGsE(struct _M0TPB5ArrayGsE*);

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*
);

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0MPC15array5Array6bufferGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*
);

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView,
  int32_t,
  int64_t
);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*,
  int64_t
);

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(struct _M0TWEOc*, int64_t);

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

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE*
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

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView,
  int32_t,
  int64_t
);

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t);

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t);

int32_t _M0IPC14byte4BytePB3Sub3sub(int32_t, int32_t);

int32_t _M0IPC14byte4BytePB3Mod3mod(int32_t, int32_t);

int32_t _M0IPC14byte4BytePB3Div3div(int32_t, int32_t);

int32_t _M0IPC14byte4BytePB3Add3add(int32_t, int32_t);

int32_t _M0MPC16string6String16unsafe__char__at(moonbit_string_t, int32_t);

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

int32_t _M0MPC13int3Int16unsafe__to__char(int32_t);

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

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TP26biolab8bio__seq3Seq** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq**,
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

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TP26biolab8bio__seq3Seq** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq**,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case**,
  int32_t,
  struct _M0TP46biolab8bio__seq3cmd4main4Case**,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq**,
  int32_t,
  struct _M0TP26biolab8bio__seq3Seq**,
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP46biolab8bio__seq3cmd4main4CaseEE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case**,
  int32_t,
  struct _M0TP46biolab8bio__seq3cmd4main4Case**,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq3SeqEE(
  struct _M0TP26biolab8bio__seq3Seq**,
  int32_t,
  struct _M0TP26biolab8bio__seq3Seq**,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t*);

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(struct _M0TUsiE**);

int32_t _M0MPB18UninitializedArray6lengthGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case**
);

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq**
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

int32_t _M0FPC15abort5abortGiE(moonbit_string_t);

moonbit_string_t _M0FPC15abort5abortGsE(moonbit_string_t);

struct _M0TPC16string10StringView _M0FPC15abort5abortGRPC16string10StringViewE(
  moonbit_string_t
);

moonbit_string_t* _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(
  moonbit_string_t
);

struct _M0TUsiE** _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(
  moonbit_string_t
);

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0FPC15abort5abortGRPB18UninitializedArrayGRP46biolab8bio__seq3cmd4main4CaseEE(
  moonbit_string_t
);

int64_t _M0FPC15abort5abortGOiE(moonbit_string_t);

struct _M0TP26biolab8bio__seq3Seq** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq3SeqEE(
  moonbit_string_t
);

int32_t _M0FP017____moonbit__initN3regS1004(
  struct _M0TPB3MapGicE*,
  int32_t,
  int32_t,
  int32_t,
  int32_t
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

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_50 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 99, 111, 
    117, 110, 116, 95, 65, 84, 65, 95, 65, 84, 65, 84, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[1]; 
} const moonbit_string_literal_0 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 0, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_84 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 115, 112, 
    108, 105, 116, 95, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_53 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 102, 105, 
    110, 100, 95, 71, 67, 67, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_37 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 99, 100, 115, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_87 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 101, 113, 
    95, 115, 97, 109, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[53]; 
} const moonbit_string_literal_133 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 52, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_48 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 65, 84, 
    65, 84, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_114 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 110, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_98 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 39, 32, 
    105, 115, 32, 110, 111, 116, 32, 97, 32, 115, 116, 97, 114, 116, 
    32, 99, 111, 100, 111, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_75 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 115, 116, 
    114, 105, 112, 95, 100, 97, 115, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[104]; 
} const moonbit_string_literal_131 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 103, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 109, 97, 105, 110, 95, 98, 108, 97, 99, 107, 98, 111, 120, 
    95, 116, 101, 115, 116, 46, 77, 111, 111, 110, 66, 105, 116, 84, 
    101, 115, 116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 
    110, 97, 108, 74, 115, 69, 114, 114, 111, 114, 46, 77, 111, 111, 
    110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 114, 
    73, 110, 116, 101, 114, 110, 97, 108, 74, 115, 69, 114, 114, 111, 
    114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_112 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_72 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 115, 108, 
    105, 99, 101, 95, 48, 95, 51, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_61 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 101, 110, 
    100, 115, 119, 105, 116, 104, 95, 84, 65, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[28]; 
} const moonbit_string_literal_100 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 27, 32, 105, 
    115, 32, 110, 111, 116, 32, 97, 32, 109, 117, 108, 116, 105, 112, 
    108, 101, 32, 111, 102, 32, 116, 104, 114, 101, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_31 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 112, 97, 114, 116, 105, 97, 
    108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_105 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 67, 111, 
    100, 111, 110, 32, 39, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_128 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 41, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_94 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 65, 65, 
    65, 67, 67, 67, 84, 65, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_126 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 46, 108, 101, 110, 103, 116, 104, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_27 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 100, 101, 102, 97, 117, 108, 
    116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_13 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 99, 111, 
    109, 112, 108, 101, 109, 101, 110, 116, 95, 67, 71, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_123 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_77 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 114, 115, 
    116, 114, 105, 112, 95, 100, 97, 115, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_76 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 108, 115, 
    116, 114, 105, 112, 95, 100, 97, 115, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_30 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 65, 84, 
    71, 71, 67, 67, 65, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_59 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 115, 116, 
    97, 114, 116, 115, 119, 105, 116, 104, 95, 65, 84, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_21 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 114, 101, 
    118, 101, 114, 115, 101, 95, 99, 111, 109, 112, 108, 101, 109, 101, 
    110, 116, 95, 114, 110, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_82 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 65, 67, 
    71, 84, 65, 67, 71, 84, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_64 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 99, 111, 
    110, 116, 97, 105, 110, 115, 95, 65, 65, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_56 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 102, 105, 
    110, 100, 95, 88, 89, 90, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_15 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 99, 111, 
    109, 112, 108, 101, 109, 101, 110, 116, 95, 67, 71, 65, 85, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_117 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 116, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_115 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 114, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_49 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 65, 84, 65, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_93 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 84, 114, 
    97, 110, 115, 108, 97, 116, 105, 111, 110, 69, 114, 114, 111, 114, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_74 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 45, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_22 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 65, 67, 
    71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[106]; 
} const moonbit_string_literal_129 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 105, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 109, 97, 105, 110, 95, 98, 108, 97, 99, 107, 98, 111, 120, 
    95, 116, 101, 115, 116, 46, 77, 111, 111, 110, 66, 105, 116, 84, 
    101, 115, 116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 
    110, 97, 108, 83, 107, 105, 112, 84, 101, 115, 116, 46, 77, 111, 
    111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 
    114, 73, 110, 116, 101, 114, 110, 97, 108, 83, 107, 105, 112, 84, 
    101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_20 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 114, 101, 
    118, 101, 114, 115, 101, 95, 99, 111, 109, 112, 108, 101, 109, 101, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_91 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 84, 65, 63, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_101 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 70, 105, 
    110, 97, 108, 32, 99, 111, 100, 111, 110, 32, 39, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_65 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 99, 111, 
    110, 116, 97, 105, 110, 115, 95, 88, 89, 90, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[40]; 
} const moonbit_string_literal_9 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 39, 65, 84, 
    71, 71, 67, 67, 65, 84, 84, 71, 84, 65, 65, 84, 71, 71, 71, 67, 67, 
    71, 67, 84, 71, 65, 65, 65, 71, 71, 71, 84, 71, 67, 67, 67, 71, 65, 
    84, 65, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_55 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 88, 89, 90, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_28 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 116, 111, 95, 115, 116, 111, 
    112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_70 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 65, 84, 
    71, 71, 67, 67, 65, 84, 84, 71, 84, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_62 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 71, 84, 71, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_79 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 84, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_54 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 114, 102, 
    105, 110, 100, 95, 71, 67, 67, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_10 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 116, 111, 
    95, 115, 116, 114, 105, 110, 103, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_132 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 83, 101, 
    113, 69, 114, 114, 111, 114, 46, 83, 101, 113, 69, 114, 114, 111, 
    114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[40]; 
} const moonbit_string_literal_25 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 39, 65, 85, 
    71, 71, 67, 67, 65, 85, 85, 71, 85, 65, 65, 85, 71, 71, 71, 67, 67, 
    71, 67, 85, 71, 65, 65, 65, 71, 71, 71, 85, 71, 67, 67, 67, 71, 65, 
    85, 65, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_12 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 67, 71, 65, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[20]; 
} const moonbit_string_literal_109 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 73, 110, 
    118, 97, 108, 105, 100, 32, 115, 116, 97, 114, 116, 32, 105, 110, 
    100, 101, 120, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_35 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 78, 78, 78, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_57 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 114, 102, 
    105, 110, 100, 95, 88, 89, 90, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_51 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 99, 111, 
    117, 110, 116, 95, 111, 118, 101, 114, 108, 97, 112, 95, 65, 84, 
    65, 95, 65, 84, 65, 84, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_47 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 99, 111, 
    117, 110, 116, 95, 111, 118, 101, 114, 108, 97, 112, 95, 65, 65, 
    65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_69 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 109, 117, 
    108, 51, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_42 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 117, 112, 
    112, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_41 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 97, 99, 
    103, 116, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_66 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 65, 67, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_32 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 84, 65, 78, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_46 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 99, 111, 
    117, 110, 116, 95, 65, 65, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_71 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 115, 108, 
    105, 99, 101, 95, 49, 95, 52, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_1 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 115, 107, 
    105, 112, 112, 101, 100, 32, 116, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_96 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 70, 97, 
    108, 115, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_18 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 99, 111, 
    109, 112, 108, 101, 109, 101, 110, 116, 95, 108, 111, 119, 101, 114, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_119 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 117, 114, 114, 111, 103, 97, 116, 
    101, 32, 112, 97, 105, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_83 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 115, 112, 
    108, 105, 116, 95, 99, 111, 117, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_43 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 108, 111, 
    119, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[31]; 
} const moonbit_string_literal_111 =
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

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_92 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 78, 79, 
    95, 69, 82, 82, 79, 82, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_67 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 71, 84, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[28]; 
} const moonbit_string_literal_39 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 27, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 99, 100, 115, 95, 105, 110, 
    116, 101, 114, 110, 97, 108, 95, 115, 116, 111, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_14 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 67, 71, 
    65, 85, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_125 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 44, 32, 
    108, 101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_122 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 111, 
    117, 110, 100, 115, 32, 99, 104, 101, 99, 107, 32, 102, 97, 105, 
    108, 101, 100, 58, 32, 97, 108, 108, 111, 99, 97, 116, 101, 95, 108, 
    101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_38 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 99, 100, 115, 95, 98, 97, 100, 
    95, 115, 116, 97, 114, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[56]; 
} const moonbit_string_literal_120 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 55, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 116, 97, 114, 116, 32, 111, 114, 
    32, 101, 110, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 114, 
    32, 83, 116, 114, 105, 110, 103, 58, 58, 99, 111, 100, 101, 112, 
    111, 105, 110, 116, 95, 108, 101, 110, 103, 116, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_44 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 65, 84, 
    65, 84, 71, 65, 65, 65, 84, 84, 84, 71, 65, 65, 65, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[28]; 
} const moonbit_string_literal_103 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 27, 69, 120, 
    116, 114, 97, 32, 105, 110, 32, 102, 114, 97, 109, 101, 32, 115, 
    116, 111, 112, 32, 99, 111, 100, 111, 110, 32, 39, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_33 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 84, 65, 78, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_19 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 99, 111, 
    109, 112, 108, 101, 109, 101, 110, 116, 95, 114, 110, 97, 95, 108, 
    111, 119, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_97 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 70, 105, 
    114, 115, 116, 32, 99, 111, 100, 111, 110, 32, 39, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_89 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 101, 113, 
    95, 100, 105, 102, 102, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_63 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 115, 116, 
    97, 114, 116, 115, 119, 105, 116, 104, 95, 71, 84, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_134 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 65, 67, 
    71, 84, 85, 77, 82, 87, 83, 89, 75, 86, 72, 68, 66, 88, 78, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_73 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 45, 45, 
    45, 65, 67, 71, 84, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_113 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 48, 49, 
    50, 51, 52, 53, 54, 55, 56, 57, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 119, 120, 121, 122, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_102 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 39, 32, 
    105, 115, 32, 110, 111, 116, 32, 97, 32, 115, 116, 111, 112, 32, 
    99, 111, 100, 111, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[51]; 
} const moonbit_string_literal_130 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 50, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 73, 110, 115, 112, 101, 
    99, 116, 69, 114, 114, 111, 114, 46, 73, 110, 115, 112, 101, 99, 
    116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_85 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 84, 84, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_68 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 97, 100, 100, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_7 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 45, 45, 
    45, 45, 45, 32, 69, 78, 68, 32, 77, 79, 79, 78, 32, 84, 69, 83, 84, 
    32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_8 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 115, 116, 97, 114, 116, 34, 44, 34, 
    102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_86 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 106, 111, 
    105, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_4 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 34, 
    105, 110, 100, 101, 120, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_90 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 69, 82, 
    82, 79, 82, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_127 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 70, 97, 
    105, 108, 117, 114, 101, 40, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_24 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 116, 114, 
    97, 110, 115, 99, 114, 105, 98, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_108 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 114, 101, 
    112, 101, 97, 116, 32, 114, 101, 115, 117, 108, 116, 32, 116, 111, 
    111, 32, 108, 97, 114, 103, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_6 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 125, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_99 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 83, 101, 
    113, 117, 101, 110, 99, 101, 32, 108, 101, 110, 103, 116, 104, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_88 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 84, 84, 
    84, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_29 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 115, 116, 111, 112, 95, 97, 
    116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_23 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 114, 101, 
    118, 101, 114, 115, 101, 95, 65, 67, 71, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[35]; 
} const moonbit_string_literal_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 34, 45, 45, 
    45, 45, 45, 32, 66, 69, 71, 73, 78, 32, 77, 79, 79, 78, 32, 84, 69, 
    83, 84, 32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_124 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    100, 115, 116, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_121 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 105, 110, 
    118, 97, 108, 105, 100, 32, 99, 111, 100, 101, 32, 112, 111, 105, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_40 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 116, 114, 
    97, 110, 115, 108, 97, 116, 101, 95, 105, 110, 118, 97, 108, 105, 
    100, 95, 99, 111, 100, 111, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_34 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 78, 78, 78, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_5 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 44, 34, 
    109, 101, 115, 115, 97, 103, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_110 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 73, 110, 
    118, 97, 108, 105, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 
    114, 32, 86, 105, 101, 119, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_95 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 84, 114, 
    117, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_106 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 39, 32, 
    105, 115, 32, 105, 110, 118, 97, 108, 105, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_118 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 92, 117, 123, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[34]; 
} const moonbit_string_literal_36 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 33, 65, 84, 
    71, 71, 67, 67, 65, 84, 84, 71, 71, 67, 67, 84, 71, 65, 65, 65, 71, 
    71, 71, 84, 71, 67, 67, 67, 67, 67, 67, 71, 84, 65, 71, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_16 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 99, 111, 
    109, 112, 108, 101, 109, 101, 110, 116, 95, 114, 110, 97, 95, 67, 
    71, 65, 85, 84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_11 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 108, 101, 
    110, 103, 116, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_80 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 85, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_60 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 84, 65, 71, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_81 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 114, 101, 
    112, 108, 97, 99, 101, 95, 84, 85, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_45 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 65, 65, 65, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_26 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 98, 97, 
    99, 107, 95, 116, 114, 97, 110, 115, 99, 114, 105, 98, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_17 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 99, 103, 
    97, 117, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_116 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_104 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 39, 32, 
    102, 111, 117, 110, 100, 46, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_107 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 110, 101, 
    103, 97, 116, 105, 118, 101, 32, 114, 101, 112, 101, 97, 116, 32, 
    99, 111, 117, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_78 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 65, 84, 
    71, 71, 67, 67, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_58 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 65, 84, 71, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_52 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 71, 67, 67, 0};

struct moonbit_object const moonbit_constant_constructor_0 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0)
  };

struct { int32_t rc; uint32_t meta; struct _M0TWuEu data; 
} const _M0MPC16string6String9to__lowerC2377l1736$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0MPC16string6String9to__lowerC2377l1736
  };

struct { int32_t rc; uint32_t meta; struct _M0TWRPC15error5ErrorEs data; 
} const _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1187$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1187
  };

struct { int32_t rc; uint32_t meta; struct _M0TWuEu data; 
} const _M0MPC16string6String9to__upperC2411l1791$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0MPC16string6String9to__upperC2411l1791
  };

struct {
  int32_t rc;
  uint32_t meta;
  struct _M0TWcERPC16string10StringView data;
  
} const _M0MPC16string10StringView5splitC2289l1145$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0MPC16string10StringView5splitC2289l1145
  };

uint32_t const moonbit_layout_table_data[87] =
  {
    sizeof(struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175)
    / 4, 1,
    offsetof(struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175, $1)
    / 4,
    sizeof(struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180)
    / 4, 1,
    offsetof(struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180, $1)
    / 4,
    sizeof(struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest, $0)
    / 4, sizeof(struct _M0TPB5ArrayGUsiEE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGUsiEE, $0) / 4, sizeof(struct _M0TUsiE) / 4,
    1, offsetof(struct _M0TUsiE, $0) / 4, sizeof(struct _M0TPB5ArrayGsE) / 4,
    1, offsetof(struct _M0TPB5ArrayGsE, $0) / 4,
    sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE, $0) / 4,
    sizeof(struct _M0TP46biolab8bio__seq3cmd4main4Case) / 4, 2,
    offsetof(struct _M0TP46biolab8bio__seq3cmd4main4Case, $0) / 4,
    offsetof(struct _M0TP46biolab8bio__seq3cmd4main4Case, $1) / 4,
    sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError, $0)
    / 4, sizeof(struct _M0TP26biolab8bio__seq3Seq) / 4, 1,
    offsetof(struct _M0TP26biolab8bio__seq3Seq, $0) / 4,
    sizeof(struct _M0TPB5EntryGicE) / 4, 1,
    offsetof(struct _M0TPB5EntryGicE, $1) / 4,
    sizeof(struct _M0TPB3MapGicE) / 4, 2,
    offsetof(struct _M0TPB3MapGicE, $0) / 4,
    offsetof(struct _M0TPB3MapGicE, $5) / 4,
    sizeof(struct _M0TPC16string10StringView) / 4, 1,
    offsetof(struct _M0TPC16string10StringView, $0) / 4,
    sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some) / 4,
    1,
    (offsetof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some, $0)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE) / 4, 
    1, offsetof(struct _M0TPB8MutLocalGORPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__) / 4, 
    2,
    offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__, $0)
    / 4,
    (offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__, $1)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4,
    sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__)
    / 4, 2,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__, $0)
    / 4,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__, $1)
    / 4, sizeof(struct _M0TPB4IterGRPC16string10StringViewE) / 4, 1,
    offsetof(struct _M0TPB4IterGRPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R38String_3a_3aiter_2eanon__u2264__l256__) / 4, 
    2, offsetof(struct _M0R38String_3a_3aiter_2eanon__u2264__l256__, $0) / 4,
    offsetof(struct _M0R38String_3a_3aiter_2eanon__u2264__l256__, $1) / 4,
    sizeof(struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE) / 4, 
    1,
    offsetof(struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE, $0) / 4,
    sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__) / 4, 
    2,
    offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__, $0) / 4,
    (offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__, $2)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB4IterGcE) / 4, 1,
    offsetof(struct _M0TPB4IterGcE, $0) / 4,
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

int64_t _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GcE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GcE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GcE = 0ll;

int64_t _M0FPB28boyer__moore__horspool__findN6constrS9990 = 0ll;

int64_t _M0FPB18brute__force__findN6constrS9991 = 0ll;

uint16_t* _M0FP26biolab8bio__seq20dna__complement__arr;

uint16_t* _M0FP26biolab8bio__seq20rna__complement__arr;

struct _M0TPB3MapGicE* _M0FP26biolab8bio__seq21codon__forward__table;

uint8_t* _M0FP26biolab8bio__seq20valid__letter__table;

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1208,
  moonbit_string_t _M0L8filenameS1177,
  int32_t _M0L5indexS1179
) {
  struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175* _closure_3226;
  struct _M0TWEOc* _M0L13handle__startS1175;
  struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180* _closure_3227;
  struct _M0TWssbEu* _M0L14handle__resultS1180;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1187;
  void* _M0L11_2atry__errS1202;
  struct moonbit_result_0 _tmp_3229;
  int32_t _handle__error__result_3230;
  int32_t _M0L6_2atmpS3046;
  void* _M0L3errS1203;
  moonbit_string_t _M0L4nameS1205;
  struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest* _M0L36_2aMoonBitTestDriverInternalSkipTestS1206;
  moonbit_string_t _M0L8_2afieldS3058;
  int32_t _M0L6_2acntS3220;
  moonbit_string_t _M0L7_2anameS1207;
  #line 515 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  moonbit_incref(_M0L8filenameS1177);
  _closure_3226
  = (struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175*)moonbit_malloc(sizeof(struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175));
  Moonbit_object_header(_closure_3226)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 0, 0);
  _closure_3226->code
  = &_M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1175;
  _closure_3226->$0 = _M0L5indexS1179;
  _closure_3226->$1 = _M0L8filenameS1177;
  _M0L13handle__startS1175 = (struct _M0TWEOc*)_closure_3226;
  moonbit_incref(_M0L8filenameS1177);
  _closure_3227
  = (struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180*)moonbit_malloc(sizeof(struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180));
  Moonbit_object_header(_closure_3227)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 3, 0);
  _closure_3227->code
  = &_M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1180;
  _closure_3227->$0 = _M0L5indexS1179;
  _closure_3227->$1 = _M0L8filenameS1177;
  _M0L14handle__resultS1180 = (struct _M0TWssbEu*)_closure_3227;
  _M0L17error__to__stringS1187
  = (struct _M0TWRPC15error5ErrorEs*)&_M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1187$closure.data;
  #line 557 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _tmp_3229
  = _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(_M0L12async__testsS1208, _M0L8filenameS1177, _M0L5indexS1179, _M0L13handle__startS1175, _M0L14handle__resultS1180, _M0L17error__to__stringS1187);
  if (_tmp_3229.tag) {
    int32_t const _M0L5_2aokS3055 = _tmp_3229.data.ok;
    _handle__error__result_3230 = _M0L5_2aokS3055;
  } else {
    void* const _M0L6_2aerrS3056 = _tmp_3229.data.err;
    moonbit_decref(_M0L17error__to__stringS1187);
    moonbit_decref(_M0L13handle__startS1175);
    _M0L11_2atry__errS1202 = _M0L6_2aerrS3056;
    goto join_1201;
  }
  if (_handle__error__result_3230) {
    moonbit_decref(_M0L17error__to__stringS1187);
    moonbit_decref(_M0L13handle__startS1175);
    _M0L6_2atmpS3046 = 1;
  } else {
    struct moonbit_result_0 _tmp_3231;
    int32_t _handle__error__result_3232;
    #line 560 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
    _tmp_3231
    = _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(_M0L12async__testsS1208, _M0L8filenameS1177, _M0L5indexS1179, _M0L13handle__startS1175, _M0L14handle__resultS1180, _M0L17error__to__stringS1187);
    if (_tmp_3231.tag) {
      int32_t const _M0L5_2aokS3053 = _tmp_3231.data.ok;
      _handle__error__result_3232 = _M0L5_2aokS3053;
    } else {
      void* const _M0L6_2aerrS3054 = _tmp_3231.data.err;
      moonbit_decref(_M0L17error__to__stringS1187);
      moonbit_decref(_M0L13handle__startS1175);
      _M0L11_2atry__errS1202 = _M0L6_2aerrS3054;
      goto join_1201;
    }
    if (_handle__error__result_3232) {
      moonbit_decref(_M0L17error__to__stringS1187);
      moonbit_decref(_M0L13handle__startS1175);
      _M0L6_2atmpS3046 = 1;
    } else {
      struct moonbit_result_0 _tmp_3233;
      int32_t _handle__error__result_3234;
      #line 563 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _tmp_3233
      = _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(_M0L12async__testsS1208, _M0L8filenameS1177, _M0L5indexS1179, _M0L13handle__startS1175, _M0L14handle__resultS1180, _M0L17error__to__stringS1187);
      if (_tmp_3233.tag) {
        int32_t const _M0L5_2aokS3051 = _tmp_3233.data.ok;
        _handle__error__result_3234 = _M0L5_2aokS3051;
      } else {
        void* const _M0L6_2aerrS3052 = _tmp_3233.data.err;
        moonbit_decref(_M0L17error__to__stringS1187);
        moonbit_decref(_M0L13handle__startS1175);
        _M0L11_2atry__errS1202 = _M0L6_2aerrS3052;
        goto join_1201;
      }
      if (_handle__error__result_3234) {
        moonbit_decref(_M0L17error__to__stringS1187);
        moonbit_decref(_M0L13handle__startS1175);
        _M0L6_2atmpS3046 = 1;
      } else {
        struct moonbit_result_0 _tmp_3235;
        int32_t _handle__error__result_3236;
        #line 566 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        _tmp_3235
        = _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(_M0L12async__testsS1208, _M0L8filenameS1177, _M0L5indexS1179, _M0L13handle__startS1175, _M0L14handle__resultS1180, _M0L17error__to__stringS1187);
        if (_tmp_3235.tag) {
          int32_t const _M0L5_2aokS3049 = _tmp_3235.data.ok;
          _handle__error__result_3236 = _M0L5_2aokS3049;
        } else {
          void* const _M0L6_2aerrS3050 = _tmp_3235.data.err;
          moonbit_decref(_M0L17error__to__stringS1187);
          moonbit_decref(_M0L13handle__startS1175);
          _M0L11_2atry__errS1202 = _M0L6_2aerrS3050;
          goto join_1201;
        }
        if (_handle__error__result_3236) {
          moonbit_decref(_M0L17error__to__stringS1187);
          moonbit_decref(_M0L13handle__startS1175);
          _M0L6_2atmpS3046 = 1;
        } else {
          struct moonbit_result_0 _tmp_3237;
          #line 569 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
          _tmp_3237
          = _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(_M0L12async__testsS1208, _M0L8filenameS1177, _M0L5indexS1179, _M0L13handle__startS1175, _M0L14handle__resultS1180, _M0L17error__to__stringS1187);
          moonbit_decref(_M0L13handle__startS1175);
          moonbit_decref(_M0L17error__to__stringS1187);
          if (_tmp_3237.tag) {
            int32_t const _M0L5_2aokS3047 = _tmp_3237.data.ok;
            _M0L6_2atmpS3046 = _M0L5_2aokS3047;
          } else {
            void* const _M0L6_2aerrS3048 = _tmp_3237.data.err;
            _M0L11_2atry__errS1202 = _M0L6_2aerrS3048;
            goto join_1201;
          }
        }
      }
    }
  }
  if (!_M0L6_2atmpS3046) {
    void* _M0L118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3057 =
      (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
    Moonbit_object_header(_M0L118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3057)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 2);
    ((struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3057)->$0
    = (moonbit_string_t)moonbit_string_literal_0.data;
    _M0L11_2atry__errS1202
    = _M0L118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3057;
    goto join_1201;
  } else {
    moonbit_decref(_M0L14handle__resultS1180);
  }
  goto joinlet_3228;
  join_1201:;
  _M0L3errS1203 = _M0L11_2atry__errS1202;
  _M0L36_2aMoonBitTestDriverInternalSkipTestS1206
  = (struct _M0DTPC15error5Error118biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L3errS1203;
  _M0L8_2afieldS3058 = _M0L36_2aMoonBitTestDriverInternalSkipTestS1206->$0;
  _M0L6_2acntS3220
  = Moonbit_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1206));
  if (_M0L6_2acntS3220 > 1) {
    int32_t _M0L11_2anew__cntS3221 = _M0L6_2acntS3220 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1206), _M0L11_2anew__cntS3221);
    moonbit_incref(_M0L8_2afieldS3058);
  } else if (_M0L6_2acntS3220 == 1) {
    #line 576 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
    moonbit_free(_M0L36_2aMoonBitTestDriverInternalSkipTestS1206);
  }
  _M0L7_2anameS1207 = _M0L8_2afieldS3058;
  _M0L4nameS1205 = _M0L7_2anameS1207;
  goto join_1204;
  goto joinlet_3238;
  join_1204:;
  #line 577 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1180(_M0L14handle__resultS1180, _M0L4nameS1205, (moonbit_string_t)moonbit_string_literal_1.data, 1);
  moonbit_decref(_M0L14handle__resultS1180);
  moonbit_decref(_M0L4nameS1205);
  joinlet_3238:;
  joinlet_3228:;
  return 0;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1187(
  struct _M0TWRPC15error5ErrorEs* _M0L6_2aenvS3045,
  void* _M0L3errS1188
) {
  void* _M0L1eS1190;
  moonbit_string_t _M0L1eS1192;
  moonbit_string_t _result_3241;
  #line 546 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  switch (Moonbit_object_tag(_M0L3errS1188)) {
    case 0: {
      struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS1193 =
        (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L3errS1188;
      moonbit_string_t _M0L4_2aeS1194 = _M0L10_2aFailureS1193->$0;
      moonbit_incref(_M0L4_2aeS1194);
      _M0L1eS1192 = _M0L4_2aeS1194;
      goto join_1191;
      break;
    }
    
    case 3: {
      struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError* _M0L15_2aInspectErrorS1195 =
        (struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError*)_M0L3errS1188;
      moonbit_string_t _M0L4_2aeS1196 = _M0L15_2aInspectErrorS1195->$0;
      moonbit_incref(_M0L4_2aeS1196);
      _M0L1eS1192 = _M0L4_2aeS1196;
      goto join_1191;
      break;
    }
    
    case 4: {
      struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError* _M0L16_2aSnapshotErrorS1197 =
        (struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError*)_M0L3errS1188;
      moonbit_string_t _M0L4_2aeS1198 = _M0L16_2aSnapshotErrorS1197->$0;
      moonbit_incref(_M0L4_2aeS1198);
      _M0L1eS1192 = _M0L4_2aeS1198;
      goto join_1191;
      break;
    }
    
    case 5: {
      struct _M0DTPC15error5Error116biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError* _M0L35_2aMoonBitTestDriverInternalJsErrorS1199 =
        (struct _M0DTPC15error5Error116biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError*)_M0L3errS1188;
      moonbit_string_t _M0L4_2aeS1200 =
        _M0L35_2aMoonBitTestDriverInternalJsErrorS1199->$0;
      moonbit_incref(_M0L4_2aeS1200);
      _M0L1eS1192 = _M0L4_2aeS1200;
      goto join_1191;
      break;
    }
    default: {
      moonbit_incref(_M0L3errS1188);
      _M0L1eS1190 = _M0L3errS1188;
      goto join_1189;
      break;
    }
  }
  join_1191:;
  return _M0L1eS1192;
  join_1189:;
  #line 552 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3241 = _M0FP15Error10to__string(_M0L1eS1190);
  moonbit_decref(_M0L1eS1190);
  return _result_3241;
}

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1180(
  struct _M0TWssbEu* _M0L6_2aenvS3042,
  moonbit_string_t _M0L10__testnameS1181,
  moonbit_string_t _M0L7messageS1182,
  int32_t _M0L7skippedS1183
) {
  struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180* _M0L14_2acasted__envS3043;
  moonbit_string_t _M0L8filenameS1177;
  int32_t _M0L5indexS1179;
  int32_t _if__result_3242;
  moonbit_string_t _M0L10file__nameS1184;
  moonbit_string_t _M0L7messageS1185;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1186;
  moonbit_string_t _M0L6_2atmpS3044;
  #line 531 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L14_2acasted__envS3043
  = (struct _M0R120_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1180*)_M0L6_2aenvS3042;
  _M0L8filenameS1177 = _M0L14_2acasted__envS3043->$1;
  _M0L5indexS1179 = _M0L14_2acasted__envS3043->$0;
  if (!_M0L7skippedS1183) {
    _if__result_3242 = 1;
  } else {
    _if__result_3242 = 0;
  }
  if (_if__result_3242) {
    
  }
  #line 537 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L10file__nameS1184
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1177, 1);
  #line 538 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L7messageS1185
  = _M0MPC16string6String14escape_2einner(_M0L7messageS1182, 1);
  #line 539 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L18_2astring__builderS1186
  = _M0MPB13StringBuilder21StringBuilder_2einner(45);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1186, (moonbit_string_t)moonbit_string_literal_3.data);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1186, _M0L10file__nameS1184);
  moonbit_decref(_M0L10file__nameS1184);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1186, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1186, _M0L5indexS1179);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1186, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1186, _M0L7messageS1185);
  moonbit_decref(_M0L7messageS1185);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1186, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 541 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS3044
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1186);
  moonbit_decref(_M0L18_2astring__builderS1186);
  #line 540 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS3044);
  moonbit_decref(_M0L6_2atmpS3044);
  #line 543 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1175(
  struct _M0TWEOc* _M0L6_2aenvS3039
) {
  struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175* _M0L14_2acasted__envS3040;
  moonbit_string_t _M0L8filenameS1177;
  int32_t _M0L5indexS1179;
  moonbit_string_t _M0L10file__nameS1176;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1178;
  moonbit_string_t _M0L6_2atmpS3041;
  #line 522 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L14_2acasted__envS3040
  = (struct _M0R119_24biolab_2fbio__seq_2fcmd_2fmain__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1175*)_M0L6_2aenvS3039;
  _M0L8filenameS1177 = _M0L14_2acasted__envS3040->$1;
  _M0L5indexS1179 = _M0L14_2acasted__envS3040->$0;
  #line 523 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L10file__nameS1176
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1177, 1);
  #line 524 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L18_2astring__builderS1178
  = _M0MPB13StringBuilder21StringBuilder_2einner(33);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1178, (moonbit_string_t)moonbit_string_literal_8.data);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1178, _M0L10file__nameS1176);
  moonbit_decref(_M0L10file__nameS1176);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1178, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1178, _M0L5indexS1179);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1178, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 526 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS3041
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1178);
  moonbit_decref(_M0L18_2astring__builderS1178);
  #line 525 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS3041);
  moonbit_decref(_M0L6_2atmpS3041);
  #line 528 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__args(
  
) {
  int32_t _M0L45moonbit__test__driver__internal__parse__int__S1132;
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139;
  int32_t _M0L57moonbit__test__driver__internal__get__cli__args__internalS1146;
  int32_t _M0L51moonbit__test__driver__internal__split__mbt__stringS1153;
  struct _M0TUsiE** _M0L6_2atmpS3038;
  struct _M0TPB5ArrayGUsiEE* _M0L16file__and__indexS1160;
  struct _M0TPB5ArrayGsE* _M0L9cli__argsS1161;
  moonbit_string_t _M0L6_2atmpS3037;
  struct _M0TPB5ArrayGsE* _M0L10test__argsS1162;
  int32_t _M0L7_2abindS1163;
  int32_t _M0L2__S1164;
  #line 194 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L45moonbit__test__driver__internal__parse__int__S1132 = 0;
  _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139
  = 0;
  _M0L57moonbit__test__driver__internal__get__cli__args__internalS1146
  = _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139;
  _M0L51moonbit__test__driver__internal__split__mbt__stringS1153 = 0;
  _M0L6_2atmpS3038 = (struct _M0TUsiE**)moonbit_empty_ref_array;
  _M0L16file__and__indexS1160
  = (struct _M0TPB5ArrayGUsiEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGUsiEE));
  Moonbit_object_header(_M0L16file__and__indexS1160)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 9, 0);
  _M0L16file__and__indexS1160->$0 = _M0L6_2atmpS3038;
  _M0L16file__and__indexS1160->$1 = 0;
  #line 283 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L9cli__argsS1161
  = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1146(_M0L57moonbit__test__driver__internal__get__cli__args__internalS1146);
  #line 285 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS3037 = _M0MPC15array5Array2atGsE(_M0L9cli__argsS1161, 1);
  moonbit_decref(_M0L9cli__argsS1161);
  #line 284 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L10test__argsS1162
  = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1153(_M0L51moonbit__test__driver__internal__split__mbt__stringS1153, _M0L6_2atmpS3037, 47);
  moonbit_decref(_M0L6_2atmpS3037);
  _M0L7_2abindS1163 = _M0L10test__argsS1162->$1;
  _M0L2__S1164 = 0;
  while (1) {
    if (_M0L2__S1164 < _M0L7_2abindS1163) {
      moonbit_string_t* _M0L3bufS3036 = _M0L10test__argsS1162->$0;
      moonbit_string_t _M0L3argS1165 =
        (moonbit_string_t)_M0L3bufS3036[_M0L2__S1164];
      struct _M0TPB5ArrayGsE* _M0L16file__and__rangeS1166;
      moonbit_string_t _M0L4fileS1167;
      moonbit_string_t _M0L5rangeS1168;
      struct _M0TPB5ArrayGsE* _M0L15start__and__endS1169;
      moonbit_string_t _M0L6_2atmpS3034;
      int32_t _M0L5startS1170;
      moonbit_string_t _M0L6_2atmpS3033;
      int32_t _M0L3endS1171;
      int32_t _M0L1iS1172;
      int32_t _M0L6_2atmpS3035;
      moonbit_incref(_M0L3argS1165);
      #line 289 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L16file__and__rangeS1166
      = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1153(_M0L51moonbit__test__driver__internal__split__mbt__stringS1153, _M0L3argS1165, 58);
      moonbit_decref(_M0L3argS1165);
      #line 290 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L4fileS1167
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1166, 0);
      #line 291 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L5rangeS1168
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1166, 1);
      moonbit_decref(_M0L16file__and__rangeS1166);
      #line 292 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L15start__and__endS1169
      = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1153(_M0L51moonbit__test__driver__internal__split__mbt__stringS1153, _M0L5rangeS1168, 45);
      moonbit_decref(_M0L5rangeS1168);
      #line 295 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS3034
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1169, 0);
      #line 295 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L5startS1170
      = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1132(_M0L45moonbit__test__driver__internal__parse__int__S1132, _M0L6_2atmpS3034);
      moonbit_decref(_M0L6_2atmpS3034);
      #line 296 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS3033
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1169, 1);
      moonbit_decref(_M0L15start__and__endS1169);
      #line 296 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L3endS1171
      = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1132(_M0L45moonbit__test__driver__internal__parse__int__S1132, _M0L6_2atmpS3033);
      moonbit_decref(_M0L6_2atmpS3033);
      _M0L1iS1172 = _M0L5startS1170;
      while (1) {
        if (_M0L1iS1172 < _M0L3endS1171) {
          struct _M0TUsiE* _M0L8_2atupleS3031;
          int32_t _M0L6_2atmpS3032;
          moonbit_incref(_M0L4fileS1167);
          _M0L8_2atupleS3031
          = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
          Moonbit_object_header(_M0L8_2atupleS3031)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 12, 0);
          _M0L8_2atupleS3031->$0 = _M0L4fileS1167;
          _M0L8_2atupleS3031->$1 = _M0L1iS1172;
          #line 298 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
          _M0MPC15array5Array4pushGUsiEE(_M0L16file__and__indexS1160, _M0L8_2atupleS3031);
          moonbit_decref(_M0L8_2atupleS3031);
          _M0L6_2atmpS3032 = _M0L1iS1172 + 1;
          _M0L1iS1172 = _M0L6_2atmpS3032;
          continue;
        } else {
          moonbit_decref(_M0L4fileS1167);
        }
        break;
      }
      _M0L6_2atmpS3035 = _M0L2__S1164 + 1;
      _M0L2__S1164 = _M0L6_2atmpS3035;
      continue;
    } else {
      moonbit_decref(_M0L10test__argsS1162);
    }
    break;
  }
  return _M0L16file__and__indexS1160;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1153(
  int32_t _M0L6_2aenvS3012,
  moonbit_string_t _M0L1sS1154,
  int32_t _M0L3sepS1155
) {
  moonbit_string_t* _M0L6_2atmpS3030;
  struct _M0TPB5ArrayGsE* _M0L3resS1156;
  struct _M0TPB8MutLocalGiE* _M0L1iS1157;
  struct _M0TPB8MutLocalGiE* _M0L5startS1158;
  int32_t _M0L3valS3025;
  int32_t _M0L6_2atmpS3026;
  #line 262 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS3030 = (moonbit_string_t*)moonbit_empty_ref_array;
  _M0L3resS1156
  = (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
  Moonbit_object_header(_M0L3resS1156)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
  _M0L3resS1156->$0 = _M0L6_2atmpS3030;
  _M0L3resS1156->$1 = 0;
  _M0L1iS1157
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1157)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1157->$0 = 0;
  _M0L5startS1158
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS1158)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS1158->$0 = 0;
  while (1) {
    int32_t _M0L3valS3013 = _M0L1iS1157->$0;
    int32_t _M0L6_2atmpS3014 = Moonbit_array_length(_M0L1sS1154);
    if (_M0L3valS3013 < _M0L6_2atmpS3014) {
      int32_t _M0L3valS3017 = _M0L1iS1157->$0;
      int32_t _M0L6_2atmpS3016;
      int32_t _M0L6_2atmpS3015;
      int32_t _M0L3valS3024;
      int32_t _M0L6_2atmpS3023;
      if (
        _M0L3valS3017 < 0
        || _M0L3valS3017 >= Moonbit_array_length(_M0L1sS1154)
      ) {
        #line 270 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3016 = _M0L1sS1154[_M0L3valS3017];
      _M0L6_2atmpS3015 = _M0L6_2atmpS3016;
      if (_M0L6_2atmpS3015 == _M0L3sepS1155) {
        int32_t _M0L3valS3019 = _M0L5startS1158->$0;
        int32_t _M0L3valS3020 = _M0L1iS1157->$0;
        moonbit_string_t _M0L6_2atmpS3018;
        int32_t _M0L3valS3022;
        int32_t _M0L6_2atmpS3021;
        #line 271 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        _M0L6_2atmpS3018
        = _M0MPC16string6String17unsafe__substring(_M0L1sS1154, _M0L3valS3019, _M0L3valS3020);
        #line 271 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        _M0MPC15array5Array4pushGsE(_M0L3resS1156, _M0L6_2atmpS3018);
        moonbit_decref(_M0L6_2atmpS3018);
        _M0L3valS3022 = _M0L1iS1157->$0;
        _M0L6_2atmpS3021 = _M0L3valS3022 + 1;
        _M0L5startS1158->$0 = _M0L6_2atmpS3021;
      }
      _M0L3valS3024 = _M0L1iS1157->$0;
      _M0L6_2atmpS3023 = _M0L3valS3024 + 1;
      _M0L1iS1157->$0 = _M0L6_2atmpS3023;
      continue;
    } else {
      moonbit_decref(_M0L1iS1157);
    }
    break;
  }
  _M0L3valS3025 = _M0L5startS1158->$0;
  _M0L6_2atmpS3026 = Moonbit_array_length(_M0L1sS1154);
  if (_M0L3valS3025 < _M0L6_2atmpS3026) {
    int32_t _M0L3valS3028 = _M0L5startS1158->$0;
    int32_t _M0L6_2atmpS3029;
    moonbit_string_t _M0L6_2atmpS3027;
    moonbit_decref(_M0L5startS1158);
    _M0L6_2atmpS3029 = Moonbit_array_length(_M0L1sS1154);
    #line 277 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
    _M0L6_2atmpS3027
    = _M0MPC16string6String17unsafe__substring(_M0L1sS1154, _M0L3valS3028, _M0L6_2atmpS3029);
    #line 277 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
    _M0MPC15array5Array4pushGsE(_M0L3resS1156, _M0L6_2atmpS3027);
    moonbit_decref(_M0L6_2atmpS3027);
  } else {
    moonbit_decref(_M0L5startS1158);
  }
  return _M0L3resS1156;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1146(
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139
) {
  moonbit_bytes_t* _M0L3tmpS1147;
  int32_t _M0L6_2atmpS3011;
  struct _M0TPB5ArrayGsE* _M0L3resS1148;
  int32_t _M0L7_2abindS1149;
  int32_t _M0L7_2abindS1150;
  int32_t _M0L1iS1151;
  #line 251 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  #line 254 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L3tmpS1147
  = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__get__cli__args__ffi();
  _M0L6_2atmpS3011 = Moonbit_array_length(_M0L3tmpS1147);
  #line 255 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1148 = _M0MPC15array5Array11new_2einnerGsE(_M0L6_2atmpS3011);
  _M0L7_2abindS1149 = 0;
  _M0L7_2abindS1150 = Moonbit_array_length(_M0L3tmpS1147);
  _M0L1iS1151 = _M0L7_2abindS1149;
  while (1) {
    if (_M0L1iS1151 < _M0L7_2abindS1150) {
      moonbit_bytes_t _M0L6_2atmpS3009;
      moonbit_string_t _M0L6_2atmpS3008;
      int32_t _M0L6_2atmpS3010;
      if (
        _M0L1iS1151 < 0 || _M0L1iS1151 >= Moonbit_array_length(_M0L3tmpS1147)
      ) {
        #line 257 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3009 = (moonbit_bytes_t)_M0L3tmpS1147[_M0L1iS1151];
      moonbit_incref(_M0L6_2atmpS3009);
      #line 257 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS3008
      = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139(_M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139, _M0L6_2atmpS3009);
      moonbit_decref(_M0L6_2atmpS3009);
      #line 257 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0MPC15array5Array4pushGsE(_M0L3resS1148, _M0L6_2atmpS3008);
      moonbit_decref(_M0L6_2atmpS3008);
      _M0L6_2atmpS3010 = _M0L1iS1151 + 1;
      _M0L1iS1151 = _M0L6_2atmpS3010;
      continue;
    } else {
      moonbit_decref(_M0L3tmpS1147);
    }
    break;
  }
  return _M0L3resS1148;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1139(
  int32_t _M0L6_2aenvS2922,
  moonbit_bytes_t _M0L5bytesS1140
) {
  struct _M0TPB13StringBuilder* _M0L3resS1141;
  int32_t _M0L3lenS1142;
  struct _M0TPB8MutLocalGiE* _M0L1iS1143;
  moonbit_string_t _result_3248;
  #line 207 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  #line 210 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1141 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3lenS1142 = Moonbit_array_length(_M0L5bytesS1140);
  _M0L1iS1143
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1143)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1143->$0 = 0;
  while (1) {
    int32_t _M0L3valS2923 = _M0L1iS1143->$0;
    if (_M0L3valS2923 < _M0L3lenS1142) {
      int32_t _M0L3valS3007 = _M0L1iS1143->$0;
      int32_t _M0L6_2atmpS3006;
      int32_t _M0L6_2atmpS3005;
      struct _M0TPB8MutLocalGiE* _M0L1cS1144;
      int32_t _M0L3valS2924;
      if (
        _M0L3valS3007 < 0
        || _M0L3valS3007 >= Moonbit_array_length(_M0L5bytesS1140)
      ) {
        #line 214 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3006 = _M0L5bytesS1140[_M0L3valS3007];
      _M0L6_2atmpS3005 = (int32_t)_M0L6_2atmpS3006;
      _M0L1cS1144
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1cS1144)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1cS1144->$0 = _M0L6_2atmpS3005;
      _M0L3valS2924 = _M0L1cS1144->$0;
      if (_M0L3valS2924 < 128) {
        int32_t _M0L3valS2926 = _M0L1cS1144->$0;
        int32_t _M0L6_2atmpS2925;
        int32_t _M0L3valS2928;
        int32_t _M0L6_2atmpS2927;
        moonbit_decref(_M0L1cS1144);
        _M0L6_2atmpS2925 = _M0L3valS2926;
        #line 216 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1141, _M0L6_2atmpS2925);
        _M0L3valS2928 = _M0L1iS1143->$0;
        _M0L6_2atmpS2927 = _M0L3valS2928 + 1;
        _M0L1iS1143->$0 = _M0L6_2atmpS2927;
      } else {
        int32_t _M0L3valS2929 = _M0L1cS1144->$0;
        if (_M0L3valS2929 < 224) {
          int32_t _M0L3valS2931 = _M0L1iS1143->$0;
          int32_t _M0L6_2atmpS2930 = _M0L3valS2931 + 1;
          int32_t _M0L3valS2940;
          int32_t _M0L6_2atmpS2939;
          int32_t _M0L6_2atmpS2933;
          int32_t _M0L3valS2938;
          int32_t _M0L6_2atmpS2937;
          int32_t _M0L6_2atmpS2936;
          int32_t _M0L6_2atmpS2935;
          int32_t _M0L6_2atmpS2934;
          int32_t _M0L6_2atmpS2932;
          int32_t _M0L3valS2942;
          int32_t _M0L6_2atmpS2941;
          int32_t _M0L3valS2944;
          int32_t _M0L6_2atmpS2943;
          if (_M0L6_2atmpS2930 >= _M0L3lenS1142) {
            moonbit_decref(_M0L1cS1144);
            moonbit_decref(_M0L1iS1143);
            break;
          }
          _M0L3valS2940 = _M0L1cS1144->$0;
          _M0L6_2atmpS2939 = _M0L3valS2940 & 31;
          _M0L6_2atmpS2933 = _M0L6_2atmpS2939 << 6;
          _M0L3valS2938 = _M0L1iS1143->$0;
          _M0L6_2atmpS2937 = _M0L3valS2938 + 1;
          if (
            _M0L6_2atmpS2937 < 0
            || _M0L6_2atmpS2937 >= Moonbit_array_length(_M0L5bytesS1140)
          ) {
            #line 222 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2936 = _M0L5bytesS1140[_M0L6_2atmpS2937];
          _M0L6_2atmpS2935 = (int32_t)_M0L6_2atmpS2936;
          _M0L6_2atmpS2934 = _M0L6_2atmpS2935 & 63;
          _M0L6_2atmpS2932 = _M0L6_2atmpS2933 | _M0L6_2atmpS2934;
          _M0L1cS1144->$0 = _M0L6_2atmpS2932;
          _M0L3valS2942 = _M0L1cS1144->$0;
          moonbit_decref(_M0L1cS1144);
          _M0L6_2atmpS2941 = _M0L3valS2942;
          #line 223 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1141, _M0L6_2atmpS2941);
          _M0L3valS2944 = _M0L1iS1143->$0;
          _M0L6_2atmpS2943 = _M0L3valS2944 + 2;
          _M0L1iS1143->$0 = _M0L6_2atmpS2943;
        } else {
          int32_t _M0L3valS2945 = _M0L1cS1144->$0;
          if (_M0L3valS2945 < 240) {
            int32_t _M0L3valS2947 = _M0L1iS1143->$0;
            int32_t _M0L6_2atmpS2946 = _M0L3valS2947 + 2;
            int32_t _M0L3valS2963;
            int32_t _M0L6_2atmpS2962;
            int32_t _M0L6_2atmpS2955;
            int32_t _M0L3valS2961;
            int32_t _M0L6_2atmpS2960;
            int32_t _M0L6_2atmpS2959;
            int32_t _M0L6_2atmpS2958;
            int32_t _M0L6_2atmpS2957;
            int32_t _M0L6_2atmpS2956;
            int32_t _M0L6_2atmpS2949;
            int32_t _M0L3valS2954;
            int32_t _M0L6_2atmpS2953;
            int32_t _M0L6_2atmpS2952;
            int32_t _M0L6_2atmpS2951;
            int32_t _M0L6_2atmpS2950;
            int32_t _M0L6_2atmpS2948;
            int32_t _M0L3valS2965;
            int32_t _M0L6_2atmpS2964;
            int32_t _M0L3valS2967;
            int32_t _M0L6_2atmpS2966;
            if (_M0L6_2atmpS2946 >= _M0L3lenS1142) {
              moonbit_decref(_M0L1cS1144);
              moonbit_decref(_M0L1iS1143);
              break;
            }
            _M0L3valS2963 = _M0L1cS1144->$0;
            _M0L6_2atmpS2962 = _M0L3valS2963 & 15;
            _M0L6_2atmpS2955 = _M0L6_2atmpS2962 << 12;
            _M0L3valS2961 = _M0L1iS1143->$0;
            _M0L6_2atmpS2960 = _M0L3valS2961 + 1;
            if (
              _M0L6_2atmpS2960 < 0
              || _M0L6_2atmpS2960 >= Moonbit_array_length(_M0L5bytesS1140)
            ) {
              #line 230 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2959 = _M0L5bytesS1140[_M0L6_2atmpS2960];
            _M0L6_2atmpS2958 = (int32_t)_M0L6_2atmpS2959;
            _M0L6_2atmpS2957 = _M0L6_2atmpS2958 & 63;
            _M0L6_2atmpS2956 = _M0L6_2atmpS2957 << 6;
            _M0L6_2atmpS2949 = _M0L6_2atmpS2955 | _M0L6_2atmpS2956;
            _M0L3valS2954 = _M0L1iS1143->$0;
            _M0L6_2atmpS2953 = _M0L3valS2954 + 2;
            if (
              _M0L6_2atmpS2953 < 0
              || _M0L6_2atmpS2953 >= Moonbit_array_length(_M0L5bytesS1140)
            ) {
              #line 231 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2952 = _M0L5bytesS1140[_M0L6_2atmpS2953];
            _M0L6_2atmpS2951 = (int32_t)_M0L6_2atmpS2952;
            _M0L6_2atmpS2950 = _M0L6_2atmpS2951 & 63;
            _M0L6_2atmpS2948 = _M0L6_2atmpS2949 | _M0L6_2atmpS2950;
            _M0L1cS1144->$0 = _M0L6_2atmpS2948;
            _M0L3valS2965 = _M0L1cS1144->$0;
            moonbit_decref(_M0L1cS1144);
            _M0L6_2atmpS2964 = _M0L3valS2965;
            #line 232 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1141, _M0L6_2atmpS2964);
            _M0L3valS2967 = _M0L1iS1143->$0;
            _M0L6_2atmpS2966 = _M0L3valS2967 + 3;
            _M0L1iS1143->$0 = _M0L6_2atmpS2966;
          } else {
            int32_t _M0L3valS2969 = _M0L1iS1143->$0;
            int32_t _M0L6_2atmpS2968 = _M0L3valS2969 + 3;
            int32_t _M0L3valS2992;
            int32_t _M0L6_2atmpS2991;
            int32_t _M0L6_2atmpS2984;
            int32_t _M0L3valS2990;
            int32_t _M0L6_2atmpS2989;
            int32_t _M0L6_2atmpS2988;
            int32_t _M0L6_2atmpS2987;
            int32_t _M0L6_2atmpS2986;
            int32_t _M0L6_2atmpS2985;
            int32_t _M0L6_2atmpS2977;
            int32_t _M0L3valS2983;
            int32_t _M0L6_2atmpS2982;
            int32_t _M0L6_2atmpS2981;
            int32_t _M0L6_2atmpS2980;
            int32_t _M0L6_2atmpS2979;
            int32_t _M0L6_2atmpS2978;
            int32_t _M0L6_2atmpS2971;
            int32_t _M0L3valS2976;
            int32_t _M0L6_2atmpS2975;
            int32_t _M0L6_2atmpS2974;
            int32_t _M0L6_2atmpS2973;
            int32_t _M0L6_2atmpS2972;
            int32_t _M0L6_2atmpS2970;
            int32_t _M0L3valS2994;
            int32_t _M0L6_2atmpS2993;
            int32_t _M0L3valS2998;
            int32_t _M0L6_2atmpS2997;
            int32_t _M0L6_2atmpS2996;
            int32_t _M0L6_2atmpS2995;
            int32_t _M0L3valS3002;
            int32_t _M0L6_2atmpS3001;
            int32_t _M0L6_2atmpS3000;
            int32_t _M0L6_2atmpS2999;
            int32_t _M0L3valS3004;
            int32_t _M0L6_2atmpS3003;
            if (_M0L6_2atmpS2968 >= _M0L3lenS1142) {
              moonbit_decref(_M0L1cS1144);
              moonbit_decref(_M0L1iS1143);
              break;
            }
            _M0L3valS2992 = _M0L1cS1144->$0;
            _M0L6_2atmpS2991 = _M0L3valS2992 & 7;
            _M0L6_2atmpS2984 = _M0L6_2atmpS2991 << 18;
            _M0L3valS2990 = _M0L1iS1143->$0;
            _M0L6_2atmpS2989 = _M0L3valS2990 + 1;
            if (
              _M0L6_2atmpS2989 < 0
              || _M0L6_2atmpS2989 >= Moonbit_array_length(_M0L5bytesS1140)
            ) {
              #line 239 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2988 = _M0L5bytesS1140[_M0L6_2atmpS2989];
            _M0L6_2atmpS2987 = (int32_t)_M0L6_2atmpS2988;
            _M0L6_2atmpS2986 = _M0L6_2atmpS2987 & 63;
            _M0L6_2atmpS2985 = _M0L6_2atmpS2986 << 12;
            _M0L6_2atmpS2977 = _M0L6_2atmpS2984 | _M0L6_2atmpS2985;
            _M0L3valS2983 = _M0L1iS1143->$0;
            _M0L6_2atmpS2982 = _M0L3valS2983 + 2;
            if (
              _M0L6_2atmpS2982 < 0
              || _M0L6_2atmpS2982 >= Moonbit_array_length(_M0L5bytesS1140)
            ) {
              #line 240 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2981 = _M0L5bytesS1140[_M0L6_2atmpS2982];
            _M0L6_2atmpS2980 = (int32_t)_M0L6_2atmpS2981;
            _M0L6_2atmpS2979 = _M0L6_2atmpS2980 & 63;
            _M0L6_2atmpS2978 = _M0L6_2atmpS2979 << 6;
            _M0L6_2atmpS2971 = _M0L6_2atmpS2977 | _M0L6_2atmpS2978;
            _M0L3valS2976 = _M0L1iS1143->$0;
            _M0L6_2atmpS2975 = _M0L3valS2976 + 3;
            if (
              _M0L6_2atmpS2975 < 0
              || _M0L6_2atmpS2975 >= Moonbit_array_length(_M0L5bytesS1140)
            ) {
              #line 241 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2974 = _M0L5bytesS1140[_M0L6_2atmpS2975];
            _M0L6_2atmpS2973 = (int32_t)_M0L6_2atmpS2974;
            _M0L6_2atmpS2972 = _M0L6_2atmpS2973 & 63;
            _M0L6_2atmpS2970 = _M0L6_2atmpS2971 | _M0L6_2atmpS2972;
            _M0L1cS1144->$0 = _M0L6_2atmpS2970;
            _M0L3valS2994 = _M0L1cS1144->$0;
            _M0L6_2atmpS2993 = _M0L3valS2994 - 65536;
            _M0L1cS1144->$0 = _M0L6_2atmpS2993;
            _M0L3valS2998 = _M0L1cS1144->$0;
            _M0L6_2atmpS2997 = _M0L3valS2998 >> 10;
            _M0L6_2atmpS2996 = _M0L6_2atmpS2997 + 55296;
            _M0L6_2atmpS2995 = _M0L6_2atmpS2996;
            #line 243 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1141, _M0L6_2atmpS2995);
            _M0L3valS3002 = _M0L1cS1144->$0;
            moonbit_decref(_M0L1cS1144);
            _M0L6_2atmpS3001 = _M0L3valS3002 & 1023;
            _M0L6_2atmpS3000 = _M0L6_2atmpS3001 + 56320;
            _M0L6_2atmpS2999 = _M0L6_2atmpS3000;
            #line 244 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1141, _M0L6_2atmpS2999);
            _M0L3valS3004 = _M0L1iS1143->$0;
            _M0L6_2atmpS3003 = _M0L3valS3004 + 4;
            _M0L1iS1143->$0 = _M0L6_2atmpS3003;
          }
        }
      }
      continue;
    } else {
      moonbit_decref(_M0L1iS1143);
    }
    break;
  }
  #line 248 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3248 = _M0MPB13StringBuilder10to__string(_M0L3resS1141);
  moonbit_decref(_M0L3resS1141);
  return _result_3248;
}

int32_t _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1132(
  int32_t _M0L6_2aenvS2915,
  moonbit_string_t _M0L1sS1133
) {
  struct _M0TPB8MutLocalGiE* _M0L3resS1134;
  int32_t _M0L3lenS1135;
  int32_t _M0L7_2abindS1136;
  int32_t _M0L1iS1137;
  int32_t _result_3250;
  #line 198 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1134
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3resS1134)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3resS1134->$0 = 0;
  _M0L3lenS1135 = Moonbit_array_length(_M0L1sS1133);
  _M0L7_2abindS1136 = 0;
  _M0L1iS1137 = _M0L7_2abindS1136;
  while (1) {
    if (_M0L1iS1137 < _M0L3lenS1135) {
      int32_t _M0L3valS2920 = _M0L3resS1134->$0;
      int32_t _M0L6_2atmpS2917 = _M0L3valS2920 * 10;
      int32_t _M0L6_2atmpS2919;
      int32_t _M0L6_2atmpS2918;
      int32_t _M0L6_2atmpS2916;
      int32_t _M0L6_2atmpS2921;
      if (
        _M0L1iS1137 < 0 || _M0L1iS1137 >= Moonbit_array_length(_M0L1sS1133)
      ) {
        #line 202 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2919 = _M0L1sS1133[_M0L1iS1137];
      _M0L6_2atmpS2918 = _M0L6_2atmpS2919 - 48;
      _M0L6_2atmpS2916 = _M0L6_2atmpS2917 + _M0L6_2atmpS2918;
      _M0L3resS1134->$0 = _M0L6_2atmpS2916;
      _M0L6_2atmpS2921 = _M0L1iS1137 + 1;
      _M0L1iS1137 = _M0L6_2atmpS2921;
      continue;
    }
    break;
  }
  _result_3250 = _M0L3resS1134->$0;
  moonbit_decref(_M0L3resS1134);
  return _result_3250;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1102,
  moonbit_string_t _M0L12_2adiscard__S1103,
  int32_t _M0L12_2adiscard__S1104,
  struct _M0TWEOc* _M0L12_2adiscard__S1105,
  struct _M0TWssbEu* _M0L12_2adiscard__S1106,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1107
) {
  struct moonbit_result_0 _result_3251;
  #line 35 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3251.tag = 1;
  _result_3251.data.ok = 0;
  return _result_3251;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1108,
  moonbit_string_t _M0L12_2adiscard__S1109,
  int32_t _M0L12_2adiscard__S1110,
  struct _M0TWEOc* _M0L12_2adiscard__S1111,
  struct _M0TWssbEu* _M0L12_2adiscard__S1112,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1113
) {
  struct moonbit_result_0 _result_3252;
  #line 35 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3252.tag = 1;
  _result_3252.data.ok = 0;
  return _result_3252;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1114,
  moonbit_string_t _M0L12_2adiscard__S1115,
  int32_t _M0L12_2adiscard__S1116,
  struct _M0TWEOc* _M0L12_2adiscard__S1117,
  struct _M0TWssbEu* _M0L12_2adiscard__S1118,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1119
) {
  struct moonbit_result_0 _result_3253;
  #line 35 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3253.tag = 1;
  _result_3253.data.ok = 0;
  return _result_3253;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1120,
  moonbit_string_t _M0L12_2adiscard__S1121,
  int32_t _M0L12_2adiscard__S1122,
  struct _M0TWEOc* _M0L12_2adiscard__S1123,
  struct _M0TWssbEu* _M0L12_2adiscard__S1124,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1125
) {
  struct moonbit_result_0 _result_3254;
  #line 35 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3254.tag = 1;
  _result_3254.data.ok = 0;
  return _result_3254;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd20main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1126,
  moonbit_string_t _M0L12_2adiscard__S1127,
  int32_t _M0L12_2adiscard__S1128,
  struct _M0TWEOc* _M0L12_2adiscard__S1129,
  struct _M0TWssbEu* _M0L12_2adiscard__S1130,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1131
) {
  struct moonbit_result_0 _result_3255;
  #line 35 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _result_3255.tag = 1;
  _result_3255.data.ok = 0;
  return _result_3255;
}

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd20main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1101
) {
  #line 12 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  return 0;
}

struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0FP46biolab8bio__seq3cmd4main10run__cases(
  
) {
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L5casesS1086;
  struct _M0TP26biolab8bio__seq3Seq* _M0L1sS1087;
  moonbit_string_t _M0L6_2atmpS2802;
  int32_t _M0L6_2atmpS2804;
  moonbit_string_t _M0L6_2atmpS2803;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2807;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2806;
  moonbit_string_t _M0L6_2atmpS2805;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2810;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2809;
  moonbit_string_t _M0L6_2atmpS2808;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2813;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2812;
  moonbit_string_t _M0L6_2atmpS2811;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2816;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2815;
  moonbit_string_t _M0L6_2atmpS2814;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2819;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2818;
  moonbit_string_t _M0L6_2atmpS2817;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2821;
  moonbit_string_t _M0L6_2atmpS2820;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2823;
  moonbit_string_t _M0L6_2atmpS2822;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2826;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2825;
  moonbit_string_t _M0L6_2atmpS2824;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2828;
  moonbit_string_t _M0L6_2atmpS2827;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3rnaS1088;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2830;
  moonbit_string_t _M0L6_2atmpS2829;
  moonbit_string_t _M0L6_2atmpS2831;
  moonbit_string_t _M0L6_2atmpS2832;
  moonbit_string_t _M0L6_2atmpS2833;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2835;
  moonbit_string_t _M0L6_2atmpS2834;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2837;
  moonbit_string_t _M0L6_2atmpS2836;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2839;
  moonbit_string_t _M0L6_2atmpS2838;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3cdsS1089;
  moonbit_string_t _M0L6_2atmpS2840;
  moonbit_string_t _M0L6_2atmpS2841;
  moonbit_string_t _M0L6_2atmpS2842;
  moonbit_string_t _M0L6_2atmpS2843;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2846;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2845;
  moonbit_string_t _M0L6_2atmpS2844;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2849;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2848;
  moonbit_string_t _M0L6_2atmpS2847;
  struct _M0TP26biolab8bio__seq3Seq* _M0L1cS1090;
  int32_t _M0L6_2atmpS2851;
  moonbit_string_t _M0L6_2atmpS2850;
  int32_t _M0L6_2atmpS2853;
  moonbit_string_t _M0L6_2atmpS2852;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2856;
  int32_t _M0L6_2atmpS2855;
  moonbit_string_t _M0L6_2atmpS2854;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2859;
  int32_t _M0L6_2atmpS2858;
  moonbit_string_t _M0L6_2atmpS2857;
  int32_t _M0L6_2atmpS2861;
  moonbit_string_t _M0L6_2atmpS2860;
  int32_t _M0L6_2atmpS2863;
  moonbit_string_t _M0L6_2atmpS2862;
  int32_t _M0L6_2atmpS2865;
  moonbit_string_t _M0L6_2atmpS2864;
  int32_t _M0L6_2atmpS2867;
  moonbit_string_t _M0L6_2atmpS2866;
  int32_t _M0L6_2atmpS2868;
  int32_t _M0L6_2atmpS2869;
  int32_t _M0L6_2atmpS2870;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2872;
  int32_t _M0L6_2atmpS2871;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2874;
  int32_t _M0L6_2atmpS2873;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2877;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2878;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2876;
  moonbit_string_t _M0L6_2atmpS2875;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2881;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2880;
  moonbit_string_t _M0L6_2atmpS2879;
  struct _M0TP26biolab8bio__seq3Seq* _M0L2slS1091;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2883;
  moonbit_string_t _M0L6_2atmpS2882;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2885;
  moonbit_string_t _M0L6_2atmpS2884;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2888;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2887;
  moonbit_string_t _M0L6_2atmpS2886;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2891;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2890;
  moonbit_string_t _M0L6_2atmpS2889;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2894;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2893;
  moonbit_string_t _M0L6_2atmpS2892;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2897;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2896;
  moonbit_string_t _M0L6_2atmpS2895;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2914;
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L5partsS1092;
  int32_t _M0L6_2atmpS2899;
  moonbit_string_t _M0L6_2atmpS2898;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2901;
  moonbit_string_t _M0L6_2atmpS2900;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6spacerS1093;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2911;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2912;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2913;
  struct _M0TP26biolab8bio__seq3Seq** _M0L6_2atmpS2910;
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L11join__partsS1094;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2903;
  moonbit_string_t _M0L6_2atmpS2902;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2905;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2906;
  int32_t _M0L6_2atmpS2904;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2908;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2909;
  int32_t _M0L6_2atmpS2907;
  #line 29 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 30 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L5casesS1086
  = _M0MPC15array5Array11new_2einnerGRP46biolab8bio__seq3cmd4main4CaseE(0);
  #line 31 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L1sS1087
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_9.data);
  #line 34 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2802 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L1sS1087);
  #line 34 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_10.data, _M0L6_2atmpS2802);
  moonbit_decref(_M0L6_2atmpS2802);
  #line 35 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2804 = _M0MP26biolab8bio__seq3Seq6length(_M0L1sS1087);
  #line 35 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2803
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2804, 10);
  #line 35 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_11.data, _M0L6_2atmpS2803);
  moonbit_decref(_M0L6_2atmpS2803);
  #line 38 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2807
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_12.data);
  #line 38 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2806 = _M0MP26biolab8bio__seq3Seq10complement(_M0L6_2atmpS2807);
  moonbit_decref(_M0L6_2atmpS2807);
  #line 38 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2805 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2806);
  moonbit_decref(_M0L6_2atmpS2806);
  #line 38 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_13.data, _M0L6_2atmpS2805);
  moonbit_decref(_M0L6_2atmpS2805);
  #line 42 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2810
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_14.data);
  #line 42 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2809 = _M0MP26biolab8bio__seq3Seq10complement(_M0L6_2atmpS2810);
  moonbit_decref(_M0L6_2atmpS2810);
  #line 42 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2808 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2809);
  moonbit_decref(_M0L6_2atmpS2809);
  #line 39 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_15.data, _M0L6_2atmpS2808);
  moonbit_decref(_M0L6_2atmpS2808);
  #line 47 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2813
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_14.data);
  #line 47 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2812
  = _M0MP26biolab8bio__seq3Seq15complement__rna(_M0L6_2atmpS2813);
  moonbit_decref(_M0L6_2atmpS2813);
  #line 47 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2811 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2812);
  moonbit_decref(_M0L6_2atmpS2812);
  #line 44 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_16.data, _M0L6_2atmpS2811);
  moonbit_decref(_M0L6_2atmpS2811);
  #line 52 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2816
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_17.data);
  #line 52 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2815 = _M0MP26biolab8bio__seq3Seq10complement(_M0L6_2atmpS2816);
  moonbit_decref(_M0L6_2atmpS2816);
  #line 52 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2814 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2815);
  moonbit_decref(_M0L6_2atmpS2815);
  #line 49 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_18.data, _M0L6_2atmpS2814);
  moonbit_decref(_M0L6_2atmpS2814);
  #line 57 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2819
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_17.data);
  #line 57 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2818
  = _M0MP26biolab8bio__seq3Seq15complement__rna(_M0L6_2atmpS2819);
  moonbit_decref(_M0L6_2atmpS2819);
  #line 57 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2817 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2818);
  moonbit_decref(_M0L6_2atmpS2818);
  #line 54 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_19.data, _M0L6_2atmpS2817);
  moonbit_decref(_M0L6_2atmpS2817);
  #line 61 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2821
  = _M0MP26biolab8bio__seq3Seq19reverse__complement(_M0L1sS1087);
  #line 61 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2820 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2821);
  moonbit_decref(_M0L6_2atmpS2821);
  #line 61 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_20.data, _M0L6_2atmpS2820);
  moonbit_decref(_M0L6_2atmpS2820);
  #line 65 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2823
  = _M0MP26biolab8bio__seq3Seq24reverse__complement__rna(_M0L1sS1087);
  #line 65 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2822 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2823);
  moonbit_decref(_M0L6_2atmpS2823);
  #line 62 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_21.data, _M0L6_2atmpS2822);
  moonbit_decref(_M0L6_2atmpS2822);
  #line 67 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2826
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_22.data);
  #line 67 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2825 = _M0MP26biolab8bio__seq3Seq7reverse(_M0L6_2atmpS2826);
  moonbit_decref(_M0L6_2atmpS2826);
  #line 67 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2824 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2825);
  moonbit_decref(_M0L6_2atmpS2825);
  #line 67 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_23.data, _M0L6_2atmpS2824);
  moonbit_decref(_M0L6_2atmpS2824);
  #line 70 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2828 = _M0MP26biolab8bio__seq3Seq10transcribe(_M0L1sS1087);
  #line 70 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2827 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2828);
  moonbit_decref(_M0L6_2atmpS2828);
  #line 70 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_24.data, _M0L6_2atmpS2827);
  moonbit_decref(_M0L6_2atmpS2827);
  #line 71 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L3rnaS1088
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_25.data);
  #line 72 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2830
  = _M0MP26biolab8bio__seq3Seq16back__transcribe(_M0L3rnaS1088);
  moonbit_decref(_M0L3rnaS1088);
  #line 72 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2829 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2830);
  moonbit_decref(_M0L6_2atmpS2830);
  #line 72 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_26.data, _M0L6_2atmpS2829);
  moonbit_decref(_M0L6_2atmpS2829);
  #line 75 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2831
  = _M0FP46biolab8bio__seq3cmd4main18translate__default(_M0L1sS1087);
  #line 75 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_27.data, _M0L6_2atmpS2831);
  moonbit_decref(_M0L6_2atmpS2831);
  #line 76 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2832
  = _M0FP46biolab8bio__seq3cmd4main19translate__to__stop(_M0L1sS1087);
  #line 76 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_28.data, _M0L6_2atmpS2832);
  moonbit_decref(_M0L6_2atmpS2832);
  #line 77 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2833
  = _M0FP46biolab8bio__seq3cmd4main19translate__stop__at(_M0L1sS1087, 64);
  #line 77 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_29.data, _M0L6_2atmpS2833);
  moonbit_decref(_M0L6_2atmpS2833);
  #line 81 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2835
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_30.data);
  #line 81 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2834
  = _M0FP46biolab8bio__seq3cmd4main18translate__default(_M0L6_2atmpS2835);
  moonbit_decref(_M0L6_2atmpS2835);
  #line 78 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_31.data, _M0L6_2atmpS2834);
  moonbit_decref(_M0L6_2atmpS2834);
  #line 83 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2837
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_32.data);
  #line 83 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2836
  = _M0FP46biolab8bio__seq3cmd4main18translate__default(_M0L6_2atmpS2837);
  moonbit_decref(_M0L6_2atmpS2837);
  #line 83 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_33.data, _M0L6_2atmpS2836);
  moonbit_decref(_M0L6_2atmpS2836);
  #line 84 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2839
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_34.data);
  #line 84 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2838
  = _M0FP46biolab8bio__seq3cmd4main18translate__default(_M0L6_2atmpS2839);
  moonbit_decref(_M0L6_2atmpS2839);
  #line 84 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_35.data, _M0L6_2atmpS2838);
  moonbit_decref(_M0L6_2atmpS2838);
  #line 87 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L3cdsS1089
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_36.data);
  #line 88 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2840
  = _M0FP46biolab8bio__seq3cmd4main14translate__cds(_M0L3cdsS1089);
  moonbit_decref(_M0L3cdsS1089);
  #line 88 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_37.data, _M0L6_2atmpS2840);
  moonbit_decref(_M0L6_2atmpS2840);
  #line 92 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2841 = _M0FP46biolab8bio__seq3cmd4main15cds__bad__start();
  #line 92 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_38.data, _M0L6_2atmpS2841);
  moonbit_decref(_M0L6_2atmpS2841);
  #line 93 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2842
  = _M0FP46biolab8bio__seq3cmd4main19cds__internal__stop(_M0L1sS1087);
  #line 93 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_39.data, _M0L6_2atmpS2842);
  moonbit_decref(_M0L6_2atmpS2842);
  #line 94 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2843 = _M0FP46biolab8bio__seq3cmd4main14invalid__codon();
  #line 94 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_40.data, _M0L6_2atmpS2843);
  moonbit_decref(_M0L6_2atmpS2843);
  #line 97 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2846
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_41.data);
  #line 97 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2845 = _M0MP26biolab8bio__seq3Seq5upper(_M0L6_2atmpS2846);
  moonbit_decref(_M0L6_2atmpS2846);
  #line 97 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2844 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2845);
  moonbit_decref(_M0L6_2atmpS2845);
  #line 97 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_42.data, _M0L6_2atmpS2844);
  moonbit_decref(_M0L6_2atmpS2844);
  #line 98 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2849
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_41.data);
  #line 98 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2848 = _M0MP26biolab8bio__seq3Seq5lower(_M0L6_2atmpS2849);
  moonbit_decref(_M0L6_2atmpS2849);
  #line 98 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2847 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2848);
  moonbit_decref(_M0L6_2atmpS2848);
  #line 98 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_43.data, _M0L6_2atmpS2847);
  moonbit_decref(_M0L6_2atmpS2847);
  #line 101 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L1cS1090
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_44.data);
  #line 102 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2851
  = _M0MP26biolab8bio__seq3Seq5count(_M0L1cS1090, (moonbit_string_t)moonbit_string_literal_45.data);
  #line 102 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2850
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2851, 10);
  #line 102 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_46.data, _M0L6_2atmpS2850);
  moonbit_decref(_M0L6_2atmpS2850);
  #line 103 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2853
  = _M0MP26biolab8bio__seq3Seq14count__overlap(_M0L1cS1090, (moonbit_string_t)moonbit_string_literal_45.data);
  moonbit_decref(_M0L1cS1090);
  #line 103 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2852
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2853, 10);
  #line 103 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_47.data, _M0L6_2atmpS2852);
  moonbit_decref(_M0L6_2atmpS2852);
  #line 107 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2856
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_48.data);
  #line 107 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2855
  = _M0MP26biolab8bio__seq3Seq5count(_M0L6_2atmpS2856, (moonbit_string_t)moonbit_string_literal_49.data);
  moonbit_decref(_M0L6_2atmpS2856);
  #line 107 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2854
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2855, 10);
  #line 104 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_50.data, _M0L6_2atmpS2854);
  moonbit_decref(_M0L6_2atmpS2854);
  #line 112 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2859
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_48.data);
  #line 112 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2858
  = _M0MP26biolab8bio__seq3Seq14count__overlap(_M0L6_2atmpS2859, (moonbit_string_t)moonbit_string_literal_49.data);
  moonbit_decref(_M0L6_2atmpS2859);
  #line 112 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2857
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2858, 10);
  #line 109 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_51.data, _M0L6_2atmpS2857);
  moonbit_decref(_M0L6_2atmpS2857);
  #line 116 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2861
  = _M0MP26biolab8bio__seq3Seq4find(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_52.data);
  #line 116 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2860
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2861, 10);
  #line 116 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_53.data, _M0L6_2atmpS2860);
  moonbit_decref(_M0L6_2atmpS2860);
  #line 117 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2863
  = _M0MP26biolab8bio__seq3Seq5rfind(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_52.data);
  #line 117 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2862
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2863, 10);
  #line 117 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_54.data, _M0L6_2atmpS2862);
  moonbit_decref(_M0L6_2atmpS2862);
  #line 118 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2865
  = _M0MP26biolab8bio__seq3Seq4find(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_55.data);
  #line 118 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2864
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2865, 10);
  #line 118 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_56.data, _M0L6_2atmpS2864);
  moonbit_decref(_M0L6_2atmpS2864);
  #line 119 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2867
  = _M0MP26biolab8bio__seq3Seq5rfind(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_55.data);
  #line 119 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2866
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2867, 10);
  #line 119 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_57.data, _M0L6_2atmpS2866);
  moonbit_decref(_M0L6_2atmpS2866);
  #line 122 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2868
  = _M0MP26biolab8bio__seq3Seq10startswith(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_58.data);
  #line 122 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_59.data, _M0L6_2atmpS2868);
  #line 123 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2869
  = _M0MP26biolab8bio__seq3Seq8endswith(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_60.data);
  #line 123 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_61.data, _M0L6_2atmpS2869);
  #line 124 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2870
  = _M0MP26biolab8bio__seq3Seq10startswith(_M0L1sS1087, (moonbit_string_t)moonbit_string_literal_62.data);
  moonbit_decref(_M0L1sS1087);
  #line 124 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_63.data, _M0L6_2atmpS2870);
  #line 130 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2872
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_44.data);
  #line 130 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2871
  = _M0MP26biolab8bio__seq3Seq8contains(_M0L6_2atmpS2872, (moonbit_string_t)moonbit_string_literal_45.data);
  moonbit_decref(_M0L6_2atmpS2872);
  #line 127 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_64.data, _M0L6_2atmpS2871);
  #line 135 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2874
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_44.data);
  #line 135 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2873
  = _M0MP26biolab8bio__seq3Seq8contains(_M0L6_2atmpS2874, (moonbit_string_t)moonbit_string_literal_55.data);
  moonbit_decref(_M0L6_2atmpS2874);
  #line 132 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_65.data, _M0L6_2atmpS2873);
  #line 139 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2877
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_66.data);
  #line 139 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2878
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_67.data);
  #line 139 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2876
  = _M0IP26biolab8bio__seq3SeqPB3Add3add(_M0L6_2atmpS2877, _M0L6_2atmpS2878);
  moonbit_decref(_M0L6_2atmpS2877);
  moonbit_decref(_M0L6_2atmpS2878);
  #line 139 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2875 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2876);
  moonbit_decref(_M0L6_2atmpS2876);
  #line 139 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_68.data, _M0L6_2atmpS2875);
  moonbit_decref(_M0L6_2atmpS2875);
  #line 140 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2881
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_66.data);
  #line 140 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2880 = _M0MP26biolab8bio__seq3Seq3mul(_M0L6_2atmpS2881, 3);
  moonbit_decref(_M0L6_2atmpS2881);
  #line 140 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2879 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2880);
  moonbit_decref(_M0L6_2atmpS2880);
  #line 140 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_69.data, _M0L6_2atmpS2879);
  moonbit_decref(_M0L6_2atmpS2879);
  #line 143 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L2slS1091
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_70.data);
  #line 144 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2883 = _M0MP26biolab8bio__seq3Seq5slice(_M0L2slS1091, 1, 4);
  #line 144 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2882 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2883);
  moonbit_decref(_M0L6_2atmpS2883);
  #line 144 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_71.data, _M0L6_2atmpS2882);
  moonbit_decref(_M0L6_2atmpS2882);
  #line 145 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2885 = _M0MP26biolab8bio__seq3Seq5slice(_M0L2slS1091, 0, 3);
  moonbit_decref(_M0L2slS1091);
  #line 145 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2884 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2885);
  moonbit_decref(_M0L6_2atmpS2885);
  #line 145 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_72.data, _M0L6_2atmpS2884);
  moonbit_decref(_M0L6_2atmpS2884);
  #line 151 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2888
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_73.data);
  #line 151 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2887
  = _M0MP26biolab8bio__seq3Seq13strip_2einner(_M0L6_2atmpS2888, (moonbit_string_t)moonbit_string_literal_74.data);
  moonbit_decref(_M0L6_2atmpS2888);
  #line 151 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2886 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2887);
  moonbit_decref(_M0L6_2atmpS2887);
  #line 148 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_75.data, _M0L6_2atmpS2886);
  moonbit_decref(_M0L6_2atmpS2886);
  #line 156 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2891
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_73.data);
  #line 156 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2890
  = _M0MP26biolab8bio__seq3Seq14lstrip_2einner(_M0L6_2atmpS2891, (moonbit_string_t)moonbit_string_literal_74.data);
  moonbit_decref(_M0L6_2atmpS2891);
  #line 156 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2889 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2890);
  moonbit_decref(_M0L6_2atmpS2890);
  #line 153 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_76.data, _M0L6_2atmpS2889);
  moonbit_decref(_M0L6_2atmpS2889);
  #line 161 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2894
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_73.data);
  #line 161 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2893
  = _M0MP26biolab8bio__seq3Seq14rstrip_2einner(_M0L6_2atmpS2894, (moonbit_string_t)moonbit_string_literal_74.data);
  moonbit_decref(_M0L6_2atmpS2894);
  #line 161 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2892 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2893);
  moonbit_decref(_M0L6_2atmpS2893);
  #line 158 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_77.data, _M0L6_2atmpS2892);
  moonbit_decref(_M0L6_2atmpS2892);
  #line 168 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2897
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_78.data);
  #line 168 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2896
  = _M0MP26biolab8bio__seq3Seq7replace(_M0L6_2atmpS2897, (moonbit_string_t)moonbit_string_literal_79.data, (moonbit_string_t)moonbit_string_literal_80.data);
  moonbit_decref(_M0L6_2atmpS2897);
  #line 168 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2895 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2896);
  moonbit_decref(_M0L6_2atmpS2896);
  #line 165 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_81.data, _M0L6_2atmpS2895);
  moonbit_decref(_M0L6_2atmpS2895);
  #line 172 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2914
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_82.data);
  #line 172 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L5partsS1092
  = _M0MP26biolab8bio__seq3Seq5split(_M0L6_2atmpS2914, (moonbit_string_t)moonbit_string_literal_66.data);
  moonbit_decref(_M0L6_2atmpS2914);
  #line 173 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2899
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq3SeqE(_M0L5partsS1092);
  #line 173 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2898
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2899, 10);
  #line 173 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_83.data, _M0L6_2atmpS2898);
  moonbit_decref(_M0L6_2atmpS2898);
  #line 174 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2901
  = _M0MPC15array5Array2atGRP26biolab8bio__seq3SeqE(_M0L5partsS1092, 1);
  moonbit_decref(_M0L5partsS1092);
  #line 174 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2900 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2901);
  moonbit_decref(_M0L6_2atmpS2901);
  #line 174 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_84.data, _M0L6_2atmpS2900);
  moonbit_decref(_M0L6_2atmpS2900);
  #line 177 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6spacerS1093
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_74.data);
  #line 179 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2911
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_66.data);
  #line 180 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2912
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_67.data);
  #line 181 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2913
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_85.data);
  _M0L6_2atmpS2910
  = (struct _M0TP26biolab8bio__seq3Seq**)moonbit_make_ref_array_raw(3);
  _M0L6_2atmpS2910[0] = _M0L6_2atmpS2911;
  _M0L6_2atmpS2910[1] = _M0L6_2atmpS2912;
  _M0L6_2atmpS2910[2] = _M0L6_2atmpS2913;
  _M0L11join__partsS1094
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE));
  Moonbit_object_header(_M0L11join__partsS1094)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L11join__partsS1094->$0 = _M0L6_2atmpS2910;
  _M0L11join__partsS1094->$1 = 3;
  #line 183 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2903
  = _M0MP26biolab8bio__seq3Seq4join(_M0L6spacerS1093, _M0L11join__partsS1094);
  moonbit_decref(_M0L6spacerS1093);
  moonbit_decref(_M0L11join__partsS1094);
  #line 183 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2902 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2903);
  moonbit_decref(_M0L6_2atmpS2903);
  #line 183 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_86.data, _M0L6_2atmpS2902);
  moonbit_decref(_M0L6_2atmpS2902);
  #line 186 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2905
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_22.data);
  #line 186 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2906
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_22.data);
  #line 186 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2904
  = _M0IP26biolab8bio__seq3SeqPB2Eq5equal(_M0L6_2atmpS2905, _M0L6_2atmpS2906);
  moonbit_decref(_M0L6_2atmpS2905);
  moonbit_decref(_M0L6_2atmpS2906);
  #line 186 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_87.data, _M0L6_2atmpS2904);
  #line 187 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2908
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_22.data);
  #line 187 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2909
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_88.data);
  #line 187 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2907
  = _M0IP26biolab8bio__seq3SeqPB2Eq5equal(_M0L6_2atmpS2908, _M0L6_2atmpS2909);
  moonbit_decref(_M0L6_2atmpS2908);
  moonbit_decref(_M0L6_2atmpS2909);
  #line 187 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main12record__bool(_M0L5casesS1086, (moonbit_string_t)moonbit_string_literal_89.data, _M0L6_2atmpS2907);
  return _M0L5casesS1086;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main19translate__to__stop(
  struct _M0TP26biolab8bio__seq3Seq* _M0L1sS1085
) {
  void* _M0L11_2atry__errS1084;
  struct moonbit_result_1 _tmp_3257;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2799;
  moonbit_string_t _result_3258;
  #line 200 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 201 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3257
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L1sS1085, 42, 1, 0, 45);
  if (_tmp_3257.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2800 =
      _tmp_3257.data.ok;
    _M0L6_2atmpS2799 = _M0L5_2aokS2800;
  } else {
    void* const _M0L6_2aerrS2801 = _tmp_3257.data.err;
    _M0L11_2atry__errS1084 = _M0L6_2aerrS2801;
    goto join_1083;
  }
  #line 201 "/home/developer/Documents/1/cmd/main/main.mbt"
  _result_3258 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2799);
  moonbit_decref(_M0L6_2atmpS2799);
  return _result_3258;
  join_1083:;
  moonbit_decref(_M0L11_2atry__errS1084);
  return (moonbit_string_t)moonbit_string_literal_90.data;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main19translate__stop__at(
  struct _M0TP26biolab8bio__seq3Seq* _M0L1sS1081,
  int32_t _M0L3symS1082
) {
  void* _M0L11_2atry__errS1080;
  struct moonbit_result_1 _tmp_3260;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2796;
  moonbit_string_t _result_3261;
  #line 207 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 208 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3260
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L1sS1081, _M0L3symS1082, 0, 0, 45);
  if (_tmp_3260.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2797 =
      _tmp_3260.data.ok;
    _M0L6_2atmpS2796 = _M0L5_2aokS2797;
  } else {
    void* const _M0L6_2aerrS2798 = _tmp_3260.data.err;
    _M0L11_2atry__errS1080 = _M0L6_2aerrS2798;
    goto join_1079;
  }
  #line 208 "/home/developer/Documents/1/cmd/main/main.mbt"
  _result_3261 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2796);
  moonbit_decref(_M0L6_2atmpS2796);
  return _result_3261;
  join_1079:;
  moonbit_decref(_M0L11_2atry__errS1080);
  return (moonbit_string_t)moonbit_string_literal_90.data;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main18translate__default(
  struct _M0TP26biolab8bio__seq3Seq* _M0L1sS1078
) {
  void* _M0L11_2atry__errS1077;
  struct moonbit_result_1 _tmp_3263;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2793;
  moonbit_string_t _result_3264;
  #line 193 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 194 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3263
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L1sS1078, 42, 0, 0, 45);
  if (_tmp_3263.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2794 =
      _tmp_3263.data.ok;
    _M0L6_2atmpS2793 = _M0L5_2aokS2794;
  } else {
    void* const _M0L6_2aerrS2795 = _tmp_3263.data.err;
    _M0L11_2atry__errS1077 = _M0L6_2aerrS2795;
    goto join_1076;
  }
  #line 194 "/home/developer/Documents/1/cmd/main/main.mbt"
  _result_3264 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2793);
  moonbit_decref(_M0L6_2atmpS2793);
  return _result_3264;
  join_1076:;
  moonbit_decref(_M0L11_2atry__errS1077);
  return (moonbit_string_t)moonbit_string_literal_90.data;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main14translate__cds(
  struct _M0TP26biolab8bio__seq3Seq* _M0L3cdsS1075
) {
  void* _M0L11_2atry__errS1074;
  struct moonbit_result_1 _tmp_3266;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2790;
  moonbit_string_t _result_3267;
  #line 214 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 215 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3266
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L3cdsS1075, 42, 0, 1, 45);
  if (_tmp_3266.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2791 =
      _tmp_3266.data.ok;
    _M0L6_2atmpS2790 = _M0L5_2aokS2791;
  } else {
    void* const _M0L6_2aerrS2792 = _tmp_3266.data.err;
    _M0L11_2atry__errS1074 = _M0L6_2aerrS2792;
    goto join_1073;
  }
  #line 215 "/home/developer/Documents/1/cmd/main/main.mbt"
  _result_3267 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L6_2atmpS2790);
  moonbit_decref(_M0L6_2atmpS2790);
  return _result_3267;
  join_1073:;
  moonbit_decref(_M0L11_2atry__errS1074);
  return (moonbit_string_t)moonbit_string_literal_90.data;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main14invalid__codon() {
  void* _M0L11_2atry__errS1072;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2787;
  struct moonbit_result_1 _tmp_3269;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2786;
  #line 247 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 249 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2787
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_91.data);
  #line 249 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3269
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L6_2atmpS2787, 42, 0, 0, 45);
  moonbit_decref(_M0L6_2atmpS2787);
  if (_tmp_3269.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2788 =
      _tmp_3269.data.ok;
    _M0L6_2atmpS2786 = _M0L5_2aokS2788;
  } else {
    void* const _M0L6_2aerrS2789 = _tmp_3269.data.err;
    _M0L11_2atry__errS1072 = _M0L6_2aerrS2789;
    goto join_1071;
  }
  moonbit_decref(_M0L6_2atmpS2786);
  return (moonbit_string_t)moonbit_string_literal_92.data;
  join_1071:;
  moonbit_decref(_M0L11_2atry__errS1072);
  return (moonbit_string_t)moonbit_string_literal_93.data;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main19cds__internal__stop(
  struct _M0TP26biolab8bio__seq3Seq* _M0L1sS1070
) {
  void* _M0L11_2atry__errS1069;
  struct moonbit_result_1 _tmp_3271;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2783;
  #line 236 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 238 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3271
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L1sS1070, 42, 0, 1, 45);
  if (_tmp_3271.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2784 =
      _tmp_3271.data.ok;
    _M0L6_2atmpS2783 = _M0L5_2aokS2784;
  } else {
    void* const _M0L6_2aerrS2785 = _tmp_3271.data.err;
    _M0L11_2atry__errS1069 = _M0L6_2aerrS2785;
    goto join_1068;
  }
  moonbit_decref(_M0L6_2atmpS2783);
  return (moonbit_string_t)moonbit_string_literal_92.data;
  join_1068:;
  moonbit_decref(_M0L11_2atry__errS1069);
  return (moonbit_string_t)moonbit_string_literal_93.data;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd4main15cds__bad__start() {
  void* _M0L11_2atry__errS1067;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2780;
  struct moonbit_result_1 _tmp_3273;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2779;
  #line 224 "/home/developer/Documents/1/cmd/main/main.mbt"
  #line 226 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0L6_2atmpS2780
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_94.data);
  #line 226 "/home/developer/Documents/1/cmd/main/main.mbt"
  _tmp_3273
  = _M0MP26biolab8bio__seq3Seq17translate_2einner(_M0L6_2atmpS2780, 42, 0, 1, 45);
  moonbit_decref(_M0L6_2atmpS2780);
  if (_tmp_3273.tag) {
    struct _M0TP26biolab8bio__seq3Seq* const _M0L5_2aokS2781 =
      _tmp_3273.data.ok;
    _M0L6_2atmpS2779 = _M0L5_2aokS2781;
  } else {
    void* const _M0L6_2aerrS2782 = _tmp_3273.data.err;
    _M0L11_2atry__errS1067 = _M0L6_2aerrS2782;
    goto join_1066;
  }
  moonbit_decref(_M0L6_2atmpS2779);
  return (moonbit_string_t)moonbit_string_literal_92.data;
  join_1066:;
  moonbit_decref(_M0L11_2atry__errS1067);
  return (moonbit_string_t)moonbit_string_literal_93.data;
}

int32_t _M0FP46biolab8bio__seq3cmd4main12record__bool(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L5casesS1063,
  moonbit_string_t _M0L4nameS1064,
  int32_t _M0L6resultS1065
) {
  moonbit_string_t _M0L6_2atmpS2778;
  #line 22 "/home/developer/Documents/1/cmd/main/main.mbt"
  if (_M0L6resultS1065) {
    _M0L6_2atmpS2778 = (moonbit_string_t)moonbit_string_literal_95.data;
  } else {
    _M0L6_2atmpS2778 = (moonbit_string_t)moonbit_string_literal_96.data;
  }
  #line 23 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0FP46biolab8bio__seq3cmd4main6record(_M0L5casesS1063, _M0L4nameS1064, _M0L6_2atmpS2778);
  moonbit_decref(_M0L6_2atmpS2778);
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd4main6record(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L5casesS1060,
  moonbit_string_t _M0L4nameS1061,
  moonbit_string_t _M0L6resultS1062
) {
  struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L6_2atmpS2777;
  #line 16 "/home/developer/Documents/1/cmd/main/main.mbt"
  moonbit_incref(_M0L4nameS1061);
  moonbit_incref(_M0L6resultS1062);
  _M0L6_2atmpS2777
  = (struct _M0TP46biolab8bio__seq3cmd4main4Case*)moonbit_malloc(sizeof(struct _M0TP46biolab8bio__seq3cmd4main4Case));
  Moonbit_object_header(_M0L6_2atmpS2777)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
  _M0L6_2atmpS2777->$0 = _M0L4nameS1061;
  _M0L6_2atmpS2777->$1 = _M0L6resultS1062;
  #line 17 "/home/developer/Documents/1/cmd/main/main.mbt"
  _M0MPC15array5Array4pushGRP46biolab8bio__seq3cmd4main4CaseE(_M0L5casesS1060, _M0L6_2atmpS2777);
  moonbit_decref(_M0L6_2atmpS2777);
  return 0;
}

struct moonbit_result_1 _M0MP26biolab8bio__seq3Seq17translate_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS1018,
  int32_t _M0L12stop__symbolS1024,
  int32_t _M0L8to__stopS1051,
  int32_t _M0L3cdsS1029,
  int32_t _M0L3gapS1026
) {
  moonbit_string_t _M0L4dataS2776;
  moonbit_string_t _M0L3seqS1017;
  int32_t _M0L1nS1019;
  int32_t _M0L10max__aminoS1020;
  uint16_t* _M0L10amino__bufS1021;
  struct _M0TPB8MutLocalGiE* _M0L10amino__lenS1022;
  int32_t _M0L6_2atmpS2775;
  int32_t _M0L10stop__codeS1023;
  int32_t _M0L6_2atmpS2774;
  int32_t _M0L9gap__codeS1025;
  int32_t _M0L6_2atmpS2773;
  int32_t _M0L7m__codeS1027;
  int32_t _M0L6_2atmpS2772;
  int32_t _M0L7x__codeS1028;
  struct _M0TPB8MutLocalGiE* _M0L1iS1045;
  int32_t _M0L3valS2771;
  uint16_t* _M0L10final__bufS1057;
  int32_t _M0L1kS1058;
  moonbit_string_t _M0L6_2atmpS2770;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2769;
  struct moonbit_result_1 _result_3289;
  #line 407 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2776 = _M0L4selfS1018->$0;
  moonbit_incref(_M0L4dataS2776);
  #line 414 "/home/developer/Documents/1/seq.mbt"
  _M0L3seqS1017 = _M0MPC16string6String9to__upper(_M0L4dataS2776);
  moonbit_decref(_M0L4dataS2776);
  _M0L1nS1019 = Moonbit_array_length(_M0L3seqS1017);
  _M0L10max__aminoS1020 = _M0L1nS1019 / 3;
  _M0L10amino__bufS1021
  = (uint16_t*)moonbit_make_string(_M0L10max__aminoS1020, 0);
  _M0L10amino__lenS1022
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L10amino__lenS1022)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L10amino__lenS1022->$0 = 0;
  _M0L6_2atmpS2775 = _M0L12stop__symbolS1024;
  _M0L10stop__codeS1023 = (uint16_t)_M0L6_2atmpS2775;
  _M0L6_2atmpS2774 = _M0L3gapS1026;
  _M0L9gap__codeS1025 = (uint16_t)_M0L6_2atmpS2774;
  _M0L6_2atmpS2773 = 77;
  _M0L7m__codeS1027 = (uint16_t)_M0L6_2atmpS2773;
  _M0L6_2atmpS2772 = 88;
  _M0L7x__codeS1028 = (uint16_t)_M0L6_2atmpS2772;
  if (_M0L3cdsS1029) {
    int32_t _M0L6_2atmpS2686;
    int32_t _M0L6_2atmpS2690;
    int32_t _M0L6_2atmpS2696;
    int32_t _M0L6_2atmpS2695;
    struct _M0TPB8MutLocalGiE* _M0L1iS1032;
    int32_t _M0L4stopS1033;
    int32_t _M0L3valS2732;
    uint16_t* _M0L10final__bufS1042;
    int32_t _M0L1kS1043;
    moonbit_string_t _M0L6_2atmpS2731;
    struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2730;
    struct moonbit_result_1 _result_3282;
    #line 425 "/home/developer/Documents/1/seq.mbt"
    _M0L6_2atmpS2686
    = _M0FP26biolab8bio__seq20is__start__codon__at(_M0L3seqS1017, 0);
    if (!_M0L6_2atmpS2686) {
      moonbit_string_t _M0L12first__codonS1030;
      moonbit_string_t _M0L6_2atmpS2689;
      moonbit_string_t _M0L6_2atmpS2688;
      void* _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2687;
      struct moonbit_result_1 _result_3274;
      moonbit_decref(_M0L10amino__lenS1022);
      moonbit_decref(_M0L10amino__bufS1021);
      #line 426 "/home/developer/Documents/1/seq.mbt"
      _M0L12first__codonS1030
      = _M0FP26biolab8bio__seq17codon__to__string(_M0L3seqS1017, 0);
      moonbit_decref(_M0L3seqS1017);
      #line 427 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2689
      = moonbit_add_string((moonbit_string_t)moonbit_string_literal_97.data, _M0L12first__codonS1030);
      moonbit_decref(_M0L12first__codonS1030);
      #line 427 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2688
      = moonbit_add_string(_M0L6_2atmpS2689, (moonbit_string_t)moonbit_string_literal_98.data);
      moonbit_decref(_M0L6_2atmpS2689);
      _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2687
      = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError));
      Moonbit_object_header(_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2687)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 1);
      ((struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError*)_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2687)->$0
      = _M0L6_2atmpS2688;
      _result_3274.tag = 0;
      _result_3274.data.err
      = _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2687;
      return _result_3274;
    }
    _M0L6_2atmpS2690 = _M0L1nS1019 % 3;
    if (_M0L6_2atmpS2690 != 0) {
      moonbit_string_t _M0L6_2atmpS2694;
      moonbit_string_t _M0L6_2atmpS2693;
      moonbit_string_t _M0L6_2atmpS2692;
      void* _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2691;
      struct moonbit_result_1 _result_3275;
      moonbit_decref(_M0L10amino__lenS1022);
      moonbit_decref(_M0L10amino__bufS1021);
      moonbit_decref(_M0L3seqS1017);
      #line 431 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2694 = _M0MPC13int3Int18to__string_2einner(_M0L1nS1019, 10);
      #line 431 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2693
      = moonbit_add_string((moonbit_string_t)moonbit_string_literal_99.data, _M0L6_2atmpS2694);
      moonbit_decref(_M0L6_2atmpS2694);
      #line 431 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2692
      = moonbit_add_string(_M0L6_2atmpS2693, (moonbit_string_t)moonbit_string_literal_100.data);
      moonbit_decref(_M0L6_2atmpS2693);
      _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2691
      = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError));
      Moonbit_object_header(_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2691)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 1);
      ((struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError*)_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2691)->$0
      = _M0L6_2atmpS2692;
      _result_3275.tag = 0;
      _result_3275.data.err
      = _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2691;
      return _result_3275;
    }
    _M0L6_2atmpS2696 = _M0L1nS1019 - 3;
    #line 434 "/home/developer/Documents/1/seq.mbt"
    _M0L6_2atmpS2695
    = _M0FP26biolab8bio__seq19is__stop__codon__at(_M0L3seqS1017, _M0L6_2atmpS2696);
    if (!_M0L6_2atmpS2695) {
      int32_t _M0L6_2atmpS2700;
      moonbit_string_t _M0L11last__codonS1031;
      moonbit_string_t _M0L6_2atmpS2699;
      moonbit_string_t _M0L6_2atmpS2698;
      void* _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2697;
      struct moonbit_result_1 _result_3276;
      moonbit_decref(_M0L10amino__lenS1022);
      moonbit_decref(_M0L10amino__bufS1021);
      _M0L6_2atmpS2700 = _M0L1nS1019 - 3;
      #line 435 "/home/developer/Documents/1/seq.mbt"
      _M0L11last__codonS1031
      = _M0FP26biolab8bio__seq17codon__to__string(_M0L3seqS1017, _M0L6_2atmpS2700);
      moonbit_decref(_M0L3seqS1017);
      #line 436 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2699
      = moonbit_add_string((moonbit_string_t)moonbit_string_literal_101.data, _M0L11last__codonS1031);
      moonbit_decref(_M0L11last__codonS1031);
      #line 436 "/home/developer/Documents/1/seq.mbt"
      _M0L6_2atmpS2698
      = moonbit_add_string(_M0L6_2atmpS2699, (moonbit_string_t)moonbit_string_literal_102.data);
      moonbit_decref(_M0L6_2atmpS2699);
      _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2697
      = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError));
      Moonbit_object_header(_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2697)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 1);
      ((struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError*)_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2697)->$0
      = _M0L6_2atmpS2698;
      _result_3276.tag = 0;
      _result_3276.data.err
      = _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2697;
      return _result_3276;
    }
    if (0 < 0 || 0 >= Moonbit_array_length(_M0L10amino__bufS1021)) {
      #line 438 "/home/developer/Documents/1/seq.mbt"
      moonbit_panic();
    }
    _M0L10amino__bufS1021[0] = _M0L7m__codeS1027;
    _M0L10amino__lenS1022->$0 = 1;
    _M0L1iS1032
    = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
    Moonbit_object_header(_M0L1iS1032)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L1iS1032->$0 = 3;
    _M0L4stopS1033 = _M0L1nS1019 - 3;
    while (1) {
      int32_t _M0L3valS2701 = _M0L1iS1032->$0;
      if (_M0L3valS2701 < _M0L4stopS1033) {
        int32_t _M0L3valS2702 = _M0L1iS1032->$0;
        int32_t _M0L2aaS1036;
        int32_t _M0L3valS2724;
        int32_t _M0L7_2abindS1037;
        int32_t _M0L3valS2707;
        int32_t _M0L6_2atmpS2709;
        int32_t _M0L6_2atmpS2708;
        int32_t _M0L3valS2711;
        int32_t _M0L6_2atmpS2710;
        int32_t _M0L3valS2726;
        int32_t _M0L6_2atmpS2725;
        #line 443 "/home/developer/Documents/1/seq.mbt"
        if (
          _M0FP26biolab8bio__seq19is__stop__codon__at(_M0L3seqS1017, _M0L3valS2702)
        ) {
          int32_t _M0L3valS2706;
          moonbit_string_t _M0L5codonS1034;
          moonbit_string_t _M0L6_2atmpS2705;
          moonbit_string_t _M0L6_2atmpS2704;
          void* _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2703;
          struct moonbit_result_1 _result_3278;
          moonbit_decref(_M0L10amino__lenS1022);
          moonbit_decref(_M0L10amino__bufS1021);
          _M0L3valS2706 = _M0L1iS1032->$0;
          moonbit_decref(_M0L1iS1032);
          #line 444 "/home/developer/Documents/1/seq.mbt"
          _M0L5codonS1034
          = _M0FP26biolab8bio__seq17codon__to__string(_M0L3seqS1017, _M0L3valS2706);
          moonbit_decref(_M0L3seqS1017);
          #line 445 "/home/developer/Documents/1/seq.mbt"
          _M0L6_2atmpS2705
          = moonbit_add_string((moonbit_string_t)moonbit_string_literal_103.data, _M0L5codonS1034);
          moonbit_decref(_M0L5codonS1034);
          #line 445 "/home/developer/Documents/1/seq.mbt"
          _M0L6_2atmpS2704
          = moonbit_add_string(_M0L6_2atmpS2705, (moonbit_string_t)moonbit_string_literal_104.data);
          moonbit_decref(_M0L6_2atmpS2705);
          _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2703
          = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError));
          Moonbit_object_header(_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2703)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 1);
          ((struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError*)_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2703)->$0
          = _M0L6_2atmpS2704;
          _result_3278.tag = 0;
          _result_3278.data.err
          = _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2703;
          return _result_3278;
        }
        _M0L3valS2724 = _M0L1iS1032->$0;
        #line 447 "/home/developer/Documents/1/seq.mbt"
        _M0L7_2abindS1037
        = _M0FP26biolab8bio__seq13lookup__codon(_M0L3seqS1017, _M0L3valS2724);
        if (_M0L7_2abindS1037 == -1) {
          int32_t _M0L3valS2712 = _M0L1iS1032->$0;
          #line 453 "/home/developer/Documents/1/seq.mbt"
          if (
            _M0FP26biolab8bio__seq23all__valid__letters__at(_M0L3seqS1017, _M0L3valS2712)
          ) {
            int32_t _M0L3valS2713 = _M0L10amino__lenS1022->$0;
            int32_t _M0L3valS2715;
            int32_t _M0L6_2atmpS2714;
            if (
              _M0L3valS2713 < 0
              || _M0L3valS2713 >= Moonbit_array_length(_M0L10amino__bufS1021)
            ) {
              #line 454 "/home/developer/Documents/1/seq.mbt"
              moonbit_panic();
            }
            _M0L10amino__bufS1021[_M0L3valS2713] = _M0L7x__codeS1028;
            _M0L3valS2715 = _M0L10amino__lenS1022->$0;
            _M0L6_2atmpS2714 = _M0L3valS2715 + 1;
            _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2714;
          } else {
            int32_t _M0L3valS2716 = _M0L1iS1032->$0;
            #line 456 "/home/developer/Documents/1/seq.mbt"
            if (
              _M0FP26biolab8bio__seq18is__gap__codon__at(_M0L3seqS1017, _M0L3valS2716, _M0L3gapS1026)
            ) {
              int32_t _M0L3valS2717 = _M0L10amino__lenS1022->$0;
              int32_t _M0L3valS2719;
              int32_t _M0L6_2atmpS2718;
              if (
                _M0L3valS2717 < 0
                || _M0L3valS2717
                   >= Moonbit_array_length(_M0L10amino__bufS1021)
              ) {
                #line 457 "/home/developer/Documents/1/seq.mbt"
                moonbit_panic();
              }
              _M0L10amino__bufS1021[_M0L3valS2717] = _M0L9gap__codeS1025;
              _M0L3valS2719 = _M0L10amino__lenS1022->$0;
              _M0L6_2atmpS2718 = _M0L3valS2719 + 1;
              _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2718;
            } else {
              int32_t _M0L3valS2723;
              moonbit_string_t _M0L5codonS1040;
              moonbit_string_t _M0L6_2atmpS2722;
              moonbit_string_t _M0L6_2atmpS2721;
              void* _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2720;
              struct moonbit_result_1 _result_3280;
              moonbit_decref(_M0L10amino__lenS1022);
              moonbit_decref(_M0L10amino__bufS1021);
              _M0L3valS2723 = _M0L1iS1032->$0;
              moonbit_decref(_M0L1iS1032);
              #line 460 "/home/developer/Documents/1/seq.mbt"
              _M0L5codonS1040
              = _M0FP26biolab8bio__seq17codon__to__string(_M0L3seqS1017, _M0L3valS2723);
              moonbit_decref(_M0L3seqS1017);
              #line 461 "/home/developer/Documents/1/seq.mbt"
              _M0L6_2atmpS2722
              = moonbit_add_string((moonbit_string_t)moonbit_string_literal_105.data, _M0L5codonS1040);
              moonbit_decref(_M0L5codonS1040);
              #line 461 "/home/developer/Documents/1/seq.mbt"
              _M0L6_2atmpS2721
              = moonbit_add_string(_M0L6_2atmpS2722, (moonbit_string_t)moonbit_string_literal_106.data);
              moonbit_decref(_M0L6_2atmpS2722);
              _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2720
              = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError));
              Moonbit_object_header(_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2720)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 1);
              ((struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError*)_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2720)->$0
              = _M0L6_2atmpS2721;
              _result_3280.tag = 0;
              _result_3280.data.err
              = _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2720;
              return _result_3280;
            }
          }
        } else {
          int32_t _M0L7_2aSomeS1038 = _M0L7_2abindS1037;
          int32_t _M0L5_2aaaS1039 = _M0L7_2aSomeS1038;
          _M0L2aaS1036 = _M0L5_2aaaS1039;
          goto join_1035;
        }
        goto joinlet_3279;
        join_1035:;
        _M0L3valS2707 = _M0L10amino__lenS1022->$0;
        _M0L6_2atmpS2709 = _M0L2aaS1036;
        _M0L6_2atmpS2708 = (uint16_t)_M0L6_2atmpS2709;
        if (
          _M0L3valS2707 < 0
          || _M0L3valS2707 >= Moonbit_array_length(_M0L10amino__bufS1021)
        ) {
          #line 449 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L10amino__bufS1021[_M0L3valS2707] = _M0L6_2atmpS2708;
        _M0L3valS2711 = _M0L10amino__lenS1022->$0;
        _M0L6_2atmpS2710 = _M0L3valS2711 + 1;
        _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2710;
        joinlet_3279:;
        _M0L3valS2726 = _M0L1iS1032->$0;
        _M0L6_2atmpS2725 = _M0L3valS2726 + 3;
        _M0L1iS1032->$0 = _M0L6_2atmpS2725;
        continue;
      } else {
        moonbit_decref(_M0L1iS1032);
        moonbit_decref(_M0L3seqS1017);
      }
      break;
    }
    _M0L3valS2732 = _M0L10amino__lenS1022->$0;
    _M0L10final__bufS1042 = (uint16_t*)moonbit_make_string(_M0L3valS2732, 0);
    _M0L1kS1043 = 0;
    while (1) {
      int32_t _M0L3valS2727 = _M0L10amino__lenS1022->$0;
      if (_M0L1kS1043 < _M0L3valS2727) {
        int32_t _M0L6_2atmpS2728;
        int32_t _M0L6_2atmpS2729;
        if (
          _M0L1kS1043 < 0
          || _M0L1kS1043 >= Moonbit_array_length(_M0L10amino__bufS1021)
        ) {
          #line 468 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS2728 = (int32_t)_M0L10amino__bufS1021[_M0L1kS1043];
        if (
          _M0L1kS1043 < 0
          || _M0L1kS1043 >= Moonbit_array_length(_M0L10final__bufS1042)
        ) {
          #line 468 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L10final__bufS1042[_M0L1kS1043] = _M0L6_2atmpS2728;
        _M0L6_2atmpS2729 = _M0L1kS1043 + 1;
        _M0L1kS1043 = _M0L6_2atmpS2729;
        continue;
      } else {
        moonbit_decref(_M0L10amino__lenS1022);
        moonbit_decref(_M0L10amino__bufS1021);
      }
      break;
    }
    _M0L6_2atmpS2731 = _M0L10final__bufS1042;
    _M0L6_2atmpS2730
    = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
    Moonbit_object_header(_M0L6_2atmpS2730)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
    _M0L6_2atmpS2730->$0 = _M0L6_2atmpS2731;
    _result_3282.tag = 1;
    _result_3282.data.ok = _M0L6_2atmpS2730;
    return _result_3282;
  }
  _M0L1iS1045
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1045)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1045->$0 = 0;
  while (1) {
    int32_t _M0L3valS2733 = _M0L1iS1045->$0;
    int32_t _M0L6_2atmpS2735 = _M0L1nS1019 % 3;
    int32_t _M0L6_2atmpS2734 = _M0L1nS1019 - _M0L6_2atmpS2735;
    if (_M0L3valS2733 < _M0L6_2atmpS2734) {
      int32_t _M0L2aaS1047;
      int32_t _M0L3valS2763 = _M0L1iS1045->$0;
      int32_t _M0L7_2abindS1048;
      int32_t _M0L3valS2736;
      int32_t _M0L6_2atmpS2738;
      int32_t _M0L6_2atmpS2737;
      int32_t _M0L3valS2740;
      int32_t _M0L6_2atmpS2739;
      int32_t _M0L3valS2765;
      int32_t _M0L6_2atmpS2764;
      #line 474 "/home/developer/Documents/1/seq.mbt"
      _M0L7_2abindS1048
      = _M0FP26biolab8bio__seq13lookup__codon(_M0L3seqS1017, _M0L3valS2763);
      if (_M0L7_2abindS1048 == -1) {
        int32_t _M0L3valS2741 = _M0L1iS1045->$0;
        #line 480 "/home/developer/Documents/1/seq.mbt"
        if (
          _M0FP26biolab8bio__seq19is__stop__codon__at(_M0L3seqS1017, _M0L3valS2741)
        ) {
          int32_t _M0L3valS2748;
          int32_t _M0L3valS2750;
          int32_t _M0L6_2atmpS2749;
          if (_M0L8to__stopS1051) {
            int32_t _M0L3valS2747;
            uint16_t* _M0L10final__bufS1052;
            int32_t _M0L1kS1053;
            moonbit_string_t _M0L6_2atmpS2746;
            struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2745;
            struct moonbit_result_1 _result_3286;
            moonbit_decref(_M0L1iS1045);
            moonbit_decref(_M0L3seqS1017);
            _M0L3valS2747 = _M0L10amino__lenS1022->$0;
            _M0L10final__bufS1052
            = (uint16_t*)moonbit_make_string(_M0L3valS2747, 0);
            _M0L1kS1053 = 0;
            while (1) {
              int32_t _M0L3valS2742 = _M0L10amino__lenS1022->$0;
              if (_M0L1kS1053 < _M0L3valS2742) {
                int32_t _M0L6_2atmpS2743;
                int32_t _M0L6_2atmpS2744;
                if (
                  _M0L1kS1053 < 0
                  || _M0L1kS1053
                     >= Moonbit_array_length(_M0L10amino__bufS1021)
                ) {
                  #line 484 "/home/developer/Documents/1/seq.mbt"
                  moonbit_panic();
                }
                _M0L6_2atmpS2743
                = (int32_t)_M0L10amino__bufS1021[_M0L1kS1053];
                if (
                  _M0L1kS1053 < 0
                  || _M0L1kS1053
                     >= Moonbit_array_length(_M0L10final__bufS1052)
                ) {
                  #line 484 "/home/developer/Documents/1/seq.mbt"
                  moonbit_panic();
                }
                _M0L10final__bufS1052[_M0L1kS1053] = _M0L6_2atmpS2743;
                _M0L6_2atmpS2744 = _M0L1kS1053 + 1;
                _M0L1kS1053 = _M0L6_2atmpS2744;
                continue;
              } else {
                moonbit_decref(_M0L10amino__lenS1022);
                moonbit_decref(_M0L10amino__bufS1021);
              }
              break;
            }
            _M0L6_2atmpS2746 = _M0L10final__bufS1052;
            _M0L6_2atmpS2745
            = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
            Moonbit_object_header(_M0L6_2atmpS2745)->meta
            = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
            _M0L6_2atmpS2745->$0 = _M0L6_2atmpS2746;
            _result_3286.tag = 1;
            _result_3286.data.ok = _M0L6_2atmpS2745;
            return _result_3286;
          }
          _M0L3valS2748 = _M0L10amino__lenS1022->$0;
          if (
            _M0L3valS2748 < 0
            || _M0L3valS2748 >= Moonbit_array_length(_M0L10amino__bufS1021)
          ) {
            #line 488 "/home/developer/Documents/1/seq.mbt"
            moonbit_panic();
          }
          _M0L10amino__bufS1021[_M0L3valS2748] = _M0L10stop__codeS1023;
          _M0L3valS2750 = _M0L10amino__lenS1022->$0;
          _M0L6_2atmpS2749 = _M0L3valS2750 + 1;
          _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2749;
        } else {
          int32_t _M0L3valS2751 = _M0L1iS1045->$0;
          #line 490 "/home/developer/Documents/1/seq.mbt"
          if (
            _M0FP26biolab8bio__seq23all__valid__letters__at(_M0L3seqS1017, _M0L3valS2751)
          ) {
            int32_t _M0L3valS2752 = _M0L10amino__lenS1022->$0;
            int32_t _M0L3valS2754;
            int32_t _M0L6_2atmpS2753;
            if (
              _M0L3valS2752 < 0
              || _M0L3valS2752 >= Moonbit_array_length(_M0L10amino__bufS1021)
            ) {
              #line 491 "/home/developer/Documents/1/seq.mbt"
              moonbit_panic();
            }
            _M0L10amino__bufS1021[_M0L3valS2752] = _M0L7x__codeS1028;
            _M0L3valS2754 = _M0L10amino__lenS1022->$0;
            _M0L6_2atmpS2753 = _M0L3valS2754 + 1;
            _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2753;
          } else {
            int32_t _M0L3valS2755 = _M0L1iS1045->$0;
            #line 493 "/home/developer/Documents/1/seq.mbt"
            if (
              _M0FP26biolab8bio__seq18is__gap__codon__at(_M0L3seqS1017, _M0L3valS2755, _M0L3gapS1026)
            ) {
              int32_t _M0L3valS2756 = _M0L10amino__lenS1022->$0;
              int32_t _M0L3valS2758;
              int32_t _M0L6_2atmpS2757;
              if (
                _M0L3valS2756 < 0
                || _M0L3valS2756
                   >= Moonbit_array_length(_M0L10amino__bufS1021)
              ) {
                #line 494 "/home/developer/Documents/1/seq.mbt"
                moonbit_panic();
              }
              _M0L10amino__bufS1021[_M0L3valS2756] = _M0L9gap__codeS1025;
              _M0L3valS2758 = _M0L10amino__lenS1022->$0;
              _M0L6_2atmpS2757 = _M0L3valS2758 + 1;
              _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2757;
            } else {
              int32_t _M0L3valS2762;
              moonbit_string_t _M0L5codonS1055;
              moonbit_string_t _M0L6_2atmpS2761;
              moonbit_string_t _M0L6_2atmpS2760;
              void* _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2759;
              struct moonbit_result_1 _result_3287;
              moonbit_decref(_M0L10amino__lenS1022);
              moonbit_decref(_M0L10amino__bufS1021);
              _M0L3valS2762 = _M0L1iS1045->$0;
              moonbit_decref(_M0L1iS1045);
              #line 497 "/home/developer/Documents/1/seq.mbt"
              _M0L5codonS1055
              = _M0FP26biolab8bio__seq17codon__to__string(_M0L3seqS1017, _M0L3valS2762);
              moonbit_decref(_M0L3seqS1017);
              #line 498 "/home/developer/Documents/1/seq.mbt"
              _M0L6_2atmpS2761
              = moonbit_add_string((moonbit_string_t)moonbit_string_literal_105.data, _M0L5codonS1055);
              moonbit_decref(_M0L5codonS1055);
              #line 498 "/home/developer/Documents/1/seq.mbt"
              _M0L6_2atmpS2760
              = moonbit_add_string(_M0L6_2atmpS2761, (moonbit_string_t)moonbit_string_literal_106.data);
              moonbit_decref(_M0L6_2atmpS2761);
              _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2759
              = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError));
              Moonbit_object_header(_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2759)->meta
              = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 1);
              ((struct _M0DTPC15error5Error39biolab_2fbio__seq_2eSeqError_2eSeqError*)_M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2759)->$0
              = _M0L6_2atmpS2760;
              _result_3287.tag = 0;
              _result_3287.data.err
              = _M0L39biolab_2fbio__seq_2eSeqError_2eSeqErrorS2759;
              return _result_3287;
            }
          }
        }
      } else {
        int32_t _M0L7_2aSomeS1049 = _M0L7_2abindS1048;
        int32_t _M0L5_2aaaS1050 = _M0L7_2aSomeS1049;
        _M0L2aaS1047 = _M0L5_2aaaS1050;
        goto join_1046;
      }
      goto joinlet_3284;
      join_1046:;
      _M0L3valS2736 = _M0L10amino__lenS1022->$0;
      _M0L6_2atmpS2738 = _M0L2aaS1047;
      _M0L6_2atmpS2737 = (uint16_t)_M0L6_2atmpS2738;
      if (
        _M0L3valS2736 < 0
        || _M0L3valS2736 >= Moonbit_array_length(_M0L10amino__bufS1021)
      ) {
        #line 476 "/home/developer/Documents/1/seq.mbt"
        moonbit_panic();
      }
      _M0L10amino__bufS1021[_M0L3valS2736] = _M0L6_2atmpS2737;
      _M0L3valS2740 = _M0L10amino__lenS1022->$0;
      _M0L6_2atmpS2739 = _M0L3valS2740 + 1;
      _M0L10amino__lenS1022->$0 = _M0L6_2atmpS2739;
      joinlet_3284:;
      _M0L3valS2765 = _M0L1iS1045->$0;
      _M0L6_2atmpS2764 = _M0L3valS2765 + 3;
      _M0L1iS1045->$0 = _M0L6_2atmpS2764;
      continue;
    } else {
      moonbit_decref(_M0L1iS1045);
      moonbit_decref(_M0L3seqS1017);
    }
    break;
  }
  _M0L3valS2771 = _M0L10amino__lenS1022->$0;
  _M0L10final__bufS1057 = (uint16_t*)moonbit_make_string(_M0L3valS2771, 0);
  _M0L1kS1058 = 0;
  while (1) {
    int32_t _M0L3valS2766 = _M0L10amino__lenS1022->$0;
    if (_M0L1kS1058 < _M0L3valS2766) {
      int32_t _M0L6_2atmpS2767;
      int32_t _M0L6_2atmpS2768;
      if (
        _M0L1kS1058 < 0
        || _M0L1kS1058 >= Moonbit_array_length(_M0L10amino__bufS1021)
      ) {
        #line 505 "/home/developer/Documents/1/seq.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2767 = (int32_t)_M0L10amino__bufS1021[_M0L1kS1058];
      if (
        _M0L1kS1058 < 0
        || _M0L1kS1058 >= Moonbit_array_length(_M0L10final__bufS1057)
      ) {
        #line 505 "/home/developer/Documents/1/seq.mbt"
        moonbit_panic();
      }
      _M0L10final__bufS1057[_M0L1kS1058] = _M0L6_2atmpS2767;
      _M0L6_2atmpS2768 = _M0L1kS1058 + 1;
      _M0L1kS1058 = _M0L6_2atmpS2768;
      continue;
    } else {
      moonbit_decref(_M0L10amino__lenS1022);
      moonbit_decref(_M0L10amino__bufS1021);
    }
    break;
  }
  _M0L6_2atmpS2770 = _M0L10final__bufS1057;
  _M0L6_2atmpS2769
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_M0L6_2atmpS2769)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _M0L6_2atmpS2769->$0 = _M0L6_2atmpS2770;
  _result_3289.tag = 1;
  _result_3289.data.ok = _M0L6_2atmpS2769;
  return _result_3289;
}

int32_t _M0FP26biolab8bio__seq13lookup__codon(
  moonbit_string_t _M0L3seqS1015,
  int32_t _M0L1iS1016
) {
  int32_t _M0L3keyS1014;
  #line 554 "/home/developer/Documents/1/seq.mbt"
  #line 555 "/home/developer/Documents/1/seq.mbt"
  _M0L3keyS1014
  = _M0FP26biolab8bio__seq10codon__key(_M0L3seqS1015, _M0L1iS1016);
  #line 556 "/home/developer/Documents/1/seq.mbt"
  return _M0MPB3Map3getGicE(_M0FP26biolab8bio__seq21codon__forward__table, _M0L3keyS1014);
}

int32_t _M0FP26biolab8bio__seq10codon__key(
  moonbit_string_t _M0L3seqS1010,
  int32_t _M0L1iS1011
) {
  int32_t _M0L6_2atmpS2685;
  int32_t _M0L2c0S1009;
  int32_t _M0L6_2atmpS2684;
  int32_t _M0L6_2atmpS2683;
  int32_t _M0L2c1S1012;
  int32_t _M0L6_2atmpS2682;
  int32_t _M0L6_2atmpS2681;
  int32_t _M0L2c2S1013;
  int32_t _M0L6_2atmpS2679;
  int32_t _M0L6_2atmpS2680;
  int32_t _M0L6_2atmpS2678;
  #line 561 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2685 = _M0L3seqS1010[_M0L1iS1011];
  _M0L2c0S1009 = (int32_t)_M0L6_2atmpS2685;
  _M0L6_2atmpS2684 = _M0L1iS1011 + 1;
  _M0L6_2atmpS2683 = _M0L3seqS1010[_M0L6_2atmpS2684];
  _M0L2c1S1012 = (int32_t)_M0L6_2atmpS2683;
  _M0L6_2atmpS2682 = _M0L1iS1011 + 2;
  _M0L6_2atmpS2681 = _M0L3seqS1010[_M0L6_2atmpS2682];
  _M0L2c2S1013 = (int32_t)_M0L6_2atmpS2681;
  _M0L6_2atmpS2679 = _M0L2c0S1009 * 65536;
  _M0L6_2atmpS2680 = _M0L2c1S1012 * 256;
  _M0L6_2atmpS2678 = _M0L6_2atmpS2679 + _M0L6_2atmpS2680;
  return _M0L6_2atmpS2678 + _M0L2c2S1013;
}

int32_t _M0FP26biolab8bio__seq19is__stop__codon__at(
  moonbit_string_t _M0L3seqS995,
  int32_t _M0L1iS996
) {
  int32_t _M0L6_2atmpS2677;
  int32_t _M0L2c0S994;
  int32_t _M0L6_2atmpS2676;
  int32_t _M0L6_2atmpS2675;
  int32_t _M0L2c1S997;
  int32_t _M0L6_2atmpS2674;
  int32_t _M0L6_2atmpS2673;
  int32_t _M0L2c2S998;
  int32_t _M0L1tS999;
  int32_t _M0L1aS1000;
  int32_t _M0L1gS1001;
  int32_t _if__result_3290;
  #line 524 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2677 = _M0L3seqS995[_M0L1iS996];
  _M0L2c0S994 = (int32_t)_M0L6_2atmpS2677;
  _M0L6_2atmpS2676 = _M0L1iS996 + 1;
  _M0L6_2atmpS2675 = _M0L3seqS995[_M0L6_2atmpS2676];
  _M0L2c1S997 = (int32_t)_M0L6_2atmpS2675;
  _M0L6_2atmpS2674 = _M0L1iS996 + 2;
  _M0L6_2atmpS2673 = _M0L3seqS995[_M0L6_2atmpS2674];
  _M0L2c2S998 = (int32_t)_M0L6_2atmpS2673;
  _M0L1tS999 = 84;
  _M0L1aS1000 = 65;
  _M0L1gS1001 = 71;
  if (_M0L2c0S994 == _M0L1tS999) {
    if (_M0L2c1S997 == _M0L1aS1000) {
      if (_M0L2c2S998 == _M0L1aS1000) {
        _if__result_3290 = 1;
      } else {
        _if__result_3290 = _M0L2c2S998 == _M0L1gS1001;
      }
    } else {
      _if__result_3290 = 0;
    }
  } else {
    _if__result_3290 = 0;
  }
  if (_if__result_3290) {
    return 1;
  } else if (_M0L2c0S994 == _M0L1tS999) {
    if (_M0L2c1S997 == _M0L1gS1001) {
      return _M0L2c2S998 == _M0L1aS1000;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

int32_t _M0FP26biolab8bio__seq20is__start__codon__at(
  moonbit_string_t _M0L3seqS986,
  int32_t _M0L1iS987
) {
  int32_t _M0L6_2atmpS2672;
  int32_t _M0L2c0S985;
  int32_t _M0L6_2atmpS2671;
  int32_t _M0L6_2atmpS2670;
  int32_t _M0L2c1S988;
  int32_t _M0L6_2atmpS2669;
  int32_t _M0L6_2atmpS2668;
  int32_t _M0L2c2S989;
  int32_t _M0L1aS990;
  int32_t _M0L1cS991;
  int32_t _M0L1tS992;
  int32_t _M0L1gS993;
  int32_t _if__result_3291;
  #line 539 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2672 = _M0L3seqS986[_M0L1iS987];
  _M0L2c0S985 = (int32_t)_M0L6_2atmpS2672;
  _M0L6_2atmpS2671 = _M0L1iS987 + 1;
  _M0L6_2atmpS2670 = _M0L3seqS986[_M0L6_2atmpS2671];
  _M0L2c1S988 = (int32_t)_M0L6_2atmpS2670;
  _M0L6_2atmpS2669 = _M0L1iS987 + 2;
  _M0L6_2atmpS2668 = _M0L3seqS986[_M0L6_2atmpS2669];
  _M0L2c2S989 = (int32_t)_M0L6_2atmpS2668;
  _M0L1aS990 = 65;
  _M0L1cS991 = 67;
  _M0L1tS992 = 84;
  _M0L1gS993 = 71;
  if (_M0L2c0S985 == _M0L1aS990) {
    _if__result_3291 = 1;
  } else if (_M0L2c0S985 == _M0L1cS991) {
    _if__result_3291 = 1;
  } else {
    _if__result_3291 = _M0L2c0S985 == _M0L1tS992;
  }
  if (_if__result_3291) {
    if (_M0L2c1S988 == _M0L1tS992) {
      return _M0L2c2S989 == _M0L1gS993;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

int32_t _M0FP26biolab8bio__seq18is__gap__codon__at(
  moonbit_string_t _M0L3seqS983,
  int32_t _M0L1iS984,
  int32_t _M0L3gapS982
) {
  int32_t _M0L1gS981;
  int32_t _M0L6_2atmpS2667;
  int32_t _M0L6_2atmpS2666;
  #line 579 "/home/developer/Documents/1/seq.mbt"
  _M0L1gS981 = _M0L3gapS982;
  _M0L6_2atmpS2667 = _M0L3seqS983[_M0L1iS984];
  _M0L6_2atmpS2666 = (int32_t)_M0L6_2atmpS2667;
  if (_M0L6_2atmpS2666 == _M0L1gS981) {
    int32_t _M0L6_2atmpS2665 = _M0L1iS984 + 1;
    int32_t _M0L6_2atmpS2664 = _M0L3seqS983[_M0L6_2atmpS2665];
    int32_t _M0L6_2atmpS2663 = (int32_t)_M0L6_2atmpS2664;
    if (_M0L6_2atmpS2663 == _M0L1gS981) {
      int32_t _M0L6_2atmpS2662 = _M0L1iS984 + 2;
      int32_t _M0L6_2atmpS2661 = _M0L3seqS983[_M0L6_2atmpS2662];
      int32_t _M0L6_2atmpS2660 = (int32_t)_M0L6_2atmpS2661;
      return _M0L6_2atmpS2660 == _M0L1gS981;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

moonbit_string_t _M0FP26biolab8bio__seq17codon__to__string(
  moonbit_string_t _M0L3seqS979,
  int32_t _M0L1iS980
) {
  struct _M0TPB13StringBuilder* _M0L3bufS978;
  int32_t _M0L6_2atmpS2653;
  int32_t _M0L6_2atmpS2652;
  int32_t _M0L6_2atmpS2656;
  int32_t _M0L6_2atmpS2655;
  int32_t _M0L6_2atmpS2654;
  int32_t _M0L6_2atmpS2659;
  int32_t _M0L6_2atmpS2658;
  int32_t _M0L6_2atmpS2657;
  moonbit_string_t _result_3292;
  #line 513 "/home/developer/Documents/1/seq.mbt"
  #line 514 "/home/developer/Documents/1/seq.mbt"
  _M0L3bufS978 = _M0MPB13StringBuilder21StringBuilder_2einner(3);
  _M0L6_2atmpS2653 = _M0L3seqS979[_M0L1iS980];
  #line 515 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2652
  = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2653);
  #line 515 "/home/developer/Documents/1/seq.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS978, _M0L6_2atmpS2652);
  _M0L6_2atmpS2656 = _M0L1iS980 + 1;
  _M0L6_2atmpS2655 = _M0L3seqS979[_M0L6_2atmpS2656];
  #line 516 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2654
  = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2655);
  #line 516 "/home/developer/Documents/1/seq.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS978, _M0L6_2atmpS2654);
  _M0L6_2atmpS2659 = _M0L1iS980 + 2;
  _M0L6_2atmpS2658 = _M0L3seqS979[_M0L6_2atmpS2659];
  #line 517 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2657
  = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2658);
  #line 517 "/home/developer/Documents/1/seq.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS978, _M0L6_2atmpS2657);
  #line 518 "/home/developer/Documents/1/seq.mbt"
  _result_3292 = _M0MPB13StringBuilder10to__string(_M0L3bufS978);
  moonbit_decref(_M0L3bufS978);
  return _result_3292;
}

int32_t _M0FP26biolab8bio__seq23all__valid__letters__at(
  moonbit_string_t _M0L3seqS976,
  int32_t _M0L1iS977
) {
  int32_t _M0L6_2atmpS2651;
  int32_t _M0L6_2atmpS2650;
  #line 571 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2651 = _M0L3seqS976[_M0L1iS977];
  _M0L6_2atmpS2650 = (int32_t)_M0L6_2atmpS2651;
  #line 572 "/home/developer/Documents/1/seq.mbt"
  if (_M0FP26biolab8bio__seq13valid__letter(_M0L6_2atmpS2650)) {
    int32_t _M0L6_2atmpS2649 = _M0L1iS977 + 1;
    int32_t _M0L6_2atmpS2648 = _M0L3seqS976[_M0L6_2atmpS2649];
    int32_t _M0L6_2atmpS2647 = (int32_t)_M0L6_2atmpS2648;
    #line 573 "/home/developer/Documents/1/seq.mbt"
    if (_M0FP26biolab8bio__seq13valid__letter(_M0L6_2atmpS2647)) {
      int32_t _M0L6_2atmpS2646 = _M0L1iS977 + 2;
      int32_t _M0L6_2atmpS2645 = _M0L3seqS976[_M0L6_2atmpS2646];
      int32_t _M0L6_2atmpS2644 = (int32_t)_M0L6_2atmpS2645;
      #line 574 "/home/developer/Documents/1/seq.mbt"
      return _M0FP26biolab8bio__seq13valid__letter(_M0L6_2atmpS2644);
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

int32_t _M0FP26biolab8bio__seq13valid__letter(int32_t _M0L1cS975) {
  #line 119 "/home/developer/Documents/1/codon_table.mbt"
  if (_M0L1cS975 < 128) {
    if (
      _M0L1cS975 < 0
      || _M0L1cS975
         >= Moonbit_array_length(_M0FP26biolab8bio__seq20valid__letter__table)
    ) {
      #line 121 "/home/developer/Documents/1/codon_table.mbt"
      moonbit_panic();
    }
    return (int32_t)_M0FP26biolab8bio__seq20valid__letter__table[_M0L1cS975];
  } else {
    return 0;
  }
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq16back__transcribe(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS957
) {
  moonbit_string_t _M0L3srcS956;
  int32_t _M0L1nS958;
  uint16_t* _M0L3bufS959;
  int32_t _M0L6_2atmpS2643;
  int32_t _M0L1uS960;
  int32_t _M0L6_2atmpS2642;
  int32_t _M0L2ulS961;
  int32_t _M0L6_2atmpS2641;
  int32_t _M0L1tS962;
  int32_t _M0L6_2atmpS2640;
  int32_t _M0L2tlS963;
  int32_t _M0L1iS964;
  moonbit_string_t _M0L6_2atmpS2639;
  struct _M0TP26biolab8bio__seq3Seq* _block_3294;
  #line 373 "/home/developer/Documents/1/seq.mbt"
  _M0L3srcS956 = _M0L4selfS957->$0;
  _M0L1nS958 = Moonbit_array_length(_M0L3srcS956);
  _M0L3bufS959 = (uint16_t*)moonbit_make_string(_M0L1nS958, 0);
  _M0L6_2atmpS2643 = 85;
  _M0L1uS960 = (uint16_t)_M0L6_2atmpS2643;
  _M0L6_2atmpS2642 = 117;
  _M0L2ulS961 = (uint16_t)_M0L6_2atmpS2642;
  _M0L6_2atmpS2641 = 84;
  _M0L1tS962 = (uint16_t)_M0L6_2atmpS2641;
  _M0L6_2atmpS2640 = 116;
  _M0L2tlS963 = (uint16_t)_M0L6_2atmpS2640;
  moonbit_incref(_M0L3srcS956);
  _M0L1iS964 = 0;
  while (1) {
    if (_M0L1iS964 < _M0L1nS958) {
      int32_t _M0L2ciS965 = _M0L3srcS956[_M0L1iS964];
      int32_t _M0L6_2atmpS2638;
      if (_M0L2ciS965 == _M0L1uS960) {
        if (
          _M0L1iS964 < 0 || _M0L1iS964 >= Moonbit_array_length(_M0L3bufS959)
        ) {
          #line 384 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L3bufS959[_M0L1iS964] = _M0L1tS962;
      } else if (_M0L2ciS965 == _M0L2ulS961) {
        if (
          _M0L1iS964 < 0 || _M0L1iS964 >= Moonbit_array_length(_M0L3bufS959)
        ) {
          #line 386 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L3bufS959[_M0L1iS964] = _M0L2tlS963;
      } else {
        if (
          _M0L1iS964 < 0 || _M0L1iS964 >= Moonbit_array_length(_M0L3bufS959)
        ) {
          #line 388 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L3bufS959[_M0L1iS964] = _M0L2ciS965;
      }
      _M0L6_2atmpS2638 = _M0L1iS964 + 1;
      _M0L1iS964 = _M0L6_2atmpS2638;
      continue;
    } else {
      moonbit_decref(_M0L3srcS956);
    }
    break;
  }
  _M0L6_2atmpS2639 = _M0L3bufS959;
  _block_3294
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3294)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3294->$0 = _M0L6_2atmpS2639;
  return _block_3294;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq10transcribe(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS946
) {
  moonbit_string_t _M0L3srcS945;
  int32_t _M0L1nS947;
  uint16_t* _M0L3bufS948;
  int32_t _M0L6_2atmpS2637;
  int32_t _M0L1tS949;
  int32_t _M0L6_2atmpS2636;
  int32_t _M0L2tlS950;
  int32_t _M0L6_2atmpS2635;
  int32_t _M0L1uS951;
  int32_t _M0L6_2atmpS2634;
  int32_t _M0L2ulS952;
  int32_t _M0L1iS953;
  moonbit_string_t _M0L6_2atmpS2633;
  struct _M0TP26biolab8bio__seq3Seq* _block_3296;
  #line 349 "/home/developer/Documents/1/seq.mbt"
  _M0L3srcS945 = _M0L4selfS946->$0;
  _M0L1nS947 = Moonbit_array_length(_M0L3srcS945);
  _M0L3bufS948 = (uint16_t*)moonbit_make_string(_M0L1nS947, 0);
  _M0L6_2atmpS2637 = 84;
  _M0L1tS949 = (uint16_t)_M0L6_2atmpS2637;
  _M0L6_2atmpS2636 = 116;
  _M0L2tlS950 = (uint16_t)_M0L6_2atmpS2636;
  _M0L6_2atmpS2635 = 85;
  _M0L1uS951 = (uint16_t)_M0L6_2atmpS2635;
  _M0L6_2atmpS2634 = 117;
  _M0L2ulS952 = (uint16_t)_M0L6_2atmpS2634;
  moonbit_incref(_M0L3srcS945);
  _M0L1iS953 = 0;
  while (1) {
    if (_M0L1iS953 < _M0L1nS947) {
      int32_t _M0L2ciS954 = _M0L3srcS945[_M0L1iS953];
      int32_t _M0L6_2atmpS2632;
      if (_M0L2ciS954 == _M0L1tS949) {
        if (
          _M0L1iS953 < 0 || _M0L1iS953 >= Moonbit_array_length(_M0L3bufS948)
        ) {
          #line 360 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L3bufS948[_M0L1iS953] = _M0L1uS951;
      } else if (_M0L2ciS954 == _M0L2tlS950) {
        if (
          _M0L1iS953 < 0 || _M0L1iS953 >= Moonbit_array_length(_M0L3bufS948)
        ) {
          #line 362 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L3bufS948[_M0L1iS953] = _M0L2ulS952;
      } else {
        if (
          _M0L1iS953 < 0 || _M0L1iS953 >= Moonbit_array_length(_M0L3bufS948)
        ) {
          #line 364 "/home/developer/Documents/1/seq.mbt"
          moonbit_panic();
        }
        _M0L3bufS948[_M0L1iS953] = _M0L2ciS954;
      }
      _M0L6_2atmpS2632 = _M0L1iS953 + 1;
      _M0L1iS953 = _M0L6_2atmpS2632;
      continue;
    } else {
      moonbit_decref(_M0L3srcS945);
    }
    break;
  }
  _M0L6_2atmpS2633 = _M0L3bufS948;
  _block_3296
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3296)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3296->$0 = _M0L6_2atmpS2633;
  return _block_3296;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq24reverse__complement__rna(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS944
) {
  moonbit_string_t _M0L4dataS2631;
  moonbit_string_t _M0L6_2atmpS2630;
  struct _M0TP26biolab8bio__seq3Seq* _block_3297;
  #line 341 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2631 = _M0L4selfS944->$0;
  moonbit_incref(_M0L4dataS2631);
  #line 342 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2630
  = _M0FP26biolab8bio__seq26reverse__complement__chars(_M0L4dataS2631, _M0FP26biolab8bio__seq20rna__complement__arr);
  moonbit_decref(_M0L4dataS2631);
  _block_3297
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3297)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3297->$0 = _M0L6_2atmpS2630;
  return _block_3297;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq7reverse(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS939
) {
  moonbit_string_t _M0L3srcS938;
  int32_t _M0L1nS940;
  uint16_t* _M0L3bufS941;
  int32_t _M0L1iS942;
  moonbit_string_t _M0L6_2atmpS2629;
  struct _M0TP26biolab8bio__seq3Seq* _block_3299;
  #line 321 "/home/developer/Documents/1/seq.mbt"
  _M0L3srcS938 = _M0L4selfS939->$0;
  _M0L1nS940 = Moonbit_array_length(_M0L3srcS938);
  _M0L3bufS941 = (uint16_t*)moonbit_make_string(_M0L1nS940, 0);
  moonbit_incref(_M0L3srcS938);
  _M0L1iS942 = 0;
  while (1) {
    if (_M0L1iS942 < _M0L1nS940) {
      int32_t _M0L6_2atmpS2627 = _M0L1nS940 - 1;
      int32_t _M0L6_2atmpS2626 = _M0L6_2atmpS2627 - _M0L1iS942;
      int32_t _M0L6_2atmpS2625 = _M0L3srcS938[_M0L6_2atmpS2626];
      int32_t _M0L6_2atmpS2628;
      if (_M0L1iS942 < 0 || _M0L1iS942 >= Moonbit_array_length(_M0L3bufS941)) {
        #line 326 "/home/developer/Documents/1/seq.mbt"
        moonbit_panic();
      }
      _M0L3bufS941[_M0L1iS942] = _M0L6_2atmpS2625;
      _M0L6_2atmpS2628 = _M0L1iS942 + 1;
      _M0L1iS942 = _M0L6_2atmpS2628;
      continue;
    } else {
      moonbit_decref(_M0L3srcS938);
    }
    break;
  }
  _M0L6_2atmpS2629 = _M0L3bufS941;
  _block_3299
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3299)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3299->$0 = _M0L6_2atmpS2629;
  return _block_3299;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq15complement__rna(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS937
) {
  moonbit_string_t _M0L4dataS2624;
  moonbit_string_t _M0L6_2atmpS2623;
  struct _M0TP26biolab8bio__seq3Seq* _block_3300;
  #line 313 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2624 = _M0L4selfS937->$0;
  moonbit_incref(_M0L4dataS2624);
  #line 314 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2623
  = _M0FP26biolab8bio__seq17complement__chars(_M0L4dataS2624, _M0FP26biolab8bio__seq20rna__complement__arr);
  moonbit_decref(_M0L4dataS2624);
  _block_3300
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3300)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3300->$0 = _M0L6_2atmpS2623;
  return _block_3300;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq10complement(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS935
) {
  moonbit_string_t _M0L4dataS2622;
  moonbit_string_t _M0L6_2atmpS2621;
  struct _M0TP26biolab8bio__seq3Seq* _block_3301;
  #line 306 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2622 = _M0L4selfS935->$0;
  moonbit_incref(_M0L4dataS2622);
  #line 307 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2621
  = _M0FP26biolab8bio__seq17complement__chars(_M0L4dataS2622, _M0FP26biolab8bio__seq20dna__complement__arr);
  moonbit_decref(_M0L4dataS2622);
  _block_3301
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3301)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3301->$0 = _M0L6_2atmpS2621;
  return _block_3301;
}

moonbit_string_t _M0FP26biolab8bio__seq17complement__chars(
  moonbit_string_t _M0L1sS928,
  uint16_t* _M0L5tableS933
) {
  int32_t _M0L1nS927;
  uint16_t* _M0L3bufS929;
  int32_t _M0L1iS930;
  #line 115 "/home/developer/Documents/1/complement.mbt"
  _M0L1nS927 = Moonbit_array_length(_M0L1sS928);
  _M0L3bufS929 = (uint16_t*)moonbit_make_string(_M0L1nS927, 0);
  _M0L1iS930 = 0;
  while (1) {
    if (_M0L1iS930 < _M0L1nS927) {
      int32_t _M0L6_2atmpS2619 = _M0L1sS928[_M0L1iS930];
      int32_t _M0L2ciS931 = (int32_t)_M0L6_2atmpS2619;
      int32_t _M0L6_2atmpS2620;
      if (_M0L2ciS931 < 128) {
        int32_t _M0L4compS932;
        if (
          _M0L2ciS931 < 0
          || _M0L2ciS931 >= Moonbit_array_length(_M0L5tableS933)
        ) {
          #line 121 "/home/developer/Documents/1/complement.mbt"
          moonbit_panic();
        }
        _M0L4compS932 = (int32_t)_M0L5tableS933[_M0L2ciS931];
        if (_M0L4compS932 != 0) {
          if (
            _M0L1iS930 < 0
            || _M0L1iS930 >= Moonbit_array_length(_M0L3bufS929)
          ) {
            #line 123 "/home/developer/Documents/1/complement.mbt"
            moonbit_panic();
          }
          _M0L3bufS929[_M0L1iS930] = _M0L4compS932;
        } else {
          int32_t _M0L6_2atmpS2617 = _M0L1sS928[_M0L1iS930];
          if (
            _M0L1iS930 < 0
            || _M0L1iS930 >= Moonbit_array_length(_M0L3bufS929)
          ) {
            #line 125 "/home/developer/Documents/1/complement.mbt"
            moonbit_panic();
          }
          _M0L3bufS929[_M0L1iS930] = _M0L6_2atmpS2617;
        }
      } else {
        int32_t _M0L6_2atmpS2618 = _M0L1sS928[_M0L1iS930];
        if (
          _M0L1iS930 < 0 || _M0L1iS930 >= Moonbit_array_length(_M0L3bufS929)
        ) {
          #line 128 "/home/developer/Documents/1/complement.mbt"
          moonbit_panic();
        }
        _M0L3bufS929[_M0L1iS930] = _M0L6_2atmpS2618;
      }
      _M0L6_2atmpS2620 = _M0L1iS930 + 1;
      _M0L1iS930 = _M0L6_2atmpS2620;
      continue;
    }
    break;
  }
  return _M0L3bufS929;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq4join(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS925,
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L6othersS922
) {
  struct _M0TPB13StringBuilder* _M0L3bufS920;
  int32_t _M0L7_2abindS921;
  int32_t _M0L1iS923;
  moonbit_string_t _M0L6_2atmpS2616;
  struct _M0TP26biolab8bio__seq3Seq* _block_3304;
  #line 288 "/home/developer/Documents/1/seq.mbt"
  #line 289 "/home/developer/Documents/1/seq.mbt"
  _M0L3bufS920 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS921 = _M0L6othersS922->$1;
  _M0L1iS923 = 0;
  while (1) {
    if (_M0L1iS923 < _M0L7_2abindS921) {
      struct _M0TP26biolab8bio__seq3Seq** _M0L3bufS2615 = _M0L6othersS922->$0;
      struct _M0TP26biolab8bio__seq3Seq* _M0L1sS924 =
        (struct _M0TP26biolab8bio__seq3Seq*)_M0L3bufS2615[_M0L1iS923];
      moonbit_string_t _M0L8_2afieldS3069;
      int32_t _M0L6_2acntS3222;
      moonbit_string_t _M0L4dataS2613;
      int32_t _M0L6_2atmpS2614;
      if (_M0L1iS923 > 0) {
        moonbit_string_t _M0L4dataS2612 = _M0L4selfS925->$0;
        moonbit_incref(_M0L4dataS2612);
        moonbit_incref(_M0L1sS924);
        #line 292 "/home/developer/Documents/1/seq.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS920, _M0L4dataS2612);
        moonbit_decref(_M0L4dataS2612);
      } else {
        moonbit_incref(_M0L1sS924);
      }
      _M0L8_2afieldS3069 = _M0L1sS924->$0;
      _M0L6_2acntS3222 = Moonbit_rc_count(Moonbit_object_header(_M0L1sS924));
      if (_M0L6_2acntS3222 > 1) {
        int32_t _M0L11_2anew__cntS3223 = _M0L6_2acntS3222 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L1sS924), _M0L11_2anew__cntS3223);
        moonbit_incref(_M0L8_2afieldS3069);
      } else if (_M0L6_2acntS3222 == 1) {
        #line 294 "/home/developer/Documents/1/seq.mbt"
        moonbit_free(_M0L1sS924);
      }
      _M0L4dataS2613 = _M0L8_2afieldS3069;
      #line 294 "/home/developer/Documents/1/seq.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS920, _M0L4dataS2613);
      moonbit_decref(_M0L4dataS2613);
      _M0L6_2atmpS2614 = _M0L1iS923 + 1;
      _M0L1iS923 = _M0L6_2atmpS2614;
      continue;
    }
    break;
  }
  #line 296 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2616 = _M0MPB13StringBuilder10to__string(_M0L3bufS920);
  moonbit_decref(_M0L3bufS920);
  _block_3304
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3304)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3304->$0 = _M0L6_2atmpS2616;
  return _block_3304;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0MP26biolab8bio__seq3Seq5split(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS912,
  moonbit_string_t _M0L3sepS913
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L3outS910;
  moonbit_string_t _M0L4dataS2609;
  int32_t _M0L6_2atmpS2611;
  struct _M0TPC16string10StringView _M0L6_2atmpS2610;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L5_2aitS911;
  #line 277 "/home/developer/Documents/1/seq.mbt"
  #line 278 "/home/developer/Documents/1/seq.mbt"
  _M0L3outS910 = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq3SeqE(0);
  _M0L4dataS2609 = _M0L4selfS912->$0;
  _M0L6_2atmpS2611 = Moonbit_array_length(_M0L3sepS913);
  moonbit_incref(_M0L3sepS913);
  _M0L6_2atmpS2610
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L3sepS913, .$1 = 0, .$2 = _M0L6_2atmpS2611
  };
  moonbit_incref(_M0L4dataS2609);
  #line 279 "/home/developer/Documents/1/seq.mbt"
  _M0L5_2aitS911
  = _M0MPC16string6String5split(_M0L4dataS2609, _M0L6_2atmpS2610);
  moonbit_decref(_M0L4dataS2609);
  moonbit_decref(_M0L6_2atmpS2610.$0);
  while (1) {
    struct _M0TPC16string10StringView _M0L4partS915;
    void* _M0L7_2abindS917;
    moonbit_string_t _M0L6_2atmpS2608;
    struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2607;
    #line 279 "/home/developer/Documents/1/seq.mbt"
    _M0L7_2abindS917
    = _M0MPB4Iter4nextGRPC16string10StringViewE(_M0L5_2aitS911);
    switch (Moonbit_object_tag(_M0L7_2abindS917)) {
      case 1: {
        struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS918 =
          (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS917;
        struct _M0TPC16string10StringView _M0L8_2afieldS3073 =
          _M0L7_2aSomeS918->$0;
        int32_t _M0L6_2acntS3224 =
          Moonbit_rc_count(Moonbit_object_header(_M0L7_2aSomeS918));
        struct _M0TPC16string10StringView _M0L7_2apartS919;
        if (_M0L6_2acntS3224 > 1) {
          int32_t _M0L11_2anew__cntS3225 = _M0L6_2acntS3224 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2aSomeS918), _M0L11_2anew__cntS3225);
          moonbit_incref(_M0L8_2afieldS3073.$0);
        } else if (_M0L6_2acntS3224 == 1) {
          #line 279 "/home/developer/Documents/1/seq.mbt"
          moonbit_free(_M0L7_2aSomeS918);
        }
        _M0L7_2apartS919 = _M0L8_2afieldS3073;
        _M0L4partS915 = _M0L7_2apartS919;
        goto join_914;
        break;
      }
      default: {
        moonbit_decref(_M0L7_2abindS917);
        moonbit_decref(_M0L5_2aitS911);
        break;
      }
    }
    goto joinlet_3306;
    join_914:;
    #line 280 "/home/developer/Documents/1/seq.mbt"
    _M0L6_2atmpS2608 = _M0MPC16string10StringView9to__owned(_M0L4partS915);
    moonbit_decref(_M0L4partS915.$0);
    _M0L6_2atmpS2607
    = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
    Moonbit_object_header(_M0L6_2atmpS2607)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
    _M0L6_2atmpS2607->$0 = _M0L6_2atmpS2608;
    #line 280 "/home/developer/Documents/1/seq.mbt"
    _M0MPC15array5Array4pushGRP26biolab8bio__seq3SeqE(_M0L3outS910, _M0L6_2atmpS2607);
    moonbit_decref(_M0L6_2atmpS2607);
    continue;
    joinlet_3306:;
    break;
  }
  return _M0L3outS910;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq7replace(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS907,
  moonbit_string_t _M0L3oldS908,
  moonbit_string_t _M0L3newS909
) {
  moonbit_string_t _M0L4dataS2602;
  int32_t _M0L6_2atmpS2606;
  struct _M0TPC16string10StringView _M0L6_2atmpS2603;
  int32_t _M0L6_2atmpS2605;
  struct _M0TPC16string10StringView _M0L6_2atmpS2604;
  moonbit_string_t _M0L6_2atmpS2601;
  struct _M0TP26biolab8bio__seq3Seq* _block_3307;
  #line 271 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2602 = _M0L4selfS907->$0;
  _M0L6_2atmpS2606 = Moonbit_array_length(_M0L3oldS908);
  moonbit_incref(_M0L3oldS908);
  _M0L6_2atmpS2603
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L3oldS908, .$1 = 0, .$2 = _M0L6_2atmpS2606
  };
  _M0L6_2atmpS2605 = Moonbit_array_length(_M0L3newS909);
  moonbit_incref(_M0L3newS909);
  _M0L6_2atmpS2604
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L3newS909, .$1 = 0, .$2 = _M0L6_2atmpS2605
  };
  moonbit_incref(_M0L4dataS2602);
  #line 272 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2601
  = _M0MPC16string6String12replace__all(_M0L4dataS2602, _M0L6_2atmpS2603, _M0L6_2atmpS2604);
  moonbit_decref(_M0L4dataS2602);
  moonbit_decref(_M0L6_2atmpS2603.$0);
  moonbit_decref(_M0L6_2atmpS2604.$0);
  _block_3307
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3307)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3307->$0 = _M0L6_2atmpS2601;
  return _block_3307;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq14rstrip_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS905,
  moonbit_string_t _M0L5charsS906
) {
  moonbit_string_t _M0L4dataS2598;
  int32_t _M0L6_2atmpS2600;
  struct _M0TPC16string10StringView _M0L6_2atmpS2599;
  struct _M0TPC16string10StringView _M0L6_2atmpS2597;
  moonbit_string_t _M0L6_2atmpS2596;
  struct _M0TP26biolab8bio__seq3Seq* _block_3308;
  #line 265 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2598 = _M0L4selfS905->$0;
  _M0L6_2atmpS2600 = Moonbit_array_length(_M0L5charsS906);
  moonbit_incref(_M0L5charsS906);
  _M0L6_2atmpS2599
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L5charsS906, .$1 = 0, .$2 = _M0L6_2atmpS2600
  };
  moonbit_incref(_M0L4dataS2598);
  #line 266 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2597
  = _M0MPC16string6String17trim__end_2einner(_M0L4dataS2598, _M0L6_2atmpS2599);
  moonbit_decref(_M0L4dataS2598);
  moonbit_decref(_M0L6_2atmpS2599.$0);
  #line 266 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2596 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2597);
  moonbit_decref(_M0L6_2atmpS2597.$0);
  _block_3308
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3308)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3308->$0 = _M0L6_2atmpS2596;
  return _block_3308;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq14lstrip_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS903,
  moonbit_string_t _M0L5charsS904
) {
  moonbit_string_t _M0L4dataS2593;
  int32_t _M0L6_2atmpS2595;
  struct _M0TPC16string10StringView _M0L6_2atmpS2594;
  struct _M0TPC16string10StringView _M0L6_2atmpS2592;
  moonbit_string_t _M0L6_2atmpS2591;
  struct _M0TP26biolab8bio__seq3Seq* _block_3309;
  #line 259 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2593 = _M0L4selfS903->$0;
  _M0L6_2atmpS2595 = Moonbit_array_length(_M0L5charsS904);
  moonbit_incref(_M0L5charsS904);
  _M0L6_2atmpS2594
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L5charsS904, .$1 = 0, .$2 = _M0L6_2atmpS2595
  };
  moonbit_incref(_M0L4dataS2593);
  #line 260 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2592
  = _M0MPC16string6String19trim__start_2einner(_M0L4dataS2593, _M0L6_2atmpS2594);
  moonbit_decref(_M0L4dataS2593);
  moonbit_decref(_M0L6_2atmpS2594.$0);
  #line 260 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2591 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2592);
  moonbit_decref(_M0L6_2atmpS2592.$0);
  _block_3309
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3309)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3309->$0 = _M0L6_2atmpS2591;
  return _block_3309;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq13strip_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS901,
  moonbit_string_t _M0L5charsS902
) {
  moonbit_string_t _M0L4dataS2588;
  int32_t _M0L6_2atmpS2590;
  struct _M0TPC16string10StringView _M0L6_2atmpS2589;
  struct _M0TPC16string10StringView _M0L6_2atmpS2587;
  moonbit_string_t _M0L6_2atmpS2586;
  struct _M0TP26biolab8bio__seq3Seq* _block_3310;
  #line 253 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2588 = _M0L4selfS901->$0;
  _M0L6_2atmpS2590 = Moonbit_array_length(_M0L5charsS902);
  moonbit_incref(_M0L5charsS902);
  _M0L6_2atmpS2589
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L5charsS902, .$1 = 0, .$2 = _M0L6_2atmpS2590
  };
  moonbit_incref(_M0L4dataS2588);
  #line 254 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2587
  = _M0MPC16string6String12trim_2einner(_M0L4dataS2588, _M0L6_2atmpS2589);
  moonbit_decref(_M0L4dataS2588);
  moonbit_decref(_M0L6_2atmpS2589.$0);
  #line 254 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2586 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2587);
  moonbit_decref(_M0L6_2atmpS2587.$0);
  _block_3310
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3310)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3310->$0 = _M0L6_2atmpS2586;
  return _block_3310;
}

int32_t _M0MP26biolab8bio__seq3Seq8contains(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS899,
  moonbit_string_t _M0L3subS900
) {
  moonbit_string_t _M0L4dataS2583;
  int32_t _M0L6_2atmpS2585;
  struct _M0TPC16string10StringView _M0L6_2atmpS2584;
  int32_t _result_3311;
  #line 246 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2583 = _M0L4selfS899->$0;
  _M0L6_2atmpS2585 = Moonbit_array_length(_M0L3subS900);
  moonbit_incref(_M0L3subS900);
  _M0L6_2atmpS2584
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L3subS900, .$1 = 0, .$2 = _M0L6_2atmpS2585
  };
  moonbit_incref(_M0L4dataS2583);
  #line 247 "/home/developer/Documents/1/seq.mbt"
  _result_3311
  = _M0MPC16string6String8contains(_M0L4dataS2583, _M0L6_2atmpS2584);
  moonbit_decref(_M0L4dataS2583);
  moonbit_decref(_M0L6_2atmpS2584.$0);
  return _result_3311;
}

int32_t _M0MP26biolab8bio__seq3Seq8endswith(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS894,
  moonbit_string_t _M0L6suffixS892
) {
  int32_t _M0L4slenS891;
  moonbit_string_t _M0L4dataS893;
  int32_t _M0L4dlenS895;
  int32_t _M0L6offsetS896;
  int32_t _M0L1iS897;
  #line 228 "/home/developer/Documents/1/seq.mbt"
  _M0L4slenS891 = Moonbit_array_length(_M0L6suffixS892);
  _M0L4dataS893 = _M0L4selfS894->$0;
  _M0L4dlenS895 = Moonbit_array_length(_M0L4dataS893);
  if (_M0L4dlenS895 < _M0L4slenS891) {
    return 0;
  } else {
    moonbit_incref(_M0L4dataS893);
  }
  _M0L6offsetS896 = _M0L4dlenS895 - _M0L4slenS891;
  _M0L1iS897 = 0;
  while (1) {
    if (_M0L1iS897 < _M0L4slenS891) {
      int32_t _M0L6_2atmpS2581 = _M0L6offsetS896 + _M0L1iS897;
      int32_t _M0L6_2atmpS2579 = _M0L4dataS893[_M0L6_2atmpS2581];
      int32_t _M0L6_2atmpS2580 = _M0L6suffixS892[_M0L1iS897];
      int32_t _M0L6_2atmpS2582;
      if (_M0L6_2atmpS2579 != _M0L6_2atmpS2580) {
        moonbit_decref(_M0L4dataS893);
        return 0;
      }
      _M0L6_2atmpS2582 = _M0L1iS897 + 1;
      _M0L1iS897 = _M0L6_2atmpS2582;
      continue;
    } else {
      moonbit_decref(_M0L4dataS893);
    }
    break;
  }
  return 1;
}

int32_t _M0MP26biolab8bio__seq3Seq10startswith(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS888,
  moonbit_string_t _M0L6prefixS886
) {
  int32_t _M0L4plenS885;
  moonbit_string_t _M0L4dataS887;
  int32_t _M0L6_2atmpS2575;
  int32_t _M0L1iS889;
  #line 210 "/home/developer/Documents/1/seq.mbt"
  _M0L4plenS885 = Moonbit_array_length(_M0L6prefixS886);
  _M0L4dataS887 = _M0L4selfS888->$0;
  _M0L6_2atmpS2575 = Moonbit_array_length(_M0L4dataS887);
  if (_M0L6_2atmpS2575 < _M0L4plenS885) {
    return 0;
  } else {
    moonbit_incref(_M0L4dataS887);
  }
  _M0L1iS889 = 0;
  while (1) {
    if (_M0L1iS889 < _M0L4plenS885) {
      int32_t _M0L6_2atmpS2576 = _M0L4dataS887[_M0L1iS889];
      int32_t _M0L6_2atmpS2577 = _M0L6prefixS886[_M0L1iS889];
      int32_t _M0L6_2atmpS2578;
      if (_M0L6_2atmpS2576 != _M0L6_2atmpS2577) {
        moonbit_decref(_M0L4dataS887);
        return 0;
      }
      _M0L6_2atmpS2578 = _M0L1iS889 + 1;
      _M0L1iS889 = _M0L6_2atmpS2578;
      continue;
    } else {
      moonbit_decref(_M0L4dataS887);
    }
    break;
  }
  return 1;
}

int32_t _M0MP26biolab8bio__seq3Seq5rfind(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS881,
  moonbit_string_t _M0L3subS882
) {
  int32_t _M0L1pS879;
  moonbit_string_t _M0L4dataS2572;
  int32_t _M0L6_2atmpS2574;
  struct _M0TPC16string10StringView _M0L6_2atmpS2573;
  int64_t _M0L7_2abindS880;
  #line 199 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2572 = _M0L4selfS881->$0;
  _M0L6_2atmpS2574 = Moonbit_array_length(_M0L3subS882);
  moonbit_incref(_M0L3subS882);
  _M0L6_2atmpS2573
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L3subS882, .$1 = 0, .$2 = _M0L6_2atmpS2574
  };
  moonbit_incref(_M0L4dataS2572);
  #line 200 "/home/developer/Documents/1/seq.mbt"
  _M0L7_2abindS880
  = _M0MPC16string6String9rev__find(_M0L4dataS2572, _M0L6_2atmpS2573);
  moonbit_decref(_M0L4dataS2572);
  moonbit_decref(_M0L6_2atmpS2573.$0);
  if (_M0L7_2abindS880 == 4294967296ll) {
    return -1;
  } else {
    int64_t _M0L7_2aSomeS883 = _M0L7_2abindS880;
    int32_t _M0L4_2apS884 = (int32_t)_M0L7_2aSomeS883;
    _M0L1pS879 = _M0L4_2apS884;
    goto join_878;
  }
  join_878:;
  return _M0L1pS879;
}

int32_t _M0MP26biolab8bio__seq3Seq4find(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS874,
  moonbit_string_t _M0L3subS875
) {
  int32_t _M0L1pS872;
  moonbit_string_t _M0L4dataS2569;
  int32_t _M0L6_2atmpS2571;
  struct _M0TPC16string10StringView _M0L6_2atmpS2570;
  int64_t _M0L7_2abindS873;
  #line 190 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2569 = _M0L4selfS874->$0;
  _M0L6_2atmpS2571 = Moonbit_array_length(_M0L3subS875);
  moonbit_incref(_M0L3subS875);
  _M0L6_2atmpS2570
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L3subS875, .$1 = 0, .$2 = _M0L6_2atmpS2571
  };
  moonbit_incref(_M0L4dataS2569);
  #line 191 "/home/developer/Documents/1/seq.mbt"
  _M0L7_2abindS873
  = _M0MPC16string6String4find(_M0L4dataS2569, _M0L6_2atmpS2570);
  moonbit_decref(_M0L4dataS2569);
  moonbit_decref(_M0L6_2atmpS2570.$0);
  if (_M0L7_2abindS873 == 4294967296ll) {
    return -1;
  } else {
    int64_t _M0L7_2aSomeS876 = _M0L7_2abindS873;
    int32_t _M0L4_2apS877 = (int32_t)_M0L7_2aSomeS876;
    _M0L1pS872 = _M0L4_2apS877;
    goto join_871;
  }
  join_871:;
  return _M0L1pS872;
}

int32_t _M0MP26biolab8bio__seq3Seq14count__overlap(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS860,
  moonbit_string_t _M0L3subS859
) {
  int32_t _M0L8sub__lenS858;
  struct _M0TPB8MutLocalGiE* _M0L5totalS861;
  moonbit_string_t _M0L4dataS862;
  struct _M0TPB8MutLocalGiE* _M0L5startS863;
  int32_t _M0L1nS864;
  int32_t _result_3318;
  #line 166 "/home/developer/Documents/1/seq.mbt"
  _M0L8sub__lenS858 = Moonbit_array_length(_M0L3subS859);
  if (_M0L8sub__lenS858 == 0) {
    moonbit_string_t _M0L8_2afieldS3085 = _M0L4selfS860->$0;
    moonbit_string_t _M0L4dataS2563;
    int32_t _M0L6_2atmpS2562;
    moonbit_incref(_M0L8_2afieldS3085);
    _M0L4dataS2563 = _M0L8_2afieldS3085;
    _M0L6_2atmpS2562 = Moonbit_array_length(_M0L4dataS2563);
    moonbit_decref(_M0L4dataS2563);
    return _M0L6_2atmpS2562 + 1;
  }
  _M0L5totalS861
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5totalS861)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5totalS861->$0 = 0;
  _M0L4dataS862 = _M0L4selfS860->$0;
  _M0L5startS863
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS863)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS863->$0 = 0;
  _M0L1nS864 = Moonbit_array_length(_M0L4dataS862);
  moonbit_incref(_M0L4dataS862);
  while (1) {
    int32_t _M0L3valS2564 = _M0L5startS863->$0;
    if (_M0L3valS2564 <= _M0L1nS864) {
      int32_t _M0L3posS866;
      int32_t _M0L3valS2568 = _M0L5startS863->$0;
      int64_t _M0L7_2abindS867;
      int32_t _M0L3valS2566;
      int32_t _M0L6_2atmpS2565;
      int32_t _M0L6_2atmpS2567;
      #line 176 "/home/developer/Documents/1/seq.mbt"
      _M0L7_2abindS867
      = _M0FP26biolab8bio__seq10find__from(_M0L4dataS862, _M0L3subS859, _M0L3valS2568);
      if (_M0L7_2abindS867 == 4294967296ll) {
        moonbit_decref(_M0L5startS863);
        moonbit_decref(_M0L4dataS862);
        break;
      } else {
        int64_t _M0L7_2aSomeS868 = _M0L7_2abindS867;
        int32_t _M0L6_2aposS869 = (int32_t)_M0L7_2aSomeS868;
        _M0L3posS866 = _M0L6_2aposS869;
        goto join_865;
      }
      goto joinlet_3317;
      join_865:;
      _M0L3valS2566 = _M0L5totalS861->$0;
      _M0L6_2atmpS2565 = _M0L3valS2566 + 1;
      _M0L5totalS861->$0 = _M0L6_2atmpS2565;
      _M0L6_2atmpS2567 = _M0L3posS866 + 1;
      _M0L5startS863->$0 = _M0L6_2atmpS2567;
      joinlet_3317:;
      continue;
    } else {
      moonbit_decref(_M0L5startS863);
      moonbit_decref(_M0L4dataS862);
    }
    break;
  }
  _result_3318 = _M0L5totalS861->$0;
  moonbit_decref(_M0L5totalS861);
  return _result_3318;
}

int32_t _M0MP26biolab8bio__seq3Seq5count(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS847,
  moonbit_string_t _M0L3subS846
) {
  int32_t _M0L8sub__lenS845;
  struct _M0TPB8MutLocalGiE* _M0L5totalS848;
  moonbit_string_t _M0L4dataS849;
  struct _M0TPB8MutLocalGiE* _M0L5startS850;
  int32_t _M0L1nS851;
  int32_t _result_3321;
  #line 142 "/home/developer/Documents/1/seq.mbt"
  _M0L8sub__lenS845 = Moonbit_array_length(_M0L3subS846);
  if (_M0L8sub__lenS845 == 0) {
    moonbit_string_t _M0L8_2afieldS3087 = _M0L4selfS847->$0;
    moonbit_string_t _M0L4dataS2556;
    int32_t _M0L6_2atmpS2555;
    moonbit_incref(_M0L8_2afieldS3087);
    _M0L4dataS2556 = _M0L8_2afieldS3087;
    _M0L6_2atmpS2555 = Moonbit_array_length(_M0L4dataS2556);
    moonbit_decref(_M0L4dataS2556);
    return _M0L6_2atmpS2555 + 1;
  }
  _M0L5totalS848
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5totalS848)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5totalS848->$0 = 0;
  _M0L4dataS849 = _M0L4selfS847->$0;
  _M0L5startS850
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS850)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS850->$0 = 0;
  _M0L1nS851 = Moonbit_array_length(_M0L4dataS849);
  moonbit_incref(_M0L4dataS849);
  while (1) {
    int32_t _M0L3valS2557 = _M0L5startS850->$0;
    if (_M0L3valS2557 <= _M0L1nS851) {
      int32_t _M0L3posS853;
      int32_t _M0L3valS2561 = _M0L5startS850->$0;
      int64_t _M0L7_2abindS854;
      int32_t _M0L3valS2559;
      int32_t _M0L6_2atmpS2558;
      int32_t _M0L6_2atmpS2560;
      #line 152 "/home/developer/Documents/1/seq.mbt"
      _M0L7_2abindS854
      = _M0FP26biolab8bio__seq10find__from(_M0L4dataS849, _M0L3subS846, _M0L3valS2561);
      if (_M0L7_2abindS854 == 4294967296ll) {
        moonbit_decref(_M0L5startS850);
        moonbit_decref(_M0L4dataS849);
        break;
      } else {
        int64_t _M0L7_2aSomeS855 = _M0L7_2abindS854;
        int32_t _M0L6_2aposS856 = (int32_t)_M0L7_2aSomeS855;
        _M0L3posS853 = _M0L6_2aposS856;
        goto join_852;
      }
      goto joinlet_3320;
      join_852:;
      _M0L3valS2559 = _M0L5totalS848->$0;
      _M0L6_2atmpS2558 = _M0L3valS2559 + 1;
      _M0L5totalS848->$0 = _M0L6_2atmpS2558;
      _M0L6_2atmpS2560 = _M0L3posS853 + _M0L8sub__lenS845;
      _M0L5startS850->$0 = _M0L6_2atmpS2560;
      joinlet_3320:;
      continue;
    } else {
      moonbit_decref(_M0L5startS850);
      moonbit_decref(_M0L4dataS849);
    }
    break;
  }
  _result_3321 = _M0L5totalS848->$0;
  moonbit_decref(_M0L5totalS848);
  return _result_3321;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq5lower(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS844
) {
  moonbit_string_t _M0L4dataS2554;
  moonbit_string_t _M0L6_2atmpS2553;
  struct _M0TP26biolab8bio__seq3Seq* _block_3322;
  #line 135 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2554 = _M0L4selfS844->$0;
  moonbit_incref(_M0L4dataS2554);
  #line 136 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2553 = _M0MPC16string6String9to__lower(_M0L4dataS2554);
  moonbit_decref(_M0L4dataS2554);
  _block_3322
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3322)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3322->$0 = _M0L6_2atmpS2553;
  return _block_3322;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq5upper(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS843
) {
  moonbit_string_t _M0L4dataS2552;
  moonbit_string_t _M0L6_2atmpS2551;
  struct _M0TP26biolab8bio__seq3Seq* _block_3323;
  #line 129 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2552 = _M0L4selfS843->$0;
  moonbit_incref(_M0L4dataS2552);
  #line 130 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2551 = _M0MPC16string6String9to__upper(_M0L4dataS2552);
  moonbit_decref(_M0L4dataS2552);
  _block_3323
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3323)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3323->$0 = _M0L6_2atmpS2551;
  return _block_3323;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq5slice(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS840,
  int32_t _M0L5startS841,
  int32_t _M0L3endS842
) {
  moonbit_string_t _M0L4dataS2549;
  int64_t _M0L6_2atmpS2550;
  struct _M0TPC16string10StringView _M0L6_2atmpS2548;
  moonbit_string_t _M0L6_2atmpS2547;
  struct _M0TP26biolab8bio__seq3Seq* _block_3324;
  #line 117 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2549 = _M0L4selfS840->$0;
  _M0L6_2atmpS2550 = (int64_t)_M0L3endS842;
  moonbit_incref(_M0L4dataS2549);
  #line 118 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2548
  = _M0MPC16string6String11sub_2einner(_M0L4dataS2549, _M0L5startS841, _M0L6_2atmpS2550);
  moonbit_decref(_M0L4dataS2549);
  #line 118 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2547 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2548);
  moonbit_decref(_M0L6_2atmpS2548.$0);
  _block_3324
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3324)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3324->$0 = _M0L6_2atmpS2547;
  return _block_3324;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3mul(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS839,
  int32_t _M0L1nS838
) {
  #line 106 "/home/developer/Documents/1/seq.mbt"
  if (_M0L1nS838 <= 0) {
    struct _M0TP26biolab8bio__seq3Seq* _block_3325 =
      (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
    Moonbit_object_header(_block_3325)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
    _block_3325->$0 = (moonbit_string_t)moonbit_string_literal_0.data;
    return _block_3325;
  } else {
    moonbit_string_t _M0L4dataS2546 = _M0L4selfS839->$0;
    moonbit_string_t _M0L6_2atmpS2545;
    struct _M0TP26biolab8bio__seq3Seq* _block_3326;
    moonbit_incref(_M0L4dataS2546);
    #line 110 "/home/developer/Documents/1/seq.mbt"
    _M0L6_2atmpS2545
    = _M0MPC16string6String6repeat(_M0L4dataS2546, _M0L1nS838);
    moonbit_decref(_M0L4dataS2546);
    _block_3326
    = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
    Moonbit_object_header(_block_3326)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
    _block_3326->$0 = _M0L6_2atmpS2545;
    return _block_3326;
  }
}

struct _M0TP26biolab8bio__seq3Seq* _M0IP26biolab8bio__seq3SeqPB3Add3add(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS836,
  struct _M0TP26biolab8bio__seq3Seq* _M0L5otherS837
) {
  moonbit_string_t _M0L4dataS2543;
  moonbit_string_t _M0L4dataS2544;
  moonbit_string_t _M0L6_2atmpS2542;
  struct _M0TP26biolab8bio__seq3Seq* _block_3327;
  #line 93 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2543 = _M0L4selfS836->$0;
  _M0L4dataS2544 = _M0L5otherS837->$0;
  #line 94 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2542 = moonbit_add_string(_M0L4dataS2543, _M0L4dataS2544);
  _block_3327
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3327)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3327->$0 = _M0L6_2atmpS2542;
  return _block_3327;
}

int64_t _M0FP26biolab8bio__seq10find__from(
  moonbit_string_t _M0L1sS825,
  moonbit_string_t _M0L3subS827,
  int32_t _M0L4fromS828
) {
  int32_t _M0L1nS824;
  int32_t _M0L1mS826;
  int32_t _if__result_3328;
  int32_t _M0L5firstS829;
  struct _M0TPB8MutLocalGiE* _M0L1iS830;
  #line 23 "/home/developer/Documents/1/seq.mbt"
  _M0L1nS824 = Moonbit_array_length(_M0L1sS825);
  _M0L1mS826 = Moonbit_array_length(_M0L3subS827);
  if (_M0L1mS826 == 0) {
    return (int64_t)_M0L4fromS828;
  }
  if (_M0L4fromS828 < 0) {
    _if__result_3328 = 1;
  } else {
    int32_t _M0L6_2atmpS2527 = _M0L4fromS828 + _M0L1mS826;
    _if__result_3328 = _M0L6_2atmpS2527 > _M0L1nS824;
  }
  if (_if__result_3328) {
    return 4294967296ll;
  }
  _M0L5firstS829 = _M0L3subS827[0];
  _M0L1iS830
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS830)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS830->$0 = _M0L4fromS828;
  while (1) {
    int32_t _M0L3valS2529 = _M0L1iS830->$0;
    int32_t _M0L6_2atmpS2528 = _M0L3valS2529 + _M0L1mS826;
    if (_M0L6_2atmpS2528 <= _M0L1nS824) {
      int32_t _M0L3valS2531 = _M0L1iS830->$0;
      int32_t _M0L6_2atmpS2530 = _M0L1sS825[_M0L3valS2531];
      struct _M0TPB8MutLocalGbE* _M0L7matchedS832;
      int32_t _M0L7_2abindS833;
      int32_t _M0L1jS834;
      int32_t _result_3331;
      int32_t _M0L3valS2541;
      int32_t _M0L6_2atmpS2540;
      if (_M0L6_2atmpS2530 != _M0L5firstS829) {
        int32_t _M0L3valS2533 = _M0L1iS830->$0;
        int32_t _M0L6_2atmpS2532 = _M0L3valS2533 + 1;
        _M0L1iS830->$0 = _M0L6_2atmpS2532;
        continue;
      }
      _M0L7matchedS832
      = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
      Moonbit_object_header(_M0L7matchedS832)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L7matchedS832->$0 = 1;
      _M0L7_2abindS833 = 1;
      _M0L1jS834 = _M0L7_2abindS833;
      while (1) {
        if (_M0L1jS834 < _M0L1mS826) {
          int32_t _M0L3valS2537 = _M0L1iS830->$0;
          int32_t _M0L6_2atmpS2536 = _M0L3valS2537 + _M0L1jS834;
          int32_t _M0L6_2atmpS2534 = _M0L1sS825[_M0L6_2atmpS2536];
          int32_t _M0L6_2atmpS2535 = _M0L3subS827[_M0L1jS834];
          int32_t _M0L6_2atmpS2538;
          if (_M0L6_2atmpS2534 != _M0L6_2atmpS2535) {
            _M0L7matchedS832->$0 = 0;
            break;
          }
          _M0L6_2atmpS2538 = _M0L1jS834 + 1;
          _M0L1jS834 = _M0L6_2atmpS2538;
          continue;
        }
        break;
      }
      _result_3331 = _M0L7matchedS832->$0;
      moonbit_decref(_M0L7matchedS832);
      if (_result_3331) {
        int32_t _M0L3valS2539 = _M0L1iS830->$0;
        moonbit_decref(_M0L1iS830);
        return (int64_t)_M0L3valS2539;
      }
      _M0L3valS2541 = _M0L1iS830->$0;
      _M0L6_2atmpS2540 = _M0L3valS2541 + 1;
      _M0L1iS830->$0 = _M0L6_2atmpS2540;
      continue;
    } else {
      moonbit_decref(_M0L1iS830);
    }
    break;
  }
  return 4294967296ll;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq19reverse__complement(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS823
) {
  moonbit_string_t _M0L4dataS2526;
  moonbit_string_t _M0L6_2atmpS2525;
  struct _M0TP26biolab8bio__seq3Seq* _block_3332;
  #line 334 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2526 = _M0L4selfS823->$0;
  moonbit_incref(_M0L4dataS2526);
  #line 335 "/home/developer/Documents/1/seq.mbt"
  _M0L6_2atmpS2525
  = _M0FP26biolab8bio__seq26reverse__complement__chars(_M0L4dataS2526, _M0FP26biolab8bio__seq20dna__complement__arr);
  moonbit_decref(_M0L4dataS2526);
  _block_3332
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3332)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3332->$0 = _M0L6_2atmpS2525;
  return _block_3332;
}

moonbit_string_t _M0FP26biolab8bio__seq26reverse__complement__chars(
  moonbit_string_t _M0L1sS815,
  uint16_t* _M0L5tableS821
) {
  int32_t _M0L1nS814;
  uint16_t* _M0L3bufS816;
  int32_t _M0L6_2atmpS2524;
  int32_t _M0L1iS817;
  #line 138 "/home/developer/Documents/1/complement.mbt"
  _M0L1nS814 = Moonbit_array_length(_M0L1sS815);
  _M0L3bufS816 = (uint16_t*)moonbit_make_string(_M0L1nS814, 0);
  _M0L6_2atmpS2524 = _M0L1nS814 - 1;
  _M0L1iS817 = _M0L6_2atmpS2524;
  while (1) {
    if (_M0L1iS817 >= 0) {
      int32_t _M0L6_2atmpS2522 = _M0L1sS815[_M0L1iS817];
      int32_t _M0L2ciS818 = (int32_t)_M0L6_2atmpS2522;
      int32_t _M0L6_2atmpS2521 = _M0L1nS814 - 1;
      int32_t _M0L8out__idxS819 = _M0L6_2atmpS2521 - _M0L1iS817;
      int32_t _M0L6_2atmpS2523;
      if (_M0L2ciS818 < 128) {
        int32_t _M0L4compS820;
        if (
          _M0L2ciS818 < 0
          || _M0L2ciS818 >= Moonbit_array_length(_M0L5tableS821)
        ) {
          #line 145 "/home/developer/Documents/1/complement.mbt"
          moonbit_panic();
        }
        _M0L4compS820 = (int32_t)_M0L5tableS821[_M0L2ciS818];
        if (_M0L4compS820 != 0) {
          if (
            _M0L8out__idxS819 < 0
            || _M0L8out__idxS819 >= Moonbit_array_length(_M0L3bufS816)
          ) {
            #line 147 "/home/developer/Documents/1/complement.mbt"
            moonbit_panic();
          }
          _M0L3bufS816[_M0L8out__idxS819] = _M0L4compS820;
        } else {
          int32_t _M0L6_2atmpS2519 = _M0L1sS815[_M0L1iS817];
          if (
            _M0L8out__idxS819 < 0
            || _M0L8out__idxS819 >= Moonbit_array_length(_M0L3bufS816)
          ) {
            #line 149 "/home/developer/Documents/1/complement.mbt"
            moonbit_panic();
          }
          _M0L3bufS816[_M0L8out__idxS819] = _M0L6_2atmpS2519;
        }
      } else {
        int32_t _M0L6_2atmpS2520 = _M0L1sS815[_M0L1iS817];
        if (
          _M0L8out__idxS819 < 0
          || _M0L8out__idxS819 >= Moonbit_array_length(_M0L3bufS816)
        ) {
          #line 152 "/home/developer/Documents/1/complement.mbt"
          moonbit_panic();
        }
        _M0L3bufS816[_M0L8out__idxS819] = _M0L6_2atmpS2520;
      }
      _M0L6_2atmpS2523 = _M0L1iS817 - 1;
      _M0L1iS817 = _M0L6_2atmpS2523;
      continue;
    }
    break;
  }
  return _M0L3bufS816;
}

int32_t _M0MP26biolab8bio__seq3Seq6length(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS812
) {
  moonbit_string_t _M0L4dataS2518;
  int32_t _result_3334;
  #line 86 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2518 = _M0L4selfS812->$0;
  moonbit_incref(_M0L4dataS2518);
  #line 87 "/home/developer/Documents/1/seq.mbt"
  _result_3334
  = _M0MPC16string6String20char__length_2einner(_M0L4dataS2518, 0, 4294967296ll);
  moonbit_decref(_M0L4dataS2518);
  return _result_3334;
}

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS811
) {
  moonbit_string_t _M0L8_2afieldS3096;
  #line 74 "/home/developer/Documents/1/seq.mbt"
  _M0L8_2afieldS3096 = _M0L4selfS811->$0;
  moonbit_incref(_M0L8_2afieldS3096);
  return _M0L8_2afieldS3096;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t _M0L1sS810
) {
  struct _M0TP26biolab8bio__seq3Seq* _block_3335;
  #line 62 "/home/developer/Documents/1/seq.mbt"
  moonbit_incref(_M0L1sS810);
  _block_3335
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3335)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 28, 0);
  _block_3335->$0 = _M0L1sS810;
  return _block_3335;
}

int32_t _M0IP26biolab8bio__seq3SeqPB2Eq5equal(
  struct _M0TP26biolab8bio__seq3Seq* _M0L10_2ax__1931S808,
  struct _M0TP26biolab8bio__seq3Seq* _M0L10_2ax__1932S809
) {
  moonbit_string_t _M0L4dataS2516;
  moonbit_string_t _M0L4dataS2517;
  #line 56 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2516 = _M0L10_2ax__1931S808->$0;
  _M0L4dataS2517 = _M0L10_2ax__1932S809->$0;
  #line 56 "/home/developer/Documents/1/seq.mbt"
  return _M0L4dataS2516 == _M0L4dataS2517
         || Moonbit_array_length(_M0L4dataS2516)
            == Moonbit_array_length(_M0L4dataS2517)
            && 0
               == memcmp(_M0L4dataS2516, _M0L4dataS2517, Moonbit_array_length(_M0L4dataS2516) * 2);
}

struct _M0TP26biolab8bio__seq3Seq* _M0MPC15array5Array2atGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L4selfS803,
  int32_t _M0L5indexS804
) {
  int32_t _M0L3lenS802;
  int32_t _if__result_3336;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS802 = _M0L4selfS803->$1;
  if (_M0L5indexS804 >= 0) {
    _if__result_3336 = _M0L5indexS804 < _M0L3lenS802;
  } else {
    _if__result_3336 = 0;
  }
  if (_if__result_3336) {
    struct _M0TP26biolab8bio__seq3Seq** _M0L6_2atmpS2514;
    struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS3097;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2514
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq3SeqE(_M0L4selfS803);
    _M0L6_2atmpS3097
    = (struct _M0TP26biolab8bio__seq3Seq*)_M0L6_2atmpS2514[_M0L5indexS804];
    if (_M0L6_2atmpS3097) {
      moonbit_incref(_M0L6_2atmpS3097);
    }
    moonbit_decref(_M0L6_2atmpS2514);
    return _M0L6_2atmpS3097;
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
  int32_t _if__result_3337;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS805 = _M0L4selfS806->$1;
  if (_M0L5indexS807 >= 0) {
    _if__result_3337 = _M0L5indexS807 < _M0L3lenS805;
  } else {
    _if__result_3337 = 0;
  }
  if (_if__result_3337) {
    moonbit_string_t* _M0L6_2atmpS2515;
    moonbit_string_t _M0L6_2atmpS3098;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2515 = _M0MPC15array5Array6bufferGsE(_M0L4selfS806);
    _M0L6_2atmpS3098 = (moonbit_string_t)_M0L6_2atmpS2515[_M0L5indexS807];
    moonbit_incref(_M0L6_2atmpS3098);
    moonbit_decref(_M0L6_2atmpS2515);
    return _M0L6_2atmpS3098;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB7printlnGsE(moonbit_string_t _M0L5inputS801) {
  moonbit_string_t _M0L6_2atmpS2513;
  #line 36 "/home/developer/.moon/lib/core/builtin/console.mbt"
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  _M0L6_2atmpS2513 = _M0IPC16string6StringPB4Show10to__string(_M0L5inputS801);
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  moonbit_println(_M0L6_2atmpS2513);
  moonbit_decref(_M0L6_2atmpS2513);
  return 0;
}

int32_t _M0IPC13int3IntPB4Hash4hash(int32_t _M0L4selfS800) {
  int32_t _tmp_3338;
  uint32_t _M0L6_2atmpS2512;
  uint32_t _M0L6_2atmpS2511;
  uint32_t _M0L3accS799;
  uint32_t _M0L6_2atmpS2510;
  uint32_t _M0L6_2atmpS2509;
  #line 588 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_3338 = 0;
  _M0L6_2atmpS2512 = *(uint32_t*)&_tmp_3338;
  _M0L6_2atmpS2511 = _M0L6_2atmpS2512 + 374761393u;
  _M0L3accS799 = _M0L6_2atmpS2511 + 4u;
  _M0L6_2atmpS2510 = *(uint32_t*)&_M0L4selfS800;
  #line 590 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS2509 = _M0FPB13consume4__acc(_M0L3accS799, _M0L6_2atmpS2510);
  #line 590 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS2509);
}

int32_t _M0MPB3Map3getGicE(
  struct _M0TPB3MapGicE* _M0L4selfS795,
  int32_t _M0L3keyS791
) {
  int32_t _M0L4hashS790;
  int32_t _M0L14capacity__maskS2508;
  int32_t _M0L6_2atmpS2507;
  int32_t _M0L1iS792;
  int32_t _M0L3idxS793;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS790 = _M0IPC13int3IntPB4Hash4hash(_M0L3keyS791);
  _M0L14capacity__maskS2508 = _M0L4selfS795->$3;
  _M0L6_2atmpS2507 = _M0L4hashS790 & _M0L14capacity__maskS2508;
  _M0L1iS792 = 0;
  _M0L3idxS793 = _M0L6_2atmpS2507;
  while (1) {
    struct _M0TPB5EntryGicE** _M0L7entriesS2506 = _M0L4selfS795->$0;
    struct _M0TPB5EntryGicE* _M0L7_2abindS794;
    if (
      _M0L3idxS793 < 0
      || _M0L3idxS793 >= Moonbit_array_length(_M0L7entriesS2506)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS794
    = (struct _M0TPB5EntryGicE*)_M0L7entriesS2506[_M0L3idxS793];
    if (_M0L7_2abindS794 == 0) {
      return -1;
    } else {
      struct _M0TPB5EntryGicE* _M0L7_2aSomeS796 = _M0L7_2abindS794;
      struct _M0TPB5EntryGicE* _M0L8_2aentryS797 = _M0L7_2aSomeS796;
      int32_t _M0L4hashS2498 = _M0L8_2aentryS797->$3;
      int32_t _if__result_3340;
      int32_t _M0L3pslS2501;
      int32_t _M0L6_2atmpS2502;
      int32_t _M0L6_2atmpS2504;
      int32_t _M0L14capacity__maskS2505;
      int32_t _M0L6_2atmpS2503;
      if (_M0L4hashS2498 == _M0L4hashS790) {
        int32_t _M0L3keyS2497 = _M0L8_2aentryS797->$4;
        _if__result_3340 = _M0L3keyS2497 == _M0L3keyS791;
      } else {
        _if__result_3340 = 0;
      }
      if (_if__result_3340) {
        int32_t _M0L5valueS2500 = _M0L8_2aentryS797->$5;
        int32_t _M0L6_2atmpS2499 = _M0L5valueS2500;
        return _M0L6_2atmpS2499;
      } else {
        moonbit_incref(_M0L8_2aentryS797);
      }
      _M0L3pslS2501 = _M0L8_2aentryS797->$2;
      moonbit_decref(_M0L8_2aentryS797);
      if (_M0L1iS792 > _M0L3pslS2501) {
        return -1;
      }
      _M0L6_2atmpS2502 = _M0L1iS792 + 1;
      _M0L6_2atmpS2504 = _M0L3idxS793 + 1;
      _M0L14capacity__maskS2505 = _M0L4selfS795->$3;
      _M0L6_2atmpS2503 = _M0L6_2atmpS2504 & _M0L14capacity__maskS2505;
      _M0L1iS792 = _M0L6_2atmpS2502;
      _M0L3idxS793 = _M0L6_2atmpS2503;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGicE* _M0MPB3Map3MapGicE(
  struct _M0TPB9ArrayViewGUicEE _M0L3arrS780,
  int64_t _M0L8capacityS782
) {
  int32_t _M0L3endS2495;
  int32_t _M0L5startS2496;
  int32_t _M0L6lengthS779;
  int32_t _M0L8capacityS781;
  struct _M0TPB3MapGicE* _M0L1mS785;
  int32_t _M0L3endS2492;
  int32_t _M0L5startS2493;
  int32_t _M0L7_2abindS786;
  int32_t _M0L2__S787;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2495 = _M0L3arrS780.$2;
  _M0L5startS2496 = _M0L3arrS780.$1;
  _M0L6lengthS779 = _M0L3endS2495 - _M0L5startS2496;
  if (_M0L8capacityS782 == 4294967296ll) {
    if (_M0L6lengthS779 == 0) {
      _M0L8capacityS781 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS781 = _M0FPB21capacity__for__length(_M0L6lengthS779);
    }
  } else {
    int64_t _M0L7_2aSomeS783 = _M0L8capacityS782;
    int32_t _M0L11_2acapacityS784 = (int32_t)_M0L7_2aSomeS783;
    int32_t _M0L6_2atmpS2494;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2494 = _M0FPB21capacity__for__length(_M0L6lengthS779);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS781
    = _M0MPC13int3Int3max(_M0L11_2acapacityS784, _M0L6_2atmpS2494);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS785 = _M0FPB8new__mapGicE(_M0L8capacityS781);
  _M0L3endS2492 = _M0L3arrS780.$2;
  _M0L5startS2493 = _M0L3arrS780.$1;
  _M0L7_2abindS786 = _M0L3endS2492 - _M0L5startS2493;
  _M0L2__S787 = 0;
  while (1) {
    if (_M0L2__S787 < _M0L7_2abindS786) {
      struct _M0TUicE** _M0L3bufS2489 = _M0L3arrS780.$0;
      int32_t _M0L5startS2491 = _M0L3arrS780.$1;
      int32_t _M0L6_2atmpS2490 = _M0L5startS2491 + _M0L2__S787;
      struct _M0TUicE* _M0L1eS788 =
        (struct _M0TUicE*)_M0L3bufS2489[_M0L6_2atmpS2490];
      int32_t _M0L6_2atmpS2486 = _M0L1eS788->$0;
      int32_t _M0L6_2atmpS2487 = _M0L1eS788->$1;
      int32_t _M0L6_2atmpS2488;
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGicE(_M0L1mS785, _M0L6_2atmpS2486, _M0L6_2atmpS2487);
      _M0L6_2atmpS2488 = _M0L2__S787 + 1;
      _M0L2__S787 = _M0L6_2atmpS2488;
      continue;
    }
    break;
  }
  return _M0L1mS785;
}

int32_t _M0MPB3Map3setGicE(
  struct _M0TPB3MapGicE* _M0L4selfS776,
  int32_t _M0L3keyS777,
  int32_t _M0L5valueS778
) {
  int32_t _M0L6_2atmpS2485;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2485 = _M0IPC13int3IntPB4Hash4hash(_M0L3keyS777);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGicE(_M0L4selfS776, _M0L3keyS777, _M0L5valueS778, _M0L6_2atmpS2485);
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGicE(
  struct _M0TPB3MapGicE* _M0L4selfS763,
  int32_t _M0L3keyS769,
  int32_t _M0L5valueS770,
  int32_t _M0L4hashS765
) {
  int32_t _M0L14capacity__maskS2484;
  int32_t _M0L6_2atmpS2483;
  int32_t _M0L3pslS760;
  int32_t _M0L3idxS761;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2484 = _M0L4selfS763->$3;
  _M0L6_2atmpS2483 = _M0L4hashS765 & _M0L14capacity__maskS2484;
  _M0L3pslS760 = 0;
  _M0L3idxS761 = _M0L6_2atmpS2483;
  while (1) {
    struct _M0TPB5EntryGicE** _M0L7entriesS2482 = _M0L4selfS763->$0;
    struct _M0TPB5EntryGicE* _M0L7_2abindS762;
    if (
      _M0L3idxS761 < 0
      || _M0L3idxS761 >= Moonbit_array_length(_M0L7entriesS2482)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS762
    = (struct _M0TPB5EntryGicE*)_M0L7entriesS2482[_M0L3idxS761];
    if (_M0L7_2abindS762 == 0) {
      int32_t _M0L4sizeS2467 = _M0L4selfS763->$1;
      int32_t _M0L8grow__atS2468 = _M0L4selfS763->$4;
      int32_t _M0L7_2abindS766;
      struct _M0TPB5EntryGicE* _M0L7_2abindS767;
      struct _M0TPB5EntryGicE* _M0L5entryS768;
      if (_M0L4sizeS2467 >= _M0L8grow__atS2468) {
        int32_t _M0L14capacity__maskS2470;
        int32_t _M0L6_2atmpS2469;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGicE(_M0L4selfS763);
        _M0L14capacity__maskS2470 = _M0L4selfS763->$3;
        _M0L6_2atmpS2469 = _M0L4hashS765 & _M0L14capacity__maskS2470;
        _M0L3pslS760 = 0;
        _M0L3idxS761 = _M0L6_2atmpS2469;
        continue;
      }
      _M0L7_2abindS766 = _M0L4selfS763->$6;
      _M0L7_2abindS767 = 0;
      _M0L5entryS768
      = (struct _M0TPB5EntryGicE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGicE));
      Moonbit_object_header(_M0L5entryS768)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 0);
      _M0L5entryS768->$0 = _M0L7_2abindS766;
      _M0L5entryS768->$1 = _M0L7_2abindS767;
      _M0L5entryS768->$2 = _M0L3pslS760;
      _M0L5entryS768->$3 = _M0L4hashS765;
      _M0L5entryS768->$4 = _M0L3keyS769;
      _M0L5entryS768->$5 = _M0L5valueS770;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGicE(_M0L4selfS763, _M0L3idxS761, _M0L5entryS768);
      moonbit_decref(_M0L5entryS768);
      return 0;
    } else {
      struct _M0TPB5EntryGicE* _M0L7_2aSomeS771 = _M0L7_2abindS762;
      struct _M0TPB5EntryGicE* _M0L14_2acurr__entryS772 = _M0L7_2aSomeS771;
      int32_t _M0L4hashS2472 = _M0L14_2acurr__entryS772->$3;
      int32_t _if__result_3343;
      int32_t _M0L3pslS2473;
      int32_t _M0L6_2atmpS2478;
      int32_t _M0L6_2atmpS2480;
      int32_t _M0L14capacity__maskS2481;
      int32_t _M0L6_2atmpS2479;
      if (_M0L4hashS2472 == _M0L4hashS765) {
        int32_t _M0L3keyS2471 = _M0L14_2acurr__entryS772->$4;
        _if__result_3343 = _M0L3keyS2471 == _M0L3keyS769;
      } else {
        _if__result_3343 = 0;
      }
      if (_if__result_3343) {
        _M0L14_2acurr__entryS772->$5 = _M0L5valueS770;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS772);
      }
      _M0L3pslS2473 = _M0L14_2acurr__entryS772->$2;
      if (_M0L3pslS760 > _M0L3pslS2473) {
        int32_t _M0L4sizeS2474 = _M0L4selfS763->$1;
        int32_t _M0L8grow__atS2475 = _M0L4selfS763->$4;
        int32_t _M0L7_2abindS773;
        struct _M0TPB5EntryGicE* _M0L7_2abindS774;
        struct _M0TPB5EntryGicE* _M0L5entryS775;
        if (_M0L4sizeS2474 >= _M0L8grow__atS2475) {
          int32_t _M0L14capacity__maskS2477;
          int32_t _M0L6_2atmpS2476;
          moonbit_decref(_M0L14_2acurr__entryS772);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGicE(_M0L4selfS763);
          _M0L14capacity__maskS2477 = _M0L4selfS763->$3;
          _M0L6_2atmpS2476 = _M0L4hashS765 & _M0L14capacity__maskS2477;
          _M0L3pslS760 = 0;
          _M0L3idxS761 = _M0L6_2atmpS2476;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGicE(_M0L4selfS763, _M0L3idxS761, _M0L14_2acurr__entryS772);
        moonbit_decref(_M0L14_2acurr__entryS772);
        _M0L7_2abindS773 = _M0L4selfS763->$6;
        _M0L7_2abindS774 = 0;
        _M0L5entryS775
        = (struct _M0TPB5EntryGicE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGicE));
        Moonbit_object_header(_M0L5entryS775)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 31, 0);
        _M0L5entryS775->$0 = _M0L7_2abindS773;
        _M0L5entryS775->$1 = _M0L7_2abindS774;
        _M0L5entryS775->$2 = _M0L3pslS760;
        _M0L5entryS775->$3 = _M0L4hashS765;
        _M0L5entryS775->$4 = _M0L3keyS769;
        _M0L5entryS775->$5 = _M0L5valueS770;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGicE(_M0L4selfS763, _M0L3idxS761, _M0L5entryS775);
        moonbit_decref(_M0L5entryS775);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS772);
      }
      _M0L6_2atmpS2478 = _M0L3pslS760 + 1;
      _M0L6_2atmpS2480 = _M0L3idxS761 + 1;
      _M0L14capacity__maskS2481 = _M0L4selfS763->$3;
      _M0L6_2atmpS2479 = _M0L6_2atmpS2480 & _M0L14capacity__maskS2481;
      _M0L3pslS760 = _M0L6_2atmpS2478;
      _M0L3idxS761 = _M0L6_2atmpS2479;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGicE(struct _M0TPB3MapGicE* _M0L4selfS753) {
  struct _M0TPB5EntryGicE* _M0L9old__headS752;
  int32_t _M0L8capacityS2466;
  int32_t _M0L13new__capacityS754;
  struct _M0TPB5EntryGicE* _M0L6_2atmpS2460;
  struct _M0TPB5EntryGicE** _M0L6_2atmpS2459;
  struct _M0TPB5EntryGicE** _M0L6_2aoldS3108;
  int32_t _M0L6_2atmpS2461;
  int32_t _M0L8capacityS2463;
  int32_t _M0L6_2atmpS2462;
  struct _M0TPB5EntryGicE* _M0L6_2atmpS2464;
  struct _M0TPB5EntryGicE* _M0L6_2aoldS3107;
  struct _M0TPB5EntryGicE* _M0L1xS755;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS752 = _M0L4selfS753->$5;
  _M0L8capacityS2466 = _M0L4selfS753->$2;
  _M0L13new__capacityS754 = _M0L8capacityS2466 << 1;
  _M0L6_2atmpS2460 = 0;
  _M0L6_2atmpS2459
  = (struct _M0TPB5EntryGicE**)moonbit_make_ref_array(_M0L13new__capacityS754, _M0L6_2atmpS2460);
  _M0L6_2aoldS3108 = _M0L4selfS753->$0;
  if (_M0L9old__headS752) {
    moonbit_incref(_M0L9old__headS752);
  }
  moonbit_decref(_M0L6_2aoldS3108);
  _M0L4selfS753->$0 = _M0L6_2atmpS2459;
  _M0L4selfS753->$2 = _M0L13new__capacityS754;
  _M0L6_2atmpS2461 = _M0L13new__capacityS754 - 1;
  _M0L4selfS753->$3 = _M0L6_2atmpS2461;
  _M0L8capacityS2463 = _M0L4selfS753->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2462 = _M0FPB21calc__grow__threshold(_M0L8capacityS2463);
  _M0L4selfS753->$4 = _M0L6_2atmpS2462;
  _M0L4selfS753->$1 = 0;
  _M0L6_2atmpS2464 = 0;
  _M0L6_2aoldS3107 = _M0L4selfS753->$5;
  if (_M0L6_2aoldS3107) {
    moonbit_decref(_M0L6_2aoldS3107);
  }
  _M0L4selfS753->$5 = _M0L6_2atmpS2464;
  _M0L4selfS753->$6 = -1;
  _M0L1xS755 = _M0L9old__headS752;
  while (1) {
    if (_M0L1xS755 == 0) {
      if (_M0L1xS755) {
        moonbit_decref(_M0L1xS755);
      }
      break;
    } else {
      struct _M0TPB5EntryGicE* _M0L7_2aSomeS757 = _M0L1xS755;
      struct _M0TPB5EntryGicE* _M0L4_2aeS758 = _M0L7_2aSomeS757;
      struct _M0TPB5EntryGicE* _M0L15next__in__chainS759 = _M0L4_2aeS758->$1;
      struct _M0TPB5EntryGicE* _M0L6_2atmpS2465 = 0;
      struct _M0TPB5EntryGicE* _M0L6_2aoldS3105 = _M0L4_2aeS758->$1;
      if (_M0L15next__in__chainS759) {
        moonbit_incref(_M0L15next__in__chainS759);
      }
      if (_M0L6_2aoldS3105) {
        moonbit_decref(_M0L6_2aoldS3105);
      }
      _M0L4_2aeS758->$1 = _M0L6_2atmpS2465;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGicE(_M0L4selfS753, _M0L4_2aeS758);
      moonbit_decref(_M0L4_2aeS758);
      _M0L1xS755 = _M0L15next__in__chainS759;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGicE(
  struct _M0TPB3MapGicE* _M0L4selfS748,
  struct _M0TPB5EntryGicE* _M0L5outerS744
) {
  int32_t _M0L4hashS743;
  int32_t _M0L14capacity__maskS2458;
  int32_t _M0L6_2atmpS2457;
  int32_t _M0L3pslS745;
  int32_t _M0L3idxS746;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS743 = _M0L5outerS744->$3;
  _M0L14capacity__maskS2458 = _M0L4selfS748->$3;
  _M0L6_2atmpS2457 = _M0L4hashS743 & _M0L14capacity__maskS2458;
  _M0L3pslS745 = 0;
  _M0L3idxS746 = _M0L6_2atmpS2457;
  while (1) {
    struct _M0TPB5EntryGicE** _M0L7entriesS2456 = _M0L4selfS748->$0;
    struct _M0TPB5EntryGicE* _M0L7_2abindS747;
    if (
      _M0L3idxS746 < 0
      || _M0L3idxS746 >= Moonbit_array_length(_M0L7entriesS2456)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS747
    = (struct _M0TPB5EntryGicE*)_M0L7entriesS2456[_M0L3idxS746];
    if (_M0L7_2abindS747 == 0) {
      int32_t _M0L4tailS2449;
      _M0L5outerS744->$2 = _M0L3pslS745;
      _M0L4tailS2449 = _M0L4selfS748->$6;
      _M0L5outerS744->$0 = _M0L4tailS2449;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGicE(_M0L4selfS748, _M0L3idxS746, _M0L5outerS744);
      return 0;
    } else {
      struct _M0TPB5EntryGicE* _M0L7_2aSomeS749 = _M0L7_2abindS747;
      struct _M0TPB5EntryGicE* _M0L7_2acurrS750 = _M0L7_2aSomeS749;
      int32_t _M0L3pslS2450 = _M0L7_2acurrS750->$2;
      if (_M0L3pslS745 > _M0L3pslS2450) {
        int32_t _M0L4tailS2451;
        moonbit_incref(_M0L7_2acurrS750);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGicE(_M0L4selfS748, _M0L3idxS746, _M0L7_2acurrS750);
        moonbit_decref(_M0L7_2acurrS750);
        _M0L5outerS744->$2 = _M0L3pslS745;
        _M0L4tailS2451 = _M0L4selfS748->$6;
        _M0L5outerS744->$0 = _M0L4tailS2451;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGicE(_M0L4selfS748, _M0L3idxS746, _M0L5outerS744);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2452 = _M0L3pslS745 + 1;
        int32_t _M0L6_2atmpS2454 = _M0L3idxS746 + 1;
        int32_t _M0L14capacity__maskS2455 = _M0L4selfS748->$3;
        int32_t _M0L6_2atmpS2453 =
          _M0L6_2atmpS2454 & _M0L14capacity__maskS2455;
        _M0L3pslS745 = _M0L6_2atmpS2452;
        _M0L3idxS746 = _M0L6_2atmpS2453;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGicE(
  struct _M0TPB3MapGicE* _M0L4selfS737,
  int32_t _M0L3idxS742,
  struct _M0TPB5EntryGicE* _M0L5entryS741
) {
  int32_t _M0L3pslS2448;
  int32_t _M0L6_2atmpS2444;
  int32_t _M0L6_2atmpS2446;
  int32_t _M0L14capacity__maskS2447;
  int32_t _M0L6_2atmpS2445;
  int32_t _M0L3pslS733;
  int32_t _M0L3idxS734;
  struct _M0TPB5EntryGicE* _M0L5entryS735;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2448 = _M0L5entryS741->$2;
  _M0L6_2atmpS2444 = _M0L3pslS2448 + 1;
  _M0L6_2atmpS2446 = _M0L3idxS742 + 1;
  _M0L14capacity__maskS2447 = _M0L4selfS737->$3;
  _M0L6_2atmpS2445 = _M0L6_2atmpS2446 & _M0L14capacity__maskS2447;
  moonbit_incref(_M0L5entryS741);
  _M0L3pslS733 = _M0L6_2atmpS2444;
  _M0L3idxS734 = _M0L6_2atmpS2445;
  _M0L5entryS735 = _M0L5entryS741;
  while (1) {
    struct _M0TPB5EntryGicE** _M0L7entriesS2443 = _M0L4selfS737->$0;
    struct _M0TPB5EntryGicE* _M0L7_2abindS736;
    if (
      _M0L3idxS734 < 0
      || _M0L3idxS734 >= Moonbit_array_length(_M0L7entriesS2443)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS736
    = (struct _M0TPB5EntryGicE*)_M0L7entriesS2443[_M0L3idxS734];
    if (_M0L7_2abindS736 == 0) {
      _M0L5entryS735->$2 = _M0L3pslS733;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGicE(_M0L4selfS737, _M0L5entryS735, _M0L3idxS734);
      moonbit_decref(_M0L5entryS735);
      break;
    } else {
      struct _M0TPB5EntryGicE* _M0L7_2aSomeS739 = _M0L7_2abindS736;
      struct _M0TPB5EntryGicE* _M0L14_2acurr__entryS740 = _M0L7_2aSomeS739;
      int32_t _M0L3pslS2433 = _M0L14_2acurr__entryS740->$2;
      if (_M0L3pslS733 > _M0L3pslS2433) {
        int32_t _M0L3pslS2438;
        int32_t _M0L6_2atmpS2434;
        int32_t _M0L6_2atmpS2436;
        int32_t _M0L14capacity__maskS2437;
        int32_t _M0L6_2atmpS2435;
        _M0L5entryS735->$2 = _M0L3pslS733;
        moonbit_incref(_M0L14_2acurr__entryS740);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGicE(_M0L4selfS737, _M0L5entryS735, _M0L3idxS734);
        moonbit_decref(_M0L5entryS735);
        _M0L3pslS2438 = _M0L14_2acurr__entryS740->$2;
        _M0L6_2atmpS2434 = _M0L3pslS2438 + 1;
        _M0L6_2atmpS2436 = _M0L3idxS734 + 1;
        _M0L14capacity__maskS2437 = _M0L4selfS737->$3;
        _M0L6_2atmpS2435 = _M0L6_2atmpS2436 & _M0L14capacity__maskS2437;
        _M0L3pslS733 = _M0L6_2atmpS2434;
        _M0L3idxS734 = _M0L6_2atmpS2435;
        _M0L5entryS735 = _M0L14_2acurr__entryS740;
        continue;
      } else {
        int32_t _M0L6_2atmpS2439 = _M0L3pslS733 + 1;
        int32_t _M0L6_2atmpS2441 = _M0L3idxS734 + 1;
        int32_t _M0L14capacity__maskS2442 = _M0L4selfS737->$3;
        int32_t _M0L6_2atmpS2440 =
          _M0L6_2atmpS2441 & _M0L14capacity__maskS2442;
        struct _M0TPB5EntryGicE* _tmp_3347 = _M0L5entryS735;
        _M0L3pslS733 = _M0L6_2atmpS2439;
        _M0L3idxS734 = _M0L6_2atmpS2440;
        _M0L5entryS735 = _tmp_3347;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGicE(
  struct _M0TPB3MapGicE* _M0L4selfS727,
  struct _M0TPB5EntryGicE* _M0L5entryS729,
  int32_t _M0L8new__idxS728
) {
  struct _M0TPB5EntryGicE** _M0L7entriesS2431;
  struct _M0TPB5EntryGicE* _M0L6_2atmpS2432;
  struct _M0TPB5EntryGicE* _M0L6_2aoldS3115;
  struct _M0TPB5EntryGicE* _M0L7_2abindS730;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2431 = _M0L4selfS727->$0;
  _M0L6_2atmpS2432 = _M0L5entryS729;
  if (
    _M0L8new__idxS728 < 0
    || _M0L8new__idxS728 >= Moonbit_array_length(_M0L7entriesS2431)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3115
  = (struct _M0TPB5EntryGicE*)_M0L7entriesS2431[_M0L8new__idxS728];
  if (_M0L6_2atmpS2432) {
    moonbit_incref(_M0L6_2atmpS2432);
  }
  if (_M0L6_2aoldS3115) {
    moonbit_decref(_M0L6_2aoldS3115);
  }
  _M0L7entriesS2431[_M0L8new__idxS728] = _M0L6_2atmpS2432;
  _M0L7_2abindS730 = _M0L5entryS729->$1;
  if (_M0L7_2abindS730 == 0) {
    _M0L4selfS727->$6 = _M0L8new__idxS728;
  } else {
    struct _M0TPB5EntryGicE* _M0L7_2aSomeS731 = _M0L7_2abindS730;
    struct _M0TPB5EntryGicE* _M0L7_2anextS732 = _M0L7_2aSomeS731;
    _M0L7_2anextS732->$0 = _M0L8new__idxS728;
  }
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGicE(
  struct _M0TPB3MapGicE* _M0L4selfS724,
  int32_t _M0L3idxS726,
  struct _M0TPB5EntryGicE* _M0L5entryS725
) {
  int32_t _M0L7_2abindS723;
  struct _M0TPB5EntryGicE** _M0L7entriesS2427;
  struct _M0TPB5EntryGicE* _M0L6_2atmpS2428;
  struct _M0TPB5EntryGicE* _M0L6_2aoldS3117;
  int32_t _M0L4sizeS2430;
  int32_t _M0L6_2atmpS2429;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS723 = _M0L4selfS724->$6;
  switch (_M0L7_2abindS723) {
    case -1: {
      struct _M0TPB5EntryGicE* _M0L6_2atmpS2422 = _M0L5entryS725;
      struct _M0TPB5EntryGicE* _M0L6_2aoldS3119 = _M0L4selfS724->$5;
      if (_M0L6_2atmpS2422) {
        moonbit_incref(_M0L6_2atmpS2422);
      }
      if (_M0L6_2aoldS3119) {
        moonbit_decref(_M0L6_2aoldS3119);
      }
      _M0L4selfS724->$5 = _M0L6_2atmpS2422;
      break;
    }
    default: {
      struct _M0TPB5EntryGicE** _M0L7entriesS2426 = _M0L4selfS724->$0;
      struct _M0TPB5EntryGicE* _M0L6_2atmpS2425;
      struct _M0TPB5EntryGicE* _M0L6_2atmpS2423;
      struct _M0TPB5EntryGicE* _M0L6_2atmpS2424;
      struct _M0TPB5EntryGicE* _M0L6_2aoldS3120;
      if (
        _M0L7_2abindS723 < 0
        || _M0L7_2abindS723 >= Moonbit_array_length(_M0L7entriesS2426)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2425
      = (struct _M0TPB5EntryGicE*)_M0L7entriesS2426[_M0L7_2abindS723];
      if (_M0L6_2atmpS2425) {
        moonbit_incref(_M0L6_2atmpS2425);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2423
      = _M0MPC16option6Option6unwrapGRPB5EntryGicEE(_M0L6_2atmpS2425);
      if (_M0L6_2atmpS2425) {
        moonbit_decref(_M0L6_2atmpS2425);
      }
      _M0L6_2atmpS2424 = _M0L5entryS725;
      _M0L6_2aoldS3120 = _M0L6_2atmpS2423->$1;
      if (_M0L6_2atmpS2424) {
        moonbit_incref(_M0L6_2atmpS2424);
      }
      if (_M0L6_2aoldS3120) {
        moonbit_decref(_M0L6_2aoldS3120);
      }
      _M0L6_2atmpS2423->$1 = _M0L6_2atmpS2424;
      moonbit_decref(_M0L6_2atmpS2423);
      break;
    }
  }
  _M0L4selfS724->$6 = _M0L3idxS726;
  _M0L7entriesS2427 = _M0L4selfS724->$0;
  _M0L6_2atmpS2428 = _M0L5entryS725;
  if (
    _M0L3idxS726 < 0
    || _M0L3idxS726 >= Moonbit_array_length(_M0L7entriesS2427)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3117
  = (struct _M0TPB5EntryGicE*)_M0L7entriesS2427[_M0L3idxS726];
  if (_M0L6_2atmpS2428) {
    moonbit_incref(_M0L6_2atmpS2428);
  }
  if (_M0L6_2aoldS3117) {
    moonbit_decref(_M0L6_2aoldS3117);
  }
  _M0L7entriesS2427[_M0L3idxS726] = _M0L6_2atmpS2428;
  _M0L4sizeS2430 = _M0L4selfS724->$1;
  _M0L6_2atmpS2429 = _M0L4sizeS2430 + 1;
  _M0L4selfS724->$1 = _M0L6_2atmpS2429;
  return 0;
}

int32_t _M0MPC13int3Int3max(int32_t _M0L4selfS721, int32_t _M0L5otherS722) {
  #line 75 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS721 > _M0L5otherS722) {
    return _M0L4selfS721;
  } else {
    return _M0L5otherS722;
  }
}

int32_t _M0FPB21capacity__for__length(int32_t _M0L6lengthS720) {
  int32_t _M0Lm8capacityS719;
  int32_t _M0L6_2atmpS2420;
  int32_t _M0L6_2atmpS2419;
  #line 71 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 72 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0Lm8capacityS719 = _M0MPC13int3Int20next__power__of__two(_M0L6lengthS720);
  _M0L6_2atmpS2420 = _M0Lm8capacityS719;
  #line 73 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2419 = _M0FPB21calc__grow__threshold(_M0L6_2atmpS2420);
  if (_M0L6lengthS720 > _M0L6_2atmpS2419) {
    int32_t _M0L6_2atmpS2421 = _M0Lm8capacityS719;
    _M0Lm8capacityS719 = _M0L6_2atmpS2421 * 2;
  }
  return _M0Lm8capacityS719;
}

struct _M0TPB3MapGicE* _M0FPB8new__mapGicE(int32_t _M0L8capacityS714) {
  int32_t _M0L8capacityS713;
  int32_t _M0L7_2abindS715;
  int32_t _M0L7_2abindS716;
  struct _M0TPB5EntryGicE* _M0L6_2atmpS2418;
  struct _M0TPB5EntryGicE** _M0L7_2abindS717;
  struct _M0TPB5EntryGicE* _M0L7_2abindS718;
  struct _M0TPB3MapGicE* _block_3348;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS713
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS714);
  _M0L7_2abindS715 = _M0L8capacityS713 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS716 = _M0FPB21calc__grow__threshold(_M0L8capacityS713);
  _M0L6_2atmpS2418 = 0;
  _M0L7_2abindS717
  = (struct _M0TPB5EntryGicE**)moonbit_make_ref_array(_M0L8capacityS713, _M0L6_2atmpS2418);
  _M0L7_2abindS718 = 0;
  _block_3348
  = (struct _M0TPB3MapGicE*)moonbit_malloc(sizeof(struct _M0TPB3MapGicE));
  Moonbit_object_header(_block_3348)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 34, 0);
  _block_3348->$0 = _M0L7_2abindS717;
  _block_3348->$1 = 0;
  _block_3348->$2 = _M0L8capacityS713;
  _block_3348->$3 = _M0L7_2abindS715;
  _block_3348->$4 = _M0L7_2abindS716;
  _block_3348->$5 = _M0L7_2abindS718;
  _block_3348->$6 = -1;
  return _block_3348;
}

int32_t _M0MPC13int3Int20next__power__of__two(int32_t _M0L4selfS712) {
  #line 33 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS712 >= 0) {
    int32_t _M0L6_2atmpS2417;
    int32_t _M0L6_2atmpS2416;
    int32_t _M0L6_2atmpS2415;
    int32_t _M0L6_2atmpS2414;
    if (_M0L4selfS712 <= 1) {
      return 1;
    }
    if (_M0L4selfS712 > 1073741824) {
      return 1073741824;
    }
    _M0L6_2atmpS2417 = _M0L4selfS712 - 1;
    #line 44 "/home/developer/.moon/lib/core/builtin/int.mbt"
    _M0L6_2atmpS2416 = moonbit_clz32(_M0L6_2atmpS2417);
    _M0L6_2atmpS2415 = _M0L6_2atmpS2416 - 1;
    _M0L6_2atmpS2414 = 2147483647 >> (_M0L6_2atmpS2415 & 31);
    return _M0L6_2atmpS2414 + 1;
  } else {
    #line 34 "/home/developer/.moon/lib/core/builtin/int.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB21calc__grow__threshold(int32_t _M0L8capacityS711) {
  int32_t _M0L6_2atmpS2413;
  #line 610 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2413 = _M0L8capacityS711 * 13;
  return _M0L6_2atmpS2413 / 16;
}

struct _M0TPB5EntryGicE* _M0MPC16option6Option6unwrapGRPB5EntryGicEE(
  struct _M0TPB5EntryGicE* _M0L4selfS709
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS709 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGicE* _M0L7_2aSomeS710 = _M0L4selfS709;
    if (_M0L7_2aSomeS710) {
      moonbit_incref(_M0L7_2aSomeS710);
    }
    return _M0L7_2aSomeS710;
  }
}

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t _M0L4selfS708) {
  #line 35 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 36 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC13int3Int18to__string_2einner(_M0L4selfS708, 10);
}

moonbit_string_t _M0MPC16string6String9to__upper(
  moonbit_string_t _M0L4selfS691
) {
  struct _M0TWuEu* _M0L6_2atmpS2410;
  int64_t _M0L7_2abindS690;
  #line 1789 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2410
  = (struct _M0TWuEu*)&_M0MPC16string6String9to__upperC2411l1791$closure.data;
  #line 1791 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS690
  = _M0MPC16string6String8find__by(_M0L4selfS691, _M0L6_2atmpS2410);
  moonbit_decref(_M0L6_2atmpS2410);
  if (_M0L7_2abindS690 == 4294967296ll) {
    moonbit_incref(_M0L4selfS691);
    return _M0L4selfS691;
  } else {
    int64_t _M0L7_2aSomeS693 = _M0L7_2abindS690;
    int32_t _M0L6_2aidxS694 = (int32_t)_M0L7_2aSomeS693;
    int32_t _M0L6_2atmpS2409 = Moonbit_array_length(_M0L4selfS691);
    struct _M0TPB13StringBuilder* _M0L3bufS695;
    int64_t _M0L6_2atmpS2408;
    struct _M0TPC16string10StringView _M0L4headS696;
    moonbit_string_t _M0L6_2atmpS2379;
    int32_t _M0L6_2atmpS2380;
    int32_t _M0L3endS2382;
    int32_t _M0L5startS2383;
    int32_t _M0L6_2atmpS2381;
    struct _M0TPC16string10StringView _M0L7_2abindS697;
    moonbit_string_t _M0L7_2abindS698;
    int32_t _M0L7_2abindS699;
    int32_t _M0L7_2abindS700;
    int32_t _M0L16_2astring__indexS701;
    moonbit_string_t _result_3354;
    #line 1794 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L3bufS695
    = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L6_2atmpS2409);
    _M0L6_2atmpS2408 = (int64_t)_M0L6_2aidxS694;
    #line 1795 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L4headS696
    = _M0MPC16string6String12view_2einner(_M0L4selfS691, 0, _M0L6_2atmpS2408);
    #line 1796 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS2379 = _M0MPC16string10StringView4data(_M0L4headS696);
    #line 1796 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS2380
    = _M0MPC16string10StringView13start__offset(_M0L4headS696);
    _M0L3endS2382 = _M0L4headS696.$2;
    _M0L5startS2383 = _M0L4headS696.$1;
    moonbit_decref(_M0L4headS696.$0);
    _M0L6_2atmpS2381 = _M0L3endS2382 - _M0L5startS2383;
    #line 1796 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L3bufS695, _M0L6_2atmpS2379, _M0L6_2atmpS2380, _M0L6_2atmpS2381);
    moonbit_decref(_M0L6_2atmpS2379);
    #line 1797 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L7_2abindS697
    = _M0MPC16string6String12view_2einner(_M0L4selfS691, _M0L6_2aidxS694, 4294967296ll);
    _M0L7_2abindS698 = _M0L7_2abindS697.$0;
    _M0L7_2abindS699 = _M0L7_2abindS697.$1;
    _M0L7_2abindS700 = _M0L7_2abindS697.$2;
    _M0L16_2astring__indexS701 = _M0L7_2abindS699;
    while (1) {
      if (_M0L16_2astring__indexS701 < _M0L7_2abindS700) {
        int32_t _M0L31_2adecoded__next__string__indexS703;
        int32_t _M0L16_2adecoded__charS704;
        int32_t _M0L7_2abindS706 =
          _M0L7_2abindS698[_M0L16_2astring__indexS701];
        int32_t _M0L6_2atmpS2389 = (int32_t)_M0L7_2abindS706;
        int32_t _if__result_3351;
        int32_t _if__result_3352;
        if (_M0L6_2atmpS2389 >= 55296) {
          int32_t _M0L6_2atmpS2388 = (int32_t)_M0L7_2abindS706;
          _if__result_3351 = _M0L6_2atmpS2388 <= 56319;
        } else {
          _if__result_3351 = 0;
        }
        if (_if__result_3351) {
          int32_t _M0L6_2atmpS2387 = _M0L16_2astring__indexS701 + 1;
          _if__result_3352 = _M0L6_2atmpS2387 < _M0L7_2abindS700;
        } else {
          _if__result_3352 = 0;
        }
        if (_if__result_3352) {
          int32_t _M0L6_2atmpS2404 = _M0L16_2astring__indexS701 + 1;
          int32_t _M0L7_2abindS707 = _M0L7_2abindS698[_M0L6_2atmpS2404];
          int32_t _M0L6_2atmpS2391 = (int32_t)_M0L7_2abindS707;
          int32_t _if__result_3353;
          if (_M0L6_2atmpS2391 >= 56320) {
            int32_t _M0L6_2atmpS2390 = (int32_t)_M0L7_2abindS707;
            _if__result_3353 = _M0L6_2atmpS2390 <= 57343;
          } else {
            _if__result_3353 = 0;
          }
          if (_if__result_3353) {
            int32_t _M0L6_2atmpS2392 = _M0L16_2astring__indexS701 + 2;
            int32_t _M0L6_2atmpS2400 = (int32_t)_M0L7_2abindS706;
            int32_t _M0L6_2atmpS2399 = _M0L6_2atmpS2400 - 55296;
            int32_t _M0L6_2atmpS2397 = _M0L6_2atmpS2399 * 1024;
            int32_t _M0L6_2atmpS2398 = (int32_t)_M0L7_2abindS707;
            int32_t _M0L6_2atmpS2396 = _M0L6_2atmpS2397 + _M0L6_2atmpS2398;
            int32_t _M0L6_2atmpS2395 = _M0L6_2atmpS2396 - 56320;
            int32_t _M0L6_2atmpS2394 = _M0L6_2atmpS2395 + 65536;
            int32_t _M0L6_2atmpS2393;
            #line 1797 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2393
            = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2394);
            _M0L31_2adecoded__next__string__indexS703 = _M0L6_2atmpS2392;
            _M0L16_2adecoded__charS704 = _M0L6_2atmpS2393;
            goto join_702;
          } else {
            int32_t _M0L6_2atmpS2401 = _M0L16_2astring__indexS701 + 1;
            int32_t _M0L6_2atmpS2403 = (int32_t)_M0L7_2abindS706;
            int32_t _M0L6_2atmpS2402;
            #line 1797 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2402
            = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2403);
            _M0L31_2adecoded__next__string__indexS703 = _M0L6_2atmpS2401;
            _M0L16_2adecoded__charS704 = _M0L6_2atmpS2402;
            goto join_702;
          }
        } else {
          int32_t _M0L6_2atmpS2405 = _M0L16_2astring__indexS701 + 1;
          int32_t _M0L6_2atmpS2407 = (int32_t)_M0L7_2abindS706;
          int32_t _M0L6_2atmpS2406;
          #line 1797 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L6_2atmpS2406
          = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2407);
          _M0L31_2adecoded__next__string__indexS703 = _M0L6_2atmpS2405;
          _M0L16_2adecoded__charS704 = _M0L6_2atmpS2406;
          goto join_702;
        }
        goto joinlet_3350;
        join_702:;
        #line 1798 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        if (
          _M0MPC14char4Char20is__ascii__lowercase(_M0L16_2adecoded__charS704)
        ) {
          int32_t _M0L6_2atmpS2386 = _M0L16_2adecoded__charS704;
          int32_t _M0L6_2atmpS2385 = _M0L6_2atmpS2386 - 32;
          int32_t _M0L6_2atmpS2384 = _M0L6_2atmpS2385;
          #line 1799 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS695, _M0L6_2atmpS2384);
        } else {
          #line 1801 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS695, _M0L16_2adecoded__charS704);
        }
        _M0L16_2astring__indexS701
        = _M0L31_2adecoded__next__string__indexS703;
        continue;
        joinlet_3350:;
      } else {
        moonbit_decref(_M0L7_2abindS698);
      }
      break;
    }
    #line 1804 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_3354 = _M0MPB13StringBuilder10to__string(_M0L3bufS695);
    moonbit_decref(_M0L3bufS695);
    return _result_3354;
  }
}

int32_t _M0MPC16string6String9to__upperC2411l1791(
  struct _M0TWuEu* _M0L6_2aenvS2412,
  int32_t _M0L1cS692
) {
  #line 1791 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1791 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  return _M0MPC14char4Char20is__ascii__lowercase(_M0L1cS692);
}

int32_t _M0MPC14char4Char20is__ascii__lowercase(int32_t _M0L4selfS689) {
  #line 103 "/home/developer/.moon/lib/core/builtin/char.mbt"
  return _M0L4selfS689 >= 97 && _M0L4selfS689 <= 122 || 0;
}

moonbit_string_t _M0MPC16string6String9to__lower(
  moonbit_string_t _M0L4selfS672
) {
  struct _M0TWuEu* _M0L6_2atmpS2376;
  int64_t _M0L7_2abindS671;
  #line 1734 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2376
  = (struct _M0TWuEu*)&_M0MPC16string6String9to__lowerC2377l1736$closure.data;
  #line 1736 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS671
  = _M0MPC16string6String8find__by(_M0L4selfS672, _M0L6_2atmpS2376);
  moonbit_decref(_M0L6_2atmpS2376);
  if (_M0L7_2abindS671 == 4294967296ll) {
    moonbit_incref(_M0L4selfS672);
    return _M0L4selfS672;
  } else {
    int64_t _M0L7_2aSomeS674 = _M0L7_2abindS671;
    int32_t _M0L6_2aidxS675 = (int32_t)_M0L7_2aSomeS674;
    int32_t _M0L6_2atmpS2375 = Moonbit_array_length(_M0L4selfS672);
    struct _M0TPB13StringBuilder* _M0L3bufS676;
    int64_t _M0L6_2atmpS2374;
    struct _M0TPC16string10StringView _M0L4headS677;
    moonbit_string_t _M0L6_2atmpS2345;
    int32_t _M0L6_2atmpS2346;
    int32_t _M0L3endS2348;
    int32_t _M0L5startS2349;
    int32_t _M0L6_2atmpS2347;
    struct _M0TPC16string10StringView _M0L7_2abindS678;
    moonbit_string_t _M0L7_2abindS679;
    int32_t _M0L7_2abindS680;
    int32_t _M0L7_2abindS681;
    int32_t _M0L16_2astring__indexS682;
    moonbit_string_t _result_3360;
    #line 1739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L3bufS676
    = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L6_2atmpS2375);
    _M0L6_2atmpS2374 = (int64_t)_M0L6_2aidxS675;
    #line 1740 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L4headS677
    = _M0MPC16string6String12view_2einner(_M0L4selfS672, 0, _M0L6_2atmpS2374);
    #line 1741 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS2345 = _M0MPC16string10StringView4data(_M0L4headS677);
    #line 1741 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS2346
    = _M0MPC16string10StringView13start__offset(_M0L4headS677);
    _M0L3endS2348 = _M0L4headS677.$2;
    _M0L5startS2349 = _M0L4headS677.$1;
    moonbit_decref(_M0L4headS677.$0);
    _M0L6_2atmpS2347 = _M0L3endS2348 - _M0L5startS2349;
    #line 1741 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L3bufS676, _M0L6_2atmpS2345, _M0L6_2atmpS2346, _M0L6_2atmpS2347);
    moonbit_decref(_M0L6_2atmpS2345);
    #line 1742 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L7_2abindS678
    = _M0MPC16string6String12view_2einner(_M0L4selfS672, _M0L6_2aidxS675, 4294967296ll);
    _M0L7_2abindS679 = _M0L7_2abindS678.$0;
    _M0L7_2abindS680 = _M0L7_2abindS678.$1;
    _M0L7_2abindS681 = _M0L7_2abindS678.$2;
    _M0L16_2astring__indexS682 = _M0L7_2abindS680;
    while (1) {
      if (_M0L16_2astring__indexS682 < _M0L7_2abindS681) {
        int32_t _M0L31_2adecoded__next__string__indexS684;
        int32_t _M0L16_2adecoded__charS685;
        int32_t _M0L7_2abindS687 =
          _M0L7_2abindS679[_M0L16_2astring__indexS682];
        int32_t _M0L6_2atmpS2355 = (int32_t)_M0L7_2abindS687;
        int32_t _if__result_3357;
        int32_t _if__result_3358;
        if (_M0L6_2atmpS2355 >= 55296) {
          int32_t _M0L6_2atmpS2354 = (int32_t)_M0L7_2abindS687;
          _if__result_3357 = _M0L6_2atmpS2354 <= 56319;
        } else {
          _if__result_3357 = 0;
        }
        if (_if__result_3357) {
          int32_t _M0L6_2atmpS2353 = _M0L16_2astring__indexS682 + 1;
          _if__result_3358 = _M0L6_2atmpS2353 < _M0L7_2abindS681;
        } else {
          _if__result_3358 = 0;
        }
        if (_if__result_3358) {
          int32_t _M0L6_2atmpS2370 = _M0L16_2astring__indexS682 + 1;
          int32_t _M0L7_2abindS688 = _M0L7_2abindS679[_M0L6_2atmpS2370];
          int32_t _M0L6_2atmpS2357 = (int32_t)_M0L7_2abindS688;
          int32_t _if__result_3359;
          if (_M0L6_2atmpS2357 >= 56320) {
            int32_t _M0L6_2atmpS2356 = (int32_t)_M0L7_2abindS688;
            _if__result_3359 = _M0L6_2atmpS2356 <= 57343;
          } else {
            _if__result_3359 = 0;
          }
          if (_if__result_3359) {
            int32_t _M0L6_2atmpS2358 = _M0L16_2astring__indexS682 + 2;
            int32_t _M0L6_2atmpS2366 = (int32_t)_M0L7_2abindS687;
            int32_t _M0L6_2atmpS2365 = _M0L6_2atmpS2366 - 55296;
            int32_t _M0L6_2atmpS2363 = _M0L6_2atmpS2365 * 1024;
            int32_t _M0L6_2atmpS2364 = (int32_t)_M0L7_2abindS688;
            int32_t _M0L6_2atmpS2362 = _M0L6_2atmpS2363 + _M0L6_2atmpS2364;
            int32_t _M0L6_2atmpS2361 = _M0L6_2atmpS2362 - 56320;
            int32_t _M0L6_2atmpS2360 = _M0L6_2atmpS2361 + 65536;
            int32_t _M0L6_2atmpS2359;
            #line 1742 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2359
            = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2360);
            _M0L31_2adecoded__next__string__indexS684 = _M0L6_2atmpS2358;
            _M0L16_2adecoded__charS685 = _M0L6_2atmpS2359;
            goto join_683;
          } else {
            int32_t _M0L6_2atmpS2367 = _M0L16_2astring__indexS682 + 1;
            int32_t _M0L6_2atmpS2369 = (int32_t)_M0L7_2abindS687;
            int32_t _M0L6_2atmpS2368;
            #line 1742 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2368
            = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2369);
            _M0L31_2adecoded__next__string__indexS684 = _M0L6_2atmpS2367;
            _M0L16_2adecoded__charS685 = _M0L6_2atmpS2368;
            goto join_683;
          }
        } else {
          int32_t _M0L6_2atmpS2371 = _M0L16_2astring__indexS682 + 1;
          int32_t _M0L6_2atmpS2373 = (int32_t)_M0L7_2abindS687;
          int32_t _M0L6_2atmpS2372;
          #line 1742 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L6_2atmpS2372
          = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2373);
          _M0L31_2adecoded__next__string__indexS684 = _M0L6_2atmpS2371;
          _M0L16_2adecoded__charS685 = _M0L6_2atmpS2372;
          goto join_683;
        }
        goto joinlet_3356;
        join_683:;
        #line 1743 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        if (
          _M0MPC14char4Char20is__ascii__uppercase(_M0L16_2adecoded__charS685)
        ) {
          int32_t _M0L6_2atmpS2352 = _M0L16_2adecoded__charS685;
          int32_t _M0L6_2atmpS2351 = _M0L6_2atmpS2352 + 32;
          int32_t _M0L6_2atmpS2350 = _M0L6_2atmpS2351;
          #line 1745 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS676, _M0L6_2atmpS2350);
        } else {
          #line 1747 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS676, _M0L16_2adecoded__charS685);
        }
        _M0L16_2astring__indexS682
        = _M0L31_2adecoded__next__string__indexS684;
        continue;
        joinlet_3356:;
      } else {
        moonbit_decref(_M0L7_2abindS679);
      }
      break;
    }
    #line 1750 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_3360 = _M0MPB13StringBuilder10to__string(_M0L3bufS676);
    moonbit_decref(_M0L3bufS676);
    return _result_3360;
  }
}

int32_t _M0MPC16string6String9to__lowerC2377l1736(
  struct _M0TWuEu* _M0L6_2aenvS2378,
  int32_t _M0L1xS673
) {
  #line 1736 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1736 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  return _M0MPC14char4Char20is__ascii__uppercase(_M0L1xS673);
}

int32_t _M0MPC14char4Char20is__ascii__uppercase(int32_t _M0L4selfS670) {
  #line 131 "/home/developer/.moon/lib/core/builtin/char.mbt"
  return _M0L4selfS670 >= 65 && _M0L4selfS670 <= 90 || 0;
}

moonbit_string_t _M0MPC16string6String12replace__all(
  moonbit_string_t _M0L4selfS645,
  struct _M0TPC16string10StringView _M0L3oldS648,
  struct _M0TPC16string10StringView _M0L3newS650
) {
  int32_t _M0L3lenS644;
  struct _M0TPB13StringBuilder* _M0L3bufS646;
  int32_t _M0L3endS2343;
  int32_t _M0L5startS2344;
  int32_t _M0L8old__lenS647;
  moonbit_string_t _M0L3newS649;
  #line 1563 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3lenS644 = Moonbit_array_length(_M0L4selfS645);
  #line 1569 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3bufS646 = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L3lenS644);
  _M0L3endS2343 = _M0L3oldS648.$2;
  _M0L5startS2344 = _M0L3oldS648.$1;
  _M0L8old__lenS647 = _M0L3endS2343 - _M0L5startS2344;
  #line 1571 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3newS649 = _M0MPC16string10StringView9to__owned(_M0L3newS650);
  if (_M0L8old__lenS647 == 0) {
    int32_t _M0L7_2abindS651;
    int32_t _M0L16_2astring__indexS652;
    moonbit_string_t _result_3366;
    #line 1574 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS646, _M0L3newS649);
    _M0L7_2abindS651 = Moonbit_array_length(_M0L4selfS645);
    _M0L16_2astring__indexS652 = 0;
    while (1) {
      if (_M0L16_2astring__indexS652 < _M0L7_2abindS651) {
        int32_t _M0L31_2adecoded__next__string__indexS654;
        int32_t _M0L16_2adecoded__charS655;
        int32_t _M0L7_2abindS657 = _M0L4selfS645[_M0L16_2astring__indexS652];
        int32_t _M0L6_2atmpS2309 = (int32_t)_M0L7_2abindS657;
        int32_t _if__result_3363;
        int32_t _if__result_3364;
        if (_M0L6_2atmpS2309 >= 55296) {
          int32_t _M0L6_2atmpS2308 = (int32_t)_M0L7_2abindS657;
          _if__result_3363 = _M0L6_2atmpS2308 <= 56319;
        } else {
          _if__result_3363 = 0;
        }
        if (_if__result_3363) {
          int32_t _M0L6_2atmpS2307 = _M0L16_2astring__indexS652 + 1;
          _if__result_3364 = _M0L6_2atmpS2307 < _M0L7_2abindS651;
        } else {
          _if__result_3364 = 0;
        }
        if (_if__result_3364) {
          int32_t _M0L6_2atmpS2324 = _M0L16_2astring__indexS652 + 1;
          int32_t _M0L7_2abindS658 = _M0L4selfS645[_M0L6_2atmpS2324];
          int32_t _M0L6_2atmpS2311 = (int32_t)_M0L7_2abindS658;
          int32_t _if__result_3365;
          if (_M0L6_2atmpS2311 >= 56320) {
            int32_t _M0L6_2atmpS2310 = (int32_t)_M0L7_2abindS658;
            _if__result_3365 = _M0L6_2atmpS2310 <= 57343;
          } else {
            _if__result_3365 = 0;
          }
          if (_if__result_3365) {
            int32_t _M0L6_2atmpS2312 = _M0L16_2astring__indexS652 + 2;
            int32_t _M0L6_2atmpS2320 = (int32_t)_M0L7_2abindS657;
            int32_t _M0L6_2atmpS2319 = _M0L6_2atmpS2320 - 55296;
            int32_t _M0L6_2atmpS2317 = _M0L6_2atmpS2319 * 1024;
            int32_t _M0L6_2atmpS2318 = (int32_t)_M0L7_2abindS658;
            int32_t _M0L6_2atmpS2316 = _M0L6_2atmpS2317 + _M0L6_2atmpS2318;
            int32_t _M0L6_2atmpS2315 = _M0L6_2atmpS2316 - 56320;
            int32_t _M0L6_2atmpS2314 = _M0L6_2atmpS2315 + 65536;
            int32_t _M0L6_2atmpS2313;
            #line 1575 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2313
            = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2314);
            _M0L31_2adecoded__next__string__indexS654 = _M0L6_2atmpS2312;
            _M0L16_2adecoded__charS655 = _M0L6_2atmpS2313;
            goto join_653;
          } else {
            int32_t _M0L6_2atmpS2321 = _M0L16_2astring__indexS652 + 1;
            int32_t _M0L6_2atmpS2323 = (int32_t)_M0L7_2abindS657;
            int32_t _M0L6_2atmpS2322;
            #line 1575 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2322
            = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2323);
            _M0L31_2adecoded__next__string__indexS654 = _M0L6_2atmpS2321;
            _M0L16_2adecoded__charS655 = _M0L6_2atmpS2322;
            goto join_653;
          }
        } else {
          int32_t _M0L6_2atmpS2325 = _M0L16_2astring__indexS652 + 1;
          int32_t _M0L6_2atmpS2327 = (int32_t)_M0L7_2abindS657;
          int32_t _M0L6_2atmpS2326;
          #line 1575 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L6_2atmpS2326
          = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2327);
          _M0L31_2adecoded__next__string__indexS654 = _M0L6_2atmpS2325;
          _M0L16_2adecoded__charS655 = _M0L6_2atmpS2326;
          goto join_653;
        }
        goto joinlet_3362;
        join_653:;
        #line 1576 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS646, _M0L16_2adecoded__charS655);
        #line 1577 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS646, _M0L3newS649);
        _M0L16_2astring__indexS652
        = _M0L31_2adecoded__next__string__indexS654;
        continue;
        joinlet_3362:;
      } else {
        moonbit_decref(_M0L3newS649);
      }
      break;
    }
    #line 1579 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_3366 = _M0MPB13StringBuilder10to__string(_M0L3bufS646);
    moonbit_decref(_M0L3bufS646);
    return _result_3366;
  } else {
    int64_t _M0L10first__endS659;
    #line 1581 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L10first__endS659
    = _M0MPC16string6String4find(_M0L4selfS645, _M0L3oldS648);
    if (_M0L10first__endS659 == 4294967296ll) {
      moonbit_decref(_M0L3newS649);
      moonbit_decref(_M0L3bufS646);
      moonbit_incref(_M0L4selfS645);
      return _M0L4selfS645;
    } else {
      int64_t _M0L7_2aSomeS660 = _M0L10first__endS659;
      int32_t _M0L6_2aendS661 = (int32_t)_M0L7_2aSomeS660;
      int32_t _M0L6_2atmpS2342 = Moonbit_array_length(_M0L4selfS645);
      struct _M0TPC16string10StringView _M0L6_2atmpS2341;
      struct _M0TPC16string10StringView _M0L4viewS662;
      int32_t _M0L3endS663;
      moonbit_string_t _result_3368;
      moonbit_incref(_M0L4selfS645);
      _M0L6_2atmpS2341
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L4selfS645, .$1 = 0, .$2 = _M0L6_2atmpS2342
      };
      _M0L4viewS662 = _M0L6_2atmpS2341;
      _M0L3endS663 = _M0L6_2aendS661;
      while (1) {
        int64_t _M0L6_2atmpS2340 = (int64_t)_M0L3endS663;
        struct _M0TPC16string10StringView _M0L3segS664;
        moonbit_string_t _M0L6_2atmpS2328;
        int32_t _M0L6_2atmpS2329;
        int32_t _M0L3endS2331;
        int32_t _M0L5startS2332;
        int32_t _M0L6_2atmpS2330;
        int32_t _M0L6_2atmpS2333;
        #line 1584 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L3segS664
        = _M0MPC16string10StringView12view_2einner(_M0L4viewS662, 0, _M0L6_2atmpS2340);
        #line 1585 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2328 = _M0MPC16string10StringView4data(_M0L3segS664);
        #line 1585 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2329
        = _M0MPC16string10StringView13start__offset(_M0L3segS664);
        _M0L3endS2331 = _M0L3segS664.$2;
        _M0L5startS2332 = _M0L3segS664.$1;
        moonbit_decref(_M0L3segS664.$0);
        _M0L6_2atmpS2330 = _M0L3endS2331 - _M0L5startS2332;
        #line 1585 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L3bufS646, _M0L6_2atmpS2328, _M0L6_2atmpS2329, _M0L6_2atmpS2330);
        moonbit_decref(_M0L6_2atmpS2328);
        #line 1586 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS646, _M0L3newS649);
        _M0L6_2atmpS2333 = _M0L3endS663 + _M0L8old__lenS647;
        if (_M0L6_2atmpS2333 <= _M0L3lenS644) {
          int32_t _M0L6_2atmpS2339 = _M0L3endS663 + _M0L8old__lenS647;
          struct _M0TPC16string10StringView _M0L10next__viewS665;
          int64_t _M0L7_2abindS666;
          #line 1589 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L10next__viewS665
          = _M0MPC16string10StringView12view_2einner(_M0L4viewS662, _M0L6_2atmpS2339, 4294967296ll);
          moonbit_decref(_M0L4viewS662.$0);
          #line 1590 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L7_2abindS666
          = _M0MPC16string10StringView4find(_M0L10next__viewS665, _M0L3oldS648);
          if (_M0L7_2abindS666 == 4294967296ll) {
            moonbit_string_t _M0L6_2atmpS2334;
            int32_t _M0L6_2atmpS2335;
            int32_t _M0L3endS2337;
            int32_t _M0L5startS2338;
            int32_t _M0L6_2atmpS2336;
            moonbit_decref(_M0L3newS649);
            #line 1592 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2334
            = _M0MPC16string10StringView4data(_M0L10next__viewS665);
            #line 1593 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0L6_2atmpS2335
            = _M0MPC16string10StringView13start__offset(_M0L10next__viewS665);
            _M0L3endS2337 = _M0L10next__viewS665.$2;
            _M0L5startS2338 = _M0L10next__viewS665.$1;
            moonbit_decref(_M0L10next__viewS665.$0);
            _M0L6_2atmpS2336 = _M0L3endS2337 - _M0L5startS2338;
            #line 1591 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L3bufS646, _M0L6_2atmpS2334, _M0L6_2atmpS2335, _M0L6_2atmpS2336);
            moonbit_decref(_M0L6_2atmpS2334);
            break;
          } else {
            int64_t _M0L7_2aSomeS667 = _M0L7_2abindS666;
            int32_t _M0L12_2anext__endS668 = (int32_t)_M0L7_2aSomeS667;
            _M0L4viewS662 = _M0L10next__viewS665;
            _M0L3endS663 = _M0L12_2anext__endS668;
            continue;
          }
        } else {
          moonbit_decref(_M0L4viewS662.$0);
          moonbit_decref(_M0L3newS649);
          break;
        }
        break;
      }
      #line 1600 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _result_3368 = _M0MPB13StringBuilder10to__string(_M0L3bufS646);
      moonbit_decref(_M0L3bufS646);
      return _result_3368;
    }
  }
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string6String5split(
  moonbit_string_t _M0L4selfS642,
  struct _M0TPC16string10StringView _M0L3sepS643
) {
  int32_t _M0L6_2atmpS2306;
  struct _M0TPC16string10StringView _M0L6_2atmpS2305;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_3369;
  #line 1168 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2306 = Moonbit_array_length(_M0L4selfS642);
  moonbit_incref(_M0L4selfS642);
  _M0L6_2atmpS2305
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS642, .$1 = 0, .$2 = _M0L6_2atmpS2306
  };
  #line 1169 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3369
  = _M0MPC16string10StringView5split(_M0L6_2atmpS2305, _M0L3sepS643);
  moonbit_decref(_M0L6_2atmpS2305.$0);
  return _result_3369;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string10StringView5split(
  struct _M0TPC16string10StringView _M0L4selfS633,
  struct _M0TPC16string10StringView _M0L3sepS632
) {
  int32_t _M0L3endS2303;
  int32_t _M0L5startS2304;
  int32_t _M0L8sep__lenS631;
  void* _M0L4SomeS2302;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS635;
  struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__* _closure_3371;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS2292;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_3372;
  #line 1139 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2303 = _M0L3sepS632.$2;
  _M0L5startS2304 = _M0L3sepS632.$1;
  _M0L8sep__lenS631 = _M0L3endS2303 - _M0L5startS2304;
  if (_M0L8sep__lenS631 == 0) {
    struct _M0TPB4IterGcE* _M0L6_2atmpS2287;
    struct _M0TWcERPC16string10StringView* _M0L6_2atmpS2288;
    struct _M0TPB4IterGRPC16string10StringViewE* _result_3370;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS2287 = _M0MPC16string10StringView4iter(_M0L4selfS633);
    _M0L6_2atmpS2288
    = (struct _M0TWcERPC16string10StringView*)&_M0MPC16string10StringView5splitC2289l1145$closure.data;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_3370
    = _M0MPB4Iter3mapGcRPC16string10StringViewE(_M0L6_2atmpS2287, _M0L6_2atmpS2288);
    moonbit_decref(_M0L6_2atmpS2287);
    moonbit_decref(_M0L6_2atmpS2288);
    return _result_3370;
  }
  moonbit_incref(_M0L4selfS633.$0);
  _M0L4SomeS2302
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
  Moonbit_object_header(_M0L4SomeS2302)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 41, 1);
  ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS2302)->$0
  = _M0L4selfS633;
  _M0L9remainingS635
  = (struct _M0TPB8MutLocalGORPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE));
  Moonbit_object_header(_M0L9remainingS635)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 44, 0);
  _M0L9remainingS635->$0 = _M0L4SomeS2302;
  moonbit_incref(_M0L3sepS632.$0);
  _closure_3371
  = (struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__*)moonbit_malloc(sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__));
  Moonbit_object_header(_closure_3371)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 47, 0);
  _closure_3371->code = &_M0MPC16string10StringView5splitC2293l1148;
  _closure_3371->$0 = _M0L9remainingS635;
  _closure_3371->$1 = _M0L3sepS632;
  _closure_3371->$2 = _M0L8sep__lenS631;
  _M0L6_2atmpS2292
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_3371;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3372
  = _M0MPB4Iter3newGRPC16string10StringViewE(_M0L6_2atmpS2292, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS2292);
  return _result_3372;
}

void* _M0MPC16string10StringView5splitC2293l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS2294
) {
  struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__* _M0L14_2acasted__envS2295;
  int32_t _M0L8sep__lenS631;
  struct _M0TPC16string10StringView _M0L3sepS632;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS635;
  void* _M0L7_2abindS636;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L14_2acasted__envS2295
  = (struct _M0R44StringView_3a_3asplit_2eanon__u2293__l1148__*)_M0L6_2aenvS2294;
  _M0L8sep__lenS631 = _M0L14_2acasted__envS2295->$2;
  _M0L3sepS632 = _M0L14_2acasted__envS2295->$1;
  _M0L9remainingS635 = _M0L14_2acasted__envS2295->$0;
  _M0L7_2abindS636 = _M0L9remainingS635->$0;
  switch (Moonbit_object_tag(_M0L7_2abindS636)) {
    case 1: {
      struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS637 =
        (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS636;
      struct _M0TPC16string10StringView _M0L7_2aviewS638 =
        _M0L7_2aSomeS637->$0;
      int64_t _M0L7_2abindS639;
      moonbit_incref(_M0L7_2aviewS638.$0);
      #line 1150 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS639
      = _M0MPC16string10StringView4find(_M0L7_2aviewS638, _M0L3sepS632);
      if (_M0L7_2abindS639 == 4294967296ll) {
        void* _M0L4NoneS2296 =
          (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
        void* _M0L6_2aoldS3125 = _M0L9remainingS635->$0;
        void* _block_3373;
        moonbit_decref(_M0L6_2aoldS3125);
        _M0L9remainingS635->$0 = _M0L4NoneS2296;
        _block_3373
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_3373)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 41, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3373)->$0
        = _M0L7_2aviewS638;
        return _block_3373;
      } else {
        int64_t _M0L7_2aSomeS640 = _M0L7_2abindS639;
        int32_t _M0L6_2aendS641 = (int32_t)_M0L7_2aSomeS640;
        int32_t _M0L6_2atmpS2299 = _M0L6_2aendS641 + _M0L8sep__lenS631;
        struct _M0TPC16string10StringView _M0L6_2atmpS2298;
        void* _M0L4SomeS2297;
        void* _M0L6_2aoldS3126;
        int64_t _M0L6_2atmpS2301;
        struct _M0TPC16string10StringView _M0L6_2atmpS2300;
        void* _block_3374;
        #line 1154 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2298
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS638, _M0L6_2atmpS2299, 4294967296ll);
        _M0L4SomeS2297
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_M0L4SomeS2297)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 41, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS2297)->$0
        = _M0L6_2atmpS2298;
        _M0L6_2aoldS3126 = _M0L9remainingS635->$0;
        moonbit_decref(_M0L6_2aoldS3126);
        _M0L9remainingS635->$0 = _M0L4SomeS2297;
        _M0L6_2atmpS2301 = (int64_t)_M0L6_2aendS641;
        #line 1155 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2300
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS638, 0, _M0L6_2atmpS2301);
        moonbit_decref(_M0L7_2aviewS638.$0);
        _block_3374
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_3374)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 41, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3374)->$0
        = _M0L6_2atmpS2300;
        return _block_3374;
      }
      break;
    }
    default: {
      return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
      break;
    }
  }
}

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC2289l1145(
  struct _M0TWcERPC16string10StringView* _M0L6_2aenvS2290,
  int32_t _M0L1cS634
) {
  moonbit_string_t _M0L6_2atmpS2291;
  struct _M0TPC16string10StringView _result_3375;
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2291 = _M0IPC14char4CharPB4Show10to__string(_M0L1cS634);
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3375
  = _M0MPC16string6String12view_2einner(_M0L6_2atmpS2291, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS2291);
  return _result_3375;
}

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t _M0L4selfS630) {
  #line 446 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 447 "/home/developer/.moon/lib/core/builtin/char.mbt"
  return _M0FPB16char__to__string(_M0L4selfS630);
}

moonbit_string_t _M0FPB16char__to__string(int32_t _M0L4charS629) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS628;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS2286;
  moonbit_string_t _result_3376;
  #line 452 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0L7_2aselfS628 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS628, _M0L4charS629);
  _M0L6_2atmpS2286 = _M0L7_2aselfS628;
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _result_3376 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS2286);
  moonbit_decref(_M0L6_2atmpS2286);
  return _result_3376;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE* _M0L4selfS624,
  struct _M0TWcERPC16string10StringView* _M0L1fS627
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__* _closure_3377;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS2280;
  int64_t _M0L10size__hintS2281;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_3378;
  #line 389 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  moonbit_incref(_M0L1fS627);
  moonbit_incref(_M0L4selfS624);
  _closure_3377
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__*)moonbit_malloc(sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__));
  Moonbit_object_header(_closure_3377)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 51, 0);
  _closure_3377->code = &_M0MPB4Iter3mapGcRPC16string10StringViewEC2282l391;
  _closure_3377->$0 = _M0L1fS627;
  _closure_3377->$1 = _M0L4selfS624;
  _M0L6_2atmpS2280
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_3377;
  _M0L10size__hintS2281 = _M0L4selfS624->$1;
  _block_3378
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_3378)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 55, 0);
  _block_3378->$0 = _M0L6_2atmpS2280;
  _block_3378->$1 = _M0L10size__hintS2281;
  return _block_3378;
}

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC2282l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS2283
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__* _M0L14_2acasted__envS2284;
  struct _M0TPB4IterGcE* _M0L4selfS624;
  struct _M0TWcERPC16string10StringView* _M0L1fS627;
  int32_t _M0L7_2abindS623;
  #line 391 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L14_2acasted__envS2284
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u2282__l391__*)_M0L6_2aenvS2283;
  _M0L4selfS624 = _M0L14_2acasted__envS2284->$1;
  _M0L1fS627 = _M0L14_2acasted__envS2284->$0;
  #line 392 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS623 = _M0MPB4Iter4nextGcE(_M0L4selfS624);
  if (_M0L7_2abindS623 == -1) {
    return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  } else {
    int32_t _M0L7_2aSomeS625 = _M0L7_2abindS623;
    int32_t _M0L4_2axS626 = _M0L7_2aSomeS625;
    struct _M0TPC16string10StringView _M0L6_2atmpS2285;
    void* _block_3379;
    #line 393 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6_2atmpS2285 = _M0L1fS627->code(_M0L1fS627, _M0L4_2axS626);
    _block_3379
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
    Moonbit_object_header(_block_3379)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 41, 1);
    ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3379)->$0
    = _M0L6_2atmpS2285;
    return _block_3379;
  }
}

struct _M0TPB4IterGcE* _M0MPC16string6String4iter(
  moonbit_string_t _M0L4selfS618
) {
  int32_t _M0L3lenS617;
  struct _M0TPB8MutLocalGiE* _M0L5indexS619;
  struct _M0R38String_3a_3aiter_2eanon__u2264__l256__* _closure_3380;
  struct _M0TWEOc* _M0L6_2atmpS2263;
  struct _M0TPB4IterGcE* _result_3381;
  #line 252 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L3lenS617 = Moonbit_array_length(_M0L4selfS618);
  _M0L5indexS619
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS619)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS619->$0 = 0;
  moonbit_incref(_M0L4selfS618);
  _closure_3380
  = (struct _M0R38String_3a_3aiter_2eanon__u2264__l256__*)moonbit_malloc(sizeof(struct _M0R38String_3a_3aiter_2eanon__u2264__l256__));
  Moonbit_object_header(_closure_3380)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 58, 0);
  _closure_3380->code = &_M0MPC16string6String4iterC2264l256;
  _closure_3380->$0 = _M0L5indexS619;
  _closure_3380->$1 = _M0L4selfS618;
  _closure_3380->$2 = _M0L3lenS617;
  _M0L6_2atmpS2263 = (struct _M0TWEOc*)_closure_3380;
  #line 256 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_3381 = _M0MPB4Iter3newGcE(_M0L6_2atmpS2263, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS2263);
  return _result_3381;
}

int32_t _M0MPC16string6String4iterC2264l256(
  struct _M0TWEOc* _M0L6_2aenvS2265
) {
  struct _M0R38String_3a_3aiter_2eanon__u2264__l256__* _M0L14_2acasted__envS2266;
  int32_t _M0L3lenS617;
  moonbit_string_t _M0L4selfS618;
  struct _M0TPB8MutLocalGiE* _M0L5indexS619;
  int32_t _M0L3valS2267;
  #line 256 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L14_2acasted__envS2266
  = (struct _M0R38String_3a_3aiter_2eanon__u2264__l256__*)_M0L6_2aenvS2265;
  _M0L3lenS617 = _M0L14_2acasted__envS2266->$2;
  _M0L4selfS618 = _M0L14_2acasted__envS2266->$1;
  _M0L5indexS619 = _M0L14_2acasted__envS2266->$0;
  _M0L3valS2267 = _M0L5indexS619->$0;
  if (_M0L3valS2267 < _M0L3lenS617) {
    int32_t _M0L3valS2279 = _M0L5indexS619->$0;
    int32_t _M0L2c1S620 = _M0L4selfS618[_M0L3valS2279];
    int32_t _if__result_3382;
    int32_t _M0L3valS2277;
    int32_t _M0L6_2atmpS2276;
    int32_t _M0L6_2atmpS2278;
    #line 259 "/home/developer/.moon/lib/core/builtin/string.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S620)) {
      int32_t _M0L3valS2269 = _M0L5indexS619->$0;
      int32_t _M0L6_2atmpS2268 = _M0L3valS2269 + 1;
      _if__result_3382 = _M0L6_2atmpS2268 < _M0L3lenS617;
    } else {
      _if__result_3382 = 0;
    }
    if (_if__result_3382) {
      int32_t _M0L3valS2275 = _M0L5indexS619->$0;
      int32_t _M0L6_2atmpS2274 = _M0L3valS2275 + 1;
      int32_t _M0L2c2S621 = _M0L4selfS618[_M0L6_2atmpS2274];
      #line 261 "/home/developer/.moon/lib/core/builtin/string.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S621)) {
        int32_t _M0L6_2atmpS2272 = (int32_t)_M0L2c1S620;
        int32_t _M0L6_2atmpS2273 = (int32_t)_M0L2c2S621;
        int32_t _M0L1cS622;
        int32_t _M0L3valS2271;
        int32_t _M0L6_2atmpS2270;
        #line 262 "/home/developer/.moon/lib/core/builtin/string.mbt"
        _M0L1cS622
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS2272, _M0L6_2atmpS2273);
        _M0L3valS2271 = _M0L5indexS619->$0;
        _M0L6_2atmpS2270 = _M0L3valS2271 + 2;
        _M0L5indexS619->$0 = _M0L6_2atmpS2270;
        return _M0L1cS622;
      }
    }
    _M0L3valS2277 = _M0L5indexS619->$0;
    _M0L6_2atmpS2276 = _M0L3valS2277 + 1;
    _M0L5indexS619->$0 = _M0L6_2atmpS2276;
    #line 269 "/home/developer/.moon/lib/core/builtin/string.mbt"
    _M0L6_2atmpS2278 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S620);
    return _M0L6_2atmpS2278;
  } else {
    return -1;
  }
}

struct _M0TPC16string10StringView _M0MPC16string6String12trim_2einner(
  moonbit_string_t _M0L4selfS615,
  struct _M0TPC16string10StringView _M0L5charsS616
) {
  int32_t _M0L6_2atmpS2262;
  struct _M0TPC16string10StringView _M0L6_2atmpS2261;
  struct _M0TPC16string10StringView _result_3383;
  #line 794 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2262 = Moonbit_array_length(_M0L4selfS615);
  moonbit_incref(_M0L4selfS615);
  _M0L6_2atmpS2261
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS615, .$1 = 0, .$2 = _M0L6_2atmpS2262
  };
  #line 799 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3383
  = _M0MPC16string10StringView12trim_2einner(_M0L6_2atmpS2261, _M0L5charsS616);
  moonbit_decref(_M0L6_2atmpS2261.$0);
  return _result_3383;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView12trim_2einner(
  struct _M0TPC16string10StringView _M0L4selfS613,
  struct _M0TPC16string10StringView _M0L5charsS614
) {
  struct _M0TPC16string10StringView _M0L6_2atmpS2260;
  struct _M0TPC16string10StringView _result_3384;
  #line 783 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 788 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2260
  = _M0MPC16string10StringView19trim__start_2einner(_M0L4selfS613, _M0L5charsS614);
  #line 788 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3384
  = _M0MPC16string10StringView17trim__end_2einner(_M0L6_2atmpS2260, _M0L5charsS614);
  moonbit_decref(_M0L6_2atmpS2260.$0);
  return _result_3384;
}

struct _M0TPC16string10StringView _M0MPC16string6String17trim__end_2einner(
  moonbit_string_t _M0L4selfS611,
  struct _M0TPC16string10StringView _M0L5charsS612
) {
  int32_t _M0L6_2atmpS2259;
  struct _M0TPC16string10StringView _M0L6_2atmpS2258;
  struct _M0TPC16string10StringView _result_3385;
  #line 757 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2259 = Moonbit_array_length(_M0L4selfS611);
  moonbit_incref(_M0L4selfS611);
  _M0L6_2atmpS2258
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS611, .$1 = 0, .$2 = _M0L6_2atmpS2259
  };
  #line 762 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3385
  = _M0MPC16string10StringView17trim__end_2einner(_M0L6_2atmpS2258, _M0L5charsS612);
  moonbit_decref(_M0L6_2atmpS2258.$0);
  return _result_3385;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView17trim__end_2einner(
  struct _M0TPC16string10StringView _M0L4selfS610,
  struct _M0TPC16string10StringView _M0L5charsS609
) {
  struct _M0TPC16string10StringView _M0L1xS605;
  #line 734 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  moonbit_incref(_M0L4selfS610.$0);
  _M0L1xS605 = _M0L4selfS610;
  while (1) {
    int32_t _M0L3endS2241 = _M0L1xS605.$2;
    int32_t _M0L5startS2242 = _M0L1xS605.$1;
    int32_t _M0L6_2atmpS2240 = _M0L3endS2241 - _M0L5startS2242;
    if (_M0L6_2atmpS2240 == 0) {
      return _M0L1xS605;
    } else {
      moonbit_string_t _M0L3strS2251 = _M0L1xS605.$0;
      moonbit_string_t _M0L3strS2254 = _M0L1xS605.$0;
      int32_t _M0L5startS2255 = _M0L1xS605.$1;
      int32_t _M0L3endS2257 = _M0L1xS605.$2;
      int64_t _M0L6_2atmpS2256 = (int64_t)_M0L3endS2257;
      int64_t _M0L6_2atmpS2253;
      int32_t _M0L6_2atmpS2252;
      int32_t _M0L4_2acS607;
      moonbit_string_t _M0L3strS2243;
      int32_t _M0L5startS2244;
      moonbit_string_t _M0L3strS2247;
      int32_t _M0L5startS2248;
      int32_t _M0L3endS2250;
      int64_t _M0L6_2atmpS2249;
      int64_t _M0L6_2atmpS2246;
      int32_t _M0L6_2atmpS2245;
      struct _M0TPC16string10StringView _M0L4_2axS608;
      moonbit_incref(_M0L3strS2254);
      moonbit_incref(_M0L3strS2251);
      #line 739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L6_2atmpS2253
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2254, -1, _M0L5startS2255, _M0L6_2atmpS2256);
      moonbit_decref(_M0L3strS2254);
      _M0L6_2atmpS2252 = (int32_t)_M0L6_2atmpS2253;
      #line 739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L4_2acS607
      = _M0MPC16string6String16unsafe__char__at(_M0L3strS2251, _M0L6_2atmpS2252);
      moonbit_decref(_M0L3strS2251);
      _M0L3strS2243 = _M0L1xS605.$0;
      _M0L5startS2244 = _M0L1xS605.$1;
      _M0L3strS2247 = _M0L1xS605.$0;
      _M0L5startS2248 = _M0L1xS605.$1;
      _M0L3endS2250 = _M0L1xS605.$2;
      _M0L6_2atmpS2249 = (int64_t)_M0L3endS2250;
      moonbit_incref(_M0L3strS2247);
      moonbit_incref(_M0L3strS2243);
      #line 739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L6_2atmpS2246
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2247, -1, _M0L5startS2248, _M0L6_2atmpS2249);
      moonbit_decref(_M0L3strS2247);
      _M0L6_2atmpS2245 = (int32_t)_M0L6_2atmpS2246;
      _M0L4_2axS608
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L3strS2243, .$1 = _M0L5startS2244, .$2 = _M0L6_2atmpS2245
      };
      #line 743 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      if (
        _M0MPC16string10StringView14contains__char(_M0L5charsS609, _M0L4_2acS607)
      ) {
        moonbit_decref(_M0L1xS605.$0);
        _M0L1xS605 = _M0L4_2axS608;
        continue;
      } else {
        moonbit_decref(_M0L4_2axS608.$0);
        return _M0L1xS605;
      }
    }
    break;
  }
}

struct _M0TPC16string10StringView _M0MPC16string6String19trim__start_2einner(
  moonbit_string_t _M0L4selfS603,
  struct _M0TPC16string10StringView _M0L5charsS604
) {
  int32_t _M0L6_2atmpS2239;
  struct _M0TPC16string10StringView _M0L6_2atmpS2238;
  struct _M0TPC16string10StringView _result_3387;
  #line 707 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2239 = Moonbit_array_length(_M0L4selfS603);
  moonbit_incref(_M0L4selfS603);
  _M0L6_2atmpS2238
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS603, .$1 = 0, .$2 = _M0L6_2atmpS2239
  };
  #line 712 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3387
  = _M0MPC16string10StringView19trim__start_2einner(_M0L6_2atmpS2238, _M0L5charsS604);
  moonbit_decref(_M0L6_2atmpS2238.$0);
  return _result_3387;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView19trim__start_2einner(
  struct _M0TPC16string10StringView _M0L4selfS602,
  struct _M0TPC16string10StringView _M0L5charsS601
) {
  struct _M0TPC16string10StringView _M0L1xS596;
  #line 686 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  moonbit_incref(_M0L4selfS602.$0);
  _M0L1xS596 = _M0L4selfS602;
  while (1) {
    int32_t _M0L3endS2222 = _M0L1xS596.$2;
    int32_t _M0L5startS2223 = _M0L1xS596.$1;
    int32_t _M0L6_2atmpS2221 = _M0L3endS2222 - _M0L5startS2223;
    if (_M0L6_2atmpS2221 == 0) {
      return _M0L1xS596;
    } else {
      moonbit_string_t _M0L3strS2231 = _M0L1xS596.$0;
      moonbit_string_t _M0L3strS2234 = _M0L1xS596.$0;
      int32_t _M0L5startS2235 = _M0L1xS596.$1;
      int32_t _M0L3endS2237 = _M0L1xS596.$2;
      int64_t _M0L6_2atmpS2236 = (int64_t)_M0L3endS2237;
      int64_t _M0L6_2atmpS2233;
      int32_t _M0L6_2atmpS2232;
      int32_t _M0L4_2acS598;
      moonbit_string_t _M0L3strS2224;
      moonbit_string_t _M0L3strS2227;
      int32_t _M0L5startS2228;
      int32_t _M0L3endS2230;
      int64_t _M0L6_2atmpS2229;
      int64_t _M0L7_2abindS1216;
      int32_t _M0L6_2atmpS2225;
      int32_t _M0L3endS2226;
      struct _M0TPC16string10StringView _M0L4_2axS599;
      moonbit_incref(_M0L3strS2234);
      moonbit_incref(_M0L3strS2231);
      #line 691 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L6_2atmpS2233
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2234, 0, _M0L5startS2235, _M0L6_2atmpS2236);
      moonbit_decref(_M0L3strS2234);
      _M0L6_2atmpS2232 = (int32_t)_M0L6_2atmpS2233;
      #line 691 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L4_2acS598
      = _M0MPC16string6String16unsafe__char__at(_M0L3strS2231, _M0L6_2atmpS2232);
      moonbit_decref(_M0L3strS2231);
      _M0L3strS2224 = _M0L1xS596.$0;
      _M0L3strS2227 = _M0L1xS596.$0;
      _M0L5startS2228 = _M0L1xS596.$1;
      _M0L3endS2230 = _M0L1xS596.$2;
      _M0L6_2atmpS2229 = (int64_t)_M0L3endS2230;
      moonbit_incref(_M0L3strS2227);
      moonbit_incref(_M0L3strS2224);
      #line 691 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS1216
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2227, 1, _M0L5startS2228, _M0L6_2atmpS2229);
      moonbit_decref(_M0L3strS2227);
      if (_M0L7_2abindS1216 == 4294967296ll) {
        _M0L6_2atmpS2225 = _M0L1xS596.$2;
      } else {
        int64_t _M0L7_2aSomeS600 = _M0L7_2abindS1216;
        _M0L6_2atmpS2225 = (int32_t)_M0L7_2aSomeS600;
      }
      _M0L3endS2226 = _M0L1xS596.$2;
      _M0L4_2axS599
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L3strS2224, .$1 = _M0L6_2atmpS2225, .$2 = _M0L3endS2226
      };
      #line 695 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      if (
        _M0MPC16string10StringView14contains__char(_M0L5charsS601, _M0L4_2acS598)
      ) {
        moonbit_decref(_M0L1xS596.$0);
        _M0L1xS596 = _M0L4_2axS599;
        continue;
      } else {
        moonbit_decref(_M0L4_2axS599.$0);
        return _M0L1xS596;
      }
    }
    break;
  }
}

int32_t _M0MPC16string10StringView14contains__char(
  struct _M0TPC16string10StringView _M0L4selfS587,
  int32_t _M0L1cS589
) {
  int32_t _M0L3endS2219;
  int32_t _M0L5startS2220;
  int32_t _M0L3lenS586;
  #line 622 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2219 = _M0L4selfS587.$2;
  _M0L5startS2220 = _M0L4selfS587.$1;
  _M0L3lenS586 = _M0L3endS2219 - _M0L5startS2220;
  if (_M0L3lenS586 > 0) {
    int32_t _M0L1cS588 = _M0L1cS589;
    int32_t _if__result_3389;
    if (_M0L1cS588 >= 0) {
      _if__result_3389 = _M0L1cS588 <= 65535;
    } else {
      _if__result_3389 = 0;
    }
    if (_if__result_3389) {
      int32_t _M0L6_2atmpS2203 = (uint16_t)_M0L1cS588;
      #line 629 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0MPC16string10StringView20contains__code__unit(_M0L4selfS587, _M0L6_2atmpS2203);
    } else if (_M0L1cS588 < 0) {
      return 0;
    } else if (_M0L3lenS586 >= 2) {
      int32_t _M0L3adjS590 = _M0L1cS588 - 65536;
      int32_t _M0L6_2atmpS2218 = _M0L3adjS590 >> 10;
      int32_t _M0L4highS591 = 55296 + _M0L6_2atmpS2218;
      if (_M0L4highS591 <= 65535) {
        int32_t _M0L4highS592 = (uint16_t)_M0L4highS591;
        int32_t _M0L6_2atmpS2217 = _M0L3adjS590 & 1023;
        int32_t _M0L6_2atmpS2216 = 56320 + _M0L6_2atmpS2217;
        int32_t _M0L3lowS593 = (uint16_t)_M0L6_2atmpS2216;
        int32_t _M0L1iS594 = 0;
        while (1) {
          int32_t _M0L6_2atmpS2204 = _M0L3lenS586 - 1;
          if (_M0L1iS594 < _M0L6_2atmpS2204) {
            moonbit_string_t _M0L3strS2206 = _M0L4selfS587.$0;
            int32_t _M0L5startS2208 = _M0L4selfS587.$1;
            int32_t _M0L6_2atmpS2207 = _M0L5startS2208 + _M0L1iS594;
            int32_t _M0L6_2atmpS2205 = _M0L3strS2206[_M0L6_2atmpS2207];
            int32_t _M0L6_2atmpS2215;
            if (_M0L6_2atmpS2205 == _M0L4highS592) {
              moonbit_string_t _M0L3strS2210 = _M0L4selfS587.$0;
              int32_t _M0L5startS2212 = _M0L4selfS587.$1;
              int32_t _M0L6_2atmpS2213 = _M0L1iS594 + 1;
              int32_t _M0L6_2atmpS2211 = _M0L5startS2212 + _M0L6_2atmpS2213;
              int32_t _M0L6_2atmpS2209 = _M0L3strS2210[_M0L6_2atmpS2211];
              int32_t _M0L6_2atmpS2214;
              if (_M0L6_2atmpS2209 == _M0L3lowS593) {
                return 1;
              }
              _M0L6_2atmpS2214 = _M0L1iS594 + 2;
              _M0L1iS594 = _M0L6_2atmpS2214;
              continue;
            }
            _M0L6_2atmpS2215 = _M0L1iS594 + 1;
            _M0L1iS594 = _M0L6_2atmpS2215;
            continue;
          }
          break;
        }
      } else {
        return 0;
      }
    } else {
      return 0;
    }
    return 0;
  } else {
    return 0;
  }
}

int32_t _M0MPC16string6String8contains(
  moonbit_string_t _M0L4selfS584,
  struct _M0TPC16string10StringView _M0L3strS585
) {
  int32_t _M0L6_2atmpS2202;
  struct _M0TPC16string10StringView _M0L6_2atmpS2201;
  int32_t _result_3391;
  #line 545 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2202 = Moonbit_array_length(_M0L4selfS584);
  moonbit_incref(_M0L4selfS584);
  _M0L6_2atmpS2201
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS584, .$1 = 0, .$2 = _M0L6_2atmpS2202
  };
  #line 546 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3391
  = _M0MPC16string10StringView8contains(_M0L6_2atmpS2201, _M0L3strS585);
  moonbit_decref(_M0L6_2atmpS2201.$0);
  return _result_3391;
}

int32_t _M0MPC16string10StringView8contains(
  struct _M0TPC16string10StringView _M0L4selfS582,
  struct _M0TPC16string10StringView _M0L3strS581
) {
  int32_t _M0L3endS2199;
  int32_t _M0L5startS2200;
  int32_t _M0L7_2abindS580;
  #line 535 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2199 = _M0L3strS581.$2;
  _M0L5startS2200 = _M0L3strS581.$1;
  _M0L7_2abindS580 = _M0L3endS2199 - _M0L5startS2200;
  switch (_M0L7_2abindS580) {
    case 0: {
      return 1;
      break;
    }
    
    case 1: {
      moonbit_string_t _M0L3strS2196 = _M0L3strS581.$0;
      int32_t _M0L5startS2197 = _M0L3strS581.$1;
      int32_t _M0L6_2atmpS2195 = _M0L3strS2196[_M0L5startS2197];
      #line 538 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0MPC16string10StringView20contains__code__unit(_M0L4selfS582, _M0L6_2atmpS2195);
      break;
    }
    default: {
      int64_t _M0L7_2abindS583;
      int32_t _M0L6_2atmpS2198;
      #line 539 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS583
      = _M0MPC16string10StringView4find(_M0L4selfS582, _M0L3strS581);
      _M0L6_2atmpS2198 = _M0L7_2abindS583 == 4294967296ll;
      return !_M0L6_2atmpS2198;
      break;
    }
  }
}

int32_t _M0MPC16string10StringView20contains__code__unit(
  struct _M0TPC16string10StringView _M0L4selfS576,
  int32_t _M0L4codeS578
) {
  int32_t _M0L3endS2193;
  int32_t _M0L5startS2194;
  int32_t _M0L7_2abindS575;
  int32_t _M0L1iS577;
  #line 524 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2193 = _M0L4selfS576.$2;
  _M0L5startS2194 = _M0L4selfS576.$1;
  _M0L7_2abindS575 = _M0L3endS2193 - _M0L5startS2194;
  _M0L1iS577 = 0;
  while (1) {
    if (_M0L1iS577 < _M0L7_2abindS575) {
      moonbit_string_t _M0L3strS2189 = _M0L4selfS576.$0;
      int32_t _M0L5startS2191 = _M0L4selfS576.$1;
      int32_t _M0L6_2atmpS2190 = _M0L5startS2191 + _M0L1iS577;
      int32_t _M0L6_2atmpS2188 = _M0L3strS2189[_M0L6_2atmpS2190];
      int32_t _M0L6_2atmpS2192;
      if (_M0L6_2atmpS2188 == _M0L4codeS578) {
        return 1;
      }
      _M0L6_2atmpS2192 = _M0L1iS577 + 1;
      _M0L1iS577 = _M0L6_2atmpS2192;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS563,
  moonbit_string_t _M0L5valueS565
) {
  int32_t _M0L3lenS2168;
  moonbit_string_t* _M0L6_2atmpS2170;
  int32_t _M0L6_2atmpS2169;
  int32_t _M0L6lengthS564;
  moonbit_string_t* _M0L3bufS2171;
  moonbit_string_t _M0L6_2aoldS3141;
  int32_t _M0L6_2atmpS2172;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2168 = _M0L4selfS563->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2170 = _M0MPC15array5Array6bufferGsE(_M0L4selfS563);
  _M0L6_2atmpS2169 = Moonbit_array_length(_M0L6_2atmpS2170);
  moonbit_decref(_M0L6_2atmpS2170);
  if (_M0L3lenS2168 == _M0L6_2atmpS2169) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGsE(_M0L4selfS563);
  }
  _M0L6lengthS564 = _M0L4selfS563->$1;
  _M0L3bufS2171 = _M0L4selfS563->$0;
  _M0L6_2aoldS3141 = (moonbit_string_t)_M0L3bufS2171[_M0L6lengthS564];
  moonbit_incref(_M0L5valueS565);
  moonbit_decref(_M0L6_2aoldS3141);
  _M0L3bufS2171[_M0L6lengthS564] = _M0L5valueS565;
  _M0L6_2atmpS2172 = _M0L6lengthS564 + 1;
  _M0L4selfS563->$1 = _M0L6_2atmpS2172;
  return 0;
}

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS566,
  struct _M0TUsiE* _M0L5valueS568
) {
  int32_t _M0L3lenS2173;
  struct _M0TUsiE** _M0L6_2atmpS2175;
  int32_t _M0L6_2atmpS2174;
  int32_t _M0L6lengthS567;
  struct _M0TUsiE** _M0L3bufS2176;
  struct _M0TUsiE* _M0L6_2aoldS3143;
  int32_t _M0L6_2atmpS2177;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2173 = _M0L4selfS566->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2175 = _M0MPC15array5Array6bufferGUsiEE(_M0L4selfS566);
  _M0L6_2atmpS2174 = Moonbit_array_length(_M0L6_2atmpS2175);
  moonbit_decref(_M0L6_2atmpS2175);
  if (_M0L3lenS2173 == _M0L6_2atmpS2174) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGUsiEE(_M0L4selfS566);
  }
  _M0L6lengthS567 = _M0L4selfS566->$1;
  _M0L3bufS2176 = _M0L4selfS566->$0;
  _M0L6_2aoldS3143 = (struct _M0TUsiE*)_M0L3bufS2176[_M0L6lengthS567];
  moonbit_incref(_M0L5valueS568);
  if (_M0L6_2aoldS3143) {
    moonbit_decref(_M0L6_2aoldS3143);
  }
  _M0L3bufS2176[_M0L6lengthS567] = _M0L5valueS568;
  _M0L6_2atmpS2177 = _M0L6lengthS567 + 1;
  _M0L4selfS566->$1 = _M0L6_2atmpS2177;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L4selfS569,
  struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L5valueS571
) {
  int32_t _M0L3lenS2178;
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L6_2atmpS2180;
  int32_t _M0L6_2atmpS2179;
  int32_t _M0L6lengthS570;
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3bufS2181;
  struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L6_2aoldS3145;
  int32_t _M0L6_2atmpS2182;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2178 = _M0L4selfS569->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2180
  = _M0MPC15array5Array6bufferGRP46biolab8bio__seq3cmd4main4CaseE(_M0L4selfS569);
  _M0L6_2atmpS2179 = Moonbit_array_length(_M0L6_2atmpS2180);
  moonbit_decref(_M0L6_2atmpS2180);
  if (_M0L3lenS2178 == _M0L6_2atmpS2179) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRP46biolab8bio__seq3cmd4main4CaseE(_M0L4selfS569);
  }
  _M0L6lengthS570 = _M0L4selfS569->$1;
  _M0L3bufS2181 = _M0L4selfS569->$0;
  _M0L6_2aoldS3145
  = (struct _M0TP46biolab8bio__seq3cmd4main4Case*)_M0L3bufS2181[
      _M0L6lengthS570
    ];
  moonbit_incref(_M0L5valueS571);
  if (_M0L6_2aoldS3145) {
    moonbit_decref(_M0L6_2aoldS3145);
  }
  _M0L3bufS2181[_M0L6lengthS570] = _M0L5valueS571;
  _M0L6_2atmpS2182 = _M0L6lengthS570 + 1;
  _M0L4selfS569->$1 = _M0L6_2atmpS2182;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L4selfS572,
  struct _M0TP26biolab8bio__seq3Seq* _M0L5valueS574
) {
  int32_t _M0L3lenS2183;
  struct _M0TP26biolab8bio__seq3Seq** _M0L6_2atmpS2185;
  int32_t _M0L6_2atmpS2184;
  int32_t _M0L6lengthS573;
  struct _M0TP26biolab8bio__seq3Seq** _M0L3bufS2186;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2aoldS3147;
  int32_t _M0L6_2atmpS2187;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2183 = _M0L4selfS572->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2185
  = _M0MPC15array5Array6bufferGRP26biolab8bio__seq3SeqE(_M0L4selfS572);
  _M0L6_2atmpS2184 = Moonbit_array_length(_M0L6_2atmpS2185);
  moonbit_decref(_M0L6_2atmpS2185);
  if (_M0L3lenS2183 == _M0L6_2atmpS2184) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRP26biolab8bio__seq3SeqE(_M0L4selfS572);
  }
  _M0L6lengthS573 = _M0L4selfS572->$1;
  _M0L3bufS2186 = _M0L4selfS572->$0;
  _M0L6_2aoldS3147
  = (struct _M0TP26biolab8bio__seq3Seq*)_M0L3bufS2186[_M0L6lengthS573];
  moonbit_incref(_M0L5valueS574);
  if (_M0L6_2aoldS3147) {
    moonbit_decref(_M0L6_2aoldS3147);
  }
  _M0L3bufS2186[_M0L6lengthS573] = _M0L5valueS574;
  _M0L6_2atmpS2187 = _M0L6lengthS573 + 1;
  _M0L4selfS572->$1 = _M0L6_2atmpS2187;
  return 0;
}

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE* _M0L4selfS552) {
  int32_t _M0L8old__capS551;
  int32_t _M0L8new__capS553;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS551 = _M0L4selfS552->$1;
  if (_M0L8old__capS551 == 0) {
    _M0L8new__capS553 = 8;
  } else {
    _M0L8new__capS553 = _M0L8old__capS551 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGsE(_M0L4selfS552, _M0L8new__capS553);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS555
) {
  int32_t _M0L8old__capS554;
  int32_t _M0L8new__capS556;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS554 = _M0L4selfS555->$1;
  if (_M0L8old__capS554 == 0) {
    _M0L8new__capS556 = 8;
  } else {
    _M0L8new__capS556 = _M0L8old__capS554 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGUsiEE(_M0L4selfS555, _M0L8new__capS556);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L4selfS558
) {
  int32_t _M0L8old__capS557;
  int32_t _M0L8new__capS559;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS557 = _M0L4selfS558->$1;
  if (_M0L8old__capS557 == 0) {
    _M0L8new__capS559 = 8;
  } else {
    _M0L8new__capS559 = _M0L8old__capS557 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRP46biolab8bio__seq3cmd4main4CaseE(_M0L4selfS558, _M0L8new__capS559);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L4selfS561
) {
  int32_t _M0L8old__capS560;
  int32_t _M0L8new__capS562;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS560 = _M0L4selfS561->$1;
  if (_M0L8old__capS560 == 0) {
    _M0L8new__capS562 = 8;
  } else {
    _M0L8new__capS562 = _M0L8old__capS560 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq3SeqE(_M0L4selfS561, _M0L8new__capS562);
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS528,
  int32_t _M0L13new__capacityS531
) {
  moonbit_string_t* _M0L8old__bufS527;
  int32_t _M0L8old__capS529;
  int32_t _M0L9copy__lenS530;
  moonbit_string_t* _M0L8new__bufS532;
  moonbit_string_t* _M0L6_2aoldS3149;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS527 = _M0L4selfS528->$0;
  _M0L8old__capS529 = Moonbit_array_length(_M0L8old__bufS527);
  if (_M0L8old__capS529 < _M0L13new__capacityS531) {
    _M0L9copy__lenS530 = _M0L8old__capS529;
  } else {
    _M0L9copy__lenS530 = _M0L13new__capacityS531;
  }
  moonbit_incref(_M0L8old__bufS527);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS532
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(_M0L8old__bufS527, _M0L13new__capacityS531, _M0L9copy__lenS530, 0, 0);
  moonbit_decref(_M0L8old__bufS527);
  _M0L6_2aoldS3149 = _M0L4selfS528->$0;
  moonbit_decref(_M0L6_2aoldS3149);
  _M0L4selfS528->$0 = _M0L8new__bufS532;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS534,
  int32_t _M0L13new__capacityS537
) {
  struct _M0TUsiE** _M0L8old__bufS533;
  int32_t _M0L8old__capS535;
  int32_t _M0L9copy__lenS536;
  struct _M0TUsiE** _M0L8new__bufS538;
  struct _M0TUsiE** _M0L6_2aoldS3151;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS533 = _M0L4selfS534->$0;
  _M0L8old__capS535 = Moonbit_array_length(_M0L8old__bufS533);
  if (_M0L8old__capS535 < _M0L13new__capacityS537) {
    _M0L9copy__lenS536 = _M0L8old__capS535;
  } else {
    _M0L9copy__lenS536 = _M0L13new__capacityS537;
  }
  moonbit_incref(_M0L8old__bufS533);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS538
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(_M0L8old__bufS533, _M0L13new__capacityS537, _M0L9copy__lenS536, 0, 0);
  moonbit_decref(_M0L8old__bufS533);
  _M0L6_2aoldS3151 = _M0L4selfS534->$0;
  moonbit_decref(_M0L6_2aoldS3151);
  _M0L4selfS534->$0 = _M0L8new__bufS538;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L4selfS540,
  int32_t _M0L13new__capacityS543
) {
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L8old__bufS539;
  int32_t _M0L8old__capS541;
  int32_t _M0L9copy__lenS542;
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L8new__bufS544;
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L6_2aoldS3153;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS539 = _M0L4selfS540->$0;
  _M0L8old__capS541 = Moonbit_array_length(_M0L8old__bufS539);
  if (_M0L8old__capS541 < _M0L13new__capacityS543) {
    _M0L9copy__lenS542 = _M0L8old__capS541;
  } else {
    _M0L9copy__lenS542 = _M0L13new__capacityS543;
  }
  moonbit_incref(_M0L8old__bufS539);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS544
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRP46biolab8bio__seq3cmd4main4CaseE(_M0L8old__bufS539, _M0L13new__capacityS543, _M0L9copy__lenS542, 0, 0);
  moonbit_decref(_M0L8old__bufS539);
  _M0L6_2aoldS3153 = _M0L4selfS540->$0;
  moonbit_decref(_M0L6_2aoldS3153);
  _M0L4selfS540->$0 = _M0L8new__bufS544;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L4selfS546,
  int32_t _M0L13new__capacityS549
) {
  struct _M0TP26biolab8bio__seq3Seq** _M0L8old__bufS545;
  int32_t _M0L8old__capS547;
  int32_t _M0L9copy__lenS548;
  struct _M0TP26biolab8bio__seq3Seq** _M0L8new__bufS550;
  struct _M0TP26biolab8bio__seq3Seq** _M0L6_2aoldS3155;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS545 = _M0L4selfS546->$0;
  _M0L8old__capS547 = Moonbit_array_length(_M0L8old__bufS545);
  if (_M0L8old__capS547 < _M0L13new__capacityS549) {
    _M0L9copy__lenS548 = _M0L8old__capS547;
  } else {
    _M0L9copy__lenS548 = _M0L13new__capacityS549;
  }
  moonbit_incref(_M0L8old__bufS545);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS550
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq3SeqE(_M0L8old__bufS545, _M0L13new__capacityS549, _M0L9copy__lenS548, 0, 0);
  moonbit_decref(_M0L8old__bufS545);
  _M0L6_2aoldS3155 = _M0L4selfS546->$0;
  moonbit_decref(_M0L6_2aoldS3155);
  _M0L4selfS546->$0 = _M0L8new__bufS550;
  return 0;
}

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L4selfS526
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS526->$1;
}

struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0MPC15array5Array11new_2einnerGRP46biolab8bio__seq3cmd4main4CaseE(
  int32_t _M0L8capacityS523
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS523 == 0) {
    struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L6_2atmpS2162 =
      (struct _M0TP46biolab8bio__seq3cmd4main4Case**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _block_3393 =
      (struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE));
    Moonbit_object_header(_block_3393)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 62, 0);
    _block_3393->$0 = _M0L6_2atmpS2162;
    _block_3393->$1 = 0;
    return _block_3393;
  } else {
    struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L6_2atmpS2163 =
      (struct _M0TP46biolab8bio__seq3cmd4main4Case**)moonbit_make_ref_array(_M0L8capacityS523, 0);
    struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _block_3394 =
      (struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE));
    Moonbit_object_header(_block_3394)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 62, 0);
    _block_3394->$0 = _M0L6_2atmpS2163;
    _block_3394->$1 = 0;
    return _block_3394;
  }
}

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(
  int32_t _M0L8capacityS524
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS524 == 0) {
    moonbit_string_t* _M0L6_2atmpS2164 =
      (moonbit_string_t*)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGsE* _block_3395 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_3395)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
    _block_3395->$0 = _M0L6_2atmpS2164;
    _block_3395->$1 = 0;
    return _block_3395;
  } else {
    moonbit_string_t* _M0L6_2atmpS2165 =
      (moonbit_string_t*)moonbit_make_ref_array(_M0L8capacityS524, (moonbit_string_t)moonbit_string_literal_0.data);
    struct _M0TPB5ArrayGsE* _block_3396 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_3396)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
    _block_3396->$0 = _M0L6_2atmpS2165;
    _block_3396->$1 = 0;
    return _block_3396;
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq3SeqE(
  int32_t _M0L8capacityS525
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS525 == 0) {
    struct _M0TP26biolab8bio__seq3Seq** _M0L6_2atmpS2166 =
      (struct _M0TP26biolab8bio__seq3Seq**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _block_3397 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE));
    Moonbit_object_header(_block_3397)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
    _block_3397->$0 = _M0L6_2atmpS2166;
    _block_3397->$1 = 0;
    return _block_3397;
  } else {
    struct _M0TP26biolab8bio__seq3Seq** _M0L6_2atmpS2167 =
      (struct _M0TP26biolab8bio__seq3Seq**)moonbit_make_ref_array(_M0L8capacityS525, 0);
    struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _block_3398 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE));
    Moonbit_object_header(_block_3398)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
    _block_3398->$0 = _M0L6_2atmpS2167;
    _block_3398->$1 = 0;
    return _block_3398;
  }
}

int64_t _M0MPC16string6String9rev__find(
  moonbit_string_t _M0L4selfS521,
  struct _M0TPC16string10StringView _M0L3strS522
) {
  int32_t _M0L6_2atmpS2161;
  struct _M0TPC16string10StringView _M0L6_2atmpS2160;
  int64_t _result_3399;
  #line 233 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2161 = Moonbit_array_length(_M0L4selfS521);
  moonbit_incref(_M0L4selfS521);
  _M0L6_2atmpS2160
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS521, .$1 = 0, .$2 = _M0L6_2atmpS2161
  };
  #line 234 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3399
  = _M0MPC16string10StringView9rev__find(_M0L6_2atmpS2160, _M0L3strS522);
  moonbit_decref(_M0L6_2atmpS2160.$0);
  return _result_3399;
}

int64_t _M0MPC16string10StringView9rev__find(
  struct _M0TPC16string10StringView _M0L4selfS520,
  struct _M0TPC16string10StringView _M0L3strS519
) {
  int32_t _M0L3endS2158;
  int32_t _M0L5startS2159;
  int32_t _M0L6_2atmpS2157;
  #line 165 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2158 = _M0L3strS519.$2;
  _M0L5startS2159 = _M0L3strS519.$1;
  _M0L6_2atmpS2157 = _M0L3endS2158 - _M0L5startS2159;
  if (_M0L6_2atmpS2157 <= 4) {
    #line 167 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB23brute__force__rev__find(_M0L4selfS520, _M0L3strS519);
  } else {
    #line 169 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB33boyer__moore__horspool__rev__find(_M0L4selfS520, _M0L3strS519);
  }
}

int64_t _M0FPB23brute__force__rev__find(
  struct _M0TPC16string10StringView _M0L8haystackS509,
  struct _M0TPC16string10StringView _M0L6needleS511
) {
  int32_t _M0L3endS2155;
  int32_t _M0L5startS2156;
  int32_t _M0L13haystack__lenS508;
  int32_t _M0L3endS2153;
  int32_t _M0L5startS2154;
  int32_t _M0L11needle__lenS510;
  #line 178 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2155 = _M0L8haystackS509.$2;
  _M0L5startS2156 = _M0L8haystackS509.$1;
  _M0L13haystack__lenS508 = _M0L3endS2155 - _M0L5startS2156;
  _M0L3endS2153 = _M0L6needleS511.$2;
  _M0L5startS2154 = _M0L6needleS511.$1;
  _M0L11needle__lenS510 = _M0L3endS2153 - _M0L5startS2154;
  if (_M0L11needle__lenS510 > 0) {
    if (_M0L13haystack__lenS508 >= _M0L11needle__lenS510) {
      moonbit_string_t _M0L3strS2151 = _M0L6needleS511.$0;
      int32_t _M0L5startS2152 = _M0L6needleS511.$1;
      int32_t _M0L13needle__firstS512 = _M0L3strS2151[_M0L5startS2152];
      int32_t _M0L7_2abindS513 =
        _M0L13haystack__lenS508 - _M0L11needle__lenS510;
      int32_t _M0L1iS514 = _M0L7_2abindS513;
      while (1) {
        if (_M0L1iS514 >= 0) {
          moonbit_string_t _M0L3strS2138 = _M0L8haystackS509.$0;
          int32_t _M0L5startS2140 = _M0L8haystackS509.$1;
          int32_t _M0L6_2atmpS2139 = _M0L5startS2140 + _M0L1iS514;
          int32_t _M0L6_2atmpS2137 = _M0L3strS2138[_M0L6_2atmpS2139];
          int32_t _M0L1jS517;
          int32_t _M0L6_2atmpS2136;
          if (_M0L6_2atmpS2137 != _M0L13needle__firstS512) {
            goto join_515;
          }
          _M0L1jS517 = 1;
          while (1) {
            if (_M0L1jS517 < _M0L11needle__lenS510) {
              moonbit_string_t _M0L3strS2146 = _M0L8haystackS509.$0;
              int32_t _M0L5startS2148 = _M0L8haystackS509.$1;
              int32_t _M0L6_2atmpS2149 = _M0L1iS514 + _M0L1jS517;
              int32_t _M0L6_2atmpS2147 = _M0L5startS2148 + _M0L6_2atmpS2149;
              int32_t _M0L6_2atmpS2141 = _M0L3strS2146[_M0L6_2atmpS2147];
              moonbit_string_t _M0L3strS2143 = _M0L6needleS511.$0;
              int32_t _M0L5startS2145 = _M0L6needleS511.$1;
              int32_t _M0L6_2atmpS2144 = _M0L5startS2145 + _M0L1jS517;
              int32_t _M0L6_2atmpS2142 = _M0L3strS2143[_M0L6_2atmpS2144];
              int32_t _M0L6_2atmpS2150;
              if (_M0L6_2atmpS2141 != _M0L6_2atmpS2142) {
                break;
              }
              _M0L6_2atmpS2150 = _M0L1jS517 + 1;
              _M0L1jS517 = _M0L6_2atmpS2150;
              continue;
            } else {
              return (int64_t)_M0L1iS514;
            }
            break;
          }
          goto join_515;
          goto joinlet_3401;
          join_515:;
          _M0L6_2atmpS2136 = _M0L1iS514 - 1;
          _M0L1iS514 = _M0L6_2atmpS2136;
          continue;
          joinlet_3401:;
        }
        break;
      }
      return 4294967296ll;
    } else {
      return 4294967296ll;
    }
  } else {
    return (int64_t)_M0L13haystack__lenS508;
  }
}

int64_t _M0MPC16string6String8find__by(
  moonbit_string_t _M0L4selfS506,
  struct _M0TWuEu* _M0L4predS507
) {
  int32_t _M0L6_2atmpS2135;
  struct _M0TPC16string10StringView _M0L6_2atmpS2134;
  int64_t _result_3403;
  #line 144 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2135 = Moonbit_array_length(_M0L4selfS506);
  moonbit_incref(_M0L4selfS506);
  _M0L6_2atmpS2134
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS506, .$1 = 0, .$2 = _M0L6_2atmpS2135
  };
  #line 145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3403
  = _M0MPC16string10StringView8find__by(_M0L6_2atmpS2134, _M0L4predS507);
  moonbit_decref(_M0L6_2atmpS2134.$0);
  return _result_3403;
}

int64_t _M0MPC16string10StringView8find__by(
  struct _M0TPC16string10StringView _M0L4selfS493,
  struct _M0TWuEu* _M0L4predS502
) {
  moonbit_string_t _M0L7_2abindS492;
  int32_t _M0L7_2abindS494;
  int32_t _M0L7_2abindS495;
  int32_t _M0L16_2astring__indexS496;
  int32_t _M0L1iS497;
  #line 131 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS492 = _M0L4selfS493.$0;
  _M0L7_2abindS494 = _M0L4selfS493.$1;
  _M0L7_2abindS495 = _M0L4selfS493.$2;
  moonbit_incref(_M0L7_2abindS492);
  _M0L16_2astring__indexS496 = _M0L7_2abindS494;
  _M0L1iS497 = 0;
  while (1) {
    if (_M0L16_2astring__indexS496 < _M0L7_2abindS495) {
      int32_t _M0L31_2adecoded__next__string__indexS499;
      int32_t _M0L16_2adecoded__charS500;
      int32_t _M0L7_2abindS504 = _M0L7_2abindS492[_M0L16_2astring__indexS496];
      int32_t _M0L6_2atmpS2115 = (int32_t)_M0L7_2abindS504;
      int32_t _if__result_3406;
      int32_t _if__result_3407;
      int32_t _M0L20_2anext__char__indexS501;
      if (_M0L6_2atmpS2115 >= 55296) {
        int32_t _M0L6_2atmpS2114 = (int32_t)_M0L7_2abindS504;
        _if__result_3406 = _M0L6_2atmpS2114 <= 56319;
      } else {
        _if__result_3406 = 0;
      }
      if (_if__result_3406) {
        int32_t _M0L6_2atmpS2113 = _M0L16_2astring__indexS496 + 1;
        _if__result_3407 = _M0L6_2atmpS2113 < _M0L7_2abindS495;
      } else {
        _if__result_3407 = 0;
      }
      if (_if__result_3407) {
        int32_t _M0L6_2atmpS2130 = _M0L16_2astring__indexS496 + 1;
        int32_t _M0L7_2abindS505 = _M0L7_2abindS492[_M0L6_2atmpS2130];
        int32_t _M0L6_2atmpS2117 = (int32_t)_M0L7_2abindS505;
        int32_t _if__result_3408;
        if (_M0L6_2atmpS2117 >= 56320) {
          int32_t _M0L6_2atmpS2116 = (int32_t)_M0L7_2abindS505;
          _if__result_3408 = _M0L6_2atmpS2116 <= 57343;
        } else {
          _if__result_3408 = 0;
        }
        if (_if__result_3408) {
          int32_t _M0L6_2atmpS2118 = _M0L16_2astring__indexS496 + 2;
          int32_t _M0L6_2atmpS2126 = (int32_t)_M0L7_2abindS504;
          int32_t _M0L6_2atmpS2125 = _M0L6_2atmpS2126 - 55296;
          int32_t _M0L6_2atmpS2123 = _M0L6_2atmpS2125 * 1024;
          int32_t _M0L6_2atmpS2124 = (int32_t)_M0L7_2abindS505;
          int32_t _M0L6_2atmpS2122 = _M0L6_2atmpS2123 + _M0L6_2atmpS2124;
          int32_t _M0L6_2atmpS2121 = _M0L6_2atmpS2122 - 56320;
          int32_t _M0L6_2atmpS2120 = _M0L6_2atmpS2121 + 65536;
          int32_t _M0L6_2atmpS2119;
          #line 133 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L6_2atmpS2119
          = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2120);
          _M0L31_2adecoded__next__string__indexS499 = _M0L6_2atmpS2118;
          _M0L16_2adecoded__charS500 = _M0L6_2atmpS2119;
          goto join_498;
        } else {
          int32_t _M0L6_2atmpS2127 = _M0L16_2astring__indexS496 + 1;
          int32_t _M0L6_2atmpS2129 = (int32_t)_M0L7_2abindS504;
          int32_t _M0L6_2atmpS2128;
          #line 133 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L6_2atmpS2128
          = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2129);
          _M0L31_2adecoded__next__string__indexS499 = _M0L6_2atmpS2127;
          _M0L16_2adecoded__charS500 = _M0L6_2atmpS2128;
          goto join_498;
        }
      } else {
        int32_t _M0L6_2atmpS2131 = _M0L16_2astring__indexS496 + 1;
        int32_t _M0L6_2atmpS2133 = (int32_t)_M0L7_2abindS504;
        int32_t _M0L6_2atmpS2132;
        #line 133 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS2132
        = _M0MPC13int3Int16unsafe__to__char(_M0L6_2atmpS2133);
        _M0L31_2adecoded__next__string__indexS499 = _M0L6_2atmpS2131;
        _M0L16_2adecoded__charS500 = _M0L6_2atmpS2132;
        goto join_498;
      }
      goto joinlet_3405;
      join_498:;
      _M0L20_2anext__char__indexS501 = _M0L1iS497 + 1;
      #line 134 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      if (_M0L4predS502->code(_M0L4predS502, _M0L16_2adecoded__charS500)) {
        moonbit_decref(_M0L7_2abindS492);
        return (int64_t)_M0L1iS497;
      }
      _M0L16_2astring__indexS496 = _M0L31_2adecoded__next__string__indexS499;
      _M0L1iS497 = _M0L20_2anext__char__indexS501;
      continue;
      joinlet_3405:;
    } else {
      moonbit_decref(_M0L7_2abindS492);
    }
    break;
  }
  return 4294967296ll;
}

moonbit_string_t _M0MPC16string6String6repeat(
  moonbit_string_t _M0L4selfS485,
  int32_t _M0L1nS484
) {
  #line 1049 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  if (_M0L1nS484 < 0) {
    #line 1051 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_107.data);
  } else if (_M0L1nS484 == 0) {
    return (moonbit_string_t)moonbit_string_literal_0.data;
  } else if (_M0L1nS484 == 1) {
    moonbit_incref(_M0L4selfS485);
    return _M0L4selfS485;
  } else {
    int32_t _M0L3lenS486 = Moonbit_array_length(_M0L4selfS485);
    int32_t _M0L5totalS487 = _M0L3lenS486 * _M0L1nS484;
    int32_t _if__result_3409;
    if (_M0L3lenS486 == 0) {
      _if__result_3409 = 1;
    } else {
      int32_t _M0L6_2atmpS2111 = _M0L5totalS487 / _M0L1nS484;
      _if__result_3409 = _M0L6_2atmpS2111 == _M0L3lenS486;
    }
    if (_if__result_3409) {
      struct _M0TPB13StringBuilder* _M0L3bufS488;
      moonbit_string_t _M0L3strS489;
      int32_t _M0L2__S490;
      moonbit_string_t _result_3411;
      #line 1060 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3bufS488
      = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L5totalS487);
      #line 1061 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3strS489 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS485);
      _M0L2__S490 = 0;
      while (1) {
        if (_M0L2__S490 < _M0L1nS484) {
          int32_t _M0L6_2atmpS2112;
          #line 1063 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS488, _M0L3strS489);
          _M0L6_2atmpS2112 = _M0L2__S490 + 1;
          _M0L2__S490 = _M0L6_2atmpS2112;
          continue;
        } else {
          moonbit_decref(_M0L3strS489);
        }
        break;
      }
      #line 1065 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _result_3411 = _M0MPB13StringBuilder10to__string(_M0L3bufS488);
      moonbit_decref(_M0L3bufS488);
      return _result_3411;
    } else {
      #line 1058 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_108.data);
    }
  }
}

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(
  moonbit_string_t _M0L4selfS483
) {
  #line 222 "/home/developer/.moon/lib/core/builtin/show.mbt"
  moonbit_incref(_M0L4selfS483);
  return _M0L4selfS483;
}

int64_t _M0MPC16string6String4find(
  moonbit_string_t _M0L4selfS481,
  struct _M0TPC16string10StringView _M0L3strS482
) {
  int32_t _M0L6_2atmpS2110;
  struct _M0TPC16string10StringView _M0L6_2atmpS2109;
  int64_t _result_3412;
  #line 101 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2110 = Moonbit_array_length(_M0L4selfS481);
  moonbit_incref(_M0L4selfS481);
  _M0L6_2atmpS2109
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS481, .$1 = 0, .$2 = _M0L6_2atmpS2110
  };
  #line 102 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3412
  = _M0MPC16string10StringView4find(_M0L6_2atmpS2109, _M0L3strS482);
  moonbit_decref(_M0L6_2atmpS2109.$0);
  return _result_3412;
}

int64_t _M0FPB33boyer__moore__horspool__rev__find(
  struct _M0TPC16string10StringView _M0L8haystackS471,
  struct _M0TPC16string10StringView _M0L6needleS473
) {
  int32_t _M0L3endS2107;
  int32_t _M0L5startS2108;
  int32_t _M0L13haystack__lenS470;
  int32_t _M0L3endS2105;
  int32_t _M0L5startS2106;
  int32_t _M0L11needle__lenS472;
  #line 203 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2107 = _M0L8haystackS471.$2;
  _M0L5startS2108 = _M0L8haystackS471.$1;
  _M0L13haystack__lenS470 = _M0L3endS2107 - _M0L5startS2108;
  _M0L3endS2105 = _M0L6needleS473.$2;
  _M0L5startS2106 = _M0L6needleS473.$1;
  _M0L11needle__lenS472 = _M0L3endS2105 - _M0L5startS2106;
  if (_M0L11needle__lenS472 > 0) {
    if (_M0L13haystack__lenS470 >= _M0L11needle__lenS472) {
      int32_t* _M0L11skip__tableS474 =
        (int32_t*)moonbit_make_int32_array(256, _M0L11needle__lenS472);
      int32_t _M0L6_2atmpS2085 = _M0L11needle__lenS472 - 1;
      int32_t _M0L1iS475 = _M0L6_2atmpS2085;
      int32_t _M0L6_2atmpS2104;
      int32_t _M0L1iS477;
      while (1) {
        if (_M0L1iS475 >= 1) {
          moonbit_string_t _M0L3strS2081 = _M0L6needleS473.$0;
          int32_t _M0L5startS2083 = _M0L6needleS473.$1;
          int32_t _M0L6_2atmpS2082 = _M0L5startS2083 + _M0L1iS475;
          int32_t _M0L6_2atmpS2080 = _M0L3strS2081[_M0L6_2atmpS2082];
          int32_t _M0L6_2atmpS2079 = (int32_t)_M0L6_2atmpS2080;
          int32_t _M0L6_2atmpS2078 = _M0L6_2atmpS2079 & 255;
          int32_t _M0L6_2atmpS2084;
          if (
            _M0L6_2atmpS2078 < 0
            || _M0L6_2atmpS2078
               >= Moonbit_array_length(_M0L11skip__tableS474)
          ) {
            #line 213 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L11skip__tableS474[_M0L6_2atmpS2078] = _M0L1iS475;
          _M0L6_2atmpS2084 = _M0L1iS475 - 1;
          _M0L1iS475 = _M0L6_2atmpS2084;
          continue;
        }
        break;
      }
      _M0L6_2atmpS2104 = _M0L13haystack__lenS470 - _M0L11needle__lenS472;
      _M0L1iS477 = _M0L6_2atmpS2104;
      while (1) {
        if (_M0L1iS477 >= 0) {
          int32_t _M0L1jS478 = 0;
          moonbit_string_t _M0L3strS2101;
          int32_t _M0L5startS2103;
          int32_t _M0L6_2atmpS2102;
          int32_t _M0L6_2atmpS2100;
          int32_t _M0L6_2atmpS2099;
          int32_t _M0L6_2atmpS2098;
          int32_t _M0L6_2atmpS2097;
          int32_t _M0L6_2atmpS2096;
          while (1) {
            if (_M0L1jS478 < _M0L11needle__lenS472) {
              moonbit_string_t _M0L3strS2091 = _M0L8haystackS471.$0;
              int32_t _M0L5startS2093 = _M0L8haystackS471.$1;
              int32_t _M0L6_2atmpS2094 = _M0L1iS477 + _M0L1jS478;
              int32_t _M0L6_2atmpS2092 = _M0L5startS2093 + _M0L6_2atmpS2094;
              int32_t _M0L6_2atmpS2086 = _M0L3strS2091[_M0L6_2atmpS2092];
              moonbit_string_t _M0L3strS2088 = _M0L6needleS473.$0;
              int32_t _M0L5startS2090 = _M0L6needleS473.$1;
              int32_t _M0L6_2atmpS2089 = _M0L5startS2090 + _M0L1jS478;
              int32_t _M0L6_2atmpS2087 = _M0L3strS2088[_M0L6_2atmpS2089];
              int32_t _M0L6_2atmpS2095;
              if (_M0L6_2atmpS2086 != _M0L6_2atmpS2087) {
                break;
              }
              _M0L6_2atmpS2095 = _M0L1jS478 + 1;
              _M0L1jS478 = _M0L6_2atmpS2095;
              continue;
            } else {
              moonbit_decref(_M0L11skip__tableS474);
              return (int64_t)_M0L1iS477;
            }
            break;
          }
          _M0L3strS2101 = _M0L8haystackS471.$0;
          _M0L5startS2103 = _M0L8haystackS471.$1;
          _M0L6_2atmpS2102 = _M0L5startS2103 + _M0L1iS477;
          _M0L6_2atmpS2100 = _M0L3strS2101[_M0L6_2atmpS2102];
          _M0L6_2atmpS2099 = (int32_t)_M0L6_2atmpS2100;
          _M0L6_2atmpS2098 = _M0L6_2atmpS2099 & 255;
          if (
            _M0L6_2atmpS2098 < 0
            || _M0L6_2atmpS2098
               >= Moonbit_array_length(_M0L11skip__tableS474)
          ) {
            #line 217 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2097 = (int32_t)_M0L11skip__tableS474[_M0L6_2atmpS2098];
          _M0L6_2atmpS2096 = _M0L1iS477 - _M0L6_2atmpS2097;
          _M0L1iS477 = _M0L6_2atmpS2096;
          continue;
        } else {
          moonbit_decref(_M0L11skip__tableS474);
        }
        break;
      }
      return 4294967296ll;
    } else {
      return 4294967296ll;
    }
  } else {
    return (int64_t)_M0L13haystack__lenS470;
  }
}

int64_t _M0MPC16string10StringView4find(
  struct _M0TPC16string10StringView _M0L4selfS469,
  struct _M0TPC16string10StringView _M0L3strS468
) {
  int32_t _M0L3endS2076;
  int32_t _M0L5startS2077;
  int32_t _M0L6_2atmpS2075;
  #line 18 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2076 = _M0L3strS468.$2;
  _M0L5startS2077 = _M0L3strS468.$1;
  _M0L6_2atmpS2075 = _M0L3endS2076 - _M0L5startS2077;
  if (_M0L6_2atmpS2075 <= 4) {
    #line 20 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB18brute__force__find(_M0L4selfS469, _M0L3strS468);
  } else {
    #line 22 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB28boyer__moore__horspool__find(_M0L4selfS469, _M0L3strS468);
  }
}

int64_t _M0FPB18brute__force__find(
  struct _M0TPC16string10StringView _M0L8haystackS458,
  struct _M0TPC16string10StringView _M0L6needleS460
) {
  int32_t _M0L3endS2073;
  int32_t _M0L5startS2074;
  int32_t _M0L13haystack__lenS457;
  int32_t _M0L3endS2071;
  int32_t _M0L5startS2072;
  int32_t _M0L11needle__lenS459;
  #line 31 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2073 = _M0L8haystackS458.$2;
  _M0L5startS2074 = _M0L8haystackS458.$1;
  _M0L13haystack__lenS457 = _M0L3endS2073 - _M0L5startS2074;
  _M0L3endS2071 = _M0L6needleS460.$2;
  _M0L5startS2072 = _M0L6needleS460.$1;
  _M0L11needle__lenS459 = _M0L3endS2071 - _M0L5startS2072;
  if (_M0L11needle__lenS459 > 0) {
    if (_M0L13haystack__lenS457 >= _M0L11needle__lenS459) {
      moonbit_string_t _M0L3strS2069 = _M0L6needleS460.$0;
      int32_t _M0L5startS2070 = _M0L6needleS460.$1;
      int32_t _M0L13needle__firstS461 = _M0L3strS2069[_M0L5startS2070];
      int32_t _M0L12forward__lenS462 =
        _M0L13haystack__lenS457 - _M0L11needle__lenS459;
      int32_t _M0L1iS463 = 0;
      while (1) {
        if (_M0L1iS463 <= _M0L12forward__lenS462) {
          moonbit_string_t _M0L3strS2056 = _M0L8haystackS458.$0;
          int32_t _M0L5startS2058 = _M0L8haystackS458.$1;
          int32_t _M0L6_2atmpS2057 = _M0L5startS2058 + _M0L1iS463;
          int32_t _M0L6_2atmpS2055 = _M0L3strS2056[_M0L6_2atmpS2057];
          int32_t _M0L1jS466;
          int32_t _M0L6_2atmpS2054;
          if (_M0L6_2atmpS2055 != _M0L13needle__firstS461) {
            goto join_464;
          }
          _M0L1jS466 = 1;
          while (1) {
            if (_M0L1jS466 < _M0L11needle__lenS459) {
              moonbit_string_t _M0L3strS2064 = _M0L8haystackS458.$0;
              int32_t _M0L5startS2066 = _M0L8haystackS458.$1;
              int32_t _M0L6_2atmpS2067 = _M0L1iS463 + _M0L1jS466;
              int32_t _M0L6_2atmpS2065 = _M0L5startS2066 + _M0L6_2atmpS2067;
              int32_t _M0L6_2atmpS2059 = _M0L3strS2064[_M0L6_2atmpS2065];
              moonbit_string_t _M0L3strS2061 = _M0L6needleS460.$0;
              int32_t _M0L5startS2063 = _M0L6needleS460.$1;
              int32_t _M0L6_2atmpS2062 = _M0L5startS2063 + _M0L1jS466;
              int32_t _M0L6_2atmpS2060 = _M0L3strS2061[_M0L6_2atmpS2062];
              int32_t _M0L6_2atmpS2068;
              if (_M0L6_2atmpS2059 != _M0L6_2atmpS2060) {
                break;
              }
              _M0L6_2atmpS2068 = _M0L1jS466 + 1;
              _M0L1jS466 = _M0L6_2atmpS2068;
              continue;
            } else {
              return (int64_t)_M0L1iS463;
            }
            break;
          }
          goto join_464;
          goto joinlet_3417;
          join_464:;
          _M0L6_2atmpS2054 = _M0L1iS463 + 1;
          _M0L1iS463 = _M0L6_2atmpS2054;
          continue;
          joinlet_3417:;
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
  struct _M0TPC16string10StringView _M0L8haystackS445,
  struct _M0TPC16string10StringView _M0L6needleS447
) {
  int32_t _M0L3endS2052;
  int32_t _M0L5startS2053;
  int32_t _M0L13haystack__lenS444;
  int32_t _M0L3endS2050;
  int32_t _M0L5startS2051;
  int32_t _M0L11needle__lenS446;
  #line 57 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2052 = _M0L8haystackS445.$2;
  _M0L5startS2053 = _M0L8haystackS445.$1;
  _M0L13haystack__lenS444 = _M0L3endS2052 - _M0L5startS2053;
  _M0L3endS2050 = _M0L6needleS447.$2;
  _M0L5startS2051 = _M0L6needleS447.$1;
  _M0L11needle__lenS446 = _M0L3endS2050 - _M0L5startS2051;
  if (_M0L11needle__lenS446 > 0) {
    if (_M0L13haystack__lenS444 >= _M0L11needle__lenS446) {
      int32_t* _M0L11skip__tableS448 =
        (int32_t*)moonbit_make_int32_array(256, _M0L11needle__lenS446);
      int32_t _M0L7_2abindS449 = _M0L11needle__lenS446 - 1;
      int32_t _M0L1iS450 = 0;
      int32_t _M0L1iS452;
      while (1) {
        if (_M0L1iS450 < _M0L7_2abindS449) {
          moonbit_string_t _M0L3strS2025 = _M0L6needleS447.$0;
          int32_t _M0L5startS2027 = _M0L6needleS447.$1;
          int32_t _M0L6_2atmpS2026 = _M0L5startS2027 + _M0L1iS450;
          int32_t _M0L6_2atmpS2024 = _M0L3strS2025[_M0L6_2atmpS2026];
          int32_t _M0L6_2atmpS2023 = (int32_t)_M0L6_2atmpS2024;
          int32_t _M0L6_2atmpS2020 = _M0L6_2atmpS2023 & 255;
          int32_t _M0L6_2atmpS2022 = _M0L11needle__lenS446 - 1;
          int32_t _M0L6_2atmpS2021 = _M0L6_2atmpS2022 - _M0L1iS450;
          int32_t _M0L6_2atmpS2028;
          if (
            _M0L6_2atmpS2020 < 0
            || _M0L6_2atmpS2020
               >= Moonbit_array_length(_M0L11skip__tableS448)
          ) {
            #line 68 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L11skip__tableS448[_M0L6_2atmpS2020] = _M0L6_2atmpS2021;
          _M0L6_2atmpS2028 = _M0L1iS450 + 1;
          _M0L1iS450 = _M0L6_2atmpS2028;
          continue;
        }
        break;
      }
      _M0L1iS452 = 0;
      while (1) {
        int32_t _M0L6_2atmpS2029 =
          _M0L13haystack__lenS444 - _M0L11needle__lenS446;
        if (_M0L1iS452 <= _M0L6_2atmpS2029) {
          int32_t _M0L7_2abindS453 = _M0L11needle__lenS446 - 1;
          int32_t _M0L1jS454 = 0;
          moonbit_string_t _M0L3strS2045;
          int32_t _M0L5startS2047;
          int32_t _M0L6_2atmpS2049;
          int32_t _M0L6_2atmpS2048;
          int32_t _M0L6_2atmpS2046;
          int32_t _M0L6_2atmpS2044;
          int32_t _M0L6_2atmpS2043;
          int32_t _M0L6_2atmpS2042;
          int32_t _M0L6_2atmpS2041;
          int32_t _M0L6_2atmpS2040;
          while (1) {
            if (_M0L1jS454 <= _M0L7_2abindS453) {
              moonbit_string_t _M0L3strS2035 = _M0L8haystackS445.$0;
              int32_t _M0L5startS2037 = _M0L8haystackS445.$1;
              int32_t _M0L6_2atmpS2038 = _M0L1iS452 + _M0L1jS454;
              int32_t _M0L6_2atmpS2036 = _M0L5startS2037 + _M0L6_2atmpS2038;
              int32_t _M0L6_2atmpS2030 = _M0L3strS2035[_M0L6_2atmpS2036];
              moonbit_string_t _M0L3strS2032 = _M0L6needleS447.$0;
              int32_t _M0L5startS2034 = _M0L6needleS447.$1;
              int32_t _M0L6_2atmpS2033 = _M0L5startS2034 + _M0L1jS454;
              int32_t _M0L6_2atmpS2031 = _M0L3strS2032[_M0L6_2atmpS2033];
              int32_t _M0L6_2atmpS2039;
              if (_M0L6_2atmpS2030 != _M0L6_2atmpS2031) {
                break;
              }
              _M0L6_2atmpS2039 = _M0L1jS454 + 1;
              _M0L1jS454 = _M0L6_2atmpS2039;
              continue;
            } else {
              moonbit_decref(_M0L11skip__tableS448);
              return (int64_t)_M0L1iS452;
            }
            break;
          }
          _M0L3strS2045 = _M0L8haystackS445.$0;
          _M0L5startS2047 = _M0L8haystackS445.$1;
          _M0L6_2atmpS2049 = _M0L1iS452 + _M0L11needle__lenS446;
          _M0L6_2atmpS2048 = _M0L6_2atmpS2049 - 1;
          _M0L6_2atmpS2046 = _M0L5startS2047 + _M0L6_2atmpS2048;
          _M0L6_2atmpS2044 = _M0L3strS2045[_M0L6_2atmpS2046];
          _M0L6_2atmpS2043 = (int32_t)_M0L6_2atmpS2044;
          _M0L6_2atmpS2042 = _M0L6_2atmpS2043 & 255;
          if (
            _M0L6_2atmpS2042 < 0
            || _M0L6_2atmpS2042
               >= Moonbit_array_length(_M0L11skip__tableS448)
          ) {
            #line 73 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2041 = (int32_t)_M0L11skip__tableS448[_M0L6_2atmpS2042];
          _M0L6_2atmpS2040 = _M0L1iS452 + _M0L6_2atmpS2041;
          _M0L1iS452 = _M0L6_2atmpS2040;
          continue;
        } else {
          moonbit_decref(_M0L11skip__tableS448);
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
  struct _M0TPB13StringBuilder* _M0L4selfS443,
  struct _M0TPC16string10StringView _M0L3strS442
) {
  int32_t _M0L3endS2018;
  int32_t _M0L5startS2019;
  int32_t _M0L8str__lenS441;
  int32_t _M0L3lenS2011;
  int32_t _M0L6_2atmpS2010;
  uint16_t* _M0L4dataS2012;
  int32_t _M0L3lenS2013;
  moonbit_string_t _M0L6_2atmpS2014;
  int32_t _M0L6_2atmpS2015;
  int32_t _M0L3lenS2017;
  int32_t _M0L6_2atmpS2016;
  #line 131 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3endS2018 = _M0L3strS442.$2;
  _M0L5startS2019 = _M0L3strS442.$1;
  _M0L8str__lenS441 = _M0L3endS2018 - _M0L5startS2019;
  _M0L3lenS2011 = _M0L4selfS443->$1;
  _M0L6_2atmpS2010 = _M0L3lenS2011 + _M0L8str__lenS441;
  #line 136 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS443, _M0L6_2atmpS2010);
  _M0L4dataS2012 = _M0L4selfS443->$0;
  _M0L3lenS2013 = _M0L4selfS443->$1;
  moonbit_incref(_M0L4dataS2012);
  #line 139 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS2014 = _M0MPC16string10StringView4data(_M0L3strS442);
  #line 140 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS2015 = _M0MPC16string10StringView13start__offset(_M0L3strS442);
  #line 137 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS2012, _M0L3lenS2013, _M0L6_2atmpS2014, _M0L6_2atmpS2015, _M0L8str__lenS441);
  moonbit_decref(_M0L4dataS2012);
  moonbit_decref(_M0L6_2atmpS2014);
  _M0L3lenS2017 = _M0L4selfS443->$1;
  _M0L6_2atmpS2016 = _M0L3lenS2017 + _M0L8str__lenS441;
  _M0L4selfS443->$1 = _M0L6_2atmpS2016;
  return 0;
}

int64_t _M0MPC16string6String29offset__of__nth__char_2einner(
  moonbit_string_t _M0L4selfS438,
  int32_t _M0L1iS439,
  int32_t _M0L13start__offsetS440,
  int64_t _M0L11end__offsetS436
) {
  int32_t _M0L11end__offsetS435;
  #line 446 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L11end__offsetS436 == 4294967296ll) {
    _M0L11end__offsetS435 = Moonbit_array_length(_M0L4selfS438);
  } else {
    int64_t _M0L7_2aSomeS437 = _M0L11end__offsetS436;
    _M0L11end__offsetS435 = (int32_t)_M0L7_2aSomeS437;
  }
  if (_M0L1iS439 >= 0) {
    #line 455 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0MPC16string6String30offset__of__nth__char__forward(_M0L4selfS438, _M0L1iS439, _M0L13start__offsetS440, _M0L11end__offsetS435);
  } else {
    int32_t _M0L6_2atmpS2009 = -_M0L1iS439;
    #line 458 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0MPC16string6String31offset__of__nth__char__backward(_M0L4selfS438, _M0L6_2atmpS2009, _M0L13start__offsetS440, _M0L11end__offsetS435);
  }
}

int64_t _M0MPC16string6String30offset__of__nth__char__forward(
  moonbit_string_t _M0L4selfS433,
  int32_t _M0L1nS431,
  int32_t _M0L13start__offsetS427,
  int32_t _M0L11end__offsetS428
) {
  int32_t _if__result_3422;
  #line 375 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L13start__offsetS427 >= 0) {
    _if__result_3422 = _M0L13start__offsetS427 <= _M0L11end__offsetS428;
  } else {
    _if__result_3422 = 0;
  }
  if (_if__result_3422) {
    int32_t _M0L13utf16__offsetS429 = _M0L13start__offsetS427;
    int32_t _M0L11char__countS430 = 0;
    while (1) {
      int32_t _if__result_3424;
      if (_M0L13utf16__offsetS429 < _M0L11end__offsetS428) {
        _if__result_3424 = _M0L11char__countS430 < _M0L1nS431;
      } else {
        _if__result_3424 = 0;
      }
      if (_if__result_3424) {
        int32_t _M0L1cS432 = _M0L4selfS433[_M0L13utf16__offsetS429];
        #line 389 "/home/developer/.moon/lib/core/builtin/string.mbt"
        if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L1cS432)) {
          int32_t _M0L6_2atmpS2005 = _M0L13utf16__offsetS429 + 2;
          int32_t _M0L6_2atmpS2006 = _M0L11char__countS430 + 1;
          _M0L13utf16__offsetS429 = _M0L6_2atmpS2005;
          _M0L11char__countS430 = _M0L6_2atmpS2006;
          continue;
        } else {
          int32_t _M0L6_2atmpS2007 = _M0L13utf16__offsetS429 + 1;
          int32_t _M0L6_2atmpS2008 = _M0L11char__countS430 + 1;
          _M0L13utf16__offsetS429 = _M0L6_2atmpS2007;
          _M0L11char__countS430 = _M0L6_2atmpS2008;
          continue;
        }
      } else {
        int32_t _if__result_3425;
        if (_M0L11char__countS430 < _M0L1nS431) {
          _if__result_3425 = 1;
        } else {
          _if__result_3425 = _M0L13utf16__offsetS429 >= _M0L11end__offsetS428;
        }
        if (_if__result_3425) {
          return 4294967296ll;
        } else {
          return (int64_t)_M0L13utf16__offsetS429;
        }
      }
      break;
    }
  } else {
    #line 382 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0FPC15abort5abortGOiE((moonbit_string_t)moonbit_string_literal_109.data);
  }
}

int64_t _M0MPC16string6String31offset__of__nth__char__backward(
  moonbit_string_t _M0L4selfS424,
  int32_t _M0L1nS422,
  int32_t _M0L13start__offsetS421,
  int32_t _M0L11end__offsetS426
) {
  int32_t _M0L13utf16__offsetS419;
  int32_t _M0L11char__countS420;
  #line 411 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L13utf16__offsetS419 = _M0L11end__offsetS426;
  _M0L11char__countS420 = 0;
  while (1) {
    int32_t _M0L6_2atmpS1999 = _M0L13utf16__offsetS419 - 1;
    int32_t _if__result_3427;
    if (_M0L6_2atmpS1999 >= _M0L13start__offsetS421) {
      _if__result_3427 = _M0L11char__countS420 < _M0L1nS422;
    } else {
      _if__result_3427 = 0;
    }
    if (_if__result_3427) {
      int32_t _M0L6_2atmpS2004 = _M0L13utf16__offsetS419 - 1;
      int32_t _M0L1cS423 = _M0L4selfS424[_M0L6_2atmpS2004];
      #line 424 "/home/developer/.moon/lib/core/builtin/string.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L1cS423)) {
        int32_t _M0L6_2atmpS2000 = _M0L13utf16__offsetS419 - 2;
        int32_t _M0L6_2atmpS2001 = _M0L11char__countS420 + 1;
        _M0L13utf16__offsetS419 = _M0L6_2atmpS2000;
        _M0L11char__countS420 = _M0L6_2atmpS2001;
        continue;
      } else {
        int32_t _M0L6_2atmpS2002 = _M0L13utf16__offsetS419 - 1;
        int32_t _M0L6_2atmpS2003 = _M0L11char__countS420 + 1;
        _M0L13utf16__offsetS419 = _M0L6_2atmpS2002;
        _M0L11char__countS420 = _M0L6_2atmpS2003;
        continue;
      }
    } else {
      int32_t _if__result_3428;
      if (_M0L11char__countS420 < _M0L1nS422) {
        _if__result_3428 = 1;
      } else {
        _if__result_3428 = _M0L13utf16__offsetS419 < _M0L13start__offsetS421;
      }
      if (_if__result_3428) {
        return 4294967296ll;
      } else {
        return (int64_t)_M0L13utf16__offsetS419;
      }
    }
    break;
  }
}

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t _M0L4selfS417,
  int32_t _M0L13start__offsetS418,
  int64_t _M0L11end__offsetS415
) {
  int32_t _M0L11end__offsetS414;
  int32_t _if__result_3429;
  #line 614 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS415 == 4294967296ll) {
    _M0L11end__offsetS414 = Moonbit_array_length(_M0L4selfS417);
  } else {
    int64_t _M0L7_2aSomeS416 = _M0L11end__offsetS415;
    _M0L11end__offsetS414 = (int32_t)_M0L7_2aSomeS416;
  }
  if (_M0L13start__offsetS418 >= 0) {
    if (_M0L13start__offsetS418 <= _M0L11end__offsetS414) {
      int32_t _M0L6_2atmpS1998 = Moonbit_array_length(_M0L4selfS417);
      _if__result_3429 = _M0L11end__offsetS414 <= _M0L6_2atmpS1998;
    } else {
      _if__result_3429 = 0;
    }
  } else {
    _if__result_3429 = 0;
  }
  if (_if__result_3429) {
    moonbit_incref(_M0L4selfS417);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS417,
                                                 .$1 = _M0L13start__offsetS418,
                                                 .$2 = _M0L11end__offsetS414};
  } else {
    #line 623 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_110.data);
  }
}

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView _M0L4selfS409
) {
  int32_t _M0L5startS408;
  int32_t _M0L3endS410;
  struct _M0TPB8MutLocalGiE* _M0L5indexS411;
  struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__* _closure_3430;
  struct _M0TWEOc* _M0L6_2atmpS1977;
  struct _M0TPB4IterGcE* _result_3431;
  #line 214 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L5startS408 = _M0L4selfS409.$1;
  _M0L3endS410 = _M0L4selfS409.$2;
  _M0L5indexS411
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS411)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS411->$0 = _M0L5startS408;
  moonbit_incref(_M0L4selfS409.$0);
  _closure_3430
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__*)moonbit_malloc(sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__));
  Moonbit_object_header(_closure_3430)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 65, 0);
  _closure_3430->code = &_M0MPC16string10StringView4iterC1978l219;
  _closure_3430->$0 = _M0L5indexS411;
  _closure_3430->$1 = _M0L3endS410;
  _closure_3430->$2 = _M0L4selfS409;
  _M0L6_2atmpS1977 = (struct _M0TWEOc*)_closure_3430;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_3431 = _M0MPB4Iter3newGcE(_M0L6_2atmpS1977, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1977);
  return _result_3431;
}

int32_t _M0MPC16string10StringView4iterC1978l219(
  struct _M0TWEOc* _M0L6_2aenvS1979
) {
  struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__* _M0L14_2acasted__envS1980;
  struct _M0TPC16string10StringView _M0L4selfS409;
  int32_t _M0L3endS410;
  struct _M0TPB8MutLocalGiE* _M0L5indexS411;
  int32_t _M0L3valS1981;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L14_2acasted__envS1980
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1978__l219__*)_M0L6_2aenvS1979;
  _M0L4selfS409 = _M0L14_2acasted__envS1980->$2;
  _M0L3endS410 = _M0L14_2acasted__envS1980->$1;
  _M0L5indexS411 = _M0L14_2acasted__envS1980->$0;
  _M0L3valS1981 = _M0L5indexS411->$0;
  if (_M0L3valS1981 < _M0L3endS410) {
    moonbit_string_t _M0L3strS1996 = _M0L4selfS409.$0;
    int32_t _M0L3valS1997 = _M0L5indexS411->$0;
    int32_t _M0L2c1S412 = _M0L3strS1996[_M0L3valS1997];
    int32_t _if__result_3432;
    int32_t _M0L3valS1994;
    int32_t _M0L6_2atmpS1993;
    int32_t _M0L6_2atmpS1995;
    #line 222 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S412)) {
      int32_t _M0L3valS1984 = _M0L5indexS411->$0;
      int32_t _M0L6_2atmpS1982 = _M0L3valS1984 + 1;
      int32_t _M0L3endS1983 = _M0L4selfS409.$2;
      _if__result_3432 = _M0L6_2atmpS1982 < _M0L3endS1983;
    } else {
      _if__result_3432 = 0;
    }
    if (_if__result_3432) {
      moonbit_string_t _M0L3strS1990 = _M0L4selfS409.$0;
      int32_t _M0L3valS1992 = _M0L5indexS411->$0;
      int32_t _M0L6_2atmpS1991 = _M0L3valS1992 + 1;
      int32_t _M0L2c2S413 = _M0L3strS1990[_M0L6_2atmpS1991];
      #line 224 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S413)) {
        int32_t _M0L3valS1986 = _M0L5indexS411->$0;
        int32_t _M0L6_2atmpS1985 = _M0L3valS1986 + 2;
        int32_t _M0L6_2atmpS1988;
        int32_t _M0L6_2atmpS1989;
        int32_t _M0L6_2atmpS1987;
        _M0L5indexS411->$0 = _M0L6_2atmpS1985;
        _M0L6_2atmpS1988 = (int32_t)_M0L2c1S412;
        _M0L6_2atmpS1989 = (int32_t)_M0L2c2S413;
        #line 226 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        _M0L6_2atmpS1987
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1988, _M0L6_2atmpS1989);
        return _M0L6_2atmpS1987;
      }
    }
    _M0L3valS1994 = _M0L5indexS411->$0;
    _M0L6_2atmpS1993 = _M0L3valS1994 + 1;
    _M0L5indexS411->$0 = _M0L6_2atmpS1993;
    #line 230 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    _M0L6_2atmpS1995 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S412);
    return _M0L6_2atmpS1995;
  } else {
    return -1;
  }
}

moonbit_string_t _M0MPC16string10StringView9to__owned(
  struct _M0TPC16string10StringView _M0L4selfS407
) {
  moonbit_string_t _M0L3strS1974;
  int32_t _M0L5startS1975;
  int32_t _M0L3endS1976;
  moonbit_string_t _result_3433;
  #line 196 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1974 = _M0L4selfS407.$0;
  _M0L5startS1975 = _M0L4selfS407.$1;
  _M0L3endS1976 = _M0L4selfS407.$2;
  moonbit_incref(_M0L3strS1974);
  #line 199 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_3433
  = _M0MPC16string6String17unsafe__substring(_M0L3strS1974, _M0L5startS1975, _M0L3endS1976);
  moonbit_decref(_M0L3strS1974);
  return _result_3433;
}

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t _M0L3strS404,
  int32_t _M0L5startS402,
  int32_t _M0L3endS403
) {
  int32_t _if__result_3434;
  int32_t _M0L3lenS405;
  int32_t _M0L6_2atmpS1972;
  int32_t _M0L6_2atmpS1973;
  moonbit_bytes_t _M0L5bytesS406;
  moonbit_bytes_t _M0L6_2atmpS1971;
  moonbit_string_t _result_3435;
  #line 91 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L5startS402 == 0) {
    int32_t _M0L6_2atmpS1970 = Moonbit_array_length(_M0L3strS404);
    _if__result_3434 = _M0L3endS403 == _M0L6_2atmpS1970;
  } else {
    _if__result_3434 = 0;
  }
  if (_if__result_3434) {
    moonbit_incref(_M0L3strS404);
    return _M0L3strS404;
  }
  _M0L3lenS405 = _M0L3endS403 - _M0L5startS402;
  _M0L6_2atmpS1972 = _M0L3lenS405 * 2;
  #line 101 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1973 = _M0IPC14byte4BytePB7Default7default();
  _M0L5bytesS406
  = (moonbit_bytes_t)moonbit_make_bytes(_M0L6_2atmpS1972, _M0L6_2atmpS1973);
  #line 102 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0MPC15array10FixedArray18blit__from__string(_M0L5bytesS406, 0, _M0L3strS404, _M0L5startS402, _M0L3lenS405);
  _M0L6_2atmpS1971 = _M0L5bytesS406;
  #line 103 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_3435
  = _M0MPC15bytes5Bytes29to__unchecked__string_2einner(_M0L6_2atmpS1971, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1971);
  return _result_3435;
}

int32_t _M0IPC14byte4BytePB7Default7default() {
  #line 231 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  return 0;
}

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t _M0L4selfS397,
  int32_t _M0L6offsetS401,
  int64_t _M0L6lengthS399
) {
  int32_t _M0L3lenS396;
  int32_t _M0L6lengthS398;
  int32_t _if__result_3436;
  #line 77 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L3lenS396 = Moonbit_array_length(_M0L4selfS397);
  if (_M0L6lengthS399 == 4294967296ll) {
    _M0L6lengthS398 = _M0L3lenS396 - _M0L6offsetS401;
  } else {
    int64_t _M0L7_2aSomeS400 = _M0L6lengthS399;
    _M0L6lengthS398 = (int32_t)_M0L7_2aSomeS400;
  }
  if (_M0L6offsetS401 >= 0) {
    if (_M0L6lengthS398 >= 0) {
      int32_t _M0L6_2atmpS1969 = _M0L6offsetS401 + _M0L6lengthS398;
      _if__result_3436 = _M0L6_2atmpS1969 <= _M0L3lenS396;
    } else {
      _if__result_3436 = 0;
    }
  } else {
    _if__result_3436 = 0;
  }
  if (_if__result_3436) {
    moonbit_incref(_M0L4selfS397);
    #line 85 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    return _M0FPB19unsafe__sub__string(_M0L4selfS397, _M0L6offsetS401, _M0L6lengthS398);
  } else {
    #line 84 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t _M0L4selfS388,
  int32_t _M0L13bytes__offsetS383,
  moonbit_string_t _M0L3strS390,
  int32_t _M0L11str__offsetS386,
  int32_t _M0L6lengthS384
) {
  int32_t _M0L6_2atmpS1968;
  int32_t _M0L6_2atmpS1967;
  int32_t _M0L2e1S382;
  int32_t _M0L6_2atmpS1966;
  int32_t _M0L2e2S385;
  int32_t _M0L4len1S387;
  int32_t _M0L4len2S389;
  int32_t _if__result_3437;
  #line 125 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L6_2atmpS1968 = _M0L6lengthS384 * 2;
  _M0L6_2atmpS1967 = _M0L13bytes__offsetS383 + _M0L6_2atmpS1968;
  _M0L2e1S382 = _M0L6_2atmpS1967 - 1;
  _M0L6_2atmpS1966 = _M0L11str__offsetS386 + _M0L6lengthS384;
  _M0L2e2S385 = _M0L6_2atmpS1966 - 1;
  _M0L4len1S387 = Moonbit_array_length(_M0L4selfS388);
  _M0L4len2S389 = Moonbit_array_length(_M0L3strS390);
  if (_M0L6lengthS384 >= 0) {
    if (_M0L13bytes__offsetS383 >= 0) {
      if (_M0L2e1S382 < _M0L4len1S387) {
        if (_M0L11str__offsetS386 >= 0) {
          _if__result_3437 = _M0L2e2S385 < _M0L4len2S389;
        } else {
          _if__result_3437 = 0;
        }
      } else {
        _if__result_3437 = 0;
      }
    } else {
      _if__result_3437 = 0;
    }
  } else {
    _if__result_3437 = 0;
  }
  if (_if__result_3437) {
    int32_t _M0L16end__str__offsetS391 =
      _M0L11str__offsetS386 + _M0L6lengthS384;
    int32_t _M0L1iS392 = _M0L11str__offsetS386;
    int32_t _M0L1jS393 = _M0L13bytes__offsetS383;
    while (1) {
      if (_M0L1iS392 < _M0L16end__str__offsetS391) {
        int32_t _M0L6_2atmpS1963 = _M0L3strS390[_M0L1iS392];
        int32_t _M0L6_2atmpS1962 = (int32_t)_M0L6_2atmpS1963;
        uint32_t _M0L1cS394 = *(uint32_t*)&_M0L6_2atmpS1962;
        uint32_t _M0L6_2atmpS1958 = _M0L1cS394 & 255u;
        int32_t _M0L6_2atmpS1957;
        int32_t _M0L6_2atmpS1959;
        uint32_t _M0L6_2atmpS1961;
        int32_t _M0L6_2atmpS1960;
        int32_t _M0L6_2atmpS1964;
        int32_t _M0L6_2atmpS1965;
        #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1957 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1958);
        if (
          _M0L1jS393 < 0 || _M0L1jS393 >= Moonbit_array_length(_M0L4selfS388)
        ) {
          #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS388[_M0L1jS393] = _M0L6_2atmpS1957;
        _M0L6_2atmpS1959 = _M0L1jS393 + 1;
        _M0L6_2atmpS1961 = _M0L1cS394 >> 8;
        #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1960 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1961);
        if (
          _M0L6_2atmpS1959 < 0
          || _M0L6_2atmpS1959 >= Moonbit_array_length(_M0L4selfS388)
        ) {
          #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS388[_M0L6_2atmpS1959] = _M0L6_2atmpS1960;
        _M0L6_2atmpS1964 = _M0L1iS392 + 1;
        _M0L6_2atmpS1965 = _M0L1jS393 + 2;
        _M0L1iS392 = _M0L6_2atmpS1964;
        _M0L1jS393 = _M0L6_2atmpS1965;
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

int32_t _M0MPC14uint4UInt8to__byte(uint32_t _M0L4selfS381) {
  int32_t _M0L6_2atmpS1956;
  #line 2519 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1956 = *(int32_t*)&_M0L4selfS381;
  return _M0L6_2atmpS1956 & 0xff;
}

struct _M0TP26biolab8bio__seq3Seq** _M0MPC15array5Array6bufferGRP26biolab8bio__seq3SeqE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq3SeqE* _M0L4selfS377
) {
  struct _M0TP26biolab8bio__seq3Seq** _M0L8_2afieldS3178;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3178 = _M0L4selfS377->$0;
  moonbit_incref(_M0L8_2afieldS3178);
  return _M0L8_2afieldS3178;
}

moonbit_string_t* _M0MPC15array5Array6bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS378
) {
  moonbit_string_t* _M0L8_2afieldS3179;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3179 = _M0L4selfS378->$0;
  moonbit_incref(_M0L8_2afieldS3179);
  return _M0L8_2afieldS3179;
}

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS379
) {
  struct _M0TUsiE** _M0L8_2afieldS3180;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3180 = _M0L4selfS379->$0;
  moonbit_incref(_M0L8_2afieldS3180);
  return _M0L8_2afieldS3180;
}

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0MPC15array5Array6bufferGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TPB5ArrayGRP46biolab8bio__seq3cmd4main4CaseE* _M0L4selfS380
) {
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L8_2afieldS3181;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3181 = _M0L4selfS380->$0;
  moonbit_incref(_M0L8_2afieldS3181);
  return _M0L8_2afieldS3181;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView _M0L4selfS375,
  int32_t _M0L13start__offsetS376,
  int64_t _M0L11end__offsetS373
) {
  int32_t _M0L11end__offsetS372;
  int32_t _if__result_3439;
  #line 105 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS373 == 4294967296ll) {
    int32_t _M0L3endS1954 = _M0L4selfS375.$2;
    int32_t _M0L5startS1955 = _M0L4selfS375.$1;
    _M0L11end__offsetS372 = _M0L3endS1954 - _M0L5startS1955;
  } else {
    int64_t _M0L7_2aSomeS374 = _M0L11end__offsetS373;
    _M0L11end__offsetS372 = (int32_t)_M0L7_2aSomeS374;
  }
  if (_M0L13start__offsetS376 >= 0) {
    if (_M0L13start__offsetS376 <= _M0L11end__offsetS372) {
      int32_t _M0L3endS1947 = _M0L4selfS375.$2;
      int32_t _M0L5startS1948 = _M0L4selfS375.$1;
      int32_t _M0L6_2atmpS1946 = _M0L3endS1947 - _M0L5startS1948;
      _if__result_3439 = _M0L11end__offsetS372 <= _M0L6_2atmpS1946;
    } else {
      _if__result_3439 = 0;
    }
  } else {
    _if__result_3439 = 0;
  }
  if (_if__result_3439) {
    moonbit_string_t _M0L3strS1949 = _M0L4selfS375.$0;
    int32_t _M0L5startS1953 = _M0L4selfS375.$1;
    int32_t _M0L6_2atmpS1950 = _M0L5startS1953 + _M0L13start__offsetS376;
    int32_t _M0L5startS1952 = _M0L4selfS375.$1;
    int32_t _M0L6_2atmpS1951 = _M0L5startS1952 + _M0L11end__offsetS372;
    moonbit_incref(_M0L3strS1949);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1949,
                                                 .$1 = _M0L6_2atmpS1950,
                                                 .$2 = _M0L6_2atmpS1951};
  } else {
    #line 114 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_110.data);
  }
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L1fS366,
  int64_t _M0L10size__hintS363
) {
  int64_t _M0L10size__hintS362;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_3440;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS363 == 4294967296ll) {
    _M0L10size__hintS362 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS364 = _M0L10size__hintS363;
    int32_t _M0L4_2anS365 = (int32_t)_M0L7_2aSomeS364;
    if (_M0L4_2anS365 > 0) {
      _M0L10size__hintS362 = (int64_t)_M0L4_2anS365;
    } else {
      _M0L10size__hintS362
      = _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE;
    }
  }
  moonbit_incref(_M0L1fS366);
  _block_3440
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_3440)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 55, 0);
  _block_3440->$0 = _M0L1fS366;
  _block_3440->$1 = _M0L10size__hintS362;
  return _block_3440;
}

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(
  struct _M0TWEOc* _M0L1fS371,
  int64_t _M0L10size__hintS368
) {
  int64_t _M0L10size__hintS367;
  struct _M0TPB4IterGcE* _block_3441;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS368 == 4294967296ll) {
    _M0L10size__hintS367 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS369 = _M0L10size__hintS368;
    int32_t _M0L4_2anS370 = (int32_t)_M0L7_2aSomeS369;
    if (_M0L4_2anS370 > 0) {
      _M0L10size__hintS367 = (int64_t)_M0L4_2anS370;
    } else {
      _M0L10size__hintS367 = _M0MPB4Iter3newN6constrS9988GcE;
    }
  }
  moonbit_incref(_M0L1fS371);
  _block_3441
  = (struct _M0TPB4IterGcE*)moonbit_malloc(sizeof(struct _M0TPB4IterGcE));
  Moonbit_object_header(_block_3441)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 69, 0);
  _block_3441->$0 = _M0L1fS371;
  _block_3441->$1 = _M0L10size__hintS367;
  return _block_3441;
}

moonbit_string_t _M0MPC13int3Int18to__string_2einner(
  int32_t _M0L4selfS346,
  int32_t _M0L5radixS345
) {
  int32_t _if__result_3442;
  int32_t _M0L12is__negativeS347;
  uint32_t _M0L3numS348;
  uint16_t* _M0L6bufferS349;
  #line 209 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS345 < 2) {
    _if__result_3442 = 1;
  } else {
    _if__result_3442 = _M0L5radixS345 > 36;
  }
  if (_if__result_3442) {
    #line 213 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_111.data);
  }
  if (_M0L4selfS346 == 0) {
    return (moonbit_string_t)moonbit_string_literal_112.data;
  }
  _M0L12is__negativeS347 = _M0L4selfS346 < 0;
  if (_M0L12is__negativeS347) {
    int32_t _M0L6_2atmpS1945 = -_M0L4selfS346;
    _M0L3numS348 = *(uint32_t*)&_M0L6_2atmpS1945;
  } else {
    _M0L3numS348 = *(uint32_t*)&_M0L4selfS346;
  }
  switch (_M0L5radixS345) {
    case 10: {
      int32_t _M0L10digit__lenS350;
      int32_t _M0L6_2atmpS1942;
      int32_t _M0L10total__lenS351;
      uint16_t* _M0L6bufferS352;
      int32_t _M0L12digit__startS353;
      #line 235 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS350 = _M0FPB12dec__count32(_M0L3numS348);
      if (_M0L12is__negativeS347) {
        _M0L6_2atmpS1942 = 1;
      } else {
        _M0L6_2atmpS1942 = 0;
      }
      _M0L10total__lenS351 = _M0L10digit__lenS350 + _M0L6_2atmpS1942;
      _M0L6bufferS352
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS351, 0);
      if (_M0L12is__negativeS347) {
        _M0L12digit__startS353 = 1;
      } else {
        _M0L12digit__startS353 = 0;
      }
      #line 239 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__dec(_M0L6bufferS352, _M0L3numS348, _M0L12digit__startS353, _M0L10total__lenS351);
      _M0L6bufferS349 = _M0L6bufferS352;
      break;
    }
    
    case 16: {
      int32_t _M0L10digit__lenS354;
      int32_t _M0L6_2atmpS1943;
      int32_t _M0L10total__lenS355;
      uint16_t* _M0L6bufferS356;
      int32_t _M0L12digit__startS357;
      #line 243 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS354 = _M0FPB12hex__count32(_M0L3numS348);
      if (_M0L12is__negativeS347) {
        _M0L6_2atmpS1943 = 1;
      } else {
        _M0L6_2atmpS1943 = 0;
      }
      _M0L10total__lenS355 = _M0L10digit__lenS354 + _M0L6_2atmpS1943;
      _M0L6bufferS356
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS355, 0);
      if (_M0L12is__negativeS347) {
        _M0L12digit__startS357 = 1;
      } else {
        _M0L12digit__startS357 = 0;
      }
      #line 247 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__hex(_M0L6bufferS356, _M0L3numS348, _M0L12digit__startS357, _M0L10total__lenS355);
      _M0L6bufferS349 = _M0L6bufferS356;
      break;
    }
    default: {
      int32_t _M0L10digit__lenS358;
      int32_t _M0L6_2atmpS1944;
      int32_t _M0L10total__lenS359;
      uint16_t* _M0L6bufferS360;
      int32_t _M0L12digit__startS361;
      #line 251 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS358
      = _M0FPB14radix__count32(_M0L3numS348, _M0L5radixS345);
      if (_M0L12is__negativeS347) {
        _M0L6_2atmpS1944 = 1;
      } else {
        _M0L6_2atmpS1944 = 0;
      }
      _M0L10total__lenS359 = _M0L10digit__lenS358 + _M0L6_2atmpS1944;
      _M0L6bufferS360
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS359, 0);
      if (_M0L12is__negativeS347) {
        _M0L12digit__startS361 = 1;
      } else {
        _M0L12digit__startS361 = 0;
      }
      #line 255 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB24int__to__string__generic(_M0L6bufferS360, _M0L3numS348, _M0L12digit__startS361, _M0L10total__lenS359, _M0L5radixS345);
      _M0L6bufferS349 = _M0L6bufferS360;
      break;
    }
  }
  if (_M0L12is__negativeS347) {
    _M0L6bufferS349[0] = 45;
  }
  return _M0L6bufferS349;
}

int32_t _M0FPB14radix__count32(
  uint32_t _M0L5valueS339,
  int32_t _M0L5radixS341
) {
  uint32_t _M0L4baseS340;
  uint32_t _M0L3numS342;
  int32_t _M0L5countS343;
  #line 189 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS339 == 0u) {
    return 1;
  }
  _M0L4baseS340 = *(uint32_t*)&_M0L5radixS341;
  _M0L3numS342 = _M0L5valueS339;
  _M0L5countS343 = 0;
  while (1) {
    if (_M0L3numS342 > 0u) {
      uint32_t _M0L6_2atmpS1940 = _M0L3numS342 / _M0L4baseS340;
      int32_t _M0L6_2atmpS1941 = _M0L5countS343 + 1;
      _M0L3numS342 = _M0L6_2atmpS1940;
      _M0L5countS343 = _M0L6_2atmpS1941;
      continue;
    } else {
      return _M0L5countS343;
    }
    break;
  }
}

int32_t _M0FPB12hex__count32(uint32_t _M0L5valueS337) {
  #line 177 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS337 == 0u) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS338;
    int32_t _M0L6_2atmpS1939;
    int32_t _M0L6_2atmpS1938;
    #line 182 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS338 = moonbit_clz32(_M0L5valueS337);
    _M0L6_2atmpS1939 = 31 - _M0L14leading__zerosS338;
    _M0L6_2atmpS1938 = _M0L6_2atmpS1939 / 4;
    return _M0L6_2atmpS1938 + 1;
  }
}

int32_t _M0FPB12dec__count32(uint32_t _M0L5valueS336) {
  #line 143 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS336 >= 100000u) {
    if (_M0L5valueS336 >= 10000000u) {
      if (_M0L5valueS336 >= 1000000000u) {
        return 10;
      } else if (_M0L5valueS336 >= 100000000u) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS336 >= 1000000u) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS336 >= 1000u) {
    if (_M0L5valueS336 >= 10000u) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS336 >= 100u) {
    return 3;
  } else if (_M0L5valueS336 >= 10u) {
    return 2;
  } else {
    return 1;
  }
}

int32_t _M0FPB20int__to__string__dec(
  uint16_t* _M0L6bufferS322,
  uint32_t _M0L3numS334,
  int32_t _M0L12digit__startS323,
  int32_t _M0L10total__lenS335
) {
  int32_t _M0L6_2atmpS1937;
  uint32_t _M0L3numS312;
  int32_t _M0L6offsetS313;
  #line 88 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1937 = _M0L10total__lenS335 - _M0L12digit__startS323;
  _M0L3numS312 = _M0L3numS334;
  _M0L6offsetS313 = _M0L6_2atmpS1937;
  while (1) {
    if (_M0L3numS312 >= 10000u) {
      uint32_t _M0L1tS314 = _M0L3numS312 / 10000u;
      uint32_t _M0L6_2atmpS1914 = _M0L3numS312 % 10000u;
      int32_t _M0L1rS315 = *(int32_t*)&_M0L6_2atmpS1914;
      int32_t _M0L2d1S316 = _M0L1rS315 / 100;
      int32_t _M0L2d2S317 = _M0L1rS315 % 100;
      int32_t _M0L6_2atmpS1913 = _M0L2d1S316 / 10;
      int32_t _M0L6_2atmpS1912 = 48 + _M0L6_2atmpS1913;
      int32_t _M0L6d1__hiS318 = (uint16_t)_M0L6_2atmpS1912;
      int32_t _M0L6_2atmpS1911 = _M0L2d1S316 % 10;
      int32_t _M0L6_2atmpS1910 = 48 + _M0L6_2atmpS1911;
      int32_t _M0L6d1__loS319 = (uint16_t)_M0L6_2atmpS1910;
      int32_t _M0L6_2atmpS1909 = _M0L2d2S317 / 10;
      int32_t _M0L6_2atmpS1908 = 48 + _M0L6_2atmpS1909;
      int32_t _M0L6d2__hiS320 = (uint16_t)_M0L6_2atmpS1908;
      int32_t _M0L6_2atmpS1907 = _M0L2d2S317 % 10;
      int32_t _M0L6_2atmpS1906 = 48 + _M0L6_2atmpS1907;
      int32_t _M0L6d2__loS321 = (uint16_t)_M0L6_2atmpS1906;
      int32_t _M0L6_2atmpS1898 = _M0L12digit__startS323 + _M0L6offsetS313;
      int32_t _M0L6_2atmpS1897 = _M0L6_2atmpS1898 - 4;
      int32_t _M0L6_2atmpS1900;
      int32_t _M0L6_2atmpS1899;
      int32_t _M0L6_2atmpS1902;
      int32_t _M0L6_2atmpS1901;
      int32_t _M0L6_2atmpS1904;
      int32_t _M0L6_2atmpS1903;
      int32_t _M0L6_2atmpS1905;
      _M0L6bufferS322[_M0L6_2atmpS1897] = _M0L6d1__hiS318;
      _M0L6_2atmpS1900 = _M0L12digit__startS323 + _M0L6offsetS313;
      _M0L6_2atmpS1899 = _M0L6_2atmpS1900 - 3;
      _M0L6bufferS322[_M0L6_2atmpS1899] = _M0L6d1__loS319;
      _M0L6_2atmpS1902 = _M0L12digit__startS323 + _M0L6offsetS313;
      _M0L6_2atmpS1901 = _M0L6_2atmpS1902 - 2;
      _M0L6bufferS322[_M0L6_2atmpS1901] = _M0L6d2__hiS320;
      _M0L6_2atmpS1904 = _M0L12digit__startS323 + _M0L6offsetS313;
      _M0L6_2atmpS1903 = _M0L6_2atmpS1904 - 1;
      _M0L6bufferS322[_M0L6_2atmpS1903] = _M0L6d2__loS321;
      _M0L6_2atmpS1905 = _M0L6offsetS313 - 4;
      _M0L3numS312 = _M0L1tS314;
      _M0L6offsetS313 = _M0L6_2atmpS1905;
      continue;
    } else {
      int32_t _M0L6_2atmpS1936 = *(int32_t*)&_M0L3numS312;
      int32_t _M0L9remainingS325 = _M0L6_2atmpS1936;
      int32_t _M0L6offsetS326 = _M0L6offsetS313;
      while (1) {
        if (_M0L9remainingS325 >= 100) {
          int32_t _M0L1tS327 = _M0L9remainingS325 / 100;
          int32_t _M0L1dS328 = _M0L9remainingS325 % 100;
          int32_t _M0L6_2atmpS1923 = _M0L1dS328 / 10;
          int32_t _M0L6_2atmpS1922 = 48 + _M0L6_2atmpS1923;
          int32_t _M0L5d__hiS329 = (uint16_t)_M0L6_2atmpS1922;
          int32_t _M0L6_2atmpS1921 = _M0L1dS328 % 10;
          int32_t _M0L6_2atmpS1920 = 48 + _M0L6_2atmpS1921;
          int32_t _M0L5d__loS330 = (uint16_t)_M0L6_2atmpS1920;
          int32_t _M0L6_2atmpS1916 = _M0L12digit__startS323 + _M0L6offsetS326;
          int32_t _M0L6_2atmpS1915 = _M0L6_2atmpS1916 - 2;
          int32_t _M0L6_2atmpS1918;
          int32_t _M0L6_2atmpS1917;
          int32_t _M0L6_2atmpS1919;
          _M0L6bufferS322[_M0L6_2atmpS1915] = _M0L5d__hiS329;
          _M0L6_2atmpS1918 = _M0L12digit__startS323 + _M0L6offsetS326;
          _M0L6_2atmpS1917 = _M0L6_2atmpS1918 - 1;
          _M0L6bufferS322[_M0L6_2atmpS1917] = _M0L5d__loS330;
          _M0L6_2atmpS1919 = _M0L6offsetS326 - 2;
          _M0L9remainingS325 = _M0L1tS327;
          _M0L6offsetS326 = _M0L6_2atmpS1919;
          continue;
        } else if (_M0L9remainingS325 >= 10) {
          int32_t _M0L6_2atmpS1931 = _M0L9remainingS325 / 10;
          int32_t _M0L6_2atmpS1930 = 48 + _M0L6_2atmpS1931;
          int32_t _M0L5d__hiS332 = (uint16_t)_M0L6_2atmpS1930;
          int32_t _M0L6_2atmpS1929 = _M0L9remainingS325 % 10;
          int32_t _M0L6_2atmpS1928 = 48 + _M0L6_2atmpS1929;
          int32_t _M0L5d__loS333 = (uint16_t)_M0L6_2atmpS1928;
          int32_t _M0L6_2atmpS1925 = _M0L12digit__startS323 + _M0L6offsetS326;
          int32_t _M0L6_2atmpS1924 = _M0L6_2atmpS1925 - 2;
          int32_t _M0L6_2atmpS1927;
          int32_t _M0L6_2atmpS1926;
          _M0L6bufferS322[_M0L6_2atmpS1924] = _M0L5d__hiS332;
          _M0L6_2atmpS1927 = _M0L12digit__startS323 + _M0L6offsetS326;
          _M0L6_2atmpS1926 = _M0L6_2atmpS1927 - 1;
          _M0L6bufferS322[_M0L6_2atmpS1926] = _M0L5d__loS333;
        } else {
          int32_t _M0L6_2atmpS1935 = _M0L12digit__startS323 + _M0L6offsetS326;
          int32_t _M0L6_2atmpS1932 = _M0L6_2atmpS1935 - 1;
          int32_t _M0L6_2atmpS1934 = 48 + _M0L9remainingS325;
          int32_t _M0L6_2atmpS1933 = (uint16_t)_M0L6_2atmpS1934;
          _M0L6bufferS322[_M0L6_2atmpS1932] = _M0L6_2atmpS1933;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB24int__to__string__generic(
  uint16_t* _M0L6bufferS302,
  uint32_t _M0L3numS306,
  int32_t _M0L12digit__startS303,
  int32_t _M0L10total__lenS305,
  int32_t _M0L5radixS296
) {
  uint32_t _M0L4baseS295;
  int32_t _M0L6_2atmpS1882;
  int32_t _M0L6_2atmpS1881;
  #line 57 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS295 = *(uint32_t*)&_M0L5radixS296;
  _M0L6_2atmpS1882 = _M0L5radixS296 - 1;
  _M0L6_2atmpS1881 = _M0L5radixS296 & _M0L6_2atmpS1882;
  if (_M0L6_2atmpS1881 == 0) {
    int32_t _M0L5shiftS297;
    uint32_t _M0L4maskS298;
    int32_t _M0L6_2atmpS1889;
    int32_t _M0L6offsetS299;
    uint32_t _M0L1nS300;
    #line 68 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS297 = moonbit_ctz32(_M0L5radixS296);
    _M0L4maskS298 = _M0L4baseS295 - 1u;
    _M0L6_2atmpS1889 = _M0L10total__lenS305 - _M0L12digit__startS303;
    _M0L6offsetS299 = _M0L6_2atmpS1889;
    _M0L1nS300 = _M0L3numS306;
    while (1) {
      if (_M0L1nS300 > 0u) {
        uint32_t _M0L6_2atmpS1888 = _M0L1nS300 & _M0L4maskS298;
        int32_t _M0L5digitS301 = *(int32_t*)&_M0L6_2atmpS1888;
        int32_t _M0L6_2atmpS1885 = _M0L12digit__startS303 + _M0L6offsetS299;
        int32_t _M0L6_2atmpS1883 = _M0L6_2atmpS1885 - 1;
        int32_t _M0L6_2atmpS1884 =
          ((moonbit_string_t)moonbit_string_literal_113.data)[_M0L5digitS301];
        int32_t _M0L6_2atmpS1886;
        uint32_t _M0L6_2atmpS1887;
        _M0L6bufferS302[_M0L6_2atmpS1883] = _M0L6_2atmpS1884;
        _M0L6_2atmpS1886 = _M0L6offsetS299 - 1;
        _M0L6_2atmpS1887 = _M0L1nS300 >> (_M0L5shiftS297 & 31);
        _M0L6offsetS299 = _M0L6_2atmpS1886;
        _M0L1nS300 = _M0L6_2atmpS1887;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1896 = _M0L10total__lenS305 - _M0L12digit__startS303;
    int32_t _M0L6offsetS307 = _M0L6_2atmpS1896;
    uint32_t _M0L1nS308 = _M0L3numS306;
    while (1) {
      if (_M0L1nS308 > 0u) {
        uint32_t _M0L1qS309 = _M0L1nS308 / _M0L4baseS295;
        uint32_t _M0L6_2atmpS1895 = _M0L1qS309 * _M0L4baseS295;
        uint32_t _M0L6_2atmpS1894 = _M0L1nS308 - _M0L6_2atmpS1895;
        int32_t _M0L5digitS310 = *(int32_t*)&_M0L6_2atmpS1894;
        int32_t _M0L6_2atmpS1892 = _M0L12digit__startS303 + _M0L6offsetS307;
        int32_t _M0L6_2atmpS1890 = _M0L6_2atmpS1892 - 1;
        int32_t _M0L6_2atmpS1891 =
          ((moonbit_string_t)moonbit_string_literal_113.data)[_M0L5digitS310];
        int32_t _M0L6_2atmpS1893;
        _M0L6bufferS302[_M0L6_2atmpS1890] = _M0L6_2atmpS1891;
        _M0L6_2atmpS1893 = _M0L6offsetS307 - 1;
        _M0L6offsetS307 = _M0L6_2atmpS1893;
        _M0L1nS308 = _M0L1qS309;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB20int__to__string__hex(
  uint16_t* _M0L6bufferS289,
  uint32_t _M0L3numS294,
  int32_t _M0L12digit__startS290,
  int32_t _M0L10total__lenS293
) {
  int32_t _M0L6_2atmpS1880;
  int32_t _M0L6offsetS284;
  uint32_t _M0L1nS285;
  #line 29 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1880 = _M0L10total__lenS293 - _M0L12digit__startS290;
  _M0L6offsetS284 = _M0L6_2atmpS1880;
  _M0L1nS285 = _M0L3numS294;
  while (1) {
    if (_M0L6offsetS284 >= 2) {
      uint32_t _M0L6_2atmpS1877 = _M0L1nS285 & 255u;
      int32_t _M0L9byte__valS286 = *(int32_t*)&_M0L6_2atmpS1877;
      int32_t _M0L2hiS287 = _M0L9byte__valS286 / 16;
      int32_t _M0L2loS288 = _M0L9byte__valS286 % 16;
      int32_t _M0L6_2atmpS1871 = _M0L12digit__startS290 + _M0L6offsetS284;
      int32_t _M0L6_2atmpS1869 = _M0L6_2atmpS1871 - 2;
      int32_t _M0L6_2atmpS1870 =
        ((moonbit_string_t)moonbit_string_literal_113.data)[_M0L2hiS287];
      int32_t _M0L6_2atmpS1874;
      int32_t _M0L6_2atmpS1872;
      int32_t _M0L6_2atmpS1873;
      int32_t _M0L6_2atmpS1875;
      uint32_t _M0L6_2atmpS1876;
      _M0L6bufferS289[_M0L6_2atmpS1869] = _M0L6_2atmpS1870;
      _M0L6_2atmpS1874 = _M0L12digit__startS290 + _M0L6offsetS284;
      _M0L6_2atmpS1872 = _M0L6_2atmpS1874 - 1;
      _M0L6_2atmpS1873
      = ((moonbit_string_t)moonbit_string_literal_113.data)[
        _M0L2loS288
      ];
      _M0L6bufferS289[_M0L6_2atmpS1872] = _M0L6_2atmpS1873;
      _M0L6_2atmpS1875 = _M0L6offsetS284 - 2;
      _M0L6_2atmpS1876 = _M0L1nS285 >> 8;
      _M0L6offsetS284 = _M0L6_2atmpS1875;
      _M0L1nS285 = _M0L6_2atmpS1876;
      continue;
    } else if (_M0L6offsetS284 == 1) {
      uint32_t _M0L6_2atmpS1879 = _M0L1nS285 & 15u;
      int32_t _M0L6nibbleS292 = *(int32_t*)&_M0L6_2atmpS1879;
      int32_t _M0L6_2atmpS1878 =
        ((moonbit_string_t)moonbit_string_literal_113.data)[_M0L6nibbleS292];
      _M0L6bufferS289[_M0L12digit__startS290] = _M0L6_2atmpS1878;
    }
    break;
  }
  return 0;
}

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS273
) {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L7_2afuncS272;
  void* _M0L6resultS274;
  int64_t _M0L7_2abindS275;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS272 = _M0L4selfS273->$0;
  moonbit_incref(_M0L7_2afuncS272);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS274 = _M0L7_2afuncS272->code(_M0L7_2afuncS272);
  moonbit_decref(_M0L7_2afuncS272);
  _M0L7_2abindS275 = _M0L4selfS273->$1;
  switch (Moonbit_object_tag(_M0L6resultS274)) {
    case 1: {
      if (_M0L7_2abindS275 == 4294967296ll) {
        
      } else {
        int64_t _M0L7_2aSomeS276 = _M0L7_2abindS275;
        int32_t _M0L4_2anS277 = (int32_t)_M0L7_2aSomeS276;
        int64_t _M0L6_2atmpS1865;
        if (_M0L4_2anS277 > 0) {
          int32_t _M0L6_2atmpS1866 = _M0L4_2anS277 - 1;
          _M0L6_2atmpS1865 = (int64_t)_M0L6_2atmpS1866;
        } else {
          _M0L6_2atmpS1865
          = _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE;
        }
        _M0L4selfS273->$1 = _M0L6_2atmpS1865;
      }
      break;
    }
    default: {
      _M0L4selfS273->$1
      = _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE;
      break;
    }
  }
  return _M0L6resultS274;
}

int32_t _M0MPB4Iter4nextGcE(struct _M0TPB4IterGcE* _M0L4selfS279) {
  struct _M0TWEOc* _M0L7_2afuncS278;
  int32_t _M0L6resultS280;
  int64_t _M0L7_2abindS281;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS278 = _M0L4selfS279->$0;
  moonbit_incref(_M0L7_2afuncS278);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS280 = _M0L7_2afuncS278->code(_M0L7_2afuncS278);
  moonbit_decref(_M0L7_2afuncS278);
  _M0L7_2abindS281 = _M0L4selfS279->$1;
  if (_M0L6resultS280 == -1) {
    _M0L4selfS279->$1 = _M0MPB4Iter4nextN6constrS9981GcE;
  } else if (_M0L7_2abindS281 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS282 = _M0L7_2abindS281;
    int32_t _M0L4_2anS283 = (int32_t)_M0L7_2aSomeS282;
    int64_t _M0L6_2atmpS1867;
    if (_M0L4_2anS283 > 0) {
      int32_t _M0L6_2atmpS1868 = _M0L4_2anS283 - 1;
      _M0L6_2atmpS1867 = (int64_t)_M0L6_2atmpS1868;
    } else {
      _M0L6_2atmpS1867 = _M0MPB4Iter4nextN6constrS9980GcE;
    }
    _M0L4selfS279->$1 = _M0L6_2atmpS1867;
  }
  return _M0L6resultS280;
}

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void* _M0L4selfS271
) {
  struct _M0TPB13StringBuilder* _M0L6loggerS270;
  struct _M0TPB6Logger _M0L6_2atmpS1864;
  moonbit_string_t _result_3449;
  #line 165 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 166 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS270 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  moonbit_incref(_M0L6loggerS270);
  _M0L6_2atmpS1864
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L6loggerS270
  };
  #line 167 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB7FailurePB4Show6output(_M0L4selfS271, _M0L6_2atmpS1864);
  if (_M0L6_2atmpS1864.$1) {
    moonbit_decref(_M0L6_2atmpS1864.$1);
  }
  #line 168 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _result_3449 = _M0MPB13StringBuilder10to__string(_M0L6loggerS270);
  moonbit_decref(_M0L6loggerS270);
  return _result_3449;
}

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t _M0L4selfS267,
  struct _M0TPB6Logger _M0L6loggerS266
) {
  moonbit_string_t _M0L6_2atmpS1862;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1862 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS267);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS266.$0->$method_0(_M0L6loggerS266.$1, _M0L6_2atmpS1862);
  moonbit_decref(_M0L6_2atmpS1862);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t _M0L4selfS269,
  struct _M0TPB6Logger _M0L6loggerS268
) {
  moonbit_string_t _M0L6_2atmpS1863;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1863 = _M0IPC13int3IntPB4Show10to__string(_M0L4selfS269);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS268.$0->$method_0(_M0L6loggerS268.$1, _M0L6_2atmpS1863);
  moonbit_decref(_M0L6_2atmpS1863);
  return 0;
}

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView _M0L4selfS265
) {
  #line 99 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  return _M0L4selfS265.$1;
}

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView _M0L4selfS264
) {
  moonbit_string_t _M0L8_2afieldS3185;
  #line 92 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L8_2afieldS3185 = _M0L4selfS264.$0;
  moonbit_incref(_M0L8_2afieldS3185);
  return _M0L8_2afieldS3185;
}

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS260,
  moonbit_string_t _M0L5valueS261,
  int32_t _M0L5startS262,
  int32_t _M0L3lenS263
) {
  int32_t _M0L6_2atmpS1861;
  int64_t _M0L6_2atmpS1860;
  struct _M0TPC16string10StringView _M0L6_2atmpS1859;
  #line 122 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1861 = _M0L5startS262 + _M0L3lenS263;
  _M0L6_2atmpS1860 = (int64_t)_M0L6_2atmpS1861;
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1859
  = _M0MPC16string6String11sub_2einner(_M0L5valueS261, _M0L5startS262, _M0L6_2atmpS1860);
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L4selfS260, _M0L6_2atmpS1859);
  moonbit_decref(_M0L6_2atmpS1859.$0);
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t _M0L4selfS253,
  int32_t _M0L5startS259,
  int64_t _M0L3endS255
) {
  int32_t _M0L3lenS252;
  int32_t _M0L3endS254;
  int32_t _M0L5startS258;
  int32_t _if__result_3450;
  #line 755 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3lenS252 = Moonbit_array_length(_M0L4selfS253);
  if (_M0L3endS255 == 4294967296ll) {
    _M0L3endS254 = _M0L3lenS252;
  } else {
    int64_t _M0L7_2aSomeS256 = _M0L3endS255;
    int32_t _M0L6_2aendS257 = (int32_t)_M0L7_2aSomeS256;
    if (_M0L6_2aendS257 < 0) {
      _M0L3endS254 = _M0L3lenS252 + _M0L6_2aendS257;
    } else {
      _M0L3endS254 = _M0L6_2aendS257;
    }
  }
  if (_M0L5startS259 < 0) {
    _M0L5startS258 = _M0L3lenS252 + _M0L5startS259;
  } else {
    _M0L5startS258 = _M0L5startS259;
  }
  if (_M0L5startS258 >= 0) {
    if (_M0L5startS258 <= _M0L3endS254) {
      _if__result_3450 = _M0L3endS254 <= _M0L3lenS252;
    } else {
      _if__result_3450 = 0;
    }
  } else {
    _if__result_3450 = 0;
  }
  if (_if__result_3450) {
    if (_M0L5startS258 < _M0L3lenS252) {
      int32_t _M0L6_2atmpS1856 = _M0L4selfS253[_M0L5startS258];
      int32_t _M0L6_2atmpS1855;
      #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1855
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1856);
      if (!_M0L6_2atmpS1855) {
        
      } else {
        #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L3endS254 < _M0L3lenS252) {
      int32_t _M0L6_2atmpS1858 = _M0L4selfS253[_M0L3endS254];
      int32_t _M0L6_2atmpS1857;
      #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1857
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1858);
      if (!_M0L6_2atmpS1857) {
        
      } else {
        #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    moonbit_incref(_M0L4selfS253);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS253,
                                                 .$1 = _M0L5startS258,
                                                 .$2 = _M0L3endS254};
  } else {
    #line 763 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS251,
  struct _M0TPB4Show _M0L4showS250
) {
  struct _M0TPB6Logger _M0L6_2atmpS1854;
  #line 116 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS251);
  _M0L6_2atmpS1854
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS251
  };
  #line 117 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS250.$0->$method_0(_M0L4showS250.$1, _M0L6_2atmpS1854);
  if (_M0L6_2atmpS1854.$1) {
    moonbit_decref(_M0L6_2atmpS1854.$1);
  }
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS249,
  struct _M0TPB4Show _M0L4showS248
) {
  struct _M0TPB6Logger _M0L6_2atmpS1853;
  #line 111 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS249);
  _M0L6_2atmpS1853
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS249
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS248.$0->$method_0(_M0L4showS248.$1, _M0L6_2atmpS1853);
  if (_M0L6_2atmpS1853.$1) {
    moonbit_decref(_M0L6_2atmpS1853.$1);
  }
  return 0;
}

int32_t _M0FPB13finalize__acc(uint32_t _M0L3accS247) {
  uint32_t _M0L6_2atmpS1852;
  #line 444 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  #line 445 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1852 = _M0FPB14avalanche__acc(_M0L3accS247);
  return *(int32_t*)&_M0L6_2atmpS1852;
}

uint32_t _M0FPB14avalanche__acc(uint32_t _M0L3accS246) {
  uint32_t _M0Lm3accS245;
  uint32_t _M0L6_2atmpS1841;
  uint32_t _M0L6_2atmpS1843;
  uint32_t _M0L6_2atmpS1842;
  uint32_t _M0L6_2atmpS1844;
  uint32_t _M0L6_2atmpS1845;
  uint32_t _M0L6_2atmpS1847;
  uint32_t _M0L6_2atmpS1846;
  uint32_t _M0L6_2atmpS1848;
  uint32_t _M0L6_2atmpS1849;
  uint32_t _M0L6_2atmpS1851;
  uint32_t _M0L6_2atmpS1850;
  #line 449 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0Lm3accS245 = _M0L3accS246;
  _M0L6_2atmpS1841 = _M0Lm3accS245;
  _M0L6_2atmpS1843 = _M0Lm3accS245;
  _M0L6_2atmpS1842 = _M0L6_2atmpS1843 >> 15;
  _M0Lm3accS245 = _M0L6_2atmpS1841 ^ _M0L6_2atmpS1842;
  _M0L6_2atmpS1844 = _M0Lm3accS245;
  _M0Lm3accS245 = _M0L6_2atmpS1844 * 2246822519u;
  _M0L6_2atmpS1845 = _M0Lm3accS245;
  _M0L6_2atmpS1847 = _M0Lm3accS245;
  _M0L6_2atmpS1846 = _M0L6_2atmpS1847 >> 13;
  _M0Lm3accS245 = _M0L6_2atmpS1845 ^ _M0L6_2atmpS1846;
  _M0L6_2atmpS1848 = _M0Lm3accS245;
  _M0Lm3accS245 = _M0L6_2atmpS1848 * 3266489917u;
  _M0L6_2atmpS1849 = _M0Lm3accS245;
  _M0L6_2atmpS1851 = _M0Lm3accS245;
  _M0L6_2atmpS1850 = _M0L6_2atmpS1851 >> 16;
  _M0Lm3accS245 = _M0L6_2atmpS1849 ^ _M0L6_2atmpS1850;
  return _M0Lm3accS245;
}

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder* _M0L4selfS244,
  moonbit_string_t _M0L3strS243
) {
  int32_t _M0L8str__lenS242;
  int32_t _M0L3lenS1836;
  int32_t _M0L6_2atmpS1835;
  uint16_t* _M0L4dataS1837;
  int32_t _M0L3lenS1838;
  int32_t _M0L3lenS1840;
  int32_t _M0L6_2atmpS1839;
  #line 86 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L8str__lenS242 = Moonbit_array_length(_M0L3strS243);
  _M0L3lenS1836 = _M0L4selfS244->$1;
  _M0L6_2atmpS1835 = _M0L3lenS1836 + _M0L8str__lenS242;
  #line 88 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS244, _M0L6_2atmpS1835);
  _M0L4dataS1837 = _M0L4selfS244->$0;
  _M0L3lenS1838 = _M0L4selfS244->$1;
  moonbit_incref(_M0L4dataS1837);
  #line 89 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1837, _M0L3lenS1838, _M0L3strS243, 0, _M0L8str__lenS242);
  moonbit_decref(_M0L4dataS1837);
  _M0L3lenS1840 = _M0L4selfS244->$1;
  _M0L6_2atmpS1839 = _M0L3lenS1840 + _M0L8str__lenS242;
  _M0L4selfS244->$1 = _M0L6_2atmpS1839;
  return 0;
}

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t* _M0L4selfS238,
  int32_t _M0L11dst__offsetS241,
  moonbit_string_t _M0L3strS239,
  int32_t _M0L11str__offsetS234,
  int32_t _M0L3lenS235
) {
  int32_t _M0L16end__str__offsetS233;
  int32_t _M0L1iS236;
  int32_t _M0L1jS237;
  #line 71 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L16end__str__offsetS233 = _M0L11str__offsetS234 + _M0L3lenS235;
  _M0L1iS236 = _M0L11str__offsetS234;
  _M0L1jS237 = _M0L11dst__offsetS241;
  while (1) {
    if (_M0L1iS236 < _M0L16end__str__offsetS233) {
      int32_t _M0L6_2atmpS1832 = _M0L3strS239[_M0L1iS236];
      int32_t _M0L6_2atmpS1833;
      int32_t _M0L6_2atmpS1834;
      if (
        _M0L1jS237 < 0 || _M0L1jS237 >= Moonbit_array_length(_M0L4selfS238)
      ) {
        #line 80 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
        moonbit_panic();
      }
      _M0L4selfS238[_M0L1jS237] = _M0L6_2atmpS1832;
      _M0L6_2atmpS1833 = _M0L1iS236 + 1;
      _M0L6_2atmpS1834 = _M0L1jS237 + 1;
      _M0L1iS236 = _M0L6_2atmpS1833;
      _M0L1jS237 = _M0L6_2atmpS1834;
      continue;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t _M0L4selfS231,
  int32_t _M0L5quoteS232
) {
  struct _M0TPB13StringBuilder* _M0L3bufS230;
  int32_t _M0L6_2atmpS1831;
  struct _M0TPC16string10StringView _M0L6_2atmpS1829;
  struct _M0TPB6Logger _M0L6_2atmpS1830;
  moonbit_string_t _result_3452;
  #line 110 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 111 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L3bufS230 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS1831 = Moonbit_array_length(_M0L4selfS231);
  moonbit_incref(_M0L4selfS231);
  _M0L6_2atmpS1829
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS231, .$1 = 0, .$2 = _M0L6_2atmpS1831
  };
  moonbit_incref(_M0L3bufS230);
  _M0L6_2atmpS1830
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS230
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0MPC16string10StringView18escape__to_2einner(_M0L6_2atmpS1829, _M0L6_2atmpS1830, _M0L5quoteS232);
  moonbit_decref(_M0L6_2atmpS1829.$0);
  if (_M0L6_2atmpS1830.$1) {
    moonbit_decref(_M0L6_2atmpS1830.$1);
  }
  #line 113 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_3452 = _M0MPB13StringBuilder10to__string(_M0L3bufS230);
  moonbit_decref(_M0L3bufS230);
  return _result_3452;
}

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView _M0L4selfS222,
  struct _M0TPB6Logger _M0L6loggerS220,
  int32_t _M0L5quoteS219
) {
  int32_t _M0L3endS1827;
  int32_t _M0L5startS1828;
  int32_t _M0L3lenS221;
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS223;
  int32_t _M0L1iS224;
  int32_t _M0L3segS225;
  #line 144 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L5quoteS219) {
    #line 150 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS220.$0->$method_3(_M0L6loggerS220.$1, 34);
  }
  _M0L3endS1827 = _M0L4selfS222.$2;
  _M0L5startS1828 = _M0L4selfS222.$1;
  _M0L3lenS221 = _M0L3endS1827 - _M0L5startS1828;
  moonbit_incref(_M0L4selfS222.$0);
  if (_M0L6loggerS220.$1) {
    moonbit_incref(_M0L6loggerS220.$1);
  }
  _M0L6_2aenvS223
  = (struct _M0TURPC16string10StringViewRPB6LoggerE*)moonbit_malloc(sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE));
  Moonbit_object_header(_M0L6_2aenvS223)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 76, 0);
  _M0L6_2aenvS223->$0 = _M0L4selfS222;
  _M0L6_2aenvS223->$1 = _M0L6loggerS220;
  _M0L1iS224 = 0;
  _M0L3segS225 = 0;
  _2afor_226:;
  while (1) {
    moonbit_string_t _M0L3strS1824;
    int32_t _M0L5startS1826;
    int32_t _M0L6_2atmpS1825;
    int32_t _M0L4codeS227;
    int32_t _M0L1cS229;
    int32_t _M0L6_2atmpS1808;
    int32_t _M0L6_2atmpS1809;
    int32_t _M0L6_2atmpS1810;
    int32_t _tmp_3456;
    int32_t _tmp_3457;
    if (_M0L1iS224 >= _M0L3lenS221) {
      #line 160 "/home/developer/.moon/lib/core/builtin/show.mbt"
      _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
      moonbit_decref(_M0L6_2aenvS223);
      break;
    }
    _M0L3strS1824 = _M0L4selfS222.$0;
    _M0L5startS1826 = _M0L4selfS222.$1;
    _M0L6_2atmpS1825 = _M0L5startS1826 + _M0L1iS224;
    _M0L4codeS227 = _M0L3strS1824[_M0L6_2atmpS1825];
    switch (_M0L4codeS227) {
      case 34: {
        _M0L1cS229 = _M0L4codeS227;
        goto join_228;
        break;
      }
      
      case 92: {
        _M0L1cS229 = _M0L4codeS227;
        goto join_228;
        break;
      }
      
      case 10: {
        int32_t _M0L6_2atmpS1811;
        int32_t _M0L6_2atmpS1812;
        #line 172 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
        #line 173 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS220.$0->$method_0(_M0L6loggerS220.$1, (moonbit_string_t)moonbit_string_literal_114.data);
        _M0L6_2atmpS1811 = _M0L1iS224 + 1;
        _M0L6_2atmpS1812 = _M0L1iS224 + 1;
        _M0L1iS224 = _M0L6_2atmpS1811;
        _M0L3segS225 = _M0L6_2atmpS1812;
        goto _2afor_226;
        break;
      }
      
      case 13: {
        int32_t _M0L6_2atmpS1813;
        int32_t _M0L6_2atmpS1814;
        #line 177 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
        #line 178 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS220.$0->$method_0(_M0L6loggerS220.$1, (moonbit_string_t)moonbit_string_literal_115.data);
        _M0L6_2atmpS1813 = _M0L1iS224 + 1;
        _M0L6_2atmpS1814 = _M0L1iS224 + 1;
        _M0L1iS224 = _M0L6_2atmpS1813;
        _M0L3segS225 = _M0L6_2atmpS1814;
        goto _2afor_226;
        break;
      }
      
      case 8: {
        int32_t _M0L6_2atmpS1815;
        int32_t _M0L6_2atmpS1816;
        #line 182 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
        #line 183 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS220.$0->$method_0(_M0L6loggerS220.$1, (moonbit_string_t)moonbit_string_literal_116.data);
        _M0L6_2atmpS1815 = _M0L1iS224 + 1;
        _M0L6_2atmpS1816 = _M0L1iS224 + 1;
        _M0L1iS224 = _M0L6_2atmpS1815;
        _M0L3segS225 = _M0L6_2atmpS1816;
        goto _2afor_226;
        break;
      }
      
      case 9: {
        int32_t _M0L6_2atmpS1817;
        int32_t _M0L6_2atmpS1818;
        #line 187 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
        #line 188 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS220.$0->$method_0(_M0L6loggerS220.$1, (moonbit_string_t)moonbit_string_literal_117.data);
        _M0L6_2atmpS1817 = _M0L1iS224 + 1;
        _M0L6_2atmpS1818 = _M0L1iS224 + 1;
        _M0L1iS224 = _M0L6_2atmpS1817;
        _M0L3segS225 = _M0L6_2atmpS1818;
        goto _2afor_226;
        break;
      }
      default: {
        if (_M0L4codeS227 < 32) {
          int32_t _M0L6_2atmpS1820;
          moonbit_string_t _M0L6_2atmpS1819;
          int32_t _M0L6_2atmpS1821;
          int32_t _M0L6_2atmpS1822;
          #line 193 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
          #line 194 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS220.$0->$method_0(_M0L6loggerS220.$1, (moonbit_string_t)moonbit_string_literal_118.data);
          _M0L6_2atmpS1820 = _M0L4codeS227 & 0xff;
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6_2atmpS1819 = _M0MPC14byte4Byte7to__hex(_M0L6_2atmpS1820);
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS220.$0->$method_0(_M0L6loggerS220.$1, _M0L6_2atmpS1819);
          moonbit_decref(_M0L6_2atmpS1819);
          #line 196 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS220.$0->$method_3(_M0L6loggerS220.$1, 125);
          _M0L6_2atmpS1821 = _M0L1iS224 + 1;
          _M0L6_2atmpS1822 = _M0L1iS224 + 1;
          _M0L1iS224 = _M0L6_2atmpS1821;
          _M0L3segS225 = _M0L6_2atmpS1822;
          goto _2afor_226;
        } else {
          int32_t _M0L6_2atmpS1823 = _M0L1iS224 + 1;
          int32_t _tmp_3455 = _M0L3segS225;
          _M0L1iS224 = _M0L6_2atmpS1823;
          _M0L3segS225 = _tmp_3455;
          goto _2afor_226;
        }
        break;
      }
    }
    goto joinlet_3454;
    join_228:;
    #line 166 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS223, _M0L3segS225, _M0L1iS224);
    #line 167 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS220.$0->$method_3(_M0L6loggerS220.$1, 92);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1808 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS229);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS220.$0->$method_3(_M0L6loggerS220.$1, _M0L6_2atmpS1808);
    _M0L6_2atmpS1809 = _M0L1iS224 + 1;
    _M0L6_2atmpS1810 = _M0L1iS224 + 1;
    _M0L1iS224 = _M0L6_2atmpS1809;
    _M0L3segS225 = _M0L6_2atmpS1810;
    continue;
    joinlet_3454:;
    _tmp_3456 = _M0L1iS224;
    _tmp_3457 = _M0L3segS225;
    _M0L1iS224 = _tmp_3456;
    _M0L3segS225 = _tmp_3457;
    continue;
    break;
  }
  if (_M0L5quoteS219) {
    #line 204 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS220.$0->$method_3(_M0L6loggerS220.$1, 34);
  }
  return 0;
}

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS215,
  int32_t _M0L3segS218,
  int32_t _M0L1iS217
) {
  struct _M0TPB6Logger _M0L6loggerS214;
  struct _M0TPC16string10StringView _M0L4selfS216;
  #line 153 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6loggerS214 = _M0L6_2aenvS215->$1;
  _M0L4selfS216 = _M0L6_2aenvS215->$0;
  if (_M0L1iS217 > _M0L3segS218) {
    int64_t _M0L6_2atmpS1807 = (int64_t)_M0L1iS217;
    struct _M0TPC16string10StringView _M0L6_2atmpS1806;
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1806
    = _M0MPC16string10StringView11sub_2einner(_M0L4selfS216, _M0L3segS218, _M0L6_2atmpS1807);
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS214.$0->$method_2(_M0L6loggerS214.$1, _M0L6_2atmpS1806);
    moonbit_decref(_M0L6_2atmpS1806.$0);
  }
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView _M0L4selfS207,
  int32_t _M0L5startS213,
  int64_t _M0L3endS209
) {
  moonbit_string_t _M0L3strS1805;
  int32_t _M0L8str__lenS206;
  int32_t _M0L8abs__endS208;
  int32_t _M0L10abs__startS212;
  int32_t _M0L5startS1793;
  int32_t _if__result_3458;
  #line 814 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1805 = _M0L4selfS207.$0;
  _M0L8str__lenS206 = Moonbit_array_length(_M0L3strS1805);
  if (_M0L3endS209 == 4294967296ll) {
    _M0L8abs__endS208 = _M0L4selfS207.$2;
  } else {
    int64_t _M0L7_2aSomeS210 = _M0L3endS209;
    int32_t _M0L6_2aendS211 = (int32_t)_M0L7_2aSomeS210;
    if (_M0L6_2aendS211 < 0) {
      int32_t _M0L3endS1803 = _M0L4selfS207.$2;
      _M0L8abs__endS208 = _M0L3endS1803 + _M0L6_2aendS211;
    } else {
      int32_t _M0L5startS1804 = _M0L4selfS207.$1;
      _M0L8abs__endS208 = _M0L5startS1804 + _M0L6_2aendS211;
    }
  }
  if (_M0L5startS213 < 0) {
    int32_t _M0L3endS1801 = _M0L4selfS207.$2;
    _M0L10abs__startS212 = _M0L3endS1801 + _M0L5startS213;
  } else {
    int32_t _M0L5startS1802 = _M0L4selfS207.$1;
    _M0L10abs__startS212 = _M0L5startS1802 + _M0L5startS213;
  }
  _M0L5startS1793 = _M0L4selfS207.$1;
  if (_M0L10abs__startS212 >= _M0L5startS1793) {
    if (_M0L10abs__startS212 <= _M0L8abs__endS208) {
      int32_t _M0L3endS1792 = _M0L4selfS207.$2;
      _if__result_3458 = _M0L8abs__endS208 <= _M0L3endS1792;
    } else {
      _if__result_3458 = 0;
    }
  } else {
    _if__result_3458 = 0;
  }
  if (_if__result_3458) {
    moonbit_string_t _M0L3strS1800;
    if (_M0L10abs__startS212 < _M0L8str__lenS206) {
      moonbit_string_t _M0L3strS1796 = _M0L4selfS207.$0;
      int32_t _M0L6_2atmpS1795 = _M0L3strS1796[_M0L10abs__startS212];
      int32_t _M0L6_2atmpS1794;
      #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1794
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1795);
      if (!_M0L6_2atmpS1794) {
        
      } else {
        #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L8abs__endS208 < _M0L8str__lenS206) {
      moonbit_string_t _M0L3strS1799 = _M0L4selfS207.$0;
      int32_t _M0L6_2atmpS1798 = _M0L3strS1799[_M0L8abs__endS208];
      int32_t _M0L6_2atmpS1797;
      #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1797
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1798);
      if (!_M0L6_2atmpS1797) {
        
      } else {
        #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    _M0L3strS1800 = _M0L4selfS207.$0;
    moonbit_incref(_M0L3strS1800);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1800,
                                                 .$1 = _M0L10abs__startS212,
                                                 .$2 = _M0L8abs__endS208};
  } else {
    #line 834 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t _M0L1bS205) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS204;
  int32_t _M0L6_2atmpS1789;
  int32_t _M0L6_2atmpS1788;
  int32_t _M0L6_2atmpS1791;
  int32_t _M0L6_2atmpS1790;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1787;
  moonbit_string_t _result_3459;
  #line 74 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L7_2aselfS204 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1789 = _M0IPC14byte4BytePB3Div3div(_M0L1bS205, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1788
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1789);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS204, _M0L6_2atmpS1788);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1791 = _M0IPC14byte4BytePB3Mod3mod(_M0L1bS205, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1790
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1791);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS204, _M0L6_2atmpS1790);
  _M0L6_2atmpS1787 = _M0L7_2aselfS204;
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_3459 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1787);
  moonbit_decref(_M0L6_2atmpS1787);
  return _result_3459;
}

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t _M0L1iS203) {
  #line 75 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L1iS203 < 10) {
    int32_t _M0L6_2atmpS1784;
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1784 = _M0IPC14byte4BytePB3Add3add(_M0L1iS203, 48);
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1784);
  } else {
    int32_t _M0L6_2atmpS1786;
    int32_t _M0L6_2atmpS1785;
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1786 = _M0IPC14byte4BytePB3Add3add(_M0L1iS203, 97);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1785 = _M0IPC14byte4BytePB3Sub3sub(_M0L6_2atmpS1786, 10);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1785);
  }
}

int32_t _M0IPC14byte4BytePB3Sub3sub(
  int32_t _M0L4selfS201,
  int32_t _M0L4thatS202
) {
  int32_t _M0L6_2atmpS1782;
  int32_t _M0L6_2atmpS1783;
  int32_t _M0L6_2atmpS1781;
  #line 120 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1782 = (int32_t)_M0L4selfS201;
  _M0L6_2atmpS1783 = (int32_t)_M0L4thatS202;
  _M0L6_2atmpS1781 = _M0L6_2atmpS1782 - _M0L6_2atmpS1783;
  return _M0L6_2atmpS1781 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Mod3mod(
  int32_t _M0L4selfS199,
  int32_t _M0L4thatS200
) {
  int32_t _M0L6_2atmpS1779;
  int32_t _M0L6_2atmpS1780;
  int32_t _M0L6_2atmpS1778;
  #line 67 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1779 = (int32_t)_M0L4selfS199;
  _M0L6_2atmpS1780 = (int32_t)_M0L4thatS200;
  _M0L6_2atmpS1778 = _M0L6_2atmpS1779 % _M0L6_2atmpS1780;
  return _M0L6_2atmpS1778 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Div3div(
  int32_t _M0L4selfS197,
  int32_t _M0L4thatS198
) {
  int32_t _M0L6_2atmpS1776;
  int32_t _M0L6_2atmpS1777;
  int32_t _M0L6_2atmpS1775;
  #line 62 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1776 = (int32_t)_M0L4selfS197;
  _M0L6_2atmpS1777 = (int32_t)_M0L4thatS198;
  _M0L6_2atmpS1775 = _M0L6_2atmpS1776 / _M0L6_2atmpS1777;
  return _M0L6_2atmpS1775 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Add3add(
  int32_t _M0L4selfS195,
  int32_t _M0L4thatS196
) {
  int32_t _M0L6_2atmpS1773;
  int32_t _M0L6_2atmpS1774;
  int32_t _M0L6_2atmpS1772;
  #line 106 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1773 = (int32_t)_M0L4selfS195;
  _M0L6_2atmpS1774 = (int32_t)_M0L4thatS196;
  _M0L6_2atmpS1772 = _M0L6_2atmpS1773 + _M0L6_2atmpS1774;
  return _M0L6_2atmpS1772 & 0xff;
}

int32_t _M0MPC16string6String16unsafe__char__at(
  moonbit_string_t _M0L4selfS192,
  int32_t _M0L5indexS193
) {
  int32_t _M0L2c1S191;
  #line 102 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
  _M0L2c1S191 = _M0L4selfS192[_M0L5indexS193];
  #line 105 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
  if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S191)) {
    int32_t _M0L6_2atmpS1771 = _M0L5indexS193 + 1;
    int32_t _M0L2c2S194 = _M0L4selfS192[_M0L6_2atmpS1771];
    int32_t _M0L6_2atmpS1769 = (int32_t)_M0L2c1S191;
    int32_t _M0L6_2atmpS1770 = (int32_t)_M0L2c2S194;
    #line 107 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
    return _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1769, _M0L6_2atmpS1770);
  } else {
    #line 109 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
    return _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S191);
  }
}

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t _M0L4selfS190) {
  int32_t _M0L6_2atmpS1768;
  #line 68 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  _M0L6_2atmpS1768 = (int32_t)_M0L4selfS190;
  return _M0L6_2atmpS1768;
}

int32_t _M0FPB32code__point__of__surrogate__pair(
  int32_t _M0L7leadingS188,
  int32_t _M0L8trailingS189
) {
  int32_t _M0L6_2atmpS1767;
  int32_t _M0L6_2atmpS1766;
  int32_t _M0L6_2atmpS1765;
  int32_t _M0L6_2atmpS1764;
  int32_t _M0L6_2atmpS1763;
  #line 40 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1767 = _M0L7leadingS188 - 55296;
  _M0L6_2atmpS1766 = _M0L6_2atmpS1767 * 1024;
  _M0L6_2atmpS1765 = _M0L6_2atmpS1766 + _M0L8trailingS189;
  _M0L6_2atmpS1764 = _M0L6_2atmpS1765 - 56320;
  _M0L6_2atmpS1763 = _M0L6_2atmpS1764 + 65536;
  return _M0L6_2atmpS1763;
}

int32_t _M0MPC16string6String20char__length_2einner(
  moonbit_string_t _M0L4selfS181,
  int32_t _M0L13start__offsetS182,
  int64_t _M0L11end__offsetS179
) {
  int32_t _M0L11end__offsetS178;
  int32_t _if__result_3460;
  #line 60 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L11end__offsetS179 == 4294967296ll) {
    _M0L11end__offsetS178 = Moonbit_array_length(_M0L4selfS181);
  } else {
    int64_t _M0L7_2aSomeS180 = _M0L11end__offsetS179;
    _M0L11end__offsetS178 = (int32_t)_M0L7_2aSomeS180;
  }
  if (_M0L13start__offsetS182 >= 0) {
    if (_M0L13start__offsetS182 <= _M0L11end__offsetS178) {
      int32_t _M0L6_2atmpS1756 = Moonbit_array_length(_M0L4selfS181);
      _if__result_3460 = _M0L11end__offsetS178 <= _M0L6_2atmpS1756;
    } else {
      _if__result_3460 = 0;
    }
  } else {
    _if__result_3460 = 0;
  }
  if (_if__result_3460) {
    int32_t _M0L12utf16__indexS183 = _M0L13start__offsetS182;
    int32_t _M0L11char__countS184 = 0;
    while (1) {
      if (_M0L12utf16__indexS183 < _M0L11end__offsetS178) {
        int32_t _M0L2c1S185 = _M0L4selfS181[_M0L12utf16__indexS183];
        int32_t _if__result_3462;
        int32_t _M0L6_2atmpS1761;
        int32_t _M0L6_2atmpS1762;
        #line 76 "/home/developer/.moon/lib/core/builtin/string.mbt"
        if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S185)) {
          int32_t _M0L6_2atmpS1757 = _M0L12utf16__indexS183 + 1;
          _if__result_3462 = _M0L6_2atmpS1757 < _M0L11end__offsetS178;
        } else {
          _if__result_3462 = 0;
        }
        if (_if__result_3462) {
          int32_t _M0L6_2atmpS1760 = _M0L12utf16__indexS183 + 1;
          int32_t _M0L2c2S186 = _M0L4selfS181[_M0L6_2atmpS1760];
          #line 78 "/home/developer/.moon/lib/core/builtin/string.mbt"
          if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S186)) {
            int32_t _M0L6_2atmpS1758 = _M0L12utf16__indexS183 + 2;
            int32_t _M0L6_2atmpS1759 = _M0L11char__countS184 + 1;
            _M0L12utf16__indexS183 = _M0L6_2atmpS1758;
            _M0L11char__countS184 = _M0L6_2atmpS1759;
            continue;
          } else {
            #line 81 "/home/developer/.moon/lib/core/builtin/string.mbt"
            _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_119.data);
          }
        }
        _M0L6_2atmpS1761 = _M0L12utf16__indexS183 + 1;
        _M0L6_2atmpS1762 = _M0L11char__countS184 + 1;
        _M0L12utf16__indexS183 = _M0L6_2atmpS1761;
        _M0L11char__countS184 = _M0L6_2atmpS1762;
        continue;
      } else {
        return _M0L11char__countS184;
      }
      break;
    }
  } else {
    #line 70 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0FPC15abort5abortGiE((moonbit_string_t)moonbit_string_literal_120.data);
  }
}

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t _M0L4selfS177) {
  #line 45 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS177 >= 56320) {
    return _M0L4selfS177 <= 57343;
  } else {
    return 0;
  }
}

int32_t _M0MPC16uint166UInt1622is__leading__surrogate(int32_t _M0L4selfS176) {
  #line 28 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS176 >= 55296) {
    return _M0L4selfS176 <= 56319;
  } else {
    return 0;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder* _M0L4selfS174,
  int32_t _M0L2chS173
) {
  uint32_t _M0L4codeS172;
  #line 95 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  #line 96 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4codeS172 = _M0MPC14char4Char8to__uint(_M0L2chS173);
  if (_M0L4codeS172 <= 65535u) {
    int32_t _M0L3lenS1735 = _M0L4selfS174->$1;
    int32_t _M0L6_2atmpS1734 = _M0L3lenS1735 + 1;
    uint16_t* _M0L4dataS1736;
    int32_t _M0L3lenS1737;
    int32_t _M0L6_2atmpS1738;
    int32_t _M0L3lenS1740;
    int32_t _M0L6_2atmpS1739;
    #line 98 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS174, _M0L6_2atmpS1734);
    _M0L4dataS1736 = _M0L4selfS174->$0;
    _M0L3lenS1737 = _M0L4selfS174->$1;
    moonbit_incref(_M0L4dataS1736);
    #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1738 = _M0MPC14uint4UInt10to__uint16(_M0L4codeS172);
    if (
      _M0L3lenS1737 < 0
      || _M0L3lenS1737 >= Moonbit_array_length(_M0L4dataS1736)
    ) {
      #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1736[_M0L3lenS1737] = _M0L6_2atmpS1738;
    moonbit_decref(_M0L4dataS1736);
    _M0L3lenS1740 = _M0L4selfS174->$1;
    _M0L6_2atmpS1739 = _M0L3lenS1740 + 1;
    _M0L4selfS174->$1 = _M0L6_2atmpS1739;
  } else if (_M0L4codeS172 <= 1114111u) {
    int32_t _M0L3lenS1742 = _M0L4selfS174->$1;
    int32_t _M0L6_2atmpS1741 = _M0L3lenS1742 + 2;
    uint32_t _M0L4codeS175;
    uint16_t* _M0L4dataS1743;
    int32_t _M0L3lenS1744;
    uint32_t _M0L6_2atmpS1747;
    uint32_t _M0L6_2atmpS1746;
    int32_t _M0L6_2atmpS1745;
    uint16_t* _M0L4dataS1748;
    int32_t _M0L3lenS1753;
    int32_t _M0L6_2atmpS1749;
    uint32_t _M0L6_2atmpS1752;
    uint32_t _M0L6_2atmpS1751;
    int32_t _M0L6_2atmpS1750;
    int32_t _M0L3lenS1755;
    int32_t _M0L6_2atmpS1754;
    #line 102 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS174, _M0L6_2atmpS1741);
    _M0L4codeS175 = _M0L4codeS172 - 65536u;
    _M0L4dataS1743 = _M0L4selfS174->$0;
    _M0L3lenS1744 = _M0L4selfS174->$1;
    _M0L6_2atmpS1747 = _M0L4codeS175 >> 10;
    _M0L6_2atmpS1746 = 55296u + _M0L6_2atmpS1747;
    moonbit_incref(_M0L4dataS1743);
    #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1745 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1746);
    if (
      _M0L3lenS1744 < 0
      || _M0L3lenS1744 >= Moonbit_array_length(_M0L4dataS1743)
    ) {
      #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1743[_M0L3lenS1744] = _M0L6_2atmpS1745;
    moonbit_decref(_M0L4dataS1743);
    _M0L4dataS1748 = _M0L4selfS174->$0;
    _M0L3lenS1753 = _M0L4selfS174->$1;
    _M0L6_2atmpS1749 = _M0L3lenS1753 + 1;
    _M0L6_2atmpS1752 = _M0L4codeS175 & 1023u;
    _M0L6_2atmpS1751 = 56320u + _M0L6_2atmpS1752;
    moonbit_incref(_M0L4dataS1748);
    #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1750 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1751);
    if (
      _M0L6_2atmpS1749 < 0
      || _M0L6_2atmpS1749 >= Moonbit_array_length(_M0L4dataS1748)
    ) {
      #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1748[_M0L6_2atmpS1749] = _M0L6_2atmpS1750;
    moonbit_decref(_M0L4dataS1748);
    _M0L3lenS1755 = _M0L4selfS174->$1;
    _M0L6_2atmpS1754 = _M0L3lenS1755 + 2;
    _M0L4selfS174->$1 = _M0L6_2atmpS1754;
  } else {
    #line 108 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_121.data);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder* _M0L4selfS166,
  int32_t _M0L8requiredS167
) {
  uint16_t* _M0L4dataS1733;
  int32_t _M0L12current__lenS165;
  int32_t _M0L13enough__spaceS168;
  int32_t _M0L13enough__spaceS169;
  uint16_t* _M0L4dataS1729;
  int32_t _M0L6_2atmpS1730;
  int32_t _M0L3lenS1731;
  uint16_t* _M0L9new__dataS171;
  uint16_t* _M0L6_2aoldS3195;
  #line 46 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4dataS1733 = _M0L4selfS166->$0;
  _M0L12current__lenS165 = Moonbit_array_length(_M0L4dataS1733);
  if (_M0L8requiredS167 <= _M0L12current__lenS165) {
    return 0;
  }
  _M0L13enough__spaceS169 = _M0L12current__lenS165;
  while (1) {
    if (_M0L13enough__spaceS169 < _M0L8requiredS167) {
      int32_t _M0L6_2atmpS1732 = _M0L13enough__spaceS169 * 2;
      _M0L13enough__spaceS169 = _M0L6_2atmpS1732;
      continue;
    } else {
      _M0L13enough__spaceS168 = _M0L13enough__spaceS169;
    }
    break;
  }
  _M0L4dataS1729 = _M0L4selfS166->$0;
  moonbit_incref(_M0L4dataS1729);
  #line 64 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1730 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L3lenS1731 = _M0L4selfS166->$1;
  #line 61 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L9new__dataS171
  = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1729, _M0L13enough__spaceS168, _M0L6_2atmpS1730, _M0L3lenS1731, 0, 0);
  moonbit_decref(_M0L4dataS1729);
  _M0L6_2aoldS3195 = _M0L4selfS166->$0;
  moonbit_decref(_M0L6_2aoldS3195);
  _M0L4selfS166->$0 = _M0L9new__dataS171;
  return 0;
}

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t _M0L4selfS164) {
  int32_t _M0L6_2atmpS1728;
  #line 2676 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1728 = *(int32_t*)&_M0L4selfS164;
  return (uint16_t)_M0L6_2atmpS1728;
}

uint32_t _M0MPC14char4Char8to__uint(int32_t _M0L4selfS163) {
  int32_t _M0L6_2atmpS1727;
  #line 1254 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1727 = _M0L4selfS163;
  return *(uint32_t*)&_M0L6_2atmpS1727;
}

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder* _M0L4selfS161
) {
  int32_t _M0L3lenS1719;
  uint16_t* _M0L4dataS1721;
  int32_t _M0L6_2atmpS1720;
  #line 148 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3lenS1719 = _M0L4selfS161->$1;
  _M0L4dataS1721 = _M0L4selfS161->$0;
  _M0L6_2atmpS1720 = Moonbit_array_length(_M0L4dataS1721);
  if (_M0L3lenS1719 == _M0L6_2atmpS1720) {
    uint16_t* _M0L4dataS1722 = _M0L4selfS161->$0;
    moonbit_incref(_M0L4dataS1722);
    return _M0L4dataS1722;
  } else {
    uint16_t* _M0L4dataS1723 = _M0L4selfS161->$0;
    int32_t _M0L3lenS1724 = _M0L4selfS161->$1;
    int32_t _M0L6_2atmpS1725;
    int32_t _M0L3lenS1726;
    uint16_t* _M0L4dataS162;
    moonbit_incref(_M0L4dataS1723);
    #line 155 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1725 = _M0IPC16uint166UInt16PB7Default7default();
    _M0L3lenS1726 = _M0L4selfS161->$1;
    #line 152 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L4dataS162
    = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1723, _M0L3lenS1724, _M0L6_2atmpS1725, _M0L3lenS1726, 0, 0);
    moonbit_decref(_M0L4dataS1723);
    return _M0L4dataS162;
  }
}

int32_t _M0IPC16uint166UInt16PB7Default7default() {
  #line 176 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  return 0;
}

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t* _M0L3srcS158,
  int32_t _M0L13allocate__lenS154,
  int32_t _M0L4initS159,
  int32_t _M0L3lenS155,
  int32_t _M0L11src__offsetS156,
  int32_t _M0L11dst__offsetS157
) {
  int32_t _if__result_3464;
  #line 97 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L13allocate__lenS154 >= 0) {
    if (_M0L3lenS155 >= 0) {
      if (_M0L11src__offsetS156 >= 0) {
        if (_M0L11dst__offsetS157 >= 0) {
          int32_t _M0L6_2atmpS1715 = _M0L11src__offsetS156 + _M0L3lenS155;
          int32_t _M0L6_2atmpS1716 = Moonbit_array_length(_M0L3srcS158);
          if (_M0L6_2atmpS1715 <= _M0L6_2atmpS1716) {
            int32_t _M0L6_2atmpS1714 = _M0L11dst__offsetS157 + _M0L3lenS155;
            _if__result_3464 = _M0L6_2atmpS1714 <= _M0L13allocate__lenS154;
          } else {
            _if__result_3464 = 0;
          }
        } else {
          _if__result_3464 = 0;
        }
      } else {
        _if__result_3464 = 0;
      }
    } else {
      _if__result_3464 = 0;
    }
  } else {
    _if__result_3464 = 0;
  }
  if (_if__result_3464) {
    moonbit_incref(_M0L3srcS158);
    #line 115 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    return _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(_M0L3srcS158, _M0L13allocate__lenS154, _M0L4initS159, _M0L11src__offsetS156, _M0L11dst__offsetS157, _M0L3lenS155);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS160;
    int32_t _M0L6_2atmpS1718;
    moonbit_string_t _M0L6_2atmpS1717;
    uint16_t* _result_3465;
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L18_2astring__builderS160
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS160, (moonbit_string_t)moonbit_string_literal_122.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS160, _M0L13allocate__lenS154);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS160, (moonbit_string_t)moonbit_string_literal_123.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS160, _M0L11src__offsetS156);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS160, (moonbit_string_t)moonbit_string_literal_124.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS160, _M0L11dst__offsetS157);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS160, (moonbit_string_t)moonbit_string_literal_125.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS160, _M0L3lenS155);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS160, (moonbit_string_t)moonbit_string_literal_126.data);
    _M0L6_2atmpS1718 = Moonbit_array_length(_M0L3srcS158);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS160, _M0L6_2atmpS1718);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L6_2atmpS1717
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS160);
    moonbit_decref(_M0L18_2astring__builderS160);
    #line 111 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _result_3465 = _M0FPC15abort5abortGAkE(_M0L6_2atmpS1717);
    moonbit_decref(_M0L6_2atmpS1717);
    return _result_3465;
  }
}

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t* _M0L3srcS151,
  int32_t _M0L13allocate__lenS148,
  int32_t _M0L4initS149,
  int32_t _M0L11src__offsetS152,
  int32_t _M0L11dst__offsetS150,
  int32_t _M0L9blit__lenS153
) {
  uint16_t* _M0L3dstS147;
  #line 79 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  _M0L3dstS147
  = (uint16_t*)moonbit_make_string(_M0L13allocate__lenS148, _M0L4initS149);
  moonbit_incref(_M0L3dstS147);
  #line 90 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS147, _M0L11dst__offsetS150, _M0L3srcS151, _M0L11src__offsetS152, _M0L9blit__lenS153, sizeof(uint16_t));
  return _M0L3dstS147;
}

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t _M0L10size__hintS145
) {
  int32_t _M0L7initialS144;
  uint16_t* _M0L4dataS146;
  struct _M0TPB13StringBuilder* _block_3466;
  #line 32 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  if (_M0L10size__hintS145 < 1) {
    _M0L7initialS144 = 1;
  } else {
    int32_t _M0L6_2atmpS1713 = _M0L10size__hintS145 + 1;
    _M0L7initialS144 = _M0L6_2atmpS1713 / 2;
  }
  _M0L4dataS146 = (uint16_t*)moonbit_make_string(_M0L7initialS144, 0);
  _block_3466
  = (struct _M0TPB13StringBuilder*)moonbit_malloc(sizeof(struct _M0TPB13StringBuilder));
  Moonbit_object_header(_block_3466)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 81, 0);
  _block_3466->$0 = _M0L4dataS146;
  _block_3466->$1 = 0;
  return _block_3466;
}

int32_t _M0MPC14byte4Byte8to__char(int32_t _M0L4selfS143) {
  int32_t _M0L6_2atmpS1712;
  #line 1868 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1712 = (int32_t)_M0L4selfS143;
  return _M0L6_2atmpS1712;
}

int32_t _M0MPC13int3Int16unsafe__to__char(int32_t _M0L4selfS142) {
  #line 1532 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  return _M0L4selfS142;
}

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t* _M0L3srcS122,
  int32_t _M0L13allocate__lenS118,
  int32_t _M0L3lenS119,
  int32_t _M0L11src__offsetS120,
  int32_t _M0L11dst__offsetS121
) {
  int32_t _if__result_3467;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS118 >= 0) {
    if (_M0L3lenS119 >= 0) {
      if (_M0L11src__offsetS120 >= 0) {
        if (_M0L11dst__offsetS121 >= 0) {
          int32_t _M0L6_2atmpS1693 = _M0L11src__offsetS120 + _M0L3lenS119;
          int32_t _M0L6_2atmpS1694;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1694
          = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS122);
          if (_M0L6_2atmpS1693 <= _M0L6_2atmpS1694) {
            int32_t _M0L6_2atmpS1692 = _M0L11dst__offsetS121 + _M0L3lenS119;
            _if__result_3467 = _M0L6_2atmpS1692 <= _M0L13allocate__lenS118;
          } else {
            _if__result_3467 = 0;
          }
        } else {
          _if__result_3467 = 0;
        }
      } else {
        _if__result_3467 = 0;
      }
    } else {
      _if__result_3467 = 0;
    }
  } else {
    _if__result_3467 = 0;
  }
  if (_if__result_3467) {
    moonbit_incref(_M0L3srcS122);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (moonbit_string_t*)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS118, (moonbit_string_t)moonbit_string_literal_0.data, _M0L3srcS122, _M0L11src__offsetS120, _M0L11dst__offsetS121, _M0L3lenS119);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS123;
    int32_t _M0L6_2atmpS1696;
    moonbit_string_t _M0L6_2atmpS1695;
    moonbit_string_t* _result_3468;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS123
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS123, (moonbit_string_t)moonbit_string_literal_122.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS123, _M0L13allocate__lenS118);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS123, (moonbit_string_t)moonbit_string_literal_123.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS123, _M0L11src__offsetS120);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS123, (moonbit_string_t)moonbit_string_literal_124.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS123, _M0L11dst__offsetS121);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS123, (moonbit_string_t)moonbit_string_literal_125.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS123, _M0L3lenS119);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS123, (moonbit_string_t)moonbit_string_literal_126.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1696 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS122);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS123, _M0L6_2atmpS1696);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1695
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS123);
    moonbit_decref(_M0L18_2astring__builderS123);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3468
    = _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(_M0L6_2atmpS1695);
    moonbit_decref(_M0L6_2atmpS1695);
    return _result_3468;
  }
}

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE** _M0L3srcS128,
  int32_t _M0L13allocate__lenS124,
  int32_t _M0L3lenS125,
  int32_t _M0L11src__offsetS126,
  int32_t _M0L11dst__offsetS127
) {
  int32_t _if__result_3469;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS124 >= 0) {
    if (_M0L3lenS125 >= 0) {
      if (_M0L11src__offsetS126 >= 0) {
        if (_M0L11dst__offsetS127 >= 0) {
          int32_t _M0L6_2atmpS1698 = _M0L11src__offsetS126 + _M0L3lenS125;
          int32_t _M0L6_2atmpS1699;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1699
          = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS128);
          if (_M0L6_2atmpS1698 <= _M0L6_2atmpS1699) {
            int32_t _M0L6_2atmpS1697 = _M0L11dst__offsetS127 + _M0L3lenS125;
            _if__result_3469 = _M0L6_2atmpS1697 <= _M0L13allocate__lenS124;
          } else {
            _if__result_3469 = 0;
          }
        } else {
          _if__result_3469 = 0;
        }
      } else {
        _if__result_3469 = 0;
      }
    } else {
      _if__result_3469 = 0;
    }
  } else {
    _if__result_3469 = 0;
  }
  if (_if__result_3469) {
    moonbit_incref(_M0L3srcS128);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TUsiE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS124, 0, _M0L3srcS128, _M0L11src__offsetS126, _M0L11dst__offsetS127, _M0L3lenS125);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS129;
    int32_t _M0L6_2atmpS1701;
    moonbit_string_t _M0L6_2atmpS1700;
    struct _M0TUsiE** _result_3470;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS129
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS129, (moonbit_string_t)moonbit_string_literal_122.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS129, _M0L13allocate__lenS124);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS129, (moonbit_string_t)moonbit_string_literal_123.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS129, _M0L11src__offsetS126);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS129, (moonbit_string_t)moonbit_string_literal_124.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS129, _M0L11dst__offsetS127);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS129, (moonbit_string_t)moonbit_string_literal_125.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS129, _M0L3lenS125);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS129, (moonbit_string_t)moonbit_string_literal_126.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1701 = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS128);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS129, _M0L6_2atmpS1701);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1700
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS129);
    moonbit_decref(_M0L18_2astring__builderS129);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3470
    = _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(_M0L6_2atmpS1700);
    moonbit_decref(_M0L6_2atmpS1700);
    return _result_3470;
  }
}

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3srcS134,
  int32_t _M0L13allocate__lenS130,
  int32_t _M0L3lenS131,
  int32_t _M0L11src__offsetS132,
  int32_t _M0L11dst__offsetS133
) {
  int32_t _if__result_3471;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS130 >= 0) {
    if (_M0L3lenS131 >= 0) {
      if (_M0L11src__offsetS132 >= 0) {
        if (_M0L11dst__offsetS133 >= 0) {
          int32_t _M0L6_2atmpS1703 = _M0L11src__offsetS132 + _M0L3lenS131;
          int32_t _M0L6_2atmpS1704;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1704
          = _M0MPB18UninitializedArray6lengthGRP46biolab8bio__seq3cmd4main4CaseE(_M0L3srcS134);
          if (_M0L6_2atmpS1703 <= _M0L6_2atmpS1704) {
            int32_t _M0L6_2atmpS1702 = _M0L11dst__offsetS133 + _M0L3lenS131;
            _if__result_3471 = _M0L6_2atmpS1702 <= _M0L13allocate__lenS130;
          } else {
            _if__result_3471 = 0;
          }
        } else {
          _if__result_3471 = 0;
        }
      } else {
        _if__result_3471 = 0;
      }
    } else {
      _if__result_3471 = 0;
    }
  } else {
    _if__result_3471 = 0;
  }
  if (_if__result_3471) {
    moonbit_incref(_M0L3srcS134);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TP46biolab8bio__seq3cmd4main4Case**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS130, 0, _M0L3srcS134, _M0L11src__offsetS132, _M0L11dst__offsetS133, _M0L3lenS131);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS135;
    int32_t _M0L6_2atmpS1706;
    moonbit_string_t _M0L6_2atmpS1705;
    struct _M0TP46biolab8bio__seq3cmd4main4Case** _result_3472;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS135
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS135, (moonbit_string_t)moonbit_string_literal_122.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS135, _M0L13allocate__lenS130);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS135, (moonbit_string_t)moonbit_string_literal_123.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS135, _M0L11src__offsetS132);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS135, (moonbit_string_t)moonbit_string_literal_124.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS135, _M0L11dst__offsetS133);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS135, (moonbit_string_t)moonbit_string_literal_125.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS135, _M0L3lenS131);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS135, (moonbit_string_t)moonbit_string_literal_126.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1706
    = _M0MPB18UninitializedArray6lengthGRP46biolab8bio__seq3cmd4main4CaseE(_M0L3srcS134);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS135, _M0L6_2atmpS1706);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1705
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS135);
    moonbit_decref(_M0L18_2astring__builderS135);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3472
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRP46biolab8bio__seq3cmd4main4CaseEE(_M0L6_2atmpS1705);
    moonbit_decref(_M0L6_2atmpS1705);
    return _result_3472;
  }
}

struct _M0TP26biolab8bio__seq3Seq** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq** _M0L3srcS140,
  int32_t _M0L13allocate__lenS136,
  int32_t _M0L3lenS137,
  int32_t _M0L11src__offsetS138,
  int32_t _M0L11dst__offsetS139
) {
  int32_t _if__result_3473;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS136 >= 0) {
    if (_M0L3lenS137 >= 0) {
      if (_M0L11src__offsetS138 >= 0) {
        if (_M0L11dst__offsetS139 >= 0) {
          int32_t _M0L6_2atmpS1708 = _M0L11src__offsetS138 + _M0L3lenS137;
          int32_t _M0L6_2atmpS1709;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1709
          = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq3SeqE(_M0L3srcS140);
          if (_M0L6_2atmpS1708 <= _M0L6_2atmpS1709) {
            int32_t _M0L6_2atmpS1707 = _M0L11dst__offsetS139 + _M0L3lenS137;
            _if__result_3473 = _M0L6_2atmpS1707 <= _M0L13allocate__lenS136;
          } else {
            _if__result_3473 = 0;
          }
        } else {
          _if__result_3473 = 0;
        }
      } else {
        _if__result_3473 = 0;
      }
    } else {
      _if__result_3473 = 0;
    }
  } else {
    _if__result_3473 = 0;
  }
  if (_if__result_3473) {
    moonbit_incref(_M0L3srcS140);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TP26biolab8bio__seq3Seq**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS136, 0, _M0L3srcS140, _M0L11src__offsetS138, _M0L11dst__offsetS139, _M0L3lenS137);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS141;
    int32_t _M0L6_2atmpS1711;
    moonbit_string_t _M0L6_2atmpS1710;
    struct _M0TP26biolab8bio__seq3Seq** _result_3474;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS141
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS141, (moonbit_string_t)moonbit_string_literal_122.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS141, _M0L13allocate__lenS136);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS141, (moonbit_string_t)moonbit_string_literal_123.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS141, _M0L11src__offsetS138);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS141, (moonbit_string_t)moonbit_string_literal_124.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS141, _M0L11dst__offsetS139);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS141, (moonbit_string_t)moonbit_string_literal_125.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS141, _M0L3lenS137);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS141, (moonbit_string_t)moonbit_string_literal_126.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1711
    = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq3SeqE(_M0L3srcS140);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS141, _M0L6_2atmpS1711);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1710
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS141);
    moonbit_decref(_M0L18_2astring__builderS141);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3474
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq3SeqEE(_M0L6_2atmpS1710);
    moonbit_decref(_M0L6_2atmpS1710);
    return _result_3474;
  }
}

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder* _M0L4selfS115,
  moonbit_string_t _M0L3objS114
) {
  struct _M0TPB6Logger _M0L6_2atmpS1690;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS115);
  _M0L6_2atmpS1690
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS115
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS114, _M0L6_2atmpS1690);
  if (_M0L6_2atmpS1690.$1) {
    moonbit_decref(_M0L6_2atmpS1690.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder* _M0L4selfS117,
  int32_t _M0L3objS116
) {
  struct _M0TPB6Logger _M0L6_2atmpS1691;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS117);
  _M0L6_2atmpS1691
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS117
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L3objS116, _M0L6_2atmpS1691);
  if (_M0L6_2atmpS1691.$1) {
    moonbit_decref(_M0L6_2atmpS1691.$1);
  }
  return 0;
}

moonbit_string_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGsE(
  moonbit_string_t* _M0L3srcS93,
  int32_t _M0L13allocate__lenS91,
  int32_t _M0L11src__offsetS94,
  int32_t _M0L11dst__offsetS92,
  int32_t _M0L9blit__lenS95
) {
  moonbit_string_t* _M0L3dstS90;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS90
  = (moonbit_string_t*)moonbit_make_ref_array(_M0L13allocate__lenS91, (moonbit_string_t)moonbit_string_literal_0.data);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGsE(_M0L3dstS90, _M0L11dst__offsetS92, _M0L3srcS93, _M0L11src__offsetS94, _M0L9blit__lenS95);
  moonbit_decref(_M0L3srcS93);
  return _M0L3dstS90;
}

struct _M0TUsiE** _M0MPB18UninitializedArray23unsafe__make__and__blitGUsiEE(
  struct _M0TUsiE** _M0L3srcS99,
  int32_t _M0L13allocate__lenS97,
  int32_t _M0L11src__offsetS100,
  int32_t _M0L11dst__offsetS98,
  int32_t _M0L9blit__lenS101
) {
  struct _M0TUsiE** _M0L3dstS96;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS96
  = (struct _M0TUsiE**)moonbit_make_ref_array(_M0L13allocate__lenS97, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGUsiEE(_M0L3dstS96, _M0L11dst__offsetS98, _M0L3srcS99, _M0L11src__offsetS100, _M0L9blit__lenS101);
  moonbit_decref(_M0L3srcS99);
  return _M0L3dstS96;
}

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3srcS105,
  int32_t _M0L13allocate__lenS103,
  int32_t _M0L11src__offsetS106,
  int32_t _M0L11dst__offsetS104,
  int32_t _M0L9blit__lenS107
) {
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3dstS102;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS102
  = (struct _M0TP46biolab8bio__seq3cmd4main4Case**)moonbit_make_ref_array(_M0L13allocate__lenS103, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRP46biolab8bio__seq3cmd4main4CaseE(_M0L3dstS102, _M0L11dst__offsetS104, _M0L3srcS105, _M0L11src__offsetS106, _M0L9blit__lenS107);
  moonbit_decref(_M0L3srcS105);
  return _M0L3dstS102;
}

struct _M0TP26biolab8bio__seq3Seq** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq** _M0L3srcS111,
  int32_t _M0L13allocate__lenS109,
  int32_t _M0L11src__offsetS112,
  int32_t _M0L11dst__offsetS110,
  int32_t _M0L9blit__lenS113
) {
  struct _M0TP26biolab8bio__seq3Seq** _M0L3dstS108;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS108
  = (struct _M0TP26biolab8bio__seq3Seq**)moonbit_make_ref_array(_M0L13allocate__lenS109, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq3SeqE(_M0L3dstS108, _M0L11dst__offsetS110, _M0L3srcS111, _M0L11src__offsetS112, _M0L9blit__lenS113);
  moonbit_decref(_M0L3srcS111);
  return _M0L3dstS108;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGsE(
  moonbit_string_t* _M0L3dstS70,
  int32_t _M0L11dst__offsetS71,
  moonbit_string_t* _M0L3srcS72,
  int32_t _M0L11src__offsetS73,
  int32_t _M0L3lenS74
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS72);
  moonbit_incref(_M0L3dstS70);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS70, _M0L11dst__offsetS71, _M0L3srcS72, _M0L11src__offsetS73, _M0L3lenS74);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGUsiEE(
  struct _M0TUsiE** _M0L3dstS75,
  int32_t _M0L11dst__offsetS76,
  struct _M0TUsiE** _M0L3srcS77,
  int32_t _M0L11src__offsetS78,
  int32_t _M0L3lenS79
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS77);
  moonbit_incref(_M0L3dstS75);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS75, _M0L11dst__offsetS76, _M0L3srcS77, _M0L11src__offsetS78, _M0L3lenS79);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3dstS80,
  int32_t _M0L11dst__offsetS81,
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3srcS82,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq** _M0L3dstS85,
  int32_t _M0L11dst__offsetS86,
  struct _M0TP26biolab8bio__seq3Seq** _M0L3srcS87,
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGkE(
  uint16_t* _M0L3dstS25,
  int32_t _M0L11dst__offsetS27,
  uint16_t* _M0L3srcS26,
  int32_t _M0L11src__offsetS28,
  int32_t _M0L3lenS30
) {
  int32_t _if__result_3475;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS25 == _M0L3srcS26) {
    _if__result_3475 = _M0L11dst__offsetS27 < _M0L11src__offsetS28;
  } else {
    _if__result_3475 = 0;
  }
  if (_if__result_3475) {
    int32_t _M0L1iS29 = 0;
    while (1) {
      if (_M0L1iS29 < _M0L3lenS30) {
        int32_t _M0L6_2atmpS1645 = _M0L11dst__offsetS27 + _M0L1iS29;
        int32_t _M0L6_2atmpS1647 = _M0L11src__offsetS28 + _M0L1iS29;
        int32_t _M0L6_2atmpS1646;
        int32_t _M0L6_2atmpS1648;
        if (
          _M0L6_2atmpS1647 < 0
          || _M0L6_2atmpS1647 >= Moonbit_array_length(_M0L3srcS26)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1646 = (int32_t)_M0L3srcS26[_M0L6_2atmpS1647];
        if (
          _M0L6_2atmpS1645 < 0
          || _M0L6_2atmpS1645 >= Moonbit_array_length(_M0L3dstS25)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS25[_M0L6_2atmpS1645] = _M0L6_2atmpS1646;
        _M0L6_2atmpS1648 = _M0L1iS29 + 1;
        _M0L1iS29 = _M0L6_2atmpS1648;
        continue;
      } else {
        moonbit_decref(_M0L3srcS26);
        moonbit_decref(_M0L3dstS25);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1653 = _M0L3lenS30 - 1;
    int32_t _M0L1iS32 = _M0L6_2atmpS1653;
    while (1) {
      if (_M0L1iS32 >= 0) {
        int32_t _M0L6_2atmpS1649 = _M0L11dst__offsetS27 + _M0L1iS32;
        int32_t _M0L6_2atmpS1651 = _M0L11src__offsetS28 + _M0L1iS32;
        int32_t _M0L6_2atmpS1650;
        int32_t _M0L6_2atmpS1652;
        if (
          _M0L6_2atmpS1651 < 0
          || _M0L6_2atmpS1651 >= Moonbit_array_length(_M0L3srcS26)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1650 = (int32_t)_M0L3srcS26[_M0L6_2atmpS1651];
        if (
          _M0L6_2atmpS1649 < 0
          || _M0L6_2atmpS1649 >= Moonbit_array_length(_M0L3dstS25)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS25[_M0L6_2atmpS1649] = _M0L6_2atmpS1650;
        _M0L6_2atmpS1652 = _M0L1iS32 - 1;
        _M0L1iS32 = _M0L6_2atmpS1652;
        continue;
      } else {
        moonbit_decref(_M0L3srcS26);
        moonbit_decref(_M0L3dstS25);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGsEE(
  moonbit_string_t* _M0L3dstS34,
  int32_t _M0L11dst__offsetS36,
  moonbit_string_t* _M0L3srcS35,
  int32_t _M0L11src__offsetS37,
  int32_t _M0L3lenS39
) {
  int32_t _if__result_3478;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS34 == _M0L3srcS35) {
    _if__result_3478 = _M0L11dst__offsetS36 < _M0L11src__offsetS37;
  } else {
    _if__result_3478 = 0;
  }
  if (_if__result_3478) {
    int32_t _M0L1iS38 = 0;
    while (1) {
      if (_M0L1iS38 < _M0L3lenS39) {
        int32_t _M0L6_2atmpS1654 = _M0L11dst__offsetS36 + _M0L1iS38;
        int32_t _M0L6_2atmpS1656 = _M0L11src__offsetS37 + _M0L1iS38;
        moonbit_string_t _M0L6_2atmpS1655;
        moonbit_string_t _M0L6_2aoldS3201;
        int32_t _M0L6_2atmpS1657;
        if (
          _M0L6_2atmpS1656 < 0
          || _M0L6_2atmpS1656 >= Moonbit_array_length(_M0L3srcS35)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1655 = (moonbit_string_t)_M0L3srcS35[_M0L6_2atmpS1656];
        if (
          _M0L6_2atmpS1654 < 0
          || _M0L6_2atmpS1654 >= Moonbit_array_length(_M0L3dstS34)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3201 = (moonbit_string_t)_M0L3dstS34[_M0L6_2atmpS1654];
        moonbit_incref(_M0L6_2atmpS1655);
        moonbit_decref(_M0L6_2aoldS3201);
        _M0L3dstS34[_M0L6_2atmpS1654] = _M0L6_2atmpS1655;
        _M0L6_2atmpS1657 = _M0L1iS38 + 1;
        _M0L1iS38 = _M0L6_2atmpS1657;
        continue;
      } else {
        moonbit_decref(_M0L3srcS35);
        moonbit_decref(_M0L3dstS34);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1662 = _M0L3lenS39 - 1;
    int32_t _M0L1iS41 = _M0L6_2atmpS1662;
    while (1) {
      if (_M0L1iS41 >= 0) {
        int32_t _M0L6_2atmpS1658 = _M0L11dst__offsetS36 + _M0L1iS41;
        int32_t _M0L6_2atmpS1660 = _M0L11src__offsetS37 + _M0L1iS41;
        moonbit_string_t _M0L6_2atmpS1659;
        moonbit_string_t _M0L6_2aoldS3203;
        int32_t _M0L6_2atmpS1661;
        if (
          _M0L6_2atmpS1660 < 0
          || _M0L6_2atmpS1660 >= Moonbit_array_length(_M0L3srcS35)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1659 = (moonbit_string_t)_M0L3srcS35[_M0L6_2atmpS1660];
        if (
          _M0L6_2atmpS1658 < 0
          || _M0L6_2atmpS1658 >= Moonbit_array_length(_M0L3dstS34)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3203 = (moonbit_string_t)_M0L3dstS34[_M0L6_2atmpS1658];
        moonbit_incref(_M0L6_2atmpS1659);
        moonbit_decref(_M0L6_2aoldS3203);
        _M0L3dstS34[_M0L6_2atmpS1658] = _M0L6_2atmpS1659;
        _M0L6_2atmpS1661 = _M0L1iS41 - 1;
        _M0L1iS41 = _M0L6_2atmpS1661;
        continue;
      } else {
        moonbit_decref(_M0L3srcS35);
        moonbit_decref(_M0L3dstS34);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGUsiEEE(
  struct _M0TUsiE** _M0L3dstS43,
  int32_t _M0L11dst__offsetS45,
  struct _M0TUsiE** _M0L3srcS44,
  int32_t _M0L11src__offsetS46,
  int32_t _M0L3lenS48
) {
  int32_t _if__result_3481;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS43 == _M0L3srcS44) {
    _if__result_3481 = _M0L11dst__offsetS45 < _M0L11src__offsetS46;
  } else {
    _if__result_3481 = 0;
  }
  if (_if__result_3481) {
    int32_t _M0L1iS47 = 0;
    while (1) {
      if (_M0L1iS47 < _M0L3lenS48) {
        int32_t _M0L6_2atmpS1663 = _M0L11dst__offsetS45 + _M0L1iS47;
        int32_t _M0L6_2atmpS1665 = _M0L11src__offsetS46 + _M0L1iS47;
        struct _M0TUsiE* _M0L6_2atmpS1664;
        struct _M0TUsiE* _M0L6_2aoldS3205;
        int32_t _M0L6_2atmpS1666;
        if (
          _M0L6_2atmpS1665 < 0
          || _M0L6_2atmpS1665 >= Moonbit_array_length(_M0L3srcS44)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1664 = (struct _M0TUsiE*)_M0L3srcS44[_M0L6_2atmpS1665];
        if (
          _M0L6_2atmpS1663 < 0
          || _M0L6_2atmpS1663 >= Moonbit_array_length(_M0L3dstS43)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3205 = (struct _M0TUsiE*)_M0L3dstS43[_M0L6_2atmpS1663];
        if (_M0L6_2atmpS1664) {
          moonbit_incref(_M0L6_2atmpS1664);
        }
        if (_M0L6_2aoldS3205) {
          moonbit_decref(_M0L6_2aoldS3205);
        }
        _M0L3dstS43[_M0L6_2atmpS1663] = _M0L6_2atmpS1664;
        _M0L6_2atmpS1666 = _M0L1iS47 + 1;
        _M0L1iS47 = _M0L6_2atmpS1666;
        continue;
      } else {
        moonbit_decref(_M0L3srcS44);
        moonbit_decref(_M0L3dstS43);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1671 = _M0L3lenS48 - 1;
    int32_t _M0L1iS50 = _M0L6_2atmpS1671;
    while (1) {
      if (_M0L1iS50 >= 0) {
        int32_t _M0L6_2atmpS1667 = _M0L11dst__offsetS45 + _M0L1iS50;
        int32_t _M0L6_2atmpS1669 = _M0L11src__offsetS46 + _M0L1iS50;
        struct _M0TUsiE* _M0L6_2atmpS1668;
        struct _M0TUsiE* _M0L6_2aoldS3207;
        int32_t _M0L6_2atmpS1670;
        if (
          _M0L6_2atmpS1669 < 0
          || _M0L6_2atmpS1669 >= Moonbit_array_length(_M0L3srcS44)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1668 = (struct _M0TUsiE*)_M0L3srcS44[_M0L6_2atmpS1669];
        if (
          _M0L6_2atmpS1667 < 0
          || _M0L6_2atmpS1667 >= Moonbit_array_length(_M0L3dstS43)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3207 = (struct _M0TUsiE*)_M0L3dstS43[_M0L6_2atmpS1667];
        if (_M0L6_2atmpS1668) {
          moonbit_incref(_M0L6_2atmpS1668);
        }
        if (_M0L6_2aoldS3207) {
          moonbit_decref(_M0L6_2aoldS3207);
        }
        _M0L3dstS43[_M0L6_2atmpS1667] = _M0L6_2atmpS1668;
        _M0L6_2atmpS1670 = _M0L1iS50 - 1;
        _M0L1iS50 = _M0L6_2atmpS1670;
        continue;
      } else {
        moonbit_decref(_M0L3srcS44);
        moonbit_decref(_M0L3dstS43);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP46biolab8bio__seq3cmd4main4CaseEE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3dstS52,
  int32_t _M0L11dst__offsetS54,
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L3srcS53,
  int32_t _M0L11src__offsetS55,
  int32_t _M0L3lenS57
) {
  int32_t _if__result_3484;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS52 == _M0L3srcS53) {
    _if__result_3484 = _M0L11dst__offsetS54 < _M0L11src__offsetS55;
  } else {
    _if__result_3484 = 0;
  }
  if (_if__result_3484) {
    int32_t _M0L1iS56 = 0;
    while (1) {
      if (_M0L1iS56 < _M0L3lenS57) {
        int32_t _M0L6_2atmpS1672 = _M0L11dst__offsetS54 + _M0L1iS56;
        int32_t _M0L6_2atmpS1674 = _M0L11src__offsetS55 + _M0L1iS56;
        struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L6_2atmpS1673;
        struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L6_2aoldS3209;
        int32_t _M0L6_2atmpS1675;
        if (
          _M0L6_2atmpS1674 < 0
          || _M0L6_2atmpS1674 >= Moonbit_array_length(_M0L3srcS53)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1673
        = (struct _M0TP46biolab8bio__seq3cmd4main4Case*)_M0L3srcS53[
            _M0L6_2atmpS1674
          ];
        if (
          _M0L6_2atmpS1672 < 0
          || _M0L6_2atmpS1672 >= Moonbit_array_length(_M0L3dstS52)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3209
        = (struct _M0TP46biolab8bio__seq3cmd4main4Case*)_M0L3dstS52[
            _M0L6_2atmpS1672
          ];
        if (_M0L6_2atmpS1673) {
          moonbit_incref(_M0L6_2atmpS1673);
        }
        if (_M0L6_2aoldS3209) {
          moonbit_decref(_M0L6_2aoldS3209);
        }
        _M0L3dstS52[_M0L6_2atmpS1672] = _M0L6_2atmpS1673;
        _M0L6_2atmpS1675 = _M0L1iS56 + 1;
        _M0L1iS56 = _M0L6_2atmpS1675;
        continue;
      } else {
        moonbit_decref(_M0L3srcS53);
        moonbit_decref(_M0L3dstS52);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1680 = _M0L3lenS57 - 1;
    int32_t _M0L1iS59 = _M0L6_2atmpS1680;
    while (1) {
      if (_M0L1iS59 >= 0) {
        int32_t _M0L6_2atmpS1676 = _M0L11dst__offsetS54 + _M0L1iS59;
        int32_t _M0L6_2atmpS1678 = _M0L11src__offsetS55 + _M0L1iS59;
        struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L6_2atmpS1677;
        struct _M0TP46biolab8bio__seq3cmd4main4Case* _M0L6_2aoldS3211;
        int32_t _M0L6_2atmpS1679;
        if (
          _M0L6_2atmpS1678 < 0
          || _M0L6_2atmpS1678 >= Moonbit_array_length(_M0L3srcS53)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1677
        = (struct _M0TP46biolab8bio__seq3cmd4main4Case*)_M0L3srcS53[
            _M0L6_2atmpS1678
          ];
        if (
          _M0L6_2atmpS1676 < 0
          || _M0L6_2atmpS1676 >= Moonbit_array_length(_M0L3dstS52)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3211
        = (struct _M0TP46biolab8bio__seq3cmd4main4Case*)_M0L3dstS52[
            _M0L6_2atmpS1676
          ];
        if (_M0L6_2atmpS1677) {
          moonbit_incref(_M0L6_2atmpS1677);
        }
        if (_M0L6_2aoldS3211) {
          moonbit_decref(_M0L6_2aoldS3211);
        }
        _M0L3dstS52[_M0L6_2atmpS1676] = _M0L6_2atmpS1677;
        _M0L6_2atmpS1679 = _M0L1iS59 - 1;
        _M0L1iS59 = _M0L6_2atmpS1679;
        continue;
      } else {
        moonbit_decref(_M0L3srcS53);
        moonbit_decref(_M0L3dstS52);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq3SeqEE(
  struct _M0TP26biolab8bio__seq3Seq** _M0L3dstS61,
  int32_t _M0L11dst__offsetS63,
  struct _M0TP26biolab8bio__seq3Seq** _M0L3srcS62,
  int32_t _M0L11src__offsetS64,
  int32_t _M0L3lenS66
) {
  int32_t _if__result_3487;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS61 == _M0L3srcS62) {
    _if__result_3487 = _M0L11dst__offsetS63 < _M0L11src__offsetS64;
  } else {
    _if__result_3487 = 0;
  }
  if (_if__result_3487) {
    int32_t _M0L1iS65 = 0;
    while (1) {
      if (_M0L1iS65 < _M0L3lenS66) {
        int32_t _M0L6_2atmpS1681 = _M0L11dst__offsetS63 + _M0L1iS65;
        int32_t _M0L6_2atmpS1683 = _M0L11src__offsetS64 + _M0L1iS65;
        struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS1682;
        struct _M0TP26biolab8bio__seq3Seq* _M0L6_2aoldS3213;
        int32_t _M0L6_2atmpS1684;
        if (
          _M0L6_2atmpS1683 < 0
          || _M0L6_2atmpS1683 >= Moonbit_array_length(_M0L3srcS62)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1682
        = (struct _M0TP26biolab8bio__seq3Seq*)_M0L3srcS62[_M0L6_2atmpS1683];
        if (
          _M0L6_2atmpS1681 < 0
          || _M0L6_2atmpS1681 >= Moonbit_array_length(_M0L3dstS61)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3213
        = (struct _M0TP26biolab8bio__seq3Seq*)_M0L3dstS61[_M0L6_2atmpS1681];
        if (_M0L6_2atmpS1682) {
          moonbit_incref(_M0L6_2atmpS1682);
        }
        if (_M0L6_2aoldS3213) {
          moonbit_decref(_M0L6_2aoldS3213);
        }
        _M0L3dstS61[_M0L6_2atmpS1681] = _M0L6_2atmpS1682;
        _M0L6_2atmpS1684 = _M0L1iS65 + 1;
        _M0L1iS65 = _M0L6_2atmpS1684;
        continue;
      } else {
        moonbit_decref(_M0L3srcS62);
        moonbit_decref(_M0L3dstS61);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1689 = _M0L3lenS66 - 1;
    int32_t _M0L1iS68 = _M0L6_2atmpS1689;
    while (1) {
      if (_M0L1iS68 >= 0) {
        int32_t _M0L6_2atmpS1685 = _M0L11dst__offsetS63 + _M0L1iS68;
        int32_t _M0L6_2atmpS1687 = _M0L11src__offsetS64 + _M0L1iS68;
        struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS1686;
        struct _M0TP26biolab8bio__seq3Seq* _M0L6_2aoldS3215;
        int32_t _M0L6_2atmpS1688;
        if (
          _M0L6_2atmpS1687 < 0
          || _M0L6_2atmpS1687 >= Moonbit_array_length(_M0L3srcS62)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1686
        = (struct _M0TP26biolab8bio__seq3Seq*)_M0L3srcS62[_M0L6_2atmpS1687];
        if (
          _M0L6_2atmpS1685 < 0
          || _M0L6_2atmpS1685 >= Moonbit_array_length(_M0L3dstS61)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3215
        = (struct _M0TP26biolab8bio__seq3Seq*)_M0L3dstS61[_M0L6_2atmpS1685];
        if (_M0L6_2atmpS1686) {
          moonbit_incref(_M0L6_2atmpS1686);
        }
        if (_M0L6_2aoldS3215) {
          moonbit_decref(_M0L6_2aoldS3215);
        }
        _M0L3dstS61[_M0L6_2atmpS1685] = _M0L6_2atmpS1686;
        _M0L6_2atmpS1688 = _M0L1iS68 - 1;
        _M0L1iS68 = _M0L6_2atmpS1688;
        continue;
      } else {
        moonbit_decref(_M0L3srcS62);
        moonbit_decref(_M0L3dstS61);
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

int32_t _M0MPB18UninitializedArray6lengthGRP46biolab8bio__seq3cmd4main4CaseE(
  struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0L4selfS23
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS23);
}

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq3SeqE(
  struct _M0TP26biolab8bio__seq3Seq** _M0L4selfS24
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS24);
}

uint32_t _M0FPB13consume4__acc(uint32_t _M0L3accS19, uint32_t _M0L5inputS20) {
  uint32_t _M0L6_2atmpS1644;
  uint32_t _M0L6_2atmpS1643;
  uint32_t _M0L6_2atmpS1642;
  #line 465 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1644 = _M0L5inputS20 * 3266489917u;
  _M0L6_2atmpS1643 = _M0L3accS19 + _M0L6_2atmpS1644;
  #line 466 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1642 = _M0FPB4rotl(_M0L6_2atmpS1643, 17);
  return _M0L6_2atmpS1642 * 668265263u;
}

uint32_t _M0FPB4rotl(uint32_t _M0L1xS17, int32_t _M0L1rS18) {
  uint32_t _M0L6_2atmpS1639;
  int32_t _M0L6_2atmpS1641;
  uint32_t _M0L6_2atmpS1640;
  #line 475 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1639 = _M0L1xS17 << (_M0L1rS18 & 31);
  _M0L6_2atmpS1641 = 32 - _M0L1rS18;
  _M0L6_2atmpS1640 = _M0L1xS17 >> (_M0L6_2atmpS1641 & 31);
  return _M0L6_2atmpS1639 | _M0L6_2atmpS1640;
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
  _M0L10_2ax__5722S16.$0->$method_0(_M0L10_2ax__5722S16.$1, (moonbit_string_t)moonbit_string_literal_127.data);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB6Logger13write__objectGsE(_M0L10_2ax__5722S16, _M0L15_2a_2aarg__5723S15);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S16.$0->$method_0(_M0L10_2ax__5722S16.$1, (moonbit_string_t)moonbit_string_literal_128.data);
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

int32_t _M0FPC15abort5abortGiE(moonbit_string_t _M0L3msgS3) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS3);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

moonbit_string_t _M0FPC15abort5abortGsE(moonbit_string_t _M0L3msgS4) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS4);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TPC16string10StringView _M0FPC15abort5abortGRPC16string10StringViewE(
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

struct _M0TP46biolab8bio__seq3cmd4main4Case** _M0FPC15abort5abortGRPB18UninitializedArrayGRP46biolab8bio__seq3cmd4main4CaseEE(
  moonbit_string_t _M0L3msgS8
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS8);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int64_t _M0FPC15abort5abortGOiE(moonbit_string_t _M0L3msgS9) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS9);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TP26biolab8bio__seq3Seq** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq3SeqEE(
  moonbit_string_t _M0L3msgS10
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS10);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t _M0FP017____moonbit__initN3regS1004(
  struct _M0TPB3MapGicE* _M0L1mS1002,
  int32_t _M0L2c0S1005,
  int32_t _M0L2c1S1006,
  int32_t _M0L2c2S1007,
  int32_t _M0L2aaS1008
) {
  int32_t _M0L6_2atmpS1451;
  int32_t _M0L6_2atmpS1452;
  int32_t _M0L6_2atmpS1450;
  int32_t _M0L6_2atmpS1449;
  #line 17 "/home/developer/Documents/1/codon_table.mbt"
  _M0L6_2atmpS1451 = _M0L2c0S1005 * 65536;
  _M0L6_2atmpS1452 = _M0L2c1S1006 * 256;
  _M0L6_2atmpS1450 = _M0L6_2atmpS1451 + _M0L6_2atmpS1452;
  _M0L6_2atmpS1449 = _M0L6_2atmpS1450 + _M0L2c2S1007;
  #line 18 "/home/developer/Documents/1/codon_table.mbt"
  _M0MPB3Map3setGicE(_M0L1mS1002, _M0L6_2atmpS1449, _M0L2aaS1008);
  return 0;
}

moonbit_string_t _M0FP15Error10to__string(void* _M0L4_2aeS1215) {
  switch (Moonbit_object_tag(_M0L4_2aeS1215)) {
    case 2: {
      return (moonbit_string_t)moonbit_string_literal_129.data;
      break;
    }
    
    case 3: {
      return (moonbit_string_t)moonbit_string_literal_130.data;
      break;
    }
    
    case 5: {
      return (moonbit_string_t)moonbit_string_literal_131.data;
      break;
    }
    
    case 1: {
      return (moonbit_string_t)moonbit_string_literal_132.data;
      break;
    }
    
    case 0: {
      return _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(_M0L4_2aeS1215);
      break;
    }
    default: {
      return (moonbit_string_t)moonbit_string_literal_133.data;
      break;
    }
  }
}

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1239,
  struct _M0TPB4Show _M0L8_2aparamS1238
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1237 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1239;
  _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(_M0L7_2aselfS1237, _M0L8_2aparamS1238);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1236,
  struct _M0TPB4Show _M0L8_2aparamS1235
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1234 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1236;
  _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(_M0L7_2aselfS1234, _M0L8_2aparamS1235);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1233,
  int32_t _M0L8_2aparamS1232
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1231 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1233;
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS1231, _M0L8_2aparamS1232);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1230,
  struct _M0TPC16string10StringView _M0L8_2aparamS1229
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1228 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1230;
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L7_2aselfS1228, _M0L8_2aparamS1229);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1227,
  moonbit_string_t _M0L8_2aparamS1224,
  int32_t _M0L8_2aparamS1225,
  int32_t _M0L8_2aparamS1226
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1223 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1227;
  _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L7_2aselfS1223, _M0L8_2aparamS1224, _M0L8_2aparamS1225, _M0L8_2aparamS1226);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1222,
  moonbit_string_t _M0L8_2aparamS1221
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1220 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1222;
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L7_2aselfS1220, _M0L8_2aparamS1221);
  return 0;
}

void moonbit_init() {
  moonbit_layout_table = moonbit_layout_table_data;
  uint16_t* _M0L1tS813 = (uint16_t*)moonbit_make_string(128, 0);
  int32_t _M0L6_2atmpS1245 = 65;
  int32_t _M0L6_2atmpS1247 = 84;
  int32_t _M0L6_2atmpS1246 = (uint16_t)_M0L6_2atmpS1247;
  int32_t _M0L6_2atmpS1248;
  int32_t _M0L6_2atmpS1250;
  int32_t _M0L6_2atmpS1249;
  int32_t _M0L6_2atmpS1251;
  int32_t _M0L6_2atmpS1253;
  int32_t _M0L6_2atmpS1252;
  int32_t _M0L6_2atmpS1254;
  int32_t _M0L6_2atmpS1256;
  int32_t _M0L6_2atmpS1255;
  int32_t _M0L6_2atmpS1257;
  int32_t _M0L6_2atmpS1259;
  int32_t _M0L6_2atmpS1258;
  int32_t _M0L6_2atmpS1260;
  int32_t _M0L6_2atmpS1262;
  int32_t _M0L6_2atmpS1261;
  int32_t _M0L6_2atmpS1263;
  int32_t _M0L6_2atmpS1265;
  int32_t _M0L6_2atmpS1264;
  int32_t _M0L6_2atmpS1266;
  int32_t _M0L6_2atmpS1268;
  int32_t _M0L6_2atmpS1267;
  int32_t _M0L6_2atmpS1269;
  int32_t _M0L6_2atmpS1271;
  int32_t _M0L6_2atmpS1270;
  int32_t _M0L6_2atmpS1272;
  int32_t _M0L6_2atmpS1274;
  int32_t _M0L6_2atmpS1273;
  int32_t _M0L6_2atmpS1275;
  int32_t _M0L6_2atmpS1277;
  int32_t _M0L6_2atmpS1276;
  int32_t _M0L6_2atmpS1278;
  int32_t _M0L6_2atmpS1280;
  int32_t _M0L6_2atmpS1279;
  int32_t _M0L6_2atmpS1281;
  int32_t _M0L6_2atmpS1283;
  int32_t _M0L6_2atmpS1282;
  int32_t _M0L6_2atmpS1284;
  int32_t _M0L6_2atmpS1286;
  int32_t _M0L6_2atmpS1285;
  int32_t _M0L6_2atmpS1287;
  int32_t _M0L6_2atmpS1289;
  int32_t _M0L6_2atmpS1288;
  int32_t _M0L6_2atmpS1290;
  int32_t _M0L6_2atmpS1292;
  int32_t _M0L6_2atmpS1291;
  int32_t _M0L6_2atmpS1293;
  int32_t _M0L6_2atmpS1295;
  int32_t _M0L6_2atmpS1294;
  int32_t _M0L6_2atmpS1296;
  int32_t _M0L6_2atmpS1298;
  int32_t _M0L6_2atmpS1297;
  int32_t _M0L6_2atmpS1299;
  int32_t _M0L6_2atmpS1301;
  int32_t _M0L6_2atmpS1300;
  int32_t _M0L6_2atmpS1302;
  int32_t _M0L6_2atmpS1304;
  int32_t _M0L6_2atmpS1303;
  int32_t _M0L6_2atmpS1305;
  int32_t _M0L6_2atmpS1307;
  int32_t _M0L6_2atmpS1306;
  int32_t _M0L6_2atmpS1308;
  int32_t _M0L6_2atmpS1310;
  int32_t _M0L6_2atmpS1309;
  int32_t _M0L6_2atmpS1311;
  int32_t _M0L6_2atmpS1313;
  int32_t _M0L6_2atmpS1312;
  int32_t _M0L6_2atmpS1314;
  int32_t _M0L6_2atmpS1316;
  int32_t _M0L6_2atmpS1315;
  int32_t _M0L6_2atmpS1317;
  int32_t _M0L6_2atmpS1319;
  int32_t _M0L6_2atmpS1318;
  int32_t _M0L6_2atmpS1320;
  int32_t _M0L6_2atmpS1322;
  int32_t _M0L6_2atmpS1321;
  int32_t _M0L6_2atmpS1323;
  int32_t _M0L6_2atmpS1325;
  int32_t _M0L6_2atmpS1324;
  int32_t _M0L6_2atmpS1326;
  int32_t _M0L6_2atmpS1328;
  int32_t _M0L6_2atmpS1327;
  int32_t _M0L6_2atmpS1329;
  int32_t _M0L6_2atmpS1331;
  int32_t _M0L6_2atmpS1330;
  int32_t _M0L6_2atmpS1332;
  int32_t _M0L6_2atmpS1334;
  int32_t _M0L6_2atmpS1333;
  int32_t _M0L6_2atmpS1335;
  int32_t _M0L6_2atmpS1337;
  int32_t _M0L6_2atmpS1336;
  int32_t _M0L6_2atmpS1338;
  int32_t _M0L6_2atmpS1340;
  int32_t _M0L6_2atmpS1339;
  int32_t _M0L6_2atmpS1341;
  int32_t _M0L6_2atmpS1343;
  int32_t _M0L6_2atmpS1342;
  int32_t _M0L6_2atmpS1344;
  int32_t _M0L6_2atmpS1346;
  int32_t _M0L6_2atmpS1345;
  uint16_t* _M0L1tS936;
  int32_t _M0L6_2atmpS1347;
  int32_t _M0L6_2atmpS1349;
  int32_t _M0L6_2atmpS1348;
  int32_t _M0L6_2atmpS1350;
  int32_t _M0L6_2atmpS1352;
  int32_t _M0L6_2atmpS1351;
  int32_t _M0L6_2atmpS1353;
  int32_t _M0L6_2atmpS1355;
  int32_t _M0L6_2atmpS1354;
  int32_t _M0L6_2atmpS1356;
  int32_t _M0L6_2atmpS1358;
  int32_t _M0L6_2atmpS1357;
  int32_t _M0L6_2atmpS1359;
  int32_t _M0L6_2atmpS1361;
  int32_t _M0L6_2atmpS1360;
  int32_t _M0L6_2atmpS1362;
  int32_t _M0L6_2atmpS1364;
  int32_t _M0L6_2atmpS1363;
  int32_t _M0L6_2atmpS1365;
  int32_t _M0L6_2atmpS1367;
  int32_t _M0L6_2atmpS1366;
  int32_t _M0L6_2atmpS1368;
  int32_t _M0L6_2atmpS1370;
  int32_t _M0L6_2atmpS1369;
  int32_t _M0L6_2atmpS1371;
  int32_t _M0L6_2atmpS1373;
  int32_t _M0L6_2atmpS1372;
  int32_t _M0L6_2atmpS1374;
  int32_t _M0L6_2atmpS1376;
  int32_t _M0L6_2atmpS1375;
  int32_t _M0L6_2atmpS1377;
  int32_t _M0L6_2atmpS1379;
  int32_t _M0L6_2atmpS1378;
  int32_t _M0L6_2atmpS1380;
  int32_t _M0L6_2atmpS1382;
  int32_t _M0L6_2atmpS1381;
  int32_t _M0L6_2atmpS1383;
  int32_t _M0L6_2atmpS1385;
  int32_t _M0L6_2atmpS1384;
  int32_t _M0L6_2atmpS1386;
  int32_t _M0L6_2atmpS1388;
  int32_t _M0L6_2atmpS1387;
  int32_t _M0L6_2atmpS1389;
  int32_t _M0L6_2atmpS1391;
  int32_t _M0L6_2atmpS1390;
  int32_t _M0L6_2atmpS1392;
  int32_t _M0L6_2atmpS1394;
  int32_t _M0L6_2atmpS1393;
  int32_t _M0L6_2atmpS1395;
  int32_t _M0L6_2atmpS1397;
  int32_t _M0L6_2atmpS1396;
  int32_t _M0L6_2atmpS1398;
  int32_t _M0L6_2atmpS1400;
  int32_t _M0L6_2atmpS1399;
  int32_t _M0L6_2atmpS1401;
  int32_t _M0L6_2atmpS1403;
  int32_t _M0L6_2atmpS1402;
  int32_t _M0L6_2atmpS1404;
  int32_t _M0L6_2atmpS1406;
  int32_t _M0L6_2atmpS1405;
  int32_t _M0L6_2atmpS1407;
  int32_t _M0L6_2atmpS1409;
  int32_t _M0L6_2atmpS1408;
  int32_t _M0L6_2atmpS1410;
  int32_t _M0L6_2atmpS1412;
  int32_t _M0L6_2atmpS1411;
  int32_t _M0L6_2atmpS1413;
  int32_t _M0L6_2atmpS1415;
  int32_t _M0L6_2atmpS1414;
  int32_t _M0L6_2atmpS1416;
  int32_t _M0L6_2atmpS1418;
  int32_t _M0L6_2atmpS1417;
  int32_t _M0L6_2atmpS1419;
  int32_t _M0L6_2atmpS1421;
  int32_t _M0L6_2atmpS1420;
  int32_t _M0L6_2atmpS1422;
  int32_t _M0L6_2atmpS1424;
  int32_t _M0L6_2atmpS1423;
  int32_t _M0L6_2atmpS1425;
  int32_t _M0L6_2atmpS1427;
  int32_t _M0L6_2atmpS1426;
  int32_t _M0L6_2atmpS1428;
  int32_t _M0L6_2atmpS1430;
  int32_t _M0L6_2atmpS1429;
  int32_t _M0L6_2atmpS1431;
  int32_t _M0L6_2atmpS1433;
  int32_t _M0L6_2atmpS1432;
  int32_t _M0L6_2atmpS1434;
  int32_t _M0L6_2atmpS1436;
  int32_t _M0L6_2atmpS1435;
  int32_t _M0L6_2atmpS1437;
  int32_t _M0L6_2atmpS1439;
  int32_t _M0L6_2atmpS1438;
  int32_t _M0L6_2atmpS1440;
  int32_t _M0L6_2atmpS1442;
  int32_t _M0L6_2atmpS1441;
  int32_t _M0L6_2atmpS1443;
  int32_t _M0L6_2atmpS1445;
  int32_t _M0L6_2atmpS1444;
  int32_t _M0L6_2atmpS1446;
  int32_t _M0L6_2atmpS1448;
  int32_t _M0L6_2atmpS1447;
  struct _M0TUicE** _M0L7_2abindS1003;
  struct _M0TUicE** _M0L6_2atmpS1637;
  struct _M0TPB9ArrayViewGUicEE _M0L6_2atmpS1636;
  struct _M0TPB3MapGicE* _M0L1mS1002;
  struct _M0TPB3MapGicE* _M0L3regS1004;
  int32_t _M0L6_2atmpS1453;
  int32_t _M0L6_2atmpS1454;
  int32_t _M0L6_2atmpS1455;
  int32_t _M0L6_2atmpS1456;
  int32_t _M0L6_2atmpS1457;
  int32_t _M0L6_2atmpS1458;
  int32_t _M0L6_2atmpS1459;
  int32_t _M0L6_2atmpS1460;
  int32_t _M0L6_2atmpS1461;
  int32_t _M0L6_2atmpS1462;
  int32_t _M0L6_2atmpS1463;
  int32_t _M0L6_2atmpS1464;
  int32_t _M0L6_2atmpS1465;
  int32_t _M0L6_2atmpS1466;
  int32_t _M0L6_2atmpS1467;
  int32_t _M0L6_2atmpS1468;
  int32_t _M0L6_2atmpS1469;
  int32_t _M0L6_2atmpS1470;
  int32_t _M0L6_2atmpS1471;
  int32_t _M0L6_2atmpS1472;
  int32_t _M0L6_2atmpS1473;
  int32_t _M0L6_2atmpS1474;
  int32_t _M0L6_2atmpS1475;
  int32_t _M0L6_2atmpS1476;
  int32_t _M0L6_2atmpS1477;
  int32_t _M0L6_2atmpS1478;
  int32_t _M0L6_2atmpS1479;
  int32_t _M0L6_2atmpS1480;
  int32_t _M0L6_2atmpS1481;
  int32_t _M0L6_2atmpS1482;
  int32_t _M0L6_2atmpS1483;
  int32_t _M0L6_2atmpS1484;
  int32_t _M0L6_2atmpS1485;
  int32_t _M0L6_2atmpS1486;
  int32_t _M0L6_2atmpS1487;
  int32_t _M0L6_2atmpS1488;
  int32_t _M0L6_2atmpS1489;
  int32_t _M0L6_2atmpS1490;
  int32_t _M0L6_2atmpS1491;
  int32_t _M0L6_2atmpS1492;
  int32_t _M0L6_2atmpS1493;
  int32_t _M0L6_2atmpS1494;
  int32_t _M0L6_2atmpS1495;
  int32_t _M0L6_2atmpS1496;
  int32_t _M0L6_2atmpS1497;
  int32_t _M0L6_2atmpS1498;
  int32_t _M0L6_2atmpS1499;
  int32_t _M0L6_2atmpS1500;
  int32_t _M0L6_2atmpS1501;
  int32_t _M0L6_2atmpS1502;
  int32_t _M0L6_2atmpS1503;
  int32_t _M0L6_2atmpS1504;
  int32_t _M0L6_2atmpS1505;
  int32_t _M0L6_2atmpS1506;
  int32_t _M0L6_2atmpS1507;
  int32_t _M0L6_2atmpS1508;
  int32_t _M0L6_2atmpS1509;
  int32_t _M0L6_2atmpS1510;
  int32_t _M0L6_2atmpS1511;
  int32_t _M0L6_2atmpS1512;
  int32_t _M0L6_2atmpS1513;
  int32_t _M0L6_2atmpS1514;
  int32_t _M0L6_2atmpS1515;
  int32_t _M0L6_2atmpS1516;
  int32_t _M0L6_2atmpS1517;
  int32_t _M0L6_2atmpS1518;
  int32_t _M0L6_2atmpS1519;
  int32_t _M0L6_2atmpS1520;
  int32_t _M0L6_2atmpS1521;
  int32_t _M0L6_2atmpS1522;
  int32_t _M0L6_2atmpS1523;
  int32_t _M0L6_2atmpS1524;
  int32_t _M0L6_2atmpS1525;
  int32_t _M0L6_2atmpS1526;
  int32_t _M0L6_2atmpS1527;
  int32_t _M0L6_2atmpS1528;
  int32_t _M0L6_2atmpS1529;
  int32_t _M0L6_2atmpS1530;
  int32_t _M0L6_2atmpS1531;
  int32_t _M0L6_2atmpS1532;
  int32_t _M0L6_2atmpS1533;
  int32_t _M0L6_2atmpS1534;
  int32_t _M0L6_2atmpS1535;
  int32_t _M0L6_2atmpS1536;
  int32_t _M0L6_2atmpS1537;
  int32_t _M0L6_2atmpS1538;
  int32_t _M0L6_2atmpS1539;
  int32_t _M0L6_2atmpS1540;
  int32_t _M0L6_2atmpS1541;
  int32_t _M0L6_2atmpS1542;
  int32_t _M0L6_2atmpS1543;
  int32_t _M0L6_2atmpS1544;
  int32_t _M0L6_2atmpS1545;
  int32_t _M0L6_2atmpS1546;
  int32_t _M0L6_2atmpS1547;
  int32_t _M0L6_2atmpS1548;
  int32_t _M0L6_2atmpS1549;
  int32_t _M0L6_2atmpS1550;
  int32_t _M0L6_2atmpS1551;
  int32_t _M0L6_2atmpS1552;
  int32_t _M0L6_2atmpS1553;
  int32_t _M0L6_2atmpS1554;
  int32_t _M0L6_2atmpS1555;
  int32_t _M0L6_2atmpS1556;
  int32_t _M0L6_2atmpS1557;
  int32_t _M0L6_2atmpS1558;
  int32_t _M0L6_2atmpS1559;
  int32_t _M0L6_2atmpS1560;
  int32_t _M0L6_2atmpS1561;
  int32_t _M0L6_2atmpS1562;
  int32_t _M0L6_2atmpS1563;
  int32_t _M0L6_2atmpS1564;
  int32_t _M0L6_2atmpS1565;
  int32_t _M0L6_2atmpS1566;
  int32_t _M0L6_2atmpS1567;
  int32_t _M0L6_2atmpS1568;
  int32_t _M0L6_2atmpS1569;
  int32_t _M0L6_2atmpS1570;
  int32_t _M0L6_2atmpS1571;
  int32_t _M0L6_2atmpS1572;
  int32_t _M0L6_2atmpS1573;
  int32_t _M0L6_2atmpS1574;
  int32_t _M0L6_2atmpS1575;
  int32_t _M0L6_2atmpS1576;
  int32_t _M0L6_2atmpS1577;
  int32_t _M0L6_2atmpS1578;
  int32_t _M0L6_2atmpS1579;
  int32_t _M0L6_2atmpS1580;
  int32_t _M0L6_2atmpS1581;
  int32_t _M0L6_2atmpS1582;
  int32_t _M0L6_2atmpS1583;
  int32_t _M0L6_2atmpS1584;
  int32_t _M0L6_2atmpS1585;
  int32_t _M0L6_2atmpS1586;
  int32_t _M0L6_2atmpS1587;
  int32_t _M0L6_2atmpS1588;
  int32_t _M0L6_2atmpS1589;
  int32_t _M0L6_2atmpS1590;
  int32_t _M0L6_2atmpS1591;
  int32_t _M0L6_2atmpS1592;
  int32_t _M0L6_2atmpS1593;
  int32_t _M0L6_2atmpS1594;
  int32_t _M0L6_2atmpS1595;
  int32_t _M0L6_2atmpS1596;
  int32_t _M0L6_2atmpS1597;
  int32_t _M0L6_2atmpS1598;
  int32_t _M0L6_2atmpS1599;
  int32_t _M0L6_2atmpS1600;
  int32_t _M0L6_2atmpS1601;
  int32_t _M0L6_2atmpS1602;
  int32_t _M0L6_2atmpS1603;
  int32_t _M0L6_2atmpS1604;
  int32_t _M0L6_2atmpS1605;
  int32_t _M0L6_2atmpS1606;
  int32_t _M0L6_2atmpS1607;
  int32_t _M0L6_2atmpS1608;
  int32_t _M0L6_2atmpS1609;
  int32_t _M0L6_2atmpS1610;
  int32_t _M0L6_2atmpS1611;
  int32_t _M0L6_2atmpS1612;
  int32_t _M0L6_2atmpS1613;
  int32_t _M0L6_2atmpS1614;
  int32_t _M0L6_2atmpS1615;
  int32_t _M0L6_2atmpS1616;
  int32_t _M0L6_2atmpS1617;
  int32_t _M0L6_2atmpS1618;
  int32_t _M0L6_2atmpS1619;
  int32_t _M0L6_2atmpS1620;
  int32_t _M0L6_2atmpS1621;
  int32_t _M0L6_2atmpS1622;
  int32_t _M0L6_2atmpS1623;
  int32_t _M0L6_2atmpS1624;
  int32_t _M0L6_2atmpS1625;
  int32_t _M0L6_2atmpS1626;
  int32_t _M0L6_2atmpS1627;
  int32_t _M0L6_2atmpS1628;
  int32_t _M0L6_2atmpS1629;
  int32_t _M0L6_2atmpS1630;
  int32_t _M0L6_2atmpS1631;
  int32_t _M0L6_2atmpS1632;
  int32_t _M0L6_2atmpS1633;
  int32_t _M0L6_2atmpS1634;
  int32_t _M0L6_2atmpS1635;
  uint8_t* _M0L1tS967;
  struct _M0TPB4IterGcE* _M0L5_2aitS968;
  if (
    _M0L6_2atmpS1245 < 0
    || _M0L6_2atmpS1245 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 27 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1245] = _M0L6_2atmpS1246;
  _M0L6_2atmpS1248 = 67;
  _M0L6_2atmpS1250 = 71;
  _M0L6_2atmpS1249 = (uint16_t)_M0L6_2atmpS1250;
  if (
    _M0L6_2atmpS1248 < 0
    || _M0L6_2atmpS1248 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 28 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1248] = _M0L6_2atmpS1249;
  _M0L6_2atmpS1251 = 71;
  _M0L6_2atmpS1253 = 67;
  _M0L6_2atmpS1252 = (uint16_t)_M0L6_2atmpS1253;
  if (
    _M0L6_2atmpS1251 < 0
    || _M0L6_2atmpS1251 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 29 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1251] = _M0L6_2atmpS1252;
  _M0L6_2atmpS1254 = 84;
  _M0L6_2atmpS1256 = 65;
  _M0L6_2atmpS1255 = (uint16_t)_M0L6_2atmpS1256;
  if (
    _M0L6_2atmpS1254 < 0
    || _M0L6_2atmpS1254 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 30 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1254] = _M0L6_2atmpS1255;
  _M0L6_2atmpS1257 = 85;
  _M0L6_2atmpS1259 = 65;
  _M0L6_2atmpS1258 = (uint16_t)_M0L6_2atmpS1259;
  if (
    _M0L6_2atmpS1257 < 0
    || _M0L6_2atmpS1257 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 31 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1257] = _M0L6_2atmpS1258;
  _M0L6_2atmpS1260 = 77;
  _M0L6_2atmpS1262 = 75;
  _M0L6_2atmpS1261 = (uint16_t)_M0L6_2atmpS1262;
  if (
    _M0L6_2atmpS1260 < 0
    || _M0L6_2atmpS1260 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 32 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1260] = _M0L6_2atmpS1261;
  _M0L6_2atmpS1263 = 82;
  _M0L6_2atmpS1265 = 89;
  _M0L6_2atmpS1264 = (uint16_t)_M0L6_2atmpS1265;
  if (
    _M0L6_2atmpS1263 < 0
    || _M0L6_2atmpS1263 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 33 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1263] = _M0L6_2atmpS1264;
  _M0L6_2atmpS1266 = 87;
  _M0L6_2atmpS1268 = 87;
  _M0L6_2atmpS1267 = (uint16_t)_M0L6_2atmpS1268;
  if (
    _M0L6_2atmpS1266 < 0
    || _M0L6_2atmpS1266 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 34 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1266] = _M0L6_2atmpS1267;
  _M0L6_2atmpS1269 = 83;
  _M0L6_2atmpS1271 = 83;
  _M0L6_2atmpS1270 = (uint16_t)_M0L6_2atmpS1271;
  if (
    _M0L6_2atmpS1269 < 0
    || _M0L6_2atmpS1269 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 35 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1269] = _M0L6_2atmpS1270;
  _M0L6_2atmpS1272 = 89;
  _M0L6_2atmpS1274 = 82;
  _M0L6_2atmpS1273 = (uint16_t)_M0L6_2atmpS1274;
  if (
    _M0L6_2atmpS1272 < 0
    || _M0L6_2atmpS1272 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 36 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1272] = _M0L6_2atmpS1273;
  _M0L6_2atmpS1275 = 75;
  _M0L6_2atmpS1277 = 77;
  _M0L6_2atmpS1276 = (uint16_t)_M0L6_2atmpS1277;
  if (
    _M0L6_2atmpS1275 < 0
    || _M0L6_2atmpS1275 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 37 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1275] = _M0L6_2atmpS1276;
  _M0L6_2atmpS1278 = 86;
  _M0L6_2atmpS1280 = 66;
  _M0L6_2atmpS1279 = (uint16_t)_M0L6_2atmpS1280;
  if (
    _M0L6_2atmpS1278 < 0
    || _M0L6_2atmpS1278 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 38 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1278] = _M0L6_2atmpS1279;
  _M0L6_2atmpS1281 = 72;
  _M0L6_2atmpS1283 = 68;
  _M0L6_2atmpS1282 = (uint16_t)_M0L6_2atmpS1283;
  if (
    _M0L6_2atmpS1281 < 0
    || _M0L6_2atmpS1281 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 39 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1281] = _M0L6_2atmpS1282;
  _M0L6_2atmpS1284 = 68;
  _M0L6_2atmpS1286 = 72;
  _M0L6_2atmpS1285 = (uint16_t)_M0L6_2atmpS1286;
  if (
    _M0L6_2atmpS1284 < 0
    || _M0L6_2atmpS1284 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 40 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1284] = _M0L6_2atmpS1285;
  _M0L6_2atmpS1287 = 66;
  _M0L6_2atmpS1289 = 86;
  _M0L6_2atmpS1288 = (uint16_t)_M0L6_2atmpS1289;
  if (
    _M0L6_2atmpS1287 < 0
    || _M0L6_2atmpS1287 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 41 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1287] = _M0L6_2atmpS1288;
  _M0L6_2atmpS1290 = 88;
  _M0L6_2atmpS1292 = 88;
  _M0L6_2atmpS1291 = (uint16_t)_M0L6_2atmpS1292;
  if (
    _M0L6_2atmpS1290 < 0
    || _M0L6_2atmpS1290 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 42 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1290] = _M0L6_2atmpS1291;
  _M0L6_2atmpS1293 = 78;
  _M0L6_2atmpS1295 = 78;
  _M0L6_2atmpS1294 = (uint16_t)_M0L6_2atmpS1295;
  if (
    _M0L6_2atmpS1293 < 0
    || _M0L6_2atmpS1293 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 43 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1293] = _M0L6_2atmpS1294;
  _M0L6_2atmpS1296 = 97;
  _M0L6_2atmpS1298 = 116;
  _M0L6_2atmpS1297 = (uint16_t)_M0L6_2atmpS1298;
  if (
    _M0L6_2atmpS1296 < 0
    || _M0L6_2atmpS1296 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 45 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1296] = _M0L6_2atmpS1297;
  _M0L6_2atmpS1299 = 99;
  _M0L6_2atmpS1301 = 103;
  _M0L6_2atmpS1300 = (uint16_t)_M0L6_2atmpS1301;
  if (
    _M0L6_2atmpS1299 < 0
    || _M0L6_2atmpS1299 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 46 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1299] = _M0L6_2atmpS1300;
  _M0L6_2atmpS1302 = 103;
  _M0L6_2atmpS1304 = 99;
  _M0L6_2atmpS1303 = (uint16_t)_M0L6_2atmpS1304;
  if (
    _M0L6_2atmpS1302 < 0
    || _M0L6_2atmpS1302 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 47 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1302] = _M0L6_2atmpS1303;
  _M0L6_2atmpS1305 = 116;
  _M0L6_2atmpS1307 = 97;
  _M0L6_2atmpS1306 = (uint16_t)_M0L6_2atmpS1307;
  if (
    _M0L6_2atmpS1305 < 0
    || _M0L6_2atmpS1305 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 48 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1305] = _M0L6_2atmpS1306;
  _M0L6_2atmpS1308 = 117;
  _M0L6_2atmpS1310 = 97;
  _M0L6_2atmpS1309 = (uint16_t)_M0L6_2atmpS1310;
  if (
    _M0L6_2atmpS1308 < 0
    || _M0L6_2atmpS1308 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 49 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1308] = _M0L6_2atmpS1309;
  _M0L6_2atmpS1311 = 109;
  _M0L6_2atmpS1313 = 107;
  _M0L6_2atmpS1312 = (uint16_t)_M0L6_2atmpS1313;
  if (
    _M0L6_2atmpS1311 < 0
    || _M0L6_2atmpS1311 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 50 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1311] = _M0L6_2atmpS1312;
  _M0L6_2atmpS1314 = 114;
  _M0L6_2atmpS1316 = 121;
  _M0L6_2atmpS1315 = (uint16_t)_M0L6_2atmpS1316;
  if (
    _M0L6_2atmpS1314 < 0
    || _M0L6_2atmpS1314 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 51 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1314] = _M0L6_2atmpS1315;
  _M0L6_2atmpS1317 = 119;
  _M0L6_2atmpS1319 = 119;
  _M0L6_2atmpS1318 = (uint16_t)_M0L6_2atmpS1319;
  if (
    _M0L6_2atmpS1317 < 0
    || _M0L6_2atmpS1317 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 52 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1317] = _M0L6_2atmpS1318;
  _M0L6_2atmpS1320 = 115;
  _M0L6_2atmpS1322 = 115;
  _M0L6_2atmpS1321 = (uint16_t)_M0L6_2atmpS1322;
  if (
    _M0L6_2atmpS1320 < 0
    || _M0L6_2atmpS1320 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 53 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1320] = _M0L6_2atmpS1321;
  _M0L6_2atmpS1323 = 121;
  _M0L6_2atmpS1325 = 114;
  _M0L6_2atmpS1324 = (uint16_t)_M0L6_2atmpS1325;
  if (
    _M0L6_2atmpS1323 < 0
    || _M0L6_2atmpS1323 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 54 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1323] = _M0L6_2atmpS1324;
  _M0L6_2atmpS1326 = 107;
  _M0L6_2atmpS1328 = 109;
  _M0L6_2atmpS1327 = (uint16_t)_M0L6_2atmpS1328;
  if (
    _M0L6_2atmpS1326 < 0
    || _M0L6_2atmpS1326 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 55 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1326] = _M0L6_2atmpS1327;
  _M0L6_2atmpS1329 = 118;
  _M0L6_2atmpS1331 = 98;
  _M0L6_2atmpS1330 = (uint16_t)_M0L6_2atmpS1331;
  if (
    _M0L6_2atmpS1329 < 0
    || _M0L6_2atmpS1329 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 56 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1329] = _M0L6_2atmpS1330;
  _M0L6_2atmpS1332 = 104;
  _M0L6_2atmpS1334 = 100;
  _M0L6_2atmpS1333 = (uint16_t)_M0L6_2atmpS1334;
  if (
    _M0L6_2atmpS1332 < 0
    || _M0L6_2atmpS1332 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 57 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1332] = _M0L6_2atmpS1333;
  _M0L6_2atmpS1335 = 100;
  _M0L6_2atmpS1337 = 104;
  _M0L6_2atmpS1336 = (uint16_t)_M0L6_2atmpS1337;
  if (
    _M0L6_2atmpS1335 < 0
    || _M0L6_2atmpS1335 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 58 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1335] = _M0L6_2atmpS1336;
  _M0L6_2atmpS1338 = 98;
  _M0L6_2atmpS1340 = 118;
  _M0L6_2atmpS1339 = (uint16_t)_M0L6_2atmpS1340;
  if (
    _M0L6_2atmpS1338 < 0
    || _M0L6_2atmpS1338 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 59 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1338] = _M0L6_2atmpS1339;
  _M0L6_2atmpS1341 = 120;
  _M0L6_2atmpS1343 = 120;
  _M0L6_2atmpS1342 = (uint16_t)_M0L6_2atmpS1343;
  if (
    _M0L6_2atmpS1341 < 0
    || _M0L6_2atmpS1341 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 60 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1341] = _M0L6_2atmpS1342;
  _M0L6_2atmpS1344 = 110;
  _M0L6_2atmpS1346 = 110;
  _M0L6_2atmpS1345 = (uint16_t)_M0L6_2atmpS1346;
  if (
    _M0L6_2atmpS1344 < 0
    || _M0L6_2atmpS1344 >= Moonbit_array_length(_M0L1tS813)
  ) {
    #line 61 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS813[_M0L6_2atmpS1344] = _M0L6_2atmpS1345;
  _M0FP26biolab8bio__seq20dna__complement__arr = _M0L1tS813;
  _M0L1tS936 = (uint16_t*)moonbit_make_string(128, 0);
  _M0L6_2atmpS1347 = 65;
  _M0L6_2atmpS1349 = 85;
  _M0L6_2atmpS1348 = (uint16_t)_M0L6_2atmpS1349;
  if (
    _M0L6_2atmpS1347 < 0
    || _M0L6_2atmpS1347 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 72 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1347] = _M0L6_2atmpS1348;
  _M0L6_2atmpS1350 = 67;
  _M0L6_2atmpS1352 = 71;
  _M0L6_2atmpS1351 = (uint16_t)_M0L6_2atmpS1352;
  if (
    _M0L6_2atmpS1350 < 0
    || _M0L6_2atmpS1350 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 73 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1350] = _M0L6_2atmpS1351;
  _M0L6_2atmpS1353 = 71;
  _M0L6_2atmpS1355 = 67;
  _M0L6_2atmpS1354 = (uint16_t)_M0L6_2atmpS1355;
  if (
    _M0L6_2atmpS1353 < 0
    || _M0L6_2atmpS1353 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 74 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1353] = _M0L6_2atmpS1354;
  _M0L6_2atmpS1356 = 85;
  _M0L6_2atmpS1358 = 65;
  _M0L6_2atmpS1357 = (uint16_t)_M0L6_2atmpS1358;
  if (
    _M0L6_2atmpS1356 < 0
    || _M0L6_2atmpS1356 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 75 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1356] = _M0L6_2atmpS1357;
  _M0L6_2atmpS1359 = 84;
  _M0L6_2atmpS1361 = 65;
  _M0L6_2atmpS1360 = (uint16_t)_M0L6_2atmpS1361;
  if (
    _M0L6_2atmpS1359 < 0
    || _M0L6_2atmpS1359 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 76 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1359] = _M0L6_2atmpS1360;
  _M0L6_2atmpS1362 = 77;
  _M0L6_2atmpS1364 = 75;
  _M0L6_2atmpS1363 = (uint16_t)_M0L6_2atmpS1364;
  if (
    _M0L6_2atmpS1362 < 0
    || _M0L6_2atmpS1362 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 77 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1362] = _M0L6_2atmpS1363;
  _M0L6_2atmpS1365 = 82;
  _M0L6_2atmpS1367 = 89;
  _M0L6_2atmpS1366 = (uint16_t)_M0L6_2atmpS1367;
  if (
    _M0L6_2atmpS1365 < 0
    || _M0L6_2atmpS1365 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 78 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1365] = _M0L6_2atmpS1366;
  _M0L6_2atmpS1368 = 87;
  _M0L6_2atmpS1370 = 87;
  _M0L6_2atmpS1369 = (uint16_t)_M0L6_2atmpS1370;
  if (
    _M0L6_2atmpS1368 < 0
    || _M0L6_2atmpS1368 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 79 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1368] = _M0L6_2atmpS1369;
  _M0L6_2atmpS1371 = 83;
  _M0L6_2atmpS1373 = 83;
  _M0L6_2atmpS1372 = (uint16_t)_M0L6_2atmpS1373;
  if (
    _M0L6_2atmpS1371 < 0
    || _M0L6_2atmpS1371 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 80 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1371] = _M0L6_2atmpS1372;
  _M0L6_2atmpS1374 = 89;
  _M0L6_2atmpS1376 = 82;
  _M0L6_2atmpS1375 = (uint16_t)_M0L6_2atmpS1376;
  if (
    _M0L6_2atmpS1374 < 0
    || _M0L6_2atmpS1374 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 81 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1374] = _M0L6_2atmpS1375;
  _M0L6_2atmpS1377 = 75;
  _M0L6_2atmpS1379 = 77;
  _M0L6_2atmpS1378 = (uint16_t)_M0L6_2atmpS1379;
  if (
    _M0L6_2atmpS1377 < 0
    || _M0L6_2atmpS1377 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 82 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1377] = _M0L6_2atmpS1378;
  _M0L6_2atmpS1380 = 86;
  _M0L6_2atmpS1382 = 66;
  _M0L6_2atmpS1381 = (uint16_t)_M0L6_2atmpS1382;
  if (
    _M0L6_2atmpS1380 < 0
    || _M0L6_2atmpS1380 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 83 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1380] = _M0L6_2atmpS1381;
  _M0L6_2atmpS1383 = 72;
  _M0L6_2atmpS1385 = 68;
  _M0L6_2atmpS1384 = (uint16_t)_M0L6_2atmpS1385;
  if (
    _M0L6_2atmpS1383 < 0
    || _M0L6_2atmpS1383 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 84 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1383] = _M0L6_2atmpS1384;
  _M0L6_2atmpS1386 = 68;
  _M0L6_2atmpS1388 = 72;
  _M0L6_2atmpS1387 = (uint16_t)_M0L6_2atmpS1388;
  if (
    _M0L6_2atmpS1386 < 0
    || _M0L6_2atmpS1386 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 85 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1386] = _M0L6_2atmpS1387;
  _M0L6_2atmpS1389 = 66;
  _M0L6_2atmpS1391 = 86;
  _M0L6_2atmpS1390 = (uint16_t)_M0L6_2atmpS1391;
  if (
    _M0L6_2atmpS1389 < 0
    || _M0L6_2atmpS1389 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 86 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1389] = _M0L6_2atmpS1390;
  _M0L6_2atmpS1392 = 88;
  _M0L6_2atmpS1394 = 88;
  _M0L6_2atmpS1393 = (uint16_t)_M0L6_2atmpS1394;
  if (
    _M0L6_2atmpS1392 < 0
    || _M0L6_2atmpS1392 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 87 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1392] = _M0L6_2atmpS1393;
  _M0L6_2atmpS1395 = 78;
  _M0L6_2atmpS1397 = 78;
  _M0L6_2atmpS1396 = (uint16_t)_M0L6_2atmpS1397;
  if (
    _M0L6_2atmpS1395 < 0
    || _M0L6_2atmpS1395 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 88 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1395] = _M0L6_2atmpS1396;
  _M0L6_2atmpS1398 = 97;
  _M0L6_2atmpS1400 = 117;
  _M0L6_2atmpS1399 = (uint16_t)_M0L6_2atmpS1400;
  if (
    _M0L6_2atmpS1398 < 0
    || _M0L6_2atmpS1398 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 90 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1398] = _M0L6_2atmpS1399;
  _M0L6_2atmpS1401 = 99;
  _M0L6_2atmpS1403 = 103;
  _M0L6_2atmpS1402 = (uint16_t)_M0L6_2atmpS1403;
  if (
    _M0L6_2atmpS1401 < 0
    || _M0L6_2atmpS1401 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 91 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1401] = _M0L6_2atmpS1402;
  _M0L6_2atmpS1404 = 103;
  _M0L6_2atmpS1406 = 99;
  _M0L6_2atmpS1405 = (uint16_t)_M0L6_2atmpS1406;
  if (
    _M0L6_2atmpS1404 < 0
    || _M0L6_2atmpS1404 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 92 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1404] = _M0L6_2atmpS1405;
  _M0L6_2atmpS1407 = 117;
  _M0L6_2atmpS1409 = 97;
  _M0L6_2atmpS1408 = (uint16_t)_M0L6_2atmpS1409;
  if (
    _M0L6_2atmpS1407 < 0
    || _M0L6_2atmpS1407 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 93 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1407] = _M0L6_2atmpS1408;
  _M0L6_2atmpS1410 = 116;
  _M0L6_2atmpS1412 = 97;
  _M0L6_2atmpS1411 = (uint16_t)_M0L6_2atmpS1412;
  if (
    _M0L6_2atmpS1410 < 0
    || _M0L6_2atmpS1410 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 94 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1410] = _M0L6_2atmpS1411;
  _M0L6_2atmpS1413 = 109;
  _M0L6_2atmpS1415 = 107;
  _M0L6_2atmpS1414 = (uint16_t)_M0L6_2atmpS1415;
  if (
    _M0L6_2atmpS1413 < 0
    || _M0L6_2atmpS1413 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 95 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1413] = _M0L6_2atmpS1414;
  _M0L6_2atmpS1416 = 114;
  _M0L6_2atmpS1418 = 121;
  _M0L6_2atmpS1417 = (uint16_t)_M0L6_2atmpS1418;
  if (
    _M0L6_2atmpS1416 < 0
    || _M0L6_2atmpS1416 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 96 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1416] = _M0L6_2atmpS1417;
  _M0L6_2atmpS1419 = 119;
  _M0L6_2atmpS1421 = 119;
  _M0L6_2atmpS1420 = (uint16_t)_M0L6_2atmpS1421;
  if (
    _M0L6_2atmpS1419 < 0
    || _M0L6_2atmpS1419 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 97 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1419] = _M0L6_2atmpS1420;
  _M0L6_2atmpS1422 = 115;
  _M0L6_2atmpS1424 = 115;
  _M0L6_2atmpS1423 = (uint16_t)_M0L6_2atmpS1424;
  if (
    _M0L6_2atmpS1422 < 0
    || _M0L6_2atmpS1422 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 98 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1422] = _M0L6_2atmpS1423;
  _M0L6_2atmpS1425 = 121;
  _M0L6_2atmpS1427 = 114;
  _M0L6_2atmpS1426 = (uint16_t)_M0L6_2atmpS1427;
  if (
    _M0L6_2atmpS1425 < 0
    || _M0L6_2atmpS1425 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 99 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1425] = _M0L6_2atmpS1426;
  _M0L6_2atmpS1428 = 107;
  _M0L6_2atmpS1430 = 109;
  _M0L6_2atmpS1429 = (uint16_t)_M0L6_2atmpS1430;
  if (
    _M0L6_2atmpS1428 < 0
    || _M0L6_2atmpS1428 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 100 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1428] = _M0L6_2atmpS1429;
  _M0L6_2atmpS1431 = 118;
  _M0L6_2atmpS1433 = 98;
  _M0L6_2atmpS1432 = (uint16_t)_M0L6_2atmpS1433;
  if (
    _M0L6_2atmpS1431 < 0
    || _M0L6_2atmpS1431 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 101 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1431] = _M0L6_2atmpS1432;
  _M0L6_2atmpS1434 = 104;
  _M0L6_2atmpS1436 = 100;
  _M0L6_2atmpS1435 = (uint16_t)_M0L6_2atmpS1436;
  if (
    _M0L6_2atmpS1434 < 0
    || _M0L6_2atmpS1434 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 102 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1434] = _M0L6_2atmpS1435;
  _M0L6_2atmpS1437 = 100;
  _M0L6_2atmpS1439 = 104;
  _M0L6_2atmpS1438 = (uint16_t)_M0L6_2atmpS1439;
  if (
    _M0L6_2atmpS1437 < 0
    || _M0L6_2atmpS1437 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 103 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1437] = _M0L6_2atmpS1438;
  _M0L6_2atmpS1440 = 98;
  _M0L6_2atmpS1442 = 118;
  _M0L6_2atmpS1441 = (uint16_t)_M0L6_2atmpS1442;
  if (
    _M0L6_2atmpS1440 < 0
    || _M0L6_2atmpS1440 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 104 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1440] = _M0L6_2atmpS1441;
  _M0L6_2atmpS1443 = 120;
  _M0L6_2atmpS1445 = 120;
  _M0L6_2atmpS1444 = (uint16_t)_M0L6_2atmpS1445;
  if (
    _M0L6_2atmpS1443 < 0
    || _M0L6_2atmpS1443 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 105 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1443] = _M0L6_2atmpS1444;
  _M0L6_2atmpS1446 = 110;
  _M0L6_2atmpS1448 = 110;
  _M0L6_2atmpS1447 = (uint16_t)_M0L6_2atmpS1448;
  if (
    _M0L6_2atmpS1446 < 0
    || _M0L6_2atmpS1446 >= Moonbit_array_length(_M0L1tS936)
  ) {
    #line 106 "/home/developer/Documents/1/complement.mbt"
    moonbit_panic();
  }
  _M0L1tS936[_M0L6_2atmpS1446] = _M0L6_2atmpS1447;
  _M0FP26biolab8bio__seq20rna__complement__arr = _M0L1tS936;
  _M0L7_2abindS1003 = (struct _M0TUicE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1637 = _M0L7_2abindS1003;
  _M0L6_2atmpS1636
  = (struct _M0TPB9ArrayViewGUicEE){
    .$0 = _M0L6_2atmpS1637, .$1 = 0, .$2 = 0
  };
  #line 15 "/home/developer/Documents/1/codon_table.mbt"
  _M0L1mS1002 = _M0MPB3Map3MapGicE(_M0L6_2atmpS1636, 64ll);
  moonbit_decref(_M0L6_2atmpS1636.$0);
  _M0L3regS1004 = _M0L1mS1002;
  _M0L6_2atmpS1453 = 84;
  _M0L6_2atmpS1454 = 84;
  _M0L6_2atmpS1455 = 84;
  moonbit_incref(_M0L3regS1004);
  #line 21 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1453, _M0L6_2atmpS1454, _M0L6_2atmpS1455, 70);
  _M0L6_2atmpS1456 = 84;
  _M0L6_2atmpS1457 = 84;
  _M0L6_2atmpS1458 = 67;
  #line 22 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1456, _M0L6_2atmpS1457, _M0L6_2atmpS1458, 70);
  _M0L6_2atmpS1459 = 84;
  _M0L6_2atmpS1460 = 84;
  _M0L6_2atmpS1461 = 65;
  #line 24 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1459, _M0L6_2atmpS1460, _M0L6_2atmpS1461, 76);
  _M0L6_2atmpS1462 = 84;
  _M0L6_2atmpS1463 = 84;
  _M0L6_2atmpS1464 = 71;
  #line 25 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1462, _M0L6_2atmpS1463, _M0L6_2atmpS1464, 76);
  _M0L6_2atmpS1465 = 67;
  _M0L6_2atmpS1466 = 84;
  _M0L6_2atmpS1467 = 84;
  #line 26 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1465, _M0L6_2atmpS1466, _M0L6_2atmpS1467, 76);
  _M0L6_2atmpS1468 = 67;
  _M0L6_2atmpS1469 = 84;
  _M0L6_2atmpS1470 = 67;
  #line 27 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1468, _M0L6_2atmpS1469, _M0L6_2atmpS1470, 76);
  _M0L6_2atmpS1471 = 67;
  _M0L6_2atmpS1472 = 84;
  _M0L6_2atmpS1473 = 65;
  #line 28 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1471, _M0L6_2atmpS1472, _M0L6_2atmpS1473, 76);
  _M0L6_2atmpS1474 = 67;
  _M0L6_2atmpS1475 = 84;
  _M0L6_2atmpS1476 = 71;
  #line 29 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1474, _M0L6_2atmpS1475, _M0L6_2atmpS1476, 76);
  _M0L6_2atmpS1477 = 65;
  _M0L6_2atmpS1478 = 84;
  _M0L6_2atmpS1479 = 84;
  #line 31 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1477, _M0L6_2atmpS1478, _M0L6_2atmpS1479, 73);
  _M0L6_2atmpS1480 = 65;
  _M0L6_2atmpS1481 = 84;
  _M0L6_2atmpS1482 = 67;
  #line 32 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1480, _M0L6_2atmpS1481, _M0L6_2atmpS1482, 73);
  _M0L6_2atmpS1483 = 65;
  _M0L6_2atmpS1484 = 84;
  _M0L6_2atmpS1485 = 65;
  #line 33 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1483, _M0L6_2atmpS1484, _M0L6_2atmpS1485, 73);
  _M0L6_2atmpS1486 = 65;
  _M0L6_2atmpS1487 = 84;
  _M0L6_2atmpS1488 = 71;
  #line 35 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1486, _M0L6_2atmpS1487, _M0L6_2atmpS1488, 77);
  _M0L6_2atmpS1489 = 71;
  _M0L6_2atmpS1490 = 84;
  _M0L6_2atmpS1491 = 84;
  #line 37 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1489, _M0L6_2atmpS1490, _M0L6_2atmpS1491, 86);
  _M0L6_2atmpS1492 = 71;
  _M0L6_2atmpS1493 = 84;
  _M0L6_2atmpS1494 = 67;
  #line 38 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1492, _M0L6_2atmpS1493, _M0L6_2atmpS1494, 86);
  _M0L6_2atmpS1495 = 71;
  _M0L6_2atmpS1496 = 84;
  _M0L6_2atmpS1497 = 65;
  #line 39 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1495, _M0L6_2atmpS1496, _M0L6_2atmpS1497, 86);
  _M0L6_2atmpS1498 = 71;
  _M0L6_2atmpS1499 = 84;
  _M0L6_2atmpS1500 = 71;
  #line 40 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1498, _M0L6_2atmpS1499, _M0L6_2atmpS1500, 86);
  _M0L6_2atmpS1501 = 84;
  _M0L6_2atmpS1502 = 67;
  _M0L6_2atmpS1503 = 84;
  #line 42 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1501, _M0L6_2atmpS1502, _M0L6_2atmpS1503, 83);
  _M0L6_2atmpS1504 = 84;
  _M0L6_2atmpS1505 = 67;
  _M0L6_2atmpS1506 = 67;
  #line 43 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1504, _M0L6_2atmpS1505, _M0L6_2atmpS1506, 83);
  _M0L6_2atmpS1507 = 84;
  _M0L6_2atmpS1508 = 67;
  _M0L6_2atmpS1509 = 65;
  #line 44 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1507, _M0L6_2atmpS1508, _M0L6_2atmpS1509, 83);
  _M0L6_2atmpS1510 = 84;
  _M0L6_2atmpS1511 = 67;
  _M0L6_2atmpS1512 = 71;
  #line 45 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1510, _M0L6_2atmpS1511, _M0L6_2atmpS1512, 83);
  _M0L6_2atmpS1513 = 67;
  _M0L6_2atmpS1514 = 67;
  _M0L6_2atmpS1515 = 84;
  #line 47 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1513, _M0L6_2atmpS1514, _M0L6_2atmpS1515, 80);
  _M0L6_2atmpS1516 = 67;
  _M0L6_2atmpS1517 = 67;
  _M0L6_2atmpS1518 = 67;
  #line 48 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1516, _M0L6_2atmpS1517, _M0L6_2atmpS1518, 80);
  _M0L6_2atmpS1519 = 67;
  _M0L6_2atmpS1520 = 67;
  _M0L6_2atmpS1521 = 65;
  #line 49 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1519, _M0L6_2atmpS1520, _M0L6_2atmpS1521, 80);
  _M0L6_2atmpS1522 = 67;
  _M0L6_2atmpS1523 = 67;
  _M0L6_2atmpS1524 = 71;
  #line 50 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1522, _M0L6_2atmpS1523, _M0L6_2atmpS1524, 80);
  _M0L6_2atmpS1525 = 65;
  _M0L6_2atmpS1526 = 67;
  _M0L6_2atmpS1527 = 84;
  #line 52 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1525, _M0L6_2atmpS1526, _M0L6_2atmpS1527, 84);
  _M0L6_2atmpS1528 = 65;
  _M0L6_2atmpS1529 = 67;
  _M0L6_2atmpS1530 = 67;
  #line 53 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1528, _M0L6_2atmpS1529, _M0L6_2atmpS1530, 84);
  _M0L6_2atmpS1531 = 65;
  _M0L6_2atmpS1532 = 67;
  _M0L6_2atmpS1533 = 65;
  #line 54 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1531, _M0L6_2atmpS1532, _M0L6_2atmpS1533, 84);
  _M0L6_2atmpS1534 = 65;
  _M0L6_2atmpS1535 = 67;
  _M0L6_2atmpS1536 = 71;
  #line 55 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1534, _M0L6_2atmpS1535, _M0L6_2atmpS1536, 84);
  _M0L6_2atmpS1537 = 71;
  _M0L6_2atmpS1538 = 67;
  _M0L6_2atmpS1539 = 84;
  #line 57 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1537, _M0L6_2atmpS1538, _M0L6_2atmpS1539, 65);
  _M0L6_2atmpS1540 = 71;
  _M0L6_2atmpS1541 = 67;
  _M0L6_2atmpS1542 = 67;
  #line 58 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1540, _M0L6_2atmpS1541, _M0L6_2atmpS1542, 65);
  _M0L6_2atmpS1543 = 71;
  _M0L6_2atmpS1544 = 67;
  _M0L6_2atmpS1545 = 65;
  #line 59 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1543, _M0L6_2atmpS1544, _M0L6_2atmpS1545, 65);
  _M0L6_2atmpS1546 = 71;
  _M0L6_2atmpS1547 = 67;
  _M0L6_2atmpS1548 = 71;
  #line 60 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1546, _M0L6_2atmpS1547, _M0L6_2atmpS1548, 65);
  _M0L6_2atmpS1549 = 84;
  _M0L6_2atmpS1550 = 65;
  _M0L6_2atmpS1551 = 84;
  #line 62 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1549, _M0L6_2atmpS1550, _M0L6_2atmpS1551, 89);
  _M0L6_2atmpS1552 = 84;
  _M0L6_2atmpS1553 = 65;
  _M0L6_2atmpS1554 = 67;
  #line 63 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1552, _M0L6_2atmpS1553, _M0L6_2atmpS1554, 89);
  _M0L6_2atmpS1555 = 67;
  _M0L6_2atmpS1556 = 65;
  _M0L6_2atmpS1557 = 84;
  #line 65 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1555, _M0L6_2atmpS1556, _M0L6_2atmpS1557, 72);
  _M0L6_2atmpS1558 = 67;
  _M0L6_2atmpS1559 = 65;
  _M0L6_2atmpS1560 = 67;
  #line 66 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1558, _M0L6_2atmpS1559, _M0L6_2atmpS1560, 72);
  _M0L6_2atmpS1561 = 67;
  _M0L6_2atmpS1562 = 65;
  _M0L6_2atmpS1563 = 65;
  #line 68 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1561, _M0L6_2atmpS1562, _M0L6_2atmpS1563, 81);
  _M0L6_2atmpS1564 = 67;
  _M0L6_2atmpS1565 = 65;
  _M0L6_2atmpS1566 = 71;
  #line 69 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1564, _M0L6_2atmpS1565, _M0L6_2atmpS1566, 81);
  _M0L6_2atmpS1567 = 65;
  _M0L6_2atmpS1568 = 65;
  _M0L6_2atmpS1569 = 84;
  #line 71 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1567, _M0L6_2atmpS1568, _M0L6_2atmpS1569, 78);
  _M0L6_2atmpS1570 = 65;
  _M0L6_2atmpS1571 = 65;
  _M0L6_2atmpS1572 = 67;
  #line 72 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1570, _M0L6_2atmpS1571, _M0L6_2atmpS1572, 78);
  _M0L6_2atmpS1573 = 65;
  _M0L6_2atmpS1574 = 65;
  _M0L6_2atmpS1575 = 65;
  #line 74 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1573, _M0L6_2atmpS1574, _M0L6_2atmpS1575, 75);
  _M0L6_2atmpS1576 = 65;
  _M0L6_2atmpS1577 = 65;
  _M0L6_2atmpS1578 = 71;
  #line 75 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1576, _M0L6_2atmpS1577, _M0L6_2atmpS1578, 75);
  _M0L6_2atmpS1579 = 71;
  _M0L6_2atmpS1580 = 65;
  _M0L6_2atmpS1581 = 84;
  #line 77 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1579, _M0L6_2atmpS1580, _M0L6_2atmpS1581, 68);
  _M0L6_2atmpS1582 = 71;
  _M0L6_2atmpS1583 = 65;
  _M0L6_2atmpS1584 = 67;
  #line 78 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1582, _M0L6_2atmpS1583, _M0L6_2atmpS1584, 68);
  _M0L6_2atmpS1585 = 71;
  _M0L6_2atmpS1586 = 65;
  _M0L6_2atmpS1587 = 65;
  #line 80 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1585, _M0L6_2atmpS1586, _M0L6_2atmpS1587, 69);
  _M0L6_2atmpS1588 = 71;
  _M0L6_2atmpS1589 = 65;
  _M0L6_2atmpS1590 = 71;
  #line 81 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1588, _M0L6_2atmpS1589, _M0L6_2atmpS1590, 69);
  _M0L6_2atmpS1591 = 84;
  _M0L6_2atmpS1592 = 71;
  _M0L6_2atmpS1593 = 84;
  #line 83 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1591, _M0L6_2atmpS1592, _M0L6_2atmpS1593, 67);
  _M0L6_2atmpS1594 = 84;
  _M0L6_2atmpS1595 = 71;
  _M0L6_2atmpS1596 = 67;
  #line 84 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1594, _M0L6_2atmpS1595, _M0L6_2atmpS1596, 67);
  _M0L6_2atmpS1597 = 84;
  _M0L6_2atmpS1598 = 71;
  _M0L6_2atmpS1599 = 71;
  #line 86 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1597, _M0L6_2atmpS1598, _M0L6_2atmpS1599, 87);
  _M0L6_2atmpS1600 = 67;
  _M0L6_2atmpS1601 = 71;
  _M0L6_2atmpS1602 = 84;
  #line 88 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1600, _M0L6_2atmpS1601, _M0L6_2atmpS1602, 82);
  _M0L6_2atmpS1603 = 67;
  _M0L6_2atmpS1604 = 71;
  _M0L6_2atmpS1605 = 67;
  #line 89 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1603, _M0L6_2atmpS1604, _M0L6_2atmpS1605, 82);
  _M0L6_2atmpS1606 = 67;
  _M0L6_2atmpS1607 = 71;
  _M0L6_2atmpS1608 = 65;
  #line 90 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1606, _M0L6_2atmpS1607, _M0L6_2atmpS1608, 82);
  _M0L6_2atmpS1609 = 67;
  _M0L6_2atmpS1610 = 71;
  _M0L6_2atmpS1611 = 71;
  #line 91 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1609, _M0L6_2atmpS1610, _M0L6_2atmpS1611, 82);
  _M0L6_2atmpS1612 = 65;
  _M0L6_2atmpS1613 = 71;
  _M0L6_2atmpS1614 = 65;
  #line 92 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1612, _M0L6_2atmpS1613, _M0L6_2atmpS1614, 82);
  _M0L6_2atmpS1615 = 65;
  _M0L6_2atmpS1616 = 71;
  _M0L6_2atmpS1617 = 71;
  #line 93 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1615, _M0L6_2atmpS1616, _M0L6_2atmpS1617, 82);
  _M0L6_2atmpS1618 = 71;
  _M0L6_2atmpS1619 = 71;
  _M0L6_2atmpS1620 = 84;
  #line 95 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1618, _M0L6_2atmpS1619, _M0L6_2atmpS1620, 71);
  _M0L6_2atmpS1621 = 71;
  _M0L6_2atmpS1622 = 71;
  _M0L6_2atmpS1623 = 67;
  #line 96 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1621, _M0L6_2atmpS1622, _M0L6_2atmpS1623, 71);
  _M0L6_2atmpS1624 = 71;
  _M0L6_2atmpS1625 = 71;
  _M0L6_2atmpS1626 = 65;
  #line 97 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1624, _M0L6_2atmpS1625, _M0L6_2atmpS1626, 71);
  _M0L6_2atmpS1627 = 71;
  _M0L6_2atmpS1628 = 71;
  _M0L6_2atmpS1629 = 71;
  #line 98 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1627, _M0L6_2atmpS1628, _M0L6_2atmpS1629, 71);
  _M0L6_2atmpS1630 = 65;
  _M0L6_2atmpS1631 = 71;
  _M0L6_2atmpS1632 = 84;
  #line 100 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1630, _M0L6_2atmpS1631, _M0L6_2atmpS1632, 83);
  _M0L6_2atmpS1633 = 65;
  _M0L6_2atmpS1634 = 71;
  _M0L6_2atmpS1635 = 67;
  #line 101 "/home/developer/Documents/1/codon_table.mbt"
  _M0FP017____moonbit__initN3regS1004(_M0L3regS1004, _M0L6_2atmpS1633, _M0L6_2atmpS1634, _M0L6_2atmpS1635, 83);
  moonbit_decref(_M0L3regS1004);
  _M0FP26biolab8bio__seq21codon__forward__table = _M0L1mS1002;
  _M0L1tS967 = (uint8_t*)moonbit_make_bytes(128, 0);
  #line 111 "/home/developer/Documents/1/codon_table.mbt"
  _M0L5_2aitS968
  = _M0MPC16string6String4iter((moonbit_string_t)moonbit_string_literal_134.data);
  while (1) {
    int32_t _M0L1cS970;
    int32_t _M0L7_2abindS972;
    int32_t _M0L6_2atmpS1638;
    #line 111 "/home/developer/Documents/1/codon_table.mbt"
    _M0L7_2abindS972 = _M0MPB4Iter4nextGcE(_M0L5_2aitS968);
    if (_M0L7_2abindS972 == -1) {
      moonbit_decref(_M0L5_2aitS968);
    } else {
      int32_t _M0L7_2aSomeS973 = _M0L7_2abindS972;
      int32_t _M0L4_2acS974 = _M0L7_2aSomeS973;
      _M0L1cS970 = _M0L4_2acS974;
      goto join_969;
    }
    goto joinlet_3491;
    join_969:;
    _M0L6_2atmpS1638 = _M0L1cS970;
    if (
      _M0L6_2atmpS1638 < 0
      || _M0L6_2atmpS1638 >= Moonbit_array_length(_M0L1tS967)
    ) {
      #line 112 "/home/developer/Documents/1/codon_table.mbt"
      moonbit_panic();
    }
    _M0L1tS967[_M0L6_2atmpS1638] = 1;
    continue;
    joinlet_3491:;
    break;
  }
  _M0FP26biolab8bio__seq20valid__letter__table = _M0L1tS967;
}

int main(int argc, char** argv) {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** _M0L6_2atmpS1244;
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1209;
  struct _M0TPB5ArrayGUsiEE* _M0L7_2abindS1210;
  int32_t _M0L7_2abindS1211;
  int32_t _M0L2__S1212;
  moonbit_runtime_init(argc, argv);
  moonbit_init();
  _M0L6_2atmpS1244
  = (struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error**)moonbit_empty_ref_array;
  _M0L12async__testsS1209
  = (struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE));
  Moonbit_object_header(_M0L12async__testsS1209)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 0);
  _M0L12async__testsS1209->$0 = _M0L6_2atmpS1244;
  _M0L12async__testsS1209->$1 = 0;
  #line 399 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0L7_2abindS1210
  = _M0FP46biolab8bio__seq3cmd20main__blackbox__test52moonbit__test__driver__internal__native__parse__args();
  _M0L7_2abindS1211 = _M0L7_2abindS1210->$1;
  _M0L2__S1212 = 0;
  while (1) {
    if (_M0L2__S1212 < _M0L7_2abindS1211) {
      struct _M0TUsiE** _M0L3bufS1243 = _M0L7_2abindS1210->$0;
      struct _M0TUsiE* _M0L3argS1213 =
        (struct _M0TUsiE*)_M0L3bufS1243[_M0L2__S1212];
      moonbit_string_t _M0L6_2atmpS1240 = _M0L3argS1213->$0;
      int32_t _M0L6_2atmpS1241 = _M0L3argS1213->$1;
      int32_t _M0L6_2atmpS1242;
      moonbit_incref(_M0L6_2atmpS1240);
      #line 400 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
      _M0FP46biolab8bio__seq3cmd20main__blackbox__test44moonbit__test__driver__internal__do__execute(_M0L12async__testsS1209, _M0L6_2atmpS1240, _M0L6_2atmpS1241);
      moonbit_decref(_M0L6_2atmpS1240);
      _M0L6_2atmpS1242 = _M0L2__S1212 + 1;
      _M0L2__S1212 = _M0L6_2atmpS1242;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1210);
    }
    break;
  }
  #line 402 "/home/developer/Documents/1/cmd/main/__generated_driver_for_blackbox_test.mbt"
  _M0IP016_24default__implP46biolab8bio__seq3cmd20main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd20main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(_M0L12async__testsS1209);
  moonbit_decref(_M0L12async__testsS1209);
  return 0;
}