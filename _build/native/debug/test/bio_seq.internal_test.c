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
struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE;

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure;

struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidFormat;

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError;

struct _M0TPB4IterGsE;

struct _M0TWRPC15error5ErrorEs;

struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest;

struct _M0TWssbEu;

struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidNumber;

struct _M0TUsiE;

struct _M0TPB13StringBuilder;

struct _M0TPB17FloatingDecimal64;

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some;

struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE;

struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533;

struct _M0TUdiE;

struct _M0TWEuQRPC15error5Error;

struct _M0BTPB6Logger;

struct _M0DTPC16result6ResultGbRP26biolab8bio__seq33MoonBitTestDriverInternalSkipTestE3Err;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq8TreeNodeRP26biolab8bio__seq10ParseErrorE3Err;

struct _M0DTPC16result6ResultGUdiERP26biolab8bio__seq10ParseErrorE3Err;

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE;

struct _M0TPB6Logger;

struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE;

struct _M0TPB5ArrayGUsiEE;

struct _M0TPB8MutLocalGOsE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok;

struct _M0TP26biolab8bio__seq8TreeNode;

struct _M0TPB5ArrayGiE;

struct _M0TPB19MulShiftAll64Result;

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE;

struct _M0DTPC16result6ResultGUdiERP26biolab8bio__seq10ParseErrorE2Ok;

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__;

struct _M0DTPC16option6OptionGOsE4Some;

struct _M0TUiUWEuQRPC15error5ErrorNsEE;

struct _M0DTPC16result6ResultGuRPB7FailureE3Err;

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError;

struct _M0TWRPC15error5ErrorEu;

struct _M0TURPC16string10StringViewRPB6LoggerE;

struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError;

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__;

struct _M0TPB8MutLocalGiE;

struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538;

struct _M0TPB4Show;

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE;

struct _M0DTPC16result6ResultGuRPC15error5ErrorE2Ok;

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err;

struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE;

struct _M0TPB8MutLocalGdE;

struct _M0BTPB4Show;

struct _M0TWuEu;

struct _M0TPC16string10StringView;

struct _M0TPB8MutLocalGbE;

struct _M0TPB8MutLocalGOdE;

struct _M0DTPC16result6ResultGuRPB7FailureE2Ok;

struct _M0DTPC16result6ResultGURP26biolab8bio__seq8TreeNodeiERP26biolab8bio__seq10ParseErrorE2Ok;

struct _M0DTPC16result6ResultGURP26biolab8bio__seq8TreeNodeiERP26biolab8bio__seq10ParseErrorE3Err;

struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__;

struct _M0KTPB6LoggerTPB13StringBuilder;

struct _M0DTPC16result6ResultGbRP26biolab8bio__seq33MoonBitTestDriverInternalSkipTestE2Ok;

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq8TreeNodeRP26biolab8bio__seq10ParseErrorE2Ok;

struct _M0TPB5ArrayGsE;

struct _M0TUWEuQRPC15error5ErrorNsE;

struct _M0TURP26biolab8bio__seq8TreeNodeiE;

struct _M0DTPC16option6OptionGdE4Some;

struct _M0TWEu;

struct _M0TPB9ArrayViewGsE;

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE;

struct _M0TWEOs;

struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE;

struct _M0DTPC16result6ResultGuRPC15error5ErrorE3Err;

struct _M0TPB7Umul128;

struct _M0TPB8Pow5Pair;

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE {
  struct _M0TP26biolab8bio__seq8TreeNode** $0;
  int32_t $1;
  
};

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure {
  moonbit_string_t $0;
  
};

struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidFormat {
  moonbit_string_t $0;
  
};

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError {
  moonbit_string_t $0;
  
};

struct _M0TPB4IterGsE {
  struct _M0TWEOs* $0;
  int64_t $1;
  
};

struct _M0TWRPC15error5ErrorEs {
  moonbit_string_t(* code)(struct _M0TWRPC15error5ErrorEs*, void*);
  
};

struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest {
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

struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidNumber {
  moonbit_string_t $0;
  
};

struct _M0TUsiE {
  moonbit_string_t $0;
  int32_t $1;
  
};

struct _M0TPB13StringBuilder {
  uint16_t* $0;
  int32_t $1;
  
};

struct _M0TPB17FloatingDecimal64 {
  uint64_t $0;
  int32_t $1;
  
};

struct _M0TPC16string10StringView {
  moonbit_string_t $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some {
  struct _M0TPC16string10StringView $0;
  
};

struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** $0;
  int32_t $1;
  
};

struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533 {
  int32_t(* code)(struct _M0TWEu*);
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TUdiE {
  double $0;
  int32_t $1;
  
};

struct _M0TWEuQRPC15error5Error {
  struct moonbit_result_0(* code)(struct _M0TWEuQRPC15error5Error*);
  
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

struct _M0DTPC16result6ResultGRP26biolab8bio__seq8TreeNodeRP26biolab8bio__seq10ParseErrorE3Err {
  void* $0;
  
};

struct _M0DTPC16result6ResultGUdiERP26biolab8bio__seq10ParseErrorE3Err {
  void* $0;
  
};

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE {
  int32_t $0;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* $5;
  
};

struct _M0TPB6Logger {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE {
  moonbit_string_t $0;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* $1;
  
};

struct _M0TPB5ArrayGUsiEE {
  struct _M0TUsiE** $0;
  int32_t $1;
  
};

struct _M0TPB8MutLocalGOsE {
  moonbit_string_t $0;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok {
  int32_t $0;
  
};

struct _M0TP26biolab8bio__seq8TreeNode {
  moonbit_string_t $0;
  void* $1;
  void* $2;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* $3;
  
};

struct _M0TPB5ArrayGiE {
  int32_t* $0;
  int32_t $1;
  
};

struct _M0TPB19MulShiftAll64Result {
  uint64_t $0;
  uint64_t $1;
  uint64_t $2;
  
};

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** $0;
  int32_t $1;
  
};

struct _M0DTPC16result6ResultGUdiERP26biolab8bio__seq10ParseErrorE2Ok {
  struct _M0TUdiE* $0;
  
};

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__ {
  int32_t(* code)(struct _M0TWRPC15error5ErrorEu*, void*);
  struct _M0TWRPC15error5ErrorEs* $0;
  struct _M0TWssbEu* $1;
  moonbit_string_t $2;
  
};

struct _M0DTPC16option6OptionGOsE4Some {
  moonbit_string_t $0;
  
};

struct _M0TUiUWEuQRPC15error5ErrorNsEE {
  int32_t $0;
  struct _M0TUWEuQRPC15error5ErrorNsE* $1;
  
};

struct _M0DTPC16result6ResultGuRPB7FailureE3Err {
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

struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError {
  moonbit_string_t $0;
  
};

struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__ {
  int32_t(* code)(struct _M0TWEu*);
  struct _M0TWssbEu* $0;
  moonbit_string_t $1;
  
};

struct _M0TPB8MutLocalGiE {
  int32_t $0;
  
};

struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538 {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TPB4Show {
  struct _M0BTPB4Show* $0;
  void* $1;
  
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

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err {
  void* $0;
  
};

struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE {
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0TPB8MutLocalGdE {
  double $0;
  
};

struct _M0BTPB4Show {
  int32_t(* $method_0)(void*, struct _M0TPB6Logger);
  moonbit_string_t(* $method_1)(void*);
  
};

struct _M0TWuEu {
  int32_t(* code)(struct _M0TWuEu*, int32_t);
  
};

struct _M0TPB8MutLocalGbE {
  int32_t $0;
  
};

struct _M0TPB8MutLocalGOdE {
  void* $0;
  
};

struct _M0DTPC16result6ResultGuRPB7FailureE2Ok {
  int32_t $0;
  
};

struct _M0DTPC16result6ResultGURP26biolab8bio__seq8TreeNodeiERP26biolab8bio__seq10ParseErrorE2Ok {
  struct _M0TURP26biolab8bio__seq8TreeNodeiE* $0;
  
};

struct _M0DTPC16result6ResultGURP26biolab8bio__seq8TreeNodeiERP26biolab8bio__seq10ParseErrorE3Err {
  void* $0;
  
};

struct _M0TPB9ArrayViewGsE {
  moonbit_string_t* $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__ {
  moonbit_string_t(* code)(struct _M0TWEOs*);
  struct _M0TPB9ArrayViewGsE $0;
  int32_t $1;
  struct _M0TPB8MutLocalGiE* $2;
  
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

struct _M0DTPC16result6ResultGRP26biolab8bio__seq8TreeNodeRP26biolab8bio__seq10ParseErrorE2Ok {
  struct _M0TP26biolab8bio__seq8TreeNode* $0;
  
};

struct _M0TPB5ArrayGsE {
  moonbit_string_t* $0;
  int32_t $1;
  
};

struct _M0TUWEuQRPC15error5ErrorNsE {
  struct _M0TWEuQRPC15error5Error* $0;
  moonbit_string_t* $1;
  
};

struct _M0TURP26biolab8bio__seq8TreeNodeiE {
  struct _M0TP26biolab8bio__seq8TreeNode* $0;
  int32_t $1;
  
};

struct _M0DTPC16option6OptionGdE4Some {
  double $0;
  
};

struct _M0TWEu {
  int32_t(* code)(struct _M0TWEu*);
  
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

struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE {
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16result6ResultGuRPC15error5ErrorE3Err {
  void* $0;
  
};

struct _M0TPB7Umul128 {
  uint64_t $0;
  uint64_t $1;
  
};

struct _M0TPB8Pow5Pair {
  uint64_t $0;
  uint64_t $1;
  
};

struct moonbit_result_2 {
  int tag;
  union { struct _M0TURP26biolab8bio__seq8TreeNodeiE* ok; void* err;  } data;
  
};

struct moonbit_result_1 {
  int tag;
  union { struct _M0TP26biolab8bio__seq8TreeNode* ok; void* err;  } data;
  
};

struct moonbit_result_3 {
  int tag;
  union { struct _M0TUdiE* ok; void* err;  } data;
  
};

struct moonbit_result_0 {
  int tag;
  union { int32_t ok; void* err;  } data;
  
};

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__5_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__12_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__0_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__15_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__7_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__17_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__10_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__2_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__3_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__11_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__6_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__9_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__1_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__8_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__14_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__16_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__4_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__13_2edyncall(
  struct _M0TWEuQRPC15error5Error*
);

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t
);

moonbit_string_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1545(
  struct _M0TWRPC15error5ErrorEs*,
  void*
);

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1538(
  struct _M0TWssbEu*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN13handle__startS1533(
  struct _M0TWEu*
);

struct moonbit_result_0 _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__test(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3458l490(
  struct _M0TWEu*
);

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3454l491(
  struct _M0TWRPC15error5ErrorEu*,
  void*
);

int32_t _M0FP26biolab8bio__seq45moonbit__test__driver__internal__catch__error(
  struct _M0TWEuQRPC15error5Error*,
  struct _M0TWEu*,
  struct _M0TWRPC15error5ErrorEu*
);

struct _M0TPB5ArrayGUsiEE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__args(
  
);

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1444(
  int32_t,
  moonbit_string_t,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1437(
  int32_t
);

moonbit_string_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430(
  int32_t,
  moonbit_bytes_t
);

int32_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1423(
  int32_t,
  moonbit_string_t
);

#define _M0FP26biolab8bio__seq52moonbit__test__driver__internal__get__cli__args__ffi moonbit_get_cli_args

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP016_24default__implP26biolab8bio__seq28MoonBit__Async__Test__Driver17run__async__testsGRP26biolab8bio__seq34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__17(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__16(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__15(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__14(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__13(
  
);

double _M0MP26biolab8bio__seq8TreeNode13total__length(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__12(
  
);

void* _M0MP26biolab8bio__seq8TreeNode8distance(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  moonbit_string_t,
  moonbit_string_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode10find__path(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  moonbit_string_t
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__11(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__10(
  
);

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode3lca(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  struct _M0TPB5ArrayGsE*
);

int32_t _M0FP26biolab8bio__seq21find__path__backtrack(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  moonbit_string_t,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__9(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__8(
  
);

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode4find(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  moonbit_string_t
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__7(
  
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode10levelorder(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__6(
  
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode9postorder(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

int32_t _M0MP26biolab8bio__seq8TreeNode17postorder__helper(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__5(
  
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode8preorder(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

int32_t _M0MP26biolab8bio__seq8TreeNode16preorder__helper(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__4(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__3(
  
);

moonbit_string_t _M0MP26biolab8bio__seq8TreeNode10to__newick(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

int32_t _M0MP26biolab8bio__seq8TreeNode18to__newick__helper(
  struct _M0TP26biolab8bio__seq8TreeNode*,
  struct _M0TPB13StringBuilder*
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__2(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__1(
  
);

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__0(
  
);

struct moonbit_result_1 _M0MP26biolab8bio__seq8TreeNode12from__newick(
  moonbit_string_t
);

struct moonbit_result_1 _M0FP26biolab8bio__seq31parse__tree__node__from__newick(
  moonbit_string_t
);

struct moonbit_result_2 _M0FP26biolab8bio__seq16parse__node__rec(
  moonbit_string_t,
  int32_t,
  int32_t
);

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode3new(
  void*,
  void*,
  void*,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode11new_2einner(
  moonbit_string_t,
  void*,
  void*,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

int32_t _M0FP26biolab8bio__seq8skip__ws(moonbit_string_t, int32_t, int32_t);

struct moonbit_result_3 _M0FP26biolab8bio__seq18parse__number__rec(
  moonbit_string_t,
  int32_t,
  int32_t
);

struct _M0TUsiE* _M0FP26biolab8bio__seq17parse__label__rec(
  moonbit_string_t,
  int32_t,
  int32_t
);

int32_t _M0FP26biolab8bio__seq14is__whitespace(int32_t);

int32_t _M0MP26biolab8bio__seq8TreeNode11count__tips(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

int32_t _M0MP26biolab8bio__seq8TreeNode7is__tip(
  struct _M0TP26biolab8bio__seq8TreeNode*
);

moonbit_string_t _M0FP26biolab8bio__seq24gen__caterpillar__newick(int32_t);

moonbit_string_t _M0FP26biolab8bio__seq21gen__balanced__newick(int32_t);

int32_t _M0FP26biolab8bio__seq18gen__balanced__rec(
  struct _M0TPB13StringBuilder*,
  int32_t,
  struct _M0TPB5ArrayGiE*
);

int32_t _M0MPC15array5Array3setGiE(struct _M0TPB5ArrayGiE*, int32_t, int32_t);

struct _M0TP26biolab8bio__seq8TreeNode* _M0MPC15array5Array6removeGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*,
  int32_t
);

moonbit_string_t _M0MPC15array5Array2atGsE(struct _M0TPB5ArrayGsE*, int32_t);

struct _M0TP26biolab8bio__seq8TreeNode* _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*,
  int32_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*,
  int32_t
);

int32_t _M0MPC15array5Array2atGiE(struct _M0TPB5ArrayGiE*, int32_t);

int32_t _M0MPC15array5Array28unsafe__truncate__to__lengthGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*,
  int32_t
);

int32_t _M0IPB9SourceLocPB4Show6output(
  moonbit_string_t,
  struct _M0TPB6Logger
);

int32_t _M0FPB7printlnGsE(moonbit_string_t);

moonbit_string_t _M0MPC16double6Double10to__string(double);

moonbit_string_t _M0FPB15ryu__to__string(double);

struct _M0TPB17FloatingDecimal64* _M0FPB15d2d__small__int(uint64_t, int32_t);

moonbit_string_t _M0FPB9to__chars(struct _M0TPB17FloatingDecimal64*, int32_t);

struct _M0TPB17FloatingDecimal64* _M0FPB3d2d(uint64_t, uint32_t);

uint64_t _M0MPC14bool4Bool10to__uint64(int32_t);

int64_t _M0MPC14bool4Bool9to__int64(int32_t);

int32_t _M0MPC14bool4Bool7to__int(int32_t);

int32_t _M0FPB17decimal__length17(uint64_t);

struct _M0TPB8Pow5Pair _M0FPB22double__computeInvPow5(int32_t);

struct _M0TPB8Pow5Pair _M0FPB19double__computePow5(int32_t);

struct _M0TPB19MulShiftAll64Result _M0FPB13mulShiftAll64(
  uint64_t,
  struct _M0TPB8Pow5Pair,
  int32_t,
  int32_t
);

int32_t _M0FPB18multipleOfPowerOf2(uint64_t, int32_t);

int32_t _M0FPB18multipleOfPowerOf5(uint64_t, int32_t);

int32_t _M0FPB10pow5Factor(uint64_t);

uint64_t _M0FPB13shiftright128(uint64_t, uint64_t, int32_t);

struct _M0TPB7Umul128 _M0FPB7umul128(uint64_t, uint64_t);

moonbit_string_t _M0FPB19string__from__bytes(
  moonbit_bytes_t,
  int32_t,
  int32_t
);

int32_t _M0FPB9log10Pow2(int32_t);

int32_t _M0FPB9log10Pow5(int32_t);

moonbit_string_t _M0FPB18copy__special__str(int32_t, int32_t, int32_t);

int32_t _M0FPB8pow5bits(int32_t);

int32_t _M0IPC13int3IntPB4Hash4hash(int32_t);

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t);

struct _M0TUWEuQRPC15error5ErrorNsE* _M0MPB3Map3getGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  int32_t
);

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3getGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
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

int32_t _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

int32_t _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

int32_t _M0MPB3Map20rehash__place__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*
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

int32_t _M0MPC13int3Int3max(int32_t, int32_t);

int32_t _M0FPB21capacity__for__length(int32_t);

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0FPB8new__mapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  int32_t
);

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0FPB8new__mapGiUWEuQRPC15error5ErrorNsEE(
  int32_t
);

int32_t _M0MPC13int3Int20next__power__of__two(int32_t);

int32_t _M0FPB21calc__grow__threshold(int32_t);

double _M0MPC16option6Option10unwrap__orGdE(void*, double);

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE(
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*
);

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0MPC16option6Option6unwrapGRPB5EntryGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*
);

int32_t _M0IPC16option6OptionPB2Eq5equalGsE(
  moonbit_string_t,
  moonbit_string_t
);

struct _M0TPB4IterGsE* _M0MPC15array13ReadOnlyArray4iterGsE(
  moonbit_string_t*
);

uint64_t _M0MPC15array13ReadOnlyArray2atGmE(uint64_t*, int32_t);

uint32_t _M0MPC15array13ReadOnlyArray2atGjE(uint32_t*, int32_t);

struct _M0TPB4IterGsE* _M0MPC15array10FixedArray4iterGsE(moonbit_string_t*);

struct _M0TPB4IterGsE* _M0MPC15array9ArrayView4iterGsE(
  struct _M0TPB9ArrayViewGsE
);

moonbit_string_t _M0MPC15array9ArrayView4iterGsEC2338l729(struct _M0TWEOs*);

moonbit_string_t _M0IPC16uint646UInt64PB4Show10to__string(uint64_t);

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t);

moonbit_string_t _M0IPC14bool4BoolPB4Show10to__string(int32_t);

uint64_t _M0MPC14uint4UInt10to__uint64(uint32_t);

struct _M0TPC16string10StringView _M0MPC16string6String4trim(
  moonbit_string_t,
  void*
);

struct _M0TPC16string10StringView _M0MPC16string6String12trim_2einner(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string10StringView12trim_2einner(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

struct _M0TPC16string10StringView _M0MPC16string10StringView17trim__end_2einner(
  struct _M0TPC16string10StringView,
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

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*,
  struct _M0TP26biolab8bio__seq8TreeNode*
);

int32_t _M0MPC15array5Array4pushGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array7reallocGUsiEE(struct _M0TPB5ArrayGUsiEE*);

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

int32_t _M0MPC15array5Array7reallocGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*
);

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*,
  int32_t
);

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

int32_t _M0MPC15array5Array6lengthGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*
);

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(int32_t);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(
  int32_t
);

struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0MPC15array5Array11new_2einnerGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  int32_t
);

int32_t _M0MPC16string6String11has__suffix(
  moonbit_string_t,
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView11has__suffix(
  struct _M0TPC16string10StringView,
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

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(moonbit_string_t);

int64_t _M0FPB33boyer__moore__horspool__rev__find(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

struct moonbit_result_0 _M0FPB12assert__true(
  int32_t,
  void*,
  moonbit_string_t
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

struct _M0TP26biolab8bio__seq8TreeNode** _M0MPC15array5Array6bufferGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0MPC15array5Array6bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*
);

int32_t* _M0MPC15array5Array6bufferGiE(struct _M0TPB5ArrayGiE*);

struct _M0TPB4IterGsE* _M0MPB4Iter3newGsE(struct _M0TWEOs*, int64_t);

struct moonbit_result_0 _M0FPB10assert__eqGiE(
  int32_t,
  int32_t,
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_0 _M0FPB4failGuE(
  struct _M0TPC16string10StringView,
  moonbit_string_t
);

moonbit_string_t _M0FPB13debug__stringGiE(int32_t);

moonbit_string_t _M0MPC16uint646UInt6418to__string_2einner(uint64_t, int32_t);

int32_t _M0FPB22int64__to__string__dec(uint16_t*, uint64_t, int32_t, int32_t);

int32_t _M0FPB26int64__to__string__generic(
  uint16_t*,
  uint64_t,
  int32_t,
  int32_t,
  int32_t
);

int32_t _M0FPB22int64__to__string__hex(uint16_t*, uint64_t, int32_t, int32_t);

int32_t _M0FPB14radix__count64(uint64_t, int32_t);

int32_t _M0FPB12hex__count64(uint64_t);

int32_t _M0FPB12dec__count64(uint64_t);

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

int32_t _M0IP016_24default__implPB4Show6outputGmE(
  uint64_t,
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

int32_t _M0IP016_24default__implPB2Eq10not__equalGOsE(
  moonbit_string_t,
  moonbit_string_t
);

uint64_t _M0MPC13int3Int10to__uint64(int32_t);

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

struct _M0TP26biolab8bio__seq8TreeNode** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0MPB18UninitializedArray23make__and__blit_2einnerGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**,
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

int32_t _M0MPB13StringBuilder13write__objectGmE(
  struct _M0TPB13StringBuilder*,
  uint64_t
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

struct _M0TP26biolab8bio__seq8TreeNode** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode**,
  int32_t,
  int32_t,
  int32_t,
  int32_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0MPB18UninitializedArray23unsafe__make__and__blitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode**,
  int32_t,
  struct _M0TP26biolab8bio__seq8TreeNode**,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray12unsafe__blitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**,
  int32_t,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**,
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TP26biolab8bio__seq8TreeNode**,
  int32_t,
  struct _M0TP26biolab8bio__seq8TreeNode**,
  int32_t,
  int32_t
);

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**,
  int32_t,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**,
  int32_t,
  int32_t
);

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t*);

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(struct _M0TUsiE**);

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode**
);

int32_t _M0MPB18UninitializedArray6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**
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

struct _M0TP26biolab8bio__seq8TreeNode* _M0FPC15abort5abortGRP26biolab8bio__seq8TreeNodeE(
  moonbit_string_t
);

struct _M0TP26biolab8bio__seq8TreeNode** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq8TreeNodeEE(
  moonbit_string_t
);

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0FPC15abort5abortGRPB18UninitializedArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEEE(
  moonbit_string_t
);

int32_t _M0FPC15abort5abortGiE(moonbit_string_t);

int64_t _M0FPC15abort5abortGOiE(moonbit_string_t);

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

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_109 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 112, 104, 
    121, 108, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_96 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 98, 101, 
    110, 99, 104, 95, 99, 111, 117, 110, 116, 95, 116, 105, 112, 115, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_11 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 116, 105, 
    112, 52, 57, 56, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 115, 107, 
    105, 112, 112, 101, 100, 32, 116, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[1]; 
} const moonbit_string_literal_1 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 0, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_111 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 112, 100, 
    98, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[40]; 
} const moonbit_string_literal_74 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 39, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 80, 97, 
    114, 115, 101, 69, 114, 114, 111, 114, 46, 73, 110, 118, 97, 108, 
    105, 100, 78, 117, 109, 98, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_16 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 56, 52, 58, 51, 45, 49, 56, 52, 58, 50, 
    55, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[53]; 
} const moonbit_string_literal_81 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 52, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_12 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 50, 49, 56, 58, 51, 45, 50, 49, 56, 58, 51, 
    49, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_122 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 97, 108, 
    105, 103, 110, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_61 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 110, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[31]; 
} const moonbit_string_literal_59 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 30, 114, 97, 
    100, 105, 120, 32, 109, 117, 115, 116, 32, 98, 101, 32, 98, 101, 
    116, 119, 101, 101, 110, 32, 50, 32, 97, 110, 100, 32, 51, 54, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_50 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 73, 110, 
    102, 105, 110, 105, 116, 121, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_48 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 78, 97, 78, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_4 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 114, 101, 115, 117, 108, 116, 34, 
    44, 34, 102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_58 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 32, 70, 
    65, 73, 76, 69, 68, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[47]; 
} const moonbit_string_literal_34 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 46, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 55, 50, 58, 51, 45, 55, 50, 58, 51, 55, 64, 
    98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[47]; 
} const moonbit_string_literal_31 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 46, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 57, 56, 58, 51, 45, 57, 56, 58, 51, 56, 64, 
    98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_46 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[28]; 
} const moonbit_string_literal_36 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 27, 69, 120, 
    116, 114, 97, 32, 99, 104, 97, 114, 97, 99, 116, 101, 114, 115, 32, 
    97, 102, 116, 101, 114, 32, 114, 111, 111, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_120 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 98, 105, 
    111, 95, 115, 101, 113, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_115 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 102, 97, 
    115, 116, 113, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_70 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 44, 32, 
    108, 101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_53 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    108, 115, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_118 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 99, 111, 
    100, 111, 110, 95, 116, 97, 98, 108, 101, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_116 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 102, 97, 
    115, 116, 97, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_67 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 111, 
    117, 110, 100, 115, 32, 99, 104, 101, 99, 107, 32, 102, 97, 105, 
    108, 101, 100, 58, 32, 97, 108, 108, 111, 99, 97, 116, 101, 95, 108, 
    101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_87 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 98, 101, 
    110, 99, 104, 95, 112, 114, 101, 111, 114, 100, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_89 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 98, 101, 
    110, 99, 104, 95, 108, 101, 118, 101, 108, 111, 114, 100, 101, 114, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_20 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 116, 105, 
    112, 52, 48, 57, 53, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_10 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 116, 105, 
    112, 48, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_0 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 9, 10, 13, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_27 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 50, 50, 58, 51, 45, 49, 50, 50, 58, 51, 
    52, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_104 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 115, 101, 
    113, 117, 116, 105, 108, 115, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_86 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 98, 101, 
    110, 99, 104, 95, 116, 111, 95, 110, 101, 119, 105, 99, 107, 95, 
    108, 97, 114, 103, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_84 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 98, 101, 
    110, 99, 104, 95, 110, 101, 119, 105, 99, 107, 95, 112, 97, 114, 
    115, 101, 95, 108, 97, 114, 103, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_73 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 41, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_41 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 58, 48, 
    46, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_108 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 115, 101, 
    113, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_91 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 98, 101, 
    110, 99, 104, 95, 102, 105, 110, 100, 95, 108, 97, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_60 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 48, 49, 
    50, 51, 52, 53, 54, 55, 56, 57, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 119, 120, 121, 122, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[47]; 
} const moonbit_string_literal_35 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 46, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 55, 54, 58, 51, 45, 55, 54, 58, 51, 57, 64, 
    98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_107 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 115, 101, 
    113, 95, 114, 101, 99, 111, 114, 100, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_44 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 105, 110, 
    100, 101, 120, 32, 111, 117, 116, 32, 111, 102, 32, 98, 111, 117, 
    110, 100, 115, 58, 32, 116, 104, 101, 32, 108, 101, 110, 32, 105, 
    115, 32, 102, 114, 111, 109, 32, 48, 32, 116, 111, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_18 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 116, 105, 
    112, 50, 48, 48, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_15 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 57, 49, 58, 51, 45, 49, 57, 49, 58, 51, 
    55, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[20]; 
} const moonbit_string_literal_102 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 115, 107, 
    98, 105, 111, 95, 100, 105, 118, 101, 114, 115, 105, 116, 121, 46, 
    109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[51]; 
} const moonbit_string_literal_75 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 50, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 73, 110, 115, 112, 101, 
    99, 116, 69, 114, 114, 111, 114, 46, 73, 110, 115, 112, 101, 99, 
    116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_71 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 46, 108, 101, 110, 103, 116, 104, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_8 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 45, 45, 
    45, 45, 45, 32, 69, 78, 68, 32, 77, 79, 79, 78, 32, 84, 69, 83, 84, 
    32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_98 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 98, 101, 
    110, 99, 104, 95, 99, 97, 116, 101, 114, 112, 105, 108, 108, 97, 
    114, 95, 108, 99, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_68 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[47]; 
} const moonbit_string_literal_32 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 46, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 57, 48, 58, 51, 45, 57, 48, 58, 51, 55, 64, 
    98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_9 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 123, 34, 
    116, 121, 112, 101, 34, 58, 34, 115, 116, 97, 114, 116, 34, 44, 34, 
    102, 105, 108, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_119 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 99, 108, 
    117, 115, 116, 97, 108, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_21 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 54, 56, 58, 51, 45, 49, 54, 56, 58, 51, 
    53, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_14 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 57, 57, 58, 51, 45, 49, 57, 57, 58, 51, 
    50, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_82 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 98, 101, 
    110, 99, 104, 95, 115, 101, 116, 117, 112, 95, 99, 104, 101, 99, 
    107, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_42 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 116, 105, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_5 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 34, 
    105, 110, 100, 101, 120, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[81]; 
} const moonbit_string_literal_76 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 80, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 77, 111, 
    111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 101, 
    114, 73, 110, 116, 101, 114, 110, 97, 108, 74, 115, 69, 114, 114, 
    111, 114, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 
    68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 
    74, 115, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[26]; 
} const moonbit_string_literal_47 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 25, 73, 108, 
    108, 101, 103, 97, 108, 65, 114, 103, 117, 109, 101, 110, 116, 69, 
    120, 99, 101, 112, 116, 105, 111, 110, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_45 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 32, 98, 
    117, 116, 32, 116, 104, 101, 32, 105, 110, 100, 101, 120, 32, 105, 
    115, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_72 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 70, 97, 
    105, 108, 117, 114, 101, 40, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_117 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 99, 111, 
    109, 112, 108, 101, 109, 101, 110, 116, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_7 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 125, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_112 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 112, 100, 
    98, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_97 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 98, 101, 
    110, 99, 104, 95, 99, 97, 116, 101, 114, 112, 105, 108, 108, 97, 
    114, 95, 102, 105, 110, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_90 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 98, 101, 
    110, 99, 104, 95, 102, 105, 110, 100, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[47]; 
} const moonbit_string_literal_33 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 46, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 56, 51, 58, 51, 45, 56, 51, 58, 51, 54, 64, 
    98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_30 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 48, 54, 58, 51, 45, 49, 48, 54, 58, 51, 
    56, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[48]; 
} const moonbit_string_literal_80 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 47, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 80, 97, 
    114, 115, 101, 69, 114, 114, 111, 114, 46, 85, 110, 98, 97, 108, 
    97, 110, 99, 101, 100, 80, 97, 114, 101, 110, 116, 104, 101, 115, 
    101, 115, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[35]; 
} const moonbit_string_literal_3 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 34, 45, 45, 
    45, 45, 45, 32, 66, 69, 71, 73, 78, 32, 77, 79, 79, 78, 32, 84, 69, 
    83, 84, 32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[27]; 
} const moonbit_string_literal_99 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 98, 101, 
    110, 99, 104, 95, 99, 97, 116, 101, 114, 112, 105, 108, 108, 97, 
    114, 95, 100, 105, 115, 116, 97, 110, 99, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_64 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 116, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_62 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 114, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_55 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 96, 32, 
    105, 115, 32, 110, 111, 116, 32, 116, 114, 117, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_40 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 58, 48, 
    46, 49, 44, 116, 105, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_69 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    100, 115, 116, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_66 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 105, 110, 
    118, 97, 108, 105, 100, 32, 99, 111, 100, 101, 32, 112, 111, 105, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_49 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 45, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_23 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 52, 54, 58, 51, 45, 49, 52, 54, 58, 51, 
    50, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_6 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 44, 34, 
    109, 101, 115, 115, 97, 103, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_19 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 55, 54, 58, 51, 45, 49, 55, 54, 58, 51, 
    49, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_110 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 112, 104, 
    121, 108, 105, 112, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_106 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 115, 101, 
    113, 102, 101, 97, 116, 117, 114, 101, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[20]; 
} const moonbit_string_literal_103 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 115, 107, 
    98, 105, 111, 95, 97, 108, 105, 103, 110, 109, 101, 110, 116, 46, 
    109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_95 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 98, 101, 
    110, 99, 104, 95, 116, 111, 116, 97, 108, 95, 108, 101, 110, 103, 
    116, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[25]; 
} const moonbit_string_literal_83 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 24, 98, 101, 
    110, 99, 104, 95, 110, 101, 119, 105, 99, 107, 95, 112, 97, 114, 
    115, 101, 95, 115, 109, 97, 108, 108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_25 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 51, 56, 58, 51, 45, 49, 51, 56, 58, 51, 
    50, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_54 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 96, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_92 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 98, 101, 
    110, 99, 104, 95, 108, 99, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[40]; 
} const moonbit_string_literal_79 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 39, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 80, 97, 
    114, 115, 101, 69, 114, 114, 111, 114, 46, 73, 110, 118, 97, 108, 
    105, 100, 70, 111, 114, 109, 97, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_65 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 92, 117, 123, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_94 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 98, 101, 
    110, 99, 104, 95, 100, 105, 115, 116, 97, 110, 99, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_85 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 98, 101, 
    110, 99, 104, 95, 116, 111, 95, 110, 101, 119, 105, 99, 107, 95, 
    115, 109, 97, 108, 108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_52 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 116, 114, 
    117, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_93 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 98, 101, 
    110, 99, 104, 95, 108, 99, 97, 95, 100, 101, 101, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[20]; 
} const moonbit_string_literal_56 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 73, 110, 
    118, 97, 108, 105, 100, 32, 115, 116, 97, 114, 116, 32, 105, 110, 
    100, 101, 120, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_37 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 116, 105, 
    112, 48, 58, 48, 46, 49, 59, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_113 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 110, 101, 
    119, 105, 99, 107, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_101 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_100 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_63 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_22 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 53, 55, 58, 51, 45, 49, 53, 55, 58, 51, 
    53, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_13 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 50, 49, 48, 58, 51, 45, 50, 49, 48, 58, 51, 
    53, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_121 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 97, 108, 
    105, 103, 110, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_88 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 98, 101, 
    110, 99, 104, 95, 112, 111, 115, 116, 111, 114, 100, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_38 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 116, 105, 
    112, 48, 58, 48, 46, 49, 44, 116, 105, 112, 49, 58, 48, 46, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_26 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 51, 48, 58, 51, 45, 49, 51, 48, 58, 51, 
    52, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[43]; 
} const moonbit_string_literal_78 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 42, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 80, 97, 
    114, 115, 101, 69, 114, 114, 111, 114, 46, 77, 105, 115, 115, 105, 
    110, 103, 83, 101, 109, 105, 99, 111, 108, 111, 110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_24 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 116, 105, 
    112, 49, 48, 48, 48, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_39 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 105, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_17 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 116, 105, 
    112, 49, 48, 48, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[83]; 
} const moonbit_string_literal_77 =
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
} const moonbit_string_literal_57 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 32, 33, 
    61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_43 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 48, 46, 49, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[49]; 
} const moonbit_string_literal_28 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 48, 115, 107, 
    98, 105, 111, 95, 116, 114, 101, 101, 95, 98, 101, 110, 99, 104, 
    46, 109, 98, 116, 58, 49, 49, 52, 58, 51, 45, 49, 49, 52, 58, 51, 
    52, 64, 98, 105, 111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 
    113, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_114 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 103, 101, 
    110, 98, 97, 110, 107, 95, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_105 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 115, 101, 
    113, 105, 111, 46, 109, 98, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_51 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 48, 46, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_29 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 59, 0};

struct moonbit_object const moonbit_constant_constructor_0 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0)
  };

struct moonbit_object const moonbit_constant_constructor_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 2)
  };

struct moonbit_object const moonbit_constant_constructor_3 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 3)
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__0_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__0_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__8_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__8_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__4_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__4_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__7_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__7_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__9_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__9_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__12_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__12_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__1_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__1_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__15_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__15_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__10_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__10_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__14_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__14_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__13_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__13_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__3_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__3_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__2_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__2_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__17_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__17_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__6_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__6_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWRPC15error5ErrorEs data; 
} const _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1545$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1545
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__11_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__11_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__16_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__16_2edyncall
  };

struct { int32_t rc; uint32_t meta; struct _M0TWEuQRPC15error5Error data; 
} const _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__5_2edyncall$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__5_2edyncall
  };

uint32_t const moonbit_layout_table_data[119] =
  {
    sizeof(struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533)
    / 4, 1,
    offsetof(struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533, $1)
    / 4,
    sizeof(struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538)
    / 4, 1,
    offsetof(struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538, $1)
    / 4,
    sizeof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest, $0)
    / 4,
    sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__)
    / 4, 2,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__, $0)
    / 4,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__, $1)
    / 4,
    sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__)
    / 4, 3,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__, $0)
    / 4,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__, $1)
    / 4,
    offsetof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__, $2)
    / 4, sizeof(struct _M0TPB5ArrayGUsiEE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGUsiEE, $0) / 4, sizeof(struct _M0TUsiE) / 4,
    1, offsetof(struct _M0TUsiE, $0) / 4, sizeof(struct _M0TPB5ArrayGsE) / 4,
    1, offsetof(struct _M0TPB5ArrayGsE, $0) / 4,
    sizeof(struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidFormat)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidFormat, $0)
    / 4, sizeof(struct _M0TPB8MutLocalGOsE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGOsE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGOdE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGOdE, $0) / 4,
    sizeof(struct _M0DTPC16option6OptionGOsE4Some) / 4, 1,
    offsetof(struct _M0DTPC16option6OptionGOsE4Some, $0) / 4,
    sizeof(struct _M0TURP26biolab8bio__seq8TreeNodeiE) / 4, 1,
    offsetof(struct _M0TURP26biolab8bio__seq8TreeNodeiE, $0) / 4,
    sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE, $0) / 4,
    sizeof(struct _M0TP26biolab8bio__seq8TreeNode) / 4, 4,
    offsetof(struct _M0TP26biolab8bio__seq8TreeNode, $0) / 4,
    offsetof(struct _M0TP26biolab8bio__seq8TreeNode, $1) / 4,
    offsetof(struct _M0TP26biolab8bio__seq8TreeNode, $2) / 4,
    offsetof(struct _M0TP26biolab8bio__seq8TreeNode, $3) / 4,
    sizeof(struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidNumber)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidNumber, $0)
    / 4, sizeof(struct _M0TPB5ArrayGiE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGiE, $0) / 4,
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
    sizeof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE) / 4, 
    2,
    offsetof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE, $5) / 4,
    sizeof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE) / 4, 2,
    offsetof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE, $0) / 4,
    offsetof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE, $5) / 4,
    sizeof(struct _M0TPB9ArrayViewGsE) / 4, 1,
    offsetof(struct _M0TPB9ArrayViewGsE, $0) / 4,
    sizeof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__)
    / 4, 2,
    (offsetof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__, $0)
     + offsetof(struct _M0TPB9ArrayViewGsE, $0))
    / 4,
    offsetof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__, $2)
    / 4,
    sizeof(struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE) / 4,
    1,
    offsetof(struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE, $0)
    / 4, sizeof(struct _M0TPB4IterGsE) / 4, 1,
    offsetof(struct _M0TPB4IterGsE, $0) / 4,
    sizeof(struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure, $0)
    / 4, sizeof(struct _M0TPC16string10StringView) / 4, 1,
    offsetof(struct _M0TPC16string10StringView, $0) / 4,
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

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__13_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__13_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__4_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__4_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__16_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__16_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__14_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__14_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__8_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__8_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__1_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__1_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__9_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__9_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__6_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__6_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__11_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__11_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__3_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__3_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__2_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__2_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__10_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__10_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__17_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__17_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__7_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__7_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__15_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__15_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__0_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__0_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__12_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__12_2edyncall$closure.data;

struct _M0TWEuQRPC15error5Error* _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__5_2eclo =
  (struct _M0TWEuQRPC15error5Error*)&_M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__5_2edyncall$closure.data;

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

moonbit_string_t _M0MPC16string6String4trimN7_2abindS7043 =
  (moonbit_string_t)moonbit_string_literal_0.data;

struct { int32_t rc; uint32_t meta; uint64_t data[30]; 
} _M0FPB26gDOUBLE__POW5__INV__SPLIT2$object =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 30, 1ull,
    2305843009213693952ull, 5955668970331000884ull, 1784059615882449851ull,
    8982663654677661702ull, 1380349269358112757ull, 7286864317269821294ull,
    2135987035920910082ull, 7005857020398200553ull, 1652639921975621497ull,
    17965325103354776697ull, 1278668206209430417ull, 8928596168509315048ull,
    1978643211784836272ull, 10075671573058298858ull, 1530901034580419511ull,
    597001226353042382ull, 1184477304306571148ull, 1527430471115325346ull,
    1832889850782397517ull, 12533209867169019542ull, 1418129833677084982ull,
    5577825024675947042ull, 2194449627517475473ull, 11006974540203867551ull,
    1697873161311732311ull, 10313493231639821582ull, 1313665730009899186ull,
    12701016819766672773ull, 2032799256770390445ull
  };

uint64_t* _M0FPB26gDOUBLE__POW5__INV__SPLIT2 =
  _M0FPB26gDOUBLE__POW5__INV__SPLIT2$object.data;

struct { int32_t rc; uint32_t meta; uint32_t data[19]; 
} _M0FPB19gPOW5__INV__OFFSETS$object =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 1414808916u,
    67458373u, 268701696u, 4195348u, 1073807360u, 1091917141u, 1108u, 
    65604u, 1073741824u, 1140850753u, 1346716752u, 1431634004u, 1365595476u,
    1073758208u, 16777217u, 66816u, 1364284433u, 89478484u, 0u
  };

uint32_t* _M0FPB19gPOW5__INV__OFFSETS =
  _M0FPB19gPOW5__INV__OFFSETS$object.data;

struct { int32_t rc; uint32_t meta; uint64_t data[26]; 
} _M0FPB21gDOUBLE__POW5__SPLIT2$object =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 0ull,
    1152921504606846976ull, 0ull, 1490116119384765625ull,
    1032610780636961552ull, 1925929944387235853ull, 7910200175544436838ull,
    1244603055572228341ull, 16941905809032713930ull, 1608611746708759036ull,
    13024893955298202172ull, 2079081953128979843ull, 6607496772837067824ull,
    1343575221513417750ull, 17332926989895652603ull, 1736530273035216783ull,
    13037379183483547984ull, 2244412773384604712ull, 1605989338741628675ull,
    1450417759929778918ull, 9630225068416591280ull, 1874621017369538693ull,
    665883850346957067ull, 1211445438634777304ull, 14931890668723713708ull,
    1565756531257009982ull
  };

uint64_t* _M0FPB21gDOUBLE__POW5__SPLIT2 =
  _M0FPB21gDOUBLE__POW5__SPLIT2$object.data;

struct { int32_t rc; uint32_t meta; uint32_t data[21]; 
} _M0FPB14gPOW5__OFFSETS$object =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 0u, 0u, 
    0u, 0u, 1073741824u, 1500076437u, 1431590229u, 1448432917u, 1091896580u,
    1079333904u, 1146442053u, 1146111296u, 1163220304u, 1073758208u,
    2521039936u, 1431721317u, 1413824581u, 1075134801u, 1431671125u,
    1363170645u, 261u
  };

uint32_t* _M0FPB14gPOW5__OFFSETS = _M0FPB14gPOW5__OFFSETS$object.data;

struct { int32_t rc; uint32_t meta; uint64_t data[26]; 
} _M0FPB20gDOUBLE__POW5__TABLE$object =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 26, 1ull, 5ull,
    25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull,
    1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull,
    6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull,
    3814697265625ull, 19073486328125ull, 95367431640625ull,
    476837158203125ull, 2384185791015625ull, 11920928955078125ull,
    59604644775390625ull, 298023223876953125ull
  };

uint64_t* _M0FPB20gDOUBLE__POW5__TABLE =
  _M0FPB20gDOUBLE__POW5__TABLE$object.data;

int64_t _M0MPB4Iter4nextN6constrS9980GsE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GsE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GsE = 0ll;

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0FP26biolab8bio__seq48moonbit__test__driver__internal__no__args__tests;

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__5_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3498
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__5();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__12_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3497
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__12();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__0_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3496
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__0();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__15_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3495
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__15();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__7_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3494
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__7();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__17_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3493
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__17();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__10_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3492
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__10();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__2_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3491
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__2();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__3_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3490
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__3();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__11_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3489
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__11();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__6_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3488
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__6();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__9_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3487
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__9();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__1_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3486
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__1();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__8_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3485
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__8();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__14_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3484
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__14();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__16_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3483
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__16();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq63____test__736b62696f5f747265655f62656e63682e6d6274__4_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3482
) {
  return _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__4();
}

struct moonbit_result_0 _M0FP26biolab8bio__seq64____test__736b62696f5f747265655f62656e63682e6d6274__13_2edyncall(
  struct _M0TWEuQRPC15error5Error* _M0L6_2aenvS3481
) {
  return _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__13();
}

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1566,
  moonbit_string_t _M0L8filenameS1535,
  int32_t _M0L5indexS1537
) {
  struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533* _closure_3746;
  struct _M0TWEu* _M0L13handle__startS1533;
  struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538* _closure_3747;
  struct _M0TWssbEu* _M0L14handle__resultS1538;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1545;
  void* _M0L11_2atry__errS1560;
  struct moonbit_result_0 _tmp_3749;
  int32_t _handle__error__result_3750;
  int32_t _M0L6_2atmpS3469;
  void* _M0L3errS1561;
  moonbit_string_t _M0L4nameS1563;
  struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest* _M0L36_2aMoonBitTestDriverInternalSkipTestS1564;
  moonbit_string_t _M0L8_2afieldS3499;
  int32_t _M0L6_2acntS3700;
  moonbit_string_t _M0L7_2anameS1565;
  #line 619 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  moonbit_incref(_M0L8filenameS1535);
  _closure_3746
  = (struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533*)moonbit_malloc(sizeof(struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533));
  Moonbit_object_header(_closure_3746)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 0, 0);
  _closure_3746->code
  = &_M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN13handle__startS1533;
  _closure_3746->$0 = _M0L5indexS1537;
  _closure_3746->$1 = _M0L8filenameS1535;
  _M0L13handle__startS1533 = (struct _M0TWEu*)_closure_3746;
  moonbit_incref(_M0L8filenameS1535);
  _closure_3747
  = (struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538*)moonbit_malloc(sizeof(struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538));
  Moonbit_object_header(_closure_3747)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 3, 0);
  _closure_3747->code
  = &_M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1538;
  _closure_3747->$0 = _M0L5indexS1537;
  _closure_3747->$1 = _M0L8filenameS1535;
  _M0L14handle__resultS1538 = (struct _M0TWssbEu*)_closure_3747;
  _M0L17error__to__stringS1545
  = (struct _M0TWRPC15error5ErrorEs*)&_M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1545$closure.data;
  #line 661 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _tmp_3749
  = _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__test(_M0L12async__testsS1566, _M0L8filenameS1535, _M0L5indexS1537, _M0L13handle__startS1533, _M0L14handle__resultS1538, _M0L17error__to__stringS1545);
  if (_tmp_3749.tag) {
    int32_t const _M0L5_2aokS3478 = _tmp_3749.data.ok;
    _handle__error__result_3750 = _M0L5_2aokS3478;
  } else {
    void* const _M0L6_2aerrS3479 = _tmp_3749.data.err;
    moonbit_decref(_M0L17error__to__stringS1545);
    moonbit_decref(_M0L13handle__startS1533);
    _M0L11_2atry__errS1560 = _M0L6_2aerrS3479;
    goto join_1559;
  }
  if (_handle__error__result_3750) {
    moonbit_decref(_M0L17error__to__stringS1545);
    moonbit_decref(_M0L13handle__startS1533);
    _M0L6_2atmpS3469 = 1;
  } else {
    struct moonbit_result_0 _tmp_3751;
    int32_t _handle__error__result_3752;
    #line 664 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
    _tmp_3751
    = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq43MoonBit__Test__Driver__Internal__With__ArgsE(_M0L12async__testsS1566, _M0L8filenameS1535, _M0L5indexS1537, _M0L13handle__startS1533, _M0L14handle__resultS1538, _M0L17error__to__stringS1545);
    if (_tmp_3751.tag) {
      int32_t const _M0L5_2aokS3476 = _tmp_3751.data.ok;
      _handle__error__result_3752 = _M0L5_2aokS3476;
    } else {
      void* const _M0L6_2aerrS3477 = _tmp_3751.data.err;
      moonbit_decref(_M0L17error__to__stringS1545);
      moonbit_decref(_M0L13handle__startS1533);
      _M0L11_2atry__errS1560 = _M0L6_2aerrS3477;
      goto join_1559;
    }
    if (_handle__error__result_3752) {
      moonbit_decref(_M0L17error__to__stringS1545);
      moonbit_decref(_M0L13handle__startS1533);
      _M0L6_2atmpS3469 = 1;
    } else {
      struct moonbit_result_0 _tmp_3753;
      int32_t _handle__error__result_3754;
      #line 667 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _tmp_3753
      = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq48MoonBit__Test__Driver__Internal__Async__No__ArgsE(_M0L12async__testsS1566, _M0L8filenameS1535, _M0L5indexS1537, _M0L13handle__startS1533, _M0L14handle__resultS1538, _M0L17error__to__stringS1545);
      if (_tmp_3753.tag) {
        int32_t const _M0L5_2aokS3474 = _tmp_3753.data.ok;
        _handle__error__result_3754 = _M0L5_2aokS3474;
      } else {
        void* const _M0L6_2aerrS3475 = _tmp_3753.data.err;
        moonbit_decref(_M0L17error__to__stringS1545);
        moonbit_decref(_M0L13handle__startS1533);
        _M0L11_2atry__errS1560 = _M0L6_2aerrS3475;
        goto join_1559;
      }
      if (_handle__error__result_3754) {
        moonbit_decref(_M0L17error__to__stringS1545);
        moonbit_decref(_M0L13handle__startS1533);
        _M0L6_2atmpS3469 = 1;
      } else {
        struct moonbit_result_0 _tmp_3755;
        int32_t _handle__error__result_3756;
        #line 670 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        _tmp_3755
        = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__Async__With__ArgsE(_M0L12async__testsS1566, _M0L8filenameS1535, _M0L5indexS1537, _M0L13handle__startS1533, _M0L14handle__resultS1538, _M0L17error__to__stringS1545);
        if (_tmp_3755.tag) {
          int32_t const _M0L5_2aokS3472 = _tmp_3755.data.ok;
          _handle__error__result_3756 = _M0L5_2aokS3472;
        } else {
          void* const _M0L6_2aerrS3473 = _tmp_3755.data.err;
          moonbit_decref(_M0L17error__to__stringS1545);
          moonbit_decref(_M0L13handle__startS1533);
          _M0L11_2atry__errS1560 = _M0L6_2aerrS3473;
          goto join_1559;
        }
        if (_handle__error__result_3756) {
          moonbit_decref(_M0L17error__to__stringS1545);
          moonbit_decref(_M0L13handle__startS1533);
          _M0L6_2atmpS3469 = 1;
        } else {
          struct moonbit_result_0 _tmp_3757;
          #line 673 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
          _tmp_3757
          = _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(_M0L12async__testsS1566, _M0L8filenameS1535, _M0L5indexS1537, _M0L13handle__startS1533, _M0L14handle__resultS1538, _M0L17error__to__stringS1545);
          moonbit_decref(_M0L13handle__startS1533);
          moonbit_decref(_M0L17error__to__stringS1545);
          if (_tmp_3757.tag) {
            int32_t const _M0L5_2aokS3470 = _tmp_3757.data.ok;
            _M0L6_2atmpS3469 = _M0L5_2aokS3470;
          } else {
            void* const _M0L6_2aerrS3471 = _tmp_3757.data.err;
            _M0L11_2atry__errS1560 = _M0L6_2aerrS3471;
            goto join_1559;
          }
        }
      }
    }
  }
  if (!_M0L6_2atmpS3469) {
    void* _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3480 =
      (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
    Moonbit_object_header(_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3480)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 5);
    ((struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3480)->$0
    = (moonbit_string_t)moonbit_string_literal_1.data;
    _M0L11_2atry__errS1560
    = _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3480;
    goto join_1559;
  } else {
    moonbit_decref(_M0L14handle__resultS1538);
  }
  goto joinlet_3748;
  join_1559:;
  _M0L3errS1561 = _M0L11_2atry__errS1560;
  _M0L36_2aMoonBitTestDriverInternalSkipTestS1564
  = (struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L3errS1561;
  _M0L8_2afieldS3499 = _M0L36_2aMoonBitTestDriverInternalSkipTestS1564->$0;
  _M0L6_2acntS3700
  = Moonbit_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1564));
  if (_M0L6_2acntS3700 > 1) {
    int32_t _M0L11_2anew__cntS3701 = _M0L6_2acntS3700 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1564), _M0L11_2anew__cntS3701);
    moonbit_incref(_M0L8_2afieldS3499);
  } else if (_M0L6_2acntS3700 == 1) {
    #line 680 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
    moonbit_free(_M0L36_2aMoonBitTestDriverInternalSkipTestS1564);
  }
  _M0L7_2anameS1565 = _M0L8_2afieldS3499;
  _M0L4nameS1563 = _M0L7_2anameS1565;
  goto join_1562;
  goto joinlet_3758;
  join_1562:;
  #line 681 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1538(_M0L14handle__resultS1538, _M0L4nameS1563, (moonbit_string_t)moonbit_string_literal_2.data, 1);
  moonbit_decref(_M0L14handle__resultS1538);
  moonbit_decref(_M0L4nameS1563);
  joinlet_3758:;
  joinlet_3748:;
  return 0;
}

moonbit_string_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN17error__to__stringS1545(
  struct _M0TWRPC15error5ErrorEs* _M0L6_2aenvS3468,
  void* _M0L3errS1546
) {
  void* _M0L1eS1548;
  moonbit_string_t _M0L1eS1550;
  moonbit_string_t _result_3761;
  #line 650 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  switch (Moonbit_object_tag(_M0L3errS1546)) {
    case 0: {
      struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS1551 =
        (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L3errS1546;
      moonbit_string_t _M0L4_2aeS1552 = _M0L10_2aFailureS1551->$0;
      moonbit_incref(_M0L4_2aeS1552);
      _M0L1eS1550 = _M0L4_2aeS1552;
      goto join_1549;
      break;
    }
    
    case 6: {
      struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError* _M0L15_2aInspectErrorS1553 =
        (struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError*)_M0L3errS1546;
      moonbit_string_t _M0L4_2aeS1554 = _M0L15_2aInspectErrorS1553->$0;
      moonbit_incref(_M0L4_2aeS1554);
      _M0L1eS1550 = _M0L4_2aeS1554;
      goto join_1549;
      break;
    }
    
    case 7: {
      struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError* _M0L16_2aSnapshotErrorS1555 =
        (struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError*)_M0L3errS1546;
      moonbit_string_t _M0L4_2aeS1556 = _M0L16_2aSnapshotErrorS1555->$0;
      moonbit_incref(_M0L4_2aeS1556);
      _M0L1eS1550 = _M0L4_2aeS1556;
      goto join_1549;
      break;
    }
    
    case 8: {
      struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError* _M0L35_2aMoonBitTestDriverInternalJsErrorS1557 =
        (struct _M0DTPC15error5Error87biolab_2fbio__seq_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError*)_M0L3errS1546;
      moonbit_string_t _M0L4_2aeS1558 =
        _M0L35_2aMoonBitTestDriverInternalJsErrorS1557->$0;
      moonbit_incref(_M0L4_2aeS1558);
      _M0L1eS1550 = _M0L4_2aeS1558;
      goto join_1549;
      break;
    }
    default: {
      moonbit_incref(_M0L3errS1546);
      _M0L1eS1548 = _M0L3errS1546;
      goto join_1547;
      break;
    }
  }
  join_1549:;
  return _M0L1eS1550;
  join_1547:;
  #line 656 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _result_3761 = _M0FP15Error10to__string(_M0L1eS1548);
  moonbit_decref(_M0L1eS1548);
  return _result_3761;
}

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN14handle__resultS1538(
  struct _M0TWssbEu* _M0L6_2aenvS3465,
  moonbit_string_t _M0L10__testnameS1539,
  moonbit_string_t _M0L7messageS1540,
  int32_t _M0L7skippedS1541
) {
  struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538* _M0L14_2acasted__envS3466;
  moonbit_string_t _M0L8filenameS1535;
  int32_t _M0L5indexS1537;
  int32_t _if__result_3762;
  moonbit_string_t _M0L10file__nameS1542;
  moonbit_string_t _M0L7messageS1543;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1544;
  moonbit_string_t _M0L6_2atmpS3467;
  #line 635 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L14_2acasted__envS3466
  = (struct _M0R91_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1538*)_M0L6_2aenvS3465;
  _M0L8filenameS1535 = _M0L14_2acasted__envS3466->$1;
  _M0L5indexS1537 = _M0L14_2acasted__envS3466->$0;
  if (!_M0L7skippedS1541) {
    _if__result_3762 = 1;
  } else {
    _if__result_3762 = 0;
  }
  if (_if__result_3762) {
    
  }
  #line 641 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L10file__nameS1542
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1535, 1);
  #line 642 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L7messageS1543
  = _M0MPC16string6String14escape_2einner(_M0L7messageS1540, 1);
  #line 643 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_3.data);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L18_2astring__builderS1544
  = _M0MPB13StringBuilder21StringBuilder_2einner(45);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1544, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1544, _M0L10file__nameS1542);
  moonbit_decref(_M0L10file__nameS1542);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1544, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1544, _M0L5indexS1537);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1544, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1544, _M0L7messageS1543);
  moonbit_decref(_M0L7messageS1543);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1544, (moonbit_string_t)moonbit_string_literal_7.data);
  #line 645 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS3467
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1544);
  moonbit_decref(_M0L18_2astring__builderS1544);
  #line 644 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS3467);
  moonbit_decref(_M0L6_2atmpS3467);
  #line 647 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_8.data);
  return 0;
}

int32_t _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__executeN13handle__startS1533(
  struct _M0TWEu* _M0L6_2aenvS3462
) {
  struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533* _M0L14_2acasted__envS3463;
  moonbit_string_t _M0L8filenameS1535;
  int32_t _M0L5indexS1537;
  moonbit_string_t _M0L10file__nameS1534;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1536;
  moonbit_string_t _M0L6_2atmpS3464;
  #line 626 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L14_2acasted__envS3463
  = (struct _M0R90_24biolab_2fbio__seq_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1533*)_M0L6_2aenvS3462;
  _M0L8filenameS1535 = _M0L14_2acasted__envS3463->$1;
  _M0L5indexS1537 = _M0L14_2acasted__envS3463->$0;
  #line 627 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L10file__nameS1534
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1535, 1);
  #line 628 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_3.data);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L18_2astring__builderS1536
  = _M0MPB13StringBuilder21StringBuilder_2einner(33);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1536, (moonbit_string_t)moonbit_string_literal_9.data);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1536, _M0L10file__nameS1534);
  moonbit_decref(_M0L10file__nameS1534);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1536, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1536, _M0L5indexS1537);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1536, (moonbit_string_t)moonbit_string_literal_7.data);
  #line 630 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS3464
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1536);
  moonbit_decref(_M0L18_2astring__builderS1536);
  #line 629 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS3464);
  moonbit_decref(_M0L6_2atmpS3464);
  #line 632 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_8.data);
  return 0;
}

struct moonbit_result_0 _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__test(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1532,
  moonbit_string_t _M0L8filenameS1529,
  int32_t _M0L5indexS1523,
  struct _M0TWEu* _M0L13handle__startS1518,
  struct _M0TWssbEu* _M0L14handle__resultS1519,
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1521
) {
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L10index__mapS1498;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS1528;
  struct _M0TWEuQRPC15error5Error* _M0L1fS1500;
  moonbit_string_t* _M0L5attrsS1501;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L7_2abindS1522;
  moonbit_string_t _M0L4nameS1504;
  moonbit_string_t _M0L4nameS1502;
  int32_t _M0L6_2atmpS3461;
  struct _M0TPB4IterGsE* _M0L5_2aitS1506;
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__* _closure_3771;
  struct _M0TWEu* _M0L6_2atmpS3452;
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__* _closure_3772;
  struct _M0TWRPC15error5ErrorEu* _M0L6_2atmpS3453;
  struct moonbit_result_0 _result_3773;
  #line 468 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  #line 476 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L7_2abindS1528
  = _M0MPB3Map3getGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0FP26biolab8bio__seq48moonbit__test__driver__internal__no__args__tests, _M0L8filenameS1529);
  if (_M0L7_2abindS1528 == 0) {
    struct moonbit_result_0 _result_3764;
    if (_M0L7_2abindS1528) {
      moonbit_decref(_M0L7_2abindS1528);
    }
    _result_3764.tag = 1;
    _result_3764.data.ok = 0;
    return _result_3764;
  } else {
    struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS1530 =
      _M0L7_2abindS1528;
    struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L13_2aindex__mapS1531 =
      _M0L7_2aSomeS1530;
    _M0L10index__mapS1498 = _M0L13_2aindex__mapS1531;
    goto join_1497;
  }
  join_1497:;
  #line 478 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L7_2abindS1522
  = _M0MPB3Map3getGiUWEuQRPC15error5ErrorNsEE(_M0L10index__mapS1498, _M0L5indexS1523);
  moonbit_decref(_M0L10index__mapS1498);
  if (_M0L7_2abindS1522 == 0) {
    struct moonbit_result_0 _result_3766;
    if (_M0L7_2abindS1522) {
      moonbit_decref(_M0L7_2abindS1522);
    }
    _result_3766.tag = 1;
    _result_3766.data.ok = 0;
    return _result_3766;
  } else {
    struct _M0TUWEuQRPC15error5ErrorNsE* _M0L7_2aSomeS1524 =
      _M0L7_2abindS1522;
    struct _M0TUWEuQRPC15error5ErrorNsE* _M0L4_2axS1525 = _M0L7_2aSomeS1524;
    struct _M0TWEuQRPC15error5Error* _M0L4_2afS1526 = _M0L4_2axS1525->$0;
    moonbit_string_t* _M0L8_2afieldS3501 = _M0L4_2axS1525->$1;
    int32_t _M0L6_2acntS3702 =
      Moonbit_rc_count(Moonbit_object_header(_M0L4_2axS1525));
    moonbit_string_t* _M0L8_2aattrsS1527;
    if (_M0L6_2acntS3702 > 1) {
      int32_t _M0L11_2anew__cntS3703 = _M0L6_2acntS3702 - 1;
      Moonbit_set_rc_count(Moonbit_object_header(_M0L4_2axS1525), _M0L11_2anew__cntS3703);
      moonbit_incref(_M0L8_2afieldS3501);
      moonbit_incref(_M0L4_2afS1526);
    } else if (_M0L6_2acntS3702 == 1) {
      #line 476 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      moonbit_free(_M0L4_2axS1525);
    }
    _M0L8_2aattrsS1527 = _M0L8_2afieldS3501;
    _M0L1fS1500 = _M0L4_2afS1526;
    _M0L5attrsS1501 = _M0L8_2aattrsS1527;
    goto join_1499;
  }
  join_1499:;
  _M0L6_2atmpS3461 = Moonbit_array_length(_M0L5attrsS1501);
  if (_M0L6_2atmpS3461 >= 1) {
    moonbit_string_t _M0L7_2anameS1505 = (moonbit_string_t)_M0L5attrsS1501[0];
    moonbit_incref(_M0L7_2anameS1505);
    _M0L4nameS1504 = _M0L7_2anameS1505;
    goto join_1503;
  } else {
    _M0L4nameS1502 = (moonbit_string_t)moonbit_string_literal_1.data;
  }
  goto joinlet_3767;
  join_1503:;
  _M0L4nameS1502 = _M0L4nameS1504;
  joinlet_3767:;
  #line 479 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L5_2aitS1506 = _M0MPC15array13ReadOnlyArray4iterGsE(_M0L5attrsS1501);
  moonbit_decref(_M0L5attrsS1501);
  while (1) {
    moonbit_string_t _M0L4attrS1508;
    moonbit_string_t _M0L7_2abindS1515;
    int32_t _M0L6_2atmpS3445;
    #line 481 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
    _M0L7_2abindS1515 = _M0MPB4Iter4nextGsE(_M0L5_2aitS1506);
    if (_M0L7_2abindS1515 == 0) {
      if (_M0L7_2abindS1515) {
        moonbit_decref(_M0L7_2abindS1515);
      }
      moonbit_decref(_M0L5_2aitS1506);
    } else {
      moonbit_string_t _M0L7_2aSomeS1516 = _M0L7_2abindS1515;
      moonbit_string_t _M0L7_2aattrS1517 = _M0L7_2aSomeS1516;
      _M0L4attrS1508 = _M0L7_2aattrS1517;
      goto join_1507;
    }
    goto joinlet_3769;
    join_1507:;
    _M0L6_2atmpS3445 = Moonbit_array_length(_M0L4attrS1508);
    if (_M0L6_2atmpS3445 >= 5) {
      int32_t _M0L6_2atmpS3451 = _M0L4attrS1508[0];
      int32_t _M0L4_2axS1509 = _M0L6_2atmpS3451;
      if (_M0L4_2axS1509 == 112) {
        int32_t _M0L6_2atmpS3450 = _M0L4attrS1508[1];
        int32_t _M0L4_2axS1510 = _M0L6_2atmpS3450;
        if (_M0L4_2axS1510 == 97) {
          int32_t _M0L6_2atmpS3449 = _M0L4attrS1508[2];
          int32_t _M0L4_2axS1511 = _M0L6_2atmpS3449;
          if (_M0L4_2axS1511 == 110) {
            int32_t _M0L6_2atmpS3448 = _M0L4attrS1508[3];
            int32_t _M0L4_2axS1512 = _M0L6_2atmpS3448;
            if (_M0L4_2axS1512 == 105) {
              int32_t _M0L6_2atmpS3447 = _M0L4attrS1508[4];
              int32_t _M0L4_2axS1513;
              moonbit_decref(_M0L4attrS1508);
              _M0L4_2axS1513 = _M0L6_2atmpS3447;
              if (_M0L4_2axS1513 == 99) {
                void* _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3446;
                struct moonbit_result_0 _result_3770;
                moonbit_decref(_M0L5_2aitS1506);
                moonbit_decref(_M0L1fS1500);
                _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3446
                = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
                Moonbit_object_header(_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3446)->meta
                = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 5);
                ((struct _M0DTPC15error5Error89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3446)->$0
                = _M0L4nameS1502;
                _result_3770.tag = 0;
                _result_3770.data.err
                = _M0L89biolab_2fbio__seq_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS3446;
                return _result_3770;
              }
            } else {
              moonbit_decref(_M0L4attrS1508);
            }
          } else {
            moonbit_decref(_M0L4attrS1508);
          }
        } else {
          moonbit_decref(_M0L4attrS1508);
        }
      } else {
        moonbit_decref(_M0L4attrS1508);
      }
    } else {
      moonbit_decref(_M0L4attrS1508);
    }
    continue;
    joinlet_3769:;
    break;
  }
  #line 487 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L13handle__startS1518->code(_M0L13handle__startS1518);
  moonbit_incref(_M0L14handle__resultS1519);
  moonbit_incref(_M0L4nameS1502);
  _closure_3771
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__*)moonbit_malloc(sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__));
  Moonbit_object_header(_closure_3771)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 9, 0);
  _closure_3771->code
  = &_M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3458l490;
  _closure_3771->$0 = _M0L14handle__resultS1519;
  _closure_3771->$1 = _M0L4nameS1502;
  _M0L6_2atmpS3452 = (struct _M0TWEu*)_closure_3771;
  moonbit_incref(_M0L17error__to__stringS1521);
  moonbit_incref(_M0L14handle__resultS1519);
  _closure_3772
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__*)moonbit_malloc(sizeof(struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__));
  Moonbit_object_header(_closure_3772)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 13, 0);
  _closure_3772->code
  = &_M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3454l491;
  _closure_3772->$0 = _M0L17error__to__stringS1521;
  _closure_3772->$1 = _M0L14handle__resultS1519;
  _closure_3772->$2 = _M0L4nameS1502;
  _M0L6_2atmpS3453 = (struct _M0TWRPC15error5ErrorEu*)_closure_3772;
  #line 488 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FP26biolab8bio__seq45moonbit__test__driver__internal__catch__error(_M0L1fS1500, _M0L6_2atmpS3452, _M0L6_2atmpS3453);
  moonbit_decref(_M0L1fS1500);
  moonbit_decref(_M0L6_2atmpS3452);
  moonbit_decref(_M0L6_2atmpS3453);
  _result_3773.tag = 1;
  _result_3773.data.ok = 1;
  return _result_3773;
}

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3458l490(
  struct _M0TWEu* _M0L6_2aenvS3459
) {
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__* _M0L14_2acasted__envS3460;
  moonbit_string_t _M0L4nameS1502;
  struct _M0TWssbEu* _M0L14handle__resultS1519;
  #line 490 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L14_2acasted__envS3460
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3458__l490__*)_M0L6_2aenvS3459;
  _M0L4nameS1502 = _M0L14_2acasted__envS3460->$1;
  _M0L14handle__resultS1519 = _M0L14_2acasted__envS3460->$0;
  #line 490 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L14handle__resultS1519->code(_M0L14handle__resultS1519, _M0L4nameS1502, (moonbit_string_t)moonbit_string_literal_1.data, 0);
  return 0;
}

int32_t _M0IP26biolab8bio__seq41MoonBit__Test__Driver__Internal__No__ArgsP26biolab8bio__seq21MoonBit__Test__Driver9run__testC3454l491(
  struct _M0TWRPC15error5ErrorEu* _M0L6_2aenvS3455,
  void* _M0L3errS1520
) {
  struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__* _M0L14_2acasted__envS3456;
  moonbit_string_t _M0L4nameS1502;
  struct _M0TWssbEu* _M0L14handle__resultS1519;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1521;
  moonbit_string_t _M0L6_2atmpS3457;
  #line 491 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L14_2acasted__envS3456
  = (struct _M0R151_40biolab_2fbio__seq_2eMoonBit__Test__Driver_3a_3arun__test_7c_40biolab_2fbio__seq_2eMoonBit__Test__Driver__Internal__No__Args_7c_2eanon__u3454__l491__*)_M0L6_2aenvS3455;
  _M0L4nameS1502 = _M0L14_2acasted__envS3456->$2;
  _M0L14handle__resultS1519 = _M0L14_2acasted__envS3456->$1;
  _M0L17error__to__stringS1521 = _M0L14_2acasted__envS3456->$0;
  #line 491 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS3457
  = _M0L17error__to__stringS1521->code(_M0L17error__to__stringS1521, _M0L3errS1520);
  #line 491 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L14handle__resultS1519->code(_M0L14handle__resultS1519, _M0L4nameS1502, _M0L6_2atmpS3457, 0);
  moonbit_decref(_M0L6_2atmpS3457);
  return 0;
}

int32_t _M0FP26biolab8bio__seq45moonbit__test__driver__internal__catch__error(
  struct _M0TWEuQRPC15error5Error* _M0L1fS1471,
  struct _M0TWEu* _M0L6on__okS1472,
  struct _M0TWRPC15error5ErrorEu* _M0L7on__errS1469
) {
  void* _M0L11_2atry__errS1467;
  struct moonbit_result_0 _tmp_3775;
  void* _M0L3errS1468;
  #line 377 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  #line 384 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _tmp_3775 = _M0L1fS1471->code(_M0L1fS1471);
  if (_tmp_3775.tag) {
    int32_t const _M0L5_2aokS3443 = _tmp_3775.data.ok;
  } else {
    void* const _M0L6_2aerrS3444 = _tmp_3775.data.err;
    _M0L11_2atry__errS1467 = _M0L6_2aerrS3444;
    goto join_1466;
  }
  #line 384 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6on__okS1472->code(_M0L6on__okS1472);
  goto joinlet_3774;
  join_1466:;
  _M0L3errS1468 = _M0L11_2atry__errS1467;
  #line 385 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L7on__errS1469->code(_M0L7on__errS1469, _M0L3errS1468);
  moonbit_decref(_M0L3errS1468);
  joinlet_3774:;
  return 0;
}

struct _M0TPB5ArrayGUsiEE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__args(
  
) {
  int32_t _M0L45moonbit__test__driver__internal__parse__int__S1423;
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430;
  int32_t _M0L57moonbit__test__driver__internal__get__cli__args__internalS1437;
  int32_t _M0L51moonbit__test__driver__internal__split__mbt__stringS1444;
  struct _M0TUsiE** _M0L6_2atmpS3442;
  struct _M0TPB5ArrayGUsiEE* _M0L16file__and__indexS1451;
  struct _M0TPB5ArrayGsE* _M0L9cli__argsS1452;
  moonbit_string_t _M0L6_2atmpS3441;
  struct _M0TPB5ArrayGsE* _M0L10test__argsS1453;
  int32_t _M0L7_2abindS1454;
  int32_t _M0L2__S1455;
  #line 194 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L45moonbit__test__driver__internal__parse__int__S1423 = 0;
  _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430
  = 0;
  _M0L57moonbit__test__driver__internal__get__cli__args__internalS1437
  = _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430;
  _M0L51moonbit__test__driver__internal__split__mbt__stringS1444 = 0;
  _M0L6_2atmpS3442 = (struct _M0TUsiE**)moonbit_empty_ref_array;
  _M0L16file__and__indexS1451
  = (struct _M0TPB5ArrayGUsiEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGUsiEE));
  Moonbit_object_header(_M0L16file__and__indexS1451)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L16file__and__indexS1451->$0 = _M0L6_2atmpS3442;
  _M0L16file__and__indexS1451->$1 = 0;
  #line 283 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L9cli__argsS1452
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1437(_M0L57moonbit__test__driver__internal__get__cli__args__internalS1437);
  #line 285 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS3441 = _M0MPC15array5Array2atGsE(_M0L9cli__argsS1452, 1);
  moonbit_decref(_M0L9cli__argsS1452);
  #line 284 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L10test__argsS1453
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1444(_M0L51moonbit__test__driver__internal__split__mbt__stringS1444, _M0L6_2atmpS3441, 47);
  moonbit_decref(_M0L6_2atmpS3441);
  _M0L7_2abindS1454 = _M0L10test__argsS1453->$1;
  _M0L2__S1455 = 0;
  while (1) {
    if (_M0L2__S1455 < _M0L7_2abindS1454) {
      moonbit_string_t* _M0L3bufS3440 = _M0L10test__argsS1453->$0;
      moonbit_string_t _M0L3argS1456 =
        (moonbit_string_t)_M0L3bufS3440[_M0L2__S1455];
      struct _M0TPB5ArrayGsE* _M0L16file__and__rangeS1457;
      moonbit_string_t _M0L4fileS1458;
      moonbit_string_t _M0L5rangeS1459;
      struct _M0TPB5ArrayGsE* _M0L15start__and__endS1460;
      moonbit_string_t _M0L6_2atmpS3438;
      int32_t _M0L5startS1461;
      moonbit_string_t _M0L6_2atmpS3437;
      int32_t _M0L3endS1462;
      int32_t _M0L1iS1463;
      int32_t _M0L6_2atmpS3439;
      moonbit_incref(_M0L3argS1456);
      #line 289 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L16file__and__rangeS1457
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1444(_M0L51moonbit__test__driver__internal__split__mbt__stringS1444, _M0L3argS1456, 58);
      moonbit_decref(_M0L3argS1456);
      #line 290 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L4fileS1458
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1457, 0);
      #line 291 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L5rangeS1459
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1457, 1);
      moonbit_decref(_M0L16file__and__rangeS1457);
      #line 292 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L15start__and__endS1460
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1444(_M0L51moonbit__test__driver__internal__split__mbt__stringS1444, _M0L5rangeS1459, 45);
      moonbit_decref(_M0L5rangeS1459);
      #line 295 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L6_2atmpS3438
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1460, 0);
      #line 295 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L5startS1461
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1423(_M0L45moonbit__test__driver__internal__parse__int__S1423, _M0L6_2atmpS3438);
      moonbit_decref(_M0L6_2atmpS3438);
      #line 296 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L6_2atmpS3437
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1460, 1);
      moonbit_decref(_M0L15start__and__endS1460);
      #line 296 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L3endS1462
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1423(_M0L45moonbit__test__driver__internal__parse__int__S1423, _M0L6_2atmpS3437);
      moonbit_decref(_M0L6_2atmpS3437);
      _M0L1iS1463 = _M0L5startS1461;
      while (1) {
        if (_M0L1iS1463 < _M0L3endS1462) {
          struct _M0TUsiE* _M0L8_2atupleS3435;
          int32_t _M0L6_2atmpS3436;
          moonbit_incref(_M0L4fileS1458);
          _M0L8_2atupleS3435
          = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
          Moonbit_object_header(_M0L8_2atupleS3435)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
          _M0L8_2atupleS3435->$0 = _M0L4fileS1458;
          _M0L8_2atupleS3435->$1 = _M0L1iS1463;
          #line 298 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
          _M0MPC15array5Array4pushGUsiEE(_M0L16file__and__indexS1451, _M0L8_2atupleS3435);
          moonbit_decref(_M0L8_2atupleS3435);
          _M0L6_2atmpS3436 = _M0L1iS1463 + 1;
          _M0L1iS1463 = _M0L6_2atmpS3436;
          continue;
        } else {
          moonbit_decref(_M0L4fileS1458);
        }
        break;
      }
      _M0L6_2atmpS3439 = _M0L2__S1455 + 1;
      _M0L2__S1455 = _M0L6_2atmpS3439;
      continue;
    } else {
      moonbit_decref(_M0L10test__argsS1453);
    }
    break;
  }
  return _M0L16file__and__indexS1451;
}

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1444(
  int32_t _M0L6_2aenvS3416,
  moonbit_string_t _M0L1sS1445,
  int32_t _M0L3sepS1446
) {
  moonbit_string_t* _M0L6_2atmpS3434;
  struct _M0TPB5ArrayGsE* _M0L3resS1447;
  struct _M0TPB8MutLocalGiE* _M0L1iS1448;
  struct _M0TPB8MutLocalGiE* _M0L5startS1449;
  int32_t _M0L3valS3429;
  int32_t _M0L6_2atmpS3430;
  #line 262 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS3434 = (moonbit_string_t*)moonbit_empty_ref_array;
  _M0L3resS1447
  = (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
  Moonbit_object_header(_M0L3resS1447)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
  _M0L3resS1447->$0 = _M0L6_2atmpS3434;
  _M0L3resS1447->$1 = 0;
  _M0L1iS1448
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1448)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1448->$0 = 0;
  _M0L5startS1449
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS1449)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS1449->$0 = 0;
  while (1) {
    int32_t _M0L3valS3417 = _M0L1iS1448->$0;
    int32_t _M0L6_2atmpS3418 = Moonbit_array_length(_M0L1sS1445);
    if (_M0L3valS3417 < _M0L6_2atmpS3418) {
      int32_t _M0L3valS3421 = _M0L1iS1448->$0;
      int32_t _M0L6_2atmpS3420;
      int32_t _M0L6_2atmpS3419;
      int32_t _M0L3valS3428;
      int32_t _M0L6_2atmpS3427;
      if (
        _M0L3valS3421 < 0
        || _M0L3valS3421 >= Moonbit_array_length(_M0L1sS1445)
      ) {
        #line 270 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3420 = _M0L1sS1445[_M0L3valS3421];
      _M0L6_2atmpS3419 = _M0L6_2atmpS3420;
      if (_M0L6_2atmpS3419 == _M0L3sepS1446) {
        int32_t _M0L3valS3423 = _M0L5startS1449->$0;
        int32_t _M0L3valS3424 = _M0L1iS1448->$0;
        moonbit_string_t _M0L6_2atmpS3422;
        int32_t _M0L3valS3426;
        int32_t _M0L6_2atmpS3425;
        #line 271 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        _M0L6_2atmpS3422
        = _M0MPC16string6String17unsafe__substring(_M0L1sS1445, _M0L3valS3423, _M0L3valS3424);
        #line 271 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        _M0MPC15array5Array4pushGsE(_M0L3resS1447, _M0L6_2atmpS3422);
        moonbit_decref(_M0L6_2atmpS3422);
        _M0L3valS3426 = _M0L1iS1448->$0;
        _M0L6_2atmpS3425 = _M0L3valS3426 + 1;
        _M0L5startS1449->$0 = _M0L6_2atmpS3425;
      }
      _M0L3valS3428 = _M0L1iS1448->$0;
      _M0L6_2atmpS3427 = _M0L3valS3428 + 1;
      _M0L1iS1448->$0 = _M0L6_2atmpS3427;
      continue;
    } else {
      moonbit_decref(_M0L1iS1448);
    }
    break;
  }
  _M0L3valS3429 = _M0L5startS1449->$0;
  _M0L6_2atmpS3430 = Moonbit_array_length(_M0L1sS1445);
  if (_M0L3valS3429 < _M0L6_2atmpS3430) {
    int32_t _M0L3valS3432 = _M0L5startS1449->$0;
    int32_t _M0L6_2atmpS3433;
    moonbit_string_t _M0L6_2atmpS3431;
    moonbit_decref(_M0L5startS1449);
    _M0L6_2atmpS3433 = Moonbit_array_length(_M0L1sS1445);
    #line 277 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
    _M0L6_2atmpS3431
    = _M0MPC16string6String17unsafe__substring(_M0L1sS1445, _M0L3valS3432, _M0L6_2atmpS3433);
    #line 277 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
    _M0MPC15array5Array4pushGsE(_M0L3resS1447, _M0L6_2atmpS3431);
    moonbit_decref(_M0L6_2atmpS3431);
  } else {
    moonbit_decref(_M0L5startS1449);
  }
  return _M0L3resS1447;
}

struct _M0TPB5ArrayGsE* _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1437(
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430
) {
  moonbit_bytes_t* _M0L3tmpS1438;
  int32_t _M0L6_2atmpS3415;
  struct _M0TPB5ArrayGsE* _M0L3resS1439;
  int32_t _M0L7_2abindS1440;
  int32_t _M0L7_2abindS1441;
  int32_t _M0L1iS1442;
  #line 251 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  #line 254 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L3tmpS1438
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__get__cli__args__ffi();
  _M0L6_2atmpS3415 = Moonbit_array_length(_M0L3tmpS1438);
  #line 255 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L3resS1439 = _M0MPC15array5Array11new_2einnerGsE(_M0L6_2atmpS3415);
  _M0L7_2abindS1440 = 0;
  _M0L7_2abindS1441 = Moonbit_array_length(_M0L3tmpS1438);
  _M0L1iS1442 = _M0L7_2abindS1440;
  while (1) {
    if (_M0L1iS1442 < _M0L7_2abindS1441) {
      moonbit_bytes_t _M0L6_2atmpS3413;
      moonbit_string_t _M0L6_2atmpS3412;
      int32_t _M0L6_2atmpS3414;
      if (
        _M0L1iS1442 < 0 || _M0L1iS1442 >= Moonbit_array_length(_M0L3tmpS1438)
      ) {
        #line 257 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3413 = (moonbit_bytes_t)_M0L3tmpS1438[_M0L1iS1442];
      moonbit_incref(_M0L6_2atmpS3413);
      #line 257 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0L6_2atmpS3412
      = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430(_M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430, _M0L6_2atmpS3413);
      moonbit_decref(_M0L6_2atmpS3413);
      #line 257 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0MPC15array5Array4pushGsE(_M0L3resS1439, _M0L6_2atmpS3412);
      moonbit_decref(_M0L6_2atmpS3412);
      _M0L6_2atmpS3414 = _M0L1iS1442 + 1;
      _M0L1iS1442 = _M0L6_2atmpS3414;
      continue;
    } else {
      moonbit_decref(_M0L3tmpS1438);
    }
    break;
  }
  return _M0L3resS1439;
}

moonbit_string_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1430(
  int32_t _M0L6_2aenvS3326,
  moonbit_bytes_t _M0L5bytesS1431
) {
  struct _M0TPB13StringBuilder* _M0L3resS1432;
  int32_t _M0L3lenS1433;
  struct _M0TPB8MutLocalGiE* _M0L1iS1434;
  moonbit_string_t _result_3781;
  #line 207 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  #line 210 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L3resS1432 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3lenS1433 = Moonbit_array_length(_M0L5bytesS1431);
  _M0L1iS1434
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1434)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1434->$0 = 0;
  while (1) {
    int32_t _M0L3valS3327 = _M0L1iS1434->$0;
    if (_M0L3valS3327 < _M0L3lenS1433) {
      int32_t _M0L3valS3411 = _M0L1iS1434->$0;
      int32_t _M0L6_2atmpS3410;
      int32_t _M0L6_2atmpS3409;
      struct _M0TPB8MutLocalGiE* _M0L1cS1435;
      int32_t _M0L3valS3328;
      if (
        _M0L3valS3411 < 0
        || _M0L3valS3411 >= Moonbit_array_length(_M0L5bytesS1431)
      ) {
        #line 214 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3410 = _M0L5bytesS1431[_M0L3valS3411];
      _M0L6_2atmpS3409 = (int32_t)_M0L6_2atmpS3410;
      _M0L1cS1435
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1cS1435)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1cS1435->$0 = _M0L6_2atmpS3409;
      _M0L3valS3328 = _M0L1cS1435->$0;
      if (_M0L3valS3328 < 128) {
        int32_t _M0L3valS3330 = _M0L1cS1435->$0;
        int32_t _M0L6_2atmpS3329;
        int32_t _M0L3valS3332;
        int32_t _M0L6_2atmpS3331;
        moonbit_decref(_M0L1cS1435);
        _M0L6_2atmpS3329 = _M0L3valS3330;
        #line 216 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1432, _M0L6_2atmpS3329);
        _M0L3valS3332 = _M0L1iS1434->$0;
        _M0L6_2atmpS3331 = _M0L3valS3332 + 1;
        _M0L1iS1434->$0 = _M0L6_2atmpS3331;
      } else {
        int32_t _M0L3valS3333 = _M0L1cS1435->$0;
        if (_M0L3valS3333 < 224) {
          int32_t _M0L3valS3335 = _M0L1iS1434->$0;
          int32_t _M0L6_2atmpS3334 = _M0L3valS3335 + 1;
          int32_t _M0L3valS3344;
          int32_t _M0L6_2atmpS3343;
          int32_t _M0L6_2atmpS3337;
          int32_t _M0L3valS3342;
          int32_t _M0L6_2atmpS3341;
          int32_t _M0L6_2atmpS3340;
          int32_t _M0L6_2atmpS3339;
          int32_t _M0L6_2atmpS3338;
          int32_t _M0L6_2atmpS3336;
          int32_t _M0L3valS3346;
          int32_t _M0L6_2atmpS3345;
          int32_t _M0L3valS3348;
          int32_t _M0L6_2atmpS3347;
          if (_M0L6_2atmpS3334 >= _M0L3lenS1433) {
            moonbit_decref(_M0L1cS1435);
            moonbit_decref(_M0L1iS1434);
            break;
          }
          _M0L3valS3344 = _M0L1cS1435->$0;
          _M0L6_2atmpS3343 = _M0L3valS3344 & 31;
          _M0L6_2atmpS3337 = _M0L6_2atmpS3343 << 6;
          _M0L3valS3342 = _M0L1iS1434->$0;
          _M0L6_2atmpS3341 = _M0L3valS3342 + 1;
          if (
            _M0L6_2atmpS3341 < 0
            || _M0L6_2atmpS3341 >= Moonbit_array_length(_M0L5bytesS1431)
          ) {
            #line 222 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS3340 = _M0L5bytesS1431[_M0L6_2atmpS3341];
          _M0L6_2atmpS3339 = (int32_t)_M0L6_2atmpS3340;
          _M0L6_2atmpS3338 = _M0L6_2atmpS3339 & 63;
          _M0L6_2atmpS3336 = _M0L6_2atmpS3337 | _M0L6_2atmpS3338;
          _M0L1cS1435->$0 = _M0L6_2atmpS3336;
          _M0L3valS3346 = _M0L1cS1435->$0;
          moonbit_decref(_M0L1cS1435);
          _M0L6_2atmpS3345 = _M0L3valS3346;
          #line 223 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1432, _M0L6_2atmpS3345);
          _M0L3valS3348 = _M0L1iS1434->$0;
          _M0L6_2atmpS3347 = _M0L3valS3348 + 2;
          _M0L1iS1434->$0 = _M0L6_2atmpS3347;
        } else {
          int32_t _M0L3valS3349 = _M0L1cS1435->$0;
          if (_M0L3valS3349 < 240) {
            int32_t _M0L3valS3351 = _M0L1iS1434->$0;
            int32_t _M0L6_2atmpS3350 = _M0L3valS3351 + 2;
            int32_t _M0L3valS3367;
            int32_t _M0L6_2atmpS3366;
            int32_t _M0L6_2atmpS3359;
            int32_t _M0L3valS3365;
            int32_t _M0L6_2atmpS3364;
            int32_t _M0L6_2atmpS3363;
            int32_t _M0L6_2atmpS3362;
            int32_t _M0L6_2atmpS3361;
            int32_t _M0L6_2atmpS3360;
            int32_t _M0L6_2atmpS3353;
            int32_t _M0L3valS3358;
            int32_t _M0L6_2atmpS3357;
            int32_t _M0L6_2atmpS3356;
            int32_t _M0L6_2atmpS3355;
            int32_t _M0L6_2atmpS3354;
            int32_t _M0L6_2atmpS3352;
            int32_t _M0L3valS3369;
            int32_t _M0L6_2atmpS3368;
            int32_t _M0L3valS3371;
            int32_t _M0L6_2atmpS3370;
            if (_M0L6_2atmpS3350 >= _M0L3lenS1433) {
              moonbit_decref(_M0L1cS1435);
              moonbit_decref(_M0L1iS1434);
              break;
            }
            _M0L3valS3367 = _M0L1cS1435->$0;
            _M0L6_2atmpS3366 = _M0L3valS3367 & 15;
            _M0L6_2atmpS3359 = _M0L6_2atmpS3366 << 12;
            _M0L3valS3365 = _M0L1iS1434->$0;
            _M0L6_2atmpS3364 = _M0L3valS3365 + 1;
            if (
              _M0L6_2atmpS3364 < 0
              || _M0L6_2atmpS3364 >= Moonbit_array_length(_M0L5bytesS1431)
            ) {
              #line 230 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3363 = _M0L5bytesS1431[_M0L6_2atmpS3364];
            _M0L6_2atmpS3362 = (int32_t)_M0L6_2atmpS3363;
            _M0L6_2atmpS3361 = _M0L6_2atmpS3362 & 63;
            _M0L6_2atmpS3360 = _M0L6_2atmpS3361 << 6;
            _M0L6_2atmpS3353 = _M0L6_2atmpS3359 | _M0L6_2atmpS3360;
            _M0L3valS3358 = _M0L1iS1434->$0;
            _M0L6_2atmpS3357 = _M0L3valS3358 + 2;
            if (
              _M0L6_2atmpS3357 < 0
              || _M0L6_2atmpS3357 >= Moonbit_array_length(_M0L5bytesS1431)
            ) {
              #line 231 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3356 = _M0L5bytesS1431[_M0L6_2atmpS3357];
            _M0L6_2atmpS3355 = (int32_t)_M0L6_2atmpS3356;
            _M0L6_2atmpS3354 = _M0L6_2atmpS3355 & 63;
            _M0L6_2atmpS3352 = _M0L6_2atmpS3353 | _M0L6_2atmpS3354;
            _M0L1cS1435->$0 = _M0L6_2atmpS3352;
            _M0L3valS3369 = _M0L1cS1435->$0;
            moonbit_decref(_M0L1cS1435);
            _M0L6_2atmpS3368 = _M0L3valS3369;
            #line 232 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1432, _M0L6_2atmpS3368);
            _M0L3valS3371 = _M0L1iS1434->$0;
            _M0L6_2atmpS3370 = _M0L3valS3371 + 3;
            _M0L1iS1434->$0 = _M0L6_2atmpS3370;
          } else {
            int32_t _M0L3valS3373 = _M0L1iS1434->$0;
            int32_t _M0L6_2atmpS3372 = _M0L3valS3373 + 3;
            int32_t _M0L3valS3396;
            int32_t _M0L6_2atmpS3395;
            int32_t _M0L6_2atmpS3388;
            int32_t _M0L3valS3394;
            int32_t _M0L6_2atmpS3393;
            int32_t _M0L6_2atmpS3392;
            int32_t _M0L6_2atmpS3391;
            int32_t _M0L6_2atmpS3390;
            int32_t _M0L6_2atmpS3389;
            int32_t _M0L6_2atmpS3381;
            int32_t _M0L3valS3387;
            int32_t _M0L6_2atmpS3386;
            int32_t _M0L6_2atmpS3385;
            int32_t _M0L6_2atmpS3384;
            int32_t _M0L6_2atmpS3383;
            int32_t _M0L6_2atmpS3382;
            int32_t _M0L6_2atmpS3375;
            int32_t _M0L3valS3380;
            int32_t _M0L6_2atmpS3379;
            int32_t _M0L6_2atmpS3378;
            int32_t _M0L6_2atmpS3377;
            int32_t _M0L6_2atmpS3376;
            int32_t _M0L6_2atmpS3374;
            int32_t _M0L3valS3398;
            int32_t _M0L6_2atmpS3397;
            int32_t _M0L3valS3402;
            int32_t _M0L6_2atmpS3401;
            int32_t _M0L6_2atmpS3400;
            int32_t _M0L6_2atmpS3399;
            int32_t _M0L3valS3406;
            int32_t _M0L6_2atmpS3405;
            int32_t _M0L6_2atmpS3404;
            int32_t _M0L6_2atmpS3403;
            int32_t _M0L3valS3408;
            int32_t _M0L6_2atmpS3407;
            if (_M0L6_2atmpS3372 >= _M0L3lenS1433) {
              moonbit_decref(_M0L1cS1435);
              moonbit_decref(_M0L1iS1434);
              break;
            }
            _M0L3valS3396 = _M0L1cS1435->$0;
            _M0L6_2atmpS3395 = _M0L3valS3396 & 7;
            _M0L6_2atmpS3388 = _M0L6_2atmpS3395 << 18;
            _M0L3valS3394 = _M0L1iS1434->$0;
            _M0L6_2atmpS3393 = _M0L3valS3394 + 1;
            if (
              _M0L6_2atmpS3393 < 0
              || _M0L6_2atmpS3393 >= Moonbit_array_length(_M0L5bytesS1431)
            ) {
              #line 239 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3392 = _M0L5bytesS1431[_M0L6_2atmpS3393];
            _M0L6_2atmpS3391 = (int32_t)_M0L6_2atmpS3392;
            _M0L6_2atmpS3390 = _M0L6_2atmpS3391 & 63;
            _M0L6_2atmpS3389 = _M0L6_2atmpS3390 << 12;
            _M0L6_2atmpS3381 = _M0L6_2atmpS3388 | _M0L6_2atmpS3389;
            _M0L3valS3387 = _M0L1iS1434->$0;
            _M0L6_2atmpS3386 = _M0L3valS3387 + 2;
            if (
              _M0L6_2atmpS3386 < 0
              || _M0L6_2atmpS3386 >= Moonbit_array_length(_M0L5bytesS1431)
            ) {
              #line 240 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3385 = _M0L5bytesS1431[_M0L6_2atmpS3386];
            _M0L6_2atmpS3384 = (int32_t)_M0L6_2atmpS3385;
            _M0L6_2atmpS3383 = _M0L6_2atmpS3384 & 63;
            _M0L6_2atmpS3382 = _M0L6_2atmpS3383 << 6;
            _M0L6_2atmpS3375 = _M0L6_2atmpS3381 | _M0L6_2atmpS3382;
            _M0L3valS3380 = _M0L1iS1434->$0;
            _M0L6_2atmpS3379 = _M0L3valS3380 + 3;
            if (
              _M0L6_2atmpS3379 < 0
              || _M0L6_2atmpS3379 >= Moonbit_array_length(_M0L5bytesS1431)
            ) {
              #line 241 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS3378 = _M0L5bytesS1431[_M0L6_2atmpS3379];
            _M0L6_2atmpS3377 = (int32_t)_M0L6_2atmpS3378;
            _M0L6_2atmpS3376 = _M0L6_2atmpS3377 & 63;
            _M0L6_2atmpS3374 = _M0L6_2atmpS3375 | _M0L6_2atmpS3376;
            _M0L1cS1435->$0 = _M0L6_2atmpS3374;
            _M0L3valS3398 = _M0L1cS1435->$0;
            _M0L6_2atmpS3397 = _M0L3valS3398 - 65536;
            _M0L1cS1435->$0 = _M0L6_2atmpS3397;
            _M0L3valS3402 = _M0L1cS1435->$0;
            _M0L6_2atmpS3401 = _M0L3valS3402 >> 10;
            _M0L6_2atmpS3400 = _M0L6_2atmpS3401 + 55296;
            _M0L6_2atmpS3399 = _M0L6_2atmpS3400;
            #line 243 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1432, _M0L6_2atmpS3399);
            _M0L3valS3406 = _M0L1cS1435->$0;
            moonbit_decref(_M0L1cS1435);
            _M0L6_2atmpS3405 = _M0L3valS3406 & 1023;
            _M0L6_2atmpS3404 = _M0L6_2atmpS3405 + 56320;
            _M0L6_2atmpS3403 = _M0L6_2atmpS3404;
            #line 244 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1432, _M0L6_2atmpS3403);
            _M0L3valS3408 = _M0L1iS1434->$0;
            _M0L6_2atmpS3407 = _M0L3valS3408 + 4;
            _M0L1iS1434->$0 = _M0L6_2atmpS3407;
          }
        }
      }
      continue;
    } else {
      moonbit_decref(_M0L1iS1434);
    }
    break;
  }
  #line 248 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _result_3781 = _M0MPB13StringBuilder10to__string(_M0L3resS1432);
  moonbit_decref(_M0L3resS1432);
  return _result_3781;
}

int32_t _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1423(
  int32_t _M0L6_2aenvS3319,
  moonbit_string_t _M0L1sS1424
) {
  struct _M0TPB8MutLocalGiE* _M0L3resS1425;
  int32_t _M0L3lenS1426;
  int32_t _M0L7_2abindS1427;
  int32_t _M0L1iS1428;
  int32_t _result_3783;
  #line 198 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L3resS1425
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3resS1425)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3resS1425->$0 = 0;
  _M0L3lenS1426 = Moonbit_array_length(_M0L1sS1424);
  _M0L7_2abindS1427 = 0;
  _M0L1iS1428 = _M0L7_2abindS1427;
  while (1) {
    if (_M0L1iS1428 < _M0L3lenS1426) {
      int32_t _M0L3valS3324 = _M0L3resS1425->$0;
      int32_t _M0L6_2atmpS3321 = _M0L3valS3324 * 10;
      int32_t _M0L6_2atmpS3323;
      int32_t _M0L6_2atmpS3322;
      int32_t _M0L6_2atmpS3320;
      int32_t _M0L6_2atmpS3325;
      if (
        _M0L1iS1428 < 0 || _M0L1iS1428 >= Moonbit_array_length(_M0L1sS1424)
      ) {
        #line 202 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS3323 = _M0L1sS1424[_M0L1iS1428];
      _M0L6_2atmpS3322 = _M0L6_2atmpS3323 - 48;
      _M0L6_2atmpS3320 = _M0L6_2atmpS3321 + _M0L6_2atmpS3322;
      _M0L3resS1425->$0 = _M0L6_2atmpS3320;
      _M0L6_2atmpS3325 = _M0L1iS1428 + 1;
      _M0L1iS1428 = _M0L6_2atmpS3325;
      continue;
    }
    break;
  }
  _result_3783 = _M0L3resS1425->$0;
  moonbit_decref(_M0L3resS1425);
  return _result_3783;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1399,
  moonbit_string_t _M0L12_2adiscard__S1400,
  int32_t _M0L12_2adiscard__S1401,
  struct _M0TWEu* _M0L12_2adiscard__S1402,
  struct _M0TWssbEu* _M0L12_2adiscard__S1403,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1404
) {
  struct moonbit_result_0 _result_3784;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _result_3784.tag = 1;
  _result_3784.data.ok = 0;
  return _result_3784;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1405,
  moonbit_string_t _M0L12_2adiscard__S1406,
  int32_t _M0L12_2adiscard__S1407,
  struct _M0TWEu* _M0L12_2adiscard__S1408,
  struct _M0TWssbEu* _M0L12_2adiscard__S1409,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1410
) {
  struct moonbit_result_0 _result_3785;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _result_3785.tag = 1;
  _result_3785.data.ok = 0;
  return _result_3785;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1411,
  moonbit_string_t _M0L12_2adiscard__S1412,
  int32_t _M0L12_2adiscard__S1413,
  struct _M0TWEu* _M0L12_2adiscard__S1414,
  struct _M0TWssbEu* _M0L12_2adiscard__S1415,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1416
) {
  struct moonbit_result_0 _result_3786;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _result_3786.tag = 1;
  _result_3786.data.ok = 0;
  return _result_3786;
}

struct moonbit_result_0 _M0IP016_24default__implP26biolab8bio__seq21MoonBit__Test__Driver9run__testGRP26biolab8bio__seq50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1417,
  moonbit_string_t _M0L12_2adiscard__S1418,
  int32_t _M0L12_2adiscard__S1419,
  struct _M0TWEu* _M0L12_2adiscard__S1420,
  struct _M0TWssbEu* _M0L12_2adiscard__S1421,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1422
) {
  struct moonbit_result_0 _result_3787;
  #line 35 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _result_3787.tag = 1;
  _result_3787.data.ok = 0;
  return _result_3787;
}

int32_t _M0IP016_24default__implP26biolab8bio__seq28MoonBit__Async__Test__Driver17run__async__testsGRP26biolab8bio__seq34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1398
) {
  #line 12 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  return 0;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__17(
  
) {
  moonbit_string_t _M0L6newickS1395;
  struct moonbit_result_1 _tmp_3788;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1396;
  void* _M0L4distS1397;
  int32_t _M0L6_2atmpS3315;
  void* _M0L4NoneS3316;
  struct moonbit_result_0 _result_3790;
  #line 214 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 215 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1395 = _M0FP26biolab8bio__seq24gen__caterpillar__newick(500);
  #line 216 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3788 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1395);
  moonbit_decref(_M0L6newickS1395);
  if (_tmp_3788.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3317 =
      _tmp_3788.data.ok;
    _M0L4treeS1396 = _M0L5_2aokS3317;
  } else {
    void* const _M0L6_2aerrS3318 = _tmp_3788.data.err;
    struct moonbit_result_0 _result_3789;
    _result_3789.tag = 0;
    _result_3789.data.err = _M0L6_2aerrS3318;
    return _result_3789;
  }
  #line 217 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L4distS1397
  = _M0MP26biolab8bio__seq8TreeNode8distance(_M0L4treeS1396, (moonbit_string_t)moonbit_string_literal_10.data, (moonbit_string_t)moonbit_string_literal_11.data);
  moonbit_decref(_M0L4treeS1396);
  switch (Moonbit_object_tag(_M0L4distS1397)) {
    case 1: {
      moonbit_decref(_M0L4distS1397);
      _M0L6_2atmpS3315 = 1;
      break;
    }
    default: {
      moonbit_decref(_M0L4distS1397);
      _M0L6_2atmpS3315 = 0;
      break;
    }
  }
  _M0L4NoneS3316
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 218 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3790
  = _M0FPB12assert__true(_M0L6_2atmpS3315, _M0L4NoneS3316, (moonbit_string_t)moonbit_string_literal_12.data);
  moonbit_decref(_M0L4NoneS3316);
  return _result_3790;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__16(
  
) {
  moonbit_string_t _M0L6newickS1391;
  struct moonbit_result_1 _tmp_3791;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1392;
  struct _M0TPB5ArrayGsE* _M0L5namesS1393;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L8ancestorS1394;
  int32_t _M0L6_2atmpS3312;
  int32_t _M0L6_2atmpS3310;
  void* _M0L4NoneS3311;
  struct moonbit_result_0 _result_3793;
  #line 203 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 204 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1391 = _M0FP26biolab8bio__seq24gen__caterpillar__newick(500);
  #line 205 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3791 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1391);
  moonbit_decref(_M0L6newickS1391);
  if (_tmp_3791.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3313 =
      _tmp_3791.data.ok;
    _M0L4treeS1392 = _M0L5_2aokS3313;
  } else {
    void* const _M0L6_2aerrS3314 = _tmp_3791.data.err;
    struct moonbit_result_0 _result_3792;
    _result_3792.tag = 0;
    _result_3792.data.err = _M0L6_2aerrS3314;
    return _result_3792;
  }
  #line 206 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5namesS1393 = _M0MPC15array5Array11new_2einnerGsE(0);
  #line 207 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0MPC15array5Array4pushGsE(_M0L5namesS1393, (moonbit_string_t)moonbit_string_literal_10.data);
  #line 208 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0MPC15array5Array4pushGsE(_M0L5namesS1393, (moonbit_string_t)moonbit_string_literal_11.data);
  #line 209 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L8ancestorS1394
  = _M0MP26biolab8bio__seq8TreeNode3lca(_M0L4treeS1392, _M0L5namesS1393);
  moonbit_decref(_M0L4treeS1392);
  moonbit_decref(_M0L5namesS1393);
  _M0L6_2atmpS3312 = _M0L8ancestorS1394 == 0;
  if (_M0L8ancestorS1394) {
    moonbit_decref(_M0L8ancestorS1394);
  }
  _M0L6_2atmpS3310 = !_M0L6_2atmpS3312;
  _M0L4NoneS3311
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 210 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3793
  = _M0FPB12assert__true(_M0L6_2atmpS3310, _M0L4NoneS3311, (moonbit_string_t)moonbit_string_literal_13.data);
  moonbit_decref(_M0L4NoneS3311);
  return _result_3793;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__15(
  
) {
  moonbit_string_t _M0L6newickS1388;
  struct moonbit_result_1 _tmp_3794;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1389;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L5foundS1390;
  int32_t _M0L6_2atmpS3307;
  int32_t _M0L6_2atmpS3305;
  void* _M0L4NoneS3306;
  struct moonbit_result_0 _result_3796;
  #line 195 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 196 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1388 = _M0FP26biolab8bio__seq24gen__caterpillar__newick(500);
  #line 197 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3794 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1388);
  moonbit_decref(_M0L6newickS1388);
  if (_tmp_3794.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3308 =
      _tmp_3794.data.ok;
    _M0L4treeS1389 = _M0L5_2aokS3308;
  } else {
    void* const _M0L6_2aerrS3309 = _tmp_3794.data.err;
    struct moonbit_result_0 _result_3795;
    _result_3795.tag = 0;
    _result_3795.data.err = _M0L6_2aerrS3309;
    return _result_3795;
  }
  #line 198 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5foundS1390
  = _M0MP26biolab8bio__seq8TreeNode4find(_M0L4treeS1389, (moonbit_string_t)moonbit_string_literal_11.data);
  moonbit_decref(_M0L4treeS1389);
  _M0L6_2atmpS3307 = _M0L5foundS1390 == 0;
  if (_M0L5foundS1390) {
    moonbit_decref(_M0L5foundS1390);
  }
  _M0L6_2atmpS3305 = !_M0L6_2atmpS3307;
  _M0L4NoneS3306
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 199 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3796
  = _M0FPB12assert__true(_M0L6_2atmpS3305, _M0L4NoneS3306, (moonbit_string_t)moonbit_string_literal_14.data);
  moonbit_decref(_M0L4NoneS3306);
  return _result_3796;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__14(
  
) {
  moonbit_string_t _M0L6newickS1386;
  struct moonbit_result_1 _tmp_3797;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1387;
  int32_t _M0L6_2atmpS3301;
  moonbit_string_t _M0L6_2atmpS3302;
  struct moonbit_result_0 _result_3799;
  #line 188 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 189 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1386 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 190 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3797 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1386);
  moonbit_decref(_M0L6newickS1386);
  if (_tmp_3797.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3303 =
      _tmp_3797.data.ok;
    _M0L4treeS1387 = _M0L5_2aokS3303;
  } else {
    void* const _M0L6_2aerrS3304 = _tmp_3797.data.err;
    struct moonbit_result_0 _result_3798;
    _result_3798.tag = 0;
    _result_3798.data.err = _M0L6_2aerrS3304;
    return _result_3798;
  }
  #line 191 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3301
  = _M0MP26biolab8bio__seq8TreeNode11count__tips(_M0L4treeS1387);
  moonbit_decref(_M0L4treeS1387);
  _M0L6_2atmpS3302 = 0;
  #line 191 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3799
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3301, 4096, _M0L6_2atmpS3302, (moonbit_string_t)moonbit_string_literal_15.data);
  if (_M0L6_2atmpS3302) {
    moonbit_decref(_M0L6_2atmpS3302);
  }
  return _result_3799;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__13(
  
) {
  moonbit_string_t _M0L6newickS1383;
  struct moonbit_result_1 _tmp_3800;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1384;
  double _M0L5totalS1385;
  int32_t _M0L6_2atmpS3297;
  void* _M0L4NoneS3298;
  struct moonbit_result_0 _result_3802;
  #line 180 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 181 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1383 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 182 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3800 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1383);
  moonbit_decref(_M0L6newickS1383);
  if (_tmp_3800.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3299 =
      _tmp_3800.data.ok;
    _M0L4treeS1384 = _M0L5_2aokS3299;
  } else {
    void* const _M0L6_2aerrS3300 = _tmp_3800.data.err;
    struct moonbit_result_0 _result_3801;
    _result_3801.tag = 0;
    _result_3801.data.err = _M0L6_2aerrS3300;
    return _result_3801;
  }
  #line 183 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5totalS1385
  = _M0MP26biolab8bio__seq8TreeNode13total__length(_M0L4treeS1384);
  moonbit_decref(_M0L4treeS1384);
  _M0L6_2atmpS3297 = _M0L5totalS1385 > 0x0p+0;
  _M0L4NoneS3298
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 184 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3802
  = _M0FPB12assert__true(_M0L6_2atmpS3297, _M0L4NoneS3298, (moonbit_string_t)moonbit_string_literal_16.data);
  moonbit_decref(_M0L4NoneS3298);
  return _result_3802;
}

double _M0MP26biolab8bio__seq8TreeNode13total__length(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1377
) {
  void* _M0L6lengthS3296;
  double _M0L6_2atmpS3295;
  struct _M0TPB8MutLocalGdE* _M0L5totalS1376;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1378;
  int32_t _M0L7_2abindS1379;
  int32_t _M0L2__S1380;
  double _result_3804;
  #line 261 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6lengthS3296 = _M0L4selfS1377->$1;
  moonbit_incref(_M0L6lengthS3296);
  #line 262 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3295
  = _M0MPC16option6Option10unwrap__orGdE(_M0L6lengthS3296, 0x0p+0);
  moonbit_decref(_M0L6lengthS3296);
  _M0L5totalS1376
  = (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
  Moonbit_object_header(_M0L5totalS1376)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5totalS1376->$0 = _M0L6_2atmpS3295;
  _M0L7_2abindS1378 = _M0L4selfS1377->$3;
  _M0L7_2abindS1379 = _M0L7_2abindS1378->$1;
  moonbit_incref(_M0L7_2abindS1378);
  _M0L2__S1380 = 0;
  while (1) {
    if (_M0L2__S1380 < _M0L7_2abindS1379) {
      struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3294 =
        _M0L7_2abindS1378->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1381 =
        (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3294[_M0L2__S1380];
      double _M0L3valS3291 = _M0L5totalS1376->$0;
      double _M0L6_2atmpS3292;
      double _M0L6_2atmpS3290;
      int32_t _M0L6_2atmpS3293;
      moonbit_incref(_M0L5childS1381);
      #line 264 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3292
      = _M0MP26biolab8bio__seq8TreeNode13total__length(_M0L5childS1381);
      moonbit_decref(_M0L5childS1381);
      _M0L6_2atmpS3290 = _M0L3valS3291 + _M0L6_2atmpS3292;
      _M0L5totalS1376->$0 = _M0L6_2atmpS3290;
      _M0L6_2atmpS3293 = _M0L2__S1380 + 1;
      _M0L2__S1380 = _M0L6_2atmpS3293;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1378);
    }
    break;
  }
  _result_3804 = _M0L5totalS1376->$0;
  moonbit_decref(_M0L5totalS1376);
  return _result_3804;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__12(
  
) {
  moonbit_string_t _M0L6newickS1373;
  struct moonbit_result_1 _tmp_3805;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1374;
  void* _M0L4distS1375;
  int32_t _M0L6_2atmpS3286;
  void* _M0L4NoneS3287;
  struct moonbit_result_0 _result_3807;
  #line 172 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 173 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1373 = _M0FP26biolab8bio__seq21gen__balanced__newick(10);
  #line 174 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3805 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1373);
  moonbit_decref(_M0L6newickS1373);
  if (_tmp_3805.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3288 =
      _tmp_3805.data.ok;
    _M0L4treeS1374 = _M0L5_2aokS3288;
  } else {
    void* const _M0L6_2aerrS3289 = _tmp_3805.data.err;
    struct moonbit_result_0 _result_3806;
    _result_3806.tag = 0;
    _result_3806.data.err = _M0L6_2aerrS3289;
    return _result_3806;
  }
  #line 175 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L4distS1375
  = _M0MP26biolab8bio__seq8TreeNode8distance(_M0L4treeS1374, (moonbit_string_t)moonbit_string_literal_17.data, (moonbit_string_t)moonbit_string_literal_18.data);
  moonbit_decref(_M0L4treeS1374);
  switch (Moonbit_object_tag(_M0L4distS1375)) {
    case 1: {
      moonbit_decref(_M0L4distS1375);
      _M0L6_2atmpS3286 = 1;
      break;
    }
    default: {
      moonbit_decref(_M0L4distS1375);
      _M0L6_2atmpS3286 = 0;
      break;
    }
  }
  _M0L4NoneS3287
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 176 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3807
  = _M0FPB12assert__true(_M0L6_2atmpS3286, _M0L4NoneS3287, (moonbit_string_t)moonbit_string_literal_19.data);
  moonbit_decref(_M0L4NoneS3287);
  return _result_3807;
}

void* _M0MP26biolab8bio__seq8TreeNode8distance(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1365,
  moonbit_string_t _M0L5name1S1366,
  moonbit_string_t _M0L5name2S1368
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5path1S1348;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5path2S1349;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1364;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1367;
  struct _M0TPB8MutLocalGiE* _M0L8lca__idxS1350;
  int32_t _M0L6_2atmpS3284;
  int32_t _M0L6_2atmpS3285;
  int32_t _M0L8min__lenS1351;
  int32_t _M0L7_2abindS1352;
  int32_t _M0L1iS1353;
  struct _M0TPB8MutLocalGdE* _M0L4distS1355;
  int32_t _M0L3valS3275;
  int32_t _M0L7_2abindS1356;
  int32_t _M0L7_2abindS1357;
  int32_t _M0L1iS1358;
  int32_t _M0L3valS3282;
  int32_t _M0L7_2abindS1360;
  int32_t _M0L7_2abindS1361;
  int32_t _M0L1iS1362;
  double _M0L3valS3283;
  void* _block_3814;
  #line 224 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 229 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7_2abindS1364
  = _M0MP26biolab8bio__seq8TreeNode10find__path(_M0L4selfS1365, _M0L5name1S1366);
  #line 229 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7_2abindS1367
  = _M0MP26biolab8bio__seq8TreeNode10find__path(_M0L4selfS1365, _M0L5name2S1368);
  if (_M0L7_2abindS1364 == 0) {
    if (_M0L7_2abindS1367) {
      moonbit_decref(_M0L7_2abindS1367);
    }
    if (_M0L7_2abindS1364) {
      moonbit_decref(_M0L7_2abindS1364);
    }
    goto join_1346;
  } else {
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2aSomeS1369 =
      _M0L7_2abindS1364;
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2apath1S1370 =
      _M0L7_2aSomeS1369;
    if (_M0L7_2abindS1367 == 0) {
      moonbit_decref(_M0L8_2apath1S1370);
      if (_M0L7_2abindS1367) {
        moonbit_decref(_M0L7_2abindS1367);
      }
      goto join_1346;
    } else {
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2aSomeS1371 =
        _M0L7_2abindS1367;
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2apath2S1372 =
        _M0L7_2aSomeS1371;
      _M0L5path1S1348 = _M0L8_2apath1S1370;
      _M0L5path2S1349 = _M0L8_2apath2S1372;
      goto join_1347;
    }
  }
  join_1347:;
  _M0L8lca__idxS1350
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L8lca__idxS1350)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L8lca__idxS1350->$0 = 0;
  #line 233 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3284
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5path1S1348);
  #line 233 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3285
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5path2S1349);
  if (_M0L6_2atmpS3284 < _M0L6_2atmpS3285) {
    #line 234 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L8min__lenS1351
    = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5path1S1348);
  } else {
    #line 236 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L8min__lenS1351
    = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5path2S1349);
  }
  _M0L7_2abindS1352 = 0;
  _M0L1iS1353 = _M0L7_2abindS1352;
  while (1) {
    if (_M0L1iS1353 < _M0L8min__lenS1351) {
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3267;
      moonbit_string_t _M0L8_2afieldS3513;
      int32_t _M0L6_2acntS3704;
      moonbit_string_t _M0L4nameS3264;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3266;
      moonbit_string_t _M0L8_2afieldS3512;
      int32_t _M0L6_2acntS3709;
      moonbit_string_t _M0L4nameS3265;
      int32_t _result_3811;
      int32_t _M0L6_2atmpS3268;
      #line 239 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3267
      = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L5path1S1348, _M0L1iS1353);
      _M0L8_2afieldS3513 = _M0L6_2atmpS3267->$0;
      _M0L6_2acntS3704
      = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3267));
      if (_M0L6_2acntS3704 > 1) {
        int32_t _M0L11_2anew__cntS3708 = _M0L6_2acntS3704 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3267), _M0L11_2anew__cntS3708);
        if (_M0L8_2afieldS3513) {
          moonbit_incref(_M0L8_2afieldS3513);
        }
      } else if (_M0L6_2acntS3704 == 1) {
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2afieldS3707 =
          _M0L6_2atmpS3267->$3;
        void* _M0L8_2afieldS3706;
        void* _M0L8_2afieldS3705;
        moonbit_decref(_M0L8_2afieldS3707);
        _M0L8_2afieldS3706 = _M0L6_2atmpS3267->$2;
        moonbit_decref(_M0L8_2afieldS3706);
        _M0L8_2afieldS3705 = _M0L6_2atmpS3267->$1;
        moonbit_decref(_M0L8_2afieldS3705);
        #line 239 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L6_2atmpS3267);
      }
      _M0L4nameS3264 = _M0L8_2afieldS3513;
      #line 239 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3266
      = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L5path2S1349, _M0L1iS1353);
      _M0L8_2afieldS3512 = _M0L6_2atmpS3266->$0;
      _M0L6_2acntS3709
      = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3266));
      if (_M0L6_2acntS3709 > 1) {
        int32_t _M0L11_2anew__cntS3713 = _M0L6_2acntS3709 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3266), _M0L11_2anew__cntS3713);
        if (_M0L8_2afieldS3512) {
          moonbit_incref(_M0L8_2afieldS3512);
        }
      } else if (_M0L6_2acntS3709 == 1) {
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2afieldS3712 =
          _M0L6_2atmpS3266->$3;
        void* _M0L8_2afieldS3711;
        void* _M0L8_2afieldS3710;
        moonbit_decref(_M0L8_2afieldS3712);
        _M0L8_2afieldS3711 = _M0L6_2atmpS3266->$2;
        moonbit_decref(_M0L8_2afieldS3711);
        _M0L8_2afieldS3710 = _M0L6_2atmpS3266->$1;
        moonbit_decref(_M0L8_2afieldS3710);
        #line 239 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L6_2atmpS3266);
      }
      _M0L4nameS3265 = _M0L8_2afieldS3512;
      #line 239 "/home/developer/Documents/1/skbio_tree.mbt"
      _result_3811
      = _M0IPC16option6OptionPB2Eq5equalGsE(_M0L4nameS3264, _M0L4nameS3265);
      if (_M0L4nameS3264) {
        moonbit_decref(_M0L4nameS3264);
      }
      if (_M0L4nameS3265) {
        moonbit_decref(_M0L4nameS3265);
      }
      if (_result_3811) {
        _M0L8lca__idxS1350->$0 = _M0L1iS1353;
      } else {
        break;
      }
      _M0L6_2atmpS3268 = _M0L1iS1353 + 1;
      _M0L1iS1353 = _M0L6_2atmpS3268;
      continue;
    }
    break;
  }
  _M0L4distS1355
  = (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
  Moonbit_object_header(_M0L4distS1355)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L4distS1355->$0 = 0x0p+0;
  _M0L3valS3275 = _M0L8lca__idxS1350->$0;
  _M0L7_2abindS1356 = _M0L3valS3275 + 1;
  #line 247 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7_2abindS1357
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5path1S1348);
  _M0L1iS1358 = _M0L7_2abindS1356;
  while (1) {
    if (_M0L1iS1358 < _M0L7_2abindS1357) {
      double _M0L3valS3270 = _M0L4distS1355->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3273;
      void* _M0L8_2afieldS3511;
      int32_t _M0L6_2acntS3714;
      void* _M0L6lengthS3272;
      double _M0L6_2atmpS3271;
      double _M0L6_2atmpS3269;
      int32_t _M0L6_2atmpS3274;
      #line 248 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3273
      = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L5path1S1348, _M0L1iS1358);
      _M0L8_2afieldS3511 = _M0L6_2atmpS3273->$1;
      _M0L6_2acntS3714
      = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3273));
      if (_M0L6_2acntS3714 > 1) {
        int32_t _M0L11_2anew__cntS3718 = _M0L6_2acntS3714 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3273), _M0L11_2anew__cntS3718);
        moonbit_incref(_M0L8_2afieldS3511);
      } else if (_M0L6_2acntS3714 == 1) {
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2afieldS3717 =
          _M0L6_2atmpS3273->$3;
        void* _M0L8_2afieldS3716;
        moonbit_string_t _M0L8_2afieldS3715;
        moonbit_decref(_M0L8_2afieldS3717);
        _M0L8_2afieldS3716 = _M0L6_2atmpS3273->$2;
        moonbit_decref(_M0L8_2afieldS3716);
        _M0L8_2afieldS3715 = _M0L6_2atmpS3273->$0;
        if (_M0L8_2afieldS3715) {
          moonbit_decref(_M0L8_2afieldS3715);
        }
        #line 248 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L6_2atmpS3273);
      }
      _M0L6lengthS3272 = _M0L8_2afieldS3511;
      #line 248 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3271
      = _M0MPC16option6Option10unwrap__orGdE(_M0L6lengthS3272, 0x0p+0);
      moonbit_decref(_M0L6lengthS3272);
      _M0L6_2atmpS3269 = _M0L3valS3270 + _M0L6_2atmpS3271;
      _M0L4distS1355->$0 = _M0L6_2atmpS3269;
      _M0L6_2atmpS3274 = _M0L1iS1358 + 1;
      _M0L1iS1358 = _M0L6_2atmpS3274;
      continue;
    } else {
      moonbit_decref(_M0L5path1S1348);
    }
    break;
  }
  _M0L3valS3282 = _M0L8lca__idxS1350->$0;
  moonbit_decref(_M0L8lca__idxS1350);
  _M0L7_2abindS1360 = _M0L3valS3282 + 1;
  #line 250 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7_2abindS1361
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5path2S1349);
  _M0L1iS1362 = _M0L7_2abindS1360;
  while (1) {
    if (_M0L1iS1362 < _M0L7_2abindS1361) {
      double _M0L3valS3277 = _M0L4distS1355->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3280;
      void* _M0L8_2afieldS3510;
      int32_t _M0L6_2acntS3719;
      void* _M0L6lengthS3279;
      double _M0L6_2atmpS3278;
      double _M0L6_2atmpS3276;
      int32_t _M0L6_2atmpS3281;
      #line 251 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3280
      = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L5path2S1349, _M0L1iS1362);
      _M0L8_2afieldS3510 = _M0L6_2atmpS3280->$1;
      _M0L6_2acntS3719
      = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3280));
      if (_M0L6_2acntS3719 > 1) {
        int32_t _M0L11_2anew__cntS3723 = _M0L6_2acntS3719 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3280), _M0L11_2anew__cntS3723);
        moonbit_incref(_M0L8_2afieldS3510);
      } else if (_M0L6_2acntS3719 == 1) {
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2afieldS3722 =
          _M0L6_2atmpS3280->$3;
        void* _M0L8_2afieldS3721;
        moonbit_string_t _M0L8_2afieldS3720;
        moonbit_decref(_M0L8_2afieldS3722);
        _M0L8_2afieldS3721 = _M0L6_2atmpS3280->$2;
        moonbit_decref(_M0L8_2afieldS3721);
        _M0L8_2afieldS3720 = _M0L6_2atmpS3280->$0;
        if (_M0L8_2afieldS3720) {
          moonbit_decref(_M0L8_2afieldS3720);
        }
        #line 251 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L6_2atmpS3280);
      }
      _M0L6lengthS3279 = _M0L8_2afieldS3510;
      #line 251 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3278
      = _M0MPC16option6Option10unwrap__orGdE(_M0L6lengthS3279, 0x0p+0);
      moonbit_decref(_M0L6lengthS3279);
      _M0L6_2atmpS3276 = _M0L3valS3277 + _M0L6_2atmpS3278;
      _M0L4distS1355->$0 = _M0L6_2atmpS3276;
      _M0L6_2atmpS3281 = _M0L1iS1362 + 1;
      _M0L1iS1362 = _M0L6_2atmpS3281;
      continue;
    } else {
      moonbit_decref(_M0L5path2S1349);
    }
    break;
  }
  _M0L3valS3283 = _M0L4distS1355->$0;
  moonbit_decref(_M0L4distS1355);
  _block_3814
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGdE4Some));
  Moonbit_object_header(_block_3814)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 1);
  ((struct _M0DTPC16option6OptionGdE4Some*)_block_3814)->$0 = _M0L3valS3283;
  return _block_3814;
  join_1346:;
  return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode10find__path(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1344,
  moonbit_string_t _M0L4nameS1345
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4pathS1343;
  #line 117 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 118 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L4pathS1343
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
  #line 119 "/home/developer/Documents/1/skbio_tree.mbt"
  if (
    _M0FP26biolab8bio__seq21find__path__backtrack(_M0L4selfS1344, _M0L4nameS1345, _M0L4pathS1343)
  ) {
    return _M0L4pathS1343;
  } else {
    moonbit_decref(_M0L4pathS1343);
    return 0;
  }
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__11(
  
) {
  moonbit_string_t _M0L6newickS1339;
  struct moonbit_result_1 _tmp_3815;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1340;
  struct _M0TPB5ArrayGsE* _M0L5namesS1341;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L8ancestorS1342;
  int32_t _M0L6_2atmpS3261;
  int32_t _M0L6_2atmpS3259;
  void* _M0L4NoneS3260;
  struct moonbit_result_0 _result_3817;
  #line 161 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 162 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1339 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 163 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3815 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1339);
  moonbit_decref(_M0L6newickS1339);
  if (_tmp_3815.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3262 =
      _tmp_3815.data.ok;
    _M0L4treeS1340 = _M0L5_2aokS3262;
  } else {
    void* const _M0L6_2aerrS3263 = _tmp_3815.data.err;
    struct moonbit_result_0 _result_3816;
    _result_3816.tag = 0;
    _result_3816.data.err = _M0L6_2aerrS3263;
    return _result_3816;
  }
  #line 164 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5namesS1341 = _M0MPC15array5Array11new_2einnerGsE(0);
  #line 165 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0MPC15array5Array4pushGsE(_M0L5namesS1341, (moonbit_string_t)moonbit_string_literal_10.data);
  #line 166 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0MPC15array5Array4pushGsE(_M0L5namesS1341, (moonbit_string_t)moonbit_string_literal_20.data);
  #line 167 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L8ancestorS1342
  = _M0MP26biolab8bio__seq8TreeNode3lca(_M0L4treeS1340, _M0L5namesS1341);
  moonbit_decref(_M0L4treeS1340);
  moonbit_decref(_M0L5namesS1341);
  _M0L6_2atmpS3261 = _M0L8ancestorS1342 == 0;
  if (_M0L8ancestorS1342) {
    moonbit_decref(_M0L8ancestorS1342);
  }
  _M0L6_2atmpS3259 = !_M0L6_2atmpS3261;
  _M0L4NoneS3260
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 168 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3817
  = _M0FPB12assert__true(_M0L6_2atmpS3259, _M0L4NoneS3260, (moonbit_string_t)moonbit_string_literal_21.data);
  moonbit_decref(_M0L4NoneS3260);
  return _result_3817;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq54____test__736b62696f5f747265655f62656e63682e6d6274__10(
  
) {
  moonbit_string_t _M0L6newickS1335;
  struct moonbit_result_1 _tmp_3818;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1336;
  struct _M0TPB5ArrayGsE* _M0L5namesS1337;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L8ancestorS1338;
  int32_t _M0L6_2atmpS3256;
  int32_t _M0L6_2atmpS3254;
  void* _M0L4NoneS3255;
  struct moonbit_result_0 _result_3820;
  #line 150 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 151 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1335 = _M0FP26biolab8bio__seq21gen__balanced__newick(10);
  #line 152 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3818 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1335);
  moonbit_decref(_M0L6newickS1335);
  if (_tmp_3818.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3257 =
      _tmp_3818.data.ok;
    _M0L4treeS1336 = _M0L5_2aokS3257;
  } else {
    void* const _M0L6_2aerrS3258 = _tmp_3818.data.err;
    struct moonbit_result_0 _result_3819;
    _result_3819.tag = 0;
    _result_3819.data.err = _M0L6_2aerrS3258;
    return _result_3819;
  }
  #line 153 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5namesS1337 = _M0MPC15array5Array11new_2einnerGsE(0);
  #line 154 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0MPC15array5Array4pushGsE(_M0L5namesS1337, (moonbit_string_t)moonbit_string_literal_17.data);
  #line 155 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0MPC15array5Array4pushGsE(_M0L5namesS1337, (moonbit_string_t)moonbit_string_literal_18.data);
  #line 156 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L8ancestorS1338
  = _M0MP26biolab8bio__seq8TreeNode3lca(_M0L4treeS1336, _M0L5namesS1337);
  moonbit_decref(_M0L4treeS1336);
  moonbit_decref(_M0L5namesS1337);
  _M0L6_2atmpS3256 = _M0L8ancestorS1338 == 0;
  if (_M0L8ancestorS1338) {
    moonbit_decref(_M0L8ancestorS1338);
  }
  _M0L6_2atmpS3254 = !_M0L6_2atmpS3256;
  _M0L4NoneS3255
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 157 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3820
  = _M0FPB12assert__true(_M0L6_2atmpS3254, _M0L4NoneS3255, (moonbit_string_t)moonbit_string_literal_22.data);
  moonbit_decref(_M0L4NoneS3255);
  return _result_3820;
}

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode3lca(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1317,
  struct _M0TPB5ArrayGsE* _M0L5namesS1311
) {
  int32_t _M0L6_2atmpS3233;
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L5pathsS1312;
  int32_t _M0L7_2abindS1313;
  int32_t _M0L2__S1314;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3253;
  int32_t _M0L6_2atmpS3252;
  struct _M0TPB8MutLocalGiE* _M0L1mS1320;
  int32_t _M0L7_2abindS1321;
  int32_t _M0L7_2abindS1322;
  int32_t _M0L1iS1323;
  int32_t _M0L8min__lenS1319;
  struct _M0TPB8MutLocalGiE* _M0L8lca__idxS1325;
  int32_t _M0L7_2abindS1326;
  int32_t _M0L1iS1327;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3244;
  int32_t _M0L3valS3245;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3243;
  #line 46 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 47 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3233 = _M0MPC15array5Array6lengthGsE(_M0L5namesS1311);
  if (_M0L6_2atmpS3233 == 0) {
    return 0;
  }
  #line 52 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L5pathsS1312
  = _M0MPC15array5Array11new_2einnerGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(0);
  _M0L7_2abindS1313 = _M0L5namesS1311->$1;
  _M0L2__S1314 = 0;
  while (1) {
    if (_M0L2__S1314 < _M0L7_2abindS1313) {
      moonbit_string_t* _M0L3bufS3235 = _M0L5namesS1311->$0;
      moonbit_string_t _M0L4nameS1315 =
        (moonbit_string_t)_M0L3bufS3235[_M0L2__S1314];
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4pathS1316;
      int32_t _result_3822;
      int32_t _M0L6_2atmpS3234;
      moonbit_incref(_M0L4nameS1315);
      #line 54 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L4pathS1316
      = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
      #line 55 "/home/developer/Documents/1/skbio_tree.mbt"
      _result_3822
      = _M0FP26biolab8bio__seq21find__path__backtrack(_M0L4selfS1317, _M0L4nameS1315, _M0L4pathS1316);
      moonbit_decref(_M0L4nameS1315);
      if (_result_3822) {
        #line 56 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0MPC15array5Array4pushGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, _M0L4pathS1316);
        moonbit_decref(_M0L4pathS1316);
      } else {
        moonbit_decref(_M0L4pathS1316);
        moonbit_decref(_M0L5pathsS1312);
        return 0;
      }
      _M0L6_2atmpS3234 = _M0L2__S1314 + 1;
      _M0L2__S1314 = _M0L6_2atmpS3234;
      continue;
    }
    break;
  }
  #line 64 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3253
  = _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, 0);
  #line 64 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3252
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS3253);
  moonbit_decref(_M0L6_2atmpS3253);
  _M0L1mS1320
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1mS1320)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1mS1320->$0 = _M0L6_2atmpS3252;
  _M0L7_2abindS1321 = 1;
  #line 65 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7_2abindS1322
  = _M0MPC15array5Array6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312);
  _M0L1iS1323 = _M0L7_2abindS1321;
  while (1) {
    if (_M0L1iS1323 < _M0L7_2abindS1322) {
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3248;
      int32_t _M0L6_2atmpS3246;
      int32_t _M0L3valS3247;
      int32_t _M0L6_2atmpS3251;
      #line 66 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3248
      = _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, _M0L1iS1323);
      #line 66 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3246
      = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS3248);
      moonbit_decref(_M0L6_2atmpS3248);
      _M0L3valS3247 = _M0L1mS1320->$0;
      if (_M0L6_2atmpS3246 < _M0L3valS3247) {
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3250;
        int32_t _M0L6_2atmpS3249;
        #line 67 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L6_2atmpS3250
        = _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, _M0L1iS1323);
        #line 67 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L6_2atmpS3249
        = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS3250);
        moonbit_decref(_M0L6_2atmpS3250);
        _M0L1mS1320->$0 = _M0L6_2atmpS3249;
      }
      _M0L6_2atmpS3251 = _M0L1iS1323 + 1;
      _M0L1iS1323 = _M0L6_2atmpS3251;
      continue;
    }
    break;
  }
  _M0L8min__lenS1319 = _M0L1mS1320->$0;
  moonbit_decref(_M0L1mS1320);
  _M0L8lca__idxS1325
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L8lca__idxS1325)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L8lca__idxS1325->$0 = 0;
  _M0L7_2abindS1326 = 0;
  _M0L1iS1327 = _M0L7_2abindS1326;
  while (1) {
    if (_M0L1iS1327 < _M0L8min__lenS1319) {
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3241;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L9candidateS1328;
      struct _M0TPB8MutLocalGbE* _M0L10all__matchS1329;
      int32_t _M0L7_2abindS1330;
      int32_t _M0L7_2abindS1331;
      int32_t _M0L1jS1332;
      int32_t _result_3827;
      int32_t _M0L6_2atmpS3242;
      #line 75 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3241
      = _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, 0);
      #line 75 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L9candidateS1328
      = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS3241, _M0L1iS1327);
      moonbit_decref(_M0L6_2atmpS3241);
      _M0L10all__matchS1329
      = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
      Moonbit_object_header(_M0L10all__matchS1329)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L10all__matchS1329->$0 = 1;
      _M0L7_2abindS1330 = 1;
      #line 77 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L7_2abindS1331
      = _M0MPC15array5Array6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312);
      _M0L1jS1332 = _M0L7_2abindS1330;
      while (1) {
        if (_M0L1jS1332 < _M0L7_2abindS1331) {
          struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3239;
          struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3238;
          moonbit_string_t _M0L8_2afieldS3515;
          int32_t _M0L6_2acntS3724;
          moonbit_string_t _M0L4nameS3236;
          moonbit_string_t _M0L4nameS3237;
          int32_t _result_3826;
          int32_t _M0L6_2atmpS3240;
          #line 78 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0L6_2atmpS3239
          = _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, _M0L1jS1332);
          #line 78 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0L6_2atmpS3238
          = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS3239, _M0L1iS1327);
          moonbit_decref(_M0L6_2atmpS3239);
          _M0L8_2afieldS3515 = _M0L6_2atmpS3238->$0;
          _M0L6_2acntS3724
          = Moonbit_rc_count(Moonbit_object_header(_M0L6_2atmpS3238));
          if (_M0L6_2acntS3724 > 1) {
            int32_t _M0L11_2anew__cntS3728 = _M0L6_2acntS3724 - 1;
            Moonbit_set_rc_count(Moonbit_object_header(_M0L6_2atmpS3238), _M0L11_2anew__cntS3728);
            if (_M0L8_2afieldS3515) {
              moonbit_incref(_M0L8_2afieldS3515);
            }
          } else if (_M0L6_2acntS3724 == 1) {
            struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2afieldS3727 =
              _M0L6_2atmpS3238->$3;
            void* _M0L8_2afieldS3726;
            void* _M0L8_2afieldS3725;
            moonbit_decref(_M0L8_2afieldS3727);
            _M0L8_2afieldS3726 = _M0L6_2atmpS3238->$2;
            moonbit_decref(_M0L8_2afieldS3726);
            _M0L8_2afieldS3725 = _M0L6_2atmpS3238->$1;
            moonbit_decref(_M0L8_2afieldS3725);
            #line 78 "/home/developer/Documents/1/skbio_tree.mbt"
            moonbit_free(_M0L6_2atmpS3238);
          }
          _M0L4nameS3236 = _M0L8_2afieldS3515;
          _M0L4nameS3237 = _M0L9candidateS1328->$0;
          if (_M0L4nameS3237) {
            moonbit_incref(_M0L4nameS3237);
          }
          #line 78 "/home/developer/Documents/1/skbio_tree.mbt"
          _result_3826
          = _M0IP016_24default__implPB2Eq10not__equalGOsE(_M0L4nameS3236, _M0L4nameS3237);
          if (_M0L4nameS3236) {
            moonbit_decref(_M0L4nameS3236);
          }
          if (_M0L4nameS3237) {
            moonbit_decref(_M0L4nameS3237);
          }
          if (_result_3826) {
            moonbit_decref(_M0L9candidateS1328);
            _M0L10all__matchS1329->$0 = 0;
            break;
          }
          _M0L6_2atmpS3240 = _M0L1jS1332 + 1;
          _M0L1jS1332 = _M0L6_2atmpS3240;
          continue;
        } else {
          moonbit_decref(_M0L9candidateS1328);
        }
        break;
      }
      _result_3827 = _M0L10all__matchS1329->$0;
      moonbit_decref(_M0L10all__matchS1329);
      if (_result_3827) {
        _M0L8lca__idxS1325->$0 = _M0L1iS1327;
      } else {
        break;
      }
      _M0L6_2atmpS3242 = _M0L1iS1327 + 1;
      _M0L1iS1327 = _M0L6_2atmpS3242;
      continue;
    }
    break;
  }
  #line 90 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3244
  = _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L5pathsS1312, 0);
  moonbit_decref(_M0L5pathsS1312);
  _M0L3valS3245 = _M0L8lca__idxS1325->$0;
  moonbit_decref(_M0L8lca__idxS1325);
  #line 90 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3243
  = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS3244, _M0L3valS3245);
  moonbit_decref(_M0L6_2atmpS3244);
  return _M0L6_2atmpS3243;
}

int32_t _M0FP26biolab8bio__seq21find__path__backtrack(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4nodeS1304,
  moonbit_string_t _M0L4nameS1305,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4pathS1303
) {
  moonbit_string_t _M0L4nameS3226;
  moonbit_string_t _M0L6_2atmpS3227;
  int32_t _result_3828;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1306;
  int32_t _M0L7_2abindS1307;
  int32_t _M0L2__S1308;
  int32_t _M0L6_2atmpS3232;
  int32_t _M0L6_2atmpS3231;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3230;
  #line 96 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 101 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L4pathS1303, _M0L4nodeS1304);
  _M0L4nameS3226 = _M0L4nodeS1304->$0;
  _M0L6_2atmpS3227 = _M0L4nameS1305;
  if (_M0L4nameS3226) {
    moonbit_incref(_M0L4nameS3226);
  }
  #line 102 "/home/developer/Documents/1/skbio_tree.mbt"
  _result_3828
  = _M0IPC16option6OptionPB2Eq5equalGsE(_M0L4nameS3226, _M0L6_2atmpS3227);
  if (_M0L4nameS3226) {
    moonbit_decref(_M0L4nameS3226);
  }
  if (_result_3828) {
    return 1;
  }
  _M0L7_2abindS1306 = _M0L4nodeS1304->$3;
  _M0L7_2abindS1307 = _M0L7_2abindS1306->$1;
  moonbit_incref(_M0L7_2abindS1306);
  _M0L2__S1308 = 0;
  while (1) {
    if (_M0L2__S1308 < _M0L7_2abindS1307) {
      struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3229 =
        _M0L7_2abindS1306->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1309 =
        (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3229[_M0L2__S1308];
      int32_t _result_3830;
      int32_t _M0L6_2atmpS3228;
      moonbit_incref(_M0L5childS1309);
      #line 106 "/home/developer/Documents/1/skbio_tree.mbt"
      _result_3830
      = _M0FP26biolab8bio__seq21find__path__backtrack(_M0L5childS1309, _M0L4nameS1305, _M0L4pathS1303);
      moonbit_decref(_M0L5childS1309);
      if (_result_3830) {
        moonbit_decref(_M0L7_2abindS1306);
        return 1;
      }
      _M0L6_2atmpS3228 = _M0L2__S1308 + 1;
      _M0L2__S1308 = _M0L6_2atmpS3228;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1306);
    }
    break;
  }
  #line 110 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3232
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L4pathS1303);
  _M0L6_2atmpS3231 = _M0L6_2atmpS3232 - 1;
  #line 110 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3230
  = _M0MPC15array5Array6removeGRP26biolab8bio__seq8TreeNodeE(_M0L4pathS1303, _M0L6_2atmpS3231);
  moonbit_decref(_M0L6_2atmpS3230);
  return 0;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__9(
  
) {
  moonbit_string_t _M0L6newickS1300;
  struct moonbit_result_1 _tmp_3831;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1301;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L5foundS1302;
  int32_t _M0L6_2atmpS3223;
  int32_t _M0L6_2atmpS3221;
  void* _M0L4NoneS3222;
  struct moonbit_result_0 _result_3833;
  #line 142 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 143 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1300 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 144 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3831 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1300);
  moonbit_decref(_M0L6newickS1300);
  if (_tmp_3831.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3224 =
      _tmp_3831.data.ok;
    _M0L4treeS1301 = _M0L5_2aokS3224;
  } else {
    void* const _M0L6_2aerrS3225 = _tmp_3831.data.err;
    struct moonbit_result_0 _result_3832;
    _result_3832.tag = 0;
    _result_3832.data.err = _M0L6_2aerrS3225;
    return _result_3832;
  }
  #line 145 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5foundS1302
  = _M0MP26biolab8bio__seq8TreeNode4find(_M0L4treeS1301, (moonbit_string_t)moonbit_string_literal_20.data);
  moonbit_decref(_M0L4treeS1301);
  _M0L6_2atmpS3223 = _M0L5foundS1302 == 0;
  if (_M0L5foundS1302) {
    moonbit_decref(_M0L5foundS1302);
  }
  _M0L6_2atmpS3221 = !_M0L6_2atmpS3223;
  _M0L4NoneS3222
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 146 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3833
  = _M0FPB12assert__true(_M0L6_2atmpS3221, _M0L4NoneS3222, (moonbit_string_t)moonbit_string_literal_23.data);
  moonbit_decref(_M0L4NoneS3222);
  return _result_3833;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__8(
  
) {
  moonbit_string_t _M0L6newickS1297;
  struct moonbit_result_1 _tmp_3834;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1298;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L5foundS1299;
  int32_t _M0L6_2atmpS3218;
  int32_t _M0L6_2atmpS3216;
  void* _M0L4NoneS3217;
  struct moonbit_result_0 _result_3836;
  #line 134 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 135 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1297 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 136 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3834 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1297);
  moonbit_decref(_M0L6newickS1297);
  if (_tmp_3834.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3219 =
      _tmp_3834.data.ok;
    _M0L4treeS1298 = _M0L5_2aokS3219;
  } else {
    void* const _M0L6_2aerrS3220 = _tmp_3834.data.err;
    struct moonbit_result_0 _result_3835;
    _result_3835.tag = 0;
    _result_3835.data.err = _M0L6_2aerrS3220;
    return _result_3835;
  }
  #line 137 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5foundS1299
  = _M0MP26biolab8bio__seq8TreeNode4find(_M0L4treeS1298, (moonbit_string_t)moonbit_string_literal_24.data);
  moonbit_decref(_M0L4treeS1298);
  _M0L6_2atmpS3218 = _M0L5foundS1299 == 0;
  if (_M0L5foundS1299) {
    moonbit_decref(_M0L5foundS1299);
  }
  _M0L6_2atmpS3216 = !_M0L6_2atmpS3218;
  _M0L4NoneS3217
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 138 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3836
  = _M0FPB12assert__true(_M0L6_2atmpS3216, _M0L4NoneS3217, (moonbit_string_t)moonbit_string_literal_25.data);
  moonbit_decref(_M0L4NoneS3217);
  return _result_3836;
}

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode4find(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1285,
  moonbit_string_t _M0L4nameS1286
) {
  moonbit_string_t _M0L4nameS3212;
  moonbit_string_t _M0L6_2atmpS3213;
  int32_t _result_3837;
  #line 271 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L4nameS3212 = _M0L4selfS1285->$0;
  _M0L6_2atmpS3213 = _M0L4nameS1286;
  if (_M0L4nameS3212) {
    moonbit_incref(_M0L4nameS3212);
  }
  #line 272 "/home/developer/Documents/1/skbio_tree.mbt"
  _result_3837
  = _M0IPC16option6OptionPB2Eq5equalGsE(_M0L4nameS3212, _M0L6_2atmpS3213);
  if (_M0L4nameS3212) {
    moonbit_decref(_M0L4nameS3212);
  }
  if (_result_3837) {
    moonbit_incref(_M0L4selfS1285);
    return _M0L4selfS1285;
  } else {
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1287 =
      _M0L4selfS1285->$3;
    int32_t _M0L7_2abindS1288 = _M0L7_2abindS1287->$1;
    int32_t _M0L2__S1289;
    moonbit_incref(_M0L7_2abindS1287);
    _M0L2__S1289 = 0;
    while (1) {
      if (_M0L2__S1289 < _M0L7_2abindS1288) {
        struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3215 =
          _M0L7_2abindS1287->$0;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1290 =
          (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3215[
            _M0L2__S1289
          ];
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L5foundS1292;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L7_2abindS1293;
        int32_t _M0L6_2atmpS3214;
        moonbit_incref(_M0L5childS1290);
        #line 276 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L7_2abindS1293
        = _M0MP26biolab8bio__seq8TreeNode4find(_M0L5childS1290, _M0L4nameS1286);
        moonbit_decref(_M0L5childS1290);
        if (_M0L7_2abindS1293 == 0) {
          if (_M0L7_2abindS1293) {
            moonbit_decref(_M0L7_2abindS1293);
          }
        } else {
          struct _M0TP26biolab8bio__seq8TreeNode* _M0L7_2aSomeS1294;
          struct _M0TP26biolab8bio__seq8TreeNode* _M0L8_2afoundS1295;
          moonbit_decref(_M0L7_2abindS1287);
          _M0L7_2aSomeS1294 = _M0L7_2abindS1293;
          _M0L8_2afoundS1295 = _M0L7_2aSomeS1294;
          _M0L5foundS1292 = _M0L8_2afoundS1295;
          goto join_1291;
        }
        goto joinlet_3839;
        join_1291:;
        return _M0L5foundS1292;
        joinlet_3839:;
        _M0L6_2atmpS3214 = _M0L2__S1289 + 1;
        _M0L2__S1289 = _M0L6_2atmpS3214;
        continue;
      } else {
        moonbit_decref(_M0L7_2abindS1287);
      }
      break;
    }
    return 0;
  }
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__7(
  
) {
  moonbit_string_t _M0L6newickS1282;
  struct moonbit_result_1 _tmp_3840;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1283;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5nodesS1284;
  int32_t _M0L6_2atmpS3208;
  moonbit_string_t _M0L6_2atmpS3209;
  struct moonbit_result_0 _result_3842;
  #line 126 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 127 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1282 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 128 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3840 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1282);
  moonbit_decref(_M0L6newickS1282);
  if (_tmp_3840.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3210 =
      _tmp_3840.data.ok;
    _M0L4treeS1283 = _M0L5_2aokS3210;
  } else {
    void* const _M0L6_2aerrS3211 = _tmp_3840.data.err;
    struct moonbit_result_0 _result_3841;
    _result_3841.tag = 0;
    _result_3841.data.err = _M0L6_2aerrS3211;
    return _result_3841;
  }
  #line 129 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5nodesS1284
  = _M0MP26biolab8bio__seq8TreeNode10levelorder(_M0L4treeS1283);
  moonbit_decref(_M0L4treeS1283);
  #line 130 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3208
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5nodesS1284);
  moonbit_decref(_M0L5nodesS1284);
  _M0L6_2atmpS3209 = 0;
  #line 130 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3842
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3208, 8191, _M0L6_2atmpS3209, (moonbit_string_t)moonbit_string_literal_26.data);
  if (_M0L6_2atmpS3209) {
    moonbit_decref(_M0L6_2atmpS3209);
  }
  return _result_3842;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode10levelorder(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1273
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6resultS1271;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5queueS1272;
  struct _M0TPB8MutLocalGiE* _M0L4headS1274;
  #line 164 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 165 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6resultS1271
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
  #line 166 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L5queueS1272
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
  #line 167 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L5queueS1272, _M0L4selfS1273);
  _M0L4headS1274
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L4headS1274)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L4headS1274->$0 = 0;
  while (1) {
    int32_t _M0L3valS3201 = _M0L4headS1274->$0;
    int32_t _M0L6_2atmpS3202;
    #line 169 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3202
    = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5queueS1272);
    if (_M0L3valS3201 < _M0L6_2atmpS3202) {
      int32_t _M0L3valS3207 = _M0L4headS1274->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L4nodeS1275;
      int32_t _M0L3valS3204;
      int32_t _M0L6_2atmpS3203;
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8_2afieldS3528;
      int32_t _M0L6_2acntS3729;
      struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1276;
      int32_t _M0L7_2abindS1277;
      int32_t _M0L2__S1278;
      #line 170 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L4nodeS1275
      = _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(_M0L5queueS1272, _M0L3valS3207);
      _M0L3valS3204 = _M0L4headS1274->$0;
      _M0L6_2atmpS3203 = _M0L3valS3204 + 1;
      _M0L4headS1274->$0 = _M0L6_2atmpS3203;
      #line 172 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L6resultS1271, _M0L4nodeS1275);
      _M0L8_2afieldS3528 = _M0L4nodeS1275->$3;
      _M0L6_2acntS3729
      = Moonbit_rc_count(Moonbit_object_header(_M0L4nodeS1275));
      if (_M0L6_2acntS3729 > 1) {
        int32_t _M0L11_2anew__cntS3733 = _M0L6_2acntS3729 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L4nodeS1275), _M0L11_2anew__cntS3733);
        moonbit_incref(_M0L8_2afieldS3528);
      } else if (_M0L6_2acntS3729 == 1) {
        void* _M0L8_2afieldS3732 = _M0L4nodeS1275->$2;
        void* _M0L8_2afieldS3731;
        moonbit_string_t _M0L8_2afieldS3730;
        moonbit_decref(_M0L8_2afieldS3732);
        _M0L8_2afieldS3731 = _M0L4nodeS1275->$1;
        moonbit_decref(_M0L8_2afieldS3731);
        _M0L8_2afieldS3730 = _M0L4nodeS1275->$0;
        if (_M0L8_2afieldS3730) {
          moonbit_decref(_M0L8_2afieldS3730);
        }
        #line 173 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L4nodeS1275);
      }
      _M0L7_2abindS1276 = _M0L8_2afieldS3528;
      _M0L7_2abindS1277 = _M0L7_2abindS1276->$1;
      _M0L2__S1278 = 0;
      while (1) {
        if (_M0L2__S1278 < _M0L7_2abindS1277) {
          struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3206 =
            _M0L7_2abindS1276->$0;
          struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1279 =
            (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3206[
              _M0L2__S1278
            ];
          int32_t _M0L6_2atmpS3205;
          moonbit_incref(_M0L5childS1279);
          #line 174 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L5queueS1272, _M0L5childS1279);
          moonbit_decref(_M0L5childS1279);
          _M0L6_2atmpS3205 = _M0L2__S1278 + 1;
          _M0L2__S1278 = _M0L6_2atmpS3205;
          continue;
        } else {
          moonbit_decref(_M0L7_2abindS1276);
        }
        break;
      }
      continue;
    } else {
      moonbit_decref(_M0L4headS1274);
      moonbit_decref(_M0L5queueS1272);
    }
    break;
  }
  return _M0L6resultS1271;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__6(
  
) {
  moonbit_string_t _M0L6newickS1268;
  struct moonbit_result_1 _tmp_3845;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1269;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5nodesS1270;
  int32_t _M0L6_2atmpS3197;
  moonbit_string_t _M0L6_2atmpS3198;
  struct moonbit_result_0 _result_3847;
  #line 118 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 119 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1268 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 120 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3845 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1268);
  moonbit_decref(_M0L6newickS1268);
  if (_tmp_3845.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3199 =
      _tmp_3845.data.ok;
    _M0L4treeS1269 = _M0L5_2aokS3199;
  } else {
    void* const _M0L6_2aerrS3200 = _tmp_3845.data.err;
    struct moonbit_result_0 _result_3846;
    _result_3846.tag = 0;
    _result_3846.data.err = _M0L6_2aerrS3200;
    return _result_3846;
  }
  #line 121 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5nodesS1270 = _M0MP26biolab8bio__seq8TreeNode9postorder(_M0L4treeS1269);
  moonbit_decref(_M0L4treeS1269);
  #line 122 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3197
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5nodesS1270);
  moonbit_decref(_M0L5nodesS1270);
  _M0L6_2atmpS3198 = 0;
  #line 122 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3847
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3197, 8191, _M0L6_2atmpS3198, (moonbit_string_t)moonbit_string_literal_27.data);
  if (_M0L6_2atmpS3198) {
    moonbit_decref(_M0L6_2atmpS3198);
  }
  return _result_3847;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode9postorder(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1267
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6resultS1266;
  #line 144 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 145 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6resultS1266
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
  #line 146 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MP26biolab8bio__seq8TreeNode17postorder__helper(_M0L4selfS1267, _M0L6resultS1266);
  return _M0L6resultS1266;
}

int32_t _M0MP26biolab8bio__seq8TreeNode17postorder__helper(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1260,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6resultS1264
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1259;
  int32_t _M0L7_2abindS1261;
  int32_t _M0L2__S1262;
  #line 151 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7_2abindS1259 = _M0L4selfS1260->$3;
  _M0L7_2abindS1261 = _M0L7_2abindS1259->$1;
  moonbit_incref(_M0L7_2abindS1259);
  _M0L2__S1262 = 0;
  while (1) {
    if (_M0L2__S1262 < _M0L7_2abindS1261) {
      struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3196 =
        _M0L7_2abindS1259->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1263 =
        (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3196[_M0L2__S1262];
      int32_t _M0L6_2atmpS3195;
      moonbit_incref(_M0L5childS1263);
      #line 156 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0MP26biolab8bio__seq8TreeNode17postorder__helper(_M0L5childS1263, _M0L6resultS1264);
      moonbit_decref(_M0L5childS1263);
      _M0L6_2atmpS3195 = _M0L2__S1262 + 1;
      _M0L2__S1262 = _M0L6_2atmpS3195;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1259);
    }
    break;
  }
  #line 158 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L6resultS1264, _M0L4selfS1260);
  return 0;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__5(
  
) {
  moonbit_string_t _M0L6newickS1256;
  struct moonbit_result_1 _tmp_3849;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1257;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5nodesS1258;
  int32_t _M0L6_2atmpS3191;
  moonbit_string_t _M0L6_2atmpS3192;
  struct moonbit_result_0 _result_3851;
  #line 110 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 111 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1256 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 112 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3849 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1256);
  moonbit_decref(_M0L6newickS1256);
  if (_tmp_3849.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3193 =
      _tmp_3849.data.ok;
    _M0L4treeS1257 = _M0L5_2aokS3193;
  } else {
    void* const _M0L6_2aerrS3194 = _tmp_3849.data.err;
    struct moonbit_result_0 _result_3850;
    _result_3850.tag = 0;
    _result_3850.data.err = _M0L6_2aerrS3194;
    return _result_3850;
  }
  #line 113 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L5nodesS1258 = _M0MP26biolab8bio__seq8TreeNode8preorder(_M0L4treeS1257);
  moonbit_decref(_M0L4treeS1257);
  #line 114 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3191
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L5nodesS1258);
  moonbit_decref(_M0L5nodesS1258);
  _M0L6_2atmpS3192 = 0;
  #line 114 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3851
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3191, 8191, _M0L6_2atmpS3192, (moonbit_string_t)moonbit_string_literal_28.data);
  if (_M0L6_2atmpS3192) {
    moonbit_decref(_M0L6_2atmpS3192);
  }
  return _result_3851;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MP26biolab8bio__seq8TreeNode8preorder(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1255
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6resultS1254;
  #line 128 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 129 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6resultS1254
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
  #line 130 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MP26biolab8bio__seq8TreeNode16preorder__helper(_M0L4selfS1255, _M0L6resultS1254);
  return _M0L6resultS1254;
}

int32_t _M0MP26biolab8bio__seq8TreeNode16preorder__helper(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1248,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6resultS1247
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1249;
  int32_t _M0L7_2abindS1250;
  int32_t _M0L2__S1251;
  #line 135 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 136 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L6resultS1247, _M0L4selfS1248);
  _M0L7_2abindS1249 = _M0L4selfS1248->$3;
  _M0L7_2abindS1250 = _M0L7_2abindS1249->$1;
  moonbit_incref(_M0L7_2abindS1249);
  _M0L2__S1251 = 0;
  while (1) {
    if (_M0L2__S1251 < _M0L7_2abindS1250) {
      struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3190 =
        _M0L7_2abindS1249->$0;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1252 =
        (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3190[_M0L2__S1251];
      int32_t _M0L6_2atmpS3189;
      moonbit_incref(_M0L5childS1252);
      #line 138 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0MP26biolab8bio__seq8TreeNode16preorder__helper(_M0L5childS1252, _M0L6resultS1247);
      moonbit_decref(_M0L5childS1252);
      _M0L6_2atmpS3189 = _M0L2__S1251 + 1;
      _M0L2__S1251 = _M0L6_2atmpS3189;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1249);
    }
    break;
  }
  return 0;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__4(
  
) {
  moonbit_string_t _M0L6newickS1243;
  struct moonbit_result_1 _tmp_3853;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1244;
  moonbit_string_t _M0L6outputS1245;
  moonbit_string_t _M0L7_2abindS1246;
  int32_t _M0L6_2atmpS3186;
  struct _M0TPC16string10StringView _M0L6_2atmpS3185;
  int32_t _M0L6_2atmpS3183;
  void* _M0L4NoneS3184;
  struct moonbit_result_0 _result_3855;
  #line 102 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 103 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1243 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 104 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3853 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1243);
  moonbit_decref(_M0L6newickS1243);
  if (_tmp_3853.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3187 =
      _tmp_3853.data.ok;
    _M0L4treeS1244 = _M0L5_2aokS3187;
  } else {
    void* const _M0L6_2aerrS3188 = _tmp_3853.data.err;
    struct moonbit_result_0 _result_3854;
    _result_3854.tag = 0;
    _result_3854.data.err = _M0L6_2aerrS3188;
    return _result_3854;
  }
  #line 105 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6outputS1245
  = _M0MP26biolab8bio__seq8TreeNode10to__newick(_M0L4treeS1244);
  moonbit_decref(_M0L4treeS1244);
  _M0L7_2abindS1246 = (moonbit_string_t)moonbit_string_literal_29.data;
  _M0L6_2atmpS3186 = Moonbit_array_length(_M0L7_2abindS1246);
  _M0L6_2atmpS3185
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1246, .$1 = 0, .$2 = _M0L6_2atmpS3186
  };
  #line 106 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3183
  = _M0MPC16string6String11has__suffix(_M0L6outputS1245, _M0L6_2atmpS3185);
  moonbit_decref(_M0L6outputS1245);
  moonbit_decref(_M0L6_2atmpS3185.$0);
  _M0L4NoneS3184
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 106 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3855
  = _M0FPB12assert__true(_M0L6_2atmpS3183, _M0L4NoneS3184, (moonbit_string_t)moonbit_string_literal_30.data);
  moonbit_decref(_M0L4NoneS3184);
  return _result_3855;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__3(
  
) {
  moonbit_string_t _M0L6newickS1239;
  struct moonbit_result_1 _tmp_3856;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1240;
  moonbit_string_t _M0L6outputS1241;
  moonbit_string_t _M0L7_2abindS1242;
  int32_t _M0L6_2atmpS3180;
  struct _M0TPC16string10StringView _M0L6_2atmpS3179;
  int32_t _M0L6_2atmpS3177;
  void* _M0L4NoneS3178;
  struct moonbit_result_0 _result_3858;
  #line 94 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 95 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1239 = _M0FP26biolab8bio__seq21gen__balanced__newick(8);
  #line 96 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3856 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1239);
  moonbit_decref(_M0L6newickS1239);
  if (_tmp_3856.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3181 =
      _tmp_3856.data.ok;
    _M0L4treeS1240 = _M0L5_2aokS3181;
  } else {
    void* const _M0L6_2aerrS3182 = _tmp_3856.data.err;
    struct moonbit_result_0 _result_3857;
    _result_3857.tag = 0;
    _result_3857.data.err = _M0L6_2aerrS3182;
    return _result_3857;
  }
  #line 97 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6outputS1241
  = _M0MP26biolab8bio__seq8TreeNode10to__newick(_M0L4treeS1240);
  moonbit_decref(_M0L4treeS1240);
  _M0L7_2abindS1242 = (moonbit_string_t)moonbit_string_literal_29.data;
  _M0L6_2atmpS3180 = Moonbit_array_length(_M0L7_2abindS1242);
  _M0L6_2atmpS3179
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1242, .$1 = 0, .$2 = _M0L6_2atmpS3180
  };
  #line 98 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3177
  = _M0MPC16string6String11has__suffix(_M0L6outputS1241, _M0L6_2atmpS3179);
  moonbit_decref(_M0L6outputS1241);
  moonbit_decref(_M0L6_2atmpS3179.$0);
  _M0L4NoneS3178
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 98 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3858
  = _M0FPB12assert__true(_M0L6_2atmpS3177, _M0L4NoneS3178, (moonbit_string_t)moonbit_string_literal_31.data);
  moonbit_decref(_M0L4NoneS3178);
  return _result_3858;
}

moonbit_string_t _M0MP26biolab8bio__seq8TreeNode10to__newick(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1238
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1237;
  moonbit_string_t _result_3859;
  #line 342 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 343 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L3bufS1237 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 344 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0MP26biolab8bio__seq8TreeNode18to__newick__helper(_M0L4selfS1238, _M0L3bufS1237);
  #line 345 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1237, (moonbit_string_t)moonbit_string_literal_29.data);
  #line 346 "/home/developer/Documents/1/skbio_tree.mbt"
  _result_3859 = _M0MPB13StringBuilder10to__string(_M0L3bufS1237);
  moonbit_decref(_M0L3bufS1237);
  return _result_3859;
}

int32_t _M0MP26biolab8bio__seq8TreeNode18to__newick__helper(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1219,
  struct _M0TPB13StringBuilder* _M0L3bufS1220
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8childrenS3172;
  int32_t _M0L6_2atmpS3171;
  moonbit_string_t _M0L1nS1228;
  moonbit_string_t _M0L7_2abindS1229;
  double _M0L1lS1233;
  void* _M0L7_2abindS1234;
  moonbit_string_t _M0L6_2atmpS3176;
  #line 350 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L8childrenS3172 = _M0L4selfS1219->$3;
  moonbit_incref(_M0L8childrenS3172);
  #line 351 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3171
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L8childrenS3172);
  moonbit_decref(_M0L8childrenS3172);
  if (_M0L6_2atmpS3171 > 0) {
    struct _M0TPB8MutLocalGbE* _M0L5firstS1221;
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1222;
    int32_t _M0L7_2abindS1223;
    int32_t _M0L2__S1224;
    #line 352 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1220, 40);
    _M0L5firstS1221
    = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
    Moonbit_object_header(_M0L5firstS1221)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L5firstS1221->$0 = 1;
    _M0L7_2abindS1222 = _M0L4selfS1219->$3;
    _M0L7_2abindS1223 = _M0L7_2abindS1222->$1;
    moonbit_incref(_M0L7_2abindS1222);
    _M0L2__S1224 = 0;
    while (1) {
      if (_M0L2__S1224 < _M0L7_2abindS1223) {
        struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS3175 =
          _M0L7_2abindS1222->$0;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1225 =
          (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS3175[
            _M0L2__S1224
          ];
        int32_t _M0L3valS3173 = _M0L5firstS1221->$0;
        int32_t _M0L6_2atmpS3174;
        if (!_M0L3valS3173) {
          moonbit_incref(_M0L5childS1225);
          #line 356 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1220, 44);
        } else {
          moonbit_incref(_M0L5childS1225);
        }
        _M0L5firstS1221->$0 = 0;
        #line 359 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0MP26biolab8bio__seq8TreeNode18to__newick__helper(_M0L5childS1225, _M0L3bufS1220);
        moonbit_decref(_M0L5childS1225);
        _M0L6_2atmpS3174 = _M0L2__S1224 + 1;
        _M0L2__S1224 = _M0L6_2atmpS3174;
        continue;
      } else {
        moonbit_decref(_M0L7_2abindS1222);
        moonbit_decref(_M0L5firstS1221);
      }
      break;
    }
    #line 361 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1220, 41);
  }
  _M0L7_2abindS1229 = _M0L4selfS1219->$0;
  if (_M0L7_2abindS1229 == 0) {
    
  } else {
    moonbit_string_t _M0L7_2aSomeS1230 = _M0L7_2abindS1229;
    moonbit_string_t _M0L4_2anS1231 = _M0L7_2aSomeS1230;
    moonbit_incref(_M0L4_2anS1231);
    _M0L1nS1228 = _M0L4_2anS1231;
    goto join_1227;
  }
  goto joinlet_3861;
  join_1227:;
  #line 365 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1220, _M0L1nS1228);
  moonbit_decref(_M0L1nS1228);
  joinlet_3861:;
  _M0L7_2abindS1234 = _M0L4selfS1219->$1;
  switch (Moonbit_object_tag(_M0L7_2abindS1234)) {
    case 1: {
      struct _M0DTPC16option6OptionGdE4Some* _M0L7_2aSomeS1235 =
        (struct _M0DTPC16option6OptionGdE4Some*)_M0L7_2abindS1234;
      double _M0L4_2alS1236 = _M0L7_2aSomeS1235->$0;
      _M0L1lS1233 = _M0L4_2alS1236;
      goto join_1232;
      break;
    }
    default:
      break;
  }
  goto joinlet_3862;
  join_1232:;
  #line 371 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1220, 58);
  #line 372 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3176 = _M0MPC16double6Double10to__string(_M0L1lS1233);
  #line 372 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1220, _M0L6_2atmpS3176);
  moonbit_decref(_M0L6_2atmpS3176);
  joinlet_3862:;
  return 0;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__2(
  
) {
  moonbit_string_t _M0L6newickS1217;
  struct moonbit_result_1 _tmp_3863;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1218;
  int32_t _M0L6_2atmpS3167;
  moonbit_string_t _M0L6_2atmpS3168;
  struct moonbit_result_0 _result_3865;
  #line 87 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 88 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1217 = _M0FP26biolab8bio__seq21gen__balanced__newick(12);
  #line 89 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3863 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1217);
  moonbit_decref(_M0L6newickS1217);
  if (_tmp_3863.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3169 =
      _tmp_3863.data.ok;
    _M0L4treeS1218 = _M0L5_2aokS3169;
  } else {
    void* const _M0L6_2aerrS3170 = _tmp_3863.data.err;
    struct moonbit_result_0 _result_3864;
    _result_3864.tag = 0;
    _result_3864.data.err = _M0L6_2aerrS3170;
    return _result_3864;
  }
  #line 90 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3167
  = _M0MP26biolab8bio__seq8TreeNode11count__tips(_M0L4treeS1218);
  moonbit_decref(_M0L4treeS1218);
  _M0L6_2atmpS3168 = 0;
  #line 90 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3865
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3167, 4096, _M0L6_2atmpS3168, (moonbit_string_t)moonbit_string_literal_32.data);
  if (_M0L6_2atmpS3168) {
    moonbit_decref(_M0L6_2atmpS3168);
  }
  return _result_3865;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__1(
  
) {
  moonbit_string_t _M0L6newickS1215;
  struct moonbit_result_1 _tmp_3866;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1216;
  int32_t _M0L6_2atmpS3163;
  moonbit_string_t _M0L6_2atmpS3164;
  struct moonbit_result_0 _result_3868;
  #line 80 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 81 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6newickS1215 = _M0FP26biolab8bio__seq21gen__balanced__newick(8);
  #line 82 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3866 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6newickS1215);
  moonbit_decref(_M0L6newickS1215);
  if (_tmp_3866.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3165 =
      _tmp_3866.data.ok;
    _M0L4treeS1216 = _M0L5_2aokS3165;
  } else {
    void* const _M0L6_2aerrS3166 = _tmp_3866.data.err;
    struct moonbit_result_0 _result_3867;
    _result_3867.tag = 0;
    _result_3867.data.err = _M0L6_2aerrS3166;
    return _result_3867;
  }
  #line 83 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3163
  = _M0MP26biolab8bio__seq8TreeNode11count__tips(_M0L4treeS1216);
  moonbit_decref(_M0L4treeS1216);
  _M0L6_2atmpS3164 = 0;
  #line 83 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3868
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3163, 256, _M0L6_2atmpS3164, (moonbit_string_t)moonbit_string_literal_33.data);
  if (_M0L6_2atmpS3164) {
    moonbit_decref(_M0L6_2atmpS3164);
  }
  return _result_3868;
}

struct moonbit_result_0 _M0FP26biolab8bio__seq53____test__736b62696f5f747265655f62656e63682e6d6274__0(
  
) {
  moonbit_string_t _M0L6tree10S1211;
  struct moonbit_result_1 _tmp_3869;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4treeS1212;
  int32_t _M0L6_2atmpS3153;
  moonbit_string_t _M0L6_2atmpS3154;
  struct moonbit_result_0 _tmp_3871;
  moonbit_string_t _M0L3catS1213;
  struct moonbit_result_1 _tmp_3873;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L7cattreeS1214;
  int32_t _M0L6_2atmpS3157;
  moonbit_string_t _M0L6_2atmpS3158;
  struct moonbit_result_0 _result_3875;
  #line 69 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 70 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6tree10S1211 = _M0FP26biolab8bio__seq21gen__balanced__newick(10);
  #line 71 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3869 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L6tree10S1211);
  moonbit_decref(_M0L6tree10S1211);
  if (_tmp_3869.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3161 =
      _tmp_3869.data.ok;
    _M0L4treeS1212 = _M0L5_2aokS3161;
  } else {
    void* const _M0L6_2aerrS3162 = _tmp_3869.data.err;
    struct moonbit_result_0 _result_3870;
    _result_3870.tag = 0;
    _result_3870.data.err = _M0L6_2aerrS3162;
    return _result_3870;
  }
  #line 72 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3153
  = _M0MP26biolab8bio__seq8TreeNode11count__tips(_M0L4treeS1212);
  moonbit_decref(_M0L4treeS1212);
  _M0L6_2atmpS3154 = 0;
  #line 72 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3871
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3153, 1024, _M0L6_2atmpS3154, (moonbit_string_t)moonbit_string_literal_34.data);
  if (_M0L6_2atmpS3154) {
    moonbit_decref(_M0L6_2atmpS3154);
  }
  if (_tmp_3871.tag) {
    int32_t const _M0L5_2aokS3155 = _tmp_3871.data.ok;
  } else {
    void* const _M0L6_2aerrS3156 = _tmp_3871.data.err;
    struct moonbit_result_0 _result_3872;
    _result_3872.tag = 0;
    _result_3872.data.err = _M0L6_2aerrS3156;
    return _result_3872;
  }
  #line 74 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L3catS1213 = _M0FP26biolab8bio__seq24gen__caterpillar__newick(100);
  #line 75 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _tmp_3873 = _M0MP26biolab8bio__seq8TreeNode12from__newick(_M0L3catS1213);
  moonbit_decref(_M0L3catS1213);
  if (_tmp_3873.tag) {
    struct _M0TP26biolab8bio__seq8TreeNode* const _M0L5_2aokS3159 =
      _tmp_3873.data.ok;
    _M0L7cattreeS1214 = _M0L5_2aokS3159;
  } else {
    void* const _M0L6_2aerrS3160 = _tmp_3873.data.err;
    struct moonbit_result_0 _result_3874;
    _result_3874.tag = 0;
    _result_3874.data.err = _M0L6_2aerrS3160;
    return _result_3874;
  }
  #line 76 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L6_2atmpS3157
  = _M0MP26biolab8bio__seq8TreeNode11count__tips(_M0L7cattreeS1214);
  moonbit_decref(_M0L7cattreeS1214);
  _M0L6_2atmpS3158 = 0;
  #line 76 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3875
  = _M0FPB10assert__eqGiE(_M0L6_2atmpS3157, 100, _M0L6_2atmpS3158, (moonbit_string_t)moonbit_string_literal_35.data);
  if (_M0L6_2atmpS3158) {
    moonbit_decref(_M0L6_2atmpS3158);
  }
  return _result_3875;
}

struct moonbit_result_1 _M0MP26biolab8bio__seq8TreeNode12from__newick(
  moonbit_string_t _M0L6newickS1210
) {
  #line 336 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 337 "/home/developer/Documents/1/skbio_tree.mbt"
  return _M0FP26biolab8bio__seq31parse__tree__node__from__newick(_M0L6newickS1210);
}

struct moonbit_result_1 _M0FP26biolab8bio__seq31parse__tree__node__from__newick(
  moonbit_string_t _M0L5inputS1200
) {
  void* _M0L4NoneS3152;
  struct _M0TPC16string10StringView _M0L6_2atmpS3151;
  moonbit_string_t _M0L7trimmedS1199;
  moonbit_string_t _M0L7_2abindS1201;
  int32_t _M0L6_2atmpS3141;
  struct _M0TPC16string10StringView _M0L6_2atmpS3140;
  int32_t _M0L6_2atmpS3139;
  int32_t _M0L6_2atmpS3150;
  int32_t _M0L6_2atmpS3149;
  int64_t _M0L6_2atmpS3148;
  struct _M0TPC16string10StringView _M0L6_2atmpS3147;
  moonbit_string_t _M0L1sS1202;
  int32_t _M0L1nS1203;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4rootS1205;
  int32_t _M0L3posS1206;
  struct moonbit_result_2 _tmp_3878;
  struct _M0TURP26biolab8bio__seq8TreeNodeiE* _M0L7_2abindS1207;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L7_2arootS1208;
  int32_t _M0L6_2aposS1209;
  int32_t _M0L6_2acntS3734;
  struct moonbit_result_1 _result_3881;
  #line 389 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L4NoneS3152
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  #line 390 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3151
  = _M0MPC16string6String4trim(_M0L5inputS1200, _M0L4NoneS3152);
  moonbit_decref(_M0L4NoneS3152);
  #line 390 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L7trimmedS1199 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3151);
  moonbit_decref(_M0L6_2atmpS3151.$0);
  _M0L7_2abindS1201 = (moonbit_string_t)moonbit_string_literal_29.data;
  _M0L6_2atmpS3141 = Moonbit_array_length(_M0L7_2abindS1201);
  _M0L6_2atmpS3140
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1201, .$1 = 0, .$2 = _M0L6_2atmpS3141
  };
  #line 392 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3139
  = _M0MPC16string6String11has__suffix(_M0L7trimmedS1199, _M0L6_2atmpS3140);
  moonbit_decref(_M0L6_2atmpS3140.$0);
  if (!_M0L6_2atmpS3139) {
    void* _M0L49biolab_2fbio__seq_2eParseError_2eMissingSemicolonS3142;
    struct moonbit_result_1 _result_3876;
    moonbit_decref(_M0L7trimmedS1199);
    _M0L49biolab_2fbio__seq_2eParseError_2eMissingSemicolonS3142
    = (struct moonbit_object*)&moonbit_constant_constructor_3 + 1;
    _result_3876.tag = 0;
    _result_3876.data.err
    = _M0L49biolab_2fbio__seq_2eParseError_2eMissingSemicolonS3142;
    return _result_3876;
  }
  _M0L6_2atmpS3150 = Moonbit_array_length(_M0L7trimmedS1199);
  _M0L6_2atmpS3149 = _M0L6_2atmpS3150 - 1;
  _M0L6_2atmpS3148 = (int64_t)_M0L6_2atmpS3149;
  #line 396 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3147
  = _M0MPC16string6String11sub_2einner(_M0L7trimmedS1199, 0, _M0L6_2atmpS3148);
  moonbit_decref(_M0L7trimmedS1199);
  #line 396 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L1sS1202 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS3147);
  moonbit_decref(_M0L6_2atmpS3147.$0);
  _M0L1nS1203 = Moonbit_array_length(_M0L1sS1202);
  #line 398 "/home/developer/Documents/1/skbio_tree.mbt"
  _tmp_3878
  = _M0FP26biolab8bio__seq16parse__node__rec(_M0L1sS1202, 0, _M0L1nS1203);
  moonbit_decref(_M0L1sS1202);
  if (_tmp_3878.tag) {
    struct _M0TURP26biolab8bio__seq8TreeNodeiE* const _M0L5_2aokS3145 =
      _tmp_3878.data.ok;
    _M0L7_2abindS1207 = _M0L5_2aokS3145;
  } else {
    void* const _M0L6_2aerrS3146 = _tmp_3878.data.err;
    struct moonbit_result_1 _result_3879;
    _result_3879.tag = 0;
    _result_3879.data.err = _M0L6_2aerrS3146;
    return _result_3879;
  }
  _M0L7_2arootS1208 = _M0L7_2abindS1207->$0;
  _M0L6_2aposS1209 = _M0L7_2abindS1207->$1;
  _M0L6_2acntS3734
  = Moonbit_rc_count(Moonbit_object_header(_M0L7_2abindS1207));
  if (_M0L6_2acntS3734 > 1) {
    int32_t _M0L11_2anew__cntS3735 = _M0L6_2acntS3734 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2abindS1207), _M0L11_2anew__cntS3735);
    moonbit_incref(_M0L7_2arootS1208);
  } else if (_M0L6_2acntS3734 == 1) {
    #line 398 "/home/developer/Documents/1/skbio_tree.mbt"
    moonbit_free(_M0L7_2abindS1207);
  }
  _M0L4rootS1205 = _M0L7_2arootS1208;
  _M0L3posS1206 = _M0L6_2aposS1209;
  goto join_1204;
  join_1204:;
  if (_M0L3posS1206 < _M0L1nS1203) {
    moonbit_string_t _M0L6_2atmpS3144;
    void* _M0L46biolab_2fbio__seq_2eParseError_2eInvalidFormatS3143;
    struct moonbit_result_1 _result_3880;
    moonbit_decref(_M0L4rootS1205);
    #line 401 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3144
    = _M0IPC16string6StringPB4Show10to__string((moonbit_string_t)moonbit_string_literal_36.data);
    _M0L46biolab_2fbio__seq_2eParseError_2eInvalidFormatS3143
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidFormat));
    Moonbit_object_header(_M0L46biolab_2fbio__seq_2eParseError_2eInvalidFormatS3143)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 4);
    ((struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidFormat*)_M0L46biolab_2fbio__seq_2eParseError_2eInvalidFormatS3143)->$0
    = _M0L6_2atmpS3144;
    _result_3880.tag = 0;
    _result_3880.data.err
    = _M0L46biolab_2fbio__seq_2eParseError_2eInvalidFormatS3143;
    return _result_3880;
  }
  _result_3881.tag = 1;
  _result_3881.data.ok = _M0L4rootS1205;
  return _result_3881;
}

struct moonbit_result_2 _M0FP26biolab8bio__seq16parse__node__rec(
  moonbit_string_t _M0L1sS1166,
  int32_t _M0L3posS1167,
  int32_t _M0L3lenS1168
) {
  int32_t _M0L6_2atmpS3138;
  struct _M0TPB8MutLocalGiE* _M0L3posS1165;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8childrenS1169;
  moonbit_string_t _M0L6_2atmpS3137;
  struct _M0TPB8MutLocalGOsE* _M0L4nameS1170;
  void* _M0L4NoneS3136;
  struct _M0TPB8MutLocalGOdE* _M0L6lengthS1171;
  int32_t _M0L3valS3073;
  int32_t _if__result_3882;
  int32_t _M0L3valS3107;
  int32_t _M0L3valS3117;
  int32_t _if__result_3896;
  moonbit_string_t _M0L8_2afieldS3543;
  int32_t _M0L6_2acntS3742;
  moonbit_string_t _M0L3valS3135;
  void* _M0L4SomeS3130;
  void* _M0L8_2afieldS3542;
  int32_t _M0L6_2acntS3744;
  void* _M0L3valS3134;
  void* _M0L6_2atmpS3131;
  void* _M0L6_2atmpS3132;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3133;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4nodeS1198;
  int32_t _M0L3valS3129;
  struct _M0TURP26biolab8bio__seq8TreeNodeiE* _M0L8_2atupleS3128;
  struct moonbit_result_2 _result_3900;
  #line 408 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 413 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS3138
  = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3posS1167, _M0L3lenS1168);
  _M0L3posS1165
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3posS1165)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3posS1165->$0 = _M0L6_2atmpS3138;
  #line 415 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L8childrenS1169
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(0);
  _M0L6_2atmpS3137 = 0;
  _M0L4nameS1170
  = (struct _M0TPB8MutLocalGOsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGOsE));
  Moonbit_object_header(_M0L4nameS1170)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
  _M0L4nameS1170->$0 = _M0L6_2atmpS3137;
  _M0L4NoneS3136
  = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  _M0L6lengthS1171
  = (struct _M0TPB8MutLocalGOdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGOdE));
  Moonbit_object_header(_M0L6lengthS1171)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 33, 0);
  _M0L6lengthS1171->$0 = _M0L4NoneS3136;
  _M0L3valS3073 = _M0L3posS1165->$0;
  if (_M0L3valS3073 < _M0L3lenS1168) {
    int32_t _M0L3valS3072 = _M0L3posS1165->$0;
    int32_t _M0L6_2atmpS3071 = _M0L1sS1166[_M0L3valS3072];
    int32_t _M0L6_2atmpS3070;
    #line 419 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3070
    = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3071);
    _if__result_3882 = _M0L6_2atmpS3070 == 40;
  } else {
    _if__result_3882 = 0;
  }
  if (_if__result_3882) {
    int32_t _M0L3valS3075 = _M0L3posS1165->$0;
    int32_t _M0L6_2atmpS3074 = _M0L3valS3075 + 1;
    int32_t _M0L3valS3077;
    int32_t _M0L6_2atmpS3076;
    int32_t _M0L3valS3081;
    int32_t _if__result_3883;
    int32_t _M0L3valS3101;
    int32_t _if__result_3892;
    int32_t _M0L3valS3104;
    int32_t _M0L6_2atmpS3103;
    int32_t _M0L3valS3106;
    int32_t _M0L6_2atmpS3105;
    _M0L3posS1165->$0 = _M0L6_2atmpS3074;
    _M0L3valS3077 = _M0L3posS1165->$0;
    #line 421 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3076
    = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3077, _M0L3lenS1168);
    _M0L3posS1165->$0 = _M0L6_2atmpS3076;
    _M0L3valS3081 = _M0L3posS1165->$0;
    if (_M0L3valS3081 < _M0L3lenS1168) {
      int32_t _M0L3valS3080 = _M0L3posS1165->$0;
      int32_t _M0L6_2atmpS3079 = _M0L1sS1166[_M0L3valS3080];
      int32_t _M0L6_2atmpS3078;
      #line 423 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3078
      = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3079);
      _if__result_3883 = _M0L6_2atmpS3078 != 41;
    } else {
      _if__result_3883 = 0;
    }
    if (_if__result_3883) {
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1173;
      int32_t _M0L8new__posS1174;
      int32_t _M0L3valS3095 = _M0L3posS1165->$0;
      struct moonbit_result_2 _tmp_3885;
      struct _M0TURP26biolab8bio__seq8TreeNodeiE* _M0L7_2abindS1182;
      struct _M0TP26biolab8bio__seq8TreeNode* _M0L8_2achildS1183;
      int32_t _M0L11_2anew__posS1184;
      int32_t _M0L6_2acntS3738;
      #line 424 "/home/developer/Documents/1/skbio_tree.mbt"
      _tmp_3885
      = _M0FP26biolab8bio__seq16parse__node__rec(_M0L1sS1166, _M0L3valS3095, _M0L3lenS1168);
      if (_tmp_3885.tag) {
        struct _M0TURP26biolab8bio__seq8TreeNodeiE* const _M0L5_2aokS3096 =
          _tmp_3885.data.ok;
        _M0L7_2abindS1182 = _M0L5_2aokS3096;
      } else {
        void* const _M0L6_2aerrS3097 = _tmp_3885.data.err;
        struct moonbit_result_2 _result_3886;
        moonbit_decref(_M0L6lengthS1171);
        moonbit_decref(_M0L4nameS1170);
        moonbit_decref(_M0L8childrenS1169);
        moonbit_decref(_M0L3posS1165);
        _result_3886.tag = 0;
        _result_3886.data.err = _M0L6_2aerrS3097;
        return _result_3886;
      }
      _M0L8_2achildS1183 = _M0L7_2abindS1182->$0;
      _M0L11_2anew__posS1184 = _M0L7_2abindS1182->$1;
      _M0L6_2acntS3738
      = Moonbit_rc_count(Moonbit_object_header(_M0L7_2abindS1182));
      if (_M0L6_2acntS3738 > 1) {
        int32_t _M0L11_2anew__cntS3739 = _M0L6_2acntS3738 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2abindS1182), _M0L11_2anew__cntS3739);
        moonbit_incref(_M0L8_2achildS1183);
      } else if (_M0L6_2acntS3738 == 1) {
        #line 424 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L7_2abindS1182);
      }
      _M0L5childS1173 = _M0L8_2achildS1183;
      _M0L8new__posS1174 = _M0L11_2anew__posS1184;
      goto join_1172;
      goto joinlet_3884;
      join_1172:;
      _M0L3posS1165->$0 = _M0L8new__posS1174;
      #line 426 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L8childrenS1169, _M0L5childS1173);
      moonbit_decref(_M0L5childS1173);
      while (1) {
        int32_t _M0L3valS3085 = _M0L3posS1165->$0;
        int32_t _if__result_3888;
        if (_M0L3valS3085 < _M0L3lenS1168) {
          int32_t _M0L3valS3084 = _M0L3posS1165->$0;
          int32_t _M0L6_2atmpS3083 = _M0L1sS1166[_M0L3valS3084];
          int32_t _M0L6_2atmpS3082;
          #line 428 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0L6_2atmpS3082
          = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3083);
          _if__result_3888 = _M0L6_2atmpS3082 == 44;
        } else {
          _if__result_3888 = 0;
        }
        if (_if__result_3888) {
          int32_t _M0L3valS3087 = _M0L3posS1165->$0;
          int32_t _M0L6_2atmpS3086 = _M0L3valS3087 + 1;
          int32_t _M0L3valS3089;
          int32_t _M0L6_2atmpS3088;
          struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1176;
          int32_t _M0L8new__posS1177;
          int32_t _M0L3valS3092;
          struct moonbit_result_2 _tmp_3890;
          struct _M0TURP26biolab8bio__seq8TreeNodeiE* _M0L7_2abindS1178;
          struct _M0TP26biolab8bio__seq8TreeNode* _M0L8_2achildS1179;
          int32_t _M0L11_2anew__posS1180;
          int32_t _M0L6_2acntS3736;
          int32_t _M0L3valS3091;
          int32_t _M0L6_2atmpS3090;
          _M0L3posS1165->$0 = _M0L6_2atmpS3086;
          _M0L3valS3089 = _M0L3posS1165->$0;
          #line 430 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0L6_2atmpS3088
          = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3089, _M0L3lenS1168);
          _M0L3posS1165->$0 = _M0L6_2atmpS3088;
          _M0L3valS3092 = _M0L3posS1165->$0;
          #line 431 "/home/developer/Documents/1/skbio_tree.mbt"
          _tmp_3890
          = _M0FP26biolab8bio__seq16parse__node__rec(_M0L1sS1166, _M0L3valS3092, _M0L3lenS1168);
          if (_tmp_3890.tag) {
            struct _M0TURP26biolab8bio__seq8TreeNodeiE* const _M0L5_2aokS3093 =
              _tmp_3890.data.ok;
            _M0L7_2abindS1178 = _M0L5_2aokS3093;
          } else {
            void* const _M0L6_2aerrS3094 = _tmp_3890.data.err;
            struct moonbit_result_2 _result_3891;
            moonbit_decref(_M0L6lengthS1171);
            moonbit_decref(_M0L4nameS1170);
            moonbit_decref(_M0L8childrenS1169);
            moonbit_decref(_M0L3posS1165);
            _result_3891.tag = 0;
            _result_3891.data.err = _M0L6_2aerrS3094;
            return _result_3891;
          }
          _M0L8_2achildS1179 = _M0L7_2abindS1178->$0;
          _M0L11_2anew__posS1180 = _M0L7_2abindS1178->$1;
          _M0L6_2acntS3736
          = Moonbit_rc_count(Moonbit_object_header(_M0L7_2abindS1178));
          if (_M0L6_2acntS3736 > 1) {
            int32_t _M0L11_2anew__cntS3737 = _M0L6_2acntS3736 - 1;
            Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2abindS1178), _M0L11_2anew__cntS3737);
            moonbit_incref(_M0L8_2achildS1179);
          } else if (_M0L6_2acntS3736 == 1) {
            #line 431 "/home/developer/Documents/1/skbio_tree.mbt"
            moonbit_free(_M0L7_2abindS1178);
          }
          _M0L5childS1176 = _M0L8_2achildS1179;
          _M0L8new__posS1177 = _M0L11_2anew__posS1180;
          goto join_1175;
          goto joinlet_3889;
          join_1175:;
          _M0L3posS1165->$0 = _M0L8new__posS1177;
          #line 433 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(_M0L8childrenS1169, _M0L5childS1176);
          moonbit_decref(_M0L5childS1176);
          _M0L3valS3091 = _M0L3posS1165->$0;
          #line 434 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0L6_2atmpS3090
          = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3091, _M0L3lenS1168);
          _M0L3posS1165->$0 = _M0L6_2atmpS3090;
          joinlet_3889:;
          continue;
        }
        break;
      }
      joinlet_3884:;
    }
    _M0L3valS3101 = _M0L3posS1165->$0;
    if (_M0L3valS3101 >= _M0L3lenS1168) {
      _if__result_3892 = 1;
    } else {
      int32_t _M0L3valS3100 = _M0L3posS1165->$0;
      int32_t _M0L6_2atmpS3099 = _M0L1sS1166[_M0L3valS3100];
      int32_t _M0L6_2atmpS3098;
      #line 438 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3098
      = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3099);
      _if__result_3892 = _M0L6_2atmpS3098 != 41;
    }
    if (_if__result_3892) {
      void* _M0L54biolab_2fbio__seq_2eParseError_2eUnbalancedParenthesesS3102;
      struct moonbit_result_2 _result_3893;
      moonbit_decref(_M0L6lengthS1171);
      moonbit_decref(_M0L4nameS1170);
      moonbit_decref(_M0L8childrenS1169);
      moonbit_decref(_M0L3posS1165);
      _M0L54biolab_2fbio__seq_2eParseError_2eUnbalancedParenthesesS3102
      = (struct moonbit_object*)&moonbit_constant_constructor_2 + 1;
      _result_3893.tag = 0;
      _result_3893.data.err
      = _M0L54biolab_2fbio__seq_2eParseError_2eUnbalancedParenthesesS3102;
      return _result_3893;
    }
    _M0L3valS3104 = _M0L3posS1165->$0;
    _M0L6_2atmpS3103 = _M0L3valS3104 + 1;
    _M0L3posS1165->$0 = _M0L6_2atmpS3103;
    _M0L3valS3106 = _M0L3posS1165->$0;
    #line 442 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3105
    = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3106, _M0L3lenS1168);
    _M0L3posS1165->$0 = _M0L6_2atmpS3105;
  }
  _M0L3valS3107 = _M0L3posS1165->$0;
  if (_M0L3valS3107 < _M0L3lenS1168) {
    int32_t _M0L3valS3113 = _M0L3posS1165->$0;
    int32_t _M0L6_2atmpS3112 = _M0L1sS1166[_M0L3valS3113];
    int32_t _M0L1cS1185;
    int32_t _if__result_3894;
    #line 446 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L1cS1185 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3112);
    if (_M0L1cS1185 != 58) {
      if (_M0L1cS1185 != 44) {
        if (_M0L1cS1185 != 41) {
          _if__result_3894 = _M0L1cS1185 != 59;
        } else {
          _if__result_3894 = 0;
        }
      } else {
        _if__result_3894 = 0;
      }
    } else {
      _if__result_3894 = 0;
    }
    if (_if__result_3894) {
      moonbit_string_t _M0L5labelS1187;
      int32_t _M0L8new__posS1188;
      int32_t _M0L3valS3111 = _M0L3posS1165->$0;
      struct _M0TUsiE* _M0L7_2abindS1189;
      moonbit_string_t _M0L8_2alabelS1190;
      int32_t _M0L11_2anew__posS1191;
      int32_t _M0L6_2acntS3740;
      moonbit_string_t _M0L6_2atmpS3108;
      moonbit_string_t _M0L6_2aoldS3545;
      int32_t _M0L3valS3110;
      int32_t _M0L6_2atmpS3109;
      #line 448 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L7_2abindS1189
      = _M0FP26biolab8bio__seq17parse__label__rec(_M0L1sS1166, _M0L3valS3111, _M0L3lenS1168);
      _M0L8_2alabelS1190 = _M0L7_2abindS1189->$0;
      _M0L11_2anew__posS1191 = _M0L7_2abindS1189->$1;
      _M0L6_2acntS3740
      = Moonbit_rc_count(Moonbit_object_header(_M0L7_2abindS1189));
      if (_M0L6_2acntS3740 > 1) {
        int32_t _M0L11_2anew__cntS3741 = _M0L6_2acntS3740 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2abindS1189), _M0L11_2anew__cntS3741);
        moonbit_incref(_M0L8_2alabelS1190);
      } else if (_M0L6_2acntS3740 == 1) {
        #line 448 "/home/developer/Documents/1/skbio_tree.mbt"
        moonbit_free(_M0L7_2abindS1189);
      }
      _M0L5labelS1187 = _M0L8_2alabelS1190;
      _M0L8new__posS1188 = _M0L11_2anew__posS1191;
      goto join_1186;
      goto joinlet_3895;
      join_1186:;
      _M0L6_2atmpS3108 = _M0L5labelS1187;
      _M0L6_2aoldS3545 = _M0L4nameS1170->$0;
      if (_M0L6_2aoldS3545) {
        moonbit_decref(_M0L6_2aoldS3545);
      }
      _M0L4nameS1170->$0 = _M0L6_2atmpS3108;
      _M0L3posS1165->$0 = _M0L8new__posS1188;
      _M0L3valS3110 = _M0L3posS1165->$0;
      #line 451 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3109
      = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3110, _M0L3lenS1168);
      _M0L3posS1165->$0 = _M0L6_2atmpS3109;
      joinlet_3895:;
    }
  }
  _M0L3valS3117 = _M0L3posS1165->$0;
  if (_M0L3valS3117 < _M0L3lenS1168) {
    int32_t _M0L3valS3116 = _M0L3posS1165->$0;
    int32_t _M0L6_2atmpS3115 = _M0L1sS1166[_M0L3valS3116];
    int32_t _M0L6_2atmpS3114;
    #line 455 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3114
    = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3115);
    _if__result_3896 = _M0L6_2atmpS3114 == 58;
  } else {
    _if__result_3896 = 0;
  }
  if (_if__result_3896) {
    int32_t _M0L3valS3119 = _M0L3posS1165->$0;
    int32_t _M0L6_2atmpS3118 = _M0L3valS3119 + 1;
    int32_t _M0L3valS3121;
    int32_t _M0L6_2atmpS3120;
    double _M0L3numS1193;
    int32_t _M0L8new__posS1194;
    int32_t _M0L3valS3125;
    struct moonbit_result_3 _tmp_3898;
    struct _M0TUdiE* _M0L7_2abindS1195;
    double _M0L6_2anumS1196;
    int32_t _M0L11_2anew__posS1197;
    void* _M0L4SomeS3122;
    void* _M0L6_2aoldS3544;
    int32_t _M0L3valS3124;
    int32_t _M0L6_2atmpS3123;
    _M0L3posS1165->$0 = _M0L6_2atmpS3118;
    _M0L3valS3121 = _M0L3posS1165->$0;
    #line 457 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3120
    = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3121, _M0L3lenS1168);
    _M0L3posS1165->$0 = _M0L6_2atmpS3120;
    _M0L3valS3125 = _M0L3posS1165->$0;
    #line 458 "/home/developer/Documents/1/skbio_tree.mbt"
    _tmp_3898
    = _M0FP26biolab8bio__seq18parse__number__rec(_M0L1sS1166, _M0L3valS3125, _M0L3lenS1168);
    if (_tmp_3898.tag) {
      struct _M0TUdiE* const _M0L5_2aokS3126 = _tmp_3898.data.ok;
      _M0L7_2abindS1195 = _M0L5_2aokS3126;
    } else {
      void* const _M0L6_2aerrS3127 = _tmp_3898.data.err;
      struct moonbit_result_2 _result_3899;
      moonbit_decref(_M0L6lengthS1171);
      moonbit_decref(_M0L4nameS1170);
      moonbit_decref(_M0L8childrenS1169);
      moonbit_decref(_M0L3posS1165);
      _result_3899.tag = 0;
      _result_3899.data.err = _M0L6_2aerrS3127;
      return _result_3899;
    }
    _M0L6_2anumS1196 = _M0L7_2abindS1195->$0;
    _M0L11_2anew__posS1197 = _M0L7_2abindS1195->$1;
    moonbit_decref(_M0L7_2abindS1195);
    _M0L3numS1193 = _M0L6_2anumS1196;
    _M0L8new__posS1194 = _M0L11_2anew__posS1197;
    goto join_1192;
    goto joinlet_3897;
    join_1192:;
    _M0L4SomeS3122
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGdE4Some));
    Moonbit_object_header(_M0L4SomeS3122)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 1);
    ((struct _M0DTPC16option6OptionGdE4Some*)_M0L4SomeS3122)->$0
    = _M0L3numS1193;
    _M0L6_2aoldS3544 = _M0L6lengthS1171->$0;
    moonbit_decref(_M0L6_2aoldS3544);
    _M0L6lengthS1171->$0 = _M0L4SomeS3122;
    _M0L3posS1165->$0 = _M0L8new__posS1194;
    _M0L3valS3124 = _M0L3posS1165->$0;
    #line 461 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3123
    = _M0FP26biolab8bio__seq8skip__ws(_M0L1sS1166, _M0L3valS3124, _M0L3lenS1168);
    _M0L3posS1165->$0 = _M0L6_2atmpS3123;
    joinlet_3897:;
  }
  _M0L8_2afieldS3543 = _M0L4nameS1170->$0;
  _M0L6_2acntS3742 = Moonbit_rc_count(Moonbit_object_header(_M0L4nameS1170));
  if (_M0L6_2acntS3742 > 1) {
    int32_t _M0L11_2anew__cntS3743 = _M0L6_2acntS3742 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L4nameS1170), _M0L11_2anew__cntS3743);
    if (_M0L8_2afieldS3543) {
      moonbit_incref(_M0L8_2afieldS3543);
    }
  } else if (_M0L6_2acntS3742 == 1) {
    #line 464 "/home/developer/Documents/1/skbio_tree.mbt"
    moonbit_free(_M0L4nameS1170);
  }
  _M0L3valS3135 = _M0L8_2afieldS3543;
  _M0L4SomeS3130
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGOsE4Some));
  Moonbit_object_header(_M0L4SomeS3130)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 36, 1);
  ((struct _M0DTPC16option6OptionGOsE4Some*)_M0L4SomeS3130)->$0
  = _M0L3valS3135;
  _M0L8_2afieldS3542 = _M0L6lengthS1171->$0;
  _M0L6_2acntS3744
  = Moonbit_rc_count(Moonbit_object_header(_M0L6lengthS1171));
  if (_M0L6_2acntS3744 > 1) {
    int32_t _M0L11_2anew__cntS3745 = _M0L6_2acntS3744 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L6lengthS1171), _M0L11_2anew__cntS3745);
    moonbit_incref(_M0L8_2afieldS3542);
  } else if (_M0L6_2acntS3744 == 1) {
    #line 464 "/home/developer/Documents/1/skbio_tree.mbt"
    moonbit_free(_M0L6lengthS1171);
  }
  _M0L3valS3134 = _M0L8_2afieldS3542;
  _M0L6_2atmpS3131 = _M0L3valS3134;
  _M0L6_2atmpS3132 = 0;
  _M0L6_2atmpS3133 = _M0L8childrenS1169;
  #line 464 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L4nodeS1198
  = _M0MP26biolab8bio__seq8TreeNode3new(_M0L4SomeS3130, _M0L6_2atmpS3131, _M0L6_2atmpS3132, _M0L6_2atmpS3133);
  moonbit_decref(_M0L4SomeS3130);
  if (_M0L6_2atmpS3131) {
    moonbit_decref(_M0L6_2atmpS3131);
  }
  if (_M0L6_2atmpS3132) {
    moonbit_decref(_M0L6_2atmpS3132);
  }
  if (_M0L6_2atmpS3133) {
    moonbit_decref(_M0L6_2atmpS3133);
  }
  _M0L3valS3129 = _M0L3posS1165->$0;
  moonbit_decref(_M0L3posS1165);
  _M0L8_2atupleS3128
  = (struct _M0TURP26biolab8bio__seq8TreeNodeiE*)moonbit_malloc(sizeof(struct _M0TURP26biolab8bio__seq8TreeNodeiE));
  Moonbit_object_header(_M0L8_2atupleS3128)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 39, 0);
  _M0L8_2atupleS3128->$0 = _M0L4nodeS1198;
  _M0L8_2atupleS3128->$1 = _M0L3valS3129;
  _result_3900.tag = 1;
  _result_3900.data.ok = _M0L8_2atupleS3128;
  return _result_3900;
}

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode3new(
  void* _M0L10name_2eoptS1154,
  void* _M0L12length_2eoptS1157,
  void* _M0L13support_2eoptS1160,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L14children_2eoptS1163
) {
  moonbit_string_t _M0L4nameS1153;
  void* _M0L6lengthS1156;
  void* _M0L7supportS1159;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8childrenS1162;
  struct _M0TP26biolab8bio__seq8TreeNode* _result_3901;
  switch (Moonbit_object_tag(_M0L10name_2eoptS1154)) {
    case 1: {
      struct _M0DTPC16option6OptionGOsE4Some* _M0L7_2aSomeS1155 =
        (struct _M0DTPC16option6OptionGOsE4Some*)_M0L10name_2eoptS1154;
      moonbit_string_t _M0L8_2afieldS3549 = _M0L7_2aSomeS1155->$0;
      if (_M0L8_2afieldS3549) {
        moonbit_incref(_M0L8_2afieldS3549);
      }
      _M0L4nameS1153 = _M0L8_2afieldS3549;
      break;
    }
    default: {
      _M0L4nameS1153 = 0;
      break;
    }
  }
  if (_M0L12length_2eoptS1157 == 0) {
    _M0L6lengthS1156
    = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  } else {
    void* _M0L7_2aSomeS1158 = _M0L12length_2eoptS1157;
    if (_M0L7_2aSomeS1158) {
      moonbit_incref(_M0L7_2aSomeS1158);
    }
    _M0L6lengthS1156 = _M0L7_2aSomeS1158;
  }
  if (_M0L13support_2eoptS1160 == 0) {
    _M0L7supportS1159
    = (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  } else {
    void* _M0L7_2aSomeS1161 = _M0L13support_2eoptS1160;
    if (_M0L7_2aSomeS1161) {
      moonbit_incref(_M0L7_2aSomeS1161);
    }
    _M0L7supportS1159 = _M0L7_2aSomeS1161;
  }
  if (_M0L14children_2eoptS1163 == 0) {
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS3069 =
      (struct _M0TP26biolab8bio__seq8TreeNode**)moonbit_empty_ref_array;
    _M0L8childrenS1162
    = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE));
    Moonbit_object_header(_M0L8childrenS1162)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 42, 0);
    _M0L8childrenS1162->$0 = _M0L6_2atmpS3069;
    _M0L8childrenS1162->$1 = 0;
  } else {
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2aSomeS1164 =
      _M0L14children_2eoptS1163;
    if (_M0L7_2aSomeS1164) {
      moonbit_incref(_M0L7_2aSomeS1164);
    }
    _M0L8childrenS1162 = _M0L7_2aSomeS1164;
  }
  _result_3901
  = _M0MP26biolab8bio__seq8TreeNode11new_2einner(_M0L4nameS1153, _M0L6lengthS1156, _M0L7supportS1159, _M0L8childrenS1162);
  if (_M0L4nameS1153) {
    moonbit_decref(_M0L4nameS1153);
  }
  moonbit_decref(_M0L6lengthS1156);
  moonbit_decref(_M0L7supportS1159);
  moonbit_decref(_M0L8childrenS1162);
  return _result_3901;
}

struct _M0TP26biolab8bio__seq8TreeNode* _M0MP26biolab8bio__seq8TreeNode11new_2einner(
  moonbit_string_t _M0L4nameS1149,
  void* _M0L6lengthS1150,
  void* _M0L7supportS1151,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8childrenS1152
) {
  struct _M0TP26biolab8bio__seq8TreeNode* _block_3902;
  #line 16 "/home/developer/Documents/1/skbio_tree.mbt"
  if (_M0L4nameS1149) {
    moonbit_incref(_M0L4nameS1149);
  }
  moonbit_incref(_M0L6lengthS1150);
  moonbit_incref(_M0L7supportS1151);
  moonbit_incref(_M0L8childrenS1152);
  _block_3902
  = (struct _M0TP26biolab8bio__seq8TreeNode*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq8TreeNode));
  Moonbit_object_header(_block_3902)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 45, 0);
  _block_3902->$0 = _M0L4nameS1149;
  _block_3902->$1 = _M0L6lengthS1150;
  _block_3902->$2 = _M0L7supportS1151;
  _block_3902->$3 = _M0L8childrenS1152;
  return _block_3902;
}

int32_t _M0FP26biolab8bio__seq8skip__ws(
  moonbit_string_t _M0L1sS1147,
  int32_t _M0L3posS1145,
  int32_t _M0L3lenS1146
) {
  struct _M0TPB8MutLocalGiE* _M0L3posS1144;
  int32_t _result_3905;
  #line 611 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L3posS1144
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3posS1144)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3posS1144->$0 = _M0L3posS1145;
  while (1) {
    int32_t _M0L3valS3066 = _M0L3posS1144->$0;
    int32_t _if__result_3904;
    if (_M0L3valS3066 < _M0L3lenS1146) {
      int32_t _M0L3valS3065 = _M0L3posS1144->$0;
      int32_t _M0L6_2atmpS3064 = _M0L1sS1147[_M0L3valS3065];
      int32_t _M0L6_2atmpS3063;
      #line 613 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS3063
      = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS3064);
      #line 613 "/home/developer/Documents/1/skbio_tree.mbt"
      _if__result_3904
      = _M0FP26biolab8bio__seq14is__whitespace(_M0L6_2atmpS3063);
    } else {
      _if__result_3904 = 0;
    }
    if (_if__result_3904) {
      int32_t _M0L3valS3068 = _M0L3posS1144->$0;
      int32_t _M0L6_2atmpS3067 = _M0L3valS3068 + 1;
      _M0L3posS1144->$0 = _M0L6_2atmpS3067;
      continue;
    }
    break;
  }
  _result_3905 = _M0L3posS1144->$0;
  moonbit_decref(_M0L3posS1144);
  return _result_3905;
}

struct moonbit_result_3 _M0FP26biolab8bio__seq18parse__number__rec(
  moonbit_string_t _M0L1sS1128,
  int32_t _M0L3posS1118,
  int32_t _M0L3lenS1127
) {
  struct _M0TPB8MutLocalGiE* _M0L3posS1117;
  struct _M0TPB8MutLocalGdE* _M0L6resultS1119;
  struct _M0TPB8MutLocalGdE* _M0L7decimalS1120;
  struct _M0TPB8MutLocalGiE* _M0L15decimal__placesS1121;
  struct _M0TPB8MutLocalGbE* _M0L12is__negativeS1122;
  struct _M0TPB8MutLocalGbE* _M0L14is__scientificS1123;
  struct _M0TPB8MutLocalGiE* _M0L9sci__signS1124;
  struct _M0TPB8MutLocalGiE* _M0L10sci__valueS1125;
  int32_t _M0L4zeroS1126;
  int32_t _M0L3valS2983;
  int32_t _if__result_3906;
  struct _M0TPB8MutLocalGbE* _M0L10has__digitS1129;
  int32_t _M0L3valS3031;
  int32_t _M0L3valS3034;
  double _M0L3valS3061;
  double _M0L3valS3062;
  double _M0L6_2atmpS3060;
  struct _M0TPB8MutLocalGdE* _M0L13final__resultS1140;
  int32_t _result_3916;
  int32_t _result_3918;
  double _M0L3valS3058;
  int32_t _M0L3valS3059;
  struct _M0TUdiE* _M0L8_2atupleS3057;
  struct moonbit_result_3 _result_3919;
  #line 506 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L3posS1117
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3posS1117)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3posS1117->$0 = _M0L3posS1118;
  _M0L6resultS1119
  = (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
  Moonbit_object_header(_M0L6resultS1119)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L6resultS1119->$0 = 0x0p+0;
  _M0L7decimalS1120
  = (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
  Moonbit_object_header(_M0L7decimalS1120)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L7decimalS1120->$0 = 0x0p+0;
  _M0L15decimal__placesS1121
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L15decimal__placesS1121)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L15decimal__placesS1121->$0 = 0;
  _M0L12is__negativeS1122
  = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
  Moonbit_object_header(_M0L12is__negativeS1122)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L12is__negativeS1122->$0 = 0;
  _M0L14is__scientificS1123
  = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
  Moonbit_object_header(_M0L14is__scientificS1123)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L14is__scientificS1123->$0 = 0;
  _M0L9sci__signS1124
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L9sci__signS1124)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L9sci__signS1124->$0 = 1;
  _M0L10sci__valueS1125
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L10sci__valueS1125)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L10sci__valueS1125->$0 = 0;
  _M0L4zeroS1126 = 48;
  _M0L3valS2983 = _M0L3posS1117->$0;
  if (_M0L3valS2983 < _M0L3lenS1127) {
    int32_t _M0L3valS2982 = _M0L3posS1117->$0;
    int32_t _M0L6_2atmpS2981 = _M0L1sS1128[_M0L3valS2982];
    int32_t _M0L6_2atmpS2980 = (int32_t)_M0L6_2atmpS2981;
    _if__result_3906 = _M0L6_2atmpS2980 == 45;
  } else {
    _if__result_3906 = 0;
  }
  if (_if__result_3906) {
    int32_t _M0L3valS2985;
    int32_t _M0L6_2atmpS2984;
    _M0L12is__negativeS1122->$0 = 1;
    _M0L3valS2985 = _M0L3posS1117->$0;
    _M0L6_2atmpS2984 = _M0L3valS2985 + 1;
    _M0L3posS1117->$0 = _M0L6_2atmpS2984;
  }
  _M0L10has__digitS1129
  = (struct _M0TPB8MutLocalGbE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGbE));
  Moonbit_object_header(_M0L10has__digitS1129)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L10has__digitS1129->$0 = 0;
  while (1) {
    int32_t _M0L3valS2986 = _M0L3posS1117->$0;
    if (_M0L3valS2986 < _M0L3lenS1127) {
      int32_t _M0L3valS3030 = _M0L3posS1117->$0;
      int32_t _M0L6_2atmpS3029 = _M0L1sS1128[_M0L3valS3030];
      int32_t _M0L1cS1130 = (int32_t)_M0L6_2atmpS3029;
      int32_t _if__result_3910;
      int32_t _if__result_3913;
      double _M0L3valS3026;
      double _M0L6_2atmpS3023;
      int32_t _M0L6_2atmpS3025;
      double _M0L6_2atmpS3024;
      double _M0L6_2atmpS3022;
      int32_t _M0L3valS3028;
      int32_t _M0L6_2atmpS3027;
      if (_M0L1cS1130 == 46) {
        int32_t _M0L3valS2988 = _M0L3posS1117->$0;
        int32_t _M0L6_2atmpS2987 = _M0L3valS2988 + 1;
        _M0L3posS1117->$0 = _M0L6_2atmpS2987;
        while (1) {
          int32_t _M0L3valS2989 = _M0L3posS1117->$0;
          if (_M0L3valS2989 < _M0L3lenS1127) {
            int32_t _M0L3valS3001 = _M0L3posS1117->$0;
            int32_t _M0L6_2atmpS3000 = _M0L1sS1128[_M0L3valS3001];
            int32_t _M0L2dcS1131 = (int32_t)_M0L6_2atmpS3000;
            int32_t _if__result_3909;
            double _M0L3valS2995;
            double _M0L6_2atmpS2992;
            int32_t _M0L6_2atmpS2994;
            double _M0L6_2atmpS2993;
            double _M0L6_2atmpS2991;
            int32_t _M0L3valS2997;
            int32_t _M0L6_2atmpS2996;
            int32_t _M0L3valS2999;
            int32_t _M0L6_2atmpS2998;
            if (_M0L2dcS1131 < _M0L4zeroS1126) {
              _if__result_3909 = 1;
            } else {
              int32_t _M0L6_2atmpS2990 = _M0L4zeroS1126 + 9;
              _if__result_3909 = _M0L2dcS1131 > _M0L6_2atmpS2990;
            }
            if (_if__result_3909) {
              break;
            }
            _M0L3valS2995 = _M0L7decimalS1120->$0;
            _M0L6_2atmpS2992 = _M0L3valS2995 * 0x1.4p+3;
            _M0L6_2atmpS2994 = _M0L2dcS1131 - _M0L4zeroS1126;
            _M0L6_2atmpS2993 = (double)_M0L6_2atmpS2994;
            _M0L6_2atmpS2991 = _M0L6_2atmpS2992 + _M0L6_2atmpS2993;
            _M0L7decimalS1120->$0 = _M0L6_2atmpS2991;
            _M0L3valS2997 = _M0L15decimal__placesS1121->$0;
            _M0L6_2atmpS2996 = _M0L3valS2997 + 1;
            _M0L15decimal__placesS1121->$0 = _M0L6_2atmpS2996;
            _M0L10has__digitS1129->$0 = 1;
            _M0L3valS2999 = _M0L3posS1117->$0;
            _M0L6_2atmpS2998 = _M0L3valS2999 + 1;
            _M0L3posS1117->$0 = _M0L6_2atmpS2998;
            continue;
          }
          break;
        }
        break;
      }
      if (_M0L1cS1130 == 101) {
        _if__result_3910 = 1;
      } else {
        _if__result_3910 = _M0L1cS1130 == 69;
      }
      if (_if__result_3910) {
        int32_t _M0L3valS3003;
        int32_t _M0L6_2atmpS3002;
        int32_t _M0L3valS3004;
        _M0L14is__scientificS1123->$0 = 1;
        _M0L3valS3003 = _M0L3posS1117->$0;
        _M0L6_2atmpS3002 = _M0L3valS3003 + 1;
        _M0L3posS1117->$0 = _M0L6_2atmpS3002;
        _M0L3valS3004 = _M0L3posS1117->$0;
        if (_M0L3valS3004 < _M0L3lenS1127) {
          int32_t _M0L3valS3010 = _M0L3posS1117->$0;
          int32_t _M0L6_2atmpS3009 = _M0L1sS1128[_M0L3valS3010];
          int32_t _M0L2scS1134 = (int32_t)_M0L6_2atmpS3009;
          if (_M0L2scS1134 == 45) {
            int32_t _M0L3valS3006;
            int32_t _M0L6_2atmpS3005;
            _M0L9sci__signS1124->$0 = -1;
            _M0L3valS3006 = _M0L3posS1117->$0;
            _M0L6_2atmpS3005 = _M0L3valS3006 + 1;
            _M0L3posS1117->$0 = _M0L6_2atmpS3005;
          } else if (_M0L2scS1134 == 43) {
            int32_t _M0L3valS3008 = _M0L3posS1117->$0;
            int32_t _M0L6_2atmpS3007 = _M0L3valS3008 + 1;
            _M0L3posS1117->$0 = _M0L6_2atmpS3007;
          }
        }
        while (1) {
          int32_t _M0L3valS3011 = _M0L3posS1117->$0;
          if (_M0L3valS3011 < _M0L3lenS1127) {
            int32_t _M0L3valS3020 = _M0L3posS1117->$0;
            int32_t _M0L6_2atmpS3019 = _M0L1sS1128[_M0L3valS3020];
            int32_t _M0L2dcS1135 = (int32_t)_M0L6_2atmpS3019;
            int32_t _if__result_3912;
            int32_t _M0L3valS3016;
            int32_t _M0L6_2atmpS3014;
            int32_t _M0L6_2atmpS3015;
            int32_t _M0L6_2atmpS3013;
            int32_t _M0L3valS3018;
            int32_t _M0L6_2atmpS3017;
            if (_M0L2dcS1135 < _M0L4zeroS1126) {
              _if__result_3912 = 1;
            } else {
              int32_t _M0L6_2atmpS3012 = _M0L4zeroS1126 + 9;
              _if__result_3912 = _M0L2dcS1135 > _M0L6_2atmpS3012;
            }
            if (_if__result_3912) {
              break;
            }
            _M0L3valS3016 = _M0L10sci__valueS1125->$0;
            _M0L6_2atmpS3014 = _M0L3valS3016 * 10;
            _M0L6_2atmpS3015 = _M0L2dcS1135 - _M0L4zeroS1126;
            _M0L6_2atmpS3013 = _M0L6_2atmpS3014 + _M0L6_2atmpS3015;
            _M0L10sci__valueS1125->$0 = _M0L6_2atmpS3013;
            _M0L3valS3018 = _M0L3posS1117->$0;
            _M0L6_2atmpS3017 = _M0L3valS3018 + 1;
            _M0L3posS1117->$0 = _M0L6_2atmpS3017;
            continue;
          }
          break;
        }
        break;
      }
      if (_M0L1cS1130 < _M0L4zeroS1126) {
        _if__result_3913 = 1;
      } else {
        int32_t _M0L6_2atmpS3021 = _M0L4zeroS1126 + 9;
        _if__result_3913 = _M0L1cS1130 > _M0L6_2atmpS3021;
      }
      if (_if__result_3913) {
        break;
      }
      _M0L3valS3026 = _M0L6resultS1119->$0;
      _M0L6_2atmpS3023 = _M0L3valS3026 * 0x1.4p+3;
      _M0L6_2atmpS3025 = _M0L1cS1130 - _M0L4zeroS1126;
      _M0L6_2atmpS3024 = (double)_M0L6_2atmpS3025;
      _M0L6_2atmpS3022 = _M0L6_2atmpS3023 + _M0L6_2atmpS3024;
      _M0L6resultS1119->$0 = _M0L6_2atmpS3022;
      _M0L10has__digitS1129->$0 = 1;
      _M0L3valS3028 = _M0L3posS1117->$0;
      _M0L6_2atmpS3027 = _M0L3valS3028 + 1;
      _M0L3posS1117->$0 = _M0L6_2atmpS3027;
      continue;
    }
    break;
  }
  _M0L3valS3031 = _M0L10has__digitS1129->$0;
  moonbit_decref(_M0L10has__digitS1129);
  if (!_M0L3valS3031) {
    moonbit_string_t _M0L6_2atmpS3033;
    void* _M0L46biolab_2fbio__seq_2eParseError_2eInvalidNumberS3032;
    struct moonbit_result_3 _result_3914;
    moonbit_decref(_M0L10sci__valueS1125);
    moonbit_decref(_M0L9sci__signS1124);
    moonbit_decref(_M0L14is__scientificS1123);
    moonbit_decref(_M0L12is__negativeS1122);
    moonbit_decref(_M0L15decimal__placesS1121);
    moonbit_decref(_M0L7decimalS1120);
    moonbit_decref(_M0L6resultS1119);
    moonbit_decref(_M0L3posS1117);
    #line 576 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS3033
    = _M0IPC16string6StringPB4Show10to__string((moonbit_string_t)moonbit_string_literal_1.data);
    _M0L46biolab_2fbio__seq_2eParseError_2eInvalidNumberS3032
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidNumber));
    Moonbit_object_header(_M0L46biolab_2fbio__seq_2eParseError_2eInvalidNumberS3032)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 51, 1);
    ((struct _M0DTPC15error5Error46biolab_2fbio__seq_2eParseError_2eInvalidNumber*)_M0L46biolab_2fbio__seq_2eParseError_2eInvalidNumberS3032)->$0
    = _M0L6_2atmpS3033;
    _result_3914.tag = 0;
    _result_3914.data.err
    = _M0L46biolab_2fbio__seq_2eParseError_2eInvalidNumberS3032;
    return _result_3914;
  }
  _M0L3valS3034 = _M0L15decimal__placesS1121->$0;
  if (_M0L3valS3034 > 0) {
    struct _M0TPB8MutLocalGdE* _M0L7divisorS1137 =
      (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
    int32_t _M0L1iS1138;
    double _M0L3valS3040;
    double _M0L3valS3041;
    double _M0L6_2atmpS3039;
    Moonbit_object_header(_M0L7divisorS1137)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L7divisorS1137->$0 = 0x1p+0;
    _M0L1iS1138 = 0;
    while (1) {
      int32_t _M0L3valS3035 = _M0L15decimal__placesS1121->$0;
      if (_M0L1iS1138 < _M0L3valS3035) {
        double _M0L3valS3037 = _M0L7divisorS1137->$0;
        double _M0L6_2atmpS3036 = _M0L3valS3037 * 0x1.4p+3;
        int32_t _M0L6_2atmpS3038;
        _M0L7divisorS1137->$0 = _M0L6_2atmpS3036;
        _M0L6_2atmpS3038 = _M0L1iS1138 + 1;
        _M0L1iS1138 = _M0L6_2atmpS3038;
        continue;
      } else {
        moonbit_decref(_M0L15decimal__placesS1121);
      }
      break;
    }
    _M0L3valS3040 = _M0L7decimalS1120->$0;
    _M0L3valS3041 = _M0L7divisorS1137->$0;
    moonbit_decref(_M0L7divisorS1137);
    _M0L6_2atmpS3039 = _M0L3valS3040 / _M0L3valS3041;
    _M0L7decimalS1120->$0 = _M0L6_2atmpS3039;
  } else {
    moonbit_decref(_M0L15decimal__placesS1121);
  }
  _M0L3valS3061 = _M0L6resultS1119->$0;
  moonbit_decref(_M0L6resultS1119);
  _M0L3valS3062 = _M0L7decimalS1120->$0;
  moonbit_decref(_M0L7decimalS1120);
  _M0L6_2atmpS3060 = _M0L3valS3061 + _M0L3valS3062;
  _M0L13final__resultS1140
  = (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
  Moonbit_object_header(_M0L13final__resultS1140)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L13final__resultS1140->$0 = _M0L6_2atmpS3060;
  _result_3916 = _M0L14is__scientificS1123->$0;
  moonbit_decref(_M0L14is__scientificS1123);
  if (_result_3916) {
    struct _M0TPB8MutLocalGdE* _M0L5pow10S1141 =
      (struct _M0TPB8MutLocalGdE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGdE));
    struct _M0TPB8MutLocalGiE* _M0L1jS1142;
    int32_t _M0L3valS3048;
    Moonbit_object_header(_M0L5pow10S1141)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L5pow10S1141->$0 = 0x1p+0;
    _M0L1jS1142
    = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
    Moonbit_object_header(_M0L1jS1142)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L1jS1142->$0 = 0;
    while (1) {
      int32_t _M0L3valS3042 = _M0L1jS1142->$0;
      int32_t _M0L3valS3043 = _M0L10sci__valueS1125->$0;
      if (_M0L3valS3042 < _M0L3valS3043) {
        double _M0L3valS3045 = _M0L5pow10S1141->$0;
        double _M0L6_2atmpS3044 = _M0L3valS3045 * 0x1.4p+3;
        int32_t _M0L3valS3047;
        int32_t _M0L6_2atmpS3046;
        _M0L5pow10S1141->$0 = _M0L6_2atmpS3044;
        _M0L3valS3047 = _M0L1jS1142->$0;
        _M0L6_2atmpS3046 = _M0L3valS3047 + 1;
        _M0L1jS1142->$0 = _M0L6_2atmpS3046;
        continue;
      } else {
        moonbit_decref(_M0L1jS1142);
        moonbit_decref(_M0L10sci__valueS1125);
      }
      break;
    }
    _M0L3valS3048 = _M0L9sci__signS1124->$0;
    moonbit_decref(_M0L9sci__signS1124);
    if (_M0L3valS3048 == 1) {
      double _M0L3valS3050 = _M0L13final__resultS1140->$0;
      double _M0L3valS3051 = _M0L5pow10S1141->$0;
      double _M0L6_2atmpS3049;
      moonbit_decref(_M0L5pow10S1141);
      _M0L6_2atmpS3049 = _M0L3valS3050 * _M0L3valS3051;
      _M0L13final__resultS1140->$0 = _M0L6_2atmpS3049;
    } else {
      double _M0L3valS3053 = _M0L13final__resultS1140->$0;
      double _M0L3valS3054 = _M0L5pow10S1141->$0;
      double _M0L6_2atmpS3052;
      moonbit_decref(_M0L5pow10S1141);
      _M0L6_2atmpS3052 = _M0L3valS3053 / _M0L3valS3054;
      _M0L13final__resultS1140->$0 = _M0L6_2atmpS3052;
    }
  } else {
    moonbit_decref(_M0L10sci__valueS1125);
    moonbit_decref(_M0L9sci__signS1124);
  }
  _result_3918 = _M0L12is__negativeS1122->$0;
  moonbit_decref(_M0L12is__negativeS1122);
  if (_result_3918) {
    double _M0L3valS3056 = _M0L13final__resultS1140->$0;
    double _M0L6_2atmpS3055 = -_M0L3valS3056;
    _M0L13final__resultS1140->$0 = _M0L6_2atmpS3055;
  }
  _M0L3valS3058 = _M0L13final__resultS1140->$0;
  moonbit_decref(_M0L13final__resultS1140);
  _M0L3valS3059 = _M0L3posS1117->$0;
  moonbit_decref(_M0L3posS1117);
  _M0L8_2atupleS3057
  = (struct _M0TUdiE*)moonbit_malloc(sizeof(struct _M0TUdiE));
  Moonbit_object_header(_M0L8_2atupleS3057)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L8_2atupleS3057->$0 = _M0L3valS3058;
  _M0L8_2atupleS3057->$1 = _M0L3valS3059;
  _result_3919.tag = 1;
  _result_3919.data.ok = _M0L8_2atupleS3057;
  return _result_3919;
}

struct _M0TUsiE* _M0FP26biolab8bio__seq17parse__label__rec(
  moonbit_string_t _M0L1sS1112,
  int32_t _M0L3posS1110,
  int32_t _M0L3lenS1111
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1108;
  struct _M0TPB8MutLocalGiE* _M0L3posS1109;
  int32_t _M0L3valS2948;
  int32_t _if__result_3920;
  moonbit_string_t _M0L6_2atmpS2978;
  int32_t _M0L3valS2979;
  struct _M0TUsiE* _block_3927;
  #line 470 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 471 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L3bufS1108 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3posS1109
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3posS1109)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3posS1109->$0 = _M0L3posS1110;
  _M0L3valS2948 = _M0L3posS1109->$0;
  if (_M0L3valS2948 < _M0L3lenS1111) {
    int32_t _M0L3valS2947 = _M0L3posS1109->$0;
    int32_t _M0L6_2atmpS2946 = _M0L1sS1112[_M0L3valS2947];
    int32_t _M0L6_2atmpS2945;
    #line 474 "/home/developer/Documents/1/skbio_tree.mbt"
    _M0L6_2atmpS2945
    = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2946);
    _if__result_3920 = _M0L6_2atmpS2945 == 39;
  } else {
    _if__result_3920 = 0;
  }
  if (_if__result_3920) {
    int32_t _M0L3valS2950 = _M0L3posS1109->$0;
    int32_t _M0L6_2atmpS2949 = _M0L3valS2950 + 1;
    int32_t _M0L3valS2970;
    int32_t _if__result_3924;
    _M0L3posS1109->$0 = _M0L6_2atmpS2949;
    while (1) {
      int32_t _M0L3valS2954 = _M0L3posS1109->$0;
      int32_t _if__result_3922;
      if (_M0L3valS2954 < _M0L3lenS1111) {
        int32_t _M0L3valS2953 = _M0L3posS1109->$0;
        int32_t _M0L6_2atmpS2952 = _M0L1sS1112[_M0L3valS2953];
        int32_t _M0L6_2atmpS2951;
        #line 476 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L6_2atmpS2951
        = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2952);
        _if__result_3922 = _M0L6_2atmpS2951 != 39;
      } else {
        _if__result_3922 = 0;
      }
      if (_if__result_3922) {
        int32_t _M0L3valS2966 = _M0L3posS1109->$0;
        int32_t _M0L6_2atmpS2965 = _M0L1sS1112[_M0L3valS2966];
        int32_t _M0L1cS1113;
        int32_t _if__result_3923;
        int32_t _M0L3valS2964;
        int32_t _M0L6_2atmpS2963;
        #line 477 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L1cS1113
        = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2965);
        if (_M0L1cS1113 == 39) {
          int32_t _M0L3valS2960 = _M0L3posS1109->$0;
          int32_t _M0L6_2atmpS2959 = _M0L3valS2960 + 1;
          if (_M0L6_2atmpS2959 < _M0L3lenS1111) {
            int32_t _M0L3valS2958 = _M0L3posS1109->$0;
            int32_t _M0L6_2atmpS2957 = _M0L3valS2958 + 1;
            int32_t _M0L6_2atmpS2956 = _M0L1sS1112[_M0L6_2atmpS2957];
            int32_t _M0L6_2atmpS2955;
            #line 480 "/home/developer/Documents/1/skbio_tree.mbt"
            _M0L6_2atmpS2955
            = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2956);
            _if__result_3923 = _M0L6_2atmpS2955 == 39;
          } else {
            _if__result_3923 = 0;
          }
        } else {
          _if__result_3923 = 0;
        }
        if (_if__result_3923) {
          int32_t _M0L3valS2962;
          int32_t _M0L6_2atmpS2961;
          #line 481 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1108, 39);
          _M0L3valS2962 = _M0L3posS1109->$0;
          _M0L6_2atmpS2961 = _M0L3valS2962 + 1;
          _M0L3posS1109->$0 = _M0L6_2atmpS2961;
        } else {
          #line 484 "/home/developer/Documents/1/skbio_tree.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1108, _M0L1cS1113);
        }
        _M0L3valS2964 = _M0L3posS1109->$0;
        _M0L6_2atmpS2963 = _M0L3valS2964 + 1;
        _M0L3posS1109->$0 = _M0L6_2atmpS2963;
        continue;
      }
      break;
    }
    _M0L3valS2970 = _M0L3posS1109->$0;
    if (_M0L3valS2970 < _M0L3lenS1111) {
      int32_t _M0L3valS2969 = _M0L3posS1109->$0;
      int32_t _M0L6_2atmpS2968 = _M0L1sS1112[_M0L3valS2969];
      int32_t _M0L6_2atmpS2967;
      #line 488 "/home/developer/Documents/1/skbio_tree.mbt"
      _M0L6_2atmpS2967
      = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2968);
      _if__result_3924 = _M0L6_2atmpS2967 == 39;
    } else {
      _if__result_3924 = 0;
    }
    if (_if__result_3924) {
      int32_t _M0L3valS2972 = _M0L3posS1109->$0;
      int32_t _M0L6_2atmpS2971 = _M0L3valS2972 + 1;
      _M0L3posS1109->$0 = _M0L6_2atmpS2971;
    }
  } else {
    while (1) {
      int32_t _M0L3valS2973 = _M0L3posS1109->$0;
      if (_M0L3valS2973 < _M0L3lenS1111) {
        int32_t _M0L3valS2977 = _M0L3posS1109->$0;
        int32_t _M0L6_2atmpS2976 = _M0L1sS1112[_M0L3valS2977];
        int32_t _M0L1cS1115;
        int32_t _if__result_3926;
        int32_t _M0L3valS2975;
        int32_t _M0L6_2atmpS2974;
        #line 493 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L1cS1115
        = _M0MPC16uint166UInt1616unsafe__to__char(_M0L6_2atmpS2976);
        if (_M0L1cS1115 == 58) {
          _if__result_3926 = 1;
        } else if (_M0L1cS1115 == 44) {
          _if__result_3926 = 1;
        } else if (_M0L1cS1115 == 41) {
          _if__result_3926 = 1;
        } else if (_M0L1cS1115 == 40) {
          _if__result_3926 = 1;
        } else {
          #line 494 "/home/developer/Documents/1/skbio_tree.mbt"
          _if__result_3926
          = _M0FP26biolab8bio__seq14is__whitespace(_M0L1cS1115);
        }
        if (_if__result_3926) {
          break;
        }
        #line 497 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1108, _M0L1cS1115);
        _M0L3valS2975 = _M0L3posS1109->$0;
        _M0L6_2atmpS2974 = _M0L3valS2975 + 1;
        _M0L3posS1109->$0 = _M0L6_2atmpS2974;
        continue;
      }
      break;
    }
  }
  #line 502 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS2978 = _M0MPB13StringBuilder10to__string(_M0L3bufS1108);
  moonbit_decref(_M0L3bufS1108);
  _M0L3valS2979 = _M0L3posS1109->$0;
  moonbit_decref(_M0L3posS1109);
  _block_3927 = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
  Moonbit_object_header(_block_3927)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
  _block_3927->$0 = _M0L6_2atmpS2978;
  _block_3927->$1 = _M0L3valS2979;
  return _block_3927;
}

int32_t _M0FP26biolab8bio__seq14is__whitespace(int32_t _M0L1cS1107) {
  #line 620 "/home/developer/Documents/1/skbio_tree.mbt"
  if (_M0L1cS1107 == 32) {
    return 1;
  } else if (_M0L1cS1107 == 9) {
    return 1;
  } else if (_M0L1cS1107 == 10) {
    return 1;
  } else {
    return _M0L1cS1107 == 13;
  }
}

int32_t _M0MP26biolab8bio__seq8TreeNode11count__tips(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1100
) {
  #line 312 "/home/developer/Documents/1/skbio_tree.mbt"
  #line 313 "/home/developer/Documents/1/skbio_tree.mbt"
  if (_M0MP26biolab8bio__seq8TreeNode7is__tip(_M0L4selfS1100)) {
    return 1;
  } else {
    struct _M0TPB8MutLocalGiE* _M0L5countS1101 =
      (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L7_2abindS1102;
    int32_t _M0L7_2abindS1103;
    int32_t _M0L2__S1104;
    int32_t _result_3929;
    Moonbit_object_header(_M0L5countS1101)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L5countS1101->$0 = 0;
    _M0L7_2abindS1102 = _M0L4selfS1100->$3;
    _M0L7_2abindS1103 = _M0L7_2abindS1102->$1;
    moonbit_incref(_M0L7_2abindS1102);
    _M0L2__S1104 = 0;
    while (1) {
      if (_M0L2__S1104 < _M0L7_2abindS1103) {
        struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS2944 =
          _M0L7_2abindS1102->$0;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L5childS1105 =
          (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS2944[
            _M0L2__S1104
          ];
        int32_t _M0L3valS2941 = _M0L5countS1101->$0;
        int32_t _M0L6_2atmpS2942;
        int32_t _M0L6_2atmpS2940;
        int32_t _M0L6_2atmpS2943;
        moonbit_incref(_M0L5childS1105);
        #line 318 "/home/developer/Documents/1/skbio_tree.mbt"
        _M0L6_2atmpS2942
        = _M0MP26biolab8bio__seq8TreeNode11count__tips(_M0L5childS1105);
        moonbit_decref(_M0L5childS1105);
        _M0L6_2atmpS2940 = _M0L3valS2941 + _M0L6_2atmpS2942;
        _M0L5countS1101->$0 = _M0L6_2atmpS2940;
        _M0L6_2atmpS2943 = _M0L2__S1104 + 1;
        _M0L2__S1104 = _M0L6_2atmpS2943;
        continue;
      } else {
        moonbit_decref(_M0L7_2abindS1102);
      }
      break;
    }
    _result_3929 = _M0L5countS1101->$0;
    moonbit_decref(_M0L5countS1101);
    return _result_3929;
  }
}

int32_t _M0MP26biolab8bio__seq8TreeNode7is__tip(
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L4selfS1099
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L8childrenS2939;
  int32_t _M0L6_2atmpS2938;
  #line 27 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L8childrenS2939 = _M0L4selfS1099->$3;
  moonbit_incref(_M0L8childrenS2939);
  #line 28 "/home/developer/Documents/1/skbio_tree.mbt"
  _M0L6_2atmpS2938
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L8childrenS2939);
  moonbit_decref(_M0L8childrenS2939);
  return _M0L6_2atmpS2938 == 0;
}

moonbit_string_t _M0FP26biolab8bio__seq24gen__caterpillar__newick(
  int32_t _M0L1nS1093
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1094;
  int32_t _M0L1iS1095;
  int32_t _M0L1iS1097;
  moonbit_string_t _result_3932;
  #line 42 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  if (_M0L1nS1093 < 2) {
    return (moonbit_string_t)moonbit_string_literal_37.data;
  }
  #line 46 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L3bufS1094 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L1iS1095 = 0;
  while (1) {
    int32_t _M0L6_2atmpS2932 = _M0L1nS1093 - 1;
    if (_M0L1iS1095 < _M0L6_2atmpS2932) {
      int32_t _M0L6_2atmpS2933;
      #line 49 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1094, 40);
      _M0L6_2atmpS2933 = _M0L1iS1095 + 1;
      _M0L1iS1095 = _M0L6_2atmpS2933;
      continue;
    }
    break;
  }
  #line 52 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1094, (moonbit_string_t)moonbit_string_literal_38.data);
  _M0L1iS1097 = 2;
  while (1) {
    if (_M0L1iS1097 < _M0L1nS1093) {
      int32_t _M0L6_2atmpS2935;
      moonbit_string_t _M0L6_2atmpS2934;
      moonbit_string_t _M0L6_2atmpS2936;
      int32_t _M0L6_2atmpS2937;
      #line 55 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1094, 41);
      #line 56 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1094, (moonbit_string_t)moonbit_string_literal_39.data);
      _M0L6_2atmpS2935 = _M0L1iS1097 - 2;
      #line 57 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0L6_2atmpS2934
      = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2935, 10);
      #line 57 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1094, _M0L6_2atmpS2934);
      moonbit_decref(_M0L6_2atmpS2934);
      #line 58 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1094, (moonbit_string_t)moonbit_string_literal_40.data);
      #line 59 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0L6_2atmpS2936 = _M0MPC13int3Int18to__string_2einner(_M0L1iS1097, 10);
      #line 59 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1094, _M0L6_2atmpS2936);
      moonbit_decref(_M0L6_2atmpS2936);
      #line 60 "/home/developer/Documents/1/skbio_tree_bench.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1094, (moonbit_string_t)moonbit_string_literal_41.data);
      _M0L6_2atmpS2937 = _M0L1iS1097 + 1;
      _M0L1iS1097 = _M0L6_2atmpS2937;
      continue;
    }
    break;
  }
  #line 63 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1094, 41);
  #line 64 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1094, 59);
  #line 65 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3932 = _M0MPB13StringBuilder10to__string(_M0L3bufS1094);
  moonbit_decref(_M0L3bufS1094);
  return _result_3932;
}

moonbit_string_t _M0FP26biolab8bio__seq21gen__balanced__newick(
  int32_t _M0L5depthS1092
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1090;
  int32_t* _M0L6_2atmpS2931;
  struct _M0TPB5ArrayGiE* _M0L7counterS1091;
  moonbit_string_t _result_3933;
  #line 7 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  #line 8 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0L3bufS1090 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS2931 = (int32_t*)moonbit_make_int32_array_raw(1);
  _M0L6_2atmpS2931[0] = 0;
  _M0L7counterS1091
  = (struct _M0TPB5ArrayGiE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGiE));
  Moonbit_object_header(_M0L7counterS1091)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 54, 0);
  _M0L7counterS1091->$0 = _M0L6_2atmpS2931;
  _M0L7counterS1091->$1 = 1;
  #line 10 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0FP26biolab8bio__seq18gen__balanced__rec(_M0L3bufS1090, _M0L5depthS1092, _M0L7counterS1091);
  moonbit_decref(_M0L7counterS1091);
  #line 11 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1090, 59);
  #line 12 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _result_3933 = _M0MPB13StringBuilder10to__string(_M0L3bufS1090);
  moonbit_decref(_M0L3bufS1090);
  return _result_3933;
}

int32_t _M0FP26biolab8bio__seq18gen__balanced__rec(
  struct _M0TPB13StringBuilder* _M0L3bufS1088,
  int32_t _M0L5depthS1087,
  struct _M0TPB5ArrayGiE* _M0L7counterS1089
) {
  int32_t _M0L6_2atmpS2929;
  int32_t _M0L6_2atmpS2930;
  #line 16 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  if (_M0L5depthS1087 == 0) {
    int32_t _M0L6_2atmpS2926;
    moonbit_string_t _M0L6_2atmpS2925;
    int32_t _M0L6_2atmpS2928;
    int32_t _M0L6_2atmpS2927;
    #line 22 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1088, (moonbit_string_t)moonbit_string_literal_42.data);
    #line 23 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0L6_2atmpS2926 = _M0MPC15array5Array2atGiE(_M0L7counterS1089, 0);
    #line 23 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0L6_2atmpS2925
    = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2926, 10);
    #line 23 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1088, _M0L6_2atmpS2925);
    moonbit_decref(_M0L6_2atmpS2925);
    #line 24 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0L6_2atmpS2928 = _M0MPC15array5Array2atGiE(_M0L7counterS1089, 0);
    _M0L6_2atmpS2927 = _M0L6_2atmpS2928 + 1;
    #line 24 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0MPC15array5Array3setGiE(_M0L7counterS1089, 0, _M0L6_2atmpS2927);
    #line 25 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1088, 58);
    #line 26 "/home/developer/Documents/1/skbio_tree_bench.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1088, (moonbit_string_t)moonbit_string_literal_43.data);
    return 0;
  }
  #line 29 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1088, 40);
  _M0L6_2atmpS2929 = _M0L5depthS1087 - 1;
  #line 30 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0FP26biolab8bio__seq18gen__balanced__rec(_M0L3bufS1088, _M0L6_2atmpS2929, _M0L7counterS1089);
  #line 31 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1088, 44);
  _M0L6_2atmpS2930 = _M0L5depthS1087 - 1;
  #line 32 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0FP26biolab8bio__seq18gen__balanced__rec(_M0L3bufS1088, _M0L6_2atmpS2930, _M0L7counterS1089);
  #line 33 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1088, 41);
  #line 35 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1088, 58);
  #line 36 "/home/developer/Documents/1/skbio_tree_bench.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1088, (moonbit_string_t)moonbit_string_literal_43.data);
  return 0;
}

int32_t _M0MPC15array5Array3setGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS1084,
  int32_t _M0L5indexS1085,
  int32_t _M0L5valueS1086
) {
  int32_t _M0L3lenS1083;
  int32_t _if__result_3934;
  #line 258 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1083 = _M0L4selfS1084->$1;
  if (_M0L5indexS1085 >= 0) {
    _if__result_3934 = _M0L5indexS1085 < _M0L3lenS1083;
  } else {
    _if__result_3934 = 0;
  }
  if (_if__result_3934) {
    int32_t* _M0L6_2atmpS2924;
    #line 263 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2924 = _M0MPC15array5Array6bufferGiE(_M0L4selfS1084);
    _M0L6_2atmpS2924[_M0L5indexS1085] = _M0L5valueS1086;
    moonbit_decref(_M0L6_2atmpS2924);
  } else {
    #line 262 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
  return 0;
}

struct _M0TP26biolab8bio__seq8TreeNode* _M0MPC15array5Array6removeGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS1080,
  int32_t _M0L5indexS1079
) {
  int32_t _if__result_3935;
  #line 384 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L5indexS1079 >= 0) {
    int32_t _M0L3lenS2912 = _M0L4selfS1080->$1;
    _if__result_3935 = _M0L5indexS1079 < _M0L3lenS2912;
  } else {
    _if__result_3935 = 0;
  }
  if (_if__result_3935) {
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS2921 =
      _M0L4selfS1080->$0;
    struct _M0TP26biolab8bio__seq8TreeNode* _M0L5valueS1081 =
      (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS2921[_M0L5indexS1079];
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS2913;
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS2914;
    int32_t _M0L6_2atmpS2915;
    int32_t _M0L3lenS2918;
    int32_t _M0L6_2atmpS2917;
    int32_t _M0L6_2atmpS2916;
    int32_t _M0L3lenS2920;
    int32_t _M0L6_2atmpS2919;
    moonbit_incref(_M0L5valueS1081);
    #line 392 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0L6_2atmpS2913
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS1080);
    #line 394 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0L6_2atmpS2914
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS1080);
    _M0L6_2atmpS2915 = _M0L5indexS1079 + 1;
    _M0L3lenS2918 = _M0L4selfS1080->$1;
    _M0L6_2atmpS2917 = _M0L3lenS2918 - _M0L5indexS1079;
    _M0L6_2atmpS2916 = _M0L6_2atmpS2917 - 1;
    #line 391 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS2913, _M0L5indexS1079, _M0L6_2atmpS2914, _M0L6_2atmpS2915, _M0L6_2atmpS2916);
    moonbit_decref(_M0L6_2atmpS2913);
    moonbit_decref(_M0L6_2atmpS2914);
    _M0L3lenS2920 = _M0L4selfS1080->$1;
    _M0L6_2atmpS2919 = _M0L3lenS2920 - 1;
    #line 398 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array28unsafe__truncate__to__lengthGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS1080, _M0L6_2atmpS2919);
    return _M0L5valueS1081;
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1082;
    int32_t _M0L3lenS2923;
    moonbit_string_t _M0L6_2atmpS2922;
    struct _M0TP26biolab8bio__seq8TreeNode* _result_3936;
    #line 387 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0L18_2astring__builderS1082
    = _M0MPB13StringBuilder21StringBuilder_2einner(60);
    #line 387 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1082, (moonbit_string_t)moonbit_string_literal_44.data);
    _M0L3lenS2923 = _M0L4selfS1080->$1;
    #line 387 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1082, _M0L3lenS2923);
    #line 387 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1082, (moonbit_string_t)moonbit_string_literal_45.data);
    #line 387 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1082, _M0L5indexS1079);
    #line 387 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0L6_2atmpS2922
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1082);
    moonbit_decref(_M0L18_2astring__builderS1082);
    #line 386 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _result_3936
    = _M0FPC15abort5abortGRP26biolab8bio__seq8TreeNodeE(_M0L6_2atmpS2922);
    moonbit_decref(_M0L6_2atmpS2922);
    return _result_3936;
  }
}

moonbit_string_t _M0MPC15array5Array2atGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS1068,
  int32_t _M0L5indexS1069
) {
  int32_t _M0L3lenS1067;
  int32_t _if__result_3937;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1067 = _M0L4selfS1068->$1;
  if (_M0L5indexS1069 >= 0) {
    _if__result_3937 = _M0L5indexS1069 < _M0L3lenS1067;
  } else {
    _if__result_3937 = 0;
  }
  if (_if__result_3937) {
    moonbit_string_t* _M0L6_2atmpS2908;
    moonbit_string_t _M0L6_2atmpS3556;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2908 = _M0MPC15array5Array6bufferGsE(_M0L4selfS1068);
    _M0L6_2atmpS3556 = (moonbit_string_t)_M0L6_2atmpS2908[_M0L5indexS1069];
    moonbit_incref(_M0L6_2atmpS3556);
    moonbit_decref(_M0L6_2atmpS2908);
    return _M0L6_2atmpS3556;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TP26biolab8bio__seq8TreeNode* _M0MPC15array5Array2atGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS1071,
  int32_t _M0L5indexS1072
) {
  int32_t _M0L3lenS1070;
  int32_t _if__result_3938;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1070 = _M0L4selfS1071->$1;
  if (_M0L5indexS1072 >= 0) {
    _if__result_3938 = _M0L5indexS1072 < _M0L3lenS1070;
  } else {
    _if__result_3938 = 0;
  }
  if (_if__result_3938) {
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS2909;
    struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS3557;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2909
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS1071);
    _M0L6_2atmpS3557
    = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L6_2atmpS2909[
        _M0L5indexS1072
      ];
    if (_M0L6_2atmpS3557) {
      moonbit_incref(_M0L6_2atmpS3557);
    }
    moonbit_decref(_M0L6_2atmpS2909);
    return _M0L6_2atmpS3557;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MPC15array5Array2atGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L4selfS1074,
  int32_t _M0L5indexS1075
) {
  int32_t _M0L3lenS1073;
  int32_t _if__result_3939;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1073 = _M0L4selfS1074->$1;
  if (_M0L5indexS1075 >= 0) {
    _if__result_3939 = _M0L5indexS1075 < _M0L3lenS1073;
  } else {
    _if__result_3939 = 0;
  }
  if (_if__result_3939) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L6_2atmpS2910;
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS3558;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2910
    = _M0MPC15array5Array6bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L4selfS1074);
    _M0L6_2atmpS3558
    = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)_M0L6_2atmpS2910[
        _M0L5indexS1075
      ];
    if (_M0L6_2atmpS3558) {
      moonbit_incref(_M0L6_2atmpS3558);
    }
    moonbit_decref(_M0L6_2atmpS2910);
    return _M0L6_2atmpS3558;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array5Array2atGiE(
  struct _M0TPB5ArrayGiE* _M0L4selfS1077,
  int32_t _M0L5indexS1078
) {
  int32_t _M0L3lenS1076;
  int32_t _if__result_3940;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS1076 = _M0L4selfS1077->$1;
  if (_M0L5indexS1078 >= 0) {
    _if__result_3940 = _M0L5indexS1078 < _M0L3lenS1076;
  } else {
    _if__result_3940 = 0;
  }
  if (_if__result_3940) {
    int32_t* _M0L6_2atmpS2911;
    int32_t _result_3941;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2911 = _M0MPC15array5Array6bufferGiE(_M0L4selfS1077);
    _result_3941 = (int32_t)_M0L6_2atmpS2911[_M0L5indexS1078];
    moonbit_decref(_M0L6_2atmpS2911);
    return _result_3941;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array5Array28unsafe__truncate__to__lengthGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS1063,
  int32_t _M0L8new__lenS1064
) {
  int32_t _M0L3lenS1062;
  #line 167 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1062 = _M0L4selfS1063->$1;
  if (_M0L8new__lenS1064 <= _M0L3lenS1062) {
    int32_t _M0L1iS1065 = _M0L8new__lenS1064;
    while (1) {
      if (_M0L1iS1065 < _M0L3lenS1062) {
        struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS2906 =
          _M0L4selfS1063->$0;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2aoldS3559;
        int32_t _M0L6_2atmpS2907;
        if (
          _M0L1iS1065 < 0
          || _M0L1iS1065 >= Moonbit_array_length(_M0L3bufS2906)
        ) {
          #line 171 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3559
        = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS2906[_M0L1iS1065];
        if (_M0L6_2aoldS3559) {
          moonbit_decref(_M0L6_2aoldS3559);
        }
        if (
          _M0L1iS1065 < 0
          || _M0L1iS1065 >= Moonbit_array_length(_M0L3bufS2906)
        ) {
          #line 171 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
          moonbit_panic();
        }
        _M0L3bufS2906[_M0L1iS1065] = 0;
        _M0L6_2atmpS2907 = _M0L1iS1065 + 1;
        _M0L1iS1065 = _M0L6_2atmpS2907;
        continue;
      }
      break;
    }
    _M0L4selfS1063->$1 = _M0L8new__lenS1064;
  } else {
    #line 169 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    moonbit_panic();
  }
  return 0;
}

int32_t _M0IPB9SourceLocPB4Show6output(
  moonbit_string_t _M0L4selfS1061,
  struct _M0TPB6Logger _M0L6loggerS1060
) {
  moonbit_string_t _M0L6_2atmpS2905;
  #line 43 "/home/developer/.moon/lib/core/builtin/autoloc.mbt"
  moonbit_incref(_M0L4selfS1061);
  _M0L6_2atmpS2905 = _M0L4selfS1061;
  #line 44 "/home/developer/.moon/lib/core/builtin/autoloc.mbt"
  _M0L6loggerS1060.$0->$method_0(_M0L6loggerS1060.$1, _M0L6_2atmpS2905);
  moonbit_decref(_M0L6_2atmpS2905);
  return 0;
}

int32_t _M0FPB7printlnGsE(moonbit_string_t _M0L5inputS1059) {
  moonbit_string_t _M0L6_2atmpS2904;
  #line 36 "/home/developer/.moon/lib/core/builtin/console.mbt"
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  _M0L6_2atmpS2904
  = _M0IPC16string6StringPB4Show10to__string(_M0L5inputS1059);
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  moonbit_println(_M0L6_2atmpS2904);
  moonbit_decref(_M0L6_2atmpS2904);
  return 0;
}

moonbit_string_t _M0MPC16double6Double10to__string(double _M0L4selfS1058) {
  #line 282 "/home/developer/.moon/lib/core/builtin/double.mbt"
  #line 284 "/home/developer/.moon/lib/core/builtin/double.mbt"
  return _M0FPB15ryu__to__string(_M0L4selfS1058);
}

moonbit_string_t _M0FPB15ryu__to__string(double _M0L3valS1045) {
  uint64_t _M0L4bitsS1046;
  uint64_t _M0L6_2atmpS2903;
  uint64_t _M0L6_2atmpS2902;
  int32_t _M0L8ieeeSignS1047;
  uint64_t _M0L12ieeeMantissaS1048;
  uint64_t _M0L6_2atmpS2901;
  uint64_t _M0L6_2atmpS2900;
  int32_t _M0L12ieeeExponentS1049;
  int32_t _if__result_3943;
  struct _M0TPB17FloatingDecimal64* _M0L7_2abindS1050;
  struct _M0TPB17FloatingDecimal64* _M0L1vS1051;
  moonbit_string_t _result_3945;
  #line 659 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  if (_M0L3valS1045 == 0x0p+0) {
    return (moonbit_string_t)moonbit_string_literal_46.data;
  }
  _M0L4bitsS1046 = *(int64_t*)&_M0L3valS1045;
  _M0L6_2atmpS2903 = _M0L4bitsS1046 >> 63;
  _M0L6_2atmpS2902 = _M0L6_2atmpS2903 & 1ull;
  _M0L8ieeeSignS1047 = _M0L6_2atmpS2902 != 0ull;
  _M0L12ieeeMantissaS1048 = _M0L4bitsS1046 & 4503599627370495ull;
  _M0L6_2atmpS2901 = _M0L4bitsS1046 >> 52;
  _M0L6_2atmpS2900 = _M0L6_2atmpS2901 & 2047ull;
  _M0L12ieeeExponentS1049 = (int32_t)_M0L6_2atmpS2900;
  if (_M0L12ieeeExponentS1049 == 2047) {
    _if__result_3943 = 1;
  } else if (_M0L12ieeeExponentS1049 == 0) {
    _if__result_3943 = _M0L12ieeeMantissaS1048 == 0ull;
  } else {
    _if__result_3943 = 0;
  }
  if (_if__result_3943) {
    int32_t _M0L6_2atmpS2891 = _M0L12ieeeExponentS1049 != 0;
    int32_t _M0L6_2atmpS2892 = _M0L12ieeeMantissaS1048 != 0ull;
    #line 676 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    return _M0FPB18copy__special__str(_M0L8ieeeSignS1047, _M0L6_2atmpS2891, _M0L6_2atmpS2892);
  }
  #line 678 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS1050
  = _M0FPB15d2d__small__int(_M0L12ieeeMantissaS1048, _M0L12ieeeExponentS1049);
  if (_M0L7_2abindS1050 == 0) {
    uint32_t _M0L6_2atmpS2893;
    if (_M0L7_2abindS1050) {
      moonbit_decref(_M0L7_2abindS1050);
    }
    _M0L6_2atmpS2893 = *(uint32_t*)&_M0L12ieeeExponentS1049;
    #line 688 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L1vS1051 = _M0FPB3d2d(_M0L12ieeeMantissaS1048, _M0L6_2atmpS2893);
  } else {
    struct _M0TPB17FloatingDecimal64* _M0L7_2aSomeS1052 = _M0L7_2abindS1050;
    struct _M0TPB17FloatingDecimal64* _M0L4_2afS1053 = _M0L7_2aSomeS1052;
    struct _M0TPB17FloatingDecimal64* _M0L1xS1054 = _M0L4_2afS1053;
    while (1) {
      uint64_t _M0L8mantissaS2899 = _M0L1xS1054->$0;
      uint64_t _M0L1qS1055 = _M0L8mantissaS2899 / 10ull;
      uint64_t _M0L8mantissaS2897 = _M0L1xS1054->$0;
      uint64_t _M0L6_2atmpS2898 = 10ull * _M0L1qS1055;
      uint64_t _M0L1rS1056 = _M0L8mantissaS2897 - _M0L6_2atmpS2898;
      int32_t _M0L8exponentS2896;
      int32_t _M0L6_2atmpS2895;
      struct _M0TPB17FloatingDecimal64* _M0L6_2atmpS2894;
      if (_M0L1rS1056 != 0ull) {
        _M0L1vS1051 = _M0L1xS1054;
        break;
      }
      _M0L8exponentS2896 = _M0L1xS1054->$1;
      moonbit_decref(_M0L1xS1054);
      _M0L6_2atmpS2895 = _M0L8exponentS2896 + 1;
      _M0L6_2atmpS2894
      = (struct _M0TPB17FloatingDecimal64*)moonbit_malloc(sizeof(struct _M0TPB17FloatingDecimal64));
      Moonbit_object_header(_M0L6_2atmpS2894)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L6_2atmpS2894->$0 = _M0L1qS1055;
      _M0L6_2atmpS2894->$1 = _M0L6_2atmpS2895;
      _M0L1xS1054 = _M0L6_2atmpS2894;
      continue;
      break;
    }
  }
  #line 690 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _result_3945 = _M0FPB9to__chars(_M0L1vS1051, _M0L8ieeeSignS1047);
  moonbit_decref(_M0L1vS1051);
  return _result_3945;
}

struct _M0TPB17FloatingDecimal64* _M0FPB15d2d__small__int(
  uint64_t _M0L12ieeeMantissaS1040,
  int32_t _M0L12ieeeExponentS1042
) {
  uint64_t _M0L2m2S1039;
  int32_t _M0L6_2atmpS2890;
  int32_t _M0L2e2S1041;
  int32_t _M0L6_2atmpS2889;
  uint64_t _M0L6_2atmpS2888;
  uint64_t _M0L4maskS1043;
  uint64_t _M0L8fractionS1044;
  int32_t _M0L6_2atmpS2887;
  uint64_t _M0L6_2atmpS2886;
  struct _M0TPB17FloatingDecimal64* _M0L6_2atmpS2885;
  #line 637 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L2m2S1039 = 4503599627370496ull | _M0L12ieeeMantissaS1040;
  _M0L6_2atmpS2890 = _M0L12ieeeExponentS1042 - 1023;
  _M0L2e2S1041 = _M0L6_2atmpS2890 - 52;
  if (_M0L2e2S1041 > 0) {
    return 0;
  }
  if (_M0L2e2S1041 < -52) {
    return 0;
  }
  _M0L6_2atmpS2889 = -_M0L2e2S1041;
  _M0L6_2atmpS2888 = 1ull << (_M0L6_2atmpS2889 & 63);
  _M0L4maskS1043 = _M0L6_2atmpS2888 - 1ull;
  _M0L8fractionS1044 = _M0L2m2S1039 & _M0L4maskS1043;
  if (_M0L8fractionS1044 != 0ull) {
    return 0;
  }
  _M0L6_2atmpS2887 = -_M0L2e2S1041;
  _M0L6_2atmpS2886 = _M0L2m2S1039 >> (_M0L6_2atmpS2887 & 63);
  _M0L6_2atmpS2885
  = (struct _M0TPB17FloatingDecimal64*)moonbit_malloc(sizeof(struct _M0TPB17FloatingDecimal64));
  Moonbit_object_header(_M0L6_2atmpS2885)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L6_2atmpS2885->$0 = _M0L6_2atmpS2886;
  _M0L6_2atmpS2885->$1 = 0;
  return _M0L6_2atmpS2885;
}

moonbit_string_t _M0FPB9to__chars(
  struct _M0TPB17FloatingDecimal64* _M0L1vS1007,
  int32_t _M0L4signS1005
) {
  int32_t _M0L6_2atmpS2884;
  moonbit_bytes_t _M0L6resultS1003;
  int32_t _M0Lm5indexS1004;
  uint64_t _M0L6outputS1006;
  int32_t _M0L7olengthS1008;
  int32_t _M0L8exponentS2883;
  int32_t _M0L6_2atmpS2882;
  int32_t _M0Lm3expS1009;
  int32_t _M0L6_2atmpS2881;
  int32_t _M0L6_2atmpS2879;
  int32_t _M0L18scientificNotationS1010;
  #line 530 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  #line 532 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2884 = _M0IPC14byte4BytePB7Default7default();
  _M0L6resultS1003
  = (moonbit_bytes_t)moonbit_make_bytes(25, _M0L6_2atmpS2884);
  _M0Lm5indexS1004 = 0;
  if (_M0L4signS1005) {
    int32_t _M0L6_2atmpS2753 = _M0Lm5indexS1004;
    int32_t _M0L6_2atmpS2754;
    if (
      _M0L6_2atmpS2753 < 0
      || _M0L6_2atmpS2753 >= Moonbit_array_length(_M0L6resultS1003)
    ) {
      #line 535 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      moonbit_panic();
    }
    _M0L6resultS1003[_M0L6_2atmpS2753] = 45;
    _M0L6_2atmpS2754 = _M0Lm5indexS1004;
    _M0Lm5indexS1004 = _M0L6_2atmpS2754 + 1;
  }
  _M0L6outputS1006 = _M0L1vS1007->$0;
  #line 539 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7olengthS1008 = _M0FPB17decimal__length17(_M0L6outputS1006);
  _M0L8exponentS2883 = _M0L1vS1007->$1;
  _M0L6_2atmpS2882 = _M0L8exponentS2883 + _M0L7olengthS1008;
  _M0Lm3expS1009 = _M0L6_2atmpS2882 - 1;
  _M0L6_2atmpS2881 = _M0Lm3expS1009;
  if (_M0L6_2atmpS2881 >= -6) {
    int32_t _M0L6_2atmpS2880 = _M0Lm3expS1009;
    _M0L6_2atmpS2879 = _M0L6_2atmpS2880 < 21;
  } else {
    _M0L6_2atmpS2879 = 0;
  }
  _M0L18scientificNotationS1010 = !_M0L6_2atmpS2879;
  if (_M0L18scientificNotationS1010) {
    int32_t _M0L7_2abindS1011 = _M0L7olengthS1008 - 1;
    uint64_t _M0L6outputS1012;
    int32_t _M0L1iS1013 = 0;
    uint64_t _M0L6outputS1014 = _M0L6outputS1006;
    int32_t _M0L6_2atmpS2755;
    int32_t _M0L6_2atmpS2759;
    int32_t _M0L6_2atmpS2758;
    int32_t _M0L6_2atmpS2757;
    int32_t _M0L6_2atmpS2756;
    int32_t _M0L6_2atmpS2763;
    int32_t _M0L6_2atmpS2764;
    int32_t _M0L6_2atmpS2765;
    int32_t _M0L6_2atmpS2766;
    int32_t _M0L6_2atmpS2767;
    int32_t _M0L6_2atmpS2773;
    int32_t _M0L6_2atmpS2806;
    moonbit_string_t _result_3947;
    while (1) {
      if (_M0L1iS1013 < _M0L7_2abindS1011) {
        uint64_t _M0L1cS1015 = _M0L6outputS1014 % 10ull;
        int32_t _M0L6_2atmpS2812 = _M0Lm5indexS1004;
        int32_t _M0L6_2atmpS2811 = _M0L6_2atmpS2812 + _M0L7olengthS1008;
        int32_t _M0L6_2atmpS2807 = _M0L6_2atmpS2811 - _M0L1iS1013;
        int32_t _M0L6_2atmpS2810 = (int32_t)_M0L1cS1015;
        int32_t _M0L6_2atmpS2809 = 48 + _M0L6_2atmpS2810;
        int32_t _M0L6_2atmpS2808 = _M0L6_2atmpS2809 & 0xff;
        int32_t _M0L6_2atmpS2813;
        uint64_t _M0L6_2atmpS2814;
        if (
          _M0L6_2atmpS2807 < 0
          || _M0L6_2atmpS2807 >= Moonbit_array_length(_M0L6resultS1003)
        ) {
          #line 547 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
          moonbit_panic();
        }
        _M0L6resultS1003[_M0L6_2atmpS2807] = _M0L6_2atmpS2808;
        _M0L6_2atmpS2813 = _M0L1iS1013 + 1;
        _M0L6_2atmpS2814 = _M0L6outputS1014 / 10ull;
        _M0L1iS1013 = _M0L6_2atmpS2813;
        _M0L6outputS1014 = _M0L6_2atmpS2814;
        continue;
      } else {
        _M0L6outputS1012 = _M0L6outputS1014;
      }
      break;
    }
    _M0L6_2atmpS2755 = _M0Lm5indexS1004;
    _M0L6_2atmpS2759 = (int32_t)_M0L6outputS1012;
    _M0L6_2atmpS2758 = _M0L6_2atmpS2759 % 10;
    _M0L6_2atmpS2757 = 48 + _M0L6_2atmpS2758;
    _M0L6_2atmpS2756 = _M0L6_2atmpS2757 & 0xff;
    if (
      _M0L6_2atmpS2755 < 0
      || _M0L6_2atmpS2755 >= Moonbit_array_length(_M0L6resultS1003)
    ) {
      #line 552 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      moonbit_panic();
    }
    _M0L6resultS1003[_M0L6_2atmpS2755] = _M0L6_2atmpS2756;
    if (_M0L7olengthS1008 > 1) {
      int32_t _M0L6_2atmpS2761 = _M0Lm5indexS1004;
      int32_t _M0L6_2atmpS2760 = _M0L6_2atmpS2761 + 1;
      if (
        _M0L6_2atmpS2760 < 0
        || _M0L6_2atmpS2760 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 554 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2760] = 46;
    } else {
      int32_t _M0L6_2atmpS2762 = _M0Lm5indexS1004;
      _M0Lm5indexS1004 = _M0L6_2atmpS2762 - 1;
    }
    _M0L6_2atmpS2763 = _M0Lm5indexS1004;
    _M0L6_2atmpS2764 = _M0L7olengthS1008 + 1;
    _M0Lm5indexS1004 = _M0L6_2atmpS2763 + _M0L6_2atmpS2764;
    _M0L6_2atmpS2765 = _M0Lm5indexS1004;
    if (
      _M0L6_2atmpS2765 < 0
      || _M0L6_2atmpS2765 >= Moonbit_array_length(_M0L6resultS1003)
    ) {
      #line 562 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      moonbit_panic();
    }
    _M0L6resultS1003[_M0L6_2atmpS2765] = 101;
    _M0L6_2atmpS2766 = _M0Lm5indexS1004;
    _M0Lm5indexS1004 = _M0L6_2atmpS2766 + 1;
    _M0L6_2atmpS2767 = _M0Lm3expS1009;
    if (_M0L6_2atmpS2767 < 0) {
      int32_t _M0L6_2atmpS2768 = _M0Lm5indexS1004;
      int32_t _M0L6_2atmpS2769;
      int32_t _M0L6_2atmpS2770;
      if (
        _M0L6_2atmpS2768 < 0
        || _M0L6_2atmpS2768 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 565 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2768] = 45;
      _M0L6_2atmpS2769 = _M0Lm5indexS1004;
      _M0Lm5indexS1004 = _M0L6_2atmpS2769 + 1;
      _M0L6_2atmpS2770 = _M0Lm3expS1009;
      _M0Lm3expS1009 = -_M0L6_2atmpS2770;
    } else {
      int32_t _M0L6_2atmpS2771 = _M0Lm5indexS1004;
      int32_t _M0L6_2atmpS2772;
      if (
        _M0L6_2atmpS2771 < 0
        || _M0L6_2atmpS2771 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 569 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2771] = 43;
      _M0L6_2atmpS2772 = _M0Lm5indexS1004;
      _M0Lm5indexS1004 = _M0L6_2atmpS2772 + 1;
    }
    _M0L6_2atmpS2773 = _M0Lm3expS1009;
    if (_M0L6_2atmpS2773 >= 100) {
      int32_t _M0L6_2atmpS2789 = _M0Lm3expS1009;
      int32_t _M0L1aS1017 = _M0L6_2atmpS2789 / 100;
      int32_t _M0L6_2atmpS2788 = _M0Lm3expS1009;
      int32_t _M0L6_2atmpS2787 = _M0L6_2atmpS2788 / 10;
      int32_t _M0L1bS1018 = _M0L6_2atmpS2787 % 10;
      int32_t _M0L6_2atmpS2786 = _M0Lm3expS1009;
      int32_t _M0L1cS1019 = _M0L6_2atmpS2786 % 10;
      int32_t _M0L6_2atmpS2774 = _M0Lm5indexS1004;
      int32_t _M0L6_2atmpS2776 = 48 + _M0L1aS1017;
      int32_t _M0L6_2atmpS2775 = _M0L6_2atmpS2776 & 0xff;
      int32_t _M0L6_2atmpS2780;
      int32_t _M0L6_2atmpS2777;
      int32_t _M0L6_2atmpS2779;
      int32_t _M0L6_2atmpS2778;
      int32_t _M0L6_2atmpS2784;
      int32_t _M0L6_2atmpS2781;
      int32_t _M0L6_2atmpS2783;
      int32_t _M0L6_2atmpS2782;
      int32_t _M0L6_2atmpS2785;
      if (
        _M0L6_2atmpS2774 < 0
        || _M0L6_2atmpS2774 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 576 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2774] = _M0L6_2atmpS2775;
      _M0L6_2atmpS2780 = _M0Lm5indexS1004;
      _M0L6_2atmpS2777 = _M0L6_2atmpS2780 + 1;
      _M0L6_2atmpS2779 = 48 + _M0L1bS1018;
      _M0L6_2atmpS2778 = _M0L6_2atmpS2779 & 0xff;
      if (
        _M0L6_2atmpS2777 < 0
        || _M0L6_2atmpS2777 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 577 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2777] = _M0L6_2atmpS2778;
      _M0L6_2atmpS2784 = _M0Lm5indexS1004;
      _M0L6_2atmpS2781 = _M0L6_2atmpS2784 + 2;
      _M0L6_2atmpS2783 = 48 + _M0L1cS1019;
      _M0L6_2atmpS2782 = _M0L6_2atmpS2783 & 0xff;
      if (
        _M0L6_2atmpS2781 < 0
        || _M0L6_2atmpS2781 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 578 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2781] = _M0L6_2atmpS2782;
      _M0L6_2atmpS2785 = _M0Lm5indexS1004;
      _M0Lm5indexS1004 = _M0L6_2atmpS2785 + 3;
    } else {
      int32_t _M0L6_2atmpS2790 = _M0Lm3expS1009;
      if (_M0L6_2atmpS2790 >= 10) {
        int32_t _M0L6_2atmpS2800 = _M0Lm3expS1009;
        int32_t _M0L1aS1020 = _M0L6_2atmpS2800 / 10;
        int32_t _M0L6_2atmpS2799 = _M0Lm3expS1009;
        int32_t _M0L1bS1021 = _M0L6_2atmpS2799 % 10;
        int32_t _M0L6_2atmpS2791 = _M0Lm5indexS1004;
        int32_t _M0L6_2atmpS2793 = 48 + _M0L1aS1020;
        int32_t _M0L6_2atmpS2792 = _M0L6_2atmpS2793 & 0xff;
        int32_t _M0L6_2atmpS2797;
        int32_t _M0L6_2atmpS2794;
        int32_t _M0L6_2atmpS2796;
        int32_t _M0L6_2atmpS2795;
        int32_t _M0L6_2atmpS2798;
        if (
          _M0L6_2atmpS2791 < 0
          || _M0L6_2atmpS2791 >= Moonbit_array_length(_M0L6resultS1003)
        ) {
          #line 583 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
          moonbit_panic();
        }
        _M0L6resultS1003[_M0L6_2atmpS2791] = _M0L6_2atmpS2792;
        _M0L6_2atmpS2797 = _M0Lm5indexS1004;
        _M0L6_2atmpS2794 = _M0L6_2atmpS2797 + 1;
        _M0L6_2atmpS2796 = 48 + _M0L1bS1021;
        _M0L6_2atmpS2795 = _M0L6_2atmpS2796 & 0xff;
        if (
          _M0L6_2atmpS2794 < 0
          || _M0L6_2atmpS2794 >= Moonbit_array_length(_M0L6resultS1003)
        ) {
          #line 584 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
          moonbit_panic();
        }
        _M0L6resultS1003[_M0L6_2atmpS2794] = _M0L6_2atmpS2795;
        _M0L6_2atmpS2798 = _M0Lm5indexS1004;
        _M0Lm5indexS1004 = _M0L6_2atmpS2798 + 2;
      } else {
        int32_t _M0L6_2atmpS2801 = _M0Lm5indexS1004;
        int32_t _M0L6_2atmpS2804 = _M0Lm3expS1009;
        int32_t _M0L6_2atmpS2803 = 48 + _M0L6_2atmpS2804;
        int32_t _M0L6_2atmpS2802 = _M0L6_2atmpS2803 & 0xff;
        int32_t _M0L6_2atmpS2805;
        if (
          _M0L6_2atmpS2801 < 0
          || _M0L6_2atmpS2801 >= Moonbit_array_length(_M0L6resultS1003)
        ) {
          #line 587 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
          moonbit_panic();
        }
        _M0L6resultS1003[_M0L6_2atmpS2801] = _M0L6_2atmpS2802;
        _M0L6_2atmpS2805 = _M0Lm5indexS1004;
        _M0Lm5indexS1004 = _M0L6_2atmpS2805 + 1;
      }
    }
    _M0L6_2atmpS2806 = _M0Lm5indexS1004;
    #line 590 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _result_3947
    = _M0FPB19string__from__bytes(_M0L6resultS1003, 0, _M0L6_2atmpS2806);
    moonbit_decref(_M0L6resultS1003);
    return _result_3947;
  } else {
    int32_t _M0L6_2atmpS2815 = _M0Lm3expS1009;
    int32_t _M0L6_2atmpS2878;
    moonbit_string_t _result_3953;
    if (_M0L6_2atmpS2815 < 0) {
      int32_t _M0L6_2atmpS2816 = _M0Lm5indexS1004;
      int32_t _M0L6_2atmpS2818;
      int32_t _M0L6_2atmpS2817;
      int32_t _M0L6_2atmpS2819;
      int32_t _M0L1iS1022;
      int32_t _M0L6_2atmpS2834;
      int32_t _M0L6_2atmpS2836;
      int32_t _M0L6_2atmpS2835;
      int32_t _M0L7currentS1024;
      int32_t _M0L1iS1025;
      uint64_t _M0L6outputS1026;
      if (
        _M0L6_2atmpS2816 < 0
        || _M0L6_2atmpS2816 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 595 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2816] = 48;
      _M0L6_2atmpS2818 = _M0Lm5indexS1004;
      _M0L6_2atmpS2817 = _M0L6_2atmpS2818 + 1;
      if (
        _M0L6_2atmpS2817 < 0
        || _M0L6_2atmpS2817 >= Moonbit_array_length(_M0L6resultS1003)
      ) {
        #line 596 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6resultS1003[_M0L6_2atmpS2817] = 46;
      _M0L6_2atmpS2819 = _M0Lm5indexS1004;
      _M0Lm5indexS1004 = _M0L6_2atmpS2819 + 2;
      _M0L1iS1022 = -1;
      while (1) {
        int32_t _M0L6_2atmpS2820 = _M0Lm3expS1009;
        if (_M0L1iS1022 > _M0L6_2atmpS2820) {
          int32_t _M0L6_2atmpS2823 = _M0Lm5indexS1004;
          int32_t _M0L6_2atmpS2822 = _M0L6_2atmpS2823 - _M0L1iS1022;
          int32_t _M0L6_2atmpS2821 = _M0L6_2atmpS2822 - 1;
          int32_t _M0L6_2atmpS2824;
          if (
            _M0L6_2atmpS2821 < 0
            || _M0L6_2atmpS2821 >= Moonbit_array_length(_M0L6resultS1003)
          ) {
            #line 599 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
            moonbit_panic();
          }
          _M0L6resultS1003[_M0L6_2atmpS2821] = 48;
          _M0L6_2atmpS2824 = _M0L1iS1022 - 1;
          _M0L1iS1022 = _M0L6_2atmpS2824;
          continue;
        }
        break;
      }
      _M0L6_2atmpS2834 = _M0Lm5indexS1004;
      _M0L6_2atmpS2836 = _M0Lm3expS1009;
      _M0L6_2atmpS2835 = -1 - _M0L6_2atmpS2836;
      _M0L7currentS1024 = _M0L6_2atmpS2834 + _M0L6_2atmpS2835;
      _M0L1iS1025 = 0;
      _M0L6outputS1026 = _M0L6outputS1006;
      while (1) {
        if (_M0L1iS1025 < _M0L7olengthS1008) {
          int32_t _M0L6_2atmpS2831 = _M0L7currentS1024 + _M0L7olengthS1008;
          int32_t _M0L6_2atmpS2830 = _M0L6_2atmpS2831 - _M0L1iS1025;
          int32_t _M0L6_2atmpS2825 = _M0L6_2atmpS2830 - 1;
          uint64_t _M0L6_2atmpS2829 = _M0L6outputS1026 % 10ull;
          int32_t _M0L6_2atmpS2828 = (int32_t)_M0L6_2atmpS2829;
          int32_t _M0L6_2atmpS2827 = 48 + _M0L6_2atmpS2828;
          int32_t _M0L6_2atmpS2826 = _M0L6_2atmpS2827 & 0xff;
          int32_t _M0L6_2atmpS2832;
          uint64_t _M0L6_2atmpS2833;
          if (
            _M0L6_2atmpS2825 < 0
            || _M0L6_2atmpS2825 >= Moonbit_array_length(_M0L6resultS1003)
          ) {
            #line 603 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
            moonbit_panic();
          }
          _M0L6resultS1003[_M0L6_2atmpS2825] = _M0L6_2atmpS2826;
          _M0L6_2atmpS2832 = _M0L1iS1025 + 1;
          _M0L6_2atmpS2833 = _M0L6outputS1026 / 10ull;
          _M0L1iS1025 = _M0L6_2atmpS2832;
          _M0L6outputS1026 = _M0L6_2atmpS2833;
          continue;
        }
        break;
      }
      _M0Lm5indexS1004 = _M0L7currentS1024 + _M0L7olengthS1008;
    } else {
      int32_t _M0L6_2atmpS2838 = _M0Lm3expS1009;
      int32_t _M0L6_2atmpS2837 = _M0L6_2atmpS2838 + 1;
      if (_M0L6_2atmpS2837 >= _M0L7olengthS1008) {
        int32_t _M0L1iS1028 = 0;
        uint64_t _M0L6outputS1029 = _M0L6outputS1006;
        int32_t _M0L6_2atmpS2849;
        int32_t _M0L6_2atmpS2854;
        int32_t _M0L7_2abindS1031;
        int32_t _M0L1iS1032;
        int32_t _M0L6_2atmpS2855;
        int32_t _M0L6_2atmpS2858;
        int32_t _M0L6_2atmpS2857;
        int32_t _M0L6_2atmpS2856;
        while (1) {
          if (_M0L1iS1028 < _M0L7olengthS1008) {
            int32_t _M0L6_2atmpS2846 = _M0Lm5indexS1004;
            int32_t _M0L6_2atmpS2845 = _M0L6_2atmpS2846 + _M0L7olengthS1008;
            int32_t _M0L6_2atmpS2844 = _M0L6_2atmpS2845 - _M0L1iS1028;
            int32_t _M0L6_2atmpS2839 = _M0L6_2atmpS2844 - 1;
            uint64_t _M0L6_2atmpS2843 = _M0L6outputS1029 % 10ull;
            int32_t _M0L6_2atmpS2842 = (int32_t)_M0L6_2atmpS2843;
            int32_t _M0L6_2atmpS2841 = 48 + _M0L6_2atmpS2842;
            int32_t _M0L6_2atmpS2840 = _M0L6_2atmpS2841 & 0xff;
            int32_t _M0L6_2atmpS2847;
            uint64_t _M0L6_2atmpS2848;
            if (
              _M0L6_2atmpS2839 < 0
              || _M0L6_2atmpS2839 >= Moonbit_array_length(_M0L6resultS1003)
            ) {
              #line 610 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
              moonbit_panic();
            }
            _M0L6resultS1003[_M0L6_2atmpS2839] = _M0L6_2atmpS2840;
            _M0L6_2atmpS2847 = _M0L1iS1028 + 1;
            _M0L6_2atmpS2848 = _M0L6outputS1029 / 10ull;
            _M0L1iS1028 = _M0L6_2atmpS2847;
            _M0L6outputS1029 = _M0L6_2atmpS2848;
            continue;
          }
          break;
        }
        _M0L6_2atmpS2849 = _M0Lm5indexS1004;
        _M0Lm5indexS1004 = _M0L6_2atmpS2849 + _M0L7olengthS1008;
        _M0L6_2atmpS2854 = _M0Lm3expS1009;
        _M0L7_2abindS1031 = _M0L6_2atmpS2854 + 1;
        _M0L1iS1032 = _M0L7olengthS1008;
        while (1) {
          if (_M0L1iS1032 < _M0L7_2abindS1031) {
            int32_t _M0L6_2atmpS2852 = _M0Lm5indexS1004;
            int32_t _M0L6_2atmpS2851 = _M0L6_2atmpS2852 + _M0L1iS1032;
            int32_t _M0L6_2atmpS2850 = _M0L6_2atmpS2851 - _M0L7olengthS1008;
            int32_t _M0L6_2atmpS2853;
            if (
              _M0L6_2atmpS2850 < 0
              || _M0L6_2atmpS2850 >= Moonbit_array_length(_M0L6resultS1003)
            ) {
              #line 615 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
              moonbit_panic();
            }
            _M0L6resultS1003[_M0L6_2atmpS2850] = 48;
            _M0L6_2atmpS2853 = _M0L1iS1032 + 1;
            _M0L1iS1032 = _M0L6_2atmpS2853;
            continue;
          }
          break;
        }
        _M0L6_2atmpS2855 = _M0Lm5indexS1004;
        _M0L6_2atmpS2858 = _M0Lm3expS1009;
        _M0L6_2atmpS2857 = _M0L6_2atmpS2858 + 1;
        _M0L6_2atmpS2856 = _M0L6_2atmpS2857 - _M0L7olengthS1008;
        _M0Lm5indexS1004 = _M0L6_2atmpS2855 + _M0L6_2atmpS2856;
      } else {
        int32_t _M0L6_2atmpS2875 = _M0Lm5indexS1004;
        int32_t _M0L6_2atmpS2874 = _M0L6_2atmpS2875 + 1;
        int32_t _M0L1iS1034 = 0;
        int32_t _M0L7currentS1035 = _M0L6_2atmpS2874;
        uint64_t _M0L6outputS1036 = _M0L6outputS1006;
        int32_t _M0L6_2atmpS2876;
        int32_t _M0L6_2atmpS2877;
        while (1) {
          if (_M0L1iS1034 < _M0L7olengthS1008) {
            int32_t _M0L6_2atmpS2870 = _M0L7olengthS1008 - _M0L1iS1034;
            int32_t _M0L6_2atmpS2868 = _M0L6_2atmpS2870 - 1;
            int32_t _M0L6_2atmpS2869 = _M0Lm3expS1009;
            int32_t _M0L7currentS1037;
            int32_t _M0L6_2atmpS2865;
            int32_t _M0L6_2atmpS2864;
            int32_t _M0L6_2atmpS2859;
            uint64_t _M0L6_2atmpS2863;
            int32_t _M0L6_2atmpS2862;
            int32_t _M0L6_2atmpS2861;
            int32_t _M0L6_2atmpS2860;
            int32_t _M0L6_2atmpS2866;
            uint64_t _M0L6_2atmpS2867;
            if (_M0L6_2atmpS2868 == _M0L6_2atmpS2869) {
              int32_t _M0L6_2atmpS2873 =
                _M0L7currentS1035 + _M0L7olengthS1008;
              int32_t _M0L6_2atmpS2872 = _M0L6_2atmpS2873 - _M0L1iS1034;
              int32_t _M0L6_2atmpS2871 = _M0L6_2atmpS2872 - 1;
              if (
                _M0L6_2atmpS2871 < 0
                || _M0L6_2atmpS2871 >= Moonbit_array_length(_M0L6resultS1003)
              ) {
                #line 622 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
                moonbit_panic();
              }
              _M0L6resultS1003[_M0L6_2atmpS2871] = 46;
              _M0L7currentS1037 = _M0L7currentS1035 - 1;
            } else {
              _M0L7currentS1037 = _M0L7currentS1035;
            }
            _M0L6_2atmpS2865 = _M0L7currentS1037 + _M0L7olengthS1008;
            _M0L6_2atmpS2864 = _M0L6_2atmpS2865 - _M0L1iS1034;
            _M0L6_2atmpS2859 = _M0L6_2atmpS2864 - 1;
            _M0L6_2atmpS2863 = _M0L6outputS1036 % 10ull;
            _M0L6_2atmpS2862 = (int32_t)_M0L6_2atmpS2863;
            _M0L6_2atmpS2861 = 48 + _M0L6_2atmpS2862;
            _M0L6_2atmpS2860 = _M0L6_2atmpS2861 & 0xff;
            if (
              _M0L6_2atmpS2859 < 0
              || _M0L6_2atmpS2859 >= Moonbit_array_length(_M0L6resultS1003)
            ) {
              #line 627 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
              moonbit_panic();
            }
            _M0L6resultS1003[_M0L6_2atmpS2859] = _M0L6_2atmpS2860;
            _M0L6_2atmpS2866 = _M0L1iS1034 + 1;
            _M0L6_2atmpS2867 = _M0L6outputS1036 / 10ull;
            _M0L1iS1034 = _M0L6_2atmpS2866;
            _M0L7currentS1035 = _M0L7currentS1037;
            _M0L6outputS1036 = _M0L6_2atmpS2867;
            continue;
          }
          break;
        }
        _M0L6_2atmpS2876 = _M0Lm5indexS1004;
        _M0L6_2atmpS2877 = _M0L7olengthS1008 + 1;
        _M0Lm5indexS1004 = _M0L6_2atmpS2876 + _M0L6_2atmpS2877;
      }
    }
    _M0L6_2atmpS2878 = _M0Lm5indexS1004;
    #line 632 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _result_3953
    = _M0FPB19string__from__bytes(_M0L6resultS1003, 0, _M0L6_2atmpS2878);
    moonbit_decref(_M0L6resultS1003);
    return _result_3953;
  }
}

struct _M0TPB17FloatingDecimal64* _M0FPB3d2d(
  uint64_t _M0L12ieeeMantissaS949,
  uint32_t _M0L12ieeeExponentS948
) {
  int32_t _M0Lm2e2S946;
  uint64_t _M0Lm2m2S947;
  uint64_t _M0L6_2atmpS2752;
  uint64_t _M0L6_2atmpS2751;
  int32_t _M0L4evenS950;
  uint64_t _M0L6_2atmpS2750;
  uint64_t _M0L2mvS951;
  int32_t _M0L7mmShiftS952;
  uint64_t _M0Lm2vrS953;
  uint64_t _M0Lm2vpS954;
  uint64_t _M0Lm2vmS955;
  int32_t _M0Lm3e10S956;
  int32_t _M0Lm17vmIsTrailingZerosS957;
  int32_t _M0Lm17vrIsTrailingZerosS958;
  int32_t _M0L6_2atmpS2652;
  int32_t _M0Lm7removedS977;
  int32_t _M0Lm16lastRemovedDigitS978;
  uint64_t _M0Lm6outputS979;
  int32_t _M0L6_2atmpS2748;
  int32_t _M0L6_2atmpS2749;
  int32_t _M0L3expS1002;
  uint64_t _M0L6_2atmpS2747;
  struct _M0TPB17FloatingDecimal64* _block_3959;
  #line 347 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0Lm2e2S946 = 0;
  _M0Lm2m2S947 = 0ull;
  if (_M0L12ieeeExponentS948 == 0u) {
    _M0Lm2e2S946 = -1076;
    _M0Lm2m2S947 = _M0L12ieeeMantissaS949;
  } else {
    int32_t _M0L6_2atmpS2651 = *(int32_t*)&_M0L12ieeeExponentS948;
    int32_t _M0L6_2atmpS2650 = _M0L6_2atmpS2651 - 1023;
    int32_t _M0L6_2atmpS2649 = _M0L6_2atmpS2650 - 52;
    _M0Lm2e2S946 = _M0L6_2atmpS2649 - 2;
    _M0Lm2m2S947 = 4503599627370496ull | _M0L12ieeeMantissaS949;
  }
  _M0L6_2atmpS2752 = _M0Lm2m2S947;
  _M0L6_2atmpS2751 = _M0L6_2atmpS2752 & 1ull;
  _M0L4evenS950 = _M0L6_2atmpS2751 == 0ull;
  _M0L6_2atmpS2750 = _M0Lm2m2S947;
  _M0L2mvS951 = 4ull * _M0L6_2atmpS2750;
  if (_M0L12ieeeMantissaS949 != 0ull) {
    _M0L7mmShiftS952 = 1;
  } else {
    _M0L7mmShiftS952 = _M0L12ieeeExponentS948 <= 1u;
  }
  _M0Lm2vrS953 = 0ull;
  _M0Lm2vpS954 = 0ull;
  _M0Lm2vmS955 = 0ull;
  _M0Lm3e10S956 = 0;
  _M0Lm17vmIsTrailingZerosS957 = 0;
  _M0Lm17vrIsTrailingZerosS958 = 0;
  _M0L6_2atmpS2652 = _M0Lm2e2S946;
  if (_M0L6_2atmpS2652 >= 0) {
    int32_t _M0L6_2atmpS2674 = _M0Lm2e2S946;
    int32_t _M0L6_2atmpS2670;
    int32_t _M0L6_2atmpS2673;
    int32_t _M0L6_2atmpS2672;
    int32_t _M0L6_2atmpS2671;
    int32_t _M0L1qS959;
    int32_t _M0L6_2atmpS2669;
    int32_t _M0L6_2atmpS2668;
    int32_t _M0L1kS960;
    int32_t _M0L6_2atmpS2667;
    int32_t _M0L6_2atmpS2666;
    int32_t _M0L6_2atmpS2665;
    int32_t _M0L1iS961;
    struct _M0TPB8Pow5Pair _M0L4pow5S962;
    uint64_t _M0L6_2atmpS2664;
    struct _M0TPB19MulShiftAll64Result _M0L7_2abindS963;
    uint64_t _M0L8_2avrOutS964;
    uint64_t _M0L8_2avpOutS965;
    uint64_t _M0L8_2avmOutS966;
    #line 383 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2670 = _M0FPB9log10Pow2(_M0L6_2atmpS2674);
    _M0L6_2atmpS2673 = _M0Lm2e2S946;
    _M0L6_2atmpS2672 = _M0L6_2atmpS2673 > 3;
    #line 383 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2671 = _M0MPC14bool4Bool7to__int(_M0L6_2atmpS2672);
    _M0L1qS959 = _M0L6_2atmpS2670 - _M0L6_2atmpS2671;
    _M0Lm3e10S956 = _M0L1qS959;
    #line 385 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2669 = _M0FPB8pow5bits(_M0L1qS959);
    _M0L6_2atmpS2668 = 125 + _M0L6_2atmpS2669;
    _M0L1kS960 = _M0L6_2atmpS2668 - 1;
    _M0L6_2atmpS2667 = _M0Lm2e2S946;
    _M0L6_2atmpS2666 = -_M0L6_2atmpS2667;
    _M0L6_2atmpS2665 = _M0L6_2atmpS2666 + _M0L1qS959;
    _M0L1iS961 = _M0L6_2atmpS2665 + _M0L1kS960;
    #line 387 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L4pow5S962 = _M0FPB22double__computeInvPow5(_M0L1qS959);
    _M0L6_2atmpS2664 = _M0Lm2m2S947;
    #line 388 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L7_2abindS963
    = _M0FPB13mulShiftAll64(_M0L6_2atmpS2664, _M0L4pow5S962, _M0L1iS961, _M0L7mmShiftS952);
    _M0L8_2avrOutS964 = _M0L7_2abindS963.$0;
    _M0L8_2avpOutS965 = _M0L7_2abindS963.$1;
    _M0L8_2avmOutS966 = _M0L7_2abindS963.$2;
    _M0Lm2vrS953 = _M0L8_2avrOutS964;
    _M0Lm2vpS954 = _M0L8_2avpOutS965;
    _M0Lm2vmS955 = _M0L8_2avmOutS966;
    if (_M0L1qS959 <= 21) {
      int32_t _M0L6_2atmpS2660 = (int32_t)_M0L2mvS951;
      uint64_t _M0L6_2atmpS2663 = _M0L2mvS951 / 5ull;
      int32_t _M0L6_2atmpS2662 = (int32_t)_M0L6_2atmpS2663;
      int32_t _M0L6_2atmpS2661 = 5 * _M0L6_2atmpS2662;
      int32_t _M0L6mvMod5S967 = _M0L6_2atmpS2660 - _M0L6_2atmpS2661;
      if (_M0L6mvMod5S967 == 0) {
        #line 400 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        _M0Lm17vrIsTrailingZerosS958
        = _M0FPB18multipleOfPowerOf5(_M0L2mvS951, _M0L1qS959);
      } else if (_M0L4evenS950) {
        uint64_t _M0L6_2atmpS2654 = _M0L2mvS951 - 1ull;
        uint64_t _M0L6_2atmpS2655;
        uint64_t _M0L6_2atmpS2653;
        #line 406 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        _M0L6_2atmpS2655 = _M0MPC14bool4Bool10to__uint64(_M0L7mmShiftS952);
        _M0L6_2atmpS2653 = _M0L6_2atmpS2654 - _M0L6_2atmpS2655;
        #line 405 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        _M0Lm17vmIsTrailingZerosS957
        = _M0FPB18multipleOfPowerOf5(_M0L6_2atmpS2653, _M0L1qS959);
      } else {
        uint64_t _M0L6_2atmpS2656 = _M0Lm2vpS954;
        uint64_t _M0L6_2atmpS2659 = _M0L2mvS951 + 2ull;
        int32_t _M0L6_2atmpS2658;
        uint64_t _M0L6_2atmpS2657;
        #line 410 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        _M0L6_2atmpS2658
        = _M0FPB18multipleOfPowerOf5(_M0L6_2atmpS2659, _M0L1qS959);
        #line 410 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        _M0L6_2atmpS2657 = _M0MPC14bool4Bool10to__uint64(_M0L6_2atmpS2658);
        _M0Lm2vpS954 = _M0L6_2atmpS2656 - _M0L6_2atmpS2657;
      }
    }
  } else {
    int32_t _M0L6_2atmpS2688 = _M0Lm2e2S946;
    int32_t _M0L6_2atmpS2687 = -_M0L6_2atmpS2688;
    int32_t _M0L6_2atmpS2682;
    int32_t _M0L6_2atmpS2686;
    int32_t _M0L6_2atmpS2685;
    int32_t _M0L6_2atmpS2684;
    int32_t _M0L6_2atmpS2683;
    int32_t _M0L1qS968;
    int32_t _M0L6_2atmpS2675;
    int32_t _M0L6_2atmpS2681;
    int32_t _M0L6_2atmpS2680;
    int32_t _M0L1iS969;
    int32_t _M0L6_2atmpS2679;
    int32_t _M0L1kS970;
    int32_t _M0L1jS971;
    struct _M0TPB8Pow5Pair _M0L4pow5S972;
    uint64_t _M0L6_2atmpS2678;
    struct _M0TPB19MulShiftAll64Result _M0L7_2abindS973;
    uint64_t _M0L8_2avrOutS974;
    uint64_t _M0L8_2avpOutS975;
    uint64_t _M0L8_2avmOutS976;
    #line 415 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2682 = _M0FPB9log10Pow5(_M0L6_2atmpS2687);
    _M0L6_2atmpS2686 = _M0Lm2e2S946;
    _M0L6_2atmpS2685 = -_M0L6_2atmpS2686;
    _M0L6_2atmpS2684 = _M0L6_2atmpS2685 > 1;
    #line 415 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2683 = _M0MPC14bool4Bool7to__int(_M0L6_2atmpS2684);
    _M0L1qS968 = _M0L6_2atmpS2682 - _M0L6_2atmpS2683;
    _M0L6_2atmpS2675 = _M0Lm2e2S946;
    _M0Lm3e10S956 = _M0L1qS968 + _M0L6_2atmpS2675;
    _M0L6_2atmpS2681 = _M0Lm2e2S946;
    _M0L6_2atmpS2680 = -_M0L6_2atmpS2681;
    _M0L1iS969 = _M0L6_2atmpS2680 - _M0L1qS968;
    #line 418 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2679 = _M0FPB8pow5bits(_M0L1iS969);
    _M0L1kS970 = _M0L6_2atmpS2679 - 125;
    _M0L1jS971 = _M0L1qS968 - _M0L1kS970;
    #line 420 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L4pow5S972 = _M0FPB19double__computePow5(_M0L1iS969);
    _M0L6_2atmpS2678 = _M0Lm2m2S947;
    #line 421 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L7_2abindS973
    = _M0FPB13mulShiftAll64(_M0L6_2atmpS2678, _M0L4pow5S972, _M0L1jS971, _M0L7mmShiftS952);
    _M0L8_2avrOutS974 = _M0L7_2abindS973.$0;
    _M0L8_2avpOutS975 = _M0L7_2abindS973.$1;
    _M0L8_2avmOutS976 = _M0L7_2abindS973.$2;
    _M0Lm2vrS953 = _M0L8_2avrOutS974;
    _M0Lm2vpS954 = _M0L8_2avpOutS975;
    _M0Lm2vmS955 = _M0L8_2avmOutS976;
    if (_M0L1qS968 <= 1) {
      _M0Lm17vrIsTrailingZerosS958 = 1;
      if (_M0L4evenS950) {
        int32_t _M0L6_2atmpS2676;
        #line 432 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        _M0L6_2atmpS2676 = _M0MPC14bool4Bool7to__int(_M0L7mmShiftS952);
        _M0Lm17vmIsTrailingZerosS957 = _M0L6_2atmpS2676 == 1;
      } else {
        uint64_t _M0L6_2atmpS2677 = _M0Lm2vpS954;
        _M0Lm2vpS954 = _M0L6_2atmpS2677 - 1ull;
      }
    } else if (_M0L1qS968 < 63) {
      #line 437 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      _M0Lm17vrIsTrailingZerosS958
      = _M0FPB18multipleOfPowerOf2(_M0L2mvS951, _M0L1qS968);
    }
  }
  _M0Lm7removedS977 = 0;
  _M0Lm16lastRemovedDigitS978 = 0;
  _M0Lm6outputS979 = 0ull;
  if (_M0Lm17vmIsTrailingZerosS957 || _M0Lm17vrIsTrailingZerosS958) {
    int32_t _if__result_3956;
    uint64_t _M0L6_2atmpS2718;
    uint64_t _M0L6_2atmpS2724;
    uint64_t _M0L6_2atmpS2725;
    int32_t _if__result_3957;
    int32_t _M0L6_2atmpS2721;
    int64_t _M0L6_2atmpS2720;
    uint64_t _M0L6_2atmpS2719;
    while (1) {
      uint64_t _M0L6_2atmpS2701 = _M0Lm2vpS954;
      uint64_t _M0L7vpDiv10S980 = _M0L6_2atmpS2701 / 10ull;
      uint64_t _M0L6_2atmpS2700 = _M0Lm2vmS955;
      uint64_t _M0L7vmDiv10S981 = _M0L6_2atmpS2700 / 10ull;
      uint64_t _M0L6_2atmpS2699;
      int32_t _M0L6_2atmpS2696;
      int32_t _M0L6_2atmpS2698;
      int32_t _M0L6_2atmpS2697;
      int32_t _M0L7vmMod10S983;
      uint64_t _M0L6_2atmpS2695;
      uint64_t _M0L7vrDiv10S984;
      uint64_t _M0L6_2atmpS2694;
      int32_t _M0L6_2atmpS2691;
      int32_t _M0L6_2atmpS2693;
      int32_t _M0L6_2atmpS2692;
      int32_t _M0L7vrMod10S985;
      int32_t _M0L6_2atmpS2690;
      if (_M0L7vpDiv10S980 <= _M0L7vmDiv10S981) {
        break;
      }
      _M0L6_2atmpS2699 = _M0Lm2vmS955;
      _M0L6_2atmpS2696 = (int32_t)_M0L6_2atmpS2699;
      _M0L6_2atmpS2698 = (int32_t)_M0L7vmDiv10S981;
      _M0L6_2atmpS2697 = 10 * _M0L6_2atmpS2698;
      _M0L7vmMod10S983 = _M0L6_2atmpS2696 - _M0L6_2atmpS2697;
      _M0L6_2atmpS2695 = _M0Lm2vrS953;
      _M0L7vrDiv10S984 = _M0L6_2atmpS2695 / 10ull;
      _M0L6_2atmpS2694 = _M0Lm2vrS953;
      _M0L6_2atmpS2691 = (int32_t)_M0L6_2atmpS2694;
      _M0L6_2atmpS2693 = (int32_t)_M0L7vrDiv10S984;
      _M0L6_2atmpS2692 = 10 * _M0L6_2atmpS2693;
      _M0L7vrMod10S985 = _M0L6_2atmpS2691 - _M0L6_2atmpS2692;
      if (_M0Lm17vmIsTrailingZerosS957) {
        _M0Lm17vmIsTrailingZerosS957 = _M0L7vmMod10S983 == 0;
      } else {
        _M0Lm17vmIsTrailingZerosS957 = 0;
      }
      if (_M0Lm17vrIsTrailingZerosS958) {
        int32_t _M0L6_2atmpS2689 = _M0Lm16lastRemovedDigitS978;
        _M0Lm17vrIsTrailingZerosS958 = _M0L6_2atmpS2689 == 0;
      } else {
        _M0Lm17vrIsTrailingZerosS958 = 0;
      }
      _M0Lm16lastRemovedDigitS978 = _M0L7vrMod10S985;
      _M0Lm2vrS953 = _M0L7vrDiv10S984;
      _M0Lm2vpS954 = _M0L7vpDiv10S980;
      _M0Lm2vmS955 = _M0L7vmDiv10S981;
      _M0L6_2atmpS2690 = _M0Lm7removedS977;
      _M0Lm7removedS977 = _M0L6_2atmpS2690 + 1;
      continue;
      break;
    }
    if (_M0Lm17vmIsTrailingZerosS957) {
      while (1) {
        uint64_t _M0L6_2atmpS2714 = _M0Lm2vmS955;
        uint64_t _M0L7vmDiv10S986 = _M0L6_2atmpS2714 / 10ull;
        uint64_t _M0L6_2atmpS2713 = _M0Lm2vmS955;
        int32_t _M0L6_2atmpS2710 = (int32_t)_M0L6_2atmpS2713;
        int32_t _M0L6_2atmpS2712 = (int32_t)_M0L7vmDiv10S986;
        int32_t _M0L6_2atmpS2711 = 10 * _M0L6_2atmpS2712;
        int32_t _M0L7vmMod10S987 = _M0L6_2atmpS2710 - _M0L6_2atmpS2711;
        uint64_t _M0L6_2atmpS2709;
        uint64_t _M0L7vpDiv10S989;
        uint64_t _M0L6_2atmpS2708;
        uint64_t _M0L7vrDiv10S990;
        uint64_t _M0L6_2atmpS2707;
        int32_t _M0L6_2atmpS2704;
        int32_t _M0L6_2atmpS2706;
        int32_t _M0L6_2atmpS2705;
        int32_t _M0L7vrMod10S991;
        int32_t _M0L6_2atmpS2703;
        if (_M0L7vmMod10S987 != 0) {
          break;
        }
        _M0L6_2atmpS2709 = _M0Lm2vpS954;
        _M0L7vpDiv10S989 = _M0L6_2atmpS2709 / 10ull;
        _M0L6_2atmpS2708 = _M0Lm2vrS953;
        _M0L7vrDiv10S990 = _M0L6_2atmpS2708 / 10ull;
        _M0L6_2atmpS2707 = _M0Lm2vrS953;
        _M0L6_2atmpS2704 = (int32_t)_M0L6_2atmpS2707;
        _M0L6_2atmpS2706 = (int32_t)_M0L7vrDiv10S990;
        _M0L6_2atmpS2705 = 10 * _M0L6_2atmpS2706;
        _M0L7vrMod10S991 = _M0L6_2atmpS2704 - _M0L6_2atmpS2705;
        if (_M0Lm17vrIsTrailingZerosS958) {
          int32_t _M0L6_2atmpS2702 = _M0Lm16lastRemovedDigitS978;
          _M0Lm17vrIsTrailingZerosS958 = _M0L6_2atmpS2702 == 0;
        } else {
          _M0Lm17vrIsTrailingZerosS958 = 0;
        }
        _M0Lm16lastRemovedDigitS978 = _M0L7vrMod10S991;
        _M0Lm2vrS953 = _M0L7vrDiv10S990;
        _M0Lm2vpS954 = _M0L7vpDiv10S989;
        _M0Lm2vmS955 = _M0L7vmDiv10S986;
        _M0L6_2atmpS2703 = _M0Lm7removedS977;
        _M0Lm7removedS977 = _M0L6_2atmpS2703 + 1;
        continue;
        break;
      }
    }
    if (_M0Lm17vrIsTrailingZerosS958) {
      int32_t _M0L6_2atmpS2717 = _M0Lm16lastRemovedDigitS978;
      if (_M0L6_2atmpS2717 == 5) {
        uint64_t _M0L6_2atmpS2716 = _M0Lm2vrS953;
        uint64_t _M0L6_2atmpS2715 = _M0L6_2atmpS2716 % 2ull;
        _if__result_3956 = _M0L6_2atmpS2715 == 0ull;
      } else {
        _if__result_3956 = 0;
      }
    } else {
      _if__result_3956 = 0;
    }
    if (_if__result_3956) {
      _M0Lm16lastRemovedDigitS978 = 4;
    }
    _M0L6_2atmpS2718 = _M0Lm2vrS953;
    _M0L6_2atmpS2724 = _M0Lm2vrS953;
    _M0L6_2atmpS2725 = _M0Lm2vmS955;
    if (_M0L6_2atmpS2724 == _M0L6_2atmpS2725) {
      if (!_M0L4evenS950) {
        _if__result_3957 = 1;
      } else {
        int32_t _M0L6_2atmpS2723 = _M0Lm17vmIsTrailingZerosS957;
        _if__result_3957 = !_M0L6_2atmpS2723;
      }
    } else {
      _if__result_3957 = 0;
    }
    if (_if__result_3957) {
      _M0L6_2atmpS2721 = 1;
    } else {
      int32_t _M0L6_2atmpS2722 = _M0Lm16lastRemovedDigitS978;
      _M0L6_2atmpS2721 = _M0L6_2atmpS2722 >= 5;
    }
    #line 487 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2720 = _M0MPC14bool4Bool9to__int64(_M0L6_2atmpS2721);
    _M0L6_2atmpS2719 = *(uint64_t*)&_M0L6_2atmpS2720;
    _M0Lm6outputS979 = _M0L6_2atmpS2718 + _M0L6_2atmpS2719;
  } else {
    int32_t _M0Lm7roundUpS992 = 0;
    uint64_t _M0L6_2atmpS2746 = _M0Lm2vpS954;
    uint64_t _M0L8vpDiv100S993 = _M0L6_2atmpS2746 / 100ull;
    uint64_t _M0L6_2atmpS2745 = _M0Lm2vmS955;
    uint64_t _M0L8vmDiv100S994 = _M0L6_2atmpS2745 / 100ull;
    uint64_t _M0L6_2atmpS2740;
    uint64_t _M0L6_2atmpS2743;
    uint64_t _M0L6_2atmpS2744;
    int32_t _M0L6_2atmpS2742;
    uint64_t _M0L6_2atmpS2741;
    if (_M0L8vpDiv100S993 > _M0L8vmDiv100S994) {
      uint64_t _M0L6_2atmpS2731 = _M0Lm2vrS953;
      uint64_t _M0L8vrDiv100S995 = _M0L6_2atmpS2731 / 100ull;
      uint64_t _M0L6_2atmpS2730 = _M0Lm2vrS953;
      int32_t _M0L6_2atmpS2727 = (int32_t)_M0L6_2atmpS2730;
      int32_t _M0L6_2atmpS2729 = (int32_t)_M0L8vrDiv100S995;
      int32_t _M0L6_2atmpS2728 = 100 * _M0L6_2atmpS2729;
      int32_t _M0L8vrMod100S996 = _M0L6_2atmpS2727 - _M0L6_2atmpS2728;
      int32_t _M0L6_2atmpS2726;
      _M0Lm7roundUpS992 = _M0L8vrMod100S996 >= 50;
      _M0Lm2vrS953 = _M0L8vrDiv100S995;
      _M0Lm2vpS954 = _M0L8vpDiv100S993;
      _M0Lm2vmS955 = _M0L8vmDiv100S994;
      _M0L6_2atmpS2726 = _M0Lm7removedS977;
      _M0Lm7removedS977 = _M0L6_2atmpS2726 + 2;
    }
    while (1) {
      uint64_t _M0L6_2atmpS2739 = _M0Lm2vpS954;
      uint64_t _M0L7vpDiv10S997 = _M0L6_2atmpS2739 / 10ull;
      uint64_t _M0L6_2atmpS2738 = _M0Lm2vmS955;
      uint64_t _M0L7vmDiv10S998 = _M0L6_2atmpS2738 / 10ull;
      uint64_t _M0L6_2atmpS2737;
      uint64_t _M0L7vrDiv10S1000;
      uint64_t _M0L6_2atmpS2736;
      int32_t _M0L6_2atmpS2733;
      int32_t _M0L6_2atmpS2735;
      int32_t _M0L6_2atmpS2734;
      int32_t _M0L7vrMod10S1001;
      int32_t _M0L6_2atmpS2732;
      if (_M0L7vpDiv10S997 <= _M0L7vmDiv10S998) {
        break;
      }
      _M0L6_2atmpS2737 = _M0Lm2vrS953;
      _M0L7vrDiv10S1000 = _M0L6_2atmpS2737 / 10ull;
      _M0L6_2atmpS2736 = _M0Lm2vrS953;
      _M0L6_2atmpS2733 = (int32_t)_M0L6_2atmpS2736;
      _M0L6_2atmpS2735 = (int32_t)_M0L7vrDiv10S1000;
      _M0L6_2atmpS2734 = 10 * _M0L6_2atmpS2735;
      _M0L7vrMod10S1001 = _M0L6_2atmpS2733 - _M0L6_2atmpS2734;
      _M0Lm7roundUpS992 = _M0L7vrMod10S1001 >= 5;
      _M0Lm2vrS953 = _M0L7vrDiv10S1000;
      _M0Lm2vpS954 = _M0L7vpDiv10S997;
      _M0Lm2vmS955 = _M0L7vmDiv10S998;
      _M0L6_2atmpS2732 = _M0Lm7removedS977;
      _M0Lm7removedS977 = _M0L6_2atmpS2732 + 1;
      continue;
      break;
    }
    _M0L6_2atmpS2740 = _M0Lm2vrS953;
    _M0L6_2atmpS2743 = _M0Lm2vrS953;
    _M0L6_2atmpS2744 = _M0Lm2vmS955;
    _M0L6_2atmpS2742
    = _M0L6_2atmpS2743 == _M0L6_2atmpS2744 || _M0Lm7roundUpS992;
    #line 522 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0L6_2atmpS2741 = _M0MPC14bool4Bool10to__uint64(_M0L6_2atmpS2742);
    _M0Lm6outputS979 = _M0L6_2atmpS2740 + _M0L6_2atmpS2741;
  }
  _M0L6_2atmpS2748 = _M0Lm3e10S956;
  _M0L6_2atmpS2749 = _M0Lm7removedS977;
  _M0L3expS1002 = _M0L6_2atmpS2748 + _M0L6_2atmpS2749;
  _M0L6_2atmpS2747 = _M0Lm6outputS979;
  _block_3959
  = (struct _M0TPB17FloatingDecimal64*)moonbit_malloc(sizeof(struct _M0TPB17FloatingDecimal64));
  Moonbit_object_header(_block_3959)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _block_3959->$0 = _M0L6_2atmpS2747;
  _block_3959->$1 = _M0L3expS1002;
  return _block_3959;
}

uint64_t _M0MPC14bool4Bool10to__uint64(int32_t _M0L4selfS945) {
  #line 110 "/home/developer/.moon/lib/core/builtin/bool.mbt"
  if (_M0L4selfS945) {
    return 1ull;
  } else {
    return 0ull;
  }
}

int64_t _M0MPC14bool4Bool9to__int64(int32_t _M0L4selfS944) {
  #line 58 "/home/developer/.moon/lib/core/builtin/bool.mbt"
  if (_M0L4selfS944) {
    return 1ll;
  } else {
    return 0ll;
  }
}

int32_t _M0MPC14bool4Bool7to__int(int32_t _M0L4selfS943) {
  #line 32 "/home/developer/.moon/lib/core/builtin/bool.mbt"
  if (_M0L4selfS943) {
    return 1;
  } else {
    return 0;
  }
}

int32_t _M0FPB17decimal__length17(uint64_t _M0L1vS942) {
  #line 280 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  if (_M0L1vS942 >= 10000000000000000ull) {
    return 17;
  }
  if (_M0L1vS942 >= 1000000000000000ull) {
    return 16;
  }
  if (_M0L1vS942 >= 100000000000000ull) {
    return 15;
  }
  if (_M0L1vS942 >= 10000000000000ull) {
    return 14;
  }
  if (_M0L1vS942 >= 1000000000000ull) {
    return 13;
  }
  if (_M0L1vS942 >= 100000000000ull) {
    return 12;
  }
  if (_M0L1vS942 >= 10000000000ull) {
    return 11;
  }
  if (_M0L1vS942 >= 1000000000ull) {
    return 10;
  }
  if (_M0L1vS942 >= 100000000ull) {
    return 9;
  }
  if (_M0L1vS942 >= 10000000ull) {
    return 8;
  }
  if (_M0L1vS942 >= 1000000ull) {
    return 7;
  }
  if (_M0L1vS942 >= 100000ull) {
    return 6;
  }
  if (_M0L1vS942 >= 10000ull) {
    return 5;
  }
  if (_M0L1vS942 >= 1000ull) {
    return 4;
  }
  if (_M0L1vS942 >= 100ull) {
    return 3;
  }
  if (_M0L1vS942 >= 10ull) {
    return 2;
  }
  return 1;
}

struct _M0TPB8Pow5Pair _M0FPB22double__computeInvPow5(int32_t _M0L1iS925) {
  int32_t _M0L6_2atmpS2648;
  int32_t _M0L6_2atmpS2647;
  int32_t _M0L4baseS924;
  int32_t _M0L5base2S926;
  int32_t _M0L6offsetS927;
  int32_t _M0L6_2atmpS2646;
  uint64_t _M0L4mul0S928;
  int32_t _M0L6_2atmpS2645;
  int32_t _M0L6_2atmpS2644;
  uint64_t _M0L4mul1S929;
  uint64_t _M0L1mS930;
  struct _M0TPB7Umul128 _M0L7_2abindS931;
  uint64_t _M0L7_2alow1S932;
  uint64_t _M0L8_2ahigh1S933;
  struct _M0TPB7Umul128 _M0L7_2abindS934;
  uint64_t _M0L7_2alow0S935;
  uint64_t _M0L8_2ahigh0S936;
  uint64_t _M0L3sumS937;
  uint64_t _M0Lm5high1S938;
  int32_t _M0L6_2atmpS2642;
  int32_t _M0L6_2atmpS2643;
  int32_t _M0L5deltaS939;
  uint64_t _M0L6_2atmpS2641;
  uint64_t _M0L6_2atmpS2633;
  int32_t _M0L6_2atmpS2640;
  uint32_t _M0L6_2atmpS2637;
  int32_t _M0L6_2atmpS2639;
  int32_t _M0L6_2atmpS2638;
  uint32_t _M0L6_2atmpS2636;
  uint32_t _M0L6_2atmpS2635;
  uint64_t _M0L6_2atmpS2634;
  uint64_t _M0L1aS940;
  uint64_t _M0L6_2atmpS2632;
  uint64_t _M0L1bS941;
  #line 239 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2648 = _M0L1iS925 + 26;
  _M0L6_2atmpS2647 = _M0L6_2atmpS2648 - 1;
  _M0L4baseS924 = _M0L6_2atmpS2647 / 26;
  _M0L5base2S926 = _M0L4baseS924 * 26;
  _M0L6offsetS927 = _M0L5base2S926 - _M0L1iS925;
  _M0L6_2atmpS2646 = _M0L4baseS924 * 2;
  #line 243 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L4mul0S928
  = _M0MPC15array13ReadOnlyArray2atGmE(_M0FPB26gDOUBLE__POW5__INV__SPLIT2, _M0L6_2atmpS2646);
  _M0L6_2atmpS2645 = _M0L4baseS924 * 2;
  _M0L6_2atmpS2644 = _M0L6_2atmpS2645 + 1;
  #line 244 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L4mul1S929
  = _M0MPC15array13ReadOnlyArray2atGmE(_M0FPB26gDOUBLE__POW5__INV__SPLIT2, _M0L6_2atmpS2644);
  if (_M0L6offsetS927 == 0) {
    return (struct _M0TPB8Pow5Pair){.$0 = _M0L4mul0S928, .$1 = _M0L4mul1S929};
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L1mS930
  = _M0MPC15array13ReadOnlyArray2atGmE(_M0FPB20gDOUBLE__POW5__TABLE, _M0L6offsetS927);
  #line 249 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS931 = _M0FPB7umul128(_M0L1mS930, _M0L4mul1S929);
  _M0L7_2alow1S932 = _M0L7_2abindS931.$0;
  _M0L8_2ahigh1S933 = _M0L7_2abindS931.$1;
  #line 250 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS934 = _M0FPB7umul128(_M0L1mS930, _M0L4mul0S928);
  _M0L7_2alow0S935 = _M0L7_2abindS934.$0;
  _M0L8_2ahigh0S936 = _M0L7_2abindS934.$1;
  _M0L3sumS937 = _M0L8_2ahigh0S936 + _M0L7_2alow1S932;
  _M0Lm5high1S938 = _M0L8_2ahigh1S933;
  if (_M0L3sumS937 < _M0L8_2ahigh0S936) {
    uint64_t _M0L6_2atmpS2631 = _M0Lm5high1S938;
    _M0Lm5high1S938 = _M0L6_2atmpS2631 + 1ull;
  }
  #line 256 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2642 = _M0FPB8pow5bits(_M0L5base2S926);
  #line 256 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2643 = _M0FPB8pow5bits(_M0L1iS925);
  _M0L5deltaS939 = _M0L6_2atmpS2642 - _M0L6_2atmpS2643;
  #line 257 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2641
  = _M0FPB13shiftright128(_M0L7_2alow0S935, _M0L3sumS937, _M0L5deltaS939);
  _M0L6_2atmpS2633 = _M0L6_2atmpS2641 + 1ull;
  _M0L6_2atmpS2640 = _M0L1iS925 / 16;
  #line 259 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2637
  = _M0MPC15array13ReadOnlyArray2atGjE(_M0FPB19gPOW5__INV__OFFSETS, _M0L6_2atmpS2640);
  _M0L6_2atmpS2639 = _M0L1iS925 % 16;
  _M0L6_2atmpS2638 = _M0L6_2atmpS2639 << 1;
  _M0L6_2atmpS2636 = _M0L6_2atmpS2637 >> (_M0L6_2atmpS2638 & 31);
  _M0L6_2atmpS2635 = _M0L6_2atmpS2636 & 3u;
  #line 259 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2634 = _M0MPC14uint4UInt10to__uint64(_M0L6_2atmpS2635);
  _M0L1aS940 = _M0L6_2atmpS2633 + _M0L6_2atmpS2634;
  _M0L6_2atmpS2632 = _M0Lm5high1S938;
  #line 260 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L1bS941
  = _M0FPB13shiftright128(_M0L3sumS937, _M0L6_2atmpS2632, _M0L5deltaS939);
  return (struct _M0TPB8Pow5Pair){.$0 = _M0L1aS940, .$1 = _M0L1bS941};
}

struct _M0TPB8Pow5Pair _M0FPB19double__computePow5(int32_t _M0L1iS907) {
  int32_t _M0L4baseS906;
  int32_t _M0L5base2S908;
  int32_t _M0L6offsetS909;
  int32_t _M0L6_2atmpS2630;
  uint64_t _M0L4mul0S910;
  int32_t _M0L6_2atmpS2629;
  int32_t _M0L6_2atmpS2628;
  uint64_t _M0L4mul1S911;
  uint64_t _M0L1mS912;
  struct _M0TPB7Umul128 _M0L7_2abindS913;
  uint64_t _M0L7_2alow1S914;
  uint64_t _M0L8_2ahigh1S915;
  struct _M0TPB7Umul128 _M0L7_2abindS916;
  uint64_t _M0L7_2alow0S917;
  uint64_t _M0L8_2ahigh0S918;
  uint64_t _M0L3sumS919;
  uint64_t _M0Lm5high1S920;
  int32_t _M0L6_2atmpS2626;
  int32_t _M0L6_2atmpS2627;
  int32_t _M0L5deltaS921;
  uint64_t _M0L6_2atmpS2618;
  int32_t _M0L6_2atmpS2625;
  uint32_t _M0L6_2atmpS2622;
  int32_t _M0L6_2atmpS2624;
  int32_t _M0L6_2atmpS2623;
  uint32_t _M0L6_2atmpS2621;
  uint32_t _M0L6_2atmpS2620;
  uint64_t _M0L6_2atmpS2619;
  uint64_t _M0L1aS922;
  uint64_t _M0L6_2atmpS2617;
  uint64_t _M0L1bS923;
  #line 213 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L4baseS906 = _M0L1iS907 / 26;
  _M0L5base2S908 = _M0L4baseS906 * 26;
  _M0L6offsetS909 = _M0L1iS907 - _M0L5base2S908;
  _M0L6_2atmpS2630 = _M0L4baseS906 * 2;
  #line 217 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L4mul0S910
  = _M0MPC15array13ReadOnlyArray2atGmE(_M0FPB21gDOUBLE__POW5__SPLIT2, _M0L6_2atmpS2630);
  _M0L6_2atmpS2629 = _M0L4baseS906 * 2;
  _M0L6_2atmpS2628 = _M0L6_2atmpS2629 + 1;
  #line 218 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L4mul1S911
  = _M0MPC15array13ReadOnlyArray2atGmE(_M0FPB21gDOUBLE__POW5__SPLIT2, _M0L6_2atmpS2628);
  if (_M0L6offsetS909 == 0) {
    return (struct _M0TPB8Pow5Pair){.$0 = _M0L4mul0S910, .$1 = _M0L4mul1S911};
  }
  #line 222 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L1mS912
  = _M0MPC15array13ReadOnlyArray2atGmE(_M0FPB20gDOUBLE__POW5__TABLE, _M0L6offsetS909);
  #line 223 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS913 = _M0FPB7umul128(_M0L1mS912, _M0L4mul1S911);
  _M0L7_2alow1S914 = _M0L7_2abindS913.$0;
  _M0L8_2ahigh1S915 = _M0L7_2abindS913.$1;
  #line 224 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS916 = _M0FPB7umul128(_M0L1mS912, _M0L4mul0S910);
  _M0L7_2alow0S917 = _M0L7_2abindS916.$0;
  _M0L8_2ahigh0S918 = _M0L7_2abindS916.$1;
  _M0L3sumS919 = _M0L8_2ahigh0S918 + _M0L7_2alow1S914;
  _M0Lm5high1S920 = _M0L8_2ahigh1S915;
  if (_M0L3sumS919 < _M0L8_2ahigh0S918) {
    uint64_t _M0L6_2atmpS2616 = _M0Lm5high1S920;
    _M0Lm5high1S920 = _M0L6_2atmpS2616 + 1ull;
  }
  #line 230 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2626 = _M0FPB8pow5bits(_M0L1iS907);
  #line 230 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2627 = _M0FPB8pow5bits(_M0L5base2S908);
  _M0L5deltaS921 = _M0L6_2atmpS2626 - _M0L6_2atmpS2627;
  #line 231 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2618
  = _M0FPB13shiftright128(_M0L7_2alow0S917, _M0L3sumS919, _M0L5deltaS921);
  _M0L6_2atmpS2625 = _M0L1iS907 / 16;
  #line 232 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2622
  = _M0MPC15array13ReadOnlyArray2atGjE(_M0FPB14gPOW5__OFFSETS, _M0L6_2atmpS2625);
  _M0L6_2atmpS2624 = _M0L1iS907 % 16;
  _M0L6_2atmpS2623 = _M0L6_2atmpS2624 << 1;
  _M0L6_2atmpS2621 = _M0L6_2atmpS2622 >> (_M0L6_2atmpS2623 & 31);
  _M0L6_2atmpS2620 = _M0L6_2atmpS2621 & 3u;
  #line 232 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2619 = _M0MPC14uint4UInt10to__uint64(_M0L6_2atmpS2620);
  _M0L1aS922 = _M0L6_2atmpS2618 + _M0L6_2atmpS2619;
  _M0L6_2atmpS2617 = _M0Lm5high1S920;
  #line 233 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L1bS923
  = _M0FPB13shiftright128(_M0L3sumS919, _M0L6_2atmpS2617, _M0L5deltaS921);
  return (struct _M0TPB8Pow5Pair){.$0 = _M0L1aS922, .$1 = _M0L1bS923};
}

struct _M0TPB19MulShiftAll64Result _M0FPB13mulShiftAll64(
  uint64_t _M0L1mS880,
  struct _M0TPB8Pow5Pair _M0L3mulS877,
  int32_t _M0L1jS893,
  int32_t _M0L7mmShiftS895
) {
  uint64_t _M0L7_2amul0S876;
  uint64_t _M0L7_2amul1S878;
  uint64_t _M0L1mS879;
  struct _M0TPB7Umul128 _M0L7_2abindS881;
  uint64_t _M0L5_2aloS882;
  uint64_t _M0L6_2atmpS883;
  struct _M0TPB7Umul128 _M0L7_2abindS884;
  uint64_t _M0L6_2alo2S885;
  uint64_t _M0L6_2ahi2S886;
  uint64_t _M0L3midS887;
  uint64_t _M0L6_2atmpS2615;
  uint64_t _M0L2hiS888;
  uint64_t _M0L3lo2S889;
  uint64_t _M0L6_2atmpS2613;
  uint64_t _M0L6_2atmpS2614;
  uint64_t _M0L4mid2S890;
  uint64_t _M0L6_2atmpS2612;
  uint64_t _M0L3hi2S891;
  int32_t _M0L6_2atmpS2611;
  int32_t _M0L6_2atmpS2610;
  uint64_t _M0L2vpS892;
  uint64_t _M0Lm2vmS894;
  int32_t _M0L6_2atmpS2609;
  int32_t _M0L6_2atmpS2608;
  uint64_t _M0L2vrS905;
  uint64_t _M0L6_2atmpS2607;
  #line 129 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2amul0S876 = _M0L3mulS877.$0;
  _M0L7_2amul1S878 = _M0L3mulS877.$1;
  _M0L1mS879 = _M0L1mS880 << 1;
  #line 137 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS881 = _M0FPB7umul128(_M0L1mS879, _M0L7_2amul0S876);
  _M0L5_2aloS882 = _M0L7_2abindS881.$0;
  _M0L6_2atmpS883 = _M0L7_2abindS881.$1;
  #line 138 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L7_2abindS884 = _M0FPB7umul128(_M0L1mS879, _M0L7_2amul1S878);
  _M0L6_2alo2S885 = _M0L7_2abindS884.$0;
  _M0L6_2ahi2S886 = _M0L7_2abindS884.$1;
  _M0L3midS887 = _M0L6_2atmpS883 + _M0L6_2alo2S885;
  if (_M0L3midS887 < _M0L6_2atmpS883) {
    _M0L6_2atmpS2615 = 1ull;
  } else {
    _M0L6_2atmpS2615 = 0ull;
  }
  _M0L2hiS888 = _M0L6_2ahi2S886 + _M0L6_2atmpS2615;
  _M0L3lo2S889 = _M0L5_2aloS882 + _M0L7_2amul0S876;
  _M0L6_2atmpS2613 = _M0L3midS887 + _M0L7_2amul1S878;
  if (_M0L3lo2S889 < _M0L5_2aloS882) {
    _M0L6_2atmpS2614 = 1ull;
  } else {
    _M0L6_2atmpS2614 = 0ull;
  }
  _M0L4mid2S890 = _M0L6_2atmpS2613 + _M0L6_2atmpS2614;
  if (_M0L4mid2S890 < _M0L3midS887) {
    _M0L6_2atmpS2612 = 1ull;
  } else {
    _M0L6_2atmpS2612 = 0ull;
  }
  _M0L3hi2S891 = _M0L2hiS888 + _M0L6_2atmpS2612;
  _M0L6_2atmpS2611 = _M0L1jS893 - 64;
  _M0L6_2atmpS2610 = _M0L6_2atmpS2611 - 1;
  #line 144 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L2vpS892
  = _M0FPB13shiftright128(_M0L4mid2S890, _M0L3hi2S891, _M0L6_2atmpS2610);
  _M0Lm2vmS894 = 0ull;
  if (_M0L7mmShiftS895) {
    uint64_t _M0L3lo3S896 = _M0L5_2aloS882 - _M0L7_2amul0S876;
    uint64_t _M0L6_2atmpS2597 = _M0L3midS887 - _M0L7_2amul1S878;
    uint64_t _M0L6_2atmpS2598;
    uint64_t _M0L4mid3S897;
    uint64_t _M0L6_2atmpS2596;
    uint64_t _M0L3hi3S898;
    int32_t _M0L6_2atmpS2595;
    int32_t _M0L6_2atmpS2594;
    if (_M0L5_2aloS882 < _M0L3lo3S896) {
      _M0L6_2atmpS2598 = 1ull;
    } else {
      _M0L6_2atmpS2598 = 0ull;
    }
    _M0L4mid3S897 = _M0L6_2atmpS2597 - _M0L6_2atmpS2598;
    if (_M0L3midS887 < _M0L4mid3S897) {
      _M0L6_2atmpS2596 = 1ull;
    } else {
      _M0L6_2atmpS2596 = 0ull;
    }
    _M0L3hi3S898 = _M0L2hiS888 - _M0L6_2atmpS2596;
    _M0L6_2atmpS2595 = _M0L1jS893 - 64;
    _M0L6_2atmpS2594 = _M0L6_2atmpS2595 - 1;
    #line 150 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0Lm2vmS894
    = _M0FPB13shiftright128(_M0L4mid3S897, _M0L3hi3S898, _M0L6_2atmpS2594);
  } else {
    uint64_t _M0L3lo3S899 = _M0L5_2aloS882 + _M0L5_2aloS882;
    uint64_t _M0L6_2atmpS2605 = _M0L3midS887 + _M0L3midS887;
    uint64_t _M0L6_2atmpS2606;
    uint64_t _M0L4mid3S900;
    uint64_t _M0L6_2atmpS2603;
    uint64_t _M0L6_2atmpS2604;
    uint64_t _M0L3hi3S901;
    uint64_t _M0L3lo4S902;
    uint64_t _M0L6_2atmpS2601;
    uint64_t _M0L6_2atmpS2602;
    uint64_t _M0L4mid4S903;
    uint64_t _M0L6_2atmpS2600;
    uint64_t _M0L3hi4S904;
    int32_t _M0L6_2atmpS2599;
    if (_M0L3lo3S899 < _M0L5_2aloS882) {
      _M0L6_2atmpS2606 = 1ull;
    } else {
      _M0L6_2atmpS2606 = 0ull;
    }
    _M0L4mid3S900 = _M0L6_2atmpS2605 + _M0L6_2atmpS2606;
    _M0L6_2atmpS2603 = _M0L2hiS888 + _M0L2hiS888;
    if (_M0L4mid3S900 < _M0L3midS887) {
      _M0L6_2atmpS2604 = 1ull;
    } else {
      _M0L6_2atmpS2604 = 0ull;
    }
    _M0L3hi3S901 = _M0L6_2atmpS2603 + _M0L6_2atmpS2604;
    _M0L3lo4S902 = _M0L3lo3S899 - _M0L7_2amul0S876;
    _M0L6_2atmpS2601 = _M0L4mid3S900 - _M0L7_2amul1S878;
    if (_M0L3lo3S899 < _M0L3lo4S902) {
      _M0L6_2atmpS2602 = 1ull;
    } else {
      _M0L6_2atmpS2602 = 0ull;
    }
    _M0L4mid4S903 = _M0L6_2atmpS2601 - _M0L6_2atmpS2602;
    if (_M0L4mid3S900 < _M0L4mid4S903) {
      _M0L6_2atmpS2600 = 1ull;
    } else {
      _M0L6_2atmpS2600 = 0ull;
    }
    _M0L3hi4S904 = _M0L3hi3S901 - _M0L6_2atmpS2600;
    _M0L6_2atmpS2599 = _M0L1jS893 - 64;
    #line 158 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _M0Lm2vmS894
    = _M0FPB13shiftright128(_M0L4mid4S903, _M0L3hi4S904, _M0L6_2atmpS2599);
  }
  _M0L6_2atmpS2609 = _M0L1jS893 - 64;
  _M0L6_2atmpS2608 = _M0L6_2atmpS2609 - 1;
  #line 160 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L2vrS905
  = _M0FPB13shiftright128(_M0L3midS887, _M0L2hiS888, _M0L6_2atmpS2608);
  _M0L6_2atmpS2607 = _M0Lm2vmS894;
  return (struct _M0TPB19MulShiftAll64Result){.$0 = _M0L2vrS905,
                                                .$1 = _M0L2vpS892,
                                                .$2 = _M0L6_2atmpS2607};
}

int32_t _M0FPB18multipleOfPowerOf2(
  uint64_t _M0L5valueS874,
  int32_t _M0L1pS875
) {
  uint64_t _M0L6_2atmpS2593;
  uint64_t _M0L6_2atmpS2592;
  uint64_t _M0L6_2atmpS2591;
  #line 124 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2593 = 1ull << (_M0L1pS875 & 63);
  _M0L6_2atmpS2592 = _M0L6_2atmpS2593 - 1ull;
  _M0L6_2atmpS2591 = _M0L5valueS874 & _M0L6_2atmpS2592;
  return _M0L6_2atmpS2591 == 0ull;
}

int32_t _M0FPB18multipleOfPowerOf5(
  uint64_t _M0L5valueS872,
  int32_t _M0L1pS873
) {
  int32_t _M0L6_2atmpS2590;
  #line 119 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  #line 120 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2590 = _M0FPB10pow5Factor(_M0L5valueS872);
  return _M0L6_2atmpS2590 >= _M0L1pS873;
}

int32_t _M0FPB10pow5Factor(uint64_t _M0L5valueS867) {
  uint64_t _M0L6_2atmpS2581;
  uint64_t _M0L6_2atmpS2582;
  uint64_t _M0L6_2atmpS2583;
  uint64_t _M0L6_2atmpS2584;
  uint64_t _M0L6_2atmpS2589;
  int32_t _M0L5countS868;
  uint64_t _M0L1vS869;
  #line 94 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2581 = _M0L5valueS867 % 5ull;
  if (_M0L6_2atmpS2581 != 0ull) {
    return 0;
  }
  _M0L6_2atmpS2582 = _M0L5valueS867 % 25ull;
  if (_M0L6_2atmpS2582 != 0ull) {
    return 1;
  }
  _M0L6_2atmpS2583 = _M0L5valueS867 % 125ull;
  if (_M0L6_2atmpS2583 != 0ull) {
    return 2;
  }
  _M0L6_2atmpS2584 = _M0L5valueS867 % 625ull;
  if (_M0L6_2atmpS2584 != 0ull) {
    return 3;
  }
  _M0L6_2atmpS2589 = _M0L5valueS867 / 625ull;
  _M0L5countS868 = 4;
  _M0L1vS869 = _M0L6_2atmpS2589;
  while (1) {
    if (_M0L1vS869 > 0ull) {
      uint64_t _M0L6_2atmpS2585 = _M0L1vS869 % 5ull;
      int32_t _M0L6_2atmpS2586;
      uint64_t _M0L6_2atmpS2587;
      if (_M0L6_2atmpS2585 != 0ull) {
        return _M0L5countS868;
      }
      _M0L6_2atmpS2586 = _M0L5countS868 + 1;
      _M0L6_2atmpS2587 = _M0L1vS869 / 5ull;
      _M0L5countS868 = _M0L6_2atmpS2586;
      _M0L1vS869 = _M0L6_2atmpS2587;
      continue;
    } else {
      struct _M0TPB13StringBuilder* _M0L18_2astring__builderS871;
      moonbit_string_t _M0L6_2atmpS2588;
      int32_t _result_3961;
      #line 114 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      _M0L18_2astring__builderS871
      = _M0MPB13StringBuilder21StringBuilder_2einner(25);
      #line 114 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS871, (moonbit_string_t)moonbit_string_literal_47.data);
      #line 114 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      _M0MPB13StringBuilder13write__objectGmE(_M0L18_2astring__builderS871, _M0L5valueS867);
      #line 114 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      _M0L6_2atmpS2588
      = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS871);
      moonbit_decref(_M0L18_2astring__builderS871);
      #line 114 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
      _result_3961 = _M0FPC15abort5abortGiE(_M0L6_2atmpS2588);
      moonbit_decref(_M0L6_2atmpS2588);
      return _result_3961;
    }
    break;
  }
}

uint64_t _M0FPB13shiftright128(
  uint64_t _M0L2loS866,
  uint64_t _M0L2hiS864,
  int32_t _M0L4distS865
) {
  int32_t _M0L6_2atmpS2580;
  uint64_t _M0L6_2atmpS2578;
  uint64_t _M0L6_2atmpS2579;
  #line 89 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2580 = 64 - _M0L4distS865;
  _M0L6_2atmpS2578 = _M0L2hiS864 << (_M0L6_2atmpS2580 & 63);
  _M0L6_2atmpS2579 = _M0L2loS866 >> (_M0L4distS865 & 63);
  return _M0L6_2atmpS2578 | _M0L6_2atmpS2579;
}

struct _M0TPB7Umul128 _M0FPB7umul128(
  uint64_t _M0L1aS854,
  uint64_t _M0L1bS857
) {
  uint64_t _M0L3aLoS853;
  uint64_t _M0L3aHiS855;
  uint64_t _M0L3bLoS856;
  uint64_t _M0L3bHiS858;
  uint64_t _M0L1xS859;
  uint64_t _M0L6_2atmpS2576;
  uint64_t _M0L6_2atmpS2577;
  uint64_t _M0L1yS860;
  uint64_t _M0L6_2atmpS2574;
  uint64_t _M0L6_2atmpS2575;
  uint64_t _M0L1zS861;
  uint64_t _M0L6_2atmpS2572;
  uint64_t _M0L6_2atmpS2573;
  uint64_t _M0L6_2atmpS2570;
  uint64_t _M0L6_2atmpS2571;
  uint64_t _M0L1wS862;
  uint64_t _M0L2loS863;
  #line 74 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L3aLoS853 = _M0L1aS854 & 4294967295ull;
  _M0L3aHiS855 = _M0L1aS854 >> 32;
  _M0L3bLoS856 = _M0L1bS857 & 4294967295ull;
  _M0L3bHiS858 = _M0L1bS857 >> 32;
  _M0L1xS859 = _M0L3aLoS853 * _M0L3bLoS856;
  _M0L6_2atmpS2576 = _M0L3aHiS855 * _M0L3bLoS856;
  _M0L6_2atmpS2577 = _M0L1xS859 >> 32;
  _M0L1yS860 = _M0L6_2atmpS2576 + _M0L6_2atmpS2577;
  _M0L6_2atmpS2574 = _M0L3aLoS853 * _M0L3bHiS858;
  _M0L6_2atmpS2575 = _M0L1yS860 & 4294967295ull;
  _M0L1zS861 = _M0L6_2atmpS2574 + _M0L6_2atmpS2575;
  _M0L6_2atmpS2572 = _M0L3aHiS855 * _M0L3bHiS858;
  _M0L6_2atmpS2573 = _M0L1yS860 >> 32;
  _M0L6_2atmpS2570 = _M0L6_2atmpS2572 + _M0L6_2atmpS2573;
  _M0L6_2atmpS2571 = _M0L1zS861 >> 32;
  _M0L1wS862 = _M0L6_2atmpS2570 + _M0L6_2atmpS2571;
  _M0L2loS863 = _M0L1aS854 * _M0L1bS857;
  return (struct _M0TPB7Umul128){.$0 = _M0L2loS863, .$1 = _M0L1wS862};
}

moonbit_string_t _M0FPB19string__from__bytes(
  moonbit_bytes_t _M0L5bytesS851,
  int32_t _M0L4fromS848,
  int32_t _M0L2toS847
) {
  int32_t _M0L3lenS846;
  int32_t _M0L6_2atmpS2569;
  uint16_t* _M0L6bufferS849;
  int32_t _M0L1iS850;
  #line 52 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L3lenS846 = _M0L2toS847 - _M0L4fromS848;
  #line 54 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2569 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L6bufferS849
  = (uint16_t*)moonbit_make_string(_M0L3lenS846, _M0L6_2atmpS2569);
  _M0L1iS850 = 0;
  while (1) {
    if (_M0L1iS850 < _M0L3lenS846) {
      int32_t _M0L6_2atmpS2567 = _M0L4fromS848 + _M0L1iS850;
      int32_t _M0L6_2atmpS2566;
      int32_t _M0L6_2atmpS2565;
      int32_t _M0L6_2atmpS2568;
      if (
        _M0L6_2atmpS2567 < 0
        || _M0L6_2atmpS2567 >= Moonbit_array_length(_M0L5bytesS851)
      ) {
        #line 56 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2566 = (int32_t)_M0L5bytesS851[_M0L6_2atmpS2567];
      _M0L6_2atmpS2565 = (uint16_t)_M0L6_2atmpS2566;
      if (
        _M0L1iS850 < 0 || _M0L1iS850 >= Moonbit_array_length(_M0L6bufferS849)
      ) {
        #line 56 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
        moonbit_panic();
      }
      _M0L6bufferS849[_M0L1iS850] = _M0L6_2atmpS2565;
      _M0L6_2atmpS2568 = _M0L1iS850 + 1;
      _M0L1iS850 = _M0L6_2atmpS2568;
      continue;
    }
    break;
  }
  return _M0L6bufferS849;
}

int32_t _M0FPB9log10Pow2(int32_t _M0L1eS845) {
  int32_t _M0L6_2atmpS2564;
  uint32_t _M0L6_2atmpS2563;
  uint32_t _M0L6_2atmpS2562;
  #line 44 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2564 = _M0L1eS845 * 78913;
  _M0L6_2atmpS2563 = *(uint32_t*)&_M0L6_2atmpS2564;
  _M0L6_2atmpS2562 = _M0L6_2atmpS2563 >> 18;
  return *(int32_t*)&_M0L6_2atmpS2562;
}

int32_t _M0FPB9log10Pow5(int32_t _M0L1eS844) {
  int32_t _M0L6_2atmpS2561;
  uint32_t _M0L6_2atmpS2560;
  uint32_t _M0L6_2atmpS2559;
  #line 37 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2561 = _M0L1eS844 * 732923;
  _M0L6_2atmpS2560 = *(uint32_t*)&_M0L6_2atmpS2561;
  _M0L6_2atmpS2559 = _M0L6_2atmpS2560 >> 20;
  return *(int32_t*)&_M0L6_2atmpS2559;
}

moonbit_string_t _M0FPB18copy__special__str(
  int32_t _M0L4signS842,
  int32_t _M0L8exponentS843,
  int32_t _M0L8mantissaS840
) {
  moonbit_string_t _M0L1sS841;
  moonbit_string_t _result_3964;
  #line 23 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  if (_M0L8mantissaS840) {
    return (moonbit_string_t)moonbit_string_literal_48.data;
  }
  if (_M0L4signS842) {
    _M0L1sS841 = (moonbit_string_t)moonbit_string_literal_49.data;
  } else {
    _M0L1sS841 = (moonbit_string_t)moonbit_string_literal_1.data;
  }
  if (_M0L8exponentS843) {
    moonbit_string_t _result_3963;
    #line 29 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
    _result_3963
    = moonbit_add_string(_M0L1sS841, (moonbit_string_t)moonbit_string_literal_50.data);
    moonbit_decref(_M0L1sS841);
    return _result_3963;
  }
  #line 31 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _result_3964
  = moonbit_add_string(_M0L1sS841, (moonbit_string_t)moonbit_string_literal_51.data);
  moonbit_decref(_M0L1sS841);
  return _result_3964;
}

int32_t _M0FPB8pow5bits(int32_t _M0L1eS839) {
  int32_t _M0L6_2atmpS2558;
  uint32_t _M0L6_2atmpS2557;
  uint32_t _M0L6_2atmpS2556;
  int32_t _M0L6_2atmpS2555;
  #line 18 "/home/developer/.moon/lib/core/builtin/double_ryu_nonjs.mbt"
  _M0L6_2atmpS2558 = _M0L1eS839 * 1217359;
  _M0L6_2atmpS2557 = *(uint32_t*)&_M0L6_2atmpS2558;
  _M0L6_2atmpS2556 = _M0L6_2atmpS2557 >> 19;
  _M0L6_2atmpS2555 = *(int32_t*)&_M0L6_2atmpS2556;
  return _M0L6_2atmpS2555 + 1;
}

int32_t _M0IPC13int3IntPB4Hash4hash(int32_t _M0L4selfS838) {
  int32_t _tmp_3965;
  uint32_t _M0L6_2atmpS2554;
  uint32_t _M0L6_2atmpS2553;
  uint32_t _M0L3accS837;
  uint32_t _M0L6_2atmpS2552;
  uint32_t _M0L6_2atmpS2551;
  #line 588 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_3965 = 0;
  _M0L6_2atmpS2554 = *(uint32_t*)&_tmp_3965;
  _M0L6_2atmpS2553 = _M0L6_2atmpS2554 + 374761393u;
  _M0L3accS837 = _M0L6_2atmpS2553 + 4u;
  _M0L6_2atmpS2552 = *(uint32_t*)&_M0L4selfS838;
  #line 590 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS2551 = _M0FPB13consume4__acc(_M0L3accS837, _M0L6_2atmpS2552);
  #line 590 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS2551);
}

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t _M0L4selfS833) {
  int32_t _tmp_3966;
  uint32_t _M0L6_2atmpS2550;
  uint32_t _M0Lm3accS831;
  int32_t _M0L7_2abindS832;
  int32_t _M0L1iS834;
  uint32_t _M0L6_2atmpS2549;
  #line 522 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_3966 = 0;
  _M0L6_2atmpS2550 = *(uint32_t*)&_tmp_3966;
  _M0Lm3accS831 = _M0L6_2atmpS2550 + 374761393u;
  _M0L7_2abindS832 = Moonbit_array_length(_M0L4selfS833);
  _M0L1iS834 = 0;
  while (1) {
    if (_M0L1iS834 < _M0L7_2abindS832) {
      uint32_t _M0L6_2atmpS2544 = _M0Lm3accS831;
      int32_t _M0L6_2atmpS2547;
      int32_t _M0L6_2atmpS2546;
      uint32_t _M0L1vS835;
      uint32_t _M0L6_2atmpS2545;
      int32_t _M0L6_2atmpS2548;
      _M0Lm3accS831 = _M0L6_2atmpS2544 + 4u;
      _M0L6_2atmpS2547 = _M0L4selfS833[_M0L1iS834];
      _M0L6_2atmpS2546 = (int32_t)_M0L6_2atmpS2547;
      _M0L1vS835 = *(uint32_t*)&_M0L6_2atmpS2546;
      _M0L6_2atmpS2545 = _M0Lm3accS831;
      #line 527 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
      _M0Lm3accS831 = _M0FPB13consume4__acc(_M0L6_2atmpS2545, _M0L1vS835);
      _M0L6_2atmpS2548 = _M0L1iS834 + 1;
      _M0L1iS834 = _M0L6_2atmpS2548;
      continue;
    }
    break;
  }
  _M0L6_2atmpS2549 = _M0Lm3accS831;
  #line 529 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS2549);
}

struct _M0TUWEuQRPC15error5ErrorNsE* _M0MPB3Map3getGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS818,
  int32_t _M0L3keyS814
) {
  int32_t _M0L4hashS813;
  int32_t _M0L14capacity__maskS2529;
  int32_t _M0L6_2atmpS2528;
  int32_t _M0L1iS815;
  int32_t _M0L3idxS816;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS813 = _M0IPC13int3IntPB4Hash4hash(_M0L3keyS814);
  _M0L14capacity__maskS2529 = _M0L4selfS818->$3;
  _M0L6_2atmpS2528 = _M0L4hashS813 & _M0L14capacity__maskS2529;
  _M0L1iS815 = 0;
  _M0L3idxS816 = _M0L6_2atmpS2528;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2527 =
      _M0L4selfS818->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS817;
    if (
      _M0L3idxS816 < 0
      || _M0L3idxS816 >= Moonbit_array_length(_M0L7entriesS2527)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS817
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2527[
        _M0L3idxS816
      ];
    if (_M0L7_2abindS817 == 0) {
      struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2516 = 0;
      return _M0L6_2atmpS2516;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS819 =
        _M0L7_2abindS817;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L8_2aentryS820 =
        _M0L7_2aSomeS819;
      int32_t _M0L4hashS2518 = _M0L8_2aentryS820->$3;
      int32_t _if__result_3969;
      int32_t _M0L3pslS2521;
      int32_t _M0L6_2atmpS2523;
      int32_t _M0L6_2atmpS2525;
      int32_t _M0L14capacity__maskS2526;
      int32_t _M0L6_2atmpS2524;
      if (_M0L4hashS2518 == _M0L4hashS813) {
        int32_t _M0L3keyS2517 = _M0L8_2aentryS820->$4;
        _if__result_3969 = _M0L3keyS2517 == _M0L3keyS814;
      } else {
        _if__result_3969 = 0;
      }
      if (_if__result_3969) {
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L5valueS2520 =
          _M0L8_2aentryS820->$5;
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2519;
        moonbit_incref(_M0L5valueS2520);
        _M0L6_2atmpS2519 = _M0L5valueS2520;
        return _M0L6_2atmpS2519;
      } else {
        moonbit_incref(_M0L8_2aentryS820);
      }
      _M0L3pslS2521 = _M0L8_2aentryS820->$2;
      moonbit_decref(_M0L8_2aentryS820);
      if (_M0L1iS815 > _M0L3pslS2521) {
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2522 = 0;
        return _M0L6_2atmpS2522;
      }
      _M0L6_2atmpS2523 = _M0L1iS815 + 1;
      _M0L6_2atmpS2525 = _M0L3idxS816 + 1;
      _M0L14capacity__maskS2526 = _M0L4selfS818->$3;
      _M0L6_2atmpS2524 = _M0L6_2atmpS2525 & _M0L14capacity__maskS2526;
      _M0L1iS815 = _M0L6_2atmpS2523;
      _M0L3idxS816 = _M0L6_2atmpS2524;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3getGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS827,
  moonbit_string_t _M0L3keyS823
) {
  int32_t _M0L4hashS822;
  int32_t _M0L14capacity__maskS2543;
  int32_t _M0L6_2atmpS2542;
  int32_t _M0L1iS824;
  int32_t _M0L3idxS825;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS822 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS823);
  _M0L14capacity__maskS2543 = _M0L4selfS827->$3;
  _M0L6_2atmpS2542 = _M0L4hashS822 & _M0L14capacity__maskS2543;
  _M0L1iS824 = 0;
  _M0L3idxS825 = _M0L6_2atmpS2542;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2541 =
      _M0L4selfS827->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS826;
    if (
      _M0L3idxS825 < 0
      || _M0L3idxS825 >= Moonbit_array_length(_M0L7entriesS2541)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS826
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2541[
        _M0L3idxS825
      ];
    if (_M0L7_2abindS826 == 0) {
      struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2530 = 0;
      return _M0L6_2atmpS2530;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS828 =
        _M0L7_2abindS826;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2aentryS829 =
        _M0L7_2aSomeS828;
      int32_t _M0L4hashS2532 = _M0L8_2aentryS829->$3;
      int32_t _if__result_3971;
      int32_t _M0L3pslS2535;
      int32_t _M0L6_2atmpS2537;
      int32_t _M0L6_2atmpS2539;
      int32_t _M0L14capacity__maskS2540;
      int32_t _M0L6_2atmpS2538;
      if (_M0L4hashS2532 == _M0L4hashS822) {
        moonbit_string_t _M0L3keyS2531 = _M0L8_2aentryS829->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3971
        = _M0L3keyS2531 == _M0L3keyS823
          || Moonbit_array_length(_M0L3keyS2531)
             == Moonbit_array_length(_M0L3keyS823)
             && 0
                == memcmp(_M0L3keyS2531, _M0L3keyS823, Moonbit_array_length(_M0L3keyS2531) * 2);
      } else {
        _if__result_3971 = 0;
      }
      if (_if__result_3971) {
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L5valueS2534 =
          _M0L8_2aentryS829->$5;
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2533;
        moonbit_incref(_M0L5valueS2534);
        _M0L6_2atmpS2533 = _M0L5valueS2534;
        return _M0L6_2atmpS2533;
      } else {
        moonbit_incref(_M0L8_2aentryS829);
      }
      _M0L3pslS2535 = _M0L8_2aentryS829->$2;
      moonbit_decref(_M0L8_2aentryS829);
      if (_M0L1iS824 > _M0L3pslS2535) {
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2536 = 0;
        return _M0L6_2atmpS2536;
      }
      _M0L6_2atmpS2537 = _M0L1iS824 + 1;
      _M0L6_2atmpS2539 = _M0L3idxS825 + 1;
      _M0L14capacity__maskS2540 = _M0L4selfS827->$3;
      _M0L6_2atmpS2538 = _M0L6_2atmpS2539 & _M0L14capacity__maskS2540;
      _M0L1iS824 = _M0L6_2atmpS2537;
      _M0L3idxS825 = _M0L6_2atmpS2538;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPB3Map3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE _M0L3arrS792,
  int64_t _M0L8capacityS794
) {
  int32_t _M0L3endS2503;
  int32_t _M0L5startS2504;
  int32_t _M0L6lengthS791;
  int32_t _M0L8capacityS793;
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L1mS797;
  int32_t _M0L3endS2500;
  int32_t _M0L5startS2501;
  int32_t _M0L7_2abindS798;
  int32_t _M0L2__S799;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2503 = _M0L3arrS792.$2;
  _M0L5startS2504 = _M0L3arrS792.$1;
  _M0L6lengthS791 = _M0L3endS2503 - _M0L5startS2504;
  if (_M0L8capacityS794 == 4294967296ll) {
    if (_M0L6lengthS791 == 0) {
      _M0L8capacityS793 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS793 = _M0FPB21capacity__for__length(_M0L6lengthS791);
    }
  } else {
    int64_t _M0L7_2aSomeS795 = _M0L8capacityS794;
    int32_t _M0L11_2acapacityS796 = (int32_t)_M0L7_2aSomeS795;
    int32_t _M0L6_2atmpS2502;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2502 = _M0FPB21capacity__for__length(_M0L6lengthS791);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS793
    = _M0MPC13int3Int3max(_M0L11_2acapacityS796, _M0L6_2atmpS2502);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS797
  = _M0FPB8new__mapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L8capacityS793);
  _M0L3endS2500 = _M0L3arrS792.$2;
  _M0L5startS2501 = _M0L3arrS792.$1;
  _M0L7_2abindS798 = _M0L3endS2500 - _M0L5startS2501;
  _M0L2__S799 = 0;
  while (1) {
    if (_M0L2__S799 < _M0L7_2abindS798) {
      struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L3bufS2497 =
        _M0L3arrS792.$0;
      int32_t _M0L5startS2499 = _M0L3arrS792.$1;
      int32_t _M0L6_2atmpS2498 = _M0L5startS2499 + _M0L2__S799;
      struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L1eS800 =
        (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L3bufS2497[
          _M0L6_2atmpS2498
        ];
      moonbit_string_t _M0L6_2atmpS2494 = _M0L1eS800->$0;
      struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2495 =
        _M0L1eS800->$1;
      int32_t _M0L6_2atmpS2496;
      moonbit_incref(_M0L6_2atmpS2495);
      moonbit_incref(_M0L6_2atmpS2494);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L1mS797, _M0L6_2atmpS2494, _M0L6_2atmpS2495);
      moonbit_decref(_M0L6_2atmpS2494);
      moonbit_decref(_M0L6_2atmpS2495);
      _M0L6_2atmpS2496 = _M0L2__S799 + 1;
      _M0L2__S799 = _M0L6_2atmpS2496;
      continue;
    }
    break;
  }
  return _M0L1mS797;
}

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L3arrS803,
  int64_t _M0L8capacityS805
) {
  int32_t _M0L3endS2514;
  int32_t _M0L5startS2515;
  int32_t _M0L6lengthS802;
  int32_t _M0L8capacityS804;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L1mS808;
  int32_t _M0L3endS2511;
  int32_t _M0L5startS2512;
  int32_t _M0L7_2abindS809;
  int32_t _M0L2__S810;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2514 = _M0L3arrS803.$2;
  _M0L5startS2515 = _M0L3arrS803.$1;
  _M0L6lengthS802 = _M0L3endS2514 - _M0L5startS2515;
  if (_M0L8capacityS805 == 4294967296ll) {
    if (_M0L6lengthS802 == 0) {
      _M0L8capacityS804 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS804 = _M0FPB21capacity__for__length(_M0L6lengthS802);
    }
  } else {
    int64_t _M0L7_2aSomeS806 = _M0L8capacityS805;
    int32_t _M0L11_2acapacityS807 = (int32_t)_M0L7_2aSomeS806;
    int32_t _M0L6_2atmpS2513;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2513 = _M0FPB21capacity__for__length(_M0L6lengthS802);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS804
    = _M0MPC13int3Int3max(_M0L11_2acapacityS807, _M0L6_2atmpS2513);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS808 = _M0FPB8new__mapGiUWEuQRPC15error5ErrorNsEE(_M0L8capacityS804);
  _M0L3endS2511 = _M0L3arrS803.$2;
  _M0L5startS2512 = _M0L3arrS803.$1;
  _M0L7_2abindS809 = _M0L3endS2511 - _M0L5startS2512;
  _M0L2__S810 = 0;
  while (1) {
    if (_M0L2__S810 < _M0L7_2abindS809) {
      struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L3bufS2508 =
        _M0L3arrS803.$0;
      int32_t _M0L5startS2510 = _M0L3arrS803.$1;
      int32_t _M0L6_2atmpS2509 = _M0L5startS2510 + _M0L2__S810;
      struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L1eS811 =
        (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)_M0L3bufS2508[
          _M0L6_2atmpS2509
        ];
      int32_t _M0L6_2atmpS2505 = _M0L1eS811->$0;
      struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2atmpS2506 = _M0L1eS811->$1;
      int32_t _M0L6_2atmpS2507;
      moonbit_incref(_M0L6_2atmpS2506);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGiUWEuQRPC15error5ErrorNsEE(_M0L1mS808, _M0L6_2atmpS2505, _M0L6_2atmpS2506);
      moonbit_decref(_M0L6_2atmpS2506);
      _M0L6_2atmpS2507 = _M0L2__S810 + 1;
      _M0L2__S810 = _M0L6_2atmpS2507;
      continue;
    }
    break;
  }
  return _M0L1mS808;
}

int32_t _M0MPB3Map3setGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS785,
  moonbit_string_t _M0L3keyS786,
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L5valueS787
) {
  int32_t _M0L6_2atmpS2492;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2492 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS786);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS785, _M0L3keyS786, _M0L5valueS787, _M0L6_2atmpS2492);
  return 0;
}

int32_t _M0MPB3Map3setGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS788,
  int32_t _M0L3keyS789,
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L5valueS790
) {
  int32_t _M0L6_2atmpS2493;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2493 = _M0IPC13int3IntPB4Hash4hash(_M0L3keyS789);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS788, _M0L3keyS789, _M0L5valueS790, _M0L6_2atmpS2493);
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS756,
  moonbit_string_t _M0L3keyS762,
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L5valueS763,
  int32_t _M0L4hashS758
) {
  int32_t _M0L14capacity__maskS2473;
  int32_t _M0L6_2atmpS2472;
  int32_t _M0L3pslS753;
  int32_t _M0L3idxS754;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2473 = _M0L4selfS756->$3;
  _M0L6_2atmpS2472 = _M0L4hashS758 & _M0L14capacity__maskS2473;
  _M0L3pslS753 = 0;
  _M0L3idxS754 = _M0L6_2atmpS2472;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2471 =
      _M0L4selfS756->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS755;
    if (
      _M0L3idxS754 < 0
      || _M0L3idxS754 >= Moonbit_array_length(_M0L7entriesS2471)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS755
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2471[
        _M0L3idxS754
      ];
    if (_M0L7_2abindS755 == 0) {
      int32_t _M0L4sizeS2456 = _M0L4selfS756->$1;
      int32_t _M0L8grow__atS2457 = _M0L4selfS756->$4;
      int32_t _M0L7_2abindS759;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS760;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS761;
      if (_M0L4sizeS2456 >= _M0L8grow__atS2457) {
        int32_t _M0L14capacity__maskS2459;
        int32_t _M0L6_2atmpS2458;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS756);
        _M0L14capacity__maskS2459 = _M0L4selfS756->$3;
        _M0L6_2atmpS2458 = _M0L4hashS758 & _M0L14capacity__maskS2459;
        _M0L3pslS753 = 0;
        _M0L3idxS754 = _M0L6_2atmpS2458;
        continue;
      }
      _M0L7_2abindS759 = _M0L4selfS756->$6;
      _M0L7_2abindS760 = 0;
      moonbit_incref(_M0L3keyS762);
      moonbit_incref(_M0L5valueS763);
      _M0L5entryS761
      = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
      Moonbit_object_header(_M0L5entryS761)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 57, 0);
      _M0L5entryS761->$0 = _M0L7_2abindS759;
      _M0L5entryS761->$1 = _M0L7_2abindS760;
      _M0L5entryS761->$2 = _M0L3pslS753;
      _M0L5entryS761->$3 = _M0L4hashS758;
      _M0L5entryS761->$4 = _M0L3keyS762;
      _M0L5entryS761->$5 = _M0L5valueS763;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS756, _M0L3idxS754, _M0L5entryS761);
      moonbit_decref(_M0L5entryS761);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS764 =
        _M0L7_2abindS755;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L14_2acurr__entryS765 =
        _M0L7_2aSomeS764;
      int32_t _M0L4hashS2461 = _M0L14_2acurr__entryS765->$3;
      int32_t _if__result_3975;
      int32_t _M0L3pslS2462;
      int32_t _M0L6_2atmpS2467;
      int32_t _M0L6_2atmpS2469;
      int32_t _M0L14capacity__maskS2470;
      int32_t _M0L6_2atmpS2468;
      if (_M0L4hashS2461 == _M0L4hashS758) {
        moonbit_string_t _M0L3keyS2460 = _M0L14_2acurr__entryS765->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3975
        = _M0L3keyS2460 == _M0L3keyS762
          || Moonbit_array_length(_M0L3keyS2460)
             == Moonbit_array_length(_M0L3keyS762)
             && 0
                == memcmp(_M0L3keyS2460, _M0L3keyS762, Moonbit_array_length(_M0L3keyS2460) * 2);
      } else {
        _if__result_3975 = 0;
      }
      if (_if__result_3975) {
        struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3575 =
          _M0L14_2acurr__entryS765->$5;
        moonbit_incref(_M0L5valueS763);
        moonbit_decref(_M0L6_2aoldS3575);
        _M0L14_2acurr__entryS765->$5 = _M0L5valueS763;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS765);
      }
      _M0L3pslS2462 = _M0L14_2acurr__entryS765->$2;
      if (_M0L3pslS753 > _M0L3pslS2462) {
        int32_t _M0L4sizeS2463 = _M0L4selfS756->$1;
        int32_t _M0L8grow__atS2464 = _M0L4selfS756->$4;
        int32_t _M0L7_2abindS766;
        struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS767;
        struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS768;
        if (_M0L4sizeS2463 >= _M0L8grow__atS2464) {
          int32_t _M0L14capacity__maskS2466;
          int32_t _M0L6_2atmpS2465;
          moonbit_decref(_M0L14_2acurr__entryS765);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS756);
          _M0L14capacity__maskS2466 = _M0L4selfS756->$3;
          _M0L6_2atmpS2465 = _M0L4hashS758 & _M0L14capacity__maskS2466;
          _M0L3pslS753 = 0;
          _M0L3idxS754 = _M0L6_2atmpS2465;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS756, _M0L3idxS754, _M0L14_2acurr__entryS765);
        moonbit_decref(_M0L14_2acurr__entryS765);
        _M0L7_2abindS766 = _M0L4selfS756->$6;
        _M0L7_2abindS767 = 0;
        moonbit_incref(_M0L3keyS762);
        moonbit_incref(_M0L5valueS763);
        _M0L5entryS768
        = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
        Moonbit_object_header(_M0L5entryS768)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 57, 0);
        _M0L5entryS768->$0 = _M0L7_2abindS766;
        _M0L5entryS768->$1 = _M0L7_2abindS767;
        _M0L5entryS768->$2 = _M0L3pslS753;
        _M0L5entryS768->$3 = _M0L4hashS758;
        _M0L5entryS768->$4 = _M0L3keyS762;
        _M0L5entryS768->$5 = _M0L5valueS763;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS756, _M0L3idxS754, _M0L5entryS768);
        moonbit_decref(_M0L5entryS768);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS765);
      }
      _M0L6_2atmpS2467 = _M0L3pslS753 + 1;
      _M0L6_2atmpS2469 = _M0L3idxS754 + 1;
      _M0L14capacity__maskS2470 = _M0L4selfS756->$3;
      _M0L6_2atmpS2468 = _M0L6_2atmpS2469 & _M0L14capacity__maskS2470;
      _M0L3pslS753 = _M0L6_2atmpS2467;
      _M0L3idxS754 = _M0L6_2atmpS2468;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS772,
  int32_t _M0L3keyS778,
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L5valueS779,
  int32_t _M0L4hashS774
) {
  int32_t _M0L14capacity__maskS2491;
  int32_t _M0L6_2atmpS2490;
  int32_t _M0L3pslS769;
  int32_t _M0L3idxS770;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2491 = _M0L4selfS772->$3;
  _M0L6_2atmpS2490 = _M0L4hashS774 & _M0L14capacity__maskS2491;
  _M0L3pslS769 = 0;
  _M0L3idxS770 = _M0L6_2atmpS2490;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2489 =
      _M0L4selfS772->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS771;
    if (
      _M0L3idxS770 < 0
      || _M0L3idxS770 >= Moonbit_array_length(_M0L7entriesS2489)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS771
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2489[
        _M0L3idxS770
      ];
    if (_M0L7_2abindS771 == 0) {
      int32_t _M0L4sizeS2474 = _M0L4selfS772->$1;
      int32_t _M0L8grow__atS2475 = _M0L4selfS772->$4;
      int32_t _M0L7_2abindS775;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS776;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS777;
      if (_M0L4sizeS2474 >= _M0L8grow__atS2475) {
        int32_t _M0L14capacity__maskS2477;
        int32_t _M0L6_2atmpS2476;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS772);
        _M0L14capacity__maskS2477 = _M0L4selfS772->$3;
        _M0L6_2atmpS2476 = _M0L4hashS774 & _M0L14capacity__maskS2477;
        _M0L3pslS769 = 0;
        _M0L3idxS770 = _M0L6_2atmpS2476;
        continue;
      }
      _M0L7_2abindS775 = _M0L4selfS772->$6;
      _M0L7_2abindS776 = 0;
      moonbit_incref(_M0L5valueS779);
      _M0L5entryS777
      = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE));
      Moonbit_object_header(_M0L5entryS777)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 62, 0);
      _M0L5entryS777->$0 = _M0L7_2abindS775;
      _M0L5entryS777->$1 = _M0L7_2abindS776;
      _M0L5entryS777->$2 = _M0L3pslS769;
      _M0L5entryS777->$3 = _M0L4hashS774;
      _M0L5entryS777->$4 = _M0L3keyS778;
      _M0L5entryS777->$5 = _M0L5valueS779;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS772, _M0L3idxS770, _M0L5entryS777);
      moonbit_decref(_M0L5entryS777);
      return 0;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS780 =
        _M0L7_2abindS771;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L14_2acurr__entryS781 =
        _M0L7_2aSomeS780;
      int32_t _M0L4hashS2479 = _M0L14_2acurr__entryS781->$3;
      int32_t _if__result_3977;
      int32_t _M0L3pslS2480;
      int32_t _M0L6_2atmpS2485;
      int32_t _M0L6_2atmpS2487;
      int32_t _M0L14capacity__maskS2488;
      int32_t _M0L6_2atmpS2486;
      if (_M0L4hashS2479 == _M0L4hashS774) {
        int32_t _M0L3keyS2478 = _M0L14_2acurr__entryS781->$4;
        _if__result_3977 = _M0L3keyS2478 == _M0L3keyS778;
      } else {
        _if__result_3977 = 0;
      }
      if (_if__result_3977) {
        struct _M0TUWEuQRPC15error5ErrorNsE* _M0L6_2aoldS3579 =
          _M0L14_2acurr__entryS781->$5;
        moonbit_incref(_M0L5valueS779);
        moonbit_decref(_M0L6_2aoldS3579);
        _M0L14_2acurr__entryS781->$5 = _M0L5valueS779;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS781);
      }
      _M0L3pslS2480 = _M0L14_2acurr__entryS781->$2;
      if (_M0L3pslS769 > _M0L3pslS2480) {
        int32_t _M0L4sizeS2481 = _M0L4selfS772->$1;
        int32_t _M0L8grow__atS2482 = _M0L4selfS772->$4;
        int32_t _M0L7_2abindS782;
        struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS783;
        struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS784;
        if (_M0L4sizeS2481 >= _M0L8grow__atS2482) {
          int32_t _M0L14capacity__maskS2484;
          int32_t _M0L6_2atmpS2483;
          moonbit_decref(_M0L14_2acurr__entryS781);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS772);
          _M0L14capacity__maskS2484 = _M0L4selfS772->$3;
          _M0L6_2atmpS2483 = _M0L4hashS774 & _M0L14capacity__maskS2484;
          _M0L3pslS769 = 0;
          _M0L3idxS770 = _M0L6_2atmpS2483;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS772, _M0L3idxS770, _M0L14_2acurr__entryS781);
        moonbit_decref(_M0L14_2acurr__entryS781);
        _M0L7_2abindS782 = _M0L4selfS772->$6;
        _M0L7_2abindS783 = 0;
        moonbit_incref(_M0L5valueS779);
        _M0L5entryS784
        = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE));
        Moonbit_object_header(_M0L5entryS784)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 62, 0);
        _M0L5entryS784->$0 = _M0L7_2abindS782;
        _M0L5entryS784->$1 = _M0L7_2abindS783;
        _M0L5entryS784->$2 = _M0L3pslS769;
        _M0L5entryS784->$3 = _M0L4hashS774;
        _M0L5entryS784->$4 = _M0L3keyS778;
        _M0L5entryS784->$5 = _M0L5valueS779;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS772, _M0L3idxS770, _M0L5entryS784);
        moonbit_decref(_M0L5entryS784);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS781);
      }
      _M0L6_2atmpS2485 = _M0L3pslS769 + 1;
      _M0L6_2atmpS2487 = _M0L3idxS770 + 1;
      _M0L14capacity__maskS2488 = _M0L4selfS772->$3;
      _M0L6_2atmpS2486 = _M0L6_2atmpS2487 & _M0L14capacity__maskS2488;
      _M0L3pslS769 = _M0L6_2atmpS2485;
      _M0L3idxS770 = _M0L6_2atmpS2486;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS738
) {
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L9old__headS737;
  int32_t _M0L8capacityS2447;
  int32_t _M0L13new__capacityS739;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2441;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L6_2atmpS2440;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L6_2aoldS3585;
  int32_t _M0L6_2atmpS2442;
  int32_t _M0L8capacityS2444;
  int32_t _M0L6_2atmpS2443;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2445;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3584;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L1xS740;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS737 = _M0L4selfS738->$5;
  _M0L8capacityS2447 = _M0L4selfS738->$2;
  _M0L13new__capacityS739 = _M0L8capacityS2447 << 1;
  _M0L6_2atmpS2441 = 0;
  _M0L6_2atmpS2440
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE**)moonbit_make_ref_array(_M0L13new__capacityS739, _M0L6_2atmpS2441);
  _M0L6_2aoldS3585 = _M0L4selfS738->$0;
  if (_M0L9old__headS737) {
    moonbit_incref(_M0L9old__headS737);
  }
  moonbit_decref(_M0L6_2aoldS3585);
  _M0L4selfS738->$0 = _M0L6_2atmpS2440;
  _M0L4selfS738->$2 = _M0L13new__capacityS739;
  _M0L6_2atmpS2442 = _M0L13new__capacityS739 - 1;
  _M0L4selfS738->$3 = _M0L6_2atmpS2442;
  _M0L8capacityS2444 = _M0L4selfS738->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2443 = _M0FPB21calc__grow__threshold(_M0L8capacityS2444);
  _M0L4selfS738->$4 = _M0L6_2atmpS2443;
  _M0L4selfS738->$1 = 0;
  _M0L6_2atmpS2445 = 0;
  _M0L6_2aoldS3584 = _M0L4selfS738->$5;
  if (_M0L6_2aoldS3584) {
    moonbit_decref(_M0L6_2aoldS3584);
  }
  _M0L4selfS738->$5 = _M0L6_2atmpS2445;
  _M0L4selfS738->$6 = -1;
  _M0L1xS740 = _M0L9old__headS737;
  while (1) {
    if (_M0L1xS740 == 0) {
      if (_M0L1xS740) {
        moonbit_decref(_M0L1xS740);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS742 =
        _M0L1xS740;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4_2aeS743 =
        _M0L7_2aSomeS742;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L15next__in__chainS744 =
        _M0L4_2aeS743->$1;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2446 =
        0;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3582 =
        _M0L4_2aeS743->$1;
      if (_M0L15next__in__chainS744) {
        moonbit_incref(_M0L15next__in__chainS744);
      }
      if (_M0L6_2aoldS3582) {
        moonbit_decref(_M0L6_2aoldS3582);
      }
      _M0L4_2aeS743->$1 = _M0L6_2atmpS2446;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS738, _M0L4_2aeS743);
      moonbit_decref(_M0L4_2aeS743);
      _M0L1xS740 = _M0L15next__in__chainS744;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS746
) {
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L9old__headS745;
  int32_t _M0L8capacityS2455;
  int32_t _M0L13new__capacityS747;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2449;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS2448;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L6_2aoldS3590;
  int32_t _M0L6_2atmpS2450;
  int32_t _M0L8capacityS2452;
  int32_t _M0L6_2atmpS2451;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2453;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3589;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L1xS748;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS745 = _M0L4selfS746->$5;
  _M0L8capacityS2455 = _M0L4selfS746->$2;
  _M0L13new__capacityS747 = _M0L8capacityS2455 << 1;
  _M0L6_2atmpS2449 = 0;
  _M0L6_2atmpS2448
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE**)moonbit_make_ref_array(_M0L13new__capacityS747, _M0L6_2atmpS2449);
  _M0L6_2aoldS3590 = _M0L4selfS746->$0;
  if (_M0L9old__headS745) {
    moonbit_incref(_M0L9old__headS745);
  }
  moonbit_decref(_M0L6_2aoldS3590);
  _M0L4selfS746->$0 = _M0L6_2atmpS2448;
  _M0L4selfS746->$2 = _M0L13new__capacityS747;
  _M0L6_2atmpS2450 = _M0L13new__capacityS747 - 1;
  _M0L4selfS746->$3 = _M0L6_2atmpS2450;
  _M0L8capacityS2452 = _M0L4selfS746->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2451 = _M0FPB21calc__grow__threshold(_M0L8capacityS2452);
  _M0L4selfS746->$4 = _M0L6_2atmpS2451;
  _M0L4selfS746->$1 = 0;
  _M0L6_2atmpS2453 = 0;
  _M0L6_2aoldS3589 = _M0L4selfS746->$5;
  if (_M0L6_2aoldS3589) {
    moonbit_decref(_M0L6_2aoldS3589);
  }
  _M0L4selfS746->$5 = _M0L6_2atmpS2453;
  _M0L4selfS746->$6 = -1;
  _M0L1xS748 = _M0L9old__headS745;
  while (1) {
    if (_M0L1xS748 == 0) {
      if (_M0L1xS748) {
        moonbit_decref(_M0L1xS748);
      }
      break;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS750 =
        _M0L1xS748;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L4_2aeS751 =
        _M0L7_2aSomeS750;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L15next__in__chainS752 =
        _M0L4_2aeS751->$1;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2454 = 0;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3587 =
        _M0L4_2aeS751->$1;
      if (_M0L15next__in__chainS752) {
        moonbit_incref(_M0L15next__in__chainS752);
      }
      if (_M0L6_2aoldS3587) {
        moonbit_decref(_M0L6_2aoldS3587);
      }
      _M0L4_2aeS751->$1 = _M0L6_2atmpS2454;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS746, _M0L4_2aeS751);
      moonbit_decref(_M0L4_2aeS751);
      _M0L1xS748 = _M0L15next__in__chainS752;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS724,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5outerS720
) {
  int32_t _M0L4hashS719;
  int32_t _M0L14capacity__maskS2429;
  int32_t _M0L6_2atmpS2428;
  int32_t _M0L3pslS721;
  int32_t _M0L3idxS722;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS719 = _M0L5outerS720->$3;
  _M0L14capacity__maskS2429 = _M0L4selfS724->$3;
  _M0L6_2atmpS2428 = _M0L4hashS719 & _M0L14capacity__maskS2429;
  _M0L3pslS721 = 0;
  _M0L3idxS722 = _M0L6_2atmpS2428;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2427 =
      _M0L4selfS724->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS723;
    if (
      _M0L3idxS722 < 0
      || _M0L3idxS722 >= Moonbit_array_length(_M0L7entriesS2427)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS723
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2427[
        _M0L3idxS722
      ];
    if (_M0L7_2abindS723 == 0) {
      int32_t _M0L4tailS2420;
      _M0L5outerS720->$2 = _M0L3pslS721;
      _M0L4tailS2420 = _M0L4selfS724->$6;
      _M0L5outerS720->$0 = _M0L4tailS2420;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS724, _M0L3idxS722, _M0L5outerS720);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS725 =
        _M0L7_2abindS723;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2acurrS726 =
        _M0L7_2aSomeS725;
      int32_t _M0L3pslS2421 = _M0L7_2acurrS726->$2;
      if (_M0L3pslS721 > _M0L3pslS2421) {
        int32_t _M0L4tailS2422;
        moonbit_incref(_M0L7_2acurrS726);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS724, _M0L3idxS722, _M0L7_2acurrS726);
        moonbit_decref(_M0L7_2acurrS726);
        _M0L5outerS720->$2 = _M0L3pslS721;
        _M0L4tailS2422 = _M0L4selfS724->$6;
        _M0L5outerS720->$0 = _M0L4tailS2422;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS724, _M0L3idxS722, _M0L5outerS720);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2423 = _M0L3pslS721 + 1;
        int32_t _M0L6_2atmpS2425 = _M0L3idxS722 + 1;
        int32_t _M0L14capacity__maskS2426 = _M0L4selfS724->$3;
        int32_t _M0L6_2atmpS2424 =
          _M0L6_2atmpS2425 & _M0L14capacity__maskS2426;
        _M0L3pslS721 = _M0L6_2atmpS2423;
        _M0L3idxS722 = _M0L6_2atmpS2424;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS733,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5outerS729
) {
  int32_t _M0L4hashS728;
  int32_t _M0L14capacity__maskS2439;
  int32_t _M0L6_2atmpS2438;
  int32_t _M0L3pslS730;
  int32_t _M0L3idxS731;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS728 = _M0L5outerS729->$3;
  _M0L14capacity__maskS2439 = _M0L4selfS733->$3;
  _M0L6_2atmpS2438 = _M0L4hashS728 & _M0L14capacity__maskS2439;
  _M0L3pslS730 = 0;
  _M0L3idxS731 = _M0L6_2atmpS2438;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2437 =
      _M0L4selfS733->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS732;
    if (
      _M0L3idxS731 < 0
      || _M0L3idxS731 >= Moonbit_array_length(_M0L7entriesS2437)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS732
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2437[
        _M0L3idxS731
      ];
    if (_M0L7_2abindS732 == 0) {
      int32_t _M0L4tailS2430;
      _M0L5outerS729->$2 = _M0L3pslS730;
      _M0L4tailS2430 = _M0L4selfS733->$6;
      _M0L5outerS729->$0 = _M0L4tailS2430;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS733, _M0L3idxS731, _M0L5outerS729);
      return 0;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS734 =
        _M0L7_2abindS732;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2acurrS735 =
        _M0L7_2aSomeS734;
      int32_t _M0L3pslS2431 = _M0L7_2acurrS735->$2;
      if (_M0L3pslS730 > _M0L3pslS2431) {
        int32_t _M0L4tailS2432;
        moonbit_incref(_M0L7_2acurrS735);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS733, _M0L3idxS731, _M0L7_2acurrS735);
        moonbit_decref(_M0L7_2acurrS735);
        _M0L5outerS729->$2 = _M0L3pslS730;
        _M0L4tailS2432 = _M0L4selfS733->$6;
        _M0L5outerS729->$0 = _M0L4tailS2432;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS733, _M0L3idxS731, _M0L5outerS729);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2433 = _M0L3pslS730 + 1;
        int32_t _M0L6_2atmpS2435 = _M0L3idxS731 + 1;
        int32_t _M0L14capacity__maskS2436 = _M0L4selfS733->$3;
        int32_t _M0L6_2atmpS2434 =
          _M0L6_2atmpS2435 & _M0L14capacity__maskS2436;
        _M0L3pslS730 = _M0L6_2atmpS2433;
        _M0L3idxS731 = _M0L6_2atmpS2434;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS703,
  int32_t _M0L3idxS708,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS707
) {
  int32_t _M0L3pslS2403;
  int32_t _M0L6_2atmpS2399;
  int32_t _M0L6_2atmpS2401;
  int32_t _M0L14capacity__maskS2402;
  int32_t _M0L6_2atmpS2400;
  int32_t _M0L3pslS699;
  int32_t _M0L3idxS700;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS701;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2403 = _M0L5entryS707->$2;
  _M0L6_2atmpS2399 = _M0L3pslS2403 + 1;
  _M0L6_2atmpS2401 = _M0L3idxS708 + 1;
  _M0L14capacity__maskS2402 = _M0L4selfS703->$3;
  _M0L6_2atmpS2400 = _M0L6_2atmpS2401 & _M0L14capacity__maskS2402;
  moonbit_incref(_M0L5entryS707);
  _M0L3pslS699 = _M0L6_2atmpS2399;
  _M0L3idxS700 = _M0L6_2atmpS2400;
  _M0L5entryS701 = _M0L5entryS707;
  while (1) {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2398 =
      _M0L4selfS703->$0;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS702;
    if (
      _M0L3idxS700 < 0
      || _M0L3idxS700 >= Moonbit_array_length(_M0L7entriesS2398)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS702
    = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2398[
        _M0L3idxS700
      ];
    if (_M0L7_2abindS702 == 0) {
      _M0L5entryS701->$2 = _M0L3pslS699;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS703, _M0L5entryS701, _M0L3idxS700);
      moonbit_decref(_M0L5entryS701);
      break;
    } else {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS705 =
        _M0L7_2abindS702;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L14_2acurr__entryS706 =
        _M0L7_2aSomeS705;
      int32_t _M0L3pslS2388 = _M0L14_2acurr__entryS706->$2;
      if (_M0L3pslS699 > _M0L3pslS2388) {
        int32_t _M0L3pslS2393;
        int32_t _M0L6_2atmpS2389;
        int32_t _M0L6_2atmpS2391;
        int32_t _M0L14capacity__maskS2392;
        int32_t _M0L6_2atmpS2390;
        _M0L5entryS701->$2 = _M0L3pslS699;
        moonbit_incref(_M0L14_2acurr__entryS706);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L4selfS703, _M0L5entryS701, _M0L3idxS700);
        moonbit_decref(_M0L5entryS701);
        _M0L3pslS2393 = _M0L14_2acurr__entryS706->$2;
        _M0L6_2atmpS2389 = _M0L3pslS2393 + 1;
        _M0L6_2atmpS2391 = _M0L3idxS700 + 1;
        _M0L14capacity__maskS2392 = _M0L4selfS703->$3;
        _M0L6_2atmpS2390 = _M0L6_2atmpS2391 & _M0L14capacity__maskS2392;
        _M0L3pslS699 = _M0L6_2atmpS2389;
        _M0L3idxS700 = _M0L6_2atmpS2390;
        _M0L5entryS701 = _M0L14_2acurr__entryS706;
        continue;
      } else {
        int32_t _M0L6_2atmpS2394 = _M0L3pslS699 + 1;
        int32_t _M0L6_2atmpS2396 = _M0L3idxS700 + 1;
        int32_t _M0L14capacity__maskS2397 = _M0L4selfS703->$3;
        int32_t _M0L6_2atmpS2395 =
          _M0L6_2atmpS2396 & _M0L14capacity__maskS2397;
        struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _tmp_3983 =
          _M0L5entryS701;
        _M0L3pslS699 = _M0L6_2atmpS2394;
        _M0L3idxS700 = _M0L6_2atmpS2395;
        _M0L5entryS701 = _tmp_3983;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS713,
  int32_t _M0L3idxS718,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS717
) {
  int32_t _M0L3pslS2419;
  int32_t _M0L6_2atmpS2415;
  int32_t _M0L6_2atmpS2417;
  int32_t _M0L14capacity__maskS2418;
  int32_t _M0L6_2atmpS2416;
  int32_t _M0L3pslS709;
  int32_t _M0L3idxS710;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS711;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2419 = _M0L5entryS717->$2;
  _M0L6_2atmpS2415 = _M0L3pslS2419 + 1;
  _M0L6_2atmpS2417 = _M0L3idxS718 + 1;
  _M0L14capacity__maskS2418 = _M0L4selfS713->$3;
  _M0L6_2atmpS2416 = _M0L6_2atmpS2417 & _M0L14capacity__maskS2418;
  moonbit_incref(_M0L5entryS717);
  _M0L3pslS709 = _M0L6_2atmpS2415;
  _M0L3idxS710 = _M0L6_2atmpS2416;
  _M0L5entryS711 = _M0L5entryS717;
  while (1) {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2414 =
      _M0L4selfS713->$0;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS712;
    if (
      _M0L3idxS710 < 0
      || _M0L3idxS710 >= Moonbit_array_length(_M0L7entriesS2414)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS712
    = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2414[
        _M0L3idxS710
      ];
    if (_M0L7_2abindS712 == 0) {
      _M0L5entryS711->$2 = _M0L3pslS709;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS713, _M0L5entryS711, _M0L3idxS710);
      moonbit_decref(_M0L5entryS711);
      break;
    } else {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS715 =
        _M0L7_2abindS712;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L14_2acurr__entryS716 =
        _M0L7_2aSomeS715;
      int32_t _M0L3pslS2404 = _M0L14_2acurr__entryS716->$2;
      if (_M0L3pslS709 > _M0L3pslS2404) {
        int32_t _M0L3pslS2409;
        int32_t _M0L6_2atmpS2405;
        int32_t _M0L6_2atmpS2407;
        int32_t _M0L14capacity__maskS2408;
        int32_t _M0L6_2atmpS2406;
        _M0L5entryS711->$2 = _M0L3pslS709;
        moonbit_incref(_M0L14_2acurr__entryS716);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(_M0L4selfS713, _M0L5entryS711, _M0L3idxS710);
        moonbit_decref(_M0L5entryS711);
        _M0L3pslS2409 = _M0L14_2acurr__entryS716->$2;
        _M0L6_2atmpS2405 = _M0L3pslS2409 + 1;
        _M0L6_2atmpS2407 = _M0L3idxS710 + 1;
        _M0L14capacity__maskS2408 = _M0L4selfS713->$3;
        _M0L6_2atmpS2406 = _M0L6_2atmpS2407 & _M0L14capacity__maskS2408;
        _M0L3pslS709 = _M0L6_2atmpS2405;
        _M0L3idxS710 = _M0L6_2atmpS2406;
        _M0L5entryS711 = _M0L14_2acurr__entryS716;
        continue;
      } else {
        int32_t _M0L6_2atmpS2410 = _M0L3pslS709 + 1;
        int32_t _M0L6_2atmpS2412 = _M0L3idxS710 + 1;
        int32_t _M0L14capacity__maskS2413 = _M0L4selfS713->$3;
        int32_t _M0L6_2atmpS2411 =
          _M0L6_2atmpS2412 & _M0L14capacity__maskS2413;
        struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _tmp_3985 =
          _M0L5entryS711;
        _M0L3pslS709 = _M0L6_2atmpS2410;
        _M0L3idxS710 = _M0L6_2atmpS2411;
        _M0L5entryS711 = _tmp_3985;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS687,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS689,
  int32_t _M0L8new__idxS688
) {
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2384;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2385;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3601;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS690;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2384 = _M0L4selfS687->$0;
  _M0L6_2atmpS2385 = _M0L5entryS689;
  if (
    _M0L8new__idxS688 < 0
    || _M0L8new__idxS688 >= Moonbit_array_length(_M0L7entriesS2384)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3601
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2384[
      _M0L8new__idxS688
    ];
  if (_M0L6_2atmpS2385) {
    moonbit_incref(_M0L6_2atmpS2385);
  }
  if (_M0L6_2aoldS3601) {
    moonbit_decref(_M0L6_2aoldS3601);
  }
  _M0L7entriesS2384[_M0L8new__idxS688] = _M0L6_2atmpS2385;
  _M0L7_2abindS690 = _M0L5entryS689->$1;
  if (_M0L7_2abindS690 == 0) {
    _M0L4selfS687->$6 = _M0L8new__idxS688;
  } else {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS691 =
      _M0L7_2abindS690;
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2anextS692 =
      _M0L7_2aSomeS691;
    _M0L7_2anextS692->$0 = _M0L8new__idxS688;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS693,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS695,
  int32_t _M0L8new__idxS694
) {
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2386;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2387;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3604;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS696;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2386 = _M0L4selfS693->$0;
  _M0L6_2atmpS2387 = _M0L5entryS695;
  if (
    _M0L8new__idxS694 < 0
    || _M0L8new__idxS694 >= Moonbit_array_length(_M0L7entriesS2386)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3604
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2386[
      _M0L8new__idxS694
    ];
  if (_M0L6_2atmpS2387) {
    moonbit_incref(_M0L6_2atmpS2387);
  }
  if (_M0L6_2aoldS3604) {
    moonbit_decref(_M0L6_2aoldS3604);
  }
  _M0L7entriesS2386[_M0L8new__idxS694] = _M0L6_2atmpS2387;
  _M0L7_2abindS696 = _M0L5entryS695->$1;
  if (_M0L7_2abindS696 == 0) {
    _M0L4selfS693->$6 = _M0L8new__idxS694;
  } else {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS697 =
      _M0L7_2abindS696;
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2anextS698 =
      _M0L7_2aSomeS697;
    _M0L7_2anextS698->$0 = _M0L8new__idxS694;
  }
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS680,
  int32_t _M0L3idxS682,
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L5entryS681
) {
  int32_t _M0L7_2abindS679;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2371;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2372;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3606;
  int32_t _M0L4sizeS2374;
  int32_t _M0L6_2atmpS2373;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS679 = _M0L4selfS680->$6;
  switch (_M0L7_2abindS679) {
    case -1: {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2366 =
        _M0L5entryS681;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3608 =
        _M0L4selfS680->$5;
      if (_M0L6_2atmpS2366) {
        moonbit_incref(_M0L6_2atmpS2366);
      }
      if (_M0L6_2aoldS3608) {
        moonbit_decref(_M0L6_2aoldS3608);
      }
      _M0L4selfS680->$5 = _M0L6_2atmpS2366;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7entriesS2370 =
        _M0L4selfS680->$0;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2369;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2367;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2368;
      struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2aoldS3609;
      if (
        _M0L7_2abindS679 < 0
        || _M0L7_2abindS679 >= Moonbit_array_length(_M0L7entriesS2370)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2369
      = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2370[
          _M0L7_2abindS679
        ];
      if (_M0L6_2atmpS2369) {
        moonbit_incref(_M0L6_2atmpS2369);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2367
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE(_M0L6_2atmpS2369);
      if (_M0L6_2atmpS2369) {
        moonbit_decref(_M0L6_2atmpS2369);
      }
      _M0L6_2atmpS2368 = _M0L5entryS681;
      _M0L6_2aoldS3609 = _M0L6_2atmpS2367->$1;
      if (_M0L6_2atmpS2368) {
        moonbit_incref(_M0L6_2atmpS2368);
      }
      if (_M0L6_2aoldS3609) {
        moonbit_decref(_M0L6_2aoldS3609);
      }
      _M0L6_2atmpS2367->$1 = _M0L6_2atmpS2368;
      moonbit_decref(_M0L6_2atmpS2367);
      break;
    }
  }
  _M0L4selfS680->$6 = _M0L3idxS682;
  _M0L7entriesS2371 = _M0L4selfS680->$0;
  _M0L6_2atmpS2372 = _M0L5entryS681;
  if (
    _M0L3idxS682 < 0
    || _M0L3idxS682 >= Moonbit_array_length(_M0L7entriesS2371)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3606
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)_M0L7entriesS2371[
      _M0L3idxS682
    ];
  if (_M0L6_2atmpS2372) {
    moonbit_incref(_M0L6_2atmpS2372);
  }
  if (_M0L6_2aoldS3606) {
    moonbit_decref(_M0L6_2aoldS3606);
  }
  _M0L7entriesS2371[_M0L3idxS682] = _M0L6_2atmpS2372;
  _M0L4sizeS2374 = _M0L4selfS680->$1;
  _M0L6_2atmpS2373 = _M0L4sizeS2374 + 1;
  _M0L4selfS680->$1 = _M0L6_2atmpS2373;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGiUWEuQRPC15error5ErrorNsEE(
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS684,
  int32_t _M0L3idxS686,
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L5entryS685
) {
  int32_t _M0L7_2abindS683;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2380;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2381;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3612;
  int32_t _M0L4sizeS2383;
  int32_t _M0L6_2atmpS2382;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS683 = _M0L4selfS684->$6;
  switch (_M0L7_2abindS683) {
    case -1: {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2375 =
        _M0L5entryS685;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3614 =
        _M0L4selfS684->$5;
      if (_M0L6_2atmpS2375) {
        moonbit_incref(_M0L6_2atmpS2375);
      }
      if (_M0L6_2aoldS3614) {
        moonbit_decref(_M0L6_2aoldS3614);
      }
      _M0L4selfS684->$5 = _M0L6_2atmpS2375;
      break;
    }
    default: {
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7entriesS2379 =
        _M0L4selfS684->$0;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2378;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2376;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2377;
      struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2aoldS3615;
      if (
        _M0L7_2abindS683 < 0
        || _M0L7_2abindS683 >= Moonbit_array_length(_M0L7entriesS2379)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2378
      = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2379[
          _M0L7_2abindS683
        ];
      if (_M0L6_2atmpS2378) {
        moonbit_incref(_M0L6_2atmpS2378);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS2376
      = _M0MPC16option6Option6unwrapGRPB5EntryGiUWEuQRPC15error5ErrorNsEEE(_M0L6_2atmpS2378);
      if (_M0L6_2atmpS2378) {
        moonbit_decref(_M0L6_2atmpS2378);
      }
      _M0L6_2atmpS2377 = _M0L5entryS685;
      _M0L6_2aoldS3615 = _M0L6_2atmpS2376->$1;
      if (_M0L6_2atmpS2377) {
        moonbit_incref(_M0L6_2atmpS2377);
      }
      if (_M0L6_2aoldS3615) {
        moonbit_decref(_M0L6_2aoldS3615);
      }
      _M0L6_2atmpS2376->$1 = _M0L6_2atmpS2377;
      moonbit_decref(_M0L6_2atmpS2376);
      break;
    }
  }
  _M0L4selfS684->$6 = _M0L3idxS686;
  _M0L7entriesS2380 = _M0L4selfS684->$0;
  _M0L6_2atmpS2381 = _M0L5entryS685;
  if (
    _M0L3idxS686 < 0
    || _M0L3idxS686 >= Moonbit_array_length(_M0L7entriesS2380)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3612
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE*)_M0L7entriesS2380[
      _M0L3idxS686
    ];
  if (_M0L6_2atmpS2381) {
    moonbit_incref(_M0L6_2atmpS2381);
  }
  if (_M0L6_2aoldS3612) {
    moonbit_decref(_M0L6_2aoldS3612);
  }
  _M0L7entriesS2380[_M0L3idxS686] = _M0L6_2atmpS2381;
  _M0L4sizeS2383 = _M0L4selfS684->$1;
  _M0L6_2atmpS2382 = _M0L4sizeS2383 + 1;
  _M0L4selfS684->$1 = _M0L6_2atmpS2382;
  return 0;
}

int32_t _M0MPC13int3Int3max(int32_t _M0L4selfS677, int32_t _M0L5otherS678) {
  #line 75 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS677 > _M0L5otherS678) {
    return _M0L4selfS677;
  } else {
    return _M0L5otherS678;
  }
}

int32_t _M0FPB21capacity__for__length(int32_t _M0L6lengthS676) {
  int32_t _M0Lm8capacityS675;
  int32_t _M0L6_2atmpS2364;
  int32_t _M0L6_2atmpS2363;
  #line 71 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 72 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0Lm8capacityS675 = _M0MPC13int3Int20next__power__of__two(_M0L6lengthS676);
  _M0L6_2atmpS2364 = _M0Lm8capacityS675;
  #line 73 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2363 = _M0FPB21calc__grow__threshold(_M0L6_2atmpS2364);
  if (_M0L6lengthS676 > _M0L6_2atmpS2363) {
    int32_t _M0L6_2atmpS2365 = _M0Lm8capacityS675;
    _M0Lm8capacityS675 = _M0L6_2atmpS2365 * 2;
  }
  return _M0Lm8capacityS675;
}

struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0FPB8new__mapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(
  int32_t _M0L8capacityS664
) {
  int32_t _M0L8capacityS663;
  int32_t _M0L7_2abindS665;
  int32_t _M0L7_2abindS666;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L6_2atmpS2361;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7_2abindS667;
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2abindS668;
  struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _block_3986;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS663
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS664);
  _M0L7_2abindS665 = _M0L8capacityS663 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS666 = _M0FPB21calc__grow__threshold(_M0L8capacityS663);
  _M0L6_2atmpS2361 = 0;
  _M0L7_2abindS667
  = (struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE**)moonbit_make_ref_array(_M0L8capacityS663, _M0L6_2atmpS2361);
  _M0L7_2abindS668 = 0;
  _block_3986
  = (struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_block_3986)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 66, 0);
  _block_3986->$0 = _M0L7_2abindS667;
  _block_3986->$1 = 0;
  _block_3986->$2 = _M0L8capacityS663;
  _block_3986->$3 = _M0L7_2abindS665;
  _block_3986->$4 = _M0L7_2abindS666;
  _block_3986->$5 = _M0L7_2abindS668;
  _block_3986->$6 = -1;
  return _block_3986;
}

struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0FPB8new__mapGiUWEuQRPC15error5ErrorNsEE(
  int32_t _M0L8capacityS670
) {
  int32_t _M0L8capacityS669;
  int32_t _M0L7_2abindS671;
  int32_t _M0L7_2abindS672;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS2362;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS673;
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2abindS674;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _block_3987;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS669
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS670);
  _M0L7_2abindS671 = _M0L8capacityS669 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS672 = _M0FPB21calc__grow__threshold(_M0L8capacityS669);
  _M0L6_2atmpS2362 = 0;
  _M0L7_2abindS673
  = (struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE**)moonbit_make_ref_array(_M0L8capacityS669, _M0L6_2atmpS2362);
  _M0L7_2abindS674 = 0;
  _block_3987
  = (struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_block_3987)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 70, 0);
  _block_3987->$0 = _M0L7_2abindS673;
  _block_3987->$1 = 0;
  _block_3987->$2 = _M0L8capacityS669;
  _block_3987->$3 = _M0L7_2abindS671;
  _block_3987->$4 = _M0L7_2abindS672;
  _block_3987->$5 = _M0L7_2abindS674;
  _block_3987->$6 = -1;
  return _block_3987;
}

int32_t _M0MPC13int3Int20next__power__of__two(int32_t _M0L4selfS662) {
  #line 33 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS662 >= 0) {
    int32_t _M0L6_2atmpS2360;
    int32_t _M0L6_2atmpS2359;
    int32_t _M0L6_2atmpS2358;
    int32_t _M0L6_2atmpS2357;
    if (_M0L4selfS662 <= 1) {
      return 1;
    }
    if (_M0L4selfS662 > 1073741824) {
      return 1073741824;
    }
    _M0L6_2atmpS2360 = _M0L4selfS662 - 1;
    #line 44 "/home/developer/.moon/lib/core/builtin/int.mbt"
    _M0L6_2atmpS2359 = moonbit_clz32(_M0L6_2atmpS2360);
    _M0L6_2atmpS2358 = _M0L6_2atmpS2359 - 1;
    _M0L6_2atmpS2357 = 2147483647 >> (_M0L6_2atmpS2358 & 31);
    return _M0L6_2atmpS2357 + 1;
  } else {
    #line 34 "/home/developer/.moon/lib/core/builtin/int.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB21calc__grow__threshold(int32_t _M0L8capacityS661) {
  int32_t _M0L6_2atmpS2356;
  #line 610 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2356 = _M0L8capacityS661 * 13;
  return _M0L6_2atmpS2356 / 16;
}

double _M0MPC16option6Option10unwrap__orGdE(
  void* _M0L4selfS658,
  double _M0L7defaultS659
) {
  #line 47 "/home/developer/.moon/lib/core/builtin/option.mbt"
  switch (Moonbit_object_tag(_M0L4selfS658)) {
    case 0: {
      return _M0L7defaultS659;
      break;
    }
    default: {
      struct _M0DTPC16option6OptionGdE4Some* _M0L7_2aSomeS660 =
        (struct _M0DTPC16option6OptionGdE4Some*)_M0L4selfS658;
      return _M0L7_2aSomeS660->$0;
      break;
    }
  }
}

struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE(
  struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L4selfS654
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS654 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L7_2aSomeS655 =
      _M0L4selfS654;
    if (_M0L7_2aSomeS655) {
      moonbit_incref(_M0L7_2aSomeS655);
    }
    return _M0L7_2aSomeS655;
  }
}

struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0MPC16option6Option6unwrapGRPB5EntryGiUWEuQRPC15error5ErrorNsEEE(
  struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L4selfS656
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS656 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGiUWEuQRPC15error5ErrorNsEE* _M0L7_2aSomeS657 =
      _M0L4selfS656;
    if (_M0L7_2aSomeS657) {
      moonbit_incref(_M0L7_2aSomeS657);
    }
    return _M0L7_2aSomeS657;
  }
}

int32_t _M0IPC16option6OptionPB2Eq5equalGsE(
  moonbit_string_t _M0L4selfS648,
  moonbit_string_t _M0L5otherS649
) {
  #line 16 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS648 == 0) {
    return _M0L5otherS649 == 0;
  } else {
    moonbit_string_t _M0L7_2aSomeS650 = _M0L4selfS648;
    moonbit_string_t _M0L4_2axS651 = _M0L7_2aSomeS650;
    if (_M0L5otherS649 == 0) {
      return 0;
    } else {
      moonbit_string_t _M0L7_2aSomeS652 = _M0L5otherS649;
      moonbit_string_t _M0L4_2ayS653 = _M0L7_2aSomeS652;
      #line 19 "/home/developer/.moon/lib/core/builtin/option.mbt"
      return _M0L4_2axS651 == _M0L4_2ayS653
             || Moonbit_array_length(_M0L4_2axS651)
                == Moonbit_array_length(_M0L4_2ayS653)
                && 0
                   == memcmp(_M0L4_2axS651, _M0L4_2ayS653, Moonbit_array_length(_M0L4_2axS651) * 2);
    }
  }
}

struct _M0TPB4IterGsE* _M0MPC15array13ReadOnlyArray4iterGsE(
  moonbit_string_t* _M0L4selfS647
) {
  moonbit_string_t* _M0L6_2atmpS2355;
  struct _M0TPB4IterGsE* _result_3988;
  #line 194 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
  moonbit_incref(_M0L4selfS647);
  _M0L6_2atmpS2355 = _M0L4selfS647;
  #line 196 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
  _result_3988 = _M0MPC15array10FixedArray4iterGsE(_M0L6_2atmpS2355);
  moonbit_decref(_M0L6_2atmpS2355);
  return _result_3988;
}

uint64_t _M0MPC15array13ReadOnlyArray2atGmE(
  uint64_t* _M0L4selfS643,
  int32_t _M0L5indexS644
) {
  uint64_t* _M0L6_2atmpS2353;
  uint64_t _result_3989;
  #line 38 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
  moonbit_incref(_M0L4selfS643);
  _M0L6_2atmpS2353 = _M0L4selfS643;
  if (
    _M0L5indexS644 < 0
    || _M0L5indexS644 >= Moonbit_array_length(_M0L6_2atmpS2353)
  ) {
    #line 40 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
    moonbit_panic();
  }
  _result_3989 = (uint64_t)_M0L6_2atmpS2353[_M0L5indexS644];
  moonbit_decref(_M0L6_2atmpS2353);
  return _result_3989;
}

uint32_t _M0MPC15array13ReadOnlyArray2atGjE(
  uint32_t* _M0L4selfS645,
  int32_t _M0L5indexS646
) {
  uint32_t* _M0L6_2atmpS2354;
  uint32_t _result_3990;
  #line 38 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
  moonbit_incref(_M0L4selfS645);
  _M0L6_2atmpS2354 = _M0L4selfS645;
  if (
    _M0L5indexS646 < 0
    || _M0L5indexS646 >= Moonbit_array_length(_M0L6_2atmpS2354)
  ) {
    #line 40 "/home/developer/.moon/lib/core/builtin/readonlyarray.mbt"
    moonbit_panic();
  }
  _result_3990 = (uint32_t)_M0L6_2atmpS2354[_M0L5indexS646];
  moonbit_decref(_M0L6_2atmpS2354);
  return _result_3990;
}

struct _M0TPB4IterGsE* _M0MPC15array10FixedArray4iterGsE(
  moonbit_string_t* _M0L4selfS642
) {
  moonbit_string_t* _M0L6_2atmpS2351;
  int32_t _M0L6_2atmpS2352;
  struct _M0TPB9ArrayViewGsE _M0L6_2atmpS2350;
  struct _M0TPB4IterGsE* _result_3991;
  #line 1566 "/home/developer/.moon/lib/core/builtin/fixedarray.mbt"
  _M0L6_2atmpS2351 = _M0L4selfS642;
  _M0L6_2atmpS2352 = Moonbit_array_length(_M0L4selfS642);
  moonbit_incref(_M0L6_2atmpS2351);
  _M0L6_2atmpS2350
  = (struct _M0TPB9ArrayViewGsE){
    .$0 = _M0L6_2atmpS2351, .$1 = 0, .$2 = _M0L6_2atmpS2352
  };
  #line 1568 "/home/developer/.moon/lib/core/builtin/fixedarray.mbt"
  _result_3991 = _M0MPC15array9ArrayView4iterGsE(_M0L6_2atmpS2350);
  moonbit_decref(_M0L6_2atmpS2350.$0);
  return _result_3991;
}

struct _M0TPB4IterGsE* _M0MPC15array9ArrayView4iterGsE(
  struct _M0TPB9ArrayViewGsE _M0L4selfS640
) {
  struct _M0TPB8MutLocalGiE* _M0L1iS638;
  int32_t _M0L3endS2348;
  int32_t _M0L5startS2349;
  int32_t _M0L3lenS639;
  struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__* _closure_3992;
  struct _M0TWEOs* _M0L6_2atmpS2336;
  int64_t _M0L6_2atmpS2337;
  struct _M0TPB4IterGsE* _result_3993;
  #line 724 "/home/developer/.moon/lib/core/builtin/arrayview.mbt"
  _M0L1iS638
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS638)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS638->$0 = 0;
  _M0L3endS2348 = _M0L4selfS640.$2;
  _M0L5startS2349 = _M0L4selfS640.$1;
  _M0L3lenS639 = _M0L3endS2348 - _M0L5startS2349;
  moonbit_incref(_M0L4selfS640.$0);
  _closure_3992
  = (struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__*)moonbit_malloc(sizeof(struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__));
  Moonbit_object_header(_closure_3992)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 77, 0);
  _closure_3992->code = &_M0MPC15array9ArrayView4iterGsEC2338l729;
  _closure_3992->$0 = _M0L4selfS640;
  _closure_3992->$1 = _M0L3lenS639;
  _closure_3992->$2 = _M0L1iS638;
  _M0L6_2atmpS2336 = (struct _M0TWEOs*)_closure_3992;
  _M0L6_2atmpS2337 = (int64_t)_M0L3lenS639;
  #line 728 "/home/developer/.moon/lib/core/builtin/arrayview.mbt"
  _result_3993 = _M0MPB4Iter3newGsE(_M0L6_2atmpS2336, _M0L6_2atmpS2337);
  moonbit_decref(_M0L6_2atmpS2336);
  return _result_3993;
}

moonbit_string_t _M0MPC15array9ArrayView4iterGsEC2338l729(
  struct _M0TWEOs* _M0L6_2aenvS2339
) {
  struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__* _M0L14_2acasted__envS2340;
  struct _M0TPB8MutLocalGiE* _M0L1iS638;
  int32_t _M0L3lenS639;
  struct _M0TPB9ArrayViewGsE _M0L4selfS640;
  int32_t _M0L3valS2341;
  #line 729 "/home/developer/.moon/lib/core/builtin/arrayview.mbt"
  _M0L14_2acasted__envS2340
  = (struct _M0R59ArrayView_3a_3aiter_7c_5bString_5d_7c_2eanon__u2338__l729__*)_M0L6_2aenvS2339;
  _M0L1iS638 = _M0L14_2acasted__envS2340->$2;
  _M0L3lenS639 = _M0L14_2acasted__envS2340->$1;
  _M0L4selfS640 = _M0L14_2acasted__envS2340->$0;
  _M0L3valS2341 = _M0L1iS638->$0;
  if (_M0L3valS2341 < _M0L3lenS639) {
    moonbit_string_t* _M0L3bufS2344 = _M0L4selfS640.$0;
    int32_t _M0L5startS2346 = _M0L4selfS640.$1;
    int32_t _M0L3valS2347 = _M0L1iS638->$0;
    int32_t _M0L6_2atmpS2345 = _M0L5startS2346 + _M0L3valS2347;
    moonbit_string_t _M0L4elemS641 =
      (moonbit_string_t)_M0L3bufS2344[_M0L6_2atmpS2345];
    int32_t _M0L3valS2343 = _M0L1iS638->$0;
    int32_t _M0L6_2atmpS2342 = _M0L3valS2343 + 1;
    _M0L1iS638->$0 = _M0L6_2atmpS2342;
    moonbit_incref(_M0L4elemS641);
    return _M0L4elemS641;
  } else {
    return 0;
  }
}

moonbit_string_t _M0IPC16uint646UInt64PB4Show10to__string(
  uint64_t _M0L4selfS637
) {
  #line 50 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 51 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC16uint646UInt6418to__string_2einner(_M0L4selfS637, 10);
}

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t _M0L4selfS636) {
  #line 35 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 36 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC13int3Int18to__string_2einner(_M0L4selfS636, 10);
}

moonbit_string_t _M0IPC14bool4BoolPB4Show10to__string(int32_t _M0L4selfS635) {
  #line 26 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L4selfS635) {
    return (moonbit_string_t)moonbit_string_literal_52.data;
  } else {
    return (moonbit_string_t)moonbit_string_literal_53.data;
  }
}

uint64_t _M0MPC14uint4UInt10to__uint64(uint32_t _M0L4selfS634) {
  #line 2494 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  return (uint64_t)_M0L4selfS634;
}

struct _M0TPC16string10StringView _M0MPC16string6String4trim(
  moonbit_string_t _M0L4selfS633,
  void* _M0L11chars_2eoptS631
) {
  struct _M0TPC16string10StringView _M0L5charsS630;
  struct _M0TPC16string10StringView _result_3994;
  switch (Moonbit_object_tag(_M0L11chars_2eoptS631)) {
    case 1: {
      struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS632 =
        (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L11chars_2eoptS631;
      struct _M0TPC16string10StringView _M0L8_2afieldS3620 =
        _M0L7_2aSomeS632->$0;
      moonbit_incref(_M0L8_2afieldS3620.$0);
      _M0L5charsS630 = _M0L8_2afieldS3620;
      break;
    }
    default: {
      int32_t _M0L6_2atmpS2335 =
        Moonbit_array_length(_M0MPC16string6String4trimN7_2abindS7043);
      moonbit_incref(_M0MPC16string6String4trimN7_2abindS7043);
      _M0L5charsS630
      = (struct _M0TPC16string10StringView){
        .$0 = _M0MPC16string6String4trimN7_2abindS7043,
          .$1 = 0,
          .$2 = _M0L6_2atmpS2335
      };
      break;
    }
  }
  _result_3994
  = _M0MPC16string6String12trim_2einner(_M0L4selfS633, _M0L5charsS630);
  moonbit_decref(_M0L5charsS630.$0);
  return _result_3994;
}

struct _M0TPC16string10StringView _M0MPC16string6String12trim_2einner(
  moonbit_string_t _M0L4selfS628,
  struct _M0TPC16string10StringView _M0L5charsS629
) {
  int32_t _M0L6_2atmpS2334;
  struct _M0TPC16string10StringView _M0L6_2atmpS2333;
  struct _M0TPC16string10StringView _result_3995;
  #line 794 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2334 = Moonbit_array_length(_M0L4selfS628);
  moonbit_incref(_M0L4selfS628);
  _M0L6_2atmpS2333
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS628, .$1 = 0, .$2 = _M0L6_2atmpS2334
  };
  #line 799 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3995
  = _M0MPC16string10StringView12trim_2einner(_M0L6_2atmpS2333, _M0L5charsS629);
  moonbit_decref(_M0L6_2atmpS2333.$0);
  return _result_3995;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView12trim_2einner(
  struct _M0TPC16string10StringView _M0L4selfS626,
  struct _M0TPC16string10StringView _M0L5charsS627
) {
  struct _M0TPC16string10StringView _M0L6_2atmpS2332;
  struct _M0TPC16string10StringView _result_3996;
  #line 783 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 788 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2332
  = _M0MPC16string10StringView19trim__start_2einner(_M0L4selfS626, _M0L5charsS627);
  #line 788 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3996
  = _M0MPC16string10StringView17trim__end_2einner(_M0L6_2atmpS2332, _M0L5charsS627);
  moonbit_decref(_M0L6_2atmpS2332.$0);
  return _result_3996;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView17trim__end_2einner(
  struct _M0TPC16string10StringView _M0L4selfS625,
  struct _M0TPC16string10StringView _M0L5charsS624
) {
  struct _M0TPC16string10StringView _M0L1xS620;
  #line 734 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  moonbit_incref(_M0L4selfS625.$0);
  _M0L1xS620 = _M0L4selfS625;
  while (1) {
    int32_t _M0L3endS2315 = _M0L1xS620.$2;
    int32_t _M0L5startS2316 = _M0L1xS620.$1;
    int32_t _M0L6_2atmpS2314 = _M0L3endS2315 - _M0L5startS2316;
    if (_M0L6_2atmpS2314 == 0) {
      return _M0L1xS620;
    } else {
      moonbit_string_t _M0L3strS2325 = _M0L1xS620.$0;
      moonbit_string_t _M0L3strS2328 = _M0L1xS620.$0;
      int32_t _M0L5startS2329 = _M0L1xS620.$1;
      int32_t _M0L3endS2331 = _M0L1xS620.$2;
      int64_t _M0L6_2atmpS2330 = (int64_t)_M0L3endS2331;
      int64_t _M0L6_2atmpS2327;
      int32_t _M0L6_2atmpS2326;
      int32_t _M0L4_2acS622;
      moonbit_string_t _M0L3strS2317;
      int32_t _M0L5startS2318;
      moonbit_string_t _M0L3strS2321;
      int32_t _M0L5startS2322;
      int32_t _M0L3endS2324;
      int64_t _M0L6_2atmpS2323;
      int64_t _M0L6_2atmpS2320;
      int32_t _M0L6_2atmpS2319;
      struct _M0TPC16string10StringView _M0L4_2axS623;
      moonbit_incref(_M0L3strS2328);
      moonbit_incref(_M0L3strS2325);
      #line 739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L6_2atmpS2327
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2328, -1, _M0L5startS2329, _M0L6_2atmpS2330);
      moonbit_decref(_M0L3strS2328);
      _M0L6_2atmpS2326 = (int32_t)_M0L6_2atmpS2327;
      #line 739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L4_2acS622
      = _M0MPC16string6String16unsafe__char__at(_M0L3strS2325, _M0L6_2atmpS2326);
      moonbit_decref(_M0L3strS2325);
      _M0L3strS2317 = _M0L1xS620.$0;
      _M0L5startS2318 = _M0L1xS620.$1;
      _M0L3strS2321 = _M0L1xS620.$0;
      _M0L5startS2322 = _M0L1xS620.$1;
      _M0L3endS2324 = _M0L1xS620.$2;
      _M0L6_2atmpS2323 = (int64_t)_M0L3endS2324;
      moonbit_incref(_M0L3strS2321);
      moonbit_incref(_M0L3strS2317);
      #line 739 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L6_2atmpS2320
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2321, -1, _M0L5startS2322, _M0L6_2atmpS2323);
      moonbit_decref(_M0L3strS2321);
      _M0L6_2atmpS2319 = (int32_t)_M0L6_2atmpS2320;
      _M0L4_2axS623
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L3strS2317, .$1 = _M0L5startS2318, .$2 = _M0L6_2atmpS2319
      };
      #line 743 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      if (
        _M0MPC16string10StringView14contains__char(_M0L5charsS624, _M0L4_2acS622)
      ) {
        moonbit_decref(_M0L1xS620.$0);
        _M0L1xS620 = _M0L4_2axS623;
        continue;
      } else {
        moonbit_decref(_M0L4_2axS623.$0);
        return _M0L1xS620;
      }
    }
    break;
  }
}

struct _M0TPC16string10StringView _M0MPC16string10StringView19trim__start_2einner(
  struct _M0TPC16string10StringView _M0L4selfS619,
  struct _M0TPC16string10StringView _M0L5charsS618
) {
  struct _M0TPC16string10StringView _M0L1xS613;
  #line 686 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  moonbit_incref(_M0L4selfS619.$0);
  _M0L1xS613 = _M0L4selfS619;
  while (1) {
    int32_t _M0L3endS2298 = _M0L1xS613.$2;
    int32_t _M0L5startS2299 = _M0L1xS613.$1;
    int32_t _M0L6_2atmpS2297 = _M0L3endS2298 - _M0L5startS2299;
    if (_M0L6_2atmpS2297 == 0) {
      return _M0L1xS613;
    } else {
      moonbit_string_t _M0L3strS2307 = _M0L1xS613.$0;
      moonbit_string_t _M0L3strS2310 = _M0L1xS613.$0;
      int32_t _M0L5startS2311 = _M0L1xS613.$1;
      int32_t _M0L3endS2313 = _M0L1xS613.$2;
      int64_t _M0L6_2atmpS2312 = (int64_t)_M0L3endS2313;
      int64_t _M0L6_2atmpS2309;
      int32_t _M0L6_2atmpS2308;
      int32_t _M0L4_2acS615;
      moonbit_string_t _M0L3strS2300;
      moonbit_string_t _M0L3strS2303;
      int32_t _M0L5startS2304;
      int32_t _M0L3endS2306;
      int64_t _M0L6_2atmpS2305;
      int64_t _M0L7_2abindS1574;
      int32_t _M0L6_2atmpS2301;
      int32_t _M0L3endS2302;
      struct _M0TPC16string10StringView _M0L4_2axS616;
      moonbit_incref(_M0L3strS2310);
      moonbit_incref(_M0L3strS2307);
      #line 691 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L6_2atmpS2309
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2310, 0, _M0L5startS2311, _M0L6_2atmpS2312);
      moonbit_decref(_M0L3strS2310);
      _M0L6_2atmpS2308 = (int32_t)_M0L6_2atmpS2309;
      #line 691 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L4_2acS615
      = _M0MPC16string6String16unsafe__char__at(_M0L3strS2307, _M0L6_2atmpS2308);
      moonbit_decref(_M0L3strS2307);
      _M0L3strS2300 = _M0L1xS613.$0;
      _M0L3strS2303 = _M0L1xS613.$0;
      _M0L5startS2304 = _M0L1xS613.$1;
      _M0L3endS2306 = _M0L1xS613.$2;
      _M0L6_2atmpS2305 = (int64_t)_M0L3endS2306;
      moonbit_incref(_M0L3strS2303);
      moonbit_incref(_M0L3strS2300);
      #line 691 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS1574
      = _M0MPC16string6String29offset__of__nth__char_2einner(_M0L3strS2303, 1, _M0L5startS2304, _M0L6_2atmpS2305);
      moonbit_decref(_M0L3strS2303);
      if (_M0L7_2abindS1574 == 4294967296ll) {
        _M0L6_2atmpS2301 = _M0L1xS613.$2;
      } else {
        int64_t _M0L7_2aSomeS617 = _M0L7_2abindS1574;
        _M0L6_2atmpS2301 = (int32_t)_M0L7_2aSomeS617;
      }
      _M0L3endS2302 = _M0L1xS613.$2;
      _M0L4_2axS616
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L3strS2300, .$1 = _M0L6_2atmpS2301, .$2 = _M0L3endS2302
      };
      #line 695 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      if (
        _M0MPC16string10StringView14contains__char(_M0L5charsS618, _M0L4_2acS615)
      ) {
        moonbit_decref(_M0L1xS613.$0);
        _M0L1xS613 = _M0L4_2axS616;
        continue;
      } else {
        moonbit_decref(_M0L4_2axS616.$0);
        return _M0L1xS613;
      }
    }
    break;
  }
}

int32_t _M0MPC16string10StringView14contains__char(
  struct _M0TPC16string10StringView _M0L4selfS604,
  int32_t _M0L1cS606
) {
  int32_t _M0L3endS2295;
  int32_t _M0L5startS2296;
  int32_t _M0L3lenS603;
  #line 622 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2295 = _M0L4selfS604.$2;
  _M0L5startS2296 = _M0L4selfS604.$1;
  _M0L3lenS603 = _M0L3endS2295 - _M0L5startS2296;
  if (_M0L3lenS603 > 0) {
    int32_t _M0L1cS605 = _M0L1cS606;
    int32_t _if__result_3999;
    if (_M0L1cS605 >= 0) {
      _if__result_3999 = _M0L1cS605 <= 65535;
    } else {
      _if__result_3999 = 0;
    }
    if (_if__result_3999) {
      int32_t _M0L6_2atmpS2279 = (uint16_t)_M0L1cS605;
      #line 629 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0MPC16string10StringView20contains__code__unit(_M0L4selfS604, _M0L6_2atmpS2279);
    } else if (_M0L1cS605 < 0) {
      return 0;
    } else if (_M0L3lenS603 >= 2) {
      int32_t _M0L3adjS607 = _M0L1cS605 - 65536;
      int32_t _M0L6_2atmpS2294 = _M0L3adjS607 >> 10;
      int32_t _M0L4highS608 = 55296 + _M0L6_2atmpS2294;
      if (_M0L4highS608 <= 65535) {
        int32_t _M0L4highS609 = (uint16_t)_M0L4highS608;
        int32_t _M0L6_2atmpS2293 = _M0L3adjS607 & 1023;
        int32_t _M0L6_2atmpS2292 = 56320 + _M0L6_2atmpS2293;
        int32_t _M0L3lowS610 = (uint16_t)_M0L6_2atmpS2292;
        int32_t _M0L1iS611 = 0;
        while (1) {
          int32_t _M0L6_2atmpS2280 = _M0L3lenS603 - 1;
          if (_M0L1iS611 < _M0L6_2atmpS2280) {
            moonbit_string_t _M0L3strS2282 = _M0L4selfS604.$0;
            int32_t _M0L5startS2284 = _M0L4selfS604.$1;
            int32_t _M0L6_2atmpS2283 = _M0L5startS2284 + _M0L1iS611;
            int32_t _M0L6_2atmpS2281 = _M0L3strS2282[_M0L6_2atmpS2283];
            int32_t _M0L6_2atmpS2291;
            if (_M0L6_2atmpS2281 == _M0L4highS609) {
              moonbit_string_t _M0L3strS2286 = _M0L4selfS604.$0;
              int32_t _M0L5startS2288 = _M0L4selfS604.$1;
              int32_t _M0L6_2atmpS2289 = _M0L1iS611 + 1;
              int32_t _M0L6_2atmpS2287 = _M0L5startS2288 + _M0L6_2atmpS2289;
              int32_t _M0L6_2atmpS2285 = _M0L3strS2286[_M0L6_2atmpS2287];
              int32_t _M0L6_2atmpS2290;
              if (_M0L6_2atmpS2285 == _M0L3lowS610) {
                return 1;
              }
              _M0L6_2atmpS2290 = _M0L1iS611 + 2;
              _M0L1iS611 = _M0L6_2atmpS2290;
              continue;
            }
            _M0L6_2atmpS2291 = _M0L1iS611 + 1;
            _M0L1iS611 = _M0L6_2atmpS2291;
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

int32_t _M0MPC16string10StringView20contains__code__unit(
  struct _M0TPC16string10StringView _M0L4selfS599,
  int32_t _M0L4codeS601
) {
  int32_t _M0L3endS2277;
  int32_t _M0L5startS2278;
  int32_t _M0L7_2abindS598;
  int32_t _M0L1iS600;
  #line 524 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2277 = _M0L4selfS599.$2;
  _M0L5startS2278 = _M0L4selfS599.$1;
  _M0L7_2abindS598 = _M0L3endS2277 - _M0L5startS2278;
  _M0L1iS600 = 0;
  while (1) {
    if (_M0L1iS600 < _M0L7_2abindS598) {
      moonbit_string_t _M0L3strS2273 = _M0L4selfS599.$0;
      int32_t _M0L5startS2275 = _M0L4selfS599.$1;
      int32_t _M0L6_2atmpS2274 = _M0L5startS2275 + _M0L1iS600;
      int32_t _M0L6_2atmpS2272 = _M0L3strS2273[_M0L6_2atmpS2274];
      int32_t _M0L6_2atmpS2276;
      if (_M0L6_2atmpS2272 == _M0L4codeS601) {
        return 1;
      }
      _M0L6_2atmpS2276 = _M0L1iS600 + 1;
      _M0L1iS600 = _M0L6_2atmpS2276;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS586,
  moonbit_string_t _M0L5valueS588
) {
  int32_t _M0L3lenS2252;
  moonbit_string_t* _M0L6_2atmpS2254;
  int32_t _M0L6_2atmpS2253;
  int32_t _M0L6lengthS587;
  moonbit_string_t* _M0L3bufS2255;
  moonbit_string_t _M0L6_2aoldS3632;
  int32_t _M0L6_2atmpS2256;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2252 = _M0L4selfS586->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2254 = _M0MPC15array5Array6bufferGsE(_M0L4selfS586);
  _M0L6_2atmpS2253 = Moonbit_array_length(_M0L6_2atmpS2254);
  moonbit_decref(_M0L6_2atmpS2254);
  if (_M0L3lenS2252 == _M0L6_2atmpS2253) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGsE(_M0L4selfS586);
  }
  _M0L6lengthS587 = _M0L4selfS586->$1;
  _M0L3bufS2255 = _M0L4selfS586->$0;
  _M0L6_2aoldS3632 = (moonbit_string_t)_M0L3bufS2255[_M0L6lengthS587];
  moonbit_incref(_M0L5valueS588);
  moonbit_decref(_M0L6_2aoldS3632);
  _M0L3bufS2255[_M0L6lengthS587] = _M0L5valueS588;
  _M0L6_2atmpS2256 = _M0L6lengthS587 + 1;
  _M0L4selfS586->$1 = _M0L6_2atmpS2256;
  return 0;
}

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS589,
  struct _M0TUsiE* _M0L5valueS591
) {
  int32_t _M0L3lenS2257;
  struct _M0TUsiE** _M0L6_2atmpS2259;
  int32_t _M0L6_2atmpS2258;
  int32_t _M0L6lengthS590;
  struct _M0TUsiE** _M0L3bufS2260;
  struct _M0TUsiE* _M0L6_2aoldS3634;
  int32_t _M0L6_2atmpS2261;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2257 = _M0L4selfS589->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2259 = _M0MPC15array5Array6bufferGUsiEE(_M0L4selfS589);
  _M0L6_2atmpS2258 = Moonbit_array_length(_M0L6_2atmpS2259);
  moonbit_decref(_M0L6_2atmpS2259);
  if (_M0L3lenS2257 == _M0L6_2atmpS2258) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGUsiEE(_M0L4selfS589);
  }
  _M0L6lengthS590 = _M0L4selfS589->$1;
  _M0L3bufS2260 = _M0L4selfS589->$0;
  _M0L6_2aoldS3634 = (struct _M0TUsiE*)_M0L3bufS2260[_M0L6lengthS590];
  moonbit_incref(_M0L5valueS591);
  if (_M0L6_2aoldS3634) {
    moonbit_decref(_M0L6_2aoldS3634);
  }
  _M0L3bufS2260[_M0L6lengthS590] = _M0L5valueS591;
  _M0L6_2atmpS2261 = _M0L6lengthS590 + 1;
  _M0L4selfS589->$1 = _M0L6_2atmpS2261;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS592,
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L5valueS594
) {
  int32_t _M0L3lenS2262;
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS2264;
  int32_t _M0L6_2atmpS2263;
  int32_t _M0L6lengthS593;
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3bufS2265;
  struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2aoldS3636;
  int32_t _M0L6_2atmpS2266;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2262 = _M0L4selfS592->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2264
  = _M0MPC15array5Array6bufferGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS592);
  _M0L6_2atmpS2263 = Moonbit_array_length(_M0L6_2atmpS2264);
  moonbit_decref(_M0L6_2atmpS2264);
  if (_M0L3lenS2262 == _M0L6_2atmpS2263) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS592);
  }
  _M0L6lengthS593 = _M0L4selfS592->$1;
  _M0L3bufS2265 = _M0L4selfS592->$0;
  _M0L6_2aoldS3636
  = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3bufS2265[_M0L6lengthS593];
  moonbit_incref(_M0L5valueS594);
  if (_M0L6_2aoldS3636) {
    moonbit_decref(_M0L6_2aoldS3636);
  }
  _M0L3bufS2265[_M0L6lengthS593] = _M0L5valueS594;
  _M0L6_2atmpS2266 = _M0L6lengthS593 + 1;
  _M0L4selfS592->$1 = _M0L6_2atmpS2266;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L4selfS595,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L5valueS597
) {
  int32_t _M0L3lenS2267;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L6_2atmpS2269;
  int32_t _M0L6_2atmpS2268;
  int32_t _M0L6lengthS596;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3bufS2270;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2aoldS3638;
  int32_t _M0L6_2atmpS2271;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS2267 = _M0L4selfS595->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS2269
  = _M0MPC15array5Array6bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L4selfS595);
  _M0L6_2atmpS2268 = Moonbit_array_length(_M0L6_2atmpS2269);
  moonbit_decref(_M0L6_2atmpS2269);
  if (_M0L3lenS2267 == _M0L6_2atmpS2268) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L4selfS595);
  }
  _M0L6lengthS596 = _M0L4selfS595->$1;
  _M0L3bufS2270 = _M0L4selfS595->$0;
  _M0L6_2aoldS3638
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)_M0L3bufS2270[
      _M0L6lengthS596
    ];
  moonbit_incref(_M0L5valueS597);
  if (_M0L6_2aoldS3638) {
    moonbit_decref(_M0L6_2aoldS3638);
  }
  _M0L3bufS2270[_M0L6lengthS596] = _M0L5valueS597;
  _M0L6_2atmpS2271 = _M0L6lengthS596 + 1;
  _M0L4selfS595->$1 = _M0L6_2atmpS2271;
  return 0;
}

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE* _M0L4selfS575) {
  int32_t _M0L8old__capS574;
  int32_t _M0L8new__capS576;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS574 = _M0L4selfS575->$1;
  if (_M0L8old__capS574 == 0) {
    _M0L8new__capS576 = 8;
  } else {
    _M0L8new__capS576 = _M0L8old__capS574 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGsE(_M0L4selfS575, _M0L8new__capS576);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS578
) {
  int32_t _M0L8old__capS577;
  int32_t _M0L8new__capS579;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS577 = _M0L4selfS578->$1;
  if (_M0L8old__capS577 == 0) {
    _M0L8new__capS579 = 8;
  } else {
    _M0L8new__capS579 = _M0L8old__capS577 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGUsiEE(_M0L4selfS578, _M0L8new__capS579);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS581
) {
  int32_t _M0L8old__capS580;
  int32_t _M0L8new__capS582;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS580 = _M0L4selfS581->$1;
  if (_M0L8old__capS580 == 0) {
    _M0L8new__capS582 = 8;
  } else {
    _M0L8new__capS582 = _M0L8old__capS580 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq8TreeNodeE(_M0L4selfS581, _M0L8new__capS582);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L4selfS584
) {
  int32_t _M0L8old__capS583;
  int32_t _M0L8new__capS585;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS583 = _M0L4selfS584->$1;
  if (_M0L8old__capS583 == 0) {
    _M0L8new__capS585 = 8;
  } else {
    _M0L8new__capS585 = _M0L8old__capS583 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L4selfS584, _M0L8new__capS585);
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS551,
  int32_t _M0L13new__capacityS554
) {
  moonbit_string_t* _M0L8old__bufS550;
  int32_t _M0L8old__capS552;
  int32_t _M0L9copy__lenS553;
  moonbit_string_t* _M0L8new__bufS555;
  moonbit_string_t* _M0L6_2aoldS3640;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS550 = _M0L4selfS551->$0;
  _M0L8old__capS552 = Moonbit_array_length(_M0L8old__bufS550);
  if (_M0L8old__capS552 < _M0L13new__capacityS554) {
    _M0L9copy__lenS553 = _M0L8old__capS552;
  } else {
    _M0L9copy__lenS553 = _M0L13new__capacityS554;
  }
  moonbit_incref(_M0L8old__bufS550);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS555
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(_M0L8old__bufS550, _M0L13new__capacityS554, _M0L9copy__lenS553, 0, 0);
  moonbit_decref(_M0L8old__bufS550);
  _M0L6_2aoldS3640 = _M0L4selfS551->$0;
  moonbit_decref(_M0L6_2aoldS3640);
  _M0L4selfS551->$0 = _M0L8new__bufS555;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS557,
  int32_t _M0L13new__capacityS560
) {
  struct _M0TUsiE** _M0L8old__bufS556;
  int32_t _M0L8old__capS558;
  int32_t _M0L9copy__lenS559;
  struct _M0TUsiE** _M0L8new__bufS561;
  struct _M0TUsiE** _M0L6_2aoldS3642;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS556 = _M0L4selfS557->$0;
  _M0L8old__capS558 = Moonbit_array_length(_M0L8old__bufS556);
  if (_M0L8old__capS558 < _M0L13new__capacityS560) {
    _M0L9copy__lenS559 = _M0L8old__capS558;
  } else {
    _M0L9copy__lenS559 = _M0L13new__capacityS560;
  }
  moonbit_incref(_M0L8old__bufS556);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS561
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(_M0L8old__bufS556, _M0L13new__capacityS560, _M0L9copy__lenS559, 0, 0);
  moonbit_decref(_M0L8old__bufS556);
  _M0L6_2aoldS3642 = _M0L4selfS557->$0;
  moonbit_decref(_M0L6_2aoldS3642);
  _M0L4selfS557->$0 = _M0L8new__bufS561;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS563,
  int32_t _M0L13new__capacityS566
) {
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L8old__bufS562;
  int32_t _M0L8old__capS564;
  int32_t _M0L9copy__lenS565;
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L8new__bufS567;
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2aoldS3644;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS562 = _M0L4selfS563->$0;
  _M0L8old__capS564 = Moonbit_array_length(_M0L8old__bufS562);
  if (_M0L8old__capS564 < _M0L13new__capacityS566) {
    _M0L9copy__lenS565 = _M0L8old__capS564;
  } else {
    _M0L9copy__lenS565 = _M0L13new__capacityS566;
  }
  moonbit_incref(_M0L8old__bufS562);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS567
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq8TreeNodeE(_M0L8old__bufS562, _M0L13new__capacityS566, _M0L9copy__lenS565, 0, 0);
  moonbit_decref(_M0L8old__bufS562);
  _M0L6_2aoldS3644 = _M0L4selfS563->$0;
  moonbit_decref(_M0L6_2aoldS3644);
  _M0L4selfS563->$0 = _M0L8new__bufS567;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L4selfS569,
  int32_t _M0L13new__capacityS572
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L8old__bufS568;
  int32_t _M0L8old__capS570;
  int32_t _M0L9copy__lenS571;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L8new__bufS573;
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L6_2aoldS3646;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS568 = _M0L4selfS569->$0;
  _M0L8old__capS570 = Moonbit_array_length(_M0L8old__bufS568);
  if (_M0L8old__capS570 < _M0L13new__capacityS572) {
    _M0L9copy__lenS571 = _M0L8old__capS570;
  } else {
    _M0L9copy__lenS571 = _M0L13new__capacityS572;
  }
  moonbit_incref(_M0L8old__bufS568);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS573
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L8old__bufS568, _M0L13new__capacityS572, _M0L9copy__lenS571, 0, 0);
  moonbit_decref(_M0L8old__bufS568);
  _M0L6_2aoldS3646 = _M0L4selfS569->$0;
  moonbit_decref(_M0L6_2aoldS3646);
  _M0L4selfS569->$0 = _M0L8new__bufS573;
  return 0;
}

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS547
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS547->$1;
}

int32_t _M0MPC15array5Array6lengthGsE(struct _M0TPB5ArrayGsE* _M0L4selfS548) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS548->$1;
}

int32_t _M0MPC15array5Array6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L4selfS549
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS549->$1;
}

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(
  int32_t _M0L8capacityS544
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS544 == 0) {
    moonbit_string_t* _M0L6_2atmpS2246 =
      (moonbit_string_t*)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGsE* _block_4002 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_4002)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
    _block_4002->$0 = _M0L6_2atmpS2246;
    _block_4002->$1 = 0;
    return _block_4002;
  } else {
    moonbit_string_t* _M0L6_2atmpS2247 =
      (moonbit_string_t*)moonbit_make_ref_array(_M0L8capacityS544, (moonbit_string_t)moonbit_string_literal_1.data);
    struct _M0TPB5ArrayGsE* _block_4003 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_4003)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 0);
    _block_4003->$0 = _M0L6_2atmpS2247;
    _block_4003->$1 = 0;
    return _block_4003;
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq8TreeNodeE(
  int32_t _M0L8capacityS545
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS545 == 0) {
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS2248 =
      (struct _M0TP26biolab8bio__seq8TreeNode**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _block_4004 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE));
    Moonbit_object_header(_block_4004)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 42, 0);
    _block_4004->$0 = _M0L6_2atmpS2248;
    _block_4004->$1 = 0;
    return _block_4004;
  } else {
    struct _M0TP26biolab8bio__seq8TreeNode** _M0L6_2atmpS2249 =
      (struct _M0TP26biolab8bio__seq8TreeNode**)moonbit_make_ref_array(_M0L8capacityS545, 0);
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _block_4005 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE));
    Moonbit_object_header(_block_4005)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 42, 0);
    _block_4005->$0 = _M0L6_2atmpS2249;
    _block_4005->$1 = 0;
    return _block_4005;
  }
}

struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0MPC15array5Array11new_2einnerGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  int32_t _M0L8capacityS546
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS546 == 0) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L6_2atmpS2250 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _block_4006 =
      (struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE));
    Moonbit_object_header(_block_4006)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 81, 0);
    _block_4006->$0 = _M0L6_2atmpS2250;
    _block_4006->$1 = 0;
    return _block_4006;
  } else {
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L6_2atmpS2251 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**)moonbit_make_ref_array(_M0L8capacityS546, 0);
    struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _block_4007 =
      (struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE));
    Moonbit_object_header(_block_4007)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 81, 0);
    _block_4007->$0 = _M0L6_2atmpS2251;
    _block_4007->$1 = 0;
    return _block_4007;
  }
}

int32_t _M0MPC16string6String11has__suffix(
  moonbit_string_t _M0L4selfS542,
  struct _M0TPC16string10StringView _M0L3strS543
) {
  int32_t _M0L6_2atmpS2245;
  struct _M0TPC16string10StringView _M0L6_2atmpS2244;
  int32_t _result_4008;
  #line 269 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS2245 = Moonbit_array_length(_M0L4selfS542);
  moonbit_incref(_M0L4selfS542);
  _M0L6_2atmpS2244
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS542, .$1 = 0, .$2 = _M0L6_2atmpS2245
  };
  #line 271 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_4008
  = _M0MPC16string10StringView11has__suffix(_M0L6_2atmpS2244, _M0L3strS543);
  moonbit_decref(_M0L6_2atmpS2244.$0);
  return _result_4008;
}

int32_t _M0MPC16string10StringView11has__suffix(
  struct _M0TPC16string10StringView _M0L4selfS538,
  struct _M0TPC16string10StringView _M0L3strS539
) {
  int64_t _M0L7_2abindS537;
  #line 262 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 264 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS537
  = _M0MPC16string10StringView9rev__find(_M0L4selfS538, _M0L3strS539);
  if (_M0L7_2abindS537 == 4294967296ll) {
    return 0;
  } else {
    int64_t _M0L7_2aSomeS540 = _M0L7_2abindS537;
    int32_t _M0L4_2aiS541 = (int32_t)_M0L7_2aSomeS540;
    int32_t _M0L3endS2242 = _M0L4selfS538.$2;
    int32_t _M0L5startS2243 = _M0L4selfS538.$1;
    int32_t _M0L6_2atmpS2238 = _M0L3endS2242 - _M0L5startS2243;
    int32_t _M0L3endS2240 = _M0L3strS539.$2;
    int32_t _M0L5startS2241 = _M0L3strS539.$1;
    int32_t _M0L6_2atmpS2239 = _M0L3endS2240 - _M0L5startS2241;
    int32_t _M0L6_2atmpS2237 = _M0L6_2atmpS2238 - _M0L6_2atmpS2239;
    return _M0L4_2aiS541 == _M0L6_2atmpS2237;
  }
}

int64_t _M0MPC16string10StringView9rev__find(
  struct _M0TPC16string10StringView _M0L4selfS536,
  struct _M0TPC16string10StringView _M0L3strS535
) {
  int32_t _M0L3endS2235;
  int32_t _M0L5startS2236;
  int32_t _M0L6_2atmpS2234;
  #line 165 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2235 = _M0L3strS535.$2;
  _M0L5startS2236 = _M0L3strS535.$1;
  _M0L6_2atmpS2234 = _M0L3endS2235 - _M0L5startS2236;
  if (_M0L6_2atmpS2234 <= 4) {
    #line 167 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB23brute__force__rev__find(_M0L4selfS536, _M0L3strS535);
  } else {
    #line 169 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB33boyer__moore__horspool__rev__find(_M0L4selfS536, _M0L3strS535);
  }
}

int64_t _M0FPB23brute__force__rev__find(
  struct _M0TPC16string10StringView _M0L8haystackS525,
  struct _M0TPC16string10StringView _M0L6needleS527
) {
  int32_t _M0L3endS2232;
  int32_t _M0L5startS2233;
  int32_t _M0L13haystack__lenS524;
  int32_t _M0L3endS2230;
  int32_t _M0L5startS2231;
  int32_t _M0L11needle__lenS526;
  #line 178 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2232 = _M0L8haystackS525.$2;
  _M0L5startS2233 = _M0L8haystackS525.$1;
  _M0L13haystack__lenS524 = _M0L3endS2232 - _M0L5startS2233;
  _M0L3endS2230 = _M0L6needleS527.$2;
  _M0L5startS2231 = _M0L6needleS527.$1;
  _M0L11needle__lenS526 = _M0L3endS2230 - _M0L5startS2231;
  if (_M0L11needle__lenS526 > 0) {
    if (_M0L13haystack__lenS524 >= _M0L11needle__lenS526) {
      moonbit_string_t _M0L3strS2228 = _M0L6needleS527.$0;
      int32_t _M0L5startS2229 = _M0L6needleS527.$1;
      int32_t _M0L13needle__firstS528 = _M0L3strS2228[_M0L5startS2229];
      int32_t _M0L7_2abindS529 =
        _M0L13haystack__lenS524 - _M0L11needle__lenS526;
      int32_t _M0L1iS530 = _M0L7_2abindS529;
      while (1) {
        if (_M0L1iS530 >= 0) {
          moonbit_string_t _M0L3strS2215 = _M0L8haystackS525.$0;
          int32_t _M0L5startS2217 = _M0L8haystackS525.$1;
          int32_t _M0L6_2atmpS2216 = _M0L5startS2217 + _M0L1iS530;
          int32_t _M0L6_2atmpS2214 = _M0L3strS2215[_M0L6_2atmpS2216];
          int32_t _M0L1jS533;
          int32_t _M0L6_2atmpS2213;
          if (_M0L6_2atmpS2214 != _M0L13needle__firstS528) {
            goto join_531;
          }
          _M0L1jS533 = 1;
          while (1) {
            if (_M0L1jS533 < _M0L11needle__lenS526) {
              moonbit_string_t _M0L3strS2223 = _M0L8haystackS525.$0;
              int32_t _M0L5startS2225 = _M0L8haystackS525.$1;
              int32_t _M0L6_2atmpS2226 = _M0L1iS530 + _M0L1jS533;
              int32_t _M0L6_2atmpS2224 = _M0L5startS2225 + _M0L6_2atmpS2226;
              int32_t _M0L6_2atmpS2218 = _M0L3strS2223[_M0L6_2atmpS2224];
              moonbit_string_t _M0L3strS2220 = _M0L6needleS527.$0;
              int32_t _M0L5startS2222 = _M0L6needleS527.$1;
              int32_t _M0L6_2atmpS2221 = _M0L5startS2222 + _M0L1jS533;
              int32_t _M0L6_2atmpS2219 = _M0L3strS2220[_M0L6_2atmpS2221];
              int32_t _M0L6_2atmpS2227;
              if (_M0L6_2atmpS2218 != _M0L6_2atmpS2219) {
                break;
              }
              _M0L6_2atmpS2227 = _M0L1jS533 + 1;
              _M0L1jS533 = _M0L6_2atmpS2227;
              continue;
            } else {
              return (int64_t)_M0L1iS530;
            }
            break;
          }
          goto join_531;
          goto joinlet_4010;
          join_531:;
          _M0L6_2atmpS2213 = _M0L1iS530 - 1;
          _M0L1iS530 = _M0L6_2atmpS2213;
          continue;
          joinlet_4010:;
        }
        break;
      }
      return 4294967296ll;
    } else {
      return 4294967296ll;
    }
  } else {
    return (int64_t)_M0L13haystack__lenS524;
  }
}

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(
  moonbit_string_t _M0L4selfS523
) {
  #line 222 "/home/developer/.moon/lib/core/builtin/show.mbt"
  moonbit_incref(_M0L4selfS523);
  return _M0L4selfS523;
}

int64_t _M0FPB33boyer__moore__horspool__rev__find(
  struct _M0TPC16string10StringView _M0L8haystackS513,
  struct _M0TPC16string10StringView _M0L6needleS515
) {
  int32_t _M0L3endS2211;
  int32_t _M0L5startS2212;
  int32_t _M0L13haystack__lenS512;
  int32_t _M0L3endS2209;
  int32_t _M0L5startS2210;
  int32_t _M0L11needle__lenS514;
  #line 203 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS2211 = _M0L8haystackS513.$2;
  _M0L5startS2212 = _M0L8haystackS513.$1;
  _M0L13haystack__lenS512 = _M0L3endS2211 - _M0L5startS2212;
  _M0L3endS2209 = _M0L6needleS515.$2;
  _M0L5startS2210 = _M0L6needleS515.$1;
  _M0L11needle__lenS514 = _M0L3endS2209 - _M0L5startS2210;
  if (_M0L11needle__lenS514 > 0) {
    if (_M0L13haystack__lenS512 >= _M0L11needle__lenS514) {
      int32_t* _M0L11skip__tableS516 =
        (int32_t*)moonbit_make_int32_array(256, _M0L11needle__lenS514);
      int32_t _M0L6_2atmpS2189 = _M0L11needle__lenS514 - 1;
      int32_t _M0L1iS517 = _M0L6_2atmpS2189;
      int32_t _M0L6_2atmpS2208;
      int32_t _M0L1iS519;
      while (1) {
        if (_M0L1iS517 >= 1) {
          moonbit_string_t _M0L3strS2185 = _M0L6needleS515.$0;
          int32_t _M0L5startS2187 = _M0L6needleS515.$1;
          int32_t _M0L6_2atmpS2186 = _M0L5startS2187 + _M0L1iS517;
          int32_t _M0L6_2atmpS2184 = _M0L3strS2185[_M0L6_2atmpS2186];
          int32_t _M0L6_2atmpS2183 = (int32_t)_M0L6_2atmpS2184;
          int32_t _M0L6_2atmpS2182 = _M0L6_2atmpS2183 & 255;
          int32_t _M0L6_2atmpS2188;
          if (
            _M0L6_2atmpS2182 < 0
            || _M0L6_2atmpS2182
               >= Moonbit_array_length(_M0L11skip__tableS516)
          ) {
            #line 213 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L11skip__tableS516[_M0L6_2atmpS2182] = _M0L1iS517;
          _M0L6_2atmpS2188 = _M0L1iS517 - 1;
          _M0L1iS517 = _M0L6_2atmpS2188;
          continue;
        }
        break;
      }
      _M0L6_2atmpS2208 = _M0L13haystack__lenS512 - _M0L11needle__lenS514;
      _M0L1iS519 = _M0L6_2atmpS2208;
      while (1) {
        if (_M0L1iS519 >= 0) {
          int32_t _M0L1jS520 = 0;
          moonbit_string_t _M0L3strS2205;
          int32_t _M0L5startS2207;
          int32_t _M0L6_2atmpS2206;
          int32_t _M0L6_2atmpS2204;
          int32_t _M0L6_2atmpS2203;
          int32_t _M0L6_2atmpS2202;
          int32_t _M0L6_2atmpS2201;
          int32_t _M0L6_2atmpS2200;
          while (1) {
            if (_M0L1jS520 < _M0L11needle__lenS514) {
              moonbit_string_t _M0L3strS2195 = _M0L8haystackS513.$0;
              int32_t _M0L5startS2197 = _M0L8haystackS513.$1;
              int32_t _M0L6_2atmpS2198 = _M0L1iS519 + _M0L1jS520;
              int32_t _M0L6_2atmpS2196 = _M0L5startS2197 + _M0L6_2atmpS2198;
              int32_t _M0L6_2atmpS2190 = _M0L3strS2195[_M0L6_2atmpS2196];
              moonbit_string_t _M0L3strS2192 = _M0L6needleS515.$0;
              int32_t _M0L5startS2194 = _M0L6needleS515.$1;
              int32_t _M0L6_2atmpS2193 = _M0L5startS2194 + _M0L1jS520;
              int32_t _M0L6_2atmpS2191 = _M0L3strS2192[_M0L6_2atmpS2193];
              int32_t _M0L6_2atmpS2199;
              if (_M0L6_2atmpS2190 != _M0L6_2atmpS2191) {
                break;
              }
              _M0L6_2atmpS2199 = _M0L1jS520 + 1;
              _M0L1jS520 = _M0L6_2atmpS2199;
              continue;
            } else {
              moonbit_decref(_M0L11skip__tableS516);
              return (int64_t)_M0L1iS519;
            }
            break;
          }
          _M0L3strS2205 = _M0L8haystackS513.$0;
          _M0L5startS2207 = _M0L8haystackS513.$1;
          _M0L6_2atmpS2206 = _M0L5startS2207 + _M0L1iS519;
          _M0L6_2atmpS2204 = _M0L3strS2205[_M0L6_2atmpS2206];
          _M0L6_2atmpS2203 = (int32_t)_M0L6_2atmpS2204;
          _M0L6_2atmpS2202 = _M0L6_2atmpS2203 & 255;
          if (
            _M0L6_2atmpS2202 < 0
            || _M0L6_2atmpS2202
               >= Moonbit_array_length(_M0L11skip__tableS516)
          ) {
            #line 217 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2201 = (int32_t)_M0L11skip__tableS516[_M0L6_2atmpS2202];
          _M0L6_2atmpS2200 = _M0L1iS519 - _M0L6_2atmpS2201;
          _M0L1iS519 = _M0L6_2atmpS2200;
          continue;
        } else {
          moonbit_decref(_M0L11skip__tableS516);
        }
        break;
      }
      return 4294967296ll;
    } else {
      return 4294967296ll;
    }
  } else {
    return (int64_t)_M0L13haystack__lenS512;
  }
}

struct moonbit_result_0 _M0FPB12assert__true(
  int32_t _M0L1xS505,
  void* _M0L3msgS507,
  moonbit_string_t _M0L3locS511
) {
  #line 123 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  if (!_M0L1xS505) {
    struct _M0TPC16string10StringView _M0L9fail__msgS506;
    struct moonbit_result_0 _result_4015;
    switch (Moonbit_object_tag(_M0L3msgS507)) {
      case 1: {
        struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS508 =
          (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L3msgS507;
        struct _M0TPC16string10StringView _M0L8_2afieldS3656 =
          _M0L7_2aSomeS508->$0;
        moonbit_incref(_M0L8_2afieldS3656.$0);
        _M0L9fail__msgS506 = _M0L8_2afieldS3656;
        break;
      }
      default: {
        struct _M0TPB13StringBuilder* _M0L18_2astring__builderS509;
        moonbit_string_t _M0L7_2abindS510;
        int32_t _M0L6_2atmpS2180;
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0L18_2astring__builderS509
        = _M0MPB13StringBuilder21StringBuilder_2einner(14);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS509, (moonbit_string_t)moonbit_string_literal_54.data);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0MPB13StringBuilder13write__objectGbE(_M0L18_2astring__builderS509, _M0L1xS505);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS509, (moonbit_string_t)moonbit_string_literal_55.data);
        #line 129 "/home/developer/.moon/lib/core/builtin/assert.mbt"
        _M0L7_2abindS510
        = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS509);
        moonbit_decref(_M0L18_2astring__builderS509);
        _M0L6_2atmpS2180 = Moonbit_array_length(_M0L7_2abindS510);
        _M0L9fail__msgS506
        = (struct _M0TPC16string10StringView){
          .$0 = _M0L7_2abindS510, .$1 = 0, .$2 = _M0L6_2atmpS2180
        };
        break;
      }
    }
    #line 131 "/home/developer/.moon/lib/core/builtin/assert.mbt"
    _result_4015 = _M0FPB4failGuE(_M0L9fail__msgS506, _M0L3locS511);
    moonbit_decref(_M0L9fail__msgS506.$0);
    return _result_4015;
  } else {
    int32_t _M0L6_2atmpS2181 = 0;
    struct moonbit_result_0 _result_4016;
    _result_4016.tag = 1;
    _result_4016.data.ok = _M0L6_2atmpS2181;
    return _result_4016;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__view(
  struct _M0TPB13StringBuilder* _M0L4selfS504,
  struct _M0TPC16string10StringView _M0L3strS503
) {
  int32_t _M0L3endS2178;
  int32_t _M0L5startS2179;
  int32_t _M0L8str__lenS502;
  int32_t _M0L3lenS2171;
  int32_t _M0L6_2atmpS2170;
  uint16_t* _M0L4dataS2172;
  int32_t _M0L3lenS2173;
  moonbit_string_t _M0L6_2atmpS2174;
  int32_t _M0L6_2atmpS2175;
  int32_t _M0L3lenS2177;
  int32_t _M0L6_2atmpS2176;
  #line 131 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3endS2178 = _M0L3strS503.$2;
  _M0L5startS2179 = _M0L3strS503.$1;
  _M0L8str__lenS502 = _M0L3endS2178 - _M0L5startS2179;
  _M0L3lenS2171 = _M0L4selfS504->$1;
  _M0L6_2atmpS2170 = _M0L3lenS2171 + _M0L8str__lenS502;
  #line 136 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS504, _M0L6_2atmpS2170);
  _M0L4dataS2172 = _M0L4selfS504->$0;
  _M0L3lenS2173 = _M0L4selfS504->$1;
  moonbit_incref(_M0L4dataS2172);
  #line 139 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS2174 = _M0MPC16string10StringView4data(_M0L3strS503);
  #line 140 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS2175 = _M0MPC16string10StringView13start__offset(_M0L3strS503);
  #line 137 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS2172, _M0L3lenS2173, _M0L6_2atmpS2174, _M0L6_2atmpS2175, _M0L8str__lenS502);
  moonbit_decref(_M0L4dataS2172);
  moonbit_decref(_M0L6_2atmpS2174);
  _M0L3lenS2177 = _M0L4selfS504->$1;
  _M0L6_2atmpS2176 = _M0L3lenS2177 + _M0L8str__lenS502;
  _M0L4selfS504->$1 = _M0L6_2atmpS2176;
  return 0;
}

int64_t _M0MPC16string6String29offset__of__nth__char_2einner(
  moonbit_string_t _M0L4selfS499,
  int32_t _M0L1iS500,
  int32_t _M0L13start__offsetS501,
  int64_t _M0L11end__offsetS497
) {
  int32_t _M0L11end__offsetS496;
  #line 446 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L11end__offsetS497 == 4294967296ll) {
    _M0L11end__offsetS496 = Moonbit_array_length(_M0L4selfS499);
  } else {
    int64_t _M0L7_2aSomeS498 = _M0L11end__offsetS497;
    _M0L11end__offsetS496 = (int32_t)_M0L7_2aSomeS498;
  }
  if (_M0L1iS500 >= 0) {
    #line 455 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0MPC16string6String30offset__of__nth__char__forward(_M0L4selfS499, _M0L1iS500, _M0L13start__offsetS501, _M0L11end__offsetS496);
  } else {
    int32_t _M0L6_2atmpS2169 = -_M0L1iS500;
    #line 458 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0MPC16string6String31offset__of__nth__char__backward(_M0L4selfS499, _M0L6_2atmpS2169, _M0L13start__offsetS501, _M0L11end__offsetS496);
  }
}

int64_t _M0MPC16string6String30offset__of__nth__char__forward(
  moonbit_string_t _M0L4selfS494,
  int32_t _M0L1nS492,
  int32_t _M0L13start__offsetS488,
  int32_t _M0L11end__offsetS489
) {
  int32_t _if__result_4017;
  #line 375 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L13start__offsetS488 >= 0) {
    _if__result_4017 = _M0L13start__offsetS488 <= _M0L11end__offsetS489;
  } else {
    _if__result_4017 = 0;
  }
  if (_if__result_4017) {
    int32_t _M0L13utf16__offsetS490 = _M0L13start__offsetS488;
    int32_t _M0L11char__countS491 = 0;
    while (1) {
      int32_t _if__result_4019;
      if (_M0L13utf16__offsetS490 < _M0L11end__offsetS489) {
        _if__result_4019 = _M0L11char__countS491 < _M0L1nS492;
      } else {
        _if__result_4019 = 0;
      }
      if (_if__result_4019) {
        int32_t _M0L1cS493 = _M0L4selfS494[_M0L13utf16__offsetS490];
        #line 389 "/home/developer/.moon/lib/core/builtin/string.mbt"
        if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L1cS493)) {
          int32_t _M0L6_2atmpS2165 = _M0L13utf16__offsetS490 + 2;
          int32_t _M0L6_2atmpS2166 = _M0L11char__countS491 + 1;
          _M0L13utf16__offsetS490 = _M0L6_2atmpS2165;
          _M0L11char__countS491 = _M0L6_2atmpS2166;
          continue;
        } else {
          int32_t _M0L6_2atmpS2167 = _M0L13utf16__offsetS490 + 1;
          int32_t _M0L6_2atmpS2168 = _M0L11char__countS491 + 1;
          _M0L13utf16__offsetS490 = _M0L6_2atmpS2167;
          _M0L11char__countS491 = _M0L6_2atmpS2168;
          continue;
        }
      } else {
        int32_t _if__result_4020;
        if (_M0L11char__countS491 < _M0L1nS492) {
          _if__result_4020 = 1;
        } else {
          _if__result_4020 = _M0L13utf16__offsetS490 >= _M0L11end__offsetS489;
        }
        if (_if__result_4020) {
          return 4294967296ll;
        } else {
          return (int64_t)_M0L13utf16__offsetS490;
        }
      }
      break;
    }
  } else {
    #line 382 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0FPC15abort5abortGOiE((moonbit_string_t)moonbit_string_literal_56.data);
  }
}

int64_t _M0MPC16string6String31offset__of__nth__char__backward(
  moonbit_string_t _M0L4selfS485,
  int32_t _M0L1nS483,
  int32_t _M0L13start__offsetS482,
  int32_t _M0L11end__offsetS487
) {
  int32_t _M0L13utf16__offsetS480;
  int32_t _M0L11char__countS481;
  #line 411 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L13utf16__offsetS480 = _M0L11end__offsetS487;
  _M0L11char__countS481 = 0;
  while (1) {
    int32_t _M0L6_2atmpS2159 = _M0L13utf16__offsetS480 - 1;
    int32_t _if__result_4022;
    if (_M0L6_2atmpS2159 >= _M0L13start__offsetS482) {
      _if__result_4022 = _M0L11char__countS481 < _M0L1nS483;
    } else {
      _if__result_4022 = 0;
    }
    if (_if__result_4022) {
      int32_t _M0L6_2atmpS2164 = _M0L13utf16__offsetS480 - 1;
      int32_t _M0L1cS484 = _M0L4selfS485[_M0L6_2atmpS2164];
      #line 424 "/home/developer/.moon/lib/core/builtin/string.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L1cS484)) {
        int32_t _M0L6_2atmpS2160 = _M0L13utf16__offsetS480 - 2;
        int32_t _M0L6_2atmpS2161 = _M0L11char__countS481 + 1;
        _M0L13utf16__offsetS480 = _M0L6_2atmpS2160;
        _M0L11char__countS481 = _M0L6_2atmpS2161;
        continue;
      } else {
        int32_t _M0L6_2atmpS2162 = _M0L13utf16__offsetS480 - 1;
        int32_t _M0L6_2atmpS2163 = _M0L11char__countS481 + 1;
        _M0L13utf16__offsetS480 = _M0L6_2atmpS2162;
        _M0L11char__countS481 = _M0L6_2atmpS2163;
        continue;
      }
    } else {
      int32_t _if__result_4023;
      if (_M0L11char__countS481 < _M0L1nS483) {
        _if__result_4023 = 1;
      } else {
        _if__result_4023 = _M0L13utf16__offsetS480 < _M0L13start__offsetS482;
      }
      if (_if__result_4023) {
        return 4294967296ll;
      } else {
        return (int64_t)_M0L13utf16__offsetS480;
      }
    }
    break;
  }
}

int32_t _M0IPC16string10StringViewPB4Show6output(
  struct _M0TPC16string10StringView _M0L4selfS479,
  struct _M0TPB6Logger _M0L6loggerS478
) {
  #line 203 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  #line 204 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L6loggerS478.$0->$method_2(_M0L6loggerS478.$1, _M0L4selfS479);
  return 0;
}

moonbit_string_t _M0MPC16string10StringView9to__owned(
  struct _M0TPC16string10StringView _M0L4selfS477
) {
  moonbit_string_t _M0L3strS2156;
  int32_t _M0L5startS2157;
  int32_t _M0L3endS2158;
  moonbit_string_t _result_4024;
  #line 196 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS2156 = _M0L4selfS477.$0;
  _M0L5startS2157 = _M0L4selfS477.$1;
  _M0L3endS2158 = _M0L4selfS477.$2;
  moonbit_incref(_M0L3strS2156);
  #line 199 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_4024
  = _M0MPC16string6String17unsafe__substring(_M0L3strS2156, _M0L5startS2157, _M0L3endS2158);
  moonbit_decref(_M0L3strS2156);
  return _result_4024;
}

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t _M0L3strS474,
  int32_t _M0L5startS472,
  int32_t _M0L3endS473
) {
  int32_t _if__result_4025;
  int32_t _M0L3lenS475;
  int32_t _M0L6_2atmpS2154;
  int32_t _M0L6_2atmpS2155;
  moonbit_bytes_t _M0L5bytesS476;
  moonbit_bytes_t _M0L6_2atmpS2153;
  moonbit_string_t _result_4026;
  #line 91 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L5startS472 == 0) {
    int32_t _M0L6_2atmpS2152 = Moonbit_array_length(_M0L3strS474);
    _if__result_4025 = _M0L3endS473 == _M0L6_2atmpS2152;
  } else {
    _if__result_4025 = 0;
  }
  if (_if__result_4025) {
    moonbit_incref(_M0L3strS474);
    return _M0L3strS474;
  }
  _M0L3lenS475 = _M0L3endS473 - _M0L5startS472;
  _M0L6_2atmpS2154 = _M0L3lenS475 * 2;
  #line 101 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS2155 = _M0IPC14byte4BytePB7Default7default();
  _M0L5bytesS476
  = (moonbit_bytes_t)moonbit_make_bytes(_M0L6_2atmpS2154, _M0L6_2atmpS2155);
  #line 102 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0MPC15array10FixedArray18blit__from__string(_M0L5bytesS476, 0, _M0L3strS474, _M0L5startS472, _M0L3lenS475);
  _M0L6_2atmpS2153 = _M0L5bytesS476;
  #line 103 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_4026
  = _M0MPC15bytes5Bytes29to__unchecked__string_2einner(_M0L6_2atmpS2153, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS2153);
  return _result_4026;
}

int32_t _M0IPC14byte4BytePB7Default7default() {
  #line 231 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  return 0;
}

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t _M0L4selfS467,
  int32_t _M0L6offsetS471,
  int64_t _M0L6lengthS469
) {
  int32_t _M0L3lenS466;
  int32_t _M0L6lengthS468;
  int32_t _if__result_4027;
  #line 77 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L3lenS466 = Moonbit_array_length(_M0L4selfS467);
  if (_M0L6lengthS469 == 4294967296ll) {
    _M0L6lengthS468 = _M0L3lenS466 - _M0L6offsetS471;
  } else {
    int64_t _M0L7_2aSomeS470 = _M0L6lengthS469;
    _M0L6lengthS468 = (int32_t)_M0L7_2aSomeS470;
  }
  if (_M0L6offsetS471 >= 0) {
    if (_M0L6lengthS468 >= 0) {
      int32_t _M0L6_2atmpS2151 = _M0L6offsetS471 + _M0L6lengthS468;
      _if__result_4027 = _M0L6_2atmpS2151 <= _M0L3lenS466;
    } else {
      _if__result_4027 = 0;
    }
  } else {
    _if__result_4027 = 0;
  }
  if (_if__result_4027) {
    moonbit_incref(_M0L4selfS467);
    #line 85 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    return _M0FPB19unsafe__sub__string(_M0L4selfS467, _M0L6offsetS471, _M0L6lengthS468);
  } else {
    #line 84 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t _M0L4selfS458,
  int32_t _M0L13bytes__offsetS453,
  moonbit_string_t _M0L3strS460,
  int32_t _M0L11str__offsetS456,
  int32_t _M0L6lengthS454
) {
  int32_t _M0L6_2atmpS2150;
  int32_t _M0L6_2atmpS2149;
  int32_t _M0L2e1S452;
  int32_t _M0L6_2atmpS2148;
  int32_t _M0L2e2S455;
  int32_t _M0L4len1S457;
  int32_t _M0L4len2S459;
  int32_t _if__result_4028;
  #line 125 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L6_2atmpS2150 = _M0L6lengthS454 * 2;
  _M0L6_2atmpS2149 = _M0L13bytes__offsetS453 + _M0L6_2atmpS2150;
  _M0L2e1S452 = _M0L6_2atmpS2149 - 1;
  _M0L6_2atmpS2148 = _M0L11str__offsetS456 + _M0L6lengthS454;
  _M0L2e2S455 = _M0L6_2atmpS2148 - 1;
  _M0L4len1S457 = Moonbit_array_length(_M0L4selfS458);
  _M0L4len2S459 = Moonbit_array_length(_M0L3strS460);
  if (_M0L6lengthS454 >= 0) {
    if (_M0L13bytes__offsetS453 >= 0) {
      if (_M0L2e1S452 < _M0L4len1S457) {
        if (_M0L11str__offsetS456 >= 0) {
          _if__result_4028 = _M0L2e2S455 < _M0L4len2S459;
        } else {
          _if__result_4028 = 0;
        }
      } else {
        _if__result_4028 = 0;
      }
    } else {
      _if__result_4028 = 0;
    }
  } else {
    _if__result_4028 = 0;
  }
  if (_if__result_4028) {
    int32_t _M0L16end__str__offsetS461 =
      _M0L11str__offsetS456 + _M0L6lengthS454;
    int32_t _M0L1iS462 = _M0L11str__offsetS456;
    int32_t _M0L1jS463 = _M0L13bytes__offsetS453;
    while (1) {
      if (_M0L1iS462 < _M0L16end__str__offsetS461) {
        int32_t _M0L6_2atmpS2145 = _M0L3strS460[_M0L1iS462];
        int32_t _M0L6_2atmpS2144 = (int32_t)_M0L6_2atmpS2145;
        uint32_t _M0L1cS464 = *(uint32_t*)&_M0L6_2atmpS2144;
        uint32_t _M0L6_2atmpS2140 = _M0L1cS464 & 255u;
        int32_t _M0L6_2atmpS2139;
        int32_t _M0L6_2atmpS2141;
        uint32_t _M0L6_2atmpS2143;
        int32_t _M0L6_2atmpS2142;
        int32_t _M0L6_2atmpS2146;
        int32_t _M0L6_2atmpS2147;
        #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS2139 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS2140);
        if (
          _M0L1jS463 < 0 || _M0L1jS463 >= Moonbit_array_length(_M0L4selfS458)
        ) {
          #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS458[_M0L1jS463] = _M0L6_2atmpS2139;
        _M0L6_2atmpS2141 = _M0L1jS463 + 1;
        _M0L6_2atmpS2143 = _M0L1cS464 >> 8;
        #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS2142 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS2143);
        if (
          _M0L6_2atmpS2141 < 0
          || _M0L6_2atmpS2141 >= Moonbit_array_length(_M0L4selfS458)
        ) {
          #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS458[_M0L6_2atmpS2141] = _M0L6_2atmpS2142;
        _M0L6_2atmpS2146 = _M0L1iS462 + 1;
        _M0L6_2atmpS2147 = _M0L1jS463 + 2;
        _M0L1iS462 = _M0L6_2atmpS2146;
        _M0L1jS463 = _M0L6_2atmpS2147;
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

int32_t _M0MPC14uint4UInt8to__byte(uint32_t _M0L4selfS451) {
  int32_t _M0L6_2atmpS2138;
  #line 2519 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS2138 = *(int32_t*)&_M0L4selfS451;
  return _M0L6_2atmpS2138 & 0xff;
}

moonbit_string_t* _M0MPC15array5Array6bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS446
) {
  moonbit_string_t* _M0L8_2afieldS3659;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3659 = _M0L4selfS446->$0;
  moonbit_incref(_M0L8_2afieldS3659);
  return _M0L8_2afieldS3659;
}

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS447
) {
  struct _M0TUsiE** _M0L8_2afieldS3660;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3660 = _M0L4selfS447->$0;
  moonbit_incref(_M0L8_2afieldS3660);
  return _M0L8_2afieldS3660;
}

struct _M0TP26biolab8bio__seq8TreeNode** _M0MPC15array5Array6bufferGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L4selfS448
) {
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L8_2afieldS3661;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3661 = _M0L4selfS448->$0;
  moonbit_incref(_M0L8_2afieldS3661);
  return _M0L8_2afieldS3661;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0MPC15array5Array6bufferGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE* _M0L4selfS449
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L8_2afieldS3662;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3662 = _M0L4selfS449->$0;
  moonbit_incref(_M0L8_2afieldS3662);
  return _M0L8_2afieldS3662;
}

int32_t* _M0MPC15array5Array6bufferGiE(struct _M0TPB5ArrayGiE* _M0L4selfS450) {
  int32_t* _M0L8_2afieldS3663;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3663 = _M0L4selfS450->$0;
  moonbit_incref(_M0L8_2afieldS3663);
  return _M0L8_2afieldS3663;
}

struct _M0TPB4IterGsE* _M0MPB4Iter3newGsE(
  struct _M0TWEOs* _M0L1fS445,
  int64_t _M0L10size__hintS442
) {
  int64_t _M0L10size__hintS441;
  struct _M0TPB4IterGsE* _block_4030;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS442 == 4294967296ll) {
    _M0L10size__hintS441 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS443 = _M0L10size__hintS442;
    int32_t _M0L4_2anS444 = (int32_t)_M0L7_2aSomeS443;
    if (_M0L4_2anS444 > 0) {
      _M0L10size__hintS441 = (int64_t)_M0L4_2anS444;
    } else {
      _M0L10size__hintS441 = _M0MPB4Iter3newN6constrS9988GsE;
    }
  }
  moonbit_incref(_M0L1fS445);
  _block_4030
  = (struct _M0TPB4IterGsE*)moonbit_malloc(sizeof(struct _M0TPB4IterGsE));
  Moonbit_object_header(_block_4030)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 84, 0);
  _block_4030->$0 = _M0L1fS445;
  _block_4030->$1 = _M0L10size__hintS441;
  return _block_4030;
}

struct moonbit_result_0 _M0FPB10assert__eqGiE(
  int32_t _M0L1aS434,
  int32_t _M0L1bS435,
  moonbit_string_t _M0L3msgS437,
  moonbit_string_t _M0L3locS440
) {
  #line 45 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  if (_M0L1aS434 != _M0L1bS435) {
    moonbit_string_t _M0L9fail__msgS436;
    int32_t _M0L6_2atmpS2134;
    struct _M0TPC16string10StringView _M0L6_2atmpS2133;
    struct moonbit_result_0 _result_4031;
    if (_M0L3msgS437 == 0) {
      struct _M0TPB13StringBuilder* _M0L18_2astring__builderS439;
      moonbit_string_t _M0L6_2atmpS2135;
      moonbit_string_t _M0L6_2atmpS2136;
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L18_2astring__builderS439
      = _M0MPB13StringBuilder21StringBuilder_2einner(6);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS439, (moonbit_string_t)moonbit_string_literal_54.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L6_2atmpS2135 = _M0FPB13debug__stringGiE(_M0L1aS434);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS439, _M0L6_2atmpS2135);
      moonbit_decref(_M0L6_2atmpS2135);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS439, (moonbit_string_t)moonbit_string_literal_57.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L6_2atmpS2136 = _M0FPB13debug__stringGiE(_M0L1bS435);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS439, _M0L6_2atmpS2136);
      moonbit_decref(_M0L6_2atmpS2136);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS439, (moonbit_string_t)moonbit_string_literal_54.data);
      #line 57 "/home/developer/.moon/lib/core/builtin/assert.mbt"
      _M0L9fail__msgS436
      = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS439);
      moonbit_decref(_M0L18_2astring__builderS439);
    } else {
      moonbit_string_t _M0L7_2aSomeS438 = _M0L3msgS437;
      if (_M0L7_2aSomeS438) {
        moonbit_incref(_M0L7_2aSomeS438);
      }
      _M0L9fail__msgS436 = _M0L7_2aSomeS438;
    }
    _M0L6_2atmpS2134 = Moonbit_array_length(_M0L9fail__msgS436);
    _M0L6_2atmpS2133
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L9fail__msgS436, .$1 = 0, .$2 = _M0L6_2atmpS2134
    };
    #line 59 "/home/developer/.moon/lib/core/builtin/assert.mbt"
    _result_4031 = _M0FPB4failGuE(_M0L6_2atmpS2133, _M0L3locS440);
    moonbit_decref(_M0L6_2atmpS2133.$0);
    return _result_4031;
  } else {
    int32_t _M0L6_2atmpS2137 = 0;
    struct moonbit_result_0 _result_4032;
    _result_4032.tag = 1;
    _result_4032.data.ok = _M0L6_2atmpS2137;
    return _result_4032;
  }
}

struct moonbit_result_0 _M0FPB4failGuE(
  struct _M0TPC16string10StringView _M0L3msgS433,
  moonbit_string_t _M0L3locS432
) {
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS431;
  moonbit_string_t _M0L6_2atmpS2132;
  void* _M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS2131;
  struct moonbit_result_0 _result_4033;
  #line 56 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L18_2astring__builderS431
  = _M0MPB13StringBuilder21StringBuilder_2einner(9);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB13StringBuilder13write__objectGRPB9SourceLocE(_M0L18_2astring__builderS431, _M0L3locS432);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS431, (moonbit_string_t)moonbit_string_literal_58.data);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB13StringBuilder13write__objectGRPC16string10StringViewE(_M0L18_2astring__builderS431, _M0L3msgS433);
  #line 58 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L6_2atmpS2132
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS431);
  moonbit_decref(_M0L18_2astring__builderS431);
  _M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS2131
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure));
  Moonbit_object_header(_M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS2131)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 87, 0);
  ((struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS2131)->$0
  = _M0L6_2atmpS2132;
  _result_4033.tag = 0;
  _result_4033.data.err
  = _M0L48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailureS2131;
  return _result_4033;
}

moonbit_string_t _M0FPB13debug__stringGiE(int32_t _M0L1tS430) {
  struct _M0TPB13StringBuilder* _M0L3bufS429;
  struct _M0TPB6Logger _M0L6_2atmpS2130;
  moonbit_string_t _result_4034;
  #line 16 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  #line 17 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _M0L3bufS429 = _M0MPB13StringBuilder21StringBuilder_2einner(50);
  moonbit_incref(_M0L3bufS429);
  _M0L6_2atmpS2130
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS429
  };
  #line 18 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L1tS430, _M0L6_2atmpS2130);
  if (_M0L6_2atmpS2130.$1) {
    moonbit_decref(_M0L6_2atmpS2130.$1);
  }
  #line 19 "/home/developer/.moon/lib/core/builtin/assert.mbt"
  _result_4034 = _M0MPB13StringBuilder10to__string(_M0L3bufS429);
  moonbit_decref(_M0L3bufS429);
  return _result_4034;
}

moonbit_string_t _M0MPC16uint646UInt6418to__string_2einner(
  uint64_t _M0L4selfS421,
  int32_t _M0L5radixS420
) {
  int32_t _if__result_4035;
  uint16_t* _M0L6bufferS422;
  #line 607 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS420 < 2) {
    _if__result_4035 = 1;
  } else {
    _if__result_4035 = _M0L5radixS420 > 36;
  }
  if (_if__result_4035) {
    #line 611 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_59.data);
  }
  if (_M0L4selfS421 == 0ull) {
    return (moonbit_string_t)moonbit_string_literal_46.data;
  }
  switch (_M0L5radixS420) {
    case 10: {
      int32_t _M0L3lenS423;
      uint16_t* _M0L6bufferS424;
      #line 622 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L3lenS423 = _M0FPB12dec__count64(_M0L4selfS421);
      _M0L6bufferS424 = (uint16_t*)moonbit_make_string(_M0L3lenS423, 0);
      #line 624 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB22int64__to__string__dec(_M0L6bufferS424, _M0L4selfS421, 0, _M0L3lenS423);
      _M0L6bufferS422 = _M0L6bufferS424;
      break;
    }
    
    case 16: {
      int32_t _M0L3lenS425;
      uint16_t* _M0L6bufferS426;
      #line 628 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L3lenS425 = _M0FPB12hex__count64(_M0L4selfS421);
      _M0L6bufferS426 = (uint16_t*)moonbit_make_string(_M0L3lenS425, 0);
      #line 630 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB22int64__to__string__hex(_M0L6bufferS426, _M0L4selfS421, 0, _M0L3lenS425);
      _M0L6bufferS422 = _M0L6bufferS426;
      break;
    }
    default: {
      int32_t _M0L3lenS427;
      uint16_t* _M0L6bufferS428;
      #line 634 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L3lenS427 = _M0FPB14radix__count64(_M0L4selfS421, _M0L5radixS420);
      _M0L6bufferS428 = (uint16_t*)moonbit_make_string(_M0L3lenS427, 0);
      #line 636 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB26int64__to__string__generic(_M0L6bufferS428, _M0L4selfS421, 0, _M0L3lenS427, _M0L5radixS420);
      _M0L6bufferS422 = _M0L6bufferS428;
      break;
    }
  }
  return _M0L6bufferS422;
}

int32_t _M0FPB22int64__to__string__dec(
  uint16_t* _M0L6bufferS406,
  uint64_t _M0L3numS418,
  int32_t _M0L12digit__startS407,
  int32_t _M0L10total__lenS419
) {
  int32_t _M0L6_2atmpS2129;
  uint64_t _M0L3numS396;
  int32_t _M0L6offsetS397;
  #line 493 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS2129 = _M0L10total__lenS419 - _M0L12digit__startS407;
  _M0L3numS396 = _M0L3numS418;
  _M0L6offsetS397 = _M0L6_2atmpS2129;
  while (1) {
    if (_M0L3numS396 >= 10000ull) {
      uint64_t _M0L1tS398 = _M0L3numS396 / 10000ull;
      uint64_t _M0L6_2atmpS2106 = _M0L3numS396 % 10000ull;
      int32_t _M0L1rS399 = (int32_t)_M0L6_2atmpS2106;
      int32_t _M0L2d1S400 = _M0L1rS399 / 100;
      int32_t _M0L2d2S401 = _M0L1rS399 % 100;
      int32_t _M0L6_2atmpS2105 = _M0L2d1S400 / 10;
      int32_t _M0L6_2atmpS2104 = 48 + _M0L6_2atmpS2105;
      int32_t _M0L6d1__hiS402 = (uint16_t)_M0L6_2atmpS2104;
      int32_t _M0L6_2atmpS2103 = _M0L2d1S400 % 10;
      int32_t _M0L6_2atmpS2102 = 48 + _M0L6_2atmpS2103;
      int32_t _M0L6d1__loS403 = (uint16_t)_M0L6_2atmpS2102;
      int32_t _M0L6_2atmpS2101 = _M0L2d2S401 / 10;
      int32_t _M0L6_2atmpS2100 = 48 + _M0L6_2atmpS2101;
      int32_t _M0L6d2__hiS404 = (uint16_t)_M0L6_2atmpS2100;
      int32_t _M0L6_2atmpS2099 = _M0L2d2S401 % 10;
      int32_t _M0L6_2atmpS2098 = 48 + _M0L6_2atmpS2099;
      int32_t _M0L6d2__loS405 = (uint16_t)_M0L6_2atmpS2098;
      int32_t _M0L6_2atmpS2090 = _M0L12digit__startS407 + _M0L6offsetS397;
      int32_t _M0L6_2atmpS2089 = _M0L6_2atmpS2090 - 4;
      int32_t _M0L6_2atmpS2092;
      int32_t _M0L6_2atmpS2091;
      int32_t _M0L6_2atmpS2094;
      int32_t _M0L6_2atmpS2093;
      int32_t _M0L6_2atmpS2096;
      int32_t _M0L6_2atmpS2095;
      int32_t _M0L6_2atmpS2097;
      _M0L6bufferS406[_M0L6_2atmpS2089] = _M0L6d1__hiS402;
      _M0L6_2atmpS2092 = _M0L12digit__startS407 + _M0L6offsetS397;
      _M0L6_2atmpS2091 = _M0L6_2atmpS2092 - 3;
      _M0L6bufferS406[_M0L6_2atmpS2091] = _M0L6d1__loS403;
      _M0L6_2atmpS2094 = _M0L12digit__startS407 + _M0L6offsetS397;
      _M0L6_2atmpS2093 = _M0L6_2atmpS2094 - 2;
      _M0L6bufferS406[_M0L6_2atmpS2093] = _M0L6d2__hiS404;
      _M0L6_2atmpS2096 = _M0L12digit__startS407 + _M0L6offsetS397;
      _M0L6_2atmpS2095 = _M0L6_2atmpS2096 - 1;
      _M0L6bufferS406[_M0L6_2atmpS2095] = _M0L6d2__loS405;
      _M0L6_2atmpS2097 = _M0L6offsetS397 - 4;
      _M0L3numS396 = _M0L1tS398;
      _M0L6offsetS397 = _M0L6_2atmpS2097;
      continue;
    } else {
      int32_t _M0L6_2atmpS2128 = (int32_t)_M0L3numS396;
      int32_t _M0L9remainingS409 = _M0L6_2atmpS2128;
      int32_t _M0L6offsetS410 = _M0L6offsetS397;
      while (1) {
        if (_M0L9remainingS409 >= 100) {
          int32_t _M0L1tS411 = _M0L9remainingS409 / 100;
          int32_t _M0L1dS412 = _M0L9remainingS409 % 100;
          int32_t _M0L6_2atmpS2115 = _M0L1dS412 / 10;
          int32_t _M0L6_2atmpS2114 = 48 + _M0L6_2atmpS2115;
          int32_t _M0L5d__hiS413 = (uint16_t)_M0L6_2atmpS2114;
          int32_t _M0L6_2atmpS2113 = _M0L1dS412 % 10;
          int32_t _M0L6_2atmpS2112 = 48 + _M0L6_2atmpS2113;
          int32_t _M0L5d__loS414 = (uint16_t)_M0L6_2atmpS2112;
          int32_t _M0L6_2atmpS2108 = _M0L12digit__startS407 + _M0L6offsetS410;
          int32_t _M0L6_2atmpS2107 = _M0L6_2atmpS2108 - 2;
          int32_t _M0L6_2atmpS2110;
          int32_t _M0L6_2atmpS2109;
          int32_t _M0L6_2atmpS2111;
          _M0L6bufferS406[_M0L6_2atmpS2107] = _M0L5d__hiS413;
          _M0L6_2atmpS2110 = _M0L12digit__startS407 + _M0L6offsetS410;
          _M0L6_2atmpS2109 = _M0L6_2atmpS2110 - 1;
          _M0L6bufferS406[_M0L6_2atmpS2109] = _M0L5d__loS414;
          _M0L6_2atmpS2111 = _M0L6offsetS410 - 2;
          _M0L9remainingS409 = _M0L1tS411;
          _M0L6offsetS410 = _M0L6_2atmpS2111;
          continue;
        } else if (_M0L9remainingS409 >= 10) {
          int32_t _M0L6_2atmpS2123 = _M0L9remainingS409 / 10;
          int32_t _M0L6_2atmpS2122 = 48 + _M0L6_2atmpS2123;
          int32_t _M0L5d__hiS416 = (uint16_t)_M0L6_2atmpS2122;
          int32_t _M0L6_2atmpS2121 = _M0L9remainingS409 % 10;
          int32_t _M0L6_2atmpS2120 = 48 + _M0L6_2atmpS2121;
          int32_t _M0L5d__loS417 = (uint16_t)_M0L6_2atmpS2120;
          int32_t _M0L6_2atmpS2117 = _M0L12digit__startS407 + _M0L6offsetS410;
          int32_t _M0L6_2atmpS2116 = _M0L6_2atmpS2117 - 2;
          int32_t _M0L6_2atmpS2119;
          int32_t _M0L6_2atmpS2118;
          _M0L6bufferS406[_M0L6_2atmpS2116] = _M0L5d__hiS416;
          _M0L6_2atmpS2119 = _M0L12digit__startS407 + _M0L6offsetS410;
          _M0L6_2atmpS2118 = _M0L6_2atmpS2119 - 1;
          _M0L6bufferS406[_M0L6_2atmpS2118] = _M0L5d__loS417;
        } else {
          int32_t _M0L6_2atmpS2127 = _M0L12digit__startS407 + _M0L6offsetS410;
          int32_t _M0L6_2atmpS2124 = _M0L6_2atmpS2127 - 1;
          int32_t _M0L6_2atmpS2126 = 48 + _M0L9remainingS409;
          int32_t _M0L6_2atmpS2125 = (uint16_t)_M0L6_2atmpS2126;
          _M0L6bufferS406[_M0L6_2atmpS2124] = _M0L6_2atmpS2125;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB26int64__to__string__generic(
  uint16_t* _M0L6bufferS386,
  uint64_t _M0L3numS390,
  int32_t _M0L12digit__startS387,
  int32_t _M0L10total__lenS389,
  int32_t _M0L5radixS380
) {
  uint64_t _M0L4baseS379;
  int32_t _M0L6_2atmpS2074;
  int32_t _M0L6_2atmpS2073;
  #line 462 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  #line 470 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS379 = _M0MPC13int3Int10to__uint64(_M0L5radixS380);
  _M0L6_2atmpS2074 = _M0L5radixS380 - 1;
  _M0L6_2atmpS2073 = _M0L5radixS380 & _M0L6_2atmpS2074;
  if (_M0L6_2atmpS2073 == 0) {
    int32_t _M0L5shiftS381;
    uint64_t _M0L4maskS382;
    int32_t _M0L6_2atmpS2081;
    int32_t _M0L6offsetS383;
    uint64_t _M0L1nS384;
    #line 473 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS381 = moonbit_ctz32(_M0L5radixS380);
    _M0L4maskS382 = _M0L4baseS379 - 1ull;
    _M0L6_2atmpS2081 = _M0L10total__lenS389 - _M0L12digit__startS387;
    _M0L6offsetS383 = _M0L6_2atmpS2081;
    _M0L1nS384 = _M0L3numS390;
    while (1) {
      if (_M0L1nS384 > 0ull) {
        uint64_t _M0L6_2atmpS2080 = _M0L1nS384 & _M0L4maskS382;
        int32_t _M0L5digitS385 = (int32_t)_M0L6_2atmpS2080;
        int32_t _M0L6_2atmpS2077 = _M0L12digit__startS387 + _M0L6offsetS383;
        int32_t _M0L6_2atmpS2075 = _M0L6_2atmpS2077 - 1;
        int32_t _M0L6_2atmpS2076 =
          ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L5digitS385];
        int32_t _M0L6_2atmpS2078;
        uint64_t _M0L6_2atmpS2079;
        _M0L6bufferS386[_M0L6_2atmpS2075] = _M0L6_2atmpS2076;
        _M0L6_2atmpS2078 = _M0L6offsetS383 - 1;
        _M0L6_2atmpS2079 = _M0L1nS384 >> (_M0L5shiftS381 & 63);
        _M0L6offsetS383 = _M0L6_2atmpS2078;
        _M0L1nS384 = _M0L6_2atmpS2079;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS2088 = _M0L10total__lenS389 - _M0L12digit__startS387;
    int32_t _M0L6offsetS391 = _M0L6_2atmpS2088;
    uint64_t _M0L1nS392 = _M0L3numS390;
    while (1) {
      if (_M0L1nS392 > 0ull) {
        uint64_t _M0L1qS393 = _M0L1nS392 / _M0L4baseS379;
        uint64_t _M0L6_2atmpS2087 = _M0L1qS393 * _M0L4baseS379;
        uint64_t _M0L6_2atmpS2086 = _M0L1nS392 - _M0L6_2atmpS2087;
        int32_t _M0L5digitS394 = (int32_t)_M0L6_2atmpS2086;
        int32_t _M0L6_2atmpS2084 = _M0L12digit__startS387 + _M0L6offsetS391;
        int32_t _M0L6_2atmpS2082 = _M0L6_2atmpS2084 - 1;
        int32_t _M0L6_2atmpS2083 =
          ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L5digitS394];
        int32_t _M0L6_2atmpS2085;
        _M0L6bufferS386[_M0L6_2atmpS2082] = _M0L6_2atmpS2083;
        _M0L6_2atmpS2085 = _M0L6offsetS391 - 1;
        _M0L6offsetS391 = _M0L6_2atmpS2085;
        _M0L1nS392 = _M0L1qS393;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB22int64__to__string__hex(
  uint16_t* _M0L6bufferS373,
  uint64_t _M0L3numS378,
  int32_t _M0L12digit__startS374,
  int32_t _M0L10total__lenS377
) {
  int32_t _M0L6_2atmpS2072;
  int32_t _M0L6offsetS368;
  uint64_t _M0L1nS369;
  #line 434 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS2072 = _M0L10total__lenS377 - _M0L12digit__startS374;
  _M0L6offsetS368 = _M0L6_2atmpS2072;
  _M0L1nS369 = _M0L3numS378;
  while (1) {
    if (_M0L6offsetS368 >= 2) {
      uint64_t _M0L6_2atmpS2069 = _M0L1nS369 & 255ull;
      int32_t _M0L9byte__valS370 = (int32_t)_M0L6_2atmpS2069;
      int32_t _M0L2hiS371 = _M0L9byte__valS370 / 16;
      int32_t _M0L2loS372 = _M0L9byte__valS370 % 16;
      int32_t _M0L6_2atmpS2063 = _M0L12digit__startS374 + _M0L6offsetS368;
      int32_t _M0L6_2atmpS2061 = _M0L6_2atmpS2063 - 2;
      int32_t _M0L6_2atmpS2062 =
        ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L2hiS371];
      int32_t _M0L6_2atmpS2066;
      int32_t _M0L6_2atmpS2064;
      int32_t _M0L6_2atmpS2065;
      int32_t _M0L6_2atmpS2067;
      uint64_t _M0L6_2atmpS2068;
      _M0L6bufferS373[_M0L6_2atmpS2061] = _M0L6_2atmpS2062;
      _M0L6_2atmpS2066 = _M0L12digit__startS374 + _M0L6offsetS368;
      _M0L6_2atmpS2064 = _M0L6_2atmpS2066 - 1;
      _M0L6_2atmpS2065
      = ((moonbit_string_t)moonbit_string_literal_60.data)[
        _M0L2loS372
      ];
      _M0L6bufferS373[_M0L6_2atmpS2064] = _M0L6_2atmpS2065;
      _M0L6_2atmpS2067 = _M0L6offsetS368 - 2;
      _M0L6_2atmpS2068 = _M0L1nS369 >> 8;
      _M0L6offsetS368 = _M0L6_2atmpS2067;
      _M0L1nS369 = _M0L6_2atmpS2068;
      continue;
    } else if (_M0L6offsetS368 == 1) {
      uint64_t _M0L6_2atmpS2071 = _M0L1nS369 & 15ull;
      int32_t _M0L6nibbleS376 = (int32_t)_M0L6_2atmpS2071;
      int32_t _M0L6_2atmpS2070 =
        ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L6nibbleS376];
      _M0L6bufferS373[_M0L12digit__startS374] = _M0L6_2atmpS2070;
    }
    break;
  }
  return 0;
}

int32_t _M0FPB14radix__count64(
  uint64_t _M0L5valueS362,
  int32_t _M0L5radixS364
) {
  uint64_t _M0L4baseS363;
  uint64_t _M0L3numS365;
  int32_t _M0L5countS366;
  #line 419 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS362 == 0ull) {
    return 1;
  }
  #line 424 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS363 = _M0MPC13int3Int10to__uint64(_M0L5radixS364);
  _M0L3numS365 = _M0L5valueS362;
  _M0L5countS366 = 0;
  while (1) {
    if (_M0L3numS365 > 0ull) {
      uint64_t _M0L6_2atmpS2059 = _M0L3numS365 / _M0L4baseS363;
      int32_t _M0L6_2atmpS2060 = _M0L5countS366 + 1;
      _M0L3numS365 = _M0L6_2atmpS2059;
      _M0L5countS366 = _M0L6_2atmpS2060;
      continue;
    } else {
      return _M0L5countS366;
    }
    break;
  }
}

int32_t _M0FPB12hex__count64(uint64_t _M0L5valueS360) {
  #line 407 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS360 == 0ull) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS361;
    int32_t _M0L6_2atmpS2058;
    int32_t _M0L6_2atmpS2057;
    #line 412 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS361 = moonbit_clz64(_M0L5valueS360);
    _M0L6_2atmpS2058 = 63 - _M0L14leading__zerosS361;
    _M0L6_2atmpS2057 = _M0L6_2atmpS2058 / 4;
    return _M0L6_2atmpS2057 + 1;
  }
}

int32_t _M0FPB12dec__count64(uint64_t _M0L5valueS359) {
  #line 343 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS359 >= 10000000000ull) {
    if (_M0L5valueS359 >= 100000000000000ull) {
      if (_M0L5valueS359 >= 10000000000000000ull) {
        if (_M0L5valueS359 >= 1000000000000000000ull) {
          if (_M0L5valueS359 >= 10000000000000000000ull) {
            return 20;
          } else {
            return 19;
          }
        } else if (_M0L5valueS359 >= 100000000000000000ull) {
          return 18;
        } else {
          return 17;
        }
      } else if (_M0L5valueS359 >= 1000000000000000ull) {
        return 16;
      } else {
        return 15;
      }
    } else if (_M0L5valueS359 >= 1000000000000ull) {
      if (_M0L5valueS359 >= 10000000000000ull) {
        return 14;
      } else {
        return 13;
      }
    } else if (_M0L5valueS359 >= 100000000000ull) {
      return 12;
    } else {
      return 11;
    }
  } else if (_M0L5valueS359 >= 100000ull) {
    if (_M0L5valueS359 >= 10000000ull) {
      if (_M0L5valueS359 >= 1000000000ull) {
        return 10;
      } else if (_M0L5valueS359 >= 100000000ull) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS359 >= 1000000ull) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS359 >= 1000ull) {
    if (_M0L5valueS359 >= 10000ull) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS359 >= 100ull) {
    return 3;
  } else if (_M0L5valueS359 >= 10ull) {
    return 2;
  } else {
    return 1;
  }
}

moonbit_string_t _M0MPC13int3Int18to__string_2einner(
  int32_t _M0L4selfS343,
  int32_t _M0L5radixS342
) {
  int32_t _if__result_4042;
  int32_t _M0L12is__negativeS344;
  uint32_t _M0L3numS345;
  uint16_t* _M0L6bufferS346;
  #line 209 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS342 < 2) {
    _if__result_4042 = 1;
  } else {
    _if__result_4042 = _M0L5radixS342 > 36;
  }
  if (_if__result_4042) {
    #line 213 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_59.data);
  }
  if (_M0L4selfS343 == 0) {
    return (moonbit_string_t)moonbit_string_literal_46.data;
  }
  _M0L12is__negativeS344 = _M0L4selfS343 < 0;
  if (_M0L12is__negativeS344) {
    int32_t _M0L6_2atmpS2056 = -_M0L4selfS343;
    _M0L3numS345 = *(uint32_t*)&_M0L6_2atmpS2056;
  } else {
    _M0L3numS345 = *(uint32_t*)&_M0L4selfS343;
  }
  switch (_M0L5radixS342) {
    case 10: {
      int32_t _M0L10digit__lenS347;
      int32_t _M0L6_2atmpS2053;
      int32_t _M0L10total__lenS348;
      uint16_t* _M0L6bufferS349;
      int32_t _M0L12digit__startS350;
      #line 235 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS347 = _M0FPB12dec__count32(_M0L3numS345);
      if (_M0L12is__negativeS344) {
        _M0L6_2atmpS2053 = 1;
      } else {
        _M0L6_2atmpS2053 = 0;
      }
      _M0L10total__lenS348 = _M0L10digit__lenS347 + _M0L6_2atmpS2053;
      _M0L6bufferS349
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS348, 0);
      if (_M0L12is__negativeS344) {
        _M0L12digit__startS350 = 1;
      } else {
        _M0L12digit__startS350 = 0;
      }
      #line 239 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__dec(_M0L6bufferS349, _M0L3numS345, _M0L12digit__startS350, _M0L10total__lenS348);
      _M0L6bufferS346 = _M0L6bufferS349;
      break;
    }
    
    case 16: {
      int32_t _M0L10digit__lenS351;
      int32_t _M0L6_2atmpS2054;
      int32_t _M0L10total__lenS352;
      uint16_t* _M0L6bufferS353;
      int32_t _M0L12digit__startS354;
      #line 243 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS351 = _M0FPB12hex__count32(_M0L3numS345);
      if (_M0L12is__negativeS344) {
        _M0L6_2atmpS2054 = 1;
      } else {
        _M0L6_2atmpS2054 = 0;
      }
      _M0L10total__lenS352 = _M0L10digit__lenS351 + _M0L6_2atmpS2054;
      _M0L6bufferS353
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS352, 0);
      if (_M0L12is__negativeS344) {
        _M0L12digit__startS354 = 1;
      } else {
        _M0L12digit__startS354 = 0;
      }
      #line 247 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__hex(_M0L6bufferS353, _M0L3numS345, _M0L12digit__startS354, _M0L10total__lenS352);
      _M0L6bufferS346 = _M0L6bufferS353;
      break;
    }
    default: {
      int32_t _M0L10digit__lenS355;
      int32_t _M0L6_2atmpS2055;
      int32_t _M0L10total__lenS356;
      uint16_t* _M0L6bufferS357;
      int32_t _M0L12digit__startS358;
      #line 251 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS355
      = _M0FPB14radix__count32(_M0L3numS345, _M0L5radixS342);
      if (_M0L12is__negativeS344) {
        _M0L6_2atmpS2055 = 1;
      } else {
        _M0L6_2atmpS2055 = 0;
      }
      _M0L10total__lenS356 = _M0L10digit__lenS355 + _M0L6_2atmpS2055;
      _M0L6bufferS357
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS356, 0);
      if (_M0L12is__negativeS344) {
        _M0L12digit__startS358 = 1;
      } else {
        _M0L12digit__startS358 = 0;
      }
      #line 255 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB24int__to__string__generic(_M0L6bufferS357, _M0L3numS345, _M0L12digit__startS358, _M0L10total__lenS356, _M0L5radixS342);
      _M0L6bufferS346 = _M0L6bufferS357;
      break;
    }
  }
  if (_M0L12is__negativeS344) {
    _M0L6bufferS346[0] = 45;
  }
  return _M0L6bufferS346;
}

int32_t _M0FPB14radix__count32(
  uint32_t _M0L5valueS336,
  int32_t _M0L5radixS338
) {
  uint32_t _M0L4baseS337;
  uint32_t _M0L3numS339;
  int32_t _M0L5countS340;
  #line 189 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS336 == 0u) {
    return 1;
  }
  _M0L4baseS337 = *(uint32_t*)&_M0L5radixS338;
  _M0L3numS339 = _M0L5valueS336;
  _M0L5countS340 = 0;
  while (1) {
    if (_M0L3numS339 > 0u) {
      uint32_t _M0L6_2atmpS2051 = _M0L3numS339 / _M0L4baseS337;
      int32_t _M0L6_2atmpS2052 = _M0L5countS340 + 1;
      _M0L3numS339 = _M0L6_2atmpS2051;
      _M0L5countS340 = _M0L6_2atmpS2052;
      continue;
    } else {
      return _M0L5countS340;
    }
    break;
  }
}

int32_t _M0FPB12hex__count32(uint32_t _M0L5valueS334) {
  #line 177 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS334 == 0u) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS335;
    int32_t _M0L6_2atmpS2050;
    int32_t _M0L6_2atmpS2049;
    #line 182 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS335 = moonbit_clz32(_M0L5valueS334);
    _M0L6_2atmpS2050 = 31 - _M0L14leading__zerosS335;
    _M0L6_2atmpS2049 = _M0L6_2atmpS2050 / 4;
    return _M0L6_2atmpS2049 + 1;
  }
}

int32_t _M0FPB12dec__count32(uint32_t _M0L5valueS333) {
  #line 143 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS333 >= 100000u) {
    if (_M0L5valueS333 >= 10000000u) {
      if (_M0L5valueS333 >= 1000000000u) {
        return 10;
      } else if (_M0L5valueS333 >= 100000000u) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS333 >= 1000000u) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS333 >= 1000u) {
    if (_M0L5valueS333 >= 10000u) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS333 >= 100u) {
    return 3;
  } else if (_M0L5valueS333 >= 10u) {
    return 2;
  } else {
    return 1;
  }
}

int32_t _M0FPB20int__to__string__dec(
  uint16_t* _M0L6bufferS319,
  uint32_t _M0L3numS331,
  int32_t _M0L12digit__startS320,
  int32_t _M0L10total__lenS332
) {
  int32_t _M0L6_2atmpS2048;
  uint32_t _M0L3numS309;
  int32_t _M0L6offsetS310;
  #line 88 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS2048 = _M0L10total__lenS332 - _M0L12digit__startS320;
  _M0L3numS309 = _M0L3numS331;
  _M0L6offsetS310 = _M0L6_2atmpS2048;
  while (1) {
    if (_M0L3numS309 >= 10000u) {
      uint32_t _M0L1tS311 = _M0L3numS309 / 10000u;
      uint32_t _M0L6_2atmpS2025 = _M0L3numS309 % 10000u;
      int32_t _M0L1rS312 = *(int32_t*)&_M0L6_2atmpS2025;
      int32_t _M0L2d1S313 = _M0L1rS312 / 100;
      int32_t _M0L2d2S314 = _M0L1rS312 % 100;
      int32_t _M0L6_2atmpS2024 = _M0L2d1S313 / 10;
      int32_t _M0L6_2atmpS2023 = 48 + _M0L6_2atmpS2024;
      int32_t _M0L6d1__hiS315 = (uint16_t)_M0L6_2atmpS2023;
      int32_t _M0L6_2atmpS2022 = _M0L2d1S313 % 10;
      int32_t _M0L6_2atmpS2021 = 48 + _M0L6_2atmpS2022;
      int32_t _M0L6d1__loS316 = (uint16_t)_M0L6_2atmpS2021;
      int32_t _M0L6_2atmpS2020 = _M0L2d2S314 / 10;
      int32_t _M0L6_2atmpS2019 = 48 + _M0L6_2atmpS2020;
      int32_t _M0L6d2__hiS317 = (uint16_t)_M0L6_2atmpS2019;
      int32_t _M0L6_2atmpS2018 = _M0L2d2S314 % 10;
      int32_t _M0L6_2atmpS2017 = 48 + _M0L6_2atmpS2018;
      int32_t _M0L6d2__loS318 = (uint16_t)_M0L6_2atmpS2017;
      int32_t _M0L6_2atmpS2009 = _M0L12digit__startS320 + _M0L6offsetS310;
      int32_t _M0L6_2atmpS2008 = _M0L6_2atmpS2009 - 4;
      int32_t _M0L6_2atmpS2011;
      int32_t _M0L6_2atmpS2010;
      int32_t _M0L6_2atmpS2013;
      int32_t _M0L6_2atmpS2012;
      int32_t _M0L6_2atmpS2015;
      int32_t _M0L6_2atmpS2014;
      int32_t _M0L6_2atmpS2016;
      _M0L6bufferS319[_M0L6_2atmpS2008] = _M0L6d1__hiS315;
      _M0L6_2atmpS2011 = _M0L12digit__startS320 + _M0L6offsetS310;
      _M0L6_2atmpS2010 = _M0L6_2atmpS2011 - 3;
      _M0L6bufferS319[_M0L6_2atmpS2010] = _M0L6d1__loS316;
      _M0L6_2atmpS2013 = _M0L12digit__startS320 + _M0L6offsetS310;
      _M0L6_2atmpS2012 = _M0L6_2atmpS2013 - 2;
      _M0L6bufferS319[_M0L6_2atmpS2012] = _M0L6d2__hiS317;
      _M0L6_2atmpS2015 = _M0L12digit__startS320 + _M0L6offsetS310;
      _M0L6_2atmpS2014 = _M0L6_2atmpS2015 - 1;
      _M0L6bufferS319[_M0L6_2atmpS2014] = _M0L6d2__loS318;
      _M0L6_2atmpS2016 = _M0L6offsetS310 - 4;
      _M0L3numS309 = _M0L1tS311;
      _M0L6offsetS310 = _M0L6_2atmpS2016;
      continue;
    } else {
      int32_t _M0L6_2atmpS2047 = *(int32_t*)&_M0L3numS309;
      int32_t _M0L9remainingS322 = _M0L6_2atmpS2047;
      int32_t _M0L6offsetS323 = _M0L6offsetS310;
      while (1) {
        if (_M0L9remainingS322 >= 100) {
          int32_t _M0L1tS324 = _M0L9remainingS322 / 100;
          int32_t _M0L1dS325 = _M0L9remainingS322 % 100;
          int32_t _M0L6_2atmpS2034 = _M0L1dS325 / 10;
          int32_t _M0L6_2atmpS2033 = 48 + _M0L6_2atmpS2034;
          int32_t _M0L5d__hiS326 = (uint16_t)_M0L6_2atmpS2033;
          int32_t _M0L6_2atmpS2032 = _M0L1dS325 % 10;
          int32_t _M0L6_2atmpS2031 = 48 + _M0L6_2atmpS2032;
          int32_t _M0L5d__loS327 = (uint16_t)_M0L6_2atmpS2031;
          int32_t _M0L6_2atmpS2027 = _M0L12digit__startS320 + _M0L6offsetS323;
          int32_t _M0L6_2atmpS2026 = _M0L6_2atmpS2027 - 2;
          int32_t _M0L6_2atmpS2029;
          int32_t _M0L6_2atmpS2028;
          int32_t _M0L6_2atmpS2030;
          _M0L6bufferS319[_M0L6_2atmpS2026] = _M0L5d__hiS326;
          _M0L6_2atmpS2029 = _M0L12digit__startS320 + _M0L6offsetS323;
          _M0L6_2atmpS2028 = _M0L6_2atmpS2029 - 1;
          _M0L6bufferS319[_M0L6_2atmpS2028] = _M0L5d__loS327;
          _M0L6_2atmpS2030 = _M0L6offsetS323 - 2;
          _M0L9remainingS322 = _M0L1tS324;
          _M0L6offsetS323 = _M0L6_2atmpS2030;
          continue;
        } else if (_M0L9remainingS322 >= 10) {
          int32_t _M0L6_2atmpS2042 = _M0L9remainingS322 / 10;
          int32_t _M0L6_2atmpS2041 = 48 + _M0L6_2atmpS2042;
          int32_t _M0L5d__hiS329 = (uint16_t)_M0L6_2atmpS2041;
          int32_t _M0L6_2atmpS2040 = _M0L9remainingS322 % 10;
          int32_t _M0L6_2atmpS2039 = 48 + _M0L6_2atmpS2040;
          int32_t _M0L5d__loS330 = (uint16_t)_M0L6_2atmpS2039;
          int32_t _M0L6_2atmpS2036 = _M0L12digit__startS320 + _M0L6offsetS323;
          int32_t _M0L6_2atmpS2035 = _M0L6_2atmpS2036 - 2;
          int32_t _M0L6_2atmpS2038;
          int32_t _M0L6_2atmpS2037;
          _M0L6bufferS319[_M0L6_2atmpS2035] = _M0L5d__hiS329;
          _M0L6_2atmpS2038 = _M0L12digit__startS320 + _M0L6offsetS323;
          _M0L6_2atmpS2037 = _M0L6_2atmpS2038 - 1;
          _M0L6bufferS319[_M0L6_2atmpS2037] = _M0L5d__loS330;
        } else {
          int32_t _M0L6_2atmpS2046 = _M0L12digit__startS320 + _M0L6offsetS323;
          int32_t _M0L6_2atmpS2043 = _M0L6_2atmpS2046 - 1;
          int32_t _M0L6_2atmpS2045 = 48 + _M0L9remainingS322;
          int32_t _M0L6_2atmpS2044 = (uint16_t)_M0L6_2atmpS2045;
          _M0L6bufferS319[_M0L6_2atmpS2043] = _M0L6_2atmpS2044;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB24int__to__string__generic(
  uint16_t* _M0L6bufferS299,
  uint32_t _M0L3numS303,
  int32_t _M0L12digit__startS300,
  int32_t _M0L10total__lenS302,
  int32_t _M0L5radixS293
) {
  uint32_t _M0L4baseS292;
  int32_t _M0L6_2atmpS1993;
  int32_t _M0L6_2atmpS1992;
  #line 57 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS292 = *(uint32_t*)&_M0L5radixS293;
  _M0L6_2atmpS1993 = _M0L5radixS293 - 1;
  _M0L6_2atmpS1992 = _M0L5radixS293 & _M0L6_2atmpS1993;
  if (_M0L6_2atmpS1992 == 0) {
    int32_t _M0L5shiftS294;
    uint32_t _M0L4maskS295;
    int32_t _M0L6_2atmpS2000;
    int32_t _M0L6offsetS296;
    uint32_t _M0L1nS297;
    #line 68 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS294 = moonbit_ctz32(_M0L5radixS293);
    _M0L4maskS295 = _M0L4baseS292 - 1u;
    _M0L6_2atmpS2000 = _M0L10total__lenS302 - _M0L12digit__startS300;
    _M0L6offsetS296 = _M0L6_2atmpS2000;
    _M0L1nS297 = _M0L3numS303;
    while (1) {
      if (_M0L1nS297 > 0u) {
        uint32_t _M0L6_2atmpS1999 = _M0L1nS297 & _M0L4maskS295;
        int32_t _M0L5digitS298 = *(int32_t*)&_M0L6_2atmpS1999;
        int32_t _M0L6_2atmpS1996 = _M0L12digit__startS300 + _M0L6offsetS296;
        int32_t _M0L6_2atmpS1994 = _M0L6_2atmpS1996 - 1;
        int32_t _M0L6_2atmpS1995 =
          ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L5digitS298];
        int32_t _M0L6_2atmpS1997;
        uint32_t _M0L6_2atmpS1998;
        _M0L6bufferS299[_M0L6_2atmpS1994] = _M0L6_2atmpS1995;
        _M0L6_2atmpS1997 = _M0L6offsetS296 - 1;
        _M0L6_2atmpS1998 = _M0L1nS297 >> (_M0L5shiftS294 & 31);
        _M0L6offsetS296 = _M0L6_2atmpS1997;
        _M0L1nS297 = _M0L6_2atmpS1998;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS2007 = _M0L10total__lenS302 - _M0L12digit__startS300;
    int32_t _M0L6offsetS304 = _M0L6_2atmpS2007;
    uint32_t _M0L1nS305 = _M0L3numS303;
    while (1) {
      if (_M0L1nS305 > 0u) {
        uint32_t _M0L1qS306 = _M0L1nS305 / _M0L4baseS292;
        uint32_t _M0L6_2atmpS2006 = _M0L1qS306 * _M0L4baseS292;
        uint32_t _M0L6_2atmpS2005 = _M0L1nS305 - _M0L6_2atmpS2006;
        int32_t _M0L5digitS307 = *(int32_t*)&_M0L6_2atmpS2005;
        int32_t _M0L6_2atmpS2003 = _M0L12digit__startS300 + _M0L6offsetS304;
        int32_t _M0L6_2atmpS2001 = _M0L6_2atmpS2003 - 1;
        int32_t _M0L6_2atmpS2002 =
          ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L5digitS307];
        int32_t _M0L6_2atmpS2004;
        _M0L6bufferS299[_M0L6_2atmpS2001] = _M0L6_2atmpS2002;
        _M0L6_2atmpS2004 = _M0L6offsetS304 - 1;
        _M0L6offsetS304 = _M0L6_2atmpS2004;
        _M0L1nS305 = _M0L1qS306;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB20int__to__string__hex(
  uint16_t* _M0L6bufferS286,
  uint32_t _M0L3numS291,
  int32_t _M0L12digit__startS287,
  int32_t _M0L10total__lenS290
) {
  int32_t _M0L6_2atmpS1991;
  int32_t _M0L6offsetS281;
  uint32_t _M0L1nS282;
  #line 29 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1991 = _M0L10total__lenS290 - _M0L12digit__startS287;
  _M0L6offsetS281 = _M0L6_2atmpS1991;
  _M0L1nS282 = _M0L3numS291;
  while (1) {
    if (_M0L6offsetS281 >= 2) {
      uint32_t _M0L6_2atmpS1988 = _M0L1nS282 & 255u;
      int32_t _M0L9byte__valS283 = *(int32_t*)&_M0L6_2atmpS1988;
      int32_t _M0L2hiS284 = _M0L9byte__valS283 / 16;
      int32_t _M0L2loS285 = _M0L9byte__valS283 % 16;
      int32_t _M0L6_2atmpS1982 = _M0L12digit__startS287 + _M0L6offsetS281;
      int32_t _M0L6_2atmpS1980 = _M0L6_2atmpS1982 - 2;
      int32_t _M0L6_2atmpS1981 =
        ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L2hiS284];
      int32_t _M0L6_2atmpS1985;
      int32_t _M0L6_2atmpS1983;
      int32_t _M0L6_2atmpS1984;
      int32_t _M0L6_2atmpS1986;
      uint32_t _M0L6_2atmpS1987;
      _M0L6bufferS286[_M0L6_2atmpS1980] = _M0L6_2atmpS1981;
      _M0L6_2atmpS1985 = _M0L12digit__startS287 + _M0L6offsetS281;
      _M0L6_2atmpS1983 = _M0L6_2atmpS1985 - 1;
      _M0L6_2atmpS1984
      = ((moonbit_string_t)moonbit_string_literal_60.data)[
        _M0L2loS285
      ];
      _M0L6bufferS286[_M0L6_2atmpS1983] = _M0L6_2atmpS1984;
      _M0L6_2atmpS1986 = _M0L6offsetS281 - 2;
      _M0L6_2atmpS1987 = _M0L1nS282 >> 8;
      _M0L6offsetS281 = _M0L6_2atmpS1986;
      _M0L1nS282 = _M0L6_2atmpS1987;
      continue;
    } else if (_M0L6offsetS281 == 1) {
      uint32_t _M0L6_2atmpS1990 = _M0L1nS282 & 15u;
      int32_t _M0L6nibbleS289 = *(int32_t*)&_M0L6_2atmpS1990;
      int32_t _M0L6_2atmpS1989 =
        ((moonbit_string_t)moonbit_string_literal_60.data)[_M0L6nibbleS289];
      _M0L6bufferS286[_M0L12digit__startS287] = _M0L6_2atmpS1989;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPB4Iter4nextGsE(struct _M0TPB4IterGsE* _M0L4selfS276) {
  struct _M0TWEOs* _M0L7_2afuncS275;
  moonbit_string_t _M0L6resultS277;
  int64_t _M0L7_2abindS278;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS275 = _M0L4selfS276->$0;
  moonbit_incref(_M0L7_2afuncS275);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS277 = _M0L7_2afuncS275->code(_M0L7_2afuncS275);
  moonbit_decref(_M0L7_2afuncS275);
  _M0L7_2abindS278 = _M0L4selfS276->$1;
  if (_M0L6resultS277 == 0) {
    _M0L4selfS276->$1 = _M0MPB4Iter4nextN6constrS9981GsE;
  } else if (_M0L7_2abindS278 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS279 = _M0L7_2abindS278;
    int32_t _M0L4_2anS280 = (int32_t)_M0L7_2aSomeS279;
    int64_t _M0L6_2atmpS1978;
    if (_M0L4_2anS280 > 0) {
      int32_t _M0L6_2atmpS1979 = _M0L4_2anS280 - 1;
      _M0L6_2atmpS1978 = (int64_t)_M0L6_2atmpS1979;
    } else {
      _M0L6_2atmpS1978 = _M0MPB4Iter4nextN6constrS9980GsE;
    }
    _M0L4selfS276->$1 = _M0L6_2atmpS1978;
  }
  return _M0L6resultS277;
}

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void* _M0L4selfS274
) {
  struct _M0TPB13StringBuilder* _M0L6loggerS273;
  struct _M0TPB6Logger _M0L6_2atmpS1977;
  moonbit_string_t _result_4049;
  #line 165 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 166 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS273 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  moonbit_incref(_M0L6loggerS273);
  _M0L6_2atmpS1977
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L6loggerS273
  };
  #line 167 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB7FailurePB4Show6output(_M0L4selfS274, _M0L6_2atmpS1977);
  if (_M0L6_2atmpS1977.$1) {
    moonbit_decref(_M0L6_2atmpS1977.$1);
  }
  #line 168 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _result_4049 = _M0MPB13StringBuilder10to__string(_M0L6loggerS273);
  moonbit_decref(_M0L6loggerS273);
  return _result_4049;
}

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t _M0L4selfS266,
  struct _M0TPB6Logger _M0L6loggerS265
) {
  moonbit_string_t _M0L6_2atmpS1973;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1973 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS266);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS265.$0->$method_0(_M0L6loggerS265.$1, _M0L6_2atmpS1973);
  moonbit_decref(_M0L6_2atmpS1973);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t _M0L4selfS268,
  struct _M0TPB6Logger _M0L6loggerS267
) {
  moonbit_string_t _M0L6_2atmpS1974;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1974 = _M0IPC13int3IntPB4Show10to__string(_M0L4selfS268);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS267.$0->$method_0(_M0L6loggerS267.$1, _M0L6_2atmpS1974);
  moonbit_decref(_M0L6_2atmpS1974);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGbE(
  int32_t _M0L4selfS270,
  struct _M0TPB6Logger _M0L6loggerS269
) {
  moonbit_string_t _M0L6_2atmpS1975;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1975 = _M0IPC14bool4BoolPB4Show10to__string(_M0L4selfS270);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS269.$0->$method_0(_M0L6loggerS269.$1, _M0L6_2atmpS1975);
  moonbit_decref(_M0L6_2atmpS1975);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGmE(
  uint64_t _M0L4selfS272,
  struct _M0TPB6Logger _M0L6loggerS271
) {
  moonbit_string_t _M0L6_2atmpS1976;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1976 = _M0IPC16uint646UInt64PB4Show10to__string(_M0L4selfS272);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS271.$0->$method_0(_M0L6loggerS271.$1, _M0L6_2atmpS1976);
  moonbit_decref(_M0L6_2atmpS1976);
  return 0;
}

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView _M0L4selfS264
) {
  #line 99 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  return _M0L4selfS264.$1;
}

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView _M0L4selfS263
) {
  moonbit_string_t _M0L8_2afieldS3665;
  #line 92 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L8_2afieldS3665 = _M0L4selfS263.$0;
  moonbit_incref(_M0L8_2afieldS3665);
  return _M0L8_2afieldS3665;
}

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS259,
  moonbit_string_t _M0L5valueS260,
  int32_t _M0L5startS261,
  int32_t _M0L3lenS262
) {
  int32_t _M0L6_2atmpS1972;
  int64_t _M0L6_2atmpS1971;
  struct _M0TPC16string10StringView _M0L6_2atmpS1970;
  #line 122 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1972 = _M0L5startS261 + _M0L3lenS262;
  _M0L6_2atmpS1971 = (int64_t)_M0L6_2atmpS1972;
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1970
  = _M0MPC16string6String11sub_2einner(_M0L5valueS260, _M0L5startS261, _M0L6_2atmpS1971);
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L4selfS259, _M0L6_2atmpS1970);
  moonbit_decref(_M0L6_2atmpS1970.$0);
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t _M0L4selfS252,
  int32_t _M0L5startS258,
  int64_t _M0L3endS254
) {
  int32_t _M0L3lenS251;
  int32_t _M0L3endS253;
  int32_t _M0L5startS257;
  int32_t _if__result_4050;
  #line 755 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3lenS251 = Moonbit_array_length(_M0L4selfS252);
  if (_M0L3endS254 == 4294967296ll) {
    _M0L3endS253 = _M0L3lenS251;
  } else {
    int64_t _M0L7_2aSomeS255 = _M0L3endS254;
    int32_t _M0L6_2aendS256 = (int32_t)_M0L7_2aSomeS255;
    if (_M0L6_2aendS256 < 0) {
      _M0L3endS253 = _M0L3lenS251 + _M0L6_2aendS256;
    } else {
      _M0L3endS253 = _M0L6_2aendS256;
    }
  }
  if (_M0L5startS258 < 0) {
    _M0L5startS257 = _M0L3lenS251 + _M0L5startS258;
  } else {
    _M0L5startS257 = _M0L5startS258;
  }
  if (_M0L5startS257 >= 0) {
    if (_M0L5startS257 <= _M0L3endS253) {
      _if__result_4050 = _M0L3endS253 <= _M0L3lenS251;
    } else {
      _if__result_4050 = 0;
    }
  } else {
    _if__result_4050 = 0;
  }
  if (_if__result_4050) {
    if (_M0L5startS257 < _M0L3lenS251) {
      int32_t _M0L6_2atmpS1967 = _M0L4selfS252[_M0L5startS257];
      int32_t _M0L6_2atmpS1966;
      #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1966
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1967);
      if (!_M0L6_2atmpS1966) {
        
      } else {
        #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L3endS253 < _M0L3lenS251) {
      int32_t _M0L6_2atmpS1969 = _M0L4selfS252[_M0L3endS253];
      int32_t _M0L6_2atmpS1968;
      #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1968
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1969);
      if (!_M0L6_2atmpS1968) {
        
      } else {
        #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    moonbit_incref(_M0L4selfS252);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS252,
                                                 .$1 = _M0L5startS257,
                                                 .$2 = _M0L3endS253};
  } else {
    #line 763 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS250,
  struct _M0TPB4Show _M0L4showS249
) {
  struct _M0TPB6Logger _M0L6_2atmpS1965;
  #line 116 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS250);
  _M0L6_2atmpS1965
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS250
  };
  #line 117 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS249.$0->$method_0(_M0L4showS249.$1, _M0L6_2atmpS1965);
  if (_M0L6_2atmpS1965.$1) {
    moonbit_decref(_M0L6_2atmpS1965.$1);
  }
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS248,
  struct _M0TPB4Show _M0L4showS247
) {
  struct _M0TPB6Logger _M0L6_2atmpS1964;
  #line 111 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS248);
  _M0L6_2atmpS1964
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS248
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS247.$0->$method_0(_M0L4showS247.$1, _M0L6_2atmpS1964);
  if (_M0L6_2atmpS1964.$1) {
    moonbit_decref(_M0L6_2atmpS1964.$1);
  }
  return 0;
}

int32_t _M0FPB13finalize__acc(uint32_t _M0L3accS246) {
  uint32_t _M0L6_2atmpS1963;
  #line 444 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  #line 445 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1963 = _M0FPB14avalanche__acc(_M0L3accS246);
  return *(int32_t*)&_M0L6_2atmpS1963;
}

uint32_t _M0FPB14avalanche__acc(uint32_t _M0L3accS245) {
  uint32_t _M0Lm3accS244;
  uint32_t _M0L6_2atmpS1952;
  uint32_t _M0L6_2atmpS1954;
  uint32_t _M0L6_2atmpS1953;
  uint32_t _M0L6_2atmpS1955;
  uint32_t _M0L6_2atmpS1956;
  uint32_t _M0L6_2atmpS1958;
  uint32_t _M0L6_2atmpS1957;
  uint32_t _M0L6_2atmpS1959;
  uint32_t _M0L6_2atmpS1960;
  uint32_t _M0L6_2atmpS1962;
  uint32_t _M0L6_2atmpS1961;
  #line 449 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0Lm3accS244 = _M0L3accS245;
  _M0L6_2atmpS1952 = _M0Lm3accS244;
  _M0L6_2atmpS1954 = _M0Lm3accS244;
  _M0L6_2atmpS1953 = _M0L6_2atmpS1954 >> 15;
  _M0Lm3accS244 = _M0L6_2atmpS1952 ^ _M0L6_2atmpS1953;
  _M0L6_2atmpS1955 = _M0Lm3accS244;
  _M0Lm3accS244 = _M0L6_2atmpS1955 * 2246822519u;
  _M0L6_2atmpS1956 = _M0Lm3accS244;
  _M0L6_2atmpS1958 = _M0Lm3accS244;
  _M0L6_2atmpS1957 = _M0L6_2atmpS1958 >> 13;
  _M0Lm3accS244 = _M0L6_2atmpS1956 ^ _M0L6_2atmpS1957;
  _M0L6_2atmpS1959 = _M0Lm3accS244;
  _M0Lm3accS244 = _M0L6_2atmpS1959 * 3266489917u;
  _M0L6_2atmpS1960 = _M0Lm3accS244;
  _M0L6_2atmpS1962 = _M0Lm3accS244;
  _M0L6_2atmpS1961 = _M0L6_2atmpS1962 >> 16;
  _M0Lm3accS244 = _M0L6_2atmpS1960 ^ _M0L6_2atmpS1961;
  return _M0Lm3accS244;
}

int32_t _M0IP016_24default__implPB2Eq10not__equalGOsE(
  moonbit_string_t _M0L1xS242,
  moonbit_string_t _M0L1yS243
) {
  int32_t _M0L6_2atmpS1951;
  #line 23 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 24 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1951
  = _M0IPC16option6OptionPB2Eq5equalGsE(_M0L1xS242, _M0L1yS243);
  return !_M0L6_2atmpS1951;
}

uint64_t _M0MPC13int3Int10to__uint64(int32_t _M0L4selfS241) {
  int64_t _M0L6_2atmpS1950;
  #line 907 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1950 = (int64_t)_M0L4selfS241;
  return *(uint64_t*)&_M0L6_2atmpS1950;
}

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder* _M0L4selfS240,
  moonbit_string_t _M0L3strS239
) {
  int32_t _M0L8str__lenS238;
  int32_t _M0L3lenS1945;
  int32_t _M0L6_2atmpS1944;
  uint16_t* _M0L4dataS1946;
  int32_t _M0L3lenS1947;
  int32_t _M0L3lenS1949;
  int32_t _M0L6_2atmpS1948;
  #line 86 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L8str__lenS238 = Moonbit_array_length(_M0L3strS239);
  _M0L3lenS1945 = _M0L4selfS240->$1;
  _M0L6_2atmpS1944 = _M0L3lenS1945 + _M0L8str__lenS238;
  #line 88 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS240, _M0L6_2atmpS1944);
  _M0L4dataS1946 = _M0L4selfS240->$0;
  _M0L3lenS1947 = _M0L4selfS240->$1;
  moonbit_incref(_M0L4dataS1946);
  #line 89 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1946, _M0L3lenS1947, _M0L3strS239, 0, _M0L8str__lenS238);
  moonbit_decref(_M0L4dataS1946);
  _M0L3lenS1949 = _M0L4selfS240->$1;
  _M0L6_2atmpS1948 = _M0L3lenS1949 + _M0L8str__lenS238;
  _M0L4selfS240->$1 = _M0L6_2atmpS1948;
  return 0;
}

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t* _M0L4selfS234,
  int32_t _M0L11dst__offsetS237,
  moonbit_string_t _M0L3strS235,
  int32_t _M0L11str__offsetS230,
  int32_t _M0L3lenS231
) {
  int32_t _M0L16end__str__offsetS229;
  int32_t _M0L1iS232;
  int32_t _M0L1jS233;
  #line 71 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L16end__str__offsetS229 = _M0L11str__offsetS230 + _M0L3lenS231;
  _M0L1iS232 = _M0L11str__offsetS230;
  _M0L1jS233 = _M0L11dst__offsetS237;
  while (1) {
    if (_M0L1iS232 < _M0L16end__str__offsetS229) {
      int32_t _M0L6_2atmpS1941 = _M0L3strS235[_M0L1iS232];
      int32_t _M0L6_2atmpS1942;
      int32_t _M0L6_2atmpS1943;
      if (
        _M0L1jS233 < 0 || _M0L1jS233 >= Moonbit_array_length(_M0L4selfS234)
      ) {
        #line 80 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
        moonbit_panic();
      }
      _M0L4selfS234[_M0L1jS233] = _M0L6_2atmpS1941;
      _M0L6_2atmpS1942 = _M0L1iS232 + 1;
      _M0L6_2atmpS1943 = _M0L1jS233 + 1;
      _M0L1iS232 = _M0L6_2atmpS1942;
      _M0L1jS233 = _M0L6_2atmpS1943;
      continue;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t _M0L4selfS227,
  int32_t _M0L5quoteS228
) {
  struct _M0TPB13StringBuilder* _M0L3bufS226;
  int32_t _M0L6_2atmpS1940;
  struct _M0TPC16string10StringView _M0L6_2atmpS1938;
  struct _M0TPB6Logger _M0L6_2atmpS1939;
  moonbit_string_t _result_4052;
  #line 110 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 111 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L3bufS226 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS1940 = Moonbit_array_length(_M0L4selfS227);
  moonbit_incref(_M0L4selfS227);
  _M0L6_2atmpS1938
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS227, .$1 = 0, .$2 = _M0L6_2atmpS1940
  };
  moonbit_incref(_M0L3bufS226);
  _M0L6_2atmpS1939
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS226
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0MPC16string10StringView18escape__to_2einner(_M0L6_2atmpS1938, _M0L6_2atmpS1939, _M0L5quoteS228);
  moonbit_decref(_M0L6_2atmpS1938.$0);
  if (_M0L6_2atmpS1939.$1) {
    moonbit_decref(_M0L6_2atmpS1939.$1);
  }
  #line 113 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_4052 = _M0MPB13StringBuilder10to__string(_M0L3bufS226);
  moonbit_decref(_M0L3bufS226);
  return _result_4052;
}

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView _M0L4selfS218,
  struct _M0TPB6Logger _M0L6loggerS216,
  int32_t _M0L5quoteS215
) {
  int32_t _M0L3endS1936;
  int32_t _M0L5startS1937;
  int32_t _M0L3lenS217;
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS219;
  int32_t _M0L1iS220;
  int32_t _M0L3segS221;
  #line 144 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L5quoteS215) {
    #line 150 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS216.$0->$method_3(_M0L6loggerS216.$1, 34);
  }
  _M0L3endS1936 = _M0L4selfS218.$2;
  _M0L5startS1937 = _M0L4selfS218.$1;
  _M0L3lenS217 = _M0L3endS1936 - _M0L5startS1937;
  moonbit_incref(_M0L4selfS218.$0);
  if (_M0L6loggerS216.$1) {
    moonbit_incref(_M0L6loggerS216.$1);
  }
  _M0L6_2aenvS219
  = (struct _M0TURPC16string10StringViewRPB6LoggerE*)moonbit_malloc(sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE));
  Moonbit_object_header(_M0L6_2aenvS219)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 97, 0);
  _M0L6_2aenvS219->$0 = _M0L4selfS218;
  _M0L6_2aenvS219->$1 = _M0L6loggerS216;
  _M0L1iS220 = 0;
  _M0L3segS221 = 0;
  _2afor_222:;
  while (1) {
    moonbit_string_t _M0L3strS1933;
    int32_t _M0L5startS1935;
    int32_t _M0L6_2atmpS1934;
    int32_t _M0L4codeS223;
    int32_t _M0L1cS225;
    int32_t _M0L6_2atmpS1917;
    int32_t _M0L6_2atmpS1918;
    int32_t _M0L6_2atmpS1919;
    int32_t _tmp_4056;
    int32_t _tmp_4057;
    if (_M0L1iS220 >= _M0L3lenS217) {
      #line 160 "/home/developer/.moon/lib/core/builtin/show.mbt"
      _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
      moonbit_decref(_M0L6_2aenvS219);
      break;
    }
    _M0L3strS1933 = _M0L4selfS218.$0;
    _M0L5startS1935 = _M0L4selfS218.$1;
    _M0L6_2atmpS1934 = _M0L5startS1935 + _M0L1iS220;
    _M0L4codeS223 = _M0L3strS1933[_M0L6_2atmpS1934];
    switch (_M0L4codeS223) {
      case 34: {
        _M0L1cS225 = _M0L4codeS223;
        goto join_224;
        break;
      }
      
      case 92: {
        _M0L1cS225 = _M0L4codeS223;
        goto join_224;
        break;
      }
      
      case 10: {
        int32_t _M0L6_2atmpS1920;
        int32_t _M0L6_2atmpS1921;
        #line 172 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
        #line 173 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS216.$0->$method_0(_M0L6loggerS216.$1, (moonbit_string_t)moonbit_string_literal_61.data);
        _M0L6_2atmpS1920 = _M0L1iS220 + 1;
        _M0L6_2atmpS1921 = _M0L1iS220 + 1;
        _M0L1iS220 = _M0L6_2atmpS1920;
        _M0L3segS221 = _M0L6_2atmpS1921;
        goto _2afor_222;
        break;
      }
      
      case 13: {
        int32_t _M0L6_2atmpS1922;
        int32_t _M0L6_2atmpS1923;
        #line 177 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
        #line 178 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS216.$0->$method_0(_M0L6loggerS216.$1, (moonbit_string_t)moonbit_string_literal_62.data);
        _M0L6_2atmpS1922 = _M0L1iS220 + 1;
        _M0L6_2atmpS1923 = _M0L1iS220 + 1;
        _M0L1iS220 = _M0L6_2atmpS1922;
        _M0L3segS221 = _M0L6_2atmpS1923;
        goto _2afor_222;
        break;
      }
      
      case 8: {
        int32_t _M0L6_2atmpS1924;
        int32_t _M0L6_2atmpS1925;
        #line 182 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
        #line 183 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS216.$0->$method_0(_M0L6loggerS216.$1, (moonbit_string_t)moonbit_string_literal_63.data);
        _M0L6_2atmpS1924 = _M0L1iS220 + 1;
        _M0L6_2atmpS1925 = _M0L1iS220 + 1;
        _M0L1iS220 = _M0L6_2atmpS1924;
        _M0L3segS221 = _M0L6_2atmpS1925;
        goto _2afor_222;
        break;
      }
      
      case 9: {
        int32_t _M0L6_2atmpS1926;
        int32_t _M0L6_2atmpS1927;
        #line 187 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
        #line 188 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS216.$0->$method_0(_M0L6loggerS216.$1, (moonbit_string_t)moonbit_string_literal_64.data);
        _M0L6_2atmpS1926 = _M0L1iS220 + 1;
        _M0L6_2atmpS1927 = _M0L1iS220 + 1;
        _M0L1iS220 = _M0L6_2atmpS1926;
        _M0L3segS221 = _M0L6_2atmpS1927;
        goto _2afor_222;
        break;
      }
      default: {
        if (_M0L4codeS223 < 32) {
          int32_t _M0L6_2atmpS1929;
          moonbit_string_t _M0L6_2atmpS1928;
          int32_t _M0L6_2atmpS1930;
          int32_t _M0L6_2atmpS1931;
          #line 193 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
          #line 194 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS216.$0->$method_0(_M0L6loggerS216.$1, (moonbit_string_t)moonbit_string_literal_65.data);
          _M0L6_2atmpS1929 = _M0L4codeS223 & 0xff;
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6_2atmpS1928 = _M0MPC14byte4Byte7to__hex(_M0L6_2atmpS1929);
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS216.$0->$method_0(_M0L6loggerS216.$1, _M0L6_2atmpS1928);
          moonbit_decref(_M0L6_2atmpS1928);
          #line 196 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS216.$0->$method_3(_M0L6loggerS216.$1, 125);
          _M0L6_2atmpS1930 = _M0L1iS220 + 1;
          _M0L6_2atmpS1931 = _M0L1iS220 + 1;
          _M0L1iS220 = _M0L6_2atmpS1930;
          _M0L3segS221 = _M0L6_2atmpS1931;
          goto _2afor_222;
        } else {
          int32_t _M0L6_2atmpS1932 = _M0L1iS220 + 1;
          int32_t _tmp_4055 = _M0L3segS221;
          _M0L1iS220 = _M0L6_2atmpS1932;
          _M0L3segS221 = _tmp_4055;
          goto _2afor_222;
        }
        break;
      }
    }
    goto joinlet_4054;
    join_224:;
    #line 166 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS219, _M0L3segS221, _M0L1iS220);
    #line 167 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS216.$0->$method_3(_M0L6loggerS216.$1, 92);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1917 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS225);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS216.$0->$method_3(_M0L6loggerS216.$1, _M0L6_2atmpS1917);
    _M0L6_2atmpS1918 = _M0L1iS220 + 1;
    _M0L6_2atmpS1919 = _M0L1iS220 + 1;
    _M0L1iS220 = _M0L6_2atmpS1918;
    _M0L3segS221 = _M0L6_2atmpS1919;
    continue;
    joinlet_4054:;
    _tmp_4056 = _M0L1iS220;
    _tmp_4057 = _M0L3segS221;
    _M0L1iS220 = _tmp_4056;
    _M0L3segS221 = _tmp_4057;
    continue;
    break;
  }
  if (_M0L5quoteS215) {
    #line 204 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS216.$0->$method_3(_M0L6loggerS216.$1, 34);
  }
  return 0;
}

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS211,
  int32_t _M0L3segS214,
  int32_t _M0L1iS213
) {
  struct _M0TPB6Logger _M0L6loggerS210;
  struct _M0TPC16string10StringView _M0L4selfS212;
  #line 153 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6loggerS210 = _M0L6_2aenvS211->$1;
  _M0L4selfS212 = _M0L6_2aenvS211->$0;
  if (_M0L1iS213 > _M0L3segS214) {
    int64_t _M0L6_2atmpS1916 = (int64_t)_M0L1iS213;
    struct _M0TPC16string10StringView _M0L6_2atmpS1915;
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1915
    = _M0MPC16string10StringView11sub_2einner(_M0L4selfS212, _M0L3segS214, _M0L6_2atmpS1916);
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS210.$0->$method_2(_M0L6loggerS210.$1, _M0L6_2atmpS1915);
    moonbit_decref(_M0L6_2atmpS1915.$0);
  }
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView _M0L4selfS203,
  int32_t _M0L5startS209,
  int64_t _M0L3endS205
) {
  moonbit_string_t _M0L3strS1914;
  int32_t _M0L8str__lenS202;
  int32_t _M0L8abs__endS204;
  int32_t _M0L10abs__startS208;
  int32_t _M0L5startS1902;
  int32_t _if__result_4058;
  #line 814 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1914 = _M0L4selfS203.$0;
  _M0L8str__lenS202 = Moonbit_array_length(_M0L3strS1914);
  if (_M0L3endS205 == 4294967296ll) {
    _M0L8abs__endS204 = _M0L4selfS203.$2;
  } else {
    int64_t _M0L7_2aSomeS206 = _M0L3endS205;
    int32_t _M0L6_2aendS207 = (int32_t)_M0L7_2aSomeS206;
    if (_M0L6_2aendS207 < 0) {
      int32_t _M0L3endS1912 = _M0L4selfS203.$2;
      _M0L8abs__endS204 = _M0L3endS1912 + _M0L6_2aendS207;
    } else {
      int32_t _M0L5startS1913 = _M0L4selfS203.$1;
      _M0L8abs__endS204 = _M0L5startS1913 + _M0L6_2aendS207;
    }
  }
  if (_M0L5startS209 < 0) {
    int32_t _M0L3endS1910 = _M0L4selfS203.$2;
    _M0L10abs__startS208 = _M0L3endS1910 + _M0L5startS209;
  } else {
    int32_t _M0L5startS1911 = _M0L4selfS203.$1;
    _M0L10abs__startS208 = _M0L5startS1911 + _M0L5startS209;
  }
  _M0L5startS1902 = _M0L4selfS203.$1;
  if (_M0L10abs__startS208 >= _M0L5startS1902) {
    if (_M0L10abs__startS208 <= _M0L8abs__endS204) {
      int32_t _M0L3endS1901 = _M0L4selfS203.$2;
      _if__result_4058 = _M0L8abs__endS204 <= _M0L3endS1901;
    } else {
      _if__result_4058 = 0;
    }
  } else {
    _if__result_4058 = 0;
  }
  if (_if__result_4058) {
    moonbit_string_t _M0L3strS1909;
    if (_M0L10abs__startS208 < _M0L8str__lenS202) {
      moonbit_string_t _M0L3strS1905 = _M0L4selfS203.$0;
      int32_t _M0L6_2atmpS1904 = _M0L3strS1905[_M0L10abs__startS208];
      int32_t _M0L6_2atmpS1903;
      #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1903
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1904);
      if (!_M0L6_2atmpS1903) {
        
      } else {
        #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L8abs__endS204 < _M0L8str__lenS202) {
      moonbit_string_t _M0L3strS1908 = _M0L4selfS203.$0;
      int32_t _M0L6_2atmpS1907 = _M0L3strS1908[_M0L8abs__endS204];
      int32_t _M0L6_2atmpS1906;
      #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1906
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1907);
      if (!_M0L6_2atmpS1906) {
        
      } else {
        #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    _M0L3strS1909 = _M0L4selfS203.$0;
    moonbit_incref(_M0L3strS1909);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1909,
                                                 .$1 = _M0L10abs__startS208,
                                                 .$2 = _M0L8abs__endS204};
  } else {
    #line 834 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t _M0L1bS201) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS200;
  int32_t _M0L6_2atmpS1898;
  int32_t _M0L6_2atmpS1897;
  int32_t _M0L6_2atmpS1900;
  int32_t _M0L6_2atmpS1899;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1896;
  moonbit_string_t _result_4059;
  #line 74 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L7_2aselfS200 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1898 = _M0IPC14byte4BytePB3Div3div(_M0L1bS201, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1897
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1898);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS200, _M0L6_2atmpS1897);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1900 = _M0IPC14byte4BytePB3Mod3mod(_M0L1bS201, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1899
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1900);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS200, _M0L6_2atmpS1899);
  _M0L6_2atmpS1896 = _M0L7_2aselfS200;
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_4059 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1896);
  moonbit_decref(_M0L6_2atmpS1896);
  return _result_4059;
}

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t _M0L1iS199) {
  #line 75 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L1iS199 < 10) {
    int32_t _M0L6_2atmpS1893;
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1893 = _M0IPC14byte4BytePB3Add3add(_M0L1iS199, 48);
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1893);
  } else {
    int32_t _M0L6_2atmpS1895;
    int32_t _M0L6_2atmpS1894;
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1895 = _M0IPC14byte4BytePB3Add3add(_M0L1iS199, 97);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1894 = _M0IPC14byte4BytePB3Sub3sub(_M0L6_2atmpS1895, 10);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1894);
  }
}

int32_t _M0IPC14byte4BytePB3Sub3sub(
  int32_t _M0L4selfS197,
  int32_t _M0L4thatS198
) {
  int32_t _M0L6_2atmpS1891;
  int32_t _M0L6_2atmpS1892;
  int32_t _M0L6_2atmpS1890;
  #line 120 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1891 = (int32_t)_M0L4selfS197;
  _M0L6_2atmpS1892 = (int32_t)_M0L4thatS198;
  _M0L6_2atmpS1890 = _M0L6_2atmpS1891 - _M0L6_2atmpS1892;
  return _M0L6_2atmpS1890 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Mod3mod(
  int32_t _M0L4selfS195,
  int32_t _M0L4thatS196
) {
  int32_t _M0L6_2atmpS1888;
  int32_t _M0L6_2atmpS1889;
  int32_t _M0L6_2atmpS1887;
  #line 67 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1888 = (int32_t)_M0L4selfS195;
  _M0L6_2atmpS1889 = (int32_t)_M0L4thatS196;
  _M0L6_2atmpS1887 = _M0L6_2atmpS1888 % _M0L6_2atmpS1889;
  return _M0L6_2atmpS1887 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Div3div(
  int32_t _M0L4selfS193,
  int32_t _M0L4thatS194
) {
  int32_t _M0L6_2atmpS1885;
  int32_t _M0L6_2atmpS1886;
  int32_t _M0L6_2atmpS1884;
  #line 62 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1885 = (int32_t)_M0L4selfS193;
  _M0L6_2atmpS1886 = (int32_t)_M0L4thatS194;
  _M0L6_2atmpS1884 = _M0L6_2atmpS1885 / _M0L6_2atmpS1886;
  return _M0L6_2atmpS1884 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Add3add(
  int32_t _M0L4selfS191,
  int32_t _M0L4thatS192
) {
  int32_t _M0L6_2atmpS1882;
  int32_t _M0L6_2atmpS1883;
  int32_t _M0L6_2atmpS1881;
  #line 106 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1882 = (int32_t)_M0L4selfS191;
  _M0L6_2atmpS1883 = (int32_t)_M0L4thatS192;
  _M0L6_2atmpS1881 = _M0L6_2atmpS1882 + _M0L6_2atmpS1883;
  return _M0L6_2atmpS1881 & 0xff;
}

int32_t _M0MPC16string6String16unsafe__char__at(
  moonbit_string_t _M0L4selfS188,
  int32_t _M0L5indexS189
) {
  int32_t _M0L2c1S187;
  #line 102 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
  _M0L2c1S187 = _M0L4selfS188[_M0L5indexS189];
  #line 105 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
  if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S187)) {
    int32_t _M0L6_2atmpS1880 = _M0L5indexS189 + 1;
    int32_t _M0L2c2S190 = _M0L4selfS188[_M0L6_2atmpS1880];
    int32_t _M0L6_2atmpS1878 = (int32_t)_M0L2c1S187;
    int32_t _M0L6_2atmpS1879 = (int32_t)_M0L2c2S190;
    #line 107 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
    return _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1878, _M0L6_2atmpS1879);
  } else {
    #line 109 "/home/developer/.moon/lib/core/builtin/deprecated.mbt"
    return _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S187);
  }
}

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t _M0L4selfS186) {
  int32_t _M0L6_2atmpS1877;
  #line 68 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  _M0L6_2atmpS1877 = (int32_t)_M0L4selfS186;
  return _M0L6_2atmpS1877;
}

int32_t _M0FPB32code__point__of__surrogate__pair(
  int32_t _M0L7leadingS184,
  int32_t _M0L8trailingS185
) {
  int32_t _M0L6_2atmpS1876;
  int32_t _M0L6_2atmpS1875;
  int32_t _M0L6_2atmpS1874;
  int32_t _M0L6_2atmpS1873;
  int32_t _M0L6_2atmpS1872;
  #line 40 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1876 = _M0L7leadingS184 - 55296;
  _M0L6_2atmpS1875 = _M0L6_2atmpS1876 * 1024;
  _M0L6_2atmpS1874 = _M0L6_2atmpS1875 + _M0L8trailingS185;
  _M0L6_2atmpS1873 = _M0L6_2atmpS1874 - 56320;
  _M0L6_2atmpS1872 = _M0L6_2atmpS1873 + 65536;
  return _M0L6_2atmpS1872;
}

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t _M0L4selfS183) {
  #line 45 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS183 >= 56320) {
    return _M0L4selfS183 <= 57343;
  } else {
    return 0;
  }
}

int32_t _M0MPC16uint166UInt1622is__leading__surrogate(int32_t _M0L4selfS182) {
  #line 28 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS182 >= 55296) {
    return _M0L4selfS182 <= 56319;
  } else {
    return 0;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder* _M0L4selfS180,
  int32_t _M0L2chS179
) {
  uint32_t _M0L4codeS178;
  #line 95 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  #line 96 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4codeS178 = _M0MPC14char4Char8to__uint(_M0L2chS179);
  if (_M0L4codeS178 <= 65535u) {
    int32_t _M0L3lenS1851 = _M0L4selfS180->$1;
    int32_t _M0L6_2atmpS1850 = _M0L3lenS1851 + 1;
    uint16_t* _M0L4dataS1852;
    int32_t _M0L3lenS1853;
    int32_t _M0L6_2atmpS1854;
    int32_t _M0L3lenS1856;
    int32_t _M0L6_2atmpS1855;
    #line 98 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS180, _M0L6_2atmpS1850);
    _M0L4dataS1852 = _M0L4selfS180->$0;
    _M0L3lenS1853 = _M0L4selfS180->$1;
    moonbit_incref(_M0L4dataS1852);
    #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1854 = _M0MPC14uint4UInt10to__uint16(_M0L4codeS178);
    if (
      _M0L3lenS1853 < 0
      || _M0L3lenS1853 >= Moonbit_array_length(_M0L4dataS1852)
    ) {
      #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1852[_M0L3lenS1853] = _M0L6_2atmpS1854;
    moonbit_decref(_M0L4dataS1852);
    _M0L3lenS1856 = _M0L4selfS180->$1;
    _M0L6_2atmpS1855 = _M0L3lenS1856 + 1;
    _M0L4selfS180->$1 = _M0L6_2atmpS1855;
  } else if (_M0L4codeS178 <= 1114111u) {
    int32_t _M0L3lenS1858 = _M0L4selfS180->$1;
    int32_t _M0L6_2atmpS1857 = _M0L3lenS1858 + 2;
    uint32_t _M0L4codeS181;
    uint16_t* _M0L4dataS1859;
    int32_t _M0L3lenS1860;
    uint32_t _M0L6_2atmpS1863;
    uint32_t _M0L6_2atmpS1862;
    int32_t _M0L6_2atmpS1861;
    uint16_t* _M0L4dataS1864;
    int32_t _M0L3lenS1869;
    int32_t _M0L6_2atmpS1865;
    uint32_t _M0L6_2atmpS1868;
    uint32_t _M0L6_2atmpS1867;
    int32_t _M0L6_2atmpS1866;
    int32_t _M0L3lenS1871;
    int32_t _M0L6_2atmpS1870;
    #line 102 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS180, _M0L6_2atmpS1857);
    _M0L4codeS181 = _M0L4codeS178 - 65536u;
    _M0L4dataS1859 = _M0L4selfS180->$0;
    _M0L3lenS1860 = _M0L4selfS180->$1;
    _M0L6_2atmpS1863 = _M0L4codeS181 >> 10;
    _M0L6_2atmpS1862 = 55296u + _M0L6_2atmpS1863;
    moonbit_incref(_M0L4dataS1859);
    #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1861 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1862);
    if (
      _M0L3lenS1860 < 0
      || _M0L3lenS1860 >= Moonbit_array_length(_M0L4dataS1859)
    ) {
      #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1859[_M0L3lenS1860] = _M0L6_2atmpS1861;
    moonbit_decref(_M0L4dataS1859);
    _M0L4dataS1864 = _M0L4selfS180->$0;
    _M0L3lenS1869 = _M0L4selfS180->$1;
    _M0L6_2atmpS1865 = _M0L3lenS1869 + 1;
    _M0L6_2atmpS1868 = _M0L4codeS181 & 1023u;
    _M0L6_2atmpS1867 = 56320u + _M0L6_2atmpS1868;
    moonbit_incref(_M0L4dataS1864);
    #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1866 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1867);
    if (
      _M0L6_2atmpS1865 < 0
      || _M0L6_2atmpS1865 >= Moonbit_array_length(_M0L4dataS1864)
    ) {
      #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1864[_M0L6_2atmpS1865] = _M0L6_2atmpS1866;
    moonbit_decref(_M0L4dataS1864);
    _M0L3lenS1871 = _M0L4selfS180->$1;
    _M0L6_2atmpS1870 = _M0L3lenS1871 + 2;
    _M0L4selfS180->$1 = _M0L6_2atmpS1870;
  } else {
    #line 108 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_66.data);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder* _M0L4selfS172,
  int32_t _M0L8requiredS173
) {
  uint16_t* _M0L4dataS1849;
  int32_t _M0L12current__lenS171;
  int32_t _M0L13enough__spaceS174;
  int32_t _M0L13enough__spaceS175;
  uint16_t* _M0L4dataS1845;
  int32_t _M0L6_2atmpS1846;
  int32_t _M0L3lenS1847;
  uint16_t* _M0L9new__dataS177;
  uint16_t* _M0L6_2aoldS3675;
  #line 46 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4dataS1849 = _M0L4selfS172->$0;
  _M0L12current__lenS171 = Moonbit_array_length(_M0L4dataS1849);
  if (_M0L8requiredS173 <= _M0L12current__lenS171) {
    return 0;
  }
  _M0L13enough__spaceS175 = _M0L12current__lenS171;
  while (1) {
    if (_M0L13enough__spaceS175 < _M0L8requiredS173) {
      int32_t _M0L6_2atmpS1848 = _M0L13enough__spaceS175 * 2;
      _M0L13enough__spaceS175 = _M0L6_2atmpS1848;
      continue;
    } else {
      _M0L13enough__spaceS174 = _M0L13enough__spaceS175;
    }
    break;
  }
  _M0L4dataS1845 = _M0L4selfS172->$0;
  moonbit_incref(_M0L4dataS1845);
  #line 64 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1846 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L3lenS1847 = _M0L4selfS172->$1;
  #line 61 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L9new__dataS177
  = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1845, _M0L13enough__spaceS174, _M0L6_2atmpS1846, _M0L3lenS1847, 0, 0);
  moonbit_decref(_M0L4dataS1845);
  _M0L6_2aoldS3675 = _M0L4selfS172->$0;
  moonbit_decref(_M0L6_2aoldS3675);
  _M0L4selfS172->$0 = _M0L9new__dataS177;
  return 0;
}

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t _M0L4selfS170) {
  int32_t _M0L6_2atmpS1844;
  #line 2676 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1844 = *(int32_t*)&_M0L4selfS170;
  return (uint16_t)_M0L6_2atmpS1844;
}

uint32_t _M0MPC14char4Char8to__uint(int32_t _M0L4selfS169) {
  int32_t _M0L6_2atmpS1843;
  #line 1254 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1843 = _M0L4selfS169;
  return *(uint32_t*)&_M0L6_2atmpS1843;
}

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder* _M0L4selfS167
) {
  int32_t _M0L3lenS1835;
  uint16_t* _M0L4dataS1837;
  int32_t _M0L6_2atmpS1836;
  #line 148 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3lenS1835 = _M0L4selfS167->$1;
  _M0L4dataS1837 = _M0L4selfS167->$0;
  _M0L6_2atmpS1836 = Moonbit_array_length(_M0L4dataS1837);
  if (_M0L3lenS1835 == _M0L6_2atmpS1836) {
    uint16_t* _M0L4dataS1838 = _M0L4selfS167->$0;
    moonbit_incref(_M0L4dataS1838);
    return _M0L4dataS1838;
  } else {
    uint16_t* _M0L4dataS1839 = _M0L4selfS167->$0;
    int32_t _M0L3lenS1840 = _M0L4selfS167->$1;
    int32_t _M0L6_2atmpS1841;
    int32_t _M0L3lenS1842;
    uint16_t* _M0L4dataS168;
    moonbit_incref(_M0L4dataS1839);
    #line 155 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1841 = _M0IPC16uint166UInt16PB7Default7default();
    _M0L3lenS1842 = _M0L4selfS167->$1;
    #line 152 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L4dataS168
    = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1839, _M0L3lenS1840, _M0L6_2atmpS1841, _M0L3lenS1842, 0, 0);
    moonbit_decref(_M0L4dataS1839);
    return _M0L4dataS168;
  }
}

int32_t _M0IPC16uint166UInt16PB7Default7default() {
  #line 176 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  return 0;
}

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t* _M0L3srcS164,
  int32_t _M0L13allocate__lenS160,
  int32_t _M0L4initS165,
  int32_t _M0L3lenS161,
  int32_t _M0L11src__offsetS162,
  int32_t _M0L11dst__offsetS163
) {
  int32_t _if__result_4061;
  #line 97 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L13allocate__lenS160 >= 0) {
    if (_M0L3lenS161 >= 0) {
      if (_M0L11src__offsetS162 >= 0) {
        if (_M0L11dst__offsetS163 >= 0) {
          int32_t _M0L6_2atmpS1831 = _M0L11src__offsetS162 + _M0L3lenS161;
          int32_t _M0L6_2atmpS1832 = Moonbit_array_length(_M0L3srcS164);
          if (_M0L6_2atmpS1831 <= _M0L6_2atmpS1832) {
            int32_t _M0L6_2atmpS1830 = _M0L11dst__offsetS163 + _M0L3lenS161;
            _if__result_4061 = _M0L6_2atmpS1830 <= _M0L13allocate__lenS160;
          } else {
            _if__result_4061 = 0;
          }
        } else {
          _if__result_4061 = 0;
        }
      } else {
        _if__result_4061 = 0;
      }
    } else {
      _if__result_4061 = 0;
    }
  } else {
    _if__result_4061 = 0;
  }
  if (_if__result_4061) {
    moonbit_incref(_M0L3srcS164);
    #line 115 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    return _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(_M0L3srcS164, _M0L13allocate__lenS160, _M0L4initS165, _M0L11src__offsetS162, _M0L11dst__offsetS163, _M0L3lenS161);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS166;
    int32_t _M0L6_2atmpS1834;
    moonbit_string_t _M0L6_2atmpS1833;
    uint16_t* _result_4062;
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L18_2astring__builderS166
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS166, (moonbit_string_t)moonbit_string_literal_67.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS166, _M0L13allocate__lenS160);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS166, (moonbit_string_t)moonbit_string_literal_68.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS166, _M0L11src__offsetS162);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS166, (moonbit_string_t)moonbit_string_literal_69.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS166, _M0L11dst__offsetS163);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS166, (moonbit_string_t)moonbit_string_literal_70.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS166, _M0L3lenS161);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS166, (moonbit_string_t)moonbit_string_literal_71.data);
    _M0L6_2atmpS1834 = Moonbit_array_length(_M0L3srcS164);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS166, _M0L6_2atmpS1834);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L6_2atmpS1833
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS166);
    moonbit_decref(_M0L18_2astring__builderS166);
    #line 111 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _result_4062 = _M0FPC15abort5abortGAkE(_M0L6_2atmpS1833);
    moonbit_decref(_M0L6_2atmpS1833);
    return _result_4062;
  }
}

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t* _M0L3srcS157,
  int32_t _M0L13allocate__lenS154,
  int32_t _M0L4initS155,
  int32_t _M0L11src__offsetS158,
  int32_t _M0L11dst__offsetS156,
  int32_t _M0L9blit__lenS159
) {
  uint16_t* _M0L3dstS153;
  #line 79 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  _M0L3dstS153
  = (uint16_t*)moonbit_make_string(_M0L13allocate__lenS154, _M0L4initS155);
  moonbit_incref(_M0L3dstS153);
  #line 90 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS153, _M0L11dst__offsetS156, _M0L3srcS157, _M0L11src__offsetS158, _M0L9blit__lenS159, sizeof(uint16_t));
  return _M0L3dstS153;
}

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t _M0L10size__hintS151
) {
  int32_t _M0L7initialS150;
  uint16_t* _M0L4dataS152;
  struct _M0TPB13StringBuilder* _block_4063;
  #line 32 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  if (_M0L10size__hintS151 < 1) {
    _M0L7initialS150 = 1;
  } else {
    int32_t _M0L6_2atmpS1829 = _M0L10size__hintS151 + 1;
    _M0L7initialS150 = _M0L6_2atmpS1829 / 2;
  }
  _M0L4dataS152 = (uint16_t*)moonbit_make_string(_M0L7initialS150, 0);
  _block_4063
  = (struct _M0TPB13StringBuilder*)moonbit_malloc(sizeof(struct _M0TPB13StringBuilder));
  Moonbit_object_header(_block_4063)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 102, 0);
  _block_4063->$0 = _M0L4dataS152;
  _block_4063->$1 = 0;
  return _block_4063;
}

int32_t _M0MPC14byte4Byte8to__char(int32_t _M0L4selfS149) {
  int32_t _M0L6_2atmpS1828;
  #line 1868 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1828 = (int32_t)_M0L4selfS149;
  return _M0L6_2atmpS1828;
}

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t* _M0L3srcS129,
  int32_t _M0L13allocate__lenS125,
  int32_t _M0L3lenS126,
  int32_t _M0L11src__offsetS127,
  int32_t _M0L11dst__offsetS128
) {
  int32_t _if__result_4064;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS125 >= 0) {
    if (_M0L3lenS126 >= 0) {
      if (_M0L11src__offsetS127 >= 0) {
        if (_M0L11dst__offsetS128 >= 0) {
          int32_t _M0L6_2atmpS1809 = _M0L11src__offsetS127 + _M0L3lenS126;
          int32_t _M0L6_2atmpS1810;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1810
          = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS129);
          if (_M0L6_2atmpS1809 <= _M0L6_2atmpS1810) {
            int32_t _M0L6_2atmpS1808 = _M0L11dst__offsetS128 + _M0L3lenS126;
            _if__result_4064 = _M0L6_2atmpS1808 <= _M0L13allocate__lenS125;
          } else {
            _if__result_4064 = 0;
          }
        } else {
          _if__result_4064 = 0;
        }
      } else {
        _if__result_4064 = 0;
      }
    } else {
      _if__result_4064 = 0;
    }
  } else {
    _if__result_4064 = 0;
  }
  if (_if__result_4064) {
    moonbit_incref(_M0L3srcS129);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (moonbit_string_t*)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS125, (moonbit_string_t)moonbit_string_literal_1.data, _M0L3srcS129, _M0L11src__offsetS127, _M0L11dst__offsetS128, _M0L3lenS126);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS130;
    int32_t _M0L6_2atmpS1812;
    moonbit_string_t _M0L6_2atmpS1811;
    moonbit_string_t* _result_4065;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS130
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS130, (moonbit_string_t)moonbit_string_literal_67.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS130, _M0L13allocate__lenS125);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS130, (moonbit_string_t)moonbit_string_literal_68.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS130, _M0L11src__offsetS127);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS130, (moonbit_string_t)moonbit_string_literal_69.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS130, _M0L11dst__offsetS128);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS130, (moonbit_string_t)moonbit_string_literal_70.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS130, _M0L3lenS126);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS130, (moonbit_string_t)moonbit_string_literal_71.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1812 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS129);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS130, _M0L6_2atmpS1812);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1811
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS130);
    moonbit_decref(_M0L18_2astring__builderS130);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4065
    = _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(_M0L6_2atmpS1811);
    moonbit_decref(_M0L6_2atmpS1811);
    return _result_4065;
  }
}

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE** _M0L3srcS135,
  int32_t _M0L13allocate__lenS131,
  int32_t _M0L3lenS132,
  int32_t _M0L11src__offsetS133,
  int32_t _M0L11dst__offsetS134
) {
  int32_t _if__result_4066;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS131 >= 0) {
    if (_M0L3lenS132 >= 0) {
      if (_M0L11src__offsetS133 >= 0) {
        if (_M0L11dst__offsetS134 >= 0) {
          int32_t _M0L6_2atmpS1814 = _M0L11src__offsetS133 + _M0L3lenS132;
          int32_t _M0L6_2atmpS1815;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1815
          = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS135);
          if (_M0L6_2atmpS1814 <= _M0L6_2atmpS1815) {
            int32_t _M0L6_2atmpS1813 = _M0L11dst__offsetS134 + _M0L3lenS132;
            _if__result_4066 = _M0L6_2atmpS1813 <= _M0L13allocate__lenS131;
          } else {
            _if__result_4066 = 0;
          }
        } else {
          _if__result_4066 = 0;
        }
      } else {
        _if__result_4066 = 0;
      }
    } else {
      _if__result_4066 = 0;
    }
  } else {
    _if__result_4066 = 0;
  }
  if (_if__result_4066) {
    moonbit_incref(_M0L3srcS135);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TUsiE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS131, 0, _M0L3srcS135, _M0L11src__offsetS133, _M0L11dst__offsetS134, _M0L3lenS132);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS136;
    int32_t _M0L6_2atmpS1817;
    moonbit_string_t _M0L6_2atmpS1816;
    struct _M0TUsiE** _result_4067;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS136
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS136, (moonbit_string_t)moonbit_string_literal_67.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS136, _M0L13allocate__lenS131);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS136, (moonbit_string_t)moonbit_string_literal_68.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS136, _M0L11src__offsetS133);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS136, (moonbit_string_t)moonbit_string_literal_69.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS136, _M0L11dst__offsetS134);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS136, (moonbit_string_t)moonbit_string_literal_70.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS136, _M0L3lenS132);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS136, (moonbit_string_t)moonbit_string_literal_71.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1817 = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS135);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS136, _M0L6_2atmpS1817);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1816
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS136);
    moonbit_decref(_M0L18_2astring__builderS136);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4067
    = _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(_M0L6_2atmpS1816);
    moonbit_decref(_M0L6_2atmpS1816);
    return _result_4067;
  }
}

struct _M0TP26biolab8bio__seq8TreeNode** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3srcS141,
  int32_t _M0L13allocate__lenS137,
  int32_t _M0L3lenS138,
  int32_t _M0L11src__offsetS139,
  int32_t _M0L11dst__offsetS140
) {
  int32_t _if__result_4068;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS137 >= 0) {
    if (_M0L3lenS138 >= 0) {
      if (_M0L11src__offsetS139 >= 0) {
        if (_M0L11dst__offsetS140 >= 0) {
          int32_t _M0L6_2atmpS1819 = _M0L11src__offsetS139 + _M0L3lenS138;
          int32_t _M0L6_2atmpS1820;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1820
          = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L3srcS141);
          if (_M0L6_2atmpS1819 <= _M0L6_2atmpS1820) {
            int32_t _M0L6_2atmpS1818 = _M0L11dst__offsetS140 + _M0L3lenS138;
            _if__result_4068 = _M0L6_2atmpS1818 <= _M0L13allocate__lenS137;
          } else {
            _if__result_4068 = 0;
          }
        } else {
          _if__result_4068 = 0;
        }
      } else {
        _if__result_4068 = 0;
      }
    } else {
      _if__result_4068 = 0;
    }
  } else {
    _if__result_4068 = 0;
  }
  if (_if__result_4068) {
    moonbit_incref(_M0L3srcS141);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TP26biolab8bio__seq8TreeNode**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS137, 0, _M0L3srcS141, _M0L11src__offsetS139, _M0L11dst__offsetS140, _M0L3lenS138);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS142;
    int32_t _M0L6_2atmpS1822;
    moonbit_string_t _M0L6_2atmpS1821;
    struct _M0TP26biolab8bio__seq8TreeNode** _result_4069;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS142
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS142, (moonbit_string_t)moonbit_string_literal_67.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS142, _M0L13allocate__lenS137);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS142, (moonbit_string_t)moonbit_string_literal_68.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS142, _M0L11src__offsetS139);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS142, (moonbit_string_t)moonbit_string_literal_69.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS142, _M0L11dst__offsetS140);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS142, (moonbit_string_t)moonbit_string_literal_70.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS142, _M0L3lenS138);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS142, (moonbit_string_t)moonbit_string_literal_71.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1822
    = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq8TreeNodeE(_M0L3srcS141);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS142, _M0L6_2atmpS1822);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1821
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS142);
    moonbit_decref(_M0L18_2astring__builderS142);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4069
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L6_2atmpS1821);
    moonbit_decref(_M0L6_2atmpS1821);
    return _result_4069;
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0MPB18UninitializedArray23make__and__blit_2einnerGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3srcS147,
  int32_t _M0L13allocate__lenS143,
  int32_t _M0L3lenS144,
  int32_t _M0L11src__offsetS145,
  int32_t _M0L11dst__offsetS146
) {
  int32_t _if__result_4070;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS143 >= 0) {
    if (_M0L3lenS144 >= 0) {
      if (_M0L11src__offsetS145 >= 0) {
        if (_M0L11dst__offsetS146 >= 0) {
          int32_t _M0L6_2atmpS1824 = _M0L11src__offsetS145 + _M0L3lenS144;
          int32_t _M0L6_2atmpS1825;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1825
          = _M0MPB18UninitializedArray6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L3srcS147);
          if (_M0L6_2atmpS1824 <= _M0L6_2atmpS1825) {
            int32_t _M0L6_2atmpS1823 = _M0L11dst__offsetS146 + _M0L3lenS144;
            _if__result_4070 = _M0L6_2atmpS1823 <= _M0L13allocate__lenS143;
          } else {
            _if__result_4070 = 0;
          }
        } else {
          _if__result_4070 = 0;
        }
      } else {
        _if__result_4070 = 0;
      }
    } else {
      _if__result_4070 = 0;
    }
  } else {
    _if__result_4070 = 0;
  }
  if (_if__result_4070) {
    moonbit_incref(_M0L3srcS147);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS143, 0, _M0L3srcS147, _M0L11src__offsetS145, _M0L11dst__offsetS146, _M0L3lenS144);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS148;
    int32_t _M0L6_2atmpS1827;
    moonbit_string_t _M0L6_2atmpS1826;
    struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _result_4071;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS148
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS148, (moonbit_string_t)moonbit_string_literal_67.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS148, _M0L13allocate__lenS143);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS148, (moonbit_string_t)moonbit_string_literal_68.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS148, _M0L11src__offsetS145);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS148, (moonbit_string_t)moonbit_string_literal_69.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS148, _M0L11dst__offsetS146);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS148, (moonbit_string_t)moonbit_string_literal_70.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS148, _M0L3lenS144);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS148, (moonbit_string_t)moonbit_string_literal_71.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1827
    = _M0MPB18UninitializedArray6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L3srcS147);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS148, _M0L6_2atmpS1827);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1826
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS148);
    moonbit_decref(_M0L18_2astring__builderS148);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_4071
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEEE(_M0L6_2atmpS1826);
    moonbit_decref(_M0L6_2atmpS1826);
    return _result_4071;
  }
}

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder* _M0L4selfS114,
  moonbit_string_t _M0L3objS113
) {
  struct _M0TPB6Logger _M0L6_2atmpS1802;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS114);
  _M0L6_2atmpS1802
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS114
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS113, _M0L6_2atmpS1802);
  if (_M0L6_2atmpS1802.$1) {
    moonbit_decref(_M0L6_2atmpS1802.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder* _M0L4selfS116,
  int32_t _M0L3objS115
) {
  struct _M0TPB6Logger _M0L6_2atmpS1803;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS116);
  _M0L6_2atmpS1803
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS116
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L3objS115, _M0L6_2atmpS1803);
  if (_M0L6_2atmpS1803.$1) {
    moonbit_decref(_M0L6_2atmpS1803.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGbE(
  struct _M0TPB13StringBuilder* _M0L4selfS118,
  int32_t _M0L3objS117
) {
  struct _M0TPB6Logger _M0L6_2atmpS1804;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS118);
  _M0L6_2atmpS1804
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS118
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGbE(_M0L3objS117, _M0L6_2atmpS1804);
  if (_M0L6_2atmpS1804.$1) {
    moonbit_decref(_M0L6_2atmpS1804.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGRPB9SourceLocE(
  struct _M0TPB13StringBuilder* _M0L4selfS120,
  moonbit_string_t _M0L3objS119
) {
  struct _M0TPB6Logger _M0L6_2atmpS1805;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS120);
  _M0L6_2atmpS1805
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS120
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IPB9SourceLocPB4Show6output(_M0L3objS119, _M0L6_2atmpS1805);
  if (_M0L6_2atmpS1805.$1) {
    moonbit_decref(_M0L6_2atmpS1805.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGRPC16string10StringViewE(
  struct _M0TPB13StringBuilder* _M0L4selfS122,
  struct _M0TPC16string10StringView _M0L3objS121
) {
  struct _M0TPB6Logger _M0L6_2atmpS1806;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS122);
  _M0L6_2atmpS1806
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS122
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IPC16string10StringViewPB4Show6output(_M0L3objS121, _M0L6_2atmpS1806);
  if (_M0L6_2atmpS1806.$1) {
    moonbit_decref(_M0L6_2atmpS1806.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGmE(
  struct _M0TPB13StringBuilder* _M0L4selfS124,
  uint64_t _M0L3objS123
) {
  struct _M0TPB6Logger _M0L6_2atmpS1807;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS124);
  _M0L6_2atmpS1807
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS124
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGmE(_M0L3objS123, _M0L6_2atmpS1807);
  if (_M0L6_2atmpS1807.$1) {
    moonbit_decref(_M0L6_2atmpS1807.$1);
  }
  return 0;
}

moonbit_string_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGsE(
  moonbit_string_t* _M0L3srcS92,
  int32_t _M0L13allocate__lenS90,
  int32_t _M0L11src__offsetS93,
  int32_t _M0L11dst__offsetS91,
  int32_t _M0L9blit__lenS94
) {
  moonbit_string_t* _M0L3dstS89;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS89
  = (moonbit_string_t*)moonbit_make_ref_array(_M0L13allocate__lenS90, (moonbit_string_t)moonbit_string_literal_1.data);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGsE(_M0L3dstS89, _M0L11dst__offsetS91, _M0L3srcS92, _M0L11src__offsetS93, _M0L9blit__lenS94);
  moonbit_decref(_M0L3srcS92);
  return _M0L3dstS89;
}

struct _M0TUsiE** _M0MPB18UninitializedArray23unsafe__make__and__blitGUsiEE(
  struct _M0TUsiE** _M0L3srcS98,
  int32_t _M0L13allocate__lenS96,
  int32_t _M0L11src__offsetS99,
  int32_t _M0L11dst__offsetS97,
  int32_t _M0L9blit__lenS100
) {
  struct _M0TUsiE** _M0L3dstS95;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS95
  = (struct _M0TUsiE**)moonbit_make_ref_array(_M0L13allocate__lenS96, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGUsiEE(_M0L3dstS95, _M0L11dst__offsetS97, _M0L3srcS98, _M0L11src__offsetS99, _M0L9blit__lenS100);
  moonbit_decref(_M0L3srcS98);
  return _M0L3dstS95;
}

struct _M0TP26biolab8bio__seq8TreeNode** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3srcS104,
  int32_t _M0L13allocate__lenS102,
  int32_t _M0L11src__offsetS105,
  int32_t _M0L11dst__offsetS103,
  int32_t _M0L9blit__lenS106
) {
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3dstS101;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS101
  = (struct _M0TP26biolab8bio__seq8TreeNode**)moonbit_make_ref_array(_M0L13allocate__lenS102, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq8TreeNodeE(_M0L3dstS101, _M0L11dst__offsetS103, _M0L3srcS104, _M0L11src__offsetS105, _M0L9blit__lenS106);
  moonbit_decref(_M0L3srcS104);
  return _M0L3dstS101;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0MPB18UninitializedArray23unsafe__make__and__blitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3srcS110,
  int32_t _M0L13allocate__lenS108,
  int32_t _M0L11src__offsetS111,
  int32_t _M0L11dst__offsetS109,
  int32_t _M0L9blit__lenS112
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3dstS107;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS107
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE**)moonbit_make_ref_array(_M0L13allocate__lenS108, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(_M0L3dstS107, _M0L11dst__offsetS109, _M0L3srcS110, _M0L11src__offsetS111, _M0L9blit__lenS112);
  moonbit_decref(_M0L3srcS110);
  return _M0L3dstS107;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGsE(
  moonbit_string_t* _M0L3dstS69,
  int32_t _M0L11dst__offsetS70,
  moonbit_string_t* _M0L3srcS71,
  int32_t _M0L11src__offsetS72,
  int32_t _M0L3lenS73
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS71);
  moonbit_incref(_M0L3dstS69);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS69, _M0L11dst__offsetS70, _M0L3srcS71, _M0L11src__offsetS72, _M0L3lenS73);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGUsiEE(
  struct _M0TUsiE** _M0L3dstS74,
  int32_t _M0L11dst__offsetS75,
  struct _M0TUsiE** _M0L3srcS76,
  int32_t _M0L11src__offsetS77,
  int32_t _M0L3lenS78
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS76);
  moonbit_incref(_M0L3dstS74);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS74, _M0L11dst__offsetS75, _M0L3srcS76, _M0L11src__offsetS77, _M0L3lenS78);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3dstS79,
  int32_t _M0L11dst__offsetS80,
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3srcS81,
  int32_t _M0L11src__offsetS82,
  int32_t _M0L3lenS83
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS81);
  moonbit_incref(_M0L3dstS79);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS79, _M0L11dst__offsetS80, _M0L3srcS81, _M0L11src__offsetS82, _M0L3lenS83);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3dstS84,
  int32_t _M0L11dst__offsetS85,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3srcS86,
  int32_t _M0L11src__offsetS87,
  int32_t _M0L3lenS88
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS86);
  moonbit_incref(_M0L3dstS84);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS84, _M0L11dst__offsetS85, _M0L3srcS86, _M0L11src__offsetS87, _M0L3lenS88);
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGkE(
  uint16_t* _M0L3dstS24,
  int32_t _M0L11dst__offsetS26,
  uint16_t* _M0L3srcS25,
  int32_t _M0L11src__offsetS27,
  int32_t _M0L3lenS29
) {
  int32_t _if__result_4072;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS24 == _M0L3srcS25) {
    _if__result_4072 = _M0L11dst__offsetS26 < _M0L11src__offsetS27;
  } else {
    _if__result_4072 = 0;
  }
  if (_if__result_4072) {
    int32_t _M0L1iS28 = 0;
    while (1) {
      if (_M0L1iS28 < _M0L3lenS29) {
        int32_t _M0L6_2atmpS1757 = _M0L11dst__offsetS26 + _M0L1iS28;
        int32_t _M0L6_2atmpS1759 = _M0L11src__offsetS27 + _M0L1iS28;
        int32_t _M0L6_2atmpS1758;
        int32_t _M0L6_2atmpS1760;
        if (
          _M0L6_2atmpS1759 < 0
          || _M0L6_2atmpS1759 >= Moonbit_array_length(_M0L3srcS25)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1758 = (int32_t)_M0L3srcS25[_M0L6_2atmpS1759];
        if (
          _M0L6_2atmpS1757 < 0
          || _M0L6_2atmpS1757 >= Moonbit_array_length(_M0L3dstS24)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS24[_M0L6_2atmpS1757] = _M0L6_2atmpS1758;
        _M0L6_2atmpS1760 = _M0L1iS28 + 1;
        _M0L1iS28 = _M0L6_2atmpS1760;
        continue;
      } else {
        moonbit_decref(_M0L3srcS25);
        moonbit_decref(_M0L3dstS24);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1765 = _M0L3lenS29 - 1;
    int32_t _M0L1iS31 = _M0L6_2atmpS1765;
    while (1) {
      if (_M0L1iS31 >= 0) {
        int32_t _M0L6_2atmpS1761 = _M0L11dst__offsetS26 + _M0L1iS31;
        int32_t _M0L6_2atmpS1763 = _M0L11src__offsetS27 + _M0L1iS31;
        int32_t _M0L6_2atmpS1762;
        int32_t _M0L6_2atmpS1764;
        if (
          _M0L6_2atmpS1763 < 0
          || _M0L6_2atmpS1763 >= Moonbit_array_length(_M0L3srcS25)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1762 = (int32_t)_M0L3srcS25[_M0L6_2atmpS1763];
        if (
          _M0L6_2atmpS1761 < 0
          || _M0L6_2atmpS1761 >= Moonbit_array_length(_M0L3dstS24)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS24[_M0L6_2atmpS1761] = _M0L6_2atmpS1762;
        _M0L6_2atmpS1764 = _M0L1iS31 - 1;
        _M0L1iS31 = _M0L6_2atmpS1764;
        continue;
      } else {
        moonbit_decref(_M0L3srcS25);
        moonbit_decref(_M0L3dstS24);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGsEE(
  moonbit_string_t* _M0L3dstS33,
  int32_t _M0L11dst__offsetS35,
  moonbit_string_t* _M0L3srcS34,
  int32_t _M0L11src__offsetS36,
  int32_t _M0L3lenS38
) {
  int32_t _if__result_4075;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS33 == _M0L3srcS34) {
    _if__result_4075 = _M0L11dst__offsetS35 < _M0L11src__offsetS36;
  } else {
    _if__result_4075 = 0;
  }
  if (_if__result_4075) {
    int32_t _M0L1iS37 = 0;
    while (1) {
      if (_M0L1iS37 < _M0L3lenS38) {
        int32_t _M0L6_2atmpS1766 = _M0L11dst__offsetS35 + _M0L1iS37;
        int32_t _M0L6_2atmpS1768 = _M0L11src__offsetS36 + _M0L1iS37;
        moonbit_string_t _M0L6_2atmpS1767;
        moonbit_string_t _M0L6_2aoldS3681;
        int32_t _M0L6_2atmpS1769;
        if (
          _M0L6_2atmpS1768 < 0
          || _M0L6_2atmpS1768 >= Moonbit_array_length(_M0L3srcS34)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1767 = (moonbit_string_t)_M0L3srcS34[_M0L6_2atmpS1768];
        if (
          _M0L6_2atmpS1766 < 0
          || _M0L6_2atmpS1766 >= Moonbit_array_length(_M0L3dstS33)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3681 = (moonbit_string_t)_M0L3dstS33[_M0L6_2atmpS1766];
        moonbit_incref(_M0L6_2atmpS1767);
        moonbit_decref(_M0L6_2aoldS3681);
        _M0L3dstS33[_M0L6_2atmpS1766] = _M0L6_2atmpS1767;
        _M0L6_2atmpS1769 = _M0L1iS37 + 1;
        _M0L1iS37 = _M0L6_2atmpS1769;
        continue;
      } else {
        moonbit_decref(_M0L3srcS34);
        moonbit_decref(_M0L3dstS33);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1774 = _M0L3lenS38 - 1;
    int32_t _M0L1iS40 = _M0L6_2atmpS1774;
    while (1) {
      if (_M0L1iS40 >= 0) {
        int32_t _M0L6_2atmpS1770 = _M0L11dst__offsetS35 + _M0L1iS40;
        int32_t _M0L6_2atmpS1772 = _M0L11src__offsetS36 + _M0L1iS40;
        moonbit_string_t _M0L6_2atmpS1771;
        moonbit_string_t _M0L6_2aoldS3683;
        int32_t _M0L6_2atmpS1773;
        if (
          _M0L6_2atmpS1772 < 0
          || _M0L6_2atmpS1772 >= Moonbit_array_length(_M0L3srcS34)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1771 = (moonbit_string_t)_M0L3srcS34[_M0L6_2atmpS1772];
        if (
          _M0L6_2atmpS1770 < 0
          || _M0L6_2atmpS1770 >= Moonbit_array_length(_M0L3dstS33)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3683 = (moonbit_string_t)_M0L3dstS33[_M0L6_2atmpS1770];
        moonbit_incref(_M0L6_2atmpS1771);
        moonbit_decref(_M0L6_2aoldS3683);
        _M0L3dstS33[_M0L6_2atmpS1770] = _M0L6_2atmpS1771;
        _M0L6_2atmpS1773 = _M0L1iS40 - 1;
        _M0L1iS40 = _M0L6_2atmpS1773;
        continue;
      } else {
        moonbit_decref(_M0L3srcS34);
        moonbit_decref(_M0L3dstS33);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGUsiEEE(
  struct _M0TUsiE** _M0L3dstS42,
  int32_t _M0L11dst__offsetS44,
  struct _M0TUsiE** _M0L3srcS43,
  int32_t _M0L11src__offsetS45,
  int32_t _M0L3lenS47
) {
  int32_t _if__result_4078;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS42 == _M0L3srcS43) {
    _if__result_4078 = _M0L11dst__offsetS44 < _M0L11src__offsetS45;
  } else {
    _if__result_4078 = 0;
  }
  if (_if__result_4078) {
    int32_t _M0L1iS46 = 0;
    while (1) {
      if (_M0L1iS46 < _M0L3lenS47) {
        int32_t _M0L6_2atmpS1775 = _M0L11dst__offsetS44 + _M0L1iS46;
        int32_t _M0L6_2atmpS1777 = _M0L11src__offsetS45 + _M0L1iS46;
        struct _M0TUsiE* _M0L6_2atmpS1776;
        struct _M0TUsiE* _M0L6_2aoldS3685;
        int32_t _M0L6_2atmpS1778;
        if (
          _M0L6_2atmpS1777 < 0
          || _M0L6_2atmpS1777 >= Moonbit_array_length(_M0L3srcS43)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1776 = (struct _M0TUsiE*)_M0L3srcS43[_M0L6_2atmpS1777];
        if (
          _M0L6_2atmpS1775 < 0
          || _M0L6_2atmpS1775 >= Moonbit_array_length(_M0L3dstS42)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3685 = (struct _M0TUsiE*)_M0L3dstS42[_M0L6_2atmpS1775];
        if (_M0L6_2atmpS1776) {
          moonbit_incref(_M0L6_2atmpS1776);
        }
        if (_M0L6_2aoldS3685) {
          moonbit_decref(_M0L6_2aoldS3685);
        }
        _M0L3dstS42[_M0L6_2atmpS1775] = _M0L6_2atmpS1776;
        _M0L6_2atmpS1778 = _M0L1iS46 + 1;
        _M0L1iS46 = _M0L6_2atmpS1778;
        continue;
      } else {
        moonbit_decref(_M0L3srcS43);
        moonbit_decref(_M0L3dstS42);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1783 = _M0L3lenS47 - 1;
    int32_t _M0L1iS49 = _M0L6_2atmpS1783;
    while (1) {
      if (_M0L1iS49 >= 0) {
        int32_t _M0L6_2atmpS1779 = _M0L11dst__offsetS44 + _M0L1iS49;
        int32_t _M0L6_2atmpS1781 = _M0L11src__offsetS45 + _M0L1iS49;
        struct _M0TUsiE* _M0L6_2atmpS1780;
        struct _M0TUsiE* _M0L6_2aoldS3687;
        int32_t _M0L6_2atmpS1782;
        if (
          _M0L6_2atmpS1781 < 0
          || _M0L6_2atmpS1781 >= Moonbit_array_length(_M0L3srcS43)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1780 = (struct _M0TUsiE*)_M0L3srcS43[_M0L6_2atmpS1781];
        if (
          _M0L6_2atmpS1779 < 0
          || _M0L6_2atmpS1779 >= Moonbit_array_length(_M0L3dstS42)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3687 = (struct _M0TUsiE*)_M0L3dstS42[_M0L6_2atmpS1779];
        if (_M0L6_2atmpS1780) {
          moonbit_incref(_M0L6_2atmpS1780);
        }
        if (_M0L6_2aoldS3687) {
          moonbit_decref(_M0L6_2aoldS3687);
        }
        _M0L3dstS42[_M0L6_2atmpS1779] = _M0L6_2atmpS1780;
        _M0L6_2atmpS1782 = _M0L1iS49 - 1;
        _M0L1iS49 = _M0L6_2atmpS1782;
        continue;
      } else {
        moonbit_decref(_M0L3srcS43);
        moonbit_decref(_M0L3dstS42);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3dstS51,
  int32_t _M0L11dst__offsetS53,
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L3srcS52,
  int32_t _M0L11src__offsetS54,
  int32_t _M0L3lenS56
) {
  int32_t _if__result_4081;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS51 == _M0L3srcS52) {
    _if__result_4081 = _M0L11dst__offsetS53 < _M0L11src__offsetS54;
  } else {
    _if__result_4081 = 0;
  }
  if (_if__result_4081) {
    int32_t _M0L1iS55 = 0;
    while (1) {
      if (_M0L1iS55 < _M0L3lenS56) {
        int32_t _M0L6_2atmpS1784 = _M0L11dst__offsetS53 + _M0L1iS55;
        int32_t _M0L6_2atmpS1786 = _M0L11src__offsetS54 + _M0L1iS55;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS1785;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2aoldS3689;
        int32_t _M0L6_2atmpS1787;
        if (
          _M0L6_2atmpS1786 < 0
          || _M0L6_2atmpS1786 >= Moonbit_array_length(_M0L3srcS52)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1785
        = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3srcS52[
            _M0L6_2atmpS1786
          ];
        if (
          _M0L6_2atmpS1784 < 0
          || _M0L6_2atmpS1784 >= Moonbit_array_length(_M0L3dstS51)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3689
        = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3dstS51[
            _M0L6_2atmpS1784
          ];
        if (_M0L6_2atmpS1785) {
          moonbit_incref(_M0L6_2atmpS1785);
        }
        if (_M0L6_2aoldS3689) {
          moonbit_decref(_M0L6_2aoldS3689);
        }
        _M0L3dstS51[_M0L6_2atmpS1784] = _M0L6_2atmpS1785;
        _M0L6_2atmpS1787 = _M0L1iS55 + 1;
        _M0L1iS55 = _M0L6_2atmpS1787;
        continue;
      } else {
        moonbit_decref(_M0L3srcS52);
        moonbit_decref(_M0L3dstS51);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1792 = _M0L3lenS56 - 1;
    int32_t _M0L1iS58 = _M0L6_2atmpS1792;
    while (1) {
      if (_M0L1iS58 >= 0) {
        int32_t _M0L6_2atmpS1788 = _M0L11dst__offsetS53 + _M0L1iS58;
        int32_t _M0L6_2atmpS1790 = _M0L11src__offsetS54 + _M0L1iS58;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2atmpS1789;
        struct _M0TP26biolab8bio__seq8TreeNode* _M0L6_2aoldS3691;
        int32_t _M0L6_2atmpS1791;
        if (
          _M0L6_2atmpS1790 < 0
          || _M0L6_2atmpS1790 >= Moonbit_array_length(_M0L3srcS52)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1789
        = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3srcS52[
            _M0L6_2atmpS1790
          ];
        if (
          _M0L6_2atmpS1788 < 0
          || _M0L6_2atmpS1788 >= Moonbit_array_length(_M0L3dstS51)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3691
        = (struct _M0TP26biolab8bio__seq8TreeNode*)_M0L3dstS51[
            _M0L6_2atmpS1788
          ];
        if (_M0L6_2atmpS1789) {
          moonbit_incref(_M0L6_2atmpS1789);
        }
        if (_M0L6_2aoldS3691) {
          moonbit_decref(_M0L6_2aoldS3691);
        }
        _M0L3dstS51[_M0L6_2atmpS1788] = _M0L6_2atmpS1789;
        _M0L6_2atmpS1791 = _M0L1iS58 - 1;
        _M0L1iS58 = _M0L6_2atmpS1791;
        continue;
      } else {
        moonbit_decref(_M0L3srcS52);
        moonbit_decref(_M0L3dstS51);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3dstS60,
  int32_t _M0L11dst__offsetS62,
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L3srcS61,
  int32_t _M0L11src__offsetS63,
  int32_t _M0L3lenS65
) {
  int32_t _if__result_4084;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS60 == _M0L3srcS61) {
    _if__result_4084 = _M0L11dst__offsetS62 < _M0L11src__offsetS63;
  } else {
    _if__result_4084 = 0;
  }
  if (_if__result_4084) {
    int32_t _M0L1iS64 = 0;
    while (1) {
      if (_M0L1iS64 < _M0L3lenS65) {
        int32_t _M0L6_2atmpS1793 = _M0L11dst__offsetS62 + _M0L1iS64;
        int32_t _M0L6_2atmpS1795 = _M0L11src__offsetS63 + _M0L1iS64;
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS1794;
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2aoldS3693;
        int32_t _M0L6_2atmpS1796;
        if (
          _M0L6_2atmpS1795 < 0
          || _M0L6_2atmpS1795 >= Moonbit_array_length(_M0L3srcS61)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1794
        = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)_M0L3srcS61[
            _M0L6_2atmpS1795
          ];
        if (
          _M0L6_2atmpS1793 < 0
          || _M0L6_2atmpS1793 >= Moonbit_array_length(_M0L3dstS60)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3693
        = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)_M0L3dstS60[
            _M0L6_2atmpS1793
          ];
        if (_M0L6_2atmpS1794) {
          moonbit_incref(_M0L6_2atmpS1794);
        }
        if (_M0L6_2aoldS3693) {
          moonbit_decref(_M0L6_2aoldS3693);
        }
        _M0L3dstS60[_M0L6_2atmpS1793] = _M0L6_2atmpS1794;
        _M0L6_2atmpS1796 = _M0L1iS64 + 1;
        _M0L1iS64 = _M0L6_2atmpS1796;
        continue;
      } else {
        moonbit_decref(_M0L3srcS61);
        moonbit_decref(_M0L3dstS60);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1801 = _M0L3lenS65 - 1;
    int32_t _M0L1iS67 = _M0L6_2atmpS1801;
    while (1) {
      if (_M0L1iS67 >= 0) {
        int32_t _M0L6_2atmpS1797 = _M0L11dst__offsetS62 + _M0L1iS67;
        int32_t _M0L6_2atmpS1799 = _M0L11src__offsetS63 + _M0L1iS67;
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2atmpS1798;
        struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE* _M0L6_2aoldS3695;
        int32_t _M0L6_2atmpS1800;
        if (
          _M0L6_2atmpS1799 < 0
          || _M0L6_2atmpS1799 >= Moonbit_array_length(_M0L3srcS61)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1798
        = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)_M0L3srcS61[
            _M0L6_2atmpS1799
          ];
        if (
          _M0L6_2atmpS1797 < 0
          || _M0L6_2atmpS1797 >= Moonbit_array_length(_M0L3dstS60)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3695
        = (struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE*)_M0L3dstS60[
            _M0L6_2atmpS1797
          ];
        if (_M0L6_2atmpS1798) {
          moonbit_incref(_M0L6_2atmpS1798);
        }
        if (_M0L6_2aoldS3695) {
          moonbit_decref(_M0L6_2aoldS3695);
        }
        _M0L3dstS60[_M0L6_2atmpS1797] = _M0L6_2atmpS1798;
        _M0L6_2atmpS1800 = _M0L1iS67 - 1;
        _M0L1iS67 = _M0L6_2atmpS1800;
        continue;
      } else {
        moonbit_decref(_M0L3srcS61);
        moonbit_decref(_M0L3dstS60);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t* _M0L4selfS20) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS20);
}

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(
  struct _M0TUsiE** _M0L4selfS21
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS21);
}

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq8TreeNodeE(
  struct _M0TP26biolab8bio__seq8TreeNode** _M0L4selfS22
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS22);
}

int32_t _M0MPB18UninitializedArray6lengthGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0L4selfS23
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS23);
}

uint32_t _M0FPB13consume4__acc(uint32_t _M0L3accS18, uint32_t _M0L5inputS19) {
  uint32_t _M0L6_2atmpS1756;
  uint32_t _M0L6_2atmpS1755;
  uint32_t _M0L6_2atmpS1754;
  #line 465 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1756 = _M0L5inputS19 * 3266489917u;
  _M0L6_2atmpS1755 = _M0L3accS18 + _M0L6_2atmpS1756;
  #line 466 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1754 = _M0FPB4rotl(_M0L6_2atmpS1755, 17);
  return _M0L6_2atmpS1754 * 668265263u;
}

uint32_t _M0FPB4rotl(uint32_t _M0L1xS16, int32_t _M0L1rS17) {
  uint32_t _M0L6_2atmpS1751;
  int32_t _M0L6_2atmpS1753;
  uint32_t _M0L6_2atmpS1752;
  #line 475 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1751 = _M0L1xS16 << (_M0L1rS17 & 31);
  _M0L6_2atmpS1753 = 32 - _M0L1rS17;
  _M0L6_2atmpS1752 = _M0L1xS16 >> (_M0L6_2atmpS1753 & 31);
  return _M0L6_2atmpS1751 | _M0L6_2atmpS1752;
}

int32_t _M0IPB7FailurePB4Show6output(
  void* _M0L10_2ax__5721S12,
  struct _M0TPB6Logger _M0L10_2ax__5722S15
) {
  struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS13;
  moonbit_string_t _M0L15_2a_2aarg__5723S14;
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2aFailureS13
  = (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L10_2ax__5721S12;
  _M0L15_2a_2aarg__5723S14 = _M0L10_2aFailureS13->$0;
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S15.$0->$method_0(_M0L10_2ax__5722S15.$1, (moonbit_string_t)moonbit_string_literal_72.data);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB6Logger13write__objectGsE(_M0L10_2ax__5722S15, _M0L15_2a_2aarg__5723S14);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S15.$0->$method_0(_M0L10_2ax__5722S15.$1, (moonbit_string_t)moonbit_string_literal_73.data);
  return 0;
}

int32_t _M0MPB6Logger13write__objectGsE(
  struct _M0TPB6Logger _M0L4selfS11,
  moonbit_string_t _M0L3objS10
) {
  #line 173 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 174 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS10, _M0L4selfS11);
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

struct _M0TP26biolab8bio__seq8TreeNode* _M0FPC15abort5abortGRP26biolab8bio__seq8TreeNodeE(
  moonbit_string_t _M0L3msgS5
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS5);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TP26biolab8bio__seq8TreeNode** _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq8TreeNodeEE(
  moonbit_string_t _M0L3msgS6
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS6);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

struct _M0TPB5ArrayGRP26biolab8bio__seq8TreeNodeE** _M0FPC15abort5abortGRPB18UninitializedArrayGRPB5ArrayGRP26biolab8bio__seq8TreeNodeEEE(
  moonbit_string_t _M0L3msgS7
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS7);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

int32_t _M0FPC15abort5abortGiE(moonbit_string_t _M0L3msgS8) {
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

moonbit_string_t _M0FP15Error10to__string(void* _M0L4_2aeS1573) {
  switch (Moonbit_object_tag(_M0L4_2aeS1573)) {
    case 1: {
      return (moonbit_string_t)moonbit_string_literal_74.data;
      break;
    }
    
    case 6: {
      return (moonbit_string_t)moonbit_string_literal_75.data;
      break;
    }
    
    case 8: {
      return (moonbit_string_t)moonbit_string_literal_76.data;
      break;
    }
    
    case 5: {
      return (moonbit_string_t)moonbit_string_literal_77.data;
      break;
    }
    
    case 3: {
      return (moonbit_string_t)moonbit_string_literal_78.data;
      break;
    }
    
    case 0: {
      return _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(_M0L4_2aeS1573);
      break;
    }
    
    case 4: {
      return (moonbit_string_t)moonbit_string_literal_79.data;
      break;
    }
    
    case 2: {
      return (moonbit_string_t)moonbit_string_literal_80.data;
      break;
    }
    default: {
      return (moonbit_string_t)moonbit_string_literal_81.data;
      break;
    }
  }
}

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1597,
  struct _M0TPB4Show _M0L8_2aparamS1596
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1595 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1597;
  _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(_M0L7_2aselfS1595, _M0L8_2aparamS1596);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1594,
  struct _M0TPB4Show _M0L8_2aparamS1593
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1592 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1594;
  _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(_M0L7_2aselfS1592, _M0L8_2aparamS1593);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1591,
  int32_t _M0L8_2aparamS1590
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1589 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1591;
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS1589, _M0L8_2aparamS1590);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1588,
  struct _M0TPC16string10StringView _M0L8_2aparamS1587
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1586 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1588;
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L7_2aselfS1586, _M0L8_2aparamS1587);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1585,
  moonbit_string_t _M0L8_2aparamS1582,
  int32_t _M0L8_2aparamS1583,
  int32_t _M0L8_2aparamS1584
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1581 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1585;
  _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L7_2aselfS1581, _M0L8_2aparamS1582, _M0L8_2aparamS1583, _M0L8_2aparamS1584);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1580,
  moonbit_string_t _M0L8_2aparamS1579
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1578 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1580;
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L7_2aselfS1578, _M0L8_2aparamS1579);
  return 0;
}

void moonbit_init() {
  moonbit_layout_table = moonbit_layout_table_data;
  moonbit_string_t* _M0L6_2atmpS1750 =
    (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1749;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1697;
  moonbit_string_t* _M0L6_2atmpS1748;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1747;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1698;
  moonbit_string_t* _M0L6_2atmpS1746;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1745;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1699;
  moonbit_string_t* _M0L6_2atmpS1744;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1743;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1700;
  moonbit_string_t* _M0L6_2atmpS1742;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1741;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1701;
  moonbit_string_t* _M0L6_2atmpS1740;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1739;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1702;
  moonbit_string_t* _M0L6_2atmpS1738;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1737;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1703;
  moonbit_string_t* _M0L6_2atmpS1736;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1735;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1704;
  moonbit_string_t* _M0L6_2atmpS1734;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1733;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1705;
  moonbit_string_t* _M0L6_2atmpS1732;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1731;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1706;
  moonbit_string_t* _M0L6_2atmpS1730;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1729;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1707;
  moonbit_string_t* _M0L6_2atmpS1728;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1727;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1708;
  moonbit_string_t* _M0L6_2atmpS1726;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1725;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1709;
  moonbit_string_t* _M0L6_2atmpS1724;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1723;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1710;
  moonbit_string_t* _M0L6_2atmpS1722;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1721;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1711;
  moonbit_string_t* _M0L6_2atmpS1720;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1719;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1712;
  moonbit_string_t* _M0L6_2atmpS1718;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1717;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1713;
  moonbit_string_t* _M0L6_2atmpS1716;
  struct _M0TUWEuQRPC15error5ErrorNsE* _M0L8_2atupleS1715;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE* _M0L8_2atupleS1714;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1474;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1696;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1695;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1694;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1605;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1475;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1693;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1692;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1691;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1606;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1476;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1690;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1689;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1688;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1607;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1477;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1687;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1686;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1685;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1608;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1478;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1684;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1683;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1682;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1609;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1479;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1681;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1680;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1679;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1610;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1480;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1678;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1677;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1676;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1611;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1481;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1675;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1674;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1673;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1612;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1482;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1672;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1671;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1670;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1613;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1483;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1669;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1668;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1667;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1614;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1484;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1666;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1665;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1664;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1615;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1485;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1663;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1662;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1661;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1616;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1486;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1660;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1659;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1658;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1617;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1487;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1657;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1656;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1655;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1618;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1488;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1654;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1653;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1652;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1619;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1489;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1651;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1650;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1649;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1620;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1490;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1648;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1647;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1646;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1621;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1491;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1645;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1644;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1643;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1622;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1492;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1642;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1641;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1640;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1623;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1493;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1639;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1638;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1637;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1624;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1494;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1636;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1635;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1634;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1625;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1495;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1633;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1632;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1631;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1626;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L7_2abindS1496;
  struct _M0TUiUWEuQRPC15error5ErrorNsEE** _M0L6_2atmpS1630;
  struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE _M0L6_2atmpS1629;
  struct _M0TPB3MapGiUWEuQRPC15error5ErrorNsEE* _M0L6_2atmpS1628;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE* _M0L8_2atupleS1627;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L7_2abindS1473;
  struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE** _M0L6_2atmpS1604;
  struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE _M0L6_2atmpS1603;
  _M0L6_2atmpS1750[0] = (moonbit_string_t)moonbit_string_literal_82.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__0_2eclo);
  _M0L8_2atupleS1749
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1749)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1749->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__0_2eclo;
  _M0L8_2atupleS1749->$1 = _M0L6_2atmpS1750;
  _M0L8_2atupleS1697
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1697)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1697->$0 = 0;
  _M0L8_2atupleS1697->$1 = _M0L8_2atupleS1749;
  _M0L6_2atmpS1748 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1748[0] = (moonbit_string_t)moonbit_string_literal_83.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__1_2eclo);
  _M0L8_2atupleS1747
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1747)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1747->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__1_2eclo;
  _M0L8_2atupleS1747->$1 = _M0L6_2atmpS1748;
  _M0L8_2atupleS1698
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1698)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1698->$0 = 1;
  _M0L8_2atupleS1698->$1 = _M0L8_2atupleS1747;
  _M0L6_2atmpS1746 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1746[0] = (moonbit_string_t)moonbit_string_literal_84.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__2_2eclo);
  _M0L8_2atupleS1745
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1745)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1745->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__2_2eclo;
  _M0L8_2atupleS1745->$1 = _M0L6_2atmpS1746;
  _M0L8_2atupleS1699
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1699)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1699->$0 = 2;
  _M0L8_2atupleS1699->$1 = _M0L8_2atupleS1745;
  _M0L6_2atmpS1744 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1744[0] = (moonbit_string_t)moonbit_string_literal_85.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__3_2eclo);
  _M0L8_2atupleS1743
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1743)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1743->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__3_2eclo;
  _M0L8_2atupleS1743->$1 = _M0L6_2atmpS1744;
  _M0L8_2atupleS1700
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1700)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1700->$0 = 3;
  _M0L8_2atupleS1700->$1 = _M0L8_2atupleS1743;
  _M0L6_2atmpS1742 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1742[0] = (moonbit_string_t)moonbit_string_literal_86.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__4_2eclo);
  _M0L8_2atupleS1741
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1741)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1741->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__4_2eclo;
  _M0L8_2atupleS1741->$1 = _M0L6_2atmpS1742;
  _M0L8_2atupleS1701
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1701)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1701->$0 = 4;
  _M0L8_2atupleS1701->$1 = _M0L8_2atupleS1741;
  _M0L6_2atmpS1740 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1740[0] = (moonbit_string_t)moonbit_string_literal_87.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__5_2eclo);
  _M0L8_2atupleS1739
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1739)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1739->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__5_2eclo;
  _M0L8_2atupleS1739->$1 = _M0L6_2atmpS1740;
  _M0L8_2atupleS1702
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1702)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1702->$0 = 5;
  _M0L8_2atupleS1702->$1 = _M0L8_2atupleS1739;
  _M0L6_2atmpS1738 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1738[0] = (moonbit_string_t)moonbit_string_literal_88.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__6_2eclo);
  _M0L8_2atupleS1737
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1737)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1737->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__6_2eclo;
  _M0L8_2atupleS1737->$1 = _M0L6_2atmpS1738;
  _M0L8_2atupleS1703
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1703)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1703->$0 = 6;
  _M0L8_2atupleS1703->$1 = _M0L8_2atupleS1737;
  _M0L6_2atmpS1736 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1736[0] = (moonbit_string_t)moonbit_string_literal_89.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__7_2eclo);
  _M0L8_2atupleS1735
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1735)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1735->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__7_2eclo;
  _M0L8_2atupleS1735->$1 = _M0L6_2atmpS1736;
  _M0L8_2atupleS1704
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1704)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1704->$0 = 7;
  _M0L8_2atupleS1704->$1 = _M0L8_2atupleS1735;
  _M0L6_2atmpS1734 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1734[0] = (moonbit_string_t)moonbit_string_literal_90.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__8_2eclo);
  _M0L8_2atupleS1733
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1733)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1733->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__8_2eclo;
  _M0L8_2atupleS1733->$1 = _M0L6_2atmpS1734;
  _M0L8_2atupleS1705
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1705)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1705->$0 = 8;
  _M0L8_2atupleS1705->$1 = _M0L8_2atupleS1733;
  _M0L6_2atmpS1732 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1732[0] = (moonbit_string_t)moonbit_string_literal_91.data;
  moonbit_incref(_M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__9_2eclo);
  _M0L8_2atupleS1731
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1731)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1731->$0
  = _M0FP26biolab8bio__seq59____test__736b62696f5f747265655f62656e63682e6d6274__9_2eclo;
  _M0L8_2atupleS1731->$1 = _M0L6_2atmpS1732;
  _M0L8_2atupleS1706
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1706)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1706->$0 = 9;
  _M0L8_2atupleS1706->$1 = _M0L8_2atupleS1731;
  _M0L6_2atmpS1730 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1730[0] = (moonbit_string_t)moonbit_string_literal_92.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__10_2eclo);
  _M0L8_2atupleS1729
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1729)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1729->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__10_2eclo;
  _M0L8_2atupleS1729->$1 = _M0L6_2atmpS1730;
  _M0L8_2atupleS1707
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1707)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1707->$0 = 10;
  _M0L8_2atupleS1707->$1 = _M0L8_2atupleS1729;
  _M0L6_2atmpS1728 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1728[0] = (moonbit_string_t)moonbit_string_literal_93.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__11_2eclo);
  _M0L8_2atupleS1727
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1727)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1727->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__11_2eclo;
  _M0L8_2atupleS1727->$1 = _M0L6_2atmpS1728;
  _M0L8_2atupleS1708
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1708)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1708->$0 = 11;
  _M0L8_2atupleS1708->$1 = _M0L8_2atupleS1727;
  _M0L6_2atmpS1726 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1726[0] = (moonbit_string_t)moonbit_string_literal_94.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__12_2eclo);
  _M0L8_2atupleS1725
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1725)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1725->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__12_2eclo;
  _M0L8_2atupleS1725->$1 = _M0L6_2atmpS1726;
  _M0L8_2atupleS1709
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1709)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1709->$0 = 12;
  _M0L8_2atupleS1709->$1 = _M0L8_2atupleS1725;
  _M0L6_2atmpS1724 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1724[0] = (moonbit_string_t)moonbit_string_literal_95.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__13_2eclo);
  _M0L8_2atupleS1723
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1723)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1723->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__13_2eclo;
  _M0L8_2atupleS1723->$1 = _M0L6_2atmpS1724;
  _M0L8_2atupleS1710
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1710)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1710->$0 = 13;
  _M0L8_2atupleS1710->$1 = _M0L8_2atupleS1723;
  _M0L6_2atmpS1722 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1722[0] = (moonbit_string_t)moonbit_string_literal_96.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__14_2eclo);
  _M0L8_2atupleS1721
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1721)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1721->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__14_2eclo;
  _M0L8_2atupleS1721->$1 = _M0L6_2atmpS1722;
  _M0L8_2atupleS1711
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1711)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1711->$0 = 14;
  _M0L8_2atupleS1711->$1 = _M0L8_2atupleS1721;
  _M0L6_2atmpS1720 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1720[0] = (moonbit_string_t)moonbit_string_literal_97.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__15_2eclo);
  _M0L8_2atupleS1719
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1719)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1719->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__15_2eclo;
  _M0L8_2atupleS1719->$1 = _M0L6_2atmpS1720;
  _M0L8_2atupleS1712
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1712)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1712->$0 = 15;
  _M0L8_2atupleS1712->$1 = _M0L8_2atupleS1719;
  _M0L6_2atmpS1718 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1718[0] = (moonbit_string_t)moonbit_string_literal_98.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__16_2eclo);
  _M0L8_2atupleS1717
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1717)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1717->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__16_2eclo;
  _M0L8_2atupleS1717->$1 = _M0L6_2atmpS1718;
  _M0L8_2atupleS1713
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1713)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1713->$0 = 16;
  _M0L8_2atupleS1713->$1 = _M0L8_2atupleS1717;
  _M0L6_2atmpS1716 = (moonbit_string_t*)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS1716[0] = (moonbit_string_t)moonbit_string_literal_99.data;
  moonbit_incref(_M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__17_2eclo);
  _M0L8_2atupleS1715
  = (struct _M0TUWEuQRPC15error5ErrorNsE*)moonbit_malloc(sizeof(struct _M0TUWEuQRPC15error5ErrorNsE));
  Moonbit_object_header(_M0L8_2atupleS1715)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 105, 0);
  _M0L8_2atupleS1715->$0
  = _M0FP26biolab8bio__seq60____test__736b62696f5f747265655f62656e63682e6d6274__17_2eclo;
  _M0L8_2atupleS1715->$1 = _M0L6_2atmpS1716;
  _M0L8_2atupleS1714
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE*)moonbit_malloc(sizeof(struct _M0TUiUWEuQRPC15error5ErrorNsEE));
  Moonbit_object_header(_M0L8_2atupleS1714)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 109, 0);
  _M0L8_2atupleS1714->$0 = 17;
  _M0L8_2atupleS1714->$1 = _M0L8_2atupleS1715;
  _M0L7_2abindS1474
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_make_ref_array_raw(18);
  _M0L7_2abindS1474[0] = _M0L8_2atupleS1697;
  _M0L7_2abindS1474[1] = _M0L8_2atupleS1698;
  _M0L7_2abindS1474[2] = _M0L8_2atupleS1699;
  _M0L7_2abindS1474[3] = _M0L8_2atupleS1700;
  _M0L7_2abindS1474[4] = _M0L8_2atupleS1701;
  _M0L7_2abindS1474[5] = _M0L8_2atupleS1702;
  _M0L7_2abindS1474[6] = _M0L8_2atupleS1703;
  _M0L7_2abindS1474[7] = _M0L8_2atupleS1704;
  _M0L7_2abindS1474[8] = _M0L8_2atupleS1705;
  _M0L7_2abindS1474[9] = _M0L8_2atupleS1706;
  _M0L7_2abindS1474[10] = _M0L8_2atupleS1707;
  _M0L7_2abindS1474[11] = _M0L8_2atupleS1708;
  _M0L7_2abindS1474[12] = _M0L8_2atupleS1709;
  _M0L7_2abindS1474[13] = _M0L8_2atupleS1710;
  _M0L7_2abindS1474[14] = _M0L8_2atupleS1711;
  _M0L7_2abindS1474[15] = _M0L8_2atupleS1712;
  _M0L7_2abindS1474[16] = _M0L8_2atupleS1713;
  _M0L7_2abindS1474[17] = _M0L8_2atupleS1714;
  _M0L6_2atmpS1696 = _M0L7_2abindS1474;
  _M0L6_2atmpS1695
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1696, .$1 = 0, .$2 = 18
  };
  #line 400 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1694
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1695, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1695.$0);
  _M0L8_2atupleS1605
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1605)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1605->$0 = (moonbit_string_t)moonbit_string_literal_100.data;
  _M0L8_2atupleS1605->$1 = _M0L6_2atmpS1694;
  _M0L7_2abindS1475
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1693 = _M0L7_2abindS1475;
  _M0L6_2atmpS1692
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1693, .$1 = 0, .$2 = 0
  };
  #line 420 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1691
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1692, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1692.$0);
  _M0L8_2atupleS1606
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1606)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1606->$0 = (moonbit_string_t)moonbit_string_literal_101.data;
  _M0L8_2atupleS1606->$1 = _M0L6_2atmpS1691;
  _M0L7_2abindS1476
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1690 = _M0L7_2abindS1476;
  _M0L6_2atmpS1689
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1690, .$1 = 0, .$2 = 0
  };
  #line 422 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1688
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1689, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1689.$0);
  _M0L8_2atupleS1607
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1607)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1607->$0 = (moonbit_string_t)moonbit_string_literal_102.data;
  _M0L8_2atupleS1607->$1 = _M0L6_2atmpS1688;
  _M0L7_2abindS1477
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1687 = _M0L7_2abindS1477;
  _M0L6_2atmpS1686
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1687, .$1 = 0, .$2 = 0
  };
  #line 424 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1685
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1686, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1686.$0);
  _M0L8_2atupleS1608
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1608)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1608->$0 = (moonbit_string_t)moonbit_string_literal_103.data;
  _M0L8_2atupleS1608->$1 = _M0L6_2atmpS1685;
  _M0L7_2abindS1478
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1684 = _M0L7_2abindS1478;
  _M0L6_2atmpS1683
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1684, .$1 = 0, .$2 = 0
  };
  #line 426 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1682
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1683, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1683.$0);
  _M0L8_2atupleS1609
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1609)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1609->$0 = (moonbit_string_t)moonbit_string_literal_104.data;
  _M0L8_2atupleS1609->$1 = _M0L6_2atmpS1682;
  _M0L7_2abindS1479
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1681 = _M0L7_2abindS1479;
  _M0L6_2atmpS1680
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1681, .$1 = 0, .$2 = 0
  };
  #line 428 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1679
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1680, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1680.$0);
  _M0L8_2atupleS1610
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1610)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1610->$0 = (moonbit_string_t)moonbit_string_literal_105.data;
  _M0L8_2atupleS1610->$1 = _M0L6_2atmpS1679;
  _M0L7_2abindS1480
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1678 = _M0L7_2abindS1480;
  _M0L6_2atmpS1677
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1678, .$1 = 0, .$2 = 0
  };
  #line 430 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1676
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1677, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1677.$0);
  _M0L8_2atupleS1611
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1611)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1611->$0 = (moonbit_string_t)moonbit_string_literal_106.data;
  _M0L8_2atupleS1611->$1 = _M0L6_2atmpS1676;
  _M0L7_2abindS1481
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1675 = _M0L7_2abindS1481;
  _M0L6_2atmpS1674
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1675, .$1 = 0, .$2 = 0
  };
  #line 432 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1673
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1674, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1674.$0);
  _M0L8_2atupleS1612
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1612)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1612->$0 = (moonbit_string_t)moonbit_string_literal_107.data;
  _M0L8_2atupleS1612->$1 = _M0L6_2atmpS1673;
  _M0L7_2abindS1482
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1672 = _M0L7_2abindS1482;
  _M0L6_2atmpS1671
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1672, .$1 = 0, .$2 = 0
  };
  #line 434 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1670
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1671, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1671.$0);
  _M0L8_2atupleS1613
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1613)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1613->$0 = (moonbit_string_t)moonbit_string_literal_108.data;
  _M0L8_2atupleS1613->$1 = _M0L6_2atmpS1670;
  _M0L7_2abindS1483
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1669 = _M0L7_2abindS1483;
  _M0L6_2atmpS1668
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1669, .$1 = 0, .$2 = 0
  };
  #line 436 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1667
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1668, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1668.$0);
  _M0L8_2atupleS1614
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1614)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1614->$0 = (moonbit_string_t)moonbit_string_literal_109.data;
  _M0L8_2atupleS1614->$1 = _M0L6_2atmpS1667;
  _M0L7_2abindS1484
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1666 = _M0L7_2abindS1484;
  _M0L6_2atmpS1665
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1666, .$1 = 0, .$2 = 0
  };
  #line 438 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1664
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1665, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1665.$0);
  _M0L8_2atupleS1615
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1615)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1615->$0 = (moonbit_string_t)moonbit_string_literal_110.data;
  _M0L8_2atupleS1615->$1 = _M0L6_2atmpS1664;
  _M0L7_2abindS1485
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1663 = _M0L7_2abindS1485;
  _M0L6_2atmpS1662
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1663, .$1 = 0, .$2 = 0
  };
  #line 440 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1661
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1662, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1662.$0);
  _M0L8_2atupleS1616
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1616)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1616->$0 = (moonbit_string_t)moonbit_string_literal_111.data;
  _M0L8_2atupleS1616->$1 = _M0L6_2atmpS1661;
  _M0L7_2abindS1486
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1660 = _M0L7_2abindS1486;
  _M0L6_2atmpS1659
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1660, .$1 = 0, .$2 = 0
  };
  #line 442 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1658
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1659, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1659.$0);
  _M0L8_2atupleS1617
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1617)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1617->$0 = (moonbit_string_t)moonbit_string_literal_112.data;
  _M0L8_2atupleS1617->$1 = _M0L6_2atmpS1658;
  _M0L7_2abindS1487
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1657 = _M0L7_2abindS1487;
  _M0L6_2atmpS1656
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1657, .$1 = 0, .$2 = 0
  };
  #line 444 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1655
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1656, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1656.$0);
  _M0L8_2atupleS1618
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1618)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1618->$0 = (moonbit_string_t)moonbit_string_literal_113.data;
  _M0L8_2atupleS1618->$1 = _M0L6_2atmpS1655;
  _M0L7_2abindS1488
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1654 = _M0L7_2abindS1488;
  _M0L6_2atmpS1653
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1654, .$1 = 0, .$2 = 0
  };
  #line 446 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1652
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1653, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1653.$0);
  _M0L8_2atupleS1619
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1619)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1619->$0 = (moonbit_string_t)moonbit_string_literal_114.data;
  _M0L8_2atupleS1619->$1 = _M0L6_2atmpS1652;
  _M0L7_2abindS1489
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1651 = _M0L7_2abindS1489;
  _M0L6_2atmpS1650
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1651, .$1 = 0, .$2 = 0
  };
  #line 448 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1649
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1650, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1650.$0);
  _M0L8_2atupleS1620
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1620)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1620->$0 = (moonbit_string_t)moonbit_string_literal_115.data;
  _M0L8_2atupleS1620->$1 = _M0L6_2atmpS1649;
  _M0L7_2abindS1490
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1648 = _M0L7_2abindS1490;
  _M0L6_2atmpS1647
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1648, .$1 = 0, .$2 = 0
  };
  #line 450 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1646
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1647, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1647.$0);
  _M0L8_2atupleS1621
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1621)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1621->$0 = (moonbit_string_t)moonbit_string_literal_116.data;
  _M0L8_2atupleS1621->$1 = _M0L6_2atmpS1646;
  _M0L7_2abindS1491
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1645 = _M0L7_2abindS1491;
  _M0L6_2atmpS1644
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1645, .$1 = 0, .$2 = 0
  };
  #line 452 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1643
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1644, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1644.$0);
  _M0L8_2atupleS1622
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1622)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1622->$0 = (moonbit_string_t)moonbit_string_literal_117.data;
  _M0L8_2atupleS1622->$1 = _M0L6_2atmpS1643;
  _M0L7_2abindS1492
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1642 = _M0L7_2abindS1492;
  _M0L6_2atmpS1641
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1642, .$1 = 0, .$2 = 0
  };
  #line 454 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1640
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1641, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1641.$0);
  _M0L8_2atupleS1623
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1623)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1623->$0 = (moonbit_string_t)moonbit_string_literal_118.data;
  _M0L8_2atupleS1623->$1 = _M0L6_2atmpS1640;
  _M0L7_2abindS1493
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1639 = _M0L7_2abindS1493;
  _M0L6_2atmpS1638
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1639, .$1 = 0, .$2 = 0
  };
  #line 456 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1637
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1638, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1638.$0);
  _M0L8_2atupleS1624
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1624)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1624->$0 = (moonbit_string_t)moonbit_string_literal_119.data;
  _M0L8_2atupleS1624->$1 = _M0L6_2atmpS1637;
  _M0L7_2abindS1494
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1636 = _M0L7_2abindS1494;
  _M0L6_2atmpS1635
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1636, .$1 = 0, .$2 = 0
  };
  #line 458 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1634
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1635, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1635.$0);
  _M0L8_2atupleS1625
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1625)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1625->$0 = (moonbit_string_t)moonbit_string_literal_120.data;
  _M0L8_2atupleS1625->$1 = _M0L6_2atmpS1634;
  _M0L7_2abindS1495
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1633 = _M0L7_2abindS1495;
  _M0L6_2atmpS1632
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1633, .$1 = 0, .$2 = 0
  };
  #line 460 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1631
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1632, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1632.$0);
  _M0L8_2atupleS1626
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1626)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1626->$0 = (moonbit_string_t)moonbit_string_literal_121.data;
  _M0L8_2atupleS1626->$1 = _M0L6_2atmpS1631;
  _M0L7_2abindS1496
  = (struct _M0TUiUWEuQRPC15error5ErrorNsEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS1630 = _M0L7_2abindS1496;
  _M0L6_2atmpS1629
  = (struct _M0TPB9ArrayViewGUiUWEuQRPC15error5ErrorNsEEE){
    .$0 = _M0L6_2atmpS1630, .$1 = 0, .$2 = 0
  };
  #line 462 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS1628
  = _M0MPB3Map3MapGiUWEuQRPC15error5ErrorNsEE(_M0L6_2atmpS1629, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1629.$0);
  _M0L8_2atupleS1627
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE*)moonbit_malloc(sizeof(struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE));
  Moonbit_object_header(_M0L8_2atupleS1627)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 112, 0);
  _M0L8_2atupleS1627->$0 = (moonbit_string_t)moonbit_string_literal_122.data;
  _M0L8_2atupleS1627->$1 = _M0L6_2atmpS1628;
  _M0L7_2abindS1473
  = (struct _M0TUsRPB3MapGiUWEuQRPC15error5ErrorNsEEE**)moonbit_make_ref_array_raw(23);
  _M0L7_2abindS1473[0] = _M0L8_2atupleS1605;
  _M0L7_2abindS1473[1] = _M0L8_2atupleS1606;
  _M0L7_2abindS1473[2] = _M0L8_2atupleS1607;
  _M0L7_2abindS1473[3] = _M0L8_2atupleS1608;
  _M0L7_2abindS1473[4] = _M0L8_2atupleS1609;
  _M0L7_2abindS1473[5] = _M0L8_2atupleS1610;
  _M0L7_2abindS1473[6] = _M0L8_2atupleS1611;
  _M0L7_2abindS1473[7] = _M0L8_2atupleS1612;
  _M0L7_2abindS1473[8] = _M0L8_2atupleS1613;
  _M0L7_2abindS1473[9] = _M0L8_2atupleS1614;
  _M0L7_2abindS1473[10] = _M0L8_2atupleS1615;
  _M0L7_2abindS1473[11] = _M0L8_2atupleS1616;
  _M0L7_2abindS1473[12] = _M0L8_2atupleS1617;
  _M0L7_2abindS1473[13] = _M0L8_2atupleS1618;
  _M0L7_2abindS1473[14] = _M0L8_2atupleS1619;
  _M0L7_2abindS1473[15] = _M0L8_2atupleS1620;
  _M0L7_2abindS1473[16] = _M0L8_2atupleS1621;
  _M0L7_2abindS1473[17] = _M0L8_2atupleS1622;
  _M0L7_2abindS1473[18] = _M0L8_2atupleS1623;
  _M0L7_2abindS1473[19] = _M0L8_2atupleS1624;
  _M0L7_2abindS1473[20] = _M0L8_2atupleS1625;
  _M0L7_2abindS1473[21] = _M0L8_2atupleS1626;
  _M0L7_2abindS1473[22] = _M0L8_2atupleS1627;
  _M0L6_2atmpS1604 = _M0L7_2abindS1473;
  _M0L6_2atmpS1603
  = (struct _M0TPB9ArrayViewGUsRPB3MapGiUWEuQRPC15error5ErrorNsEEEE){
    .$0 = _M0L6_2atmpS1604, .$1 = 0, .$2 = 23
  };
  #line 399 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0FP26biolab8bio__seq48moonbit__test__driver__internal__no__args__tests
  = _M0MPB3Map3MapGsRPB3MapGiUWEuQRPC15error5ErrorNsEEE(_M0L6_2atmpS1603, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1603.$0);
}

int main(int argc, char** argv) {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** _M0L6_2atmpS1602;
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1567;
  struct _M0TPB5ArrayGUsiEE* _M0L7_2abindS1568;
  int32_t _M0L7_2abindS1569;
  int32_t _M0L2__S1570;
  moonbit_runtime_init(argc, argv);
  moonbit_init();
  _M0L6_2atmpS1602
  = (struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error**)moonbit_empty_ref_array;
  _M0L12async__testsS1567
  = (struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE));
  Moonbit_object_header(_M0L12async__testsS1567)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 116, 0);
  _M0L12async__testsS1567->$0 = _M0L6_2atmpS1602;
  _M0L12async__testsS1567->$1 = 0;
  #line 503 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0L7_2abindS1568
  = _M0FP26biolab8bio__seq52moonbit__test__driver__internal__native__parse__args();
  _M0L7_2abindS1569 = _M0L7_2abindS1568->$1;
  _M0L2__S1570 = 0;
  while (1) {
    if (_M0L2__S1570 < _M0L7_2abindS1569) {
      struct _M0TUsiE** _M0L3bufS1601 = _M0L7_2abindS1568->$0;
      struct _M0TUsiE* _M0L3argS1571 =
        (struct _M0TUsiE*)_M0L3bufS1601[_M0L2__S1570];
      moonbit_string_t _M0L6_2atmpS1598 = _M0L3argS1571->$0;
      int32_t _M0L6_2atmpS1599 = _M0L3argS1571->$1;
      int32_t _M0L6_2atmpS1600;
      moonbit_incref(_M0L6_2atmpS1598);
      #line 504 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
      _M0FP26biolab8bio__seq44moonbit__test__driver__internal__do__execute(_M0L12async__testsS1567, _M0L6_2atmpS1598, _M0L6_2atmpS1599);
      moonbit_decref(_M0L6_2atmpS1598);
      _M0L6_2atmpS1600 = _M0L2__S1570 + 1;
      _M0L2__S1570 = _M0L6_2atmpS1600;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1568);
    }
    break;
  }
  #line 506 "/home/developer/Documents/1/__generated_driver_for_internal_test.mbt"
  _M0IP016_24default__implP26biolab8bio__seq28MoonBit__Async__Test__Driver17run__async__testsGRP26biolab8bio__seq34MoonBit__Async__Test__Driver__ImplE(_M0L12async__testsS1567);
  moonbit_decref(_M0L12async__testsS1567);
  return 0;
}