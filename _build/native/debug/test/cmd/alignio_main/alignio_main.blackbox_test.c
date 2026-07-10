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
struct _M0TPB5EntryGsRPB13StringBuilderE;

struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE;

struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError;

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentERP26biolab8bio__seq12AlignIOErrorE3Err;

struct _M0TPB5ArrayGRPC16string10StringViewE;

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure;

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError;

struct _M0TPB8MutLocalGORPC16string10StringViewE;

struct _M0TWRPC15error5ErrorEs;

struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest;

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__;

struct _M0TWssbEu;

struct _M0TUsiE;

struct _M0TPB13StringBuilder;

struct _M0TPB8MutLocalGRPB13StringBuilderE;

struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360;

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some;

struct _M0TPB5EntryGssE;

struct _M0TUssE;

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq12AlignIOErrorE2Ok;

struct _M0TPB3MapGsRPB5ArrayGiEE;

struct _M0TPB3MapGsRPB5ArrayGsEE;

struct _M0TP26biolab8bio__seq3Seq;

struct _M0BTPB6Logger;

struct _M0DTPC16result6ResultGuRP26biolab8bio__seq12AlignIOErrorE3Err;

struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE;

struct _M0TPB4IterGcE;

struct _M0TPB6Logger;

struct _M0TWcERPC16string10StringView;

struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__;

struct _M0TPB5ArrayGUsiEE;

struct _M0TPB9ArrayViewGUsRPB5ArrayGsEEE;

struct _M0TPB8MutLocalGOsE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok;

struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__;

struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365;

struct _M0TPB5ArrayGiE;

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE;

struct _M0TWEOUsRPB13StringBuilderE;

struct _M0TUsRPB13StringBuilderE;

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError;

struct _M0TWRPC15error5ErrorEu;

struct _M0TURPC16string10StringViewRPB6LoggerE;

struct _M0TPB4IterGRPC16string10StringViewE;

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq12AlignIOErrorE3Err;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd29alignio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE2Ok;

struct _M0TPB5EntryGsRPB5ArrayGiEE;

struct _M0TUsRPB5ArrayGsEE;

struct _M0TPB8MutLocalGiE;

struct _M0TWERPC16option6OptionGRPC16string10StringViewE;

struct _M0TPB3MapGssE;

struct _M0TPB3MapGsRPB13StringBuilderE;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq20MultipleSeqAlignmentRP26biolab8bio__seq12AlignIOErrorE3Err;

struct _M0TPB4Show;

struct _M0TWEOc;

struct _M0TUsRPB5ArrayGiEE;

struct _M0TP26biolab8bio__seq9SeqRecord;

struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__;

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error;

struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE;

struct _M0R38String_3a_3aiter_2eanon__u1894__l256__;

struct _M0TPB5EntryGsRPB5ArrayGsEE;

struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err;

struct _M0BTPB4Show;

struct _M0TWuEu;

struct _M0TPC16string10StringView;

struct _M0KTPB6LoggerTPB13StringBuilder;

struct _M0TPB5ArrayGsE;

struct _M0TPB9ArrayViewGUssEE;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd29alignio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE3Err;

struct _M0TPB4IterGUsRPB13StringBuilderEE;

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE;

struct _M0DTPC16result6ResultGRP26biolab8bio__seq20MultipleSeqAlignmentRP26biolab8bio__seq12AlignIOErrorE2Ok;

struct _M0DTPC16result6ResultGuRP26biolab8bio__seq12AlignIOErrorE2Ok;

struct _M0TP26biolab8bio__seq20MultipleSeqAlignment;

struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError;

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentERP26biolab8bio__seq12AlignIOErrorE2Ok;

struct _M0TPB5EntryGsRPB13StringBuilderE {
  int32_t $0;
  struct _M0TPB5EntryGsRPB13StringBuilderE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TPB13StringBuilder* $5;
  
};

struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE {
  struct _M0TPB5EntryGsRPB13StringBuilderE* $0;
  
};

struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError {
  moonbit_string_t $0;
  
};

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentERP26biolab8bio__seq12AlignIOErrorE3Err {
  void* $0;
  
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

struct _M0TPB8MutLocalGORPC16string10StringViewE {
  void* $0;
  
};

struct _M0TWRPC15error5ErrorEs {
  moonbit_string_t(* code)(struct _M0TWRPC15error5ErrorEs*, void*);
  
};

struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest {
  moonbit_string_t $0;
  
};

struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TWcERPC16string10StringView* $0;
  struct _M0TPB4IterGcE* $1;
  
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

struct _M0TPB8MutLocalGRPB13StringBuilderE {
  struct _M0TPB13StringBuilder* $0;
  
};

struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360 {
  int32_t(* code)(struct _M0TWEOc*);
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TPC16string10StringView {
  moonbit_string_t $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some {
  struct _M0TPC16string10StringView $0;
  
};

struct _M0TPB5EntryGssE {
  int32_t $0;
  struct _M0TPB5EntryGssE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  moonbit_string_t $5;
  
};

struct _M0TUssE {
  moonbit_string_t $0;
  moonbit_string_t $1;
  
};

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq12AlignIOErrorE2Ok {
  moonbit_string_t $0;
  
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

struct _M0TPB3MapGsRPB5ArrayGsEE {
  struct _M0TPB5EntryGsRPB5ArrayGsEE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* $5;
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

struct _M0DTPC16result6ResultGuRP26biolab8bio__seq12AlignIOErrorE3Err {
  void* $0;
  
};

struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE {
  struct _M0TUsRPB13StringBuilderE** $0;
  int32_t $1;
  int32_t $2;
  
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

struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__ {
  struct _M0TUsRPB13StringBuilderE*(* code)(
    struct _M0TWEOUsRPB13StringBuilderE*
  );
  struct _M0TPB8MutLocalGiE* $0;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE* $1;
  
};

struct _M0TPB5ArrayGUsiEE {
  struct _M0TUsiE** $0;
  int32_t $1;
  
};

struct _M0TPB9ArrayViewGUsRPB5ArrayGsEEE {
  struct _M0TUsRPB5ArrayGsEE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0TPB8MutLocalGOsE {
  moonbit_string_t $0;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok {
  int32_t $0;
  
};

struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__ {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  struct _M0TPB8MutLocalGORPC16string10StringViewE* $0;
  struct _M0TPC16string10StringView $1;
  int32_t $2;
  
};

struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365 {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0TPB5ArrayGiE {
  int32_t* $0;
  int32_t $1;
  
};

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** $0;
  int32_t $1;
  
};

struct _M0TWEOUsRPB13StringBuilderE {
  struct _M0TUsRPB13StringBuilderE*(* code)(
    struct _M0TWEOUsRPB13StringBuilderE*
  );
  
};

struct _M0TUsRPB13StringBuilderE {
  moonbit_string_t $0;
  struct _M0TPB13StringBuilder* $1;
  
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

struct _M0TPB4IterGRPC16string10StringViewE {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* $0;
  int64_t $1;
  
};

struct _M0DTPC16result6ResultGsRP26biolab8bio__seq12AlignIOErrorE3Err {
  void* $0;
  
};

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd29alignio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE2Ok {
  int32_t $0;
  
};

struct _M0TPB5EntryGsRPB5ArrayGiEE {
  int32_t $0;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TPB5ArrayGiE* $5;
  
};

struct _M0TUsRPB5ArrayGsEE {
  moonbit_string_t $0;
  struct _M0TPB5ArrayGsE* $1;
  
};

struct _M0TPB8MutLocalGiE {
  int32_t $0;
  
};

struct _M0TWERPC16option6OptionGRPC16string10StringViewE {
  void*(* code)(struct _M0TWERPC16option6OptionGRPC16string10StringViewE*);
  
};

struct _M0TPB3MapGssE {
  struct _M0TPB5EntryGssE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGssE* $5;
  int32_t $6;
  
};

struct _M0TPB3MapGsRPB13StringBuilderE {
  struct _M0TPB5EntryGsRPB13StringBuilderE** $0;
  int32_t $1;
  int32_t $2;
  int32_t $3;
  int32_t $4;
  struct _M0TPB5EntryGsRPB13StringBuilderE* $5;
  int32_t $6;
  
};

struct _M0DTPC16result6ResultGRP26biolab8bio__seq20MultipleSeqAlignmentRP26biolab8bio__seq12AlignIOErrorE3Err {
  void* $0;
  
};

struct _M0TPB4Show {
  struct _M0BTPB4Show* $0;
  void* $1;
  
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

struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  int32_t $1;
  struct _M0TPC16string10StringView $2;
  
};

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error {
  struct moonbit_result_0(* code)(
    struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error*,
    struct _M0TWuEu*,
    struct _M0TWRPC15error5ErrorEu*
  );
  
};

struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE {
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** $0;
  int32_t $1;
  
};

struct _M0R38String_3a_3aiter_2eanon__u1894__l256__ {
  int32_t(* code)(struct _M0TWEOc*);
  struct _M0TPB8MutLocalGiE* $0;
  moonbit_string_t $1;
  int32_t $2;
  
};

struct _M0TPB5EntryGsRPB5ArrayGsEE {
  int32_t $0;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* $1;
  int32_t $2;
  int32_t $3;
  moonbit_string_t $4;
  struct _M0TPB5ArrayGsE* $5;
  
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

struct _M0KTPB6LoggerTPB13StringBuilder {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0TPB5ArrayGsE {
  moonbit_string_t* $0;
  int32_t $1;
  
};

struct _M0TPB9ArrayViewGUssEE {
  struct _M0TUssE** $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd29alignio__main__blackbox__test33MoonBitTestDriverInternalSkipTestE3Err {
  void* $0;
  
};

struct _M0TPB4IterGUsRPB13StringBuilderEE {
  struct _M0TWEOUsRPB13StringBuilderE* $0;
  int64_t $1;
  
};

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE {
  struct _M0TP26biolab8bio__seq9SeqRecord** $0;
  int32_t $1;
  
};

struct _M0DTPC16result6ResultGRP26biolab8bio__seq20MultipleSeqAlignmentRP26biolab8bio__seq12AlignIOErrorE2Ok {
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* $0;
  
};

struct _M0DTPC16result6ResultGuRP26biolab8bio__seq12AlignIOErrorE2Ok {
  int32_t $0;
  
};

struct _M0TP26biolab8bio__seq20MultipleSeqAlignment {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* $0;
  struct _M0TPB3MapGssE* $1;
  struct _M0TPB3MapGsRPB5ArrayGsEE* $2;
  
};

struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError {
  moonbit_string_t $0;
  
};

struct _M0DTPC16result6ResultGRPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentERP26biolab8bio__seq12AlignIOErrorE2Ok {
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* $0;
  
};

struct moonbit_result_3 {
  int tag;
  union {
    struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* ok;
    void* err;
    
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

struct moonbit_result_1 {
  int tag;
  union { struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* ok; void* err; 
  } data;
  
};

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1372(
  struct _M0TWRPC15error5ErrorEs*,
  void*
);

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1365(
  struct _M0TWssbEu*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1360(
  struct _M0TWEOc*
);

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__args(
  
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1338(
  int32_t,
  moonbit_string_t,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1331(
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324(
  int32_t,
  moonbit_bytes_t
);

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1317(
  int32_t,
  moonbit_string_t
);

#define _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__get__cli__args__ffi moonbit_get_cli_args

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEOc*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*
);

int32_t _M0FP46biolab8bio__seq3cmd13alignio__main13alignio__test();

struct moonbit_result_0 _M0FP46biolab8bio__seq3cmd13alignio__main18run__alignio__test(
  
);

int32_t _M0FP46biolab8bio__seq3cmd13alignio__main4emit(
  moonbit_string_t,
  moonbit_string_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main19alignment__to__json(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main17records__to__json(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main16record__to__json(
  struct _M0TP26biolab8bio__seq9SeqRecord*
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(
  moonbit_string_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main12json__escape(
  moonbit_string_t
);

moonbit_string_t _M0MP26biolab8bio__seq20MultipleSeqAlignment6column(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*,
  int32_t
);

int32_t _M0MP26biolab8bio__seq20MultipleSeqAlignment12num__records(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

struct moonbit_result_2 _M0FP26biolab8bio__seq14alignio__write(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*,
  moonbit_string_t
);

moonbit_string_t _M0FP26biolab8bio__seq35write__phylip__sequential__multiple(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*
);

moonbit_string_t _M0FP26biolab8bio__seq23write__phylip__multiple(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*
);

moonbit_string_t _M0FP26biolab8bio__seq27write__fasta__as__alignment(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*
);

moonbit_string_t _M0FP26biolab8bio__seq24write__clustal__multiple(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*
);

struct moonbit_result_1 _M0FP26biolab8bio__seq13alignio__read(
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_3 _M0FP26biolab8bio__seq14alignio__parse(
  moonbit_string_t,
  moonbit_string_t
);

struct moonbit_result_1 _M0FP26biolab8bio__seq27parse__fasta__as__alignment(
  moonbit_string_t
);

moonbit_string_t _M0FP26biolab8bio__seq14write__clustal(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

struct moonbit_result_1 _M0FP26biolab8bio__seq14parse__clustal(
  moonbit_string_t
);

moonbit_string_t _M0FP26biolab8bio__seq25write__phylip__sequential(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

moonbit_string_t _M0FP26biolab8bio__seq13write__phylip(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

int32_t _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*
);

struct moonbit_result_1 _M0FP26biolab8bio__seq25parse__phylip__sequential(
  moonbit_string_t
);

struct moonbit_result_1 _M0FP26biolab8bio__seq13parse__phylip(
  moonbit_string_t
);

struct moonbit_result_1 _M0FP26biolab8bio__seq22parse__phylip__generic(
  moonbit_string_t,
  int32_t
);

struct moonbit_result_1 _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TPB3MapGssE*,
  struct _M0TPB3MapGsRPB5ArrayGsEE*
);

struct moonbit_result_1 _M0MP26biolab8bio__seq20MultipleSeqAlignment11new_2einner(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  struct _M0TPB3MapGssE*,
  struct _M0TPB3MapGsRPB5ArrayGsEE*
);

int64_t _M0FP26biolab8bio__seq10parse__int(moonbit_string_t);

moonbit_string_t _M0FP26biolab8bio__seq12write__fasta(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

moonbit_string_t _M0FP26biolab8bio__seq19build__fasta__title(
  moonbit_string_t,
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

int32_t _M0MP26biolab8bio__seq3Seq9get__char(
  struct _M0TP26biolab8bio__seq3Seq*,
  int32_t
);

int32_t _M0MP26biolab8bio__seq3Seq6length(struct _M0TP26biolab8bio__seq3Seq*);

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq*
);

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t
);

moonbit_string_t _M0MPC15array5Array2atGsE(struct _M0TPB5ArrayGsE*, int32_t);

struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0MPC15array5Array2atGRP26biolab8bio__seq20MultipleSeqAlignmentE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*,
  int32_t
);

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*,
  int32_t
);

struct _M0TPC16string10StringView _M0MPC15array5Array2atGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  int32_t
);

int32_t _M0FPB7printlnGsE(moonbit_string_t);

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t);

struct _M0TUsRPB13StringBuilderE* _M0MPB5Iter24nextGsRPB13StringBuilderE(
  struct _M0TPB4IterGUsRPB13StringBuilderEE*
);

struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0MPB3Map5iter2GsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*
);

struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0MPB3Map4iterGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*
);

struct _M0TUsRPB13StringBuilderE* _M0MPB3Map4iterGsRPB13StringBuilderEC2286l711(
  struct _M0TWEOUsRPB13StringBuilderE*
);

struct _M0TPB13StringBuilder* _M0MPB3Map3getGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  moonbit_string_t
);

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0MPB3Map3MapGsRPB5ArrayGiEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE,
  int64_t
);

struct _M0TPB3MapGssE* _M0MPB3Map3MapGssE(
  struct _M0TPB9ArrayViewGUssEE,
  int64_t
);

struct _M0TPB3MapGsRPB5ArrayGsEE* _M0MPB3Map3MapGsRPB5ArrayGsEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGsEEE,
  int64_t
);

struct _M0TPB3MapGsRPB13StringBuilderE* _M0MPB3Map3MapGsRPB13StringBuilderE(
  struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE,
  int64_t
);

int32_t _M0MPB3Map3setGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGiE*
);

int32_t _M0MPB3Map3setGssE(
  struct _M0TPB3MapGssE*,
  moonbit_string_t,
  moonbit_string_t
);

int32_t _M0MPB3Map3setGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGsE*
);

int32_t _M0MPB3Map3setGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  moonbit_string_t,
  struct _M0TPB13StringBuilder*
);

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGiE*,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGssE(
  struct _M0TPB3MapGssE*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE*,
  moonbit_string_t,
  struct _M0TPB5ArrayGsE*,
  int32_t
);

int32_t _M0MPB3Map15set__with__hashGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  moonbit_string_t,
  struct _M0TPB13StringBuilder*,
  int32_t
);

int32_t _M0MPB3Map4growGsRPB5ArrayGiEE(struct _M0TPB3MapGsRPB5ArrayGiEE*);

int32_t _M0MPB3Map4growGssE(struct _M0TPB3MapGssE*);

int32_t _M0MPB3Map4growGsRPB5ArrayGsEE(struct _M0TPB3MapGsRPB5ArrayGsEE*);

int32_t _M0MPB3Map4growGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map20rehash__place__entryGssE(
  struct _M0TPB3MapGssE*,
  struct _M0TPB5EntryGssE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE*,
  struct _M0TPB5EntryGsRPB5ArrayGsEE*
);

int32_t _M0MPB3Map20rehash__place__entryGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  struct _M0TPB5EntryGsRPB13StringBuilderE*
);

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map10push__awayGssE(
  struct _M0TPB3MapGssE*,
  int32_t,
  struct _M0TPB5EntryGssE*
);

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGsEE*
);

int32_t _M0MPB3Map10push__awayGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  int32_t,
  struct _M0TPB5EntryGsRPB13StringBuilderE*
);

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGssE(
  struct _M0TPB3MapGssE*,
  struct _M0TPB5EntryGssE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE*,
  struct _M0TPB5EntryGsRPB5ArrayGsEE*,
  int32_t
);

int32_t _M0MPB3Map10set__entryGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  struct _M0TPB5EntryGsRPB13StringBuilderE*,
  int32_t
);

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

int32_t _M0MPB3Map20add__entry__to__tailGssE(
  struct _M0TPB3MapGssE*,
  int32_t,
  struct _M0TPB5EntryGssE*
);

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE*,
  int32_t,
  struct _M0TPB5EntryGsRPB5ArrayGsEE*
);

int32_t _M0MPB3Map20add__entry__to__tailGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE*,
  int32_t,
  struct _M0TPB5EntryGsRPB13StringBuilderE*
);

int32_t _M0MPC13int3Int3max(int32_t, int32_t);

int32_t _M0FPB21capacity__for__length(int32_t);

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0FPB8new__mapGsRPB5ArrayGiEE(int32_t);

struct _M0TPB3MapGssE* _M0FPB8new__mapGssE(int32_t);

struct _M0TPB3MapGsRPB5ArrayGsEE* _M0FPB8new__mapGsRPB5ArrayGsEE(int32_t);

struct _M0TPB3MapGsRPB13StringBuilderE* _M0FPB8new__mapGsRPB13StringBuilderE(
  int32_t
);

int32_t _M0MPC13int3Int20next__power__of__two(int32_t);

int32_t _M0FPB21calc__grow__threshold(int32_t);

struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(
  struct _M0TPB5EntryGsRPB5ArrayGiEE*
);

struct _M0TPB5EntryGssE* _M0MPC16option6Option6unwrapGRPB5EntryGssEE(
  struct _M0TPB5EntryGssE*
);

struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGsEEE(
  struct _M0TPB5EntryGsRPB5ArrayGsEE*
);

struct _M0TPB5EntryGsRPB13StringBuilderE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB13StringBuilderEE(
  struct _M0TPB5EntryGsRPB13StringBuilderE*
);

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t);

int32_t _M0MPC16string6String9get__char(moonbit_string_t, int32_t);

struct _M0TPC16string10StringView _M0MPC16string10StringView7replace(
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView,
  struct _M0TPC16string10StringView
);

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

void* _M0MPC16string10StringView5splitC1923l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC1919l1145(
  struct _M0TWcERPC16string10StringView*,
  int32_t
);

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t);

moonbit_string_t _M0FPB16char__to__string(int32_t);

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE*,
  struct _M0TWcERPC16string10StringView*
);

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC1912l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*
);

struct _M0TPB4IterGcE* _M0MPC16string6String4iter(moonbit_string_t);

int32_t _M0MPC16string6String4iterC1894l256(struct _M0TWEOc*);

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

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array7reallocGUsiEE(struct _M0TPB5ArrayGUsiEE*);

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

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

int32_t _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*,
  int32_t
);

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq20MultipleSeqAlignmentE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*
);

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
);

int32_t _M0MPC15array5Array6lengthGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE*
);

int32_t _M0MPC15array5Array6lengthGsE(struct _M0TPB5ArrayGsE*);

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(int32_t);

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(
  int32_t
);

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

int32_t _M0MPB13StringBuilder11write__iter(
  struct _M0TPB13StringBuilder*,
  struct _M0TPB4IterGcE*
);

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t,
  int32_t,
  int64_t
);

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView
);

int32_t _M0MPC16string10StringView4iterC1772l219(struct _M0TWEOc*);

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

struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0MPC15array5Array6bufferGRP26biolab8bio__seq20MultipleSeqAlignmentE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*
);

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*
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

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE*,
  int64_t
);

struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0MPB4Iter3newGUsRPB13StringBuilderEE(
  struct _M0TWEOUsRPB13StringBuilderE*,
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

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE*
);

struct _M0TUsRPB13StringBuilderE* _M0MPB4Iter4nextGUsRPB13StringBuilderEE(
  struct _M0TPB4IterGUsRPB13StringBuilderEE*
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
} const moonbit_string_literal_50 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 44, 32, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[112]; 
} const moonbit_string_literal_92 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 111, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 97, 108, 105, 103, 110, 105, 111, 95, 109, 97, 105, 110, 
    95, 98, 108, 97, 99, 107, 98, 111, 120, 95, 116, 101, 115, 116, 46, 
    77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 
    118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 74, 115, 69, 
    114, 114, 111, 114, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 
    115, 116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 
    97, 108, 74, 115, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[41]; 
} const moonbit_string_literal_91 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 40, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 46, 65, 108, 
    105, 103, 110, 73, 79, 69, 114, 114, 111, 114, 46, 65, 108, 105, 
    103, 110, 73, 79, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_66 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 69, 109, 
    112, 116, 121, 32, 80, 72, 89, 76, 73, 80, 32, 102, 105, 108, 101, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_36 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 123, 34, 
    110, 117, 109, 95, 114, 101, 99, 111, 114, 100, 115, 34, 58, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_43 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 123, 34, 
    101, 114, 114, 111, 114, 34, 58, 32, 110, 117, 108, 108, 125, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_78 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 117, 114, 114, 111, 103, 97, 116, 
    101, 32, 112, 97, 105, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[53]; 
} const moonbit_string_literal_90 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 52, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[114]; 
} const moonbit_string_literal_89 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 113, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 97, 108, 105, 103, 110, 105, 111, 95, 109, 97, 105, 110, 
    95, 98, 108, 97, 99, 107, 98, 111, 120, 95, 116, 101, 115, 116, 46, 
    77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 
    118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 83, 107, 105, 
    112, 84, 101, 115, 116, 46, 77, 111, 111, 110, 66, 105, 116, 84, 
    101, 115, 116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 
    110, 97, 108, 83, 107, 105, 112, 84, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_59 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 110, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_25 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 65, 70, 
    49, 57, 49, 54, 54, 48, 46, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[31]; 
} const moonbit_string_literal_72 =
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

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_37 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 44, 32, 
    34, 97, 108, 105, 103, 110, 109, 101, 110, 116, 95, 108, 101, 110, 
    103, 116, 104, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_38 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 97, 108, 
    105, 103, 110, 109, 101, 110, 116, 95, 108, 101, 110, 103, 116, 104, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_73 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_46 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 61, 61, 
    61, 32, 67, 65, 83, 69, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_13 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 99, 108, 
    117, 115, 116, 97, 108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_84 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 44, 32, 
    108, 101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[6]; 
} const moonbit_string_literal_19 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 5, 102, 97, 
    115, 116, 97, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_81 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 111, 
    117, 110, 100, 115, 32, 99, 104, 101, 99, 107, 32, 102, 97, 105, 
    108, 101, 100, 58, 32, 97, 108, 108, 111, 99, 97, 116, 101, 95, 108, 
    101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_49 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 44, 32, 
    34, 114, 101, 99, 111, 114, 100, 115, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_41 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 97, 108, 
    105, 103, 110, 109, 101, 110, 116, 95, 99, 111, 108, 117, 109, 110, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[48]; 
} const moonbit_string_literal_21 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 47, 84, 65, 
    84, 65, 67, 65, 84, 84, 65, 65, 65, 71, 65, 65, 71, 71, 71, 71, 71, 
    65, 84, 71, 67, 71, 71, 65, 84, 65, 65, 65, 84, 71, 71, 65, 65, 65, 
    71, 71, 67, 71, 65, 65, 65, 71, 65, 71, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[56]; 
} const moonbit_string_literal_79 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 55, 105, 110, 
    118, 97, 108, 105, 100, 32, 115, 116, 97, 114, 116, 32, 111, 114, 
    32, 101, 110, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 114, 
    32, 83, 116, 114, 105, 110, 103, 58, 58, 99, 111, 100, 101, 112, 
    111, 105, 110, 116, 95, 108, 101, 110, 103, 116, 104, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_31 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 119, 114, 
    105, 116, 101, 95, 112, 104, 121, 108, 105, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[17]; 
} const moonbit_string_literal_61 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 16, 85, 110, 
    107, 110, 111, 119, 110, 32, 102, 111, 114, 109, 97, 116, 58, 32, 
    0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[175]; 
} const moonbit_string_literal_12 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 174, 62, 115, 
    101, 113, 49, 32, 68, 101, 115, 99, 114, 105, 112, 116, 105, 111, 
    110, 32, 49, 10, 65, 84, 71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 
    67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 
    71, 65, 84, 67, 71, 65, 84, 67, 10, 62, 115, 101, 113, 50, 32, 68, 
    101, 115, 99, 114, 105, 112, 116, 105, 111, 110, 32, 50, 10, 65, 
    84, 71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 
    65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 
    84, 67, 10, 62, 115, 101, 113, 51, 32, 68, 101, 115, 99, 114, 105, 
    112, 116, 105, 111, 110, 32, 51, 10, 65, 84, 71, 67, 67, 71, 65, 
    84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 
    67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_87 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 41, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[15]; 
} const moonbit_string_literal_40 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 14, 44, 32, 
    34, 99, 111, 108, 117, 109, 110, 95, 53, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_74 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 48, 49, 
    50, 51, 52, 53, 54, 55, 56, 57, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 119, 120, 121, 122, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_65 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 32, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[78]; 
} const moonbit_string_literal_44 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 77, 123, 34, 
    101, 114, 114, 111, 114, 34, 58, 32, 34, 86, 97, 108, 117, 101, 69, 
    114, 114, 111, 114, 34, 44, 32, 34, 109, 101, 115, 115, 97, 103, 
    101, 34, 58, 32, 34, 77, 111, 114, 101, 32, 116, 104, 97, 110, 32, 
    111, 110, 101, 32, 97, 108, 105, 103, 110, 109, 101, 110, 116, 32, 
    102, 111, 117, 110, 100, 32, 105, 110, 32, 104, 97, 110, 100, 108, 
    101, 34, 125, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[51]; 
} const moonbit_string_literal_88 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 50, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 73, 110, 115, 112, 101, 
    99, 116, 69, 114, 114, 111, 114, 46, 73, 110, 115, 112, 101, 99, 
    116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_85 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 46, 108, 101, 110, 103, 116, 104, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_42 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 10, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_53 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 44, 32, 
    34, 100, 101, 115, 99, 114, 105, 112, 116, 105, 111, 110, 34, 58, 
    32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_9 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 84, 101, 
    115, 116, 32, 101, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_7 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 45, 45, 
    45, 45, 45, 32, 69, 78, 68, 32, 77, 79, 79, 78, 32, 84, 69, 83, 84, 
    32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_82 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_30 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 115, 101, 
    113, 51, 0
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
} const moonbit_string_literal_48 =
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

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_39 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 123, 34, 
    99, 111, 108, 117, 109, 110, 95, 48, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_52 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 44, 32, 
    34, 110, 97, 109, 101, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_86 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 70, 97, 
    105, 108, 117, 114, 101, 40, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[40]; 
} const moonbit_string_literal_63 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 39, 77, 111, 
    114, 101, 32, 116, 104, 97, 110, 32, 111, 110, 101, 32, 97, 108, 
    105, 103, 110, 109, 101, 110, 116, 32, 102, 111, 117, 110, 100, 32, 
    105, 110, 32, 104, 97, 110, 100, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[8]; 
} const moonbit_string_literal_51 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 7, 123, 34, 
    105, 100, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_70 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 114, 101, 
    112, 101, 97, 116, 32, 114, 101, 115, 117, 108, 116, 32, 116, 111, 
    111, 32, 108, 97, 114, 103, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_6 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 125, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_17 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 112, 104, 
    121, 108, 105, 112, 45, 115, 101, 113, 117, 101, 110, 116, 105, 97, 
    108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[35]; 
} const moonbit_string_literal_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 34, 45, 45, 
    45, 45, 45, 32, 66, 69, 71, 73, 78, 32, 77, 79, 79, 78, 32, 84, 69, 
    83, 84, 32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_75 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 114, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[56]; 
} const moonbit_string_literal_68 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 55, 65, 108, 
    108, 32, 115, 101, 113, 117, 101, 110, 99, 101, 115, 32, 109, 117, 
    115, 116, 32, 104, 97, 118, 101, 32, 116, 104, 101, 32, 115, 97, 
    109, 101, 32, 108, 101, 110, 103, 116, 104, 32, 105, 110, 32, 97, 
    110, 32, 97, 108, 105, 103, 110, 109, 101, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_60 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 116, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_83 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    100, 115, 116, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_80 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 105, 110, 
    118, 97, 108, 105, 100, 32, 99, 111, 100, 101, 32, 112, 111, 105, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_33 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 115, 101, 
    113, 49, 32, 100, 101, 115, 99, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_16 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 112, 97, 
    114, 115, 101, 95, 112, 104, 121, 108, 105, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[13]; 
} const moonbit_string_literal_55 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 12, 44, 32, 
    34, 108, 101, 110, 103, 116, 104, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[20]; 
} const moonbit_string_literal_45 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 19, 114, 101, 
    97, 100, 95, 109, 117, 108, 116, 105, 112, 108, 101, 95, 101, 114, 
    114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_22 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 65, 70, 
    49, 57, 49, 54, 53, 57, 46, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_5 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 44, 34, 
    109, 101, 115, 115, 97, 103, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[48]; 
} const moonbit_string_literal_64 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 47, 67, 76, 
    85, 83, 84, 65, 76, 32, 88, 32, 40, 49, 46, 56, 49, 41, 32, 109, 
    117, 108, 116, 105, 112, 108, 101, 32, 115, 101, 113, 117, 101, 110, 
    99, 101, 32, 97, 108, 105, 103, 110, 109, 101, 110, 116, 10, 10, 
    10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_35 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 119, 114, 
    105, 116, 101, 95, 102, 97, 115, 116, 97, 95, 97, 108, 105, 103, 
    110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_28 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 115, 101, 
    113, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[23]; 
} const moonbit_string_literal_71 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 22, 73, 110, 
    118, 97, 108, 105, 100, 32, 105, 110, 100, 101, 120, 32, 102, 111, 
    114, 32, 86, 105, 101, 119, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[18]; 
} const moonbit_string_literal_20 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 17, 112, 97, 
    114, 115, 101, 95, 102, 97, 115, 116, 97, 95, 97, 108, 105, 103, 
    110, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[11]; 
} const moonbit_string_literal_23 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 10, 65, 70, 
    49, 57, 49, 54, 53, 56, 46, 49, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[21]; 
} const moonbit_string_literal_27 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 20, 65, 84, 
    71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 
    84, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_34 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 115, 101, 
    113, 50, 32, 100, 101, 115, 99, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[29]; 
} const moonbit_string_literal_62 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 28, 78, 111, 
    32, 97, 108, 105, 103, 110, 109, 101, 110, 116, 32, 102, 111, 117, 
    110, 100, 32, 105, 110, 32, 104, 97, 110, 100, 108, 101, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[38]; 
} const moonbit_string_literal_32 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 37, 65, 84, 
    71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 
    84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 
    67, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_77 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 92, 117, 123, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_76 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_56 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 34, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_54 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 32, 
    34, 115, 101, 113, 34, 58, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[7]; 
} const moonbit_string_literal_15 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 6, 112, 104, 
    121, 108, 105, 112, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_14 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 112, 97, 
    114, 115, 101, 95, 99, 108, 117, 115, 116, 97, 108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_69 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 110, 101, 
    103, 97, 116, 105, 118, 101, 32, 114, 101, 112, 101, 97, 116, 32, 
    99, 111, 117, 110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[22]; 
} const moonbit_string_literal_67 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 21, 73, 110, 
    118, 97, 108, 105, 100, 32, 80, 72, 89, 76, 73, 80, 32, 104, 101, 
    97, 100, 101, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_57 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 92, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_29 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 115, 101, 
    113, 50, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[5]; 
} const moonbit_string_literal_47 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 4, 32, 61, 
    61, 61, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[14]; 
} const moonbit_string_literal_26 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 13, 119, 114, 
    105, 116, 101, 95, 99, 108, 117, 115, 116, 97, 108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[48]; 
} const moonbit_string_literal_24 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 47, 84, 65, 
    84, 65, 67, 65, 84, 65, 65, 65, 65, 71, 65, 65, 71, 71, 71, 71, 71, 
    65, 84, 71, 67, 71, 71, 65, 84, 65, 65, 65, 84, 71, 71, 65, 65, 65, 
    71, 71, 67, 71, 65, 65, 65, 71, 65, 71, 65, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[372]; 
} const moonbit_string_literal_10 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 371, 67, 76, 
    85, 83, 84, 65, 76, 32, 88, 32, 40, 49, 46, 56, 49, 41, 32, 109, 
    117, 108, 116, 105, 112, 108, 101, 32, 115, 101, 113, 117, 101, 110, 
    99, 101, 32, 97, 108, 105, 103, 110, 109, 101, 110, 116, 10, 10, 
    65, 70, 49, 57, 49, 54, 53, 57, 46, 49, 32, 32, 32, 32, 32, 32, 84, 
    65, 84, 65, 67, 65, 84, 84, 65, 65, 65, 71, 65, 65, 71, 71, 71, 71, 
    71, 65, 84, 71, 67, 71, 71, 65, 84, 65, 65, 65, 84, 71, 71, 65, 65, 
    65, 71, 71, 67, 71, 65, 65, 65, 71, 10, 65, 70, 49, 57, 49, 54, 53, 
    56, 46, 49, 32, 32, 32, 32, 32, 32, 84, 65, 84, 65, 67, 65, 84, 84, 
    65, 65, 65, 71, 65, 65, 71, 71, 71, 71, 71, 65, 84, 71, 67, 71, 71, 
    65, 84, 65, 65, 65, 84, 71, 71, 65, 65, 65, 71, 71, 67, 71, 65, 65, 
    65, 71, 10, 65, 70, 49, 57, 49, 54, 54, 49, 46, 49, 32, 32, 32, 32, 
    32, 32, 84, 65, 84, 65, 67, 65, 84, 84, 65, 65, 65, 71, 65, 65, 71, 
    71, 71, 71, 71, 65, 84, 71, 67, 71, 71, 65, 84, 65, 65, 65, 84, 71, 
    71, 65, 65, 65, 71, 71, 67, 71, 65, 65, 65, 71, 10, 65, 70, 49, 57, 
    49, 54, 54, 48, 46, 49, 32, 32, 32, 32, 32, 32, 84, 65, 84, 65, 67, 
    65, 84, 65, 65, 65, 65, 71, 65, 65, 71, 71, 71, 71, 71, 65, 84, 71, 
    67, 71, 71, 65, 84, 65, 65, 65, 84, 71, 71, 65, 65, 65, 71, 71, 67, 
    71, 65, 65, 65, 71, 10, 10, 65, 70, 49, 57, 49, 54, 53, 57, 46, 49, 
    32, 32, 32, 32, 32, 32, 65, 71, 65, 10, 65, 70, 49, 57, 49, 54, 53, 
    56, 46, 49, 32, 32, 32, 32, 32, 32, 65, 71, 65, 10, 65, 70, 49, 57, 
    49, 54, 54, 49, 46, 49, 32, 32, 32, 32, 32, 32, 65, 71, 65, 10, 65, 
    70, 49, 57, 49, 54, 54, 48, 46, 49, 32, 32, 32, 32, 32, 32, 65, 71, 
    65, 10, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[24]; 
} const moonbit_string_literal_18 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 23, 112, 97, 
    114, 115, 101, 95, 112, 104, 121, 108, 105, 112, 95, 115, 101, 113, 
    117, 101, 110, 116, 105, 97, 108, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[105]; 
} const moonbit_string_literal_11 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 104, 51, 32, 
    50, 48, 10, 115, 101, 113, 49, 32, 32, 32, 32, 32, 32, 32, 32, 65, 
    84, 71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 
    65, 84, 10, 115, 101, 113, 50, 32, 32, 32, 32, 32, 32, 32, 32, 65, 
    84, 71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 
    65, 84, 10, 115, 101, 113, 51, 32, 32, 32, 32, 32, 32, 32, 32, 65, 
    84, 71, 67, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 65, 84, 67, 71, 
    65, 84, 10, 0
  };

struct moonbit_object const moonbit_constant_constructor_0 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0)
  };

struct {
  int32_t rc;
  uint32_t meta;
  struct _M0TWcERPC16string10StringView data;
  
} const _M0MPC16string10StringView5splitC1919l1145$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0MPC16string10StringView5splitC1919l1145
  };

struct { int32_t rc; uint32_t meta; struct _M0TWRPC15error5ErrorEs data; 
} const _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1372$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1372
  };

uint32_t const moonbit_layout_table_data[147] =
  {
    sizeof(struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360)
    / 4, 1,
    offsetof(struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360, $1)
    / 4,
    sizeof(struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365)
    / 4, 1,
    offsetof(struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365, $1)
    / 4,
    sizeof(struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest, $0)
    / 4, sizeof(struct _M0TPB5ArrayGUsiEE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGUsiEE, $0) / 4, sizeof(struct _M0TUsiE) / 4,
    1, offsetof(struct _M0TUsiE, $0) / 4, sizeof(struct _M0TPB5ArrayGsE) / 4,
    1, offsetof(struct _M0TPB5ArrayGsE, $0) / 4,
    sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE, $0) / 4,
    sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE)
    / 4, 1,
    offsetof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE, $0)
    / 4,
    sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError, $0)
    / 4, sizeof(struct _M0TP26biolab8bio__seq9SeqRecord) / 4, 5,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $0) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $1) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $2) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $3) / 4,
    offsetof(struct _M0TP26biolab8bio__seq9SeqRecord, $4) / 4,
    sizeof(struct _M0TP26biolab8bio__seq20MultipleSeqAlignment) / 4, 
    3, offsetof(struct _M0TP26biolab8bio__seq20MultipleSeqAlignment, $0) / 4,
    offsetof(struct _M0TP26biolab8bio__seq20MultipleSeqAlignment, $1) / 4,
    offsetof(struct _M0TP26biolab8bio__seq20MultipleSeqAlignment, $2) / 4,
    sizeof(struct _M0TPB8MutLocalGOsE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGOsE, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGRPB13StringBuilderE) / 4, 1,
    offsetof(struct _M0TPB8MutLocalGRPB13StringBuilderE, $0) / 4,
    sizeof(struct _M0TP26biolab8bio__seq3Seq) / 4, 1,
    offsetof(struct _M0TP26biolab8bio__seq3Seq, $0) / 4,
    sizeof(struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE) / 4, 
    1,
    offsetof(struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE, $0) / 4,
    sizeof(struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__)
    / 4, 2,
    offsetof(struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__, $0)
    / 4,
    offsetof(struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__, $1)
    / 4, sizeof(struct _M0TUsRPB13StringBuilderE) / 4, 2,
    offsetof(struct _M0TUsRPB13StringBuilderE, $0) / 4,
    offsetof(struct _M0TUsRPB13StringBuilderE, $1) / 4,
    sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE) / 4, 3,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGiEE, $5) / 4,
    sizeof(struct _M0TPB5EntryGssE) / 4, 3,
    offsetof(struct _M0TPB5EntryGssE, $1) / 4,
    offsetof(struct _M0TPB5EntryGssE, $4) / 4,
    offsetof(struct _M0TPB5EntryGssE, $5) / 4,
    sizeof(struct _M0TPB5EntryGsRPB5ArrayGsEE) / 4, 3,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGsEE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGsEE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRPB5ArrayGsEE, $5) / 4,
    sizeof(struct _M0TPB5EntryGsRPB13StringBuilderE) / 4, 3,
    offsetof(struct _M0TPB5EntryGsRPB13StringBuilderE, $1) / 4,
    offsetof(struct _M0TPB5EntryGsRPB13StringBuilderE, $4) / 4,
    offsetof(struct _M0TPB5EntryGsRPB13StringBuilderE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRPB5ArrayGiEE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGiEE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGiEE, $5) / 4,
    sizeof(struct _M0TPB3MapGssE) / 4, 2,
    offsetof(struct _M0TPB3MapGssE, $0) / 4,
    offsetof(struct _M0TPB3MapGssE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRPB5ArrayGsEE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGsEE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB5ArrayGsEE, $5) / 4,
    sizeof(struct _M0TPB3MapGsRPB13StringBuilderE) / 4, 2,
    offsetof(struct _M0TPB3MapGsRPB13StringBuilderE, $0) / 4,
    offsetof(struct _M0TPB3MapGsRPB13StringBuilderE, $5) / 4,
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
    sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__) / 4, 
    2,
    offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__, $0)
    / 4,
    (offsetof(struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__, $1)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4,
    sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__)
    / 4, 2,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__, $0)
    / 4,
    offsetof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__, $1)
    / 4, sizeof(struct _M0TPB4IterGRPC16string10StringViewE) / 4, 1,
    offsetof(struct _M0TPB4IterGRPC16string10StringViewE, $0) / 4,
    sizeof(struct _M0R38String_3a_3aiter_2eanon__u1894__l256__) / 4, 
    2, offsetof(struct _M0R38String_3a_3aiter_2eanon__u1894__l256__, $0) / 4,
    offsetof(struct _M0R38String_3a_3aiter_2eanon__u1894__l256__, $1) / 4,
    sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__) / 4, 
    2,
    offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__, $0) / 4,
    (offsetof(struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__, $2)
     + offsetof(struct _M0TPC16string10StringView, $0))
    / 4, sizeof(struct _M0TPB4IterGcE) / 4, 1,
    offsetof(struct _M0TPB4IterGcE, $0) / 4,
    sizeof(struct _M0TPB4IterGUsRPB13StringBuilderEE) / 4, 1,
    offsetof(struct _M0TPB4IterGUsRPB13StringBuilderEE, $0) / 4,
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

int64_t _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9980GUsRPB13StringBuilderEE = 0ll;

int64_t _M0MPB4Iter4nextN6constrS9981GUsRPB13StringBuilderEE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GcE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE = 0ll;

int64_t _M0MPB4Iter3newN6constrS9988GUsRPB13StringBuilderEE = 0ll;

int64_t _M0FPB28boyer__moore__horspool__findN6constrS9990 = 0ll;

int64_t _M0FPB18brute__force__findN6constrS9991 = 0ll;

int32_t _M0FP26biolab8bio__seq18fasta__line__width = 60;

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1393,
  moonbit_string_t _M0L8filenameS1362,
  int32_t _M0L5indexS1364
) {
  struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360* _closure_3246;
  struct _M0TWEOc* _M0L13handle__startS1360;
  struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365* _closure_3247;
  struct _M0TWssbEu* _M0L14handle__resultS1365;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS1372;
  void* _M0L11_2atry__errS1387;
  struct moonbit_result_0 _tmp_3249;
  int32_t _handle__error__result_3250;
  int32_t _M0L6_2atmpS2931;
  void* _M0L3errS1388;
  moonbit_string_t _M0L4nameS1390;
  struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest* _M0L36_2aMoonBitTestDriverInternalSkipTestS1391;
  moonbit_string_t _M0L8_2afieldS2943;
  int32_t _M0L6_2acntS3212;
  moonbit_string_t _M0L7_2anameS1392;
  #line 515 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  moonbit_incref(_M0L8filenameS1362);
  _closure_3246
  = (struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360*)moonbit_malloc(sizeof(struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360));
  Moonbit_object_header(_closure_3246)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 0, 0);
  _closure_3246->code
  = &_M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1360;
  _closure_3246->$0 = _M0L5indexS1364;
  _closure_3246->$1 = _M0L8filenameS1362;
  _M0L13handle__startS1360 = (struct _M0TWEOc*)_closure_3246;
  moonbit_incref(_M0L8filenameS1362);
  _closure_3247
  = (struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365*)moonbit_malloc(sizeof(struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365));
  Moonbit_object_header(_closure_3247)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 3, 0);
  _closure_3247->code
  = &_M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1365;
  _closure_3247->$0 = _M0L5indexS1364;
  _closure_3247->$1 = _M0L8filenameS1362;
  _M0L14handle__resultS1365 = (struct _M0TWssbEu*)_closure_3247;
  _M0L17error__to__stringS1372
  = (struct _M0TWRPC15error5ErrorEs*)&_M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1372$closure.data;
  #line 557 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _tmp_3249
  = _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(_M0L12async__testsS1393, _M0L8filenameS1362, _M0L5indexS1364, _M0L13handle__startS1360, _M0L14handle__resultS1365, _M0L17error__to__stringS1372);
  if (_tmp_3249.tag) {
    int32_t const _M0L5_2aokS2940 = _tmp_3249.data.ok;
    _handle__error__result_3250 = _M0L5_2aokS2940;
  } else {
    void* const _M0L6_2aerrS2941 = _tmp_3249.data.err;
    moonbit_decref(_M0L17error__to__stringS1372);
    moonbit_decref(_M0L13handle__startS1360);
    _M0L11_2atry__errS1387 = _M0L6_2aerrS2941;
    goto join_1386;
  }
  if (_handle__error__result_3250) {
    moonbit_decref(_M0L17error__to__stringS1372);
    moonbit_decref(_M0L13handle__startS1360);
    _M0L6_2atmpS2931 = 1;
  } else {
    struct moonbit_result_0 _tmp_3251;
    int32_t _handle__error__result_3252;
    #line 560 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
    _tmp_3251
    = _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(_M0L12async__testsS1393, _M0L8filenameS1362, _M0L5indexS1364, _M0L13handle__startS1360, _M0L14handle__resultS1365, _M0L17error__to__stringS1372);
    if (_tmp_3251.tag) {
      int32_t const _M0L5_2aokS2938 = _tmp_3251.data.ok;
      _handle__error__result_3252 = _M0L5_2aokS2938;
    } else {
      void* const _M0L6_2aerrS2939 = _tmp_3251.data.err;
      moonbit_decref(_M0L17error__to__stringS1372);
      moonbit_decref(_M0L13handle__startS1360);
      _M0L11_2atry__errS1387 = _M0L6_2aerrS2939;
      goto join_1386;
    }
    if (_handle__error__result_3252) {
      moonbit_decref(_M0L17error__to__stringS1372);
      moonbit_decref(_M0L13handle__startS1360);
      _M0L6_2atmpS2931 = 1;
    } else {
      struct moonbit_result_0 _tmp_3253;
      int32_t _handle__error__result_3254;
      #line 563 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _tmp_3253
      = _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(_M0L12async__testsS1393, _M0L8filenameS1362, _M0L5indexS1364, _M0L13handle__startS1360, _M0L14handle__resultS1365, _M0L17error__to__stringS1372);
      if (_tmp_3253.tag) {
        int32_t const _M0L5_2aokS2936 = _tmp_3253.data.ok;
        _handle__error__result_3254 = _M0L5_2aokS2936;
      } else {
        void* const _M0L6_2aerrS2937 = _tmp_3253.data.err;
        moonbit_decref(_M0L17error__to__stringS1372);
        moonbit_decref(_M0L13handle__startS1360);
        _M0L11_2atry__errS1387 = _M0L6_2aerrS2937;
        goto join_1386;
      }
      if (_handle__error__result_3254) {
        moonbit_decref(_M0L17error__to__stringS1372);
        moonbit_decref(_M0L13handle__startS1360);
        _M0L6_2atmpS2931 = 1;
      } else {
        struct moonbit_result_0 _tmp_3255;
        int32_t _handle__error__result_3256;
        #line 566 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        _tmp_3255
        = _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(_M0L12async__testsS1393, _M0L8filenameS1362, _M0L5indexS1364, _M0L13handle__startS1360, _M0L14handle__resultS1365, _M0L17error__to__stringS1372);
        if (_tmp_3255.tag) {
          int32_t const _M0L5_2aokS2934 = _tmp_3255.data.ok;
          _handle__error__result_3256 = _M0L5_2aokS2934;
        } else {
          void* const _M0L6_2aerrS2935 = _tmp_3255.data.err;
          moonbit_decref(_M0L17error__to__stringS1372);
          moonbit_decref(_M0L13handle__startS1360);
          _M0L11_2atry__errS1387 = _M0L6_2aerrS2935;
          goto join_1386;
        }
        if (_handle__error__result_3256) {
          moonbit_decref(_M0L17error__to__stringS1372);
          moonbit_decref(_M0L13handle__startS1360);
          _M0L6_2atmpS2931 = 1;
        } else {
          struct moonbit_result_0 _tmp_3257;
          #line 569 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
          _tmp_3257
          = _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(_M0L12async__testsS1393, _M0L8filenameS1362, _M0L5indexS1364, _M0L13handle__startS1360, _M0L14handle__resultS1365, _M0L17error__to__stringS1372);
          moonbit_decref(_M0L13handle__startS1360);
          moonbit_decref(_M0L17error__to__stringS1372);
          if (_tmp_3257.tag) {
            int32_t const _M0L5_2aokS2932 = _tmp_3257.data.ok;
            _M0L6_2atmpS2931 = _M0L5_2aokS2932;
          } else {
            void* const _M0L6_2aerrS2933 = _tmp_3257.data.err;
            _M0L11_2atry__errS1387 = _M0L6_2aerrS2933;
            goto join_1386;
          }
        }
      }
    }
  }
  if (!_M0L6_2atmpS2931) {
    void* _M0L127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2942 =
      (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
    Moonbit_object_header(_M0L127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2942)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 2);
    ((struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2942)->$0
    = (moonbit_string_t)moonbit_string_literal_0.data;
    _M0L11_2atry__errS1387
    = _M0L127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS2942;
    goto join_1386;
  } else {
    moonbit_decref(_M0L14handle__resultS1365);
  }
  goto joinlet_3248;
  join_1386:;
  _M0L3errS1388 = _M0L11_2atry__errS1387;
  _M0L36_2aMoonBitTestDriverInternalSkipTestS1391
  = (struct _M0DTPC15error5Error127biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L3errS1388;
  _M0L8_2afieldS2943 = _M0L36_2aMoonBitTestDriverInternalSkipTestS1391->$0;
  _M0L6_2acntS3212
  = Moonbit_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1391));
  if (_M0L6_2acntS3212 > 1) {
    int32_t _M0L11_2anew__cntS3213 = _M0L6_2acntS3212 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS1391), _M0L11_2anew__cntS3213);
    moonbit_incref(_M0L8_2afieldS2943);
  } else if (_M0L6_2acntS3212 == 1) {
    #line 576 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
    moonbit_free(_M0L36_2aMoonBitTestDriverInternalSkipTestS1391);
  }
  _M0L7_2anameS1392 = _M0L8_2afieldS2943;
  _M0L4nameS1390 = _M0L7_2anameS1392;
  goto join_1389;
  goto joinlet_3258;
  join_1389:;
  #line 577 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1365(_M0L14handle__resultS1365, _M0L4nameS1390, (moonbit_string_t)moonbit_string_literal_1.data, 1);
  moonbit_decref(_M0L14handle__resultS1365);
  moonbit_decref(_M0L4nameS1390);
  joinlet_3258:;
  joinlet_3248:;
  return 0;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN17error__to__stringS1372(
  struct _M0TWRPC15error5ErrorEs* _M0L6_2aenvS2930,
  void* _M0L3errS1373
) {
  void* _M0L1eS1375;
  moonbit_string_t _M0L1eS1377;
  moonbit_string_t _result_3261;
  #line 546 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  switch (Moonbit_object_tag(_M0L3errS1373)) {
    case 0: {
      struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS1378 =
        (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L3errS1373;
      moonbit_string_t _M0L4_2aeS1379 = _M0L10_2aFailureS1378->$0;
      moonbit_incref(_M0L4_2aeS1379);
      _M0L1eS1377 = _M0L4_2aeS1379;
      goto join_1376;
      break;
    }
    
    case 3: {
      struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError* _M0L15_2aInspectErrorS1380 =
        (struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError*)_M0L3errS1373;
      moonbit_string_t _M0L4_2aeS1381 = _M0L15_2aInspectErrorS1380->$0;
      moonbit_incref(_M0L4_2aeS1381);
      _M0L1eS1377 = _M0L4_2aeS1381;
      goto join_1376;
      break;
    }
    
    case 4: {
      struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError* _M0L16_2aSnapshotErrorS1382 =
        (struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError*)_M0L3errS1373;
      moonbit_string_t _M0L4_2aeS1383 = _M0L16_2aSnapshotErrorS1382->$0;
      moonbit_incref(_M0L4_2aeS1383);
      _M0L1eS1377 = _M0L4_2aeS1383;
      goto join_1376;
      break;
    }
    
    case 5: {
      struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError* _M0L35_2aMoonBitTestDriverInternalJsErrorS1384 =
        (struct _M0DTPC15error5Error125biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError*)_M0L3errS1373;
      moonbit_string_t _M0L4_2aeS1385 =
        _M0L35_2aMoonBitTestDriverInternalJsErrorS1384->$0;
      moonbit_incref(_M0L4_2aeS1385);
      _M0L1eS1377 = _M0L4_2aeS1385;
      goto join_1376;
      break;
    }
    default: {
      moonbit_incref(_M0L3errS1373);
      _M0L1eS1375 = _M0L3errS1373;
      goto join_1374;
      break;
    }
  }
  join_1376:;
  return _M0L1eS1377;
  join_1374:;
  #line 552 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3261 = _M0FP15Error10to__string(_M0L1eS1375);
  moonbit_decref(_M0L1eS1375);
  return _result_3261;
}

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN14handle__resultS1365(
  struct _M0TWssbEu* _M0L6_2aenvS2927,
  moonbit_string_t _M0L10__testnameS1366,
  moonbit_string_t _M0L7messageS1367,
  int32_t _M0L7skippedS1368
) {
  struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365* _M0L14_2acasted__envS2928;
  moonbit_string_t _M0L8filenameS1362;
  int32_t _M0L5indexS1364;
  int32_t _if__result_3262;
  moonbit_string_t _M0L10file__nameS1369;
  moonbit_string_t _M0L7messageS1370;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1371;
  moonbit_string_t _M0L6_2atmpS2929;
  #line 531 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L14_2acasted__envS2928
  = (struct _M0R129_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c1365*)_M0L6_2aenvS2927;
  _M0L8filenameS1362 = _M0L14_2acasted__envS2928->$1;
  _M0L5indexS1364 = _M0L14_2acasted__envS2928->$0;
  if (!_M0L7skippedS1368) {
    _if__result_3262 = 1;
  } else {
    _if__result_3262 = 0;
  }
  if (_if__result_3262) {
    
  }
  #line 537 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L10file__nameS1369
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1362, 1);
  #line 538 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L7messageS1370
  = _M0MPC16string6String14escape_2einner(_M0L7messageS1367, 1);
  #line 539 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L18_2astring__builderS1371
  = _M0MPB13StringBuilder21StringBuilder_2einner(45);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1371, (moonbit_string_t)moonbit_string_literal_3.data);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1371, _M0L10file__nameS1369);
  moonbit_decref(_M0L10file__nameS1369);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1371, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1371, _M0L5indexS1364);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1371, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1371, _M0L7messageS1370);
  moonbit_decref(_M0L7messageS1370);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1371, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 541 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2929
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1371);
  moonbit_decref(_M0L18_2astring__builderS1371);
  #line 540 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS2929);
  moonbit_decref(_M0L6_2atmpS2929);
  #line 543 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__executeN13handle__startS1360(
  struct _M0TWEOc* _M0L6_2aenvS2924
) {
  struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360* _M0L14_2acasted__envS2925;
  moonbit_string_t _M0L8filenameS1362;
  int32_t _M0L5indexS1364;
  moonbit_string_t _M0L10file__nameS1361;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS1363;
  moonbit_string_t _M0L6_2atmpS2926;
  #line 522 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L14_2acasted__envS2925
  = (struct _M0R128_24biolab_2fbio__seq_2fcmd_2falignio__main__blackbox__test_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c1360*)_M0L6_2aenvS2924;
  _M0L8filenameS1362 = _M0L14_2acasted__envS2925->$1;
  _M0L5indexS1364 = _M0L14_2acasted__envS2925->$0;
  #line 523 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L10file__nameS1361
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS1362, 1);
  #line 524 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L18_2astring__builderS1363
  = _M0MPB13StringBuilder21StringBuilder_2einner(33);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1363, (moonbit_string_t)moonbit_string_literal_8.data);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS1363, _M0L10file__nameS1361);
  moonbit_decref(_M0L10file__nameS1361);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1363, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS1363, _M0L5indexS1364);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS1363, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 526 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2926
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS1363);
  moonbit_decref(_M0L18_2astring__builderS1363);
  #line 525 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS2926);
  moonbit_decref(_M0L6_2atmpS2926);
  #line 528 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__args(
  
) {
  int32_t _M0L45moonbit__test__driver__internal__parse__int__S1317;
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324;
  int32_t _M0L57moonbit__test__driver__internal__get__cli__args__internalS1331;
  int32_t _M0L51moonbit__test__driver__internal__split__mbt__stringS1338;
  struct _M0TUsiE** _M0L6_2atmpS2923;
  struct _M0TPB5ArrayGUsiEE* _M0L16file__and__indexS1345;
  struct _M0TPB5ArrayGsE* _M0L9cli__argsS1346;
  moonbit_string_t _M0L6_2atmpS2922;
  struct _M0TPB5ArrayGsE* _M0L10test__argsS1347;
  int32_t _M0L7_2abindS1348;
  int32_t _M0L2__S1349;
  #line 194 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L45moonbit__test__driver__internal__parse__int__S1317 = 0;
  _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324
  = 0;
  _M0L57moonbit__test__driver__internal__get__cli__args__internalS1331
  = _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324;
  _M0L51moonbit__test__driver__internal__split__mbt__stringS1338 = 0;
  _M0L6_2atmpS2923 = (struct _M0TUsiE**)moonbit_empty_ref_array;
  _M0L16file__and__indexS1345
  = (struct _M0TPB5ArrayGUsiEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGUsiEE));
  Moonbit_object_header(_M0L16file__and__indexS1345)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 9, 0);
  _M0L16file__and__indexS1345->$0 = _M0L6_2atmpS2923;
  _M0L16file__and__indexS1345->$1 = 0;
  #line 283 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L9cli__argsS1346
  = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1331(_M0L57moonbit__test__driver__internal__get__cli__args__internalS1331);
  #line 285 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2922 = _M0MPC15array5Array2atGsE(_M0L9cli__argsS1346, 1);
  moonbit_decref(_M0L9cli__argsS1346);
  #line 284 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L10test__argsS1347
  = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1338(_M0L51moonbit__test__driver__internal__split__mbt__stringS1338, _M0L6_2atmpS2922, 47);
  moonbit_decref(_M0L6_2atmpS2922);
  _M0L7_2abindS1348 = _M0L10test__argsS1347->$1;
  _M0L2__S1349 = 0;
  while (1) {
    if (_M0L2__S1349 < _M0L7_2abindS1348) {
      moonbit_string_t* _M0L3bufS2921 = _M0L10test__argsS1347->$0;
      moonbit_string_t _M0L3argS1350 =
        (moonbit_string_t)_M0L3bufS2921[_M0L2__S1349];
      struct _M0TPB5ArrayGsE* _M0L16file__and__rangeS1351;
      moonbit_string_t _M0L4fileS1352;
      moonbit_string_t _M0L5rangeS1353;
      struct _M0TPB5ArrayGsE* _M0L15start__and__endS1354;
      moonbit_string_t _M0L6_2atmpS2919;
      int32_t _M0L5startS1355;
      moonbit_string_t _M0L6_2atmpS2918;
      int32_t _M0L3endS1356;
      int32_t _M0L1iS1357;
      int32_t _M0L6_2atmpS2920;
      moonbit_incref(_M0L3argS1350);
      #line 289 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L16file__and__rangeS1351
      = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1338(_M0L51moonbit__test__driver__internal__split__mbt__stringS1338, _M0L3argS1350, 58);
      moonbit_decref(_M0L3argS1350);
      #line 290 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L4fileS1352
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1351, 0);
      #line 291 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L5rangeS1353
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS1351, 1);
      moonbit_decref(_M0L16file__and__rangeS1351);
      #line 292 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L15start__and__endS1354
      = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1338(_M0L51moonbit__test__driver__internal__split__mbt__stringS1338, _M0L5rangeS1353, 45);
      moonbit_decref(_M0L5rangeS1353);
      #line 295 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS2919
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1354, 0);
      #line 295 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L5startS1355
      = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1317(_M0L45moonbit__test__driver__internal__parse__int__S1317, _M0L6_2atmpS2919);
      moonbit_decref(_M0L6_2atmpS2919);
      #line 296 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS2918
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS1354, 1);
      moonbit_decref(_M0L15start__and__endS1354);
      #line 296 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L3endS1356
      = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1317(_M0L45moonbit__test__driver__internal__parse__int__S1317, _M0L6_2atmpS2918);
      moonbit_decref(_M0L6_2atmpS2918);
      _M0L1iS1357 = _M0L5startS1355;
      while (1) {
        if (_M0L1iS1357 < _M0L3endS1356) {
          struct _M0TUsiE* _M0L8_2atupleS2916;
          int32_t _M0L6_2atmpS2917;
          moonbit_incref(_M0L4fileS1352);
          _M0L8_2atupleS2916
          = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
          Moonbit_object_header(_M0L8_2atupleS2916)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 12, 0);
          _M0L8_2atupleS2916->$0 = _M0L4fileS1352;
          _M0L8_2atupleS2916->$1 = _M0L1iS1357;
          #line 298 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
          _M0MPC15array5Array4pushGUsiEE(_M0L16file__and__indexS1345, _M0L8_2atupleS2916);
          moonbit_decref(_M0L8_2atupleS2916);
          _M0L6_2atmpS2917 = _M0L1iS1357 + 1;
          _M0L1iS1357 = _M0L6_2atmpS2917;
          continue;
        } else {
          moonbit_decref(_M0L4fileS1352);
        }
        break;
      }
      _M0L6_2atmpS2920 = _M0L2__S1349 + 1;
      _M0L2__S1349 = _M0L6_2atmpS2920;
      continue;
    } else {
      moonbit_decref(_M0L10test__argsS1347);
    }
    break;
  }
  return _M0L16file__and__indexS1345;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS1338(
  int32_t _M0L6_2aenvS2897,
  moonbit_string_t _M0L1sS1339,
  int32_t _M0L3sepS1340
) {
  moonbit_string_t* _M0L6_2atmpS2915;
  struct _M0TPB5ArrayGsE* _M0L3resS1341;
  struct _M0TPB8MutLocalGiE* _M0L1iS1342;
  struct _M0TPB8MutLocalGiE* _M0L5startS1343;
  int32_t _M0L3valS2910;
  int32_t _M0L6_2atmpS2911;
  #line 262 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L6_2atmpS2915 = (moonbit_string_t*)moonbit_empty_ref_array;
  _M0L3resS1341
  = (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
  Moonbit_object_header(_M0L3resS1341)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
  _M0L3resS1341->$0 = _M0L6_2atmpS2915;
  _M0L3resS1341->$1 = 0;
  _M0L1iS1342
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1342)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1342->$0 = 0;
  _M0L5startS1343
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS1343)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS1343->$0 = 0;
  while (1) {
    int32_t _M0L3valS2898 = _M0L1iS1342->$0;
    int32_t _M0L6_2atmpS2899 = Moonbit_array_length(_M0L1sS1339);
    if (_M0L3valS2898 < _M0L6_2atmpS2899) {
      int32_t _M0L3valS2902 = _M0L1iS1342->$0;
      int32_t _M0L6_2atmpS2901;
      int32_t _M0L6_2atmpS2900;
      int32_t _M0L3valS2909;
      int32_t _M0L6_2atmpS2908;
      if (
        _M0L3valS2902 < 0
        || _M0L3valS2902 >= Moonbit_array_length(_M0L1sS1339)
      ) {
        #line 270 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2901 = _M0L1sS1339[_M0L3valS2902];
      _M0L6_2atmpS2900 = _M0L6_2atmpS2901;
      if (_M0L6_2atmpS2900 == _M0L3sepS1340) {
        int32_t _M0L3valS2904 = _M0L5startS1343->$0;
        int32_t _M0L3valS2905 = _M0L1iS1342->$0;
        moonbit_string_t _M0L6_2atmpS2903;
        int32_t _M0L3valS2907;
        int32_t _M0L6_2atmpS2906;
        #line 271 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        _M0L6_2atmpS2903
        = _M0MPC16string6String17unsafe__substring(_M0L1sS1339, _M0L3valS2904, _M0L3valS2905);
        #line 271 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        _M0MPC15array5Array4pushGsE(_M0L3resS1341, _M0L6_2atmpS2903);
        moonbit_decref(_M0L6_2atmpS2903);
        _M0L3valS2907 = _M0L1iS1342->$0;
        _M0L6_2atmpS2906 = _M0L3valS2907 + 1;
        _M0L5startS1343->$0 = _M0L6_2atmpS2906;
      }
      _M0L3valS2909 = _M0L1iS1342->$0;
      _M0L6_2atmpS2908 = _M0L3valS2909 + 1;
      _M0L1iS1342->$0 = _M0L6_2atmpS2908;
      continue;
    } else {
      moonbit_decref(_M0L1iS1342);
    }
    break;
  }
  _M0L3valS2910 = _M0L5startS1343->$0;
  _M0L6_2atmpS2911 = Moonbit_array_length(_M0L1sS1339);
  if (_M0L3valS2910 < _M0L6_2atmpS2911) {
    int32_t _M0L3valS2913 = _M0L5startS1343->$0;
    int32_t _M0L6_2atmpS2914;
    moonbit_string_t _M0L6_2atmpS2912;
    moonbit_decref(_M0L5startS1343);
    _M0L6_2atmpS2914 = Moonbit_array_length(_M0L1sS1339);
    #line 277 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
    _M0L6_2atmpS2912
    = _M0MPC16string6String17unsafe__substring(_M0L1sS1339, _M0L3valS2913, _M0L6_2atmpS2914);
    #line 277 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
    _M0MPC15array5Array4pushGsE(_M0L3resS1341, _M0L6_2atmpS2912);
    moonbit_decref(_M0L6_2atmpS2912);
  } else {
    moonbit_decref(_M0L5startS1343);
  }
  return _M0L3resS1341;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS1331(
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324
) {
  moonbit_bytes_t* _M0L3tmpS1332;
  int32_t _M0L6_2atmpS2896;
  struct _M0TPB5ArrayGsE* _M0L3resS1333;
  int32_t _M0L7_2abindS1334;
  int32_t _M0L7_2abindS1335;
  int32_t _M0L1iS1336;
  #line 251 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  #line 254 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3tmpS1332
  = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__get__cli__args__ffi();
  _M0L6_2atmpS2896 = Moonbit_array_length(_M0L3tmpS1332);
  #line 255 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1333 = _M0MPC15array5Array11new_2einnerGsE(_M0L6_2atmpS2896);
  _M0L7_2abindS1334 = 0;
  _M0L7_2abindS1335 = Moonbit_array_length(_M0L3tmpS1332);
  _M0L1iS1336 = _M0L7_2abindS1334;
  while (1) {
    if (_M0L1iS1336 < _M0L7_2abindS1335) {
      moonbit_bytes_t _M0L6_2atmpS2894;
      moonbit_string_t _M0L6_2atmpS2893;
      int32_t _M0L6_2atmpS2895;
      if (
        _M0L1iS1336 < 0 || _M0L1iS1336 >= Moonbit_array_length(_M0L3tmpS1332)
      ) {
        #line 257 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2894 = (moonbit_bytes_t)_M0L3tmpS1332[_M0L1iS1336];
      moonbit_incref(_M0L6_2atmpS2894);
      #line 257 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0L6_2atmpS2893
      = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324(_M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324, _M0L6_2atmpS2894);
      moonbit_decref(_M0L6_2atmpS2894);
      #line 257 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0MPC15array5Array4pushGsE(_M0L3resS1333, _M0L6_2atmpS2893);
      moonbit_decref(_M0L6_2atmpS2893);
      _M0L6_2atmpS2895 = _M0L1iS1336 + 1;
      _M0L1iS1336 = _M0L6_2atmpS2895;
      continue;
    } else {
      moonbit_decref(_M0L3tmpS1332);
    }
    break;
  }
  return _M0L3resS1333;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS1324(
  int32_t _M0L6_2aenvS2807,
  moonbit_bytes_t _M0L5bytesS1325
) {
  struct _M0TPB13StringBuilder* _M0L3resS1326;
  int32_t _M0L3lenS1327;
  struct _M0TPB8MutLocalGiE* _M0L1iS1328;
  moonbit_string_t _result_3268;
  #line 207 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  #line 210 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1326 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3lenS1327 = Moonbit_array_length(_M0L5bytesS1325);
  _M0L1iS1328
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1328)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1328->$0 = 0;
  while (1) {
    int32_t _M0L3valS2808 = _M0L1iS1328->$0;
    if (_M0L3valS2808 < _M0L3lenS1327) {
      int32_t _M0L3valS2892 = _M0L1iS1328->$0;
      int32_t _M0L6_2atmpS2891;
      int32_t _M0L6_2atmpS2890;
      struct _M0TPB8MutLocalGiE* _M0L1cS1329;
      int32_t _M0L3valS2809;
      if (
        _M0L3valS2892 < 0
        || _M0L3valS2892 >= Moonbit_array_length(_M0L5bytesS1325)
      ) {
        #line 214 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2891 = _M0L5bytesS1325[_M0L3valS2892];
      _M0L6_2atmpS2890 = (int32_t)_M0L6_2atmpS2891;
      _M0L1cS1329
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1cS1329)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1cS1329->$0 = _M0L6_2atmpS2890;
      _M0L3valS2809 = _M0L1cS1329->$0;
      if (_M0L3valS2809 < 128) {
        int32_t _M0L3valS2811 = _M0L1cS1329->$0;
        int32_t _M0L6_2atmpS2810;
        int32_t _M0L3valS2813;
        int32_t _M0L6_2atmpS2812;
        moonbit_decref(_M0L1cS1329);
        _M0L6_2atmpS2810 = _M0L3valS2811;
        #line 216 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1326, _M0L6_2atmpS2810);
        _M0L3valS2813 = _M0L1iS1328->$0;
        _M0L6_2atmpS2812 = _M0L3valS2813 + 1;
        _M0L1iS1328->$0 = _M0L6_2atmpS2812;
      } else {
        int32_t _M0L3valS2814 = _M0L1cS1329->$0;
        if (_M0L3valS2814 < 224) {
          int32_t _M0L3valS2816 = _M0L1iS1328->$0;
          int32_t _M0L6_2atmpS2815 = _M0L3valS2816 + 1;
          int32_t _M0L3valS2825;
          int32_t _M0L6_2atmpS2824;
          int32_t _M0L6_2atmpS2818;
          int32_t _M0L3valS2823;
          int32_t _M0L6_2atmpS2822;
          int32_t _M0L6_2atmpS2821;
          int32_t _M0L6_2atmpS2820;
          int32_t _M0L6_2atmpS2819;
          int32_t _M0L6_2atmpS2817;
          int32_t _M0L3valS2827;
          int32_t _M0L6_2atmpS2826;
          int32_t _M0L3valS2829;
          int32_t _M0L6_2atmpS2828;
          if (_M0L6_2atmpS2815 >= _M0L3lenS1327) {
            moonbit_decref(_M0L1cS1329);
            moonbit_decref(_M0L1iS1328);
            break;
          }
          _M0L3valS2825 = _M0L1cS1329->$0;
          _M0L6_2atmpS2824 = _M0L3valS2825 & 31;
          _M0L6_2atmpS2818 = _M0L6_2atmpS2824 << 6;
          _M0L3valS2823 = _M0L1iS1328->$0;
          _M0L6_2atmpS2822 = _M0L3valS2823 + 1;
          if (
            _M0L6_2atmpS2822 < 0
            || _M0L6_2atmpS2822 >= Moonbit_array_length(_M0L5bytesS1325)
          ) {
            #line 222 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS2821 = _M0L5bytesS1325[_M0L6_2atmpS2822];
          _M0L6_2atmpS2820 = (int32_t)_M0L6_2atmpS2821;
          _M0L6_2atmpS2819 = _M0L6_2atmpS2820 & 63;
          _M0L6_2atmpS2817 = _M0L6_2atmpS2818 | _M0L6_2atmpS2819;
          _M0L1cS1329->$0 = _M0L6_2atmpS2817;
          _M0L3valS2827 = _M0L1cS1329->$0;
          moonbit_decref(_M0L1cS1329);
          _M0L6_2atmpS2826 = _M0L3valS2827;
          #line 223 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1326, _M0L6_2atmpS2826);
          _M0L3valS2829 = _M0L1iS1328->$0;
          _M0L6_2atmpS2828 = _M0L3valS2829 + 2;
          _M0L1iS1328->$0 = _M0L6_2atmpS2828;
        } else {
          int32_t _M0L3valS2830 = _M0L1cS1329->$0;
          if (_M0L3valS2830 < 240) {
            int32_t _M0L3valS2832 = _M0L1iS1328->$0;
            int32_t _M0L6_2atmpS2831 = _M0L3valS2832 + 2;
            int32_t _M0L3valS2848;
            int32_t _M0L6_2atmpS2847;
            int32_t _M0L6_2atmpS2840;
            int32_t _M0L3valS2846;
            int32_t _M0L6_2atmpS2845;
            int32_t _M0L6_2atmpS2844;
            int32_t _M0L6_2atmpS2843;
            int32_t _M0L6_2atmpS2842;
            int32_t _M0L6_2atmpS2841;
            int32_t _M0L6_2atmpS2834;
            int32_t _M0L3valS2839;
            int32_t _M0L6_2atmpS2838;
            int32_t _M0L6_2atmpS2837;
            int32_t _M0L6_2atmpS2836;
            int32_t _M0L6_2atmpS2835;
            int32_t _M0L6_2atmpS2833;
            int32_t _M0L3valS2850;
            int32_t _M0L6_2atmpS2849;
            int32_t _M0L3valS2852;
            int32_t _M0L6_2atmpS2851;
            if (_M0L6_2atmpS2831 >= _M0L3lenS1327) {
              moonbit_decref(_M0L1cS1329);
              moonbit_decref(_M0L1iS1328);
              break;
            }
            _M0L3valS2848 = _M0L1cS1329->$0;
            _M0L6_2atmpS2847 = _M0L3valS2848 & 15;
            _M0L6_2atmpS2840 = _M0L6_2atmpS2847 << 12;
            _M0L3valS2846 = _M0L1iS1328->$0;
            _M0L6_2atmpS2845 = _M0L3valS2846 + 1;
            if (
              _M0L6_2atmpS2845 < 0
              || _M0L6_2atmpS2845 >= Moonbit_array_length(_M0L5bytesS1325)
            ) {
              #line 230 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2844 = _M0L5bytesS1325[_M0L6_2atmpS2845];
            _M0L6_2atmpS2843 = (int32_t)_M0L6_2atmpS2844;
            _M0L6_2atmpS2842 = _M0L6_2atmpS2843 & 63;
            _M0L6_2atmpS2841 = _M0L6_2atmpS2842 << 6;
            _M0L6_2atmpS2834 = _M0L6_2atmpS2840 | _M0L6_2atmpS2841;
            _M0L3valS2839 = _M0L1iS1328->$0;
            _M0L6_2atmpS2838 = _M0L3valS2839 + 2;
            if (
              _M0L6_2atmpS2838 < 0
              || _M0L6_2atmpS2838 >= Moonbit_array_length(_M0L5bytesS1325)
            ) {
              #line 231 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2837 = _M0L5bytesS1325[_M0L6_2atmpS2838];
            _M0L6_2atmpS2836 = (int32_t)_M0L6_2atmpS2837;
            _M0L6_2atmpS2835 = _M0L6_2atmpS2836 & 63;
            _M0L6_2atmpS2833 = _M0L6_2atmpS2834 | _M0L6_2atmpS2835;
            _M0L1cS1329->$0 = _M0L6_2atmpS2833;
            _M0L3valS2850 = _M0L1cS1329->$0;
            moonbit_decref(_M0L1cS1329);
            _M0L6_2atmpS2849 = _M0L3valS2850;
            #line 232 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1326, _M0L6_2atmpS2849);
            _M0L3valS2852 = _M0L1iS1328->$0;
            _M0L6_2atmpS2851 = _M0L3valS2852 + 3;
            _M0L1iS1328->$0 = _M0L6_2atmpS2851;
          } else {
            int32_t _M0L3valS2854 = _M0L1iS1328->$0;
            int32_t _M0L6_2atmpS2853 = _M0L3valS2854 + 3;
            int32_t _M0L3valS2877;
            int32_t _M0L6_2atmpS2876;
            int32_t _M0L6_2atmpS2869;
            int32_t _M0L3valS2875;
            int32_t _M0L6_2atmpS2874;
            int32_t _M0L6_2atmpS2873;
            int32_t _M0L6_2atmpS2872;
            int32_t _M0L6_2atmpS2871;
            int32_t _M0L6_2atmpS2870;
            int32_t _M0L6_2atmpS2862;
            int32_t _M0L3valS2868;
            int32_t _M0L6_2atmpS2867;
            int32_t _M0L6_2atmpS2866;
            int32_t _M0L6_2atmpS2865;
            int32_t _M0L6_2atmpS2864;
            int32_t _M0L6_2atmpS2863;
            int32_t _M0L6_2atmpS2856;
            int32_t _M0L3valS2861;
            int32_t _M0L6_2atmpS2860;
            int32_t _M0L6_2atmpS2859;
            int32_t _M0L6_2atmpS2858;
            int32_t _M0L6_2atmpS2857;
            int32_t _M0L6_2atmpS2855;
            int32_t _M0L3valS2879;
            int32_t _M0L6_2atmpS2878;
            int32_t _M0L3valS2883;
            int32_t _M0L6_2atmpS2882;
            int32_t _M0L6_2atmpS2881;
            int32_t _M0L6_2atmpS2880;
            int32_t _M0L3valS2887;
            int32_t _M0L6_2atmpS2886;
            int32_t _M0L6_2atmpS2885;
            int32_t _M0L6_2atmpS2884;
            int32_t _M0L3valS2889;
            int32_t _M0L6_2atmpS2888;
            if (_M0L6_2atmpS2853 >= _M0L3lenS1327) {
              moonbit_decref(_M0L1cS1329);
              moonbit_decref(_M0L1iS1328);
              break;
            }
            _M0L3valS2877 = _M0L1cS1329->$0;
            _M0L6_2atmpS2876 = _M0L3valS2877 & 7;
            _M0L6_2atmpS2869 = _M0L6_2atmpS2876 << 18;
            _M0L3valS2875 = _M0L1iS1328->$0;
            _M0L6_2atmpS2874 = _M0L3valS2875 + 1;
            if (
              _M0L6_2atmpS2874 < 0
              || _M0L6_2atmpS2874 >= Moonbit_array_length(_M0L5bytesS1325)
            ) {
              #line 239 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2873 = _M0L5bytesS1325[_M0L6_2atmpS2874];
            _M0L6_2atmpS2872 = (int32_t)_M0L6_2atmpS2873;
            _M0L6_2atmpS2871 = _M0L6_2atmpS2872 & 63;
            _M0L6_2atmpS2870 = _M0L6_2atmpS2871 << 12;
            _M0L6_2atmpS2862 = _M0L6_2atmpS2869 | _M0L6_2atmpS2870;
            _M0L3valS2868 = _M0L1iS1328->$0;
            _M0L6_2atmpS2867 = _M0L3valS2868 + 2;
            if (
              _M0L6_2atmpS2867 < 0
              || _M0L6_2atmpS2867 >= Moonbit_array_length(_M0L5bytesS1325)
            ) {
              #line 240 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2866 = _M0L5bytesS1325[_M0L6_2atmpS2867];
            _M0L6_2atmpS2865 = (int32_t)_M0L6_2atmpS2866;
            _M0L6_2atmpS2864 = _M0L6_2atmpS2865 & 63;
            _M0L6_2atmpS2863 = _M0L6_2atmpS2864 << 6;
            _M0L6_2atmpS2856 = _M0L6_2atmpS2862 | _M0L6_2atmpS2863;
            _M0L3valS2861 = _M0L1iS1328->$0;
            _M0L6_2atmpS2860 = _M0L3valS2861 + 3;
            if (
              _M0L6_2atmpS2860 < 0
              || _M0L6_2atmpS2860 >= Moonbit_array_length(_M0L5bytesS1325)
            ) {
              #line 241 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS2859 = _M0L5bytesS1325[_M0L6_2atmpS2860];
            _M0L6_2atmpS2858 = (int32_t)_M0L6_2atmpS2859;
            _M0L6_2atmpS2857 = _M0L6_2atmpS2858 & 63;
            _M0L6_2atmpS2855 = _M0L6_2atmpS2856 | _M0L6_2atmpS2857;
            _M0L1cS1329->$0 = _M0L6_2atmpS2855;
            _M0L3valS2879 = _M0L1cS1329->$0;
            _M0L6_2atmpS2878 = _M0L3valS2879 - 65536;
            _M0L1cS1329->$0 = _M0L6_2atmpS2878;
            _M0L3valS2883 = _M0L1cS1329->$0;
            _M0L6_2atmpS2882 = _M0L3valS2883 >> 10;
            _M0L6_2atmpS2881 = _M0L6_2atmpS2882 + 55296;
            _M0L6_2atmpS2880 = _M0L6_2atmpS2881;
            #line 243 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1326, _M0L6_2atmpS2880);
            _M0L3valS2887 = _M0L1cS1329->$0;
            moonbit_decref(_M0L1cS1329);
            _M0L6_2atmpS2886 = _M0L3valS2887 & 1023;
            _M0L6_2atmpS2885 = _M0L6_2atmpS2886 + 56320;
            _M0L6_2atmpS2884 = _M0L6_2atmpS2885;
            #line 244 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS1326, _M0L6_2atmpS2884);
            _M0L3valS2889 = _M0L1iS1328->$0;
            _M0L6_2atmpS2888 = _M0L3valS2889 + 4;
            _M0L1iS1328->$0 = _M0L6_2atmpS2888;
          }
        }
      }
      continue;
    } else {
      moonbit_decref(_M0L1iS1328);
    }
    break;
  }
  #line 248 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3268 = _M0MPB13StringBuilder10to__string(_M0L3resS1326);
  moonbit_decref(_M0L3resS1326);
  return _result_3268;
}

int32_t _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S1317(
  int32_t _M0L6_2aenvS2800,
  moonbit_string_t _M0L1sS1318
) {
  struct _M0TPB8MutLocalGiE* _M0L3resS1319;
  int32_t _M0L3lenS1320;
  int32_t _M0L7_2abindS1321;
  int32_t _M0L1iS1322;
  int32_t _result_3270;
  #line 198 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L3resS1319
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3resS1319)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3resS1319->$0 = 0;
  _M0L3lenS1320 = Moonbit_array_length(_M0L1sS1318);
  _M0L7_2abindS1321 = 0;
  _M0L1iS1322 = _M0L7_2abindS1321;
  while (1) {
    if (_M0L1iS1322 < _M0L3lenS1320) {
      int32_t _M0L3valS2805 = _M0L3resS1319->$0;
      int32_t _M0L6_2atmpS2802 = _M0L3valS2805 * 10;
      int32_t _M0L6_2atmpS2804;
      int32_t _M0L6_2atmpS2803;
      int32_t _M0L6_2atmpS2801;
      int32_t _M0L6_2atmpS2806;
      if (
        _M0L1iS1322 < 0 || _M0L1iS1322 >= Moonbit_array_length(_M0L1sS1318)
      ) {
        #line 202 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2804 = _M0L1sS1318[_M0L1iS1322];
      _M0L6_2atmpS2803 = _M0L6_2atmpS2804 - 48;
      _M0L6_2atmpS2801 = _M0L6_2atmpS2802 + _M0L6_2atmpS2803;
      _M0L3resS1319->$0 = _M0L6_2atmpS2801;
      _M0L6_2atmpS2806 = _M0L1iS1322 + 1;
      _M0L1iS1322 = _M0L6_2atmpS2806;
      continue;
    }
    break;
  }
  _result_3270 = _M0L3resS1319->$0;
  moonbit_decref(_M0L3resS1319);
  return _result_3270;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1287,
  moonbit_string_t _M0L12_2adiscard__S1288,
  int32_t _M0L12_2adiscard__S1289,
  struct _M0TWEOc* _M0L12_2adiscard__S1290,
  struct _M0TWssbEu* _M0L12_2adiscard__S1291,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1292
) {
  struct moonbit_result_0 _result_3271;
  #line 35 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3271.tag = 1;
  _result_3271.data.ok = 0;
  return _result_3271;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1293,
  moonbit_string_t _M0L12_2adiscard__S1294,
  int32_t _M0L12_2adiscard__S1295,
  struct _M0TWEOc* _M0L12_2adiscard__S1296,
  struct _M0TWssbEu* _M0L12_2adiscard__S1297,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1298
) {
  struct moonbit_result_0 _result_3272;
  #line 35 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3272.tag = 1;
  _result_3272.data.ok = 0;
  return _result_3272;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1299,
  moonbit_string_t _M0L12_2adiscard__S1300,
  int32_t _M0L12_2adiscard__S1301,
  struct _M0TWEOc* _M0L12_2adiscard__S1302,
  struct _M0TWssbEu* _M0L12_2adiscard__S1303,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1304
) {
  struct moonbit_result_0 _result_3273;
  #line 35 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3273.tag = 1;
  _result_3273.data.ok = 0;
  return _result_3273;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1305,
  moonbit_string_t _M0L12_2adiscard__S1306,
  int32_t _M0L12_2adiscard__S1307,
  struct _M0TWEOc* _M0L12_2adiscard__S1308,
  struct _M0TWssbEu* _M0L12_2adiscard__S1309,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1310
) {
  struct moonbit_result_0 _result_3274;
  #line 35 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3274.tag = 1;
  _result_3274.data.ok = 0;
  return _result_3274;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1311,
  moonbit_string_t _M0L12_2adiscard__S1312,
  int32_t _M0L12_2adiscard__S1313,
  struct _M0TWEOc* _M0L12_2adiscard__S1314,
  struct _M0TWssbEu* _M0L12_2adiscard__S1315,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S1316
) {
  struct moonbit_result_0 _result_3275;
  #line 35 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _result_3275.tag = 1;
  _result_3275.data.ok = 0;
  return _result_3275;
}

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S1286
) {
  #line 12 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd13alignio__main13alignio__test() {
  void* _M0L11_2atry__errS1285;
  struct moonbit_result_0 _tmp_3277;
  #line 199 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  #line 201 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3277 = _M0FP46biolab8bio__seq3cmd13alignio__main18run__alignio__test();
  if (_tmp_3277.tag) {
    int32_t const _M0L5_2aokS2798 = _tmp_3277.data.ok;
  } else {
    void* const _M0L6_2aerrS2799 = _tmp_3277.data.err;
    _M0L11_2atry__errS1285 = _M0L6_2aerrS2799;
    goto join_1284;
  }
  goto joinlet_3276;
  join_1284:;
  moonbit_decref(_M0L11_2atry__errS1285);
  #line 200 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_9.data);
  joinlet_3276:;
  return 0;
}

struct moonbit_result_0 _M0FP46biolab8bio__seq3cmd13alignio__main18run__alignio__test(
  
) {
  moonbit_string_t _M0L15clustal__sampleS1252;
  moonbit_string_t _M0L14phylip__sampleS1253;
  moonbit_string_t _M0L19phylip__seq__sampleS1254;
  moonbit_string_t _M0L20fasta__align__sampleS1255;
  struct moonbit_result_1 _tmp_3278;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L14clustal__alignS1256;
  moonbit_string_t _M0L6_2atmpS2720;
  struct moonbit_result_1 _tmp_3280;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L13phylip__alignS1257;
  moonbit_string_t _M0L6_2atmpS2721;
  struct moonbit_result_1 _tmp_3282;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L18phylip__seq__alignS1258;
  moonbit_string_t _M0L6_2atmpS2722;
  struct moonbit_result_1 _tmp_3284;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L12fasta__alignS1259;
  moonbit_string_t _M0L6_2atmpS2723;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2789;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2784;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2788;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2785;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2787;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2786;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2783;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L20write__clustal__recsS1260;
  struct _M0TPB3MapGssE* _M0L6_2atmpS2779;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L6_2atmpS2780;
  struct moonbit_result_1 _tmp_3286;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L21write__clustal__alignS1261;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2776;
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2775;
  struct moonbit_result_2 _tmp_3288;
  moonbit_string_t _M0L12clustal__outS1262;
  moonbit_string_t _M0L6_2atmpS2724;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2774;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2769;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2773;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2770;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2772;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2771;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2768;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L19write__phylip__recsS1263;
  struct _M0TPB3MapGssE* _M0L6_2atmpS2764;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L6_2atmpS2765;
  struct moonbit_result_1 _tmp_3290;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L20write__phylip__alignS1264;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2761;
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2760;
  struct moonbit_result_2 _tmp_3292;
  moonbit_string_t _M0L11phylip__outS1265;
  moonbit_string_t _M0L6_2atmpS2725;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2759;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2756;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2758;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2757;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2755;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L18write__fasta__recsS1266;
  struct _M0TPB3MapGssE* _M0L6_2atmpS2751;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L6_2atmpS2752;
  struct moonbit_result_1 _tmp_3294;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L19write__fasta__alignS1267;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2748;
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2747;
  struct moonbit_result_2 _tmp_3296;
  moonbit_string_t _M0L10fasta__outS1268;
  moonbit_string_t _M0L6_2atmpS2726;
  int32_t _M0L6_2atmpS2734;
  moonbit_string_t _M0L6_2atmpS2733;
  moonbit_string_t _M0L6_2atmpS2732;
  moonbit_string_t _M0L6_2atmpS2729;
  int32_t _M0L6_2atmpS2731;
  moonbit_string_t _M0L6_2atmpS2730;
  moonbit_string_t _M0L6_2atmpS2728;
  moonbit_string_t _M0L6_2atmpS2727;
  moonbit_string_t _M0L1cS1271;
  moonbit_string_t _M0L4col0S1269;
  moonbit_string_t _M0L7_2abindS1272;
  moonbit_string_t _M0L1cS1277;
  moonbit_string_t _M0L4col5S1275;
  moonbit_string_t _M0L7_2abindS1278;
  moonbit_string_t _M0L6_2atmpS2740;
  moonbit_string_t _M0L6_2atmpS2739;
  moonbit_string_t _M0L6_2atmpS2737;
  moonbit_string_t _M0L6_2atmpS2738;
  moonbit_string_t _M0L6_2atmpS2736;
  moonbit_string_t _M0L6_2atmpS2735;
  void* _M0L11_2atry__errS1283;
  moonbit_string_t _M0L12err__payloadS1281;
  moonbit_string_t _M0L6_2atmpS2744;
  moonbit_string_t _M0L6_2atmpS2743;
  struct moonbit_result_1 _tmp_3301;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2742;
  int32_t _M0L6_2atmpS2741;
  struct moonbit_result_0 _result_3302;
  #line 73 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L15clustal__sampleS1252
  = (moonbit_string_t)moonbit_string_literal_10.data;
  _M0L14phylip__sampleS1253
  = (moonbit_string_t)moonbit_string_literal_11.data;
  _M0L19phylip__seq__sampleS1254
  = (moonbit_string_t)moonbit_string_literal_11.data;
  _M0L20fasta__align__sampleS1255
  = (moonbit_string_t)moonbit_string_literal_12.data;
  #line 82 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3278
  = _M0FP26biolab8bio__seq13alignio__read(_M0L15clustal__sampleS1252, (moonbit_string_t)moonbit_string_literal_13.data);
  if (_tmp_3278.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2796 =
      _tmp_3278.data.ok;
    _M0L14clustal__alignS1256 = _M0L5_2aokS2796;
  } else {
    void* const _M0L6_2aerrS2797 = _tmp_3278.data.err;
    struct moonbit_result_0 _result_3279;
    moonbit_decref(_M0L20fasta__align__sampleS1255);
    moonbit_decref(_M0L19phylip__seq__sampleS1254);
    moonbit_decref(_M0L14phylip__sampleS1253);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3279.tag = 0;
    _result_3279.data.err = _M0L6_2aerrS2797;
    return _result_3279;
  }
  #line 83 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2720
  = _M0FP46biolab8bio__seq3cmd13alignio__main19alignment__to__json(_M0L14clustal__alignS1256);
  #line 83 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_14.data, _M0L6_2atmpS2720);
  moonbit_decref(_M0L6_2atmpS2720);
  #line 85 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3280
  = _M0FP26biolab8bio__seq13alignio__read(_M0L14phylip__sampleS1253, (moonbit_string_t)moonbit_string_literal_15.data);
  moonbit_decref(_M0L14phylip__sampleS1253);
  if (_tmp_3280.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2794 =
      _tmp_3280.data.ok;
    _M0L13phylip__alignS1257 = _M0L5_2aokS2794;
  } else {
    void* const _M0L6_2aerrS2795 = _tmp_3280.data.err;
    struct moonbit_result_0 _result_3281;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L20fasta__align__sampleS1255);
    moonbit_decref(_M0L19phylip__seq__sampleS1254);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3281.tag = 0;
    _result_3281.data.err = _M0L6_2aerrS2795;
    return _result_3281;
  }
  #line 86 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2721
  = _M0FP46biolab8bio__seq3cmd13alignio__main19alignment__to__json(_M0L13phylip__alignS1257);
  moonbit_decref(_M0L13phylip__alignS1257);
  #line 86 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_16.data, _M0L6_2atmpS2721);
  moonbit_decref(_M0L6_2atmpS2721);
  #line 88 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3282
  = _M0FP26biolab8bio__seq13alignio__read(_M0L19phylip__seq__sampleS1254, (moonbit_string_t)moonbit_string_literal_17.data);
  moonbit_decref(_M0L19phylip__seq__sampleS1254);
  if (_tmp_3282.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2792 =
      _tmp_3282.data.ok;
    _M0L18phylip__seq__alignS1258 = _M0L5_2aokS2792;
  } else {
    void* const _M0L6_2aerrS2793 = _tmp_3282.data.err;
    struct moonbit_result_0 _result_3283;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L20fasta__align__sampleS1255);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3283.tag = 0;
    _result_3283.data.err = _M0L6_2aerrS2793;
    return _result_3283;
  }
  #line 91 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2722
  = _M0FP46biolab8bio__seq3cmd13alignio__main19alignment__to__json(_M0L18phylip__seq__alignS1258);
  moonbit_decref(_M0L18phylip__seq__alignS1258);
  #line 91 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_18.data, _M0L6_2atmpS2722);
  moonbit_decref(_M0L6_2atmpS2722);
  #line 93 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3284
  = _M0FP26biolab8bio__seq13alignio__read(_M0L20fasta__align__sampleS1255, (moonbit_string_t)moonbit_string_literal_19.data);
  moonbit_decref(_M0L20fasta__align__sampleS1255);
  if (_tmp_3284.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2790 =
      _tmp_3284.data.ok;
    _M0L12fasta__alignS1259 = _M0L5_2aokS2790;
  } else {
    void* const _M0L6_2aerrS2791 = _tmp_3284.data.err;
    struct moonbit_result_0 _result_3285;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3285.tag = 0;
    _result_3285.data.err = _M0L6_2aerrS2791;
    return _result_3285;
  }
  #line 94 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2723
  = _M0FP46biolab8bio__seq3cmd13alignio__main19alignment__to__json(_M0L12fasta__alignS1259);
  moonbit_decref(_M0L12fasta__alignS1259);
  #line 94 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_20.data, _M0L6_2atmpS2723);
  moonbit_decref(_M0L6_2atmpS2723);
  #line 98 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2789
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_21.data);
  #line 97 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2784
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2789, (moonbit_string_t)moonbit_string_literal_22.data, (moonbit_string_t)moonbit_string_literal_22.data, (moonbit_string_t)moonbit_string_literal_22.data);
  moonbit_decref(_M0L6_2atmpS2789);
  #line 104 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2788
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_21.data);
  #line 103 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2785
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2788, (moonbit_string_t)moonbit_string_literal_23.data, (moonbit_string_t)moonbit_string_literal_23.data, (moonbit_string_t)moonbit_string_literal_23.data);
  moonbit_decref(_M0L6_2atmpS2788);
  #line 110 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2787
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_24.data);
  #line 109 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2786
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2787, (moonbit_string_t)moonbit_string_literal_25.data, (moonbit_string_t)moonbit_string_literal_25.data, (moonbit_string_t)moonbit_string_literal_25.data);
  moonbit_decref(_M0L6_2atmpS2787);
  _M0L6_2atmpS2783
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(3);
  _M0L6_2atmpS2783[0] = _M0L6_2atmpS2784;
  _M0L6_2atmpS2783[1] = _M0L6_2atmpS2785;
  _M0L6_2atmpS2783[2] = _M0L6_2atmpS2786;
  _M0L20write__clustal__recsS1260
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L20write__clustal__recsS1260)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L20write__clustal__recsS1260->$0 = _M0L6_2atmpS2783;
  _M0L20write__clustal__recsS1260->$1 = 3;
  _M0L6_2atmpS2779 = 0;
  _M0L6_2atmpS2780 = 0;
  #line 116 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3286
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(_M0L20write__clustal__recsS1260, _M0L6_2atmpS2779, _M0L6_2atmpS2780);
  moonbit_decref(_M0L20write__clustal__recsS1260);
  if (_M0L6_2atmpS2779) {
    moonbit_decref(_M0L6_2atmpS2779);
  }
  if (_M0L6_2atmpS2780) {
    moonbit_decref(_M0L6_2atmpS2780);
  }
  if (_tmp_3286.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2781 =
      _tmp_3286.data.ok;
    _M0L21write__clustal__alignS1261 = _M0L5_2aokS2781;
  } else {
    void* const _M0L6_2aerrS2782 = _tmp_3286.data.err;
    struct moonbit_result_0 _result_3287;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3287.tag = 0;
    _result_3287.data.err = _M0L6_2aerrS2782;
    return _result_3287;
  }
  _M0L6_2atmpS2776
  = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS2776[0] = _M0L21write__clustal__alignS1261;
  _M0L6_2atmpS2775
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
  Moonbit_object_header(_M0L6_2atmpS2775)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
  _M0L6_2atmpS2775->$0 = _M0L6_2atmpS2776;
  _M0L6_2atmpS2775->$1 = 1;
  #line 117 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3288
  = _M0FP26biolab8bio__seq14alignio__write(_M0L6_2atmpS2775, (moonbit_string_t)moonbit_string_literal_13.data);
  moonbit_decref(_M0L6_2atmpS2775);
  if (_tmp_3288.tag) {
    moonbit_string_t const _M0L5_2aokS2777 = _tmp_3288.data.ok;
    _M0L12clustal__outS1262 = _M0L5_2aokS2777;
  } else {
    void* const _M0L6_2aerrS2778 = _tmp_3288.data.err;
    struct moonbit_result_0 _result_3289;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3289.tag = 0;
    _result_3289.data.err = _M0L6_2aerrS2778;
    return _result_3289;
  }
  #line 118 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2724
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L12clustal__outS1262);
  moonbit_decref(_M0L12clustal__outS1262);
  #line 118 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_26.data, _M0L6_2atmpS2724);
  moonbit_decref(_M0L6_2atmpS2724);
  #line 122 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2774
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_27.data);
  #line 121 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2769
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2774, (moonbit_string_t)moonbit_string_literal_28.data, (moonbit_string_t)moonbit_string_literal_28.data, (moonbit_string_t)moonbit_string_literal_28.data);
  moonbit_decref(_M0L6_2atmpS2774);
  #line 128 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2773
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_27.data);
  #line 127 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2770
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2773, (moonbit_string_t)moonbit_string_literal_29.data, (moonbit_string_t)moonbit_string_literal_29.data, (moonbit_string_t)moonbit_string_literal_29.data);
  moonbit_decref(_M0L6_2atmpS2773);
  #line 134 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2772
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_27.data);
  #line 133 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2771
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2772, (moonbit_string_t)moonbit_string_literal_30.data, (moonbit_string_t)moonbit_string_literal_30.data, (moonbit_string_t)moonbit_string_literal_30.data);
  moonbit_decref(_M0L6_2atmpS2772);
  _M0L6_2atmpS2768
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(3);
  _M0L6_2atmpS2768[0] = _M0L6_2atmpS2769;
  _M0L6_2atmpS2768[1] = _M0L6_2atmpS2770;
  _M0L6_2atmpS2768[2] = _M0L6_2atmpS2771;
  _M0L19write__phylip__recsS1263
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L19write__phylip__recsS1263)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L19write__phylip__recsS1263->$0 = _M0L6_2atmpS2768;
  _M0L19write__phylip__recsS1263->$1 = 3;
  _M0L6_2atmpS2764 = 0;
  _M0L6_2atmpS2765 = 0;
  #line 140 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3290
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(_M0L19write__phylip__recsS1263, _M0L6_2atmpS2764, _M0L6_2atmpS2765);
  moonbit_decref(_M0L19write__phylip__recsS1263);
  if (_M0L6_2atmpS2764) {
    moonbit_decref(_M0L6_2atmpS2764);
  }
  if (_M0L6_2atmpS2765) {
    moonbit_decref(_M0L6_2atmpS2765);
  }
  if (_tmp_3290.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2766 =
      _tmp_3290.data.ok;
    _M0L20write__phylip__alignS1264 = _M0L5_2aokS2766;
  } else {
    void* const _M0L6_2aerrS2767 = _tmp_3290.data.err;
    struct moonbit_result_0 _result_3291;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3291.tag = 0;
    _result_3291.data.err = _M0L6_2aerrS2767;
    return _result_3291;
  }
  _M0L6_2atmpS2761
  = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS2761[0] = _M0L20write__phylip__alignS1264;
  _M0L6_2atmpS2760
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
  Moonbit_object_header(_M0L6_2atmpS2760)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
  _M0L6_2atmpS2760->$0 = _M0L6_2atmpS2761;
  _M0L6_2atmpS2760->$1 = 1;
  #line 141 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3292
  = _M0FP26biolab8bio__seq14alignio__write(_M0L6_2atmpS2760, (moonbit_string_t)moonbit_string_literal_15.data);
  moonbit_decref(_M0L6_2atmpS2760);
  if (_tmp_3292.tag) {
    moonbit_string_t const _M0L5_2aokS2762 = _tmp_3292.data.ok;
    _M0L11phylip__outS1265 = _M0L5_2aokS2762;
  } else {
    void* const _M0L6_2aerrS2763 = _tmp_3292.data.err;
    struct moonbit_result_0 _result_3293;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3293.tag = 0;
    _result_3293.data.err = _M0L6_2aerrS2763;
    return _result_3293;
  }
  #line 142 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2725
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L11phylip__outS1265);
  moonbit_decref(_M0L11phylip__outS1265);
  #line 142 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_31.data, _M0L6_2atmpS2725);
  moonbit_decref(_M0L6_2atmpS2725);
  #line 146 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2759
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_32.data);
  #line 145 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2756
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2759, (moonbit_string_t)moonbit_string_literal_28.data, (moonbit_string_t)moonbit_string_literal_28.data, (moonbit_string_t)moonbit_string_literal_33.data);
  moonbit_decref(_M0L6_2atmpS2759);
  #line 152 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2758
  = _M0MP26biolab8bio__seq3Seq3new((moonbit_string_t)moonbit_string_literal_32.data);
  #line 151 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2757
  = _M0MP26biolab8bio__seq9SeqRecord11new_2einner(_M0L6_2atmpS2758, (moonbit_string_t)moonbit_string_literal_29.data, (moonbit_string_t)moonbit_string_literal_29.data, (moonbit_string_t)moonbit_string_literal_34.data);
  moonbit_decref(_M0L6_2atmpS2758);
  _M0L6_2atmpS2755
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_raw(2);
  _M0L6_2atmpS2755[0] = _M0L6_2atmpS2756;
  _M0L6_2atmpS2755[1] = _M0L6_2atmpS2757;
  _M0L18write__fasta__recsS1266
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
  Moonbit_object_header(_M0L18write__fasta__recsS1266)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
  _M0L18write__fasta__recsS1266->$0 = _M0L6_2atmpS2755;
  _M0L18write__fasta__recsS1266->$1 = 2;
  _M0L6_2atmpS2751 = 0;
  _M0L6_2atmpS2752 = 0;
  #line 158 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3294
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(_M0L18write__fasta__recsS1266, _M0L6_2atmpS2751, _M0L6_2atmpS2752);
  moonbit_decref(_M0L18write__fasta__recsS1266);
  if (_M0L6_2atmpS2751) {
    moonbit_decref(_M0L6_2atmpS2751);
  }
  if (_M0L6_2atmpS2752) {
    moonbit_decref(_M0L6_2atmpS2752);
  }
  if (_tmp_3294.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2753 =
      _tmp_3294.data.ok;
    _M0L19write__fasta__alignS1267 = _M0L5_2aokS2753;
  } else {
    void* const _M0L6_2aerrS2754 = _tmp_3294.data.err;
    struct moonbit_result_0 _result_3295;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3295.tag = 0;
    _result_3295.data.err = _M0L6_2aerrS2754;
    return _result_3295;
  }
  _M0L6_2atmpS2748
  = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
  _M0L6_2atmpS2748[0] = _M0L19write__fasta__alignS1267;
  _M0L6_2atmpS2747
  = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
  Moonbit_object_header(_M0L6_2atmpS2747)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
  _M0L6_2atmpS2747->$0 = _M0L6_2atmpS2748;
  _M0L6_2atmpS2747->$1 = 1;
  #line 159 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3296
  = _M0FP26biolab8bio__seq14alignio__write(_M0L6_2atmpS2747, (moonbit_string_t)moonbit_string_literal_19.data);
  moonbit_decref(_M0L6_2atmpS2747);
  if (_tmp_3296.tag) {
    moonbit_string_t const _M0L5_2aokS2749 = _tmp_3296.data.ok;
    _M0L10fasta__outS1268 = _M0L5_2aokS2749;
  } else {
    void* const _M0L6_2aerrS2750 = _tmp_3296.data.err;
    struct moonbit_result_0 _result_3297;
    moonbit_decref(_M0L14clustal__alignS1256);
    moonbit_decref(_M0L15clustal__sampleS1252);
    _result_3297.tag = 0;
    _result_3297.data.err = _M0L6_2aerrS2750;
    return _result_3297;
  }
  #line 160 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2726
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L10fasta__outS1268);
  moonbit_decref(_M0L10fasta__outS1268);
  #line 160 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_35.data, _M0L6_2atmpS2726);
  moonbit_decref(_M0L6_2atmpS2726);
  #line 165 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2734
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment12num__records(_M0L14clustal__alignS1256);
  #line 165 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2733
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2734, 10);
  #line 164 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2732
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_36.data, _M0L6_2atmpS2733);
  moonbit_decref(_M0L6_2atmpS2733);
  #line 164 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2729
  = moonbit_add_string(_M0L6_2atmpS2732, (moonbit_string_t)moonbit_string_literal_37.data);
  moonbit_decref(_M0L6_2atmpS2732);
  #line 167 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2731
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(_M0L14clustal__alignS1256);
  #line 167 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2730
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2731, 10);
  #line 164 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2728 = moonbit_add_string(_M0L6_2atmpS2729, _M0L6_2atmpS2730);
  moonbit_decref(_M0L6_2atmpS2730);
  moonbit_decref(_M0L6_2atmpS2729);
  #line 164 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2727
  = moonbit_add_string(_M0L6_2atmpS2728, (moonbit_string_t)moonbit_string_literal_6.data);
  moonbit_decref(_M0L6_2atmpS2728);
  #line 162 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_38.data, _M0L6_2atmpS2727);
  moonbit_decref(_M0L6_2atmpS2727);
  #line 171 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L7_2abindS1272
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment6column(_M0L14clustal__alignS1256, 0);
  if (_M0L7_2abindS1272 == 0) {
    if (_M0L7_2abindS1272) {
      moonbit_decref(_M0L7_2abindS1272);
    }
    _M0L4col0S1269 = (moonbit_string_t)moonbit_string_literal_0.data;
  } else {
    moonbit_string_t _M0L7_2aSomeS1273 = _M0L7_2abindS1272;
    moonbit_string_t _M0L4_2acS1274 = _M0L7_2aSomeS1273;
    _M0L1cS1271 = _M0L4_2acS1274;
    goto join_1270;
  }
  goto joinlet_3298;
  join_1270:;
  _M0L4col0S1269 = _M0L1cS1271;
  joinlet_3298:;
  #line 175 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L7_2abindS1278
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment6column(_M0L14clustal__alignS1256, 5);
  moonbit_decref(_M0L14clustal__alignS1256);
  if (_M0L7_2abindS1278 == 0) {
    if (_M0L7_2abindS1278) {
      moonbit_decref(_M0L7_2abindS1278);
    }
    _M0L4col5S1275 = (moonbit_string_t)moonbit_string_literal_0.data;
  } else {
    moonbit_string_t _M0L7_2aSomeS1279 = _M0L7_2abindS1278;
    moonbit_string_t _M0L4_2acS1280 = _M0L7_2aSomeS1279;
    _M0L1cS1277 = _M0L4_2acS1280;
    goto join_1276;
  }
  goto joinlet_3299;
  join_1276:;
  _M0L4col5S1275 = _M0L1cS1277;
  joinlet_3299:;
  #line 182 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2740
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L4col0S1269);
  moonbit_decref(_M0L4col0S1269);
  #line 181 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2739
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_39.data, _M0L6_2atmpS2740);
  moonbit_decref(_M0L6_2atmpS2740);
  #line 181 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2737
  = moonbit_add_string(_M0L6_2atmpS2739, (moonbit_string_t)moonbit_string_literal_40.data);
  moonbit_decref(_M0L6_2atmpS2739);
  #line 184 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2738
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L4col5S1275);
  moonbit_decref(_M0L4col5S1275);
  #line 181 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2736 = moonbit_add_string(_M0L6_2atmpS2737, _M0L6_2atmpS2738);
  moonbit_decref(_M0L6_2atmpS2738);
  moonbit_decref(_M0L6_2atmpS2737);
  #line 181 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2735
  = moonbit_add_string(_M0L6_2atmpS2736, (moonbit_string_t)moonbit_string_literal_6.data);
  moonbit_decref(_M0L6_2atmpS2736);
  #line 179 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_41.data, _M0L6_2atmpS2735);
  moonbit_decref(_M0L6_2atmpS2735);
  #line 189 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2744
  = moonbit_add_string(_M0L15clustal__sampleS1252, (moonbit_string_t)moonbit_string_literal_42.data);
  #line 189 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2743
  = moonbit_add_string(_M0L6_2atmpS2744, _M0L15clustal__sampleS1252);
  moonbit_decref(_M0L15clustal__sampleS1252);
  moonbit_decref(_M0L6_2atmpS2744);
  #line 189 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _tmp_3301
  = _M0FP26biolab8bio__seq13alignio__read(_M0L6_2atmpS2743, (moonbit_string_t)moonbit_string_literal_13.data);
  moonbit_decref(_M0L6_2atmpS2743);
  if (_tmp_3301.tag) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2745 =
      _tmp_3301.data.ok;
    _M0L6_2atmpS2742 = _M0L5_2aokS2745;
  } else {
    void* const _M0L6_2aerrS2746 = _tmp_3301.data.err;
    _M0L11_2atry__errS1283 = _M0L6_2aerrS2746;
    goto join_1282;
  }
  moonbit_decref(_M0L6_2atmpS2742);
  _M0L12err__payloadS1281 = (moonbit_string_t)moonbit_string_literal_43.data;
  goto joinlet_3300;
  join_1282:;
  moonbit_decref(_M0L11_2atry__errS1283);
  _M0L12err__payloadS1281 = (moonbit_string_t)moonbit_string_literal_44.data;
  joinlet_3300:;
  #line 195 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2741
  = _M0FP46biolab8bio__seq3cmd13alignio__main4emit((moonbit_string_t)moonbit_string_literal_45.data, _M0L12err__payloadS1281);
  moonbit_decref(_M0L12err__payloadS1281);
  _result_3302.tag = 1;
  _result_3302.data.ok = _M0L6_2atmpS2741;
  return _result_3302;
}

int32_t _M0FP46biolab8bio__seq3cmd13alignio__main4emit(
  moonbit_string_t _M0L4nameS1250,
  moonbit_string_t _M0L7payloadS1251
) {
  moonbit_string_t _M0L6_2atmpS2719;
  moonbit_string_t _M0L6_2atmpS2718;
  #line 66 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  #line 67 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2719
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_46.data, _M0L4nameS1250);
  #line 67 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2718
  = moonbit_add_string(_M0L6_2atmpS2719, (moonbit_string_t)moonbit_string_literal_47.data);
  moonbit_decref(_M0L6_2atmpS2719);
  #line 67 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS2718);
  moonbit_decref(_M0L6_2atmpS2718);
  #line 68 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FPB7printlnGsE(_M0L7payloadS1251);
  #line 69 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_48.data);
  return 0;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main19alignment__to__json(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1249
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1248;
  int32_t _M0L6_2atmpS2713;
  moonbit_string_t _M0L6_2atmpS2712;
  int32_t _M0L6_2atmpS2715;
  moonbit_string_t _M0L6_2atmpS2714;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L6_2atmpS2717;
  moonbit_string_t _M0L6_2atmpS2716;
  moonbit_string_t _result_3303;
  #line 53 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  #line 54 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L3bufS1248 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 55 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1248, (moonbit_string_t)moonbit_string_literal_36.data);
  #line 56 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2713
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment12num__records(_M0L5alignS1249);
  #line 56 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2712
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2713, 10);
  #line 56 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1248, _M0L6_2atmpS2712);
  moonbit_decref(_M0L6_2atmpS2712);
  #line 57 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1248, (moonbit_string_t)moonbit_string_literal_37.data);
  #line 58 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2715
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(_M0L5alignS1249);
  #line 58 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2714
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2715, 10);
  #line 58 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1248, _M0L6_2atmpS2714);
  moonbit_decref(_M0L6_2atmpS2714);
  #line 59 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1248, (moonbit_string_t)moonbit_string_literal_49.data);
  #line 60 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2717
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(_M0L5alignS1249);
  #line 60 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2716
  = _M0FP46biolab8bio__seq3cmd13alignio__main17records__to__json(_M0L6_2atmpS2717);
  moonbit_decref(_M0L6_2atmpS2717);
  #line 60 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1248, _M0L6_2atmpS2716);
  moonbit_decref(_M0L6_2atmpS2716);
  #line 61 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1248, 125);
  #line 62 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _result_3303 = _M0MPB13StringBuilder10to__string(_M0L3bufS1248);
  moonbit_decref(_M0L3bufS1248);
  return _result_3303;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main17records__to__json(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1244
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1242;
  int32_t _M0L7_2abindS1243;
  int32_t _M0L1iS1245;
  moonbit_string_t _result_3305;
  #line 39 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  #line 40 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L3bufS1242 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 41 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1242, 91);
  _M0L7_2abindS1243 = _M0L7recordsS1244->$1;
  _M0L1iS1245 = 0;
  while (1) {
    if (_M0L1iS1245 < _M0L7_2abindS1243) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2711 =
        _M0L7recordsS1244->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1246 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2711[_M0L1iS1245];
      moonbit_string_t _M0L6_2atmpS2709;
      int32_t _M0L6_2atmpS2710;
      if (_M0L1iS1245 > 0) {
        moonbit_incref(_M0L3recS1246);
        #line 44 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1242, (moonbit_string_t)moonbit_string_literal_50.data);
      } else {
        moonbit_incref(_M0L3recS1246);
      }
      #line 46 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0L6_2atmpS2709
      = _M0FP46biolab8bio__seq3cmd13alignio__main16record__to__json(_M0L3recS1246);
      moonbit_decref(_M0L3recS1246);
      #line 46 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1242, _M0L6_2atmpS2709);
      moonbit_decref(_M0L6_2atmpS2709);
      _M0L6_2atmpS2710 = _M0L1iS1245 + 1;
      _M0L1iS1245 = _M0L6_2atmpS2710;
      continue;
    }
    break;
  }
  #line 48 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1242, 93);
  #line 49 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _result_3305 = _M0MPB13StringBuilder10to__string(_M0L3bufS1242);
  moonbit_decref(_M0L3bufS1242);
  return _result_3305;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main16record__to__json(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1241
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1240;
  moonbit_string_t _M0L2idS2699;
  moonbit_string_t _M0L6_2atmpS2698;
  moonbit_string_t _M0L4nameS2701;
  moonbit_string_t _M0L6_2atmpS2700;
  moonbit_string_t _M0L11descriptionS2703;
  moonbit_string_t _M0L6_2atmpS2702;
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2706;
  moonbit_string_t _M0L6_2atmpS2705;
  moonbit_string_t _M0L6_2atmpS2704;
  int32_t _M0L6_2atmpS2708;
  moonbit_string_t _M0L6_2atmpS2707;
  moonbit_string_t _result_3306;
  #line 22 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  #line 23 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L3bufS1240 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 24 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, (moonbit_string_t)moonbit_string_literal_51.data);
  _M0L2idS2699 = _M0L3recS1241->$1;
  moonbit_incref(_M0L2idS2699);
  #line 25 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2698
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L2idS2699);
  moonbit_decref(_M0L2idS2699);
  #line 25 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, _M0L6_2atmpS2698);
  moonbit_decref(_M0L6_2atmpS2698);
  #line 26 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, (moonbit_string_t)moonbit_string_literal_52.data);
  _M0L4nameS2701 = _M0L3recS1241->$2;
  moonbit_incref(_M0L4nameS2701);
  #line 27 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2700
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L4nameS2701);
  moonbit_decref(_M0L4nameS2701);
  #line 27 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, _M0L6_2atmpS2700);
  moonbit_decref(_M0L6_2atmpS2700);
  #line 28 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, (moonbit_string_t)moonbit_string_literal_53.data);
  _M0L11descriptionS2703 = _M0L3recS1241->$3;
  moonbit_incref(_M0L11descriptionS2703);
  #line 29 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2702
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L11descriptionS2703);
  moonbit_decref(_M0L11descriptionS2703);
  #line 29 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, _M0L6_2atmpS2702);
  moonbit_decref(_M0L6_2atmpS2702);
  #line 30 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, (moonbit_string_t)moonbit_string_literal_54.data);
  _M0L3seqS2706 = _M0L3recS1241->$0;
  moonbit_incref(_M0L3seqS2706);
  #line 31 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2705 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2706);
  moonbit_decref(_M0L3seqS2706);
  #line 31 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2704
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(_M0L6_2atmpS2705);
  moonbit_decref(_M0L6_2atmpS2705);
  #line 31 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, _M0L6_2atmpS2704);
  moonbit_decref(_M0L6_2atmpS2704);
  #line 32 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, (moonbit_string_t)moonbit_string_literal_55.data);
  #line 33 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2708 = _M0MP26biolab8bio__seq9SeqRecord6length(_M0L3recS1241);
  #line 33 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2707
  = _M0MPC13int3Int18to__string_2einner(_M0L6_2atmpS2708, 10);
  #line 33 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1240, _M0L6_2atmpS2707);
  moonbit_decref(_M0L6_2atmpS2707);
  #line 34 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1240, 125);
  #line 35 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _result_3306 = _M0MPB13StringBuilder10to__string(_M0L3bufS1240);
  moonbit_decref(_M0L3bufS1240);
  return _result_3306;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main12json__string(
  moonbit_string_t _M0L1sS1239
) {
  moonbit_string_t _M0L6_2atmpS2697;
  moonbit_string_t _M0L6_2atmpS2696;
  moonbit_string_t _result_3307;
  #line 17 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  #line 18 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2697
  = _M0FP46biolab8bio__seq3cmd13alignio__main12json__escape(_M0L1sS1239);
  #line 18 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2696
  = moonbit_add_string((moonbit_string_t)moonbit_string_literal_56.data, _M0L6_2atmpS2697);
  moonbit_decref(_M0L6_2atmpS2697);
  #line 18 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _result_3307
  = moonbit_add_string(_M0L6_2atmpS2696, (moonbit_string_t)moonbit_string_literal_56.data);
  moonbit_decref(_M0L6_2atmpS2696);
  return _result_3307;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd13alignio__main12json__escape(
  moonbit_string_t _M0L1sS1231
) {
  int32_t _M0L6_2atmpS2695;
  int32_t _M0L6_2atmpS2694;
  struct _M0TPB13StringBuilder* _M0L3bufS1230;
  struct _M0TPB4IterGcE* _M0L5_2aitS1232;
  moonbit_string_t _result_3310;
  #line 2 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L6_2atmpS2695 = Moonbit_array_length(_M0L1sS1231);
  _M0L6_2atmpS2694 = _M0L6_2atmpS2695 + 2;
  #line 3 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L3bufS1230
  = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L6_2atmpS2694);
  #line 4 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _M0L5_2aitS1232 = _M0MPC16string6String4iter(_M0L1sS1231);
  while (1) {
    int32_t _M0L1cS1234;
    int32_t _M0L7_2abindS1236;
    #line 4 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
    _M0L7_2abindS1236 = _M0MPB4Iter4nextGcE(_M0L5_2aitS1232);
    if (_M0L7_2abindS1236 == -1) {
      moonbit_decref(_M0L5_2aitS1232);
    } else {
      int32_t _M0L7_2aSomeS1237 = _M0L7_2abindS1236;
      int32_t _M0L4_2acS1238 = _M0L7_2aSomeS1237;
      _M0L1cS1234 = _M0L4_2acS1238;
      goto join_1233;
    }
    goto joinlet_3309;
    join_1233:;
    if (_M0L1cS1234 == 92) {
      #line 6 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1230, (moonbit_string_t)moonbit_string_literal_57.data);
    } else if (_M0L1cS1234 == 34) {
      #line 7 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1230, (moonbit_string_t)moonbit_string_literal_58.data);
    } else if (_M0L1cS1234 == 10) {
      #line 8 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1230, (moonbit_string_t)moonbit_string_literal_59.data);
    } else if (_M0L1cS1234 == 9) {
      #line 9 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1230, (moonbit_string_t)moonbit_string_literal_60.data);
    } else {
      #line 10 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1230, _M0L1cS1234);
    }
    continue;
    joinlet_3309:;
    break;
  }
  #line 13 "/home/developer/Documents/1/cmd/alignio_main/alignio_tests.mbt"
  _result_3310 = _M0MPB13StringBuilder10to__string(_M0L3bufS1230);
  moonbit_decref(_M0L3bufS1230);
  return _result_3310;
}

moonbit_string_t _M0MP26biolab8bio__seq20MultipleSeqAlignment6column(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L4selfS1217,
  int32_t _M0L3colS1218
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS2686;
  int32_t _M0L6_2atmpS2685;
  int32_t _if__result_3311;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS2693;
  int32_t _M0L6_2atmpS2692;
  struct _M0TPB13StringBuilder* _M0L3bufS1219;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS1220;
  int32_t _M0L7_2abindS1221;
  int32_t _M0L2__S1222;
  moonbit_string_t _M0L6_2atmpS2691;
  #line 81 "/home/developer/Documents/1/align.mbt"
  _M0L7recordsS2686 = _M0L4selfS1217->$0;
  moonbit_incref(_M0L7recordsS2686);
  #line 85 "/home/developer/Documents/1/align.mbt"
  _M0L6_2atmpS2685
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS2686);
  moonbit_decref(_M0L7recordsS2686);
  if (_M0L6_2atmpS2685 == 0) {
    return 0;
  }
  if (_M0L3colS1218 < 0) {
    _if__result_3311 = 1;
  } else {
    int32_t _M0L6_2atmpS2687;
    #line 88 "/home/developer/Documents/1/align.mbt"
    _M0L6_2atmpS2687
    = _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(_M0L4selfS1217);
    _if__result_3311 = _M0L3colS1218 >= _M0L6_2atmpS2687;
  }
  if (_if__result_3311) {
    return 0;
  }
  _M0L7recordsS2693 = _M0L4selfS1217->$0;
  moonbit_incref(_M0L7recordsS2693);
  #line 91 "/home/developer/Documents/1/align.mbt"
  _M0L6_2atmpS2692
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS2693);
  moonbit_decref(_M0L7recordsS2693);
  #line 91 "/home/developer/Documents/1/align.mbt"
  _M0L3bufS1219
  = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L6_2atmpS2692);
  _M0L7_2abindS1220 = _M0L4selfS1217->$0;
  _M0L7_2abindS1221 = _M0L7_2abindS1220->$1;
  moonbit_incref(_M0L7_2abindS1220);
  _M0L2__S1222 = 0;
  while (1) {
    if (_M0L2__S1222 < _M0L7_2abindS1221) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2690 =
        _M0L7_2abindS1220->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1223 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2690[_M0L2__S1222];
      int32_t _M0L1cS1225;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2688 = _M0L3recS1223->$0;
      int32_t _M0L7_2abindS1226;
      int32_t _M0L6_2atmpS2689;
      moonbit_incref(_M0L3seqS2688);
      #line 93 "/home/developer/Documents/1/align.mbt"
      _M0L7_2abindS1226
      = _M0MP26biolab8bio__seq3Seq9get__char(_M0L3seqS2688, _M0L3colS1218);
      moonbit_decref(_M0L3seqS2688);
      if (_M0L7_2abindS1226 == -1) {
        #line 95 "/home/developer/Documents/1/align.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1219, 63);
      } else {
        int32_t _M0L7_2aSomeS1227 = _M0L7_2abindS1226;
        int32_t _M0L4_2acS1228 = _M0L7_2aSomeS1227;
        _M0L1cS1225 = _M0L4_2acS1228;
        goto join_1224;
      }
      goto joinlet_3313;
      join_1224:;
      #line 94 "/home/developer/Documents/1/align.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1219, _M0L1cS1225);
      joinlet_3313:;
      _M0L6_2atmpS2689 = _M0L2__S1222 + 1;
      _M0L2__S1222 = _M0L6_2atmpS2689;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1220);
    }
    break;
  }
  #line 98 "/home/developer/Documents/1/align.mbt"
  _M0L6_2atmpS2691 = _M0MPB13StringBuilder10to__string(_M0L3bufS1219);
  moonbit_decref(_M0L3bufS1219);
  return _M0L6_2atmpS2691;
}

int32_t _M0MP26biolab8bio__seq20MultipleSeqAlignment12num__records(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L4selfS1216
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS2684;
  int32_t _result_3314;
  #line 50 "/home/developer/Documents/1/align.mbt"
  _M0L7recordsS2684 = _M0L4selfS1216->$0;
  moonbit_incref(_M0L7recordsS2684);
  #line 51 "/home/developer/Documents/1/align.mbt"
  _result_3314
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS2684);
  moonbit_decref(_M0L7recordsS2684);
  return _result_3314;
}

struct moonbit_result_2 _M0FP26biolab8bio__seq14alignio__write(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L10alignmentsS1215,
  moonbit_string_t _M0L6formatS1214
) {
  #line 46 "/home/developer/Documents/1/alignio.mbt"
  if (
    _M0L6formatS1214 == (moonbit_string_t)moonbit_string_literal_13.data
    || Moonbit_array_length(_M0L6formatS1214) == 7
       && 0
          == memcmp(_M0L6formatS1214, (moonbit_string_t)moonbit_string_literal_13.data, 14)
  ) {
    moonbit_string_t _M0L6_2atmpS2678;
    struct moonbit_result_2 _result_3319;
    #line 51 "/home/developer/Documents/1/alignio.mbt"
    _M0L6_2atmpS2678
    = _M0FP26biolab8bio__seq24write__clustal__multiple(_M0L10alignmentsS1215);
    _result_3319.tag = 1;
    _result_3319.data.ok = _M0L6_2atmpS2678;
    return _result_3319;
  } else if (
           _M0L6formatS1214
           == (moonbit_string_t)moonbit_string_literal_15.data
           || Moonbit_array_length(_M0L6formatS1214) == 6
              && 0
                 == memcmp(_M0L6formatS1214, (moonbit_string_t)moonbit_string_literal_15.data, 12)
         ) {
    moonbit_string_t _M0L6_2atmpS2679;
    struct moonbit_result_2 _result_3318;
    #line 52 "/home/developer/Documents/1/alignio.mbt"
    _M0L6_2atmpS2679
    = _M0FP26biolab8bio__seq23write__phylip__multiple(_M0L10alignmentsS1215);
    _result_3318.tag = 1;
    _result_3318.data.ok = _M0L6_2atmpS2679;
    return _result_3318;
  } else if (
           _M0L6formatS1214
           == (moonbit_string_t)moonbit_string_literal_17.data
           || Moonbit_array_length(_M0L6formatS1214) == 17
              && 0
                 == memcmp(_M0L6formatS1214, (moonbit_string_t)moonbit_string_literal_17.data, 34)
         ) {
    moonbit_string_t _M0L6_2atmpS2680;
    struct moonbit_result_2 _result_3317;
    #line 53 "/home/developer/Documents/1/alignio.mbt"
    _M0L6_2atmpS2680
    = _M0FP26biolab8bio__seq35write__phylip__sequential__multiple(_M0L10alignmentsS1215);
    _result_3317.tag = 1;
    _result_3317.data.ok = _M0L6_2atmpS2680;
    return _result_3317;
  } else if (
           _M0L6formatS1214
           == (moonbit_string_t)moonbit_string_literal_19.data
           || Moonbit_array_length(_M0L6formatS1214) == 5
              && 0
                 == memcmp(_M0L6formatS1214, (moonbit_string_t)moonbit_string_literal_19.data, 10)
         ) {
    moonbit_string_t _M0L6_2atmpS2681;
    struct moonbit_result_2 _result_3316;
    #line 54 "/home/developer/Documents/1/alignio.mbt"
    _M0L6_2atmpS2681
    = _M0FP26biolab8bio__seq27write__fasta__as__alignment(_M0L10alignmentsS1215);
    _result_3316.tag = 1;
    _result_3316.data.ok = _M0L6_2atmpS2681;
    return _result_3316;
  } else {
    moonbit_string_t _M0L6_2atmpS2683;
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2682;
    struct moonbit_result_2 _result_3315;
    #line 55 "/home/developer/Documents/1/alignio.mbt"
    _M0L6_2atmpS2683
    = moonbit_add_string((moonbit_string_t)moonbit_string_literal_61.data, _M0L6formatS1214);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2682
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2682)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2682)->$0
    = _M0L6_2atmpS2683;
    _result_3315.tag = 0;
    _result_3315.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2682;
    return _result_3315;
  }
}

moonbit_string_t _M0FP26biolab8bio__seq35write__phylip__sequential__multiple(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L10alignmentsS1210
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1208;
  int32_t _M0L7_2abindS1209;
  int32_t _M0L1iS1211;
  moonbit_string_t _result_3321;
  #line 108 "/home/developer/Documents/1/alignio.mbt"
  #line 111 "/home/developer/Documents/1/alignio.mbt"
  _M0L3bufS1208 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS1209 = _M0L10alignmentsS1210->$1;
  _M0L1iS1211 = 0;
  while (1) {
    if (_M0L1iS1211 < _M0L7_2abindS1209) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L3bufS2677 =
        _M0L10alignmentsS1210->$0;
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1212 =
        (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*)_M0L3bufS2677[
          _M0L1iS1211
        ];
      moonbit_string_t _M0L6_2atmpS2675;
      int32_t _M0L6_2atmpS2676;
      if (_M0L1iS1211 > 0) {
        moonbit_incref(_M0L5alignS1212);
        #line 114 "/home/developer/Documents/1/alignio.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1208, 10);
      } else {
        moonbit_incref(_M0L5alignS1212);
      }
      #line 116 "/home/developer/Documents/1/alignio.mbt"
      _M0L6_2atmpS2675
      = _M0FP26biolab8bio__seq25write__phylip__sequential(_M0L5alignS1212);
      moonbit_decref(_M0L5alignS1212);
      #line 116 "/home/developer/Documents/1/alignio.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1208, _M0L6_2atmpS2675);
      moonbit_decref(_M0L6_2atmpS2675);
      _M0L6_2atmpS2676 = _M0L1iS1211 + 1;
      _M0L1iS1211 = _M0L6_2atmpS2676;
      continue;
    }
    break;
  }
  #line 118 "/home/developer/Documents/1/alignio.mbt"
  _result_3321 = _M0MPB13StringBuilder10to__string(_M0L3bufS1208);
  moonbit_decref(_M0L3bufS1208);
  return _result_3321;
}

moonbit_string_t _M0FP26biolab8bio__seq23write__phylip__multiple(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L10alignmentsS1204
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1202;
  int32_t _M0L7_2abindS1203;
  int32_t _M0L1iS1205;
  moonbit_string_t _result_3323;
  #line 95 "/home/developer/Documents/1/alignio.mbt"
  #line 96 "/home/developer/Documents/1/alignio.mbt"
  _M0L3bufS1202 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS1203 = _M0L10alignmentsS1204->$1;
  _M0L1iS1205 = 0;
  while (1) {
    if (_M0L1iS1205 < _M0L7_2abindS1203) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L3bufS2674 =
        _M0L10alignmentsS1204->$0;
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1206 =
        (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*)_M0L3bufS2674[
          _M0L1iS1205
        ];
      moonbit_string_t _M0L6_2atmpS2672;
      int32_t _M0L6_2atmpS2673;
      if (_M0L1iS1205 > 0) {
        moonbit_incref(_M0L5alignS1206);
        #line 99 "/home/developer/Documents/1/alignio.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1202, 10);
      } else {
        moonbit_incref(_M0L5alignS1206);
      }
      #line 101 "/home/developer/Documents/1/alignio.mbt"
      _M0L6_2atmpS2672
      = _M0FP26biolab8bio__seq13write__phylip(_M0L5alignS1206);
      moonbit_decref(_M0L5alignS1206);
      #line 101 "/home/developer/Documents/1/alignio.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1202, _M0L6_2atmpS2672);
      moonbit_decref(_M0L6_2atmpS2672);
      _M0L6_2atmpS2673 = _M0L1iS1205 + 1;
      _M0L1iS1205 = _M0L6_2atmpS2673;
      continue;
    }
    break;
  }
  #line 103 "/home/developer/Documents/1/alignio.mbt"
  _result_3323 = _M0MPB13StringBuilder10to__string(_M0L3bufS1202);
  moonbit_decref(_M0L3bufS1202);
  return _result_3323;
}

moonbit_string_t _M0FP26biolab8bio__seq27write__fasta__as__alignment(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L10alignmentsS1193
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L12all__recordsS1191;
  int32_t _M0L7_2abindS1192;
  int32_t _M0L2__S1194;
  moonbit_string_t _result_3326;
  #line 70 "/home/developer/Documents/1/alignio.mbt"
  #line 71 "/home/developer/Documents/1/alignio.mbt"
  _M0L12all__recordsS1191
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS1192 = _M0L10alignmentsS1193->$1;
  _M0L2__S1194 = 0;
  while (1) {
    if (_M0L2__S1194 < _M0L7_2abindS1192) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L3bufS2671 =
        _M0L10alignmentsS1193->$0;
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1195 =
        (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*)_M0L3bufS2671[
          _M0L2__S1194
        ];
      struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7_2abindS1196;
      int32_t _M0L7_2abindS1197;
      int32_t _M0L2__S1198;
      int32_t _M0L6_2atmpS2670;
      moonbit_incref(_M0L5alignS1195);
      #line 73 "/home/developer/Documents/1/alignio.mbt"
      _M0L7_2abindS1196
      = _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(_M0L5alignS1195);
      moonbit_decref(_M0L5alignS1195);
      _M0L7_2abindS1197 = _M0L7_2abindS1196->$1;
      _M0L2__S1198 = 0;
      while (1) {
        if (_M0L2__S1198 < _M0L7_2abindS1197) {
          struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2669 =
            _M0L7_2abindS1196->$0;
          struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1199 =
            (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2669[
              _M0L2__S1198
            ];
          int32_t _M0L6_2atmpS2668;
          moonbit_incref(_M0L3recS1199);
          #line 74 "/home/developer/Documents/1/alignio.mbt"
          _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L12all__recordsS1191, _M0L3recS1199);
          moonbit_decref(_M0L3recS1199);
          _M0L6_2atmpS2668 = _M0L2__S1198 + 1;
          _M0L2__S1198 = _M0L6_2atmpS2668;
          continue;
        } else {
          moonbit_decref(_M0L7_2abindS1196);
        }
        break;
      }
      _M0L6_2atmpS2670 = _M0L2__S1194 + 1;
      _M0L2__S1194 = _M0L6_2atmpS2670;
      continue;
    }
    break;
  }
  #line 77 "/home/developer/Documents/1/alignio.mbt"
  _result_3326
  = _M0FP26biolab8bio__seq12write__fasta(_M0L12all__recordsS1191);
  moonbit_decref(_M0L12all__recordsS1191);
  return _result_3326;
}

moonbit_string_t _M0FP26biolab8bio__seq24write__clustal__multiple(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L10alignmentsS1187
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1185;
  int32_t _M0L7_2abindS1186;
  int32_t _M0L1iS1188;
  moonbit_string_t _result_3328;
  #line 82 "/home/developer/Documents/1/alignio.mbt"
  #line 83 "/home/developer/Documents/1/alignio.mbt"
  _M0L3bufS1185 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS1186 = _M0L10alignmentsS1187->$1;
  _M0L1iS1188 = 0;
  while (1) {
    if (_M0L1iS1188 < _M0L7_2abindS1186) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L3bufS2667 =
        _M0L10alignmentsS1187->$0;
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1189 =
        (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*)_M0L3bufS2667[
          _M0L1iS1188
        ];
      moonbit_string_t _M0L6_2atmpS2665;
      int32_t _M0L6_2atmpS2666;
      if (_M0L1iS1188 > 0) {
        moonbit_incref(_M0L5alignS1189);
        #line 86 "/home/developer/Documents/1/alignio.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1185, 10);
      } else {
        moonbit_incref(_M0L5alignS1189);
      }
      #line 88 "/home/developer/Documents/1/alignio.mbt"
      _M0L6_2atmpS2665
      = _M0FP26biolab8bio__seq14write__clustal(_M0L5alignS1189);
      moonbit_decref(_M0L5alignS1189);
      #line 88 "/home/developer/Documents/1/alignio.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1185, _M0L6_2atmpS2665);
      moonbit_decref(_M0L6_2atmpS2665);
      _M0L6_2atmpS2666 = _M0L1iS1188 + 1;
      _M0L1iS1188 = _M0L6_2atmpS2666;
      continue;
    }
    break;
  }
  #line 90 "/home/developer/Documents/1/alignio.mbt"
  _result_3328 = _M0MPB13StringBuilder10to__string(_M0L3bufS1185);
  moonbit_decref(_M0L3bufS1185);
  return _result_3328;
}

struct moonbit_result_1 _M0FP26biolab8bio__seq13alignio__read(
  moonbit_string_t _M0L7contentS1183,
  moonbit_string_t _M0L6formatS1184
) {
  struct moonbit_result_3 _tmp_3329;
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L10alignmentsS1182;
  int32_t _M0L6_2atmpS2658;
  int32_t _M0L6_2atmpS2660;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2662;
  struct moonbit_result_1 _result_3333;
  #line 28 "/home/developer/Documents/1/alignio.mbt"
  #line 32 "/home/developer/Documents/1/alignio.mbt"
  _tmp_3329
  = _M0FP26biolab8bio__seq14alignio__parse(_M0L7contentS1183, _M0L6formatS1184);
  if (_tmp_3329.tag) {
    struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* const _M0L5_2aokS2663 =
      _tmp_3329.data.ok;
    _M0L10alignmentsS1182 = _M0L5_2aokS2663;
  } else {
    void* const _M0L6_2aerrS2664 = _tmp_3329.data.err;
    struct moonbit_result_1 _result_3330;
    _result_3330.tag = 0;
    _result_3330.data.err = _M0L6_2aerrS2664;
    return _result_3330;
  }
  #line 33 "/home/developer/Documents/1/alignio.mbt"
  _M0L6_2atmpS2658
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq20MultipleSeqAlignmentE(_M0L10alignmentsS1182);
  if (_M0L6_2atmpS2658 == 0) {
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2659;
    struct moonbit_result_1 _result_3331;
    moonbit_decref(_M0L10alignmentsS1182);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2659
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2659)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2659)->$0
    = (moonbit_string_t)moonbit_string_literal_62.data;
    _result_3331.tag = 0;
    _result_3331.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2659;
    return _result_3331;
  }
  #line 36 "/home/developer/Documents/1/alignio.mbt"
  _M0L6_2atmpS2660
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq20MultipleSeqAlignmentE(_M0L10alignmentsS1182);
  if (_M0L6_2atmpS2660 > 1) {
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2661;
    struct moonbit_result_1 _result_3332;
    moonbit_decref(_M0L10alignmentsS1182);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2661
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2661)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2661)->$0
    = (moonbit_string_t)moonbit_string_literal_63.data;
    _result_3332.tag = 0;
    _result_3332.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2661;
    return _result_3332;
  }
  #line 39 "/home/developer/Documents/1/alignio.mbt"
  _M0L6_2atmpS2662
  = _M0MPC15array5Array2atGRP26biolab8bio__seq20MultipleSeqAlignmentE(_M0L10alignmentsS1182, 0);
  moonbit_decref(_M0L10alignmentsS1182);
  _result_3333.tag = 1;
  _result_3333.data.ok = _M0L6_2atmpS2662;
  return _result_3333;
}

struct moonbit_result_3 _M0FP26biolab8bio__seq14alignio__parse(
  moonbit_string_t _M0L7contentS1181,
  moonbit_string_t _M0L6formatS1180
) {
  #line 10 "/home/developer/Documents/1/alignio.mbt"
  if (
    _M0L6formatS1180 == (moonbit_string_t)moonbit_string_literal_13.data
    || Moonbit_array_length(_M0L6formatS1180) == 7
       && 0
          == memcmp(_M0L6formatS1180, (moonbit_string_t)moonbit_string_literal_13.data, 14)
  ) {
    struct moonbit_result_1 _tmp_3344;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2638;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2637;
    struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2636;
    struct moonbit_result_3 _result_3346;
    #line 15 "/home/developer/Documents/1/alignio.mbt"
    _tmp_3344 = _M0FP26biolab8bio__seq14parse__clustal(_M0L7contentS1181);
    if (_tmp_3344.tag) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2639 =
        _tmp_3344.data.ok;
      _M0L6_2atmpS2638 = _M0L5_2aokS2639;
    } else {
      void* const _M0L6_2aerrS2640 = _tmp_3344.data.err;
      struct moonbit_result_3 _result_3345;
      _result_3345.tag = 0;
      _result_3345.data.err = _M0L6_2aerrS2640;
      return _result_3345;
    }
    _M0L6_2atmpS2637
    = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
    _M0L6_2atmpS2637[0] = _M0L6_2atmpS2638;
    _M0L6_2atmpS2636
    = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
    Moonbit_object_header(_M0L6_2atmpS2636)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
    _M0L6_2atmpS2636->$0 = _M0L6_2atmpS2637;
    _M0L6_2atmpS2636->$1 = 1;
    _result_3346.tag = 1;
    _result_3346.data.ok = _M0L6_2atmpS2636;
    return _result_3346;
  } else if (
           _M0L6formatS1180
           == (moonbit_string_t)moonbit_string_literal_15.data
           || Moonbit_array_length(_M0L6formatS1180) == 6
              && 0
                 == memcmp(_M0L6formatS1180, (moonbit_string_t)moonbit_string_literal_15.data, 12)
         ) {
    struct moonbit_result_1 _tmp_3341;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2643;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2642;
    struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2641;
    struct moonbit_result_3 _result_3343;
    #line 16 "/home/developer/Documents/1/alignio.mbt"
    _tmp_3341 = _M0FP26biolab8bio__seq13parse__phylip(_M0L7contentS1181);
    if (_tmp_3341.tag) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2644 =
        _tmp_3341.data.ok;
      _M0L6_2atmpS2643 = _M0L5_2aokS2644;
    } else {
      void* const _M0L6_2aerrS2645 = _tmp_3341.data.err;
      struct moonbit_result_3 _result_3342;
      _result_3342.tag = 0;
      _result_3342.data.err = _M0L6_2aerrS2645;
      return _result_3342;
    }
    _M0L6_2atmpS2642
    = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
    _M0L6_2atmpS2642[0] = _M0L6_2atmpS2643;
    _M0L6_2atmpS2641
    = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
    Moonbit_object_header(_M0L6_2atmpS2641)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
    _M0L6_2atmpS2641->$0 = _M0L6_2atmpS2642;
    _M0L6_2atmpS2641->$1 = 1;
    _result_3343.tag = 1;
    _result_3343.data.ok = _M0L6_2atmpS2641;
    return _result_3343;
  } else if (
           _M0L6formatS1180
           == (moonbit_string_t)moonbit_string_literal_17.data
           || Moonbit_array_length(_M0L6formatS1180) == 17
              && 0
                 == memcmp(_M0L6formatS1180, (moonbit_string_t)moonbit_string_literal_17.data, 34)
         ) {
    struct moonbit_result_1 _tmp_3338;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2648;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2647;
    struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2646;
    struct moonbit_result_3 _result_3340;
    #line 17 "/home/developer/Documents/1/alignio.mbt"
    _tmp_3338
    = _M0FP26biolab8bio__seq25parse__phylip__sequential(_M0L7contentS1181);
    if (_tmp_3338.tag) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2649 =
        _tmp_3338.data.ok;
      _M0L6_2atmpS2648 = _M0L5_2aokS2649;
    } else {
      void* const _M0L6_2aerrS2650 = _tmp_3338.data.err;
      struct moonbit_result_3 _result_3339;
      _result_3339.tag = 0;
      _result_3339.data.err = _M0L6_2aerrS2650;
      return _result_3339;
    }
    _M0L6_2atmpS2647
    = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
    _M0L6_2atmpS2647[0] = _M0L6_2atmpS2648;
    _M0L6_2atmpS2646
    = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
    Moonbit_object_header(_M0L6_2atmpS2646)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
    _M0L6_2atmpS2646->$0 = _M0L6_2atmpS2647;
    _M0L6_2atmpS2646->$1 = 1;
    _result_3340.tag = 1;
    _result_3340.data.ok = _M0L6_2atmpS2646;
    return _result_3340;
  } else if (
           _M0L6formatS1180
           == (moonbit_string_t)moonbit_string_literal_19.data
           || Moonbit_array_length(_M0L6formatS1180) == 5
              && 0
                 == memcmp(_M0L6formatS1180, (moonbit_string_t)moonbit_string_literal_19.data, 10)
         ) {
    struct moonbit_result_1 _tmp_3335;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2653;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2652;
    struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L6_2atmpS2651;
    struct moonbit_result_3 _result_3337;
    #line 18 "/home/developer/Documents/1/alignio.mbt"
    _tmp_3335
    = _M0FP26biolab8bio__seq27parse__fasta__as__alignment(_M0L7contentS1181);
    if (_tmp_3335.tag) {
      struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* const _M0L5_2aokS2654 =
        _tmp_3335.data.ok;
      _M0L6_2atmpS2653 = _M0L5_2aokS2654;
    } else {
      void* const _M0L6_2aerrS2655 = _tmp_3335.data.err;
      struct moonbit_result_3 _result_3336;
      _result_3336.tag = 0;
      _result_3336.data.err = _M0L6_2aerrS2655;
      return _result_3336;
    }
    _M0L6_2atmpS2652
    = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment**)moonbit_make_ref_array_raw(1);
    _M0L6_2atmpS2652[0] = _M0L6_2atmpS2653;
    _M0L6_2atmpS2651
    = (struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE));
    Moonbit_object_header(_M0L6_2atmpS2651)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 21, 0);
    _M0L6_2atmpS2651->$0 = _M0L6_2atmpS2652;
    _M0L6_2atmpS2651->$1 = 1;
    _result_3337.tag = 1;
    _result_3337.data.ok = _M0L6_2atmpS2651;
    return _result_3337;
  } else {
    moonbit_string_t _M0L6_2atmpS2657;
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2656;
    struct moonbit_result_3 _result_3334;
    #line 19 "/home/developer/Documents/1/alignio.mbt"
    _M0L6_2atmpS2657
    = moonbit_add_string((moonbit_string_t)moonbit_string_literal_61.data, _M0L6formatS1180);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2656
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2656)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2656)->$0
    = _M0L6_2atmpS2657;
    _result_3334.tag = 0;
    _result_3334.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2656;
    return _result_3334;
  }
}

struct moonbit_result_1 _M0FP26biolab8bio__seq27parse__fasta__as__alignment(
  moonbit_string_t _M0L7contentS1179
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1178;
  struct _M0TPB3MapGssE* _M0L6_2atmpS2634;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L6_2atmpS2635;
  struct moonbit_result_1 _result_3347;
  #line 61 "/home/developer/Documents/1/alignio.mbt"
  #line 64 "/home/developer/Documents/1/alignio.mbt"
  _M0L7recordsS1178 = _M0FP26biolab8bio__seq12parse__fasta(_M0L7contentS1179);
  _M0L6_2atmpS2634 = 0;
  _M0L6_2atmpS2635 = 0;
  #line 65 "/home/developer/Documents/1/alignio.mbt"
  _result_3347
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(_M0L7recordsS1178, _M0L6_2atmpS2634, _M0L6_2atmpS2635);
  moonbit_decref(_M0L7recordsS1178);
  if (_M0L6_2atmpS2634) {
    moonbit_decref(_M0L6_2atmpS2634);
  }
  if (_M0L6_2atmpS2635) {
    moonbit_decref(_M0L6_2atmpS2635);
  }
  return _result_3347;
}

moonbit_string_t _M0FP26biolab8bio__seq14write__clustal(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1166
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1164;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1165;
  int32_t _M0L7n__colsS1167;
  int32_t _M0L11line__widthS1168;
  int32_t _M0L11id__paddingS1169;
  struct _M0TPB8MutLocalGiE* _M0L3posS1170;
  moonbit_string_t _result_3350;
  #line 104 "/home/developer/Documents/1/clustal_io.mbt"
  #line 105 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L3bufS1164 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 106 "/home/developer/Documents/1/clustal_io.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1164, (moonbit_string_t)moonbit_string_literal_64.data);
  #line 108 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L7recordsS1165
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(_M0L5alignS1166);
  #line 109 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L7n__colsS1167
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(_M0L5alignS1166);
  _M0L11line__widthS1168 = 60;
  _M0L11id__paddingS1169 = 36;
  _M0L3posS1170
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3posS1170)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3posS1170->$0 = 0;
  while (1) {
    int32_t _M0L3valS2618 = _M0L3posS1170->$0;
    if (_M0L3valS2618 < _M0L7n__colsS1167) {
      int32_t _M0L3valS2632 = _M0L3posS1170->$0;
      int32_t _M0L6_2atmpS2631 = _M0L3valS2632 + _M0L11line__widthS1168;
      int32_t _M0L3endS1171;
      int32_t _M0L7_2abindS1172;
      int32_t _M0L2__S1173;
      if (_M0L6_2atmpS2631 < _M0L7n__colsS1167) {
        int32_t _M0L3valS2633 = _M0L3posS1170->$0;
        _M0L3endS1171 = _M0L3valS2633 + _M0L11line__widthS1168;
      } else {
        _M0L3endS1171 = _M0L7n__colsS1167;
      }
      _M0L7_2abindS1172 = _M0L7recordsS1165->$1;
      _M0L2__S1173 = 0;
      while (1) {
        if (_M0L2__S1173 < _M0L7_2abindS1172) {
          struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2630 =
            _M0L7recordsS1165->$0;
          struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1174 =
            (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2630[
              _M0L2__S1173
            ];
          moonbit_string_t _M0L2idS2619 = _M0L3recS1174->$1;
          moonbit_string_t _M0L2idS2628;
          int32_t _M0L6_2atmpS2627;
          int32_t _M0L7paddingS1175;
          struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2970;
          int32_t _M0L6_2acntS3214;
          struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2626;
          moonbit_string_t _M0L6_2atmpS2623;
          int32_t _M0L3valS2624;
          int64_t _M0L6_2atmpS2625;
          struct _M0TPC16string10StringView _M0L6_2atmpS2622;
          moonbit_string_t _M0L6_2atmpS2621;
          int32_t _M0L6_2atmpS2629;
          moonbit_incref(_M0L2idS2619);
          moonbit_incref(_M0L3recS1174);
          #line 117 "/home/developer/Documents/1/clustal_io.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1164, _M0L2idS2619);
          moonbit_decref(_M0L2idS2619);
          _M0L2idS2628 = _M0L3recS1174->$1;
          _M0L6_2atmpS2627 = Moonbit_array_length(_M0L2idS2628);
          _M0L7paddingS1175 = _M0L11id__paddingS1169 - _M0L6_2atmpS2627;
          if (_M0L7paddingS1175 > 0) {
            moonbit_string_t _M0L6_2atmpS2620;
            #line 120 "/home/developer/Documents/1/clustal_io.mbt"
            _M0L6_2atmpS2620
            = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_65.data, _M0L7paddingS1175);
            #line 120 "/home/developer/Documents/1/clustal_io.mbt"
            _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1164, _M0L6_2atmpS2620);
            moonbit_decref(_M0L6_2atmpS2620);
          } else {
            #line 122 "/home/developer/Documents/1/clustal_io.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1164, 32);
          }
          _M0L8_2afieldS2970 = _M0L3recS1174->$0;
          _M0L6_2acntS3214
          = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1174));
          if (_M0L6_2acntS3214 > 1) {
            int32_t _M0L11_2anew__cntS3219 = _M0L6_2acntS3214 - 1;
            Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1174), _M0L11_2anew__cntS3219);
            moonbit_incref(_M0L8_2afieldS2970);
          } else if (_M0L6_2acntS3214 == 1) {
            struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3218 =
              _M0L3recS1174->$4;
            moonbit_string_t _M0L8_2afieldS3217;
            moonbit_string_t _M0L8_2afieldS3216;
            moonbit_string_t _M0L8_2afieldS3215;
            moonbit_decref(_M0L8_2afieldS3218);
            _M0L8_2afieldS3217 = _M0L3recS1174->$3;
            moonbit_decref(_M0L8_2afieldS3217);
            _M0L8_2afieldS3216 = _M0L3recS1174->$2;
            moonbit_decref(_M0L8_2afieldS3216);
            _M0L8_2afieldS3215 = _M0L3recS1174->$1;
            moonbit_decref(_M0L8_2afieldS3215);
            #line 124 "/home/developer/Documents/1/clustal_io.mbt"
            moonbit_free(_M0L3recS1174);
          }
          _M0L3seqS2626 = _M0L8_2afieldS2970;
          #line 124 "/home/developer/Documents/1/clustal_io.mbt"
          _M0L6_2atmpS2623
          = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2626);
          moonbit_decref(_M0L3seqS2626);
          _M0L3valS2624 = _M0L3posS1170->$0;
          _M0L6_2atmpS2625 = (int64_t)_M0L3endS1171;
          #line 124 "/home/developer/Documents/1/clustal_io.mbt"
          _M0L6_2atmpS2622
          = _M0MPC16string6String11sub_2einner(_M0L6_2atmpS2623, _M0L3valS2624, _M0L6_2atmpS2625);
          moonbit_decref(_M0L6_2atmpS2623);
          #line 124 "/home/developer/Documents/1/clustal_io.mbt"
          _M0L6_2atmpS2621
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2622);
          moonbit_decref(_M0L6_2atmpS2622.$0);
          #line 124 "/home/developer/Documents/1/clustal_io.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1164, _M0L6_2atmpS2621);
          moonbit_decref(_M0L6_2atmpS2621);
          #line 125 "/home/developer/Documents/1/clustal_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1164, 10);
          _M0L6_2atmpS2629 = _M0L2__S1173 + 1;
          _M0L2__S1173 = _M0L6_2atmpS2629;
          continue;
        }
        break;
      }
      #line 127 "/home/developer/Documents/1/clustal_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1164, 10);
      _M0L3posS1170->$0 = _M0L3endS1171;
      continue;
    } else {
      moonbit_decref(_M0L3posS1170);
      moonbit_decref(_M0L7recordsS1165);
    }
    break;
  }
  #line 130 "/home/developer/Documents/1/clustal_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1164, 10);
  #line 131 "/home/developer/Documents/1/clustal_io.mbt"
  _result_3350 = _M0MPB13StringBuilder10to__string(_M0L3bufS1164);
  moonbit_decref(_M0L3bufS1164);
  return _result_3350;
}

struct moonbit_result_1 _M0FP26biolab8bio__seq14parse__clustal(
  moonbit_string_t _M0L7contentS1125
) {
  moonbit_string_t _M0L7_2abindS1126;
  int32_t _M0L6_2atmpS2617;
  struct _M0TPC16string10StringView _M0L6_2atmpS2616;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2615;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS1124;
  int32_t _M0L1nS1127;
  struct _M0TPB8MutLocalGiE* _M0L1iS1128;
  struct _M0TUsRPB13StringBuilderE** _M0L7_2abindS1130;
  struct _M0TUsRPB13StringBuilderE** _M0L6_2atmpS2614;
  struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE _M0L6_2atmpS2613;
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L7recordsS1129;
  int32_t _M0L6_2atmpS2612;
  int32_t _M0L2gtS1131;
  int32_t _M0L6_2atmpS2611;
  int32_t _M0L2crS1132;
  int32_t _M0L6_2atmpS2610;
  int32_t _M0L2spS1133;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L12seq__recordsS1151;
  struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0L5_2aitS1152;
  struct _M0TPB3MapGssE* _M0L6_2atmpS2608;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L6_2atmpS2609;
  struct moonbit_result_1 _result_3362;
  #line 19 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L7_2abindS1126 = (moonbit_string_t)moonbit_string_literal_42.data;
  _M0L6_2atmpS2617 = Moonbit_array_length(_M0L7_2abindS1126);
  _M0L6_2atmpS2616
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1126, .$1 = 0, .$2 = _M0L6_2atmpS2617
  };
  #line 22 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L6_2atmpS2615
  = _M0MPC16string6String5split(_M0L7contentS1125, _M0L6_2atmpS2616);
  moonbit_decref(_M0L6_2atmpS2616.$0);
  #line 22 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L5linesS1124
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2615);
  moonbit_decref(_M0L6_2atmpS2615);
  #line 23 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L1nS1127
  = _M0MPC15array5Array6lengthGRPC16string10StringViewE(_M0L5linesS1124);
  _M0L1iS1128
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1128)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1128->$0 = 0;
  _M0L7_2abindS1130
  = (struct _M0TUsRPB13StringBuilderE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2614 = _M0L7_2abindS1130;
  _M0L6_2atmpS2613
  = (struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE){
    .$0 = _M0L6_2atmpS2614, .$1 = 0, .$2 = 0
  };
  #line 25 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L7recordsS1129
  = _M0MPB3Map3MapGsRPB13StringBuilderE(_M0L6_2atmpS2613, 16ll);
  moonbit_decref(_M0L6_2atmpS2613.$0);
  _M0L6_2atmpS2612 = 62;
  _M0L2gtS1131 = (uint16_t)_M0L6_2atmpS2612;
  _M0L6_2atmpS2611 = 13;
  _M0L2crS1132 = (uint16_t)_M0L6_2atmpS2611;
  _M0L6_2atmpS2610 = 32;
  _M0L2spS1133 = (uint16_t)_M0L6_2atmpS2610;
  while (1) {
    int32_t _M0L3valS2536 = _M0L1iS1128->$0;
    if (_M0L3valS2536 < _M0L1nS1127) {
      int32_t _M0L3valS2602 = _M0L1iS1128->$0;
      struct _M0TPC16string10StringView _M0L4lineS1134;
      int32_t _M0L3lenS1135;
      struct _M0TPB8MutLocalGiE* _M0L7trimmedS1136;
      int32_t _if__result_3352;
      int32_t _M0L3valS2540;
      int32_t _M0L5firstS1138;
      int32_t _M0L6_2atmpS2550;
      int32_t _M0L6_2atmpS2549;
      int32_t _if__result_3353;
      int32_t _M0L3valS2574;
      int32_t _if__result_3354;
      struct _M0TPB8MutLocalGiE* _M0L1jS1139;
      int32_t _M0L3valS2583;
      int32_t _M0L3valS2601;
      int64_t _M0L6_2atmpS2600;
      struct _M0TPC16string10StringView _M0L6_2atmpS2599;
      moonbit_string_t _M0L4nameS1141;
      int32_t _M0L3valS2598;
      struct _M0TPB8MutLocalGiE* _M0L1kS1142;
      int32_t _M0L3valS2595;
      int32_t _M0L3valS2597;
      int64_t _M0L6_2atmpS2596;
      struct _M0TPC16string10StringView _M0L6_2atmpS2594;
      moonbit_string_t _M0L9seq__partS1144;
      struct _M0TPB13StringBuilder* _M0L2sbS1146;
      struct _M0TPB13StringBuilder* _M0L7_2abindS1147;
      int32_t _M0L3valS2593;
      int32_t _M0L6_2atmpS2592;
      #line 31 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L4lineS1134
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1124, _M0L3valS2602);
      #line 32 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L3lenS1135 = _M0MPC16string10StringView6length(_M0L4lineS1134);
      _M0L7trimmedS1136
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L7trimmedS1136)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L7trimmedS1136->$0 = _M0L3lenS1135;
      if (_M0L3lenS1135 > 0) {
        int32_t _M0L6_2atmpS2538 = _M0L3lenS1135 - 1;
        int32_t _M0L6_2atmpS2537;
        #line 34 "/home/developer/Documents/1/clustal_io.mbt"
        _M0L6_2atmpS2537
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, _M0L6_2atmpS2538);
        _if__result_3352 = _M0L6_2atmpS2537 == _M0L2crS1132;
      } else {
        _if__result_3352 = 0;
      }
      if (_if__result_3352) {
        int32_t _M0L6_2atmpS2539 = _M0L3lenS1135 - 1;
        _M0L7trimmedS1136->$0 = _M0L6_2atmpS2539;
      }
      _M0L3valS2540 = _M0L7trimmedS1136->$0;
      if (_M0L3valS2540 == 0) {
        int32_t _M0L3valS2542;
        int32_t _M0L6_2atmpS2541;
        moonbit_decref(_M0L7trimmedS1136);
        moonbit_decref(_M0L4lineS1134.$0);
        _M0L3valS2542 = _M0L1iS1128->$0;
        _M0L6_2atmpS2541 = _M0L3valS2542 + 1;
        _M0L1iS1128->$0 = _M0L6_2atmpS2541;
        continue;
      }
      #line 41 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L5firstS1138
      = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 0);
      if (_M0L5firstS1138 == _M0L2gtS1131) {
        moonbit_decref(_M0L7trimmedS1136);
        moonbit_decref(_M0L4lineS1134.$0);
        moonbit_decref(_M0L1iS1128);
        moonbit_decref(_M0L5linesS1124);
        break;
      }
      _M0L6_2atmpS2550 = 42;
      _M0L6_2atmpS2549 = (uint16_t)_M0L6_2atmpS2550;
      if (_M0L5firstS1138 == _M0L6_2atmpS2549) {
        _if__result_3353 = 1;
      } else {
        int32_t _M0L6_2atmpS2548 = 58;
        int32_t _M0L6_2atmpS2547 = (uint16_t)_M0L6_2atmpS2548;
        if (_M0L5firstS1138 == _M0L6_2atmpS2547) {
          _if__result_3353 = 1;
        } else {
          int32_t _M0L6_2atmpS2546 = 46;
          int32_t _M0L6_2atmpS2545 = (uint16_t)_M0L6_2atmpS2546;
          if (_M0L5firstS1138 == _M0L6_2atmpS2545) {
            _if__result_3353 = 1;
          } else {
            int32_t _M0L6_2atmpS2544 = 32;
            int32_t _M0L6_2atmpS2543 = (uint16_t)_M0L6_2atmpS2544;
            _if__result_3353 = _M0L5firstS1138 == _M0L6_2atmpS2543;
          }
        }
      }
      if (_if__result_3353) {
        int32_t _M0L3valS2552;
        int32_t _M0L6_2atmpS2551;
        moonbit_decref(_M0L7trimmedS1136);
        moonbit_decref(_M0L4lineS1134.$0);
        _M0L3valS2552 = _M0L1iS1128->$0;
        _M0L6_2atmpS2551 = _M0L3valS2552 + 1;
        _M0L1iS1128->$0 = _M0L6_2atmpS2551;
        continue;
      }
      _M0L3valS2574 = _M0L7trimmedS1136->$0;
      if (_M0L3valS2574 >= 7) {
        int32_t _M0L6_2atmpS2571;
        int32_t _M0L6_2atmpS2573;
        int32_t _M0L6_2atmpS2572;
        #line 53 "/home/developer/Documents/1/clustal_io.mbt"
        _M0L6_2atmpS2571
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 0);
        _M0L6_2atmpS2573 = 67;
        _M0L6_2atmpS2572 = (uint16_t)_M0L6_2atmpS2573;
        if (_M0L6_2atmpS2571 == _M0L6_2atmpS2572) {
          int32_t _M0L6_2atmpS2568;
          int32_t _M0L6_2atmpS2570;
          int32_t _M0L6_2atmpS2569;
          #line 54 "/home/developer/Documents/1/clustal_io.mbt"
          _M0L6_2atmpS2568
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 1);
          _M0L6_2atmpS2570 = 76;
          _M0L6_2atmpS2569 = (uint16_t)_M0L6_2atmpS2570;
          if (_M0L6_2atmpS2568 == _M0L6_2atmpS2569) {
            int32_t _M0L6_2atmpS2565;
            int32_t _M0L6_2atmpS2567;
            int32_t _M0L6_2atmpS2566;
            #line 55 "/home/developer/Documents/1/clustal_io.mbt"
            _M0L6_2atmpS2565
            = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 2);
            _M0L6_2atmpS2567 = 85;
            _M0L6_2atmpS2566 = (uint16_t)_M0L6_2atmpS2567;
            if (_M0L6_2atmpS2565 == _M0L6_2atmpS2566) {
              int32_t _M0L6_2atmpS2562;
              int32_t _M0L6_2atmpS2564;
              int32_t _M0L6_2atmpS2563;
              #line 56 "/home/developer/Documents/1/clustal_io.mbt"
              _M0L6_2atmpS2562
              = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 3);
              _M0L6_2atmpS2564 = 83;
              _M0L6_2atmpS2563 = (uint16_t)_M0L6_2atmpS2564;
              if (_M0L6_2atmpS2562 == _M0L6_2atmpS2563) {
                int32_t _M0L6_2atmpS2559;
                int32_t _M0L6_2atmpS2561;
                int32_t _M0L6_2atmpS2560;
                #line 57 "/home/developer/Documents/1/clustal_io.mbt"
                _M0L6_2atmpS2559
                = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 4);
                _M0L6_2atmpS2561 = 84;
                _M0L6_2atmpS2560 = (uint16_t)_M0L6_2atmpS2561;
                if (_M0L6_2atmpS2559 == _M0L6_2atmpS2560) {
                  int32_t _M0L6_2atmpS2556;
                  int32_t _M0L6_2atmpS2558;
                  int32_t _M0L6_2atmpS2557;
                  #line 58 "/home/developer/Documents/1/clustal_io.mbt"
                  _M0L6_2atmpS2556
                  = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 5);
                  _M0L6_2atmpS2558 = 65;
                  _M0L6_2atmpS2557 = (uint16_t)_M0L6_2atmpS2558;
                  if (_M0L6_2atmpS2556 == _M0L6_2atmpS2557) {
                    int32_t _M0L6_2atmpS2553;
                    int32_t _M0L6_2atmpS2555;
                    int32_t _M0L6_2atmpS2554;
                    #line 59 "/home/developer/Documents/1/clustal_io.mbt"
                    _M0L6_2atmpS2553
                    = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, 6);
                    _M0L6_2atmpS2555 = 76;
                    _M0L6_2atmpS2554 = (uint16_t)_M0L6_2atmpS2555;
                    _if__result_3354 = _M0L6_2atmpS2553 == _M0L6_2atmpS2554;
                  } else {
                    _if__result_3354 = 0;
                  }
                } else {
                  _if__result_3354 = 0;
                }
              } else {
                _if__result_3354 = 0;
              }
            } else {
              _if__result_3354 = 0;
            }
          } else {
            _if__result_3354 = 0;
          }
        } else {
          _if__result_3354 = 0;
        }
      } else {
        _if__result_3354 = 0;
      }
      if (_if__result_3354) {
        int32_t _M0L3valS2576;
        int32_t _M0L6_2atmpS2575;
        moonbit_decref(_M0L7trimmedS1136);
        moonbit_decref(_M0L4lineS1134.$0);
        _M0L3valS2576 = _M0L1iS1128->$0;
        _M0L6_2atmpS2575 = _M0L3valS2576 + 1;
        _M0L1iS1128->$0 = _M0L6_2atmpS2575;
        continue;
      }
      _M0L1jS1139
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1jS1139)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1jS1139->$0 = 0;
      while (1) {
        int32_t _M0L3valS2579 = _M0L1jS1139->$0;
        int32_t _M0L3valS2580 = _M0L7trimmedS1136->$0;
        int32_t _if__result_3356;
        if (_M0L3valS2579 < _M0L3valS2580) {
          int32_t _M0L3valS2578 = _M0L1jS1139->$0;
          int32_t _M0L6_2atmpS2577;
          #line 64 "/home/developer/Documents/1/clustal_io.mbt"
          _M0L6_2atmpS2577
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, _M0L3valS2578);
          _if__result_3356 = _M0L6_2atmpS2577 != _M0L2spS1133;
        } else {
          _if__result_3356 = 0;
        }
        if (_if__result_3356) {
          int32_t _M0L3valS2582 = _M0L1jS1139->$0;
          int32_t _M0L6_2atmpS2581 = _M0L3valS2582 + 1;
          _M0L1jS1139->$0 = _M0L6_2atmpS2581;
          continue;
        }
        break;
      }
      _M0L3valS2583 = _M0L1jS1139->$0;
      if (_M0L3valS2583 == 0) {
        int32_t _M0L3valS2585;
        int32_t _M0L6_2atmpS2584;
        moonbit_decref(_M0L1jS1139);
        moonbit_decref(_M0L7trimmedS1136);
        moonbit_decref(_M0L4lineS1134.$0);
        _M0L3valS2585 = _M0L1iS1128->$0;
        _M0L6_2atmpS2584 = _M0L3valS2585 + 1;
        _M0L1iS1128->$0 = _M0L6_2atmpS2584;
        continue;
      }
      _M0L3valS2601 = _M0L1jS1139->$0;
      _M0L6_2atmpS2600 = (int64_t)_M0L3valS2601;
      #line 71 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L6_2atmpS2599
      = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1134, 0, _M0L6_2atmpS2600);
      #line 71 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L4nameS1141 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2599);
      moonbit_decref(_M0L6_2atmpS2599.$0);
      _M0L3valS2598 = _M0L1jS1139->$0;
      moonbit_decref(_M0L1jS1139);
      _M0L1kS1142
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1kS1142)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1kS1142->$0 = _M0L3valS2598;
      while (1) {
        int32_t _M0L3valS2588 = _M0L1kS1142->$0;
        int32_t _M0L3valS2589 = _M0L7trimmedS1136->$0;
        int32_t _if__result_3358;
        if (_M0L3valS2588 < _M0L3valS2589) {
          int32_t _M0L3valS2587 = _M0L1kS1142->$0;
          int32_t _M0L6_2atmpS2586;
          #line 73 "/home/developer/Documents/1/clustal_io.mbt"
          _M0L6_2atmpS2586
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1134, _M0L3valS2587);
          _if__result_3358 = _M0L6_2atmpS2586 == _M0L2spS1133;
        } else {
          _if__result_3358 = 0;
        }
        if (_if__result_3358) {
          int32_t _M0L3valS2591 = _M0L1kS1142->$0;
          int32_t _M0L6_2atmpS2590 = _M0L3valS2591 + 1;
          _M0L1kS1142->$0 = _M0L6_2atmpS2590;
          continue;
        }
        break;
      }
      _M0L3valS2595 = _M0L1kS1142->$0;
      moonbit_decref(_M0L1kS1142);
      _M0L3valS2597 = _M0L7trimmedS1136->$0;
      moonbit_decref(_M0L7trimmedS1136);
      _M0L6_2atmpS2596 = (int64_t)_M0L3valS2597;
      #line 76 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L6_2atmpS2594
      = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1134, _M0L3valS2595, _M0L6_2atmpS2596);
      moonbit_decref(_M0L4lineS1134.$0);
      #line 76 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L9seq__partS1144
      = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2594);
      moonbit_decref(_M0L6_2atmpS2594.$0);
      #line 77 "/home/developer/Documents/1/clustal_io.mbt"
      _M0L7_2abindS1147
      = _M0MPB3Map3getGsRPB13StringBuilderE(_M0L7recordsS1129, _M0L4nameS1141);
      if (_M0L7_2abindS1147 == 0) {
        struct _M0TPB13StringBuilder* _M0L2sbS1150;
        if (_M0L7_2abindS1147) {
          moonbit_decref(_M0L7_2abindS1147);
        }
        #line 80 "/home/developer/Documents/1/clustal_io.mbt"
        _M0L2sbS1150 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
        #line 81 "/home/developer/Documents/1/clustal_io.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L2sbS1150, _M0L9seq__partS1144);
        moonbit_decref(_M0L9seq__partS1144);
        #line 82 "/home/developer/Documents/1/clustal_io.mbt"
        _M0MPB3Map3setGsRPB13StringBuilderE(_M0L7recordsS1129, _M0L4nameS1141, _M0L2sbS1150);
        moonbit_decref(_M0L4nameS1141);
        moonbit_decref(_M0L2sbS1150);
      } else {
        struct _M0TPB13StringBuilder* _M0L7_2aSomeS1148;
        struct _M0TPB13StringBuilder* _M0L5_2asbS1149;
        moonbit_decref(_M0L4nameS1141);
        _M0L7_2aSomeS1148 = _M0L7_2abindS1147;
        _M0L5_2asbS1149 = _M0L7_2aSomeS1148;
        _M0L2sbS1146 = _M0L5_2asbS1149;
        goto join_1145;
      }
      goto joinlet_3359;
      join_1145:;
      #line 78 "/home/developer/Documents/1/clustal_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L2sbS1146, _M0L9seq__partS1144);
      moonbit_decref(_M0L2sbS1146);
      moonbit_decref(_M0L9seq__partS1144);
      joinlet_3359:;
      _M0L3valS2593 = _M0L1iS1128->$0;
      _M0L6_2atmpS2592 = _M0L3valS2593 + 1;
      _M0L1iS1128->$0 = _M0L6_2atmpS2592;
      continue;
    } else {
      moonbit_decref(_M0L1iS1128);
      moonbit_decref(_M0L5linesS1124);
    }
    break;
  }
  #line 88 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L12seq__recordsS1151
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  #line 88 "/home/developer/Documents/1/clustal_io.mbt"
  _M0L5_2aitS1152 = _M0MPB3Map5iter2GsRPB13StringBuilderE(_M0L7recordsS1129);
  moonbit_decref(_M0L7recordsS1129);
  while (1) {
    moonbit_string_t _M0L4nameS1154;
    struct _M0TPB13StringBuilder* _M0L2sbS1155;
    struct _M0TUsRPB13StringBuilderE* _M0L7_2abindS1159;
    moonbit_string_t _M0L6_2atmpS2607;
    struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS1156;
    struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1157;
    struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2606;
    struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2605;
    struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2604;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2603;
    #line 89 "/home/developer/Documents/1/clustal_io.mbt"
    _M0L7_2abindS1159
    = _M0MPB5Iter24nextGsRPB13StringBuilderE(_M0L5_2aitS1152);
    if (_M0L7_2abindS1159 == 0) {
      if (_M0L7_2abindS1159) {
        moonbit_decref(_M0L7_2abindS1159);
      }
      moonbit_decref(_M0L5_2aitS1152);
    } else {
      struct _M0TUsRPB13StringBuilderE* _M0L7_2aSomeS1160 = _M0L7_2abindS1159;
      struct _M0TUsRPB13StringBuilderE* _M0L4_2axS1161 = _M0L7_2aSomeS1160;
      moonbit_string_t _M0L7_2anameS1162 = _M0L4_2axS1161->$0;
      struct _M0TPB13StringBuilder* _M0L8_2afieldS2975 = _M0L4_2axS1161->$1;
      int32_t _M0L6_2acntS3220 =
        Moonbit_rc_count(Moonbit_object_header(_M0L4_2axS1161));
      struct _M0TPB13StringBuilder* _M0L5_2asbS1163;
      if (_M0L6_2acntS3220 > 1) {
        int32_t _M0L11_2anew__cntS3221 = _M0L6_2acntS3220 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L4_2axS1161), _M0L11_2anew__cntS3221);
        moonbit_incref(_M0L8_2afieldS2975);
        moonbit_incref(_M0L7_2anameS1162);
      } else if (_M0L6_2acntS3220 == 1) {
        #line 89 "/home/developer/Documents/1/clustal_io.mbt"
        moonbit_free(_M0L4_2axS1161);
      }
      _M0L5_2asbS1163 = _M0L8_2afieldS2975;
      _M0L4nameS1154 = _M0L7_2anameS1162;
      _M0L2sbS1155 = _M0L5_2asbS1163;
      goto join_1153;
    }
    goto joinlet_3361;
    join_1153:;
    #line 90 "/home/developer/Documents/1/clustal_io.mbt"
    _M0L6_2atmpS2607 = _M0MPB13StringBuilder10to__string(_M0L2sbS1155);
    moonbit_decref(_M0L2sbS1155);
    #line 90 "/home/developer/Documents/1/clustal_io.mbt"
    _M0L3seqS1156 = _M0MP26biolab8bio__seq3Seq3new(_M0L6_2atmpS2607);
    moonbit_decref(_M0L6_2atmpS2607);
    _M0L7_2abindS1157 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
    _M0L6_2atmpS2606 = _M0L7_2abindS1157;
    _M0L6_2atmpS2605
    = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
      .$0 = _M0L6_2atmpS2606, .$1 = 0, .$2 = 0
    };
    #line 96 "/home/developer/Documents/1/clustal_io.mbt"
    _M0L6_2atmpS2604 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2605, 0ll);
    moonbit_decref(_M0L6_2atmpS2605.$0);
    moonbit_incref(_M0L4nameS1154);
    moonbit_incref(_M0L4nameS1154);
    _M0L6_2atmpS2603
    = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
    Moonbit_object_header(_M0L6_2atmpS2603)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
    _M0L6_2atmpS2603->$0 = _M0L3seqS1156;
    _M0L6_2atmpS2603->$1 = _M0L4nameS1154;
    _M0L6_2atmpS2603->$2 = _M0L4nameS1154;
    _M0L6_2atmpS2603->$3 = _M0L4nameS1154;
    _M0L6_2atmpS2603->$4 = _M0L6_2atmpS2604;
    #line 91 "/home/developer/Documents/1/clustal_io.mbt"
    _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L12seq__recordsS1151, _M0L6_2atmpS2603);
    moonbit_decref(_M0L6_2atmpS2603);
    continue;
    joinlet_3361:;
    break;
  }
  _M0L6_2atmpS2608 = 0;
  _M0L6_2atmpS2609 = 0;
  #line 99 "/home/developer/Documents/1/clustal_io.mbt"
  _result_3362
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(_M0L12seq__recordsS1151, _M0L6_2atmpS2608, _M0L6_2atmpS2609);
  moonbit_decref(_M0L12seq__recordsS1151);
  if (_M0L6_2atmpS2608) {
    moonbit_decref(_M0L6_2atmpS2608);
  }
  if (_M0L6_2atmpS2609) {
    moonbit_decref(_M0L6_2atmpS2609);
  }
  return _result_3362;
}

moonbit_string_t _M0FP26biolab8bio__seq25write__phylip__sequential(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1115
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1113;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1114;
  int32_t _M0L7n__seqsS1116;
  int32_t _M0L7n__colsS1117;
  moonbit_string_t _M0L6_2atmpS2520;
  moonbit_string_t _M0L6_2atmpS2521;
  int32_t _M0L11name__widthS1118;
  int32_t _M0L7_2abindS1119;
  int32_t _M0L2__S1120;
  moonbit_string_t _result_3364;
  #line 227 "/home/developer/Documents/1/phylip_io.mbt"
  #line 228 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L3bufS1113 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 229 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7recordsS1114
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(_M0L5alignS1115);
  #line 230 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7n__seqsS1116
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1114);
  #line 231 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7n__colsS1117
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(_M0L5alignS1115);
  #line 233 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2520
  = _M0MPC13int3Int18to__string_2einner(_M0L7n__seqsS1116, 10);
  #line 233 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1113, _M0L6_2atmpS2520);
  moonbit_decref(_M0L6_2atmpS2520);
  #line 234 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1113, 32);
  #line 235 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2521
  = _M0MPC13int3Int18to__string_2einner(_M0L7n__colsS1117, 10);
  #line 235 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1113, _M0L6_2atmpS2521);
  moonbit_decref(_M0L6_2atmpS2521);
  #line 236 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1113, 10);
  _M0L11name__widthS1118 = 10;
  _M0L7_2abindS1119 = _M0L7recordsS1114->$1;
  _M0L2__S1120 = 0;
  while (1) {
    if (_M0L2__S1120 < _M0L7_2abindS1119) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2535 =
        _M0L7recordsS1114->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1121 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2535[_M0L2__S1120];
      moonbit_string_t _M0L2idS2525 = _M0L3recS1121->$1;
      int32_t _M0L6_2atmpS2524 = Moonbit_array_length(_M0L2idS2525);
      moonbit_string_t _M0L4nameS1122;
      struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2977;
      int32_t _M0L6_2acntS3222;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2523;
      moonbit_string_t _M0L6_2atmpS2522;
      int32_t _M0L6_2atmpS2534;
      if (_M0L6_2atmpS2524 >= _M0L11name__widthS1118) {
        moonbit_string_t _M0L2idS2527 = _M0L3recS1121->$1;
        int64_t _M0L6_2atmpS2528 = (int64_t)_M0L11name__widthS1118;
        struct _M0TPC16string10StringView _M0L6_2atmpS2526;
        moonbit_incref(_M0L2idS2527);
        moonbit_incref(_M0L3recS1121);
        #line 241 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2526
        = _M0MPC16string6String11sub_2einner(_M0L2idS2527, 0, _M0L6_2atmpS2528);
        moonbit_decref(_M0L2idS2527);
        #line 241 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L4nameS1122
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2526);
        moonbit_decref(_M0L6_2atmpS2526.$0);
      } else {
        moonbit_string_t _M0L2idS2529 = _M0L3recS1121->$1;
        moonbit_string_t _M0L2idS2533 = _M0L3recS1121->$1;
        int32_t _M0L6_2atmpS2532 = Moonbit_array_length(_M0L2idS2533);
        int32_t _M0L6_2atmpS2531 = _M0L11name__widthS1118 - _M0L6_2atmpS2532;
        moonbit_string_t _M0L6_2atmpS2530;
        moonbit_incref(_M0L2idS2529);
        moonbit_incref(_M0L3recS1121);
        #line 243 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2530
        = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_65.data, _M0L6_2atmpS2531);
        #line 243 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L4nameS1122 = moonbit_add_string(_M0L2idS2529, _M0L6_2atmpS2530);
        moonbit_decref(_M0L6_2atmpS2530);
        moonbit_decref(_M0L2idS2529);
      }
      #line 245 "/home/developer/Documents/1/phylip_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1113, _M0L4nameS1122);
      moonbit_decref(_M0L4nameS1122);
      #line 246 "/home/developer/Documents/1/phylip_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1113, 32);
      _M0L8_2afieldS2977 = _M0L3recS1121->$0;
      _M0L6_2acntS3222
      = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1121));
      if (_M0L6_2acntS3222 > 1) {
        int32_t _M0L11_2anew__cntS3227 = _M0L6_2acntS3222 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1121), _M0L11_2anew__cntS3227);
        moonbit_incref(_M0L8_2afieldS2977);
      } else if (_M0L6_2acntS3222 == 1) {
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3226 =
          _M0L3recS1121->$4;
        moonbit_string_t _M0L8_2afieldS3225;
        moonbit_string_t _M0L8_2afieldS3224;
        moonbit_string_t _M0L8_2afieldS3223;
        moonbit_decref(_M0L8_2afieldS3226);
        _M0L8_2afieldS3225 = _M0L3recS1121->$3;
        moonbit_decref(_M0L8_2afieldS3225);
        _M0L8_2afieldS3224 = _M0L3recS1121->$2;
        moonbit_decref(_M0L8_2afieldS3224);
        _M0L8_2afieldS3223 = _M0L3recS1121->$1;
        moonbit_decref(_M0L8_2afieldS3223);
        #line 247 "/home/developer/Documents/1/phylip_io.mbt"
        moonbit_free(_M0L3recS1121);
      }
      _M0L3seqS2523 = _M0L8_2afieldS2977;
      #line 247 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2522
      = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2523);
      moonbit_decref(_M0L3seqS2523);
      #line 247 "/home/developer/Documents/1/phylip_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1113, _M0L6_2atmpS2522);
      moonbit_decref(_M0L6_2atmpS2522);
      #line 248 "/home/developer/Documents/1/phylip_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1113, 10);
      _M0L6_2atmpS2534 = _M0L2__S1120 + 1;
      _M0L2__S1120 = _M0L6_2atmpS2534;
      continue;
    } else {
      moonbit_decref(_M0L7recordsS1114);
    }
    break;
  }
  #line 250 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1113, 10);
  #line 251 "/home/developer/Documents/1/phylip_io.mbt"
  _result_3364 = _M0MPB13StringBuilder10to__string(_M0L3bufS1113);
  moonbit_decref(_M0L3bufS1113);
  return _result_3364;
}

moonbit_string_t _M0FP26biolab8bio__seq13write__phylip(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L5alignS1095
) {
  struct _M0TPB13StringBuilder* _M0L3bufS1093;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1094;
  int32_t _M0L7n__seqsS1096;
  int32_t _M0L7n__colsS1097;
  moonbit_string_t _M0L6_2atmpS2490;
  moonbit_string_t _M0L6_2atmpS2491;
  int32_t _M0L11line__widthS1098;
  int32_t _M0L11name__widthS1099;
  int32_t _M0L11group__sizeS1100;
  struct _M0TPB8MutLocalGiE* _M0L3posS1101;
  moonbit_string_t _result_3368;
  #line 177 "/home/developer/Documents/1/phylip_io.mbt"
  #line 178 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L3bufS1093 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 179 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7recordsS1094
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(_M0L5alignS1095);
  #line 180 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7n__seqsS1096
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1094);
  #line 181 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7n__colsS1097
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(_M0L5alignS1095);
  #line 183 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 32);
  #line 184 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2490
  = _M0MPC13int3Int18to__string_2einner(_M0L7n__seqsS1096, 10);
  #line 184 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1093, _M0L6_2atmpS2490);
  moonbit_decref(_M0L6_2atmpS2490);
  #line 185 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 32);
  #line 186 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2491
  = _M0MPC13int3Int18to__string_2einner(_M0L7n__colsS1097, 10);
  #line 186 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1093, _M0L6_2atmpS2491);
  moonbit_decref(_M0L6_2atmpS2491);
  #line 187 "/home/developer/Documents/1/phylip_io.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 10);
  _M0L11line__widthS1098 = 60;
  _M0L11name__widthS1099 = 10;
  _M0L11group__sizeS1100 = 10;
  _M0L3posS1101
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3posS1101)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3posS1101->$0 = 0;
  while (1) {
    int32_t _M0L3valS2492 = _M0L3posS1101->$0;
    if (_M0L3valS2492 < _M0L7n__colsS1097) {
      int32_t _M0L3valS2518 = _M0L3posS1101->$0;
      int32_t _M0L6_2atmpS2517 = _M0L3valS2518 + _M0L11line__widthS1098;
      int32_t _M0L3endS1102;
      int32_t _M0L7_2abindS1103;
      int32_t _M0L2__S1104;
      int32_t _M0L3valS2516;
      int32_t _M0L6_2atmpS2515;
      if (_M0L6_2atmpS2517 < _M0L7n__colsS1097) {
        int32_t _M0L3valS2519 = _M0L3posS1101->$0;
        _M0L3endS1102 = _M0L3valS2519 + _M0L11line__widthS1098;
      } else {
        _M0L3endS1102 = _M0L7n__colsS1097;
      }
      _M0L7_2abindS1103 = _M0L7recordsS1094->$1;
      _M0L2__S1104 = 0;
      while (1) {
        if (_M0L2__S1104 < _M0L7_2abindS1103) {
          struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2514 =
            _M0L7recordsS1094->$0;
          struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1105 =
            (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2514[
              _M0L2__S1104
            ];
          moonbit_string_t _M0L2idS2504 = _M0L3recS1105->$1;
          int32_t _M0L6_2atmpS2503 = Moonbit_array_length(_M0L2idS2504);
          moonbit_string_t _M0L4nameS1106;
          struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2984;
          int32_t _M0L6_2acntS3228;
          struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2502;
          moonbit_string_t _M0L8seq__strS1107;
          int32_t _M0L3valS2501;
          struct _M0TPB8MutLocalGiE* _M0L1iS1108;
          int32_t _M0L6_2atmpS2513;
          if (_M0L6_2atmpS2503 >= _M0L11name__widthS1099) {
            moonbit_string_t _M0L2idS2506 = _M0L3recS1105->$1;
            int64_t _M0L6_2atmpS2507 = (int64_t)_M0L11name__widthS1099;
            struct _M0TPC16string10StringView _M0L6_2atmpS2505;
            moonbit_incref(_M0L2idS2506);
            moonbit_incref(_M0L3recS1105);
            #line 198 "/home/developer/Documents/1/phylip_io.mbt"
            _M0L6_2atmpS2505
            = _M0MPC16string6String11sub_2einner(_M0L2idS2506, 0, _M0L6_2atmpS2507);
            moonbit_decref(_M0L2idS2506);
            #line 198 "/home/developer/Documents/1/phylip_io.mbt"
            _M0L4nameS1106
            = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2505);
            moonbit_decref(_M0L6_2atmpS2505.$0);
          } else {
            moonbit_string_t _M0L2idS2508 = _M0L3recS1105->$1;
            moonbit_string_t _M0L2idS2512 = _M0L3recS1105->$1;
            int32_t _M0L6_2atmpS2511 = Moonbit_array_length(_M0L2idS2512);
            int32_t _M0L6_2atmpS2510 =
              _M0L11name__widthS1099 - _M0L6_2atmpS2511;
            moonbit_string_t _M0L6_2atmpS2509;
            moonbit_incref(_M0L2idS2508);
            moonbit_incref(_M0L3recS1105);
            #line 200 "/home/developer/Documents/1/phylip_io.mbt"
            _M0L6_2atmpS2509
            = _M0MPC16string6String6repeat((moonbit_string_t)moonbit_string_literal_65.data, _M0L6_2atmpS2510);
            #line 200 "/home/developer/Documents/1/phylip_io.mbt"
            _M0L4nameS1106
            = moonbit_add_string(_M0L2idS2508, _M0L6_2atmpS2509);
            moonbit_decref(_M0L6_2atmpS2509);
            moonbit_decref(_M0L2idS2508);
          }
          #line 202 "/home/developer/Documents/1/phylip_io.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1093, _M0L4nameS1106);
          moonbit_decref(_M0L4nameS1106);
          #line 203 "/home/developer/Documents/1/phylip_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 32);
          _M0L8_2afieldS2984 = _M0L3recS1105->$0;
          _M0L6_2acntS3228
          = Moonbit_rc_count(Moonbit_object_header(_M0L3recS1105));
          if (_M0L6_2acntS3228 > 1) {
            int32_t _M0L11_2anew__cntS3233 = _M0L6_2acntS3228 - 1;
            Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS1105), _M0L11_2anew__cntS3233);
            moonbit_incref(_M0L8_2afieldS2984);
          } else if (_M0L6_2acntS3228 == 1) {
            struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3232 =
              _M0L3recS1105->$4;
            moonbit_string_t _M0L8_2afieldS3231;
            moonbit_string_t _M0L8_2afieldS3230;
            moonbit_string_t _M0L8_2afieldS3229;
            moonbit_decref(_M0L8_2afieldS3232);
            _M0L8_2afieldS3231 = _M0L3recS1105->$3;
            moonbit_decref(_M0L8_2afieldS3231);
            _M0L8_2afieldS3230 = _M0L3recS1105->$2;
            moonbit_decref(_M0L8_2afieldS3230);
            _M0L8_2afieldS3229 = _M0L3recS1105->$1;
            moonbit_decref(_M0L8_2afieldS3229);
            #line 204 "/home/developer/Documents/1/phylip_io.mbt"
            moonbit_free(_M0L3recS1105);
          }
          _M0L3seqS2502 = _M0L8_2afieldS2984;
          #line 204 "/home/developer/Documents/1/phylip_io.mbt"
          _M0L8seq__strS1107
          = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2502);
          moonbit_decref(_M0L3seqS2502);
          _M0L3valS2501 = _M0L3posS1101->$0;
          _M0L1iS1108
          = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
          Moonbit_object_header(_M0L1iS1108)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
          _M0L1iS1108->$0 = _M0L3valS2501;
          while (1) {
            int32_t _M0L3valS2493 = _M0L1iS1108->$0;
            if (_M0L3valS2493 < _M0L3endS1102) {
              int32_t _M0L3valS2499 = _M0L1iS1108->$0;
              int32_t _M0L6_2atmpS2498 =
                _M0L3valS2499 + _M0L11group__sizeS1100;
              int32_t _M0L4gendS1109;
              int32_t _M0L3valS2496;
              int64_t _M0L6_2atmpS2497;
              struct _M0TPC16string10StringView _M0L6_2atmpS2495;
              moonbit_string_t _M0L6_2atmpS2494;
              if (_M0L6_2atmpS2498 < _M0L3endS1102) {
                int32_t _M0L3valS2500 = _M0L1iS1108->$0;
                _M0L4gendS1109 = _M0L3valS2500 + _M0L11group__sizeS1100;
              } else {
                _M0L4gendS1109 = _M0L3endS1102;
              }
              _M0L3valS2496 = _M0L1iS1108->$0;
              _M0L6_2atmpS2497 = (int64_t)_M0L4gendS1109;
              #line 208 "/home/developer/Documents/1/phylip_io.mbt"
              _M0L6_2atmpS2495
              = _M0MPC16string6String11sub_2einner(_M0L8seq__strS1107, _M0L3valS2496, _M0L6_2atmpS2497);
              #line 208 "/home/developer/Documents/1/phylip_io.mbt"
              _M0L6_2atmpS2494
              = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2495);
              moonbit_decref(_M0L6_2atmpS2495.$0);
              #line 208 "/home/developer/Documents/1/phylip_io.mbt"
              _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS1093, _M0L6_2atmpS2494);
              moonbit_decref(_M0L6_2atmpS2494);
              if (_M0L4gendS1109 < _M0L3endS1102) {
                #line 210 "/home/developer/Documents/1/phylip_io.mbt"
                _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 32);
              }
              _M0L1iS1108->$0 = _M0L4gendS1109;
              continue;
            } else {
              moonbit_decref(_M0L1iS1108);
              moonbit_decref(_M0L8seq__strS1107);
            }
            break;
          }
          #line 214 "/home/developer/Documents/1/phylip_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 32);
          #line 215 "/home/developer/Documents/1/phylip_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 10);
          _M0L6_2atmpS2513 = _M0L2__S1104 + 1;
          _M0L2__S1104 = _M0L6_2atmpS2513;
          continue;
        }
        break;
      }
      _M0L3valS2516 = _M0L3posS1101->$0;
      _M0L6_2atmpS2515 = _M0L3valS2516 + _M0L11line__widthS1098;
      if (_M0L6_2atmpS2515 < _M0L7n__colsS1097) {
        #line 218 "/home/developer/Documents/1/phylip_io.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS1093, 10);
      }
      _M0L3posS1101->$0 = _M0L3endS1102;
      continue;
    } else {
      moonbit_decref(_M0L3posS1101);
      moonbit_decref(_M0L7recordsS1094);
    }
    break;
  }
  #line 222 "/home/developer/Documents/1/phylip_io.mbt"
  _result_3368 = _M0MPB13StringBuilder10to__string(_M0L3bufS1093);
  moonbit_decref(_M0L3bufS1093);
  return _result_3368;
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MP26biolab8bio__seq20MultipleSeqAlignment4iter(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L4selfS1092
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L8_2afieldS2991;
  #line 208 "/home/developer/Documents/1/align.mbt"
  _M0L8_2afieldS2991 = _M0L4selfS1092->$0;
  moonbit_incref(_M0L8_2afieldS2991);
  return _M0L8_2afieldS2991;
}

int32_t _M0MP26biolab8bio__seq20MultipleSeqAlignment22get__alignment__length(
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L4selfS1091
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS2487;
  int32_t _M0L6_2atmpS2486;
  #line 56 "/home/developer/Documents/1/align.mbt"
  _M0L7recordsS2487 = _M0L4selfS1091->$0;
  moonbit_incref(_M0L7recordsS2487);
  #line 59 "/home/developer/Documents/1/align.mbt"
  _M0L6_2atmpS2486
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS2487);
  moonbit_decref(_M0L7recordsS2487);
  if (_M0L6_2atmpS2486 == 0) {
    return 0;
  } else {
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS2489 =
      _M0L4selfS1091->$0;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2488;
    int32_t _result_3369;
    moonbit_incref(_M0L7recordsS2489);
    #line 62 "/home/developer/Documents/1/align.mbt"
    _M0L6_2atmpS2488
    = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS2489, 0);
    moonbit_decref(_M0L7recordsS2489);
    #line 62 "/home/developer/Documents/1/align.mbt"
    _result_3369 = _M0MP26biolab8bio__seq9SeqRecord6length(_M0L6_2atmpS2488);
    moonbit_decref(_M0L6_2atmpS2488);
    return _result_3369;
  }
}

struct moonbit_result_1 _M0FP26biolab8bio__seq25parse__phylip__sequential(
  moonbit_string_t _M0L7contentS1090
) {
  #line 49 "/home/developer/Documents/1/phylip_io.mbt"
  #line 52 "/home/developer/Documents/1/phylip_io.mbt"
  return _M0FP26biolab8bio__seq22parse__phylip__generic(_M0L7contentS1090, 1);
}

struct moonbit_result_1 _M0FP26biolab8bio__seq13parse__phylip(
  moonbit_string_t _M0L7contentS1089
) {
  #line 41 "/home/developer/Documents/1/phylip_io.mbt"
  #line 44 "/home/developer/Documents/1/phylip_io.mbt"
  return _M0FP26biolab8bio__seq22parse__phylip__generic(_M0L7contentS1089, 0);
}

struct moonbit_result_1 _M0FP26biolab8bio__seq22parse__phylip__generic(
  moonbit_string_t _M0L7contentS1024,
  int32_t _M0L10sequentialS1076
) {
  moonbit_string_t _M0L7_2abindS1025;
  int32_t _M0L6_2atmpS2485;
  struct _M0TPC16string10StringView _M0L6_2atmpS2484;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2483;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS1023;
  int32_t _M0L1nS1026;
  struct _M0TPB8MutLocalGiE* _M0L1iS1027;
  int32_t _M0L6_2atmpS2482;
  int32_t _M0L2crS1028;
  int32_t _M0L6_2atmpS2481;
  int32_t _M0L2spS1029;
  int32_t _M0L3valS2398;
  int32_t _M0L3valS2480;
  struct _M0TPC16string10StringView _M0L6headerS1034;
  int32_t _M0L4hlenS1035;
  struct _M0TPB8MutLocalGiE* _M0L8htrimmedS1036;
  int32_t _if__result_3373;
  struct _M0TPB8MutLocalGiE* _M0L10space__posS1037;
  int32_t _M0L3valS2479;
  int64_t _M0L6_2atmpS2478;
  struct _M0TPC16string10StringView _M0L6_2atmpS2477;
  moonbit_string_t _M0L14num__seqs__strS1039;
  int32_t _M0L1nS1042;
  int32_t _M0L9num__seqsS1040;
  int64_t _M0L7_2abindS1043;
  int32_t _M0L3valS2475;
  struct _M0TPB8MutLocalGiE* _M0L8col__posS1046;
  int32_t _M0L3valS2472;
  int32_t _M0L3valS2474;
  int64_t _M0L6_2atmpS2473;
  struct _M0TPC16string10StringView _M0L6_2atmpS2471;
  moonbit_string_t _M0L14num__cols__strS1048;
  int32_t _M0L1nS1051;
  int64_t _M0L7_2abindS1052;
  int32_t _M0L3valS2416;
  int32_t _M0L6_2atmpS2415;
  struct _M0TUsRPB13StringBuilderE** _M0L7_2abindS1056;
  struct _M0TUsRPB13StringBuilderE** _M0L6_2atmpS2469;
  struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE _M0L6_2atmpS2467;
  int64_t _M0L6_2atmpS2468;
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L7recordsS1055;
  struct _M0TPB5ArrayGsE* _M0L14current__namesS1057;
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L12seq__recordsS1077;
  int32_t _M0L7_2abindS1078;
  int32_t _M0L2__S1079;
  struct _M0TPB3MapGssE* _M0L6_2atmpS2465;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L6_2atmpS2466;
  struct moonbit_result_1 _result_3392;
  #line 57 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7_2abindS1025 = (moonbit_string_t)moonbit_string_literal_42.data;
  _M0L6_2atmpS2485 = Moonbit_array_length(_M0L7_2abindS1025);
  _M0L6_2atmpS2484
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS1025, .$1 = 0, .$2 = _M0L6_2atmpS2485
  };
  #line 61 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2483
  = _M0MPC16string6String5split(_M0L7contentS1024, _M0L6_2atmpS2484);
  moonbit_decref(_M0L6_2atmpS2484.$0);
  #line 61 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L5linesS1023
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2483);
  moonbit_decref(_M0L6_2atmpS2483);
  #line 62 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L1nS1026
  = _M0MPC15array5Array6lengthGRPC16string10StringViewE(_M0L5linesS1023);
  _M0L1iS1027
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS1027)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS1027->$0 = 0;
  _M0L6_2atmpS2482 = 13;
  _M0L2crS1028 = (uint16_t)_M0L6_2atmpS2482;
  _M0L6_2atmpS2481 = 32;
  _M0L2spS1029 = (uint16_t)_M0L6_2atmpS2481;
  while (1) {
    int32_t _M0L3valS2390 = _M0L1iS1027->$0;
    if (_M0L3valS2390 < _M0L1nS1026) {
      int32_t _M0L3valS2397 = _M0L1iS1027->$0;
      struct _M0TPC16string10StringView _M0L4lineS1030;
      int32_t _M0L3lenS1031;
      struct _M0TPB8MutLocalGiE* _M0L7trimmedS1032;
      int32_t _if__result_3371;
      int32_t _M0L3valS2394;
      int32_t _M0L3valS2396;
      int32_t _M0L6_2atmpS2395;
      #line 68 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L4lineS1030
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1023, _M0L3valS2397);
      #line 69 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L3lenS1031 = _M0MPC16string10StringView6length(_M0L4lineS1030);
      _M0L7trimmedS1032
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L7trimmedS1032)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L7trimmedS1032->$0 = _M0L3lenS1031;
      if (_M0L3lenS1031 > 0) {
        int32_t _M0L6_2atmpS2392 = _M0L3lenS1031 - 1;
        int32_t _M0L6_2atmpS2391;
        #line 71 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2391
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1030, _M0L6_2atmpS2392);
        moonbit_decref(_M0L4lineS1030.$0);
        _if__result_3371 = _M0L6_2atmpS2391 == _M0L2crS1028;
      } else {
        moonbit_decref(_M0L4lineS1030.$0);
        _if__result_3371 = 0;
      }
      if (_if__result_3371) {
        int32_t _M0L6_2atmpS2393 = _M0L3lenS1031 - 1;
        _M0L7trimmedS1032->$0 = _M0L6_2atmpS2393;
      }
      _M0L3valS2394 = _M0L7trimmedS1032->$0;
      moonbit_decref(_M0L7trimmedS1032);
      if (_M0L3valS2394 > 0) {
        break;
      }
      _M0L3valS2396 = _M0L1iS1027->$0;
      _M0L6_2atmpS2395 = _M0L3valS2396 + 1;
      _M0L1iS1027->$0 = _M0L6_2atmpS2395;
      continue;
    }
    break;
  }
  _M0L3valS2398 = _M0L1iS1027->$0;
  if (_M0L3valS2398 >= _M0L1nS1026) {
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2399;
    struct moonbit_result_1 _result_3372;
    moonbit_decref(_M0L1iS1027);
    moonbit_decref(_M0L5linesS1023);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2399
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2399)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2399)->$0
    = (moonbit_string_t)moonbit_string_literal_66.data;
    _result_3372.tag = 0;
    _result_3372.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2399;
    return _result_3372;
  }
  _M0L3valS2480 = _M0L1iS1027->$0;
  #line 84 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6headerS1034
  = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1023, _M0L3valS2480);
  #line 85 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L4hlenS1035 = _M0MPC16string10StringView6length(_M0L6headerS1034);
  _M0L8htrimmedS1036
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L8htrimmedS1036)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L8htrimmedS1036->$0 = _M0L4hlenS1035;
  if (_M0L4hlenS1035 > 0) {
    int32_t _M0L6_2atmpS2401 = _M0L4hlenS1035 - 1;
    int32_t _M0L6_2atmpS2400;
    #line 87 "/home/developer/Documents/1/phylip_io.mbt"
    _M0L6_2atmpS2400
    = _M0MPC16string10StringView11unsafe__get(_M0L6headerS1034, _M0L6_2atmpS2401);
    _if__result_3373 = _M0L6_2atmpS2400 == _M0L2crS1028;
  } else {
    _if__result_3373 = 0;
  }
  if (_if__result_3373) {
    int32_t _M0L6_2atmpS2402 = _M0L4hlenS1035 - 1;
    _M0L8htrimmedS1036->$0 = _M0L6_2atmpS2402;
  }
  _M0L10space__posS1037
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L10space__posS1037)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L10space__posS1037->$0 = 0;
  while (1) {
    int32_t _M0L3valS2405 = _M0L10space__posS1037->$0;
    int32_t _M0L3valS2406 = _M0L8htrimmedS1036->$0;
    int32_t _if__result_3375;
    if (_M0L3valS2405 < _M0L3valS2406) {
      int32_t _M0L3valS2404 = _M0L10space__posS1037->$0;
      int32_t _M0L6_2atmpS2403;
      #line 92 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2403
      = _M0MPC16string10StringView11unsafe__get(_M0L6headerS1034, _M0L3valS2404);
      _if__result_3375 = _M0L6_2atmpS2403 != _M0L2spS1029;
    } else {
      _if__result_3375 = 0;
    }
    if (_if__result_3375) {
      int32_t _M0L3valS2408 = _M0L10space__posS1037->$0;
      int32_t _M0L6_2atmpS2407 = _M0L3valS2408 + 1;
      _M0L10space__posS1037->$0 = _M0L6_2atmpS2407;
      continue;
    }
    break;
  }
  _M0L3valS2479 = _M0L10space__posS1037->$0;
  _M0L6_2atmpS2478 = (int64_t)_M0L3valS2479;
  #line 95 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2477
  = _M0MPC16string10StringView11sub_2einner(_M0L6headerS1034, 0, _M0L6_2atmpS2478);
  #line 95 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L14num__seqs__strS1039
  = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2477);
  moonbit_decref(_M0L6_2atmpS2477.$0);
  #line 96 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7_2abindS1043
  = _M0FP26biolab8bio__seq10parse__int(_M0L14num__seqs__strS1039);
  moonbit_decref(_M0L14num__seqs__strS1039);
  if (_M0L7_2abindS1043 == 4294967296ll) {
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2476;
    struct moonbit_result_1 _result_3377;
    moonbit_decref(_M0L10space__posS1037);
    moonbit_decref(_M0L8htrimmedS1036);
    moonbit_decref(_M0L6headerS1034.$0);
    moonbit_decref(_M0L1iS1027);
    moonbit_decref(_M0L5linesS1023);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2476
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2476)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2476)->$0
    = (moonbit_string_t)moonbit_string_literal_67.data;
    _result_3377.tag = 0;
    _result_3377.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2476;
    return _result_3377;
  } else {
    int64_t _M0L7_2aSomeS1044 = _M0L7_2abindS1043;
    int32_t _M0L4_2anS1045 = (int32_t)_M0L7_2aSomeS1044;
    _M0L1nS1042 = _M0L4_2anS1045;
    goto join_1041;
  }
  goto joinlet_3376;
  join_1041:;
  _M0L9num__seqsS1040 = _M0L1nS1042;
  joinlet_3376:;
  _M0L3valS2475 = _M0L10space__posS1037->$0;
  moonbit_decref(_M0L10space__posS1037);
  _M0L8col__posS1046
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L8col__posS1046)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L8col__posS1046->$0 = _M0L3valS2475;
  while (1) {
    int32_t _M0L3valS2411 = _M0L8col__posS1046->$0;
    int32_t _M0L3valS2412 = _M0L8htrimmedS1036->$0;
    int32_t _if__result_3379;
    if (_M0L3valS2411 < _M0L3valS2412) {
      int32_t _M0L3valS2410 = _M0L8col__posS1046->$0;
      int32_t _M0L6_2atmpS2409;
      #line 102 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2409
      = _M0MPC16string10StringView11unsafe__get(_M0L6headerS1034, _M0L3valS2410);
      _if__result_3379 = _M0L6_2atmpS2409 == _M0L2spS1029;
    } else {
      _if__result_3379 = 0;
    }
    if (_if__result_3379) {
      int32_t _M0L3valS2414 = _M0L8col__posS1046->$0;
      int32_t _M0L6_2atmpS2413 = _M0L3valS2414 + 1;
      _M0L8col__posS1046->$0 = _M0L6_2atmpS2413;
      continue;
    }
    break;
  }
  _M0L3valS2472 = _M0L8col__posS1046->$0;
  moonbit_decref(_M0L8col__posS1046);
  _M0L3valS2474 = _M0L8htrimmedS1036->$0;
  moonbit_decref(_M0L8htrimmedS1036);
  _M0L6_2atmpS2473 = (int64_t)_M0L3valS2474;
  #line 105 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L6_2atmpS2471
  = _M0MPC16string10StringView11sub_2einner(_M0L6headerS1034, _M0L3valS2472, _M0L6_2atmpS2473);
  moonbit_decref(_M0L6headerS1034.$0);
  #line 105 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L14num__cols__strS1048
  = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2471);
  moonbit_decref(_M0L6_2atmpS2471.$0);
  #line 106 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7_2abindS1052
  = _M0FP26biolab8bio__seq10parse__int(_M0L14num__cols__strS1048);
  moonbit_decref(_M0L14num__cols__strS1048);
  if (_M0L7_2abindS1052 == 4294967296ll) {
    void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2470;
    struct moonbit_result_1 _result_3381;
    moonbit_decref(_M0L1iS1027);
    moonbit_decref(_M0L5linesS1023);
    _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2470
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
    Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2470)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
    ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2470)->$0
    = (moonbit_string_t)moonbit_string_literal_67.data;
    _result_3381.tag = 0;
    _result_3381.data.err
    = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2470;
    return _result_3381;
  } else {
    int64_t _M0L7_2aSomeS1053 = _M0L7_2abindS1052;
    int32_t _M0L4_2anS1054 = (int32_t)_M0L7_2aSomeS1053;
    _M0L1nS1051 = _M0L4_2anS1054;
    goto join_1050;
  }
  goto joinlet_3380;
  join_1050:;
  joinlet_3380:;
  _M0L3valS2416 = _M0L1iS1027->$0;
  _M0L6_2atmpS2415 = _M0L3valS2416 + 1;
  _M0L1iS1027->$0 = _M0L6_2atmpS2415;
  _M0L7_2abindS1056
  = (struct _M0TUsRPB13StringBuilderE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2469 = _M0L7_2abindS1056;
  _M0L6_2atmpS2467
  = (struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE){
    .$0 = _M0L6_2atmpS2469, .$1 = 0, .$2 = 0
  };
  _M0L6_2atmpS2468 = (int64_t)_M0L9num__seqsS1040;
  #line 112 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L7recordsS1055
  = _M0MPB3Map3MapGsRPB13StringBuilderE(_M0L6_2atmpS2467, _M0L6_2atmpS2468);
  moonbit_decref(_M0L6_2atmpS2467.$0);
  #line 113 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L14current__namesS1057 = _M0MPC15array5Array11new_2einnerGsE(0);
  while (1) {
    int32_t _M0L3valS2417 = _M0L1iS1027->$0;
    if (_M0L3valS2417 < _M0L1nS1026) {
      int32_t _M0L3valS2457 = _M0L1iS1027->$0;
      struct _M0TPC16string10StringView _M0L4lineS1058;
      int32_t _M0L3lenS1059;
      struct _M0TPB8MutLocalGiE* _M0L7trimmedS1060;
      int32_t _if__result_3383;
      int32_t _M0L3valS2421;
      struct _M0TPB8MutLocalGiE* _M0L1jS1062;
      int32_t _M0L3valS2430;
      int32_t _M0L3valS2452;
      moonbit_string_t _M0L4nameS1064;
      int32_t _M0L3valS2451;
      struct _M0TPB8MutLocalGiE* _M0L1kS1065;
      int32_t _M0L3valS2448;
      int32_t _M0L3valS2450;
      int64_t _M0L6_2atmpS2449;
      struct _M0TPC16string10StringView _M0L6_2atmpS2443;
      moonbit_string_t _M0L7_2abindS1068;
      int32_t _M0L6_2atmpS2447;
      struct _M0TPC16string10StringView _M0L6_2atmpS2444;
      moonbit_string_t _M0L7_2abindS1069;
      int32_t _M0L6_2atmpS2446;
      struct _M0TPC16string10StringView _M0L6_2atmpS2445;
      struct _M0TPC16string10StringView _M0L6_2atmpS2442;
      moonbit_string_t _M0L9seq__partS1067;
      struct _M0TPB13StringBuilder* _M0L2sbS1071;
      struct _M0TPB13StringBuilder* _M0L7_2abindS1072;
      int32_t _if__result_3389;
      int32_t _M0L3valS2441;
      int32_t _M0L6_2atmpS2440;
      #line 116 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L4lineS1058
      = _M0MPC15array5Array2atGRPC16string10StringViewE(_M0L5linesS1023, _M0L3valS2457);
      #line 117 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L3lenS1059 = _M0MPC16string10StringView6length(_M0L4lineS1058);
      _M0L7trimmedS1060
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L7trimmedS1060)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L7trimmedS1060->$0 = _M0L3lenS1059;
      if (_M0L3lenS1059 > 0) {
        int32_t _M0L6_2atmpS2419 = _M0L3lenS1059 - 1;
        int32_t _M0L6_2atmpS2418;
        #line 119 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2418
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1058, _M0L6_2atmpS2419);
        _if__result_3383 = _M0L6_2atmpS2418 == _M0L2crS1028;
      } else {
        _if__result_3383 = 0;
      }
      if (_if__result_3383) {
        int32_t _M0L6_2atmpS2420 = _M0L3lenS1059 - 1;
        _M0L7trimmedS1060->$0 = _M0L6_2atmpS2420;
      }
      _M0L3valS2421 = _M0L7trimmedS1060->$0;
      if (_M0L3valS2421 == 0) {
        int32_t _M0L3valS2423;
        int32_t _M0L6_2atmpS2422;
        moonbit_decref(_M0L7trimmedS1060);
        moonbit_decref(_M0L4lineS1058.$0);
        _M0L3valS2423 = _M0L1iS1027->$0;
        _M0L6_2atmpS2422 = _M0L3valS2423 + 1;
        _M0L1iS1027->$0 = _M0L6_2atmpS2422;
        continue;
      }
      _M0L1jS1062
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1jS1062)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1jS1062->$0 = 0;
      while (1) {
        int32_t _M0L3valS2426 = _M0L1jS1062->$0;
        int32_t _M0L3valS2427 = _M0L7trimmedS1060->$0;
        int32_t _if__result_3385;
        if (_M0L3valS2426 < _M0L3valS2427) {
          int32_t _M0L3valS2425 = _M0L1jS1062->$0;
          int32_t _M0L6_2atmpS2424;
          #line 127 "/home/developer/Documents/1/phylip_io.mbt"
          _M0L6_2atmpS2424
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1058, _M0L3valS2425);
          _if__result_3385 = _M0L6_2atmpS2424 != _M0L2spS1029;
        } else {
          _if__result_3385 = 0;
        }
        if (_if__result_3385) {
          int32_t _M0L3valS2429 = _M0L1jS1062->$0;
          int32_t _M0L6_2atmpS2428 = _M0L3valS2429 + 1;
          _M0L1jS1062->$0 = _M0L6_2atmpS2428;
          continue;
        }
        break;
      }
      _M0L3valS2430 = _M0L1jS1062->$0;
      if (_M0L3valS2430 == 0) {
        int32_t _M0L3valS2432;
        int32_t _M0L6_2atmpS2431;
        moonbit_decref(_M0L1jS1062);
        moonbit_decref(_M0L7trimmedS1060);
        moonbit_decref(_M0L4lineS1058.$0);
        _M0L3valS2432 = _M0L1iS1027->$0;
        _M0L6_2atmpS2431 = _M0L3valS2432 + 1;
        _M0L1iS1027->$0 = _M0L6_2atmpS2431;
        continue;
      }
      _M0L3valS2452 = _M0L1jS1062->$0;
      if (_M0L3valS2452 > 10) {
        struct _M0TPC16string10StringView _M0L6_2atmpS2453;
        #line 134 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2453
        = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1058, 0, 10ll);
        #line 134 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L4nameS1064
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2453);
        moonbit_decref(_M0L6_2atmpS2453.$0);
      } else {
        int32_t _M0L3valS2456 = _M0L1jS1062->$0;
        int64_t _M0L6_2atmpS2455 = (int64_t)_M0L3valS2456;
        struct _M0TPC16string10StringView _M0L6_2atmpS2454;
        #line 134 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2454
        = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1058, 0, _M0L6_2atmpS2455);
        #line 134 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L4nameS1064
        = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2454);
        moonbit_decref(_M0L6_2atmpS2454.$0);
      }
      _M0L3valS2451 = _M0L1jS1062->$0;
      moonbit_decref(_M0L1jS1062);
      _M0L1kS1065
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1kS1065)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1kS1065->$0 = _M0L3valS2451;
      while (1) {
        int32_t _M0L3valS2435 = _M0L1kS1065->$0;
        int32_t _M0L3valS2436 = _M0L7trimmedS1060->$0;
        int32_t _if__result_3387;
        if (_M0L3valS2435 < _M0L3valS2436) {
          int32_t _M0L3valS2434 = _M0L1kS1065->$0;
          int32_t _M0L6_2atmpS2433;
          #line 136 "/home/developer/Documents/1/phylip_io.mbt"
          _M0L6_2atmpS2433
          = _M0MPC16string10StringView11unsafe__get(_M0L4lineS1058, _M0L3valS2434);
          _if__result_3387 = _M0L6_2atmpS2433 == _M0L2spS1029;
        } else {
          _if__result_3387 = 0;
        }
        if (_if__result_3387) {
          int32_t _M0L3valS2438 = _M0L1kS1065->$0;
          int32_t _M0L6_2atmpS2437 = _M0L3valS2438 + 1;
          _M0L1kS1065->$0 = _M0L6_2atmpS2437;
          continue;
        }
        break;
      }
      _M0L3valS2448 = _M0L1kS1065->$0;
      moonbit_decref(_M0L1kS1065);
      _M0L3valS2450 = _M0L7trimmedS1060->$0;
      moonbit_decref(_M0L7trimmedS1060);
      _M0L6_2atmpS2449 = (int64_t)_M0L3valS2450;
      #line 139 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2443
      = _M0MPC16string10StringView11sub_2einner(_M0L4lineS1058, _M0L3valS2448, _M0L6_2atmpS2449);
      moonbit_decref(_M0L4lineS1058.$0);
      _M0L7_2abindS1068 = (moonbit_string_t)moonbit_string_literal_65.data;
      _M0L6_2atmpS2447 = Moonbit_array_length(_M0L7_2abindS1068);
      _M0L6_2atmpS2444
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L7_2abindS1068, .$1 = 0, .$2 = _M0L6_2atmpS2447
      };
      _M0L7_2abindS1069 = (moonbit_string_t)moonbit_string_literal_0.data;
      _M0L6_2atmpS2446 = Moonbit_array_length(_M0L7_2abindS1069);
      _M0L6_2atmpS2445
      = (struct _M0TPC16string10StringView){
        .$0 = _M0L7_2abindS1069, .$1 = 0, .$2 = _M0L6_2atmpS2446
      };
      #line 139 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2442
      = _M0MPC16string10StringView7replace(_M0L6_2atmpS2443, _M0L6_2atmpS2444, _M0L6_2atmpS2445);
      moonbit_decref(_M0L6_2atmpS2443.$0);
      moonbit_decref(_M0L6_2atmpS2444.$0);
      moonbit_decref(_M0L6_2atmpS2445.$0);
      #line 139 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L9seq__partS1067
      = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2442);
      moonbit_decref(_M0L6_2atmpS2442.$0);
      #line 141 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L7_2abindS1072
      = _M0MPB3Map3getGsRPB13StringBuilderE(_M0L7recordsS1055, _M0L4nameS1064);
      if (_M0L7_2abindS1072 == 0) {
        struct _M0TPB13StringBuilder* _M0L2sbS1075;
        if (_M0L7_2abindS1072) {
          moonbit_decref(_M0L7_2abindS1072);
        }
        #line 144 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L2sbS1075 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
        #line 145 "/home/developer/Documents/1/phylip_io.mbt"
        _M0IPB13StringBuilderPB6Logger13write__string(_M0L2sbS1075, _M0L9seq__partS1067);
        moonbit_decref(_M0L9seq__partS1067);
        #line 146 "/home/developer/Documents/1/phylip_io.mbt"
        _M0MPB3Map3setGsRPB13StringBuilderE(_M0L7recordsS1055, _M0L4nameS1064, _M0L2sbS1075);
        moonbit_decref(_M0L2sbS1075);
        #line 147 "/home/developer/Documents/1/phylip_io.mbt"
        _M0MPC15array5Array4pushGsE(_M0L14current__namesS1057, _M0L4nameS1064);
        moonbit_decref(_M0L4nameS1064);
      } else {
        struct _M0TPB13StringBuilder* _M0L7_2aSomeS1073;
        struct _M0TPB13StringBuilder* _M0L5_2asbS1074;
        moonbit_decref(_M0L4nameS1064);
        _M0L7_2aSomeS1073 = _M0L7_2abindS1072;
        _M0L5_2asbS1074 = _M0L7_2aSomeS1073;
        _M0L2sbS1071 = _M0L5_2asbS1074;
        goto join_1070;
      }
      goto joinlet_3388;
      join_1070:;
      #line 142 "/home/developer/Documents/1/phylip_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L2sbS1071, _M0L9seq__partS1067);
      moonbit_decref(_M0L2sbS1071);
      moonbit_decref(_M0L9seq__partS1067);
      joinlet_3388:;
      if (_M0L10sequentialS1076) {
        int32_t _M0L6_2atmpS2439;
        #line 150 "/home/developer/Documents/1/phylip_io.mbt"
        _M0L6_2atmpS2439
        = _M0MPC15array5Array6lengthGsE(_M0L14current__namesS1057);
        _if__result_3389 = _M0L6_2atmpS2439 >= _M0L9num__seqsS1040;
      } else {
        _if__result_3389 = 0;
      }
      if (_if__result_3389) {
        moonbit_decref(_M0L1iS1027);
        moonbit_decref(_M0L5linesS1023);
        break;
      }
      _M0L3valS2441 = _M0L1iS1027->$0;
      _M0L6_2atmpS2440 = _M0L3valS2441 + 1;
      _M0L1iS1027->$0 = _M0L6_2atmpS2440;
      continue;
    } else {
      moonbit_decref(_M0L1iS1027);
      moonbit_decref(_M0L5linesS1023);
    }
    break;
  }
  #line 156 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L12seq__recordsS1077
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS1078 = _M0L14current__namesS1057->$1;
  _M0L2__S1079 = 0;
  while (1) {
    if (_M0L2__S1079 < _M0L7_2abindS1078) {
      moonbit_string_t* _M0L3bufS2464 = _M0L14current__namesS1057->$0;
      moonbit_string_t _M0L4nameS1080 =
        (moonbit_string_t)_M0L3bufS2464[_M0L2__S1079];
      struct _M0TPB13StringBuilder* _M0L2sbS1082;
      struct _M0TPB13StringBuilder* _M0L7_2abindS1085;
      moonbit_string_t _M0L6_2atmpS2462;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS1083;
      struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS1084;
      struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2461;
      struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2460;
      struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2459;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2458;
      int32_t _M0L6_2atmpS2463;
      moonbit_incref(_M0L4nameS1080);
      #line 158 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L7_2abindS1085
      = _M0MPB3Map3getGsRPB13StringBuilderE(_M0L7recordsS1055, _M0L4nameS1080);
      if (_M0L7_2abindS1085 == 0) {
        if (_M0L7_2abindS1085) {
          moonbit_decref(_M0L7_2abindS1085);
        }
        moonbit_decref(_M0L4nameS1080);
      } else {
        struct _M0TPB13StringBuilder* _M0L7_2aSomeS1086 = _M0L7_2abindS1085;
        struct _M0TPB13StringBuilder* _M0L5_2asbS1087 = _M0L7_2aSomeS1086;
        _M0L2sbS1082 = _M0L5_2asbS1087;
        goto join_1081;
      }
      goto joinlet_3391;
      join_1081:;
      #line 160 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2462 = _M0MPB13StringBuilder10to__string(_M0L2sbS1082);
      moonbit_decref(_M0L2sbS1082);
      #line 160 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L3seqS1083 = _M0MP26biolab8bio__seq3Seq3new(_M0L6_2atmpS2462);
      moonbit_decref(_M0L6_2atmpS2462);
      _M0L7_2abindS1084
      = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
      _M0L6_2atmpS2461 = _M0L7_2abindS1084;
      _M0L6_2atmpS2460
      = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
        .$0 = _M0L6_2atmpS2461, .$1 = 0, .$2 = 0
      };
      #line 166 "/home/developer/Documents/1/phylip_io.mbt"
      _M0L6_2atmpS2459 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2460, 0ll);
      moonbit_decref(_M0L6_2atmpS2460.$0);
      moonbit_incref(_M0L4nameS1080);
      moonbit_incref(_M0L4nameS1080);
      _M0L6_2atmpS2458
      = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
      Moonbit_object_header(_M0L6_2atmpS2458)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
      _M0L6_2atmpS2458->$0 = _M0L3seqS1083;
      _M0L6_2atmpS2458->$1 = _M0L4nameS1080;
      _M0L6_2atmpS2458->$2 = _M0L4nameS1080;
      _M0L6_2atmpS2458->$3 = _M0L4nameS1080;
      _M0L6_2atmpS2458->$4 = _M0L6_2atmpS2459;
      #line 161 "/home/developer/Documents/1/phylip_io.mbt"
      _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L12seq__recordsS1077, _M0L6_2atmpS2458);
      moonbit_decref(_M0L6_2atmpS2458);
      joinlet_3391:;
      _M0L6_2atmpS2463 = _M0L2__S1079 + 1;
      _M0L2__S1079 = _M0L6_2atmpS2463;
      continue;
    } else {
      moonbit_decref(_M0L14current__namesS1057);
      moonbit_decref(_M0L7recordsS1055);
    }
    break;
  }
  _M0L6_2atmpS2465 = 0;
  _M0L6_2atmpS2466 = 0;
  #line 172 "/home/developer/Documents/1/phylip_io.mbt"
  _result_3392
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(_M0L12seq__recordsS1077, _M0L6_2atmpS2465, _M0L6_2atmpS2466);
  moonbit_decref(_M0L12seq__recordsS1077);
  if (_M0L6_2atmpS2465) {
    moonbit_decref(_M0L6_2atmpS2465);
  }
  if (_M0L6_2atmpS2466) {
    moonbit_decref(_M0L6_2atmpS2466);
  }
  return _result_3392;
}

struct moonbit_result_1 _M0MP26biolab8bio__seq20MultipleSeqAlignment3new(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1022,
  struct _M0TPB3MapGssE* _M0L17annotations_2eoptS1015,
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L25column__annotations_2eoptS1019
) {
  struct _M0TPB3MapGssE* _M0L11annotationsS1014;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L19column__annotationsS1018;
  struct moonbit_result_1 _result_3393;
  if (_M0L17annotations_2eoptS1015 == 0) {
    struct _M0TUssE** _M0L7_2abindS1017 =
      (struct _M0TUssE**)moonbit_empty_ref_array;
    struct _M0TUssE** _M0L6_2atmpS2389 = _M0L7_2abindS1017;
    struct _M0TPB9ArrayViewGUssEE _M0L6_2atmpS2388 =
      (struct _M0TPB9ArrayViewGUssEE){.$0 = _M0L6_2atmpS2389,
                                        .$1 = 0,
                                        .$2 = 0};
    #line 32 "/home/developer/Documents/1/align.mbt"
    _M0L11annotationsS1014 = _M0MPB3Map3MapGssE(_M0L6_2atmpS2388, 0ll);
    moonbit_decref(_M0L6_2atmpS2388.$0);
  } else {
    struct _M0TPB3MapGssE* _M0L7_2aSomeS1016 = _M0L17annotations_2eoptS1015;
    if (_M0L7_2aSomeS1016) {
      moonbit_incref(_M0L7_2aSomeS1016);
    }
    _M0L11annotationsS1014 = _M0L7_2aSomeS1016;
  }
  if (_M0L25column__annotations_2eoptS1019 == 0) {
    struct _M0TUsRPB5ArrayGsEE** _M0L7_2abindS1021 =
      (struct _M0TUsRPB5ArrayGsEE**)moonbit_empty_ref_array;
    struct _M0TUsRPB5ArrayGsEE** _M0L6_2atmpS2387 = _M0L7_2abindS1021;
    struct _M0TPB9ArrayViewGUsRPB5ArrayGsEEE _M0L6_2atmpS2386 =
      (struct _M0TPB9ArrayViewGUsRPB5ArrayGsEEE){.$0 = _M0L6_2atmpS2387,
                                                   .$1 = 0,
                                                   .$2 = 0};
    #line 33 "/home/developer/Documents/1/align.mbt"
    _M0L19column__annotationsS1018
    = _M0MPB3Map3MapGsRPB5ArrayGsEE(_M0L6_2atmpS2386, 0ll);
    moonbit_decref(_M0L6_2atmpS2386.$0);
  } else {
    struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L7_2aSomeS1020 =
      _M0L25column__annotations_2eoptS1019;
    if (_M0L7_2aSomeS1020) {
      moonbit_incref(_M0L7_2aSomeS1020);
    }
    _M0L19column__annotationsS1018 = _M0L7_2aSomeS1020;
  }
  _result_3393
  = _M0MP26biolab8bio__seq20MultipleSeqAlignment11new_2einner(_M0L7recordsS1022, _M0L11annotationsS1014, _M0L19column__annotationsS1018);
  moonbit_decref(_M0L11annotationsS1014);
  moonbit_decref(_M0L19column__annotationsS1018);
  return _result_3393;
}

struct moonbit_result_1 _M0MP26biolab8bio__seq20MultipleSeqAlignment11new_2einner(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS1006,
  struct _M0TPB3MapGssE* _M0L11annotationsS1012,
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L19column__annotationsS1013
) {
  int32_t _M0L6_2atmpS2379;
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS2385;
  struct moonbit_result_1 _result_3396;
  #line 30 "/home/developer/Documents/1/align.mbt"
  #line 35 "/home/developer/Documents/1/align.mbt"
  _M0L6_2atmpS2379
  = _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1006);
  if (_M0L6_2atmpS2379 > 0) {
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2384;
    int32_t _M0L10first__lenS1007;
    int32_t _M0L7_2abindS1008;
    int32_t _M0L2__S1009;
    #line 36 "/home/developer/Documents/1/align.mbt"
    _M0L6_2atmpS2384
    = _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS1006, 0);
    #line 36 "/home/developer/Documents/1/align.mbt"
    _M0L10first__lenS1007
    = _M0MP26biolab8bio__seq9SeqRecord6length(_M0L6_2atmpS2384);
    moonbit_decref(_M0L6_2atmpS2384);
    _M0L7_2abindS1008 = _M0L7recordsS1006->$1;
    _M0L2__S1009 = 0;
    while (1) {
      if (_M0L2__S1009 < _M0L7_2abindS1008) {
        struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2383 =
          _M0L7recordsS1006->$0;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS1010 =
          (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2383[
            _M0L2__S1009
          ];
        int32_t _M0L6_2atmpS2380;
        int32_t _M0L6_2atmpS2382;
        moonbit_incref(_M0L3recS1010);
        #line 38 "/home/developer/Documents/1/align.mbt"
        _M0L6_2atmpS2380
        = _M0MP26biolab8bio__seq9SeqRecord6length(_M0L3recS1010);
        moonbit_decref(_M0L3recS1010);
        if (_M0L6_2atmpS2380 != _M0L10first__lenS1007) {
          void* _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2381 =
            (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError));
          struct moonbit_result_1 _result_3395;
          Moonbit_object_header(_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2381)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 24, 1);
          ((struct _M0DTPC15error5Error47biolab_2fbio__seq_2eAlignIOError_2eAlignIOError*)_M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2381)->$0
          = (moonbit_string_t)moonbit_string_literal_68.data;
          _result_3395.tag = 0;
          _result_3395.data.err
          = _M0L47biolab_2fbio__seq_2eAlignIOError_2eAlignIOErrorS2381;
          return _result_3395;
        }
        _M0L6_2atmpS2382 = _M0L2__S1009 + 1;
        _M0L2__S1009 = _M0L6_2atmpS2382;
        continue;
      }
      break;
    }
  }
  moonbit_incref(_M0L7recordsS1006);
  moonbit_incref(_M0L11annotationsS1012);
  moonbit_incref(_M0L19column__annotationsS1013);
  _M0L6_2atmpS2385
  = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq20MultipleSeqAlignment));
  Moonbit_object_header(_M0L6_2atmpS2385)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 34, 0);
  _M0L6_2atmpS2385->$0 = _M0L7recordsS1006;
  _M0L6_2atmpS2385->$1 = _M0L11annotationsS1012;
  _M0L6_2atmpS2385->$2 = _M0L19column__annotationsS1013;
  _result_3396.tag = 1;
  _result_3396.data.ok = _M0L6_2atmpS2385;
  return _result_3396;
}

int64_t _M0FP26biolab8bio__seq10parse__int(moonbit_string_t _M0L1sS1000) {
  int32_t _M0L1nS999;
  #line 21 "/home/developer/Documents/1/phylip_io.mbt"
  _M0L1nS999 = Moonbit_array_length(_M0L1sS1000);
  if (_M0L1nS999 == 0) {
    return 4294967296ll;
  } else {
    struct _M0TPB8MutLocalGiE* _M0L6resultS1001 =
      (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
    int32_t _M0L4zeroS1002;
    int32_t _M0L1iS1003;
    int32_t _M0L3valS2378;
    Moonbit_object_header(_M0L6resultS1001)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
    _M0L6resultS1001->$0 = 0;
    _M0L4zeroS1002 = 48;
    _M0L1iS1003 = 0;
    while (1) {
      if (_M0L1iS1003 < _M0L1nS999) {
        int32_t _M0L6_2atmpS2376 = _M0L1sS1000[_M0L1iS1003];
        int32_t _M0L1cS1004 = (int32_t)_M0L6_2atmpS2376;
        int32_t _if__result_3398;
        int32_t _M0L3valS2375;
        int32_t _M0L6_2atmpS2373;
        int32_t _M0L6_2atmpS2374;
        int32_t _M0L6_2atmpS2372;
        int32_t _M0L6_2atmpS2377;
        if (_M0L1cS1004 < _M0L4zeroS1002) {
          _if__result_3398 = 1;
        } else {
          int32_t _M0L6_2atmpS2371 = _M0L4zeroS1002 + 9;
          _if__result_3398 = _M0L1cS1004 > _M0L6_2atmpS2371;
        }
        if (_if__result_3398) {
          moonbit_decref(_M0L6resultS1001);
          return 4294967296ll;
        }
        _M0L3valS2375 = _M0L6resultS1001->$0;
        _M0L6_2atmpS2373 = _M0L3valS2375 * 10;
        _M0L6_2atmpS2374 = _M0L1cS1004 - _M0L4zeroS1002;
        _M0L6_2atmpS2372 = _M0L6_2atmpS2373 + _M0L6_2atmpS2374;
        _M0L6resultS1001->$0 = _M0L6_2atmpS2372;
        _M0L6_2atmpS2377 = _M0L1iS1003 + 1;
        _M0L1iS1003 = _M0L6_2atmpS2377;
        continue;
      }
      break;
    }
    _M0L3valS2378 = _M0L6resultS1001->$0;
    moonbit_decref(_M0L6resultS1001);
    return (int64_t)_M0L3valS2378;
  }
}

moonbit_string_t _M0FP26biolab8bio__seq12write__fasta(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS989
) {
  struct _M0TPB13StringBuilder* _M0L3bufS987;
  int32_t _M0L7_2abindS988;
  int32_t _M0L2__S990;
  moonbit_string_t _result_3401;
  #line 99 "/home/developer/Documents/1/fasta_io.mbt"
  #line 100 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L3bufS987 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L7_2abindS988 = _M0L7recordsS989->$1;
  _M0L2__S990 = 0;
  while (1) {
    if (_M0L2__S990 < _M0L7_2abindS988) {
      struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS2370 =
        _M0L7recordsS989->$0;
      struct _M0TP26biolab8bio__seq9SeqRecord* _M0L3recS991 =
        (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS2370[_M0L2__S990];
      moonbit_string_t _M0L2idS2367 = _M0L3recS991->$1;
      moonbit_string_t _M0L11descriptionS2368 = _M0L3recS991->$3;
      moonbit_string_t _M0L5titleS992;
      struct _M0TP26biolab8bio__seq3Seq* _M0L8_2afieldS2998;
      int32_t _M0L6_2acntS3234;
      struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2366;
      moonbit_string_t _M0L3seqS993;
      int32_t _M0L1nS994;
      struct _M0TPB8MutLocalGiE* _M0L1iS995;
      int32_t _M0L6_2atmpS2369;
      moonbit_incref(_M0L11descriptionS2368);
      moonbit_incref(_M0L2idS2367);
      moonbit_incref(_M0L3recS991);
      #line 107 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L5titleS992
      = _M0FP26biolab8bio__seq19build__fasta__title(_M0L2idS2367, _M0L11descriptionS2368);
      moonbit_decref(_M0L2idS2367);
      moonbit_decref(_M0L11descriptionS2368);
      #line 108 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS987, 62);
      #line 109 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS987, _M0L5titleS992);
      moonbit_decref(_M0L5titleS992);
      #line 110 "/home/developer/Documents/1/fasta_io.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS987, 10);
      _M0L8_2afieldS2998 = _M0L3recS991->$0;
      _M0L6_2acntS3234
      = Moonbit_rc_count(Moonbit_object_header(_M0L3recS991));
      if (_M0L6_2acntS3234 > 1) {
        int32_t _M0L11_2anew__cntS3239 = _M0L6_2acntS3234 - 1;
        Moonbit_set_rc_count(Moonbit_object_header(_M0L3recS991), _M0L11_2anew__cntS3239);
        moonbit_incref(_M0L8_2afieldS2998);
      } else if (_M0L6_2acntS3234 == 1) {
        struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L8_2afieldS3238 =
          _M0L3recS991->$4;
        moonbit_string_t _M0L8_2afieldS3237;
        moonbit_string_t _M0L8_2afieldS3236;
        moonbit_string_t _M0L8_2afieldS3235;
        moonbit_decref(_M0L8_2afieldS3238);
        _M0L8_2afieldS3237 = _M0L3recS991->$3;
        moonbit_decref(_M0L8_2afieldS3237);
        _M0L8_2afieldS3236 = _M0L3recS991->$2;
        moonbit_decref(_M0L8_2afieldS3236);
        _M0L8_2afieldS3235 = _M0L3recS991->$1;
        moonbit_decref(_M0L8_2afieldS3235);
        #line 112 "/home/developer/Documents/1/fasta_io.mbt"
        moonbit_free(_M0L3recS991);
      }
      _M0L3seqS2366 = _M0L8_2afieldS2998;
      #line 112 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L3seqS993 = _M0MP26biolab8bio__seq3Seq10to__string(_M0L3seqS2366);
      moonbit_decref(_M0L3seqS2366);
      _M0L1nS994 = Moonbit_array_length(_M0L3seqS993);
      _M0L1iS995
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1iS995)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1iS995->$0 = 0;
      while (1) {
        int32_t _M0L3valS2358 = _M0L1iS995->$0;
        if (_M0L3valS2358 < _M0L1nS994) {
          int32_t _M0L3valS2364 = _M0L1iS995->$0;
          int32_t _M0L6_2atmpS2363 =
            _M0L3valS2364 + _M0FP26biolab8bio__seq18fasta__line__width;
          int32_t _M0L3endS996;
          int32_t _M0L3valS2361;
          int64_t _M0L6_2atmpS2362;
          struct _M0TPC16string10StringView _M0L6_2atmpS2360;
          moonbit_string_t _M0L6_2atmpS2359;
          if (_M0L6_2atmpS2363 < _M0L1nS994) {
            int32_t _M0L3valS2365 = _M0L1iS995->$0;
            _M0L3endS996
            = _M0L3valS2365 + _M0FP26biolab8bio__seq18fasta__line__width;
          } else {
            _M0L3endS996 = _M0L1nS994;
          }
          _M0L3valS2361 = _M0L1iS995->$0;
          _M0L6_2atmpS2362 = (int64_t)_M0L3endS996;
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2360
          = _M0MPC16string6String11sub_2einner(_M0L3seqS993, _M0L3valS2361, _M0L6_2atmpS2362);
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2359
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2360);
          moonbit_decref(_M0L6_2atmpS2360.$0);
          #line 117 "/home/developer/Documents/1/fasta_io.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS987, _M0L6_2atmpS2359);
          moonbit_decref(_M0L6_2atmpS2359);
          #line 118 "/home/developer/Documents/1/fasta_io.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3bufS987, 10);
          _M0L1iS995->$0 = _M0L3endS996;
          continue;
        } else {
          moonbit_decref(_M0L1iS995);
          moonbit_decref(_M0L3seqS993);
        }
        break;
      }
      _M0L6_2atmpS2369 = _M0L2__S990 + 1;
      _M0L2__S990 = _M0L6_2atmpS2369;
      continue;
    }
    break;
  }
  #line 122 "/home/developer/Documents/1/fasta_io.mbt"
  _result_3401 = _M0MPB13StringBuilder10to__string(_M0L3bufS987);
  moonbit_decref(_M0L3bufS987);
  return _result_3401;
}

moonbit_string_t _M0FP26biolab8bio__seq19build__fasta__title(
  moonbit_string_t _M0L2idS986,
  moonbit_string_t _M0L11descriptionS985
) {
  int32_t _M0L6_2atmpS2354;
  #line 129 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2354 = Moonbit_array_length(_M0L11descriptionS985);
  if (_M0L6_2atmpS2354 == 0) {
    moonbit_incref(_M0L2idS986);
    return _M0L2idS986;
  } else {
    int32_t _M0L6_2atmpS2356 = Moonbit_array_length(_M0L2idS986);
    struct _M0TPC16string10StringView _M0L6_2atmpS2355;
    int32_t _result_3402;
    moonbit_incref(_M0L2idS986);
    _M0L6_2atmpS2355
    = (struct _M0TPC16string10StringView){
      .$0 = _M0L2idS986, .$1 = 0, .$2 = _M0L6_2atmpS2356
    };
    #line 132 "/home/developer/Documents/1/fasta_io.mbt"
    _result_3402
    = _M0MPC16string6String11has__prefix(_M0L11descriptionS985, _M0L6_2atmpS2355);
    moonbit_decref(_M0L6_2atmpS2355.$0);
    if (_result_3402) {
      moonbit_incref(_M0L11descriptionS985);
      return _M0L11descriptionS985;
    } else {
      moonbit_string_t _M0L6_2atmpS2357;
      moonbit_string_t _result_3403;
      #line 135 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L6_2atmpS2357
      = moonbit_add_string(_M0L2idS986, (moonbit_string_t)moonbit_string_literal_65.data);
      #line 135 "/home/developer/Documents/1/fasta_io.mbt"
      _result_3403
      = moonbit_add_string(_M0L6_2atmpS2357, _M0L11descriptionS985);
      moonbit_decref(_M0L6_2atmpS2357);
      return _result_3403;
    }
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0FP26biolab8bio__seq12parse__fasta(
  moonbit_string_t _M0L7contentS957
) {
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L7recordsS955;
  moonbit_string_t _M0L7_2abindS958;
  int32_t _M0L6_2atmpS2353;
  struct _M0TPC16string10StringView _M0L6_2atmpS2352;
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L6_2atmpS2351;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L5linesS956;
  moonbit_string_t _M0L6_2atmpS2350;
  struct _M0TPB8MutLocalGOsE* _M0L14current__titleS959;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS2349;
  struct _M0TPB8MutLocalGRPB13StringBuilderE* _M0L12current__seqS960;
  int32_t _M0L6_2atmpS2348;
  int32_t _M0L2gtS961;
  int32_t _M0L6_2atmpS2347;
  int32_t _M0L2crS962;
  int32_t _M0L6_2atmpS2346;
  int32_t _M0L2spS963;
  int32_t _M0L6_2atmpS2345;
  int32_t _M0L3tabS964;
  int32_t _M0L7_2abindS965;
  int32_t _M0L2__S966;
  moonbit_string_t _M0L5titleS981;
  moonbit_string_t _M0L8_2afieldS3004;
  int32_t _M0L6_2acntS3242;
  moonbit_string_t _M0L7_2abindS982;
  struct _M0TPB13StringBuilder* _M0L8_2afieldS3003;
  int32_t _M0L6_2acntS3240;
  struct _M0TPB13StringBuilder* _M0L3valS2344;
  moonbit_string_t _M0L6_2atmpS2343;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2342;
  #line 24 "/home/developer/Documents/1/fasta_io.mbt"
  #line 25 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7recordsS955
  = _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(0);
  _M0L7_2abindS958 = (moonbit_string_t)moonbit_string_literal_42.data;
  _M0L6_2atmpS2353 = Moonbit_array_length(_M0L7_2abindS958);
  _M0L6_2atmpS2352
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS958, .$1 = 0, .$2 = _M0L6_2atmpS2353
  };
  #line 26 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2351
  = _M0MPC16string6String5split(_M0L7contentS957, _M0L6_2atmpS2352);
  moonbit_decref(_M0L6_2atmpS2352.$0);
  #line 26 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L5linesS956
  = _M0MPB4Iter9to__arrayGRPC16string10StringViewE(_M0L6_2atmpS2351);
  moonbit_decref(_M0L6_2atmpS2351);
  _M0L6_2atmpS2350 = 0;
  _M0L14current__titleS959
  = (struct _M0TPB8MutLocalGOsE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGOsE));
  Moonbit_object_header(_M0L14current__titleS959)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 39, 0);
  _M0L14current__titleS959->$0 = _M0L6_2atmpS2350;
  #line 28 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2349 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L12current__seqS960
  = (struct _M0TPB8MutLocalGRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGRPB13StringBuilderE));
  Moonbit_object_header(_M0L12current__seqS960)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 42, 0);
  _M0L12current__seqS960->$0 = _M0L6_2atmpS2349;
  _M0L6_2atmpS2348 = 62;
  _M0L2gtS961 = (uint16_t)_M0L6_2atmpS2348;
  _M0L6_2atmpS2347 = 13;
  _M0L2crS962 = (uint16_t)_M0L6_2atmpS2347;
  _M0L6_2atmpS2346 = 32;
  _M0L2spS963 = (uint16_t)_M0L6_2atmpS2346;
  _M0L6_2atmpS2345 = 9;
  _M0L3tabS964 = (uint16_t)_M0L6_2atmpS2345;
  _M0L7_2abindS965 = _M0L5linesS956->$1;
  _M0L2__S966 = 0;
  while (1) {
    if (_M0L2__S966 < _M0L7_2abindS965) {
      struct _M0TPC16string10StringView* _M0L3bufS2341 = _M0L5linesS956->$0;
      struct _M0TPC16string10StringView _M0L4lineS967 =
        _M0L3bufS2341[_M0L2__S966];
      int32_t _M0L1nS970;
      struct _M0TPB8MutLocalGiE* _M0L12trimmed__lenS971;
      int32_t _if__result_3406;
      int32_t _M0L3valS2325;
      int32_t _M0L6_2atmpS2326;
      int32_t _M0L6_2atmpS2321;
      moonbit_incref(_M0L4lineS967.$0);
      #line 35 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L1nS970 = _M0MPC16string10StringView6length(_M0L4lineS967);
      _M0L12trimmed__lenS971
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L12trimmed__lenS971)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L12trimmed__lenS971->$0 = _M0L1nS970;
      if (_M0L1nS970 > 0) {
        int32_t _M0L6_2atmpS2323 = _M0L1nS970 - 1;
        int32_t _M0L6_2atmpS2322;
        #line 37 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2322
        = _M0MPC16string10StringView11unsafe__get(_M0L4lineS967, _M0L6_2atmpS2323);
        _if__result_3406 = _M0L6_2atmpS2322 == _M0L2crS962;
      } else {
        _if__result_3406 = 0;
      }
      if (_if__result_3406) {
        int32_t _M0L6_2atmpS2324 = _M0L1nS970 - 1;
        _M0L12trimmed__lenS971->$0 = _M0L6_2atmpS2324;
      }
      _M0L3valS2325 = _M0L12trimmed__lenS971->$0;
      if (_M0L3valS2325 == 0) {
        moonbit_decref(_M0L12trimmed__lenS971);
        moonbit_decref(_M0L4lineS967.$0);
        goto join_968;
      }
      #line 43 "/home/developer/Documents/1/fasta_io.mbt"
      _M0L6_2atmpS2326
      = _M0MPC16string10StringView11unsafe__get(_M0L4lineS967, 0);
      if (_M0L6_2atmpS2326 == _M0L2gtS961) {
        moonbit_string_t _M0L5titleS973;
        moonbit_string_t _M0L7_2abindS974 = _M0L14current__titleS959->$0;
        struct _M0TPB13StringBuilder* _M0L3valS2329;
        moonbit_string_t _M0L6_2atmpS2328;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS2327;
        int32_t _M0L3valS2332;
        moonbit_string_t _M0L6_2atmpS2331;
        moonbit_string_t _M0L6_2atmpS2330;
        moonbit_string_t _M0L6_2aoldS3006;
        struct _M0TPB13StringBuilder* _M0L6_2atmpS2336;
        struct _M0TPB13StringBuilder* _M0L6_2aoldS3005;
        if (_M0L7_2abindS974 == 0) {
          
        } else {
          moonbit_string_t _M0L7_2aSomeS975 = _M0L7_2abindS974;
          moonbit_string_t _M0L8_2atitleS976 = _M0L7_2aSomeS975;
          moonbit_incref(_M0L8_2atitleS976);
          _M0L5titleS973 = _M0L8_2atitleS976;
          goto join_972;
        }
        goto joinlet_3407;
        join_972:;
        _M0L3valS2329 = _M0L12current__seqS960->$0;
        moonbit_incref(_M0L3valS2329);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2328 = _M0MPB13StringBuilder10to__string(_M0L3valS2329);
        moonbit_decref(_M0L3valS2329);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2327
        = _M0FP26biolab8bio__seq24fasta__title__to__record(_M0L5titleS973, _M0L6_2atmpS2328);
        moonbit_decref(_M0L5titleS973);
        moonbit_decref(_M0L6_2atmpS2328);
        #line 46 "/home/developer/Documents/1/fasta_io.mbt"
        _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS955, _M0L6_2atmpS2327);
        moonbit_decref(_M0L6_2atmpS2327);
        joinlet_3407:;
        _M0L3valS2332 = _M0L12trimmed__lenS971->$0;
        if (_M0L3valS2332 > 1) {
          int32_t _M0L3valS2335 = _M0L12trimmed__lenS971->$0;
          int64_t _M0L6_2atmpS2334;
          struct _M0TPC16string10StringView _M0L6_2atmpS2333;
          moonbit_decref(_M0L12trimmed__lenS971);
          _M0L6_2atmpS2334 = (int64_t)_M0L3valS2335;
          #line 51 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2333
          = _M0MPC16string10StringView11sub_2einner(_M0L4lineS967, 1, _M0L6_2atmpS2334);
          moonbit_decref(_M0L4lineS967.$0);
          #line 51 "/home/developer/Documents/1/fasta_io.mbt"
          _M0L6_2atmpS2331
          = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2333);
          moonbit_decref(_M0L6_2atmpS2333.$0);
        } else {
          moonbit_decref(_M0L12trimmed__lenS971);
          moonbit_decref(_M0L4lineS967.$0);
          _M0L6_2atmpS2331 = (moonbit_string_t)moonbit_string_literal_0.data;
        }
        _M0L6_2atmpS2330 = _M0L6_2atmpS2331;
        _M0L6_2aoldS3006 = _M0L14current__titleS959->$0;
        if (_M0L6_2aoldS3006) {
          moonbit_decref(_M0L6_2aoldS3006);
        }
        _M0L14current__titleS959->$0 = _M0L6_2atmpS2330;
        #line 56 "/home/developer/Documents/1/fasta_io.mbt"
        _M0L6_2atmpS2336 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
        _M0L6_2aoldS3005 = _M0L12current__seqS960->$0;
        moonbit_decref(_M0L6_2aoldS3005);
        _M0L12current__seqS960->$0 = _M0L6_2atmpS2336;
      } else {
        int32_t _M0L1iS977 = 0;
        while (1) {
          int32_t _M0L3valS2337 = _M0L12trimmed__lenS971->$0;
          if (_M0L1iS977 < _M0L3valS2337) {
            int32_t _M0L1cS978;
            int32_t _if__result_3409;
            int32_t _M0L6_2atmpS2340;
            #line 59 "/home/developer/Documents/1/fasta_io.mbt"
            _M0L1cS978
            = _M0MPC16string10StringView11unsafe__get(_M0L4lineS967, _M0L1iS977);
            if (_M0L1cS978 != _M0L2spS963) {
              _if__result_3409 = _M0L1cS978 != _M0L3tabS964;
            } else {
              _if__result_3409 = 0;
            }
            if (_if__result_3409) {
              struct _M0TPB13StringBuilder* _M0L3valS2338 =
                _M0L12current__seqS960->$0;
              int32_t _M0L6_2atmpS2339;
              moonbit_incref(_M0L3valS2338);
              #line 61 "/home/developer/Documents/1/fasta_io.mbt"
              _M0L6_2atmpS2339
              = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS978);
              #line 61 "/home/developer/Documents/1/fasta_io.mbt"
              _M0IPB13StringBuilderPB6Logger11write__char(_M0L3valS2338, _M0L6_2atmpS2339);
              moonbit_decref(_M0L3valS2338);
            }
            _M0L6_2atmpS2340 = _M0L1iS977 + 1;
            _M0L1iS977 = _M0L6_2atmpS2340;
            continue;
          } else {
            moonbit_decref(_M0L12trimmed__lenS971);
            moonbit_decref(_M0L4lineS967.$0);
          }
          break;
        }
      }
      goto join_968;
      goto joinlet_3405;
      join_968:;
      _M0L6_2atmpS2321 = _M0L2__S966 + 1;
      _M0L2__S966 = _M0L6_2atmpS2321;
      continue;
      joinlet_3405:;
    } else {
      moonbit_decref(_M0L5linesS956);
    }
    break;
  }
  _M0L8_2afieldS3004 = _M0L14current__titleS959->$0;
  _M0L6_2acntS3242
  = Moonbit_rc_count(Moonbit_object_header(_M0L14current__titleS959));
  if (_M0L6_2acntS3242 > 1) {
    int32_t _M0L11_2anew__cntS3243 = _M0L6_2acntS3242 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L14current__titleS959), _M0L11_2anew__cntS3243);
    if (_M0L8_2afieldS3004) {
      moonbit_incref(_M0L8_2afieldS3004);
    }
  } else if (_M0L6_2acntS3242 == 1) {
    #line 66 "/home/developer/Documents/1/fasta_io.mbt"
    moonbit_free(_M0L14current__titleS959);
  }
  _M0L7_2abindS982 = _M0L8_2afieldS3004;
  if (_M0L7_2abindS982 == 0) {
    if (_M0L7_2abindS982) {
      moonbit_decref(_M0L7_2abindS982);
    }
    moonbit_decref(_M0L12current__seqS960);
  } else {
    moonbit_string_t _M0L7_2aSomeS983 = _M0L7_2abindS982;
    moonbit_string_t _M0L8_2atitleS984 = _M0L7_2aSomeS983;
    _M0L5titleS981 = _M0L8_2atitleS984;
    goto join_980;
  }
  goto joinlet_3410;
  join_980:;
  _M0L8_2afieldS3003 = _M0L12current__seqS960->$0;
  _M0L6_2acntS3240
  = Moonbit_rc_count(Moonbit_object_header(_M0L12current__seqS960));
  if (_M0L6_2acntS3240 > 1) {
    int32_t _M0L11_2anew__cntS3241 = _M0L6_2acntS3240 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L12current__seqS960), _M0L11_2anew__cntS3241);
    moonbit_incref(_M0L8_2afieldS3003);
  } else if (_M0L6_2acntS3240 == 1) {
    #line 68 "/home/developer/Documents/1/fasta_io.mbt"
    moonbit_free(_M0L12current__seqS960);
  }
  _M0L3valS2344 = _M0L8_2afieldS3003;
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2343 = _M0MPB13StringBuilder10to__string(_M0L3valS2344);
  moonbit_decref(_M0L3valS2344);
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2342
  = _M0FP26biolab8bio__seq24fasta__title__to__record(_M0L5titleS981, _M0L6_2atmpS2343);
  moonbit_decref(_M0L5titleS981);
  moonbit_decref(_M0L6_2atmpS2343);
  #line 68 "/home/developer/Documents/1/fasta_io.mbt"
  _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(_M0L7recordsS955, _M0L6_2atmpS2342);
  moonbit_decref(_M0L6_2atmpS2342);
  joinlet_3410:;
  return _M0L7recordsS955;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0FP26biolab8bio__seq24fasta__title__to__record(
  moonbit_string_t _M0L5titleS950,
  moonbit_string_t _M0L8seq__strS946
) {
  moonbit_string_t _M0L2idS944;
  moonbit_string_t _M0L11descriptionS945;
  int32_t _M0L3posS949;
  moonbit_string_t _M0L7_2abindS952;
  int32_t _M0L6_2atmpS2320;
  struct _M0TPC16string10StringView _M0L6_2atmpS2319;
  int64_t _M0L7_2abindS951;
  int64_t _M0L6_2atmpS2318;
  struct _M0TPC16string10StringView _M0L6_2atmpS2317;
  moonbit_string_t _M0L6_2atmpS2316;
  struct _M0TP26biolab8bio__seq3Seq* _M0L6_2atmpS2312;
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS947;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2315;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2314;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2313;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_3413;
  #line 78 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7_2abindS952 = (moonbit_string_t)moonbit_string_literal_65.data;
  _M0L6_2atmpS2320 = Moonbit_array_length(_M0L7_2abindS952);
  _M0L6_2atmpS2319
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L7_2abindS952, .$1 = 0, .$2 = _M0L6_2atmpS2320
  };
  #line 79 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L7_2abindS951
  = _M0MPC16string6String4find(_M0L5titleS950, _M0L6_2atmpS2319);
  moonbit_decref(_M0L6_2atmpS2319.$0);
  if (_M0L7_2abindS951 == 4294967296ll) {
    moonbit_incref(_M0L5titleS950);
    moonbit_incref(_M0L5titleS950);
    _M0L2idS944 = _M0L5titleS950;
    _M0L11descriptionS945 = _M0L5titleS950;
    goto join_943;
  } else {
    int64_t _M0L7_2aSomeS953 = _M0L7_2abindS951;
    int32_t _M0L6_2aposS954 = (int32_t)_M0L7_2aSomeS953;
    _M0L3posS949 = _M0L6_2aposS954;
    goto join_948;
  }
  join_948:;
  _M0L6_2atmpS2318 = (int64_t)_M0L3posS949;
  #line 80 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2317
  = _M0MPC16string6String11sub_2einner(_M0L5titleS950, 0, _M0L6_2atmpS2318);
  #line 80 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2316 = _M0MPC16string10StringView9to__owned(_M0L6_2atmpS2317);
  moonbit_decref(_M0L6_2atmpS2317.$0);
  moonbit_incref(_M0L5titleS950);
  _M0L2idS944 = _M0L6_2atmpS2316;
  _M0L11descriptionS945 = _M0L5titleS950;
  goto join_943;
  join_943:;
  #line 84 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2312 = _M0MP26biolab8bio__seq3Seq3new(_M0L8seq__strS946);
  _M0L7_2abindS947 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2315 = _M0L7_2abindS947;
  _M0L6_2atmpS2314
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2315, .$1 = 0, .$2 = 0
  };
  #line 88 "/home/developer/Documents/1/fasta_io.mbt"
  _M0L6_2atmpS2313 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2314, 0ll);
  moonbit_decref(_M0L6_2atmpS2314.$0);
  moonbit_incref(_M0L2idS944);
  _block_3413
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_3413)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _block_3413->$0 = _M0L6_2atmpS2312;
  _block_3413->$1 = _M0L2idS944;
  _block_3413->$2 = _M0L2idS944;
  _block_3413->$3 = _M0L11descriptionS945;
  _block_3413->$4 = _M0L6_2atmpS2313;
  return _block_3413;
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MP26biolab8bio__seq9SeqRecord11new_2einner(
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS938,
  moonbit_string_t _M0L2idS939,
  moonbit_string_t _M0L4nameS940,
  moonbit_string_t _M0L11descriptionS941
) {
  struct _M0TUsRPB5ArrayGiEE** _M0L7_2abindS942;
  struct _M0TUsRPB5ArrayGiEE** _M0L6_2atmpS2311;
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L6_2atmpS2310;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L6_2atmpS2309;
  struct _M0TP26biolab8bio__seq9SeqRecord* _block_3414;
  #line 31 "/home/developer/Documents/1/seq_record.mbt"
  _M0L7_2abindS942 = (struct _M0TUsRPB5ArrayGiEE**)moonbit_empty_ref_array;
  _M0L6_2atmpS2311 = _M0L7_2abindS942;
  _M0L6_2atmpS2310
  = (struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE){
    .$0 = _M0L6_2atmpS2311, .$1 = 0, .$2 = 0
  };
  #line 42 "/home/developer/Documents/1/seq_record.mbt"
  _M0L6_2atmpS2309 = _M0MPB3Map3MapGsRPB5ArrayGiEE(_M0L6_2atmpS2310, 0ll);
  moonbit_decref(_M0L6_2atmpS2310.$0);
  moonbit_incref(_M0L3seqS938);
  moonbit_incref(_M0L2idS939);
  moonbit_incref(_M0L4nameS940);
  moonbit_incref(_M0L11descriptionS941);
  _block_3414
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq9SeqRecord));
  Moonbit_object_header(_block_3414)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 27, 0);
  _block_3414->$0 = _M0L3seqS938;
  _block_3414->$1 = _M0L2idS939;
  _block_3414->$2 = _M0L4nameS940;
  _block_3414->$3 = _M0L11descriptionS941;
  _block_3414->$4 = _M0L6_2atmpS2309;
  return _block_3414;
}

int32_t _M0MP26biolab8bio__seq9SeqRecord6length(
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L4selfS937
) {
  struct _M0TP26biolab8bio__seq3Seq* _M0L3seqS2308;
  int32_t _result_3415;
  #line 84 "/home/developer/Documents/1/seq_record.mbt"
  _M0L3seqS2308 = _M0L4selfS937->$0;
  moonbit_incref(_M0L3seqS2308);
  #line 85 "/home/developer/Documents/1/seq_record.mbt"
  _result_3415 = _M0MP26biolab8bio__seq3Seq6length(_M0L3seqS2308);
  moonbit_decref(_M0L3seqS2308);
  return _result_3415;
}

int32_t _M0MP26biolab8bio__seq3Seq9get__char(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS935,
  int32_t _M0L1iS936
) {
  moonbit_string_t _M0L4dataS2307;
  int32_t _result_3416;
  #line 123 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2307 = _M0L4selfS935->$0;
  moonbit_incref(_M0L4dataS2307);
  #line 124 "/home/developer/Documents/1/seq.mbt"
  _result_3416 = _M0MPC16string6String9get__char(_M0L4dataS2307, _M0L1iS936);
  moonbit_decref(_M0L4dataS2307);
  return _result_3416;
}

int32_t _M0MP26biolab8bio__seq3Seq6length(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS934
) {
  moonbit_string_t _M0L4dataS2306;
  int32_t _result_3417;
  #line 86 "/home/developer/Documents/1/seq.mbt"
  _M0L4dataS2306 = _M0L4selfS934->$0;
  moonbit_incref(_M0L4dataS2306);
  #line 87 "/home/developer/Documents/1/seq.mbt"
  _result_3417
  = _M0MPC16string6String20char__length_2einner(_M0L4dataS2306, 0, 4294967296ll);
  moonbit_decref(_M0L4dataS2306);
  return _result_3417;
}

moonbit_string_t _M0MP26biolab8bio__seq3Seq10to__string(
  struct _M0TP26biolab8bio__seq3Seq* _M0L4selfS933
) {
  moonbit_string_t _M0L8_2afieldS3015;
  #line 74 "/home/developer/Documents/1/seq.mbt"
  _M0L8_2afieldS3015 = _M0L4selfS933->$0;
  moonbit_incref(_M0L8_2afieldS3015);
  return _M0L8_2afieldS3015;
}

struct _M0TP26biolab8bio__seq3Seq* _M0MP26biolab8bio__seq3Seq3new(
  moonbit_string_t _M0L1sS932
) {
  struct _M0TP26biolab8bio__seq3Seq* _block_3418;
  #line 62 "/home/developer/Documents/1/seq.mbt"
  moonbit_incref(_M0L1sS932);
  _block_3418
  = (struct _M0TP26biolab8bio__seq3Seq*)moonbit_malloc(sizeof(struct _M0TP26biolab8bio__seq3Seq));
  Moonbit_object_header(_block_3418)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 45, 0);
  _block_3418->$0 = _M0L1sS932;
  return _block_3418;
}

moonbit_string_t _M0MPC15array5Array2atGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS921,
  int32_t _M0L5indexS922
) {
  int32_t _M0L3lenS920;
  int32_t _if__result_3419;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS920 = _M0L4selfS921->$1;
  if (_M0L5indexS922 >= 0) {
    _if__result_3419 = _M0L5indexS922 < _M0L3lenS920;
  } else {
    _if__result_3419 = 0;
  }
  if (_if__result_3419) {
    moonbit_string_t* _M0L6_2atmpS2302;
    moonbit_string_t _M0L6_2atmpS3016;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2302 = _M0MPC15array5Array6bufferGsE(_M0L4selfS921);
    _M0L6_2atmpS3016 = (moonbit_string_t)_M0L6_2atmpS2302[_M0L5indexS922];
    moonbit_incref(_M0L6_2atmpS3016);
    moonbit_decref(_M0L6_2atmpS2302);
    return _M0L6_2atmpS3016;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0MPC15array5Array2atGRP26biolab8bio__seq20MultipleSeqAlignmentE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L4selfS924,
  int32_t _M0L5indexS925
) {
  int32_t _M0L3lenS923;
  int32_t _if__result_3420;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS923 = _M0L4selfS924->$1;
  if (_M0L5indexS925 >= 0) {
    _if__result_3420 = _M0L5indexS925 < _M0L3lenS923;
  } else {
    _if__result_3420 = 0;
  }
  if (_if__result_3420) {
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L6_2atmpS2303;
    struct _M0TP26biolab8bio__seq20MultipleSeqAlignment* _M0L6_2atmpS3017;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2303
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq20MultipleSeqAlignmentE(_M0L4selfS924);
    _M0L6_2atmpS3017
    = (struct _M0TP26biolab8bio__seq20MultipleSeqAlignment*)_M0L6_2atmpS2303[
        _M0L5indexS925
      ];
    if (_M0L6_2atmpS3017) {
      moonbit_incref(_M0L6_2atmpS3017);
    }
    moonbit_decref(_M0L6_2atmpS2303);
    return _M0L6_2atmpS3017;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TP26biolab8bio__seq9SeqRecord* _M0MPC15array5Array2atGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS927,
  int32_t _M0L5indexS928
) {
  int32_t _M0L3lenS926;
  int32_t _if__result_3421;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS926 = _M0L4selfS927->$1;
  if (_M0L5indexS928 >= 0) {
    _if__result_3421 = _M0L5indexS928 < _M0L3lenS926;
  } else {
    _if__result_3421 = 0;
  }
  if (_if__result_3421) {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS2304;
    struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS3018;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2304
    = _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS927);
    _M0L6_2atmpS3018
    = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L6_2atmpS2304[
        _M0L5indexS928
      ];
    if (_M0L6_2atmpS3018) {
      moonbit_incref(_M0L6_2atmpS3018);
    }
    moonbit_decref(_M0L6_2atmpS2304);
    return _M0L6_2atmpS3018;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

struct _M0TPC16string10StringView _M0MPC15array5Array2atGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS930,
  int32_t _M0L5indexS931
) {
  int32_t _M0L3lenS929;
  int32_t _if__result_3422;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS929 = _M0L4selfS930->$1;
  if (_M0L5indexS931 >= 0) {
    _if__result_3422 = _M0L5indexS931 < _M0L3lenS929;
  } else {
    _if__result_3422 = 0;
  }
  if (_if__result_3422) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS2305;
    struct _M0TPC16string10StringView _M0L6_2atmpS3019;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS2305
    = _M0MPC15array5Array6bufferGRPC16string10StringViewE(_M0L4selfS930);
    _M0L6_2atmpS3019 = _M0L6_2atmpS2305[_M0L5indexS931];
    moonbit_incref(_M0L6_2atmpS3019.$0);
    moonbit_decref(_M0L6_2atmpS2305);
    return _M0L6_2atmpS3019;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB7printlnGsE(moonbit_string_t _M0L5inputS919) {
  moonbit_string_t _M0L6_2atmpS2301;
  #line 36 "/home/developer/.moon/lib/core/builtin/console.mbt"
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  _M0L6_2atmpS2301 = _M0IPC16string6StringPB4Show10to__string(_M0L5inputS919);
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  moonbit_println(_M0L6_2atmpS2301);
  moonbit_decref(_M0L6_2atmpS2301);
  return 0;
}

int32_t _M0IPC16string6StringPB4Hash4hash(moonbit_string_t _M0L4selfS915) {
  int32_t _tmp_3423;
  uint32_t _M0L6_2atmpS2300;
  uint32_t _M0Lm3accS913;
  int32_t _M0L7_2abindS914;
  int32_t _M0L1iS916;
  uint32_t _M0L6_2atmpS2299;
  #line 522 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _tmp_3423 = 0;
  _M0L6_2atmpS2300 = *(uint32_t*)&_tmp_3423;
  _M0Lm3accS913 = _M0L6_2atmpS2300 + 374761393u;
  _M0L7_2abindS914 = Moonbit_array_length(_M0L4selfS915);
  _M0L1iS916 = 0;
  while (1) {
    if (_M0L1iS916 < _M0L7_2abindS914) {
      uint32_t _M0L6_2atmpS2294 = _M0Lm3accS913;
      int32_t _M0L6_2atmpS2297;
      int32_t _M0L6_2atmpS2296;
      uint32_t _M0L1vS917;
      uint32_t _M0L6_2atmpS2295;
      int32_t _M0L6_2atmpS2298;
      _M0Lm3accS913 = _M0L6_2atmpS2294 + 4u;
      _M0L6_2atmpS2297 = _M0L4selfS915[_M0L1iS916];
      _M0L6_2atmpS2296 = (int32_t)_M0L6_2atmpS2297;
      _M0L1vS917 = *(uint32_t*)&_M0L6_2atmpS2296;
      _M0L6_2atmpS2295 = _M0Lm3accS913;
      #line 527 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
      _M0Lm3accS913 = _M0FPB13consume4__acc(_M0L6_2atmpS2295, _M0L1vS917);
      _M0L6_2atmpS2298 = _M0L1iS916 + 1;
      _M0L1iS916 = _M0L6_2atmpS2298;
      continue;
    }
    break;
  }
  _M0L6_2atmpS2299 = _M0Lm3accS913;
  #line 529 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  return _M0FPB13finalize__acc(_M0L6_2atmpS2299);
}

struct _M0TUsRPB13StringBuilderE* _M0MPB5Iter24nextGsRPB13StringBuilderE(
  struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0L4selfS912
) {
  #line 1099 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  #line 1100 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  return _M0MPB4Iter4nextGUsRPB13StringBuilderEE(_M0L4selfS912);
}

struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0MPB3Map5iter2GsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS911
) {
  #line 725 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 727 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  return _M0MPB3Map4iterGsRPB13StringBuilderE(_M0L4selfS911);
}

struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0MPB3Map4iterGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS901
) {
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L4headS2293;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE* _M0L11curr__entryS900;
  int32_t _M0L3lenS902;
  struct _M0TPB8MutLocalGiE* _M0L9remainingS903;
  struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__* _closure_3425;
  struct _M0TWEOUsRPB13StringBuilderE* _M0L6_2atmpS2284;
  int64_t _M0L6_2atmpS2285;
  struct _M0TPB4IterGUsRPB13StringBuilderEE* _result_3426;
  #line 705 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4headS2293 = _M0L4selfS901->$5;
  if (_M0L4headS2293) {
    moonbit_incref(_M0L4headS2293);
  }
  _M0L11curr__entryS900
  = (struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE));
  Moonbit_object_header(_M0L11curr__entryS900)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 48, 0);
  _M0L11curr__entryS900->$0 = _M0L4headS2293;
  _M0L3lenS902 = _M0L4selfS901->$1;
  _M0L9remainingS903
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L9remainingS903)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L9remainingS903->$0 = _M0L3lenS902;
  _closure_3425
  = (struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__*)moonbit_malloc(sizeof(struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__));
  Moonbit_object_header(_closure_3425)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 51, 0);
  _closure_3425->code = &_M0MPB3Map4iterGsRPB13StringBuilderEC2286l711;
  _closure_3425->$0 = _M0L9remainingS903;
  _closure_3425->$1 = _M0L11curr__entryS900;
  _M0L6_2atmpS2284 = (struct _M0TWEOUsRPB13StringBuilderE*)_closure_3425;
  _M0L6_2atmpS2285 = (int64_t)_M0L3lenS902;
  #line 710 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _result_3426
  = _M0MPB4Iter3newGUsRPB13StringBuilderEE(_M0L6_2atmpS2284, _M0L6_2atmpS2285);
  moonbit_decref(_M0L6_2atmpS2284);
  return _result_3426;
}

struct _M0TUsRPB13StringBuilderE* _M0MPB3Map4iterGsRPB13StringBuilderEC2286l711(
  struct _M0TWEOUsRPB13StringBuilderE* _M0L6_2aenvS2287
) {
  struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__* _M0L14_2acasted__envS2288;
  struct _M0TPB8MutLocalGORPB5EntryGsRPB13StringBuilderEE* _M0L11curr__entryS900;
  struct _M0TPB8MutLocalGiE* _M0L9remainingS903;
  int32_t _M0L3valS2289;
  #line 711 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14_2acasted__envS2288
  = (struct _M0R103Map_3a_3aiter_7c_5bString_2c_20moonbitlang_2fcore_2fbuiltin_2fStringBuilder_5d_7c_2eanon__u2286__l711__*)_M0L6_2aenvS2287;
  _M0L11curr__entryS900 = _M0L14_2acasted__envS2288->$1;
  _M0L9remainingS903 = _M0L14_2acasted__envS2288->$0;
  _M0L3valS2289 = _M0L9remainingS903->$0;
  if (_M0L3valS2289 > 0) {
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS905 =
      _M0L11curr__entryS900->$0;
    if (_M0L7_2abindS905 == 0) {
      goto join_904;
    } else {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS906 =
        _M0L7_2abindS905;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L4_2axS907 =
        _M0L7_2aSomeS906;
      moonbit_string_t _M0L6_2akeyS908 = _M0L4_2axS907->$4;
      struct _M0TPB13StringBuilder* _M0L8_2avalueS909 = _M0L4_2axS907->$5;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2anextS910 =
        _M0L4_2axS907->$1;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3021 =
        _M0L11curr__entryS900->$0;
      int32_t _M0L3valS2291;
      int32_t _M0L6_2atmpS2290;
      struct _M0TUsRPB13StringBuilderE* _M0L8_2atupleS2292;
      if (_M0L7_2anextS910) {
        moonbit_incref(_M0L7_2anextS910);
      }
      moonbit_incref(_M0L8_2avalueS909);
      moonbit_incref(_M0L6_2akeyS908);
      if (_M0L6_2aoldS3021) {
        moonbit_decref(_M0L6_2aoldS3021);
      }
      _M0L11curr__entryS900->$0 = _M0L7_2anextS910;
      _M0L3valS2291 = _M0L9remainingS903->$0;
      _M0L6_2atmpS2290 = _M0L3valS2291 - 1;
      _M0L9remainingS903->$0 = _M0L6_2atmpS2290;
      _M0L8_2atupleS2292
      = (struct _M0TUsRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TUsRPB13StringBuilderE));
      Moonbit_object_header(_M0L8_2atupleS2292)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 55, 0);
      _M0L8_2atupleS2292->$0 = _M0L6_2akeyS908;
      _M0L8_2atupleS2292->$1 = _M0L8_2avalueS909;
      return _M0L8_2atupleS2292;
    }
  } else {
    goto join_904;
  }
  join_904:;
  return 0;
}

struct _M0TPB13StringBuilder* _M0MPB3Map3getGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS896,
  moonbit_string_t _M0L3keyS892
) {
  int32_t _M0L4hashS891;
  int32_t _M0L14capacity__maskS2283;
  int32_t _M0L6_2atmpS2282;
  int32_t _M0L1iS893;
  int32_t _M0L3idxS894;
  #line 236 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 237 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS891 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS892);
  _M0L14capacity__maskS2283 = _M0L4selfS896->$3;
  _M0L6_2atmpS2282 = _M0L4hashS891 & _M0L14capacity__maskS2283;
  _M0L1iS893 = 0;
  _M0L3idxS894 = _M0L6_2atmpS2282;
  while (1) {
    struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2281 =
      _M0L4selfS896->$0;
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS895;
    if (
      _M0L3idxS894 < 0
      || _M0L3idxS894 >= Moonbit_array_length(_M0L7entriesS2281)
    ) {
      #line 239 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS895
    = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2281[
        _M0L3idxS894
      ];
    if (_M0L7_2abindS895 == 0) {
      struct _M0TPB13StringBuilder* _M0L6_2atmpS2270 = 0;
      return _M0L6_2atmpS2270;
    } else {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS897 =
        _M0L7_2abindS895;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L8_2aentryS898 =
        _M0L7_2aSomeS897;
      int32_t _M0L4hashS2272 = _M0L8_2aentryS898->$3;
      int32_t _if__result_3429;
      int32_t _M0L3pslS2275;
      int32_t _M0L6_2atmpS2277;
      int32_t _M0L6_2atmpS2279;
      int32_t _M0L14capacity__maskS2280;
      int32_t _M0L6_2atmpS2278;
      if (_M0L4hashS2272 == _M0L4hashS891) {
        moonbit_string_t _M0L3keyS2271 = _M0L8_2aentryS898->$4;
        #line 240 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3429
        = _M0L3keyS2271 == _M0L3keyS892
          || Moonbit_array_length(_M0L3keyS2271)
             == Moonbit_array_length(_M0L3keyS892)
             && 0
                == memcmp(_M0L3keyS2271, _M0L3keyS892, Moonbit_array_length(_M0L3keyS2271) * 2);
      } else {
        _if__result_3429 = 0;
      }
      if (_if__result_3429) {
        struct _M0TPB13StringBuilder* _M0L5valueS2274 = _M0L8_2aentryS898->$5;
        struct _M0TPB13StringBuilder* _M0L6_2atmpS2273;
        moonbit_incref(_M0L5valueS2274);
        _M0L6_2atmpS2273 = _M0L5valueS2274;
        return _M0L6_2atmpS2273;
      } else {
        moonbit_incref(_M0L8_2aentryS898);
      }
      _M0L3pslS2275 = _M0L8_2aentryS898->$2;
      moonbit_decref(_M0L8_2aentryS898);
      if (_M0L1iS893 > _M0L3pslS2275) {
        struct _M0TPB13StringBuilder* _M0L6_2atmpS2276 = 0;
        return _M0L6_2atmpS2276;
      }
      _M0L6_2atmpS2277 = _M0L1iS893 + 1;
      _M0L6_2atmpS2279 = _M0L3idxS894 + 1;
      _M0L14capacity__maskS2280 = _M0L4selfS896->$3;
      _M0L6_2atmpS2278 = _M0L6_2atmpS2279 & _M0L14capacity__maskS2280;
      _M0L1iS893 = _M0L6_2atmpS2277;
      _M0L3idxS894 = _M0L6_2atmpS2278;
      continue;
    }
    break;
  }
}

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0MPB3Map3MapGsRPB5ArrayGiEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGiEEE _M0L3arrS848,
  int64_t _M0L8capacityS850
) {
  int32_t _M0L3endS2235;
  int32_t _M0L5startS2236;
  int32_t _M0L6lengthS847;
  int32_t _M0L8capacityS849;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L1mS853;
  int32_t _M0L3endS2232;
  int32_t _M0L5startS2233;
  int32_t _M0L7_2abindS854;
  int32_t _M0L2__S855;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2235 = _M0L3arrS848.$2;
  _M0L5startS2236 = _M0L3arrS848.$1;
  _M0L6lengthS847 = _M0L3endS2235 - _M0L5startS2236;
  if (_M0L8capacityS850 == 4294967296ll) {
    if (_M0L6lengthS847 == 0) {
      _M0L8capacityS849 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS849 = _M0FPB21capacity__for__length(_M0L6lengthS847);
    }
  } else {
    int64_t _M0L7_2aSomeS851 = _M0L8capacityS850;
    int32_t _M0L11_2acapacityS852 = (int32_t)_M0L7_2aSomeS851;
    int32_t _M0L6_2atmpS2234;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2234 = _M0FPB21capacity__for__length(_M0L6lengthS847);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS849
    = _M0MPC13int3Int3max(_M0L11_2acapacityS852, _M0L6_2atmpS2234);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS853 = _M0FPB8new__mapGsRPB5ArrayGiEE(_M0L8capacityS849);
  _M0L3endS2232 = _M0L3arrS848.$2;
  _M0L5startS2233 = _M0L3arrS848.$1;
  _M0L7_2abindS854 = _M0L3endS2232 - _M0L5startS2233;
  _M0L2__S855 = 0;
  while (1) {
    if (_M0L2__S855 < _M0L7_2abindS854) {
      struct _M0TUsRPB5ArrayGiEE** _M0L3bufS2229 = _M0L3arrS848.$0;
      int32_t _M0L5startS2231 = _M0L3arrS848.$1;
      int32_t _M0L6_2atmpS2230 = _M0L5startS2231 + _M0L2__S855;
      struct _M0TUsRPB5ArrayGiEE* _M0L1eS856 =
        (struct _M0TUsRPB5ArrayGiEE*)_M0L3bufS2229[_M0L6_2atmpS2230];
      moonbit_string_t _M0L6_2atmpS2226 = _M0L1eS856->$0;
      struct _M0TPB5ArrayGiE* _M0L6_2atmpS2227 = _M0L1eS856->$1;
      int32_t _M0L6_2atmpS2228;
      moonbit_incref(_M0L6_2atmpS2227);
      moonbit_incref(_M0L6_2atmpS2226);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB5ArrayGiEE(_M0L1mS853, _M0L6_2atmpS2226, _M0L6_2atmpS2227);
      moonbit_decref(_M0L6_2atmpS2226);
      moonbit_decref(_M0L6_2atmpS2227);
      _M0L6_2atmpS2228 = _M0L2__S855 + 1;
      _M0L2__S855 = _M0L6_2atmpS2228;
      continue;
    }
    break;
  }
  return _M0L1mS853;
}

struct _M0TPB3MapGssE* _M0MPB3Map3MapGssE(
  struct _M0TPB9ArrayViewGUssEE _M0L3arrS859,
  int64_t _M0L8capacityS861
) {
  int32_t _M0L3endS2246;
  int32_t _M0L5startS2247;
  int32_t _M0L6lengthS858;
  int32_t _M0L8capacityS860;
  struct _M0TPB3MapGssE* _M0L1mS864;
  int32_t _M0L3endS2243;
  int32_t _M0L5startS2244;
  int32_t _M0L7_2abindS865;
  int32_t _M0L2__S866;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2246 = _M0L3arrS859.$2;
  _M0L5startS2247 = _M0L3arrS859.$1;
  _M0L6lengthS858 = _M0L3endS2246 - _M0L5startS2247;
  if (_M0L8capacityS861 == 4294967296ll) {
    if (_M0L6lengthS858 == 0) {
      _M0L8capacityS860 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS860 = _M0FPB21capacity__for__length(_M0L6lengthS858);
    }
  } else {
    int64_t _M0L7_2aSomeS862 = _M0L8capacityS861;
    int32_t _M0L11_2acapacityS863 = (int32_t)_M0L7_2aSomeS862;
    int32_t _M0L6_2atmpS2245;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2245 = _M0FPB21capacity__for__length(_M0L6lengthS858);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS860
    = _M0MPC13int3Int3max(_M0L11_2acapacityS863, _M0L6_2atmpS2245);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS864 = _M0FPB8new__mapGssE(_M0L8capacityS860);
  _M0L3endS2243 = _M0L3arrS859.$2;
  _M0L5startS2244 = _M0L3arrS859.$1;
  _M0L7_2abindS865 = _M0L3endS2243 - _M0L5startS2244;
  _M0L2__S866 = 0;
  while (1) {
    if (_M0L2__S866 < _M0L7_2abindS865) {
      struct _M0TUssE** _M0L3bufS2240 = _M0L3arrS859.$0;
      int32_t _M0L5startS2242 = _M0L3arrS859.$1;
      int32_t _M0L6_2atmpS2241 = _M0L5startS2242 + _M0L2__S866;
      struct _M0TUssE* _M0L1eS867 =
        (struct _M0TUssE*)_M0L3bufS2240[_M0L6_2atmpS2241];
      moonbit_string_t _M0L6_2atmpS2237 = _M0L1eS867->$0;
      moonbit_string_t _M0L6_2atmpS2238 = _M0L1eS867->$1;
      int32_t _M0L6_2atmpS2239;
      moonbit_incref(_M0L6_2atmpS2238);
      moonbit_incref(_M0L6_2atmpS2237);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGssE(_M0L1mS864, _M0L6_2atmpS2237, _M0L6_2atmpS2238);
      moonbit_decref(_M0L6_2atmpS2237);
      moonbit_decref(_M0L6_2atmpS2238);
      _M0L6_2atmpS2239 = _M0L2__S866 + 1;
      _M0L2__S866 = _M0L6_2atmpS2239;
      continue;
    }
    break;
  }
  return _M0L1mS864;
}

struct _M0TPB3MapGsRPB5ArrayGsEE* _M0MPB3Map3MapGsRPB5ArrayGsEE(
  struct _M0TPB9ArrayViewGUsRPB5ArrayGsEEE _M0L3arrS870,
  int64_t _M0L8capacityS872
) {
  int32_t _M0L3endS2257;
  int32_t _M0L5startS2258;
  int32_t _M0L6lengthS869;
  int32_t _M0L8capacityS871;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L1mS875;
  int32_t _M0L3endS2254;
  int32_t _M0L5startS2255;
  int32_t _M0L7_2abindS876;
  int32_t _M0L2__S877;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2257 = _M0L3arrS870.$2;
  _M0L5startS2258 = _M0L3arrS870.$1;
  _M0L6lengthS869 = _M0L3endS2257 - _M0L5startS2258;
  if (_M0L8capacityS872 == 4294967296ll) {
    if (_M0L6lengthS869 == 0) {
      _M0L8capacityS871 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS871 = _M0FPB21capacity__for__length(_M0L6lengthS869);
    }
  } else {
    int64_t _M0L7_2aSomeS873 = _M0L8capacityS872;
    int32_t _M0L11_2acapacityS874 = (int32_t)_M0L7_2aSomeS873;
    int32_t _M0L6_2atmpS2256;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2256 = _M0FPB21capacity__for__length(_M0L6lengthS869);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS871
    = _M0MPC13int3Int3max(_M0L11_2acapacityS874, _M0L6_2atmpS2256);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS875 = _M0FPB8new__mapGsRPB5ArrayGsEE(_M0L8capacityS871);
  _M0L3endS2254 = _M0L3arrS870.$2;
  _M0L5startS2255 = _M0L3arrS870.$1;
  _M0L7_2abindS876 = _M0L3endS2254 - _M0L5startS2255;
  _M0L2__S877 = 0;
  while (1) {
    if (_M0L2__S877 < _M0L7_2abindS876) {
      struct _M0TUsRPB5ArrayGsEE** _M0L3bufS2251 = _M0L3arrS870.$0;
      int32_t _M0L5startS2253 = _M0L3arrS870.$1;
      int32_t _M0L6_2atmpS2252 = _M0L5startS2253 + _M0L2__S877;
      struct _M0TUsRPB5ArrayGsEE* _M0L1eS878 =
        (struct _M0TUsRPB5ArrayGsEE*)_M0L3bufS2251[_M0L6_2atmpS2252];
      moonbit_string_t _M0L6_2atmpS2248 = _M0L1eS878->$0;
      struct _M0TPB5ArrayGsE* _M0L6_2atmpS2249 = _M0L1eS878->$1;
      int32_t _M0L6_2atmpS2250;
      moonbit_incref(_M0L6_2atmpS2249);
      moonbit_incref(_M0L6_2atmpS2248);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB5ArrayGsEE(_M0L1mS875, _M0L6_2atmpS2248, _M0L6_2atmpS2249);
      moonbit_decref(_M0L6_2atmpS2248);
      moonbit_decref(_M0L6_2atmpS2249);
      _M0L6_2atmpS2250 = _M0L2__S877 + 1;
      _M0L2__S877 = _M0L6_2atmpS2250;
      continue;
    }
    break;
  }
  return _M0L1mS875;
}

struct _M0TPB3MapGsRPB13StringBuilderE* _M0MPB3Map3MapGsRPB13StringBuilderE(
  struct _M0TPB9ArrayViewGUsRPB13StringBuilderEE _M0L3arrS881,
  int64_t _M0L8capacityS883
) {
  int32_t _M0L3endS2268;
  int32_t _M0L5startS2269;
  int32_t _M0L6lengthS880;
  int32_t _M0L8capacityS882;
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L1mS886;
  int32_t _M0L3endS2265;
  int32_t _M0L5startS2266;
  int32_t _M0L7_2abindS887;
  int32_t _M0L2__S888;
  #line 83 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3endS2268 = _M0L3arrS881.$2;
  _M0L5startS2269 = _M0L3arrS881.$1;
  _M0L6lengthS880 = _M0L3endS2268 - _M0L5startS2269;
  if (_M0L8capacityS883 == 4294967296ll) {
    if (_M0L6lengthS880 == 0) {
      _M0L8capacityS882 = 8;
    } else {
      #line 95 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L8capacityS882 = _M0FPB21capacity__for__length(_M0L6lengthS880);
    }
  } else {
    int64_t _M0L7_2aSomeS884 = _M0L8capacityS883;
    int32_t _M0L11_2acapacityS885 = (int32_t)_M0L7_2aSomeS884;
    int32_t _M0L6_2atmpS2267;
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L6_2atmpS2267 = _M0FPB21capacity__for__length(_M0L6lengthS880);
    #line 90 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    _M0L8capacityS882
    = _M0MPC13int3Int3max(_M0L11_2acapacityS885, _M0L6_2atmpS2267);
  }
  #line 98 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L1mS886 = _M0FPB8new__mapGsRPB13StringBuilderE(_M0L8capacityS882);
  _M0L3endS2265 = _M0L3arrS881.$2;
  _M0L5startS2266 = _M0L3arrS881.$1;
  _M0L7_2abindS887 = _M0L3endS2265 - _M0L5startS2266;
  _M0L2__S888 = 0;
  while (1) {
    if (_M0L2__S888 < _M0L7_2abindS887) {
      struct _M0TUsRPB13StringBuilderE** _M0L3bufS2262 = _M0L3arrS881.$0;
      int32_t _M0L5startS2264 = _M0L3arrS881.$1;
      int32_t _M0L6_2atmpS2263 = _M0L5startS2264 + _M0L2__S888;
      struct _M0TUsRPB13StringBuilderE* _M0L1eS889 =
        (struct _M0TUsRPB13StringBuilderE*)_M0L3bufS2262[_M0L6_2atmpS2263];
      moonbit_string_t _M0L6_2atmpS2259 = _M0L1eS889->$0;
      struct _M0TPB13StringBuilder* _M0L6_2atmpS2260 = _M0L1eS889->$1;
      int32_t _M0L6_2atmpS2261;
      moonbit_incref(_M0L6_2atmpS2260);
      moonbit_incref(_M0L6_2atmpS2259);
      #line 100 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map3setGsRPB13StringBuilderE(_M0L1mS886, _M0L6_2atmpS2259, _M0L6_2atmpS2260);
      moonbit_decref(_M0L6_2atmpS2259);
      moonbit_decref(_M0L6_2atmpS2260);
      _M0L6_2atmpS2261 = _M0L2__S888 + 1;
      _M0L2__S888 = _M0L6_2atmpS2261;
      continue;
    }
    break;
  }
  return _M0L1mS886;
}

int32_t _M0MPB3Map3setGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS835,
  moonbit_string_t _M0L3keyS836,
  struct _M0TPB5ArrayGiE* _M0L5valueS837
) {
  int32_t _M0L6_2atmpS2222;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2222 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS836);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(_M0L4selfS835, _M0L3keyS836, _M0L5valueS837, _M0L6_2atmpS2222);
  return 0;
}

int32_t _M0MPB3Map3setGssE(
  struct _M0TPB3MapGssE* _M0L4selfS838,
  moonbit_string_t _M0L3keyS839,
  moonbit_string_t _M0L5valueS840
) {
  int32_t _M0L6_2atmpS2223;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2223 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS839);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGssE(_M0L4selfS838, _M0L3keyS839, _M0L5valueS840, _M0L6_2atmpS2223);
  return 0;
}

int32_t _M0MPB3Map3setGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS841,
  moonbit_string_t _M0L3keyS842,
  struct _M0TPB5ArrayGsE* _M0L5valueS843
) {
  int32_t _M0L6_2atmpS2224;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2224 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS842);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB5ArrayGsEE(_M0L4selfS841, _M0L3keyS842, _M0L5valueS843, _M0L6_2atmpS2224);
  return 0;
}

int32_t _M0MPB3Map3setGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS844,
  moonbit_string_t _M0L3keyS845,
  struct _M0TPB13StringBuilder* _M0L5valueS846
) {
  int32_t _M0L6_2atmpS2225;
  #line 127 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2225 = _M0IPC16string6StringPB4Hash4hash(_M0L3keyS845);
  #line 129 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0MPB3Map15set__with__hashGsRPB13StringBuilderE(_M0L4selfS844, _M0L3keyS845, _M0L5valueS846, _M0L6_2atmpS2225);
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS774,
  moonbit_string_t _M0L3keyS780,
  struct _M0TPB5ArrayGiE* _M0L5valueS781,
  int32_t _M0L4hashS776
) {
  int32_t _M0L14capacity__maskS2167;
  int32_t _M0L6_2atmpS2166;
  int32_t _M0L3pslS771;
  int32_t _M0L3idxS772;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2167 = _M0L4selfS774->$3;
  _M0L6_2atmpS2166 = _M0L4hashS776 & _M0L14capacity__maskS2167;
  _M0L3pslS771 = 0;
  _M0L3idxS772 = _M0L6_2atmpS2166;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2165 =
      _M0L4selfS774->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS773;
    if (
      _M0L3idxS772 < 0
      || _M0L3idxS772 >= Moonbit_array_length(_M0L7entriesS2165)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS773
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2165[_M0L3idxS772];
    if (_M0L7_2abindS773 == 0) {
      int32_t _M0L4sizeS2150 = _M0L4selfS774->$1;
      int32_t _M0L8grow__atS2151 = _M0L4selfS774->$4;
      int32_t _M0L7_2abindS777;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS778;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS779;
      if (_M0L4sizeS2150 >= _M0L8grow__atS2151) {
        int32_t _M0L14capacity__maskS2153;
        int32_t _M0L6_2atmpS2152;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB5ArrayGiEE(_M0L4selfS774);
        _M0L14capacity__maskS2153 = _M0L4selfS774->$3;
        _M0L6_2atmpS2152 = _M0L4hashS776 & _M0L14capacity__maskS2153;
        _M0L3pslS771 = 0;
        _M0L3idxS772 = _M0L6_2atmpS2152;
        continue;
      }
      _M0L7_2abindS777 = _M0L4selfS774->$6;
      _M0L7_2abindS778 = 0;
      moonbit_incref(_M0L3keyS780);
      moonbit_incref(_M0L5valueS781);
      _M0L5entryS779
      = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE));
      Moonbit_object_header(_M0L5entryS779)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 59, 0);
      _M0L5entryS779->$0 = _M0L7_2abindS777;
      _M0L5entryS779->$1 = _M0L7_2abindS778;
      _M0L5entryS779->$2 = _M0L3pslS771;
      _M0L5entryS779->$3 = _M0L4hashS776;
      _M0L5entryS779->$4 = _M0L3keyS780;
      _M0L5entryS779->$5 = _M0L5valueS781;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS774, _M0L3idxS772, _M0L5entryS779);
      moonbit_decref(_M0L5entryS779);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS782 = _M0L7_2abindS773;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L14_2acurr__entryS783 =
        _M0L7_2aSomeS782;
      int32_t _M0L4hashS2155 = _M0L14_2acurr__entryS783->$3;
      int32_t _if__result_3435;
      int32_t _M0L3pslS2156;
      int32_t _M0L6_2atmpS2161;
      int32_t _M0L6_2atmpS2163;
      int32_t _M0L14capacity__maskS2164;
      int32_t _M0L6_2atmpS2162;
      if (_M0L4hashS2155 == _M0L4hashS776) {
        moonbit_string_t _M0L3keyS2154 = _M0L14_2acurr__entryS783->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3435
        = _M0L3keyS2154 == _M0L3keyS780
          || Moonbit_array_length(_M0L3keyS2154)
             == Moonbit_array_length(_M0L3keyS780)
             && 0
                == memcmp(_M0L3keyS2154, _M0L3keyS780, Moonbit_array_length(_M0L3keyS2154) * 2);
      } else {
        _if__result_3435 = 0;
      }
      if (_if__result_3435) {
        struct _M0TPB5ArrayGiE* _M0L6_2aoldS3046 =
          _M0L14_2acurr__entryS783->$5;
        moonbit_incref(_M0L5valueS781);
        moonbit_decref(_M0L6_2aoldS3046);
        _M0L14_2acurr__entryS783->$5 = _M0L5valueS781;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS783);
      }
      _M0L3pslS2156 = _M0L14_2acurr__entryS783->$2;
      if (_M0L3pslS771 > _M0L3pslS2156) {
        int32_t _M0L4sizeS2157 = _M0L4selfS774->$1;
        int32_t _M0L8grow__atS2158 = _M0L4selfS774->$4;
        int32_t _M0L7_2abindS784;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS785;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS786;
        if (_M0L4sizeS2157 >= _M0L8grow__atS2158) {
          int32_t _M0L14capacity__maskS2160;
          int32_t _M0L6_2atmpS2159;
          moonbit_decref(_M0L14_2acurr__entryS783);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB5ArrayGiEE(_M0L4selfS774);
          _M0L14capacity__maskS2160 = _M0L4selfS774->$3;
          _M0L6_2atmpS2159 = _M0L4hashS776 & _M0L14capacity__maskS2160;
          _M0L3pslS771 = 0;
          _M0L3idxS772 = _M0L6_2atmpS2159;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGiEE(_M0L4selfS774, _M0L3idxS772, _M0L14_2acurr__entryS783);
        moonbit_decref(_M0L14_2acurr__entryS783);
        _M0L7_2abindS784 = _M0L4selfS774->$6;
        _M0L7_2abindS785 = 0;
        moonbit_incref(_M0L3keyS780);
        moonbit_incref(_M0L5valueS781);
        _M0L5entryS786
        = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGiEE));
        Moonbit_object_header(_M0L5entryS786)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 59, 0);
        _M0L5entryS786->$0 = _M0L7_2abindS784;
        _M0L5entryS786->$1 = _M0L7_2abindS785;
        _M0L5entryS786->$2 = _M0L3pslS771;
        _M0L5entryS786->$3 = _M0L4hashS776;
        _M0L5entryS786->$4 = _M0L3keyS780;
        _M0L5entryS786->$5 = _M0L5valueS781;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS774, _M0L3idxS772, _M0L5entryS786);
        moonbit_decref(_M0L5entryS786);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS783);
      }
      _M0L6_2atmpS2161 = _M0L3pslS771 + 1;
      _M0L6_2atmpS2163 = _M0L3idxS772 + 1;
      _M0L14capacity__maskS2164 = _M0L4selfS774->$3;
      _M0L6_2atmpS2162 = _M0L6_2atmpS2163 & _M0L14capacity__maskS2164;
      _M0L3pslS771 = _M0L6_2atmpS2161;
      _M0L3idxS772 = _M0L6_2atmpS2162;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGssE(
  struct _M0TPB3MapGssE* _M0L4selfS790,
  moonbit_string_t _M0L3keyS796,
  moonbit_string_t _M0L5valueS797,
  int32_t _M0L4hashS792
) {
  int32_t _M0L14capacity__maskS2185;
  int32_t _M0L6_2atmpS2184;
  int32_t _M0L3pslS787;
  int32_t _M0L3idxS788;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2185 = _M0L4selfS790->$3;
  _M0L6_2atmpS2184 = _M0L4hashS792 & _M0L14capacity__maskS2185;
  _M0L3pslS787 = 0;
  _M0L3idxS788 = _M0L6_2atmpS2184;
  while (1) {
    struct _M0TPB5EntryGssE** _M0L7entriesS2183 = _M0L4selfS790->$0;
    struct _M0TPB5EntryGssE* _M0L7_2abindS789;
    if (
      _M0L3idxS788 < 0
      || _M0L3idxS788 >= Moonbit_array_length(_M0L7entriesS2183)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS789
    = (struct _M0TPB5EntryGssE*)_M0L7entriesS2183[_M0L3idxS788];
    if (_M0L7_2abindS789 == 0) {
      int32_t _M0L4sizeS2168 = _M0L4selfS790->$1;
      int32_t _M0L8grow__atS2169 = _M0L4selfS790->$4;
      int32_t _M0L7_2abindS793;
      struct _M0TPB5EntryGssE* _M0L7_2abindS794;
      struct _M0TPB5EntryGssE* _M0L5entryS795;
      if (_M0L4sizeS2168 >= _M0L8grow__atS2169) {
        int32_t _M0L14capacity__maskS2171;
        int32_t _M0L6_2atmpS2170;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGssE(_M0L4selfS790);
        _M0L14capacity__maskS2171 = _M0L4selfS790->$3;
        _M0L6_2atmpS2170 = _M0L4hashS792 & _M0L14capacity__maskS2171;
        _M0L3pslS787 = 0;
        _M0L3idxS788 = _M0L6_2atmpS2170;
        continue;
      }
      _M0L7_2abindS793 = _M0L4selfS790->$6;
      _M0L7_2abindS794 = 0;
      moonbit_incref(_M0L3keyS796);
      moonbit_incref(_M0L5valueS797);
      _M0L5entryS795
      = (struct _M0TPB5EntryGssE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGssE));
      Moonbit_object_header(_M0L5entryS795)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 64, 0);
      _M0L5entryS795->$0 = _M0L7_2abindS793;
      _M0L5entryS795->$1 = _M0L7_2abindS794;
      _M0L5entryS795->$2 = _M0L3pslS787;
      _M0L5entryS795->$3 = _M0L4hashS792;
      _M0L5entryS795->$4 = _M0L3keyS796;
      _M0L5entryS795->$5 = _M0L5valueS797;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGssE(_M0L4selfS790, _M0L3idxS788, _M0L5entryS795);
      moonbit_decref(_M0L5entryS795);
      return 0;
    } else {
      struct _M0TPB5EntryGssE* _M0L7_2aSomeS798 = _M0L7_2abindS789;
      struct _M0TPB5EntryGssE* _M0L14_2acurr__entryS799 = _M0L7_2aSomeS798;
      int32_t _M0L4hashS2173 = _M0L14_2acurr__entryS799->$3;
      int32_t _if__result_3437;
      int32_t _M0L3pslS2174;
      int32_t _M0L6_2atmpS2179;
      int32_t _M0L6_2atmpS2181;
      int32_t _M0L14capacity__maskS2182;
      int32_t _M0L6_2atmpS2180;
      if (_M0L4hashS2173 == _M0L4hashS792) {
        moonbit_string_t _M0L3keyS2172 = _M0L14_2acurr__entryS799->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3437
        = _M0L3keyS2172 == _M0L3keyS796
          || Moonbit_array_length(_M0L3keyS2172)
             == Moonbit_array_length(_M0L3keyS796)
             && 0
                == memcmp(_M0L3keyS2172, _M0L3keyS796, Moonbit_array_length(_M0L3keyS2172) * 2);
      } else {
        _if__result_3437 = 0;
      }
      if (_if__result_3437) {
        moonbit_string_t _M0L6_2aoldS3050 = _M0L14_2acurr__entryS799->$5;
        moonbit_incref(_M0L5valueS797);
        moonbit_decref(_M0L6_2aoldS3050);
        _M0L14_2acurr__entryS799->$5 = _M0L5valueS797;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS799);
      }
      _M0L3pslS2174 = _M0L14_2acurr__entryS799->$2;
      if (_M0L3pslS787 > _M0L3pslS2174) {
        int32_t _M0L4sizeS2175 = _M0L4selfS790->$1;
        int32_t _M0L8grow__atS2176 = _M0L4selfS790->$4;
        int32_t _M0L7_2abindS800;
        struct _M0TPB5EntryGssE* _M0L7_2abindS801;
        struct _M0TPB5EntryGssE* _M0L5entryS802;
        if (_M0L4sizeS2175 >= _M0L8grow__atS2176) {
          int32_t _M0L14capacity__maskS2178;
          int32_t _M0L6_2atmpS2177;
          moonbit_decref(_M0L14_2acurr__entryS799);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGssE(_M0L4selfS790);
          _M0L14capacity__maskS2178 = _M0L4selfS790->$3;
          _M0L6_2atmpS2177 = _M0L4hashS792 & _M0L14capacity__maskS2178;
          _M0L3pslS787 = 0;
          _M0L3idxS788 = _M0L6_2atmpS2177;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGssE(_M0L4selfS790, _M0L3idxS788, _M0L14_2acurr__entryS799);
        moonbit_decref(_M0L14_2acurr__entryS799);
        _M0L7_2abindS800 = _M0L4selfS790->$6;
        _M0L7_2abindS801 = 0;
        moonbit_incref(_M0L3keyS796);
        moonbit_incref(_M0L5valueS797);
        _M0L5entryS802
        = (struct _M0TPB5EntryGssE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGssE));
        Moonbit_object_header(_M0L5entryS802)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 64, 0);
        _M0L5entryS802->$0 = _M0L7_2abindS800;
        _M0L5entryS802->$1 = _M0L7_2abindS801;
        _M0L5entryS802->$2 = _M0L3pslS787;
        _M0L5entryS802->$3 = _M0L4hashS792;
        _M0L5entryS802->$4 = _M0L3keyS796;
        _M0L5entryS802->$5 = _M0L5valueS797;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGssE(_M0L4selfS790, _M0L3idxS788, _M0L5entryS802);
        moonbit_decref(_M0L5entryS802);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS799);
      }
      _M0L6_2atmpS2179 = _M0L3pslS787 + 1;
      _M0L6_2atmpS2181 = _M0L3idxS788 + 1;
      _M0L14capacity__maskS2182 = _M0L4selfS790->$3;
      _M0L6_2atmpS2180 = _M0L6_2atmpS2181 & _M0L14capacity__maskS2182;
      _M0L3pslS787 = _M0L6_2atmpS2179;
      _M0L3idxS788 = _M0L6_2atmpS2180;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS806,
  moonbit_string_t _M0L3keyS812,
  struct _M0TPB5ArrayGsE* _M0L5valueS813,
  int32_t _M0L4hashS808
) {
  int32_t _M0L14capacity__maskS2203;
  int32_t _M0L6_2atmpS2202;
  int32_t _M0L3pslS803;
  int32_t _M0L3idxS804;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2203 = _M0L4selfS806->$3;
  _M0L6_2atmpS2202 = _M0L4hashS808 & _M0L14capacity__maskS2203;
  _M0L3pslS803 = 0;
  _M0L3idxS804 = _M0L6_2atmpS2202;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7entriesS2201 =
      _M0L4selfS806->$0;
    struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS805;
    if (
      _M0L3idxS804 < 0
      || _M0L3idxS804 >= Moonbit_array_length(_M0L7entriesS2201)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS805
    = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)_M0L7entriesS2201[_M0L3idxS804];
    if (_M0L7_2abindS805 == 0) {
      int32_t _M0L4sizeS2186 = _M0L4selfS806->$1;
      int32_t _M0L8grow__atS2187 = _M0L4selfS806->$4;
      int32_t _M0L7_2abindS809;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS810;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5entryS811;
      if (_M0L4sizeS2186 >= _M0L8grow__atS2187) {
        int32_t _M0L14capacity__maskS2189;
        int32_t _M0L6_2atmpS2188;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB5ArrayGsEE(_M0L4selfS806);
        _M0L14capacity__maskS2189 = _M0L4selfS806->$3;
        _M0L6_2atmpS2188 = _M0L4hashS808 & _M0L14capacity__maskS2189;
        _M0L3pslS803 = 0;
        _M0L3idxS804 = _M0L6_2atmpS2188;
        continue;
      }
      _M0L7_2abindS809 = _M0L4selfS806->$6;
      _M0L7_2abindS810 = 0;
      moonbit_incref(_M0L3keyS812);
      moonbit_incref(_M0L5valueS813);
      _M0L5entryS811
      = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGsEE));
      Moonbit_object_header(_M0L5entryS811)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 69, 0);
      _M0L5entryS811->$0 = _M0L7_2abindS809;
      _M0L5entryS811->$1 = _M0L7_2abindS810;
      _M0L5entryS811->$2 = _M0L3pslS803;
      _M0L5entryS811->$3 = _M0L4hashS808;
      _M0L5entryS811->$4 = _M0L3keyS812;
      _M0L5entryS811->$5 = _M0L5valueS813;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGsEE(_M0L4selfS806, _M0L3idxS804, _M0L5entryS811);
      moonbit_decref(_M0L5entryS811);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2aSomeS814 = _M0L7_2abindS805;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L14_2acurr__entryS815 =
        _M0L7_2aSomeS814;
      int32_t _M0L4hashS2191 = _M0L14_2acurr__entryS815->$3;
      int32_t _if__result_3439;
      int32_t _M0L3pslS2192;
      int32_t _M0L6_2atmpS2197;
      int32_t _M0L6_2atmpS2199;
      int32_t _M0L14capacity__maskS2200;
      int32_t _M0L6_2atmpS2198;
      if (_M0L4hashS2191 == _M0L4hashS808) {
        moonbit_string_t _M0L3keyS2190 = _M0L14_2acurr__entryS815->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3439
        = _M0L3keyS2190 == _M0L3keyS812
          || Moonbit_array_length(_M0L3keyS2190)
             == Moonbit_array_length(_M0L3keyS812)
             && 0
                == memcmp(_M0L3keyS2190, _M0L3keyS812, Moonbit_array_length(_M0L3keyS2190) * 2);
      } else {
        _if__result_3439 = 0;
      }
      if (_if__result_3439) {
        struct _M0TPB5ArrayGsE* _M0L6_2aoldS3054 =
          _M0L14_2acurr__entryS815->$5;
        moonbit_incref(_M0L5valueS813);
        moonbit_decref(_M0L6_2aoldS3054);
        _M0L14_2acurr__entryS815->$5 = _M0L5valueS813;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS815);
      }
      _M0L3pslS2192 = _M0L14_2acurr__entryS815->$2;
      if (_M0L3pslS803 > _M0L3pslS2192) {
        int32_t _M0L4sizeS2193 = _M0L4selfS806->$1;
        int32_t _M0L8grow__atS2194 = _M0L4selfS806->$4;
        int32_t _M0L7_2abindS816;
        struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS817;
        struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5entryS818;
        if (_M0L4sizeS2193 >= _M0L8grow__atS2194) {
          int32_t _M0L14capacity__maskS2196;
          int32_t _M0L6_2atmpS2195;
          moonbit_decref(_M0L14_2acurr__entryS815);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB5ArrayGsEE(_M0L4selfS806);
          _M0L14capacity__maskS2196 = _M0L4selfS806->$3;
          _M0L6_2atmpS2195 = _M0L4hashS808 & _M0L14capacity__maskS2196;
          _M0L3pslS803 = 0;
          _M0L3idxS804 = _M0L6_2atmpS2195;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGsEE(_M0L4selfS806, _M0L3idxS804, _M0L14_2acurr__entryS815);
        moonbit_decref(_M0L14_2acurr__entryS815);
        _M0L7_2abindS816 = _M0L4selfS806->$6;
        _M0L7_2abindS817 = 0;
        moonbit_incref(_M0L3keyS812);
        moonbit_incref(_M0L5valueS813);
        _M0L5entryS818
        = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB5ArrayGsEE));
        Moonbit_object_header(_M0L5entryS818)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 69, 0);
        _M0L5entryS818->$0 = _M0L7_2abindS816;
        _M0L5entryS818->$1 = _M0L7_2abindS817;
        _M0L5entryS818->$2 = _M0L3pslS803;
        _M0L5entryS818->$3 = _M0L4hashS808;
        _M0L5entryS818->$4 = _M0L3keyS812;
        _M0L5entryS818->$5 = _M0L5valueS813;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGsEE(_M0L4selfS806, _M0L3idxS804, _M0L5entryS818);
        moonbit_decref(_M0L5entryS818);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS815);
      }
      _M0L6_2atmpS2197 = _M0L3pslS803 + 1;
      _M0L6_2atmpS2199 = _M0L3idxS804 + 1;
      _M0L14capacity__maskS2200 = _M0L4selfS806->$3;
      _M0L6_2atmpS2198 = _M0L6_2atmpS2199 & _M0L14capacity__maskS2200;
      _M0L3pslS803 = _M0L6_2atmpS2197;
      _M0L3idxS804 = _M0L6_2atmpS2198;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map15set__with__hashGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS822,
  moonbit_string_t _M0L3keyS828,
  struct _M0TPB13StringBuilder* _M0L5valueS829,
  int32_t _M0L4hashS824
) {
  int32_t _M0L14capacity__maskS2221;
  int32_t _M0L6_2atmpS2220;
  int32_t _M0L3pslS819;
  int32_t _M0L3idxS820;
  #line 133 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L14capacity__maskS2221 = _M0L4selfS822->$3;
  _M0L6_2atmpS2220 = _M0L4hashS824 & _M0L14capacity__maskS2221;
  _M0L3pslS819 = 0;
  _M0L3idxS820 = _M0L6_2atmpS2220;
  while (1) {
    struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2219 =
      _M0L4selfS822->$0;
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS821;
    if (
      _M0L3idxS820 < 0
      || _M0L3idxS820 >= Moonbit_array_length(_M0L7entriesS2219)
    ) {
      #line 141 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS821
    = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2219[
        _M0L3idxS820
      ];
    if (_M0L7_2abindS821 == 0) {
      int32_t _M0L4sizeS2204 = _M0L4selfS822->$1;
      int32_t _M0L8grow__atS2205 = _M0L4selfS822->$4;
      int32_t _M0L7_2abindS825;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS826;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5entryS827;
      if (_M0L4sizeS2204 >= _M0L8grow__atS2205) {
        int32_t _M0L14capacity__maskS2207;
        int32_t _M0L6_2atmpS2206;
        #line 145 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map4growGsRPB13StringBuilderE(_M0L4selfS822);
        _M0L14capacity__maskS2207 = _M0L4selfS822->$3;
        _M0L6_2atmpS2206 = _M0L4hashS824 & _M0L14capacity__maskS2207;
        _M0L3pslS819 = 0;
        _M0L3idxS820 = _M0L6_2atmpS2206;
        continue;
      }
      _M0L7_2abindS825 = _M0L4selfS822->$6;
      _M0L7_2abindS826 = 0;
      moonbit_incref(_M0L3keyS828);
      moonbit_incref(_M0L5valueS829);
      _M0L5entryS827
      = (struct _M0TPB5EntryGsRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB13StringBuilderE));
      Moonbit_object_header(_M0L5entryS827)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 74, 0);
      _M0L5entryS827->$0 = _M0L7_2abindS825;
      _M0L5entryS827->$1 = _M0L7_2abindS826;
      _M0L5entryS827->$2 = _M0L3pslS819;
      _M0L5entryS827->$3 = _M0L4hashS824;
      _M0L5entryS827->$4 = _M0L3keyS828;
      _M0L5entryS827->$5 = _M0L5valueS829;
      #line 150 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB13StringBuilderE(_M0L4selfS822, _M0L3idxS820, _M0L5entryS827);
      moonbit_decref(_M0L5entryS827);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS830 =
        _M0L7_2abindS821;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L14_2acurr__entryS831 =
        _M0L7_2aSomeS830;
      int32_t _M0L4hashS2209 = _M0L14_2acurr__entryS831->$3;
      int32_t _if__result_3441;
      int32_t _M0L3pslS2210;
      int32_t _M0L6_2atmpS2215;
      int32_t _M0L6_2atmpS2217;
      int32_t _M0L14capacity__maskS2218;
      int32_t _M0L6_2atmpS2216;
      if (_M0L4hashS2209 == _M0L4hashS824) {
        moonbit_string_t _M0L3keyS2208 = _M0L14_2acurr__entryS831->$4;
        #line 154 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _if__result_3441
        = _M0L3keyS2208 == _M0L3keyS828
          || Moonbit_array_length(_M0L3keyS2208)
             == Moonbit_array_length(_M0L3keyS828)
             && 0
                == memcmp(_M0L3keyS2208, _M0L3keyS828, Moonbit_array_length(_M0L3keyS2208) * 2);
      } else {
        _if__result_3441 = 0;
      }
      if (_if__result_3441) {
        struct _M0TPB13StringBuilder* _M0L6_2aoldS3058 =
          _M0L14_2acurr__entryS831->$5;
        moonbit_incref(_M0L5valueS829);
        moonbit_decref(_M0L6_2aoldS3058);
        _M0L14_2acurr__entryS831->$5 = _M0L5valueS829;
        return 0;
      } else {
        moonbit_incref(_M0L14_2acurr__entryS831);
      }
      _M0L3pslS2210 = _M0L14_2acurr__entryS831->$2;
      if (_M0L3pslS819 > _M0L3pslS2210) {
        int32_t _M0L4sizeS2211 = _M0L4selfS822->$1;
        int32_t _M0L8grow__atS2212 = _M0L4selfS822->$4;
        int32_t _M0L7_2abindS832;
        struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS833;
        struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5entryS834;
        if (_M0L4sizeS2211 >= _M0L8grow__atS2212) {
          int32_t _M0L14capacity__maskS2214;
          int32_t _M0L6_2atmpS2213;
          moonbit_decref(_M0L14_2acurr__entryS831);
          #line 162 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
          _M0MPB3Map4growGsRPB13StringBuilderE(_M0L4selfS822);
          _M0L14capacity__maskS2214 = _M0L4selfS822->$3;
          _M0L6_2atmpS2213 = _M0L4hashS824 & _M0L14capacity__maskS2214;
          _M0L3pslS819 = 0;
          _M0L3idxS820 = _M0L6_2atmpS2213;
          continue;
        }
        #line 166 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB13StringBuilderE(_M0L4selfS822, _M0L3idxS820, _M0L14_2acurr__entryS831);
        moonbit_decref(_M0L14_2acurr__entryS831);
        _M0L7_2abindS832 = _M0L4selfS822->$6;
        _M0L7_2abindS833 = 0;
        moonbit_incref(_M0L3keyS828);
        moonbit_incref(_M0L5valueS829);
        _M0L5entryS834
        = (struct _M0TPB5EntryGsRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TPB5EntryGsRPB13StringBuilderE));
        Moonbit_object_header(_M0L5entryS834)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 74, 0);
        _M0L5entryS834->$0 = _M0L7_2abindS832;
        _M0L5entryS834->$1 = _M0L7_2abindS833;
        _M0L5entryS834->$2 = _M0L3pslS819;
        _M0L5entryS834->$3 = _M0L4hashS824;
        _M0L5entryS834->$4 = _M0L3keyS828;
        _M0L5entryS834->$5 = _M0L5valueS829;
        #line 168 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB13StringBuilderE(_M0L4selfS822, _M0L3idxS820, _M0L5entryS834);
        moonbit_decref(_M0L5entryS834);
        return 0;
      } else {
        moonbit_decref(_M0L14_2acurr__entryS831);
      }
      _M0L6_2atmpS2215 = _M0L3pslS819 + 1;
      _M0L6_2atmpS2217 = _M0L3idxS820 + 1;
      _M0L14capacity__maskS2218 = _M0L4selfS822->$3;
      _M0L6_2atmpS2216 = _M0L6_2atmpS2217 & _M0L14capacity__maskS2218;
      _M0L3pslS819 = _M0L6_2atmpS2215;
      _M0L3idxS820 = _M0L6_2atmpS2216;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS740
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L9old__headS739;
  int32_t _M0L8capacityS2125;
  int32_t _M0L13new__capacityS741;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2119;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L6_2atmpS2118;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L6_2aoldS3065;
  int32_t _M0L6_2atmpS2120;
  int32_t _M0L8capacityS2122;
  int32_t _M0L6_2atmpS2121;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2123;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3064;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L1xS742;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS739 = _M0L4selfS740->$5;
  _M0L8capacityS2125 = _M0L4selfS740->$2;
  _M0L13new__capacityS741 = _M0L8capacityS2125 << 1;
  _M0L6_2atmpS2119 = 0;
  _M0L6_2atmpS2118
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE**)moonbit_make_ref_array(_M0L13new__capacityS741, _M0L6_2atmpS2119);
  _M0L6_2aoldS3065 = _M0L4selfS740->$0;
  if (_M0L9old__headS739) {
    moonbit_incref(_M0L9old__headS739);
  }
  moonbit_decref(_M0L6_2aoldS3065);
  _M0L4selfS740->$0 = _M0L6_2atmpS2118;
  _M0L4selfS740->$2 = _M0L13new__capacityS741;
  _M0L6_2atmpS2120 = _M0L13new__capacityS741 - 1;
  _M0L4selfS740->$3 = _M0L6_2atmpS2120;
  _M0L8capacityS2122 = _M0L4selfS740->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2121 = _M0FPB21calc__grow__threshold(_M0L8capacityS2122);
  _M0L4selfS740->$4 = _M0L6_2atmpS2121;
  _M0L4selfS740->$1 = 0;
  _M0L6_2atmpS2123 = 0;
  _M0L6_2aoldS3064 = _M0L4selfS740->$5;
  if (_M0L6_2aoldS3064) {
    moonbit_decref(_M0L6_2aoldS3064);
  }
  _M0L4selfS740->$5 = _M0L6_2atmpS2123;
  _M0L4selfS740->$6 = -1;
  _M0L1xS742 = _M0L9old__headS739;
  while (1) {
    if (_M0L1xS742 == 0) {
      if (_M0L1xS742) {
        moonbit_decref(_M0L1xS742);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS744 = _M0L1xS742;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4_2aeS745 = _M0L7_2aSomeS744;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L15next__in__chainS746 =
        _M0L4_2aeS745->$1;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2124 = 0;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3062 =
        _M0L4_2aeS745->$1;
      if (_M0L15next__in__chainS746) {
        moonbit_incref(_M0L15next__in__chainS746);
      }
      if (_M0L6_2aoldS3062) {
        moonbit_decref(_M0L6_2aoldS3062);
      }
      _M0L4_2aeS745->$1 = _M0L6_2atmpS2124;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(_M0L4selfS740, _M0L4_2aeS745);
      moonbit_decref(_M0L4_2aeS745);
      _M0L1xS742 = _M0L15next__in__chainS746;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGssE(struct _M0TPB3MapGssE* _M0L4selfS748) {
  struct _M0TPB5EntryGssE* _M0L9old__headS747;
  int32_t _M0L8capacityS2133;
  int32_t _M0L13new__capacityS749;
  struct _M0TPB5EntryGssE* _M0L6_2atmpS2127;
  struct _M0TPB5EntryGssE** _M0L6_2atmpS2126;
  struct _M0TPB5EntryGssE** _M0L6_2aoldS3070;
  int32_t _M0L6_2atmpS2128;
  int32_t _M0L8capacityS2130;
  int32_t _M0L6_2atmpS2129;
  struct _M0TPB5EntryGssE* _M0L6_2atmpS2131;
  struct _M0TPB5EntryGssE* _M0L6_2aoldS3069;
  struct _M0TPB5EntryGssE* _M0L1xS750;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS747 = _M0L4selfS748->$5;
  _M0L8capacityS2133 = _M0L4selfS748->$2;
  _M0L13new__capacityS749 = _M0L8capacityS2133 << 1;
  _M0L6_2atmpS2127 = 0;
  _M0L6_2atmpS2126
  = (struct _M0TPB5EntryGssE**)moonbit_make_ref_array(_M0L13new__capacityS749, _M0L6_2atmpS2127);
  _M0L6_2aoldS3070 = _M0L4selfS748->$0;
  if (_M0L9old__headS747) {
    moonbit_incref(_M0L9old__headS747);
  }
  moonbit_decref(_M0L6_2aoldS3070);
  _M0L4selfS748->$0 = _M0L6_2atmpS2126;
  _M0L4selfS748->$2 = _M0L13new__capacityS749;
  _M0L6_2atmpS2128 = _M0L13new__capacityS749 - 1;
  _M0L4selfS748->$3 = _M0L6_2atmpS2128;
  _M0L8capacityS2130 = _M0L4selfS748->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2129 = _M0FPB21calc__grow__threshold(_M0L8capacityS2130);
  _M0L4selfS748->$4 = _M0L6_2atmpS2129;
  _M0L4selfS748->$1 = 0;
  _M0L6_2atmpS2131 = 0;
  _M0L6_2aoldS3069 = _M0L4selfS748->$5;
  if (_M0L6_2aoldS3069) {
    moonbit_decref(_M0L6_2aoldS3069);
  }
  _M0L4selfS748->$5 = _M0L6_2atmpS2131;
  _M0L4selfS748->$6 = -1;
  _M0L1xS750 = _M0L9old__headS747;
  while (1) {
    if (_M0L1xS750 == 0) {
      if (_M0L1xS750) {
        moonbit_decref(_M0L1xS750);
      }
      break;
    } else {
      struct _M0TPB5EntryGssE* _M0L7_2aSomeS752 = _M0L1xS750;
      struct _M0TPB5EntryGssE* _M0L4_2aeS753 = _M0L7_2aSomeS752;
      struct _M0TPB5EntryGssE* _M0L15next__in__chainS754 = _M0L4_2aeS753->$1;
      struct _M0TPB5EntryGssE* _M0L6_2atmpS2132 = 0;
      struct _M0TPB5EntryGssE* _M0L6_2aoldS3067 = _M0L4_2aeS753->$1;
      if (_M0L15next__in__chainS754) {
        moonbit_incref(_M0L15next__in__chainS754);
      }
      if (_M0L6_2aoldS3067) {
        moonbit_decref(_M0L6_2aoldS3067);
      }
      _M0L4_2aeS753->$1 = _M0L6_2atmpS2132;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGssE(_M0L4selfS748, _M0L4_2aeS753);
      moonbit_decref(_M0L4_2aeS753);
      _M0L1xS750 = _M0L15next__in__chainS754;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS756
) {
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L9old__headS755;
  int32_t _M0L8capacityS2141;
  int32_t _M0L13new__capacityS757;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS2135;
  struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L6_2atmpS2134;
  struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L6_2aoldS3075;
  int32_t _M0L6_2atmpS2136;
  int32_t _M0L8capacityS2138;
  int32_t _M0L6_2atmpS2137;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS2139;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2aoldS3074;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L1xS758;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS755 = _M0L4selfS756->$5;
  _M0L8capacityS2141 = _M0L4selfS756->$2;
  _M0L13new__capacityS757 = _M0L8capacityS2141 << 1;
  _M0L6_2atmpS2135 = 0;
  _M0L6_2atmpS2134
  = (struct _M0TPB5EntryGsRPB5ArrayGsEE**)moonbit_make_ref_array(_M0L13new__capacityS757, _M0L6_2atmpS2135);
  _M0L6_2aoldS3075 = _M0L4selfS756->$0;
  if (_M0L9old__headS755) {
    moonbit_incref(_M0L9old__headS755);
  }
  moonbit_decref(_M0L6_2aoldS3075);
  _M0L4selfS756->$0 = _M0L6_2atmpS2134;
  _M0L4selfS756->$2 = _M0L13new__capacityS757;
  _M0L6_2atmpS2136 = _M0L13new__capacityS757 - 1;
  _M0L4selfS756->$3 = _M0L6_2atmpS2136;
  _M0L8capacityS2138 = _M0L4selfS756->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2137 = _M0FPB21calc__grow__threshold(_M0L8capacityS2138);
  _M0L4selfS756->$4 = _M0L6_2atmpS2137;
  _M0L4selfS756->$1 = 0;
  _M0L6_2atmpS2139 = 0;
  _M0L6_2aoldS3074 = _M0L4selfS756->$5;
  if (_M0L6_2aoldS3074) {
    moonbit_decref(_M0L6_2aoldS3074);
  }
  _M0L4selfS756->$5 = _M0L6_2atmpS2139;
  _M0L4selfS756->$6 = -1;
  _M0L1xS758 = _M0L9old__headS755;
  while (1) {
    if (_M0L1xS758 == 0) {
      if (_M0L1xS758) {
        moonbit_decref(_M0L1xS758);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2aSomeS760 = _M0L1xS758;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L4_2aeS761 = _M0L7_2aSomeS760;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L15next__in__chainS762 =
        _M0L4_2aeS761->$1;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS2140 = 0;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2aoldS3072 =
        _M0L4_2aeS761->$1;
      if (_M0L15next__in__chainS762) {
        moonbit_incref(_M0L15next__in__chainS762);
      }
      if (_M0L6_2aoldS3072) {
        moonbit_decref(_M0L6_2aoldS3072);
      }
      _M0L4_2aeS761->$1 = _M0L6_2atmpS2140;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB5ArrayGsEE(_M0L4selfS756, _M0L4_2aeS761);
      moonbit_decref(_M0L4_2aeS761);
      _M0L1xS758 = _M0L15next__in__chainS762;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map4growGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS764
) {
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L9old__headS763;
  int32_t _M0L8capacityS2149;
  int32_t _M0L13new__capacityS765;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS2143;
  struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L6_2atmpS2142;
  struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L6_2aoldS3080;
  int32_t _M0L6_2atmpS2144;
  int32_t _M0L8capacityS2146;
  int32_t _M0L6_2atmpS2145;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS2147;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3079;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L1xS766;
  #line 561 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L9old__headS763 = _M0L4selfS764->$5;
  _M0L8capacityS2149 = _M0L4selfS764->$2;
  _M0L13new__capacityS765 = _M0L8capacityS2149 << 1;
  _M0L6_2atmpS2143 = 0;
  _M0L6_2atmpS2142
  = (struct _M0TPB5EntryGsRPB13StringBuilderE**)moonbit_make_ref_array(_M0L13new__capacityS765, _M0L6_2atmpS2143);
  _M0L6_2aoldS3080 = _M0L4selfS764->$0;
  if (_M0L9old__headS763) {
    moonbit_incref(_M0L9old__headS763);
  }
  moonbit_decref(_M0L6_2aoldS3080);
  _M0L4selfS764->$0 = _M0L6_2atmpS2142;
  _M0L4selfS764->$2 = _M0L13new__capacityS765;
  _M0L6_2atmpS2144 = _M0L13new__capacityS765 - 1;
  _M0L4selfS764->$3 = _M0L6_2atmpS2144;
  _M0L8capacityS2146 = _M0L4selfS764->$2;
  #line 567 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS2145 = _M0FPB21calc__grow__threshold(_M0L8capacityS2146);
  _M0L4selfS764->$4 = _M0L6_2atmpS2145;
  _M0L4selfS764->$1 = 0;
  _M0L6_2atmpS2147 = 0;
  _M0L6_2aoldS3079 = _M0L4selfS764->$5;
  if (_M0L6_2aoldS3079) {
    moonbit_decref(_M0L6_2aoldS3079);
  }
  _M0L4selfS764->$5 = _M0L6_2atmpS2147;
  _M0L4selfS764->$6 = -1;
  _M0L1xS766 = _M0L9old__headS763;
  while (1) {
    if (_M0L1xS766 == 0) {
      if (_M0L1xS766) {
        moonbit_decref(_M0L1xS766);
      }
      break;
    } else {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS768 = _M0L1xS766;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L4_2aeS769 =
        _M0L7_2aSomeS768;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L15next__in__chainS770 =
        _M0L4_2aeS769->$1;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS2148 = 0;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3077 =
        _M0L4_2aeS769->$1;
      if (_M0L15next__in__chainS770) {
        moonbit_incref(_M0L15next__in__chainS770);
      }
      if (_M0L6_2aoldS3077) {
        moonbit_decref(_M0L6_2aoldS3077);
      }
      _M0L4_2aeS769->$1 = _M0L6_2atmpS2148;
      #line 577 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20rehash__place__entryGsRPB13StringBuilderE(_M0L4selfS764, _M0L4_2aeS769);
      moonbit_decref(_M0L4_2aeS769);
      _M0L1xS766 = _M0L15next__in__chainS770;
      continue;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS708,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5outerS704
) {
  int32_t _M0L4hashS703;
  int32_t _M0L14capacity__maskS2087;
  int32_t _M0L6_2atmpS2086;
  int32_t _M0L3pslS705;
  int32_t _M0L3idxS706;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS703 = _M0L5outerS704->$3;
  _M0L14capacity__maskS2087 = _M0L4selfS708->$3;
  _M0L6_2atmpS2086 = _M0L4hashS703 & _M0L14capacity__maskS2087;
  _M0L3pslS705 = 0;
  _M0L3idxS706 = _M0L6_2atmpS2086;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2085 =
      _M0L4selfS708->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS707;
    if (
      _M0L3idxS706 < 0
      || _M0L3idxS706 >= Moonbit_array_length(_M0L7entriesS2085)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS707
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2085[_M0L3idxS706];
    if (_M0L7_2abindS707 == 0) {
      int32_t _M0L4tailS2078;
      _M0L5outerS704->$2 = _M0L3pslS705;
      _M0L4tailS2078 = _M0L4selfS708->$6;
      _M0L5outerS704->$0 = _M0L4tailS2078;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS708, _M0L3idxS706, _M0L5outerS704);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS709 = _M0L7_2abindS707;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2acurrS710 = _M0L7_2aSomeS709;
      int32_t _M0L3pslS2079 = _M0L7_2acurrS710->$2;
      if (_M0L3pslS705 > _M0L3pslS2079) {
        int32_t _M0L4tailS2080;
        moonbit_incref(_M0L7_2acurrS710);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGiEE(_M0L4selfS708, _M0L3idxS706, _M0L7_2acurrS710);
        moonbit_decref(_M0L7_2acurrS710);
        _M0L5outerS704->$2 = _M0L3pslS705;
        _M0L4tailS2080 = _M0L4selfS708->$6;
        _M0L5outerS704->$0 = _M0L4tailS2080;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(_M0L4selfS708, _M0L3idxS706, _M0L5outerS704);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2081 = _M0L3pslS705 + 1;
        int32_t _M0L6_2atmpS2083 = _M0L3idxS706 + 1;
        int32_t _M0L14capacity__maskS2084 = _M0L4selfS708->$3;
        int32_t _M0L6_2atmpS2082 =
          _M0L6_2atmpS2083 & _M0L14capacity__maskS2084;
        _M0L3pslS705 = _M0L6_2atmpS2081;
        _M0L3idxS706 = _M0L6_2atmpS2082;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGssE(
  struct _M0TPB3MapGssE* _M0L4selfS717,
  struct _M0TPB5EntryGssE* _M0L5outerS713
) {
  int32_t _M0L4hashS712;
  int32_t _M0L14capacity__maskS2097;
  int32_t _M0L6_2atmpS2096;
  int32_t _M0L3pslS714;
  int32_t _M0L3idxS715;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS712 = _M0L5outerS713->$3;
  _M0L14capacity__maskS2097 = _M0L4selfS717->$3;
  _M0L6_2atmpS2096 = _M0L4hashS712 & _M0L14capacity__maskS2097;
  _M0L3pslS714 = 0;
  _M0L3idxS715 = _M0L6_2atmpS2096;
  while (1) {
    struct _M0TPB5EntryGssE** _M0L7entriesS2095 = _M0L4selfS717->$0;
    struct _M0TPB5EntryGssE* _M0L7_2abindS716;
    if (
      _M0L3idxS715 < 0
      || _M0L3idxS715 >= Moonbit_array_length(_M0L7entriesS2095)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS716
    = (struct _M0TPB5EntryGssE*)_M0L7entriesS2095[_M0L3idxS715];
    if (_M0L7_2abindS716 == 0) {
      int32_t _M0L4tailS2088;
      _M0L5outerS713->$2 = _M0L3pslS714;
      _M0L4tailS2088 = _M0L4selfS717->$6;
      _M0L5outerS713->$0 = _M0L4tailS2088;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGssE(_M0L4selfS717, _M0L3idxS715, _M0L5outerS713);
      return 0;
    } else {
      struct _M0TPB5EntryGssE* _M0L7_2aSomeS718 = _M0L7_2abindS716;
      struct _M0TPB5EntryGssE* _M0L7_2acurrS719 = _M0L7_2aSomeS718;
      int32_t _M0L3pslS2089 = _M0L7_2acurrS719->$2;
      if (_M0L3pslS714 > _M0L3pslS2089) {
        int32_t _M0L4tailS2090;
        moonbit_incref(_M0L7_2acurrS719);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGssE(_M0L4selfS717, _M0L3idxS715, _M0L7_2acurrS719);
        moonbit_decref(_M0L7_2acurrS719);
        _M0L5outerS713->$2 = _M0L3pslS714;
        _M0L4tailS2090 = _M0L4selfS717->$6;
        _M0L5outerS713->$0 = _M0L4tailS2090;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGssE(_M0L4selfS717, _M0L3idxS715, _M0L5outerS713);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2091 = _M0L3pslS714 + 1;
        int32_t _M0L6_2atmpS2093 = _M0L3idxS715 + 1;
        int32_t _M0L14capacity__maskS2094 = _M0L4selfS717->$3;
        int32_t _M0L6_2atmpS2092 =
          _M0L6_2atmpS2093 & _M0L14capacity__maskS2094;
        _M0L3pslS714 = _M0L6_2atmpS2091;
        _M0L3idxS715 = _M0L6_2atmpS2092;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS726,
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5outerS722
) {
  int32_t _M0L4hashS721;
  int32_t _M0L14capacity__maskS2107;
  int32_t _M0L6_2atmpS2106;
  int32_t _M0L3pslS723;
  int32_t _M0L3idxS724;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS721 = _M0L5outerS722->$3;
  _M0L14capacity__maskS2107 = _M0L4selfS726->$3;
  _M0L6_2atmpS2106 = _M0L4hashS721 & _M0L14capacity__maskS2107;
  _M0L3pslS723 = 0;
  _M0L3idxS724 = _M0L6_2atmpS2106;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7entriesS2105 =
      _M0L4selfS726->$0;
    struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS725;
    if (
      _M0L3idxS724 < 0
      || _M0L3idxS724 >= Moonbit_array_length(_M0L7entriesS2105)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS725
    = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)_M0L7entriesS2105[_M0L3idxS724];
    if (_M0L7_2abindS725 == 0) {
      int32_t _M0L4tailS2098;
      _M0L5outerS722->$2 = _M0L3pslS723;
      _M0L4tailS2098 = _M0L4selfS726->$6;
      _M0L5outerS722->$0 = _M0L4tailS2098;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGsEE(_M0L4selfS726, _M0L3idxS724, _M0L5outerS722);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2aSomeS727 = _M0L7_2abindS725;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2acurrS728 = _M0L7_2aSomeS727;
      int32_t _M0L3pslS2099 = _M0L7_2acurrS728->$2;
      if (_M0L3pslS723 > _M0L3pslS2099) {
        int32_t _M0L4tailS2100;
        moonbit_incref(_M0L7_2acurrS728);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB5ArrayGsEE(_M0L4selfS726, _M0L3idxS724, _M0L7_2acurrS728);
        moonbit_decref(_M0L7_2acurrS728);
        _M0L5outerS722->$2 = _M0L3pslS723;
        _M0L4tailS2100 = _M0L4selfS726->$6;
        _M0L5outerS722->$0 = _M0L4tailS2100;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGsEE(_M0L4selfS726, _M0L3idxS724, _M0L5outerS722);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2101 = _M0L3pslS723 + 1;
        int32_t _M0L6_2atmpS2103 = _M0L3idxS724 + 1;
        int32_t _M0L14capacity__maskS2104 = _M0L4selfS726->$3;
        int32_t _M0L6_2atmpS2102 =
          _M0L6_2atmpS2103 & _M0L14capacity__maskS2104;
        _M0L3pslS723 = _M0L6_2atmpS2101;
        _M0L3idxS724 = _M0L6_2atmpS2102;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map20rehash__place__entryGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS735,
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5outerS731
) {
  int32_t _M0L4hashS730;
  int32_t _M0L14capacity__maskS2117;
  int32_t _M0L6_2atmpS2116;
  int32_t _M0L3pslS732;
  int32_t _M0L3idxS733;
  #line 585 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L4hashS730 = _M0L5outerS731->$3;
  _M0L14capacity__maskS2117 = _M0L4selfS735->$3;
  _M0L6_2atmpS2116 = _M0L4hashS730 & _M0L14capacity__maskS2117;
  _M0L3pslS732 = 0;
  _M0L3idxS733 = _M0L6_2atmpS2116;
  while (1) {
    struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2115 =
      _M0L4selfS735->$0;
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS734;
    if (
      _M0L3idxS733 < 0
      || _M0L3idxS733 >= Moonbit_array_length(_M0L7entriesS2115)
    ) {
      #line 588 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS734
    = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2115[
        _M0L3idxS733
      ];
    if (_M0L7_2abindS734 == 0) {
      int32_t _M0L4tailS2108;
      _M0L5outerS731->$2 = _M0L3pslS732;
      _M0L4tailS2108 = _M0L4selfS735->$6;
      _M0L5outerS731->$0 = _M0L4tailS2108;
      #line 592 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map20add__entry__to__tailGsRPB13StringBuilderE(_M0L4selfS735, _M0L3idxS733, _M0L5outerS731);
      return 0;
    } else {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS736 =
        _M0L7_2abindS734;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2acurrS737 =
        _M0L7_2aSomeS736;
      int32_t _M0L3pslS2109 = _M0L7_2acurrS737->$2;
      if (_M0L3pslS732 > _M0L3pslS2109) {
        int32_t _M0L4tailS2110;
        moonbit_incref(_M0L7_2acurrS737);
        #line 597 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10push__awayGsRPB13StringBuilderE(_M0L4selfS735, _M0L3idxS733, _M0L7_2acurrS737);
        moonbit_decref(_M0L7_2acurrS737);
        _M0L5outerS731->$2 = _M0L3pslS732;
        _M0L4tailS2110 = _M0L4selfS735->$6;
        _M0L5outerS731->$0 = _M0L4tailS2110;
        #line 600 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map20add__entry__to__tailGsRPB13StringBuilderE(_M0L4selfS735, _M0L3idxS733, _M0L5outerS731);
        return 0;
      } else {
        int32_t _M0L6_2atmpS2111 = _M0L3pslS732 + 1;
        int32_t _M0L6_2atmpS2113 = _M0L3idxS733 + 1;
        int32_t _M0L14capacity__maskS2114 = _M0L4selfS735->$3;
        int32_t _M0L6_2atmpS2112 =
          _M0L6_2atmpS2113 & _M0L14capacity__maskS2114;
        _M0L3pslS732 = _M0L6_2atmpS2111;
        _M0L3idxS733 = _M0L6_2atmpS2112;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS667,
  int32_t _M0L3idxS672,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS671
) {
  int32_t _M0L3pslS2029;
  int32_t _M0L6_2atmpS2025;
  int32_t _M0L6_2atmpS2027;
  int32_t _M0L14capacity__maskS2028;
  int32_t _M0L6_2atmpS2026;
  int32_t _M0L3pslS663;
  int32_t _M0L3idxS664;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS665;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2029 = _M0L5entryS671->$2;
  _M0L6_2atmpS2025 = _M0L3pslS2029 + 1;
  _M0L6_2atmpS2027 = _M0L3idxS672 + 1;
  _M0L14capacity__maskS2028 = _M0L4selfS667->$3;
  _M0L6_2atmpS2026 = _M0L6_2atmpS2027 & _M0L14capacity__maskS2028;
  moonbit_incref(_M0L5entryS671);
  _M0L3pslS663 = _M0L6_2atmpS2025;
  _M0L3idxS664 = _M0L6_2atmpS2026;
  _M0L5entryS665 = _M0L5entryS671;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2024 =
      _M0L4selfS667->$0;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS666;
    if (
      _M0L3idxS664 < 0
      || _M0L3idxS664 >= Moonbit_array_length(_M0L7entriesS2024)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS666
    = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2024[_M0L3idxS664];
    if (_M0L7_2abindS666 == 0) {
      _M0L5entryS665->$2 = _M0L3pslS663;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB5ArrayGiEE(_M0L4selfS667, _M0L5entryS665, _M0L3idxS664);
      moonbit_decref(_M0L5entryS665);
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS669 = _M0L7_2abindS666;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L14_2acurr__entryS670 =
        _M0L7_2aSomeS669;
      int32_t _M0L3pslS2014 = _M0L14_2acurr__entryS670->$2;
      if (_M0L3pslS663 > _M0L3pslS2014) {
        int32_t _M0L3pslS2019;
        int32_t _M0L6_2atmpS2015;
        int32_t _M0L6_2atmpS2017;
        int32_t _M0L14capacity__maskS2018;
        int32_t _M0L6_2atmpS2016;
        _M0L5entryS665->$2 = _M0L3pslS663;
        moonbit_incref(_M0L14_2acurr__entryS670);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB5ArrayGiEE(_M0L4selfS667, _M0L5entryS665, _M0L3idxS664);
        moonbit_decref(_M0L5entryS665);
        _M0L3pslS2019 = _M0L14_2acurr__entryS670->$2;
        _M0L6_2atmpS2015 = _M0L3pslS2019 + 1;
        _M0L6_2atmpS2017 = _M0L3idxS664 + 1;
        _M0L14capacity__maskS2018 = _M0L4selfS667->$3;
        _M0L6_2atmpS2016 = _M0L6_2atmpS2017 & _M0L14capacity__maskS2018;
        _M0L3pslS663 = _M0L6_2atmpS2015;
        _M0L3idxS664 = _M0L6_2atmpS2016;
        _M0L5entryS665 = _M0L14_2acurr__entryS670;
        continue;
      } else {
        int32_t _M0L6_2atmpS2020 = _M0L3pslS663 + 1;
        int32_t _M0L6_2atmpS2022 = _M0L3idxS664 + 1;
        int32_t _M0L14capacity__maskS2023 = _M0L4selfS667->$3;
        int32_t _M0L6_2atmpS2021 =
          _M0L6_2atmpS2022 & _M0L14capacity__maskS2023;
        struct _M0TPB5EntryGsRPB5ArrayGiEE* _tmp_3451 = _M0L5entryS665;
        _M0L3pslS663 = _M0L6_2atmpS2020;
        _M0L3idxS664 = _M0L6_2atmpS2021;
        _M0L5entryS665 = _tmp_3451;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGssE(
  struct _M0TPB3MapGssE* _M0L4selfS677,
  int32_t _M0L3idxS682,
  struct _M0TPB5EntryGssE* _M0L5entryS681
) {
  int32_t _M0L3pslS2045;
  int32_t _M0L6_2atmpS2041;
  int32_t _M0L6_2atmpS2043;
  int32_t _M0L14capacity__maskS2044;
  int32_t _M0L6_2atmpS2042;
  int32_t _M0L3pslS673;
  int32_t _M0L3idxS674;
  struct _M0TPB5EntryGssE* _M0L5entryS675;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2045 = _M0L5entryS681->$2;
  _M0L6_2atmpS2041 = _M0L3pslS2045 + 1;
  _M0L6_2atmpS2043 = _M0L3idxS682 + 1;
  _M0L14capacity__maskS2044 = _M0L4selfS677->$3;
  _M0L6_2atmpS2042 = _M0L6_2atmpS2043 & _M0L14capacity__maskS2044;
  moonbit_incref(_M0L5entryS681);
  _M0L3pslS673 = _M0L6_2atmpS2041;
  _M0L3idxS674 = _M0L6_2atmpS2042;
  _M0L5entryS675 = _M0L5entryS681;
  while (1) {
    struct _M0TPB5EntryGssE** _M0L7entriesS2040 = _M0L4selfS677->$0;
    struct _M0TPB5EntryGssE* _M0L7_2abindS676;
    if (
      _M0L3idxS674 < 0
      || _M0L3idxS674 >= Moonbit_array_length(_M0L7entriesS2040)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS676
    = (struct _M0TPB5EntryGssE*)_M0L7entriesS2040[_M0L3idxS674];
    if (_M0L7_2abindS676 == 0) {
      _M0L5entryS675->$2 = _M0L3pslS673;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGssE(_M0L4selfS677, _M0L5entryS675, _M0L3idxS674);
      moonbit_decref(_M0L5entryS675);
      break;
    } else {
      struct _M0TPB5EntryGssE* _M0L7_2aSomeS679 = _M0L7_2abindS676;
      struct _M0TPB5EntryGssE* _M0L14_2acurr__entryS680 = _M0L7_2aSomeS679;
      int32_t _M0L3pslS2030 = _M0L14_2acurr__entryS680->$2;
      if (_M0L3pslS673 > _M0L3pslS2030) {
        int32_t _M0L3pslS2035;
        int32_t _M0L6_2atmpS2031;
        int32_t _M0L6_2atmpS2033;
        int32_t _M0L14capacity__maskS2034;
        int32_t _M0L6_2atmpS2032;
        _M0L5entryS675->$2 = _M0L3pslS673;
        moonbit_incref(_M0L14_2acurr__entryS680);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGssE(_M0L4selfS677, _M0L5entryS675, _M0L3idxS674);
        moonbit_decref(_M0L5entryS675);
        _M0L3pslS2035 = _M0L14_2acurr__entryS680->$2;
        _M0L6_2atmpS2031 = _M0L3pslS2035 + 1;
        _M0L6_2atmpS2033 = _M0L3idxS674 + 1;
        _M0L14capacity__maskS2034 = _M0L4selfS677->$3;
        _M0L6_2atmpS2032 = _M0L6_2atmpS2033 & _M0L14capacity__maskS2034;
        _M0L3pslS673 = _M0L6_2atmpS2031;
        _M0L3idxS674 = _M0L6_2atmpS2032;
        _M0L5entryS675 = _M0L14_2acurr__entryS680;
        continue;
      } else {
        int32_t _M0L6_2atmpS2036 = _M0L3pslS673 + 1;
        int32_t _M0L6_2atmpS2038 = _M0L3idxS674 + 1;
        int32_t _M0L14capacity__maskS2039 = _M0L4selfS677->$3;
        int32_t _M0L6_2atmpS2037 =
          _M0L6_2atmpS2038 & _M0L14capacity__maskS2039;
        struct _M0TPB5EntryGssE* _tmp_3453 = _M0L5entryS675;
        _M0L3pslS673 = _M0L6_2atmpS2036;
        _M0L3idxS674 = _M0L6_2atmpS2037;
        _M0L5entryS675 = _tmp_3453;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS687,
  int32_t _M0L3idxS692,
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5entryS691
) {
  int32_t _M0L3pslS2061;
  int32_t _M0L6_2atmpS2057;
  int32_t _M0L6_2atmpS2059;
  int32_t _M0L14capacity__maskS2060;
  int32_t _M0L6_2atmpS2058;
  int32_t _M0L3pslS683;
  int32_t _M0L3idxS684;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5entryS685;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2061 = _M0L5entryS691->$2;
  _M0L6_2atmpS2057 = _M0L3pslS2061 + 1;
  _M0L6_2atmpS2059 = _M0L3idxS692 + 1;
  _M0L14capacity__maskS2060 = _M0L4selfS687->$3;
  _M0L6_2atmpS2058 = _M0L6_2atmpS2059 & _M0L14capacity__maskS2060;
  moonbit_incref(_M0L5entryS691);
  _M0L3pslS683 = _M0L6_2atmpS2057;
  _M0L3idxS684 = _M0L6_2atmpS2058;
  _M0L5entryS685 = _M0L5entryS691;
  while (1) {
    struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7entriesS2056 =
      _M0L4selfS687->$0;
    struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS686;
    if (
      _M0L3idxS684 < 0
      || _M0L3idxS684 >= Moonbit_array_length(_M0L7entriesS2056)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS686
    = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)_M0L7entriesS2056[_M0L3idxS684];
    if (_M0L7_2abindS686 == 0) {
      _M0L5entryS685->$2 = _M0L3pslS683;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB5ArrayGsEE(_M0L4selfS687, _M0L5entryS685, _M0L3idxS684);
      moonbit_decref(_M0L5entryS685);
      break;
    } else {
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2aSomeS689 = _M0L7_2abindS686;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L14_2acurr__entryS690 =
        _M0L7_2aSomeS689;
      int32_t _M0L3pslS2046 = _M0L14_2acurr__entryS690->$2;
      if (_M0L3pslS683 > _M0L3pslS2046) {
        int32_t _M0L3pslS2051;
        int32_t _M0L6_2atmpS2047;
        int32_t _M0L6_2atmpS2049;
        int32_t _M0L14capacity__maskS2050;
        int32_t _M0L6_2atmpS2048;
        _M0L5entryS685->$2 = _M0L3pslS683;
        moonbit_incref(_M0L14_2acurr__entryS690);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB5ArrayGsEE(_M0L4selfS687, _M0L5entryS685, _M0L3idxS684);
        moonbit_decref(_M0L5entryS685);
        _M0L3pslS2051 = _M0L14_2acurr__entryS690->$2;
        _M0L6_2atmpS2047 = _M0L3pslS2051 + 1;
        _M0L6_2atmpS2049 = _M0L3idxS684 + 1;
        _M0L14capacity__maskS2050 = _M0L4selfS687->$3;
        _M0L6_2atmpS2048 = _M0L6_2atmpS2049 & _M0L14capacity__maskS2050;
        _M0L3pslS683 = _M0L6_2atmpS2047;
        _M0L3idxS684 = _M0L6_2atmpS2048;
        _M0L5entryS685 = _M0L14_2acurr__entryS690;
        continue;
      } else {
        int32_t _M0L6_2atmpS2052 = _M0L3pslS683 + 1;
        int32_t _M0L6_2atmpS2054 = _M0L3idxS684 + 1;
        int32_t _M0L14capacity__maskS2055 = _M0L4selfS687->$3;
        int32_t _M0L6_2atmpS2053 =
          _M0L6_2atmpS2054 & _M0L14capacity__maskS2055;
        struct _M0TPB5EntryGsRPB5ArrayGsEE* _tmp_3455 = _M0L5entryS685;
        _M0L3pslS683 = _M0L6_2atmpS2052;
        _M0L3idxS684 = _M0L6_2atmpS2053;
        _M0L5entryS685 = _tmp_3455;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10push__awayGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS697,
  int32_t _M0L3idxS702,
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5entryS701
) {
  int32_t _M0L3pslS2077;
  int32_t _M0L6_2atmpS2073;
  int32_t _M0L6_2atmpS2075;
  int32_t _M0L14capacity__maskS2076;
  int32_t _M0L6_2atmpS2074;
  int32_t _M0L3pslS693;
  int32_t _M0L3idxS694;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5entryS695;
  #line 178 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L3pslS2077 = _M0L5entryS701->$2;
  _M0L6_2atmpS2073 = _M0L3pslS2077 + 1;
  _M0L6_2atmpS2075 = _M0L3idxS702 + 1;
  _M0L14capacity__maskS2076 = _M0L4selfS697->$3;
  _M0L6_2atmpS2074 = _M0L6_2atmpS2075 & _M0L14capacity__maskS2076;
  moonbit_incref(_M0L5entryS701);
  _M0L3pslS693 = _M0L6_2atmpS2073;
  _M0L3idxS694 = _M0L6_2atmpS2074;
  _M0L5entryS695 = _M0L5entryS701;
  while (1) {
    struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2072 =
      _M0L4selfS697->$0;
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS696;
    if (
      _M0L3idxS694 < 0
      || _M0L3idxS694 >= Moonbit_array_length(_M0L7entriesS2072)
    ) {
      #line 184 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      moonbit_panic();
    }
    _M0L7_2abindS696
    = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2072[
        _M0L3idxS694
      ];
    if (_M0L7_2abindS696 == 0) {
      _M0L5entryS695->$2 = _M0L3pslS693;
      #line 187 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0MPB3Map10set__entryGsRPB13StringBuilderE(_M0L4selfS697, _M0L5entryS695, _M0L3idxS694);
      moonbit_decref(_M0L5entryS695);
      break;
    } else {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS699 =
        _M0L7_2abindS696;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L14_2acurr__entryS700 =
        _M0L7_2aSomeS699;
      int32_t _M0L3pslS2062 = _M0L14_2acurr__entryS700->$2;
      if (_M0L3pslS693 > _M0L3pslS2062) {
        int32_t _M0L3pslS2067;
        int32_t _M0L6_2atmpS2063;
        int32_t _M0L6_2atmpS2065;
        int32_t _M0L14capacity__maskS2066;
        int32_t _M0L6_2atmpS2064;
        _M0L5entryS695->$2 = _M0L3pslS693;
        moonbit_incref(_M0L14_2acurr__entryS700);
        #line 193 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        _M0MPB3Map10set__entryGsRPB13StringBuilderE(_M0L4selfS697, _M0L5entryS695, _M0L3idxS694);
        moonbit_decref(_M0L5entryS695);
        _M0L3pslS2067 = _M0L14_2acurr__entryS700->$2;
        _M0L6_2atmpS2063 = _M0L3pslS2067 + 1;
        _M0L6_2atmpS2065 = _M0L3idxS694 + 1;
        _M0L14capacity__maskS2066 = _M0L4selfS697->$3;
        _M0L6_2atmpS2064 = _M0L6_2atmpS2065 & _M0L14capacity__maskS2066;
        _M0L3pslS693 = _M0L6_2atmpS2063;
        _M0L3idxS694 = _M0L6_2atmpS2064;
        _M0L5entryS695 = _M0L14_2acurr__entryS700;
        continue;
      } else {
        int32_t _M0L6_2atmpS2068 = _M0L3pslS693 + 1;
        int32_t _M0L6_2atmpS2070 = _M0L3idxS694 + 1;
        int32_t _M0L14capacity__maskS2071 = _M0L4selfS697->$3;
        int32_t _M0L6_2atmpS2069 =
          _M0L6_2atmpS2070 & _M0L14capacity__maskS2071;
        struct _M0TPB5EntryGsRPB13StringBuilderE* _tmp_3457 = _M0L5entryS695;
        _M0L3pslS693 = _M0L6_2atmpS2068;
        _M0L3idxS694 = _M0L6_2atmpS2069;
        _M0L5entryS695 = _tmp_3457;
        continue;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS639,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS641,
  int32_t _M0L8new__idxS640
) {
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS2006;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS2007;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3099;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS642;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2006 = _M0L4selfS639->$0;
  _M0L6_2atmpS2007 = _M0L5entryS641;
  if (
    _M0L8new__idxS640 < 0
    || _M0L8new__idxS640 >= Moonbit_array_length(_M0L7entriesS2006)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3099
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS2006[_M0L8new__idxS640];
  if (_M0L6_2atmpS2007) {
    moonbit_incref(_M0L6_2atmpS2007);
  }
  if (_M0L6_2aoldS3099) {
    moonbit_decref(_M0L6_2aoldS3099);
  }
  _M0L7entriesS2006[_M0L8new__idxS640] = _M0L6_2atmpS2007;
  _M0L7_2abindS642 = _M0L5entryS641->$1;
  if (_M0L7_2abindS642 == 0) {
    _M0L4selfS639->$6 = _M0L8new__idxS640;
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS643 = _M0L7_2abindS642;
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2anextS644 = _M0L7_2aSomeS643;
    _M0L7_2anextS644->$0 = _M0L8new__idxS640;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGssE(
  struct _M0TPB3MapGssE* _M0L4selfS645,
  struct _M0TPB5EntryGssE* _M0L5entryS647,
  int32_t _M0L8new__idxS646
) {
  struct _M0TPB5EntryGssE** _M0L7entriesS2008;
  struct _M0TPB5EntryGssE* _M0L6_2atmpS2009;
  struct _M0TPB5EntryGssE* _M0L6_2aoldS3102;
  struct _M0TPB5EntryGssE* _M0L7_2abindS648;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2008 = _M0L4selfS645->$0;
  _M0L6_2atmpS2009 = _M0L5entryS647;
  if (
    _M0L8new__idxS646 < 0
    || _M0L8new__idxS646 >= Moonbit_array_length(_M0L7entriesS2008)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3102
  = (struct _M0TPB5EntryGssE*)_M0L7entriesS2008[_M0L8new__idxS646];
  if (_M0L6_2atmpS2009) {
    moonbit_incref(_M0L6_2atmpS2009);
  }
  if (_M0L6_2aoldS3102) {
    moonbit_decref(_M0L6_2aoldS3102);
  }
  _M0L7entriesS2008[_M0L8new__idxS646] = _M0L6_2atmpS2009;
  _M0L7_2abindS648 = _M0L5entryS647->$1;
  if (_M0L7_2abindS648 == 0) {
    _M0L4selfS645->$6 = _M0L8new__idxS646;
  } else {
    struct _M0TPB5EntryGssE* _M0L7_2aSomeS649 = _M0L7_2abindS648;
    struct _M0TPB5EntryGssE* _M0L7_2anextS650 = _M0L7_2aSomeS649;
    _M0L7_2anextS650->$0 = _M0L8new__idxS646;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS651,
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5entryS653,
  int32_t _M0L8new__idxS652
) {
  struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7entriesS2010;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS2011;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2aoldS3105;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS654;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2010 = _M0L4selfS651->$0;
  _M0L6_2atmpS2011 = _M0L5entryS653;
  if (
    _M0L8new__idxS652 < 0
    || _M0L8new__idxS652 >= Moonbit_array_length(_M0L7entriesS2010)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3105
  = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)_M0L7entriesS2010[_M0L8new__idxS652];
  if (_M0L6_2atmpS2011) {
    moonbit_incref(_M0L6_2atmpS2011);
  }
  if (_M0L6_2aoldS3105) {
    moonbit_decref(_M0L6_2aoldS3105);
  }
  _M0L7entriesS2010[_M0L8new__idxS652] = _M0L6_2atmpS2011;
  _M0L7_2abindS654 = _M0L5entryS653->$1;
  if (_M0L7_2abindS654 == 0) {
    _M0L4selfS651->$6 = _M0L8new__idxS652;
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2aSomeS655 = _M0L7_2abindS654;
    struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2anextS656 = _M0L7_2aSomeS655;
    _M0L7_2anextS656->$0 = _M0L8new__idxS652;
  }
  return 0;
}

int32_t _M0MPB3Map10set__entryGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS657,
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5entryS659,
  int32_t _M0L8new__idxS658
) {
  struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2012;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS2013;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3108;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS660;
  #line 205 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7entriesS2012 = _M0L4selfS657->$0;
  _M0L6_2atmpS2013 = _M0L5entryS659;
  if (
    _M0L8new__idxS658 < 0
    || _M0L8new__idxS658 >= Moonbit_array_length(_M0L7entriesS2012)
  ) {
    #line 210 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3108
  = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2012[
      _M0L8new__idxS658
    ];
  if (_M0L6_2atmpS2013) {
    moonbit_incref(_M0L6_2atmpS2013);
  }
  if (_M0L6_2aoldS3108) {
    moonbit_decref(_M0L6_2aoldS3108);
  }
  _M0L7entriesS2012[_M0L8new__idxS658] = _M0L6_2atmpS2013;
  _M0L7_2abindS660 = _M0L5entryS659->$1;
  if (_M0L7_2abindS660 == 0) {
    _M0L4selfS657->$6 = _M0L8new__idxS658;
  } else {
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS661 =
      _M0L7_2abindS660;
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2anextS662 =
      _M0L7_2aSomeS661;
    _M0L7_2anextS662->$0 = _M0L8new__idxS658;
  }
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGiEE(
  struct _M0TPB3MapGsRPB5ArrayGiEE* _M0L4selfS624,
  int32_t _M0L3idxS626,
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L5entryS625
) {
  int32_t _M0L7_2abindS623;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1975;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1976;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3110;
  int32_t _M0L4sizeS1978;
  int32_t _M0L6_2atmpS1977;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS623 = _M0L4selfS624->$6;
  switch (_M0L7_2abindS623) {
    case -1: {
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1970 = _M0L5entryS625;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3112 =
        _M0L4selfS624->$5;
      if (_M0L6_2atmpS1970) {
        moonbit_incref(_M0L6_2atmpS1970);
      }
      if (_M0L6_2aoldS3112) {
        moonbit_decref(_M0L6_2aoldS3112);
      }
      _M0L4selfS624->$5 = _M0L6_2atmpS1970;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7entriesS1974 =
        _M0L4selfS624->$0;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1973;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1971;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1972;
      struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2aoldS3113;
      if (
        _M0L7_2abindS623 < 0
        || _M0L7_2abindS623 >= Moonbit_array_length(_M0L7entriesS1974)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS1973
      = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1974[
          _M0L7_2abindS623
        ];
      if (_M0L6_2atmpS1973) {
        moonbit_incref(_M0L6_2atmpS1973);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS1971
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(_M0L6_2atmpS1973);
      if (_M0L6_2atmpS1973) {
        moonbit_decref(_M0L6_2atmpS1973);
      }
      _M0L6_2atmpS1972 = _M0L5entryS625;
      _M0L6_2aoldS3113 = _M0L6_2atmpS1971->$1;
      if (_M0L6_2atmpS1972) {
        moonbit_incref(_M0L6_2atmpS1972);
      }
      if (_M0L6_2aoldS3113) {
        moonbit_decref(_M0L6_2aoldS3113);
      }
      _M0L6_2atmpS1971->$1 = _M0L6_2atmpS1972;
      moonbit_decref(_M0L6_2atmpS1971);
      break;
    }
  }
  _M0L4selfS624->$6 = _M0L3idxS626;
  _M0L7entriesS1975 = _M0L4selfS624->$0;
  _M0L6_2atmpS1976 = _M0L5entryS625;
  if (
    _M0L3idxS626 < 0
    || _M0L3idxS626 >= Moonbit_array_length(_M0L7entriesS1975)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3110
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE*)_M0L7entriesS1975[_M0L3idxS626];
  if (_M0L6_2atmpS1976) {
    moonbit_incref(_M0L6_2atmpS1976);
  }
  if (_M0L6_2aoldS3110) {
    moonbit_decref(_M0L6_2aoldS3110);
  }
  _M0L7entriesS1975[_M0L3idxS626] = _M0L6_2atmpS1976;
  _M0L4sizeS1978 = _M0L4selfS624->$1;
  _M0L6_2atmpS1977 = _M0L4sizeS1978 + 1;
  _M0L4selfS624->$1 = _M0L6_2atmpS1977;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGssE(
  struct _M0TPB3MapGssE* _M0L4selfS628,
  int32_t _M0L3idxS630,
  struct _M0TPB5EntryGssE* _M0L5entryS629
) {
  int32_t _M0L7_2abindS627;
  struct _M0TPB5EntryGssE** _M0L7entriesS1984;
  struct _M0TPB5EntryGssE* _M0L6_2atmpS1985;
  struct _M0TPB5EntryGssE* _M0L6_2aoldS3116;
  int32_t _M0L4sizeS1987;
  int32_t _M0L6_2atmpS1986;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS627 = _M0L4selfS628->$6;
  switch (_M0L7_2abindS627) {
    case -1: {
      struct _M0TPB5EntryGssE* _M0L6_2atmpS1979 = _M0L5entryS629;
      struct _M0TPB5EntryGssE* _M0L6_2aoldS3118 = _M0L4selfS628->$5;
      if (_M0L6_2atmpS1979) {
        moonbit_incref(_M0L6_2atmpS1979);
      }
      if (_M0L6_2aoldS3118) {
        moonbit_decref(_M0L6_2aoldS3118);
      }
      _M0L4selfS628->$5 = _M0L6_2atmpS1979;
      break;
    }
    default: {
      struct _M0TPB5EntryGssE** _M0L7entriesS1983 = _M0L4selfS628->$0;
      struct _M0TPB5EntryGssE* _M0L6_2atmpS1982;
      struct _M0TPB5EntryGssE* _M0L6_2atmpS1980;
      struct _M0TPB5EntryGssE* _M0L6_2atmpS1981;
      struct _M0TPB5EntryGssE* _M0L6_2aoldS3119;
      if (
        _M0L7_2abindS627 < 0
        || _M0L7_2abindS627 >= Moonbit_array_length(_M0L7entriesS1983)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS1982
      = (struct _M0TPB5EntryGssE*)_M0L7entriesS1983[_M0L7_2abindS627];
      if (_M0L6_2atmpS1982) {
        moonbit_incref(_M0L6_2atmpS1982);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS1980
      = _M0MPC16option6Option6unwrapGRPB5EntryGssEE(_M0L6_2atmpS1982);
      if (_M0L6_2atmpS1982) {
        moonbit_decref(_M0L6_2atmpS1982);
      }
      _M0L6_2atmpS1981 = _M0L5entryS629;
      _M0L6_2aoldS3119 = _M0L6_2atmpS1980->$1;
      if (_M0L6_2atmpS1981) {
        moonbit_incref(_M0L6_2atmpS1981);
      }
      if (_M0L6_2aoldS3119) {
        moonbit_decref(_M0L6_2aoldS3119);
      }
      _M0L6_2atmpS1980->$1 = _M0L6_2atmpS1981;
      moonbit_decref(_M0L6_2atmpS1980);
      break;
    }
  }
  _M0L4selfS628->$6 = _M0L3idxS630;
  _M0L7entriesS1984 = _M0L4selfS628->$0;
  _M0L6_2atmpS1985 = _M0L5entryS629;
  if (
    _M0L3idxS630 < 0
    || _M0L3idxS630 >= Moonbit_array_length(_M0L7entriesS1984)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3116
  = (struct _M0TPB5EntryGssE*)_M0L7entriesS1984[_M0L3idxS630];
  if (_M0L6_2atmpS1985) {
    moonbit_incref(_M0L6_2atmpS1985);
  }
  if (_M0L6_2aoldS3116) {
    moonbit_decref(_M0L6_2aoldS3116);
  }
  _M0L7entriesS1984[_M0L3idxS630] = _M0L6_2atmpS1985;
  _M0L4sizeS1987 = _M0L4selfS628->$1;
  _M0L6_2atmpS1986 = _M0L4sizeS1987 + 1;
  _M0L4selfS628->$1 = _M0L6_2atmpS1986;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB5ArrayGsEE(
  struct _M0TPB3MapGsRPB5ArrayGsEE* _M0L4selfS632,
  int32_t _M0L3idxS634,
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L5entryS633
) {
  int32_t _M0L7_2abindS631;
  struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7entriesS1993;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS1994;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2aoldS3122;
  int32_t _M0L4sizeS1996;
  int32_t _M0L6_2atmpS1995;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS631 = _M0L4selfS632->$6;
  switch (_M0L7_2abindS631) {
    case -1: {
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS1988 = _M0L5entryS633;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2aoldS3124 =
        _M0L4selfS632->$5;
      if (_M0L6_2atmpS1988) {
        moonbit_incref(_M0L6_2atmpS1988);
      }
      if (_M0L6_2aoldS3124) {
        moonbit_decref(_M0L6_2aoldS3124);
      }
      _M0L4selfS632->$5 = _M0L6_2atmpS1988;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7entriesS1992 =
        _M0L4selfS632->$0;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS1991;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS1989;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS1990;
      struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2aoldS3125;
      if (
        _M0L7_2abindS631 < 0
        || _M0L7_2abindS631 >= Moonbit_array_length(_M0L7entriesS1992)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS1991
      = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)_M0L7entriesS1992[
          _M0L7_2abindS631
        ];
      if (_M0L6_2atmpS1991) {
        moonbit_incref(_M0L6_2atmpS1991);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS1989
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGsEEE(_M0L6_2atmpS1991);
      if (_M0L6_2atmpS1991) {
        moonbit_decref(_M0L6_2atmpS1991);
      }
      _M0L6_2atmpS1990 = _M0L5entryS633;
      _M0L6_2aoldS3125 = _M0L6_2atmpS1989->$1;
      if (_M0L6_2atmpS1990) {
        moonbit_incref(_M0L6_2atmpS1990);
      }
      if (_M0L6_2aoldS3125) {
        moonbit_decref(_M0L6_2aoldS3125);
      }
      _M0L6_2atmpS1989->$1 = _M0L6_2atmpS1990;
      moonbit_decref(_M0L6_2atmpS1989);
      break;
    }
  }
  _M0L4selfS632->$6 = _M0L3idxS634;
  _M0L7entriesS1993 = _M0L4selfS632->$0;
  _M0L6_2atmpS1994 = _M0L5entryS633;
  if (
    _M0L3idxS634 < 0
    || _M0L3idxS634 >= Moonbit_array_length(_M0L7entriesS1993)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3122
  = (struct _M0TPB5EntryGsRPB5ArrayGsEE*)_M0L7entriesS1993[_M0L3idxS634];
  if (_M0L6_2atmpS1994) {
    moonbit_incref(_M0L6_2atmpS1994);
  }
  if (_M0L6_2aoldS3122) {
    moonbit_decref(_M0L6_2aoldS3122);
  }
  _M0L7entriesS1993[_M0L3idxS634] = _M0L6_2atmpS1994;
  _M0L4sizeS1996 = _M0L4selfS632->$1;
  _M0L6_2atmpS1995 = _M0L4sizeS1996 + 1;
  _M0L4selfS632->$1 = _M0L6_2atmpS1995;
  return 0;
}

int32_t _M0MPB3Map20add__entry__to__tailGsRPB13StringBuilderE(
  struct _M0TPB3MapGsRPB13StringBuilderE* _M0L4selfS636,
  int32_t _M0L3idxS638,
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L5entryS637
) {
  int32_t _M0L7_2abindS635;
  struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2002;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS2003;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3128;
  int32_t _M0L4sizeS2005;
  int32_t _M0L6_2atmpS2004;
  #line 516 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS635 = _M0L4selfS636->$6;
  switch (_M0L7_2abindS635) {
    case -1: {
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS1997 =
        _M0L5entryS637;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3130 =
        _M0L4selfS636->$5;
      if (_M0L6_2atmpS1997) {
        moonbit_incref(_M0L6_2atmpS1997);
      }
      if (_M0L6_2aoldS3130) {
        moonbit_decref(_M0L6_2aoldS3130);
      }
      _M0L4selfS636->$5 = _M0L6_2atmpS1997;
      break;
    }
    default: {
      struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7entriesS2001 =
        _M0L4selfS636->$0;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS2000;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS1998;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS1999;
      struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2aoldS3131;
      if (
        _M0L7_2abindS635 < 0
        || _M0L7_2abindS635 >= Moonbit_array_length(_M0L7entriesS2001)
      ) {
        #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS2000
      = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2001[
          _M0L7_2abindS635
        ];
      if (_M0L6_2atmpS2000) {
        moonbit_incref(_M0L6_2atmpS2000);
      }
      #line 523 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
      _M0L6_2atmpS1998
      = _M0MPC16option6Option6unwrapGRPB5EntryGsRPB13StringBuilderEE(_M0L6_2atmpS2000);
      if (_M0L6_2atmpS2000) {
        moonbit_decref(_M0L6_2atmpS2000);
      }
      _M0L6_2atmpS1999 = _M0L5entryS637;
      _M0L6_2aoldS3131 = _M0L6_2atmpS1998->$1;
      if (_M0L6_2atmpS1999) {
        moonbit_incref(_M0L6_2atmpS1999);
      }
      if (_M0L6_2aoldS3131) {
        moonbit_decref(_M0L6_2aoldS3131);
      }
      _M0L6_2atmpS1998->$1 = _M0L6_2atmpS1999;
      moonbit_decref(_M0L6_2atmpS1998);
      break;
    }
  }
  _M0L4selfS636->$6 = _M0L3idxS638;
  _M0L7entriesS2002 = _M0L4selfS636->$0;
  _M0L6_2atmpS2003 = _M0L5entryS637;
  if (
    _M0L3idxS638 < 0
    || _M0L3idxS638 >= Moonbit_array_length(_M0L7entriesS2002)
  ) {
    #line 526 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
    moonbit_panic();
  }
  _M0L6_2aoldS3128
  = (struct _M0TPB5EntryGsRPB13StringBuilderE*)_M0L7entriesS2002[
      _M0L3idxS638
    ];
  if (_M0L6_2atmpS2003) {
    moonbit_incref(_M0L6_2atmpS2003);
  }
  if (_M0L6_2aoldS3128) {
    moonbit_decref(_M0L6_2aoldS3128);
  }
  _M0L7entriesS2002[_M0L3idxS638] = _M0L6_2atmpS2003;
  _M0L4sizeS2005 = _M0L4selfS636->$1;
  _M0L6_2atmpS2004 = _M0L4sizeS2005 + 1;
  _M0L4selfS636->$1 = _M0L6_2atmpS2004;
  return 0;
}

int32_t _M0MPC13int3Int3max(int32_t _M0L4selfS621, int32_t _M0L5otherS622) {
  #line 75 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS621 > _M0L5otherS622) {
    return _M0L4selfS621;
  } else {
    return _M0L5otherS622;
  }
}

int32_t _M0FPB21capacity__for__length(int32_t _M0L6lengthS620) {
  int32_t _M0Lm8capacityS619;
  int32_t _M0L6_2atmpS1968;
  int32_t _M0L6_2atmpS1967;
  #line 71 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 72 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0Lm8capacityS619 = _M0MPC13int3Int20next__power__of__two(_M0L6lengthS620);
  _M0L6_2atmpS1968 = _M0Lm8capacityS619;
  #line 73 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1967 = _M0FPB21calc__grow__threshold(_M0L6_2atmpS1968);
  if (_M0L6lengthS620 > _M0L6_2atmpS1967) {
    int32_t _M0L6_2atmpS1969 = _M0Lm8capacityS619;
    _M0Lm8capacityS619 = _M0L6_2atmpS1969 * 2;
  }
  return _M0Lm8capacityS619;
}

struct _M0TPB3MapGsRPB5ArrayGiEE* _M0FPB8new__mapGsRPB5ArrayGiEE(
  int32_t _M0L8capacityS596
) {
  int32_t _M0L8capacityS595;
  int32_t _M0L7_2abindS597;
  int32_t _M0L7_2abindS598;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L6_2atmpS1963;
  struct _M0TPB5EntryGsRPB5ArrayGiEE** _M0L7_2abindS599;
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2abindS600;
  struct _M0TPB3MapGsRPB5ArrayGiEE* _block_3458;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS595
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS596);
  _M0L7_2abindS597 = _M0L8capacityS595 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS598 = _M0FPB21calc__grow__threshold(_M0L8capacityS595);
  _M0L6_2atmpS1963 = 0;
  _M0L7_2abindS599
  = (struct _M0TPB5EntryGsRPB5ArrayGiEE**)moonbit_make_ref_array(_M0L8capacityS595, _M0L6_2atmpS1963);
  _M0L7_2abindS600 = 0;
  _block_3458
  = (struct _M0TPB3MapGsRPB5ArrayGiEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB5ArrayGiEE));
  Moonbit_object_header(_block_3458)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 79, 0);
  _block_3458->$0 = _M0L7_2abindS599;
  _block_3458->$1 = 0;
  _block_3458->$2 = _M0L8capacityS595;
  _block_3458->$3 = _M0L7_2abindS597;
  _block_3458->$4 = _M0L7_2abindS598;
  _block_3458->$5 = _M0L7_2abindS600;
  _block_3458->$6 = -1;
  return _block_3458;
}

struct _M0TPB3MapGssE* _M0FPB8new__mapGssE(int32_t _M0L8capacityS602) {
  int32_t _M0L8capacityS601;
  int32_t _M0L7_2abindS603;
  int32_t _M0L7_2abindS604;
  struct _M0TPB5EntryGssE* _M0L6_2atmpS1964;
  struct _M0TPB5EntryGssE** _M0L7_2abindS605;
  struct _M0TPB5EntryGssE* _M0L7_2abindS606;
  struct _M0TPB3MapGssE* _block_3459;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS601
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS602);
  _M0L7_2abindS603 = _M0L8capacityS601 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS604 = _M0FPB21calc__grow__threshold(_M0L8capacityS601);
  _M0L6_2atmpS1964 = 0;
  _M0L7_2abindS605
  = (struct _M0TPB5EntryGssE**)moonbit_make_ref_array(_M0L8capacityS601, _M0L6_2atmpS1964);
  _M0L7_2abindS606 = 0;
  _block_3459
  = (struct _M0TPB3MapGssE*)moonbit_malloc(sizeof(struct _M0TPB3MapGssE));
  Moonbit_object_header(_block_3459)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 83, 0);
  _block_3459->$0 = _M0L7_2abindS605;
  _block_3459->$1 = 0;
  _block_3459->$2 = _M0L8capacityS601;
  _block_3459->$3 = _M0L7_2abindS603;
  _block_3459->$4 = _M0L7_2abindS604;
  _block_3459->$5 = _M0L7_2abindS606;
  _block_3459->$6 = -1;
  return _block_3459;
}

struct _M0TPB3MapGsRPB5ArrayGsEE* _M0FPB8new__mapGsRPB5ArrayGsEE(
  int32_t _M0L8capacityS608
) {
  int32_t _M0L8capacityS607;
  int32_t _M0L7_2abindS609;
  int32_t _M0L7_2abindS610;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L6_2atmpS1965;
  struct _M0TPB5EntryGsRPB5ArrayGsEE** _M0L7_2abindS611;
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2abindS612;
  struct _M0TPB3MapGsRPB5ArrayGsEE* _block_3460;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS607
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS608);
  _M0L7_2abindS609 = _M0L8capacityS607 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS610 = _M0FPB21calc__grow__threshold(_M0L8capacityS607);
  _M0L6_2atmpS1965 = 0;
  _M0L7_2abindS611
  = (struct _M0TPB5EntryGsRPB5ArrayGsEE**)moonbit_make_ref_array(_M0L8capacityS607, _M0L6_2atmpS1965);
  _M0L7_2abindS612 = 0;
  _block_3460
  = (struct _M0TPB3MapGsRPB5ArrayGsEE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB5ArrayGsEE));
  Moonbit_object_header(_block_3460)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 87, 0);
  _block_3460->$0 = _M0L7_2abindS611;
  _block_3460->$1 = 0;
  _block_3460->$2 = _M0L8capacityS607;
  _block_3460->$3 = _M0L7_2abindS609;
  _block_3460->$4 = _M0L7_2abindS610;
  _block_3460->$5 = _M0L7_2abindS612;
  _block_3460->$6 = -1;
  return _block_3460;
}

struct _M0TPB3MapGsRPB13StringBuilderE* _M0FPB8new__mapGsRPB13StringBuilderE(
  int32_t _M0L8capacityS614
) {
  int32_t _M0L8capacityS613;
  int32_t _M0L7_2abindS615;
  int32_t _M0L7_2abindS616;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L6_2atmpS1966;
  struct _M0TPB5EntryGsRPB13StringBuilderE** _M0L7_2abindS617;
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2abindS618;
  struct _M0TPB3MapGsRPB13StringBuilderE* _block_3461;
  #line 57 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  #line 58 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L8capacityS613
  = _M0MPC13int3Int20next__power__of__two(_M0L8capacityS614);
  _M0L7_2abindS615 = _M0L8capacityS613 - 1;
  #line 63 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L7_2abindS616 = _M0FPB21calc__grow__threshold(_M0L8capacityS613);
  _M0L6_2atmpS1966 = 0;
  _M0L7_2abindS617
  = (struct _M0TPB5EntryGsRPB13StringBuilderE**)moonbit_make_ref_array(_M0L8capacityS613, _M0L6_2atmpS1966);
  _M0L7_2abindS618 = 0;
  _block_3461
  = (struct _M0TPB3MapGsRPB13StringBuilderE*)moonbit_malloc(sizeof(struct _M0TPB3MapGsRPB13StringBuilderE));
  Moonbit_object_header(_block_3461)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 91, 0);
  _block_3461->$0 = _M0L7_2abindS617;
  _block_3461->$1 = 0;
  _block_3461->$2 = _M0L8capacityS613;
  _block_3461->$3 = _M0L7_2abindS615;
  _block_3461->$4 = _M0L7_2abindS616;
  _block_3461->$5 = _M0L7_2abindS618;
  _block_3461->$6 = -1;
  return _block_3461;
}

int32_t _M0MPC13int3Int20next__power__of__two(int32_t _M0L4selfS594) {
  #line 33 "/home/developer/.moon/lib/core/builtin/int.mbt"
  if (_M0L4selfS594 >= 0) {
    int32_t _M0L6_2atmpS1962;
    int32_t _M0L6_2atmpS1961;
    int32_t _M0L6_2atmpS1960;
    int32_t _M0L6_2atmpS1959;
    if (_M0L4selfS594 <= 1) {
      return 1;
    }
    if (_M0L4selfS594 > 1073741824) {
      return 1073741824;
    }
    _M0L6_2atmpS1962 = _M0L4selfS594 - 1;
    #line 44 "/home/developer/.moon/lib/core/builtin/int.mbt"
    _M0L6_2atmpS1961 = moonbit_clz32(_M0L6_2atmpS1962);
    _M0L6_2atmpS1960 = _M0L6_2atmpS1961 - 1;
    _M0L6_2atmpS1959 = 2147483647 >> (_M0L6_2atmpS1960 & 31);
    return _M0L6_2atmpS1959 + 1;
  } else {
    #line 34 "/home/developer/.moon/lib/core/builtin/int.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB21calc__grow__threshold(int32_t _M0L8capacityS593) {
  int32_t _M0L6_2atmpS1958;
  #line 610 "/home/developer/.moon/lib/core/builtin/linked_hash_map.mbt"
  _M0L6_2atmpS1958 = _M0L8capacityS593 * 13;
  return _M0L6_2atmpS1958 / 16;
}

struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGiEEE(
  struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L4selfS585
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS585 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGiEE* _M0L7_2aSomeS586 = _M0L4selfS585;
    if (_M0L7_2aSomeS586) {
      moonbit_incref(_M0L7_2aSomeS586);
    }
    return _M0L7_2aSomeS586;
  }
}

struct _M0TPB5EntryGssE* _M0MPC16option6Option6unwrapGRPB5EntryGssEE(
  struct _M0TPB5EntryGssE* _M0L4selfS587
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS587 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGssE* _M0L7_2aSomeS588 = _M0L4selfS587;
    if (_M0L7_2aSomeS588) {
      moonbit_incref(_M0L7_2aSomeS588);
    }
    return _M0L7_2aSomeS588;
  }
}

struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB5ArrayGsEEE(
  struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L4selfS589
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS589 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB5ArrayGsEE* _M0L7_2aSomeS590 = _M0L4selfS589;
    if (_M0L7_2aSomeS590) {
      moonbit_incref(_M0L7_2aSomeS590);
    }
    return _M0L7_2aSomeS590;
  }
}

struct _M0TPB5EntryGsRPB13StringBuilderE* _M0MPC16option6Option6unwrapGRPB5EntryGsRPB13StringBuilderEE(
  struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L4selfS591
) {
  #line 38 "/home/developer/.moon/lib/core/builtin/option.mbt"
  if (_M0L4selfS591 == 0) {
    #line 40 "/home/developer/.moon/lib/core/builtin/option.mbt"
    moonbit_panic();
  } else {
    struct _M0TPB5EntryGsRPB13StringBuilderE* _M0L7_2aSomeS592 =
      _M0L4selfS591;
    if (_M0L7_2aSomeS592) {
      moonbit_incref(_M0L7_2aSomeS592);
    }
    return _M0L7_2aSomeS592;
  }
}

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t _M0L4selfS584) {
  #line 35 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 36 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC13int3Int18to__string_2einner(_M0L4selfS584, 10);
}

int32_t _M0MPC16string6String9get__char(
  moonbit_string_t _M0L4selfS581,
  int32_t _M0L3idxS580
) {
  int32_t _if__result_3462;
  #line 1986 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  if (_M0L3idxS580 >= 0) {
    int32_t _M0L6_2atmpS1950 = Moonbit_array_length(_M0L4selfS581);
    _if__result_3462 = _M0L3idxS580 < _M0L6_2atmpS1950;
  } else {
    _if__result_3462 = 0;
  }
  if (_if__result_3462) {
    int32_t _M0L1cS582 = _M0L4selfS581[_M0L3idxS580];
    #line 1989 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L1cS582)) {
      int32_t _M0L6_2atmpS1951 = _M0L3idxS580 + 1;
      int32_t _M0L6_2atmpS1952 = Moonbit_array_length(_M0L4selfS581);
      if (_M0L6_2atmpS1951 < _M0L6_2atmpS1952) {
        int32_t _M0L6_2atmpS1956 = _M0L3idxS580 + 1;
        int32_t _M0L4nextS583 = _M0L4selfS581[_M0L6_2atmpS1956];
        #line 1992 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L4nextS583)) {
          int32_t _M0L6_2atmpS1954 = (int32_t)_M0L1cS582;
          int32_t _M0L6_2atmpS1955 = (int32_t)_M0L4nextS583;
          int32_t _M0L6_2atmpS1953;
          #line 1993 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0L6_2atmpS1953
          = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1954, _M0L6_2atmpS1955);
          return _M0L6_2atmpS1953;
        } else {
          return -1;
        }
      } else {
        return -1;
      }
    } else {
      #line 1997 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L1cS582)) {
        return -1;
      } else {
        int32_t _M0L6_2atmpS1957;
        #line 2000 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS1957
        = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS582);
        return _M0L6_2atmpS1957;
      }
    }
  } else {
    return -1;
  }
}

struct _M0TPC16string10StringView _M0MPC16string10StringView7replace(
  struct _M0TPC16string10StringView _M0L4selfS573,
  struct _M0TPC16string10StringView _M0L3oldS574,
  struct _M0TPC16string10StringView _M0L3newS579
) {
  int64_t _M0L7_2abindS572;
  #line 1459 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1464 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS572
  = _M0MPC16string10StringView4find(_M0L4selfS573, _M0L3oldS574);
  if (_M0L7_2abindS572 == 4294967296ll) {
    moonbit_incref(_M0L4selfS573.$0);
    return _M0L4selfS573;
  } else {
    int64_t _M0L7_2aSomeS575 = _M0L7_2abindS572;
    int32_t _M0L6_2aendS576 = (int32_t)_M0L7_2aSomeS575;
    struct _M0TPB13StringBuilder* _M0L7_2aselfS578;
    int64_t _M0L6_2atmpS1942;
    struct _M0TPC16string10StringView _M0L6_2atmpS1941;
    struct _M0TPB4IterGcE* _M0L6_2atmpS1940;
    struct _M0TPB4IterGcE* _M0L6_2atmpS1943;
    int32_t _M0L3endS1948;
    int32_t _M0L5startS1949;
    int32_t _M0L6_2atmpS1947;
    int32_t _M0L6_2atmpS1946;
    struct _M0TPC16string10StringView _M0L6_2atmpS1945;
    struct _M0TPB4IterGcE* _M0L6_2atmpS1944;
    struct _M0TPB13StringBuilder* _M0L6_2atmpS1939;
    moonbit_string_t _M0L7_2abindS577;
    int32_t _M0L6_2atmpS1938;
    #line 1466 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L7_2aselfS578 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
    _M0L6_2atmpS1942 = (int64_t)_M0L6_2aendS576;
    #line 1467 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1941
    = _M0MPC16string10StringView12view_2einner(_M0L4selfS573, 0, _M0L6_2atmpS1942);
    #line 1467 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1940 = _M0MPC16string10StringView4iter(_M0L6_2atmpS1941);
    moonbit_decref(_M0L6_2atmpS1941.$0);
    #line 1466 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0MPB13StringBuilder11write__iter(_M0L7_2aselfS578, _M0L6_2atmpS1940);
    moonbit_decref(_M0L6_2atmpS1940);
    #line 1468 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1943 = _M0MPC16string10StringView4iter(_M0L3newS579);
    #line 1466 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0MPB13StringBuilder11write__iter(_M0L7_2aselfS578, _M0L6_2atmpS1943);
    moonbit_decref(_M0L6_2atmpS1943);
    _M0L3endS1948 = _M0L3oldS574.$2;
    _M0L5startS1949 = _M0L3oldS574.$1;
    _M0L6_2atmpS1947 = _M0L3endS1948 - _M0L5startS1949;
    _M0L6_2atmpS1946 = _M0L6_2aendS576 + _M0L6_2atmpS1947;
    #line 1469 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1945
    = _M0MPC16string10StringView12view_2einner(_M0L4selfS573, _M0L6_2atmpS1946, 4294967296ll);
    #line 1469 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1944 = _M0MPC16string10StringView4iter(_M0L6_2atmpS1945);
    moonbit_decref(_M0L6_2atmpS1945.$0);
    #line 1466 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0MPB13StringBuilder11write__iter(_M0L7_2aselfS578, _M0L6_2atmpS1944);
    moonbit_decref(_M0L6_2atmpS1944);
    _M0L6_2atmpS1939 = _M0L7_2aselfS578;
    #line 1466 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L7_2abindS577 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1939);
    moonbit_decref(_M0L6_2atmpS1939);
    _M0L6_2atmpS1938 = Moonbit_array_length(_M0L7_2abindS577);
    return (struct _M0TPC16string10StringView){.$0 = _M0L7_2abindS577,
                                                 .$1 = 0,
                                                 .$2 = _M0L6_2atmpS1938};
  }
}

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPB4Iter9to__arrayGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS564
) {
  int64_t _M0L7_2abindS563;
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L6resultS565;
  #line 841 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS563 = _M0L4selfS564->$1;
  if (_M0L7_2abindS563 == 4294967296ll) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS1937 =
      (struct _M0TPC16string10StringView*)moonbit_empty_ref_valtype_array;
    _M0L6resultS565
    = (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_M0L6resultS565)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 95, 0);
    _M0L6resultS565->$0 = _M0L6_2atmpS1937;
    _M0L6resultS565->$1 = 0;
  } else {
    int64_t _M0L7_2aSomeS566 = _M0L7_2abindS563;
    int32_t _M0L4_2anS567 = (int32_t)_M0L7_2aSomeS566;
    #line 844 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6resultS565
    = _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(_M0L4_2anS567);
  }
  _2awhile_571:;
  while (1) {
    void* _M0L7_2abindS568;
    #line 847 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L7_2abindS568
    = _M0MPB4Iter4nextGRPC16string10StringViewE(_M0L4selfS564);
    switch (Moonbit_object_tag(_M0L7_2abindS568)) {
      case 1: {
        struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS569 =
          (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS568;
        struct _M0TPC16string10StringView _M0L8_2afieldS3134 =
          _M0L7_2aSomeS569->$0;
        int32_t _M0L6_2acntS3244 =
          Moonbit_rc_count(Moonbit_object_header(_M0L7_2aSomeS569));
        struct _M0TPC16string10StringView _M0L4_2axS570;
        if (_M0L6_2acntS3244 > 1) {
          int32_t _M0L11_2anew__cntS3245 = _M0L6_2acntS3244 - 1;
          Moonbit_set_rc_count(Moonbit_object_header(_M0L7_2aSomeS569), _M0L11_2anew__cntS3245);
          moonbit_incref(_M0L8_2afieldS3134.$0);
        } else if (_M0L6_2acntS3244 == 1) {
          #line 847 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
          moonbit_free(_M0L7_2aSomeS569);
        }
        _M0L4_2axS570 = _M0L8_2afieldS3134;
        #line 848 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
        _M0MPC15array5Array4pushGRPC16string10StringViewE(_M0L6resultS565, _M0L4_2axS570);
        moonbit_decref(_M0L4_2axS570.$0);
        goto _2awhile_571;
        break;
      }
      default: {
        moonbit_decref(_M0L7_2abindS568);
        break;
      }
    }
    break;
  }
  return _M0L6resultS565;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string6String5split(
  moonbit_string_t _M0L4selfS561,
  struct _M0TPC16string10StringView _M0L3sepS562
) {
  int32_t _M0L6_2atmpS1936;
  struct _M0TPC16string10StringView _M0L6_2atmpS1935;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_3464;
  #line 1168 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1936 = Moonbit_array_length(_M0L4selfS561);
  moonbit_incref(_M0L4selfS561);
  _M0L6_2atmpS1935
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS561, .$1 = 0, .$2 = _M0L6_2atmpS1936
  };
  #line 1169 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3464
  = _M0MPC16string10StringView5split(_M0L6_2atmpS1935, _M0L3sepS562);
  moonbit_decref(_M0L6_2atmpS1935.$0);
  return _result_3464;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPC16string10StringView5split(
  struct _M0TPC16string10StringView _M0L4selfS552,
  struct _M0TPC16string10StringView _M0L3sepS551
) {
  int32_t _M0L3endS1933;
  int32_t _M0L5startS1934;
  int32_t _M0L8sep__lenS550;
  void* _M0L4SomeS1932;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS554;
  struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__* _closure_3466;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS1922;
  struct _M0TPB4IterGRPC16string10StringViewE* _result_3467;
  #line 1139 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1933 = _M0L3sepS551.$2;
  _M0L5startS1934 = _M0L3sepS551.$1;
  _M0L8sep__lenS550 = _M0L3endS1933 - _M0L5startS1934;
  if (_M0L8sep__lenS550 == 0) {
    struct _M0TPB4IterGcE* _M0L6_2atmpS1917;
    struct _M0TWcERPC16string10StringView* _M0L6_2atmpS1918;
    struct _M0TPB4IterGRPC16string10StringViewE* _result_3465;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _M0L6_2atmpS1917 = _M0MPC16string10StringView4iter(_M0L4selfS552);
    _M0L6_2atmpS1918
    = (struct _M0TWcERPC16string10StringView*)&_M0MPC16string10StringView5splitC1919l1145$closure.data;
    #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    _result_3465
    = _M0MPB4Iter3mapGcRPC16string10StringViewE(_M0L6_2atmpS1917, _M0L6_2atmpS1918);
    moonbit_decref(_M0L6_2atmpS1917);
    moonbit_decref(_M0L6_2atmpS1918);
    return _result_3465;
  }
  moonbit_incref(_M0L4selfS552.$0);
  _M0L4SomeS1932
  = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
  Moonbit_object_header(_M0L4SomeS1932)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 101, 1);
  ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS1932)->$0
  = _M0L4selfS552;
  _M0L9remainingS554
  = (struct _M0TPB8MutLocalGORPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGORPC16string10StringViewE));
  Moonbit_object_header(_M0L9remainingS554)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 104, 0);
  _M0L9remainingS554->$0 = _M0L4SomeS1932;
  moonbit_incref(_M0L3sepS551.$0);
  _closure_3466
  = (struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__*)moonbit_malloc(sizeof(struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__));
  Moonbit_object_header(_closure_3466)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 107, 0);
  _closure_3466->code = &_M0MPC16string10StringView5splitC1923l1148;
  _closure_3466->$0 = _M0L9remainingS554;
  _closure_3466->$1 = _M0L3sepS551;
  _closure_3466->$2 = _M0L8sep__lenS550;
  _M0L6_2atmpS1922
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_3466;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3467
  = _M0MPB4Iter3newGRPC16string10StringViewE(_M0L6_2atmpS1922, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1922);
  return _result_3467;
}

void* _M0MPC16string10StringView5splitC1923l1148(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS1924
) {
  struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__* _M0L14_2acasted__envS1925;
  int32_t _M0L8sep__lenS550;
  struct _M0TPC16string10StringView _M0L3sepS551;
  struct _M0TPB8MutLocalGORPC16string10StringViewE* _M0L9remainingS554;
  void* _M0L7_2abindS555;
  #line 1148 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L14_2acasted__envS1925
  = (struct _M0R44StringView_3a_3asplit_2eanon__u1923__l1148__*)_M0L6_2aenvS1924;
  _M0L8sep__lenS550 = _M0L14_2acasted__envS1925->$2;
  _M0L3sepS551 = _M0L14_2acasted__envS1925->$1;
  _M0L9remainingS554 = _M0L14_2acasted__envS1925->$0;
  _M0L7_2abindS555 = _M0L9remainingS554->$0;
  switch (Moonbit_object_tag(_M0L7_2abindS555)) {
    case 1: {
      struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some* _M0L7_2aSomeS556 =
        (struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L7_2abindS555;
      struct _M0TPC16string10StringView _M0L7_2aviewS557 =
        _M0L7_2aSomeS556->$0;
      int64_t _M0L7_2abindS558;
      moonbit_incref(_M0L7_2aviewS557.$0);
      #line 1150 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L7_2abindS558
      = _M0MPC16string10StringView4find(_M0L7_2aviewS557, _M0L3sepS551);
      if (_M0L7_2abindS558 == 4294967296ll) {
        void* _M0L4NoneS1926 =
          (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
        void* _M0L6_2aoldS3135 = _M0L9remainingS554->$0;
        void* _block_3468;
        moonbit_decref(_M0L6_2aoldS3135);
        _M0L9remainingS554->$0 = _M0L4NoneS1926;
        _block_3468
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_3468)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 101, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3468)->$0
        = _M0L7_2aviewS557;
        return _block_3468;
      } else {
        int64_t _M0L7_2aSomeS559 = _M0L7_2abindS558;
        int32_t _M0L6_2aendS560 = (int32_t)_M0L7_2aSomeS559;
        int32_t _M0L6_2atmpS1929 = _M0L6_2aendS560 + _M0L8sep__lenS550;
        struct _M0TPC16string10StringView _M0L6_2atmpS1928;
        void* _M0L4SomeS1927;
        void* _M0L6_2aoldS3136;
        int64_t _M0L6_2atmpS1931;
        struct _M0TPC16string10StringView _M0L6_2atmpS1930;
        void* _block_3469;
        #line 1154 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS1928
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS557, _M0L6_2atmpS1929, 4294967296ll);
        _M0L4SomeS1927
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_M0L4SomeS1927)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 101, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_M0L4SomeS1927)->$0
        = _M0L6_2atmpS1928;
        _M0L6_2aoldS3136 = _M0L9remainingS554->$0;
        moonbit_decref(_M0L6_2aoldS3136);
        _M0L9remainingS554->$0 = _M0L4SomeS1927;
        _M0L6_2atmpS1931 = (int64_t)_M0L6_2aendS560;
        #line 1155 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
        _M0L6_2atmpS1930
        = _M0MPC16string10StringView12view_2einner(_M0L7_2aviewS557, 0, _M0L6_2atmpS1931);
        moonbit_decref(_M0L7_2aviewS557.$0);
        _block_3469
        = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
        Moonbit_object_header(_block_3469)->meta
        = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 101, 1);
        ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3469)->$0
        = _M0L6_2atmpS1930;
        return _block_3469;
      }
      break;
    }
    default: {
      return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
      break;
    }
  }
}

struct _M0TPC16string10StringView _M0MPC16string10StringView5splitC1919l1145(
  struct _M0TWcERPC16string10StringView* _M0L6_2aenvS1920,
  int32_t _M0L1cS553
) {
  moonbit_string_t _M0L6_2atmpS1921;
  struct _M0TPC16string10StringView _result_3470;
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1921 = _M0IPC14char4CharPB4Show10to__string(_M0L1cS553);
  #line 1145 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3470
  = _M0MPC16string6String12view_2einner(_M0L6_2atmpS1921, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1921);
  return _result_3470;
}

moonbit_string_t _M0IPC14char4CharPB4Show10to__string(int32_t _M0L4selfS549) {
  #line 446 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 447 "/home/developer/.moon/lib/core/builtin/char.mbt"
  return _M0FPB16char__to__string(_M0L4selfS549);
}

moonbit_string_t _M0FPB16char__to__string(int32_t _M0L4charS548) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS547;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1916;
  moonbit_string_t _result_3471;
  #line 452 "/home/developer/.moon/lib/core/builtin/char.mbt"
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0L7_2aselfS547 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS547, _M0L4charS548);
  _M0L6_2atmpS1916 = _M0L7_2aselfS547;
  #line 454 "/home/developer/.moon/lib/core/builtin/char.mbt"
  _result_3471 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1916);
  moonbit_decref(_M0L6_2atmpS1916);
  return _result_3471;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3mapGcRPC16string10StringViewE(
  struct _M0TPB4IterGcE* _M0L4selfS543,
  struct _M0TWcERPC16string10StringView* _M0L1fS546
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__* _closure_3472;
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2atmpS1910;
  int64_t _M0L10size__hintS1911;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_3473;
  #line 389 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  moonbit_incref(_M0L1fS546);
  moonbit_incref(_M0L4selfS543);
  _closure_3472
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__*)moonbit_malloc(sizeof(struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__));
  Moonbit_object_header(_closure_3472)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 111, 0);
  _closure_3472->code = &_M0MPB4Iter3mapGcRPC16string10StringViewEC1912l391;
  _closure_3472->$0 = _M0L1fS546;
  _closure_3472->$1 = _M0L4selfS543;
  _M0L6_2atmpS1910
  = (struct _M0TWERPC16option6OptionGRPC16string10StringViewE*)_closure_3472;
  _M0L10size__hintS1911 = _M0L4selfS543->$1;
  _block_3473
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_3473)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 115, 0);
  _block_3473->$0 = _M0L6_2atmpS1910;
  _block_3473->$1 = _M0L10size__hintS1911;
  return _block_3473;
}

void* _M0MPB4Iter3mapGcRPC16string10StringViewEC1912l391(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L6_2aenvS1913
) {
  struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__* _M0L14_2acasted__envS1914;
  struct _M0TPB4IterGcE* _M0L4selfS543;
  struct _M0TWcERPC16string10StringView* _M0L1fS546;
  int32_t _M0L7_2abindS542;
  #line 391 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L14_2acasted__envS1914
  = (struct _M0R97Iter_3a_3amap_7c_5bChar_2c_20moonbitlang_2fcore_2fstring_2fStringView_5d_7c_2eanon__u1912__l391__*)_M0L6_2aenvS1913;
  _M0L4selfS543 = _M0L14_2acasted__envS1914->$1;
  _M0L1fS546 = _M0L14_2acasted__envS1914->$0;
  #line 392 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2abindS542 = _M0MPB4Iter4nextGcE(_M0L4selfS543);
  if (_M0L7_2abindS542 == -1) {
    return (struct moonbit_object*)&moonbit_constant_constructor_0 + 1;
  } else {
    int32_t _M0L7_2aSomeS544 = _M0L7_2abindS542;
    int32_t _M0L4_2axS545 = _M0L7_2aSomeS544;
    struct _M0TPC16string10StringView _M0L6_2atmpS1915;
    void* _block_3474;
    #line 393 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
    _M0L6_2atmpS1915 = _M0L1fS546->code(_M0L1fS546, _M0L4_2axS545);
    _block_3474
    = (void*)moonbit_malloc(sizeof(struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some));
    Moonbit_object_header(_block_3474)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 101, 1);
    ((struct _M0DTPC16option6OptionGRPC16string10StringViewE4Some*)_block_3474)->$0
    = _M0L6_2atmpS1915;
    return _block_3474;
  }
}

struct _M0TPB4IterGcE* _M0MPC16string6String4iter(
  moonbit_string_t _M0L4selfS537
) {
  int32_t _M0L3lenS536;
  struct _M0TPB8MutLocalGiE* _M0L5indexS538;
  struct _M0R38String_3a_3aiter_2eanon__u1894__l256__* _closure_3475;
  struct _M0TWEOc* _M0L6_2atmpS1893;
  struct _M0TPB4IterGcE* _result_3476;
  #line 252 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L3lenS536 = Moonbit_array_length(_M0L4selfS537);
  _M0L5indexS538
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS538)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS538->$0 = 0;
  moonbit_incref(_M0L4selfS537);
  _closure_3475
  = (struct _M0R38String_3a_3aiter_2eanon__u1894__l256__*)moonbit_malloc(sizeof(struct _M0R38String_3a_3aiter_2eanon__u1894__l256__));
  Moonbit_object_header(_closure_3475)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 118, 0);
  _closure_3475->code = &_M0MPC16string6String4iterC1894l256;
  _closure_3475->$0 = _M0L5indexS538;
  _closure_3475->$1 = _M0L4selfS537;
  _closure_3475->$2 = _M0L3lenS536;
  _M0L6_2atmpS1893 = (struct _M0TWEOc*)_closure_3475;
  #line 256 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_3476 = _M0MPB4Iter3newGcE(_M0L6_2atmpS1893, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1893);
  return _result_3476;
}

int32_t _M0MPC16string6String4iterC1894l256(
  struct _M0TWEOc* _M0L6_2aenvS1895
) {
  struct _M0R38String_3a_3aiter_2eanon__u1894__l256__* _M0L14_2acasted__envS1896;
  int32_t _M0L3lenS536;
  moonbit_string_t _M0L4selfS537;
  struct _M0TPB8MutLocalGiE* _M0L5indexS538;
  int32_t _M0L3valS1897;
  #line 256 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L14_2acasted__envS1896
  = (struct _M0R38String_3a_3aiter_2eanon__u1894__l256__*)_M0L6_2aenvS1895;
  _M0L3lenS536 = _M0L14_2acasted__envS1896->$2;
  _M0L4selfS537 = _M0L14_2acasted__envS1896->$1;
  _M0L5indexS538 = _M0L14_2acasted__envS1896->$0;
  _M0L3valS1897 = _M0L5indexS538->$0;
  if (_M0L3valS1897 < _M0L3lenS536) {
    int32_t _M0L3valS1909 = _M0L5indexS538->$0;
    int32_t _M0L2c1S539 = _M0L4selfS537[_M0L3valS1909];
    int32_t _if__result_3477;
    int32_t _M0L3valS1907;
    int32_t _M0L6_2atmpS1906;
    int32_t _M0L6_2atmpS1908;
    #line 259 "/home/developer/.moon/lib/core/builtin/string.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S539)) {
      int32_t _M0L3valS1899 = _M0L5indexS538->$0;
      int32_t _M0L6_2atmpS1898 = _M0L3valS1899 + 1;
      _if__result_3477 = _M0L6_2atmpS1898 < _M0L3lenS536;
    } else {
      _if__result_3477 = 0;
    }
    if (_if__result_3477) {
      int32_t _M0L3valS1905 = _M0L5indexS538->$0;
      int32_t _M0L6_2atmpS1904 = _M0L3valS1905 + 1;
      int32_t _M0L2c2S540 = _M0L4selfS537[_M0L6_2atmpS1904];
      #line 261 "/home/developer/.moon/lib/core/builtin/string.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S540)) {
        int32_t _M0L6_2atmpS1902 = (int32_t)_M0L2c1S539;
        int32_t _M0L6_2atmpS1903 = (int32_t)_M0L2c2S540;
        int32_t _M0L1cS541;
        int32_t _M0L3valS1901;
        int32_t _M0L6_2atmpS1900;
        #line 262 "/home/developer/.moon/lib/core/builtin/string.mbt"
        _M0L1cS541
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1902, _M0L6_2atmpS1903);
        _M0L3valS1901 = _M0L5indexS538->$0;
        _M0L6_2atmpS1900 = _M0L3valS1901 + 2;
        _M0L5indexS538->$0 = _M0L6_2atmpS1900;
        return _M0L1cS541;
      }
    }
    _M0L3valS1907 = _M0L5indexS538->$0;
    _M0L6_2atmpS1906 = _M0L3valS1907 + 1;
    _M0L5indexS538->$0 = _M0L6_2atmpS1906;
    #line 269 "/home/developer/.moon/lib/core/builtin/string.mbt"
    _M0L6_2atmpS1908 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S539);
    return _M0L6_2atmpS1908;
  } else {
    return -1;
  }
}

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS524,
  moonbit_string_t _M0L5valueS526
) {
  int32_t _M0L3lenS1873;
  moonbit_string_t* _M0L6_2atmpS1875;
  int32_t _M0L6_2atmpS1874;
  int32_t _M0L6lengthS525;
  moonbit_string_t* _M0L3bufS1876;
  moonbit_string_t _M0L6_2aoldS3139;
  int32_t _M0L6_2atmpS1877;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1873 = _M0L4selfS524->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1875 = _M0MPC15array5Array6bufferGsE(_M0L4selfS524);
  _M0L6_2atmpS1874 = Moonbit_array_length(_M0L6_2atmpS1875);
  moonbit_decref(_M0L6_2atmpS1875);
  if (_M0L3lenS1873 == _M0L6_2atmpS1874) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGsE(_M0L4selfS524);
  }
  _M0L6lengthS525 = _M0L4selfS524->$1;
  _M0L3bufS1876 = _M0L4selfS524->$0;
  _M0L6_2aoldS3139 = (moonbit_string_t)_M0L3bufS1876[_M0L6lengthS525];
  moonbit_incref(_M0L5valueS526);
  moonbit_decref(_M0L6_2aoldS3139);
  _M0L3bufS1876[_M0L6lengthS525] = _M0L5valueS526;
  _M0L6_2atmpS1877 = _M0L6lengthS525 + 1;
  _M0L4selfS524->$1 = _M0L6_2atmpS1877;
  return 0;
}

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS527,
  struct _M0TUsiE* _M0L5valueS529
) {
  int32_t _M0L3lenS1878;
  struct _M0TUsiE** _M0L6_2atmpS1880;
  int32_t _M0L6_2atmpS1879;
  int32_t _M0L6lengthS528;
  struct _M0TUsiE** _M0L3bufS1881;
  struct _M0TUsiE* _M0L6_2aoldS3141;
  int32_t _M0L6_2atmpS1882;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1878 = _M0L4selfS527->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1880 = _M0MPC15array5Array6bufferGUsiEE(_M0L4selfS527);
  _M0L6_2atmpS1879 = Moonbit_array_length(_M0L6_2atmpS1880);
  moonbit_decref(_M0L6_2atmpS1880);
  if (_M0L3lenS1878 == _M0L6_2atmpS1879) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGUsiEE(_M0L4selfS527);
  }
  _M0L6lengthS528 = _M0L4selfS527->$1;
  _M0L3bufS1881 = _M0L4selfS527->$0;
  _M0L6_2aoldS3141 = (struct _M0TUsiE*)_M0L3bufS1881[_M0L6lengthS528];
  moonbit_incref(_M0L5valueS529);
  if (_M0L6_2aoldS3141) {
    moonbit_decref(_M0L6_2aoldS3141);
  }
  _M0L3bufS1881[_M0L6lengthS528] = _M0L5valueS529;
  _M0L6_2atmpS1882 = _M0L6lengthS528 + 1;
  _M0L4selfS527->$1 = _M0L6_2atmpS1882;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS530,
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L5valueS532
) {
  int32_t _M0L3lenS1883;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS1885;
  int32_t _M0L6_2atmpS1884;
  int32_t _M0L6lengthS531;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3bufS1886;
  struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3143;
  int32_t _M0L6_2atmpS1887;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1883 = _M0L4selfS530->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1885
  = _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS530);
  _M0L6_2atmpS1884 = Moonbit_array_length(_M0L6_2atmpS1885);
  moonbit_decref(_M0L6_2atmpS1885);
  if (_M0L3lenS1883 == _M0L6_2atmpS1884) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS530);
  }
  _M0L6lengthS531 = _M0L4selfS530->$1;
  _M0L3bufS1886 = _M0L4selfS530->$0;
  _M0L6_2aoldS3143
  = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3bufS1886[_M0L6lengthS531];
  moonbit_incref(_M0L5valueS532);
  if (_M0L6_2aoldS3143) {
    moonbit_decref(_M0L6_2aoldS3143);
  }
  _M0L3bufS1886[_M0L6lengthS531] = _M0L5valueS532;
  _M0L6_2atmpS1887 = _M0L6lengthS531 + 1;
  _M0L4selfS530->$1 = _M0L6_2atmpS1887;
  return 0;
}

int32_t _M0MPC15array5Array4pushGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS533,
  struct _M0TPC16string10StringView _M0L5valueS535
) {
  int32_t _M0L3lenS1888;
  struct _M0TPC16string10StringView* _M0L6_2atmpS1890;
  int32_t _M0L6_2atmpS1889;
  int32_t _M0L6lengthS534;
  struct _M0TPC16string10StringView* _M0L3bufS1891;
  struct _M0TPC16string10StringView _M0L6_2aoldS3145;
  int32_t _M0L6_2atmpS1892;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS1888 = _M0L4selfS533->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS1890
  = _M0MPC15array5Array6bufferGRPC16string10StringViewE(_M0L4selfS533);
  _M0L6_2atmpS1889 = Moonbit_array_length(_M0L6_2atmpS1890);
  moonbit_decref(_M0L6_2atmpS1890);
  if (_M0L3lenS1888 == _M0L6_2atmpS1889) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGRPC16string10StringViewE(_M0L4selfS533);
  }
  _M0L6lengthS534 = _M0L4selfS533->$1;
  _M0L3bufS1891 = _M0L4selfS533->$0;
  _M0L6_2aoldS3145 = _M0L3bufS1891[_M0L6lengthS534];
  moonbit_incref(_M0L5valueS535.$0);
  moonbit_decref(_M0L6_2aoldS3145.$0);
  _M0L3bufS1891[_M0L6lengthS534] = _M0L5valueS535;
  _M0L6_2atmpS1892 = _M0L6lengthS534 + 1;
  _M0L4selfS533->$1 = _M0L6_2atmpS1892;
  return 0;
}

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE* _M0L4selfS513) {
  int32_t _M0L8old__capS512;
  int32_t _M0L8new__capS514;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS512 = _M0L4selfS513->$1;
  if (_M0L8old__capS512 == 0) {
    _M0L8new__capS514 = 8;
  } else {
    _M0L8new__capS514 = _M0L8old__capS512 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGsE(_M0L4selfS513, _M0L8new__capS514);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS516
) {
  int32_t _M0L8old__capS515;
  int32_t _M0L8new__capS517;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS515 = _M0L4selfS516->$1;
  if (_M0L8old__capS515 == 0) {
    _M0L8new__capS517 = 8;
  } else {
    _M0L8new__capS517 = _M0L8old__capS515 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGUsiEE(_M0L4selfS516, _M0L8new__capS517);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS519
) {
  int32_t _M0L8old__capS518;
  int32_t _M0L8new__capS520;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS518 = _M0L4selfS519->$1;
  if (_M0L8old__capS518 == 0) {
    _M0L8new__capS520 = 8;
  } else {
    _M0L8new__capS520 = _M0L8old__capS518 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(_M0L4selfS519, _M0L8new__capS520);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS522
) {
  int32_t _M0L8old__capS521;
  int32_t _M0L8new__capS523;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS521 = _M0L4selfS522->$1;
  if (_M0L8old__capS521 == 0) {
    _M0L8new__capS523 = 8;
  } else {
    _M0L8new__capS523 = _M0L8old__capS521 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(_M0L4selfS522, _M0L8new__capS523);
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS489,
  int32_t _M0L13new__capacityS492
) {
  moonbit_string_t* _M0L8old__bufS488;
  int32_t _M0L8old__capS490;
  int32_t _M0L9copy__lenS491;
  moonbit_string_t* _M0L8new__bufS493;
  moonbit_string_t* _M0L6_2aoldS3147;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS488 = _M0L4selfS489->$0;
  _M0L8old__capS490 = Moonbit_array_length(_M0L8old__bufS488);
  if (_M0L8old__capS490 < _M0L13new__capacityS492) {
    _M0L9copy__lenS491 = _M0L8old__capS490;
  } else {
    _M0L9copy__lenS491 = _M0L13new__capacityS492;
  }
  moonbit_incref(_M0L8old__bufS488);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS493
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(_M0L8old__bufS488, _M0L13new__capacityS492, _M0L9copy__lenS491, 0, 0);
  moonbit_decref(_M0L8old__bufS488);
  _M0L6_2aoldS3147 = _M0L4selfS489->$0;
  moonbit_decref(_M0L6_2aoldS3147);
  _M0L4selfS489->$0 = _M0L8new__bufS493;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS495,
  int32_t _M0L13new__capacityS498
) {
  struct _M0TUsiE** _M0L8old__bufS494;
  int32_t _M0L8old__capS496;
  int32_t _M0L9copy__lenS497;
  struct _M0TUsiE** _M0L8new__bufS499;
  struct _M0TUsiE** _M0L6_2aoldS3149;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS494 = _M0L4selfS495->$0;
  _M0L8old__capS496 = Moonbit_array_length(_M0L8old__bufS494);
  if (_M0L8old__capS496 < _M0L13new__capacityS498) {
    _M0L9copy__lenS497 = _M0L8old__capS496;
  } else {
    _M0L9copy__lenS497 = _M0L13new__capacityS498;
  }
  moonbit_incref(_M0L8old__bufS494);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS499
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(_M0L8old__bufS494, _M0L13new__capacityS498, _M0L9copy__lenS497, 0, 0);
  moonbit_decref(_M0L8old__bufS494);
  _M0L6_2aoldS3149 = _M0L4selfS495->$0;
  moonbit_decref(_M0L6_2aoldS3149);
  _M0L4selfS495->$0 = _M0L8new__bufS499;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS501,
  int32_t _M0L13new__capacityS504
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8old__bufS500;
  int32_t _M0L8old__capS502;
  int32_t _M0L9copy__lenS503;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8new__bufS505;
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2aoldS3151;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS500 = _M0L4selfS501->$0;
  _M0L8old__capS502 = Moonbit_array_length(_M0L8old__bufS500);
  if (_M0L8old__capS502 < _M0L13new__capacityS504) {
    _M0L9copy__lenS503 = _M0L8old__capS502;
  } else {
    _M0L9copy__lenS503 = _M0L13new__capacityS504;
  }
  moonbit_incref(_M0L8old__bufS500);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS505
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(_M0L8old__bufS500, _M0L13new__capacityS504, _M0L9copy__lenS503, 0, 0);
  moonbit_decref(_M0L8old__bufS500);
  _M0L6_2aoldS3151 = _M0L4selfS501->$0;
  moonbit_decref(_M0L6_2aoldS3151);
  _M0L4selfS501->$0 = _M0L8new__bufS505;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS507,
  int32_t _M0L13new__capacityS510
) {
  struct _M0TPC16string10StringView* _M0L8old__bufS506;
  int32_t _M0L8old__capS508;
  int32_t _M0L9copy__lenS509;
  struct _M0TPC16string10StringView* _M0L8new__bufS511;
  struct _M0TPC16string10StringView* _M0L6_2aoldS3153;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS506 = _M0L4selfS507->$0;
  _M0L8old__capS508 = Moonbit_array_length(_M0L8old__bufS506);
  if (_M0L8old__capS508 < _M0L13new__capacityS510) {
    _M0L9copy__lenS509 = _M0L8old__capS508;
  } else {
    _M0L9copy__lenS509 = _M0L13new__capacityS510;
  }
  moonbit_incref(_M0L8old__bufS506);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS511
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(_M0L8old__bufS506, _M0L13new__capacityS510, _M0L9copy__lenS509, 0, 0);
  moonbit_decref(_M0L8old__bufS506);
  _M0L6_2aoldS3153 = _M0L4selfS507->$0;
  moonbit_decref(_M0L6_2aoldS3153);
  _M0L4selfS507->$0 = _M0L8new__bufS511;
  return 0;
}

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq20MultipleSeqAlignmentE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L4selfS484
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS484->$1;
}

int32_t _M0MPC15array5Array6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS485
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS485->$1;
}

int32_t _M0MPC15array5Array6lengthGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS486
) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS486->$1;
}

int32_t _M0MPC15array5Array6lengthGsE(struct _M0TPB5ArrayGsE* _M0L4selfS487) {
  #line 140 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  return _M0L4selfS487->$1;
}

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(
  int32_t _M0L8capacityS481
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS481 == 0) {
    moonbit_string_t* _M0L6_2atmpS1867 =
      (moonbit_string_t*)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGsE* _block_3478 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_3478)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
    _block_3478->$0 = _M0L6_2atmpS1867;
    _block_3478->$1 = 0;
    return _block_3478;
  } else {
    moonbit_string_t* _M0L6_2atmpS1868 =
      (moonbit_string_t*)moonbit_make_ref_array(_M0L8capacityS481, (moonbit_string_t)moonbit_string_literal_0.data);
    struct _M0TPB5ArrayGsE* _block_3479 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_3479)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 0);
    _block_3479->$0 = _M0L6_2atmpS1868;
    _block_3479->$1 = 0;
    return _block_3479;
  }
}

struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0MPC15array5Array11new_2einnerGRP26biolab8bio__seq9SeqRecordE(
  int32_t _M0L8capacityS482
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS482 == 0) {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS1869 =
      (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _block_3480 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
    Moonbit_object_header(_block_3480)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
    _block_3480->$0 = _M0L6_2atmpS1869;
    _block_3480->$1 = 0;
    return _block_3480;
  } else {
    struct _M0TP26biolab8bio__seq9SeqRecord** _M0L6_2atmpS1870 =
      (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array(_M0L8capacityS482, 0);
    struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _block_3481 =
      (struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE));
    Moonbit_object_header(_block_3481)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 18, 0);
    _block_3481->$0 = _M0L6_2atmpS1870;
    _block_3481->$1 = 0;
    return _block_3481;
  }
}

struct _M0TPB5ArrayGRPC16string10StringViewE* _M0MPC15array5Array11new_2einnerGRPC16string10StringViewE(
  int32_t _M0L8capacityS483
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS483 == 0) {
    struct _M0TPC16string10StringView* _M0L6_2atmpS1871 =
      (struct _M0TPC16string10StringView*)moonbit_empty_ref_valtype_array;
    struct _M0TPB5ArrayGRPC16string10StringViewE* _block_3482 =
      (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_block_3482)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 95, 0);
    _block_3482->$0 = _M0L6_2atmpS1871;
    _block_3482->$1 = 0;
    return _block_3482;
  } else {
    struct _M0TPC16string10StringView* _M0L6_2atmpS1872 =
      (struct _M0TPC16string10StringView*)moonbit_make_ref_valtype_array(_M0L8capacityS483, sizeof(struct _M0TPC16string10StringView), Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 98, 0), &(struct _M0TPC16string10StringView){.$0 = (moonbit_string_t)moonbit_string_literal_0.data, .$1 = 0, .$2 = 0});
    struct _M0TPB5ArrayGRPC16string10StringViewE* _block_3483 =
      (struct _M0TPB5ArrayGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGRPC16string10StringViewE));
    Moonbit_object_header(_block_3483)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 95, 0);
    _block_3483->$0 = _M0L6_2atmpS1872;
    _block_3483->$1 = 0;
    return _block_3483;
  }
}

int32_t _M0MPC16string6String11has__prefix(
  moonbit_string_t _M0L4selfS479,
  struct _M0TPC16string10StringView _M0L3strS480
) {
  int32_t _M0L6_2atmpS1866;
  struct _M0TPC16string10StringView _M0L6_2atmpS1865;
  int32_t _result_3484;
  #line 297 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1866 = Moonbit_array_length(_M0L4selfS479);
  moonbit_incref(_M0L4selfS479);
  _M0L6_2atmpS1865
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS479, .$1 = 0, .$2 = _M0L6_2atmpS1866
  };
  #line 299 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3484
  = _M0MPC16string10StringView11has__prefix(_M0L6_2atmpS1865, _M0L3strS480);
  moonbit_decref(_M0L6_2atmpS1865.$0);
  return _result_3484;
}

int32_t _M0MPC16string10StringView11has__prefix(
  struct _M0TPC16string10StringView _M0L4selfS475,
  struct _M0TPC16string10StringView _M0L3strS476
) {
  int64_t _M0L7_2abindS474;
  #line 290 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  #line 292 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L7_2abindS474
  = _M0MPC16string10StringView4find(_M0L4selfS475, _M0L3strS476);
  if (_M0L7_2abindS474 == 4294967296ll) {
    return 0;
  } else {
    int64_t _M0L7_2aSomeS477 = _M0L7_2abindS474;
    int32_t _M0L4_2aiS478 = (int32_t)_M0L7_2aSomeS477;
    return _M0L4_2aiS478 == 0;
  }
}

moonbit_string_t _M0MPC16string6String6repeat(
  moonbit_string_t _M0L4selfS467,
  int32_t _M0L1nS466
) {
  #line 1049 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  if (_M0L1nS466 < 0) {
    #line 1051 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_69.data);
  } else if (_M0L1nS466 == 0) {
    return (moonbit_string_t)moonbit_string_literal_0.data;
  } else if (_M0L1nS466 == 1) {
    moonbit_incref(_M0L4selfS467);
    return _M0L4selfS467;
  } else {
    int32_t _M0L3lenS468 = Moonbit_array_length(_M0L4selfS467);
    int32_t _M0L5totalS469 = _M0L3lenS468 * _M0L1nS466;
    int32_t _if__result_3485;
    if (_M0L3lenS468 == 0) {
      _if__result_3485 = 1;
    } else {
      int32_t _M0L6_2atmpS1863 = _M0L5totalS469 / _M0L1nS466;
      _if__result_3485 = _M0L6_2atmpS1863 == _M0L3lenS468;
    }
    if (_if__result_3485) {
      struct _M0TPB13StringBuilder* _M0L3bufS470;
      moonbit_string_t _M0L3strS471;
      int32_t _M0L2__S472;
      moonbit_string_t _result_3487;
      #line 1060 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3bufS470
      = _M0MPB13StringBuilder21StringBuilder_2einner(_M0L5totalS469);
      #line 1061 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _M0L3strS471 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS467);
      _M0L2__S472 = 0;
      while (1) {
        if (_M0L2__S472 < _M0L1nS466) {
          int32_t _M0L6_2atmpS1864;
          #line 1063 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
          _M0IPB13StringBuilderPB6Logger13write__string(_M0L3bufS470, _M0L3strS471);
          _M0L6_2atmpS1864 = _M0L2__S472 + 1;
          _M0L2__S472 = _M0L6_2atmpS1864;
          continue;
        } else {
          moonbit_decref(_M0L3strS471);
        }
        break;
      }
      #line 1065 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      _result_3487 = _M0MPB13StringBuilder10to__string(_M0L3bufS470);
      moonbit_decref(_M0L3bufS470);
      return _result_3487;
    } else {
      #line 1058 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
      return _M0FPC15abort5abortGsE((moonbit_string_t)moonbit_string_literal_70.data);
    }
  }
}

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(
  moonbit_string_t _M0L4selfS465
) {
  #line 222 "/home/developer/.moon/lib/core/builtin/show.mbt"
  moonbit_incref(_M0L4selfS465);
  return _M0L4selfS465;
}

int64_t _M0MPC16string6String4find(
  moonbit_string_t _M0L4selfS463,
  struct _M0TPC16string10StringView _M0L3strS464
) {
  int32_t _M0L6_2atmpS1862;
  struct _M0TPC16string10StringView _M0L6_2atmpS1861;
  int64_t _result_3488;
  #line 101 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L6_2atmpS1862 = Moonbit_array_length(_M0L4selfS463);
  moonbit_incref(_M0L4selfS463);
  _M0L6_2atmpS1861
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS463, .$1 = 0, .$2 = _M0L6_2atmpS1862
  };
  #line 102 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _result_3488
  = _M0MPC16string10StringView4find(_M0L6_2atmpS1861, _M0L3strS464);
  moonbit_decref(_M0L6_2atmpS1861.$0);
  return _result_3488;
}

int64_t _M0MPC16string10StringView4find(
  struct _M0TPC16string10StringView _M0L4selfS462,
  struct _M0TPC16string10StringView _M0L3strS461
) {
  int32_t _M0L3endS1859;
  int32_t _M0L5startS1860;
  int32_t _M0L6_2atmpS1858;
  #line 18 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1859 = _M0L3strS461.$2;
  _M0L5startS1860 = _M0L3strS461.$1;
  _M0L6_2atmpS1858 = _M0L3endS1859 - _M0L5startS1860;
  if (_M0L6_2atmpS1858 <= 4) {
    #line 20 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB18brute__force__find(_M0L4selfS462, _M0L3strS461);
  } else {
    #line 22 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
    return _M0FPB28boyer__moore__horspool__find(_M0L4selfS462, _M0L3strS461);
  }
}

int64_t _M0FPB18brute__force__find(
  struct _M0TPC16string10StringView _M0L8haystackS451,
  struct _M0TPC16string10StringView _M0L6needleS453
) {
  int32_t _M0L3endS1856;
  int32_t _M0L5startS1857;
  int32_t _M0L13haystack__lenS450;
  int32_t _M0L3endS1854;
  int32_t _M0L5startS1855;
  int32_t _M0L11needle__lenS452;
  #line 31 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1856 = _M0L8haystackS451.$2;
  _M0L5startS1857 = _M0L8haystackS451.$1;
  _M0L13haystack__lenS450 = _M0L3endS1856 - _M0L5startS1857;
  _M0L3endS1854 = _M0L6needleS453.$2;
  _M0L5startS1855 = _M0L6needleS453.$1;
  _M0L11needle__lenS452 = _M0L3endS1854 - _M0L5startS1855;
  if (_M0L11needle__lenS452 > 0) {
    if (_M0L13haystack__lenS450 >= _M0L11needle__lenS452) {
      moonbit_string_t _M0L3strS1852 = _M0L6needleS453.$0;
      int32_t _M0L5startS1853 = _M0L6needleS453.$1;
      int32_t _M0L13needle__firstS454 = _M0L3strS1852[_M0L5startS1853];
      int32_t _M0L12forward__lenS455 =
        _M0L13haystack__lenS450 - _M0L11needle__lenS452;
      int32_t _M0L1iS456 = 0;
      while (1) {
        if (_M0L1iS456 <= _M0L12forward__lenS455) {
          moonbit_string_t _M0L3strS1839 = _M0L8haystackS451.$0;
          int32_t _M0L5startS1841 = _M0L8haystackS451.$1;
          int32_t _M0L6_2atmpS1840 = _M0L5startS1841 + _M0L1iS456;
          int32_t _M0L6_2atmpS1838 = _M0L3strS1839[_M0L6_2atmpS1840];
          int32_t _M0L1jS459;
          int32_t _M0L6_2atmpS1837;
          if (_M0L6_2atmpS1838 != _M0L13needle__firstS454) {
            goto join_457;
          }
          _M0L1jS459 = 1;
          while (1) {
            if (_M0L1jS459 < _M0L11needle__lenS452) {
              moonbit_string_t _M0L3strS1847 = _M0L8haystackS451.$0;
              int32_t _M0L5startS1849 = _M0L8haystackS451.$1;
              int32_t _M0L6_2atmpS1850 = _M0L1iS456 + _M0L1jS459;
              int32_t _M0L6_2atmpS1848 = _M0L5startS1849 + _M0L6_2atmpS1850;
              int32_t _M0L6_2atmpS1842 = _M0L3strS1847[_M0L6_2atmpS1848];
              moonbit_string_t _M0L3strS1844 = _M0L6needleS453.$0;
              int32_t _M0L5startS1846 = _M0L6needleS453.$1;
              int32_t _M0L6_2atmpS1845 = _M0L5startS1846 + _M0L1jS459;
              int32_t _M0L6_2atmpS1843 = _M0L3strS1844[_M0L6_2atmpS1845];
              int32_t _M0L6_2atmpS1851;
              if (_M0L6_2atmpS1842 != _M0L6_2atmpS1843) {
                break;
              }
              _M0L6_2atmpS1851 = _M0L1jS459 + 1;
              _M0L1jS459 = _M0L6_2atmpS1851;
              continue;
            } else {
              return (int64_t)_M0L1iS456;
            }
            break;
          }
          goto join_457;
          goto joinlet_3490;
          join_457:;
          _M0L6_2atmpS1837 = _M0L1iS456 + 1;
          _M0L1iS456 = _M0L6_2atmpS1837;
          continue;
          joinlet_3490:;
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
  struct _M0TPC16string10StringView _M0L8haystackS438,
  struct _M0TPC16string10StringView _M0L6needleS440
) {
  int32_t _M0L3endS1835;
  int32_t _M0L5startS1836;
  int32_t _M0L13haystack__lenS437;
  int32_t _M0L3endS1833;
  int32_t _M0L5startS1834;
  int32_t _M0L11needle__lenS439;
  #line 57 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
  _M0L3endS1835 = _M0L8haystackS438.$2;
  _M0L5startS1836 = _M0L8haystackS438.$1;
  _M0L13haystack__lenS437 = _M0L3endS1835 - _M0L5startS1836;
  _M0L3endS1833 = _M0L6needleS440.$2;
  _M0L5startS1834 = _M0L6needleS440.$1;
  _M0L11needle__lenS439 = _M0L3endS1833 - _M0L5startS1834;
  if (_M0L11needle__lenS439 > 0) {
    if (_M0L13haystack__lenS437 >= _M0L11needle__lenS439) {
      int32_t* _M0L11skip__tableS441 =
        (int32_t*)moonbit_make_int32_array(256, _M0L11needle__lenS439);
      int32_t _M0L7_2abindS442 = _M0L11needle__lenS439 - 1;
      int32_t _M0L1iS443 = 0;
      int32_t _M0L1iS445;
      while (1) {
        if (_M0L1iS443 < _M0L7_2abindS442) {
          moonbit_string_t _M0L3strS1808 = _M0L6needleS440.$0;
          int32_t _M0L5startS1810 = _M0L6needleS440.$1;
          int32_t _M0L6_2atmpS1809 = _M0L5startS1810 + _M0L1iS443;
          int32_t _M0L6_2atmpS1807 = _M0L3strS1808[_M0L6_2atmpS1809];
          int32_t _M0L6_2atmpS1806 = (int32_t)_M0L6_2atmpS1807;
          int32_t _M0L6_2atmpS1803 = _M0L6_2atmpS1806 & 255;
          int32_t _M0L6_2atmpS1805 = _M0L11needle__lenS439 - 1;
          int32_t _M0L6_2atmpS1804 = _M0L6_2atmpS1805 - _M0L1iS443;
          int32_t _M0L6_2atmpS1811;
          if (
            _M0L6_2atmpS1803 < 0
            || _M0L6_2atmpS1803
               >= Moonbit_array_length(_M0L11skip__tableS441)
          ) {
            #line 68 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L11skip__tableS441[_M0L6_2atmpS1803] = _M0L6_2atmpS1804;
          _M0L6_2atmpS1811 = _M0L1iS443 + 1;
          _M0L1iS443 = _M0L6_2atmpS1811;
          continue;
        }
        break;
      }
      _M0L1iS445 = 0;
      while (1) {
        int32_t _M0L6_2atmpS1812 =
          _M0L13haystack__lenS437 - _M0L11needle__lenS439;
        if (_M0L1iS445 <= _M0L6_2atmpS1812) {
          int32_t _M0L7_2abindS446 = _M0L11needle__lenS439 - 1;
          int32_t _M0L1jS447 = 0;
          moonbit_string_t _M0L3strS1828;
          int32_t _M0L5startS1830;
          int32_t _M0L6_2atmpS1832;
          int32_t _M0L6_2atmpS1831;
          int32_t _M0L6_2atmpS1829;
          int32_t _M0L6_2atmpS1827;
          int32_t _M0L6_2atmpS1826;
          int32_t _M0L6_2atmpS1825;
          int32_t _M0L6_2atmpS1824;
          int32_t _M0L6_2atmpS1823;
          while (1) {
            if (_M0L1jS447 <= _M0L7_2abindS446) {
              moonbit_string_t _M0L3strS1818 = _M0L8haystackS438.$0;
              int32_t _M0L5startS1820 = _M0L8haystackS438.$1;
              int32_t _M0L6_2atmpS1821 = _M0L1iS445 + _M0L1jS447;
              int32_t _M0L6_2atmpS1819 = _M0L5startS1820 + _M0L6_2atmpS1821;
              int32_t _M0L6_2atmpS1813 = _M0L3strS1818[_M0L6_2atmpS1819];
              moonbit_string_t _M0L3strS1815 = _M0L6needleS440.$0;
              int32_t _M0L5startS1817 = _M0L6needleS440.$1;
              int32_t _M0L6_2atmpS1816 = _M0L5startS1817 + _M0L1jS447;
              int32_t _M0L6_2atmpS1814 = _M0L3strS1815[_M0L6_2atmpS1816];
              int32_t _M0L6_2atmpS1822;
              if (_M0L6_2atmpS1813 != _M0L6_2atmpS1814) {
                break;
              }
              _M0L6_2atmpS1822 = _M0L1jS447 + 1;
              _M0L1jS447 = _M0L6_2atmpS1822;
              continue;
            } else {
              moonbit_decref(_M0L11skip__tableS441);
              return (int64_t)_M0L1iS445;
            }
            break;
          }
          _M0L3strS1828 = _M0L8haystackS438.$0;
          _M0L5startS1830 = _M0L8haystackS438.$1;
          _M0L6_2atmpS1832 = _M0L1iS445 + _M0L11needle__lenS439;
          _M0L6_2atmpS1831 = _M0L6_2atmpS1832 - 1;
          _M0L6_2atmpS1829 = _M0L5startS1830 + _M0L6_2atmpS1831;
          _M0L6_2atmpS1827 = _M0L3strS1828[_M0L6_2atmpS1829];
          _M0L6_2atmpS1826 = (int32_t)_M0L6_2atmpS1827;
          _M0L6_2atmpS1825 = _M0L6_2atmpS1826 & 255;
          if (
            _M0L6_2atmpS1825 < 0
            || _M0L6_2atmpS1825
               >= Moonbit_array_length(_M0L11skip__tableS441)
          ) {
            #line 73 "/home/developer/.moon/lib/core/builtin/string_methods.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS1824 = (int32_t)_M0L11skip__tableS441[_M0L6_2atmpS1825];
          _M0L6_2atmpS1823 = _M0L1iS445 + _M0L6_2atmpS1824;
          _M0L1iS445 = _M0L6_2atmpS1823;
          continue;
        } else {
          moonbit_decref(_M0L11skip__tableS441);
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
  struct _M0TPB13StringBuilder* _M0L4selfS436,
  struct _M0TPC16string10StringView _M0L3strS435
) {
  int32_t _M0L3endS1801;
  int32_t _M0L5startS1802;
  int32_t _M0L8str__lenS434;
  int32_t _M0L3lenS1794;
  int32_t _M0L6_2atmpS1793;
  uint16_t* _M0L4dataS1795;
  int32_t _M0L3lenS1796;
  moonbit_string_t _M0L6_2atmpS1797;
  int32_t _M0L6_2atmpS1798;
  int32_t _M0L3lenS1800;
  int32_t _M0L6_2atmpS1799;
  #line 131 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3endS1801 = _M0L3strS435.$2;
  _M0L5startS1802 = _M0L3strS435.$1;
  _M0L8str__lenS434 = _M0L3endS1801 - _M0L5startS1802;
  _M0L3lenS1794 = _M0L4selfS436->$1;
  _M0L6_2atmpS1793 = _M0L3lenS1794 + _M0L8str__lenS434;
  #line 136 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS436, _M0L6_2atmpS1793);
  _M0L4dataS1795 = _M0L4selfS436->$0;
  _M0L3lenS1796 = _M0L4selfS436->$1;
  moonbit_incref(_M0L4dataS1795);
  #line 139 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1797 = _M0MPC16string10StringView4data(_M0L3strS435);
  #line 140 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1798 = _M0MPC16string10StringView13start__offset(_M0L3strS435);
  #line 137 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1795, _M0L3lenS1796, _M0L6_2atmpS1797, _M0L6_2atmpS1798, _M0L8str__lenS434);
  moonbit_decref(_M0L4dataS1795);
  moonbit_decref(_M0L6_2atmpS1797);
  _M0L3lenS1800 = _M0L4selfS436->$1;
  _M0L6_2atmpS1799 = _M0L3lenS1800 + _M0L8str__lenS434;
  _M0L4selfS436->$1 = _M0L6_2atmpS1799;
  return 0;
}

int32_t _M0MPB13StringBuilder11write__iter(
  struct _M0TPB13StringBuilder* _M0L4selfS432,
  struct _M0TPB4IterGcE* _M0L4iterS429
) {
  #line 44 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  while (1) {
    int32_t _M0L7_2abindS428;
    #line 48 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
    _M0L7_2abindS428 = _M0MPB4Iter4nextGcE(_M0L4iterS429);
    if (_M0L7_2abindS428 == -1) {
      
    } else {
      int32_t _M0L7_2aSomeS430 = _M0L7_2abindS428;
      int32_t _M0L5_2achS431 = _M0L7_2aSomeS430;
      #line 49 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
      _M0IPB13StringBuilderPB6Logger11write__char(_M0L4selfS432, _M0L5_2achS431);
      continue;
    }
    break;
  }
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String12view_2einner(
  moonbit_string_t _M0L4selfS426,
  int32_t _M0L13start__offsetS427,
  int64_t _M0L11end__offsetS424
) {
  int32_t _M0L11end__offsetS423;
  int32_t _if__result_3496;
  #line 614 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS424 == 4294967296ll) {
    _M0L11end__offsetS423 = Moonbit_array_length(_M0L4selfS426);
  } else {
    int64_t _M0L7_2aSomeS425 = _M0L11end__offsetS424;
    _M0L11end__offsetS423 = (int32_t)_M0L7_2aSomeS425;
  }
  if (_M0L13start__offsetS427 >= 0) {
    if (_M0L13start__offsetS427 <= _M0L11end__offsetS423) {
      int32_t _M0L6_2atmpS1792 = Moonbit_array_length(_M0L4selfS426);
      _if__result_3496 = _M0L11end__offsetS423 <= _M0L6_2atmpS1792;
    } else {
      _if__result_3496 = 0;
    }
  } else {
    _if__result_3496 = 0;
  }
  if (_if__result_3496) {
    moonbit_incref(_M0L4selfS426);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS426,
                                                 .$1 = _M0L13start__offsetS427,
                                                 .$2 = _M0L11end__offsetS423};
  } else {
    #line 623 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_71.data);
  }
}

struct _M0TPB4IterGcE* _M0MPC16string10StringView4iter(
  struct _M0TPC16string10StringView _M0L4selfS418
) {
  int32_t _M0L5startS417;
  int32_t _M0L3endS419;
  struct _M0TPB8MutLocalGiE* _M0L5indexS420;
  struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__* _closure_3497;
  struct _M0TWEOc* _M0L6_2atmpS1771;
  struct _M0TPB4IterGcE* _result_3498;
  #line 214 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L5startS417 = _M0L4selfS418.$1;
  _M0L3endS419 = _M0L4selfS418.$2;
  _M0L5indexS420
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5indexS420)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5indexS420->$0 = _M0L5startS417;
  moonbit_incref(_M0L4selfS418.$0);
  _closure_3497
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__*)moonbit_malloc(sizeof(struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__));
  Moonbit_object_header(_closure_3497)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 122, 0);
  _closure_3497->code = &_M0MPC16string10StringView4iterC1772l219;
  _closure_3497->$0 = _M0L5indexS420;
  _closure_3497->$1 = _M0L3endS419;
  _closure_3497->$2 = _M0L4selfS418;
  _M0L6_2atmpS1771 = (struct _M0TWEOc*)_closure_3497;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_3498 = _M0MPB4Iter3newGcE(_M0L6_2atmpS1771, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1771);
  return _result_3498;
}

int32_t _M0MPC16string10StringView4iterC1772l219(
  struct _M0TWEOc* _M0L6_2aenvS1773
) {
  struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__* _M0L14_2acasted__envS1774;
  struct _M0TPC16string10StringView _M0L4selfS418;
  int32_t _M0L3endS419;
  struct _M0TPB8MutLocalGiE* _M0L5indexS420;
  int32_t _M0L3valS1775;
  #line 219 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L14_2acasted__envS1774
  = (struct _M0R42StringView_3a_3aiter_2eanon__u1772__l219__*)_M0L6_2aenvS1773;
  _M0L4selfS418 = _M0L14_2acasted__envS1774->$2;
  _M0L3endS419 = _M0L14_2acasted__envS1774->$1;
  _M0L5indexS420 = _M0L14_2acasted__envS1774->$0;
  _M0L3valS1775 = _M0L5indexS420->$0;
  if (_M0L3valS1775 < _M0L3endS419) {
    moonbit_string_t _M0L3strS1790 = _M0L4selfS418.$0;
    int32_t _M0L3valS1791 = _M0L5indexS420->$0;
    int32_t _M0L2c1S421 = _M0L3strS1790[_M0L3valS1791];
    int32_t _if__result_3499;
    int32_t _M0L3valS1788;
    int32_t _M0L6_2atmpS1787;
    int32_t _M0L6_2atmpS1789;
    #line 222 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S421)) {
      int32_t _M0L3valS1778 = _M0L5indexS420->$0;
      int32_t _M0L6_2atmpS1776 = _M0L3valS1778 + 1;
      int32_t _M0L3endS1777 = _M0L4selfS418.$2;
      _if__result_3499 = _M0L6_2atmpS1776 < _M0L3endS1777;
    } else {
      _if__result_3499 = 0;
    }
    if (_if__result_3499) {
      moonbit_string_t _M0L3strS1784 = _M0L4selfS418.$0;
      int32_t _M0L3valS1786 = _M0L5indexS420->$0;
      int32_t _M0L6_2atmpS1785 = _M0L3valS1786 + 1;
      int32_t _M0L2c2S422 = _M0L3strS1784[_M0L6_2atmpS1785];
      #line 224 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S422)) {
        int32_t _M0L3valS1780 = _M0L5indexS420->$0;
        int32_t _M0L6_2atmpS1779 = _M0L3valS1780 + 2;
        int32_t _M0L6_2atmpS1782;
        int32_t _M0L6_2atmpS1783;
        int32_t _M0L6_2atmpS1781;
        _M0L5indexS420->$0 = _M0L6_2atmpS1779;
        _M0L6_2atmpS1782 = (int32_t)_M0L2c1S421;
        _M0L6_2atmpS1783 = (int32_t)_M0L2c2S422;
        #line 226 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        _M0L6_2atmpS1781
        = _M0FPB32code__point__of__surrogate__pair(_M0L6_2atmpS1782, _M0L6_2atmpS1783);
        return _M0L6_2atmpS1781;
      }
    }
    _M0L3valS1788 = _M0L5indexS420->$0;
    _M0L6_2atmpS1787 = _M0L3valS1788 + 1;
    _M0L5indexS420->$0 = _M0L6_2atmpS1787;
    #line 230 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    _M0L6_2atmpS1789 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L2c1S421);
    return _M0L6_2atmpS1789;
  } else {
    return -1;
  }
}

moonbit_string_t _M0MPC16string10StringView9to__owned(
  struct _M0TPC16string10StringView _M0L4selfS416
) {
  moonbit_string_t _M0L3strS1768;
  int32_t _M0L5startS1769;
  int32_t _M0L3endS1770;
  moonbit_string_t _result_3500;
  #line 196 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1768 = _M0L4selfS416.$0;
  _M0L5startS1769 = _M0L4selfS416.$1;
  _M0L3endS1770 = _M0L4selfS416.$2;
  moonbit_incref(_M0L3strS1768);
  #line 199 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _result_3500
  = _M0MPC16string6String17unsafe__substring(_M0L3strS1768, _M0L5startS1769, _M0L3endS1770);
  moonbit_decref(_M0L3strS1768);
  return _result_3500;
}

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t _M0L3strS413,
  int32_t _M0L5startS411,
  int32_t _M0L3endS412
) {
  int32_t _if__result_3501;
  int32_t _M0L3lenS414;
  int32_t _M0L6_2atmpS1766;
  int32_t _M0L6_2atmpS1767;
  moonbit_bytes_t _M0L5bytesS415;
  moonbit_bytes_t _M0L6_2atmpS1765;
  moonbit_string_t _result_3502;
  #line 91 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L5startS411 == 0) {
    int32_t _M0L6_2atmpS1764 = Moonbit_array_length(_M0L3strS413);
    _if__result_3501 = _M0L3endS412 == _M0L6_2atmpS1764;
  } else {
    _if__result_3501 = 0;
  }
  if (_if__result_3501) {
    moonbit_incref(_M0L3strS413);
    return _M0L3strS413;
  }
  _M0L3lenS414 = _M0L3endS412 - _M0L5startS411;
  _M0L6_2atmpS1766 = _M0L3lenS414 * 2;
  #line 101 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1767 = _M0IPC14byte4BytePB7Default7default();
  _M0L5bytesS415
  = (moonbit_bytes_t)moonbit_make_bytes(_M0L6_2atmpS1766, _M0L6_2atmpS1767);
  #line 102 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0MPC15array10FixedArray18blit__from__string(_M0L5bytesS415, 0, _M0L3strS413, _M0L5startS411, _M0L3lenS414);
  _M0L6_2atmpS1765 = _M0L5bytesS415;
  #line 103 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_3502
  = _M0MPC15bytes5Bytes29to__unchecked__string_2einner(_M0L6_2atmpS1765, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS1765);
  return _result_3502;
}

int32_t _M0IPC14byte4BytePB7Default7default() {
  #line 231 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  return 0;
}

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t _M0L4selfS406,
  int32_t _M0L6offsetS410,
  int64_t _M0L6lengthS408
) {
  int32_t _M0L3lenS405;
  int32_t _M0L6lengthS407;
  int32_t _if__result_3503;
  #line 77 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L3lenS405 = Moonbit_array_length(_M0L4selfS406);
  if (_M0L6lengthS408 == 4294967296ll) {
    _M0L6lengthS407 = _M0L3lenS405 - _M0L6offsetS410;
  } else {
    int64_t _M0L7_2aSomeS409 = _M0L6lengthS408;
    _M0L6lengthS407 = (int32_t)_M0L7_2aSomeS409;
  }
  if (_M0L6offsetS410 >= 0) {
    if (_M0L6lengthS407 >= 0) {
      int32_t _M0L6_2atmpS1763 = _M0L6offsetS410 + _M0L6lengthS407;
      _if__result_3503 = _M0L6_2atmpS1763 <= _M0L3lenS405;
    } else {
      _if__result_3503 = 0;
    }
  } else {
    _if__result_3503 = 0;
  }
  if (_if__result_3503) {
    moonbit_incref(_M0L4selfS406);
    #line 85 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    return _M0FPB19unsafe__sub__string(_M0L4selfS406, _M0L6offsetS410, _M0L6lengthS407);
  } else {
    #line 84 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t _M0L4selfS397,
  int32_t _M0L13bytes__offsetS392,
  moonbit_string_t _M0L3strS399,
  int32_t _M0L11str__offsetS395,
  int32_t _M0L6lengthS393
) {
  int32_t _M0L6_2atmpS1762;
  int32_t _M0L6_2atmpS1761;
  int32_t _M0L2e1S391;
  int32_t _M0L6_2atmpS1760;
  int32_t _M0L2e2S394;
  int32_t _M0L4len1S396;
  int32_t _M0L4len2S398;
  int32_t _if__result_3504;
  #line 125 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L6_2atmpS1762 = _M0L6lengthS393 * 2;
  _M0L6_2atmpS1761 = _M0L13bytes__offsetS392 + _M0L6_2atmpS1762;
  _M0L2e1S391 = _M0L6_2atmpS1761 - 1;
  _M0L6_2atmpS1760 = _M0L11str__offsetS395 + _M0L6lengthS393;
  _M0L2e2S394 = _M0L6_2atmpS1760 - 1;
  _M0L4len1S396 = Moonbit_array_length(_M0L4selfS397);
  _M0L4len2S398 = Moonbit_array_length(_M0L3strS399);
  if (_M0L6lengthS393 >= 0) {
    if (_M0L13bytes__offsetS392 >= 0) {
      if (_M0L2e1S391 < _M0L4len1S396) {
        if (_M0L11str__offsetS395 >= 0) {
          _if__result_3504 = _M0L2e2S394 < _M0L4len2S398;
        } else {
          _if__result_3504 = 0;
        }
      } else {
        _if__result_3504 = 0;
      }
    } else {
      _if__result_3504 = 0;
    }
  } else {
    _if__result_3504 = 0;
  }
  if (_if__result_3504) {
    int32_t _M0L16end__str__offsetS400 =
      _M0L11str__offsetS395 + _M0L6lengthS393;
    int32_t _M0L1iS401 = _M0L11str__offsetS395;
    int32_t _M0L1jS402 = _M0L13bytes__offsetS392;
    while (1) {
      if (_M0L1iS401 < _M0L16end__str__offsetS400) {
        int32_t _M0L6_2atmpS1757 = _M0L3strS399[_M0L1iS401];
        int32_t _M0L6_2atmpS1756 = (int32_t)_M0L6_2atmpS1757;
        uint32_t _M0L1cS403 = *(uint32_t*)&_M0L6_2atmpS1756;
        uint32_t _M0L6_2atmpS1752 = _M0L1cS403 & 255u;
        int32_t _M0L6_2atmpS1751;
        int32_t _M0L6_2atmpS1753;
        uint32_t _M0L6_2atmpS1755;
        int32_t _M0L6_2atmpS1754;
        int32_t _M0L6_2atmpS1758;
        int32_t _M0L6_2atmpS1759;
        #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1751 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1752);
        if (
          _M0L1jS402 < 0 || _M0L1jS402 >= Moonbit_array_length(_M0L4selfS397)
        ) {
          #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS397[_M0L1jS402] = _M0L6_2atmpS1751;
        _M0L6_2atmpS1753 = _M0L1jS402 + 1;
        _M0L6_2atmpS1755 = _M0L1cS403 >> 8;
        #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS1754 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS1755);
        if (
          _M0L6_2atmpS1753 < 0
          || _M0L6_2atmpS1753 >= Moonbit_array_length(_M0L4selfS397)
        ) {
          #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS397[_M0L6_2atmpS1753] = _M0L6_2atmpS1754;
        _M0L6_2atmpS1758 = _M0L1iS401 + 1;
        _M0L6_2atmpS1759 = _M0L1jS402 + 2;
        _M0L1iS401 = _M0L6_2atmpS1758;
        _M0L1jS402 = _M0L6_2atmpS1759;
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

int32_t _M0MPC14uint4UInt8to__byte(uint32_t _M0L4selfS390) {
  int32_t _M0L6_2atmpS1750;
  #line 2519 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1750 = *(int32_t*)&_M0L4selfS390;
  return _M0L6_2atmpS1750 & 0xff;
}

moonbit_string_t* _M0MPC15array5Array6bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS385
) {
  moonbit_string_t* _M0L8_2afieldS3167;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3167 = _M0L4selfS385->$0;
  moonbit_incref(_M0L8_2afieldS3167);
  return _M0L8_2afieldS3167;
}

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS386
) {
  struct _M0TUsiE** _M0L8_2afieldS3168;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3168 = _M0L4selfS386->$0;
  moonbit_incref(_M0L8_2afieldS3168);
  return _M0L8_2afieldS3168;
}

struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0MPC15array5Array6bufferGRP26biolab8bio__seq20MultipleSeqAlignmentE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq20MultipleSeqAlignmentE* _M0L4selfS387
) {
  struct _M0TP26biolab8bio__seq20MultipleSeqAlignment** _M0L8_2afieldS3169;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3169 = _M0L4selfS387->$0;
  moonbit_incref(_M0L8_2afieldS3169);
  return _M0L8_2afieldS3169;
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPC15array5Array6bufferGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TPB5ArrayGRP26biolab8bio__seq9SeqRecordE* _M0L4selfS388
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L8_2afieldS3170;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3170 = _M0L4selfS388->$0;
  moonbit_incref(_M0L8_2afieldS3170);
  return _M0L8_2afieldS3170;
}

struct _M0TPC16string10StringView* _M0MPC15array5Array6bufferGRPC16string10StringViewE(
  struct _M0TPB5ArrayGRPC16string10StringViewE* _M0L4selfS389
) {
  struct _M0TPC16string10StringView* _M0L8_2afieldS3171;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS3171 = _M0L4selfS389->$0;
  moonbit_incref(_M0L8_2afieldS3171);
  return _M0L8_2afieldS3171;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView12view_2einner(
  struct _M0TPC16string10StringView _M0L4selfS383,
  int32_t _M0L13start__offsetS384,
  int64_t _M0L11end__offsetS381
) {
  int32_t _M0L11end__offsetS380;
  int32_t _if__result_3506;
  #line 105 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  if (_M0L11end__offsetS381 == 4294967296ll) {
    int32_t _M0L3endS1748 = _M0L4selfS383.$2;
    int32_t _M0L5startS1749 = _M0L4selfS383.$1;
    _M0L11end__offsetS380 = _M0L3endS1748 - _M0L5startS1749;
  } else {
    int64_t _M0L7_2aSomeS382 = _M0L11end__offsetS381;
    _M0L11end__offsetS380 = (int32_t)_M0L7_2aSomeS382;
  }
  if (_M0L13start__offsetS384 >= 0) {
    if (_M0L13start__offsetS384 <= _M0L11end__offsetS380) {
      int32_t _M0L3endS1741 = _M0L4selfS383.$2;
      int32_t _M0L5startS1742 = _M0L4selfS383.$1;
      int32_t _M0L6_2atmpS1740 = _M0L3endS1741 - _M0L5startS1742;
      _if__result_3506 = _M0L11end__offsetS380 <= _M0L6_2atmpS1740;
    } else {
      _if__result_3506 = 0;
    }
  } else {
    _if__result_3506 = 0;
  }
  if (_if__result_3506) {
    moonbit_string_t _M0L3strS1743 = _M0L4selfS383.$0;
    int32_t _M0L5startS1747 = _M0L4selfS383.$1;
    int32_t _M0L6_2atmpS1744 = _M0L5startS1747 + _M0L13start__offsetS384;
    int32_t _M0L5startS1746 = _M0L4selfS383.$1;
    int32_t _M0L6_2atmpS1745 = _M0L5startS1746 + _M0L11end__offsetS380;
    moonbit_incref(_M0L3strS1743);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1743,
                                                 .$1 = _M0L6_2atmpS1744,
                                                 .$2 = _M0L6_2atmpS1745};
  } else {
    #line 114 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    return _M0FPC15abort5abortGRPC16string10StringViewE((moonbit_string_t)moonbit_string_literal_71.data);
  }
}

struct _M0TPB4IterGcE* _M0MPB4Iter3newGcE(
  struct _M0TWEOc* _M0L1fS369,
  int64_t _M0L10size__hintS366
) {
  int64_t _M0L10size__hintS365;
  struct _M0TPB4IterGcE* _block_3507;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS366 == 4294967296ll) {
    _M0L10size__hintS365 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS367 = _M0L10size__hintS366;
    int32_t _M0L4_2anS368 = (int32_t)_M0L7_2aSomeS367;
    if (_M0L4_2anS368 > 0) {
      _M0L10size__hintS365 = (int64_t)_M0L4_2anS368;
    } else {
      _M0L10size__hintS365 = _M0MPB4Iter3newN6constrS9988GcE;
    }
  }
  moonbit_incref(_M0L1fS369);
  _block_3507
  = (struct _M0TPB4IterGcE*)moonbit_malloc(sizeof(struct _M0TPB4IterGcE));
  Moonbit_object_header(_block_3507)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 126, 0);
  _block_3507->$0 = _M0L1fS369;
  _block_3507->$1 = _M0L10size__hintS365;
  return _block_3507;
}

struct _M0TPB4IterGRPC16string10StringViewE* _M0MPB4Iter3newGRPC16string10StringViewE(
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L1fS374,
  int64_t _M0L10size__hintS371
) {
  int64_t _M0L10size__hintS370;
  struct _M0TPB4IterGRPC16string10StringViewE* _block_3508;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS371 == 4294967296ll) {
    _M0L10size__hintS370 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS372 = _M0L10size__hintS371;
    int32_t _M0L4_2anS373 = (int32_t)_M0L7_2aSomeS372;
    if (_M0L4_2anS373 > 0) {
      _M0L10size__hintS370 = (int64_t)_M0L4_2anS373;
    } else {
      _M0L10size__hintS370
      = _M0MPB4Iter3newN6constrS9988GRPC16string10StringViewE;
    }
  }
  moonbit_incref(_M0L1fS374);
  _block_3508
  = (struct _M0TPB4IterGRPC16string10StringViewE*)moonbit_malloc(sizeof(struct _M0TPB4IterGRPC16string10StringViewE));
  Moonbit_object_header(_block_3508)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 115, 0);
  _block_3508->$0 = _M0L1fS374;
  _block_3508->$1 = _M0L10size__hintS370;
  return _block_3508;
}

struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0MPB4Iter3newGUsRPB13StringBuilderEE(
  struct _M0TWEOUsRPB13StringBuilderE* _M0L1fS379,
  int64_t _M0L10size__hintS376
) {
  int64_t _M0L10size__hintS375;
  struct _M0TPB4IterGUsRPB13StringBuilderEE* _block_3509;
  #line 270 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  if (_M0L10size__hintS376 == 4294967296ll) {
    _M0L10size__hintS375 = 4294967296ll;
  } else {
    int64_t _M0L7_2aSomeS377 = _M0L10size__hintS376;
    int32_t _M0L4_2anS378 = (int32_t)_M0L7_2aSomeS377;
    if (_M0L4_2anS378 > 0) {
      _M0L10size__hintS375 = (int64_t)_M0L4_2anS378;
    } else {
      _M0L10size__hintS375
      = _M0MPB4Iter3newN6constrS9988GUsRPB13StringBuilderEE;
    }
  }
  moonbit_incref(_M0L1fS379);
  _block_3509
  = (struct _M0TPB4IterGUsRPB13StringBuilderEE*)moonbit_malloc(sizeof(struct _M0TPB4IterGUsRPB13StringBuilderEE));
  Moonbit_object_header(_block_3509)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 129, 0);
  _block_3509->$0 = _M0L1fS379;
  _block_3509->$1 = _M0L10size__hintS375;
  return _block_3509;
}

moonbit_string_t _M0MPC13int3Int18to__string_2einner(
  int32_t _M0L4selfS349,
  int32_t _M0L5radixS348
) {
  int32_t _if__result_3510;
  int32_t _M0L12is__negativeS350;
  uint32_t _M0L3numS351;
  uint16_t* _M0L6bufferS352;
  #line 209 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS348 < 2) {
    _if__result_3510 = 1;
  } else {
    _if__result_3510 = _M0L5radixS348 > 36;
  }
  if (_if__result_3510) {
    #line 213 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_72.data);
  }
  if (_M0L4selfS349 == 0) {
    return (moonbit_string_t)moonbit_string_literal_73.data;
  }
  _M0L12is__negativeS350 = _M0L4selfS349 < 0;
  if (_M0L12is__negativeS350) {
    int32_t _M0L6_2atmpS1739 = -_M0L4selfS349;
    _M0L3numS351 = *(uint32_t*)&_M0L6_2atmpS1739;
  } else {
    _M0L3numS351 = *(uint32_t*)&_M0L4selfS349;
  }
  switch (_M0L5radixS348) {
    case 10: {
      int32_t _M0L10digit__lenS353;
      int32_t _M0L6_2atmpS1736;
      int32_t _M0L10total__lenS354;
      uint16_t* _M0L6bufferS355;
      int32_t _M0L12digit__startS356;
      #line 235 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS353 = _M0FPB12dec__count32(_M0L3numS351);
      if (_M0L12is__negativeS350) {
        _M0L6_2atmpS1736 = 1;
      } else {
        _M0L6_2atmpS1736 = 0;
      }
      _M0L10total__lenS354 = _M0L10digit__lenS353 + _M0L6_2atmpS1736;
      _M0L6bufferS355
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS354, 0);
      if (_M0L12is__negativeS350) {
        _M0L12digit__startS356 = 1;
      } else {
        _M0L12digit__startS356 = 0;
      }
      #line 239 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__dec(_M0L6bufferS355, _M0L3numS351, _M0L12digit__startS356, _M0L10total__lenS354);
      _M0L6bufferS352 = _M0L6bufferS355;
      break;
    }
    
    case 16: {
      int32_t _M0L10digit__lenS357;
      int32_t _M0L6_2atmpS1737;
      int32_t _M0L10total__lenS358;
      uint16_t* _M0L6bufferS359;
      int32_t _M0L12digit__startS360;
      #line 243 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS357 = _M0FPB12hex__count32(_M0L3numS351);
      if (_M0L12is__negativeS350) {
        _M0L6_2atmpS1737 = 1;
      } else {
        _M0L6_2atmpS1737 = 0;
      }
      _M0L10total__lenS358 = _M0L10digit__lenS357 + _M0L6_2atmpS1737;
      _M0L6bufferS359
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS358, 0);
      if (_M0L12is__negativeS350) {
        _M0L12digit__startS360 = 1;
      } else {
        _M0L12digit__startS360 = 0;
      }
      #line 247 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__hex(_M0L6bufferS359, _M0L3numS351, _M0L12digit__startS360, _M0L10total__lenS358);
      _M0L6bufferS352 = _M0L6bufferS359;
      break;
    }
    default: {
      int32_t _M0L10digit__lenS361;
      int32_t _M0L6_2atmpS1738;
      int32_t _M0L10total__lenS362;
      uint16_t* _M0L6bufferS363;
      int32_t _M0L12digit__startS364;
      #line 251 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS361
      = _M0FPB14radix__count32(_M0L3numS351, _M0L5radixS348);
      if (_M0L12is__negativeS350) {
        _M0L6_2atmpS1738 = 1;
      } else {
        _M0L6_2atmpS1738 = 0;
      }
      _M0L10total__lenS362 = _M0L10digit__lenS361 + _M0L6_2atmpS1738;
      _M0L6bufferS363
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS362, 0);
      if (_M0L12is__negativeS350) {
        _M0L12digit__startS364 = 1;
      } else {
        _M0L12digit__startS364 = 0;
      }
      #line 255 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB24int__to__string__generic(_M0L6bufferS363, _M0L3numS351, _M0L12digit__startS364, _M0L10total__lenS362, _M0L5radixS348);
      _M0L6bufferS352 = _M0L6bufferS363;
      break;
    }
  }
  if (_M0L12is__negativeS350) {
    _M0L6bufferS352[0] = 45;
  }
  return _M0L6bufferS352;
}

int32_t _M0FPB14radix__count32(
  uint32_t _M0L5valueS342,
  int32_t _M0L5radixS344
) {
  uint32_t _M0L4baseS343;
  uint32_t _M0L3numS345;
  int32_t _M0L5countS346;
  #line 189 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS342 == 0u) {
    return 1;
  }
  _M0L4baseS343 = *(uint32_t*)&_M0L5radixS344;
  _M0L3numS345 = _M0L5valueS342;
  _M0L5countS346 = 0;
  while (1) {
    if (_M0L3numS345 > 0u) {
      uint32_t _M0L6_2atmpS1734 = _M0L3numS345 / _M0L4baseS343;
      int32_t _M0L6_2atmpS1735 = _M0L5countS346 + 1;
      _M0L3numS345 = _M0L6_2atmpS1734;
      _M0L5countS346 = _M0L6_2atmpS1735;
      continue;
    } else {
      return _M0L5countS346;
    }
    break;
  }
}

int32_t _M0FPB12hex__count32(uint32_t _M0L5valueS340) {
  #line 177 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS340 == 0u) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS341;
    int32_t _M0L6_2atmpS1733;
    int32_t _M0L6_2atmpS1732;
    #line 182 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS341 = moonbit_clz32(_M0L5valueS340);
    _M0L6_2atmpS1733 = 31 - _M0L14leading__zerosS341;
    _M0L6_2atmpS1732 = _M0L6_2atmpS1733 / 4;
    return _M0L6_2atmpS1732 + 1;
  }
}

int32_t _M0FPB12dec__count32(uint32_t _M0L5valueS339) {
  #line 143 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS339 >= 100000u) {
    if (_M0L5valueS339 >= 10000000u) {
      if (_M0L5valueS339 >= 1000000000u) {
        return 10;
      } else if (_M0L5valueS339 >= 100000000u) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS339 >= 1000000u) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS339 >= 1000u) {
    if (_M0L5valueS339 >= 10000u) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS339 >= 100u) {
    return 3;
  } else if (_M0L5valueS339 >= 10u) {
    return 2;
  } else {
    return 1;
  }
}

int32_t _M0FPB20int__to__string__dec(
  uint16_t* _M0L6bufferS325,
  uint32_t _M0L3numS337,
  int32_t _M0L12digit__startS326,
  int32_t _M0L10total__lenS338
) {
  int32_t _M0L6_2atmpS1731;
  uint32_t _M0L3numS315;
  int32_t _M0L6offsetS316;
  #line 88 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1731 = _M0L10total__lenS338 - _M0L12digit__startS326;
  _M0L3numS315 = _M0L3numS337;
  _M0L6offsetS316 = _M0L6_2atmpS1731;
  while (1) {
    if (_M0L3numS315 >= 10000u) {
      uint32_t _M0L1tS317 = _M0L3numS315 / 10000u;
      uint32_t _M0L6_2atmpS1708 = _M0L3numS315 % 10000u;
      int32_t _M0L1rS318 = *(int32_t*)&_M0L6_2atmpS1708;
      int32_t _M0L2d1S319 = _M0L1rS318 / 100;
      int32_t _M0L2d2S320 = _M0L1rS318 % 100;
      int32_t _M0L6_2atmpS1707 = _M0L2d1S319 / 10;
      int32_t _M0L6_2atmpS1706 = 48 + _M0L6_2atmpS1707;
      int32_t _M0L6d1__hiS321 = (uint16_t)_M0L6_2atmpS1706;
      int32_t _M0L6_2atmpS1705 = _M0L2d1S319 % 10;
      int32_t _M0L6_2atmpS1704 = 48 + _M0L6_2atmpS1705;
      int32_t _M0L6d1__loS322 = (uint16_t)_M0L6_2atmpS1704;
      int32_t _M0L6_2atmpS1703 = _M0L2d2S320 / 10;
      int32_t _M0L6_2atmpS1702 = 48 + _M0L6_2atmpS1703;
      int32_t _M0L6d2__hiS323 = (uint16_t)_M0L6_2atmpS1702;
      int32_t _M0L6_2atmpS1701 = _M0L2d2S320 % 10;
      int32_t _M0L6_2atmpS1700 = 48 + _M0L6_2atmpS1701;
      int32_t _M0L6d2__loS324 = (uint16_t)_M0L6_2atmpS1700;
      int32_t _M0L6_2atmpS1692 = _M0L12digit__startS326 + _M0L6offsetS316;
      int32_t _M0L6_2atmpS1691 = _M0L6_2atmpS1692 - 4;
      int32_t _M0L6_2atmpS1694;
      int32_t _M0L6_2atmpS1693;
      int32_t _M0L6_2atmpS1696;
      int32_t _M0L6_2atmpS1695;
      int32_t _M0L6_2atmpS1698;
      int32_t _M0L6_2atmpS1697;
      int32_t _M0L6_2atmpS1699;
      _M0L6bufferS325[_M0L6_2atmpS1691] = _M0L6d1__hiS321;
      _M0L6_2atmpS1694 = _M0L12digit__startS326 + _M0L6offsetS316;
      _M0L6_2atmpS1693 = _M0L6_2atmpS1694 - 3;
      _M0L6bufferS325[_M0L6_2atmpS1693] = _M0L6d1__loS322;
      _M0L6_2atmpS1696 = _M0L12digit__startS326 + _M0L6offsetS316;
      _M0L6_2atmpS1695 = _M0L6_2atmpS1696 - 2;
      _M0L6bufferS325[_M0L6_2atmpS1695] = _M0L6d2__hiS323;
      _M0L6_2atmpS1698 = _M0L12digit__startS326 + _M0L6offsetS316;
      _M0L6_2atmpS1697 = _M0L6_2atmpS1698 - 1;
      _M0L6bufferS325[_M0L6_2atmpS1697] = _M0L6d2__loS324;
      _M0L6_2atmpS1699 = _M0L6offsetS316 - 4;
      _M0L3numS315 = _M0L1tS317;
      _M0L6offsetS316 = _M0L6_2atmpS1699;
      continue;
    } else {
      int32_t _M0L6_2atmpS1730 = *(int32_t*)&_M0L3numS315;
      int32_t _M0L9remainingS328 = _M0L6_2atmpS1730;
      int32_t _M0L6offsetS329 = _M0L6offsetS316;
      while (1) {
        if (_M0L9remainingS328 >= 100) {
          int32_t _M0L1tS330 = _M0L9remainingS328 / 100;
          int32_t _M0L1dS331 = _M0L9remainingS328 % 100;
          int32_t _M0L6_2atmpS1717 = _M0L1dS331 / 10;
          int32_t _M0L6_2atmpS1716 = 48 + _M0L6_2atmpS1717;
          int32_t _M0L5d__hiS332 = (uint16_t)_M0L6_2atmpS1716;
          int32_t _M0L6_2atmpS1715 = _M0L1dS331 % 10;
          int32_t _M0L6_2atmpS1714 = 48 + _M0L6_2atmpS1715;
          int32_t _M0L5d__loS333 = (uint16_t)_M0L6_2atmpS1714;
          int32_t _M0L6_2atmpS1710 = _M0L12digit__startS326 + _M0L6offsetS329;
          int32_t _M0L6_2atmpS1709 = _M0L6_2atmpS1710 - 2;
          int32_t _M0L6_2atmpS1712;
          int32_t _M0L6_2atmpS1711;
          int32_t _M0L6_2atmpS1713;
          _M0L6bufferS325[_M0L6_2atmpS1709] = _M0L5d__hiS332;
          _M0L6_2atmpS1712 = _M0L12digit__startS326 + _M0L6offsetS329;
          _M0L6_2atmpS1711 = _M0L6_2atmpS1712 - 1;
          _M0L6bufferS325[_M0L6_2atmpS1711] = _M0L5d__loS333;
          _M0L6_2atmpS1713 = _M0L6offsetS329 - 2;
          _M0L9remainingS328 = _M0L1tS330;
          _M0L6offsetS329 = _M0L6_2atmpS1713;
          continue;
        } else if (_M0L9remainingS328 >= 10) {
          int32_t _M0L6_2atmpS1725 = _M0L9remainingS328 / 10;
          int32_t _M0L6_2atmpS1724 = 48 + _M0L6_2atmpS1725;
          int32_t _M0L5d__hiS335 = (uint16_t)_M0L6_2atmpS1724;
          int32_t _M0L6_2atmpS1723 = _M0L9remainingS328 % 10;
          int32_t _M0L6_2atmpS1722 = 48 + _M0L6_2atmpS1723;
          int32_t _M0L5d__loS336 = (uint16_t)_M0L6_2atmpS1722;
          int32_t _M0L6_2atmpS1719 = _M0L12digit__startS326 + _M0L6offsetS329;
          int32_t _M0L6_2atmpS1718 = _M0L6_2atmpS1719 - 2;
          int32_t _M0L6_2atmpS1721;
          int32_t _M0L6_2atmpS1720;
          _M0L6bufferS325[_M0L6_2atmpS1718] = _M0L5d__hiS335;
          _M0L6_2atmpS1721 = _M0L12digit__startS326 + _M0L6offsetS329;
          _M0L6_2atmpS1720 = _M0L6_2atmpS1721 - 1;
          _M0L6bufferS325[_M0L6_2atmpS1720] = _M0L5d__loS336;
        } else {
          int32_t _M0L6_2atmpS1729 = _M0L12digit__startS326 + _M0L6offsetS329;
          int32_t _M0L6_2atmpS1726 = _M0L6_2atmpS1729 - 1;
          int32_t _M0L6_2atmpS1728 = 48 + _M0L9remainingS328;
          int32_t _M0L6_2atmpS1727 = (uint16_t)_M0L6_2atmpS1728;
          _M0L6bufferS325[_M0L6_2atmpS1726] = _M0L6_2atmpS1727;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB24int__to__string__generic(
  uint16_t* _M0L6bufferS305,
  uint32_t _M0L3numS309,
  int32_t _M0L12digit__startS306,
  int32_t _M0L10total__lenS308,
  int32_t _M0L5radixS299
) {
  uint32_t _M0L4baseS298;
  int32_t _M0L6_2atmpS1676;
  int32_t _M0L6_2atmpS1675;
  #line 57 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS298 = *(uint32_t*)&_M0L5radixS299;
  _M0L6_2atmpS1676 = _M0L5radixS299 - 1;
  _M0L6_2atmpS1675 = _M0L5radixS299 & _M0L6_2atmpS1676;
  if (_M0L6_2atmpS1675 == 0) {
    int32_t _M0L5shiftS300;
    uint32_t _M0L4maskS301;
    int32_t _M0L6_2atmpS1683;
    int32_t _M0L6offsetS302;
    uint32_t _M0L1nS303;
    #line 68 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS300 = moonbit_ctz32(_M0L5radixS299);
    _M0L4maskS301 = _M0L4baseS298 - 1u;
    _M0L6_2atmpS1683 = _M0L10total__lenS308 - _M0L12digit__startS306;
    _M0L6offsetS302 = _M0L6_2atmpS1683;
    _M0L1nS303 = _M0L3numS309;
    while (1) {
      if (_M0L1nS303 > 0u) {
        uint32_t _M0L6_2atmpS1682 = _M0L1nS303 & _M0L4maskS301;
        int32_t _M0L5digitS304 = *(int32_t*)&_M0L6_2atmpS1682;
        int32_t _M0L6_2atmpS1679 = _M0L12digit__startS306 + _M0L6offsetS302;
        int32_t _M0L6_2atmpS1677 = _M0L6_2atmpS1679 - 1;
        int32_t _M0L6_2atmpS1678 =
          ((moonbit_string_t)moonbit_string_literal_74.data)[_M0L5digitS304];
        int32_t _M0L6_2atmpS1680;
        uint32_t _M0L6_2atmpS1681;
        _M0L6bufferS305[_M0L6_2atmpS1677] = _M0L6_2atmpS1678;
        _M0L6_2atmpS1680 = _M0L6offsetS302 - 1;
        _M0L6_2atmpS1681 = _M0L1nS303 >> (_M0L5shiftS300 & 31);
        _M0L6offsetS302 = _M0L6_2atmpS1680;
        _M0L1nS303 = _M0L6_2atmpS1681;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1690 = _M0L10total__lenS308 - _M0L12digit__startS306;
    int32_t _M0L6offsetS310 = _M0L6_2atmpS1690;
    uint32_t _M0L1nS311 = _M0L3numS309;
    while (1) {
      if (_M0L1nS311 > 0u) {
        uint32_t _M0L1qS312 = _M0L1nS311 / _M0L4baseS298;
        uint32_t _M0L6_2atmpS1689 = _M0L1qS312 * _M0L4baseS298;
        uint32_t _M0L6_2atmpS1688 = _M0L1nS311 - _M0L6_2atmpS1689;
        int32_t _M0L5digitS313 = *(int32_t*)&_M0L6_2atmpS1688;
        int32_t _M0L6_2atmpS1686 = _M0L12digit__startS306 + _M0L6offsetS310;
        int32_t _M0L6_2atmpS1684 = _M0L6_2atmpS1686 - 1;
        int32_t _M0L6_2atmpS1685 =
          ((moonbit_string_t)moonbit_string_literal_74.data)[_M0L5digitS313];
        int32_t _M0L6_2atmpS1687;
        _M0L6bufferS305[_M0L6_2atmpS1684] = _M0L6_2atmpS1685;
        _M0L6_2atmpS1687 = _M0L6offsetS310 - 1;
        _M0L6offsetS310 = _M0L6_2atmpS1687;
        _M0L1nS311 = _M0L1qS312;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB20int__to__string__hex(
  uint16_t* _M0L6bufferS292,
  uint32_t _M0L3numS297,
  int32_t _M0L12digit__startS293,
  int32_t _M0L10total__lenS296
) {
  int32_t _M0L6_2atmpS1674;
  int32_t _M0L6offsetS287;
  uint32_t _M0L1nS288;
  #line 29 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS1674 = _M0L10total__lenS296 - _M0L12digit__startS293;
  _M0L6offsetS287 = _M0L6_2atmpS1674;
  _M0L1nS288 = _M0L3numS297;
  while (1) {
    if (_M0L6offsetS287 >= 2) {
      uint32_t _M0L6_2atmpS1671 = _M0L1nS288 & 255u;
      int32_t _M0L9byte__valS289 = *(int32_t*)&_M0L6_2atmpS1671;
      int32_t _M0L2hiS290 = _M0L9byte__valS289 / 16;
      int32_t _M0L2loS291 = _M0L9byte__valS289 % 16;
      int32_t _M0L6_2atmpS1665 = _M0L12digit__startS293 + _M0L6offsetS287;
      int32_t _M0L6_2atmpS1663 = _M0L6_2atmpS1665 - 2;
      int32_t _M0L6_2atmpS1664 =
        ((moonbit_string_t)moonbit_string_literal_74.data)[_M0L2hiS290];
      int32_t _M0L6_2atmpS1668;
      int32_t _M0L6_2atmpS1666;
      int32_t _M0L6_2atmpS1667;
      int32_t _M0L6_2atmpS1669;
      uint32_t _M0L6_2atmpS1670;
      _M0L6bufferS292[_M0L6_2atmpS1663] = _M0L6_2atmpS1664;
      _M0L6_2atmpS1668 = _M0L12digit__startS293 + _M0L6offsetS287;
      _M0L6_2atmpS1666 = _M0L6_2atmpS1668 - 1;
      _M0L6_2atmpS1667
      = ((moonbit_string_t)moonbit_string_literal_74.data)[
        _M0L2loS291
      ];
      _M0L6bufferS292[_M0L6_2atmpS1666] = _M0L6_2atmpS1667;
      _M0L6_2atmpS1669 = _M0L6offsetS287 - 2;
      _M0L6_2atmpS1670 = _M0L1nS288 >> 8;
      _M0L6offsetS287 = _M0L6_2atmpS1669;
      _M0L1nS288 = _M0L6_2atmpS1670;
      continue;
    } else if (_M0L6offsetS287 == 1) {
      uint32_t _M0L6_2atmpS1673 = _M0L1nS288 & 15u;
      int32_t _M0L6nibbleS295 = *(int32_t*)&_M0L6_2atmpS1673;
      int32_t _M0L6_2atmpS1672 =
        ((moonbit_string_t)moonbit_string_literal_74.data)[_M0L6nibbleS295];
      _M0L6bufferS292[_M0L12digit__startS293] = _M0L6_2atmpS1672;
    }
    break;
  }
  return 0;
}

int32_t _M0MPB4Iter4nextGcE(struct _M0TPB4IterGcE* _M0L4selfS270) {
  struct _M0TWEOc* _M0L7_2afuncS269;
  int32_t _M0L6resultS271;
  int64_t _M0L7_2abindS272;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS269 = _M0L4selfS270->$0;
  moonbit_incref(_M0L7_2afuncS269);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS271 = _M0L7_2afuncS269->code(_M0L7_2afuncS269);
  moonbit_decref(_M0L7_2afuncS269);
  _M0L7_2abindS272 = _M0L4selfS270->$1;
  if (_M0L6resultS271 == -1) {
    _M0L4selfS270->$1 = _M0MPB4Iter4nextN6constrS9981GcE;
  } else if (_M0L7_2abindS272 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS273 = _M0L7_2abindS272;
    int32_t _M0L4_2anS274 = (int32_t)_M0L7_2aSomeS273;
    int64_t _M0L6_2atmpS1657;
    if (_M0L4_2anS274 > 0) {
      int32_t _M0L6_2atmpS1658 = _M0L4_2anS274 - 1;
      _M0L6_2atmpS1657 = (int64_t)_M0L6_2atmpS1658;
    } else {
      _M0L6_2atmpS1657 = _M0MPB4Iter4nextN6constrS9980GcE;
    }
    _M0L4selfS270->$1 = _M0L6_2atmpS1657;
  }
  return _M0L6resultS271;
}

void* _M0MPB4Iter4nextGRPC16string10StringViewE(
  struct _M0TPB4IterGRPC16string10StringViewE* _M0L4selfS276
) {
  struct _M0TWERPC16option6OptionGRPC16string10StringViewE* _M0L7_2afuncS275;
  void* _M0L6resultS277;
  int64_t _M0L7_2abindS278;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS275 = _M0L4selfS276->$0;
  moonbit_incref(_M0L7_2afuncS275);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS277 = _M0L7_2afuncS275->code(_M0L7_2afuncS275);
  moonbit_decref(_M0L7_2afuncS275);
  _M0L7_2abindS278 = _M0L4selfS276->$1;
  switch (Moonbit_object_tag(_M0L6resultS277)) {
    case 1: {
      if (_M0L7_2abindS278 == 4294967296ll) {
        
      } else {
        int64_t _M0L7_2aSomeS279 = _M0L7_2abindS278;
        int32_t _M0L4_2anS280 = (int32_t)_M0L7_2aSomeS279;
        int64_t _M0L6_2atmpS1659;
        if (_M0L4_2anS280 > 0) {
          int32_t _M0L6_2atmpS1660 = _M0L4_2anS280 - 1;
          _M0L6_2atmpS1659 = (int64_t)_M0L6_2atmpS1660;
        } else {
          _M0L6_2atmpS1659
          = _M0MPB4Iter4nextN6constrS9980GRPC16string10StringViewE;
        }
        _M0L4selfS276->$1 = _M0L6_2atmpS1659;
      }
      break;
    }
    default: {
      _M0L4selfS276->$1
      = _M0MPB4Iter4nextN6constrS9981GRPC16string10StringViewE;
      break;
    }
  }
  return _M0L6resultS277;
}

struct _M0TUsRPB13StringBuilderE* _M0MPB4Iter4nextGUsRPB13StringBuilderEE(
  struct _M0TPB4IterGUsRPB13StringBuilderEE* _M0L4selfS282
) {
  struct _M0TWEOUsRPB13StringBuilderE* _M0L7_2afuncS281;
  struct _M0TUsRPB13StringBuilderE* _M0L6resultS283;
  int64_t _M0L7_2abindS284;
  #line 38 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L7_2afuncS281 = _M0L4selfS282->$0;
  moonbit_incref(_M0L7_2afuncS281);
  #line 41 "/home/developer/.moon/lib/core/builtin/iterator.mbt"
  _M0L6resultS283 = _M0L7_2afuncS281->code(_M0L7_2afuncS281);
  moonbit_decref(_M0L7_2afuncS281);
  _M0L7_2abindS284 = _M0L4selfS282->$1;
  if (_M0L6resultS283 == 0) {
    _M0L4selfS282->$1 = _M0MPB4Iter4nextN6constrS9981GUsRPB13StringBuilderEE;
  } else if (_M0L7_2abindS284 == 4294967296ll) {
    
  } else {
    int64_t _M0L7_2aSomeS285 = _M0L7_2abindS284;
    int32_t _M0L4_2anS286 = (int32_t)_M0L7_2aSomeS285;
    int64_t _M0L6_2atmpS1661;
    if (_M0L4_2anS286 > 0) {
      int32_t _M0L6_2atmpS1662 = _M0L4_2anS286 - 1;
      _M0L6_2atmpS1661 = (int64_t)_M0L6_2atmpS1662;
    } else {
      _M0L6_2atmpS1661 = _M0MPB4Iter4nextN6constrS9980GUsRPB13StringBuilderEE;
    }
    _M0L4selfS282->$1 = _M0L6_2atmpS1661;
  }
  return _M0L6resultS283;
}

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void* _M0L4selfS268
) {
  struct _M0TPB13StringBuilder* _M0L6loggerS267;
  struct _M0TPB6Logger _M0L6_2atmpS1656;
  moonbit_string_t _result_3517;
  #line 165 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 166 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS267 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  moonbit_incref(_M0L6loggerS267);
  _M0L6_2atmpS1656
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L6loggerS267
  };
  #line 167 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB7FailurePB4Show6output(_M0L4selfS268, _M0L6_2atmpS1656);
  if (_M0L6_2atmpS1656.$1) {
    moonbit_decref(_M0L6_2atmpS1656.$1);
  }
  #line 168 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _result_3517 = _M0MPB13StringBuilder10to__string(_M0L6loggerS267);
  moonbit_decref(_M0L6loggerS267);
  return _result_3517;
}

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t _M0L4selfS264,
  struct _M0TPB6Logger _M0L6loggerS263
) {
  moonbit_string_t _M0L6_2atmpS1654;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1654 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS264);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS263.$0->$method_0(_M0L6loggerS263.$1, _M0L6_2atmpS1654);
  moonbit_decref(_M0L6_2atmpS1654);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t _M0L4selfS266,
  struct _M0TPB6Logger _M0L6loggerS265
) {
  moonbit_string_t _M0L6_2atmpS1655;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1655 = _M0IPC13int3IntPB4Show10to__string(_M0L4selfS266);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS265.$0->$method_0(_M0L6loggerS265.$1, _M0L6_2atmpS1655);
  moonbit_decref(_M0L6_2atmpS1655);
  return 0;
}

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView _M0L4selfS262
) {
  #line 99 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  return _M0L4selfS262.$1;
}

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView _M0L4selfS261
) {
  moonbit_string_t _M0L8_2afieldS3176;
  #line 92 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L8_2afieldS3176 = _M0L4selfS261.$0;
  moonbit_incref(_M0L8_2afieldS3176);
  return _M0L8_2afieldS3176;
}

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS257,
  moonbit_string_t _M0L5valueS258,
  int32_t _M0L5startS259,
  int32_t _M0L3lenS260
) {
  int32_t _M0L6_2atmpS1653;
  int64_t _M0L6_2atmpS1652;
  struct _M0TPC16string10StringView _M0L6_2atmpS1651;
  #line 122 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1653 = _M0L5startS259 + _M0L3lenS260;
  _M0L6_2atmpS1652 = (int64_t)_M0L6_2atmpS1653;
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS1651
  = _M0MPC16string6String11sub_2einner(_M0L5valueS258, _M0L5startS259, _M0L6_2atmpS1652);
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L4selfS257, _M0L6_2atmpS1651);
  moonbit_decref(_M0L6_2atmpS1651.$0);
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t _M0L4selfS250,
  int32_t _M0L5startS256,
  int64_t _M0L3endS252
) {
  int32_t _M0L3lenS249;
  int32_t _M0L3endS251;
  int32_t _M0L5startS255;
  int32_t _if__result_3518;
  #line 755 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3lenS249 = Moonbit_array_length(_M0L4selfS250);
  if (_M0L3endS252 == 4294967296ll) {
    _M0L3endS251 = _M0L3lenS249;
  } else {
    int64_t _M0L7_2aSomeS253 = _M0L3endS252;
    int32_t _M0L6_2aendS254 = (int32_t)_M0L7_2aSomeS253;
    if (_M0L6_2aendS254 < 0) {
      _M0L3endS251 = _M0L3lenS249 + _M0L6_2aendS254;
    } else {
      _M0L3endS251 = _M0L6_2aendS254;
    }
  }
  if (_M0L5startS256 < 0) {
    _M0L5startS255 = _M0L3lenS249 + _M0L5startS256;
  } else {
    _M0L5startS255 = _M0L5startS256;
  }
  if (_M0L5startS255 >= 0) {
    if (_M0L5startS255 <= _M0L3endS251) {
      _if__result_3518 = _M0L3endS251 <= _M0L3lenS249;
    } else {
      _if__result_3518 = 0;
    }
  } else {
    _if__result_3518 = 0;
  }
  if (_if__result_3518) {
    if (_M0L5startS255 < _M0L3lenS249) {
      int32_t _M0L6_2atmpS1648 = _M0L4selfS250[_M0L5startS255];
      int32_t _M0L6_2atmpS1647;
      #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1647
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1648);
      if (!_M0L6_2atmpS1647) {
        
      } else {
        #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L3endS251 < _M0L3lenS249) {
      int32_t _M0L6_2atmpS1650 = _M0L4selfS250[_M0L3endS251];
      int32_t _M0L6_2atmpS1649;
      #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1649
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1650);
      if (!_M0L6_2atmpS1649) {
        
      } else {
        #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    moonbit_incref(_M0L4selfS250);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS250,
                                                 .$1 = _M0L5startS255,
                                                 .$2 = _M0L3endS251};
  } else {
    #line 763 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS248,
  struct _M0TPB4Show _M0L4showS247
) {
  struct _M0TPB6Logger _M0L6_2atmpS1646;
  #line 116 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS248);
  _M0L6_2atmpS1646
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS248
  };
  #line 117 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS247.$0->$method_0(_M0L4showS247.$1, _M0L6_2atmpS1646);
  if (_M0L6_2atmpS1646.$1) {
    moonbit_decref(_M0L6_2atmpS1646.$1);
  }
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS246,
  struct _M0TPB4Show _M0L4showS245
) {
  struct _M0TPB6Logger _M0L6_2atmpS1645;
  #line 111 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS246);
  _M0L6_2atmpS1645
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS246
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS245.$0->$method_0(_M0L4showS245.$1, _M0L6_2atmpS1645);
  if (_M0L6_2atmpS1645.$1) {
    moonbit_decref(_M0L6_2atmpS1645.$1);
  }
  return 0;
}

int32_t _M0FPB13finalize__acc(uint32_t _M0L3accS244) {
  uint32_t _M0L6_2atmpS1644;
  #line 444 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  #line 445 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1644 = _M0FPB14avalanche__acc(_M0L3accS244);
  return *(int32_t*)&_M0L6_2atmpS1644;
}

uint32_t _M0FPB14avalanche__acc(uint32_t _M0L3accS243) {
  uint32_t _M0Lm3accS242;
  uint32_t _M0L6_2atmpS1633;
  uint32_t _M0L6_2atmpS1635;
  uint32_t _M0L6_2atmpS1634;
  uint32_t _M0L6_2atmpS1636;
  uint32_t _M0L6_2atmpS1637;
  uint32_t _M0L6_2atmpS1639;
  uint32_t _M0L6_2atmpS1638;
  uint32_t _M0L6_2atmpS1640;
  uint32_t _M0L6_2atmpS1641;
  uint32_t _M0L6_2atmpS1643;
  uint32_t _M0L6_2atmpS1642;
  #line 449 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0Lm3accS242 = _M0L3accS243;
  _M0L6_2atmpS1633 = _M0Lm3accS242;
  _M0L6_2atmpS1635 = _M0Lm3accS242;
  _M0L6_2atmpS1634 = _M0L6_2atmpS1635 >> 15;
  _M0Lm3accS242 = _M0L6_2atmpS1633 ^ _M0L6_2atmpS1634;
  _M0L6_2atmpS1636 = _M0Lm3accS242;
  _M0Lm3accS242 = _M0L6_2atmpS1636 * 2246822519u;
  _M0L6_2atmpS1637 = _M0Lm3accS242;
  _M0L6_2atmpS1639 = _M0Lm3accS242;
  _M0L6_2atmpS1638 = _M0L6_2atmpS1639 >> 13;
  _M0Lm3accS242 = _M0L6_2atmpS1637 ^ _M0L6_2atmpS1638;
  _M0L6_2atmpS1640 = _M0Lm3accS242;
  _M0Lm3accS242 = _M0L6_2atmpS1640 * 3266489917u;
  _M0L6_2atmpS1641 = _M0Lm3accS242;
  _M0L6_2atmpS1643 = _M0Lm3accS242;
  _M0L6_2atmpS1642 = _M0L6_2atmpS1643 >> 16;
  _M0Lm3accS242 = _M0L6_2atmpS1641 ^ _M0L6_2atmpS1642;
  return _M0Lm3accS242;
}

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder* _M0L4selfS241,
  moonbit_string_t _M0L3strS240
) {
  int32_t _M0L8str__lenS239;
  int32_t _M0L3lenS1628;
  int32_t _M0L6_2atmpS1627;
  uint16_t* _M0L4dataS1629;
  int32_t _M0L3lenS1630;
  int32_t _M0L3lenS1632;
  int32_t _M0L6_2atmpS1631;
  #line 86 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L8str__lenS239 = Moonbit_array_length(_M0L3strS240);
  _M0L3lenS1628 = _M0L4selfS241->$1;
  _M0L6_2atmpS1627 = _M0L3lenS1628 + _M0L8str__lenS239;
  #line 88 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS241, _M0L6_2atmpS1627);
  _M0L4dataS1629 = _M0L4selfS241->$0;
  _M0L3lenS1630 = _M0L4selfS241->$1;
  moonbit_incref(_M0L4dataS1629);
  #line 89 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS1629, _M0L3lenS1630, _M0L3strS240, 0, _M0L8str__lenS239);
  moonbit_decref(_M0L4dataS1629);
  _M0L3lenS1632 = _M0L4selfS241->$1;
  _M0L6_2atmpS1631 = _M0L3lenS1632 + _M0L8str__lenS239;
  _M0L4selfS241->$1 = _M0L6_2atmpS1631;
  return 0;
}

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t* _M0L4selfS235,
  int32_t _M0L11dst__offsetS238,
  moonbit_string_t _M0L3strS236,
  int32_t _M0L11str__offsetS231,
  int32_t _M0L3lenS232
) {
  int32_t _M0L16end__str__offsetS230;
  int32_t _M0L1iS233;
  int32_t _M0L1jS234;
  #line 71 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L16end__str__offsetS230 = _M0L11str__offsetS231 + _M0L3lenS232;
  _M0L1iS233 = _M0L11str__offsetS231;
  _M0L1jS234 = _M0L11dst__offsetS238;
  while (1) {
    if (_M0L1iS233 < _M0L16end__str__offsetS230) {
      int32_t _M0L6_2atmpS1624 = _M0L3strS236[_M0L1iS233];
      int32_t _M0L6_2atmpS1625;
      int32_t _M0L6_2atmpS1626;
      if (
        _M0L1jS234 < 0 || _M0L1jS234 >= Moonbit_array_length(_M0L4selfS235)
      ) {
        #line 80 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
        moonbit_panic();
      }
      _M0L4selfS235[_M0L1jS234] = _M0L6_2atmpS1624;
      _M0L6_2atmpS1625 = _M0L1iS233 + 1;
      _M0L6_2atmpS1626 = _M0L1jS234 + 1;
      _M0L1iS233 = _M0L6_2atmpS1625;
      _M0L1jS234 = _M0L6_2atmpS1626;
      continue;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t _M0L4selfS228,
  int32_t _M0L5quoteS229
) {
  struct _M0TPB13StringBuilder* _M0L3bufS227;
  int32_t _M0L6_2atmpS1623;
  struct _M0TPC16string10StringView _M0L6_2atmpS1621;
  struct _M0TPB6Logger _M0L6_2atmpS1622;
  moonbit_string_t _result_3520;
  #line 110 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 111 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L3bufS227 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS1623 = Moonbit_array_length(_M0L4selfS228);
  moonbit_incref(_M0L4selfS228);
  _M0L6_2atmpS1621
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS228, .$1 = 0, .$2 = _M0L6_2atmpS1623
  };
  moonbit_incref(_M0L3bufS227);
  _M0L6_2atmpS1622
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS227
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0MPC16string10StringView18escape__to_2einner(_M0L6_2atmpS1621, _M0L6_2atmpS1622, _M0L5quoteS229);
  moonbit_decref(_M0L6_2atmpS1621.$0);
  if (_M0L6_2atmpS1622.$1) {
    moonbit_decref(_M0L6_2atmpS1622.$1);
  }
  #line 113 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_3520 = _M0MPB13StringBuilder10to__string(_M0L3bufS227);
  moonbit_decref(_M0L3bufS227);
  return _result_3520;
}

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView _M0L4selfS219,
  struct _M0TPB6Logger _M0L6loggerS217,
  int32_t _M0L5quoteS216
) {
  int32_t _M0L3endS1619;
  int32_t _M0L5startS1620;
  int32_t _M0L3lenS218;
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS220;
  int32_t _M0L1iS221;
  int32_t _M0L3segS222;
  #line 144 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L5quoteS216) {
    #line 150 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS217.$0->$method_3(_M0L6loggerS217.$1, 34);
  }
  _M0L3endS1619 = _M0L4selfS219.$2;
  _M0L5startS1620 = _M0L4selfS219.$1;
  _M0L3lenS218 = _M0L3endS1619 - _M0L5startS1620;
  moonbit_incref(_M0L4selfS219.$0);
  if (_M0L6loggerS217.$1) {
    moonbit_incref(_M0L6loggerS217.$1);
  }
  _M0L6_2aenvS220
  = (struct _M0TURPC16string10StringViewRPB6LoggerE*)moonbit_malloc(sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE));
  Moonbit_object_header(_M0L6_2aenvS220)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 136, 0);
  _M0L6_2aenvS220->$0 = _M0L4selfS219;
  _M0L6_2aenvS220->$1 = _M0L6loggerS217;
  _M0L1iS221 = 0;
  _M0L3segS222 = 0;
  _2afor_223:;
  while (1) {
    moonbit_string_t _M0L3strS1616;
    int32_t _M0L5startS1618;
    int32_t _M0L6_2atmpS1617;
    int32_t _M0L4codeS224;
    int32_t _M0L1cS226;
    int32_t _M0L6_2atmpS1600;
    int32_t _M0L6_2atmpS1601;
    int32_t _M0L6_2atmpS1602;
    int32_t _tmp_3524;
    int32_t _tmp_3525;
    if (_M0L1iS221 >= _M0L3lenS218) {
      #line 160 "/home/developer/.moon/lib/core/builtin/show.mbt"
      _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
      moonbit_decref(_M0L6_2aenvS220);
      break;
    }
    _M0L3strS1616 = _M0L4selfS219.$0;
    _M0L5startS1618 = _M0L4selfS219.$1;
    _M0L6_2atmpS1617 = _M0L5startS1618 + _M0L1iS221;
    _M0L4codeS224 = _M0L3strS1616[_M0L6_2atmpS1617];
    switch (_M0L4codeS224) {
      case 34: {
        _M0L1cS226 = _M0L4codeS224;
        goto join_225;
        break;
      }
      
      case 92: {
        _M0L1cS226 = _M0L4codeS224;
        goto join_225;
        break;
      }
      
      case 10: {
        int32_t _M0L6_2atmpS1603;
        int32_t _M0L6_2atmpS1604;
        #line 172 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
        #line 173 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS217.$0->$method_0(_M0L6loggerS217.$1, (moonbit_string_t)moonbit_string_literal_59.data);
        _M0L6_2atmpS1603 = _M0L1iS221 + 1;
        _M0L6_2atmpS1604 = _M0L1iS221 + 1;
        _M0L1iS221 = _M0L6_2atmpS1603;
        _M0L3segS222 = _M0L6_2atmpS1604;
        goto _2afor_223;
        break;
      }
      
      case 13: {
        int32_t _M0L6_2atmpS1605;
        int32_t _M0L6_2atmpS1606;
        #line 177 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
        #line 178 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS217.$0->$method_0(_M0L6loggerS217.$1, (moonbit_string_t)moonbit_string_literal_75.data);
        _M0L6_2atmpS1605 = _M0L1iS221 + 1;
        _M0L6_2atmpS1606 = _M0L1iS221 + 1;
        _M0L1iS221 = _M0L6_2atmpS1605;
        _M0L3segS222 = _M0L6_2atmpS1606;
        goto _2afor_223;
        break;
      }
      
      case 8: {
        int32_t _M0L6_2atmpS1607;
        int32_t _M0L6_2atmpS1608;
        #line 182 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
        #line 183 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS217.$0->$method_0(_M0L6loggerS217.$1, (moonbit_string_t)moonbit_string_literal_76.data);
        _M0L6_2atmpS1607 = _M0L1iS221 + 1;
        _M0L6_2atmpS1608 = _M0L1iS221 + 1;
        _M0L1iS221 = _M0L6_2atmpS1607;
        _M0L3segS222 = _M0L6_2atmpS1608;
        goto _2afor_223;
        break;
      }
      
      case 9: {
        int32_t _M0L6_2atmpS1609;
        int32_t _M0L6_2atmpS1610;
        #line 187 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
        #line 188 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS217.$0->$method_0(_M0L6loggerS217.$1, (moonbit_string_t)moonbit_string_literal_60.data);
        _M0L6_2atmpS1609 = _M0L1iS221 + 1;
        _M0L6_2atmpS1610 = _M0L1iS221 + 1;
        _M0L1iS221 = _M0L6_2atmpS1609;
        _M0L3segS222 = _M0L6_2atmpS1610;
        goto _2afor_223;
        break;
      }
      default: {
        if (_M0L4codeS224 < 32) {
          int32_t _M0L6_2atmpS1612;
          moonbit_string_t _M0L6_2atmpS1611;
          int32_t _M0L6_2atmpS1613;
          int32_t _M0L6_2atmpS1614;
          #line 193 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
          #line 194 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS217.$0->$method_0(_M0L6loggerS217.$1, (moonbit_string_t)moonbit_string_literal_77.data);
          _M0L6_2atmpS1612 = _M0L4codeS224 & 0xff;
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6_2atmpS1611 = _M0MPC14byte4Byte7to__hex(_M0L6_2atmpS1612);
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS217.$0->$method_0(_M0L6loggerS217.$1, _M0L6_2atmpS1611);
          moonbit_decref(_M0L6_2atmpS1611);
          #line 196 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS217.$0->$method_3(_M0L6loggerS217.$1, 125);
          _M0L6_2atmpS1613 = _M0L1iS221 + 1;
          _M0L6_2atmpS1614 = _M0L1iS221 + 1;
          _M0L1iS221 = _M0L6_2atmpS1613;
          _M0L3segS222 = _M0L6_2atmpS1614;
          goto _2afor_223;
        } else {
          int32_t _M0L6_2atmpS1615 = _M0L1iS221 + 1;
          int32_t _tmp_3523 = _M0L3segS222;
          _M0L1iS221 = _M0L6_2atmpS1615;
          _M0L3segS222 = _tmp_3523;
          goto _2afor_223;
        }
        break;
      }
    }
    goto joinlet_3522;
    join_225:;
    #line 166 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS220, _M0L3segS222, _M0L1iS221);
    #line 167 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS217.$0->$method_3(_M0L6loggerS217.$1, 92);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1600 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS226);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS217.$0->$method_3(_M0L6loggerS217.$1, _M0L6_2atmpS1600);
    _M0L6_2atmpS1601 = _M0L1iS221 + 1;
    _M0L6_2atmpS1602 = _M0L1iS221 + 1;
    _M0L1iS221 = _M0L6_2atmpS1601;
    _M0L3segS222 = _M0L6_2atmpS1602;
    continue;
    joinlet_3522:;
    _tmp_3524 = _M0L1iS221;
    _tmp_3525 = _M0L3segS222;
    _M0L1iS221 = _tmp_3524;
    _M0L3segS222 = _tmp_3525;
    continue;
    break;
  }
  if (_M0L5quoteS216) {
    #line 204 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS217.$0->$method_3(_M0L6loggerS217.$1, 34);
  }
  return 0;
}

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS212,
  int32_t _M0L3segS215,
  int32_t _M0L1iS214
) {
  struct _M0TPB6Logger _M0L6loggerS211;
  struct _M0TPC16string10StringView _M0L4selfS213;
  #line 153 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6loggerS211 = _M0L6_2aenvS212->$1;
  _M0L4selfS213 = _M0L6_2aenvS212->$0;
  if (_M0L1iS214 > _M0L3segS215) {
    int64_t _M0L6_2atmpS1599 = (int64_t)_M0L1iS214;
    struct _M0TPC16string10StringView _M0L6_2atmpS1598;
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1598
    = _M0MPC16string10StringView11sub_2einner(_M0L4selfS213, _M0L3segS215, _M0L6_2atmpS1599);
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS211.$0->$method_2(_M0L6loggerS211.$1, _M0L6_2atmpS1598);
    moonbit_decref(_M0L6_2atmpS1598.$0);
  }
  return 0;
}

int32_t _M0MPC16string10StringView11unsafe__get(
  struct _M0TPC16string10StringView _M0L4selfS209,
  int32_t _M0L5indexS210
) {
  moonbit_string_t _M0L3strS1595;
  int32_t _M0L5startS1597;
  int32_t _M0L6_2atmpS1596;
  #line 128 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1595 = _M0L4selfS209.$0;
  _M0L5startS1597 = _M0L4selfS209.$1;
  _M0L6_2atmpS1596 = _M0L5startS1597 + _M0L5indexS210;
  return _M0L3strS1595[_M0L6_2atmpS1596];
}

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView _M0L4selfS202,
  int32_t _M0L5startS208,
  int64_t _M0L3endS204
) {
  moonbit_string_t _M0L3strS1594;
  int32_t _M0L8str__lenS201;
  int32_t _M0L8abs__endS203;
  int32_t _M0L10abs__startS207;
  int32_t _M0L5startS1582;
  int32_t _if__result_3526;
  #line 814 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS1594 = _M0L4selfS202.$0;
  _M0L8str__lenS201 = Moonbit_array_length(_M0L3strS1594);
  if (_M0L3endS204 == 4294967296ll) {
    _M0L8abs__endS203 = _M0L4selfS202.$2;
  } else {
    int64_t _M0L7_2aSomeS205 = _M0L3endS204;
    int32_t _M0L6_2aendS206 = (int32_t)_M0L7_2aSomeS205;
    if (_M0L6_2aendS206 < 0) {
      int32_t _M0L3endS1592 = _M0L4selfS202.$2;
      _M0L8abs__endS203 = _M0L3endS1592 + _M0L6_2aendS206;
    } else {
      int32_t _M0L5startS1593 = _M0L4selfS202.$1;
      _M0L8abs__endS203 = _M0L5startS1593 + _M0L6_2aendS206;
    }
  }
  if (_M0L5startS208 < 0) {
    int32_t _M0L3endS1590 = _M0L4selfS202.$2;
    _M0L10abs__startS207 = _M0L3endS1590 + _M0L5startS208;
  } else {
    int32_t _M0L5startS1591 = _M0L4selfS202.$1;
    _M0L10abs__startS207 = _M0L5startS1591 + _M0L5startS208;
  }
  _M0L5startS1582 = _M0L4selfS202.$1;
  if (_M0L10abs__startS207 >= _M0L5startS1582) {
    if (_M0L10abs__startS207 <= _M0L8abs__endS203) {
      int32_t _M0L3endS1581 = _M0L4selfS202.$2;
      _if__result_3526 = _M0L8abs__endS203 <= _M0L3endS1581;
    } else {
      _if__result_3526 = 0;
    }
  } else {
    _if__result_3526 = 0;
  }
  if (_if__result_3526) {
    moonbit_string_t _M0L3strS1589;
    if (_M0L10abs__startS207 < _M0L8str__lenS201) {
      moonbit_string_t _M0L3strS1585 = _M0L4selfS202.$0;
      int32_t _M0L6_2atmpS1584 = _M0L3strS1585[_M0L10abs__startS207];
      int32_t _M0L6_2atmpS1583;
      #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1583
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1584);
      if (!_M0L6_2atmpS1583) {
        
      } else {
        #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L8abs__endS203 < _M0L8str__lenS201) {
      moonbit_string_t _M0L3strS1588 = _M0L4selfS202.$0;
      int32_t _M0L6_2atmpS1587 = _M0L3strS1588[_M0L8abs__endS203];
      int32_t _M0L6_2atmpS1586;
      #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS1586
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS1587);
      if (!_M0L6_2atmpS1586) {
        
      } else {
        #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    _M0L3strS1589 = _M0L4selfS202.$0;
    moonbit_incref(_M0L3strS1589);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS1589,
                                                 .$1 = _M0L10abs__startS207,
                                                 .$2 = _M0L8abs__endS203};
  } else {
    #line 834 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC16string10StringView6length(
  struct _M0TPC16string10StringView _M0L4selfS200
) {
  int32_t _M0L3endS1579;
  int32_t _M0L5startS1580;
  #line 49 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3endS1579 = _M0L4selfS200.$2;
  _M0L5startS1580 = _M0L4selfS200.$1;
  return _M0L3endS1579 - _M0L5startS1580;
}

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t _M0L1bS199) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS198;
  int32_t _M0L6_2atmpS1576;
  int32_t _M0L6_2atmpS1575;
  int32_t _M0L6_2atmpS1578;
  int32_t _M0L6_2atmpS1577;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS1574;
  moonbit_string_t _result_3527;
  #line 74 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L7_2aselfS198 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1576 = _M0IPC14byte4BytePB3Div3div(_M0L1bS199, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1575
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1576);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS198, _M0L6_2atmpS1575);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1578 = _M0IPC14byte4BytePB3Mod3mod(_M0L1bS199, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS1577
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS1578);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS198, _M0L6_2atmpS1577);
  _M0L6_2atmpS1574 = _M0L7_2aselfS198;
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_3527 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS1574);
  moonbit_decref(_M0L6_2atmpS1574);
  return _result_3527;
}

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t _M0L1iS197) {
  #line 75 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L1iS197 < 10) {
    int32_t _M0L6_2atmpS1571;
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1571 = _M0IPC14byte4BytePB3Add3add(_M0L1iS197, 48);
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1571);
  } else {
    int32_t _M0L6_2atmpS1573;
    int32_t _M0L6_2atmpS1572;
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1573 = _M0IPC14byte4BytePB3Add3add(_M0L1iS197, 97);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS1572 = _M0IPC14byte4BytePB3Sub3sub(_M0L6_2atmpS1573, 10);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS1572);
  }
}

int32_t _M0IPC14byte4BytePB3Sub3sub(
  int32_t _M0L4selfS195,
  int32_t _M0L4thatS196
) {
  int32_t _M0L6_2atmpS1569;
  int32_t _M0L6_2atmpS1570;
  int32_t _M0L6_2atmpS1568;
  #line 120 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1569 = (int32_t)_M0L4selfS195;
  _M0L6_2atmpS1570 = (int32_t)_M0L4thatS196;
  _M0L6_2atmpS1568 = _M0L6_2atmpS1569 - _M0L6_2atmpS1570;
  return _M0L6_2atmpS1568 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Mod3mod(
  int32_t _M0L4selfS193,
  int32_t _M0L4thatS194
) {
  int32_t _M0L6_2atmpS1566;
  int32_t _M0L6_2atmpS1567;
  int32_t _M0L6_2atmpS1565;
  #line 67 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1566 = (int32_t)_M0L4selfS193;
  _M0L6_2atmpS1567 = (int32_t)_M0L4thatS194;
  _M0L6_2atmpS1565 = _M0L6_2atmpS1566 % _M0L6_2atmpS1567;
  return _M0L6_2atmpS1565 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Div3div(
  int32_t _M0L4selfS191,
  int32_t _M0L4thatS192
) {
  int32_t _M0L6_2atmpS1563;
  int32_t _M0L6_2atmpS1564;
  int32_t _M0L6_2atmpS1562;
  #line 62 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1563 = (int32_t)_M0L4selfS191;
  _M0L6_2atmpS1564 = (int32_t)_M0L4thatS192;
  _M0L6_2atmpS1562 = _M0L6_2atmpS1563 / _M0L6_2atmpS1564;
  return _M0L6_2atmpS1562 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Add3add(
  int32_t _M0L4selfS189,
  int32_t _M0L4thatS190
) {
  int32_t _M0L6_2atmpS1560;
  int32_t _M0L6_2atmpS1561;
  int32_t _M0L6_2atmpS1559;
  #line 106 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS1560 = (int32_t)_M0L4selfS189;
  _M0L6_2atmpS1561 = (int32_t)_M0L4thatS190;
  _M0L6_2atmpS1559 = _M0L6_2atmpS1560 + _M0L6_2atmpS1561;
  return _M0L6_2atmpS1559 & 0xff;
}

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t _M0L4selfS188) {
  int32_t _M0L6_2atmpS1558;
  #line 68 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  _M0L6_2atmpS1558 = (int32_t)_M0L4selfS188;
  return _M0L6_2atmpS1558;
}

int32_t _M0FPB32code__point__of__surrogate__pair(
  int32_t _M0L7leadingS186,
  int32_t _M0L8trailingS187
) {
  int32_t _M0L6_2atmpS1557;
  int32_t _M0L6_2atmpS1556;
  int32_t _M0L6_2atmpS1555;
  int32_t _M0L6_2atmpS1554;
  int32_t _M0L6_2atmpS1553;
  #line 40 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS1557 = _M0L7leadingS186 - 55296;
  _M0L6_2atmpS1556 = _M0L6_2atmpS1557 * 1024;
  _M0L6_2atmpS1555 = _M0L6_2atmpS1556 + _M0L8trailingS187;
  _M0L6_2atmpS1554 = _M0L6_2atmpS1555 - 56320;
  _M0L6_2atmpS1553 = _M0L6_2atmpS1554 + 65536;
  return _M0L6_2atmpS1553;
}

int32_t _M0MPC16string6String20char__length_2einner(
  moonbit_string_t _M0L4selfS179,
  int32_t _M0L13start__offsetS180,
  int64_t _M0L11end__offsetS177
) {
  int32_t _M0L11end__offsetS176;
  int32_t _if__result_3528;
  #line 60 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L11end__offsetS177 == 4294967296ll) {
    _M0L11end__offsetS176 = Moonbit_array_length(_M0L4selfS179);
  } else {
    int64_t _M0L7_2aSomeS178 = _M0L11end__offsetS177;
    _M0L11end__offsetS176 = (int32_t)_M0L7_2aSomeS178;
  }
  if (_M0L13start__offsetS180 >= 0) {
    if (_M0L13start__offsetS180 <= _M0L11end__offsetS176) {
      int32_t _M0L6_2atmpS1546 = Moonbit_array_length(_M0L4selfS179);
      _if__result_3528 = _M0L11end__offsetS176 <= _M0L6_2atmpS1546;
    } else {
      _if__result_3528 = 0;
    }
  } else {
    _if__result_3528 = 0;
  }
  if (_if__result_3528) {
    int32_t _M0L12utf16__indexS181 = _M0L13start__offsetS180;
    int32_t _M0L11char__countS182 = 0;
    while (1) {
      if (_M0L12utf16__indexS181 < _M0L11end__offsetS176) {
        int32_t _M0L2c1S183 = _M0L4selfS179[_M0L12utf16__indexS181];
        int32_t _if__result_3530;
        int32_t _M0L6_2atmpS1551;
        int32_t _M0L6_2atmpS1552;
        #line 76 "/home/developer/.moon/lib/core/builtin/string.mbt"
        if (_M0MPC16uint166UInt1622is__leading__surrogate(_M0L2c1S183)) {
          int32_t _M0L6_2atmpS1547 = _M0L12utf16__indexS181 + 1;
          _if__result_3530 = _M0L6_2atmpS1547 < _M0L11end__offsetS176;
        } else {
          _if__result_3530 = 0;
        }
        if (_if__result_3530) {
          int32_t _M0L6_2atmpS1550 = _M0L12utf16__indexS181 + 1;
          int32_t _M0L2c2S184 = _M0L4selfS179[_M0L6_2atmpS1550];
          #line 78 "/home/developer/.moon/lib/core/builtin/string.mbt"
          if (_M0MPC16uint166UInt1623is__trailing__surrogate(_M0L2c2S184)) {
            int32_t _M0L6_2atmpS1548 = _M0L12utf16__indexS181 + 2;
            int32_t _M0L6_2atmpS1549 = _M0L11char__countS182 + 1;
            _M0L12utf16__indexS181 = _M0L6_2atmpS1548;
            _M0L11char__countS182 = _M0L6_2atmpS1549;
            continue;
          } else {
            #line 81 "/home/developer/.moon/lib/core/builtin/string.mbt"
            _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_78.data);
          }
        }
        _M0L6_2atmpS1551 = _M0L12utf16__indexS181 + 1;
        _M0L6_2atmpS1552 = _M0L11char__countS182 + 1;
        _M0L12utf16__indexS181 = _M0L6_2atmpS1551;
        _M0L11char__countS182 = _M0L6_2atmpS1552;
        continue;
      } else {
        return _M0L11char__countS182;
      }
      break;
    }
  } else {
    #line 70 "/home/developer/.moon/lib/core/builtin/string.mbt"
    return _M0FPC15abort5abortGiE((moonbit_string_t)moonbit_string_literal_79.data);
  }
}

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t _M0L4selfS175) {
  #line 45 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS175 >= 56320) {
    return _M0L4selfS175 <= 57343;
  } else {
    return 0;
  }
}

int32_t _M0MPC16uint166UInt1622is__leading__surrogate(int32_t _M0L4selfS174) {
  #line 28 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS174 >= 55296) {
    return _M0L4selfS174 <= 56319;
  } else {
    return 0;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder* _M0L4selfS172,
  int32_t _M0L2chS171
) {
  uint32_t _M0L4codeS170;
  #line 95 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  #line 96 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4codeS170 = _M0MPC14char4Char8to__uint(_M0L2chS171);
  if (_M0L4codeS170 <= 65535u) {
    int32_t _M0L3lenS1525 = _M0L4selfS172->$1;
    int32_t _M0L6_2atmpS1524 = _M0L3lenS1525 + 1;
    uint16_t* _M0L4dataS1526;
    int32_t _M0L3lenS1527;
    int32_t _M0L6_2atmpS1528;
    int32_t _M0L3lenS1530;
    int32_t _M0L6_2atmpS1529;
    #line 98 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS172, _M0L6_2atmpS1524);
    _M0L4dataS1526 = _M0L4selfS172->$0;
    _M0L3lenS1527 = _M0L4selfS172->$1;
    moonbit_incref(_M0L4dataS1526);
    #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1528 = _M0MPC14uint4UInt10to__uint16(_M0L4codeS170);
    if (
      _M0L3lenS1527 < 0
      || _M0L3lenS1527 >= Moonbit_array_length(_M0L4dataS1526)
    ) {
      #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1526[_M0L3lenS1527] = _M0L6_2atmpS1528;
    moonbit_decref(_M0L4dataS1526);
    _M0L3lenS1530 = _M0L4selfS172->$1;
    _M0L6_2atmpS1529 = _M0L3lenS1530 + 1;
    _M0L4selfS172->$1 = _M0L6_2atmpS1529;
  } else if (_M0L4codeS170 <= 1114111u) {
    int32_t _M0L3lenS1532 = _M0L4selfS172->$1;
    int32_t _M0L6_2atmpS1531 = _M0L3lenS1532 + 2;
    uint32_t _M0L4codeS173;
    uint16_t* _M0L4dataS1533;
    int32_t _M0L3lenS1534;
    uint32_t _M0L6_2atmpS1537;
    uint32_t _M0L6_2atmpS1536;
    int32_t _M0L6_2atmpS1535;
    uint16_t* _M0L4dataS1538;
    int32_t _M0L3lenS1543;
    int32_t _M0L6_2atmpS1539;
    uint32_t _M0L6_2atmpS1542;
    uint32_t _M0L6_2atmpS1541;
    int32_t _M0L6_2atmpS1540;
    int32_t _M0L3lenS1545;
    int32_t _M0L6_2atmpS1544;
    #line 102 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS172, _M0L6_2atmpS1531);
    _M0L4codeS173 = _M0L4codeS170 - 65536u;
    _M0L4dataS1533 = _M0L4selfS172->$0;
    _M0L3lenS1534 = _M0L4selfS172->$1;
    _M0L6_2atmpS1537 = _M0L4codeS173 >> 10;
    _M0L6_2atmpS1536 = 55296u + _M0L6_2atmpS1537;
    moonbit_incref(_M0L4dataS1533);
    #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1535 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1536);
    if (
      _M0L3lenS1534 < 0
      || _M0L3lenS1534 >= Moonbit_array_length(_M0L4dataS1533)
    ) {
      #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1533[_M0L3lenS1534] = _M0L6_2atmpS1535;
    moonbit_decref(_M0L4dataS1533);
    _M0L4dataS1538 = _M0L4selfS172->$0;
    _M0L3lenS1543 = _M0L4selfS172->$1;
    _M0L6_2atmpS1539 = _M0L3lenS1543 + 1;
    _M0L6_2atmpS1542 = _M0L4codeS173 & 1023u;
    _M0L6_2atmpS1541 = 56320u + _M0L6_2atmpS1542;
    moonbit_incref(_M0L4dataS1538);
    #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1540 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS1541);
    if (
      _M0L6_2atmpS1539 < 0
      || _M0L6_2atmpS1539 >= Moonbit_array_length(_M0L4dataS1538)
    ) {
      #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS1538[_M0L6_2atmpS1539] = _M0L6_2atmpS1540;
    moonbit_decref(_M0L4dataS1538);
    _M0L3lenS1545 = _M0L4selfS172->$1;
    _M0L6_2atmpS1544 = _M0L3lenS1545 + 2;
    _M0L4selfS172->$1 = _M0L6_2atmpS1544;
  } else {
    #line 108 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_80.data);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder* _M0L4selfS164,
  int32_t _M0L8requiredS165
) {
  uint16_t* _M0L4dataS1523;
  int32_t _M0L12current__lenS163;
  int32_t _M0L13enough__spaceS166;
  int32_t _M0L13enough__spaceS167;
  uint16_t* _M0L4dataS1519;
  int32_t _M0L6_2atmpS1520;
  int32_t _M0L3lenS1521;
  uint16_t* _M0L9new__dataS169;
  uint16_t* _M0L6_2aoldS3187;
  #line 46 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4dataS1523 = _M0L4selfS164->$0;
  _M0L12current__lenS163 = Moonbit_array_length(_M0L4dataS1523);
  if (_M0L8requiredS165 <= _M0L12current__lenS163) {
    return 0;
  }
  _M0L13enough__spaceS167 = _M0L12current__lenS163;
  while (1) {
    if (_M0L13enough__spaceS167 < _M0L8requiredS165) {
      int32_t _M0L6_2atmpS1522 = _M0L13enough__spaceS167 * 2;
      _M0L13enough__spaceS167 = _M0L6_2atmpS1522;
      continue;
    } else {
      _M0L13enough__spaceS166 = _M0L13enough__spaceS167;
    }
    break;
  }
  _M0L4dataS1519 = _M0L4selfS164->$0;
  moonbit_incref(_M0L4dataS1519);
  #line 64 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS1520 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L3lenS1521 = _M0L4selfS164->$1;
  #line 61 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L9new__dataS169
  = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1519, _M0L13enough__spaceS166, _M0L6_2atmpS1520, _M0L3lenS1521, 0, 0);
  moonbit_decref(_M0L4dataS1519);
  _M0L6_2aoldS3187 = _M0L4selfS164->$0;
  moonbit_decref(_M0L6_2aoldS3187);
  _M0L4selfS164->$0 = _M0L9new__dataS169;
  return 0;
}

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t _M0L4selfS162) {
  int32_t _M0L6_2atmpS1518;
  #line 2676 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1518 = *(int32_t*)&_M0L4selfS162;
  return (uint16_t)_M0L6_2atmpS1518;
}

uint32_t _M0MPC14char4Char8to__uint(int32_t _M0L4selfS161) {
  int32_t _M0L6_2atmpS1517;
  #line 1254 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1517 = _M0L4selfS161;
  return *(uint32_t*)&_M0L6_2atmpS1517;
}

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder* _M0L4selfS159
) {
  int32_t _M0L3lenS1509;
  uint16_t* _M0L4dataS1511;
  int32_t _M0L6_2atmpS1510;
  #line 148 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3lenS1509 = _M0L4selfS159->$1;
  _M0L4dataS1511 = _M0L4selfS159->$0;
  _M0L6_2atmpS1510 = Moonbit_array_length(_M0L4dataS1511);
  if (_M0L3lenS1509 == _M0L6_2atmpS1510) {
    uint16_t* _M0L4dataS1512 = _M0L4selfS159->$0;
    moonbit_incref(_M0L4dataS1512);
    return _M0L4dataS1512;
  } else {
    uint16_t* _M0L4dataS1513 = _M0L4selfS159->$0;
    int32_t _M0L3lenS1514 = _M0L4selfS159->$1;
    int32_t _M0L6_2atmpS1515;
    int32_t _M0L3lenS1516;
    uint16_t* _M0L4dataS160;
    moonbit_incref(_M0L4dataS1513);
    #line 155 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS1515 = _M0IPC16uint166UInt16PB7Default7default();
    _M0L3lenS1516 = _M0L4selfS159->$1;
    #line 152 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L4dataS160
    = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS1513, _M0L3lenS1514, _M0L6_2atmpS1515, _M0L3lenS1516, 0, 0);
    moonbit_decref(_M0L4dataS1513);
    return _M0L4dataS160;
  }
}

int32_t _M0IPC16uint166UInt16PB7Default7default() {
  #line 176 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  return 0;
}

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t* _M0L3srcS156,
  int32_t _M0L13allocate__lenS152,
  int32_t _M0L4initS157,
  int32_t _M0L3lenS153,
  int32_t _M0L11src__offsetS154,
  int32_t _M0L11dst__offsetS155
) {
  int32_t _if__result_3532;
  #line 97 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L13allocate__lenS152 >= 0) {
    if (_M0L3lenS153 >= 0) {
      if (_M0L11src__offsetS154 >= 0) {
        if (_M0L11dst__offsetS155 >= 0) {
          int32_t _M0L6_2atmpS1505 = _M0L11src__offsetS154 + _M0L3lenS153;
          int32_t _M0L6_2atmpS1506 = Moonbit_array_length(_M0L3srcS156);
          if (_M0L6_2atmpS1505 <= _M0L6_2atmpS1506) {
            int32_t _M0L6_2atmpS1504 = _M0L11dst__offsetS155 + _M0L3lenS153;
            _if__result_3532 = _M0L6_2atmpS1504 <= _M0L13allocate__lenS152;
          } else {
            _if__result_3532 = 0;
          }
        } else {
          _if__result_3532 = 0;
        }
      } else {
        _if__result_3532 = 0;
      }
    } else {
      _if__result_3532 = 0;
    }
  } else {
    _if__result_3532 = 0;
  }
  if (_if__result_3532) {
    moonbit_incref(_M0L3srcS156);
    #line 115 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    return _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(_M0L3srcS156, _M0L13allocate__lenS152, _M0L4initS157, _M0L11src__offsetS154, _M0L11dst__offsetS155, _M0L3lenS153);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS158;
    int32_t _M0L6_2atmpS1508;
    moonbit_string_t _M0L6_2atmpS1507;
    uint16_t* _result_3533;
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L18_2astring__builderS158
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS158, (moonbit_string_t)moonbit_string_literal_81.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS158, _M0L13allocate__lenS152);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS158, (moonbit_string_t)moonbit_string_literal_82.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS158, _M0L11src__offsetS154);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS158, (moonbit_string_t)moonbit_string_literal_83.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS158, _M0L11dst__offsetS155);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS158, (moonbit_string_t)moonbit_string_literal_84.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS158, _M0L3lenS153);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS158, (moonbit_string_t)moonbit_string_literal_85.data);
    _M0L6_2atmpS1508 = Moonbit_array_length(_M0L3srcS156);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS158, _M0L6_2atmpS1508);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L6_2atmpS1507
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS158);
    moonbit_decref(_M0L18_2astring__builderS158);
    #line 111 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _result_3533 = _M0FPC15abort5abortGAkE(_M0L6_2atmpS1507);
    moonbit_decref(_M0L6_2atmpS1507);
    return _result_3533;
  }
}

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t* _M0L3srcS149,
  int32_t _M0L13allocate__lenS146,
  int32_t _M0L4initS147,
  int32_t _M0L11src__offsetS150,
  int32_t _M0L11dst__offsetS148,
  int32_t _M0L9blit__lenS151
) {
  uint16_t* _M0L3dstS145;
  #line 79 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  _M0L3dstS145
  = (uint16_t*)moonbit_make_string(_M0L13allocate__lenS146, _M0L4initS147);
  moonbit_incref(_M0L3dstS145);
  #line 90 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS145, _M0L11dst__offsetS148, _M0L3srcS149, _M0L11src__offsetS150, _M0L9blit__lenS151, sizeof(uint16_t));
  return _M0L3dstS145;
}

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t _M0L10size__hintS143
) {
  int32_t _M0L7initialS142;
  uint16_t* _M0L4dataS144;
  struct _M0TPB13StringBuilder* _block_3534;
  #line 32 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  if (_M0L10size__hintS143 < 1) {
    _M0L7initialS142 = 1;
  } else {
    int32_t _M0L6_2atmpS1503 = _M0L10size__hintS143 + 1;
    _M0L7initialS142 = _M0L6_2atmpS1503 / 2;
  }
  _M0L4dataS144 = (uint16_t*)moonbit_make_string(_M0L7initialS142, 0);
  _block_3534
  = (struct _M0TPB13StringBuilder*)moonbit_malloc(sizeof(struct _M0TPB13StringBuilder));
  Moonbit_object_header(_block_3534)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 141, 0);
  _block_3534->$0 = _M0L4dataS144;
  _block_3534->$1 = 0;
  return _block_3534;
}

int32_t _M0MPC14byte4Byte8to__char(int32_t _M0L4selfS141) {
  int32_t _M0L6_2atmpS1502;
  #line 1868 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS1502 = (int32_t)_M0L4selfS141;
  return _M0L6_2atmpS1502;
}

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t* _M0L3srcS121,
  int32_t _M0L13allocate__lenS117,
  int32_t _M0L3lenS118,
  int32_t _M0L11src__offsetS119,
  int32_t _M0L11dst__offsetS120
) {
  int32_t _if__result_3535;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS117 >= 0) {
    if (_M0L3lenS118 >= 0) {
      if (_M0L11src__offsetS119 >= 0) {
        if (_M0L11dst__offsetS120 >= 0) {
          int32_t _M0L6_2atmpS1483 = _M0L11src__offsetS119 + _M0L3lenS118;
          int32_t _M0L6_2atmpS1484;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1484
          = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS121);
          if (_M0L6_2atmpS1483 <= _M0L6_2atmpS1484) {
            int32_t _M0L6_2atmpS1482 = _M0L11dst__offsetS120 + _M0L3lenS118;
            _if__result_3535 = _M0L6_2atmpS1482 <= _M0L13allocate__lenS117;
          } else {
            _if__result_3535 = 0;
          }
        } else {
          _if__result_3535 = 0;
        }
      } else {
        _if__result_3535 = 0;
      }
    } else {
      _if__result_3535 = 0;
    }
  } else {
    _if__result_3535 = 0;
  }
  if (_if__result_3535) {
    moonbit_incref(_M0L3srcS121);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (moonbit_string_t*)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS117, (moonbit_string_t)moonbit_string_literal_0.data, _M0L3srcS121, _M0L11src__offsetS119, _M0L11dst__offsetS120, _M0L3lenS118);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS122;
    int32_t _M0L6_2atmpS1486;
    moonbit_string_t _M0L6_2atmpS1485;
    moonbit_string_t* _result_3536;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS122
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS122, (moonbit_string_t)moonbit_string_literal_81.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS122, _M0L13allocate__lenS117);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS122, (moonbit_string_t)moonbit_string_literal_82.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS122, _M0L11src__offsetS119);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS122, (moonbit_string_t)moonbit_string_literal_83.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS122, _M0L11dst__offsetS120);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS122, (moonbit_string_t)moonbit_string_literal_84.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS122, _M0L3lenS118);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS122, (moonbit_string_t)moonbit_string_literal_85.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1486 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS121);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS122, _M0L6_2atmpS1486);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1485
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS122);
    moonbit_decref(_M0L18_2astring__builderS122);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3536
    = _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(_M0L6_2atmpS1485);
    moonbit_decref(_M0L6_2atmpS1485);
    return _result_3536;
  }
}

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE** _M0L3srcS127,
  int32_t _M0L13allocate__lenS123,
  int32_t _M0L3lenS124,
  int32_t _M0L11src__offsetS125,
  int32_t _M0L11dst__offsetS126
) {
  int32_t _if__result_3537;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS123 >= 0) {
    if (_M0L3lenS124 >= 0) {
      if (_M0L11src__offsetS125 >= 0) {
        if (_M0L11dst__offsetS126 >= 0) {
          int32_t _M0L6_2atmpS1488 = _M0L11src__offsetS125 + _M0L3lenS124;
          int32_t _M0L6_2atmpS1489;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1489
          = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS127);
          if (_M0L6_2atmpS1488 <= _M0L6_2atmpS1489) {
            int32_t _M0L6_2atmpS1487 = _M0L11dst__offsetS126 + _M0L3lenS124;
            _if__result_3537 = _M0L6_2atmpS1487 <= _M0L13allocate__lenS123;
          } else {
            _if__result_3537 = 0;
          }
        } else {
          _if__result_3537 = 0;
        }
      } else {
        _if__result_3537 = 0;
      }
    } else {
      _if__result_3537 = 0;
    }
  } else {
    _if__result_3537 = 0;
  }
  if (_if__result_3537) {
    moonbit_incref(_M0L3srcS127);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TUsiE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS123, 0, _M0L3srcS127, _M0L11src__offsetS125, _M0L11dst__offsetS126, _M0L3lenS124);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS128;
    int32_t _M0L6_2atmpS1491;
    moonbit_string_t _M0L6_2atmpS1490;
    struct _M0TUsiE** _result_3538;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS128
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS128, (moonbit_string_t)moonbit_string_literal_81.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS128, _M0L13allocate__lenS123);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS128, (moonbit_string_t)moonbit_string_literal_82.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS128, _M0L11src__offsetS125);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS128, (moonbit_string_t)moonbit_string_literal_83.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS128, _M0L11dst__offsetS126);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS128, (moonbit_string_t)moonbit_string_literal_84.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS128, _M0L3lenS124);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS128, (moonbit_string_t)moonbit_string_literal_85.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1491 = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS127);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS128, _M0L6_2atmpS1491);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1490
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS128);
    moonbit_decref(_M0L18_2astring__builderS128);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3538
    = _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(_M0L6_2atmpS1490);
    moonbit_decref(_M0L6_2atmpS1490);
    return _result_3538;
  }
}

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23make__and__blit_2einnerGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS133,
  int32_t _M0L13allocate__lenS129,
  int32_t _M0L3lenS130,
  int32_t _M0L11src__offsetS131,
  int32_t _M0L11dst__offsetS132
) {
  int32_t _if__result_3539;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS129 >= 0) {
    if (_M0L3lenS130 >= 0) {
      if (_M0L11src__offsetS131 >= 0) {
        if (_M0L11dst__offsetS132 >= 0) {
          int32_t _M0L6_2atmpS1493 = _M0L11src__offsetS131 + _M0L3lenS130;
          int32_t _M0L6_2atmpS1494;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1494
          = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L3srcS133);
          if (_M0L6_2atmpS1493 <= _M0L6_2atmpS1494) {
            int32_t _M0L6_2atmpS1492 = _M0L11dst__offsetS132 + _M0L3lenS130;
            _if__result_3539 = _M0L6_2atmpS1492 <= _M0L13allocate__lenS129;
          } else {
            _if__result_3539 = 0;
          }
        } else {
          _if__result_3539 = 0;
        }
      } else {
        _if__result_3539 = 0;
      }
    } else {
      _if__result_3539 = 0;
    }
  } else {
    _if__result_3539 = 0;
  }
  if (_if__result_3539) {
    moonbit_incref(_M0L3srcS133);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS129, 0, _M0L3srcS133, _M0L11src__offsetS131, _M0L11dst__offsetS132, _M0L3lenS130);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS134;
    int32_t _M0L6_2atmpS1496;
    moonbit_string_t _M0L6_2atmpS1495;
    struct _M0TP26biolab8bio__seq9SeqRecord** _result_3540;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS134
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS134, (moonbit_string_t)moonbit_string_literal_81.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS134, _M0L13allocate__lenS129);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS134, (moonbit_string_t)moonbit_string_literal_82.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS134, _M0L11src__offsetS131);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS134, (moonbit_string_t)moonbit_string_literal_83.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS134, _M0L11dst__offsetS132);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS134, (moonbit_string_t)moonbit_string_literal_84.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS134, _M0L3lenS130);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS134, (moonbit_string_t)moonbit_string_literal_85.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1496
    = _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(_M0L3srcS133);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS134, _M0L6_2atmpS1496);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1495
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS134);
    moonbit_decref(_M0L18_2astring__builderS134);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3540
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRP26biolab8bio__seq9SeqRecordEE(_M0L6_2atmpS1495);
    moonbit_decref(_M0L6_2atmpS1495);
    return _result_3540;
  }
}

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23make__and__blit_2einnerGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3srcS139,
  int32_t _M0L13allocate__lenS135,
  int32_t _M0L3lenS136,
  int32_t _M0L11src__offsetS137,
  int32_t _M0L11dst__offsetS138
) {
  int32_t _if__result_3541;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS135 >= 0) {
    if (_M0L3lenS136 >= 0) {
      if (_M0L11src__offsetS137 >= 0) {
        if (_M0L11dst__offsetS138 >= 0) {
          int32_t _M0L6_2atmpS1498 = _M0L11src__offsetS137 + _M0L3lenS136;
          int32_t _M0L6_2atmpS1499;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS1499
          = _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(_M0L3srcS139);
          if (_M0L6_2atmpS1498 <= _M0L6_2atmpS1499) {
            int32_t _M0L6_2atmpS1497 = _M0L11dst__offsetS138 + _M0L3lenS136;
            _if__result_3541 = _M0L6_2atmpS1497 <= _M0L13allocate__lenS135;
          } else {
            _if__result_3541 = 0;
          }
        } else {
          _if__result_3541 = 0;
        }
      } else {
        _if__result_3541 = 0;
      }
    } else {
      _if__result_3541 = 0;
    }
  } else {
    _if__result_3541 = 0;
  }
  if (_if__result_3541) {
    moonbit_incref(_M0L3srcS139);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return _M0MPB18UninitializedArray23unsafe__make__and__blitGRPC16string10StringViewE(_M0L3srcS139, _M0L13allocate__lenS135, _M0L11src__offsetS137, _M0L11dst__offsetS138, _M0L3lenS136);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS140;
    int32_t _M0L6_2atmpS1501;
    moonbit_string_t _M0L6_2atmpS1500;
    struct _M0TPC16string10StringView* _result_3542;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS140
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS140, (moonbit_string_t)moonbit_string_literal_81.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS140, _M0L13allocate__lenS135);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS140, (moonbit_string_t)moonbit_string_literal_82.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS140, _M0L11src__offsetS137);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS140, (moonbit_string_t)moonbit_string_literal_83.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS140, _M0L11dst__offsetS138);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS140, (moonbit_string_t)moonbit_string_literal_84.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS140, _M0L3lenS136);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS140, (moonbit_string_t)moonbit_string_literal_85.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1501
    = _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(_M0L3srcS139);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS140, _M0L6_2atmpS1501);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS1500
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS140);
    moonbit_decref(_M0L18_2astring__builderS140);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_3542
    = _M0FPC15abort5abortGRPB18UninitializedArrayGRPC16string10StringViewEE(_M0L6_2atmpS1500);
    moonbit_decref(_M0L6_2atmpS1500);
    return _result_3542;
  }
}

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder* _M0L4selfS114,
  moonbit_string_t _M0L3objS113
) {
  struct _M0TPB6Logger _M0L6_2atmpS1480;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS114);
  _M0L6_2atmpS1480
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS114
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS113, _M0L6_2atmpS1480);
  if (_M0L6_2atmpS1480.$1) {
    moonbit_decref(_M0L6_2atmpS1480.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder* _M0L4selfS116,
  int32_t _M0L3objS115
) {
  struct _M0TPB6Logger _M0L6_2atmpS1481;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS116);
  _M0L6_2atmpS1481
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS116
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L3objS115, _M0L6_2atmpS1481);
  if (_M0L6_2atmpS1481.$1) {
    moonbit_decref(_M0L6_2atmpS1481.$1);
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
  = (moonbit_string_t*)moonbit_make_ref_array(_M0L13allocate__lenS90, (moonbit_string_t)moonbit_string_literal_0.data);
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

struct _M0TP26biolab8bio__seq9SeqRecord** _M0MPB18UninitializedArray23unsafe__make__and__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS104,
  int32_t _M0L13allocate__lenS102,
  int32_t _M0L11src__offsetS105,
  int32_t _M0L11dst__offsetS103,
  int32_t _M0L9blit__lenS106
) {
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS101;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS101
  = (struct _M0TP26biolab8bio__seq9SeqRecord**)moonbit_make_ref_array(_M0L13allocate__lenS102, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(_M0L3dstS101, _M0L11dst__offsetS103, _M0L3srcS104, _M0L11src__offsetS105, _M0L9blit__lenS106);
  moonbit_decref(_M0L3srcS104);
  return _M0L3dstS101;
}

struct _M0TPC16string10StringView* _M0MPB18UninitializedArray23unsafe__make__and__blitGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3srcS110,
  int32_t _M0L13allocate__lenS108,
  int32_t _M0L11src__offsetS111,
  int32_t _M0L11dst__offsetS109,
  int32_t _M0L9blit__lenS112
) {
  struct _M0TPC16string10StringView* _M0L3dstS107;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS107
  = (struct _M0TPC16string10StringView*)moonbit_make_ref_valtype_array(_M0L13allocate__lenS108, sizeof(struct _M0TPC16string10StringView), Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 98, 0), &(struct _M0TPC16string10StringView){.$0 = (moonbit_string_t)moonbit_string_literal_0.data, .$1 = 0, .$2 = 0});
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGRPC16string10StringViewE(_M0L3dstS107, _M0L11dst__offsetS109, _M0L3srcS110, _M0L11src__offsetS111, _M0L9blit__lenS112);
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS79,
  int32_t _M0L11dst__offsetS80,
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS81,
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

int32_t _M0MPB18UninitializedArray12unsafe__blitGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L3dstS84,
  int32_t _M0L11dst__offsetS85,
  struct _M0TPC16string10StringView* _M0L3srcS86,
  int32_t _M0L11src__offsetS87,
  int32_t _M0L3lenS88
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS86);
  moonbit_incref(_M0L3dstS84);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPC16string10StringViewEE(_M0L3dstS84, _M0L11dst__offsetS85, _M0L3srcS86, _M0L11src__offsetS87, _M0L3lenS88);
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGkE(
  uint16_t* _M0L3dstS24,
  int32_t _M0L11dst__offsetS26,
  uint16_t* _M0L3srcS25,
  int32_t _M0L11src__offsetS27,
  int32_t _M0L3lenS29
) {
  int32_t _if__result_3543;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS24 == _M0L3srcS25) {
    _if__result_3543 = _M0L11dst__offsetS26 < _M0L11src__offsetS27;
  } else {
    _if__result_3543 = 0;
  }
  if (_if__result_3543) {
    int32_t _M0L1iS28 = 0;
    while (1) {
      if (_M0L1iS28 < _M0L3lenS29) {
        int32_t _M0L6_2atmpS1435 = _M0L11dst__offsetS26 + _M0L1iS28;
        int32_t _M0L6_2atmpS1437 = _M0L11src__offsetS27 + _M0L1iS28;
        int32_t _M0L6_2atmpS1436;
        int32_t _M0L6_2atmpS1438;
        if (
          _M0L6_2atmpS1437 < 0
          || _M0L6_2atmpS1437 >= Moonbit_array_length(_M0L3srcS25)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1436 = (int32_t)_M0L3srcS25[_M0L6_2atmpS1437];
        if (
          _M0L6_2atmpS1435 < 0
          || _M0L6_2atmpS1435 >= Moonbit_array_length(_M0L3dstS24)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS24[_M0L6_2atmpS1435] = _M0L6_2atmpS1436;
        _M0L6_2atmpS1438 = _M0L1iS28 + 1;
        _M0L1iS28 = _M0L6_2atmpS1438;
        continue;
      } else {
        moonbit_decref(_M0L3srcS25);
        moonbit_decref(_M0L3dstS24);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1443 = _M0L3lenS29 - 1;
    int32_t _M0L1iS31 = _M0L6_2atmpS1443;
    while (1) {
      if (_M0L1iS31 >= 0) {
        int32_t _M0L6_2atmpS1439 = _M0L11dst__offsetS26 + _M0L1iS31;
        int32_t _M0L6_2atmpS1441 = _M0L11src__offsetS27 + _M0L1iS31;
        int32_t _M0L6_2atmpS1440;
        int32_t _M0L6_2atmpS1442;
        if (
          _M0L6_2atmpS1441 < 0
          || _M0L6_2atmpS1441 >= Moonbit_array_length(_M0L3srcS25)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1440 = (int32_t)_M0L3srcS25[_M0L6_2atmpS1441];
        if (
          _M0L6_2atmpS1439 < 0
          || _M0L6_2atmpS1439 >= Moonbit_array_length(_M0L3dstS24)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS24[_M0L6_2atmpS1439] = _M0L6_2atmpS1440;
        _M0L6_2atmpS1442 = _M0L1iS31 - 1;
        _M0L1iS31 = _M0L6_2atmpS1442;
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
  int32_t _if__result_3546;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS33 == _M0L3srcS34) {
    _if__result_3546 = _M0L11dst__offsetS35 < _M0L11src__offsetS36;
  } else {
    _if__result_3546 = 0;
  }
  if (_if__result_3546) {
    int32_t _M0L1iS37 = 0;
    while (1) {
      if (_M0L1iS37 < _M0L3lenS38) {
        int32_t _M0L6_2atmpS1444 = _M0L11dst__offsetS35 + _M0L1iS37;
        int32_t _M0L6_2atmpS1446 = _M0L11src__offsetS36 + _M0L1iS37;
        moonbit_string_t _M0L6_2atmpS1445;
        moonbit_string_t _M0L6_2aoldS3193;
        int32_t _M0L6_2atmpS1447;
        if (
          _M0L6_2atmpS1446 < 0
          || _M0L6_2atmpS1446 >= Moonbit_array_length(_M0L3srcS34)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1445 = (moonbit_string_t)_M0L3srcS34[_M0L6_2atmpS1446];
        if (
          _M0L6_2atmpS1444 < 0
          || _M0L6_2atmpS1444 >= Moonbit_array_length(_M0L3dstS33)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3193 = (moonbit_string_t)_M0L3dstS33[_M0L6_2atmpS1444];
        moonbit_incref(_M0L6_2atmpS1445);
        moonbit_decref(_M0L6_2aoldS3193);
        _M0L3dstS33[_M0L6_2atmpS1444] = _M0L6_2atmpS1445;
        _M0L6_2atmpS1447 = _M0L1iS37 + 1;
        _M0L1iS37 = _M0L6_2atmpS1447;
        continue;
      } else {
        moonbit_decref(_M0L3srcS34);
        moonbit_decref(_M0L3dstS33);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1452 = _M0L3lenS38 - 1;
    int32_t _M0L1iS40 = _M0L6_2atmpS1452;
    while (1) {
      if (_M0L1iS40 >= 0) {
        int32_t _M0L6_2atmpS1448 = _M0L11dst__offsetS35 + _M0L1iS40;
        int32_t _M0L6_2atmpS1450 = _M0L11src__offsetS36 + _M0L1iS40;
        moonbit_string_t _M0L6_2atmpS1449;
        moonbit_string_t _M0L6_2aoldS3195;
        int32_t _M0L6_2atmpS1451;
        if (
          _M0L6_2atmpS1450 < 0
          || _M0L6_2atmpS1450 >= Moonbit_array_length(_M0L3srcS34)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1449 = (moonbit_string_t)_M0L3srcS34[_M0L6_2atmpS1450];
        if (
          _M0L6_2atmpS1448 < 0
          || _M0L6_2atmpS1448 >= Moonbit_array_length(_M0L3dstS33)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3195 = (moonbit_string_t)_M0L3dstS33[_M0L6_2atmpS1448];
        moonbit_incref(_M0L6_2atmpS1449);
        moonbit_decref(_M0L6_2aoldS3195);
        _M0L3dstS33[_M0L6_2atmpS1448] = _M0L6_2atmpS1449;
        _M0L6_2atmpS1451 = _M0L1iS40 - 1;
        _M0L1iS40 = _M0L6_2atmpS1451;
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
  int32_t _if__result_3549;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS42 == _M0L3srcS43) {
    _if__result_3549 = _M0L11dst__offsetS44 < _M0L11src__offsetS45;
  } else {
    _if__result_3549 = 0;
  }
  if (_if__result_3549) {
    int32_t _M0L1iS46 = 0;
    while (1) {
      if (_M0L1iS46 < _M0L3lenS47) {
        int32_t _M0L6_2atmpS1453 = _M0L11dst__offsetS44 + _M0L1iS46;
        int32_t _M0L6_2atmpS1455 = _M0L11src__offsetS45 + _M0L1iS46;
        struct _M0TUsiE* _M0L6_2atmpS1454;
        struct _M0TUsiE* _M0L6_2aoldS3197;
        int32_t _M0L6_2atmpS1456;
        if (
          _M0L6_2atmpS1455 < 0
          || _M0L6_2atmpS1455 >= Moonbit_array_length(_M0L3srcS43)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1454 = (struct _M0TUsiE*)_M0L3srcS43[_M0L6_2atmpS1455];
        if (
          _M0L6_2atmpS1453 < 0
          || _M0L6_2atmpS1453 >= Moonbit_array_length(_M0L3dstS42)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3197 = (struct _M0TUsiE*)_M0L3dstS42[_M0L6_2atmpS1453];
        if (_M0L6_2atmpS1454) {
          moonbit_incref(_M0L6_2atmpS1454);
        }
        if (_M0L6_2aoldS3197) {
          moonbit_decref(_M0L6_2aoldS3197);
        }
        _M0L3dstS42[_M0L6_2atmpS1453] = _M0L6_2atmpS1454;
        _M0L6_2atmpS1456 = _M0L1iS46 + 1;
        _M0L1iS46 = _M0L6_2atmpS1456;
        continue;
      } else {
        moonbit_decref(_M0L3srcS43);
        moonbit_decref(_M0L3dstS42);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1461 = _M0L3lenS47 - 1;
    int32_t _M0L1iS49 = _M0L6_2atmpS1461;
    while (1) {
      if (_M0L1iS49 >= 0) {
        int32_t _M0L6_2atmpS1457 = _M0L11dst__offsetS44 + _M0L1iS49;
        int32_t _M0L6_2atmpS1459 = _M0L11src__offsetS45 + _M0L1iS49;
        struct _M0TUsiE* _M0L6_2atmpS1458;
        struct _M0TUsiE* _M0L6_2aoldS3199;
        int32_t _M0L6_2atmpS1460;
        if (
          _M0L6_2atmpS1459 < 0
          || _M0L6_2atmpS1459 >= Moonbit_array_length(_M0L3srcS43)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1458 = (struct _M0TUsiE*)_M0L3srcS43[_M0L6_2atmpS1459];
        if (
          _M0L6_2atmpS1457 < 0
          || _M0L6_2atmpS1457 >= Moonbit_array_length(_M0L3dstS42)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3199 = (struct _M0TUsiE*)_M0L3dstS42[_M0L6_2atmpS1457];
        if (_M0L6_2atmpS1458) {
          moonbit_incref(_M0L6_2atmpS1458);
        }
        if (_M0L6_2aoldS3199) {
          moonbit_decref(_M0L6_2aoldS3199);
        }
        _M0L3dstS42[_M0L6_2atmpS1457] = _M0L6_2atmpS1458;
        _M0L6_2atmpS1460 = _M0L1iS49 - 1;
        _M0L1iS49 = _M0L6_2atmpS1460;
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRP26biolab8bio__seq9SeqRecordEE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3dstS51,
  int32_t _M0L11dst__offsetS53,
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L3srcS52,
  int32_t _M0L11src__offsetS54,
  int32_t _M0L3lenS56
) {
  int32_t _if__result_3552;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS51 == _M0L3srcS52) {
    _if__result_3552 = _M0L11dst__offsetS53 < _M0L11src__offsetS54;
  } else {
    _if__result_3552 = 0;
  }
  if (_if__result_3552) {
    int32_t _M0L1iS55 = 0;
    while (1) {
      if (_M0L1iS55 < _M0L3lenS56) {
        int32_t _M0L6_2atmpS1462 = _M0L11dst__offsetS53 + _M0L1iS55;
        int32_t _M0L6_2atmpS1464 = _M0L11src__offsetS54 + _M0L1iS55;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1463;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3201;
        int32_t _M0L6_2atmpS1465;
        if (
          _M0L6_2atmpS1464 < 0
          || _M0L6_2atmpS1464 >= Moonbit_array_length(_M0L3srcS52)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1463
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3srcS52[
            _M0L6_2atmpS1464
          ];
        if (
          _M0L6_2atmpS1462 < 0
          || _M0L6_2atmpS1462 >= Moonbit_array_length(_M0L3dstS51)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3201
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3dstS51[
            _M0L6_2atmpS1462
          ];
        if (_M0L6_2atmpS1463) {
          moonbit_incref(_M0L6_2atmpS1463);
        }
        if (_M0L6_2aoldS3201) {
          moonbit_decref(_M0L6_2aoldS3201);
        }
        _M0L3dstS51[_M0L6_2atmpS1462] = _M0L6_2atmpS1463;
        _M0L6_2atmpS1465 = _M0L1iS55 + 1;
        _M0L1iS55 = _M0L6_2atmpS1465;
        continue;
      } else {
        moonbit_decref(_M0L3srcS52);
        moonbit_decref(_M0L3dstS51);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1470 = _M0L3lenS56 - 1;
    int32_t _M0L1iS58 = _M0L6_2atmpS1470;
    while (1) {
      if (_M0L1iS58 >= 0) {
        int32_t _M0L6_2atmpS1466 = _M0L11dst__offsetS53 + _M0L1iS58;
        int32_t _M0L6_2atmpS1468 = _M0L11src__offsetS54 + _M0L1iS58;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2atmpS1467;
        struct _M0TP26biolab8bio__seq9SeqRecord* _M0L6_2aoldS3203;
        int32_t _M0L6_2atmpS1469;
        if (
          _M0L6_2atmpS1468 < 0
          || _M0L6_2atmpS1468 >= Moonbit_array_length(_M0L3srcS52)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1467
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3srcS52[
            _M0L6_2atmpS1468
          ];
        if (
          _M0L6_2atmpS1466 < 0
          || _M0L6_2atmpS1466 >= Moonbit_array_length(_M0L3dstS51)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3203
        = (struct _M0TP26biolab8bio__seq9SeqRecord*)_M0L3dstS51[
            _M0L6_2atmpS1466
          ];
        if (_M0L6_2atmpS1467) {
          moonbit_incref(_M0L6_2atmpS1467);
        }
        if (_M0L6_2aoldS3203) {
          moonbit_decref(_M0L6_2aoldS3203);
        }
        _M0L3dstS51[_M0L6_2atmpS1466] = _M0L6_2atmpS1467;
        _M0L6_2atmpS1469 = _M0L1iS58 - 1;
        _M0L1iS58 = _M0L6_2atmpS1469;
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

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGRPC16string10StringViewEE(
  struct _M0TPC16string10StringView* _M0L3dstS60,
  int32_t _M0L11dst__offsetS62,
  struct _M0TPC16string10StringView* _M0L3srcS61,
  int32_t _M0L11src__offsetS63,
  int32_t _M0L3lenS65
) {
  int32_t _if__result_3555;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS60 == _M0L3srcS61) {
    _if__result_3555 = _M0L11dst__offsetS62 < _M0L11src__offsetS63;
  } else {
    _if__result_3555 = 0;
  }
  if (_if__result_3555) {
    int32_t _M0L1iS64 = 0;
    while (1) {
      if (_M0L1iS64 < _M0L3lenS65) {
        int32_t _M0L6_2atmpS1471 = _M0L11dst__offsetS62 + _M0L1iS64;
        int32_t _M0L6_2atmpS1473 = _M0L11src__offsetS63 + _M0L1iS64;
        struct _M0TPC16string10StringView _M0L6_2atmpS1472;
        struct _M0TPC16string10StringView _M0L6_2aoldS3205;
        int32_t _M0L6_2atmpS1474;
        if (
          _M0L6_2atmpS1473 < 0
          || _M0L6_2atmpS1473 >= Moonbit_array_length(_M0L3srcS61)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1472 = _M0L3srcS61[_M0L6_2atmpS1473];
        if (
          _M0L6_2atmpS1471 < 0
          || _M0L6_2atmpS1471 >= Moonbit_array_length(_M0L3dstS60)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3205 = _M0L3dstS60[_M0L6_2atmpS1471];
        moonbit_incref(_M0L6_2atmpS1472.$0);
        moonbit_decref(_M0L6_2aoldS3205.$0);
        _M0L3dstS60[_M0L6_2atmpS1471] = _M0L6_2atmpS1472;
        _M0L6_2atmpS1474 = _M0L1iS64 + 1;
        _M0L1iS64 = _M0L6_2atmpS1474;
        continue;
      } else {
        moonbit_decref(_M0L3srcS61);
        moonbit_decref(_M0L3dstS60);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS1479 = _M0L3lenS65 - 1;
    int32_t _M0L1iS67 = _M0L6_2atmpS1479;
    while (1) {
      if (_M0L1iS67 >= 0) {
        int32_t _M0L6_2atmpS1475 = _M0L11dst__offsetS62 + _M0L1iS67;
        int32_t _M0L6_2atmpS1477 = _M0L11src__offsetS63 + _M0L1iS67;
        struct _M0TPC16string10StringView _M0L6_2atmpS1476;
        struct _M0TPC16string10StringView _M0L6_2aoldS3207;
        int32_t _M0L6_2atmpS1478;
        if (
          _M0L6_2atmpS1477 < 0
          || _M0L6_2atmpS1477 >= Moonbit_array_length(_M0L3srcS61)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS1476 = _M0L3srcS61[_M0L6_2atmpS1477];
        if (
          _M0L6_2atmpS1475 < 0
          || _M0L6_2atmpS1475 >= Moonbit_array_length(_M0L3dstS60)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS3207 = _M0L3dstS60[_M0L6_2atmpS1475];
        moonbit_incref(_M0L6_2atmpS1476.$0);
        moonbit_decref(_M0L6_2aoldS3207.$0);
        _M0L3dstS60[_M0L6_2atmpS1475] = _M0L6_2atmpS1476;
        _M0L6_2atmpS1478 = _M0L1iS67 - 1;
        _M0L1iS67 = _M0L6_2atmpS1478;
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

int32_t _M0MPB18UninitializedArray6lengthGRP26biolab8bio__seq9SeqRecordE(
  struct _M0TP26biolab8bio__seq9SeqRecord** _M0L4selfS22
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS22);
}

int32_t _M0MPB18UninitializedArray6lengthGRPC16string10StringViewE(
  struct _M0TPC16string10StringView* _M0L4selfS23
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS23);
}

uint32_t _M0FPB13consume4__acc(uint32_t _M0L3accS18, uint32_t _M0L5inputS19) {
  uint32_t _M0L6_2atmpS1434;
  uint32_t _M0L6_2atmpS1433;
  uint32_t _M0L6_2atmpS1432;
  #line 465 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1434 = _M0L5inputS19 * 3266489917u;
  _M0L6_2atmpS1433 = _M0L3accS18 + _M0L6_2atmpS1434;
  #line 466 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1432 = _M0FPB4rotl(_M0L6_2atmpS1433, 17);
  return _M0L6_2atmpS1432 * 668265263u;
}

uint32_t _M0FPB4rotl(uint32_t _M0L1xS16, int32_t _M0L1rS17) {
  uint32_t _M0L6_2atmpS1429;
  int32_t _M0L6_2atmpS1431;
  uint32_t _M0L6_2atmpS1430;
  #line 475 "/home/developer/.moon/lib/core/builtin/hasher.mbt"
  _M0L6_2atmpS1429 = _M0L1xS16 << (_M0L1rS17 & 31);
  _M0L6_2atmpS1431 = 32 - _M0L1rS17;
  _M0L6_2atmpS1430 = _M0L1xS16 >> (_M0L6_2atmpS1431 & 31);
  return _M0L6_2atmpS1429 | _M0L6_2atmpS1430;
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
  _M0L10_2ax__5722S15.$0->$method_0(_M0L10_2ax__5722S15.$1, (moonbit_string_t)moonbit_string_literal_86.data);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB6Logger13write__objectGsE(_M0L10_2ax__5722S15, _M0L15_2a_2aarg__5723S14);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S15.$0->$method_0(_M0L10_2ax__5722S15.$1, (moonbit_string_t)moonbit_string_literal_87.data);
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

struct _M0TPC16string10StringView* _M0FPC15abort5abortGRPB18UninitializedArrayGRPC16string10StringViewEE(
  moonbit_string_t _M0L3msgS9
) {
  #line 47 "/home/developer/.moon/lib/core/abort/abort.mbt"
  #line 49 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_println(_M0L3msgS9);
  #line 50 "/home/developer/.moon/lib/core/abort/abort.mbt"
  moonbit_panic();
}

moonbit_string_t _M0FP15Error10to__string(void* _M0L4_2aeS1400) {
  switch (Moonbit_object_tag(_M0L4_2aeS1400)) {
    case 3: {
      return (moonbit_string_t)moonbit_string_literal_88.data;
      break;
    }
    
    case 2: {
      return (moonbit_string_t)moonbit_string_literal_89.data;
      break;
    }
    
    case 0: {
      return _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(_M0L4_2aeS1400);
      break;
    }
    
    case 4: {
      return (moonbit_string_t)moonbit_string_literal_90.data;
      break;
    }
    
    case 1: {
      return (moonbit_string_t)moonbit_string_literal_91.data;
      break;
    }
    default: {
      return (moonbit_string_t)moonbit_string_literal_92.data;
      break;
    }
  }
}

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1423,
  struct _M0TPB4Show _M0L8_2aparamS1422
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1421 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1423;
  _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(_M0L7_2aselfS1421, _M0L8_2aparamS1422);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1420,
  struct _M0TPB4Show _M0L8_2aparamS1419
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1418 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1420;
  _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(_M0L7_2aselfS1418, _M0L8_2aparamS1419);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1417,
  int32_t _M0L8_2aparamS1416
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1415 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1417;
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS1415, _M0L8_2aparamS1416);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1414,
  struct _M0TPC16string10StringView _M0L8_2aparamS1413
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1412 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1414;
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L7_2aselfS1412, _M0L8_2aparamS1413);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS1411,
  moonbit_string_t _M0L8_2aparamS1408,
  int32_t _M0L8_2aparamS1409,
  int32_t _M0L8_2aparamS1410
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1407 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1411;
  _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L7_2aselfS1407, _M0L8_2aparamS1408, _M0L8_2aparamS1409, _M0L8_2aparamS1410);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS1406,
  moonbit_string_t _M0L8_2aparamS1405
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS1404 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS1406;
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L7_2aselfS1404, _M0L8_2aparamS1405);
  return 0;
}

void moonbit_init() {
  moonbit_layout_table = moonbit_layout_table_data;
}

int main(int argc, char** argv) {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** _M0L6_2atmpS1428;
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS1394;
  struct _M0TPB5ArrayGUsiEE* _M0L7_2abindS1395;
  int32_t _M0L7_2abindS1396;
  int32_t _M0L2__S1397;
  moonbit_runtime_init(argc, argv);
  moonbit_init();
  _M0L6_2atmpS1428
  = (struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error**)moonbit_empty_ref_array;
  _M0L12async__testsS1394
  = (struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE));
  Moonbit_object_header(_M0L12async__testsS1394)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 144, 0);
  _M0L12async__testsS1394->$0 = _M0L6_2atmpS1428;
  _M0L12async__testsS1394->$1 = 0;
  #line 399 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0L7_2abindS1395
  = _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test52moonbit__test__driver__internal__native__parse__args();
  _M0L7_2abindS1396 = _M0L7_2abindS1395->$1;
  _M0L2__S1397 = 0;
  while (1) {
    if (_M0L2__S1397 < _M0L7_2abindS1396) {
      struct _M0TUsiE** _M0L3bufS1427 = _M0L7_2abindS1395->$0;
      struct _M0TUsiE* _M0L3argS1398 =
        (struct _M0TUsiE*)_M0L3bufS1427[_M0L2__S1397];
      moonbit_string_t _M0L6_2atmpS1424 = _M0L3argS1398->$0;
      int32_t _M0L6_2atmpS1425 = _M0L3argS1398->$1;
      int32_t _M0L6_2atmpS1426;
      moonbit_incref(_M0L6_2atmpS1424);
      #line 400 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
      _M0FP46biolab8bio__seq3cmd29alignio__main__blackbox__test44moonbit__test__driver__internal__do__execute(_M0L12async__testsS1394, _M0L6_2atmpS1424, _M0L6_2atmpS1425);
      moonbit_decref(_M0L6_2atmpS1424);
      _M0L6_2atmpS1426 = _M0L2__S1397 + 1;
      _M0L2__S1397 = _M0L6_2atmpS1426;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS1395);
    }
    break;
  }
  #line 402 "/home/developer/Documents/1/cmd/alignio_main/__generated_driver_for_blackbox_test.mbt"
  _M0IP016_24default__implP46biolab8bio__seq3cmd29alignio__main__blackbox__test28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd29alignio__main__blackbox__test34MoonBit__Async__Test__Driver__ImplE(_M0L12async__testsS1394);
  moonbit_decref(_M0L12async__testsS1394);
  return 0;
}