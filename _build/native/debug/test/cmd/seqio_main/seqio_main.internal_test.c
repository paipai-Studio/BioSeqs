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
struct _M0TURPC16string10StringViewRPB6LoggerE;

struct _M0TPB8MutLocalGiE;

struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362;

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure;

struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest;

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError;

struct _M0TWRPC15error5ErrorEs;

struct _M0TPB4Show;

struct _M0TWssbEu;

struct _M0TUsiE;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd11seqio__main33MoonBitTestDriverInternalSkipTestE3Err;

struct _M0TPB13StringBuilder;

struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357;

struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error;

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd11seqio__main33MoonBitTestDriverInternalSkipTestE2Ok;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err;

struct _M0BTPB6Logger;

struct _M0BTPB4Show;

struct _M0TWuEu;

struct _M0TPC16string10StringView;

struct _M0KTPB6LoggerTPB13StringBuilder;

struct _M0TPB6Logger;

struct _M0TPB5ArrayGUsiEE;

struct _M0TPB5ArrayGsE;

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok;

struct _M0TWEu;

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE;

struct _M0DTPC15error5Error107biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError;

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError;

struct _M0TWRPC15error5ErrorEu;

struct _M0TPB6Logger {
  struct _M0BTPB6Logger* $0;
  void* $1;
  
};

struct _M0TPC16string10StringView {
  moonbit_string_t $0;
  int32_t $1;
  int32_t $2;
  
};

struct _M0TURPC16string10StringViewRPB6LoggerE {
  struct _M0TPC16string10StringView $0;
  struct _M0TPB6Logger $1;
  
};

struct _M0TPB8MutLocalGiE {
  int32_t $0;
  
};

struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362 {
  int32_t(* code)(
    struct _M0TWssbEu*,
    moonbit_string_t,
    moonbit_string_t,
    int32_t
  );
  int32_t $0;
  moonbit_string_t $1;
  
};

struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure {
  moonbit_string_t $0;
  
};

struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest {
  moonbit_string_t $0;
  
};

struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError {
  moonbit_string_t $0;
  
};

struct _M0TWRPC15error5ErrorEs {
  moonbit_string_t(* code)(struct _M0TWRPC15error5ErrorEs*, void*);
  
};

struct _M0TPB4Show {
  struct _M0BTPB4Show* $0;
  void* $1;
  
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

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd11seqio__main33MoonBitTestDriverInternalSkipTestE3Err {
  void* $0;
  
};

struct _M0TPB13StringBuilder {
  uint16_t* $0;
  int32_t $1;
  
};

struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357 {
  int32_t(* code)(struct _M0TWEu*);
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

struct _M0DTPC16result6ResultGbRP46biolab8bio__seq3cmd11seqio__main33MoonBitTestDriverInternalSkipTestE2Ok {
  int32_t $0;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE3Err {
  void* $0;
  
};

struct _M0BTPB6Logger {
  int32_t(* $method_0)(void*, moonbit_string_t);
  int32_t(* $method_1)(void*, moonbit_string_t, int32_t, int32_t);
  int32_t(* $method_2)(void*, struct _M0TPC16string10StringView);
  int32_t(* $method_3)(void*, int32_t);
  int32_t(* $method_4)(void*, struct _M0TPB4Show);
  int32_t(* $method_5)(void*, struct _M0TPB4Show);
  
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

struct _M0TPB5ArrayGUsiEE {
  struct _M0TUsiE** $0;
  int32_t $1;
  
};

struct _M0TPB5ArrayGsE {
  moonbit_string_t* $0;
  int32_t $1;
  
};

struct _M0DTPC16result6ResultGOuRPC15error5ErrorE2Ok {
  int32_t $0;
  
};

struct _M0TWEu {
  int32_t(* code)(struct _M0TWEu*);
  
};

struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** $0;
  int32_t $1;
  
};

struct _M0DTPC15error5Error107biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError {
  moonbit_string_t $0;
  
};

struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError {
  moonbit_string_t $0;
  
};

struct _M0TWRPC15error5ErrorEu {
  int32_t(* code)(struct _M0TWRPC15error5ErrorEu*, void*);
  
};

struct moonbit_result_0 {
  int tag;
  union { int32_t ok; void* err;  } data;
  
};

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd11seqio__main34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*
);

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__args(
  
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS412(
  int32_t,
  moonbit_string_t,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS405(
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398(
  int32_t,
  moonbit_bytes_t
);

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S391(
  int32_t,
  moonbit_string_t
);

#define _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__get__cli__args__ffi moonbit_get_cli_args

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t
);

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN17error__to__stringS369(
  struct _M0TWRPC15error5ErrorEs*,
  void*
);

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN14handle__resultS362(
  struct _M0TWssbEu*,
  moonbit_string_t,
  moonbit_string_t,
  int32_t
);

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN13handle__startS357(
  struct _M0TWEu*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*,
  moonbit_string_t,
  int32_t,
  struct _M0TWEu*,
  struct _M0TWssbEu*,
  struct _M0TWRPC15error5ErrorEs*
);

moonbit_string_t _M0MPC15array5Array2atGsE(struct _M0TPB5ArrayGsE*, int32_t);

int32_t _M0FPB7printlnGsE(moonbit_string_t);

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t);

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE*,
  moonbit_string_t
);

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE*,
  struct _M0TUsiE*
);

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE*);

int32_t _M0MPC15array5Array7reallocGUsiEE(struct _M0TPB5ArrayGUsiEE*);

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE*,
  int32_t
);

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE*,
  int32_t
);

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(int32_t);

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(moonbit_string_t);

int32_t _M0IPB13StringBuilderPB6Logger11write__view(
  struct _M0TPB13StringBuilder*,
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

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t);

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t);

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

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t*);

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(struct _M0TUsiE**);

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

struct { int32_t rc; uint32_t meta; uint16_t const data[35]; 
} const moonbit_string_literal_2 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 34, 45, 45, 
    45, 45, 45, 32, 66, 69, 71, 73, 78, 32, 77, 79, 79, 78, 32, 84, 69, 
    83, 84, 32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
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
} const moonbit_string_literal_15 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 116, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_13 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 114, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_20 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    100, 115, 116, 95, 111, 102, 102, 115, 101, 116, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[19]; 
} const moonbit_string_literal_17 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 18, 105, 110, 
    118, 97, 108, 105, 100, 32, 99, 111, 100, 101, 32, 112, 111, 105, 
    110, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[12]; 
} const moonbit_string_literal_5 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 11, 44, 34, 
    109, 101, 115, 115, 97, 103, 101, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[53]; 
} const moonbit_string_literal_26 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 52, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 46, 83, 110, 97, 112, 115, 
    104, 111, 116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_12 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 110, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[31]; 
} const moonbit_string_literal_9 =
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

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_10 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 48, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_21 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 44, 32, 
    108, 101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_18 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 98, 111, 
    117, 110, 100, 115, 32, 99, 104, 101, 99, 107, 32, 102, 97, 105, 
    108, 101, 100, 58, 32, 97, 108, 108, 111, 99, 97, 116, 101, 95, 108, 
    101, 110, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[4]; 
} const moonbit_string_literal_16 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 3, 92, 117, 123, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[96]; 
} const moonbit_string_literal_27 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 95, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 115, 101, 113, 105, 111, 95, 109, 97, 105, 110, 46, 77, 
    111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 
    101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 74, 115, 69, 114, 
    114, 111, 114, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 115, 
    116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 97, 
    108, 74, 115, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_24 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 41, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[37]; 
} const moonbit_string_literal_11 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 36, 48, 49, 
    50, 51, 52, 53, 54, 55, 56, 57, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 
    118, 119, 120, 121, 122, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[3]; 
} const moonbit_string_literal_14 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 2, 92, 98, 0};

struct { int32_t rc; uint32_t meta; uint16_t const data[51]; 
} const moonbit_string_literal_25 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 50, 109, 111, 
    111, 110, 98, 105, 116, 108, 97, 110, 103, 47, 99, 111, 114, 101, 
    47, 98, 117, 105, 108, 116, 105, 110, 46, 73, 110, 115, 112, 101, 
    99, 116, 69, 114, 114, 111, 114, 46, 73, 110, 115, 112, 101, 99, 
    116, 69, 114, 114, 111, 114, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_22 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 15, 44, 32, 
    115, 114, 99, 46, 108, 101, 110, 103, 116, 104, 32, 61, 32, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[33]; 
} const moonbit_string_literal_7 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 32, 45, 45, 
    45, 45, 45, 32, 69, 78, 68, 32, 77, 79, 79, 78, 32, 84, 69, 83, 84, 
    32, 82, 69, 83, 85, 76, 84, 32, 45, 45, 45, 45, 45, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[98]; 
} const moonbit_string_literal_28 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 97, 98, 105, 
    111, 108, 97, 98, 47, 98, 105, 111, 95, 115, 101, 113, 47, 99, 109, 
    100, 47, 115, 101, 113, 105, 111, 95, 109, 97, 105, 110, 46, 77, 
    111, 111, 110, 66, 105, 116, 84, 101, 115, 116, 68, 114, 105, 118, 
    101, 114, 73, 110, 116, 101, 114, 110, 97, 108, 83, 107, 105, 112, 
    84, 101, 115, 116, 46, 77, 111, 111, 110, 66, 105, 116, 84, 101, 
    115, 116, 68, 114, 105, 118, 101, 114, 73, 110, 116, 101, 114, 110, 
    97, 108, 83, 107, 105, 112, 84, 101, 115, 116, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[16]; 
} const moonbit_string_literal_19 =
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

struct { int32_t rc; uint32_t meta; uint16_t const data[10]; 
} const moonbit_string_literal_4 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 9, 44, 34, 
    105, 110, 100, 101, 120, 34, 58, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[9]; 
} const moonbit_string_literal_23 =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 8, 70, 97, 
    105, 108, 117, 114, 101, 40, 0
  };

struct { int32_t rc; uint32_t meta; uint16_t const data[2]; 
} const moonbit_string_literal_6 =
  { Moonbit_make_static_rc(moonbit_BLOCK_KIND_VAL_ARRAY), 1, 125, 0};

struct { int32_t rc; uint32_t meta; struct _M0TWRPC15error5ErrorEs data; 
} const _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN17error__to__stringS369$closure =
  {
    Moonbit_make_static_rc(moonbit_BLOCK_KIND_REGULAR),
    Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0),
    _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN17error__to__stringS369
  };

uint32_t const moonbit_layout_table_data[36] =
  {
    sizeof(struct _M0TPB5ArrayGUsiEE) / 4, 1,
    offsetof(struct _M0TPB5ArrayGUsiEE, $0) / 4, sizeof(struct _M0TUsiE) / 4,
    1, offsetof(struct _M0TUsiE, $0) / 4, sizeof(struct _M0TPB5ArrayGsE) / 4,
    1, offsetof(struct _M0TPB5ArrayGsE, $0) / 4,
    sizeof(struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357)
    / 4, 1,
    offsetof(struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357, $1)
    / 4,
    sizeof(struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362)
    / 4, 1,
    offsetof(struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362, $1)
    / 4,
    sizeof(struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest)
    / 4, 1,
    offsetof(struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest, $0)
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

int32_t _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd11seqio__main34MoonBit__Async__Test__Driver__ImplE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S434
) {
  #line 12 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  return 0;
}

struct _M0TPB5ArrayGUsiEE* _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__args(
  
) {
  int32_t _M0L45moonbit__test__driver__internal__parse__int__S391;
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398;
  int32_t _M0L57moonbit__test__driver__internal__get__cli__args__internalS405;
  int32_t _M0L51moonbit__test__driver__internal__split__mbt__stringS412;
  struct _M0TUsiE** _M0L6_2atmpS896;
  struct _M0TPB5ArrayGUsiEE* _M0L16file__and__indexS419;
  struct _M0TPB5ArrayGsE* _M0L9cli__argsS420;
  moonbit_string_t _M0L6_2atmpS895;
  struct _M0TPB5ArrayGsE* _M0L10test__argsS421;
  int32_t _M0L7_2abindS422;
  int32_t _M0L2__S423;
  #line 194 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L45moonbit__test__driver__internal__parse__int__S391 = 0;
  _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398 = 0;
  _M0L57moonbit__test__driver__internal__get__cli__args__internalS405
  = _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398;
  _M0L51moonbit__test__driver__internal__split__mbt__stringS412 = 0;
  _M0L6_2atmpS896 = (struct _M0TUsiE**)moonbit_empty_ref_array;
  _M0L16file__and__indexS419
  = (struct _M0TPB5ArrayGUsiEE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGUsiEE));
  Moonbit_object_header(_M0L16file__and__indexS419)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 0, 0);
  _M0L16file__and__indexS419->$0 = _M0L6_2atmpS896;
  _M0L16file__and__indexS419->$1 = 0;
  #line 283 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L9cli__argsS420
  = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS405(_M0L57moonbit__test__driver__internal__get__cli__args__internalS405);
  #line 285 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS895 = _M0MPC15array5Array2atGsE(_M0L9cli__argsS420, 1);
  moonbit_decref(_M0L9cli__argsS420);
  #line 284 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L10test__argsS421
  = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS412(_M0L51moonbit__test__driver__internal__split__mbt__stringS412, _M0L6_2atmpS895, 47);
  moonbit_decref(_M0L6_2atmpS895);
  _M0L7_2abindS422 = _M0L10test__argsS421->$1;
  _M0L2__S423 = 0;
  while (1) {
    if (_M0L2__S423 < _M0L7_2abindS422) {
      moonbit_string_t* _M0L3bufS894 = _M0L10test__argsS421->$0;
      moonbit_string_t _M0L3argS424 =
        (moonbit_string_t)_M0L3bufS894[_M0L2__S423];
      struct _M0TPB5ArrayGsE* _M0L16file__and__rangeS425;
      moonbit_string_t _M0L4fileS426;
      moonbit_string_t _M0L5rangeS427;
      struct _M0TPB5ArrayGsE* _M0L15start__and__endS428;
      moonbit_string_t _M0L6_2atmpS892;
      int32_t _M0L5startS429;
      moonbit_string_t _M0L6_2atmpS891;
      int32_t _M0L3endS430;
      int32_t _M0L1iS431;
      int32_t _M0L6_2atmpS893;
      moonbit_incref(_M0L3argS424);
      #line 289 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L16file__and__rangeS425
      = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS412(_M0L51moonbit__test__driver__internal__split__mbt__stringS412, _M0L3argS424, 58);
      moonbit_decref(_M0L3argS424);
      #line 290 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L4fileS426
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS425, 0);
      #line 291 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L5rangeS427
      = _M0MPC15array5Array2atGsE(_M0L16file__and__rangeS425, 1);
      moonbit_decref(_M0L16file__and__rangeS425);
      #line 292 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L15start__and__endS428
      = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS412(_M0L51moonbit__test__driver__internal__split__mbt__stringS412, _M0L5rangeS427, 45);
      moonbit_decref(_M0L5rangeS427);
      #line 295 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L6_2atmpS892
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS428, 0);
      #line 295 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L5startS429
      = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S391(_M0L45moonbit__test__driver__internal__parse__int__S391, _M0L6_2atmpS892);
      moonbit_decref(_M0L6_2atmpS892);
      #line 296 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L6_2atmpS891
      = _M0MPC15array5Array2atGsE(_M0L15start__and__endS428, 1);
      moonbit_decref(_M0L15start__and__endS428);
      #line 296 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L3endS430
      = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S391(_M0L45moonbit__test__driver__internal__parse__int__S391, _M0L6_2atmpS891);
      moonbit_decref(_M0L6_2atmpS891);
      _M0L1iS431 = _M0L5startS429;
      while (1) {
        if (_M0L1iS431 < _M0L3endS430) {
          struct _M0TUsiE* _M0L8_2atupleS889;
          int32_t _M0L6_2atmpS890;
          moonbit_incref(_M0L4fileS426);
          _M0L8_2atupleS889
          = (struct _M0TUsiE*)moonbit_malloc(sizeof(struct _M0TUsiE));
          Moonbit_object_header(_M0L8_2atupleS889)->meta
          = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 3, 0);
          _M0L8_2atupleS889->$0 = _M0L4fileS426;
          _M0L8_2atupleS889->$1 = _M0L1iS431;
          #line 298 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
          _M0MPC15array5Array4pushGUsiEE(_M0L16file__and__indexS419, _M0L8_2atupleS889);
          moonbit_decref(_M0L8_2atupleS889);
          _M0L6_2atmpS890 = _M0L1iS431 + 1;
          _M0L1iS431 = _M0L6_2atmpS890;
          continue;
        } else {
          moonbit_decref(_M0L4fileS426);
        }
        break;
      }
      _M0L6_2atmpS893 = _M0L2__S423 + 1;
      _M0L2__S423 = _M0L6_2atmpS893;
      continue;
    } else {
      moonbit_decref(_M0L10test__argsS421);
    }
    break;
  }
  return _M0L16file__and__indexS419;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN51moonbit__test__driver__internal__split__mbt__stringS412(
  int32_t _M0L6_2aenvS870,
  moonbit_string_t _M0L1sS413,
  int32_t _M0L3sepS414
) {
  moonbit_string_t* _M0L6_2atmpS888;
  struct _M0TPB5ArrayGsE* _M0L3resS415;
  struct _M0TPB8MutLocalGiE* _M0L1iS416;
  struct _M0TPB8MutLocalGiE* _M0L5startS417;
  int32_t _M0L3valS883;
  int32_t _M0L6_2atmpS884;
  #line 262 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS888 = (moonbit_string_t*)moonbit_empty_ref_array;
  _M0L3resS415
  = (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
  Moonbit_object_header(_M0L3resS415)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 0);
  _M0L3resS415->$0 = _M0L6_2atmpS888;
  _M0L3resS415->$1 = 0;
  _M0L1iS416
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS416)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS416->$0 = 0;
  _M0L5startS417
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L5startS417)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L5startS417->$0 = 0;
  while (1) {
    int32_t _M0L3valS871 = _M0L1iS416->$0;
    int32_t _M0L6_2atmpS872 = Moonbit_array_length(_M0L1sS413);
    if (_M0L3valS871 < _M0L6_2atmpS872) {
      int32_t _M0L3valS875 = _M0L1iS416->$0;
      int32_t _M0L6_2atmpS874;
      int32_t _M0L6_2atmpS873;
      int32_t _M0L3valS882;
      int32_t _M0L6_2atmpS881;
      if (
        _M0L3valS875 < 0 || _M0L3valS875 >= Moonbit_array_length(_M0L1sS413)
      ) {
        #line 270 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS874 = _M0L1sS413[_M0L3valS875];
      _M0L6_2atmpS873 = _M0L6_2atmpS874;
      if (_M0L6_2atmpS873 == _M0L3sepS414) {
        int32_t _M0L3valS877 = _M0L5startS417->$0;
        int32_t _M0L3valS878 = _M0L1iS416->$0;
        moonbit_string_t _M0L6_2atmpS876;
        int32_t _M0L3valS880;
        int32_t _M0L6_2atmpS879;
        #line 271 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        _M0L6_2atmpS876
        = _M0MPC16string6String17unsafe__substring(_M0L1sS413, _M0L3valS877, _M0L3valS878);
        #line 271 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        _M0MPC15array5Array4pushGsE(_M0L3resS415, _M0L6_2atmpS876);
        moonbit_decref(_M0L6_2atmpS876);
        _M0L3valS880 = _M0L1iS416->$0;
        _M0L6_2atmpS879 = _M0L3valS880 + 1;
        _M0L5startS417->$0 = _M0L6_2atmpS879;
      }
      _M0L3valS882 = _M0L1iS416->$0;
      _M0L6_2atmpS881 = _M0L3valS882 + 1;
      _M0L1iS416->$0 = _M0L6_2atmpS881;
      continue;
    } else {
      moonbit_decref(_M0L1iS416);
    }
    break;
  }
  _M0L3valS883 = _M0L5startS417->$0;
  _M0L6_2atmpS884 = Moonbit_array_length(_M0L1sS413);
  if (_M0L3valS883 < _M0L6_2atmpS884) {
    int32_t _M0L3valS886 = _M0L5startS417->$0;
    int32_t _M0L6_2atmpS887;
    moonbit_string_t _M0L6_2atmpS885;
    moonbit_decref(_M0L5startS417);
    _M0L6_2atmpS887 = Moonbit_array_length(_M0L1sS413);
    #line 277 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
    _M0L6_2atmpS885
    = _M0MPC16string6String17unsafe__substring(_M0L1sS413, _M0L3valS886, _M0L6_2atmpS887);
    #line 277 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
    _M0MPC15array5Array4pushGsE(_M0L3resS415, _M0L6_2atmpS885);
    moonbit_decref(_M0L6_2atmpS885);
  } else {
    moonbit_decref(_M0L5startS417);
  }
  return _M0L3resS415;
}

struct _M0TPB5ArrayGsE* _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN57moonbit__test__driver__internal__get__cli__args__internalS405(
  int32_t _M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398
) {
  moonbit_bytes_t* _M0L3tmpS406;
  int32_t _M0L6_2atmpS869;
  struct _M0TPB5ArrayGsE* _M0L3resS407;
  int32_t _M0L7_2abindS408;
  int32_t _M0L7_2abindS409;
  int32_t _M0L1iS410;
  #line 251 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  #line 254 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L3tmpS406
  = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__get__cli__args__ffi();
  _M0L6_2atmpS869 = Moonbit_array_length(_M0L3tmpS406);
  #line 255 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L3resS407 = _M0MPC15array5Array11new_2einnerGsE(_M0L6_2atmpS869);
  _M0L7_2abindS408 = 0;
  _M0L7_2abindS409 = Moonbit_array_length(_M0L3tmpS406);
  _M0L1iS410 = _M0L7_2abindS408;
  while (1) {
    if (_M0L1iS410 < _M0L7_2abindS409) {
      moonbit_bytes_t _M0L6_2atmpS867;
      moonbit_string_t _M0L6_2atmpS866;
      int32_t _M0L6_2atmpS868;
      if (_M0L1iS410 < 0 || _M0L1iS410 >= Moonbit_array_length(_M0L3tmpS406)) {
        #line 257 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS867 = (moonbit_bytes_t)_M0L3tmpS406[_M0L1iS410];
      moonbit_incref(_M0L6_2atmpS867);
      #line 257 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0L6_2atmpS866
      = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398(_M0L61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398, _M0L6_2atmpS867);
      moonbit_decref(_M0L6_2atmpS867);
      #line 257 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0MPC15array5Array4pushGsE(_M0L3resS407, _M0L6_2atmpS866);
      moonbit_decref(_M0L6_2atmpS866);
      _M0L6_2atmpS868 = _M0L1iS410 + 1;
      _M0L1iS410 = _M0L6_2atmpS868;
      continue;
    } else {
      moonbit_decref(_M0L3tmpS406);
    }
    break;
  }
  return _M0L3resS407;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN61moonbit__test__driver__internal__utf8__bytes__to__mbt__stringS398(
  int32_t _M0L6_2aenvS780,
  moonbit_bytes_t _M0L5bytesS399
) {
  struct _M0TPB13StringBuilder* _M0L3resS400;
  int32_t _M0L3lenS401;
  struct _M0TPB8MutLocalGiE* _M0L1iS402;
  moonbit_string_t _result_947;
  #line 207 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  #line 210 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L3resS400 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L3lenS401 = Moonbit_array_length(_M0L5bytesS399);
  _M0L1iS402
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L1iS402)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L1iS402->$0 = 0;
  while (1) {
    int32_t _M0L3valS781 = _M0L1iS402->$0;
    if (_M0L3valS781 < _M0L3lenS401) {
      int32_t _M0L3valS865 = _M0L1iS402->$0;
      int32_t _M0L6_2atmpS864;
      int32_t _M0L6_2atmpS863;
      struct _M0TPB8MutLocalGiE* _M0L1cS403;
      int32_t _M0L3valS782;
      if (
        _M0L3valS865 < 0
        || _M0L3valS865 >= Moonbit_array_length(_M0L5bytesS399)
      ) {
        #line 214 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS864 = _M0L5bytesS399[_M0L3valS865];
      _M0L6_2atmpS863 = (int32_t)_M0L6_2atmpS864;
      _M0L1cS403
      = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
      Moonbit_object_header(_M0L1cS403)->meta
      = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
      _M0L1cS403->$0 = _M0L6_2atmpS863;
      _M0L3valS782 = _M0L1cS403->$0;
      if (_M0L3valS782 < 128) {
        int32_t _M0L3valS784 = _M0L1cS403->$0;
        int32_t _M0L6_2atmpS783;
        int32_t _M0L3valS786;
        int32_t _M0L6_2atmpS785;
        moonbit_decref(_M0L1cS403);
        _M0L6_2atmpS783 = _M0L3valS784;
        #line 216 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS400, _M0L6_2atmpS783);
        _M0L3valS786 = _M0L1iS402->$0;
        _M0L6_2atmpS785 = _M0L3valS786 + 1;
        _M0L1iS402->$0 = _M0L6_2atmpS785;
      } else {
        int32_t _M0L3valS787 = _M0L1cS403->$0;
        if (_M0L3valS787 < 224) {
          int32_t _M0L3valS789 = _M0L1iS402->$0;
          int32_t _M0L6_2atmpS788 = _M0L3valS789 + 1;
          int32_t _M0L3valS798;
          int32_t _M0L6_2atmpS797;
          int32_t _M0L6_2atmpS791;
          int32_t _M0L3valS796;
          int32_t _M0L6_2atmpS795;
          int32_t _M0L6_2atmpS794;
          int32_t _M0L6_2atmpS793;
          int32_t _M0L6_2atmpS792;
          int32_t _M0L6_2atmpS790;
          int32_t _M0L3valS800;
          int32_t _M0L6_2atmpS799;
          int32_t _M0L3valS802;
          int32_t _M0L6_2atmpS801;
          if (_M0L6_2atmpS788 >= _M0L3lenS401) {
            moonbit_decref(_M0L1cS403);
            moonbit_decref(_M0L1iS402);
            break;
          }
          _M0L3valS798 = _M0L1cS403->$0;
          _M0L6_2atmpS797 = _M0L3valS798 & 31;
          _M0L6_2atmpS791 = _M0L6_2atmpS797 << 6;
          _M0L3valS796 = _M0L1iS402->$0;
          _M0L6_2atmpS795 = _M0L3valS796 + 1;
          if (
            _M0L6_2atmpS795 < 0
            || _M0L6_2atmpS795 >= Moonbit_array_length(_M0L5bytesS399)
          ) {
            #line 222 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
            moonbit_panic();
          }
          _M0L6_2atmpS794 = _M0L5bytesS399[_M0L6_2atmpS795];
          _M0L6_2atmpS793 = (int32_t)_M0L6_2atmpS794;
          _M0L6_2atmpS792 = _M0L6_2atmpS793 & 63;
          _M0L6_2atmpS790 = _M0L6_2atmpS791 | _M0L6_2atmpS792;
          _M0L1cS403->$0 = _M0L6_2atmpS790;
          _M0L3valS800 = _M0L1cS403->$0;
          moonbit_decref(_M0L1cS403);
          _M0L6_2atmpS799 = _M0L3valS800;
          #line 223 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
          _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS400, _M0L6_2atmpS799);
          _M0L3valS802 = _M0L1iS402->$0;
          _M0L6_2atmpS801 = _M0L3valS802 + 2;
          _M0L1iS402->$0 = _M0L6_2atmpS801;
        } else {
          int32_t _M0L3valS803 = _M0L1cS403->$0;
          if (_M0L3valS803 < 240) {
            int32_t _M0L3valS805 = _M0L1iS402->$0;
            int32_t _M0L6_2atmpS804 = _M0L3valS805 + 2;
            int32_t _M0L3valS821;
            int32_t _M0L6_2atmpS820;
            int32_t _M0L6_2atmpS813;
            int32_t _M0L3valS819;
            int32_t _M0L6_2atmpS818;
            int32_t _M0L6_2atmpS817;
            int32_t _M0L6_2atmpS816;
            int32_t _M0L6_2atmpS815;
            int32_t _M0L6_2atmpS814;
            int32_t _M0L6_2atmpS807;
            int32_t _M0L3valS812;
            int32_t _M0L6_2atmpS811;
            int32_t _M0L6_2atmpS810;
            int32_t _M0L6_2atmpS809;
            int32_t _M0L6_2atmpS808;
            int32_t _M0L6_2atmpS806;
            int32_t _M0L3valS823;
            int32_t _M0L6_2atmpS822;
            int32_t _M0L3valS825;
            int32_t _M0L6_2atmpS824;
            if (_M0L6_2atmpS804 >= _M0L3lenS401) {
              moonbit_decref(_M0L1cS403);
              moonbit_decref(_M0L1iS402);
              break;
            }
            _M0L3valS821 = _M0L1cS403->$0;
            _M0L6_2atmpS820 = _M0L3valS821 & 15;
            _M0L6_2atmpS813 = _M0L6_2atmpS820 << 12;
            _M0L3valS819 = _M0L1iS402->$0;
            _M0L6_2atmpS818 = _M0L3valS819 + 1;
            if (
              _M0L6_2atmpS818 < 0
              || _M0L6_2atmpS818 >= Moonbit_array_length(_M0L5bytesS399)
            ) {
              #line 230 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS817 = _M0L5bytesS399[_M0L6_2atmpS818];
            _M0L6_2atmpS816 = (int32_t)_M0L6_2atmpS817;
            _M0L6_2atmpS815 = _M0L6_2atmpS816 & 63;
            _M0L6_2atmpS814 = _M0L6_2atmpS815 << 6;
            _M0L6_2atmpS807 = _M0L6_2atmpS813 | _M0L6_2atmpS814;
            _M0L3valS812 = _M0L1iS402->$0;
            _M0L6_2atmpS811 = _M0L3valS812 + 2;
            if (
              _M0L6_2atmpS811 < 0
              || _M0L6_2atmpS811 >= Moonbit_array_length(_M0L5bytesS399)
            ) {
              #line 231 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS810 = _M0L5bytesS399[_M0L6_2atmpS811];
            _M0L6_2atmpS809 = (int32_t)_M0L6_2atmpS810;
            _M0L6_2atmpS808 = _M0L6_2atmpS809 & 63;
            _M0L6_2atmpS806 = _M0L6_2atmpS807 | _M0L6_2atmpS808;
            _M0L1cS403->$0 = _M0L6_2atmpS806;
            _M0L3valS823 = _M0L1cS403->$0;
            moonbit_decref(_M0L1cS403);
            _M0L6_2atmpS822 = _M0L3valS823;
            #line 232 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS400, _M0L6_2atmpS822);
            _M0L3valS825 = _M0L1iS402->$0;
            _M0L6_2atmpS824 = _M0L3valS825 + 3;
            _M0L1iS402->$0 = _M0L6_2atmpS824;
          } else {
            int32_t _M0L3valS827 = _M0L1iS402->$0;
            int32_t _M0L6_2atmpS826 = _M0L3valS827 + 3;
            int32_t _M0L3valS850;
            int32_t _M0L6_2atmpS849;
            int32_t _M0L6_2atmpS842;
            int32_t _M0L3valS848;
            int32_t _M0L6_2atmpS847;
            int32_t _M0L6_2atmpS846;
            int32_t _M0L6_2atmpS845;
            int32_t _M0L6_2atmpS844;
            int32_t _M0L6_2atmpS843;
            int32_t _M0L6_2atmpS835;
            int32_t _M0L3valS841;
            int32_t _M0L6_2atmpS840;
            int32_t _M0L6_2atmpS839;
            int32_t _M0L6_2atmpS838;
            int32_t _M0L6_2atmpS837;
            int32_t _M0L6_2atmpS836;
            int32_t _M0L6_2atmpS829;
            int32_t _M0L3valS834;
            int32_t _M0L6_2atmpS833;
            int32_t _M0L6_2atmpS832;
            int32_t _M0L6_2atmpS831;
            int32_t _M0L6_2atmpS830;
            int32_t _M0L6_2atmpS828;
            int32_t _M0L3valS852;
            int32_t _M0L6_2atmpS851;
            int32_t _M0L3valS856;
            int32_t _M0L6_2atmpS855;
            int32_t _M0L6_2atmpS854;
            int32_t _M0L6_2atmpS853;
            int32_t _M0L3valS860;
            int32_t _M0L6_2atmpS859;
            int32_t _M0L6_2atmpS858;
            int32_t _M0L6_2atmpS857;
            int32_t _M0L3valS862;
            int32_t _M0L6_2atmpS861;
            if (_M0L6_2atmpS826 >= _M0L3lenS401) {
              moonbit_decref(_M0L1cS403);
              moonbit_decref(_M0L1iS402);
              break;
            }
            _M0L3valS850 = _M0L1cS403->$0;
            _M0L6_2atmpS849 = _M0L3valS850 & 7;
            _M0L6_2atmpS842 = _M0L6_2atmpS849 << 18;
            _M0L3valS848 = _M0L1iS402->$0;
            _M0L6_2atmpS847 = _M0L3valS848 + 1;
            if (
              _M0L6_2atmpS847 < 0
              || _M0L6_2atmpS847 >= Moonbit_array_length(_M0L5bytesS399)
            ) {
              #line 239 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS846 = _M0L5bytesS399[_M0L6_2atmpS847];
            _M0L6_2atmpS845 = (int32_t)_M0L6_2atmpS846;
            _M0L6_2atmpS844 = _M0L6_2atmpS845 & 63;
            _M0L6_2atmpS843 = _M0L6_2atmpS844 << 12;
            _M0L6_2atmpS835 = _M0L6_2atmpS842 | _M0L6_2atmpS843;
            _M0L3valS841 = _M0L1iS402->$0;
            _M0L6_2atmpS840 = _M0L3valS841 + 2;
            if (
              _M0L6_2atmpS840 < 0
              || _M0L6_2atmpS840 >= Moonbit_array_length(_M0L5bytesS399)
            ) {
              #line 240 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS839 = _M0L5bytesS399[_M0L6_2atmpS840];
            _M0L6_2atmpS838 = (int32_t)_M0L6_2atmpS839;
            _M0L6_2atmpS837 = _M0L6_2atmpS838 & 63;
            _M0L6_2atmpS836 = _M0L6_2atmpS837 << 6;
            _M0L6_2atmpS829 = _M0L6_2atmpS835 | _M0L6_2atmpS836;
            _M0L3valS834 = _M0L1iS402->$0;
            _M0L6_2atmpS833 = _M0L3valS834 + 3;
            if (
              _M0L6_2atmpS833 < 0
              || _M0L6_2atmpS833 >= Moonbit_array_length(_M0L5bytesS399)
            ) {
              #line 241 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
              moonbit_panic();
            }
            _M0L6_2atmpS832 = _M0L5bytesS399[_M0L6_2atmpS833];
            _M0L6_2atmpS831 = (int32_t)_M0L6_2atmpS832;
            _M0L6_2atmpS830 = _M0L6_2atmpS831 & 63;
            _M0L6_2atmpS828 = _M0L6_2atmpS829 | _M0L6_2atmpS830;
            _M0L1cS403->$0 = _M0L6_2atmpS828;
            _M0L3valS852 = _M0L1cS403->$0;
            _M0L6_2atmpS851 = _M0L3valS852 - 65536;
            _M0L1cS403->$0 = _M0L6_2atmpS851;
            _M0L3valS856 = _M0L1cS403->$0;
            _M0L6_2atmpS855 = _M0L3valS856 >> 10;
            _M0L6_2atmpS854 = _M0L6_2atmpS855 + 55296;
            _M0L6_2atmpS853 = _M0L6_2atmpS854;
            #line 243 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS400, _M0L6_2atmpS853);
            _M0L3valS860 = _M0L1cS403->$0;
            moonbit_decref(_M0L1cS403);
            _M0L6_2atmpS859 = _M0L3valS860 & 1023;
            _M0L6_2atmpS858 = _M0L6_2atmpS859 + 56320;
            _M0L6_2atmpS857 = _M0L6_2atmpS858;
            #line 244 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
            _M0IPB13StringBuilderPB6Logger11write__char(_M0L3resS400, _M0L6_2atmpS857);
            _M0L3valS862 = _M0L1iS402->$0;
            _M0L6_2atmpS861 = _M0L3valS862 + 4;
            _M0L1iS402->$0 = _M0L6_2atmpS861;
          }
        }
      }
      continue;
    } else {
      moonbit_decref(_M0L1iS402);
    }
    break;
  }
  #line 248 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_947 = _M0MPB13StringBuilder10to__string(_M0L3resS400);
  moonbit_decref(_M0L3resS400);
  return _result_947;
}

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__argsN45moonbit__test__driver__internal__parse__int__S391(
  int32_t _M0L6_2aenvS773,
  moonbit_string_t _M0L1sS392
) {
  struct _M0TPB8MutLocalGiE* _M0L3resS393;
  int32_t _M0L3lenS394;
  int32_t _M0L7_2abindS395;
  int32_t _M0L1iS396;
  int32_t _result_949;
  #line 198 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L3resS393
  = (struct _M0TPB8MutLocalGiE*)moonbit_malloc(sizeof(struct _M0TPB8MutLocalGiE));
  Moonbit_object_header(_M0L3resS393)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_SCALAR, 0, 0);
  _M0L3resS393->$0 = 0;
  _M0L3lenS394 = Moonbit_array_length(_M0L1sS392);
  _M0L7_2abindS395 = 0;
  _M0L1iS396 = _M0L7_2abindS395;
  while (1) {
    if (_M0L1iS396 < _M0L3lenS394) {
      int32_t _M0L3valS778 = _M0L3resS393->$0;
      int32_t _M0L6_2atmpS775 = _M0L3valS778 * 10;
      int32_t _M0L6_2atmpS777;
      int32_t _M0L6_2atmpS776;
      int32_t _M0L6_2atmpS774;
      int32_t _M0L6_2atmpS779;
      if (_M0L1iS396 < 0 || _M0L1iS396 >= Moonbit_array_length(_M0L1sS392)) {
        #line 202 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        moonbit_panic();
      }
      _M0L6_2atmpS777 = _M0L1sS392[_M0L1iS396];
      _M0L6_2atmpS776 = _M0L6_2atmpS777 - 48;
      _M0L6_2atmpS774 = _M0L6_2atmpS775 + _M0L6_2atmpS776;
      _M0L3resS393->$0 = _M0L6_2atmpS774;
      _M0L6_2atmpS779 = _M0L1iS396 + 1;
      _M0L1iS396 = _M0L6_2atmpS779;
      continue;
    }
    break;
  }
  _result_949 = _M0L3resS393->$0;
  moonbit_decref(_M0L3resS393);
  return _result_949;
}

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__execute(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS390,
  moonbit_string_t _M0L8filenameS359,
  int32_t _M0L5indexS361
) {
  struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357* _closure_950;
  struct _M0TWEu* _M0L13handle__startS357;
  struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362* _closure_951;
  struct _M0TWssbEu* _M0L14handle__resultS362;
  struct _M0TWRPC15error5ErrorEs* _M0L17error__to__stringS369;
  void* _M0L11_2atry__errS384;
  struct moonbit_result_0 _tmp_953;
  int32_t _handle__error__result_954;
  int32_t _M0L6_2atmpS761;
  void* _M0L3errS385;
  moonbit_string_t _M0L4nameS387;
  struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest* _M0L36_2aMoonBitTestDriverInternalSkipTestS388;
  moonbit_string_t _M0L8_2afieldS900;
  int32_t _M0L6_2acntS940;
  moonbit_string_t _M0L7_2anameS389;
  #line 515 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  moonbit_incref(_M0L8filenameS359);
  _closure_950
  = (struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357*)moonbit_malloc(sizeof(struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357));
  Moonbit_object_header(_closure_950)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 9, 0);
  _closure_950->code
  = &_M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN13handle__startS357;
  _closure_950->$0 = _M0L5indexS361;
  _closure_950->$1 = _M0L8filenameS359;
  _M0L13handle__startS357 = (struct _M0TWEu*)_closure_950;
  moonbit_incref(_M0L8filenameS359);
  _closure_951
  = (struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362*)moonbit_malloc(sizeof(struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362));
  Moonbit_object_header(_closure_951)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 12, 0);
  _closure_951->code
  = &_M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN14handle__resultS362;
  _closure_951->$0 = _M0L5indexS361;
  _closure_951->$1 = _M0L8filenameS359;
  _M0L14handle__resultS362 = (struct _M0TWssbEu*)_closure_951;
  _M0L17error__to__stringS369
  = (struct _M0TWRPC15error5ErrorEs*)&_M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN17error__to__stringS369$closure.data;
  #line 557 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _tmp_953
  = _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main41MoonBit__Test__Driver__Internal__No__ArgsE(_M0L12async__testsS390, _M0L8filenameS359, _M0L5indexS361, _M0L13handle__startS357, _M0L14handle__resultS362, _M0L17error__to__stringS369);
  if (_tmp_953.tag) {
    int32_t const _M0L5_2aokS770 = _tmp_953.data.ok;
    _handle__error__result_954 = _M0L5_2aokS770;
  } else {
    void* const _M0L6_2aerrS771 = _tmp_953.data.err;
    moonbit_decref(_M0L17error__to__stringS369);
    moonbit_decref(_M0L13handle__startS357);
    _M0L11_2atry__errS384 = _M0L6_2aerrS771;
    goto join_383;
  }
  if (_handle__error__result_954) {
    moonbit_decref(_M0L17error__to__stringS369);
    moonbit_decref(_M0L13handle__startS357);
    _M0L6_2atmpS761 = 1;
  } else {
    struct moonbit_result_0 _tmp_955;
    int32_t _handle__error__result_956;
    #line 560 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
    _tmp_955
    = _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main43MoonBit__Test__Driver__Internal__With__ArgsE(_M0L12async__testsS390, _M0L8filenameS359, _M0L5indexS361, _M0L13handle__startS357, _M0L14handle__resultS362, _M0L17error__to__stringS369);
    if (_tmp_955.tag) {
      int32_t const _M0L5_2aokS768 = _tmp_955.data.ok;
      _handle__error__result_956 = _M0L5_2aokS768;
    } else {
      void* const _M0L6_2aerrS769 = _tmp_955.data.err;
      moonbit_decref(_M0L17error__to__stringS369);
      moonbit_decref(_M0L13handle__startS357);
      _M0L11_2atry__errS384 = _M0L6_2aerrS769;
      goto join_383;
    }
    if (_handle__error__result_956) {
      moonbit_decref(_M0L17error__to__stringS369);
      moonbit_decref(_M0L13handle__startS357);
      _M0L6_2atmpS761 = 1;
    } else {
      struct moonbit_result_0 _tmp_957;
      int32_t _handle__error__result_958;
      #line 563 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _tmp_957
      = _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main48MoonBit__Test__Driver__Internal__Async__No__ArgsE(_M0L12async__testsS390, _M0L8filenameS359, _M0L5indexS361, _M0L13handle__startS357, _M0L14handle__resultS362, _M0L17error__to__stringS369);
      if (_tmp_957.tag) {
        int32_t const _M0L5_2aokS766 = _tmp_957.data.ok;
        _handle__error__result_958 = _M0L5_2aokS766;
      } else {
        void* const _M0L6_2aerrS767 = _tmp_957.data.err;
        moonbit_decref(_M0L17error__to__stringS369);
        moonbit_decref(_M0L13handle__startS357);
        _M0L11_2atry__errS384 = _M0L6_2aerrS767;
        goto join_383;
      }
      if (_handle__error__result_958) {
        moonbit_decref(_M0L17error__to__stringS369);
        moonbit_decref(_M0L13handle__startS357);
        _M0L6_2atmpS761 = 1;
      } else {
        struct moonbit_result_0 _tmp_959;
        int32_t _handle__error__result_960;
        #line 566 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
        _tmp_959
        = _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main50MoonBit__Test__Driver__Internal__Async__With__ArgsE(_M0L12async__testsS390, _M0L8filenameS359, _M0L5indexS361, _M0L13handle__startS357, _M0L14handle__resultS362, _M0L17error__to__stringS369);
        if (_tmp_959.tag) {
          int32_t const _M0L5_2aokS764 = _tmp_959.data.ok;
          _handle__error__result_960 = _M0L5_2aokS764;
        } else {
          void* const _M0L6_2aerrS765 = _tmp_959.data.err;
          moonbit_decref(_M0L17error__to__stringS369);
          moonbit_decref(_M0L13handle__startS357);
          _M0L11_2atry__errS384 = _M0L6_2aerrS765;
          goto join_383;
        }
        if (_handle__error__result_960) {
          moonbit_decref(_M0L17error__to__stringS369);
          moonbit_decref(_M0L13handle__startS357);
          _M0L6_2atmpS761 = 1;
        } else {
          struct moonbit_result_0 _tmp_961;
          #line 569 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
          _tmp_961
          = _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(_M0L12async__testsS390, _M0L8filenameS359, _M0L5indexS361, _M0L13handle__startS357, _M0L14handle__resultS362, _M0L17error__to__stringS369);
          moonbit_decref(_M0L13handle__startS357);
          moonbit_decref(_M0L17error__to__stringS369);
          if (_tmp_961.tag) {
            int32_t const _M0L5_2aokS762 = _tmp_961.data.ok;
            _M0L6_2atmpS761 = _M0L5_2aokS762;
          } else {
            void* const _M0L6_2aerrS763 = _tmp_961.data.err;
            _M0L11_2atry__errS384 = _M0L6_2aerrS763;
            goto join_383;
          }
        }
      }
    }
  }
  if (!_M0L6_2atmpS761) {
    void* _M0L109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS772 =
      (void*)moonbit_malloc(sizeof(struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest));
    Moonbit_object_header(_M0L109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS772)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 15, 1);
    ((struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS772)->$0
    = (moonbit_string_t)moonbit_string_literal_0.data;
    _M0L11_2atry__errS384
    = _M0L109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTestS772;
    goto join_383;
  } else {
    moonbit_decref(_M0L14handle__resultS362);
  }
  goto joinlet_952;
  join_383:;
  _M0L3errS385 = _M0L11_2atry__errS384;
  _M0L36_2aMoonBitTestDriverInternalSkipTestS388
  = (struct _M0DTPC15error5Error109biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalSkipTest_2eMoonBitTestDriverInternalSkipTest*)_M0L3errS385;
  _M0L8_2afieldS900 = _M0L36_2aMoonBitTestDriverInternalSkipTestS388->$0;
  _M0L6_2acntS940
  = Moonbit_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS388));
  if (_M0L6_2acntS940 > 1) {
    int32_t _M0L11_2anew__cntS941 = _M0L6_2acntS940 - 1;
    Moonbit_set_rc_count(Moonbit_object_header(_M0L36_2aMoonBitTestDriverInternalSkipTestS388), _M0L11_2anew__cntS941);
    moonbit_incref(_M0L8_2afieldS900);
  } else if (_M0L6_2acntS940 == 1) {
    #line 576 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
    moonbit_free(_M0L36_2aMoonBitTestDriverInternalSkipTestS388);
  }
  _M0L7_2anameS389 = _M0L8_2afieldS900;
  _M0L4nameS387 = _M0L7_2anameS389;
  goto join_386;
  goto joinlet_962;
  join_386:;
  #line 577 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN14handle__resultS362(_M0L14handle__resultS362, _M0L4nameS387, (moonbit_string_t)moonbit_string_literal_1.data, 1);
  moonbit_decref(_M0L14handle__resultS362);
  moonbit_decref(_M0L4nameS387);
  joinlet_962:;
  joinlet_952:;
  return 0;
}

moonbit_string_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN17error__to__stringS369(
  struct _M0TWRPC15error5ErrorEs* _M0L6_2aenvS760,
  void* _M0L3errS370
) {
  void* _M0L1eS372;
  moonbit_string_t _M0L1eS374;
  moonbit_string_t _result_965;
  #line 546 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  switch (Moonbit_object_tag(_M0L3errS370)) {
    case 0: {
      struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS375 =
        (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L3errS370;
      moonbit_string_t _M0L4_2aeS376 = _M0L10_2aFailureS375->$0;
      moonbit_incref(_M0L4_2aeS376);
      _M0L1eS374 = _M0L4_2aeS376;
      goto join_373;
      break;
    }
    
    case 2: {
      struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError* _M0L15_2aInspectErrorS377 =
        (struct _M0DTPC15error5Error58moonbitlang_2fcore_2fbuiltin_2eInspectError_2eInspectError*)_M0L3errS370;
      moonbit_string_t _M0L4_2aeS378 = _M0L15_2aInspectErrorS377->$0;
      moonbit_incref(_M0L4_2aeS378);
      _M0L1eS374 = _M0L4_2aeS378;
      goto join_373;
      break;
    }
    
    case 3: {
      struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError* _M0L16_2aSnapshotErrorS379 =
        (struct _M0DTPC15error5Error60moonbitlang_2fcore_2fbuiltin_2eSnapshotError_2eSnapshotError*)_M0L3errS370;
      moonbit_string_t _M0L4_2aeS380 = _M0L16_2aSnapshotErrorS379->$0;
      moonbit_incref(_M0L4_2aeS380);
      _M0L1eS374 = _M0L4_2aeS380;
      goto join_373;
      break;
    }
    
    case 4: {
      struct _M0DTPC15error5Error107biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError* _M0L35_2aMoonBitTestDriverInternalJsErrorS381 =
        (struct _M0DTPC15error5Error107biolab_2fbio__seq_2fcmd_2fseqio__main_2eMoonBitTestDriverInternalJsError_2eMoonBitTestDriverInternalJsError*)_M0L3errS370;
      moonbit_string_t _M0L4_2aeS382 =
        _M0L35_2aMoonBitTestDriverInternalJsErrorS381->$0;
      moonbit_incref(_M0L4_2aeS382);
      _M0L1eS374 = _M0L4_2aeS382;
      goto join_373;
      break;
    }
    default: {
      moonbit_incref(_M0L3errS370);
      _M0L1eS372 = _M0L3errS370;
      goto join_371;
      break;
    }
  }
  join_373:;
  return _M0L1eS374;
  join_371:;
  #line 552 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_965 = _M0FP15Error10to__string(_M0L1eS372);
  moonbit_decref(_M0L1eS372);
  return _result_965;
}

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN14handle__resultS362(
  struct _M0TWssbEu* _M0L6_2aenvS757,
  moonbit_string_t _M0L10__testnameS363,
  moonbit_string_t _M0L7messageS364,
  int32_t _M0L7skippedS365
) {
  struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362* _M0L14_2acasted__envS758;
  moonbit_string_t _M0L8filenameS359;
  int32_t _M0L5indexS361;
  int32_t _if__result_966;
  moonbit_string_t _M0L10file__nameS366;
  moonbit_string_t _M0L7messageS367;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS368;
  moonbit_string_t _M0L6_2atmpS759;
  #line 531 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L14_2acasted__envS758
  = (struct _M0R110_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__result_7c362*)_M0L6_2aenvS757;
  _M0L8filenameS359 = _M0L14_2acasted__envS758->$1;
  _M0L5indexS361 = _M0L14_2acasted__envS758->$0;
  if (!_M0L7skippedS365) {
    _if__result_966 = 1;
  } else {
    _if__result_966 = 0;
  }
  if (_if__result_966) {
    
  }
  #line 537 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L10file__nameS366
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS359, 1);
  #line 538 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L7messageS367
  = _M0MPC16string6String14escape_2einner(_M0L7messageS364, 1);
  #line 539 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L18_2astring__builderS368
  = _M0MPB13StringBuilder21StringBuilder_2einner(45);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS368, (moonbit_string_t)moonbit_string_literal_3.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS368, _M0L10file__nameS366);
  moonbit_decref(_M0L10file__nameS366);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS368, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS368, _M0L5indexS361);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS368, (moonbit_string_t)moonbit_string_literal_5.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS368, _M0L7messageS367);
  moonbit_decref(_M0L7messageS367);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS368, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 541 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS759
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS368);
  moonbit_decref(_M0L18_2astring__builderS368);
  #line 540 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS759);
  moonbit_decref(_M0L6_2atmpS759);
  #line 543 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

int32_t _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__executeN13handle__startS357(
  struct _M0TWEu* _M0L6_2aenvS754
) {
  struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357* _M0L14_2acasted__envS755;
  moonbit_string_t _M0L8filenameS359;
  int32_t _M0L5indexS361;
  moonbit_string_t _M0L10file__nameS358;
  struct _M0TPB13StringBuilder* _M0L18_2astring__builderS360;
  moonbit_string_t _M0L6_2atmpS756;
  #line 522 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L14_2acasted__envS755
  = (struct _M0R109_24biolab_2fbio__seq_2fcmd_2fseqio__main_2emoonbit__test__driver__internal__do__execute_2ehandle__start_7c357*)_M0L6_2aenvS754;
  _M0L8filenameS359 = _M0L14_2acasted__envS755->$1;
  _M0L5indexS361 = _M0L14_2acasted__envS755->$0;
  #line 523 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L10file__nameS358
  = _M0MPC16string6String14escape_2einner(_M0L8filenameS359, 1);
  #line 524 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_2.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L18_2astring__builderS360
  = _M0MPB13StringBuilder21StringBuilder_2einner(33);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS360, (moonbit_string_t)moonbit_string_literal_8.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGsE(_M0L18_2astring__builderS360, _M0L10file__nameS358);
  moonbit_decref(_M0L10file__nameS358);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS360, (moonbit_string_t)moonbit_string_literal_4.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS360, _M0L5indexS361);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS360, (moonbit_string_t)moonbit_string_literal_6.data);
  #line 526 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L6_2atmpS756
  = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS360);
  moonbit_decref(_M0L18_2astring__builderS360);
  #line 525 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE(_M0L6_2atmpS756);
  moonbit_decref(_M0L6_2atmpS756);
  #line 528 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0FPB7printlnGsE((moonbit_string_t)moonbit_string_literal_7.data);
  return 0;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main41MoonBit__Test__Driver__Internal__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S327,
  moonbit_string_t _M0L12_2adiscard__S328,
  int32_t _M0L12_2adiscard__S329,
  struct _M0TWEu* _M0L12_2adiscard__S330,
  struct _M0TWssbEu* _M0L12_2adiscard__S331,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S332
) {
  struct moonbit_result_0 _result_967;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_967.tag = 1;
  _result_967.data.ok = 0;
  return _result_967;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main43MoonBit__Test__Driver__Internal__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S333,
  moonbit_string_t _M0L12_2adiscard__S334,
  int32_t _M0L12_2adiscard__S335,
  struct _M0TWEu* _M0L12_2adiscard__S336,
  struct _M0TWssbEu* _M0L12_2adiscard__S337,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S338
) {
  struct moonbit_result_0 _result_968;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_968.tag = 1;
  _result_968.data.ok = 0;
  return _result_968;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main48MoonBit__Test__Driver__Internal__Async__No__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S339,
  moonbit_string_t _M0L12_2adiscard__S340,
  int32_t _M0L12_2adiscard__S341,
  struct _M0TWEu* _M0L12_2adiscard__S342,
  struct _M0TWssbEu* _M0L12_2adiscard__S343,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S344
) {
  struct moonbit_result_0 _result_969;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_969.tag = 1;
  _result_969.data.ok = 0;
  return _result_969;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main50MoonBit__Test__Driver__Internal__Async__With__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S345,
  moonbit_string_t _M0L12_2adiscard__S346,
  int32_t _M0L12_2adiscard__S347,
  struct _M0TWEu* _M0L12_2adiscard__S348,
  struct _M0TWssbEu* _M0L12_2adiscard__S349,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S350
) {
  struct moonbit_result_0 _result_970;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_970.tag = 1;
  _result_970.data.ok = 0;
  return _result_970;
}

struct moonbit_result_0 _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main21MoonBit__Test__Driver9run__testGRP46biolab8bio__seq3cmd11seqio__main50MoonBit__Test__Driver__Internal__With__Bench__ArgsE(
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12_2adiscard__S351,
  moonbit_string_t _M0L12_2adiscard__S352,
  int32_t _M0L12_2adiscard__S353,
  struct _M0TWEu* _M0L12_2adiscard__S354,
  struct _M0TWssbEu* _M0L12_2adiscard__S355,
  struct _M0TWRPC15error5ErrorEs* _M0L12_2adiscard__S356
) {
  struct moonbit_result_0 _result_971;
  #line 35 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _result_971.tag = 1;
  _result_971.data.ok = 0;
  return _result_971;
}

moonbit_string_t _M0MPC15array5Array2atGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS325,
  int32_t _M0L5indexS326
) {
  int32_t _M0L3lenS324;
  int32_t _if__result_972;
  #line 181 "/home/developer/.moon/lib/core/builtin/array.mbt"
  _M0L3lenS324 = _M0L4selfS325->$1;
  if (_M0L5indexS326 >= 0) {
    _if__result_972 = _M0L5indexS326 < _M0L3lenS324;
  } else {
    _if__result_972 = 0;
  }
  if (_if__result_972) {
    moonbit_string_t* _M0L6_2atmpS753;
    moonbit_string_t _M0L6_2atmpS901;
    #line 186 "/home/developer/.moon/lib/core/builtin/array.mbt"
    _M0L6_2atmpS753 = _M0MPC15array5Array6bufferGsE(_M0L4selfS325);
    _M0L6_2atmpS901 = (moonbit_string_t)_M0L6_2atmpS753[_M0L5indexS326];
    moonbit_incref(_M0L6_2atmpS901);
    moonbit_decref(_M0L6_2atmpS753);
    return _M0L6_2atmpS901;
  } else {
    #line 185 "/home/developer/.moon/lib/core/builtin/array.mbt"
    moonbit_panic();
  }
}

int32_t _M0FPB7printlnGsE(moonbit_string_t _M0L5inputS323) {
  moonbit_string_t _M0L6_2atmpS752;
  #line 36 "/home/developer/.moon/lib/core/builtin/console.mbt"
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  _M0L6_2atmpS752 = _M0IPC16string6StringPB4Show10to__string(_M0L5inputS323);
  #line 37 "/home/developer/.moon/lib/core/builtin/console.mbt"
  moonbit_println(_M0L6_2atmpS752);
  moonbit_decref(_M0L6_2atmpS752);
  return 0;
}

moonbit_string_t _M0IPC13int3IntPB4Show10to__string(int32_t _M0L4selfS322) {
  #line 35 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 36 "/home/developer/.moon/lib/core/builtin/show.mbt"
  return _M0MPC13int3Int18to__string_2einner(_M0L4selfS322, 10);
}

int32_t _M0MPC15array5Array4pushGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS316,
  moonbit_string_t _M0L5valueS318
) {
  int32_t _M0L3lenS742;
  moonbit_string_t* _M0L6_2atmpS744;
  int32_t _M0L6_2atmpS743;
  int32_t _M0L6lengthS317;
  moonbit_string_t* _M0L3bufS745;
  moonbit_string_t _M0L6_2aoldS902;
  int32_t _M0L6_2atmpS746;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS742 = _M0L4selfS316->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS744 = _M0MPC15array5Array6bufferGsE(_M0L4selfS316);
  _M0L6_2atmpS743 = Moonbit_array_length(_M0L6_2atmpS744);
  moonbit_decref(_M0L6_2atmpS744);
  if (_M0L3lenS742 == _M0L6_2atmpS743) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGsE(_M0L4selfS316);
  }
  _M0L6lengthS317 = _M0L4selfS316->$1;
  _M0L3bufS745 = _M0L4selfS316->$0;
  _M0L6_2aoldS902 = (moonbit_string_t)_M0L3bufS745[_M0L6lengthS317];
  moonbit_incref(_M0L5valueS318);
  moonbit_decref(_M0L6_2aoldS902);
  _M0L3bufS745[_M0L6lengthS317] = _M0L5valueS318;
  _M0L6_2atmpS746 = _M0L6lengthS317 + 1;
  _M0L4selfS316->$1 = _M0L6_2atmpS746;
  return 0;
}

int32_t _M0MPC15array5Array4pushGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS319,
  struct _M0TUsiE* _M0L5valueS321
) {
  int32_t _M0L3lenS747;
  struct _M0TUsiE** _M0L6_2atmpS749;
  int32_t _M0L6_2atmpS748;
  int32_t _M0L6lengthS320;
  struct _M0TUsiE** _M0L3bufS750;
  struct _M0TUsiE* _M0L6_2aoldS904;
  int32_t _M0L6_2atmpS751;
  #line 305 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L3lenS747 = _M0L4selfS319->$1;
  #line 306 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L6_2atmpS749 = _M0MPC15array5Array6bufferGUsiEE(_M0L4selfS319);
  _M0L6_2atmpS748 = Moonbit_array_length(_M0L6_2atmpS749);
  moonbit_decref(_M0L6_2atmpS749);
  if (_M0L3lenS747 == _M0L6_2atmpS748) {
    #line 307 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
    _M0MPC15array5Array7reallocGUsiEE(_M0L4selfS319);
  }
  _M0L6lengthS320 = _M0L4selfS319->$1;
  _M0L3bufS750 = _M0L4selfS319->$0;
  _M0L6_2aoldS904 = (struct _M0TUsiE*)_M0L3bufS750[_M0L6lengthS320];
  moonbit_incref(_M0L5valueS321);
  if (_M0L6_2aoldS904) {
    moonbit_decref(_M0L6_2aoldS904);
  }
  _M0L3bufS750[_M0L6lengthS320] = _M0L5valueS321;
  _M0L6_2atmpS751 = _M0L6lengthS320 + 1;
  _M0L4selfS319->$1 = _M0L6_2atmpS751;
  return 0;
}

int32_t _M0MPC15array5Array7reallocGsE(struct _M0TPB5ArrayGsE* _M0L4selfS311) {
  int32_t _M0L8old__capS310;
  int32_t _M0L8new__capS312;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS310 = _M0L4selfS311->$1;
  if (_M0L8old__capS310 == 0) {
    _M0L8new__capS312 = 8;
  } else {
    _M0L8new__capS312 = _M0L8old__capS310 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGsE(_M0L4selfS311, _M0L8new__capS312);
  return 0;
}

int32_t _M0MPC15array5Array7reallocGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS314
) {
  int32_t _M0L8old__capS313;
  int32_t _M0L8new__capS315;
  #line 245 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__capS313 = _M0L4selfS314->$1;
  if (_M0L8old__capS313 == 0) {
    _M0L8new__capS315 = 8;
  } else {
    _M0L8new__capS315 = _M0L8old__capS313 * 2;
  }
  #line 248 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0MPC15array5Array14resize__bufferGUsiEE(_M0L4selfS314, _M0L8new__capS315);
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS299,
  int32_t _M0L13new__capacityS302
) {
  moonbit_string_t* _M0L8old__bufS298;
  int32_t _M0L8old__capS300;
  int32_t _M0L9copy__lenS301;
  moonbit_string_t* _M0L8new__bufS303;
  moonbit_string_t* _M0L6_2aoldS906;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS298 = _M0L4selfS299->$0;
  _M0L8old__capS300 = Moonbit_array_length(_M0L8old__bufS298);
  if (_M0L8old__capS300 < _M0L13new__capacityS302) {
    _M0L9copy__lenS301 = _M0L8old__capS300;
  } else {
    _M0L9copy__lenS301 = _M0L13new__capacityS302;
  }
  moonbit_incref(_M0L8old__bufS298);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS303
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(_M0L8old__bufS298, _M0L13new__capacityS302, _M0L9copy__lenS301, 0, 0);
  moonbit_decref(_M0L8old__bufS298);
  _M0L6_2aoldS906 = _M0L4selfS299->$0;
  moonbit_decref(_M0L6_2aoldS906);
  _M0L4selfS299->$0 = _M0L8new__bufS303;
  return 0;
}

int32_t _M0MPC15array5Array14resize__bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS305,
  int32_t _M0L13new__capacityS308
) {
  struct _M0TUsiE** _M0L8old__bufS304;
  int32_t _M0L8old__capS306;
  int32_t _M0L9copy__lenS307;
  struct _M0TUsiE** _M0L8new__bufS309;
  struct _M0TUsiE** _M0L6_2aoldS908;
  #line 189 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8old__bufS304 = _M0L4selfS305->$0;
  _M0L8old__capS306 = Moonbit_array_length(_M0L8old__bufS304);
  if (_M0L8old__capS306 < _M0L13new__capacityS308) {
    _M0L9copy__lenS307 = _M0L8old__capS306;
  } else {
    _M0L9copy__lenS307 = _M0L13new__capacityS308;
  }
  moonbit_incref(_M0L8old__bufS304);
  #line 193 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8new__bufS309
  = _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(_M0L8old__bufS304, _M0L13new__capacityS308, _M0L9copy__lenS307, 0, 0);
  moonbit_decref(_M0L8old__bufS304);
  _M0L6_2aoldS908 = _M0L4selfS305->$0;
  moonbit_decref(_M0L6_2aoldS908);
  _M0L4selfS305->$0 = _M0L8new__bufS309;
  return 0;
}

struct _M0TPB5ArrayGsE* _M0MPC15array5Array11new_2einnerGsE(
  int32_t _M0L8capacityS297
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  if (_M0L8capacityS297 == 0) {
    moonbit_string_t* _M0L6_2atmpS740 =
      (moonbit_string_t*)moonbit_empty_ref_array;
    struct _M0TPB5ArrayGsE* _block_973 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_973)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 0);
    _block_973->$0 = _M0L6_2atmpS740;
    _block_973->$1 = 0;
    return _block_973;
  } else {
    moonbit_string_t* _M0L6_2atmpS741 =
      (moonbit_string_t*)moonbit_make_ref_array(_M0L8capacityS297, (moonbit_string_t)moonbit_string_literal_0.data);
    struct _M0TPB5ArrayGsE* _block_974 =
      (struct _M0TPB5ArrayGsE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGsE));
    Moonbit_object_header(_block_974)->meta
    = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 6, 0);
    _block_974->$0 = _M0L6_2atmpS741;
    _block_974->$1 = 0;
    return _block_974;
  }
}

moonbit_string_t _M0IPC16string6StringPB4Show10to__string(
  moonbit_string_t _M0L4selfS296
) {
  #line 222 "/home/developer/.moon/lib/core/builtin/show.mbt"
  moonbit_incref(_M0L4selfS296);
  return _M0L4selfS296;
}

int32_t _M0IPB13StringBuilderPB6Logger11write__view(
  struct _M0TPB13StringBuilder* _M0L4selfS295,
  struct _M0TPC16string10StringView _M0L3strS294
) {
  int32_t _M0L3endS738;
  int32_t _M0L5startS739;
  int32_t _M0L8str__lenS293;
  int32_t _M0L3lenS731;
  int32_t _M0L6_2atmpS730;
  uint16_t* _M0L4dataS732;
  int32_t _M0L3lenS733;
  moonbit_string_t _M0L6_2atmpS734;
  int32_t _M0L6_2atmpS735;
  int32_t _M0L3lenS737;
  int32_t _M0L6_2atmpS736;
  #line 131 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3endS738 = _M0L3strS294.$2;
  _M0L5startS739 = _M0L3strS294.$1;
  _M0L8str__lenS293 = _M0L3endS738 - _M0L5startS739;
  _M0L3lenS731 = _M0L4selfS295->$1;
  _M0L6_2atmpS730 = _M0L3lenS731 + _M0L8str__lenS293;
  #line 136 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS295, _M0L6_2atmpS730);
  _M0L4dataS732 = _M0L4selfS295->$0;
  _M0L3lenS733 = _M0L4selfS295->$1;
  moonbit_incref(_M0L4dataS732);
  #line 139 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS734 = _M0MPC16string10StringView4data(_M0L3strS294);
  #line 140 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS735 = _M0MPC16string10StringView13start__offset(_M0L3strS294);
  #line 137 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS732, _M0L3lenS733, _M0L6_2atmpS734, _M0L6_2atmpS735, _M0L8str__lenS293);
  moonbit_decref(_M0L4dataS732);
  moonbit_decref(_M0L6_2atmpS734);
  _M0L3lenS737 = _M0L4selfS295->$1;
  _M0L6_2atmpS736 = _M0L3lenS737 + _M0L8str__lenS293;
  _M0L4selfS295->$1 = _M0L6_2atmpS736;
  return 0;
}

moonbit_string_t _M0MPC16string6String17unsafe__substring(
  moonbit_string_t _M0L3strS290,
  int32_t _M0L5startS288,
  int32_t _M0L3endS289
) {
  int32_t _if__result_975;
  int32_t _M0L3lenS291;
  int32_t _M0L6_2atmpS728;
  int32_t _M0L6_2atmpS729;
  moonbit_bytes_t _M0L5bytesS292;
  moonbit_bytes_t _M0L6_2atmpS727;
  moonbit_string_t _result_976;
  #line 91 "/home/developer/.moon/lib/core/builtin/string.mbt"
  if (_M0L5startS288 == 0) {
    int32_t _M0L6_2atmpS726 = Moonbit_array_length(_M0L3strS290);
    _if__result_975 = _M0L3endS289 == _M0L6_2atmpS726;
  } else {
    _if__result_975 = 0;
  }
  if (_if__result_975) {
    moonbit_incref(_M0L3strS290);
    return _M0L3strS290;
  }
  _M0L3lenS291 = _M0L3endS289 - _M0L5startS288;
  _M0L6_2atmpS728 = _M0L3lenS291 * 2;
  #line 101 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0L6_2atmpS729 = _M0IPC14byte4BytePB7Default7default();
  _M0L5bytesS292
  = (moonbit_bytes_t)moonbit_make_bytes(_M0L6_2atmpS728, _M0L6_2atmpS729);
  #line 102 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _M0MPC15array10FixedArray18blit__from__string(_M0L5bytesS292, 0, _M0L3strS290, _M0L5startS288, _M0L3lenS291);
  _M0L6_2atmpS727 = _M0L5bytesS292;
  #line 103 "/home/developer/.moon/lib/core/builtin/string.mbt"
  _result_976
  = _M0MPC15bytes5Bytes29to__unchecked__string_2einner(_M0L6_2atmpS727, 0, 4294967296ll);
  moonbit_decref(_M0L6_2atmpS727);
  return _result_976;
}

int32_t _M0IPC14byte4BytePB7Default7default() {
  #line 231 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  return 0;
}

moonbit_string_t _M0MPC15bytes5Bytes29to__unchecked__string_2einner(
  moonbit_bytes_t _M0L4selfS283,
  int32_t _M0L6offsetS287,
  int64_t _M0L6lengthS285
) {
  int32_t _M0L3lenS282;
  int32_t _M0L6lengthS284;
  int32_t _if__result_977;
  #line 77 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L3lenS282 = Moonbit_array_length(_M0L4selfS283);
  if (_M0L6lengthS285 == 4294967296ll) {
    _M0L6lengthS284 = _M0L3lenS282 - _M0L6offsetS287;
  } else {
    int64_t _M0L7_2aSomeS286 = _M0L6lengthS285;
    _M0L6lengthS284 = (int32_t)_M0L7_2aSomeS286;
  }
  if (_M0L6offsetS287 >= 0) {
    if (_M0L6lengthS284 >= 0) {
      int32_t _M0L6_2atmpS725 = _M0L6offsetS287 + _M0L6lengthS284;
      _if__result_977 = _M0L6_2atmpS725 <= _M0L3lenS282;
    } else {
      _if__result_977 = 0;
    }
  } else {
    _if__result_977 = 0;
  }
  if (_if__result_977) {
    moonbit_incref(_M0L4selfS283);
    #line 85 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    return _M0FPB19unsafe__sub__string(_M0L4selfS283, _M0L6offsetS287, _M0L6lengthS284);
  } else {
    #line 84 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
    moonbit_panic();
  }
}

int32_t _M0MPC15array10FixedArray18blit__from__string(
  moonbit_bytes_t _M0L4selfS274,
  int32_t _M0L13bytes__offsetS269,
  moonbit_string_t _M0L3strS276,
  int32_t _M0L11str__offsetS272,
  int32_t _M0L6lengthS270
) {
  int32_t _M0L6_2atmpS724;
  int32_t _M0L6_2atmpS723;
  int32_t _M0L2e1S268;
  int32_t _M0L6_2atmpS722;
  int32_t _M0L2e2S271;
  int32_t _M0L4len1S273;
  int32_t _M0L4len2S275;
  int32_t _if__result_978;
  #line 125 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
  _M0L6_2atmpS724 = _M0L6lengthS270 * 2;
  _M0L6_2atmpS723 = _M0L13bytes__offsetS269 + _M0L6_2atmpS724;
  _M0L2e1S268 = _M0L6_2atmpS723 - 1;
  _M0L6_2atmpS722 = _M0L11str__offsetS272 + _M0L6lengthS270;
  _M0L2e2S271 = _M0L6_2atmpS722 - 1;
  _M0L4len1S273 = Moonbit_array_length(_M0L4selfS274);
  _M0L4len2S275 = Moonbit_array_length(_M0L3strS276);
  if (_M0L6lengthS270 >= 0) {
    if (_M0L13bytes__offsetS269 >= 0) {
      if (_M0L2e1S268 < _M0L4len1S273) {
        if (_M0L11str__offsetS272 >= 0) {
          _if__result_978 = _M0L2e2S271 < _M0L4len2S275;
        } else {
          _if__result_978 = 0;
        }
      } else {
        _if__result_978 = 0;
      }
    } else {
      _if__result_978 = 0;
    }
  } else {
    _if__result_978 = 0;
  }
  if (_if__result_978) {
    int32_t _M0L16end__str__offsetS277 =
      _M0L11str__offsetS272 + _M0L6lengthS270;
    int32_t _M0L1iS278 = _M0L11str__offsetS272;
    int32_t _M0L1jS279 = _M0L13bytes__offsetS269;
    while (1) {
      if (_M0L1iS278 < _M0L16end__str__offsetS277) {
        int32_t _M0L6_2atmpS719 = _M0L3strS276[_M0L1iS278];
        int32_t _M0L6_2atmpS718 = (int32_t)_M0L6_2atmpS719;
        uint32_t _M0L1cS280 = *(uint32_t*)&_M0L6_2atmpS718;
        uint32_t _M0L6_2atmpS714 = _M0L1cS280 & 255u;
        int32_t _M0L6_2atmpS713;
        int32_t _M0L6_2atmpS715;
        uint32_t _M0L6_2atmpS717;
        int32_t _M0L6_2atmpS716;
        int32_t _M0L6_2atmpS720;
        int32_t _M0L6_2atmpS721;
        #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS713 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS714);
        if (
          _M0L1jS279 < 0 || _M0L1jS279 >= Moonbit_array_length(_M0L4selfS274)
        ) {
          #line 142 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS274[_M0L1jS279] = _M0L6_2atmpS713;
        _M0L6_2atmpS715 = _M0L1jS279 + 1;
        _M0L6_2atmpS717 = _M0L1cS280 >> 8;
        #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
        _M0L6_2atmpS716 = _M0MPC14uint4UInt8to__byte(_M0L6_2atmpS717);
        if (
          _M0L6_2atmpS715 < 0
          || _M0L6_2atmpS715 >= Moonbit_array_length(_M0L4selfS274)
        ) {
          #line 143 "/home/developer/.moon/lib/core/builtin/bytes.mbt"
          moonbit_panic();
        }
        _M0L4selfS274[_M0L6_2atmpS715] = _M0L6_2atmpS716;
        _M0L6_2atmpS720 = _M0L1iS278 + 1;
        _M0L6_2atmpS721 = _M0L1jS279 + 2;
        _M0L1iS278 = _M0L6_2atmpS720;
        _M0L1jS279 = _M0L6_2atmpS721;
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

int32_t _M0MPC14uint4UInt8to__byte(uint32_t _M0L4selfS267) {
  int32_t _M0L6_2atmpS712;
  #line 2519 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS712 = *(int32_t*)&_M0L4selfS267;
  return _M0L6_2atmpS712 & 0xff;
}

moonbit_string_t* _M0MPC15array5Array6bufferGsE(
  struct _M0TPB5ArrayGsE* _M0L4selfS265
) {
  moonbit_string_t* _M0L8_2afieldS911;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS911 = _M0L4selfS265->$0;
  moonbit_incref(_M0L8_2afieldS911);
  return _M0L8_2afieldS911;
}

struct _M0TUsiE** _M0MPC15array5Array6bufferGUsiEE(
  struct _M0TPB5ArrayGUsiEE* _M0L4selfS266
) {
  struct _M0TUsiE** _M0L8_2afieldS912;
  #line 184 "/home/developer/.moon/lib/core/builtin/arraycore_nonjs.mbt"
  _M0L8_2afieldS912 = _M0L4selfS266->$0;
  moonbit_incref(_M0L8_2afieldS912);
  return _M0L8_2afieldS912;
}

moonbit_string_t _M0MPC13int3Int18to__string_2einner(
  int32_t _M0L4selfS249,
  int32_t _M0L5radixS248
) {
  int32_t _if__result_980;
  int32_t _M0L12is__negativeS250;
  uint32_t _M0L3numS251;
  uint16_t* _M0L6bufferS252;
  #line 209 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5radixS248 < 2) {
    _if__result_980 = 1;
  } else {
    _if__result_980 = _M0L5radixS248 > 36;
  }
  if (_if__result_980) {
    #line 213 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_9.data);
  }
  if (_M0L4selfS249 == 0) {
    return (moonbit_string_t)moonbit_string_literal_10.data;
  }
  _M0L12is__negativeS250 = _M0L4selfS249 < 0;
  if (_M0L12is__negativeS250) {
    int32_t _M0L6_2atmpS711 = -_M0L4selfS249;
    _M0L3numS251 = *(uint32_t*)&_M0L6_2atmpS711;
  } else {
    _M0L3numS251 = *(uint32_t*)&_M0L4selfS249;
  }
  switch (_M0L5radixS248) {
    case 10: {
      int32_t _M0L10digit__lenS253;
      int32_t _M0L6_2atmpS708;
      int32_t _M0L10total__lenS254;
      uint16_t* _M0L6bufferS255;
      int32_t _M0L12digit__startS256;
      #line 235 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS253 = _M0FPB12dec__count32(_M0L3numS251);
      if (_M0L12is__negativeS250) {
        _M0L6_2atmpS708 = 1;
      } else {
        _M0L6_2atmpS708 = 0;
      }
      _M0L10total__lenS254 = _M0L10digit__lenS253 + _M0L6_2atmpS708;
      _M0L6bufferS255
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS254, 0);
      if (_M0L12is__negativeS250) {
        _M0L12digit__startS256 = 1;
      } else {
        _M0L12digit__startS256 = 0;
      }
      #line 239 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__dec(_M0L6bufferS255, _M0L3numS251, _M0L12digit__startS256, _M0L10total__lenS254);
      _M0L6bufferS252 = _M0L6bufferS255;
      break;
    }
    
    case 16: {
      int32_t _M0L10digit__lenS257;
      int32_t _M0L6_2atmpS709;
      int32_t _M0L10total__lenS258;
      uint16_t* _M0L6bufferS259;
      int32_t _M0L12digit__startS260;
      #line 243 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS257 = _M0FPB12hex__count32(_M0L3numS251);
      if (_M0L12is__negativeS250) {
        _M0L6_2atmpS709 = 1;
      } else {
        _M0L6_2atmpS709 = 0;
      }
      _M0L10total__lenS258 = _M0L10digit__lenS257 + _M0L6_2atmpS709;
      _M0L6bufferS259
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS258, 0);
      if (_M0L12is__negativeS250) {
        _M0L12digit__startS260 = 1;
      } else {
        _M0L12digit__startS260 = 0;
      }
      #line 247 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB20int__to__string__hex(_M0L6bufferS259, _M0L3numS251, _M0L12digit__startS260, _M0L10total__lenS258);
      _M0L6bufferS252 = _M0L6bufferS259;
      break;
    }
    default: {
      int32_t _M0L10digit__lenS261;
      int32_t _M0L6_2atmpS710;
      int32_t _M0L10total__lenS262;
      uint16_t* _M0L6bufferS263;
      int32_t _M0L12digit__startS264;
      #line 251 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0L10digit__lenS261
      = _M0FPB14radix__count32(_M0L3numS251, _M0L5radixS248);
      if (_M0L12is__negativeS250) {
        _M0L6_2atmpS710 = 1;
      } else {
        _M0L6_2atmpS710 = 0;
      }
      _M0L10total__lenS262 = _M0L10digit__lenS261 + _M0L6_2atmpS710;
      _M0L6bufferS263
      = (uint16_t*)moonbit_make_string(_M0L10total__lenS262, 0);
      if (_M0L12is__negativeS250) {
        _M0L12digit__startS264 = 1;
      } else {
        _M0L12digit__startS264 = 0;
      }
      #line 255 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
      _M0FPB24int__to__string__generic(_M0L6bufferS263, _M0L3numS251, _M0L12digit__startS264, _M0L10total__lenS262, _M0L5radixS248);
      _M0L6bufferS252 = _M0L6bufferS263;
      break;
    }
  }
  if (_M0L12is__negativeS250) {
    _M0L6bufferS252[0] = 45;
  }
  return _M0L6bufferS252;
}

int32_t _M0FPB14radix__count32(
  uint32_t _M0L5valueS242,
  int32_t _M0L5radixS244
) {
  uint32_t _M0L4baseS243;
  uint32_t _M0L3numS245;
  int32_t _M0L5countS246;
  #line 189 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS242 == 0u) {
    return 1;
  }
  _M0L4baseS243 = *(uint32_t*)&_M0L5radixS244;
  _M0L3numS245 = _M0L5valueS242;
  _M0L5countS246 = 0;
  while (1) {
    if (_M0L3numS245 > 0u) {
      uint32_t _M0L6_2atmpS706 = _M0L3numS245 / _M0L4baseS243;
      int32_t _M0L6_2atmpS707 = _M0L5countS246 + 1;
      _M0L3numS245 = _M0L6_2atmpS706;
      _M0L5countS246 = _M0L6_2atmpS707;
      continue;
    } else {
      return _M0L5countS246;
    }
    break;
  }
}

int32_t _M0FPB12hex__count32(uint32_t _M0L5valueS240) {
  #line 177 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS240 == 0u) {
    return 1;
  } else {
    int32_t _M0L14leading__zerosS241;
    int32_t _M0L6_2atmpS705;
    int32_t _M0L6_2atmpS704;
    #line 182 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L14leading__zerosS241 = moonbit_clz32(_M0L5valueS240);
    _M0L6_2atmpS705 = 31 - _M0L14leading__zerosS241;
    _M0L6_2atmpS704 = _M0L6_2atmpS705 / 4;
    return _M0L6_2atmpS704 + 1;
  }
}

int32_t _M0FPB12dec__count32(uint32_t _M0L5valueS239) {
  #line 143 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  if (_M0L5valueS239 >= 100000u) {
    if (_M0L5valueS239 >= 10000000u) {
      if (_M0L5valueS239 >= 1000000000u) {
        return 10;
      } else if (_M0L5valueS239 >= 100000000u) {
        return 9;
      } else {
        return 8;
      }
    } else if (_M0L5valueS239 >= 1000000u) {
      return 7;
    } else {
      return 6;
    }
  } else if (_M0L5valueS239 >= 1000u) {
    if (_M0L5valueS239 >= 10000u) {
      return 5;
    } else {
      return 4;
    }
  } else if (_M0L5valueS239 >= 100u) {
    return 3;
  } else if (_M0L5valueS239 >= 10u) {
    return 2;
  } else {
    return 1;
  }
}

int32_t _M0FPB20int__to__string__dec(
  uint16_t* _M0L6bufferS225,
  uint32_t _M0L3numS237,
  int32_t _M0L12digit__startS226,
  int32_t _M0L10total__lenS238
) {
  int32_t _M0L6_2atmpS703;
  uint32_t _M0L3numS215;
  int32_t _M0L6offsetS216;
  #line 88 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS703 = _M0L10total__lenS238 - _M0L12digit__startS226;
  _M0L3numS215 = _M0L3numS237;
  _M0L6offsetS216 = _M0L6_2atmpS703;
  while (1) {
    if (_M0L3numS215 >= 10000u) {
      uint32_t _M0L1tS217 = _M0L3numS215 / 10000u;
      uint32_t _M0L6_2atmpS680 = _M0L3numS215 % 10000u;
      int32_t _M0L1rS218 = *(int32_t*)&_M0L6_2atmpS680;
      int32_t _M0L2d1S219 = _M0L1rS218 / 100;
      int32_t _M0L2d2S220 = _M0L1rS218 % 100;
      int32_t _M0L6_2atmpS679 = _M0L2d1S219 / 10;
      int32_t _M0L6_2atmpS678 = 48 + _M0L6_2atmpS679;
      int32_t _M0L6d1__hiS221 = (uint16_t)_M0L6_2atmpS678;
      int32_t _M0L6_2atmpS677 = _M0L2d1S219 % 10;
      int32_t _M0L6_2atmpS676 = 48 + _M0L6_2atmpS677;
      int32_t _M0L6d1__loS222 = (uint16_t)_M0L6_2atmpS676;
      int32_t _M0L6_2atmpS675 = _M0L2d2S220 / 10;
      int32_t _M0L6_2atmpS674 = 48 + _M0L6_2atmpS675;
      int32_t _M0L6d2__hiS223 = (uint16_t)_M0L6_2atmpS674;
      int32_t _M0L6_2atmpS673 = _M0L2d2S220 % 10;
      int32_t _M0L6_2atmpS672 = 48 + _M0L6_2atmpS673;
      int32_t _M0L6d2__loS224 = (uint16_t)_M0L6_2atmpS672;
      int32_t _M0L6_2atmpS664 = _M0L12digit__startS226 + _M0L6offsetS216;
      int32_t _M0L6_2atmpS663 = _M0L6_2atmpS664 - 4;
      int32_t _M0L6_2atmpS666;
      int32_t _M0L6_2atmpS665;
      int32_t _M0L6_2atmpS668;
      int32_t _M0L6_2atmpS667;
      int32_t _M0L6_2atmpS670;
      int32_t _M0L6_2atmpS669;
      int32_t _M0L6_2atmpS671;
      _M0L6bufferS225[_M0L6_2atmpS663] = _M0L6d1__hiS221;
      _M0L6_2atmpS666 = _M0L12digit__startS226 + _M0L6offsetS216;
      _M0L6_2atmpS665 = _M0L6_2atmpS666 - 3;
      _M0L6bufferS225[_M0L6_2atmpS665] = _M0L6d1__loS222;
      _M0L6_2atmpS668 = _M0L12digit__startS226 + _M0L6offsetS216;
      _M0L6_2atmpS667 = _M0L6_2atmpS668 - 2;
      _M0L6bufferS225[_M0L6_2atmpS667] = _M0L6d2__hiS223;
      _M0L6_2atmpS670 = _M0L12digit__startS226 + _M0L6offsetS216;
      _M0L6_2atmpS669 = _M0L6_2atmpS670 - 1;
      _M0L6bufferS225[_M0L6_2atmpS669] = _M0L6d2__loS224;
      _M0L6_2atmpS671 = _M0L6offsetS216 - 4;
      _M0L3numS215 = _M0L1tS217;
      _M0L6offsetS216 = _M0L6_2atmpS671;
      continue;
    } else {
      int32_t _M0L6_2atmpS702 = *(int32_t*)&_M0L3numS215;
      int32_t _M0L9remainingS228 = _M0L6_2atmpS702;
      int32_t _M0L6offsetS229 = _M0L6offsetS216;
      while (1) {
        if (_M0L9remainingS228 >= 100) {
          int32_t _M0L1tS230 = _M0L9remainingS228 / 100;
          int32_t _M0L1dS231 = _M0L9remainingS228 % 100;
          int32_t _M0L6_2atmpS689 = _M0L1dS231 / 10;
          int32_t _M0L6_2atmpS688 = 48 + _M0L6_2atmpS689;
          int32_t _M0L5d__hiS232 = (uint16_t)_M0L6_2atmpS688;
          int32_t _M0L6_2atmpS687 = _M0L1dS231 % 10;
          int32_t _M0L6_2atmpS686 = 48 + _M0L6_2atmpS687;
          int32_t _M0L5d__loS233 = (uint16_t)_M0L6_2atmpS686;
          int32_t _M0L6_2atmpS682 = _M0L12digit__startS226 + _M0L6offsetS229;
          int32_t _M0L6_2atmpS681 = _M0L6_2atmpS682 - 2;
          int32_t _M0L6_2atmpS684;
          int32_t _M0L6_2atmpS683;
          int32_t _M0L6_2atmpS685;
          _M0L6bufferS225[_M0L6_2atmpS681] = _M0L5d__hiS232;
          _M0L6_2atmpS684 = _M0L12digit__startS226 + _M0L6offsetS229;
          _M0L6_2atmpS683 = _M0L6_2atmpS684 - 1;
          _M0L6bufferS225[_M0L6_2atmpS683] = _M0L5d__loS233;
          _M0L6_2atmpS685 = _M0L6offsetS229 - 2;
          _M0L9remainingS228 = _M0L1tS230;
          _M0L6offsetS229 = _M0L6_2atmpS685;
          continue;
        } else if (_M0L9remainingS228 >= 10) {
          int32_t _M0L6_2atmpS697 = _M0L9remainingS228 / 10;
          int32_t _M0L6_2atmpS696 = 48 + _M0L6_2atmpS697;
          int32_t _M0L5d__hiS235 = (uint16_t)_M0L6_2atmpS696;
          int32_t _M0L6_2atmpS695 = _M0L9remainingS228 % 10;
          int32_t _M0L6_2atmpS694 = 48 + _M0L6_2atmpS695;
          int32_t _M0L5d__loS236 = (uint16_t)_M0L6_2atmpS694;
          int32_t _M0L6_2atmpS691 = _M0L12digit__startS226 + _M0L6offsetS229;
          int32_t _M0L6_2atmpS690 = _M0L6_2atmpS691 - 2;
          int32_t _M0L6_2atmpS693;
          int32_t _M0L6_2atmpS692;
          _M0L6bufferS225[_M0L6_2atmpS690] = _M0L5d__hiS235;
          _M0L6_2atmpS693 = _M0L12digit__startS226 + _M0L6offsetS229;
          _M0L6_2atmpS692 = _M0L6_2atmpS693 - 1;
          _M0L6bufferS225[_M0L6_2atmpS692] = _M0L5d__loS236;
        } else {
          int32_t _M0L6_2atmpS701 = _M0L12digit__startS226 + _M0L6offsetS229;
          int32_t _M0L6_2atmpS698 = _M0L6_2atmpS701 - 1;
          int32_t _M0L6_2atmpS700 = 48 + _M0L9remainingS228;
          int32_t _M0L6_2atmpS699 = (uint16_t)_M0L6_2atmpS700;
          _M0L6bufferS225[_M0L6_2atmpS698] = _M0L6_2atmpS699;
        }
        break;
      }
    }
    break;
  }
  return 0;
}

int32_t _M0FPB24int__to__string__generic(
  uint16_t* _M0L6bufferS205,
  uint32_t _M0L3numS209,
  int32_t _M0L12digit__startS206,
  int32_t _M0L10total__lenS208,
  int32_t _M0L5radixS199
) {
  uint32_t _M0L4baseS198;
  int32_t _M0L6_2atmpS648;
  int32_t _M0L6_2atmpS647;
  #line 57 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L4baseS198 = *(uint32_t*)&_M0L5radixS199;
  _M0L6_2atmpS648 = _M0L5radixS199 - 1;
  _M0L6_2atmpS647 = _M0L5radixS199 & _M0L6_2atmpS648;
  if (_M0L6_2atmpS647 == 0) {
    int32_t _M0L5shiftS200;
    uint32_t _M0L4maskS201;
    int32_t _M0L6_2atmpS655;
    int32_t _M0L6offsetS202;
    uint32_t _M0L1nS203;
    #line 68 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
    _M0L5shiftS200 = moonbit_ctz32(_M0L5radixS199);
    _M0L4maskS201 = _M0L4baseS198 - 1u;
    _M0L6_2atmpS655 = _M0L10total__lenS208 - _M0L12digit__startS206;
    _M0L6offsetS202 = _M0L6_2atmpS655;
    _M0L1nS203 = _M0L3numS209;
    while (1) {
      if (_M0L1nS203 > 0u) {
        uint32_t _M0L6_2atmpS654 = _M0L1nS203 & _M0L4maskS201;
        int32_t _M0L5digitS204 = *(int32_t*)&_M0L6_2atmpS654;
        int32_t _M0L6_2atmpS651 = _M0L12digit__startS206 + _M0L6offsetS202;
        int32_t _M0L6_2atmpS649 = _M0L6_2atmpS651 - 1;
        int32_t _M0L6_2atmpS650 =
          ((moonbit_string_t)moonbit_string_literal_11.data)[_M0L5digitS204];
        int32_t _M0L6_2atmpS652;
        uint32_t _M0L6_2atmpS653;
        _M0L6bufferS205[_M0L6_2atmpS649] = _M0L6_2atmpS650;
        _M0L6_2atmpS652 = _M0L6offsetS202 - 1;
        _M0L6_2atmpS653 = _M0L1nS203 >> (_M0L5shiftS200 & 31);
        _M0L6offsetS202 = _M0L6_2atmpS652;
        _M0L1nS203 = _M0L6_2atmpS653;
        continue;
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS662 = _M0L10total__lenS208 - _M0L12digit__startS206;
    int32_t _M0L6offsetS210 = _M0L6_2atmpS662;
    uint32_t _M0L1nS211 = _M0L3numS209;
    while (1) {
      if (_M0L1nS211 > 0u) {
        uint32_t _M0L1qS212 = _M0L1nS211 / _M0L4baseS198;
        uint32_t _M0L6_2atmpS661 = _M0L1qS212 * _M0L4baseS198;
        uint32_t _M0L6_2atmpS660 = _M0L1nS211 - _M0L6_2atmpS661;
        int32_t _M0L5digitS213 = *(int32_t*)&_M0L6_2atmpS660;
        int32_t _M0L6_2atmpS658 = _M0L12digit__startS206 + _M0L6offsetS210;
        int32_t _M0L6_2atmpS656 = _M0L6_2atmpS658 - 1;
        int32_t _M0L6_2atmpS657 =
          ((moonbit_string_t)moonbit_string_literal_11.data)[_M0L5digitS213];
        int32_t _M0L6_2atmpS659;
        _M0L6bufferS205[_M0L6_2atmpS656] = _M0L6_2atmpS657;
        _M0L6_2atmpS659 = _M0L6offsetS210 - 1;
        _M0L6offsetS210 = _M0L6_2atmpS659;
        _M0L1nS211 = _M0L1qS212;
        continue;
      }
      break;
    }
  }
  return 0;
}

int32_t _M0FPB20int__to__string__hex(
  uint16_t* _M0L6bufferS192,
  uint32_t _M0L3numS197,
  int32_t _M0L12digit__startS193,
  int32_t _M0L10total__lenS196
) {
  int32_t _M0L6_2atmpS646;
  int32_t _M0L6offsetS187;
  uint32_t _M0L1nS188;
  #line 29 "/home/developer/.moon/lib/core/builtin/to_string.mbt"
  _M0L6_2atmpS646 = _M0L10total__lenS196 - _M0L12digit__startS193;
  _M0L6offsetS187 = _M0L6_2atmpS646;
  _M0L1nS188 = _M0L3numS197;
  while (1) {
    if (_M0L6offsetS187 >= 2) {
      uint32_t _M0L6_2atmpS643 = _M0L1nS188 & 255u;
      int32_t _M0L9byte__valS189 = *(int32_t*)&_M0L6_2atmpS643;
      int32_t _M0L2hiS190 = _M0L9byte__valS189 / 16;
      int32_t _M0L2loS191 = _M0L9byte__valS189 % 16;
      int32_t _M0L6_2atmpS637 = _M0L12digit__startS193 + _M0L6offsetS187;
      int32_t _M0L6_2atmpS635 = _M0L6_2atmpS637 - 2;
      int32_t _M0L6_2atmpS636 =
        ((moonbit_string_t)moonbit_string_literal_11.data)[_M0L2hiS190];
      int32_t _M0L6_2atmpS640;
      int32_t _M0L6_2atmpS638;
      int32_t _M0L6_2atmpS639;
      int32_t _M0L6_2atmpS641;
      uint32_t _M0L6_2atmpS642;
      _M0L6bufferS192[_M0L6_2atmpS635] = _M0L6_2atmpS636;
      _M0L6_2atmpS640 = _M0L12digit__startS193 + _M0L6offsetS187;
      _M0L6_2atmpS638 = _M0L6_2atmpS640 - 1;
      _M0L6_2atmpS639
      = ((moonbit_string_t)moonbit_string_literal_11.data)[
        _M0L2loS191
      ];
      _M0L6bufferS192[_M0L6_2atmpS638] = _M0L6_2atmpS639;
      _M0L6_2atmpS641 = _M0L6offsetS187 - 2;
      _M0L6_2atmpS642 = _M0L1nS188 >> 8;
      _M0L6offsetS187 = _M0L6_2atmpS641;
      _M0L1nS188 = _M0L6_2atmpS642;
      continue;
    } else if (_M0L6offsetS187 == 1) {
      uint32_t _M0L6_2atmpS645 = _M0L1nS188 & 15u;
      int32_t _M0L6nibbleS195 = *(int32_t*)&_M0L6_2atmpS645;
      int32_t _M0L6_2atmpS644 =
        ((moonbit_string_t)moonbit_string_literal_11.data)[_M0L6nibbleS195];
      _M0L6bufferS192[_M0L12digit__startS193] = _M0L6_2atmpS644;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(
  void* _M0L4selfS186
) {
  struct _M0TPB13StringBuilder* _M0L6loggerS185;
  struct _M0TPB6Logger _M0L6_2atmpS634;
  moonbit_string_t _result_987;
  #line 165 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 166 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS185 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  moonbit_incref(_M0L6loggerS185);
  _M0L6_2atmpS634
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L6loggerS185
  };
  #line 167 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB7FailurePB4Show6output(_M0L4selfS186, _M0L6_2atmpS634);
  if (_M0L6_2atmpS634.$1) {
    moonbit_decref(_M0L6_2atmpS634.$1);
  }
  #line 168 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _result_987 = _M0MPB13StringBuilder10to__string(_M0L6loggerS185);
  moonbit_decref(_M0L6loggerS185);
  return _result_987;
}

int32_t _M0IP016_24default__implPB4Show6outputGsE(
  moonbit_string_t _M0L4selfS182,
  struct _M0TPB6Logger _M0L6loggerS181
) {
  moonbit_string_t _M0L6_2atmpS632;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS632 = _M0IPC16string6StringPB4Show10to__string(_M0L4selfS182);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS181.$0->$method_0(_M0L6loggerS181.$1, _M0L6_2atmpS632);
  moonbit_decref(_M0L6_2atmpS632);
  return 0;
}

int32_t _M0IP016_24default__implPB4Show6outputGiE(
  int32_t _M0L4selfS184,
  struct _M0TPB6Logger _M0L6loggerS183
) {
  moonbit_string_t _M0L6_2atmpS633;
  #line 159 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS633 = _M0IPC13int3IntPB4Show10to__string(_M0L4selfS184);
  #line 160 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6loggerS183.$0->$method_0(_M0L6loggerS183.$1, _M0L6_2atmpS633);
  moonbit_decref(_M0L6_2atmpS633);
  return 0;
}

int32_t _M0MPC16string10StringView13start__offset(
  struct _M0TPC16string10StringView _M0L4selfS180
) {
  #line 99 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  return _M0L4selfS180.$1;
}

moonbit_string_t _M0MPC16string10StringView4data(
  struct _M0TPC16string10StringView _M0L4selfS179
) {
  moonbit_string_t _M0L8_2afieldS913;
  #line 92 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L8_2afieldS913 = _M0L4selfS179.$0;
  moonbit_incref(_M0L8_2afieldS913);
  return _M0L8_2afieldS913;
}

int32_t _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS175,
  moonbit_string_t _M0L5valueS176,
  int32_t _M0L5startS177,
  int32_t _M0L3lenS178
) {
  int32_t _M0L6_2atmpS631;
  int64_t _M0L6_2atmpS630;
  struct _M0TPC16string10StringView _M0L6_2atmpS629;
  #line 122 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS631 = _M0L5startS177 + _M0L3lenS178;
  _M0L6_2atmpS630 = (int64_t)_M0L6_2atmpS631;
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L6_2atmpS629
  = _M0MPC16string6String11sub_2einner(_M0L5valueS176, _M0L5startS177, _M0L6_2atmpS630);
  #line 123 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L4selfS175, _M0L6_2atmpS629);
  moonbit_decref(_M0L6_2atmpS629.$0);
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string6String11sub_2einner(
  moonbit_string_t _M0L4selfS168,
  int32_t _M0L5startS174,
  int64_t _M0L3endS170
) {
  int32_t _M0L3lenS167;
  int32_t _M0L3endS169;
  int32_t _M0L5startS173;
  int32_t _if__result_988;
  #line 755 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3lenS167 = Moonbit_array_length(_M0L4selfS168);
  if (_M0L3endS170 == 4294967296ll) {
    _M0L3endS169 = _M0L3lenS167;
  } else {
    int64_t _M0L7_2aSomeS171 = _M0L3endS170;
    int32_t _M0L6_2aendS172 = (int32_t)_M0L7_2aSomeS171;
    if (_M0L6_2aendS172 < 0) {
      _M0L3endS169 = _M0L3lenS167 + _M0L6_2aendS172;
    } else {
      _M0L3endS169 = _M0L6_2aendS172;
    }
  }
  if (_M0L5startS174 < 0) {
    _M0L5startS173 = _M0L3lenS167 + _M0L5startS174;
  } else {
    _M0L5startS173 = _M0L5startS174;
  }
  if (_M0L5startS173 >= 0) {
    if (_M0L5startS173 <= _M0L3endS169) {
      _if__result_988 = _M0L3endS169 <= _M0L3lenS167;
    } else {
      _if__result_988 = 0;
    }
  } else {
    _if__result_988 = 0;
  }
  if (_if__result_988) {
    if (_M0L5startS173 < _M0L3lenS167) {
      int32_t _M0L6_2atmpS626 = _M0L4selfS168[_M0L5startS173];
      int32_t _M0L6_2atmpS625;
      #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS625
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS626);
      if (!_M0L6_2atmpS625) {
        
      } else {
        #line 765 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L3endS169 < _M0L3lenS167) {
      int32_t _M0L6_2atmpS628 = _M0L4selfS168[_M0L3endS169];
      int32_t _M0L6_2atmpS627;
      #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS627
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS628);
      if (!_M0L6_2atmpS627) {
        
      } else {
        #line 768 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    moonbit_incref(_M0L4selfS168);
    return (struct _M0TPC16string10StringView){.$0 = _M0L4selfS168,
                                                 .$1 = _M0L5startS173,
                                                 .$2 = _M0L3endS169};
  } else {
    #line 763 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

int32_t _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS166,
  struct _M0TPB4Show _M0L4showS165
) {
  struct _M0TPB6Logger _M0L6_2atmpS624;
  #line 116 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS166);
  _M0L6_2atmpS624
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS166
  };
  #line 117 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS165.$0->$method_0(_M0L4showS165.$1, _M0L6_2atmpS624);
  if (_M0L6_2atmpS624.$1) {
    moonbit_decref(_M0L6_2atmpS624.$1);
  }
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(
  struct _M0TPB13StringBuilder* _M0L4selfS164,
  struct _M0TPB4Show _M0L4showS163
) {
  struct _M0TPB6Logger _M0L6_2atmpS623;
  #line 111 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  moonbit_incref(_M0L4selfS164);
  _M0L6_2atmpS623
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS164
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0L4showS163.$0->$method_0(_M0L4showS163.$1, _M0L6_2atmpS623);
  if (_M0L6_2atmpS623.$1) {
    moonbit_decref(_M0L6_2atmpS623.$1);
  }
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger13write__string(
  struct _M0TPB13StringBuilder* _M0L4selfS162,
  moonbit_string_t _M0L3strS161
) {
  int32_t _M0L8str__lenS160;
  int32_t _M0L3lenS618;
  int32_t _M0L6_2atmpS617;
  uint16_t* _M0L4dataS619;
  int32_t _M0L3lenS620;
  int32_t _M0L3lenS622;
  int32_t _M0L6_2atmpS621;
  #line 86 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L8str__lenS160 = Moonbit_array_length(_M0L3strS161);
  _M0L3lenS618 = _M0L4selfS162->$1;
  _M0L6_2atmpS617 = _M0L3lenS618 + _M0L8str__lenS160;
  #line 88 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS162, _M0L6_2atmpS617);
  _M0L4dataS619 = _M0L4selfS162->$0;
  _M0L3lenS620 = _M0L4selfS162->$1;
  moonbit_incref(_M0L4dataS619);
  #line 89 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0MPC15array10FixedArray26unsafe__blit__from__string(_M0L4dataS619, _M0L3lenS620, _M0L3strS161, 0, _M0L8str__lenS160);
  moonbit_decref(_M0L4dataS619);
  _M0L3lenS622 = _M0L4selfS162->$1;
  _M0L6_2atmpS621 = _M0L3lenS622 + _M0L8str__lenS160;
  _M0L4selfS162->$1 = _M0L6_2atmpS621;
  return 0;
}

int32_t _M0MPC15array10FixedArray26unsafe__blit__from__string(
  uint16_t* _M0L4selfS156,
  int32_t _M0L11dst__offsetS159,
  moonbit_string_t _M0L3strS157,
  int32_t _M0L11str__offsetS152,
  int32_t _M0L3lenS153
) {
  int32_t _M0L16end__str__offsetS151;
  int32_t _M0L1iS154;
  int32_t _M0L1jS155;
  #line 71 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L16end__str__offsetS151 = _M0L11str__offsetS152 + _M0L3lenS153;
  _M0L1iS154 = _M0L11str__offsetS152;
  _M0L1jS155 = _M0L11dst__offsetS159;
  while (1) {
    if (_M0L1iS154 < _M0L16end__str__offsetS151) {
      int32_t _M0L6_2atmpS614 = _M0L3strS157[_M0L1iS154];
      int32_t _M0L6_2atmpS615;
      int32_t _M0L6_2atmpS616;
      if (
        _M0L1jS155 < 0 || _M0L1jS155 >= Moonbit_array_length(_M0L4selfS156)
      ) {
        #line 80 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
        moonbit_panic();
      }
      _M0L4selfS156[_M0L1jS155] = _M0L6_2atmpS614;
      _M0L6_2atmpS615 = _M0L1iS154 + 1;
      _M0L6_2atmpS616 = _M0L1jS155 + 1;
      _M0L1iS154 = _M0L6_2atmpS615;
      _M0L1jS155 = _M0L6_2atmpS616;
      continue;
    }
    break;
  }
  return 0;
}

moonbit_string_t _M0MPC16string6String14escape_2einner(
  moonbit_string_t _M0L4selfS149,
  int32_t _M0L5quoteS150
) {
  struct _M0TPB13StringBuilder* _M0L3bufS148;
  int32_t _M0L6_2atmpS613;
  struct _M0TPC16string10StringView _M0L6_2atmpS611;
  struct _M0TPB6Logger _M0L6_2atmpS612;
  moonbit_string_t _result_990;
  #line 110 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 111 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L3bufS148 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  _M0L6_2atmpS613 = Moonbit_array_length(_M0L4selfS149);
  moonbit_incref(_M0L4selfS149);
  _M0L6_2atmpS611
  = (struct _M0TPC16string10StringView){
    .$0 = _M0L4selfS149, .$1 = 0, .$2 = _M0L6_2atmpS613
  };
  moonbit_incref(_M0L3bufS148);
  _M0L6_2atmpS612
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L3bufS148
  };
  #line 112 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0MPC16string10StringView18escape__to_2einner(_M0L6_2atmpS611, _M0L6_2atmpS612, _M0L5quoteS150);
  moonbit_decref(_M0L6_2atmpS611.$0);
  if (_M0L6_2atmpS612.$1) {
    moonbit_decref(_M0L6_2atmpS612.$1);
  }
  #line 113 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_990 = _M0MPB13StringBuilder10to__string(_M0L3bufS148);
  moonbit_decref(_M0L3bufS148);
  return _result_990;
}

int32_t _M0MPC16string10StringView18escape__to_2einner(
  struct _M0TPC16string10StringView _M0L4selfS140,
  struct _M0TPB6Logger _M0L6loggerS138,
  int32_t _M0L5quoteS137
) {
  int32_t _M0L3endS609;
  int32_t _M0L5startS610;
  int32_t _M0L3lenS139;
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS141;
  int32_t _M0L1iS142;
  int32_t _M0L3segS143;
  #line 144 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L5quoteS137) {
    #line 150 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS138.$0->$method_3(_M0L6loggerS138.$1, 34);
  }
  _M0L3endS609 = _M0L4selfS140.$2;
  _M0L5startS610 = _M0L4selfS140.$1;
  _M0L3lenS139 = _M0L3endS609 - _M0L5startS610;
  moonbit_incref(_M0L4selfS140.$0);
  if (_M0L6loggerS138.$1) {
    moonbit_incref(_M0L6loggerS138.$1);
  }
  _M0L6_2aenvS141
  = (struct _M0TURPC16string10StringViewRPB6LoggerE*)moonbit_malloc(sizeof(struct _M0TURPC16string10StringViewRPB6LoggerE));
  Moonbit_object_header(_M0L6_2aenvS141)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 25, 0);
  _M0L6_2aenvS141->$0 = _M0L4selfS140;
  _M0L6_2aenvS141->$1 = _M0L6loggerS138;
  _M0L1iS142 = 0;
  _M0L3segS143 = 0;
  _2afor_144:;
  while (1) {
    moonbit_string_t _M0L3strS606;
    int32_t _M0L5startS608;
    int32_t _M0L6_2atmpS607;
    int32_t _M0L4codeS145;
    int32_t _M0L1cS147;
    int32_t _M0L6_2atmpS590;
    int32_t _M0L6_2atmpS591;
    int32_t _M0L6_2atmpS592;
    int32_t _tmp_994;
    int32_t _tmp_995;
    if (_M0L1iS142 >= _M0L3lenS139) {
      #line 160 "/home/developer/.moon/lib/core/builtin/show.mbt"
      _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
      moonbit_decref(_M0L6_2aenvS141);
      break;
    }
    _M0L3strS606 = _M0L4selfS140.$0;
    _M0L5startS608 = _M0L4selfS140.$1;
    _M0L6_2atmpS607 = _M0L5startS608 + _M0L1iS142;
    _M0L4codeS145 = _M0L3strS606[_M0L6_2atmpS607];
    switch (_M0L4codeS145) {
      case 34: {
        _M0L1cS147 = _M0L4codeS145;
        goto join_146;
        break;
      }
      
      case 92: {
        _M0L1cS147 = _M0L4codeS145;
        goto join_146;
        break;
      }
      
      case 10: {
        int32_t _M0L6_2atmpS593;
        int32_t _M0L6_2atmpS594;
        #line 172 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
        #line 173 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS138.$0->$method_0(_M0L6loggerS138.$1, (moonbit_string_t)moonbit_string_literal_12.data);
        _M0L6_2atmpS593 = _M0L1iS142 + 1;
        _M0L6_2atmpS594 = _M0L1iS142 + 1;
        _M0L1iS142 = _M0L6_2atmpS593;
        _M0L3segS143 = _M0L6_2atmpS594;
        goto _2afor_144;
        break;
      }
      
      case 13: {
        int32_t _M0L6_2atmpS595;
        int32_t _M0L6_2atmpS596;
        #line 177 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
        #line 178 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS138.$0->$method_0(_M0L6loggerS138.$1, (moonbit_string_t)moonbit_string_literal_13.data);
        _M0L6_2atmpS595 = _M0L1iS142 + 1;
        _M0L6_2atmpS596 = _M0L1iS142 + 1;
        _M0L1iS142 = _M0L6_2atmpS595;
        _M0L3segS143 = _M0L6_2atmpS596;
        goto _2afor_144;
        break;
      }
      
      case 8: {
        int32_t _M0L6_2atmpS597;
        int32_t _M0L6_2atmpS598;
        #line 182 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
        #line 183 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS138.$0->$method_0(_M0L6loggerS138.$1, (moonbit_string_t)moonbit_string_literal_14.data);
        _M0L6_2atmpS597 = _M0L1iS142 + 1;
        _M0L6_2atmpS598 = _M0L1iS142 + 1;
        _M0L1iS142 = _M0L6_2atmpS597;
        _M0L3segS143 = _M0L6_2atmpS598;
        goto _2afor_144;
        break;
      }
      
      case 9: {
        int32_t _M0L6_2atmpS599;
        int32_t _M0L6_2atmpS600;
        #line 187 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
        #line 188 "/home/developer/.moon/lib/core/builtin/show.mbt"
        _M0L6loggerS138.$0->$method_0(_M0L6loggerS138.$1, (moonbit_string_t)moonbit_string_literal_15.data);
        _M0L6_2atmpS599 = _M0L1iS142 + 1;
        _M0L6_2atmpS600 = _M0L1iS142 + 1;
        _M0L1iS142 = _M0L6_2atmpS599;
        _M0L3segS143 = _M0L6_2atmpS600;
        goto _2afor_144;
        break;
      }
      default: {
        if (_M0L4codeS145 < 32) {
          int32_t _M0L6_2atmpS602;
          moonbit_string_t _M0L6_2atmpS601;
          int32_t _M0L6_2atmpS603;
          int32_t _M0L6_2atmpS604;
          #line 193 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
          #line 194 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS138.$0->$method_0(_M0L6loggerS138.$1, (moonbit_string_t)moonbit_string_literal_16.data);
          _M0L6_2atmpS602 = _M0L4codeS145 & 0xff;
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6_2atmpS601 = _M0MPC14byte4Byte7to__hex(_M0L6_2atmpS602);
          #line 195 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS138.$0->$method_0(_M0L6loggerS138.$1, _M0L6_2atmpS601);
          moonbit_decref(_M0L6_2atmpS601);
          #line 196 "/home/developer/.moon/lib/core/builtin/show.mbt"
          _M0L6loggerS138.$0->$method_3(_M0L6loggerS138.$1, 125);
          _M0L6_2atmpS603 = _M0L1iS142 + 1;
          _M0L6_2atmpS604 = _M0L1iS142 + 1;
          _M0L1iS142 = _M0L6_2atmpS603;
          _M0L3segS143 = _M0L6_2atmpS604;
          goto _2afor_144;
        } else {
          int32_t _M0L6_2atmpS605 = _M0L1iS142 + 1;
          int32_t _tmp_993 = _M0L3segS143;
          _M0L1iS142 = _M0L6_2atmpS605;
          _M0L3segS143 = _tmp_993;
          goto _2afor_144;
        }
        break;
      }
    }
    goto joinlet_992;
    join_146:;
    #line 166 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(_M0L6_2aenvS141, _M0L3segS143, _M0L1iS142);
    #line 167 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS138.$0->$method_3(_M0L6loggerS138.$1, 92);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS590 = _M0MPC16uint166UInt1616unsafe__to__char(_M0L1cS147);
    #line 168 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS138.$0->$method_3(_M0L6loggerS138.$1, _M0L6_2atmpS590);
    _M0L6_2atmpS591 = _M0L1iS142 + 1;
    _M0L6_2atmpS592 = _M0L1iS142 + 1;
    _M0L1iS142 = _M0L6_2atmpS591;
    _M0L3segS143 = _M0L6_2atmpS592;
    continue;
    joinlet_992:;
    _tmp_994 = _M0L1iS142;
    _tmp_995 = _M0L3segS143;
    _M0L1iS142 = _tmp_994;
    _M0L3segS143 = _tmp_995;
    continue;
    break;
  }
  if (_M0L5quoteS137) {
    #line 204 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS138.$0->$method_3(_M0L6loggerS138.$1, 34);
  }
  return 0;
}

int32_t _M0MPC16string10StringView18escape__to_2einnerN14flush__segmentS4141(
  struct _M0TURPC16string10StringViewRPB6LoggerE* _M0L6_2aenvS133,
  int32_t _M0L3segS136,
  int32_t _M0L1iS135
) {
  struct _M0TPB6Logger _M0L6loggerS132;
  struct _M0TPC16string10StringView _M0L4selfS134;
  #line 153 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6loggerS132 = _M0L6_2aenvS133->$1;
  _M0L4selfS134 = _M0L6_2aenvS133->$0;
  if (_M0L1iS135 > _M0L3segS136) {
    int64_t _M0L6_2atmpS589 = (int64_t)_M0L1iS135;
    struct _M0TPC16string10StringView _M0L6_2atmpS588;
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS588
    = _M0MPC16string10StringView11sub_2einner(_M0L4selfS134, _M0L3segS136, _M0L6_2atmpS589);
    #line 155 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6loggerS132.$0->$method_2(_M0L6loggerS132.$1, _M0L6_2atmpS588);
    moonbit_decref(_M0L6_2atmpS588.$0);
  }
  return 0;
}

struct _M0TPC16string10StringView _M0MPC16string10StringView11sub_2einner(
  struct _M0TPC16string10StringView _M0L4selfS125,
  int32_t _M0L5startS131,
  int64_t _M0L3endS127
) {
  moonbit_string_t _M0L3strS587;
  int32_t _M0L8str__lenS124;
  int32_t _M0L8abs__endS126;
  int32_t _M0L10abs__startS130;
  int32_t _M0L5startS575;
  int32_t _if__result_996;
  #line 814 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
  _M0L3strS587 = _M0L4selfS125.$0;
  _M0L8str__lenS124 = Moonbit_array_length(_M0L3strS587);
  if (_M0L3endS127 == 4294967296ll) {
    _M0L8abs__endS126 = _M0L4selfS125.$2;
  } else {
    int64_t _M0L7_2aSomeS128 = _M0L3endS127;
    int32_t _M0L6_2aendS129 = (int32_t)_M0L7_2aSomeS128;
    if (_M0L6_2aendS129 < 0) {
      int32_t _M0L3endS585 = _M0L4selfS125.$2;
      _M0L8abs__endS126 = _M0L3endS585 + _M0L6_2aendS129;
    } else {
      int32_t _M0L5startS586 = _M0L4selfS125.$1;
      _M0L8abs__endS126 = _M0L5startS586 + _M0L6_2aendS129;
    }
  }
  if (_M0L5startS131 < 0) {
    int32_t _M0L3endS583 = _M0L4selfS125.$2;
    _M0L10abs__startS130 = _M0L3endS583 + _M0L5startS131;
  } else {
    int32_t _M0L5startS584 = _M0L4selfS125.$1;
    _M0L10abs__startS130 = _M0L5startS584 + _M0L5startS131;
  }
  _M0L5startS575 = _M0L4selfS125.$1;
  if (_M0L10abs__startS130 >= _M0L5startS575) {
    if (_M0L10abs__startS130 <= _M0L8abs__endS126) {
      int32_t _M0L3endS574 = _M0L4selfS125.$2;
      _if__result_996 = _M0L8abs__endS126 <= _M0L3endS574;
    } else {
      _if__result_996 = 0;
    }
  } else {
    _if__result_996 = 0;
  }
  if (_if__result_996) {
    moonbit_string_t _M0L3strS582;
    if (_M0L10abs__startS130 < _M0L8str__lenS124) {
      moonbit_string_t _M0L3strS578 = _M0L4selfS125.$0;
      int32_t _M0L6_2atmpS577 = _M0L3strS578[_M0L10abs__startS130];
      int32_t _M0L6_2atmpS576;
      #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS576
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS577);
      if (!_M0L6_2atmpS576) {
        
      } else {
        #line 840 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    if (_M0L8abs__endS126 < _M0L8str__lenS124) {
      moonbit_string_t _M0L3strS581 = _M0L4selfS125.$0;
      int32_t _M0L6_2atmpS580 = _M0L3strS581[_M0L8abs__endS126];
      int32_t _M0L6_2atmpS579;
      #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
      _M0L6_2atmpS579
      = _M0MPC16uint166UInt1623is__trailing__surrogate(_M0L6_2atmpS580);
      if (!_M0L6_2atmpS579) {
        
      } else {
        #line 843 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
        moonbit_panic();
      }
    }
    _M0L3strS582 = _M0L4selfS125.$0;
    moonbit_incref(_M0L3strS582);
    return (struct _M0TPC16string10StringView){.$0 = _M0L3strS582,
                                                 .$1 = _M0L10abs__startS130,
                                                 .$2 = _M0L8abs__endS126};
  } else {
    #line 834 "/home/developer/.moon/lib/core/builtin/stringview.mbt"
    moonbit_panic();
  }
}

moonbit_string_t _M0MPC14byte4Byte7to__hex(int32_t _M0L1bS123) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS122;
  int32_t _M0L6_2atmpS571;
  int32_t _M0L6_2atmpS570;
  int32_t _M0L6_2atmpS573;
  int32_t _M0L6_2atmpS572;
  struct _M0TPB13StringBuilder* _M0L6_2atmpS569;
  moonbit_string_t _result_997;
  #line 74 "/home/developer/.moon/lib/core/builtin/show.mbt"
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L7_2aselfS122 = _M0MPB13StringBuilder21StringBuilder_2einner(0);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS571 = _M0IPC14byte4BytePB3Div3div(_M0L1bS123, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS570
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS571);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS122, _M0L6_2atmpS570);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS573 = _M0IPC14byte4BytePB3Mod3mod(_M0L1bS123, 16);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0L6_2atmpS572
  = _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(_M0L6_2atmpS573);
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS122, _M0L6_2atmpS572);
  _M0L6_2atmpS569 = _M0L7_2aselfS122;
  #line 83 "/home/developer/.moon/lib/core/builtin/show.mbt"
  _result_997 = _M0MPB13StringBuilder10to__string(_M0L6_2atmpS569);
  moonbit_decref(_M0L6_2atmpS569);
  return _result_997;
}

int32_t _M0MPC14byte4Byte7to__hexN14to__hex__digitS4156(int32_t _M0L1iS121) {
  #line 75 "/home/developer/.moon/lib/core/builtin/show.mbt"
  if (_M0L1iS121 < 10) {
    int32_t _M0L6_2atmpS566;
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS566 = _M0IPC14byte4BytePB3Add3add(_M0L1iS121, 48);
    #line 77 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS566);
  } else {
    int32_t _M0L6_2atmpS568;
    int32_t _M0L6_2atmpS567;
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS568 = _M0IPC14byte4BytePB3Add3add(_M0L1iS121, 97);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    _M0L6_2atmpS567 = _M0IPC14byte4BytePB3Sub3sub(_M0L6_2atmpS568, 10);
    #line 79 "/home/developer/.moon/lib/core/builtin/show.mbt"
    return _M0MPC14byte4Byte8to__char(_M0L6_2atmpS567);
  }
}

int32_t _M0IPC14byte4BytePB3Sub3sub(
  int32_t _M0L4selfS119,
  int32_t _M0L4thatS120
) {
  int32_t _M0L6_2atmpS564;
  int32_t _M0L6_2atmpS565;
  int32_t _M0L6_2atmpS563;
  #line 120 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS564 = (int32_t)_M0L4selfS119;
  _M0L6_2atmpS565 = (int32_t)_M0L4thatS120;
  _M0L6_2atmpS563 = _M0L6_2atmpS564 - _M0L6_2atmpS565;
  return _M0L6_2atmpS563 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Mod3mod(
  int32_t _M0L4selfS117,
  int32_t _M0L4thatS118
) {
  int32_t _M0L6_2atmpS561;
  int32_t _M0L6_2atmpS562;
  int32_t _M0L6_2atmpS560;
  #line 67 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS561 = (int32_t)_M0L4selfS117;
  _M0L6_2atmpS562 = (int32_t)_M0L4thatS118;
  _M0L6_2atmpS560 = _M0L6_2atmpS561 % _M0L6_2atmpS562;
  return _M0L6_2atmpS560 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Div3div(
  int32_t _M0L4selfS115,
  int32_t _M0L4thatS116
) {
  int32_t _M0L6_2atmpS558;
  int32_t _M0L6_2atmpS559;
  int32_t _M0L6_2atmpS557;
  #line 62 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS558 = (int32_t)_M0L4selfS115;
  _M0L6_2atmpS559 = (int32_t)_M0L4thatS116;
  _M0L6_2atmpS557 = _M0L6_2atmpS558 / _M0L6_2atmpS559;
  return _M0L6_2atmpS557 & 0xff;
}

int32_t _M0IPC14byte4BytePB3Add3add(
  int32_t _M0L4selfS113,
  int32_t _M0L4thatS114
) {
  int32_t _M0L6_2atmpS555;
  int32_t _M0L6_2atmpS556;
  int32_t _M0L6_2atmpS554;
  #line 106 "/home/developer/.moon/lib/core/builtin/byte.mbt"
  _M0L6_2atmpS555 = (int32_t)_M0L4selfS113;
  _M0L6_2atmpS556 = (int32_t)_M0L4thatS114;
  _M0L6_2atmpS554 = _M0L6_2atmpS555 + _M0L6_2atmpS556;
  return _M0L6_2atmpS554 & 0xff;
}

int32_t _M0MPC16uint166UInt1616unsafe__to__char(int32_t _M0L4selfS112) {
  int32_t _M0L6_2atmpS553;
  #line 68 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  _M0L6_2atmpS553 = (int32_t)_M0L4selfS112;
  return _M0L6_2atmpS553;
}

int32_t _M0MPC16uint166UInt1623is__trailing__surrogate(int32_t _M0L4selfS111) {
  #line 45 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  if (_M0L4selfS111 >= 56320) {
    return _M0L4selfS111 <= 57343;
  } else {
    return 0;
  }
}

int32_t _M0IPB13StringBuilderPB6Logger11write__char(
  struct _M0TPB13StringBuilder* _M0L4selfS109,
  int32_t _M0L2chS108
) {
  uint32_t _M0L4codeS107;
  #line 95 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  #line 96 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4codeS107 = _M0MPC14char4Char8to__uint(_M0L2chS108);
  if (_M0L4codeS107 <= 65535u) {
    int32_t _M0L3lenS532 = _M0L4selfS109->$1;
    int32_t _M0L6_2atmpS531 = _M0L3lenS532 + 1;
    uint16_t* _M0L4dataS533;
    int32_t _M0L3lenS534;
    int32_t _M0L6_2atmpS535;
    int32_t _M0L3lenS537;
    int32_t _M0L6_2atmpS536;
    #line 98 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS109, _M0L6_2atmpS531);
    _M0L4dataS533 = _M0L4selfS109->$0;
    _M0L3lenS534 = _M0L4selfS109->$1;
    moonbit_incref(_M0L4dataS533);
    #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS535 = _M0MPC14uint4UInt10to__uint16(_M0L4codeS107);
    if (
      _M0L3lenS534 < 0 || _M0L3lenS534 >= Moonbit_array_length(_M0L4dataS533)
    ) {
      #line 99 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS533[_M0L3lenS534] = _M0L6_2atmpS535;
    moonbit_decref(_M0L4dataS533);
    _M0L3lenS537 = _M0L4selfS109->$1;
    _M0L6_2atmpS536 = _M0L3lenS537 + 1;
    _M0L4selfS109->$1 = _M0L6_2atmpS536;
  } else if (_M0L4codeS107 <= 1114111u) {
    int32_t _M0L3lenS539 = _M0L4selfS109->$1;
    int32_t _M0L6_2atmpS538 = _M0L3lenS539 + 2;
    uint32_t _M0L4codeS110;
    uint16_t* _M0L4dataS540;
    int32_t _M0L3lenS541;
    uint32_t _M0L6_2atmpS544;
    uint32_t _M0L6_2atmpS543;
    int32_t _M0L6_2atmpS542;
    uint16_t* _M0L4dataS545;
    int32_t _M0L3lenS550;
    int32_t _M0L6_2atmpS546;
    uint32_t _M0L6_2atmpS549;
    uint32_t _M0L6_2atmpS548;
    int32_t _M0L6_2atmpS547;
    int32_t _M0L3lenS552;
    int32_t _M0L6_2atmpS551;
    #line 102 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0MPB13StringBuilder19grow__if__necessary(_M0L4selfS109, _M0L6_2atmpS538);
    _M0L4codeS110 = _M0L4codeS107 - 65536u;
    _M0L4dataS540 = _M0L4selfS109->$0;
    _M0L3lenS541 = _M0L4selfS109->$1;
    _M0L6_2atmpS544 = _M0L4codeS110 >> 10;
    _M0L6_2atmpS543 = 55296u + _M0L6_2atmpS544;
    moonbit_incref(_M0L4dataS540);
    #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS542 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS543);
    if (
      _M0L3lenS541 < 0 || _M0L3lenS541 >= Moonbit_array_length(_M0L4dataS540)
    ) {
      #line 104 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS540[_M0L3lenS541] = _M0L6_2atmpS542;
    moonbit_decref(_M0L4dataS540);
    _M0L4dataS545 = _M0L4selfS109->$0;
    _M0L3lenS550 = _M0L4selfS109->$1;
    _M0L6_2atmpS546 = _M0L3lenS550 + 1;
    _M0L6_2atmpS549 = _M0L4codeS110 & 1023u;
    _M0L6_2atmpS548 = 56320u + _M0L6_2atmpS549;
    moonbit_incref(_M0L4dataS545);
    #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS547 = _M0MPC14uint4UInt10to__uint16(_M0L6_2atmpS548);
    if (
      _M0L6_2atmpS546 < 0
      || _M0L6_2atmpS546 >= Moonbit_array_length(_M0L4dataS545)
    ) {
      #line 105 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
      moonbit_panic();
    }
    _M0L4dataS545[_M0L6_2atmpS546] = _M0L6_2atmpS547;
    moonbit_decref(_M0L4dataS545);
    _M0L3lenS552 = _M0L4selfS109->$1;
    _M0L6_2atmpS551 = _M0L3lenS552 + 2;
    _M0L4selfS109->$1 = _M0L6_2atmpS551;
  } else {
    #line 108 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0FPC15abort5abortGuE((moonbit_string_t)moonbit_string_literal_17.data);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder19grow__if__necessary(
  struct _M0TPB13StringBuilder* _M0L4selfS101,
  int32_t _M0L8requiredS102
) {
  uint16_t* _M0L4dataS530;
  int32_t _M0L12current__lenS100;
  int32_t _M0L13enough__spaceS103;
  int32_t _M0L13enough__spaceS104;
  uint16_t* _M0L4dataS526;
  int32_t _M0L6_2atmpS527;
  int32_t _M0L3lenS528;
  uint16_t* _M0L9new__dataS106;
  uint16_t* _M0L6_2aoldS923;
  #line 46 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L4dataS530 = _M0L4selfS101->$0;
  _M0L12current__lenS100 = Moonbit_array_length(_M0L4dataS530);
  if (_M0L8requiredS102 <= _M0L12current__lenS100) {
    return 0;
  }
  _M0L13enough__spaceS104 = _M0L12current__lenS100;
  while (1) {
    if (_M0L13enough__spaceS104 < _M0L8requiredS102) {
      int32_t _M0L6_2atmpS529 = _M0L13enough__spaceS104 * 2;
      _M0L13enough__spaceS104 = _M0L6_2atmpS529;
      continue;
    } else {
      _M0L13enough__spaceS103 = _M0L13enough__spaceS104;
    }
    break;
  }
  _M0L4dataS526 = _M0L4selfS101->$0;
  moonbit_incref(_M0L4dataS526);
  #line 64 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L6_2atmpS527 = _M0IPC16uint166UInt16PB7Default7default();
  _M0L3lenS528 = _M0L4selfS101->$1;
  #line 61 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L9new__dataS106
  = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS526, _M0L13enough__spaceS103, _M0L6_2atmpS527, _M0L3lenS528, 0, 0);
  moonbit_decref(_M0L4dataS526);
  _M0L6_2aoldS923 = _M0L4selfS101->$0;
  moonbit_decref(_M0L6_2aoldS923);
  _M0L4selfS101->$0 = _M0L9new__dataS106;
  return 0;
}

int32_t _M0MPC14uint4UInt10to__uint16(uint32_t _M0L4selfS99) {
  int32_t _M0L6_2atmpS525;
  #line 2676 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS525 = *(int32_t*)&_M0L4selfS99;
  return (uint16_t)_M0L6_2atmpS525;
}

uint32_t _M0MPC14char4Char8to__uint(int32_t _M0L4selfS98) {
  int32_t _M0L6_2atmpS524;
  #line 1254 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS524 = _M0L4selfS98;
  return *(uint32_t*)&_M0L6_2atmpS524;
}

moonbit_string_t _M0MPB13StringBuilder10to__string(
  struct _M0TPB13StringBuilder* _M0L4selfS96
) {
  int32_t _M0L3lenS516;
  uint16_t* _M0L4dataS518;
  int32_t _M0L6_2atmpS517;
  #line 148 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  _M0L3lenS516 = _M0L4selfS96->$1;
  _M0L4dataS518 = _M0L4selfS96->$0;
  _M0L6_2atmpS517 = Moonbit_array_length(_M0L4dataS518);
  if (_M0L3lenS516 == _M0L6_2atmpS517) {
    uint16_t* _M0L4dataS519 = _M0L4selfS96->$0;
    moonbit_incref(_M0L4dataS519);
    return _M0L4dataS519;
  } else {
    uint16_t* _M0L4dataS520 = _M0L4selfS96->$0;
    int32_t _M0L3lenS521 = _M0L4selfS96->$1;
    int32_t _M0L6_2atmpS522;
    int32_t _M0L3lenS523;
    uint16_t* _M0L4dataS97;
    moonbit_incref(_M0L4dataS520);
    #line 155 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L6_2atmpS522 = _M0IPC16uint166UInt16PB7Default7default();
    _M0L3lenS523 = _M0L4selfS96->$1;
    #line 152 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
    _M0L4dataS97
    = _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(_M0L4dataS520, _M0L3lenS521, _M0L6_2atmpS522, _M0L3lenS523, 0, 0);
    moonbit_decref(_M0L4dataS520);
    return _M0L4dataS97;
  }
}

int32_t _M0IPC16uint166UInt16PB7Default7default() {
  #line 176 "/home/developer/.moon/lib/core/builtin/uint16_char.mbt"
  return 0;
}

uint16_t* _M0MPC15array10FixedArray23make__and__blit_2einnerGkE(
  uint16_t* _M0L3srcS93,
  int32_t _M0L13allocate__lenS89,
  int32_t _M0L4initS94,
  int32_t _M0L3lenS90,
  int32_t _M0L11src__offsetS91,
  int32_t _M0L11dst__offsetS92
) {
  int32_t _if__result_999;
  #line 97 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L13allocate__lenS89 >= 0) {
    if (_M0L3lenS90 >= 0) {
      if (_M0L11src__offsetS91 >= 0) {
        if (_M0L11dst__offsetS92 >= 0) {
          int32_t _M0L6_2atmpS512 = _M0L11src__offsetS91 + _M0L3lenS90;
          int32_t _M0L6_2atmpS513 = Moonbit_array_length(_M0L3srcS93);
          if (_M0L6_2atmpS512 <= _M0L6_2atmpS513) {
            int32_t _M0L6_2atmpS511 = _M0L11dst__offsetS92 + _M0L3lenS90;
            _if__result_999 = _M0L6_2atmpS511 <= _M0L13allocate__lenS89;
          } else {
            _if__result_999 = 0;
          }
        } else {
          _if__result_999 = 0;
        }
      } else {
        _if__result_999 = 0;
      }
    } else {
      _if__result_999 = 0;
    }
  } else {
    _if__result_999 = 0;
  }
  if (_if__result_999) {
    moonbit_incref(_M0L3srcS93);
    #line 115 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    return _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(_M0L3srcS93, _M0L13allocate__lenS89, _M0L4initS94, _M0L11src__offsetS91, _M0L11dst__offsetS92, _M0L3lenS90);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS95;
    int32_t _M0L6_2atmpS515;
    moonbit_string_t _M0L6_2atmpS514;
    uint16_t* _result_1000;
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L18_2astring__builderS95
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS95, (moonbit_string_t)moonbit_string_literal_18.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS95, _M0L13allocate__lenS89);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS95, (moonbit_string_t)moonbit_string_literal_19.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS95, _M0L11src__offsetS91);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS95, (moonbit_string_t)moonbit_string_literal_20.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS95, _M0L11dst__offsetS92);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS95, (moonbit_string_t)moonbit_string_literal_21.data);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS95, _M0L3lenS90);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS95, (moonbit_string_t)moonbit_string_literal_22.data);
    _M0L6_2atmpS515 = Moonbit_array_length(_M0L3srcS93);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS95, _M0L6_2atmpS515);
    #line 112 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _M0L6_2atmpS514
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS95);
    moonbit_decref(_M0L18_2astring__builderS95);
    #line 111 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
    _result_1000 = _M0FPC15abort5abortGAkE(_M0L6_2atmpS514);
    moonbit_decref(_M0L6_2atmpS514);
    return _result_1000;
  }
}

uint16_t* _M0MPC15array10FixedArray23unsafe__make__and__blitGkE(
  uint16_t* _M0L3srcS86,
  int32_t _M0L13allocate__lenS83,
  int32_t _M0L4initS84,
  int32_t _M0L11src__offsetS87,
  int32_t _M0L11dst__offsetS85,
  int32_t _M0L9blit__lenS88
) {
  uint16_t* _M0L3dstS82;
  #line 79 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  _M0L3dstS82
  = (uint16_t*)moonbit_make_string(_M0L13allocate__lenS83, _M0L4initS84);
  moonbit_incref(_M0L3dstS82);
  #line 90 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  moonbit_unsafe_val_array_blit(_M0L3dstS82, _M0L11dst__offsetS85, _M0L3srcS86, _M0L11src__offsetS87, _M0L9blit__lenS88, sizeof(uint16_t));
  return _M0L3dstS82;
}

struct _M0TPB13StringBuilder* _M0MPB13StringBuilder21StringBuilder_2einner(
  int32_t _M0L10size__hintS80
) {
  int32_t _M0L7initialS79;
  uint16_t* _M0L4dataS81;
  struct _M0TPB13StringBuilder* _block_1001;
  #line 32 "/home/developer/.moon/lib/core/builtin/stringbuilder_buffer.mbt"
  if (_M0L10size__hintS80 < 1) {
    _M0L7initialS79 = 1;
  } else {
    int32_t _M0L6_2atmpS510 = _M0L10size__hintS80 + 1;
    _M0L7initialS79 = _M0L6_2atmpS510 / 2;
  }
  _M0L4dataS81 = (uint16_t*)moonbit_make_string(_M0L7initialS79, 0);
  _block_1001
  = (struct _M0TPB13StringBuilder*)moonbit_malloc(sizeof(struct _M0TPB13StringBuilder));
  Moonbit_object_header(_block_1001)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 30, 0);
  _block_1001->$0 = _M0L4dataS81;
  _block_1001->$1 = 0;
  return _block_1001;
}

int32_t _M0MPC14byte4Byte8to__char(int32_t _M0L4selfS78) {
  int32_t _M0L6_2atmpS509;
  #line 1868 "/home/developer/.moon/lib/core/builtin/intrinsics.mbt"
  _M0L6_2atmpS509 = (int32_t)_M0L4selfS78;
  return _M0L6_2atmpS509;
}

moonbit_string_t* _M0MPB18UninitializedArray23make__and__blit_2einnerGsE(
  moonbit_string_t* _M0L3srcS70,
  int32_t _M0L13allocate__lenS66,
  int32_t _M0L3lenS67,
  int32_t _M0L11src__offsetS68,
  int32_t _M0L11dst__offsetS69
) {
  int32_t _if__result_1002;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS66 >= 0) {
    if (_M0L3lenS67 >= 0) {
      if (_M0L11src__offsetS68 >= 0) {
        if (_M0L11dst__offsetS69 >= 0) {
          int32_t _M0L6_2atmpS500 = _M0L11src__offsetS68 + _M0L3lenS67;
          int32_t _M0L6_2atmpS501;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS501 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS70);
          if (_M0L6_2atmpS500 <= _M0L6_2atmpS501) {
            int32_t _M0L6_2atmpS499 = _M0L11dst__offsetS69 + _M0L3lenS67;
            _if__result_1002 = _M0L6_2atmpS499 <= _M0L13allocate__lenS66;
          } else {
            _if__result_1002 = 0;
          }
        } else {
          _if__result_1002 = 0;
        }
      } else {
        _if__result_1002 = 0;
      }
    } else {
      _if__result_1002 = 0;
    }
  } else {
    _if__result_1002 = 0;
  }
  if (_if__result_1002) {
    moonbit_incref(_M0L3srcS70);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (moonbit_string_t*)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS66, (moonbit_string_t)moonbit_string_literal_0.data, _M0L3srcS70, _M0L11src__offsetS68, _M0L11dst__offsetS69, _M0L3lenS67);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS71;
    int32_t _M0L6_2atmpS503;
    moonbit_string_t _M0L6_2atmpS502;
    moonbit_string_t* _result_1003;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS71
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS71, (moonbit_string_t)moonbit_string_literal_18.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS71, _M0L13allocate__lenS66);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS71, (moonbit_string_t)moonbit_string_literal_19.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS71, _M0L11src__offsetS68);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS71, (moonbit_string_t)moonbit_string_literal_20.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS71, _M0L11dst__offsetS69);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS71, (moonbit_string_t)moonbit_string_literal_21.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS71, _M0L3lenS67);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS71, (moonbit_string_t)moonbit_string_literal_22.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS503 = _M0MPB18UninitializedArray6lengthGsE(_M0L3srcS70);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS71, _M0L6_2atmpS503);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS502
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS71);
    moonbit_decref(_M0L18_2astring__builderS71);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_1003
    = _M0FPC15abort5abortGRPB18UninitializedArrayGsEE(_M0L6_2atmpS502);
    moonbit_decref(_M0L6_2atmpS502);
    return _result_1003;
  }
}

struct _M0TUsiE** _M0MPB18UninitializedArray23make__and__blit_2einnerGUsiEE(
  struct _M0TUsiE** _M0L3srcS76,
  int32_t _M0L13allocate__lenS72,
  int32_t _M0L3lenS73,
  int32_t _M0L11src__offsetS74,
  int32_t _M0L11dst__offsetS75
) {
  int32_t _if__result_1004;
  #line 168 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  if (_M0L13allocate__lenS72 >= 0) {
    if (_M0L3lenS73 >= 0) {
      if (_M0L11src__offsetS74 >= 0) {
        if (_M0L11dst__offsetS75 >= 0) {
          int32_t _M0L6_2atmpS505 = _M0L11src__offsetS74 + _M0L3lenS73;
          int32_t _M0L6_2atmpS506;
          #line 179 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
          _M0L6_2atmpS506
          = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS76);
          if (_M0L6_2atmpS505 <= _M0L6_2atmpS506) {
            int32_t _M0L6_2atmpS504 = _M0L11dst__offsetS75 + _M0L3lenS73;
            _if__result_1004 = _M0L6_2atmpS504 <= _M0L13allocate__lenS72;
          } else {
            _if__result_1004 = 0;
          }
        } else {
          _if__result_1004 = 0;
        }
      } else {
        _if__result_1004 = 0;
      }
    } else {
      _if__result_1004 = 0;
    }
  } else {
    _if__result_1004 = 0;
  }
  if (_if__result_1004) {
    moonbit_incref(_M0L3srcS76);
    #line 185 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    return (struct _M0TUsiE**)moonbit_make_ref_array_with_blit(_M0L13allocate__lenS72, 0, _M0L3srcS76, _M0L11src__offsetS74, _M0L11dst__offsetS75, _M0L3lenS73);
  } else {
    struct _M0TPB13StringBuilder* _M0L18_2astring__builderS77;
    int32_t _M0L6_2atmpS508;
    moonbit_string_t _M0L6_2atmpS507;
    struct _M0TUsiE** _result_1005;
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L18_2astring__builderS77
    = _M0MPB13StringBuilder21StringBuilder_2einner(89);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS77, (moonbit_string_t)moonbit_string_literal_18.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS77, _M0L13allocate__lenS72);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS77, (moonbit_string_t)moonbit_string_literal_19.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS77, _M0L11src__offsetS74);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS77, (moonbit_string_t)moonbit_string_literal_20.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS77, _M0L11dst__offsetS75);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS77, (moonbit_string_t)moonbit_string_literal_21.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS77, _M0L3lenS73);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0IPB13StringBuilderPB6Logger13write__string(_M0L18_2astring__builderS77, (moonbit_string_t)moonbit_string_literal_22.data);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS508 = _M0MPB18UninitializedArray6lengthGUsiEE(_M0L3srcS76);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0MPB13StringBuilder13write__objectGiE(_M0L18_2astring__builderS77, _M0L6_2atmpS508);
    #line 182 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _M0L6_2atmpS507
    = _M0MPB13StringBuilder10to__string(_M0L18_2astring__builderS77);
    moonbit_decref(_M0L18_2astring__builderS77);
    #line 181 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
    _result_1005
    = _M0FPC15abort5abortGRPB18UninitializedArrayGUsiEEE(_M0L6_2atmpS507);
    moonbit_decref(_M0L6_2atmpS507);
    return _result_1005;
  }
}

int32_t _M0MPB13StringBuilder13write__objectGsE(
  struct _M0TPB13StringBuilder* _M0L4selfS63,
  moonbit_string_t _M0L3objS62
) {
  struct _M0TPB6Logger _M0L6_2atmpS497;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS63);
  _M0L6_2atmpS497
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS63
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS62, _M0L6_2atmpS497);
  if (_M0L6_2atmpS497.$1) {
    moonbit_decref(_M0L6_2atmpS497.$1);
  }
  return 0;
}

int32_t _M0MPB13StringBuilder13write__objectGiE(
  struct _M0TPB13StringBuilder* _M0L4selfS65,
  int32_t _M0L3objS64
) {
  struct _M0TPB6Logger _M0L6_2atmpS498;
  #line 17 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  moonbit_incref(_M0L4selfS65);
  _M0L6_2atmpS498
  = (struct _M0TPB6Logger){
    _M0FP0119moonbitlang_2fcore_2fbuiltin_2fStringBuilder_2eas___40moonbitlang_2fcore_2fbuiltin_2eLogger_2estatic__method__table__id,
      _M0L4selfS65
  };
  #line 23 "/home/developer/.moon/lib/core/builtin/stringbuilder.mbt"
  _M0IP016_24default__implPB4Show6outputGiE(_M0L3objS64, _M0L6_2atmpS498);
  if (_M0L6_2atmpS498.$1) {
    moonbit_decref(_M0L6_2atmpS498.$1);
  }
  return 0;
}

moonbit_string_t* _M0MPB18UninitializedArray23unsafe__make__and__blitGsE(
  moonbit_string_t* _M0L3srcS53,
  int32_t _M0L13allocate__lenS51,
  int32_t _M0L11src__offsetS54,
  int32_t _M0L11dst__offsetS52,
  int32_t _M0L9blit__lenS55
) {
  moonbit_string_t* _M0L3dstS50;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS50
  = (moonbit_string_t*)moonbit_make_ref_array(_M0L13allocate__lenS51, (moonbit_string_t)moonbit_string_literal_0.data);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGsE(_M0L3dstS50, _M0L11dst__offsetS52, _M0L3srcS53, _M0L11src__offsetS54, _M0L9blit__lenS55);
  moonbit_decref(_M0L3srcS53);
  return _M0L3dstS50;
}

struct _M0TUsiE** _M0MPB18UninitializedArray23unsafe__make__and__blitGUsiEE(
  struct _M0TUsiE** _M0L3srcS59,
  int32_t _M0L13allocate__lenS57,
  int32_t _M0L11src__offsetS60,
  int32_t _M0L11dst__offsetS58,
  int32_t _M0L9blit__lenS61
) {
  struct _M0TUsiE** _M0L3dstS56;
  #line 133 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0L3dstS56
  = (struct _M0TUsiE**)moonbit_make_ref_array(_M0L13allocate__lenS57, 0);
  #line 143 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  _M0MPB18UninitializedArray12unsafe__blitGUsiEE(_M0L3dstS56, _M0L11dst__offsetS58, _M0L3srcS59, _M0L11src__offsetS60, _M0L9blit__lenS61);
  moonbit_decref(_M0L3srcS59);
  return _M0L3dstS56;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGsE(
  moonbit_string_t* _M0L3dstS40,
  int32_t _M0L11dst__offsetS41,
  moonbit_string_t* _M0L3srcS42,
  int32_t _M0L11src__offsetS43,
  int32_t _M0L3lenS44
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS42);
  moonbit_incref(_M0L3dstS40);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS40, _M0L11dst__offsetS41, _M0L3srcS42, _M0L11src__offsetS43, _M0L3lenS44);
  return 0;
}

int32_t _M0MPB18UninitializedArray12unsafe__blitGUsiEE(
  struct _M0TUsiE** _M0L3dstS45,
  int32_t _M0L11dst__offsetS46,
  struct _M0TUsiE** _M0L3srcS47,
  int32_t _M0L11src__offsetS48,
  int32_t _M0L3lenS49
) {
  #line 119 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_incref(_M0L3srcS47);
  moonbit_incref(_M0L3dstS45);
  #line 128 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  moonbit_unsafe_ref_array_blit(_M0L3dstS45, _M0L11dst__offsetS46, _M0L3srcS47, _M0L11src__offsetS48, _M0L3lenS49);
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGkE(
  uint16_t* _M0L3dstS13,
  int32_t _M0L11dst__offsetS15,
  uint16_t* _M0L3srcS14,
  int32_t _M0L11src__offsetS16,
  int32_t _M0L3lenS18
) {
  int32_t _if__result_1006;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS13 == _M0L3srcS14) {
    _if__result_1006 = _M0L11dst__offsetS15 < _M0L11src__offsetS16;
  } else {
    _if__result_1006 = 0;
  }
  if (_if__result_1006) {
    int32_t _M0L1iS17 = 0;
    while (1) {
      if (_M0L1iS17 < _M0L3lenS18) {
        int32_t _M0L6_2atmpS470 = _M0L11dst__offsetS15 + _M0L1iS17;
        int32_t _M0L6_2atmpS472 = _M0L11src__offsetS16 + _M0L1iS17;
        int32_t _M0L6_2atmpS471;
        int32_t _M0L6_2atmpS473;
        if (
          _M0L6_2atmpS472 < 0
          || _M0L6_2atmpS472 >= Moonbit_array_length(_M0L3srcS14)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS471 = (int32_t)_M0L3srcS14[_M0L6_2atmpS472];
        if (
          _M0L6_2atmpS470 < 0
          || _M0L6_2atmpS470 >= Moonbit_array_length(_M0L3dstS13)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS13[_M0L6_2atmpS470] = _M0L6_2atmpS471;
        _M0L6_2atmpS473 = _M0L1iS17 + 1;
        _M0L1iS17 = _M0L6_2atmpS473;
        continue;
      } else {
        moonbit_decref(_M0L3srcS14);
        moonbit_decref(_M0L3dstS13);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS478 = _M0L3lenS18 - 1;
    int32_t _M0L1iS20 = _M0L6_2atmpS478;
    while (1) {
      if (_M0L1iS20 >= 0) {
        int32_t _M0L6_2atmpS474 = _M0L11dst__offsetS15 + _M0L1iS20;
        int32_t _M0L6_2atmpS476 = _M0L11src__offsetS16 + _M0L1iS20;
        int32_t _M0L6_2atmpS475;
        int32_t _M0L6_2atmpS477;
        if (
          _M0L6_2atmpS476 < 0
          || _M0L6_2atmpS476 >= Moonbit_array_length(_M0L3srcS14)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS475 = (int32_t)_M0L3srcS14[_M0L6_2atmpS476];
        if (
          _M0L6_2atmpS474 < 0
          || _M0L6_2atmpS474 >= Moonbit_array_length(_M0L3dstS13)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L3dstS13[_M0L6_2atmpS474] = _M0L6_2atmpS475;
        _M0L6_2atmpS477 = _M0L1iS20 - 1;
        _M0L1iS20 = _M0L6_2atmpS477;
        continue;
      } else {
        moonbit_decref(_M0L3srcS14);
        moonbit_decref(_M0L3dstS13);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGsEE(
  moonbit_string_t* _M0L3dstS22,
  int32_t _M0L11dst__offsetS24,
  moonbit_string_t* _M0L3srcS23,
  int32_t _M0L11src__offsetS25,
  int32_t _M0L3lenS27
) {
  int32_t _if__result_1009;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS22 == _M0L3srcS23) {
    _if__result_1009 = _M0L11dst__offsetS24 < _M0L11src__offsetS25;
  } else {
    _if__result_1009 = 0;
  }
  if (_if__result_1009) {
    int32_t _M0L1iS26 = 0;
    while (1) {
      if (_M0L1iS26 < _M0L3lenS27) {
        int32_t _M0L6_2atmpS479 = _M0L11dst__offsetS24 + _M0L1iS26;
        int32_t _M0L6_2atmpS481 = _M0L11src__offsetS25 + _M0L1iS26;
        moonbit_string_t _M0L6_2atmpS480;
        moonbit_string_t _M0L6_2aoldS929;
        int32_t _M0L6_2atmpS482;
        if (
          _M0L6_2atmpS481 < 0
          || _M0L6_2atmpS481 >= Moonbit_array_length(_M0L3srcS23)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS480 = (moonbit_string_t)_M0L3srcS23[_M0L6_2atmpS481];
        if (
          _M0L6_2atmpS479 < 0
          || _M0L6_2atmpS479 >= Moonbit_array_length(_M0L3dstS22)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS929 = (moonbit_string_t)_M0L3dstS22[_M0L6_2atmpS479];
        moonbit_incref(_M0L6_2atmpS480);
        moonbit_decref(_M0L6_2aoldS929);
        _M0L3dstS22[_M0L6_2atmpS479] = _M0L6_2atmpS480;
        _M0L6_2atmpS482 = _M0L1iS26 + 1;
        _M0L1iS26 = _M0L6_2atmpS482;
        continue;
      } else {
        moonbit_decref(_M0L3srcS23);
        moonbit_decref(_M0L3dstS22);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS487 = _M0L3lenS27 - 1;
    int32_t _M0L1iS29 = _M0L6_2atmpS487;
    while (1) {
      if (_M0L1iS29 >= 0) {
        int32_t _M0L6_2atmpS483 = _M0L11dst__offsetS24 + _M0L1iS29;
        int32_t _M0L6_2atmpS485 = _M0L11src__offsetS25 + _M0L1iS29;
        moonbit_string_t _M0L6_2atmpS484;
        moonbit_string_t _M0L6_2aoldS931;
        int32_t _M0L6_2atmpS486;
        if (
          _M0L6_2atmpS485 < 0
          || _M0L6_2atmpS485 >= Moonbit_array_length(_M0L3srcS23)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS484 = (moonbit_string_t)_M0L3srcS23[_M0L6_2atmpS485];
        if (
          _M0L6_2atmpS483 < 0
          || _M0L6_2atmpS483 >= Moonbit_array_length(_M0L3dstS22)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS931 = (moonbit_string_t)_M0L3dstS22[_M0L6_2atmpS483];
        moonbit_incref(_M0L6_2atmpS484);
        moonbit_decref(_M0L6_2aoldS931);
        _M0L3dstS22[_M0L6_2atmpS483] = _M0L6_2atmpS484;
        _M0L6_2atmpS486 = _M0L1iS29 - 1;
        _M0L1iS29 = _M0L6_2atmpS486;
        continue;
      } else {
        moonbit_decref(_M0L3srcS23);
        moonbit_decref(_M0L3dstS22);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPC15array10FixedArray12unsafe__blitGRPB17UnsafeMaybeUninitGUsiEEE(
  struct _M0TUsiE** _M0L3dstS31,
  int32_t _M0L11dst__offsetS33,
  struct _M0TUsiE** _M0L3srcS32,
  int32_t _M0L11src__offsetS34,
  int32_t _M0L3lenS36
) {
  int32_t _if__result_1012;
  #line 38 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
  if (_M0L3dstS31 == _M0L3srcS32) {
    _if__result_1012 = _M0L11dst__offsetS33 < _M0L11src__offsetS34;
  } else {
    _if__result_1012 = 0;
  }
  if (_if__result_1012) {
    int32_t _M0L1iS35 = 0;
    while (1) {
      if (_M0L1iS35 < _M0L3lenS36) {
        int32_t _M0L6_2atmpS488 = _M0L11dst__offsetS33 + _M0L1iS35;
        int32_t _M0L6_2atmpS490 = _M0L11src__offsetS34 + _M0L1iS35;
        struct _M0TUsiE* _M0L6_2atmpS489;
        struct _M0TUsiE* _M0L6_2aoldS933;
        int32_t _M0L6_2atmpS491;
        if (
          _M0L6_2atmpS490 < 0
          || _M0L6_2atmpS490 >= Moonbit_array_length(_M0L3srcS32)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS489 = (struct _M0TUsiE*)_M0L3srcS32[_M0L6_2atmpS490];
        if (
          _M0L6_2atmpS488 < 0
          || _M0L6_2atmpS488 >= Moonbit_array_length(_M0L3dstS31)
        ) {
          #line 50 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS933 = (struct _M0TUsiE*)_M0L3dstS31[_M0L6_2atmpS488];
        if (_M0L6_2atmpS489) {
          moonbit_incref(_M0L6_2atmpS489);
        }
        if (_M0L6_2aoldS933) {
          moonbit_decref(_M0L6_2aoldS933);
        }
        _M0L3dstS31[_M0L6_2atmpS488] = _M0L6_2atmpS489;
        _M0L6_2atmpS491 = _M0L1iS35 + 1;
        _M0L1iS35 = _M0L6_2atmpS491;
        continue;
      } else {
        moonbit_decref(_M0L3srcS32);
        moonbit_decref(_M0L3dstS31);
      }
      break;
    }
  } else {
    int32_t _M0L6_2atmpS496 = _M0L3lenS36 - 1;
    int32_t _M0L1iS38 = _M0L6_2atmpS496;
    while (1) {
      if (_M0L1iS38 >= 0) {
        int32_t _M0L6_2atmpS492 = _M0L11dst__offsetS33 + _M0L1iS38;
        int32_t _M0L6_2atmpS494 = _M0L11src__offsetS34 + _M0L1iS38;
        struct _M0TUsiE* _M0L6_2atmpS493;
        struct _M0TUsiE* _M0L6_2aoldS935;
        int32_t _M0L6_2atmpS495;
        if (
          _M0L6_2atmpS494 < 0
          || _M0L6_2atmpS494 >= Moonbit_array_length(_M0L3srcS32)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2atmpS493 = (struct _M0TUsiE*)_M0L3srcS32[_M0L6_2atmpS494];
        if (
          _M0L6_2atmpS492 < 0
          || _M0L6_2atmpS492 >= Moonbit_array_length(_M0L3dstS31)
        ) {
          #line 54 "/home/developer/.moon/lib/core/builtin/fixedarray_block.mbt"
          moonbit_panic();
        }
        _M0L6_2aoldS935 = (struct _M0TUsiE*)_M0L3dstS31[_M0L6_2atmpS492];
        if (_M0L6_2atmpS493) {
          moonbit_incref(_M0L6_2atmpS493);
        }
        if (_M0L6_2aoldS935) {
          moonbit_decref(_M0L6_2aoldS935);
        }
        _M0L3dstS31[_M0L6_2atmpS492] = _M0L6_2atmpS493;
        _M0L6_2atmpS495 = _M0L1iS38 - 1;
        _M0L1iS38 = _M0L6_2atmpS495;
        continue;
      } else {
        moonbit_decref(_M0L3srcS32);
        moonbit_decref(_M0L3dstS31);
      }
      break;
    }
  }
  return 0;
}

int32_t _M0MPB18UninitializedArray6lengthGsE(moonbit_string_t* _M0L4selfS11) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS11);
}

int32_t _M0MPB18UninitializedArray6lengthGUsiEE(
  struct _M0TUsiE** _M0L4selfS12
) {
  #line 113 "/home/developer/.moon/lib/core/builtin/uninitialized_array.mbt"
  return Moonbit_array_length(_M0L4selfS12);
}

int32_t _M0IPB7FailurePB4Show6output(
  void* _M0L10_2ax__5721S7,
  struct _M0TPB6Logger _M0L10_2ax__5722S10
) {
  struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure* _M0L10_2aFailureS8;
  moonbit_string_t _M0L15_2a_2aarg__5723S9;
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2aFailureS8
  = (struct _M0DTPC15error5Error48moonbitlang_2fcore_2fbuiltin_2eFailure_2eFailure*)_M0L10_2ax__5721S7;
  _M0L15_2a_2aarg__5723S9 = _M0L10_2aFailureS8->$0;
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S10.$0->$method_0(_M0L10_2ax__5722S10.$1, (moonbit_string_t)moonbit_string_literal_23.data);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0MPB6Logger13write__objectGsE(_M0L10_2ax__5722S10, _M0L15_2a_2aarg__5723S9);
  #line 37 "/home/developer/.moon/lib/core/builtin/failure.mbt"
  _M0L10_2ax__5722S10.$0->$method_0(_M0L10_2ax__5722S10.$1, (moonbit_string_t)moonbit_string_literal_24.data);
  return 0;
}

int32_t _M0MPB6Logger13write__objectGsE(
  struct _M0TPB6Logger _M0L4selfS6,
  moonbit_string_t _M0L3objS5
) {
  #line 173 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  #line 174 "/home/developer/.moon/lib/core/builtin/traits.mbt"
  _M0IP016_24default__implPB4Show6outputGsE(_M0L3objS5, _M0L4selfS6);
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

moonbit_string_t _M0FP15Error10to__string(void* _M0L4_2aeS441) {
  switch (Moonbit_object_tag(_M0L4_2aeS441)) {
    case 2: {
      return (moonbit_string_t)moonbit_string_literal_25.data;
      break;
    }
    
    case 0: {
      return _M0IP016_24default__implPB4Show10to__stringGRPB7FailureE(_M0L4_2aeS441);
      break;
    }
    
    case 3: {
      return (moonbit_string_t)moonbit_string_literal_26.data;
      break;
    }
    
    case 4: {
      return (moonbit_string_t)moonbit_string_literal_27.data;
      break;
    }
    default: {
      return (moonbit_string_t)moonbit_string_literal_28.data;
      break;
    }
  }
}

int32_t _M0IP016_24default__implPB6Logger61write_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS464,
  struct _M0TPB4Show _M0L8_2aparamS463
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS462 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS464;
  _M0IP016_24default__implPB6Logger5writeGRPB13StringBuilderE(_M0L7_2aselfS462, _M0L8_2aparamS463);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger84write__string__interpolation_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS461,
  struct _M0TPB4Show _M0L8_2aparamS460
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS459 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS461;
  _M0IP016_24default__implPB6Logger28write__string__interpolationGRPB13StringBuilderE(_M0L7_2aselfS459, _M0L8_2aparamS460);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__char_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS458,
  int32_t _M0L8_2aparamS457
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS456 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS458;
  _M0IPB13StringBuilderPB6Logger11write__char(_M0L7_2aselfS456, _M0L8_2aparamS457);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger67write__view_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS455,
  struct _M0TPC16string10StringView _M0L8_2aparamS454
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS453 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS455;
  _M0IPB13StringBuilderPB6Logger11write__view(_M0L7_2aselfS453, _M0L8_2aparamS454);
  return 0;
}

int32_t _M0IP016_24default__implPB6Logger72write__substring_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLoggerGRPB13StringBuilderE(
  void* _M0L11_2aobj__ptrS452,
  moonbit_string_t _M0L8_2aparamS449,
  int32_t _M0L8_2aparamS450,
  int32_t _M0L8_2aparamS451
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS448 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS452;
  _M0IP016_24default__implPB6Logger16write__substringGRPB13StringBuilderE(_M0L7_2aselfS448, _M0L8_2aparamS449, _M0L8_2aparamS450, _M0L8_2aparamS451);
  return 0;
}

int32_t _M0IPB13StringBuilderPB6Logger69write__string_2edyncall__as___40moonbitlang_2fcore_2fbuiltin_2eLogger(
  void* _M0L11_2aobj__ptrS447,
  moonbit_string_t _M0L8_2aparamS446
) {
  struct _M0TPB13StringBuilder* _M0L7_2aselfS445 =
    (struct _M0TPB13StringBuilder*)_M0L11_2aobj__ptrS447;
  _M0IPB13StringBuilderPB6Logger13write__string(_M0L7_2aselfS445, _M0L8_2aparamS446);
  return 0;
}

void moonbit_init() {
  moonbit_layout_table = moonbit_layout_table_data;
}

int main(int argc, char** argv) {
  struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error** _M0L6_2atmpS469;
  struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE* _M0L12async__testsS435;
  struct _M0TPB5ArrayGUsiEE* _M0L7_2abindS436;
  int32_t _M0L7_2abindS437;
  int32_t _M0L2__S438;
  moonbit_runtime_init(argc, argv);
  moonbit_init();
  _M0L6_2atmpS469
  = (struct _M0TWWuEuWRPC15error5ErrorEuEOuQRPC15error5Error**)moonbit_empty_ref_array;
  _M0L12async__testsS435
  = (struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE*)moonbit_malloc(sizeof(struct _M0TPB5ArrayGVWEuQRPC15error5ErrorE));
  Moonbit_object_header(_M0L12async__testsS435)->meta
  = Moonbit_make_regular_object_header(MOONBIT_REGULAR_LAYOUT_CLASS_INDEXED, 33, 0);
  _M0L12async__testsS435->$0 = _M0L6_2atmpS469;
  _M0L12async__testsS435->$1 = 0;
  #line 399 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0L7_2abindS436
  = _M0FP46biolab8bio__seq3cmd11seqio__main52moonbit__test__driver__internal__native__parse__args();
  _M0L7_2abindS437 = _M0L7_2abindS436->$1;
  _M0L2__S438 = 0;
  while (1) {
    if (_M0L2__S438 < _M0L7_2abindS437) {
      struct _M0TUsiE** _M0L3bufS468 = _M0L7_2abindS436->$0;
      struct _M0TUsiE* _M0L3argS439 =
        (struct _M0TUsiE*)_M0L3bufS468[_M0L2__S438];
      moonbit_string_t _M0L6_2atmpS465 = _M0L3argS439->$0;
      int32_t _M0L6_2atmpS466 = _M0L3argS439->$1;
      int32_t _M0L6_2atmpS467;
      moonbit_incref(_M0L6_2atmpS465);
      #line 400 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
      _M0FP46biolab8bio__seq3cmd11seqio__main44moonbit__test__driver__internal__do__execute(_M0L12async__testsS435, _M0L6_2atmpS465, _M0L6_2atmpS466);
      moonbit_decref(_M0L6_2atmpS465);
      _M0L6_2atmpS467 = _M0L2__S438 + 1;
      _M0L2__S438 = _M0L6_2atmpS467;
      continue;
    } else {
      moonbit_decref(_M0L7_2abindS436);
    }
    break;
  }
  #line 402 "/home/developer/Documents/1/cmd/seqio_main/__generated_driver_for_internal_test.mbt"
  _M0IP016_24default__implP46biolab8bio__seq3cmd11seqio__main28MoonBit__Async__Test__Driver17run__async__testsGRP46biolab8bio__seq3cmd11seqio__main34MoonBit__Async__Test__Driver__ImplE(_M0L12async__testsS435);
  moonbit_decref(_M0L12async__testsS435);
  return 0;
}