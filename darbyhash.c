/* 
 * This is a modification of Ley Yureiv's t1ha hash https://github.com/leo-yuriev/t1ha
 *
 * The modifications include assuming AVX and AES instruction sets as well as optimizing it for large strings
 */
#include "darbyhash.h"

/* assume unaligned reads are ok. address sanitizers and valgrind don't like this but it will work for most CPUs */
#define DH_UNALIGNED_OK 1

/* 'magic' primes */
static const uint64_t p0 = 17048867929148541611ull;
static const uint64_t p1 = 9386433910765580089ull;
static const uint64_t p2 = 15343884574428479051ull;
static const uint64_t p3 = 13662985319504319857ull;
static const uint64_t p4 = 11242949449147999147ull;
static const uint64_t p5 = 13862205317416547141ull;
static const uint64_t p6 = 14653293970879851569ull;

/* rotations */
static const unsigned s0 = 41;
static const unsigned s1 = 17;

static __inline uint64_t dh_rot64(uint64_t v, unsigned s) {
  return (v >> s) | (v << (64 - s));
}

static __inline uint64_t dh_mix(uint64_t v, uint64_t p) {
  v *= p;
  return v ^ dh_rot64(v, s0);
}

/* xor high and low parts of full 128-bit product */
static __inline uint64_t dh_mux64(uint64_t v, uint64_t p) {
  __uint128_t r = (__uint128_t)v * (__uint128_t)p;
  return ((uint64_t)r) ^ (r >> 64);
}

static __inline uint64_t fetch64_le(const void *v) {
  return *(const uint64_t *)v;
}

static __inline uint32_t fetch32_le(const void *v) {
  return *(const uint32_t *)v;
}

static __inline uint16_t fetch16_le(const void *v) {
  return *(const uint16_t *)v;
}

static __inline uint64_t dh_tail64_le(const void *v, size_t tail) {
  const uint8_t *p = (const uint8_t *)v;
  uint64_t r = 0;
  switch (tail & 7) {
#if DH_UNALIGNED_OK
  /* For most CPUs this code is better when not needed
   * copying for alignment or byte reordering.
   */
  case 0:
    return fetch64_le(p);
  case 7:
    r = (uint64_t)p[6] << 8;
  case 6:
    r += p[5];
    r <<= 8;
  case 5:
    r += p[4];
    r <<= 32;
  case 4:
    return r + fetch32_le(p);
  case 3:
    r = (uint64_t)p[2] << 16;
  case 2:
    return r + fetch16_le(p);
  case 1:
    return p[0];
#else
  /* For most CPUs this code is better than a
   * copying for alignment and/or byte reordering. */
  case 0:
    r = p[7] << 8;
  case 7:
    r += p[6];
    r <<= 8;
  case 6:
    r += p[5];
    r <<= 8;
  case 5:
    r += p[4];
    r <<= 8;
  case 4:
    r += p[3];
    r <<= 8;
  case 3:
    r += p[2];
    r <<= 8;
  case 2:
    r += p[1];
    r <<= 8;
  case 1:
    return r + p[0];
#endif
  }
  // not reachable
  return 0;
}

uint64_t darbyhash_aes_avx(const void *data, size_t len, uint64_t seed) {
  uint64_t a = seed;
  uint64_t b = len;

  if (len > 32) {
    __m128i x = _mm_set_epi64x(a, b);
    __m128i y = _mm_aesenc_si128(x, _mm_set_epi64x(p0, p1));

    const __m128i *v = (const __m128i *)data;
    const __m128i *const detent =
        (const __m128i *)((const uint8_t *)data + (len & ~15ul));
    data = detent;

    if (len & 16) {
      x = _mm_add_epi64(x, _mm_loadu_si128(v++));
      y = _mm_aesenc_si128(x, y);
    }
    len &= 15;

    if (v + 7 < detent) {
      __m128i salt = y;
      do {
        __m128i t = _mm_aesenc_si128(_mm_loadu_si128(v+0), salt);
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+1));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+2));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+3));

        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+4));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+5));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+6));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v+7));

        salt = _mm_add_epi64(salt, _mm_set_epi64x(p2, p3));
        t = _mm_aesenc_si128(x, t);
        x = _mm_add_epi64(y, x);
        y = t;
        v+=8;
      } while (v + 7 < detent);
    }

    __m128i v0y;
    __m128i v1x;
    switch (detent - v)
    {
    case 8:
    case 7:
      v0y = _mm_add_epi64(y, _mm_loadu_si128(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si128(v++));
      x = _mm_aesdec_si128(x, v0y);
      y = _mm_aesdec_si128(y, v1x);
    case 6:
    case 5:
      v0y = _mm_add_epi64(y, _mm_loadu_si128(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si128(v++));
      x = _mm_aesdec_si128(x, v0y);
      y = _mm_aesdec_si128(y, v1x);
    case 4:
    case 3:
      v0y = _mm_add_epi64(y, _mm_loadu_si128(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si128(v++));
      x = _mm_aesdec_si128(x, v0y);
      y = _mm_aesdec_si128(y, v1x);
    case 2:
    case 1:
      v0y = _mm_add_epi64(y, _mm_loadu_si128(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si128(v++));
      x = _mm_aesdec_si128(x, v0y);
      y = _mm_aesdec_si128(y, v1x);
    }

    x = _mm_add_epi64(_mm_aesdec_si128(x, _mm_aesenc_si128(y, x)), y);
    a = _mm_cvtsi128_si64(x);
    b = _mm_extract_epi64(x, 1);
  }

  /* Finish last 32 bytes with a more standard multipicative hash function */
  const uint64_t *v = (const uint64_t *)data;
  switch (len) {
  default:
    b += dh_mux64(*v++, p4);
  /* fall through */
  case 24:
  case 23:
  case 22:
  case 21:
  case 20:
  case 19:
  case 18:
  case 17:
    a += dh_mux64(*v++, p3);
  /* fall through */
  case 16:
  case 15:
  case 14:
  case 13:
  case 12:
  case 11:
  case 10:
  case 9:
    b += dh_mux64(*v++, p2);
  /* fall through */
  case 8:
  case 7:
  case 6:
  case 5:
  case 4:
  case 3:
  case 2:
  case 1:
    a += dh_mux64(dh_tail64_le(v, len), p1);
  /* fall through */
  case 0:
    return dh_mux64(dh_rot64(a + b, s1), p4) + dh_mix(a ^ b, p0);
  }
}
uint64_t darbyhash(const void *data, size_t len, uint64_t seed) {
  uint64_t a = seed;
  uint64_t b = len;

  if (len > 32) {
    __m256i x = _mm256_set_epi64x(0, 0, a, b);
    __m256i y = _mm_aesenc_si256(x, _mm256_set_epi64x(p0, p1, p5, p6));

    const __m256i *v = (const __m256i *)data;
    const __m256i *const detent =
        (const __m256i *)((const uint8_t *)data + (len & ~15ul));
    data = detent;

    if (len & 16) {
      x = _mm_add_epi64(x, _mm_loadu_si256(v++));
      y = _mm_aesenc_si256(x, y);
    }
    len &= 15;

    if (v + 7 < detent) {
      __m256i salt = y;
      do {
        __m256i t = _mm_aesenc_si256(_mm_loadu_si256(v+0), salt);
        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+1));
        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+2));
        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+3));

        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+4));
        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+5));
        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+6));
        t = _mm_aesdec_si256(t, _mm_loadu_si256(v+7));

        salt = _mm_add_epi64(salt, _mm256_set_epi64x(p2, p3, p4, p0));
        t = _mm_aesenc_si256(x, t);
        x = _mm_add_epi64(y, x);
        y = t;
        v+=8;
      } while (v + 7 < detent);
    }

    __m256i v0y;
    __m256i v1x;
    switch (detent - v)
    {
    case 8:
    case 7:
      v0y = _mm_add_epi64(y, _mm_loadu_si256(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si256(v++));
      x = _mm_aesdec_si256(x, v0y);
      y = _mm_aesdec_si256(y, v1x);
    case 6:
    case 5:
      v0y = _mm_add_epi64(y, _mm_loadu_si256(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si256(v++));
      x = _mm_aesdec_si256(x, v0y);
      y = _mm_aesdec_si256(y, v1x);
    case 4:
    case 3:
      v0y = _mm_add_epi64(y, _mm_loadu_si256(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si256(v++));
      x = _mm_aesdec_si256(x, v0y);
      y = _mm_aesdec_si256(y, v1x);
    case 2:
    case 1:
      v0y = _mm_add_epi64(y, _mm_loadu_si256(v++));
      v1x = _mm_sub_epi64(x, _mm_loadu_si256(v++));
      x = _mm_aesdec_si256(x, v0y);
      y = _mm_aesdec_si256(y, v1x);
    }

    x = _mm_add_epi64(_mm_aesdec_si256(x, _mm_aesenc_si256(y, x)), y);
    a = _mm_cvtsi256_si64(x);
    b = _mm_extract_epi64(x, 1);
  }

  /* Finish last 32 bytes with a more standard multipicative hash function */
  const uint64_t *v = (const uint64_t *)data;
  switch (len) {
  default:
    b += dh_mux64(*v++, p4);
  /* fall through */
  case 24:
  case 23:
  case 22:
  case 21:
  case 20:
  case 19:
  case 18:
  case 17:
    a += dh_mux64(*v++, p3);
  /* fall through */
  case 16:
  case 15:
  case 14:
  case 13:
  case 12:
  case 11:
  case 10:
  case 9:
    b += dh_mux64(*v++, p2);
  /* fall through */
  case 8:
  case 7:
  case 6:
  case 5:
  case 4:
  case 3:
  case 2:
  case 1:
    a += dh_mux64(dh_tail64_le(v, len), p1);
  /* fall through */
  case 0:
    return dh_mux64(dh_rot64(a + b, s1), p4) + dh_mix(a ^ b, p0);
  }
}

uint64_t darbyhash_avx(const void *data, size_t len, uint64_t seed) {
  uint64_t a = seed;
  uint64_t b = len;

  if (len > 32) {
    __m128i x = _mm_set_epi64x(a, b);
    __m128i y = _mm_aesenc_si128(x, _mm_set_epi64x(p0, p1));

    const __m128i *v = (const __m128i *)data;
    const __m128i *const detent =
        (const __m128i *)((const uint8_t *)data + (len & ~15ul));
    data = detent;

    if (len & 16) {
      x = _mm_add_epi64(x, _mm_loadu_si128(v++));
      y = _mm_aesenc_si128(x, y);
    }
    len &= 15;

    if (v + 7 < detent) {
      __m128i salt = y;
      do {
        __m128i t = _mm_aesenc_si128(_mm_loadu_si128(v++), salt);
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));

        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));
        t = _mm_aesdec_si128(t, _mm_loadu_si128(v++));

        salt = _mm_add_epi64(salt, _mm_set_epi64x(p2, p3));
        t = _mm_aesenc_si128(x, t);
        x = _mm_add_epi64(y, x);
        y = t;
      } while (v + 7 < detent);
    }

    while (v < detent) {
      __m128i v0y = _mm_add_epi64(y, _mm_loadu_si128(v++));
      __m128i v1x = _mm_sub_epi64(x, _mm_loadu_si128(v++));
      x = _mm_aesdec_si128(x, v0y);
      y = _mm_aesdec_si128(y, v1x);
    }

    x = _mm_add_epi64(_mm_aesdec_si128(x, _mm_aesenc_si128(y, x)), y);
    a = _mm_cvtsi128_si64(x);
    b = _mm_extract_epi64(x, 1);
  }

  const uint64_t *v = (const uint64_t *)data;
  switch (len) {
  default:
    b += dh_mux64(*v++, p4);
  case 24:
  case 23:
  case 22:
  case 21:
  case 20:
  case 19:
  case 18:
  case 17:
    a += dh_mux64(*v++, p3);
  case 16:
  case 15:
  case 14:
  case 13:
  case 12:
  case 11:
  case 10:
  case 9:
    b += dh_mux64(*v++, p2);
  case 8:
  case 7:
  case 6:
  case 5:
  case 4:
  case 3:
  case 2:
  case 1:
    a += dh_mux64(dh_tail64_le(v, len), p1);
  case 0:
    return dh_mux64(dh_rot64(a + b, s1), p4) + dh_mix(a ^ b, p0);
  }
}
