#include <stdint.h>

#include "gf2_16.h"
#include "gf2_32.h"


#include <stdio.h>

// In this file we construct the Galois field GF(2^32)
// as a finite field extension of GF(2^16).  The
// field is constructed using the following irreducible
// polynomial, whose coefficents should be understood
// as the binary representation of elements of GF(2^16):
//
//    x^2 + x + 8192
//
// Because the constant term is the only coefficent not
// equal to 0 or 1, modular reduction is fairly fast, requiring
// only one additional multipliation in GF(2^16).
//
// Moreover, we get very lucky with this polynomial because
// the Itoh-Tsujii method for calculating multiplicative inverses
// works especially well for this polynomial.  The main task of
// raising to the r = 2^16 power turns out to be almost trivial.


uint32_t gf2_32_mult( uint32_t a, uint32_t b ) {
  uint16_t a0, a1, b0, b1, za0, za1, zb0, zb1;
  uint32_t alog0, alog1, blog0, blog1;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);

  alog0 = gf2_16_log_table[a0];
  alog1 = gf2_16_log_table[a1];
  blog0 = gf2_16_log_table[b0];
  blog1 = gf2_16_log_table[b1];

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );

  uint16_t c0, c1, c2;
  uint16_t z, t;

  t = gf2_16_exp_table[ alog0 + blog0 ];
  z = za0 | zb0;
  c0 = (z & t) ^ t;

  t = gf2_16_exp_table[ alog0 + blog1 ];
  z = za0 | zb1;
  c1 = (z & t) ^ t;

  t = gf2_16_exp_table[ alog1 + blog0 ];
  z = za1 | zb0;
  c1 ^= (z & t) ^ t;

  t = gf2_16_exp_table[ alog1 + blog1 ];
  z = za1 | zb1;
  c2 = (z & t) ^ t;

  // Now reduce by the field extension irreducible
  // polynomial x^2 + x + 8192

  //uint16_t vlog = gf2_16_log_table[8192];
  // this magic number is equal to gf2_16_log_table[8192];
  const uint16_t vlog = 0x000d;

  c1 ^= c2;

  //z = zeroMask( c2 );
  //t = gf2_16_exp_table[ gf2_16_log_table[ c2 ] + vlog ];
  t = gf2_16_exp_table[ alog1 + blog1 + vlog ];
  c0 ^= (z & t) ^ t;

  // polynomial is now reduced, output the result
  uint32_t res = (((uint32_t) c1) << 16) | ((uint32_t) c0);
  return res;
}

uint32_t gf2_32_square( uint32_t a ) {
  uint16_t a0, a1, b0, b2, c1, c0, t, z;
  uint32_t alog0, alog1;

  // a1(x) + a0
  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);

  alog0 = gf2_16_log_table[a0];
  alog1 = gf2_16_log_table[a1];

  // (a1(x) + a0)^2 = b2(x^2) + b0
  z = zeroMask( a0 );
  t = gf2_16_exp_table[ alog0 + alog0 ];
  b0 = (z & t) ^ t;

  z = zeroMask( a1 );
  t = gf2_16_exp_table[ alog1 + alog1 ];
  b2 = (z & t) ^ t;

  // Now reduce by the field extension irreducible
  // polynomial x^2 + x + 8192

  //uint16_t vlog = gf2_16_log_table[8192];
  // this magic number is equal to gf2_16_log_table[8192];
  const uint16_t vlog = 0x000d;

  c1 = b2;

  //z = zeroMask( b2 );
  //t = gf2_16_exp_table[ gf2_16_log_table[ b2 ] + vlog ];
  t = gf2_16_exp_table[ alog1 + alog1 + vlog ];


  c0 = b0;
  c0 ^= (z & t) ^ t;

  uint32_t res = (((uint32_t) c1) << 16) | ((uint32_t) c0);
  return res;
}

uint32_t gf2_32_inv( uint32_t a ) {
  // Let q = 2^16.  let r = (q^2 - 1) / (q - 1) = 2^16 + 1
  // This rather special number has to do with caluclating finite field norms.
  // For all x in GF((2^16)^2), x^r yields a value in GF(2^16); that is, for which
  // the high coefficent is 0.  We exploit this fact to perform fast inversions
  // in GF(2^32) by reducing them to inversion in GF(2^16).
  //
  // We are going to calculate a^(-1) = a^(-r) * a^(r-1).  The algorithm below goes
  // roughly as follows:
  //
  //   s = a^(r-1) = a^(2^16)
  //   t = a * s = a^r
  //   b = t^(-1) * s = a^(-1)

  uint16_t s0, s1;
  uint16_t t0, t1;

  // split a into high and low coefficents
  uint16_t a0 = (uint16_t) a;
  uint16_t a1 = (uint16_t) (a >> 16);

  // This is some crazy magic! We calculate s = a^(r-1) = a^(2^16) by the
  //  following formula, which seems too insanely simple to be correct,
  //  but somehow works out.
  s1 = a1;
  s0 = a0 ^ a1;

  // Now we calculate t = s * a = a^r.  This is guaranteed, by a theorem
  // of finite fields, to be a constant polynomial, i.e., the high coefficent
  // of t must be equal to 0.  We exploit this to do only the minimal amount
  // of calculation necessary to find the low coefficent of t.  This saves
  // several GF(2^16) multiplies and adds.
  t0 = gf2_16_mult( s0, a0 );

  // Multiply the high coefficents and perform polynomial reduction all at once...
  const uint16_t vlog = 0x000d;
  uint32_t logsum = ((uint32_t) gf2_16_log_table[s1]) + ((uint32_t) gf2_16_log_table[a1] + vlog);
  t1 = gf2_16_exp_table[ logsum ];
  //t1 = gf2_16_exp_table[ ((uint32_t) gf2_16_log_table[ t1 ]) + vlog ];
  uint16_t z = zeroMask( s1 ) | zeroMask( a1 );
  t1 = (z & t1) ^ t1;
  t0 ^= t1;

  // Now, we can find the inverse of t by just taking the inverse of t0 in
  // GF(2^16), which we can easily do via our lookup tables.  However, we
  // can save a few table lookups by instead just directly taking the log of t0,
  // and using that to perform the necessary divisions below.
  //uint16_t t0_inv = gf2_16_inv( t0 );
  uint32_t t0_inv_log = Q - gf2_16_log_table[ t0 ];

  // Now we calculate b = t^(-1) * s = a^(-r) * a^(r-1) = a^(-1).  Because the
  // high coefficent of t is 0, we simply multiply the inverse of t0 by the coefficents
  // of s.  We do this by using the t0_log value we calculated just above.
  // We know t0 != 0, so we don't have to mask for that case.

  z = zeroMask( s0 );
  uint16_t b0 = gf2_16_exp_table[ ((uint32_t) gf2_16_log_table[s0]) + t0_inv_log ];
  b0 = (b0 & z) ^ b0;
  //uint16_t b0 = gf2_16_mult( s0, t0_inv );

  z = zeroMask( s1 );
  uint16_t b1 = gf2_16_exp_table[ ((uint32_t) gf2_16_log_table[s1]) + t0_inv_log ];
  b1 = (b1 & z) ^ b1;
  //uint16_t b1 = gf2_16_mult( s1, t0_inv );

  // Package the result and ship it
  uint32_t b = (((uint32_t) b1) << 16) | ((uint32_t) b0);
  return b;
}

// r = (2^32 - 1) / (2^16 - 1) = 65537 = 2^16 + 1
// r-1 = 65536 = 2^16
