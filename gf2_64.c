#include <stdint.h>

#include "gf2_16.h"
#include "gf2_64.h"

// In this file, we construct the Galois field GF(2^64)
// as a finite field extension of GF(2^16). The
// field is constructed using the following irreducible
// polynomial, whose coefficents should be understood
// as the binary representation of elements of GF(2^16):
//
//    x^4 + x^2 + 2x + 1
//


inline uint64_t gf2_64_mult( uint64_t a, uint64_t b ) {
  uint16_t a0, a1, a2, a3;
  uint16_t za0, za1, za2, za3;
  uint32_t alog0, alog1, alog2, alog3;

  uint16_t b0, b1, b2, b3;
  uint16_t zb0, zb1, zb2, zb3;
  uint32_t blog0, blog1, blog2, blog3;

  uint16_t t, z;
  uint16_t c0, c1, c2, c3, c4, c5, c6;
  uint16_t d0, d1, d2, d3;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  a2 = (uint16_t) (a >> 32);
  a3 = (uint16_t) (a >> 48);

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = gf2_16_log_table[a0];
  alog1 = gf2_16_log_table[a1];
  alog2 = gf2_16_log_table[a2];
  alog3 = gf2_16_log_table[a3];

  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);
  b2 = (uint16_t) (b >> 32);
  b3 = (uint16_t) (b >> 48);

  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );
  zb2 = zeroMask( b2 );
  zb3 = zeroMask( b3 );

  blog0 = gf2_16_log_table[b0];
  blog1 = gf2_16_log_table[b1];
  blog2 = gf2_16_log_table[b2];
  blog3 = gf2_16_log_table[b3];

  // Do the raw coefficent multiplications
  // via GF(2^16) tables

  t   = gf2_16_exp_table[ alog0 + blog0 ];
  z   = za0 | zb0;
  c0  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog0 + blog1 ];
  z   = za0 | zb1;
  c1  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog1 + blog0 ];
  z   = za1 | zb0;
  c1 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog1 + blog1 ];
  z   = za1 | zb1;
  c2  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog2 + blog0 ];
  z   = za2 | zb0;
  c2 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog0 + blog2 ];
  z   = za0 | zb2;
  c2 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog1 + blog2 ];
  z   = za1 | zb2;
  c3  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog2 + blog1 ];
  z   = za2 | zb1;
  c3 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + blog0 ];
  z   = za3 | zb0;
  c3 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog0 + blog3 ];
  z   = za0 | zb3;
  c3 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog1 + blog3 ];
  z   = za1 | zb3;
  c4  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog2 + blog2 ];
  z   = za2 | zb2;
  c4 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + blog1 ];
  z   = za3 | zb1;
  c4 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog2 + blog3 ];
  z   = za2 | zb3;
  c5  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + blog2 ];
  z   = za3 | zb2;
  c5 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + blog3 ];
  z   = za3 | zb3;
  c6  = (z & t) ^ t;

  // Now perform the modular reduction
  uint16_t c4x2 = gf2_16_mult( c4, 2 );
  uint16_t c5x2 = gf2_16_mult( c5, 2 );
  uint16_t c6x2 = gf2_16_mult( c6, 2 );

  d3  = c3 ^ c5 ^ c6x2;
  d2  = c2 ^ c4 ^ c5x2;
  d1  = c1 ^ c4x2 ^ c5 ^ c6x2;
  d0  = c0 ^ c4 ^ c6;

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}

// Multiply a*b, but only calculate the low coefficent.
inline uint16_t gf2_64_mult_low_coeff( uint64_t a, uint64_t b ) {
  uint16_t a0, a1, a2, a3;
  uint16_t za0, za1, za2, za3;
  uint32_t alog0, alog1, alog2, alog3;

  uint16_t b0, b1, b2, b3;
  uint16_t zb0, zb1, zb2, zb3;
  uint32_t blog0, blog1, blog2, blog3;

  uint16_t t, z;
  uint16_t c0, c4, c6;
  uint16_t d0;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  a2 = (uint16_t) (a >> 32);
  a3 = (uint16_t) (a >> 48);

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = gf2_16_log_table[a0];
  alog1 = gf2_16_log_table[a1];
  alog2 = gf2_16_log_table[a2];
  alog3 = gf2_16_log_table[a3];

  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);
  b2 = (uint16_t) (b >> 32);
  b3 = (uint16_t) (b >> 48);

  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );
  zb2 = zeroMask( b2 );
  zb3 = zeroMask( b3 );

  blog0 = gf2_16_log_table[b0];
  blog1 = gf2_16_log_table[b1];
  blog2 = gf2_16_log_table[b2];
  blog3 = gf2_16_log_table[b3];

  // Do the raw coefficent multiplications
  // via GF(2^16) tables

  t   = gf2_16_exp_table[ alog0 + blog0 ];
  z   = za0 | zb0;
  c0  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog1 + blog3 ];
  z   = za1 | zb3;
  c4  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog2 + blog2 ];
  z   = za2 | zb2;
  c4 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + blog1 ];
  z   = za3 | zb1;
  c4 ^= (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + blog3 ];
  z   = za3 | zb3;
  c6  = (z & t) ^ t;

  d0  = c0 ^ c4 ^ c6;

  return d0;
}

// This rather magical-looking function computes
// the value a^(2^16).  Because the field has characteristic
// two, squaring commutes with addition, as does any number
// of iterations of squaring.  Moreover, raising the coefficents
// (which lie in GF(2^16)) to the 2^16 power is the identity
// function.  Thus, the overall effect of raising to the 2^16
// power is a linear transformation of the original coefficents.
//
// The calculations below implement multiplication by the following
// matrix, which represents the action of the linear transformation.
//
//  0x0001  0x0000  0x0000  0x0000
//  0x0fd3  0x9b04  0x1f3f  0x0000
//  0xc870  0x393c  0x9b04  0x0000
//  0x29d2  0x6d0b  0x36ef  0x0001
//
inline uint64_t gf2_64_square16( uint64_t a ) {
  uint16_t a0, a1, a2, a3;
  uint16_t d0, d1, d2, d3;

  a3 = (uint16_t) (a >> 48);
  a2 = (uint16_t) (a >> 32);
  a1 = (uint16_t) (a >> 16);
  a0 = (uint16_t) a;

  d3 = a3;
  d2 = gf2_16_mult( 0x0fd3, a3 ) ^
       gf2_16_mult( 0x9b04, a2 ) ^
       gf2_16_mult( 0x1f3f, a1 );
  d1 = gf2_16_mult( 0xc870, a3 ) ^
       gf2_16_mult( 0x393c, a2 ) ^
       gf2_16_mult( 0x9b04, a1 );
  d0 = gf2_16_mult( 0x29d2, a3 ) ^
       gf2_16_mult( 0x6d0b, a2 ) ^
       gf2_16_mult( 0x36ef, a1 ) ^
       a0;

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}


inline uint64_t gf2_64_square( uint64_t a ) {
  uint16_t a0, a1, a2, a3;
  uint16_t za0, za1, za2, za3;
  uint32_t alog0, alog1, alog2, alog3;

  uint16_t t, z;
  uint16_t c0, c1, c2, c3, c4, c5, c6;
  uint16_t d0, d1, d2, d3;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  a2 = (uint16_t) (a >> 32);
  a3 = (uint16_t) (a >> 48);

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = gf2_16_log_table[a0];
  alog1 = gf2_16_log_table[a1];
  alog2 = gf2_16_log_table[a2];
  alog3 = gf2_16_log_table[a3];

  // Do the raw coefficent multiplications
  // via GF(2^16) tables

  t   = gf2_16_exp_table[ alog0 + alog0 ];
  z   = za0;
  c0  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog1 + alog1 ];
  z   = za1;
  c2  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog2 + alog2 ];
  z   = za2;
  c4  = (z & t) ^ t;

  t   = gf2_16_exp_table[ alog3 + alog3 ];
  z   = za3;
  c6  = (z & t) ^ t;

  // Now perform the modular reduction
  uint16_t c4x2 = gf2_16_mult( c4, 2 );
  uint16_t c6x2 = gf2_16_mult( c6, 2 );

  d3  = c6x2;
  d2  = c2 ^ c4;
  d1  = c4x2 ^ c6x2;
  d0  = c0 ^ c4 ^ c6;

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}

inline uint64_t gf2_64_pointwise_mult( uint16_t x, uint64_t a ) {
  uint16_t a0, a1, a2, a3;
  uint16_t d0, d1, d2, d3;

  a3 = (uint16_t) (a >> 48);
  a2 = (uint16_t) (a >> 32);
  a1 = (uint16_t) (a >> 16);
  a0 = (uint16_t) a;

  d3 = gf2_16_mult( x, a3 );
  d2 = gf2_16_mult( x, a2 );
  d1 = gf2_16_mult( x, a1 );
  d0 = gf2_16_mult( x, a0 );

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}

uint64_t gf2_64_inv( uint64_t a ) {
  // Let q = 2^16.  Let r = (q^4 - 1)/(q - 1) = 2^48 + 2^32 + 2^16 + 1
  // This rather special number has to do with caluclating finite field norms.
  // For all x in GF((2^16)^4), x^r yields a value in GF(2^16); that is, for which
  // the high coefficents are 0.  We exploit this fact to perform fast inversions
  // in GF(2^64) by reducing them to inversion in GF(2^16).
  //
  // We are going to calculate a^(-1) = a^(-r) * a^(r-1).  The algorithm below goes
  // roughly as follows:
  //
  //   s = a^(r-1)
  //     = a^(q^3 + q^2 + q)
  //     = a^(q(q^2 + q + 1))
  //     = a^(q(q(q+1)+1))
  //   t = a * s = a^r
  //   b = t^(-1) * s = a^(-1)

  // Compute s = a^(r-1)
  uint64_t s;
  s = gf2_64_square16( a );
  s = gf2_64_mult( s, a );
  s = gf2_64_square16( s );
  s = gf2_64_mult( s, a );
  s = gf2_64_square16( s );

  // t = s * a = a^r
  // Because we know t in GF(2^16), we can save
  // some work by only calculating the low coefficent
  // of this multiplication.
  uint16_t t0 = gf2_64_mult_low_coeff( s, a );

  // Now invert t0
  uint16_t t0_inv = gf2_16_inv( t0 );

  // b = t^(-1) * s = a^(-r) * a^(r-1) = a(-1)
  uint64_t b = gf2_64_pointwise_mult( t0_inv, s );
  return b;
}
