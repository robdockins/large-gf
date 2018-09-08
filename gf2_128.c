#include <stdint.h>

#include "gf2_16.h"
#include "gf2_128.h"

// In this file, we construct the Galois field GF(2^128)
// as a finite field extension of GF(2^16). The
// field is constructed using the following irreducible
// polynomial, whose coefficents should be understood
// as the binary representation of elements of GF(2^16):
//
//    x^8 + x^3 + x + 8
//

inline uint128_t gf2_128_mult( uint128_t a, uint128_t b ) {
  uint16_t a0, a1, a2, a3, a4, a5, a6, a7;
  uint16_t b0, b1, b2, b3, b4, b5, b6, b7;
  uint16_t za0, za1, za2, za3, za4, za5, za6, za7;
  uint16_t zb0, zb1, zb2, zb3, zb4, zb5, zb6, zb7;

  uint32_t alog0, alog1, alog2, alog3;
  uint32_t alog4, alog5, alog6, alog7;
  uint32_t blog0, blog1, blog2, blog3;
  uint32_t blog4, blog5, blog6, blog7;

  uint16_t c0, c1, c2, c3, c4, c5, c6, c7;
  uint16_t c8, c9, c10, c11, c12, c13, c14;

  uint16_t d0, d1, d2, d3, d4, d5, d6, d7;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  a2 = (uint16_t) (a >> 32);
  a3 = (uint16_t) (a >> 48);
  a4 = (uint16_t) (a >> 64);
  a5 = (uint16_t) (a >> 80);
  a6 = (uint16_t) (a >> 96);
  a7 = (uint16_t) (a >> 112);

  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);
  b2 = (uint16_t) (b >> 32);
  b3 = (uint16_t) (b >> 48);
  b4 = (uint16_t) (b >> 64);
  b5 = (uint16_t) (b >> 80);
  b6 = (uint16_t) (b >> 96);
  b7 = (uint16_t) (b >> 112);

  alog0 = gf2_16_log_table[a0];
  alog1 = gf2_16_log_table[a1];
  alog2 = gf2_16_log_table[a2];
  alog3 = gf2_16_log_table[a3];
  alog4 = gf2_16_log_table[a4];
  alog5 = gf2_16_log_table[a5];
  alog6 = gf2_16_log_table[a6];
  alog7 = gf2_16_log_table[a7];

  blog0 = gf2_16_log_table[b0];
  blog1 = gf2_16_log_table[b1];
  blog2 = gf2_16_log_table[b2];
  blog3 = gf2_16_log_table[b3];
  blog4 = gf2_16_log_table[b4];
  blog5 = gf2_16_log_table[b5];
  blog6 = gf2_16_log_table[b6];
  blog7 = gf2_16_log_table[b7];

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );
  za4 = zeroMask( a4 );
  za5 = zeroMask( a5 );
  za6 = zeroMask( a6 );
  za7 = zeroMask( a7 );

  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );
  zb2 = zeroMask( b2 );
  zb3 = zeroMask( b3 );
  zb4 = zeroMask( b4 );
  zb5 = zeroMask( b5 );
  zb6 = zeroMask( b6 );
  zb7 = zeroMask( b7 );

  c0  = gf2_16_expadd( za0|zb0, alog0, blog0 );

  c1  = gf2_16_expadd( za0|zb1, alog0, blog1 );
  c1 ^= gf2_16_expadd( za1|zb0, alog1, blog0 );

  c2  = gf2_16_expadd( za0|zb2, alog0, blog2 );
  c2 ^= gf2_16_expadd( za1|zb1, alog1, blog1 );
  c2 ^= gf2_16_expadd( za2|zb0, alog2, blog0 );

  c3  = gf2_16_expadd( za0|zb3, alog0, blog3 );
  c3 ^= gf2_16_expadd( za1|zb2, alog1, blog2 );
  c3 ^= gf2_16_expadd( za2|zb1, alog2, blog1 );
  c3 ^= gf2_16_expadd( za3|zb0, alog3, blog0 );

  c4  = gf2_16_expadd( za0|zb4, alog0, blog4 );
  c4 ^= gf2_16_expadd( za1|zb3, alog1, blog3 );
  c4 ^= gf2_16_expadd( za2|zb2, alog2, blog2 );
  c4 ^= gf2_16_expadd( za3|zb1, alog3, blog1 );
  c4 ^= gf2_16_expadd( za4|zb0, alog4, blog0 );

  c5  = gf2_16_expadd( za0|zb5, alog0, blog5 );
  c5 ^= gf2_16_expadd( za1|zb4, alog1, blog4 );
  c5 ^= gf2_16_expadd( za2|zb3, alog2, blog3 );
  c5 ^= gf2_16_expadd( za3|zb2, alog3, blog2 );
  c5 ^= gf2_16_expadd( za4|zb1, alog4, blog1 );
  c5 ^= gf2_16_expadd( za5|zb0, alog5, blog0 );

  c6  = gf2_16_expadd( za0|zb6, alog0, blog6 );
  c6 ^= gf2_16_expadd( za1|zb5, alog1, blog5 );
  c6 ^= gf2_16_expadd( za2|zb4, alog2, blog4 );
  c6 ^= gf2_16_expadd( za3|zb3, alog3, blog3 );
  c6 ^= gf2_16_expadd( za4|zb2, alog4, blog2 );
  c6 ^= gf2_16_expadd( za5|zb1, alog5, blog1 );
  c6 ^= gf2_16_expadd( za6|zb0, alog6, blog0 );

  c7  = gf2_16_expadd( za0|zb7, alog0, blog7 );
  c7 ^= gf2_16_expadd( za1|zb6, alog1, blog6 );
  c7 ^= gf2_16_expadd( za2|zb5, alog2, blog5 );
  c7 ^= gf2_16_expadd( za3|zb4, alog3, blog4 );
  c7 ^= gf2_16_expadd( za4|zb3, alog4, blog3 );
  c7 ^= gf2_16_expadd( za5|zb2, alog5, blog2 );
  c7 ^= gf2_16_expadd( za6|zb1, alog6, blog1 );
  c7 ^= gf2_16_expadd( za7|zb0, alog7, blog0 );

  c8  = gf2_16_expadd( za1|zb7, alog1, blog7 );
  c8 ^= gf2_16_expadd( za2|zb6, alog2, blog6 );
  c8 ^= gf2_16_expadd( za3|zb5, alog3, blog5 );
  c8 ^= gf2_16_expadd( za4|zb4, alog4, blog4 );
  c8 ^= gf2_16_expadd( za5|zb3, alog5, blog3 );
  c8 ^= gf2_16_expadd( za6|zb2, alog6, blog2 );
  c8 ^= gf2_16_expadd( za7|zb1, alog7, blog1 );

  c9  = gf2_16_expadd( za2|zb7, alog2, blog7 );
  c9 ^= gf2_16_expadd( za3|zb6, alog3, blog6 );
  c9 ^= gf2_16_expadd( za4|zb5, alog4, blog5 );
  c9 ^= gf2_16_expadd( za5|zb4, alog5, blog4 );
  c9 ^= gf2_16_expadd( za6|zb3, alog6, blog3 );
  c9 ^= gf2_16_expadd( za7|zb2, alog7, blog2 );

  c10 = gf2_16_expadd( za3|zb7, alog3, blog7 );
  c10^= gf2_16_expadd( za4|zb6, alog4, blog6 );
  c10^= gf2_16_expadd( za5|zb5, alog5, blog5 );
  c10^= gf2_16_expadd( za6|zb4, alog6, blog4 );
  c10^= gf2_16_expadd( za7|zb3, alog7, blog3 );

  c11 = gf2_16_expadd( za4|zb7, alog4, blog7 );
  c11^= gf2_16_expadd( za5|zb6, alog5, blog6 );
  c11^= gf2_16_expadd( za6|zb5, alog6, blog5 );
  c11^= gf2_16_expadd( za7|zb4, alog7, blog4 );

  c12 = gf2_16_expadd( za5|zb7, alog5, blog7 );
  c12^= gf2_16_expadd( za6|zb6, alog6, blog6 );
  c12^= gf2_16_expadd( za7|zb5, alog7, blog5 );

  c13 = gf2_16_expadd( za6|zb7, alog6, blog7 );
  c13^= gf2_16_expadd( za7|zb6, alog7, blog6 );

  c14 = gf2_16_expadd( za7|zb7, alog7, blog7 );

  // Now, modular reduction
  uint16_t log8  = 3; // gf2_16_log_table[8];

  uint16_t c14x8 = gf2_16_expadd( zeroMask(c14), gf2_16_log_table[c14], log8 );
  uint16_t c13x8 = gf2_16_expadd( zeroMask(c13), gf2_16_log_table[c13], log8 );
  uint16_t c12x8 = gf2_16_expadd( zeroMask(c12), gf2_16_log_table[c12], log8 );
  uint16_t c11x8 = gf2_16_expadd( zeroMask(c11), gf2_16_log_table[c11], log8 );
  uint16_t c10x8 = gf2_16_expadd( zeroMask(c10), gf2_16_log_table[c10], log8 );
  uint16_t  c9x8 = gf2_16_expadd( zeroMask( c9), gf2_16_log_table[ c9], log8 );
  uint16_t  c8x8 = gf2_16_expadd( zeroMask( c8), gf2_16_log_table[ c8], log8 );

  d7 = c14   ^ c12   ^ c7;
  d6 = c14x8 ^ c13   ^ c11  ^ c6;
  d5 = c13x8 ^ c12   ^ c10  ^ c5;
  d4 = c14   ^ c12x8 ^ c11  ^ c9 ^ c4;
  d3 = c13   ^ c11x8 ^ c10  ^ c8 ^ c3;
  d2 = c14   ^ c10x8 ^ c9   ^ c2;
  d1 = c14x8 ^ c13   ^ c9x8 ^ c8  ^ c1;
  d0 = c13x8 ^ c8x8  ^ c0;

  uint128_t d =
    (((uint128_t) d7) << 112) |
    (((uint128_t) d6) <<  96) |
    (((uint128_t) d5) <<  80) |
    (((uint128_t) d4) <<  64) |
    (((uint128_t) d3) <<  48) |
    (((uint128_t) d2) <<  32) |
    (((uint128_t) d1) <<  16) |
     ((uint128_t) d0);

  return d;
}
