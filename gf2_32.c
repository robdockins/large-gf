#include <stdint.h>
#include "gf2_16.h"
#include "gf2_32.h"

#include <stdio.h>


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

  z = zeroMask( c2 );
  t = gf2_16_exp_table[ gf2_16_log_table[ c2 ] + vlog ];
  c0 ^= (z & t) ^ t;

  uint32_t res = (((uint32_t) c1) << 16) | ((uint32_t) c0);
  return res;
}
