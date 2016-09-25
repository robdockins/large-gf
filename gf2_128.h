#ifndef GF2_128_H
#define GF2_128_H

#include <stdint.h>

typedef unsigned __int128 uint128_t;

inline uint128_t gf2_128_add( uint128_t x, uint128_t y) {
  return (x ^ y);
}

uint128_t gf2_128_mult( uint128_t, uint128_t );
//uint128_t gf2_128_square( uint128_t );
//uint128_t gf2_128_inv( uint128_t );
//uint128_t gf2_128_square16( uint128_t );

#endif
