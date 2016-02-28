#ifndef GF2_16_H
#define GF2_16_H

#include <stdint.h>

#define FIELD_WIDTH 16
#define FIELD_SIZE (1 << FIELD_WIDTH)
#define Q (FIELD_SIZE - 1)
#define gf2_16_irreducible ((1 << 12) | (1<<3) | (1<<1) | (1<<0))

// irreducible polynomial is
//   x^17 + x^12 + x^3 + x + 1

extern uint16_t gf2_16_log_table[FIELD_SIZE];
extern uint16_t gf2_16_exp_table[2*FIELD_SIZE];

// This function computes, in a branch-free way, the value:
//    (x == 0) ? 0xFFFF : 0x0000
inline uint16_t zeroMask( uint16_t x ) {
  uint32_t y = (uint32_t) x;
  uint32_t z = (y - 1) >> 16;
  return (uint16_t) z;
}

// Compute polynomial addition in GF(2^16)
inline uint16_t gf2_16_add( uint16_t x, uint16_t y ) {
  return (x ^ y);
}

inline uint16_t logsum_modQ( uint16_t a, uint16_t b ) {
  uint32_t sum = ((uint32_t) a) + ((uint32_t) b);
  uint16_t d = ((uint16_t) sum) + ((uint16_t) (sum >> 16));
  return d;
}

// Compute exp( (a + b)%Q )
inline uint16_t gf2_16_expadd( uint16_t zmask, uint16_t a, uint16_t b ) {
  uint32_t sum = ((uint32_t) a) + ((uint32_t) b);
  uint16_t t = gf2_16_exp_table[ sum ];
  return (zmask & t) ^ t;
}

// Compute exp( (a + b + c)%Q )
inline uint16_t gf2_16_expadd3( uint16_t zmask, uint16_t a, uint16_t b, uint16_t c ) {
  uint32_t sum = ((uint32_t) logsum_modQ( a, b )) + ((uint32_t) c);
  uint16_t t = gf2_16_exp_table[ sum ];
  return (zmask & t) ^ t;
}

// Compute polynomial multiplication in GF(2^16)
inline uint16_t gf2_16_mult( uint16_t x, uint16_t y ) {
  uint32_t a = (uint32_t) gf2_16_log_table[x];
  uint32_t b = (uint32_t) gf2_16_log_table[y];
  uint16_t z = zeroMask(x) | zeroMask(y);
  return gf2_16_expadd( z, a, b );
}

// Compute polynomial division in GF(2^16).
// In other words, find q such that
//      mult(q, y) = x  (mod irreducible)
//
// Precondition: y != 0
inline uint16_t gf2_16_div( uint16_t x, uint16_t y ) {
  uint32_t a = (uint32_t) gf2_16_log_table[x];
  uint32_t b = (uint32_t) gf2_16_log_table[y];
  uint32_t c = a + Q - b;
  uint16_t d = gf2_16_exp_table[c];
  uint16_t z = zeroMask(x);
  return (z & d) ^ d;
}

// precondition x != 0
inline uint16_t gf2_16_inv( uint16_t x ) {
  uint32_t a = (uint32_t) gf2_16_log_table[x];
  uint32_t b = Q - a;
  return gf2_16_exp_table[b];
}

inline uint16_t gf2_16_exp( uint16_t x ) {
  return gf2_16_exp_table[x];
}

inline uint16_t gf2_16_log( uint16_t x ) {
  return gf2_16_log_table[x];
}


#endif
