#ifndef GF2_16_H
#define GF2_16_H

#include <stdint.h>

#define FIELD_WIDTH 16
#define FIELD_SIZE (1 << FIELD_WIDTH)
#define Q (FIELD_SIZE - 1)
#define gf2_16_irreducible ((1 << 12) | (1<<3) | (1<<1) | (1<<0))

extern uint16_t gf2_16_log_table[FIELD_SIZE];
extern uint16_t gf2_16_exp_table[2*FIELD_SIZE + 0x0f];

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

// Compute polynomial multiplication in GF(2^16)
inline uint16_t gf2_16_mult( uint16_t x, uint16_t y ) {
  uint32_t a = (uint32_t) gf2_16_log_table[x];
  uint32_t b = (uint32_t) gf2_16_log_table[y];
  uint32_t c = a + b;
  uint16_t d = gf2_16_exp_table[c];
  uint16_t z = zeroMask(x) | zeroMask(y);
  return (z & d) ^ d;
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
