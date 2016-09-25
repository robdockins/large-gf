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
  uint16_t a = gf2_16_log_table[x];
  uint16_t b = gf2_16_log_table[y];
  uint16_t z = zeroMask(x) | zeroMask(y);
  return gf2_16_expadd( z, a, b );
}

// Compute polynomial division in GF(2^16).
// In other words, find q such that
//      mult(q, y) = x  (mod irreducible)
//
// Precondition: y != 0
inline uint16_t gf2_16_div( uint16_t x, uint16_t y ) {
  uint16_t a = gf2_16_log_table[x];
  uint16_t b = gf2_16_log_table[y];
  uint16_t z = zeroMask(x);
  return gf2_16_expadd( z, a, Q - b );
}

// precondition x != 0
inline uint16_t gf2_16_inv( uint16_t x ) {
  uint16_t a = gf2_16_log_table[x];
  uint16_t b = Q - a;
  return gf2_16_exp_table[b];
}

inline uint16_t gf2_16_exp( uint16_t x ) {
  return gf2_16_exp_table[x];
}

inline uint16_t gf2_16_log( uint16_t x ) {
  return gf2_16_log_table[x];
}

// We want to calculate a^x, for 0 <= x < 2^16, and a != 0.
//
// Let Q = 2^16-1, which is the order of the mutiplicative
// group of GF(2^16). It is a fact of this multiplicative group
// that a^Q = 1.  Now, let l = log(a), then
// a^x = exp(log(a))^x = exp(l)^x = exp( l*x ).
// Moreover, because of the previous fact,
// this is equal to exp( (l*x)%Q ).
//
// Finally, it is another amazing fact of modular arithmetic
// that we can calculate (z%Q) very easily when z
// is a 32-bit number: we can simply add the high
// 16-bits to the low-16 bits.  This result may carry into
// a 17-bit number.  This is fine, as our exponential table
// is defined for twice the field size, i.e., for 17-bit indices.
inline uint16_t gf2_16_pow( uint16_t a, uint16_t x ) {
  uint32_t log;

  // First calculate the multiplied logarithm.
  log = (uint32_t) gf2_16_log_table[a];
  log = log * (uint32_t) x;

  // Now take the modulus by (2^16-1).  The following
  // works by some magic of modular arithmetic.
  log = (log >> 16) + (log & 0xFFFF);

  return gf2_16_exp_table[ log ];
}

#endif
