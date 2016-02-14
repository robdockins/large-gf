#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define FIELD_WIDTH 16
#define FIELD_SIZE (1 << FIELD_WIDTH)
#define Q (FIELD_SIZE - 1)

uint16_t gf2_16_log_table[FIELD_SIZE];
uint16_t gf2_16_exp_table[2*FIELD_SIZE];

// Binary representation of the irreducible polynomial:
//   x^16 + x^12 + x^3 + x + 1
//
// The leading x^16 coefficent is implicit and elided.
const uint16_t irreducible = (1 << 12) | (1<<3) | (1<<1) | (1<<0);

inline uint16_t nextPower( uint16_t b ) {
  return
    (b >> 15) ?
      ((b << 1) ^ irreducible) :
      (b << 1);
}

// This function computes, in a branch-free way, the value:
//    (x == 0) ? 0xFFFF : 0x0000
inline uint16_t zeroMask( uint16_t x ) {
  uint32_t y = (uint32_t) x;
  uint32_t z = (y - 1) >> 16;
  return (uint16_t) z;
}

void __attribute__ ((constructor)) init_tables() {
  // avoid multiple init
  if( gf2_16_exp_table[0] != 0 ) return;

  unsigned int i = 0;
  uint16_t b = 1;

  for( ; i < Q; i++ ) {
    gf2_16_exp_table[i] = b;
    gf2_16_log_table[b] = i;
    b = nextPower( b );
  }

  for( ; i < 2*FIELD_SIZE; i++ ) {
    gf2_16_exp_table[i] = b;
    b = nextPower( b );
  }
}

// Compute polynomial addition in GF(2^16)
uint16_t gf2_16_add( uint16_t x, uint16_t y ) {
  return (x ^ y);
}

// Compute polynomial multiplication in GF(2^16)
uint16_t gf2_16_mult( uint16_t x, uint16_t y ) {
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
uint16_t gf2_16_div( uint16_t x, uint16_t y ) {
  uint32_t a = (uint32_t) gf2_16_log_table[x];
  uint32_t b = (uint32_t) gf2_16_log_table[y];
  uint32_t c = a + Q - b;
  uint16_t d = gf2_16_exp_table[c];
  uint16_t z = zeroMask(x);
  return (z & d) ^ d;
}

uint16_t gf2_16_exp( uint16_t x ) {
  return gf2_16_exp_table[x];
}

uint16_t gf2_16_log( uint16_t x ) {
  return gf2_16_log_table[x];
}

int main() {
  printf("main start\n");
  
  for( int i=0; i<20; i++ ) {
    printf( "EXP[%d] = 0x%.4x\n", i, gf2_16_exp_table[i] );
  }
  for( int i=0; i<20; i++ ) {
    printf( "LOG[%d] = 0x%.4x\n", i, gf2_16_log_table[i] );
  }

  uint16_t x, y;
  x = 0x1234;
  y = 0xabcd;
  uint16_t res = gf2_16_mult( x, y );

  printf( "mult( 0x%.4x, 0x%.4x ) = 0x%.4x\n", x, y, res );

  return 0;
}
  
