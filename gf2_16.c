#include <stdint.h>

#include "gf2_16.h"

uint16_t gf2_16_log_table[FIELD_SIZE];
uint16_t gf2_16_exp_table[2*FIELD_SIZE];

inline uint16_t nextPower( uint16_t b ) {
  return
    (b >> 15) ?
      ((b << 1) ^ gf2_16_irreducible) :
      (b << 1);
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
