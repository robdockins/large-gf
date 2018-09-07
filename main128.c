#include <stdlib.h>
#include <stdio.h>

#include "gf2_16.h"
#include "gf2_128.h"

const unsigned long MAX_ROUNDS = 100000000;

int main(void) {
  unsigned int randreg;
  FILE* urand = fopen("/dev/urandom", "r");
  fread( &randreg, sizeof(randreg), 1, urand );
  fclose(urand);
  //randreg = 0x1650338b;
  printf( "randreg = %#.8x\n", randreg );

  uint64_t i;
  for( i = 0; i < MAX_ROUNDS; i++ ) {
    uint128_t x, y, z;

    x = (uint128_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint128_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint128_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint128_t) rand_r( &randreg );

    y = (uint128_t) rand_r( &randreg );
    y <<= 32;
    y |= (uint128_t) rand_r( &randreg );
    y <<= 32;
    y |= (uint128_t) rand_r( &randreg );
    y <<= 32;
    y |= (uint128_t) rand_r( &randreg );

    z = (uint128_t) rand_r( &randreg );
    z <<= 32;
    z |= (uint128_t) rand_r( &randreg );
    z <<= 32;
    z |= (uint128_t) rand_r( &randreg );
    z <<= 32;
    z |= (uint128_t) rand_r( &randreg );

    // Associativity test
    if( gf2_128_mult( x, gf2_128_mult( y, z ) ) ==
        gf2_128_mult( gf2_128_mult( x, y ), z ) ) {
    } else {
      printf( "Associativity fail:\n" );
      printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
      printf( "  y = 0x%.16lx%.16lx\n", (uint64_t) (y>>64), (uint64_t) y );
      printf( "  z = 0x%.16lx%.16lx\n", (uint64_t) (z>>64), (uint64_t) z );
    }

  }

  printf( "%lld Tests completed\n", MAX_ROUNDS );
  return 0;
}
