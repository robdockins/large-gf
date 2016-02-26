#include <stdlib.h>
#include <stdio.h>

#include "gf2_16.h"
#include "gf2_64.h"

//const unsigned long MAX_ROUNDS = 1000000000;
const unsigned long MAX_ROUNDS = 10000000;

int main() {
  unsigned int randreg;
  FILE* urand = fopen("/dev/urandom", "r");
  fread( &randreg, sizeof(randreg), 1, urand );
  fclose(urand);
  printf( "randreg = %#.8x\n", randreg );

  unsigned long i;
  for( i = 0; i < MAX_ROUNDS; i++ ) {
    uint64_t x, y, z;

    x = (uint64_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint64_t) rand_r( &randreg );

    y = (uint64_t) rand_r( &randreg );
    y <<= 32;
    y |= (uint64_t) rand_r( &randreg );

    z = (uint64_t) rand_r( &randreg );
    z <<= 32;
    z |= (uint64_t) rand_r( &randreg );

    // Distrbutivity test
    if( gf2_64_mult( x, y ^ z ) ==
        gf2_64_mult( x, y ) ^ gf2_64_mult( x, z ) ) {
    } else {
      printf( "Distributivity fail:\n" );
      printf( "  x = 0x%.16lx\n", x );
      printf( "  y = 0x%.16lx\n", y );
      printf( "  z = 0x%.16lx\n", z );
    }

    // Associativity test
    if( gf2_64_mult( x, gf2_64_mult( y, z ) ) ==
        gf2_64_mult( gf2_64_mult( x, y ), z ) ) {
    } else {
      printf( "Associativity fail:\n" );
      printf( "  x = 0x%.16lx\n", x );
      printf( "  y = 0x%.16lx\n", y );
      printf( "  z = 0x%.16lx\n", z );
    }

    // Commutitivity test
    if( gf2_64_mult( x, y ) == gf2_64_mult( y, x ) ) {
    } else {
      printf( "Commutivity fail:\n" );
      printf( "  x = 0x%.16lx\n", x );
      printf( "  y = 0x%.16lx\n", y );
    }

    // squaring test
    if( gf2_64_mult( x, x ) == gf2_64_square( x ) ) {
    } else {
      printf( "Squaring fail:\n" );
      printf( "  x = 0x%.16lx\n", x );
    }

    // inverse test
    if( x == 0 || (gf2_64_mult( x, gf2_64_inv( x ) ) == 1) ) {
    } else {
      printf( "Inversion fail:\n" );
      printf( "  x = 0x%.16lx\n", x );
    }
  }

  printf( "%d Tests completed\n", MAX_ROUNDS );
}
