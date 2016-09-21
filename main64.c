#include <stdlib.h>
#include <stdio.h>

#include "gf2_16.h"
#include "gf2_64.h"

//const unsigned long MAX_ROUNDS = 1000000000;
const uint64_t MAX_ROUNDS = 100000000;

const uint64_t w64_prim = 0b11011L;

uint64_t advance( uint64_t x ) {
  if( x & (1L << 63) ) {
    return (x << 1) ^ w64_prim;
  } else {
    return x << 1;
  }
}



int main_asdf() {
  uint64_t x = 0x1;
  uint64_t y = 0x1;

  for(int i=0; i<128; i++) {
    uint64_t fx = gf2_64_iso(x);
    if( fx == y ) {
    } else {
      printf( "Failed at %d\n", i );
      printf( "  x   = 0x%.16llx\n", x );
      printf( "  y   = 0x%.16llx\n", y );
      printf( " f(x) = 0x%.16llx\n", fx );
    }
    
    x = advance(x);
    y = gf2_64_mult( y, gf2_64_generator );
  }

  return 0;
}


/*   for(uint64_t i=0;i<64; i++) { */
/*     uint64_t x = gf2_64_pow( gf2_64_generator, i ); */
/*     printf( "  x = 0x%.16llx\n", x ); */
/*   } */
/* } */

#if 1
int main() {
  unsigned int randreg;
  FILE* urand = fopen("/dev/urandom", "r");
  fread( &randreg, sizeof(randreg), 1, urand );
  fclose(urand);
  printf( "randreg = %#.8x\n", randreg );

  uint64_t i;
  for( i = 0; i < MAX_ROUNDS; i++ ) {
    uint64_t x, y, z;

    x = (uint64_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint64_t) rand_r( &randreg );

    //y = (uint64_t) rand_r( &randreg );
    //y <<= 32;
    //y |= (uint64_t) rand_r( &randreg );

    //z = (uint64_t) rand_r( &randreg );
    //z <<= 32;
    //z |= (uint64_t) rand_r( &randreg );

    z  = gf2_64_pow( x, 64 );
    z ^= gf2_64_pow( x,  4 );
    z ^= gf2_64_pow( x,  3 );
    z ^= x;

    if( z == 1 ) {
      printf( "Isomorphism generator?:\n" );
      printf( "  z = 0x%.16llx\n", z );
    }

    /*
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
    */

  }

  printf( "%lld Tests completed\n", MAX_ROUNDS );
  return 0;
}
#endif
