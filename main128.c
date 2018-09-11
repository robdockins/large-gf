#include <stdlib.h>
#include <stdio.h>

#include "gf2_16.h"
#include "gf2_128.h"

//const unsigned long MAX_ROUNDS = 100000000;
const unsigned long MAX_ROUNDS = 10000000;
//const unsigned long MAX_ROUNDS = 100;

int main(void) {
  unsigned int randreg;
  FILE* urand = fopen("/dev/urandom", "r");
  fread( &randreg, sizeof(randreg), 1, urand );
  fclose(urand);
  //randreg = 0x1650338b;
  printf( "randreg = %#.8x\n", randreg );

  uint128_t x, y, z, w, v;

  uint64_t i;
  for( i = 0; i < MAX_ROUNDS; i++ ) {

    x = (uint128_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint128_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint128_t) rand_r( &randreg );
    x <<= 32;
    x |= (uint128_t) rand_r( &randreg );

    /* y = (uint128_t) rand_r( &randreg ); */
    /* y <<= 32; */
    /* y |= (uint128_t) rand_r( &randreg ); */
    /* y <<= 32; */
    /* y |= (uint128_t) rand_r( &randreg ); */
    /* y <<= 32; */
    /* y |= (uint128_t) rand_r( &randreg ); */

    /* z = (uint128_t) rand_r( &randreg ); */
    /* z <<= 32; */
    /* z |= (uint128_t) rand_r( &randreg ); */
    /* z <<= 32; */
    /* z |= (uint128_t) rand_r( &randreg ); */
    /* z <<= 32; */
    /* z |= (uint128_t) rand_r( &randreg ); */

#if 0
    // Commutivity test
    w = gf2_128_mult( x, y );
    v = gf2_128_mult( y, x );
    if( w == v ) {
    } else {
      printf( "Commutivity fail:\n" );
      printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
      printf( "  y = 0x%.16lx%.16lx\n", (uint64_t) (y>>64), (uint64_t) y );
      printf( "  w = 0x%.16lx%.16lx\n", (uint64_t) (w>>64), (uint64_t) w );
      printf( "  v = 0x%.16lx%.16lx\n", (uint64_t) (v>>64), (uint64_t) v );
    }

    // Associativity test
    w = gf2_128_mult( x, gf2_128_mult( y, z ) );
    v = gf2_128_mult( gf2_128_mult( x, y ), z );
    if( w == v ) {
    } else {
      printf( "Associativity fail:\n" );
      printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
      printf( "  y = 0x%.16lx%.16lx\n", (uint64_t) (y>>64), (uint64_t) y );
      printf( "  z = 0x%.16lx%.16lx\n", (uint64_t) (z>>64), (uint64_t) z );
      printf( "  w = 0x%.16lx%.16lx\n", (uint64_t) (w>>64), (uint64_t) w );
      printf( "  v = 0x%.16lx%.16lx\n", (uint64_t) (v>>64), (uint64_t) v );
    }

    // Distrbutivity test
    w = gf2_128_mult( x, y ^ z);
    v = gf2_128_mult( x, y ) ^ gf2_128_mult( x, z );
    if( w == v ) {
    } else {
      printf( "Distributivity fail:\n" );
      printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
      printf( "  y = 0x%.16lx%.16lx\n", (uint64_t) (y>>64), (uint64_t) y );
      printf( "  z = 0x%.16lx%.16lx\n", (uint64_t) (z>>64), (uint64_t) z );
      printf( "  w = 0x%.16lx%.16lx\n", (uint64_t) (w>>64), (uint64_t) w );
      printf( "  v = 0x%.16lx%.16lx\n", (uint64_t) (v>>64), (uint64_t) v );
    }

    // Square test
    w = gf2_128_mult( x, x );
    v = gf2_128_square( x );
    if( w == v ) {
    } else {
      printf( "Square fail:\n" );
      printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
      printf( "  w = 0x%.16lx%.16lx\n", (uint64_t) (w>>64), (uint64_t) w );
      printf( "  v = 0x%.16lx%.16lx\n", (uint64_t) (v>>64), (uint64_t) v );
    }

    // Square^16 test
    w = x;
    for (int j=0; j<16; j++) {
      w = gf2_128_square( w );
    }
    v = gf2_128_square16( x );
    if( w == v ) {
    } else {
      printf( "Square16 fail:\n" );
      printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
      printf( "  w = 0x%.16lx%.16lx\n", (uint64_t) (w>>64), (uint64_t) w );
      printf( "  v = 0x%.16lx%.16lx\n", (uint64_t) (v>>64), (uint64_t) v );
    }
#endif

    // inverse test
    if( x != 0 ) {
      w = gf2_128_inv( x );
      v = gf2_128_mult( x, w );
      if( v == 0x1 ) {
      } else {
        printf( "Inversion fail:\n" );
        printf( "  x = 0x%.16lx%.16lx\n", (uint64_t) (x>>64), (uint64_t) x );
        printf( "  w = 0x%.16lx%.16lx\n", (uint64_t) (w>>64), (uint64_t) w );
        printf( "  v = 0x%.16lx%.16lx\n", (uint64_t) (v>>64), (uint64_t) v );
      }
    }

  }

  printf( "%lld Tests completed\n", MAX_ROUNDS );

  return 0;
}
