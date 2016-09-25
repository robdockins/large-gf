#include <stdlib.h>
#include <stdio.h>

#include "gf_complete.h"

const unsigned long MAX_ROUNDS = 100000000;

int main() {
  unsigned int randreg;
  //FILE* urand = fopen("/dev/urandom", "r");
  //fread( &randreg, sizeof(randreg), 1, urand );
  //fclose(urand);
  randreg = 0x1650338b;
  printf( "randreg = %#.8x\n", randreg );

  gf_t gf;
  gf_init_easy(&gf, 64);

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

    /* // Distrbutivity test */
    /* if( gf.multiply.w64( &gf, x, y ^ z ) == */
    /*     gf.multiply.w64( &gf, x, y ) ^ gf.multiply.w64( &gf, x, z ) ) { */
    /* } else { */
    /*   printf( "Distributivity fail:\n" ); */
    /*   printf( "  x = 0x%.16lx\n", x ); */
    /*   printf( "  y = 0x%.16lx\n", y ); */
    /*   printf( "  z = 0x%.16lx\n", z ); */
    /* } */

    /* // Associativity test */
    /* if( gf.multiply.w64( &gf, x, gf.multiply.w64( &gf, y, z ) ) == */
    /*     gf.multiply.w64( &gf, gf.multiply.w64( &gf, x, y ), z ) ) { */
    /* } else { */
    /*   printf( "Associativity fail:\n" ); */
    /*   printf( "  x = 0x%.16lx\n", x ); */
    /*   printf( "  y = 0x%.16lx\n", y ); */
    /*   printf( "  z = 0x%.16lx\n", z ); */
    /* } */

    // inverse test
    if( x == 0 || gf.inverse.w64( &gf, x ) != 0) {
    } else {
      printf( "Inversion fail:\n" );
      printf( "  x = 0x%.16lx\n", x );
    }

    /* // inverse test */
    /* if( x == 0 || (gf.multiply.w64( &gf, x, gf.inverse.w64( &gf, x ) ) == 1) ) { */
    /* } else { */
    /*   printf( "Inversion fail:\n" ); */
    /*   printf( "  x = 0x%.16lx\n", x ); */
    /* } */

  }

  printf( "%d Tests completed\n", MAX_ROUNDS );
 
}
