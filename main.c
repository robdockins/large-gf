#include <stdlib.h>
#include <stdio.h>

#include "gf2_16.h"
#include "gf2_32.h"

//const unsigned long MAX_ROUNDS = 1000000000;
const unsigned long MAX_ROUNDS = 100000000;

int main() {
  //printf( "gf2_16_log(8192) = 0x%.4x\n", gf2_16_log_table[8192] );

  unsigned int randreg;
  FILE* urand = fopen("/dev/urandom", "r");
  fread( &randreg, sizeof(randreg), 1, urand );
  fclose(urand);
  printf( "randreg = %#.8x\n", randreg );

  unsigned long i;
  for( i = 0; i < MAX_ROUNDS; i++ ) {
    uint32_t x = (uint32_t) rand_r( &randreg );
    uint32_t y = (uint32_t) rand_r( &randreg );
    uint32_t z = (uint32_t) rand_r( &randreg );

    // Distrbutivity test
    if( gf2_32_mult( x, y ^ z) ==
        (gf2_32_mult( x, y ) ^ gf2_32_mult( x, z )) ) {
    } else {
      printf( "Distributivity fail:\n" );
      printf( "  x = 0x%.8x\n", x );
      printf( "  y = 0x%.8x\n", y );
      printf( "  z = 0x%.8x\n", z );
    }

    // Associativity test
    if( gf2_32_mult( x, gf2_32_mult( y, z ) ) ==
        gf2_32_mult( gf2_32_mult( x, y ), z ) ) {
    } else {
      printf( "Associativity fail:\n" );
      printf( "  x = 0x%.8x\n", x );
      printf( "  y = 0x%.8x\n", y );
      printf( "  z = 0x%.8x\n", z );
    }

    // Commutitivity test
    if( gf2_32_mult( x, y ) == gf2_32_mult( y, x ) ) {
    } else {
      printf( "Commutivity fail:\n" );
      printf( "  x = 0x%.8x\n", x );
      printf( "  y = 0x%.8x\n", y );
    }

    // squaring test
    if( gf2_32_mult( x, x ) == gf2_32_square( x ) ) {
    } else {
      printf( "Squaring fail:\n" );
      printf( "  x = 0x%.8x\n", x );
    }

    // inverse test
    if( x == 0 || (gf2_32_mult( x, gf2_32_inv( x ) ) == 1) ) {
    } else {
      printf( "Inversion fail:\n" );
      printf( "  x = 0x%.8x\n", x );
    }
  }

  printf( "%d Tests completed\n", MAX_ROUNDS );
}


void print16logtable()
{
  for( int i = 0; i < FIELD_SIZE; i+=16 ) {
    printf( "  , 0x%.4x , 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x, 0x%.4x\n",
            gf2_16_log_table[i+0],
            gf2_16_log_table[i+1],
            gf2_16_log_table[i+2],
            gf2_16_log_table[i+3],
            gf2_16_log_table[i+4],
            gf2_16_log_table[i+5],
            gf2_16_log_table[i+6],
            gf2_16_log_table[i+7],
            gf2_16_log_table[i+8],
            gf2_16_log_table[i+9],
            gf2_16_log_table[i+10],
            gf2_16_log_table[i+11],
            gf2_16_log_table[i+12],
            gf2_16_log_table[i+13],
            gf2_16_log_table[i+14],
            gf2_16_log_table[i+15] );
  }
}
