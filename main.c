#include <stdlib.h>
#include <stdio.h>

#include "gf2_16.h"
#include "gf2_32.h"

const unsigned long MAX_ROUNDS = 10000000;

int main() {
  unsigned int randreg = 0xa1b2c3d4U;
    //0x12345678U;

  unsigned long i;
  for( i = 0; i < MAX_ROUNDS; i++ ) {
    uint32_t x = (uint32_t) rand_r( &randreg );
    uint32_t y = (uint32_t) rand_r( &randreg );
    uint32_t z = (uint32_t) rand_r( &randreg );

    uint32_t ans1 = gf2_32_mult( x, y ^ z);
    uint32_t ans2 = gf2_32_mult( x, y ) ^ gf2_32_mult( x, z );

    /*
    printf("Test values:");
    printf( "  x = 0x%.8x", x );
    printf( "  y = 0x%.8x", y );
    printf( "  z = 0x%.8x\n", z );
    printf( "  ans1 = 0x%.8x\n", ans1 );
    printf( "  ans2 = 0x%.8x\n", ans2 );
*/

    // Distrbutivity test
    if( ans1 == ans2 ) {
    } else {
      printf( "Distributivity fail:\n" );
      printf( "  x = 0x%.8x\n", x );
      printf( "  y = 0x%.8x\n", y );
      printf( "  z = 0x%.8x\n", z );
    }
  
    // Associativity test
    if( gf2_32_mult( x, gf2_32_mult( y, z) ) ==
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
