#ifndef GF2_64_H
#define GF2_64_H

inline uint32_t gf2_32_add( uint64_t x, uint64_t y) {
  return (x ^ y);
}

uint64_t gf2_64_mult( uint64_t, uint64_t );
uint64_t gf2_64_square( uint64_t );
uint64_t gf2_64_inv( uint64_t );
uint64_t gf2_64_square16( uint64_t );
uint64_t gf2_64_pow( uint64_t, uint64_t );
uint64_t gf2_64_pow_alternate( uint64_t, uint64_t );
int gf2_64_generator( uint64_t a );

uint64_t gf2_64_iso( uint64_t );

#endif
