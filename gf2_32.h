#ifndef GF2_32_H
#define GF2_32_H

inline uint32_t gf2_32_add( uint32_t x, uint32_t y) {
  return (x ^ y);
}
  
uint32_t gf2_32_mult( uint32_t, uint32_t );

#endif
