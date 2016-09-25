typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long uint64_t;

#define FIELD_WIDTH 16
#define FIELD_SIZE (1 << FIELD_WIDTH)
#define Q (FIELD_SIZE - 1)
#define gf2_16_irreducible ((1 << 12) | (1<<3) | (1<<1) | (1<<0))

//__constant unsigned long MEM_SIZE = 128 * 10000;
__constant unsigned long TASK_SIZE = 10000;


// irreducible polynomial is
//   x^17 + x^12 + x^3 + x + 1

// This function computes, in a branch-free way, the value:
//    (x == 0) ? 0xFFFF : 0x0000
inline uint16_t zeroMask( uint16_t x ) {
  uint32_t y = (uint32_t) x;
  uint32_t z = (y - 1) >> 16;
  return (uint16_t) z;
}

inline uint16_t logsum_modQ( uint16_t a, uint16_t b ) {
  uint32_t sum = ((uint32_t) a) + ((uint32_t) b);
  uint16_t d = ((uint16_t) sum) + ((uint16_t) (sum >> 16));
  return d;
}

// Compute exp( (a + b)%Q )
inline uint16_t gf2_16_expadd( __global uint16_t* exp_table, uint16_t zmask, uint16_t a, uint16_t b ) {
  uint32_t sum = ((uint32_t) a) + ((uint32_t) b);
  uint16_t t = exp_table[ sum ];
  return (zmask & t) ^ t;
}

// Compute exp( (a + b + c)%Q )
inline uint16_t gf2_16_expadd3( __global uint16_t* exp_table, uint16_t zmask, uint16_t a, uint16_t b, uint16_t c ) {
  uint32_t sum = ((uint32_t) logsum_modQ( a, b )) + ((uint32_t) c);
  uint16_t t = exp_table[ sum ];
  return (zmask & t) ^ t;
}

__kernel
void gf2_16_mult( __global uint16_t* log_table,
                  __global uint16_t* exp_table,
                  __global uint16_t* x,
                  __global uint16_t* y,
                  __global uint16_t* out )
{
  size_t base = TASK_SIZE * get_global_id(0);
  for(long i=0; i<TASK_SIZE; i++) {
    uint16_t xv = x[base + i];
    uint16_t yv = y[base + i];
    uint16_t a = log_table[xv];
    uint16_t b = log_table[yv];
    uint16_t z = zeroMask(xv) | zeroMask(yv);
    out[base + i] = gf2_16_expadd( exp_table, z, a, b );
  }
}

inline
uint32_t do_gf2_32_mult( __global uint16_t* log_table,
                      __global uint16_t* exp_table,
                      uint32_t a, uint32_t b ) {
  uint16_t a0, a1, b0, b1, za0, za1, zb0, zb1;
  uint32_t alog0, alog1, blog0, blog1;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);

  alog0 = log_table[a0];
  alog1 = log_table[a1];
  blog0 = log_table[b0];
  blog1 = log_table[b1];

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );

  uint16_t c0, c1, c2;

  c0  = gf2_16_expadd( exp_table, za0|zb0, alog0, blog0 );

  c1  = gf2_16_expadd( exp_table, za0|zb1, alog0, blog1 );
  c1 ^= gf2_16_expadd( exp_table, za1|zb0, alog1, blog0 );

  c2  = gf2_16_expadd( exp_table, za1|zb1, alog1, blog1 );

  // Now reduce by the field extension irreducible
  // polynomial x^2 + x + 8192

  //uint16_t vlog = gf2_16_log_table[8192];
  // this magic number is equal to gf2_16_log_table[8192];
  const uint16_t vlog = 0x000d;

  c1 ^= c2;
  c0 ^= gf2_16_expadd3( exp_table, zeroMask(c2), alog1, blog1, vlog );

  // polynomial is now reduced, output the result
  uint32_t res = (((uint32_t) c1) << 16) | ((uint32_t) c0);
  return res;
}



__kernel
void gf2_32_mult( __global uint16_t* log_table,
                  __global uint16_t* exp_table,
                  __global uint32_t* x,
                  __global uint32_t* y,
                  __global uint32_t* out )
{
  size_t base = TASK_SIZE * get_global_id(0);
  for(long i=0; i<TASK_SIZE; i++) {
    out[base + i] = do_gf2_32_mult( log_table, exp_table, x[base+i], y[base+i] );
  }
}


inline uint64_t do_gf2_64_mult( __global uint16_t* log_table,
                                __global uint16_t* exp_table,
                                uint64_t a, uint64_t b ) {
  uint16_t a0, a1, a2, a3;
  uint16_t za0, za1, za2, za3;
  uint16_t alog0, alog1, alog2, alog3;

  uint16_t b0, b1, b2, b3;
  uint16_t zb0, zb1, zb2, zb3;
  uint16_t blog0, blog1, blog2, blog3;

  uint16_t c0, c1, c2, c3, c4, c5, c6;
  uint16_t d0, d1, d2, d3;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  a2 = (uint16_t) (a >> 32);
  a3 = (uint16_t) (a >> 48);

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = log_table[a0];
  alog1 = log_table[a1];
  alog2 = log_table[a2];
  alog3 = log_table[a3];

  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);
  b2 = (uint16_t) (b >> 32);
  b3 = (uint16_t) (b >> 48);

  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );
  zb2 = zeroMask( b2 );
  zb3 = zeroMask( b3 );

  blog0 = log_table[b0];
  blog1 = log_table[b1];
  blog2 = log_table[b2];
  blog3 = log_table[b3];

  // Do the raw coefficent multiplications
  // via GF(2^16) tables

  c0  = gf2_16_expadd( exp_table, za0|zb0, alog0, blog0 );

  c1  = gf2_16_expadd( exp_table, za0|zb1, alog0, blog1 );
  c1 ^= gf2_16_expadd( exp_table, za1|zb0, alog1, blog0 );

  c2  = gf2_16_expadd( exp_table, za1|zb1, alog1, blog1 );
  c2 ^= gf2_16_expadd( exp_table, za2|zb0, alog2, blog0 );
  c2 ^= gf2_16_expadd( exp_table, za0|zb2, alog0, blog2 );

  c3  = gf2_16_expadd( exp_table, za1|zb2, alog1, blog2 );
  c3 ^= gf2_16_expadd( exp_table, za2|zb1, alog2, blog1 );
  c3 ^= gf2_16_expadd( exp_table, za3|zb0, alog3, blog0 );
  c3 ^= gf2_16_expadd( exp_table, za0|zb3, alog0, blog3 );

  c4  = gf2_16_expadd( exp_table, za1|zb3, alog1, blog3 );
  c4 ^= gf2_16_expadd( exp_table, za2|zb2, alog2, blog2 );
  c4 ^= gf2_16_expadd( exp_table, za3|zb1, alog3, blog1 );

  c5  = gf2_16_expadd( exp_table, za2|zb3, alog2, blog3 );
  c5 ^= gf2_16_expadd( exp_table, za3|zb2, alog3, blog2 );

  c6  = gf2_16_expadd( exp_table, za3|zb3, alog3, blog3 );

  uint16_t log2  = log_table[2];
  uint16_t c4log = log_table[c4];
  uint16_t c5log = log_table[c5];
  uint16_t c6log = log_table[c6];

  // Now perform the modular reduction
  uint16_t c4x2 = gf2_16_expadd( exp_table, zeroMask(c4), c4log, log2 );
  uint16_t c5x2 = gf2_16_expadd( exp_table, zeroMask(c5), c5log, log2 );
  uint16_t c6x2 = gf2_16_expadd( exp_table, zeroMask(c6), c6log, log2 );

  d3  = c3 ^ c5 ^ c6x2;
  d2  = c2 ^ c4 ^ c5x2;
  d1  = c1 ^ c4x2 ^ c5 ^ c6x2;
  d0  = c0 ^ c4 ^ c6;

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}



__kernel
void gf2_64_mult( __global uint16_t* log_table,
                  __global uint16_t* exp_table,
                  __global uint64_t* x,
                  __global uint64_t* y,
                  __global uint64_t* out )
{
  size_t base = TASK_SIZE * get_global_id(0);
  for(long i=0; i<TASK_SIZE; i++) {
    out[base + i] = do_gf2_64_mult( log_table, exp_table, x[base+i], y[base+i] );
  }
}

// Multiply a*b, but only calculate the low coefficent.
inline uint16_t gf2_64_mult_low_coeff( __global uint16_t* log_table,
                                       __global uint16_t* exp_table,
                                       uint64_t a, uint64_t b ) {
  uint16_t a0, a1, a2, a3;
  uint16_t za0, za1, za2, za3;
  uint16_t alog0, alog1, alog2, alog3;

  uint16_t b0, b1, b2, b3;
  uint16_t zb0, zb1, zb2, zb3;
  uint16_t blog0, blog1, blog2, blog3;

  uint16_t t, z;
  uint16_t c0, c4, c6;
  uint16_t d0;

  a0 = (uint16_t) a;
  a1 = (uint16_t) (a >> 16);
  a2 = (uint16_t) (a >> 32);
  a3 = (uint16_t) (a >> 48);

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = log_table[a0];
  alog1 = log_table[a1];
  alog2 = log_table[a2];
  alog3 = log_table[a3];

  b0 = (uint16_t) b;
  b1 = (uint16_t) (b >> 16);
  b2 = (uint16_t) (b >> 32);
  b3 = (uint16_t) (b >> 48);

  zb0 = zeroMask( b0 );
  zb1 = zeroMask( b1 );
  zb2 = zeroMask( b2 );
  zb3 = zeroMask( b3 );

  blog0 = log_table[b0];
  blog1 = log_table[b1];
  blog2 = log_table[b2];
  blog3 = log_table[b3];

  // Do the raw coefficent multiplications
  // via GF(2^16) tables

  c0  = gf2_16_expadd( exp_table, za0|zb0, alog0, blog0 );

  c4  = gf2_16_expadd( exp_table, za1|zb3, alog1, blog3 );
  c4 ^= gf2_16_expadd( exp_table, za2|zb2, alog2, blog2 );
  c4 ^= gf2_16_expadd( exp_table, za3|zb1, alog3, blog1 );

  c6  = gf2_16_expadd( exp_table, za3|zb3, alog3, blog3 );

  d0  = c0 ^ c4 ^ c6;

  return d0;
}

// This rather magical-looking function computes
// the value a^(2^16).  Because the field has characteristic
// two, squaring commutes with addition, as does any number
// of iterations of squaring.  Moreover, raising the coefficents
// (which lie in GF(2^16)) to the 2^16 power is the identity
// function.  Thus, the overall effect of raising to the 2^16
// power is a linear transformation of the original coefficents.
//
// The calculations below implement multiplication by the following
// matrix, which represents the action of the linear transformation.
//
//  0x0001  0x0000  0x0000  0x0000
//  0x0fd3  0x9b04  0x1f3f  0x0000
//  0xc870  0x393c  0x9b04  0x0000
//  0x29d2  0x6d0b  0x36ef  0x0001
//
inline uint64_t gf2_64_square16( __global uint16_t* log_table,
                                 __global uint16_t* exp_table,
                                 uint64_t a ) {
  uint16_t a0, a1, a2, a3;
  uint16_t d0, d1, d2, d3;
  uint16_t za0, za1, za2, za3;
  uint16_t alog0, alog1, alog2, alog3;

  a3 = (uint16_t) (a >> 48);
  a2 = (uint16_t) (a >> 32);
  a1 = (uint16_t) (a >> 16);
  a0 = (uint16_t) a;

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = log_table[a0];
  alog1 = log_table[a1];
  alog2 = log_table[a2];
  alog3 = log_table[a3];

  /* d3 = a3; */
  /* d2 = gf2_16_expadd( za3, alog3, gf2_16_log_table[0x0fd3] ) ^ */
  /*      gf2_16_expadd( za2, alog2, gf2_16_log_table[0x9b04] ) ^ */
  /*      gf2_16_expadd( za1, alog1, gf2_16_log_table[0x1f3f] ); */
  /* d1 = gf2_16_expadd( za3, alog3, gf2_16_log_table[0xc870] ) ^ */
  /*      gf2_16_expadd( za2, alog2, gf2_16_log_table[0x393c] ) ^ */
  /*      gf2_16_expadd( za1, alog1, gf2_16_log_table[0x9b04] ); */
  /* d0 = gf2_16_expadd( za3, alog3, gf2_16_log_table[0x29d2] ) ^ */
  /*      gf2_16_expadd( za2, alog2, gf2_16_log_table[0x6d0b] ) ^ */
  /*      gf2_16_expadd( za1, alog1, gf2_16_log_table[0x36ef] ) ^ */
  /*      a0; */

  // The below magic values are the GF(2^16) logarithms of the a
  // above magic numbers, which represent the matrix defining the
  // linear transformation, a^(2^16).
  d3 = a3;
  d2 = gf2_16_expadd( exp_table, za3, alog3, 0x6cd0 ) ^
       gf2_16_expadd( exp_table, za2, alog2, 0x3190 ) ^
       gf2_16_expadd( exp_table, za1, alog1, 0x6321 );
  d1 = gf2_16_expadd( exp_table, za3, alog3, 0x2514 ) ^
       gf2_16_expadd( exp_table, za2, alog2, 0xc643 ) ^
       gf2_16_expadd( exp_table, za1, alog1, 0x3190 );
  d0 = gf2_16_expadd( exp_table, za3, alog3, 0x49fd ) ^
       gf2_16_expadd( exp_table, za2, alog2, 0xd8dd ) ^
       gf2_16_expadd( exp_table, za1, alog1, 0x53f1 ) ^
       a0;

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}

inline uint64_t gf2_64_pointwise_mult( __global uint16_t* log_table,
                                       __global uint16_t* exp_table,
                                       uint16_t xlog, uint64_t a ) {
  uint16_t a0, a1, a2, a3;
  uint16_t d0, d1, d2, d3;

  uint16_t za0, za1, za2, za3;
  uint16_t alog0, alog1, alog2, alog3;

  a3 = (uint16_t) (a >> 48);
  a2 = (uint16_t) (a >> 32);
  a1 = (uint16_t) (a >> 16);
  a0 = (uint16_t) a;

  za0 = zeroMask( a0 );
  za1 = zeroMask( a1 );
  za2 = zeroMask( a2 );
  za3 = zeroMask( a3 );

  alog0 = log_table[a0];
  alog1 = log_table[a1];
  alog2 = log_table[a2];
  alog3 = log_table[a3];

  d3 = gf2_16_expadd( exp_table, za3, alog3, xlog );
  d2 = gf2_16_expadd( exp_table, za2, alog2, xlog );
  d1 = gf2_16_expadd( exp_table, za1, alog1, xlog );
  d0 = gf2_16_expadd( exp_table, za0, alog0, xlog );

  uint64_t d =
    (((uint64_t) d3) << 48) |
    (((uint64_t) d2) << 32) |
    (((uint64_t) d1) << 16) |
    ((uint64_t) d0);

  return d;
}

uint64_t do_gf2_64_inv( __global uint16_t* log_table,
                        __global uint16_t* exp_table,
                        uint64_t a ) {

  // Let q = 2^16.  Let r = (q^4 - 1)/(q - 1) = 2^48 + 2^32 + 2^16 + 1
  // This rather special number has to do with caluclating finite field norms.
  // For all x in GF((2^16)^4), x^r yields a value in GF(2^16); that is, for which
  // the high coefficents are 0.  We exploit this fact to perform fast inversions
  // in GF(2^64) by reducing them to inversion in GF(2^16).
  //
  // We are going to calculate a^(-1) = a^(-r) * a^(r-1).  The algorithm below goes
  // roughly as follows:
  //
  //   s = a^(r-1)
  //     = a^(q^3 + q^2 + q)
  //     = a^(q(q^2 + q + 1))
  //     = a^(q(q(q+1)+1))
  //     = (((a^q * a)^q)*a)^q
  //   t = a * s = a^r
  //   b = t^(-1) * s = a^(-1)

  // Compute s = a^(r-1)
  uint64_t s = a;
  for(int i=0;;i++) {
    s = gf2_64_square16( log_table, exp_table, s );
    if(i>=2) break;
    s = do_gf2_64_mult( log_table, exp_table, s, a );
  }

  // t = s * a = a^r
  // Because we know t in GF(2^16), we can save
  // some work by only calculating the low coefficent
  // of this multiplication.
  uint16_t t0 = gf2_64_mult_low_coeff( log_table, exp_table, s, a );

  // Now invert t0
  //uint16_t t0_inv = gf2_16_inv( t0 );
  uint16_t t0_inv_log = Q - log_table[ t0 ];

  // b = t^(-1) * s = a^(-r) * a^(r-1) = a(-1)
  uint64_t b = gf2_64_pointwise_mult( log_table, exp_table, t0_inv_log, s );
  return b;
}

__kernel
void gf2_64_inv( __global uint16_t* log_table,
                 __global uint16_t* exp_table,
                 __global uint64_t* x,
                 __global uint64_t* out )
{
  size_t base = TASK_SIZE * get_global_id(0);
  for(long i=0; i<TASK_SIZE; i++) {
    out[base + i] = do_gf2_64_inv( log_table, exp_table, x[base+i] );
  }
}
