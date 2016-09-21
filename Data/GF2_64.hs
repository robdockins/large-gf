module Data.GF2_64 (
  GF2_64(..)
, pow
, square
, frobenius
) where

import Data.Bits
import Data.Ratio
import Data.Word
import Numeric
import Test.QuickCheck.Arbitrary

foreign import ccall unsafe
  "gf2_64_mult_noinline" gf2_64_mult
  :: Word64 -> Word64 -> Word64

foreign import ccall unsafe
  "gf2_64_square_noinline" gf2_64_square
  :: Word64 -> Word64

foreign import ccall unsafe
  "gf2_64_square16" gf2_64_square16
  :: Word64 -> Word64

foreign import ccall unsafe
  "gf2_64_pow" gf2_64_pow
  :: Word64 -> Word64 -> Word64

foreign import ccall unsafe
  "gf2_64_inv" gf2_64_inv
  :: Word64 -> Word64

newtype GF2_64 = GF2_64 { unGF2_64 :: Word64 }
  deriving (Eq, Ord)

instance Show GF2_64 where
  showsPrec _d (GF2_64 x) = ("0x"++) . showHex x

instance Bounded GF2_64 where
  minBound = GF2_64 0x0
  maxBound = GF2_64 (complement 0x0)

instance Num GF2_64 where
  fromInteger = GF2_64 . fromInteger
  (GF2_64 x) + (GF2_64 y) = GF2_64 (x `xor` y)
  (GF2_64 x) * (GF2_64 y) = GF2_64 (gf2_64_mult x y)
  signum (GF2_64 x) = if x == 0 then GF2_64 0x0 else GF2_64 0x1
  negate x = x
  abs x = x

instance Fractional GF2_64 where
  recip (GF2_64 x) = GF2_64 (gf2_64_inv x)
  fromRational r =
    fromIntegral (numerator r) /
    fromIntegral (denominator r)

instance Arbitrary GF2_64 where
  arbitrary = fmap GF2_64 arbitraryBoundedRandom
  shrink (GF2_64 x) = fmap GF2_64 (shrink x)

instance CoArbitrary GF2_64 where
  coarbitrary (GF2_64 x) = coarbitrary x

pow64 :: GF2_64 -> Word64 -> GF2_64
pow64 (GF2_64 x) i = GF2_64 (gf2_64_pow x i)

pow :: Integral i => GF2_64 -> i -> GF2_64
pow x i
  | i < 0 = recip (pow64 x (fromIntegral (negate i)))
  | otherwise = pow64 x (fromIntegral i)

square :: GF2_64 -> GF2_64
square (GF2_64 x) = GF2_64 (gf2_64_square x)

frobenius8 :: GF2_64 -> Word8 -> GF2_64
frobenius8 (GF2_64 x) i = GF2_64 yfinal
  where sixteens = (i .&. 0x30) `shiftR` 4
        ones     = i .&. 0x0F
        xs       = iterate gf2_64_square16 x
        xfinal   = xs !! (fromIntegral sixteens)
        ys       = iterate gf2_64_square xfinal
        yfinal   = ys !! (fromIntegral ones)

frobenius :: Integral i => GF2_64 -> i -> GF2_64
frobenius x i = frobenius8 x i'
 where
  i' | i >= 0    = fromIntegral i
     | otherwise = 64 - (0x3F .&. (fromIntegral (negate i)))
