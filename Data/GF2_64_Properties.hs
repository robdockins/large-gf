{-# LANGUAGE ParallelListComp #-}
module Data.GF2_64_Properties where

import Data.Bits
import Data.Word
import Test.QuickCheck

import Data.GF2_64
import Data.GF2_64_Poly

add_comm :: GF2_64 -> GF2_64 -> Bool
add_comm x y = x+y == y+x

add_assoc :: GF2_64 -> GF2_64 -> GF2_64 -> Bool
add_assoc x y z = x+(y+z) == (x+y)+z

add_unit :: GF2_64 -> Bool
add_unit x = x+0 == x

add_inv :: GF2_64 -> Bool
add_inv x = x + negate x == 0

mult_comm :: GF2_64 -> GF2_64 -> Bool
mult_comm x y = x*y == y*x

mult_assoc :: GF2_64 -> GF2_64 -> GF2_64 -> Bool
mult_assoc x y z = x*(y*z) == (x*y)*z

mult_unit :: GF2_64 -> Bool
mult_unit x = x*1 == x

mult_inv :: GF2_64 -> Property
mult_inv x = x /= 0 ==> x * recip x == 1

mult_distrib :: GF2_64 -> GF2_64 -> GF2_64 -> Bool
mult_distrib x y z = x*(y+z) == (x*y) + (x*z)

square_correct :: GF2_64 -> Bool
square_correct x = square x == x*x


pow_correct :: GF2_64 -> Int -> Bool
pow_correct x i = pow x i == (if i >= 0 then y else recip y)
  where
   i'      = abs i
   y       = foldr (*) 1 factors
   factors = map snd $ filter (testBit i' . fst) squares
   squares = [ (d, y)
             | d <- [0..63]
             | y <- iterate square x
             ]

frobenius_correct :: GF2_64 -> Word8 -> Bool
frobenius_correct x i = frobenius x i == squares !! (fromIntegral i)
  where
    squares = iterate square x


multAddPoly_distrib :: Poly -> Poly -> Poly -> Bool
multAddPoly_distrib x y z =
  multPoly x (addPoly y z) ==
  addPoly (multPoly x y) (multPoly x z)

divModPoly_correct :: Poly -> Poly -> Bool
divModPoly_correct x y = addPoly (multPoly q y) r == x
  where (q,r) = divModPoly x y

egcd_correct :: Poly -> Poly -> Bool
egcd_correct x y =
   multPoly u x `addPoly` multPoly v y == g &&
   multPoly a g == x &&
   multPoly b g == y &&
   deg g <= deg a &&
   deg g <= deg b &&
   (deg b == 0 || (deg u < deg b - deg g)) &&
   (deg a == 0 || (deg v < deg a - deg g))

 where
 (g,u,v,a,b) = egcd x y
