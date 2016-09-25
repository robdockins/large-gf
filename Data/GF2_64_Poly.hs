{-# LANGUAGE ParallelListComp #-}
{-# LANGUAGE MultiWayIf #-}

module Data.GF2_64_Poly
( Poly
, deg
, constant
, term
, addPoly
, negatePoly
, multPoly
, poly
, divModPoly
, egcd
, isZeroPoly
, formal_derivitave
, squarePoly
, unsquarePoly
, scalarMult
, berlekamp
, traceMod
, factorBerlekamp
, roots
, reconstruct
, isoMatrix
, gaussInv
) where

import qualified Data.Map.Strict as IMap
import Data.Foldable
import Data.List
import Data.Maybe
import Data.Word
import Data.Bits

import Test.QuickCheck

import Data.GF2_64

newtype Poly = Poly (IMap.Map Integer GF2_64)
 deriving (Eq, Ord)

instance Show Poly where
  show p@(Poly xs) | deg p == 0 = show (IMap.findWithDefault 0 0 xs)
  show (Poly xs) = concat (intersperse " + " (reverse terms))
    where
    terms = [ show c ++ "·x^" ++ show d
            | (d,c) <- IMap.assocs xs
--            , c /= 0
            ]

instance Arbitrary Poly where
  arbitrary = do
    d <- choose (0,100)
    xs <- vectorOf (d+1) arbitrary
    return $ mkPoly (IMap.fromList (zip [0..] xs))
  shrink (Poly xs) = fmap Poly (shrink xs)

poly :: [(GF2_64, Integer)] -> Poly
poly ts = foldr addPoly (Poly IMap.empty)
            [ term coeff power | (coeff,power) <- ts ]

deg :: Poly -> Integer
deg (Poly xs) | IMap.null xs = 0
              | otherwise    = fst $ IMap.findMax xs

constant :: GF2_64 -> Poly
constant b = term b 0

term :: GF2_64 -> Integer -> Poly
term     0 power = Poly IMap.empty
term coeff power = Poly (IMap.singleton power coeff)

mkPoly :: IMap.Map Integer GF2_64 -> Poly
mkPoly xs = Poly (IMap.filter (/= 0) xs)

addPoly :: Poly -> Poly -> Poly
addPoly (Poly xs) (Poly ys) = Poly zs
 where
  zs = IMap.mergeWithKey f id id xs ys
  f _k x y = let z = x + y in
             if z /= 0 then Just z else Nothing

-- Negation in a field of chararistic 2 is the identity
negatePoly :: Poly -> Poly
negatePoly p = p

multPoly :: Poly -> Poly -> Poly
multPoly (Poly xs) (Poly ys) = zs
 where
  zs = foldr addPoly (Poly IMap.empty) terms
  terms = [ Poly (IMap.mapKeysMonotonic (+i) (fmap (*y) xs))
          | (i, y) <- IMap.assocs ys
          , y /= 0
          ]

-- precondition: input poly is not the constant 0 polynomial
asMonic :: IMap.Map Integer GF2_64 -> (Integer, GF2_64, GF2_64, IMap.Map Integer GF2_64)
asMonic xs = (d,a,a',fmap (*a') xs')
  where ((d,a),xs') = IMap.deleteFindMax xs
        a' = recip a

-- precondition: y is not the constant 0 polynomial
divModPoly :: Poly -> Poly -> (Poly, Poly)
divModPoly p1@(Poly xs) (Poly ys)
  | d1 < d2   = (constant 0, p1)
  | otherwise = (mkPoly q, mkPoly r')
  
 where
 d1 = deg p1
 (d2, a, a', ys') = asMonic ys
 xs' = fmap (*a') xs
 (q, r) = loop IMap.empty xs'
 r' = fmap (*a) r

 loop q x | IMap.null x = (q, x)
          | i < d2      = (q, x)
          | c == 0      = loop q x'
          | otherwise   = loop (IMap.insert j c q) x''
  where
    ((i,c),x') = IMap.deleteFindMax x
    j    = i - d2
    ys'' = IMap.mapKeysMonotonic (+j) (fmap (*c) ys')
    x''  = IMap.mergeWithKey (\_k a b -> Just $ a - b) id (fmap negate) x' ys''

egcd :: Poly -> Poly -> (Poly, Poly, Poly, Poly, Poly)
egcd x y
  | deg x >= deg y = egcd_loop x y (constant 1) (constant 0)
                                   (constant 0) (constant 1)
  | otherwise      = (gcd, v, u, b, a)
     where (gcd, u, v, a, b) = egcd_loop y x (constant 1) (constant 0)
                                             (constant 0) (constant 1)

isZeroPoly :: Poly -> Bool
isZeroPoly (Poly xs) = all (==0) (IMap.elems xs)


scalarMult :: GF2_64 -> Poly -> Poly
scalarMult c (Poly xs) = Poly (fmap (*c) xs)

normalize_egcd :: (Poly, Poly, Poly, Poly, Poly)
               -> (Poly, Poly, Poly, Poly, Poly)
normalize_egcd x@(g, _, _, _, _) | isZeroPoly g = x
normalize_egcd (g, u, v, a, b) = (g', u', v', a', b')
  where
  Poly gs = g
  (i,c,c',gs') = asMonic gs
  g'   = Poly (IMap.insert i 1 gs')
  u'   = scalarMult c' u
  v'   = scalarMult c' v
  a'   = scalarMult c a
  b'   = scalarMult c b

-- precondition: deg x >= deg y
egcd_loop :: Poly -> Poly
          -> Poly -> Poly
          -> Poly -> Poly
          -> (Poly, Poly, Poly, Poly, Poly)
egcd_loop r r' s s' t t' =
   if isZeroPoly r' then
      normalize_egcd (r, s, t, t', s')
   else
      egcd_loop r' r'' s' s'' t' t''
  where
  (q,r'')  = divModPoly r r'
  s''      = addPoly s (negatePoly (multPoly q s'))
  t''      = addPoly t (negatePoly (multPoly q t'))


repeatedAddition :: Integer -> GF2_64 -> Maybe GF2_64
repeatedAddition i x
  | odd i     = Just x
  | otherwise = Nothing

formal_derivitave :: Poly -> Poly
formal_derivitave (Poly xs) = Poly ys
  where
  ys = IMap.mapKeysMonotonic (subtract 1) $
       IMap.mapMaybeWithKey repeatedAddition xs

squarePoly :: Poly -> Poly
squarePoly (Poly xs) = Poly ys
 where
 ys = IMap.mapKeysMonotonic (2*) $
      fmap square xs

unsquarePoly :: Poly -> Poly
unsquarePoly (Poly xs) = Poly ys
 where
 ys = IMap.mapKeysMonotonic (`div` 2) $
      IMap.mapMaybeWithKey f xs
 f i x
   | even i    = Just (frobenius x (-1))
   | otherwise = Nothing

traceMod :: Poly -> Poly -> Poly
traceMod x f = foldr addPoly (constant 0) bs
 where
   bs = berlekamp x f

berlekamp :: Poly -> Poly -> [Poly]
berlekamp x0 f = take 64 xs
  where 
   xs = x0 :
        [ r
        | x <- xs
        , let (_,r) = divModPoly (squarePoly x) f
        ]

factorBerlekamp :: Poly -> [Poly]
factorBerlekamp f
  | f == constant 1 = []
  | deg f < 1       = [f]
  | null xs         = [f]
  | otherwise       = 
      (if deg g1 < deg f then factorBerlekamp g1 else [g1]) ++
      (if deg g2 < deg f then factorBerlekamp g2 else [g2])
  where
  (g1,_,_,g2,_) = egcd f (xs!!0)

  xs = filter (\x -> deg x > 0)
        [ traceMod (term (pow 0x10001 i) 1) f
        | i <- [0..63]
        ]

{-

0x33d22d5be78d6cc is a root of x^64 + x^4 + x^3 + x^1 + x^0

-}

isoMatrix :: [GF2_64]
isoMatrix = take 64 xs
 where
   a  = 0x33d22d5be78d6cc :: GF2_64
   xs = 0x1 : [ a * x | x <- xs ]

gaussInv :: [Word64] -> ([Word64],[Word64])
gaussInv start = reduce 0 start ident
 where
 ident = [ 1 `shiftL` i
         | i <- [0..63]
         ]

 reduce i m v 
   | i >= 64 = elim 0 m v
   | otherwise = reduce (i+1) m'' v''
       where
       idx = findRow i m
       m'  = swapRow i idx m
       v'  = swapRow i idx v
       m'' = reduceRow i m'
       v'' = reduceRow i v'

 swapRow i j m =
   [ if | k == i -> m!!j
        | k == j -> m!!i
        | otherwise -> m!!k
   | k <- [0..63]
   ]

 reduceRow i m = take i m ++ r : xs
   where
   r  = m!!i
   xs = [ x `xor` r
        | x <- drop (i+1) m
        ]

 findRow i m = case xs of
                 (j:_) -> j
                 [] -> error $ unwords ["findRow:", show i, show (map GF2_64 m)]

  where xs = [ j
             | (j,r) <- zip [i..63] (drop i m)
             , testBit r j
             ]

 elim i m v = (m,v)


reconstruct :: Poly
reconstruct =
  foldr multPoly (constant 1)
     [ poly [(1,1),(r,0)]
     | r <- roots
     ]

roots :: [GF2_64]
roots = [
  0x2ac17bd7d83af675,
  0x5814dfabd3629f1f,
  0x8b1e521e802f86a3,
  0x294104f0842da8ee,
  0x858aab7c79da600a,
  0xaaec8347554e149a,
  0xab7babde214b16e1,
  0x8b1e05fe096e4c45,
  0xf7f37fd1825f11b1,
  0xf7f30216066d1d57,
  0x858aaa45fdeca9c6,
  0x739ed608ae358fc8,
  0x7d07092d0c05c83c,
  0xdf3bfc55d350c437,
  0x5814d2a0f0758fe5,
  0x739eab4ab52efd0c,
  0x8c177f2c68a9bed9,
  0xaaecde6034e80ee8,
  0xa83d64b90871d30,
  0x2ac171866b907e36,
  0xab7b7748cffe4f2e,
  0x33d22d5be78d6cc,
  0xdf3b5ae4c7f30608,
  0xa832b866fdd9f81,
  0x7d07dbabc28d0bb3,
  0xf7f37d68925a7d92,
  0x8c17f133423d8440,
  0x33d7cc43f69205d,
  0x29418b82e2362395,
  0xf7f300af16688e99,
  0x2ac113578109217c,
  0x739e6c8d553e1a70,
  0xdf3bcd2b88225b3c,
  0x8c17b4dd7992de2c,
  0x8f359951536f7d78,
  0xab7b1b360ae8bde3,
  0x8b1e1b3903b4dc1a,
  0x785a186f75b77b5d,
  0x29414d48f9a9a31f,
  0x8c173ac25306ec90,
  0x785ae4605235b48f,
  0x7d076e2cf09c320a,
  0x8f351dd05a0c3f1c,
  0xaaecbfcd03b1f660,
  0x8b1e4cd98af510cb,
  0x8f3562b5319c14b5,
  0xab7bc7a0e45da2d1,
  0x785a3836eaca3275,
  0xa8347e66bf5fea6,
  0x858a3c876fd1295c,
  0xdf3b6b9a9c81377e,
  0x858a3dbeebe7fb8f,
  0x739e11cf4e258f88,
  0xaaece2ea6217a9c1,
  0x2ac1190632a3fcbd,
  0x5814bf17bb5164c9,
  0x785ac439cd480d13,
  0xa83ba2b94af6911,
  0x7d07bcaa3e140b8b,
  0x2941c23a9fb27ae6,
  0x33dba7712f3f70b,
  0x8f35e63438ff58b0,
  0x33de46693e207e0,
  0x5814b21c9846c41b]



{-
squareFreeFactors :: Poly -> [(Poly,Int)]
squareFreeFactors p = sff (xs Seq.>< 1) 1 (constant 1)
  where (a,a',xs) = asMonic p

sff :: Seq.Seq GF2_64 -> Int -> Poly -> [(Poly,Int)]
sff f i r = 
  where
   g = formal_derivitave f
   (c,_,_,w0,_) = egcd f g
   loop w
-}