module Normal (dnorm, dnormal, dStdNormal, pnorm, pnormal , pStdNormal,
    qnorm, qnormal, qStdNormal, rnorm, rnormal, rStdNormal) where
{-| This library implements some fairly standard elementary functions dealing
with normal distributions. For each function, I also created a short-hand that
matches the function names in R.
-}

import Random exposing (Generator)
import Result exposing (Result)

{-| dnormal calculates the pdf at the point for a normal with a given mean and
standard deviation. For example, 
    dnormal 2 1 1 == Ok 0.24197072451914337

If the standard deviation is less than 0, then it returns an error. For standard
deviation equal to 0, then if the mean is not equal to the point x it returns 0
and returns and error of infinity otherwise.

dStdNormal unwraps the Result constructor since calculations can always take
place. 
-}
dnormal : Float -> Float -> Float -> Result String Float
dnormal mean std x =
    if | std < 0   -> Err "Standard deviation cannot be less than 0."
       | std == 0  ->
           if x == mean then Err "Infinity."
                        else Ok 0
       | otherwise -> Ok ((1/(std*(sqrt (2*pi))))*(e^(-(((x-mean)^2)/(2*(
            std^2))))))

dnorm = dnormal

dStdNormal : Float -> Float
dStdNormal x =
    let y  = dnormal 0 1 x
    in
       case y of
           Ok y' -> y'

{-| Approximates normal cdf function using the method outlined in West 2004. For
example, if you wanted to know the cdf value at 6 for a normal distribution with
mean 0 and standard deviation 2 you would use 
    pnormal 0 2 6 == Ok 0.9986501019683699

If std is less than 0, it returns an error. If it equals 0, then it returns Ok 0
if x is less than the mean and Ok 1 otherwise.

REFERENCE 
    West, Graeme. "Better approximations to cumulative normal functions."
    Wilmott Magazine 9 (2005): 70-76.
-}
pnormal : Float -> Float -> Float -> Result String Float
pnormal mean std x =
    if | std < 0   -> Err "Standard deviation cannot be less than 0."
       | std == 0  ->
           if x < mean then Ok 0 else Ok 1
       | otherwise -> Ok (pStdNormal ((x - mean)/std))

pnorm = pnormal

pStdNormal : Float -> Float
pStdNormal x =
    let cumNorm z =
        if | ((abs z) > 37)||(abs z == 37) -> 0
           | (abs z < 7.07106781186547)||(abs z == 7.07106781186547) ->
               let exp = e^((-((abs z)^2))/2)
                   b1 = ((((((0.0352624965998911*(abs z)) +
                   0.700383064443688)*(abs z) + 6.37396220353165)*(abs z) +
                   33.912866078383)*(abs z) + 112.079291497871)*(abs z) +
                   221.213596169931)*(abs z) + 220.206867912376
                   b2 = (((((((0.0883883476483184*(abs x)) +
                   1.75566716318264)*(abs z) + 16.064177579207)*(abs z) +
                   86.7807322029461)*(abs z) + 296.564248779674)*(abs z) +
                   637.333633378831)*(abs z) + 793.826512519948)*(abs z) +
                   440.413735824752
               in
                (exp*b1)/b2
           | ((abs z) > 7.07106781186547)&&(abs z < 37) ->
               let exp = e^((-((abs z)^2))/2)
                   b = (abs z) + (1/((abs z) + (2/((abs z) + (3/((abs z) +
                   (4/((abs z) + 0.65))))))))
               in
                (exp/b)/(2.506628274631)
    in
       if x > 0 then 1- (cumNorm x)
                else cumNorm x

{-| Calculates the quantile function for the Normal distribution using Algorithm
AS 241 found in Wichura 1977. For values outside of (0,1), the function returns
Nothing since the quantile function is not defined outside of those values. This
function might be useful if you want to find a cut-off z-score for something
normally distributed with mean m standard deviation s and significant
level p (modulo some arithmetic). For example,
    
    qnormal 0 1 0.975 == Ok 1.9628974824967147

This implementation should theoretically be accurate up to 1 part in 10^(-16).

If the p is less than or equal to 0, it returns an error for
negative infinity. Similarly, if p is greater than or equal to 1, it returns an
error for infinity.

REFERENCE 
    Wichura, Michael J. "Algorithm AS 241: The percentage points of the normal
    distribution." Applied Statistics (1988): 477-484.
-}
qnormal : Float -> Float -> Float -> Result String Float
qnormal mean std p =
    let z = qStdNormal p
    in
       case z of
           Ok y  -> Ok (mean + std*y)
           Err x -> Err x

qnorm = qnormal

qStdNormal : Float -> Result String Float
qStdNormal p =
    let q = p-0.5
    in
       if | (p < 0)||(p == 0) -> Err "Negative infinity."
          | (p == 1)||(p > 1) -> Err "Infinity."
          | (abs q < 0.425)||(abs q == 0.425) ->
                let r = 0.180625 - q*q
                in
                    Ok (q * (((((((r * 2509.0809287301226727
                    +33430.575583588128105) * r + 67265.770927008700853) * r
                    +45921.953931549871457) * r + 13731.693765509461125) * r
                    +1971.5909503065514427) * r + 133.14166789178437745) * r
                    +3.387132872796366608)/ (((((((r * 5226.495278852854561
                    +28729.085735721942674) * r + 39307.89580009271061) * r
                                        +21213.794301586595867) * r + 5394.1960214247511077) * r
                    +687.1870074920579083) * r + 42.313330701600911252) * r +
                    1))
          | otherwise ->
            if (q < 0)
              then
                let r' = sqrt (-(logBase e p))
                in
                   if (r' > 5.0)
                      then
                        let r = r'-5.0
                        in
                           Ok (-((((((((r * 2.01033439929228813265e-7
                           + 2.71155556874348757815e-5) * r
                           + 0.0012426609473880784386) * r +
                           0.026532189526576123093) * r +
                           0.29656057182850489123) * r + 1.7848265399172913358)
                           * r + 5.4637849111641143699) *r +
                           6.6579046435011037772)/(((((((r
                           * 2.04426310338993978564e-15 +
                           1.4215117583164458887e-7) * r +
                           1.8463183175100546818e-5) * r +
                           7.868691311456132591e-4) * r +
                           0.0148753612908506148525)* r +
                           0.13692988092273580531) * r + 0.59983220655588793769)
                           * r + 1)))
                      else
                        let r = r'-1.6
                        in
                           Ok (-((((((((r * 7.7454501427834140764e-4 +
                           0.0227238449892691845833) * r +
                           0.24178072517745061177) * r + 1.27045825245236838258)
                           * r + 3.64784832476320460504) * r +
                           5.7694972214606914055) * r + 4.6303378461565452959) *
                           r + 1.42343711074968357734)/ (((((((r *
                           1.05075007164441684324e-9 +
                           5.475938084995344946e-4) *r +
                           0.0151986665636164571966) * r +
                           0.14810397642748007459) * r + 0.68976733498510000455)
                           * r + 1.6763848301838038494) * r +
                           2.05319162663775882187) * r + 1)))
              else
                 let r' = sqrt (-(logBase e (1-p)))
                 in
                    if (r' > 5.0)
                      then
                        let r = r'-5.0
                        in
                           Ok ((((((((r * 2.01033439929228813265e-7 +
                           2.71155556874348757815e-5) * r +
                           0.0012426609473880784386) * r +
                           0.026532189526576123093) * r +
                           0.29656057182850489123) * r +1.7848265399172913358) *
                           r + 5.4637849111641143699) * r +
                           6.6579046435011037772)/(((((((r *
                           2.04426310338993978564e-15 +
                           1.4215117583164458887e-7)*r +
                           1.8463183175100546818e-5) * r +
                           7.868691311456132591e-4) * r +
                           0.0148753612908506148525)* r +
                           0.13692988092273580531) * r + 0.59983220655588793769)
                           * r + 1))
                      else
                        let r = r'-1.6
                        in
                           Ok ((((((((r * 7.7454501427834140764e-4 +
                           0.0227238449892691845833) * r +
                           0.24178072517745061177) * r + 1.27045825245236838258)
                           * r + 3.64784832476320460504) * r +
                           5.7694972214606914055) * r + 4.6303378461565452959) *
                           r + 1.42343711074968357734)/(((((((r *
                           1.05075007164441684324e-9 +
                           5.475938084995344946e-4) * r +
                           0.0151986665636164571966) * r +
                           0.14810397642748007459) * r + 0.68976733498510000455)
                           * r + 1.6763848301838038494) * r +
                           2.05319162663775882187) * r + 1))

{-| Sample from a normal distribution using the Marsaglia polar method. This
creates a generator to be used with the Random module to collect normally
distribued samples. 

If std is less than 0, it returns an error. If it equals 0, then morally we
should have the Dirac delta but that's not easily implementable as far as I
know.
-}
rnormal : Float -> Float -> Result String (Generator Float)
rnormal mean std =
    if | std < 0 -> Err "Standard deviation cannot be less than 0."
       | std == 0 -> Err "Dirac delta."
       | otherwise ->
            let normGen s =
                let
                    (x, s') = Random.generate (Random.float -1 1) s
                    (y, s'') = Random.generate (Random.float -1 1) s'
                    z = x^2 + y^2
                in
                    if z < 1
                      then (mean + std*(x*(sqrt ((-2*(logBase e z))/z))) ,s'')
                      else normGen s''
            in
                Ok (Random.customGenerator normGen)

rnorm = rnormal

rStdNormal : Generator Float
rStdNormal =
    let r = rnormal 0 1
    in
      case r of
        Ok r' -> r'
