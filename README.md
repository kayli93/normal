# Normal
Elm library providing simple tools for working with normal distributions. The four main functions are dnormal, pnormal, qnormal, and rnormal, which I've aliased to dnorm, pnorm, qnorm, and rnorm respectively. This is intended to be evocative of the function naming in R. 

# dnormal
dnormal takes a mean, a standard deviation, and a value, and it returns the height of the distribution that value. 

If you try to use a negative standard deviation, dnormal will give you an error. If you try to use a standard deviation of 0, then it will return 0 except if the value equals the mean in which case it returns an error of "Infinity."

dStdNormal does the same thing as dnormal but since it uses the standard normal distribution we don't have to worry about error catching.

Example usage: 
  dnormal 2 1 1 == Ok 0.24197072451914337

# pnormal
pnormal approximates the cumulative distribution function for a normal distribution with a given mean and standard deviation. It will return an error if given a negative standard deviation. If the standard deviation is zero, it will return the indicator function for a given value being greater than or equal to the mean. 

pStdNormal is pnormal for mean 0 and standard deviation 1. We don't need to worry about error catching in this case.

Example usage:
  pnormal 0 2 6 == Ok 0.9986501019683699

The algorithm for this function is pulled from a paper by Graeme West (2004) referenced below. 

# qnormal
qnormal approximates the quantile function for a normal distribution with a given mean and standard deviation. It will return an error of "Negative infinity" if the quantile is less than or equal to 0 and "Infinity" if it is greater than or equal to 1. 

qStdNormal is qnormal for mean 0 and standard deviation 1. We have to catch the same errors here.

Example usage:
  qnormal 0 1 0.975 == Ok 1.96289748249T67147
  
The algorithm for this function is pulled from a paper by Michael Wichura (1977) referenced below. It is, in theory, accurate up to 1 part in 10^(-16). 
  
# rnormal
rnormal creates a generator using the Marsaglia polar method. This creates a generator to be used with the Random module to collect normally distribued samples. 

If the standard deviation is less than 0, it returns an error. If it equals 0, then morally we should have the Dirac delta but that's not easily implementable as far as I know.

# References
  West, Graeme. "Better approximations to cumulative normal functions." Wilmott Magazine 9 (2004): 70-76.
  
Wichura, Michael J. "Algorithm AS 241: The percentage points of the normal distribution." Applied Statistics (1988): 477-484.
