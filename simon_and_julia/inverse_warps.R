library(tidyverse)
library(fdasrvf)


inverse_warps = function(tstar, gam){
  gam_inverse = approx(gam*(max(tstar) - min(tstar)) + min(tstar), tstar, xout = tstar)$y
  
  gam_inverse
}






