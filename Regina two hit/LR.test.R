#Likelihood Ratio Test

# Deviance = -2log(L_model/L_null)
# -- Deviance ~ Chisq(r), r = n.pars increase from null

LR.test = function(par.vec, nll.vec, nll.null, n.pars){

  ### Calculate LR confidence interval
  if(which.min(nll.vec) > 1){ #If there is a midpoint
  lower = 1:which(nll.vec == min(nll.vec))
  upper = which(nll.vec == min(nll.vec)):length(nll.vec)
  
  #Vector of pars and likelihoods below the min to feed to approx
  pars.lower = par.vec[lower]
  lk.lower = nll.vec[lower]
  
  #Vector of pars and likelihoods above the min to feed to approx
  pars.upper = par.vec[upper]
  lk.upper = nll.vec[upper]
  
  #Use approx to find the point, (y), and likelihood, (x), where the profile
  # crosses the point equal to -log(L0)+Chisq(.95, df)
  lower.CI = approx(lk.lower, pars.lower, xout = (nll.null+qchisq(.95, df = n.pars)/2) ); lower.CI
  upper.CI = approx(lk.upper, pars.upper, xout = (nll.null+qchisq(.95, df = n.pars)/2) ); upper.CI
  
  out = list(lower = lower.CI, upper = upper.CI)

  }else{
    upper.CI = approx(nll.vec, par.vec, xout = (nll.null+qchisq(.95, df = n.pars)/2) )
    out = upper.CI

}
out
}


LR.Threshold = function(NLL_best, df){
  threshold = qchisq(.95, df)/2+NLL_best
  threshold
}

LR.CI = function(threshold, nll.vec, pars.vec){
  if(any(which(nll.vec > threshold) < which.min(nll.vec))  &  any(which(nll.vec > threshold) > which.min(nll.vec)) ){ #If the minimum is not an end point
    #Find string before and after the min value
  lower = nll.vec[1:which.min(nll.vec)]
  upper = nll.vec[(which.min(nll.vec)-1):length(nll.vec)]
  #Extract low CI from first string, and upper CI from upper string
  CI = c(pars.vec[which.min(abs(lower-threshold))], pars.vec[length(lower)+which.min(abs(upper-threshold))])
  }else{
    #If the first value is the minimum
    if(any(which(nll.vec > threshold) < which.min(nll.vec)) == FALSE){ CI = c(pars.vec[1], pars.vec[which.min(abs(nll.vec-threshold))]) }
    if(any(which(nll.vec > threshold) > which.min(nll.vec)) == FALSE){ CI = c(pars.vec[which.min(abs(nll.vec-threshold))], pars.vec[length(nll.vec)])}
  }
  CI
}
  