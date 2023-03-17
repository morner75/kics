
cash_flow <- function(tenor,ytm,freq=2){
  dt <- 1/freq
  if (tenor < dt){
    t <- tenor
    CF <- 1
  }else{
    t <- setdiff(seq(tenor%%dt, tenor, by = dt),0)
    CF <- c(rep(ytm/freq,length(t)-1), 1 + ytm/freq)
  }
  names(CF) <- t
  CF
}



#' cash_flow_matrix
#'
#' Construct cash flow matrix
#'
#' @param tenor Vector of bond maturities
#' @param ytm Vector of yield to maturities
#' @param freq Compounding frequency per year (default is 2)
#'
#' @return A matrix with cash flows for different tenors and yield to maturities
#' @export
cash_flow_mat <- function(tenor,ytm,freq=2){

  dt <- 1/freq
  cum_t <- NULL
  for(tr in tenor){
    temp_t <- setdiff(seq(tr%%dt, tr, by = dt),0)
    cum_t <- union(cum_t, temp_t)
  }

  C_mat <- matrix(0, ncol=length(tenor),nrow=length(cum_t))
  rownames(C_mat) <- as.character(cum_t)
  colnames(C_mat) <- as.character(tenor)


  for(i in seq_along(tenor)){
    cf <- cash_flow(tenor=tenor[i],ytm=ytm[i])
    C_mat[match(names(cf),rownames(C_mat)),match(tenor[i],colnames(C_mat))] <- cf
  }
  C_mat
}

ytm2price0 <- function (tenor, ytm, freq) {

  if(ytm > 1) ytm <- ytm/100

  dt <- 1/freq
  if (freq == 0)  return(1/(1 + ytm)^tenor)
  if (tenor < dt)  return(1/(1 +ytm*tenor))

  t <- setdiff(seq(tenor%%dt, tenor, by = dt),0)


  CF <- c(rep(ytm/freq,length(t)-1), 1 + ytm/freq)
  DF = 1/(t%%dt*ytm+1)*(1 + ytm/freq)^(-(t %/%dt))

  sum(CF * DF)
}

#' Yield to maturity to market price
#'
#' Computes the price of a bond given its tenor, yield to maturity, and coupon frequency.
#'
#' @param tenor length of time until maturity of the bond, in years.
#' @param ytm the yield to maturity of the bond.
#' @param freq the number of coupon payments per year.
#'
#' @return The price of the bond.
#'
#'
#' @export
ytm2price<- function (tenor, ytm, freq) {

  if(length(tenor)!=length(ytm)) stop("tenor and ytm have the same length")

  n <- length(tenor)

  if(n==1){
    return(ytm2price0(tenor=tenor,ytm=ytm,freq=freq))
  }else{
    price <- vector(mode="double",length=n)
    for(i in seq_len(n)){
      price[i] <- ytm2price0(tenor=tenor[i],ytm=ytm[i],freq=freq)
    }
    return(price)
  }

}

#' wilston W
#'
#' Calculates the Wilson's W function with given parameters
#'
#' @param u A vector of time values, observed durations to maturity
#' @param v A vector of maturity values
#' @param a A numeric value, speed of convergence
#' @param ufr A numeric value, ultimate forward intensity(continuous rate)
#'
#' @return Returns a matrix of Wilston's function values
#'
#' @examples
#' WilsonW(u = c(1, 2, 3), v = c(2, 3, 4), a = 0.05, ufr = 0.03)
#'
#' @export
wilsonW <- function(u,v,a,ufr) {

  W <- diag(exp(-ufr*u))%*%WilsonH(u,v,a)%*%diag(exp(-ufr*v))

  return(W)
}



#' wilston H
#'
#' Calculates the Wilson's H function with given parameters
#'
#' @param u A vector of time values, observed durations to maturity
#' @param v A vector of maturity values
#' @param a A numeric value, speed of convergence
#'
#' @return Returns a matrix of Wilston's function values
#'
#' @export
wilsonH <- function(u,v,a) {

  H <- matrix(NA,length(u),length(v))

  for (i in 1:length(u)) {
    for (k in 1:length(v)) {
      H[i,k] <- a*min(u[i],v[k])+(exp(-a*(u[i]+v[k]))-exp(-a*abs(u[i]-v[k])))/2
    }}
  return(H)
}



