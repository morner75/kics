
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

#' Calculate market Bond Price from Yield-to-Maturity
#'
#' This function calculates the price of a bond with a given tenor, yield-to-maturity, and coupon frequency.
#'
#' @param tenor The number of years remaining until the bond reaches maturity.
#' @param ytm The yield-to-maturity of the bond.
#' @param freq The frequency of coupon payments per year. If freq=0, this indicates a zero-coupon bond.
#'
#' @return The price of the bond.
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
#' @export
wilsonW <- function(u,v,a,ufr) {

  W <- diag(exp(-ufr*u))%*%wilsonH(u,v,a)%*%diag(exp(-ufr*v))

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

#' Wilson G Function
#'
#' Computes the values of the Wilson G function given u and v.
#'
#' @param u A vector of numeric values.
#' @param v A vector of numeric values.
#' @param a A numeric value for parameter a.
#' @return Returns a matrix of numeric values representing the Wilson G function with dimensions length(u) by length(v).
#' @export
wilsonG <- function(u,v,a){

  G <- matrix(NA,nrow=length(u),ncol=length(v))

  for (i in 1:length(u)) {
    for (k in 1:length(v)) {
      if(v[k]<=u[i]){
        G[i,k] <- a - a*exp(-a*u[i])*cosh(a*v[k])
      }else if(u[i]<=v[k]){
        G[i,k] <- a*exp(-a*v[k])*sinh(a*u[i])
      }
    }}
  return(G)

}


#' Discount pricing function with Wilson Function and UFR
#'
#' This function calculates the discounted prices for a given vector of tenors and yields,
#' taking into account the Ultimate Forward Rate (UFR) and the Wilson function.
#'
#' @param v a vector of tenors for which to calculate the discounted prices.
#' @param tenor a vector of the number of years remaining until the bond reaches maturity.
#' @param ytm a vector of yields to maturity of observed tenors.
#' @param a a parameter for the Wilson function.
#' @param ufr the Ultimate Forward Rate.
#' @param freq The frequency of coupon payments per year. If freq=0, this indicates a zero-coupon bond.
#' @return a vector of discounted prices for the given tenors.
#' @export
discount_pricing <- function(v,tenor,ytm,a,ufr,freq=2){

  # market value recalculated from ytm
  p <- ytm2price(tenor,ytm, freq)

  # cash flows from observed tenors and ytms
  Cmat <- cash_flow_mat(tenor,ytm,freq)

  # durations which produce a cash flow
  u <- as.numeric(rownames(Cmat))

  # discount factor with ufr
  d <- exp(-ufr*u)

  # ufr-discounted cash flows
  Q <- diag(as.vector(d))%*%Cmat

  # solution for auxiliary vector b from equation
  b <- solve(t(Q)%*%wilsonH(u,u,a)%*%Q)%*%(p-colSums(Q))

  as.vector(exp(-ufr*v)+exp(-ufr*v)*wilsonH(v,u,a)%*%Q%*%b)
}

#' Discrete spot rates from Smith-Wilson methods
#'
#' This function calculates the spot rates for a given vector of tenors and yields,
#' taking into account the Ultimate Forward Rate (UFR) and the Wilson function.
#'
#' @param v a vector of tenors for which to calculate the discounted prices.
#' @param tenor a vector of the number of years remaining until the bond reaches maturity.
#' @param ytm a vector of yields to maturity of observed tenors.
#' @param a a parameter for the Wilson function.
#' @param ufr the Ultimate Forward Rate.
#' @param freq The frequency of coupon payments per year. If freq=0, this indicates a zero-coupon bond.
#' @return a vector of discounted prices for the given tenors.
#' @export
compute_spots_sw <- function(v,tenor,ytm,a,ufr,freq){

  intensity <- vector(mode="numeric",length=length(v))

  for( i in seq_along(v)){
    if(v[i]==0){

      # market value recalculated from ytm
      p <- ytm2price(tenor,ytm, freq)

      # cash flows from observed tenors and ytms
      Cmat <- cash_flow_mat(tenor,ytm,freq)

      # durations which produce a cash flow
      u <- as.numeric(rownames(Cmat))

      # discount factor with ufr
      d <- exp(-ufr*u)

      # ufr-discounted cash flows
      Q <- diag(as.vector(d))%*%Cmat

      # solution for auxiliary vector b from equation
      b <- solve(t(Q)%*%wilsonH(u,u,a)%*%Q)%*%(p-colSums(Q))

      intensity[i] <- ufr-a*colSums(Q%*%b)+a*exp(-a*u)%*%Q%*%b


    } else{
       Pv <- discount_pricing(v[i],tenor,ytm,a,ufr,freq)
       intensity[i] <- -log(Pv)/v[i]

    }
  }

  res <- cont2disc(intensity,convert="discrete")
  names(res) <- format(v,digits=2)
  res
}


#' Calculate Forward Rates from Spot Rates
#'
#' Given a vector of spot rates, this function calculates the corresponding forward rates.
#'
#' @param spot_rates A numeric vector of spot rates.
#' @return A numeric vector of forward rates.
#' @export
compute_forwards <- function(spot_rates){
  pv <- (1+spot_rates)^(seq_along(spot_rates)-1)
  forward_rates <- `[`(pv/lag(pv)-1,-1)
  if(!is.null(names(spot_rates))) names(forward_rates) <- names(spot_rates[-length(spot_rates)])
  forward_rates
}

