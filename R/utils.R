
#' 지수평활화가중평균
#'
#' Compute Exponetially Weighted Moving Average.
#'
#' @param vector A numeric vector, 평균을 구할 시계열.
#' @param lambda A numeric 평활화 승수값.
#' @param round Logical 5bp 단위 반올림?
#'
#' @return A numeric.
#' @export
#'
#' @examples
#' x <- 1:20
#' expSmooth(x,lambda=0.85)
expSmooth <- function(vector,lambda=NULL, round=TRUE){

  if(is.null(lambda)) lambda=2/(length(vector)+1)

  n <- length(vector)
  weight <- c(1 ,rep(lambda,times=n-1))*(1-lambda)^rev(seq_along(vector)-1)
  avg <- sum(weight*vector)

  if(round) return(round(20*avg)/20)
  return(avg)
}

#' Update Long Term Forward Rate
#'
#' 장기선도금리 업데이트
#'
#' @param data A data.frame which contains yearly CD rates(INT_CD91_Y) and CPI(CPI_YG) series.
#' @param update_year A numeric 업데이트할 연도.
#' @param last_LTFR  A numeric 마지막 장기선도금리.
#' @param BOK_target A numeric 한국은행 물가안정목표.
#' @param var_names A characteㄱ vector of length 3, 연도, 무위험지표금리, 인플레이션율 변수명을 지정하는 벡터
#'
#' @return A numeric 다음년도 장기선도금리.
#' @export
#'
#' @examples
#' print("추후에...")
update_LTFR <- function(data,update_year,last_LTFR,BOK_target=2,
                        var_names=c(year="year",rfr="INT_CD91_Y",inflation="CPI_YG")){

  if(!is.data.frame(data)) stop("data must be a data.frame")
  if(sum(var_names %in% colnames(data))<2) stop("data must have year, rfr, inflation columns")

  data %>%
    select(all_of(var_names)) %>%
    set_names(names(var_names)) %>%
    filter(year>=1991,year<=update_year) %>%
    mutate(INT_REAL=(rfr-inflation)/(1+inflation/100)) %>%
    summarise(LTFR=expSmooth(INT_REAL)+BOK_target) %>%
    mutate(LTFR=case_when(LTFR-last_LTFR > 0.15 ~ last_LTFR + 0.15,
                          LTFR-last_LTFR >= -0.15 ~ last_LTFR,
                          TRUE ~ last_LTFR-0.15)) %>%  pull()

}

#' Effective interest rate and nominal interest rate converter
#'
#' This function converts either the nominal or effective interest rate to the other one.
#'
#' @param int A decimal number indicating the interest rate. This must be less than or equal to 1.
#' @param conversion An integer indicating the number of compounding periods in a year. Default is 2.
#' @param convert A string indicating which type of interest rates to convert to. Either "effective" or "nominal". Default is "effective".
#'
#' @return A decimal number representing either the effective interest rate or nominal interest rate.
#' @examples
#' nom2eff(0.12,4,"effective")
#' nom2eff(0.06,2,"nominal")
#' @export
nom2eff <- function(int,conversion=2,convert=c("effective","nominal")) {

  if(int>1) stop("interest rates must be inputed in a decimal not in a percentage")
  freq <- conversion
  convert <- match.arg(arg=convert, choices=c("effective","nominal"))
  switch(convert,
         effective= (1+int/freq)^freq-1,
         nominal=freq*((1+int)^(1/freq)-1)
  )
}


#' Convert Interest Rates from Continuous to Discrete (or vice versa)
#'
#' This function allows users to convert interest rates between continuous and discrete forms.
#'
#' @param int A numeric value representing the interest rate to be converted.
#' @param convert A character string specifying the direction of the conversion. Acceptable values are "discrete" and "continuous".
#' @return A numeric value representing the converted interest rate.
#' @examples
#' cont2disc(0.05,"discrete")
#' cont2disc(0.126,"continuous")
#' @export
cont2disc <- function(int,convert=c("discrete","continuous")) {

  if(int>1) stop("interest rates must be inputed in a decimal not in a percentage")

  convert <- match.arg(arg=convert, choices=c("discrete","continuous"))
  switch(convert,
         discrete= exp(int)-1,
         continuous=log(1+int)
  )
}


#' Calculate market Bond Price from Yield-to-Maturity
#'
#' This function calculates the price of a bond with a given tenor, yield-to-maturity, and coupon frequency.
#'
#' @param Tenor The number of years remaining until the bond reaches maturity.
#' @param YTM The yield-to-maturity of the bond.
#' @param Freq The frequency of coupon payments per year. If Freq=0, this indicates a zero-coupon bond.
#'
#' @return The price of the bond.
#'
#' @examples
#' YTM2Price(5, 0.03, 2)
#' YTM2Price(10, 0.05, 4)
#'
#' @export
YTM2Price <- function(Tenor, YTM, Freq){

  # If the frequency is zero, return the present value of the bond
  if(Freq == 0) return(1 / (1 + YTM)^Tenor)

  # Calculate the time steps
  dt = 1 / Freq

  # Create a sequence of time steps
  t <- seq(dt, Tenor, by = dt)

  # Calculate the cash flow of the bond
  CF <- c(rep(YTM / Freq, length(t) - 1), 1 + YTM / Freq)

  # Calculate the discount factor of the bond
  DF = ifelse(t %% 1 == 0, (1 + YTM / Freq) ^ (-t * Freq), 1 / (1 + YTM * t))

  # Return the present value of the bond
  return(sum(CF * DF))
}


#' WilstonW
#'
#' Calculates the Wilson's function with given parameters
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
WilsonW<-function(u,v,a,ufr) {

  W <- matrix(NA,length(u),length(v))

  for (i in 1:length(u)) {
    for (k in 1:length(v)) {
      W[i,k] <- exp(-ufr*(u[i]+v[k]))*(a*min(u[i],v[k])
                                       -0.5*exp(-a*max(u[i],v[k]))*(exp(a*min(u[i],v[k]))
                                                                    -exp(-a*min(u[i],v[k]))))
    }}
  return(W)
}

