
expSmooth <- function(vector,lambda=NULL, round=TRUE){

  if(is.null(lambda)) lambda=2/(length(vector)+1)

  n <- length(vector)
  weight <- c(1 ,rep(lambda,times=n-1))*(1-lambda)^rev(seq_along(vector)-1)
  avg <- sum(weight*vector)

  if(round) return(round(20*avg)/20)
  return(avg)
}

update_LTFR <- function(data,update_year,last_LTFR,BOK_target=2){

  if(!is.data.frame(data)) stop("data must be a data.frame")
  if(sum(c("INT_CD91_Y","CPI_YG")  %in% colnames(data))<2) stop("data must have 'year','INT_CD91_Y','CPI_YG' columns")

  data %>%
    select(year, INT_CD91_Y,CPI_YG) %>%
    filter(year>=1991,year<=update_year) %>%
    mutate(INT_REAL=(INT_CD91_Y-CPI_YG)/(1+CPI_YG/100)) %>%
    summarise(LTFR=expSmooth(INT_REAL)+BOK_target) %>%
    mutate(LTFR=case_when(LTFR-last_LTFR > 0.15 ~ last_LTFR + 0.15,
                          LTFR-last_LTFR >= -0.15 ~ last_LTFR,
                          TRUE ~ last_LTFR-0.15)) %>% pull

}
