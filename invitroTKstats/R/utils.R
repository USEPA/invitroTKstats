#' Formatting function for X-axis in log10-scale
#'
#' @param x
#'
#' @return text with desired expression
#'
#' @import scales
#'
#' @export calc_fup_uc
scientific_10 <- function(x) {
  out <- gsub("1e", "10^", scientific_format()(x))
  out <- gsub("\\+","",out)
  out <- gsub("10\\^01","10",out)
  out <- parse(text=gsub("10\\^00","1",out))
}

#' Heaviside
#'
#' @param x a vector
#'
#' @param threshold defaults to 0
#'
#' @return  a vector of 1 and 0, 1 indicates input element is above threshold
#'
#'
#' @export Heaviside
Heaviside <- function(x, threshold=0)
{
  out <- rep(0,length(x))
  out[x >= threshold] <- 1
  return(out)
}

#' Convert a runjags object to a list
#'
#' @param runjagsdata.in result from autorun.jags(), an object of class runjags
#'
#' @return a list contains data from the input runjags object
#'
#'
#' @export runjagsdata.to.list
runjagsdata.to.list <- function(runjagsdata.in)
{
  temp <- strsplit(runjagsdata.in,"\n")[[1]]
  list.out <- list()
  for (i in 1:length(temp))
  {
    temp2 <- strsplit(temp[i]," <- ")[[1]]
    list.out[[gsub("\"","",temp2[1])]] <- eval(parse(text=temp2[2]))
  }
  return(list.out)
}
