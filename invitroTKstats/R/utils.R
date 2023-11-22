#' Formatting function for X-axis in log10-scale
#'
#' @param x (Character) String to be formatted. 
#'
#' @return Text with desired expression. Replace any scientific e notation to ten notation, 
#' simplify 10^01 to 10 and 10^0 to 1.
#'
#' @importFrom scales scientific_format
#'
scientific_10 <- function(x) {
  out <- gsub("1e", "10^", scientific_format()(x))
  out <- gsub("\\+","",out)
  out <- gsub("10\\^01","10",out)
  out <- parse(text=gsub("10\\^00","1",out))
}

#' Heaviside
#'
#' @param x (Numeric) A numeric vector.
#'
#' @param threshold (Numeric) A threshold value used to compare to elements in \code{x}. (Defaults to 0.)
#'
#' @return A vector of 1 and 0 in which 1 indicates the element in \code{x} is larger or equal to the threshold.
#'
#'
Heaviside <- function(x, threshold=0)
{
  out <- rep(0,length(x))
  out[x >= threshold] <- 1
  return(out)
}

#' Convert a runjags-class object to a list
#'
#' @param runjagsdata.in (\code{runjags} Object) MCMC results from \code{autorun.jags}.
#'
#' @return A list object containing MCMC results from the provided runjags object.
#'
#'
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
