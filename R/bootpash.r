library(pash)  #remove when integrated with pash
library(roxygen2)
#source("./R/input_correction.r") #remove when integrated with pash, and Inputlx corrected

#' Calculate midd-classes of an age interval
#' 
#' @description 
#' Calculate midd-classes of an age interval. \cr\cr
#' \emph{\bold{Internal function}}
#' @param x Vector with begining of age classes.
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
getInterval<-function(x) x+c(diff(x),diff(x)[length(x)-1])/2

#' Converting dx to lx assuming that all events are exactly observed (no censoring)
#'
#' @description 
#' Converting dx to lx assuming that all events are exactly observed (no censoring). \cr\cr
#' \emph{\bold{Internal function}}
#' @param dx Vector with death counts
#' @return A vector with population size at the beginnig of interval.
#' @details 
#' The presence or absence of open interval does not matter here.
#' @seealso \code{\link{lx2dx}}, \code{\link{dx2age}}, and \code{\link{age2dx}}.
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
dx2lx<-function(dx) sum(dx)-c(0,cumsum(dx))[seq_along(dx)]

#' Converting lx to dx assuming that all events are exactly observed (no censoring)
#'
#' @description Converting lx to dx assuming that all events are exactly observed (no censoring). \cr\cr
#' \emph{\bold{Internal function}}
#' @param lx population size at the beginning of the age interval.
#' @return A vector with death counts. 
#' @seealso \code{\link{dx2lx}}, \code{\link{dx2age}}, and \code{\link{age2dx}}.
#' @keywords internal
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
lx2dx<-function(lx) c(-diff(lx),lx[length(lx)])

#' Constructing vector of age at deaths assuming that each death occurs in the middle of an age interval
#'
#' @description
#' The function converts counts in vector dx into into individual age/time at deaths.
#' Notice that dx does not contain information about censoring so resulting vector has only exactly observed events.\cr\cr
#' \emph{\bold{Internal function}}
#' @param dx Death counts vector.
#' @param x Beginning of the Age/time class. A vector of the same length as \code{ndx}.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @seealso \code{\link{lx2dx}}, \code{\link{dx2lx}}, and \code{\link{age2dx}}.
#' @examples
#' \dontrun{
#' # Data:
#' x=c(0,0.5,2,5,10,13,15)
#' dx=c(1,2,6,15,22,6,1)
#'
#' # Should give the same vector as in dx.
#' age2dx(dx2age(dx,x),x)
#'
#' }
#' @keywords internal
dx2age<-function(dx,x) unlist(sapply(seq_along(x),function (kk) rep(getInterval(x)[kk],dx[kk])))

#' Constructing dx from individual ages/times at deaths
#'
#' @description
#' The function converts ages/times at death into dx according.
#' Notice that age at deaths are assumed to be exactly observed.\cr\cr
#' \emph{\bold{Internal function}}
#' @param times vector with age/times at deaths.
#' @param x Beginning of the Age/time class. The resulting vector of \code{dx} will have the same length.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @seealso \code{\link{lx2dx}}, \code{\link{dx2lx}}, and \code{\link{dx2age}}.
#' @examples
#' \dontrun{
#' # Data:
#' x=c(0,0.5,2,5,10,13,15)
#' dx=c(1,2,6,15,22,6,1)
#'
#' # Should give the same vector as in dx.
#' age2dx(dx2age(dx,x),x)
#'
#' }
#' @keywords internal
age2dx <- function(times,x) (hist(x = times, plot = FALSE, right = FALSE, breaks = c(x, x[length(x)] + 1000))$counts)

#' Geting pace and shape measures for bootstrap computations
#'
#' @description Geting pace and shape measures for bootstrap computations. \cr\cr
#' \emph{\bold{Internal function}}
#' @param dx A vector with integer counts (deaths).
#' @param x Beginning of the Age/time classes. Vector with the same length as \code{dx}
#' @param pash.parent A parent pash object.
#' @param pace.type Which pace measure should be returned (default "all")?
#' Use "none" if you don't wat to return any pace measure. See \code{\link{GetPace}} for details.
#' @param shape.type Which shape measure should be returned (default "all")?
#' Use "none" if you don't wat to return any shape measure. See \code{\link{GetShape}} for details.
#' @param q GetPace parameter. Quantile specification for age where q percent of the life-table population is still alive (defaults to median).
#' See \code{\link{GetPace}} for details.
#' @param harmonized GetShape parameter. Should the harmonized version of the shape measures be returned (default \code{TRUE})?
#' See \code{\link{GetShape}} for details.
#' @keywords internal
getpash <- function(dx, x, pash.parent, pace.type = "all", shape.type = "all", q = 0.5, harmonized = TRUE){
  if (missing(x)) x <- pash.parent$lt$x
  nax <- pash.parent$lt$nax
  lx <- dx2lx(dx)
  obj <- Inputlx(x = x, lx = lx, last_open = attr(pash.parent,'last_open'), time_unit = attr(pash.parent,'time_unit'), messages = FALSE)
  if (pace.type == "none") p <- NULL else p <- GetPace(obj, q = q, type = pace.type)
  if (shape.type == 'none') s <- NULL else s <- GetShape(obj, harmonized = harmonized, type = shape.type)
  return(c(p,s))
}

#' Fast JackKnife method performed on \code{dx}
#' 
#' @description Fast JackKnife method performed on \code{dx}.\cr\cr
#' \emph{\bold{Internal function}}
#' @param dx A vector with integer counts (deaths).
#' @param x Beginning of the Age/time classes. Vector with the same length as \code{dx}.
#' @param ... Other parameters passed to \code{\link{getpash}} function.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall.
#' @keywords internal
JackKnife_dx <- function(dx, x, ...) sapply(x[dx>0], function(kk) getpash(dx = dx - (x == kk), x = x, ...))

#' Fast Bootstrap method performed on \code{dx}
#' @description Fast Bootstrap method performed on \code{dx}.\cr\cr
#' \emph{\bold{Internal function}}
#' @param x Beginning of the Age/time class.
#' @param dx A vector with integer counts (deaths).
#' @param pash.parent A parent pash object.
#' @param N Number of bootstrap replicates.
#' @param pace.type Which pace measure should be returned (default "all")?
#' Use "none" if you don't wat to return any pace measure. See \code{\link{GetPace}} for details.
#' @param shape.type Which shape measure should be returned (default "all")?
#' Use "none" if you don't wat to return any shape measure. See \code{\link{GetShape}} for details.
#' @param q GetPace parameter. Quantile specification for age where q percent of the life-table population is still alive (defaults to median).
#' See \code{\link{GetPace}} for details.
#' @param harmonized GetShape parameter. Should the harmonized version of the shape measures be returned (default \code{TRUE})?
#' See \code{\link{GetShape}} for details.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall.
#' @keywords internal
Bootstrap_dx <- function(dx, 
                         x, 
                         pash.parent, 
                         N = 1000, 
                         pace.type = "all", 
                         shape.type = "all", 
                         q = 0.5, 
                         harmonized = TRUE) 
  replicate(n = N, expr = suppressWarnings(getpash(dx = age2dx(sample(x = x, 
                                                     size = sum(dx), 
                                                     replace = TRUE, 
                                                     prob = dx / sum(dx)), 
                                              x = x),
                                  x = x,
                                  pash.parent = pash.parent,
                                  pace.type = pace.type,
                                  shape.type = shape.type,
                                  q = q,
                                  harmonized = harmonized)))

#' Calculation of JeckKnife acceleration parameter for BCA method
#' 
#' @description Calculation of JeckKnife acceleration parameter for BCA method. \cr\cr
#' \emph{\bold{Internal function}}
#' @param JackKnifeMat JackKnife matrix with one row per each considered pash measure. See \code{\link{JackKnife_dx}}.
#' @param dx A vector with integer counts (deaths).
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall. Page 186.
#' @examples
#' \dontrun{
#' 
#' #Get some data
#' population.size <- 10000
#' obj <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,nax = australia_10y$nax, nx = australia_10y$nx, last_open = TRUE)
#' dx <- lx2dx(round(object$lt$lx*population.size))
#' x <- obj$lt$x
#' 
#' #Acceleration parameters calculated by slow method
#' times <- dx2age(dx,x)
#' gp <- function(times, x, pash.parent) getpash(dx = age2dx(times = times, x=x), x=x, pash.parent = pash.parent)
#' J1 <- as.matrix(sapply(seq_along(times), function(k) gp(times = times[-k], x = x, pash.parent = obj)))
#' a.slow <- JackKnifeAcc(J1)
#' 
#' #Acceleration parameters calculated by fast method
#' J2 <- JackKnife_dx(dx = org.dx, x = org.x, pash.parent = obj)
#' a.fast <- JackKnifeAcc(J2,dx)
#' 
#' #Compare results
#' cbind(a.slow = a.slow, a.fast = a.fast, diff = round(a.slow - a.fast, 12))
#' 
#' }
#' @keywords internal
JackKnifeAcc <- function(JackKnifeMat, dx){
  JackKnifeMat <- as.matrix(JackKnifeMat)
  if (missing(dx)){ 
    # Classic approach
    # Jackknife constructed from times, very slow method for big populations
    L <- rowSums(JackKnifeMat, na.rm = TRUE) / dim(JackKnifeMat)[2] - JackKnifeMat
    a <- rowSums(L^3, na.rm = TRUE) / (6 * rowSums(L^2, na.rm = TRUE)^1.5)
  } else { 
    # Eficient method to calculate Jacknife a
    # jackknife constructed from dx 
    dx <- dx[dx > 0] #JackKniefe omitted this values during its construction
    DX <- t(matrix(dx, length(dx), dim(JackKnifeMat)[1]))
    L <- rowSums(DX * JackKnifeMat, na.rm = TRUE) / sum(dx) - JackKnifeMat
    Lsq <- DX*(L^2)
    Lcu <- DX*(L^3)
    a = rowSums(Lcu, na.rm = TRUE) / (6 * rowSums(Lsq, na.rm = TRUE)^1.5)
  }
  a[!is.finite(a)] <- 0
  a
}

#' Checking if all elements in the vector are equal
#' @description Checking if all elements in the vector are equal. .\cr\cr
#' \emph{\bold{Internal function}}
#' @param x A vector be analized.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @keywords internal
all.elements.equal <- function(x) round(sum(abs(diff(x))),floor(-log10(.Machine$double.eps^0.8))) == 0

#' Calcualtion of the BCA type of confidence intervals
#' 
#' @description Calcualtion of the BCA type of confidence intervals. \cr\cr
#' \emph{\bold{Internal function}}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall.
#' @keywords internal
BCA.CI <- function(OrgEst, BootEst, JackKnifeMat, Orgdx, alpha = 0.05) {
  av <- JackKnifeAcc(JackKnifeMat = JackKnifeMat, dx = Orgdx)
  rna <- rownames(BootEst)
  result <- matrix(NA, length(av), 2)
  for (j in seq_along(av)){
    a <- av[j]
    be <- BootEst[j,]
    be <- be[!is.na(be)]
    if (!all.elements.equal(be)) {
      oe <- OrgEst[j]
      B <- length(be)
      zalpha <- stats:::qnorm(alpha / 2)
      z0 <- stats:::qnorm(sum(be <= oe)/B)
      a1 <- stats:::pnorm(z0 + (z0 + zalpha) / (1 - a * (z0 + zalpha)))
      a2 <- stats:::pnorm(z0 + (z0 - zalpha) / (1 - a * (z0 - zalpha)))
      CI <- stats:::quantile(be, c(a1,a2), na.rm = T)
    } else CI <- c(NA, NA)
    result[j,] <- CI
  }
  ind <- result[,1] > OrgEst
  ind[is.na(ind)] <- FALSE
  if (sum(ind) > 0) warning(paste('BCa confidence intervals do not cover original estimate in ', rna[ind], '. The lower bound was adjusted.\n',sep=''))
  result[ind,1] <- OrgEst[ind]
  ind <- result[,2]< OrgEst
  ind[is.na(ind)] <- FALSE
  if (sum(ind) > 0) warning(paste('BCa confidence intervals do not cover original estimate in ', rna[ind], '. The upper bound was adjusted.\n',sep=''))
  result[ind,2] <- OrgEst[ind]
  rownames(result) <- rna
  colnames(result) <- c('Lower bound','Upper bound')
  as.data.frame(result)
}

#' Calcualtion of the Percentile type of confidence intervals
#'
#' @description Calcualtion of the Percentile type of confidence intervals. \cr\cr
#' \emph{\bold{Internal function}}
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall.
#' @keywords internal
PER.CI <- function(OrgEst, BootEst, alpha=0.05) {
  rna <- rownames(BootEst)
  result <- matrix(NA, dim(BootEst)[1], 2)
  for (j in seq_len(dim(BootEst)[1])){
    be <- BootEst[j,]
    be <- be[!is.na(be)]
    if (!all.elements.equal(be)) {
      conf <- 1-alpha
      zalpha <- (1 + c(-conf, conf))/2
      CI <- stats:::quantile(be, c(zalpha[1], zalpha[2]), na.rm = TRUE)
    } else CI <- c(NA, NA)
    result[j,1:2] <- CI
  }
  ind <- result[,1]>OrgEst
  ind[is.na(ind)] <- FALSE
  if (sum(ind)>0) warning(paste('Percentile confidence intervals do not cover original estimate in ',rna[ind],'. The lower bound was adjusted.\n',sep=''))
  
  result[ind,1] <- OrgEst[ind]
  ind <- result[,2] < OrgEst
  ind[is.na(ind)] <- FALSE
  if (sum(ind)>0) warning(paste('Percentile confidence intervals do not cover original estimate in ',rna[ind],'. The upper bound was adjusted.\n',sep=''))
  
  result[ind,2] <- OrgEst[ind]
  rownames(result) <- rna
  colnames(result) <- c('Lower bound','Upper bound')
  as.data.frame(result)
}

#' Default procedure to perform bootstrap and calculate confidence intervals of pash measures
#' 
#' @description Default procedure to perform bootstrap and calculate confidence intervals of pash measures. \cr\cr
#' \emph{\bold{Internal function}}
#' @param x Beginning of the Age/time class.
#' @param dx A vector with integer counts (deaths).
#' @param pash.parent A parent pash object.
#' @param N Number of bootstrap replicates.
#' @param js.er Maximal acceptable fraction of JackKnife errors.
#' @param bs.er Maximal acceptable fraction of Bootstrap errors.
#' @param pace.type Which pace measure should be returned (default "all")?
#' Use "none" if you don't wat to return any pace measure. See \code{\link{GetPace}} for details.
#' @param shape.type Which shape measure should be returned (default "all")?
#' Use "none" if you don't wat to return any shape measure. See \code{\link{GetShape}} for details.
#' @param q GetPace parameter. Quantile specification for age where q percent of the life-table population is still alive (defaults to median).
#' See \code{\link{GetPace}} for details.
#' @param harmonized GetShape parameter. Should the harmonized version of the shape measures be returned (default \code{TRUE})?
#' See \code{\link{GetShape}} for details.
#' @param trace Logical indicating if to show summary of performed bootstrap.
#' @param alpha Significance level.
#' @author Maciej J. Danko <\email{danko@demogr.mpg.de}> <\email{maciej.danko@gmail.com}>
#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall.
#' @examples 
#' \dontrun{
#' 
#' #Get some data
#' population.size=105
#' obj <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,nax = australia_10y$nax, nx = australia_10y$nx, last_open = TRUE)
#' dx <- lx2dx(round(object$lt$lx*population.size))
#' x <- obj$lt$x
#' 
#' #Calcualte confidence intervals usung bootstrap
#' z <- boot.default(dx = dx, x = x, pash.parent = obj, trace = TRUE)
#' z
#' 
#' }
#' @keywords internal
boot.default <- function(dx, 
                         x, 
                         pash.parent, 
                         N = 1000, 
                         bs.er = 0.25, 
                         jk.er = 0.10, 
                         pace.type = "all", 
                         shape.type = "all", 
                         q = 0.5, 
                         harmonized = TRUE,
                         trace = TRUE,
                         alpha = 0.05){
  
  # Detect mistakes
  .AnalizeBoot <- function(jn, bt, trace, ind0){
    proc_fails_jn <- rowSums(is.na(jn)) / dim(jn)[2]
    proc_fails_bt <- rowSums(is.na(bt)) / dim(bt)[2]
    ind <- (proc_fails_jn >= jk.er) | (proc_fails_bt >= bs.er) | ind0
    ind2 <- (proc_fails_bt >= bs.er) | ind0
    vec <- seq_len(dim(bt)[1])
    if (trace) {
      printo <- round(100 * rbind(JackKnife = proc_fails_jn, Bootstrap = proc_fails_bt), 2)
      cat('\nBootstrap with', format(N, scientific = FALSE), 'replicates\n')
      cat('\nFraction [%] of cases where a particular measure could not be evaluated:\n')
      print(printo)
      if (sum(ind) > 0) cat('\nMeasures for witch BCa CI will be not evaluated:',
                            paste(rownames(bt)[vec[ind]], sep = '', collpase = ','), '\n')
      if (sum(ind2) > 0) cat('Measures for witch percentile CI will be not evaluated:',
                             paste(rownames(bt)[vec[ind2]], sep = '', collpase = ','), '\n')
    }
    ind <- vec[ind]
    ind2 <- vec[ind2]
    list(BCa = ind, P = ind2)
  }
  
  OrgPash <- getpash(dx = dx, 
                       x = x,
                       pash.parent = pash.parent,
                       pace.type = pace.type,
                       shape.type = shape.type,
                       q = q,
                       harmonized = harmonized)
  ind0 <- is.na(OrgPash)
  if (N > 1){
    jn <- JackKnife_dx(dx = dx, 
                    x = x,
                    pash.parent = pash.parent,
                    pace.type = pace.type,
                    shape.type = shape.type,
                    q = q,
                    harmonized = harmonized)
    bt <- Bootstrap_dx(dx=dx, 
                   x = x,
                   pash.parent = pash.parent,
                   pace.type = pace.type,
                   shape.type = shape.type,
                   q = q,
                   harmonized = harmonized)
 
    z <- .AnalizeBoot(jn = jn, bt = bt, trace = trace, ind0 = ind0)
    CI.Percentile <- as.matrix(PER.CI(OrgEst = OrgPash, 
                                      BootEst = bt,
                                      alpha = alpha))
    CI.BCa <- as.matrix(BCA.CI(OrgEst = OrgPash, 
                               BootEst = bt, 
                               JackKnifeMat = jn, 
                               Orgdx = dx,
                               alpha = alpha))
    CI.BCa[z$BCa,] <- NA
    CI.Percentile[z$P,] <- NA
    CI.BCa <- as.data.frame(CI.BCa)
    CI.Percentile <- as.data.frame(CI.Percentile)
  } else {
    CI.BCa <- NULL
    CI.Percentile <- NULL
  }
  return(list(Pash = OrgPash,
       CI.BCa = CI.BCa,
       CI.Percentile = CI.Percentile))
}


#' Claculates condfidence intervals for pash measures
#' 
#' @description Claculates condfidence intervals around pash measures for a given \code{pash} object. \cr\cr Generic function.
#' @param object \code{Pash} object.
#' @param population.size The size of population for the constructed pash object. If not given then it will be read from the \code{object}.
#' @param N Number of bootstrap replicates.
#' @param pace.type Which pace measure should be returned (default "all")?
#' Use "none" if you don't wat to return any pace measure. See \code{\link{GetPace}} for details.
#' @param shape.type Which shape measure should be returned (default "all")?
#' Use "none" if you don't wat to return any shape measure. See \code{\link{GetShape}} for details.
#' @param q GetPace parameter. Quantile specification for age where q percent of the life-table population is still alive (defaults to median).
#' See \code{\link{GetPace}} for details.
#' @param harmonized GetShape parameter. Should the harmonized version of the shape measures be returned (default \code{TRUE})?
#' See \code{\link{GetShape}} for details.
#' @param trace Logical indicating if to show summary of performed bootstrap.
#' @param alpha Significance level.
#' @param bs.er Maximal acceptable fraction of Bootstrap errors 
#' (a measure does not exists for a certain bootstraped data).
#' @param js.er Maximal acceptable fraction of JackKnife errors 
#' (a measure does not exists for a certain JackKnife cases).

#' @references 
#' Efron, B., & Tibshirani, R. J. (1993). An introduction to the bootstrap. New York: Chapman & Hall.
#' @examples 
#' \dontrun{
#' 
#' #Get a data
#' obj <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,nax = australia_10y$nax, nx = australia_10y$nx, last_open = TRUE)
#' 
#' ci1 <- confint(object = obj, population.size = 1000, trace = TRUE)
#' print(ci1, CI.type = 'BCa')
#' print(ci1, CI.type = 'Percentile')
#'
#' #Small population size
#' obj <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,nax = australia_10y$nax, last_open = FALSE) 
#' ci2 <- confint(object = obj, population.size = 300, trace = TRUE)
#' ci2
#' print(ci2)
#' 
#' #Very small population size
#' # There is a mistake in pash package. The linear extrapolation doesn't work for some sets with open last intervals.
#' obj <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,nax = australia_10y$nax, nx = australia_10y$nx, last_open = TRUE)
#' ci3 <- confint(object = obj, population.size = 10, shape.type = 'none', trace = TRUE)
#' 
#' }
#' @export
confint.pash<-function(object, 
                   population.size, 
                   N = 1000, 
                   pace.type = "all", 
                   shape.type = "all", 
                   q = 0.5, 
                   harmonized = TRUE,
                   trace = TRUE,
                   alpha = 0.05,
                   bs.er = 0.25, 
                   jk.er = 0.10){
  i.ndx <- attributes(object)$source$input$dx
  i.lx <- attributes(object)$source$input$lx
  i.x <- attributes(object)$source$input$x
  if ((length(i.ndx) == 0) && (length(i.lx) != 0)) i.ndx <- -diff(c(i.lx, 0))
  # If population size not given try to read from the pash attribute
  if (missing(population.size) && length(attr(object, 'population.size')>0)) population.size <- attr(object, 'population.size')
  # If population size is still unknown try to get it from source attribute
  if (missing(population.size)) {
    mps <- TRUE
    message('The parameter "population.size" was not given. It was retrieved from "source" attribute of pash object. Notice that it may not reflect truth (e.g. it could be equal to standardized population size). Population size is crucial for PCLM estimation.\n')
    if (length(i.ndx) == 0) stop('Population size cannot be retrieved from pash object.')
    population.size <- sum(i.ndx)
    if (population.size <= 1.01) stop('Population size cannot be retrieved from pash object.')
    message('Population size set to ', population.size, '\n')
  } else {
    mps <- FALSE
    if (!is.numeric(population.size)) stop('Please give a correct population size.')
  }
  
  res<-boot.default(dx = lx2dx(round(object$lt$lx*population.size)), 
                    x = object$lt$x, 
                    pash.parent=object, 
                    N = N, 
                    bs.er = bs.er, 
                    jk.er = jk.er, 
                    pace.type = pace.type,
                    shape.type = shape.type,
                    q = q,
                    harmonized = harmonized,
                    trace = trace,
                    alpha = alpha)
  res$pash.parent <- object
  res$population.size <- population.size
  res$q <-q
  res$harmonized <- harmonized
  res$alpha <- alpha
  class(res) <- 'bootpash'
  res
}

#' Printing the confidence intervals for pash object
#' 
#' @description  
#' Generic function to print the confidence intervals for pash object.
#' @param oject An \code{bootpash} object with fitted confidence interval for \code{pash} object.
#' @param CI.type The type of printed confidence intervals. One from "BCa" or "Percentile".
#' @export
print.bootpash <- function(object, CI.type = c('BCa', 'Percentile'), digits = 4){
  if (tolower(CI.type[1]) == 'bca') CI <- object$CI.BCa else if (tolower(CI.type[1]) == 'percentile') CI <- object$CI.Percentile else stop('Unknown CI type.')
  Mat <- cbind(lower = CI[,1], pash = object$Pash, upper = CI[,2])
  print(round(Mat, digits), quote = FALSE)
  invisible(Mat)
}
  
