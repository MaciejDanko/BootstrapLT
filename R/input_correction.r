# Input lx ----------------------------------------------------------------

#' Convert a Life-table Survivorship Function into a Pace-Shape Object
#'
#' Given an age vector and corresponding survival probabilities a complete
#' life-table is calculated and a pace-shape object constructed.
#'
#' @param x Start of the age interval.
#' @param lx Probability to survive up until age x.
#' @param nax Subject-time alive in [x, x+n) for those who die in same interval
#'   (either numeric scalar, numeric vector or one of \code{c("udd", "cfm")}).
#' @param nx Width of age interval [x, x+n) (either numeric scalar, numeric
#'   vector or \code{"auto"}).
#' @param last_open Is the last age group open (TRUE) or closed (FALSE, default)?
#' @param time_unit The unit of the ages (by default "years").
#' @param messages Mesages/warnings ON (TRUE) or OFF (FALSE).
#' @section nx handling: For nx you may provide either a numeric scalar, a
#'   numeric vector, or let the function determine the width for you
#'   (\code{"auto"}, default). A scalar will be recycled for each age group. A
#'   vector must be as long as the age vector and allows you to specify the
#'   width of each age group separately. By default, the width of the age
#'   groups are calculated from differencing the age vector.
#'
#' @section nax handling: nax may be provided as either a numeric scalar, a
#'   numeric vector, or calculated via the uniform distribution of deaths
#'   (\code{udd}) method (default) or the constant force of mortality assumption
#'   (option \code{cfm}).
#'
#'   The (\code{udd}) method assumes a linear decline of the l(x) function over
#'   the width of an age group, implying that those who die in that age group
#'   die on average halfway into it (also known as the "midpoint" assumption):
#'
#'   nax = n/2 (see Preston et al. 2001, p. 46)
#'
#'   Assuming the mortality rate during age interval [x, x+n) to be constant
#'   (\code{cfm} method) implies an exponentially declining l(x) function
#'   within [x, x+n) and will produce nax values smaller than those calculated
#'   via the \code{udd} method. Preston et al. (2001, p. 46) provide an expression
#'   for nax given the assumption of constant mortality. Restating this expression
#'   in terms of nqx and npx leads to:
#'
#'   nax = -n/nqx - n/log(npx) + n
#'
#'   If the last age group is open and the \code{udd} or \code{cfm} method is
#'   used, then the last nmx value is log-linearly extrapolated based on the
#'   preceding two nmx and the nax, and the ex for the last age group id
#'   calculated using the constant hazard assumption.
#'
#' @source Preston, Samuel H., Patrick Heuveline, and Michel Guillot (2001).
#'   Demography: Measuring and modeling population processes. Oxford: Blackwell.
#'
#' @return A pace-shape object.
#'
#' @examples
#' Inputlx(x = prestons_lx$x, lx = prestons_lx$lx)
#' # open last age group
#' Inputlx(x = prestons_lx$x, lx = prestons_lx$lx, last_open = TRUE)
#' # different nax assumptions
#' Inputlx(x = prestons_lx$x, lx = prestons_lx$lx, nax = "cfm")
#'
#' @export
Inputlx <- function (x, lx,
                     nax = "udd",
                     nx = "auto",
                     last_open = FALSE,
                     time_unit = "years",
                     messages = TRUE) {

  # Input validation --------------------------------------------------------

  # k: number of age-groups
  k = length(x)

  # validate parameters
  pash:::ValidateAge(x)
  val_nx = pash:::Validatenx(nx, x, last_open)
  nx_ = val_nx[["nx"]]
  val_nax = pash:::Validatenax(nax, x, nx, last_open)
  nax_ = val_nax[["nax"]]

  # validate data
  pash:::Validatelx(lx)

  # Build life-table --------------------------------------------------------

  # set radix to 1
  lx_ = lx / lx[1L]

  # ndx: life-table deaths in age group [x, x+n)
  ndx = c(lx_[-k] - lx_[-1L], lx_[k])
  # nqx: life-table probability of death in age group [x, x+n)
  # given survival to age x
  nqx = ndx/lx_
  # npx: life-table probability of surviving age group [x, x+n)
  # given survival to age x
  npx = 1-nqx

  # nax: amount of subject-time at risk in age group [x, x+n)
  # contributed by those who die in that age group
  if (identical(val_nax[["nax_mode"]], "udd")) {
    nax_ = pash:::naxUDD(nx_, k, last_open)
  }
  if (identical(val_nax[["nax_mode"]], "cfm")) {
    nax_ = pash:::naxCFMfromnqx(x, nx_, nqx, npx, k, last_open)
  }

  # nLx: amount of subject-time at risk in age group [x, x+n)
  nLx = nx_*(lx_-ndx) + nax_*ndx
  # nmx: life-table mortality rate in age group [x, x+n]
  nmx = ndx/nLx
  if (identical(last_open, TRUE)) {
    # in case of an open age group nmx will be NA for this age group (as the
    # width of the open age group is unknown). therefore we log-linearly
    # extrapolate nmx based on the preceding two nmx and calculate nLx of the
    # last age group using the constant hazard assumption.
    if (val_nax[["nax_mode"]] %in% c("udd", "cfm")) {
      nmx[k] = pash:::LinearExtrapolation(x = x[c(k-2, k-1)], y = nmx[c(k-2, k-1)],
                                   xextra = x[k], loga = TRUE)
      nax_[k] = 1/nmx[k]
      nLx[k] = nax_[k]*lx_[k]
      if (messages) message("Inputlx() and last_open = TRUE: nmx of open age group log-linearly extrapolated based on preceding two nmx.")
    }
    if (val_nax[["nax_mode"]] %in% c("scalar", "vector")) {
      nLx[k] = nax_[k]*lx_[k]
      nmx[k] = ndx[k]/nLx[k]
    }
  }
  # Tx: amount of subject-time at risk above age group [x, x+n)
  Tx = rev(cumsum(rev(nLx)))
  # ex: life expectancy at age x
  ex = Tx/lx_
  # when lx becomes 0 ex becomes NaN. set it to 0
  ex[is.nan(ex)] = 0

  # Construct pash object ---------------------------------------------------

  # construct the pace-shape object, a validated life-table
  pash =  pash:::ConstructPash(
    x = x, nx = nx_, nmx = nmx, nax = nax_,
    nqx = nqx, npx = npx, lx = lx_, ndx = ndx,
    nLx = nLx, Tx = Tx, ex = ex,
    time_unit = time_unit, nax_mode = val_nax[["nax_mode"]], last_open = last_open,
    type = "lx",
    input = list(x = x, lx = lx, nax = nax, nx = nx)
  )

  return(pash)

}
