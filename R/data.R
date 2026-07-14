# Example data documentation --------------------------------------------------

#' Weekly English Channel crossings (UK example)
#'
#' A weekly count series of irregular small-boat crossings of the English
#' Channel towards the United Kingdom, beginning in early 2018. The series has
#' frequent zeros in its early weeks, which makes it a useful example for
#' zero inflation with the Poisson model. The series is used in the
#' irregular-migration application of Zens and Bijak (2026).
#'
#' @format A data frame with 376 rows and 2 variables:
#' \describe{
#'   \item{week}{ISO week label, e.g. `"2018-W01"`.}
#'   \item{count}{Non-negative integer count of weekly crossings.}
#' }
#' @source Weekly aggregates of detected irregular English Channel crossings,
#'   compiled from operational/agency records as described in Zens and Bijak
#'   (2026), \doi{10.1214/26-AOAS2171}.
#' @references
#' Zens, G. and Bijak, J. (2026). Dynamic Count Models with Flexible Innovation
#' Processes for Irregular Maritime Migration. \emph{The Annals of Applied
#' Statistics}. \doi{10.1214/26-AOAS2171}.
#' @examples
#' data(uk_weekly)
#' plot(uk_weekly$count, type = "h")
#' mean(uk_weekly$count == 0)
"uk_weekly"

#' Weekly Mediterranean crossings (Mediterranean example)
#'
#' A longer weekly count series of irregular sea crossings on the Mediterranean
#' route, beginning in autumn 2015. The counts are larger in magnitude with
#' fewer zeros than [uk_weekly]. The series is used in the irregular-migration
#' application of Zens and Bijak (2026).
#'
#' @format A data frame with 494 rows and 2 variables:
#' \describe{
#'   \item{week}{ISO week label, e.g. `"2015-W40"`.}
#'   \item{count}{Non-negative integer count of weekly crossings.}
#' }
#' @source Weekly aggregates of detected irregular Mediterranean crossings,
#'   compiled from operational/agency records as described in Zens and Bijak
#'   (2026), \doi{10.1214/26-AOAS2171}.
#' @references
#' Zens, G. and Bijak, J. (2026). Dynamic Count Models with Flexible Innovation
#' Processes for Irregular Maritime Migration. \emph{The Annals of Applied
#' Statistics}. \doi{10.1214/26-AOAS2171}.
#' @examples
#' data(med_weekly)
#' summary(med_weekly$count)
"med_weekly"
