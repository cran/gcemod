#' Head-and-neck cancer competing-risks example dataset
#'
#' An example head-and-neck cancer competing-risks dataset used to illustrate the
#' \pkg{gcemod} functions. The data are provided for illustration only and do not
#' represent real patients. The event of interest is cancer recurrence and the
#' competing event is death without recurrence; recurrence is the more frequent
#' event (recurrence-dominant, \eqn{\omega^{+} > 1}). Contrast with
#' \code{\link{prostate}}.
#'
#' @format A data frame with 1000 rows and 11 variables:
#' \describe{
#'   \item{id}{sequential subject identifier}
#'   \item{time}{follow-up time in years}
#'   \item{status}{competing-risks status: 0 = censored, 1 = cancer recurrence
#'     (event of interest), 2 = death without recurrence (competing event)}
#'   \item{age}{age at diagnosis (years)}
#'   \item{ps}{ECOG performance status (factor: 0, 1, 2)}
#'   \item{female}{sex indicator (1 = female)}
#'   \item{smoker}{smoking indicator (1 = more than 10 pack-years)}
#'   \item{t_cat}{T category (factor: 1-4)}
#'   \item{n_cat}{N category (factor: 0-3)}
#'   \item{p16}{p16 status (1 = positive)}
#'   \item{site}{primary site (factor: oropharynx, larynx, hypopharynx, oralcavity)}
#' }
#' @keywords datasets
"hn"

#' Prostate cancer competing-risks example dataset
#'
#' An example prostate cancer competing-risks dataset used to illustrate the
#' \pkg{gcemod} functions. The data are provided for illustration only and do not
#' represent real patients. The event of interest is distant metastasis or
#' prostate-cancer death and the competing event is death from other causes; the
#' competing event is the more frequent (competing-death-dominant,
#' \eqn{\omega^{+} < 1}). Contrast with \code{\link{hn}}.
#'
#' @format A data frame with 1000 rows and 9 variables:
#' \describe{
#'   \item{id}{sequential subject identifier}
#'   \item{time}{follow-up time in years}
#'   \item{status}{competing-risks status: 0 = censored, 1 = distant metastasis or
#'     prostate-cancer death (event of interest), 2 = death from other causes
#'     (competing event)}
#'   \item{age}{age at diagnosis (years)}
#'   \item{ps}{performance status (factor: 0, 1)}
#'   \item{psa}{pre-treatment PSA (ng/mL, capped at 20)}
#'   \item{gleason}{Gleason grade group (factor: "<=6", "3+4", "4+3", ">=8")}
#'   \item{t2b}{stage indicator (1 = T2b or higher)}
#'   \item{comorbidity}{significant comorbidity indicator (1 = yes)}
#' }
#' @keywords datasets
"prostate"
