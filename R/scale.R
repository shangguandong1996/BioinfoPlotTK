#' scale_min_max
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mt <- t(matrix(c(rep(0, 5), rnorm(5), rep(2,5)), nrow = 5))
#' scale_min_max(mt)
#'
scale_min_max <- function(data){

    scale_function <- function(x) {
        (x - min(x)) / (max(x)-min(x)) ^ as.logical(max(x)-min(x))
    }

    scale_data <- apply(data, 1, scale_function)
    scale_data_t <- t(scale_data)

    return(scale_data_t)

}


#' scale_zscale
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' mt <- t(matrix(c(rep(0, 5), rnorm(5), rep(2,5)), nrow = 5))
#' scale_zscale(mt)
#'
scale_zscale <- function(data){

    scale_function <- function(x) {
        (x - mean(x)) / (sd(x)^as.logical(mean(x)))
    }

    scale_data <- apply(data, 1, scale_function)
    scale_data_t <- t(scale_data)

    return(scale_data_t)


}
