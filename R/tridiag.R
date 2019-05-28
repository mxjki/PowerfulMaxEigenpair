#' @title Tridiagonal matrix
#' @description Generate tridiagonal matrix Q based on three input vectors.
#'
#' @param upper The upper diagonal vector.
#' @param lower The lower diagonal vector.
#' @param main The main diagonal vector.
#' @return A tridiagonal matrix is returned.
#' 
#' @examples a = c(1:7)^2
#' b = c(1:7)^2
#' c = -c(1:8)^2
#' tridiag(b, a, c)

#' @export
tridiag <- function(upper, lower, main) {
    out <- matrix(0, length(main), length(main))
    diag(out) <- main
    indx <- seq.int(length(upper))
    out[cbind(indx + 1, indx)] <- lower
    out[cbind(indx, indx + 1)] <- upper
    return(out)
}
