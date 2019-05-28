#' @title Solve the linear equation (-Q-zI)w=v.
#' @description Construct the solution of linear equation (-Q-zI)w=v.
#'
#' @param Q The given tridiagonal matrix.
#' @param v The column vector on the right hand of  equation.
#' @param z The Rayleigh shift.
#' @return A solution sequence \eqn{w} to the equation (-Q-zI)w=v.
#' 
#' @examples
#' nn = 8
#' a = c(1:(nn - 1))^2
#' b = c(1:(nn - 1))^2
#' C = c(b[1], a[1:(nn - 2)] + b[2:(nn - 1)], a[nn - 1] + nn^2)
#' Q = tridiag(b, a, -C)
#' zstart = 6
#' thomas.tri.sol(Q, z=zstart, v=rep(1,dim(Q)[1]))

#' @export
thomas.tri.sol = function(Q, v, z) {
    
    N = dim(Q)[1] - 1
    indx <- seq.int(N)
    a = Q[cbind(indx + 1, indx)]
    b = Q[cbind(indx, indx + 1)]
    C = -diag(Q)
    
    d = rep(0, N)
    d[1] = b[1]/(z - C[1])
    for (i in 2:N) {
        d[i] = b[i]/(z - C[i] - a[i - 1] * d[i - 1])
    }
    
    psi = rep(0, N + 1)
    psi[1] = v[1]/(C[1] - z)
    for (i in 2:(N + 1)) {
        psi[i] = (v[i] + a[i - 1] * psi[i - 1])/(C[i] - z + a[i - 1] * d[i - 
            1])
    }
    
    w = rep(0, N + 1)
    w[N + 1] = psi[N + 1]
    for (i in N:1) {
        w[i] = psi[i] - d[i] * w[i + 1]
    }
    
    return(w)
}
