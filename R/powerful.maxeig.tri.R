#' @title Tridiagonal matrix maximal eigenpair
#' @description Calculate the maximal eigenpair for the tridiagonal matrix by
#' Thomas algorithm.
#'
#' @param a The lower diagonal vector.
#' @param b The upper diagonal vector.
#' @param C The main diagonal vector. 
#' @param digit.thresh The precise level of output results.
#'
#' @return A list of eigenpair object are returned, with components \eqn{z}, \eqn{v} and \eqn{iter}.
#' \item{z}{The approximating sequence of the maximal eigenvalue.}
#' \item{v}{The approximating eigenfunction of the corresponding eigenvector.}
#' \item{iter}{The number of iterations.}
#'
#' @examples
#' nn = 8
#' a = c(1:(nn - 1))^2
#' b = c(1:(nn - 1))^2
#' C = c(b[1], a[1:(nn - 2)] + b[2:(nn - 1)], a[nn - 1] + nn^2)
#' powerful.maxeig.tri(a, b, C, digit.thresh = 6)

#' @export
powerful.maxeig.tri = function(a, b, C, digit.thresh = 6) {
    
    N = length(a)
    Q = tridiag(b, a, -C)
    
    # check input vectors
    if (sum(Re(a * b) <= 0) > 0) {
        stop("The product of input vectors a and b should be all positive!")
    }
    if (sum(Im(C) != 0) > 0) {
        stop("Input vector c should be real!")
    }
    
    if (sum(C[2:N] - a[1:(N - 1)] - b[2:N] != 0) == 0 && C[1] == b[1]) {
        
        a_tilde = a
        b_tilde = b
        C_tilde = C
        m = 0
        
    } else {
        # print('Adjust Qtilde')
        m = max(max(c(-C[1] + abs(b[1]), -C[2:N] + abs(a[1:(N - 1)]) + abs(b[2:N]), 
            -C[N + 1] + abs(a[N]))), 0)
        # print(paste0('m=', m))
        u = a * b
        
        C_tilde = C + m
        
        b_tilde = rep(0, N)
        b_tilde[1] = C_tilde[1]
        
        a_tilde = rep(0, N)
        a_tilde[1] = a[1] * b[1]/b_tilde[1]
        for (k in 2:N) {
            b_tilde[k] = C_tilde[k] - a[k - 1] * b[k - 1]/b_tilde[k - 1]
            a_tilde[k] = a[k] * b[k]/b_tilde[k]
        }
        
    }
    
    if (C_tilde[N + 1] == a_tilde[N]) {
        stop(paste0("The maximal eigenpair is trivial! m=", m))
    }
    
    Q_sym = tridiag(sqrt(a * b), sqrt(a * b), -C_tilde)
    
    b_tilde = c(b_tilde[1:N], C_tilde[N + 1] - a_tilde[N])
    
    M = matrix(0, N + 1, N + 1)
    Phi = rep(0, N + 1)
    for (k in 1:N) {
        M[k, k] = 1
        M[k, (k + 1):(N + 1)] = cumprod(a_tilde[k:N]/b_tilde[k:N])
        Phi[k] = sum((M[k, ]/b_tilde)[k:(N + 1)])
    }
    M[N + 1, N + 1] = 1
    
    Phi[N + 1] = 1/b_tilde[N + 1]
    
    d = rep(1, N + 1)
    d[1] = 1
    d[2:(N + 1)] = sqrt(cumprod(a[1:k]/b[1:k]))
    
    z = list()
    rz = list()
    v = list()
    g = list()
    w = list()
    
    w[[1]] = sqrt(Phi)
    
    v[[1]] = w[[1]]/sqrt(sum(w[[1]]^2))
    
    xi = max(Re((v[[1]] %*% sqrt(M))/(b_tilde * v[[1]] - c(sqrt(a_tilde * 
        b_tilde[1:N]) * (v[[1]][2:(N + 1)]), 0))))
    
    z[[1]] = 1/xi
    
    rz[[1]] = round(z[[1]], digit.thresh)
    
    ratio_z = 1
    iter = 0
    
    while (ratio_z >= 10^(-digit.thresh)) {
        
        iter = iter + 1
        
        w = append(w, list(thomas.tri.sol(Q = Q_sym, v = v[[iter]], z = z[[iter]])))
        
        v = append(v, list(w[[iter + 1]]/sqrt(sum(w[[iter + 1]]^2))))
        
        g = append(g, list(diag(d) %*% (w[[iter + 1]]/sqrt(sum(w[[iter + 
            1]]^2)))))
        
        xi = max(Re((v[[iter + 1]] %*% sqrt(M))/(b_tilde * v[[iter + 1]] - 
            c(sqrt(a_tilde * b_tilde[1:N]) * v[[iter + 1]][2:(N + 1)], 0))))
        
        z = append(z, list(1/xi))
        
        ratio_z = abs(round(z[[iter + 1]], digit.thresh) - round(z[[iter]], 
            digit.thresh))
        
        rz[[iter + 1]] = round(z[[iter + 1]], digit.thresh)
        
    }
    
    
    if (ratio_z == 0) {
        v = v[-(iter + 1)]
        rz = rz[-(iter + 1)]
        
        iter = iter - 1
    }
    
    return(list(m = m, z = round(unlist(z), digit.thresh), v = v, g = g, 
        iter = iter))
}
