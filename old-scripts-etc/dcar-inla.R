require(INLA)

tcar_inla <- function(
                      cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                              "log.prior", "quit"),
                      theta = NULL) {

    interpret.theta <- function() {
        list(tau = exp(theta[1L]), rho = 1 / (1 + exp(-theta[2L])))
    }

    graph <- function() {
        require(Matrix)
        as(diag(1, nrow(B)) + B, 'sparseMatrix')
    }

    Q <- function() {
        require(Matrix)
        param <- interpret.theta()
        as(param$tau * (Diagonal(x = colSums(B)) - param$rho * B), 'sparseMatrix')
    }

    mu <- function() rep(0, nrow(B))

    log.norm.const <- function() {
        require(Matrix)
        param <- interpret.theta()
        Q <- as(param$tau * (Diagonal(x = colSums(B)) - param$rho * B), 'sparseMatrix')

        res <- nrow(B) * (-0.5 * log(2 * pi)) +
            0.5 * determinant(Q, logarithm = TRUE)$modulus

        return(res)
    }

    log.prior <- function() {
        param = interpret.theta()

        dgamma(param$tau, 1, 5e-05, log = TRUE) + log(param$tau) +
            log(1) + log(param$rho) + log(1 - param$rho)
    }

    initial <- function() {
        return(c(0, 0))
    }

    quit <- function() invisible()

    res <- do.call(match.arg(cmd), args = list())
    return(res)
}

## requires k, A, and dmax
dcar_inla <- function(cmd = c(
                          "graph", "Q", "mu", "initial",
                          "log.norm.const", "log.prior", "quit"
                      ),
                      theta = NULL) {

    interpret.theta <- function() {
        return(list(
            tau = exp(theta[1L]),
            rho = 1 / (1 + exp(-theta[2L])),
            x0n = 1 / (1 + exp(-theta[3L]))
        ))
    }

    graph <- function() {
        require(Matrix)

        ret <- matrix(1, nrow = nrow(A), ncol = nrow(A))
        ret[A > dmax] <- 0
        as(ret, 'sparseMatrix')
    }

    C <- function(x0n) {
        ## rad <- max((200 - 50) * x0n + 50 + x1, 0)
        rad <- (dmax - 50) * x0n + 50
        ret <- 1 - 1 / (1 + exp(-k * (A - rad)))
        diag(ret) <- 0
        ret[A > dmax] <- 0
        return(as(ret, 'sparseMatrix'))
    }

    Q <- function() {
        require(Matrix)

        param <- interpret.theta()
        return(param$tau * (Diagonal(x = colSums(C(param$x0n))) -
                            param$rho * C(param$x0n)))
    }

    mu <- function() rep(0, nrow(A))

    log.norm.const <- function() {
        param <- interpret.theta()
        Q <- param$tau * (Diagonal(x = colSums(C(param$x0n))) -
                          param$rho * C(param$x0n))

        res <- nrow(A) * (-0.5 * log(2 * pi)) +
            0.5 * determinant(Q, logarithm = TRUE)$modulus

        return(res)
    }

    log.prior <- function() {
        param <- interpret.theta()

        res <- dgamma(param$tau, 1, .00005, log = TRUE) + log(param$tau) +
            log(1) + log(param$rho) + log(1 - param$rho) +
            ## dbeta(param$x0n, 1, 2, log = TRUE) + 
            log(1) + log(param$x0n) + log(1 - param$x0n)
        ## dnorm(param$x1, 0, dmax / 4, log = TRUE)

        return(res)
    }

    initial <- function() c(0, 0, 0)

    quit <- function() invisible()

    return(do.call(match.arg(cmd), args = list()))
}

dcar_tmarg <- function(inla_fit, dmax) {
    marg_tau <- inla.tmarginal(exp, inla_fit$marginals.hyperpar[[1]])
    marg_rho <- inla.tmarginal(
        function(x) 1 / (1 + exp(-x)),
        inla_fit$marginals.hyperpar[[2]]
    )
    marg_x0 <- inla.tmarginal(
        function(x) (dmax - 50) * (1 / (1 + exp(-x))) + 50,
        inla_fit$marginals.hyperpar[[3]]
    )

    return(list(tau = marg_tau, rho = marg_rho, x0 = marg_x0))
}

tcar_tmarg <- function(inla_fit) {
    marg_tau <- inla.tmarginal(exp, inla_fit$marginals.hyperpar[[1]])
    marg_rho <- inla.tmarginal(
        function(x) 1 / (1 + exp(-x)),
        inla_fit$marginals.hyperpar[[2]]
    )

    return(list(tau = marg_tau, rho = marg_rho))
}


plot_dcar_hyper <- function(inla_fit, dmax) {
    l <- dcar_tmarg(inla_fit, dmax)

    p1 <- ggplot(as_tibble(l[[1]]), aes(x = x, y = y)) +
        geom_line()
    p2 <- ggplot(as_tibble(l[[2]]), aes(x = x, y = y)) +
        geom_line()
    p3 <- ggplot(as_tibble(l[[3]]), aes(x = x, y = y)) +
        geom_line()

    grid.arrange(p1, p2, p3, nrow = 1)
}

plot_tcar_hyper <- function(inla_fit) {
    l <- tcar_tmarg(inla_fit)

    p1 <- ggplot(as_tibble(l[[1]]), aes(x = x, y = y)) +
        geom_line()
    p2 <- ggplot(as_tibble(l[[2]]), aes(x = x, y = y)) +
        geom_line()

    grid.arrange(p1, p2, nrow = 1)
}
