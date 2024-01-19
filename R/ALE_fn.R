#' Accumulated Local Effects workhorse
#' 
#' @description Calculates the Accumulated Local Effects (ALE) from a given data.frame, model, predictions, and covariate. This riffs on the ALEPlot function available in the ALEPlot package. 
#'
#' @param X the data.frame to get the covariate from
#' @param X.model the name of the response variable
#' @param pred.fun a function to calculate new predictions from the model
#' @param J the column index of the covariate of interest
#' @param K an integer value that determines the number of "windows" or breaks to calculate the model predictions over. More increase computational time but serves smooths the ALE predictions.
#' @param type is either (response) or (prob) from predict.randomForest; if (prob) then n sets of predictions are returned for the n levels in var- if "response" then the factorized predicted response values are returned
#' @param multi is a logical for either multivariate factor as the response variable (TRUE) or not (FALSE- the default)
#' 
#' @return A list that contains:
#' \itemize{
#'      \item \strong{K}: the number of realized breaks
#'      \item \strong{x.values}: the break values trialed
#'      \item \strong{class}: the class of the covariate
#'      \item \strong{quantile}: the quantile of the breaks
#'      \item \strong{fJ}: the ALEs evaluated at a given x value
#' }
#' 
#' @export
#' 
#' 
ALE_fn <- function (X, X.model, pred.fun, J, K = 40, type='response', multi=FALSE){
    N = dim(X)[1]
    d = dim(X)[2]
    if (is(X[, J], "factor")) {
        X[, J] <- droplevels(X[, J])
        x.count <- as.numeric(table(X[, J]))
        x.prob <- x.count/sum(x.count)
        K <- nlevels(X[, J])
        D.cum <- matrix(0, K, K)
        D <- matrix(0, K, K)
        for (j in setdiff(1:d, J)) {
            if (is(X[, j], "factor")) {
                A = table(X[, J], X[, j])
                A = A/x.count
                for (i in 1:(K - 1)) {
                  for (k in (i + 1):K) {
                    D[i, k] = sum(abs(A[i, ] - A[k, ]))/2
                    D[k, i] = D[i, k]
                  }
                }
                D.cum <- D.cum + D
            }
            else {
                q.x.all <- quantile(X[, j], probs = seq(0, 1, 
                  length.out = 100), na.rm = TRUE, names = FALSE)
                x.ecdf = tapply(X[, j], X[, J], ecdf)
                for (i in 1:(K - 1)) {
                  for (k in (i + 1):K) {
                    D[i, k] = max(abs(x.ecdf[[i]](q.x.all) - 
                      x.ecdf[[k]](q.x.all)))
                    D[k, i] = D[i, k]
                  }
                }
                D.cum <- D.cum + D
            }
        }
        D1D <- cmdscale(D.cum, k = 1)
        ind.ord <- sort(D1D, index.return = T)$ix
        ord.ind <- sort(ind.ord, index.return = T)$ix
        levs.orig <- levels(X[, J])
        levs.ord <- levs.orig[ind.ord]
        x.ord <- ord.ind[as.numeric(X[, J])]
        row.ind.plus <- (1:N)[x.ord < K]
        row.ind.neg <- (1:N)[x.ord > 1]
        X.plus <- X
        X.neg <- X
        X.plus[row.ind.plus, J] <- levs.ord[x.ord[row.ind.plus] + 
            1]
        X.neg[row.ind.neg, J] <- levs.ord[x.ord[row.ind.neg] - 
            1]
        y.hat <- pred.fun(X.model = X.model, newdata = X, type = type)
        y.hat.plus <- pred.fun(X.model = X.model, newdata = X.plus[row.ind.plus, 
            ], type = type)
        y.hat.neg <- pred.fun(X.model = X.model, newdata = X.neg[row.ind.neg, 
            ], type = type)
        Delta.plus <- y.hat.plus - y.hat[row.ind.plus]
        Delta.neg <- y.hat[row.ind.neg] - y.hat.neg
        if (!multi) {
            if (type == "prob") {
                Delta <- as.numeric(tapply(c(Delta.plus[, 2], 
                  Delta.neg[, 2]), c(x.ord[row.ind.plus], x.ord[row.ind.neg] - 
                  1), mean))
            }
            else {
                Delta <- as.numeric(tapply(c(Delta.plus, Delta.neg), 
                  c(x.ord[row.ind.plus], x.ord[row.ind.neg] - 
                    1), mean))
            }
            fJ <- c(0, cumsum(Delta))
            fJ = fJ - sum(fJ * x.prob[ind.ord])
            x <- levs.ord
            q <- rep(NA, length(x))
            class <- rep("factor", length(x))
        }
        else {
            if (type == "prob") {
                Delta <- sapply(1:ncol(Delta.plus), function(x) as.numeric(tapply(c(Delta.plus[, 
                  x], Delta.neg[, x]), c(x.ord[row.ind.plus], 
                  x.ord[row.ind.neg] - 1), mean)))
                if(is.null(nrow(Delta))){
                    fJ <- rbind(0, cumsum(Delta))
                    fJ = fJ - sum(fJ * x.prob[ind.ord])
                }else{
                    fJ = apply(Delta, 2, function(x) c(0, cumsum(x)))
                    fJ = apply(fJ, 2, function(x) x - sum(x * 
                      x.prob[ind.ord]))
                }                
                x <- levs.ord
                q <- rep(NA, length(x))
                class <- rep("factor", length(x))
            }
            else {
                Delta <- as.numeric(tapply(rbind(Delta.plus, Delta.neg), 
                  c(x.ord[row.ind.plus], x.ord[row.ind.neg] - 
                    1), mean))
                fJ <- c(0, cumsum(Delta))
                fJ = fJ - sum(fJ * x.prob[ind.ord])
                x <- levs.ord
                q <- rep(NA, length(x))
                class <- rep("factor", length(x))
            }
        }
    }
    else if (is(X[, J], "numeric") | is(X[, J], "integer")) {
        z = c(min(X[, J]), as.numeric(quantile(X[, J], seq(1/K, 
            1, length.out = K), type = 1)))
        f <- ecdf(X[, J])
        z = unique(z)
        q <- f(z)
        K = length(z) - 1
        fJ = numeric(K)
        a1 = as.numeric(cut(X[, J], breaks = z, include.lowest = TRUE))
        X1 = X
        X2 = X
        X1[, J] = z[a1]
        X2[, J] = z[a1 + 1]
        y.hat1 = pred.fun(X.model = X.model, newdata = X1, type = type)
        y.hat2 = pred.fun(X.model = X.model, newdata = X2, type = type)
        Delta = y.hat2 - y.hat1
        if (!multi) {
            if (type == "prob") {
                Delta <- as.numeric(tapply(c(Delta[, 2]), a1, 
                  mean))
            }
            else {
                Delta <- as.numeric(tapply(Delta, a1, mean))
            }
            fJ = c(0, cumsum(Delta))
            b1 <- as.numeric(table(a1))
            fJ = fJ - sum((fJ[1:K] + fJ[2:(K + 1)])/2 * b1)/sum(b1)
            x <- z
            class <- rep("numeric", length(x))
        }
        else {
            if (type == "prob") {
                Delta <- sapply(1:ncol(Delta), function(x) as.numeric(tapply(c(Delta[, 
                  x]), a1, mean)))
                fJ = apply(Delta, 2, function(x) c(0, cumsum(x)))
                b1 <- as.numeric(table(a1))
                fJ = apply(fJ, 2, function(fJ) fJ - sum((fJ[1:K] + 
                  fJ[2:(K + 1)])/2 * b1)/sum(b1))
                x <- z
                class <- rep("numeric", length(x))
            }
            else {
                Delta <- as.numeric(tapply(Delta, a1, mean))
                fJ = c(0, cumsum(Delta))
                b1 <- as.numeric(table(a1))
                fJ = fJ - sum((fJ[1:K] + fJ[2:(K + 1)])/2 * b1)/sum(b1)
                x <- z
                class <- rep("numeric", length(x))
            }
        }
    }
    list(K = K, x.values = x, class = class, quantile = q, f.values = fJ)
}