robfmtest <-
function (formula, power = 2:3, type = c("regressor"), 
            data = list(), x.weights = c("HAT", "MCD"), testtype = "Wald", ...) 
  {
    dname <- paste(deparse(substitute(formula)))
    if (!inherits(formula, "formula")) {
      mf <- model.frame(formula)
      X <- if (is.matrix(formula$x)) 
        formula$x
      else model.matrix(terms(formula), mf)
      y <- if (is.vector(formula$y)) 
        formula$y
      else model.response(mf)
    }
    else {
      mf <- model.frame(formula, data = data)
      y <- model.response(mf)
      X <- model.matrix(formula, data = data)
    }
    y <- scale(y, center = attr(terms(mf), "intercept") > 
                 0L)
    k <- ncol(X)
    n <- nrow(X)
    type <- match.arg(type)
    switch(type, regressor = {
      Z <- as.matrix(mf[, which(!sapply(mf, is.factor))[-1]])
      Z <- matrix(as.vector(t(sapply(as.vector(Z), "^", 
                                     power))), nrow = n)
    })
    x.weights <- match.arg(x.weights)
    if (testtype == "Wald") {
      XZ <- cbind(X[, -1], Z)
      q <- ncol(Z)
      switch(x.weights, HAT =  {xwght <- sqrt(1 - hat(X[,-1]))},
       MCD = {
        tmc <- qchisq(0.95, k)
        mcd.est <- cov.rob(X[, -1], method="mcd")
        dist <- mahalanobis(as.matrix(X[, -1]), mcd.est$center, mcd.est$cov) 
        xwght <- ifelse(dist < tmc, 1, tmc/dist)
      })
      rob.fit <- rlm(y ~ XZ, method="M", weights = xwght)
      gmmnum <- seq(k+1, k+q)
      gmr <- rob.fit$coefficients[gmmnum]
      VVr <- vcov(rob.fit)[gmmnum, gmmnum]
      reset <- gmr%*%solve(VVr)%*%gmr
      df <- q
      pval <- as.vector(pchisq(reset, df, lower.tail = FALSE))
    }
    RVAL <- list(statistic = drop(reset), dof = df, method = "FMrob test", 
                 p.value = pval, data.name = dname)
    class(RVAL) <- "robfmtest"
    return(RVAL)
  }
