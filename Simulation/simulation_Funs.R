calculate_REs <- function (lmeObject, newdata) {
    data <- lmeObject$data
    formYx <- formula(lmeObject)
    mfX <- model.frame(terms(formYx), data = data)
    TermsX <- attr(mfX, "terms")
    mfX_new <- model.frame(TermsX, data = newdata)
    X_new <- model.matrix(formYx, mfX_new)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    mfZ_new <- model.frame(TermsZ, data = newdata)
    Z_new <- model.matrix(formYz, mfZ_new)
    y_new <- model.response(mfX_new, "numeric")
    idVar <- names(lmeObject$modelStruct$reStruct)
    if (length(idVar) > 1)
        stop("the current version of the function only works with a single grouping variable.\n")
    if (is.null(newdata[[idVar]]))
        stop("subject id variable not in newdata.")
    id <- match(newdata[[idVar]], unique(newdata[[idVar]]))
    n <- length(unique(id))
    betas <- fixef(lmeObject)
    D <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*",
                lmeObject$sigma^2)[[1]]
    modes <- matrix(0.0, n, ncol(Z_new))
    for (i in seq_len(n)) {
        id_i <- id == i
        X_new_id <- X_new[id_i, , drop = FALSE]
        Z_new_id <- Z_new[id_i, , drop = FALSE]
        Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) + 
                            lmeObject$sigma^2 * diag(sum(id_i)))
        DZtVinv <- tcrossprod(D, Z_new_id) %*% Vi_inv
        modes[i, ] <- c(DZtVinv %*% (y_new[id_i] - X_new_id %*% betas))
    }
    modes
}

true_pi <- function (u, t, trueValues, test_data, keep, scenario) {
    betas <- trueValues[["betas"]]
    gammas <- trueValues[["gammas"]]
    if (scenario %in% c("I", "III")) {
        alpha <- trueValues[["alpha"]]
    } else {
        alpha1 <- trueValues[["alpha1"]]
        alpha2 <- trueValues[["alpha2"]]
    }
    phi <- trueValues[["phi"]]
    b <- trueValues[["b_test"]][keep, ]
    Bkn <- trueValues[["Bkn"]]
    kn <- trueValues[["kn"]]
    W <- model.matrix(~ group, data = test_data[!duplicated(test_data$id), ])
    group <- W[, 2]
    eta.t <- c(W %*% gammas)
    S <- function (t, i) {
        h <- function (s) {
            group0 <- 1 - group[i]
            group1 <- group[i]
            if (scenario == "I") {
                BS <- ns(s, knots = kn, Boundary.knots = Bkn)
                XX <- cbind(group0, group1, group0*BS[, 1], group1*BS[, 1],
                            group0*BS[, 2], group1*BS[, 2], group0*BS[, 3], group1*BS[, 3])
                ZZ <- cbind(1, BS)
                f <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
                exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f * alpha)
            } else if (scenario == "II") {
                BS <- ns(s, knots = kn, Boundary.knots = Bkn)
                XX <- cbind(group0, group1, group0*BS[, 1], group1*BS[, 1],
                            group0*BS[, 2], group1*BS[, 2], group0*BS[, 3], group1*BS[, 3])
                ZZ <- cbind(1, BS)
                f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
                dBS <- dns(s, knots = kn, Boundary.knots = Bkn)
                XXd <- cbind(group0*dBS[, 1], group1*dBS[, 1],
                             group0*dBS[, 2], group1*dBS[, 2], group0*dBS[, 3], group1*dBS[, 3])
                ZZd <- dBS
                f2 <- as.vector(XXd %*% betas[3:8] + rowSums(ZZd * b[rep(i, nrow(ZZd)), 2:4]))
                exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha1 + f2 * alpha2)
            } else {
                iBS <- ins(s, knots = kn, Boundary.knots = Bkn)
                XXi <- cbind(group0*s, group1*s, group0*iBS[, 1], group1*iBS[, 1],
                             group0*iBS[, 2], group1*iBS[, 2], group0*iBS[, 3], group1*iBS[, 3])
                ZZi <- cbind(s, iBS)
                f <- as.vector(XXi %*% betas + rowSums(ZZi * b[rep(i, nrow(ZZi)), ]))
                exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f * alpha)
            }
        }
        exp(-integrate(h, lower = 0, upper = t)$value)
    }
    n <- length(eta.t)
    if (length(u) < n)
        u <- rep(u, length = n)
    if (length(t) < n)
        t <- rep(t, length = n)
    out <- numeric(n)
    for (i in seq_len(n)) {
        out[i] <- S(u[i], i) / S(t[i], i) 
    }
    names(out) <- names(t)
    out
}

aucJM2 <- function (object, newdata, Tstart, ...) {
    UseMethod("aucJM2")
}

aucJM2.JMbayes <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, 
                            idVar = "id", simulate = FALSE, M = 100, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata'.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    Thoriz <- Thoriz + 1e-07
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$Terms$termsT
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    is_counting <- attr(SurvT, "type") == "counting"
    Time <- if (is_counting) {
        ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
    } else {
        SurvT[, 1]
    }
    timeVar <- object$timeVar
    ordTime <- order(Time)
    newdata2 <- newdata[ordTime, ]
    newdata2 <- newdata2[Time[ordTime] > Tstart, ]
    newdata2 <- newdata2[newdata2[[timeVar]] <= Tstart, ]
    pi.u.t <- if (is_counting) {
        survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, 
                  simulate = simulate, M = M, LeftTrunc_var = all.vars(TermsT)[1L])
    } else {
        survfitJM(object, newdata = newdata2, idVar = idVar, survTimes = Thoriz, 
                  simulate = simulate, M = M)
    }
    pi.u.t <- sapply(pi.u.t$summaries, "[", 1, 2)
    # find comparable subjects
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    if (is_counting) {
        f <- factor(id, levels = unique(id))
        Time <- tapply(SurvT[, 2], f, tail, 1)
        event <- tapply(SurvT[, 3], f, tail, 1)
    } else{
        Time <- SurvT[!duplicated(id), 1]
        event <- SurvT[!duplicated(id), 2]
    }
    names(Time) <- names(event) <- as.character(unique(id))
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + 1e-07
    }
    if (!all(names(pi.u.t) == names(Time)))
        stop("mismatch between 'Time' variable names and survival probabilities names.")
    auc <- if (length(Time) > 1) {
        pairs <- combn(as.character(unique(id)), 2)
        Ti <- Time[pairs[1, ]]
        Tj <- Time[pairs[2, ]]
        di <- event[pairs[1, ]]
        dj <- event[pairs[2, ]]
        pi.u.t.i <- pi.u.t[pairs[1, ]]
        pi.u.t.j <- pi.u.t[pairs[2, ]]
        ind1 <- (Ti <= Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti <= Thoriz & di == 0) & Tj > Thoriz
        ind3 <- (Ti <= Thoriz & di == 1) & (Tj <= Thoriz & dj == 0)
        ind4 <- (Ti <= Thoriz & di == 0) & (Tj <= Thoriz & dj == 0)
        names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- paste(names(Ti), 
                                                                          names(Tj), 
                                                                          sep = "_")
        ind <- ind1 | ind2 | ind3 | ind4
        if (any(ind2)) {
            nams <- strsplit(names(ind2[ind2]), "_")
            nams_i <- sapply(nams, "[", 1)
            unq_nams_i <- unique(nams_i)
            newdata2_i <- newdata2[id %in% unq_nams_i, ]
            keep_i <- unique(newdata2_i[[idVar]])
            pi2 <- 1 - true_pi(Thoriz, Time[unq_nams_i], trueValues, newdata2_i, keep_i, 
                               scenario)
            ind[ind2] <- ind[ind2] * pi2[nams_i]
        }
        if (any(ind3)) {
            nams <- strsplit(names(ind3[ind3]), "_")
            nams_j <- sapply(nams, "[", 2)
            unq_nams_j <- unique(nams_j)
            newdata2_j <- newdata2[id %in% unq_nams_j, ]
            keep_j <- unique(newdata2_j[[idVar]])
            pi3 <- true_pi(Thoriz, Time[unq_nams_j], trueValues, newdata2_j, keep_j, 
                           scenario)
            ind[ind3] <- ind[ind3] * pi3[nams_j]
        }
        if (any(ind4)) {
            nams <- strsplit(names(ind4[ind4]), "_")
            nams_i <- sapply(nams, "[", 1)
            nams_j <- sapply(nams, "[", 2)
            unq_nams_i <- unique(nams_i)
            unq_nams_j <- unique(nams_j)
            newdata2_i <- newdata2[id %in% unq_nams_i, ]
            keep_i <- unique(newdata2_i[[idVar]])
            newdata2_j <- newdata2[id %in% unq_nams_j, ]
            keep_j <- unique(newdata2_j[[idVar]])
            pi4_i <- 1 - true_pi(Thoriz, Time[unq_nams_i], trueValues, newdata2_i, keep_i, 
                                 scenario)
            pi4_j <- true_pi(Thoriz, Time[unq_nams_j], trueValues, newdata2_j, keep_j, 
                             scenario)
            
            ind[ind4] <- ind[ind4] * pi4_i[nams_i] * pi4_j[nams_j]
        }
        sum((pi.u.t.i < pi.u.t.j) * c(ind), na.rm = TRUE) / sum(ind, na.rm = TRUE)
    } else {
        NA
    }
    out <- list(auc = auc, Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)), 
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "aucJM"
    out
}

aucJM2.coxph <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, idVar = "id", 
                          respVar = "y", timeVar = "time", evTimeVar = "Time",
                          summary = c("value", "slope", "area"), 
                          tranfFun = function (x) x, ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    Thoriz <- Thoriz + 1e-07
    summary <- match.arg(summary)
    if (summary %in% c("slope", "area"))
        newdata$area <- newdata$slope <- 0
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$terms
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    ordTime <- order(Time)
    newdata2 <- newdata[ordTime, ]
    newdata2 <- dataLM(newdata2, Tstart, idVar, respVar, timeVar, evTimeVar, summary, 
                       tranfFun)
    pi.u.t <- c(summary(survfit(object, newdata = newdata2), times = Thoriz)$surv)
    # find comparable subjects
    id <- newdata2[[idVar]]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    Time <- SurvT[!duplicated(id), 1]
    event <- SurvT[!duplicated(id), 2]
    names(pi.u.t) <- names(Time) <- names(event) <- as.character(unique(id))
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + 1e-07
    }
    auc <- if (length(Time) > 1) {
        pairs <- combn(as.character(unique(id)), 2)
        Ti <- Time[pairs[1, ]]
        Tj <- Time[pairs[2, ]]
        di <- event[pairs[1, ]]
        dj <- event[pairs[2, ]]
        pi.u.t.i <- pi.u.t[pairs[1, ]]
        pi.u.t.j <- pi.u.t[pairs[2, ]]
        ind1 <- (Ti <= Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti <= Thoriz & di == 0) & Tj > Thoriz
        ind3 <- (Ti <= Thoriz & di == 1) & (Tj <= Thoriz & dj == 0)
        ind4 <- (Ti <= Thoriz & di == 0) & (Tj <= Thoriz & dj == 0)
        names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <- paste(names(Ti), 
                                                                          names(Tj), 
                                                                          sep = "_")
        ind <- ind1 | ind2 | ind3 | ind4
        if (any(ind2)) {
            nams <- strsplit(names(ind2[ind2]), "_")
            nams_i <- sapply(nams, "[", 1)
            unq_nams_i <- unique(nams_i)
            ND <- newdata2[id %in% unq_nams_i, ]
            tt <- model.response(model.frame(TermsT, ND))[, 1]
            keep_i <- unique(ND[[idVar]])
            pi2 <- 1 - true_pi(Thoriz, tt, trueValues, ND, keep_i, scenario)
            ind[ind2] <- ind[ind2] * pi2[nams_i]
        }
        if (any(ind3)) {
            nams <- strsplit(names(ind3[ind3]), "_")
            nams_j <- sapply(nams, "[", 2)
            unq_nams_j <- unique(nams_j)
            ND <- newdata2[id %in% unq_nams_j, ]
            tt <- model.response(model.frame(TermsT, ND))[, 1]
            keep_j <- unique(ND[[idVar]])
            pi3 <- true_pi(Thoriz, tt, trueValues, ND, keep_j, scenario)
            ind[ind3] <- ind[ind3] * pi3[nams_j]
        }
        if (any(ind4)) {
            nams <- strsplit(names(ind4[ind4]), "_")
            nams_i <- sapply(nams, "[", 1)
            nams_j <- sapply(nams, "[", 2)
            unq_nams_i <- unique(nams_i)
            unq_nams_j <- unique(nams_j)
            ND_i <- newdata2[id %in% unq_nams_i, ]
            ND_j <- newdata2[id %in% unq_nams_j, ]
            tt_i <- model.response(model.frame(TermsT, ND_i))[, 1]
            tt_j <- model.response(model.frame(TermsT, ND_j))[, 1]
            keep_i <- unique(ND_i[[idVar]])
            keep_j <- unique(ND_j[[idVar]])
            pi4_i <- 1 - true_pi(Thoriz, tt_i, trueValues, ND_i, keep_i, scenario)
            pi4_j <- true_pi(Thoriz, tt_j, trueValues, ND_j, keep_j, scenario)
            ind[ind4] <- ind[ind4] * pi4_i[nams_i] * pi4_j[nams_j]
        }
        sum((pi.u.t.i < pi.u.t.j) * c(ind), na.rm = TRUE) / sum(ind, na.rm = TRUE)
    } else {
        NA
    }    
    out <- list(auc = auc, Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)), 
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "aucJM"
    out
}

prederrJM2 <- function (object, newdata, Tstart, Thoriz, ...) {
    UseMethod("prederrJM2")
}

prederrJM2.JMbayes <- function (object, newdata, Tstart, Thoriz, 
                                lossFun = c("square", "absolute"), interval = FALSE, 
                                idVar = "id", simulate = FALSE, M = 100, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    lossFun <- if (is.function(lossFun)) {
        lf <- lossFun
        match.fun(lossFun)
    } else {
        lf <- match.arg(lossFun)
        if (lf == "absolute") function (x) abs(x) else function (x) x*x
    }
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$Terms$termsT
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    is_counting <- attr(SurvT, "type") == "counting"
    Time <- if (is_counting) {
        ave(SurvT[, 2], id, FUN = function (x) tail(x, 1))
    } else {
        SurvT[, 1]
    }
    timeVar <- object$timeVar
    newdata2 <- newdata[Time > Tstart, ]
    SurvT <- model.response(model.frame(TermsT, newdata2))
    if (is_counting) {
        id2 <- newdata2[[idVar]]
        f <- factor(id2, levels = unique(id2))
        Time <- ave(SurvT[, 2], f, FUN = function (x) tail(x, 1))
        delta <- ave(SurvT[, 3], f, FUN = function (x) tail(x, 1))
    } else {
        Time <- SurvT[, 1]
        delta <- SurvT[, 2]
    }
    timesInd <- newdata2[[timeVar]] <= Tstart
    aliveThoriz <- newdata2[Time > Thoriz & timesInd, ]
    deadThoriz <- newdata2[Time <= Thoriz & delta == 1 & timesInd, ]
    indCens <- Time < Thoriz & delta == 0 & timesInd
    censThoriz <- newdata2[indCens, ]
    nr <- length(unique(newdata2[[idVar]]))
    idalive <- unique(aliveThoriz[[idVar]])
    iddead <- unique(deadThoriz[[idVar]])
    idcens <- unique(censThoriz[[idVar]])
    prederr <- if (length(unique(Time)) > 1 && nrow(aliveThoriz) > 1 &&
                   nrow(deadThoriz) > 1) {
        Surv.aliveThoriz <- if (is_counting) {
            survfitJM(object, newdata = aliveThoriz, idVar = idVar, simulate = simulate, 
                      M = M, survTimes = Thoriz, last.time = rep(Tstart, length(idalive)),
                      LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
            survfitJM(object, newdata = aliveThoriz, idVar = idVar, simulate = simulate, 
                      M = M, survTimes = Thoriz, last.time = rep(Tstart, length(idalive)))
        }
        Surv.deadThoriz <- if (is_counting) {
            survfitJM(object, newdata = deadThoriz, idVar = idVar, simulate = simulate,
                      survTimes = Thoriz, last.time = rep(Tstart, length(iddead)),
                      LeftTrunc_var = all.vars(TermsT)[1L])
        } else {
            survfitJM(object, newdata = deadThoriz, idVar = idVar, simulate = simulate,
                      survTimes = Thoriz, last.time = rep(Tstart, length(iddead)))
        }
        Surv.aliveThoriz <- sapply(Surv.aliveThoriz$summaries, "[", 2)
        Surv.deadThoriz <- sapply(Surv.deadThoriz$summaries, "[", 2)
        if (nrow(censThoriz)) {
            Surv.censThoriz <- if (is_counting) {
                survfitJM(object, newdata = censThoriz, idVar = idVar, simulate = simulate, 
                          M = M, survTimes = Thoriz, last.time = rep(Tstart, length(idcens)),
                          LeftTrunc_var = all.vars(TermsT)[1L])
            } else {
                survfitJM(object, newdata = censThoriz, idVar = idVar, simulate = simulate, 
                          M = M, survTimes = Thoriz, last.time = rep(Tstart, length(idcens)))
            }
            Surv.censThoriz <- sapply(Surv.censThoriz$summaries, "[", 2)
            tt <- Time[indCens]
            keep <- unique(censThoriz[[idVar]])
            weights <- true_pi(Thoriz, tt[!duplicated(censThoriz[[idVar]])], 
                                trueValues, censThoriz, keep, scenario)
        } else {
            Surv.censThoriz <- weights <- NA
        }
        if (!interval) {
            (1/nr) * sum(lossFun(1 - Surv.aliveThoriz), lossFun(0 - Surv.deadThoriz),
                         weights * lossFun(1 - Surv.censThoriz) + 
                             (1 - weights) * lossFun(0 - Surv.censThoriz), na.rm = TRUE)
        } else {
            TimeCens <- object$y$Time
            deltaCens <- 1 - object$y$event
            KMcens <- survfit(Surv(TimeCens, deltaCens) ~ 1)
            times <- TimeCens[TimeCens > Tstart & TimeCens < Thoriz & !deltaCens]
            times <- sort(unique(times))
            k <- as.numeric(table(times))
            w <- summary(KMcens, times = Tstart)$surv / summary(KMcens, times = times)$surv
            prederr.times <- sapply(times, 
                                    function (t) prederrJM(object, newdata, Tstart, t,
                                                           interval = FALSE, idVar = idVar, 
                                                           simulate = simulate)$prederr)
            num <- sum(prederr.times * w * k, na.rm = TRUE)
            den <- sum(w * k, na.rm = TRUE)
            num / den
        }
    } else {
        nr <- NA
        NA
    }
    out <- list(prederr = prederr, nr = nr, Tstart = Tstart, Thoriz = Thoriz, 
                interval = interval, classObject = class(object), 
                nameObject = deparse(substitute(object)), lossFun = lf)
    class(out) <- "prederrJM"
    out
}

prederrJM2.coxph <- function (object, newdata, Tstart, Thoriz, 
                              lossFun = c("absolute", "square"), interval = FALSE, 
                              idVar = "id", timeVar = "time", respVar = "y", 
                              evTimeVar = "Time", summary = c("value", "slope", "area"), 
                              tranfFun = function (x) x, ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata'.\n")
    lossFun <- if (is.function(lossFun)) {
        lf <- lossFun
        match.fun(lossFun)
    } else {
        lf <- match.arg(lossFun)
        if (lf == "absolute") function (x) abs(x) else function (x) x * x
    }
    summary <- match.arg(summary)
    if (summary %in% c("slope", "area"))
        newdata$area <- newdata$slope <- 0
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$terms
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    newdata2 <- dataLM(newdata, Tstart, idVar, respVar, timeVar, evTimeVar, summary, 
                       tranfFun)
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    Time <- SurvT[, 1]
    delta <- SurvT[, 2]
    indCens <- Time < Thoriz & delta == 0
    nr <- nrow(newdata2)
    aliveThoriz.id <- newdata2[Time > Thoriz, ]
    deadThoriz.id <- newdata2[Time <= Thoriz & delta == 1, ]
    prederr <- if (length(unique(Time)) > 1 && nrow(aliveThoriz.id) > 1 &&
                   nrow(deadThoriz.id) > 1) {
        Surv.aliveThoriz <- c(summary(survfit(object, newdata = aliveThoriz.id), 
                                      times = Thoriz)$surv)
        Surv.deadThoriz <- c(summary(survfit(object, newdata = deadThoriz.id), 
                                     times = Thoriz)$surv)
        if (sum(indCens) > 1) {
            censThoriz.id <- newdata2[indCens, ]
            Surv.censThoriz <- c(summary(survfit(object, newdata = censThoriz.id), 
                                         times = Thoriz)$surv)
            tt <- model.response(model.frame(TermsT, censThoriz.id))[, 1]
            nn <- length(tt)
            keep <- unique(censThoriz.id[[idVar]])
            weights <- true_pi(Thoriz, tt, trueValues, censThoriz.id, keep, scenario)
        } else {
            Surv.censThoriz <- weights <- NA
        }
        if (!interval) {
            (1/nr) * sum(lossFun(1 - Surv.aliveThoriz), lossFun(0 - Surv.deadThoriz),
                         weights * lossFun(1 - Surv.censThoriz) + 
                             (1 - weights) * lossFun(0 - Surv.censThoriz))
        } else {
            TimeCens <- model.response(model.frame(TermsT, newdata))[, 1]
            deltaCens <- 1 - model.response(model.frame(TermsT, newdata))[, 2]
            KMcens <- survfit(Surv(TimeCens, deltaCens) ~ 1)
            times <- TimeCens[TimeCens > Tstart & TimeCens <= Thoriz & !deltaCens]
            times <- sort(unique(times))
            k <- as.numeric(table(times))
            w <- summary(KMcens, times = Tstart)$surv / summary(KMcens, times = times)$surv
            prederr.times <- sapply(times, 
                                    function (t) prederrJM(object, newdata, Tstart, t,
                                                           interval = FALSE, idVar = idVar, 
                                                           timeVar = timeVar,
                                                           respVar = respVar, 
                                                           evTimeVar = evTimeVar, 
                                                           summary = summary, 
                                                           tranfFun = tranfFun)$prederr)
            num <- sum(prederr.times * w * k, na.rm = TRUE)
            den <- sum(w * k, na.rm = TRUE)
            num / den
        }
    } else {
        nr <- NA
        NA
    }
    out <- list(prederr = prederr, nr = nr, Tstart = Tstart, Thoriz = Thoriz, 
                interval = interval, classObject = class(object), 
                nameObject = deparse(substitute(object)), lossFun = lf)
    class(out) <- "prederrJM"
    out
}

