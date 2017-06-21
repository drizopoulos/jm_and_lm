##########################################################################################
# Aim: Code to run the cross-validation analysis of the Aortic Valve dataset (Table 1)   #
#                                                                                        #
# Required packages: development version of JMbayes from GitHub                          #
#     (https://github.com/drizopoulos/JMbayes), and parallel and xtable from CRAN        #
#                                                                                        #
# Note: The code requires considerable amount of time to run (depending also on your     #
# machine's configuration)                                                               #
#                                                                                        #
# Author: Dimitris Rizopoulos                                                            #
##########################################################################################


library("parallel")
library("xtable")
n <- 500 # number of subjects
V <- 5 # number of folds in the cross-validation
M <- 20 # number of times to replicate the cross-validation procedure
Res <- vector("list", M)
for (mm in seq_len(M)) {
    
    set.seed(123 + mm)
    splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
    
    runCV <- function (i) {
        # load data
        con <- url("https://raw.github.com/drizopoulos/jm_and_lm/master/case_study/simulated_AoValv.RData")
        load(con)
        close(con)
        train_data <- AoValv[!AoValv$id %in% i, ]
        train_data$id <- match(train_data$id, unique(train_data$id))
        test_data <- AoValv[AoValv$id %in% i, ]
        test_data$id <- match(test_data$id, unique(test_data$id))
        #################################################################################
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
        CV_AoValv <- function (train_data, test_data, followup_times = c(5.5, 7.5, 9.5)) {
            test <- try({
                ##################
                # Joint Modeling #
                ##################
                
                library("JMbayes")
                library("splines")
                # Fit Joint Models
                lmeFit <- lme(sqrt(AoGradient) ~ 0 + TypeOp + TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                              data = train_data, 
                              random = list(id = pdDiag(form = ~ ns(time, k = c(2.5, 6), B = c(0.5, 13)))))
                survFit <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex, 
                                 data = train_data[!duplicated(train_data$id), ], x = TRUE)
                # current value
                jointFit1 <- jointModelBayes(lmeFit, survFit, timeVar = "time", 
                                             n.iter = 100000L, 
                                             priors = list(priorA.tauBs = 1e-01, priorB.tauBs = 1e-01))
                
                # current value + slope
                dForm <- list(fixed = ~ 0 + TypeOp:dns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                              indFixed = 3:8, 
                              random = ~ 0 + dns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                              indRandom = 2:4)
                jointFit2 <- update(jointFit1, param = "td-both", extraForm = dForm)
                # cumulative
                iForm <- list(fixed = ~ 0 + TypeOp:time + TypeOp:ins(time, k = c(2.5, 6), B = c(0.5, 13)),
                              indFixed = 1:8, 
                              random = ~ 0 + time + ins(time, k = c(2.5, 6), B = c(0.5, 13)), 
                              indRandom = 1:4)
                jointFit3 <- update(jointFit1, param = "td-extra", extraForm = iForm)
                
                ######################################################################################
                
                # Calculate Performance Measures
                JM_Models <- list("JM_value" = jointFit1, "JM_value-slope" = jointFit2, 
                                  "JM_cumulative" = jointFit3)
                combos <- expand.grid("model_name" = names(JM_Models), "time" = followup_times)
                # AUC
                auc_fun <- function (model_name, time) {
                    aucJM(JM_Models[[model_name]], newdata = test_data, Tstart = time, Dt = 2)$auc
                }
                JM_AUCs <- mapply(auc_fun, combos$model_name, combos$time)
                # PE
                pe_fun <- function (model_name, time) {
                    prederrJM(JM_Models[[model_name]], newdata = test_data, Tstart = time, 
                              Thoriz = time + 2, lossFun = "square")$prederr
                }
                JM_PEs <- mapply(pe_fun, combos$model_name, combos$time)
                
                ######################################################################################
                ######################################################################################
                
                ###############
                # Landmarking #
                ###############
                
                dataLM <- JMbayes:::dataLM
                
                # function to fit Landmark models
                LM_models_fun <- function (time) {
                    # Landmark data sets
                    D1 <- dataLM(train_data, time, respVar = "AoGradient", timeVar = "time", 
                                 evTimeVar = "EvTime")
                    D2 <- dataLM(train_data, time, respVar = "AoGradient", timeVar = "time", 
                                 evTimeVar = "EvTime", summary = "slope", tranfFun = sqrt)
                    D3 <- dataLM(train_data, time, respVar = "AoGradient", timeVar = "time", 
                                 evTimeVar = "EvTime", summary = "area", tranfFun = sqrt)
                    # Cox models
                    CoxLM1 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + sqrt(AoGradient), 
                                    data = D1)
                    CoxLM2 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + sqrt(AoGradient) + 
                                        slope, data = D2)
                    CoxLM3 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + area, data = D3)
                    list("LM_value" = CoxLM1, "LM_value-slope" = CoxLM2, "LM_cumulative" = CoxLM3)
                }
                
                ######################################################################################
                
                # Calculate Performance Measures
                # AUC
                auc_fun <- function (time) {
                    LM_models <- LM_models_fun(time)
                    auc_objs <- mapply(aucJM, LM_models, summary = c("value", "slope", "area"),
                                       MoreArgs = list(newdata = test_data, Tstart = time, 
                                                       Dt = 2, timeVar = "time", 
                                                       respVar = "AoGradient", 
                                                       evTimeVar = "EvTime", 
                                                       tranfFun = sqrt), SIMPLIFY = FALSE)
                    sapply(auc_objs, "[[", "auc")
                }
                LM_AUCs <- sapply(followup_times, auc_fun)
                # PE
                pe_fun <- function (time) {
                    LM_models <- LM_models_fun(time)
                    pe_objs <- mapply(prederrJM, LM_models, summary = c("value", "slope", "area"), 
                                      MoreArgs = list(newdata = test_data, Tstart = time, 
                                                      Thoriz = time + 2, timeVar = "time", 
                                                      respVar = "AoGradient", 
                                                      evTimeVar = "EvTime", tranfFun = sqrt, 
                                                      lossFun = "square"), SIMPLIFY = FALSE)
                    sapply(pe_objs, "[[", "prederr")
                }
                LM_PEs <- sapply(followup_times, pe_fun)
                
                ######################################################################################
                ######################################################################################
                
                ##############################
                # Landmarking + Mixed Models #
                ##############################
                
                # function to fit Landmark mixed models
                LMmixed_models_fun <- function (time) {
                    # create LM data
                    active_train <- train_data[train_data$EvTime > time & train_data$time <= time, ]
                    active_test <- test_data[test_data$EvTime > time & test_data$time <= time, ]
                    # fit mixed model
                    lmeFit <- lme(sqrt(AoGradient) ~ 0 + TypeOp + TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                                  data = active_train, 
                                  random = list(id = pdDiag(form = ~ ns(time, k = c(2.5, 6), B = c(0.5, 13)))))
                    betas <- fixef(lmeFit)
                    b <- data.matrix(ranef(lmeFit))
                    b_test <- calculate_REs(lmeFit, active_test)
                    id <- match(active_train$id, unique(active_train$id))
                    id_test <- match(active_test$id, unique(active_test$id))
                    active_train <- active_train[tapply(row.names(active_train), id, tail, 1), ]
                    active_test <- active_test[tapply(row.names(active_test), id_test, tail, 1), ]
                    active_train[["time"]] <- time
                    active_test[["time"]] <- time
                    # calculate fitted values at event item
                    Xvalue <- model.matrix(sqrt(AoGradient) ~ 0 + TypeOp + 
                                               TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                                           data = active_train)
                    Zvalue <- model.matrix(~ ns(time, k = c(2.5, 6), B = c(0.5, 13)), data = active_train)
                    active_train$value <- c(Xvalue %*% betas) + rowSums(Zvalue * b)
                    ##
                    Xvalue_test <- model.matrix(sqrt(AoGradient) ~ 0 + TypeOp + 
                                                    TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                                                data = active_test)
                    Zvalue_test <- model.matrix(~ ns(time, k = c(2.5, 6), B = c(0.5, 13)), data = active_test)
                    active_test$value <- c(Xvalue_test %*% betas) + rowSums(Zvalue_test * b_test)
                    ##
                    Xslope <- model.matrix(dForm$fixed, data = active_train)
                    Zslope <- model.matrix(dForm$random, data = active_train)
                    active_train$slope <- c(Xslope %*% betas[dForm$indFixed]) + 
                        rowSums(Zslope * b[, dForm$indRandom])
                    ##
                    Xslope_test <- model.matrix(dForm$fixed, data = active_test)
                    Zslope_test <- model.matrix(dForm$random, data = active_test)
                    active_test$slope <- c(Xslope_test %*% betas[dForm$indFixed]) + 
                        rowSums(Zslope_test * b_test[, dForm$indRandom])
                    ##
                    Xarea <- model.matrix(iForm$fixed, data = active_train)
                    Zarea <- model.matrix(iForm$random, data = active_train)
                    active_train$area <- c(Xarea %*% betas[iForm$indFixed]) + 
                        rowSums(Zarea * b[, iForm$indRandom])
                    ##
                    Xarea_test <- model.matrix(iForm$fixed, data = active_test)
                    Zarea_test <- model.matrix(iForm$random, data = active_test)
                    active_test$area <- c(Xarea_test %*% betas[iForm$indFixed]) + 
                        rowSums(Zarea_test * b_test[, iForm$indRandom])
                    # Fit Cox models
                    CoxLM1 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + 
                                        value, data = active_train)
                    CoxLM2 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + 
                                        value + slope, 
                                    data = active_train)
                    CoxLM3 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + 
                                        area, data = active_train)
                    list(models = list("LMmixed_value" = CoxLM1, "LMmixed_value-slope" = CoxLM2, 
                                       "LMmixed_cumulative" = CoxLM3), data = active_test)
                }
                
                ######################################################################################
                
                # Calculate Performance Measures
                # AUC
                auc_fun <- function (time) {
                    run_LMmixed_models <- LMmixed_models_fun(time)
                    LM_models <- run_LMmixed_models$models
                    aucs <- numeric(length(LM_models))
                    for (i in seq_along(aucs)) {
                        test <- names(LM_models[i]) %in% c("LMmixed_value", "LMmixed_value-slope")
                        aucs[i] <- aucJM(LM_models[[i]], newdata = run_LMmixed_models$data, 
                                         Tstart = time, Dt = 2, timeVar = "time", evTimeVar = "EvTime",
                                         respVar = if (test) "value" else "area")$auc
                    }
                    aucs
                }
                LMmixed_AUCs <- sapply(followup_times, auc_fun)
                
                # PE
                pe_fun <- function (time) {
                    run_LMmixed_models <- LMmixed_models_fun(time)
                    LM_models <- run_LMmixed_models$models
                    pes <- numeric(length(LM_models))
                    for (i in seq_along(pes)) {
                        test <- names(LM_models[i]) %in% c("LMmixed_value", "LMmixed_value-slope")
                        pes[i] <- prederrJM(LM_models[[i]], newdata = run_LMmixed_models$data, 
                                            Tstart = time, Thoriz = time + 2, timeVar = "time", 
                                            evTimeVar = "EvTime", lossFun = "square", 
                                            respVar = if (test) "value" else "area")$prederr
                    }
                    pes
                }
                LMmixed_PEs <- sapply(followup_times, pe_fun)
                
                ######################################################################################
                ######################################################################################
                
                # Collect results
                
                AoValv_accuracy <- expand.grid("functional form" = c("value", "value-slope", "cumulative"), 
                                               "time" = followup_times)
                AoValv_accuracy$JM_AUC <- JM_AUCs
                AoValv_accuracy$LM_AUC <- c(LM_AUCs)
                AoValv_accuracy$LMmixed_AUC <- c(LMmixed_AUCs)
                AoValv_accuracy$JM_PE <- JM_PEs
                AoValv_accuracy$LM_PE <- c(LM_PEs)
                AoValv_accuracy$LMmixed_PE <- c(LMmixed_PEs)
                AoValv_accuracy
            }, TRUE)
            if (!inherits(test, "try-error")) {
                return(test)
            } else {
                AoValv_accuracy <- expand.grid("functional form" = c("value", "value-slope", "cumulative"), 
                                               "time" = followup_times)
                AoValv_accuracy$JM_AUC <- AoValv_accuracy$LM_AUC <- NA
                AoValv_accuracy$LMmixed_AUC <- AoValv_accuracy$JM_PE <- NA
                AoValv_accuracy$LM_PE <- AoValv_accuracy$LMmixed_PE <- NA
                return(AoValv_accuracy)
            }
        }
        CV_AoValv(train_data, test_data)
    }
    cl <- makeCluster(5)
    time <- system.time(res <- parLapply(cl, splits, runCV))
    stopCluster(cl)
    Res[[mm]] <- res
    print(mm)
}

##########################################################################################
##########################################################################################

out <- do.call("rbind", unlist(Res, recursive = FALSE))
out <- data.frame(
    "functional_form" = rep(out$`functional form`, 3),
    "time" = factor(rep(out$time, 3), labels = paste("t* =", c(5.5, 7.5, 9.5))),
    "AUC" = c(out$JM_AUC, out$LM_AUC, out$LMmixed_AUC),
    "PE" = c(out$JM_PE, out$LM_PE, out$LMmixed_PE),
    "model" = rep(c("JM", "LM", "LM\nmixed"), each = nrow(out))
)
out$functional_form <- factor(out$functional_form, 
                              labels = c("value", "value + slope", "cumulative"))

AUCs <- round(with(out, tapply(AUC, list(time, model, functional_form), mean, na.rm = TRUE)), 3)
PEs <- round(with(out, tapply(PE, list(time, model, functional_form), mean, na.rm = TRUE)), 3)

followup_times <- c(5.5, 7.5, 9.5)

AoValv_accuracy <- expand.grid("functional form" = c("value", "value+slope", "cumulative"), 
                               "time" = followup_times)

AoValv_accuracy$JM_AUC <- c(t(AUCs[, 1,]))
AoValv_accuracy$LM_AUC <- c(t(AUCs[, 2,]))
AoValv_accuracy$LMmixed_AUC <- c(t(AUCs[, 3,]))
AoValv_accuracy$JM_PE <- c(t(PEs[, 1,]))
AoValv_accuracy$LM_PE <- c(t(PEs[, 2,]))
AoValv_accuracy$LMmixed_PE <- c(t(PEs[, 3,]))

# Table 1 main paper

print(xtable(AoValv_accuracy, label = "Tab:AccMeas", caption = NULL, 
             align = c("l", rep("r", 8)), digits = c(0, 0, 1, rep(3, 6))), 
      math.style.negative = TRUE, sanitize.text.function = function (x) x,
      include.rownames = FALSE)
