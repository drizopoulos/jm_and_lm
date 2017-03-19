##########################################################################################
# Aim: Read simulation results and produce the figures in the supplementary material     #
# Input: The user needs to specify the scenario for which the plots are required; options# 
#        are Ia, Ib, IIa, IIb, IIIa, IIIb. These correspond to Scenarios I, II and II    #
#        described in the paper, and the suffix 'a' or 'b' denote whether the correct of #
#        misspecified functional form for the time effect the linear mixed model were    #
#        used.                                                                           #
# Author: Dimitris Rizopoulos                                                            #
##########################################################################################


library("JMbayes")
library("splines")
library("xtable")

con <- url("https://raw.github.com/drizopoulos/jm_and_lm/master/case_study/simulated_AoValv.RData")
load(con)
close(con)

# Joint Modeling
lmeFit <- lme(sqrt(AoGradient) ~ 0 + TypeOp + TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), 
              data = AoValv, 
              random = list(id = pdDiag(form = ~ ns(time, k = c(2.5, 6), B = c(0.5, 13)))))
survFit <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex, data = AoValv.id, 
                 x = TRUE)

# current value
jointFit1 <- jointModelBayes(lmeFit, survFit, timeVar = "time", n.iter = 100000L)

# current value + slope
dForm <- list(fixed = ~ 0 + TypeOp:dns(time, k = c(2.5, 6), B = c(0.5, 13)), indFixed = 3:8, 
              random = ~ 0 + dns(time, k = c(2.5, 6), B = c(0.5, 13)), indRandom = 2:4)
jointFit2 <- update(jointFit1, param = "td-both", extraForm = dForm)

# cumulative
iForm <- list(fixed = ~ 0 + TypeOp:time + TypeOp:ins(time, k = c(2.5, 6), B = c(0.5, 13)),
              indFixed = 1:8, random = ~ 0 + time + ins(time, k = c(2.5, 6), B = c(0.5, 13)), 
              indRandom = 1:4)
jointFit3 <- update(jointFit1, param = "td-extra", extraForm = iForm)


###################
# Extract Results #
###################

# Longitudinal Process

fCI <- function (x) {
    g <- function (y) {
        y <- sprintf("%.3f", y)
        out <- paste("(", paste(y, collapse = "; "), ")", sep = "")
        gsub("-", "$-$", out, fixed = TRUE)
    }
    if (is.matrix(x)) apply(x, 2, g) else g(x)
}
dY <- data.frame(
    Value = c(jointFit1$postMeans$betas, jointFit1$postMeans$sigma), 
    "95\\% CI" = c(fCI(jointFit1$CI$betas), fCI(jointFit1$CI$sigma)),
    #
    Value = c(jointFit2$postMeans$betas, jointFit2$postMeans$sigma), 
    "95\\% CI" = c(fCI(jointFit2$CI$betas), fCI(jointFit2$CI$sigma)),
    #
    Value = c(jointFit3$postMeans$betas, jointFit3$postMeans$sigma), 
    "95\\% CI" = c(fCI(jointFit3$CI$betas), fCI(jointFit3$CI$sigma)),
    check.names = FALSE)
row.names(dY) <- c(head(gsub("ns(time, k = c(2.5, 6), B = c(0.5, 13))", "B-spln", row.names(dY), fixed = TRUE), -1), "$\\sigma$")

print(xtable(dY, label = "Tab:Res-Y", 
             caption = "Estimated coefficients and 95\\% credibility intervals for the parameters of the longitudinal submodels.",
             align = c("l", rep("r", 6))), 
      math.style.negative = TRUE, sanitize.text.function = function (x) x)

dT <- data.frame(
    Value = c(jointFit1$postMeans$gammas, jointFit1$postMeans$alphas, NA), 
    "95\\% CI" = c(fCI(jointFit1$CI$gammas), fCI(jointFit1$CI$alphas), NA),
    #
    Value = c(jointFit2$postMeans$gammas, jointFit2$postMeans$alphas, jointFit2$postMeans$Dalphas), 
    "95\\% CI" = c(fCI(jointFit2$CI$gammas), fCI(jointFit2$CI$alphas), fCI(jointFit2$CI$Dalphas)),
    #
    Value = c(jointFit3$postMeans$gammas, jointFit3$postMeans$Dalphas, NA), 
    "95\\% CI" = c(fCI(jointFit3$CI$gammas), fCI(jointFit3$CI$Dalphas), NA),
    check.names = FALSE)

print(xtable(dT, label = "Tab:Res-T", 
             caption = "Estimated coefficients and 95\\% credibility intervals for the parameters of the survival submodels.",
             align = c("l", rep("r", 6))), 
      math.style.negative = TRUE, sanitize.text.function = function (x) x)


##########################################################################################
##########################################################################################
##########################################################################################

# Standard Landmarking

dataLM <- JMbayes:::dataLM

LM_models_fun <- function (time) {
    # Landmark data sets
    D1 <- dataLM(AoValv, time, respVar = "AoGradient", timeVar = "time", 
                 evTimeVar = "EvTime")
    D2 <- dataLM(AoValv, time, respVar = "AoGradient", timeVar = "time", 
                 evTimeVar = "EvTime", summary = "slope", tranfFun = sqrt)
    D3 <- dataLM(AoValv, time, respVar = "AoGradient", timeVar = "time", 
                 evTimeVar = "EvTime", summary = "area", tranfFun = sqrt)
    # Cox models
    CoxLM1 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + sqrt(AoGradient), 
                    data = D1)
    CoxLM2 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + sqrt(AoGradient) + 
                        slope, data = D2)
    CoxLM3 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + area, data = D3)
    list("LM_value" = CoxLM1, "LM_value-slope" = CoxLM2, "LM_cumulative" = CoxLM3)
}

###################
# Extract Results #
###################

extract_coefs <- function (time) {
    models <- LM_models_fun(time)
    fCI <- function (x) {
        g <- function (y) {
            y <- sprintf("%.3f", y)
            out <- paste("(", paste(y, collapse = "; "), ")", sep = "")
            gsub("-", "$-$", out, fixed = TRUE)
        }
        if (is.matrix(x)) apply(x, 1, g) else g(x)
    }
    data.frame(
        time = c(time, rep(NA, length(coef(models[[2]])) - 1)),
        Value = c(coef(models[[1]]), NA), 
        "95\\% CI" = c(fCI(confint(models[[1]])), NA),
        #
        Value = coef(models[[2]]), 
        "95\\% CI" = fCI(confint(models[[2]])),
        #
        Value = c(coef(models[[3]]), NA), 
        "95\\% CI" = c(fCI(confint(models[[3]])), NA),
        check.names = FALSE)
}

dT <- do.call("rbind", lapply(c(5.5, 7.5, 9.5), extract_coefs))
cap <- paste("Estimated coefficients and 95\\% confidence intervals for the parameters", 
             "of the Cox models fitted to landmark datasets at follow-up times $t = 5.5,",
             "7.5,$", "9.5 years.")
print(xtable(dT, label = "Tab:Res-LM", caption = cap, align = c("l", rep("r", 7))), 
      math.style.negative = TRUE, sanitize.text.function = function (x) x)

##########################################################################################
##########################################################################################
##########################################################################################

# Mixed Model Landmarking

dataLM <- JMbayes:::dataLM

LMmixed_models_fun <- function (time) {
    # create LM data
    active_data <- AoValv[AoValv$EvTime > time & AoValv$time <= time, ]
    # fit mixed model
    lmeFit <- lme(sqrt(AoGradient) ~ 0 + TypeOp + TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), 
                  data = active_data, 
                  random = list(id = pdDiag(form = ~ ns(time, k = c(2.5, 6), B = c(0.5, 13)))))
    betas <- fixef(lmeFit)
    b <- data.matrix(ranef(lmeFit))
    id <- match(active_data$id, unique(active_data$id))
    # calculate fitted values at event item
    Xvalue <- model.matrix(~ 0 + TypeOp + TypeOp:ns(time, k = c(2.5, 6), B = c(0.5, 13)), data = active_data)
    Zvalue <- model.matrix(~ ns(time, k = c(2.5, 6), B = c(0.5, 13)), data = active_data)
    active_data$value <- c(Xvalue %*% betas) + rowSums(Zvalue * b[id, ])
    Xslope <- model.matrix(dForm$fixed, data = active_data)
    Zslope <- model.matrix(dForm$random, data = active_data)
    active_data$slope <- c(Xslope %*% betas[dForm$indFixed]) + 
        rowSums(Zslope * b[id, dForm$indRandom])
    Xarea <- model.matrix(iForm$fixed, data = active_data)
    Zarea <- model.matrix(iForm$random, data = active_data)
    active_data$area <- c(Xarea %*% betas[iForm$indFixed]) + 
        rowSums(Zarea * b[id, iForm$indRandom])
    # create landmark data set
    dataLM <- JMbayes:::dataLM
    DD <- dataLM(active_data, time, respVar = "value", timeVar = "time", 
                 evTimeVar = "EvTime")
    # Fit Cox models
    CoxLM1 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + value, data = DD)
    CoxLM2 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + value + slope, 
                    data = DD)
    CoxLM3 <- coxph(Surv(EvTime, event) ~ TypeOp + Age + sex + area, data = DD)
    list(models = list("LMmixed_value" = CoxLM1, "LMmixed_value-slope" = CoxLM2, 
                       "LMmixed_cumulative" = CoxLM3), data = active_data)
}

###################
# Extract Results #
###################

extract_coefs <- function (time) {
    models <- LMmixed_models_fun(time)$models
    fCI <- function (x) {
        g <- function (y) {
            y <- sprintf("%.3f", y)
            out <- paste("(", paste(y, collapse = "; "), ")", sep = "")
            gsub("-", "$-$", out, fixed = TRUE)
        }
        if (is.matrix(x)) apply(x, 1, g) else g(x)
    }
    data.frame(
        time = c(time, rep(NA, length(coef(models[[2]])) - 1)),
        Value = c(coef(models[[1]]), NA), 
        "95\\% CI" = c(fCI(confint(models[[1]])), NA),
        #
        Value = coef(models[[2]]), 
        "95\\% CI" = fCI(confint(models[[2]])),
        #
        Value = c(coef(models[[3]]), NA), 
        "95\\% CI" = c(fCI(confint(models[[3]])), NA),
        check.names = FALSE)
}

dT <- do.call("rbind", lapply(c(5.5, 7.5, 9.5), extract_coefs))
cap <- paste("Estimated coefficients and 95\\% confidence intervals for the parameters", 
             "of the Cox models fitted to landmark datasets at follow-up times $t = 5.5,",
             "7.5,$", "9.5 years.")
print(xtable(dT, label = "Tab:Res-LM", caption = cap, align = c("l", rep("r", 7))), 
      math.style.negative = TRUE, sanitize.text.function = function (x) x)


