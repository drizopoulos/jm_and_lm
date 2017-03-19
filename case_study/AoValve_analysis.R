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

con <- url("https://raw.github.com/drizopoulos/jm_and_lm/master/case_study/simulated_AoValve.RData")
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

