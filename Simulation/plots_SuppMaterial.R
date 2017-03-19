##########################################################################################
# Aim: Read simulation results and produce the figures in the supplementary material     #
# Input: The user needs to specify the scenario for which the plots are required; options# 
#        are Ia, Ib, IIa, IIb, IIIa, IIIb. These correspond to Scenarios I, II and II    #
#        described in the paper, and the suffix 'a' or 'b' denote whether the correct of #
#        misspecified functional form for the time effect the linear mixed model were    #
#        used.                                                                           #
# Author: Dimitris Rizopoulos                                                            #
##########################################################################################


plot_fun <- function (scenario, diff = TRUE) {
    require("lattice")
    nam <- paste0("combinedResults", scenario, ".RData")
    con_text <- paste0("https://raw.github.com/drizopoulos/lm_and_jm/master/Simulation/",
                       nam)
    con <- url(con_text)
    load(con)
    close(con)
    ###############################################################
    out <- do.call("rbind", res)
    nams_aucs <- c("JM_AUC", "LM_AUC", "LMmixed_AUC")
    nams_pes <- c("JM_PE", "LM_PE", "LMmixed_PE")
    AUCs <- data.matrix(out[nams_aucs])
    PEs <- data.matrix(out[nams_pes])
    out <- if (diff) {
        AUCs <- cbind(AUCs[, "JM_AUC"] - AUCs[, "LM_AUC"],
                      AUCs[, "JM_AUC"] - AUCs[, "LMmixed_AUC"])
        PEs <- cbind(PEs[, "JM_PE"] - PEs[, "LM_PE"],
                     PEs[, "JM_PE"] - PEs[, "LMmixed_PE"])
        data.frame(
            "functional_form" = rep(out$`functional form`, 2),
            "time" = factor(rep(out$time, 2), labels = paste("t =", c(5.5, 7.5, 9.5))),
            "AUC" = as.vector(AUCs),
            "PE" = as.vector(PEs),
            "model" = rep(c("JM vs LM", "JM\nvs\nLM mixed"), each = nrow(out))
        )
    } else {
        data.frame(
            "functional_form" = rep(out$`functional form`, 3),
            "time" = factor(rep(out$time, 3), labels = paste("t =", c(5.5, 7.5, 9.5))),
            "AUC" = c(AUCs),
            "PE" = c(PEs),
            "model" = rep(c("JM", "LM", "LM\nmixed"), each = nrow(out))
        )
    }
    out$functional_form <- factor(out$functional_form, 
                                  labels = c("value", "value & slope", "cumulative"))
    obj_auc <- bwplot(AUC ~ model | time * functional_form, data = out, diff = diff, 
                      panel = function (..., diff) {
                          panel.bwplot(..., col = "lightgrey", coef = 0, pch = "|")
                          if (diff) panel.abline(h = 0, lty = 2, col = 1)
                      }, as.table = TRUE, layout = c(3, 3), 
                      scales = list(y = list(relation = "free", tick.number = 4)),
                      ylab = if (diff) expression(paste(Delta, "AUC")) else "AUC",
                      par.strip.text = list(cex = 0.7),
                      par.settings = list(box.rectangle = list(fill = "lightgrey"),
                                          par.ylab.text = list(cex = 0.7), 
                                          axis.text = list(cex = 0.6)))
    obj_pe <- bwplot(PE ~ model | time * functional_form, data = out, diff = diff,
                     panel = function (...) {
                         panel.bwplot(..., col = "lightgrey", coef = 0, pch = "|")
                         if (diff) panel.abline(h = 0, lty = 2, col = 1)
                     }, as.table = TRUE, layout = c(3, 3), 
                     scales = list(y = list(relation = "free", tick.number = 4)),
                     ylab = if (diff) expression(paste(Delta, "PE")) else "PE",
                     par.strip.text = list(cex = 0.7),
                     par.settings = list(box.rectangle = list(fill = "lightgrey"),
                                         par.ylab.text = list(cex = 0.7), 
                                         axis.text = list(cex = 0.6)))
    list(AUC = obj_auc, PE = obj_pe, 
         mean_AUC = with(out, tapply(AUC, list(time, model, functional_form), 
                                     mean, na.rm = TRUE)),
         mean_PE = with(out, tapply(PE, list(time, model, functional_form), 
                                    mean, na.rm = TRUE))
    )
}

##########################################################################################

plots <- plot_fun("Ia")
plots$AUC
plots$PE

plots$mean_AUC
plots$mean_PE






