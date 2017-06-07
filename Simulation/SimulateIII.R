##########################################################################################
# Aim: Simulates data from simulation Scenrio III                                        #
#                                                                                        #
# Required packages: MASS and splines from CRAN                                          #
#                                                                                        #
# Author: Dimitris Rizopoulos                                                            #
##########################################################################################

library("MASS")
library("splines")
n <- 1000 # number of subjects
K <- 15  # number of planned repeated measurements per subject, per outcome
t.max <- 15 # maximum follow-up time

################################################

# parameters for the linear mixed effects model
betas <- c("Group0" = 3.4554, "Group1" = 2.9470, 
           "Group0:Time1" = 1.0027, "Group1:Time1" = 0.9709, 
           "Group0:Time2" = 4.1290, "Group1:Time2" = 4.0893,
           "Group0:Time3" = 6.2182, "Group1:Time3" = 6.6909)
sigma.y <- 0.564 # measurement error standard deviation


# parameters for the survival model
gammas <- c("(Intercept)" = -5.7296, "Group" = 0.46) # coefficients for baseline covariates
alpha <- 0.0365 # association parameter value
phi <- 0.9518 # shape for the Weibull baseline hazard
meanCens0 <- 10 # mean of the uniform censoring distribution for group 0
meanCens1 <- 14 # mean of the uniform censoring distribution for group 1

D <- matrix(c(0.5686193, 0.2126076, 0.1547322, 0.4354939,
              0.2126076, 1.6721086, 2.3299235, 2.1926166,
              0.1547322, 2.329923, 5.0230656, 2.8873934,
              0.4354939, 2.1926166, 2.8873934, 4.0286104), 4, 4)
D <- (D + t(D)) / 2

################################################

Bkn <- c(0, 13)
kn <- c(2.5, 6)

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
DF <- data.frame(year = times, group = factor(rep(group, each = K)))
X <- model.matrix(~ 0 + group + group:ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
Z <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group)

################################################

#simulate random effects
b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)

# simulate event times
eta.t <- as.vector(W %*% gammas)
invS <- function (t, u, i) {
    h <- function (s) {
        group0 <- 1 - group[i]
        group1 <- group[i]
        ###
        iBS <- ins(s, knots = kn, Boundary.knots = Bkn)
        XXi <- cbind(group0*s, group1*s, group0*iBS[, 1], group1*iBS[, 1],
                     group0*iBS[, 2], group1*iBS[, 2], group0*iBS[, 3], group1*iBS[, 3])
        ZZi <- cbind(s, iBS)
        f <- as.vector(XXi %*% betas + rowSums(ZZi * b[rep(i, nrow(ZZi)), ]))
        exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f * alpha)
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
        tries <- tries - 1
        Up <- Up + 200
        Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}
na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
group <- group[na.ind]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from an exponential distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- numeric(n)
Ctimes[group == 0] <- runif(sum(group == 0), 0, 2 * meanCens0)
Ctimes[group == 1] <- runif(sum(group == 1), 0, 2 * meanCens1)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################


# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
names(dat) <- c("time", "group", "id", "y", "Time", "event")

#summary(tapply(id, id, length))
#table(event)
#n
#mean(event)

set <- sample(unique(id), 400)
train_data <- dat[!dat$id %in% set, ]
train_data$id <- match(train_data$id, unique(train_data$id))
test_data <- dat[dat$id %in% set, ]
test_data$id <- match(test_data$id, unique(test_data$id))

trueValues <- list(betas = betas, phi = phi, gammas = gammas, alpha = alpha,
                   b_test = b[sort(set), ], Bkn = Bkn, kn = kn)

# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W,
   betas, sigma.y, gammas, alpha, eta.t, eta.y, phi, t.max,
   trueTimes, u, Root, invS, D, b, K, set, dat,
   times, group, i, tries, Up, Bkn, kn, DF, meanCens0, meanCens1)

