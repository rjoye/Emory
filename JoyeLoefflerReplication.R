setwd("/Users/redbirds36/Downloads/dataverse_files")

rm(list = ls())
set.seed(2074530682)
install.packages('AER')
library(AER)
install.packages('RItools')
library(RItools)
d.kn <- read.csv('study1.csv', stringsAsFactors = FALSE)
B <- 1000
cat('Note: Bootstraps are set to', B, 'replicates.\n')

#FUNCTIONS
resampleInStratum <- function(val, data, strata) {
  ok <- which(data$strata == val)
  sample(ok, length(ok), replace = TRUE)
}

resample <- function(data) {
  vals <- unique(data$strata)
  do.call(c, lapply(vals, resampleInStratum, data))
}

boot <- function(formula, FUN, data, fitted_model, B) {
  resamples <- do.call(cbind, lapply(1:B, function(b) resample(data)))
  boots <- matrix(NA, nrow = B, ncol = length(coef(fitted_model)))
  boots <- as.data.frame(boots)
  names(boots) <- names(coef(fitted_model))
  for (b in 1:B) {
    resampled_data <- data[resamples[, b],]
    coefs <- coef(FUN(formula, data = resampled_data, weights = weights))
    for (var in names(coefs)) {
      boots[b, var] <- coefs[var]
    }
  }
  boots
}

cace <- function(response, baseline, data, B = 5, strata = NULL, ipw = FALSE,
                 conditioning_covariate = NULL, incl_sameparty = FALSE) {
  
  y <- data[, paste(response, '01', sep = '.')]
  data <- data[!is.na(y), ]
  
  if (!ipw) {
    data$weights <- 1
  } else {
    data$weights <- data[, paste('ipw', response, sep = '.')]
  }
  
  if (is.null(strata)) {
    data$strata <- 1
  }  else {
    data$strata <- data[, strata]
  }
  
  if (is.null(conditioning_covariate)) {
    if (incl_sameparty) {
      formula <- paste(response, '.01',
                       '~attended.session + same.party +', baseline, 
                       '.NArecode | . - attended.session + assigned.to.session', 
                       sep = '')
    } else {
      formula <- paste(response, '.01',
                       '~attended.session +', baseline, 
                       '.NArecode | . - attended.session + assigned.to.session', 
                       sep = '')
    }
  } else {
    formula <- paste(response, '.01',
                     '~attended.session * same.party +', baseline, 
                     '.NArecode | . - attended.session * same.party + assigned.to.session * same.party', 
                     sep = '')
  }
  fitted_model <- ivreg(as.formula(formula), data = data, weights = weights)
  boots <- boot(formula, ivreg, data, fitted_model, B)
  out <- list(treat = 'attended.session', fitted_model = fitted_model, boots = boots)
  class(out) <- 'est'
  out
}

itt <- function(response, baseline, data, B = 5, strata = NULL, ipw = FALSE, incl_sameparty = FALSE) {
  
  y <- data[, paste(response, '01', sep = '.')]
  data <- data[!is.na(y), ]
  
  if (!ipw) {
    data$weights <- 1
  } else {
    data$weights <- data[, paste('ipw', response, sep = '.')]
  }
  
  if (is.null(strata)) {
    data$strata <- 1
  }  else {
    data$strata <- data[, strata]
  }
  
  if (incl_sameparty) {
    formula <- paste(response, '.01',
                     '~ assigned.to.session + same.party +', baseline, 
                     '.NArecode', sep = '')
  } else {
    formula <- paste(response, '.01',
                     '~ assigned.to.session +', baseline, 
                     '.NArecode', sep = '')
    
  }
  fitted_model <- lm(as.formula(formula), data = data, weights = weights)
  boots <- boot(formula, lm, data, fitted_model, B)
  out <- list(treat = 'assigned.to.session', fitted_model = fitted_model, boots = boots)
  class(out) <- 'est'
  out
}

summary.est <- function(object) {
  tr <- object$treat
  fm <- object$fitted_model
  boots <- object$boots[, tr]
  B <- length(boots)
  N <- nrow(fm$model)
  EST <- coef(fm)[tr]
  SE <- bootSE(boots, N, B)
  P <- bootP(boots)
  q025 <- quantile(boots, .025)
  q975 <- quantile(boots, .975)
  out <- data.frame(B = EST, SE = SE, P = P, q025 = q025, q975 = q975, N = N)  
  rownames(out) <- terms(fm)[[2]]
  round(out, 3)
}

bootSE <- function(x, n, B) {
  sq <- (x - mean(x)) ^ 2
  sqrt(n / (n - 1)) * sqrt(sum(sq) / (B - 1))
}

bootP <- function(x) {
  2 * min(mean(x < 0), mean(x > 0))
}

NArecode <- function(x) {
  out <- x
  out[is.na(x)] <- -1
  factor(out)
}

ConditionalTrimmingBounds <- function(y, x, z) {
  #  Calculates conditional trimming bounds by subsetting by a covariate,
  #    performing trimming bounds analysis for values of y, z at each
  #    covariate value, and averaging results
  #
  #  Args:
  #    y: response variable
  #    x: covariate
  #    z: treatment variable
  #
  #  Returns:
  #    a pair of bounds on the treatment effect
  x.recode <- x  # Recode missing covariate values and coerce to factor
  x.recode[is.na(x)] <- max(x, na.rm=TRUE) + 1
  x.recode <- factor(x.recode)
  
  # Function to calc trimming bounds for (y, z) at value of x.recode
  f <- function(val) TrimmingBounds(y[x.recode == val], z[x.recode == val])
  
  # Apply f at each value of x.recode, bind results
  tbs <- do.call(rbind, lapply(levels(x.recode), f))
  
  # Drop undefined bounds, calculate weighted means, return
  w <- tabulate(x.recode) / length(x.recode)  # Calculate weights
  c(mean(tbs[, 1], weights = w, na.rm = TRUE),
    mean(tbs[, 2], weights = w, na.rm = TRUE))
}

ComplierReporter <- function(response, baseline, data) {
  data <- subset(data, assigned.to.session == 1)
  y1 <- data[, response]
  y0 <- data[, baseline]
  cr <- which(data$attended.session == 1 & !is.na(y1))
  cn <- which(data$attended.session == 1 &  is.na(y1))
  nr <- which(data$attended.session == 0 & !is.na(y1))
  nn <- which(data$attended.session == 0 &  is.na(y1))
  type <- rep('', nrow(data))
  type[cr] <- 'cr'
  type[cn] <- 'cn'
  type[nr] <- 'nr'
  type[nn] <- 'nn'
  out <- data.frame(mean.cr = rep(NA, 2), 
                    mean.cn = rep(NA, 2), 
                    mean.nr = rep(NA, 2), 
                    mean.nn = rep(NA, 2), 
                    Fstat = rep(NA, 2), 
                    Pvalue = rep(NA, 2))
  out[1, 1:4] <- c(mean(y0[cr], na.rm = TRUE), 
                   mean(y0[cn], na.rm = TRUE), 
                   mean(y0[nr], na.rm = TRUE), 
                   mean(y0[nn], na.rm = TRUE))
  out[2, 1:4] <- c(length(cr), 
                   length(cn), 
                   length(nr), 
                   length(nn))
  fm <- aov(y0 ~ type)
  out[1, 5] <- summary(fm)[[1]]$'F value'[1]
  out[1, 6] <- summary(fm)[[1]]$'Pr(>F)'[1]
  out
}

TrimmingBounds <- function(y, z) {
  #  Calculates trimming bounds by replacing missing values of y with
  #    extreme (lo & hi) observed values for opposite treatment value
  #
  #  Args:
  #    y: response variable
  #    z: treatment variable
  #
  #  Returns:
  #    a pair of bounds on the treatment effect
  
  # if all observations are either treated or control, trimming bounds
  #   are not defined
  if (sum(z) %in% c(0, length(z))) {
    return(c(NA, NA))
  }
  pT <- mean(!is.na(y[z == 1]))    # percent treated that are observed
  pC <- mean(!is.na(y[z == 0]))    # percent control that are observed
  zobs <- z[!is.na(y)]             # treatments for observed
  yobs <- y[!is.na(y)]             # responses for observed
  yTobs <- sort(yobs[zobs == 1])   # observed treated y's, sorted lo values 1st
  yCobs <- sort(yobs[zobs == 0])   # observed control y's, sorted lo values 1st
  if (pC > pT) {                   # if observed more controls than treated...
    q <- (pC - pT) / pC            # ratio of extra missingness among treated
    nobs <- length(yCobs)          # number of observed control values
    yClo <- tail(yCobs, -ceiling(q * nobs))      # trim lo control y-values
    yChi <- head(yCobs, ceiling((1 - q) * nobs)) # trim hi control y-values
    lb <- mean(yTobs) - mean(yClo) # average effect trimming lo control y-values
    ub <- mean(yTobs) - mean(yChi) # average effect trimming hi control y-values
  } else {                         # if observed more controls than treated...
    q <- (pT - pC) / pT            # ratio of extra missingness among control
    nobs <- length(yTobs)          # number of observed treated values
    yTlo <- tail(yTobs, -ceiling(q * nobs))      # trim lo treated y-values
    yThi <- head(yTobs, ceiling((1 - q) * nobs)) # trim hi treated y-values
    lb <- mean(yThi) - mean(yCobs) # average effect trimming hi treated y-values
    ub <- mean(yTlo) - mean(yCobs) # average effect trimming lo treated y-values
  }
  # return bounds
  if (is.na(lb) | is.na(ub)) {
    c(NA, NA)
  } else {
    if (lb < ub) {
      c(lb, ub)
    } else {
      c(ub, lb)
    }
  }
}

summarize <- function(x) {
  c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), sum(is.na(x)))
}

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /  (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

loo <- function(response, baseline, data = d1, strata = 'congressional.district') {
  # leave one district out ivreg
  formula <- paste(response, '.01',
                   '~attended.session + same.party +', baseline, 
                   '.NArecode | . - attended.session + assigned.to.session', 
                   sep = '')
  formula <- as.formula(formula)
  fitted_model <- ivreg(formula, data = data)
  data <- data[-fitted_model$na.action, ]
  n <- nrow(data)
  if (!is.null(strata)) {
    strata <- data[, strata]
    vals <- unique(strata)
  } else {
    strata <- 1
    vals <- 1
  }
  LOOEST <- function(val) {
    loo_data <- subset(data, strata != val)
    coef(ivreg(formula, data = loo_data))
  }
  loos <- do.call(rbind, lapply(vals, LOOEST))
  counts <- do.call(c, lapply(vals, function(val) sum(strata == val)))
  list(fitted_model = fitted_model, loos = loos, counts = counts, n = n)
}

SummarizeLOO <- function(object) {
  treat <- 'attended.session'
  LOO_ESTs <- object$loos[, treat]
  EST <- coef(object$fitted_model)[treat]
  LOO_SE <- sqrt(sum((1 - object$counts / object$n) * (EST - LOO_ESTs) ^ 2)) 
  LOO_P <- 2 * (1 - pnorm(abs(EST / LOO_SE)))
  out <- data.frame(EST = EST, SE = LOO_SE, P = LOO_P, N = object$n)
  rownames(out) <- terms(object$fitted_model)[[2]]
  round(out, 3)
}

SummarizeConditionalEffects <- function(object) {
  coefs <- coef(object$fitted_model)
  N <- nrow(object$fitted_model$model)
  EST1 <- coefs['attended.session'] + coefs['attended.session:same.party']
  EST2 <- coefs['attended.session']
  boots1 <- object$boots[, 'attended.session'] + object$boots[, 'attended.session:same.party']
  boots2 <- object$boots[, 'attended.session']
  B <- length(boots1)
  SE1 <- bootSE(boots1, N, B)
  SE2 <- bootSE(boots2, N, B)
  P1 <- bootP(boots1)
  P2 <- bootP(boots2)
  q025 <- c(quantile(boots1, .025), quantile(boots2, .025))
  q975 <- c(quantile(boots1, .975), quantile(boots2, .975))
  P_compare <- 2 * min(mean(boots2 < boots1), mean(boots2 > boots1))
  results <- data.frame(EST = c(EST1, EST2), SE = c(SE1, SE2), P = c(P1, P2),
                        q025 = q025, q975 = q975, N = c(N, NA), P_compare = c(P_compare, NA))
  round(results, 3)
}

#DATA CLEANING
d1 <- data.frame(id = 10000 + as.numeric(rownames(d.kn)))
d1$congressional.district <- d.kn$geo.grp

# Willingness to participate, treatment assignment, and compliance
d1$willing.to.participate <- 1 * (d.kn$slfrpt %in% 1:2)
d1$assigned.to.info.only <- 1 * (d.kn$exp.grp %in% 4)
d1$assigned.to.session <- 1 * (d.kn$exp.grp %in% 1:2)
d1$attended.session <- 1 * (d.kn$attendance == 1)
.ok <- (d1$assigned.to.info.only == 0 & d1$assigned.to.session == 0) | 
  d1$willing.to.participate == 0
d1$assigned.to.session[.ok] <- d1$attended.session[.ok] <- NA

# Demographics
d1$Female <- 1 * (d.kn$ppgender == 2)
d1$Age <- d.kn$ppage
d1$Education <- d.kn$ppeducat # education == 1 for no HS, 2 for HS, 3 for some college, 4 for college+
d1$White <- 1 *(d.kn$ppethm == 1)
d1$Black  <- 1 *(d.kn$ppethm == 2)
d1$Latino <- 1 *(d.kn$ppethm == 4)
d1$Other <- 1 *(d.kn$ppethm %in% c(3, 5))
d1$Income <- NA
d1$Income[d.kn$ppincimp %in% 1:3] <- 1 # less than $10k
d1$Income[d.kn$ppincimp %in% 4:5] <- 2 # $10k to $15k
d1$Income[d.kn$ppincimp == 6] <- 3 # $15k to $20k
d1$Income[d.kn$ppincimp == 7] <- 4 # $20k to $25k
d1$Income[d.kn$ppincimp == 8] <- 5 # $25k to $30k
d1$Income[d.kn$ppincimp %in% 9:10] <- 6 # $30k to $40k
d1$Income[d.kn$ppincimp == 11] <- 7 # $40k to $50k
d1$Income[d.kn$ppincimp %in% 12:15] <- 8 # $50k to $100k
d1$Income[d.kn$ppincimp %in% 16:19] <- 10 # $100k or more

# Party/Position groups
d1$Party7pt <- d.kn$partyid7
.constituent.rep <- d.kn$partyid3 == 1
.constituent.dem <- d.kn$partyid3 == 3
.constituent.rep[is.na(d.kn$partyid3)] <- 0
.constituent.dem[is.na(d.kn$partyid3)] <- 0
.member.rep <- d.kn$geo.grp %in% c(3, 7, 8, 12, 16)
.member.dem <- !.member.rep
d1$same.party <- 1 * ((.constituent.rep & .member.rep) | (.constituent.dem & .member.dem))
.anti.immig <- d.kn$geo.grp %in% c(3, 7, 12, 13, 16)
.pro.immig <- !.anti.immig

# Path to citizenship
d1$path.to.citizenship.baseline <- with(d.kn,
                                        6 * (bas32 == 1 & bas32c == 2) +
                                          5 * (bas32 == 1 & bas32c == 1) +
                                          4 * (bas32 == 3 & bas32b == 1) +
                                          3 * (bas32 == 3 & bas32b == 3) +
                                          2 * (bas32 == 3 & bas32b == 2) +
                                          1 * (bas32 == 2 & bas32c == 1) +
                                          0 * (bas32 == 2 & bas32c == 2))
d1$path.to.citizenship.baseline[.anti.immig] <- 6 - d1$path.to.citizenship.baseline[.anti.immig]
d1$path.to.citizenship.followup <- with(d.kn,
                                        6 * (fol32 == 1 & fol32c == 2) +
                                          5 * (fol32 == 1 & fol32c == 1) +
                                          4 * (fol32 == 3 & fol32b == 1) +
                                          3 * (fol32 == 3 & fol32b == 3) +
                                          2 * (fol32 == 3 & fol32b == 2) +
                                          1 * (fol32 == 2 & fol32c == 1) +
                                          0 * (fol32 == 2 & fol32c == 2))
d1$path.to.citizenship.followup[.anti.immig] <- 6 - d1$path.to.citizenship.followup[.anti.immig]

# Agree on Increase Legal Immigration
d1$legal.immigration.baseline <- with(d.kn,
                                      4 * (bas30a == 1 & bas30b == 2) +
                                        3 * (bas30a == 1 & bas30b == 1) +
                                        2 * (bas30a == 2) +
                                        1 * (bas30a == 3 & bas30b == 1) +
                                        0 * (bas30a == 3 & bas30b == 2))
d1$legal.immigration.baseline[.anti.immig] <- 4 - d1$legal.immigration.baseline[.anti.immig]
d1$legal.immigration.followup <- with(d.kn,
                                      4 * (fol30a == 1 & fol30b == 2) +
                                        3 * (fol30a == 1 & fol30b == 1) +
                                        2 * (fol30a == 2) +
                                        1 * (fol30a == 3 & fol30b == 1) +
                                        0 * (fol30a == 3 & fol30b == 2))
d1$legal.immigration.followup[.anti.immig] <- 4 - d1$legal.immigration.followup[.anti.immig]

# Trust MOC
d1$trust.baseline <- 3 - (d.kn$bas17 - 1)
d1$trust.baseline[d.kn$bas17 == 5] <- NA  # set DK's to unobserved
d1$trust.followup <- 3 - (d.kn$fol17 - 1)
d1$trust.followup[d.kn$fol17 == 5] <- NA  # set DK's to unobserved

# Approve MOC
d1$approval.baseline <- 4 - (d.kn$bas22 - 1)
d1$approval.baseline[d.kn$bas22 == 7] <- NA # set DK's to unobserved
d1$approval.followup <- 4 - (d.kn$fol22 - 1)
d1$approval.followup[d.kn$fol22 == 7] <- NA # set DK's to unobserved

# Intent to Vote
d1$turnout.baseline <- 4 - (d.kn$bas20 - 1)
d1$turnout.baseline[d.kn$bas20 %in% 6:7] <- 0
d1$vote.baseline <- 4 - (d.kn$bas21 - 1)
d1$vote.baseline[d.kn$bas21 %in% 6:7] <- 0
d1$vote.followup <- 4 - (d.kn$fol21 - 1)
d1$vote.followup[d.kn$fol21 %in% 6:7] <- 0

# Actual Vote
d1$vote.november <- 1 * (d.kn$nov3 == 2 & d.kn$nov5 == 1)
d1$turnout.november <- 1 * (d.kn$nov3 == 2)

# Drop 'true controls' who were not assigned to either attend session or receive information
d1 <- subset(d1, !is.na(attended.session))

# recode baseline on 0-1 scale (to estimate Null SDs for Cohen's d's)
d1$path.to.citizenship.baseline.01 <- scale01(d1$path.to.citizenship.baseline)
d1$legal.immigration.baseline.01 <- scale01(d1$legal.immigration.baseline)
d1$trust.baseline.01 <- scale01(d1$trust.baseline)
d1$approval.baseline.01 <- scale01(d1$approval.baseline)
d1$vote.baseline.01 <- scale01(d1$vote.baseline)
d1$turnout.baseline.01 <- scale01(d1$turnout.baseline)

# recode baseline to include missingness, then create factor
d1$path.to.citizenship.baseline.NArecode <- NArecode(d1$path.to.citizenship.baseline)
d1$legal.immigration.baseline.NArecode <- NArecode(d1$legal.immigration.baseline)
d1$trust.baseline.NArecode <- NArecode(d1$trust.baseline)
d1$approval.baseline.NArecode <- NArecode(d1$approval.baseline)
d1$vote.baseline.NArecode <- NArecode(d1$vote.baseline)
d1$turnout.baseline.NArecode <- NArecode(d1$turnout.baseline)

# put responses on 0-1 scale
d1$path.to.citizenship.followup.01 <- scale01(d1$path.to.citizenship.followup)
d1$legal.immigration.followup.01 <- scale01(d1$legal.immigration.followup)
d1$trust.followup.01 <- scale01(d1$trust.followup)
d1$approval.followup.01 <- scale01(d1$approval.followup)
d1$vote.followup.01 <- scale01(d1$vote.followup)
d1$vote.november.01 <- scale01(d1$vote.november)
d1$turnout.november.01 <- scale01(d1$turnout.november)

# Calculate inverse probability of response weights
prob <- glm(is.na(path.to.citizenship.followup) ~ 
              path.to.citizenship.baseline.NArecode + 
              assigned.to.session + same.party, 
            data = d1, family = binomial)$fitted.values
d1$ipw.path.to.citizenship.followup <- 1 / prob
prob <- glm(is.na(legal.immigration.followup) ~ 
              legal.immigration.baseline.NArecode + 
              assigned.to.session + same.party, 
            data = d1, family = binomial)$fitted.values
d1$ipw.legal.immigration.followup <- 1 / prob
prob <- glm(is.na(trust.followup) ~ 
              trust.baseline.NArecode + 
              assigned.to.session + same.party, 
            data = d1, family = binomial)$fitted.values
d1$ipw.trust.followup <- 1 / prob
prob <- glm(is.na(approval.followup) ~ 
              approval.baseline.NArecode + 
              assigned.to.session + same.party, 
            data = d1, family = binomial)$fitted.values
d1$ipw.approval.followup <- 1 / prob
prob <- glm(is.na(vote.followup) ~ 
              vote.baseline.NArecode + 
              assigned.to.session + same.party, 
            data = d1, family = binomial)$fitted.values
d1$ipw.vote.followup <- 1 / prob
prob <- glm(is.na(vote.november) ~ vote.baseline.NArecode + assigned.to.session 
            + same.party, data = d1, family = binomial)$fitted.values
d1$ipw.vote.november <- 1 / prob
prob <- glm(is.na(turnout.november) ~ turnout.baseline.NArecode + assigned.to.session 
            + same.party, data = d1, family = binomial)$fitted.values
d1$ipw.turnout.november <- 1 / prob


# main analysis

#CACE
data <- d1
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)

#COMPLIER REPORTER
tableS2 <- data.frame(Complier_Reporter = rep(NA, 12),
                      Complier_Nonreporter = rep(NA, 12),
                      Noncomplier_Reporter = rep(NA, 12),
                      Noncomplier_Nonreporter = rep(NA, 12),
                      F = rep(NA, 12),
                      P = rep(NA, 12))
tableS2[ 1: 2, ] <- ComplierReporter('path.to.citizenship.followup', 
                                     'path.to.citizenship.baseline.01', d1)
tableS2[ 3: 4, ] <- ComplierReporter('legal.immigration.followup', 
                                     'legal.immigration.baseline.01', d1)
tableS2[ 5: 6, ] <- ComplierReporter('trust.followup', 'trust.baseline.01', d1)
tableS2[ 7: 8, ] <- ComplierReporter('approval.followup', 
                                     'approval.baseline.01', d1)
tableS2[ 9:10, ] <- ComplierReporter('vote.followup', 'vote.baseline.01', d1)
tableS2[11:12, ] <- ComplierReporter('vote.november', 'vote.baseline.01', d1)

# reporting
#TABLE 2
table2 <- data.frame(CACE = rep(NA, 6), 
                     SE = rep(NA, 6),  
                     P = rep(NA, 6),
                     q025 = rep(NA, 6),  
                     q975 = rep(NA, 6),  
                     N = rep(NA, 6))
table2[1, ] <- summary(CACE_S1_path)
table2[2, ] <- summary(CACE_S1_immg)
table2[3, ] <- summary(CACE_S1_trus)
table2[4, ] <- summary(CACE_S1_appr)
table2[5, ] <- summary(CACE_S1_vint)
table2[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2$d <- round(table2$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                      'Approval', 'Vote Intent', 'Actual Vote')
print(table2[, -(4:5)])

rownames(tableS2) <- c('Path to Citizenship', '1', 'Legal Immigration', '2', 'Trust', '3', 
                       'Approval', '4', 'Vote Intent', '5', 'Actual Vote', '6')

cat('\n\nTable S2\n')
print(round(tableS2, 3), digits = 3)

#Extension
#WOMEN
data <- subset(d1, Female==1)
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
table2women <- data.frame(CACE = rep(NA, 6), 
                     SE = rep(NA, 6),  
                     P = rep(NA, 6),
                     q025 = rep(NA, 6),  
                     q975 = rep(NA, 6),  
                     N = rep(NA, 6))
table2women[1, ] <- summary(CACE_S1_path)
table2women[2, ] <- summary(CACE_S1_immg)
table2women[3, ] <- summary(CACE_S1_trus)
table2women[4, ] <- summary(CACE_S1_appr)
table2women[5, ] <- summary(CACE_S1_vint)
table2women[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2women$d <- round(table2women$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2women) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                      'Approval', 'Vote Intent', 'Actual Vote')
print(table2women[, -(4:5)])

#MEN
data <- subset(d1, Female==0)
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
table2men <- data.frame(CACE = rep(NA, 6), 
                          SE = rep(NA, 6),  
                          P = rep(NA, 6),
                          q025 = rep(NA, 6),  
                          q975 = rep(NA, 6),  
                          N = rep(NA, 6))
table2men[1, ] <- summary(CACE_S1_path)
table2men[2, ] <- summary(CACE_S1_immg)
table2men[3, ] <- summary(CACE_S1_trus)
table2men[4, ] <- summary(CACE_S1_appr)
table2men[5, ] <- summary(CACE_S1_vint)
table2men[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2men$d <- round(table2men$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2men) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                           'Approval', 'Vote Intent', 'Actual Vote')
print(table2men[, -(4:5)])

#HIGH INCOME
summary(d1$Income)
data <- subset(d1, Income>7)
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
table2high <- data.frame(CACE = rep(NA, 6), 
                          SE = rep(NA, 6),  
                          P = rep(NA, 6),
                          q025 = rep(NA, 6),  
                          q975 = rep(NA, 6),  
                          N = rep(NA, 6))
table2high[1, ] <- summary(CACE_S1_path)
table2high[2, ] <- summary(CACE_S1_immg)
table2high[3, ] <- summary(CACE_S1_trus)
table2high[4, ] <- summary(CACE_S1_appr)
table2high[5, ] <- summary(CACE_S1_vint)
table2high[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2high$d <- round(table2high$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2high) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                           'Approval', 'Vote Intent', 'Actual Vote')
print(table2high[, -(4:5)])

#LOW INCOME
data <- subset(d1, Income<=7)
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
table2low <- data.frame(CACE = rep(NA, 6), 
                         SE = rep(NA, 6),  
                         P = rep(NA, 6),
                         q025 = rep(NA, 6),  
                         q975 = rep(NA, 6),  
                         N = rep(NA, 6))
table2low[1, ] <- summary(CACE_S1_path)
table2low[2, ] <- summary(CACE_S1_immg)
table2low[3, ] <- summary(CACE_S1_trus)
table2low[4, ] <- summary(CACE_S1_appr)
table2low[5, ] <- summary(CACE_S1_vint)
table2low[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2low$d <- round(table2low$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2low) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                          'Approval', 'Vote Intent', 'Actual Vote')
print(table2low[, -(4:5)])


summary(d1$Age)
#OLDER THAN 44.5
data <- subset(d1, Age>44.5)
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
table2old <- data.frame(CACE = rep(NA, 6), 
                         SE = rep(NA, 6),  
                         P = rep(NA, 6),
                         q025 = rep(NA, 6),  
                         q975 = rep(NA, 6),  
                         N = rep(NA, 6))
table2old[1, ] <- summary(CACE_S1_path)
table2old[2, ] <- summary(CACE_S1_immg)
table2old[3, ] <- summary(CACE_S1_trus)
table2old[4, ] <- summary(CACE_S1_appr)
table2old[5, ] <- summary(CACE_S1_vint)
table2old[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2old$d <- round(table2old$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2old) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                          'Approval', 'Vote Intent', 'Actual Vote')
print(table2old[, -(4:5)])

#44.5 OR YOUNGER
data <- subset(d1, Age<=44.5)
strata <- 'congressional.district'
CACE_S1_path <- cace('path.to.citizenship.followup', 
                     'path.to.citizenship.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_immg <- cace('legal.immigration.followup', 
                     'legal.immigration.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_trus <- cace('trust.followup', 
                     'trust.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_appr <- cace('approval.followup', 
                     'approval.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_vint <- cace('vote.followup', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_actv <- cace('vote.november', 
                     'vote.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
CACE_S1_turn <- cace('turnout.november', 
                     'turnout.baseline', 
                     data, B, strata, incl_sameparty = TRUE)
table2young <- data.frame(CACE = rep(NA, 6), 
                        SE = rep(NA, 6),  
                        P = rep(NA, 6),
                        q025 = rep(NA, 6),  
                        q975 = rep(NA, 6),  
                        N = rep(NA, 6))
table2young[1, ] <- summary(CACE_S1_path)
table2young[2, ] <- summary(CACE_S1_immg)
table2young[3, ] <- summary(CACE_S1_trus)
table2young[4, ] <- summary(CACE_S1_appr)
table2young[5, ] <- summary(CACE_S1_vint)
table2young[6, ] <- summary(CACE_S1_actv)

Est_Null_SD <- c(
  sd(d1$path.to.citizenship.baseline.01, na.rm = TRUE),
  sd(d1$legal.immigration.baseline.01, na.rm = TRUE),
  sd(d1$trust.baseline.01, na.rm = TRUE),
  sd(d1$approval.baseline.01, na.rm = TRUE),
  sd(d1$vote.baseline.01, na.rm = TRUE),
  sd(d1$vote.november.01[d1$assigned.to.session == 0], na.rm = TRUE))
table2young$d <- round(table2young$CACE / Est_Null_SD, 3)

cat('\n\nTable 2\n')
rownames(table2young) <- c('Path to Citizenship', 'Legal Immigration', 'Trust', 
                         'Approval', 'Vote Intent', 'Actual Vote')
print(table2young[, -(4:5)])