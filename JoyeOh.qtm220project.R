#Replication of "Adolescent marijuana use higher in states with medical 
#marijuana laws?" by Wall et. al (2011), with an extension of repeating the 
#analysis including additional data from 2009 to 2011

#Ryan Joye and Chris Oh
#December 19, 2015

# Loop for assigning MML passed
for (i in 1:357) {
  if (nsduh02to08$state[i] == "Alaska" || nsduh02to08$state[i] == "California" || 
      nsduh02to08$state[i] == "Colorado" || nsduh02to08$state[i] == "Hawaii" ||
      nsduh02to08$state[i] == "Maine" || nsduh02to08$state[i] == "Nevada" ||
      nsduh02to08$state[i] == "Oregon" || nsduh02to08$state[i] == "Washington") {
    nsduh02to08$mml[i] = 1
  }
  if (nsduh02to08$state[i] == "Michigan" && nsduh02to08$year[i] >= 2008) {
    nsduh02to08$mml[i] = 1
  }
  if (nsduh02to08$state[i] == "Montana" && nsduh02to08$year[i] >= 2004) {
    nsduh02to08$mml[i] = 1
  }
  if (nsduh02to08$state[i] == "New Mexico" && nsduh02to08$year[i] >= 2007) {
    nsduh02to08$mml[i] = 1
  }
  if (nsduh02to08$state[i] == "Rhode Island" && nsduh02to08$year[i] >= 2006) {
    nsduh02to08$mml[i] = 1
  }
  if (nsduh02to08$state[i] == "Vermont" && nsduh02to08$year[i] >= 2004) {
    nsduh02to08$mml[i] = 1
  }
}

# Diff in diff regression
attach(nsduh02to08)
lm.fit = lm(prevalence~mml+as.factor(state)+as.factor(year), data=nsduh02to08)
summary(lm.fit)

#WLS
lm.fit2 = lm(prevalence~mml+as.factor(state)+as.factor(year), 
             data=nsduh02to08, weights=nsduh02to08$n)
summary(lm.fit2)

#2002-2011
for (i in 1:510) {
  if (nsduh02to11$state[i] == "Alaska" || nsduh02to11$state[i] == "California" || 
        nsduh02to11$state[i] == "Colorado" || nsduh02to11$state[i] == "Hawaii" ||
        nsduh02to11$state[i] == "Maine" || nsduh02to11$state[i] == "Nevada" ||
        nsduh02to11$state[i] == "Oregon" || nsduh02to11$state[i] == "Washington") {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "Michigan" && nsduh02to11$year[i] >= 2008) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "Montana" && nsduh02to11$year[i] >= 2004) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "New Mexico" && nsduh02to11$year[i] >= 2007) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "Rhode Island" && nsduh02to11$year[i] >= 2006) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "Vermont" && nsduh02to11$year[i] >= 2004) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "Arizona" && nsduh02to11$year[i] >= 2010) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "District of Columbia" && nsduh02to11$year[i] >= 2010) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "Delaware" && nsduh02to11$year[i] >= 2011) {
    nsduh02to11$mml[i] = 1
  }
  if (nsduh02to11$state[i] == "New Jersey" && nsduh02to11$year[i] >= 2010) {
    nsduh02to11$mml[i] = 1
  }
}

# Diff in diff regression
attach(nsduh02to11)
lm.fit3 = lm(prevalence~mml+as.factor(state)+as.factor(year), data=nsduh02to11)
summary(lm.fit3)

#Weighted Least Squares
lm.fit4 = lm(prevalence~mml+as.factor(state)+as.factor(year), 
             data=nsduh02to11, weights=nsduh02to11$n)
summary(lm.fit4)

#Synthetic Control Methods
install.packages("Synth")
library("Synth")

#MI
nsduh02to11$state = as.character(nsduh02to11$state)
dataprep.out_MI = dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2008,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "Michigan",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2008,
  time.plot = 2002:2011)

dataprep.out_MI$X1
dataprep.out_MI$X0
dataprep.out_MI$Z1
dataprep.out_MI$Z0

synth.out_MI = synth(data.prep.obj = dataprep.out_MI)

gaps_MI = dataprep.out_MI$Y1plot - (dataprep.out_MI$Y0plot %*% synth.out_MI$solution.w)
gaps_MI

synth.tables_MI = synth.tab(dataprep.res = dataprep.out_MI, 
                             synth.res = synth.out_MI)
names(synth.tables_MI)
synth.tables_MI$tab.pred
synth.tables_MI$tab.v
synth.tables_MI$tab.w
synth.tables_MI$tab.loss

path.plot(synth.res = synth.out_MI, dataprep.res = dataprep.out_MI,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("Michigan",
          "synthetic Michigan"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_MI, dataprep.res = dataprep.out_MI,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "Michigan", tr.intake=2008)

#Montana
dataprep.out_MT = dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2004,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "Montana",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2004,
  time.plot = 2002:2011)

dataprep.out_MT$X1
dataprep.out_MT$X0
dataprep.out_MT$Z1
dataprep.out_MT$Z0

synth.out_MT = synth(data.prep.obj = dataprep.out_MT)

gaps_MT = dataprep.out_MT$Y1plot - (dataprep.out_MT$Y0plot %*% synth.out_MT$solution.w)
gaps_MT

synth.tables_MT = synth.tab(dataprep.res = dataprep.out_MT, 
                            synth.res = synth.out_MT)
synth.tables_MT$tab.pred
synth.tables_MT$tab.v
synth.tables_MT$tab.w
synth.tables_MT$tab.loss

path.plot(synth.res = synth.out_MT, dataprep.res = dataprep.out_MT,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("Montana",
                                      "synthetic Montana"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_MT, dataprep.res = dataprep.out_MT,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "Montana", tr.intake=2004)

#New Mexico
dataprep.out_NM = dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2007,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "New Mexico",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2007,
  time.plot = 2002:2011)

synth.out_NM = synth(data.prep.obj = dataprep.out_NM)

gaps_NM = dataprep.out_MT$Y1plot - (dataprep.out_NM$Y0plot %*% synth.out_NM$solution.w)
gaps_NM

synth.tables_NM = synth.tab(dataprep.res = dataprep.out_NM, 
                            synth.res = synth.out_NM)
synth.tables_NM$tab.pred
synth.tables_NM$tab.v
synth.tables_NM$tab.w
synth.tables_NM$tab.loss

path.plot(synth.res = synth.out_NM, dataprep.res = dataprep.out_NM,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("New Mexico",
                                      "synthetic New Mexico"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_NM, dataprep.res = dataprep.out_NM,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "New Mexico", tr.intake=2007)

#Rhode Island
dataprep.out_RI = dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2006,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "Rhode Island",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2006,
  time.plot = 2002:2011)

synth.out_RI = synth(data.prep.obj = dataprep.out_RI)

path.plot(synth.res = synth.out_RI, dataprep.res = dataprep.out_RI,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("Rhode Island",
                                      "synthetic Rhode Island"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_RI, dataprep.res = dataprep.out_RI,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "Rhode Island", tr.intake=2006)

#Vermont
dataprep.out_VT= dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2004,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "Vermont",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2004,
  time.plot = 2002:2011)

synth.out_VT = synth(data.prep.obj = dataprep.out_VT)

path.plot(synth.res = synth.out_VT, dataprep.res = dataprep.out_VT,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("Vermont",
                                      "synthetic Vermont"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_VT, dataprep.res = dataprep.out_VT,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "Vermont", tr.intake=2004)

#Arizona
dataprep.out_AZ= dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2010,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "Arizona",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2010,
  time.plot = 2002:2011)

synth.out_AZ = synth(data.prep.obj = dataprep.out_AZ)

path.plot(synth.res = synth.out_AZ, dataprep.res = dataprep.out_AZ,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("Arizona",
                                      "synthetic Arizona"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_AZ, dataprep.res = dataprep.out_AZ,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "Arizona", tr.intake=2010)

#District of Columbia
dataprep.out_DC= dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2010,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "District of Columbia",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2010,
  time.plot = 2002:2011)

synth.out_DC = synth(data.prep.obj = dataprep.out_DC)

path.plot(synth.res = synth.out_DC, dataprep.res = dataprep.out_DC,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("District of Columbia",
                                      "synthetic District of Columbia"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_DC, dataprep.res = dataprep.out_DC,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "District of Columbia", tr.intake=2010)

#Delaware
dataprep.out_DE= dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2011,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "Delaware",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2011,
  time.plot = 2002:2011)

synth.out_DE = synth(data.prep.obj = dataprep.out_DE)

path.plot(synth.res = synth.out_DE, dataprep.res = dataprep.out_DE,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("Delaware",
                                      "synthetic Delaware"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_DE, dataprep.res = dataprep.out_DE,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "Delaware", tr.intake=2011)

#New Jersey
dataprep.out_NJ= dataprep(
  foo = nsduh02to11,
  predictors = "prevalence",
  predictors.op = "mean",
  time.predictors.prior = 2002:2010,
  dependent = "prevalence",
  unit.names.variable = "state",
  unit.variable = "state.no",
  time.variable = "year",
  treatment.identifier = "New Jersey",
  controls.identifier = c("Alabama","Arkansas","Connecticut","Florida","Georgia","Idaho",
                          "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana",
                          "Maryland","Massachusetts","Minnesota","Mississippi","Missouri",
                          "Nebraska","New Hampshire","New York","North Carolina",
                          "North Dakota","Ohio","Oklahoma","Pennsylvania","South Carolina",
                          "South Dakota","Tennessee","Texas","Utah","Virginia",
                          "West Virginia","Wisconsin","Wyoming"),
  time.optimize.ssr = 2002:2010,
  time.plot = 2002:2011)

dataprep.out_NJ
synth.out_NJ = synth(data.prep.obj = dataprep.out_NJ)

path.plot(synth.res = synth.out_NJ, dataprep.res = dataprep.out_NJ,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(0, 15), Legend = c("New Jersey",
                                      "synthetic New Jersey"), Legend.position = "bottomright")
gaps.plot(synth.res = synth.out_NJ, dataprep.res = dataprep.out_NJ,
          Ylab = "prevalence of marijuana usage among 12-17 year-olds", Xlab = "year",
          Ylim = c(-2.5, 2.5), Main = "New Jersey", tr.intake = 2010)
