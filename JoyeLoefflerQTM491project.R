#Power Analysis for a die dropping experiment

# Height dropped is 9.25 inches (the hight of our econometrics book). Dropped on wood panel floor.
test = c(3,5,5,5,3,3,5,1,3,3,2,2,6,3,4,2,1,2,5,6,1,1,4,1,1,2,5,5,4,6,6,2,2,4,2,6,1,2,1,5,3,2,6,5,6,1,6,2,1,1,2,5,5,2,6,2,4,4,1,2,1,2,2,4,6,4,1,2,5,6,6,3,5,5,2,5,2,5,3,2,6,4,6,3,6,4,6,1,5,2,5,2,2,5,4,1,5,3,3,6,1,4,5,5,2,5,5,2,1,3,1,6,5,1,2,5,2,1,4,4,4,1,4,5,5,6,6,2,5,2,1,2,5,5,4,1,5,4,2,2,3,5,2,3,4,4,3,2,1,1,1,4,2,4,1,5,3,2,5,4,5,2,1,5,5,2,5,2,1,1,5,5,3,5,2,6,6,4,3,5,4,1,3,1,1,2,6,3,1,4,2,4,1,5,6,3,6,3,6,3,4,3,4,3,4,1,3,2,2,5,6,6,2,5,4,3,5,3,4,4,4,3,4,2,5,5,2,4,2,3,6,5,5,1,3,6,6,2,1,1,3,3,3,6,6,3,3,3,6,6,1,5,6,5,4,1,4,3,2,3,5,3,6,2,4,4,3,1,3,5,4,1,6,5,6,1,3,3,1,1,6,5,2,5,3,1,1,2,1,1,5,6,1,6,6,3,3,2,5,3,3,3,4,4,3,5,6,5,1,2,1,1,1,1,4,4,1,4,6,6,2,2,3,6,3,3,4,2,4,4,6,4,6,2,1,2,3,1,4,3,3,3,5,1,4,3,2,3,4,3,6,6,2,5,5,1,4,4,6,3,6,5,5,2,5,4,2,6,4,6,1,1,5,4,6,4,6,3,2,3,5,6,6,5,4,6,1,6,6,1,2,2,1,5,2,4,4,1,4,5,2,4,1,4,1,6,3,1,4,1,5,1,1,3,1,2,3,6,3,3,2,1,2,5,6,4,5,6,3,1,2,6,2,5,3,1,2,4,3,3,6,6,2,4,6,2,3,1,6,5,5,3,3,1,5,6,2,5,3,2,5,3,1,2,2,4,2,3,4,1,4,5,4,4,3,6,3,2,2,6,6,3,6,4,6,1,4,5,5,3,6,1,5,4,1,1,3,4,6,2,3,1,1,3,5,3,4,4,4,2,4,4,4,6,2,1,3,1,5,1,5,4,3,5,6,1,2,3,5,2,1,6,2,5,3,3,4,3,3,4,5,1,2,5,4,4,6,6,2,3,6,5,4,2,5,2,4,1,5,1,3,2,6,1,5,4,3,2,6,5,2,2,6,4,1,3,6,4,4,4,4,5,2,6,5,1,1,6,1,4,3,3,1,5,1,5,6,4,3,6)
length(test)
ftable(test)
dif=mean(test)-3.5
delta=3-3.5
#n needed for power=.8 for found effect size
power1 = power.t.test(delta=dif,sd=sd(test),power=.8,sig.level=.05,
               type="two.sample",alternative="two.sided")  
#n needed for power=.8 for effect size of .5
power2 = power.t.test(delta=delta,sd=sd(test),power=.8,sig.level=.05,
             type="two.sample",alternative="two.sided")

t1 = t.test(x=test,y=NULL,alternative="two.sided",mu=3,conf.level=.95)
t2 = t.test(x=test,y=NULL,alternative="great",mu=3,conf.level=.95)
t3 = t.test(x=test,y=NULL,alternative='less',mu=3.5,conf.level=.95)

install.packages("pwr")
library(pwr)

prop6 = 93/600
h = ES.h(1/12,1/6)
h2 = ES.h(.15,1/6)

#n needed for power=.8 for effect size of halving 1/6
power3 = pwr.p.test(h=h,sig.level=.05,power=.8,alternative="two.sided")
#n needed for power=.8 for found effect size 
power4 = pwr.p.test(h=prop6-(1/6),sig.level=.05,power=.8,alternative="two.sided")
#n needed for power=.8 for effect size of decreasing 1/6 by 10%
power5 = pwr.p.test(h=h2,sig.level=.05,power=.8,alternative="two.sided")
#power of n=600 and for effect size of decreasing 1/6 by 50%
power6 = pwr.p.test(h=(1/12)-(1/6),sig.level=.05,n=600,alternative="two.sided")
#power of n=600 and for effect size of decreasing 1/6 by 10%
power7 = pwr.p.test(h=h2,sig.level=.05,n=600,alternative="two.sided")

prop.test(x=93,n=600,p=1/12,alternative="two.sided",conf.level=.95)
prop.test(x=93,n=600,p=1/12,alternative="great",conf.level=.95)
prop.test(x=93,n=600,p=1/6,alternative="less",conf.level=.95)
prop.test(x=93,n=600,p=.15,alternative="great",conf.level=.95)

sigma=sqrt(((1/12)*(11/12))/600)
#power curve for proportion of 6s
curve(pnorm(((sqrt(600)*(x-.15))/sigma)-qnorm(1-.05)),
      from=1/12,
      to=1/6,
      main="Power Curve for Proportion",
      xlab="Proportion of 6s",
      ylab="Power")

sigma1=sqrt((.15*.85)/600)
#power curve for proportion of 6s
curve(pnorm(((sqrt(600)*(x-.15))/sigma1)-qnorm(1-.05)),
      from=.15,
      to=1/6,
      main="Power Curve for Proportion",
      xlab="Proportion of 6s",
      ylab="Power")

#power curve for mean die roll
curve(pnorm(((sqrt(600)*(x-3))/sd(test))-qnorm(1-.05)),
      from=3,
      to=3.5,
      main="Power Curve for Mean",
      xlab="Average Die Roll",
      ylab="Power")


