# ======================================
# Simulate for class practical
# ======================================
library(secr)

# get sundarbans fits
fits = readRDS("./inst/code4data/sundarbansFits.Rds")

# extract the models
fit0=fits$fit0
fit.D=fits$fit.D
fit.D.d2water=fits$fit.D.d2water
fit.d2water=fits$fit.d2water
fit.d2water_s=fits$fit.d2water_s
fit.d2water_l=fits$fit.d2water_l
fit0.h2.both=fits$fit0.h2.both
fit0.h2.sig=fits$fit0.h2.sig
fit0.h2.lam=fits$fit0.h2.lam

# look at AICs and choose best
AIC(fit0,fit.D,fit.D.d2water,fit.d2water,fit.d2water_s,fit.d2water_l,fit0.h2.both,fit0.h2.sig,fit0.h2.lam)

# create scr objects and simulate population
mask = fit0.h2.both$mask
Dsurf = predictDsurface(fit0.h2.both,mask)
#splotcovariate(Dsurf,covariate="D.0")
covariates(Dsurf)$D.0 = covariates(Dsurf)$D.0*2.75
set.seed(12345)
pop = sim.popn("D.0",core=Dsurf,model2D="IHP",buffer=0)
N = nrow(pop);N
plot(mask)
plot(pop,frame=FALSE,add=TRUE)

# Thin traps and simulate capture history
traps = as.matrix(traps(fit0.h2.both$capthist))
keep = c(1:29)*4
trapdf = as.data.frame(traps[keep,])
simtraps = read.traps(data=trapdf,detector="multi")
fitdetpars = detectpar(fit0.h2.both)
g0 = 0.2
sigma = 3000
set.seed(1)
simch = sim.capthist(simtraps,popn=pop,detectfn="HHN",noccasions=5,detectpar=list(g0=g0,sigma=sigma))
plot(simch,tracks=TRUE)
testfit = secr.fit(simch,mask=mask,detectfn="HHN")
N;region.N(testfit)

# turn secr capture history into object for use with mt5751a
crch = apply(simch,c(1,2),sum) # accumulate across traps (remove space from history)
summary(crch)
dim(crch)

# fit some models to check:
est.M0 = fit.cr(crch, model="M0")
est.Mt = fit.cr(crch, model="Mt")
est.Mb = fit.cr(crch, model="Mb")
est.Mtb = fit.cr(crch, model="Mtb")
est.Mh2 = fit.cr(crch, model="Mh2",start.mu1=0.8,start.mu2=0.25)

aics = c(est.M0$AIC,est.Mt$AIC,est.Mb$AIC,est.Mtb$AIC,est.Mh2$AIC)
names(aics) = c("M0","Mt","Mb","Mtb","Mh2")
sort(aics-min(aics))

est.M0$Nhat
est.Mb$Nhat
est.Mt$Nhat
est.Mtb$Nhat
est.Mh2$Nhat


est.M0$phat
est.Mb$phat
est.Mt$phat
est.Mh2$phat

library(ggcorrplot)
ggcorrplot(est.M0$beta.corrmatrix, type="lower",lab=TRUE)
ggcorrplot(est.Mt$beta.corrmatrix, type="lower",lab=TRUE)
ggcorrplot(est.Mb$beta.corrmatrix, type="lower",lab=TRUE)

# save as Rdata object for mt5751a
tigerch = crch
save(tigerch,file="./data/tigerch.rda")

# save as Rdata object for secr
tigersecrch = simch
tigermask = mask
tigerpop = pop
save(tigersecrch,tigermask,tigerpop,file="./data/tigersecr.rda")
