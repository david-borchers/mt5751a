library(RMark)
data("edwards.eberhardt")
eedat = edwards.eberhardt

nocc = length(strsplit(eedat[1,1],split="")[[1]])
n = nrow(eedat)
ch = data.frame(matrix(rep(NA,length=(n*nocc)),nrow=n))
colnames(ch) = paste(rep("occ",nocc),1:nocc,sep="")
for(i in 1:n) {
  ch[i,] = as.integer(strsplit(eedat[i,1],split="")[[1]][1:nocc])
}


require(ggcorrplot)

est.M0 = fit.cr(ch, model="M0",start.n0=60,start.p=0.1)
#est.Mt = fit.cr(ch, model="Mt",start.n0=60,start.p=0.1)
est.Mb = fit.cr(ch, model="Mb",start.n0=60,start.p=0.1)
#est.Mtb = fit.cr(ch, model="Mtb",start.n0=60,start.p=0.1)
est.Mh2 = fit.cr(ch, model="Mh2",start.mu1=0.04,start.mu2=0.2,start.phi=0.85)

aics = c(est.M0$AIC,est.Mb$AIC,est.Mh2$AIC)
names(aics) = c("M0","Mb","Mh2")
sort(aics-min(aics))
est.M0$Nhat
est.Mb$Nhat
est.Mh2$Nhat

est.M0$phat
est.Mb$phat
est.Mh2$phat

ggcorrplot(est.M0$beta.corrmatrix, type="lower",lab=TRUE)
ggcorrplot(est.Mb$beta.corrmatrix, type="lower",lab=TRUE)
ggcorrplot(est.Mh2$beta.corrmatrix, type="lower",lab=TRUE)

edwards.eberhardt = ch
save(edwards.eberhardt,file="./data/edwards_eberhardt.rda")
