Ffliter =10
dat<- exp_dat_ld
N= dat[1,"samplesize.exposure"]
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
dat=transform(dat,F=(N-2)*R2/(1-R2))
outTab=dat[dat$F>Ffliter,]
write.csv(outTab,file = "GFR_ld.csv",row.names = F)

rm(N)
rm(Ffliter)
rm(dat)
rm(outTab)