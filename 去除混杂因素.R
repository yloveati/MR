#去除混杂因素，即筛选出来同样作用于结局的SNP
library(MendelianRandomization)
dat<- read.csv("GFR_ld.csv",header = T,sep = ",",check.names = F)
snpld<- dat$SNP
y<- seq_along(snpld)
chunks<- split(snpld,ceiling(y/100))
outTab= data.frame()
for(i in names(chunks)){
  confounder=phenoscanner(
    snpquery = chunks[[i]],
    catalogue = "GWAS",
    pvalue = 1e-05,
    proxies = "None",
    r2 = 0.8,
    build = 37)
  outTab=rbind(outTab,
               confounder$results)
}
delSnp=c("rs13078960","rs2030323")
dat=dat[!dat$SNP%in%delSnp,]