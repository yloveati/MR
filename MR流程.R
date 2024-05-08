library(TwoSampleMR)
library(VariantAnnotation)
library(ieugwasr)
library(gwasglue)
vcf<- readVcf("ieu-a-1284.vcf.gz")
vcf<- gwasvcf_to_TwoSampleMR(vcf = vcf,type = "exposure")
library(CMplot)
data<-vcf[,c("SNP","chr.exposure","pos.exposure","pval.exposure")]
colnames(data)<-c("SNP","CHR","BP","pvalue")
#绘制线性曼哈顿图
CMplot(data,plot.type = "m",LOG10 =TRUE, threshold = 5e-8, threshold.lwd = 3, threshold.lty = 1, signal.cex = 0.2, chr.den.col = NULL,cex = 0.2, bin.size = 1e5, ylim = c(0,50), file = "pdf", file.output = TRUE, width = 15, height = 9,verbose = TRUE)
#绘制圈图
CMplot(data,plot.type = "c", LOG10 = TRUE,threshold = 5e-8, threshold.lwd = 3, threshold.lty = 1, signal.cex = 0.2, chr.den.col = NULL,cex = 0.2, bin.size = 1e5, ylim = c(0,100), file = "pdf", file.output = TRUE, width = 7, height = 7,verbose = TRUE)
#根据p<5e-8筛选
exp<- subset(vcf,pval.exposure<5e-08)
#作为csv写出，并进行下一步LD
write.csv(exp,file = "GFR_exp.csv",row.names = F)
exp_dat<- read_exposure_data(filename = "GFR_exp.csv",
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             effect_allele_col = "effect_allele.exposure",
                             other_allele_col = "other_allele.exposure",
                             pval_col = "pval.exposure",
                             eaf_col = "eaf.exposure",
                             samplesize_col = "samplesize.exposure",
                             clump = F)
exp_dat_ld<- clump_data(exp_dat,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 1,clump_p2 = 1,pop = "EAS")
exp_dat_ld<- clump_data(dat = exp_dat,clump_kb = 10000,clump_r2 = 0.001,clump_p1 = 1,clump_p2 = 1,plink_bin = "plink.exe",bfile = "D:/human/EAS")
#计算F值
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

#导入暴露因素
exp_dat<- read_exposure_data(filename = "GFR_ld.csv",
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = "beta.exposure",
                             se_col = "se.exposure",
                             effect_allele_col = "effect_allele.exposure",
                             other_allele_col = "other_allele.exposure",
                             pval_col = "pval.exposure",
                             eaf_col = "eaf.exposure",
                             clump = F)
#导入结局因素
vcf<- readVcf("bbj-a-137.vcf.gz")
out<- gwasvcf_to_TwoSampleMR(vcf = vcf,type = "outcome")

intersection_dat<- merge(exp_dat,out,by.x="SNP",by.y="SNP")
write.csv(intersection_dat,file = "OP_out.csv")
#剔除工具变量与结局结果相关的因素p>0.05
intersection_dat_select<- subset(intersection_dat,pval.outcome>0.05)
write.csv(intersection_dat,file = "OP_out_select.csv")
out_dat<- read_outcome_data(snps = exp_dat$SNP,
                            filename = "OP_out.csv",
                            snp_col = "SNP",
                            beta_col = "beta.outcome",
                            sep = ",",
                            se_col = "se.outcome",
                            effect_allele_col = "effect_allele.outcome",
                            other_allele_col = "other_allele.outcome",
                            pval_col = "pval.outcome",
                            eaf_col = "eaf.outcome")

out_dat_select<- read_outcome_data(snps = exp_dat$SNP,
                            filename = "OP_out_select.csv",
                            snp_col = "SNP",
                            beta_col = "beta.outcome",
                            sep = ",",
                            se_col = "se.outcome",
                            effect_allele_col = "effect_allele.outcome",
                            other_allele_col = "other_allele.outcome",
                            pval_col = "pval.outcome",
                            eaf_col = "eaf.outcome")
#整合暴露数据与结局数据
res<- harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat)
res_sel<- harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat_select)
#输出用于MR的工具变量
outTab<-res[res$mr_keep.outcome=="TRUE",]
write.csv(outTab,file = "table.SNP.csv",row.names = F)
#检测异质性
heterRes<- mr_heterogeneity(res)
write.csv(heterRes,file = "table.heterogenety.csv",row.names = F)
#若结果p<0.05，则需要presso，NbDistribution越大，精度越高，计算机计算越多，min3000
presso=run_mr_presso(dat = res,NbDistribution = 3000)
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`,file = "table.MR-PRESSO.csv")
#异质性漏斗图
mr_funnel_plot(singlesnp_results =mr_singlesnp(res))

#多效性检测（SNP通过其他暴露因素影响到结局）：若存在多效性，则通过看文献了解存在哪些暴露因素，用phenoscanner网站剔除与暴露因素有关的SNP
pleioTab<- mr_pleiotropy_test(res)#若p<0.05说明存在多效性（有截距），则不做了
write.csv(pleioTab,file = "table.pleiotropy.csv",row.names = F)


#MR分析

mrRes<- mr(res, method_list = c("mr_ivw",
                                "mr_egger_regression",
                                "mr_simple_median",
                                "mr_penalised_weighted_median",
                                "mr_two_sample_ml"
                                ))
write.csv(mrRes,file = "table.mrRes.csv",row.names = F)
# OR值版本
mrTab<- generate_odds_ratios(mr_res = mrRes)
write.csv(mrTab,file = "table.mrTab.csv",row.names = F)
#绘制散点图
mr_scatter_plot(mr_results = mrRes,dat = res)
#森林图
mr_forest_plot(mr_singlesnp(res))
#leaveoneout
mr_leaveoneout(res)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(res))
