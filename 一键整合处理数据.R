library(TwoSampleMR)
library(gwasglue)
library(ieugwasr)
library(data.table)
exp_file_path<- "E:/EXP_DAT/"
out_file_path<- "D:/xlxz/jingshen/"
storge_file_path<- "D:/GutandDepression/"
a<- "/"
files<- list.files(exp_file_path)
files2<- list.files(out_file_path)
ord_mr<- ".mr.csv"
ord_or<- ".or.csv"
for (f in files2) {
  out<- fread(paste0(out_file_path,f))
  for (file in files) {
    exp<- fread(paste0(exp_file_path,file),header = T)
    exp<- subset(exp,P.weightedSumZ<1e-5)
    write.csv(exp,file = "exp.csv")
    exp<- read_exposure_data("exp.csv",
                             clump = F,
                             sep = ",",
                             snp_col = "rsID",
                             beta_col = "beta",
                             se_col = "SE",
                             pval_col = "P.weightedSumZ",
                             effect_allele_col = "eff.allele",
                             other_allele_col = "ref.allele")
    exp<- clump_data(dat = exp,
                     clump_kb = 10000,
                     clump_r2 = 0.001,
                     clump_p1 = 1,
                     clump_p2 = 1,
                     pop = "EUR")
    exp<- snp_add_eaf(exp,build = "37",pop = "EUR")
    out_dat<- merge(exp,out,by.x = "SNP",by.y = "rsids")
    write.csv(out_dat,file = "out.csv")
    out_dat<- read_outcome_data(filename = "out.csv",
                                snps = exp$SNP,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta",
                                se_col = "sebeta",
                                eaf_col = "af_alt",
                                effect_allele_col = "alt",
                                other_allele_col = "ref",
                                pval_col = "pval")
    res<- harmonise_data(exposure_dat = exp,outcome_dat = out_dat)
    mr_res<- mr(dat = res)
    mrTab<- generate_odds_ratios(mr_res)
    write.csv(mr_res,file = paste0(storge_file_path,f,a,file,ord_mr))
    write.csv(mrTab,file = paste0(storge_file_path,f,a,file,ord_or))
  }
}








