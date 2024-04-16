# Replacing missing instruments from outcome GWAS with proxies
if(length(missing_IVs) == 0) {
  
  print("All exposure IVs found in outcome GWAS.")
  
} else {
  
  print("Some exposure IVs missing from outcome GWAS.")
  #out_full <- fread(paste0(deparse(substitute(out)), "_outcome.txt", sep = ""))
  
  for (i in 1:length(missing_IVs)) {
    
    proxies <- LDproxy(snp = missing_IVs[i], pop = "EUR", r2d = "r2", token = "6fb632e022ef", file = FALSE)
    proxies <- proxies[proxies$R2 > 0.8, ]
    proxy_present = FALSE
    
    if(length(proxies$RS_Number) == 0){
      
      print(paste0("No proxy SNP available for ", missing_IVs[i]))
      
    } else {
      
      for (j in 1:length(proxies$RS_Number)) {
        
        proxy_present <- proxies$RS_Number[j] %in% out_full$SNP
        
        if (proxy_present) {
          proxy_SNP = proxies$RS_Number[j]
          proxy_SNP_allele_1 = str_sub(proxies$Alleles[j], 2, 2)
          proxy_SNP_allele_2 = str_sub(proxies$Alleles[j], 4, 4)
          original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
          original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
          break
        }
      }
    }
    
    if(proxy_present == TRUE) {
      print(paste0("Proxy SNP found. ", missing_IVs[i], " replaced with ", proxy_SNP))
      proxy_row <- out_dat[1, ]
      proxy_row$SNP = missing_IVs[i]
      proxy_row$beta.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "BETA"])
      proxy_row$se.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "SE"])
      if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome = original_SNP_allele_1
      if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome = original_SNP_allele_2
      if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome = original_SNP_allele_1
      if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome = original_SNP_allele_2
      proxy_row$pval.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "P"])
      proxy_row$samplesize.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N"])
      if("N_case" %in% colnames(out_full)) proxy_row$ncase.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N_case"])
      if("N_control" %in% colnames(out_full))proxy_row$ncontrol.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N_control"])
      proxy_row$chr.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "chr.exposure"])
      proxy_row$pos.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "pos.exposure"])
      if("AF1" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "AF1"])
      out_dat <- rbind(out_dat, proxy_row)
    }
    
    if(proxy_present == FALSE) {
      print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
    }
  }
  
}