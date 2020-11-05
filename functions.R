extract_good_probe_list <- function(dataset) {
  # dataset - expression dataset as dataframe
  # probe_percentile_threshold - percentile at which to use as cut-off for detected probes 
  # number of samples in which probe must be expressed in - fixed at 0.8 - i.e 80% of samples
  # calculate quantile threshold for each sample
  # samples as columns and genes as rows
  dataset <- t(dataset)
  sample_quantiles <- apply(dataset, 2, quantile, probs = 0.10 ) # Find the 90th percentile of for each sample
  sample_quantiles <- replicate(nrow(dataset), sample_quantiles)
  sample_quantiles <- t(sample_quantiles)
  # convert to dataframe
  dataset_check     <- dataset- sample_quantiles
  dataset_check     <- as.data.frame(dataset_check)
  dataset_replace   <- dataset_check %>% mutate_each(funs(replace(.,.<0,NA)))
  
  rownames(dataset_replace) <- rownames(dataset)
  count_na_row              <- apply(dataset_replace,1,function(x){sum(is.na(x))})
  dataset_check_2           <- as.data.frame(cbind(dataset, count_na_row))
  dataset_check_2$pct       <- dataset_check_2$count_na_row > 0.8 * ncol(dataset)
  dataset_check_2           <- subset(dataset_check_2, dataset_check_2$pct == FALSE)
  
  # subset good probese
  good_probes <- rownames( dataset_check_2)
  
  return(good_probes)
}
filter_by_variance <- function(dataset){
  tdataset     <- t(dataset)
  tdataset     <- as.data.frame(tdataset)
  rownames(tdataset) <- colnames(dataset)
  colnames(tdataset) <- rownames(dataset)
  tdataset$varexp  <- apply(dataset,2,var)
  tdataset$varexp  <- as.numeric(tdataset$varexp)
  x <- quantile(tdataset$varexp, prob = 0.5)
  tdataset <- subset(tdataset, tdataset$varexp > x)
  return(as.data.frame(t(tdataset)))
  
}


linear_reg_genes <- function(var_name){
  exp_mat        <- clean_mat_ad
  exp_mat$Y      <- exp_mat[,var_name]
  exp_mat$case <- ifelse(exp_mat$case %in% c("normal","Normal","control","Control","con"),0,1)
  ln_regs <- lm(Y ~ case, data = exp_mat)
  return(ln_regs)
} 

sva_cleaning_mat <- function(edat){
  # edat is the expression matrix with sex and case as the last two columns 
  mod_data   <- model.matrix( ~1 + case +age+sex, data = edat)
  mod0       <- model.matrix(~1 +age+sex, data = edat)
  n_cols     <- ncol(edat) - 3
  p_cols     <- ncol(edat) - 1
  expr_ae_ad <- apply(edat[,1:n_cols],2,as.numeric)
  #expr_ae_ad <- as.matrix(eset_ae_ad_pb[,1:2063])
  #expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)
  
  n.survar   <- num.sv( t(expr_ae_ad),mod_data, method = "leek") 
  
  svobj  <- sva(t(expr_ae_ad), mod_data, mod0, n.sv=n.survar, method="two-step")
  if(svobj$n.sv >1 ){
    # adjust for sva
    X     <- cbind(mod_data, svobj$sv)
    fit   <- lmFit(t(expr_ae_ad), X )
    contrasts.matrix <- cbind("C1" = c(1,0,0,0,replicate(svobj$n.sv,0)),"C2" = c(1,1,0,0,replicate(svobj$n.sv,0)))
    fit_contrats <-  contrasts.fit(fit, contrasts.matrix)
    eB           <- eBayes(fit)
    return_table <- topTable(eB,2, adjust.method = "fdr",confint = TRUE, number = n_cols )
    
  }
  else{
    
    X <- mod_data
    fit <- lmFit(t(expr_ae_ad),X)
    contrasts.matrix <- cbind("C1"= c(1,0,0,0), "C2" = c(1,1,0,0)
    )
    fit_contrasts <- contrasts.fit(fit,contrasts.matrix)
    eB            <- eBayes(fit)
    
    return_table <- topTable(eB,2,  adjust.method = "fdr", confint = TRUE, number = n_cols)}
  return(list(return_table,svobj))
} 
edat <- exp_df_fem

sva_cleaning_sex <- function(edat){
  # edat is a combination of tranpose of the expression matrix with phenotype data(case, sex and age)
  edat$case  <- relevel(edat$case, ref = "Normal")
  mod_data   <- model.matrix( ~1 + case +age , data = edat)
  mod0       <- model.matrix(~1 +age, data = edat)
  n_cols     <- ncol(edat) - 3
  p_cols     <- ncol(edat) - 1
  expr_ae_ad_sv <- apply(edat[,1:n_cols],2,as.numeric)
  
  #expr_ae_ad <- as.matrix(eset_ae_ad_pb[,1:2063])
  #expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)
  
  n.survar   <- num.sv( t(expr_ae_ad_sv),mod_data, method = "be") 
  
  svobj  <- sva(t(expr_ae_ad_sv), mod_data, mod0, n.sv=n.survar)
  if(svobj$n.sv >1 ){
    # adjust for sva
    X     <- cbind(mod_data, svobj$sv)
    fit   <- lmFit(t(expr_ae_ad_sv), X )
    
    eB           <- eBayes(fit)
    return_table <- topTable(eB,2, adjust.method = "fdr",confint = TRUE, number = n_cols )
    
  }
  else{
    
    X <- mod_data
    fit <- lmFit(t(expr_ae_ad_sv),X)
    eB            <- eBayes(fit)
    
    return_table <- topTable(eB,2,  adjust.method = "fdr", confint = TRUE, number = n_cols)}
  return(list(return_table,svobj))
  }