clean_data_fin <- apply(clean_data_fin,2,as.double)  # 
probes_names  <- extract_good_probe_list(clean_data_fin)

# function is used for filtering: the two thresholds we have used so far are probs = 0.9 and probs = 0.10

eset_ae_ad_p  <- clean_data_final[, colnames(clean_data_final) %in% c(probes_names, "case", "sex","age")]
eset_ae_ad_p  <- eset_ae_ad_p[complete.cases(eset_ae_ad_p),]



exp_df_male  <- subset(eset_ae_ad_p, eset_ae_ad_p$sex %in% c( "male","MALE","M","m",",male","Male"))
exp_df_fem   <- subset(eset_ae_ad_p, eset_ae_ad_p$sex %in% c("female", "FEMALE","F","f","Female"))

# Sex Stratified
results_limma_female_age   <- sva_cleaning_sex(exp_df_fem) 
results_limma_male_age     <- sva_cleaning_sex(exp_df_male)
results_limma_comb         <- sva_cleaning_mat(eset_ae_ad_p)


sv_male      <- results_limma_male_age[[2]]
male_limma   <- results_limma_male_age[[1]]
sv_comb      <- results_limma_comb[[2]]
comb_limma   <- results_limma_comb[[1]]
sv_female      <- results_limma_female_age[[2]]
male_limma   <- results_limma_female_age[[1]]

male_sv_df   <- as.data.frame(cbind(as.character(exp_df_male$case),
                                    exp_df_male$age, sv_male$sv))
comb_sv_df  <- as.data.frame(cbind(as.character(eset_ae_ad_p$case),eset_ae_ad_p$sex,
                                   eset_ae_ad_p$age, results_limma_comb[[2]]$sv))
female_sv_df   <- as.data.frame(cbind(as.character(exp_df_fem$case),
                                    exp_df_fem$age, sv_female$sv))


male_sv_names <- sprintf("sv-[%s]", seq(1:sv_male$n.sv))
colnames(sv_male$sv) <- male_sv_names
colnames(male_sv_df) <-c("case","age", colnames(sv_male$sv))

female_sv_names <- sprintf("sv-[%s]", seq(1:sv_female$n.sv))
colnames(sv_female$sv) <- female_sv_names
colnames(female_sv_df) <-c("case","age", colnames(sv_female$sv))


comb_sv_names <- sprintf("sv-[%s]", seq(1:sv_comb$n.sv))
colnames(sv_comb$sv) <- comb_sv_names

colnames(comb_sv_df) <- c("case","sex","age",colnames(sv_comb$sv))

write.csv(male_sv_df, "G118553_EC_sv_male.csv")
write.csv(comb_sv_df, "G188553_EC_sv_comb.csv")
write.csv(female_sv_df, "G118553 _EC_sv_fem.csv")

male_pheno <- subset(pheno_data_tissue, pheno_data_tissue$sex == "MALE")
female_pheno <- subset(pheno_data_tissue, pheno_data_tissue$sex =="FEMALE")

write.csv(male_pheno, "male_pheno_G11853_ec.csv")
write.csv(pheno_data_tissue, "comb_pheno_G11853_ec.csv")
write.csv(female_pheno, "female_pheno_G11853_ec.csv")