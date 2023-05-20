# This code was used to fit the maturation curves of nodal CT with age

library(mgcv)
Data <- read.csv("F:/Data/CBDP.csv") 
Measures <- read.csv("F:/Data/CBDP/cortical_thickness.csv", header = FALSE) 

scan_age = Data$scan_age
Data$sub_id <- as.factor(Data$sub_id)
subject=Data$sub_id
Data$sex <- as.factor(Data$sex)
sex=Data$sex

age_term_all<-matrix(NA,nrow=521,ncol=1000)
Pvalue_all <- matrix(NA,nrow=1000,ncol=1)
# outpath = 'F:/data/CBDP/results/'
for(roi in 1:1000){
  Data$Measure <-(Measures[,roi]) 
  Data$Measure <- as.numeric(Data$Measure)
  Measure=Data$Measure
  fit_G4 <- gam(Measure ~ s(scan_age, bs="tp", k=3) + sex + s(subject, bs="re"), data=Data, method="REML")
  # save(fit_G4,file=paste0(outpath,'fit_G',roi,'.Rdata'))
  
  summary_model <- summary(fit_G4)
  Pvalue_all[roi,1] = summary_model$s.table[1,4]
  age_term = predict(fit_G4, type = 'terms')[, 's(scan_age)']
  age_term_all[,roi]=age_term
 
}
write.table(age_term_all,"F:/data/CBDP/results/GAM_CT_age.csv",row.names=F, col.names=F)
write.table(Pvalue_all,"F:/data/CBDP/results/Pvalue.csv",row.names=F, col.names=F)





