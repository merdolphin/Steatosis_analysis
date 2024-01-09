library(tidyverse)
library(readxl)
library(compareGroups)
library(VIM) ## library(VIM) aggr(df1, plot = F, numbers = T, prop = T)
library(FSA)
library(ggsci)
library(gridExtra)
library(ggpubr)
library("plotrix")
library(tidyverse)
library(ggplot2)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(caret)
library(party)
library(randomForest)
library(mice)
library(org.Hs.eg.db)
library(ggpubr)
library(biomaRt)
library(tximport)
library(clusterProfiler)
library(BiocParallel)

setwd("/home/lina/scratch/18.NASH/Rscripts")

ClinicalData <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                                         sheet = "Sheet1")
aggr(ClinicalData, plot=F)
colnames(ClinicalData) -> tmp1
write.csv(tmp1,"colnames.1.csv")
colnames.1 <- read.table("colnames.3.csv", sep=",")
colnames(ClinicalData) <- colnames.1
colnames(ClinicalData)
colSums(is.na(ClinicalData)) 

ClinicalData %>% dplyr::select(Diabetes, T2DM) %>% filter(Diabetes == 0) %>% as.data.frame() 

ClinicalData %>% dplyr::select(Metformin_treatment) %>% summary()

ClinicalData %>% dplyr::select(`Diabetes (0=No, 1=Yes, 2= Ex diabetic)`) %>% 
  mutate_if(is.numeric, as.factor) %>% group_by(`Diabetes (0=No, 1=Yes, 2= Ex diabetic)`) %>% 
  summary()
########add antilipid drug usage######




ClinicalData[,colSums(is.na(ClinicalData))<25] -> CD.1
CD.1 %>% colnames()
CD.1 %>% mutate(Gender = factor(Gender, levels=c(0,1), labels=c("Female","Male")), 
                Diabetes = as.factor(Diabetes),
                        Previous.gestational.DM = as.factor(Previous.gestational.DM),
                        Metformin_treatment = factor(Metformin_treatment, levels=c(0,1), labels=c("No","Yes")), 
                T2DM = as.factor(T2DM)
                        ) -> CD.2

CD.2 %>% dplyr::select(Age, Diabetes, Gender, `Weight`, `BMI(kg/m2)`,
                       MetSscore, MetS_Trig, MetS_BP, MetS_BSL) -> CD.2.anthropometricParm

CD.2 %>% dplyr::select(Diabetes, ALT, AST, `ALT/AST`, GGT, ALP, TotalBili, Albumin, 
                       Urea,SerumCreatine, eGFR ) -> CD.2.liver.renal

CD.2 %>% dplyr::select(Diabetes, TotalChol, HDL, LDL, `NonHDL-Chol`, Trig
                       ) -> CD.2.lipid.profile

CD.2.lipid.profile

CD.2 %>% dplyr::select(Diabetes, T2DM, Previous.gestational.DM, `HbA1c(%)`,
                       IFG, FBG, Cpeptide, Metformin_treatment,
                       Insulin, `HOMA2%B`, `HOMA2%S` ) -> CD.2.glucose.para

CD.2 %>% dplyr::select(Diabetes, NAFL, Histo_STEATOSIS, HISTO_NAS,
                       Histo_INFLAM, Histo_BALLOON, Histo_FIBROSIS, HISTO_NAS,
                       LiverFibrosisgrouping,
                       LiverGrossgrouping) %>% 
  mutate_if(is.numeric, as.factor) -> CD.2.NAFLD


meth = 1
createTable( compareGroups(Diabetes ~ Age + Gender +  `Weight` + 
                             `BMI(kg/m2)` + MetSscore + MetS_Trig + MetS_BP + MetS_BSL, 
                           data = CD.2.anthropometricParm, 
                           method = meth),  
             show.p.mul = TRUE, hide = c(Gender = "Male")) -> CD.2.anthropometricParm.T

CD.2.anthropometricParm.T

createTable( compareGroups(Diabetes ~ ALT + AST + #`ALT/AST` + 
                             GGT + ALP + TotalBili + Albumin + 
                           Urea + SerumCreatine + eGFR, 
                           data = CD.2.liver.renal, 
                           method = meth),  
             show.p.mul = TRUE) -> CD.2.liver.renal.T

createTable( compareGroups(Diabetes ~ T2DM + `HbA1c(%)` +
                             IFG + FBG + Cpeptide + Metformin_treatment +
                           Insulin + `HOMA2%B` + `HOMA2%S`, 
                           data = CD.2.glucose.para, 
                           method = meth),  
             show.p.mul = TRUE) -> CD.2.glucose.para.T

createTable( compareGroups(Diabetes ~ NAFL + Histo_STEATOSIS + HISTO_NAS + 
                             Histo_INFLAM + Histo_BALLOON + Histo_FIBROSIS +
                             HISTO_NAS + LiverFibrosisgrouping +
                             LiverGrossgrouping, 
                           data= CD.2.NAFLD,
                           method = meth),
             show.p.mul = TRUE) -> CD.2.NAFLD.T


createTable( compareGroups(Diabetes ~  TotalChol + HDL + `NonHDL-Chol` + LDL +
                              + Trig,
                           data= CD.2.lipid.profile,
                           method = 2),
             show.p.mul = TRUE) -> CD.2.lipid.profile.T

rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
)

rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
      )-> restab

export2csv(restab, file='table1.csv')

help("compareGroups")
help(createTable)

CD.2$LiverFibrosisgrouping


CD.2 %>% dplyr::select(Diabetes, LiverFibrosisgrouping) %>% 
  mutate(LiverFibrosisgrouping = as.factor(LiverFibrosisgrouping)) %>%
  group_by(Diabetes, LiverFibrosisgrouping) %>% na.omit() %>% summarise(n=n())

CD.2 %>% dplyr::select(Diabetes, Histo_STEATOSIS) %>% 
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>%
  group_by(Diabetes, Histo_STEATOSIS) %>% na.omit() %>% summarise(n=n())

CD.2 %>% dplyr::select(Diabetes, HISTO_NAS) %>% 
  mutate(HISTO_NAS = as.factor(HISTO_NAS)) %>%
  group_by(Diabetes, HISTO_NAS) %>% na.omit() %>% summarise(n=n())

################################################################
################################################################
##################
################## based on NAFL grouping
##################
###############################################################

### Ex-Diabetes as 1
ClinicalData$Diabetes[ClinicalData$Diabetes == "2"] <- "1"
ClinicalData$T2DM[ClinicalData$T2DM == "2"] <- "1"

ClinicalData

ClinicalData[,colSums(is.na(ClinicalData))<25] -> CD.1
CD.1 %>% duplicated()
CD.1 %>% mutate(Gender = factor(Gender, levels=c(0,1), labels=c("Female","Male")), 
                Diabetes = as.factor(Diabetes),
                Previous.gestational.DM = as.factor(Previous.gestational.DM),
                Metformin_treatment = factor(Metformin_treatment, levels=c(0,1), labels=c("No","Yes")), 
                T2DM = as.factor(T2DM)
) -> CD.2

CD.2 %>% mutate(NAFL = str_replace(NAFL,"No Path F1", "No Path")) %>% 
  mutate(NAFL = str_replace(NAFL, "No Path Inflam 1", "No Path")) %>%
  #mutate(NAFL = str_replace(NAFL, "NAFL", "NAFL+NASH")) %>%
  mutate(NAFL = str_replace(NAFL, "NASH", "NAFL")) %>%
  filter(! is.na(NAFL)) -> CD.2


CD.2 %>% dplyr::select(Age, NAFL, Gender, `Weight(kg)`, `BMI(kg/m2)`,
                       MetSscore, MetS_Trig, MetS_BP, MetS_BSL) -> CD.2.anthropometricParm

CD.2 %>% dplyr::select(NAFL, ALT, AST, `ALT/AST`, GGT, ALP, TotalBili, Albumin, 
                       Urea,SerumCreatine, eGFR ) -> CD.2.liver.renal

CD.2 %>% dplyr::select(NAFL, TotalChol, HDL, LDL, `NonHDL-Chol`, Trig
) -> CD.2.lipid.profile

CD.2.lipid.profile

CD.2 %>% dplyr::select(NAFL, Diabetes, T2DM, Previous.gestational.DM, 
                       IFG, FBG, Cpeptide, Metformin_treatment,
                       Insulin, `HOMA2%B`, `HOMA2%S` ) -> CD.2.glucose.para

CD.2 %>% dplyr::select(NAFL, NAFL, Histo_STEATOSIS, HISTO_NAS,
                       Histo_INFLAM, Histo_BALLOON, Histo_FIBROSIS, HISTO_NAS,
                       LiverFibrosisgrouping,
                       LiverGrossgrouping) %>% 
  mutate_if(is.numeric, as.factor) -> CD.2.NAFLD


meth = 1
createTable( compareGroups(NAFL ~ Age + Gender +  `Weight(kg)` + 
                             `BMI(kg/m2)` + MetSscore + MetS_Trig + MetS_BP + MetS_BSL, 
                           data = CD.2.anthropometricParm, 
                           method = meth),  
             show.p.mul = TRUE, hide = c(Gender = "Male")) -> CD.2.anthropometricParm.T

CD.2.anthropometricParm.T

createTable( compareGroups(NAFL ~ ALT + AST + #`ALT/AST` + 
                             GGT + ALP + TotalBili + Albumin + 
                             Urea + SerumCreatine + eGFR, 
                           data = CD.2.liver.renal, 
                           method = meth),  
             show.p.mul = TRUE) -> CD.2.liver.renal.T

createTable( compareGroups(NAFL ~ T2DM + 
                             IFG + FBG + Cpeptide + Metformin_treatment +
                             Insulin + `HOMA2%B` + `HOMA2%S`, 
                           data = CD.2.glucose.para, 
                           method = meth),  
             show.p.mul = TRUE) -> CD.2.glucose.para.T

createTable( compareGroups(NAFL ~ NAFL + Histo_STEATOSIS + HISTO_NAS + 
                             Histo_INFLAM + Histo_BALLOON + Histo_FIBROSIS +
                             HISTO_NAS + LiverFibrosisgrouping +
                             LiverGrossgrouping, 
                           data= CD.2.NAFLD,
                           method = meth),
             show.p.mul = TRUE) -> CD.2.NAFLD.T


createTable( compareGroups(NAFL ~  TotalChol + HDL + `NonHDL-Chol` + LDL +
                             + Trig,
                           data= CD.2.lipid.profile,
                           method = 2),
             show.p.mul = TRUE) -> CD.2.lipid.profile.T


rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
)










################################################################
################################################################
##################
################## based on Inflammation
##################
###############################################################


ClinicalData <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                          sheet = "Sheet1")

colnames.1 <- read.table("colnames.3.csv", sep=",")
colnames(ClinicalData) <- colnames.1
colnames(ClinicalData)
colSums(is.na(ClinicalData)) 


ClinicalData %>% dplyr::select(contains("Histo")) %>% mutate_if(is.numeric, as.factor) %>%
  summary()

ClinicalData$Diabetes[ClinicalData$Diabetes == "2"] <- "1"
ClinicalData$T2DM[ClinicalData$T2DM == "2"] <- "1"


ClinicalData[,colSums(is.na(ClinicalData))<25] -> CD.1


CD.1 %>% mutate(Gender = factor(Gender, levels=c(0,1), labels=c("Female","Male")), 
                Diabetes = as.factor(Diabetes),
                Previous.gestational.DM = as.factor(Previous.gestational.DM),
                Metformin_treatment = factor(Metformin_treatment, levels=c(0,1), labels=c("No","Yes")), 
                T2DM = as.factor(T2DM) 
                #Chol_Rx_name = as.factor(Chol_Rx_name)
) -> CD.2


CD.2 %>% gather(key, value, -Histo_INFLAM) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_INFLAM))) %>% arrange(p.value) %>% as.data.frame()

CD.2 %>% mutate(NAFL = str_replace(NAFL,"No Path F1", "No Path")) %>% 
  mutate(NAFL = str_replace(NAFL, "No Path Inflam 1", "No Path")) %>%
  #mutate(NAFL = str_replace(NAFL, "NAFL", "NAFL+NASH")) %>%
  #mutate(NAFL = str_replace(NAFL, "NASH", "NAFL")) %>%
  mutate(Histo_INFLAM=str_replace(Histo_INFLAM, "3","2")) %>% 
  mutate(Histo_FIBROSIS = str_replace(Histo_FIBROSIS,"1a","1")) %>%
  mutate(Histo_FIBROSIS = str_replace(Histo_FIBROSIS,"1b","1")) %>%
  mutate(Histo_FIBROSIS = str_replace(Histo_FIBROSIS,"1c","1")) %>%
  filter(! is.na(Histo_INFLAM)) -> CD.2
CD.2$b

CD.2 %>% dplyr::select(Age, Histo_INFLAM, Gender, `Weight`, `BMI(kg/m2)`,
                       MetSscore, MetS_Trig, MetS_BP, MetS_BSL) -> CD.2.anthropometricParm

CD.2 %>% dplyr::select(Histo_INFLAM, ALT, AST, `ALT/AST`, GGT, ALP, TotalBili, Albumin, Ferritin,
                       Urea,SerumCreatine, eGFR ) -> CD.2.liver.renal

CD.2 %>% dplyr::select(Histo_INFLAM, TotalChol, HDL, LDL, `NonHDL-Chol`, Trig#, Chol_Rx_name
) -> CD.2.lipid.profile

CD.2.anthropometricParm %>% dplyr::select(`BMI(kg/m2)`) %>% na.omit() %>% 
  summarize(mean=mean(`BMI(kg/m2)`), se=std.error(`BMI(kg/m2)`))



###################

CD.2 %>% dplyr::select(Histo_INFLAM, Diabetes, T2DM, Previous.gestational.DM, `HbA1c(%)`,
                       IFG, FBG, Cpeptide, Metformin_treatment,
                       Insulin, `HOMA2%B`, `HOMA2%S`, HOMA2IR ) %>%  mutate(`HOMA2%B` = as.numeric(`HOMA2%B`)) -> 
  CD.2.glucose.para

CD.2.glucose.para %>% group_by(Histo_INFLAM) %>% mutate(`HOMA2%B` = as.numeric(`HOMA2%B`)) %>%
  summarise_at(vars(`HOMA2%B`,`HOMA2%S`), list(name=mean, sd),na.rm=TRUE)


help("compareGroups")


CD.2 %>% dplyr::select(Histo_INFLAM, NAFL,  HISTO_NAS,
                       Histo_STEATOSIS, Histo_BALLOON, Histo_FIBROSIS, HISTO_NAS,
                       LiverFibrosisgrouping,
                       LiverGrossgrouping) %>% 
  mutate_if(is.numeric, as.factor) %>% 
  mutate(HISTO_NAS = as.numeric(HISTO_NAS)) %>%
  mutate(LiverGrossgrouping=as.numeric(LiverGrossgrouping))-> CD.2.NAFLD

CD.2.NAFLD

meth = 4
createTable( compareGroups(Histo_INFLAM ~ Age + Gender +  `Weight` + 
                             `BMI(kg/m2)` + MetSscore + MetS_Trig + MetS_BP + MetS_BSL, 
                           data = CD.2.anthropometricParm, 
                           method = meth), hide = c(Gender = "Male")) -> CD.2.anthropometricParm.T

CD.2.anthropometricParm.T

createTable( compareGroups(Histo_INFLAM ~ ALT + AST + #`ALT/AST` + 
                             GGT + ALP + TotalBili + Albumin + Ferritin +
                             Urea + SerumCreatine + eGFR, 
                           data = CD.2.liver.renal, 
                           method = meth)) -> CD.2.liver.renal.T

CD.2.liver.renal.T

createTable( compareGroups(Histo_INFLAM ~ Diabetes + T2DM + `HbA1c(%)` +
                             IFG + FBG + Cpeptide + Metformin_treatment +
                             Insulin + `HOMA2%B` + `HOMA2%S` + HOMA2IR, 
                           data = CD.2.glucose.para, 
                           method = meth)) -> CD.2.glucose.para.T

CD.2.glucose.para.T



compareGroups(Histo_INFLAM ~ Diabetes + T2DM + `HbA1c(%)` +
                IFG + FBG + Cpeptide + Metformin_treatment +
                Insulin + `HOMA2%B` + `HOMA2%S`, 
              data = CD.2.glucose.para, subset = Diabetes == "1",
              method = meth) %>% createTable()


createTable( compareGroups(Histo_INFLAM ~ NAFL + HISTO_NAS + 
                             Histo_INFLAM + Histo_BALLOON + Histo_FIBROSIS +
                             #LiverFibrosisgrouping +
                             LiverGrossgrouping, 
                           data= CD.2.NAFLD,
                           method = meth)) -> CD.2.NAFLD.T

CD.2.NAFLD %>% mutate(HISTO_NAS = as.numeric(HISTO_NAS)) %>% group_by(Histo_INFLAM) %>%
  summarize(mean=mean(HISTO_NAS), sd = sd(HISTO_NAS))

createTable( compareGroups( Histo_INFLAM ~  TotalChol + HDL + `NonHDL-Chol` + LDL +
                              + Trig,
                            data= CD.2.lipid.profile,
                            method = meth)) -> CD.2.lipid.profile.T


rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
)

rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
) -> restab




















################################################################
################################################################
##################
################## based on Steatosis
##################
###############################################################

ClinicalData <- read_xlsx("/home/lina/scratch/18.NASH/1.sample/HumanPatients_230216.xlsx",
                          sheet = "Sheet1")

colnames.1 <- read.table("colnames.3.csv", sep=",")
colnames(ClinicalData) <- colnames.1
colnames(ClinicalData)
colSums(is.na(ClinicalData)) 


ClinicalData %>% dplyr::select(contains("Histo")) %>% mutate_if(is.numeric, as.factor) %>%
  summary()

ClinicalData$Diabetes[ClinicalData$Diabetes == "2"] <- "1"
ClinicalData$T2DM[ClinicalData$T2DM == "2"] <- "1"

ClinicalData$Chol_Rx_name

#ClinicalData$Chol_Rx_name[is.na(ClinicalData$Chol_Rx_name)] <- "No"
#ClinicalData$Chol_Rx_name[ClinicalData$Chol_Rx_name==0] <- "No"
#ClinicalData$Chol_Rx_name[ClinicalData$Chol_Rx_name!="No"] <- "Yes"



ClinicalData[,colSums(is.na(ClinicalData))<25] -> CD.1


CD.1 %>% mutate(Gender = factor(Gender, levels=c(0,1), labels=c("Female","Male")), 
                Diabetes = as.factor(Diabetes),
                Previous.gestational.DM = as.factor(Previous.gestational.DM),
                Metformin_treatment = factor(Metformin_treatment, levels=c(0,1), labels=c("No","Yes")), 
                T2DM = as.factor(T2DM) 
                #Chol_Rx_name = as.factor(Chol_Rx_name)
) -> CD.2


CD.2 %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS))) %>% arrange(p.value) %>% as.data.frame()

CD.2 %>% mutate(NAFL = str_replace(NAFL,"No Path F1", "No Path")) %>% 
  mutate(NAFL = str_replace(NAFL, "No Path Inflam 1", "No Path")) %>%
  #mutate(NAFL = str_replace(NAFL, "NAFL", "NAFL+NASH")) %>%
  #mutate(NAFL = str_replace(NAFL, "NASH", "NAFL")) %>%
  mutate(Histo_STEATOSIS=str_replace(Histo_STEATOSIS, "3","2")) %>% 
  mutate(Histo_FIBROSIS = str_replace(Histo_FIBROSIS,"1a","1")) %>%
  mutate(Histo_FIBROSIS = str_replace(Histo_FIBROSIS,"1b","1")) %>%
  mutate(Histo_FIBROSIS = str_replace(Histo_FIBROSIS,"1c","1")) %>%
  filter(! is.na(Histo_STEATOSIS)) -> CD.2


CD.2 %>% dplyr::select(Age, Histo_STEATOSIS, Gender, `Weight(kg)`, `BMI(kg/m2)`,
                       MetSscore, MetS_Trig, MetS_BP, MetS_BSL) -> CD.2.anthropometricParm

CD.2 %>% dplyr::select(Histo_STEATOSIS, ALT, AST, `ALT/AST`, GGT, ALP, TotalBili, Albumin, Ferritin,
                       Urea,SerumCreatine, eGFR ) -> CD.2.liver.renal

CD.2 %>% dplyr::select(Histo_STEATOSIS, TotalChol, HDL, LDL, `NonHDL-Chol`, Trig#, Chol_Rx_name
                       ) -> CD.2.lipid.profile

CD.2.anthropometricParm %>% dplyr::select(`BMI(kg/m2)`) %>% na.omit() %>% 
  summarize(mean=mean(`BMI(kg/m2)`), se=std.error(`BMI(kg/m2)`))

shapiro.test()


###################

CD.2 %>% dplyr::select(Histo_STEATOSIS, Diabetes, T2DM, Previous.gestational.DM, `HbA1c(%)`,
                       IFG, FBG, Cpeptide, Metformin_treatment,
                       Insulin, `HOMA2%B`, `HOMA2%S`, HOMA2IR ) %>%  mutate(`HOMA2%B` = as.numeric(`HOMA2%B`)) -> 
  CD.2.glucose.para

CD.2.glucose.para %>% group_by(Histo_STEATOSIS) %>% mutate(`HOMA2%B` = as.numeric(`HOMA2%B`)) %>%
  summarise_at(vars(`HOMA2%B`,`HOMA2%S`), list(name=mean, sd),na.rm=TRUE)
  

help("compareGroups")


CD.2 %>% dplyr::select(Histo_STEATOSIS, NAFL,  HISTO_NAS,
                       Histo_INFLAM, Histo_BALLOON, Histo_FIBROSIS, HISTO_NAS,
                       LiverFibrosisgrouping,
                       LiverGrossgrouping) %>% 
  mutate_if(is.numeric, as.factor) %>% 
  mutate(HISTO_NAS = as.numeric(HISTO_NAS)) %>%
  mutate(LiverGrossgrouping=as.numeric(LiverGrossgrouping))-> CD.2.NAFLD

CD.2.NAFLD

meth = 4
createTable( compareGroups(Histo_STEATOSIS ~ Age + Gender +  `Weight(kg)` + 
                             `BMI(kg/m2)` + MetSscore + MetS_Trig + MetS_BP + MetS_BSL, 
                           data = CD.2.anthropometricParm, 
                           method = meth), hide = c(Gender = "Male")) -> CD.2.anthropometricParm.T

CD.2.anthropometricParm.T

createTable( compareGroups(Histo_STEATOSIS ~ ALT + AST + #`ALT/AST` + 
                             GGT + ALP + TotalBili + Albumin + Ferritin +
                             Urea + SerumCreatine + eGFR, 
                           data = CD.2.liver.renal, 
                           method = meth)) -> CD.2.liver.renal.T

CD.2.liver.renal.T

createTable( compareGroups(Histo_STEATOSIS ~ Diabetes + T2DM + `HbA1c(%)` +
                             IFG + FBG + Cpeptide + Metformin_treatment +
                             Insulin + `HOMA2%B` + `HOMA2%S` + HOMA2IR, 
                           data = CD.2.glucose.para, 
                           method = meth)) -> CD.2.glucose.para.T

CD.2.glucose.para.T



compareGroups(Histo_STEATOSIS ~ Diabetes + T2DM + `HbA1c(%)` +
                IFG + FBG + Cpeptide + Metformin_treatment +
                Insulin + `HOMA2%B` + `HOMA2%S`, 
              data = CD.2.glucose.para, subset = Diabetes == "1",
              method = meth) %>% createTable()


createTable( compareGroups(Histo_STEATOSIS ~ NAFL + HISTO_NAS + 
                             Histo_INFLAM + Histo_BALLOON + Histo_FIBROSIS +
                              #LiverFibrosisgrouping +
                             LiverGrossgrouping, 
                           data= CD.2.NAFLD,
                           method = meth)) -> CD.2.NAFLD.T

CD.2.NAFLD %>% mutate(HISTO_NAS = as.numeric(HISTO_NAS)) %>% group_by(Histo_STEATOSIS) %>%
  summarize(mean=mean(HISTO_NAS), sd = sd(HISTO_NAS))

createTable( compareGroups( Histo_STEATOSIS ~  TotalChol + HDL + `NonHDL-Chol` + LDL +
                             + Trig,
                           data= CD.2.lipid.profile,
                           method = meth)) -> CD.2.lipid.profile.T


rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
)

rbind("NAFLD severity" = CD.2.NAFLD.T,
      "Anthropometric parameters" = CD.2.anthropometricParm.T, 
      "Liver enzyme and renal function" = CD.2.liver.renal.T,
      "Lipid profiles" = CD.2.lipid.profile.T,
      "Dibates-related factors" = CD.2.glucose.para.T
) -> restab
 
export2csv(restab, file='table2.csv')

library(lme4)


CD.2 %>% colnames()
#lmer(Age ~ Histo_FIBROSIS + (1|ID), data=CD.2)

model_lm <- lm(ALT ~ Histo_FIBROSIS, data=CD.2)
#model_lmm <- lmer(ALT ~ Histo_FIBROSIS + (1 | ID), data=CD.2)
summary(model_lmm)


####################
library(broom)


###########
clean_background.strong <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("black"),
                          axis.text = element_text(size = 14, color = "black"),
                          axis.title = element_text(color = "black"),
                          axis.text.x = element_text(size=14, face="bold"),
                          axis.text.y = element_text(size=14, face="bold"),
                          legend.text = element_text(size = 14),
                          legend.key = element_rect("white"))
pal <- c("lightsalmon1", "gold1", "palegreen4")

CD.2.liver.renal %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(shapiro.test(x= .$value))) 

CD.2.NAFLD %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS))) 

CD.2.anthropometricParm %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS))) 

CD.2.liver.renal %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS))) 

CD.2.lipid.profile %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS))) 

CD.2.glucose.para %>% gather(key, value, -Histo_STEATOSIS) %>% group_by(key) %>%
  do(tidy(kruskal.test(x= .$value, g = .$Histo_STEATOSIS))) 

help(mutate_at)
CD.2.liver.renal %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
  dplyr::select(Histo_STEATOSIS, ALT) %>%
  group_by(Histo_STEATOSIS) %>% na.omit() %>% summarise_at(vars(ALT), list(mean=mean,sd=sd)) %>%
  ggplot(aes(x=Histo_STEATOSIS, y=mean, group = 1), log="y")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) +
  geom_line() +
  geom_point(size=8) +
  clean_background +
  scale_fill_nejm() +
  scale_color_nejm() 


library(superb)

CD.1 %>% dplyr::select(Histo_STEATOSIS, ALT) %>% str()

kruskal.test(ALT ~ Histo_STEATOSIS, data = CD.1)
dunnTest(ALT ~ Histo_STEATOSIS, data = CD.2,method="bonferroni")

m1 <- glm(ALT ~ Histo_STEATOSIS, data = CD.2)

 
 CD.2.liver.renal %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
   dplyr::select(Histo_STEATOSIS, ALT) %>% na.omit() %>%
   group_by(Histo_STEATOSIS) -> df_ALT
 
 CD.2.liver.renal %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
   dplyr::select(Histo_STEATOSIS, AST) %>% na.omit() %>%
   group_by(Histo_STEATOSIS) -> df_AST
 
 CD.2.liver.renal %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
   dplyr::select(Histo_STEATOSIS, GGT) %>% na.omit() %>%
   group_by(Histo_STEATOSIS) -> df_GGT
 
 CD.2.liver.renal %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
   dplyr::select(Histo_STEATOSIS, TotalBili) %>% na.omit() %>%
   group_by(Histo_STEATOSIS) -> df_TotalBili
 
 CD.2.liver.renal %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
   dplyr::select(Histo_STEATOSIS, Ferritin) %>% na.omit() %>%
   group_by(Histo_STEATOSIS) -> df_Ferritin
 
 my_comparisons <- list( c("0", "1"), c("1", "2"), c("0", "2") )
 my_comparisons.1 <- list(  c("1", "2"), c("0", "2") )  
 my_comparisons.2 <- list(  c("0", "1"), c("0", "2") ) 
 my_comparisons.3 <- list( c("0", "2") )  
 
pdf("try.pdf")

p_ALT <- ggboxplot(df_ALT, x="Histo_STEATOSIS", y="ALT", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
            alpha=0.2, outlier.shape = NA, add = "jitter",
            add.params = list(alpha = 0.5),width=0.8) +
    stat_compare_means(label.y=150)+
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                       method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                       label.y = c(60,80,100))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(0,150) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "ALT (U/L)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )

p_ALT

p_AST <- ggboxplot(df_AST, x="Histo_STEATOSIS", y="AST", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y=100)+
  stat_compare_means(comparisons = my_comparisons.1, label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     label.y = c(60,70))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(0,100) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "AST (U/L)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )


p_AST

p_GGT <- ggboxplot(df_GGT, x="Histo_STEATOSIS", y="GGT", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y=100)+
  stat_compare_means(comparisons = my_comparisons.1, label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     label.y = c(50,60))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(0,100) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "GGT (U/L)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )


p_GGT

p_TBIL <- ggboxplot(df_TotalBili, x="Histo_STEATOSIS", y="TotalBili", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y=40)+
  stat_compare_means(comparisons = my_comparisons.2, label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     label.y = c(15,20))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(0,40) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "TBIL (\u00B5 mol/L)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )

p_TBIL

p_Ferritin <- ggboxplot(df_Ferritin, x="Histo_STEATOSIS", y="Ferritin", fill="Histo_STEATOSIS", 
          color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y=600)+
  stat_compare_means(comparisons = my_comparisons.3, label = "p.signif", 
                     method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
                     label.y = c(500))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(0,600) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Ferritin(ng/mL)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )

pdf("try.pdf")
p_Ferritin
dev.off()

pdf("LiverEnzyme.pdf", width=15, height=4.5)
ggarrange(p_ALT, p_AST, p_GGT, p_TBIL, p_Ferritin,
          labels=c("C","D","E","F","G"),
          nrow=1)
dev.off()


ClinicalData %>% dplyr::select(Age) %>% na.omit() %>% summarize(mean=mean(Age), sd=sd(Age))

ClinicalData %>% dplyr::select(`BMI(kg/m2)`) %>% na.omit() %>% summarize(mean=sd(`BMI(kg/m2)`))

CD.2 %>% dplyr::select(TotalChol, Histo_STEATOSIS) %>% mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) -> df_TotalChol

ggboxplot(df_TotalChol, x="Histo_STEATOSIS", y="TotalChol", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y=7)+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
  #                   method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
  #                   label.y = c(5,6,7))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(1,8) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Total Cholesterol (mmol/L)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )


CD.2 %>% dplyr::select(TotalChol, HDL, `NonHDL-Chol`, Trig, Histo_STEATOSIS, LDL) %>% 
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
  gather(key="Lipids", value="value", -Histo_STEATOSIS) %>%
  mutate(value=value, Lipids = factor(Lipids, levels=c("TotalChol", "NonHDL-Chol","LDL", "Trig","HDL"))) %>% 
  group_by(Histo_STEATOSIS, Lipids) %>% na.omit() %>%
  summarize(mean=mean(value), se=std.error(value)) -> df_lipid
df_lipid

dev.off()

pdf("try.pdf")

df_lipid %>% ggplot(aes(x=Histo_STEATOSIS, y=mean, color=Lipids)) +
  geom_line(aes(group=Lipids, color=Lipids), size=2) +
  geom_point(aes(group=Lipids), shape=19, size=3, stroke=2, position=position_dodge(0.05)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, linewidth=1, position=position_dodge(0.05)) +
  #scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  annotate(geom="text",x=2.5,y=3.9, label="Total cholesterol", size=4.5) +
  annotate(geom="text",x=2.5,y=1.44, label="NonHDL-cholesterol", size=4.5) +
  annotate(geom="text",x=2.8,y=0.7, label="LDL", size=4.5) +
  annotate(geom="text",x=2.3,y=0.306, label="Triglycerides", size=4.5) +
  annotate(geom="text",x=2.3,y=0.306, label="HDL", size=4.5) +
  #annotate(geom="text",x=3.3,y=3.56, label="***", size=4.5) +
  annotate(geom="text",x=3.3,y=3.16, label="*", size=4.5) +
  annotate(geom="text",x=3.3,y=1.68, label="*", size=4.5) +
  annotate(geom="text",x=3.3,y=0.969, label="*", size=4.5) +
  #ylim(-2,2) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Lipid profile (mmol/L)") + # guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
       # legend.position = right, 
        legend.background = element_rect(fill = alpha("darkkhaki",.2))
  )+scale_color_discrete(name = "", labels = c("Total cholesterol","NonHDL-cholesterol","LDL","Triglycerides","HDL"))

df_lipid

dev.off()

CD.2 %>% dplyr::select(Histo_STEATOSIS,  HISTO_NAS,
                       Histo_INFLAM, Histo_BALLOON, Histo_FIBROSIS, HISTO_NAS,
                       LiverFibrosisgrouping,
                       LiverGrossgrouping) %>% 
  mutate(LiverGrossgrouping=as.numeric(LiverGrossgrouping))-> df.NAFLD


df.NAFLD %>% 
  mutate(Histo_STEATOSIS = as.factor(Histo_STEATOSIS)) %>% 
  gather(key="NAFLD", value="value", -c(Histo_STEATOSIS, LiverGrossgrouping, LiverFibrosisgrouping)) %>%
  #mutate(value=value, NAFLD = factor(NAFLD, levels=c("HISTO_NAS", "Histo","Trig","HDL", "LDL"))) %>% 
  group_by(Histo_STEATOSIS, NAFLD) %>% mutate(value=as.numeric(value)) %>%
  summarize(mean=mean(value), se=std.error(value), SD=sd(value)) -> df_NAFLD.1

df_NAFLD.1
help("annotate")
df_NAFLD.1 %>% ggplot(aes(x=Histo_STEATOSIS, y=mean, color=NAFLD)) +
  geom_line(aes(group=NAFLD, color=NAFLD), size=2) +
  geom_point(aes(group=NAFLD), shape=19, size=3, stroke=2, position=position_dodge(0.05)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, linewidth=1, position=position_dodge(0.05)) +
  #scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  annotate(geom="text",x=2.7,y=3.56, label="NAS", size=4.5) +
  annotate(geom="text",x=2.5,y=1.44, label="Fibrosis", size=4.5) +
  annotate(geom="text",x=2.8,y=0.7, label="Inflammation", size=4.5) +
  annotate(geom="text",x=2.3,y=0.306, label="Balloning", size=4.5) +
  annotate(geom="text",x=3.3,y=3.56, label="***", size=4.5) +
  annotate(geom="text",x=3.3,y=1.44, label="***", size=4.5) +
  annotate(geom="text",x=3.3,y=0.732, label="***", size=4.5) +
  annotate(geom="text",x=3.3,y=0.306, label="***", size=4.5) +
  #geom_text(aes(label = names(NAFLD), x = 3, colour = names(NAFLD), y = c(NAFLD), hjust = -.02))+
  #xlim(0,5) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Grade") + # guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        # legend.position = right, 
        legend.background = element_rect(fill = alpha("darkkhaki",.2))
  )#+scale_color_discrete(name = "", labels = c("Ballooning","Fibrosis","Inflammation","NAS"))


help(cor)
help("cor.test")
cor.test( as.numeric(df.NAFLD$HISTO_NAS), as.numeric(df.NAFLD$Histo_STEATOSIS),method="kendall")
cor.test( as.numeric(df.NAFLD$Histo_FIBROSIS), as.numeric(df.NAFLD$Histo_STEATOSIS),method="kendall")

ClinicalData %>% dplyr::select(HISTO_NAS, Histo_STEATOSIS) %>% na.omit() %>%
  mutate(HISTO_NAS = as.numeric(HISTO_NAS)) %>% group_by(Histo_STEATOSIS) %>%
  summarize(mean=mean(HISTO_NAS), sd = sd(HISTO_NAS))

ggboxplot(df_TotalChol, x="Histo_STEATOSIS", y="TotalChol", fill="Histo_STEATOSIS", color="Histo_STEATOSIS",
          alpha=0.2, outlier.shape = NA, add = "jitter",
          add.params = list(alpha = 0.5),width=0.8) +
  stat_compare_means(label.y=7)+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
  #                   method = "wilcox.test", paired = F, hide.ns = T, show.legend = F, size = 8,
  #                   label.y = c(5,6,7))+ # Add pairwise comparisons p-value
  scale_fill_nejm() +
  scale_color_nejm() +
  clean_background.strong +
  ylim(1,8) +
  scale_x_discrete(labels=c("S0","S1","S2+S3")) +
  scale_fill_discrete(name = "Steatosis", labels  = c("S0 (n = 37)", "S1 (n = 41)", "S2+S3: (n = 36)")) +
  labs(x = "Steatosis level",
       y = "Total Cholesterol (mmol/L)") +  guides(color=FALSE, alpha = FALSE)+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 16),
        legend.position = "bottom"
  )


dev.off()


pdf("try.pdf")


dev.off()
