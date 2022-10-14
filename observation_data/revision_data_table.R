setwd("c:/users/jwambaug/git/EcoSEEM-Consensus-Model-for-Chemicals-in-Surface-Water/observation_data/")
dir()
final.chems <- read.csv("all_chem_res_2022_10_13.csv")
dim(final.chems)
4
423+172
citation()
R.version
ls()
head(final.chems)
km <- read.excel("all_chem_res_km.xlsx")
library(readXL)
library(readxl)
km <- read.excel("all_chem_res_km.xlsx")
km <- read_excel("all_chem_res_km.xlsx")
head(km)
km <- as.data.frame(km)
dim(final.chems)
all.data <- merge(final.chems,km, all.x=TRUE,all.y=TRUE)
colnames(all.data)
head(all.data)
head(all.data)[,1:5]
head(all.data)[,6:10]
colnames(all.data)
colnames(km)
colnames(all.data)
all.data$Conc.Median.ugpL <- all.data$km_50th
all.data$Conc.95th.ugpL <- all.data$km_95th
all.data$Censor.Median.ugpL <- all.data$bloq_med_cen_level_dslv
all.data$Censor.Max.ugpL <- all.data$bloq_max_cen_level_dslv
colnames(all.data)
colnames(all.data)[1:10]
colnames(all.data)[1:15]
colnames(all.data)[1:10]
all.data <- all.data[,c(1:9,128:131,10:127)]
colnames(all.data)[1:10]
head(all.data)[,1:15]
all.data[,9]
all.data[,10]
all.data[,11]
all.data[1:10,1:15]
sum(!is.na(all.data[,10]))
sum(!is.na(all.data[,11]))
sum(!is.na(all.data[,12]))
sum(!is.na(all.data[,13]))
write.csv(all.data,row.names=FALSE,file="all_chem_results.csv")