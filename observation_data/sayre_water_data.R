library("EnvStats")
library('tidyr')
library('dplyr')
library('ggplot2')
library('readxl')
library('NADA2') #! only for cboxplot
library('forcats')
library('viridis')

## A. From Method 1. Exclude chemicals with fewer than 50 uncensored (detected) concentration values

df_d <- read.csv("C:\\Users\\rrace\\Documents\\ECOSEEM\\smaller_df_21july2020.csv")

df_d <- within(df_d, bloq[res_val == 0] <- 1)
df_d$obs <- with(df_d, ifelse(bloq==1, loq_val, res_val))
df_d$ln_obs <- log(df_d$obs)
DF <- df_d[!is.na(df_d$ln_obs),]
DF$cen <- with(DF, ifelse(bloq==1, TRUE, FALSE))

#no_obs <- DF %>%
#  group_by(chem_id, bloq) %>%
#  summarise(mean = mean(obs), n = n()) %>% filter(n()==1) %>% filter(bloq == 1)
#chem_ids_no_obs <- unique(no_obs$chem_id)
#remove(no_obs)

gt50_obs <- DF[which(DF$bloq == 0),] %>%
  group_by(chem_id) %>%
  filter(n() > 50)
chem_ids <- unique(gt50_obs$chem_id)
remove(gt50_obs)

gt50_recs <- DF[which(!DF$chem_id %in% chem_ids),] %>%
  group_by(chem_id) %>%
  filter(n() > 50)
bloq_chem_ids <- unique(gt50_recs$chem_id)
#

DF_smaller <- DF[which(DF$chem_id %in% chem_ids),]
DF_smaller$year = substr(DF_smaller$ActivityStartDate, start = 1, stop = 4)
DF_w_HUC <- merge(DF_smaller,station_unique,by='MonitoringLocationIdentifier', all.x=TRUE)
# saved so I don't have to redo this slow step
#DF_w_HUC <- read.csv("C:\\Users\\rrace\\Documents\\ECOSEEM\\DF_smaller_w_HUC_15may2021.csv")
#! can remove other unnecessary fields
DF_w_HUC <- subset(DF_w_HUC, select = -c(OrganizationFormalName,MonitoringLocationDescriptionText))
DF_w_HUC$HUC02 <- substr(DF_w_HUC$HUCEightDigitCode, start = 1, stop = 2)
chem_ids <- unique(DF_w_HUC$chem_id)

plot_it <- DF_w_HUC %>% 
  group_by(chem_id) %>% 
  summarise(n=n_distinct(HUC02))

plot_it <- arrange(plot_it,n)

plot_it %>% 
  group_by(n) %>% 
  summarise(n_chem=n_distinct(chem_id)) %>% 
  ggplot(., aes(x=n, y=n_chem)) +
  geom_bar(stat='identity') +
  theme_bw() +
  xlab('Count of HUCs sampled per chemical') +
  ylab('Count of chemicals with at least 50 quantified values')

#cboxplot(this_res$obs, this_res$cen, this_res$season, show=TRUE, LOG=TRUE)

## B. From Method 2. Compare concentration ranges of technical min and reporting min nondetect samples (with detects)
#https://search.r-project.org/CRAN/refmans/EnvStats/html/twoSampleLinearRankTestCensored.html
#null hypothesis is that both are the same

limit_comparisons = list()
for (i in 1:length(chem_ids)){
  chem = chem_ids[i]
  this_res <- DF_w_HUC[which(DF_w_HUC$chem_id == chem),] 
  
  tryCatch({
    x <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$loq_type == 1 | this_res$cen == FALSE),]$obs, #tech min
                                         x.censored = this_res[which(this_res$loq_type == 1 | this_res$cen == FALSE),]$cen, 
                                         y = this_res[which(this_res$loq_type == 3 | this_res$cen == FALSE),]$obs, # reporting min
                                         y.censored = this_res[which(this_res$loq_type == 3 | this_res$cen == FALSE),]$cen,
                                         test="peto-peto"
    )
    if(x$sample.size[[1]] > 0 & x$sample.size[[2]] > 0){
      limit_comparisons[[i]] <- c(chem,x$p.value)      
    }

  },
  error=function(cond) {})
}

df <- as.data.frame(as.matrix(limit_comparisons))
wider_limit_df <- df %>% unnest_wider(V1)
names(wider_limit_df)[1] <- "chem_id"
names(wider_limit_df)[2] <- "limit_pval"
n_significant_per_limit <- sum(wider_limit_df$limit_pval < 0.05 & wider_limit_df$limit_pval > 0, na.rm=TRUE)
n_notsignificant_per_limit <- sum(wider_limit_df$limit_pval > 0.05, na.rm=TRUE)
prop_chem_diff_limit = round(n_significant_per_limit/(n_notsignificant_per_limit+n_significant_per_limit),3)
paste("Result B: For", prop_chem_diff_limit*100, "% of chemicals, per-chemical concentration ranges were shown to differ by limit type.", sep = " ")

# how much does excluding the reporting min decrease the chemical space?

## C. From Method 2. Compare concentration ranges of bulk and dissolved samples

phase_comparisons = list()
for (i in 1:length(chem_ids)){
  chem=chem_ids[i]
  this_res <- DF_w_HUC[which(DF_w_HUC$chem_id == chem),] 
  
  tryCatch({
    x <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$phase == 2),]$obs, #dissolved
                                         x.censored = this_res[which(this_res$phase == 2),]$cen, 
                                         y = this_res[which(this_res$phase == 3),]$obs, #bulk
                                         y.censored = this_res[which(this_res$phase == 3),]$cen,
                                         test="peto-peto")
    if(x$sample.size[[1]] > 0 & x$sample.size[[2]] > 0){
      phase_comparisons[[i]] <- c(chem,x$p.value)      
    }
  },
  error=function(cond) {})
}     

df <- as.data.frame(as.matrix(phase_comparisons))
wider_phase_df <- df %>% unnest_wider(V1)
names(wider_phase_df)[1] <- "chem_id"
names(wider_phase_df)[2] <- "phase_pval"
n_significant_per_phase <- sum(wider_phase_df$phase_pval < 0.05 & wider_phase_df$phase_pval > 0, na.rm=TRUE)
n_notsignificant_per_phase <- sum(wider_phase_df$phase_pval > 0.05, na.rm=TRUE)
prop_chem_diff_phase = round(n_significant_per_phase/(n_notsignificant_per_phase+n_significant_per_phase),3)
paste("Result C: For", prop_chem_same_phase*100, "% of chemicals, per-chemical concentration ranges were shown to differ by phase.", sep = " ")

# Perform a logistic regression to evaluate if physicochemical properties can 
#   contribute to an explanation of chemicals for which concentrations differ by phase

# Add Bool column to indicate Peto-Peto test results
wider_phase_df$phase_diff <- ifelse(wider_phase_df$phase_pval < 0.05 & wider_phase_df$phase_pval > 0,1,ifelse(wider_phase_df$phase_pval > 0.05,0,NA))

# Join physicochemical properties
pcp <- read.csv("C:\\Users\\rrace\\Documents\\ECOSEEM\\dtxsids_physchemprop_24jan2021_dedup.csv")
#! it is my coding style to write over temporary tables with a generic identifier
df <- merge(x = wider_phase_df, y = pcp, by = "chem_id", all.x=TRUE)

# Check p-values of logistic regression on test results vs. water-relevant physicochemical properties
summary(glm(phase_diff ~ as.numeric(OCTANOL_WATER_PARTITION_LOGP_OPERA_PRED) +
              as.numeric(VAPOR_PRESSURE_MMHG_OPERA_PRED) +
              as.numeric(WATER_SOLUBILITY_MOL.L_OPERA_PRED) +
              as.numeric(log.KAW),
            data = df[!is.na(df$phase_diff), ], family = binomial))
# The result: none of the p-values were compelling

#! Plot conc differences between bulk and dissolved again

# Because differences in concentration were statistically different, separate 
#   dissolved for better comparison with toxicity values
#! can remove this step
dslv_DF <- DF_w_HUC[which(DF_w_HUC$phase == 2),]
bulk_DF <- DF_w_HUC[which(DF_w_HUC$phase == 3),]

#remove(DF_w_HUC)

# Find new subset of chemicals with gt50 obs
gt50_obs <- dslv_DF[which(dslv_DF$bloq == 0),] %>%
  group_by(chem_id) %>%
  filter(n() > 50)
new_chem_ids <- unique(gt50_obs$chem_id)
remove(gt50_obs)

# Find new subset of chemicals with gt50 obs
gt50_bulk_obs <- bulk_DF[which(bulk_DF$bloq == 0),] %>%
  group_by(chem_id) %>%
  filter(n() > 50)
b_chem_ids <- unique(gt50_obs$chem_id)

dslv_DF_smaller <- dslv_DF[which(dslv_DF$chem_id %in% new_chem_ids),]
bulk_DF_smaller <- bulk_DF[which(bulk_DF$chem_id %in% b_chem_ids),]
remove(dslv_DF)

gt50_dslv_recs <- gt50_recs[which(gt50_recs$phase == 2),] %>%
  group_by(chem_id) %>%
  filter(n() > 50)
bloq_dslv_chem_ids <- unique(gt50_dslv_recs$chem_id)

gt50_bulk_recs <- gt50_recs[which(gt50_recs$phase == 3),] %>%
  group_by(chem_id) %>%
  filter(n() > 50)
bloq_bulk_chem_ids <- unique(gt50_bulk_recs$chem_id)

## D. From Method 2. Compare concentration ranges of samples from each season pair

season_comparisons = list()
for (i in 1:length(new_chem_ids)){
  chem = new_chem_ids[i]
  this_res <- dslv_DF_smaller[which(dslv_DF_smaller$chem_id == chem),] 
  
  tryCatch({
    a <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$season == 0),]$obs,
                                         x.censored = this_res[which(this_res$season == 0),]$cen, 
                                         y = this_res[which(this_res$season == 1),]$obs,
                                         y.censored = this_res[which(this_res$season == 1),]$cen,
                                         test="peto-peto"
    )
    b <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$season == 0),]$obs,
                                         x.censored = this_res[which(this_res$season == 0),]$cen, 
                                         y = this_res[which(this_res$season == 2),]$obs,
                                         y.censored = this_res[which(this_res$season == 2),]$cen,
                                         test="peto-peto"
    )
    c <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$season == 0),]$obs,
                                         x.censored = this_res[which(this_res$season == 0),]$cen, 
                                         y = this_res[which(this_res$season == 3),]$obs,
                                         y.censored = this_res[which(this_res$season == 3),]$cen,
                                         test="peto-peto"
    )
    d <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$season == 1),]$obs,
                                         x.censored = this_res[which(this_res$season == 1),]$cen, 
                                         y = this_res[which(this_res$season == 2),]$obs,
                                         y.censored = this_res[which(this_res$season == 2),]$cen,
                                         test="peto-peto"
    )
    e <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$season == 1),]$obs,
                                         x.censored = this_res[which(this_res$season == 1),]$cen, 
                                         y = this_res[which(this_res$season == 3),]$obs,
                                         y.censored = this_res[which(this_res$season == 3),]$cen,
                                         test="peto-peto"
    )
    f <- twoSampleLinearRankTestCensored(x = this_res[which(this_res$season == 2),]$obs,
                                         x.censored = this_res[which(this_res$season == 2),]$cen, 
                                         y = this_res[which(this_res$season == 3),]$obs,
                                         y.censored = this_res[which(this_res$season == 3),]$cen,
                                         test="peto-peto"
    )
    if(a$sample.size[[1]] > 0 & a$sample.size[[2]] > 0 & b$sample.size[[1]] > 0 & b$sample.size[[2]] > 0 & c$sample.size[[1]] > 0 & c$sample.size[[2]] > 0 & d$sample.size[[1]] > 0 & d$sample.size[[2]] > 0 & e$sample.size[[1]] > 0 & e$sample.size[[2]] > 0 & f$sample.size[[1]] > 0 & f$sample.size[[2]] > 0){
        season_comparisons[[i]] <- c(chem,a$p.value,b$p.value,c$p.value,d$p.value,e$p.value,f$p.value)
    }
  },
  error=function(cond) {})
}

df <- as.data.frame(as.matrix(season_comparisons))
wider_season_df <- df %>% unnest_wider(V1)
names(wider_season_df)[1] <- "chem_id"
names(wider_season_df)[2] <- "season01_pval"
names(wider_season_df)[3] <- "season02_pval"
names(wider_season_df)[4] <- "season03_pval"
names(wider_season_df)[5] <- "season12_pval"
names(wider_season_df)[6] <- "season13_pval"
names(wider_season_df)[7] <- "season23_pval"
wider_season_df <- wider_season_df %>%
  mutate(n_significant_per_chem = rowSums(select(wider_season_df, `season01_pval`:`season23_pval`) < 0.05))
n_significant_per_season01 <- sum(wider_season_df$season01_pval < 0.05 & wider_season_df$season01_pval > 0, na.rm=TRUE)
n_significant_per_season02 <- sum(wider_season_df$season02_pval < 0.05 & wider_season_df$season01_pval > 0, na.rm=TRUE)
n_significant_per_season03 <- sum(wider_season_df$season03_pval < 0.05 & wider_season_df$season01_pval > 0, na.rm=TRUE)
n_significant_per_season12 <- sum(wider_season_df$season12_pval < 0.05 & wider_season_df$season01_pval > 0, na.rm=TRUE)
n_significant_per_season13 <- sum(wider_season_df$season13_pval < 0.05 & wider_season_df$season01_pval > 0, na.rm=TRUE)
n_significant_per_season23 <- sum(wider_season_df$season23_pval < 0.05 & wider_season_df$season01_pval > 0, na.rm=TRUE)
n_notsignificant_per_season01 <- sum(wider_season_df$season01_pval > 0.05, na.rm=TRUE)
n_notsignificant_per_season02 <- sum(wider_season_df$season02_pval > 0.05, na.rm=TRUE)
n_notsignificant_per_season03 <- sum(wider_season_df$season03_pval > 0.05, na.rm=TRUE)
n_notsignificant_per_season12 <- sum(wider_season_df$season12_pval > 0.05, na.rm=TRUE)
n_notsignificant_per_season13 <- sum(wider_season_df$season13_pval > 0.05, na.rm=TRUE)
n_notsignificant_per_season23 <- sum(wider_season_df$season23_pval > 0.05, na.rm=TRUE)
prop_chem_diff_season01 = round(n_significant_per_season01/(n_notsignificant_per_season01+n_significant_per_season01),3)
prop_chem_diff_season02 = round(n_significant_per_season02/(n_notsignificant_per_season02+n_significant_per_season02),3)
prop_chem_diff_season03 = round(n_significant_per_season03/(n_notsignificant_per_season03+n_significant_per_season03),3)
prop_chem_diff_season12 = round(n_significant_per_season12/(n_notsignificant_per_season12+n_significant_per_season12),3)
prop_chem_diff_season13 = round(n_significant_per_season13/(n_notsignificant_per_season13+n_significant_per_season13),3)
prop_chem_diff_season23 = round(n_significant_per_season23/(n_notsignificant_per_season23+n_significant_per_season23),3)
df <- merge(x = wider_season_df, y = pcp, by = "chem_id", all.x=TRUE)
df$EPAPCS[df$EPAPCS == 0] <- "Non-pesticide"
df$EPAPCS[df$EPAPCS == 1] <- "Pesticide"

# figure 3
ggplot(df[!is.na(df$EPAPCS), ],
       aes(x = n_significant_per_chem, fill=EPAPCS)) + 
  geom_histogram(position="stack", alpha=0.5) +
  facet_grid(EPAPCS ~ .) + 
  theme_bw() + 
  xlab("Season pair differences") + 
  ylab("Count of chemicals") +
  scale_fill_manual(values = c("#34f900", "#644cff"))

## E. From Method 3. Calculate lognormal Maximum Likelihood Estimate parameters

mle_envstats = list()

for (i in 1:length(new_chem_ids)) {
  chem = new_chem_ids[i]
  this_res <- dslv_DF_smaller[which(dslv_DF_smaller$chem_id == chem),]  
  tryCatch({
    env_mle_05 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.05,
                                  #ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    env_mle_25 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.25,
                                  #ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    env_mle_50 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.5,
                                  ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    env_mle_75 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.75,
                                  #ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    env_mle_95 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.95,
                                  #ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    env_mle_99 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.99,
                                  #ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.99)
    n_cen = length(env_mle_50$censoring.levels)
    mle_envstats[[i]] <- c(chem,
                           env_mle_05$sample.size,
                           env_mle_05$percent.censored,
                           env_mle_05$quantiles[[1]], 
                           env_mle_25$quantiles[[1]], 
                           env_mle_50$quantiles[[1]], 
                           env_mle_50$interval$limits[[1]], #lcl
                           env_mle_50$interval$limits[[2]], #ucl
                           env_mle_75$quantiles[[1]], 
                           env_mle_95$quantiles[[1]],
                           env_mle_99$quantiles[[1]],
                           env_mle_50$parameters[[2]],
                           env_mle_50$censoring.levels[1],
                           env_mle_50$censoring.levels[n_cen]
    )
    
  },
  error=function(cond) {})
  
} 

df <- as.data.frame(as.matrix(mle_envstats))
wider_mle_df <- df %>% unnest_wider(V1)
names(wider_mle_df)[1] <- "chem_id"
names(wider_mle_df)[2] <- "sample_size_dslv"
names(wider_mle_df)[3] <- "pct_cen_dslv"
names(wider_mle_df)[4] <- "mle_5th_dslv"
names(wider_mle_df)[5] <- "mle_25th_dslv"
names(wider_mle_df)[6] <- "mle_50th_dslv"
names(wider_mle_df)[7] <- "mle_50th_lcl_dslv"
names(wider_mle_df)[8] <- "mle_50th_ucl_dslv"
names(wider_mle_df)[9] <- "mle_75th_dslv"
names(wider_mle_df)[10] <- "mle_95th_dslv"
names(wider_mle_df)[11] <- "mle_99th_dslv"
names(wider_mle_df)[12] <- "mle_sdlog_dslv"
names(wider_mle_df)[13] <- "min_cen_level_dslv"
names(wider_mle_df)[14] <- "max_cen_level_dslv"

write.csv(wider_mle_df,"nov1dsvl.csv")

#! create bulk df smaller
mle_envstats_b = list()

for (i in 1:length(b_chem_ids)) {
  chem = b_chem_ids[i]
  this_res <- bulk_DF_smaller[which(bulk_DF_smaller$chem_id == chem),]  
  tryCatch({
    env_mle_50 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.5,
                                  ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    mle_envstats_b[[i]] <- c(chem,
                           env_mle_50$sample.size,
                           env_mle_50$percent.censored,
                           env_mle_50$quantiles[[1]], 
                           env_mle_50$interval$limits[[1]], #lcl
                           env_mle_50$interval$limits[[2]], #ucl
                           env_mle_50$parameters[[2]], #sdlog     
                           env_mle_50$censoring.levels[n_cen]
                          )
    
  },
  error=function(cond) {})
  
} 


df <- as.data.frame(as.matrix(mle_envstats_b))
wider_bulk_mle_df <- df %>% unnest_wider(V1)
names(wider_bulk_mle_df)[1] <- "chem_id"
names(wider_bulk_mle_df)[2] <- "sample_size_bulk"
names(wider_bulk_mle_df)[3] <- "pct_cen_bulk"
names(wider_bulk_mle_df)[4] <- "mle_50th_bulk"
names(wider_bulk_mle_df)[5] <- "mle_50th_lcl_bulk"
names(wider_bulk_mle_df)[6] <- "mle_50th_ucl_bulk"
names(wider_bulk_mle_df)[7] <- "mle_sdlog_bulk"
names(wider_bulk_mle_df)[8] <- "max_cen_level_bulk"

write.csv(wider_bulk_mle_df,"nov1bulk.csv")

mle_envstats_bloq = list()

for (i in 1:length(bloq_dslv_chem_ids)) {
  chem = bloq_dslv_chem_ids[i]
  this_res <- gt50_dslv_recs[which(gt50_dslv_recs$chem_id == chem),]  
  tryCatch({
    env_mle_50 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.5,
                                  ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    n_cen = length(env_mle_50$censoring.levels)
    med_limit <- median(this_res[which(this_res$bloq == 1),]$obs)
    mle_envstats_bloq[[i]] <- c(chem,
                             env_mle_50$sample.size,
                             env_mle_50$percent.censored,
                             env_mle_50$quantiles[[1]], 
                             env_mle_50$interval$limits[[1]], #lcl
                             env_mle_50$interval$limits[[2]], #ucl
                             env_mle_50$parameters[[2]], #sdlog
                             env_mle_50$censoring.levels[n_cen],
                             med_limit
    )
    
  },
  error=function(cond) {})
  
} 

df <- as.data.frame(as.matrix(mle_envstats_bloq))
wider_bloq_mle_df <- df %>% unnest_wider(V1)
names(wider_bloq_mle_df)[1] <- "chem_id"
names(wider_bloq_mle_df)[2] <- "bloq_sample_size_dslv"
names(wider_bloq_mle_df)[3] <- "bloq_pct_cen_dslv"
names(wider_bloq_mle_df)[4] <- "bloq_mle_50th_dslv"
names(wider_bloq_mle_df)[5] <- "bloq_mle_50th_lcl_dslv"
names(wider_bloq_mle_df)[6] <- "bloq_mle_50th_ucl_dslv"
names(wider_bloq_mle_df)[7] <- "bloq_mle_sdlog_dslv"
names(wider_bloq_mle_df)[8] <- "bloq_max_cen_level_dslv"
names(wider_bloq_mle_df)[9] <- "bloq_med_cen_level_dslv"

mle_envstats_bloq_b = list()

for (i in 1:length(bloq_bulk_chem_ids)) {
  chem = bloq_bulk_chem_ids[i]
  this_res <- gt50_bulk_recs[which(gt50_bulk_recs$chem_id == chem),]  
  tryCatch({
    env_mle_50 <- eqlnormCensored(this_res$obs, as.logical(this_res$cen), method = "mle",
                                  censoring.side = "left", 
                                  p = 0.5,
                                  ci = TRUE, ci.method = "exact.for.complete", ci.type = "two-sided",
                                  conf.level = 0.95)
    n_cen = length(env_mle_50$censoring.levels)
    med_limit <- median(this_res[which(this_res$bloq == 1),]$obs)
    mle_envstats_bloq_b[[i]] <- c(chem,
                                env_mle_50$sample.size,
                                env_mle_50$percent.censored,
                                env_mle_50$quantiles[[1]], 
                                env_mle_50$interval$limits[[1]], #lcl
                                env_mle_50$interval$limits[[2]], #ucl
                                env_mle_50$parameters[[2]], #sdlog
                                env_mle_50$censoring.levels[n_cen],
                                med_limit
    )
    
  },
  error=function(cond) {})
  
} 

df <- as.data.frame(as.matrix(mle_envstats_bloq_b))
wider_bloq_bulk_mle_df <- df %>% unnest_wider(V1)
names(wider_bloq_bulk_mle_df)[1] <- "chem_id"
names(wider_bloq_bulk_mle_df)[2] <- "bloq_sample_size_bulk"
names(wider_bloq_bulk_mle_df)[3] <- "bloq_pct_cen_bulk"
names(wider_bloq_bulk_mle_df)[4] <- "bloq_mle_50th_bulk"
names(wider_bloq_bulk_mle_df)[5] <- "bloq_mle_50th_lcl_bulk"
names(wider_bloq_bulk_mle_df)[6] <- "bloq_mle_50th_ucl_bulk"
names(wider_bloq_bulk_mle_df)[7] <- "bloq_mle_sdlog_bulk"
names(wider_bloq_bulk_mle_df)[8] <- "bloq_max_cen_level_bulk"
names(wider_bloq_bulk_mle_df)[9] <- "bloq_med_cen_level_bulk"

## F. From Method 3. Calculate Kaplan-Meier empirical 95th percentile

km_envstats = list()
for (i in 1:length(new_chem_ids)) {
  chem = new_chem_ids[i]
  this_res <- dslv_DF_smaller[which(dslv_DF_smaller$chem_id == chem),] 
  tryCatch({  
    km <- ppointsCensored(this_res$obs, as.logical(this_res$cen), censoring.side = "left")
    idx_km50 = which.min(abs(km$Cumulative.Probabilities - 0.50))
    idx_km95 = which.min(abs(km$Cumulative.Probabilities - 0.95))
    idx_km99 = which.min(abs(km$Cumulative.Probabilities - 0.99))
    km_envstats[[i]] <- c(chem,
                          km$Order.Statistics[idx_km50],
                          km$Order.Statistics[idx_km95],
                          km$Order.Statistics[idx_km99])
  },
  error=function(cond) {})
}

df <- as.data.frame(as.matrix(km_envstats))
wider_km_df <- df %>% unnest_wider(V1)
names(wider_km_df)[1] <- "chem_id"
names(wider_km_df)[2] <- "km_50th"
names(wider_km_df)[3] <- "km_95th"
names(wider_km_df)[4] <- "km_99th"

df <- merge(wider_km_df,wider_mle_df,by = "chem_id", all.x=TRUE)

## G. From Method 4. Compare SSDs to representative values

# G.1. Load SSDs
#! the original publication supplement has been altered to add dtxsids
ssd_values <- read_excel("C:\\Users\\rrace\\Documents\\ECOSEEM\\etc4373-sup-0002-supmat.xlsx", 2)

# Join physicochemical properties
pcp <- read.csv("C:\\Users\\rrace\\Documents\\ECOSEEM\\dtxsids_physchemprop_24jan2021_dedup.csv")
repval_w_pcp <- merge(df,pcp,by = "chem_id", all.x=TRUE)
repval_w_pcp_ssd <- merge(repval_w_pcp,ssd_values,by = "DTXSID", all.x=TRUE)

# Filter out representative values with no SSDs (96 removed)
df <- repval_w_pcp_ssd[!is.na(repval_w_pcp_ssd$`QualityScore-Chronic NOEC`), ]

# Filter out SSDs with lower quality flags (2 removed)
final_df <- df[!grepl("^(2)", df$`QualityScore-Chronic NOEC`),]
df <- final_df[order(final_df$`10LogSSDMedianConcentration(ug/L)-MuChronic NOEC`),] 
df <- df[!is.na(df$chem_id), ]
df <- df[!is.na(df$km_95th), ]
df <- df[order(df$pct_cen),] 
df$idx <- seq.int(nrow(df))

# figure 5
ggplot(data=df) + 
  geom_segment(aes(x = idx, xend = idx, y = km_50th, yend = km_95th),
               size=1, lty="11", color = "#ff8000"
  ) +
  geom_segment(aes(x = idx, xend = idx, y = mle_50th, yend = mle_95th),
               position = position_nudge(x = 0.4),  size=1, lty="31",color = "#1ab3ff") +  #, 
  #geom_point(aes(x = idx, y = min_cen_level), shape = 6, color = 'grey30') +
  geom_point(aes(x = idx, y = max_cen_level), shape = 2, color = 'grey30') +
  #"#0072B2"
  scale_y_log10() + 
  theme_bw() + 
  theme(axis.text.x=element_blank()) + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(x = "Evaluation set chemicals from lowest to highest percent censoring",
       y = "Concentrations (log \u03bcg/L)")


# G. Case 1: The least risky condition (high observed values compared with lowest part of the SSD)
#df$mu <- 10^df$`10LogSSDMedianConcentration(ug/L)-MuChronic NOEC`
#df$sigma <- 10^df$`10LogSSDSlope(ug/L)-SigmaChronic NOEC`

# replace slope when indicated
#"Slope too low to be realistic, replace by 0.7"
replace_slope_idx <- grepl('^Slope', df$`Remark on chronic SSD slope`)
df$`10LogSSDSlope(ug/L)-SigmaChronic NOEC`[replace_slope_idx] <- 0.7 

df$ssd_5th = df$`10LogSSDMedianConcentration(ug/L)-MuChronic NOEC` - 1.645*df$`10LogSSDSlope(ug/L)-SigmaChronic NOEC`
df$ssd_1st = df$`10LogSSDMedianConcentration(ug/L)-MuChronic NOEC` - 2.33*df$`10LogSSDSlope(ug/L)-SigmaChronic NOEC`
#df$color = ifelse(log10(df$km_50th)>df$ssd_5th, "red", ifelse(log10(df$km_99th)>df$ssd_1st,"yellow","green"))
df$category = ifelse(log10(df$km_50th)>df$ssd_5th, "a",
                  ifelse(log10(df$km_50th)>df$ssd_1st,"b",
                         ifelse(log10(df$km_95th)>df$ssd_5th,"c",
                                ifelse(log10(df$km_95th)>df$ssd_1st,"d",
                                       ifelse(log10(df$km_99th)>df$ssd_5th,"e",
                                              ifelse(log10(df$km_99th)>df$ssd_1st,"f","g"
                                              )
                                       )
                                )
                         )
            )
)


df <- df[order(df$pct_cen),] 
df$idx <- seq.int(nrow(df))


final_chem_ids <- unique(df$chem_id)



# old figure
ggplot(data = df,
       aes(x=ssd_5th,
           y=log10(km_95th))) + 
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  geom_text(data = filter(df, log10(km_95th)>ssd_5th),aes(label=PREFERRED_NAME), vjust = -0.5, hjust = 0) +
  #scale_x_log10() + 
  #scale_y_log10() + 
  theme_bw()

df$ranking <- df$ssd_1st - log10(df$km_99th)
df <- df[order(df$category,df$ranking),] 
#df$idx <- seq.int(nrow(df))

#https://github.com/tidyverse/forcats/issues/16
fct_reordern <- function(.f, ..., .desc=FALSE, ordered=NA) {
  stopifnot(length(.desc) %in% c(1, ...length()))
  stopifnot(all(.desc %in% c(TRUE, FALSE)))
  f <- forcats:::check_factor(.f)
  # The radix method is used to support a vector of .desc (other methods only
  # support scalar values for .desc).
  new_order <- base::order(..., method="radix", decreasing=.desc)
  f_sorted <- f[new_order]
  fct_inorder(f=f_sorted, ordered=ordered)
}

# figure 6
ggplot(data=df) + 
  geom_segment(aes(x = fct_reordern(DTXSID, category, ranking),
                   xend = fct_reordern(DTXSID, category, ranking),
                   y = log10(km_50th),
                   yend = log10(km_99th),color = category),
               size=1, lty="11", alpha = 0.9
  ) +
  geom_segment(aes(x = fct_reordern(DTXSID, category, ranking),
                   xend = fct_reordern(DTXSID, category, ranking),
                   y = `10LogSSDMedianConcentration(ug/L)-MuChronic NOEC`,
                   yend = ssd_1st,color = category),
               position = position_nudge(x = 0.4),size=1, lty="31", alpha = 0.9) +  #, 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust=1, size=5),
        legend.title = element_blank(),legend.position="none")+
  labs(x = "Chemicals from highest to lowest bioactivity exposure overlap",
       y = "Concentrations (log \u03bcg/L)") + 
  scale_color_viridis(option='plasma',discrete=TRUE)#+
#scale_color_manual(values = c("#009E73",'#e31a1d','#ff7f00'))

ggplot(data=df[which(df$color %in% c('red','yellow')),]) + 
  geom_segment(aes(x = fct_reorder(INPUT, ranking), xend = fct_reorder(INPUT, ranking), y = log10(km_50th), yend = log10(km_99th),color = color),
               size=1, lty="11", alpha = 0.8
  ) +
  geom_segment(aes(x = fct_reorder(INPUT, ranking), xend = fct_reorder(INPUT, ranking), y = `10LogSSDMedianConcentration(ug/L)-MuChronic NOEC`, yend = ssd_1st,color = color),
               position = position_nudge(x = 0.4), size=1, lty="31", alpha = 0.4) +  #,\
  geom_point(aes(x = fct_reorder(INPUT, ranking), y=log10(km_99th)), shape=2) +
  geom_point(aes(x = fct_reorder(INPUT, ranking), y=ssd_1st)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust=1, size=5),
        legend.title = element_blank(),legend.position="none")+
  labs(x = "Chemicals from highest to lowest bioactivity exposure overlap",
       y = "Concentrations (log \u03bcg/L)") 
  
  
  +
  scale_color_manual(values = c("#009E73",'#e31a1d','#ff7f00'))


# Two chemicals are prioritized on the basis of this lowest risk condition: acetic acid and bifenthrin.
# Their highest observed environmental concentrations exceed estimated no-effect
# concentrations in the most sensitive species.

# G. Case 2: high estimated values compared with lowest part of the SSD

ggplot(data = df) + #,aes(x=ssd_5th, y=log10(mle_95th)
  geom_segment(aes(x = ssd_1st, xend = ssd_5th, y = log(km_95th), yend = log(km_99th), color=color),
  ) +
  #geom_point() +
  #geom_abline(slope=1, intercept=0) +
  #geom_text(data = filter(df, log10(mle_95th)>ssd_5th),aes(label=PREFERRED_NAME), vjust = -0.5, hjust = 0) +
  theme_bw() +
  scale_color_manual(values = c("green", "red", "yellow"))+ 
  labs(x = "Lowest predicted no-effect concentrations across freshwater species (log \u03bcg/L)",
       y = "Highest observed concentrations (log \u03bcg/L)")

# One additional chemical is prioritized on the basis of this risk condition: 
# chloroform. 

# old figure
ggplot(data = df,
       aes(x=ssd_5th,
           y=log10(mle_95th-mle_50th+mle_50th_ucl))) + 
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  geom_text(data = filter(df, log10(mle_95th-mle_50th+mle_50th_ucl)>ssd_5th),aes(label=PREFERRED_NAME), vjust = -0.5, hjust = 0) +
  theme_bw()

library("readr")
# hist of bloq
pct <- dslv_DF_smaller[which(dslv_DF_smaller$chem_id %in% final_chem_ids),]  %>%
  group_by(chem_id, bloq) %>%
  summarise(n = n())
data_wide <- spread(pct, bloq, n)
names(data_wide)[names(data_wide) == '0'] <- "result"
names(data_wide)[names(data_wide) == '1'] <- "bloq"
        
data_wide$pct_bloq <- data_wide$bloq/(data_wide$result+data_wide$bloq)
data_wide$log10_count <- log10(data_wide$result+data_wide$bloq)
write.csv(data_wide, "pct_12nov2021.csv")

bloq_avg <- DF_w_HUC[which(DF_w_HUC$loq_type %in% c(1,3) & DF_w_HUC$bloq == 1),] %>%
  group_by(chem_id, loq_type) %>%
  summarize(mean_conc = median(obs, na.rm = TRUE))

res_avg <- DF_w_HUC[which(DF_w_HUC$bloq == 0),] %>%
  group_by(chem_id) %>%
  summarize(mean_conc = median(obs, na.rm = TRUE))

wide = bloq_avg %>% spread(loq_type, mean_conc)

plot_limits <- merge(res_avg, wide, by="chem_id")
plot_limits$rpt_gt_tech <- ifelse(plot_limits$`3` > plot_limits$`1`,1,0)
plot_limits$tech_gt_mean <- ifelse(plot_limits$rpt_gt_tech == 0,
                                   ifelse(plot_limits$`1` > plot_limits$mean_conc, 1, 0),
                                   NA)

# figure 2
ggplot(data = plot_limits[which(plot_limits$chem_id %in% final_chem_ids),]) + 
  geom_point(aes(x=mean_conc,y=`1`), color="#ff8000", alpha = 0.9, shape = 17) + #technical
  geom_point(aes(x=mean_conc,y=`3`), color="#644cff", alpha = 0.9, shape = 15) + #reporting
  geom_segment(aes(x = mean_conc, xend = mean_conc, y = `1`, yend = `3`),
               color = 'azure4', alpha = 0.9) +
  geom_abline(slope=1, intercept = 0, linetype='dashed') + 
  scale_y_log10() + 
  scale_x_log10() + 
  theme_bw() + 
  labs(x = "Median measured value per chemical (log \u03bcg/L)",
       y= "Median left-censoring limit per chemical (log \u03bcg/L)")

sum(na.omit(plot_limits[which(plot_limits$chem_id %in% final_chem_ids),])$rpt_gt_tech)
length(na.omit(plot_limits[which(plot_limits$chem_id %in% final_chem_ids),])$rpt_gt_tech)

#! export for paper 2
first_res <- merge(x = wider_mle_df, y = wider_bulk_mle_df, by = "chem_id", all = TRUE)
second_res <- merge(x = first_res, y = wider_bloq_bulk_mle_df, by = "chem_id", all = TRUE)
second_res <- second_res[!is.na(second_res$chem_id),]
third_res <- merge(x = second_res, y = wider_bloq_mle_df, by = "chem_id", all = TRUE)
third_res <- third_res[!is.na(third_res$chem_id),]
write.csv(third_res,"all_chem_res_2021dec05.csv")