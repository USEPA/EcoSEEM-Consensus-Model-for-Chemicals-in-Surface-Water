#https://www4.stat.ncsu.edu/~reich/st590/code/

library('RMySQL')
library('dbConnect')
library('RISmed')
library('data.table')
library('dplyr')
library('rjags')
library('ggplot2')
library('readxl')
library('knitr')
library('mcmcplots')
library('purrr')

# Connect
mydb = dbConnect(MySQL(), user='rsayre01', password=password, dbname='sbox_rsayre01_ecoseem_v2', host='mysql-res1.epa.gov')

# Retrieve info
query = sprintf("SELECT DTXSID, ResultMeasureValue_ug_L, LOQ_ug_L FROM result WHERE DTXSID IS NOT NULL AND DTXSID != '-' AND id NOT IN (select id from result WHERE ActivityMediaSubdivisionName IN ('Interstitial','Elutriation','Finished Water','Industrial Waste','Wastewater Treatment Plant Effluent', 'Groundwater', 'Municipal Waste', 'Wet Fall Material')) AND id NOT IN (select id from result WHERE ActivityTypeCode IN ('Quality Control Sample-Field Spike', 'Quality Control Sample-Field Blank', 'Quality Control Sample-Lab Matrix Spike', 'Quality Control Sample-Trip Blank', 'Quality Control Sample-Spike Solution')) AND id NOT IN (select id from result WHERE ResultDetectionConditionText = 'Systematic Contamination')")
res_set = dbGetQuery(mydb, query)

# Convert zero, blank string to null
res_set$ResultMeasureValue_ug_L[res_set$ResultMeasureValue_ug_L == '0'] <- NA
res_set$LOQ_ug_L[res_set$LOQ_ug_L == '0'] <- NA


#log??
# Summarize
#summary_fxns = list(ln_geom_mean = function(x) mean(log(as.numeric(x), na.rm = TRUE)),
#                    ln_geom_std = function(x) sd(log(as.numeric(x), na.rm = TRUE))
#                    )
#sapply(summary_fxns, function(fn){res_set %>% summarise_all(fn)})
                    
loqs <- res_set %>%
   group_by(DTXSID) %>%
   summarize(ln_geom_mean_loq = mean(log(as.numeric(LOQ_ug_L)), na.rm = TRUE))

rmvs <- res_set %>%
   group_by(DTXSID) %>%
   summarize(ln_geom_mean_rmv = mean(log(as.numeric(ResultMeasureValue_ug_L)), na.rm = TRUE))

rmv_sds <- res_set %>%
   group_by(DTXSID) %>%
   summarize(ln_geom_sd_rmv = sd(log(as.numeric(ResultMeasureValue_ug_L)), na.rm = TRUE))

df <- merge(x=rmvs, y=loqs, by='DTXSID')

# log things

# Set ID as factor
df$DTXSID <- factor(df$DTXSID, levels = df$DTXSID[order(df$ln_geom_mean_rmv)])

# Plot RMV vs LOQ
ggplot(df, aes(x=DTXSID,y=ln_geom_mean_rmv)) +
   geom_point(color='#304389') +
   geom_segment(aes(xend=DTXSID, yend=ln_geom_mean_loq), color='azure3', alpha=0.5) +
   geom_point(aes(x=DTXSID,y=ln_geom_mean_loq), color='#5b87d9') +
   xlab('DTXSID') +
   ylab('Geometric mean of concentration (ln ug/L)') 

model_res = read_excel('L:\\Lab\\NCCT_ExpoCast\\ExpoCast2019\\ECOSEEM1\\ecoseem_13august2019.xlsx')
model_rmv_df <- merge(x=df, y=model_res, by='DTXSID')

# Set null NPV to min "max_millbs" value
model_rmv_df$NPV_max_millbs[is.na(model_rmv_df$NPV_max_millbs)] <- '0.025' 
model_rmv_df$ln_NPV = log(as.numeric(model_rmv_df$NPV_max_millbs))

# Scale to zero mean, unit variance
s_rmv <- scale(model_rmv_df$ln_geom_mean_rmv)
rmv_mean <- attr(s_rmv, "scaled:center")
rmv_sd <- attr(s_rmv, "scaled:scale")
model_rmv_df$ln_geom_mean_rmv_0m <- s_rmv
model_rmv_df$ln_geom_mean_loq_0m <- (model_rmv_df$ln_geom_mean_loq - rmv_mean)/rmv_sd
model_rmv_df$ln_USEtox_0m <- scale(model_rmv_df$ln_USEtox)
model_rmv_df$ln_RAIDAR_0m <- scale(model_rmv_df$ln_RAIDAR)
model_rmv_df$ln_EXAMS_0m <- scale(model_rmv_df$ln_EXAMS_mg_l_kg_hr)
model_rmv_df$ln_NPV_0m <- scale(model_rmv_df$ln_NPV)

# Remove rows with no measured and no LOQ
model_rmv_df <- model_rmv_df[!(is.nan(model_rmv_df$ln_geom_mean_rmv_0m) & is.nan(model_rmv_df$ln_geom_mean_loq_0m)), ]

model_rmv_df$ln_geom_mean_loq_0m[is.nan(model_rmv_df$ln_geom_mean_loq_0m)] <- NA
model_rmv_df$ln_geom_mean_rmv_0m[is.nan(model_rmv_df$ln_geom_mean_rmv_0m)] <- NA
model_rmv_df$isAboveLOQ = as.numeric( model_rmv_df$ln_geom_mean_rmv_0m > model_rmv_df$ln_geom_mean_loq_0m )
#model_rmv_df$ln_geom_mean_loq_0m[is.na(model_rmv_df$ln_geom_mean_loq_0m)] <- 0.9*model_rmv_df$ln_geom_mean_rmv_0m[is.na(model_rmv_df$ln_geom_mean_loq_0m)]
model_rmv_df$isAboveLOQ[is.na(model_rmv_df$ln_geom_mean_loq_0m)] <- 1
model_rmv_df$isAboveLOQ[is.na(model_rmv_df$isAboveLOQ)] <- 0

# Fewer columns
complete <- model_rmv_df[, c(1,49,50,55,51,52,53,54)]

# Train test split
nchem <- length(complete$DTXSID)
set.seed(42)
test <- sample(1:nchem, nchem * 0.25, replace=T)

Y     <- as.data.frame(complete[, 2])
X     <- complete[, 5:8]
X_U   <- as.data.frame(complete[, 5])
X_R   <- as.data.frame(complete[, 6])
X_E   <- as.data.frame(complete[, 7])
X_N   <- as.data.frame(complete[, 8])

Yo    <- Y[-test,]    # Observed data
LOQo  <- as.data.frame(complete[, 3])[-test,]
LOQo[is.na(LOQo )] <- Yo[is.na(LOQo )]*0.9
isAboveLOQ  <- as.data.frame(complete[, 4])[-test,]
Xo    <- X[-test,]

Yp    <- Y[test,]     # Set aside for prediction
Xp    <- X[test,]

no    <- length(Yo)
np    <- length(Yp)
p     <- ncol(Xo)

# JAGS Gaussian model
model_string <- "model{

# Likelihood
for(i in 1:no){
  isAboveLOQ[i] ~ dinterval(Yo[i], LOQo[i])
  Yo[i]   ~ dnorm(muo[i],inv.var)
  muo[i] <- alpha + inprod(Xo[i,],beta[])
}

# Prediction
for(i in 1:np){
# We don't need to add a term for normally distributed error because we are not
# calculating a likelihood for these data.
  Yp[i] <- alpha + inprod(Xp[i,],beta[])
}

# Priors
for(j in 1:p){
  beta[j] ~ dnorm(0,0.0001)
}
alpha     ~ dnorm(0, 0.01)
inv.var   ~ dgamma(0.001, 0.001)
sigma     <- 1/sqrt(inv.var)
}"

model <- jags.model(textConnection(model_string), 
                    n.chains=4,
                    data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,np=np,p=p,Xo=Xo,Xp=Xp))
# burn in
update(model, 10000, progress.bar="none")
samp <- coda.samples(model, 
                     variable.names=c("beta","sigma","Yp","alpha","muo"), 
                     n.iter=20000, progress.bar="none")

Ypred_g <- samp[[1]][,1:np]
coeff_g <- samp[[1]][,(np+1):(np+1+length(X))]
muo_g <- samp[[1]][,(np+2+length(X)):(ncol(samp[[1]])-1)]

rsq_g <- vector()
fstat_g <- vector()
for(j in 1:20000){
   this_vec <- muo_g[j,]
   this_rsq <- summary(lm(Yo ~ this_vec))$r.squared
   this_fstat <- summary(lm(Yo ~ this_vec))$fstatistic[[1]] 
   rsq_g <- c(rsq_g, this_rsq)
   fstat_g <- c(fstat_g, this_fstat)
}

# JAGS BLASSO model
model_string_b <- "model{

# Likelihood
for(i in 1:no){
  isAboveLOQ[i] ~ dinterval(Yo[i], LOQo[i])
  Yo[i]   ~ dnorm(muo[i],inv.var)
  muo[i] <- alpha + inprod(Xo[i,],beta[])
}

# Prediction
for(i in 1:np){
# We don't need to add a term for normally distributed error because we are not
# calculating a likelihood for these data.
  Yp[i] <- alpha + inprod(Xp[i,],beta[])
}

# Priors
for(j in 1:p){
  beta[j] ~ ddexp(0,inv.var.b)
}
alpha     ~ dnorm(0, 0.01)
inv.var   ~ dgamma(0.01, 0.01)
inv.var.b ~ dgamma(0.01, 0.01)
}"

model_b <- jags.model(textConnection(model_string_b),
                      n.chains=4,
                      data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,np=np,p=p,Xo=Xo,Xp=Xp))

update(model_b, 10000, progress.bar="none")
samp_b <- coda.samples(model_b, 
                     variable.names=c("beta","Yp","alpha","muo"), 
                     n.iter=20000, progress.bar="none")
#dic_b  <- dic.samples(model_b, 
#                    variable.names=c("beta", "alpha"), 
#                    n.iter=20000, progress.bar="none")


Ypred_b <- samp_b[[1]][,1:np]
coeff_b <- samp_b[[1]][,(np+1):(np+1+length(X))]
muo_b <- samp_b[[1]][,(np+2+length(X)):(ncol(samp_b[[1]]))]

rsq_b <- vector()
fstat_b <- vector()
for(j in 1:20000){
   this_vec <- muo_b[j,]
   this_rsq <- summary(lm(Yo ~ this_vec))$r.squared
   this_fstat <- summary(lm(Yo ~ this_vec))$fstatistic[[1]] 
   rsq_b <- c(rsq_b, this_rsq)
   fstat_b <- c(fstat_b, this_fstat)
}

# Gaussian priors with a prior on the variance of beta, ??j???Normal(0,??2b) where ??2b???InvGamma(0.01,0.01).
model_string_g2 <- "model{

# Likelihood
for(i in 1:no){
isAboveLOQ[i] ~ dinterval(Yo[i], LOQo[i])
Yo[i]   ~ dnorm(muo[i],inv.var)
muo[i] <- alpha + inprod(Xo[i,],beta[])
}

# Prediction
for(i in 1:np){
# We don't need to add a term for normally distributed error because we are not
# calculating a likelihood for these data.
Yp[i] <- alpha + inprod(Xp[i,],beta[])
}

# Prior for beta
for(j in 1:p){
beta[j] ~ dnorm(0,inv.var.b)
}

# Prior for the inverse variance
inv.var   ~ dgamma(0.01, 0.01)
inv.var.b ~ dgamma(0.01, 0.01)
alpha     ~ dnorm(0, 0.01)

}"

model_g2 <- jags.model(textConnection(model_string_g2), 
                       n.chains=4,
                       data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,np=np,p=p,Xo=Xo,Xp=Xp))

update(model_g2, 10000, progress.bar="none")
samp_g2 <- coda.samples(model_g2, 
                        variable.names=c("beta","Yp","alpha","muo"), 
                        n.iter=20000, progress.bar="none")
#dic_g2  <- dic.samples(model_g2, 
#                    variable.names=c("beta", "alpha"), 
#                    n.iter=20000, progress.bar="none")


Ypred_g2 <- samp_g2[[1]][,1:np]
coeff_g2 <- samp_g2[[1]][,(np+1):(np+1+length(X))]
muo_g2 <- samp_g2[[1]][,(np+2+length(X)):(ncol(samp_g2[[1]]))]

rsq_g2 <- vector()
fstat_g2 <- vector()
for(j in 1:20000){
   this_vec <- muo_g2[j,]
   this_rsq <- summary(lm(Yo ~ this_vec), na.rm=TRUE)$r.squared
   this_fstat <- summary(lm(Yo ~ this_vec))$fstatistic[[1]] 
   rsq_g2 <- c(rsq_g2, this_rsq)
   fstat_g2 <- c(fstat_g2, this_fstat)
}

for(index in 1:p){    
   d1 <- density(coeff_g[,index+1])
   d2 <- density(coeff_g2[,index+1])
   d4 <- density(coeff_b[,index+1])
   mx <- max(d1$y,d2$y,d3$y,d4$y)
   
   plot(d1,ylim=c(0,mx),xlab="Beta",ylab="Posterior density",main=colnames(X)[index])
   lines(d2,col=2)
   lines(d4,col=3)
   legend("topright",c("Gaussian 1", "Gaussian 2", "LASSO"),lty=1,col=1:3,inset=0.05)
}

post_mn1   <- apply(Ypred_g,2,mean)
post_sd1   <- apply(Ypred_g,2,sd)
post_low1  <- apply(Ypred_g,2,quantile,0.05)
post_high1 <- apply(Ypred_g,2,quantile,0.95)

MSE1   <- mean((post_mn1-Yp)^2, na.rm=TRUE)
BIAS1  <- mean(post_mn1-Yp, na.rm=TRUE)
AVESD1 <- mean(post_sd1, na.rm=TRUE)
COV1   <- mean(Yp>post_low1 & Yp<post_high1, na.rm=TRUE)

post_mn2   <- apply(Ypred_g2 ,2,mean)
post_sd2   <- apply(Ypred_g2 ,2,sd)
post_low2  <- apply(Ypred_g2 ,2,quantile,0.05)
post_high2 <- apply(Ypred_g2 ,2,quantile,0.95)

MSE2   <- mean((post_mn2-Yp)^2, na.rm=TRUE)
BIAS2  <- mean(post_mn2-Yp, na.rm=TRUE)
AVESD2 <- mean(post_sd2, na.rm=TRUE)
COV2   <- mean(Yp>post_low2 & Yp<post_high2, na.rm=TRUE)

post_mn4   <- apply(Ypred_b ,2,mean)
post_sd4   <- apply(Ypred_b ,2,sd)
post_low4  <- apply(Ypred_b ,2,quantile,0.05)
post_high4 <- apply(Ypred_b ,2,quantile,0.95)

MSE4   <- mean((post_mn4-Yp)^2, na.rm=TRUE)
BIAS4  <- mean(post_mn4-Yp, na.rm=TRUE)
AVESD4 <- mean(post_sd4, na.rm=TRUE)
COV4   <- mean(Yp>post_low4 & Yp<post_high4, na.rm=TRUE)

Ypred_0 <- rep(mean(Yp, na.rm=TRUE), times=4640000, na.rm=TRUE)

post_mn0   <- apply(Ypred_0 ,2,mean)
post_sd0   <- apply(Ypred_0 ,2,sd)
post_low0  <- apply(Ypred_0 ,2,quantile,0.05)
post_high0 <- apply(Ypred_0 ,2,quantile,0.95)

MSE0   <- mean((mean(Ypred_0)-Yp)^2, na.rm=TRUE)
BIAS0  <- mean(mean(Ypred_0)-Yp, na.rm=TRUE)
AVESD0 <- mean(sd(Ypred_0), na.rm=TRUE)
COV0   <- mean(Yp>post_low0 & Yp<post_high0, na.rm=TRUE)

MSE    <- c(MSE1,MSE2,MSE4)
BIAS   <- c(BIAS1,BIAS2,BIAS4)
AVESD  <- c(AVESD1,AVESD2,AVESD4)
COV    <- c(COV1,COV2,COV4)


x = list(Ypred_g, Ypred_g2, Ypred_b)

OUTPUT <- cbind(MSE,BIAS,AVESD,COV)
kable(OUTPUT,digits=3)

mm = read_excel('L:\\Lab\\NCCT_ExpoCast\\ExpoCast2019\\ECOSEEM1\\molar_mass_for_ecoseem.xls')

Xo=X_U[-test,];Xp=X_U[test,]

model_U <- jags.model(textConnection(model_string_U),
                      n.chains=4,
                      data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,Xo=Xo))
update(model_U, 10000, progress.bar="none")

samp_U <- coda.samples(model_U, 
                       variable.names=c("beta","alpha","muo"), 
                       n.iter=20000, progress.bar="none")

coeff_U <- samp_U[[1]][,1:2]
muo_U <- samp_U[[1]][,3:ncol(samp_U[[1]])]

rsq_U  <- summary(lm(Yo ~ muo_U[j,]))$r.squared
fstat_U <- summary(lm(Yo ~ muo_U[j,]))$fstatistic[[1]] 

# JAGS BLASSO model
model_string_R <- "model{

# Likelihood
for(i in 1:no){
isAboveLOQ[i] ~ dinterval(Yo[i], LOQo[i])
Yo[i]   ~ dnorm(muo[i],inv.var)
muo[i] <- alpha + beta*Xo[i]
}

# Prior
beta ~ ddexp(0,inv.var.b)

alpha     ~ dnorm(0, 0.01)
inv.var   ~ dgamma(0.01, 0.01)
inv.var.b ~ dgamma(0.01, 0.01)

}"

Xo=X_R[-test,];Xp=X_R[test,]

model_R <- jags.model(textConnection(model_string_R),
                      n.chains=4,
                      data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,Xo=Xo))
update(model_R, 10000, progress.bar="none")

samp_R <- coda.samples(model_R, 
                       variable.names=c("beta","alpha","muo"), 
                       n.iter=20000, progress.bar="none")

coeff_R <- samp_R[[1]][,1:2]
muo_R <- samp_R[[1]][,3:ncol(samp_R[[1]])]

rsq_R  <- summary(lm(Yo ~ muo_R[j,]))$r.squared
fstat_R <- summary(lm(Yo ~ muo_R[j,]))$fstatistic[[1]] 

# JAGS BLASSO model
model_string_E <- "model{

# Likelihood
for(i in 1:no){
isAboveLOQ[i] ~ dinterval(Yo[i], LOQo[i])
Yo[i]   ~ dnorm(muo[i],inv.var)
muo[i] <- alpha + beta*Xo[i]
}

# Prior
beta ~ ddexp(0,inv.var.b)

alpha     ~ dnorm(0, 0.01)
inv.var   ~ dgamma(0.01, 0.01)
inv.var.b ~ dgamma(0.01, 0.01)

}"

Xo=X_E[-test,];Xp=X_E[test,]

model_E <- jags.model(textConnection(model_string_E),
                      n.chains=4,
                      data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,Xo=Xo))
update(model_E, 10000, progress.bar="none")

samp_E <- coda.samples(model_E, 
                       variable.names=c("beta","alpha","muo"), 
                       n.iter=20000, progress.bar="none")

coeff_E <- samp_E[[1]][,1:2]
muo_E <- samp_E[[1]][,3:ncol(samp_E[[1]])]

rsq_E  <- summary(lm(Yo ~ muo_E[j,]))$r.squared
fstat_E <- summary(lm(Yo ~ muo_E[j,]))$fstatistic[[1]] 

# JAGS BLASSO model
model_string_N <- "model{

# Likelihood
for(i in 1:no){
isAboveLOQ[i] ~ dinterval(Yo[i], LOQo[i])
Yo[i]   ~ dnorm(muo[i],inv.var)
muo[i] <- alpha + beta*Xo[i]
}

# Prior
beta ~ ddexp(0,inv.var.b)

alpha     ~ dnorm(0, 0.01)
inv.var   ~ dgamma(0.01, 0.01)
inv.var.b ~ dgamma(0.01, 0.01)

}"

Xo=X_N[-test,];Xp=X_N[test,]

model_N <- jags.model(textConnection(model_string_N),
                      n.chains=4,
                      data = list(Yo=Yo,LOQo=LOQo,isAboveLOQ=isAboveLOQ,no=no,Xo=Xo))
update(model_N, 10000, progress.bar="none")

samp_N <- coda.samples(model_N, 
                       variable.names=c("beta","alpha","muo"), 
                       n.iter=20000, progress.bar="none")

coeff_N <- samp_N[[1]][,1:2]
muo_N <- samp_N[[1]][,3:ncol(samp_N[[1]])]

rsq_N  <- summary(lm(Yo ~ muo_N[j,]))$r.squared
fstat_N <- summary(lm(Yo ~ muo_N[j,]))$fstatistic[[1]] 

# Add SD back in
df <- merge(x=as.data.frame(complete[, 1:3])[~test,], y=rmv_sds, by='DTXSID')

df[, c(4)] <- (df[, c(4)] - rmv_mean)/rmv_sd #! this should be the same scale as RMV - fix this
df = merge(x=df, y=model_res, on='DTXSID')

plots <- NULL
for (i in 1:length(Yo)){
   plt <- ggplot(NULL, aes(x=Yo[i],y=mean(muo_b[,i])))+
      geom_segment(aes(x=Yo[i],
                       xend=Yo[i],
                       y=quantile(muo_b[,i],0.05),
                       yend=quantile(muo_b[,i],0.95)))
   plots <- rbind(plots, plt)
}
                                                                                          yend=quantile(muo_b[,i],0.95)))}
post_mn   <- apply(muo_g,2,mean)
MSE   <- mean((post_mn-Yo)^2, na.rm=TRUE)
# Plot predictions on validation set vs withheld values
ggplot(NULL, aes(x=preds$statistics.Mean[1:191],y=df$ln_geom_mean_rmv_0m)) +
   geom_point(aes(x=Ytest_B[1:191],y=df$ln_geom_mean_rmv_0m), shape=3) +
   geom_point() +
   geom_abline(intercept=0, slope=1, linetype=2) +
   geom_segment(aes(x=preds$statistics.Mean[1:191]-preds$statistics.SD[1:191],
                    xend=preds$statistics.Mean[1:191]+preds$statistics.SD[1:191],
                    yend=df$ln_geom_mean_rmv_0m)) +
   geom_segment(aes(xend=preds$statistics.Mean[1:191],
                    y=df$ln_geom_mean_rmv_0m-df$ln_geom_sd_rmv,
                    yend=df$ln_geom_mean_rmv_0m+df$ln_geom_sd_rmv)) +
   title('Comparison of ECOSEEM results to validation set') +
   xlab('ECOSEEM regression result (ln ug/L, scaled)') +
   ylab('Geometric mean of measured concentration (ln ug/L, scaled)')



badly_underpred <- NULL
for (i in 1:length(Yo)) {
   lt_sd <- Yo[i] - sd_df[i,9]
   mu <- mean(muo_b[,i])
   if(! is.na(lt_sd) && lt_sd > mu) {
      badly_underpred <- c(badly_underpred, i)
   } else {
      next}
}

#Association between physchemprop and predictivity
model_rmv_df[,17] <- as.numeric(model_rmv_df[,17])
model_rmv_df[,18] <- as.numeric(model_rmv_df[,18])
model_rmv_df[,19] <- as.numeric(model_rmv_df[,19])
model_rmv_df[,21] <- as.numeric(model_rmv_df[,21])
model_rmv_df[,25] <- as.numeric(model_rmv_df[,25])
model_rmv_df[,28] <- as.numeric(model_rmv_df[,28])
model_rmv_df[,31] <- as.numeric(model_rmv_df[,31])
model_rmv_df[,32] <- as.numeric(model_rmv_df[,32])
model_rmv_df[,33] <- as.numeric(model_rmv_df[,33])
model_rmv_df[,34] <- as.numeric(model_rmv_df[,34])
model_rmv_df[,35] <- as.numeric(model_rmv_df[,35])
model_rmv_df[,36] <- as.numeric(model_rmv_df[,36])
model_rmv_df[,37] <- as.numeric(model_rmv_df[,37])
model_rmv_df[,39] <- as.numeric(model_rmv_df[,39])
model_rmv_df[,40] <- as.numeric(model_rmv_df[,40])
model_rmv_df[,41] <- as.numeric(model_rmv_df[,41])
model_rmv_df[,42] <- as.numeric(model_rmv_df[,42])
model_rmv_df[,43] <- as.numeric(model_rmv_df[,43])
model_rmv_df[,44] <- as.numeric(model_rmv_df[,44])
model_rmv_df[,45] <- as.numeric(model_rmv_df[,45])

plot_physchemprop <- model_rmv_df[-test,c(1,2,17,18,19,21,25,28,31,32,33,34,35,36,37,39,40,41,42,43,44,45)]
plot_physchemprop$lower_ci <- lower_ci
plot_physchemprop$upper_ci <- upper_ci
