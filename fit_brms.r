library(brms)	#package to fit hierarchical models in Stan with simplefied syntax; will also load rstan								
options(mc.cores = parallel::detectCores())     

dat <- read.csv('data.csv')                 
cat <- read.csv('categories.csv')            

studies <- unique(dat$study.sub)             
 
dat[,4:9]     <- NA                              
colnames(dat) <- c('ESD','SV','study','coastal_open_lab','situ_not_flow_indirect','aggregate_fecal','diatom','ballast','warm_cold')

for(i in 1:length(studies)){                                              
  dat[dat$study==studies[i],4:9]<- cat[cat$study.sub==studies[i],2:7]   
} 
dat$y <- log(dat$SV)   
dat$x <- log(dat$ESD)

dat_no_fecal <- dat[dat$aggregate_fecal==1,]
#############################################################
## FIT MODELS ###############################################
#############################################################
niter   <- 5000
nchains <- 4

##--GLOBAL--#############
prior_global          <- get_prior(y ~ x + (x|study),data=dat)
prior_global$prior[2] <- 'uniform(1E-10,2)'
prior_global$prior[9] <- 'uniform(1E-10,2)'
fit_global            <- brm(y ~ x + (x|study),open_progress=TRUE,
						                       data=dat,
									           iter=niter,
											   chains=nchains, 
											   prior=prior_global)

##--COASTAL_OPEN_LAB--#############
prior_coastal_open_lab           <- get_prior(y ~ x + (x|coastal_open_lab/study), data=dat)
prior_coastal_open_lab$prior[2]  <- 'uniform(1E-10,2)'
prior_coastal_open_lab$prior[10] <- 'uniform(1E-10,2)'
prior_coastal_open_lab$prior[13] <- 'uniform(1E-10,1)'
fit_coastal_open_lab             <- brm(y ~ x + (x|coastal_open_lab/study), open_progress=TRUE,
																			data=dat,
																			iter=niter,
																			chains=nchains, 
																			prior=prior_coastal_open_lab)

##--IN-SITU--#############
prior_insitu_not_flow_in           <- get_prior(y ~ x + (x|situ_not_flow_indirect/study),data=dat)
prior_insitu_not_flow_in$prior[2]  <- 'uniform(1E-10,2)'
prior_insitu_not_flow_in$prior[10] <- 'uniform(1E-10,2)'
prior_insitu_not_flow_in$prior[13] <- 'uniform(1E-10,1)'
fit_insitu_not_flow_in             <- brm(y ~ x + (x|situ_not_flow_indirect/study), open_progress=TRUE,
																					data=dat,
																					iter=niter,
																				    chains=nchains,
																					prior=prior_insitu_not_flow_in)

##--AGGREGATES--#############
prior_agg_fecal           <- get_prior(y ~ x + (x|aggregate_fecal/study), data=dat)
prior_agg_fecal$prior[2]  <- 'uniform(1E-10,2)'
prior_agg_fecal$prior[10] <- 'uniform(1E-10,2)'
prior_agg_fecal$prior[13] <- 'uniform(1E-10,1)'
fit_agg_fecal             <- brm(y ~ x + (x|aggregate_fecal/study), open_progress=TRUE,
   	 															    data=dat,
																    iter=niter,
																	chains=nchains,
																    prior=prior_agg_fecal)

##--DIATOM--#############
prior_diatom_non           <- get_prior(y ~ x + (x|diatom/study),  data=dat_no_fecal)
prior_diatom_non$prior[2]  <- 'uniform(1E-10,2)'
prior_diatom_non$prior[10] <- 'uniform(1E-10,2)'
prior_diatom_non$prior[13] <- 'uniform(1E-10,1)'
fit_diatom_non             <- brm(y ~ x + (x|diatom/study), open_progress=TRUE, 
							         			            data=dat_no_fecal,
											                iter=niter,
												            chains=nchains,
					                                        prior=prior_diatom_non)

##--BALLAST--#############
prior_ballast_non           <- get_prior(y ~ x + (x|ballast/study),data=dat_no_fecal)
prior_ballast_non$prior[2]  <- 'uniform(1E-10,2)'
prior_ballast_non$prior[10] <- 'uniform(1E-10,2)'
prior_ballast_non$prior[13] <- 'uniform(1E-10,1)'
fit_ballast_non             <- brm(y ~ x + (x|ballast/study), open_progress=TRUE,
							         			              data=dat_no_fecal,
												              iter=niter,
												              chains=nchains,
												              prior=prior_ballast_non)

##--WARM/COLD--#############
prior_warm_cold           <- get_prior(y ~ x + (x|warm_cold/study),data=dat)
prior_warm_cold$prior[2]  <- 'uniform(1E-10,2)'
prior_warm_cold$prior[9]  <- 'uniform(1E-10,2)'
prior_warm_cold$prior[10] <- 'uniform(1E-10,2)'
prior_warm_cold$prior[12] <- 'uniform(1E-10,2)'
prior_warm_cold$prior[13] <- 'uniform(1E-10,1)'
fit_warm_cold             <- brm(y ~ x + (x|warm_cold/study), open_progress=TRUE, 
															  data=dat,
															  iter=niter,
															  chains=nchains, 
					                                          prior=prior_warm_cold)

