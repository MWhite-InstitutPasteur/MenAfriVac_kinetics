


#####################################
#####################################
##                                 ## 
##  ######  ####  #####     ####   ##
##    ##   ##  ## ##  ##   ##  ##  ##
##    ##   ###### #####       ##   ##
##    ##   ##  ## ##  ##   ##  ##  ##
##    ##   ##  ## #####     ####   ##
##                                 ##
#####################################
#####################################


TAB3 = matrix(NA, nrow=4, ncol=63)

rownames(TAB3) = c("PsA-TT-002_MenA", "PsA-TT-002_logELISA", "PsA-TT-003_MenA", "PsA-TT-003_logELISA")
colnames(TAB3) = c("beta_prim_PsA_med",     "beta_prim_PsA_low",     "beta_prim_PsA_high", 
                   "beta_boost_PsA_med",    "beta_boost_PsA_low",    "beta_boost_PsA_high", 
		       "beta_Hib_med",          "beta_Hib_low",          "beta_Hib_high", 
                   "d_cs_med",              "d_cs_low",              "d_cs_low_high", 
                   "d_cl_med",              "d_cl_low",              "d_cl_high", 
                   "d_a_med",               "d_a_low",               "d_a_high", 
		       "rho_prim_PsA_med",      "rho_prim_PsA_low",      "rho_prim_PsA_high", 
                   "rho_boost_PsA_med",     "rho_boost_PsA_low",     "rho_boost_PsA_high", 
		       "rho_Hib_med",           "rho_Hib_low",           "rho_Hib_high", 
                   "sd_beta_prim_PsA_med",  "sd_beta_prim_PsA_low",  "sd_beta_prim_PsA_high", 
                   "sd_beta_boost_PsA_med", "sd_beta_boost_PsA_low", "sd_beta_boost_PsA_high",  
		       "sd_beta_Hib_med",       "sd_beta_Hib_low",       "sd_beta_Hib_high", 
                   "sd_d_cs_med",           "sd_d_cs_low",           "sd_d_cs_high", 
                   "sd_d_cl_med",           "sd_d_cl_low",           "sd_d_cl_high", 
                   "sd_d_a_med",            "sd_d_a_low",            "sd_d_a_high", 
		       "sd_rho_prim_PsA_med",   "sd_rho_prim_PsA_low",   "sd_rho_prim_PsA_high", 
                   "sd_rho_boost_PsA_med",  "sd_rho_boost_PsA_low",  "sd_rho_boost_PsA_high", 
		       "sd_rho_Hib_med",        "sd_rho_Hib_low",        "sd_rho_Hib_high", 
                   "P_obs_med",             "P_obs_low",             "P_obs_high", 
                   "likelihood_med",        "likelihood_low",        "likelihood_high", 
                   "prior_med",             "prior_low",             "prior_high")

TAB3_median = matrix(NA, nrow=4, ncol=27)

rownames(TAB3_median) = c("PsA-TT-002_MenA", "PsA-TT-002_logELISA", "PsA-TT-003_MenA", "PsA-TT-003_logELISA")
colnames(TAB3_median) = c("beta_prim_PsA_med",     "beta_prim_PsA_low",     "beta_prim_PsA_high", 
                         "beta_boost_PsA_med",    "beta_boost_PsA_low",    "beta_boost_PsA_high", 
		             "beta_Hib_med",          "beta_Hib_low",          "beta_Hib_high", 
                         "d_cs_med",              "d_cs_low",              "d_cs_low_high", 
                         "d_cl_med",              "d_cl_low",              "d_cl_high", 
                         "d_a_med",               "d_a_low",               "d_a_high", 
		             "rho_prim_PsA_med",      "rho_prim_PsA_low",      "rho_prim_PsA_high", 
                         "rho_boost_PsA_med",     "rho_boost_PsA_low",     "rho_boost_PsA_high", 
		             "rho_Hib_med",           "rho_Hib_low",           "rho_Hib_high")



#######################
#######################
##                   ##
##  PsA-TT-002_MenA  ##
##                   ##
#######################
#######################

MCMC_pop_1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_002_mod2_2rho/Output/global_rep1.txt" )

MCMC_pop_2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_002_mod2_2rho/Output/global_rep2.txt" )


MCMC_pop = rbind( MCMC_pop_1[floor(0.5*nrow(MCMC_pop_1)):(nrow(MCMC_pop_1)-1),],
                  MCMC_pop_2[floor(0.5*nrow(MCMC_pop_2)):(nrow(MCMC_pop_2)-1),] )
colnames(MCMC_pop) = c("beta_prim_PsA", "beta_boost_PsA", "beta_Hib", 
		           "d_cs", "d_cl", "d_a", 
		           "rho_prim_PsA", "rho_boost_PsA", "rho_Hib", 
		           "sd_beta_prim_PsA", "sd_beta_boost_PsA", "sd_beta_Hib",  
		           "sd_d_cs", "sd_d_cl", "sd_d_a", 
		           "sd_rho_prim_PsA", "sd_rho_boost_PsA", "sd_rho_Hib", 
		           "P_obs", "likelihood", "prior")



N_glob_par = 9

TAB3_median[1,1:3]   = quantile( exp(MCMC_pop[,1]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,4:6]   = quantile( exp(MCMC_pop[,2]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,7:9]   = quantile( exp(MCMC_pop[,3]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,10:12] = quantile( exp(MCMC_pop[,4]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,13:15] = quantile( exp(MCMC_pop[,5]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,16:18] = quantile( exp(MCMC_pop[,6]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,19:21] = quantile( exp(MCMC_pop[,7])/(1+exp(MCMC_pop[,7])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,22:24] = quantile( exp(MCMC_pop[,8])/(1+exp(MCMC_pop[,8])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[1,25:27] = quantile( exp(MCMC_pop[,9])/(1+exp(MCMC_pop[,9])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )




for(p in 1:N_glob_par)
{
	MCMC_pop[,N_glob_par+p] = 1/sqrt(MCMC_pop[,N_glob_par+p])
}

for(p in 1:6)
{
	MCMC_pop[,p]            = exp( MCMC_pop[,p] + 0.5*MCMC_pop[,N_glob_par+p]^2 )
	MCMC_pop[,N_glob_par+p] = sqrt( exp(MCMC_pop[,N_glob_par+p]^2) -1 )*MCMC_pop[,p]
}


for(p in 7:9)
{
	for(k in 1:nrow(MCMC_pop))
	{
		f_m <- function(x){ 
			0.3989423*exp( -0.5*( ( log(x/(1-x))-MCMC_pop[k,p] )/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}

		f_m2 <- function(x){ 
			x*0.3989423*exp( -0.5*( (log(x/(1-x))-MCMC_pop[k,p])/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}
	
		tryCatch(
		{
			moment_1 <- integrate(f_m, lower=0, upper=1)$value
			moment_2 <- integrate(f_m2, lower=0, upper=1)$value
		}, error=function(e){ NULL }
		)

		MCMC_pop[k,p] <- moment_1
		MCMC_pop[k,N_glob_par+p] <- sqrt( moment_2 - moment_1^2 )
	}
}




for(j in 1:21)
{
	TAB3[1,(3*(j-1)+1):(3*j)] = quantile( MCMC_pop[,j], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
}





###########################
###########################
##                       ##
##  PsA-TT-002_logELISA  ##
##                       ##
###########################
###########################


MCMC_pop_1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_002_mod2_rawELISA/Output/global_rep1.txt" )

MCMC_pop_2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_002_mod2_rawELISA/Output/global_rep2.txt" )


MCMC_pop = rbind( MCMC_pop_1[floor(0.5*nrow(MCMC_pop_1)):(nrow(MCMC_pop_1)-1),],
                  MCMC_pop_2[floor(0.5*nrow(MCMC_pop_2)):(nrow(MCMC_pop_2)-1),] )

colnames(MCMC_pop) =  c("beta_prim_PsA", "beta_boost_PsA",  
		           "d_cs", "d_cl", "d_a", 
		           "rho_prim_PsA", "rho_boost_PsA",  
		           "sd_beta_prim_PsA", "sd_beta_boost_PsA",  
		           "sd_d_cs", "sd_d_cl", "sd_d_a", 
		           "sd_rho_prim_PsA", "sd_rho_boost_PsA",  
		           "P_obs", "likelihood", "prior")


N_glob_par = 7


TAB3_median[2,1:3]   = quantile( exp(MCMC_pop[,1]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[2,4:6]   = quantile( exp(MCMC_pop[,2]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[2,10:12] = quantile( exp(MCMC_pop[,3]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[2,13:15] = quantile( exp(MCMC_pop[,4]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[2,16:18] = quantile( exp(MCMC_pop[,5]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[2,19:21] = quantile( exp(MCMC_pop[,6])/(1+exp(MCMC_pop[,6])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[2,22:24] = quantile( exp(MCMC_pop[,7])/(1+exp(MCMC_pop[,7])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )





for(p in 1:N_glob_par)
{
	MCMC_pop[,N_glob_par+p] = 1/sqrt(MCMC_pop[,N_glob_par+p])
}

for(p in 1:5)
{
	MCMC_pop[,p]            = exp( MCMC_pop[,p] + 0.5*MCMC_pop[,N_glob_par+p]^2 )
	MCMC_pop[,N_glob_par+p] = sqrt( exp(MCMC_pop[,N_glob_par+p]^2) -1 )*MCMC_pop[,p]
}


for(p in 6:7)
{
	for(k in 1:nrow(MCMC_pop))
	{
		f_m <- function(x){ 
			0.3989423*exp( -0.5*( ( log(x/(1-x))-MCMC_pop[k,p] )/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}

		f_m2 <- function(x){ 
			x*0.3989423*exp( -0.5*( (log(x/(1-x))-MCMC_pop[k,p])/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}
	
		tryCatch(
		{
			moment_1 <- integrate(f_m, lower=0, upper=1)$value
			moment_2 <- integrate(f_m2, lower=0, upper=1)$value
		}, error=function(e){ NULL }
		)

		MCMC_pop[k,p] <- moment_1
		MCMC_pop[k,N_glob_par+p] <- sqrt( moment_2 - moment_1^2 )
	}
}



TAB3[2,1:3]   = quantile( MCMC_pop[,1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,4:6]   = quantile( MCMC_pop[,2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[2,10:12] = quantile( MCMC_pop[,3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,13:15] = quantile( MCMC_pop[,4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,16:18] = quantile( MCMC_pop[,5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,19:21] = quantile( MCMC_pop[,6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,22:24] = quantile( MCMC_pop[,7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[2,28:30] = quantile( MCMC_pop[,8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,31:33] = quantile( MCMC_pop[,9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[2,37:39] = quantile( MCMC_pop[,10], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,40:42] = quantile( MCMC_pop[,11], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,43:45] = quantile( MCMC_pop[,12], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,46:48] = quantile( MCMC_pop[,13], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,49:51] = quantile( MCMC_pop[,14], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[2,55:57] = quantile( MCMC_pop[,15], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,58:60] = quantile( MCMC_pop[,16], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[2,61:63] = quantile( MCMC_pop[,17], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )





#######################
#######################
##                   ##
##  PsA-TT-003_MenA  ##
##                   ##
#######################
#######################

MCMC_pop_1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_003_mod2/Output/global_rep1.txt" )

MCMC_pop_2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_003_mod2/Output/global_rep2.txt" )


MCMC_pop = rbind( MCMC_pop_1[floor(0.5*nrow(MCMC_pop_1)):(nrow(MCMC_pop_1)-1),],
                  MCMC_pop_2[floor(0.5*nrow(MCMC_pop_2)):(nrow(MCMC_pop_2)-1),] )
colnames(MCMC_pop) = c("beta_prim_PsA",
		           "d_cs", "d_cl", "d_a", 
		           "rho_prim_PsA",  
		           "sd_beta_prim_PsA",   
		           "sd_d_cs", "sd_d_cl", "sd_d_a", 
		           "sd_rho_prim_PsA",  
		           "P_obs", "likelihood", "prior")


N_glob_par = 5



TAB3_median[3,1:3]   = quantile( exp(MCMC_pop[,1]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[3,10:12] = quantile( exp(MCMC_pop[,2]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[3,13:15] = quantile( exp(MCMC_pop[,3]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[3,16:18] = quantile( exp(MCMC_pop[,4]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[3,19:21] = quantile( exp(MCMC_pop[,5])/(1+exp(MCMC_pop[,5])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )





for(p in 1:N_glob_par)
{
	MCMC_pop[,N_glob_par+p] = 1/sqrt(MCMC_pop[,N_glob_par+p])
}

for(p in 1:4)
{
	MCMC_pop[,p]            = exp( MCMC_pop[,p] + 0.5*MCMC_pop[,N_glob_par+p]^2 )
	MCMC_pop[,N_glob_par+p] = sqrt( exp(MCMC_pop[,N_glob_par+p]^2) -1 )*MCMC_pop[,p]
}


for(p in 5)
{
	for(k in 1:nrow(MCMC_pop))
	{
		f_m <- function(x){ 
			0.3989423*exp( -0.5*( ( log(x/(1-x))-MCMC_pop[k,p] )/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}

		f_m2 <- function(x){ 
			x*0.3989423*exp( -0.5*( (log(x/(1-x))-MCMC_pop[k,p])/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}
	
		tryCatch(
		{
			moment_1 <- integrate(f_m, lower=0, upper=1)$value
			moment_2 <- integrate(f_m2, lower=0, upper=1)$value
		}, error=function(e){ NULL }
		)

		MCMC_pop[k,p] <- moment_1
		MCMC_pop[k,N_glob_par+p] <- sqrt( moment_2 - moment_1^2 )
	}
}






TAB3[3,1:3]   = quantile( MCMC_pop[,1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[3,10:12] = quantile( MCMC_pop[,2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,13:15] = quantile( MCMC_pop[,3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,16:18] = quantile( MCMC_pop[,4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,19:21] = quantile( MCMC_pop[,5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[3,28:30] = quantile( MCMC_pop[,6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[3,37:39] = quantile( MCMC_pop[,7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,40:42] = quantile( MCMC_pop[,8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,43:45] = quantile( MCMC_pop[,9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,46:48] = quantile( MCMC_pop[,10], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[3,55:57] = quantile( MCMC_pop[,11], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,58:60] = quantile( MCMC_pop[,12], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[3,61:63] = quantile( MCMC_pop[,13], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )







###########################
###########################
##                       ##
##  PsA-TT-003_logELISA  ##
##                       ##
###########################
###########################

MCMC_pop_1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_003_mod2_rawELISA/Output/global_rep1.txt" )

MCMC_pop_2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_003_mod2_rawELISA/Output/global_rep2.txt" )


MCMC_pop = rbind( MCMC_pop_1[floor(0.5*nrow(MCMC_pop_1)):(nrow(MCMC_pop_1)-1),],
                  MCMC_pop_2[floor(0.5*nrow(MCMC_pop_2)):(nrow(MCMC_pop_2)-1),] )
colnames(MCMC_pop) = c("mean_beta_prim", "mean_t_short", "mean_t_long",  "mean_t_IgG", "mean_rho_prim", 
                       "sd_beta_prim", "sd_t_short", "sd_t_long", "sd_t_IgG", "sd_rho_prim",
     			     "sig_obs", "likelihood", "prior")



N_glob_par = 5



TAB3_median[4,1:3]   = quantile( exp(MCMC_pop[,1]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[4,10:12] = quantile( exp(MCMC_pop[,2]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[4,13:15] = quantile( exp(MCMC_pop[,3]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[4,16:18] = quantile( exp(MCMC_pop[,4]), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3_median[4,19:21] = quantile( exp(MCMC_pop[,5])/(1+exp(MCMC_pop[,5])), prob=c(0.5, 0.025, 0.975), na.rm=TRUE )




for(p in 1:N_glob_par)
{
	MCMC_pop[,N_glob_par+p] = 1/sqrt(MCMC_pop[,N_glob_par+p])
}

for(p in 1:4)
{
	MCMC_pop[,p]            = exp( MCMC_pop[,p] + 0.5*MCMC_pop[,N_glob_par+p]^2 )
	MCMC_pop[,N_glob_par+p] = sqrt( exp(MCMC_pop[,N_glob_par+p]^2) -1 )*MCMC_pop[,p]
}


for(p in 5)
{
	for(k in 1:nrow(MCMC_pop))
	{
		f_m <- function(x){ 
			0.3989423*exp( -0.5*( ( log(x/(1-x))-MCMC_pop[k,p] )/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}

		f_m2 <- function(x){ 
			x*0.3989423*exp( -0.5*( (log(x/(1-x))-MCMC_pop[k,p])/MCMC_pop[k,N_glob_par+p] )^2 )/( MCMC_pop[k,N_glob_par+p]*(1-x) )
		}
	
		tryCatch(
		{
			moment_1 <- integrate(f_m, lower=0, upper=1)$value
			moment_2 <- integrate(f_m2, lower=0, upper=1)$value
		}, error=function(e){ NULL }
		)

		MCMC_pop[k,p] <- moment_1
		MCMC_pop[k,N_glob_par+p] <- sqrt( moment_2 - moment_1^2 )
	}
}






TAB3[4,1:3]   = quantile( MCMC_pop[,1], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[4,10:12] = quantile( MCMC_pop[,2], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,13:15] = quantile( MCMC_pop[,3], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,16:18] = quantile( MCMC_pop[,4], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,19:21] = quantile( MCMC_pop[,5], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[4,28:30] = quantile( MCMC_pop[,6], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[4,37:39] = quantile( MCMC_pop[,7], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,40:42] = quantile( MCMC_pop[,8], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,43:45] = quantile( MCMC_pop[,9], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,46:48] = quantile( MCMC_pop[,10], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

TAB3[4,55:57] = quantile( MCMC_pop[,11], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,58:60] = quantile( MCMC_pop[,12], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )
TAB3[4,61:63] = quantile( MCMC_pop[,13], prob=c(0.5, 0.025, 0.975), na.rm=TRUE )





 