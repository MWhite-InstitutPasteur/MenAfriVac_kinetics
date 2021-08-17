##################################
##################################
##################################
###                            ###
###   Data: 002                ###
###                            ###
##################################
##################################
##################################


#####################################
## Read in data file

AB_data_002 = read.csv( "C:/U/Ab_dynamics/MenAfriVac/Data/Data_002/Data_002_process.csv" )

N_part_002 = nrow(AB_data_002)


AB_Cpp_002 = read.table( "C:/U/Ab_dynamics/MenAfriVac/Data/Data_002/Data_002_rawELISA_Cpp.txt" )



#################################
#################################
#################################
###                           ###
###   MCMC output: MenA_002   ###
###                           ###
#################################
#################################
#################################


#####################################
## Read in individual level MCMC output

MCMC_MenA_002_rep1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_002_mod2_2rho/Output/local_rep1.txt" )

MCMC_MenA_002_rep2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_002_mod2_2rho/Output/local_rep2.txt" )

MCMC_MenA_002 = rbind( MCMC_MenA_002_rep1[(0.5*nrow(MCMC_MenA_002_rep1)):nrow(MCMC_MenA_002_rep1),],
                       MCMC_MenA_002_rep2[(0.5*nrow(MCMC_MenA_002_rep2)):nrow(MCMC_MenA_002_rep2),] ) 



rm(MCMC_MenA_002_rep1)
rm(MCMC_MenA_002_rep2)



#####################################
#####################################
#####################################
###                               ###
###   MCMC output: logelisa_002   ###
###                               ###
#####################################
#####################################
#####################################


#####################################
## Read in individual level MCMC output

MCMC_logelisa_002_rep1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_002_mod2_rawELISA/Output/local_rep1.txt" )

MCMC_logelisa_002_rep2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_002_mod2_rawELISA/Output/local_rep2.txt" )

MCMC_logelisa_002 = rbind( MCMC_logelisa_002_rep1[(0.5*nrow(MCMC_logelisa_002_rep1)):nrow(MCMC_logelisa_002_rep1),],
                           MCMC_logelisa_002_rep2[(0.5*nrow(MCMC_logelisa_002_rep2)):nrow(MCMC_logelisa_002_rep2),] ) 


rm(MCMC_logelisa_002_rep1)
rm(MCMC_logelisa_002_rep2)




##################################
##################################
##################################
###                            ###
###   Data: 003                ###
###                            ###
##################################
##################################
##################################


#####################################
## Read in data file

AB_data_003 = read.csv( "C:/U/Ab_dynamics/MenAfriVac/Data/Data_003/Data_003_process.csv" )

N_part_003 = nrow(AB_data_003)


AB_Cpp_003 = read.table( "C:/U/Ab_dynamics/MenAfriVac/Data/Data_003/Data_003_rawELISA_Cpp.txt" )



#################################
#################################
#################################
###                           ###
###   MCMC output: MenA_003   ###
###                           ###
#################################
#################################
#################################


#####################################
## Read in individual level MCMC output

MCMC_MenA_003_rep1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_003_mod2/Output/local_rep1.txt" )

MCMC_MenA_003_rep2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/MenA_003_mod2/Output/local_rep2.txt" )

MCMC_MenA_003 = rbind( MCMC_MenA_003_rep1[(0.5*nrow(MCMC_MenA_003_rep1)):nrow(MCMC_MenA_003_rep1),],
                       MCMC_MenA_003_rep2[(0.5*nrow(MCMC_MenA_003_rep2)):nrow(MCMC_MenA_003_rep2),] ) 

rm(MCMC_MenA_003_rep1)
rm(MCMC_MenA_003_rep2)



#####################################
#####################################
#####################################
###                               ###
###   MCMC output: logelisa_003   ###
###                               ###
#####################################
#####################################
#####################################


#####################################
## Read in individual level MCMC output

MCMC_logelisa_003_rep1 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_003_mod2_rawELISA/Output/local_rep1.txt" )

MCMC_logelisa_003_rep2 = read.table( "C:/U/Ab_dynamics/MenAfriVac/AB_models_Gibbs/logelisa_003_mod2_rawELISA/Output/local_rep2.txt" )

MCMC_logelisa_003 = rbind( MCMC_logelisa_003_rep1[(0.5*nrow(MCMC_logelisa_003_rep1)):nrow(MCMC_logelisa_003_rep1),],
                           MCMC_logelisa_003_rep2[(0.5*nrow(MCMC_logelisa_003_rep2)):nrow(MCMC_logelisa_003_rep2),] ) 

rm(MCMC_logelisa_003_rep1)
rm(MCMC_logelisa_003_rep2)




##################################
##################################
##################################
###                            ###
###   Data frame               ###
###                            ###
##################################
##################################
##################################

AB_df_002 = matrix(NA, nrow=N_part_002, ncol=18)
colnames(AB_df_002) = c("ID", "country", "age", "sex", "prim_vac", "boost_vac", "SBA_base", "SBA_peak", "SBA_6m", "SBA_12m", "SBA_24m", "elisa_base", "elisa_peak", "elisa_6m", "elisa_12m", "elisa_24m", "height", "weight") 
AB_df_002 = as.data.frame(AB_df_002)

for(i in 1:N_part_002)
{
	AB_df_002[i,1] = paste( AB_data_002[i,1], "_002", sep="" )
}

AB_df_002[,c(2,3,4,5,6,7,8,12,13,17,18)] = AB_data_002[,c(2,4,5,10,11,30,32,21,23,6,7)]

AB_df_002[,3] = AB_df_002[,3]/12

for(i in 1:N_part_002)
{
	if( AB_df_002[i,5]=="Hib-TT" && AB_df_002[i,6]=="PsA-TT" )
	{
		AB_df_002[i,3]  = AB_df_002[i,3] + (AB_data_002[i,17] - AB_data_002[i,14])/365  
		AB_df_002[i,8]  = AB_data_002[i,35]		
		AB_df_002[i,13] = AB_data_002[i,26]	
	}
}	



AB_df_003 = matrix(NA, nrow=N_part_003, ncol=18)
colnames(AB_df_003) = c("ID", "country", "age", "sex", "prim_vac", "boost_vac", "SBA_base", "SBA_peak", "SBA_6m", "SBA_12m", "SBA_24m", "elisa_base", "elisa_peak", "elisa_6m", "elisa_12m", "elisa_24m", "height", "weight") 
AB_df_003 = as.data.frame(AB_df_003)

for(i in 1:N_part_003)
{
	AB_df_003[i,1] = paste( AB_data_003[i,1], "_003", sep="" )
}

AB_df_003[,c(2,3,4,5,7,8,12,13,17,18)] = AB_data_003[,c(2,4,5,9,20,22,15,17,7,8)]


AB_df_all = rbind( AB_df_002, AB_df_003 )


for(i in 1:N_part_002)
{
	if( AB_df_002[i,5]=="Hib-TT" && AB_df_002[i,6]=="PsA-TT" )
	{
		AB_df_002[i,3]  = AB_df_002[i,3] + (AB_data_002[i,17] - AB_data_002[i,14])/365  
		AB_df_002[i,8]  = AB_data_002[i,35]		
		AB_df_002[i,13] = AB_data_002[i,26]	
	}
}	

AB_df_all = rbind( AB_df_002, AB_df_003 )



#########################################
#########################################
##                                     ##
##  #     #  ####  ####   ##### ##     ##
##  ##   ## ##  ## ## ##  ##    ##     ## 
##  ####### ##  ## ##  ## ####  ##     ##
##  ## # ## ##  ## ## ##  ##    ##     ##
##  ##   ##  ####  ####   ##### #####  ##
##                                     ##
#########################################
#########################################

x_1 = 0.1
y_1 = 2

x_2 = 2
y_2 = 128

ELISA_axis_transform = function( AB_elisa )
{
	exp( log(y_1) + ( (log(y_2) - log(y_1) )/( log(x_2) - log(x_1) ))*( log(AB_elisa) - log(x_1) ) )
}


#######################################
#######################################
#######################################
###                                 ###
###   Model predict: MenA_002       ###
###                                 ###
#######################################
#######################################
#######################################

MenA_002_pred <- function( n )
{
	###################################
	## Data for participant n

	tt = as.numeric(AB_Cpp_002[n,11:17])
	AB = as.numeric(AB_Cpp_002[n,18:24])

	tt = tt[ which(AB > -0.5) ]
	AB = AB[ which(AB > -0.5) ]


	tt  = tt - tt[1]

	if( is.na(AB_data_002[n,12])==FALSE )
	{
		t_vac_boost = AB_data_002[n,15] - AB_data_002[n,12]
	}else{
		t_vac_boost = AB_data_002[n,15] - median(AB_data_002[,12],na.rm=TRUE)
	}

	###################################
	## Model prediction for participant i

	N_loc_par = 10
	MCMC_ind_n = MCMC_MenA_002[,((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]


	AB_mod = matrix(NA, nrow=nrow(MCMC_MenA_002), ncol=length(tt_plot))
	
	for(i in 1:nrow(AB_mod))
	{
		beta_prim_PsA  = MCMC_ind_n[i,1]
		beta_boost_PsA = MCMC_ind_n[i,2]
		beta_Hib       = MCMC_ind_n[i,3]
		d_cs           = MCMC_ind_n[i,4]
		d_cl           = MCMC_ind_n[i,5]
		d_a            = MCMC_ind_n[i,6]
		rho_prim_PsA   = MCMC_ind_n[i,7]
		rho_boost_PsA  = MCMC_ind_n[i,8]
		rho_Hib        = MCMC_ind_n[i,9]
		Ab_0           = MCMC_ind_n[i,10]

		r_cs = log(2)/d_cs
		r_cl = log(2)/d_cl
		r_a  = log(2)/d_a



		if( AB_Cpp_002[n,7] == 1 )
		{
			AB_tt = Ab_0*exp(-r_cl*tt_plot) 
	
			AB_tt = AB_tt + beta_prim_PsA*(       (rho_prim_PsA/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                          ((1 - rho_prim_PsA)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )
		}

		if( AB_Cpp_002[n,7] == 2 )
		{
			AB_tt = Ab_0*exp(-r_cl*(tt_plot+t_vac_boost)) 

			AB_tt = AB_tt + beta_Hib*(       (rho_Hib/(r_a - r_cs))*( exp(-r_cs*(tt_plot+t_vac_boost)) - exp(-r_a*(tt_plot+t_vac_boost)) ) +
				                     ((1 - rho_Hib)/(r_a - r_cl))*( exp(-r_cl*(tt_plot+t_vac_boost)) - exp(-r_a*(tt_plot+t_vac_boost))) )
	
			AB_tt = AB_tt + beta_prim_PsA*(       (rho_prim_PsA/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                          ((1 - rho_prim_PsA)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )
		}



		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975) )
	}


	###################
	## Return output

	AB_quant[2,]
}




#######################################
#######################################
#######################################
###                                 ###
###   Model predict: logelisa_002   ###
###                                 ###
#######################################
#######################################
#######################################

logelisa_002_pred <- function( n )
{
	###################################
	## Data for participant n


	tt = as.numeric(AB_Cpp_002[n,11:17])
	AB = as.numeric(AB_Cpp_002[n,25:31])

	tt = tt[ which(AB > -0.5) ]
	AB = AB[ which(AB > -0.5) ]



	tt  = tt - tt[1]

	if( is.na(AB_data_002[n,12])==FALSE )
	{
		t_vac_boost = AB_data_002[n,15] - AB_data_002[n,12]
	}else{
		t_vac_boost = AB_data_002[n,15] - median(AB_data_002[,12],na.rm=TRUE)
	}

	###################################
	## Model prediction for participant i

	N_loc_par = 8
	MCMC_ind_n = MCMC_logelisa_002[,((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]

	AB_mod = matrix(NA, nrow=nrow(MCMC_logelisa_002), ncol=length(tt_plot))
	
	for(i in 1:nrow(AB_mod))
	{
		beta_prim  = MCMC_ind_n[i,1]
		beta_boost = MCMC_ind_n[i,2]
		d_cs       = MCMC_ind_n[i,3]
		d_cl       = MCMC_ind_n[i,4]
		d_a        = MCMC_ind_n[i,5]
		rho_prim   = MCMC_ind_n[i,6]
		rho_boost  = MCMC_ind_n[i,7]
		Ab_0       = MCMC_ind_n[i,8]

		r_cs = log(2)/d_cs
		r_cl = log(2)/d_cl
		r_a  = log(2)/d_a


		if( AB_Cpp_002[n,7] == 1 )
		{
			AB_tt = Ab_0*exp(-r_cl*tt_plot) 
	
			AB_tt = AB_tt + beta_prim*(       (rho_prim/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                      ((1 - rho_prim)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )
		}

		if( AB_Cpp_002[n,7] == 2 )
		{
			AB_tt = Ab_0*exp(-r_cl*(tt_plot+t_vac_boost)) 

			AB_tt = AB_tt + beta_prim*(       (rho_prim/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                      ((1 - rho_prim)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )
		}



		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975) )
	}

	
	###################
	## Return output

	AB_quant[2,]
}





#######################################
#######################################
#######################################
###                                 ###
###   Model predict: MenA_003       ###
###                                 ###
#######################################
#######################################
#######################################

MenA_003_pred <- function( n )
{
	###################################
	## Data for participant n

	tt = as.numeric(AB_Cpp_003[n,11:14])
	AB = as.numeric(AB_Cpp_003[n,15:18])

	tt = tt[ which(AB > -0.5) ]
	AB = AB[ which(AB > -0.5) ]


	tt  = tt - tt[1]


	###################################
	## Model prediction for participant i

	N_loc_par = 6
	MCMC_ind_n = MCMC_MenA_003[,((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]

	AB_mod = matrix(NA, nrow=nrow(MCMC_ind_n), ncol=length(tt_plot))
	
	for(i in 1:nrow(AB_mod))
	{
		beta_prim  = MCMC_ind_n[i,1]
		d_cs       = MCMC_ind_n[i,2]
		d_cl       = MCMC_ind_n[i,3]
		d_a        = MCMC_ind_n[i,4]
		rho_prim   = MCMC_ind_n[i,5]
		AB_0       = MCMC_ind_n[i,6]

		r_cs = log(2)/d_cs
		r_cl = log(2)/d_cl
		r_a  = log(2)/d_a

		AB_tt = AB_0*exp(-tt_plot*r_cl) +
                    beta_prim*( (rho_prim/(r_a-r_cs))*( exp(-tt_plot*r_cs) - exp(-tt_plot*r_a)) +
                                ((1-rho_prim)/(r_a-r_cl))*(exp(-tt_plot*r_cl) - exp(-tt_plot*r_a)) )

		AB_tt[which(tt_plot<0)] = AB_0

		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975) )
	}


	###################
	## Return output

	AB_quant[2,]
}



#######################################
#######################################
#######################################
###                                 ###
###   Model predict: logelisa_003   ###
###                                 ###
#######################################
#######################################
#######################################

logelisa_003_pred <- function( n )
{
	###################################
	## Data for participant n

	tt = as.numeric(AB_Cpp_003[n,11:14])
	AB = as.numeric(AB_Cpp_003[n,19:22])

	tt = tt[ which(AB > -0.5) ]
	AB = AB[ which(AB > -0.5) ]


	tt  = tt - tt[1]


	###################################
	## Model prediction for participant i

	N_loc_par = 6
	MCMC_ind_n = MCMC_logelisa_003[,((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]

	AB_mod = matrix(NA, nrow=nrow(MCMC_ind_n), ncol=length(tt_plot))
	
	for(i in 1:nrow(AB_mod))
	{
		beta_prim  = MCMC_ind_n[i,1]
		d_cs       = MCMC_ind_n[i,2]
		d_cl       = MCMC_ind_n[i,3]
		d_a        = MCMC_ind_n[i,4]
		rho_prim   = MCMC_ind_n[i,5]
		AB_0       = MCMC_ind_n[i,6]

		r_cs = log(2)/d_cs
		r_cl = log(2)/d_cl
		r_a  = log(2)/d_a

		AB_tt = AB_0*exp(-tt_plot*r_cl) +
                    beta_prim*( (rho_prim/(r_a-r_cs))*( exp(-tt_plot*r_cs) - exp(-tt_plot*r_a)) +
                                ((1-rho_prim)/(r_a-r_cl))*(exp(-tt_plot*r_cl) - exp(-tt_plot*r_a)) )

		AB_tt[which(tt_plot<0)] = AB_0

		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975) )
	}


	###################
	## Return output

	AB_quant[2,]
}



###########################
###########################
##                       ##
##   ####  #### #     #  ##
##  ##      ##  ##   ##  ##
##   ####   ##  #######  ##
##      ##  ##  ## # ##  ##
##   ####  #### ##   ##  ##
##                       ##
###########################
###########################

tt_plot = c(0.5*365, 1*365, 2*365)

for(i in 1:N_part_002)
{
	AB_df_all[i,9:11] = MenA_002_pred( i )

	AB_df_all[i,14:16] = logelisa_002_pred( i )
}


for(i in 1:N_part_003)
{
	AB_df_all[N_part_002+i,9:11] = MenA_003_pred( i )

	AB_df_all[N_part_002+i,14:16] = logelisa_003_pred( i )
}


site_col = rep("forestgreen", N_part_002 + N_part_003)

site_col[which(AB_df_all$country=="Gambia")]  = "gold"
site_col[which(AB_df_all$country=="Senegal")] = "royalblue"


site_shape = rep(17, N_part_002 + N_part_003)
site_shape[1:N_part_002] = 19



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


AB_df_SBA = cbind( AB_df_all, c( rep("two", N_part_002), rep("three", N_part_003) ) )
colnames(AB_df_SBA)[19] = "trial"

AB_df_SBA = AB_df_SBA[which( AB_df_all$SBA_12m <= AB_df_all$SBA_peak ),]



#############################################
## Relevel factors

AB_df_SBA$trial <- relevel( AB_df_SBA$trial, ref="two")

AB_df_SBA$country <- relevel( AB_df_SBA$country, ref="Gambia")

AB_df_SBA$sex <- relevel( AB_df_SBA$sex, ref="female")





###################
## Percentage reduction
## Linear model SBA_12m

AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ trial*age + country + log10(SBA_base) + log10(SBA_peak),  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)








AB_df_elisa = cbind( AB_df_all, c( rep("two", N_part_002), rep("three", N_part_003) ) )
colnames(AB_df_elisa)[19] = "trial"

AB_df_elisa = AB_df_elisa[which( AB_df_all$elisa_12m <= AB_df_all$elisa_peak ),]



#############################################
## Relevel factors

AB_df_elisa$trial <- relevel( AB_df_elisa$trial, ref="two")

AB_df_elisa$country <- relevel( AB_df_elisa$country, ref="Gambia")

AB_df_elisa$sex <- relevel( AB_df_elisa$sex, ref="female")



###################
## Percentage reduction
## Linear model elisa_12m

AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ trial*age + country + log10(elisa_base) + log10(elisa_peak),  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)







###################
## 12 month levels
## Linear model SBA_12m

AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ trial*age + country + log10(SBA_base) + log10(SBA_peak),  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)





###################
## 12 month levels
## Linear model elisa_12m

AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ trial*age + country  + log10(elisa_base) + log10(elisa_peak),  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)















###########################################################
###########################################################
##                                                       ##
##   ####  ##  ## #####   ######  ####  #####       ##   ##
##  ##     ##  ## ##  ##    ##   ##  ## ##  ##     ##    ##
##   ####  ##  ## #####     ##   ###### #####     ####   ##
##      ## ##  ## ##        ##   ##  ## ##  ##   ##  ##  ##
##   ####   ####  ##        ##   ##  ## #####     ####   ##
##                                                       ##
###########################################################
###########################################################


AB_df_SBA = cbind( AB_df_all, c( rep("two", N_part_002), rep("three", N_part_003) ) )
colnames(AB_df_SBA)[19] = "trial"

AB_df_SBA = AB_df_SBA[which( AB_df_all$SBA_12m <= AB_df_all$SBA_peak ),]



#############################################
## Relevel factors

AB_df_SBA$trial <- relevel( AB_df_SBA$trial, ref="two")

AB_df_SBA$country <- relevel( AB_df_SBA$country, ref="Gambia")

AB_df_SBA$sex <- relevel( AB_df_SBA$sex, ref="female")





###################
## Percentage reduction
## Linear model SBA_12m

AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ trial*age + country + sex + height + weight + log10(SBA_base) + log10(SBA_peak),  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)








AB_df_elisa = cbind( AB_df_all, c( rep("two", N_part_002), rep("three", N_part_003) ) )
colnames(AB_df_elisa)[19] = "trial"

AB_df_elisa = AB_df_elisa[which( AB_df_all$elisa_12m <= AB_df_all$elisa_peak ),]



#############################################
## Relevel factors

AB_df_elisa$trial <- relevel( AB_df_elisa$trial, ref="two")

AB_df_elisa$country <- relevel( AB_df_elisa$country, ref="Gambia")

AB_df_elisa$sex <- relevel( AB_df_elisa$sex, ref="female")



###################
## Percentage reduction
## Linear model elisa_12m

AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ trial*age + country + sex + height + weight + log10(elisa_base) + log10(elisa_peak),  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)







###################
## 12 month levels
## Linear model SBA_12m

AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ trial*age + country + sex + height + weight + log10(SBA_base) + log10(SBA_peak),  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)





###################
## 12 month levels
## Linear model elisa_12m

AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ trial*age + country + sex + height + weight  + log10(elisa_base) + log10(elisa_peak),  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)









###########################################################
###########################################################
##                                                       ##
##   ####  ##  ## #####   ######  ####  #####    #####   ##
##  ##     ##  ## ##  ##    ##   ##  ## ##  ##   ##      ##
##   ####  ##  ## #####     ##   ###### #####     ####   ##
##      ## ##  ## ##        ##   ##  ## ##  ##       ##  ##
##   ####   ####  ##        ##   ##  ## #####     ####   ##
##                                                       ##
###########################################################
###########################################################

###########################################################
##  SBA 1 year red



AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ trial,  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)





AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ country,  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)





AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ age,  
                         data=AB_df_SBA[which(AB_df_SBA$trial=="two"),], family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)







AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ age,  
                         data=AB_df_SBA[which(AB_df_SBA$trial=="three"),], family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)






AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ sex ,  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)





AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ height,  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)





AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ weight ,  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)





AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~ log10(SBA_base),  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)





AB_mod_SBA_12mred <- glm( (1 - SBA_12m/SBA_peak) ~  log10(SBA_peak),  
                         data=AB_df_SBA, family="binomial")

summary(AB_mod_SBA_12mred)



cbind(
	summary(AB_mod_SBA_12mred)$coef[,1], 
	summary(AB_mod_SBA_12mred)$coef[,1] - 1.96*summary(AB_mod_SBA_12mred)$coef[,2],   
	summary(AB_mod_SBA_12mred)$coef[,1] + 1.96*summary(AB_mod_SBA_12mred)$coef[,2]
)



###########################################################
##  SBA 1 year level



AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ trial,  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)






AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ country ,  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)








AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ age,  
                         data=AB_df_SBA[which(AB_df_SBA$trial=="two"),])

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)





AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ age,  
                         data=AB_df_SBA[which(AB_df_SBA$trial=="three"),])

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)








AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ sex ,  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)







AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ height ,  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)







AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ weight,  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)







AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ log10(SBA_base),  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)






AB_mod_SBA_12m <- lm( log10(SBA_12m) ~ log10(SBA_peak),  
                         data=AB_df_SBA)

summary(AB_mod_SBA_12m)



cbind(
	summary(AB_mod_SBA_12m)$coef[,1], 
	summary(AB_mod_SBA_12m)$coef[,1] - 1.96*summary(AB_mod_SBA_12m)$coef[,2],   
	summary(AB_mod_SBA_12m)$coef[,1] + 1.96*summary(AB_mod_SBA_12m)$coef[,2]
)




###########################################################
##  Elisa reduction after 1 year



AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ trial,  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)





AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ country ,  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)






AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ age,  
                            data=AB_df_elisa[which(AB_df_elisa$trial=="two"),], family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)




AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ age,  
                            data=AB_df_elisa[which(AB_df_elisa$trial=="three"),], family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)






AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ sex ,  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)





AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ height ,  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)





AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ weight ,  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)






AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~ log10(elisa_base),  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)






AB_mod_elisa_12mred <- glm( (1 - elisa_12m/elisa_peak) ~  log10(elisa_peak),  
                            data=AB_df_elisa, family="binomial")

summary(AB_mod_elisa_12mred)



cbind(
	summary(AB_mod_elisa_12mred)$coef[,1], 
	summary(AB_mod_elisa_12mred)$coef[,1] - 1.96*summary(AB_mod_elisa_12mred)$coef[,2],   
	summary(AB_mod_elisa_12mred)$coef[,1] + 1.96*summary(AB_mod_elisa_12mred)$coef[,2]
)






###########################################################
##  Elisa level after 1 year


AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ trial,  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)






AB_mod_elisa_12m <- lm( log10(elisa_12m) ~  country ,  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)






AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ age,  
                         data=AB_df_elisa[which(AB_df_elisa$trial=="two"),])

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)





AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ age,  
                         data=AB_df_elisa[which(AB_df_elisa$trial=="three"),])

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)







AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ sex ,  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)






AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ height ,  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)






AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ weight,  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)






AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ log10(elisa_base),  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)






AB_mod_elisa_12m <- lm( log10(elisa_12m) ~ log10(elisa_peak),  
                         data=AB_df_elisa)

summary(AB_mod_elisa_12m)



cbind(
	summary(AB_mod_elisa_12m)$coef[,1], 
	summary(AB_mod_elisa_12m)$coef[,1] - 1.96*summary(AB_mod_elisa_12m)$coef[,2],   
	summary(AB_mod_elisa_12m)$coef[,1] + 1.96*summary(AB_mod_elisa_12m)$coef[,2]
)


