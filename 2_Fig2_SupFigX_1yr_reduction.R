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

AB_df_002 = matrix(NA, nrow=N_part_002, ncol=16)
colnames(AB_df_002) = c("ID", "country", "age", "sex", "prim_vac", "boost_vac", "SBA_base", "SBA_peak", "SBA_6m", "SBA_12m", "SBA_24m", "elisa_base", "elisa_peak", "elisa_6m", "elisa_12m", "elisa_24m") 
AB_df_002 = as.data.frame(AB_df_002)

for(i in 1:N_part_002)
{
	AB_df_002[i,1] = paste( AB_data_002[i,1], "_002", sep="" )
}

AB_df_002[,c(2,3,4,5,6,7,8,12,13)] = AB_data_002[,c(2,4,5,10,11,30,32,21,23)]

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



AB_df_003 = matrix(NA, nrow=N_part_003, ncol=16)
colnames(AB_df_003) = c("ID", "country", "age", "sex", "prim_vac", "boost_vac", "SBA_base", "SBA_peak", "SBA_6m", "SBA_12m", "SBA_24m", "elisa_base", "elisa_peak", "elisa_6m", "elisa_12m", "elisa_24m") 
AB_df_003 = as.data.frame(AB_df_003)

for(i in 1:N_part_003)
{
	AB_df_003[i,1] = paste( AB_data_003[i,1], "_003", sep="" )
}

AB_df_003[,c(2,3,4,5,7,8,12,13)] = AB_data_003[,c(2,4,5,9,20,22,15,17)]


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


##################################
##################################
##                              ##
##  ##### ####  ####     ####   ##
##  ##     ##  ##       ##  ##  ##
##  ####   ##  ## ###      ##   ## 
##  ##     ##  ##  ##     ##    ## 
##  ##    ####  ####     #####  ##
##                              ##
################################## 
##################################


main.size = 0.9
axis.size = 0.75
lab.size  = 1.0


tiff( file="2_Fig2_1yr_reduction.tif", width=24, height=18, units="cm", res=500)

lay.mat = rbind( c( 1,  2,  3  ),
                 c( 4,  5,  6  ),
                 c( 7,  7,  7  ) )

layout(lay.mat, heights=c(10,10,1))
layout.show(7)



par(mar = c(2.5,2.5,1,1.75))
par(mgp=c(1.5,0.5,0))



############################
##                        ## 
##  PANEL A               ##
##  SBA reduction vs age  ##
##                        ##
############################


line_seq_x <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

line_seq_y <- c(0.5, 1, 2, 5, 10, 20, 50)



plot( x=1e10, y=1e10,  
xlim=c(0.5,50), ylim=c(-0.01,1.01), log="x", 
xlab="age (years)", ylab="percentage reduction in first year",
main="(A) SBA titer reduction vs age",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$age, y=1 - AB_df_all$SBA_12m/AB_df_all$SBA_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1, at=c(0.5, 1, 2, 5, 10, 20, 50), labels=c(0.5, 1, 2, 5, 10, 20, 50), cex.axis=axis.size)

axis(2,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%",  "40%", "60%", "80%", "100%" ), 
cex.axis=axis.size)



AB_red <- 1 - AB_df_all$SBA_12m/AB_df_all$SBA_peak
age_x  <- log(AB_df_all$age)

AB_red[which(AB_red < 0)] <- NA

age_seq <- seq(from=log(0.5), to=log(50), by=0.01)

logit_mod <- glm( AB_red ~ age_x, family="binomial" )

logit_mod_fit <- predict.glm( logit_mod, data.frame(age_x=age_seq), interval="prediction", se.fit=TRUE)

logit_mod_med = exp(logit_mod_fit$fit)/(1 + exp(logit_mod_fit$fit))

logit_mod_low = exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit))

logit_mod_high = exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit))



points(x=exp(age_seq), y=logit_mod_med,  type='l', lwd=2 )

points(x=exp(age_seq), y=logit_mod_low,  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(age_seq), y=logit_mod_high, type='l', 
lty="longdash", col="grey", lwd=2 )




#################################
##                             ## 
##  PANEL B                    ##
##  SBA reduction vs SBA_base  ##
##                             ##
#################################


line_seq_x <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
line_seq_y <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

plot( x=1e10, y=1e10,  
xlim=c(2^(0),2^19), ylim=c(-0.01,1.01), log="x", 
xlab="pre-MenAfriVac SBA titer", ylab="percentage reduction in first year",
main="(B) SBA titer reduction vs pre-MenAfriVac SBA titer",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$SBA_base, y=1 - AB_df_all$SBA_12m/AB_df_all$SBA_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)

axis(2,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%",  "40%", "60%", "80%", "100%" ), 
cex.axis=axis.size)





AB_red <- 1 - AB_df_all$SBA_12m/AB_df_all$SBA_peak
l_AB_base <- log(AB_df_all$SBA_base)

AB_red[which(AB_red < 0)] <- NA

l_AB_base_seq <- seq(from=log(0.5), to=log(2^20), by=0.01)

logit_mod <- glm( AB_red ~ l_AB_base, family="binomial" )

logit_mod_fit <- predict.glm( logit_mod, data.frame(l_AB_base=l_AB_base_seq), interval="prediction", se.fit=TRUE)

logit_mod_med = exp(logit_mod_fit$fit)/(1 + exp(logit_mod_fit$fit))

logit_mod_low = exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit))

logit_mod_high = exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit))



points(x=exp(l_AB_base_seq), y=logit_mod_med,  type='l', lwd=2 )

points(x=exp(l_AB_base_seq), y=logit_mod_low,  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_base_seq), y=logit_mod_high, type='l', 
lty="longdash", col="grey", lwd=2 )








#################################
##                             ## 
##  PANEL C                    ##
##  SBA reduction vs SBA_peak  ##
##                             ##
#################################


line_seq_x <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
line_seq_y <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

plot( x=1e10, y=1e10,  
xlim=c(2^(0),2^19), ylim=c(-0.01,1.01), log="x", 
xlab="post-MenAfriVac SBA titer", ylab="percentage reduction in first year",
main="(C) SBA titer reduction vs post-MenAfriVac SBA titer",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$SBA_peak, y=1 - AB_df_all$SBA_12m/AB_df_all$SBA_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)

axis(2,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%",  "40%", "60%", "80%", "100%" ), 
cex.axis=axis.size)






AB_red <- 1 - AB_df_all$SBA_12m/AB_df_all$SBA_peak
l_AB_peak <- log(AB_df_all$SBA_peak)

AB_red[which(AB_red < 0)] <- NA

l_AB_peak_seq <- seq(from=log(0.5), to=log(2^20), by=0.01)

logit_mod <- glm( AB_red ~ l_AB_peak, family="binomial" )

logit_mod_fit <- predict.glm( logit_mod, data.frame(l_AB_peak=l_AB_peak_seq), interval="prediction", se.fit=TRUE)

logit_mod_med = exp(logit_mod_fit$fit)/(1 + exp(logit_mod_fit$fit))

logit_mod_low = exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit))

logit_mod_high = exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit))



points(x=exp(l_AB_peak_seq), y=logit_mod_med,  type='l', lwd=2 )

points(x=exp(l_AB_peak_seq), y=logit_mod_low,  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_peak_seq), y=logit_mod_high, type='l', 
lty="longdash", col="grey", lwd=2 )






##############################
##                          ## 
##  PANEL D                 ##
##  ELISA reduction vs age  ##
##                          ##
##############################


line_seq_x <- c(0, 0.2, 0.4, 0.6, 0.8, 1)

line_seq_y <- c(0.5, 1, 2, 5, 10, 20, 50)



plot( x=1e10, y=1e10,  
xlim=c(0.5,50), ylim=c(-0.01,1.01), log="x", 
xlab="age (years)", ylab="percentage reduction in first year",
main="(D) IgG ELISA reduction vs age",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$age, y=1 - AB_df_all$elisa_12m/AB_df_all$elisa_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1, at=c(0.5, 1, 2, 5, 10, 20, 50), labels=c(0.5, 1, 2, 5, 10, 20, 50), cex.axis=axis.size)

axis(2,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%",  "40%", "60%", "80%", "100%" ), 
cex.axis=axis.size)





AB_red <- 1 - AB_df_all$elisa_12m/AB_df_all$elisa_peak
age_x  <- log(AB_df_all$age)

AB_red[which(AB_red < 0)] <- NA

age_seq <- seq(from=log(0.5), to=log(50), by=0.01)

logit_mod <- glm( AB_red ~ age_x, family="binomial" )

logit_mod_fit <- predict.glm( logit_mod, data.frame(age_x=age_seq), interval="prediction", se.fit=TRUE)

logit_mod_med = exp(logit_mod_fit$fit)/(1 + exp(logit_mod_fit$fit))

logit_mod_low = exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit))

logit_mod_high = exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit))



points(x=exp(age_seq), y=logit_mod_med,  type='l', lwd=2 )

points(x=exp(age_seq), y=logit_mod_low,  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(age_seq), y=logit_mod_high, type='l', 
lty="longdash", col="grey", lwd=2 )





#####################################
##                                 ## 
##  PANEL E                        ##
##  ELISA reduction vs ELISA_base  ##
##                                 ##
#####################################


line_seq_x <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
line_seq_y <- c(0.1, 1, 10, 100, 1000, 3000)

plot( x=1e10, y=1e10,  
xlim=c(0.08,3000), ylim=c(-0.01,1.01), log="x", 
xlab=expression(paste( "pre-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
ylab="percentage reduction in first year",
main="(E) IgG ELISA reduction vs pre-MenAfriVac IgG ELISA",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$elisa_base, y=1 - AB_df_all$elisa_12m/AB_df_all$elisa_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)

axis(2,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%",  "40%", "60%", "80%", "100%" ), 
cex.axis=axis.size)




AB_red <- 1 - AB_df_all$elisa_12m/AB_df_all$elisa_peak
l_AB_base <- log(AB_df_all$elisa_base)

AB_red[which(AB_red < 0)] <- NA

l_AB_base_seq <- seq(from=log(0.05), to=log(2^20), by=0.01)

logit_mod <- glm( AB_red ~ l_AB_base, family="binomial" )

logit_mod_fit <- predict.glm( logit_mod, data.frame(l_AB_base=l_AB_base_seq), interval="prediction", se.fit=TRUE)

logit_mod_med = exp(logit_mod_fit$fit)/(1 + exp(logit_mod_fit$fit))

logit_mod_low = exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit))

logit_mod_high = exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit))



points(x=exp(l_AB_base_seq), y=logit_mod_med,  type='l', lwd=2 )

points(x=exp(l_AB_base_seq), y=logit_mod_low,  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_base_seq), y=logit_mod_high, type='l', 
lty="longdash", col="grey", lwd=2 )






###################################
##                               ## 
##  PANEL F                      ##
##  ELISA reduction vs SBA_peak  ##
##                               ##
###################################


line_seq_x <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
line_seq_y <- c(0.1, 1, 10, 100, 1000, 3000)


plot( x=1e10, y=1e10,  
xlim=c(0.08,3000), ylim=c(-0.01,1.01), log="x", 
xlab=expression(paste( "post-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
ylab="percentage reduction in first year",
main="(F) IgG ELISA reduction vs post-MenAfriVac IgG ELISA",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$elisa_peak, y=1 - AB_df_all$elisa_12m/AB_df_all$elisa_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)

axis(2,  at = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
     labels = c( "0%", "20%",  "40%", "60%", "80%", "100%" ), 
cex.axis=axis.size)





AB_red <- 1 - AB_df_all$elisa_12m/AB_df_all$elisa_peak
l_AB_peak <- log(AB_df_all$elisa_peak)

AB_red[which(AB_red < 0)] <- NA

l_AB_peak_seq <- seq(from=log(0.05), to=log(2^20), by=0.01)

logit_mod <- glm( AB_red ~ l_AB_peak, family="binomial" )

logit_mod_fit <- predict.glm( logit_mod, data.frame(l_AB_peak=l_AB_peak_seq), interval="prediction", se.fit=TRUE)

logit_mod_med = exp(logit_mod_fit$fit)/(1 + exp(logit_mod_fit$fit))

logit_mod_low = exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit - 1.96*logit_mod_fit$se.fit))

logit_mod_high = exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit)/(1 + exp(logit_mod_fit$fit + 1.96*logit_mod_fit$se.fit))



points(x=exp(l_AB_peak_seq), y=logit_mod_med,  type='l', lwd=2 )

points(x=exp(l_AB_peak_seq), y=logit_mod_low,  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_peak_seq), y=logit_mod_high, type='l', 
lty="longdash", col="grey", lwd=2 )





###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("Gambia (PsA-TT-002)", "Gambia (PsA-TT-003)", "Mali (PsA-TT-002)", "Mali (PsA-TT-003)", "Senegal (PsA-TT-003)"), 
       col = c("gold", "gold", "forestgreen", "forestgreen", "royalblue"), 
       pch=c(19, 17, 19, 17, 17),
       ncol=5, cex=1.25, bty="n" )


dev.off()

	




#########################################################
#########################################################
##                                                     ##
##   ####  ##  ## #####    ##### ####  ####       ##   ##
##  ##     ##  ## ##  ##   ##     ##  ##         ##    ##  
##   ####  ##  ## #####    ####   ##  ## ###    ####   ##
##      ## ##  ## ##       ##     ##  ##  ##   ##  ##  ## 
##   ####   ####  ##       ##    ####  ####     ####   ##
##                                                     ##
#########################################################
#########################################################


main.size = 0.9
axis.size = 0.75
lab.size  = 1.0


tiff( file="2_SupFig14_1yr_Ab.tif", width=24, height=18, units="cm", res=500)

lay.mat = rbind( c( 1,  2,  3  ),
                 c( 4,  5,  6  ),
                 c( 7,  7,  7  ) )

layout(lay.mat, heights=c(10,10,1))
layout.show(7)



par(mar = c(2.5,2.5,1,1.75))
par(mgp=c(1.5,0.5,0))



############################
##                        ## 
##  PANEL A               ##
##  SBA reduction vs age  ##
##                        ##
############################


line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

line_seq_y <- c(0.5, 1, 2, 5, 10, 20, 50)



plot( x=1e10, y=1e10,  
xlim=c(0.5,50), ylim=c(2^(-0.25),2^19), log="xy", 
xlab="age (years)", 
ylab="SBA titer after 1 year",
main="(A) SBA titer vs age",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$age, y=AB_df_all$SBA_12m, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1, at=c(0.5, 1, 2, 5, 10, 20, 50), labels=c(0.5, 1, 2, 5, 10, 20, 50), cex.axis=axis.size)

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$SBA_12m)
age_x  <- log(AB_df_all$age)

AB_red[which(AB_red < 0)] <- NA

age_seq <- seq(from=log(0.5), to=log(50), by=0.01)

lin_mod <- lm( AB_red ~ age_x )

lin_mod_fit <- predict( lin_mod, data.frame(age_x=age_seq), interval="confidence")

points(x=exp(age_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )



#################################
##                             ## 
##  PANEL B                    ##
##  SBA reduction vs SBA_base  ##
##                             ##
#################################


line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)
line_seq_y <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

plot( x=1e10, y=1e10,  
xlim=c(2^(-0.25),2^19), ylim=c(2^(-0.25),2^19), log="xy", 
xlab="pre-MenAfriVac SBA titer", 
ylab="SBA titer after 1 year",
main="(B) SBA titer vs pre-MenAfriVac SBA titer",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$SBA_base, y=AB_df_all$SBA_12m, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$SBA_12m)
l_AB_base <- log(AB_df_all$SBA_base)

AB_red[which(AB_red < 0)] <- NA

l_AB_base_seq <- seq(from=log(0.5), to=log(2^20), by=0.01)

lin_mod <- lm( AB_red ~ l_AB_base )

lin_mod_fit <- predict( lin_mod, data.frame(l_AB_base=l_AB_base_seq), interval="confidence")

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )






#################################
##                             ## 
##  PANEL C                    ##
##  SBA reduction vs SBA_peak  ##
##                             ##
#################################


line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)
line_seq_y <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

plot( x=1e10, y=1e10,  
xlim=c(2^(-0.25),2^19), ylim=c(2^(-0.25),2^19), log="xy", 
xlab="post-MenAfriVac SBA titer",
ylab="SBA titer after 1 year",
main="(C) SBA titer vs post-MenAfriVac SBA titer",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$SBA_peak, y=AB_df_all$SBA_12m, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$SBA_12m)
l_AB_peak <- log(AB_df_all$SBA_peak)

AB_red[which(AB_red < 0)] <- NA

l_AB_peak_seq <- seq(from=log(0.5), to=log(2^20), by=0.01)

lin_mod <- lm( AB_red ~ l_AB_peak )

lin_mod_fit <- predict( lin_mod, data.frame(l_AB_peak=l_AB_peak_seq), interval="confidence")

points(x=exp(l_AB_peak_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(l_AB_peak_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_peak_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )





##############################
##                          ## 
##  PANEL D                 ##
##  ELISA reduction vs age  ##
##                          ##
##############################


line_seq_x <- c(0.1, 1, 10, 100, 1000, 3000)

line_seq_y <- c(0.5, 1, 2, 5, 10, 20, 50)



plot( x=1e10, y=1e10,  
xlim=c(0.5,50), ylim=c(0.08,3000), log="xy", 
xlab="age (years)", 
ylab=expression(paste( "IgG ELISA after 1 year (", mu, "g/mL)", sep="" )),
main="(D) IgG ELISA vs age",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$age, y=AB_df_all$elisa_12m, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1, at=c(0.5, 1, 2, 5, 10, 20, 50), labels=c(0.5, 1, 2, 5, 10, 20, 50), cex.axis=axis.size)

axis(2,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$elisa_12m)
age_x  <- log(AB_df_all$age)

AB_red[which(AB_red < 0)] <- NA

age_seq <- seq(from=log(0.5), to=log(50), by=0.01)

lin_mod <- lm( AB_red ~ age_x )

lin_mod_fit <- predict( lin_mod, data.frame(age_x=age_seq), interval="confidence")

points(x=exp(age_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )



#####################################
##                                 ## 
##  PANEL E                        ##
##  ELISA reduction vs ELISA_base  ##
##                                 ##
#####################################


line_seq_x <- c(0.1, 1, 10, 100, 1000, 3000)
line_seq_y <- c(0.1, 1, 10, 100, 1000, 3000)

plot( x=1e10, y=1e10,  
xlim=c(0.08,3000), ylim=c(0.08,3000), log="xy", 
xlab=expression(paste( "pre-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
ylab=expression(paste( "IgG ELISA after 1 year (", mu, "g/mL)", sep="" )),
main="(E) IgG ELISA vs pre-MenAfriVac IgG ELISA",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$elisa_base, y=AB_df_all$elisa_12m, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)

axis(2,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$elisa_12m)
l_AB_base <- log(AB_df_all$elisa_base)

AB_red[which(AB_red < 0)] <- NA

l_AB_base_seq <- seq(from=log(0.005), to=log(2^20), by=0.01)

lin_mod <- lm( AB_red ~ l_AB_base )

lin_mod_fit <- predict( lin_mod, data.frame(l_AB_base=l_AB_base_seq), interval="confidence")

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )






###################################
##                               ## 
##  PANEL F                      ##
##  ELISA reduction vs SBA_peak  ##
##                               ##
###################################


line_seq_x <- c(0.1, 1, 10, 100, 1000, 3000)
line_seq_y <- c(0.1, 1, 10, 100, 1000, 3000)


plot( x=1e10, y=1e10,  
xlim=c(0.08,3000), ylim=c(0.08,3000), log="xy", 
xlab=expression(paste( "post-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
ylab=expression(paste( "IgG ELISA after 1 year (", mu, "g/mL)", sep="" )),
main="(F) IgG ELISA vs post-MenAfriVac IgG ELISA",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$elisa_peak, y=AB_df_all$elisa_12m, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)

axis(2,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$elisa_12m)
l_AB_peak <- log(AB_df_all$elisa_peak)

AB_red[which(AB_red < 0)] <- NA

l_AB_peak_seq <- seq(from=log(0.005), to=log(2^20), by=0.01)

lin_mod <- lm( AB_red ~ l_AB_peak )

lin_mod_fit <- predict( lin_mod, data.frame(l_AB_peak=l_AB_peak_seq), interval="confidence")

points(x=exp(l_AB_peak_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(l_AB_peak_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_peak_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )




###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("Gambia (PsA-TT-002)", "Gambia (PsA-TT-003)", "Mali (PsA-TT-002)", "Mali (PsA-TT-003)", "Senegal (PsA-TT-003)"), 
       col = c("gold", "gold", "forestgreen", "forestgreen", "royalblue"), 
       pch=c(19, 17, 19, 17, 17),
       ncol=5, cex=1.25, bty="n" )


dev.off()

	














########################################################
########################################################
##                                                    ##
##   ####  ##  ## #####    ##### ####  ####    #####  ##
##  ##     ##  ## ##  ##   ##     ##  ##          ##  ##  
##   ####  ##  ## #####    ####   ##  ## ###     ##   ##
##      ## ##  ## ##       ##     ##  ##  ##    ##    ## 
##   ####   ####  ##       ##    ####  ####    ##     ##
##                                                    ##
########################################################
########################################################


main.size = 0.9
axis.size = 0.75
lab.size  = 1.0


tiff( file="2_SupFig4_Ab.tif", width=16, height=18, units="cm", res=500)

lay.mat = rbind( c( 1,  2  ),
                 c( 3,  4  ),
                 c( 5,  5  ) )

layout(lay.mat, heights=c(10,10,1))
layout.show(5)



par(mar = c(2.5,2.5,1,1.75))
par(mgp=c(1.5,0.5,0))



############################
##                        ## 
##  PANEL A               ##
##  SBA reduction vs age  ##
##                        ##
############################


line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

line_seq_y <- c(0.5, 1, 2, 5, 10, 20, 50)



plot( x=1e10, y=1e10,  
xlim=c(0.5,50), ylim=c(2^(-0.25),2^19), log="xy", 
xlab="age (years)", 
ylab="post-MenAfriVac SBA titer",
main="(A) post-MenAfriVac SBA titer vs age",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$age, y=AB_df_all$SBA_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1, at=c(0.5, 1, 2, 5, 10, 20, 50), labels=c(0.5, 1, 2, 5, 10, 20, 50), cex.axis=axis.size)

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$SBA_peak)
age_x  <- log(AB_df_all$age)

AB_red[which(AB_red < 0)] <- NA

age_seq <- seq(from=log(0.5), to=log(50), by=0.01)

lin_mod <- lm( AB_red ~ age_x )

lin_mod_fit <- predict( lin_mod, data.frame(age_x=age_seq), interval="confidence")

points(x=exp(age_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )



#################################
##                             ## 
##  PANEL B                    ##
##  SBA reduction vs SBA_base  ##
##                             ##
#################################


line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)
line_seq_y <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)

plot( x=1e10, y=1e10,  
xlim=c(2^(-0.25),2^19), ylim=c(2^(-0.25),2^19), log="xy", 
xlab="pre-MenAfriVac SBA titer", 
ylab="post-MenAfriVac SBA titer",
main="(B) post-MenAfriVac SBA titer vs pre-MenAfriVac SBA titer",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$SBA_base, y=AB_df_all$SBA_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "524288" ), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$SBA_peak)
l_AB_base <- log(AB_df_all$SBA_base)

AB_red[which(AB_red < 0)] <- NA

l_AB_base_seq <- seq(from=log(0.5), to=log(2^20), by=0.01)

lin_mod <- lm( AB_red ~ l_AB_base )

lin_mod_fit <- predict( lin_mod, data.frame(l_AB_base=l_AB_base_seq), interval="confidence")

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )








##############################
##                          ## 
##  PANEL C                 ##
##  ELISA reduction vs age  ##
##                          ##
##############################


line_seq_x <- c(0.1, 1, 10, 100, 1000, 3000)

line_seq_y <- c(0.5, 1, 2, 5, 10, 20, 50)



plot( x=1e10, y=1e10,  
xlim=c(0.5,50), ylim=c(0.08,3000), log="xy", 
xlab="age (years)", 
ylab=expression(paste( "post-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
main="(C) post-MenAfriVac IgG ELISA vs age",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$age, y=AB_df_all$elisa_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1, at=c(0.5, 1, 2, 5, 10, 20, 50), labels=c(0.5, 1, 2, 5, 10, 20, 50), cex.axis=axis.size)

axis(2,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$elisa_peak)
age_x  <- log(AB_df_all$age)

AB_red[which(AB_red < 0)] <- NA

age_seq <- seq(from=log(0.5), to=log(50), by=0.01)

lin_mod <- lm( AB_red ~ age_x )

lin_mod_fit <- predict( lin_mod, data.frame(age_x=age_seq), interval="confidence")

points(x=exp(age_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(age_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )



#####################################
##                                 ## 
##  PANEL D                        ##
##  ELISA reduction vs ELISA_base  ##
##                                 ##
#####################################


line_seq_x <- c(0.1, 1, 10, 100, 1000, 3000)
line_seq_y <- c(0.1, 1, 10, 100, 1000, 3000)

plot( x=1e10, y=1e10,  
xlim=c(0.08,3000), ylim=c(0.08,3000), log="xy", 
xlab=expression(paste( "pre-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
ylab=expression(paste( "post-MenAfriVac IgG ELISA (", mu, "g/mL)", sep="" )),
main="(D) post-MenAfriVac IgG ELISA vs pre-MenAfriVac IgG ELISA",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


points(x=AB_df_all$elisa_base, y=AB_df_all$elisa_peak, pch=site_shape, col=site_col)


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


axis(1,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)

axis(2,  at = c( 0.1, 1, 10, 100, 1000, 3000),  labels = c( "0.1", "1", "10", "100", "1000", ""), 
cex.axis=axis.size)



AB_red <- log(AB_df_all$elisa_peak)
l_AB_base <- log(AB_df_all$elisa_base)

AB_red[which(AB_red < 0)] <- NA

l_AB_base_seq <- seq(from=log(0.005), to=log(2^20), by=0.01)

lin_mod <- lm( AB_red ~ l_AB_base )

lin_mod_fit <- predict( lin_mod, data.frame(l_AB_base=l_AB_base_seq), interval="confidence")

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,1]),  type='l', lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,2]),  type='l', 
lty="longdash", col="grey", lwd=2 )

points(x=exp(l_AB_base_seq), y=exp(lin_mod_fit[,3]), type='l', 
lty="longdash", col="grey", lwd=2 )






###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("Gambia (PsA-TT-002)", "Gambia (PsA-TT-003)", "Mali (PsA-TT-002)", "Mali (PsA-TT-003)", "Senegal (PsA-TT-003)"), 
       col = c("gold", "gold", "forestgreen", "forestgreen", "royalblue"), 
       pch=c(19, 17, 19, 17, 17),
       ncol=5, cex=0.9, bty="n" )


dev.off()

	























