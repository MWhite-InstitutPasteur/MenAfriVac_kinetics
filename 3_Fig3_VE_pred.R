################################
################################
################################
###                          ###
###    ####   ####   ####    ###
###   ##  ## ##  ## ##  ##   ###
###   ##  ## ##  ##    ##    ###
###   ##  ## ##  ##   ##     ###
###    ####   ####   #####   ### 
###                          ###
################################
################################
################################


VE_beta  = 128

VE_alpha_1 = 2000

VE_alpha_1 = 4


N_sam = 100


x_1 = 0.1
y_1 = 2

x_2 = 2
y_2 = 128

ELISA_axis_transform = function( AB_elisa )
{
	exp( log(y_1) + ( (log(y_2) - log(y_1) )/( log(x_2) - log(x_1) ))*( log(AB_elisa) - log(x_1) ) )
}




tt_plot = seq(from=0, to=20*365, by=5)
log2 = log(2)
N_tt = length(tt_plot)




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
	MCMC_ind_n = MCMC_MenA_002[seq(from=1, to=nrow(MCMC_MenA_002), length=N_sam),((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]


	AB_mod = matrix(NA, nrow=N_sam, ncol=length(tt_plot))
	
	for(i in 1:N_sam)
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


		AB_tt = Ab_0*exp(-r_cl*tt_plot) 


		if( (AB_Cpp_002[n,7] == 1) && (AB_Cpp_002[n,8] == 1) )
		{
			AB_tt = AB_tt + beta_prim_PsA*(       (rho_prim_PsA/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                          ((1 - rho_prim_PsA)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )

			tt_temp = tt_plot[which(tt_plot>t_vac_boost)] - t_vac_boost

			AB_tt[which(tt_plot>t_vac_boost)] = AB_tt[which(tt_plot>t_vac_boost)] +
     		                                                   beta_boost_PsA*(        (rho_boost_PsA/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
					                                                   ((1 - rho_boost_PsA)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp)) )
		}


		if( (AB_Cpp_002[n,7] == 1) && (AB_Cpp_002[n,8] == 2) )
		{
			AB_tt = AB_tt + beta_prim_PsA*(       (rho_prim_PsA/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                          ((1 - rho_prim_PsA)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )

			tt_temp = tt_plot[which(tt_plot>t_vac_boost)] - t_vac_boost

			AB_tt[which(tt_plot>t_vac_boost)] = AB_tt[which(tt_plot>t_vac_boost)] +
     		                                                   beta_Hib*(        (rho_Hib/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
					                                                   ((1 - rho_Hib)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp)) )
		}


		if( (AB_Cpp_002[n,7] == 2) && (AB_Cpp_002[n,8] == 1) )
		{
			AB_tt = AB_tt + beta_Hib*(       (rho_Hib/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                     ((1 - rho_Hib)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )

			tt_temp = tt_plot[which(tt_plot>t_vac_boost)] - t_vac_boost

			AB_tt[which(tt_plot>t_vac_boost)] = AB_tt[which(tt_plot>t_vac_boost)] +
     		                                                   beta_prim_PsA*(        (rho_prim_PsA/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
					                                                   ((1 - rho_prim_PsA)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp)) )
		}

		AB_mod[i,] = AB_tt
	}

	###################
	## Return output

	OUTPUT <- list()

	OUTPUT[[1]] = AB_mod
	OUTPUT[[2]] = AB
	OUTPUT[[3]] = tt

	OUTPUT
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
	MCMC_ind_n = MCMC_logelisa_002[seq(from=1, to=nrow(MCMC_logelisa_002), length=N_sam),((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]

	AB_mod = matrix(NA, nrow=N_sam, ncol=length(tt_plot))
	
	for(i in 1:N_sam)
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

		AB_tt = Ab_0*exp(-r_cl*tt_plot) 

		if( (AB_Cpp_002[n,7] == 1) && (AB_Cpp_002[n,8] == 2) )
		{
			AB_tt = AB_tt + beta_prim*(       (rho_prim/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                      ((1 - rho_prim)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )
		}


		if( (AB_Cpp_002[n,7] == 2) && (AB_Cpp_002[n,8] == 1) )
		{
			tt_temp = tt_plot[which(tt_plot>t_vac_boost)] - t_vac_boost

			AB_tt[which(tt_plot>t_vac_boost)] = AB_tt[which(tt_plot>t_vac_boost)] +
	     	                                                   beta_prim*(         (rho_prim/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
					                                               ((1 - rho_prim)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp)) )
		}


		if( (AB_Cpp_002[n,7] == 1) && (AB_Cpp_002[n,8] == 1) )
		{
			AB_tt = AB_tt + beta_prim*(       (rho_prim/(r_a - r_cs))*( exp(-r_cs*tt_plot) - exp(-r_a*tt_plot) ) +
				                      ((1 - rho_prim)/(r_a - r_cl))*( exp(-r_cl*tt_plot) - exp(-r_a*tt_plot)) )

			tt_temp = tt_plot[which(tt_plot>t_vac_boost)] - t_vac_boost

			AB_tt[which(tt_plot>t_vac_boost)] = AB_tt[which(tt_plot>t_vac_boost)] +
	     	                                                   beta_boost*(        (rho_boost/(r_a - r_cs))*( exp(-r_cs*tt_temp) - exp(-r_a*tt_temp) ) +
					                                               ((1 - rho_boost)/(r_a - r_cl))*( exp(-r_cl*tt_temp) - exp(-r_a*tt_temp)) )
		}

		AB_mod[i,] = AB_tt
	}
	
	###################
	## Return output

	OUTPUT <- list()

	OUTPUT[[1]] = ELISA_axis_transform( AB_mod )
	OUTPUT[[2]] = ELISA_axis_transform( AB )
	OUTPUT[[3]] = tt

	OUTPUT
}







PsA_PsA_002_index = intersect( which(AB_data_002$prim_vac=="PsA-TT"), which(AB_data_002$boost_vac=="PsA-TT") )  
N_PsA_PsA_002 = length(PsA_PsA_002_index)


PsA_Hib_002_index = intersect( which(AB_data_002$prim_vac=="PsA-TT"), which(AB_data_002$boost_vac=="Hib-TT") )  
N_PsA_Hib_002 = length(PsA_Hib_002_index)


Hib_PsA_002_index = intersect( which(AB_data_002$prim_vac=="Hib-TT"), which(AB_data_002$boost_vac=="PsA-TT") )  
N_Hib_PsA_002 = length(Hib_PsA_002_index)


PsA_PsA_MenA_002_array_1 = array(NA, dim=c(N_PsA_PsA_002, N_sam, N_tt) )
PsA_Hib_MenA_002_array_1 = array(NA, dim=c(N_PsA_Hib_002, N_sam, N_tt) )
Hib_PsA_MenA_002_array_1 = array(NA, dim=c(N_Hib_PsA_002, N_sam, N_tt) )


PsA_PsA_logelisa_002_array_1 = array(NA, dim=c(N_PsA_PsA_002, N_sam, N_tt) )
PsA_Hib_logelisa_002_array_1 = array(NA, dim=c(N_PsA_Hib_002, N_sam, N_tt) )
Hib_PsA_logelisa_002_array_1 = array(NA, dim=c(N_Hib_PsA_002, N_sam, N_tt) )




PsA_PsA_MenA_002_array_2 = array(NA, dim=c(N_PsA_PsA_002, N_sam, N_tt) )
PsA_Hib_MenA_002_array_2 = array(NA, dim=c(N_PsA_Hib_002, N_sam, N_tt) )
Hib_PsA_MenA_002_array_2 = array(NA, dim=c(N_Hib_PsA_002, N_sam, N_tt) )


PsA_PsA_logelisa_002_array_2 = array(NA, dim=c(N_PsA_PsA_002, N_sam, N_tt) )
PsA_Hib_logelisa_002_array_2 = array(NA, dim=c(N_PsA_Hib_002, N_sam, N_tt) )
Hib_PsA_logelisa_002_array_2 = array(NA, dim=c(N_Hib_PsA_002, N_sam, N_tt) )



for(i in 1:N_PsA_PsA_002)
{
	PsA_PsA_MenA_002_array_1[i,,]     = MenA_002_pred( PsA_PsA_002_index[i] )[[1]]
	PsA_PsA_logelisa_002_array_1[i,,] = logelisa_002_pred( PsA_PsA_002_index[i] )[[1]]


	PsA_PsA_MenA_002_array_2[i,,]     = PsA_PsA_MenA_002_array_1[i,,]
	PsA_PsA_logelisa_002_array_2[i,,] = PsA_PsA_logelisa_002_array_1[i,,]
}

PsA_PsA_MenA_002_array_1     = 1 - 1/( 1 + (PsA_PsA_MenA_002_array_1/VE_beta)^VE_alpha_1 )
PsA_PsA_logelisa_002_array_1 = 1 - 1/( 1 + (PsA_PsA_logelisa_002_array_1/VE_beta)^VE_alpha_1 )

PsA_PsA_MenA_002_VE_1     = colSums(PsA_PsA_MenA_002_array_1)/N_PsA_PsA_002
PsA_PsA_logelisa_002_VE_1 = colSums(PsA_PsA_logelisa_002_array_1)/N_PsA_PsA_002

rm(PsA_PsA_MenA_002_array_1)
rm(PsA_PsA_logelisa_002_array_1)



PsA_PsA_MenA_002_array_2     = 1 - 1/( 1 + (PsA_PsA_MenA_002_array_2/VE_beta)^VE_alpha_2 )
PsA_PsA_logelisa_002_array_2 = 1 - 1/( 1 + (PsA_PsA_logelisa_002_array_2/VE_beta)^VE_alpha_2 )

PsA_PsA_MenA_002_VE_2     = colSums(PsA_PsA_MenA_002_array_2)/N_PsA_PsA_002
PsA_PsA_logelisa_002_VE_2 = colSums(PsA_PsA_logelisa_002_array_2)/N_PsA_PsA_002

rm(PsA_PsA_MenA_002_array_2)
rm(PsA_PsA_logelisa_002_array_2)



for(i in 1:N_PsA_Hib_002)
{
	PsA_Hib_MenA_002_array_1[i,,]     = MenA_002_pred( PsA_Hib_002_index[i] )[[1]]
	PsA_Hib_logelisa_002_array_1[i,,] = logelisa_002_pred( PsA_Hib_002_index[i] )[[1]]


	PsA_Hib_MenA_002_array_2[i,,]     = PsA_Hib_MenA_002_array_1[i,,]
	PsA_Hib_logelisa_002_array_2[i,,] = PsA_Hib_logelisa_002_array_1[i,,]
}

PsA_Hib_MenA_002_array_1     = 1 - 1/( 1 + (PsA_Hib_MenA_002_array_1/VE_beta)^VE_alpha_1 )
PsA_Hib_logelisa_002_array_1 = 1 - 1/( 1 + (PsA_Hib_logelisa_002_array_1/VE_beta)^VE_alpha_1 )

PsA_Hib_MenA_002_VE_1     = colSums(PsA_Hib_MenA_002_array_1)/N_PsA_Hib_002
PsA_Hib_logelisa_002_VE_1 = colSums(PsA_Hib_logelisa_002_array_1)/N_PsA_Hib_002

rm(PsA_Hib_MenA_002_array_1)
rm(PsA_Hib_logelisa_002_array_1)


PsA_Hib_MenA_002_array_2     = 1 - 1/( 1 + (PsA_Hib_MenA_002_array_2/VE_beta)^VE_alpha_2 )
PsA_Hib_logelisa_002_array_2 = 1 - 1/( 1 + (PsA_Hib_logelisa_002_array_2/VE_beta)^VE_alpha_2 )

PsA_Hib_MenA_002_VE_2     = colSums(PsA_Hib_MenA_002_array_2)/N_PsA_Hib_002
PsA_Hib_logelisa_002_VE_2 = colSums(PsA_Hib_logelisa_002_array_2)/N_PsA_Hib_002

rm(PsA_Hib_MenA_002_array_2)
rm(PsA_Hib_logelisa_002_array_2)




for(i in 1:N_Hib_PsA_002)
{
	Hib_PsA_MenA_002_array_1[i,,]     = MenA_002_pred( Hib_PsA_002_index[i] )[[1]]
	Hib_PsA_logelisa_002_array_1[i,,] = logelisa_002_pred( Hib_PsA_002_index[i] )[[1]]


	Hib_PsA_MenA_002_array_2[i,,]     = Hib_PsA_MenA_002_array_1[i,,]
	Hib_PsA_logelisa_002_array_2[i,,] = Hib_PsA_logelisa_002_array_1[i,,]
}

Hib_PsA_MenA_002_array_1     = 1 - 1/( 1 + (Hib_PsA_MenA_002_array_1/VE_beta)^VE_alpha_1 )
Hib_PsA_logelisa_002_array_1 = 1 - 1/( 1 + (Hib_PsA_logelisa_002_array_1/VE_beta)^VE_alpha_1 )

Hib_PsA_MenA_002_VE_1     = colSums(Hib_PsA_MenA_002_array_1)/N_Hib_PsA_002
Hib_PsA_logelisa_002_VE_1 = colSums(Hib_PsA_logelisa_002_array_1)/N_Hib_PsA_002

rm(Hib_PsA_MenA_002_array_1)
rm(Hib_PsA_logelisa_002_array_1)


Hib_PsA_MenA_002_array_2     = 1 - 1/( 1 + (Hib_PsA_MenA_002_array_2/VE_beta)^VE_alpha_2 )
Hib_PsA_logelisa_002_array_2 = 1 - 1/( 1 + (Hib_PsA_logelisa_002_array_2/VE_beta)^VE_alpha_2 )

Hib_PsA_MenA_002_VE_2     = colSums(Hib_PsA_MenA_002_array_2)/N_Hib_PsA_002
Hib_PsA_logelisa_002_VE_2 = colSums(Hib_PsA_logelisa_002_array_2)/N_Hib_PsA_002

rm(Hib_PsA_MenA_002_array_2)
rm(Hib_PsA_logelisa_002_array_2)











PsA_PsA_MenA_002_VE_1_quant = matrix(NA, nrow=5, ncol=N_tt)
PsA_Hib_MenA_002_VE_1_quant = matrix(NA, nrow=5, ncol=N_tt)
Hib_PsA_MenA_002_VE_1_quant = matrix(NA, nrow=5, ncol=N_tt)



PsA_PsA_logelisa_002_VE_1_quant = matrix(NA, nrow=5, ncol=N_tt)
PsA_Hib_logelisa_002_VE_1_quant = matrix(NA, nrow=5, ncol=N_tt)
Hib_PsA_logelisa_002_VE_1_quant = matrix(NA, nrow=5, ncol=N_tt)



for(j in 1:N_tt)
{
	PsA_PsA_MenA_002_VE_1_quant[,j] <- quantile( PsA_PsA_MenA_002_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	PsA_Hib_MenA_002_VE_1_quant[,j] <- quantile( PsA_Hib_MenA_002_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	Hib_PsA_MenA_002_VE_1_quant[,j] <- quantile( Hib_PsA_MenA_002_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )


	PsA_PsA_logelisa_002_VE_1_quant[,j] <- quantile( PsA_PsA_logelisa_002_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	PsA_Hib_logelisa_002_VE_1_quant[,j] <- quantile( PsA_Hib_logelisa_002_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	Hib_PsA_logelisa_002_VE_1_quant[,j] <- quantile( Hib_PsA_logelisa_002_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
}






PsA_PsA_MenA_002_VE_2_quant = matrix(NA, nrow=5, ncol=N_tt)
PsA_Hib_MenA_002_VE_2_quant = matrix(NA, nrow=5, ncol=N_tt)
Hib_PsA_MenA_002_VE_2_quant = matrix(NA, nrow=5, ncol=N_tt)


PsA_PsA_logelisa_002_VE_2_quant = matrix(NA, nrow=5, ncol=N_tt)
PsA_Hib_logelisa_002_VE_2_quant = matrix(NA, nrow=5, ncol=N_tt)
Hib_PsA_logelisa_002_VE_2_quant = matrix(NA, nrow=5, ncol=N_tt)



for(j in 1:N_tt)
{
	PsA_PsA_MenA_002_VE_2_quant[,j] <- quantile( PsA_PsA_MenA_002_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	PsA_Hib_MenA_002_VE_2_quant[,j] <- quantile( PsA_Hib_MenA_002_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	Hib_PsA_MenA_002_VE_2_quant[,j] <- quantile( Hib_PsA_MenA_002_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )


	PsA_PsA_logelisa_002_VE_2_quant[,j] <- quantile( PsA_PsA_logelisa_002_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	PsA_Hib_logelisa_002_VE_2_quant[,j] <- quantile( PsA_Hib_logelisa_002_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	Hib_PsA_logelisa_002_VE_2_quant[,j] <- quantile( Hib_PsA_logelisa_002_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
}









################################
################################
################################
###                          ###
###    ####   ####   ####    ###
###   ##  ## ##  ## ##  ##   ###
###   ##  ## ##  ##    ##    ###
###   ##  ## ##  ## ##  ##   ###
###    ####   ####   ####    ### 
###                          ###
################################
################################
################################

N_sam = 100

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
	MCMC_ind_n = MCMC_MenA_003[seq(from=1, to=nrow(MCMC_MenA_003), length=N_sam),((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]

	AB_mod = matrix(NA, nrow=N_sam, ncol=length(tt_plot))
	
	for(i in 1:N_sam)
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

		AB_mod[i,] = AB_tt
	}

	###################
	## Return output

	OUTPUT <- list()

	OUTPUT[[1]] = AB_mod
	OUTPUT[[2]] = AB
	OUTPUT[[3]] = tt

	OUTPUT
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
	MCMC_ind_n = MCMC_logelisa_003[seq(from=1, to=nrow(MCMC_logelisa_003), length=N_sam),((N_loc_par+1)*(n-1)+1):((N_loc_par+1)*n)]

	AB_mod = matrix(NA, nrow=N_sam, ncol=length(tt_plot))
	
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

		AB_mod[i,] = AB_tt
	}


	###################
	## Return output

	OUTPUT <- list()

	OUTPUT[[1]] = ELISA_axis_transform( AB_mod )
	OUTPUT[[2]] = ELISA_axis_transform( AB )
	OUTPUT[[3]] = tt

	OUTPUT
}







PsA_003_index = 1:nrow(AB_data_003) 
N_PsA_003 = length(PsA_003_index)





PsA_MenA_003_array_1  = array(NA, dim=c(N_PsA_003, N_sam, N_tt) )

PsA_logelisa_003_array_1  = array(NA, dim=c(N_PsA_003, N_sam, N_tt) )


##PsA_MenA_003_array_2  = array(NA, dim=c(N_PsA_003, N_sam, N_tt) )
##
##PsA_logelisa_003_array_2  = array(NA, dim=c(N_PsA_003, N_sam, N_tt) )



for(i in 1:N_PsA_003)
{
	PsA_MenA_003_array_1[i,,]     = MenA_003_pred( PsA_003_index[i] )[[1]]
	PsA_logelisa_003_array_1[i,,] = logelisa_003_pred( PsA_003_index[i] )[[1]]


	#PsA_MenA_003_array_2[i,,]     = PsA_MenA_003_array_1[i,,]
	#PsA_logelisa_003_array_2[i,,] = PsA_logelisa_003_array_1[i,,]
}

PsA_MenA_003_array_1     = 1 - 1/( 1 + (PsA_MenA_003_array_1/VE_beta)^VE_alpha_1 )
PsA_logelisa_003_array_1 = 1 - 1/( 1 + (PsA_logelisa_003_array_1/VE_beta)^VE_alpha_1 )

PsA_MenA_003_VE_1     = colSums(PsA_MenA_003_array_1)/N_PsA_003
PsA_logelisa_003_VE_1 = colSums(PsA_logelisa_003_array_1)/N_PsA_003

rm(PsA_MenA_003_array_1)
rm(PsA_logelisa_003_array_1)



##PsA_MenA_003_array_2     = 1 - 1/( 1 + (PsA_MenA_003_array_2/VE_beta)^VE_alpha_2 )
##PsA_logelisa_003_array_2 = 1 - 1/( 1 + (PsA_logelisa_003_array_2/VE_beta)^VE_alpha_2 )
##
##PsA_MenA_003_VE_2     = colSums(PsA_MenA_003_array_2)/N_PsA_003
##PsA_logelisa_003_VE_2 = colSums(PsA_logelisa_003_array_2)/N_PsA_003
##
##rm(PsA_MenA_003_array_2)
##rm(PsA_logelisa_003_array_2)



PsA_MenA_003_VE_1_quant  = matrix(NA, nrow=5, ncol=N_tt)


PsA_logelisa_003_VE_1_quant  = matrix(NA, nrow=5, ncol=N_tt)


for(j in 1:N_tt)
{
	PsA_MenA_003_VE_1_quant[,j]  <- quantile( PsA_MenA_003_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )

	PsA_logelisa_003_VE_1_quant[,j]  <- quantile( PsA_logelisa_003_VE_1[,j], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )
}




##PsA_MenA_003_VE_2_quant  = matrix(NA, nrow=5, ncol=N_tt)
##
##PsA_logelisa_003_VE_2_quant  = matrix(NA, nrow=5, ncol=N_tt)
##
##
##for(j in 1:N_tt)
##{
##	PsA_MenA_003_VE_2_quant[,j]  <- quantile( PsA_MenA_003_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )
##
##	PsA_logelisa_003_VE_2_quant[,j]  <- quantile( PsA_logelisa_003_VE_2[,j], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )
##}















#################################
#################################
##                             ##
##  ##### ####  ####    ####   ##
##  ##     ##  ##      ##  ##  ##
##  ####   ##  ## ###     ##   ## 
##  ##     ##  ##  ##  ##  ##  ## 
##  ##    ####  ####    ####   ##
##                             ##
################################# 
#################################

line_seq_x <- c(0, 0.25, 0.5, 0.75, 1)

line_seq_y <- 365*c(0, 5, 10, 15, 20)




tiff( file="3_Fig3_VE_pred.tif", width=24, height=16, units="cm", res=500)

lay.mat = rbind( c( 1,   2,  3 ),
                 c( 4,   5,  6 ) )

layout(lay.mat, heights=c(10,10))
layout.show(6)





main.size = 1.2
axis.size = 0.8
lab.size  = 1.2


par(mar = c(3, 3 ,1.5, 0.5))
par(mgp=c(1.5,0.6,0))



############################
##                        ## 
##  PANEL 1               ##
##  PsA_PsA_002           ##
##                        ##
############################


plot( x=-100, y=-100,  
xlim=c(-50,20*365), ylim=c(0,1.02), 
xlab="time (years)", ylab="vaccine efficacy",
main="(A) PsA-TT-002: MenAfriVac/MenAfriVac",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}




polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_PsA_MenA_002_VE_1_quant[1,], rev(PsA_PsA_MenA_002_VE_1_quant[5,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_PsA_MenA_002_VE_1_quant[2,], rev(PsA_PsA_MenA_002_VE_1_quant[4,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)


points(x=tt_plot, y=PsA_PsA_MenA_002_VE_1_quant[3,], type='l', lwd=2, col="firebrick1")



axis(1, at=365*c(0, 5, 10, 15, 20), labels=c(0, 5, 10, 15, 20), cex.axis=axis.size)

axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)



############################
##                        ## 
##  PANEL 2               ##
##  Gambia_PsA_Hib_002    ##
##                        ##
############################


plot( x=-100, y=-100,  
xlim=c(-50,20*365), ylim=c(0,1.02), 
xlab="time (years)", ylab="vaccine efficacy",
main="(B) PsA-TT-002: MenAfriVac/Hib-TT",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}



polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_Hib_MenA_002_VE_1_quant[1,], rev(PsA_Hib_MenA_002_VE_1_quant[5,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_Hib_MenA_002_VE_1_quant[2,], rev(PsA_Hib_MenA_002_VE_1_quant[4,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)


points(x=tt_plot, y=PsA_Hib_MenA_002_VE_1_quant[3,], type='l', lwd=2, col="firebrick1")




axis(1, at=365*c(0, 5, 10, 15, 20), labels=c(0, 5, 10, 15, 20), cex.axis=axis.size)

axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)



############################
##                        ## 
##  PANEL 3               ##
##  Hib_PsA_002           ##
##                        ##
############################


plot( x=-100, y=-100,  
xlim=c(-50,20*365), ylim=c(0,1.02), 
xlab="time (years)", ylab="vaccine efficacy",
main="(C) PsA-TT-002: Hib-TT/MenAfriVac",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}






polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( Hib_PsA_MenA_002_VE_1_quant[1,], rev(Hib_PsA_MenA_002_VE_1_quant[5,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( Hib_PsA_MenA_002_VE_1_quant[2,], rev(Hib_PsA_MenA_002_VE_1_quant[4,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)



points(x=tt_plot, y=Hib_PsA_MenA_002_VE_1_quant[3,], type='l', lwd=2, col="firebrick1")



axis(1, at=365*c(0, 5, 10, 15, 20), labels=c(0, 5, 10, 15, 20), cex.axis=axis.size)

axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)




############################
##                        ## 
##  PANEL 4               ##
##  Gambia_003            ##
##                        ##
############################


plot( x=-100, y=-100,  
xlim=c(-50,20*365), ylim=c(0,1.02), 
xlab="time (years)", ylab="vaccine efficacy",
main="(D) PsA-TT-003: MenAfriVac",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_MenA_003_VE_1_quant[1,], rev(PsA_MenA_003_VE_1_quant[5,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_MenA_003_VE_1_quant[2,], rev(PsA_MenA_003_VE_1_quant[4,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

points(x=tt_plot, y=PsA_MenA_003_VE_1_quant[3,], type='l', lwd=2, col="firebrick1")



axis(1, at=365*c(0, 5, 10, 15, 20), labels=c(0, 5, 10, 15, 20), cex.axis=axis.size)

axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)





############################
##                        ## 
##  PANEL 5               ##
##  Mali_003              ##
##                        ##
############################

line_seq_x <- c(0, 0.25, 0.5, 0.75, 1)


line_seq_y <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)



plot( x=1e10, y=1e10,  
xlim=c(2^(0),2^20), ylim=c(0,1.02), log="x", 
xlab="SBA titer", ylab="vaccine efficacy",
main="(E) Assumed dose-response relationship",
xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size )


for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}


MenA_seq = exp(seq(from=log(1), to=log(1e10), length=10000))

VE_DR = 1 - 1/( 1 + (MenA_seq/VE_beta)^VE_alpha_1 )

points(x=MenA_seq, y=VE_DR, type='l', lwd=2, col="black")



axis(1,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=axis.size)

axis(2, at=c(0.0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)






############################
##                        ## 
##  PANEL 6               ##
##                        ##
############################

plot.new()




dev.off()

	






