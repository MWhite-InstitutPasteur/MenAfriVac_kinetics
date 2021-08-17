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

		AB_tt[which(tt_plot<0)] = Ab_0

		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975) )
	}


	###################
	## Return output

	OUTPUT <- list()

	OUTPUT[[1]] = AB_quant
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

		AB_tt[which(tt_plot<0)] = Ab_0

		AB_mod[i,] = AB_tt
	}

	AB_quant <- matrix(NA, nrow=3, ncol=length(tt_plot))
	for(j in 1:length(tt_plot))
	{
		AB_quant[,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975) )
	}

	
	###################
	## Return output

	OUTPUT <- list()

	OUTPUT[[1]] = ELISA_axis_transform( AB_quant )
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






tt_plot = seq(from=0, to=5.2*365, by=1)
log2 = log(2)
N_tt = length(tt_plot)


PsA_PsA_MenA_002_mat = matrix(NA, nrow=N_PsA_PsA_002, ncol=N_tt)
PsA_Hib_MenA_002_mat = matrix(NA, nrow=N_PsA_Hib_002, ncol=N_tt)
Hib_PsA_MenA_002_mat = matrix(NA, nrow=N_Hib_PsA_002, ncol=N_tt)



PsA_PsA_logelisa_002_mat = matrix(NA, nrow=N_PsA_PsA_002, ncol=N_tt)
PsA_Hib_logelisa_002_mat = matrix(NA, nrow=N_PsA_Hib_002, ncol=N_tt)
Hib_PsA_logelisa_002_mat = matrix(NA, nrow=N_Hib_PsA_002, ncol=N_tt)



for(i in 1:N_PsA_PsA_002)
{
	PsA_PsA_MenA_002_mat[i,]     = MenA_002_pred( PsA_PsA_002_index[i] )[[1]][2,]
	PsA_PsA_logelisa_002_mat[i,] = logelisa_002_pred( PsA_PsA_002_index[i] )[[1]][2,]
}

for(i in 1:N_PsA_Hib_002)
{
	PsA_Hib_MenA_002_mat[i,]     = MenA_002_pred( PsA_Hib_002_index[i] )[[1]][2,]
	PsA_Hib_logelisa_002_mat[i,] = logelisa_002_pred( PsA_Hib_002_index[i] )[[1]][2,]
}

for(i in 1:N_Hib_PsA_002)
{
	Hib_PsA_MenA_002_mat[i,]     = MenA_002_pred( Hib_PsA_002_index[i] )[[1]][2,]
	Hib_PsA_logelisa_002_mat[i,] = logelisa_002_pred( Hib_PsA_002_index[i] )[[1]][2,]
}


PsA_PsA_MenA_002_quant = matrix(NA, nrow=5, ncol=N_tt)
PsA_Hib_MenA_002_quant = matrix(NA, nrow=5, ncol=N_tt)
Hib_PsA_MenA_002_quant = matrix(NA, nrow=5, ncol=N_tt)

PsA_PsA_logelisa_002_quant = matrix(NA, nrow=5, ncol=N_tt)
PsA_Hib_logelisa_002_quant = matrix(NA, nrow=5, ncol=N_tt)
Hib_PsA_logelisa_002_quant = matrix(NA, nrow=5, ncol=N_tt)

for(j in 1:N_tt)
{
	PsA_PsA_MenA_002_quant[,j] <- quantile( PsA_PsA_MenA_002_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	PsA_Hib_MenA_002_quant[,j] <- quantile( PsA_Hib_MenA_002_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	Hib_PsA_MenA_002_quant[,j] <- quantile( Hib_PsA_MenA_002_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )

	PsA_PsA_logelisa_002_quant[,j] <- quantile( PsA_PsA_logelisa_002_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	PsA_Hib_logelisa_002_quant[,j] <- quantile( PsA_Hib_logelisa_002_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
	Hib_PsA_logelisa_002_quant[,j] <- quantile( Hib_PsA_logelisa_002_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975) )
}



PsA_PsA_MenA_002_95CI = matrix(NA, nrow=3, ncol=N_tt)
PsA_Hib_MenA_002_95CI = matrix(NA, nrow=3, ncol=N_tt)
Hib_PsA_MenA_002_95CI = matrix(NA, nrow=3, ncol=N_tt)

PsA_PsA_logelisa_002_95CI = matrix(NA, nrow=3, ncol=N_tt)
PsA_Hib_logelisa_002_95CI = matrix(NA, nrow=3, ncol=N_tt)
Hib_PsA_logelisa_002_95CI = matrix(NA, nrow=3, ncol=N_tt)

for(j in 2:N_tt)
{
	PsA_PsA_MenA_temp = PsA_PsA_MenA_002_mat[,j]
	PsA_PsA_MenA_temp[which(PsA_PsA_MenA_temp < 2, arr.ind=TRUE)] = 2

	PsA_PsA_MenA_002_95CI[1,j] = exp(mean(log(PsA_PsA_MenA_temp)))
	PsA_PsA_MenA_002_95CI[2,j] = exp(t.test( log(PsA_PsA_MenA_temp) )$conf.in[1])
	PsA_PsA_MenA_002_95CI[3,j] = exp(t.test( log(PsA_PsA_MenA_temp) )$conf.in[2])


	PsA_Hib_MenA_temp = PsA_Hib_MenA_002_mat[,j]
	PsA_Hib_MenA_temp[which(PsA_Hib_MenA_temp < 2, arr.ind=TRUE)] = 2

	PsA_Hib_MenA_002_95CI[1,j] = exp(mean(log(PsA_Hib_MenA_temp)))
	PsA_Hib_MenA_002_95CI[2,j] = exp(t.test( log(PsA_Hib_MenA_temp) )$conf.in[1])
	PsA_Hib_MenA_002_95CI[3,j] = exp(t.test( log(PsA_Hib_MenA_temp) )$conf.in[2])


	Hib_PsA_MenA_temp = Hib_PsA_MenA_002_mat[,j]
	Hib_PsA_MenA_temp[which(Hib_PsA_MenA_temp < 2, arr.ind=TRUE)] = 2

	Hib_PsA_MenA_002_95CI[1,j] = exp(mean(log(Hib_PsA_MenA_temp)))
	Hib_PsA_MenA_002_95CI[2,j] = exp(t.test( log(Hib_PsA_MenA_temp) )$conf.in[1])
	Hib_PsA_MenA_002_95CI[3,j] = exp(t.test( log(Hib_PsA_MenA_temp) )$conf.in[2])


	PsA_PsA_logelisa_temp = PsA_PsA_logelisa_002_mat[,j]
	PsA_PsA_logelisa_temp[which(PsA_PsA_logelisa_temp < 0.1, arr.ind=TRUE)] = 0.1

	PsA_PsA_logelisa_002_95CI[1,j] = exp(mean(log(PsA_PsA_logelisa_temp)))
	PsA_PsA_logelisa_002_95CI[2,j] = exp(t.test( log(PsA_PsA_logelisa_temp) )$conf.in[1])
	PsA_PsA_logelisa_002_95CI[3,j] = exp(t.test( log(PsA_PsA_logelisa_temp) )$conf.in[2])


	PsA_Hib_logelisa_temp = PsA_Hib_logelisa_002_mat[,j]
	PsA_Hib_logelisa_temp[which(PsA_Hib_logelisa_temp < 0.1, arr.ind=TRUE)] = 0.1

	PsA_Hib_logelisa_002_95CI[1,j] = exp(mean(log(PsA_Hib_logelisa_temp)))
	PsA_Hib_logelisa_002_95CI[2,j] = exp(t.test( log(PsA_Hib_logelisa_temp) )$conf.in[1])
	PsA_Hib_logelisa_002_95CI[3,j] = exp(t.test( log(PsA_Hib_logelisa_temp) )$conf.in[2])


	Hib_PsA_logelisa_temp = Hib_PsA_logelisa_002_mat[,j]
	Hib_PsA_logelisa_temp[which(Hib_PsA_logelisa_temp < 0.1, arr.ind=TRUE)] = 0.1

	Hib_PsA_logelisa_002_95CI[1,j] = exp(mean(log(Hib_PsA_logelisa_temp)))
	Hib_PsA_logelisa_002_95CI[2,j] = exp(t.test( log(Hib_PsA_logelisa_temp) )$conf.in[1])
	Hib_PsA_logelisa_002_95CI[3,j] = exp(t.test( log(Hib_PsA_logelisa_temp) )$conf.in[2])
}


#################################
#################################
##                             ##
##  INDIVIDUAL SIMULATIONS     ##
##                             ##
################################# 
#################################




Gambia_PsA_PsA_002_index = intersect( which(AB_data_002$country=="Gambia"), intersect( which(AB_data_002$prim_vac=="PsA-TT"), which(AB_data_002$boost_vac=="PsA-TT") ) ) 

Gambia_PsA_Hib_002_index = intersect( which(AB_data_002$country=="Gambia"), intersect( which(AB_data_002$prim_vac=="PsA-TT"), which(AB_data_002$boost_vac=="Hib-TT") ) ) 

Gambia_Hib_PsA_002_index = intersect( which(AB_data_002$country=="Gambia"), intersect( which(AB_data_002$prim_vac=="Hib-TT"), which(AB_data_002$boost_vac=="PsA-TT") ) ) 

Mali_PsA_PsA_002_index = intersect( which(AB_data_002$country=="Mali"), intersect( which(AB_data_002$prim_vac=="PsA-TT"), which(AB_data_002$boost_vac=="PsA-TT") ) ) 

Mali_PsA_Hib_002_index = intersect( which(AB_data_002$country=="Mali"), intersect( which(AB_data_002$prim_vac=="PsA-TT"), which(AB_data_002$boost_vac=="Hib-TT") ) ) 

Mali_Hib_PsA_002_index = intersect( which(AB_data_002$country=="Mali"), intersect( which(AB_data_002$prim_vac=="Hib-TT"), which(AB_data_002$boost_vac=="PsA-TT") ) ) 





PsA_PsA_002_1_MenA     = MenA_002_pred( Gambia_PsA_PsA_002_index[3] )
PsA_PsA_002_1_logelisa = logelisa_002_pred( Gambia_PsA_PsA_002_index[3] )

PsA_PsA_002_2_MenA     = MenA_002_pred( Mali_PsA_PsA_002_index[4] )
PsA_PsA_002_2_logelisa = logelisa_002_pred( Mali_PsA_PsA_002_index[4] )


PsA_Hib_002_1_MenA     = MenA_002_pred( Gambia_PsA_Hib_002_index[5] )
PsA_Hib_002_1_logelisa = logelisa_002_pred( Gambia_PsA_Hib_002_index[5] )

PsA_Hib_002_2_MenA     = MenA_002_pred( Mali_PsA_Hib_002_index[1] )
PsA_Hib_002_2_logelisa = logelisa_002_pred( Mali_PsA_Hib_002_index[1] )


Hib_PsA_002_1_MenA     = MenA_002_pred( Gambia_Hib_PsA_002_index[1] )
Hib_PsA_002_1_logelisa = logelisa_002_pred( Gambia_Hib_PsA_002_index[1] )

Hib_PsA_002_2_MenA     = MenA_002_pred( Mali_Hib_PsA_002_index[1] )
Hib_PsA_002_2_logelisa = logelisa_002_pred( Mali_Hib_PsA_002_index[1] )







#################################
#################################
##                             ##
##  VE thresholds              ##
##                             ##
################################# 
#################################


PsA_PsA_MenA_002_VE128 = matrix(NA, nrow=N_PsA_PsA_002, ncol=N_tt)
PsA_Hib_MenA_002_VE128 = matrix(NA, nrow=N_PsA_Hib_002, ncol=N_tt)
Hib_PsA_MenA_002_VE128 = matrix(NA, nrow=N_Hib_PsA_002, ncol=N_tt)

PsA_PsA_logelisa_002_VE128 = matrix(NA, nrow=N_PsA_PsA_002, ncol=N_tt)
PsA_Hib_logelisa_002_VE128 = matrix(NA, nrow=N_PsA_Hib_002, ncol=N_tt)
Hib_PsA_logelisa_002_VE128 = matrix(NA, nrow=N_Hib_PsA_002, ncol=N_tt)


for(i in 1:N_PsA_PsA_002)
{
	PsA_PsA_MenA_002_VE128[i,]     = as.numeric(PsA_PsA_MenA_002_mat[i,] > 128)
	PsA_PsA_logelisa_002_VE128[i,] = as.numeric(PsA_PsA_logelisa_002_mat[i,] > 128)
}


for(i in 1:N_PsA_Hib_002)
{
	PsA_Hib_MenA_002_VE128[i,]     = as.numeric(PsA_Hib_MenA_002_mat[i,] > 128)
	PsA_Hib_logelisa_002_VE128[i,] = as.numeric(PsA_Hib_logelisa_002_mat[i,] > 128)
}


for(i in 1:N_Hib_PsA_002)
{
	Hib_PsA_MenA_002_VE128[i,]     = as.numeric(Hib_PsA_MenA_002_mat[i,] > 128)
	Hib_PsA_logelisa_002_VE128[i,] = as.numeric(Hib_PsA_logelisa_002_mat[i,] > 128)
}





PsA_PsA_MenA_002_VE128 = apply(X=PsA_PsA_MenA_002_VE128, MARGIN=2, FUN=mean)
PsA_Hib_MenA_002_VE128 = apply(X=PsA_Hib_MenA_002_VE128, MARGIN=2, FUN=mean)
Hib_PsA_MenA_002_VE128 = apply(X=Hib_PsA_MenA_002_VE128, MARGIN=2, FUN=mean)



PsA_PsA_logelisa_002_VE128 = apply(X=PsA_PsA_logelisa_002_VE128, MARGIN=2, FUN=mean)
PsA_Hib_logelisa_002_VE128 = apply(X=PsA_Hib_logelisa_002_VE128, MARGIN=2, FUN=mean)
Hib_PsA_logelisa_002_VE128 = apply(X=Hib_PsA_logelisa_002_VE128, MARGIN=2, FUN=mean)





PsA_PsA_MenA_002_VE1024 = matrix(NA, nrow=N_PsA_PsA_002, ncol=N_tt)
PsA_Hib_MenA_002_VE1024 = matrix(NA, nrow=N_PsA_Hib_002, ncol=N_tt)
Hib_PsA_MenA_002_VE1024 = matrix(NA, nrow=N_Hib_PsA_002, ncol=N_tt)

PsA_PsA_logelisa_002_VE1024 = matrix(NA, nrow=N_PsA_PsA_002, ncol=N_tt)
PsA_Hib_logelisa_002_VE1024 = matrix(NA, nrow=N_PsA_Hib_002, ncol=N_tt)
Hib_PsA_logelisa_002_VE1024 = matrix(NA, nrow=N_Hib_PsA_002, ncol=N_tt)




for(i in 1:N_PsA_PsA_002)
{
	PsA_PsA_MenA_002_VE1024[i,]     = as.numeric(PsA_PsA_MenA_002_mat[i,] > 1024)
	PsA_PsA_logelisa_002_VE1024[i,] = as.numeric(PsA_PsA_logelisa_002_mat[i,] > 1024)
}


for(i in 1:N_PsA_Hib_002)
{
	PsA_Hib_MenA_002_VE1024[i,]     = as.numeric(PsA_Hib_MenA_002_mat[i,] > 1024)
	PsA_Hib_logelisa_002_VE1024[i,] = as.numeric(PsA_Hib_logelisa_002_mat[i,] > 1024)
}


for(i in 1:N_Hib_PsA_002)
{
	Hib_PsA_MenA_002_VE1024[i,]     = as.numeric(Hib_PsA_MenA_002_mat[i,] > 1024)
	Hib_PsA_logelisa_002_VE1024[i,] = as.numeric(Hib_PsA_logelisa_002_mat[i,] > 1024)
}



PsA_PsA_MenA_002_VE1024 = apply(X=PsA_PsA_MenA_002_VE1024, MARGIN=2, FUN=mean)
PsA_Hib_MenA_002_VE1024 = apply(X=PsA_Hib_MenA_002_VE1024, MARGIN=2, FUN=mean)
Hib_PsA_MenA_002_VE1024 = apply(X=Hib_PsA_MenA_002_VE1024, MARGIN=2, FUN=mean)



PsA_PsA_logelisa_002_VE1024 = apply(X=PsA_PsA_logelisa_002_VE1024, MARGIN=2, FUN=mean)
PsA_Hib_logelisa_002_VE1024 = apply(X=PsA_Hib_logelisa_002_VE1024, MARGIN=2, FUN=mean)
Hib_PsA_logelisa_002_VE1024 = apply(X=Hib_PsA_logelisa_002_VE1024, MARGIN=2, FUN=mean)


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

	OUTPUT <- list()

	OUTPUT[[1]] = AB_quant
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

	OUTPUT <- list()

	OUTPUT[[1]] = ELISA_axis_transform( AB_quant )
	OUTPUT[[2]] = ELISA_axis_transform( AB )
	OUTPUT[[3]] = tt

	OUTPUT
}







N_PsA_003 = nrow(AB_data_003)




tt_plot = seq(from=0, to=5.2*365, by=1)
log2 = log(2)
N_tt = length(tt_plot)





MenA_003_mat  = matrix(NA, nrow=N_PsA_003, ncol=N_tt)


logelisa_003_mat  = matrix(NA, nrow=N_PsA_003, ncol=N_tt)




for(i in 1:N_PsA_003)
{
	MenA_003_mat[i,]     = MenA_003_pred( i )[[1]][2,]
	logelisa_003_mat[i,] = logelisa_003_pred( i )[[1]][2,]
}




MenA_003_quant  = matrix(NA, nrow=5, ncol=N_tt)

logelisa_003_quant  = matrix(NA, nrow=5, ncol=N_tt)



for(j in 1:N_tt)
{
	MenA_003_quant[,j]  <- quantile( MenA_003_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )

	logelisa_003_quant[,j]  <- quantile( logelisa_003_mat[,j], prob=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE )
}








MenA_003_95CI = matrix(NA, nrow=3, ncol=N_tt)

logelisa_003_95CI = matrix(NA, nrow=3, ncol=N_tt)

for(j in 2:N_tt)
{
	MenA_temp = MenA_003_mat[,j]
	MenA_temp[which(MenA_temp < 2, arr.ind=TRUE)] = 2


	if( sd(MenA_003_mat[,j]) > 0 )
	{
		MenA_003_95CI[1,j] = exp(mean(log(MenA_temp)))
		MenA_003_95CI[2,j] = exp(t.test( log(MenA_temp) )$conf.in[1])
		MenA_003_95CI[3,j] = exp(t.test( log(MenA_temp) )$conf.in[2])
	}else{
		MenA_003_95CI[1,j] = mean(MenA_temp)
		MenA_003_95CI[2,j] = mean(MenA_temp)
		MenA_003_95CI[3,j] = mean(MenA_temp)
	}


	logelisa_temp = logelisa_003_mat[,j]
	logelisa_temp[which(logelisa_temp < 0.1, arr.ind=TRUE)] = 0.1


	logelisa_003_95CI[1,j] = exp(mean(log(logelisa_temp)))
	logelisa_003_95CI[2,j] = exp(t.test( log(logelisa_temp) )$conf.in[1])
	logelisa_003_95CI[3,j] = exp(t.test( log(logelisa_temp) )$conf.in[2])
}



#################################
#################################
##                             ##
##  INDIVIDUAL SIMULATIONS     ##
##                             ##
################################# 
#################################




Gambia_003_index = which(AB_data_003$country=="Gambia") 

Mali_003_index = which(AB_data_003$country=="Mali") 

Senegal_003_index = which(AB_data_003$country=="Senegal") 




PsA_003_1_MenA     = MenA_003_pred( Gambia_003_index[3] )
PsA_003_1_logelisa = logelisa_003_pred( Gambia_003_index[3] )


PsA_003_2_MenA     = MenA_003_pred( Senegal_003_index[1] )
PsA_003_2_logelisa = logelisa_003_pred( Senegal_003_index[1] )





#################################
#################################
##                             ##
##  VE thresholds              ##
##                             ##
################################# 
#################################


MenA_003_VE128 = matrix(NA, nrow=N_PsA_003, ncol=N_tt)


logelisa_003_VE128 = matrix(NA, nrow=N_PsA_003, ncol=N_tt)




for(i in 1:N_PsA_003)
{
	MenA_003_VE128[i,]     = as.numeric(MenA_003_mat[i,] > 128)
	logelisa_003_VE128[i,] = as.numeric(logelisa_003_mat[i,] > 128)
}



MenA_003_VE128  = apply(X=MenA_003_VE128,  MARGIN=2, FUN=mean, na.rm=TRUE)

logelisa_003_VE128  = apply(X=logelisa_003_VE128,  MARGIN=2, FUN=mean, na.rm=TRUE)







MenA_003_VE1024  = matrix(NA, nrow=N_PsA_003, ncol=N_tt)

logelisa_003_VE1024  = matrix(NA, nrow=N_PsA_003, ncol=N_tt)



for(i in 1:N_PsA_003)
{
	MenA_003_VE1024[i,]     = as.numeric(MenA_003_mat[i,] > 1024)
	logelisa_003_VE1024[i,] = as.numeric(logelisa_003_mat[i,] > 1024)
}


MenA_003_VE1024  = apply(X=MenA_003_VE1024,  MARGIN=2, FUN=mean, na.rm=TRUE)

logelisa_003_VE1024  = apply(X=logelisa_003_VE1024,  MARGIN=2, FUN=mean, na.rm=TRUE)





PsA_PsA_MenA_002_GMT  = rep(NA, N_tt)
Hib_PsA_MenA_002_GMT  = rep(NA, N_tt)
PsA_Hib_MenA_002_GMT  = rep(NA, N_tt)

PsA_PsA_logelisa_002_GMT  = rep(NA, N_tt)
Hib_PsA_logelisa_002_GMT  = rep(NA, N_tt)
PsA_Hib_logelisa_002_GMT  = rep(NA, N_tt)


PsA_PsA_MenA_002_mat[which(PsA_PsA_MenA_002_mat < 2, arr.ind=TRUE)] = 2
PsA_Hib_MenA_002_mat[which(PsA_Hib_MenA_002_mat < 2, arr.ind=TRUE)] = 2
Hib_PsA_MenA_002_mat[which(Hib_PsA_MenA_002_mat < 2, arr.ind=TRUE)] = 2

PsA_PsA_logelisa_002_mat[which(PsA_PsA_logelisa_002_mat < 0.1, arr.ind=TRUE)] = 0.1
PsA_Hib_logelisa_002_mat[which(PsA_Hib_logelisa_002_mat < 0.1, arr.ind=TRUE)] = 0.1
Hib_PsA_logelisa_002_mat[which(Hib_PsA_logelisa_002_mat < 0.1, arr.ind=TRUE)] = 0.1

for(j in 1:N_tt)
{
	PsA_PsA_MenA_002_GMT[j] = exp(mean(log(PsA_PsA_MenA_002_mat[,j]), na.rm=TRUE))
	Hib_PsA_MenA_002_GMT[j] = exp(mean(log(Hib_PsA_MenA_002_mat[,j]), na.rm=TRUE))
	PsA_Hib_MenA_002_GMT[j] = exp(mean(log(PsA_Hib_MenA_002_mat[,j]), na.rm=TRUE))

	PsA_PsA_logelisa_002_GMT[j] = exp(mean(log(PsA_PsA_logelisa_002_mat[,j]), na.rm=TRUE))
	Hib_PsA_logelisa_002_GMT[j] = exp(mean(log(Hib_PsA_logelisa_002_mat[,j]), na.rm=TRUE))
	PsA_Hib_logelisa_002_GMT[j] = exp(mean(log(PsA_Hib_logelisa_002_mat[,j]), na.rm=TRUE))
}



MenA_003_GMT  = rep(NA, N_tt)

logelisa_003_GMT  = rep(NA, N_tt)


MenA_003_mat[which(MenA_003_mat < 2, arr.ind=TRUE)] = 2


logelisa_003_mat[which(logelisa_003_mat < 0.1, arr.ind=TRUE)] = 0.1


for(j in 1:N_tt)
{
	MenA_003_GMT[j] = exp(mean(log(MenA_003_mat[,j]), na.rm=TRUE))

	logelisa_003_GMT[j] = exp(mean(log(logelisa_003_mat[,j]), na.rm=TRUE))
}




###################################
###################################
##                               ##
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
###################################

Diallo_SBA = matrix(NA, nrow=5, ncol=4)
colnames(Diallo_SBA) = c("time", "GMC", "GMC_low", "GMC_high")

Diallo_SBA[,1] = c( 0, 4*7, 26*7, 52*7, 4*52*7)
Diallo_SBA[,2] = c( 223.5, 4712.6, 3040.8, 2889.4, 2281.1 )
Diallo_SBA[,3] = c( 181.3, 4336.0, 2758.6, 2643.3, 1910.5 )
Diallo_SBA[,4] = c( 274.9, 5122.0, 3351.8, 3158.5, 2723.5 )


Diallo_IgG = matrix(NA, nrow=5, ncol=4)
colnames(Diallo_IgG) = c("time", "GMC", "GMC_low", "GMC_high")

Diallo_IgG[,1] = c( 0, 4*7, 26*7, 52*7, 4*52*7)
Diallo_IgG[,2] = c( 2.1, 65.6, 22.3, 15.3, 12.1 )
Diallo_IgG[,3] = c( 1.9, 60.0, 20.1, 13.7, 10.5 )
Diallo_IgG[,4] = c( 2.5, 71.6, 24.7, 17.0, 14.1 )


Tapia_SBA_PsAPsA = matrix(NA, nrow=7, ncol=4)
colnames(Tapia_SBA_PsAPsA) = c("time", "GMC", "GMC_low", "GMC_high")

Tapia_SBA_PsAPsA[,1] = c( 0, 4*7, 38*7, 42*7, 52*7, 104*7, 5*52*7 )
Tapia_SBA_PsAPsA[,2] = c( 14.3, 6234.5, 1130.6, 10037.4, 4485.8, 2720.8, 3781.2 )
Tapia_SBA_PsAPsA[,3] = c( 9.9, 4947.9, 666.0, 7884.5, 3579.6, 1960.1, 2410.4 )
Tapia_SBA_PsAPsA[,4] = c( 20.7, 7855.7, 1919.2, 12778.2, 5621.4, 3776.8, 5931.4 )


Tapia_SBA_PsAHib = matrix(NA, nrow=7, ncol=4)
colnames(Tapia_SBA_PsAHib) = c("time", "GMC", "GMC_low", "GMC_high")

Tapia_SBA_PsAHib[,1] = c( 0, 4*7, 38*7, 42*7, 52*7, 104*7, 5*52*7 )
Tapia_SBA_PsAHib[,2] = c( 14.3, 6234.5, 801.3, 1649.1, 1035.0, 1313.7, 2363.8  )
Tapia_SBA_PsAHib[,3] = c( 9.9, 4947.9, 457.6, 1022.7, 678.4, 875.1, 1256.7 )
Tapia_SBA_PsAHib[,4] = c( 20.7, 7855.7, 1403.1, 2659.3, 1578.9, 1971.9, 4446.1 )


Tapia_SBA_HibPsA = matrix(NA, nrow=7, ncol=4)
colnames(Tapia_SBA_HibPsA) = c("time", "GMC", "GMC_low", "GMC_high")

Tapia_SBA_HibPsA[,1] = c( 0,    4*7,  38*7, 42*7,    52*7,   104*7,  5*52*7 )
Tapia_SBA_HibPsA[,2] = c( 12.6, 60.9, 42.6, 9342.9,  3722.5, 2527.6, 2681.7  )
Tapia_SBA_HibPsA[,3] = c( 8.7,  39.8, 20.3, 7043.8,  2551.4, 1796.9, 2325.0 )
Tapia_SBA_HibPsA[,4] = c( 18.2, 93.2, 89.4, 12392.4, 5431.3, 3555.5, 5830.0 )




Tapia_IgG_PsAPsA = matrix(NA, nrow=7, ncol=4)
colnames(Tapia_IgG_PsAPsA) = c("time", "GMC", "GMC_low", "GMC_high")

Tapia_IgG_PsAPsA[,1] = c( 0, 4*7, 38*7, 42*7, 52*7, 104*7, 5*52*7 )
Tapia_IgG_PsAPsA[,2] = c( 0.1, 18.2, 1.0, 38.1, 7.3, 3.7, 2.0 )
Tapia_IgG_PsAPsA[,3] = c( 0.1, 16.0, 0.8, 25.5, 5.4, 2.7, 1.4 )
Tapia_IgG_PsAPsA[,4] = c( 0.1, 20.7, 1.4, 57.2, 9.8, 5.0, 2.8 )


Tapia_IgG_PsAHib = matrix(NA, nrow=7, ncol=4)
colnames(Tapia_IgG_PsAHib) = c("time", "GMC", "GMC_low", "GMC_high")

Tapia_IgG_PsAHib[,1] = c( 0, 4*7, 38*7, 42*7, 52*7, 104*7, 5*52*7 )
Tapia_IgG_PsAHib[,2] = c( 0.1, 18.2, 1.0, 1.1, 0.9, 0.9, 1.1  )
Tapia_IgG_PsAHib[,3] = c( 0.1, 16.0, 0.7, 0.7, 0.6, 0.7, 0.7 )
Tapia_IgG_PsAHib[,4] = c( 0.1, 20.7, 1.4, 1.6, 1.3, 1.3, 1.6 )


Tapia_IgG_HibPsA = matrix(NA, nrow=7, ncol=4)
colnames(Tapia_IgG_HibPsA) = c("time", "GMC", "GMC_low", "GMC_high")

Tapia_IgG_HibPsA[,1] = c( 0,    4*7,  38*7, 42*7,    52*7,   104*7,  5*52*7 )
Tapia_IgG_HibPsA[,2] = c( 0.1, 0.1, 0.1, 15.4, 1.6, 1.2, 1.2  )
Tapia_IgG_HibPsA[,3] = c( 0.1, 0.1, 0.1, 11.7, 1.2, 0.9, 0.9 )
Tapia_IgG_HibPsA[,4] = c( 0.2, 0.1, 0.2, 20.2, 2.1, 1.7, 1.6 )







##################################
##################################
##                              ##
##  #####  ##     ####  ######  ##
##  ##  ## ##    ##  ##   ##    ##
##  #####  ##    ##  ##   ##    ## 
##  ##     ##    ##  ##   ##    ## 
##  ##     #####  ####    ##    ##
##                              ##
################################## 
##################################




line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)


line_seq_y <- 365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)




main.size = 0.85
axis.size = 0.5
lab.size  = 1

tiff( file="1_Fig1_Ab_dynamics.tif", width=25, height=20, units="cm", res=500)

lay.mat = rbind( c( 1,  2,  3,  4  ),
                 c( 5,  7,  9,  17 ),
                 c( 6,  8,  10, 18 ),
                 c( 11, 12, 13, 19 ),
                 c( 14, 15, 16, 20 ),
                 c( 21, 21, 21, 21 ) )

layout(lay.mat, heights=c(1,10,10,10,10,1))
layout.show(21)






############################
## Labels on top          ## 
############################

par(mar = c(0,0,0,0))


plot.new()
title( "PsA-TT-002: MenAfriVac/MenAfriVac", cex.main=1.0, line=-1)

plot.new()
title( "PsA-TT-002: MenAfriVac/Hib-TT", cex.main=1.0, line=-1)

plot.new()
title( "PsA-TT-002: Hib-TT/MenAfriVac", cex.main=1.0, line=-1)


plot.new()
title( "PsA-TT-003: MenAfriVac", cex.main=1.0, line=-1)



par(mar = c(2.5,2,1,1.75))
par(mgp=c(1.3,0.5,0))




############################
##                        ## 
##  PANEL 1               ##
##  PsA_PsA_002_1         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(A) Individual G002_1; age = ", round(AB_data_002[Gambia_PsA_PsA_002_index[3],4]/12,2), " years", sep=""),
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
	  y=c( PsA_PsA_002_1_logelisa[[1]][1,], rev(PsA_PsA_002_1_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_PsA_002_1_MenA[[1]][1,], rev(PsA_PsA_002_1_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)



points(x=tt_plot, y=PsA_PsA_002_1_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_PsA_002_1_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")


points(x=PsA_PsA_002_1_logelisa[[3]], y=PsA_PsA_002_1_logelisa[[2]], pch=19, col="dodgerblue")

for(j in 1:length(PsA_PsA_002_1_MenA[[2]]))
{
	arrows(y0=PsA_PsA_002_1_MenA[[2]][j], x0=PsA_PsA_002_1_MenA[[3]][j], 
    	       y1=0.5*PsA_PsA_002_1_MenA[[2]][j], x1=PsA_PsA_002_1_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}



text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=430, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))


############################
##                        ## 
##  PANEL 2               ##
##  PsA_PsA_002_2         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(E) Individual M002_1; age = ", round(AB_data_002[Mali_PsA_PsA_002_index[4],4]/12,2), " years", sep=""),
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
	  y=c( PsA_PsA_002_2_logelisa[[1]][1,], rev(PsA_PsA_002_2_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_PsA_002_2_MenA[[1]][1,], rev(PsA_PsA_002_2_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)


points(x=tt_plot, y=PsA_PsA_002_2_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_PsA_002_2_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")


points(x=PsA_PsA_002_2_logelisa[[3]], y=PsA_PsA_002_2_logelisa[[2]], pch=19, col="dodgerblue")


for(j in 1:length(PsA_PsA_002_2_MenA[[2]]))
{
	arrows(y0=PsA_PsA_002_2_MenA[[2]][j], x0=PsA_PsA_002_2_MenA[[3]][j], 
    	       y1=0.5*PsA_PsA_002_2_MenA[[2]][j], x1=PsA_PsA_002_2_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}


text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=430, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))








############################
##                        ## 
##  PANEL 3               ##
##  PsA_Hib_002_1         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(B) Individual G002_3; age = ", round(AB_data_002[Gambia_PsA_Hib_002_index[5],4]/12,2), " years", sep=""),
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
	  y=c( PsA_Hib_002_1_logelisa[[1]][1,], rev(PsA_Hib_002_1_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_Hib_002_1_MenA[[1]][1,], rev(PsA_Hib_002_1_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

points(x=tt_plot, y=PsA_Hib_002_1_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_Hib_002_1_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")

points(x=PsA_Hib_002_1_logelisa[[3]], y=PsA_Hib_002_1_logelisa[[2]], pch=19, col="dodgerblue")


for(j in 1:length(PsA_Hib_002_1_MenA[[2]]))
{
	arrows(y0=PsA_Hib_002_1_MenA[[2]][j], x0=PsA_Hib_002_1_MenA[[3]][j], 
    	       y1=0.5*PsA_Hib_002_1_MenA[[2]][j], x1=PsA_Hib_002_1_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}




text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=410, y=600000, labels="Hib-TT", cex=0.7 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))






############################
##                        ## 
##  PANEL 4               ##
##  PsA_Hib_002_2         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(F) Individual M002_3; age = ", round(AB_data_002[Mali_PsA_Hib_002_index[1],4]/12,2), " years", sep=""),
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
	  y=c( PsA_Hib_002_2_logelisa[[1]][1,], rev(PsA_Hib_002_2_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_Hib_002_2_MenA[[1]][1,], rev(PsA_Hib_002_2_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

points(x=tt_plot, y=PsA_Hib_002_2_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_Hib_002_2_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")

points(x=PsA_Hib_002_2_logelisa[[3]], y=PsA_Hib_002_2_logelisa[[2]], pch=19, col="dodgerblue")


for(j in 1:length(PsA_Hib_002_2_MenA[[2]]))
{
	arrows(y0=PsA_Hib_002_2_MenA[[2]][j], x0=PsA_Hib_002_2_MenA[[3]][j], 
    	       y1=0.5*PsA_Hib_002_2_MenA[[2]][j], x1=PsA_Hib_002_2_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}


text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=410, y=600000, labels="Hib-TT", cex=0.7 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))






############################
##                        ## 
##  PANEL 5               ##
##  Hib_PsA_002_1         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(C) Individual G002_5; age = ", round(AB_data_002[Gambia_Hib_PsA_002_index[1],4]/12,2), " years", sep=""),
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
	  y=c( Hib_PsA_002_1_logelisa[[1]][1,], rev(Hib_PsA_002_1_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( Hib_PsA_002_1_MenA[[1]][1,], rev(Hib_PsA_002_1_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

points(x=tt_plot, y=Hib_PsA_002_1_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=Hib_PsA_002_1_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")

points(x=Hib_PsA_002_1_logelisa[[3]], y=Hib_PsA_002_1_logelisa[[2]], pch=19, col="dodgerblue")

for(j in 1:length(Hib_PsA_002_1_MenA[[2]]))
{
	arrows(y0=Hib_PsA_002_1_MenA[[2]][j], x0=Hib_PsA_002_1_MenA[[3]][j], 
    	       y1=0.5*Hib_PsA_002_1_MenA[[2]][j], x1=Hib_PsA_002_1_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}





text(x=120, y=600000, labels="Hib-TT", cex=0.7 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=430, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))






############################
##                        ## 
##  PANEL 6               ##
##  Hib_PsA_002_2         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(G) Individual M002_5; age = ", round(AB_data_002[Mali_Hib_PsA_002_index[1],4]/12,2), " years", sep=""),
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
	  y=c( Hib_PsA_002_2_logelisa[[1]][1,], rev(Hib_PsA_002_2_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( Hib_PsA_002_2_MenA[[1]][1,], rev(Hib_PsA_002_2_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

points(x=tt_plot, y=Hib_PsA_002_2_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=Hib_PsA_002_2_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")

points(x=Hib_PsA_002_2_logelisa[[3]], y=Hib_PsA_002_2_logelisa[[2]], pch=19, col="dodgerblue")


for(j in 1:length(Hib_PsA_002_2_MenA[[2]]))
{
	arrows(y0=Hib_PsA_002_2_MenA[[2]][j], x0=Hib_PsA_002_2_MenA[[3]][j], 
    	       y1=0.5*Hib_PsA_002_2_MenA[[2]][j], x1=Hib_PsA_002_2_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}



text(x=120, y=600000, labels="Hib-TT", cex=0.7 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=430, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))










############################
##                        ## 
##  PANEL 7               ##
##  PsA_PsA_002_1         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main="(I) Population average (n = 64)",
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
	  y=c( PsA_PsA_logelisa_002_95CI[2,], rev(PsA_PsA_logelisa_002_95CI[3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_PsA_MenA_002_95CI[2,], rev(PsA_PsA_MenA_002_95CI[3,]) ),
	  col=rgb(255/256,69/256,0/256,0.25), border=NA)

points(x=tt_plot, y=PsA_PsA_logelisa_002_GMT, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_PsA_MenA_002_GMT, type='l', lwd=2, col="firebrick1")




text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=430, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))



for(j in 1:nrow(Tapia_SBA_PsAPsA))
{
	points( x=Tapia_SBA_PsAPsA[j,1], y=Tapia_SBA_PsAPsA[j,2], pch=19, col="firebrick1")

	arrows( y0=Tapia_SBA_PsAPsA[j,3], x0=Tapia_SBA_PsAPsA[j,1], 
    	        y1=Tapia_SBA_PsAPsA[j,4], x1=Tapia_SBA_PsAPsA[j,1], 
     	        length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}



for(j in 1:nrow(Tapia_IgG_PsAPsA))
{
	points( x=Tapia_IgG_PsAPsA[j,1], y=ELISA_axis_transform(Tapia_IgG_PsAPsA[j,2]), pch=19, col="dodgerblue")

	arrows( y0=ELISA_axis_transform(Tapia_IgG_PsAPsA[j,3]), x0=Tapia_IgG_PsAPsA[j,1], 
    	        y1=ELISA_axis_transform(Tapia_IgG_PsAPsA[j,4]), x1=Tapia_IgG_PsAPsA[j,1], 
     	        length=0.03, angle=90, code=3, col="dodgerblue", lwd=1)	
}




############################
##                        ## 
##  PANEL 8               ##
##  PsA_Hib_002_1         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main="(J) Population average (n = 66)",
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
	  y=c( PsA_Hib_logelisa_002_95CI[2,], rev(PsA_Hib_logelisa_002_95CI[3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_Hib_MenA_002_95CI[2,], rev(PsA_Hib_MenA_002_95CI[3,]) ),
	  col=rgb(255/256,69/256,0/256,0.25), border=NA)



points(x=tt_plot, y=PsA_Hib_logelisa_002_GMT, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_Hib_MenA_002_GMT, type='l', lwd=2, col="firebrick1")




text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=410, y=600000, labels="Hib-TT", cex=0.7 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))




for(j in 1:nrow(Tapia_SBA_PsAHib))
{
	points( x=Tapia_SBA_PsAHib[j,1], y=Tapia_SBA_PsAHib[j,2], pch=19, col="firebrick1")

	arrows( y0=Tapia_SBA_PsAHib[j,3], x0=Tapia_SBA_PsAHib[j,1], 
    	        y1=Tapia_SBA_PsAHib[j,4], x1=Tapia_SBA_PsAHib[j,1], 
     	        length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}



for(j in 1:nrow(Tapia_IgG_PsAHib))
{
	points( x=Tapia_IgG_PsAHib[j,1], y=ELISA_axis_transform(Tapia_IgG_PsAHib[j,2]), pch=19, col="dodgerblue")

	arrows( y0=ELISA_axis_transform(Tapia_IgG_PsAHib[j,3]), x0=Tapia_IgG_PsAHib[j,1], 
    	        y1=ELISA_axis_transform(Tapia_IgG_PsAHib[j,4]), x1=Tapia_IgG_PsAHib[j,1], 
     	        length=0.03, angle=90, code=3, col="dodgerblue", lwd=1)	
}









############################
##                        ## 
##  PANEL 9               ##
##  Hib_PsA_002_1         ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main="(K) Population average (n = 63)",
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
	  y=c( Hib_PsA_logelisa_002_95CI[2,], rev(Hib_PsA_logelisa_002_95CI[3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( Hib_PsA_MenA_002_95CI[2,], rev(Hib_PsA_MenA_002_95CI[3,]) ),
	  col=rgb(255/256,69/256,0/256,0.25), border=NA)

points(x=tt_plot, y=Hib_PsA_logelisa_002_GMT, type='l', lwd=2, col="dodgerblue")


points(x=tt_plot, y=Hib_PsA_MenA_002_GMT, type='l', lwd=2, col="firebrick1")




text(x=120, y=600000, labels="Hib-TT", cex=0.7 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

text(x=430, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=275, y0=300000, x1=275, y1=100000, length=0.02)



axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))





for(j in 1:nrow(Tapia_SBA_HibPsA))
{
	points( x=Tapia_SBA_HibPsA[j,1], y=Tapia_SBA_HibPsA[j,2], pch=19, col="firebrick1")

	arrows( y0=Tapia_SBA_HibPsA[j,3], x0=Tapia_SBA_HibPsA[j,1], 
    	        y1=Tapia_SBA_HibPsA[j,4], x1=Tapia_SBA_HibPsA[j,1], 
     	        length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}



for(j in 1:nrow(Tapia_IgG_HibPsA))
{
	points( x=Tapia_IgG_HibPsA[j,1], y=ELISA_axis_transform(Tapia_IgG_HibPsA[j,2]), pch=19, col="dodgerblue")

	arrows( y0=ELISA_axis_transform(Tapia_IgG_HibPsA[j,3]), x0=Tapia_IgG_HibPsA[j,1], 
    	        y1=ELISA_axis_transform(Tapia_IgG_HibPsA[j,4]), x1=Tapia_IgG_HibPsA[j,1], 
     	        length=0.03, angle=90, code=3, col="dodgerblue", lwd=1)	
}











############################
##                        ## 
##  PANEL 10              ##
##  PsA_PsA_002           ##
##                        ##
############################

line_seq_x <- c(0, 0.25, 0.5, 0.75, 1)


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(0,1.02), 
xlab="time (years)", ylab="",
main="(M) Proportion greater than threshold",
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

points(x=tt_plot, y=PsA_PsA_logelisa_002_VE128, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_PsA_MenA_002_VE128, type='l', lwd=2, col="firebrick1")




points(x=tt_plot, y=PsA_PsA_logelisa_002_VE1024, type='l', lwd=2, lty="dashed", col="dodgerblue")


points(x=tt_plot, y=PsA_PsA_MenA_002_VE1024, type='l', lwd=2, lty="dashed", col="firebrick1")

legend(x="bottomright",
col="black", 
lty=c(1,3), lwd=2,
legend=c("128 / 2", "1024 / 8.9"),
bty='n', cex=0.8 )


par(mgp=c(1.3,0.25,0))


axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)



par(mgp=c(1.3,0.5,0))

mtext(side = 2, line = 1, 
cex=axis.size, col="black",
text=expression(paste( "proportion", sep="" )))


axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)






############################
##                        ## 
##  PANEL 11              ##
##  PsA_Hib_002           ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(0,1.02), 
xlab="time (years)", ylab="",
main="(N) Proportion greater than threshold",
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


points(x=tt_plot, y=PsA_Hib_logelisa_002_VE128, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_Hib_MenA_002_VE128, type='l', lwd=2, col="firebrick1")




points(x=tt_plot, y=PsA_Hib_logelisa_002_VE1024, type='l',  lwd=2, lty="dashed", col="dodgerblue")

points(x=tt_plot, y=PsA_Hib_MenA_002_VE1024, type='l', lwd=2, lty="dashed", col="firebrick1")

legend(x="bottomright",
col="black", 
lty=c(1,3), lwd=2,
legend=c("128 / 2", "1024 / 8.9"),
bty='n', cex=0.8 )


par(mgp=c(1.3,0.25,0))


axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)



par(mgp=c(1.3,0.5,0))

mtext(side = 2, line = 1, 
cex=axis.size, col="black",
text=expression(paste( "proportion", sep="" )))


axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)











############################
##                        ## 
##  PANEL 12              ##
##  Hib_PsA_002           ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(0,1.02), 
xlab="time (years)", ylab="",
main="(O) Proportion greater than threshold",
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


points(x=tt_plot, y=Hib_PsA_logelisa_002_VE128, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=Hib_PsA_MenA_002_VE128, type='l', lwd=2, col="firebrick1")




points(x=tt_plot, y=Hib_PsA_logelisa_002_VE1024, type='l', lwd=2, lty="dashed", col="dodgerblue")

points(x=tt_plot, y=Hib_PsA_MenA_002_VE1024, type='l', lwd=2, lty="dashed", col="firebrick1")


legend(x="bottomright",
col="black", 
lty=c(1,3), lwd=2,
legend=c("128 / 2", "1024 / 8.9"),
bty='n', cex=0.8 )


par(mgp=c(1.3,0.25,0))


axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)



par(mgp=c(1.3,0.5,0))

mtext(side = 2, line = 1, 
cex=axis.size, col="black",
text=expression(paste( "proportion", sep="" )))


axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)






############################
##                        ## 
##  PANEL 1               ##
##  PsA_003_1             ##
##                        ##
############################


line_seq_x <- 2^c(-2, 1, 4, 7, 10, 13, 16, 19)


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(D) Individual G003_1; age = ", round(AB_data_003[Gambia_003_index[3],4],2), " years", sep=""),
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
	  y=c( PsA_003_1_logelisa[[1]][1,], rev(PsA_003_1_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_003_1_MenA[[1]][1,], rev(PsA_003_1_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)


points(x=tt_plot, y=PsA_003_1_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_003_1_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")



points(x=PsA_003_1_logelisa[[3]], y=PsA_003_1_logelisa[[2]], pch=19, col="dodgerblue")



for(j in 1:length(PsA_003_1_MenA[[2]]))
{
	arrows(y0=PsA_003_1_MenA[[2]][j], x0=PsA_003_1_MenA[[3]][j], 
    	       y1=0.5*PsA_003_1_MenA[[2]][j], x1=PsA_003_1_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}






text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))




############################
##                        ## 
##  PANEL 2               ##
##  PsA_003_2             ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main=paste( "(H) Individual S003_1; age = ", round(AB_data_003[Senegal_003_index[1],4],2), " years", sep=""),
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
	  y=c( PsA_003_2_logelisa[[1]][1,], rev(PsA_003_2_logelisa[[1]][3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)

polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( PsA_003_2_MenA[[1]][1,], rev(PsA_003_2_MenA[[1]][3,]) ),
	  col=rgb(255/256,48/256,48/256,0.25), border=NA)

points(x=tt_plot, y=PsA_003_2_logelisa[[1]][2,], type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=PsA_003_2_MenA[[1]][2,], type='l', lwd=2, col="firebrick1")

points(x=PsA_003_2_logelisa[[3]], y=PsA_003_2_logelisa[[2]], pch=19, col="dodgerblue")


for(j in 1:length(PsA_003_2_MenA[[2]]))
{
	arrows(y0=PsA_003_2_MenA[[2]][j], x0=PsA_003_2_MenA[[3]][j], 
    	       y1=0.5*PsA_003_2_MenA[[2]][j], x1=PsA_003_2_MenA[[3]][j], 
     	       length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}




text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))








############################
##                        ## 
##  PANEL 7               ##
##  PsA_003               ##
##                        ##
############################


plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(2^(0),2^20), log="y",
xlab="time (years)", ylab="",
main="(L) Population average (n = 604)",
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
	  y=c( logelisa_003_95CI[2,], rev(logelisa_003_95CI[3,]) ),
	  col=rgb(30/256,144/256,255/256,0.25), border=NA)


polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( MenA_003_95CI[2,], rev(MenA_003_95CI[3,]) ),
	  col=rgb(255/256,69/256,0/256,0.25), border=NA)



points(x=tt_plot, y=logelisa_003_GMT, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=MenA_003_GMT, type='l', lwd=2, col="firebrick1")




text(x=120, y=600000, labels="MenAfriVac", cex=0.5 )
arrows(x0=0, y0=300000, x1=0, y1=100000, length=0.02)

axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

axis(2,  at = c( 1,   2,    16,   128,    1024,   8192,    65536,     524288, 1048576), 
     labels = c( "", "2",  "16", "128",  "1024", "8192",  "65536",  "",      "2^20" ), 
cex.axis=0.5)

par(mgp=c(1.3,0.25,0))

axis(4,  at = ELISA_axis_transform( c( 0.1, 1, 10, 100, 1000) ), 
     labels=c( "0.1", "1",  "10", "100",  "1000" ), 
cex.axis=0.5)


mtext(side = 2, line = 1, 
cex=axis.size, col="firebrick1",
text=expression(paste( "SBA titer", sep="" )))

mtext(side = 4, line = 1, 
cex=axis.size, col="dodgerblue",
text=expression(paste( "IgG ELISA (", mu, "g/mL)", sep="" )))






for(j in 1:nrow(Diallo_SBA))
{
	points( x=Diallo_SBA[j,1], y=Diallo_SBA[j,2], pch=19, col="firebrick1")

	arrows( y0=Diallo_SBA[j,3], x0=Diallo_SBA[j,1], 
    	        y1=Diallo_SBA[j,4], x1=Diallo_SBA[j,1], 
     	        length=0.03, angle=90, code=3, col="firebrick1", lwd=1)	
}



for(j in 1:nrow(Diallo_IgG))
{
	points( x=Diallo_IgG[j,1], y=ELISA_axis_transform(Diallo_IgG[j,2]), pch=19, col="dodgerblue")

	arrows( y0=ELISA_axis_transform(Diallo_IgG[j,3]), x0=Diallo_IgG[j,1], 
    	        y1=ELISA_axis_transform(Diallo_IgG[j,4]), x1=Diallo_IgG[j,1], 
     	        length=0.03, angle=90, code=3, col="dodgerblue", lwd=1)	
}





############################
##                        ## 
##  PANEL 10              ##
##  VE_002                ##
##                        ##
############################

line_seq_x <- c(0, 0.25, 0.5, 0.75, 1)



plot( x=1e10, y=1e10,  
xlim=c(-50,1898), ylim=c(0,1.02), 
xlab="time (years)", ylab="",
main="(P) Proportion greater than threshold",
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



points(x=tt_plot, y=logelisa_003_VE128, type='l', lwd=2, col="dodgerblue")

points(x=tt_plot, y=MenA_003_VE128, type='l', lwd=2, col="firebrick1")



points(x=tt_plot, y=logelisa_003_VE1024, type='l', lwd=2, lty="dashed", col="dodgerblue")

points(x=tt_plot, y=MenA_003_VE1024, type='l', lwd=2, lty="dashed", col="firebrick1")


legend(x="bottomright",
col="black", 
lty=c(1,3), lwd=2,
legend=c("128 / 2", "1024 / 8.9"),
bty='n', cex=0.8 )


par(mgp=c(1.3,0.25,0))


axis(1, at=365*c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), labels=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.2), cex.axis=axis.size)


par(mgp=c(1.3,0.5,0))

mtext(side = 2, line = 1, 
cex=axis.size, col="black",
text=expression(paste( "proportion", sep="" )))


axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("0%", "25%", "50%", "75%", "100%"), cex.axis=axis.size)






###############	
##           ##
##  LEGEND   ##
##           ##
###############


legend(x='center', 
       legend = c("SBA titer   ", "IgG ELISA   "), 
       col = c("firebrick1", "dodgerblue"), 
       lty=c(1,1), lwd=2,
       ncol=2, cex=1.5, bty="n" )



dev.off()

	

