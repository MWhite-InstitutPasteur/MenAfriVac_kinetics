#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include "randlib.h"
#include <omp.h>
#include <vector>
#include <algorithm>

using namespace std;


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//          //                                        //
//   ####   //  Setting up objects,                   //
//  ##  ##  //  declaring array sizes                 //
//  ##  ##  //                                        //
//  ##  ##  //                                        //
//   ####   //                                        //
//          //                                        //
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
// 0.1 Declare global variables of sizes of arrays

#define N_mcmc 3000000       // Number of MCMC iterations: indexed by mc

#define N_adapt 600000         // Number of MCMC iterations in adaptive phase
#define N_tune_start 10000    // Number of MCMC iterations in adaptive phase
#define N_tune_end 500000

#define N_data_cols 31       // Number of columns in data frame
#define N_part 193           // Number of participants to read in data from (2071)
#define N_t 7                // Maximum number of data points per participant
#define N_loc_par 8          // Number of individual-level parameters
#define N_glob_par 6         // Number of population-level parameters
#define LARGE 1e12           // large number needed for priors

#define log2 0.69314718055995


/////////////////////////////////////////////////////////////////	
// 0.2 Create structure to hold data for participant n
//     and local parameter estimates

struct part_n
{
	//////////////////////////////////////
	// Covariate information

	int ID;                      // individual's ID
	int country;                 // country (Mali = 1; Gambia = 2)
	double age;                  // age
	int sex;                     // gender (female = 1; male = 2)
	double height;               // height (visit 1)
	double weight;               // weigth (visit 1)


	//////////////////////////////////////
	// Vaccine info

	int prim_vac;                // primary vaccine (PsA-TT = 1; Hib-TT = 2)
	int boost_vac;               // primary vaccine (PsA-TT = 1; Hib-TT = 2)

	double t_prim_vac;           // time of primary vaccine dose
	double t_boost_vac;          // time of booster vaccine dose


	//////////////////////////////////////
	// Antibody data

	int N_sam;                   // Number of samples of antibody data

	vector<double> AB;
	vector<double> tt;

	vector<double> lAB;


	//////////////////////////////////////
	// Individual-level parameters

	double beta_prim;     // initial drug concentration
	double beta_boost;    // initial drug concentration
	double t_short;       // half-life of short-lived component of drugs
	double t_long;        // half-life of short-lived component of drugs
	double t_IgG;         // half-life of short-lived component of drugs
	double rho_prim;      // half-life of short-lived component of drugs
	double rho_boost;     // half-life of short-lived component of drugs
	double Ab_0;

	double lbeta_prim;      // log(initial drug concentration)
	double lbeta_boost;     // log(initial drug concentration)
	double lt_short;        // half-life of short-lived component of drugs
	double lt_long;         // half-life of short-lived component of drugs
	double lt_IgG;          // half-life of short-lived component of drugs
	double logitrho_prim;
	double logitrho_boost;
	double lAb_0;           // log(initial drug concentration)

	double r_short;       // drug decay rate
	double r_long;        // drug decay rate
	double r_IgG;         // drug decay rate



	//////////////////////////////////////
	// Likelihood

	double SSq;           // sum of squares
	double data_like;     // data likelihood
	double mix_like;      // mixed-effects likelihood

	double lhood;         // individual-level likelihood
};


/////////////////////////////////////////////////////////////////
// 0.3 Create structure for global parameters to be estimated

struct params
{
	/////////////////////////////////////////////////////////////
	// Population-level parameters describing mixed effects

	double mu_par[N_glob_par];
	double tau_par[N_glob_par];


	/////////////////////////////////////////////////////////////
	// Parameter for observational error

	double sig_obs;        // precision of observational error
	double log_sig_obs;

	double sig_obs_scale;

	/////////////////////////////////////////////////////////////
	// Log likelihood and prior

	double SSq_total;

	double loglike;
	double prior;


	/////////////////////////////////////////////////////////////
	// Prior distributions

	double prior_MM[N_glob_par];
	double prior_MM_CV[N_glob_par];
	double prior_SIG[N_glob_par];
	double prior_SIG_CV[N_glob_par];

	double prior_LN_MM[N_glob_par];
	double prior_LN_MM_CV[N_glob_par];
	double prior_LN_SIG[N_glob_par];
	double prior_LN_SIG_CV[N_glob_par];

	double prior_mu[N_glob_par];
	double prior_tau[N_glob_par];
	double prior_k[N_glob_par];
	double prior_theta[N_glob_par];


	/////////////////////////////////////////////////////////////
	// Individual-level parameter book-keeping

	double Y_par[N_glob_par];
	double Ymu2_par[N_glob_par];
};


/////////////////////////////////////////////////////////////////
// 0.4 Individual-level structure for MCMC tuning

struct part_n_MCMC
{
	float par_vec[N_loc_par];                             // Parameter vector (in float format for setgmn) (lAB_0, rr)

	float par_vec_test[N_loc_par];                        // Test parameter vector for MCMC update (in float format for setgmn)
	float work[N_loc_par];                                // Scratch vector for setgmn

	double par_S1[N_loc_par];                             // Sum of parameters
	double par_S2[N_loc_par][N_loc_par];                  // Sum of product of pairs

	float COV_MAT[N_loc_par][N_loc_par];                  // covariance matrix (in float format for setgmn)
	float COV_MAT_dummy[N_loc_par][N_loc_par];            // dummy covariance matrix: setgmn gives back sqrt(COV_MAT) or similar so we feed it a dummy

	float GMN_parm[(N_loc_par)*(N_loc_par + 3) / 2 + 1];  // array for setgmn output

	int denom;                                            // denominator for tracking SD calculations

	double step_scale;                                    // scalar for tuning acceptance rate

	int accept;                                           // number of accepted steps
};


////////////////////////////////////////////////////
// 0.5 Initialise functions

double SSq_n(part_n* p);
double data_like_n(part_n* p, params* theta);
double mix_like_n(part_n* p, params* theta);
double global_prior(params* priors);
double local_prior(part_n* p);
double rm_scale(double step_scale, int step, int N_step_adapt, double log_prob);
double gammln(const double xx);



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//        //                                            //
//   ##   //  Initialise main object, read in data and  //
//  ###   //  fill out objects                          //
//   ##   //                                            //
//   ##   //                                            // 
//  ####  //                                            //
//        //                                            //
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// 1.1 Initialise main object - need to choose whether
//     to run in console or as a .exe


int main(int argc, char** argv)
{

	// do we have the correct command line?
	if (argc != 4)
	{
		std::cout << "Incorrect command line.\n";
		return 0;
	}

	char* AB_input_File      = argv[1];
	char* global_output_File = argv[2];
	char* local_output_File  = argv[3];


	//////////////////////////////////////////////////////
	// 1.2 Declare seed, buffer for writing to and clock

	setall(time(NULL), 7);

	int cl = clock();


	///////////////////////////////////////////////////////
	// 1.3 Read in MenAfriVac data 

	std::ifstream AB_Stream(AB_input_File);

	if (AB_Stream.fail())
	{
		std::cout << "Failure reading in data." << endl;
	}

	vector<vector<double>> AB_data_read;

	AB_data_read.resize(N_part);
	for (int i = 0; i < N_part; i++)
	{
		AB_data_read[i].resize(N_data_cols);
	}

	for (int i = 0; i<N_part; i++)
	{
		for (int j = 0; j<N_data_cols; j++)
		{
			AB_Stream >> AB_data_read[i][j];
		}
	}

	AB_Stream.close();


	//////////////////////////////////////////////////////
	// 1.4 Create global parameter objects

	params theta, theta_p1;


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// Priors on global parameter

	////////////////////////////
	// beta_0; initial drug concentration

	theta.prior_MM[0]     = 1000.0;      // Mean value of the parameter in the population
	theta.prior_MM_CV[0]  = 0.25;         // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[0]    = 1000.0;       // Between Person standard deviation
	theta.prior_SIG_CV[0] = 0.33;         // Coefficient of variation in Between Person sd


	////////////////////////////
	// beta_boost; initial drug concentration

	theta.prior_MM[1]     = 1000.0;     // Mean value of the parameter in the population
	theta.prior_MM_CV[1]  = 0.25;        // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[1]    = 1000.0;      // Between Person standard deviation
	theta.prior_SIG_CV[1] = 0.33;        // Coefficient of variation in Between Person sd



	////////////////////////////
	// t_short;  half-life of short-lived component

	theta.prior_MM[2]     = 10.0;      // Mean value of the parameter in the population
	theta.prior_MM_CV[2]  = 0.25;       // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[2]    = 10.0;       // Between Person standard deviation
	theta.prior_SIG_CV[2] = 0.33;       // Coefficient of variation in Between Person sd


	////////////////////////////
	// t_long;  half-life of long-lived component

	theta.prior_MM[3]     = 1000.0;   // Mean value of the parameter in the population
	theta.prior_MM_CV[3]  = 0.25;      // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[3]    = 1000.0;    // Between Person standard deviation
	theta.prior_SIG_CV[3] = 0.33;      // Coefficient of variation in Between Person sd


	////////////////////////////
	// t_IgG;  half-life of IgG

	theta.prior_MM[4]     = 21.0;     // Mean value of the parameter in the population
	theta.prior_MM_CV[4]  = 0.1;     // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[4]    = 5.0;      // Between Person standard deviation
	theta.prior_SIG_CV[4] = 0.33;      // Coefficient of variation in Between Person sd


	////////////////////////////
	// rho_prim;  proportion of short-lived component at start

	theta.prior_MM[5]     = 0.9;      // Mean value of the parameter in the population
	theta.prior_MM_CV[5]  = 0.1;      // Coefficient of variation in estimate of population - level parameter

	theta.prior_SIG[5]    = 0.1;      // Between Person standard deviation
	theta.prior_SIG_CV[5] = 0.25;      // Coefficient of variation in Between Person sd


	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	// Transformation of priors 

	for (int p = 0; p < 5; p++)
	{
		theta.prior_LN_MM[p]     = log( theta.prior_MM[p] / sqrt(1.0 + pow(theta.prior_SIG[p] / theta.prior_MM[p], 2.0)) );   
		theta.prior_LN_MM_CV[p]  = theta.prior_MM_CV[p];

		theta.prior_LN_SIG[p]    = sqrt(log( 1.0 + pow(theta.prior_SIG[p] / theta.prior_MM[p], 2.0) ));          
		theta.prior_LN_SIG_CV[p] = theta.prior_SIG_CV[p];                                                       
	}

	for (int p = 5; p < N_glob_par; p++)
	{
		theta.prior_LN_MM[p] = 2.625251;
		theta.prior_LN_MM_CV[p] = 2.5*theta.prior_MM_CV[p];

		theta.prior_LN_SIG[p] = 1.096994;
		theta.prior_LN_SIG_CV[p] = 2.5*theta.prior_SIG_CV[p];
	}


	for (int p = 0; p < N_glob_par; p++)
	{
		theta.prior_mu[p]  = theta.prior_LN_MM[p];
		theta.prior_tau[p] = 1.0/pow(theta.prior_LN_MM[p]* theta.prior_LN_MM_CV[p], 2.0);

		theta.prior_k[p]     = 1.0/pow( 2.0*theta.prior_LN_SIG_CV[p], 2.0);
		theta.prior_theta[p] = pow( 2.0*theta.prior_LN_SIG_CV[p]/ theta.prior_LN_SIG[p], 2.0);
	}


	for (int p = 0; p < N_glob_par; p++)
	{
		theta.mu_par[p]  = genunf(0.9, 1.1)*theta.prior_mu[p];
		theta.tau_par[p] = 0.1*genunf(0.9, 1.1)*theta.prior_k[p]* theta.prior_theta[p];
	}


	theta.sig_obs = genunf(0.9, 1.1);      // Precision of observational error

	theta.log_sig_obs = log(theta.sig_obs);

	theta.sig_obs_scale = 1.0;


	//////////////////////////////////////////////////////////
	// 1.5 Create individual-level objects for participant n

	part_n* part;
	part = new part_n[N_part];

	for (int n = 0; n<N_part; n++)
	{
		/////////////////////////////////////////////////
		// Fill out covariate data

		part[n].ID      = AB_data_read[n][0];
		part[n].country = AB_data_read[n][1];
		part[n].age     = AB_data_read[n][2];
		part[n].sex     = AB_data_read[n][3];
		part[n].height  = AB_data_read[n][4];
		part[n].weight  = AB_data_read[n][5];


		/////////////////////////////////////////////////
		// Fill out vaccine data

		part[n].prim_vac = AB_data_read[n][6];
		part[n].boost_vac = AB_data_read[n][7];

		part[n].t_prim_vac = AB_data_read[n][8];
		part[n].t_boost_vac = AB_data_read[n][9];


		/////////////////////////////////////////////////
		// Fill antibody data

		part[n].N_sam = 0;

		for (int j = 0; j<N_t; j++)
		{
			if (AB_data_read[n][24 + j] > -0.5)
			{
				part[n].tt.push_back(AB_data_read[n][10 + j]);
				part[n].AB.push_back(AB_data_read[n][24 + j]);

				part[n].lAB.push_back(log(AB_data_read[n][24 + j]));

				part[n].N_sam = part[n].N_sam + 1;
			}
		}


		/////////////////////////////////////////////////
		// Randomly assign individual-level parameters

		part[n].beta_prim  = genunf(100.0, 200.0);
		part[n].beta_boost = genunf(100.0, 200.0);
		part[n].t_short    = genunf(5.0, 10.0);
		part[n].t_long     = genunf(500.0, 1500.0);
		part[n].t_IgG      = genunf(10.0, 30.0);
		part[n].rho_prim   = genunf(0.8, 0.9);
		part[n].rho_boost  = genunf(0.8, 0.9);
		part[n].Ab_0       = genunf(0.0, 0.01);

		part[n].lbeta_prim  = log(part[n].beta_prim);
		part[n].lbeta_boost = log(part[n].beta_boost);

		part[n].lt_short = log(part[n].t_short);
		part[n].lt_long  = log(part[n].t_long);
		part[n].lt_IgG   = log(part[n].t_IgG);

		part[n].r_short = log2/part[n].t_short;
		part[n].r_long  = log2/part[n].t_long;
		part[n].r_IgG   = log2 / part[n].t_IgG;

		part[n].logitrho_prim  = log(part[n].rho_prim / (1.0 - part[n].rho_prim));
		part[n].logitrho_boost = log(part[n].rho_boost / (1.0 - part[n].rho_boost));

		part[n].lAb_0 = log(part[n].Ab_0);


		/////////////////////////////////////////////////
		// Calculate individual-level likelihood

		part[n].SSq       = SSq_n(&part[n]);

		part[n].data_like = data_like_n(&part[n], &theta);

		part[n].mix_like  = mix_like_n(&part[n], &theta);

		part[n].lhood     = part[n].data_like + part[n].mix_like;
	}

	AB_data_read.clear();


	//////////////////////////////////////////////////////
	// 1.8 Initialise adaptive MCMC object for individual-level parameters
	//     One object for each participant.

	part_n_MCMC* part_MCMC;
	part_MCMC = new part_n_MCMC[N_part];

	for (int n = 0; n<N_part; n++)
	{
		///////////////////////////////////
		// Parameter vector for MVN update

		part_MCMC[n].par_vec[0] = part[n].lbeta_prim;
		part_MCMC[n].par_vec[1] = part[n].lbeta_boost;
		part_MCMC[n].par_vec[2] = part[n].lt_short;
		part_MCMC[n].par_vec[3] = part[n].lt_long;
		part_MCMC[n].par_vec[4] = part[n].lt_IgG;
		part_MCMC[n].par_vec[5] = part[n].logitrho_prim;
		part_MCMC[n].par_vec[6] = part[n].logitrho_boost;
		part_MCMC[n].par_vec[7] = part[n].lAb_0;


		/////////////////////////////
		// Initialise diagonal covariance matrix

		for (int p = 0; p<N_loc_par; p++)
		{
			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].COV_MAT[p][q] = 0.0;
			}
		}

		part_MCMC[n].COV_MAT[0][0] = 0.2*0.2;
		part_MCMC[n].COV_MAT[1][1] = 0.2*0.2;
		part_MCMC[n].COV_MAT[2][2] = 0.2*0.2;
		part_MCMC[n].COV_MAT[3][3] = 0.2*0.2;
		part_MCMC[n].COV_MAT[4][4] = 0.2*0.2;
		part_MCMC[n].COV_MAT[5][5] = 0.2*0.2;
		part_MCMC[n].COV_MAT[6][6] = 0.2*0.2;
		part_MCMC[n].COV_MAT[7][7] = 0.2*0.2;


		/////////////////////////////
		// Counting moments

		for (int p = 0; p<N_loc_par; p++)
		{
			part_MCMC[n].par_S1[p] = part_MCMC[n].par_vec[p];

			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].par_S2[p][q] = part_MCMC[n].par_vec[p] * part_MCMC[n].par_vec[q];
			}
		}

		part_MCMC[n].denom = 1;


		/////////////////////////////
		// Set up dummy covariance matrix including
		// step-size scaling

		part_MCMC[n].step_scale = 0.1;

		for (int p = 0; p<N_loc_par; p++)
		{
			for (int q = 0; q<N_loc_par; q++)
			{
				part_MCMC[n].COV_MAT_dummy[p][q] = part_MCMC[n].step_scale*part_MCMC[n].COV_MAT[p][q];
			}
		}

		part_MCMC[n].accept = 0.0;
	}


	////////////////////////////////////////////////////////
	// 1.6 Book-keeping

	for (int p = 0; p < N_glob_par; p++)
	{
		theta.Y_par[p] = 0.0;
		theta.Ymu2_par[p] = 0.0;

		for (int n = 0; n < N_part; n++)
		{
			theta.Y_par[p]    = theta.Y_par[p] + part_MCMC[n].par_vec[p];
			theta.Ymu2_par[p] = theta.Ymu2_par[p] + (part_MCMC[n].par_vec[p] - theta.mu_par[p])*(part_MCMC[n].par_vec[p] - theta.mu_par[p]);
		}
	}

	theta.SSq_total = 0.0;
	for (int n = 0; n < N_part; n++)
	{
		theta.SSq_total = theta.SSq_total + part[n].SSq;
	}

	theta_p1 = theta;


	////////////////////////////////////////////////////////
	// 1.7 Create objects for updating local parameters

	part_n* part_p1;
	part_p1 = new part_n[N_part];

	for (int n = 0; n<N_part; n++)
	{
		part_p1[n] = part[n];
	}


	///////////////////////////////////////////////////////////////////////////
	// 1.9 Test output of likelihood

	for (int n = 0; n<N_part; n++)
	{
		cout << n << "\t" << "beta_prim: " << part[n].lbeta_prim << "\t" << "beta_boost: " << part[n].lbeta_boost << "\t" << "t_short: " << part[n].lt_short << "\t" << "t_long: " << part[n].lt_long << "\t" 
			              << "t_IgG: " << part[n].lt_IgG << "\t" << "rho_prim: " << part[n].logitrho_prim << "\t" << "rho_boost: " << part[n].logitrho_boost << "\t" << "Ab_0: " << part[n].lAb_0 << "\t" << "logL: " << part[n].lhood << endl;
	}


	///////////////////////////////////////////////////////////////////////////
	// 1.10 Initialise parameters for MCMC likelihood, Robbins-Munro 
	//      acceptance and output

	double loglike = global_prior(&theta);
	for (int n = 0; n<N_part; n++)
	{
		loglike = loglike + part[n].lhood;
	}

	double log_prob, loglike_p1;

	double log_loc_prob[N_part];

	int glob_out = max(2, (int)((int)N_mcmc) / 10000);
	int loc_out = max(2, (int)((int)N_mcmc) / 1000);


	vector<double> randomU(N_part);

	vector<double> loglike_vec_p1(N_part);


	///////////////////////////////////////////////////////////////////////////
	// 1.11 Open file for output and write first line

	std::cout << "START MCMC" << endl;

	cout << 0 << "\t";
	for (int p = 0; p < N_glob_par; p++)
	{
		cout << theta.mu_par[p] << "\t";
	}
	for (int p = 0; p < N_glob_par; p++)
	{
		cout << theta.tau_par[p] << "\t";
	}
	cout << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;


	std::ofstream global_MCMC_Stream(global_output_File);
	
	for (int p = 0; p < N_glob_par; p++)
	{
		global_MCMC_Stream << theta.mu_par[p] << "\t";
	}
	for (int p = 0; p < N_glob_par; p++)
	{
		global_MCMC_Stream << theta.tau_par[p] << "\t";
	}
	global_MCMC_Stream << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;


	std::ofstream local_MCMC_Stream(local_output_File);

	for (int n = 0; n<N_part; n++)
	{
		local_MCMC_Stream << part[n].beta_prim << "\t" << part[n].beta_boost << "\t" << part[n].t_short << "\t" << part[n].t_long << "\t" << part[n].t_IgG << "\t" << part[n].rho_prim << "\t" << part[n].rho_boost << "\t" << part[n].Ab_0 << "\t" << part[n].lhood << "\t";
	}
	local_MCMC_Stream << endl;


	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	//          //                                         //
	//   ####   //  Begin MCMC fitting procedure           //
	//  ##  ##  //                                         //
	//     ##   //                                         //
	//    ##    //                                         //
	//   #####  //                                         //
	//          //                                         //
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////

	for (int mc = 1; mc<N_mcmc; mc++)
	{
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////
		//       //                                    //
		//  2.1  //  UPDATE STAGE 1: INDIVIDUAL-LEVEL  //
		//       //  Metropolis-Hastings sampler       //
		//       //                                    //
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////

		/////////////////////////////////////////////
		// 2.1.1. Proposal step

		for (int n = 0; n < N_part; n++)
		{
			////////////////////////////////////////////////
			// Update COV_MAT_dummay

			for (int p = 0; p < N_loc_par; p++)
			{
				for (int q = 0; q < N_loc_par; q++)
				{
					part_MCMC[n].COV_MAT_dummy[p][q] = part_MCMC[n].step_scale*part_MCMC[n].COV_MAT[p][q];
				}
			}


			///////////////////////////////////////////////
			// Multi-variate Normal proposal step

			setgmn(part_MCMC[n].par_vec, *part_MCMC[n].COV_MAT_dummy, N_loc_par, part_MCMC[n].GMN_parm);

			genmn(part_MCMC[n].GMN_parm, part_MCMC[n].par_vec_test, part_MCMC[n].work);

			part_p1[n].lbeta_prim     = part_MCMC[n].par_vec_test[0];
			part_p1[n].lbeta_boost    = part_MCMC[n].par_vec_test[1];
			part_p1[n].lt_short       = part_MCMC[n].par_vec_test[2];
			part_p1[n].lt_long        = part_MCMC[n].par_vec_test[3];
			part_p1[n].lt_IgG         = part_MCMC[n].par_vec_test[4];
			part_p1[n].logitrho_prim  = part_MCMC[n].par_vec_test[5];
			part_p1[n].logitrho_boost = part_MCMC[n].par_vec_test[6];
			part_p1[n].lAb_0          = part_MCMC[n].par_vec_test[7];

			randomU[n] = genunf(0.0, 1.0);
		}


		/////////////////////////////////////////////
		// 2.1.2. Update step
		
		//#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < N_part; n++)
		{
			////////////////////////////////////////////////////////
			// 2.1.2.1. Only proceed if allowable parameters proposed

			if (local_prior(&part_p1[n]) > -0.5*LARGE)
			{
				part_p1[n].beta_prim = exp(part_p1[n].lbeta_prim);
				part_p1[n].beta_boost = exp(part_p1[n].lbeta_boost);

				part_p1[n].t_short = exp(part_p1[n].lt_short);
				part_p1[n].r_short = log2 / part_p1[n].t_short;

				part_p1[n].t_long = exp(part_p1[n].lt_long);
				part_p1[n].r_long = log2 / part_p1[n].t_long;

				part_p1[n].t_IgG = exp(part_p1[n].lt_IgG);
				part_p1[n].r_IgG = log2 / part_p1[n].t_IgG;

				part_p1[n].rho_prim = exp(part_p1[n].logitrho_prim) / (1.0 + exp(part_p1[n].logitrho_prim));
				part_p1[n].rho_boost = exp(part_p1[n].logitrho_boost) / (1.0 + exp(part_p1[n].logitrho_boost));

				part_p1[n].Ab_0 = exp(part_p1[n].lAb_0);

				part_p1[n].SSq = SSq_n(&part_p1[n]);
				part_p1[n].data_like = data_like_n(&part_p1[n], &theta);

				part_p1[n].mix_like = mix_like_n(&part_p1[n], &theta);

				part_p1[n].lhood = part_p1[n].data_like + part_p1[n].mix_like;


				double log_prob_n = part_p1[n].lhood - part[n].lhood;

				log_loc_prob[n] = _finite(log_prob_n) ? std::min(log_prob_n, 0.0) : -LARGE;


				////////////////////////////////////////
				// 2.1.2.2. Update if necessary

				if (log(randomU[n]) < log_loc_prob[n])
				{
					part[n] = part_p1[n];

					for (int p = 0; p < N_loc_par; p++)
					{
						part_MCMC[n].par_vec[p] = part_MCMC[n].par_vec_test[p];
					}

					part_MCMC[n].accept = part_MCMC[n].accept + 1;
				}


				////////////////////////////////////////
				// 2.1.3. Adjust step-size with Robbins-Monro
				//        Only do this for a local step within allowed range

				if (mc < N_adapt)
				{
					part_MCMC[n].step_scale = rm_scale(part_MCMC[n].step_scale, mc, N_adapt, log_loc_prob[n]);
				}
			}


			////////////////////////////////////////////////////////////
			// Running account of sums and sums of squares

			if (mc < N_tune_end)
			{
				for (int p = 0; p < N_loc_par; p++)
				{
					part_MCMC[n].par_S1[p] = part_MCMC[n].par_S1[p] + part_MCMC[n].par_vec[p];
					
					for (int q = 0; q < N_loc_par; q++)
					{
						part_MCMC[n].par_S2[p][q] = part_MCMC[n].par_S2[p][q] + part_MCMC[n].par_vec[p] * part_MCMC[n].par_vec[q];
					}
				}

				part_MCMC[n].denom = part_MCMC[n].denom + 1;
			}


			////////////////////////////////////////////////////////////
			// 2.1.4. TUNING STAGE 1

			/////////////////////////////////
			// Update covariance matrix

			if ((mc >= N_tune_start) && (mc < N_tune_end))
			{
				for (int p = 0; p < N_loc_par; p++)
				{
					for (int q = 0; q < N_loc_par; q++)
					{
						if (part_MCMC[n].accept / part_MCMC[n].denom > 0.01)
						{
							part_MCMC[n].COV_MAT[p][q] = part_MCMC[n].par_S2[p][q] / (part_MCMC[n].denom) - part_MCMC[n].par_S1[p] * part_MCMC[n].par_S1[q] / (part_MCMC[n].denom*part_MCMC[n].denom);
						}
					}
				}
			}

		}


		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);
		
		for (int n = 0; n < N_part; n++)
		{
			loglike = loglike + part[n].lhood;
		}



		///////////////////////////////////////////////////
		///////////////////////////////////////////////////
		//       //                                      //
		//  2.2  //  UPDATE STAGE 2 - POPULATION-LEVEL   //
		//       //  Gibbs sampler                       // 
		//       //                                      //
		///////////////////////////////////////////////////
		///////////////////////////////////////////////////

		for (int p = 0; p < N_glob_par; p++)
		{
			theta.Y_par[p] = 0.0;

			for (int n = 0; n < N_part; n++)
			{
				theta.Y_par[p] = theta.Y_par[p] + part_MCMC[n].par_vec[p];
			}

			theta.mu_par[p] = gennor( (theta.prior_tau[p] * theta.prior_mu[p] + theta.tau_par[p] * theta.Y_par[p]) / (theta.prior_tau[p] + N_part*theta.tau_par[p]),
					                  1.0 / sqrt(theta.prior_tau[p] + N_part*theta.tau_par[p]));

			theta.Ymu2_par[p] = 0.0;

			for (int n = 0; n < N_part; n++)
			{
				theta.Ymu2_par[p] = theta.Ymu2_par[p] + (part_MCMC[n].par_vec[p] - theta.mu_par[p])*(part_MCMC[n].par_vec[p] - theta.mu_par[p]);
			}


			theta.tau_par[p] = gengam( 1.0/theta.prior_theta[p] + 0.5*theta.Ymu2_par[p] ,
				                       0.5*N_part + theta.prior_k[p]);


		}

		theta_p1 = theta;


		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < N_part; n++)
		{
			//part[n].data_like = data_like_n(&part[n], &theta);

			part[n].mix_like = mix_like_n(&part[n], &theta);

			part[n].lhood = part[n].data_like + part[n].mix_like;

			loglike = loglike + part[n].lhood;
		}

		

		///////////////////////////////////////////////////
		///////////////////////////////////////////////////
		//       //                                      //
		//  2.3  //  UPDATE STAGE 3 - OBS-VARIANCE       //
		//       //  Metropolis-Hastings update          // 
		//       //                                      //
		///////////////////////////////////////////////////
		///////////////////////////////////////////////////

		///////////////////////////////////////////////
		// 2.3.1 Propose new parameters

		#pragma omp parallel for schedule(dynamic,4)
		for (int n = 0; n < N_part; n++)
		{
			part_p1[n] = part[n];
		}

		///////////////////////////////////////////////
		// Normal proposal step

		theta_p1.sig_obs = gennor(theta.sig_obs, 0.1*theta.sig_obs_scale);
	
		//////////////////////////////////////////////
		// 2.3.2 If parameters are acceptable, attempt update

		if (global_prior(&theta_p1) > -0.5*LARGE)
		{
			theta_p1.log_sig_obs = log(theta_p1.sig_obs);

			#pragma omp parallel for schedule(dynamic,4)
			for (int n = 0; n<N_part; n++)
			{
				part_p1[n].data_like = data_like_n(&part_p1[n], &theta_p1);
				
				part_p1[n].mix_like = mix_like_n(&part_p1[n], &theta_p1);

				part_p1[n].lhood = part_p1[n].data_like + part_p1[n].mix_like;

				loglike_vec_p1[n] = part_p1[n].lhood;
			}


			loglike_p1 = global_prior(&theta_p1);
			for (int n = 0; n<N_part; n++)
			{
				loglike_p1 = loglike_p1 + loglike_vec_p1[n];
			}

			const double log_prob0 = loglike_p1 - loglike;
			log_prob = _finite(log_prob0) ? std::min(log_prob0, 0.0) : -LARGE;

			if (log(genunf(0, 1)) < log_prob)
			{
				loglike = loglike_p1;
				theta = theta_p1;

				for (int n = 0; n<N_part; n++)
				{
					//part[n].lhood = loglike_vec_p1[n];

					part[n] = part_p1[n];
				}
			}


			////////////////////////////////////////
			// 2.3.4 Adjust step-size with Robbins-Monro

			if (mc < N_adapt)
			{
				theta.sig_obs_scale = rm_scale(theta.sig_obs_scale, mc + 1, N_adapt, log_prob);
			}

		}


		//////////////////////////////////////////////////////////////
		// 2.1.5. Update the total likelihood given the local updates

		loglike = global_prior(&theta);

		for (int n = 0; n < N_part; n++)
		{
			loglike = loglike + part[n].lhood;
		}


		/////////////////////////////////////////
		/////////////////////////////////////////
		//       //                            //
		//  2.3  //  Output results to file    //
		//       //                            //
		/////////////////////////////////////////
		/////////////////////////////////////////

		if (mc%glob_out == 0)
		{
			cout << 100 * ((double)mc) / ((double)N_mcmc) << "% complete." << "\t";
			for (int p = 0; p < N_glob_par; p++)
			{
				cout << theta.mu_par[p] << "\t";
			}
			for (int p = 0; p < N_glob_par; p++)
			{
				cout << theta.tau_par[p] << "\t";
			}
			cout << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;


			cout << "scalar    " << theta.sig_obs_scale << endl;


			for (int p = 0; p < N_glob_par; p++)
			{
				global_MCMC_Stream << theta.mu_par[p] << "\t";
			}
			for (int p = 0; p < N_glob_par; p++)
			{
				global_MCMC_Stream << theta.tau_par[p] << "\t";
			}
			global_MCMC_Stream << theta.sig_obs << "\t" << loglike << "\t" << global_prior(&theta) << endl;
		}



		if (mc%loc_out == 0)
		{
			for (int n = 0; n<N_part; n++)
			{
				local_MCMC_Stream << part[n].beta_prim << "\t" << part[n].beta_boost << "\t" << part[n].t_short << "\t" << part[n].t_long << "\t" << part[n].t_IgG << "\t" << part[n].rho_prim << "\t" << part[n].rho_boost << "\t" << part[n].Ab_0 << "\t" << part[n].lhood << "\t";
			}
			local_MCMC_Stream << endl;
		}

	}


	//////////////////////////////////////
	// 2.4. Calculate and output acceptance rates

	cout << "local acceptance:  " << part_MCMC[0].accept << "\t" << (double(part_MCMC[0].accept)) / (double(N_mcmc)) << endl;

	cout << "Time taken: " << cl << endl;

	return 0;
}



///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//          //                                       //
//   ####   //  Functions                            //
//  ##  ##  //                                       //
//     ##   //                                       //
//  ##  ##  //                                       //
//   ####   //                                       //
//          //                                       // 
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
// 3.1 Individual-level likelihood

double SSq_n(part_n* p)
{
	///////////////////////////////////////
	// Model-predicted antibody titre

	vector<double> model_AB(p->N_sam), model_lAB(p->N_sam);

	for (int j = 0; j < p->N_sam; j++)
	{
		model_AB[j] = p->Ab_0*exp(-p->r_long * (p->tt[j] - p->t_prim_vac));
	}

	if ( (p->prim_vac == 1) && (p->boost_vac==2) )
	{
		for (int j = 0; j < p->N_sam; j++)
		{
			model_AB[j] = model_AB[j] + p->beta_prim*(           (p->rho_prim / (p->r_short - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_prim_vac)) - exp(-p->r_short * (p->tt[j] - p->t_prim_vac)))
				                                       + ((1.0 - p->rho_prim) / (p->r_long - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_prim_vac)) - exp(-p->r_long * (p->tt[j] - p->t_prim_vac))));
		}
	}


	if ((p->prim_vac == 2) && (p->boost_vac == 1))
	{
		for (int j = 0; j < p->N_sam; j++)
		{
			if (p->tt[j] > p->t_boost_vac)
			{
				model_AB[j] = model_AB[j] + p->beta_prim*(           (p->rho_prim / (p->r_short - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_boost_vac)) - exp(-p->r_short * (p->tt[j] - p->t_boost_vac)))
				         	                               + ((1.0 - p->rho_prim) / (p->r_long - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_boost_vac)) - exp(-p->r_long * (p->tt[j] - p->t_boost_vac))));
			}
		}
	}


	if ((p->prim_vac == 1) && (p->boost_vac == 1))
	{
		for (int j = 0; j < p->N_sam; j++)
		{
			model_AB[j] = model_AB[j] + p->beta_prim*(           (p->rho_prim / (p->r_short - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_prim_vac)) - exp(-p->r_short * (p->tt[j] - p->t_prim_vac)))
				                                       + ((1.0 - p->rho_prim) / (p->r_long - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_prim_vac)) - exp(-p->r_long * (p->tt[j] - p->t_prim_vac)))    );

			if (p->tt[j] > p->t_boost_vac)
			{
				model_AB[j] = model_AB[j] + p->beta_boost*(           (p->rho_boost / (p->r_short - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_boost_vac)) - exp(-p->r_short * (p->tt[j] - p->t_boost_vac)))
					                                        + ((1.0 - p->rho_boost) / (p->r_long - p->r_IgG))*(exp(-p->r_IgG * (p->tt[j] - p->t_boost_vac)) - exp(-p->r_long * (p->tt[j] - p->t_boost_vac))));
			}
		}
	}


	for (int j = 0; j < p->N_sam; j++)
	{
		model_lAB[j] = log( model_AB[j] );
	}


	///////////////////////////////////////
	// Sum of squares

	double SSq = 0.0;

	for (int j = 0; j < p->N_sam; j++)
	{
		SSq = SSq + (p->lAB[j] - model_lAB[j])*(p->lAB[j] - model_lAB[j]);
	}

	return SSq;
}


double data_like_n(part_n* p, params* theta)
{
	double loglike = 0.0;

	for (int j = 0; j < p->N_sam; j++)
	{
		loglike = loglike - 0.9189385 - theta->log_sig_obs;// -p->lAB[j];
	}

	loglike = loglike - 0.5*p->SSq / (theta->sig_obs*theta->sig_obs);

	return loglike;

//	return ((double) p->N_sam)*( - 0.9189385 - theta->log_sig_obs ) - 0.5*p->SSq/(theta->sig_obs*theta->sig_obs);
}



double mix_like_n(part_n* p, params* theta)
{
	double mixlike = 0.0;

	///////////////////////////////////////
	// Mixed-effects component of likelihood

	mixlike = -6.43257 + 0.5*log(theta->tau_par[0]*theta->tau_par[1] *theta->tau_par[2] *theta->tau_par[3] * theta->tau_par[4] * theta->tau_par[5] * theta->tau_par[5])
		               - 0.5*theta->tau_par[0] *(p->lbeta_prim - theta->mu_par[0])*(p->lbeta_prim - theta->mu_par[0])
		               - 0.5*theta->tau_par[1] *(p->lbeta_boost - theta->mu_par[1])*(p->lbeta_boost - theta->mu_par[1])
		               - 0.5*theta->tau_par[2] *(p->lt_short - theta->mu_par[2])*(p->lt_short - theta->mu_par[2])
		               - 0.5*theta->tau_par[3] *(p->lt_long - theta->mu_par[3])*(p->lt_long - theta->mu_par[3])
		               - 0.5*theta->tau_par[4] *(p->lt_IgG - theta->mu_par[4])*(p->lt_IgG - theta->mu_par[4])
		               - 0.5*theta->tau_par[5] *(p->logitrho_prim - theta->mu_par[5])*(p->logitrho_prim - theta->mu_par[5])
	                   - 0.5*theta->tau_par[5] *(p->logitrho_boost - theta->mu_par[5])*(p->logitrho_boost - theta->mu_par[5]);

	return mixlike;
}



///////////////////////////////////////////////////////
// 3.2 Individual-level prior likelihood: excludes ineligible steps

double local_prior(part_n* p)
{
	if (p->beta_prim < 0.0) { return -LARGE; }
	if (p->beta_boost < 0.0) { return -LARGE; }

	if (p->lt_short > p->lt_long) { return -LARGE; }
	if (p->lt_short > p->lt_IgG) { return -LARGE; }


//	if (p->lt_IgG > 5.298317 ) { return -LARGE; }


	if (p->rho_prim < 0.0) { return -LARGE; }
	if (p->rho_prim > 1.0) { return -LARGE; }

	if (p->rho_boost < 0.0) { return -LARGE; }
	if (p->rho_boost > 1.0) { return -LARGE; }

	if (p->Ab_0 < 0.0) { return -LARGE; }

	return 0.0;
}


///////////////////////////////////////////////////////
// 3.3 Population-level prior likelihood: excludes ineligible steps

double global_prior(params* priors)
{
	double logprior = 0.0;

	for (int p = 0; p < N_glob_par; p++)
	{
		///////////////////////////
		// Normal prior on mean_par[p] 

		logprior = logprior - 0.9189385 + 0.5*log(priors->prior_tau[p]) - 0.5*priors->prior_tau[p] * (priors->mu_par[p] - priors->prior_mu[p])*(priors->mu_par[p] - priors->prior_mu[p]);


		///////////////////////////
		// Gamma prior on tau_par[p]

		logprior = logprior + (priors->prior_k[p] - 1.0)*log(priors->tau_par[p]) - priors->tau_par[p] / priors->prior_theta[p] - priors->prior_k[p] *log(priors->prior_theta[p]) - gammln(priors->prior_k[p]);
	}


	///////////////////////////
	// Gamma prior on tau_obs

	if (priors->sig_obs < 0.0) { return -LARGE; }
	if (priors->sig_obs > 2.0) { return -LARGE; }

	double mean_sig_obs = 2.0;
	double sd_sig_obs = 1.0;
	double k_sig_obs = (mean_sig_obs / sd_sig_obs)*(mean_sig_obs / sd_sig_obs);
	double theta_sig_obs = sd_sig_obs*sd_sig_obs / mean_sig_obs;

	double log_prior_sig_obs = (k_sig_obs - 1)*log(priors->sig_obs) - priors->sig_obs / theta_sig_obs - k_sig_obs*log(theta_sig_obs) - gammln(k_sig_obs);



//	double log_prior_sig_obs = (priors->prior_k_tau_obs - 1.0)*log(priors->tau_obs) - priors->tau_obs / priors->prior_theta_tau_obs - priors->prior_k_tau_obs*log(priors->prior_theta_tau_obs) - gammln(priors->prior_k_tau_obs);

	return logprior + log_prior_sig_obs;
}


/////////////////////////////////////////////////////////
// 3.4 Robbins-Monro algorithm for setting acceptance rate.
//
// Robbins-Munroe stochastic approximation adjustment
//	return adjustment
//	i iteration number
//	a scaling factor for adjustment size
//	m number of iterations when adjustment halves from value at i=0
//	d outcome at iteration i (e.g. calculated acceptance probability (min of MH ratio and 1),
//								or acceptance or no acceptance)
//	p desired mean for d (e.g. desired acceptance probability)

double rm_scale(double step_scale, int step, int N_step_adapt, double log_prob)
{
	double dd = exp(log_prob);
	if (dd < -30) { dd = 0.0; }
	dd = std::min(dd, 1.0);

	double rm_temp = (dd - 0.23) / ((double(step) + 1) / (0.01*(double(N_step_adapt)) + 1));

	double out = step_scale*exp(rm_temp);

	out = std::max(out, 0.02);
	out = std::min(out, 5.0);

	return out;
}


/////////////////////////////////////////////////////////
// 3.6 Log gamma function, based on gamma.h from NRC3

double gammln(const double xx) {
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = { 57.1562356658629235, -59.5979603554754912,
		14.1360979747417471, -0.491913816097620199, 0.339946499848118887e-4,
		0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3,
		-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3,
		0.844182239838527433e-4, -0.261908384015814087e-4, 0.368991826595316234e-5 };
	if (xx <= 0) throw("bad arg in gammln");
	y = x = xx;
	tmp = x + 5.242187500000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j<14; j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005*ser / x);
}

