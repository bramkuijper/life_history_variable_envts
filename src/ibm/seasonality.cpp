
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "bramauxiliary.h"

//#define NDEBUG
//#define DISTRIBUTION

using namespace std;

const int NumGen = 80000*12;
const size_t Npop = 2000;
const size_t skip = 5 * 12;
bool do_stats = 0;

size_t seed = 0;

// mutation rates and variances
double sdmu = 0;
double mu_a_tau = 0;
double mu_a_eps = 0;

// mortality prob reproductives
double m_r = 0;

// mortality prob nonreproductives
double m_n = 0;

// current timestep
int time_x = 0;

// current environment
double environment = 0;

// duration of a breeding bout
size_t b = 10;

// how total amount of envt translates to clutch
double clutch_multiplier = 1.0;

// minimum number of offspring produced by everybody
double min_clutch = 1.0;

// amplitude of sinusoidal process
double ampl = 1.0;

// degree of stochasticity
double stoch = 1.0;

// standard deviation of stochasticity
double sd_envt = 1.0;

// statistics

// tracking number of reproductives and nonreproductives
size_t Nrep = 0;
size_t Nnonrep = 0;

// tracking number of kids
size_t Noff = 0;



// gsl initialization
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual struct
struct Individual
{
    double a_tau[2]; // intercept for time sensitivity
    double b_tau[2]; // slope for time sensitivity
    double a_eps[2]; // intercept for sensitivity to environmental threshold
    double b_eps[2]; // intercept for sensitivity to environmental threshold

    size_t btime; // the amount of breeding time
    size_t ntime; // the amount of nonbreeding time
    double sum_envt; // cumulative amount of environment
};

typedef Individual Population[Npop];
typedef Individual Population2[Npop*500];
Population Rep, NonRep;
Population2 Offspring;

string filename("sim_seasonal");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// get some of the parameter values
// as command line arguments
void init_arguments(int argc, char *argv[])
{
    mu_a_tau = atof(argv[2]);
    mu_a_eps = atof(argv[3]);
    sdmu = atof(argv[4]);
    m_r = atof(argv[5]);
    m_n = atof(argv[6]);
    b = atoi(argv[7]);
    clutch_multiplier = atof(argv[8]);
    min_clutch = atof(argv[9]);
    ampl = atof(argv[10]);
    stoch = atof(argv[11]);
    sd_envt = atof(argv[12]);
}

void Mutate(Individual &ind)
{
    if (gsl_rng_uniform(r) < mu_a_tau)
    {
        ind.a_tau[0] += gsl_ran_gaussian(r,sdmu);
    }
    
    if (gsl_rng_uniform(r) < mu_a_tau)
    {
        ind.a_tau[1] += gsl_ran_gaussian(r,sdmu);
    }
    
    if (gsl_rng_uniform(r) < mu_a_tau)
    {
        ind.b_tau[0] += gsl_ran_gaussian(r,sdmu);
    }
    
    if (gsl_rng_uniform(r) < mu_a_tau)
    {
        ind.b_tau[1] += gsl_ran_gaussian(r,sdmu);
    }


    // mutate environmental sensitivity 
    if (gsl_rng_uniform(r) < mu_a_eps)
    {
        ind.a_eps[0] += gsl_ran_gaussian(r,sdmu);
    }

    if (gsl_rng_uniform(r) < mu_a_eps)
    {
        ind.a_eps[1] += gsl_ran_gaussian(r,sdmu);
    }
    
    if (gsl_rng_uniform(r) < mu_a_eps)
    {
        ind.b_eps[0] += gsl_ran_gaussian(r,sdmu);
    }

    if (gsl_rng_uniform(r) < mu_a_eps)
    {
        ind.b_eps[1] += gsl_ran_gaussian(r,sdmu);
    }
}


void write_parameters()
{
	DataFile << endl
		<< endl
        << "mu_a_eps;" << mu_a_eps << ";" << endl
        << "mu_a_tau;" << mu_a_tau << ";" << endl
        << "sdmu;" << sdmu << ";" << endl
        << "m_r;" << m_r << ";" << endl
        << "m_n;" << m_n << ";" << endl
        << "b;" << b << ";" << endl
        << "ampl;" << ampl << ";" << endl
        << "stoch;" << stoch << ";" << endl
        << "sd_envt;" << sd_envt << ";" << endl
        << "clutch_multiplier;" << clutch_multiplier << ";" << endl
        << "min_clutch;" << min_clutch << ";" << endl
		<< "seed;" << seed << ";"<< endl;
}


void write_data_headers()
{
    DataFile << "time;environment;mean_a_tau;var_a_tau;mean_b_tau;var_b_tau;mean_a_eps;var_a_eps;mean_b_eps;var_b_eps;mean_p;var_p;mean_ntime;Nrep;Nnonrep;Noff;" << endl;
}

void write_data()
{
    double mean_a_tau =0;
    double ss_a_tau = 0;
    double mean_a_eps =0;
    double ss_a_eps = 0;

    double mean_b_tau =0;
    double ss_b_tau = 0;
    double mean_b_eps =0;
    double ss_b_eps = 0;

    double mean_p =0;
    double ss_p = 0;

    double p = 0;
    double mean_ntime = 0;

    // get stats from the population
    for (size_t i =  0; i < Nnonrep; ++i)
    {
        mean_a_tau += NonRep[i].a_tau[0] + NonRep[i].a_tau[1];
        ss_a_tau += (NonRep[i].a_tau[0] + NonRep[i].a_tau[1]) * (NonRep[i].a_tau[0] + NonRep[i].a_tau[1]);
        
        mean_a_eps += NonRep[i].a_eps[0] + NonRep[i].a_eps[1];
        ss_a_eps += (NonRep[i].a_eps[0] + NonRep[i].a_eps[1]) * (NonRep[i].a_eps[0] + NonRep[i].a_eps[1]);
        
        mean_b_tau += NonRep[i].b_tau[0] + NonRep[i].b_tau[1];
        ss_b_tau += (NonRep[i].b_tau[0] + NonRep[i].b_tau[1]) * (NonRep[i].b_tau[0] + NonRep[i].b_tau[1]);
        
        mean_b_eps += NonRep[i].b_eps[0] + NonRep[i].b_eps[1];
        ss_b_eps += (NonRep[i].b_eps[0] + NonRep[i].b_eps[1]) * (NonRep[i].b_eps[0] + NonRep[i].b_eps[1]);

        
        p = 1.0 / (1.0 + exp(
                    - (NonRep[i].a_tau[0] + NonRep[i].a_tau[1] + (NonRep[i].b_tau[0] + NonRep[i].b_tau[1]) * NonRep[i].ntime) 
                    - (NonRep[i].a_eps[0] + NonRep[i].a_eps[1] + (NonRep[i].b_eps[0] + NonRep[i].b_eps[1]) * environment) 
                    )
                );

        mean_p += p;

        ss_p += p * p;

        mean_ntime += NonRep[i].ntime;
    }
    
    // get stats from the population
    for (size_t i =  0; i < Nrep; ++i)
    {
        mean_a_tau += Rep[i].a_tau[0] + Rep[i].a_tau[1];
        ss_a_tau += (Rep[i].a_tau[0] + Rep[i].a_tau[1]) * (Rep[i].a_tau[0] + Rep[i].a_tau[1]);
        
        mean_a_eps += Rep[i].a_eps[0] + Rep[i].a_eps[1];
        ss_a_eps += (Rep[i].a_eps[0] + Rep[i].a_eps[1]) * (Rep[i].a_eps[0] + Rep[i].a_eps[1]);
        
        mean_b_tau += Rep[i].b_tau[0] + Rep[i].b_tau[1];
        ss_b_tau += (Rep[i].b_tau[0] + Rep[i].b_tau[1]) * (Rep[i].b_tau[0] + Rep[i].b_tau[1]);
        
        mean_b_eps += Rep[i].b_eps[0] + Rep[i].b_eps[1];
        ss_b_eps += (Rep[i].b_eps[0] + Rep[i].b_eps[1]) * (Rep[i].b_eps[0] + Rep[i].b_eps[1]);
        
        p = 1.0 / (1.0 + exp(
                    - (Rep[i].a_tau[0] + Rep[i].a_tau[1] + (Rep[i].b_tau[0] + Rep[i].b_tau[1]) * Rep[i].ntime) 
                    - (Rep[i].a_eps[0] + Rep[i].a_eps[1] + (Rep[i].b_eps[0] + Rep[i].b_eps[1]) * environment) 
                    )
                );

        mean_p += p; 

        ss_p += p * p;
       
        // should work as ntime is not incremented whilst being a breeder
        mean_ntime += Rep[i].ntime;
    }

    int Ntot = Nrep + Nnonrep;

    DataFile << time_x << ";"
            << environment << ";"
            << (mean_a_tau/Ntot) << ";"
            << (ss_a_tau/Ntot - pow(mean_a_tau/Ntot,2.0)) << ";"
            << (mean_b_tau/Ntot) << ";"
            << (ss_b_tau/Ntot - pow(mean_b_tau/Ntot,2.0)) << ";"
            << (mean_a_eps/Ntot) << ";"
            << (ss_a_eps/Ntot - pow(mean_a_eps/Ntot,2.0)) << ";"
            << (mean_b_eps/Ntot) << ";"
            << (ss_b_eps/Ntot - pow(mean_b_eps/Ntot,2.0)) << ";"
            << (mean_p/Ntot) << ";"
            << (ss_p/Ntot - pow(mean_p/Ntot,2.0)) << ";"
            << (mean_ntime/Ntot) << ";"
            << Nrep << ";"
            << Nnonrep << ";"
            << Noff << ";"
            << endl;
}

// initialize the simulation
// by giving all the individuals 
// genotypic values
//
// and doing some other stuff (e.g., random seed)
void init_pop()
{
    // get the timestamp (with nanosecs)
    // to initialize the seed

	seed = get_nanoseconds();
    
    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);


	// initialize the population
	for (size_t i = 0; i < Npop/2; ++i)
	{
        for (size_t j = 0; j < 2; ++j)
        {
            NonRep[i].a_tau[j] = 0;
            NonRep[i].a_eps[j] = 0;
            NonRep[i].b_tau[j] = 0;
            NonRep[i].b_eps[j] = 0;

            Rep[i].a_tau[j] = 0;
            Rep[i].a_eps[j] = 0;
            Rep[i].b_tau[j] = 0;
            Rep[i].b_eps[j] = 0;
        }

        NonRep[i].ntime = 0;
        NonRep[i].btime = 0;
        Rep[i].ntime = 0;
        Rep[i].btime = 0;
	}

    Nnonrep = Npop/2;
    Nrep = Npop/2;

    environment = 0.0;
}

// create an offspring
void create_kid(Individual &mother, Individual &father, Individual &kid)
{
    kid.a_tau[0] = mother.a_tau[gsl_rng_uniform_int(r, 2)];
    kid.a_tau[1] = father.a_tau[gsl_rng_uniform_int(r, 2)];
    
    kid.b_tau[0] = mother.b_tau[gsl_rng_uniform_int(r, 2)];
    kid.b_tau[1] = father.b_tau[gsl_rng_uniform_int(r, 2)];

    kid.a_eps[0] = mother.a_eps[gsl_rng_uniform_int(r, 2)];
    kid.a_eps[1] = father.a_eps[gsl_rng_uniform_int(r, 2)];
    
    kid.b_eps[0] = mother.b_eps[gsl_rng_uniform_int(r, 2)];
    kid.b_eps[1] = father.b_eps[gsl_rng_uniform_int(r, 2)];
    
    kid.ntime = 0;
    kid.btime = 0;

    Mutate(kid);
}

// Survival of juveniles to reproductive adults
void reproduce()
{
    if (Nnonrep && Nrep == 0)
    {
        write_parameters();
        exit(1);
    }

    int clutch;

    Noff = 0;
    double p;

    // recruit reproductives
    for (int nonrep_i = 0; nonrep_i < Nnonrep; ++nonrep_i)
    {
        p = 1.0/(1.0 + 
                exp( -(NonRep[nonrep_i].a_tau[0] + NonRep[nonrep_i].a_tau[1] + (NonRep[nonrep_i].b_tau[0] + NonRep[nonrep_i].b_tau[1]) * NonRep[nonrep_i].ntime) 
                    - (NonRep[nonrep_i].a_eps[0] + NonRep[nonrep_i].a_eps[1] + (NonRep[nonrep_i].b_eps[0] + NonRep[nonrep_i].b_eps[1]) * environment)
                    )
                );

        // nonreproductive becoming reproductive
        if (gsl_rng_uniform(r) < p)
        {
            // set breeding time to 0
            NonRep[nonrep_i].btime = 0;
            NonRep[nonrep_i].ntime = 0;
            NonRep[nonrep_i].sum_envt = 0;

            Rep[Nrep++] = NonRep[nonrep_i];
            NonRep[nonrep_i] = NonRep[Nnonrep-1];
            --nonrep_i;
            --Nnonrep;
        }
        else
        {
            ++NonRep[nonrep_i].ntime;
        }
    }

    size_t father = 0;

    // now let the nonreproductives die or reproduce
    for (int rep_i = 0; rep_i < Nrep; ++rep_i)
    {
        // increment time counter
        ++Rep[rep_i].btime;
        ++Rep[rep_i].ntime;

        // increment environmental values (additive model)
        Rep[rep_i].sum_envt += environment;

        // end of breeding bout
        if (Rep[rep_i].btime >= b)
        {
            clutch = floor(clutch_multiplier * Rep[rep_i].sum_envt);

            if (gsl_rng_uniform(r) < (clutch_multiplier * Rep[rep_i].sum_envt - clutch))
            {
                ++clutch;
            }

            clutch = min_clutch + (clutch < 0 ? 0 : clutch);

            if (Nrep > 1)
            {
                // create offspring
                for (int off_i = 0; off_i < clutch; ++off_i)
                {
                    // pick a father from the other nonbreeders
                    do {
                            father = gsl_rng_uniform_int(r, Nrep);
                    }
                    while(father == rep_i);

                    Individual Kid;
                    create_kid(Rep[rep_i], Rep[father], Kid);
                    
                    // add kid to offspring stack
                    Offspring[Noff++] = Kid;
                }
            }
            
            // kids made. Now let's remove this individual from the breeding population
            NonRep[Nnonrep] = Rep[rep_i];
            ++Nnonrep;

            Rep[rep_i] = Rep[Nrep-1];
            --rep_i;
            --Nrep;

            assert(Nrep >= 0);
            assert(Nrep < Npop);
            assert(rep_i >= -1);
            assert(Nnonrep >= 0);
            assert(Nnonrep < Npop);
        }
    }
}

void replace_adults()
{
    // mortality of nonreproductives
    for (size_t nonrep_i = 0; nonrep_i < Nnonrep; ++nonrep_i)
    {
        if (gsl_rng_uniform(r) < m_n)
        {
            NonRep[nonrep_i] = NonRep[Nnonrep-1];
            --nonrep_i;
            --Nnonrep;
        }
    }

    // mortality of reproductives
    for (size_t rep_i = 0; rep_i < Nrep; ++rep_i)
    {
        if (gsl_rng_uniform(r) < m_r)
        {
            Rep[rep_i] = Rep[Nrep-1];
            --rep_i;
            --Nrep;
        }
    }

    size_t deficit = Npop - Nnonrep - Nrep;

    for (size_t i = 0; i < deficit; ++i)
    {
        if (Noff == 0)
        {
            break;
        }

        NonRep[Nnonrep++] = Offspring[gsl_rng_uniform_int(r, Noff)];
        --Noff;
    }
    
    // update the environment
    environment = ampl * sin(2.0 * M_PI / 12 * time_x) + stoch * gsl_ran_gaussian(r, sd_envt);
}

int main(int argc, char ** argv)
{
	init_arguments(argc, argv);
	write_data_headers();
	init_pop();

	for (time_x = 0; time_x <= NumGen; ++time_x)
	{
        do_stats = time_x % skip == 0 || time_x > NumGen - 250;

		reproduce();

        replace_adults();

        if (do_stats)
		{
			write_data();
		}
	}

	write_parameters();
}
