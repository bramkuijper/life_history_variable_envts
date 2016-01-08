
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

const int NumGen = 100000;
const int Npop = 5000;
const int skip = 10;
int Noffspring = 0;
int N = 0;
int Nsurv = 0;
int total_offspring = 0;
bool do_stats = 0;

double epsilon = 0; // the value of the environment
double epsilon_sens = 0; // the perceived value of the environment

double size_min = 0.1;
double init_size = 0.1;

int offspring_control = 0;

double theta = 0; 
double omega_b_2 = 0; // width of the selection function
double omega_m_2 = 0; // width of the selection function
double omega_size2 = 0; // width of the selection function
double wmin = 1.0; // width of the selection function
double	A = 0; // baseline of the optimum size
double	B = 0; // amplitude of fluctuations
double sigma_e = 1.0; // developmental noise
double sigma_ksi = 0.1; // variance of the autocorrelated process
double rho_t = 0.5; // temporal autocorrelation
double freq = 0; // temporal autocorrelation
double mu_g 	  = 0.05;            // mutation rate
double sdmu_g         = 0.4;			 // standard deviation mutation size
double mu_m 	  = 0.05;            // mutation rate
double sdmu_m         = 0.4;			 // standard deviation mutation size
double mu_b 	  = 0.05;            // mutation rate
double sdmu_b         = 0.4;			 // standard deviation mutation size
double ksi = 0;			 // standard deviation mutation size
double ampl = 2.0;
double tau = 0.0;
double mean_surv = 0.0;
double total_surv = 0.0;

const int nloci_b = 2; // number of alleles underlying genetic architecture
const int nloci_g = 2; // number of alleles underlying genetic architecture
const int nloci_m = 2; // number of alleles underlying genetic architecture
const int nloci_s = 2; // number of alleles underlying genetic architecture

int generation = 0;
unsigned seed = 0;



// gsl initialization
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual struct
struct Individual
{
    // assuming haploid inheritance 
    // and 50 loci coding for g to get proper values of 
    // the genetic variance covariance matrix
    double g[nloci_g]; 
    double b[nloci_b]; 
    double m[nloci_m];
    double s[nloci_s];
    double phen;
    double phen_m;
    double phen_b;
    double phen_g;
    double phen_s;
    double cumsurv;
    double surv;
    int generation;
};

typedef Individual Population[Npop];
typedef Individual Population2[Npop*500];
Population Pop;
Population2 Offspring;

string filename("sim_evolving_m");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

#ifdef DISTRIBUTION
string filename_new2(create_filename("sim_evolving_m_dist"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION

// get some of the parameter values
// as command line arguments
void initArguments(int argc, char *argv[])
{
	mu_g = atof(argv[1]);
	mu_m = atof(argv[2]);
	mu_b = atof(argv[3]);
	mu_s = atof(argv[3]);
	sdmu_g = atof(argv[4]);
	sdmu_m = atof(argv[5]);
	sdmu_b = atof(argv[6]);
	sdmu_s = atof(argv[6]);
	freq = atof(argv[7]);
	sigma_e = sqrt(atof(argv[8]));
	sigma_ksi = sqrt(atof(argv[9]));
	wmin = atof(argv[10]);
	rho_t = atof(argv[11]);
	omega_b_2 = atof(argv[12]);
	omega_m_2 = atof(argv[13]);
	tau = atof(argv[14]);
	A = atof(argv[15]);
	B = atof(argv[16]);
	omega_size2 = atof(argv[17]);
    init_size = atof(argv[18]);
    offspring_control = atoi(argv[19]);
}

void MutateG(double &G)
{
	G += gsl_rng_uniform(r)<mu_g ? gsl_ran_gaussian(r, sdmu_g) : 0;
}

void MutateM(double &G)
{
	G += gsl_rng_uniform(r)<mu_m ? gsl_ran_gaussian(r,sdmu_m) : 0;
}

void MutateB(double &G)
{
	G += gsl_rng_uniform(r)<mu_b ? gsl_ran_gaussian(r,sdmu_b) : 0;
}

// standard dev only evolves 
void MutateS(double &G)
{
    if (gsl_rng_uniform(r) < mu_s)
    {
        G += gsl_ran_gaussian(r,sdmu_s);

        if (G < 0)
        {
            G = 0;
        }
    }
}

void WriteParameters()
{

	DataFile << endl
		<< endl
		<< "type:;" << "evol_m_size" << ";" << endl
		<< "control:;" << offspring_control << ";" << endl
        << "size_min:;" << size_min << ";" << endl
        << "mu_g:;" << mu_g << ";" << endl
        << "mu_m:;" << mu_m << ";" << endl
        << "mu_b:;" << mu_b << ";" << endl
        << "mu_s:;" << mu_s << ";" << endl
        << "sdmu_g:;" << sdmu_g << ";" << endl
        << "sdmu_m:;" << sdmu_m << ";" << endl
        << "sdmu_b:;" << sdmu_b << ";" << endl
        << "sdmu_s:;" << sdmu_s << ";" << endl
        << "wmin:;" << wmin << ";" << endl
        << "freq:;" << freq << ";" << endl
        << "sigma_e:;" << sigma_e << ";" << endl
        << "sigma_ksi:;" << sigma_ksi << ";" << endl
        << "rho_t:;" << rho_t << ";" << endl
        << "tau:;" << tau << ";" << endl
        << "A:;" << A << ";" << endl
        << "B:;" << B << ";" << endl
		<< "seed:;" << seed << ";"<< endl;
}


void write_data()
{
    double meanphen = 0;
    double meanphen_m = 0;
    double meang = 0;
    double ssg = 0;
    double meanm = 0;
    double ssm = 0;
    double meanb = 0;
    double ssb = 0;
    double means = 0;
    double sss = 0;
    double ssphen = 0;

    // get stats from the population
    for (int i =  0; i < N; ++i)
    {
        // stats for m
        meang += Pop[i].phen_g;
        ssg += Pop[i].phen_g * Pop[i].phen_g;

        // stats for m
        meanm += Pop[i].phen_m;
        ssm += Pop[i].phen_m * Pop[i].phen_m;

        meanb += Pop[i].phen_b;
        ssb += Pop[i].phen_b * Pop[i].phen_b;
        
        means += Pop[i].phen_s;
        sss += Pop[i].phen_s * Pop[i].phen_s;

        meanphen += Pop[i].phen;
        meanphen_m += Pop[i].phen_m;

        ssphen += Pop[i].phen * Pop[i].phen;
    }

    DataFile << generation << ";" << epsilon << ";" << total_offspring << ";" << Nsurv << ";" << mean_surv << ";" << ksi << ";";

    DataFile 
            << (meanphen/N) << ";"
            << (ssphen/N - pow(meanphen/N,2.0)) << ";"
            << (meanphen_m/N) << ";"
            << (meang/(N)) << ";"
            << (ssg/(N) - pow(meang/N,2.0)) << ";"
            << (means/(N)) << ";"
            << (sss/(N) - pow(means/N,2.0)) << ";"
            << (meanm/(N)) << ";"
            << (ssm/N - pow(meanm/N,2.0)) << ";" 
            << (meanb/N) << ";"
            << (ssb/N - pow(meanb/N,2.0)) << ";"  << endl;
}

// initialize the simulation
// by giving all the individuals 
// genotypic values
//
// and doing some other stuff (e.g., random seed)
void Init()
{
    // get the timestamp (with nanosecs)
    // to initialize the seed

	seed = get_nanoseconds();
    
    // set the seed to the random number generator
    // stupidly enough, for gsl this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));


    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);


	// initialize the whole populatin
	for (int i = 0; i < Npop; ++i)
	{
        Pop[i].phen = init_size;
        Pop[i].phen_m = init_size;

        for (int j = 0; j < nloci_g; ++j)
        {
            Pop[i].g[j] = .5*init_size;
        }

        for (int j = 0; j < nloci_b; ++j)
        {
            Pop[i].b[j] = 0;
        }

        for (int j = 0; j < nloci_m; ++j)
        {
            Pop[i].m[j] = 0;
        }
        
        for (int j = 0; j < nloci_s; ++j)
        {
            Pop[i].m[j] = 0;
        }

        Pop[i].generation = generation - 1;
	}

    // burn in, first let envt adapt
    epsilon = 1;
    N = Npop;
}

// create an offspring
void Create_Kid(int mother, int father, Individual &kid, double resources)
{
    double sum_g = 0; // sum over all the breeding values of the offspring coding for the actual phenotype
    double sum_b = 0; // sum over all the breeding values of the offspring coding for the norm of reaction
    double sum_m = 0; // sum over all the breeding values of the offspring coding for the maternal effect

    assert(mother >= 0 && mother <= N);
    assert(father >= 0 && father <= N);

    // we assume all loci are unlinked (i.e., parent of origin (poi)
    // of one allele is independent of poi of any other alles
    for (int i = 0; i < nloci_g;++i)
    {
        kid.g[i] = i % 2 == 0 ? Pop[mother].g[i + gsl_rng_uniform_int(r, 2)] : Pop[father].g[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateG(kid.g[i]);
        sum_g += kid.g[i];
    }

    for (int i = 0; i < nloci_b; ++i)
    {
        kid.b[i] = i % 2 == 0 ? Pop[mother].b[i + gsl_rng_uniform_int(r, 2)] : Pop[father].b[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateB(kid.b[i]);
        sum_b += kid.b[i];
    }

    for (int i = 0; i < nloci_m; ++i)
    {
        kid.m[i] = i % 2 == 0 ? Pop[mother].m[i + gsl_rng_uniform_int(r, 2)] : Pop[father].m[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateM(kid.m[i]);
        sum_m += kid.m[i];
    }
    
    for (int i = 0; i < nloci_s; ++i)
    {
        kid.s[i] = i % 2 == 0 ? Pop[mother].s[i + gsl_rng_uniform_int(r, 2)] : Pop[father].s[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateS(kid.s[i]);
        sum_s += kid.s[i];
    }

    kid.phen_m = sum_m;
    kid.phen_g = sum_g;
    kid.phen_b = sum_b;
    kid.phen_s = sum_s;
    kid.generation = generation;
    // kid's phenotype at time t is then:
    // random noise +
    // plasticity (b * b) 
    // breeding values + 
    // m(t-1) * z(t-1)
    //
    // check for any leaks in the population
    assert(Pop[mother].generation == kid.generation - 1);
    assert(Pop[father].generation == kid.generation - 1);

    // z = a + b * epsilon + m * R + x
    kid.phen = kid.phen_g  + kid.phen_b * epsilon_sens + kid.phen_m * Pop[mother].resources + gsl_ran_gaussian(r,s);

    assert(isnan(kid.phen) == 0);
}


// Survival of juveniles to reproductive adults
void reproduce()
{
    // reset the number of offspring produced
    Noffspring = 0;

    // count the total number of survivors
    total_surv = 0;

    mean_surv = 0;

    double mean_size = 0;

    size_t clutch_size;

    // calculate envtal variables as functions of epsilon
    resources = intercept_r + amplitude_r * epsilon;
    mmint = intercept_mmin + amplitude_mmin * epsilon * rho_mmin;

    // environmental predictability
    epsilon_sens = epsilon + gsl_ran_gaussian(noise_detect);

    // allow all individuals to reproduce
    // as they are hermaphrodites
    for (size_t i = 0; i < N; ++i)
    {
        size_t father;

        // choose a father (no selfing)
        do {
            father = gsl_rng_uniform_int(r, N);
        } while (father == i);

        clutch_size = floor(resources);

        // stochastic rounding of continous clutch size
        if (gsl_rng_uniform(r) < resources - clutch_size)
        {
            ++clutch_size;
        }

        // create offspring
        for (size_t kid_i = 0; kid_i < clutch_size; ++kid_i)
        {
            // create individual offspring
            Individual Kid;
            Create_Kid(i, father, Kid, resources);

            // reduce maternal resources
            if (Kid.phen <= clutch_size)
            {
                clutch_size -= Kid.phen;

                // size-dependent juv survival
                if (gsl_rng_uniform(r) < 1.0 - exp(-ct * (Kid.phen - mmint)))
                {
                    Offspring[Noffspring++] = Kid;
                    mean_size += Kid.phen;
                    mean_surv += Kid.surv;
                }
            } 
            else if (gsl_rng_uniform(r) < clutch_size/Kid.phen)
            {
                if (gsl_rng_uniform(r) < 1.0 - exp(-ct * (Kid.phen - mmint)))
                {
                    Offspring[Noffspring++] = Kid;
                    mean_size += Kid.phen;
                    mean_surv += Kid.surv;
                }
                break;
            }
            else
            {
                break;
            }
        }
    }


    mean_surv /= Noffspring;
    mean_size /= Noffspring;

    if (Noffspring <= 1 || total_surv == 0)
    {
        cout << (total_surv == 0 ? "no survivors. " : "no offspring produced. ") << endl;
        write_data();
        exit(1);
    }
}

void replace_adults()
{
    // sample from newborns with replacement
    for (size_t i = 0; i < N; ++i)
    {
        Pop[i] = Offspring[gsl_rng_uniform_int(r, Noffspring)];
    }
    
    // update the environment
    xt = nu * xt + gsl_ran_gaussian(r, sd_stoch);
    epsilon = sin(freq * (generation)) + xt;
}

void write_data_headers()
{
    DataFile << "generation;epsilon;noffspring;nsurv;mean_surv;ksi;meanz;varz;meanphen_m;meang;varg;means;vars;meanm;varm;meanb;varb;" << endl;
}

int main(int argc, char ** argv)
{
	init_arguments(argc, argv);
	write_data_headers();
	init_pop();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		reproduce();
        replace();

        if (do_stats)
		{
			write_data();
		}
	}

	write_parameters();
}
