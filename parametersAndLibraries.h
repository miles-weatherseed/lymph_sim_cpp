#ifndef __PARAMETERS_INCLUDED__
#define __PARAMETERS_INCLUDED__
#define DEFAULT_RADIUS (500)
//in micrometres
#define DEFAULT_CONTACT_RADIUS (20)
#define T_CELL_DENSITY (1E6)
//per mm^3.
#define TOTAL_TNUM (T_CELL_DENSITY*4.0*0.3333*3.1415926 * radius*radius*radius * 1E-9)
//2.6E7
#define DEFAULT_AGSPEC_FREQ (1E-5)
#define DEFAULT_DENDNUM (720)
#define DEFAULT_SEED (6)
#define DEFAULT_NUM_REPEATS (10)
#define DEFAULT_NUM_TIME_MEASUREMENTS (300)
#define DEFAULT_T_VELOCITY_MEAN (10.06354)
#define DEFAULT_T_VELOCITY_STDEV (0.5) //Only for gauss
#define T_GAMMA_SHAPE (2.28567314561) //Only for gamma
#define T_GAMMA_SCALE(v) (v/DEFAULT_T_VELOCITY_MEAN*4.40287799651)
#define DEFAULT_T_FREE_PATH_MEAN (25)
#define DEFAULT_T_FREE_PATH_STDEV (3)
#define DEFAULT_NUM_ANTIGEN_IN_CONTACT_AREA (500)
#define DEFAULT_TCELL_ACTIVATION_THRESHOLD (20)
//Density around 2.5 complexes / mumetre^2; so 8*2.5 ~ 20; 10 from "T cell activation is determined by the number of presentated ags" [2013] & ~10 from "functional anatomy of T cell activation and synapse formation" [2010]
#define DEFAULT_FIRST_DC_ARRIVAL (1080)
#define DEFAULT_DC_ARRIVAL_DURATION (360)
#define DEFAULT_ANTIGEN_DECAY_RATE (0.00135)
//Per minute, from k2=ln(2)/650mins + MHCturnover=ln(2)/36hrs. w/ MHC turnover of 7hours:  0.002716730707689163
//Set by MHC turnover.
#define DEFAULT_COGNATE_RATIO_MEAN_AT_DERMIS (0.2)
#define PROBABILITY_TOLERANCE (1e-8)
//Set probabilities below this to 0
#define TIME_TOLERANCE (1e-6)
#define POSITION_OOB_TOLERANCE (100*numDCells)
#define MINIMUM_CONTACT_GRID_MEMORY_REDUCTION_FACTOR (0.25)

//Stationary DC version
#define DEFAULT_TIMESTEP (0.05)
#define NUM_TIMESTEPS_HOUR(ts) (60 / ts)
#define NUM_TIMESTEPS_DAY(ts) (1440 / ts)
#define NUM_TIMESTEPS_2DAY(ts) (2880 / ts)

using namespace std;
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
// #include <boost/math/distributions/pareto.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/format.hpp>
#include <random> //needs c++11
#include <sys/stat.h> //For checking existence of file
#include <ctime>

#if defined(__linux__)
#include <unistd.h>
#elif defined(__WIN32) || defined(__CYGWIN__)
#include <windows.h>
		typedef unsigned int uint;
	#else
		#error "Unknown platform"
		return 1
#endif

#endif //__PARAMETERS_INCLUDED__