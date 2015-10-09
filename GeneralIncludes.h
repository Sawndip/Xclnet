#ifndef GENERAL_INCLUDES_H_
#define GENERAL_INCLUDES_H_

//#define DEBUG_MODE // not currently used
//#define DEBUG_MODE_MAIN // screen display of timestep based change of V and RHO
//#define DEBUG_MODE_NETWORK // print to file of connectivity statistics
//#define DEBUG_MODE_SPIKES // screen display of spike transfers
//#define DEBUG_MODE_SYNAPSE // screen display of synapse updates

#define ENABLE_SYNAPSE_UPDATES
//#define ENABLE_FIXED_TRANSFERS
//#define ENABLE_TRANSFER_RHO_INITIAL

//#define SYN_USE_SUPRATHRESHOLD_TIMESTEP
#define SYN_USE_FLAT_POTENTIAL
#define SYN_USE_HARD_BOUNDS

#define SYN_USE_CONST_INITIALISATION
//#define SYN_USE_RAND_UNIFORM_INITIALISATION
//#define SYN_USE_INVIVO_DOUBLE_WELL_INITIALISATION
//#define SYN_POTENTIATE_SUBSET_OF_SYNS
//#define SYN_DELAYED_POTENTIATE_SUBSET_OF_SYNS

#define MONITOR_UP_DOWN_POPS
//#define MONITOR_STIM_POPS

//#define LEARNING_CONST_EXTERN_STIM
//#define LEARNING_INDEP_POISSON_STIM
#define LEARNING_REPEATED_PATTERNED_STIM

#define FILE_NAME_LENGTH (50)
#define TEXT_BUFFER_LENGTH (100)
#define EPSILLON (0.0000001)

// Data reporters
#define RECORDER_NEURON_ID (48) /*(3)*/
#define RECORDER_SYNAPSE_ID (0) /* for modulo addressed multiple synaptic recordings this needs to be less than 400*/
//#define RECORDER_MULTI_SYNAPSE_SKIP (64000) /*(64000)*/ /*(450)*/

#define USE_GPU (1) /* 1=gpu, 0=cpu */
#define NETWORK_SEED (-14) /*(-14)*/
#define PARALLEL_SEED (2) /*keep positive for random123*/
#define GAUSSIAN_SYNAPTIC_SEED (-12)
#define UNIFORM_SYNAPTIC_SEED (-11)

#define RAN2_RESETTABLE_SEED (-1)

#define MAX_TIME_STEPS (5000) /*(30500000)*/ /*(50000000)*/ /*(12000000)*/ /*(12000000)*/ /*(300000)*/ /*no of timesteps, each of size dt*/

// Network schema
#define NO_EXC (100) /*(8000)*/ /*(400)*/ /*(10000)*/
#define NO_INH (0) /*(2000)*/
#define NO_LIFS (NO_EXC + NO_INH) 
#define CONNECTIVITY_PROBABILITY (1.1) /*(0.05)*/ /*(0.05)*/

// Time step sizes and statistical bin widths
#define LIF_DT (0.00001) /* modify refrac time and calcium delay in tandem, also MAX_TIME_STEPS */
#define SYN_DT LIF_DT /*TODO: at a later stage I will have the synapse update more slowly than the lif*/
#define BIN_SIZE (1.) /*(1.)*/ /*(0.1)*/


// Stimulation of subpopulation /* using secs despite inconsistency with other parameter units */
#define STIM_ON (1.)
#define STIM_OFF (500.) /* want stimulation for full duration of simulation */
#define J_STIM (50) /* a delta stim should cause a spike */ /*(24.55)*/ /*24.55mV approx 50Hz, 34.8mV approx 100Hz*/
#define NO_STIM_LIFS (100)
#define STIM_OFFSET (0)
#define NO_STIM_SUBSETS (1)

#define STIM_PATTERN_DURATION (50000) /* in timesteps (otherwise should use it to initialise a variable)*/
//#define STIM_PATTERN_DURATION_SECS (0.1) /* pattern duration in secs, for regular patterned stim */
#define STIM_PATTERN_AV_RATE (1.) /* for exponential drawing of interspike interval */
#define STIM_PATTERN_PAUSE_DURATION (50000) /* in timesteps: pause between pattern repeats */
#define STIM_FIXED_OFFSET_ISI (500) /* in timesteps: pause between stimuli on adjacent neurons */

// Transfer voltages
#define J_EE (2.2) /*(22.0)*/ /*(0.2)*/ /*(0.1)*/
#define J_IE (0.1)
#define J_II (-0.4)
#define J_EI (-0.4)

//#define J_EXT (11.45) /*(15.315)*/ /*(15.315)*/ /*(11.046)*/ /*(11.45) new in-vivo*/ /*(11.046) new in-vitro*/ /*(7.07)*/ /*(6.966) in-vivo*/ /*(7.07) 1hz in-vitro*/
#define J_EXT (0.0)


// LIF specific
#define LIF_V_INITIAL (-60.0)
#define LIF_V_REST (-70.0)
#define LIF_V_RESET (-70.0) /*(-60.0)*/
#define LIF_V_THRESHOLD (-50.0)
#define LIF_CM (0.001)
#define LIF_RM (20.0)
#define LIF_SIGMA (0) /*(5)*/
#define LIF_REFRAC_TIME (2000) /*200*/ /*timesteps*/

// Synapse model specific
#define SYN_RHO_INITIAL (0.0) /*arbitrary initialisation for now */ /*(0.019)*/ /*(0.406595)*/ /*(0.164840)*/ /*(0.406595)*/ /*(0.019) in-vivo*/ /*(0.164840) new in-vitro*/ /*(0.16492)*/ /*(0.203586)*/ /*(1.0)*/
#define SYN_CA_INITIAL (0.0)
#define SYN_CALCIUM_DELAY (461) /*46*/ /*4.6098ms*/ /*timesteps (needs to be modified when DT is modified above*/
#define SYN_GAMMA_P (725.085)
#define SYN_GAMMA_D (331.909)
#define SYN_THETA_P (1.3)
#define SYN_THETA_D (1.0)
#define SYN_SIGMA (0) /*(3.35)*/ /*3.35;*/ /*TODO: switch synapse noise back on*/
#define SYN_TAU (346.3615)
#define SYN_TAU_CA (0.0226936)
#define SYN_C_PRE (0.56175) /*(0.33705)*/ /*(0.56175)*/ /*(0.5617539)*/
#define SYN_C_POST (1.23964) /*(0.74378)*/ /*(1.23964)*/

#define SYN_RHO_FIXED SYN_RHO_INITIAL /*(0.5)*/



#include <stdio.h>
#include <errno.h> // Error numbers for file opening
//#include <limits.h> // System data limits
//#include <float.h> // Limits of system implementation of floating values
#include <stdlib.h> // malloc()
#include <fcntl.h> // manipulate file descriptors, opencl
//#include <strcmp.h>
#include <math.h> // supposedly for fmin()
#include <time.h> // to time main loop
#include <string.h> // strcpy() strcat()
//#include <ctype.h>
//#include <assert.h>
//#include <locale.h>
//#include <stddef.h>
#include <sys/types.h> // opencl
#include <sys/stat.h> // For mkdir()
//#include <sys/dir.h>

//from Mac OpenCL demo
#include <unistd.h>
#ifdef __APPLE__ //ifdef added from web suggestions
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#endif /*GENERAL_INCLUDES_H_*/
