/*
 *  DataReporters.h
 *  XclNet
 *
 *  Created by David Higgins on 09/04/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "GeneralIncludes.h"
#include "cl_Synapse.h"

//Debug
#include "cl_LIFNeuron.h"


// File pointers for recording data
// Raster output file
char* raster_name; //"raster.dat"
// Intracellular recorder file
char* intracellular_name; //"intracellular.dat"
// Averaged population activity file
char* average_activity_name;// = "network_activity.dat";
// Individual synaptic activity file
char* synaptic_activity_name;// = "single_synapse.dat";
// Final synaptic strengths of all dynamic synapses file
char* synaptic_strength_name;// = "final_synaptic_strength.dat";


FILE *raster_output;
FILE *intracellular_output;
FILE *average_activity_ouput; //network_activity_output
FILE *synaptic_activity_output; //single_synapse_output
FILE *synaptic_strength_output; //final_synaptic_strength_output


// Summary variables for monitoring network firing rate
unsigned int no_spiking_bins;
//CONSIDER: since we use a timestepping approach these variables could be condensed
// to single value variables and printed out during the simulation
double *summary_exc_spikes;
double *summary_inh_spikes;
// Variables for manipulating subset of neurons
double *lif_injection_spikes;
//int no_injection_lifs;

#ifdef MONITOR_UP_DOWN_POPS
//#define MONITOR_UP_DOWN_POPS /* local definition */
//Summary variablse for initially_UP pop
double *UP_pop_rho;
double *UP_pop_M;
double *UP_pop_S;
unsigned int *UP_pop_n;
float *UP_pop_max;
float *UP_pop_min;

//Summary variablse for initially_UP pop
double *DOWN_pop_rho;
double *DOWN_pop_M;
double *DOWN_pop_S;
unsigned int *DOWN_pop_n;
float *DOWN_pop_max;
float *DOWN_pop_min;
#endif /* MONITOR_UP_DOWN_POPS */

#ifdef MONITOR_STIM_POPS
//#define MONITOR_STIM_POPS /* local definition */
// Summary variables for monitoring synapses which receive only background activity
double *non_summary_rho;
double *non_summary_M;
double *non_summary_S;
unsigned int *non_summary_n;
float *non_summary_max;
float *non_summary_min;
// Summary variables for synapses which receive high stim pre and post activity
double *stim_summary_rho;
double *stim_summary_M;
double *stim_summary_S;
unsigned int *stim_summary_n;
float *stim_summary_max;
float *stim_summary_min;
// Summary variables for synapses which receive high stim pre activity
double *pre_summary_rho;
double *pre_summary_M;
double *pre_summary_S;
unsigned int *pre_summary_n;
float *pre_summary_max;
float *pre_summary_min;
// Summary variables for synapses which receive high stim post activity
double *post_summary_rho;
double *post_summary_M;
double *post_summary_S;
unsigned int *post_summary_n;
float *post_summary_max;
float *post_summary_min;
#endif /* MONITOR_STIM_POPS */

//Debugging variables
float *lif_gauss_totals;
float *lif_mean_destination;
char* lif_debug_name;
FILE *lif_debug_output;
int *lif_debug_no_EE;
int *lif_debug_no_IE;
int *lif_debug_no_EI;
int *lif_debug_no_II;
float *lif_mean_dest_EE;
float *lif_mean_dest_IE;
float *lif_mean_dest_EI;
float *lif_mean_dest_II;
int *lif_in_EE;
int *lif_in_IE;
int *lif_in_EI;
int *lif_in_II;
float *lif_currents_EE;
float *lif_currents_IE;
float *lif_currents_EI;
float *lif_currents_II;

void reporters_setup();
void reporters_close();

void print_raster_spike(int t, int lif_no, float isi);
void print_network_summary_activity();
void print_synapse_activity(int t, cl_Synapse *syn);
void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const);
void print_lif_debug(cl_LIFNeuron *lif);