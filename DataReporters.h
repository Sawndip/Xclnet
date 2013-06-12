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

//TODO: rename file pointer variable names for DataReporters
FILE *raster_output;
FILE *intracellular_output;
FILE *average_activity_ouput; //network_activity_output
FILE *synaptic_activity_output; //single_synapse_output
FILE *synaptic_strength_output; //final_synaptic_strength_output

// Summary variables for monitoring network firing rate
//CONSIDER: since we use a timestepping approach these variables could be condensed
// to single value variables and printed out during the simulation
float *summary_exc_spikes;
float *summary_inh_spikes;
unsigned int no_spiking_bins;

// Variables for manipulating subset of neurons
float *lif_injection_spikes;
int no_injection_lifs;


// Summary variables for monitoring synapses which receive only background activity
float *non_summary_rho;
float *non_summary_M;
float *non_summary_S;
unsigned int *non_summary_n;
// Summary variables for synapses which receive high stim pre and post activity
float *stim_summary_rho;
float *stim_summary_M;
float *stim_summary_S;
unsigned int *stim_summary_n;
// Summary variables for synapses which receive high stim pre activity
float *pre_summary_rho;
float *pre_summary_M;
float *pre_summary_S;
unsigned int *pre_summary_n;
// Summary variables for synapses which receive high stim post activity
float *post_summary_rho;
float *post_summary_M;
float *post_summary_S;
unsigned int *post_summary_n;


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

void print_raster_spike(int t, int lif_no);
void print_network_summary_activity();
void print_synapse_activity(int t, cl_Synapse *syn);
void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const);
void print_lif_debug(cl_LIFNeuron *lif);