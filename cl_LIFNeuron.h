#ifndef CL_LIFNEURON_H_
#define CL_LIFNEURON_H_

#include "GeneralIncludes.h"

typedef struct LIFNeuron{
    float * V;
    float * I;
	float * gauss;
    unsigned int * time_since_spike;
	
	//new, for ISI recorder
	unsigned int * time_of_last_spike;
	
	unsigned int * no_outgoing_synapses;
	unsigned int * no_outgoing_ee_synapses;
	signed int ** outgoing_synapse_index;
	
	unsigned int * no_incoming_synapses;
	signed int ** incoming_synapse_index;
	
	unsigned int * no_incoming_exc_synapses;
	signed int ** incoming_lif_index;
	
	//unsigned int * no_incoming_exc_synapses;
	signed int ** outgoing_lif_index;

	float v_rest;
	float v_reset;
	float v_threshold;
	float tau_m_e;
	float tau_m_i;
	//float r_m;
	//float c_m;
	float sigma;
	float refrac_time_exc;  //TODO: why is this a float? I think it was as I wasn't sure whether to use timesteps or seconds here.
	float refrac_time_inh;
	float dt;
	unsigned int no_lifs;
	unsigned int no_exc;
	
	unsigned int time_step; // required for the random123 number generator
	unsigned int random123_seed;
    
    unsigned int time_next_stim_on;
    unsigned int time_next_stim_off;
    unsigned int stim_repeats;
	
	//unsigned char * subpopulation_flag; // manipulations will be performed on this population
	
	// All outgoing synaptic dynamics are identical so we keep them with their associated pre-synaptic neurons
	//synaptic dynamics variables
	float tau_ampa_decay;
	float tau_nmda_decay;
	float tau_gaba_decay;
	float tau_ampa_rise;
	float tau_nmda_rise;
	float tau_gaba_rise;
	float proportion_ampa; // versus nmda
	
	unsigned int spike_delay; // no of timesteps since a spike occurred before it gets added to synaptic dynamics
	// Warning: the above variable assumes that the ISI is never less than this delay
	
	float * s_ampa;
	float * x_ampa;
	float * s_nmda;
	float * x_nmda;
    float * s_gaba;
	float * x_gaba;
	float * H_exc_spike_input;
    float * H_inh_spike_input;
} cl_LIFNeuron;


#endif /*CL_LIFNEURON_H_*/
