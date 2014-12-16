/*
 *  DataReporters.c
 *  XclNet
 *
 *  Created by David Higgins on 09/04/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "DataReporters.h"

void reporters_setup(){
	char outfile[FILE_NAME_LENGTH];
	raster_name = "raster.dat";
	intracellular_name = "intracellular.dat";
	average_activity_name = "network_activity.dat";
	synaptic_activity_name = "single_synapse.dat";
	synaptic_strength_name = "final_synaptic_strength.dat";
	
	// Make sure that directory 'output' exists
	if(mkdir("output",(S_IRUSR | S_IWUSR | S_IXUSR)) == -1){
		if (errno == EEXIST){
			printf("Directory 'output' already exists.\n");
		}
		else{
			perror("Error creating directory 'output'");
		}
	}
	
	// Raster output file
	strcpy(outfile, "output/");
	strcat(outfile, raster_name);
	//printf("DEBUG: %s\n", outfile);
	raster_output = fopen(outfile, "a");
	if(raster_output == NULL){
		perror("Error: failed to open raster output file\n");
	}
	fprintf(raster_output, "\n\n\n\n\n# Raster output (t, lif_no)\n");
	
	// Intracellular recording from a single neuron
	strcpy(outfile, "output/");
	strcat(outfile, intracellular_name);
	//printf("DEBUG: %s\n", outfile);
	intracellular_output = fopen(outfile, "a");
	if(intracellular_output == NULL){
		perror("Error: failed to open intracellular output file\n");
	}
	lif_currents_EE = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_EI = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_IE = calloc(MAX_TIME_STEPS, sizeof(float));
	lif_currents_II = calloc(MAX_TIME_STEPS, sizeof(float));
	fprintf(intracellular_output, "\n\n\n\n\n# Intracellular recorder (t, V(t), time since spike(t), Iext(t), gauss(t-1), Itot(t)), currents:{EE,IE,EI,II}\n# Neuron ID: %d\n", RECORDER_NEURON_ID);
	
	// Population spiking activity
	strcpy(outfile, "output/");
	strcat(outfile, average_activity_name);
	//printf("DEBUG: %s\n", outfile);
	average_activity_ouput = fopen(outfile, "a");
	if(average_activity_ouput == NULL){
		perror("Error: failed to open average activity output file\n");
	}
	no_spiking_bins = (LIF_DT / BIN_SIZE) * MAX_TIME_STEPS;
	// Setup the bins for recording average population spiking behaviour
	summary_exc_spikes = calloc(no_spiking_bins, sizeof(double));
	summary_inh_spikes = calloc(no_spiking_bins, sizeof(double));
	// Code for monitoring selectively manipulated neurons
	lif_injection_spikes = calloc(no_spiking_bins, sizeof(double));
	//fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity (time bin, TotSpikes/N, ExcSpikes/NE, InhSpikes/NI, InstantaneousExcRate, InstantaneousInhRate, RhoAvSubPop, RhoStdevSubPop, RhoAv, RhoStdev, SupPopSelectiveStimExcRate, NoExcSpikes)\n# all normalised to their respective population sizes\n");
	//fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity (time bin, InstantaneousExcRate, InstantaneousInhRate, InstantaneousStimRate, RhoAvNonPop, RhoStdevNonPop, RhoAvStimPop, RhoStdevStimPop, RhoAvPrePop, RhoStdevPrePop, RhoAvPostPop, RhoStdevPostPop, NupdatesNonPop, NupdatesStimPop, NupdatesPrePop, NupdatesPostPop, RhoMaxNonPop, RhoMinNonPop, RhoMaxStimPop, RhoMinStimPop, RhoMaxPrePop, RhoMinPrePop, RhoMaxPostPop, RhoMinPostPop )\n# all normalised to their respective population sizes\n");
	
	#ifdef MONITOR_UP_DOWN_POPS
		fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity not all accurate?(time bin, InstantaneousExcRate, InstantaneousInhRate, InstantaneousStimRate, RhoAvNonPop, RhoStdevNonPop, RhoAvStimPop, RhoStdevStimPop, RhoAvPrePop, RhoStdevPrePop, RhoAvPostPop, RhoStdevPostPop, NupdatesNonPop, NupdatesStimPop, NupdatesPrePop, NupdatesPostPop, RhoMaxNonPop, RhoMinNonPop, RhoMaxStimPop, RhoMinStimPop, RhoMaxPrePop, RhoMinPrePop, RhoMaxPostPop, RhoMinPostPop )\n# all normalised to their respective population sizes\n");
		//New code for monitoring initially_UP population
		UP_pop_rho = calloc(no_spiking_bins, sizeof(double));
		UP_pop_M = calloc(no_spiking_bins, sizeof(double));
		UP_pop_S = calloc(no_spiking_bins, sizeof(double));
		UP_pop_n = calloc(no_spiking_bins, sizeof(unsigned int));
		UP_pop_max = calloc(no_spiking_bins, sizeof(float));
		UP_pop_min = calloc(no_spiking_bins, sizeof(float));
	
		//New code for monitoring initially_DOWN population
		DOWN_pop_rho = calloc(no_spiking_bins, sizeof(double));
		DOWN_pop_M = calloc(no_spiking_bins, sizeof(double));
		DOWN_pop_S = calloc(no_spiking_bins, sizeof(double));
		DOWN_pop_n = calloc(no_spiking_bins, sizeof(unsigned int));
		DOWN_pop_max = calloc(no_spiking_bins, sizeof(float));
		DOWN_pop_min = calloc(no_spiking_bins, sizeof(float));
	#endif /* MONITOR_UP_DOWN_POPS */
	
	#ifdef MONITOR_STIM_POPS
		//fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity  some nonsense here in sequence(time bin, InstantaneousExcRate, InstantaneousInhRate, InstantaneousStimRate, RhoAvNonPop, RhoStdevNonPop, RhoAvStimPop, RhoStdevStimPop, RhoAvPrePop, RhoStdevPrePop, RhoAvPostPop, RhoStdevPostPop, NupdatesNonPop, NupdatesStimPop, NupdatesPrePop, NupdatesPostPop, RhoMaxNonPop, RhoMinNonPop, RhoMaxStimPop, RhoMinStimPop, RhoMaxPrePop, RhoMinPrePop, RhoMaxPostPop, RhoMinPostPop )\n# all normalised to their respective population sizes\n");
    
        fprintf(average_activity_ouput, "\n\n\n\n\n# Summary network activity (time bin, InstantaneousExcRate, InstantaneousInhRate, InstantaneousStimRate, sequence:{RhoAv, RhoStdev, Nupdates, RhoMax, RhoMin} for {Stim, Non, Pre, Post} pops )\n# all normalised to their respective population sizes\n");
    
		//TODO: new code for monitoring multiple synapses here
		non_summary_rho = calloc(no_spiking_bins, sizeof(double));
		non_summary_M = calloc(no_spiking_bins, sizeof(double));
		non_summary_S = calloc(no_spiking_bins, sizeof(double));
		non_summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
		non_summary_max = calloc(no_spiking_bins, sizeof(float));
		non_summary_min = calloc(no_spiking_bins, sizeof(float));
		//TODO: new code for monitoring main population synapses here
		stim_summary_rho = calloc(no_spiking_bins, sizeof(double));
		stim_summary_M = calloc(no_spiking_bins, sizeof(double));
		stim_summary_S = calloc(no_spiking_bins, sizeof(double));
		stim_summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
		stim_summary_max = calloc(no_spiking_bins, sizeof(float));
		stim_summary_min = calloc(no_spiking_bins, sizeof(float));
		//TODO: new code for monitoring multiple synapses here
		pre_summary_rho = calloc(no_spiking_bins, sizeof(double));
		pre_summary_M = calloc(no_spiking_bins, sizeof(double));
		pre_summary_S = calloc(no_spiking_bins, sizeof(double));
		pre_summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
		pre_summary_max = calloc(no_spiking_bins, sizeof(float));
		pre_summary_min = calloc(no_spiking_bins, sizeof(float));
		//TODO: new code for monitoring main population synapses here
		post_summary_rho = calloc(no_spiking_bins, sizeof(double));
		post_summary_M = calloc(no_spiking_bins, sizeof(double));
		post_summary_S = calloc(no_spiking_bins, sizeof(double));
		post_summary_n = calloc(no_spiking_bins, sizeof(unsigned int));
		post_summary_max = calloc(no_spiking_bins, sizeof(float));
		post_summary_min = calloc(no_spiking_bins, sizeof(float));
	#endif /* MONITOR_STIM_POPS */
    
	// Detailed recording from single synapse
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_activity_name);
	//printf("DEBUG: %s\n", outfile);
	synaptic_activity_output = fopen(outfile, "a");
	if(synaptic_activity_output == NULL){
		perror("Error: failed to open synaptic activity output file\n");
	}
	fprintf(synaptic_activity_output, "\n\n\n\n\n# Single synapse recorder (t, rho(t), ca(t), preT(t), postT(t))\n# Synapse ID: %d\n", RECORDER_SYNAPSE_ID);
	
	// Final state of all dynamic synapses
	strcpy(outfile, "output/");
	strcat(outfile, synaptic_strength_name);
	//printf("DEBUG: %s\n", outfile);
	synaptic_strength_output = fopen(outfile, "a");
	if(synaptic_strength_output == NULL){
		perror("Error: failed to open synaptic strength output file\n");
	}
	fprintf(synaptic_strength_output, "\n\n\n\n\n# Final synaptic strengths (syn_id, pre_syn_lif_id, post_syn_lif_id, rho_initial, rho_final)\n");
	//fprintf(synaptic_strength_output, "# Receives high rate stim, Pre: %d and Post: %d\n", ((*syn).pre_lif[RECORDER_SYNAPSE_ID]<100?1,0), ((*syn).post_lif[RECORDER_SYNAPSE_ID]<100?1,0));
	
	#ifdef DEBUG_MODE_NETWORK
		//Debugging stuff
		lif_debug_name = "lif_debug.dat";
		// LIF connectivity and activity output file
		strcpy(outfile, "output/");
		strcat(outfile, lif_debug_name);
		//printf("DEBUG: %s\n", outfile);
		lif_debug_output = fopen(outfile, "a");
		if(lif_debug_output == NULL){
			perror("Error: failed to open lif debug output file\n");
		}
		fprintf(lif_debug_output, "\n\n\n\n\n# LIF debug output (neuron id, mean destination id, no outgoing synapses, sum of gauss, no outgoing EE syns,\n# {no within class, mean dest within class}:{EE,EI,IE,II},\n# no incomming:{EE,EI,IE,II}), in_sub_pop\n");
		lif_mean_destination = calloc(NO_LIFS, sizeof(float));
		lif_gauss_totals = calloc(NO_LIFS, sizeof(float));
		lif_debug_no_EE = calloc(NO_LIFS, sizeof(int));
		lif_debug_no_EI = calloc(NO_LIFS, sizeof(int));
		lif_debug_no_IE = calloc(NO_LIFS, sizeof(int));
		lif_debug_no_II = calloc(NO_LIFS, sizeof(int));
		lif_mean_dest_EE = calloc(NO_LIFS, sizeof(float));
		lif_mean_dest_EI = calloc(NO_LIFS, sizeof(float));
		lif_mean_dest_IE = calloc(NO_LIFS, sizeof(float));
		lif_mean_dest_II = calloc(NO_LIFS, sizeof(float));
		lif_in_EE = calloc(NO_LIFS, sizeof(int));
		lif_in_EI = calloc(NO_LIFS, sizeof(int));
		lif_in_IE = calloc(NO_LIFS, sizeof(int));
		lif_in_II = calloc(NO_LIFS, sizeof(int));
	#endif /* DEBUG_MODE_NETWORK */
}


// New version of print_raster_spike to include ISI in output
void print_raster_spike(int t, int lif_no, float isi){
	// A spike has occurred, add its occurrence to raster file
	// print inter-spike-interval too
	// t, lif_id, isi
	fprintf(raster_output, "%d %d %f\n", t, lif_no, isi);
}
/*void print_raster_spike(int t, int lif_no){
	// A spike has occurred, add its occurrence to raster file
	// t, lif_id
	fprintf(raster_output, "%d %d\n", t, lif_no);
}*/


void print_network_summary_activity(){
	// bin_id_ms, no_spikes, no_exc_spikes, no_inh_spikes, exc_freq, inh_freq
	printf("Outputting network summary activity\n");
	for(int i = 0; i < no_spiking_bins; i++){
		//printf("DEBUG: i %d\n", i);
		//fprintf(average_activity_ouput, "%d %f %f %f %f %f %f %f %f %f %f %d\n", i, ((summary_inh_spikes[i] + summary_exc_spikes[i] + lif_injection_spikes[i]) / NO_LIFS), (summary_exc_spikes[i] / (NO_EXC - no_injection_lifs)), (summary_inh_spikes[i] / NO_INH), ((summary_exc_spikes[i] / (NO_EXC - no_injection_lifs)) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)), summary_rho[i]/summary_n[i], sqrt(summary_S[i]/(summary_n[i]-1)),  pop_summary_rho[i]/pop_summary_n[i], sqrt(pop_summary_S[i]/(pop_summary_n[i]-1)), ((lif_injection_spikes[i] / no_injection_lifs) * (1.0 / BIN_SIZE)), (int) summary_exc_spikes[i] );
		#ifdef MONITOR_UP_DOWN_POPS
			fprintf(average_activity_ouput, "%d %f %f %f %f %f %d %f %f %f %f %d %f %f\n", i, ((summary_exc_spikes[i] / (NO_EXC - NO_STIM_LIFS)) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)), ((lif_injection_spikes[i] / NO_STIM_LIFS) * (1.0 / BIN_SIZE)), UP_pop_rho[i]/UP_pop_n[i], sqrt(UP_pop_S[i]/(UP_pop_n[i]-1)), UP_pop_n[i], UP_pop_max[i], UP_pop_min[i], DOWN_pop_rho[i]/DOWN_pop_n[i], sqrt(DOWN_pop_S[i]/(DOWN_pop_n[i]-1)), DOWN_pop_n[i], DOWN_pop_max[i], DOWN_pop_min[i] );
		#endif /* MONITOR_UP_DOWN_POPS */
		#ifdef MONITOR_STIM_POPS	
			fprintf(average_activity_ouput, "%d %f %f %f %f %f %d %f %f %f %f %d %f %f %f %f %d %f %f %f %f %d %f %f\n", i, ((summary_exc_spikes[i] / (NO_EXC - NO_STIM_LIFS)) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)), ((lif_injection_spikes[i] / NO_STIM_LIFS) * (1.0 / BIN_SIZE)), stim_summary_rho[i]/stim_summary_n[i], sqrt(stim_summary_S[i]/(stim_summary_n[i]-1)), stim_summary_n[i], stim_summary_max[i], stim_summary_min[i], non_summary_rho[i]/non_summary_n[i], sqrt(non_summary_S[i]/(non_summary_n[i]-1)), non_summary_n[i], non_summary_max[i], non_summary_min[i], /**/pre_summary_rho[i]/pre_summary_n[i], sqrt(pre_summary_S[i]/(pre_summary_n[i]-1)), pre_summary_n[i], pre_summary_max[i], pre_summary_min[i], /**/ post_summary_rho[i]/post_summary_n[i], sqrt(post_summary_S[i]/(post_summary_n[i]-1)), post_summary_n[i], post_summary_max[i], post_summary_min[i] );
		#endif /* MONITOR_STIM_POPS */
		//fprintf(average_activity_ouput, "%d %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %d %f %f\n", i, ((summary_exc_spikes[i] / (NO_EXC - NO_STIM_LIFS)) * (1.0 / BIN_SIZE)), ((summary_inh_spikes[i] / NO_INH) * (1.0 / BIN_SIZE)), ((lif_injection_spikes[i] / NO_STIM_LIFS) * (1.0 / BIN_SIZE)), non_summary_rho[i]/non_summary_n[i], sqrt(non_summary_S[i]/(non_summary_n[i]-1)),  stim_summary_rho[i]/stim_summary_n[i], sqrt(stim_summary_S[i]/(stim_summary_n[i]-1)),  pre_summary_rho[i]/pre_summary_n[i], sqrt(pre_summary_S[i]/(pre_summary_n[i]-1)),  post_summary_rho[i]/post_summary_n[i], sqrt(post_summary_S[i]/(post_summary_n[i]-1)), non_summary_n[i], stim_summary_n[i], pre_summary_n[i], post_summary_n[i], non_summary_max[i], non_summary_min[i], stim_summary_max[i], stim_summary_min[i], pre_summary_max[i], pre_summary_min[i], post_summary_max[i], post_summary_min[i], UP_pop_rho[i]/UP_pop_n[i], sqrt(UP_pop_S[i]/(UP_pop_n[i]-1)), UP_pop_n[i], UP_pop_max[i], UP_pop_min[i], DOWN_pop_rho[i]/DOWN_pop_n[i], sqrt(DOWN_pop_S[i]/(DOWN_pop_n[i]-1)), DOWN_pop_n[i], DOWN_pop_max[i], DOWN_pop_min[i] );
	}
}


void print_synapse_activity(int t, cl_Synapse *syn){
	// in event-based model preT and postT here are one timestep in past wrt other variables
	// t, rho, ca, preT, postT
	fprintf(synaptic_activity_output, "%d %f %f %d %d\n", t, (*syn).rho[RECORDER_SYNAPSE_ID], (*syn).ca[RECORDER_SYNAPSE_ID], (*syn).preT[RECORDER_SYNAPSE_ID], (*syn).postT[RECORDER_SYNAPSE_ID]);
}


void print_synapses_final_state(cl_Synapse *syn, SynapseConsts *syn_const){
	// syn_id, pre_lif_id, post_lif_id, rho_initial, rho_final
	for(int i = 0; i < (*syn_const).no_syns; i++){
		//fprintf(synaptic_strength_output, "%d %d %d %f %f %d\n", i, (*syn).pre_lif[i], (*syn).post_lif[i], (*syn).rho_initial[i], (*syn).rho[i], (*syn).receives_stimulation_flag[i]);
		fprintf(synaptic_strength_output, "%d %d %d %0.10f %0.10f\n", i, (*syn).pre_lif[i], (*syn).post_lif[i], (*syn).rho_initial[i], (*syn).rho[i]);
	}
}


void print_lif_debug(cl_LIFNeuron *lif){
	// lif_id, mean_dest_id, no_out_syns, gauss_total, no_out_EE_syns, no_EE, dest_EE, no_EI, dest_EI, no_IE, dest_IE, no_II, dest_II, in_EE, in_EI, in_IE, in_II 
	printf("\nLIF Debug: saving connection statistics to file...\n");
	for(int i = 0; i < (*lif).no_lifs; i++){
		//fprintf(lif_debug_output, "%d %f %d %f %d %d %f %d %f %d %f %d %f %d %d %d %d %d\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i], (*lif).no_outgoing_ee_synapses[i], lif_debug_no_EE[i], lif_mean_dest_EE[i], lif_debug_no_EI[i], lif_mean_dest_EI[i], lif_debug_no_IE[i], lif_mean_dest_IE[i], lif_debug_no_II[i], lif_mean_dest_II[i], lif_in_EE[i], lif_in_EI[i], lif_in_IE[i], lif_in_II[i], (*lif).subpopulation_flag[i]);
		fprintf(lif_debug_output, "%d %f %d %f %d %d %f %d %f %d %f %d %f %d %d %d %d\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i], (*lif).no_outgoing_ee_synapses[i], lif_debug_no_EE[i], lif_mean_dest_EE[i], lif_debug_no_EI[i], lif_mean_dest_EI[i], lif_debug_no_IE[i], lif_mean_dest_IE[i], lif_debug_no_II[i], lif_mean_dest_II[i], lif_in_EE[i], lif_in_EI[i], lif_in_IE[i], lif_in_II[i]);
		
		//printf("%d %f %d %f\n", i, lif_mean_destination[i], (*lif).no_outgoing_synapses[i], lif_gauss_totals[i]);
	}
}


void reporters_close(){
	fclose(raster_output);
	fclose(intracellular_output);
	fclose(average_activity_ouput);
	fclose(synaptic_activity_output);
	fclose(synaptic_strength_output);
	#ifdef DEBUG_MODE_NETWORK
		fclose(lif_debug_output);
	#endif /* DEBUG_MODE_NETWORK */
}
