/*
 *  kernel.cl
 *  XclNet
 *
 *  Created by David Higgins on 26/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
 
#define PI (4.*atan(1.))

#include "Random123/philox.h"



// Leaky integrate and fire kernel
//
__kernel void lif(
	__global float* input_v, // membrane voltage
	__global float* input_i, // input current
	//__global float* input_gauss, // gaussian noise on membrane potential
	__global unsigned int* input_spike, // refractory period count up variable
	
	//State variables for random number generator
	/*__global unsigned int* d_z,
	__global unsigned int* d_w,
	__global unsigned int* d_jsr,
	__global unsigned int* d_jcong,*/
	
	const float v_rest, // resting membrane voltage
	const float v_reset, // reset membrane voltage
	const float v_threshold, // threshold voltage for spiking
	const float tau_m, // membrane time constant
	//const float c_m, // membrane capacitance
	const float sigma, // size of noise
	const float refrac_time, // duration of refractory period
	const float dt, // time step size
	const unsigned int no_lifs, // number of lifs in simulation
	
	const unsigned int time_step, // used for indexing the random number generator
	const unsigned int random_seed, // seed for the random number generator
	
	//TODO: if gauss stream permanently removed then it should be removed from here, etc.
	__global float* output_gauss 
	)
{
	int i = get_global_id(0);
	if ( i < no_lifs ){
		float v = input_v[i];
		float input_current = input_i[i];
		unsigned int time_since_spike = input_spike[i];
		
		philox2x32_key_t key;
		philox2x32_ctr_t ctr;
		philox2x32_ctr_t rand_val;
		float2 uni_rand;
		
		// Generate Gaussian(0,1) noise using Random123 library implementation		
		key.v[0] = i;
		ctr.v[0] = time_step;
		ctr.v[1] = random_seed;
		rand_val = philox2x32_R(10, ctr, key);
		// Convert to Uniform distribution (1/(2^32 +2))
		uni_rand.x = rand_val.v[0] * 2.328306435454494e-10;
		uni_rand.y = rand_val.v[1] * 2.328306435454494e-10;
		// Box-Muller transform
		float r = sqrt( (float)(-2.0*log(uni_rand.x)) );
		float theta = 2.0 * PI * uni_rand.y;
		float random_value = r * sin(theta);
		
	
		float new_v;
		float dv = 0;
		float noise = 0;
		//float tau_m = r_m * c_m;
		
	
		//REMINDER: initialise time_since_spike to refrac_time in main program,
		// otherwise system always resets to V_reset upon initialisation
		if (time_since_spike == 0){
			// A spike has just occurred, reset membrane voltage to reset potential
			v = v_reset;
		}
	
		// If refractory period is 0 OR if it's been longer than the refractory period since the last spike
		//CONSIDER: changed to >= to allow removal of logical OR which didn't work: (refrac_time==0)||
		if ( time_since_spike >= refrac_time ){
			// Apply leak current
			dv = (-(v - v_rest) / tau_m);
			// Apply the external current
			// Note: I use one input current variable (to cut down on streams to GPU)
			//  an external current/voltage should be added directly to this variable (outside the kernel)
			//  a synaptic current/voltage step should be multiplied by (tau_m/dt), for a delta spike, before adding to this variable,
			//  in order to counter rescaling which happens on next three lines of executable code.
			// input_current is treated as a voltage step, despite the variable name, hence the division by tau_m
            // Removing multiplication of input_current by dt, should also remove division in main program
			//dv += (input_current / tau_m);
			// Apply noise
			//noise = sqrt(dt / tau_m) * sigma * rnd.value;
			noise = sqrt(dt / tau_m) * sigma * random_value;
		}

		new_v = v + (dv * dt) + (input_current / tau_m) + noise;
		// Apply lower threshold to membrane voltage (no longer desired)
		/*if (new_v < v_rest){
			new_v = v_rest;
		}*/
	
		//Check if a spike has just occurred
		if (new_v > v_threshold){
			// A spike has just occurred, set time since last spike to 0
			time_since_spike = 0;
		}
		else{
			// No spike occurred, increment time since last spike
			time_since_spike++;
		}
		

		input_spike[i] = time_since_spike;
		input_v[i] = new_v;
		output_gauss[i] = random_value;
	}
}

