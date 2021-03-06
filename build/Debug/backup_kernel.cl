/*
 *  kernel.cl
 *  XclNet
 *
 *  Created by David Higgins on 26/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
/*typedef struct SynapseConstsStruct{
	float gamma_p;
	float gamma_d;
	float theta_p;
	float theta_d;
	unsigned int delay;
	float sigma;
	float tau;
	float tau_ca;
	float c_pre;
	float c_post;
	float dt;
	unsigned int no_syns;
} SynapseConsts;*/
 
#define PI (4.*atan(1.))

// Simple compute kernel which computes the square of an input array
//
__kernel void square(                                  
   __global float* input,
   __global float* output,
   const unsigned int count)
{
	int i = get_global_id(0);
	if(i < count)
		output[i] = input[i] * input[i];
}

// Second attempt at Marsaglia's generators
// This time with corrected equations!
// taken from http://www.math.niu.edu/~rusin/known-math/99/RNG
typedef struct random_struct2{
	float value;
	unsigned int d_z;
	unsigned int d_w;
	unsigned int d_jsr;
	unsigned int d_jcong;
} Random2;
void my_two(Random2 *rdm){
	// Multiply-with-carry
	(*rdm).d_z = 36969*((*rdm).d_z&65535)+((*rdm).d_z>>16);
	(*rdm).d_w = 18000*((*rdm).d_w&65535)+((*rdm).d_w>>16);
	unsigned int d_mwc = ((*rdm).d_z<<16)+(*rdm).d_w;
	
	// 3-shift-register generator
	(*rdm).d_jsr^=((*rdm).d_jsr<<17);
	(*rdm).d_jsr^=((*rdm).d_jsr>>13);
	(*rdm).d_jsr^=((*rdm).d_jsr<<5);
	//unsigned int d_shr3 = (*rdm).d_jsr;
	
	// Congruential generator
	(*rdm).d_jcong = 69069*(*rdm).d_jcong+1234567;
	//unsigned int d_cong = (*rdm).d_jcong;
	
	// KISS generator (combines above three generators)
	unsigned int d_kiss = ((d_mwc^(*rdm).d_jcong)+(*rdm).d_jsr);
	
	// Convert to Uniform(0,1) distribution
	float d_uni = (d_kiss*2.328306e-10);
	
	(*rdm).value = d_uni;
}
void my_GetNormal(Random2 *rdm)
{
	// Use Box-Muller algorithm
	my_two(rdm);
	float r = sqrt( -2.0*log((*rdm).value) );
	
	my_two(rdm); // Don't forget this second call!
	float theta = 2.0*PI*(*rdm).value;
	
	(*rdm).value = r*sin(theta);
}

typedef struct random_struct{
	float value;
	unsigned int m_z;
	unsigned int m_w;
} RandomStruct;
// Multiply-with-Carry random number generator
// static unsigned int m_w = 521288629, m_z = 362436069;
// Code basically from http://www.codeproject.com/Articles/25172/Simple-Random-Number-Generation
unsigned int GetUint(RandomStruct *rnd)
{
	(*rnd).m_z = 36969 * ((*rnd).m_z & 65535) + ((*rnd).m_z >> 16);
	(*rnd).m_w = 18000 * ((*rnd).m_w & 65535) + ((*rnd).m_w >> 16);
	return ((*rnd).m_z << 16) + (*rnd).m_w;
}
// Uniform distribution in interval (0,1)
float GetUniform(RandomStruct *rnd)
{
	// 0 <= u < 2^32
	unsigned int u = GetUint(rnd);
	// The magic number below is 1/(2^32 + 2).
	// The result is strictly between 0 and 1.
	(*rnd).value = (u + 1.0) * 2.328306435454494e-10;
	return (*rnd).value;
}
// Gaussian(0,1) distribution using Box-Muller algorithm
void GetNormal(RandomStruct *rnd)
{
	// Use Box-Muller algorithm
	float u1 = GetUniform(rnd);
	float u2 = GetUniform(rnd);
	float r = sqrt( -2.0*log(u1) );
	float theta = 2.0*PI*u2;
	(*rnd).value = r*sin(theta);
}

// Leaky integrate and fire kernel
//
__kernel void lif(
	__global float* input_v, // membrane voltage
	__global float* input_i, // input current
	__global float* input_gauss, // gaussian noise on membrane potential
	__global unsigned int* input_spike, // refractory period count down variable
	
	const float v_rest, // resting membrane voltage
	const float v_reset, // reset membrane voltage
	const float v_threshold, // threshold voltage for spiking
	const float r_m, // membrane resistance
	const float c_m, // membrane capacitance
	const float sigma, // size of noise
	const float refrac_time, // duration of refractory period
	const float dt, // time step size
	const unsigned int no_lifs // number of lifs in simulation	
	)
{
	int i = get_global_id(0);
	if ( i < no_lifs ){
		float v = input_v[i];
		float input_current = input_i[i];
		float gauss = input_gauss[i];
		unsigned int time_since_spike = input_spike[i];
	
		float new_v;
		float dv = 0;
		float noise = 0;
		float tau_m = r_m * c_m;
	
		//TODO: decide initial value for time_since_spike, otherwise system always resets to V_reset upon initialisation
		if (time_since_spike == 0){
			// A spike has just occurred, reset membrane voltage to reset potential
			v = v_reset;
		}
	
		// If refractory period is 0 OR if it's been longer than the refractory period since the last spike
		//CONSIDER: changed to >= to allow removal of logical OR which didn't work: (refrac_time==0)||
		if ( time_since_spike >= refrac_time ){
			// Apply leak current
			dv = (-(v - v_rest) / (c_m * r_m));
			// Apply the external current
			dv += (input_current / c_m);
			// Apply noise
			noise = sqrt(dt / tau_m) * sigma * gauss;
		}

		new_v = v + (dv * dt) + noise;
		// Apply lower threshold to membrane voltage
		if (new_v < v_rest){
			new_v = v_rest;
		}
	
		//Check if a spike has just occurred
		if (new_v > v_threshold){
			// A spike has just occurred, set time since last spike to 0
			time_since_spike = 0;
		}
		else{
			// No spike occurred, increment time since last spike
			time_since_spike++;
		}
	
		//new_v = uni_kiss(153512,22532,3167,43733);
		/*RandomStruct rnd;
		rnd.m_w = 521288629;
		rnd.m_z = 362436069;
		GetUniform(&rnd);
		new_v = rnd.value;*/
		input_spike[i] = time_since_spike;
		input_v[i] = new_v;
	}
}


	/*
	//TODO: opencl guarantees support for a minimum of only 8 const args, probably safer to pass params by passing a struct
	const float gamma_p, // potentiation learning rate
	const float gamma_d, // depression learning rate
	const float theta_p, // potentiation threshold
	const float theta_d, // depression threshold
	const float tau, // time constant for synaptic efficacy
	const float tau_ca, // time constant for calcium concentration
	const float c_pre, // increase in Ca by pre-synaptic spike
	const float c_post, // increase in Ca by post-synaptic spike
	const float sigma, // size of noise
	const float dt, // time step size	
	const unsigned int no_syns // number of synapses in simulation
	*/
// Graupner 2012 Synapse kernel
//
__kernel void synapse( 
	__global float* input_rho, // synaptic efficacy
	__global float* input_ca, // synaptic calcium concentration
	__global float* input_gauss, // gaussian noise on synaptic efficacy
	__global unsigned int* input_pre_spike, // number of pre-synaptic spikes occurring at time t-D
	__global unsigned int* input_post_spike, // number of post-synaptic spikes at time t
	//__const SynapseConsts syn_const
	const unsigned int no_syns
	)
{
	float gamma_p = 725.085;
	float gamma_d = 331.909;
	float theta_p = 1.3;
	float theta_d = 1.0;
	float sigma = 0; //3.35; //TODO: switch noise back on
	float tau = 346.3615;
	float tau_ca = 0.0226936;
	float c_pre = 0.5617539;
	float c_post = 1.23964;
	float dt = 0.001;
	
	int i = get_global_id(0);
	
	if (i < no_syns){
		float rho = input_rho[i];
		float ca = input_ca[i];
		float gauss = input_gauss[i];
		unsigned int pre_spike = input_pre_spike[i];
		unsigned int post_spike = input_post_spike[i];
	
		float new_rho, new_ca;
		float drho = 0, dca = 0;
		float noise = 0;
	
		unsigned int h_pot = 0;
		unsigned int h_dep = 0;
		unsigned int h_noise = 0;
	
		// Update Calcium concentration
		dca = (-ca/tau_ca) + (c_pre * pre_spike) + (c_post * post_spike);
		new_ca = ca + (dca * dt);
	
		// Update Synaptic efficacy
		if (new_ca > theta_p){
			h_pot = 1;
			h_noise = 1;
		}
		if (new_ca > theta_d){
			h_dep = 1;
			h_noise = 1;
		}

		// Calculate rho update
		drho = (-rho * (1.0 - rho) * (0.5 - rho)) + (gamma_p * (1.0 - rho) * h_pot) - (gamma_d * rho * h_dep);
		drho /= tau;
		// Calculate noise
		noise = (h_noise * sigma * sqrt(dt/tau) * gauss);
		// Calculate new rho value
		new_rho = rho + (drho * dt) + noise;
	
		input_rho[i] = new_rho;
		input_ca[i] = new_ca;
	}
}
