/*
 *  HandleOpenCL.c
 *  XclNet
 *
 *  Created by David Higgins on 19/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "HandleOpenCL.h"


char* print_cl_errstring(cl_int err) {
    switch (err) {
        case CL_SUCCESS:                          return strdup("Success!");
        case CL_DEVICE_NOT_FOUND:                 return strdup("Device not found.");
        case CL_DEVICE_NOT_AVAILABLE:             return strdup("Device not available");
        case CL_COMPILER_NOT_AVAILABLE:           return strdup("Compiler not available");
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return strdup("Memory object allocation failure");
        case CL_OUT_OF_RESOURCES:                 return strdup("Out of resources");
        case CL_OUT_OF_HOST_MEMORY:               return strdup("Out of host memory");
        case CL_PROFILING_INFO_NOT_AVAILABLE:     return strdup("Profiling information not available");
        case CL_MEM_COPY_OVERLAP:                 return strdup("Memory copy overlap");
        case CL_IMAGE_FORMAT_MISMATCH:            return strdup("Image format mismatch");
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return strdup("Image format not supported");
        case CL_BUILD_PROGRAM_FAILURE:            return strdup("Program build failure");
        case CL_MAP_FAILURE:                      return strdup("Map failure");
        case CL_INVALID_VALUE:                    return strdup("Invalid value");
        case CL_INVALID_DEVICE_TYPE:              return strdup("Invalid device type");
        case CL_INVALID_PLATFORM:                 return strdup("Invalid platform");
        case CL_INVALID_DEVICE:                   return strdup("Invalid device");
        case CL_INVALID_CONTEXT:                  return strdup("Invalid context");
        case CL_INVALID_QUEUE_PROPERTIES:         return strdup("Invalid queue properties");
        case CL_INVALID_COMMAND_QUEUE:            return strdup("Invalid command queue");
        case CL_INVALID_HOST_PTR:                 return strdup("Invalid host pointer");
        case CL_INVALID_MEM_OBJECT:               return strdup("Invalid memory object");
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return strdup("Invalid image format descriptor");
        case CL_INVALID_IMAGE_SIZE:               return strdup("Invalid image size");
        case CL_INVALID_SAMPLER:                  return strdup("Invalid sampler");
        case CL_INVALID_BINARY:                   return strdup("Invalid binary");
        case CL_INVALID_BUILD_OPTIONS:            return strdup("Invalid build options");
        case CL_INVALID_PROGRAM:                  return strdup("Invalid program");
        case CL_INVALID_PROGRAM_EXECUTABLE:       return strdup("Invalid program executable");
        case CL_INVALID_KERNEL_NAME:              return strdup("Invalid kernel name");
        case CL_INVALID_KERNEL_DEFINITION:        return strdup("Invalid kernel definition");
        case CL_INVALID_KERNEL:                   return strdup("Invalid kernel");
        case CL_INVALID_ARG_INDEX:                return strdup("Invalid argument index");
        case CL_INVALID_ARG_VALUE:                return strdup("Invalid argument value");
        case CL_INVALID_ARG_SIZE:                 return strdup("Invalid argument size");
        case CL_INVALID_KERNEL_ARGS:              return strdup("Invalid kernel arguments");
        case CL_INVALID_WORK_DIMENSION:           return strdup("Invalid work dimension");
        case CL_INVALID_WORK_GROUP_SIZE:          return strdup("Invalid work group size");
        case CL_INVALID_WORK_ITEM_SIZE:           return strdup("Invalid work item size");
        case CL_INVALID_GLOBAL_OFFSET:            return strdup("Invalid global offset");
        case CL_INVALID_EVENT_WAIT_LIST:          return strdup("Invalid event wait list");
        case CL_INVALID_EVENT:                    return strdup("Invalid event");
        case CL_INVALID_OPERATION:                return strdup("Invalid operation");
        case CL_INVALID_GL_OBJECT:                return strdup("Invalid OpenGL object");
        case CL_INVALID_BUFFER_SIZE:              return strdup("Invalid buffer size");
        case CL_INVALID_MIP_LEVEL:                return strdup("Invalid mip-map level");
        default:                                  return strdup("Unknown");
    }
}


char* readKernelSource(char * filename){
	// Try loading the kernel from a source file
	FILE *f_kernel;
	
	f_kernel = fopen(filename, "rb");
	if (f_kernel == NULL){
		perror("Error in loading kernel from file");
	}
	
	fseek(f_kernel, 0, SEEK_END);
	long pos = ftell(f_kernel);
	fseek(f_kernel, 0, SEEK_SET);
	
	char *KernelSource = malloc(pos);
	fread(KernelSource, pos, 1, f_kernel);
	fclose(f_kernel);
	//printf("Source:\n %s", KernelSource);	
	
	return KernelSource;
}

int setupCL(CL *cl){
	// Create initial OpenCL context and queue
	//
	
	if( getPlatformIDs(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	if( connectToComputeDevice(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( createComputeContext(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	if( createCommandQueue(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	
	return !(EXIT_FAILURE);
}

int getPlatformIDs(CL *cl){
	// Get platform IDs
	//
	
	
	printf("getting platform IDs...\n");
	
	//TODO: more than one platform?
	(*cl).err = clGetPlatformIDs(1, &(*cl).platform, NULL);
	//printf("DEBUG platform id: %d\n", (*cl).platform);
	
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to get platform ID!\n%s\n", print_cl_errstring((*cl).err));
		return EXIT_FAILURE;
	}
	//clGetDeviceInfo((*cl).device_id, CL_DEVICE_MAX_CONSTANT_ARGS, sizeof(cl_ulong), &dev_info, NULL);
	//printf("Max no args: %d\n", (unsigned int)dev_info);
	//printf("test %d\n", !(EXIT_FAILURE));
	return !(EXIT_FAILURE);
}

int connectToComputeDevice(CL *cl){
	// Connect to a compute device
	//
	int gpu = USE_GPU;
	//cl_uint dev_info;
	
	printf("connecting to compute device...\n");

	printf("gpu: %d\n", gpu);
	
	//TODO: more than one device?
	cl_device_id local_device_ids[10];
	unsigned int local_num_devices;
	//(*cl).err = clGetDeviceIDs((*cl).platform, (gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU), 1, &(*cl).device_id, NULL);
	(*cl).err = clGetDeviceIDs((*cl).platform, (gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU), 10, local_device_ids, &local_num_devices);
	//printf("DEBUG device id: %d\n", (*cl).device_id);
	
	
	int device = 0;

	(*cl).device_id = local_device_ids[device];
	
	printf("%d devices found, connecting to device: %d\n", local_num_devices, device);
	fflush(stdout);
	
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to create a device group!\n%s\n", print_cl_errstring((*cl).err));
		return EXIT_FAILURE;
	}

	return !(EXIT_FAILURE);
}

int createComputeContext(CL *cl){
	// Create a compute context
    //
	
	printf("creating compute context...\n");
	
    (*cl).context = clCreateContext(0, 1, &(*cl).device_id, NULL, NULL, &(*cl).err);
    if (!(*cl).context)
    {
        printf("Error: Failed to create a compute context!\n%s\n", print_cl_errstring((*cl).err));
        return EXIT_FAILURE;
    }
	return !(EXIT_FAILURE);
}

int createCommandQueue(CL *cl){
	// Create a command commands
    //
	
	printf("creating command queue...\n");

    (*cl).commands = clCreateCommandQueue((*cl).context, (*cl).device_id, 0, &(*cl).err);
    if (!(*cl).commands)
    {
        printf("Error: Failed to create a command commands!\n%s\n", print_cl_errstring((*cl).err));
        return EXIT_FAILURE;
    }
	return !(EXIT_FAILURE);
}

int makeProgram(CL *cl, char* KernelSource, char* k_name){
	// Make a compute kernel from source
	//
	
	//printf("Making program from kernel source:\n");
	//printf("%s", KernelSource);
	
	if( createProgram(cl, &KernelSource) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( buildProgram(cl) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}

	if( createKernel(cl, k_name) == EXIT_FAILURE){
		return EXIT_FAILURE;
	}
	return !(EXIT_FAILURE);
}

int createProgram(CL *cl, char ** KernelSource){
    // Create the compute program from the source buffer
    //
	
	printf("creating program from source...\n");
	
    (*cl).program = clCreateProgramWithSource((*cl).context, 1, (const char **) KernelSource, NULL, &(*cl).err);
    if (!(*cl).program)
    {
        printf("Error: Failed to create compute program!\n%s\n", print_cl_errstring((*cl).err));
        return EXIT_FAILURE;
    }
	return !(EXIT_FAILURE);
}

int buildProgram(CL *cl){
	// Build the program executable
	//
	
	printf("building program...\n");
	
	char* options = "-I /home/dhiggins/include/";
	//char* options = "";
	printf("build options: %s\n", options);
	

	(*cl).err = clBuildProgram((*cl).program, 0, NULL, options, NULL, NULL);
	if ((*cl).err != CL_SUCCESS)
	{
		size_t len;
		char buffer[2048];
		
		printf("Error: Failed to build program executable!\n%s\n", print_cl_errstring((*cl).err));
		clGetProgramBuildInfo((*cl).program, (*cl).device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		printf("Error code: %d\n", (*cl).err);
		exit(1);
	}
	return !(EXIT_FAILURE);
}

int createKernel(CL *cl, char * k_name){
	// Create the compute kernel in the program we wish to run
	//
	
	printf("creating kernel (%s)...\n", k_name);
	
	(*cl).kernel = clCreateKernel((*cl).program, k_name, &(*cl).err);
	if (!(*cl).kernel || (*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to create compute kernel!\n%s\n", print_cl_errstring((*cl).err));
		exit(1);
	}
	return !(EXIT_FAILURE);
}


int createLifIObufs(CL *cl){
	// Create IO buffers for transferring data to and from kernel
	//
	
	printf("creating LIF io buffers...\n");
	
    // Create the input and output arrays in device memory for our calculation
    //
    cl_int err1 = 0;
    cl_int err2 = 0;
    cl_int err3 = 0;
    cl_int err4 = 0;
    
	// Mapped memory with returned error values for debugging if a failure occurs
    (*cl).input_v = clCreateBuffer((*cl).context, CL_MEM_ALLOC_HOST_PTR,  sizeof(float) * (*cl).job_size, NULL, &err1); // read-write
	(*cl).input_spike = clCreateBuffer((*cl).context, CL_MEM_ALLOC_HOST_PTR,  sizeof(unsigned int) * (*cl).job_size, NULL, &err2); // read-write
    
    (*cl).input_current = clCreateBuffer((*cl).context, CL_MEM_ALLOC_HOST_PTR,  sizeof(float) * (*cl).job_size, NULL, &err3); // write only
	
	(*cl).gauss = clCreateBuffer((*cl).context,  CL_MEM_ALLOC_HOST_PTR,  sizeof(float) * (*cl).job_size, NULL, &err4); // read only
	

    if (!(*cl).input_v || !(*cl).input_current || !(*cl).gauss || !(*cl).input_spike)
    {
        printf("Error: Failed to allocate device memory!\n%s\n%s\n%s\n%s\n", print_cl_errstring(err1), print_cl_errstring(err2), print_cl_errstring(err3), print_cl_errstring(err4));
	    exit(1);
	}
	return !(EXIT_FAILURE);
}


int mapLifIObufs(CL *cl, cl_LIFNeuron *lif){
	// Create memory maps to IO buffers
	//
	
	printf("creating maps to LIF io buffers...\n");
    
    // I'm not sure that reusing err in the following commands is safe, do they reset it on each call? So I use separate variables instead.
    cl_int err1 = 0;
    cl_int err2 = 0;
    cl_int err3 = 0;
    cl_int err4 = 0;
	
    // Mapped memory is pinned (prevented from being swapped) hence faster (on occasion)
    printf("DEBUG: beginning map operation\n");

    (*lif).V = clEnqueueMapBuffer( (*cl).commands, (*cl).input_v , CL_TRUE,  (CL_MAP_READ | CL_MAP_WRITE), 0, sizeof(cl_float) * (*lif).no_lifs, 0, NULL, NULL, &err1 );
    (*lif).I = clEnqueueMapBuffer( (*cl).commands, (*cl).input_current , CL_TRUE,  (CL_MAP_WRITE), 0, sizeof(cl_float) * (*lif).no_lifs, 0, NULL, NULL, &err2);
    (*lif).time_since_spike = clEnqueueMapBuffer( (*cl).commands, (*cl).input_spike , CL_TRUE,  (CL_MAP_READ | CL_MAP_WRITE), 0, sizeof(cl_float) * (*lif).no_lifs, 0, NULL, NULL, &err3);
    (*lif).gauss = clEnqueueMapBuffer( (*cl).commands, (*cl).gauss , CL_TRUE,  (CL_MAP_READ), 0, sizeof(cl_float) * (*lif).no_lifs, 0, NULL, NULL, &err4 );
    printf("DEBUG: maps created\n");
    

    if (!(*lif).V || !(*lif).I|| !(*lif).time_since_spike || !(*lif).gauss)
    {
        printf("Error: Failed to create maps to device memory!\n%s\n%s\n%s\n%s\n", print_cl_errstring(err1), print_cl_errstring(err2), print_cl_errstring(err3), print_cl_errstring(err4));
	    exit(1);
	}
	return !(EXIT_FAILURE);
}


int enqueueLifInputBuf(CL *cl, cl_LIFNeuron *lif, cl_MarsagliaStruct *rnd){
	// Enqueue data for copying to Input buffers
	//
	
	//printf("enqueueing LIF input buffer...\n");
	
    // Write our data set into the input array in device memory
    //

	(*cl).err = clEnqueueWriteBuffer((*cl).commands, (*cl).input_current, CL_TRUE, 0, sizeof(float) * (*lif).no_lifs, (*lif).I, 0, NULL, NULL);

    
    if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to write to source array!\n%s\n", print_cl_errstring((*cl).err));
        exit(1);
    }
	return !(EXIT_FAILURE);
}


int setLifKernelArgs(CL *cl, cl_LIFNeuron *lif){
	// Set the Kernel arguments
	//
	
	//printf("setting args for LIF compute kernel...\n");
	
    // Set the arguments to our compute kernel
    //
    (*cl).err = 0;
    (*cl).err  = clSetKernelArg((*cl).kernel, 0, sizeof(cl_mem), &(*cl).input_v);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 1, sizeof(cl_mem), &(*cl).input_current);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 2, sizeof(cl_mem), &(*cl).input_spike);
	

	
	(*cl).err  |= clSetKernelArg((*cl).kernel, 3, sizeof(float), &(*lif).v_rest);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 4, sizeof(float), &(*lif).v_reset);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 5, sizeof(float), &(*lif).v_threshold);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 6, sizeof(float), &(*lif).tau_m);

	(*cl).err  |= clSetKernelArg((*cl).kernel, 7, sizeof(float), &(*lif).sigma);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 8, sizeof(float), &(*lif).refrac_time);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 9, sizeof(float), &(*lif).dt);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 10, sizeof(unsigned int), &(*lif).no_lifs);
	
	(*cl).err  |= clSetKernelArg((*cl).kernel, 11, sizeof(unsigned int), &(*lif).time_step);
	(*cl).err  |= clSetKernelArg((*cl).kernel, 12, sizeof(unsigned int), &(*lif).random123_seed);
	
	(*cl).err  |= clSetKernelArg((*cl).kernel, 13, sizeof(cl_mem), &(*cl).gauss);
	
    if ((*cl).err != CL_SUCCESS)
    {
        printf("Error: Failed to set kernel arguments!\n%s\n", print_cl_errstring((*cl).err));
        exit(1);
    }
	return !(EXIT_FAILURE);
}


int getMaxWorkSize(CL *cl){
	// Get the maximum work group size for executing the kernel on the device
	// and pad global work size such that it is a multiple of local
	//
	//size_t my_local_var;
	
	printf("getting max work group size...\n");
	
	(*cl).err = clGetKernelWorkGroupInfo((*cl).kernel, (*cl).device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof((*cl).local), &(*cl).local, NULL);

	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to retrieve kernel work group info! %d\n", (*cl).err);
		exit(1);
	}
	printf("CL_KERNEL_WORK_GROUP_SIZE: %d\n", (int)(*cl).local);

	
	(*cl).global = (*cl).job_size;
	(*cl).local = fmin((*cl).local, (*cl).global); // Copes with case global < local
	while( (*cl).global % (*cl).local != 0){
		// Pad the global number of work items such that it is divided evenly by local
		(*cl).global++;
		//printf("New value for global: %d\n", (*cl).global);
	}
	printf("Setting global work size: %d, local work group size: %d, real no jobs %d\n", (int)(*cl).global, (int)(*cl).local, (*cl).job_size);
	
	return !(EXIT_FAILURE);
}


int enqueueLifKernel(CL *cl){
	// Execute the kernel over the entire range of our 1d input data set
	// using the maximum number of work group items for this device
	//
	
	//printf("sending the LIF kernel to the process queue...\n");
	if((*cl).err){
		printf("Error already occurred\n%s\n", print_cl_errstring((*cl).err));
	}
    
    //TODO: modified enqueue kernel to return an event, which will then be waited for and cleared (hopefully a fix for nvidia)
	(*cl).err = clEnqueueNDRangeKernel((*cl).commands, (*cl).kernel, 1, NULL, &(*cl).global, &(*cl).local, 0, NULL, NULL);
	
    
	if ((*cl).err)
	{
		printf("Error: Failed to execute kernel!\n%s\n", print_cl_errstring((*cl).err));
		return EXIT_FAILURE;
	}
    
    
	return !(EXIT_FAILURE);
}


int waitForKernel(CL *cl){
	// Wait for the command commands to get serviced before reading back results
	//
	
	//printf("waiting for kernel to finish...\n");
	
	(*cl).err = clFinish((*cl).commands);
	if ((*cl).err)
	{
		printf("Error: Failed to finish command queue!\n%s\n", print_cl_errstring((*cl).err));
		return EXIT_FAILURE;
	}
	return !(EXIT_FAILURE);
}


int enqueueLifOutputBuf(CL *cl, cl_LIFNeuron *lif, cl_MarsagliaStruct *rnd){
	// Enqueue Kernel outputs for reading from buffers to system memory
	//
	

    // Read these memory buffers on each kernel run
    (*cl).err = clEnqueueReadBuffer( (*cl).commands, (*cl).input_v, CL_TRUE, 0, sizeof(float) * (*lif).no_lifs, (*lif).V, 0, NULL, NULL );
    (*cl).err |= clEnqueueReadBuffer( (*cl).commands, (*cl).input_spike, CL_TRUE, 0, sizeof(unsigned int) * (*lif).no_lifs, (*lif).time_since_spike, 0, NULL, NULL );
    
    
	if ((*cl).err != CL_SUCCESS)
	{
		printf("Error: Failed to read output array!\n%s\n", print_cl_errstring((*cl).err));
		exit(1);
	}
	return !(EXIT_FAILURE);
}


void shutdownLifKernel(CL *cl){
	// Shutdown and cleanup
	//
	
	printf("shutting down and cleaning up LIF kernel memory...");
	
	clReleaseMemObject((*cl).input_v);
	clReleaseMemObject((*cl).input_current);
	clReleaseMemObject((*cl).gauss);
	clReleaseMemObject((*cl).input_spike);
	

	clReleaseProgram((*cl).program);
	clReleaseKernel((*cl).kernel);
	clReleaseCommandQueue((*cl).commands);
	clReleaseContext((*cl).context);
	
	printf("done\n");
}



