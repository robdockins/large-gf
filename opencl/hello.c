#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <stdio.h>
#include <stdlib.h>
 
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
 
#define MEM_SIZE (128)
#define MAX_SOURCE_SIZE (0x100000)

const char * get_error_string(cl_int err) {
  switch(err){
  case 0: return "CL_SUCCESS";
  case -1: return "CL_DEVICE_NOT_FOUND";
  case -2: return "CL_DEVICE_NOT_AVAILABLE";
  case -3: return "CL_COMPILER_NOT_AVAILABLE";
  case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
  case -5: return "CL_OUT_OF_RESOURCES";
  case -6: return "CL_OUT_OF_HOST_MEMORY";
  case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
  case -8: return "CL_MEM_COPY_OVERLAP";
  case -9: return "CL_IMAGE_FORMAT_MISMATCH";
  case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
  case -11: return "CL_BUILD_PROGRAM_FAILURE";
  case -12: return "CL_MAP_FAILURE";
    
  case -30: return "CL_INVALID_VALUE";
  case -31: return "CL_INVALID_DEVICE_TYPE";
  case -32: return "CL_INVALID_PLATFORM";
  case -33: return "CL_INVALID_DEVICE";
  case -34: return "CL_INVALID_CONTEXT";
  case -35: return "CL_INVALID_QUEUE_PROPERTIES";
  case -36: return "CL_INVALID_COMMAND_QUEUE";
  case -37: return "CL_INVALID_HOST_PTR";
  case -38: return "CL_INVALID_MEM_OBJECT";
  case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
  case -40: return "CL_INVALID_IMAGE_SIZE";
  case -41: return "CL_INVALID_SAMPLER";
  case -42: return "CL_INVALID_BINARY";
  case -43: return "CL_INVALID_BUILD_OPTIONS";
  case -44: return "CL_INVALID_PROGRAM";
  case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
  case -46: return "CL_INVALID_KERNEL_NAME";
  case -47: return "CL_INVALID_KERNEL_DEFINITION";
  case -48: return "CL_INVALID_KERNEL";
  case -49: return "CL_INVALID_ARG_INDEX";
  case -50: return "CL_INVALID_ARG_VALUE";
  case -51: return "CL_INVALID_ARG_SIZE";
  case -52: return "CL_INVALID_KERNEL_ARGS";
  case -53: return "CL_INVALID_WORK_DIMENSION";
  case -54: return "CL_INVALID_WORK_GROUP_SIZE";
  case -55: return "CL_INVALID_WORK_ITEM_SIZE";
  case -56: return "CL_INVALID_GLOBAL_OFFSET";
  case -57: return "CL_INVALID_EVENT_WAIT_LIST";
  case -58: return "CL_INVALID_EVENT";
  case -59: return "CL_INVALID_OPERATION";
  case -60: return "CL_INVALID_GL_OBJECT";
  case -61: return "CL_INVALID_BUFFER_SIZE";
  case -62: return "CL_INVALID_MIP_LEVEL";
  case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
  default: return "Unknown OpenCL error";
  }
}

void check_err( cl_int err, const char* activity  ) {
  if( err == 0 ) {
  } else {
    const char* msg = get_error_string(err);
    printf("OpenCL Error while %s : %s\n", activity, msg);
    exit(err);
  }
}

int main()
{
  cl_device_id device_id = NULL;
  cl_context context = NULL;
  cl_command_queue command_queue = NULL;
  cl_mem memobj = NULL;
  cl_program program = NULL;
  cl_kernel kernel = NULL;
  cl_platform_id platform_id = NULL;
  cl_uint ret_num_devices;
  cl_uint ret_num_platforms;
  cl_int ret;
  
  char string[MEM_SIZE];
 
  FILE *fp;
  char fileName[] = "./hello.cl";
  char *source_str;
  size_t source_size;
 
  /* Load the source code containing the kernel*/
  fp = fopen(fileName, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel.\n");
    exit(1);
  }
  source_str = (char*)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);
  
  /* Get Platform and Device Info */
  ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
  printf("Found %d platforms\n", ret_num_platforms);

  ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
  
  /* Create OpenCL context */
  context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
  check_err(ret, "creating context");
 
  /* Create Command Queue */
  command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
  check_err(ret, "creating command queue");

  /* Create Memory Buffer */
  memobj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, MEM_SIZE * sizeof(char), NULL, &ret);
  check_err(ret, "creating memory buffer");

  /* Create Kernel Program from the source */
  program = clCreateProgramWithSource(context, 1, (const char **)&source_str,
                                      (const size_t *)&source_size, &ret);
  check_err(ret, "creating program");
 
  /* Build Kernel Program */
  ret = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
  if( ret == CL_BUILD_PROGRAM_FAILURE ) {
    const size_t MAX_LOG = 10000;
    char log[MAX_LOG];
    size_t log_sz;
    ret = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
                                MAX_LOG, &log, &log_sz );
    check_err(ret, "getting program build info!");
    puts(log);
    exit(-1);
  } else {
    check_err(ret, "building program");
  }
 
  /* Create OpenCL Kernel */
  kernel = clCreateKernel(program, "hello", &ret);
  check_err(ret, "building 'hello' kernel");
 
  /* Set OpenCL Kernel Parameters */
  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&memobj);
  check_err(ret, "setting kernel parameter");
 
  /* Execute OpenCL Kernel */
  ret = clEnqueueTask(command_queue, kernel, 0, NULL,NULL);
  check_err(ret, "executing kernel");
 
  /* Copy results from the memory buffer */
  ret = clEnqueueReadBuffer(command_queue, memobj, CL_TRUE, 0,
                            MEM_SIZE * sizeof(char),string, 0, NULL, NULL);
  check_err(ret, "reading result buffer");
   
  /* Display Result */
  puts(string);
 
  /* Finalization */
  ret = clFlush(command_queue);
  check_err(ret, "flushing command queue");

  ret = clFinish(command_queue);
  check_err(ret, "finishing command queue");

  ret = clReleaseKernel(kernel);
  check_err(ret, "releasing kernel");

  ret = clReleaseProgram(program);
  check_err(ret, "releasing program");

  ret = clReleaseMemObject(memobj);
  check_err(ret, "releasing memory buffer");

  ret = clReleaseCommandQueue(command_queue);
  check_err(ret, "releasing command queue");

  ret = clReleaseContext(context);
  check_err(ret, "releasing context");
  
  free(source_str);
  return 0;
}
