__global__ void kerneltest(long long *cstart, long long *cend, long long *cmemsz, long long *cmember, long long *crcw, 
	double *cinitial, double *crank, long long *rcwgraph, long long *outdeg, long long *corder, long long *ctemp, long long *ctempg)
{
	long long w = blockIdx.z*blockDim.z + threadIdx.z + (*cstart);
	long long num_threads_z = blockDim.z * gridDim.z;
	
	for(;w < (*cend);w+=num_threads_z)
	{
		long long size = cmemsz[w];
		long long j = blockIdx.x*blockDim.x+threadIdx.x;
		long long num_threads_x = blockDim.x * gridDim.x;
		for(;j<size;j+=num_threads_x){
			long long node = cmember[ctemp[w]+j];
			long long k = blockIdx.y*blockDim.y+threadIdx.y;
			long long size1 = crcw[node];
			long long num_threads_y = blockDim.y * gridDim.y;
			for(;k<size1;k+=num_threads_y){
				atomicAdd(&cinitial[node], 0.85*crank[rcwgraph[ctempg[node]+k]]/outdeg[rcwgraph[ctempg[node]+k]]);
			}
		}
	}
}

__global__ void kerneltest1(long long *cstart, long long *cend, long long *cmemsz, long long *cmember, long long *crcw, 
	double *cinitial, double *crank, long long *rcwgraph, long long *outdeg, long long *corder, long long *ctemp, long long *ctempg)
{
	long long w = (*cstart);
	long long size = cmemsz[w];
	long long j = blockIdx.x*blockDim.x+threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;j<size;j+=num_threads_x){
		long long node = cmember[ctemp[w]+j];
		long long k = blockIdx.y*blockDim.y+threadIdx.y;
		long long size1 = crcw[node];
		long long num_threads_y = blockDim.y * gridDim.y;
		for(;k<size1;k+=num_threads_y){
			atomicAdd(&cinitial[node], 0.85*crank[rcwgraph[ctempg[node]+k]]/outdeg[rcwgraph[ctempg[node]+k]]);
		}
	}
}

__global__ void kernel1test(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cparent)
{
	long long w = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;w<(*cn);w+=num_threads_x){
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;
		long long num_threads_y = blockDim.y * gridDim.y;
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

__global__ void kernel1test1(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cparent)
{
	long long w = *cn;
	long long size = csize[w];
	long long j = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;j<size;j+=num_threads_x){
		long long node = cgraph[ctemp[w]+j];
		atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
	}
}

__global__ void kernel2test(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cparent, long long *cmarked)
{
 
	long long w = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;w<(*cn);w+=num_threads_x){
		if(cmarked[w] != 0) continue; 
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;
		long long num_threads_y = blockDim.y * gridDim.y;	
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

__global__ void kernel2test1(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cparent, long long *cmarked)
{
 
	long long w = *cn;
	if(cmarked[w] == 0){
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;
		long long num_threads_y = blockDim.y * gridDim.y;	
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}
 
__global__ void kernel3test(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg)
{
	long long w = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;w<(*cn);w+=num_threads_x){
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;
		long long num_threads_y = blockDim.y * gridDim.y;
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}


__global__ void kernel3test1(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg)
{
	long long w = *cn;
	long long size = csize[w];
	long long j = blockIdx.y*blockDim.y + threadIdx.y;
	long long num_threads_y = blockDim.y * gridDim.y;
	for(;j<size;j+=num_threads_y){
		long long node = cgraph[ctemp[w]+j];
		atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
	}
}
 
__global__ void kernel4test(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cmarked)
{
	long long w = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;w<(*cn);w+=num_threads_x){
		if(cmarked[w] != 0) continue;
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;	
		long long num_threads_y = blockDim.y * gridDim.y;
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}

__global__ void kernel4test1(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cmarked)
{
	long long w = *cn;
	if(cmarked[w]==0){
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;	
		long long num_threads_y = blockDim.y * gridDim.y;
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}
