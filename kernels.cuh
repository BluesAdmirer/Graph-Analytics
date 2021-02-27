__global__ void kernel(int *cstart, int *cend, int *cmemsz, int *cmember, int *crcw, 
	double *cinitial, double *crank, int *rcwgraph, int *outdeg, int *corder, int *ctemp, int *ctempg)
{
	int w = blockIdx.z*blockDim.z + threadIdx.z + (*cstart);
	if(w < (*cend)){
		int size = cmemsz[corder[w]];
		int j = blockIdx.x*blockDim.x+threadIdx.x;
		if(j < size){
			int node = cmember[ctemp[corder[w]]+j];
			int k = blockIdx.y*blockDim.y+threadIdx.y;
			int size1 = crcw[node];
			if(k < size1){
				atomicAdd(&cinitial[node], 0.85*crank[rcwgraph[ctempg[node]+k]]/outdeg[rcwgraph[ctempg[node]+k]]);
			}
		}
	}
}

__global__ void kernel1(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg, int *cparent)
{
	int w = blockIdx.x*blockDim.x + threadIdx.x;
	if(w < (*cn)){
		int size = csize[w];
		int j = blockIdx.y*blockDim.y + threadIdx.y;
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

__global__ void kernel2(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg, int *cparent, int *cmarked)
{

	int w = blockIdx.x*blockDim.x + threadIdx.x;
	if(w < (*cn) && cmarked[w]==0){
		int size = csize[w];
		int j = blockIdx.y*blockDim.y + threadIdx.y;		
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

__global__ void kernel3(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg)
{
	int w = blockIdx.x*blockDim.x + threadIdx.x;
	if(w < (*cn)){
		int size = csize[w];
		int j = blockIdx.y*blockDim.y + threadIdx.y;
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}

__global__ void kernel4(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg, int *cmarked)
{
	int w = blockIdx.x*blockDim.x + threadIdx.x;
	if(w < (*cn) && cmarked[w]==0){
		int size = csize[w];
		int j = blockIdx.y*blockDim.y + threadIdx.y;	
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}