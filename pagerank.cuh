#include "cuda.h"

__global__ void kernel(int *cstart, int *cend, int *cmemsz, int *cmember, int *crcw, 
	double *cinitial, double *crank, int *rcwgraph, int *outdeg, int *corder, int *ctemp, int *ctempg)
{
	int w = blockIdx.x + (*cstart);
	if(w < (*cend)){
		int size = cmemsz[corder[w]];
		int j = threadIdx.x;
		if(j < size){
			int node = cmember[ctemp[corder[w]]+j];
			int k = threadIdx.y;
			int size1 = crcw[node];
			if(k < size1){
				atomicAdd(&cinitial[node], crank[rcwgraph[ctempg[node]+k]]/outdeg[rcwgraph[ctempg[node]+k]]);
			}
		}
	}
}

__global__ void update(double *cinitial, int *cn){
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int num_threads = blockDim.x * gridDim.x;
	for(int i=0;i<(*cn);i+=num_threads){
		int ver = i + tid;
		if(ver < (*cn)){
			cinitial[ver]=0.85*cinitial[ver];
			printf("%f\n",cinitial[ver] );
		}
	}
}