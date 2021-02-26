__global__ void kernel(int *cstart, int *cend, int *cmemsz, int *cmember, int *crcw, 
	double *cinitial, double *crank, int *rcwgraph, int *outdeg, int *corder, int *ctemp, int *ctempg)
{
	// cstart = par[i]
	// cend = par[i+1]
	// w = from cstart to cend 
	int w = blockIdx.x + (*cstart);
	if(w < (*cend)){
		// size = number of nodes in SCC numbered w
		int size = cmemsz[corder[w]];
		// j index of node in SCC w
		int j = threadIdx.x;
		if(j < size){
			int node = cmember[ctemp[corder[w]]+j];
			int k = threadIdx.y;
			// kth node in rcw graph of node
			int size1 = crcw[node];
			if(k < size1){
				// update cinitial
				atomicAdd(&cinitial[node], 0.85*crank[rcwgraph[ctempg[node]+k]]/outdeg[rcwgraph[ctempg[node]+k]]);
			}
		}
	}
}

__global__ void kernel1(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg, int *cparent)
{
	int w = blockIdx.x;
	if(w < (*cn)){
		int size = csize[w];
		int j = threadIdx.x;
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			// parent node ki rank value and wo node ki degree
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

__global__ void kernel2(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg, int *cparent, int *cmarked)
{
	int w = blockIdx.x;
	// cmarked == dead node or not
	if(w < (*cn) && cmarked[w]==0){
		int size = csize[w];
		int j = threadIdx.x;
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			// parent node ki rank value and wo node ki degree
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

__global__ void kernel3(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg)
{
	// simple code
	int w = blockIdx.x;
	if(w < (*cn)){
		int size = csize[w];
		int j = threadIdx.x;
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}

__global__ void kernel4(int *cn, int *csize, int *cmem, int *cgraph, 
							int *ctemp, double *ccurr, double *crank, int *coutdeg, int *cmarked)
{
	int w = blockIdx.x;
	// only dead node combination
	if(w < (*cn) && cmarked[w]==0){
		int size = csize[w];
		int j = threadIdx.x;
		if(j < size){
			int node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[node]/coutdeg[node]);
		}
	}
}
