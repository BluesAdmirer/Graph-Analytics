// cross edge (Edge which have both end nodes are in different Components) contribution
/*
	Example,
		Components (levelwise),
			0 (level0)
			1 2 (level1)
			3 4 5 (level2)
		Component edges
			0 -> 1
			0 -> 2
			1 -> 3
			1 -> 4
			2 -> 5
			1 -> 5
	For level1, cstart = 1
		    cend = 3
	let us say w = 2
	for component 2, for all the nodes which have incoming edges from component 0,
	we will initialize that nodes with contribution of the node from which edge is coming
*/
__global__ void kerneltest(long long *cstart, long long *cend, long long *cmemsz, long long *cmember, long long *crcw, 
	double *cinitial, double *crank, long long *rcwgraph, long long *outdeg, long long *corder, long long *ctemp, long long *ctempg)
{
	// w = Component number
	long long w = blockIdx.z*blockDim.z + threadIdx.z + (*cstart);
	// num_threads_z = max threads in z dimension
	long long num_threads_z = blockDim.z * gridDim.z;
	// w is on the same level i.e. w belongs to [cstart, cend)
	for(;w < (*cend);w+=num_threads_z)
	{
		// size = number of nodes in component w
		long long size = cmemsz[w];
		// j is the index of node
		long long j = blockIdx.x*blockDim.x+threadIdx.x;
		// num_thread_x = max threads in x dimension
		long long num_threads_x = blockDim.x * gridDim.x;
		for(;j<size;j+=num_threads_x){
			// node values of node
			// ctemp[w] = number of nodes skipped to get to the list of members of component w
			// ctemp[w] + j = jth node index of members of component w
			long long node = cmember[ctemp[w]+j];
			// k is the index of node from which node have incoming edge in the graph rcwgrah
			long long k = blockIdx.y*blockDim.y+threadIdx.y;
			// size1 = size of node's cross edge list 
			long long size1 = crcw[node];
			// num_threads_y = max threads in y dimension
			long long num_threads_y = blockDim.y * gridDim.y;
			for(;k<size1;k+=num_threads_y){
				// there is an cross edge (u -> node)
				// ctempg[node] = number of nodes skipped to get to the edge list of the node
				// u = rcwgraph[ctempg[node]+k] => kth node in rcw list of the node
				// node u's pagerank already been calculated so its contribution is fixed and it is stored in initial
				atomicAdd(&cinitial[node], 0.85*crank[rcwgraph[ctempg[node]+k]]/outdeg[rcwgraph[ctempg[node]+k]]);
			}
		}
	}
}


// This function works same as above kerneltest function
// One thing is changed
// Component is fixed (cstart/w)
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

// cn = pivot (node is changing)
__global__ void kernel1test(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cparent)
{
	// w = index of node
	long long w = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;w<(*cn);w+=num_threads_x){
		// size = size of adj. list of node at index w
		long long size = csize[w];
		long long j = blockIdx.y*blockDim.y + threadIdx.y;
		long long num_threads_y = blockDim.y * gridDim.y;
		for(;j<size;j+=num_threads_y){
			long long node = cgraph[ctemp[w]+j];
			atomicAdd(&ccurr[w], crank[cparent[node]]/coutdeg[node]);
		}
	}
}

// node is fixed
__global__ void kernel1test1(long long *cn, long long *csize, long long *cmem, long long *cgraph, 
							long long *ctemp, double *ccurr, double *crank, long long *coutdeg, long long *cparent)
{
	// w = index of node
	long long w = *cn;
	// size of adj. list of node at w
	long long size = csize[w];
	long long j = blockIdx.x*blockDim.x + threadIdx.x;
	long long num_threads_x = blockDim.x * gridDim.x;
	for(;j<size;j+=num_threads_x){
		// edge (u -> w) u = node 
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
