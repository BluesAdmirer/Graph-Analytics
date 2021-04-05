#include "kernels.cuh"
using namespace std;

// graph = graph with reversed edges
// parent = to store the parent of identical nodes
// left = identical nodes
// mapit = member vector = (mapit[i] = members of component[i] excluding identical nodes)
float computeparalleli(vector<vector<long long>> &graph, long long *parent, vector<long long> left, long long n, long long *outdeg, vector<long long> &mapit, double *rank,double *initial, long long nn)
{
	// total = time taken by kernel;
	float total = 0.0;
	long long i, iterations = 0;
	// thres = max allowed error
	double damp=0.85, thres=1e-10, error = 0;
	double randomp=(1-damp)/graph.size();
	long long thresh=10000; // to sort the nodes
	long long pivot=0;
	// sorting
	// o to pivot => nodes with smaller edge list ( < thres)
	// pivot to n => nodes with higher edge list ( > thres)
	for(i=0;i<n;i++)
	{
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}
	// curr ranks 
	double *curr = (double *)malloc(n*sizeof(double));
	for(long long i=0;i<n;i++){
		curr[i]=0;
	}
	// temp, sz, graphh == CSR
	// temp = prefix sum of sz
	// sz = size of graph[node]
	// graphh = edge list
	// mem = vector to int array of mapit
	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}
	
	long long *cn, *cm, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&cparent, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);

	do  
	{   
		// 0 to pivot calculation
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}

		// ccurr = array of size n(number of members of component i) initialized with 0
		// ccurr = stores the ranks computed within the component with the help of crank array
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);
		
		// cn = pivot
		// csize = size of adjacency list of nodes
		// cmem = members of the component
		// cgraph = edge list
		// ctemp = prefix sum of csize
		// coutdeg = outdegree of nodes
		// cparent = parents of identical nodes
		kernel1test<<<blockB ,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;

		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<n;i++)
		{   
			// pivot to n calculation (separate calculation for all the nodes as adj. list size is high)
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);
				
				// cm = ith node
				kernel1test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		}   

		double anse=0;
		// anse = max error
		for(i=0;i<n;i++){
			anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		}
		
		// set the rank values
		for(i=0;i<n;i++)
		{   
			{
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}   
		}
		iterations++;
		error = anse;
	}while(error > thres);
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(cparent);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(crank);
	return total;
}

void computeranki(vector < vector < long long > > & graph, long long *parent,vector < long long > left,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial, long long nn)
{
	double damp=0.85;
	double thres=1e-10;
	long long i,j;
	vector < double > curr(n);
	double error=0;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			{
				long long node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++){
					// if node is identical to another parent[node]=node
					// else parent[node] = u (u is parent)
					// as identical node's pagerank is not calculated and assigned values afterwards,
					// we use parent[graph[node][j]] instead graph[node][j] because graph[node][j] might be identical node
					ans=ans+rank[parent[graph[node][j]]]/outdeg[graph[node][j]];
				}
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
		}
		for(i=0;i<n;i++)
		{
			{
				rank[mapit[i]]=curr[i];
			}
		}
	}while(error > thres);
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	curr.clear();
}

float computeparallelid(vector < vector < long long > > & graph,long long *parent,vector < long long > & left,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial, long long nn)
{
	// It is same calculation as computeparalleli with one change
	// as dead node computation is included
	// we use marked array = to store the nodes which are dead nodes 
	// we calculate the pagerank for nodes which are not marked
	float total = 0.0;
	double thres=1e-10;
	double dis=1e-12;
	double value=((dis)*10.0)/n;
	long long i;
	double *prev = (double *)malloc(n*sizeof(double));
	double *curr = (double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	long long *marked = (long long *)malloc(n*sizeof(long long));
	memset(marked,0,n*sizeof(long long));
	double error=0;
	long long iterations=0;
	double damp = 0.85;
	double randomp=(1-damp)/graph.size();

	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long pivot=0;
	long long thresh=10000;
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}

	long long *cn, *cm, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph, *cmarked;;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&cparent, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&cmarked, n*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);

	do  
	{
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(long long), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel2test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<n;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel2test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		}  

		double anse=0;
		for(i=0;i<n;i++){
			if(!marked[i]){ // if not dead node
				anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
			}
		}
		iterations++;
		for(i=0;i<n;i++){
			if(!marked[i])   { // if not dead node
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}   
		}
		if(iterations%20==0){   
			// after every 20 every operations
			// we check whether any node has been dead or not
			// By dead, we mean that whether its pagerank is changing significantly or not
			for(i=0;i<n;i++){   
				if(!marked[i]){ // if not dead
					if(fabs(prev[i]-curr[i]) < value) // if difference between pagerank of previous iteration and next iteration is < value mark dead.
						marked[i]=1;
					else
						prev[i]=curr[i];
				}   
			}   
		}   
		error = anse;
	}while(error > thres );
	// left node computations
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(cparent);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(crank);
	cudaFree(cmarked);
	return total;
}

void computerankid(vector < vector < long long > > & graph,long long *parent,vector < long long > & left, long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial, long long nn)
{
	double damp = 0.85;
	double thres=1e-10;
	long long i,j;
	vector < double > curr(n);
	vector < double > prev(n);
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	double dis=1e-12;
	double value=(dis*10.0)/n;
	bool *marked = (bool *)malloc(n*sizeof(bool));
	memset(marked,0,n*sizeof(bool));
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			if(!marked[i])
			{
				long long node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{
					ans=ans+rank[parent[graph[node][j]]]/outdeg[graph[node][j]];
				}
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
		}
		for(i=0;i<n;i++)
			if(!marked[i]) rank[mapit[i]]=curr[i];
		iterations++;
		if(iterations%20==0)
		{
			for(i=0;i<n;i++)
			{
				if(!marked[i])
				{
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}
			}
		}
	}while(error > thres );
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
}

float computeparallel(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial, long long nn)
{
	float total = 0.0;
	double damp=0.85;
	double thres=1e-10;
	long long i;
	double *curr = (double *)malloc(n*sizeof(double));
	double error=0;
	long long iterations=0;
	double randomp=(1-damp)/graph.size();

	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}

	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}

	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long pivot=0;
	long long thresh=10000;
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}

	long long *cn, *cm, *cmem, *csize, *coutdeg, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);

	do
	{
		error=0;
		
		for(i=0;i<n;i++){
			curr[i]=0;
		}

		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel3test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);

		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
		for(i=pivot;i<n;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel3test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		}  

		double anse=0;
		for(i=0;i<n;i++){
			anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		}

		for(i=0;i<n;i++){
			rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
		}

		iterations++;
		error = anse;
	}while(error > thres );
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(crank);
	return total;
}


void computerank(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial)
{
	double damp=0.85;
	double thres=1e-10;
	long long i,j;
	vector < double > curr(n);
	double error=0;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			{
				long long node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
		}
		for(i=0;i<n;i++)
		{
			{
				rank[mapit[i]]=curr[i];
			}
		}
	}while(error > thres );
	curr.clear();
}

float computeparalleld(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial, long long nn)
{
	float total = 0.0;
	double thres=1e-10;
	double dis=1e-12;
	double value=((dis)*10.0)/n;
	long long i;
	double *curr = (double *)malloc(n*sizeof(double));
	double *prev = (double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	long long *marked = (long long *)malloc(n*sizeof(long long));
	memset(marked,0,n*sizeof(long long));
	double error=0;
	long long iterations=0;
	double damp = 0.85;
	double randomp=(1-damp)/graph.size();

	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long thresh=10000;
	long long pivot=0;
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}

	long long *cn, *cm, *cmem, *csize, *coutdeg, *ctemp, *cgraph, *cmarked;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&cmarked, n*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);

	do
	{
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(long long), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel4test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<n;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel4test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		}  

		double anse=0;
		for(i=0;i<n;i++)
			if(!marked[i])
				anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		iterations++;
		for(i=0;i<n;i++)
		{
			if(!marked[i])   
			{
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}   
		}
		if(iterations%20==0)
		{   
			for(i=0;i<n;i++)
			{   
				if(!marked[i])
				{   

					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else
						prev[i]=curr[i];
				}   
			}   
		}   
		error = anse;
	}while(error > thres );
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(cmarked);
	cudaFree(crank);
	return total;
}

void computerankd(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial)
{
	double damp = 0.85;
	double thres=1e-10;
	long long i,j;
	vector < double > curr(n);
	vector < double > prev(n);
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	double dis=1e-12;
	double value=(dis*10.0)/n;
	// double bound=1e-5;
	bool *marked = (bool *)malloc(n*sizeof(bool));
	memset(marked,0,n*sizeof(bool));
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			if(!marked[i])
			{
				long long node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
		}
		for(i=0;i<n;i++)
			if(!marked[i]) rank[mapit[i]]=curr[i];
		iterations++;
		if(iterations%20==0)
		{
			for(i=0;i<n;i++)
			{
				if(!marked[i])
				{
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}
			}
		}
	}while(error > thres );
	curr.clear();
}

float computeparallelc(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers, long long nn)
{
	float total = 0.0;
	double damp=0.85;
	double thres=1e-10;
	long long i;
	double *curr = (double *)malloc(n*sizeof(double));
	double error=0;
	long long iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	// 0 to limit are the nodes which are not part of chain
	// limit to n are the nodes which are part of chain (we have formula to calculate the pagerank of them)
	// limit to n are dealt with at the end
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(node==redir[node])
		{   
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	} 

	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long pivot=0;
	long long thresh = 10000;
	for(i=0;i<limit;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   

	long long *cn, *cm, *cmem, *csize, *coutdeg, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel3test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;
	
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<limit;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel3test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		}  

		double anse=0;
		for(i=0;i<limit;i++){
			anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		}

		for(i=0;i<limit;i++){
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
		}
		iterations++;
		error = anse; 
	}while(error > thres);
	// limit to n calculation
	for(i=limit;i<n;i++)
	{   
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(crank);
	return total;
}

void computerankc(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers)
{
	double damp=0.85;
	double thres=1e-10;
	long long i,j;
	vector < double > curr(n);
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{
		long long node=mapit[i];
		if(node==redir[node])
		{
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			long long node=mapit[i];
			double ans=0;
			for(j=0;j<graph[node].size();j++)
			{
				ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
			}
			curr[i]=randomp+damp*ans+initial[mapit[i]];
			error=max(error,fabs(curr[i]-rank[node]));
		}
		for(i=0;i<limit;i++)
			rank[mapit[i]]=curr[i];
		iterations++;
	}while(error > thres );
	for(i=limit;i<n;i++)
	{
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}
}

float computeparalleldc(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers, long long nn)
{
	float total = 0.0;
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double(n);
	long long i;
	double *curr = (double *)malloc(n*sizeof(double));
	double *prev = (double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++){
		prev[i]=1.0/n;
	}
	long long *marked = (long long *)malloc(n*sizeof(long long));
	memset(marked,0,n*sizeof(long long));
	double error=0;
	long long iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(node==redir[node])
		{   
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	}

	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long pivot=0;
	long long thresh=10000;
	for(i=0;i<limit;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   

	long long *cn, *cm, *cmem, *csize, *coutdeg, *ctemp, *cgraph, *cmarked;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&cmarked, n*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(long long), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel4test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;

		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<limit;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel4test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		} 

		double anse=0;
		for(i=0;i<limit;i++)
			if(!marked[i])
				anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		iterations++;
		for(i=0;i<limit;i++){
			if(!marked[i]){
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}
		}
		if(iterations%20==0)
		{   
			for(i=0;i<limit;i++)
			{   
				if(!marked[i])
				{   
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else
						prev[i]=curr[i];
				}   
			}   
		}
		error = anse;
	}while(error > thres);
	for(i=limit;i<n;i++)
	{   
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(cmarked);
	cudaFree(crank);
	return total;
}

void computerankdc(vector < vector < long long > > & graph,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers)
{
	double damp=0.85;
	double thres=1e-10;
	long long i, j;
	vector < double > curr(n);
	vector < double > prev(n,1.0/n);
	double value=((1e-12)*10.0)/double ( n );
	bool *marked = (bool *)malloc(n*sizeof(bool));
	memset(marked,0,n*sizeof(bool));
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{
		long long node=mapit[i];
		if(node==redir[node])
		{
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			long long node=mapit[i];
			if(!marked[i])
			{
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				}
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
		}
		for(i=0;i<limit;i++) if(!marked[i])
			rank[mapit[i]]=curr[i];
		iterations++;
		if(iterations%20==0){
			for(i=0;i<limit;i++)
			{
				if(!marked[i])
				{
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}
			}
		}
	}while(error > thres );
	for(i=limit;i<n;i++)
	{
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}
}

float computeparallelic(vector < vector < long long > > & graph,long long *parent,vector <long long > & left, long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers, long long nn)
{
	float total = 0.0;
	double damp=0.85;
	double thres=1e-10;
	long long i;
	double *curr = (double *)malloc(n*sizeof(double));
	double error=0;
	long long iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(node==redir[node])
		{   
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	}
	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++)
	{
		if(i)
		{
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}
		else
		{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long pivot=0;
	long long thresh=10000;
	for(i=0;i<limit;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   

	long long *cn, *cm, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&cparent, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));
	cudaMalloc((void**)&cparent, nn*sizeof(long long));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(long long), cudaMemcpyHostToDevice);
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel1test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;

		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<limit;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel1test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		} 

		double anse=0;
		for(i=0;i<limit;i++){
			anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		}

		for(i=0;i<limit;i++)
		{   
			{
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}   
		}
		iterations++;
		error = anse;
	}while(error > thres );
	for(i=limit;i<n;i++)
	{   
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(cparent);
	cudaFree(crank);
	return total;
}

void computerankic(vector < vector < long long > > & graph,long long *parent,vector < long long > & left,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers)
{
	double damp=0.85;
	double thres=1e-10;
	long long i, j;
	vector < double > curr(n);
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{
		long long node=mapit[i];
		if(node==redir[node])
		{
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			long long node=mapit[i];
			double ans=0;
			for(j=0;j<graph[node].size();j++)
			{
				ans=ans+rank[parent[graph[node][j]]]/outdeg[graph[node][j]];
			}
			curr[i]=randomp+damp*ans+initial[mapit[i]];
			error=max(error,fabs(curr[i]-rank[node]));
		}
		for(i=0;i<limit;i++)
			rank[mapit[i]]=curr[i];
		iterations++;
	}while(error > thres );
	for(i=limit;i<n;i++)
	{
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
}

float computeparallelidc(vector < vector < long long > > & graph, long long *parent,vector <long long> & left,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers, long long nn)
{
	float total = 0.0;
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	long long i;
	double *curr = (double *)malloc(n*sizeof(double));
	double *prev = (double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++){
		prev[i]=1.0/n;
	}
	long long *marked = (long long *)malloc(n*sizeof(long long));
	memset(marked,0,n*sizeof(long long));
	double error=0;
	long long iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{   
		long long node=mapit[i];
		if(node==redir[node])
		{   
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	}

	long long *mem = (long long *)malloc(n*sizeof(long long));
	long long *sz = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	long long *temp = (long long *)malloc(n*sizeof(long long));
	long long szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	long long *graphh = (long long *)malloc(szz*sizeof(long long));
	long long k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	long long pivot=0;
	long long thresh=10000;
	for(i=0;i<limit;i++)
	{   
		long long node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			long long temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   

	long long *cn, *cm, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph, *cmarked;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(long long));
	cudaMalloc((void**)&cm, sizeof(long long));
	cudaMalloc((void**)&cmem, n*sizeof(long long));
	cudaMalloc((void**)&csize, n*sizeof(long long));
	cudaMalloc((void**)&coutdeg, nn*sizeof(long long));
	cudaMalloc((void**)&cparent, nn*sizeof(long long));
	cudaMalloc((void**)&ctemp, n*sizeof(long long));
	cudaMalloc((void**)&cgraph, szz*sizeof(long long));
	cudaMalloc((void**)&cmarked, n*sizeof(long long));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));
	cudaMalloc((void**)&cparent, nn*sizeof(long long));

	cudaMemcpy(cn, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(long long), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(long long), cudaMemcpyHostToDevice);

		dim3 threadB(32,32);
		dim3 blockB(32,32);

		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);

		cudaEventRecord(start, 0);

		kernel2test<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);
		
		cudaDeviceSynchronize();

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
		total += elapsedTime;
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		for(i=pivot;i<limit;i++)
		{   
			{   
				cudaMemcpy(cm, &i, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kernel2test1<<<32,32>>>(cm, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);
				
				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);
			}   
		} 

		double anse=0;
		for(i=0;i<limit;i++)
			if(!marked[i])
				anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		iterations++;
		for(i=0;i<limit;i++)
		{
			if(!marked[i])   
			{
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}   
		}
		if(iterations%20==0)
		{   
			for(i=0;i<limit;i++)
			{   
				if(!marked[i])
				{   

					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else
						prev[i]=curr[i];
				}   
			}   
		}   
		error = anse;
	}while(error > thres);
	for(i=limit;i<n;i++)
	{   
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	cudaFree(cn);
	cudaFree(cm);
	cudaFree(cmem);
	cudaFree(csize);
	cudaFree(coutdeg);
	cudaFree(cparent);
	cudaFree(ctemp);
	cudaFree(cgraph);
	cudaFree(ccurr);
	cudaFree(crank);
	cudaFree(cmarked);
	return total;
}

void computerankidc(vector < vector < long long > > & graph,long long *parent,vector < long long > & left,long long n,long long *outdeg,vector < long long > &  mapit,double *rank,double *initial,long long *level,long long *redir,double *powers)
{
	double damp=0.85;
	double thres=1e-10;
	long long i, j;
	vector < double > curr(n);
	vector < double > prev(n,1.0/n);
	double value=((1e-12)*10.0)/double ( n );
	bool *marked = (bool *)malloc(n*sizeof(bool));
	memset(marked,0,n*sizeof(bool));
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	long long limit=0;
	for(i=0;i<n;i++)
	{
		long long node=mapit[i];
		if(node==redir[node])
		{
			long long temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			long long node=mapit[i];
			if(!marked[i])
			{
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{
					ans=ans+rank[parent[graph[node][j]]]/outdeg[graph[node][j]];
				}
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
		}
		for(i=0;i<limit;i++) if(!marked[i])
			rank[mapit[i]]=curr[i];
		iterations++;
		if(iterations%20==0){
			for(i=0;i<limit;i++)
			{
				if(!marked[i])
				{
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}
			}
		}
	}while(error > thres );
	for(i=limit;i<n;i++)
	{
		long long node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
}
