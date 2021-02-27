#include "kernels1.cuh"
using namespace std;

int computeparalleli(vector<vector<int>> &graph, int parent[], vector<int> left, int n, int outdeg[], vector<int> &mapit, double rank[],double initial[], int nn)
{
	int i, iterations = 0;
	double damp=0.85, thres=1e-10, error = 0;
	double randomp=(1-damp)/graph.size();
	double curr[n];
	for(int i=0;i<n;i++){
		curr[i]=0;
	}
	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}
	
	int *cn, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&cparent, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));

	cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}

		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel1<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		double anse=0;
		for(i=0;i<n;i++){
			anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		}

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
	return iterations;
}

int  computeparallelid(vector < vector < int > > & graph,int parent[],vector < int > & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[], int nn)
{
	double thres=1e-10;
	double dis=1e-12;
	double value=((dis)*10.0)/n;
	int i;
	double prev[n], curr[n];
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	int marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	int iterations=0;
	double damp = 0.85;
	double randomp=(1-damp)/graph.size();

	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph, *cmarked;;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&cparent, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&cmarked, n*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);


		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel2<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		double anse=0;
		for(i=0;i<n;i++){
			if(!marked[i]){
				anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
			}
		}
		iterations++;
		for(i=0;i<n;i++){
			if(!marked[i])   {
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
			}   
		}
		if(iterations%20==0){   
			for(i=0;i<n;i++){   
				if(!marked[i]){   
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else
						prev[i]=curr[i];
				}   
			}   
		}   
		error = anse;
	}while(error > thres );
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	return iterations;
}

int  computeparallel(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[], int nn)
{
	double damp=0.85;
	double thres=1e-10;
	int i;
	double curr[n];
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();

	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

	do
	{
		error=0;
		
		for(i=0;i<n;i++){
			curr[i]=0;
		}

		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);


		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel3<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

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
	return iterations;
}

int  computeparalleld(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[], int nn)
{
	double thres=1e-10;
	double dis=1e-12;
	double value=((dis)*10.0)/n;
	int i;
	double curr[n], prev[n];
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	int marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	int iterations=0;
	double damp = 0.85;
	double randomp=(1-damp)/graph.size();

	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *ctemp, *cgraph, *cmarked;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&cmarked, n*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

	do
	{
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel4<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

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
	return iterations;
}

int computeparallelc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[], int nn)
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	double curr[n];
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	for(i=0;i<n;i++)
	{   
		int node=mapit[i];
		if(node==redir[node])
		{   
			int temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	} 
	vector<int> spare;
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare.push_back(graph[node][j]);
	}

	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		

		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel3<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);
	
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

		double anse=0;
		for(i=0;i<limit;i++){
			anse=max(anse, fabs(randomp+initial[mapit[i]]+damp*curr[i]-rank[mapit[i]]));
		}

		for(i=0;i<limit;i++){
				rank[mapit[i]]=damp*curr[i]+randomp+initial[mapit[i]];
		}
		iterations++;
		for(j=0;j<spare.size();j++)
		{   
			double val=powers[level[spare[j]]];
			rank[spare[j]]=rank[redir[spare[j]]]*val+(1.0-val)/graph.size();
		}
		error = anse; 
	}while(error > thres);
	for(i=limit;i<n;i++)
	{   
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	return iterations;
}

int computeparalleldc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[], int nn)
{
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	int i,j;
	double curr[n], prev[n];
	for(i=0;i<n;i++){
		prev[i]=1.0/n;
	}
	int marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	for(i=0;i<n;i++)
	{   
		int node=mapit[i];
		if(node==redir[node])
		{   
			int temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	}
	vector<int> spare(12);
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare.push_back(graph[node][j]);
	}

	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *ctemp, *cgraph, *cmarked;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&cmarked, n*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));

	cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel4<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);

		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

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
		for(j=0;j<spare.size();j++)
		{   
			double val=powers[level[spare[j]]];
			rank[spare[j]]=rank[redir[spare[j]]]*val+(1.0-val)/graph.size();
		}   
		error = anse;
	}while(error > thres);
	for(i=limit;i<n;i++)
	{   
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	return iterations;
}

int computeparallelic(vector < vector < int > > & graph,int parent[],vector <int > & left, int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[], int nn)
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	double curr[n];
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	for(i=0;i<n;i++)
	{   
		int node=mapit[i];
		if(node==redir[node])
		{   
			int temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	}
	vector<int> spare;
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare.push_back(graph[node][j]);
	}
	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
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
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&cparent, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));
	cudaMalloc((void**)&cparent, nn*sizeof(int));

	cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel1<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);

		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

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
		for(j=0;j<spare.size();j++)
		{   
			double val=powers[level[spare[j]]];
			rank[spare[j]]=rank[redir[spare[j]]]*val+(1.0-val)/graph.size();
		}
		error = anse;
	}while(error > thres );
	for(i=limit;i<n;i++)
	{   
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	return iterations;
}

int computeparallelidc(vector < vector < int > > & graph, int parent[],vector <int> & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[], int nn)
{
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	int i,j;
	double curr[n], prev[n];
	for(i=0;i<n;i++){
		prev[i]=1.0/n;
	}
	int marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	for(i=0;i<n;i++)
	{   
		int node=mapit[i];
		if(node==redir[node])
		{   
			int temp=mapit[limit];
			mapit[limit]=node;
			mapit[i]=temp;
			limit++;
		}   
	}
	vector<int> spare;
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare.push_back(graph[node][j]);
	}   

	int mem[n],sz[n];
	for(i=0;i<n;i++){
		mem[i]=mapit[i];
		sz[i]=graph[mapit[i]].size();
	}
	int temp[n];
	int szz=0;
	for(i=0;i<n;i++){
		if(i){
			temp[i]=temp[i-1]+graph[mapit[i-1]].size();
		}else{
			temp[i]=0;
		}
		szz+=graph[mapit[i]].size();
	}
	int graphh[szz];
	int k=0;
	for(i=0;i<n;i++){
		for(auto c:graph[mapit[i]]){
			graphh[k++]=c;
		}
	}

	int *cn, *cmem, *csize, *coutdeg, *cparent, *ctemp, *cgraph, *cmarked;
	double *ccurr, *crank;

	cudaMalloc((void**)&cn, sizeof(int));
	cudaMalloc((void**)&cmem, n*sizeof(int));
	cudaMalloc((void**)&csize, n*sizeof(int));
	cudaMalloc((void**)&coutdeg, nn*sizeof(int));
	cudaMalloc((void**)&cparent, nn*sizeof(int));
	cudaMalloc((void**)&ctemp, n*sizeof(int));
	cudaMalloc((void**)&cgraph, szz*sizeof(int));
	cudaMalloc((void**)&cmarked, n*sizeof(int));
	cudaMalloc((void**)&ccurr, n*sizeof(double));
	cudaMalloc((void**)&crank, nn*sizeof(double));
	cudaMalloc((void**)&cparent, nn*sizeof(int));

	cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);

	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

		dim3 threadB(1024,1024);
		dim3 blockB(63555,63535);
		kernel2<<<blockB,threadB>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);
		
		cudaMemcpy(curr, ccurr, n*sizeof(double), cudaMemcpyDeviceToHost);

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
		for(j=0;j<spare.size();j++)
		{   
			double val=powers[level[spare[j]]];
			rank[spare[j]]=rank[redir[spare[j]]]*val+(1.0-val)/graph.size();
		}      
	}while(error > thres);
	for(i=limit;i<n;i++)
	{   
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	return iterations;
}