#include "kernels.cuh"
using namespace std;

int computeparalleli(vector<vector<int>> &graph, int parent[], vector<int> left, int n, int outdeg[], vector<int> &mapit, double rank[],double initial[], int nn)
{
	double damp=0.85;
	double thres=1e-10;
	int i;
	double curr[n];
	for(int i=0;i<n;i++){
		curr[i]=0;
	}
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cparent;
		cudaMalloc((void**)&cparent, nn*sizeof(int));
		cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel1<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);

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
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cparent;
		cudaMalloc((void**)&cparent, nn*sizeof(int));
		cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cmarked;
		cudaMalloc((void**)&cmarked, n*sizeof(int));
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel2<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);

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
	do
	{
		error=0;
		
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel3<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);

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
	do
	{
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cmarked;
		cudaMalloc((void**)&cmarked, n*sizeof(int));
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel4<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);

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
	// double value=((1e-12)*10.0)/double(n);
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
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel3<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg);

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
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cmarked;
		cudaMalloc((void**)&cmarked, n*sizeof(int));
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel4<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cmarked);

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
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cparent;
		cudaMalloc((void**)&cparent, nn*sizeof(int));
		cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);


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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel1<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent);

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
	do  
	{   
		error=0;
		for(i=0;i<n;i++){
			curr[i]=0;
		}
		int *cn;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &limit, sizeof(int), cudaMemcpyHostToDevice);

		int mem[n],sz[n];
		for(i=0;i<n;i++){
			mem[i]=mapit[i];
			sz[i]=graph[mapit[i]].size();
		}
		int *cmem;
		cudaMalloc((void**)&cmem, n*sizeof(int));
		cudaMemcpy(cmem, mem, n*sizeof(int), cudaMemcpyHostToDevice);

		int *csize;
		cudaMalloc((void**)&csize, n*sizeof(int));
		cudaMemcpy(csize, sz, n*sizeof(int), cudaMemcpyHostToDevice);

		double *ccurr;
		cudaMalloc((void**)&ccurr, n*sizeof(double));
		cudaMemcpy(ccurr, curr, n*sizeof(double), cudaMemcpyHostToDevice);

		double *crank;
		cudaMalloc((void**)&crank, nn*sizeof(double));
		cudaMemcpy(crank, rank, nn*sizeof(double), cudaMemcpyHostToDevice);

		int *coutdeg;
		cudaMalloc((void**)&coutdeg, nn*sizeof(int));
		cudaMemcpy(coutdeg, outdeg, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cparent;
		cudaMalloc((void**)&cparent, nn*sizeof(int));
		cudaMemcpy(cparent, parent, nn*sizeof(int), cudaMemcpyHostToDevice);

		int *cmarked;
		cudaMalloc((void**)&cmarked, n*sizeof(int));
		cudaMemcpy(cmarked, marked, n*sizeof(int), cudaMemcpyHostToDevice);

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

		int *ctemp;
		cudaMalloc((void**)&ctemp, n*sizeof(int));
		cudaMemcpy(ctemp, temp, n*sizeof(int), cudaMemcpyHostToDevice);

		int *cgraph;

		cudaMalloc((void**)&cgraph, szz*sizeof(int));
		cudaMemcpy(cgraph, graphh, szz*sizeof(int), cudaMemcpyHostToDevice);

		kernel2<<<10,10>>>(cn, csize, cmem, cgraph, ctemp, ccurr, crank, coutdeg, cparent, cmarked);

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

long long int computerank(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	// double value=1e-12/n;
	//bool marked[n];
	//memset(marked,0,sizeof(marked));
	double error=0;
	// time_t start,end;
	// double time_diff;
	//start=clock();
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	int sumiterations=0;
	do
	{
		error=0;
		//#pragma omp parallel for private(i,j) 
		for(i=0;i<n;i++)
		{
			//if(marked[i]<countit)	
			{
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
				//operations[omp_get_thread_num()]+=graph[node].size();
			}
		}
		for(i=0;i<n;i++)
		{
			//if(marked[i]<countit)
			{
				//	if(fabs(rank[mapit[i]]-curr[i]) < value )marked[i]++;
				//	else marked[i]=marked[i]/2;
				//operations[omp_get_thread_num()]+=2;
				rank[mapit[i]]=curr[i];
			}
		}
		iterations++;
	}while(error > thres );
	sumiterations=iterations*n;
	return iterations;
	curr.clear();
	return sumiterations;
}


long long int computerankd(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n);
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	double dis=1e-12;
	double value=(dis*10.0)/n;
	// double bound=1e-5;
	bool  marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	long long  iterations=0;
	double damp = 0.85;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			if(!marked[i])
			{
				int node=mapit[i];
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
	return iterations;
}

long long int computerankc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double error=0;
	long long  iterations=0;
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
	vector < int > spare;
	for(i=0;i<limit;i++)
	{
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare.push_back(graph[node][j]);
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			int node=mapit[i];
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
		for(i=0;i<spare.size();i++)
		{
			double val=powers[level[spare[i]]];
			rank[spare[i]]=rank[redir[spare[i]]]*val+(1.0-val)/graph.size();
		}
	}while(error > thres );
	for(i=limit;i<n;i++)
	{
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}
	return iterations;
}




long long int computerankdc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n,1.0/n);
	double value=((1e-12)*10.0)/double ( n );
	bool marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	long long  iterations=0;
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
	vector < int > spare;
	for(i=0;i<limit;i++)
	{
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare.push_back(graph[node][j]);
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			int node=mapit[i];
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
		if(iterations%20==0)
			for(i=0;i<limit;i++)
			{
				if(!marked[i])
				{
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}
			}
		for(i=0;i<spare.size();i++)
		{
			double val=powers[level[spare[i]]];
			rank[spare[i]]=rank[redir[spare[i]]]*val+(1.0-val)/graph.size();
		}
	}while(error > thres );
	for(i=limit;i<n;i++)
	{
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}
	return iterations;
}

long long int computeranki(vector < vector < int > > & graph, int parent[],vector < int > left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double error=0;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			{
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++){
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

		iterations++;
	}while(error > thres );
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	return iterations;
}


long long int computerankid(vector < vector < int > > & graph,int parent[],vector < int > & left, int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n);
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	double dis=1e-12;
	double value=(dis*10.0)/n;
	bool  marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	long long  iterations=0;
	double damp = 0.85;
	double randomp=(1-damp)/graph.size();
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			if(!marked[i])
			{
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{
					ans=ans+rank[parent[graph[node][j]]]/outdeg[graph[node][j]];
				}
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				error=max(error,fabs(curr[i]-rank[node]));
			}
			else{
				printf("dead node:: %d\n", mapit[i]);
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
		for(i=0;i<n;i++){
			printf("%f ", rank[mapit[i]]);
		}
		printf("\n");
	}while(error > thres );
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	return iterations;
}

long long int computerankic(vector < vector < int > > & graph,int parent[],vector < int > & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double error=0;
	long long  iterations=0;
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
	vector < int > spare;
	for(i=0;i<limit;i++)
	{
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[parent[graph[node][j]]]!=parent[graph[node][j]]) spare.push_back(parent[graph[node][j]]);
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			int node=mapit[i];
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
		for(i=0;i<spare.size();i++)
		{
			double val=powers[level[spare[i]]];
			rank[spare[i]]=rank[redir[spare[i]]]*val+(1.0-val)/graph.size();
		}
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

long long int computerankidc(vector < vector < int > > & graph,int parent[],vector < int > & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n,1.0/n);
	double value=((1e-12)*10.0)/double ( n );
	bool marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	long long  iterations=0;
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
	vector < int > spare;
	for(i=0;i<limit;i++)
	{
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[parent[graph[node][j]]]!=parent[graph[node][j]]) spare.push_back(parent[graph[node][j]]);
	}
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			int node=mapit[i];
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
		if(iterations%20==0)
			for(i=0;i<limit;i++)
			{
				if(!marked[i])
				{
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}
			}
		for(i=0;i<spare.size();i++)
		{
			double val=powers[level[spare[i]]];
			rank[spare[i]]=rank[redir[spare[i]]]*val+(1.0-val)/graph.size();
		}
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
