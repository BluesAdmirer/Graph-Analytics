#include<cstdio>
#include<iostream>
#include<vector>
#include<cstring>
#include<queue>
#include<algorithm>
#include<fstream>
#include "Final_Pagerank_cuda.cuh"

int inc=0;
using namespace std;

int find(int u,int parent[])
{
	if(parent[u]<0) return u;
	return parent[u]=find(parent[u],parent);
}

int unionit(int u,int v,int parent[])
{
	int pu=find(u,parent);
	int pv=find(v,parent);
	if(pu==pv) return 0;
	if(-parent[pu]>-parent[pv])
	{   
		parent[pu]=parent[pu]+parent[pv];
		parent[pv]=pu;
	}   
	else
	{   
		parent[pv]=parent[pu]+parent[pv];
		parent[pu]=pv;
	}   
	return 1;
}

void dfs(vector < vector < int > > & graph,int visit[],int nvisit[],int node)
{
	nvisit[node]=1;
	for(int i=0;i<graph[node].size();i++)
		if(nvisit[graph[node][i]]==-1)
			dfs(graph,visit,nvisit,graph[node][i]);
	visit[inc++]=node;
}


void rdfs(vector < vector < int > > & graph,int nvisit[],int node,int no[],int com)
{
	nvisit[node]=1;
	no[node]=com;
	for(int i=0;i<graph[node].size();i++)
		if(nvisit[graph[node][i]]==-1)
			rdfs(graph,nvisit,graph[node][i],no,com);
}


void toposort(int node,int order[],vector < vector < int > > & graph,int visit[])
{
	int i;
	visit[node]=1;
	for(i=0;i<graph[node].size();i++)
		if(!visit[graph[node][i]]) toposort(graph[node][i],order,graph,visit);
	order[inc--]=node;
}

void topobfs(vector < vector < int > > & graph,int order[],int visit[])
{
	int i,j;
	queue < int > line;
	memset(visit,-1,graph.size()*sizeof(int));
	int indegree[graph.size()];
	memset(indegree,0,sizeof(indegree));
	for(i=0;i<graph.size();i++)
		for(j=0;j<graph[i].size();j++)
			indegree[graph[i][j]]++;
	for(i=0;i<graph.size();i++)
	{
		if(!indegree[i]) {
			line.push(i);
			visit[i]=0;
			indegree[i]--;
		}
	}
	while(!line.empty())
	{
		int node=line.front();
		line.pop();
		order[inc++]=node;
		for(i=0;i<graph[node].size();i++)
		{
			indegree[graph[node][i]]--;
			if(indegree[graph[node][i]]==0)
			{
				line.push(graph[node][i]);
				visit[graph[node][i]]=visit[node]+1;
			}
		}
	}
}

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


int optchain=0,optdead=0,optident=0;

int main()
{
	time_t first = clock();
	ifstream fin;
	fin.open("test.txt");
	ofstream fout;
	fout.open("testcu.txt");
	int n,m;
	fin >> n >> m;
	int i,j;
	vector < vector < int > > graph(n);
	vector < vector < int > > rgraph(n);
	vector < vector < int > > rcgraph(n);
	vector < vector < int > > rcwgraph(n);
	int u,v;
	int outdeg[n];
	memset(outdeg,0,sizeof(outdeg));
	for(i=0;i<m;i++)
	{
		fin >> u >> v;
		u--;v--;
		graph[u].push_back(v);
		rgraph[v].push_back(u);
		outdeg[u]++;
	}
	int visit[n];
	int no[n];
	memset(no,-1,sizeof(no));
	memset(visit,-1,sizeof(visit));
	int nvisit[n];
	memset(nvisit,-1,sizeof(nvisit));
	for(i=0;i<n;i++)
		if(nvisit[i]==-1) {
			dfs(graph,visit,nvisit,i);	
		}
	memset(nvisit,-1,sizeof(nvisit));
	int com=0;
	for(i=n-1;i>=0;i--)
	{
		if(nvisit[visit[i]]==-1)
		{
			rdfs(rgraph,nvisit,visit[i],no,com);
			com++;
		}
	}
	for(i=0;i<n;i++)
		for(j=0;j<rgraph[i].size();j++)
			if(no[i]==no[rgraph[i][j]]) rcgraph[i].push_back(rgraph[i][j]);
			else rcwgraph[i].push_back(rgraph[i][j]);
	vector < vector < int > > members(com);
	vector < vector < int >  > compgr(com);
	for(i=0;i<n;i++)
		for(j=0;j<graph[i].size();j++)
			if(no[i]!=no[graph[i][j]]) compgr[no[i]].push_back(no[graph[i][j]]);
	int order[com];
	memset(nvisit,0,sizeof(nvisit));
	inc=0;
	
	topobfs(compgr,order,nvisit);

	int number[n];
	memset(number,0,sizeof(number));
	for(i=0;i<n;i++) if(rgraph[i].size()==1) number[rgraph[i][0]]++;
	int equiperc=0;
	for(i=0;i<n;i++) equiperc=equiperc+max(0,number[i]-1);
	double vai=double(equiperc)/n;
	double ratio=double(m)/n;
	if(vai>0.06 && ratio>3.0)
		optident=1;

	int parent2[n];
	memset(parent2,-1,sizeof(parent2));
	int parent1[n];
	memset(parent1,-1,sizeof(parent1));
	for(i=0;i<n;i++)
	{
		if(rgraph[i].size()>1 || graph[i].size()>1 ) continue;
		for(j=0;j<rcgraph[i].size();j++)
		{
			if(graph[rcgraph[i][j]].size()>1 || rgraph[rcgraph[i][j]].size()>1) continue;
			if(unionit(rcgraph[i][j],i,parent1))
				parent2[i]=rcgraph[i][j];
		}
	}
	int redir[n];
	int levelz[n];
	memset(levelz,0,sizeof(levelz));
	for(i=0;i<n;i++) redir[i]=i;
	double powers[n];
	powers[0]=1;
	for(i=1;i<n;i++)
		powers[i]=powers[i-1]*0.85;
	int vac=0;
	int temp=0;
	for(i=0;i<n;i++)
	{
		if(rgraph[i].size()>1 || graph[i].size()>1 ) continue;
		if(parent2[i]!=-1) continue;
		int node=i;
		int iterations=0;
		while(graph[node].size())
		{
			node=graph[node][0];
			if(no[node]!=no[i] || node==i || graph[node].size()>1 || rgraph[node].size()>1) break;
			iterations++;
			redir[node]=i;
			levelz[node]=iterations;
		}
		vac=vac+iterations;
		temp+=redir[i];
	}
	double rac=double(vac)/n;
	if(rac>0.2)
		optchain=1;

	if(optident==1 && optchain==0 && optdead==0)
	{
		int parent[n];
		vector < vector < int > > left(com);
		for(i=0;i<n;i++)
			parent[i]=i;
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
		for(i=0;i<n;i++)
		{
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) continue;
			if(rgraph[i].size()==1)
			{
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],no[i]),i));
			}
			else
			{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,no[i]),i));
			}
		}
		for(i=0;i<n;i++)
			sort(hvalues[i].begin(),hvalues[i].end());
		for(int k=0;k<n;k++)
		{
			for(i=0;i<hvalues[k].size();i++)
			{
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++)
				{
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();
		int noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		}
		
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		double initial[n];
		memset(initial,0,sizeof(initial));
		int *cn;
		double *cinitial;
		cudaMalloc((void**)&cn, sizeof(int));
		cudaMemcpy(cn, &n, sizeof(int), cudaMemcpyHostToDevice);
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);

			for(j=par[i];j<par[i+1];j++){
				long long val=computeparalleli(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,n );
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	if(optident==1 && optchain==0 && optdead==1)	
	{
		int parent[n];
		vector < vector < int > > left(com);
		for(i=0;i<n;i++)
			parent[i]=i;
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
		for(i=0;i<n;i++)
		{
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) continue;
			if(rgraph[i].size()==1)
			{
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],no[i]),i));
			}
			else
			{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,no[i]),i));
			}
		}
		for(i=0;i<n;i++)
			sort(hvalues[i].begin(),hvalues[i].end());
		for(int k=0;k<n;k++)
		{
			for(i=0;i<hvalues[k].size();i++)
			{
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++)
				{
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();
		int noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		}
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);

			for(j=par[i];j<par[i+1];j++){
				long long val=computeparallelid(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	if(optident==0 && optchain==0 && optdead==0)
	{
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<n;i++)
			members[no[i]].push_back(i);
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);

			for(j=par[i];j<par[i+1];j++){
				long long val=computeparallel(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	if(optident==0 && optchain==0 && optdead==1)
	{
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<n;i++)
			members[no[i]].push_back(i);
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);

			for(j=par[i];j<par[i+1];j++){
				long long val=computeparalleld(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	if(optident==0 && optchain==1 && optdead==0)
	{
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		for(i=0;i<n;i++)
			members[no[i]].push_back(i);
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			for(j=par[i];j<par[i+1];j++)
			{
				long long val=computeparallelc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers,n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}

	}
	if(optident==0 && optchain==1 && optdead==1)
	{

		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		for(i=0;i<n;i++)
			members[no[i]].push_back(i);
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			for(j=par[i];j<par[i+1];j++)
			{
				long long val=computeparalleldc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){ 
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	if(optident==1 && optchain==1 && optdead==0)
	{
		int parent[n];
		vector < vector < int > > left(com);
		for(i=0;i<n;i++)
			parent[i]=i;
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
		for(i=0;i<n;i++)
		{
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) continue;
			if(rgraph[i].size()==1)
			{
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],no[i]),i));
			}
			else
			{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,no[i]),i));
			}
		}
		for(i=0;i<n;i++)
			sort(hvalues[i].begin(),hvalues[i].end());
		for(int k=0;k<n;k++)
		{
			for(i=0;i<hvalues[k].size();i++)
			{
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++)
				{
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();
		int noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		}
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			
			for(j=par[i];j<par[i+1];j++)
			{
				long long val=computeparallelic(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	if(optident==1 && optchain==1 && optdead==1)
	{
		int parent[n];
		vector < vector < int > > left(com);
		for(i=0;i<n;i++)
			parent[i]=i;
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
		for(i=0;i<n;i++)
		{
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) continue;
			if(rgraph[i].size()==1)
			{
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],no[i]),i));
			}
			else
			{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,no[i]),i));
			}
		}
		for(i=0;i<n;i++)
			sort(hvalues[i].begin(),hvalues[i].end());
		for(int k=0;k<n;k++)
		{
			for(i=0;i<hvalues[k].size();i++)
			{
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++)
				{
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();
		int noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		}
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<par.size()-1;i++)
		{
			int *cstart;
			cudaMalloc((void**)&cstart, sizeof(int));
			cudaMemcpy(cstart, &par[i], sizeof(int), cudaMemcpyHostToDevice);

			int *cend;
			cudaMalloc((void**)&cend, sizeof(int));
			cudaMemcpy(cend, &par[i+1], sizeof(int), cudaMemcpyHostToDevice);

			int *corder;
			cudaMalloc((void**)&corder, com*sizeof(int));
			cudaMemcpy(corder, order, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmemsz, memsz[com];
			int szz=0;
			for(int i1=0;i1<com;i1++){
				memsz[i1]=members[order[i1]].size();
				szz+=members[order[i1]].size();
			}
			cudaMalloc((void**)&cmemsz, com*sizeof(int));
			cudaMemcpy(cmemsz, memsz, com*sizeof(int), cudaMemcpyHostToDevice);

			int temp[com];
			for(int i1=0;i1<com;i1++){
				if(i1) temp[i1]=temp[order[i1-1]]+memsz[order[i1-1]];
				else temp[i1]=0;
			}
			int kk=0;
			int mem[szz];
			for(int i1=0;i1<com;i1++){
				for(int c:members[order[i1]]){
					mem[kk++]=c;
				}
			}

			int tempg[n];
			for(int i1=0;i1<n;i1++){
				if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
				else tempg[i1]=0;
			}
			int szzz = tempg[n-1]+rcwgraph[n-1].size();
			int kkk=0;
			int edges[szzz];
			for(int i1=0;i1<n;i1++){
				for(int c:rcwgraph[i1]){
					edges[kkk++]=c;
				}
			}
			int *ctemp;
			cudaMalloc((void**)&ctemp, com*sizeof(int));
			cudaMemcpy(ctemp, temp, com*sizeof(int), cudaMemcpyHostToDevice);

			int *cmembers;
			cudaMalloc((void**)&cmembers, szz*sizeof(int));
			cudaMemcpy(cmembers, mem, szz*sizeof(int), cudaMemcpyHostToDevice);

			int *ctempg;
			cudaMalloc((void**)&ctempg, n*sizeof(int));
			cudaMemcpy(ctempg, tempg, n*sizeof(int), cudaMemcpyHostToDevice);

			int *cedges;
			cudaMalloc((void**)&cedges, szzz*sizeof(int));
			cudaMemcpy(cedges, edges, szzz*sizeof(int), cudaMemcpyHostToDevice);

			int rcw[n];
			for(int i1=0;i1<n;i1++){
				rcw[i1] = rcwgraph[i1].size();
			}
			int *crcw;
			cudaMalloc((void**)&crcw, n*sizeof(int));
			cudaMemcpy(crcw, rcw, n*sizeof(int), cudaMemcpyHostToDevice);

			double *cinitial;
			cudaMalloc((void**)&cinitial, n*sizeof(double));
			cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);

			double *crank;
			cudaMalloc((void**)&crank, n*sizeof(double));
			cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

			int *coutdeg;
			cudaMalloc((void**)&coutdeg, n*sizeof(int));
			cudaMemcpy(coutdeg, outdeg, n*sizeof(int), cudaMemcpyHostToDevice);

			dim3 threadB(10,10);
			kernel<<<10,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
							cedges, coutdeg, corder, ctemp, ctempg);


			cudaDeviceSynchronize();

			cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);

			for(j=par[i];j<par[i+1];j++)
			{
				long long val=computeparallelidc(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
		}
		double sum=0;
		for(i=0;i<n;i++){
			sum=sum+rank[i];
		}
		for(i=0;i<n;i++){
			rank[i]=rank[i]/sum;
		}
		for(i=0;i<n;i++){
			fout << rank[i] << "\n";
		}
	}
	return 0;
}

