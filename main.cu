#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include <set>
#include <map>
#include <chrono>
#include "functions.cuh"

using namespace std;
using namespace std::chrono;

double total = 0.0;

long long inc=0;

long long find(long long u,long long *parent){
	if(parent[u]<0) return u;
	return parent[u]=find(parent[u],parent);
}

long long unionit(long long u,long long v,long long *parent){
	long long pu=find(u,parent);
	long long pv=find(v,parent);
	if(pu==pv) return 0;
	if(-parent[pu]>-parent[pv]){   
		parent[pu]=parent[pu]+parent[pv];
		parent[pv]=pu;
	}   
	else{   
		parent[pv]=parent[pu]+parent[pv];
		parent[pu]=pv;
	}   
	return 1;
}

void dfs(vector < vector < long long > > & graph,long long *visit,long long *nvisit,long long node){
	nvisit[node]=1;
	for(long long i=0;i<graph[node].size();i++)
		if(nvisit[graph[node][i]]==-1)
			dfs(graph,visit,nvisit,graph[node][i]);
	visit[inc++]=node;
}

void rdfs(vector < vector < long long > > & graph,long long *nvisit,long long node,long long *component,long long com){
	nvisit[node]=1;
	component[node]=com;
	for(long long i=0;i<graph[node].size();i++)
		if(nvisit[graph[node][i]]==-1)
			rdfs(graph,nvisit,graph[node][i],component,com);
}

void topobfs(vector < vector < long long > > & graph, long long *order, long long *visit){
	long long i,j;
	queue < long long > line;
	memset(visit, -1, sizeof(long long)*graph.size());
	long long indegree[graph.size()];
	memset(indegree,0,graph.size()*sizeof(long long));
	for(i=0;i<graph.size();i++){
		for(j=0;j<graph[i].size();j++){
			indegree[graph[i][j]]++;
		}
	}
	for(i=0;i<graph.size();i++){
		if(!indegree[i]) {
			line.push(i);
			visit[i]=0;
			indegree[i]--;
		}
	}
	while(!line.empty()){
		long long node=line.front();
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

long long optchain=0, optdead=0, optident=0;

int main(){

    auto start = high_resolution_clock::now();

	ifstream fin;
	fin.open("input.txt");
	
	ofstream fout;
	fout.open("output.txt");
	
	long long n,m;
	fin >> n >> m;
	
	long long i,j;
	vector < vector < long long > > graph(n), rgraph(n), rcgraph(n), rcwgraph(n);

	long long *outdeg = (long long *)malloc(n*sizeof(long long));
	memset(outdeg,0,n*sizeof(long long));

	set<long long> s;
	vector<pair<long long,long long>> edgess;

	for(i=0;i<m;i++){
		long long u,v;
		fin >> u >> v;
		s.insert(u);
		s.insert(v);
		edgess.push_back(make_pair(u,v));
	}
	map<long long,long long> hash;
	long long cnt=0;
	for(auto k:s){
		hash[k]=cnt++;
	}

	cout << cnt << " " << n << "\n";

	for(i=0;i<m;i++){
		long long u=hash[edgess[i].first], v = hash[edgess[i].second];
		graph[u].push_back(v);
		rgraph[v].push_back(u);
		outdeg[u]++;
	}

	long long *visit = (long long *)malloc(n*sizeof(long long));
	memset(visit, -1, n*sizeof(long long));

	long long *component = (long long *)malloc(n*sizeof(long long));
	memset(component, -1, n*sizeof(long long));
	
	long long *nvisit = (long long *)malloc(n*sizeof(long long));
	memset(nvisit, -1, n*sizeof(long long));
	
	for(i=0;i<n;i++){
		if(nvisit[i]==-1) {
			dfs(graph,visit,nvisit,i);	
		}
	}

	memset(nvisit,-1,n*sizeof(long long));
	
	long long com=0;
	for(i=n-1;i>=0;i--){
		if(nvisit[visit[i]]==-1){
			rdfs(rgraph,nvisit,visit[i],component,com);
			com++;
		}
	}

	for(i=0;i<n;i++){
		for(j=0;j<rgraph[i].size();j++){
			if(component[i]==component[rgraph[i][j]]){ 
				rcgraph[i].push_back(rgraph[i][j]);
			}
			else{ 
				rcwgraph[i].push_back(rgraph[i][j]);
			}
		}
	}
	
	vector < vector < long long > > members(com), compgr(com);
	
	for(i=0;i<n;i++){
		for(j=0;j<graph[i].size();j++){
			if(component[i]!=component[graph[i][j]]){
				compgr[component[i]].push_back(component[graph[i][j]]);
			}
		}
	}

	long long *order = (long long *)malloc(com*sizeof(long long));
	memset(nvisit,0,n*sizeof(long long));
	
	inc=0;
	topobfs(compgr,order,nvisit);
	
	long long *number = (long long *)malloc(n*sizeof(long long));
	memset(number,0,n*sizeof(long long));

	for(i=0;i<n;i++){
		if(rgraph[i].size()==1){
			number[rgraph[i][0]]++;
		}
	}

	long long equiperc=0;
	for(i=0;i<n;i++){
		equiperc=equiperc+max((long long)0,number[i]-1);
	}
	
	double vai=double(equiperc)/n;
	double ratio=double(m)/n;

	if(vai>0.06 && ratio>3.0)
		optident=1;
	
	long long *parent2 = (long long *)malloc(n*sizeof(long long));
	memset(parent2,-1,n*sizeof(long long));
	
	long long *parent1 = (long long *)malloc(n*sizeof(long long));
	memset(parent1,-1,n*sizeof(long long));
	
	for(i=0;i<n;i++){
		if(rgraph[i].size()>1 || graph[i].size()>1 ) 
			continue;
		for(j=0;j<rcgraph[i].size();j++){
			if(graph[rcgraph[i][j]].size()>1 || rgraph[rcgraph[i][j]].size()>1) 
				continue;
			if(unionit(rcgraph[i][j],i,parent1)){
				parent2[i]=rcgraph[i][j];
			}
		}
	}
	
	long long *redir = (long long *)malloc(n*sizeof(long long));
	for(i=0;i<n;i++){
		redir[i]=i;
	}

	long long *levelz = (long long *)malloc(n*sizeof(long long));
	memset(levelz,0,n*sizeof(long long));

	double *powers = (double *)malloc(n*sizeof(double));
	powers[0]=1;
	for(i=1;i<n;i++){
		powers[i]=powers[i-1]*0.85;
	}
	
	long long vac=0;
	
	for(i=0;i<n;i++)
	{
		if(rgraph[i].size()>1 || graph[i].size()>1 ) continue;
		if(parent2[i]!=-1) continue;
		long long node=i;
		long long iterations=0;
		while(graph[node].size())
		{
			node=graph[node][0];
			if(component[node]!=component[i] || node==i || graph[node].size()>1 || rgraph[node].size()>1) break;
			iterations++;
			redir[node]=i;
			levelz[node]=iterations;
		}
		vac=vac+iterations;
	}

	double rac=double(vac)/n;
	if(rac>0.2)
		optchain=1;

	long long *tempg = (long long *)malloc(n*sizeof(long long));
	for(long long i1=0;i1<n;i1++){
		if(i1) tempg[i1]=tempg[i1-1]+rcwgraph[i1-1].size();
		else tempg[i1]=0;
	}
	long long szzz = tempg[n-1]+rcwgraph[n-1].size();
	long long kkk=0;
	long long *edges = (long long *)malloc(szzz*sizeof(long long));
	for(long long i1=0;i1<n;i1++){
		for(long long c:rcwgraph[i1]){
			edges[kkk++]=c;
		}
	}

	long long *rcw = (long long *)malloc(n*sizeof(long long));
	for(long long i1=0;i1<n;i1++){
		rcw[i1] = rcwgraph[i1].size();
	}

	long long *cstart, *cend, *corder, *cmemsz, *ctemp, *coutdeg,
					*cmembers, *ctempg, *cedges, *crcw;
	double *cinitial,*crank;

	cudaMalloc((void**)&cstart, sizeof(long long));
	cudaMalloc((void**)&cend, sizeof(long long));
	cudaMalloc((void**)&corder, com*sizeof(long long));
	cudaMalloc((void**)&cmemsz, com*sizeof(long long));
	cudaMalloc((void**)&ctemp, com*sizeof(long long));
	cudaMalloc((void**)&crcw, n*sizeof(long long));
	cudaMalloc((void**)&cinitial, n*sizeof(double));
	cudaMalloc((void**)&crank, n*sizeof(double));
	cudaMalloc((void**)&coutdeg, n*sizeof(long long));
	cudaMalloc((void**)&ctempg, n*sizeof(long long));
	cudaMalloc((void**)&cedges, szzz*sizeof(long long));

	cudaMemcpy(corder, order, com*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(ctempg, tempg, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(cedges, edges, szzz*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(crcw, rcw, n*sizeof(long long), cudaMemcpyHostToDevice);
	cudaMemcpy(coutdeg, outdeg, n*sizeof(long long), cudaMemcpyHostToDevice);

	double *rank = (double *)malloc(n*sizeof(double));
	for(i=0;i<n;i++){
		rank[i]=1.0/n;
	}

	cout << optident << " " << optchain << " " << optdead << "\n";

	if(optident==1 && optchain==0 && optdead==0)
	{
		long long *parent = (long long *)malloc(n*sizeof(long long));
		vector < vector < long long > > left(com);
		for(i=0;i<n;i++)
			parent[i]=i;
		vector < vector <  pair  <  pair < long long , long long > , long long >  > > hvalues(n);
		for(i=0;i<n;i++)
		{
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) continue;
			if(rgraph[i].size()==1)
			{
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],component[i]),i));
			}
			else
			{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%n ].push_back(make_pair(make_pair(val,component[i]),i));
			}
		}

		for(i=0;i<n;i++)
			sort(hvalues[i].begin(),hvalues[i].end());
		for(long long k=0;k<n;k++)
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
		long long noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i) 
			{
				members[component[i]].push_back(i);
			}
			else
			{
				left[component[i]].push_back(i);
				noo++;
			}
		}
		
		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++)
		{
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]])
				j++;
			par.push_back(j);
			i=j-1;
		}
		long long thresh=100000;
		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));
		long long w;

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}

		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}

		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));

		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}

		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);

		for(i=0;i<par.size()-1;i++)
		{
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cn, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1])
			{
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
			for(j=par[i];j<pivot;j++){
				total += computeparalleli(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computeranki(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial, n);
			}
		}
	}

	if(optident==1 && optchain==0 && optdead==1)	
	{
		long long *parent = (long long *)malloc(n*sizeof(long long));
		vector < vector < long long > > left(com);
		for(i=0;i<n;i++){
			parent[i]=i;
		}
		vector < vector <  pair  <  pair < long long , long long > , long long >  > > hvalues(n);
		for(i=0;i<n;i++){
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) 
				continue;
			if(rgraph[i].size()==1){
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],component[i]),i));
			}
			else{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,component[i]),i));
			}
		}
		for(i=0;i<n;i++){
			sort(hvalues[i].begin(),hvalues[i].end());
		}
		for(long long k=0;k<n;k++){
			for(i=0;i<hvalues[k].size();i++){
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++){
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();
		long long noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i){
				members[component[i]].push_back(i);
			}
			else{
				left[component[i]].push_back(i);
				noo++;
			}
		}

		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}
		
		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,sizeof(initial));

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}

		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);
		
		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cn, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
			for(j=par[i];j<pivot;j++){
				total += computeparallelid(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerankid(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial, n);
			}
		}
	}

	if(optident==0 && optchain==0 && optdead==0)
	{
		vector < long long > par;
		par.push_back(0);

		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}
		
		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));
		
		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;

		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}

		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);

		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);


				kerneltest1<<<blockB,threadB>>>(cn, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
			for(j=par[i];j<pivot;j++){
				total += computeparallel(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerank(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
		}
	}

	if(optident==0 && optchain==0 && optdead==1)
	{
		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));
		
		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}
		
		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);
		
		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
				
			for(j=par[i];j<pivot;j++){
				total += computeparalleld(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerankd(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
		}
	}

	if(optident==0 && optchain==1 && optdead==0)
	{
		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}

		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);
		
		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
				
			for(j=par[i];j<pivot;j++){
				total += computeparallelc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerankc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial, levelz, redir, powers);
			}
		}
	}

	if(optident==0 && optchain==1 && optdead==1)
	{
		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}
		
		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);
		
		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cn, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
			
			for(j=par[i];j<pivot;j++){
				total += computeparalleldc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerankdc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial, levelz, redir, powers);
			}
		}
	}
	
	if(optident==1 && optchain==1 && optdead==0)
	{
		long long parent[n];
		vector < vector < long long > > left(com);

		for(i=0;i<n;i++){
			parent[i]=i;
		}

		vector < vector <  pair  <  pair < long long , long long > , long long >  > > hvalues(n);

		for(i=0;i<n;i++){
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) 
				continue;
			if(rgraph[i].size()==1){
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],component[i]),i));
			}
			else{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,component[i]),i));
			}
		}

		for(i=0;i<n;i++){
			sort(hvalues[i].begin(),hvalues[i].end());
		}

		for(long long k=0;k<n;k++){
			for(i=0;i<hvalues[k].size();i++){
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++){
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}

		hvalues.clear();
		long long noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i) 
			{
				members[component[i]].push_back(i);
			}
			else
			{
				left[component[i]].push_back(i);
				noo++;
			}
		}

		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}
		
		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}

		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);
		
		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cn, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}

			for(j=par[i];j<pivot;j++){
				total += computeparallelic(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerankic(rcgraph, parent, left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial, levelz, redir, powers);
			}
		}
	}

	if(optident==1 && optchain==1 && optdead==1)
	{
		long long *parent = (long long *)malloc(n*sizeof(long long));
		vector < vector < long long > > left(com);
		for(i=0;i<n;i++){
			parent[i]=i;
		}
		vector < vector <  pair  <  pair < long long , long long > , long long >  > > hvalues(n);
		
		for(i=0;i<n;i++){
			if(rgraph[i].size()!=1 && rgraph[i].size()!=2) continue;
			if(rgraph[i].size()==1){
				hvalues[(rgraph[i][0])%n].push_back(make_pair(make_pair(rgraph[i][0],component[i]),i));
			}
			else{
				long long val=max(rgraph[i][1]+1,rgraph[i][0]+1)*(long long)(n+1)+min(rgraph[i][0]+1,rgraph[i][1]+1);
				hvalues[(val)%(long long)n ].push_back(make_pair(make_pair(val,component[i]),i));
			}
		}
	
		for(i=0;i<n;i++){
			sort(hvalues[i].begin(),hvalues[i].end());
		}
		
		for(long long k=0;k<n;k++){
			for(i=0;i<hvalues[k].size();i++){
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++){
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();

		long long noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i){
				members[component[i]].push_back(i);
			}
			else{
				left[component[i]].push_back(i);
				noo++;
			}
		}
	
		vector < long long > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			long long j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}
		
		double *initial = (double *)malloc(n*sizeof(double));
		memset(initial,0,n*sizeof(double));

		long long *memsz = (long long *)malloc(com*sizeof(long long));
		long long *temp = (long long *)malloc(com*sizeof(long long));
		long long szz=0;
		for(long long i1=0;i1<com;i1++){
			memsz[i1]=members[order[i1]].size();
			szz+=members[order[i1]].size();
		}
		for(long long i1=0;i1<com;i1++){
			if(i1) temp[i1]=temp[i1-1]+memsz[i1-1];
			else temp[i1]=0;
		}
		long long kk=0;
		long long *mem = (long long *)malloc(szz*sizeof(long long));
		for(long long i1=0;i1<com;i1++){
			for(long long c:members[order[i1]]){
				mem[kk++]=c;
			}
		}
		
		cudaMalloc((void**)&cmembers, szz*sizeof(long long));
		
		cudaMemcpy(cmemsz, memsz, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(ctemp, temp, com*sizeof(long long), cudaMemcpyHostToDevice);
		cudaMemcpy(cmembers, mem, szz*sizeof(long long), cudaMemcpyHostToDevice);
		
		long long w;
		long long thresh=100000;

		for(i=0;i<par.size()-1;i++){
			long long pivot=par[i];
			for(w=par[i];w<par[i+1];w++)
			{
				long long sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					long long temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			for(w=par[i];w<pivot;w++)
			{
				long long *cn;
				cudaMalloc((void**)&cn, sizeof(long long));
				cudaMemcpy(cn, &w, sizeof(long long), cudaMemcpyHostToDevice);

				dim3 threadB(32,32);
				dim3 blockB(32,32);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest1<<<blockB,threadB>>>(cn, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;

				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
				cudaFree(cn);
			}

			if(pivot < par[i+1]){
				cudaMemcpy(cstart, &pivot, sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cend, &par[i+1], sizeof(long long), cudaMemcpyHostToDevice);
				cudaMemcpy(cinitial, initial, n*sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(crank, rank, n*sizeof(double), cudaMemcpyHostToDevice);

				dim3 threadB(8,8,16);
				dim3 blockB(8,8,16);

				cudaEvent_t start, stop;
				cudaEventCreate(&start);
				cudaEventCreate(&stop);

				cudaEventRecord(start, 0);

				kerneltest<<<blockB,threadB>>>(cstart, cend, cmemsz, cmembers, crcw, cinitial, crank,
								cedges, coutdeg, corder, ctemp, ctempg);

				cudaDeviceSynchronize();

				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				float elapsedTime;
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cudaEventDestroy(start);
				cudaEventDestroy(stop);
				total += elapsedTime;
				cudaMemcpy(initial, cinitial, n*sizeof(double), cudaMemcpyDeviceToHost);
			}
	
			for(j=par[i];j<pivot;j++){
				total += computeparallelidc(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers, n);
			}
			for(j=pivot;j<par[i+1];j++){
				computerankidc(rcgraph, parent, left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial, levelz, redir, powers);
			}
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
	auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken: "
         << duration.count() / 1000000.0 << "seconds" << "\n";

    cout << "kernel time: " << total << "\n\n\n";
    cudaFree(cstart);
    cudaFree(cend);
    cudaFree(corder);
    cudaFree(cmemsz);
    cudaFree(ctemp);
    cudaFree(coutdeg);
    cudaFree(cmembers);
    cudaFree(ctempg);
    cudaFree(cedges);
    cudaFree(crcw);
    cudaFree(cinitial);
    cudaFree(crank);
	return 0;
}
