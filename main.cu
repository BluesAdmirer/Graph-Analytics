#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>
#include "main.h"

using namespace std;

int inc=0;

int find(int u,int parent[]){
	if(parent[u]<0) return u;
	return parent[u]=find(parent[u],parent);
}

int unionit(int u,int v,int parent[]){
	int pu=find(u,parent);
	int pv=find(v,parent);
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

void dfs(vector < vector < int > > & graph,int visit[],int nvisit[],int node){
	nvisit[node]=1;
	for(int i=0;i<graph[node].size();i++)
		if(nvisit[graph[node][i]]==-1)
			dfs(graph,visit,nvisit,graph[node][i]);
	visit[inc++]=node;
}

void rdfs(vector < vector < int > > & graph,int nvisit[],int node,int component[],int com){
	nvisit[node]=1;
	component[node]=com;
	for(int i=0;i<graph[node].size();i++)
		if(nvisit[graph[node][i]]==-1)
			rdfs(graph,nvisit,graph[node][i],component,com);
}

void topobfs(vector < vector < int > > & graph, int order[], int visit[]){
	int i,j;
	queue < int > line;
	memset(visit, -1, sizeof(int)*graph.size());
	int indegree[graph.size()];
	memset(indegree,0,sizeof(indegree));
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

int optchain=0, optdead=0, optident=0;

int main(){

	int w;

	ifstream fin;
	fin.open("input.txt");
	
	ofstream fout;
	fout.open("output.txt");
	
	int n,m;
	fin >> n >> m;
	
	int i,j,u,v;
	vector < vector < int > > graph(n), rgraph(n), rcgraph(n), rcwgraph(n);

	int outdeg[n];
	memset(outdeg,0,sizeof(outdeg));

	for(i=0;i<m;i++){
		fin >> u >> v, --u,--v;
		graph[u].push_back(v);
		rgraph[v].push_back(u);
		outdeg[u]++;
	}

	int visit[n];
	memset(visit, -1, sizeof(visit));

	int component[n];
	memset(component, -1, sizeof(component));
	
	int nvisit[n];
	memset(nvisit, -1, sizeof(nvisit));
	
	for(i=0;i<n;i++){
		if(nvisit[i]==-1) {
			dfs(graph,visit,nvisit,i);	
		}
	}

	memset(nvisit,-1,sizeof(nvisit));
	
	int com=0;
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
	
	vector < vector < int > > members(com), compgr(com);
	
	for(i=0;i<n;i++){
		for(j=0;j<graph[i].size();j++){
			if(component[i]!=component[graph[i][j]]){
				compgr[component[i]].push_back(component[graph[i][j]]);
			}
		}
	}

	int order[com];
	memset(nvisit,0,sizeof(nvisit));
	
	inc=0;
	topobfs(compgr,order,nvisit);
	
	int number[n];
	memset(number,0,sizeof(number));

	for(i=0;i<n;i++){
		if(rgraph[i].size()==1){
			number[rgraph[i][0]]++;
		}
	}

	int equiperc=0;
	for(i=0;i<n;i++){
		equiperc=equiperc+max(0,number[i]-1);
	}
	
	double vai=double(equiperc)/n;
	double ratio=double(m)/n;

	if(vai>0.06 && ratio>3.0)
		optident=1;
	
	int parent2[n];
	memset(parent2,-1,sizeof(parent2));
	
	int parent1[n];
	memset(parent1,-1,sizeof(parent1));
	
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
	
	int redir[n];
	for(i=0;i<n;i++){
		redir[i]=i;
	}

	int levelz[n];
	memset(levelz,0,sizeof(levelz));
	
	double randomp=0.15/n;

	double powers[n];
	powers[0]=1;
	for(i=1;i<n;i++){
		powers[i]=powers[i-1]*0.85;
	}
	
	int vac=0;
	
	for(i=0;i<n;i++)
	{
		if(rgraph[i].size()>1 || graph[i].size()>1 ) continue;
		if(parent2[i]!=-1) continue;
		int node=i;
		int iterations=0;
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

	// optdead =  1;

	if(optident==1 && optchain==0 && optdead==0){
		int parent[n];
		vector < vector < int > > left(com);
		for(i=0;i<n;i++){
			parent[i]=i;
		}
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
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
			if(parent[i]==i) {
				members[component[i]].push_back(i);
			}
			else{
				left[component[i]].push_back(i);
				noo++;
			}
		}

		double rank[n];
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		int thresh=100000;

		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<par.size()-1;i++){
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++){
				int sum=0;
				for(j=0;j<members[order[w]].size();j++){
					sum=sum+rgraph[members[order[w]][j]].size();
				}
				if(sum>thresh){
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++){
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					}
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++){
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					}
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
				long long val=computeparalleli(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++){
				long long val=computeranki(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
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
		for(i=0;i<n;i++){
			parent[i]=i;
		}
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
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
		for(int k=0;k<n;k++){
			for(i=0;i<hvalues[k].size();i++){
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++){
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();
		int noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i){
				members[component[i]].push_back(i);
			}
			else{
				left[component[i]].push_back(i);
				noo++;
			}
		}
		
		double rank[n];
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}

		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		int thresh=100000;
		
		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<par.size()-1;i++){
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++){
				int sum=0;
				for(j=0;j<members[order[w]].size();j++){
					sum=sum+rgraph[members[order[w]][j]].size();
				}
				if(sum>thresh){
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
				long long val=computeparallelid(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
				long long val=computerankid(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
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
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}

		vector < int > par;
		par.push_back(0);

		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		int thresh=100000;
		
		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++){
					sum=sum+rgraph[members[order[w]][j]].size();
				}
				if(sum>thresh){
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
			{
				long long val=computeparallel(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerank(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
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
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}

		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh){
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
			{
				long long val=computeparalleld(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankd(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
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
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}

		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++){
					sum=sum+rgraph[members[order[w]][j]].size();
				}
				if(sum>thresh)
				{
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
			{
				long long val=computeparallelc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
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
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}

		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		for(i=0;i<n;i++){
			members[component[i]].push_back(i);
		}

		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++){
					sum=sum+rgraph[members[order[w]][j]].size();
				}
				if(sum>thresh)
				{
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
			{
				long long val=computeparalleldc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankdc(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
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

		for(i=0;i<n;i++){
			parent[i]=i;
		}

		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);

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

		for(int k=0;k<n;k++){
			for(i=0;i<hvalues[k].size();i++){
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++){
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
				members[component[i]].push_back(i);
			}
			else
			{
				left[component[i]].push_back(i);
				noo++;
			}
		}
		
		double rank[n];
		for(i=0;i<n;i++){
			rank[i]=1.0/n;
		}

		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}

		int thresh=100000;
		
		double initial[n];
		memset(initial,0,sizeof(initial));
		
		for(i=0;i<par.size()-1;i++){
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++){
					sum=sum+rgraph[members[order[w]][j]].size();
				}
				if(sum>thresh)
				{
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
			{
				long long val=computeparallelic(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankic(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
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
		for(i=0;i<n;i++){
			parent[i]=i;
		}
		vector < vector <  pair  <  pair < long long , int > , int >  > > hvalues(n);
		
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
		
		for(int k=0;k<n;k++){
			for(i=0;i<hvalues[k].size();i++){
				for(j=i;j<hvalues[k].size() && hvalues[k][j].first==hvalues[k][i].first ;j++){
					parent[hvalues[k][j].second]=hvalues[k][i].second;
				}
				i=j-1;
			}
		}
		hvalues.clear();

		int noo=0;
		for(i=0;i<n;i++){
			if(parent[i]==i){
				members[component[i]].push_back(i);
			}
			else{
				left[component[i]].push_back(i);
				noo++;
			}
		}
		
		double rank[n];
		for(i=0;i<n;i++) rank[i]=1.0/n;
		vector < int > par;
		par.push_back(0);
		for(i=0;i<com;i++){
			int j=i;
			while(j<com && nvisit[order[j]]==nvisit[order[i]]){
				j++;
			}
			par.push_back(j);
			i=j-1;
		}
	
		int thresh=100000;
	
		double initial[n];
		memset(initial,0,sizeof(initial));

		for(i=0;i<par.size()-1;i++){
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(sum>thresh)
				{
					int temp=order[pivot];
					order[pivot]=order[w];
					order[w]=temp;
					pivot++;
				}
			}
			int k;
			for(w=par[i];w<pivot;w++)
			{
#pragma omp parallel for private(j,k)
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
#pragma omp parallel for private(w,j,k)
			for(w=pivot;w<par[i+1];w++)
			{
				for(j=0;j<members[order[w]].size();j++)
				{
					int node=members[order[w]][j];
					for(k=0;k<rcwgraph[node].size();k++)
						initial[node]+=rank[rcwgraph[node][k]]/outdeg[rcwgraph[node][k]];
					initial[node]=0.85*initial[node];
				}
			}
			for(j=par[i];j<pivot;j++)
			{
				long long val=computeparallelic(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankic(rcgraph,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
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

