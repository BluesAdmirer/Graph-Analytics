#include<cstdio>
#include<iostream>
#include<vector>
#include<cstring>
#include<queue>
#include<fstream>
#include<algorithm>
#include<omp.h>
#include<iostream>
#include<fstream>
#include "final_pagerank.h"

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

void topobfs(vector < vector < int > > & graph,int level[],int order[],int visit[])
{
	int i,j;
	queue < int > line;
	memset(visit, -1, sizeof(int)*graph.size());
	// memset(visit,-1,sizeof(visit));
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
		level[visit[node]]++;
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


int optchain=0;
int optdead=0;
int optident=0;




int main()
{
	ifstream fin;
	fin.open("test.txt");
	ofstream fout;
	fout.open("testcpp.txt");
	int n,m;
	// cin>>n>>m;
	fin >> n >> m;
	int i,j;
	vector < vector < int > > graph(n);
	vector < vector < int > > rgraph(n);
	vector < vector < int > > rcgraph(n);
	vector < vector < int > > rcwgraph(n);
	int u,v;
	int outdeg[n];
	int indeg[n];
	memset(outdeg,0,sizeof(outdeg));
	for(i=0;i<m;i++)
	{
		// cin>>u>>v;
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
	//time_t start,end;
	//double time_diff;
	//start=clock();
	time_t start,end,start1,end1;
	double time_diff,time1;
	start=clock();
	for(i=0;i<n;i++)
		if(nvisit[i]==-1) {
			dfs(graph,visit,nvisit,i);	
		}
	//cout<<"done"<<endl;
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
	end=clock();
	time_diff=((double)(end-start))/CLOCKS_PER_SEC;
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
	inc=com-1;
	int level[com];
	memset(level,0,sizeof(level));
	inc=0;
	start=clock();
	topobfs(compgr,level,order,nvisit);
	end=clock();
	time_diff=time_diff+((double)(end-start))/CLOCKS_PER_SEC;
	cout<<"After scc+topo= "<<time_diff<<endl;



	start=clock();
	int number[n];
	memset(number,0,sizeof(number));
	for(i=0;i<n;i++) if(rgraph[i].size()==1) number[rgraph[i][0]]++;
	int equiperc=0;
	for(i=0;i<n;i++) equiperc=equiperc+max(0,number[i]-1);
	double vai=double(equiperc)/n;
	double ratio=double(m)/n;
	cout<<"Percent of 1-degree ident nodes= "<<vai<<endl;
	if(vai>0.06 && ratio>3.0)
		optident=1;
	end=clock();
	time_diff=time_diff+((double)(end-start))/CLOCKS_PER_SEC;
	cout<<"After ident nodes= "<<time_diff<<endl;



	///chain preprocessing ///
	int parent2[n];
	memset(parent2,-1,sizeof(parent2));
	int parent1[n];
	memset(parent1,-1,sizeof(parent1));
	start1=clock();
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
	end1=clock();
	time1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;
	int redir[n];
	int levelz[n];
	memset(levelz,0,sizeof(levelz));
	for(i=0;i<n;i++) redir[i]=i;
	double randomp=0.15/n;
	double powers[n];
	powers[0]=1;
	for(i=1;i<n;i++)
		powers[i]=powers[i-1]*0.85;
	start1=clock();
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
			if(no[node]!=no[i] || node==i || graph[node].size()>1 || rgraph[node].size()>1) break;
			iterations++;
			redir[node]=i;
			levelz[node]=iterations;
		}
		vac=vac+iterations;
	}
	double rac=double(vac)/n;
	cout<<"percent of chain nodes= "<<rac<<endl;
	if(rac>0.2)
		optchain=1;
	end1=clock();
	time1 = time1 + ((double) (end1 - start1)) / CLOCKS_PER_SEC;
	time_diff=time_diff+time1;
	cout<<"After chain proces= "<<time_diff<<endl;



	// optchain=1;
	optident=1;
	optdead=0;
	optchain=0;
	// optident=0;
	//cout<<optident<<" "<<optdead<<" "<<optchain<<endl;
	if(optident==1 && optchain==0 && optdead==0)
	{
		////////////////////////////   equinodes preprocess
		int parent[n];
		vector < vector < int > > left(com);
		vector < vector < int > > alt(n);
		for(i=0;i<n;i++)
			alt[i].resize(rcgraph[i].size());
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
		start=clock();
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
		for(i=0;i<n;i++)
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
				for(j=0;j<rcgraph[i].size();j++)
				{
					alt[i][j]=outdeg[rcgraph[i][j]];
					rcgraph[i][j]=parent[rcgraph[i][j]];
				}
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		end=clock();
		cout<<"the noo is : "<<noo<<endl;
		cout<<"the perc is  : "<<((double)noo)/n<<" "<<n<<endl;
		time_diff=time_diff+((double)(end-start))/CLOCKS_PER_SEC;
		cout<<"After Equinode preprocessing= "<<time_diff<<endl;
		//////////////////////////////////

		double rank[n];
		//vector < double > rank1(n,1.0/n);
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
		int thresh=100000;
		end=clock();
		double initial[n];
		memset(initial,0,sizeof(initial));
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
// 			for(j=par[i];j<pivot;j++)
// 				long long val=computeparalleli(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
// #pragma omp parallel for private(j)
// 			for(j=pivot;j<par[i+1];j++)
// 				long long val=computeranki(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
		}
		double e=omp_get_wtime();
		cout<<e-s<<endl;
		time_diff = time_diff + ((double) (e - s)) ;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		for(i=0;i<n;i++){
			fout << initial[i] << "\n";
		}
		printf("The time taken is %lf\n",time_diff);
	}
	if(optident==1 && optchain==0 && optdead==1)	
	{
		////////////////////////////   equinodes preprocess
		int parent[n];
		vector < vector < int > > left(com);
		vector < vector < int > > alt(n);
		for(i=0;i<n;i++)
			alt[i].resize(rcgraph[i].size());
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
		start=clock();
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
		for(i=0;i<n;i++)
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
				for(j=0;j<rcgraph[i].size();j++)
				{
					alt[i][j]=outdeg[rcgraph[i][j]];
					rcgraph[i][j]=parent[rcgraph[i][j]];
				}
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		end=clock();
		cout<<"the noo is : "<<noo<<endl;
		cout<<"the perc is  : "<<((double)noo)/n<<" "<<n<<endl;
		time_diff=time_diff+((double)(end-start))/CLOCKS_PER_SEC;
		cout<<"After Equinode preprocessing= "<<time_diff<<endl;
		//////////////////////////////////

		double rank[n];
		//vector < double > rank1(n,1.0/n);
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
		int thresh=100000;
		end=clock();
		double initial[n];
		memset(initial,0,sizeof(initial));
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{

			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
				long long val=computeparallelid(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
				long long val=computerankid(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial);
		}
		double e=omp_get_wtime();
		cout<<e-s<<endl;
		time_diff = time_diff + ((double) (e - s)) ;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The time taken is %lf\n",time_diff);
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
		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<n;i++)
			members[no[i]].push_back(i);
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
				long long val=computeparallel(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerank(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
		}
		double e=omp_get_wtime();
		cout<<e-s<<endl;
		time_diff = time_diff + ((double) (e - s)) ;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The time taken is %lf\n",time_diff);
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
		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		for(i=0;i<n;i++)
			members[no[i]].push_back(i);
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
				long long val=computeparalleld(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankd(rcgraph,members[order[j]].size(),outdeg,members[order[j]],rank,initial);
			}
		}
		double e=omp_get_wtime();
		cout<<e-s<<endl;
		time_diff = time_diff + ((double) (e - s)) ;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The time taken is %lf\n",time_diff);
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
		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
		double e=omp_get_wtime();
		cout<<e-s<<endl;
		time_diff = time_diff+((double) (e - s)) ;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The total time taken is %lf\n",time_diff);

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
		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
		double e=omp_get_wtime();
		cout<<e-s<<endl;
		time_diff = time_diff+((double) (e - s)) ;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The total time taken is %lf\n",time_diff);

	}
	if(optident==1 && optchain==1 && optdead==0)
	{
		////////////////////////////   equinodes preprocess
		int parent[n];
		vector < vector < int > > left(com);
		vector < vector < int > > alt(n);
		for(i=0;i<n;i++)
			alt[i].resize(rcgraph[i].size());
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
		start=clock();
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
		for(i=0;i<n;i++)
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
				for(j=0;j<rcgraph[i].size();j++)
				{
					alt[i][j]=outdeg[rcgraph[i][j]];
					rcgraph[i][j]=parent[rcgraph[i][j]];
				}
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		end=clock();
		cout<<"the noo is : "<<noo<<endl;
		cout<<"the perc is  : "<<((double)noo)/n<<" "<<n<<endl;
		time_diff=time_diff+((double)(end-start))/CLOCKS_PER_SEC;
		cout<<"After Equinode preprocessing= "<<time_diff<<endl;
		//////////////////////////////////

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
		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
				long long val=computeparallelic(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankic(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
		}
		double e=omp_get_wtime();
		time_diff = time_diff + ((double) (e - s)) ;
		cout<<e-s<<endl;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The time taken is %lf\n",time_diff);
	}
	if(optident==1 && optchain==1 && optdead==1)
	{
		////////////////////////////   equinodes preprocess
		int parent[n];
		vector < vector < int > > left(com);
		vector < vector < int > > alt(n);
		for(i=0;i<n;i++)
			alt[i].resize(rcgraph[i].size());
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
		start=clock();
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
		for(i=0;i<n;i++)
			if(parent[i]==i) 
			{
				members[no[i]].push_back(i);
				for(j=0;j<rcgraph[i].size();j++)
				{
					alt[i][j]=outdeg[rcgraph[i][j]];
					rcgraph[i][j]=parent[rcgraph[i][j]];
				}
			}
			else
			{
				left[no[i]].push_back(i);
				noo++;
			}
		end=clock();
		cout<<"the noo is : "<<noo<<endl;
		cout<<"the perc is  : "<<((double)noo)/n<<" "<<n<<endl;
		time_diff=time_diff+((double)(end-start))/CLOCKS_PER_SEC;
		cout<<"After Equinode preprocessing= "<<time_diff<<endl;
		//////////////////////////////////

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
		int thresh=100000;
		double initial[n];
		memset(initial,0,sizeof(initial));
		double s=omp_get_wtime();
		long long no_iterations=0;
		int w;
		int gsumg=0;
		int prmax=-1;
		int compno=-1;
		for(i=0;i<par.size()-1;i++)
		{
			int pivot=par[i];
			for(int w=par[i];w<par[i+1];w++)
			{
				int sum=0;
				for(j=0;j<members[order[w]].size();j++)
					sum=sum+rgraph[members[order[w]][j]].size();
				if(prmax<sum)
				{
					compno=order[w];
					prmax=sum;
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
				long long val=computeparallelic(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
#pragma omp parallel for private(j)
			for(j=pivot;j<par[i+1];j++)
			{
				long long val=computerankic(rcgraph,alt,parent,left[order[j]],members[order[j]].size(),outdeg,members[order[j]],rank,initial,levelz,redir,powers);
			}
		}
		double e=omp_get_wtime();
		time_diff = time_diff + ((double) (e - s)) ;
		cout<<e-s<<endl;
		double sum=0;
		for(i=0;i<n;i++)
			sum=sum+rank[i];
		for(i=0;i<n;i++) rank[i]=rank[i]/sum;
		printf("The time taken is %lf\n",time_diff);
	}






	//for(i=0;i<12;i++)
	//	operation+=operations[i];
	//cout<<"the total operations are " <<operation<<endl;
	//cout<<"the total iterations are " <<no_iterations<<endl;
	/*
	   int cross[com];
	   memset(cross,0,sizeof(cross));
	   for(i=0;i<n;i++)
	   for(j=0;j<graph[i].size();j++)
	   {
	   if(nvisit[no[graph[i][j]]]==nvisit[no[i]]) continue; 
	   cross[no[graph[i][j]]]++;
	   }
	   int max_crossedges=0;
	   int avg_crossedges=0;
	   int sum_crossedges=0;
	   int no_levels=par.size()-1;
	   int no_components=com;
	   int maxno_crossedges=0;

	   for(i=0;i<com;i++){
	   max_crossedges=max(max_crossedges,cross[i]);
	   if(max_crossedges==cross[i]) maxno_crossedges=i;
	   sum_crossedges=sum_crossedges+cross[i];
	   }
	   avg_crossedges=sum_crossedges/com;
	   double sum=0;
	   for(i=0;i<n;i++)
	   sum=sum+rank[i];
	   for(i=0;i<n;i++) rank[i]=rank[i]/sum;
	   printf("The time taken is %lf\n",time_diff);
	   cout<<"no of vertices"<<n<<endl;
	   cout<<"no of edges"<<m<<endl;
	   cout<<"comp no max cross"<<maxno_crossedges<<endl;
	   cout<<"max cross="<<max_crossedges<<endl;
	   cout<<"avg cross="<<avg_crossedges<<endl;
	   cout<<"no com="<<no_components<<endl;
	   cout<<"no levels="<<no_levels<<endl;*/
	/////
	/*ofstream outmiddle;
	  outmiddle.open ("outmiddle0.txt");
	  for(i=0;i<n;i++)
	  outmiddle<<rank[i]<<endl;
	  outmiddle.close();*/
	/////

	//	for(i=0;i<n;i++) rank1[i]=rank[i]*1e8;
	//	sort(rank1.begin(),rank1.end());
	int count=0;
	/*for(i=0;i<n;i++)
	  {
	  int j=i;
	  while(fabs(rank1[j]-rank1[i])==0  && j<n)
	  {
	  j++;
	  count++;
	  }
	  i=j-1;
	  }*/
	//cout<<n<<" "<<count<<" "<<n/count<<endl;
	//	end=clock();
	//	time_diff = ((double) (end - start)) / CLOCKS_PER_SEC;
	//	printf("The time taken is %lf\n",time_diff);
	//for(i=0;i<n;i++) cout<<rank[i]<<endl;
	//	cout<<" total components="<<com<<endl;
	return 0;
}

