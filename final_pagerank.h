#include<cstdio>
#include<iostream>
#include<cmath>
#include<vector>
#include<cstring>
#include<omp.h>

using namespace std;


double damp=0.85;

//int countit=20;
//double value=1e-8;

//  toposeqmodif_strong.cpp  //
int  computeparallel(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double damp=0.85;
	double thres=1e-10;
	double value=1e-12/n;
	int i,j;
	vector < double > curr(n);
	//bool marked[n];
	//memset(marked,0,sizeof(marked));
	double error=0;
	time_t start,end;
	double time_diff;
	//start=clock();
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
	//omp_set_dynamic(0);
	//omp_set_num_threads(12);
	for(i=0;i<n;i++)
	{
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}
	}
	do
	{
		error=0;
#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{
			//if(marked[i]<countit)
			{
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				//operations[omp_get_thread_num()]+=graph[node].size();
				//error=error+fabs(curr[i]-rank[node]);
			}
		}
		for(i=pivot;i<n;i++)
		{
			//if(marked[i]<countit)
			{
				int node=mapit[i];
				double anse[12];
				memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
				for(j=0;j<graph[node].size();j++)
					anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/outdeg[graph[node][j]];
				double ans=0;
				//operations[0]+=graph[node].size();
				for(j=0;j<12;j++)
					ans=ans+anse[j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
				//error=error+fabs(curr[i]-rank[node]);
			}
		}
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<n;i++)
			//if(marked[i]<countit)
			anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
#pragma omp parallel for private(i) 
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
	//end=clock();
	time_diff = ((double ) (end - start)) / CLOCKS_PER_SEC;
	curr.clear();
	return iterations;
	//cout<<iterations<<endl;
	//	printf("The time taken is %0.9lf\n",time_diff);
	//for(i=0;i<n;i++)
	//	cout<<curr[i]<<endl;
}

long long int computerank(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double value=1e-12/n;
	//bool marked[n];
	//memset(marked,0,sizeof(marked));
	double error=0;
	time_t start,end;
	double time_diff;
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
	//end=clock();
	time_diff = ((double ) (end - start)) / CLOCKS_PER_SEC;
	curr.clear();
	//cout<<iterations<<endl;
	//	printf("The time taken is %0.9lf\n",time_diff);
	//for(i=0;i<n;i++)
	//	cout<<curr[i]<<endl;
	return sumiterations;
}






// mainalgo.cpp //
int  computeparalleld(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double thres=1e-10;
	double dis=1e-12;
	double value=((dis)*10.0)/n;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n);
#pragma omp parallel for private(i) 
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	bool marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
	for(i=0;i<n;i++)
	{
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}
	}
	do
	{
		error=0;
#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{
			if(!marked[i])
			{
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}
		}

		for(i=pivot;i<n;i++)
		{
			if(!marked[i])
			{
				int node=mapit[i];
				double anse[12];
				memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
				for(j=0;j<graph[node].size();j++)
					anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/outdeg[graph[node][j]];
				double ans=0;
				for(j=0;j<12;j++)
					ans=ans+anse[j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}
		}
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<n;i++)
			if(!marked[i])
				anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
		iterations++;
#pragma omp parallel for private(i) 
		for(i=0;i<n;i++)
			if(!marked[i]) rank[mapit[i]]=curr[i];
		if(iterations%20==0)
		{
#pragma omp parallel for private(i) 
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
	}while(error > thres );
	return iterations;
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
	double bound=1e-5;
	bool  marked[n];
	memset(marked,0,sizeof(marked));
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





// oritoposeqmodif_strong.cpp//
int computeparallelc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	int i,j;
	vector < double > curr(n);
	double error=0;
	time_t start,end;
	double time_diff;
	start=clock();
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
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
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   
	vector < vector < int > > spare(12);
#pragma omp parallel for private(i,j) 
	for(i=0;i<pivot;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	for(i=pivot;i<limit;i++)
	{   
		int node=mapit[i];
#pragma omp parallel for private(j) 
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	do  
	{   
		error=0;

#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{   
			int node=mapit[i];
			double ans=0;
			for(j=0;j<graph[node].size();j++)
			{   
				ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
			}   
			curr[i]=randomp+damp*ans+initial[mapit[i]];
		}   
		for(i=pivot;i<limit;i++)
		{   
			int node=mapit[i];
			double anse[12];
			memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
			for(j=0;j<graph[node].size();j++)
			{   
				anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/outdeg[graph[node][j]];
			}   
			double ans=0;
			for(j=0;j<12;j++)
				ans=ans+anse[j];
			curr[i]=randomp+damp*ans+initial[mapit[i]];
		}   
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
		{   
			anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		}   
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
			rank[mapit[i]]=curr[i];
		iterations++;
#pragma omp parallel for private(i,j) 
		for(i=0;i<12;i++)
		{   
			for(j=0;j<spare[i].size();j++)
			{   
				double val=powers[level[spare[i][j]]];
				rank[spare[i][j]]=rank[redir[spare[i][j]]]*val+(1.0-val)/graph.size();
			}   
		}   
	}while(error > thres );
#pragma omp parallel for private(i)
	for(i=limit;i<n;i++)
	{   
		int node=mapit[i];
		double val=powers[level[node]];
		rank[node]=rank[redir[node]]*val+(1.0-val)/graph.size();
	}   
	return iterations;
}



long long int computerankc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double value=((1e-12)*10.0)/double ( n );
	double error=0;
	time_t start,end;
	double time_diff;
	start=clock();
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	double s,e;
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
	int sumiterations=0;
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



// deadoritoposeqmodif_strong.cpp //


int computeparalleldc(vector < vector < int > > & graph,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n,1.0/n);
	bool marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	time_t start,end;
	double time_diff;
	start=clock();
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
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
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   
	vector < vector < int > > spare(12);
#pragma omp parallel for private(i,j) 
	for(i=0;i<pivot;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	for(i=pivot;i<limit;i++)
	{   
		int node=mapit[i];
#pragma omp parallel for private(j) 
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	do  
	{   
		error=0;

#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{   
			if(!marked[i])
			{   
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{   
					ans=ans+rank[graph[node][j]]/outdeg[graph[node][j]];
				}   
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		for(i=pivot;i<limit;i++)
		{   
			if(!marked[i])
			{   
				int node=mapit[i];
				double anse[12];
				memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
				for(j=0;j<graph[node].size();j++)
				{   
					anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/outdeg[graph[node][j]];
				}   
				double ans=0;
				for(j=0;j<12;j++)
					ans=ans+anse[j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
		{   
			if(!marked[i])
				anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		}   
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
			if(!marked[i])
				rank[mapit[i]]=curr[i];
		iterations++;
		if(iterations%20==0)
#pragma omp parallel for private(i) 
			for(i=0;i<limit;i++)
			{   
				if(!marked[i])
				{   
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}   
			}   
#pragma omp parallel for private(i,j) 
		for(i=0;i<12;i++)
		{   
			for(j=0;j<spare[i].size();j++)
			{   
				double val=powers[level[spare[i][j]]];
				rank[spare[i][j]]=rank[redir[spare[i][j]]]*val+(1.0-val)/graph.size();
			}   
		}   
	}while(error > thres );
#pragma omp parallel for private(i)
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
	time_t start,end;
	double time_diff;
	start=clock();
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	double s,e;
	s=omp_get_wtime();
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
	e=omp_get_wtime();
	int sumiterations=0;
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





// toposeqmodifident_strong.cpp //

int  computeparalleli(vector < vector < int > > & graph , vector < vector < int > > & alt,int parent[],vector <int > left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double damp=0.85;
	double thres=1e-10;
	double value=1e-12/n;
	int i,j;
	vector < double > curr(n);
	double error=0;
	time_t start,end;
	double time_diff;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
	for(i=0;i<n;i++)
	{   
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   
	do  
	{   
		error=0;
#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{   
			{   
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/alt[node][j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		for(i=pivot;i<n;i++)
		{   
			{   
				int node=mapit[i];
				double anse[12];
				memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
				for(j=0;j<graph[node].size();j++)
					anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/alt[node][j];
				double ans=0;
				for(j=0;j<12;j++)
					ans=ans+anse[j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<n;i++)
			anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
#pragma omp parallel for private(i) 
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
	time_diff = ((double ) (end - start)) / CLOCKS_PER_SEC;
	curr.clear();
	return iterations;
}

long long int computeranki(vector < vector < int > > & graph, vector < vector < int > > & alt, int parent[],vector < int > left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double value=1e-12/n;
	double error=0;
	time_t start,end;
	double time_diff;
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	int sumiterations=0;
	do
	{
		error=0;
		for(i=0;i<n;i++)
		{
			{
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/alt[node][j];
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
	sumiterations=iterations*n;
	return iterations;
	time_diff = ((double ) (end - start)) / CLOCKS_PER_SEC;
	curr.clear();
	return sumiterations;
}


// mainalgoident.cpp //


int  computeparallelid(vector < vector < int > > & graph,vector < vector < int > > & alt,int parent[],vector < int > & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double thres=1e-10;
	double dis=1e-12;
	double value=((dis)*10.0)/n;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n);
#pragma omp parallel for private(i) 
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	bool marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
	for(i=0;i<n;i++)
	{   
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   
	do  
	{   
		error=0;
#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{   
			if(!marked[i])
			{   
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
					ans=ans+rank[graph[node][j]]/alt[node][j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   

		for(i=pivot;i<n;i++)
		{   
			if(!marked[i])
			{   
				int node=mapit[i];
				double anse[12];
				memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
				for(j=0;j<graph[node].size();j++)
					anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/alt[node][j];
				double ans=0;
				for(j=0;j<12;j++)
					ans=ans+anse[j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<n;i++)
			if(!marked[i])
				anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
		iterations++;
#pragma omp parallel for private(i) 
		for(i=0;i<n;i++)
			if(!marked[i]) rank[mapit[i]]=curr[i];
		if(iterations%20==0)
		{   
#pragma omp parallel for private(i) 
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
	}while(error > thres );
	for(i=0;i<left.size();i++)
		rank[left[i]]=rank[parent[left[i]]];
	return iterations;
}

long long int computerankid(vector < vector < int > > & graph,vector < vector < int > >  & alt,int parent[],vector < int > & left, int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[])
{
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n);
	for(i=0;i<n;i++)
		prev[i]=rank[mapit[i]];
	double dis=1e-12;
	double value=(dis*10.0)/n;
	double bound=1e-5;
	bool  marked[n];
	memset(marked,0,sizeof(marked));
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
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{
					ans=ans+rank[graph[node][j]]/alt[node][j];
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
	return iterations;
}

// oritoposeqmodifident_strong.cpp //


int computeparallelic(vector < vector < int > > & graph,vector < vector < int >  > & alt,int parent[],vector <int > & left, int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	int i,j;
	vector < double > curr(n);
	double error=0;
	time_t start,end;
	double time_diff;
	start=clock();
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
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
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   
	vector < vector < int > > spare(12);
#pragma omp parallel for private(i,j) 
	for(i=0;i<pivot;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	for(i=pivot;i<limit;i++)
	{   
		int node=mapit[i];
#pragma omp parallel for private(j) 
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	do  
	{   
		error=0;

#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{   
			int node=mapit[i];
			double ans=0;
			for(j=0;j<graph[node].size();j++)
			{   
				ans=ans+rank[graph[node][j]]/alt[node][j];
			}   
			curr[i]=randomp+damp*ans+initial[mapit[i]];
		}   
		for(i=pivot;i<limit;i++)
		{   
			int node=mapit[i];
			double anse[12];
			memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
			for(j=0;j<graph[node].size();j++)
			{   
				anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/alt[node][j];
			}   
			double ans=0;
			for(j=0;j<12;j++)
				ans=ans+anse[j];
			curr[i]=randomp+damp*ans+initial[mapit[i]];
		}   
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
		{   
			anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		}   
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
			rank[mapit[i]]=curr[i];
		iterations++;
#pragma omp parallel for private(i,j) 
		for(i=0;i<12;i++)
		{   
			for(j=0;j<spare[i].size();j++)
			{   
				double val=powers[level[spare[i][j]]];
				rank[spare[i][j]]=rank[redir[spare[i][j]]]*val+(1.0-val)/graph.size();
			}   
		}   
	}while(error > thres );
#pragma omp parallel for private(i)
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

long long int computerankic(vector < vector < int > > & graph,vector < vector < int > > &  alt,int parent[],vector < int > & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	int i,j;
	vector < double > curr(n);
	double value=((1e-12)*10.0)/double ( n );
	double error=0;
	time_t start,end;
	double time_diff;
	start=clock();
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	double s,e;
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
	int sumiterations=0;
	do
	{
		error=0;
		for(i=0;i<limit;i++)
		{
			int node=mapit[i];
			double ans=0;
			for(j=0;j<graph[node].size();j++)
			{
				ans=ans+rank[graph[node][j]]/alt[node][j];
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


// deadoritoposeqmodifident_strong.cpp //


int computeparallelidc(vector < vector < int > > & graph,vector < vector < int > > & alt,int parent[],vector <int> & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
{
	double damp=0.85;
	double thres=1e-10;
	double value=((1e-12)*10.0)/double ( n );
	int i,j;
	vector < double > curr(n);
	vector < double > prev(n,1.0/n);
	bool marked[n];
	memset(marked,0,sizeof(marked));
	double error=0;
	time_t start,end;
	double time_diff;
	start=clock();
	int iterations=0;
	double randomp=(1-damp)/graph.size();
	int pivot=0;
	int thresh=10000;
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
	for(i=0;i<limit;i++)
	{   
		int node=mapit[i];
		if(graph[node].size()<thresh)
		{   
			int temp=mapit[pivot];
			mapit[pivot]=node;
			mapit[i]=temp;
			pivot++;
		}   
	}   
	vector < vector < int > > spare(12);
#pragma omp parallel for private(i,j) 
	for(i=0;i<pivot;i++)
	{   
		int node=mapit[i];
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	for(i=pivot;i<limit;i++)
	{   
		int node=mapit[i];
#pragma omp parallel for private(j) 
		for(j=0;j<graph[node].size();j++)
			if(redir[graph[node][j]]!=graph[node][j]) spare[omp_get_thread_num()].push_back(graph[node][j]);
	}   
	do  
	{   
		error=0;

#pragma omp parallel for private(i,j) 
		for(i=0;i<pivot;i++)
		{   
			if(!marked[i])
			{   
				int node=mapit[i];
				double ans=0;
				for(j=0;j<graph[node].size();j++)
				{   
					ans=ans+rank[graph[node][j]]/alt[node][j];
				}   
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		for(i=pivot;i<limit;i++)
		{   
			if(!marked[i])
			{   
				int node=mapit[i];
				double anse[12];
				memset(anse,0,sizeof(anse));
#pragma omp parallel for private(j) 
				for(j=0;j<graph[node].size();j++)
				{   
					anse[omp_get_thread_num()]=anse[omp_get_thread_num()]+rank[graph[node][j]]/alt[node][j];
				}   
				double ans=0;
				for(j=0;j<12;j++)
					ans=ans+anse[j];
				curr[i]=randomp+damp*ans+initial[mapit[i]];
			}   
		}   
		double anse[12];
		memset(anse,0,sizeof(anse));
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
		{   
			if(!marked[i])
				anse[omp_get_thread_num()]=max(anse[omp_get_thread_num()],fabs(curr[i]-rank[mapit[i]]));
		}   
		for(i=0;i<12;i++)
			error=max(error,anse[i]);
#pragma omp parallel for private(i) 
		for(i=0;i<limit;i++)
			if(!marked[i])
				rank[mapit[i]]=curr[i];
		iterations++;
		if(iterations%20==0)
#pragma omp parallel for private(i) 
			for(i=0;i<limit;i++)
			{   
				if(!marked[i])
				{   
					if(fabs(prev[i]-curr[i]) < value )marked[i]=1;
					else prev[i]=curr[i];
				}   
			}   
#pragma omp parallel for private(i,j) 
		for(i=0;i<12;i++)
		{   
			for(j=0;j<spare[i].size();j++)
			{   
				double val=powers[level[spare[i][j]]];
				rank[spare[i][j]]=rank[redir[spare[i][j]]]*val+(1.0-val)/graph.size();
			}   
		}   
	}while(error > thres );
#pragma omp parallel for private(i)
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


long long int computerankidc(vector < vector < int > > & graph,vector < vector <int > > & alt,int parent[],vector < int > & left,int n,int outdeg[],vector < int > &  mapit,double rank[],double initial[],int level[],int redir[],double powers[])
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
	time_t start,end;
	double time_diff;
	start=clock();
	long long  iterations=0;
	double randomp=(1-damp)/graph.size();
	int limit=0;
	double s,e;
	s=omp_get_wtime();
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
	e=omp_get_wtime();
	int sumiterations=0;
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
					ans=ans+rank[graph[node][j]]/alt[node][j];
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
