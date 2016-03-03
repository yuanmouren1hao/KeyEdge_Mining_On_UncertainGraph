#include "stdafx.h"

#include "graph.h"
#include <queue>
#include <vector>

using namespace std;

Graph::Graph()
:nId(0),nV(0),nE(0)
{

}

Graph::~Graph()
{
	
}

void Graph::Init()
{
	nId = 0;
	nV = 0; 
	nE = 0;
	for(int u = 1; u <  MAX; u++)
	{
		for(int v = 1; v < MAX; v++)
		{
			matrix[u][v].dP = 0;
			matrix[u][v].iC = 0;
			matrix[u][v].iLabel = 0;
		}
	}
}

void Init_F(Flow& f)
{
	for(int u = 0; u <  MAX; u++)
	{
		for(int v = 0; v < MAX; v++)
		{
			f[u][v] = 0;
		}
	}
}


void MSN_STRUCT::Init(int nVnum)
{
	nV = nVnum;
	for(int i = 0; i <= nV; i++)
	{
		nVStage[i] = -1;
	}
	for(int i = 0; i <= nV; i++)
	{
		for(int j = 0; j <= nV; j++)
		{
			nMSN[i][j] = 0;
		}
	}
}

/*��ʼ��ʣ��ͼ(���ʣ��ͼgfΪԭͼg)*/
void Init_Gf(Graph& g,GF& gf)
{
	for(int u = 1; u <= g.nV; u++)
	{
		for(int v = 1; v <= g.nV; v++)
		{
			gf[u][v] = g.matrix[u][v].iC;
		}
	}
}

/*ͨ����һ��ʣ��ͼ�͵�ǰ�ҵ������������㵱ǰ״̬ʣ��ͼ*/
void Compute_Gf(GF& gf,vector<BFLOW>& bf)
{
	vector<BFLOW>::iterator ite;
	for(ite = bf.begin(); ite != bf.end(); ite++)
	{
		gf[ite->u][ite->v] -= ite->flow;  
		gf[ite->v][ite->u] += ite->flow;
	}
} 

/*ֱ�Ӷ�ʣ��ͼGf������ȱ����õ��ֲ�ͼmsn*/
void BFS2MSN(GF& gf,int source,MSN& msn)
{
	bool isVisited[MAX];                                            /*���ʱ�־*/
	for(int i = 0; i < MAX; i++)
	{
		isVisited[i] = false;  
	} 

	msn.nVStage[source] = 0;
	isVisited[source] = true;

	queue<int> iQueue;
	int u;                                                          /*��ǰ����ڵ�*/
	iQueue.push(source);
	while(!iQueue.empty())
	{
		u = iQueue.front();
		iQueue.pop();
		for(int i = 1; i <= msn.nV; i++)
		{
			if(gf[u][i] == 0)                                       /*ʣ��ͼ�б�����Ϊ0*/
			{
				msn.nMSN[u][i] = 0;
			}
			else
			{
				if(!isVisited[i])                                   /*�ڵ�û�з���*/
				{
					isVisited[i] = true;
					msn.nVStage[i] = msn.nVStage[u]+1;
					msn.nMSN[u][i] = gf[u][i];
					iQueue.push(i);
				}
				else if(msn.nVStage[u]+1 == msn.nVStage[i])
				{   /*�ڵ�i�ѷ���,����ǰ�ڵ�u��ǡ����i��һ��*/
					msn.nMSN[u][i] = gf[u][i]; 
				}
				else
				{
					msn.nMSN[u][i] = 0;                             /*ͬһ��䲻���*/
				}
			}
		}
	}
}

/*DFS�������ֳ�����msn��������bf(ʹ��stack�����)*/
void DFS2BF(MSN& msn,int n,int source,int sink,vector<BFLOW>& bf)
{
	bool isDeleted[MAX];                                           /*�������Ӧ���ڶѿռ�����������ռ�*/
	for(int i = 0; i < MAX; i++)
	{
		isDeleted[i] = false;
	}

	int vStack[MAX];                                                /*�涨��0���ռ䲻��*/
	int top = 1;
	vStack[top] = source;
	while(top >= 1)
	{
		if(vStack[top] == sink)                                     /*�ҵ�һ����Դ�㵽���·��*/
		{
			assert(top >= 2);
			int minC = msn.nMSN[vStack[1]][vStack[2]];
			int i;
			for(i = 2; i < top; i++)                                /*�������·���ϵ���С��*/
			{
				if(msn.nMSN[vStack[i]][vStack[i+1]] < minC)
				{
					minC = msn.nMSN[vStack[i]][vStack[i+1]];
				}
			}
			for(i = top; i > 1; i--) 
			{
				vector<BFLOW>::iterator ite;
				for(ite = bf.begin(); ite != bf.end(); ite++)
				{
					if(ite->u == vStack[i-1] && ite->v == vStack[i])
					{
						ite->flow += minC;                           /*�������Ѿ�����*/
						break;
					}
				}
				if(ite == bf.end())                                  /*��û�г��ֹ�*/
				{
					BFLOW tempBF;
					tempBF.u = vStack[i-1];
					tempBF.v = vStack[i];
					tempBF.flow = minC;
					bf.push_back(tempBF);
				}
				
				msn.nMSN[vStack[i-1]][vStack[i]] -= minC;             /*����·�����ٱ�������*/
				if(msn.nMSN[vStack[i-1]][vStack[i]] == 0)             /*��������*/
				{
					top = i-1;
				}
			}
		}
		else                                                          /*���ܼ��������ȥ����ǲ�����*/
		{
			int i;
			for(i = 1; i <= n; i++)
			{
				if(!isDeleted[i] && msn.nMSN[vStack[top]][i] > 0)     /*��ָ����һ��ı�*/
				{
					vStack[++top] = i;                                /*�ڵ���ջ*/
					break;
				}
			}
			if( n+1 == i)
			{
				isDeleted[vStack[top]] = true;
				top--;
			}
		}
	}
}

/*ͨ��Դ�����������õ������ֵ*/
int MaxFlowValue(Flow& f,int source,int n)
{
	int maxFV = 0;
	for(int i = 1; i <= n; i++)
	{
		if(source != i && f[source][i] != 0)
		{
			maxFV += f[source][i];
		}
	}
	return maxFV;
}

/*Dinic �㷨O(nmn)��������f��ȡ�������ʱ��ʣ��ͼgf*/
int Dinic(Graph& g,int source,int sink,Flow& f,GF& gf)
{
	Init_F(f);                                                     /*��ʼ��f[u,v]��f[v,u];*/
	Init_Gf(g,gf);                                                 /*ֱ������ԭͼg��ʼ��ʣ��ͼgf*/

	MSN msn;
	msn.Init(g.nV);
	BFS2MSN(gf,source,msn);                                       /*�õ��ֲ�����msn*/

	vector<BFLOW> bf;
	/*ͨ��sink������ж��Ƿ��ڷֲ�������(-1)����*/
	while(msn.nVStage[sink] != -1) 
	{
		/*�����������Ѱ��������bf(���������������棩*/
		DFS2BF(msn,g.nV,source,sink,bf);

	    /*ͨ��������bf������f*/
		vector<BFLOW>::iterator ite;
		for(ite = bf.begin(); ite != bf.end(); ite++)
		{
			f[ite->u][ite->v] += ite->flow;
			f[ite->v][ite->u] = -f[ite->u][ite->v];
		}
		
		Compute_Gf(gf,bf);                                        /*���¼���ʣ��ͼ*/
		bf.clear();
		msn.Init(g.nV);
		BFS2MSN(gf,source,msn);                                   /*���¼���ʣ��ͼ��Ӧ�Ĳ��ͼ*/
	}
	return MaxFlowValue(f,source,g.nV);
}

/*ͨ��ʣ��ͼ��⵱ǰ�������Ӧ����С��(��α���)*/
void MinCut(GF& gf,int n,int source,vector<int>& S)  
{
	bool isVisited[MAX];
	for(int i = 0; i < MAX; i++)
	{
		isVisited[i] = false;
	}
	
	queue<int> iQueue;

	isVisited[source] = true;
	iQueue.push(source);
	int u = source;
	while(!iQueue.empty())
	{
		u = iQueue.front();
		S.push_back(u);
		iQueue.pop();
		for(int v = 1; v <= n; v++)
		{
			if(!isVisited[v] && gf[u][v] > 0)
			{
				isVisited[v] = true;
				iQueue.push(v);
			}
		}
	}
}