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

/*初始化剩余图(最初剩余图gf为原图g)*/
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

/*通过上一次剩余图和当前找到的阻塞流计算当前状态剩余图*/
void Compute_Gf(GF& gf,vector<BFLOW>& bf)
{
	vector<BFLOW>::iterator ite;
	for(ite = bf.begin(); ite != bf.end(); ite++)
	{
		gf[ite->u][ite->v] -= ite->flow;  
		gf[ite->v][ite->u] += ite->flow;
	}
} 

/*直接对剩余图Gf宽度优先遍历得到分层图msn*/
void BFS2MSN(GF& gf,int source,MSN& msn)
{
	bool isVisited[MAX];                                            /*访问标志*/
	for(int i = 0; i < MAX; i++)
	{
		isVisited[i] = false;  
	} 

	msn.nVStage[source] = 0;
	isVisited[source] = true;

	queue<int> iQueue;
	int u;                                                          /*当前处理节点*/
	iQueue.push(source);
	while(!iQueue.empty())
	{
		u = iQueue.front();
		iQueue.pop();
		for(int i = 1; i <= msn.nV; i++)
		{
			if(gf[u][i] == 0)                                       /*剩余图中边容量为0*/
			{
				msn.nMSN[u][i] = 0;
			}
			else
			{
				if(!isVisited[i])                                   /*节点没有访问*/
				{
					isVisited[i] = true;
					msn.nVStage[i] = msn.nVStage[u]+1;
					msn.nMSN[u][i] = gf[u][i];
					iQueue.push(i);
				}
				else if(msn.nVStage[u]+1 == msn.nVStage[i])
				{   /*节点i已访问,而当前节点u又恰好是i上一层*/
					msn.nMSN[u][i] = gf[u][i]; 
				}
				else
				{
					msn.nMSN[u][i] = 0;                             /*同一层间不设边*/
				}
			}
		}
	}
}

/*DFS方法求解分成网络msn中阻塞流bf(使用stack来解决)*/
void DFS2BF(MSN& msn,int n,int source,int sink,vector<BFLOW>& bf)
{
	bool isDeleted[MAX];                                           /*合理情况应该在堆空间中申请这个空间*/
	for(int i = 0; i < MAX; i++)
	{
		isDeleted[i] = false;
	}

	int vStack[MAX];                                                /*规定第0个空间不用*/
	int top = 1;
	vStack[top] = source;
	while(top >= 1)
	{
		if(vStack[top] == sink)                                     /*找到一条从源点到汇点路径*/
		{
			assert(top >= 2);
			int minC = msn.nMSN[vStack[1]][vStack[2]];
			int i;
			for(i = 2; i < top; i++)                                /*求解这条路径上的最小边*/
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
						ite->flow += minC;                           /*这条边已经存在*/
						break;
					}
				}
				if(ite == bf.end())                                  /*还没有出现过*/
				{
					BFLOW tempBF;
					tempBF.u = vStack[i-1];
					tempBF.v = vStack[i];
					tempBF.flow = minC;
					bf.push_back(tempBF);
				}
				
				msn.nMSN[vStack[i-1]][vStack[i]] -= minC;             /*沿着路径减少边上容量*/
				if(msn.nMSN[vStack[i-1]][vStack[i]] == 0)             /*满容量边*/
				{
					top = i-1;
				}
			}
		}
		else                                                          /*不能继续深度下去，标记并回溯*/
		{
			int i;
			for(i = 1; i <= n; i++)
			{
				if(!isDeleted[i] && msn.nMSN[vStack[top]][i] > 0)     /*有指向下一层的边*/
				{
					vStack[++top] = i;                                /*节点入栈*/
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

/*通过源点流出的流得到最大流值*/
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

/*Dinic 算法O(nmn)求解最大流f和取得最大流时的剩余图gf*/
int Dinic(Graph& g,int source,int sink,Flow& f,GF& gf)
{
	Init_F(f);                                                     /*初始化f[u,v]和f[v,u];*/
	Init_Gf(g,gf);                                                 /*直接利用原图g初始化剩余图gf*/

	MSN msn;
	msn.Init(g.nV);
	BFS2MSN(gf,source,msn);                                       /*得到分层网络msn*/

	vector<BFLOW> bf;
	/*通过sink层次来判断是否还在分层网络中(-1)不在*/
	while(msn.nVStage[sink] != -1) 
	{
		/*采用深度优先寻找阻塞流bf(采用向量容器保存）*/
		DFS2BF(msn,g.nV,source,sink,bf);

	    /*通过阻塞流bf增加流f*/
		vector<BFLOW>::iterator ite;
		for(ite = bf.begin(); ite != bf.end(); ite++)
		{
			f[ite->u][ite->v] += ite->flow;
			f[ite->v][ite->u] = -f[ite->u][ite->v];
		}
		
		Compute_Gf(gf,bf);                                        /*重新计算剩余图*/
		bf.clear();
		msn.Init(g.nV);
		BFS2MSN(gf,source,msn);                                   /*重新计算剩余图对应的层次图*/
	}
	return MaxFlowValue(f,source,g.nV);
}

/*通过剩余图求解当前最大流对应的最小割(层次遍历)*/
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