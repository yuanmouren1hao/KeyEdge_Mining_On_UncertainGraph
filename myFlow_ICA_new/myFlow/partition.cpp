#include "stdafx.h"

#include "partition.h"
#include "graph.h"
#include "resultprint.h"
#include <memory>
#include <stack>                                                                 /*使用栈来优化空间*/
#include <iostream>

/*引入状态的头文件*/
#include "state.h"

using namespace std;

StateSet::StateSet(int num):numE(num)
{	
}

StateSet::~StateSet()
{
}

/*用于改善空间的优化结构*/
/*构造函数初始化*/
Collection::Collection(int num):numE(num),I(0),j(0)
{
	for(int i = 0; i <= numE; i++)
	{
		this->lower[i] = 0;  
		this->upper[i] = 1;
		/*初始情况下Ad_C集合为空(由I保证)*/
	}
}

/*构造函数初始化，有上界和下界的构造函数*/
Collection::Collection(vector<int> down_, vector<int>top_, int num)
{
	this->numE = num;
	this->I = 0;
	this->j = 0;
	for(int i = 0; i <= num; i++)
	{
		this->lower[i] = down_[i];  
		this->upper[i] = top_[i];
	}
}

Collection::~Collection()
{
}

/*将状态向量转换成对应图g_target*/
void FromG2G(Graph& g_original,Graph& g_target,int* lu)
{
	g_target.Init();
	g_target.nV = g_original.nV;  
	g_target.nId = g_original.nId;                                                /*图的编号可以不用*/
	for(int u = 1; u <= g_original.nV; u++)
	{
		for(int v = 1; v <= g_original.nV; v++)
		{
			if((g_original.matrix[u][v].iC > 0) //边
				&& (lu[g_original.matrix[u][v].iLabel] == 1))
			{
				g_target.matrix[u][v].dP = g_original.matrix[u][v].dP;
				g_target.matrix[u][v].iC = g_original.matrix[u][v].iC;
				g_target.matrix[u][v].iLabel = g_original.matrix[u][v].iLabel;
				g_target.nE++;  
			}
		}
	}
}

/*通过网络和相应的最大流分布求解划分点x0*/
void ComputeX0(Graph& g,Flow Fd,int* x0)
{
	for(int u = 1; u <= g.nV; u++)
	{
		for(int v = 1; v <= g.nV; v++)
		{
			if(Fd[u][v] > 0) /**********/
			{
				x0[g.matrix[u][v].iLabel] = 1;
			}
		}
	}
}

/*子图l被子图空间中子图包含的概率*/
double CalculateP(Graph& g,int* l)
{
	double dp = 1.0;			
	for(int u = 1; u <= g.nV; u++)
	{
		for(int v = 1; v <= g.nV; v++)
		{
			if(g.matrix[u][v].iC > 0) 
			{
				if(l[g.matrix[u][v].iLabel] == 1)                                 /*只对l中为1的边记概率*/
				{
					dp = dp*(g.matrix[u][v].dP);
				}
			}
		}
	}
	return dp;
}

/*采用Chin-Chia Jane and Yih-Wenn Laih的方法得到C(j)(1-q)*/
void GetCurrentC_j(Collection &c,StateSet &Cj)
{
	/*分别初始化Cj的上下界*/
	memcpy(Cj.lower,c.lower,(Cj.numE+1)*sizeof(int)); 
	memcpy(Cj.upper,c.upper,(Cj.numE+1)*sizeof(int));

	/*初始(第一次)Collection={|E|,(00,...0),(11,...,1),(NULL),0,0}
	单独处理一下*/
	if(c.I == 0){ return ; }

	/*保证Ad_c[1],...,保证Ad_c[j-1]这些边都要大于划分线
	在二元情况下，相当于要大于下界(本质就是取上界)*/
	for(int i = 1; i < c.j; i++)
	{
		Cj.lower[c.Ad_C[i]] = Cj.upper[c.Ad_C[i]];
	}
	/*对第j个Ad_[j]对应的pivot要小于划分线
	在二元情况下，相当于是要小于上界(本质就是取下界)*/
	Cj.upper[c.Ad_C[c.j]] = Cj.lower[c.Ad_C[c.j]];
}

/*采用Chin-Chia Jane and Yih-Wenn Laih方法获得C0*/
void GetCurrentC0(Collection &c,StateSet &C0)
{
	memcpy(C0.lower,c.lower,(C0.numE+1)*sizeof(int));
	memcpy(C0.upper,c.upper,(C0.numE+1)*sizeof(int));
	for(int i = 1; i <= c.I; i++)
	{
		C0.lower[c.Ad_C[i]] = C0.upper[c.Ad_C[i]];
	}
}

/*通过闭合区间与划分线得到pivot集合ad(c)*/
/*函数返回ad(c)集合的大小和Ad_c*/
int GetAd_C(StateSet &Cj,int *x0,int *Ad_c) 
{
	int iSize = 0; 
	for(int i = 1; i <= Cj.numE; i++)
	{
		if(x0[i] > Cj.lower[i]) //f[i] > l[i]
		{
			iSize++;
			Ad_c[iSize] = i; 
		}
	}
	return iSize;
}

/*下面利用划分线进行划分,得到能够导出
C(1),...,C(I)的collection*/
void GetNextCollection(StateSet &Cj,int *x0,Collection &nextCollection)
{
	nextCollection.j = 1;                                                         /*设定需要处理要求的Cj*/
	nextCollection.I = GetAd_C(Cj,x0,nextCollection.Ad_C);                        /*设定I和Ad_c[]*/
	memcpy(nextCollection.lower,Cj.lower,(nextCollection.numE+1)*sizeof(int));    /*设定lower[]*/
	memcpy(nextCollection.upper,Cj.upper,(nextCollection.numE+1)*sizeof(int));    /*设定upper[]*/
}

/*保存当前能够取到最大流值的闭合区间中下界子图概率大者*/
#define SAVEM(C) \
	tmpP = CalculateP(g,(C).lower); \
	if(tmpP > dR)  \
	{  dR = tmpP; memcpy(lofMaxP,(C).lower,(nE+1)*sizeof(int)); } 



/*保存所有的状态*/
void saveAllState(StateSet c,Lower_subGraph * StateMtrix)
{
	StateMtrix->State_Num++;
	StateMtrix->Edge_Num=c.numE;
	// 将下界子图保存在状态矩阵中
	for (int ii=1;ii<=c.numE;ii++)
	{
		if (c.lower[ii] == 0 && c.upper[ii] == 0)
		{
			StateMtrix->State[StateMtrix->State_Num][ii] = 0;
		}else if (c.lower[ii] == 1 && c.upper[ii] == 1)
		{
			StateMtrix->State[StateMtrix->State_Num][ii] = 1;
		}else if (c.lower[ii] == 0 && c.upper[ii] == 1)
		{
			StateMtrix->State[StateMtrix->State_Num][ii] = 2;
		}
	}
	return;
}

	  
/*通过对可能事件模型进行划分得到最可靠最大流分布*/
double  old_GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd)
{	
	assert((source > 0 && source <= g.nV) && (sink > 0 && sink <= g.nV));         /*保证输入合法*/
	
	int nE = g.nE;                                                                /*向量大小是不变的*/
	double dR = 0.0;                                                              /*最可靠最大流分布概率*/

	Flow Fd;  
	GF gf;                                                                        /*剩余图*/
	int Fmax = Dinic(g,source,sink,Fd,gf);                                        /*得到所有边都存在的最大流*/
	assert(Fmax >= 0);                                                             /*保证最大流大于0有意义*/

	/*完成一个六元组的初始化*/
	Collection c(nE);
	stack<Collection> cStack;                                                     /*存放需要进一步划分的子图空间*/
	cStack.push(c);
	
	Graph cur_g;
	StateSet Cj(nE),C0(nE);                                                       /*存放当前需要处理的闭合区间[]*/
	Collection next_c(nE);                                                        /*保存下一个需要进栈处理的collection*/

	int tmpF = 0;
	double tmpP = 0.0;
	int * lofMaxP = new int[nE+1];                                                 /*临时存储满足要求的闭合区间中概率最大的下界*/
	int * x0 = new int[nE+1];                                                      /*划分线*/
	memset(lofMaxP,0,(nE+1)*sizeof(int));
	memset(x0,0,(nE+1)*sizeof(int));
	while(!cStack.empty())
	{	
		GetCurrentC_j(cStack.top(),Cj);                                            /*取得当前需要处理的闭合区间[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*变换成下一个待处理的状态*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G(g,cur_g,Cj.lower);                                                 /*下界向量对应的子图*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf); 
		if(tmpF >= Fmax)                                                           /*完备区间下界对应子图能够满足最大流*/
		{
			SAVEM(Cj);
		}  
		else if((FromG2G(g,cur_g,Cj.upper),                                        /*对上界对应的子图求最大流*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*同时满足F(lc) < Fmax <= F(uc)情况，
			采用Chin-Chia Jane and Yih-Wenn Laih的方法划分*/
			memset(x0,0,(nE+1)*sizeof(int));                                       /*初始化划分线x0*/
			ComputeX0(g,Fd,x0);                                                    /*获得划分线*/
			GetNextCollection(Cj,x0,next_c);
			cStack.push(next_c); 
			GetCurrentC0(next_c,C0);                                               /*通过划分知：C0必定能够取得最大流*/
			SAVEM(C0);                                                             /*直接考虑C0下界子图*/
		}
	}

	/*在最终的满足最大流且去掉任意一条边后
	都小于最大流的概率最大下界上运行最大流算法*/
	FromG2G(g,cur_g,lofMaxP); 

	maxflow = Dinic(cur_g,source,sink,resultFd,gf);                                /*resultFd保存最终结果*/

	delete[] x0;                                                                   /*释放申请的空间*/
	delete[] lofMaxP; 

	return dR;
}

/************************************************************************/
/* 计算子图区间概率之积*/
/* g为不确定图，c为子图区间
/************************************************************************/
double compute_all_subgraph(Graph &g, StateSet c)
{
	double p =1;
	for (int i=1; i<= c.numE; i++)
	{
		if (c.lower[i]==1 && c.upper[i]==1)
		{
			p *= g.AllEdge_p[i];
		}
		else if (c.lower[i]==0 && c.upper[i]==0)
		{
			p *= (1-g.AllEdge_p[i]);
		}
		else if (c.lower[i]==0 && c.upper[i]==1)
		{
			p *= 1;
		}
	}
	return p;
}



/*通过对可能事件模型进行划分得到最可靠最大流分布*/
double  GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd,Lower_subGraph * StateMtrix)
{	
	assert((source > 0 && source <= g.nV) && (sink > 0 && sink <= g.nV));         /*保证输入合法*/
	
	int nE = g.nE;                                                                /*向量大小是不变的*/
	double dR = 0.0;                                                              /*最可靠最大流分布概率*/

	Flow Fd;  
	GF gf;                                                                        /*剩余图*/
	int Fmax = Dinic(g,source,sink,Fd,gf);                                        /*得到所有边都存在的最大流*/
	assert(Fmax >= 0);                                                             /*保证最大流大于0有意义*/

	/*保存随机流网络可靠性，所有满足最大流的子概率之和，初始化为1*/
	double max_p2 = 0;

	/*完成一个六元组的初始化*/
	Collection c(nE);
	stack<Collection> cStack;                                                     /*存放需要进一步划分的子图空间*/
	cStack.push(c);
	
	Graph cur_g;
	StateSet Cj(nE),C0(nE);                                                       /*存放当前需要处理的闭合区间[]*/
	Collection next_c(nE);                                                        /*保存下一个需要进栈处理的collection*/

	int tmpF = 0;
	double tmpP = 0.0;
	int * lofMaxP = new int[nE+1];                                                 /*临时存储满足要求的闭合区间中概率最大的下界*/
	int * x0 = new int[nE+1];                                                      /*划分线*/
	memset(lofMaxP,0,(nE+1)*sizeof(int));
	memset(x0,0,(nE+1)*sizeof(int));
	while(!cStack.empty())
	{	
		GetCurrentC_j(cStack.top(),Cj);                                            /*取得当前需要处理的闭合区间[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*变换成下一个待处理的状态*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G(g,cur_g,Cj.lower);                                                 /*下界向量对应的子图*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf); 
		if(tmpF >= Fmax)                                                           /*完备区间下界对应子图能够满足最大流*/
		{
			//SAVEM(Cj);
			/*以上语句代替下一段*/
			tmpP = CalculateP(g,(Cj).lower); 
			if(tmpP > dR)
			{
				dR = tmpP; 
				memcpy(lofMaxP,(Cj).lower,(nE+1)*sizeof(int));
				for (int ii =1; ii<= g.nE; ii++)
				{   //更新并保存极小子图
					StateMtrix->State[0][ii] = Cj.lower[ii];
				}
			} 
			//尝试输出下界子图
			
			for (int ii=1;ii<=Cj.numE;ii++) {
				cout<<Cj.lower[ii]; }
			cout<<"--Cj-->";
			for (int ii=1;ii<=Cj.numE;ii++) {
				cout<<Cj.upper[ii]; }
			cout<<endl;
			
			//保存子图区间
			saveAllState(Cj,StateMtrix);
			//加上该子图区间的所有子图的概率之和
			max_p2 += compute_all_subgraph(g, Cj);
		}  
		else if((FromG2G(g,cur_g,Cj.upper),                                        /*对上界对应的子图求最大流*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*同时满足F(lc) < Fmax <= F(uc)情况，
			采用Chin-Chia Jane and Yih-Wenn Laih的方法划分*/
			memset(x0,0,(nE+1)*sizeof(int));                                       /*初始化划分线x0*/
			ComputeX0(g,Fd,x0);                                                    /*获得划分线*/
			GetNextCollection(Cj,x0,next_c);
			cStack.push(next_c); 
			GetCurrentC0(next_c,C0);                                               /*通过划分知：C0必定能够取得最大流*/
			//SAVEM(C0);                                                             /*直接考虑C0下界子图*/
			/*以上语句代替下一段*/
			tmpP = CalculateP(g,(C0).lower); 
			if(tmpP > dR)
			{
					dR = tmpP; 
					memcpy(lofMaxP,(C0).lower,(nE+1)*sizeof(int));
					for (int ii =1; ii<= g.nE; ii++)
					{   //更新并保存极小子图
						StateMtrix->State[0][ii] = C0.lower[ii];
					}
			} 

			//尝试输出下界子图
			for (int ii=1;ii<=C0.numE;ii++) {
				cout<<C0.lower[ii]; }
			cout<<"--C0-->";
			for (int ii=1;ii<=C0.numE;ii++) {
				cout<<C0.upper[ii]; }
			cout<<endl; 
			//保存子图区间
			saveAllState(C0,StateMtrix);
			//加上该子图区间的所有子图的概率之和，求容量可靠性
			max_p2 += compute_all_subgraph(g, C0);
		}
	}
	cout<<endl;

	/*在最终的满足最大流且去掉任意一条边后
	都小于最大流的概率最大下界上运行最大流算法*/
	FromG2G(g,cur_g,lofMaxP); 

	maxflow = Dinic(cur_g,source,sink,resultFd,gf);                                /*resultFd保存最终结果*/
	//保存最大流
	g.max_flow = maxflow;
	/*在图数据中保存原始的  容量可靠性   和  分布可靠性*/
	g.max_p2 = max_p2;
	g.max_p1 = dR;

	delete[] x0;                                                                   /*释放申请的空间*/
	delete[] lofMaxP;
	return dR;
}

/************************************************************************/
/* 计算子图区间概率之积*/
/* g为不确定图，c为子图区间，i为故障的边
/************************************************************************/
double compute_subgraph(Graph &g, StateSet c, int ii)
{
	double p =1;
	for (int i=1; i<= c.numE; i++)
	{
		//如果不是断掉的边，断掉的边不需要计算
		if (ii != i)
		{
			if (c.lower[i]==1 && c.upper[i]==1)
			{
				p *= g.AllEdge_p[i];
			}
			else if (c.lower[i]==0 && c.upper[i]==0)
			{
				p *= (1-g.AllEdge_p[i]);
			}
			else if (c.lower[i]==0 && c.upper[i]==1)
			{
				p *= 1;
			}
		}
	}
	return p;
}

/*通过对可能事件模型进行划分得到最可靠最大流分布*/
double base_GetMPMF(Graph& g, int source, int sink, int &maxflow, Flow& resultFd, double &max_p2, int break_edge)
{	
	assert((source > 0 && source <= g.nV) && (sink > 0 && sink <= g.nV));         /*保证输入合法*/
	
	/*这里需要修改一下边的长度，因为发生故障之后边减少一个，但是依然有一个空位*/
	int nE = g.nE+1;                                                                /*向量大小是不变的*/
	double dR = 0.0;                                                              /*最可靠最大流分布概率*/

	Flow Fd;  
	GF gf;                                                                        /*剩余图*/
	int Fmax = Dinic(g,source,sink,Fd,gf);                                        /*得到所有边都存在的最大流*/
	assert(Fmax >= 0);                                                             /*保证最大流大于0有意义*/

	max_p2 = 0;//保存随机流网络可靠性，所有满足最大流的子概率之和，初始化为1

	/*完成一个六元组的初始化*/
	Collection c(nE);
	stack<Collection> cStack;                                                     /*存放需要进一步划分的子图空间*/
	cStack.push(c);
	
	Graph cur_g;
	StateSet Cj(nE),C0(nE);                                                       /*存放当前需要处理的闭合区间[]*/
	Collection next_c(nE);                                                        /*保存下一个需要进栈处理的collection*/

	int tmpF = 0;
	double tmpP = 0.0;
	int * lofMaxP = new int[nE+1];                                                 /*临时存储满足要求的闭合区间中概率最大的下界*/
	int * x0 = new int[nE+1];                                                      /*划分线*/
	memset(lofMaxP,0,(nE+1)*sizeof(int));
	memset(x0,0,(nE+1)*sizeof(int));
	while(!cStack.empty())
	{	
		GetCurrentC_j(cStack.top(),Cj);                                            /*取得当前需要处理的闭合区间[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*变换成下一个待处理的状态*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G(g,cur_g,Cj.lower);                                                 /*下界向量对应的子图*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf); 
		if(tmpF >= Fmax)                                                           /*完备区间下界对应子图能够满足最大流*/
		{
			SAVEM(Cj);
			//加上该子图区间的所有子图的概率之和
			max_p2 += compute_subgraph(g, Cj, break_edge);
		}  
		else if((FromG2G(g,cur_g,Cj.upper),                                        /*对上界对应的子图求最大流*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*同时满足F(lc) < Fmax <= F(uc)情况，
			采用Chin-Chia Jane and Yih-Wenn Laih的方法划分*/
			memset(x0,0,(nE+1)*sizeof(int));                                       /*初始化划分线x0*/
			ComputeX0(g,Fd,x0);                                                    /*获得划分线*/
			GetNextCollection(Cj,x0,next_c);
			cStack.push(next_c); 
			GetCurrentC0(next_c,C0);                                               /*通过划分知：C0必定能够取得最大流*/
			SAVEM(C0);                                                             /*直接考虑C0下界子图*/
			//加上该子图区间的所有子图的概率之和
			max_p2 += compute_subgraph(g, C0, break_edge);
		}
	}

	/*在最终的满足最大流且去掉任意一条边后
	都小于最大流的概率最大下界上运行最大流算法*/
	FromG2G(g,cur_g,lofMaxP); 

	maxflow = Dinic(cur_g,source,sink,resultFd,gf);                                /*resultFd保存最终结果*/

	delete[] x0;                                                                   /*释放申请的空间*/
	delete[] lofMaxP; 

	return dR;
}