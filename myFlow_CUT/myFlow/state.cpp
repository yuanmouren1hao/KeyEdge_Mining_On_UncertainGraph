#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"

#include <set>
#include <stack>
#include <vector>

using namespace std;

/*判断一个子图的最大流*/
/*输入：vector<int>表示的一个子图,removedEdge表示被移除的边，如果为0的话那么就表示没有移除，返回：该子图的最大流*/
int DinicByVertor(vector<int> a, Graph& g, int removedEdge, int source, int sink)
{
	/*通过vector构造图*/
	Graph new_g;
	new_g.nV = g.nV;
	if (0==removedEdge)
	{
		new_g.nE = g.nE;
	}else{
		new_g.nE = g.nE -1;
	}

	/*copy矩阵信息*/
	int aa,bb,cc,dd;
	double ff;
	for (int i=1; i<=g.nE; i++)
	{
		if (i!=removedEdge)
		{
			if (a[i] == 0)
			{
				new_g.nE--;
			}
			else if (a[i] == 1)
			{
				aa = g.AllEdge[i][1];//起始点
				bb = g.AllEdge[i][2];//终点
				cc = g.AllEdge[i][3];//容量
				dd = g.AllEdge[i][4];//有效位（此处不使用）

				ff = g.AllEdge_p[i];//可靠性

				new_g.matrix[aa][bb].iLabel = i;
				new_g.matrix[aa][bb].iC = cc;
				new_g.matrix[aa][bb].dP = ff;
			}
		}
	}

	/*使用原生Dinic算法计算最大流*/
	Flow f;
	GF gf;
	return Dinic(new_g, source, sink, f, gf);
}


/*采用Chin-Chia Jane and Yih-Wenn Laih的方法得到C(j)(1-q)*/
void GetCurrentC_j_new(Collection &c,StateSet &Cj)
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

/*将状态向量转换成对应图g_target*/
void FromG2G_new(Graph& g_original,Graph& g_target,int* lu)
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
void ComputeX0_new(Graph& g,Flow Fd,int* x0)
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

/*通过闭合区间与划分线得到pivot集合ad(c)*/
/*函数返回ad(c)集合的大小和Ad_c*/
int GetAd_C_new(StateSet &Cj,int *x0,int *Ad_c) 
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
void GetNextCollection_new(StateSet &Cj,int *x0,Collection &nextCollection)
{
	nextCollection.j = 1;                                                         /*设定需要处理要求的Cj*/
	nextCollection.I = GetAd_C_new(Cj,x0,nextCollection.Ad_C);                        /*设定I和Ad_c[]*/
	memcpy(nextCollection.lower,Cj.lower,(nextCollection.numE+1)*sizeof(int));    /*设定lower[]*/
	memcpy(nextCollection.upper,Cj.upper,(nextCollection.numE+1)*sizeof(int));    /*设定upper[]*/
}


/*采用Chin-Chia Jane and Yih-Wenn Laih方法获得C0*/
void GetCurrentC0_new(Collection &c,StateSet &C0)
{
	memcpy(C0.lower,c.lower,(C0.numE+1)*sizeof(int));
	memcpy(C0.upper,c.upper,(C0.numE+1)*sizeof(int));
	for(int i = 1; i <= c.I; i++)
	{
		C0.lower[c.Ad_C[i]] = C0.upper[c.Ad_C[i]];
	}
}

/*子图l被子图空间中子图包含的概率*/
double CalculateP_new(Graph& g,int* l)
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

vector<int> getVertorByC(StateSet C0)
{
	vector<int> v;
	v.push_back(0);

	//将区间缩写为一个vector
	for (int i=1; i<=C0.numE; i++)
	{
		if (C0.lower[i]==0 && C0.upper[i]==0)
		{
			v.push_back(0);
		}else if (C0.lower[i]==0 && C0.upper[i]==1)
		{
			v.push_back(2);
		}else if (C0.lower[i]==1 && C0.upper[i]==1)
		{
			v.push_back(1);
		}
	}
	return v;
}

/************************************************************************/
/* 根据状态空间的上下界进行划分，获取满足最大流的区间                        */
/************************************************************************/
/*输入是上届，下界，使用vector表示的，计算获取的是一个set表示的是分割之后的子图区间，
不存在的边一定要使用0代替，set的元素是set<vector<int>> */
void StateDivision(Graph& g, vector<int> down_, vector<int> top_, int source, int sink, set<vector<int>>& newState)
{	
	int nE = g.nE;                                                                /*向量大小是不变的*/
	double dR = 0.0;                                                              /*最可靠最大流分布概率*/

	Flow Fd;  
	GF gf;/*剩余图*/
	int Fmax = g.max_flow;/*得到所有边都存在的最大流*/
	assert(Fmax >= 0);/*保证最大流大于0有意义*/

	/*完成一个六元组的初始化*/
	Collection c(down_, top_, nE);
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
		GetCurrentC_j_new(cStack.top(),Cj);                                            /*取得当前需要处理的闭合区间[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*变换成下一个待处理的状态*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G_new(g,cur_g,Cj.lower);                                                 /*下界向量对应的子图*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf);
		if(tmpF >= Fmax)                                                           /*完备区间下界对应子图能够满足最大流*/
		{
			newState.insert(getVertorByC(Cj));
		}  
		else if((FromG2G_new(g,cur_g,Cj.upper),                                        /*对上界对应的子图求最大流*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*同时满足F(lc) < Fmax <= F(uc)情况，
			采用Chin-Chia Jane and Yih-Wenn Laih的方法划分*/
			memset(x0,0,(nE+1)*sizeof(int));   /*初始化划分线x0*/
			ComputeX0_new(g,Fd,x0);/*获得划分线*/
			GetNextCollection_new(Cj,x0,next_c);
			cStack.push(next_c);
			GetCurrentC0_new(next_c,C0);/*通过划分知：C0必定能够取得最大流*/
			//以下替代宏定义部分
			//将C0插入到新的状态区间中，这些是满足最大流的区间
			newState.insert(getVertorByC(C0));
		}
	}

	delete[] x0;                                                                   /*释放申请的空间*/
	delete[] lofMaxP;

	return;
}


/************************************************************************/
/* 上面为做实验部分                                                       */
/************************************************************************/
/*移除A类中的某一条边*/
void removeA_Edge(Graph& g,int i, Edge &TempEdge)
{
	int u=g.AllEdge[i][1];
	int v=g.AllEdge[i][2];
	TempEdge.dP=g.matrix[u][v].dP;
	TempEdge.iC=g.matrix[u][v].iC;
	TempEdge.iLabel=g.matrix[u][v].iLabel;
	g.nE--;
	g.matrix[u][v].dP=0;
	g.matrix[u][v].iC=0;
	g.matrix[u][v].iLabel=0;
	return;
};

/*恢复某一条A类边*/
void restoreA_Edge(Graph& g,int i, Edge &TempEdge)
{
	int u=g.AllEdge[i][1];
	int v=g.AllEdge[i][2];
	g.nE++;
	g.matrix[u][v].dP=TempEdge.dP;
	g.matrix[u][v].iC=TempEdge.iC;
	g.matrix[u][v].iLabel=TempEdge.iLabel;
	return;
}


/*初始化临时存储边*/
void init_TempEdge(Edge &TempEdge)
{
	TempEdge.dP=0;
	TempEdge.iC=0;
	TempEdge.iLabel=0;
	return;
}

/*计算一个A类边*/
void computetSimgleA(Graph& g, int source, int sink, int EdgeNum, KeyEdgeSet& AllKeyEdge)
{
	int new_maxflow = 0; /*新  最大流*/
	double new_p1 =0;/*新  分布可靠性*/
	double new_p2 =0;/*新  容量可靠性*/
	Flow new_maxPmaxF;/*新的最大流分布*/
	Edge TempEdge;/*临时存储一条边的信息*/

	/*初始化临时存储边*/
	init_TempEdge(TempEdge);
	/*去除某一条A类边*/
	removeA_Edge(g, EdgeNum, TempEdge);
	//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
	new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, EdgeNum);

	/* 将断掉边之后的最大流，分布可靠性，容量可靠性分别保存 */
	int now_edge_num = AllKeyEdge.A_num;
	AllKeyEdge.A_EdgeInfo[now_edge_num].ChangeAmount_c = new_maxflow;
	AllKeyEdge.A_EdgeInfo[now_edge_num].ChangeAmount_p1 = new_p1;
	AllKeyEdge.A_EdgeInfo[now_edge_num].ChangeAmount_p2 = new_p2;
	AllKeyEdge.A_EdgeInfo[now_edge_num].Edge = EdgeNum;
	AllKeyEdge.A_EdgeInfo[now_edge_num].Edge_Class = 'A';

	/* 恢复某一条A类边 */
	restoreA_Edge(g,EdgeNum, TempEdge);

	return;
}

/*调换连个边的位置*/
void Exchange_KeyEdge(KeyEdgeSet& AllKeyEdge,int i,int j)
{
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.ChangeAmount_p1 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1;
	TempKeyEdge.ChangeAmount_p2 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2;
	TempKeyEdge.Edge = AllKeyEdge.A_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class = AllKeyEdge.A_EdgeInfo[i].Edge_Class;

	AllKeyEdge.A_EdgeInfo[i].ChangeAmount_c = AllKeyEdge.A_EdgeInfo[j].ChangeAmount_c;
	AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1 = AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p1;
	AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2 = AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p2;
	AllKeyEdge.A_EdgeInfo[i].Edge = AllKeyEdge.A_EdgeInfo[j].Edge;
	AllKeyEdge.A_EdgeInfo[i].Edge_Class = AllKeyEdge.A_EdgeInfo[j].Edge_Class;

	AllKeyEdge.A_EdgeInfo[j].ChangeAmount_c = TempKeyEdge.ChangeAmount_c;
	AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p1 = TempKeyEdge.ChangeAmount_p1;
	AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p2 = TempKeyEdge.ChangeAmount_p2;
	AllKeyEdge.A_EdgeInfo[j].Edge = TempKeyEdge.Edge;
	AllKeyEdge.A_EdgeInfo[j].Edge_Class = TempKeyEdge.Edge_Class;

	return;
}

/*使用冒泡法排序一段区间,根据的是如果分布可靠性一致，而流量可靠性一致
/************************************************************************/
/* AllKeyEdge保存所有关键边，i为初始位置，j为结束位置                       */
/************************************************************************/
void Bubbling_KeyEdge_2(KeyEdgeSet& AllKeyEdge, int begin_, int end_)
{
	/*需要满足i<j, 将从i到j的区间排序*/
	if (begin_ > end_)
	{
		return;
	}

	for (int i=begin_; i<end_; i++)
	{
		for (int j=begin_+1; j<=end_; j++)
		{
			if (AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p2 < AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge(AllKeyEdge, i, j);
			}
		}
	}
	return;
}

/*排序一段A类边*/
void sortPartA(KeyEdgeSet& AllKeyEdge, int _from, int _to)
{
	/*需要满足i<j, 将从i到j的区间排序*/
	if (_from > _to)
	{
		return;
	}

	/*比较分布可靠性*/
	for (int i=_from; i<_to; i++)
	{
		for (int j=_from+1; j<=_to; j++)
		{
			if (AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p1 < AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1)
			{
				Exchange_KeyEdge(AllKeyEdge, i, j);
			}
		}
	}

	/*分布可靠性一致，比较容量可靠性*/
	int iii = _from;
	int jjj = iii+1;
	while(jjj <= _to)
	{
		if (AllKeyEdge.A_EdgeInfo[jjj].ChangeAmount_p1 != AllKeyEdge.A_EdgeInfo[iii].ChangeAmount_p1)
		{
			iii++;
			jjj++;
		}else
		{
			//找到相同变化的区间
			for (int kkk =jjj; kkk <= AllKeyEdge.EdgeNum; kkk++)
			{
				if (AllKeyEdge.A_EdgeInfo[kkk].ChangeAmount_p1 != AllKeyEdge.A_EdgeInfo[iii].ChangeAmount_p1)
				{
					jjj = kkk;
					break;
				}
			}
			//使用冒泡排序一段
			Bubbling_KeyEdge_2(AllKeyEdge, iii, jjj-1);
			iii=jjj;
			jjj=iii+1;
		}

	}

	return;
}

/*对A类边重新计算*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, int EdgeFlow[][2])
{
	/*如果没有A类边的个数为0*/
	if (key_edge_set->A_num==0)
	{
		return;
	}

	/*定义一个已经排好序的关键边组合*/
	KeyEdgeSet AllKeyEdge;
	AllKeyEdge.A_num = 0;
	AllKeyEdge.EdgeNum = key_edge_set->EdgeNum;
	int a_num = key_edge_set->A_num;
	
	/*当前已经计算边的个数*/
	int i=1;
	int j=i+1;
	while(i<=a_num)
	{
		/*如果是最后一个*/
		if (i == a_num)
		{
			/*i这个边不需要计算，直接比较出来*/
			AllKeyEdge.A_num++;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge_Class = 'A';
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge = EdgeFlow[i][0];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_c = EdgeFlow[i][1];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p1 = 0;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p2 = 0;
			break;
		}

		if (EdgeFlow[i][1] < EdgeFlow[j][1])
		{
			/*i这个边不需要计算，直接比较出来*/
			AllKeyEdge.A_num++;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge_Class = 'A';
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge = EdgeFlow[i][0];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_c = EdgeFlow[i][1];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p1 = 0;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p2 = 0;

			/*i和j的步长都+1*/
			i++;
			j++;

		}else if (EdgeFlow[i][1] = EdgeFlow[j][1])
		{
			/*至少有2个连续的重复，需要计算重复的   首先找到有几个重复的*/
			while(EdgeFlow[j][1] == EdgeFlow[i][1])
			{
				j++;
			}

			/*对于每一个连续的需要计算可靠性  重复计算*/
			for (int ii =i; ii<=j-1;ii++)
			{
				/*把一个A类边的信息计算出来， 保存到AllKeyEdge中*/
				AllKeyEdge.A_num++;
				computetSimgleA(g, source, sink, EdgeFlow[ii][0], AllKeyEdge);
			}
			/*计算完成之后排序这一段,保存在AllKeyEdge*/
			sortPartA(AllKeyEdge, i, j-1);

			/*增加步长*/
			i=j;
			j++;
		}
	}

	/*将AllKeyEdge 的A部分保存到key_edge_set中*/
	for (int i=1; i<=AllKeyEdge.A_num; i++)
	{
		key_edge_set->A_EdgeInfo[i].ChangeAmount_c = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_c;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p1 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2;
		key_edge_set->A_EdgeInfo[i].Edge_Class = AllKeyEdge.A_EdgeInfo[i].Edge_Class;
		key_edge_set->A_EdgeInfo[i].Edge = AllKeyEdge.A_EdgeInfo[i].Edge;
	}

	return;
}


/////////////////////////////////////////////////////////////////////////
//所有关于B类边的操作

/*计算子图概率*/
/*输入：vector<int>表示的一个子图,remvedEdge表示被删除的边，如果为，返回：该子图的概率*/
double computeSubgraphProbabilityByVector(vector<int> a, int remvedEdge, Graph &g)
{
	double p = 1;
	for (int i=1; i<=g.nE; i++)
	{
		if (i!=remvedEdge)
		{
			if (a[i] == 0)
			{
				p *= 1-g.AllEdge_p[i];
			}else if (a[i] == 1)
			{
				p *= g.AllEdge_p[i];
			}else if (a[i] == 2)
			{
				p *= 1;
			}
		}
	}
	return p;
}


/*计算下界子图概率*/
/*输入：vector<int>表示的一个子图,remvedEdge表示被删除的边，如果为，返回：该子图的概率, 只使用下界*/
double computeSubgraphProbabilityByVector_lower(vector<int> a, int remvedEdge, Graph &g)
{
	double p = 1;
	for (int i=1; i<=g.nE; i++)
	{
		if (i!=remvedEdge)
		{
			if (a[i] == 1)
			{
				p *= g.AllEdge_p[i];
			}
		}
	}
	return p;
}

/*计算一些列子图的概率之和*/
/*输入：set<vector<int>>表示所有的子图集合,removedEdge表示需要移除的边，输出：所有子图概率之和*/
double computeProbabilitySum(set<vector<int>> new_set, Graph &g, int removedEdge)
{
	double p = 0;

	set<vector<int>>::iterator it;
	for (it=new_set.begin(); it!=new_set.end(); it++)
	{
		p += computeSubgraphProbabilityByVector(*it, 0, g);
	}
	return p;
}

/*比较两个偏序对的大小*/
/*b不大于a返回1，b>=a返回2*/
int compareSubgraph(vector<int> a, vector<int> b, int length, int removeEdge)
{
	for (int i=1; i<=length; i++)
	{
		if (i!=removeEdge)
		{
			if (b[i]<a[i])
			{
				return 1;
			}
		}
	}
	return 2;
}

/*根据所有的状态区间获取所有的满足最大流的状态*/
set<vector<int>> getAllState(Lower_subGraph* StateMtrix)
{
	/*将矩阵中的区间放到栈中*/
	stack<vector<int>> s;

	for (int i=1; i<= StateMtrix->State_Num; i++)
	{
		/*定义一个子图区间*/
		vector<int> v;
		v.push_back(0);/*为了填充第0个元素，从第一个元素开始使用*/

		for (int j=1; j<=StateMtrix->Edge_Num; j++)
		{
			if (StateMtrix->State[i][j]==0 && StateMtrix->State_upper[i][j]==0)
			{
				v.push_back(0);
			}else if (StateMtrix->State[i][j]==1 && StateMtrix->State_upper[i][j]==1){
				v.push_back(1);
			}
			else if (StateMtrix->State[i][j]==0 && StateMtrix->State_upper[i][j]==1){
				v.push_back(2);
			}
		}
		s.push(v);
	}

	/*保存所有状态*/
	set<vector<int>> allState ;
	/*将栈中的状态进行分解*/
	while (!s.empty())
	{
		vector<int> a = s.top();
		s.pop();
		
		for (int i=1; i<=StateMtrix->Edge_Num; i++)
		{
			if (a[i]==2)
			{
				/*如果出现X，也就是2，将子图区间分解，分别插入的栈中*/
				vector<int> temp_V = a;
				temp_V[i] =0;
				s.push(temp_V);
				temp_V[i] =1;
				s.push(temp_V);
				break;
			}
			
			/*如果直接没有的话，直接插入到 集合中*/
			if (i==StateMtrix->Edge_Num)
			{
				allState.insert(a);		
			}
		}
	}
	return allState;
}


/*获取所有满足最大流的子图区间中的极小子图的概率*/
double getP1ByStateSection(set<vector<int>>& allStateSection,Graph& g)
{
	double p1 = 0;

	set<vector<int>>::iterator it;
	double temp_p;
	for (it=allStateSection.begin(); it!=allStateSection.end(); it++)
	{
		temp_p= computeSubgraphProbabilityByVector_lower(*it, 0, g);
		if (temp_p > p1)
		{
			p1 = temp_p;
		}
	}

	return p1;
}

/*B类边增量计算，变量都使用引用的方式*/
void sortKeyEdge_B(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	/*如果没有B类边的个数为0*/
	if (key_edge_set->B1_num == 0)
	{
		return;
	}

	/*获取所有的满足最大流的状态的子图区间，注意是子图区间，2表示x，即0或者1*/
	set<vector<int>> allStateSection;

	for (int i =1; i<= key_edge_set->B1_num; i++)
	{
		/*断掉B1类边最大流不变*/
		key_edge_set->B1_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		/*边号*/
		int edge_num = key_edge_set->B1_EdgeInfo[i].Edge;

		/*将所有可以满足最大流区间的子图保存在allStateSection里面*/
		//遍历原来的满足最大流的子图区间
		for (int ii=1; ii<= StateMtrix->State_Num; ii++)
		{
			//如果在断掉的位置为0或则2的话，那么整个区间都满足最大流
			if (StateMtrix->State[ii][edge_num]==0)
			{
				/*定义一个子图区间*/
				vector<int> v;
				v.push_back(0);/*为了填充第0个元素，从第一个元素开始使用*/
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					//如果是断掉的边的话，直接保存为0
					if (j==edge_num)
					{
						v.push_back(0);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
					{
						v.push_back(0);
					}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
						v.push_back(1);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
					{
						v.push_back(2);
					}
				}
				allStateSection.insert(v);
			}
			//该区间部分满足最大流，需要进一步的判断
			else if (StateMtrix->State[ii][edge_num]==1)
			{
				//下界
				vector<int>TempV_low;
				TempV_low.push_back(0);
				//如果断掉之后 下界依然满足最大流，那么整个区间都满足
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_low.push_back(0);
					}
					else
					{
						TempV_low.push_back(StateMtrix->State[ii][j]);
					}
				}
				//上届
				vector<int>TempV_upper;
				TempV_upper.push_back(0);
				//如果断掉之后，区间上届不满足最大流，那么整个区间都不满足
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_upper.push_back(0);
					}
					else
					{
						TempV_upper.push_back(StateMtrix->State_upper[ii][j]);
					}
				}

				//如果下界满足最大流,将该区间直接保存在allStateSection中
				if (DinicByVertor(TempV_low,g,edge_num,source,sink) == g.max_flow)
				{
					vector<int> v;
					v.push_back(0);/*为了填充第0个元素，从第一个元素开始使用*/
					for (int j=1; j<=StateMtrix->Edge_Num; j++)
					{
						//如果是断掉的边的话，直接保存为0
						if (j==edge_num)
						{
							v.push_back(0);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
						{
							v.push_back(0);
						}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
							v.push_back(1);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
						{
							v.push_back(2);
						}
					}
					allStateSection.insert(v);
				}
				//如果上界子图不满足最大流，那么整个区间都不满足最大流，直接舍弃即可
				else if (DinicByVertor(TempV_upper,g,edge_num,source,sink) < g.max_flow)
				{
					//do nothing
					//cout<<"can subgraph allowed."<<endl;
				}
				//下界不满足最大流，上界满足最大流的话，需要重新分解子图区间，但是分解的子图区间比较小而已
				else
				{
					StateDivision(g,TempV_low,TempV_upper,source,sink,allStateSection);
				}
			}
		}

		//找到所有了满足最大流区间的话，直接获取极小子图和最大流可靠性即可
		key_edge_set->B1_EdgeInfo[i].ChangeAmount_p1 = getP1ByStateSection(allStateSection, g);
		/*最大流可靠性（容量可靠性）为所有满足最大流的概率之和*/
		key_edge_set->B1_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(allStateSection, g, edge_num);
	}
	return;
}

/************************************************************************/
/*for edge b  属于B2类边，断掉之后流量和分布可靠性不变  ，不在极小子图上     */
/************************************************************************/
/*B类边增量计算，变量都使用引用的方式*/
void sortKeyEdge_b(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	/*如果没有B2类边的个数为0*/
	if (key_edge_set->B2_num == 0)
	{
		return;
	}

	/*获取所有的满足最大流的状态*/
	set<vector<int>> allStateSection ;

	for (int i =1; i<= key_edge_set->B2_num; i++)
	{
		/*断掉B类边最大流不变，分布可靠性不变*/
		key_edge_set->B2_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		key_edge_set->B2_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
		/*记录边的序号*/
		int edge_num = key_edge_set->B2_EdgeInfo[i].Edge;


		/*将所有可以满足最大流区间的子图保存在allStateSection里面*/
		//遍历原来的满足最大流的子图区间
		for (int ii=1; ii<= StateMtrix->State_Num; ii++)
		{
			//如果在断掉的位置为0或则2的话，那么整个区间都满足最大流
			if (StateMtrix->State[ii][edge_num]==0)
			{
				/*定义一个子图区间*/
				vector<int> v;
				v.push_back(0);/*为了填充第0个元素，从第一个元素开始使用*/
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					//如果是断掉的边的话，直接保存为0
					if (j==edge_num)
					{
						v.push_back(0);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
					{
						v.push_back(0);
					}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
						v.push_back(1);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
					{
						v.push_back(2);
					}
				}
				allStateSection.insert(v);
			}
			//该区间部分满足最大流，需要进一步的判断
			else if (StateMtrix->State[ii][edge_num]==1)
			{
				//下界
				vector<int>TempV_low;
				TempV_low.push_back(0);
				//如果断掉之后 下界依然满足最大流，那么整个区间都满足
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_low.push_back(0);
					}
					else
					{
						TempV_low.push_back(StateMtrix->State[ii][j]);
					}
				}
				//上届
				vector<int>TempV_upper;
				TempV_upper.push_back(0);
				//如果断掉之后，区间上届不满足最大流，那么整个区间都不满足
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_upper.push_back(0);
					}
					else
					{
						TempV_upper.push_back(StateMtrix->State_upper[ii][j]);
					}
				}

				//如果下界满足最大流,将该区间直接保存在allStateSection中
				if (DinicByVertor(TempV_low,g,edge_num,source,sink) == g.max_flow)
				{
					vector<int> v;
					v.push_back(0);/*为了填充第0个元素，从第一个元素开始使用*/
					for (int j=1; j<=StateMtrix->Edge_Num; j++)
					{
						//如果是断掉的边的话，直接保存为0
						if (j==edge_num)
						{
							v.push_back(0);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
						{
							v.push_back(0);
						}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
							v.push_back(1);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
						{
							v.push_back(2);
						}
					}
					allStateSection.insert(v);
				}
				//如果上界子图不满足最大流，那么整个区间都不满足最大流，直接舍弃即可
				else if (DinicByVertor(TempV_upper,g,edge_num,source,sink) < g.max_flow)
				{
					//do nothing
					//cout<<"can subgraph allowed."<<endl;
				}
				//下界不满足最大流，上界满足最大流的话，需要重新分解子图区间，但是分解的子图区间比较小而已
				else
				{
					StateDivision(g,TempV_low,TempV_upper,source,sink,allStateSection);
				}
			}
		}
		/*最大流可靠性（容量可靠性）为所有满足最大流的概率之和*/
		key_edge_set->B2_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(allStateSection, g, edge_num);
	}

	return;
}



//////////////////////////////////////////////////////////////////
//C类边的操作

/*对于C类边进行计算,C类边不需要计算*/
void sortKeyEdge_C(KeyEdgeSet *key_edge_set,Graph &g)
{
	/*每一个C类边的发生故障之后的，状态不改变，直接保存即可*/
	for (int i=1; i <= key_edge_set->C_num; i++)
	{
		key_edge_set->C_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p2 = g.max_p2;
	}
	return;
}

/*对边进行的大分类*/
void sortKeyEdge(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set)
{
	//cout<<"StateMtrix->Edge_Num:"<<StateMtrix->Edge_Num<<endl;
	//cout<<"StateMtrix->State_Num:"<<StateMtrix->State_Num<<endl;
	/*初始化关键边集合*/
	key_edge_set->EdgeNum=StateMtrix->Edge_Num;
	int temp_sum=0;
	

	for (int ii=1;ii<=StateMtrix->Edge_Num;ii++)
	{
		temp_sum=0;
		/*通过求和的方式计算边的类别*/
		for (int jj=1;jj<=StateMtrix->State_Num;jj++)
		{
			temp_sum+=StateMtrix->State[jj][ii];
		}
		/*将关键边存储，ii表示的是边的编号*/
		if(temp_sum==StateMtrix->State_Num)
		{
			//cout<<ii<<" 是A类边。"<<endl;
			key_edge_set->A_num++;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge=ii;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge_Class='A';
		}
		else if (0==temp_sum)
		{
			//cout<<ii<<" 是C类边。"<<endl;
			key_edge_set->C_num++;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge=ii;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge_Class='C';
		}
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num && StateMtrix->State[0][ii]==1)
		{
			//cout<<ii<<" 是B类边"<<endl;
			key_edge_set->B1_num++;
			key_edge_set->B1_EdgeInfo[key_edge_set->B1_num].Edge=ii;
			key_edge_set->B1_EdgeInfo[key_edge_set->B1_num].Edge_Class='B';
		}
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num && StateMtrix->State[0][ii]==0)
		{
			//cout<<ii<<" 是B类边"<<endl;
			key_edge_set->B2_num++;
			key_edge_set->B2_EdgeInfo[key_edge_set->B2_num].Edge=ii;
			key_edge_set->B2_EdgeInfo[key_edge_set->B2_num].Edge_Class='b';
		}
	}
}

/************************************************************************/
/* 根据g和断掉的边获取，最大流                                             */
/************************************************************************/
int getMAXflowWithBreakEdge(Graph& g, int BreakEdge, int source, int sink)
{
	/*使用原来的图  先构造，获取最大流，然后恢复原有图的数据*/
	g.nE--;

	int aa,bb,cc,dd;//边BreakEdge的起始点，终点，流量，有效位
	double ff;//边BreakEdge的可靠性
	aa = g.AllEdge[BreakEdge][1];//起始点
	bb = g.AllEdge[BreakEdge][2];//终点
	cc = g.AllEdge[BreakEdge][3];//容量
	dd = g.AllEdge[BreakEdge][4];//有效位（此处不使用）
	ff = g.AllEdge_p[BreakEdge];//可靠性

	g.matrix[aa][bb].iLabel = 0;
	g.matrix[aa][bb].iC = 0;
	g.matrix[aa][bb].dP = 0;
	/*使用原生Dinic算法计算最大流*/
	Flow f;
	GF gf;
	int max_flow_remine = Dinic(g, source, sink, f, gf);

	/*计算完成之后，  需要恢复原有的图*/
	g.nE++;
	g.matrix[aa][bb].iLabel = BreakEdge;
	g.matrix[aa][bb].iC = cc;
	g.matrix[aa][bb].dP = ff;

	/*返回数据*/
	return max_flow_remine;
}


/************************************************************************/
/* 计算每条边断掉之后流量减少量                                            */
/************************************************************************/
void getALLEdgeFlowDecrease(Graph& g, int EdgeFlow[][2], int source, int sink)
{
	vector<int> v;
	v.push_back(0);
	for (int i=1; i<=g.nE; i++)
	{
		v.push_back(0);
	}

	for (int i=1; i<=g.nE; i++)
	{
		EdgeFlow[i][0] = i;
		EdgeFlow[i][1] = getMAXflowWithBreakEdge(g, i, source, sink);
	}

	/*将EdgeFlow进行排序*/
	for (int i =1; i<g.nE; i++)
	{
		for(int j=i+1; j<=g.nE; j++)
		{
			if (EdgeFlow[j][1] < EdgeFlow[i][1])
			{
				int temp;

				temp = EdgeFlow[j][1];
				EdgeFlow[j][1] = EdgeFlow[i][1];
				EdgeFlow[i][1] = temp;

				temp = EdgeFlow[j][0];
				EdgeFlow[j][0] = EdgeFlow[i][0];
				EdgeFlow[i][0] = temp;
			}
		}
	}

	return;
}



/*通过状态矩阵计算边的类别*/
void computeEdgeClass(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*先计算每条边断掉之后流量减少*/
	int EdgeFlow[MAX_E_NUM][2];
	getALLEdgeFlowDecrease(g, &EdgeFlow[0], source, sink);
	
	/*定性分析，将不确定图中的边分为A，B1，B2，C四类边*/
	sortKeyEdge(StateMtrix, key_edge_set);
	
	/*定量计算*/
	/*A类边重新计算*/
	sortKeyEdge_A(key_edge_set, g, source, sink, &EdgeFlow[0]);
	/*B类边增量计算*/
	sortKeyEdge_B(key_edge_set, g, source, sink, StateMtrix);
	/*for class edge b */
	sortKeyEdge_b(key_edge_set, g, source, sink, StateMtrix);
	/*C类边不需要计算*/
	sortKeyEdge_C(key_edge_set, g);
	
	return;
}














//////////////////////////////////////////////////////////////////
//下面是接口，不在该类中调用


/*将关键边输出*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"关键边输出如下："<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.A_num;i++)
	{
		out<<setw(5)<<key_edge_set.A_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.A_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.B1_num;i++)
	{
		out<<setw(5)<<key_edge_set.B1_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B1_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B1_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B1_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.B1_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.B2_num;i++)
	{
		out<<setw(5)<<key_edge_set.B2_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B2_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B2_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B2_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.B2_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.C_num;i++)
	{
		out<<setw(5)<<key_edge_set.C_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.C_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	out<<"---------------------------------"<<endl<<"---------------------------------"<<endl<<endl;
}


void init_KeyEdgeSet(KeyEdgeSet &key_edge_set,Lower_subGraph &StateMtrix)
{
	/*初始化关键边集合*/
	key_edge_set.A_num=0;
	key_edge_set.B1_num=0;
	key_edge_set.B2_num=0;
	key_edge_set.C_num=0;
	key_edge_set.EdgeNum=0;

	for (int i=0;i<MAX_E_NUM;i++)
	{
		key_edge_set.A_EdgeInfo[i].Edge=0;
		key_edge_set.B1_EdgeInfo[i].Edge=0;
		key_edge_set.B2_EdgeInfo[i].Edge=0;
		key_edge_set.C_EdgeInfo[i].Edge=0;

		key_edge_set.A_EdgeInfo[i].Edge_Class=0;
		key_edge_set.B1_EdgeInfo[i].Edge_Class=0;
		key_edge_set.B2_EdgeInfo[i].Edge_Class=0;
		key_edge_set.C_EdgeInfo[i].Edge_Class=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.B1_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.B2_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_c=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.B1_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.B2_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_p1 = 0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.B1_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.B2_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_p2 = 0;
	}
	/*初始化下界子图*/
	StateMtrix.Edge_Num=0;
	StateMtrix.State_Num=0;
	for (int i=0;i<NUM_OF_SUB_GRAPH;i++)
	{
		for (int j=0;j<MAX_E_NUM;j++)
		{
			StateMtrix.State[i][j]=0;
		}
	}
	return;
}