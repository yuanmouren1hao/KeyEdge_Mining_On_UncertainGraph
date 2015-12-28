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

/*对A类边重新计算*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*如果没有A类边的个数为0*/
	if (key_edge_set->A_num==0)
	{
		return;
	}

	int new_maxflow = 0; /*新  最大流*/
	double new_p1 =0;/*新  分布可靠性*/
	double new_p2 =0;/*新  容量可靠性*/
	Flow new_maxPmaxF;/*新的最大流分布*/
	Edge TempEdge;/*临时存储一条边的信息*/

	/*对于A类中的边，去掉一个*/
	for (int i = 1; i<=key_edge_set->A_num; i++)
	{
		/*初始化临时存储边*/
		init_TempEdge(TempEdge);
		/*去除某一条A类边*/
		removeA_Edge(g, key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
		//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
		new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->A_EdgeInfo[i].Edge);
		
		/* 将断掉边之后的最大流，分布可靠性，容量可靠性分别保存 */
		key_edge_set->A_EdgeInfo[i].ChangeAmount_c = new_maxflow;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p1 = new_p1;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = new_p2;

		/* 恢复某一条A类边 */
		restoreA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
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
		p += computeSubgraphProbabilityByVector(*it, removedEdge, g);
	}
	return p;
}


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

	/*copy举证信息*/
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




/*B类边增量计算，变量都使用引用的方式*/
void sortKeyEdge_B(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	/*如果没有B类边的个数为0*/
	if (key_edge_set->B_num == 0)
	{
		return;
	}

	/*极小子图*/
	/*原不确定图的极小子图保存在StateMtrix中的第0个子图下界中*/
	vector<int> minimumGraph;
	minimumGraph.push_back(0);/*从第一个元素开始使用*/
	for (int i=1; i<=StateMtrix->Edge_Num; i++)
	{
		minimumGraph.push_back(StateMtrix->State[0][i]);
	}

	/*获取所有的满足最大流的状态*/
	set<vector<int>> allState = getAllState(StateMtrix);

	for (int i =1; i<= key_edge_set->B_num; i++)
	{
		/*断掉B类边最大流不变*/
		key_edge_set->B_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		/*获取参与过滤的下界子图*/
		vector<int> lowerGraph;

		/*判断该边是否在极小子图中*/
		if (minimumGraph[key_edge_set->B_EdgeInfo[i].Edge]==0)
		{/*如果该边不在极小子图中*/
			
			/*分布可靠性不变*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
			/*需要过滤的图就是极小子图*/
			lowerGraph = minimumGraph;

			set<vector<int>> new_set;
			set<vector<int>>::iterator it;
			for (it = allState.begin(); it!=allState.end(); it++)
			{
				vector<int> temp_V ;
				if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==2)
				{
					/*大于筛选子图能满足最大流*/
					temp_V = *it;
					temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
					new_set.insert(temp_V);
				}else if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==1)
				{
					/*不大于表示不能确定是否满足，此时需要计算*/
					if (DinicByVertor(*it, g, key_edge_set->B_EdgeInfo[i].Edge, source, sink)==g.max_flow)
					{
						/*如果该子图能够满足最大流：保存，否则：不保存*/
						temp_V = *it;
						temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
						new_set.insert(temp_V);
					}
				}
			}
			
			/*最大流可靠性（容量可靠性）为所有满足最大流的概率之和*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(new_set, g, key_edge_set->B_EdgeInfo[i].Edge);

		}else if(minimumGraph[key_edge_set->B_EdgeInfo[i].Edge]==1)
		{/*该边在极小子图中*/

			/*记录最大分布可靠性*/
			double max_p1 =0;

			/*找到不含有该边的其他下界子图，最为过滤子图*/
			for (int ii=1; ii<= StateMtrix->State_Num; ii++)
			{
				if (StateMtrix->State[ii][key_edge_set->B_EdgeInfo[i].Edge] == 0)
				{
					/*该子图可用*/
					lowerGraph.push_back(0);
					for (int jj=1; jj<=StateMtrix->Edge_Num; jj++)
					{
						/*作为过滤子图*/
						lowerGraph.push_back(StateMtrix->State[ii][jj]);
					}
					break;
				}
			}

			/*将过滤子图作为极小子图，后续更新即可*/
			max_p1 = computeSubgraphProbabilityByVector(lowerGraph, key_edge_set->B_EdgeInfo[i].Edge, g);

			set<vector<int>> new_set;
			set<vector<int>>::iterator it;
			for (it = allState.begin(); it!=allState.end(); it++)
			{
				vector<int> temp_V ;
				if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==2)
				{
					/*大于筛选子图能满足最大流*/
					temp_V = *it;
					temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
					new_set.insert(temp_V);

					/*判断是为极小子图*/
					if (computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g)>max_p1)
					{
						lowerGraph = temp_V;
						max_p1 = computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g);
					}

				}else if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==1)
				{
					/*不大于表示不能确定是否满足，此时需要计算*/
					int new_maxflow = DinicByVertor(*it, g, key_edge_set->B_EdgeInfo[i].Edge, source, sink);
					if (new_maxflow==g.max_flow)
					{
						/*如果该子图能够满足最大流：保存，否则：不保存*/
						temp_V = *it;
						temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
						new_set.insert(temp_V);

						/*判断是否为极小子图*/
						if (computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g)>max_p1)
						{
							lowerGraph = temp_V;
							max_p1 = computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g);
						}
					}
					
				}
			}

			/*分布可靠性为不断更新的值*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p1 = max_p1;
			/*最大流可靠性（容量可靠性）为所有满足最大流的概率之和*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(new_set, g, key_edge_set->B_EdgeInfo[i].Edge);

		/*该边在极小子图中*/
		}

	/*对于所有的B类边循环结束*/
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
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num)
		{
			//cout<<ii<<" 是B类边"<<endl;
			key_edge_set->B_num++;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge=ii;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge_Class='B';
		}
	}
}

/*通过状态矩阵计算边的类别*/
void computeEdgeClass(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*定性分析*/
	sortKeyEdge(StateMtrix, key_edge_set);
	
	/*定量计算*/
	/*A类边重新计算*/
	sortKeyEdge_A(key_edge_set, g, source, sink);
	/*B类边增量计算*/
	sortKeyEdge_B(key_edge_set, g, source, sink, StateMtrix);
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
	for (int i=1;i<=key_edge_set.B_num;i++)
	{
		out<<setw(5)<<key_edge_set.B_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_p2<<endl;
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
	key_edge_set.B_num=0;
	key_edge_set.C_num=0;
	key_edge_set.EdgeNum=0;
	for (int i=0;i<MAX_E_NUM;i++)
	{
		key_edge_set.A_EdgeInfo[i].Edge=0;
		key_edge_set.B_EdgeInfo[i].Edge=0;
		key_edge_set.C_EdgeInfo[i].Edge=0;

		key_edge_set.A_EdgeInfo[i].Edge_Class=0;
		key_edge_set.B_EdgeInfo[i].Edge_Class=0;
		key_edge_set.C_EdgeInfo[i].Edge_Class=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_c=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_p1 = 0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount_p2 = 0;
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