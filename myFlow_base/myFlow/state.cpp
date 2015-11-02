#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"

/*引入外部的原始的计算最可靠最大流分布的算法*/
//#include "../Base/partition.cpp"


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

/*调换A类边的两个边的位置*/
/*
void Exchange_A(KeyEdgeSet *key_edge_set,int i,int j)
{
	//cout<<i<<"--->"<<j<<endl;
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c=key_edge_set->A_EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.Edge=key_edge_set->A_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class=key_edge_set->A_EdgeInfo[i].Edge_Class;

	key_edge_set->A_EdgeInfo[i].ChangeAmount_c=key_edge_set->A_EdgeInfo[j].ChangeAmount_c;
	key_edge_set->A_EdgeInfo[i].Edge=key_edge_set->A_EdgeInfo[j].Edge;
	key_edge_set->A_EdgeInfo[i].Edge_Class=key_edge_set->A_EdgeInfo[j].Edge_Class;

	key_edge_set->A_EdgeInfo[j].ChangeAmount_c=TempKeyEdge.ChangeAmount_c;
	key_edge_set->A_EdgeInfo[j].Edge=TempKeyEdge.Edge;
	key_edge_set->A_EdgeInfo[j].Edge_Class=TempKeyEdge.Edge_Class;

	return;
}
*/



/*对于A类边进行排序*/
void sortKeyEdge_ABC(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*测试输出*/
	cout<<"sortKeyEdge_ABC"<<endl;

	int new_maxflow = 0; /*保存新最大流*/
	double new_p1 =0;/*新的最大可靠性(最可靠最大流分布)*/
	double new_p2 =0;/*新的网络可靠性*/
	Flow new_maxPmaxF;/*新的最大流分布*/
	Edge TempEdge;/*临时存储一条边的信息*/


	/*如果有A类边的话*/
	if (key_edge_set->A_num > 0)
	{
		/*对于A类中的边，去掉一个*/
		for (int i=1; i<=key_edge_set->A_num; i++)
		{
			/*初始化临时存储边*/
			init_TempEdge(TempEdge);
			new_maxflow = 0;
			new_p1 = 0;
			new_p2 = 0;

			/*去除某一条A类边*/
			removeA_Edge(g, key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
			//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
			new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->A_EdgeInfo[i].Edge);

			/* 将改变量保存在A边信息,保存的是断掉边之后能够达到的最大流 */
			key_edge_set->A_EdgeInfo[i].ChangeAmount_c = new_maxflow;
			key_edge_set->A_EdgeInfo[i].ChangeAmount_p1 = new_p1;
			key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = new_p2;

			/* 恢复某一条A类边 */
			restoreA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
		}
	}

	/*如果有B类边的话*/
	if (key_edge_set->B_num > 0)
	{
		/*对于A类中的边，去掉一个*/
		for (int i=1; i<=key_edge_set->B_num; i++)
		{
			/*初始化临时存储边*/
			init_TempEdge(TempEdge);
			new_maxflow = 0;
			new_p1 = 0;
			new_p2 = 0;

			/*去除某一条A类边*/
			removeA_Edge(g, key_edge_set->B_EdgeInfo[i].Edge, TempEdge);
			//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
			new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->B_EdgeInfo[i].Edge);

			/* 将改变量保存在A边信息,保存的是断掉边之后能够达到的最大流 */
			key_edge_set->B_EdgeInfo[i].ChangeAmount_c = new_maxflow;
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p1 = new_p1;
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = new_p2;

			/* 恢复某一条A类边 */
			restoreA_Edge(g,key_edge_set->B_EdgeInfo[i].Edge, TempEdge);
		}
	}

	/*如果有C类边的话*/
	if (key_edge_set->C_num > 0)
	{
		/*对于A类中的边，去掉一个*/
		for (int i=1; i<=key_edge_set->C_num; i++)
		{
			/*初始化临时存储边*/
			init_TempEdge(TempEdge);
			new_maxflow = 0;
			new_p1 = 0;
			new_p2 = 0;

			/*去除某一条A类边*/
			removeA_Edge(g, key_edge_set->C_EdgeInfo[i].Edge, TempEdge);
			//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
			new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->C_EdgeInfo[i].Edge);

			/* 将改变量保存在A边信息,保存的是断掉边之后能够达到的最大流 */
			key_edge_set->C_EdgeInfo[i].ChangeAmount_c = new_maxflow;
			key_edge_set->C_EdgeInfo[i].ChangeAmount_p1 = new_p1;
			key_edge_set->C_EdgeInfo[i].ChangeAmount_p2 = new_p2;

			/* 恢复某一条A类边 */
			restoreA_Edge(g,key_edge_set->C_EdgeInfo[i].Edge, TempEdge);
		}
	}
	
	return;
}

/*初始化*/
void init_TempP(double TempP[])
{
	for (int i=0;i<MAX_E_NUM;i++)
	{
		TempP[i]=0.0;
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
	/*定性分类*/
	sortKeyEdge(StateMtrix, key_edge_set);
	
	/*对关键边进行定量计算*/
	sortKeyEdge_ABC(key_edge_set, g, source, sink);
	/*测试输出*/
	
	for (int i=1; i<=key_edge_set->A_num; i++)
	{
		cout<<setw(5)<<key_edge_set->A_EdgeInfo[i].Edge<<setw(10)<<key_edge_set->A_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set->A_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set->A_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set->A_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1; i<=key_edge_set->B_num; i++)
	{
		cout<<setw(5)<<key_edge_set->B_EdgeInfo[i].Edge<<setw(10)<<key_edge_set->B_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set->B_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set->B_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set->B_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1; i<=key_edge_set->C_num; i++)
	{
		cout<<setw(5)<<key_edge_set->C_EdgeInfo[i].Edge<<setw(10)<<key_edge_set->C_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set->C_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set->C_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set->C_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	
	cout<<"==========================="<<endl<<endl<<endl;
	return;
}


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