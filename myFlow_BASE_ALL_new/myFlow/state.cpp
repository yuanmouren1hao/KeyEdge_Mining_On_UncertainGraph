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


//////////////////////////////////////////////////////////////
//new
//////////////////////////////////////////////////////////////

/*调换2个边的位置*/
void Exchange_KeyEdge(KeyEdgeSet *key_edge_set,int i,int j)
{
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c  = key_edge_set->EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.ChangeAmount_p1 = key_edge_set->EdgeInfo[i].ChangeAmount_p1;
	TempKeyEdge.ChangeAmount_p2 = key_edge_set->EdgeInfo[i].ChangeAmount_p2;
	TempKeyEdge.Edge			= key_edge_set->EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class		= key_edge_set->EdgeInfo[i].Edge_Class;

	key_edge_set->EdgeInfo[i].ChangeAmount_c		=  key_edge_set->EdgeInfo[j].ChangeAmount_c	;
	key_edge_set->EdgeInfo[i].ChangeAmount_p1		=  key_edge_set->EdgeInfo[j].ChangeAmount_p1;
	key_edge_set->EdgeInfo[i].ChangeAmount_p2		=  key_edge_set->EdgeInfo[j].ChangeAmount_p2;
	key_edge_set->EdgeInfo[i].Edge					=  key_edge_set->EdgeInfo[j].Edge			;
	key_edge_set->EdgeInfo[i].Edge_Class			=  key_edge_set->EdgeInfo[j].Edge_Class		;
	
	key_edge_set->EdgeInfo[j].ChangeAmount_c		=  TempKeyEdge.ChangeAmount_c  ;
	key_edge_set->EdgeInfo[j].ChangeAmount_p1		=  TempKeyEdge.ChangeAmount_p1 ;
	key_edge_set->EdgeInfo[j].ChangeAmount_p2		=  TempKeyEdge.ChangeAmount_p2 ;
	key_edge_set->EdgeInfo[j].Edge					=  TempKeyEdge.Edge			   ;
	key_edge_set->EdgeInfo[j].Edge_Class			=  TempKeyEdge.Edge_Class	   ;

	return;
}

//排序一段p1
void Bubbling_KeyEdge_1(KeyEdgeSet *key_edge_set, int begin_, int end_)
{
	/*需要满足i<j, 将从i到j的区间排序*/
	if (begin_ > end_)
	{
		return;
	}

	for (int i=begin_; i<end_; i++)
	{
		for (int j=i+1; j<=end_; j++)
		{
			if (key_edge_set->EdgeInfo[j].ChangeAmount_p1 < key_edge_set->EdgeInfo[i].ChangeAmount_p1)
			{
				Exchange_KeyEdge(key_edge_set, i, j);
			}
		}
	}
	return;
}

//排序一段p2
void Bubbling_KeyEdge_2(KeyEdgeSet *key_edge_set, int begin_, int end_)
{
	/*需要满足i<j, 将从i到j的区间排序*/
	if (begin_ > end_)
	{
		return;
	}

	for (int i=begin_; i<end_; i++)
	{
		for (int j=i+1; j<=end_; j++)
		{
			if (key_edge_set->EdgeInfo[j].ChangeAmount_p2 < key_edge_set->EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge(key_edge_set, i, j);
			}
		}
	}
	return;
}


//排序，三重排序
void rankEdge(KeyEdgeSet *key_edge_set)
{
	/*首先按照从小到大排序最大流，ChangeAmount_c， 使用冒泡法排序*/
	for (int i =1; i< key_edge_set->EdgeNum; i++)
	{
		for (int j = i+1; j<= key_edge_set->EdgeNum; j++)
		{
			if (key_edge_set->EdgeInfo[j].ChangeAmount_c < key_edge_set->EdgeInfo[i].ChangeAmount_c)
			{
				Exchange_KeyEdge(key_edge_set, i,j);
			}
		}
	}

	/*在最大流ChangeAmount_c一致的情况下，排序分布可靠性ChangeAmount_p1*/
	int ii =1;
	int jj=ii+1;
	while(jj <= key_edge_set->EdgeNum)
	{
		if (key_edge_set->EdgeInfo[jj].ChangeAmount_c != key_edge_set->EdgeInfo[ii].ChangeAmount_c)
		{
			//前后不同
			ii++;
			jj++;
		}else
		{
			//找到相同变化的区间
			for (int k =jj; k <= key_edge_set->EdgeNum; k++)
			{
				if (key_edge_set->EdgeInfo[k].ChangeAmount_c != key_edge_set->EdgeInfo[ii].ChangeAmount_c )
				{
					jj = k;
					break;
				}
				if ( k == key_edge_set->EdgeNum)
				{
					jj = k+1;
					break;
				}
			}
			//使用冒泡排序一段ChangeAmount_p1
			Bubbling_KeyEdge_1(key_edge_set, ii, jj-1);
			ii=jj;
			jj=ii+1;
		}
	}

	/*在分布可靠性ChangeAmount_p1一致的情况下，排序容量可靠性ChangeAmount_p2*/
	int iii = 1;
	int jjj=iii+1;
	while(jjj<=key_edge_set->EdgeNum)
	{
		if (key_edge_set->EdgeInfo[jjj].ChangeAmount_c != key_edge_set->EdgeInfo[iii].ChangeAmount_c)
		{
			iii++;
			jjj++;
		}
		else if (key_edge_set->EdgeInfo[jjj].ChangeAmount_p1 != key_edge_set->EdgeInfo[iii].ChangeAmount_p1 && key_edge_set->EdgeInfo[jjj].ChangeAmount_c == key_edge_set->EdgeInfo[iii].ChangeAmount_c)
		{
			iii++;
			jjj++;
		}else
		{
			//找到相同变化的区间
			for (int kkk =jjj; kkk <= key_edge_set->EdgeNum; kkk++)
			{
				if (key_edge_set->EdgeInfo[kkk].ChangeAmount_p1 != key_edge_set->EdgeInfo[iii].ChangeAmount_p1)
				{
					jjj = kkk;
					break;
				}
				if ( kkk == key_edge_set->EdgeNum)
				{
					jjj = kkk+1;
					break;
				}
			}
			//使用冒泡排序一段
			Bubbling_KeyEdge_2(key_edge_set, iii, jjj-1);
			iii=jjj;
			jjj=iii+1;
		}
	}

	return;
}

////////////////////////////////////

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

/*移除A类中的某一条边*/
void remove_Edge(Graph& g,int i, Edge &TempEdge)
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

//base_all算法获取所有的最大流，分布可靠性，容量可靠性
void computeBaseAll(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	//对于其中的每一条边
	for (int i=1; i<=g.nE; i++)
	{
		int new_maxflow = 0; /*新  最大流*/
		double new_p1 =0;/*新  分布可靠性*/
		double new_p2 =0;/*新  容量可靠性*/
		Flow new_maxPmaxF;/*新的最大流分布*/
		Edge TempEdge;/*临时存储一条边的信息*/

		int EdgeNum = i;//移除的边号
		/*初始化临时存储边*/
		init_TempEdge(TempEdge);
		/*去除某一条A类边*/
		remove_Edge(g, EdgeNum, TempEdge);
		//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
		new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, EdgeNum);

		/* 将断掉边之后的最大流，分布可靠性，容量可靠性分别保存 */
		key_edge_set->EdgeNum++;
		key_edge_set->EdgeInfo[EdgeNum].ChangeAmount_c = new_maxflow;
		key_edge_set->EdgeInfo[EdgeNum].ChangeAmount_p1 = new_p1;
		key_edge_set->EdgeInfo[EdgeNum].ChangeAmount_p2 = new_p2;
		key_edge_set->EdgeInfo[EdgeNum].Edge = EdgeNum;
		key_edge_set->EdgeInfo[EdgeNum].Edge_Class = 'S';

		/* 恢复某一条A类边 */
		restoreA_Edge(g,EdgeNum, TempEdge);
	}

	//先对流量排序，之后相同的流量按照分布可靠性排序，相同的分布可靠性按照容量可靠性排序
	rankEdge(key_edge_set);

	return;
}

/*将关键边输出*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"关键边输出如下："<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.EdgeNum;i++)
	{
		out<<setw(5)<<key_edge_set.EdgeInfo[i].Edge<<setw(10)<<key_edge_set.EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	out<<"---------------------------------"<<endl<<"---------------------------------"<<endl<<endl;
}


void init_KeyEdgeSet(KeyEdgeSet &key_edge_set)
{
	/*初始化关键边集合*/
	key_edge_set.EdgeNum=0;
	for (int i=0;i<MAX_E_NUM;i++)
	{
		key_edge_set.EdgeInfo[i].Edge=0;
		key_edge_set.EdgeInfo[i].Edge_Class=0;
		key_edge_set.EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.EdgeInfo[i].ChangeAmount_p2 = 0;
	}

	return;
}