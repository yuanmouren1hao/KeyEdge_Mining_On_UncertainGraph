#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"

#include <windows.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////
//1.

/*调换2个边的位置*/
void Exchange_KeyEdge(KeyEdgeSet *key_edge_set,int i,int j)
{
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c  = key_edge_set->EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.ChangeAmount_p2 = key_edge_set->EdgeInfo[i].ChangeAmount_p2;
	TempKeyEdge.Edge			= key_edge_set->EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class		= key_edge_set->EdgeInfo[i].Edge_Class;

	key_edge_set->EdgeInfo[i].ChangeAmount_c		=  key_edge_set->EdgeInfo[j].ChangeAmount_c	;
	key_edge_set->EdgeInfo[i].ChangeAmount_p2		=  key_edge_set->EdgeInfo[j].ChangeAmount_p2;
	key_edge_set->EdgeInfo[i].Edge					=  key_edge_set->EdgeInfo[j].Edge			;
	key_edge_set->EdgeInfo[i].Edge_Class			=  key_edge_set->EdgeInfo[j].Edge_Class		;

	key_edge_set->EdgeInfo[j].ChangeAmount_c		=  TempKeyEdge.ChangeAmount_c  ;
	key_edge_set->EdgeInfo[j].ChangeAmount_p2		=  TempKeyEdge.ChangeAmount_p2 ;
	key_edge_set->EdgeInfo[j].Edge					=  TempKeyEdge.Edge			   ;
	key_edge_set->EdgeInfo[j].Edge_Class			=  TempKeyEdge.Edge_Class	   ;

	return;
}

//计算一条边断掉之后
int getMAXflowWithBreakEdge(Graph& g, int BreakEdge, int source, int sink)
{
	//使用原来的图  先构造，获取最大流，然后恢复原有图的数据
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
	//使用原生Dinic算法计算最大流
	Flow f;
	GF gf;
	int max_flow_remine = Dinic(g, source, sink, f, gf);

	//计算完成之后，  需要恢复原有的图
	g.nE++;
	g.matrix[aa][bb].iLabel = BreakEdge;
	g.matrix[aa][bb].iC = cc;
	g.matrix[aa][bb].dP = ff;

	//返回数据
	return max_flow_remine;
}

//计算每条边断掉之后的流量
void getALLEdgeFlowDecrease(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	//每一条边断掉计算流量
	key_edge_set->EdgeNum = g.nE;
	for (int i=1; i<=g.nE; i++)
	{
		key_edge_set->EdgeInfo[i].Edge = i;
		key_edge_set->EdgeInfo[i].Edge_Class = 'T';
		key_edge_set->EdgeInfo[i].ChangeAmount_c = getMAXflowWithBreakEdge(g, i ,source,sink);
	}

	//按照流量进行排序一下，冒泡，小靠前
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
	return;
}


///////////////////////////////////////////////////////////////////
//2.

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
			if (key_edge_set->EdgeInfo[j].ChangeAmount_p2 < key_edge_set->EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge(key_edge_set, i, j);
			}
		}
	}
	return;
}

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

void reCompute(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, int BreakEdge)
{
	int edgeNum = key_edge_set->EdgeInfo[BreakEdge].Edge;

	int new_maxflow = 0; /*新  最大流*/
	double new_p2 =0;/*新  容量可靠性*/
	Flow new_maxPmaxF;/*新的最大流分布*/
	Edge TempEdge;/*临时存储一条边的信息*/

	//初始化临时存储边
	init_TempEdge(TempEdge);
	//去除某一条A类边
	remove_Edge(g, edgeNum, TempEdge);
	//重新计算最大流new_maxflow，流分布可靠性new_p1，网络可靠性new_p2
	new_p2 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, edgeNum);

	//将断掉边之后的最大流，分布可靠性，容量可靠性分别保存
	key_edge_set->EdgeInfo[BreakEdge].ChangeAmount_c = new_maxflow;
	key_edge_set->EdgeInfo[BreakEdge].ChangeAmount_p2 = new_p2;
	key_edge_set->EdgeInfo[BreakEdge].Edge = edgeNum;
	key_edge_set->EdgeInfo[BreakEdge].Edge_Class = 'T';

	//恢复某一条边
	restoreA_Edge(g,edgeNum, TempEdge);

	return;
}

//根据最大流计算判断是否计算    分布可靠性    和     容量可靠性
void computeP1P2(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, ostream& out_test)
{
	if (key_edge_set->EdgeNum == 0)
	{
		return;
	}

	//根据流量遍历边
	int i=1, j=i+1;
	while(j <= key_edge_set->EdgeNum)
	{
		//流量不同，不需要计算i
		if (key_edge_set->EdgeInfo[j].ChangeAmount_c != key_edge_set->EdgeInfo[i].ChangeAmount_c)
		{
			//前后不同
			i++;
			j++;
			//分布可靠性不用计算，直接置为0
			key_edge_set->EdgeInfo[i].ChangeAmount_p2 = 0;
		}else
		{
			//i和j表示的流量一样，把i表示的先计算出来
			//重新计算   分布可靠性
			reCompute(key_edge_set, g, source, sink, i);

			//找到相同变化的区间
			for (int k =j; k <= key_edge_set->EdgeNum; k++)
			{
				if (key_edge_set->EdgeInfo[k].ChangeAmount_c != key_edge_set->EdgeInfo[i].ChangeAmount_c )
				{
					j = k;
					break;
				}
				
				//i和j一样，计算j的变化p2
				//重新计算  分布可靠性
				reCompute(key_edge_set, g, source, sink, k);

				if ( k == key_edge_set->EdgeNum)
				{
					j = k+1;
					break;
				}				
			}
			//使用冒泡排序一段ChangeAmount_p1
			Bubbling_KeyEdge_1(key_edge_set, i, j-1);
			i=j;
			j=i+1;
		}
	}
	return;
}

///////////////////////////////////////////////////////////////
void computeBase(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, ostream &out_test)
{
	//记录时间
	__int64 start = 0; /*用于测量时间(精确到1ms)*/ 
	__int64 frequency = 0; /*与机器平台相关*/ 
	__int64 counter = 0;
	double timeCost = 0.0;

	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency); 
	QueryPerformanceCounter((LARGE_INTEGER*)&start); /*记录开始时间*/
	//核心计算过程
	//1.分别计算最大流
	getALLEdgeFlowDecrease(key_edge_set, g, source, sink);
	//2.根据最大流计算判断是否计算 分布可靠性
	computeP1P2(key_edge_set, g, source, sink, out_test);

	QueryPerformanceCounter((LARGE_INTEGER*)&counter); /*记录结束时间*/ 
	timeCost = (counter - start) / double(frequency)*1000;/*返回单位是毫秒*/
	out_test<<"V"<<g.nV<<"E"<<g.nE<<",s,t,";
	out_test<<source<<","<<sink<<","<<timeCost<<",";
	return;
}

//////////////////////////////////////////////////////////////////
//以下为接口需要的函数

/*将关键边输出*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"关键边输出如下："<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.EdgeNum;i++)
	{
		out<<setw(5)<<key_edge_set.EdgeInfo[i].Edge<<setw(10)<<key_edge_set.EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_p2<<endl;
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
		key_edge_set.EdgeInfo[i].ChangeAmount_p2 = 0;
	}
	
	return;
}