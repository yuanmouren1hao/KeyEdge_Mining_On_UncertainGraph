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

/*����2���ߵ�λ��*/
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

//����һ���߶ϵ�֮��
int getMAXflowWithBreakEdge(Graph& g, int BreakEdge, int source, int sink)
{
	//ʹ��ԭ����ͼ  �ȹ��죬��ȡ�������Ȼ��ָ�ԭ��ͼ������
	g.nE--;

	int aa,bb,cc,dd;//��BreakEdge����ʼ�㣬�յ㣬��������Чλ
	double ff;//��BreakEdge�Ŀɿ���
	aa = g.AllEdge[BreakEdge][1];//��ʼ��
	bb = g.AllEdge[BreakEdge][2];//�յ�
	cc = g.AllEdge[BreakEdge][3];//����
	dd = g.AllEdge[BreakEdge][4];//��Чλ���˴���ʹ�ã�
	ff = g.AllEdge_p[BreakEdge];//�ɿ���

	g.matrix[aa][bb].iLabel = 0;
	g.matrix[aa][bb].iC = 0;
	g.matrix[aa][bb].dP = 0;
	//ʹ��ԭ��Dinic�㷨���������
	Flow f;
	GF gf;
	int max_flow_remine = Dinic(g, source, sink, f, gf);

	//�������֮��  ��Ҫ�ָ�ԭ�е�ͼ
	g.nE++;
	g.matrix[aa][bb].iLabel = BreakEdge;
	g.matrix[aa][bb].iC = cc;
	g.matrix[aa][bb].dP = ff;

	//��������
	return max_flow_remine;
}

//����ÿ���߶ϵ�֮�������
void getALLEdgeFlowDecrease(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	//ÿһ���߶ϵ���������
	key_edge_set->EdgeNum = g.nE;
	for (int i=1; i<=g.nE; i++)
	{
		key_edge_set->EdgeInfo[i].Edge = i;
		key_edge_set->EdgeInfo[i].Edge_Class = 'T';
		key_edge_set->EdgeInfo[i].ChangeAmount_c = getMAXflowWithBreakEdge(g, i ,source,sink);
	}

	//����������������һ�£�ð�ݣ�С��ǰ
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

//����һ��p1
void Bubbling_KeyEdge_1(KeyEdgeSet *key_edge_set, int begin_, int end_)
{
	/*��Ҫ����i<j, ����i��j����������*/
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

//����һ��p2
void Bubbling_KeyEdge_2(KeyEdgeSet *key_edge_set, int begin_, int end_)
{
	/*��Ҫ����i<j, ����i��j����������*/
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

/*�ָ�ĳһ��A���*/
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

/*��ʼ����ʱ�洢��*/
void init_TempEdge(Edge &TempEdge)
{
	TempEdge.dP=0;
	TempEdge.iC=0;
	TempEdge.iLabel=0;
	return;
}

/*�Ƴ�A���е�ĳһ����*/
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

	int new_maxflow = 0; /*��  �����*/
	double new_p1 =0;/*��  �ֲ��ɿ���*/
	double new_p2 =0;/*��  �����ɿ���*/
	Flow new_maxPmaxF;/*�µ�������ֲ�*/
	Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/

	//��ʼ����ʱ�洢��
	init_TempEdge(TempEdge);
	//ȥ��ĳһ��A���
	remove_Edge(g, edgeNum, TempEdge);
	//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
	new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, edgeNum);

	//���ϵ���֮�����������ֲ��ɿ��ԣ������ɿ��Էֱ𱣴�
	key_edge_set->EdgeInfo[BreakEdge].ChangeAmount_c = new_maxflow;
	key_edge_set->EdgeInfo[BreakEdge].ChangeAmount_p1 = new_p1;
	key_edge_set->EdgeInfo[BreakEdge].ChangeAmount_p2 = new_p2;
	key_edge_set->EdgeInfo[BreakEdge].Edge = edgeNum;
	key_edge_set->EdgeInfo[BreakEdge].Edge_Class = 'T';

	//�ָ�ĳһ����
	restoreA_Edge(g,edgeNum, TempEdge);

	return;
}

//��������������ж��Ƿ����    �ֲ��ɿ���    ��     �����ɿ���
void computeP1P2(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, ostream& out_test)
{
	if (key_edge_set->EdgeNum == 0)
	{
		return;
	}

	//��¼��Ҫ�ظ�����ĸ���
	int needRecomputerNum = 0;
	//ͳ��ʱ��
	__int64 start = 0; /*���ڲ���ʱ��(��ȷ��1ms)*/ 
	__int64 frequency = 0; /*�����ƽ̨���*/ 
	__int64 counter = 0;
	double timeCost = 0.0;

	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency); 
	QueryPerformanceCounter((LARGE_INTEGER*)&start); /*��¼��ʼʱ��*/

	//��������������
	int i=1, j=i+1;
	while(j <= key_edge_set->EdgeNum)
	{
		//������ͬ������Ҫ����i
		if (key_edge_set->EdgeInfo[j].ChangeAmount_c != key_edge_set->EdgeInfo[i].ChangeAmount_c)
		{
			//ǰ��ͬ
			i++;
			j++;
			//�ֲ��ɿ��Բ��ü��㣬ֱ����Ϊ0
			key_edge_set->EdgeInfo[i].ChangeAmount_p1 = 0;
			key_edge_set->EdgeInfo[i].ChangeAmount_p2 = 0;
		}else
		{
			//i��j��ʾ������һ������i��ʾ���ȼ������
			reCompute(key_edge_set, g, source, sink, i);
			needRecomputerNum++;

			//�ҵ���ͬ�仯������
			for (int k =j; k <= key_edge_set->EdgeNum; k++)
			{
				if (key_edge_set->EdgeInfo[k].ChangeAmount_c != key_edge_set->EdgeInfo[i].ChangeAmount_c )
				{
					j = k;
					break;
				}
				
				//i��jһ��������j�ı仯p1��p2
				reCompute(key_edge_set, g, source, sink, k);
				needRecomputerNum++;

				if ( k == key_edge_set->EdgeNum)
				{
					j = k+1;
					break;
				}				
			}
			//ʹ��ð������һ��ChangeAmount_p1
			Bubbling_KeyEdge_1(key_edge_set, i, j-1);
			i=j;
			j=i+1;
		}
	}
	//�����ʡ�ļ������
	out_test<<key_edge_set->EdgeNum-needRecomputerNum<<",";
	QueryPerformanceCounter((LARGE_INTEGER*)&counter); /*��¼����ʱ��*/ 
	timeCost = (counter - start) / double(frequency)*1000;/*���ص�λ�Ǻ���*/
	out_test<<timeCost/needRecomputerNum<<endl;


	//�ڷֲ��ɿ���ChangeAmount_p1һ�µ�����£����������ɿ���ChangeAmount_p2
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
			//�ҵ���ͬ�仯������
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
			//ʹ��ð������һ��
			Bubbling_KeyEdge_2(key_edge_set, iii, jjj-1);
			iii=jjj;
			jjj=iii+1;
		}
	}

	return;
}

///////////////////////////////////////////////////////////////
void computeBase(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, ostream &out_test)
{
	//��¼ʱ��
	__int64 start = 0; /*���ڲ���ʱ��(��ȷ��1ms)*/ 
	__int64 frequency = 0; /*�����ƽ̨���*/ 
	__int64 counter = 0;
	double timeCost = 0.0;

	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency); 
	QueryPerformanceCounter((LARGE_INTEGER*)&start); /*��¼��ʼʱ��*/

	//1.�ֱ���������
	getALLEdgeFlowDecrease(key_edge_set, g, source, sink);

	QueryPerformanceCounter((LARGE_INTEGER*)&counter); /*��¼����ʱ��*/ 
	timeCost = (counter - start) / double(frequency)*1000;/*���ص�λ�Ǻ���*/
	out_test<<"V"<<g.nV<<"E"<<g.nE<<",s,t,";
	out_test<<source<<","<<sink<<","<<timeCost<<",";


	//2.��������������ж��Ƿ����    �ֲ��ɿ���    ��     �����ɿ���
	computeP1P2(key_edge_set, g, source, sink, out_test);


	//3.���򣬸����������ֲ��ɿ��ԣ������ɿ���  ����ά��
	//��������һ�������Ѿ����

	return;
}

//////////////////////////////////////////////////////////////////
//����Ϊ�ӿ���Ҫ�ĺ���

/*���ؼ������*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"�ؼ���������£�"<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.EdgeNum;i++)
	{
		out<<setw(5)<<key_edge_set.EdgeInfo[i].Edge<<setw(10)<<key_edge_set.EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	out<<"---------------------------------"<<endl<<"---------------------------------"<<endl<<endl;
}


void init_KeyEdgeSet(KeyEdgeSet &key_edge_set)
{
	/*��ʼ���ؼ��߼���*/
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