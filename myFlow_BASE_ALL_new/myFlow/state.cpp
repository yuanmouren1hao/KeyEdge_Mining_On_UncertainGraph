#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"

/*�����ⲿ��ԭʼ�ļ�����ɿ�������ֲ����㷨*/
//#include "../Base/partition.cpp"

using namespace std;


//////////////////////////////////////////////////////////////
//new
//////////////////////////////////////////////////////////////

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


//������������
void rankEdge(KeyEdgeSet *key_edge_set)
{
	/*���Ȱ��մ�С���������������ChangeAmount_c�� ʹ��ð�ݷ�����*/
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

	/*�������ChangeAmount_cһ�µ�����£�����ֲ��ɿ���ChangeAmount_p1*/
	int ii =1;
	int jj=ii+1;
	while(jj <= key_edge_set->EdgeNum)
	{
		if (key_edge_set->EdgeInfo[jj].ChangeAmount_c != key_edge_set->EdgeInfo[ii].ChangeAmount_c)
		{
			//ǰ��ͬ
			ii++;
			jj++;
		}else
		{
			//�ҵ���ͬ�仯������
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
			//ʹ��ð������һ��ChangeAmount_p1
			Bubbling_KeyEdge_1(key_edge_set, ii, jj-1);
			ii=jj;
			jj=ii+1;
		}
	}

	/*�ڷֲ��ɿ���ChangeAmount_p1һ�µ�����£����������ɿ���ChangeAmount_p2*/
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

////////////////////////////////////

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

//base_all�㷨��ȡ���е���������ֲ��ɿ��ԣ������ɿ���
void computeBaseAll(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	//�������е�ÿһ����
	for (int i=1; i<=g.nE; i++)
	{
		int new_maxflow = 0; /*��  �����*/
		double new_p1 =0;/*��  �ֲ��ɿ���*/
		double new_p2 =0;/*��  �����ɿ���*/
		Flow new_maxPmaxF;/*�µ�������ֲ�*/
		Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/

		int EdgeNum = i;//�Ƴ��ıߺ�
		/*��ʼ����ʱ�洢��*/
		init_TempEdge(TempEdge);
		/*ȥ��ĳһ��A���*/
		remove_Edge(g, EdgeNum, TempEdge);
		//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
		new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, EdgeNum);

		/* ���ϵ���֮�����������ֲ��ɿ��ԣ������ɿ��Էֱ𱣴� */
		key_edge_set->EdgeNum++;
		key_edge_set->EdgeInfo[EdgeNum].ChangeAmount_c = new_maxflow;
		key_edge_set->EdgeInfo[EdgeNum].ChangeAmount_p1 = new_p1;
		key_edge_set->EdgeInfo[EdgeNum].ChangeAmount_p2 = new_p2;
		key_edge_set->EdgeInfo[EdgeNum].Edge = EdgeNum;
		key_edge_set->EdgeInfo[EdgeNum].Edge_Class = 'S';

		/* �ָ�ĳһ��A��� */
		restoreA_Edge(g,EdgeNum, TempEdge);
	}

	//�ȶ���������֮����ͬ���������շֲ��ɿ���������ͬ�ķֲ��ɿ��԰��������ɿ�������
	rankEdge(key_edge_set);

	return;
}

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