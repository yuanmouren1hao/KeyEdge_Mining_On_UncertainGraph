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

/*�Ƴ�A���е�ĳһ����*/
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

/*����A��ߵ������ߵ�λ��*/
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



/*����A��߽�������*/
void sortKeyEdge_ABC(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*�������*/
	cout<<"sortKeyEdge_ABC"<<endl;

	int new_maxflow = 0; /*�����������*/
	double new_p1 =0;/*�µ����ɿ���(��ɿ�������ֲ�)*/
	double new_p2 =0;/*�µ�����ɿ���*/
	Flow new_maxPmaxF;/*�µ�������ֲ�*/
	Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/


	/*�����A��ߵĻ�*/
	if (key_edge_set->A_num > 0)
	{
		/*����A���еıߣ�ȥ��һ��*/
		for (int i=1; i<=key_edge_set->A_num; i++)
		{
			/*��ʼ����ʱ�洢��*/
			init_TempEdge(TempEdge);
			new_maxflow = 0;
			new_p1 = 0;
			new_p2 = 0;

			/*ȥ��ĳһ��A���*/
			removeA_Edge(g, key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
			//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
			new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->A_EdgeInfo[i].Edge);

			/* ���ı���������A����Ϣ,������Ƕϵ���֮���ܹ��ﵽ������� */
			key_edge_set->A_EdgeInfo[i].ChangeAmount_c = new_maxflow;
			key_edge_set->A_EdgeInfo[i].ChangeAmount_p1 = new_p1;
			key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = new_p2;

			/* �ָ�ĳһ��A��� */
			restoreA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
		}
	}

	/*�����B��ߵĻ�*/
	if (key_edge_set->B_num > 0)
	{
		/*����A���еıߣ�ȥ��һ��*/
		for (int i=1; i<=key_edge_set->B_num; i++)
		{
			/*��ʼ����ʱ�洢��*/
			init_TempEdge(TempEdge);
			new_maxflow = 0;
			new_p1 = 0;
			new_p2 = 0;

			/*ȥ��ĳһ��A���*/
			removeA_Edge(g, key_edge_set->B_EdgeInfo[i].Edge, TempEdge);
			//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
			new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->B_EdgeInfo[i].Edge);

			/* ���ı���������A����Ϣ,������Ƕϵ���֮���ܹ��ﵽ������� */
			key_edge_set->B_EdgeInfo[i].ChangeAmount_c = new_maxflow;
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p1 = new_p1;
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = new_p2;

			/* �ָ�ĳһ��A��� */
			restoreA_Edge(g,key_edge_set->B_EdgeInfo[i].Edge, TempEdge);
		}
	}

	/*�����C��ߵĻ�*/
	if (key_edge_set->C_num > 0)
	{
		/*����A���еıߣ�ȥ��һ��*/
		for (int i=1; i<=key_edge_set->C_num; i++)
		{
			/*��ʼ����ʱ�洢��*/
			init_TempEdge(TempEdge);
			new_maxflow = 0;
			new_p1 = 0;
			new_p2 = 0;

			/*ȥ��ĳһ��A���*/
			removeA_Edge(g, key_edge_set->C_EdgeInfo[i].Edge, TempEdge);
			//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
			new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->C_EdgeInfo[i].Edge);

			/* ���ı���������A����Ϣ,������Ƕϵ���֮���ܹ��ﵽ������� */
			key_edge_set->C_EdgeInfo[i].ChangeAmount_c = new_maxflow;
			key_edge_set->C_EdgeInfo[i].ChangeAmount_p1 = new_p1;
			key_edge_set->C_EdgeInfo[i].ChangeAmount_p2 = new_p2;

			/* �ָ�ĳһ��A��� */
			restoreA_Edge(g,key_edge_set->C_EdgeInfo[i].Edge, TempEdge);
		}
	}
	
	return;
}

/*��ʼ��*/
void init_TempP(double TempP[])
{
	for (int i=0;i<MAX_E_NUM;i++)
	{
		TempP[i]=0.0;
	}
	return;
}



/*�Ա߽��еĴ����*/
void sortKeyEdge(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set)
{
	//cout<<"StateMtrix->Edge_Num:"<<StateMtrix->Edge_Num<<endl;
	//cout<<"StateMtrix->State_Num:"<<StateMtrix->State_Num<<endl;
	/*��ʼ���ؼ��߼���*/
	key_edge_set->EdgeNum=StateMtrix->Edge_Num;
	int temp_sum=0;
	
	for (int ii=1;ii<=StateMtrix->Edge_Num;ii++)
	{
		temp_sum=0;
		/*ͨ����͵ķ�ʽ����ߵ����*/
		for (int jj=1;jj<=StateMtrix->State_Num;jj++)
		{
			temp_sum+=StateMtrix->State[jj][ii];
		}
		/*���ؼ��ߴ洢��ii��ʾ���Ǳߵı��*/
		if(temp_sum==StateMtrix->State_Num)
		{
			//cout<<ii<<" ��A��ߡ�"<<endl;
			key_edge_set->A_num++;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge=ii;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge_Class='A';
		}
		else if (0==temp_sum)
		{
			//cout<<ii<<" ��C��ߡ�"<<endl;
			key_edge_set->C_num++;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge=ii;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge_Class='C';
		}
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num)
		{
			//cout<<ii<<" ��B���"<<endl;
			key_edge_set->B_num++;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge=ii;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge_Class='B';
		}
	}
}

/*ͨ��״̬�������ߵ����*/
void computeEdgeClass(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*���Է���*/
	sortKeyEdge(StateMtrix, key_edge_set);
	
	/*�Թؼ��߽��ж�������*/
	sortKeyEdge_ABC(key_edge_set, g, source, sink);
	/*�������*/
	
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


/*���ؼ������*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"�ؼ���������£�"<<endl<<"---------------------------------"<<endl;
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
	/*��ʼ���ؼ��߼���*/
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
	/*��ʼ���½���ͼ*/
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