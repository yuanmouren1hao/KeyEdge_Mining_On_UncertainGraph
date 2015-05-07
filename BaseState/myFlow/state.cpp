#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"


using namespace std;

/*�Ƴ�A���е�ĳһ����*/
void removeA_Edge(Graph& g,int i,int AllEdge[][5],Edge &TempEdge)
{
	int u=AllEdge[i][1];
	int v=AllEdge[i][2];
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
void restoreA_Edge(Graph& g,int i,int AllEdge[][5],Edge &TempEdge)
{
	int u=AllEdge[i][1];
	int v=AllEdge[i][2];
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
void Exchange_A(KeyEdgeSet *key_edge_set,int i,int j)
{
	//cout<<i<<"--->"<<j<<endl;
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount=key_edge_set->A_EdgeInfo[i].ChangeAmount;
	TempKeyEdge.Edge=key_edge_set->A_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class=key_edge_set->A_EdgeInfo[i].Edge_Class;

	key_edge_set->A_EdgeInfo[i].ChangeAmount=key_edge_set->A_EdgeInfo[j].ChangeAmount;
	key_edge_set->A_EdgeInfo[i].Edge=key_edge_set->A_EdgeInfo[j].Edge;
	key_edge_set->A_EdgeInfo[i].Edge_Class=key_edge_set->A_EdgeInfo[j].Edge_Class;

	key_edge_set->A_EdgeInfo[j].ChangeAmount=TempKeyEdge.ChangeAmount;
	key_edge_set->A_EdgeInfo[j].Edge=TempKeyEdge.Edge;
	key_edge_set->A_EdgeInfo[j].Edge_Class=TempKeyEdge.Edge_Class;

	return;
}

double getRemainA_Flow(Flow Fd,GF gf )
{
	for (int i=1;i<=7;i++)
	{
		for (int j = 1;j <= 7; j++ )
		{
			cout<<Fd[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;
	for (int i=1;i<=7;i++)
	{
		for (int j = 1;j <= 7; j++ )
		{
			cout<<gf[i][j]<<"  ";
		}
		cout<<endl;
	}
	return 0;
}


/*����A��߽�������*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set,Graph& g,int source,int sink,int AllEdge[][5])
{
	//cout<<"sortKeyEdge_A"<<endl;
	/*���û��A��ߵĸ���Ϊ0*/
	if (key_edge_set->A_num==0)
	{
		return;
	}

	//GF gf; /*ʣ��ͼ*/

	/*��ʱ�洢һ���ߵ���Ϣ*/
	Edge TempEdge;
	

	/*����A���еıߣ�ȥ��һ��*/
	for (int i=1;i<=key_edge_set->A_num;i++)
	{
		/*��ʱ�洢 ʣ�����������ɿ���*/
		double remain_flow=1;
		/*�����*/
		int Fmax;
		Flow Fd;

		/*��ʼ����ʱ�洢��*/
		init_TempEdge(TempEdge);
		/*ȥ��ĳһ��A���*/
		removeA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge,AllEdge,TempEdge);
		/*����A�߶ϵ�֮��ʣ������������������Ŀɿ���*/
		remain_flow = old_GetMPMF(g,source, sink,Fmax,Fd);
		/*���ı���������A����Ϣ*/
		key_edge_set->A_EdgeInfo[i].ChangeAmount=Fmax;
		key_edge_set->A_EdgeInfo[i].change2 = remain_flow;
		/*�ָ�ĳһ��A���*/
		restoreA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge,AllEdge,TempEdge);
	}
	/*ʹ��ð������ʽ��A����ڲ�����*/
	for (int mi=1;mi<key_edge_set->A_num;mi++)
	{
		for (int mj=mi+1;mj<=key_edge_set->A_num;mj++)
		{
			if (key_edge_set->A_EdgeInfo[mi].ChangeAmount>key_edge_set->A_EdgeInfo[mj].ChangeAmount)
			{
				Exchange_A(key_edge_set,mi,mj);
			}
		}
	}

	/*���տɿ��Զ�A�߽�������*/
	/*TODO �����Ҫ���Խ�������*/
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

/*����B��ߵ������ߵ�λ��*/
void Exchange_B(KeyEdgeSet *key_edge_set,int i,int j)
{
	//cout<<i<<"--->"<<j<<endl;
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount=key_edge_set->B_EdgeInfo[i].ChangeAmount;
	TempKeyEdge.Edge=key_edge_set->B_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class=key_edge_set->B_EdgeInfo[i].Edge_Class;
	TempKeyEdge.change2=key_edge_set->B_EdgeInfo[i].change2;

	key_edge_set->B_EdgeInfo[i].ChangeAmount=key_edge_set->B_EdgeInfo[j].ChangeAmount;
	key_edge_set->B_EdgeInfo[i].Edge=key_edge_set->B_EdgeInfo[j].Edge;
	key_edge_set->B_EdgeInfo[i].Edge_Class=key_edge_set->B_EdgeInfo[j].Edge_Class;
	key_edge_set->B_EdgeInfo[i].change2=key_edge_set->B_EdgeInfo[j].change2;

	key_edge_set->B_EdgeInfo[j].ChangeAmount=TempKeyEdge.ChangeAmount;
	key_edge_set->B_EdgeInfo[j].Edge=TempKeyEdge.Edge;
	key_edge_set->B_EdgeInfo[j].Edge_Class=TempKeyEdge.Edge_Class;
	key_edge_set->B_EdgeInfo[j].change2=TempKeyEdge.change2;

	return;
}

void Exchange_B_C(KeyEdgeSet *key_edge_set,int bb)
{
	key_edge_set->C_num++;
	key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge=key_edge_set->B_EdgeInfo[bb].Edge;
	key_edge_set->C_EdgeInfo[key_edge_set->C_num].ChangeAmount=key_edge_set->B_EdgeInfo[bb].ChangeAmount;
	key_edge_set->C_EdgeInfo[key_edge_set->C_num].change2=key_edge_set->B_EdgeInfo[bb].change2;
	key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge_Class='C';

	key_edge_set->B_EdgeInfo[bb].Edge_Class=0;
	key_edge_set->B_EdgeInfo[bb].ChangeAmount=0;
	key_edge_set->B_EdgeInfo[bb].change2=0;
	key_edge_set->B_EdgeInfo[bb].Edge=0;
	key_edge_set->B_num--;
	return;
}

/*�жϱ��Ƿ���ĳһ����ͼ��*/
int is_in_subgraph(int subgraph, int edge, Lower_subGraph *StateMtrix)
{
		if (StateMtrix->State[subgraph][edge]==1)
		{
			return true;
		}
		else
		{
			return false;
		}
}

/*����B��߽�������*/
 void sortKeyEdge_B(KeyEdgeSet *key_edge_set,double AllEdge_p[],Lower_subGraph *StateMtrix)
{
	if (0 == key_edge_set->B_num)
	{
		/*no B edge*/
		return;
	}
	//cout<<"sortKeyEdge_B"<<endl;
	double TempP[MAX_E_NUM];
	init_TempP(TempP);
	double p;
	for (int i=1;i<=StateMtrix->State_Num;i++)
	{
		p=1.0;
		/*��һ���½���ͼ�Ŀɿ���*/
		for (int j=1;j<=StateMtrix->Edge_Num;j++)
		{
			if (1==StateMtrix->State[i][j])
			{
				p*=AllEdge_p[j];
			}
		}
		StateMtrix->SubGraphP[i]=p;
		/*�����½���ͼ�б߶�Ӧ�����ɿ�*/
		for (int jj=1;jj<=StateMtrix->Edge_Num;jj++)
		{
			if (1==StateMtrix->State[i][jj]&&TempP[jj]<p)
			{
				TempP[jj]=p;
			}
		}
	}
	/*���±ߵĿɿ���*/
	for (int i=1;i<=key_edge_set->B_num;i++)
	{
		key_edge_set->B_EdgeInfo[i].ChangeAmount=TempP[key_edge_set->B_EdgeInfo[i].Edge];
	}
	/*����B����ڲ���������ð�ݷ�*/
	for (int mi=1;mi<key_edge_set->B_num;mi++)
	{
		for (int mj=mi+1;mj<=key_edge_set->B_num;mj++)
		{
			if (key_edge_set->B_EdgeInfo[mi].ChangeAmount < key_edge_set->B_EdgeInfo[mj].ChangeAmount)
			{
				Exchange_B(key_edge_set,mi,mj);
			}
		}
	}

	double a=key_edge_set->B_EdgeInfo[1].ChangeAmount;
	int j=1;
	/*����ʼB���еķ�B�߸���ΪC���*/
	for (int i=1;i<=key_edge_set->B_num;i++)
	{
		if (a != key_edge_set->B_EdgeInfo[i].ChangeAmount)
		{
			j=i;
			break;
		}	
	}

	/*��ΪB�۱ߵ������Ǳ仯�ģ������ȼ�¼����*/
	int B_num = key_edge_set->B_num;
	for (int bb=j;bb<=B_num;bb++)
	{
		Exchange_B_C(key_edge_set,bb);
	}

	/*��B�߶ϵ����ܴﵽ�����ɿ��Լ������*/
	/*�����˱ߵ����ɿ��Ե�״̬*/
	for (int i=1;i<=key_edge_set->B_num;i++)
	{
		/*����B�����ÿһ����*/
		for (int j=1; j<=StateMtrix->State_Num; j++)
		{
			if (!is_in_subgraph(j,key_edge_set->B_EdgeInfo[i].Edge,StateMtrix) && StateMtrix->SubGraphP[j]<key_edge_set->B_EdgeInfo[1].ChangeAmount)
			{
				if (key_edge_set->B_EdgeInfo[i].change2<StateMtrix->SubGraphP[j])
				{
					key_edge_set->B_EdgeInfo[i].change2=StateMtrix->SubGraphP[j];
				}
			}
		}
	}

	/*�Ը��º��B�߽���������Ҫ��change2*/
	for (int mi=1;mi<key_edge_set->B_num;mi++)
	{
		for (int mj=mi+1;mj<=key_edge_set->B_num;mj++)
		{
			if (key_edge_set->B_EdgeInfo[mi].change2 > key_edge_set->B_EdgeInfo[mj].change2)
			{
				Exchange_B(key_edge_set,mi,mj);
			}
		}
	}

	return;
}

/*����C�ı�*/
void Exchange_C(KeyEdgeSet *key_edge_set,int i,int j)
{
	//cout<<i<<"--->"<<j<<endl;
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount=key_edge_set->C_EdgeInfo[i].ChangeAmount;
	TempKeyEdge.Edge=key_edge_set->C_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class=key_edge_set->C_EdgeInfo[i].Edge_Class;
	TempKeyEdge.change2=key_edge_set->C_EdgeInfo[i].change2;

	key_edge_set->C_EdgeInfo[i].ChangeAmount=key_edge_set->C_EdgeInfo[j].ChangeAmount;
	key_edge_set->C_EdgeInfo[i].Edge=key_edge_set->C_EdgeInfo[j].Edge;
	key_edge_set->C_EdgeInfo[i].Edge_Class=key_edge_set->C_EdgeInfo[j].Edge_Class;
	key_edge_set->C_EdgeInfo[i].change2=key_edge_set->C_EdgeInfo[j].change2;

	key_edge_set->C_EdgeInfo[j].ChangeAmount=TempKeyEdge.ChangeAmount;
	key_edge_set->C_EdgeInfo[j].Edge=TempKeyEdge.Edge;
	key_edge_set->C_EdgeInfo[j].Edge_Class=TempKeyEdge.Edge_Class;
	key_edge_set->C_EdgeInfo[j].change2=TempKeyEdge.change2;

	return;
}

/*��C������������*/
void sortKeyEdge_C(KeyEdgeSet *key_edge_set)
{
	if (0 == key_edge_set->C_num)
	{
		return;
	}
	for (int mi=1;mi<key_edge_set->C_num;mi++)
	{
		for (int mj=mi+1;mj<=key_edge_set->C_num;mj++)
		{
			if (key_edge_set->C_EdgeInfo[mi].ChangeAmount < key_edge_set->C_EdgeInfo[mj].ChangeAmount)
			{
				Exchange_C(key_edge_set,mi,mj);
			}
		}
	}
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


/*����A���  ʣ�������->�仯������*/
void resetEdge_A(KeyEdgeSet *key_edge_set, int maxFlow)
{
	if (0 == key_edge_set->A_num)
	{
		return;
	}
	for (int i=1;i<=key_edge_set->A_num;i++)
	{
		key_edge_set->A_EdgeInfo[i].ChangeAmount = maxFlow - key_edge_set->A_EdgeInfo[i].ChangeAmount;
		if (key_edge_set->A_EdgeInfo[i].ChangeAmount == maxFlow && key_edge_set->A_EdgeInfo[i].change2 == 0)
		{
			key_edge_set->A_EdgeInfo[i].change2 = 1;
		}
	}
	return;
}

/*ͨ��״̬�������ߵ����*/
void computeEdgeClass(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set,Graph& g,int source,int sink,int AllEdge[][5],double AllEdge_p[], int maxFlow)
{
	/*�Ա߽��еĴ����*/
	sortKeyEdge(StateMtrix,key_edge_set);
	/*�Թؼ��߽�������*/
	sortKeyEdge_A(key_edge_set,g,source,sink,AllEdge);
	resetEdge_A(key_edge_set,maxFlow);
	sortKeyEdge_B(key_edge_set,AllEdge_p,StateMtrix);
	sortKeyEdge_C(key_edge_set);

	return;
}


/*���ؼ������*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"�ؼ���������£�"<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.A_num;i++)
	{
		out<<setw(5)<<key_edge_set.A_EdgeInfo[i].Edge<<setw(5)<<key_edge_set.A_EdgeInfo[i].Edge_Class<<setw(15)<<key_edge_set.A_EdgeInfo[i].ChangeAmount<<setw(5)<<key_edge_set.A_EdgeInfo[i].change2<<endl;
	}
	for (int i=1;i<=key_edge_set.B_num;i++)
	{
		out<<setw(5)<<key_edge_set.B_EdgeInfo[i].Edge<<setw(5)<<key_edge_set.B_EdgeInfo[i].Edge_Class<<setw(15)<<key_edge_set.B_EdgeInfo[i].ChangeAmount<<setw(5)<<key_edge_set.B_EdgeInfo[i].change2<<endl;
	}
	for (int i=1;i<=key_edge_set.C_num;i++)
	{
		out<<setw(5)<<key_edge_set.C_EdgeInfo[i].Edge<<setw(5)<<key_edge_set.C_EdgeInfo[i].Edge_Class<<setw(15)<<key_edge_set.C_EdgeInfo[i].ChangeAmount<<setw(5)<<key_edge_set.C_EdgeInfo[i].change2<<endl;
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

		key_edge_set.A_EdgeInfo[i].ChangeAmount=0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount=0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount=0;

		key_edge_set.A_EdgeInfo[i].change2=0;
		key_edge_set.B_EdgeInfo[i].change2=0;
		key_edge_set.C_EdgeInfo[i].change2=0;
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