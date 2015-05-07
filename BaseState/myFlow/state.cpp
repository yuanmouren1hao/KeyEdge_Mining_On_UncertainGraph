#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"


using namespace std;

/*移除A类中的某一条边*/
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

/*恢复某一条A类边*/
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


/*初始化临时存储边*/
void init_TempEdge(Edge &TempEdge)
{
	TempEdge.dP=0;
	TempEdge.iC=0;
	TempEdge.iLabel=0;
	return;
}

/*调换A类边的两个边的位置*/
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


/*对于A类边进行排序*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set,Graph& g,int source,int sink,int AllEdge[][5])
{
	//cout<<"sortKeyEdge_A"<<endl;
	/*如果没有A类边的个数为0*/
	if (key_edge_set->A_num==0)
	{
		return;
	}

	//GF gf; /*剩余图*/

	/*临时存储一条边的信息*/
	Edge TempEdge;
	

	/*对于A类中的边，去掉一个*/
	for (int i=1;i<=key_edge_set->A_num;i++)
	{
		/*临时存储 剩余随机流网络可靠性*/
		double remain_flow=1;
		/*最大流*/
		int Fmax;
		Flow Fd;

		/*初始化临时存储边*/
		init_TempEdge(TempEdge);
		/*去除某一条A类边*/
		removeA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge,AllEdge,TempEdge);
		/*计算A边断掉之后剩余最大流的随机流网络的可靠性*/
		remain_flow = old_GetMPMF(g,source, sink,Fmax,Fd);
		/*将改变量保存在A边信息*/
		key_edge_set->A_EdgeInfo[i].ChangeAmount=Fmax;
		key_edge_set->A_EdgeInfo[i].change2 = remain_flow;
		/*恢复某一条A类边*/
		restoreA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge,AllEdge,TempEdge);
	}
	/*使用冒泡排序方式将A类边内部排序*/
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

	/*按照可靠性对A边进行排序*/
	/*TODO 如果需要可以进行排序*/
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

/*调换B类边的两个边的位置*/
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

/*判断边是否在某一个子图中*/
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

/*对于B类边进行排序*/
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
		/*求一个下界子图的可靠性*/
		for (int j=1;j<=StateMtrix->Edge_Num;j++)
		{
			if (1==StateMtrix->State[i][j])
			{
				p*=AllEdge_p[j];
			}
		}
		StateMtrix->SubGraphP[i]=p;
		/*更新下界子图中边对应的最大可靠*/
		for (int jj=1;jj<=StateMtrix->Edge_Num;jj++)
		{
			if (1==StateMtrix->State[i][jj]&&TempP[jj]<p)
			{
				TempP[jj]=p;
			}
		}
	}
	/*更新边的可靠性*/
	for (int i=1;i<=key_edge_set->B_num;i++)
	{
		key_edge_set->B_EdgeInfo[i].ChangeAmount=TempP[key_edge_set->B_EdgeInfo[i].Edge];
	}
	/*对于B类边内部进行排序，冒泡法*/
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
	/*将初始B边中的非B边更新为C类边*/
	for (int i=1;i<=key_edge_set->B_num;i++)
	{
		if (a != key_edge_set->B_EdgeInfo[i].ChangeAmount)
		{
			j=i;
			break;
		}	
	}

	/*因为B累边的数量是变化的，所以先记录下来*/
	int B_num = key_edge_set->B_num;
	for (int bb=j;bb<=B_num;bb++)
	{
		Exchange_B_C(key_edge_set,bb);
	}

	/*将B边断掉后能达到的最大可靠性计算出来*/
	/*不含此边的最大可靠性的状态*/
	for (int i=1;i<=key_edge_set->B_num;i++)
	{
		/*对于B类边中每一条边*/
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

	/*对更新后的B边进行排序，主要是change2*/
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

/*交换C的边*/
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

/*将C类百年进行排序*/
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


/*更新A类边  剩余的流量->变化的流量*/
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

/*通过状态矩阵计算边的类别*/
void computeEdgeClass(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set,Graph& g,int source,int sink,int AllEdge[][5],double AllEdge_p[], int maxFlow)
{
	/*对边进行的大分类*/
	sortKeyEdge(StateMtrix,key_edge_set);
	/*对关键边进行排序*/
	sortKeyEdge_A(key_edge_set,g,source,sink,AllEdge);
	resetEdge_A(key_edge_set,maxFlow);
	sortKeyEdge_B(key_edge_set,AllEdge_p,StateMtrix);
	sortKeyEdge_C(key_edge_set);

	return;
}


/*将关键边输出*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"关键边输出如下："<<endl<<"---------------------------------"<<endl;
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

		key_edge_set.A_EdgeInfo[i].ChangeAmount=0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount=0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount=0;

		key_edge_set.A_EdgeInfo[i].change2=0;
		key_edge_set.B_EdgeInfo[i].change2=0;
		key_edge_set.C_EdgeInfo[i].change2=0;
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