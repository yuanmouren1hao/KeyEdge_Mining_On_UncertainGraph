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

/*����Chin-Chia Jane and Yih-Wenn Laih�ķ����õ�C(j)(1-q)*/
void GetCurrentC_j_new(Collection &c,StateSet &Cj)
{
	/*�ֱ��ʼ��Cj�����½�*/
	memcpy(Cj.lower,c.lower,(Cj.numE+1)*sizeof(int)); 
	memcpy(Cj.upper,c.upper,(Cj.numE+1)*sizeof(int));

	/*��ʼ(��һ��)Collection={|E|,(00,...0),(11,...,1),(NULL),0,0}
	��������һ��*/
	if(c.I == 0){ return ; }

	/*��֤Ad_c[1],...,��֤Ad_c[j-1]��Щ�߶�Ҫ���ڻ�����
	�ڶ�Ԫ����£��൱��Ҫ�����½�(���ʾ���ȡ�Ͻ�)*/
	for(int i = 1; i < c.j; i++)
	{
		Cj.lower[c.Ad_C[i]] = Cj.upper[c.Ad_C[i]];
	}
	/*�Ե�j��Ad_[j]��Ӧ��pivotҪС�ڻ�����
	�ڶ�Ԫ����£��൱����ҪС���Ͻ�(���ʾ���ȡ�½�)*/
	Cj.upper[c.Ad_C[c.j]] = Cj.lower[c.Ad_C[c.j]];
}

/*��״̬����ת���ɶ�Ӧͼg_target*/
void FromG2G_new(Graph& g_original,Graph& g_target,int* lu)
{
	g_target.Init();
	g_target.nV = g_original.nV;  
	g_target.nId = g_original.nId;                                                /*ͼ�ı�ſ��Բ���*/
	for(int u = 1; u <= g_original.nV; u++)
	{
		for(int v = 1; v <= g_original.nV; v++)
		{
			if((g_original.matrix[u][v].iC > 0) //��
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


/*ͨ���������Ӧ��������ֲ���⻮�ֵ�x0*/
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

/*ͨ���պ������뻮���ߵõ�pivot����ad(c)*/
/*��������ad(c)���ϵĴ�С��Ad_c*/
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

/*�������û����߽��л���,�õ��ܹ�����
C(1),...,C(I)��collection*/
void GetNextCollection_new(StateSet &Cj,int *x0,Collection &nextCollection)
{
	nextCollection.j = 1;                                                         /*�趨��Ҫ����Ҫ���Cj*/
	nextCollection.I = GetAd_C_new(Cj,x0,nextCollection.Ad_C);                        /*�趨I��Ad_c[]*/
	memcpy(nextCollection.lower,Cj.lower,(nextCollection.numE+1)*sizeof(int));    /*�趨lower[]*/
	memcpy(nextCollection.upper,Cj.upper,(nextCollection.numE+1)*sizeof(int));    /*�趨upper[]*/
}


/*����Chin-Chia Jane and Yih-Wenn Laih�������C0*/
void GetCurrentC0_new(Collection &c,StateSet &C0)
{
	memcpy(C0.lower,c.lower,(C0.numE+1)*sizeof(int));
	memcpy(C0.upper,c.upper,(C0.numE+1)*sizeof(int));
	for(int i = 1; i <= c.I; i++)
	{
		C0.lower[c.Ad_C[i]] = C0.upper[c.Ad_C[i]];
	}
}

/*��ͼl����ͼ�ռ�����ͼ�����ĸ���*/
double CalculateP_new(Graph& g,int* l)
{
	double dp = 1.0;			
	for(int u = 1; u <= g.nV; u++)
	{
		for(int v = 1; v <= g.nV; v++)
		{
			if(g.matrix[u][v].iC > 0) 
			{
				if(l[g.matrix[u][v].iLabel] == 1)                                 /*ֻ��l��Ϊ1�ı߼Ǹ���*/
				{
					dp = dp*(g.matrix[u][v].dP);
				}
			}
		}
	}
	return dp;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//���ϲ�Ҫ

//////////////////////////////////////////////////////////////////////////
//1.��������  �߶ϵ�   ����������

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

//����ͼ ���б߶ϵ�֮��������������������һ����ά����AllEdgeFlowDecrease��
void getALLEdgeFlowDecrease(int AllEdgeFlowDecrease[][2], Graph& g, int source, int sink)
{
	for (int i=1; i<=g.nE; i++)
	{
		AllEdgeFlowDecrease[i][0] = i;
		AllEdgeFlowDecrease[i][1] = getMAXflowWithBreakEdge(g, i, source, sink);
	}
	return;
}

//////////////////////////////////////////////////////////////////////////
//2.���Է���  ABC���

//���Է���
void classifyABC(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set, int AllEdgeFlowDecrease[][2])
{
	//��ʼ���ؼ��߼���
	key_edge_set->EdgeNum=StateMtrix->Edge_Num;
	int temp_sum=0;

	for (int ii=1; ii<=StateMtrix->Edge_Num; ii++)
	{
		temp_sum = 0;
		//ͨ����͵ķ�ʽ����ߵ����
		for (int jj=1; jj <= StateMtrix->State_Num; jj++)
		{
			if (StateMtrix->State[jj][ii] == 0)
			{
				temp_sum += 0;
			}
			else if (StateMtrix->State[jj][ii] == 1)
			{
				temp_sum += 1;
			}
			else if (StateMtrix->State[jj][ii] == 2)
			{
				temp_sum += 0;
			}
		}

		//���ؼ��ߴ洢��ii��ʾ���Ǳߵı��
		if(temp_sum == StateMtrix->State_Num)
		{
			//A���
			key_edge_set->A_num++;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge=ii;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge_Class='A';
			//���� �����
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].ChangeAmount_c = AllEdgeFlowDecrease[ii][1];
		}
		else if(0<temp_sum && temp_sum<StateMtrix->State_Num && StateMtrix->State[0][ii] == 1)
		{
			//B���
			key_edge_set->B_num++;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge=ii;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge_Class='B';
			//���� �����
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].ChangeAmount_c = AllEdgeFlowDecrease[ii][1];
		}
		else if(0<=temp_sum && temp_sum <StateMtrix->State_Num && StateMtrix->State[0][ii]==0)
		{
			//C���
			key_edge_set->C_num++;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge=ii;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge_Class='C';
			//���� �����
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].ChangeAmount_c = AllEdgeFlowDecrease[ii][1];
		}
	}
}


//////////////////////////////////////////////////////////////////////////

//����2���ߵ�λ��  A��
void Exchange_KeyEdge_A(KeyEdgeSet *key_edge_set,int i,int j)
{
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c  = key_edge_set->A_EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.ChangeAmount_p2 = key_edge_set->A_EdgeInfo[i].ChangeAmount_p2;
	TempKeyEdge.Edge			= key_edge_set->A_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class		= key_edge_set->A_EdgeInfo[i].Edge_Class;

	key_edge_set->A_EdgeInfo[i].ChangeAmount_c		=  key_edge_set->A_EdgeInfo[j].ChangeAmount_c	;
	key_edge_set->A_EdgeInfo[i].ChangeAmount_p2		=  key_edge_set->A_EdgeInfo[j].ChangeAmount_p2;
	key_edge_set->A_EdgeInfo[i].Edge				=  key_edge_set->A_EdgeInfo[j].Edge			;
	key_edge_set->A_EdgeInfo[i].Edge_Class			=  key_edge_set->A_EdgeInfo[j].Edge_Class		;

	key_edge_set->A_EdgeInfo[j].ChangeAmount_c		=  TempKeyEdge.ChangeAmount_c  ;
	key_edge_set->A_EdgeInfo[j].ChangeAmount_p2		=  TempKeyEdge.ChangeAmount_p2 ;
	key_edge_set->A_EdgeInfo[j].Edge				=  TempKeyEdge.Edge			   ;
	key_edge_set->A_EdgeInfo[j].Edge_Class			=  TempKeyEdge.Edge_Class	   ;

	return;
}

/*�ָ�ĳһ����*/
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

/*�Ƴ�ĳһ����*/
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

//���¼���һ��A��
void reCompute_A(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, int BreakEdge)
{
	int edgeNum = key_edge_set->A_EdgeInfo[BreakEdge].Edge;

	int new_maxflow = 0; /*��  �����*/
	double new_p2 =0;/*��  �����ɿ���*/
	Flow new_maxPmaxF;/*�µ�������ֲ�*/
	Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/

	//��ʼ����ʱ�洢��
	init_TempEdge(TempEdge);
	//ȥ��ĳһ��A���
	remove_Edge(g, edgeNum, TempEdge);
	//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
	new_p2 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, edgeNum);

	//���ϵ���֮�����������ֲ��ɿ��ԣ������ɿ��Էֱ𱣴�
	key_edge_set->A_EdgeInfo[BreakEdge].ChangeAmount_p2 = new_p2;
	//�ָ�ĳһ����
	restoreA_Edge(g,edgeNum, TempEdge);

	return;
}

//����һ��p1
void Bubbling_KeyEdge_A_1(KeyEdgeSet *key_edge_set, int begin_, int end_)
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
			if (key_edge_set->A_EdgeInfo[j].ChangeAmount_p2 < key_edge_set->A_EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge_A(key_edge_set, i, j);
			}
		}
	}
	return;
}

//����һ��p2  A��
void Bubbling_KeyEdge_A_2(KeyEdgeSet *key_edge_set, int begin_, int end_)
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
			if (key_edge_set->A_EdgeInfo[j].ChangeAmount_p2 < key_edge_set->A_EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge_A(key_edge_set, i, j);
			}
		}
	}
	return;
}


void reCompute_APA(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, int BreakEdge)
{
	int edgeNum = key_edge_set->A_EdgeInfo[BreakEdge].Edge;

	double new_p2 =1;/*��  �ֲ��ɿ���*/
	Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/
	Flow Fd;  
	GF gf;
	int g_num_of_edges = g.nE;

	//��ʼ����ʱ�洢��
	init_TempEdge(TempEdge);
	//ȥ��ĳһ��A���
	remove_Edge(g, edgeNum, TempEdge);
	//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
	int max_flow = Dinic(g, source, sink, Fd, gf);
	//new_p2 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, edgeNum);
	for (int i=1; i<=g_num_of_edges; i++)
	{
		int _begin = g.AllEdge[i][1];
		int _end = g.AllEdge[i][2];
		if (Fd[_begin][_end] > 0)
		{
			new_p2 *= g.AllEdge_p[i];
		}
	}

	//���ϵ���֮�����������ֲ��ɿ��ԣ������ɿ��Էֱ𱣴�
	key_edge_set->A_EdgeInfo[BreakEdge].ChangeAmount_p2 = new_p2;
	//�ָ�ĳһ����
	restoreA_Edge(g,edgeNum, TempEdge);

	return;
}

//A������¼���
void sortKeyEdge_A(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	if (key_edge_set->A_num <= 1)
	{
		return;
	}
	//��A����������������  ð�ݣ�С��ǰ
	for (int ia =1; ia< key_edge_set->A_num; ia++)
	{
		for (int ja = ia+1; ja<= key_edge_set->A_num; ja++)
		{
			if (key_edge_set->A_EdgeInfo[ja].ChangeAmount_c < key_edge_set->A_EdgeInfo[ia].ChangeAmount_c)
			{
				Exchange_KeyEdge_A(key_edge_set, ia,ja);
			}
		}
	}

	//��������������
	int i=1, j=i+1;
	while(j <= key_edge_set->A_num)
	{
		//��������һ��
		if (i == key_edge_set->A_num)
		{
			//�ֲ��ɿ��Բ��ü��㣬ֱ����Ϊ0
			key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = 0;
			break;
		}

		//������ͬ������Ҫ����i
		if (key_edge_set->A_EdgeInfo[j].ChangeAmount_c != key_edge_set->A_EdgeInfo[i].ChangeAmount_c)
		{
			//ǰ��ͬ
			i++;
			j++;
			//�ֲ��ɿ��Բ��ü��㣬ֱ����Ϊ0
			key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = 0;
		}else
		{
			//i��j��ʾ������һ������i��ʾ���ȼ������

			//�����i����������ֲ�
			reCompute_APA(key_edge_set, g, source, sink, i);

			//�ҵ���ͬ�仯������
			for (int k =j; k <= key_edge_set->A_num; k++)
			{
				if (key_edge_set->A_EdgeInfo[k].ChangeAmount_c != key_edge_set->A_EdgeInfo[i].ChangeAmount_c )
				{
					j = k;
					break;
				}

				//i��jһ��������j�ı仯p1��p2
				reCompute_APA(key_edge_set, g, source, sink, k);

				if ( k == key_edge_set->A_num)
				{
					j = k+1;
					break;
				}				
			}
			//ʹ��ð������һ��ChangeAmount_p1
			Bubbling_KeyEdge_A_1(key_edge_set, i, j-1);
			i=j;
			j=i+1;
		}
	}
	return;
}


//////////////////////////////////////////////////////////////////////////
//B1��߼���

//����һ��ʹ��������ʾ��ͼ�������
/*���룺vector<int>��ʾ��һ����ͼ,removedEdge��ʾ���Ƴ��ıߣ����Ϊ0�Ļ���ô�ͱ�ʾû���Ƴ������أ�����ͼ�������*/
int DinicByVertor(vector<int> a, Graph& g, int removedEdge, int source, int sink)
{
	/*ͨ��vector����ͼ*/
	Graph new_g;
	new_g.nV = g.nV;
	if (0==removedEdge)
	{
		new_g.nE = g.nE;
	}else{
		new_g.nE = g.nE -1;
	}

	/*copy������Ϣ*/
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
				aa = g.AllEdge[i][1];//��ʼ��
				bb = g.AllEdge[i][2];//�յ�
				cc = g.AllEdge[i][3];//����
				dd = g.AllEdge[i][4];//��Чλ���˴���ʹ�ã�

				ff = g.AllEdge_p[i];//�ɿ���

				new_g.matrix[aa][bb].iLabel = i;
				new_g.matrix[aa][bb].iC = cc;
				new_g.matrix[aa][bb].dP = ff;
			}
		}
	}

	/*ʹ��ԭ��Dinic�㷨���������*/
	Flow f;
	GF gf;
	return Dinic(new_g, source, sink, f, gf);
}

//����ȡ��һ�������ʾΪһ������
vector<int> getVertorByC(StateSet C0)
{
	vector<int> v;
	v.push_back(0);

	//��������дΪһ��vector
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
/* ����״̬�ռ�����½���л��֣���ȡ���������������                        */
/************************************************************************/
/*�������Ͻ죬�½磬ʹ��vector��ʾ�ģ������ȡ����һ��set��ʾ���Ƿָ�֮�����ͼ���䣬
�����ڵı�һ��Ҫʹ��0���棬set��Ԫ����set<vector<int>> */
void StateDivision(Graph& g, vector<int> down_, vector<int> top_, int source, int sink, set<vector<int>>& newState)
{	
	int nE = g.nE;                                                                /*������С�ǲ����*/
	double dR = 0.0;                                                              /*��ɿ�������ֲ�����*/

	Flow Fd;  
	GF gf;/*ʣ��ͼ*/
	int Fmax = g.max_flow;/*�õ����б߶����ڵ������*/
	assert(Fmax >= 0);/*��֤���������0������*/

	/*���һ����Ԫ��ĳ�ʼ��*/
	Collection c(down_, top_, nE);
	stack<Collection> cStack;                                                     /*�����Ҫ��һ�����ֵ���ͼ�ռ�*/
	cStack.push(c);

	Graph cur_g;
	StateSet Cj(nE),C0(nE);                                                       /*��ŵ�ǰ��Ҫ�����ıպ�����[]*/
	Collection next_c(nE);                                                        /*������һ����Ҫ��ջ������collection*/

	int tmpF = 0;
	double tmpP = 0.0;
	int * lofMaxP = new int[nE+1];                                                 /*��ʱ�洢����Ҫ��ıպ������и��������½�*/
	int * x0 = new int[nE+1];                                                      /*������*/
	memset(lofMaxP,0,(nE+1)*sizeof(int));
	memset(x0,0,(nE+1)*sizeof(int));
	while(!cStack.empty())
	{	
		GetCurrentC_j_new(cStack.top(),Cj);                                            /*ȡ�õ�ǰ��Ҫ�����ıպ�����[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*�任����һ����������״̬*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G_new(g,cur_g,Cj.lower);                                                 /*�½�������Ӧ����ͼ*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf);
		if(tmpF >= Fmax)                                                           /*�걸�����½��Ӧ��ͼ�ܹ����������*/
		{
			newState.insert(getVertorByC(Cj));
		}  
		else if((FromG2G_new(g,cur_g,Cj.upper),                                        /*���Ͻ��Ӧ����ͼ�������*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*ͬʱ����F(lc) < Fmax <= F(uc)�����
		 ����Chin-Chia Jane and Yih-Wenn Laih�ķ�������*/
			memset(x0,0,(nE+1)*sizeof(int));   /*��ʼ��������x0*/
			ComputeX0_new(g,Fd,x0);/*��û�����*/
			GetNextCollection_new(Cj,x0,next_c);
			cStack.push(next_c);
			GetCurrentC0_new(next_c,C0);/*ͨ������֪��C0�ض��ܹ�ȡ�������*/
			//��������궨�岿��
			//��C0���뵽�µ�״̬�����У���Щ�����������������
			newState.insert(getVertorByC(C0));
		}
	}

	delete[] x0;                                                                   /*�ͷ�����Ŀռ�*/
	delete[] lofMaxP;
	return;
}

//���������½���ͼ����
/*���룺vector<int>��ʾ��һ����ͼ,remvedEdge��ʾ��ɾ���ıߣ����Ϊ�����أ�����ͼ�ĸ���, ֻʹ���½�*/
double computeSubgraphProbabilityByVector_lower(vector<int> a, int remvedEdge, Graph &g)
{
	double p = 1;
	//�½���ڱߵĸ���֮��
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

/*��ȡ�����������������ͼ�����еļ�С��ͼ�ĸ���*/
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

//����һ���������ͼ����֮��
/*���룺vector<int>��ʾ��һ����ͼ,remvedEdge��ʾ��ɾ���ıߣ����Ϊ�����أ�����ͼ�ĸ���*/
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

/*����һЩ����ͼ�ĸ���֮��*/
/*���룺set<vector<int>>��ʾ���е���ͼ����,removedEdge��ʾ��Ҫ�Ƴ��ıߣ������������ͼ����֮��*/
double getP2ByStateSection(set<vector<int>> new_set, Graph &g, int removedEdge)
{
	double p = 0;

	set<vector<int>>::iterator it;
	for (it=new_set.begin(); it!=new_set.end(); it++)
	{
		p += computeSubgraphProbabilityByVector(*it, removedEdge, g);
	}
	return p;
}

//����2���ߵ�λ��  B1��
void Exchange_KeyEdge_B(KeyEdgeSet *key_edge_set,int i,int j)
{
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c  = key_edge_set->B_EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.ChangeAmount_p2 = key_edge_set->B_EdgeInfo[i].ChangeAmount_p2;
	TempKeyEdge.Edge			= key_edge_set->B_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class		= key_edge_set->B_EdgeInfo[i].Edge_Class;

	key_edge_set->B_EdgeInfo[i].ChangeAmount_c		=  key_edge_set->B_EdgeInfo[j].ChangeAmount_c	;
	key_edge_set->B_EdgeInfo[i].ChangeAmount_p2	=  key_edge_set->B_EdgeInfo[j].ChangeAmount_p2;
	key_edge_set->B_EdgeInfo[i].Edge				=  key_edge_set->B_EdgeInfo[j].Edge			;
	key_edge_set->B_EdgeInfo[i].Edge_Class			=  key_edge_set->B_EdgeInfo[j].Edge_Class		;

	key_edge_set->B_EdgeInfo[j].ChangeAmount_c		=  TempKeyEdge.ChangeAmount_c  ;
	key_edge_set->B_EdgeInfo[j].ChangeAmount_p2	=  TempKeyEdge.ChangeAmount_p2 ;
	key_edge_set->B_EdgeInfo[j].Edge				=  TempKeyEdge.Edge			   ;
	key_edge_set->B_EdgeInfo[j].Edge_Class			=  TempKeyEdge.Edge_Class	   ;

	return;
}



//B1���㣬�ڼ�����ͼ�ϵ�b���
void sortKeyEdge_B(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	//���B1Ϊ0����ֻ��һ��B1Ҳ����Ҫ����
	if (key_edge_set->B_num <= 1)
	{
		return;
	}

	for (int i =1; i<= key_edge_set->B_num; i++)
	{
		//�ߺ�
		int edge_num = key_edge_set->B_EdgeInfo[i].Edge;
		double temp_P = 0;

		/*�����п�������������������ͼ������allStateSection����*/
		//����ԭ�����������������ͼ����
		for (int ii=1; ii<= StateMtrix->State_Num; ii++)
		{
			//����ڶϵ���λ��Ϊ0����2�Ļ�����ô�������䶼���������
			if (StateMtrix->State[ii][edge_num]==0 || StateMtrix->State[ii][edge_num]==2)
			{
				//��ȡ��λ��Ϊ0����ΪX���½������������
				if (temp_P < g.lower_graph_P[ii])
				{
					temp_P = g.lower_graph_P[ii];
				}
			}
		}

		//������ɿ��ԣ������ɿ��ԣ�Ϊ��������������ĸ���֮��
		key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = temp_P;
	}

	//��B1���� ð��
	for (int i =1; i< key_edge_set->B_num; i++)
	{
		for (int j = i+1; j<= key_edge_set->B_num; j++)
		{
			if (key_edge_set->B_EdgeInfo[j].ChangeAmount_p2 < key_edge_set->B_EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge_B(key_edge_set, i,j);
			}
		}
	}

	return;
}


//////////////////////////////////////////////////////////////////////////
//B2��߼���



//////////////////////////////////////////////////////////////////////////

//C���  ����Ҫ����
void sortKeyEdge_C(KeyEdgeSet *key_edge_set, Graph &g)
{
	//ÿһ��C��ߵķ�������֮��ģ�״̬���ı䣬ֱ�ӱ��漴��
	for (int i=1; i <= key_edge_set->C_num; i++)
	{
		key_edge_set->C_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p2 = g.max_p2;
	}
	return;
}

////////////////////////////////////////////////////////////////////
// ���¹��캯��

void computeICA_APA(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	//1.�ֱ���������
	int AllEdgeFlowDecrease[MAX_E_NUM][2];
	getALLEdgeFlowDecrease(&AllEdgeFlowDecrease[0], g, source, sink);
	//2.���� ���Է���  ��ΪA��B��C�����
	classifyABC(StateMtrix, key_edge_set, &AllEdgeFlowDecrease[0]);

	//3.����   ��������    ��   ����
	//A������¼���
	sortKeyEdge_A(key_edge_set, g, source, sink);
	//B�����������
	sortKeyEdge_B(key_edge_set, g, source, sink, StateMtrix);
	//C��߲���Ҫ����
	sortKeyEdge_C(key_edge_set, g);

	return;
}


//////////////////////////////////////////////////////////////////
//�����ǽӿڣ����ڸ����е���


/*���ؼ������*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"�ؼ���������£�"<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.A_num;i++)
	{
		out<<setw(5)<<key_edge_set.A_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.A_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.B_num;i++)
	{
		out<<setw(5)<<key_edge_set.B_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.C_num;i++)
	{
		out<<setw(5)<<key_edge_set.C_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.C_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_p2<<endl;
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