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

/*�ж�һ����ͼ�������*/
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
	StateSet Cj(nE),C0(nE);                                                       /*��ŵ�ǰ��Ҫ����ıպ�����[]*/
	Collection next_c(nE);                                                        /*������һ����Ҫ��ջ�����collection*/

	int tmpF = 0;
	double tmpP = 0.0;
	int * lofMaxP = new int[nE+1];                                                 /*��ʱ�洢����Ҫ��ıպ������и��������½�*/
	int * x0 = new int[nE+1];                                                      /*������*/
	memset(lofMaxP,0,(nE+1)*sizeof(int));
	memset(x0,0,(nE+1)*sizeof(int));
	while(!cStack.empty())
	{	
		GetCurrentC_j_new(cStack.top(),Cj);                                            /*ȡ�õ�ǰ��Ҫ����ıպ�����[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*�任����һ���������״̬*/
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


/************************************************************************/
/* ����Ϊ��ʵ�鲿��                                                       */
/************************************************************************/
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

/*����һ��A���*/
void computetSimgleA(Graph& g, int source, int sink, int EdgeNum, KeyEdgeSet& AllKeyEdge)
{
	int new_maxflow = 0; /*��  �����*/
	double new_p1 =0;/*��  �ֲ��ɿ���*/
	double new_p2 =0;/*��  �����ɿ���*/
	Flow new_maxPmaxF;/*�µ�������ֲ�*/
	Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/

	/*��ʼ����ʱ�洢��*/
	init_TempEdge(TempEdge);
	/*ȥ��ĳһ��A���*/
	removeA_Edge(g, EdgeNum, TempEdge);
	//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
	new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, EdgeNum);

	/* ���ϵ���֮�����������ֲ��ɿ��ԣ������ɿ��Էֱ𱣴� */
	int now_edge_num = AllKeyEdge.A_num;
	AllKeyEdge.A_EdgeInfo[now_edge_num].ChangeAmount_c = new_maxflow;
	AllKeyEdge.A_EdgeInfo[now_edge_num].ChangeAmount_p1 = new_p1;
	AllKeyEdge.A_EdgeInfo[now_edge_num].ChangeAmount_p2 = new_p2;
	AllKeyEdge.A_EdgeInfo[now_edge_num].Edge = EdgeNum;
	AllKeyEdge.A_EdgeInfo[now_edge_num].Edge_Class = 'A';

	/* �ָ�ĳһ��A��� */
	restoreA_Edge(g,EdgeNum, TempEdge);

	return;
}

/*���������ߵ�λ��*/
void Exchange_KeyEdge(KeyEdgeSet& AllKeyEdge,int i,int j)
{
	KeyEdge TempKeyEdge;
	TempKeyEdge.ChangeAmount_c = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_c;
	TempKeyEdge.ChangeAmount_p1 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1;
	TempKeyEdge.ChangeAmount_p2 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2;
	TempKeyEdge.Edge = AllKeyEdge.A_EdgeInfo[i].Edge;
	TempKeyEdge.Edge_Class = AllKeyEdge.A_EdgeInfo[i].Edge_Class;

	AllKeyEdge.A_EdgeInfo[i].ChangeAmount_c = AllKeyEdge.A_EdgeInfo[j].ChangeAmount_c;
	AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1 = AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p1;
	AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2 = AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p2;
	AllKeyEdge.A_EdgeInfo[i].Edge = AllKeyEdge.A_EdgeInfo[j].Edge;
	AllKeyEdge.A_EdgeInfo[i].Edge_Class = AllKeyEdge.A_EdgeInfo[j].Edge_Class;

	AllKeyEdge.A_EdgeInfo[j].ChangeAmount_c = TempKeyEdge.ChangeAmount_c;
	AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p1 = TempKeyEdge.ChangeAmount_p1;
	AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p2 = TempKeyEdge.ChangeAmount_p2;
	AllKeyEdge.A_EdgeInfo[j].Edge = TempKeyEdge.Edge;
	AllKeyEdge.A_EdgeInfo[j].Edge_Class = TempKeyEdge.Edge_Class;

	return;
}

/*ʹ��ð�ݷ�����һ������,���ݵ�������ֲ��ɿ���һ�£��������ɿ���һ��
/************************************************************************/
/* AllKeyEdge�������йؼ��ߣ�iΪ��ʼλ�ã�jΪ����λ��                       */
/************************************************************************/
void Bubbling_KeyEdge_2(KeyEdgeSet& AllKeyEdge, int begin_, int end_)
{
	/*��Ҫ����i<j, ����i��j����������*/
	if (begin_ > end_)
	{
		return;
	}

	for (int i=begin_; i<end_; i++)
	{
		for (int j=begin_+1; j<=end_; j++)
		{
			if (AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p2 < AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2)
			{
				Exchange_KeyEdge(AllKeyEdge, i, j);
			}
		}
	}
	return;
}

/*����һ��A���*/
void sortPartA(KeyEdgeSet& AllKeyEdge, int _from, int _to)
{
	/*��Ҫ����i<j, ����i��j����������*/
	if (_from > _to)
	{
		return;
	}

	/*�ȽϷֲ��ɿ���*/
	for (int i=_from; i<_to; i++)
	{
		for (int j=_from+1; j<=_to; j++)
		{
			if (AllKeyEdge.A_EdgeInfo[j].ChangeAmount_p1 < AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1)
			{
				Exchange_KeyEdge(AllKeyEdge, i, j);
			}
		}
	}

	/*�ֲ��ɿ���һ�£��Ƚ������ɿ���*/
	int iii = _from;
	int jjj = iii+1;
	while(jjj <= _to)
	{
		if (AllKeyEdge.A_EdgeInfo[jjj].ChangeAmount_p1 != AllKeyEdge.A_EdgeInfo[iii].ChangeAmount_p1)
		{
			iii++;
			jjj++;
		}else
		{
			//�ҵ���ͬ�仯������
			for (int kkk =jjj; kkk <= AllKeyEdge.EdgeNum; kkk++)
			{
				if (AllKeyEdge.A_EdgeInfo[kkk].ChangeAmount_p1 != AllKeyEdge.A_EdgeInfo[iii].ChangeAmount_p1)
				{
					jjj = kkk;
					break;
				}
			}
			//ʹ��ð������һ��
			Bubbling_KeyEdge_2(AllKeyEdge, iii, jjj-1);
			iii=jjj;
			jjj=iii+1;
		}

	}

	return;
}

/*��A������¼���*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink, int EdgeFlow[][2])
{
	/*���û��A��ߵĸ���Ϊ0*/
	if (key_edge_set->A_num==0)
	{
		return;
	}

	/*����һ���Ѿ��ź���Ĺؼ������*/
	KeyEdgeSet AllKeyEdge;
	AllKeyEdge.A_num = 0;
	AllKeyEdge.EdgeNum = key_edge_set->EdgeNum;
	int a_num = key_edge_set->A_num;
	
	/*��ǰ�Ѿ�����ߵĸ���*/
	int i=1;
	int j=i+1;
	while(i<=a_num)
	{
		/*��������һ��*/
		if (i == a_num)
		{
			/*i����߲���Ҫ���㣬ֱ�ӱȽϳ���*/
			AllKeyEdge.A_num++;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge_Class = 'A';
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge = EdgeFlow[i][0];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_c = EdgeFlow[i][1];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p1 = 0;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p2 = 0;
			break;
		}

		if (EdgeFlow[i][1] < EdgeFlow[j][1])
		{
			/*i����߲���Ҫ���㣬ֱ�ӱȽϳ���*/
			AllKeyEdge.A_num++;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge_Class = 'A';
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].Edge = EdgeFlow[i][0];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_c = EdgeFlow[i][1];
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p1 = 0;
			AllKeyEdge.A_EdgeInfo[AllKeyEdge.A_num].ChangeAmount_p2 = 0;

			/*i��j�Ĳ�����+1*/
			i++;
			j++;

		}else if (EdgeFlow[i][1] = EdgeFlow[j][1])
		{
			/*������2���������ظ�����Ҫ�����ظ���   �����ҵ��м����ظ���*/
			while(EdgeFlow[j][1] == EdgeFlow[i][1])
			{
				j++;
			}

			/*����ÿһ����������Ҫ����ɿ���  �ظ�����*/
			for (int ii =i; ii<=j-1;ii++)
			{
				/*��һ��A��ߵ���Ϣ��������� ���浽AllKeyEdge��*/
				AllKeyEdge.A_num++;
				computetSimgleA(g, source, sink, EdgeFlow[ii][0], AllKeyEdge);
			}
			/*�������֮��������һ��,������AllKeyEdge*/
			sortPartA(AllKeyEdge, i, j-1);

			/*���Ӳ���*/
			i=j;
			j++;
		}
	}

	/*��AllKeyEdge ��A���ֱ��浽key_edge_set��*/
	for (int i=1; i<=AllKeyEdge.A_num; i++)
	{
		key_edge_set->A_EdgeInfo[i].ChangeAmount_c = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_c;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p1 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p1;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = AllKeyEdge.A_EdgeInfo[i].ChangeAmount_p2;
		key_edge_set->A_EdgeInfo[i].Edge_Class = AllKeyEdge.A_EdgeInfo[i].Edge_Class;
		key_edge_set->A_EdgeInfo[i].Edge = AllKeyEdge.A_EdgeInfo[i].Edge;
	}

	return;
}


/////////////////////////////////////////////////////////////////////////
//���й���B��ߵĲ���

/*������ͼ����*/
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


/*�����½���ͼ����*/
/*���룺vector<int>��ʾ��һ����ͼ,remvedEdge��ʾ��ɾ���ıߣ����Ϊ�����أ�����ͼ�ĸ���, ֻʹ���½�*/
double computeSubgraphProbabilityByVector_lower(vector<int> a, int remvedEdge, Graph &g)
{
	double p = 1;
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

/*����һЩ����ͼ�ĸ���֮��*/
/*���룺set<vector<int>>��ʾ���е���ͼ����,removedEdge��ʾ��Ҫ�Ƴ��ıߣ������������ͼ����֮��*/
double computeProbabilitySum(set<vector<int>> new_set, Graph &g, int removedEdge)
{
	double p = 0;

	set<vector<int>>::iterator it;
	for (it=new_set.begin(); it!=new_set.end(); it++)
	{
		p += computeSubgraphProbabilityByVector(*it, 0, g);
	}
	return p;
}

/*�Ƚ�����ƫ��ԵĴ�С*/
/*b������a����1��b>=a����2*/
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

/*�������е�״̬�����ȡ���е������������״̬*/
set<vector<int>> getAllState(Lower_subGraph* StateMtrix)
{
	/*�������е�����ŵ�ջ��*/
	stack<vector<int>> s;

	for (int i=1; i<= StateMtrix->State_Num; i++)
	{
		/*����һ����ͼ����*/
		vector<int> v;
		v.push_back(0);/*Ϊ������0��Ԫ�أ��ӵ�һ��Ԫ�ؿ�ʼʹ��*/

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

	/*��������״̬*/
	set<vector<int>> allState ;
	/*��ջ�е�״̬���зֽ�*/
	while (!s.empty())
	{
		vector<int> a = s.top();
		s.pop();
		
		for (int i=1; i<=StateMtrix->Edge_Num; i++)
		{
			if (a[i]==2)
			{
				/*�������X��Ҳ����2������ͼ����ֽ⣬�ֱ�����ջ��*/
				vector<int> temp_V = a;
				temp_V[i] =0;
				s.push(temp_V);
				temp_V[i] =1;
				s.push(temp_V);
				break;
			}
			
			/*���ֱ��û�еĻ���ֱ�Ӳ��뵽 ������*/
			if (i==StateMtrix->Edge_Num)
			{
				allState.insert(a);		
			}
		}
	}
	return allState;
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

/*B����������㣬������ʹ�����õķ�ʽ*/
void sortKeyEdge_B(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	/*���û��B��ߵĸ���Ϊ0*/
	if (key_edge_set->B1_num == 0)
	{
		return;
	}

	/*��ȡ���е������������״̬����ͼ���䣬ע������ͼ���䣬2��ʾx����0����1*/
	set<vector<int>> allStateSection;

	for (int i =1; i<= key_edge_set->B1_num; i++)
	{
		/*�ϵ�B1������������*/
		key_edge_set->B1_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		/*�ߺ�*/
		int edge_num = key_edge_set->B1_EdgeInfo[i].Edge;

		/*�����п�������������������ͼ������allStateSection����*/
		//����ԭ�����������������ͼ����
		for (int ii=1; ii<= StateMtrix->State_Num; ii++)
		{
			//����ڶϵ���λ��Ϊ0����2�Ļ�����ô�������䶼���������
			if (StateMtrix->State[ii][edge_num]==0)
			{
				/*����һ����ͼ����*/
				vector<int> v;
				v.push_back(0);/*Ϊ������0��Ԫ�أ��ӵ�һ��Ԫ�ؿ�ʼʹ��*/
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					//����Ƕϵ��ıߵĻ���ֱ�ӱ���Ϊ0
					if (j==edge_num)
					{
						v.push_back(0);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
					{
						v.push_back(0);
					}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
						v.push_back(1);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
					{
						v.push_back(2);
					}
				}
				allStateSection.insert(v);
			}
			//�����䲿���������������Ҫ��һ�����ж�
			else if (StateMtrix->State[ii][edge_num]==1)
			{
				//�½�
				vector<int>TempV_low;
				TempV_low.push_back(0);
				//����ϵ�֮�� �½���Ȼ�������������ô�������䶼����
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_low.push_back(0);
					}
					else
					{
						TempV_low.push_back(StateMtrix->State[ii][j]);
					}
				}
				//�Ͻ�
				vector<int>TempV_upper;
				TempV_upper.push_back(0);
				//����ϵ�֮�������Ͻ첻�������������ô�������䶼������
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_upper.push_back(0);
					}
					else
					{
						TempV_upper.push_back(StateMtrix->State_upper[ii][j]);
					}
				}

				//����½����������,��������ֱ�ӱ�����allStateSection��
				if (DinicByVertor(TempV_low,g,edge_num,source,sink) == g.max_flow)
				{
					vector<int> v;
					v.push_back(0);/*Ϊ������0��Ԫ�أ��ӵ�һ��Ԫ�ؿ�ʼʹ��*/
					for (int j=1; j<=StateMtrix->Edge_Num; j++)
					{
						//����Ƕϵ��ıߵĻ���ֱ�ӱ���Ϊ0
						if (j==edge_num)
						{
							v.push_back(0);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
						{
							v.push_back(0);
						}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
							v.push_back(1);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
						{
							v.push_back(2);
						}
					}
					allStateSection.insert(v);
				}
				//����Ͻ���ͼ���������������ô�������䶼�������������ֱ����������
				else if (DinicByVertor(TempV_upper,g,edge_num,source,sink) < g.max_flow)
				{
					//do nothing
					//cout<<"can subgraph allowed."<<endl;
				}
				//�½粻������������Ͻ�����������Ļ�����Ҫ���·ֽ���ͼ���䣬���Ƿֽ����ͼ����Ƚ�С����
				else
				{
					StateDivision(g,TempV_low,TempV_upper,source,sink,allStateSection);
				}
			}
		}

		//�ҵ��������������������Ļ���ֱ�ӻ�ȡ��С��ͼ��������ɿ��Լ���
		key_edge_set->B1_EdgeInfo[i].ChangeAmount_p1 = getP1ByStateSection(allStateSection, g);
		/*������ɿ��ԣ������ɿ��ԣ�Ϊ��������������ĸ���֮��*/
		key_edge_set->B1_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(allStateSection, g, edge_num);
	}
	return;
}

/************************************************************************/
/*for edge b  ����B2��ߣ��ϵ�֮�������ͷֲ��ɿ��Բ���  �����ڼ�С��ͼ��     */
/************************************************************************/
/*B����������㣬������ʹ�����õķ�ʽ*/
void sortKeyEdge_b(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	/*���û��B2��ߵĸ���Ϊ0*/
	if (key_edge_set->B2_num == 0)
	{
		return;
	}

	/*��ȡ���е������������״̬*/
	set<vector<int>> allStateSection ;

	for (int i =1; i<= key_edge_set->B2_num; i++)
	{
		/*�ϵ�B�����������䣬�ֲ��ɿ��Բ���*/
		key_edge_set->B2_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		key_edge_set->B2_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
		/*��¼�ߵ����*/
		int edge_num = key_edge_set->B2_EdgeInfo[i].Edge;


		/*�����п�������������������ͼ������allStateSection����*/
		//����ԭ�����������������ͼ����
		for (int ii=1; ii<= StateMtrix->State_Num; ii++)
		{
			//����ڶϵ���λ��Ϊ0����2�Ļ�����ô�������䶼���������
			if (StateMtrix->State[ii][edge_num]==0)
			{
				/*����һ����ͼ����*/
				vector<int> v;
				v.push_back(0);/*Ϊ������0��Ԫ�أ��ӵ�һ��Ԫ�ؿ�ʼʹ��*/
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					//����Ƕϵ��ıߵĻ���ֱ�ӱ���Ϊ0
					if (j==edge_num)
					{
						v.push_back(0);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
					{
						v.push_back(0);
					}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
						v.push_back(1);
					}
					else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
					{
						v.push_back(2);
					}
				}
				allStateSection.insert(v);
			}
			//�����䲿���������������Ҫ��һ�����ж�
			else if (StateMtrix->State[ii][edge_num]==1)
			{
				//�½�
				vector<int>TempV_low;
				TempV_low.push_back(0);
				//����ϵ�֮�� �½���Ȼ�������������ô�������䶼����
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_low.push_back(0);
					}
					else
					{
						TempV_low.push_back(StateMtrix->State[ii][j]);
					}
				}
				//�Ͻ�
				vector<int>TempV_upper;
				TempV_upper.push_back(0);
				//����ϵ�֮�������Ͻ첻�������������ô�������䶼������
				for (int j=1; j<=StateMtrix->Edge_Num; j++)
				{
					if (j==edge_num)
					{
						TempV_upper.push_back(0);
					}
					else
					{
						TempV_upper.push_back(StateMtrix->State_upper[ii][j]);
					}
				}

				//����½����������,��������ֱ�ӱ�����allStateSection��
				if (DinicByVertor(TempV_low,g,edge_num,source,sink) == g.max_flow)
				{
					vector<int> v;
					v.push_back(0);/*Ϊ������0��Ԫ�أ��ӵ�һ��Ԫ�ؿ�ʼʹ��*/
					for (int j=1; j<=StateMtrix->Edge_Num; j++)
					{
						//����Ƕϵ��ıߵĻ���ֱ�ӱ���Ϊ0
						if (j==edge_num)
						{
							v.push_back(0);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==0)
						{
							v.push_back(0);
						}else if (StateMtrix->State[ii][j]==1 && StateMtrix->State_upper[ii][j]==1){
							v.push_back(1);
						}
						else if (StateMtrix->State[ii][j]==0 && StateMtrix->State_upper[ii][j]==1 && j!=edge_num)
						{
							v.push_back(2);
						}
					}
					allStateSection.insert(v);
				}
				//����Ͻ���ͼ���������������ô�������䶼�������������ֱ����������
				else if (DinicByVertor(TempV_upper,g,edge_num,source,sink) < g.max_flow)
				{
					//do nothing
					//cout<<"can subgraph allowed."<<endl;
				}
				//�½粻������������Ͻ�����������Ļ�����Ҫ���·ֽ���ͼ���䣬���Ƿֽ����ͼ����Ƚ�С����
				else
				{
					StateDivision(g,TempV_low,TempV_upper,source,sink,allStateSection);
				}
			}
		}
		/*������ɿ��ԣ������ɿ��ԣ�Ϊ��������������ĸ���֮��*/
		key_edge_set->B2_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(allStateSection, g, edge_num);
	}

	return;
}



//////////////////////////////////////////////////////////////////
//C��ߵĲ���

/*����C��߽��м���,C��߲���Ҫ����*/
void sortKeyEdge_C(KeyEdgeSet *key_edge_set,Graph &g)
{
	/*ÿһ��C��ߵķ�������֮��ģ�״̬���ı䣬ֱ�ӱ��漴��*/
	for (int i=1; i <= key_edge_set->C_num; i++)
	{
		key_edge_set->C_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p2 = g.max_p2;
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
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num && StateMtrix->State[0][ii]==1)
		{
			//cout<<ii<<" ��B���"<<endl;
			key_edge_set->B1_num++;
			key_edge_set->B1_EdgeInfo[key_edge_set->B1_num].Edge=ii;
			key_edge_set->B1_EdgeInfo[key_edge_set->B1_num].Edge_Class='B';
		}
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num && StateMtrix->State[0][ii]==0)
		{
			//cout<<ii<<" ��B���"<<endl;
			key_edge_set->B2_num++;
			key_edge_set->B2_EdgeInfo[key_edge_set->B2_num].Edge=ii;
			key_edge_set->B2_EdgeInfo[key_edge_set->B2_num].Edge_Class='b';
		}
	}
}

/************************************************************************/
/* ����g�Ͷϵ��ı߻�ȡ�������                                             */
/************************************************************************/
int getMAXflowWithBreakEdge(Graph& g, int BreakEdge, int source, int sink)
{
	/*ʹ��ԭ����ͼ  �ȹ��죬��ȡ�������Ȼ��ָ�ԭ��ͼ������*/
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
	/*ʹ��ԭ��Dinic�㷨���������*/
	Flow f;
	GF gf;
	int max_flow_remine = Dinic(g, source, sink, f, gf);

	/*�������֮��  ��Ҫ�ָ�ԭ�е�ͼ*/
	g.nE++;
	g.matrix[aa][bb].iLabel = BreakEdge;
	g.matrix[aa][bb].iC = cc;
	g.matrix[aa][bb].dP = ff;

	/*��������*/
	return max_flow_remine;
}


/************************************************************************/
/* ����ÿ���߶ϵ�֮������������                                            */
/************************************************************************/
void getALLEdgeFlowDecrease(Graph& g, int EdgeFlow[][2], int source, int sink)
{
	vector<int> v;
	v.push_back(0);
	for (int i=1; i<=g.nE; i++)
	{
		v.push_back(0);
	}

	for (int i=1; i<=g.nE; i++)
	{
		EdgeFlow[i][0] = i;
		EdgeFlow[i][1] = getMAXflowWithBreakEdge(g, i, source, sink);
	}

	/*��EdgeFlow��������*/
	for (int i =1; i<g.nE; i++)
	{
		for(int j=i+1; j<=g.nE; j++)
		{
			if (EdgeFlow[j][1] < EdgeFlow[i][1])
			{
				int temp;

				temp = EdgeFlow[j][1];
				EdgeFlow[j][1] = EdgeFlow[i][1];
				EdgeFlow[i][1] = temp;

				temp = EdgeFlow[j][0];
				EdgeFlow[j][0] = EdgeFlow[i][0];
				EdgeFlow[i][0] = temp;
			}
		}
	}

	return;
}



/*ͨ��״̬�������ߵ����*/
void computeEdgeClass(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*�ȼ���ÿ���߶ϵ�֮����������*/
	int EdgeFlow[MAX_E_NUM][2];
	getALLEdgeFlowDecrease(g, &EdgeFlow[0], source, sink);
	
	/*���Է���������ȷ��ͼ�еı߷�ΪA��B1��B2��C�����*/
	sortKeyEdge(StateMtrix, key_edge_set);
	
	/*��������*/
	/*A������¼���*/
	sortKeyEdge_A(key_edge_set, g, source, sink, &EdgeFlow[0]);
	/*B�����������*/
	sortKeyEdge_B(key_edge_set, g, source, sink, StateMtrix);
	/*for class edge b */
	sortKeyEdge_b(key_edge_set, g, source, sink, StateMtrix);
	/*C��߲���Ҫ����*/
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
		out<<setw(5)<<key_edge_set.A_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.A_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.B1_num;i++)
	{
		out<<setw(5)<<key_edge_set.B1_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B1_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B1_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B1_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.B1_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.B2_num;i++)
	{
		out<<setw(5)<<key_edge_set.B2_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B2_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B2_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B2_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.B2_EdgeInfo[i].ChangeAmount_p2<<endl;
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
	key_edge_set.B1_num=0;
	key_edge_set.B2_num=0;
	key_edge_set.C_num=0;
	key_edge_set.EdgeNum=0;

	for (int i=0;i<MAX_E_NUM;i++)
	{
		key_edge_set.A_EdgeInfo[i].Edge=0;
		key_edge_set.B1_EdgeInfo[i].Edge=0;
		key_edge_set.B2_EdgeInfo[i].Edge=0;
		key_edge_set.C_EdgeInfo[i].Edge=0;

		key_edge_set.A_EdgeInfo[i].Edge_Class=0;
		key_edge_set.B1_EdgeInfo[i].Edge_Class=0;
		key_edge_set.B2_EdgeInfo[i].Edge_Class=0;
		key_edge_set.C_EdgeInfo[i].Edge_Class=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.B1_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.B2_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_c=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.B1_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.B2_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_p1 = 0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.B1_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.B2_EdgeInfo[i].ChangeAmount_p2 = 0;
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