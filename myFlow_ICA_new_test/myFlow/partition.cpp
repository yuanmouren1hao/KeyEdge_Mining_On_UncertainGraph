#include "stdafx.h"

#include "partition.h"
#include "graph.h"
#include "resultprint.h"
#include <memory>
#include <stack>                                                                 /*ʹ��ջ���Ż��ռ�*/
#include <iostream>

/*����״̬��ͷ�ļ�*/
#include "state.h"

using namespace std;

StateSet::StateSet(int num):numE(num)
{	
}

StateSet::~StateSet()
{
}

/*���ڸ��ƿռ���Ż��ṹ*/
/*���캯����ʼ��*/
Collection::Collection(int num):numE(num),I(0),j(0)
{
	for(int i = 0; i <= numE; i++)
	{
		this->lower[i] = 0;  
		this->upper[i] = 1;
		/*��ʼ�����Ad_C����Ϊ��(��I��֤)*/
	}
}

/*���캯����ʼ�������Ͻ���½�Ĺ��캯��*/
Collection::Collection(vector<int> down_, vector<int>top_, int num)
{
	this->numE = num;
	this->I = 0;
	this->j = 0;
	for(int i = 0; i <= num; i++)
	{
		this->lower[i] = down_[i];  
		this->upper[i] = top_[i];
	}
}

Collection::~Collection()
{
}

/*��״̬����ת���ɶ�Ӧͼg_target*/
void FromG2G(Graph& g_original,Graph& g_target,int* lu)
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
void ComputeX0(Graph& g,Flow Fd,int* x0)
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

/*��ͼl����ͼ�ռ�����ͼ�����ĸ���*/
double CalculateP(Graph& g,int* l)
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

/*����Chin-Chia Jane and Yih-Wenn Laih�ķ����õ�C(j)(1-q)*/
void GetCurrentC_j(Collection &c,StateSet &Cj)
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

/*����Chin-Chia Jane and Yih-Wenn Laih�������C0*/
void GetCurrentC0(Collection &c,StateSet &C0)
{
	memcpy(C0.lower,c.lower,(C0.numE+1)*sizeof(int));
	memcpy(C0.upper,c.upper,(C0.numE+1)*sizeof(int));
	for(int i = 1; i <= c.I; i++)
	{
		C0.lower[c.Ad_C[i]] = C0.upper[c.Ad_C[i]];
	}
}

/*ͨ���պ������뻮���ߵõ�pivot����ad(c)*/
/*��������ad(c)���ϵĴ�С��Ad_c*/
int GetAd_C(StateSet &Cj,int *x0,int *Ad_c) 
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
void GetNextCollection(StateSet &Cj,int *x0,Collection &nextCollection)
{
	nextCollection.j = 1;                                                         /*�趨��Ҫ����Ҫ���Cj*/
	nextCollection.I = GetAd_C(Cj,x0,nextCollection.Ad_C);                        /*�趨I��Ad_c[]*/
	memcpy(nextCollection.lower,Cj.lower,(nextCollection.numE+1)*sizeof(int));    /*�趨lower[]*/
	memcpy(nextCollection.upper,Cj.upper,(nextCollection.numE+1)*sizeof(int));    /*�趨upper[]*/
}

/*���浱ǰ�ܹ�ȡ�������ֵ�ıպ��������½���ͼ���ʴ���*/
#define SAVEM(C) \
	tmpP = CalculateP(g,(C).lower); \
	if(tmpP > dR)  \
	{  dR = tmpP; memcpy(lofMaxP,(C).lower,(nE+1)*sizeof(int)); } 



/*�������е�״̬*/
void saveAllState(StateSet c,Lower_subGraph * StateMtrix)
{
	StateMtrix->State_Num++;
	StateMtrix->Edge_Num=c.numE;
	// ���½���ͼ������״̬������
	for (int ii=1;ii<=c.numE;ii++)
	{
		if (c.lower[ii] == 0 && c.upper[ii] == 0)
		{
			StateMtrix->State[StateMtrix->State_Num][ii] = 0;
		}else if (c.lower[ii] == 1 && c.upper[ii] == 1)
		{
			StateMtrix->State[StateMtrix->State_Num][ii] = 1;
		}else if (c.lower[ii] == 0 && c.upper[ii] == 1)
		{
			StateMtrix->State[StateMtrix->State_Num][ii] = 2;
		}
	}
	return;
}

	  
/*ͨ���Կ����¼�ģ�ͽ��л��ֵõ���ɿ�������ֲ�*/
double  old_GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd)
{	
	assert((source > 0 && source <= g.nV) && (sink > 0 && sink <= g.nV));         /*��֤����Ϸ�*/
	
	int nE = g.nE;                                                                /*������С�ǲ����*/
	double dR = 0.0;                                                              /*��ɿ�������ֲ�����*/

	Flow Fd;  
	GF gf;                                                                        /*ʣ��ͼ*/
	int Fmax = Dinic(g,source,sink,Fd,gf);                                        /*�õ����б߶����ڵ������*/
	assert(Fmax >= 0);                                                             /*��֤���������0������*/

	/*���һ����Ԫ��ĳ�ʼ��*/
	Collection c(nE);
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
		GetCurrentC_j(cStack.top(),Cj);                                            /*ȡ�õ�ǰ��Ҫ����ıպ�����[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*�任����һ���������״̬*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G(g,cur_g,Cj.lower);                                                 /*�½�������Ӧ����ͼ*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf); 
		if(tmpF >= Fmax)                                                           /*�걸�����½��Ӧ��ͼ�ܹ����������*/
		{
			SAVEM(Cj);
		}  
		else if((FromG2G(g,cur_g,Cj.upper),                                        /*���Ͻ��Ӧ����ͼ�������*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*ͬʱ����F(lc) < Fmax <= F(uc)�����
			����Chin-Chia Jane and Yih-Wenn Laih�ķ�������*/
			memset(x0,0,(nE+1)*sizeof(int));                                       /*��ʼ��������x0*/
			ComputeX0(g,Fd,x0);                                                    /*��û�����*/
			GetNextCollection(Cj,x0,next_c);
			cStack.push(next_c); 
			GetCurrentC0(next_c,C0);                                               /*ͨ������֪��C0�ض��ܹ�ȡ�������*/
			SAVEM(C0);                                                             /*ֱ�ӿ���C0�½���ͼ*/
		}
	}

	/*�����յ������������ȥ������һ���ߺ�
	��С��������ĸ�������½�������������㷨*/
	FromG2G(g,cur_g,lofMaxP); 

	maxflow = Dinic(cur_g,source,sink,resultFd,gf);                                /*resultFd�������ս��*/

	delete[] x0;                                                                   /*�ͷ�����Ŀռ�*/
	delete[] lofMaxP; 

	return dR;
}

/************************************************************************/
/* ������ͼ�������֮��*/
/* gΪ��ȷ��ͼ��cΪ��ͼ����
/************************************************************************/
double compute_all_subgraph(Graph &g, StateSet c)
{
	double p =1;
	for (int i=1; i<= c.numE; i++)
	{
		if (c.lower[i]==1 && c.upper[i]==1)
		{
			p *= g.AllEdge_p[i];
		}
		else if (c.lower[i]==0 && c.upper[i]==0)
		{
			p *= (1-g.AllEdge_p[i]);
		}
		else if (c.lower[i]==0 && c.upper[i]==1)
		{
			p *= 1;
		}
	}
	return p;
}



/*ͨ���Կ����¼�ģ�ͽ��л��ֵõ���ɿ�������ֲ�*/
double  GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd,Lower_subGraph * StateMtrix)
{	
	assert((source > 0 && source <= g.nV) && (sink > 0 && sink <= g.nV));         /*��֤����Ϸ�*/
	
	int nE = g.nE;                                                                /*������С�ǲ����*/
	double dR = 0.0;                                                              /*��ɿ�������ֲ�����*/

	Flow Fd;  
	GF gf;                                                                        /*ʣ��ͼ*/
	int Fmax = Dinic(g,source,sink,Fd,gf);                                        /*�õ����б߶����ڵ������*/
	assert(Fmax >= 0);                                                             /*��֤���������0������*/

	/*�������������ɿ��ԣ�����������������Ӹ���֮�ͣ���ʼ��Ϊ1*/
	double max_p2 = 0;

	/*���һ����Ԫ��ĳ�ʼ��*/
	Collection c(nE);
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
		GetCurrentC_j(cStack.top(),Cj);                                            /*ȡ�õ�ǰ��Ҫ����ıպ�����[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*�任����һ���������״̬*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G(g,cur_g,Cj.lower);                                                 /*�½�������Ӧ����ͼ*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf); 
		if(tmpF >= Fmax)                                                           /*�걸�����½��Ӧ��ͼ�ܹ����������*/
		{
			//SAVEM(Cj);
			/*������������һ��*/
			tmpP = CalculateP(g,(Cj).lower); 
			if(tmpP > dR)
			{
				dR = tmpP; 
				memcpy(lofMaxP,(Cj).lower,(nE+1)*sizeof(int));
				for (int ii =1; ii<= g.nE; ii++)
				{   //���²����漫С��ͼ
					StateMtrix->State[0][ii] = Cj.lower[ii];
				}
			} 
			//��������½���ͼ
			
			for (int ii=1;ii<=Cj.numE;ii++) {
				cout<<Cj.lower[ii]; }
			cout<<"--Cj-->";
			for (int ii=1;ii<=Cj.numE;ii++) {
				cout<<Cj.upper[ii]; }
			cout<<endl;
			
			//������ͼ����
			saveAllState(Cj,StateMtrix);
			//���ϸ���ͼ�����������ͼ�ĸ���֮��
			max_p2 += compute_all_subgraph(g, Cj);
		}  
		else if((FromG2G(g,cur_g,Cj.upper),                                        /*���Ͻ��Ӧ����ͼ�������*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*ͬʱ����F(lc) < Fmax <= F(uc)�����
			����Chin-Chia Jane and Yih-Wenn Laih�ķ�������*/
			memset(x0,0,(nE+1)*sizeof(int));                                       /*��ʼ��������x0*/
			ComputeX0(g,Fd,x0);                                                    /*��û�����*/
			GetNextCollection(Cj,x0,next_c);
			cStack.push(next_c); 
			GetCurrentC0(next_c,C0);                                               /*ͨ������֪��C0�ض��ܹ�ȡ�������*/
			//SAVEM(C0);                                                             /*ֱ�ӿ���C0�½���ͼ*/
			/*������������һ��*/
			tmpP = CalculateP(g,(C0).lower); 
			if(tmpP > dR)
			{
					dR = tmpP; 
					memcpy(lofMaxP,(C0).lower,(nE+1)*sizeof(int));
					for (int ii =1; ii<= g.nE; ii++)
					{   //���²����漫С��ͼ
						StateMtrix->State[0][ii] = C0.lower[ii];
					}
			} 

			//��������½���ͼ
			for (int ii=1;ii<=C0.numE;ii++) {
				cout<<C0.lower[ii]; }
			cout<<"--C0-->";
			for (int ii=1;ii<=C0.numE;ii++) {
				cout<<C0.upper[ii]; }
			cout<<endl; 
			//������ͼ����
			saveAllState(C0,StateMtrix);
			//���ϸ���ͼ�����������ͼ�ĸ���֮�ͣ��������ɿ���
			max_p2 += compute_all_subgraph(g, C0);
		}
	}
	cout<<endl;

	/*�����յ������������ȥ������һ���ߺ�
	��С��������ĸ�������½�������������㷨*/
	FromG2G(g,cur_g,lofMaxP); 

	maxflow = Dinic(cur_g,source,sink,resultFd,gf);                                /*resultFd�������ս��*/
	//���������
	g.max_flow = maxflow;
	/*��ͼ�����б���ԭʼ��  �����ɿ���   ��  �ֲ��ɿ���*/
	g.max_p2 = max_p2;
	g.max_p1 = dR;

	delete[] x0;                                                                   /*�ͷ�����Ŀռ�*/
	delete[] lofMaxP;
	return dR;
}

/************************************************************************/
/* ������ͼ�������֮��*/
/* gΪ��ȷ��ͼ��cΪ��ͼ���䣬iΪ���ϵı�
/************************************************************************/
double compute_subgraph(Graph &g, StateSet c, int ii)
{
	double p =1;
	for (int i=1; i<= c.numE; i++)
	{
		//������Ƕϵ��ıߣ��ϵ��ı߲���Ҫ����
		if (ii != i)
		{
			if (c.lower[i]==1 && c.upper[i]==1)
			{
				p *= g.AllEdge_p[i];
			}
			else if (c.lower[i]==0 && c.upper[i]==0)
			{
				p *= (1-g.AllEdge_p[i]);
			}
			else if (c.lower[i]==0 && c.upper[i]==1)
			{
				p *= 1;
			}
		}
	}
	return p;
}

/*ͨ���Կ����¼�ģ�ͽ��л��ֵõ���ɿ�������ֲ�*/
double base_GetMPMF(Graph& g, int source, int sink, int &maxflow, Flow& resultFd, double &max_p2, int break_edge)
{	
	assert((source > 0 && source <= g.nV) && (sink > 0 && sink <= g.nV));         /*��֤����Ϸ�*/
	
	/*������Ҫ�޸�һ�±ߵĳ��ȣ���Ϊ��������֮��߼���һ����������Ȼ��һ����λ*/
	int nE = g.nE+1;                                                                /*������С�ǲ����*/
	double dR = 0.0;                                                              /*��ɿ�������ֲ�����*/

	Flow Fd;  
	GF gf;                                                                        /*ʣ��ͼ*/
	int Fmax = Dinic(g,source,sink,Fd,gf);                                        /*�õ����б߶����ڵ������*/
	assert(Fmax >= 0);                                                             /*��֤���������0������*/

	max_p2 = 0;//�������������ɿ��ԣ�����������������Ӹ���֮�ͣ���ʼ��Ϊ1

	/*���һ����Ԫ��ĳ�ʼ��*/
	Collection c(nE);
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
		GetCurrentC_j(cStack.top(),Cj);                                            /*ȡ�õ�ǰ��Ҫ����ıպ�����[]Cj*/
		if(cStack.top().j < cStack.top().I)
		{
			cStack.top().j++;                                                      /*�任����һ���������״̬*/
		}
		else
		{
			cStack.pop();
		}

		FromG2G(g,cur_g,Cj.lower);                                                 /*�½�������Ӧ����ͼ*/
		tmpF = Dinic(cur_g,source,sink,Fd,gf); 
		if(tmpF >= Fmax)                                                           /*�걸�����½��Ӧ��ͼ�ܹ����������*/
		{
			SAVEM(Cj);
			//���ϸ���ͼ�����������ͼ�ĸ���֮��
			max_p2 += compute_subgraph(g, Cj, break_edge);
		}  
		else if((FromG2G(g,cur_g,Cj.upper),                                        /*���Ͻ��Ӧ����ͼ�������*/
			Dinic(cur_g,source,sink,Fd,gf)) >= Fmax) 
		{/*ͬʱ����F(lc) < Fmax <= F(uc)�����
			����Chin-Chia Jane and Yih-Wenn Laih�ķ�������*/
			memset(x0,0,(nE+1)*sizeof(int));                                       /*��ʼ��������x0*/
			ComputeX0(g,Fd,x0);                                                    /*��û�����*/
			GetNextCollection(Cj,x0,next_c);
			cStack.push(next_c); 
			GetCurrentC0(next_c,C0);                                               /*ͨ������֪��C0�ض��ܹ�ȡ�������*/
			SAVEM(C0);                                                             /*ֱ�ӿ���C0�½���ͼ*/
			//���ϸ���ͼ�����������ͼ�ĸ���֮��
			max_p2 += compute_subgraph(g, C0, break_edge);
		}
	}

	/*�����յ������������ȥ������һ���ߺ�
	��С��������ĸ�������½�������������㷨*/
	FromG2G(g,cur_g,lofMaxP); 

	maxflow = Dinic(cur_g,source,sink,resultFd,gf);                                /*resultFd�������ս��*/

	delete[] x0;                                                                   /*�ͷ�����Ŀռ�*/
	delete[] lofMaxP; 

	return dR;
}