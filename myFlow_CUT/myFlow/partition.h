#ifndef SPA_H
#define SPA_H

#include "ConstDef.h"
#include "graph.h"
#include <vector>
#include <set>

/*����״̬��ͷ�ļ�*/
#include "state.h"

using namespace std;

class StateSet
{
public:
	StateSet(int num);  
	~StateSet();

public:
	int numE;                            /*������*/
	int lower[MAX_E_NUM];                /*�ռ�0����*/
	int upper[MAX_E_NUM];                /*�ռ�0����*/
};  //һ���������½��״̬����

class Collection
{
public:
	Collection(int num, set<int>& MinCutEdges);
	Collection(int num);
	~Collection();

public:
	int numE;                             /*������*/
	int lower[MAX_E_NUM];                 /*�ռ�0����*/
	int upper[MAX_E_NUM];                 /*�ռ�0����*/
	int Ad_C[MAX_E_NUM];                  /*ͼ������x0[k] > lower[k]�ı߱�ŵļ���*/
	int I;                                /*A_d_C�����бߵĸ���*/
	int j;                                /*��ǰC(j)*/
};                                        /*����һ��6Ԫ����*/


/*ͨ�����ֵķ�ʽ�õ���ɿ�������ֲ�(rusultFdΪ��ɿ�������ֲ�,����ֵΪ�����)*/
double GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd,Lower_subGraph * StateMtrix, set<int>& MinCutEdges); 
/*ͨ���Կ����¼�ģ�ͽ��л��ֵõ���ɿ�������ֲ�*/
double old_GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd);

/*ͨ���Կ����¼�ģ�ͽ��л��ֵõ���ɿ�������ֲ�*/
double base_GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd, double &max_p2, int break_edge);
#endif