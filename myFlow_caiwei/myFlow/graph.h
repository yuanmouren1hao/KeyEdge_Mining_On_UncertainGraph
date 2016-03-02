#ifndef GRAPH_H
#define GRAPH_H

#include "ConstDef.h"
#include <vector>
using namespace std;

struct Edge
{
	int iC;                                  /*��������*/
	double dP;                               /*�߿ɿ���*/
	int iLabel;                              /*�ߵı��*/

	Edge():iC(0),dP(0),iLabel(0){}
};

struct CP
{
	int iC; 
	double dP;
};                                           /*������������*/

typedef Edge AdjMatrix[MAX][MAX];            /*�������в���*/

class Graph
{
public:
	Graph();
	~Graph();

	void Init();                             /*ͼ����*/


public: 
	int nId;                                 /*ͼ���*/
	int nV;                                  /*������*/
	int nE;                                  /*����*/

public:
	AdjMatrix matrix;                        /*�ڽӾ���*/
};


typedef int Flow[MAX][MAX];                  /*������ľ�������(0��0�в���)*/
typedef int GF[MAX][MAX];                    /*ʣ��ͼ����(0��0�в���)*/

struct MSN_STRUCT                            /*�ֳ���������*/
{
	int nMSN[MAX][MAX];                      /*��ŷֲ�����ߵ�����(0��0�в���)*/
	int nVStage[MAX];                        /*��¼�ֲ������ж���Ĳ��(0��0�в���)*/
	int nV;                                  /*�������(��ԭ���綥�������ͬ)*/

	void Init(int nVnum);
}; 

typedef struct MSN_STRUCT MSN;      

struct BFnode
{
	int u;
	int v;
	int flow;
};
typedef struct BFnode BFLOW;                 /*������*/            

/*Dinic�㷨*/
int Dinic(Graph& g,int source,int sink,Flow& f,GF& gf);

/*ͨ��ʣ��ͼ�õ���Ӧ���������С��*/
void MinCut(GF& gf,int n,int source,vector<int>& S);

#endif