#ifndef GRAPH_H
#define GRAPH_H

#include "ConstDef.h"
#include <vector>
#include <set>
using namespace std;

struct Edge
{
	int iC;                                  /*边上容量*/
	double dP;                               /*边可靠度*/
	int iLabel;                              /*边的编号*/

	Edge():iC(0),dP(0),iLabel(0){}
};

struct CP
{
	int iC; 
	double dP;
};                                           /*存放容量与概率*/

typedef Edge AdjMatrix[MAX][MAX];            /*零行零列不用*/

class Graph
{
public:
	Graph();
	~Graph();

	void Init();                             /*图清零*/


public: 
	int nId;                                 /*图编号*/
	int nV;                                  /*顶点数*/
	int nE;                                  /*边数*/

public:
	AdjMatrix matrix;                        /*邻接矩阵*/
	/*存储边的信息,不使用0*/
	int AllEdge[MAX_E_NUM][5];/*1234，起点，终点，容量，边的标号*/
	double AllEdge_p[MAX_E_NUM];/**/

	int max_flow;//最大流
	double max_p1;//流分布可靠性
	double max_p2;//网络可靠性
};


typedef int Flow[MAX][MAX];                  /*存放流的矩阵类型(0行0列不用)*/
typedef int GF[MAX][MAX];                    /*剩余图类型(0行0列不用)*/

struct MSN_STRUCT                            /*分成网络类型*/
{
	int nMSN[MAX][MAX];                      /*存放分层网络边的容量(0行0列不用)*/
	int nVStage[MAX];                        /*记录分层网络中顶点的层次(0行0列不用)*/
	int nV;                                  /*顶点个数(与原网络顶点个数相同)*/

	void Init(int nVnum);
}; 

typedef struct MSN_STRUCT MSN;      

struct BFnode
{
	int u;
	int v;
	int flow;
};
typedef struct BFnode BFLOW;                 /*阻塞流*/            

/*Dinic算法*/
int Dinic(Graph& g,int source,int sink,Flow& f,GF& gf);

/*通过剩余图得到相应最大流的最小割*/
void MinCut(GF& gf,int n,int source,set<int>& S);

#endif