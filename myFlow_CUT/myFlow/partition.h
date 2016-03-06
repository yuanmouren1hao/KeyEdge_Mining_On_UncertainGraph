#ifndef SPA_H
#define SPA_H

#include "ConstDef.h"
#include "graph.h"
#include <vector>
#include <set>

/*引入状态的头文件*/
#include "state.h"
#include "cut.h"

using namespace std;

class StateSet
{
public:
	StateSet(int num);  
	~StateSet();

public:
	int numE;                            /*边数量*/
	int lower[MAX_E_NUM];                /*空间0不用*/
	int upper[MAX_E_NUM];                /*空间0不用*/
};  //一个具有上下界的状态集合

class Collection
{
public:
	Collection(int num, CertainEdge& certainEdge);
	Collection(int num);
	Collection(vector<int> down_, vector<int> top_, int num);
	~Collection();

public:
	int numE;                             /*边数量*/
	int lower[MAX_E_NUM];                 /*空间0不用*/
	int upper[MAX_E_NUM];                 /*空间0不用*/
	int Ad_C[MAX_E_NUM];                  /*图中满足x0[k] > lower[k]的边编号的集合*/
	int I;                                /*A_d_C集合中边的个数*/
	int j;                                /*当前C(j)*/
};                                        /*定义一个6元组用*/


/*通过划分的方式得到最可靠最大流分布(rusultFd为最可靠最大流分布,返回值为其概率)*/
double GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd,Lower_subGraph * StateMtrix, CertainEdge& certainEdge); 
/*通过对可能事件模型进行划分得到最可靠最大流分布*/
double old_GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd);

/*通过对可能事件模型进行划分得到最可靠最大流分布*/
double base_GetMPMF(Graph& g,int source,int sink,int &maxflow,Flow& resultFd, double &max_p2, int break_edge);
#endif