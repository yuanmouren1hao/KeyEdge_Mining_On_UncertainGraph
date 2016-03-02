#ifndef RESULT_PRINT_H
#define RESULT_PRINT_H

#include "graph.h"
#include "partition.h"

/*程序状态与结果打印*/
void PrintGraph(Graph& g);                       /*打印图g*/

void PrintFlow(ofstream &out, Flow& f,int n);    /*打印流分布*/

void PrintMSN(MSN& msn);                         /*输出分层网络*/

void PrintX(int* x,int n);                       /*打印向量X*/

void PrintRunTime(long start_t,long end_t);      /*打印运行的时间*/

void PrintFmax_Prob(ostream &out,
					int s,
					int t,
					int fmax,
					double dp);                  /*打印最大流值与最可靠最大流的概率*/
void PrintTime(ostream &out,double cost);        /*打印时间*/


void PrintN(const char * fileName,
			int nIterator,
			int nDinic,
			Graph& g,
			char* algName);                      /* 向文件打印算法迭代次数和Dinic算法调用次数*/

void PrintCollection(Collection& c);             /*打印一个Collection*/

void PrintState(StateSet &statSet);

#endif