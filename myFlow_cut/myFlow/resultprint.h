#ifndef RESULT_PRINT_H
#define RESULT_PRINT_H

#include "graph.h"
#include "partition.h"

/*����״̬������ӡ*/
void PrintGraph(Graph& g);                       /*��ӡͼg*/

void PrintFlow(ofstream &out, Flow& f,int n);    /*��ӡ���ֲ�*/

void PrintMSN(MSN& msn);                         /*����ֲ�����*/

void PrintX(int* x,int n);                       /*��ӡ����X*/

void PrintRunTime(long start_t,long end_t);      /*��ӡ���е�ʱ��*/

void PrintFmax_Prob(ostream &out,
					int s,
					int t,
					int fmax,
					double dp);                  /*��ӡ�����ֵ����ɿ�������ĸ���*/
void PrintTime(ostream &out,double cost);        /*��ӡʱ��*/


void PrintN(const char * fileName,
			int nIterator,
			int nDinic,
			Graph& g,
			char* algName);                      /* ���ļ���ӡ�㷨����������Dinic�㷨���ô���*/

void PrintCollection(Collection& c);             /*��ӡһ��Collection*/

void PrintState(StateSet &statSet);

#endif