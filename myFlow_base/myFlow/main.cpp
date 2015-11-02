// Network Flow.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "InputReader.h"
#include "graph.h"
#include "resultprint.h"
#include <iostream>
#include <string.h>
#include <windows.h>
#include  "psapi.h"                                           /*��������(���̵����������Psapi.lib)*/

/*����ʹ���½���ͼ��һЩ�ṹ�뷽��*/
#include "state.h"

using namespace std; 

/*��release�汾���ܹ��ɹ�����*/
#pragma  once
#pragma  message("Psapi.h --> Linking with Psapi.lib")
#pragma  comment(lib,"Psapi.lib")


/*ȫ�ֱ������е��½���ͼ״̬*/
Lower_subGraph StateMtrix;
/*�������еĹؼ�����Ϣ*/
KeyEdgeSet key_edge_set;
/*�洢�ߵ���Ϣ,��ʹ��0*/
int AllEdge[MAX_E_NUM][5];/*1234����㣬�յ㣬��������Чλ*/
double AllEdge_p[MAX_E_NUM];/**/


int _tmain(int argc, char* argv[])
{
	char stFileName[BUFFER_SIZE]     = WORK_SPACE;             /*Դ����*/
	char fileName[BUFFER_SIZE]       = WORK_SPACE;             /*���ݴ���ļ�*/
	char resultFileName[BUFFER_SIZE] = WORK_SPACE;             /*ʵ��������ļ�*/

	
	if (4 != argc)
	{
		cout<<"Command Params : "<<endl
			<<"\tSource_Sink File Name"<<endl
			<<"\tGraph data File Name"<<endl
			<<"\tResult File Name"<<endl;
		return 1;
	}

	strcat_s(stFileName,argv[1]);
	strcat_s(fileName,argv[2]);
	strcat_s(resultFileName,argv[3]);
	
	/*
	char argv1[BUFFER_SIZE]="data\\s-t\\test_graph_st.txt";
	char argv2[BUFFER_SIZE]="data\\mydata\\TestGraph.txt";
	char argv3[BUFFER_SIZE]="results\\new\\TestGraph.txt";
	
	strcat_s(stFileName,argv1);
	strcat_s(fileName,argv2);
	strcat_s(resultFileName,argv3);
	cout<<stFileName<<endl;
	*/
	
	ofstream out_result;
	out_result.open(resultFileName,ios::out|ios::app);         /*����ʵ����*/
	if(!out_result.is_open())
	{
		printf("open result file failed...\n");
		exit(1);
	}

	HANDLE hProcess;                                          /*���ڲ����ռ�õ��ڴ�*/
	PROCESS_MEMORY_COUNTERS pmc;
	hProcess = OpenProcess(PROCESS_QUERY_INFORMATION 
		| PROCESS_VM_READ,FALSE,GetCurrentProcessId());
	if (NULL == hProcess)
    {
        cout << "Process Hanle Error !" << endl;
		out_result.close();                                    /*��ֹ��Դй¶*/
		return -1;
    }

	double dP = 0.0;
	double timeCost = 0.0;
	__int64 start = 0;                                         /*���ڲ���ʱ��(��ȷ��1ms)*/
	__int64 frequency = 0;                                     /*�����ƽ̨���*/
	__int64 counter = 0;
	SIZE_T memsize;

	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency); 

	InputReader inReader(fileName,stFileName);	
	int s,t;
	Flow maxPmaxF;
	int maxflow;
	Graph g;
	g.Init();

	while(inReader.ReadGraph(g))                               /*��ȡ�ļ��е�ͼ*/
	{
		inReader.ReadSourceSink(s,t);
		/*��һ��ͼ���ݴ���һ��ͼ*/

		QueryPerformanceCounter((LARGE_INTEGER*)&start);      /*��¼��ʼʱ��*/
		dP = GetMPMF(g,s,t,maxflow,maxPmaxF,&StateMtrix);     /*���к����㷨*/
		//cout<< dP<< endl;

		/*ͨ��״̬�������ߵ����*/
		computeEdgeClass(&StateMtrix,&key_edge_set,g,s,t);

		QueryPerformanceCounter((LARGE_INTEGER*)&counter);    /*��¼����ʱ��*/
		timeCost = (counter - start) / double(frequency)*1000;/*���ص�λ�Ǻ���*/

		GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc));   /*��ý���ʹ�õ��ڴ�ʹ�����*/
		memsize = pmc.WorkingSetSize;                         /*��ý������ĵ��ڴ�*/

		/*�����н��������ļ�*/
		PrintFmax_Prob(out_result,s,t,maxflow,dP);            /*������ɿ�������ֲ����ʵ�����ļ�*/
		out_result<<"״̬�����㷨���ĵ��ڴ�Ϊ��"              /*����ڴ�ʹ�����������ļ�*/
			<<(double)memsize/MB<<endl;
		PrintTime(out_result,timeCost);                       /*��������ʱ�䵽����ļ�*/
		PrintFlow(out_result,maxPmaxF,g.nV);                  /*�����ɿ���������ֲ�������ļ�*/	
		/*���ؼ������*/
		printKeyEdge(out_result,key_edge_set);
		g.Init();
		/*��ʼ��*/
		init_KeyEdgeSet(key_edge_set,StateMtrix);
		
	}

	CloseHandle(hProcess);                                    /*�رս��̾��*/
	out_result.close();                                       /*�رս���ļ�*/
	//getchar();
	return 0;
}