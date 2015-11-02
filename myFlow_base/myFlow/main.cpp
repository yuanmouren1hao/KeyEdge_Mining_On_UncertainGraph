// Network Flow.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "InputReader.h"
#include "graph.h"
#include "resultprint.h"
#include <iostream>
#include <string.h>
#include <windows.h>
#include  "psapi.h"                                           /*与进程相关(工程的属性中添加Psapi.lib)*/

/*用于使用下界子图的一些结构与方法*/
#include "state.h"

using namespace std; 

/*在release版本中能够成功生成*/
#pragma  once
#pragma  message("Psapi.h --> Linking with Psapi.lib")
#pragma  comment(lib,"Psapi.lib")


/*全局保存所有的下界子图状态*/
Lower_subGraph StateMtrix;
/*保存所有的关键边信息*/
KeyEdgeSet key_edge_set;
/*存储边的信息,不使用0*/
int AllEdge[MAX_E_NUM][5];/*1234，起点，终点，容量，有效位*/
double AllEdge_p[MAX_E_NUM];/**/


int _tmain(int argc, char* argv[])
{
	char stFileName[BUFFER_SIZE]     = WORK_SPACE;             /*源点汇点*/
	char fileName[BUFFER_SIZE]       = WORK_SPACE;             /*数据存放文件*/
	char resultFileName[BUFFER_SIZE] = WORK_SPACE;             /*实验结果存放文件*/

	
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
	out_result.open(resultFileName,ios::out|ios::app);         /*保存实验结果*/
	if(!out_result.is_open())
	{
		printf("open result file failed...\n");
		exit(1);
	}

	HANDLE hProcess;                                          /*用于测进程占用的内存*/
	PROCESS_MEMORY_COUNTERS pmc;
	hProcess = OpenProcess(PROCESS_QUERY_INFORMATION 
		| PROCESS_VM_READ,FALSE,GetCurrentProcessId());
	if (NULL == hProcess)
    {
        cout << "Process Hanle Error !" << endl;
		out_result.close();                                    /*防止资源泄露*/
		return -1;
    }

	double dP = 0.0;
	double timeCost = 0.0;
	__int64 start = 0;                                         /*用于测量时间(精确到1ms)*/
	__int64 frequency = 0;                                     /*与机器平台相关*/
	__int64 counter = 0;
	SIZE_T memsize;

	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency); 

	InputReader inReader(fileName,stFileName);	
	int s,t;
	Flow maxPmaxF;
	int maxflow;
	Graph g;
	g.Init();

	while(inReader.ReadGraph(g))                               /*读取文件中的图*/
	{
		inReader.ReadSourceSink(s,t);
		/*读一个图数据处理一个图*/

		QueryPerformanceCounter((LARGE_INTEGER*)&start);      /*记录开始时间*/
		dP = GetMPMF(g,s,t,maxflow,maxPmaxF,&StateMtrix);     /*运行核心算法*/
		//cout<< dP<< endl;

		/*通过状态矩阵计算边的类别*/
		computeEdgeClass(&StateMtrix,&key_edge_set,g,s,t);

		QueryPerformanceCounter((LARGE_INTEGER*)&counter);    /*记录结束时间*/
		timeCost = (counter - start) / double(frequency)*1000;/*返回单位是毫秒*/

		GetProcessMemoryInfo( hProcess, &pmc, sizeof(pmc));   /*获得进程使用的内存使用情况*/
		memsize = pmc.WorkingSetSize;                         /*获得进程消耗的内存*/

		/*将运行结果输出到文件*/
		PrintFmax_Prob(out_result,s,t,maxflow,dP);            /*保存最可靠最大流分布概率到结果文件*/
		out_result<<"状态划分算法消耗的内存为："              /*输出内存使用情况到结果文件*/
			<<(double)memsize/MB<<endl;
		PrintTime(out_result,timeCost);                       /*保存运行时间到结果文件*/
		PrintFlow(out_result,maxPmaxF,g.nV);                  /*输出最可靠的最大流分布到结果文件*/	
		/*将关键边输出*/
		printKeyEdge(out_result,key_edge_set);
		g.Init();
		/*初始化*/
		init_KeyEdgeSet(key_edge_set,StateMtrix);
		
	}

	CloseHandle(hProcess);                                    /*关闭进程句柄*/
	out_result.close();                                       /*关闭结果文件*/
	//getchar();
	return 0;
}