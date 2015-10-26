#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include "graph.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include <time.h>
#include <windows.h>
#include <iomanip> 

using namespace std;

/* Change any of these parameters to match your needs */ 
#define POPSIZE 50 /* population size */ 
#define MAXGENS 1000 /* max. number of generations */ 
//#define NVARS 26 /* no. of problem variables TODO*/ 
const int NVARS = 64 ;
int source;
int dest;
string FileName = "data/data2/V16E64.txt";
string stFileName ="data/st/V6E10.txt";
char * resultFileName = "result/data2/V16E64.txt";

#define PXOVER 0.95 /* probability of crossover */ 
#define PMUTATION 0.02 /* probability of mutation */ 
#define TRUE 1 
#define FALSE 0 

int generation; /* current generation no. */ 
int cur_best; /* best individual */ 
FILE *galog; /* an output file */ 

/*uncertern graph*/
Graph g;
Flow flow;
GF gf;
int AllEdge[MAX_E_NUM][5];/*1234，起点，终点，容量，有效位*/
double AllEdge_p[MAX_E_NUM];/*存储边概率*/


struct genotype /* genotype (GT), a member of the population */ 
{ 
	short int gene[NVARS]; /* a string of variables */
	int flow; /*the flow that subgraph can be*/
	int is_maxFlow; /*1 for yes, 0 for no, 2 for havenot calculate yet. */
	double fitness; /* GT's fitness */ 
	//double upper[NVARS]; /* GT's variables upper bound */ 
	//double lower[NVARS]; /* GT's variables lower bound */ 
	double rfitness; /* relative fitness */ 
	double cfitness; /* cumulative fitness */ 
}; 

struct genotype population[POPSIZE+1]; /* population */ 
struct genotype newpopulation[POPSIZE+1]; /* new population; */ 

/* Declaration of procedures used by this genetic algorithm */ 

void initialize(void); 
double randval(double, double); 
void evaluate(void); 
void keep_the_best(void); 
void elitist(void); 
void select(void); 
void crossover(void); 
void Xover(int,int); 
void swap(short int *, short int *); 
void mutate(void); 
void report(void); 
void loadgraph(void);
double calculateSubgraph(int, int &, int &);

/***************************************************************/ 
/* Initialization function: Initializes the values of genes */ 
/* within the variables bounds. It also initializes (to zero) */ 
/* all fitness values for each member of the population. It */ 
/* reads upper and lower bounds of each variable from the */ 
/* input file `gadata.txt'. It randomly generates values */ 
/* between these bounds for each gene of each genotype in the */ 
/* population. The format of the input file `gadata.txt' is */ 
/* var1_lower_bound var1_upper bound */ 
/* var2_lower_bound var2_upper bound ... */ 
/***************************************************************/ 
void initialize(void) 
{ 
	int i, j; 

	/*使用一个满足最大流的个体初始化最初种群*/
	population[0].cfitness = 0;
	population[0].fitness = 0;
	population[0].flow = g.maxFlow;
	population[0].is_maxFlow = 1;
	population[0].rfitness = 0;
	int u, v;/*起始顶点*/
	for (int i = 1; i <= NVARS; i++)
	{
		u = AllEdge[i][1];
		v = AllEdge[i][2];
		if (flow[u][v] > 0)
			population[0].gene[i-1] = 1;
		else
			population[0].gene[i-1] = 0;
	}


	for (j = 1; j < POPSIZE; j++) 
	{ 
		population[j].fitness = 0; 
		population[j].rfitness = 0; 
		population[j].cfitness = 0; 
		population[j].is_maxFlow = 0;
		population[j].flow = 0;
		for (i = 0; i < NVARS; i++){
			population[j].gene[i] = rand()%2;
		}
	} 

	return;
} 
 

/*************************************************************/ 
/* Evaluation function: This takes a user defined function. */ 
/* Each time this is changed, the code has to be recompiled. */ 
/* The current function is: x[1]^2-x[1]*x[2]+x[3] */ 
/*************************************************************/ 
void evaluate(void) 
{ 
	int mem; 

	for (mem = 0; mem < POPSIZE; mem++) 
	{ 
		/*TODO 这里需要计算子图的可靠性*/
		population[mem].fitness = calculateSubgraph(mem, population[mem].is_maxFlow, population[mem].flow);
	} 
	return;
} 


/************************************************************************/
/*计算子图的 最大流 和 可靠性                                           */
/************************************************************************/
double calculateSubgraph(int mem, int & is_maxFlow, int & flow){
	int i;
	Graph sub_graph;
	int u,v;
	Flow sub_flow;
	GF sub_gf;
	double p = 1.0;

	sub_graph.nE = g.nE;
	sub_graph.nV = g.nV;
	sub_graph.nId = g.nId;

	for (i=0; i<NVARS; ++i)
	{
		if ( population[mem].gene[i] == 1 )
		{
			u = AllEdge[i+1][1];
			v = AllEdge[i+1][2];
			sub_graph.matrix[u][v].iC = AllEdge[i+1][3];
		}
	}
	/*使用dinic算法 计算子图的最大流*/
	flow = Dinic(sub_graph, source, dest, sub_flow, sub_gf);

	/*如果不能满足最大流不需要计算，直接赋值fitness为0。否则计算fitness=各边之积*/
	if (flow == g.maxFlow ){
		is_maxFlow = 1;
		for (i=0; i<NVARS; ++i)
		{
			if ( population[mem].gene[i] == 1 ){
				p *= AllEdge_p[i+1];
			}
		}
		return p;
	}else{
		is_maxFlow = 0;
		return 0.0;
	}
	
}



/***************************************************************/ 
/* Keep_the_best function: This function keeps track of the */ 
/* best member of the population. Note that the last entry in */ 
/* the array Population holds a copy of the best individual */ 
/***************************************************************/ 
/*找到最优的一个  个体（能够满足最大流，且可靠性最大） */
void keep_the_best() 
{ 
	int mem; 
	int i; 
	cur_best = 0; /* stores the index of the best individual */ 

	for (mem = 0; mem < POPSIZE; mem++) 
	{ 
		/*找到能够达到最大流，且可靠性最大的子图*/
		if (population[mem].fitness > population[POPSIZE].fitness) 
		{ 
			cur_best = mem; 
			population[POPSIZE].fitness = population[mem].fitness; 
		} 
	} 
	/* once the best member in the population is found, copy the genes */ 
	for (i = 0; i < NVARS; i++) 
		population[POPSIZE].gene[i] = population[cur_best].gene[i];
	population[POPSIZE].is_maxFlow = population[cur_best].is_maxFlow;
	population[POPSIZE].flow = population[cur_best].flow;
	population[POPSIZE].fitness = population[cur_best].fitness;
	population[POPSIZE].rfitness = population[cur_best].rfitness;
	population[POPSIZE].cfitness = population[cur_best].cfitness;
} 

/****************************************************************/ 
/* Elitist function: The best member of the previous generation */ 
/* is stored as the last in the array. If the best member of */ 
/* the current generation is worse then the best member of the */ 
/* previous generation, the latter one would replace the worst */ 
/* member of the current population */ 
/****************************************************************/ 
void elitist() 
{ 
	int i; 
	double best, worst; /* best and worst fitness values */ 
	int best_mem, worst_mem; /* indexes of the best and worst member */ 

	best = population[0].fitness; 
	worst = population[0].fitness;
	best_mem = 0;
	worst_mem = 0;

	/*找出最好和最坏的个体*/
	for (i = 0; i < POPSIZE; ++i) 
	{
		if (population[i].fitness > best)
		{
			best = population[i].fitness;
			best_mem = i;
		}

		if (population[i].fitness < worst)
		{
			worst = population[i].fitness;
			worst_mem = i;
		}
	} 
	/* if best individual from the new population is better than */ 
	/* the best individual from the previous population, then */ 
	/* copy the best from the new population; else replace the */ 
	/* worst individual from the current population with the */ 
	/* best one from the previous generation */ 
	/*判断当前一代最优种群是否   好过于  上一代的种群的最优*/
	if (best > population[POPSIZE].fitness) 
	{ 
		for (i = 0; i < NVARS; i++) 
			population[POPSIZE].gene[i] = population[best_mem].gene[i]; 
		population[POPSIZE].fitness = population[best_mem].fitness;
		population[POPSIZE].flow = population[best_mem].flow;
		population[POPSIZE].is_maxFlow = population[best_mem].is_maxFlow;
	} 
	else 
	{ 
		for (i = 0; i < NVARS; i++) 
			population[worst_mem].gene[i] = population[POPSIZE].gene[i]; 
		population[worst_mem].fitness = population[POPSIZE].fitness; 
		population[worst_mem].flow = population[POPSIZE].flow;
		population[worst_mem].is_maxFlow = population[POPSIZE].is_maxFlow;
	} 
} 


/**************************************************************/ 
/* Selection function: Standard proportional selection for */ 
/* maximization problems incorporating elitist model - makes */ 
/* sure that the best member survives */ 
/**************************************************************/ 
/*选择遗传操作*/
void select(void) 
{ 
	int mem, i, j; 
	double sum = 0; 
	double p; 

	/* find total fitness of the population */ 
	for (mem = 0; mem < POPSIZE; mem++) 
	{ 
		sum += population[mem].fitness; 
	} 

	/* calculate relative fitness */ 
	for (mem = 0; mem < POPSIZE; mem++) 
	{ 
		population[mem].rfitness = population[mem].fitness/sum; 
	} 
	population[0].cfitness = population[0].rfitness; 

	/* calculate cumulative fitness */ 
	for (mem = 1; mem < POPSIZE; mem++) 
	{ 
		population[mem].cfitness = population[mem-1].cfitness + 
			population[mem].rfitness; 
	} 

	/* finally select survivors using cumulative fitness. */ 

	for (i = 0; i < POPSIZE; i++) 
	{ 
		p = rand()%1000/1000.0; 
		if (p < population[0].cfitness) 
			newpopulation[i] = population[0]; 
		else 
		{ 
			for (j = 0; j < POPSIZE;j++) 
				if (p >= population[j].cfitness && 
					p<population[j+1].cfitness) 
					newpopulation[i] = population[j+1]; 
		} 
	} 
	/* once a new population is created, copy it back */ 

	for (i = 0; i < POPSIZE; i++) 
		population[i] = newpopulation[i]; 
} 

/***************************************************************/ 
/* Crossover selection: selects two parents that take part in */ 
/* the crossover. Implements a single point crossover */ 
/***************************************************************/ 
void crossover(void) 
 { 
	int mem, one; 
	int first = 0; /* count of the number of members chosen */ 
	double x; 

	for (mem = 0; mem < POPSIZE; ++mem) 
	{ 
		x = rand()%1000/1000.0;
		/*满足交叉率*/
		if (x < PXOVER) 
		{ 
			++first; 
			if (first % 2 == 0) 
				Xover(one, mem); 
			else 
				one = mem; 
		} 
	} 
} 
/**************************************************************/ 
/* Crossover: performs crossover of the two selected parents. */ 
/**************************************************************/ 
/*执行交差操作*/
void Xover(int one, int two) 
{ 
	int i; 
	int point; /* crossover point */ 

	/* select crossover point */ 
	/*选择交差位置*/
	point = (rand() % (NVARS - 1)) + 1; 

	for (i = 0; i < point; i++) 
		swap(&population[one].gene[i], &population[two].gene[i]); 
} 

/*************************************************************/ 
/* Swap: A swap procedure that helps in swapping 2 variables */ 
/*************************************************************/ 
void swap(short int *x,  short int *y) 
{ 
	int temp; 

	temp = *x; 
	*x = *y; 
	*y = temp; 

} 

/**************************************************************/ 
/* Mutation: Random uniform mutation. A variable selected for */ 
/* mutation is replaced by a random value between lower and */ 
/* upper bounds of this variable */ 
/**************************************************************/ 
/*执行变异操作*/
void mutate(void) 
{ 
	int i, j; 
	//double lbound, hbound; 
	double x; 

	for (i = 0; i < POPSIZE; i++) 
		for (j = 0; j < NVARS; j++) 
		{ 
			x = rand()%1000/1000.0;
			/*如果满足变异率*/
			if (x < PMUTATION) 
			{ 
				/* find the bounds on the variable to be mutated */ 
				//lbound = population[i].lower[j]; 
				//hbound = population[i].upper[j]; 
				/*如果是1变为0.如果是0变为1*/
				population[i].gene[j] = (population[i].gene[j]+1)%2; 
			} 
		} 
} 

/***************************************************************/ 
/* Report function: Reports progress of the simulation. Data */ 
/* dumped into the output file are separated by commas */ 
/***************************************************************/ 
/*生成种群的报告*/
void report(void) 
{ 
	int i; 
	double best_val; /* best population fitness */ 
	double avg; /* avg population fitness */ 
	double stddev; /* std. deviation of population fitness */ 
	double sum_square; /* sum of square for std. calc */ 
	double square_sum; /* square of sum for std. calc */ 
	double sum; /* total population fitness */ 

	sum = 0.0; 
	sum_square = 0.0; 

	for (i = 0; i < POPSIZE; i++) 
	{ 
		sum += population[i].fitness; 
		sum_square += population[i].fitness * population[i].fitness; 
	} 

	avg = sum/(double)POPSIZE; 
	square_sum = avg * avg * POPSIZE; 
	stddev = sqrt((sum_square - square_sum)/(POPSIZE - 1)); 
	best_val = population[POPSIZE].fitness; 

	fprintf(galog, "\n%5d, %f, %f, %f ", generation, 
		best_val, avg, stddev); 
} 

/************************************************************************/
/* loadgraph                                                                     */
/************************************************************************/
void loadgraph(void){
	ifstream in;          /*读graph*/
	ifstream stIn;        /*源点s汇点t*/
	char strLine[32];     /*记录in读取的内容*/
	int gId;	

	in.open(FileName);
	if(!in.is_open())
	{
		cout<<"Input file graph.txt doesn't exist!"<<endl;
		exit(1);
	}

	stIn.open(stFileName);
	if(!stIn.is_open())
	{
		cout<<"Input file st.txt doesn't exist!"<<endl;
		in.close();
		exit(1);
	}
	
	/*read st file name*/
	stIn >> source >> dest;

	/*read first line and edge and v*/
	in>>strLine;
	assert(strLine[0] == 'g');
	in>>strLine;
	assert(strLine[0] == '#');
	in>> gId;
	if(gId == 0)
	{
		return ;                                         /*已经没有图*/
	}
	g.nId = gId;
	in>>strLine;
	assert(strLine[0] == 's');                                /*读取顶点数和边数*/
	in>>g.nV>>g.nE;

	assert(g.nV < MAX);                                       /*防止邻接矩阵大小不够用*/

	/*下面读取边的信息*/
	int u,v;                                                  /*u,v是一条边的两个顶点*/
	for(int i = 1; i <= g.nE; i++)
	{
		in>>strLine;
		assert(strLine[0] == 'e');                            
		in>>u>>v;                                             /*注意这个不能与下面的写到一起*/
		in>>g.matrix[u][v].iC>>g.matrix[u][v].dP>>g.matrix[u][v].iLabel;
		/*将边的信息保存在二维数组中，方便查询*/
		AllEdge[i][1]=u;
		AllEdge[i][2]=v;
		AllEdge[i][3]=g.matrix[u][v].iC;
		AllEdge[i][4]=0;
		AllEdge_p[i]=g.matrix[u][v].dP;
	}

	return ;
}

/**************************************************************/ 
/* Main function: Each generation involves selecting the best */ 
/* members, performing crossover & mutation and then */ 
/* evaluating the resulting population, until the terminating */ 
/* condition is satisfied */ 
/**************************************************************/ 
void main2(void) 
{ 
	int i; 
	__int64 start = 0;                                         /*用于测量时间(精确到1ms)*/
	__int64 frequency = 0;                                     /*与机器平台相关*/
	__int64 counter = 0;
	double timeCost = 0.0;

	if ((galog = fopen(resultFileName, "w"))==NULL) 
	{ 
		exit(1); 
	} 
	generation = 0; 

	fprintf(galog, "\n generation best average standard \n"); 
	fprintf(galog, " number value fitness deviation \n"); 

	loadgraph();
	/*获取最大流*/
	g.maxFlow = Dinic(g, source, dest, flow, gf);
	fprintf(galog, " maxFlow:%d \n", g.maxFlow); 

	/*
	for(int i =1; i<= NVARS; i++){
		for (int j = 1; j<=NVARS; j++)
		{
			cout<<setw(4)<<flow[i][j];
		}
		cout<<endl;
	}*/

	initialize(); 
	evaluate(); 
	keep_the_best(); 


	/*随机数种子*/
	srand(time(0));
	QueryPerformanceFrequency((LARGE_INTEGER*)&frequency); 
	/*开始时间*/
	QueryPerformanceCounter((LARGE_INTEGER*)&start);      /*记录开始时间*/
	while(generation<MAXGENS) 
	{ 
		generation++; 
		select(); 
		crossover(); 
		mutate(); 
		report(); 
		evaluate(); 
		elitist(); 
	}

	/*结束时间*/
	QueryPerformanceCounter((LARGE_INTEGER*)&counter);      /*记录结束时间*/
	timeCost = (counter - start)*1000 / double(frequency);/*返回单位是毫秒*/
	fprintf(galog, " \nExecute time is:%f ms\n", timeCost); 

	fprintf(galog,"\n\n Simulation completed\n"); 
	fprintf(galog,"\n Best member: \n"); 

	for (i = 0; i < NVARS; i++) 
	{ 
		fprintf (galog,"\n var(%d) = %d",i,population[POPSIZE].gene[i]); 
	} 
	fprintf(galog,"\n\n Best fitness = %f",population[POPSIZE].fitness); 
	fclose(galog); 
	printf("Success\n");
} 
/***************************************************************/