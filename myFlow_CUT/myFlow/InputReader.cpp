#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
using namespace std;


InputReader::InputReader(const char* fileName, const char* stFileName)
{
	in.open(fileName);
	if(!in.is_open())
	{
		cout<<"Input file "<<fileName<<" doesn't exist!"<<endl;
		exit(1);
	}

	stIn.open(stFileName);
	if(!stIn.is_open())
	{
		cout<<"Input file "<<stFileName<<" doesn't exist!"<<endl;
		in.close();
		exit(1);
	}
}

InputReader::~InputReader()
{
	in.close();
	stIn.close();
}

void InputReader::ReadFirstLine()
{
	in>>strLine;
	assert(strLine[0] == 'g');
	
	in>>strLine;
	assert(strLine[0] == '#');

	in>> gId;
}


bool InputReader::ReadGraph(Graph &g)
{
	ReadFirstLine();
	if(gId == 0)
	{
		return false;                                         /*已经没有图*/
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
		/*将边的信息保存在二维数组中*/
		g.AllEdge[i][1]=u;
		g.AllEdge[i][2]=v;
		g.AllEdge[i][3]=g.matrix[u][v].iC;
		g.AllEdge[i][4]=g.matrix[u][v].iLabel;
		g.AllEdge_p[i]=g.matrix[u][v].dP;
	}

	return true;
}

void InputReader::ReadSourceSink(int &s, int &t)
{
	stIn>>s>>t;
}