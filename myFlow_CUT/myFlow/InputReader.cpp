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
		return false;                                         /*�Ѿ�û��ͼ*/
	}
	g.nId = gId;

	in>>strLine;
	assert(strLine[0] == 's');                                /*��ȡ�������ͱ���*/
	in>>g.nV>>g.nE;

	assert(g.nV < MAX);                                       /*��ֹ�ڽӾ����С������*/

    /*�����ȡ�ߵ���Ϣ*/
	int u,v;                                                  /*u,v��һ���ߵ���������*/
	for(int i = 1; i <= g.nE; i++)
	{
		in>>strLine;
		assert(strLine[0] == 'e');                            
		in>>u>>v;                                             /*ע����������������д��һ��*/
		in>>g.matrix[u][v].iC>>g.matrix[u][v].dP>>g.matrix[u][v].iLabel;
		/*���ߵ���Ϣ�����ڶ�ά������*/
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