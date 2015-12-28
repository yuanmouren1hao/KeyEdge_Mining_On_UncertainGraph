#ifndef INPUTREADER_H
#define INPUTREADER_H

#include "ConstDef.h"
#include <fstream>
using namespace std;

class Graph;

class InputReader
{
public:
	InputReader(const char* fileName, const char* stFileName);
	~InputReader();

	bool ReadGraph(Graph& g);
	void ReadSourceSink(int &s, int &t);

private:
	void ReadFirstLine();  /*读取（g # 1）*/

private:
	ifstream in;          /*读graph*/
	ifstream stIn;        /*源点s汇点t*/
	char strLine[32];     /*记录in读取的内容*/
	int gId;
};

#endif