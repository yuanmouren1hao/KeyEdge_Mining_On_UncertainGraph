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
	void ReadFirstLine();  /*��ȡ��g # 1��*/

private:
	ifstream in;          /*��graph*/
	ifstream stIn;        /*Դ��s���t*/
	char strLine[32];     /*��¼in��ȡ������*/
	int gId;
};

#endif