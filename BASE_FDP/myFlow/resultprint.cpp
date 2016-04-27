#include "stdafx.h"

#include "resultprint.h"
#include <iostream>
#include <iomanip>
#include <time.h>
#include <fstream>

using namespace std;


/*向控制台输出*/

void PrintGraph(Graph& g)
{
	cout<<g.nId<<endl;
	cout<<g.nV<<"   "<<g.nE<<endl;
	for(int i = 1; i <= g.nV; i++)
	{
		for(int j = 1; j <= g.nV; j++)
		{
			cout<<setw(4)<<g.matrix[i][j].iC<<"  ";
		}
		cout<<endl;
	}
	cout<<"---------------------------------------"<<endl;
}

void PrintFlow(ofstream &out, Flow& f, int n) 
{
	out<<"流分布如下: "<<endl;
	out<<"---------------------------------"<<endl;
	out<<setiosflags(ios::left);
	for(int i = 1; i <= n; i++)
	{
		for(int j = 1; j <= n; j++)
		{	
			out<<setw(4)<<(f[i][j] > 0 ? f[i][j] : 0);
		}
		out<<endl;
	}
	out<<"---------------------------------"<<endl;
}

void PrintMSN(MSN& msn)
{
	cout<<"分层网络中的顶点层次"<<endl;
	for(int i = 1; i <= msn.nV; i++)
	{
		cout<<"顶点"<<i<<"层次："<<msn.nVStage[i]<<endl;
	}
	cout<<"分层网络中的顶点层次"<<endl;
	for(int i = 1; i <= msn.nV; i++)
	{
		cout<<setiosflags(ios::left);
		for(int j = 1; j <= msn.nV; j++)
		{
			cout<<setw(4)<<msn.nMSN[i][j];
		}
		cout<<endl;
	}
	cout<<"--------------------------------------"<<endl;
}

void PrintX(int* x,int n)
{
	for(int i = 1; i <= n; i++)
	{
		cout<<x[i];
	}
	cout<<endl<<"---------------------------------"<<endl;
}

/*向文件中输出数据*/

void PrintN(const char * fileName,
					int nIterator,
					int nDinic,
					Graph& g,
					char* algName)
{
	ofstream fout(fileName, ios::out | ios::app);
	if(!fout)
	{
		cout<<"文件打开失败......."<<endl;
		return ;
	}
	fout<<"网络流图G("<<g.nV<<","<<g.nE<<") 使用---"
		<<algName<<"---迭代次数： "<<nIterator<<endl;
	fout<<"使用---"
		<<algName<<"---调用Dinic算法次数： "<<nDinic<<endl;
	fout.close();
}

void PrintCollection(Collection& c)
{
	ofstream fout("..\\printTemp\\collection.txt", ios::out | ios::app);
	if(!fout)
	{
		cout<<"文件打开失败......."<<endl;
		return ;
	}
	fout<<endl;
	fout<<"nE = "<<c.numE<<"\tI = "<<c.I<<"\tj = "<<c.j<<endl;
	for(int i = 1; i <= c.numE; i++)
	{
		fout<<c.lower[i]<<"  ";
	}
	fout<<endl;
	for(int i = 1; i <= c.numE; i++)
	{
		fout<<c.upper[i]<<"  ";
	}
	fout<<endl;
	for(int i = 1; i <= c.I; i++)
	{
		fout<<c.Ad_C[i]<<"  ";
	}
	fout<<endl<<"-------------------------------"<<endl;
	fout.close();
}


void PrintState(StateSet &statSet)
{
	ofstream fout("..\\printTemp\\temp1.txt", ios::out | ios::app);
	if(!fout)
	{
		cout<<"文件打开失败......."<<endl;
		return ;
	}
	for(int i = 1; i <= statSet.numE; i++)
	{
		fout<<statSet.lower[i]<<"  ";
	}
	fout<<endl;
	for(int i = 1; i <= statSet.numE; i++)
	{
		fout<<statSet.upper[i]<<"  ";
	}
	fout<<endl<<endl;
	fout.close();
}

void PrintFmax_Prob(ostream &out,int s,int t,int fmax,double dp)
{
	out<<s<<"--"<<t<<"最大流值为："<<fmax<<endl;
	out<<"最可靠最大流分布概率为 ： "<<dp<<endl;
}

void PrintTime(ostream &out,double cost)
{
	out<<"状态划分算法运行时间(ms) ： "<<cost<<endl;
	out<<"---------------------------------"<<endl;
}