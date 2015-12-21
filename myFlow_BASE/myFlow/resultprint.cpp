#include "stdafx.h"

#include "resultprint.h"
#include <iostream>
#include <iomanip>
#include <time.h>
#include <fstream>

using namespace std;


/*�����̨���*/

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
	out<<"���ֲ�����: "<<endl;
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
	cout<<"�ֲ������еĶ�����"<<endl;
	for(int i = 1; i <= msn.nV; i++)
	{
		cout<<"����"<<i<<"��Σ�"<<msn.nVStage[i]<<endl;
	}
	cout<<"�ֲ������еĶ�����"<<endl;
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

/*���ļ����������*/

void PrintN(const char * fileName,
					int nIterator,
					int nDinic,
					Graph& g,
					char* algName)
{
	ofstream fout(fileName, ios::out | ios::app);
	if(!fout)
	{
		cout<<"�ļ���ʧ��......."<<endl;
		return ;
	}
	fout<<"������ͼG("<<g.nV<<","<<g.nE<<") ʹ��---"
		<<algName<<"---���������� "<<nIterator<<endl;
	fout<<"ʹ��---"
		<<algName<<"---����Dinic�㷨������ "<<nDinic<<endl;
	fout.close();
}

void PrintCollection(Collection& c)
{
	ofstream fout("..\\printTemp\\collection.txt", ios::out | ios::app);
	if(!fout)
	{
		cout<<"�ļ���ʧ��......."<<endl;
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
		cout<<"�ļ���ʧ��......."<<endl;
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
	out<<s<<"--"<<t<<"�����ֵΪ��"<<fmax<<endl;
	out<<"��ɿ�������ֲ�����Ϊ �� "<<dp<<endl;
}

void PrintTime(ostream &out,double cost)
{
	out<<"״̬�����㷨����ʱ��(ms) �� "<<cost<<endl;
	out<<"---------------------------------"<<endl;
}