#include "stdafx.h"

#include "InputReader.h"
#include "graph.h"
#include <iostream>
#include <string>
#include "state.h"
#include <iomanip>
#include "partition.h"

#include <set>
#include <stack>
#include <vector>

using namespace std;

/*�Ƴ�A���е�ĳһ����*/
void removeA_Edge(Graph& g,int i, Edge &TempEdge)
{
	int u=g.AllEdge[i][1];
	int v=g.AllEdge[i][2];
	TempEdge.dP=g.matrix[u][v].dP;
	TempEdge.iC=g.matrix[u][v].iC;
	TempEdge.iLabel=g.matrix[u][v].iLabel;
	g.nE--;
	g.matrix[u][v].dP=0;
	g.matrix[u][v].iC=0;
	g.matrix[u][v].iLabel=0;
	return;
};

/*�ָ�ĳһ��A���*/
void restoreA_Edge(Graph& g,int i, Edge &TempEdge)
{
	int u=g.AllEdge[i][1];
	int v=g.AllEdge[i][2];
	g.nE++;
	g.matrix[u][v].dP=TempEdge.dP;
	g.matrix[u][v].iC=TempEdge.iC;
	g.matrix[u][v].iLabel=TempEdge.iLabel;
	return;
}


/*��ʼ����ʱ�洢��*/
void init_TempEdge(Edge &TempEdge)
{
	TempEdge.dP=0;
	TempEdge.iC=0;
	TempEdge.iLabel=0;
	return;
}

/*��A������¼���*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*���û��A��ߵĸ���Ϊ0*/
	if (key_edge_set->A_num==0)
	{
		return;
	}

	int new_maxflow = 0; /*��  �����*/
	double new_p1 =0;/*��  �ֲ��ɿ���*/
	double new_p2 =0;/*��  �����ɿ���*/
	Flow new_maxPmaxF;/*�µ�������ֲ�*/
	Edge TempEdge;/*��ʱ�洢һ���ߵ���Ϣ*/

	/*����A���еıߣ�ȥ��һ��*/
	for (int i = 1; i<=key_edge_set->A_num; i++)
	{
		/*��ʼ����ʱ�洢��*/
		init_TempEdge(TempEdge);
		/*ȥ��ĳһ��A���*/
		removeA_Edge(g, key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
		//���¼��������new_maxflow�����ֲ��ɿ���new_p1������ɿ���new_p2
		new_p1 = base_GetMPMF(g, source, sink, new_maxflow, new_maxPmaxF, new_p2, key_edge_set->A_EdgeInfo[i].Edge);
		
		/* ���ϵ���֮�����������ֲ��ɿ��ԣ������ɿ��Էֱ𱣴� */
		key_edge_set->A_EdgeInfo[i].ChangeAmount_c = new_maxflow;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p1 = new_p1;
		key_edge_set->A_EdgeInfo[i].ChangeAmount_p2 = new_p2;

		/* �ָ�ĳһ��A��� */
		restoreA_Edge(g,key_edge_set->A_EdgeInfo[i].Edge, TempEdge);
	}

	return;
}


/////////////////////////////////////////////////////////////////////////
//���й���B��ߵĲ���

/*������ͼ����*/
/*���룺vector<int>��ʾ��һ����ͼ,remvedEdge��ʾ��ɾ���ıߣ����Ϊ�����أ�����ͼ�ĸ���*/
double computeSubgraphProbabilityByVector(vector<int> a, int remvedEdge, Graph &g)
{
	double p = 1;
	for (int i=1; i<=g.nE; i++)
	{
		if (i!=remvedEdge)
		{
			if (a[i] == 0)
			{
				p *= 1-g.AllEdge_p[i];
			}else if (a[i] == 1)
			{
				p *= g.AllEdge_p[i];
			}
		}
	}
	return p;
}

/*����һЩ����ͼ�ĸ���֮��*/
/*���룺set<vector<int>>��ʾ���е���ͼ����,removedEdge��ʾ��Ҫ�Ƴ��ıߣ������������ͼ����֮��*/
double computeProbabilitySum(set<vector<int>> new_set, Graph &g, int removedEdge)
{
	double p = 0;

	set<vector<int>>::iterator it;
	for (it=new_set.begin(); it!=new_set.end(); it++)
	{
		p += computeSubgraphProbabilityByVector(*it, removedEdge, g);
	}
	return p;
}


/*�ж�һ����ͼ�������*/
/*���룺vector<int>��ʾ��һ����ͼ,removedEdge��ʾ���Ƴ��ıߣ����Ϊ0�Ļ���ô�ͱ�ʾû���Ƴ������أ�����ͼ�������*/
int DinicByVertor(vector<int> a, Graph& g, int removedEdge, int source, int sink)
{
	/*ͨ��vector����ͼ*/
	Graph new_g;
	new_g.nV = g.nV;
	if (0==removedEdge)
	{
		new_g.nE = g.nE;
	}else{
		new_g.nE = g.nE -1;
	}

	/*copy��֤��Ϣ*/
	int aa,bb,cc,dd;
	double ff;
	for (int i=1; i<=g.nE; i++)
	{
		if (i!=removedEdge)
		{
			if (a[i] == 0)
			{
				new_g.nE--;
			}
			else if (a[i] == 1)
			{
				aa = g.AllEdge[i][1];//��ʼ��
				bb = g.AllEdge[i][2];//�յ�
				cc = g.AllEdge[i][3];//����
				dd = g.AllEdge[i][4];//��Чλ���˴���ʹ�ã�

				ff = g.AllEdge_p[i];//�ɿ���

				new_g.matrix[aa][bb].iLabel = i;
				new_g.matrix[aa][bb].iC = cc;
				new_g.matrix[aa][bb].dP = ff;
			}
		}
	}

	/*ʹ��ԭ��Dinic�㷨���������*/
	Flow f;
	GF gf;
	return Dinic(new_g, source, sink, f, gf);
}



/*�Ƚ�����ƫ��ԵĴ�С*/
/*b������a����1��b>=a����2*/
int compareSubgraph(vector<int> a, vector<int> b, int length, int removeEdge)
{
	for (int i=1; i<=length; i++)
	{
		if (i!=removeEdge)
		{
			if (b[i]<a[i])
			{
				return 1;
			}
		}
	}
	return 2;
}

/*�������е�״̬�����ȡ���е������������״̬*/
set<vector<int>> getAllState(Lower_subGraph* StateMtrix)
{
	/*�������е�����ŵ�ջ��*/
	stack<vector<int>> s;

	for (int i=1; i<= StateMtrix->State_Num; i++)
	{
		/*����һ����ͼ����*/
		vector<int> v;
		v.push_back(0);/*Ϊ������0��Ԫ�أ��ӵ�һ��Ԫ�ؿ�ʼʹ��*/

		for (int j=1; j<=StateMtrix->Edge_Num; j++)
		{
			if (StateMtrix->State[i][j]==0 && StateMtrix->State_upper[i][j]==0)
			{
				v.push_back(0);
			}else if (StateMtrix->State[i][j]==1 && StateMtrix->State_upper[i][j]==1){
				v.push_back(1);
			}
			else if (StateMtrix->State[i][j]==0 && StateMtrix->State_upper[i][j]==1){
				v.push_back(2);
			}
		}
		s.push(v);
	}

	/*��������״̬*/
	set<vector<int>> allState ;
	/*��ջ�е�״̬���зֽ�*/
	while (!s.empty())
	{
		vector<int> a = s.top();
		s.pop();
		
		for (int i=1; i<=StateMtrix->Edge_Num; i++)
		{
			if (a[i]==2)
			{
				/*�������X��Ҳ����2������ͼ����ֽ⣬�ֱ�����ջ��*/
				vector<int> temp_V = a;
				temp_V[i] =0;
				s.push(temp_V);
				temp_V[i] =1;
				s.push(temp_V);
				break;
			}
			
			/*���ֱ��û�еĻ���ֱ�Ӳ��뵽 ������*/
			if (i==StateMtrix->Edge_Num)
			{
				allState.insert(a);		
			}
		}
	}
	return allState;
}




/*B����������㣬������ʹ�����õķ�ʽ*/
void sortKeyEdge_B(KeyEdgeSet * key_edge_set, Graph& g, int source, int sink, Lower_subGraph *StateMtrix)
{
	/*���û��B��ߵĸ���Ϊ0*/
	if (key_edge_set->B_num == 0)
	{
		return;
	}

	/*��С��ͼ*/
	/*ԭ��ȷ��ͼ�ļ�С��ͼ������StateMtrix�еĵ�0����ͼ�½���*/
	vector<int> minimumGraph;
	minimumGraph.push_back(0);/*�ӵ�һ��Ԫ�ؿ�ʼʹ��*/
	for (int i=1; i<=StateMtrix->Edge_Num; i++)
	{
		minimumGraph.push_back(StateMtrix->State[0][i]);
	}

	/*��ȡ���е������������״̬*/
	set<vector<int>> allState = getAllState(StateMtrix);

	for (int i =1; i<= key_edge_set->B_num; i++)
	{
		/*�ϵ�B������������*/
		key_edge_set->B_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		/*��ȡ������˵��½���ͼ*/
		vector<int> lowerGraph;

		/*�жϸñ��Ƿ��ڼ�С��ͼ��*/
		if (minimumGraph[key_edge_set->B_EdgeInfo[i].Edge]==0)
		{/*����ñ߲��ڼ�С��ͼ��*/
			
			/*�ֲ��ɿ��Բ���*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
			/*��Ҫ���˵�ͼ���Ǽ�С��ͼ*/
			lowerGraph = minimumGraph;

			set<vector<int>> new_set;
			set<vector<int>>::iterator it;
			for (it = allState.begin(); it!=allState.end(); it++)
			{
				vector<int> temp_V ;
				if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==2)
				{
					/*����ɸѡ��ͼ�����������*/
					temp_V = *it;
					temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
					new_set.insert(temp_V);
				}else if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==1)
				{
					/*�����ڱ�ʾ����ȷ���Ƿ����㣬��ʱ��Ҫ����*/
					if (DinicByVertor(*it, g, key_edge_set->B_EdgeInfo[i].Edge, source, sink)==g.max_flow)
					{
						/*�������ͼ�ܹ���������������棬���򣺲�����*/
						temp_V = *it;
						temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
						new_set.insert(temp_V);
					}
				}
			}
			
			/*������ɿ��ԣ������ɿ��ԣ�Ϊ��������������ĸ���֮��*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(new_set, g, key_edge_set->B_EdgeInfo[i].Edge);

		}else if(minimumGraph[key_edge_set->B_EdgeInfo[i].Edge]==1)
		{/*�ñ��ڼ�С��ͼ��*/

			/*��¼���ֲ��ɿ���*/
			double max_p1 =0;

			/*�ҵ������иñߵ������½���ͼ����Ϊ������ͼ*/
			for (int ii=1; ii<= StateMtrix->State_Num; ii++)
			{
				if (StateMtrix->State[ii][key_edge_set->B_EdgeInfo[i].Edge] == 0)
				{
					/*����ͼ����*/
					lowerGraph.push_back(0);
					for (int jj=1; jj<=StateMtrix->Edge_Num; jj++)
					{
						/*��Ϊ������ͼ*/
						lowerGraph.push_back(StateMtrix->State[ii][jj]);
					}
					break;
				}
			}

			/*��������ͼ��Ϊ��С��ͼ���������¼���*/
			max_p1 = computeSubgraphProbabilityByVector(lowerGraph, key_edge_set->B_EdgeInfo[i].Edge, g);

			set<vector<int>> new_set;
			set<vector<int>>::iterator it;
			for (it = allState.begin(); it!=allState.end(); it++)
			{
				vector<int> temp_V ;
				if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==2)
				{
					/*����ɸѡ��ͼ�����������*/
					temp_V = *it;
					temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
					new_set.insert(temp_V);

					/*�ж���Ϊ��С��ͼ*/
					if (computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g)>max_p1)
					{
						lowerGraph = temp_V;
						max_p1 = computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g);
					}

				}else if (compareSubgraph(lowerGraph, *it, StateMtrix->Edge_Num, key_edge_set->B_EdgeInfo[i].Edge)==1)
				{
					/*�����ڱ�ʾ����ȷ���Ƿ����㣬��ʱ��Ҫ����*/
					int new_maxflow = DinicByVertor(*it, g, key_edge_set->B_EdgeInfo[i].Edge, source, sink);
					if (new_maxflow==g.max_flow)
					{
						/*�������ͼ�ܹ���������������棬���򣺲�����*/
						temp_V = *it;
						temp_V[key_edge_set->B_EdgeInfo[i].Edge] =0;
						new_set.insert(temp_V);

						/*�ж��Ƿ�Ϊ��С��ͼ*/
						if (computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g)>max_p1)
						{
							lowerGraph = temp_V;
							max_p1 = computeSubgraphProbabilityByVector(temp_V, key_edge_set->B_EdgeInfo[i].Edge, g);
						}
					}
					
				}
			}

			/*�ֲ��ɿ���Ϊ���ϸ��µ�ֵ*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p1 = max_p1;
			/*������ɿ��ԣ������ɿ��ԣ�Ϊ��������������ĸ���֮��*/
			key_edge_set->B_EdgeInfo[i].ChangeAmount_p2 = computeProbabilitySum(new_set, g, key_edge_set->B_EdgeInfo[i].Edge);

		/*�ñ��ڼ�С��ͼ��*/
		}

	/*�������е�B���ѭ������*/
	}

	return;
}


//////////////////////////////////////////////////////////////////
//C��ߵĲ���

/*����C��߽��м���,C��߲���Ҫ����*/
void sortKeyEdge_C(KeyEdgeSet *key_edge_set,Graph &g)
{
	/*ÿһ��C��ߵķ�������֮��ģ�״̬���ı䣬ֱ�ӱ��漴��*/
	for (int i=1; i <= key_edge_set->C_num; i++)
	{
		key_edge_set->C_EdgeInfo[i].ChangeAmount_c = g.max_flow;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p1 = g.max_p1;
		key_edge_set->C_EdgeInfo[i].ChangeAmount_p2 = g.max_p2;
	}
	return;
}

/*�Ա߽��еĴ����*/
void sortKeyEdge(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set)
{
	//cout<<"StateMtrix->Edge_Num:"<<StateMtrix->Edge_Num<<endl;
	//cout<<"StateMtrix->State_Num:"<<StateMtrix->State_Num<<endl;
	/*��ʼ���ؼ��߼���*/
	key_edge_set->EdgeNum=StateMtrix->Edge_Num;
	int temp_sum=0;
	

	for (int ii=1;ii<=StateMtrix->Edge_Num;ii++)
	{
		temp_sum=0;
		/*ͨ����͵ķ�ʽ����ߵ����*/
		for (int jj=1;jj<=StateMtrix->State_Num;jj++)
		{
			temp_sum+=StateMtrix->State[jj][ii];
		}
		/*���ؼ��ߴ洢��ii��ʾ���Ǳߵı��*/
		if(temp_sum==StateMtrix->State_Num)
		{
			//cout<<ii<<" ��A��ߡ�"<<endl;
			key_edge_set->A_num++;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge=ii;
			key_edge_set->A_EdgeInfo[key_edge_set->A_num].Edge_Class='A';
		}
		else if (0==temp_sum)
		{
			//cout<<ii<<" ��C��ߡ�"<<endl;
			key_edge_set->C_num++;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge=ii;
			key_edge_set->C_EdgeInfo[key_edge_set->C_num].Edge_Class='C';
		}
		else if(0<temp_sum && temp_sum <StateMtrix->State_Num)
		{
			//cout<<ii<<" ��B���"<<endl;
			key_edge_set->B_num++;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge=ii;
			key_edge_set->B_EdgeInfo[key_edge_set->B_num].Edge_Class='B';
		}
	}
}

/*ͨ��״̬�������ߵ����*/
void computeEdgeClass(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink)
{
	/*���Է���*/
	sortKeyEdge(StateMtrix, key_edge_set);
	
	/*��������*/
	/*A������¼���*/
	sortKeyEdge_A(key_edge_set, g, source, sink);
	/*B�����������*/
	sortKeyEdge_B(key_edge_set, g, source, sink, StateMtrix);
	/*C��߲���Ҫ����*/
	sortKeyEdge_C(key_edge_set, g);
	
	return;
}














//////////////////////////////////////////////////////////////////
//�����ǽӿڣ����ڸ����е���


/*���ؼ������*/
void printKeyEdge(ostream &out,KeyEdgeSet &key_edge_set)
{
	out<<"�ؼ���������£�"<<endl<<"---------------------------------"<<endl;
	for (int i=1;i<=key_edge_set.A_num;i++)
	{
		out<<setw(5)<<key_edge_set.A_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.A_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.A_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.B_num;i++)
	{
		out<<setw(5)<<key_edge_set.B_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.B_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.B_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	for (int i=1;i<=key_edge_set.C_num;i++)
	{
		out<<setw(5)<<key_edge_set.C_EdgeInfo[i].Edge<<setw(10)<<key_edge_set.C_EdgeInfo[i].Edge_Class<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_c<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_p1<<setw(10)<<key_edge_set.C_EdgeInfo[i].ChangeAmount_p2<<endl;
	}
	out<<"---------------------------------"<<endl<<"---------------------------------"<<endl<<endl;
}


void init_KeyEdgeSet(KeyEdgeSet &key_edge_set,Lower_subGraph &StateMtrix)
{
	/*��ʼ���ؼ��߼���*/
	key_edge_set.A_num=0;
	key_edge_set.B_num=0;
	key_edge_set.C_num=0;
	key_edge_set.EdgeNum=0;
	for (int i=0;i<MAX_E_NUM;i++)
	{
		key_edge_set.A_EdgeInfo[i].Edge=0;
		key_edge_set.B_EdgeInfo[i].Edge=0;
		key_edge_set.C_EdgeInfo[i].Edge=0;

		key_edge_set.A_EdgeInfo[i].Edge_Class=0;
		key_edge_set.B_EdgeInfo[i].Edge_Class=0;
		key_edge_set.C_EdgeInfo[i].Edge_Class=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount_c=0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_c=0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount_p1 = 0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_p1 = 0;

		key_edge_set.A_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.B_EdgeInfo[i].ChangeAmount_p2 = 0;
		key_edge_set.C_EdgeInfo[i].ChangeAmount_p2 = 0;
	}
	/*��ʼ���½���ͼ*/
	StateMtrix.Edge_Num=0;
	StateMtrix.State_Num=0;
	for (int i=0;i<NUM_OF_SUB_GRAPH;i++)
	{
		for (int j=0;j<MAX_E_NUM;j++)
		{
			StateMtrix.State[i][j]=0;
		}
	}
	return;
}