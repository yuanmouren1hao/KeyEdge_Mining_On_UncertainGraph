#ifndef STATE_H
#define STATE_H

/*�洢�����½���ͼ��״̬*/
typedef struct {
	int State_Num;							/*�½���ͼ�ĸ���*/
	int Edge_Num;							/*�ߵĸ���*/
	int State[NUM_OF_SUB_GRAPH][MAX_E_NUM]; /*�洢�����½���ͼ��һ�б�ʾһ����ͼ��0��ʹ��*/
	int State_upper[NUM_OF_SUB_GRAPH][MAX_E_NUM];/*�洢������ͼ�Ͻ磬һ�б�ʾһ����ͼ��0��ʹ��*/
}Lower_subGraph;


/*����ؼ������ݽṹ*/
typedef struct{
	int Edge;								/*�ߵı��*/
	char Edge_Class;						/*�ߵ����ABC��*/
	double ChangeAmount_c;					/*�����ı仯��*/
	double ChangeAmount_p1;					/*���������ɿ��Ա仯��*/
	double ChangeAmount_p2;					/*���ֲ��ɿ��Ա仯��*/
}KeyEdge;
/*����ؼ��߼���*/
struct KeyEdgeSet{
	int EdgeNum;							/*�ߵĸ���*/
	KeyEdge EdgeInfo[MAX_E_NUM];			/*�洢�ߵļ�������Ϣ,0��ʹ��*/

	/*Ĭ�Ϲ��캯��*/
	KeyEdgeSet()
	{
		EdgeNum = 0;
		for (int i=0; i<MAX_E_NUM; i++)
		{
			EdgeInfo[i].ChangeAmount_c = 0;
			EdgeInfo[i].ChangeAmount_p1 = 0;
			EdgeInfo[i].ChangeAmount_p2 = 0;
			EdgeInfo[i].Edge = 0;
		}
	};

};



/*ͨ��״̬�������ߵ����*/
void computeEdgeClass(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set,Graph& g,int source,int sink);
/*�Ա߽��еĴ����*/
void sortKeyEdge(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set);
/*����A��߽�������*/
void sortKeyEdge_ABC(KeyEdgeSet *key_edge_set,Graph& g,int source,int sink);
/*���ؼ������*/
void printKeyEdge(ostream &file_out,KeyEdgeSet &key_edge_set);
/*��ʼ���ؼ��߼�������ѭ������*/
void init_KeyEdgeSet(KeyEdgeSet &key_edge_set);

//base_all���㷨����������
void computeBaseAll(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink);


#endif