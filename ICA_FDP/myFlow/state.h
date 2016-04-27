#ifndef STATE_H
#define STATE_H

/*�洢�����½���ͼ��״̬*/
typedef struct {
	int State_Num;							/*�½���ͼ�ĸ���*/
	int Edge_Num;							/*�ߵĸ���*/
	int State[NUM_OF_SUB_GRAPH][MAX_E_NUM]; /*�洢�����½���ͼ��һ�б�ʾһ����ͼ��0��ʹ��*/
}Lower_subGraph;


/*����ؼ������ݽṹ*/
typedef struct{
	int Edge;								/*�ߵı��*/
	char Edge_Class;						/*�ߵ����ABC��*/
	double ChangeAmount_c;					/*�����ı仯��*/
	double ChangeAmount_p2;					/*���ֲ��ɿ��Ա仯��*/
}KeyEdge;
/*����ؼ��߼���*/
struct KeyEdgeSet{
	int EdgeNum;							/*�ߵĸ���*/
	int A_num;								/*A��ߵĸ���*/
	int B_num;								/*B��ߵĸ���*/
	int C_num;								/*C��ߵĸ���*/
	KeyEdge A_EdgeInfo[MAX_E_NUM];			/*�洢A�ߵ���Ϣ��������,0��ʹ��*/
	KeyEdge B_EdgeInfo[MAX_E_NUM];			/*�洢B�ߵ���Ϣ��������,0��ʹ��*/
	KeyEdge C_EdgeInfo[MAX_E_NUM];			/*�洢C�ߵ���Ϣ��������,0��ʹ��*/

	/*Ĭ�Ϲ��캯��*/
	KeyEdgeSet()
	{
		EdgeNum = 0;
		A_num = 0;
		B_num = 0;
		C_num = 0;

		for (int i=0; i<MAX_E_NUM; i++)
		{
			A_EdgeInfo[i].ChangeAmount_c = 0;
			A_EdgeInfo[i].ChangeAmount_p2 = 0;
			A_EdgeInfo[i].Edge = 0;
			B_EdgeInfo[i].ChangeAmount_c = 0;
			B_EdgeInfo[i].ChangeAmount_p2 = 0;
			B_EdgeInfo[i].Edge = 0;
			C_EdgeInfo[i].ChangeAmount_c = 0;
			C_EdgeInfo[i].ChangeAmount_p2 = 0;
			C_EdgeInfo[i].Edge = 0;
		}
	};
};



//ica�����㷨
void computeICA(Lower_subGraph *StateMtrix, KeyEdgeSet *key_edge_set, Graph& g, int source, int sink);


/*���ؼ������*/
void printKeyEdge(ostream &file_out,KeyEdgeSet &key_edge_set);
/*��ʼ���ؼ��߼�������ѭ������*/
void init_KeyEdgeSet(KeyEdgeSet &key_edge_set,Lower_subGraph &StateMtrix);


#endif