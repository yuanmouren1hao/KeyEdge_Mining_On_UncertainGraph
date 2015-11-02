#ifndef STATE_H
#define STATE_H

/*存储所有下界子图的状态*/
typedef struct {
	int State_Num;							/*下界子图的个数*/
	int Edge_Num;							/*边的个数*/
	int State[NUM_OF_SUB_GRAPH][MAX_E_NUM]; /*存储所有下界子图，一行表示一个子图，0不使用*/
}Lower_subGraph;


/*定义关键边数据结构*/
typedef struct{
	int Edge;								/*边的编号*/
	char Edge_Class;						/*边的类别，ABC类*/
	double ChangeAmount;					/*改变量，如果是A类边则表示流量的变化，B为可靠性的变化*/								
}KeyEdge;
/*定义关键边集合*/
typedef struct{
	int EdgeNum;							/*边的个数*/
	int A_num;								/*A类边的个数*/
	int B_num;								/*B类边的个数*/
	int C_num;								/*C类边的个数*/
	KeyEdge A_EdgeInfo[MAX_E_NUM];			/*存储A边的信息，有排序,0不使用*/
	KeyEdge B_EdgeInfo[MAX_E_NUM];			/*存储B边的信息，有排序,0不使用*/
	KeyEdge C_EdgeInfo[MAX_E_NUM];			/*存储C边的信息，有排序,0不使用*/
}KeyEdgeSet;



/*通过状态矩阵计算边的类别*/
void computeEdgeClass(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set,Graph& g,int source,int sink,int AllEdge[][5],double AllEdge_p[]);
/*对边进行的大分类*/
void sortKeyEdge(Lower_subGraph *StateMtrix,KeyEdgeSet *key_edge_set);
/*对于A类边进行排序*/
void sortKeyEdge_A(KeyEdgeSet *key_edge_set,Graph& g,int source,int sink,int AllEdge[][5]);
/*对于B类边进行排序*/
void sortKeyEdge_B(KeyEdgeSet *key_edge_set,double AllEdge_p[],Lower_subGraph *StateMtrix);
/*将关键边输出*/
void printKeyEdge(ostream &file_out,KeyEdgeSet &key_edge_set);
/*初始化关键边集，用于循环计算*/
void init_KeyEdgeSet(KeyEdgeSet &key_edge_set,Lower_subGraph &StateMtrix);


#endif