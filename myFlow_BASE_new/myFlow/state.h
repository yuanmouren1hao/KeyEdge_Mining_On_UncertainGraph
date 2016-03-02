#ifndef STATE_H
#define STATE_H

/*存储所有下界子图的状态*/
typedef struct {
	int State_Num;							/*下界子图的个数*/
	int Edge_Num;							/*边的个数*/
	int State[NUM_OF_SUB_GRAPH][MAX_E_NUM]; /*存储所有下界子图，一行表示一个子图，0不使用*/
	int State_upper[NUM_OF_SUB_GRAPH][MAX_E_NUM];/*存储所有子图上界，一行表示一个子图，0不使用*/
}Lower_subGraph;


/*定义关键边数据结构*/
typedef struct{
	int Edge;								/*边的编号*/
	char Edge_Class;						/*边的类别，ABC类*/
	double ChangeAmount_c;					/*流量的变化量*/
	double ChangeAmount_p1;					/*随机流网络可靠性变化量*/
	double ChangeAmount_p2;					/*流分布可靠性变化量*/
}KeyEdge;
/*定义关键边集合*/
struct KeyEdgeSet{
	int EdgeNum;							/*边的个数*/
	KeyEdge EdgeInfo[MAX_E_NUM];			/*存储边的信息，0不使用*/

	/*默认构造函数*/
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


/*将关键边输出*/
void printKeyEdge(ostream &file_out,KeyEdgeSet &key_edge_set);
/*初始化关键边集，用于循环计算*/
void init_KeyEdgeSet(KeyEdgeSet &key_edge_set);

//base算法，核心计算
void computeBase(KeyEdgeSet *key_edge_set, Graph& g, int source, int sink);

#endif