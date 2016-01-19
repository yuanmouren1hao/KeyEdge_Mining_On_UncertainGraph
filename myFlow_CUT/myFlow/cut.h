#ifndef CUT_H
#define CUT_H

#include <set>

/*拆分子图SplitGraph*/
typedef struct{

}SplitGraph; 


/*确定边的集合，包括割集中的边和悬挂边*/
typedef struct{
	set<int> MinCutEdges;/*割集中的边*/
	set<int> HangEdges;/*悬挂边*/
}CertainEdge;


/*获取确定的边，包括悬挂边和割集中的边*/
/*参数说明：g:不确定图；source:顶点；sink:终点； certainEdge：确定边集合的引用 */
void getCertainEdges(Graph& g, int source, int sink, CertainEdge& certainEdge);

#endif