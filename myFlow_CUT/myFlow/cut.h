#ifndef CUT_H
#define CUT_H

#include <set>

/*拆分子图SplitGraph*/
typedef struct{

}SplitGraph; 




/*获取确定图中某个最大流对应割集中的边*/
/*参数说明：g:不确定图；source:顶点；sink:终点； MinCutEdges：割边集合的引用 */
int getMinCutEdges(Graph& g, int source, int sink, set<int>& MinCutEdges);

#endif