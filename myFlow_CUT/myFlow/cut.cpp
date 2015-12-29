#include <set>
#include <vector>

/*关于图的一些算法*/
#include "graph.h"
using namespace std;



/*获取确定图中某个最大流对应割集中的边*/
/*参数说明：g:不确定图；source:顶点；sink:终点； MinCutEdges：割边集合的引用 */
int getMinCutEdges(Graph& g, int source, int sink, set<int>& MinCutEdges)
{
	/*求得剩余图,其中gf为剩余图*/
	Flow f;
	GF gf;
	int max_flow = Dinic(g,source,sink,f,gf);

	/*源测割点集*/
	set<int>S;
	/*使用剩余图获取源测割点集合*/
	MinCut(gf,g.nV, source, S);
	/*起点在S中，终点不在S中的边为关键边*/
	for (int i=1; i<=g.nE; i++)
	{
		if (S.count(g.AllEdge[i][1]) && !S.count(g.AllEdge[i][2]))
		{
			MinCutEdges.insert(g.AllEdge[i][4]);
		}
	}

	return 0;
}


