#include <set>
#include <vector>

/*关于图的一些算法*/
#include "graph.h"
#include "cut.h"
using namespace std;



/*获取确定图中某个最大流对应割集中的边*/
/*参数说明：g:不确定图；source:顶点；sink:终点； certainEdge：确定边集合的引用 */
int getMinCutEdges(Graph& g, int source, int sink,  CertainEdge& certainEdge)
{
	/*求得剩余图,其中gf为剩余图*/
	Flow f;
	GF gf;
	int max_flow = Dinic(g,source,sink,f,gf);

	/*源测割点集*/
	set<int> S;
	/*使用剩余图获取源测割点集合*/
	MinCut(gf,g.nV, source, S);
	/*起点在S中，终点不在S中的边为关键边*/
	for (int i=1; i<=g.nE; i++)
	{
		if (S.count(g.AllEdge[i][1]) && !S.count(g.AllEdge[i][2]))
		{
			certainEdge.MinCutEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_MinCutEdges++;
		}
	}

	return 0;
}

/*获取确定的边，悬挂边*/
/*参数说明：g:不确定图；source:顶点；sink:终点； certainEdge：确定边集合的引用 */
void getHangEdge(Graph& g, int source, int sink, CertainEdge& certainEdge)
{
	set<int> in_set,out_set;/*起点in,终点out*/
	for(int i=1; i<=g.nE; i++)
	{
		in_set.insert(g.AllEdge[i][1]);
		out_set.insert(g.AllEdge[i][2]);
	}

	for(int i=1; i<=g.nE; i++)
	{
		/*源点流入的边是悬挂边*/
		if (g.AllEdge[i][2] == source)
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}
		/*汇点流出的边是悬挂边*/
		if (g.AllEdge[i][1] == sink)
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}

		/*只有流入，没有流出的点*/
		if (in_set.count(g.AllEdge[i][2])==0 && out_set.count(g.AllEdge[i][2])>0 && sink != g.AllEdge[i][2])
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}

		/*只有流出，没有流入的点*/
		if (in_set.count(g.AllEdge[i][1])>0 && out_set.count(g.AllEdge[i][1])==0 && source!=g.AllEdge[i][1])
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}

	}
	//certainEdge.HangEdges.insert(0);
	return;
}


/*获取确定的边，包括悬挂边和割集中的边*/
/*参数说明：g:不确定图；source:顶点；sink:终点； certainEdge：确定边集合的引用 */
void getCertainEdges(Graph& g, int source, int sink, CertainEdge& certainEdge)
{
	/*获取割集中的边*/
	getMinCutEdges(g,source,sink,certainEdge);
	/*获取悬挂边*/
	getHangEdge(g,source,sink,certainEdge);

	return;
}
