#include <set>
#include <vector>

/*����ͼ��һЩ�㷨*/
#include "graph.h"
using namespace std;



/*��ȡȷ��ͼ��ĳ���������Ӧ��еı�*/
/*����˵����g:��ȷ��ͼ��source:���㣻sink:�յ㣻 MinCutEdges����߼��ϵ����� */
int getMinCutEdges(Graph& g, int source, int sink, set<int>& MinCutEdges)
{
	/*���ʣ��ͼ,����gfΪʣ��ͼ*/
	Flow f;
	GF gf;
	int max_flow = Dinic(g,source,sink,f,gf);

	/*Դ���㼯*/
	set<int>S;
	/*ʹ��ʣ��ͼ��ȡԴ���㼯��*/
	MinCut(gf,g.nV, source, S);
	/*�����S�У��յ㲻��S�еı�Ϊ�ؼ���*/
	for (int i=1; i<=g.nE; i++)
	{
		if (S.count(g.AllEdge[i][1]) && !S.count(g.AllEdge[i][2]))
		{
			MinCutEdges.insert(g.AllEdge[i][4]);
		}
	}

	return 0;
}


