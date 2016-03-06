#include <set>
#include <vector>

/*����ͼ��һЩ�㷨*/
#include "graph.h"
#include "cut.h"
using namespace std;



/*��ȡȷ��ͼ��ĳ���������Ӧ��еı�*/
/*����˵����g:��ȷ��ͼ��source:���㣻sink:�յ㣻 certainEdge��ȷ���߼��ϵ����� */
int getMinCutEdges(Graph& g, int source, int sink,  CertainEdge& certainEdge)
{
	/*���ʣ��ͼ,����gfΪʣ��ͼ*/
	Flow f;
	GF gf;
	int max_flow = Dinic(g,source,sink,f,gf);

	/*Դ���㼯*/
	set<int> S;
	/*ʹ��ʣ��ͼ��ȡԴ���㼯��*/
	MinCut(gf,g.nV, source, S);
	/*�����S�У��յ㲻��S�еı�Ϊ�ؼ���*/
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

/*��ȡȷ���ıߣ����ұ�*/
/*����˵����g:��ȷ��ͼ��source:���㣻sink:�յ㣻 certainEdge��ȷ���߼��ϵ����� */
void getHangEdge(Graph& g, int source, int sink, CertainEdge& certainEdge)
{
	set<int> in_set,out_set;/*���in,�յ�out*/
	for(int i=1; i<=g.nE; i++)
	{
		in_set.insert(g.AllEdge[i][1]);
		out_set.insert(g.AllEdge[i][2]);
	}

	for(int i=1; i<=g.nE; i++)
	{
		/*Դ������ı������ұ�*/
		if (g.AllEdge[i][2] == source)
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}
		/*��������ı������ұ�*/
		if (g.AllEdge[i][1] == sink)
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}

		/*ֻ�����룬û�������ĵ�*/
		if (in_set.count(g.AllEdge[i][2])==0 && out_set.count(g.AllEdge[i][2])>0 && sink != g.AllEdge[i][2])
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}

		/*ֻ��������û������ĵ�*/
		if (in_set.count(g.AllEdge[i][1])>0 && out_set.count(g.AllEdge[i][1])==0 && source!=g.AllEdge[i][1])
		{
			certainEdge.HangEdges.insert(g.AllEdge[i][4]);
			certainEdge.num_HangEdges++;
		}

	}
	//certainEdge.HangEdges.insert(0);
	return;
}


/*��ȡȷ���ıߣ��������ұߺ͸�еı�*/
/*����˵����g:��ȷ��ͼ��source:���㣻sink:�յ㣻 certainEdge��ȷ���߼��ϵ����� */
void getCertainEdges(Graph& g, int source, int sink, CertainEdge& certainEdge)
{
	/*��ȡ��еı�*/
	getMinCutEdges(g,source,sink,certainEdge);
	/*��ȡ���ұ�*/
	getHangEdge(g,source,sink,certainEdge);

	return;
}
