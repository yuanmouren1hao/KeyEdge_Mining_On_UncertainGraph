#ifndef CUT_H
#define CUT_H

#include <set>

/*�����ͼSplitGraph*/
typedef struct{

}SplitGraph; 




/*��ȡȷ��ͼ��ĳ���������Ӧ��еı�*/
/*����˵����g:��ȷ��ͼ��source:���㣻sink:�յ㣻 MinCutEdges����߼��ϵ����� */
int getMinCutEdges(Graph& g, int source, int sink, set<int>& MinCutEdges);

#endif