#ifndef CUT_H
#define CUT_H

#include <set>

/*�����ͼSplitGraph*/
typedef struct{

}SplitGraph; 


/*ȷ���ߵļ��ϣ�������еıߺ����ұ�*/
typedef struct{
	set<int> MinCutEdges;/*��еı�*/
	set<int> HangEdges;/*���ұ�*/
}CertainEdge;


/*��ȡȷ���ıߣ��������ұߺ͸�еı�*/
/*����˵����g:��ȷ��ͼ��source:���㣻sink:�յ㣻 certainEdge��ȷ���߼��ϵ����� */
void getCertainEdges(Graph& g, int source, int sink, CertainEdge& certainEdge);

#endif