/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr> 
include colibri functions (previously in cghseg)
*/
#ifndef _Node_h_
#define _Node_h_


class Node
{
public:
	int Index;
	double Value;
	int LowIndex;
	int HighIndex;
	Node(int I, double V, int LI, int HI);
	Node();
	bool operator<(Node &Other);
	bool operator<=(Node &Other);
	Node operator=(Node &Other);
};

void Swap(Node &a, Node &b);


#endif


