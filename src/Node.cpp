/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr> 
include colibri functions (previously in cghseg)
*/
#include "Node.h"


Node::Node(int I, double V, int LI, int HI)
{
	Value = V;
	Index = I;
	LowIndex = LI;
	HighIndex = HI;
}

bool Node::operator<(Node &Other)
{
	return (Value < Other.Value);
}

bool Node::operator<=(Node &Other)
{
	return (Value <= Other.Value);
}

Node::Node()
{
	Value = 0;
	Index = 0;
	LowIndex = 0;
	HighIndex = 0;	
}


Node Node::operator=(Node &Other)
{
	Value = Other.Value;
	Index = Other.Index;
	LowIndex = Other.LowIndex;
	HighIndex = Other.HighIndex;
	return *this;
}

void Swap(Node &a, Node &b)
{
	Node c = a;
	a = b;
	b = c;
}

