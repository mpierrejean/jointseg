/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr> 
include colibri functions (previously in cghseg)
*/

#include <list>
#include <iostream>
#include <fstream>
#include "Node.h"

#ifndef _Heap_
#define _Heap_


using namespace std;

class Heap
{
public:
  Node *MyHeap;
  int HeapSize;
  //Heap(Node *TheN = NULL, int NbNodes = 0, int AllocationSize = 0);
  Heap(int AllocationSize);
  Heap();
  ~Heap();
  void AddNode(Node N);
  void RemoveHead();
private:
  //void Debug();
	int AllocatedSize;
	void ReAllocate();
};

//ostream & operator<<(ostream &s, const Heap &H);



#endif



