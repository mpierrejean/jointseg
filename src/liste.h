
/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr>
include colibri functions (previously in cghseg)
*/
#ifndef LISTEH
#define LISTEH
#include "polynome2.h"

//#include <math.h> // not needed as included by polynome2.h
#include <R_ext/Arith.h>
#define NUMLIB_POSINF R_PosInf

class Liste {
private:
  double max, min;
  Liste *next;
  Polynome2 *poly;

public:
  /* constructors and destructors */
  Liste()
  : max(0.), min(0.), next(NULL), poly(NULL) {}
  Liste(double max_, double min_)
  : max(max_), min(min_), next(NULL) , poly(NULL){}
  Liste(double max_, double min_, Polynome2 *poly_)
  : max(max_), min(min_), next(NULL), poly(poly_) {}
  Liste(Polynome2 *poly_)
  : max(0.), min(0.), next(NULL), poly(poly_) {}
  ~Liste(){
    delete next;
    delete poly;
  }
  /* fonction setter and getter */
  inline
  double getMax();
  inline
  void setMax(double max_);
  inline
  double getMin();
  inline
  void setMin(double min_);
  inline
  void setPolynome(Polynome2 * poly_);
  inline
  Polynome2 *getPolynome();
  inline
  Liste * getNext();
  inline
  void setNext(Liste * next_);
  /* Useful */
  inline
  void setToNull();
  inline
  void insert(Liste * maillon_);
  inline
  int compteIntervals();
  inline
  int compteCand(int *index, int step);
  inline
  Liste * removeDoublon();
  inline
  void checkForDoublon();
  /* show and others */
  void show();
  void showAllNext();
  /* */
  inline
  void computeRoots(double);
  inline
  void add(double, double, double);
  inline
  void computeMinOrMax(double*, int*);
  void resetMaillonBorders(Polynome2*);
  void resetAllBorders(Polynome2*);
};
/* Setter and Getter */
/* */
double Liste::getMax()
{
  return(this->max);
}
void Liste::setMax(double max_)
{
  this->max = max_;
}
double Liste::getMin()
{
  return(this->min);
}
void Liste::setMin(double min_)
{
  this->min = min_;
}
void Liste::setPolynome(Polynome2 * poly_)
{
  this->poly=poly_;
}
Polynome2* Liste::getPolynome()
{
  return(this->poly	);
}
/* */
Liste * Liste::getNext()
{
  return(this->next);
}
void Liste::setNext(Liste * next_)
{
  this->next = next_;
}
void Liste::setToNull()
{
  max=0.;
  min=0.;
  poly=NULL;
  next=NULL;
}
void Liste::insert(Liste * maillon_)
{
  maillon_->setNext(this->getNext());
  this->setNext(maillon_);
}
int Liste::compteIntervals(){
  Liste *l;
  int tmp = 0;
  l=this;
  while(l != NULL){
      tmp= tmp+1;
      l=l->getNext();
  }
  return(tmp);
}


int Liste::compteCand(int *index, int step){
  Liste *l;
  int tmp = 0;
  l=this;
  while(l != NULL){
      if(index[l->getPolynome()->getOrigine()] < step)
	{
		index[l->getPolynome()->getOrigine()] = step;
		tmp = tmp+1;
	}
      l=l->getNext();
  }
  return(tmp);
}

Liste * Liste::removeDoublon()
{
  Liste *next = this->getNext();
  if(next != NULL)
    {
      if(next->getPolynome() == this->getPolynome())
        {
          //std::cerr<<"erase"<<std::endl;
          this->setMin(next->getMin());
          this->setNext(next->getNext());
          next->setToNull();
          delete next;
          return(this);
        } else
          {
            return(next);
          }
    } else
      {
        return(NULL);
      }
}
void Liste::checkForDoublon()
{
  Liste *l = this;
  while(l != NULL)
    {
      l=l->removeDoublon();
    }
}
void Liste::computeRoots(double a0_)
{
  Liste *l;
  l=this;
  while(l != NULL)
    {
      l->getPolynome()->roots(a0_);
      l=l->getNext();
    }
}
void Liste::add(double a2_, double a1_, double a0_)
{
  Liste *l;
  l=this;
  while(l != NULL)
    {
      l->getPolynome()->add(a2_, a1_, a0_);
      l=l->getNext();
    }
}
void Liste::computeMinOrMax(double * min, int * which)
{
  Liste *l;
  double tmp = NUMLIB_POSINF;
  *min = NUMLIB_POSINF;
  * which=-1;
  l=this;
  while(l != NULL)
    {
      l->getPolynome()->minOrMax(min, &tmp, which);
      l=l->getNext();
    }
}
#endif //INLINED

