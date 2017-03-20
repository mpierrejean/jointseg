
/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr> 
include colibri functions (previously in cghseg)
*/
#include "liste.h"

void Liste::show()
{
        //std::cout << "Max : " << this->getMax() << ", Min : " << this->getMin()<< std::endl;
        this->poly->show();
}
void Liste::showAllNext()
{
        Liste *l;
        l = this;
        while(l != NULL)
        {
                l->show();
                l=l->getNext();
        }

}
void Liste::resetMaillonBorders(Polynome2 *poly_)
{
	//if(this->getPolynome()->getRacine2() == NAN)
	if(this->getPolynome()->getRacine2() == 0.)
	//if( isnan(this->getPolynome()->getRacine2()) )
	{
		this->setPolynome(poly_);
	} else if(this->getPolynome()->getRacine1() >= this->getMax())
	{
		if(this->getPolynome()->getRacine2() >= this->getMax())
		{
			this->setPolynome(poly_);
		} else if(this->getPolynome()->getRacine2() > this->getMin()) 
			{
				Liste * maillon= new Liste(this->getPolynome()->getRacine2(), this->getMin(), poly_);
				this->insert(maillon);
				this->setMin(this->getPolynome()->getRacine2());
			} else 
			{
			}
		
	
	} else if( this->getPolynome()->getRacine1() > this->getMin())
		{
		
			if(this->getPolynome()->getRacine2() > this->getMin())
			{
				Liste *maillon3 = new Liste(this->getPolynome()->getRacine2(), this->getMin(), poly_);
				Liste *maillon2 = new Liste(this->getPolynome()->getRacine1(), this->getPolynome()->getRacine2(), this->getPolynome());
				this->setMin(this->getPolynome()->getRacine1());
				this->setPolynome(poly_);
				this->insert(maillon3);
				this->insert(maillon2);
			} else
			{
				Liste *maillon2 = new Liste(this->getPolynome()->getRacine1(), this->getMin(), this->getPolynome());
				this->setMin(this->getPolynome()->getRacine1());
				this->setPolynome(poly_);
				this->insert(maillon2);
			}
		
		} else 
		{
			this->setPolynome(poly_);
		}
}
void Liste::resetAllBorders(Polynome2 *poly_)
{

	Liste *lCurrent, *lNext;
	lCurrent = this;
	lNext = this->getNext();
	lCurrent->resetMaillonBorders(poly_);
	lCurrent=lNext;
	while(lCurrent != NULL)
	{
		lNext=lCurrent->getNext();	
		lCurrent->resetMaillonBorders(poly_);
		lCurrent=lNext;
	}
}


