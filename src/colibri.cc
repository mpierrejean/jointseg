/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr>
include colibri functions (previously in cghseg)
*/

#include "colibri.h"

////////////////////

void colibri_sn_c (const double *profil, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, int *origine,
double *cout_n, double *allCost)
{
	int nb     = *nbi;
	int Kmax   = *Kmaxi;
	double min = *mini;
	double max = *maxi;

	double *minCostBefore  = new double[nb];
	double *minCostCurrent = new double[nb];

	double *tmp; //1
	int minPosition;
	double minCurrent;
	//int * origine = (int *) malloc(nb * sizeof(int));
	int i = 0;
	int i2 = 0;
	double somme = 0;
	double sommeC = 0;
	int turn = 1;
	//char c = 13;

    /* Initialisation Cout en 1 segment */
    while(i < nb)
	{
		somme = somme + profil[i];
		sommeC = sommeC + profil[i]*profil[i];
		minCostBefore[i] = sommeC - pow(somme, 2) / (i+1);
		origine[i]=0;
		allCost[i] = minCostBefore[i];
		i++;
	}

    /* Save */
    cout_n[0] = minCostBefore[nb-1];


    /* Initialisation Polynome Cost */
    //Polynome2 * p1;
	Liste * l1;
	//Polynome2 * pTest;

	Polynome2 **stock = new Polynome2* [nb]; 

        i=0;
	while(i < nb)
	{
		stock[i] = new Polynome2();
		i++;	
	}


    /* Boucle turn 1 -> Kmax -1 */
    while( turn < Kmax)
    {
      /* initalisation */
      i  = turn;
      i2 = turn+ turn*nb;
      stock[i]->reset(1.0, -2*profil[i], profil[i]*profil[i]+minCostBefore[turn -1],  turn);
      stock[i]->setStatus(2);
      l1 = new Liste(max, min, stock[i]);

      /* Min */
      l1->computeMinOrMax(&minCurrent, &minPosition);
      minCostCurrent[i] = minCurrent;
      origine[i2]       = i;
      allCost[i2]       = minCurrent;

      /* iterate */
      i++;
      i2++;

      while(i < nb)
      {
	/* Slide 1 and Prune */
	   l1->computeRoots(minCostBefore[i-1]);
	   stock[i]->reset(0.0, 0.0, minCostBefore[i-1],  i);
	   l1->resetAllBorders(stock[i]);
	   l1->checkForDoublon();
	   l1->add(1.0, -2*profil[i], profil[i]*profil[i]);

	   /* Compute Min */
  	   l1->computeMinOrMax(&minCurrent, &minPosition);
   	   minCostCurrent[i]=minCurrent;
	   origine[i2] = minPosition;
	   allCost[i2] = minCurrent;
	   /* iterate */
	   i++;	
           i2++;
      }

      /* Save */
      cout_n[turn] = minCostCurrent[nb-1];
	  
     /* */
     tmp=minCostCurrent;
     minCostCurrent=minCostBefore;
     minCostBefore=tmp;
	
	
     //delete(l1);
     /* iterate */
     turn++;

   }
	
   /* Free All */
   /* free stock */
   i=0;
   while(i < nb)
   {
     delete(stock[i]);	
     i++;
   }
   delete[] stock;  
   delete[] minCostBefore;
   delete[] minCostCurrent;
}
///////////////////////////////////////


