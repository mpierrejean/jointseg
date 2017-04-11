
/* March 2017 Guillem Rigaill <guillem.rigaill@inra.fr> 
include colibri functions (previously in cghseg)
*/
#include "colibri.h"
#include<R_ext/Arith.h>

// this function is visible by R
extern "C" {

void colibri_sn_R_c (double *profil, int *nbi, int *Kmax_, double *mini, double *maxi, int *origine,
double *cout_n, double *allCost){
	colibri_sn_c (profil, nbi, Kmax_, mini, maxi, origine, cout_n, allCost);
  }

}

