#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "functionsdef.h"
#define NMAX 100


int scala[NMAX];


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

********************* FUNZIONI PER ANALISI STATISTICA ******************** 

<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

int generascale(double erre){
int k,m,nn;
double kprova;


 if(NX!=NY) {
   fprintf(stderr, "ATTENTION: analysis is working only for NX = NY\n");
   exit(1);
 }

 kprova=1;
 m=0; k=2;
 scala[m]=k;

 while(k<NX){
   nn=0;
   while((int) kprova<=k){
     kprova=2*pow(erre,nn);
     nn++;
   }
   m++;
   k=(int) kprova;
   scala[m]=k;
 }
 if(scala[m-1]!=NX)scala[m++]=NX;

 return m;


}



void computeSAR(int **MS, double *sar, double *sarvar, int numscale, int irel, int numrel){
  int k,kx,ky,jx,jy,j,lsquares,count;
  int *speciesvector;
  char nomefi[128];
  FILE *fout;
  double scra;
  speciesvector = malloc(sizeof(int)*(NS));

  sprintf(nomefi,"SAR_NU%e_NX%d_NY%d_numrel%d.ker%d",NU,NX,NY,numrel,kernsize);
  fout =fopen(nomefi,"w");
  fprintf(fout,"##nu=%g NX=%d  NY=%d numrel=%d kernsize=%d\n",NU,NX,NY,irel, kernsize);  


   for(k=0;k<numscale;k++){   
     lsquares=NX/scala[k]; 
     for(kx=0;kx<lsquares;kx++){
       for(ky=0;ky<lsquares;ky++){

	 for(j=0;j<NS;j++)  speciesvector[j]=0;
	 
	 for(jx=kx*scala[k];jx<(kx+1)*scala[k];jx++)
	   for(jy=ky*scala[k];jy<(ky+1)*scala[k];jy++){
	     speciesvector[MS[jx][jy]]++;       
	   }    
	 count=0;
	 for(j=0;j<NS;j++) if (speciesvector[j]!=0) count++;
	 scra = (double)count;
	 sar[k]    += scra;
	 sarvar[k] += scra*scra;	 
	 
       }       
     }
     {double var=0,med=0,norm=0;
       norm = 1./((double) (lsquares*lsquares) * irel);
       med=sar[k]*norm;
       var=(sqrt( norm* sarvar[k] - med*med));
       var *=sqrt(norm);
       fprintf(fout,"%d %lf %lf\n",scala[k]*scala[k],med,var);
     }
   }
   
   fclose(fout);

}




void computeSAD(int **MS,double *NxS){
  int k,kx,ky;
  int *numberofindividuals;
 

  numberofindividuals = malloc(sizeof(int)*(NS));


  for(k=0;k<NS;k++)
    numberofindividuals[k]=0;



  /* conta numero di individui per ogni specie */
  for(kx=0;kx<NX;kx++){
    for(ky=0;ky<NY;ky++){
      numberofindividuals[MS[kx][ky]]++; 
    }
  }

  for(k=0;k<NS;k++)
    NxS[numberofindividuals[k]]++;

  /*
  for(k=0;k<NS;k++){
    fprintf(stderr,"%d %d \n",k,numberofindividuals[k]);
  }
  */
  

  free(numberofindividuals);
  
}


/* ******************** PLOT SAD HISTOGRAM ********************* */


/* commento: non e' ancora il prestonplot per ora e' un semplice istogramma binnato logaritmico */
/* TODO: implement true preston plot */

double prestonplot(double *NxS, int NREL){
  int k,nbin,ib;
  double *h;
  double base,norm,total;
  char nomefi[128];
  FILE *fplot;

  base=1.0/log(2.0);
  nbin=(int) log((double) NTOT)*base;
  norm=1./(double) NREL;


  h=(double *) malloc(sizeof(double) * (nbin+1));
  for(ib=0;ib<nbin;ib++)
    h[ib]=0.0;

  for(k=0;k<NTOT;k++){
    ib=(int) log((double) (k+1))*base;
    if(ib>=nbin+1)fprintf(stderr, "ERROR ib=%d nbin=%d\n",ib,nbin);
    h[ib]+=NxS[k]*norm;
  }


  sprintf(nomefi,"PPLOT_NU%e_NX%d_NY%d_NREL%d.ker%d",NU,NX,NY,NREL,kernsize);
  fplot =fopen(nomefi,"w");
  fprintf(fplot,"##nu=%g NX=%d  NY=%d numrel=%d kernsize=%d\n",NU,NX,NY,NREL, kernsize);  
  
  total=0;
  for(ib=0;ib<nbin;ib++){
    total += h[ib];
    fprintf(fplot,"%d %g\n",ib,h[ib]);
  }
  fclose(fplot);

  free(h);
  return total;
}

/* ******************* PRINT CONFIGURATION ********************* */

void printconf(int irel, int **matrice){
int kx,ky;
FILE *fdump;
char fname[128];

#ifdef BINARY
   if(kernsize>0)
     sprintf(fname,"BConf_NX%d_NY%d_NU%e_Kernsize%d_rel.%d",NX,NY,NU,kernsize,irel);
   else
     sprintf(fname,"BConf_NX%d_NY%d_NU%e_KernNN_rel.%d",NX,NY,NU,irel);
#else
   if(kernsize>0)
     sprintf(fname,"AConf_NX%d_NY%d_NU%e_Kernsize%d_rel.%d",NX,NY,NU,kernsize,irel);
   else
     sprintf(fname,"AConf_NX%d_NY%d_NU%e_KernNN_rel.%d",NX,NY,NU,irel);
#endif
   fdump = fopen (fname,"w");
   
#ifdef BINARY
   for(kx=0;kx<NX;kx++)
     fwrite(matrice[kx], sizeof(int), NY, fdump);   
#else
   for(kx=0;kx<NX;kx++)
     for(ky=0;ky<NY;ky++)
       fprintf(fdump, "%d %d %d\n",kx,ky,matrice[kx][ky]);
#endif
   fclose(fdump);

}


