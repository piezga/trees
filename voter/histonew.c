#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define MAXDATA 100000
int main(int argc, char *argv[]){
  int nbin,i,ib,j;
  int ndata=MAXDATA,ntmp;
  float *h,x,*X,dd;
  float min,max,med,var;
  FILE *fout,*fin;
  char fname[512],outname[512];

  if(argc==3){
    sscanf(argv[1],"%d",&nbin); 
    sscanf(argv[2],"%s",&outname); 
  }
  else{
    fprintf(stderr,"usage histo <nbin> ,<outname>\n");
    exit(99);
  }
  fprintf(stderr,"Evaluating pdf and cumulative from stdinput\n");
  h = (float *)malloc(nbin*sizeof(float));
  X = (float *)malloc(MAXDATA*sizeof(float));
  for(i=0;i<nbin;i++)
    h[i]=0.;
  
  
  min=1.0e10;
  max=-1.0e10;
  med=0.0;
  var=0.0;
  
  i=0;
  while(scanf("%g\n", &x )==1){
    if(i>ndata){
      ntmp=ndata+100000;
      X=realloc(X,ntmp*sizeof(float));
      if(X==NULL){
	fprintf(stderr,"Not enough memory to allocate X \n");
	exit(9);
      }
      ndata=ntmp;
    }
    
    X[i]=x;
    med += x;
    var += x*x;
    if(x<min)min=x;
    if(x>max)max=x;
    i++;    
  }
  ndata=i;
  fprintf(stderr,"datapoint=%d\n",ndata);
  med /= (float) ndata;
  var /= (float) ndata;
  var= sqrt(var-med*med);
 

  dd=(max-min)/(float) (nbin);
  
  fprintf(stderr,"average, variance\n");
  fprintf(stderr,"%f %f\n", med,var);
  fprintf(stderr,"max,min,bin\n");
  fprintf(stderr,"%f %f %f\n", max,min,dd);
  
  for(i=0;i<ndata;i++){
    ib=(X[i]-min)/dd;
    if(ib<0||ib>nbin)exit(99);
    h[ib]+=1.0;
  }
  
  sprintf(fname,"%s.histo",outname);
  fout = fopen(fname,"w");
  fprintf(fout,"#datapoint=%d\n",ndata);
  fprintf(fout,"#average, variance\n");
  fprintf(fout,"#%f %f\n", med,var);
  fprintf(fout,"#max,min,bin\n");
  fprintf(fout,"#%f %f %f\n", max,min,dd);
  fprintf(fout,"##############################\n");
  med=0.0;
  for(ib=0;ib<nbin;ib++){
    med+=h[ib]/(dd*ndata);
    fprintf(fout,"%g %g %g\n", min+dd*(ib)+dd*.5, h[ib]/(dd*ndata),med*dd);
  }
  fprintf(stderr,"PDF normlized %f\n", med*dd);
  fclose(fout);
  
  exit(1);
}
