#include <stdio.h>
#include <stdlib.h>
#define MAXCOL 255

int main(int argc, char *argv[]){
int NX,NY;
 int kx,ky,max;
char fname[128];
int **map;
FILE *fin;
int seme1,seme2,seme3;
int R[256],G[256],B[256];
 if(argc!=4){
   fprintf(stderr, "usage %s <nome file> <NX> <NY>\n");
   exit(9);
 }



sprintf(fname,"%s",argv[1]);
sscanf(argv[2],"%d",&NX);
sscanf(argv[3],"%d",&NY);


 map = (int **) malloc(sizeof(int *) * NX );
for (kx=0;kx<NX;kx++)
  map[kx] = (int *) malloc(sizeof(int) * NY );





fin = fopen(fname,"r");
 fprintf(stderr, "open %s\n",fname);
for(kx=0;kx<NX;kx++)
  fread(map[kx], sizeof(int), NY, fin);

fclose(fin);
 fprintf(stderr, "%s has been read\n",fname);
 max=0;
for(kx=0;kx<NX;kx++){
  for(ky=0;ky<NY;ky++){
    if(map[kx][ky]>max)max=map[kx][ky];
  }
 }
  fprintf(stderr, "max is %d\n",max);
 seme1=0.0*6832121;
 seme2=895973;
 seme3=23234;

 for(kx=0;kx<=max;kx++){
   R[kx]= (seme1*kx)%256; 
   G[kx]= (seme2*kx)%256; 
   B[kx]= (seme3*kx)%256; 
 }
 

fprintf(stdout, "P3\n");
fprintf(stdout, "%d %d\n",NX,NY);
fprintf(stdout, "%d\n",MAXCOL);

for(kx=0;kx<NX;kx++){
  for(ky=0;ky<NY;ky++){
    fprintf(stdout, " %d %d %d  \n",R[map[kx][ky]],G[map[kx][ky]],B[map[kx][ky]]);
  }
 }

  exit(1);
}
