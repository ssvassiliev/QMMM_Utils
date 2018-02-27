#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#define MAXPRM 10000
int main(int argc, char** argv) 
{
  FILE *fp;
  char line_buf[200];
  int i,j,N;
  char* ambtrs="torsions.dat";
  char atoms[MAXPRM][12];
  float IDIVF[MAXPRM]; int PN[MAXPRM];
  float PK[MAXPRM][4],PHASE[MAXPRM][4];
  float tmpPHASE, tmpPK;

/* AmbTrs CT CT N  C    0   0 180 180  0.530  0.000  0.150  0.500 1.0 */

/*           IDIVF    PK        PHASE           PN */
/* CT-CT-N -C    1    0.50        180.0            -4.         phi,psi,parm94 */
/* CT-CT-N -C    1    0.15        180.0            -3.         phi,psi,parm94 */
/* CT-CT-N -C    1    0.00          0.0            -2.         JCC,7,(1986),230 */
/* CT-CT-N -C    1    0.53          0.0             1.         phi,psi,parm94 */


  fp=fopen(ambtrs, "rt");
  if(fp==NULL)
    {printf("** E ** torsions file not found\n");return(1);}

  i=0;
  while(1)
    {
      if(fgets(line_buf,200,fp)==NULL)
	break;
      if(PN[i]>=0)
	strncpy(atoms[i],line_buf,11);

      if(atoms[i][0]=='X')
	   atoms[i][0]='*';   
      if(atoms[i][9]=='X')
	   atoms[i][9]='*';   
      for(j=0;j<11;j++)
	{
	  if(atoms[i][j]=='-')
	    atoms[i][j]=' ';   
	}
      sscanf(&line_buf[12], "%f%f%f%i",&IDIVF[i],&tmpPK,&tmpPHASE,&PN[i]);
      PK[i][(int)abs(PN[i])-1]=tmpPK;
      PHASE[i][(int)abs(PN[i])-1]=tmpPHASE;    
      if(PN[i]>0)
	i++;
    }
  N=i-1;


  fclose(fp);


  printf("!\n!Torsions\n!\n");
 for(i=0;i<N;i++)
   printf("AmbTrs %s %5i%5i%5i%5i %7.3f%7.3f%7.3f%7.3f%4.1f\n", atoms[i],\
	  (int)(PHASE[i][0]), (int)PHASE[i][1], (int)PHASE[i][2],(int)PHASE[i][3], \
	  PK[i][0], PK[i][1], PK[i][2], PK[i][3], IDIVF[i] );
 return(0);
}
