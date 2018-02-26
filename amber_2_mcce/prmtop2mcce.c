/* File prmtop2mcce.c */
/* This program reads AMBER parm7 topology file and converts it to MCCE  */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../File_Utilities/file_read_write.h"
#include "../File_Utilities/file_read_write.c"

char at_type[1000][3];
float VDW_RAD[1000];
float VDW_EPS[1000];
int n_types;

void write_atype(char *label)
{
  int i,j;
  FILE *fp;

  fp = fopen(label,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}  
  printf("Writing MCCE force field file: %s\n", label);

  for(i=0;i<N;i++)
    for(j=0;j<n_types;j++)
      {
	if(!strncmp(Type[i],at_type[j],2))
	  {
	    if(strcspn(atom_name[i]," ")==4)
	      {
		fprintf(fp,"VDW_RAD  %s01 %s %f\n",  RESIDUE_LABEL[0], atom_name[i],VDW_RAD[j]);
		fprintf(fp,"VDW_EPS  %s01 %s %f\n",  RESIDUE_LABEL[0], atom_name[i],VDW_EPS[j]);
		break;
	      }
	    else
	      {
		fprintf(fp,"VDW_RAD  %s01  %s%f\n",  RESIDUE_LABEL[0], atom_name[i],VDW_RAD[j]);
		fprintf(fp,"VDW_EPS  %s01  %s%f\n",  RESIDUE_LABEL[0], atom_name[i],VDW_EPS[j]);
		break;
	      }

	    fprintf(fp,"VDW_RAD  %s01 %s %f\n",  RESIDUE_LABEL[0], atom_name[i],VDW_RAD[j]);
	    fprintf(fp,"VDW_EPS  %s01 %s %f\n",  RESIDUE_LABEL[0], atom_name[i],VDW_EPS[j]);
	    break;
	  }
      }
  fclose(fp);
  
}



int read_amber_params(char *forcefield)
{
  FILE *fp;
  long i,j,k;
  char line_buf[200];  

  fp=fopen(forcefield, "rt");
  if(fp==NULL)
    {printf("** E ** forcefield file not found\n");return(0);}
  fgets(line_buf,82,fp);


  /* ----------------- READ NONBONDED PARAMS ---------------------*/
  while(1)
    {
      fgets(line_buf,82,fp);
      if(!strncmp(line_buf,"NONB",4))
	break;
    }
  printf("FF NONB Section found\n");
  k=0;
  while(1)
    {
      fgets(line_buf,160,fp);
      if(line_buf[0]=='\n')
	break;
	 strncpy(at_type[k],&line_buf[0],2);
	 at_type[k][2]='\0';
	 sscanf(&line_buf[3],"%f %f",&VDW_RAD[k],&VDW_EPS[k]);
	 k++;
    }

  n_types=k;

      for(i=0;i<n_types;i++)
	printf("%s %f %f\n",at_type[i],VDW_RAD[i],VDW_EPS[i]);

      return(0);
}


int main(int argc, char **argv)
{  
  char label[80], label_ff[80];
   
  printf("\n **** MD utilities, (c) S.Vasil'ev, Brock U. Biology ****\n\n");  

  if(!read_parm7(argv[1]))
    exit(0);
  read_amber_params(argv[2]);
  convert_amb_types_to_radii();
  strcpy(label,RESIDUE_LABEL[0]);
  strcpy(label_ff,RESIDUE_LABEL[0]);
  strcat(label,".tpl");
  strcat(label_ff,".ffm");
  write_mcce(label);
  write_atype(label_ff);
  return(1);
}

