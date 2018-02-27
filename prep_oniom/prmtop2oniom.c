/* File prmtop2mcce.c */
/* This program reads AMBER parm7 topology file and converts it to MCCE  */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../File_Utilities/file_read_write.h"
#include "../File_Utilities/file_read_write.c"

char params[10][100];
int prmcnt;
int h_list[MAX_NATOMS][6];
int h_num[MAX_NATOMS];
char linkAtom[MAX_NATOMS][5];
char linkType[MAX_NATOMS][5];

void read_freeze(char *restr, int *at_freeze, char *at_level, int *link)
{
  FILE *fp;
  char line_buf[200];
  char aname[5], resname[5];  
  int i,j,k,start,stop,nrest,nlink,bnd2;
  char bondedTo[5], lAtm[5], lTyp[5];


  int cf=1; char level;

  fp=fopen(restr, "rt");
  if(fp==NULL)
    {printf("** E ** restraints file not found\n");return;}

 for(i=0;i<N;i++)
   for(j=0;j<5;j++)
     if(atom_name[i][j]==' ')
       {(atom_name[i][j]='\0');break;}

  // MAKE LISTS OF BONDED HYDROGEN ATOMS
  for(i=0;i<N;i++)
    {
      h_num[i]=0;
      // BONDS WITH HYDROGEN     
      for(j=0;j<NBONH;j++)
	{
	  if(IBH[j]==i)
	    h_list[i][h_num[i]++]=JBH[j];
	  if(JBH[j]==i)
	    h_list[i][h_num[i]++]=IBH[j];
	}  
    }

  prmcnt=0;
  for(;;)
    {
      if(fgets(line_buf,82,fp)==NULL)
	break;
      nrest=0;nlink=0;     
      /*-------------------------FF PARAMETER FILES------------------------------*/
      if(!strncmp(line_buf,"PARAMS",6))
	{ sscanf(&line_buf[7], "%s", &params[prmcnt++]);continue;}

      /*-------------------------MODEL LEVELS------------------------------------*/
      if(!strncmp(line_buf,"LEVEL ALL",9))
	{
	  sscanf(&line_buf[20], "%c", &level);
	  for(i=0,nrest=0;i<N;i++)
	    {nrest++; at_level[i]=level;}
	  printf("Model level  %c set for %6i atoms: (all)\n",level,nrest);continue;
	}
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"LEVEL INDEX ",12))
	{
	  sscanf(&line_buf[20], "%i%i %c", &start, &stop, &level);
	  for(i=0,nrest=0;i<N;i++)
	    if((i>=start-1)&&(i<=stop-1))
	      {nrest++; at_level[i]=level;}
	  printf("Model level  %c set for %6i atoms: index %i-%i\n", level,nrest,start,stop);continue;
	}   
      if(!strncmp(line_buf,"LEVEL INDEXH",12))
	{
	  sscanf(&line_buf[20], "%i%i %c", &start, &stop, &level);
	  for(i=0,nrest=0;i<N;i++)
	    if((i>=start-1)&&(i<=stop-1))
	      {
		nrest++; at_level[i]=level;
		for(k=0;k<h_num[i];k++)
		  at_level[h_list[i][k]]=level;
	      }
	  printf("Model level  %c set for %6i atoms: index %i-%i\n", level,nrest,start,stop);continue;
	}   
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"LEVEL RESID ",12))
	{
	  sscanf(&line_buf[20], "%i%i %c", &start, &stop, &level);
	  for(i=0,nrest=0;i<N;i++)
	    if((resSeq[i]>=start-1)&&(resSeq[i]<=stop-1))
	      {nrest++; at_level[i]=level;}
	  printf("Model level  %c set for %6i atoms: resid %i-%i\n", level,nrest,start,stop);continue;
	}   
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"LEVEL ANAME ",12))
	{
	  sscanf(&line_buf[20], "%s %c",aname, &level);
	  for(i=0,nrest=0;i<N;i++)
	    if(!strncmp(&atom_name[i][0], aname, strlen(aname)))
	      {nrest++; at_level[i]=level;}
	  printf("Model level  %c set for %6i atoms: aname %s*\n", level, nrest, aname);continue;
	}
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"LEVEL RNAME ",12))
	{
	  sscanf(&line_buf[20], "%s %c", resname, &level);
	  for(i=0,nrest=0;i<N;i++)
	    if(!strncmp(resName[i], resname, strlen(resname)))
	      {nrest++; at_level[i]=level;}
	  printf("Model level  %c set for %6i atoms: rname %s\n", level, nrest, resname);continue;
	}
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"LEVEL RNAME:ATOM",16))
	{
	  sscanf(&line_buf[20], "%s %i %i %c", resname, &start, &stop, &level);
	  for(j=0;j<NRES;j++)
	    {
	      if(!strncmp(RESIDUE_LABEL[j],resname,strlen(resname)))
		for(i=RESIDUE_POINTER[j]-1+start-1;i<=RESIDUE_POINTER[j]-1+stop-1;i++)
		  {
		    nrest++; 
		    at_level[i]=level;
		  }
	    }
	  printf("Model level  %c set for %6i atoms: rname %s:%i-%i\n", level, nrest, resname,start,stop);
	  continue;
	}
    /*-------------------------------------------------------------*/
      nrest=0;
      if(!strncmp(line_buf,"LEVEL RNAME:ANAME",17))
	{
	  sscanf(&line_buf[20], "%s %s %c", resname, aname, &level);
	  for(j=0;j<NRES;j++)
	    {
	      if(!strncmp(RESIDUE_LABEL[j],resname,strlen(resname)))
		for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
		  if(!strncmp(&atom_name[i][0],aname,strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))
		    {
		      nrest++; 
		      at_level[i]=level;
		      for(k=0;k<h_num[i];k++)
			at_level[h_list[i][k]]=level;
		    }
	    }
	  printf("Model level  %c set for %6i atoms: rname %s/%s\n", level, nrest, resname,aname);
	  continue;
	}
    /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"LEVEL RESID:ANAME",17))
	{
	  sscanf(&line_buf[20], "%i%i%s %c", &start, &stop, aname, &level);
	  for(i=0,nrest=0;i<N;i++)
	    if((resSeq[i]>=start-1)&&(resSeq[i]<=stop-1))
	      {
		if(!strncmp(&atom_name[i][0],aname,strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))
		  {
		    nrest++; 
		    at_level[i]=level;
		      for(k=0;k<h_num[i];k++)
			at_level[h_list[i][k]]=level;
		  }
	      }
      printf("Model level  %c set for %6i atoms: resid %i-%i/%s\n",level,nrest,start,stop, aname);
      continue;
	}

      /*-------------------------LINKS------------------------------------*/

      if(!strncmp(line_buf,"LINK RNAME:ANAME",16))
	{
	  sscanf(&line_buf[20], "%s %s %s %s %s", resname, aname, lAtm, lTyp, bondedTo);
	  for(j=0;j<NRES;j++)
	    {
	      if(!strncmp(RESIDUE_LABEL[j],resname,strlen(resname)))
		{
		for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
		  if(!strncmp(&atom_name[i][0], bondedTo, strlen(bondedTo))&&(strlen(bondedTo)==strlen(&atom_name[i][0])))
		    bnd2=i+1;
		for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
		    if(!strncmp(&atom_name[i][0], aname, strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))
		      {
			nlink++; 
			link[i]=bnd2;
			strcpy(&linkAtom[i][0],lAtm);
			strcpy(&linkType[i][0],lTyp);
		      }		
		}  
	    }
	  printf("Link %s %s-%s:%i\n", resname, aname, bondedTo,  nlink);continue;
	}
      /*--------------------------------------------------------------------*/
     if(!strncmp(line_buf,"LINK RESID:ANAME",16))
	{
	  sscanf(&line_buf[20], "%i %i %s %s %s %s", &start, &stop, aname, lAtm, lTyp, bondedTo);
	  for(j=0;j<NRES;j++)
	    {
	      if((j>=start-1)&&(j<=stop-1))
		{
		for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
		  if(!strncmp(&atom_name[i][0], bondedTo, strlen(bondedTo))&&(strlen(bondedTo)==strlen(&atom_name[i][0])))
		    bnd2=i+1;
		for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
		    if(!strncmp(&atom_name[i][0], aname, strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))
		      {
			nlink++; 
			link[i]=bnd2;
			strcpy(&linkAtom[i][0],lAtm);
			strcpy(&linkType[i][0],lTyp);
		      }		
		}  
	    }
	  printf("Link resid %i-%i %s-%s:%i\n", start,stop, aname, bondedTo,  nlink);continue;
	}
     /*--------------------------------------------------------------------*/
     if(!strncmp(line_buf,"LINK RESID:INDEX",16))
	{
	  sscanf(&line_buf[20], "%i %s %s %s %i", &start,  aname, lAtm, lTyp, &bnd2);
	  j=start-1;
	  for(i=RESIDUE_POINTER[j]-1;i<RESIDUE_POINTER[j+1]-1;i++)
	    if(!strncmp(&atom_name[i][0], aname, strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))
	      {
		nlink++; 
		link[i]=bnd2;
		strcpy(&linkAtom[i][0],lAtm);
		strcpy(&linkType[i][0],lTyp);
	      }		
	  printf("Link resid %i %s-%i:%i\n", start, aname, bnd2,  nlink);continue;
	}

      /*---------------------------FREESES----------------------------------*/

      if(!strncmp(line_buf,"FREEZE ALL",10))
	{
	  sscanf(&line_buf[20], "%i", &cf);
	  for(i=0,nrest=0;i<N;i++)
	    {nrest++; at_freeze[i]=cf;}
	  printf("Freeze type %2i set for %6i atoms: (all)\n",cf,nrest);continue;
	}
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"FREEZE INDEX ",13))
	{
	  sscanf(&line_buf[20], "%i%i%i", &start, &stop, &cf);
	  for(i=0,nrest=0;i<N;i++)
	    if((i>=start-1)&&(i<=stop-1))
	      {nrest++; at_freeze[i]=cf;}
	  printf("Freeze type %2i set for %6i atoms: index %i-%i\n", cf, nrest, start,stop);continue;
	}
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"FREEZE INDEXH",13))
	{
	  sscanf(&line_buf[20], "%i%i%i", &start, &stop, &cf);
	  for(i=0,nrest=0;i<N;i++)
	    if((i>=start-1)&&(i<=stop-1))
	      {
		nrest++; at_freeze[i]=cf;
		for(k=0;k<h_num[i];k++)
		  at_freeze[h_list[i][k]]=cf;
	      }
	  printf("Freeze type %2i set for %6i atoms: index %i-%i\n", cf, nrest, start,stop);continue;
	}
    /*-------------------------------------------------------------*/
     if(!strncmp(line_buf,"FREEZE RESID ",13))
	{
	  sscanf(&line_buf[20], "%i %i %i", &start, &stop, &cf);
	  for(i=0,nrest=0;i<N;i++)
	    if((resSeq[i]>=start-1)&&(resSeq[i]<=stop-1))
	      {nrest++; at_freeze[i]=cf;}
	  printf("Freeze type %2i set for %6i atoms: resid %i-%i\n",cf,nrest,start,stop);continue;
	}
      /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"FREEZE ANAME",12))
	{
	  sscanf(&line_buf[20], "%s %i",aname, &cf);
	  for(i=0,nrest=0;i<N;i++)
	    if(!strncmp(&atom_name[i][0], aname, strlen(aname)))
	      {nrest++; at_freeze[i]=cf;}
	  printf("Freeze type %2i set for %6i atoms: aname %s*\n", cf, nrest, aname);continue;
	}
    /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"FREEZE RNAME ",13))
	{
	  sscanf(&line_buf[20], "%s%i", resname, &cf);
	  for(i=0,nrest=0;i<N;i++)
	    if(!strncmp(resName[i], resname, strlen(resname)))
	      {nrest++; at_freeze[i]=cf;}
	  printf("Freeze type %2i set for %6i atoms: rname %s\n", cf, nrest, resname);continue;
	}
    /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"FREEZE RNAME:ANAME",18))
	{
	  sscanf(&line_buf[20], "%s %s %i", resname, aname, &cf);
	  for(j=0;j<NRES;j++)
	    {
	      if(!strncmp(RESIDUE_LABEL[j],resname,strlen(resname)))
		for(i=RESIDUE_POINTER[j];i<=RESIDUE_POINTER[j+1]-1;i++)
		  if(!strncmp(&atom_name[i][0],aname,strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))	
		    {
		      nrest++; 
		      at_freeze[i]=cf;
		      for(k=0;k<h_num[i];k++)
			at_freeze[h_list[i][k]]=cf;
		    }
	    }
	  printf("Freeze type  %2i set for %6i atoms: rname %s\n", cf, nrest, resname);
	  continue;
	}
    /*-------------------------------------------------------------*/
      if(!strncmp(line_buf,"FREEZE RESID:ANAME",18))
	{
	  sscanf(&line_buf[20], "%i%i%s%i", &start, &stop, aname, &cf);
	  for(i=0,nrest=0;i<N;i++)
	    if((resSeq[i]>=start-1)&&(resSeq[i]<=stop-1))
	      {
		if(!strncmp(&atom_name[i][0],aname,strlen(aname))&&(strlen(aname)==strlen(&atom_name[i][0])))
		  {
		    nrest++; 
		    at_freeze[i]=cf;
		      for(k=0;k<h_num[i];k++)
			at_freeze[h_list[i][k]]=cf;
		  }
	      }
      printf("Freeze type %2i set for %6i atoms: resid %i-%i/%s\n",cf,nrest,start,stop, aname);
      continue;
	}
  
  if(strlen(line_buf)>5)
    printf("Unknown keyword: %s",line_buf);
}
fclose(fp);
}



void write_oniom(char *restr, char *label)
{
  long int i,j;
  int D_chrg;
  FILE *fp;
  float chrg, chrgf;
  int *at_freeze;
  char *at_level;
  int *link;
  int bnd2;
  char tmpstr[100];

  
  chrg=0;
  chrgf=0;
  
  at_freeze=malloc(N*sizeof(int));
  at_level=malloc(N*sizeof(char));
  link=malloc(N*sizeof(int));

  for(i=0;i<N;i++){
    at_freeze[i]=0.0;
    at_level[i]='L';
  }

  mass_to_element();
  read_freeze(restr, at_freeze, at_level, link);
 

  // MAKE LISTS OF BONDED ATOMS
  for(i=0;i<N;i++)
    {// read only 3 digits in charges  
      chrg+=nearbyint(charge[i]*1000)/1000;
      chrgf+=charge[i];
      connect_num[i]=0;
      // BONDS WITHOUT HYDROGEN
      for(j=0;j<NBONA;j++)
	{
	  if(IB[j]==i)
	    connect_list[i][connect_num[i]++]=JB[j];  
	  if(JB[j]==i)
	    connect_list[i][connect_num[i]++]=IB[j];
	}
      // BONDS WITH HYDROGEN     
      for(j=0;j<NBONH;j++)
	{
	  if(IBH[j]==i)
	    connect_list[i][connect_num[i]++]=JBH[j];
	  if(JBH[j]==i)
	    connect_list[i][connect_num[i]++]=IBH[j];
	}  
    }

  D_chrg = nearbyint((chrg-chrgf)*1000);

  if(nearbyint(chrgf*1000)/1000 != 0.0)
    printf("\n ** WARNING! ** Total charge is %.4f\n Check if the charge was intended to be non zero\n\n", chrgf);
  // CORRECT CHARGE
  
  printf("Total/Trunc charge: %.4f/%.3f, Diff: 0.00%i\n",chrgf, chrg, abs(D_chrg) );
  printf("Changing charges on %i atoms by 0.001\n", abs(D_chrg) );
  for(i=0;i<abs(D_chrg);i++)
    if(D_chrg>0)
      charge[N-1-i]-=0.001;
    else
      charge[N-1-i]+=0.001;
  chrg=0;
  for(i=0;i<N;i++)
    chrg+=nearbyint(charge[i]*1000)/1000;
  printf("Corrected charge:   %.4f\n", chrg);

 
  fp = fopen(label,"wt");
  if(fp == NULL){printf("** E ** could not open file %s for writing, abort\n", label); exit(1);}  
  printf("Writing ONIOM input file: %s\n", label);
  
  //Strip spaces
  for(i=0;i<N;i++)
      if(Type[i][1]==' ')Type[i][1]='\0';
  fprintf(fp,"%cmem=320MW\n",'%');
  fprintf(fp,"%cnproc=8\n",'%');
  fprintf(fp,"#T ONIOM(HF/6-31G:Amber=softonly)=EmbedCharge Geom=Connectivity Opt\n\nTitle Card\n\n0 1 0 1 0 1\n");
  for(i=0;i<N;i++)
    {
      sprintf(tmpstr,"%s-%s-%.3f",Element[i],Type[i],nearbyint(charge[i]*1000)/1000);
      fprintf(fp,"%s",tmpstr);   
      for(j=0;j<12-strlen(tmpstr);j++)
	fprintf(fp," ");  
      fprintf(fp,"%3i%11.6f%11.6f%11.6f %c ",at_freeze[i],X[i],Y[i],Z[i],at_level[i]);
      if(link[i])
	fprintf(fp,"%s-%s-0.1 %i\n",&linkAtom[i][0], &linkType[i][0],link[i]);
      else
	fprintf(fp,"\n");
    }

  // PRINT CONNECTIVITY
  fprintf(fp,"\n");
  for(i=0;i<N;i++)
    {
      fprintf(fp,"%i ",i+1);
      
      for(j=0;j<connect_num[i];j++)
	if(connect_list[i][j]>i)
	  fprintf(fp,"%i 1.0 ",connect_list[i][j]+1);		
      fprintf(fp,"\n");
    }
  fprintf(fp,"\n");
  for(i=0;i<prmcnt;i++)
  fprintf(fp,"@%s\n",params[i]);
  fprintf(fp,"\n\n");
  fclose(fp);

  //Save indexed file for debugging
  fp=fopen("debug.com","wt");
  for(i=0;i<N;i++)
    {
      chainID[i][0]=at_level[i];
      tempFactor[i]=at_freeze[i];
      sprintf(tmpstr,"%5i %s-%s-%.3f",i,Element[i],Type[i],nearbyint(charge[i]*1000)/1000);
      fprintf(fp,"%s",tmpstr);   
      for(j=0;j<18-strlen(tmpstr);j++)
	fprintf(fp," ");  
      fprintf(fp,"%3i%11.6f%11.6f%11.6f %c ",at_freeze[i],X[i],Y[i],Z[i],at_level[i]);
      if(link[i])
	{fprintf(fp,"%s-%s-0.1 %i\n",&linkAtom[i][0], &linkType[i][0],link[i]); chainID[i][0]='Q';}
      else
	fprintf(fp,"\n");
    }
  fclose(fp); 

 for(i=0;i<N;i++)
   if(link[i])
     chainID[link[i]-1][0]='q';
     
  free(at_freeze);
  free(at_level);
  free(link);
}



int main(int argc, char **argv)
{  
  int c;
  char usage[]="\n"

    " -----------------------------------------------------------------\n"
    " * This utility prepares input file for g09 oniom calculations   *\n"
    " -----------------------------------------------------------------\n"
    "\n  Command line options:\n"
    "  -c inpcrd      input  coordinates\n"
    "  -o oniom       output ONIOM in\n"
    "  -p prmtop      input AMBER7  parameters\n"
    "  -r const.in    input frozen atoms\n"
    "  -h             prints this message\n\n"
    " Keywords:\n"
    "LEVEL ALL           L|M|H\n"
    "LEVEL RESID         1 2 L|M,H\n"
    "LEVEL INDEX         1 2 L|M,H\n"
    "LEVEL INDEXH        1 2 L|M,H\n"
    "LEVEL ANAME         MG  L|M,H\n"
    "LEVEL RNAME         BCL L|M,H\n"
    "LEVEL RNAME:ATOM    BCL 1 4 L|M,H\n"
    "LEVEL RNAME:ANAME   BCL CA L|M,H\n"
    "LEVEL RESID:ANAME   1 2 CA L|M,H\n"
    "FREEZE ALL          -1|0\n"
    "FREEZE INDEX        1 2 -1|0\n"
    "FREEZE INDEXH       1 2 -1|0\n"
    "FREEZE RESID        1 2 -1|0\n" 
    "FREEZE ANAME        MG -1|0\n"
    "FREEZE RNAME        BCL -1|0\n"
    "FREEZE RNAME:ANAME  BCL CA -1|0\n" 
    "FREEZE RESID:ANAME  1 2 CA -1|0\n" 
    "LINK RNAME:ANAME    U10 C7 H HC C8\n"
    "LINK RESID:ANAME    1 10 CB H HC CA\n"
    "LINK RESID:INDEX    1 CB H HC 1221\n"
    "PARAMS /path/file\n"
    "DO NOT FORGET TO ADD 1 TO VMD INDEX!!!\n";

  char *prmtop = "prmtop";
  char *coord = "inpcrd";
  char *oniom = "oniom.com";
  char *restr = "const.in";

  opterr=0;
  while( (c=getopt(argc,argv,"p:c:r:o:h")) != -1) 
    {
      switch(c) 
	{
	case 'o': 
	  oniom=optarg;
	  break;
	case 'p': 
	  prmtop=optarg;
	  break;
	case 'c': 
	  coord=optarg;
	  break;
	case 'r': 
	  restr=optarg;
	  break;
	case 'h':
	  printf("%s",usage);
	  abort();
	}
    }
  
  if(!read_parm7(prmtop))
    exit(0);

  read_amber_coor(coord, N, X, Y, Z);
  write_oniom(restr,oniom);
  read_parm7(prmtop);
  write_one_pdb("oniom.pdb");

  return(1);
}

