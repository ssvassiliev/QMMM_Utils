float *bkbn;
int *ibkbn;

void define_bkbn(void)
{
  int i,j,k;
  double bkbn_crg, crgtmp;
  for(i=0,j=0;i<NRES;i++)
    {
      if(!strncmp(RESIDUE_LABEL[i],"ALA",3)||
	 !strncmp(RESIDUE_LABEL[i],"ARG",3)||
	 !strncmp(RESIDUE_LABEL[i],"ARN",3)||
	 !strncmp(RESIDUE_LABEL[i],"ASH",3)||
	 !strncmp(RESIDUE_LABEL[i],"ASP",3)||
	 !strncmp(RESIDUE_LABEL[i],"ASN",3)||
	 !strncmp(RESIDUE_LABEL[i],"CYM",3)||
	 !strncmp(RESIDUE_LABEL[i],"CYS",3)||
	 !strncmp(RESIDUE_LABEL[i],"CYX",3)||
	 !strncmp(RESIDUE_LABEL[i],"GLH",3)||
	 !strncmp(RESIDUE_LABEL[i],"GLN",3)||
	 !strncmp(RESIDUE_LABEL[i],"GLU",3)||
	 !strncmp(RESIDUE_LABEL[i],"HIE",3)||
	 !strncmp(RESIDUE_LABEL[i],"HID",3)||
	 !strncmp(RESIDUE_LABEL[i],"HIP",3)||
	 !strncmp(RESIDUE_LABEL[i],"ILE",3)||
	 !strncmp(RESIDUE_LABEL[i],"LEU",3)||
	 !strncmp(RESIDUE_LABEL[i],"LYN",3)||
	 !strncmp(RESIDUE_LABEL[i],"LYS",3)||
	 !strncmp(RESIDUE_LABEL[i],"MET",3)||
	 !strncmp(RESIDUE_LABEL[i],"PHE",3)||
	 !strncmp(RESIDUE_LABEL[i],"SER",3)||
	 !strncmp(RESIDUE_LABEL[i],"THR",3)||
	 !strncmp(RESIDUE_LABEL[i],"TRP",3)||
	 !strncmp(RESIDUE_LABEL[i],"TYR",3)||
	 !strncmp(RESIDUE_LABEL[i],"VAL",3))
	{
	  k=-1;
	  crgtmp=0.0;
	  while(strncmp(atom_name[RESIDUE_POINTER[i]+k],"CB",2))
	    {
	      bkbn[RESIDUE_POINTER[i]+k]=charge[RESIDUE_POINTER[i]+k];
	      ibkbn[RESIDUE_POINTER[i]+k]=1;
	      crgtmp+=bkbn[RESIDUE_POINTER[i]+k];
	      k++;
	      if(k>25)
		{printf("** ERROR ** Residue %i is not protein\n",i);break;}
	    }
	  ibkbn[RESIDUE_POINTER[i]+k]=-1; // CB
	  bkbn[RESIDUE_POINTER[i+1]-3]=charge[RESIDUE_POINTER[i+1]-3];
	  ibkbn[RESIDUE_POINTER[i+1]-3]=1;
	  bkbn[RESIDUE_POINTER[i+1]-2]=charge[RESIDUE_POINTER[i+1]-2];
	  ibkbn[RESIDUE_POINTER[i+1]-2]=1;
	  crgtmp+=bkbn[RESIDUE_POINTER[i+1]-3];
	  crgtmp+=bkbn[RESIDUE_POINTER[i+1]-2];
	  bkbn[RESIDUE_POINTER[i]+k-2]-=crgtmp;
	  j++;
	}

      if(!strncmp(RESIDUE_LABEL[i],"PRO",3))
	{
	  crgtmp=0.0;
	  bkbn[RESIDUE_POINTER[i]-1]=charge[RESIDUE_POINTER[i]-1];
	  bkbn[RESIDUE_POINTER[i+1]-5]=charge[RESIDUE_POINTER[i+1]-5];
	  bkbn[RESIDUE_POINTER[i+1]-4]=charge[RESIDUE_POINTER[i+1]-4];
	  bkbn[RESIDUE_POINTER[i+1]-3]=charge[RESIDUE_POINTER[i+1]-3];
	  bkbn[RESIDUE_POINTER[i+1]-2]=charge[RESIDUE_POINTER[i+1]-2];
	  ibkbn[RESIDUE_POINTER[i]-1]=1;
	  ibkbn[RESIDUE_POINTER[i+1]-5]=1;
	  ibkbn[RESIDUE_POINTER[i+1]-4]=1;
	  ibkbn[RESIDUE_POINTER[i+1]-3]=1;
	  ibkbn[RESIDUE_POINTER[i+1]-2]=1;
	  crgtmp+=bkbn[RESIDUE_POINTER[i]-1];
	  crgtmp+=bkbn[RESIDUE_POINTER[i+1]-5];
	  crgtmp+=bkbn[RESIDUE_POINTER[i+1]-4];
	  crgtmp+=bkbn[RESIDUE_POINTER[i+1]-3];
	  crgtmp+=bkbn[RESIDUE_POINTER[i+1]-2];
	  bkbn[RESIDUE_POINTER[i+1]-5]-=crgtmp;
	  j++; 
	}

      if(!strncmp(RESIDUE_LABEL[i],"GLY",3))
	{
	  for(k=RESIDUE_POINTER[i]-1;k<RESIDUE_POINTER[i+1]-1;k++)
	    {bkbn[k]=charge[k];ibkbn[k]=1;}
	  j++;
	}

    }
  bkbn_crg=0;
  for(i=0;i<N;i++){
    bkbn_crg+=bkbn[i];
    tempFactor[i]=ibkbn[i];
  }

  //  write_one_pdb("test.pdb");
  printf("Backbone: %i residues charge = %f\n",j,bkbn_crg);
}

