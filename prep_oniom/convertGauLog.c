#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../File_Utilities/file_read_write.h"
#include "../File_Utilities/file_read_write.c"


int main(int argc, char **argv)
{  
  int c;
  char *gaulog="oniom.log";

  char usage[]="\n"

    " -------------------------------------------------------------------------\n"
    " * This utility extraxts coordinates from gaussian log and saves as rst7  *\n"
    " --------------------------------------------------------------------------\n"
    "\n  Command line options:\n"
    "  -l oniom.log    input  gaussian log\n";


  opterr=0;
  while( (c=getopt(argc,argv,"l:h")) != -1) 
    {
      switch(c) 
	{
	case 'l': 
	  gaulog=optarg;
	  break;

	case 'h':
	  printf("%s",usage);
	  abort();
	}
    }
  
  read_gaulog(gaulog);
}
