#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hsi_index.h"


// Function List
// -------------
// + PrintUsage
// + main



///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  PrintUsage
//
int PrintUsage () {
  fprintf(stderr,"\n");
  fprintf(stderr,"  USAGE: ./sindex {OPT.S} [file]\n");
  fprintf(stderr,"  OPT.S: -alphabet [DNA/RNA/Protein]\n");
  fprintf(stderr,"         --force\n");
  fprintf(stderr,"\n");
  return 1;
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  main
//
int main (int argc, char ** argv) {

  if (argc < 2) return PrintUsage();
  
  char set_alpha = 0;
  int force = 0;

  int arg_parser = 1;
  while (arg_parser < argc-1) {
  
    if (!strcmp(argv[arg_parser],"-alphabet")) {

      // Specifying an alphabet
      set_alpha = argv[++arg_parser][0];
      if (set_alpha > 96) set_alpha -= 32;
      if (set_alpha != 'D' && set_alpha != 'R' && set_alpha != 'P') {
        fprintf(stderr,"  Error: Unrecognized alphabet '%s'\n",argv[arg_parser]);
        return 1;
      }

    } else if (!strcmp(argv[arg_parser],"--force")) {

      // Force index construction, even if one exists
      force = 1;
      
    } else {
      fprintf(stderr,"  Warning: Unrecognized option '%s' ignored\n",argv[arg_parser]);
    }

    arg_parser++;

  }

  BuildSeqIndex(argv[argc-1],HSI_STD_BENCHMARK_DIST,set_alpha,force);
  return 0;

}
