#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hsi_index.h"


// Function List
// -------------
// + PrintUsage
// + PrintSeqStats
// + main




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  PrintUsage
//
int PrintUsage () {
  fprintf(stderr,"\n  Usage: ./sstat [file.fa]\n\n");
  return 1;
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  PrintSeqStats
//
void PrintSeqStats (FILE_META * FM) {

  // Parse the metadata from our index file
  ParseHSIMetadata(FM);

  // We'll open up our index and jump ahead to the end of the metadata
  FILE * inf = fopen(FM->index_fname,"r");
  fseek(inf,112,0);

  char smallhex[9];
  char bighex[17];

  // Are there multiple alphabets in use?
  int multiple_alphas = FM->Alphabets[0] + FM->Alphabets[1] + FM->Alphabets[2] - 1;
  
  // Now we'll run through and record all the data we feel like grabbin'!
  FM->Seqs = malloc(FM->num_seqs*sizeof(SEQ_META *));
  uint32_t i,j;
  uint32_t longest_seq = 0;
  uint32_t shortest_seq = 0;
  uint64_t total_residues = 0;
  for (i=0; i<FM->num_seqs; i++) {

    FM->Seqs[i]	= malloc(sizeof(SEQ_META));
    FM->Seqs[i]->SeqName = malloc((FM->max_seq_name_len+1)*sizeof(char));

    fscanf(inf,": %s",FM->Seqs[i]->SeqName);

    char alphabet;
    fscanf(inf,"%s %c %s\n",smallhex,&alphabet,bighex);

    FM->Seqs[i]->seq_length = (uint32_t)strtol(smallhex,NULL,16);
    FM->Seqs[i]->alphabet   = alphabet;

    if (!shortest_seq || FM->Seqs[i]->seq_length < shortest_seq)
      shortest_seq = FM->Seqs[i]->seq_length;
    if (!longest_seq  || FM->Seqs[i]->seq_length > longest_seq)
      longest_seq  = FM->Seqs[i]->seq_length;

    total_residues += (uint64_t)FM->Seqs[i]->seq_length;

  }
  
  fclose(inf);

  // For formatting, we'll need to know how many digits are in the
  // longest sequence's length
  char len_str[10];
  sprintf(len_str,"%u",longest_seq);
  int max_str_len = strlen(len_str);

  char res_str[20];
  sprintf(res_str,"%llu",total_residues);
  int res_str_len = strlen(res_str);
    
  // Now it's just time to report!
  fprintf(stdout,"Number of Sequences : %u\n",FM->num_seqs);
  fprintf(stdout,"Sequence Alphabet(s):");
  if (FM->Alphabets[0]) fprintf(stdout," DNA");
  if (FM->Alphabets[1]) fprintf(stdout," RNA");
  if (FM->Alphabets[2]) fprintf(stdout," Protein");
  fprintf(stdout,"\n");

  fprintf(stdout,"Total Number of Residues: %llu\n",total_residues);

  fprintf(stdout,"Shortest Sequence Length: ");
  sprintf(res_str,"%u",shortest_seq);
  for (i=strlen(res_str); i<res_str_len; i++) fprintf(stdout," ");
  fprintf(stdout,"%u\n",shortest_seq);

  fprintf(stdout,"Longest  Sequence Length: ");
  sprintf(res_str,"%u",longest_seq);
  for (i=strlen(res_str); i<res_str_len; i++) fprintf(stdout," ");
  fprintf(stdout,"%u\n",longest_seq);

  for (i=0; i<FM->num_seqs; i++) {

    fprintf(stdout,": %s ",FM->Seqs[i]->SeqName);

    j=0;
    while (strlen(FM->Seqs[i]->SeqName)+j <= (int)FM->max_seq_name_len) {
      fprintf(stdout," ");
      j++;
    }
    
    // If there are multiple alphabet types in this file, let them know what
    // this sequence's type is.
    if (multiple_alphas) {
      if (FM->Alphabets[2]) {
	if      (FM->Seqs[i]->alphabet == 'D') fprintf(stdout,"DNA     ");
	else if (FM->Seqs[i]->alphabet == 'R') fprintf(stdout,"RNA     ");
	else                            /* P */fprintf(stdout,"Protein ");
      } else {
	if      (FM->Seqs[i]->alphabet == 'D') fprintf(stdout,"DNA ");
	else                            /* R */fprintf(stdout,"RNA ");
      }
    }

    sprintf(len_str,"%u",FM->Seqs[i]->seq_length);
    while (strlen(FM->Seqs[i]->SeqName)+j+strlen(len_str) <= (int)FM->max_seq_name_len + max_str_len) {
      fprintf(stdout," ");
      j++;
    }
    fprintf(stdout,"%s\n",len_str);
    
  }

}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  main
//
int main (int argc, char ** argv) {

  if (argc != 2) return PrintUsage();

  FILE_META * FM = malloc(sizeof(FILE_META));
  
  // If we don't have an index on this file, build one
  GetIndexFname(argv[1],FM);
  FILE * inf = fopen(FM->index_fname,"r");
  if (inf == NULL)
    BuildSeqIndex(FM->fname,HSI_STD_BENCHMARK_DIST,0,1);
  else
    fclose(inf);

  PrintSeqStats(FM);

  // Clear out our file_meta struct
  uint32_t i;
  for (i=0; i<FM->num_seqs; i++) {
    free(FM->Seqs[i]->SeqName);
    free(FM->Seqs[i]);
  }
  free(FM);

  return 0;

}
