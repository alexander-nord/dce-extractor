#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hsi_index.h"


// Function List
// -------------
// + PrintUsage
// + FindStartByte
// + ReadDNARevcomp
// + ReadRNARevcomp
// + ReadProteinRev
// + ExtractSequenceRange
// + main



///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  PrintUsage
//
//  About: Provide some basic usage instructions
//
int PrintUsage () {
  fprintf(stderr,"\n");
  fprintf(stderr,"  USAGE: ./sfetch {OPT.S} [file] [seqname]\n");
  fprintf(stderr,"  OPT.S: -range [start]..[end]\n");
  fprintf(stderr,"\n");
  return 1;
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  FindStartByte
//
//  About: Identify where in the index the first byte related to our sequence
//         of interest resides.
//
uint64_t FindStartByte (FILE_META * FM, char * SeqName) {

  // First off, what's the starting byte of our list of sequence names?
  // NOTE: The reason for the numbers is given in more detail in 'seq_index.c'
  long int list_start = 112;
  long int list_entry_len = (long int)FM->max_seq_name_len + 31; // Includes '\n'

  // Open up the index file
  FILE * inf = fopen(FM->index_fname,"r");

  // Where we'll store the next sequence name that we compare against
  char CompName[FM->max_seq_name_len];
  
  // We'll be binary-searching our way through the names
  long int search_low  = 0;
  long int search_high = FM->num_seqs;

  // Until we find our line, keep cruising along...
  while (search_low <= search_high) {

    char last_chance = 0;
    if (search_low == search_high)
      last_chance = 1;

    long int search = (search_low + search_high) / 2;
    long int pos    = list_start + list_entry_len * search;
    
    fseek(inf,pos,0);
    fscanf(inf,": %s",CompName);

    if (!strcmp(CompName,SeqName) || last_chance)
      break;

    // Do we need to look higher or lower?
    if (IsLexLessThan(CompName,SeqName))
      search_low = search;
    else
      search_high = search;
    
  }

  // If we've exhausted our list, then it seems we're on the wildest goose chase
  if (strcmp(CompName,SeqName)) {
    fclose(inf);
    fprintf(stderr,"  Error:  Failed to locate sequence '%s' in file\n",SeqName);
    exit(1);
  }

  char alphabet;
  char seq_len_hex[9];
  char start_byte_hex[17];
  fscanf(inf,"%s %c %s\n",seq_len_hex,&alphabet,start_byte_hex);
  fclose(inf);

  return strtoll(start_byte_hex,NULL,16);
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  ReadDNARevcomp
//
//  About: Read in a DNA sequence in a revcomp-y way
//
void ReadDNARevcomp (FILE * inf, char * Seq, uint32_t len) {
  while (len) {
    char c = (char)getc(inf);
    if (c > 32) {
      switch (c) {
      case 'A':
	Seq[--len] = 'T';
	break;
      case 'C':
	Seq[--len] = 'G';
	break;
      case 'G':
	Seq[--len] = 'C';
	break;
      case 'T':
	Seq[--len] = 'A';
	break;
      case 'a':
	Seq[--len] = 't';
	break;
      case 'c':
	Seq[--len] = 'g';
	break;
      case 'g':
	Seq[--len] = 'c';
	break;
      case 't':
	Seq[--len] = 'a';
	break;
      default:
	Seq[--len] = 'N';
      }
    }
  }
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  ReadRNARevcomp
//
//  About: Read in an RNA sequence in a revcomp-y way
//
void ReadRNARevcomp (FILE * inf, char * Seq, uint32_t len) {
  while (len) {
    char c = (char)getc(inf);
    if (c > 32) {
      switch (c) {
      case 'A':
	Seq[--len] = 'U';
	break;
      case 'C':
	Seq[--len] = 'G';
	break;
      case 'G':
	Seq[--len] = 'C';
	break;
      case 'U':
	Seq[--len] = 'A';
	break;
      case 'a':
	Seq[--len] = 'u';
	break;
      case 'c':
	Seq[--len] = 'g';
	break;
      case 'g':
	Seq[--len] = 'c';
	break;
      case 'u':
	Seq[--len] = 'a';
	break;
      default:
	Seq[--len] = 'N';
      }
    }
  }
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  ReadProteinRev
//
//  About: Read in a protein sequence in a reverse-y way
//
void ReadProteinRev (FILE * inf, char * Seq, uint32_t len) {
  while (len) {
    char c = (char)getc(inf);
    if (c > 32)
      Seq[--len] = c;
  }
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  ExtractSequenceRange
//
//  About: Pull a range of sequence from a fasta file
//
void ExtractSequenceRange (FILE_META * FM, char * SeqName, uint32_t start, uint32_t end) {

  // Where does the sequence of interest reside in the index file?
  uint64_t index_start = FindStartByte(FM,SeqName);

  char name[FM->max_seq_name_len+1];
  char smallhex[9];
  char bighex[17];

  // Open the index and jump to our starting position!
  FILE * inf = fopen(FM->index_fname,"r");
  fseek(inf,index_start,0);

  // Sanity check -- do we have a matching name?
  fscanf(inf,"SEQ NAME: %s\n",name);
  if (strcmp(name,SeqName)) {
    fclose(inf);
    fprintf(stderr,"  Error: Search landed in wrong position in index file!\n");
    exit(1);
  }

  // Read in the sequence length
  fscanf(inf,"SEQ LENGTH: %s\n",smallhex);
  uint32_t seq_len = strtol(smallhex,NULL,16);

  // If the user wants the full sequence, we'll need to oblige them
  uint8_t fetching_full_seq = 0;
  if (start == 0 && end == 0) {
    fetching_full_seq = 1;
    start = 1;
    end = seq_len;
  }

  // Let's also switch over to [0..seq_len-1]
  start--;
  end--;
  
  // Now we can actually assume we're gonna be checking out the sequence!
  // We'll want start to always be lower than end, so here's a way to make that happen
  int revcomp = 0;
  if (start > end) {
    uint32_t temp = start;
    start = end;
    end = temp;
    revcomp = 1;
  }
  
  // Sanity check -- is this a valid range?
  if (end >= seq_len) {
    fclose(inf);
    fprintf(stderr,"  Error: Requested end coordinate is larger than sequence length");
    fprintf(stderr," (%u)\n",seq_len);
    exit(1);
  }

  // Read in the alphabet code
  char alphabet;
  fscanf(inf,"SEQ ALPHABET: %c\n",&alphabet);

  // What's the start byte in the actual sequence file?
  fscanf(inf,"SEQ START BYTE: %s\n",bighex);
  uint64_t seq_start_byte = strtoll(bighex,NULL,16);

  // We also record the number of benchmarks for each sequence,
  // so don't forget to scan through that bit!
  fseek(inf,29,1);

  // How far do we need to scan from our current position in inf to get to our
  // benchmark? NOTE that this will be one too high, so if we're at zero that
  // means we don't go to the next benchmark!
  if (start >= FM->benchmark_dist) {
    uint32_t benchmark = 9 * (start / FM->benchmark_dist) - 9;
    fseek(inf,benchmark,1);
    fscanf(inf,"%s\n",smallhex);
    seq_start_byte += (uint64_t)strtol(smallhex,NULL,16);
  }

  // No need for the index anymore!
  fclose(inf);

  // Now that we know where we're jumping to in our sequence file (in terms of bytes),
  // how many characters do we need to gobble up before we're at the requested range?
  uint32_t offset = start - FM->benchmark_dist * (start / FM->benchmark_dist);

  // OOOOHHH, it's scannin' time!
  inf = fopen(FM->fname,"r");
  fseek(inf,seq_start_byte,0);
  char c;
  while (offset) {
    c = (char)fgetc(inf);
    if (c > 32)
      offset--;
  }

  printf(">\%s",SeqName);
  if (!fetching_full_seq) {
    if (revcomp) printf("/%d-%d",end+1,start+1);
    else         printf("/%d-%d",start+1,end+1);
  }
  printf("\n");
  
  // AWWWW YEAH! Now, we READ!
  int line_len = 0;
  int new_line = 60;
  if (revcomp) {

    uint32_t seq_len = end-start+1;
    char * Seq = malloc(seq_len * sizeof(char));
    
    if      (alphabet == 'D') ReadDNARevcomp(inf,Seq,seq_len);
    else if (alphabet == 'R') ReadRNARevcomp(inf,Seq,seq_len);
    else               /* P */ReadProteinRev(inf,Seq,seq_len);

    uint32_t pos = 0;
    while (pos < seq_len) {
      printf("%c",Seq[pos++]);
      line_len++;
      if (line_len == new_line) {
	printf("\n");
	line_len = 0;
      }
    }

    free(Seq);

  } else {
    
    while (start <= end) {
      c = (char)getc(inf);
      if (c > 32) {
	printf("%c",c);
	start++;
	line_len++;
	if (line_len == new_line) {
	  printf("\n");
	  line_len = 0;
	}
      }
    }
    
  }
  
  if (line_len) printf("\n");

  // That's it, baby!
  fclose(inf);

}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  main
// 
int main (int argc, char ** argv) {

  if (argc < 3) return PrintUsage();

  int i=0;
      
  // Are we being asked to pull a specific range?
  uint32_t start = 0;
  uint32_t end   = 0;

  // Are we being asked to index a file?
  int indexing = 0;

  // Are we being asked to index a file *and* index it assuming a particular alphabet?
  char set_alpha = 0;

  int arg_parser = 1;
  int end_arg_parse = argc-2; // Assume we're pulling a seq, until we see order to index
  while (arg_parser < end_arg_parse) {

    if (!strcmp(argv[arg_parser],"-range")) {

      // Requesting a particular range
      arg_parser++;
      while (argv[arg_parser][i] != '.') {
	start = start * 10 + (uint32_t)(argv[arg_parser][i] - '0');
	i++;
      }
      i += 2;
      while (i < strlen(argv[arg_parser])) {
	end  = end * 10 + (uint32_t)(argv[arg_parser][i] - '0');
	i++;
      }

    } else if (!strcmp(argv[arg_parser],"--index")) {

      // Indicating that we're indexing
      indexing = 1;
      end_arg_parse++;

    } else if (!strcmp(argv[arg_parser],"-alphabet")) {

      // Specifying an alphabet
      set_alpha = argv[++arg_parser][0];
      if (set_alpha > 96) set_alpha -= 32;
      if (set_alpha != 'D' && set_alpha != 'R' && set_alpha != 'P') {
	fprintf(stderr,"  Error: Unrecognized alphabet '%s'\n",argv[arg_parser]);
	return 1;
      }
      
    } else {
      fprintf(stderr,"  Warning: Unrecognized option '%s' ignored\n",argv[arg_parser]);
    }
    
    arg_parser++;
    
  }

  // We'll need a file_meta struct to store information about the seq file (duh)
  FILE_META * FM = malloc(sizeof(FILE_META));

  // If this file hasn't been indexed, build an index right quick
  GetIndexFname(argv[argc-2],FM);
  FILE * inf = fopen(FM->index_fname,"r");
  if (inf == NULL) {
    fprintf(stderr,"  Creating index file '%s' (this may take a minute)\n",FM->index_fname);
    BuildSeqIndex(FM->fname,HSI_STD_BENCHMARK_DIST,0,1);
  } else {
    fclose(inf);
  }

  // Pull in the index's metadata
  ParseHSIMetadata(FM);
  
  // Extraction time!
  ExtractSequenceRange(FM,argv[argc-1],start,end);

  return 0;

}
