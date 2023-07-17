#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#ifndef _HSI_Index_H_
#define _HSI_Index_H_

#define HSI_STD_BENCHMARK_DIST 100000

// The struct where we store the metadata related to a single sequence
typedef struct _seq_metadata_t_ {
  char * SeqName;
  uint32_t * Benchmarks;
  uint32_t num_benchmarks;
  uint32_t seq_length;
  uint64_t seq_start_byte;
  char alphabet; // (D)NA / (R)NA / (P)rotein <-[anything that isn't DNA or RNA]
} SEQ_META;


// The struct where we store the metadata related to a whole fasta file
typedef struct _file_metadata_t_ {
  char * fname;
  char * index_fname;
  //char md5[33];
  uint64_t filesize;
  uint32_t num_seqs;
  uint32_t max_seq_name_len;
  uint32_t benchmark_dist;
  SEQ_META ** Seqs;
  char Alphabets[3]; // 0=DNA / 1=RNA / 2=Protein
  char set_alpha;
} FILE_META;



/*
 *  GENERIC FUNCTIONS (not specific to fasta indexing)
 */
int  IsLexLessThan (char * A, char * B);
void LexSort (FILE_META * FM, uint32_t * Index);


/*
 *  INDEX CONSTRUCTION FUNCTIONS
 */
void GetIndexFname (char * fname, FILE_META * FM);
void GetFileSize (char * fname, uint64_t * filesize);
void GetMD5 (char * fname, char * md5);
int  IndexAlreadyExists(FILE_META * FM);
void ResizeFMSeqs (FILE_META * FM, uint32_t * max_seqs);
void ResizeBenchmarks (SEQ_META * S, uint32_t * max_benchmarks);
char DetermineAlphabet (char * Residues);
void GetSeqMetadata (FILE_META * FM);
int  GetFileMetadata (char * infname, FILE_META * FM, int force);
void WriteIndex (FILE_META * FM);
void BuildSeqIndex (char * infname, uint32_t benchmark_dist, char set_alpha, int force);


/*
 *  INDEX PARSING FUNCTIONS
 */
void ParseHSIMetadata (FILE_META * FM);


#endif
