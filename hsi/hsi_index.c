#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "hsi_index.h"


// Function List
// -------------
//
// + IsLexLessThan
// + LexSort
// 
// + GetIndexFname
// + GetMD5
// + IndexAlreadyExists
// + GetFileSize
// + ResizeFMSeqs
// + ResizeBenchmarks
// + DetermineAlphabet
// + GetSeqMetadata
// + GetFileMetadata
// + WriteIndex
// + BuildSeqIndex
// 
// + ParseHSIMetadata




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  IsLexLessThan
//
//  About: Is string A lexicographically less than string B?
//
int IsLexLessThan (char * A, char * B) {

  int A_len = strlen(A);
  int B_len = strlen(B);

  int min_len = A_len;
  if (B_len < A_len)
    min_len = B_len;

  int i=0;
  while (i<min_len && A[i] == B[i])
    i++;

  // If we've hit the end of a string, then that string wins!
  if (i==min_len) {
    if (A_len==min_len)
      return 1;
    return 0;
  }

  // Quick comparison and we're solid!
  if (A[i] < B[i])
    return 1;
  return 0;
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  LexSort
//
//  About: Produce a lexicographic sorting of sequence names
//
void LexSort (FILE_META * FM, uint32_t * Index) {

  uint32_t i,j,k,j_lim,k_lim;

  uint32_t * Read  = malloc(FM->num_seqs*sizeof(uint32_t));
  uint32_t * Write = malloc(FM->num_seqs*sizeof(uint32_t));
  uint32_t * Swap;
  
  // Initialize the index
  for (i=0; i<FM->num_seqs; i++)
    Read[i] = i;

  // Merge sort
  uint32_t blocksize  = 1;
  uint32_t num_blocks = FM->num_seqs / blocksize;
  uint32_t pos        = 0;
  while (num_blocks) {

    for (i=0; i+blocksize<FM->num_seqs; i+=2*blocksize) {

      j     = i;
      j_lim = j+blocksize;
      k     = j_lim;
      k_lim = k+blocksize;
      if (k_lim>FM->num_seqs)
	k_lim = FM->num_seqs;
      
      while (j<j_lim && k<k_lim) {
	if (IsLexLessThan(FM->Seqs[Read[j]]->SeqName,FM->Seqs[Read[k]]->SeqName))
	  Write[pos++] = Read[j++];
	else
	  Write[pos++] = Read[k++];
      }

      while (j<j_lim)
	Write[pos++] = Read[j++];

      while (k<k_lim)
	Write[pos++] = Read[k++];
      
    }

    // If we didn't have enough to fill up a full block, we'll still need to
    // copy over the final values into 'Write'
    while (i<FM->num_seqs)
      Write[pos++] = Read[i++];

    // Read becomes write, write becomes read, the first shall be last, etc!
    Swap  = Read;
    Read  = Write;
    Write = Swap;

    // Brace yourselves for another iteration!
    num_blocks /= 2;
    blocksize  *= 2;
    pos         = 0;
    
  }

  // Our final sorted index is stored in 'Read'
  for (i=0; i<FM->num_seqs; i++)
    Index[i] = Read[i];
  
  free(Read);
  free(Write);
  
  
}




//
// ------------------------------------------------------------------------------------
//




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  GetIndexFname
//
//  About: Set the name of our output index, based on the input file's name
//
void GetIndexFname (char * fname, FILE_META * FM) {

  uint32_t fname_len = (uint32_t)strlen(fname);

  FM->fname = malloc((fname_len+1) * sizeof(char));
  strcpy(FM->fname,fname);
  FM->fname[fname_len] = 0;

  FM->index_fname = malloc((fname_len+5) * sizeof(char));
  strcpy(FM->index_fname,fname);
  strcat(FM->index_fname,".hsi"); // Hexxxy Sequence Index
  FM->index_fname[fname_len+4] = 0;
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  GetMD5
//
//  About: Get the MD5 hash for the input file
//  NOTE : I'm pulling this, since it takes a loooooong time for large files
//
void GetMD5 (char * fname, char * md5) {

  int fname_len = strlen(fname);

  char cmd[fname_len+4];
  strcpy(cmd,"md5 ");
  strcat(cmd,fname);
  
  int i;
  for (i=0; i<fname_len; i++)
    cmd[i+4] = fname[i];
  
  FILE * md5output = popen(cmd,"r");
  fscanf(md5output,"MD5 %s = %s",cmd,md5);
  pclose(md5output);

  md5[32] = 0;
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  IndexAlreadyExists
//
//  About: Check whether an index file already exists by comparing the filesizes
//
int IndexAlreadyExists (FILE_META * FM) {

  FILE * inf = fopen(FM->index_fname,"r");
  if (inf) {

    // Sorry to disturb you ma'am -- just a routine file opening
    fclose(inf);

    // Instead of checking the MD5 hash matches (SLOW!) we'll compare
    // the filesizes.
    uint64_t filesize;
    GetFileSize(FM->fname,&filesize);
    if (filesize == FM->filesize)
      return 1;
    return 0;
    
    // Well, *a* file exists with that name -- does it match our md5?
    /*
    char md5[33];
    fscanf(inf,"MD5: %s",md5);
    md5[32] = 0;
    if (!strcmp(FM->md5,md5))
      return 1;
    */
    
  }

  // Today, Alex learned that trying to close a NULL file pointer
  // causes a segfault -- it seems so obvious in retrospect!
  //fclose(inf);

  // Nope! Not a chance
  return 0;

}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  GetFileSize
//
//  About: Check the size of a file
//
void GetFileSize (char * fname, uint64_t * filesize) {

  int  cmd_len = strlen(fname) + 26;
  char cmd[cmd_len];
  strcpy(cmd,"ls -l ");
  strcat(cmd,fname);
  strcat(cmd," | awk '{print $5}'");
  cmd[cmd_len-1] = 0;

  FILE * sizecheck = popen(cmd,"r");
  fscanf(sizecheck,"%llu",filesize);
  pclose(sizecheck);

}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  ResizeFMSeqs
//
//  About: Realloc the file metadata (if we read in more sequences than initially expected)
//
void ResizeFMSeqs (FILE_META * FM, uint32_t * max_seqs) {
  
  SEQ_META ** TM = malloc((*max_seqs) * sizeof(SEQ_META *));

  uint32_t i,j;
  for (i=0; i<*max_seqs; i++) {
    SEQ_META * SM = FM->Seqs[i];
    TM[i] = malloc(sizeof(SEQ_META));
    TM[i]->SeqName = malloc((strlen(SM->SeqName)+1) * sizeof(char));
    strcpy(TM[i]->SeqName,SM->SeqName);
    TM[i]->seq_length = SM->seq_length;
    TM[i]->alphabet = SM->alphabet;
    TM[i]->seq_start_byte = SM->seq_start_byte;
    TM[i]->num_benchmarks = SM->num_benchmarks;
    if (TM[i]->num_benchmarks) {
      TM[i]->Benchmarks = malloc(SM->num_benchmarks * sizeof(uint32_t));
      for (j=0; j<TM[i]->num_benchmarks; j++)
	TM[i]->Benchmarks[j] = SM->Benchmarks[j];
    }
  }

  *max_seqs *= 2;

  FM->Seqs = realloc(FM->Seqs,*max_seqs * sizeof(SEQ_META *));

  for (i=0; i<*max_seqs/2; i++)
    FM->Seqs[i] = TM[i];

  free(TM);
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION: ResizeBenchmarks
//
//  About: Resize the space we've allocated for byte benchmarks for a given sequence
//
void ResizeBenchmarks (SEQ_META * S, uint32_t * max_benchmarks) {

  uint32_t Tmp[*max_benchmarks];
  
  uint32_t i;
  for (i=0; i<*max_benchmarks; i++)
    Tmp[i] = S->Benchmarks[i];

  *max_benchmarks *= 2;
  S->Benchmarks = realloc(S->Benchmarks,*max_benchmarks*sizeof(uint32_t));

  for (i=0; i<*max_benchmarks/2; i++)
    S->Benchmarks[i] = Tmp[i];
    
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION: DetermineAlphabet
//
//  About: Determine the alphabet of a sequence according to the residues observed
//
char DetermineAlphabet (char * Residues) {

  int i;
  int num_residues = 0;
  for (i=0; i<26; i++) {
    if (Residues[i])
      num_residues++;
  }

  // Are we only using nucleotide characters? (incl. 'N' as ambiguity)
  if (Residues[0] + Residues[2] + Residues[6] + Residues[13] + Residues[19] + Residues[20] == num_residues) {
    // The ultimate question: U or T?
    if (Residues[20])
      return 'R';
    return 'D';
  }

  // Looks protein-y to me!
  return 'P'; 
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  GetSeqMetadata
//
//  About: Extract the metadata for each of the sequences in our file
//
void GetSeqMetadata (FILE_META * FM) {

  int i;

  // How many sequences should we initially allocate space for?
  uint32_t max_seqs = 20;
  FM->Seqs = malloc(max_seqs*sizeof(SEQ_META *));
  FM->num_seqs = 0;

  // Open the file
  FILE * inf = fopen(FM->fname,"r");

  // A pointer we'll use to track the current sequence's metadata
  SEQ_META * SM;

  // How many bytes deep are we into this file?
  uint64_t byte_depth = 0;

  // How many benchmarks will we start off allocating space for?
  uint32_t std_max_benchmarks = 20;

  // What's the length of the longest sequence name?
  uint32_t max_seq_name_len = 0;

  // What residues have we seen?
  char Residues[26];

  // Set the Alphabets flags to 0
  FM->Alphabets[0] = 0;
  FM->Alphabets[1] = 0;
  FM->Alphabets[2] = 0;
  
  // A bajillion variables we'll use while running through the file
  uint32_t max_benchmarks;
  uint32_t dist_to_next_bm;
  int line_len;
  size_t burner = 0;
  char * line = NULL;
  uint32_t * seq_len;
  uint32_t dist_from_start;

  // Go through line by line -- fun!
  while ((line_len = getline(&line,&burner,inf)) != -1) {

    // Is this the start of a new sequence?
    if (line[0] == '>') {

      // If we've already been parsing a sequence, figure out what alphabet
      // it belonged to and reset the Residues array
      if (FM->num_seqs) {
	char alphabet;
	if (FM->set_alpha) alphabet = FM->set_alpha;
	else               alphabet = DetermineAlphabet(Residues);
	SM->alphabet = alphabet;
	if      (alphabet == 'D') { FM->Alphabets[0] = 1; }
	else if (alphabet == 'N') { FM->Alphabets[1] = 1; }
	else               /* P */{ FM->Alphabets[2] = 1; }
      }
      for (i=0; i<26; i++)
	Residues[i] = 0;

      // Do we need to resize our file metadata (upon seeing more sequences than
      // initially expected)?
      if (FM->num_seqs == max_seqs)
	ResizeFMSeqs(FM,&max_seqs);

      // How long is this sequence's name (just name -- skip comments!)
      int name_len = 1;
      while (line[name_len+1] > 32)
	name_len++;

      // Is this the longest sequence name we've seen?
      if (name_len > max_seq_name_len)
	max_seq_name_len = (uint32_t)name_len;

      // Allocate the sequence metadata struct we'll use for this seq
      FM->Seqs[FM->num_seqs] = malloc(sizeof(SEQ_META));
      SM = FM->Seqs[FM->num_seqs];

      // Copy that dang name over!
      SM->SeqName = malloc((name_len+1)*sizeof(char));
      SM->SeqName[name_len] = 0;
      while (name_len) {
	SM->SeqName[name_len-1] = line[name_len];
	name_len--;
      }

      // Record some information to kick off our indexing
      SM->seq_start_byte = byte_depth+(uint64_t)line_len;
      SM->seq_length = 0;
      seq_len = &(SM->seq_length);

      // Allocate our initial array of benchmarks for this sequence
      max_benchmarks = std_max_benchmarks;
      SM->Benchmarks = malloc(max_benchmarks*sizeof(uint32_t));
      SM->num_benchmarks = 0;

      // How long do we have to go (in characters) until the next benchmark?
      dist_to_next_bm = FM->benchmark_dist;

      // How many bytes have we traveled from the beginning?
      dist_from_start = 0;

      // This is another sequence, alright!
      FM->num_seqs += 1;

    } else if (FM->num_seqs) {

      // We're in a sequence!
      i=0;
      while (line[i] > 32) {

	// Increase the sequence length, note that we've moved forward
	*seq_len += 1;
	dist_from_start++;

	// Note the residue ("corrected" from ASCII)
	line[i] -= 65;
	if (line[i] > 25)
	  line[i] -= 32;
	Residues[line[i]] = 1;

	// We're one step closer to the next benchmark!
	dist_to_next_bm--;

	// Uh-oh -- is it benchmark time?!
	if (!dist_to_next_bm) {

	  // Do we need to resize?
	  if (SM->num_benchmarks == max_benchmarks)
	    ResizeBenchmarks(SM,&max_benchmarks);

	  // BENCHMARK TIME!
	  SM->Benchmarks[SM->num_benchmarks] = dist_from_start;
	  SM->num_benchmarks++;

	  // Got a ways to go again :p
	  dist_to_next_bm = FM->benchmark_dist;

	}
	
	i++;

      }

      // We stepped forward to hit the (presumably newline) less-than-32 character
      //dist_from_start++;

      // In case there's some garbage spacing in the file, be sure to catch it!
      dist_from_start += line_len - i;
      
    }

    // We're really cruising!
    byte_depth += line_len;
    
  }

  // Closing up shop
  fclose(inf);

  // Record the longest sequence name's length
  FM->max_seq_name_len = max_seq_name_len;

  // Quick sanity check
  if (FM->num_seqs == 0) {
    fprintf(stderr,"  ERROR: No (FASTA-formatted) sequence data identified in '%s'\n",
	    FM->fname);
  } else {
    // Catch the final alphabet
    char alphabet;
    if (FM->set_alpha) alphabet = FM->set_alpha;
    else               alphabet = DetermineAlphabet(Residues);
    SM->alphabet = alphabet;
    if      (alphabet == 'D') { FM->Alphabets[0] = 1; }
    else if (alphabet == 'N') { FM->Alphabets[1] = 1; }
    else               /* P */{ FM->Alphabets[2] = 1; }
  }

}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION: GetFileMetadata
//
//  About: Get the metadata from our file
//
int GetFileMetadata (char * infname, FILE_META * FM, int force) {

  GetIndexFname(infname,FM);

  //GetMD5(FM->fname,FM->md5);

  // Check whether we've already indexed this file
  GetFileSize(FM->fname,&(FM->filesize));
  if (!force && IndexAlreadyExists(FM))
    return 1;

  GetSeqMetadata(FM);
  return 0;
  
}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  WriteIndex
//
//  About: Write the index information out to a file
//
void WriteIndex (FILE_META * FM) {

  uint32_t i,j;

  // DEBUGGING
  //FILE * outf = stdout;

  FILE * outf = fopen(FM->index_fname,"w");

  // What alphabets are present in this file?
  int alphabet_code = FM->Alphabets[0] + 2 * FM->Alphabets[1] + 4 * FM->Alphabets[2];

  fprintf(outf,"FILE SIZE: %016llx\n",FM->filesize);             // 12 + 16 =  28 bytes
  //fprintf(outf,"MD5: %s\n",FM->md5); // 6 + 32 =  38 bytes
  fprintf(outf,"NUM SEQS: %08x\n",FM->num_seqs);                 // 11 +  8 =  19 bytes
  fprintf(outf,"ALPHABETS: %d\n",alphabet_code);                 // 12 +  1 =  13 bytes
  fprintf(outf,"LONGEST SEQ NAME: %08x\n",FM->max_seq_name_len); // 19 +  8 =  27 bytes
  fprintf(outf,"BENCHMARK DIST: %08x\n",FM->benchmark_dist);     // 17 +  8 =  25 bytes
  //------------------------------------------------------------------------= 112 bytes

  // This 'preamble' takes 112 bytes, each listing will be
  // (2 + FM->max_seq_name_len + 1 + 8 + 1 + 1 + 1 + 16 + 1) bytes
  uint64_t seq_start_byte = 112 + (uint64_t)(FM->num_seqs * (FM->max_seq_name_len+31));

  // NOTE that we need to pre-compute these for each sequence as it appears in the
  // original (potentially unsorted order), since the listings will be alphabetically
  // sorted in the 'index' segment
  uint64_t IndexStartByte[FM->num_seqs];
  for (i=0; i<FM->num_seqs; i++) {

    IndexStartByte[i] = seq_start_byte;

    // Compute the location of the next sequence's start byte
    // (within *the index*, not the sequence file) using the
    // count that we see in the next section
    SEQ_META * SM = FM->Seqs[i];
    seq_start_byte += (uint64_t)(110 + FM->max_seq_name_len + (9 * SM->num_benchmarks));

  }

  // Before we print out the list of sequence names, we'll want to
  // organize them into a lexicographic sorting.
  uint32_t * Index = malloc(sizeof(uint32_t) * FM->num_seqs);
  LexSort(FM,Index);

  // NOW we can list off each of the sequences in sorted order,
  // along with its length and the byte index that it starts at
  // in this file!
  for (i=0; i<FM->num_seqs; i++) {
    
    SEQ_META * SM = FM->Seqs[Index[i]];
    
    fprintf(outf,": %s",SM->SeqName);
    for (j=strlen(SM->SeqName); j<FM->max_seq_name_len; j++)
      fprintf(outf," ");
    fprintf(outf," %08x %c %016llx\n",SM->seq_length,SM->alphabet,IndexStartByte[Index[i]]);

  }

  // The shortest-lived stars burn brightest
  free(Index);

  // FINALLY, the easy stuff!  Tell 'em the sequence position indices in
  // the original FASTA!
  for (i=0; i<FM->num_seqs; i++) {

    SEQ_META * SM = FM->Seqs[i];
    
    fprintf(outf,"SEQ NAME: %s",SM->SeqName);
    for (j=strlen(SM->SeqName); j<FM->max_seq_name_len; j++)
      fprintf(outf," ");
    fprintf(outf,"\n");                                            // 11 + FM->max_seq_name_len bytes
    fprintf(outf,"SEQ LENGTH: %08x\n",SM->seq_length);             // 13 +  8 = 21 bytes
    fprintf(outf,"SEQ ALPHABET: %c\n",SM->alphabet);               // 15 +  1 = 16 bytes
    fprintf(outf,"SEQ START BYTE: %016llx\n",SM->seq_start_byte);  // 17 + 16 = 33 bytes
    fprintf(outf,"SEQ NUM BENCHMARKS: %08x\n",SM->num_benchmarks); // 21 +  8 = 29 bytes
    for (j=0; j<SM->num_benchmarks; j++)
      fprintf(outf,"%08x\n",SM->Benchmarks[j]);                    // SM->num_benchmarks * (8 + 1) bytes
    //------------------------------------------------------------- 110 + FM->max_seq_name_len + (9 * SM->num_benchmarks) bytes
  }

  fclose(outf);

}




///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  BuildSeqIndex
//
//  About: Build an index on the sequences from our input file
//
void BuildSeqIndex (char * infname, uint32_t benchmark_dist, char set_alpha, int force) {
  FILE_META * FM = malloc(sizeof(FILE_META));
  FM->set_alpha = set_alpha;
  FM->benchmark_dist = benchmark_dist;
  if (GetFileMetadata(infname,FM,force))
    return;
  WriteIndex(FM);
}



//
// ------------------------------------------------------------------------------------
//



///////////////////////////////////////////////////////////////////////////////////
//
//  FUNCTION:  ParseHSIMetadata
//
//  About: Parse an HSI file's metadata component
//
void ParseHSIMetadata (FILE_META * FM) {

  char smallhex[9];
  char bighex[17];

  FILE * inf = fopen(FM->index_fname,"r");
  if (inf == NULL) {
    fprintf(stderr,"  Error:  Failed to open index file '%s'\n",FM->index_fname);
    exit(1);
  }
  
  // Grab the file size
  fscanf(inf,"FILE SIZE: %s\n",bighex);
  FM->filesize = strtoll(bighex,NULL,16);

  // If this filesize doesn't match the sequence file, create a new index
  uint64_t fsize_check;
  GetFileSize(FM->fname,&fsize_check);
  if (fsize_check != FM->filesize) {
    fclose(inf);
    BuildSeqIndex(FM->fname,HSI_STD_BENCHMARK_DIST,0,1);
    inf = fopen(FM->index_fname,"r");
    fseek(inf,28,0);
  }
  
  //fscanf(inf,"MD5: %s\n",FM->md5);

  fscanf(inf,"NUM SEQS: %s\n",smallhex);
  FM->num_seqs = strtol(smallhex,NULL,16);

  int alphabet_code;
  fscanf(inf,"ALPHABETS: %d\n",&alphabet_code);

  fscanf(inf,"LONGEST SEQ NAME: %s\n",smallhex);
  FM->max_seq_name_len = strtol(smallhex,NULL,16);

  fscanf(inf,"BENCHMARK DIST: %s\n",smallhex);
  FM->benchmark_dist = strtol(smallhex,NULL,16);

  // Parse the alphabet code
  int j=4;
  while (j) {
    if (alphabet_code >= j) {
      FM->Alphabets[j/2] = 1;
      alphabet_code -= j;
    } else {
      FM->Alphabets[j/2] = 0;
    }
    j /= 2;
  }

  fclose(inf);
  
}






// EOF



