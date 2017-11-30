#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include <zlib.h>
#define MAX_ID_LEN 255
#define MAX_SEQ_LEN 1023
#define HASH_SIZE 8589934592

static char *data[] = { "alpha", "bravo", "charlie", "delta",
			"echo", "foxtrot", "golf", "hotel", "india", "juliet",
			"kilo", "lima", "mike", "november", "oscar", "papa",
			"quebec", "romeo", "sierra", "tango", "uniform",
			"victor", "whisky", "x-ray", "yankee", "zulu" };

/* Type to hold the forward and reverse read
   of a sequence pair with quality scores */
typedef struct sqp {
  char id[MAX_ID_LEN+1];
  char seq[MAX_SEQ_LEN+1];
  char qual[MAX_SEQ_LEN+1];
} Sqp;
typedef struct sqp* SQP;

int is_gz( const char* fq_fn );
SQPDB init_SQPDB( size_t size );
int read_fastq( FILE* fastq, SQP curr_seq );
int gz_read_fastq( gzFile fastq, SQP curr_seq );
FILE * fileOpen(const char *name, char access_mode[]);

int main(int argc, char* argv[]) {

  extern char* optarg;
  extern int optin;
  char fq_fn[MAX_ID_LEN+1];
  int ich;
  int IS_GZ;
  FILE* fq_p;
  gzFile fq_gz_p;
  SQP curr_seq;

  if ( argc == 1 ) {
    help();
  }

  while( (ich=getopt( argc, argv, "f:k:q:" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_fn, optarg );
      break;
    case 'k' :
      k = atoi( optarg );
      break;
    case 'q' :
      qual_cut = atoi( optarg );
      break;
    case 'h' :
      help();
    }
  }

  /* Initialize the SQP curr_seq */
  curr_seq = (SQP)malloc(sizeof( Sqp ));

  /* Initialize the hash */
  int hcreate( HASH_SIZE );
  
  IS_GZ = is_gz( fq_fn );
  /* Open the fastq file */
  if ( IS_GZ ) {
    fq_gz_p = gzopen( fq_fn, "r" );
    if ( fq_gz_p == NULL ) {
      help();
    }
  }
  else {
    fq_p = fileOpen( fq_fn, "r" );
    if ( fq_p == NULL ) {
      help();
    }
  }
  while( read_fastq_somehow( IS_GZ, fq_p, fq_gz_p, curr_seq ) ) {
    add_fastq_kmers( curr_seq );
  }

  close_fastq_somehow( IS_GZ, fq_p, fq_gz_p );


  /* Write out something */
  ENTRY e, *ep;
  int i;
  
  hcreate(30);
  
  for (i = 0; i < 24; i++) {
    e.key = data[i];
    /* data is just an integer, instead of a
       pointer to something */
    e.data = (void *) i;
    ep = hsearch(e, ENTER);
    /* there should be no failures */
    if (ep == NULL) {
      fprintf(stderr, "entry failed\n");
      exit(EXIT_FAILURE);
    }
  }
  
  for (i = 22; i < 26; i++) {
    /* print two entries from the table, and
       show that two are not in the table */
    e.key = data[i];
    ep = hsearch(e, FIND);
    printf("%9.9s -> %9.9s:%d\n", e.key,
	   ep ? ep->key : "NULL", ep ? (int)(ep->data) : 0);
  }
  hdestroy();
  exit(EXIT_SUCCESS);
}

/* is_gz - takes pointer to a string understood to
   be a filename. Returns true (1) iff the string
   ends in ".gz"
   Returns false (0) otherwise
*/
int is_gz( const char* fq_fn ) {
  size_t fn_len;
  fn_len = strlen( fq_fn );
  if ( (fq_fn[fn_len-3] == '.') &&
       (fq_fn[fn_len-2] == 'g') &&
       (fq_fn[fn_len-1] == 'z') ) {
    return 1;
  }
  return 0;
}

/* add_fastq_kmers
   Args: SQP curr_seq
         int k
   Returns: void
   Does: adds all the kmers in the input SQP that pass
         filtering to the current hash
*/
void add_fastq_kmers( SQP curr_seq, int k ) {
  size_t len;

}

/* init_SQPDB
   Initialize and return a SQPDB
   Args: (1) size_t size - how big the array of sqps should be
*/
SQPDB init_SQPDB( size_t size ) {
  size_t i;
  SQPDB sqpdb;
  SQP first_seq;

  /* Try to allocate memory for SQPDB */
  sqpdb = (SQPDB)malloc(sizeof(Sqpdb));
  if ( sqpdb == NULL ) {
    fprintf( stderr, "Not enough memories for database" );
    return NULL;
  }

  /* Try to allocate memory for sqps array */
  first_seq = (SQP)malloc( size * sizeof(Sqp) );
  sqpdb->sqps = (SQP*)malloc( size * sizeof(SQP) );
  
  if ( (first_seq   == NULL) ||
       (sqpdb->sqps == NULL) ) {
    fprintf( stderr, "Not enough memories for database" );
    return NULL;
  } 
    
  for( i = 0; i < size; i++ ) {
    sqpdb->sqps[i] = &first_seq[i];
  }

  sqpdb->seq_len   = MAX_SEQ_LEN;
  sqpdb->num_reads = 0;
  sqpdb->size      = size;

  return sqpdb;
}



/* read_fastq
   Return 1 => more sequence to be had
          0 => EOF
 */
int read_fastq( FILE* fastq, SQP curr_seq ) {
  char c;
  size_t i;
  c = fgetc( fastq );
  if ( c == EOF ) return 0;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return 0;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return 0;
    }
    curr_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      curr_seq->id[i] = '\0';
    }
  }
  curr_seq->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
          (c != EOF) ) {
    c = fgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
  while ( (c != '\n') &&
          (c != EOF) &&
          (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      curr_seq->seq[i++] = c;
    }
    c = fgetc( fastq );
  }
  curr_seq->seq[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", curr_seq->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
         (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
  i = 0;
  while( (c != '\n') &&
         (c != EOF) &&
         (i < MAX_SEQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      curr_seq->qual[i++] = c;
    }
    c = fgetc( fastq );
  }
  curr_seq->qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     INIT_ALN_SEQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_SEQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return 0;
  }
  return 1;
}

/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}
