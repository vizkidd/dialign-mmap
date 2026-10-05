
                 /*******************\
                 *                   *
                 *     DIALIGN 2     *
                 *                   *
                 *     input.c       *
                 *                   *
                 \*******************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "define.h"
#include "dialign.h"
#include "alig_graph_closure.h"
#include "mmstore.h"
#include "pratique.h"

extern int max_dia , self_comparison ; 
extern int sim_score[21][21]; 
extern int max_sim_score ;
extern float av_sim_score_pep ;
extern float av_sim_score_nuc ;
extern char par_dir[ NAME_LEN ] ;
extern double **tp400_prot, **tp400_dna, **tp400_trans ; 
extern int *seqlen ;
extern int seqnum ; 


int word_count( char *str ) {

  /* number of blank/tab-separated words in a line (the line terminator
     counts as a separator; an empty string has no words).  The original
     looped up to strlen(str)-1, which underflows for an empty string and
     silently ignores the last character of a line without a final newline. */
  short word = 0 ; 
  size_t i ; 
  int word_len = 0 ; 

  for( i = 0 ; str[i] != '\0' ; i++ ) { 
    if( ( str[i] != ' ' ) && ( str[i] != '\t' ) && 
        ( str[i] != '\n' ) && ( str[i] != '\r' ) ) {  
      if( ! word ) { 
        word_len++ ; 
        word = 1 ; 
      }
    }
    else 
      word = 0 ; 
  }
  return( word_len ) ; 
} 


void exclude_frg_read( char *file_name , int ***exclude_list) {

  char exclude_file_name[ NAME_LEN ] ;
  FILE *fp;
  char line[ 10000 ] ;
  int i, len, beg1, beg2, seq1, seq2; 
  int sv = 0, hv, word_num  ;

  strcpy( exclude_file_name , file_name );
  strcat( exclude_file_name , ".xfr" );

  if( (fp = fopen( exclude_file_name, "r")) == NULL)
    erreur("\n\n cannot find file with excluded fragments \n\n\n");

  

  while( fgets( line , MLINE , fp ) != NULL ) {
    if( strlen( line ) > 4 ) {   
      sscanf(line,"%d %d %d %d %d", &seq1, &seq2, &beg1, &beg2 , &len  );

      if( seq1 > seqnum ){
        printf ("\n\n exclueded fragment makes no sense!\n\n");
        printf (" wrong sequence no %d in fragment\n\n", seq1 );
        printf ("%d %d %d %d %d \n\n ", seq1, seq2, beg1, beg2 , len  );
        exit(1) ;
      }
 
      if( seq2 > seqnum ){
        printf ("\n\n    excluded fragment makes no sense!\n\n");
        printf ("    wrong sequence no %d in fragment\n\n", seq2 );
        printf ("    %d %d %d %d %d \n\n", seq1, seq2, beg1, beg2 , len );
        exit(1) ;
      }

/*
      seq1 = seq1 - 1; 
      seq2 = seq2 - 1;
*/

      if( beg1 + len > seqlen[ seq1 - 1 ] + 1 ){
        printf ("\n\n    excluded fragment makes no sense!\n");
        printf ("    fragment");
        printf ("     \" %d %d %d %d %d \"\n", seq1, seq2, beg1, beg2 , len );
        printf ("    doesn't fit into sequence %d:\n", seq1 );
        printf ("    sequence %d has length =  %d\n\n", seq1 , seqlen[ seq1 - 1 ] );
        exit(1) ;
      }


 
      for( i = 0 ; i < len ; i++ ) {  
        exclude_list[ seq1 - 1 ][ seq2 - 1 ][ beg1 + i ] = beg2 + i ;
      }
    }
  }
  
} /* excluded_frg_read  */ 







void ws_remove( char *str ) {
  int pv = 0 ;

  while( ( str[ pv ] == ' ' ) || ( str[ pv ] == '\t' ) ) {
    pv++ ;
  }

  strcpy( str , str + pv );
}

void n_clean( char *str ) {
  int pv = 0 ;
  char *char_ptr ;

  while( ( str[ pv ] == ' ' ) || 
         ( str[ pv ] == '\t' ) || 
         ( str[ pv ] == '>' ) ) {
    pv++ ;
  }
  strcpy( str , str + pv ) ;

  if( ( char_ptr = strchr( str ,' ') ) != NULL)
    *char_ptr = '\0';
  if( ( char_ptr = strchr( str ,'\t') ) != NULL)
    *char_ptr = '\0';
  if( ( char_ptr = strchr( str ,'\n') ) != NULL)
    *char_ptr = '\0';




}


/*
 *  seq_read()  --  reads the multi-FASTA file.
 *
 *  Semantics are those of the original DIALIGN 2.2.1 reader:
 *    - the first non-blank line must start with '>' (after blanks/tabs),
 *    - only the letters A-Z / a-z of a sequence line are kept (upper-cased),
 *    - the name is the header up to the first blank/tab, and is cut/padded
 *      to SEQ_NAME_LEN characters for the alignment output;
 *  but it is not limited by any line length or sequence count: the file is
 *  mapped read-only and scanned twice (count, then fill), and every
 *  sequence is stored in the mmap store with exactly the room it needs.
 */
static int fasta_scan_open( char *seq_file , const char **buf , size_t *len )
{
  int fd;
  struct stat sb;
  void *p;

  if( ( fd = open( seq_file , O_RDONLY ) ) < 0 ) {
    printf("\n\n Cannot find sequence file %s \n\n\n", seq_file );
    exit(1) ;
  }
  if( fstat( fd , &sb ) != 0 || sb.st_size == 0 ) {
    close( fd );
    erreur("\n\n  file not in FASTA format  \n\n");
  }
  p = mmap( NULL , (size_t) sb.st_size , PROT_READ , MAP_PRIVATE , fd , 0 );
  close( fd );
  if( p == MAP_FAILED )
    erreur("\n\n cannot map sequence file \n\n");
#ifndef __EMSCRIPTEN__
  /* pure performance hint (return value unused, no functional effect);
     not all libc/mmap-emulation targets implement it (e.g. Emscripten's
     WASM build reading from MEMFS), so skip it there rather than risk a
     link failure over an optimization hint. */
  madvise( p , (size_t) sb.st_size , MADV_SEQUENTIAL );
#endif
  *buf = (const char *) p ;
  *len = (size_t) sb.st_size ;
  return 0;
}

/* iterate over lines; returns 0 at end.  [*ls,*le) is the line without '\n',
   with leading blanks/tabs already skipped                                  */
static int fasta_next_line( const char **pos , const char *end ,
                            const char **ls , const char **le )
{
  const char *p = *pos, *e;
  if( p >= end ) return 0;
  e = memchr( p , '\n' , (size_t)( end - p ) );
  *le  = e ? e : end ;
  *pos = e ? e + 1 : end ;
  while( p < *le && ( *p == ' ' || *p == '\t' ) ) p++ ;
  *ls = p ;
  return 1;
}

int seq_read( char *seq_file , char ***sqp , char ***sqnp , char ***fsqnp ) {

  const char *buf, *end, *pos, *ls, *le, *q;
  size_t flen;
  long long *cnt = NULL, *hoff = NULL ;
  size_t cap = 0 ;
  int sn = -1 , i, first = 1 ;
  char **sq, **sqn, **fsqn ;

  fasta_scan_open( seq_file , &buf , &flen );
  end = buf + flen ;

  /* ---- pass 1: check format, count sequences and residues ---- */
  pos = buf ;
  while( fasta_next_line( &pos , end , &ls , &le ) ) {
    if( first ) {
      if( ls == le ) continue ;                    /* blank line          */
      if( *ls != '>' )
        erreur("\n\n  file not in FASTA format  \n\n");
      first = 0 ;
    }
    if( ls < le && *ls == '>' ) {
      sn++ ;
      if( (size_t) sn >= cap ) {
        cap = cap ? cap * 2 : 64 ;
        cnt  = (long long *) realloc( cnt  , cap * sizeof(long long) );
        hoff = (long long *) realloc( hoff , cap * sizeof(long long) );
        if( !cnt || !hoff ) erreur("\n\n problems with memory allocation for the sequence index \n\n");
      }
      cnt[ sn ]  = 0 ;
      hoff[ sn ] = (long long)( ls - buf ) ;
    }
    else
      for( q = ls ; q < le ; q++ )
        if( ( *q >= 65 && *q <= 90 ) || ( *q >= 97 && *q <= 122 ) )
          cnt[ sn ]++ ;
  }
  if( first )
    erreur("\n\n  file not in FASTA format  \n\n");
  sn++ ;                                           /* number of sequences */

  if( self_comparison == 1 && sn != 1 ) {
    printf("\n\n With option \"self comparison\" input file must contain one single sequence \n\n" ); 
    exit(1) ;
  }

  sq   = (char **) mm_calloc( sn + 2 , sizeof(char *) );
  sqn  = (char **) mm_calloc( sn + 2 , sizeof(char *) );
  fsqn = (char **) mm_calloc( sn + 2 , sizeof(char *) );
  seqlen = (int *) mm_calloc( sn + 2 , sizeof(int) );
  if( !sq || !sqn || !fsqn || !seqlen )
    erreur("\n\n problems with memory allocation for `seqlen' \n\n");

  /* ---- pass 2: names ---- */
  for( i = 0 ; i < sn ; i++ ) {
    const char *nb, *ne, *hp = buf + hoff[ i ] ;
    size_t nl, crc ;
    const char *hl = memchr( hp , '\n' , (size_t)( end - hp ) );
    if( !hl ) hl = end ;
    while( hp < hl && ( *hp == ' ' || *hp == '\t' || *hp == '>' ) ) hp++ ;
    nb = hp ;
    while( hp < hl && *hp != ' ' && *hp != '\t' && *hp != '\r' ) hp++ ;
    ne = hp ;
    nl = (size_t)( ne - nb );

    fsqn[ i ] = (char *) mm_calloc( nl + 3 , sizeof(char) );
    memcpy( fsqn[ i ] , nb , nl );

    sqn[ i ] = (char *) mm_calloc( SEQ_NAME_LEN + 3 , sizeof(char) );
    for( crc = 0 ; crc < SEQ_NAME_LEN ; crc++ )
      sqn[ i ][ crc ] = ( crc < nl ) ? nb[ crc ] : ' ' ;
    sqn[ i ][ SEQ_NAME_LEN ] = '\0' ;

    if( cnt[ i ] > 2000000000LL )
      erreur("\n\n sequence too long (more than 2*10^9 residues) \n\n");
    seqlen[ i ] = (int) cnt[ i ] ;
    /* +3: room for the one-position shift done by seq_shift() and a NUL */
    sq[ i ] = (char *) mm_calloc( (size_t) cnt[ i ] + 3 , sizeof(char) );
  }

  /* ---- pass 3: residues ---- */
  {
    int cur = -1 ;
    size_t j = 0 ;
    pos = buf ;
    while( fasta_next_line( &pos , end , &ls , &le ) ) {
      if( ls < le && *ls == '>' ) { cur++ ; j = 0 ; }
      else if( cur >= 0 )
        for( q = ls ; q < le ; q++ )
          if( ( *q >= 65 && *q <= 90 ) || ( *q >= 97 && *q <= 122 ) )
            sq[ cur ][ j++ ] = toupper( (unsigned char) *q ) ;
    }
  }
  munmap( (void *) buf , flen );
  free( cnt );
  free( hoff );

  if( self_comparison ) {
    seqlen[ 1 ] = seqlen[ 0 ] ;
    sq[ 1 ]   = (char *) mm_calloc( (size_t) seqlen[ 0 ] + 3 , sizeof(char) );
    memcpy( sq[ 1 ] , sq[ 0 ] , (size_t) seqlen[ 0 ] + 1 );
    sqn[ 1 ]  = (char *) mm_calloc( SEQ_NAME_LEN + 3 , sizeof(char) );
    strcpy( sqn[ 1 ] , sqn[ 0 ] ) ;
    fsqn[ 1 ] = (char *) mm_calloc( strlen( fsqn[ 0 ] ) + 3 , sizeof(char) );
    strcpy( fsqn[ 1 ] , fsqn[ 0 ] ) ;
    sn++ ;
  }

  *sqp = sq ; *sqnp = sqn ; *fsqnp = fsqn ;
  return( sn );
}



void matrix_read( FILE *fp_mat ) {
  int i, j;
  char line[MLINE], dummy[MLINE];
 
  fgets( line , MLINE , fp_mat );
  fgets( line , MLINE , fp_mat );


  for( i = 1 ; i <= 20 ; i++ ) {
    for(j=i;j<=20;j++) {
      fscanf( fp_mat , "%d" , &sim_score[i][j]);
      sim_score[j][i] = sim_score[i][j];  
      if ( sim_score[i][j] > max_sim_score )
        max_sim_score = sim_score[i][j] ;
    }

    fscanf( fp_mat, "%s\n", dummy);
  }

  fclose(fp_mat);

  for( i = 0 ; i <= 20 ; i++ ) {
    sim_score[i][0] = 0 ;
    sim_score[0][i] = 0 ;
  }

/*
 sim_score[0][0] = max_sim_score ;
*/

}
 


void tp400_read( int w_type , double **pr_ptr ) {  
 
  /* reads probabilities from file */
   /* w_type = 0 (protein), 1 (dna w/o transl.), 2 (dna with transl.) */  

  char line[MLINE], file_name[MLINE], suffix[10], str[MLINE] ;
  int sum, len, max_sim, i ;
  double pr;

  FILE *fp;
 
  if ( w_type == 0 ) {
    strcpy( suffix , "prot" );
  }

  if ( w_type == 1 ) {
    strcpy( suffix , "dna" );
  }  

  if ( w_type == 2 ) {
    strcpy( suffix , "trans" );
  }
 
  strcpy( file_name , par_dir ); 
  strcat( file_name , "/tp400_" );
  strcat( file_name , suffix );


 if ( ( fp = fopen( file_name , "r" ) ) == NULL ) { 
   printf("\n\n Cannot find the file %s \n\n", file_name );    
   printf(" Make sure the environment variable DIALIGN2_DIR points\n");
   printf(" to a directory containing the files \n\n");
   printf("   BLOSUM \n   tp400_dna\n   tp400_prot \n   tp400_trans \n\n" );
   printf(" These files should be contained in the DIALIGN package \n\n\n" ) ;
   exit(1) ;
 }


  if ( fgets( line , MLINE , fp ) == NULL ) 
    { printf("\n\n problem with file %s  \n\n", file_name ); exit(1); }
  else
    if( w_type % 2 )  
      av_sim_score_nuc = atof( line );
    else
      av_sim_score_pep = atof( line );
     

  while( fgets( line , MLINE , fp ) != NULL )
   {
      sscanf(line,"%d %d %s", &len, &sum, str  );

      pr = atof(str);
      pr_ptr[len][sum] = pr;

    }


}    /*  tp400_read  */




