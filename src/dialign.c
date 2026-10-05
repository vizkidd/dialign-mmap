
                 /************************\
                 *                        *
                 *     DIALIGN 2.2.1      *
                 *                        *
                 *       dialign.c        * 
                 *                        *
                 *       written by       *
                 *                        *
                 *    B. Morgenstern      *
                 *                        *
                 \************************/




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#include <sys/stat.h>
#include "dialign.h"
#include "define.h"
#include "alig_graph_closure.h"
#include "mmstore.h"
#include "pratique.h"



FILE *fp_dia, *fp_dpa, *fp_frg , *fp_mot ; 
struct multi_frag *anchor_frg ;

int col_score = 0; 
int char_num[ MAX_REGEX ] ;
char *mot_char[ MAX_REGEX ] ; 
int regex_len , mot_len = 0 ; 


clock_t beg_pa , end_pa , beg_ali , end_ali , beg_ts , end_ts ;
float time_diff_pa , time_diff_ali , perc_pa_time , time_diff_srt ; 
float total_pa_time = 0 ; 


float mot_factor , mot_offset_factor , max_mot_offset ; 

int wgt_type_plot = 0 , motifs = 0 ; 
int bubblesort = 0 , cd_gobics = 0 ; 
int nas = 0 , ref_seq = 0 , i_max ; 
int speed_optimized = 0 ; 
int online = 0 ; 
int time_stamps = 0 ; 
int break1 = 0 ; 
int break2 = 0 ; 
int wgt_print = 0 ; 
int wgt_print_x = 0 ; 
short max_itnum = MAX_ITNUM ; 
int quali_num = 1 ; 
int wgt_plot = 0 ; 
int self_comparison = 0;
short exclude_frg = 0; 
int ***exclude_list ;
int max_sim_score = -2000 ; 
int sf_mat = 0 ; 
char nuc1, nuc2, nuc3 ;
short crick_strand = 0;
int frg_count = 0; 
int dna_speed = 0;
char pst_name[NAME_LEN];
int cont_it = 1 , wgt_type = 0  ;
int mask = 0, strict = 0 , textual_alignment = 1;
char prn[ NAME_LEN ] ;
int redundant, print_max_nd = 1;
int lmax = MAX_DIA;
char **arguments;
int  pr_av_nd = 0, pr_av_max_nd ;
char input_line[ NAME_LEN ];
char input_parameters[ NAME_LEN ];
int print_status = 0 ;
char clust_sim[NAME_LEN] ;
float tot_weight = 0, av_len;
int anchors = 0;
int pa_only = 0;
int dia_num = 0;
int max_dia_num = 0;
float av_dia_num = 0;
float av_max_dia_num = 0;
int afc_file = 0;
int afc_filex = 0;
int dia_pa_file = 0;
int frag_file = 0;
int argnum;
int standard_out = 0;
int plot_num = 4 ;
int default_name = 1;
int fasta_file = 0;
int cw_file = 0; 
int msf_file = 0;
char *upg_str;
int dcount = 0;
char tmp_dir[NAME_LEN] = "";  /* -tmpdir: directory of the mmap store */
int no_tree = 0;            /* -notree: do not compute the guide tree */
int no_mmap = 0;             /* -nommap: use ordinary heap memory */
int num_threads = 0;         /* -threads: 0 = OpenMP default */


int **shift; 
int   thr_sim_score = 4 ;
char **seq;              /* sequences (rows live in the mmap store) */
int sim_score[21][21];  /* similarity matrix */
float av_sim_score_pep ;
float av_sim_score_nuc ;
float **glob_sim;        /* overall similarity between any two sequences */
float **wgt_prot  ;      /* `weight' of diagonals */
float **wgt_dna   ;      /* `weight' of diagonals */
float **wgt_trans ;      /* `weight' of diagonals */
float **min_weight;      /* `weight' of diagonals */
int min_dia = MIN_DIA ;             /* minimum length of diagonals */
int max_dia = MAX_DIA ;  /* maximum length of diagonals */
int iter_cond_prob = 0;
int *seqlen;                /* lengths of sequences */
char **full_name ;
unsigned char **stop_it_p;   /* stop_it_p[i][j] = 1: no further alignment of the pair (i,j) is
                            needed in this iteration step (inverse of the original
                            short cont_it_p[i][j]; 0 = fresh zero page, no init) */ 
float score;
int maxlen;              /* maximum length of sequences */
int seqnum;              /* number of sequences */
int *num_dia_bf;         /* num_dia_bf[ istep ] = number of diagonals from
                            all pairwise alignments BEFORE FILTER
                            PROCEDURE in iteration step `istep' */     
int *num_dia_af;         /* num_dia_af[istep] = number of diagonals from
                            all pairwise alignments AFTER FILTER 
                            PROCEDURE in iteration step `it' */     
int num_dia_anc;         /* number of diagonals definde by anchored 
                            regions */
int num_all_it_dia = 0;  /* total number of diagonals in multiple alignment 
                            in all iteration steps */
float weight_sum_bf;     /* sum of weights  of diagonals in multiple 
                            alignment before filter procedure */  
float weight_sum_af;     /* sum of weights  of diagonals in multiple 
                            alignment after fliter procedure*/
float threshold = 0.0 ;  /* threshold T */
int num_dia_p;           /* number of diagonals in pairwise alignment */ 
int long_output = 0;     /* if long_output = 1, a log-file is produced.  */   
int frg_mult_file = 0 ; 
int frg_mult_file_v = 0 ; 
int overlap_weights = 1 ;  
int ow_force = 0 ;
int anc_num = 0;          /* number of anchored regions 
                            (specified in file *.anc) */
int par_count;           /* number of parameters       */
float pairalignsum;      /* sum of weights in pairwise alignment */ 
int pairalignlen;        /* sum of aligned residues in pairwise alignment */
char amino_acid[22];
int istep;  
struct multi_frag        /* pointer to first diagonal in multiple alignment */
      *this_it_dia;      /* in current iteration step */  
struct multi_frag        /* last diagonal of the previous iteration step (NULL in step 1) */
      *this_it_prev = NULL,
      *list_tail = NULL;
struct multi_frag        /* pointer to first diagonal in multiple alignment */
      *all_it_dia;       /* in all iteration step */
struct multi_frag *end_dia;  
                         /* pointer to last diagonal in multiple alignment */

char par_dir[NAME_LEN];
char **seq_name;
char mat_name[NAME_LEN];         /* name of file containing similarity matrix */
char mat_name_p[NAME_LEN];
char anc_name[NAME_LEN];  /* anchored regions */
char seq_file[NAME_LEN];
char input_name[NAME_LEN];
char tmp_str[NAME_LEN];
char output_name[NAME_LEN];
char printname[NAME_LEN];
char mot_regex[MAX_REGEX] ; 

char *par_file;

short **mot_pos ;       /* positions of pre-defined motifs */ 

unsigned char **amino;  /* amino acid residues in protein sequences or 
                          translated DNA sequences, respective */

unsigned char **amino_c;  /* amino acid residues on crick strand */ 
 
CLOSURE *clos;         /* closure data structure for GABIOS-LIB */

unsigned char **closed_pos_row;  /* closed_pos_row[i] -> N rows of (seqlen[i]+2) flags,
                          all stored in ONE block of the mmap store.
                          CLOSED_POS(i,j)[p] == 0 <=> the p-th residue of 
                          sequence i is not yet directly (by one diagonal) 
                          aligned with any residue of sequence j (this is
                          the original DIALIGN table `open_pos[i][j][p]`,
                          stored inverted: at the beginning of the first
                          iteration step every flag is 0, i.e. a fresh,
                          never-touched zero page of the backing file, so
                          no time or disk space is spent initialising the
                          O(N^2 * L) table).  In the subsequent iteration
                          steps, only those parts of the sequence are 
                          considered that are not yet aligned. */     

  
struct multi_frag *pair_dia;   /* diagonals in pairwise alignemnt */


double **tp400_prot ;    /* propability distribution for sums of similarity
                       socores in diagonals occurring in comparison matrix
                       (by random experiments and approximation  */

double **tp400_dna ;    /* propability distribution for sums of similarity
                       socores in diagonals occurring in comparison matrix
                       (by random experiments and approximation  */

double **tp400_trans ;    /* propability distribution for sums of similarity
                       socores in diagonals occurring in comparison matrix
                       (by random experiments and approximation  */


char dia_pa_name[NAME_LEN];
char frag_file_name[NAME_LEN];
char mot_file_name[NAME_LEN];


/********************************/
/* prototypes                   */
/********************************/

 extern float mot_dist_factor ( int offset , float parameter ) ;
 extern int word_count( char *seq ) ; 
 extern void subst_mat(char *file_name, int fragno , struct  multi_frag *smp );
 extern int seq_read( char *in_file , char ***sq , char ***sqn , char ***fsqn) ;
 extern int anc_read( char *file_name ) ;
 extern int multi_anc_read( char *file_name ) ;
 extern void randomize( int r_numb , FILE *fp1 );
 extern int mini2(int a, int b);
 extern int maxi2(int a, int b);
 extern int mini3(int a, int b, int c);
 extern int num_test( char *cp );
 extern void mini(int *a, int b);
 extern void maxi(int *a, int b);
 extern void filter( int *num, struct multi_frag *vector );
 extern void throw_out( float *weight_sum );
 extern void sel_test();
 extern void regex_parse( char *regex ) ;
 extern void seq_parse( char *regex ) ;
 extern void exclude_frg_read( char *file_name , int ***exclude_list ) ;
 extern void ow_all( struct multi_frag *first , int number , int nthreads );
 extern float frag_chain( int n1 , int n2 , FILE *fp , FILE *fp2, int *num , struct multi_frag **pd , int par_safe );
 extern void rel_wgt_calc( int l1 , int l2 , float **wgt , int type );
 extern void para_read( int num , char **arg ); 
 extern void frag_sort( int number , struct multi_frag *dp , int olw );
 extern void ow_frag_sort( int number , struct multi_frag *dp , int olw );
 extern void bubble_sort( int number , struct multi_frag *dp );
 extern void ow_bubble_sort( int number , struct multi_frag *dp );
 extern void seq_shift();
 extern int translate(char c1, char c2, char c3, int s , int i);
 extern char invert( char c1 ) ;
 extern int int_test(float f);
 extern int match_test( struct multi_frag *dia, int mn);
 
 extern void para_print(char *s_f, FILE *f);
 extern void ali_arrange(int fragno , struct  multi_frag *smp, FILE *fp, FILE *fp2, FILE *fp3 , FILE *fp4 , FILE *fp_csc );
 extern void print_log( struct multi_frag *d , FILE *fp_l , FILE *fp_fs);
 extern void print_fragments( struct multi_frag *d , FILE *fp_frg );
 extern void tp400_read( int wgt_type , double **pr_ptr );
 extern void ow_add(struct  multi_frag *sm1 , struct  multi_frag *sm2);
 extern void av_tree_print();
 extern void matrix_read( FILE *fp_mat ) ;
 extern void mem_alloc( ) ;
 

                    /******************************/
                    /*           main             */
                    /******************************/



/* Output-safety check (KNOWN_ISSUES.md #1): when the input file's own
   extension is .fa/.seq/.fasta, that extension is stripped to build the
   output-file stem and then re-suffixed (see the ".fa"/".ali"/".ms"/".cw"
   handling below); for -fa output this reconstructs the exact original
   input filename whenever the input was itself named "*.fa" -- which is
   the ordinary, expected naming convention for FASTA files, so this is not
   a rare edge case. fopen(...,"w") on that name would silently truncate
   the user's input sequences and replace them with the alignment output.
   Uses (st_dev,st_ino), not string comparison, so it is robust to
   relative vs. absolute paths and symlinks; if either file can't be
   stat()'d (output doesn't exist yet, e.g.) there is nothing to collide
   with. */
static int would_clobber_input( const char *out_name , const char *in_name )
{
  struct stat sb_in, sb_out;
  if( stat( in_name , &sb_in )  != 0 ) return 0;
  if( stat( out_name , &sb_out ) != 0 ) return 0;
  return ( sb_in.st_dev == sb_out.st_dev && sb_in.st_ino == sb_out.st_ino );
}

/* .ali/.ms/.cw never collide with a stripped .fa/.seq/.fasta input in
   practice (their suffixes aren't in the strip list), so this is a cheap
   defensive backstop, not a change to any normal code path: refuse rather
   than silently truncate in whatever unanticipated setup would trigger it. */
static void refuse_if_would_clobber_input( const char *out_name , const char *in_name )
{
  if( ! would_clobber_input( out_name , in_name ) ) return;
  char msg[ 2 * NAME_LEN + 200 ];
  sprintf( msg,
    "\n\n  refusing to write output file `%s': it is the same file as\n"
    "  the input `%s'. Writing it would overwrite your input sequences.\n"
    "  Rename the input file (or its extension) and re-run.\n\n",
    out_name, in_name );
  erreur( msg );
}

/* For -fa specifically, the collision is the common case (any "*.fa" input),
   so rather than refuse -fa outright for the ordinary naming convention,
   fall back to a name that can't collide and say so on stderr. out_name is
   updated in place; alt_buf must be at least NAME_LEN+32 bytes. */
static void avoid_input_collision( char *out_name , const char *in_name , char *alt_buf )
{
  if( ! would_clobber_input( out_name , in_name ) ) return;
  if( strlen( in_name ) + strlen( ".dialign-aligned.fa" ) >= NAME_LEN ) {
    /* pathological input-name length: fall back to a hard refusal rather
       than risk overflowing the NAME_LEN-sized itname2 buffer */
    refuse_if_would_clobber_input( out_name , in_name );
    return;
  }
  sprintf( alt_buf , "%s.dialign-aligned.fa" , in_name );
  fprintf( stderr,
    "\n  note: -fa output `%s' would overwrite the input file `%s';\n"
    "  writing it to `%s' instead.\n\n",
    out_name , in_name , alt_buf );
  strcpy( out_name , alt_buf );
}

int main(int argc, char **argv)
{
 int k,  anc1, dia_counter, tmpi1, tmpi2 ;

 struct multi_frag *current_dia, *diagonal1, *diagonal2, *anc_dia;  
                        /* pointers to diagonals in multiple alignment */ 

 char str[NAME_LEN], dist_name[NAME_LEN]; 
 char par_str[NAME_LEN];  
 char *char_ptr;
 char prn2[NAME_LEN];
 char logname[NAME_LEN];
 char fsm_name[NAME_LEN];
 char dia_name[NAME_LEN];
 char csc_name[NAME_LEN];
 char itname[NAME_LEN], itname2[NAME_LEN], itname3[NAME_LEN];
 char itname4[NAME_LEN];
 char dialign_dir[NAME_LEN];

 int i, j, hv, sv, fv; 


 FILE *fp_ali, *fp2, *fp3, *fp4, *fp_log, *fp_fsm, *fp_st , *fp_csc ; 
 FILE *fp_matrix ;               /* file containing similarity matrix */

 strcpy(mat_name,MATNAME);
 strcpy( clust_sim , "av" );
 
 par_file = (char *) calloc((size_t) NAME_LEN , sizeof(char) );


 if( time_stamps ) 
   beg_ali = clock() ; 

 strcpy ( dialign_dir , "DIALIGN2_DIR" );

 par_file = getenv(dialign_dir);
#ifdef DIALIGN2_DIR_DEFAULT
 /* Only defined by the WASM build (wasm/build.sh): a browser has no shell
    environment to export DIALIGN2_DIR from, and reaching Emscripten's
    emulated environment from page JavaScript varies between Emscripten
    versions, so the build bakes in where it preloads the data files.
    An explicitly set DIALIGN2_DIR still wins. Not defined for native
    builds, where this whole block compiles away. */
 if( par_file == NULL )
   par_file = DIALIGN2_DIR_DEFAULT;
#endif
 if( par_file == NULL )
   {
     printf("\n \n \n    Please set the environmentvariable DIALIGN2_DIR \n");
     printf("    as described in the README file \n"); 
     exit(1);
   }

 argnum = argc;

 strcpy( par_dir , par_file );

if(argc == 1)
  {
    printf("\n    usage: %s [ options ] <seq_file> \n\n", argv[0] );
    printf("    <seq_file> contains input sequences in FASTA format.\n"); 
    printf("    Per default, sequences are assumed to be protein sequences.\n" ) ;
    printf("    For DNA alignment, please use one of these options: \n\n");
    printf("     -n    DNA sequences; similarity calculated at the nucleotide level \n\n"); 
    printf("     -nt   DNA sequences; similarity calculated at the peptide level\n");
    printf("           (by translation using the genetic code) \n\n");
    printf("     -lgs  long genomic sequences: Both nucleotide and peptide\n");
    printf("           similarities calculated \n\n");  
    printf("    Many more options are available, please consult the \n");
    printf("    DIALIGN USER_GUIDE that should come with the DIALIGN package.\n");
    printf("    For more information on DIALIGN, please visit the DIALIGN\n"); 
    printf("    home page at BiBiServ (Bielefeld Bioinformatic Server): \n\n") ;
    printf("        http://bibiserv.techfak.uni-bielefeld.de/dialign/ \n\n");    
    exit(1) ;
  }

 arguments = ( char ** ) calloc( argnum , sizeof ( char * ) );

 for( i = 0 ; i < argnum ; i++ )
   {
     arguments[i] = ( char *)  calloc( NAME_LEN , sizeof (char) );
     strcpy( arguments[i] , argv[i] );
   }
 


 strcpy( input_name , argv[ argc - 1 ] );
  
 threshold = 0.0 ;


 para_read( argnum , arguments );
#ifdef _OPENMP
 if( num_threads > 0 ) omp_set_num_threads( num_threads );
#endif

 if( ( textual_alignment == 0 ) && ( col_score == 1 ) ) { 
   printf("\n\n   Option -csc makes sense only if \"textual alignment\"");
   printf(" is produced. \n");
   printf("   This can be enforced with option -ta \n\n");
   printf("   program terminated \n\n\n");
   exit(1) ;
 } 


 if( cd_gobics ) {
 strcpy( input_line , "program parameters:  " ) ; 
 for( i = 1 ; i < ( argnum -1 ) ; i++ ) {
     strcat( input_line , argv[i] );
     strcat( input_line , " " );
   }
 }
 else {
 strcpy( input_line , "program call:  " ) ; 
 for( i = 0 ; i < argnum ; i++ ) {
     strcat( input_line , argv[i] );
     strcat( input_line , " " );
   }
 }


 if ( wgt_type > 0 )  
   strict = 1 ; 

 strcpy( seq_file , input_name );

 if(
        ( ! strcmp( input_name + strlen( input_name ) - 4 , ".seq" ) )
     || ( ! strcmp( input_name + strlen( input_name ) - 3 , ".fa" ) )
     || ( ! strcmp( input_name + strlen( input_name ) - 6 , ".fasta" ) )
   )
 if( ( char_ptr = strrchr(input_name,'.') ) != NULL)
   *char_ptr = '\0';


 strcpy( anc_name , input_name );
 strcat( anc_name , ".anc" );

 /* ---- backing store for all large data structures ---- */
 if( no_mmap )
   mm_set_enabled( 0 );
 else {
   char mmdir[ NAME_LEN ] ;
   if( tmp_dir[0] )
     strcpy( mmdir , tmp_dir ) ;
   else if( getenv( "DIALIGN_TMPDIR" ) != NULL )
     mmdir[0] = '\0' ;                      /* mmstore reads the variable */
   else {                                   /* default: next to the input */
     char *slash ;
     strcpy( mmdir , seq_file ) ;
     if( ( slash = strrchr( mmdir , '/' ) ) != NULL ) {
       if( slash == mmdir ) slash[1] = '\0' ; else *slash = '\0' ;
     }
     else
       strcpy( mmdir , "." ) ;
   }
   if( mmdir[0] ) mm_set_dir( mmdir );
 }

 seqnum = seq_read( seq_file , &seq , &seq_name , &full_name ) ;

 if ( motifs )
   regex_parse( mot_regex ) ; 


 if( ( seqnum == 2 ) && ( iter_cond_prob == 0 ) ) 
   max_itnum = 1 ; 

 
     if(  ( ow_force == 0 ) && ( seqnum > OVERLAP_THRESHOLD )  )
       overlap_weights = 0;
     if( seqnum == 2 )
       overlap_weights = 0;

  if( seqnum < 2 ) { 

    if( cd_gobics ) {
      printf("\n\n         Something is wrong with your sequence file. Maybe you entered a\n");
      printf("         MS WORD or RFT file or your file contains only one single sequence.\n");
      printf("         Please note that our server only accepts plain text files. \n\n");  
      printf("         For more information, please consult our online manual \n");
      printf("         at the CHAOS/DIALIGN home page:\n\n");  
      printf("             http://dialign.gobics.de/chaos-dialign-manual");
    }

    else { 
      printf("\n\n         Your sequence file containes only a single sequence.\n");
      printf("         Please make sure your input file contains at least two sequences.\n\n");
      printf("         For more information, please consult the online manual \n");
      printf("         at the DIALIGN home page: \n\n");
      printf("             http://bibiserv.techfak.uni-bielefeld.de/dialign/manual.html ");
    }



    printf("\n       \n \n \n \n");
    exit(1);
  }

  maxlen = 0;

  

  stop_it_p = (unsigned char **) mm_matrix( seqnum , seqnum , sizeof(unsigned char) );





  for( i = 0 ; i < seqnum ; i++ )
   {
    av_len = av_len + seqlen[i];

    if( seqlen[i] == 0 )
      {
        printf("\n \n \n                       WARNING: \n \n");
        printf("          Sequence %d contains no residues.\n",i+1);
        printf("          Please inspect the sequence file.\n \n ");
        printf("\n \n          Program terminated \n \n \n " );     

        exit(1);
      }
 
    if(maxlen < seqlen[i])
       maxlen = seqlen[i];
   }

  av_len = av_len / seqnum;

  if ( motifs )
    seq_parse( mot_regex ) ; 
  
  seq_shift();


   glob_sim = (float **) mm_matrix( seqnum , seqnum , sizeof(float) );

   strcpy(par_str,"sdfsdf");

   if( argc > 1 )
   {
   strcpy(str,par_dir);
   strcat(str,"/");
   strcat(str,mat_name);
   strcpy(mat_name_p,str);
   
   if( (fp_matrix = fopen(mat_name_p, "r")) == NULL)
   {


   printf("\n\n Cannot find the file %s \n\n", mat_name );
   printf(" Make sure the environment variable DIALIGN2_DIR points\n");
   printf(" to a directory containing the files \n\n");
   printf("   BLOSUM \n   tp400_dna\n   tp400_prot \n   tp400_trans \n\n" );
   printf(" These files should be contained in the DIALIGN package \n\n\n" ) ;
   exit(1) ;




     printf("\n \n \n \n              ATTENTION ! \n \n");
     printf("\n   There is no similarity matrix `%s'. \n", mat_name);
     printf("   in the directory \n \n");
     printf("           %s\n \n", par_dir);
     exit(1);
   }
   }


    if( sf_mat && wgt_type == 1 ) {
      printf("\n\n  Option -mat needs peptide-level similarities; it cannot be used with -n \n\n");
      exit(1);
    }

    if( wgt_type != 1 )
      matrix_read( fp_matrix );

    mem_alloc(  );


    if( wgt_type != 1 ) {
      amino = (unsigned char **) mm_calloc( seqnum , sizeof(unsigned char *) );
      for( i = 0 ; i < seqnum ; i++ ) 
        amino[i] = (unsigned char *) mm_calloc( (size_t) seqlen[i]+5 , sizeof(unsigned char) );
    }

    if( crick_strand ) { 
      amino_c = (unsigned char **) mm_calloc( seqnum , sizeof(unsigned char *) );
      for( i = 0 ; i < seqnum ; i++ )
        amino_c[i] = (unsigned char *) mm_calloc( (size_t) seqlen[i]+5 , sizeof(unsigned char) );
    }
 

             /******************************************************  
             *                                                     *      
             *  read file, that contains data of anchored regions  *
             *                                                     *      
             ******************************************************/  



if( anchors ) {
  multi_anc_read( input_name );
}

if( exclude_frg ) { 

  exclude_list = (int ***) mm_calloc( seqnum , sizeof(int **) );
  for(i = 0 ; i < seqnum ; i++ ) 
    exclude_list[ i ] = (int **) mm_calloc( seqnum , sizeof(int *) );
  for(i = 0 ; i < seqnum ; i++ ) 
  for(j = 0 ; j < seqnum ; j++ ) 
    exclude_list[ i ][ j ]  = (int *) mm_calloc( (size_t) seqlen[ i ] + 1 , sizeof(int) );

  exclude_frg_read ( input_name , exclude_list ) ;
}



   if( wgt_type == 0 ) 
     tp400_read( 0 , tp400_prot);
   if( wgt_type % 2 )
     tp400_read( 1 , tp400_dna );
   if( wgt_type > 1 )
     tp400_read( 2 , tp400_trans );



           /****************************\
           *                            * 
           *    Name of output files    *  
           *                            * 
           \****************************/

   if( default_name )
     {
       strcpy( printname , input_name);
       strcpy( prn , printname);
     } 
   else
     {  
       strcpy( printname , output_name );
       strcpy( prn , printname);
     }
    

   strcpy(prn2 , prn); 
  
   if( default_name )
     strcat(prn,".ali");

   strcat(prn2,".fa");  
    


   strcpy(logname,printname);
   strcat(logname,".log");

   strcpy(fsm_name , printname);
   strcat(fsm_name,".fsm");

   if( print_status ) {
     strcpy( pst_name , printname );
     strcat( pst_name,".sta");
   }    

   if( afc_file )
     {
       strcpy( dia_name , printname );  
       strcat( dia_name , ".afc" );
       fp_dia = fopen( dia_name , "w" );
       fprintf(fp_dia,"\n #  %s \n\n  seq_len: " , input_line );
       for( i = 0 ; i < seqnum ; i++ )
         fprintf(fp_dia,"  %d ", seqlen[i] );
       fprintf(fp_dia,"\n\n");

     }

   if( col_score ) { 
     strcpy( csc_name , printname );  
     strcat( csc_name , ".csc" );
     fp_csc = fopen( csc_name , "w" );
   }

   if( dia_pa_file )
     {
       strcpy( dia_pa_name , printname );  
       strcat( dia_pa_name , ".fop" );

       fp_dpa = fopen( dia_pa_name , "w" );


       fprintf(fp_dpa,"\n #  %s \n\n  seq_len: " , input_line );
       for( i = 0 ; i < seqnum ; i++ ) 
         fprintf(fp_dpa,"  %d ", seqlen[i] ); 
       fprintf(fp_dpa,"\n\n");
       fclose( fp_dpa ) ;
     }


   if( motifs ) {
     strcpy( mot_file_name , printname );  
     strcat( mot_file_name , ".mot" );
     fp_mot = fopen( mot_file_name , "w" );
      
     fprintf(fp_mot,"\n #  %s \n\n   " , input_line );
     fprintf(fp_mot," motif: %s \n\n", mot_regex ); 
     fprintf(fp_mot," max offset for motifs = %d \n\n", (int) max_mot_offset ); 
     fprintf(fp_mot," the following fragments contain the motif: \n\n" ); 
     fprintf(fp_mot,"   seq1 seq2    beg1 beg1 len    wgt" ); 
     fprintf(fp_mot,"   # mot    mot_wgt  \n\n" ); 
   }


   if( frag_file ) {
     strcpy( frag_file_name , printname );  
     strcat( frag_file_name , ".frg" );
     fp_frg = fopen( frag_file_name , "w" );
      
     fprintf(fp_frg,"\n #  %s \n\n  seq_len: " , input_line );
     for( i = 0 ; i < seqnum ; i++ )
       fprintf(fp_frg,"  %d ", seqlen[i] );
     fprintf(fp_frg,"\n  sequences: " );
     for( i = 0 ; i < seqnum ; i++ )
       fprintf(fp_frg,"  %s ", seq_name[i] );

     fprintf(fp_frg ,"\n\n");
   }



  clos = newAligGraphClosure(seqnum, seqlen, 0, NULL);

  {
    /* one block for all N*N rows; sequence i owns N rows of (seqlen[i]+2) bytes */
    size_t total = 16 ;
    closed_pos_row = (unsigned char **) mm_calloc( seqnum , sizeof(unsigned char *) );
    for( i = 0 ; i < seqnum ; i++ )
      total += (size_t) seqnum * ( (size_t) seqlen[i] + 2 ) ;
    {
      unsigned char *blk = (unsigned char *) mm_calloc( total , 1 );
      size_t off = 0 ;
      for( i = 0 ; i < seqnum ; i++ ) {
        closed_pos_row[i] = blk + off ;
        off += (size_t) seqnum * ( (size_t) seqlen[i] + 2 ) ;
      }
    }
  }

   	  /**************************************
          *                                     *
          *      definition of  `amino'         *       
    	  *                                     *
          **************************************/




  if( wgt_type > 1 ) 
    for(hv=0;hv<seqnum;hv++)
    for(i=1;i<=seqlen[hv]-2;i++)
      {


        if( translate( seq[hv][i],seq[hv][i+1],seq[hv][i+2],hv,i ) == -1)
          exit(1);


        amino[hv][i] = translate( seq[hv][i],seq[hv][i+1],seq[hv][i+2],hv,i);
   
        if( crick_strand ) { 
          nuc1 = invert( seq[hv][i+2] );
          nuc2 = invert( seq[hv][i+1] );
          nuc3 = invert( seq[hv][i] );
 
          amino_c[hv][i] = translate( nuc1 , nuc2 , nuc3 , hv , i);
        }
      }


   if( wgt_type == 0 ) 
   for(hv=0;hv<seqnum;hv++)
   for(i=1;i<=seqlen[hv];i++)
    {
     if( seq[hv][i] == 'C' ) amino[hv][i] = 1;           
     if( seq[hv][i] == 'S' ) amino[hv][i] = 2;           
     if( seq[hv][i] == 'T' ) amino[hv][i] = 3;           
     if( seq[hv][i] == 'P' ) amino[hv][i] = 4;           
     if( seq[hv][i] == 'A' ) amino[hv][i] = 5;           
     if( seq[hv][i] == 'G' ) amino[hv][i] = 6;           
     if( seq[hv][i] == 'N' ) amino[hv][i] = 7;           
     if( seq[hv][i] == 'D' ) amino[hv][i] = 8;           
     if( seq[hv][i] == 'E' ) amino[hv][i] = 9;           
     if( seq[hv][i] == 'Q' ) amino[hv][i] = 10;           
     if( seq[hv][i] == 'H' ) amino[hv][i] = 11;           
     if( seq[hv][i] == 'R' ) amino[hv][i] = 12;           
     if( seq[hv][i] == 'K' ) amino[hv][i] = 13;           
     if( seq[hv][i] == 'M' ) amino[hv][i] = 14;           
     if( seq[hv][i] == 'I' ) amino[hv][i] = 15;           
     if( seq[hv][i] == 'L' ) amino[hv][i] = 16;           
     if( seq[hv][i] == 'V' ) amino[hv][i] = 17;           
     if( seq[hv][i] == 'F' ) amino[hv][i] = 18;           
     if( seq[hv][i] == 'Y' ) amino[hv][i] = 19;           
     if( seq[hv][i] == 'W' ) amino[hv][i] = 20;           
    }


     
     amino_acid[0] = 'X';           
     amino_acid[1] = 'C';           
     amino_acid[2] = 'S';           
     amino_acid[3] = 'T';           
     amino_acid[4] = 'P';           
     amino_acid[5] = 'A';           
     amino_acid[6] = 'G';           
     amino_acid[7] = 'N';           
     amino_acid[8] = 'D';           
     amino_acid[9] = 'E';           
     amino_acid[10] = 'Q';           
     amino_acid[11] = 'H';           
     amino_acid[12] = 'R';           
     amino_acid[13] = 'K';           
     amino_acid[14] = 'M';           
     amino_acid[15] = 'I';           
     amino_acid[16] = 'L';           
     amino_acid[17] = 'V';           
     amino_acid[18] = 'F';           
     amino_acid[19] = 'Y';           
     amino_acid[20] = 'W';



num_dia_anc = anc_num * (seqnum-1);




if ( anchors ) {

  if( time_stamps )
    beg_ts = clock() ;

  if ( nas == 0 ) 
    if( bubblesort ) 
      bubble_sort ( anc_num , anchor_frg ) ;  
    else
      frag_sort ( anc_num , anchor_frg , 0 ) ;


  if( time_stamps) {
    end_ts = clock() ;
    time_diff_srt = (float) ( end_ts - beg_ts ) / CLOCKS_PER_SEC ;
    if( time_stamps )
      printf (" for anc: time_diff_srt = %f \n", time_diff_srt );
  }


  filter( &anc_num , anchor_frg);
/*  exit(1) ; 
*/
}


if(long_output)
  {
    fp_log = fopen(logname,"w");
    fprintf(fp_log,"\n #  %s \n\n   " , input_line );
  }

if(frg_mult_file) {
  fp_fsm = fopen(fsm_name,"w");
  fprintf(fp_fsm,"\n #  %s \n\n" , input_line );
}





if( 
    (  num_dia_bf = (int *) calloc( ( max_itnum + 1 )  ,  sizeof( int ) )  )
      == NULL
  )
     {
          printf(" problems with memory allocation for `num_dia_bf' !  \n \n");
          exit(1);
     }


if( 
    (  num_dia_af = (int *) calloc( ( max_itnum + 1 )  ,  sizeof( int ) )  )
      == NULL
  )
     {
          printf(" problems with memory allocation for `num_dia_af' !  \n \n");
          exit(1);
     }


all_it_dia = frag_new(); 
current_dia = all_it_dia;


  strcpy(itname,printname);
  strcpy(itname2,printname);
  strcpy(itname3,printname);
  strcpy(itname4,printname);
  sprintf(str,".ali");
  
  if( default_name )
  strcat(itname,str);
  
  sprintf(str,".fa");
  strcat(itname2,str);


  if( msf_file )
  strcat(itname3,".ms");
  
  
  if( cw_file )
  strcat(itname4,".cw");
  





        if( textual_alignment ) {
          refuse_if_would_clobber_input( itname , seq_file );
          fp_ali = fopen(itname,"w");
        }
 
        if( standard_out )
          fp_ali = stdout;

   
        if( textual_alignment )
        if(fasta_file) {
          char itname2_alt[ NAME_LEN + 32 ];
          avoid_input_collision( itname2 , seq_file , itname2_alt );
          fp2 = fopen(itname2,"w");
        }

        if(msf_file) {
          refuse_if_would_clobber_input( itname3 , seq_file );
          fp3 = fopen(itname3,"w");
        }

        if(cw_file) {
          refuse_if_would_clobber_input( itname4 , seq_file );
          fp4 = fopen(itname4,"w");
        }

        if( textual_alignment )
          para_print(seq_file , fp_ali);


                    /***************************\
                    *                           * 
                    *      ITERATION START      *   
                    *                           * 
                    \***************************/









istep = 0 ; 
 while( ( cont_it == 1 ) && ( istep < max_itnum ) ) 
  {

    cont_it = 0 ;
    istep++ ; 

/* printf("\n  istep = %d \n", istep ); */


    this_it_dia = current_dia;
    this_it_prev = list_tail;           /* last node of the previous steps (NULL in step 1) */
 
    strcpy(itname,printname);
    strcpy(itname2,printname);
    strcpy(itname3,printname);
    strcpy(itname4,printname);
    sprintf(str,".ali");

        
    if( default_name )
      strcat(itname,str); 

    sprintf(str,".fa");
    strcat(itname2,str); 

    if( msf_file )
      strcat(itname3,".ms"); 
   
    
    if( cw_file )
      strcat(itname4,".cw"); 
  
    weight_sum_af = 0;
    num_dia_bf[ istep ] = 0;
 
    if( time_stamps ) 
      beg_pa = clock(); 


    if( ref_seq == 0 ) 
      i_max = seqnum ; 
    else
      i_max = 1 ; 

    /* ---------------------------------------------------------------
       pairwise alignments.  The pairs are processed in blocks; inside a
       block they run in parallel (OpenMP) when nothing forces a serial
       run, and the results are appended to the fragment list strictly in
       the original (i,j) order, so the output is independent of the
       number of threads.
       --------------------------------------------------------------- */
    {
      const int PBLOCK = 512 ;
      int par_ok = ! ( iter_cond_prob || afc_file || dia_pa_file || motifs ||
                       print_status || long_output || wgt_print || wgt_print_x ||
                       pr_av_max_nd ) ;
      int  *bi = (int *) malloc( PBLOCK * sizeof(int) ) ;
      int  *bj = (int *) malloc( PBLOCK * sizeof(int) ) ;
      int  *bn = (int *) malloc( PBLOCK * sizeof(int) ) ;
      struct multi_frag **bd = (struct multi_frag **) malloc( PBLOCK * sizeof(struct multi_frag *) ) ;
      int nb = 0, b, last_i = -1, last_j = -1 ;
      int done = 0 ;

      i = 0 ; j = 1 ;
      while( ! done )
        {
          /* collect the next block of pairs that have to be aligned */
          nb = 0 ;
          while( nb < PBLOCK )
            {
              if( i >= i_max ) { done = 1 ; break ; }
              if( j >= seqnum ) { i++ ; j = i + 1 ; continue ; }
              if( ! stop_it_p[ i ][ j ] ) { bi[nb] = i ; bj[nb] = j ; nb++ ; }
              j++ ;
            }
          if( nb == 0 ) break ;

#pragma omp parallel for schedule(dynamic,1) if( par_ok && nb > 1 )
          for( b = 0 ; b < nb ; b++ )
            {
              int np = 0 ;
              struct multi_frag *pd = NULL ;
              frag_chain( bi[b] , bj[b] , fp_ali , fp_mot , &np , &pd , par_ok ) ;
              bn[b] = np ;
              bd[b] = pd ;
            }

          for( b = 0 ; b < nb ; b++ )
            {
              int num_p = bn[b] ;
              struct multi_frag *pd = bd[b] ;

              for(k=0;k<num_p;k++)
                {
                  *current_dia = pd[k];

                  SETNX( current_dia , frag_new() );
                  end_dia = current_dia;   
                  current_dia = NX(current_dia);

                }

              num_dia_bf[ istep ] = num_dia_bf[ istep ] + num_p;

              for(hv=0; hv<num_p;hv++)
                weight_sum_af = weight_sum_af + (pd[hv]).weight;

              if(num_p)
                free(pd);
            }
          last_i = bi[nb-1] ; last_j = bj[nb-1] ;
        }

      /* The global weight tables must hold what the last aligned pair left
         in them (ow_add() reads them later); with private per-thread tables
         that has to be restored explicitly. */
      if( par_ok && last_i >= 0 ) {
        if( wgt_type == 0 )  rel_wgt_calc( seqlen[last_i] , seqlen[last_j] , wgt_prot , 0 );
        if( wgt_type % 2 )   rel_wgt_calc( seqlen[last_i] , seqlen[last_j] , wgt_dna  , 1 );
        if( wgt_type > 1 )   rel_wgt_calc( seqlen[last_i] , seqlen[last_j] , wgt_trans, 2 );
      }
      free( bi ) ; free( bj ) ; free( bn ) ; free( bd ) ;
    }


    if( time_stamps ) { 
      end_pa = clock(); 

      time_diff_pa = (float) ( end_pa - beg_pa ) / CLOCKS_PER_SEC ; 
      if( time_stamps )
        printf (" time_diff_pa = %f \n", time_diff_pa ); 
      total_pa_time = total_pa_time + time_diff_pa;
    }




    if( break1 ) {
      printf("\n  break1\n");
      exit(1) ;
    }


/*
    if( pa_only ) {
      printf("\n\n istep = %d, pa finished - exit \n\n", istep );       
      exit(1);
    }
*/

    if(overlap_weights)
      {
        diagonal1 = this_it_dia;
        dia_counter = 0;

        if( print_status ) {
          /* status file is updated every 100 diagonals: serial loop */
          if( diagonal1 != NULL )   
          while( NX(diagonal1) != NULL )   
            {
              dia_counter++;
              if( ( dia_counter % 100 ) == 0 )
                {                
                  fp_st = fopen( pst_name ,"w");

                  fprintf(fp_st," dsd  %s \n", input_line);
                  fprintf(fp_st,"\n\n\n    Status of the program run:\n");  
                  fprintf(fp_st,"    ==========================\n\n");  
                  if( seqnum > 2 ) { 
                    fprintf(fp_st,"      iteration step %d in ", istep); 
                    fprintf(fp_st,"multiple alignment\n" );
                  }
                  fprintf(fp_st,"      calculating overlap weight for diagonals\n");
                  fprintf(fp_st,"      current diagonal = %d\n\n", dia_counter );
                  fprintf(fp_st,"      total number of"); 
                  fprintf(fp_st," diagonals: %d\n\n\n\n", num_dia_bf[ istep ]);
                  fclose(fp_st);
                }

              diagonal2 = NX(diagonal1);

              while(NX(diagonal2) != NULL) 
                {
                  if( diagonal1->trans == diagonal2->trans ) 
                    ow_add(diagonal1 , diagonal2); 
                  diagonal2 = NX(diagonal2);        
                }
              diagonal1 = NX(diagonal1);        
            }
        }
        else {
          int nth = 1 ;
#ifdef _OPENMP
          nth = omp_get_max_threads() ;
#endif
          ow_all( diagonal1 , num_dia_bf[ istep ] , nth ) ;
        }
        if( bubblesort )   
          ow_bubble_sort( num_dia_bf[ istep ] , this_it_dia ); 
        else 
          frag_sort( num_dia_bf[ istep ] , this_it_dia , overlap_weights ); 
      }
    else /* no overlap_weights */ {
      beg_ts = clock() ; 

      if( bubblesort ) 
        bubble_sort( num_dia_bf[ istep ] , this_it_dia );
      else 
        frag_sort( num_dia_bf[ istep ] , this_it_dia , overlap_weights );

      end_ts = clock() ; 
      time_diff_srt = (float) ( end_ts - beg_ts ) / CLOCKS_PER_SEC ;
      if( time_stamps )
        printf (" time_diff_srt = %f \n", time_diff_srt );
  }


    num_dia_af[ istep ] = num_dia_bf[ istep ];
    weight_sum_bf = weight_sum_af;

    pairalignsum = 0;
    pairalignlen = 0;


    filter( num_dia_af + istep , this_it_dia); 
    num_all_it_dia = num_all_it_dia + num_dia_af[ istep ];


/*
    if( pa_only == 0 ) {
      printf("\n\n istep = %d, filter finished - exit \n\n", istep );       
      exit(1);
    }
*/


    weight_sum_af = 0;
       
    /* print_log() only produces output for -lo / -fsm / -fsmv (its other
       results, pairalignsum/len, are never used).  It walks the whole
       diagonal list once per sequence pair, i.e. O(N^2 * D), so it is
       skipped when it cannot write anything. */
    if( long_output || frg_mult_file )
      print_log( this_it_dia , fp_log , fp_fsm );

    if( frag_file )
      print_fragments( this_it_dia , fp_frg );

    throw_out( &weight_sum_af );

    sel_test( );

     

    threshold = threshold ;

    if( break2 ) {
      printf("\n  break2\n");
      exit(1) ;
    }


  } /* while ( cond_it == 1 ) */  


                    /***************************\
                    *                           * 
                    *       ITERATION END       *   
                    *                           * 
                    \***************************/

strcpy( dist_name , printname);
strcat(dist_name , ".dst");





if ( ref_seq == 0 ) {
  if( no_tree ) {
    upg_str = (char *) mm_calloc( 2 , sizeof(char) ) ;
  }
  else
    av_tree_print();
}



      if( standard_out )
        fp_ali = stdout;
      
if(sf_mat){
  subst_mat( input_name , num_all_it_dia , all_it_dia ) ;
}


if( textual_alignment )
  ali_arrange( num_all_it_dia , all_it_dia , fp_ali , fp2, fp3, fp4, fp_csc );


if(long_output) 
  {
/*     fprintf(fp_log "\n\n thr = %f , lmax = %d , speed = %f  */  
    fprintf(fp_log, "\n\n    total sum of weights: %f \n\n\n", tot_weight);
    fclose(fp_log);
  }




if( argnum == 1 )
  {
    printf("\n     Program terminated normally\n");
    printf("     Results are contained in file `%s' \n \n \n", itname);
  }


av_dia_num = 2 * dia_num ;
av_dia_num = av_dia_num / ( seqnum * ( seqnum - 1) ) ;

av_max_dia_num = 2 * max_dia_num ;
av_max_dia_num = av_max_dia_num / ( seqnum * ( seqnum - 1) ) ;



tmpi1 = av_dia_num ;
tmpi2 = av_max_dia_num ;

if(pr_av_nd)
  printf(" %d ", tmpi1 );

if(pr_av_max_nd)
  printf(" %d ", tmpi2 );



/* pr_av_nd/pr_av_max_nd (-pand/-o) print into fp_ali, but fp_ali is only
   ever opened when textual_alignment is set (see the fopen(itname,"w")
   above) -- some modes (-lgs/-lgs_t/-nta) turn textual_alignment off
   without disabling -pand/-o. Without this guard fp_ali is an
   uninitialized stack FILE*, and fprintf()ing through it is undefined
   behaviour: reproduced as a deterministic SEGV under a controlled build,
   and as a non-deterministic crash/hang across repeated runs of an
   unmodified build (the garbage FILE* varies run to run) -- e.g.
   `dialign2-2 -lgs_t -pand seqs.fa`. Confirmed present in the original
   2.2.1 too; this guard changes nothing when textual_alignment is set
   (the ordinary case), matching how the fclose(fp_ali) just below already
   guards on the same condition. */
if( textual_alignment ) {
if(pr_av_nd)
  fprintf(fp_ali, "    %d fragments considered for alignment \n", tmpi1 );

if(pr_av_max_nd)
  fprintf(fp_ali, "    %d fragments simultaneously stored \n\n", tmpi2 );
}

if( textual_alignment )
  fclose(fp_ali);


  if( time_stamps ){ 
    end_ali = clock() ; 
    time_diff_ali = (float) ( end_ali - beg_ali ) / CLOCKS_PER_SEC ;

    perc_pa_time = total_pa_time / time_diff_ali * 100 ; 
    printf (" time_diff_ali = %f \n", time_diff_ali ); 
    printf (" total_pa_time = %f \n", total_pa_time );
    printf (" corresponds to %f percent \n\n", perc_pa_time );
  } 
}    /* main */



