
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
#include "dialign.h"
#include "define.h"
#include "alig_graph_closure.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include <pthread.h>
#include <omp.h>

FILE *fp_dia, *fp_dpa, *fp_frg, *fp_mot;
struct multi_frag *anchor_frg;

int col_score = 0;
int char_num[MAX_REGEX];
char *mot_char[MAX_REGEX];
int regex_len, mot_len = 0;

clock_t beg_pa, end_pa, beg_ali, end_ali, beg_ts, end_ts;
double time_diff_pa, time_diff_ali, perc_pa_time, time_diff_srt;
double total_pa_time = 0;

double mot_factor, mot_offset_factor, max_mot_offset;

int wgt_type_plot = 0, motifs = 0;
int bubblesort = 0, cd_gobics = 0;
int nas = 0, ref_seq = 0, i_max;
int speed_optimized = 0;
int online = 0;
int time_stamps = 0;
int break1 = 0;
int break2 = 0;
int wgt_print = 0;
int wgt_print_x = 0;
short max_itnum = MAX_ITNUM;
int quali_num = 1;
int wgt_plot = 0;
int self_comparison = 0;
short exclude_frg = 0;
int ***exclude_list;
int max_sim_score = -2000;
int sf_mat = 0;
char nuc1, nuc2, nuc3;
short crick_strand = 0;
int frg_count = 0;
int dna_speed = 0;
char pst_name[NAME_LEN];
int cont_it = 1, wgt_type = 0;
int mask = 0, strict = 0, textual_alignment = 1;
char prn[NAME_LEN];
int redundant, print_max_nd = 1;
int lmax = MAX_DIA;
char **arguments;
int pr_av_nd = 0, pr_av_max_nd;
char input_line[NAME_LEN];
char input_parameters[NAME_LEN];
int print_status = 0;
char clust_sim[NAME_LEN];
double tot_weight = 0, av_len;
int anchors = 0;
int pa_only = 0;
int dia_num = 0;
int max_dia_num = 0;
double av_dia_num = 0;
double av_max_dia_num = 0;
int afc_file = 0;
int afc_filex = 0;
int dia_pa_file = 0;
int frag_file = 0;
int argnum;
int standard_out = 0;
int plot_num = 4;
int default_name = 1;
int fasta_file = 0;
int cw_file = 0;
int msf_file = 0;
char *upg_str;
int dcount = 0;

int **shift;
int thr_sim_score = 4;
char *seq[MAX_SEQNUM];    /* sequences */
char *newseq[MAX_SEQNUM]; /* sequences */
int sim_score[21][21];    /* similarity matrix */
double av_sim_score_pep;
double av_sim_score_nuc;
// double **glob_sim;     /* overall similarity between any two sequences */
double **wgt_prot;     /* `weight' of diagonals */
double **wgt_dna;      /* `weight' of diagonals */
double **wgt_trans;    /* `weight' of diagonals */
double **min_weight;   /* `weight' of diagonals */
int min_dia = MIN_DIA; /* minimum length of diagonals */
int max_dia = MAX_DIA; /* maximum length of diagonals */
int iter_cond_prob = 0;
int *seqlen; /* lengths of sequences */
char *full_name[MAX_SEQNUM];
double **pair_score;
short **cont_it_p;
double score;
int maxlen;             /* maximum length of sequences */
int seqnum;             /* number of sequences */
int *num_dia_bf;        /* num_dia_bf[ istep ] = number of diagonals from
                           all pairwise alignments BEFORE FILTER
                           PROCEDURE in iteration step `istep' */
int *num_dia_af;        /* num_dia_af[istep] = number of diagonals from
                           all pairwise alignments AFTER FILTER
                           PROCEDURE in iteration step `it' */
int num_dia_anc;        /* number of diagonals definde by anchored
                           regions */
int num_all_it_dia = 0; /* total number of diagonals in multiple alignment
                           in all iteration steps */
double weight_sum_bf;   /* sum of weights  of diagonals in multiple
                          alignment before filter procedure */
double weight_sum_af;   /* sum of weights  of diagonals in multiple
                          alignment after fliter procedure*/
double threshold = 0.0; /* threshold T */
int num_dia_p;          /* number of diagonals in pairwise alignment */
int long_output = 0;    /* if long_output = 1, a log-file is produced.  */
int frg_mult_file = 0;
int frg_mult_file_v = 0;
int overlap_weights = 1;
int ow_force = 0;
int total_anc_num = 0; /* number of anchored regions
                   (specified in file *.anc) */
int par_count;         /* number of parameters       */
double pairalignsum;   /* sum of weights in pairwise alignment */
int pairalignlen;      /* sum of aligned residues in pairwise alignment */
char amino_acid[22];
int istep;
struct multi_frag /* pointer to first diagonal in multiple alignment */
    *this_it_dia; /* in current iteration step */
struct multi_frag /* pointer to first diagonal in multiple alignment */
    *all_it_dia;  /* in all iteration step */
// struct multi_frag *end_dia;
/* pointer to last diagonal in multiple alignment */

char par_dir[NAME_LEN];
char *seq_name[MAX_SEQNUM];
char mat_name[NAME_LEN]; /* name of file containing similarity matrix */
char mat_name_p[NAME_LEN];
char anc_name[NAME_LEN]; /* anchored regions */
// char anc_name_sorted[NAME_LEN]; /* anchored regions */
char ow_filename[NAME_LEN]; /* overlap weights file */
char seq_file[NAME_LEN];
char input_name[NAME_LEN];
char tmp_str[NAME_LEN];
char output_name[NAME_LEN];
char printname[NAME_LEN];
char mot_regex[MAX_REGEX];

char *par_file;

short **mot_pos; /* positions of pre-defined motifs */

int **amino; /* amino acid residues in protein sequences or
                translated DNA sequences, respective */

int **amino_c; /* amino acid residues on crick strand */

CLOSURE *clos; /* closure data structure for GABIOS-LIB */

int ***open_pos; /* open_pos[i][j][p] = 1, if the p-th residue of
                sequence i is not yet directly (by one diagonal)
                aligned with any residue of sequence j and
                open_pos[i][j][r] = 0 otherwise. So, at the
                beginning of the first iteration step, all values
                are 1. In the subsequent iteration steps,
                only those parts of the sequence are considered,
                that are not yet aligned. */

struct multi_frag *pair_dia; /* diagonals in pairwise alignemnt */

double **tp400_prot; /* propability distribution for sums of similarity
                   socores in diagonals occurring in comparison matrix
                   (by random experiments and approximation  */

double **tp400_dna; /* propability distribution for sums of similarity
                   socores in diagonals occurring in comparison matrix
                   (by random experiments and approximation  */

double **tp400_trans; /* propability distribution for sums of similarity
                   socores in diagonals occurring in comparison matrix
                   (by random experiments and approximation  */

char dia_pa_name[NAME_LEN];
char frag_file_name[NAME_LEN];
char mot_file_name[NAME_LEN];

/********************************/
/* prototypes                   */
/********************************/

extern float mot_dist_factor(int offset, float parameter);
extern int word_count(char *seq);
// extern void subst_mat(char *file_name, int fragno, struct multi_frag *smp);

extern int seq_read(char *in_file, char *sq[MAX_SEQNUM], char **sqn, char **fsqn);
extern int seq_read_mmap(const char *seq_file, char **sq, char **sqn, char **fsqn);

extern mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset);
extern fasta_value *get_seq_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_value *get_seqname_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern openpos_value get_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset);
extern openpos_value set_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, int *value, size_t input_offset);
extern fasta_len_value *get_seqlen_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
// extern ow_value get_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset);
// extern ow_value set_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, double *value, size_t input_offset);
extern diag_value set_diags_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, double *value, size_t input_offset);
extern diag_value get_diagval_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, size_t input_offset);
extern int multi_anc_read_mmap(char *file_name);
extern int get_anchor_count(char *file_name);
extern void para_print_mmap(char *s_f, FILE *fpi, mmapped_file *mapped_fasta);
extern void subst_mat_mmap(char *file_name, int fragno, struct multi_frag *frg, mmapped_file *mapped_fasta);
extern void ali_arrange_mmap(int fragno, struct multi_frag *d, FILE *fp, FILE *fp2, FILE *fp3, FILE *fp4, FILE *fp_col_score, mmapped_file *mapped_fasta);
extern size_t get_seqcount_mmapped_fasta(mmapped_file *mapped_fasta);
// extern void bubble_sort_mmap(int number);
extern void bubble_sort(int number, struct multi_frag *dp, mmapped_file *mapped_anc, mmapped_file *mapped_fasta, mmapped_file *mapped_ow);
extern void bubble_sort_mmap(int number, mmapped_file *mapped_anc, mmapped_file *mapped_fasta, mmapped_file *mapped_ow);
// // extern void frag_sort_mmap(int number, int olw);
// extern void frag_sort_mmap(int number, int olw, mmapped_file *mapped_fasta, mmapped_file *mapped_ow, struct multi_frag *dp);
extern void frag_sort_mmap(int number, int olw, mmapped_file *mapped_fasta, mmapped_file *mapped_anc, mmapped_file *mapped_diags);
// extern void filter_all(int *number);
// extern void filter(int *number, int anchor_step, mmapped_file *mmapped_anc, mmapped_file *mmapped_fasta, mmapped_file *mmapped_diags);
extern void filter(int *number, int anchor_step, mmapped_file *mmapped_anc, mmapped_file *mmapped_fasta, mmapped_file *mmapped_diags, mmapped_file *mmapped_ow);
extern struct multi_frag *get_anchor_empty();
extern struct multi_frag *get_anchor_by_num(mmapped_file *mapped_fasta, mmapped_file *mapped_anc, int *anc_num, size_t offset_anc, int get_next);
extern struct multi_frag *get_anchor_seqs(mmapped_file *mapped_fasta, mmapped_file *mapped_anc, int *seqnum1, int *seqnum2, size_t fasta_offset, int get_next);
extern frag_chain_value *frag_chain_mmap(int n1, int n2, FILE *fp1, FILE *fp_m, int *number, mmapped_file *mapped_fasta, size_t *input_offset, size_t seqlen_n1, char *seq_n1, char *seqname_n1);
extern void ow_add_mmap(int total_diagonals, mmapped_file *mapped_diags, mmapped_file *mapped_fasta);
// extern void frag_sort(mmapped_file *mmapped_fasta, mmapped_file *mmapped_anc, int number, struct multi_frag *dp, int olw);
// extern void ow_add_mmap(struct multi_frag *sm1, struct multi_frag *sm2, mmapped_file *mapped_fasta, mmapped_file *mapped_ow, size_t input_offset);

extern int anc_read(char *file_name);
extern int multi_anc_read(char *file_name);
extern void randomize(int r_numb, FILE *fp1);
extern int mini2(int a, int b);
extern int maxi2(int a, int b);
extern int mini3(int a, int b, int c);
extern int num_test(char *cp);
extern void mini(int *a, int b);
extern void maxi(int *a, int b);
// extern void filter(int *num, struct multi_frag *vector);
// extern void throw_out(double *weight_sum);
extern void throw_out(double *weight_sum, mmapped_file *mapped_diags);
// extern void sel_test();
extern void sel_test(mmapped_file *mapped_diags);
// extern frag_chain_value frag_chain_mmap(int n1, int n2, FILE *fp, FILE *fp2, int *num, char *mapped_fasta, struct stat *sb, size_t input_offset);
extern void para_read(int num, char **arg);
// extern void frag_sort(int number, struct multi_frag *dp, int olw);
extern void ow_frag_sort(int number, struct multi_frag *dp, int olw);
// extern void bubble_sort(int number, struct multi_frag *dp);
extern void ow_bubble_sort(int number, struct multi_frag *dp);
extern void seq_shift();
extern int translate(char c1, char c2, char c3, int s, int i);
extern char invert(char c1);
extern int int_test(double f);
extern int match_test(struct multi_frag *dia, int mn);

// extern void para_print(char *s_f, FILE *f);
// extern void ali_arrange(int fragno, struct multi_frag *smp, FILE *fp, FILE *fp2, FILE *fp3, FILE *fp4, FILE *fp_csc);
extern void print_log(FILE *fp_l, FILE *fp_fs, mmapped_file *mapped_fasta, mmapped_file *mapped_diags);
// extern void print_fragments(struct multi_frag *d, FILE *fp_frg);
extern void print_fragments(mmapped_file *mapped_diags, FILE *fp_ff2);
extern void tp400_read(int wgt_type, double **pr_ptr);
//  extern void ow_add(struct  multi_frag *sm1 , struct  multi_frag *sm2);
// extern void av_tree_print();
// extern void av_tree_print_mmap(mmapped_file *mapped_fasta);
extern void av_tree_print_mmap(mmapped_file *mapped_fasta, mmapped_file *mapped_ow);
extern void matrix_read(FILE *fp_mat);
extern void mem_alloc();

/******************************/
/*           main             */
/******************************/

void create_openpos_file(mmapped_file *mmapped_fasta)
{
  // Writing 1's to file to remove memory limitations. The file can be mmaped later
  //  #UNCOMMENT LATER //DEBUG
  struct stat *sb = &(mmapped_fasta->sb);
  // printf("\n ptr sb: %p \n", &sb); // DEBUG
  // // printf("\n ptr sb int: %d \n", &sb);                                       // DEBUG
  // printf("\n sb.st_size: %d \n", sb->st_size);                              // DEBUG
  // printf("\n mmapped_fasta->sb.st_size: %d \n", mmapped_fasta->sb.st_size); // DEBUG
                                                                            // printf("\n mmapped_fasta->mapped_file[0]: %c \n", mmapped_fasta->mapped_file[0]); // DEBUG
                                                                            //   FILE *open_pos_f;
                                                                            //   // char *file_line = NULL;
                                                                            //   // size_t seq_len_val = -1;
                                                                            //   int seq_file_length = sizeof(seq_file) / sizeof(char);
                                                                            //   char open_pos_fname[seq_file_length]; // = seq_file;
                                                                            //   // printf("\n seq_file: %s \n", seq_file); // DEBUG
                                                                            //   // strcpy(open_pos_fname, seq_file);
                                                                            //   snprintf(open_pos_fname, sizeof(seq_file), "%s", seq_file);
                                                                            //   strcat(open_pos_fname, ".op");
                                                                            //   // printf("\n here 15.5.5 \n"); // DEBUG
                                                                            //   open_pos_f = fopen(open_pos_fname, "w");
                                                                            // #pragma omp parallel for num_threads(omp_get_max_threads()) // schedule(static) // private(file_line, seq_len_val)
                                                                            //   for (int i = 0; i < seqnum; i++)
                                                                            //   {
                                                                            // #pragma omp parallel for shared(i) num_threads(omp_get_max_threads()) // schedule(static) // private(file_line, seq_len_val)
                                                                            //     for (int j = 0; j < seqnum; j++)
                                                                            //     {
                                                                            //       // fprintf(open_pos_f, ">%d,%d:", i, j);
                                                                            //       // for (hv = 1; hv <= seqlen[i]; hv++)
                                                                            //       size_t seq_len_val = get_seqlen_mmapped_fasta(mmapped_fasta->mapped_file, &i, sb);
                                                                            //       // printf("\n seq_len_val: %d \n", seq_len_val);                                     // DEBUG
                                                                            //       char *file_line = (char *)malloc(((seq_len_val + 1) * 2) * sizeof(char));
                                                                            //       strcpy(file_line, "");
                                                                            //       for (int hv = 1; hv < seq_len_val; hv++)
                                                                            //       {
                                                                            //         //  open_pos[i][j][hv] = 1;
                                                                            //         // fprintf(open_pos_f, "%d\t%d\n", hv, 1);
                                                                            //         // fprintf(open_pos_f, "%d,", 1);
                                                                            //         strcat(file_line, "1,");
                                                                            //       }
                                                                            //       // fprintf(open_pos_f, "%d", 1);
                                                                            //       strcat(file_line, "1");
                                                                            // // printf("\n file_line:%s \n", file_line); // DEBUG
                                                                            // #pragma omp critical
                                                                            //       {
                                                                            //         fprintf(open_pos_f, ">%d,%d:%s\n", i, j, file_line);
                                                                            //       }
                                                                            //       // fprintf(open_pos_f, "\n");
                                                                            //       free(file_line);
                                                                            //     }
                                                                            //   }
                                                                            //   // free(file_line);
                                                                            //   fclose(open_pos_f);
                                                                            //   // exit(0); // DEBUG

  // Writing 1's to file to remove memory limitations. The file can be mmaped later
  //  #UNCOMMENT LATER //DEBUG
  FILE *open_pos_f;

  int seq_file_length = sizeof(seq_file) / sizeof(char);
  char open_pos_fname[seq_file_length]; // = seq_file;

  snprintf(open_pos_fname, sizeof(seq_file), "%s", seq_file);
  strcat(open_pos_fname, ".op");

  open_pos_f = fopen(open_pos_fname, "w");

  size_t input_offset = 0;

#pragma omp parallel for shared(input_offset) num_threads(omp_get_max_threads()) // schedule(dynamic)
  for (int i = 0; i < seqnum; i++)
  {
    // char *file_line = NULL;
    // char *i_file_line = NULL;
    // size_t i_seq_len_val = 0;
    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);
#pragma omp atomic write
    input_offset = len_val->line_offset; // to continue from this line next time
    size_t seq_len_val = len_val->len;

#pragma omp parallel for shared(i, seq_len_val) private(input_offset) num_threads(omp_get_max_threads()) // schedule(dynamic)
    for (int j = 0; j < seqnum; j++)
    {

      if (input_offset >= sb->st_size)
      {
        printf("Error: Could not find sequence length for sequence %d, %d\n", i, j);
        exit(1);
      }

      // fasta_len_value len_val = get_seqlen_mmapped_fasta(mmapped_fasta->mapped_file, &i, sb, input_offset);

      // size_t seq_len_val = len_val.len;

      // #pragma omp atomic write
      //       i_seq_len_val = i_seq_len_val + 5 + (seq_len_val * 2);

      // printf(" seq_len_val: %d \n", seq_len_val); // DEBUG
      // #pragma omp atomic write
      //       input_offset = len_val.line_offset; // to continue from this line next time

      char *file_line = (char *)malloc((seq_len_val + 1) * 2 * sizeof(char));
      strcpy(file_line, "");
      // fprintf(open_pos_f, ">%d,%d:", i, j);
      // sprintf(file_line, ">%d,%d:", i, j);
      for (int hv = 1; hv <= seq_len_val; hv++)
      {
        // fprintf(open_pos_f, "%d,", 1);
        strcat(file_line, "1,");
      }
      strcat(file_line, "1");
// fprintf(open_pos_f, "%d", 1);
// fprintf(open_pos_f, "\n");
#pragma omp critical
      fprintf(open_pos_f, ">%d,%d:%s\n", i, j, file_line);
      free(file_line);
    }
    // i_file_line = (char *)malloc((i_seq_len_val) * sizeof(char));
    // strcpy(i_file_line, file_line);
    // #pragma omp critical
    // fprintf(open_pos_f, "%s", i_file_line);
    // free(file_line);
    // free(i_file_line);
  }

  fclose(open_pos_f);
}

void create_overlapweights_file()
{
  // Writing 0's to overlap weights file to remove memory limitations. The file can be mmaped later
  FILE *ow_f;

  int seq_file_length = sizeof(seq_file) / sizeof(char);
  char ow_fname[seq_file_length]; // = seq_file;

  snprintf(ow_fname, sizeof(seq_file), "%s", seq_file);
  strcat(ow_fname, ".ow");

  ow_f = fopen(ow_fname, "w");

#pragma omp parallel for num_threads(omp_get_max_threads()) // schedule(dynamic)
  for (int i = 0; i < seqnum; i++)
  {
#pragma omp parallel for shared(i, seq_len_val) num_threads(omp_get_max_threads()) // schedule(dynamic)
    for (int j = 0; j < seqnum; j++)
    {
      char *file_line = (char *)malloc((4 + 1) * 2 * sizeof(char));
      strcpy(file_line, "");
      // for (int hv = 0; hv < (2 * 2) - 1; hv++) // 4 elements (0,0), (0,1), (1,0), (1,1)
      // {
      //   strcat(file_line, "0,");
      // }
      strcat(file_line, "0");
#pragma omp critical
      {
        fprintf(ow_f, ">%d,%d:%s\n", i, j, file_line);
      }
      free(file_line);
    }
  }
  fclose(ow_f);
}

main(int argc, char **argv)
{
  int k, anc1, dia_counter, tmpi1, tmpi2;

  // struct multi_frag *current_dia, *diagonal1, *diagonal2, *anc_dia;
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

  FILE *fp_ali, *fp2, *fp3, *fp4, *fp_log, *fp_fsm, *fp_st, *fp_csc;
  FILE *fp_matrix; /* file containing similarity matrix */

  strcpy(mat_name, MATNAME);
  strcpy(clust_sim, "av");

  par_file = (char *)calloc((size_t)NAME_LEN, sizeof(char));

  if (time_stamps)
    beg_ali = clock();

  strcpy(dialign_dir, "DIALIGN2_DIR");

  if ((par_file = getenv(dialign_dir)) == NULL)
  {
    printf("\n \n \n    Please set the environmentvariable DIALIGN2_DIR \n");
    printf("    as described in the README file \n");
    exit(1);
  }

  argnum = argc;

  strcpy(par_dir, par_file);

  if (argc == 1)
  {
    printf("\n    usage: %s [ options ] <seq_file> \n\n", argv[0]);
    printf("    <seq_file> contains input sequences in FASTA format.\n");
    printf("    Per default, sequences are assumedchar *test_name = get_seqname_mmapped_fasta(fasta, 0); to be protein sequences.\n");
    printf("    For DNA alignment, please use one of these options: \n\n");
    printf("     -n    DNA sequences; similarity calculated at the nucleotide level \n\n");
    printf("     -nt   DNA sequences; similarity calculated at the peptide level\n");
    printf("           (by translation using the genetic code) \n\n");
    printf("     -lgs  long genomic sequences: Both nucleotide and peptide\n");
    printf("           similarities calculated \n\n");
    printf("    Many more options are available, please consult the \n");
    printf("    DIALIGN USER_GUIDE that should come with the DIALIGN package.\n");
    printf("    For more information on DIALIGN, please visit the DIALIGN\n");
    printf("    home page at BiBiServ (Bielefeld Bioinformatic Server): \n\n");
    printf("        http://bibiserv.techfak.uni-bielefeld.de/dialign/ \n\n");
    exit(1);
  }

  arguments = (char **)calloc(argnum, sizeof(char *));

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < argnum; i++)
  {
    arguments[i] = (char *)calloc(NAME_LEN, sizeof(char));
    strcpy(arguments[i], argv[i]);
  }

  strcpy(input_name, argv[argc - 1]);

  threshold = 0.0;

  para_read(argnum, arguments);

  for (i = 0; i < argnum; i++)
  {
    free(arguments[i]);
  }
  free(arguments);

  if ((textual_alignment == 0) && (col_score == 1))
  {
    printf("\n\n   Option -csc makes sense only if \"textual alignment\"");
    printf(" is produced. \n");
    printf("   This can be enforced with option -ta \n\n");
    printf("   program terminated \n\n\n");
    exit(1);
  }

  if (cd_gobics)
  {
    strcpy(input_line, "program parameters:  ");
    for (i = 1; i < (argnum - 1); i++)
    {
      strcat(input_line, argv[i]);
      strcat(input_line, " ");
    }
  }
  else
  {
    strcpy(input_line, "program call:  ");
    for (i = 0; i < argnum; i++)
    {
      strcat(input_line, argv[i]);
      strcat(input_line, " ");
    }
  }

  if (wgt_type > 0)
    strict = 1;

  strcpy(seq_file, input_name);

  if (
      (!strcmp(input_name + strlen(input_name) - 4, ".seq")) || (!strcmp(input_name + strlen(input_name) - 3, ".fa")) || (!strcmp(input_name + strlen(input_name) - 6, ".fasta")))
    if ((char_ptr = strrchr(input_name, '.')) != NULL)
      *char_ptr = '\0';

  strcpy(anc_name, input_name);
  strcat(anc_name, ".anc");

  // strcpy(anc_name_sorted, input_name);
  // strcat(anc_name_sorted, ".sorted.anc");

  strcpy(ow_filename, seq_file);
  strcat(ow_filename, ".ow");

  mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  size_t file_size = mmapped_fasta->sb.st_size;

  mmapped_file *mmapped_anc = mmap_file(anc_name, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);
  size_t file_size_anc = mmapped_anc->sb.st_size;

  // mmapped_file *mmapped_diags;
  // mmapped_file *sorted_diags;

  int seq_file_length = sizeof(seq_file) / sizeof(char);

  char original_diag_fname[seq_file_length];
  strcpy(&original_diag_fname, &input_name);
  strcat(original_diag_fname, ".diags");
  char sorted_diag_fname[seq_file_length];
  strcpy(&sorted_diag_fname, &original_diag_fname);
  strcat(sorted_diag_fname, ".sorted");

  remove(original_diag_fname);
  remove(sorted_diag_fname);

  // fclose(fopen(diag_fname, "w"));
  // char pathtrace_fname[seq_file_length]; // = seq_file;
  // strcpy(&pathtrace_fname, &input_name);
  // strcat(pathtrace_fname, ".trace");
  // remove(pathtrace_fname);

  // exit(0);
  // int fd = open(seq_file, O_RDONLY);
  // if (fd == -1)
  // {
  //   perror("Error opening file");
  //   exit(1);
  // }

  // struct stat sb;
  // if (fstat(fd, &sb) == -1)
  // {
  //   perror("Error getting file size");
  //   exit(1);
  // }

  // // mmap the file into memory
  // char *fasta = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  // if (fasta == MAP_FAILED)
  // {
  //   perror("Error mmapping file");
  //   exit(1);
  // }

  // close(fd);

  seqnum = get_seqcount_mmapped_fasta(mmapped_fasta);
  printf(" seqnum:%d\n", seqnum); // DEBUG

  // create thread for open_pos file creation

  pthread_t open_pos_thread;
  // mmapped_file mmapped_fasta;
  // mmapped_fasta.mapped_file = fasta;
  // mmapped_fasta.sb = sb;

  int open_pos_thread_status = pthread_create(&open_pos_thread, NULL, create_openpos_file, mmapped_fasta);

  pthread_t ow_thread;
  int ow_thread_status = pthread_create(&ow_thread, NULL, create_overlapweights_file, NULL);

  // // DEBUG
  // pthread_join(ow_thread, NULL);
  // mmapped_file *mmapped_ow = mmap_file(ow_filename, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);
  // size_t file_size_ow = mmapped_ow->sb.st_size;

  // int val2 = 2;
  // int val1 = 1;
  // int val0 = 0;
  // int val3 = 3;

  // double vald1 = 0.45791;

  // double vald2 = 50.12;

  // printf("ow set value [%d,%d,%d]:%f \n", val0, val1, val0, set_ow_mmapped(mmapped_ow, &val0, &val1, &val0, &vald1, 0).value);
  // printf("ow get value [%d,%d,%d]:%g \n", val0, val1, val0, get_ow_mmapped(mmapped_ow, &val0, &val1, &val0, 0).value);
  // printf("ow get value [%d,%d,%d]:%g \n", val0, val2, val0, get_ow_mmapped(mmapped_ow, &val0, &val2, &val0, 0).value);
  // printf("ow set value [%d,%d,%d]:%f \n", val0, val2, val0, set_ow_mmapped(mmapped_ow, &val0, &val2, &val0, &vald1, 0).value);
  // printf("ow get value [%d,%d,%d]:%g \n", val0, val2, val0, get_ow_mmapped(mmapped_ow, &val0, &val2, &val0, 0).value);
  // printf("ow set value [%d,%d,%d]:%f \n", val0, val2, val0, set_ow_mmapped(mmapped_ow, &val0, &val2, &val0, &vald2, 0).value);
  // printf("ow get value [%d,%d,%d]:%g \n", val0, val2, val0, get_ow_mmapped(mmapped_ow, &val0, &val2, &val0, 0).value);
  // msync(mmapped_ow->mapped_file, mmapped_ow->sb.st_size, MS_SYNC);
  // munmap(mmapped_ow->mapped_file, mmapped_ow->sb.st_size);
  // exit(0);

  //   int val1 = 0;
  // int* test_len = get_seqlen_mmapped_fasta(mmapped_fasta->mapped_file, &val1, sb);
  // printf(" test_len:%d\n", *test_len);
  // exit(0);

  // printf(" here 8.1\n"); // DEBUG

  // pthread_join(open_pos_thread, NULL);
  // int fd1 = open("tmp.fa.op", O_RDWR);
  // if (fd1 == -1)
  // {
  //   perror("Error opening file");
  //   exit(1);
  // }
  // struct stat sb1;
  // if (fstat(fd1, &sb1) == -1)
  // {
  //   perror("Error getting file size");
  //   exit(1);
  // }
  // char *op_f1 = mmap(NULL, sb1.st_size, PROT_WRITE, MAP_SHARED, fd1, 0);
  // if (op_f1 == MAP_FAILED)
  // {
  //   perror("Error mmapping file");
  //   exit(1);
  // }

  // int val2 = 2;
  // int val1 = 1;
  // int val0 = 0;
  // int val3 = 3;
  // int val4 = 4;
  // int val5 = 5;

  // // printf("i:%d j:%d open_pos value:%d \n", val1, val3, test_op1);                            // DEBUG
  // printf("open_pos get value [%d,%d,%d]:%d \n", val1, val2, val5, get_openpos_mmapped(op_f1, &sb1, &val1, &val2, &val5, 0).value);        // DEBUG must be == 1
  // printf("open_pos set value [%d,%d,%d]:%d \n", val1, val2, val5, set_openpos_mmapped(op_f1, &sb1, &val1, &val2, &val5, &val5, 0).value); // DEBUG must be == 1
  // printf("open_pos get value [%d,%d,%d]:%d \n", val1, val2, val5, get_openpos_mmapped(op_f1, &sb1, &val1, &val2, &val5, 0).value);        // DEBUG must be == 1
  // printf("open_pos get value [%d,%d,%d]:%d \n", val0, val2, val5, get_openpos_mmapped(op_f1, &sb1, &val0, &val2, &val5, 0).value);        // DEBUG must be == 1
  // printf("open_pos set value [%d,%d,%d]:%d \n", val0, val2, val5, set_openpos_mmapped(op_f1, &sb1, &val0, &val2, &val5, &val0, 0).value); // DEBUG must be == 0
  // printf("open_pos get value [%d,%d,%d]:%d \n", val0, val0, val4, get_openpos_mmapped(op_f1, &sb1, &val0, &val0, &val4, 0).value);        // DEBUG must be == 1
  // printf("open_pos set value [%d,%d,%d]:%d \n", val0, val0, val4, set_openpos_mmapped(op_f1, &sb1, &val0, &val0, &val4, &val2, 0).value); // DEBUG must be == 2
  // printf("open_pos get value [%d,%d,%d]:%d \n", val0, val0, val4, get_openpos_mmapped(op_f1, &sb1, &val0, &val0, &val4, 0).value);        // DEBUG must be == 2
  // printf("open_pos get value [%d,%d,%d]:%d \n", val0, val2, val2, get_openpos_mmapped(op_f1, &sb1, &val0, &val2, &val2, 0).value);        // DEBUG must be == 1
  // msync(op_f1, sb1.st_size, MS_SYNC);
  // munmap(op_f1, sb1.st_size);
  // munmap(mmapped_fasta->mapped_file, sb1.st_size);
  // exit(0);

  if (motifs)
    regex_parse(mot_regex);

  if ((seqnum == 2) && (iter_cond_prob == 0))
    max_itnum = 1;

  if ((ow_force == 0) && (seqnum > OVERLAP_THRESHOLD))
    overlap_weights = 0;
  if (seqnum == 2)
    overlap_weights = 0;

  if (seqnum < 2)
  {

    if (cd_gobics)
    {
      printf("\n\n         Something is wrong with your sequence file. Maybe you entered a\n");
      printf("         MS WORD or RFT file or your file contains only one single sequence.\n");
      printf("         Please note that our server only accepts plain text files. \n\n");
      printf("         For more information, please consult our online manual \n");
      printf("         at the CHAOS/DIALIGN home page:\n\n");
      printf("             http://dialign.gobics.de/chaos-dialign-manual");
    }

    else
    {
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

  // printf("\n here 9 \n"); // DEBUG

  // if( (pair_score = (float **) calloc( seqnum , sizeof(float *) )) == NULL)
  if ((pair_score = (double **)malloc(seqnum * sizeof(double *))) == NULL)
  {
    printf(" problems with memory allocation for `pair_score' !  \n \n");
    exit(1);
  }

  // printf("\n here 10 \n"); // DEBUG

  for (i = 0; i < seqnum; i++)
    // if( (pair_score[i] = (float *) calloc( seqnum , sizeof(float) )) == NULL)
    if ((pair_score[i] = (double *)malloc(seqnum * sizeof(double))) == NULL)
    {
      printf(" problems with memory allocation for `pair_score' !  \n \n");
      exit(1);
    }

  // printf("\n here 11 \n"); // DEBUG

  // if(( cont_it_p = (short **) calloc( seqnum , sizeof( short *))) == NULL )
  if ((cont_it_p = (short **)malloc(seqnum * sizeof(short *))) == NULL)
  {
    printf(" problems with memory allocation for `cont_it_p ' !  \n \n");
    exit(1);
  }

  // printf("\n here 12 \n"); // DEBUG

  for (i = 0; i < seqnum; i++)
    // if( (cont_it_p[i] = (short *) calloc( seqnum , sizeof(short) )) == NULL)
    if ((cont_it_p[i] = (short *)malloc(seqnum * sizeof(short))) == NULL)
    {
      printf(" problems with memory allocation for `cont_it_p' !  \n \n");
      exit(1);
    }

  printf("\n here 13 \n"); // DEBUG

#pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
    for (j = 0; j < seqnum; j++)
      cont_it_p[i][j] = 1;

  printf("\n here 14 \n"); // DEBUG

  // if ((glob_sim =
  //          (double **)malloc(seqnum * sizeof(double *))) == NULL)
  // {
  //   printf("Problems with memory allocation for glob_sim\n");
  //   exit(1);
  // }

  size_t input_offset = 0;

  for (i = 0; i < seqnum; i++)
  {
    // printf("\n here 14.1 : %d \n", i); //DEBUG
    if (input_offset >= file_size)
    {
      printf("Error: Could not find sequence length for sequence %d\n", i);
    }

    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);

    int seq_len_val = len_val->len;
    input_offset = len_val->line_offset + seq_len_val + 1;
    // av_len = av_len + seqlen[i];
    // printf("\n here 14.1 :: i:%d num:%d :: seq_len: %d :: ofsset:%d file_size:%d \n", i, len_val.seq_num, seq_len_val, input_offset, sb->st_size); // DEBUG
    av_len = av_len + (double)seq_len_val;

    // printf("\n here 14.2 \n"); //DEBUG

    // if( seqlen[i] == 0 )
    if (seq_len_val == 0)
    {
      printf("\n \n \n                       WARNING: \n \n");
      printf("          Sequence %d contains no residues.\n", i + 1);
      printf("          Please inspect the sequence file.\n \n ");
      printf("\n \n          Program terminated \n \n \n ");

      exit(1);
    }

    // if ((glob_sim[i] =
    //          (double *)malloc(seqnum * sizeof(double))) == NULL)
    // {
    //   printf("Problems with memory allocation for glob_sim \n");
    //   exit(1);
    // }

    // printf("\n here 14.3 \n"); //DEBUG

    // if(maxlen < seqlen[i])
    //    maxlen = seqlen[i];
    if (maxlen < seq_len_val)
      maxlen = seq_len_val;
  }

  printf("\n here 14.4 \n"); // DEBUG

  av_len = av_len / seqnum;

  if (motifs)
    seq_parse(mot_regex);

  // seq_shift(); //Why do they do this?

  // if ((glob_sim =
  //      (float **)malloc(seqnum * sizeof(float *))) == NULL)
  // {
  //   printf("Problems with memory allocation for glob_sim\n");
  //   exit(1);
  // }

  // for (i = 0; i < seqnum; i++)
  // {

  //   if ((glob_sim[i] =
  //        (float *)malloc(seqnum * sizeof(float))) == NULL)
  //   {
  //     printf("Problems with memory allocation for glob_sim \n");
  //     exit(1);
  //   }
  // }

  strcpy(par_str, "sdfsdf");

  if (argc > 1)
  {
    strcpy(str, par_dir);
    strcat(str, "/");
    strcat(str, mat_name);
    strcpy(mat_name_p, str);

    if ((fp_matrix = fopen(mat_name_p, "r")) == NULL)
    {

      printf("\n\n Cannot find the file %s \n\n", mat_name);
      printf(" Make sure the environment variable DIALIGN2_DIR points\n");
      printf(" to a directory containing the files \n\n");
      printf("   BLOSUM \n   tp400_dna\n   tp400_prot \n   tp400_trans \n\n");
      printf(" These files should be contained in the DIALIGN package \n\n\n");
      exit(1);

      printf("\n \n \n \n              ATTENTION ! \n \n");
      printf("\n   There is no similarity matrix `%s'. \n", mat_name);
      printf("   in the directory \n \n");
      printf("           %s\n \n", par_dir);
      exit(1);
    }
  }

  if (wgt_type != 1)
    matrix_read(fp_matrix);

  mem_alloc();

  printf("\n here 14.5 \n"); // DEBUG

  if (wgt_type != 1)
  { // if( (amino = (int **) calloc( seqnum , sizeof(int *) ) ) == NULL)
    if ((amino = (int **)malloc(seqnum * sizeof(int *))) == NULL)
    {
      printf(" problems with memory allocation");
      printf(" for `amino' !  \n \n");
      exit(1);
    }
  }

  if (crick_strand)
  {
    // if( (amino_c = (int **) calloc( seqnum , sizeof(int *) ) ) == NULL)
    if ((amino_c = (int **)malloc(seqnum * sizeof(int *))) == NULL)
    {
      printf(" problems with memory allocation");
      printf(" for `amino_c' !  \n \n");
      exit(1);
    }
  }

  input_offset = 0;

  for (i = 0; i < seqnum; i++)
  {
    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);
    size_t seq_len_val = len_val->len;
    input_offset = len_val->line_offset;

    // size_t seq_len_val = get_seqlen_mmapped_fasta(mmapped_fasta->mapped_file, &i, sb);
    if (wgt_type != 1)
    {
      // if ((amino[i] = (int *)malloc((seqlen[i] + 5) * sizeof(int))) == NULL)
      if ((amino[i] = (int *)malloc((seq_len_val + 5) * sizeof(int))) == NULL)
      {
        printf(" problems with memory allocation");
        printf(" for `amino[%d]' !  \n \n", i);
        exit(1);
      }
    }
    if (crick_strand)
    {
      // if ((amino_c[i] = (int *)malloc((seqlen[i] + 5) * sizeof(int))) == NULL)
      if ((amino_c[i] = (int *)malloc((seq_len_val + 5) * sizeof(int))) == NULL)
      {
        printf(" problems with memory allocation");
        printf(" for `amino_c[%d]' !  \n \n", i);
        exit(1);
      }
    }
  }

  printf("\n here 14.5 \n"); // DEBUG

  /******************************************************
   *                                                     *
   *  read file, that contains data of anchored regions  *
   *                                                     *
   ******************************************************/

  if (anchors)
  {
    // multi_anc_read(input_name);
    //  multi_anc_read_mmap(input_name);
    total_anc_num = get_anchor_count(input_name);
  }

  printf("\n here 15 \n"); // DEBUG

  if (exclude_frg)
  {

    // if( ( exclude_list = (int ***) calloc( seqnum , sizeof(int **) )) == NULL)
    if ((exclude_list = (int ***)malloc(seqnum * sizeof(int **))) == NULL)
    {
      printf(" problems with memory allocation for 'exclude_list' \n \n");
      exit(1);
    }

    for (i = 0; i < seqnum; i++)
      // if( ( exclude_list[ i ] = (int **) calloc( seqnum , sizeof(int *) )) == NULL)
      if ((exclude_list[i] = (int **)malloc(seqnum * sizeof(int *))) == NULL)
      {
        printf(" problems with memory allocation for 'exclude_list' \n \n");
        exit(1);
      }

    input_offset = 0;

    for (i = 0; i < seqnum; i++)
      for (j = 0; j < seqnum; j++)
      {
        fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);
        size_t seq_len_val = len_val->len;
        input_offset = len_val->line_offset;

        // if ((exclude_list[i][j] = (int *)malloc((seqlen[i] + 1) * sizeof(int))) == NULL)
        if ((exclude_list[i][j] = (int *)malloc((seq_len_val + 1) * sizeof(int))) == NULL)
        {
          printf(" problems with memory allocation for 'exclude_list' \n \n");
          exit(1);
        }
      }

    exclude_frg_read(input_name, exclude_list); // TODO: CHECK this
  }

  if (wgt_type == 0)
    tp400_read(0, tp400_prot);
  if (wgt_type % 2)
    tp400_read(1, tp400_dna);
  if (wgt_type > 1)
    tp400_read(2, tp400_trans);

  /****************************\
  *                            *
  *    Name of output files    *
  *                            *
  \****************************/

  if (default_name)
  {
    strcpy(printname, input_name);
    strcpy(prn, printname);
  }
  else
  {
    strcpy(printname, output_name);
    strcpy(prn, printname);
  }

  strcpy(prn2, prn);

  if (default_name)
    strcat(prn, ".ali");

  strcat(prn2, ".fa");

  printf("\n here 15.1 \n"); // DEBUG

  strcpy(logname, printname);
  strcat(logname, ".log");

  strcpy(fsm_name, printname);
  strcat(fsm_name, ".fsm");

  if (print_status)
  {
    strcpy(pst_name, printname);
    strcat(pst_name, ".sta");
  }

  if (afc_file)
  {
    strcpy(dia_name, printname);
    strcat(dia_name, ".afc");
    fp_dia = fopen(dia_name, "w");
    fprintf(fp_dia, "\n #  %s \n\n  seq_len: ", input_line);

    input_offset = 0;

    for (i = 0; i < seqnum; i++)
    {
      // fprintf(fp_dia, "  %d ", seqlen[i]);
      fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);
      size_t seq_len_val = len_val->len;
      input_offset = len_val->line_offset;

      fprintf(fp_dia, "  %d ", seq_len_val);
    }
    fprintf(fp_dia, "\n\n");
  }

  if (col_score)
  {
    strcpy(csc_name, printname);
    strcat(csc_name, ".csc");
    fp_csc = fopen(csc_name, "w");
  }

  if (dia_pa_file)
  {
    strcpy(dia_pa_name, printname);
    strcat(dia_pa_name, ".fop");

    fp_dpa = fopen(dia_pa_name, "w");

    fprintf(fp_dpa, "\n #  %s \n\n  seq_len: ", input_line);
    input_offset = 0;
    for (i = 0; i < seqnum; i++)
    {
      // fprintf(fp_dpa, "  %d ", seqlen[i]);
      fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);
      size_t seq_len_val = len_val->len;
      input_offset = len_val->line_offset;

      fprintf(fp_dpa, "  %d ", seq_len_val);
    }
    fprintf(fp_dpa, "\n\n");
    fclose(fp_dpa);
  }

  printf("\n here 15.2 \n"); // DEBUG

  if (motifs)
  {
    strcpy(mot_file_name, printname);
    strcat(mot_file_name, ".mot");
    fp_mot = fopen(mot_file_name, "w");

    fprintf(fp_mot, "\n #  %s \n\n   ", input_line);
    fprintf(fp_mot, " motif: %s \n\n", mot_regex);
    fprintf(fp_mot, " max offset for motifs = %d \n\n", (int)max_mot_offset);
    fprintf(fp_mot, " the following fragments contain the motif: \n\n");
    fprintf(fp_mot, "   seq1 seq2    beg1 beg1 len    wgt");
    fprintf(fp_mot, "   # mot    mot_wgt  \n\n");
  }

  printf("\n here 15.3 \n"); // DEBUG

  if (frag_file)
  {
    strcpy(frag_file_name, printname);
    strcat(frag_file_name, ".frg");
    fp_frg = fopen(frag_file_name, "w");

    fprintf(fp_frg, "\n #  %s \n\n  seq_len: ", input_line);
    input_offset = 0;
    for (i = 0; i < seqnum; i++)
    {
      // fprintf(fp_frg, "  %d ", seqlen[i]);
      fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset);
      size_t seq_len_val = len_val->len;
      input_offset = len_val->line_offset;

      fprintf(fp_frg, "  %d ", seq_len_val);
    }
    fprintf(fp_frg, "\n  sequences: ");

    input_offset = 0;
    for (i = 0; i < seqnum; i++)
    {
      // fprintf(fp_frg, "  %s ", seq_name[i]);
      fasta_value *name_val = get_seqname_mmapped_fasta(mmapped_fasta, &i, input_offset);
      // name_val->data++;
      char *seq_name_val = name_val->data;
      input_offset = name_val->line_offset;

      fprintf(fp_frg, "  %s ", seq_name_val);
      free(seq_name_val);
      // free(name_val->data);
      free(name_val);
    }

    fprintf(fp_frg, "\n\n");
  }

  printf("\n here 15.4 \n"); // DEBUG

  clos = newAligGraphClosure(seqnum, 0, NULL, mmapped_fasta);
  // clos = newAligGraphClosure(seqnum, seqlen, 0, NULL);

  //   // if( (open_pos = (int *** ) calloc( seqnum , sizeof(int **))) == NULL)
  //   if( (open_pos = (int *** ) malloc( seqnum * sizeof(int **))) == NULL)
  //      {
  //        printf("Problems with memory allocation for open_pos\n");
  //        exit(1);
  //      }

  // printf("\n here 15.6 \n"); //DEBUG

  //   for(i=0;i<seqnum;i++)
  //       {
  //         if( (open_pos[i] =
  //         // (int ** ) calloc( seqnum , sizeof(int *))) == NULL)
  //         (int ** ) malloc( seqnum * sizeof(int *))) == NULL)
  //           {
  //              printf("Problems with memory allocation for open_pos\n");
  //              exit(1);
  //           }
  //       }

  // printf("\n here 15.7 \n"); //DEBUG

  //   for(i=0;i<seqnum;i++)
  //   for(j=0;j<seqnum;j++)
  //       {
  //         printf("\n here 15.7.1 : %d  %d \n", i , j); //DEBUG
  //        if( (open_pos[i][j] =
  //       //  (int * ) calloc( ( seqlen[i]+2) , sizeof(int) ) ) == NULL)
  //        (int * ) malloc( ( seqlen[i]+2) * sizeof(int) ) ) == NULL)
  //           {
  //              printf("Problems with memory allocation for open_pos\n");
  //              exit(1);
  //           }
  //       }

  // printf("\n here 15.8 \n"); //DEBUG

  //   for( i = 0 ; i <seqnum ; i++)
  //   for( j = 0 ; j <seqnum ; j++)
  //   for( hv = 1 ; hv <= seqlen[i] ; hv++)
  //      open_pos[i][j][hv] = 1;

  //   // Writing 1's to file to remove memory limitations. The file can be mmaped later
  //   //  #UNCOMMENT LATER //DEBUG
  //   FILE *open_pos_f;
  //   // char *file_line = NULL;
  //   // size_t seq_len_val = -1;
  //   int seq_file_length = sizeof(seq_file) / sizeof(char);
  //   char open_pos_fname[seq_file_length]; // = seq_file;
  //   // printf("\n seq_file: %s \n", seq_file); // DEBUG
  //   // strcpy(open_pos_fname, seq_file);
  //   snprintf(open_pos_fname, sizeof(seq_file), "%s", seq_file);
  //   strcat(open_pos_fname, ".op");
  //   // printf("\n here 15.5.5 \n"); // DEBUG
  //   open_pos_f = fopen(open_pos_fname, "w");
  // #pragma omp parallel for num_threads(omp_get_max_threads()) // schedule(static) // private(file_line, seq_len_val)
  //   for (i = 0; i < seqnum; i++)
  //   {
  // #pragma omp parallel for shared(i) num_threads(omp_get_max_threads()) // schedule(static) // private(file_line, seq_len_val)
  //     for (j = 0; j < seqnum; j++)
  //     {
  //       // fprintf(open_pos_f, ">%d,%d:", i, j);
  //       // for (hv = 1; hv <= seqlen[i]; hv++)
  //       size_t seq_len_val = get_seqlen_mmapped_fasta(mmapped_fasta->mapped_file, &i, sb);
  //       char *file_line = (char *)malloc(((seq_len_val + 1) * 2) * sizeof(char));
  //       strcpy(file_line, "");
  //       for (hv = 1; hv < seq_len_val; hv++)
  //       {
  //         //  open_pos[i][j][hv] = 1;
  //         // fprintf(open_pos_f, "%d\t%d\n", hv, 1);
  //         // fprintf(open_pos_f, "%d,", 1);
  //         strcat(file_line, "1,");
  //       }
  //       // fprintf(open_pos_f, "%d", 1);
  //       strcat(file_line, "1");
  // // printf("\n file_line:%s \n", file_line); // DEBUG
  // #pragma omp critical
  //       {
  //         fprintf(open_pos_f, ">%d,%d:%s\n", i, j, file_line);
  //       }
  //       // fprintf(open_pos_f, "\n");
  //       free(file_line);
  //     }
  //   }
  //   // free(file_line);
  //   fclose(open_pos_f);
  //   // exit(0); // DEBUG
  /**************************************
   *                                     *
   *      definition of  `amino'         *
   *                                     *
   **************************************/

  printf("\n here 15.9 \n"); // DEBUG

  input_offset = 0;

  if (wgt_type > 1)
    for (hv = 0; hv < seqnum; hv++)
      for (i = 1; i <= seqlen[hv] - 2; i++)
      {
        fasta_value *seq_val = get_seq_mmapped_fasta(mmapped_fasta, &hv, input_offset);
        // seq_val->data++;
        char *seq_hv = seq_val->data;

        input_offset = seq_val->line_offset + strlen(seq_hv) + 1;

        // if( translate( seq[hv][i],seq[hv][i+1],seq[hv][i+2],hv,i ) == -1)
        if (translate(seq_hv[i], seq_hv[i + 1], seq_hv[i + 2], hv, i) == -1)
          exit(1);

        // amino[hv][i] = translate( seq[hv][i],seq[hv][i+1],seq[hv][i+2],hv,i);
        amino[hv][i] = translate(seq_hv[i], seq_hv[i + 1], seq_hv[i + 2], hv, i);

        if (crick_strand)
        {
          // nuc1 = invert( seq[hv][i+2] );
          // nuc2 = invert( seq[hv][i+1] );
          // nuc3 = invert( seq[hv][i] );
          nuc1 = invert(seq_hv[i + 2]);
          nuc2 = invert(seq_hv[i + 1]);
          nuc3 = invert(seq_hv[i]);

          amino_c[hv][i] = translate(nuc1, nuc2, nuc3, hv, i);
        }

        free(seq_hv);
        // free(seq_val->data);
        free(seq_val);
      }

  // printf("\n here 15.6 \n"); // DEBUG

  input_offset = 0;

  if (wgt_type == 0)
    for (hv = 0; hv < seqnum; hv++)
    {
      fasta_value *seq_val = get_seq_mmapped_fasta(mmapped_fasta, &hv, input_offset);
      // seq_val->data++;
      char *seq_hv = seq_val->data;
      fasta_len_value *len_val = get_seqlen_mmapped_fasta(mmapped_fasta, &hv, input_offset);
      size_t seq_len_hv = len_val->len;

      input_offset = seq_val->line_offset + strlen(seq_hv) + 1;

      // for (i = 1; i <= seqlen[hv]; i++)
      for (i = 1; i <= seq_len_hv; i++)
      {
        // if( seq[hv][i] == 'C' ) amino[hv][i] = 1;
        // if( seq[hv][i] == 'S' ) amino[hv][i] = 2;
        // if( seq[hv][i] == 'T' ) amino[hv][i] = 3;
        // if( seq[hv][i] == 'P' ) amino[hv][i] = 4;
        // if( seq[hv][i] == 'A' ) amino[hv][i] = 5;
        // if( seq[hv][i] == 'G' ) amino[hv][i] = 6;
        // if( seq[hv][i] == 'N' ) amino[hv][i] = 7;
        // if( seq[hv][i] == 'D' ) amino[hv][i] = 8;
        // if( seq[hv][i] == 'E' ) amino[hv][i] = 9;
        // if( seq[hv][i] == 'Q' ) amino[hv][i] = 10;
        // if( seq[hv][i] == 'H' ) amino[hv][i] = 11;
        // if( seq[hv][i] == 'R' ) amino[hv][i] = 12;
        // if( seq[hv][i] == 'K' ) amino[hv][i] = 13;
        // if( seq[hv][i] == 'M' ) amino[hv][i] = 14;
        // if( seq[hv][i] == 'I' ) amino[hv][i] = 15;
        // if( seq[hv][i] == 'L' ) amino[hv][i] = 16;
        // if( seq[hv][i] == 'V' ) amino[hv][i] = 17;
        // if( seq[hv][i] == 'F' ) amino[hv][i] = 18;
        // if( seq[hv][i] == 'Y' ) amino[hv][i] = 19;
        // if( seq[hv][i] == 'W' ) amino[hv][i] = 20;

        if (seq_hv[i] == 'C')
          amino[hv][i] = 1;
        if (seq_hv[i] == 'S')
          amino[hv][i] = 2;
        if (seq_hv[i] == 'T')
          amino[hv][i] = 3;
        if (seq_hv[i] == 'P')
          amino[hv][i] = 4;
        if (seq_hv[i] == 'A')
          amino[hv][i] = 5;
        if (seq_hv[i] == 'G')
          amino[hv][i] = 6;
        if (seq_hv[i] == 'N')
          amino[hv][i] = 7;
        if (seq_hv[i] == 'D')
          amino[hv][i] = 8;
        if (seq_hv[i] == 'E')
          amino[hv][i] = 9;
        if (seq_hv[i] == 'Q')
          amino[hv][i] = 10;
        if (seq_hv[i] == 'H')
          amino[hv][i] = 11;
        if (seq_hv[i] == 'R')
          amino[hv][i] = 12;
        if (seq_hv[i] == 'K')
          amino[hv][i] = 13;
        if (seq_hv[i] == 'M')
          amino[hv][i] = 14;
        if (seq_hv[i] == 'I')
          amino[hv][i] = 15;
        if (seq_hv[i] == 'L')
          amino[hv][i] = 16;
        if (seq_hv[i] == 'V')
          amino[hv][i] = 17;
        if (seq_hv[i] == 'F')
          amino[hv][i] = 18;
        if (seq_hv[i] == 'Y')
          amino[hv][i] = 19;
        if (seq_hv[i] == 'W')
          amino[hv][i] = 20;
      }

      free(seq_hv);
      // free(seq_val->data);
      free(seq_val);
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

  num_dia_anc = total_anc_num * (seqnum - 1);

  // printf("\n here 16 \n"); // DEBUG

  pthread_join(ow_thread, NULL);
  // exit(0); // DEBUG

  // // mmapped_file *mmapped_ow = mmap_file(ow_filename, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);
  mmapped_file *mmapped_ow = mmap_file(ow_filename, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);
  size_t file_size_ow = mmapped_ow->sb.st_size;

  // wait or join open_pos file creation thread here
  pthread_join(open_pos_thread, NULL);
  // exit(0);

  if (anchors)
  {
    if (time_stamps)
      beg_ts = clock();

    if (nas == 0)
      if (bubblesort)
      { // bubble_sort(total_anc_num, anchor_frg);
        // bubble_sort_mmap(total_anc_num, mmapped_anc, mmapped_fasta, mmapped_ow);
        printf("PLACEHOLDER FOR bubble_sort_mmap()\n"); // DEBUG
        exit(0);                                        // DEBUG
      }
      else
      {
      // frag_sort(total_anc_num, anchor_frg, 0);
      //         i = 1;
      //         struct multi_frag *dp = get_anchor_by_num(mmapped_fasta, mmapped_anc, mmapped_ow, &i, 0, 0);
      // #pragma omp critical
      //         {
      //           frag_sort_mmap(total_anc_num, 0, mmapped_fasta, mmapped_ow, dp);
      //         }
      //         free(dp);

#pragma omp critical
        {
          frag_sort_mmap(total_anc_num, 0, mmapped_fasta, mmapped_anc, NULL);
        }
      }

    if (time_stamps)
    {
      end_ts = clock();
      time_diff_srt = (double)(end_ts - beg_ts) / CLOCKS_PER_SEC;
      if (time_stamps)
        printf(" for anc: time_diff_srt = %f \n", time_diff_srt);
    }

    // filter(&total_anc_num, anchor_frg);

    // filter_all(&total_anc_num);

    // munmap(mmapped_anc->mapped_file, mmapped_anc->sb.st_size);
    // mmapped_anc = mmap_file(anc_name_sorted, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);
    // file_size_anc = mmapped_anc->sb.st_size;

    // int first_anc = 1;
    // struct multi_frag *first_diag = get_anchor_by_num(mmapped_fasta, mmapped_anc, &first_anc, 0, 0);
    filter(&total_anc_num, 1, mmapped_anc, mmapped_fasta, NULL, mmapped_ow);

    /*  exit(1) ;
     */
  }

  // printf("\n here 17 \n"); // DEBUG

  if (long_output)
  {
    fp_log = fopen(logname, "w");
    fprintf(fp_log, "\n #  %s \n\n   ", input_line);
  }

  if (frg_mult_file)
  {
    fp_fsm = fopen(fsm_name, "w");
    fprintf(fp_fsm, "\n #  %s \n\n", input_line);
  }

  // printf("\n here 17.1 \n"); // DEBUG

  // if (
  //     (num_dia_bf = (int *)calloc((max_itnum + 1), sizeof(int))) == NULL)
  if (
      (num_dia_bf = (int *)malloc((max_itnum + 1) * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `num_dia_bf' !  \n \n");
    exit(1);
  }

  // printf("\n here 17.2 \n"); // DEBUG

  // if (
  //     (num_dia_af = (int *)calloc((max_itnum + 1), sizeof(int))) == NULL)
  if (
      (num_dia_af = (int *)malloc((max_itnum + 1) * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `num_dia_af' !  \n \n");
    exit(1);
  }

  // printf("\n here 17.3 \n"); // DEBUG

  // all_it_dia = (struct multi_frag *)calloc(1, sizeof(struct multi_frag));
  // all_it_dia = (struct multi_frag *)malloc(1 * sizeof(struct multi_frag));
  // current_dia = all_it_dia;

  strcpy(itname, printname);
  strcpy(itname2, printname);
  strcpy(itname3, printname);
  strcpy(itname4, printname);
  sprintf(str, ".ali");

  if (default_name)
    strcat(itname, str);

  sprintf(str, ".fa.aln");
  strcat(itname2, str);

  if (msf_file)
    strcat(itname3, ".ms");

  if (cw_file)
    strcat(itname4, ".cw");

  if (textual_alignment)
    fp_ali = fopen(itname, "w");

  if (standard_out)
    fp_ali = stdout;

  if (textual_alignment)
    if (fasta_file)
      fp2 = fopen(itname2, "w");

  if (msf_file)
    fp3 = fopen(itname3, "w");

  if (cw_file)
    fp4 = fopen(itname4, "w");

  if (textual_alignment)
    // para_print(seq_file, fp_ali);
    para_print_mmap(seq_file, fp_ali, mmapped_fasta);

  /***************************\
  *                           *
  *      ITERATION START      *
  *                           *
  \***************************/

  // FILE *diag_f;
  // diag_f = fopen(original_diag_fname, "w");

  istep = 0;
  while ((cont_it == 1) && (istep < max_itnum))
  {
    int total_diagonals = 0;

    // remove files if it is not the last iteration
    // remove(diag_fname);

    cont_it = 0;
    istep++;

    /* printf("\n  istep = %d \n", istep ); */

    // this_it_dia = current_dia;

    strcpy(itname, printname);
    strcpy(itname2, printname);
    strcpy(itname3, printname);
    strcpy(itname4, printname);
    sprintf(str, ".ali");

    if (default_name)
      strcat(itname, str);

    sprintf(str, ".fa");
    strcat(itname2, str);

    if (msf_file)
      strcat(itname3, ".ms");

    if (cw_file)
      strcat(itname4, ".cw");

    weight_sum_af = 0;
    num_dia_bf[istep] = 0;

    if (time_stamps)
      beg_pa = clock();

    if (ref_seq == 0)
      i_max = seqnum;
    else
      i_max = 1;

    // printf("\n here 17.4 \n"); // DEBUG

    size_t input_offset_i = 0;
    size_t input_offset_j = 0;
    // #pragma omp parallel for num_threads(omp_get_max_threads())
    for (i = 0; i < i_max; i++)
    {

      printf(" i_max:%d i:%d\n", i_max, i); // DEBUG
      fasta_len_value *seqlen_n1_val = get_seqlen_mmapped_fasta(mmapped_fasta, &i, input_offset_i);
      int seqlen_n1 = seqlen_n1_val->len;
      fasta_value *seqname_n1_val = get_seqname_mmapped_fasta(mmapped_fasta, &i, input_offset_i);
      // seqname_n1_val->data++;
      int seqname_n1_len = strlen(seqname_n1_val->data);
      // char seqname_n1[seqname_n1_len];
      char seqname_n1[seqname_n1_len + 1];
      strcpy(seqname_n1, seqname_n1_val->data);
      input_offset_i = seqname_n1_val->line_offset + strlen(seqname_n1) + 1;
      fasta_value *seq1_val = get_seq_mmapped_fasta(mmapped_fasta, &i, input_offset_i);
      char seq_n1[seqlen_n1 + 2];
      // printf("DEBUG seq_val:%s seqname_n1_val->data:%s seqlen_n1:%d\n", seq1_val->data, seqname_n1_val->data, seqlen_n1); // DEBUG
      // seq1_val->data++;
      strcpy(seq_n1, seq1_val->data);
#pragma omp atomic write
      input_offset_i = seq1_val->line_offset + strlen(seq1_val) + 1;
#pragma omp atomic write
      input_offset_j = input_offset_i;
// #pragma omp atomic write
//       *seq_n1 = seq_n1[1];
#pragma omp parallel for shared(i, seqlen_n1, seq_n1, seqname_n1) private(input_offset) num_threads(omp_get_max_threads())
      for (j = i + 1; j < seqnum; j++)
      {

        /****************************************\
        *                                        *
        *          PAIRWISE  ALIGNMENT           *
        *                                        *
        \****************************************/

        if (cont_it_p[i][j])
        {

          /*
                            printf("\n out of frc it %d : wgt 20 = %f \n", istep ,  wgt_dna[ 20 ][ 20 ] ) ;
          */

          // printf("\n here 17.5 \n"); //DEBUG

          // score = frag_chain( i , j , fp_ali, fp_mot, &num_dia_p );

          if (input_offset_i >= mmapped_fasta->sb.st_size)
          {
            printf("Coud not fetch fragment chain value for i = %d, j = %d offset:%ld mmapped_fasta->sb.st_size:%ld \n", i, j, input_offset_i, mmapped_fasta->sb.st_size);
            exit(1);
            // input_offset = 0;
          }

          frag_chain_value *frag_chain_val = frag_chain_mmap(i, j, fp_ali, fp_mot, &num_dia_p, mmapped_fasta, input_offset_j, seqlen_n1, &seq_n1, &seqname_n1);
          score = frag_chain_val->value;

#pragma omp atomic write
          input_offset_j = frag_chain_val->line_offset;

          free(frag_chain_val);
        }
        else
        {
          score = 0;
          num_dia_p = 0;
        }

        if (istep == 1)
        {
          pair_score[j][i] = score;
          pair_score[i][j] = score;

          // char *file_line = (char *)calloc(1, ((4 * sizeof(int)) + (2 * sizeof(double)) + (8 * sizeof(char))) + 1);
          // strcpy(file_line, "");
          // sprintf(file_line, ">%d,%d:%g\n>%d,%d:%g\n", i, j, score, j, i, score);
          // fprintf(pairscores_f, "%s", file_line);
          // free(file_line);
        }

        // for (k = 0; k < num_dia_p; k++)
        // {
        //   *current_dia = pair_dia[k];
        //   printf("k:%d num_dia_p:%d current_dia->s[0]:%d current_dia->s[1]:%d current_dia->weight:%f\n", k, num_dia_p, current_dia->s[0] + 1, current_dia->s[1] + 1, current_dia->weight); // DEBUG

        //   // printf("current_dia->next->s[0]:%d current_dia->next->s[1]:%d current_dia->next->weight:%f\n", current_dia->next->s[0] + 1, current_dia->next->s[1] + 1, current_dia->next->weight); // DEBUG

        //   current_dia->next = (struct multi_frag *)calloc(1, sizeof(struct multi_frag)); // calloc() returns null if current_dia->next is null. malloc gives garbage values :(
        //   // current_dia->next = (struct multi_frag *)malloc(1 * sizeof(struct multi_frag));
        //   // // int next_anc = current_dia->anc_num + 1;
        //   // int next_anc = current_dia->line_num + 1;
        //   // // int next_s1 = current_dia->s[0] + 1;
        //   // // int next_s2 = current_dia->s[1] + 1;
        //   // current_dia->next = get_anchor_by_num(mmapped_fasta ,mmapped_anc, &next_anc, current_dia->line_offset, 1);
        //   // // current_dia->next = get_anchor_seqs(mmapped_fasta, mmapped_anc, &next_s1, &next_s2, current_dia->fasta_offset, 1);
        //   end_dia = current_dia;
        //   // current_dia = current_dia->next;
        //   // current_dia = get_anchor_empty(); // doesn't work. it requires something other than empty anchor

        //   current_dia = current_dia->next;

        //   printf("(next) k:%d num_dia_p:%d current_dia->s[0]:%d current_dia->s[1]:%d current_dia->weight:%f\n", k, num_dia_p, current_dia->s[0] + 1, current_dia->s[1] + 1, current_dia->weight); // DEBUG
        //   current_dia->pred = end_dia;
        //   printf("(end dia) k:%d num_dia_p:%d current_dia->pred->s[0]:%d current_dia->pred->s[1]:%d current_dia->pred->weight:%f\n", k, num_dia_p, current_dia->pred->s[0] + 1, current_dia->pred->s[1] + 1, current_dia->pred->weight); // DEBUG
        // }

        // num_dia_bf[istep] = num_dia_bf[istep] + num_dia_p;

        // printf(" num_dia_bf[istep]:%d\n", num_dia_bf[istep]); // DEBUG

        for (hv = 0; hv < num_dia_p; hv++)
          weight_sum_af = weight_sum_af + (pair_dia[hv]).weight;

        if (num_dia_p)
          free(pair_dia);

        // #pragma omp critical
        //         {
        //           if (num_dia_p && pair_dia[k] != NULL)
        //             free(pair_dia[k]);
        //         }

      } /*    for(j = i+1 ; j<seqnum ; j++) */

      // free(seq1_val->data);
      free(seq1_val);
      // free(seqname_n1_val->data);
      free(seqname_n1_val);

    } /*   for(i = 0 ; i<seqnum ; i++) */

    // free(pair_dia);
    // exit(0); // DEBUG

    if (time_stamps)
    {
      end_pa = clock();

      time_diff_pa = (double)(end_pa - beg_pa) / CLOCKS_PER_SEC;
      if (time_stamps)
        printf(" time_diff_pa = %f \n", time_diff_pa);
      total_pa_time = total_pa_time + time_diff_pa;
    }

    if (break1)
    {
      printf("\n  break1\n");
      exit(1);
    }

    /*
        if( pa_only ) {
          printf("\n\n istep = %d, pa finished - exit \n\n", istep );
          exit(1);
        }
    */

    // int test_counter = 0;
    input_offset = 0;
    // size_t fasta_offset = 0;

    mmapped_file *mmapped_diags = mmap_file(&original_diag_fname, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
    total_diagonals = get_seqcount_mmapped_fasta(mmapped_diags);
    // exit(0); // DEBUG

    if (overlap_weights)
    {

      ow_add_mmap(total_diagonals, mmapped_diags, mmapped_fasta);
      printf("total diagonals:%d\n", total_diagonals); // DEBUG
      // exit(0);                                         // DEBUG

      if (bubblesort)
        // ow_bubble_sort(num_dia_bf[istep], this_it_dia);
        // ow_bubble_sort(total_diagonals, this_it_dia);
        exit(0); // DEBUG - placeholer
      else
      {
//         // frag_sort(mmapped_fasta, mmapped_anc, num_dia_bf[istep], this_it_dia, overlap_weights);
//         // frag_sort(mmapped_file *mmapped_fasta, mmapped_file *mmapped_anc, int number, struct multi_frag *dp, int olw);
// #pragma omp critical
//         {
//           frag_sort_mmap(num_dia_bf[istep], overlap_weights, mmapped_fasta, mmapped_ow, this_it_dia);
//         }
#pragma omp critical
        {
          frag_sort_mmap(total_diagonals, overlap_weights, mmapped_fasta, mmapped_anc, mmapped_diags);
        }
      }
      // munmap(mmapped_diags->mapped_file, mmapped_diags->sb.st_size);
      // free(mmapped_diags->file_name);
      // strcat(diag_fname, ".sorted");
    }
    else /* no overlap_weights */
    {
      beg_ts = clock();

      if (bubblesort)
        // // bubble_sort(num_dia_bf[istep], this_it_dia);
        // bubble_sort(num_dia_bf[istep], this_it_dia, mmapped_anc, mmapped_fasta, mmapped_ow);
        exit(0); // DEBUG - placeholder
      else
      {
        //         // frag_sort(mmapped_file *mmapped_fasta, mmapped_file *mmapped_anc, int number, struct multi_frag *dp, int olw);
        //         //  frag_sort(mmapped_fasta, mmapped_anc, num_dia_bf[istep], this_it_dia, overlap_weights);
        // #pragma omp critical
        //         {
        //           frag_sort_mmap(num_dia_bf[istep], overlap_weights, mmapped_fasta, mmapped_ow, this_it_dia);
        //         }
        // mmapped_diags = mmap_file(&diag_fname, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
        // strcpy(&diag_fname, mmapped_diags->file_name);
        // total_diagonals = get_seqcount_mmapped_fasta(mmapped_diags);
#pragma omp critical
        {
          frag_sort_mmap(total_diagonals, overlap_weights, mmapped_fasta, mmapped_anc, mmapped_diags);
        }
        // munmap(mmapped_diags->mapped_file, mmapped_diags->sb.st_size);
        // free(mmapped_diags->file_name);
        // strcat(diag_fname, ".sorted");
      }

      end_ts = clock();
      time_diff_srt = (double)(end_ts - beg_ts) / CLOCKS_PER_SEC;
      if (time_stamps)
        printf(" time_diff_srt = %f \n", time_diff_srt);
    }

    num_dia_af[istep] = num_dia_bf[istep];
    weight_sum_bf = weight_sum_af;

    pairalignsum = 0;
    pairalignlen = 0;

    // // #pragma omp critical
    // // {
    // filter(num_dia_af + istep, this_it_dia, mmapped_anc, mmapped_fasta, mmapped_ow);
    // char diags_name_sorted[NAME_LEN];
    // strcpy(diags_name_sorted, diag_fname);
    // strcat(diags_name_sorted, ".sorted");
    // if (istep > 1)
    // {          // DEBUG
    //   exit(0); // DEBUG //for 1st iteration debugging
    // } // DEBUG
    // munmap(mmapped_diags->file_name, mmapped_diags->sb.st_size);
    // free(mmapped_diags->file_name);
    // mmapped_diags = mmap_file(&original_diag_fname, mmapped_diags->file_mode, mmapped_diags->protocol, mmapped_diags->map_mode, 0);
    // printf("original_diag_fname:%s diag_fname:%s file_name:%s\n", original_diag_fname, sorted_diag_fname, mmapped_diags->file_name); // DEBUG
    // total_diagonals = get_seqcount_mmapped_fasta(mmapped_diags);
    mmapped_file *sorted_diags = mmap_file(&sorted_diag_fname, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
    int total_diagonals_sorted = get_seqcount_mmapped_fasta(sorted_diags);
    // printf("HIT PLACEHOLDER total_seqs:%d\n", get_seqcount_mmapped_fasta(mmapped_fasta)); // DEBUG
    // exit(0);
    filter(&total_diagonals_sorted, 0, mmapped_anc, mmapped_fasta, sorted_diags, mmapped_ow);
    // munmap(mmapped_diags->mapped_file, mmapped_diags->sb.st_size);
    // free(mmapped_diags->file_name);
    printf("HIT PLACEHOLDER filter()\n"); // DEBUG
    // exit(0);                              // DEBUG - placeholder
    // // }
    num_all_it_dia = num_all_it_dia + num_dia_af[istep];

    /*
        if( pa_only == 0 ) {
          printf("\n\n istep = %d, filter finished - exit \n\n", istep );
          exit(1);
        }
    */
    // printf("original_diag_fname:%s diag_fname:%s file_name:%s\n", original_diag_fname, sorted_diag_fname, mmapped_diags->file_name); // DEBUG
    weight_sum_af = 0;
    // munmap(mmapped_diags->file_name, mmapped_diags->sb.st_size);
    // mmapped_diags = mmap_file(&sorted_diag_fname, mmapped_diags->file_mode, mmapped_diags->protocol, mmapped_diags->map_mode, 0);
    // mmapped_file *all_diags = mmap_file(&original_diag_fname, O_RDONLY, PROT_READ, MAP_SHARED, 0);
#pragma omp critical
    print_log(fp_log, fp_fsm, mmapped_fasta, sorted_diags);
    // printf("HIT PLACEHOLDER print_log()\n"); // DEBUG
    // exit(0); // DEBUG - placeholder

    if (frag_file)
      print_fragments(sorted_diags, fp_frg);
    // printf("HIT PLACEHOLDER print_fragments()\n"); // DEBUG
    // exit(0);                                       // DEBUG - placeholder

    // printf("HIT PLACEHOLDER throw_out()\n"); // DEBUG
    // exit(0);                                 // DEBUG
    throw_out(&weight_sum_af, sorted_diags);

    // printf("HIT PLACEHOLDER sel_test()\n"); // DEBUG
    // exit(0);                                // DEBUG - placeholder
    // sel_test(mmapped_diags);

    threshold = threshold;

    if (break2)
    {
      printf("\n  break2\n");
      exit(1);
    }

    // munmap(all_diags->file_name, all_diags->sb.st_size);
    // free(all_diags->file_name);
    // free(all_diags);
    munmap(mmapped_diags->file_name, mmapped_diags->sb.st_size);
    free(mmapped_diags->file_name);
    free(mmapped_diags);

    munmap(sorted_diags->file_name, sorted_diags->sb.st_size);
    free(sorted_diags->file_name);
    free(sorted_diags);

    printf("max_itnum:%d\n", max_itnum); // DEBUG
    // if(istep >0){ //DEBUG
    // exit(0);                             // DEBUG //for 1st iteration debugging
    // } //DEBUG
  } /* while ( cond_it == 1 ) */
  // fclose(diag_f);
  // exit(0);
  // if (pair_dia != NULL)
  //   free(pair_dia);

  // printf("\n here 18 \n"); // DEBUG

  /***************************\
  *                           *
  *       ITERATION END       *
  *                           *
  \***************************/

  strcpy(dist_name, printname);
  strcat(dist_name, ".dst");

  if (ref_seq == 0)
  {
    // mmapped_file *all_diags = mmap_file(&original_diag_fname, O_RDONLY, PROT_READ, MAP_SHARED, 0);
    av_tree_print_mmap(mmapped_fasta, mmapped_ow);
  }
  // exit(0);
  if (standard_out)
    fp_ali = stdout;

  // printf("\n here 18.2 \n"); // DEBUG

  if (sf_mat)
  {
    if (wgt_type == 0 || wgt_type > 1)
    { // printf("\n here 18.2.-3 \n"); // DEBUG
      // subst_mat_mmap(input_name, num_all_it_dia, all_it_dia, mmapped_fasta);
      printf("HIT PLACEHOLDER subst_mat_mmap()\n"); // DEBUG
      exit(0);                                      // DEBUG - placeholder
    }
    else
    {
      perror("Substitution matrix output not supported for this wgt_type (-n)\n");
    }
  }

  printf("DEBUG num_all_it_dia:%d\n", num_all_it_dia); // DEBUG
  if (textual_alignment)
    ali_arrange_mmap(num_all_it_dia, all_it_dia, fp_ali, fp2, fp3, fp4, fp_csc, mmapped_fasta);
  printf("HIT PLACEHOLDER ali_arrange_mmap()\n"); // DEBUG
  // exit(0);                                        // DEBUG - placeholder

  if (long_output)
  {
    /*     fprintf(fp_log "\n\n thr = %f , lmax = %d , speed = %f  */
    fprintf(fp_log, "\n\n    total sum of weights: %f \n\n\n", tot_weight);
    fclose(fp_log);
  }

  if (argnum == 1)
  {
    printf("\n     Program terminated normally\n");
    printf("     Results are contained in file `%s' \n \n \n", itname);
  }

  av_dia_num = 2 * dia_num;
  av_dia_num = av_dia_num / (seqnum * (seqnum - 1));

  av_max_dia_num = 2 * max_dia_num;
  av_max_dia_num = av_max_dia_num / (seqnum * (seqnum - 1));

  tmpi1 = av_dia_num;
  tmpi2 = av_max_dia_num;

  if (pr_av_nd)
    printf(" %d ", tmpi1);

  if (pr_av_max_nd)
    printf(" %d ", tmpi2);

  if (pr_av_nd)
    fprintf(fp_ali, "    %d fragments considered for alignment \n", tmpi1);

  if (pr_av_max_nd)
    fprintf(fp_ali, "    %d fragments simultaneously stored \n\n", tmpi2);

  if (textual_alignment)
    fclose(fp_ali);

  if (time_stamps)
  {
    end_ali = clock();
    time_diff_ali = (double)(end_ali - beg_ali) / CLOCKS_PER_SEC;

    perc_pa_time = total_pa_time / time_diff_ali * 100;
    printf(" time_diff_ali = %f \n", time_diff_ali);
    printf(" total_pa_time = %f \n", total_pa_time);
    printf(" corresponds to %f percent \n\n", perc_pa_time);
  }
  munmap(mmapped_fasta->mapped_file, file_size);
  munmap(mmapped_anc->mapped_file, file_size_anc);
  munmap(mmapped_ow->mapped_file, mmapped_ow->sb.st_size);
  // munmap(mmapped_diags->mapped_file, mmapped_diags->sb.st_size);
  // free(sb);
  free(mmapped_fasta->file_name);
  free(mmapped_anc->file_name);
  free(mmapped_ow->file_name);
  free(mmapped_fasta);
  free(mmapped_anc);
  // free(mmapped_diags);
  free(mmapped_ow);
  // exit(0);

  free(num_dia_bf);
  free(num_dia_af);
  // free(all_it_dia);

  if (textual_alignment)
    if (fasta_file)
      fclose(fp2);

  if (msf_file)
    fclose(fp3);

  if (cw_file)
    fclose(fp4);

  unlink(sorted_diag_fname);
} /* main */
