
/********************\
*                    *
*    DIALIGN 2       *
*                    *
*    frag_chain.c    *
*                    *
\********************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "define.h"
#include "dialign.h"
#include "alig_graph_closure.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

extern short **mot_pos;
extern double mot_factor, mot_offset_factor, max_mot_offset;
extern int self_comparison, ***exclude_list;
extern short crick_strand, exclude_frg;
extern char pst_name[NAME_LEN];
extern int wgt_type, dna_speed;
extern double **wgt_prot, **wgt_dna, **wgt_trans;
extern int istep, lmax;
// extern char *seq[MAX_SEQNUM];
extern double av_sim_score_pep;
extern double av_sim_score_nuc;
// extern int *seqlen;

/* o.k. with read only global var. */
extern char seq_file[NAME_LEN];

extern int mot_len, motifs, iter_cond_prob, wgt_print, wgt_print_x;
extern struct multi_frag *pair_dia;
extern int afc_file, afc_filex, dia_pa_file;
extern int thr_sim_score;
extern CLOSURE *clos;
// extern int ***open_pos;
extern int sim_score[21][21];
extern int long_output;
extern int min_dia, max_dia, strict, seqnum;
extern int **amino, **amino_c;
extern char dia_pa_name[NAME_LEN];
extern char input_name[NAME_LEN];
extern FILE *fp_dia, *fp_dpa, *fp_mot;
// extern char *seq_name[MAX_SEQNUM];
extern int print_max_nd, pr_av_max_nd;
extern int dia_num, redundant, print_status, dcount, cont_it;
extern char input_line[NAME_LEN];
extern int max_dia_num;
extern double threshold;

extern double maxf2(double a, double b);
extern void rel_wgt_calc(int l1, int l2, double **wgt);
extern int mini2(int a, int b);
extern void weight_print(double **wgt);
extern int mini3(int a, int b, int c);
extern void wgt_prnt(int seqlen_l1, int seqlen_l2);
extern void wgt_prnt_x();
extern float mot_dist_factor(int offset, float parameter);

extern fasta_value *get_seq_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_len_value *get_seqlen_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_value *get_seqname_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern openpos_value get_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset);
extern openpos_value set_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, int *value, size_t input_offset);
// extern ow_value set_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, double *value, size_t input_offset);
extern mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset);

frag_chain_value *frag_chain_mmap(int n1, int n2, FILE *fp1, FILE *fp_m, int *number, mmapped_file *mapped_fasta, size_t input_offset, size_t seqlen_n1, char *seq_n1, char *seqname_n1)
{
  /* pairwise alignment */

  /* `i'  denotes positions in the 1. sequence ( seq[n1] ),
     `j'  denotes positions in the 2. sequence ( seq[n2] )  */

  //   int fd = open(seq_file, O_RDONLY);
  //   if (fd == -1) {
  //       perror("Error opening file");
  //       exit(1);
  //   }

  //   struct stat sb;
  //   if (fstat(fd, &sb) == -1) {
  //       perror("Error getting file size");
  //       exit(1);
  //   }

  //   // mmap the file into memory
  //   char *fasta = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  //   if (fasta == MAP_FAILED) {
  //       perror("Error mmapping file");
  //       exit(1);
  //   }

  //   close(fd);

  frag_chain_value *ret_val = (frag_chain_value *)malloc(1 * sizeof(frag_chain_value));
  // struct multi_frag_min *ret_val = (struct multi_frag_min *)calloc(1, sizeof(struct multi_frag_min));

  ret_val->line_offset = -1;
  ret_val->value = -1;

  // fasta_len_value seqlen_n1_val = get_seqlen_mmapped_fasta(mapped_fasta, &n1, sb, input_offset);
  // int seqlen_n1 = seqlen_n1_val.len;
  // fasta_value seq1_val = get_seq_mmapped_fasta(mapped_fasta, &n1, sb, input_offset);
  // char seq_n1[seqlen_n1];
  // strcpy(seq_n1, seq1_val.data);
  // fasta_value seqname_n1_val = get_seqname_mmapped_fasta(mapped_fasta, &n1, sb, input_offset);
  // int seqname_n1_len = strlen(seqname_n1_val.data);
  // char seqname_n1[seqname_n1_len];
  // strcpy(seqname_n1, seqname_n1_val.data);

  // printf("seq_n1: %s\n", seq_n1); // DEBUG

  // printf("seqlen_n1: %d\n", seqlen_n1); // DEBUG
  // exit(0);

  // input_offset = seq1_val.line_offset; // +strlen(seq_n1) + 1;

  // ret_val.line_offset = input_offset;

  fasta_value *seqname_n2_val = get_seqname_mmapped_fasta(mapped_fasta, &n2, input_offset);
  int seqname_n2_len = strlen(seqname_n2_val->data);
  char seqname_n2[seqname_n2_len + 1];
  // seqname_n2_val->data++;
  strcpy(seqname_n2, seqname_n2_val->data);

  input_offset = seqname_n2_val->line_offset + strlen(seqname_n2) + 1;
  fasta_len_value *seqlen_n2_val = get_seqlen_mmapped_fasta(mapped_fasta, &n2, input_offset);
  int seqlen_n2 = seqlen_n2_val->len;
  fasta_value *seq2_val = get_seq_mmapped_fasta(mapped_fasta, &n2, input_offset);
  char seq_n2[seqlen_n2 + 2];
  // seq2_val->data++;
  strcpy(seq_n2, seq2_val->data);
  // *seq_n2 = seq_n2[1];
  // input_offset = seq2_val->line_offset + strlen(seq_n1) + 1;

  ret_val->line_offset = seq2_val->line_offset + strlen(seq_n2) + 1;

  // printf("n1:%d n2:%d\n", n1, n2);        // DEBUG
  // printf("seq_n1: %s\n", seq_n1);         // DEBUG
  // printf("seq_n2: %s\n", seq_n2);         // DEBUG
  // printf("seqlen_n1: %d\n", seqlen_n1);   // DEBUG
  // printf("seqlen_n2: %d\n", seqlen_n2);   // DEBUG
  // printf("seqname_n1: %s\n", seqname_n1); // DEBUG
  // printf("seqname_n2: %s\n", seqname_n2); // DEBUG

  // exit(0);

  // FILE *open_pos_f;
  //  char *open_pos_fname; // = seq_file;
  int seq_file_length = sizeof(seq_file) / sizeof(char);
  char open_pos_fname[seq_file_length]; // = seq_file;
  // snprintf(open_pos_fname, sizeof(seq_file), "%s", seq_file);
  strcpy(open_pos_fname, seq_file);
  strcat(open_pos_fname, ".op");

  // printf("open_pos_fname: %s\n", open_pos_fname); // DEBUG

  // int op_fd = open(open_pos_fname, O_RDWR);
  // if (op_fd == -1)
  // {
  //   perror("Error opening file");
  //   exit(1);
  // }
  // struct stat op_sb;
  // if (fstat(op_fd, &op_sb) == -1)
  // {
  //   perror("Error getting file size");
  //   exit(1);
  // }
  // char *op_arr = mmap(NULL, op_sb.st_size, PROT_WRITE, MAP_SHARED, op_fd, 0);
  // if (op_arr == MAP_FAILED)
  // {
  //   perror("Error mmapping file");
  //   exit(1);
  // }

  // close(op_fd);

  mmapped_file *mmapped_op = mmap_file(open_pos_fname, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);

  // printf("&(mmapped_op->sb): %p\n", &(mmapped_op->sb)); // DEBUG

  // printf("mmapped_op->sb.st_size: %ld\n", mmapped_op->sb.st_size); // DEBUG

  // struct stat *op_sb = &(mmapped_op->sb);
  size_t file_size = mmapped_op->sb.st_size;

  int mot_match[MAX_DIA * 3];

  int mot_match_num, mot_offset;
  double mot_wgt_sum, this_mot_wgt;

  double thr; /* threshold for the weight of fragments starting
                at a given point (i,j). For any new pair (i,j),
                thr = 0. However, if a fragment with positive weight w
                is found starting at (i,j), thr is defined to be w
                and any further fragment starting at (i,j) is
                taken into consideration, only if its weight excedes
                thr. This is, because it is not meaningful to consider
                a fragment containing another fragment with
                higher weight. */
  int i, j, k, l, m, diff1, diff2, diff3, diff4, hv, hv2, numsubseq;
  int ende2;   /* denote the last position considered in the 2nd
                  sequence. Coincides with seqlen[n2], respectively,
                  exept if nucleotide diagonals are translated into
                  peptide diagonals. In this case,
                           ende2 = seqlen[n2]-2  */
  int start_a; /* diagonals begining at a position (i,j) are only
                 considered if the similarity-value at (i,j)
                 exceeds a certain threshold, respectively if
                 seq[n1][i] = seq[n2][j]. In this case the value
                 of `start_a' is 1, otherwise the value is 0  */
  int start_pep, start_pep_c, start_dna, start_dna1, trpl_start;
  int match; /* number of matches or sum of similarity values
                in a given diagonal */
  int trans, start_count, crick_wgt = 0;
  int match_p, match_p_c, match_d;
  int kmaxloc; /* maximum length of diagonals starting at a given
                  position (i,j) of the dot matrix.
                  kmaxloc = min{ max_dia, seqlen[n1]-i+1 , seqlen[n2]-j+1} */
  int lmax_real;
  int mnum = 0; /* number of current diagonal */
  int *ub_int = NULL;
  int *lb_int = NULL;
  int limit;      /* min { ub_int[i] ; ende2 }  */
  int bound_test; /* = 1 , if diagonal under consideration is consistent
                     with ub_int and lb_int.
                     = 0 , if not.   */
  int max_nd = 0, new_region = 0, current_nd = 0;
  short accepted = 0;
  char ch;
  double total_sum, wgt_k_match, wgt_k_match_c;

  double thr2, mot_wgt;

  struct pair_frag **diap = NULL;     /* diap[i] = pointer to last diagonal ending
                                       in the (i-1)-th column  */
  struct pair_frag **prec_vec = NULL; /* prec_vec[j] = pointer to diagonal with
                                       maximum sum of weights accumulated
                                       at a given position (i,j)  */
  struct pair_frag *current_dia = NULL, *hp = NULL, *cp = NULL, *cd = NULL;

  FILE *fp_st, *nd_fp;

  /*
   printf( "\n  in frag_chain: iter = %d wgt_type = %d  \n\n", istep , wgt_type );
   printf( "\n  in frag_chain: iter = %d wgt_dna 20 = %f \n\n", istep ,   wgt_dna[ 20 ][ 20 ] );
  */

  if (print_status)
    if (seqnum > 20)
    {
      fp_st = fopen(pst_name, "w");

      fprintf(fp_st, "\n\n\n    Status of the program run:\n");
      fprintf(fp_st, "    ==========================\n\n");
      fprintf(fp_st, "      %s \n\n", input_line);
      fprintf(fp_st, "      iteration step %d in multiple alignment\n\n", istep);
      fprintf(fp_st, "      aligning seq %d /", n1 + 1);
      fprintf(fp_st, " seq %d\n", n2 + 1);
      fprintf(fp_st, "      total number of");
      fprintf(fp_st, " sequences: %d\n\n", seqnum);
      fprintf(fp_st, "\n\n\n");

      fclose(fp_st);
    }

  //  if ( (ub_int = (int *) calloc( ( seqlen[n1] + 3 ) , sizeof(int) ) ) == NULL)
  if ((ub_int = (int *)calloc((seqlen_n1 + 3), sizeof(int))) == NULL)
  {
    printf("problems with memory allocation for ub_int!  \n \n");
    exit(1);
  }

  //  if ( (lb_int = (int *) calloc( (seqlen[n1]+3) , sizeof(int) ) ) == NULL)
  if ((lb_int = (int *)calloc((seqlen_n1 + 3), sizeof(int))) == NULL)
  {
    printf("problems with memory allocation for lb_int!  \n \n");
    exit(1);
  }

  //  if ( (prec_vec = (struct pair_frag **)
  //          calloc( (seqlen[n2]+3) , sizeof(struct pair_frag *) ) ) == NULL)
  if ((prec_vec = (struct pair_frag **)
           calloc((seqlen_n2 + 3), sizeof(struct pair_frag *))) == NULL)
  {
    printf("problems with memory allocation for prec_vec!  \n \n");
    exit(1);
  }

  //  if(
  //      (diap = (struct pair_frag **) calloc( (seqlen[n1] + 3) ,
  //               sizeof(struct pair_frag *) ))  ==  NULL
  //    )
  if (
      (diap = (struct pair_frag **)calloc((seqlen_n1 + 3),
                                          sizeof(struct pair_frag *))) == NULL)
  {
    printf("\n \n \n      ATTENTION: \n \n \n");
    printf("      problems with memory allocation\n");
    printf("      for diagonals! \n");
    exit(1);
  }

//  for( i = 1 ; i<= seqlen[n1] ; i++ )
#pragma omp parallel for num_threads(omp_get_max_threads())
  // for (i = 1; i <= seqlen_n1; i++)
  for (i = 1; i <= seqlen_n1; i++)
  {
    diap[i] = NULL;
    // diap[i] = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag ));
  }

  if ((diap[0] = (struct pair_frag *)
           calloc(1, sizeof(struct pair_frag))) == NULL)
  {
    printf("problems with memory allocation for diap!  \n \n");
    exit(1);
  }

  // diap[0]->b1 = diap[0]->b2 = diap[0]->ext = 0;
  // diap[0]->weight = 0.0;
  // diap[0]->trans = diap[0]->cs = 0;
  // diap[0]->prec = diap[0]->last = NULL;
  // diap[0]->sum = 0.0;

  // prec_vec[0] = NULL;
//  for( j = 1 ; j< seqlen[n2]+3 ; j++ )
#pragma omp parallel for num_threads(omp_get_max_threads())
  for (j = 1; j < seqlen_n2 + 3; j++)
  {
    // prec_vec[j] = NULL;
    // // prec_vec[j] = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag));
    // prec_vec[j]->last = prec_vec[j]->prec = NULL;
    // prec_vec[j]->sum = 0.0;
    // prec_vec[j]->b1 = prec_vec[j]->b2 = prec_vec[j]->ext = 0;
    // prec_vec[j]->weight = 0.0;
    // prec_vec[j]->trans = prec_vec[j]->cs = 0;
#pragma omp atomic read
    prec_vec[j] = diap[0];
  }

  // printf("\n here 17.6 \n"); //DEBUG

  if (dia_pa_file)
    fp_dpa = fopen(dia_pa_name, "a");

  //    ende2 = seqlen[n2];
  ende2 = seqlen_n2;

  /* Calculation of rel_weight  */

  if (iter_cond_prob == 0)
  {
    //    if( wgt_type == 0 )
    //      rel_wgt_calc( seqlen[n1] , seqlen[n2] , wgt_prot );
    //    if( wgt_type % 2 )
    //      rel_wgt_calc( seqlen[n1] , seqlen[n2] , wgt_dna);
    //    if( wgt_type > 1 )
    //      rel_wgt_calc( seqlen[n1] , seqlen[n2] , wgt_trans);
    if (wgt_type == 0)
      rel_wgt_calc(seqlen_n1, seqlen_n2, wgt_prot);
    if (wgt_type % 2)
      rel_wgt_calc(seqlen_n1, seqlen_n2, wgt_dna);
    if (wgt_type > 1)
      rel_wgt_calc(seqlen_n1, seqlen_n2, wgt_trans);

    // wgt_prnt(seqlen_n1, seqlen_n2); // DEBUG

    if (istep == 1)
      if (wgt_print || wgt_print_x)
      {
        wgt_prnt(seqlen_n1, seqlen_n2);
        if (wgt_print_x)
          exit(1);
      }

  } /* if( iter_cond_prob == 0 ) */

// printf("\n here 17.7 \n"); //DEBUG

//  for( hv = 1 ; hv <= seqlen[ n1 ] ; hv++ )
#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 1; hv <= seqlen_n1; hv++)
  // for (hv = 0; hv < seqlen_n1; hv++)
  {
    lb_int[hv] = predFrontier(clos, n1, hv, n2);
    ub_int[hv] = succFrontier(clos, n1, hv, n2);
    if (lb_int[hv] != ub_int[hv])
    {
      {
#pragma omp atomic update
        lb_int[hv]++;
#pragma omp atomic update
        ub_int[hv]--;
      }
    }
  }

  mnum = 0;

  if (iter_cond_prob || (istep == 1))
    new_region = 1;

  /* DP START */

  size_t op_offset1 = 0, op_offset2 = 0, op_offset3 = 0, op_offset4 = 0;

  // printf("seqlen_n1: %d seqlen_n2: %d\n", seqlen_n1, seqlen_n2); // DEBUG

  //  for( i = 1 ; i <= seqlen[n1] ; i++ )
  for (i = 1; i <= seqlen_n1; i++)
  // for (i = 0; i < seqlen_n1; i++)
  {
    if (op_offset1 >= file_size)
    {
      printf("Error: Could not find open position for %d, %d, %d\n", n1, n2, i);
      exit(1);
    }
    //  if( open_pos[n1][n2][i] )

    openpos_value open_pos_val = get_openpos_mmapped(mmapped_op, &n1, &n2, &i, op_offset1);
    // printf("n1:%d n2:%d i: %d open_pos_val.value: %c \n", n1, n2, i, open_pos_val.value + '0'); // DEBUG
    op_offset1 = open_pos_val.line_offset;
    // printf(" 1. open_pos[%d][%d][i:%d]:%d\n", n1, n2, i, open_pos_val.value); // DEBUG
    // printf(" seq[n1][%d]:%c seq[n2][%d]:%c seq[n1][i] == seq[n2][j]:%d \n", i, seq_n1[i], j, seq_n2[j], seq_n1[i] == seq_n2[j]); // DEBUG
    if (open_pos_val.value > 0)
    {
      if (new_region)
      {

        diff2 = (succFrontier(clos, n1, i, n2) - predFrontier(clos, n1, i, n2) - 1);

        if (diff2 < 0)
          diff2 = 0;

        diff1 = (succFrontier(clos, n2, lb_int[i], n1) - predFrontier(clos, n2, lb_int[i], n1) - 1);
        if (diff1 < 0)
          diff1 = 0;

        /*
          printf(" new region, i = %d diff = %d , %d \n", i, diff1 , diff2  );
        */

        // printf("\n here 17.7.1 \n"); //DEBUG

        if (iter_cond_prob)
          if ((diff1 > 0) && (diff2 > 0))
          {
            if (wgt_type == 0)
              rel_wgt_calc(diff1, diff2, wgt_prot);
            if (wgt_type % 2)
              rel_wgt_calc(diff1, diff2, wgt_dna);
            if (wgt_type > 1)
              rel_wgt_calc(diff1, diff2, wgt_trans);
          }
      }

      limit = mini2(ub_int[i], ende2);
      for (j = lb_int[i]; j <= limit; j++)
      {
        // printf(" seq[n1][%d]:%c seq[n2][%d]:%c seq[n1][i] == seq[n2][j]:%d \n", i, seq_n1[i], j, seq_n2[j], seq_n1[i] == seq_n2[j]); // DEBUG
        if (wgt_type != 1)
          start_pep = (sim_score[amino[n1][i]][amino[n2][j]] >= thr_sim_score);

        if (crick_strand)
          start_pep_c = (sim_score[amino_c[n1][i]][amino_c[n2][j]] >= thr_sim_score);

        if (wgt_type % 2)
        {

          //    if( strict )
          //      start_dna = ( (seq[n1][i] == seq[n2][j]) &&
          //                    ( seq[n1][i] == 'A' ||
          //                      seq[n1][i] == 'C' ||
          //                      seq[n1][i] == 'T' ||
          //                      seq[n1][i] == 'G' ||
          //                      seq[n1][i] == 'U'  ) );

          //    else
          //      start_dna = (seq[n1][i] == seq[n2][j]);

          if (strict)
          {
            start_dna = ((seq_n1[i] == seq_n2[j]) &&
                         (seq_n1[i] == 'A' ||
                          seq_n1[i] == 'C' ||
                          seq_n1[i] == 'T' ||
                          seq_n1[i] == 'G' ||
                          seq_n1[i] == 'U'));
          }
          else
          {
            start_dna = (seq_n1[i] == seq_n2[j]);
          }

          // printf("\n here 17.7.2 \n"); //DEBUG

          if (dna_speed)
            //    if( ( i < seqlen[n1] ) && ( j < limit ) )
            if ((i < seqlen_n1) && (j < limit))
            {
              //  if( strict )
              //    start_dna1 = ( (seq[n1][ i + 1 ] == seq[ n2 ][ j + 1 ]) &&
              //                  ( seq[n1][ i + 1 ] == 'A' ||
              //                    seq[n1][ i + 1 ] == 'C' ||
              //                    seq[n1][ i + 1 ] == 'T' ||
              //                    seq[n1][ i + 1 ] == 'G' ||
              //                    seq[n1][ i + 1 ] == 'U'  ) );

              //  else
              //    start_dna1 = ( seq[ n1 ][ i + 1 ] == seq[ n2 ][ j + 1 ] );
              if (strict)
                start_dna1 = ((seq_n1[i + 1] == seq_n2[j + 1]) &&
                              (seq_n1[i + 1] == 'A' ||
                               seq_n1[i + 1] == 'C' ||
                               seq_n1[i + 1] == 'T' ||
                               seq_n1[i + 1] == 'G' ||
                               seq_n1[i + 1] == 'U'));

              else
                start_dna1 = (seq_n1[i + 1] == seq_n2[j + 1]);
              start_dna = start_dna * start_dna1;
            }
        }

        // printf("\n here 17.7.3 \n"); //DEBUG

        if (wgt_type != 1)
          start_a = start_pep;
        else
          start_a = start_dna;

        if (wgt_type == 3)
          start_a = start_pep + start_dna;

        if (crick_strand)
          start_a = start_a + start_pep_c;

        if (self_comparison)
          if (i == j)
            start_a = 0;

        if (exclude_frg)
          if (j == exclude_list[n1][n2][i])
            start_a = 0;

        if (start_a)
        {

          match = 0;
          match_d = 0;
          match_p = 0;
          match_p_c = 0;
          thr = 0;
          /*
                           start_count++ ;
          */
          bound_test = 1;

          if (wgt_type > 1)
            lmax_real = lmax * 3;
          else
            lmax_real = lmax;

          //  kmaxloc =
          //    mini3( lmax_real , seqlen[n1]-i+1 , seqlen[n2]-j+1 );
          kmaxloc =
              mini3(lmax_real, seqlen_n1 - i + 1, seqlen_n2 - j + 1);

          if (motifs)
          {
            for (k = 1; k <= kmaxloc; k++)
              if ((mot_pos[n1][i + k - 1] == 1) &&
                  (mot_pos[n2][j + k - 1] == 1))
              {
                mot_match[k] = 1;
                /*     printf(" match in  %d  %d  %d \n", i, j, k ); */
              }
              else
                mot_match[k] = 0;
          }

          /*******************\
          *                   *
          *  fragments start  *
          *                   *
          \*******************/

          k = 1;
          mot_match_num = 0;
          mot_wgt_sum = 0;

          while ((k <= kmaxloc) && start_a)
          {
            if (motifs)
            {
              if (((i - j) * (i - j)) < (max_mot_offset * max_mot_offset))
                if (k >= mot_len)
                  if (mot_match[k - mot_len + 1])
                  {
                    mot_offset = (i - j);
                    this_mot_wgt = mot_dist_factor((i - j), mot_offset_factor);

                    // printf("\n here 17.7.4 \n"); //DEBUG

                    /*
                                           printf("  i - j = %d , tmw = %f \n", i - j , this_mot_wgt );
                    */
                    mot_wgt_sum = mot_wgt_sum + this_mot_wgt;
                    mot_match_num++;
                  }
            }
            int i_deferred_idx1 = i + k - 1;
            int i_deferred_idx2 = i + k;
            int i_deferred_idx3 = i + k + 1;
            //  if( open_pos[n1][n2][ i + k - 1 ] )

            openpos_value open_pos_val1 = get_openpos_mmapped(mmapped_op, &n1, &n2, &i_deferred_idx1, op_offset2);
            // printf("n1:%d n2:%d i + k - 1: %d get open_pos_val1: %c\n", n1, n2, i + k - 1, open_pos_val1.value + '0'); // DEBUG
            op_offset2 = open_pos_val1.line_offset;
            if (open_pos_val1.value > 0)
            {
              bound_test = bound_test *
                           (j + k - 1 >= lb_int[i + k - 1]);
              bound_test = bound_test *
                           (j + k - 1 <= ub_int[i + k - 1]);

              trpl_start = 0;

              if (wgt_type < 2)
                trans = 0;
              else
                trans = 1;

              if (start_pep ||
                  (crick_strand && start_pep_c))
                if ((wgt_type > 1) && ((k % 3) == 1))
                {
                  trpl_start = 1;

                  trpl_start = trpl_start *
                               (j + k >= lb_int[i + k]);
                  trpl_start = trpl_start *
                               (j + k <= ub_int[i + k]);
                  //  trpl_start = trpl_start * open_pos[ n1 ][ n2 ][ i + k ] ;
                  // printf("n1:%d n2:%d i + k: %d\n", n1, n2, i + k); // DEBUG
                  openpos_value open_pos_def1 = get_openpos_mmapped(mmapped_op, &n1, &n2, &i_deferred_idx2, op_offset3);
                  // printf("get open_pos_def1:%c i + k - 1:%d\n", open_pos_def1.value + '0', i_deferred_idx2); // DEBUG
                  op_offset3 = open_pos_def1.line_offset;
                  trpl_start = trpl_start * (open_pos_def1.value);

                  trpl_start = trpl_start *
                               (j + k + 1 >= lb_int[i + k + 1]);
                  trpl_start = trpl_start *
                               (j + k + 1 <= ub_int[i + k + 1]);
                  //  trpl_start = trpl_start * open_pos[ n1 ][ n2 ][ i + k + 1 ] ;
                  // printf("n1:%d n2:%d i + k + 1: %d\n", n1, n2, i + k + 1); // DEBUG
                  openpos_value open_pos_def2 = get_openpos_mmapped(mmapped_op, &n1, &n2, &i_deferred_idx3, op_offset4);
                  op_offset4 = open_pos_def2.line_offset;
                  // printf("get open_pos_def2:%c i + k + 1:%d\n", open_pos_def2.value + '0', i_deferred_idx3); // DEBUG
                  trpl_start = trpl_start * (open_pos_def2.value);
                }

              // printf("\n here 17.8 \n"); //DEBUG

              if (
                  bound_test &&
                  ((wgt_type != 2) || trpl_start))
              {
                if (start_pep)
                  if (
                      (wgt_type == 0) ||
                      ((wgt_type > 1) && trpl_start))
                    match_p = match_p + sim_score[amino[n1][i + k - 1]]
                                                 [amino[n2][j + k - 1]];

                if (crick_strand)
                  if (start_pep_c)
                    if (
                        (wgt_type == 0) ||
                        ((wgt_type > 1) && trpl_start))
                      match_p_c = match_p_c + sim_score[amino_c[n1][i + k - 1]]
                                                       [amino_c[n2][j + k - 1]];

                if (start_dna)
                  if (wgt_type % 2)
                  {
                    //    if( !strict ||
                    //       (seq[n1][i+k-1] == 'A' ||
                    //        seq[n1][i+k-1] == 'C' ||
                    //        seq[n1][i+k-1] == 'T' ||
                    //        seq[n1][i+k-1] == 'G' ||
                    //        seq[n1][i+k-1] == 'U'  ))
                    //      match_d = match_d +
                    //      ( seq[n1][i+k-1] == seq[n2][j+k-1] );
                    if (!strict ||
                        (seq_n1[i + k - 1] == 'A' ||
                         seq_n1[i + k - 1] == 'C' ||
                         seq_n1[i + k - 1] == 'T' ||
                         seq_n1[i + k - 1] == 'G' ||
                         seq_n1[i + k - 1] == 'U'))
                      match_d = match_d +
                                (seq_n1[i + k - 1] == seq_n2[j + k - 1]);
                  }

                wgt_k_match = 0;

                if (wgt_type == 0)
                  wgt_k_match = wgt_prot[k][match_p];
                if (wgt_type == 1)
                  wgt_k_match = wgt_dna[k][match_d];

                if (wgt_type > 1)
                {
                  if (start_pep)
                    wgt_k_match = wgt_trans[(k + 2) / 3][match_p];

                  if (crick_strand)
                    if (start_pep_c)
                    {
                      if (wgt_trans[(k + 2) / 3][match_p_c] > wgt_k_match)
                      {
                        wgt_k_match = wgt_trans[(k + 2) / 3][match_p_c];
                        crick_wgt = 1;
                      }
                      else
                        crick_wgt = 0;
                    }
                }

                if (start_dna)
                  if (wgt_type == 3)
                    if (k <= lmax)
                      if (wgt_dna[k][match_d] > wgt_k_match)
                      {
                        wgt_k_match = wgt_dna[k][match_d];
                        trans = 0;
                      }

                if (wgt_type == 0)
                  if (match_p <= (k * av_sim_score_pep))
                    start_pep = 0;

                if (wgt_type == 1)
                  if (match_d <= (k * av_sim_score_nuc))
                    start_dna = 0;

                if (start_pep)
                  if (wgt_type > 1)
                    if ((match_p * 3) <= (k * av_sim_score_pep))
                      start_pep = 0;

                if (start_pep_c)
                  if (wgt_type > 1)
                    if ((match_p_c * 3) <= (k * av_sim_score_pep))
                      start_pep_c = 0;

                if (wgt_type != 1)
                  start_a = start_pep;
                else
                  start_a = start_dna;

                if (wgt_type == 3)
                  start_a = start_pep + start_dna;

                if (crick_strand)
                  start_a = start_a + start_pep_c;

                if (exclude_frg)
                  if (exclude_list[n1][n2][i + k] == j + k)
                    start_a = 0;

                if (motifs)
                  if (mot_wgt_sum > 0)
                  {
                    fprintf(fp_mot, "  %4d %4d  ", n1 + 1, n2 + 1);
                    fprintf(fp_mot, "  %4d %4d %3d  ", i, j, k);
                    fprintf(fp_mot, "  %5.2f ", wgt_k_match);
                    mot_wgt = mot_wgt_sum * mot_factor;
                    wgt_k_match = wgt_k_match + mot_wgt;
                    fprintf(fp_mot, "  %2d ", mot_match_num);
                    fprintf(fp_mot, "     %5.2f \n", wgt_k_match);
                  }
                //  printf("\n here 17.9 \n"); //DEBUG

                /*
                   if( wgt_k_match > 0 )     printf(" k = %d min_dia = %d wgt_k_match = %f thr = %f \n", k, min_dia , wgt_k_match , thr );
                */
                if (k >= min_dia)
                  if (wgt_k_match > thr)
                  {
                    // if ((current_dia = (struct pair_frag *)
                    //          calloc(1, sizeof(struct pair_frag))) == NULL)
                    if ((current_dia = (struct pair_frag *)
                             calloc(1, sizeof(struct pair_frag))) == NULL)
                    {
                      printf("\n \n \n      ATTENTION: \n \n \n");
                      printf("      too many diagonals in\n");
                      printf("      pairwise alignment of");
                      printf(" sequences\n");
                      printf("      %s  and  ", seqname_n1);
                      printf("%s\n \n \n \n", seqname_n2);

                      fprintf(fp1, "\n \n      ATTENTION:\n \n");
                      fprintf(fp1, "      too many diagonals\n");
                      fprintf(fp1, "      in pairwise alignment");
                      fprintf(fp1, " of sequences\n");
                      fprintf(fp1, "      %s  and  ", seqname_n1);
                      fprintf(fp1, "%s\n \n \n \n", seqname_n2);

                      exit(1);
                    }

                    current_dia->b1 = i;
                    current_dia->b2 = j;
                    current_dia->ext = k + 2 * trans;
                    current_dia->weight = wgt_k_match;
                    current_dia->trans = trans;
                    current_dia->cs = crick_wgt;

                    // printf(" current_dia->weight:%f \n", current_dia->weight); // DEBUG

                    // if (prec_vec[j] == NULL)
                    // {
                    //   // printf(" prec_vec[%d] = NULL \n", j); // DEBUG
                    //   prec_vec[j] = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag *));
                    //   prec_vec[j] = diap[0];
                    // }
                    // printf(" i:%d j:%d current_dia->weight:%f (prec_vec[j])->sum:%f\n", i, j, current_dia->weight, (prec_vec[j])->sum); // DEBUG
                    current_dia->sum = current_dia->weight + (prec_vec[j])->sum;
                    // current_dia->prec = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag *));
                    current_dia->prec = prec_vec[j];
                    // current_dia->last = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag *));
                    current_dia->last = diap[i + (current_dia->ext)];
                    diap[i + (current_dia->ext)] = current_dia;

                    // printf(" diap[%d]:%p current_dia->sum:%f prec_vec[%d]:%p current_dia->prec:%p current_dia->last:%p\n", i + (current_dia->ext), diap[i + (current_dia->ext)], current_dia->sum, j, prec_vec[j], current_dia->prec, current_dia->last); // DEBUG

                    mnum++;

                    if (print_max_nd)
                    {
                      current_nd++;
                      if (current_nd > max_nd)
                        max_nd = current_nd;
                    }

                    dia_num++;
                    if (afc_file)
                    {
                      fprintf(fp_dia, "FRG %d ", dia_num);
                      fprintf(fp_dia, "name: %s ", seqname_n1);
                      fprintf(fp_dia, " %s ", seqname_n2);
                      if (seqnum > 2)
                      {
                        fprintf(fp_dia, "  seq: %d", n1 + 1);
                        fprintf(fp_dia, " %d", n2 + 1);
                      }
                      fprintf(fp_dia, "  beg: %d %d", i, j);
                      fprintf(fp_dia, "  len: %d", current_dia->ext);
                      fprintf(fp_dia, "  wgt: %6.3f", current_dia->weight);

                      /*
                                                         if( BETA )
                                                         if( iter_cond_prob )
                                                           {
                                                             fprintf(fp_dia,"   d1 = %d d2 = %d ", diff1, diff2 );
                                                           }
                      */

                      fprintf(fp_dia, "  it = %d ", istep);
                      if ((wgt_type == 3) || crick_strand)
                        if (current_dia->trans)
                          fprintf(fp_dia, " P-frg");
                        else
                          fprintf(fp_dia, " N-frg");
                      fprintf(fp_dia, "\n");
                      if (afc_filex)
                      {
                        fprintf(fp_dia, "SEG1 ");
                        for (hv = 0; hv < current_dia->ext; hv++)
                        {
                          //    ch = seq[n1][ i + hv ] ;
                          ch = seq_n1[i + hv];
                          fprintf(fp_dia, "%c", ch);
                        }
                        fprintf(fp_dia, "\n");
                        fprintf(fp_dia, "SEG2 ");
                        for (hv = 0; hv < current_dia->ext; hv++)
                        {
                          //    ch = seq[n2][ j + hv ] ;
                          ch = seq_n2[j + hv];
                          fprintf(fp_dia, "%c", ch);
                        }
                        fprintf(fp_dia, "\n\n");
                      }
                    }

                    if (!redundant)
                    {
                      thr2 = maxf2(thr, (current_dia->weight));
                      thr = thr2;
                    }

                  } /*   if( wgt[k][match] > thr  )  */
              } /*   if ( bound_test )              */
            } /*   if( open_pos ...  )  */
            k++;
          } /*   while( ( k <= kmaxloc ) && start_a ) */
        } /*   if( start_a )                    */
      } /*   for(j=lb_int[i];j<=limit;j++)   */
      new_region = 0;
    } /*   if( open_pos )                     */
    else if (iter_cond_prob)

      new_region = 1;

    if (print_status)
      //  if( ( ( seqlen[n1] + seqlen[n2] ) > 1000 ) )
      if (((seqlen_n1 + seqlen_n2) > 1000))
        if (!(i % 100))
        {
          fp_st = fopen(pst_name, "w");

          fprintf(fp_st, "\n\n\n    Status of the program run:\n");
          fprintf(fp_st, "    ==========================\n\n");
          fprintf(fp_st, "      %s \n\n", input_line);
          if (seqnum > 2)
          {
            fprintf(fp_st, "      iteration step %d in", istep);
            fprintf(fp_st, " multiple alignment\n\n");
          }
          if (seqnum > 2)
          {
            fprintf(fp_st, "      aligning seq %d /", n1 + 1);
            fprintf(fp_st, " seq %d\n", n2 + 1);
            fprintf(fp_st, "      total number of");
            fprintf(fp_st, " sequences: %d\n\n", seqnum);
          }
          fprintf(fp_st, "      current position in");
          fprintf(fp_st, " sequence %d:  %8d\n", n1 + 1, i);
          fprintf(fp_st, "      length of seq %d:", n1 + 1);

          // fprintf(fp_st,"                 %8d\n\n", seqlen[n1]);
          fprintf(fp_st, "                 %8d\n\n", seqlen_n1);

          /*
                      if( iter_cond_prob || ( istep == 1 ) )
                        {
                          if( open_pos[n1][n2][i] )
                            {
                              fprintf(fp_st,"      diff1 = %d \n", diff1 );
                              fprintf(fp_st,"      diff2 = %d \n", diff2 );
                            }
                          else
                            fprintf(fp_st,"      position already aligned");
                        }
          */

          fprintf(fp_st, "\n\n\n");

          fclose(fp_st);
        }

    // cp = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag *));
    cp = diap[i + 1];
    hp = NULL;
    // hp = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag *));
    accepted = 0;

    while (cp != NULL)
    {
      j = cp->b2 + cp->ext;
      // printf(" j:%d (prec_vec[j])->sum:%f cp->sum:%f\n", j, (prec_vec[j])->sum, cp->sum); // DEBUG
      // printf("cp->b2:%d cp->ext:%d prec_vec[%d]->sum:%d\n", cp->b2, cp->ext, j, prec_vec[j]->sum); // DEBUG
      // printf(" prec_vec[%d]->sum:%d\n", j, prec_vec[j]->sum); // DEBUG
      double prev_vec_j_sum = 0;
#pragma omp atomic read
      prev_vec_j_sum = (prec_vec[j])->sum;
      // printf(" prec_vec[%d]->sum:%d cp:%p\n", j, prev_vec_j_sum, cp); // DEBUG
      if (prev_vec_j_sum < cp->sum)
      {
        prec_vec[j] = cp;
        accepted = 1;

        hp = cp;
        cp = cp->last;
      }
      else
      {
        cp = cp->last;

        if (accepted)
        {
          free(hp->last);

          hp->last = cp;
        }
        else
        {
          free(diap[i + 1]);
          diap[i + 1] = cp;
        }

        current_nd--;
      }
    }

    //  for( hv=2 ; hv < ( seqlen[n2] + 3 ) ; hv++ )

#pragma omp parallel for num_threads(omp_get_max_threads())
    for (hv = 2; hv < (seqlen_n2 + 3); hv++)
    {
      double prev_vec_hv_sum = (prec_vec[hv])->sum;
      double prev_vec_hvm1_sum = 0;
#pragma omp atomic read
      prev_vec_hvm1_sum = (prec_vec[hv - 1])->sum;
      if (prev_vec_hv_sum < prev_vec_hvm1_sum)
      {
#pragma omp atomic read
        prec_vec[hv] = prec_vec[hv - 1];
      }
    }

  } /*   for(i= ... )                   */
  /*
  printf (" start_count = %d \n ", start_count );
   */
  if (pr_av_max_nd)
    if (istep == 1)
    {
      if ((nd_fp = fopen("nd_file", "a")) == NULL)
      {
        printf("\n\n  nd_fp could not be opened \n\n");
        exit(1);
      }

      fprintf(nd_fp, " %2d/%2d %8d  %8d \n", n1 + 1, n2 + 1, mnum, max_nd);
      fclose(nd_fp);
    }

  numsubseq = 0; /* counts diagonals in alignment */

  //  current_dia = prec_vec[ seqlen[n2] + 1 ];
  current_dia = prec_vec[seqlen_n2 + 1];
  // cd = (struct pair_frag *)malloc(1 * sizeof(struct pair_frag *));
  cd = current_dia;

  total_sum = cd->sum;

  printf(" total_sum = %f \n", total_sum); // DEBUG
  while (cd->prec != NULL)
  {
    printf(" numsubseq = %d \n", numsubseq);           // DEBUG
    printf(" cd->weight = %f \n", cd->weight);         // DEBUG
    printf(" cd->sum = %f \n", cd->sum);               // DEBUG
    printf(" before - &(cd->prec) = %p \n", cd->prec); // DEBUG
    printf(" before - &(cd) = %p \n", cd);             // DEBUG
    // #pragma omp atomic update
    numsubseq++;
    cd = cd->prec;
    printf(" after - &(cd) = %p \n", cd);             // DEBUG
    printf(" after - &(cd->prec) = %p \n", cd->prec); // DEBUG
  }

  if (numsubseq)
  {
    hv = numsubseq - 1;

    // if (
    //     (
    //         pair_dia = (struct multi_frag *)
    //             calloc((numsubseq + 1), sizeof(struct multi_frag))) == NULL)
    if (
        (
            pair_dia = (struct multi_frag *)
                calloc((numsubseq + 1), sizeof(struct multi_frag))) == NULL)
    {
      printf("problems with memory allocation for `pair_dia'!  \n \n");
      exit(1);
    }

    // Creating/appending to diags file
    FILE *diag_f;
    char diag_fname[seq_file_length];
    strcpy(&diag_fname, &input_name);
    strcat(diag_fname, ".diags");
    diag_f = fopen(diag_fname, "a+");

    // // Creating/appending to diagonal pathtrace file
    // FILE *pathtrace_f;
    // char pathtrace_fname[seq_file_length];
    // strcpy(&pathtrace_fname, &input_name);
    // strcat(pathtrace_fname, ".trace");
    // pathtrace_f = fopen(pathtrace_fname, "a+");

    // int prev_s0, prev_s1, prev_hv, prev_idx;

    while (hv >= 0)
    {
      if (dia_pa_file)
      {
        fprintf(fp_dpa, " %3d) ", ++dcount);
        if (seqnum > 2)
          fprintf(fp_dpa, "seq: %3d %3d ", n1 + 1, n2 + 1);
        fprintf(fp_dpa, " beg: %6d %6d ", current_dia->b1, current_dia->b2);
        fprintf(fp_dpa, " len: %2d ", current_dia->ext);
        fprintf(fp_dpa, " weight: %5.2f ", current_dia->weight);
        fprintf(fp_dpa, " it: %d ", istep);
        if ((wgt_type == 3) || crick_strand)
          if (current_dia->trans)
            fprintf(fp_dpa, " P-frg");
          else
            fprintf(fp_dpa, " N-frg");
        if (current_dia->trans)
          if (crick_strand)
          {
            if (current_dia->cs)
              fprintf(fp_dpa, " crick ");
            else
              fprintf(fp_dpa, " watson ");
          }
        fprintf(fp_dpa, "\n");
      }

      (pair_dia[hv]).anc_num = n1 + 1; // viz
      (pair_dia[hv]).b[0] = current_dia->b1;
      (pair_dia[hv]).b[1] = current_dia->b2;
      (pair_dia[hv]).s[0] = n1;
      (pair_dia[hv]).s[1] = n2;
      (pair_dia[hv]).sel = 1;
      (pair_dia[hv]).ext = current_dia->ext;
      (pair_dia[hv]).weight = current_dia->weight;
      // if ((pair_dia[hv]).ow == 0)
      // {
      (pair_dia[hv]).ow = current_dia->weight;
      // }
      printf("  (pair_dia[hv]).weight:%f \n", (pair_dia[hv]).weight); // DEBUG
      printf("  (pair_dia[hv]).ow:%f \n", (pair_dia[hv]).ow);         // DEBUG
      (pair_dia[hv]).trans = current_dia->trans;
      if (crick_strand)
      {
        (pair_dia[hv]).cs = current_dia->cs;
      }
      else
      {
        (pair_dia[hv]).cs = (short)0;
      }
      (pair_dia[hv]).it = istep;

      // double default_addedwghts_val = 0;
      char *file_line = (char *)calloc(1, ((3 * sizeof(short)) + (9 * sizeof(int)) + (3 * sizeof(double)) + (17 * sizeof(char))) + 1);
      strcpy(file_line, "");
      sprintf(file_line, ">%d,%d,%d,%d:%d,%d,%d,%d,%u,%u,%d,%g,%g,%g,%u\n", (pair_dia[hv]).s[0], (pair_dia[hv]).s[1], hv, (pair_dia[hv]).it, (pair_dia[hv]).anc_num, numsubseq, (pair_dia[hv]).b[0], (pair_dia[hv]).b[1], (pair_dia[hv]).sel, (pair_dia[hv]).trans, (pair_dia[hv]).ext, (pair_dia[hv]).weight, (pair_dia[hv]).ow, total_sum, (pair_dia[hv]).cs);
      // printf("file_line:%s\n", file_line); // DEBUG
      // exit(0);                             // DEBUG
#pragma omp critical
      fprintf(diag_f, "%s", file_line);
      // printf("%s", file_line); // DEBUG
      free(file_line);
      // printf(">%d,%d,%d,%d:%d,%d,%d,%d,%u,%u,%d,%g,%g,%g,%u\n", (pair_dia[hv]).s[0], (pair_dia[hv]).s[1], hv, istep, (pair_dia[hv]).anc_num, numsubseq, (pair_dia[hv]).b[0], (pair_dia[hv]).b[1], (pair_dia[hv]).sel, (pair_dia[hv]).trans, (pair_dia[hv]).ext, (pair_dia[hv]).weight, (pair_dia[hv]).ow, total_sum, (pair_dia[hv]).cs); // DEBUG

      hv--;
      current_dia = current_dia->prec;
    }
    fclose(diag_f);
    // exit(0); // DEBUG

    /*     if( dia_pa_file )
           fprintf(fp_dpa, " \n" );
    */

    /*    modified in LGI-VITRY
           if( iter_cond_prob )
    */

    cont_it = 1;

  } /* if(numsubseq) */

  *number = numsubseq;

  if (long_output)
  {
    printf("Seq. %3d -%3d: ", n1 + 1, n2 + 1);
    printf("T = %2.2f,", threshold);
    printf(" %3d D. in alignment,", *number);
    printf("%6d D. in matrix", mnum);
    printf("\n");
  }

  // int iter = 0;
  // //  for( hv=0 ; hv < seqlen[n1]+3 ; hv++ )
  for (hv = 0; hv < seqlen_n1 + 3; hv++)
  {
    current_dia = diap[hv];

    // while (!current_dia)
    while (current_dia != NULL)
    {
      hp = current_dia;
      // struct pair_frag *tmp_hp = current_dia;

      current_dia = current_dia->last;
      // current_dia = current_dia->prec;
      // printf(" freeing hp: iter:%d sum:%d weight:%d\n", iter, hp->sum, hp->weight); //DEBUG
      free(hp);
      // iter++; // DEBUG
      // free(tmp_hp);
    }
    // iter = 0; // DEBUG
  }

  if (istep == 1)
  {
    max_dia_num = max_dia_num + max_nd;
  }

  // for (size_t i = 0; i < seqlen_n2 + 3; i++)
  // {
  //   if (prec_vec[i])
  //   {
  //     printf(" free prec_vec[%d]:%p\n", i, prec_vec[i]); // DEBUG
  //     free(prec_vec[i]);
  //   }
  // }
  // free(prec_vec);

  // for (size_t i = 0; i < seqlen_n1 + 3; i++)
  // {
  //   if (diap[i] != NULL)
  //   {
  //     printf(" free diap[%d]:%p\n", i, diap[i]); // DEBUG
  //     free(diap[i]);
  //   }
  // }

  // free(seq2_val->data);
  // free(seqname_n2_val->data);
  free(seq2_val);
  free(seqname_n2_val);

  free(diap);

  free(ub_int);
  free(lb_int);
  free(prec_vec);
  if (dia_pa_file)
  {
    fclose(fp_dpa);
  }

  // free(seq_n1);
  // free(seq_n2);
  // free(seqname_n1);
  // free(seqname_n2);

  munmap(mmapped_op->mapped_file, file_size);
  // free(op_sb);
  free(mmapped_op->file_name);
  free(mmapped_op);

  // return (total_sum);
  ret_val->value = total_sum;

  // printf(" ret_val->line_offset:%d \n", ret_val->line_offset); // DEBUG
  printf("n1:%d n2:%d\n", n1, n2);        // DEBUG
  printf("seq_n1: %s\n", seq_n1);         // DEBUG
  printf("seq_n2: %s\n", seq_n2);         // DEBUG
  printf("seqlen_n1: %d\n", seqlen_n1);   // DEBUG
  printf("seqlen_n2: %d\n", seqlen_n2);   // DEBUG
  printf("seqname_n1: %s\n", seqname_n1); // DEBUG
  printf("seqname_n2: %s\n", seqname_n2); // DEBUG

  return ret_val;

} /* frag_chain */
