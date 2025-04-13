
/*************************\
*                         *
*     DIALIGN 2           *
*                         *
*     functions.c         *
*                         *
\*************************/

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

extern int iter_cond_prob, col_score;
extern short **cont_it_p;

extern char seq_file[NAME_LEN];
extern char pst_name[NAME_LEN];
extern int istep;
// extern char *seq[MAX_SEQNUM]; /* sequences */
// extern int *seqlen;           /* lengths of sequences */
extern double tot_weight;
extern int plot_num;
extern char clust_sim[NAME_LEN];
extern char *upg_str;

extern int total_anc_num;

extern char input_name[NAME_LEN];
extern int anchors;
extern int frg_mult_file;
extern int frg_mult_file_v;
extern short crick_strand;
/* extern char dia_pa_name[NAME_LEN]; */
extern int pr_av_max_nd, wgt_type;
extern int mask;
extern char prn[NAME_LEN];
extern char input_line[NAME_LEN];
extern int print_status;
extern char *seq_name[MAX_SEQNUM];
extern char printname[NAME_LEN];
extern char amino_acid[22];
extern int par_count;
extern int *num_dia_bf;
extern int *num_dia_af;

extern double pairalignsum;
extern int pairalignlen;

extern struct multi_frag *this_it_dia;
extern struct multi_frag *all_it_dia;

extern CLOSURE *clos;
// extern double **glob_sim;

// extern int ***open_pos;

extern double **wgt_prot, **wgt_dna, **wgt_trans;

extern int sim_score[21][21];

extern int min_dia, max_dia;
extern int long_output;
extern int seqnum;
extern short dots;
extern double threshold;
extern int num_all_it_dia;

extern int num_dia_p, overlap_weights;

extern int **amino;
extern int **amino_c;

extern int **shift;
extern double **tp400_prot;
extern double **tp400_dna;
extern double **tp400_trans;

// extern openpos_value get_openpos_mmapped(char* mapped_file, struct stat *sb, int *x, int* y, int* z) ;
extern ow_value get_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset);
extern ow_value set_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, double *value, size_t input_offset);
extern openpos_value set_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, int *value, size_t input_offset);
extern fasta_value *get_seq_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_len_value *get_seqlen_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_value *get_seqname_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern diag_value set_diags_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, double *value, size_t input_offset);
// extern diag_value get_diags_seq_mmapped(mmapped_file *mapped_file, int *x, int *y, int *val_index, size_t input_offset);
extern diag_value get_diagval_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, size_t input_offset);
extern diag_line get_diagvals_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, size_t input_offset);
extern diag_value get_diagval_line_mmapped(mmapped_file *mapped_file, int *line_num, int *val_index, size_t input_offset);
extern diag_line get_diagvals_line_mmapped(mmapped_file *mapped_file, int *line_num, size_t input_offset);
extern size_t get_seqcount_mmapped_fasta(mmapped_file *mapped_fasta);
extern line_value get_offset_of_line(mmapped_file *mapped_file, int *line_num, size_t input_offset);
extern line_value get_line_of_offset(mmapped_file *mapped_file, size_t input_offset);
extern mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset);
extern mmapped_file *resize_mmap_file(mmapped_file *mapped_file, size_t new_size);
extern mmapped_file *shift_mmap_contents(mmapped_file *mapped_fasta, mmapped_file *mapped_file, size_t start_offset, size_t shift_amount);
// extern struct multi_frag *get_anchor_by_num(char *mapped_anc, struct stat *sb_anc, int *anc_num, size_t offset_anc, int get_next);
extern struct multi_frag *get_anchor_by_num(mmapped_file *mapped_fasta, mmapped_file *mapped_anc, int *anc_num, size_t offset_anc, int get_next);
// extern mmapped_file *swap_multifrag_anchor(mmapped_file *mapped_file, int line_num1, int line_num2);
// extern struct multi_frag *get_anchor_by_num(char *mapped, struct stat *sb, int *anc_num, size_t input_offset, int get_next);
extern int get_anchor_count(char *file_name);
// extern size_t file_length_from_offset(mmapped_file *mapped_file, size_t input_offset);
extern size_t get_size_from_offset(mmapped_file *mapped_file, size_t input_offset);

int num_test(char *cp)
{
  int result = 1;
  int i;
  char *strng;

  strng = cp;

  for (i = 0; i < strlen(strng); i++)
    if (!isdigit(strng[i]) && (strng[i] != '.'))
    {
      result = 0;
      /*     printf("\n %c is no digit !!!\n", strng[i]);   */
    }

  return result;
}

double minf2(double a, double b)
{
  if (a < b)
    return a;
  else
    return b;
}

double maxf2(double a, double b)
{
  if (a > b)
    return a;
  else
    return b;
}

int mini2(int a, int b)
{
  if (a < b)
    return a;
  else
    return b;
}

int maxi2(int a, int b)
{
  if (a > b)
    return a;
  else
    return b;
}

int mini3(int a, int b, int c)
{
  return mini2(a, mini2(b, c));
}

void minf(double *a, double b)
{
  if (*a > b)
    *a = b;
}

void mini(int *a, int b)
{
  if (*a > b)
    *a = b;
}

void maxf(double *a, double b)
{
  if (*a < b)
    *a = b;
}

void maxi(int *a, int b)
{
  if (*a < b)
    *a = b;
}

char invert(char c1)
{
  char c2 = c1;

  if (c1 == 'T')
    c2 = 'A';
  if (c1 == 'A')
    c2 = 'T';
  if (c1 == 'C')
    c2 = 'G';
  if (c1 == 'G')
    c2 = 'C';

  return (c2);
}

int translate(char c1, char c2, char c3, int seqno, int pos)
{
  /* translation of triplets into amino acids */

  int amac; /* resulting amino acid */

  if (c1 == 'T')
  {
    if (c2 == 'T')
    {
      if (c3 == 'T')
        amac = 18;
      if (c3 == 'C')
        amac = 18;
      if (c3 == 'A')
        amac = 16;
      if (c3 == 'G')
        amac = 16;
    }
    if (c2 == 'C')
      amac = 2;
    if (c2 == 'A')
    {
      if (c3 == 'T')
        amac = 19;
      if (c3 == 'C')
        amac = 19;
      if (c3 == 'A')
        amac = 0; /* stop codon */
      if (c3 == 'G')
        amac = 0;
    }
    if (c2 == 'G')
    {
      if (c3 == 'T')
        amac = 1;
      if (c3 == 'C')
        amac = 1;
      if (c3 == 'A')
        amac = 20;

      if (c3 == 'G')
        amac = 20;
    }
  }

  if (c1 == 'C')
  {
    if (c2 == 'T')
      amac = 16;
    if (c2 == 'C')
      amac = 4;
    if (c2 == 'A')
    {
      if (c3 == 'T')
        amac = 11;
      if (c3 == 'C')
        amac = 11;
      if (c3 == 'A')
        amac = 10;
      if (c3 == 'G')
        amac = 10;
    }
    if (c2 == 'G')
      amac = 12;
  }

  if (c1 == 'A')
  {
    if (c2 == 'T')
    {
      if (c3 == 'T')
        amac = 15;
      if (c3 == 'C')
        amac = 15;
      if (c3 == 'A')
        amac = 15;
      if (c3 == 'G')
        amac = 14;
    }
    if (c2 == 'C')
      amac = 3;
    if (c2 == 'A')
    {
      if (c3 == 'T')
        amac = 7;
      if (c3 == 'C')
        amac = 7;
      if (c3 == 'A')
        amac = 13;
      if (c3 == 'G')
        amac = 13;
    }
    if (c2 == 'G')
    {
      if (c3 == 'T')
        amac = 2;
      if (c3 == 'C')
        amac = 2;
      if (c3 == 'A')
        amac = 12;
      if (c3 == 'G')
        amac = 12;
    }
  }

  if (c1 == 'G')
  {
    if (c2 == 'T')
      amac = 17;
    if (c2 == 'C')
      amac = 5;
    if (c2 == 'A')
    {
      if (c3 == 'T')
        amac = 8;
      if (c3 == 'C')
        amac = 8;
      if (c3 == 'A')
        amac = 9;
      if (c3 == 'G')
        amac = 9;
    }
    if (c2 == 'G')
      amac = 6;
  }

  if (
      (c1 != 'A' && c1 != 'T' && c1 != 'G' && c1 != 'C') ||
      (c2 != 'A' && c2 != 'T' && c2 != 'G' && c2 != 'C') ||
      (c3 != 'A' && c3 != 'T' && c3 != 'G' && c3 != 'C'))
    return (0);
  else
    return (amac);

} /*  translate */

int int_test(float f)
{
  int i = f;

  if (i == f)
    return (1);
  else
    return (0);
}

/*==========================================================
 *         OLD SORT FUNCTION (BUBBLE-SORT)
 *=========================================================*/

void change(struct multi_frag *a, struct multi_frag *b)
{
  struct multi_frag c, *an, *bn;

  c = *a;
  an = a->next;
  bn = b->next;

  *a = *b;
  *b = c;

  a->next = an;
  b->next = bn;
}

// void change_mmap(char *mapped_file, struct stat *sb, int *anc_num1, int *anc_num2, size_t line_offset1, size_t line_offset2)
void change_mmap(mmapped_file *mapped_file, struct multi_frag *anc_num1, struct multi_frag *anc_num2, size_t line_offset1, size_t line_offset2)
{
  // size_t line_offset1 = 0, line_offset2 = 0;
  // line_offset1 = get_offset_of_line(mapped_file, sb, anc_num1, 0);
  // line_offset2 = get_offset_of_line(mapped_file, sb, anc_num2, 0);

  char *line_map1 = &(mapped_file->mapped_file[line_offset1]);
  char *line_map2 = &(mapped_file->mapped_file[line_offset2]);

  int len1 = strcspn(line_map1, "\n");
  char *line1 = strndup(line_map1, len1);
  int line_len1 = strlen(line1);

  int len2 = strcspn(line_map2, "\n");
  char *line2 = strndup(line_map2, len2);
  int line_len2 = strlen(line2);

  char line1_buffer[line_len1], line2_buffer[line_len2];

  strcpy(line1_buffer, line1);
  strcpy(line2_buffer, line2);

  memmove(&(mapped_file->mapped_file[line_offset1]), line2_buffer, line_len2);
  memmove(&(mapped_file->mapped_file[line_offset2]), line1_buffer, line_len1);

  msync(&(mapped_file->mapped_file[line_offset1]), sizeof(char) * line_len2, MS_SYNC);
  msync(&(mapped_file->mapped_file[line_offset2]), sizeof(char) * line_len1, MS_SYNC);

  free(line1);
  free(line2);
}

void pair_change(struct seq_pair *a, struct seq_pair *b)
{
  struct seq_pair c;

  c = *a;
  *a = *b;
  *b = c;
}

void ow_bubble_sort(int number, struct multi_frag *dp)
{
  /* sorting diagonals in multiple alignment according to their
     overlap weights */

  struct multi_frag *hp;
  int hv1, hv2;

  FILE *fp_st;

  for (hv1 = 1; hv1 < number; hv1++)
  {
    hp = dp;

    if (print_status)
      if ((hv1 % 100) == 0)
      {
        fp_st = fopen(pst_name, "w");

        fprintf(fp_st, "\n\n\n    Status of the program run:\n");
        fprintf(fp_st, "    ==========================\n\n");
        fprintf(fp_st, "      %s \n\n", input_line);
        fprintf(fp_st, "      iteration step %d in multiple alignment\n", istep);
        fprintf(fp_st, "      overlap weight sorting of diagonals\n");
        fprintf(fp_st, "      current diagonal = %d\n\n", hv1);
        fprintf(fp_st, "      total number of");
        fprintf(fp_st, " diagonals: %d\n\n\n\n", number);
        fclose(fp_st);
      }

    for (hv2 = hv1; hv2 < number; hv2++)
    {
      if (hp->ow < (hp->next)->ow)
        change(hp, hp->next);
      hp = hp->next;
    }
  }
} /*  ow_bubble_sort */

void bubble_sort(int number, struct multi_frag *dp, mmapped_file *mapped_anc, mmapped_file *mapped_fasta)
{
  struct multi_frag *hp;
  int hv1, hv2;
  FILE *fp_st;

  for (hv1 = 1; hv1 < number; hv1++)
  {
    hp = dp;

    if (print_status)
      if ((hv1 % 100) == 0)
      {
        fp_st = fopen(pst_name, "w");

        fprintf(fp_st, "\n\n\n    Status of the program run:\n");
        fprintf(fp_st, "    ==========================\n\n");
        fprintf(fp_st, "      %s \n\n", input_line);
        fprintf(fp_st, "      iteration step %d\n", istep);
        fprintf(fp_st, "      ind. weight sorting of diagonals\n");
        fprintf(fp_st, "      current diagonal = %d\n\n", hv1);
        fprintf(fp_st, "      total number of");
        fprintf(fp_st, " diagonals: %d\n\n\n\n", number);
        fclose(fp_st);
      }

    for (hv2 = hv1; hv2 < number; hv2++)
    {
      // if (hp->weight < (hp->next)->weight)
      //   change(hp, hp->next);
      // hp = hp->next;
      // int next_anc = hv1 + 1;
      int next_anc = hp->next->anc_num;
      if (hp->weight < (hp->next)->weight)
      {
        // change(hp, hp->next);
        // int next_anc = hv1 + 1;
        change_mmap(mapped_anc, hp, hp->next, hp->line_offset, hp->next->line_offset);
      }
      // hp = hp->next;
      hp = get_anchor_by_num(mapped_fasta, mapped_anc, &next_anc, hp->next->line_offset, 1);
    }
  }
} /*bubble sort*/

void bubble_sort_mmap(int number, mmapped_file *mapped_anc, mmapped_file *mapped_fasta)
{
  /* sorting diagonals in multiple alignment according to their
       individual weights */

  printf(" here 15.10. bubblesort\n"); // DEBUG
  struct multi_frag *hp;
  int hv1, hv2;
  FILE *fp_st;
  unsigned int anc_idx = 0;

  // int fd;
  // struct stat sb;
  // char *mapped;
  size_t offset = 0;
  int sn = -1;

  // char anc_file_name[NAME_LEN];
  // strcpy(anc_file_name, input_name);
  // strcat(anc_file_name, ".anc");

  // Open the file
  // fd = open(anc_file_name, O_RDONLY);
  // if (fd == -1)
  // {
  //   perror("Error opening file");
  //   return -1;
  // }

  // // Get the size of the file
  // if (fstat(fd, &sb) == -1)
  // {
  //   perror("Error getting file size");
  //   close(fd);
  //   return -1;
  // }

  // // Map the file into memory
  // mapped = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  // if (mapped == MAP_FAILED)
  // {
  //   perror("Error mapping file");
  //   close(fd);
  //   return -1;
  // }

  // // Close the file descriptor, as it is no longer needed after mmap
  // close(fd);

  // mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
  // struct stat *anc_sb = &(mmapped_anc->sb);
  // size_t file_size = anc_sb->st_size;

  // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  // struct stat *sb_fasta = &(mmapped_fasta->sb);
  // size_t file_size_fasta = sb_fasta->st_size;

  size_t input_offset_hv1 = 0;

  for (hv1 = 1; hv1 < number; hv1++)
  {
    hp = get_anchor_by_num(mapped_fasta, mapped_anc, &hv1, input_offset_hv1, 1);
    input_offset_hv1 = hp->line_offset;

    if (print_status)
      if ((hv1 % 100) == 0)
      {
        fp_st = fopen(pst_name, "w");

        fprintf(fp_st, "\n\n\n    Status of the program run:\n");
        fprintf(fp_st, "    ==========================\n\n");
        fprintf(fp_st, "      %s \n\n", input_line);
        fprintf(fp_st, "      iteration step %d\n", istep);
        fprintf(fp_st, "      ind. weight sorting of diagonals\n");
        fprintf(fp_st, "      current diagonal = %d\n\n", hv1);
        fprintf(fp_st, "      total number of");
        fprintf(fp_st, " diagonals: %d\n\n\n\n", number);
        fclose(fp_st);
      }

    for (hv2 = hv1; hv2 < number; hv2++)
    {
      // int next_anc = hv1 + 1;
      int next_anc = hp->next->anc_num;
      if (hp->weight < (hp->next)->weight)
      {
        // change(hp, hp->next);
        // int next_anc = hv1 + 1;
        change_mmap(mapped_anc, hp, hp->next, hp->line_offset, hp->next->line_offset);
      }
      // hp = hp->next;
      hp = get_anchor_by_num(mapped_fasta, mapped_anc, &next_anc, hp->next->line_offset, 1);
      // input_offset = hp->line_offset;
    }
  }

  // Unmap the file after reading
  // munmap(mmapped_anc->mapped_file, file_size);
  // // free(anc_sb);
  // free(mmapped_anc);

} /*  bubble_sort */

/*==========================================================
 *         NEW SORT FUNCTION (QUICK-SORT)
 *=========================================================*/

/***********************************************************
 *                   change()
 ***********************************************************/
void change_struct_el_frg(struct multi_frag **a, int l, int r)
{
  struct multi_frag *dummy;
  dummy = a[l];
  a[l] = a[r];
  a[r] = dummy;
}

void change_struct_el_int(int **a, int l, int r)
{
  int *dummy;
  dummy = a[l];
  a[l] = a[r];
  a[r] = dummy;
}

// int iter = 0;

void change_struct_el_mmap(mmapped_file *mapped_anc, struct multi_frag *l_element, struct multi_frag *r_element)
{

  // iter++;

  // size_t shift_start_l = 0;
  // size_t shift_start_r = 0;

  size_t original_size = mapped_anc->sb.st_size;
  size_t l_offset = l_element->line_offset;
  size_t r_offset = r_element->line_offset;

  // Get pointers to the two lines
  char *line1 = &(mapped_anc->mapped_file)[l_offset];
  char *line2 = &(mapped_anc->mapped_file)[r_offset];

  // Determine the length of each line
  size_t len1 = strcspn(line1, "\n") + 1; // +1 to include newline
  size_t len2 = strcspn(line2, "\n") + 1; // +1 to include newline

  // Allocate temporary buffers to hold the lines
  char *temp1 = strndup(line1, len1);
  char *temp2 = strndup(line2, len2);

  printf(" Line 1: \"%.*s\"\n Line 2: \"%.*s\"\n", (int)len1, line1, (int)len2, line2);

  // strcat(temp1, "\n");
  // strcat(temp2, "\n");

  // strcat(temp1, "\n\0");
  // strcat(temp2, "\n\0");

  // printf(" (b) line1[0]:%c line2[0]:%c\n", line1[0], line2[0]); // DEBUG
  printf(" (b) temp1:%s temp2:%s\n", temp1, temp2); // DEBUG
  // printf(" (b) l_offset:%d r_offset:%d\n", l_offset, r_offset); // DEBUG

  // Calculate the difference in length between the lines
  int length_difference = 0;

  if (line1 > line2)
  {
    length_difference = len1 - len2;
  }
  else
  {
    length_difference = len2 - len1;
  }

  printf(" length_difference:%d len1:%d len2:%d\n", length_difference, len1, len2); // DEBUG

  // if (iter == 2)
  // {
  //   exit(0); // DEBUG
  // }

  // iter++; // DEBUG

  // printf(" iter:%d\n", iter); // DEBUG

  // If line lengths differ, adjust the file size and contents as needed
  if (length_difference != 0)
  {
    // // Resize the file to the new size
    size_t new_size = mapped_anc->sb.st_size + abs(length_difference);
    mapped_anc = resize_mmap_file(mapped_anc, new_size);

    // Update line1 and line2 pointers after remapping
    line1 = &(mapped_anc->mapped_file[l_offset]);
    line2 = &(mapped_anc->mapped_file[r_offset]);

    // printf(" new_size:%d original_size:%d new_size - original_size:%d\n", new_size, original_size, new_size - original_size); //DEBUG

    if (length_difference > 0)
    { // length_difference > 0
      if (r_offset > l_offset)
      {                                                            // r_offset > l_offset && length_difference > 0
        printf("r_offset > l_offset && length_difference > 0 \n"); // DEBUG
        size_t flength_l_offset = get_size_from_offset(mapped_anc, l_offset);
        memmove(&(mapped_anc->mapped_file[r_offset + len2 - length_difference]), &(mapped_anc->mapped_file[r_offset + len2]), (flength_l_offset - len1 - length_difference) * sizeof(char));
        memmove(line2, temp1, len1);
        msync(&(mapped_anc->mapped_file[r_offset]), len1, MS_SYNC);

        size_t flength_r_offset = get_size_from_offset(mapped_anc, r_offset);
        memmove(&(mapped_anc->mapped_file[l_offset + len1 + length_difference]), &(mapped_anc->mapped_file[l_offset + len1]), (flength_r_offset - len2 - length_difference) * sizeof(char));
        memmove(line1, temp2, len2);
        msync(&(mapped_anc->mapped_file[l_offset]), len2, MS_SYNC);
      } // r_offset > l_offset && length_difference > 0
      else
      {                                                            // r_offset < l_offset && length_difference > 0
        printf(" r_offset < l_offset && length_difference > 0\n"); // DEBUG
        size_t flength_r_offset = get_size_from_offset(mapped_anc, r_offset);
        memmove(&(mapped_anc->mapped_file[l_offset + len1 - length_difference]), &(mapped_anc->mapped_file[l_offset + len1]), (flength_r_offset - len2 - length_difference) * sizeof(char));
        memmove(line1, temp2, len2);
        msync(&(mapped_anc->mapped_file[l_offset]), len2, MS_SYNC);

        size_t flength_l_offset = get_size_from_offset(mapped_anc, l_offset);
        memmove(&(mapped_anc->mapped_file[r_offset + len2 + length_difference]), &(mapped_anc->mapped_file[r_offset + len2]), (flength_l_offset - len1 - length_difference) * sizeof(char));
        memmove(line2, temp1, len1);
        msync(&(mapped_anc->mapped_file[r_offset]), len1, MS_SYNC);
      } // r_offset < l_offset && length_difference > 0

      line_value l_line_val = get_line_of_offset(mapped_anc, l_offset + len1 + length_difference);
      line_value r_line_val = get_line_of_offset(mapped_anc, r_offset + len2 - length_difference);

      l_element->line_offset = l_line_val.line_offset;
      r_element->line_offset = r_line_val.line_offset;
    }
    else
    {
      // length_difference < 0
      if (r_offset > l_offset)
      {                                                           // r_offset > l_offset && length_difference < 0
        printf("r_offset > l_offset && length_difference < 0\n"); // DEBUG
        size_t flength_l_offset = get_size_from_offset(mapped_anc, l_offset);
        memmove(&(mapped_anc->mapped_file[l_offset + len1 - abs(length_difference)]), &(mapped_anc->mapped_file[l_offset + len1]), (flength_l_offset - len2) * sizeof(char));
        memmove(line1, temp2, len2);
        msync(&(mapped_anc->mapped_file[l_offset]), len2, MS_SYNC);

        size_t flength_r_offset = get_size_from_offset(mapped_anc, r_offset);
        memmove(&(mapped_anc->mapped_file[r_offset + len2]), &(mapped_anc->mapped_file[r_offset + len2 - abs(length_difference)]), (flength_r_offset - len1) * sizeof(char));
        memmove(&(mapped_anc->mapped_file[r_offset - abs(length_difference)]), temp1, len1);
        msync(&(mapped_anc->mapped_file[r_offset]), len1, MS_SYNC);
      } // r_offset > l_offset && length_difference < 0
      else
      {                                                           // r_offset < l_offset && length_difference < 0
        printf("r_offset < l_offset && length_difference < 0\n"); // DEBUG
        size_t flength_r_offset = get_size_from_offset(mapped_anc, r_offset);
        memmove(&(mapped_anc->mapped_file[r_offset + len2 - abs(length_difference)]), &(mapped_anc->mapped_file[r_offset + len2]), (flength_r_offset - len1) * sizeof(char));
        memmove(line2, temp1, len1);
        msync(&(mapped_anc->mapped_file[r_offset]), len1, MS_SYNC);

        size_t flength_l_offset = get_size_from_offset(mapped_anc, l_offset);
        memmove(&(mapped_anc->mapped_file[l_offset + len2 - abs(length_difference)]), &(mapped_anc->mapped_file[l_offset + len1 - abs(length_difference)]), (flength_l_offset - len2) * sizeof(char));
        memmove(&(mapped_anc->mapped_file[l_offset - abs(length_difference)]), temp2, len2);
        msync(&(mapped_anc->mapped_file[l_offset - abs(length_difference)]), len2, MS_SYNC);
      } // r_offset < l_offset &&  length_difference < 0

      line_value l_line_val = get_line_of_offset(mapped_anc, l_offset + len1 + abs(length_difference));
      line_value r_line_val = get_line_of_offset(mapped_anc, r_offset + len2 - abs(length_difference));

      l_element->line_offset = l_line_val.line_offset;
      r_element->line_offset = r_line_val.line_offset - abs(length_difference);
    } // length_difference < 0
  } // length_difference != 0
  else
  {                                      // length_difference == 0
    printf(" length_difference == 0\n"); // DEBUG
    memmove(line1, temp2, len2);
    memmove(line2, temp1, len1);

    msync(&(mapped_anc->mapped_file[l_offset]), len2, MS_SYNC);
    msync(&(mapped_anc->mapped_file[r_offset]), len1, MS_SYNC);
  } // length_difference == 0

  printf(" (swapped anc) anc_num_l:%d anc_num_r:%d\n", l_element->anc_num, r_element->anc_num); // DEBUG

  if (original_size != mapped_anc->sb.st_size)
  {
    // printf(" original size != new size : %d != %d\n", original_size, mapped_anc->sb.st_size);
    mapped_anc = resize_mmap_file(mapped_anc, original_size); // Shrink file again after shifting to resolve NULL bit intrusion
  }

  // Free temporary buffers
  free(temp1);
  free(temp2);
}

// void change_struct_el_mmap(char *mapped_anc, struct stat *sb_anc, size_t l_offset, size_t r_offset)
// {
//   // struct multi_frag *dummy;
//   // dummy = a[l];
//   // a[l] = a[r];
//   // a[r] = dummy;

//   // int *dummy;
//   // dummy = a[l];
//   // a[l] = a[r];
//   // a[r] = dummy;

//   // Get pointers to the two lines in the mmapped file
//   char *line1 = &mapped_anc[l_offset];
//   char *line2 = &mapped_anc[r_offset];

//   // Determine the length of each line for swapping
//   size_t len1 = strcspn(line1, "\n") + 1; //+1 for newline char
//   size_t len2 = strcspn(line2, "\n") + 1; //+1 for newline char

//   // Allocate memory to temporarily store the lines
//   char *temp1 = strndup(line1, len1);
//   char *temp2 = strndup(line2, len2);

//   // memset(temp1 + len1, '\n', 1);
//   // memset(temp2 + len2, '\n', 1); // Pad temp2 with spaces

//   printf(" (b) l_line:%s r_line:%s len1:%d len2:%d\n", temp1, temp2, len1, len2);

//   if (len1 != len2)
//     if (len1 < len2)
//     {
//       memset(temp1 + len1, ' ', len2 - len1); // Pad temp1 with spaces
//       len1 = len2;                            // Update len1 to match len2 for consistent swapping
//     }
//     else if (len2 < len1)
//     {
//       memset(temp2 + len2, ' ', len1 - len2); // Pad temp2 with spaces
//       len2 = len1;                            // Update len2 to match len1 for consistent swapping
//     }

//   printf(" (a) l_line:%s r_line:%s len1:%d len2:%d\n", temp1, temp2, len1, len2);                           // DEBUG
//   printf(" mapped_anc[l_offset]:%c mapped_anc[r_offset]:%c\n", mapped_anc[l_offset], mapped_anc[r_offset]); // DEBUG
//   // exit(0);
//   // Perform the swap
//   memmove(&mapped_anc[l_offset], temp2, len2);
//   memmove(&mapped_anc[r_offset], temp1, len1);

//   // Synchronize changes to the file
//   msync(&mapped_anc[l_offset - len1], len1 + len2, MS_SYNC);
//   msync(&mapped_anc[r_offset - len2], len2 + len1, MS_SYNC);

//   printf(" (a) mapped_anc[l_offset]:%c mapped_anc[r_offset]:%c\n", mapped_anc[l_offset], mapped_anc[r_offset]); // DEBUG

//   // Free temporary buffers
//   free(temp1);
//   free(temp2);
// }
/***********************************************************
 *                    change_first()
 ***********************************************************/
void change_first(struct multi_frag *a, struct multi_frag *b)
{
  struct multi_frag c, *an, *bn;

  if (a == b)
  { /* Change the first list-element with the second one (old first-el.). */
    c = *a;
    an = a->next;
    bn = a->next->next;

    *a = *(a->next);
    a->next = bn;

    *an = c;
    an->next = a;
  }
  else /* Change the new first list-el. with the old first-el. */
  {
    c = *a;             /* Make a copy of the new first-listelement a. */
    an = a->next;       /* Make a copy of the pointer at the second-el. */
    bn = b->next->next; /* Make a copy of the pointer old first-el. shows at. */

    *a = *(b->next); /* Whrite the value of the old first-el. on the place of the new first-el.*/
    a->next = bn;    /* Bend his "next" pointer at the next el. of the old first-el. */

    *(b->next) = c;     /* Whrite the value of the new fist-el. on the place of the old first-el. */
    b->next->next = an; /* Bend his "next" pointer at the next el. of the new first-el. */

    b->next = a;
  }
}

/***********************************************************
 *                    quicksort_ow()
 ***********************************************************/
void quicksort_ow(int **array, int left, int right, mmapped_file *mapped_diags, size_t l_offset, size_t r_offset, int *val_index, int *total_diagonals)
{

  if (left >= right)
    return;

  int zero_idx = 0;
  int l = left;
  int r = right;
  int pivot_idx = array[(left + right) / 2];

  // for (size_t i = 0; i < total_anc_num; i++)
  // {
  //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
  // }

  printf(" left:%d right:%d pivot_idx:%d\n", left, right, pivot_idx); // DEBUG

  // Get the pivot weight for comparison
  diag_value pivot = get_diagval_line_mmapped(mapped_diags, &pivot_idx, val_index, l_offset);

  while (l <= r)
  {
    // Move 'l' right until we find an element >= pivot weight
    int array_l = array[l]; // because anc_num is one indexed
    // printf(" array_l:%d\n", array_l); // DEBUG
    diag_value element_l = get_diagval_line_mmapped(mapped_diags, &array_l, val_index, l_offset);
    // l_offset = element_l.line_offset;
    while (element_l.value > pivot.value && l < *total_diagonals)
    {
      l++;
      if (l > right)
        break;
      array_l = array[l];
      // printf(" array_l:%d\n", array_l); // DEBUG
      element_l = get_diagval_line_mmapped(mapped_diags, &array_l, val_index, element_l.line_offset);
      // l_offset = element_l.line_offset;
    }
    // free(element_l);
    l_offset = element_l.line_offset;

    // Move 'r' left until we find an element <= pivot weight
    int array_r = array[r]; // because anc_num is 1 indexed
    // printf(" array_r:%d array[r]:%d\n", array_r, array[r]); // DEBUG
    diag_value element_r = get_diagval_line_mmapped(mapped_diags, &array_r, val_index, r_offset);
    // r_offset = element_r.line_offset;
    while (element_r.value < pivot.value && r > 1)
    {
      r--;
      if (r < left)
        break;
      array_r = array[r];
      // printf(" array_r:%d\n", array_r); // DEBUG
      element_r = get_diagval_line_mmapped(mapped_diags, &array_r, val_index, element_r.line_offset);
      // r_offset = element_r.line_offset;
    }
    // free(element_r);
    r_offset = element_r.line_offset;

    // If l and r have not crossed, swap the elements in the array only
    if (l <= r)
    {
      int temp = array[l];
      array[l] = array[r];
      array[r] = temp;

      l++;
      r--;
    }
  }
  // exit(0); // DEBUG
  // free(pivot);

  // Debug output to trace recursive partitioning
  // printf("Partition complete: Left: %d, Right: %d, l: %d, r: %d\n", left, right, l, r);

  // printf(" l:%d r:%d\n", l, r);                                      // DEBUG
  // printf(" array[%d]:%d, array[%d]:%d\n", l, array[l], r, array[r]); // DEBUG

  // for (size_t i = 0; i < 24; i++)
  // {
  //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
  // }

  // Recur into left and right partitions
  if (left < r)
    quicksort_ow(array, left, r, mapped_diags, l_offset, r_offset, val_index, total_diagonals);
  if (l < right)
    quicksort_ow(array, l, right, mapped_diags, l_offset, r_offset, val_index, total_diagonals);
} /* quicksort_ow */

// // void quicksort_ow(struct multi_frag **array, int left, int right)
// void quicksort_ow(int **array, int left, int right, mmapped_file *mapped_diags, size_t l_offset, size_t r_offset, int *val_index)
// {
//   // if (left >= right)
//   //   return;

//   // int l = left, r = right;
//   // // struct multi_frag *element;
//   // int *element_idx;
//   // // element_idx = array[(left + right) / 2];
//   // element_idx = array[(left + right) / 2] + 1; // because anc_num is 1 indexed

//   // struct multi_frag *element = get_anchor_by_num(mapped_fasta ,mapped_anc, &element_idx, l_offset, 0);
//   // struct multi_frag *element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &l, l_offset, 0);
//   // struct multi_frag *element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &r, r_offset, 0);

//   // do
//   // {
//   //   element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &l, element_l->line_offset, 0);
//   //   element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &r, element_r->line_offset, 0);
//   //   // while (array[l]->ow > element->ow)
//   //   while (element_l->ow > element->ow && l <= r)
//   //   {
//   //     l++;
//   //     element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &l, element_l->line_offset, 0);
//   //   }
//   //   // while (element->ow > array[r]->ow)
//   //   while (element->ow > element_r->ow && r >= 1)
//   //   {
//   //     r--;
//   //     element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &r, element_r->prev_offset, 0);
//   //   }

//   //   if (l < r)
//   //   {
//   //     change_struct_el(array, l, r);
//   //     // change_struct_el_mmap(mapped_anc, element_l->line_offset, element_r->line_offset);
//   //     // change_struct_el_mmap(mapped_anc, element_l, element_r);
//   //   }
//   //   if (l <= r)
//   //   {
//   //     l++;
//   //     r--;
//   //   }
//   // } while (l <= r);

//   // if (left < r)
//   //   quicksort_ow(array, left, r, mapped_fasta, mapped_anc, element_l->line_offset, element_r->line_offset);
//   // if (l < right)
//   //   quicksort_ow(array, l, right, mapped_fasta, mapped_anc, element_l->line_offset, element_r->line_offset);

//   if (left >= right)
//     return;

//   int l = left;
//   int r = right;
//   int pivot_idx = array[(left + right) / 2];

//   int zero_idx = 0;

//   // for (size_t i = 0; i < total_anc_num; i++)
//   // {
//   //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   // }
//   // struct multi_frag *test_frg = frg_array;
//   // while(test_frg != NULL){
//   //   printf("test_frg->s[0]:%d test_frg->s[1]:%d\n", test_frg->s[0], test_frg->s[1]); //DEBUG
//   //   test_frg = test_frg->next;
//   // }
//   // printf(" frg_array[pivot_idx]:%p \n", frg_array[pivot_idx]); //DEBUG
//   printf(" left:%d right:%d pivot_idx:%d\n", left, right, pivot_idx); // DEBUG

//   // Get the pivot weight for comparison
//   // // struct multi_frag *pivot = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &pivot_idx, l_offset, 0);
//   // struct multi_frag *pivot = frg_array;
//   diag_value pivot = get_diagval_line_mmapped(mapped_diags, &pivot_idx, val_index, l_offset);
//   // printf(" pivot:%p \n", pivot); // DEBUG
//   // exit(0);                       // DEBUG
//   // int iter_idx = 0;
//   // while (pivot != NULL && iter_idx < pivot_idx)
//   // {
//   //   // printf("pivot->s[0]:%d pivot->s[1]:%d\n", pivot->s[0], pivot->s[1]); //DEBUG
//   //   pivot = pivot->next;
//   //   iter_idx++;
//   // }

//   // int pivot_s0 = pivot.s0;
//   // int pivot_s1 = pivot.s1;
//   // // ow_value pivot_ow = get_ow_mmapped(mapped_ow, &pivot_s0, &pivot_s1, &zero_idx, 0);
//   // ow_value pivot_ow = get_diagval_mmapped(mapped_diags, &pivot_s0, &pivot_s1, &zero_idx, 0);
//   // double pivot_val = pivot.value;

//   while (l <= r)
//   {
//     // Move 'l' right until we find an element >= pivot weight
//     int array_l = array[l]; // because anc_num is one indexed
//     // printf(" array_l:%d\n", array_l); // DEBUG
//     // // struct multi_frag *element_l = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &array_l, l_offset, 0);
//     // struct multi_frag *element_l = frg_array;
//     diag_value element_l = get_diagval_line_mmapped(mapped_diags, &array_l, val_index, l_offset);
//     l_offset = element_l.line_offset;
//     // int element_l_s0 = element_l.s0;
//     // int element_l_s1 = element_l.s1;
//     // ow_value element_l_ow = get_ow_mmapped(mapped_ow, &element_l_s0, &element_l_s0, &zero_idx, 0);
//     // double element_l_val = element_l.value;
//     // while (element_l->ow > pivot->ow && l < total_anc_num)
//     while (element_l.value > pivot.value && l < total_anc_num)
//     {
//       l++;
//       if (l > right)
//         break;
//       array_l = array[l];
//       // printf(" array_l:%d\n", array_l); // DEBUG
//       // // element_l = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &array_l, l_offset, 0);
//       // element_l = frg_array;
//       element_l = get_diagval_line_mmapped(mapped_diags, &array_l, val_index, l_offset);
//       l_offset = element_l.line_offset;
//       // iter_idx = 0;
//       // while (element_l != NULL && iter_idx < array_l)
//       // {
//       //   element_l = element_l->next;
//       //   iter_idx++;
//       // }
//       // element_l_s0 = element_l->s[0];
//       // element_l_s1 = element_l->s[1];
//       // element_l_ow = get_ow_mmapped(mapped_ow, &element_l_s0, &element_l_s1, &zero_idx, 0);
//     }
//     // free(element_l);

//     // Move 'r' left until we find an element <= pivot weight
//     int array_r = array[r]; // because anc_num is 1 indexed
//     // printf(" array_r:%d array[r]:%d\n", array_r, array[r]); // DEBUG
//     // // struct multi_frag *element_r = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &array_r, r_offset, 0);
//     // struct multi_frag *element_r = frg_array;
//     diag_value element_r = get_diagval_line_mmapped(mapped_diags, &array_r, val_index, r_offset);
//     r_offset = element_r.line_offset;
//     // iter_idx = 0;
//     // while (element_r != NULL && iter_idx < array_r)
//     // {
//     //   element_r = element_r->next;
//     //   iter_idx++;
//     //   // printf("iter_idx:%d\n", iter_idx); //DEBUG
//     // }
//     // int element_r_s0 = element_r->s[0];
//     // int element_r_s1 = element_r->s[1];
//     // ow_value element_r_ow = get_ow_mmapped(mapped_ow, &element_r_s0, &element_r_s1, &zero_idx, 0);
//     // while (element_r->ow < pivot->ow && r > 1)
//     while (element_r.value < pivot.value && r > 1)
//     {
//       r--;
//       if (r < left)
//         break;
//       array_r = array[r];
//       // printf(" array_r:%d\n", array_r); // DEBUG
//       // element_r = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &array_r, r_offset, 0);
//       // element_r = frg_array;
//       element_r = get_diagval_line_mmapped(mapped_diags, &array_r, val_index, r_offset);
//       r_offset = element_r.line_offset;
//       // iter_idx = 0;
//       // while (element_r != NULL && iter_idx < array_r)
//       // {
//       //   element_r = element_r->next;
//       //   iter_idx++;
//       // }
//       // element_r_s0 = element_r->s[0];
//       // element_r_s1 = element_r->s[1];
//       // element_r_ow = get_ow_mmapped(mapped_ow, &element_r_s0, &element_r_s1, &zero_idx, 0);
//     }
//     // free(element_r);

//     // If l and r have not crossed, swap the elements in the array only
//     if (l <= r)
//     {
//       int temp = array[l];
//       array[l] = array[r];
//       array[r] = temp;

//       l++;
//       r--;
//     }
//   }

//   // free(pivot);

//   // Debug output to trace recursive partitioning
//   printf("Partition complete: Left: %d, Right: %d, l: %d, r: %d\n", left, right, l, r);

//   printf(" l:%d r:%d\n", l, r);                                      // DEBUG
//   printf(" array[%d]:%d, array[%d]:%d\n", l, array[l], r, array[r]); // DEBUG

//   // for (size_t i = 0; i < 24; i++)
//   // {
//   //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   // }

//   // Recur into left and right partitions
//   if (left < r)
//     quicksort_ow(array, left, r, mapped_diags, l_offset, r_offset, val_index);
//   if (l < right)
//     quicksort_ow(array, l, right, mapped_diags, l_offset, r_offset, val_index);

// } /*quicksort_ow*/

// void quicksort_ow(struct multi_frag **frg_array, int left, int right)
// {
//   // if (left >= right)
//   //   return;

//   // int l = left, r = right;
//   // // struct multi_frag *element;
//   // int *element_idx;
//   // // element_idx = array[(left + right) / 2];
//   // element_idx = array[(left + right) / 2] + 1; // because anc_num is 1 indexed

//   // struct multi_frag *element = get_anchor_by_num(mapped_fasta ,mapped_anc, &element_idx, l_offset, 0);
//   // struct multi_frag *element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &l, l_offset, 0);
//   // struct multi_frag *element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &r, r_offset, 0);

//   // do
//   // {
//   //   element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &l, element_l->line_offset, 0);
//   //   element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &r, element_r->line_offset, 0);
//   //   // while (array[l]->ow > element->ow)
//   //   while (element_l->ow > element->ow && l <= r)
//   //   {
//   //     l++;
//   //     element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &l, element_l->line_offset, 0);
//   //   }
//   //   // while (element->ow > array[r]->ow)
//   //   while (element->ow > element_r->ow && r >= 1)
//   //   {
//   //     r--;
//   //     element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &r, element_r->prev_offset, 0);
//   //   }

//   //   if (l < r)
//   //   {
//   //     change_struct_el(array, l, r);
//   //     // change_struct_el_mmap(mapped_anc, element_l->line_offset, element_r->line_offset);
//   //     // change_struct_el_mmap(mapped_anc, element_l, element_r);
//   //   }
//   //   if (l <= r)
//   //   {
//   //     l++;
//   //     r--;
//   //   }
//   // } while (l <= r);

//   // if (left < r)
//   //   quicksort_ow(array, left, r, mapped_fasta, mapped_anc, element_l->line_offset, element_r->line_offset);
//   // if (l < right)
//   //   quicksort_ow(array, l, right, mapped_fasta, mapped_anc, element_l->line_offset, element_r->line_offset);

//   // if (left >= right)
//   //   return;

//   int l = left;
//   int r = right;

//   int zero_idx = 0;

//   // // for (size_t i = 0; i < total_anc_num; i++)
//   // // {
//   // //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   // // }
//   // // struct multi_frag *test_frg = frg_array;
//   // // while(test_frg != NULL){
//   // //   printf("test_frg->s[0]:%d test_frg->s[1]:%d\n", test_frg->s[0], test_frg->s[1]); //DEBUG
//   // //   test_frg = test_frg->next;
//   // // }
//   // // printf(" frg_array[pivot_idx]:%p \n", frg_array[pivot_idx]); //DEBUG
//   // printf(" left:%d right:%d pivot_idx:%d\n", left, right, pivot_idx); // DEBUG

//   // Get the pivot weight for comparison
//   // struct multi_frag *pivot = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &pivot_idx, l_offset, 0);

//   struct multi_frag *pivot = frg_array[(left + right) / 2];
//   // printf(" pivot:%p \n", pivot); // DEBUG

//   int pivot_s0 = pivot->s[0];
//   int pivot_s1 = pivot->s[1];
//   // ow_value pivot_ow = get_ow_mmapped(mapped_ow, &pivot_s0, &pivot_s1, &zero_idx, 0); //TODO: fix set_ow_mmapped() & get_ow_mmapped() and then convert this code to use those functions
//   printf(" left:%d right:%d pivot_idx:%d\n", left, right, (left + right) / 2); // DEBUG
//   do
//   {
//     while (frg_array[l]->ow > pivot->ow)
//       l++;
//     while (pivot->ow > frg_array[r]->ow)
//       r--;

//     if (l < r)
//       change_struct_el_frg(frg_array, l, r);
//     if (l <= r)
//     {
//       l++;
//       r--;
//     }
//   } while (l <= r);

//   // free(pivot);

//   // Debug output to trace recursive partitioning
//   // printf("Partition complete: Left: %d, Right: %d, l: %d, r: %d\n", left, right, l, r);

//   // printf(" l:%d r:%d\n", l, r);                                      // DEBUG
//   // printf(" array[%d]:%d, array[%d]:%d\n", l, frg_array[l], r, frg_array[r]); // DEBUG

//   // for (size_t i = 0; i < 24; i++)
//   // {
//   //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   // }

//   // Recur into left and right partitions
//   if (left < r)
//     quicksort_ow(frg_array, left, r);
//   if (l < right)
//     quicksort_ow(frg_array, l, right);

// } /*quicksort_ow*/

/***********************************************************
 *                    quicksort_weight()
 ***********************************************************/
// void quicksort_weight(struct multi_frag **array, int left, int right)
void quicksort_weight(int **array, int left, int right, mmapped_file *mapped_fasta, mmapped_file *mapped_anc, size_t l_offset, size_t r_offset)
{

  if (left >= right)
    return;

  int l = left;
  int r = right;
  int pivot_idx = array[(left + right) / 2];

  // for (size_t i = 0; i < total_anc_num; i++)
  // {
  //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
  // }

  printf(" left:%d right:%d pivot_idx:%d\n", left, right, pivot_idx); // DEBUG

  // Get the pivot weight for comparison
  struct multi_frag *pivot = get_anchor_by_num(mapped_fasta, mapped_anc, &pivot_idx, l_offset, 0);

  while (l <= r)
  {
    // Move 'l' right until we find an element >= pivot weight
    int array_l = array[l];           // because anc_num is one indexed
    printf(" array_l:%d\n", array_l); // DEBUG
    struct multi_frag *element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &array_l, l_offset, 0);
    while (element_l->weight > pivot->weight && l < total_anc_num)
    {
      l++;
      if (l > right)
        break;
      array_l = array[l];
      printf(" array_l:%d\n", array_l); // DEBUG
      element_l = get_anchor_by_num(mapped_fasta, mapped_anc, &array_l, l_offset, 0);
    }
    free(element_l);

    // Move 'r' left until we find an element <= pivot weight
    int array_r = array[r];                                 // because anc_num is 1 indexed
    printf(" array_r:%d array[r]:%d\n", array_r, array[r]); // DEBUG
    struct multi_frag *element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &array_r, r_offset, 0);
    while (element_r->weight < pivot->weight && r > 1)
    {
      r--;
      if (r < left)
        break;
      array_r = array[r];
      printf(" array_r:%d\n", array_r); // DEBUG
      element_r = get_anchor_by_num(mapped_fasta, mapped_anc, &array_r, r_offset, 0);
    }
    free(element_r);

    // If l and r have not crossed, swap the elements in the array only
    if (l <= r)
    {
      int temp = array[l];
      array[l] = array[r];
      array[r] = temp;

      l++;
      r--;
    }
  }

  free(pivot);

  // Debug output to trace recursive partitioning
  printf("Partition complete: Left: %d, Right: %d, l: %d, r: %d\n", left, right, l, r);

  printf(" l:%d r:%d\n", l, r);                                      // DEBUG
  printf(" array[%d]:%d, array[%d]:%d\n", l, array[l], r, array[r]); // DEBUG

  // for (size_t i = 0; i < 24; i++)
  // {
  //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
  // }

  // Recur into left and right partitions
  if (left < r)
    quicksort_weight(array, left, r, mapped_fasta, mapped_anc, l_offset, r_offset);
  if (l < right)
    quicksort_weight(array, l, right, mapped_fasta, mapped_anc, l_offset, r_offset);

} /*quicksort_weight*/

/***********************************************************
 *                    assemble_list()
 ***********************************************************/
// void assemble_list(struct multi_frag **array, struct multi_frag *dp, int number)
// void assemble_list(int **array, struct multi_frag *dp, int number)
// {
//   int i, index = 0;
//   for (i = 0; i < number - 1; i++)
//   {
//     if (dp->anc_num == array[i])
//       index = i;
//     // array[i]->next = array[i + 1];
//     change_struct_el(array, i + 1, i);
//   }

//   // array[number - 1]->next = 0; // array[number] = 0;
//   if (dp->anc_num == array[number - 1])
//     index = number - 1;

//   if (index != 0)
//     // change_first(array[0], array[index - 1]);
//     change_struct_el(array, 0, index - 1);
// } /* assemble_list */

void assemble_list(struct multi_frag **array, struct multi_frag *dp, int number)
{
  int i, index = 0;
  for (i = 0; i < number - 1; i++)
  {
    if (dp == array[i])
      index = i;
    array[i]->next = array[i + 1];
  }

  array[number - 1]->next = 0;
  if (dp == array[number - 1])
    index = number - 1;

  if (index != 0)
    change_first(array[0], array[index - 1]);
} /* assemble_list */

/***********************************************************
 *                      frag_sort()
 ************************************************************/
// void frag_sort_mmap(int number, int olw, mmapped_file *mapped_fasta, mmapped_file *mapped_anc, mmapped_file *mapped_ow, struct multi_frag *dp)
// {
//   printf(" here 15.10 frag sort\n"); // DEBUG

//   int i = 1;

//   // char anc_file_name[NAME_LEN];
//   // strcpy(anc_file_name, input_name);
//   // strcat(anc_file_name, ".anc");

//   // mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
//   // struct stat *anc_sb = &(mmapped_anc->sb);
//   // size_t file_size = anc_sb->st_size;

//   // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
//   // struct stat *sb_fasta = &(mmapped_fasta->sb);
//   // size_t file_size_fasta = sb_fasta->st_size;

//   // struct multi_frag **array;
//   // int **array;

//   // if ((array = (struct multi_frag **)calloc(number + 1, sizeof(struct multi_frag *))) == 0)
//   // if ((array = (struct multi_frag **)malloc((number + 1) * sizeof(struct multi_frag *))) == 0)
//   // if ((array = (int **)malloc((number + 1) * sizeof(int *))) == 0)
//   // {
//   //   printf(" problems with memory allocation for `all_clades'\n \n");
//   //   exit(1);
//   // }

//   // array[0] = dp;

//   // printf(" mapped_anc[0]:%c\n", mmapped_anc->mapped_file[0]); //DEBUG

//   // int first_anc_num = 1;

//   // struct multi_frag *dp = get_anchor_by_num(mapped_anc, sb_anc, &first_anc_num, 0, 1);

//   // array[0] = dp; // get_anchor_by_num(mmapped_anc, anc_sb, &number);

//   // int anc_num = 1;
//   // int anc_num = dp->anc_num;
//   // array[0] = anc_num;
//   // // struct multi_frag *dp = get_anchor_by_num(mapped_anc, sb_anc, &anc_num, 0, 1);

//   // // printf("dp: %p\n", dp);                 // DEBUG
//   // // printf("dp->next: %p\n", dp->next);     // DEBUG
//   // // printf("dp->ow: %f\n", dp->ow);         // DEBUG
//   // // printf("dp->weight: %f\n", dp->weight); // DEBUG
//   // // printf("dp->s[0]: %d\n", dp->s[0]);     // DEBUG
//   // // printf("dp->s[1]: %d\n", dp->s[1]);     // DEBUG
//   // // printf("dp->b[0]: %d\n", dp->b[0]);     // DEBUG
//   // // printf("dp->b[1]: %d\n", dp->b[1]);     // DEBUG
//   // // printf("dp->ext: %d\n", dp->ext);       // DEBUG

//   // // while (array[i - 1]->next)
//   // size_t input_offset = 0;
//   // // while (i <= number)
//   // while (i <= number)
//   // {
//   //   // printf("array[%d] (i -1): %p\n", i - 1, array[i - 1]->next); // DEBUG
//   //   // int im1 = i - 1;
//   //   anc_num++;
//   //   array[i] = anc_num;
//   //   // array[i - 1] = get_anchor_by_num(mapped_anc, sb_anc, &i, input_offset, 1);
//   //   // // printf("array[%d]->next offset:%d\n", i - 1, array[i - 1]->next->line_offset); // DEBUG
//   //   // // printf("array[%d] (i): %p\n", i, array[i - 1]->next); // DEBUG
//   //   // array[i] = array[i - 1]->next;
//   //   // // printf(" i:%d input_offset:%d array[%d]_offset:%d array[%d]_offset:%d\n", i, input_offset, i - 1, array[i - 1]->line_offset, i, array[i]->line_offset); // DEBUG
//   //   // input_offset = array[i]->line_offset;

//   //   // printf("array[%d] offset:%d\n", i, array[i]->line_offset);             // DEBUG
//   //   // printf("array[%d]->next offset:%d\n", i, array[i]->next->line_offset); // DEBUG
//   //   i++;
//   // }

//   // // DEBUG
//   // for (size_t i = 0; i < number + 1; i++)
//   // {
//   //   printf(" array[%d] anc_num:%d ->next:%p\n", i, array[i]->anc_num, array[i]->next); // DEBUG
//   // }

//   // printf(" number:%d total_anc_num:%d\n", number, total_anc_num); // DEBUG
//   // printf(" number:%d i:%d\n", number, i);                         // DEBUG

//   // printf("array max index (frag sort):%d number:%d\n", i, number); // DEBUG
//   // exit(0); //DEBUG

//   if (olw)
//   { // quicksort_ow
//     // int **array;
//     struct multi_frag **frg_array;
//     // struct multi_frag *frg_element;
//     int frgarr_size = 0;
//     // frg_element = dp;
//     // while (frg_element != NULL)
//     // {
//     // frg_element = frg_element->next;
//     // frgarr_size++;
//     // }

//     if ((frg_array = (struct multi_frag **)calloc(number + 1, sizeof(struct multi_frag *))) == 0)
//     {
//       printf(" olw: problems with memory allocation for `all_clades'\n \n");
//       exit(1);
//     }
//     // if ((array = (int **)malloc((frgarr_size + 1) * sizeof(int *))) == 0)
//     // {
//     //   printf(" weights: problems with memory allocation for `all_clades'\n \n");
//     //   exit(1);
//     // }
//     frg_array[0] = dp;
//     // array[0] = dp->anc_num;
//     while (frg_array[i - 1]->next)
//     {
//       frg_array[i] = frg_array[i - 1]->next;
//       // array[i] = frg_array[i]->anc_num;
//       i++;
//       frgarr_size++;
//     }
//     free(frg_array);

//     printf(" number:%d total_anc_num:%d\n", frgarr_size, total_anc_num); // DEBUG
//     printf(" number:%d i:%d\n", frgarr_size, i);                         // DEBUG

//     printf("array max index (frag sort):%d number:%d\n", i, frgarr_size); // DEBUG

//     printf(" here 15.10.1 frag sort - quicksort_ow()\n"); // DEBUG
//     quicksort_ow(frg_array, 0, frgarr_size);
//     // // quicksort_ow(array, 0, number, mapped_anc, 0, 0);
//     // // quicksort_ow(array, dp, 0, number - 1, mapped_fasta, mapped_anc, mapped_ow, 0, 0);
//     // quicksort_ow(array, dp, 0, frgarr_size - 1, mapped_fasta, mapped_anc, mapped_ow, 0, 0);
//     assemble_list(frg_array, dp, frgarr_size);
//     printf(" post assemble (frag sort)\n"); // DEBUG
//     struct multi_frag *test_frg = frg_array[0];
//     int zero_idx = 0;
//     i = 0;
//     int sz, so;
//     printf(" mapped_ow:%s\n", mapped_ow->file_name); // DEBUG
//     while (test_frg != NULL)
//     {
//       sz = test_frg->s[0];
//       so = test_frg->s[1];
//       double test_ow = get_ow_mmapped(mapped_ow, &sz, &so, &zero_idx, test_frg->ow_offset).value;
//       printf(" array[%d]:%d s1:%d s2:%d weight:%f test_ow:%g ow:%g \n", i, test_frg->anc_num, sz, so, test_frg->weight, test_ow, test_frg->ow); // DEBUG
//       test_frg = test_frg->next;
//       i++;
//     }

//     // for (size_t i = 0; i < number; i++)
//     // {
//     //   int iter_test = 0;
//     //   struct multi_frag *test_frg = dp;
//     //   while (test_frg != NULL && iter_test < array[i])
//     //   {
//     //     test_frg = test_frg->next;
//     //     iter_test++;
//     //   }
//     //   printf(" array[%d]:%d s1:%d s2:%d weight:%f ow:%f\n", i, array[i], test_frg->s[0], test_frg->s[1], test_frg->weight, test_frg->ow); // DEBUG
//     // }
//     // free(frg_array);

//   } // quicksort_ow
//   else
//   { // quicksort_weight
//     int **array;
//     if ((array = (int **)malloc((number + 1) * sizeof(int *))) == 0)
//     {
//       printf(" weights: problems with memory allocation for `all_clades'\n \n");
//       exit(1);
//     }
//     int anc_num = dp->anc_num;
//     array[0] = anc_num;
//     // struct multi_frag *dp = get_anchor_by_num(mapped_anc, sb_anc, &anc_num, 0, 1);

//     // printf("dp: %p\n", dp);                 // DEBUG
//     // printf("dp->next: %p\n", dp->next);     // DEBUG
//     // printf("dp->ow: %f\n", dp->ow);         // DEBUG
//     // printf("dp->weight: %f\n", dp->weight); // DEBUG
//     // printf("dp->s[0]: %d\n", dp->s[0]);     // DEBUG
//     // printf("dp->s[1]: %d\n", dp->s[1]);     // DEBUG
//     // printf("dp->b[0]: %d\n", dp->b[0]);     // DEBUG
//     // printf("dp->b[1]: %d\n", dp->b[1]);     // DEBUG
//     // printf("dp->ext: %d\n", dp->ext);       // DEBUG

//     // while (array[i - 1]->next)
//     size_t input_offset = 0;
//     // while (i <= number)
//     while (i <= number)
//     {
//       // printf("array[%d] (i -1): %p\n", i - 1, array[i - 1]->next); // DEBUG
//       // int im1 = i - 1;
//       anc_num++;
//       array[i] = anc_num;
//       // array[i - 1] = get_anchor_by_num(mapped_anc, sb_anc, &i, input_offset, 1);
//       // // printf("array[%d]->next offset:%d\n", i - 1, array[i - 1]->next->line_offset); // DEBUG
//       // // printf("array[%d] (i): %p\n", i, array[i - 1]->next); // DEBUG
//       // array[i] = array[i - 1]->next;
//       // // printf(" i:%d input_offset:%d array[%d]_offset:%d array[%d]_offset:%d\n", i, input_offset, i - 1, array[i - 1]->line_offset, i, array[i]->line_offset); // DEBUG
//       // input_offset = array[i]->line_offset;

//       // printf("array[%d] offset:%d\n", i, array[i]->line_offset);             // DEBUG
//       // printf("array[%d]->next offset:%d\n", i, array[i]->next->line_offset); // DEBUG
//       i++;
//     }

//     printf(" number:%d total_anc_num:%d\n", number, total_anc_num); // DEBUG
//     printf(" number:%d i:%d\n", number, i);                         // DEBUG

//     printf("array max index (frag sort):%d number:%d\n", i, number); // DEBUG

//     printf(" here 15.10.2 frag sort - quicksort_weight()\n"); // DEBUG
//     // quicksort_weight(array, 0, number, mapped_anc, 0, 0);
//     quicksort_weight(array, 0, number - 1, mapped_fasta, mapped_anc, mapped_ow, 0, 0);

//     // Swapping index 0 and 1 to complete sorting
//     change_struct_el_int(array, 0, 1);

//     i = 0;
//     // int j = 0;
//     size_t i_offset = 0; //, j_offset = 0;
//     struct multi_frag *element_i;
//     // struct multi_frag *element_j;
//     char anc_name_sorted[NAME_LEN];
//     strcpy(anc_name_sorted, input_name);
//     strcat(anc_name_sorted, ".sorted.anc");

//     FILE *anc_sorted_file = fopen(anc_name_sorted, "w");

//     for (i = 0; i < number; i++)
//     {
//       int anc_i = array[i];
//       // element_i = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &anc_i, i_offset, 0); //uncoment if not working
//       // if (olw)
//       // {
//       //   int iter_index = 0;
//       //   element_i = dp;
//       //   while (element_i != NULL && iter_index < anc_i)
//       //   {
//       //     element_i = element_i->next;
//       //     iter_index++;
//       //   }
//       // }
//       // else
//       // {
//       element_i = get_anchor_by_num(mapped_fasta, mapped_anc, mapped_ow, &anc_i, i_offset, 0); // uncoment if not working
//       // }
//       // element_j = get_anchor_by_num(mapped_anc, &anc_j, j_offset, 0);
//       i_offset = element_i->line_offset;
//       // j_offset = element_j->line_offset;
//       // change_struct_el_mmap(mapped_anc, element_i, element_j);
//       printf(" array[%d]: element_i->s[0]:%d element_i->s[1]:%d\n", i, element_i->s[0] + 1, element_i->s[1] + 1); // DEBUG
//       fprintf(anc_sorted_file, "%d\t%d\t%d\t%d\t%d\t%f\n", element_i->s[0] + 1, element_i->s[1] + 1, element_i->b[0], element_i->b[1], element_i->ext, element_i->weight);
//     }
//     fclose(anc_sorted_file);

//     printf(" post assemble (frag sort)\n"); // DEBUG

//     for (size_t i = 0; i < total_anc_num; i++)
//     {
//       printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//     }

//     // free(array);

//   } // quicksort_weight

//   // printf(" pre assemble (frag sort)\n"); // DEBUG
//   // if (olw)
//   // {
//   //   for (size_t i = 0; i < number; i++)
//   //   {
//   //     int iter_test = 0;
//   //     struct multi_frag *test_frg = dp;
//   //     while (test_frg != NULL && iter_test < array[i])
//   //     {
//   //       test_frg = test_frg->next;
//   //       iter_test++;
//   //     }
//   //     printf(" array[%d]:%d s1:%d s2:%d weight:%f ow:%f\n", i, array[i], test_frg->s[0], test_frg->s[1], test_frg->weight, test_frg->ow); // DEBUG
//   //   }
//   // }
//   // else
//   // {
//   //   for (size_t i = 0; i < total_anc_num; i++)
//   //   {
//   //     printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   //   }
//   // }

//   // i = 1;
//   // struct multi_frag* dp = get_anchor_by_num(mapped_anc, &i, 0,0);
//   // // assemble_list(array, dp, number + 1);

//   // int j = 0;
//   // for(i = 0; dp->anc_num != array[i]; i++)
//   // {
//   //   j++;
//   // }

//   // //Swapping index 0 with dp->anc_num to complete sorting
//   // change_struct_el(array, 0, j);

// } /* frag_sort mmap */

void frag_sort_mmap(int number, int olw, mmapped_file *mapped_fasta, mmapped_file *mapped_anc, mmapped_file *mapped_diags)
{
  printf(" here 15.10 frag sort\n"); // DEBUG

  int i = 1;

  // struct multi_frag_min *ret_val = (struct multi_frag_min *)calloc(1, sizeof(struct multi_frag_min));

  // ret_val->line_offset = -1;
  // ret_val->total_sum = -1;

  // char anc_file_name[NAME_LEN];
  // strcpy(anc_file_name, input_name);
  // strcat(anc_file_name, ".anc");

  // mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
  // struct stat *anc_sb = &(mmapped_anc->sb);
  // size_t file_size = anc_sb->st_size;

  // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  // struct stat *sb_fasta = &(mmapped_fasta->sb);
  // size_t file_size_fasta = sb_fasta->st_size;

  // struct multi_frag **array;
  int **array;

  // if ((array = (struct multi_frag **)calloc(number + 1, sizeof(struct multi_frag *))) == 0)
  // if ((array = (struct multi_frag **)malloc((number + 1) * sizeof(struct multi_frag *))) == 0)
  if ((array = (int **)malloc((number + 1) * sizeof(int *))) == 0)
  {
    printf(" problems with memory allocation for `all_clades'\n \n");
    exit(1);
  }

  // array[0] = dp;

  // printf(" mapped_anc[0]:%c\n", mmapped_anc->mapped_file[0]); //DEBUG

  // int first_anc_num = 1;

  // struct multi_frag *dp = get_anchor_by_num(mapped_anc, sb_anc, &first_anc_num, 0, 1);

  // array[0] = dp; // get_anchor_by_num(mmapped_anc, anc_sb, &number);

  // int anc_num = 1;
  int anc_num = 1; // dp->anc_num;

  // struct multi_frag *dp = get_anchor_by_num(mapped_anc, sb_anc, &anc_num, 0, 1);

  // printf("dp: %p\n", dp);                 // DEBUG
  // printf("dp->next: %p\n", dp->next);     // DEBUG
  // printf("dp->ow: %f\n", dp->ow);         // DEBUG
  // printf("dp->weight: %f\n", dp->weight); // DEBUG
  // printf("dp->s[0]: %d\n", dp->s[0]);     // DEBUG
  // printf("dp->s[1]: %d\n", dp->s[1]);     // DEBUG
  // printf("dp->b[0]: %d\n", dp->b[0]);     // DEBUG
  // printf("dp->b[1]: %d\n", dp->b[1]);     // DEBUG
  // printf("dp->ext: %d\n", dp->ext);       // DEBUG

  // while (array[i - 1]->next)
  size_t input_offset = 0;

  // // DEBUG
  // for (size_t i = 0; i < number + 1; i++)
  // {
  //   printf(" array[%d] anc_num:%d ->next:%p\n", i, array[i]->anc_num, array[i]->next); // DEBUG
  // }

  printf(" number:%d total_anc_num:%d\n", number, total_anc_num); // DEBUG
  printf(" number:%d i:%d\n", number, i);                         // DEBUG

  printf("array max index (frag sort):%d number:%d\n", i, number); // DEBUG
  // exit(0); //DEBUG

  if (olw)
  {
    int zero_idx = 0, one_idx = 1;
    int hv_subseqs = get_diagval_mmapped(mapped_diags, &zero_idx, &one_idx, &zero_idx, &istep, &one_idx, 0).value;
    array[0] = get_diagval_mmapped(mapped_diags, &zero_idx, &one_idx, &hv_subseqs, &istep, &zero_idx, 0).line_num;
    // array[0] = anc_num;
    while (i <= number)
    {
      anc_num++;
      array[i] = anc_num;
      i++;
    }
    printf(" here 15.10.1 frag sort - quicksort_ow()\n"); // DEBUG
    // quicksort_ow(array, 0, number, mapped_anc, 0, 0);
    // quicksort_ow(array, 0, number - 1, mapped_fasta, mapped_anc, mapped_ow, 0, 0);
    int ow_idx = 8;
    int total_diagonals = get_seqcount_mmapped_fasta(mapped_diags);
    quicksort_ow(array, 0, number - 1, mapped_diags, 0, 0, &ow_idx, &total_diagonals);
  }
  else
  {
    array[0] = anc_num;
    while (i <= number)
    {
      anc_num++;
      array[i] = anc_num;
      i++;
    }
    printf(" here 15.10.2 frag sort - quicksort_weight()\n"); // DEBUG
    // quicksort_weight(array, 0, number, mapped_anc, 0, 0);
    quicksort_weight(array, 0, number - 1, mapped_fasta, mapped_anc, 0, 0);
  }

  // Swapping index 0 and 1 to complete sorting
  change_struct_el_int(array, 0, 1);

  // printf(" pre assemble (frag sort)\n"); // DEBUG
  // for (size_t i = 0; i < total_anc_num; i++)
  // {
  //   printf(" array[%d]:%d\n", i, array[i]); // DEBUG
  // }

  i = 0;
  // int j = 0;
  size_t i_offset = 0; //, j_offset = 0;

  // struct multi_frag *element_j;
  char name_sorted[NAME_LEN];
  FILE *sorted_file;
  if (olw)
  {
    // char diags_name_sorted[NAME_LEN];
    // strcpy(name_sorted, mapped_diags->file_name);
    strcpy(&name_sorted, &input_name);
    strcat(name_sorted, ".diags.sorted");

    sorted_file = fopen(name_sorted, "a+");
  }
  else
  {
    // char anc_name_sorted[NAME_LEN];
    strcpy(&name_sorted, &input_name);
    strcat(name_sorted, ".sorted.anc");

    sorted_file = fopen(name_sorted, "w");
  }

  for (i = 0; i < number; i++)
  {
    int anc_i = array[i];
    if (olw)
    {

      diag_line element_i = get_diagvals_line_mmapped(mapped_diags, &anc_i, i_offset);
      i_offset = element_i.line_offset;
      // j_offset = element_j->line_offset;
      // change_struct_el_mmap(mapped_anc, element_i, element_j);
      printf(" array[%d]: offset:%u element_i->s[0]:%d element_i->s[1]:%d\n", i, i_offset, element_i.s0, element_i.s1); // DEBUG
      char *file_line = (char *)calloc(1, ((3 * sizeof(short)) + (9 * sizeof(int)) + (3 * sizeof(double)) + (17 * sizeof(char))) + 1);
      // strcpy(file_line, "");
      sprintf(file_line, ">%d,%d,%d,%d:%d,%d,%d,%d,%u,%u,%d,%g,%g,%g,%u\n", element_i.s0, element_i.s1, element_i.hv, istep, element_i.anc_num, element_i.numsubseq, element_i.b0, element_i.b1, element_i.sel, element_i.trans, element_i.ext, element_i.weight, element_i.ow, element_i.total_sum, element_i.cs);
      // printf("file_line:%s\n", file_line); // DEBUG
      // exit(0);                             // DEBUG
      fprintf(sorted_file, "%s", file_line);
      free(file_line);
    }
    else
    {
      struct multi_frag *element_i = get_anchor_by_num(mapped_fasta, mapped_anc, &anc_i, i_offset, 0);
      // element_j = get_anchor_by_num(mapped_anc, &anc_j, j_offset, 0);
      i_offset = element_i->line_offset;
      // j_offset = element_j->line_offset;
      // change_struct_el_mmap(mapped_anc, element_i, element_j);
      printf(" file:%s array[%d]: element_i->s[0]:%d element_i->s[1]:%d element_i->wgt:%f\n", mapped_anc->file_name, i, element_i->s[0] + 1, element_i->s[1] + 1, element_i->weight); // DEBUG
      fprintf(sorted_file, "%d\t%d\t%d\t%d\t%d\t%f\n", element_i->s[0] + 1, element_i->s[1] + 1, element_i->b[0], element_i->b[1], element_i->ext, element_i->weight);
    }
  }
  fclose(sorted_file);
  // exit(0);
  // strcpy(filename, name_sorted);

  // if(olw){
  //   // iterate through array and trace the path of the diagonals
  //   i_offset = 0;
  //   for (i = 0; i < number; i++)
  //   {
  //     diag_line element_arr = get_diagvals_line_mmapped(mapped_diags, &array[i], i_offset);
  //     i_offset = element_arr.line_offset;
  //     ret_val->idx = anc_i;
  //     ret_val->hv = element_i.hv;
  //     ret_val->istep = element_i.istep;
  //     ret_val->s0 = element_i.s0;
  //     ret_val->s1 = element_i.s1;
  //   }
  // }

  // i = 1;
  // struct multi_frag* dp = get_anchor_by_num(mapped_anc, &i, 0,0);
  // // assemble_list(array, dp, number + 1);

  // int j = 0;
  // for(i = 0; dp->anc_num != array[i]; i++)
  // {
  //   j++;
  // }

  // //Swapping index 0 with dp->anc_num to complete sorting
  // change_struct_el(array, 0, j);

  free(array);

  if (olw)
  {
    // munmap(mapped_diags->mapped_file, mapped_diags->sb.st_size);
    // mapped_diags = mmap_file(name_sorted, mapped_diags->file_mode, mapped_diags->protocol, mapped_diags->map_mode, 0);
    mmapped_file *tmp_diags = mmap_file(name_sorted, mapped_diags->file_mode, mapped_diags->protocol, mapped_diags->map_mode, 0);

    printf(" post assemble (frag sort)\n"); // DEBUG
    for (size_t i = 0; i < number; i++)
    {
      diag_line test_line = get_diagvals_line_mmapped(mapped_diags, &i, 0);
      printf(" array[%d]:%d s1:%d s2:%d weight:%f ow:%f test_ow:%g\n", i, array[i], test_line.s0, test_line.s1, test_line.weight, test_line.ow); // DEBUG
    }
    munmap(tmp_diags->mapped_file, tmp_diags->sb.st_size);
    free(tmp_diags->file_name);
    free(tmp_diags);
  }
  else
  {
    munmap(mapped_anc->mapped_file, mapped_anc->sb.st_size);
    mapped_anc = mmap_file(name_sorted, mapped_anc->file_mode, mapped_anc->protocol, mapped_anc->map_mode, 0);
    printf(" post assemble (frag sort)\n"); // DEBUG
    for (size_t i = 0; i < total_anc_num; i++)
    {
      struct multi_frag *test_line = get_anchor_by_num(mapped_fasta, mapped_anc, &i, i_offset, 0);
      // printf(" test_line:%p\n", test_line);                                                                                                                        // DEBUG
      printf(" array[%d]:%d s1:%d s2:%d weight:%f anc_num:%d test_ow:%g\n", i, array[i], test_line->s[0], test_line->s[1], test_line->anc_num, test_line->weight); // DEBUG
    }
  }
} /* frag_sort mmap */

// void frag_sort(mmapped_file *mmapped_fasta, mmapped_file *mmapped_anc, mmapped_file *mmapped_ow, int number, struct multi_frag *dp, int olw)
// {
//   int i = 1;

//   // char anc_file_name[NAME_LEN];
//   // strcpy(anc_file_name, input_name);
//   // strcat(anc_file_name, ".anc");

//   // mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
//   // // struct stat *anc_sb = &(mmapped_anc->sb);
//   // size_t file_size = mmapped_anc->sb.st_size;

//   struct multi_frag **array;
//   // if ((array = (struct multi_frag **)calloc(number + 1, sizeof(struct multi_frag *))) == 0)
//   if ((array = (struct multi_frag **)malloc(number + 1 * sizeof(struct multi_frag *))) == 0)
//   {
//     printf(" problems with memory allocation for `all_clades'\n \n");
//     exit(1);
//   }

//   array[0] = dp;

//   // printf("dp: %p\n", dp);                 // DEBUG
//   // printf("dp->next: %p\n", dp->next);     // DEBUG
//   // printf("dp->ow: %f\n", dp->ow);         // DEBUG
//   // printf("dp->weight: %f\n", dp->weight); // DEBUG
//   // printf("dp->s[0]: %d\n", dp->s[0]);     // DEBUG
//   // printf("dp->s[1]: %d\n", dp->s[1]);     // DEBUG
//   // printf("dp->b[0]: %d\n", dp->b[0]);     // DEBUG
//   // printf("dp->b[1]: %d\n", dp->b[1]);     // DEBUG
//   // printf("dp->ext: %d\n", dp->ext);       // DEBUG

//   int next_anc = dp->anc_num;

//   // while (array[i - 1]->next)
//   while (i <= number)
//   {
//     // printf("array[%d] (i -1): %p\n", i - 1, array[i - 1]->next); // DEBUG
//     // int next_anc = array[i - 1]->anc_num + 1;

//     next_anc = array[i - 1]->anc_num - 1;
//     // array[i] = array[i - 1]->next;
//     array[i] = get_anchor_by_num(mmapped_fasta, mmapped_anc, mmapped_ow, &next_anc, array[i - 1]->line_offset, 1);
//     // printf("array[%d] (i): %p\n", i, array[i]->next); // DEBUG
//     i++;
//   }

//   if (olw)
//     // quicksort_ow(array, 0, number, mmapped_anc, 0, 0);
//     quicksort_ow(array, 0, number - 1, mmapped_fasta, mmapped_anc, mmapped_ow, 0, 0);
//   else
//     // quicksort_weight(array, 0, number, mmapped_anc, 0, 0);
//     quicksort_weight(array, 0, number - 1, mmapped_fasta, mmapped_anc, mmapped_ow, 0, 0);

//   // Swapping index 0 and 1 to complete sorting
//   change_struct_el(array, 0, 1);

//   // printf(" pre assemble (frag sort)\n"); //DEBUG
//   //   for (size_t i = 0; i < total_anc_num; i++)
//   //   {
//   //     printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   //   }

//   i = 0;
//   // int j = 0;
//   size_t i_offset = 0; //, j_offset = 0;
//   struct multi_frag *element_i;
//   // struct multi_frag *element_j;
//   char anc_name_sorted[NAME_LEN];
//   strcpy(anc_name_sorted, input_name);
//   strcat(anc_name_sorted, ".sorted.anc");

//   FILE *anc_sorted_file = fopen(anc_name_sorted, "w");

//   for (i = 0; i < number; i++)
//   {
//     int anc_i = array[i];
//     element_i = get_anchor_by_num(mmapped_fasta, mmapped_anc, mmapped_ow, &anc_i, i_offset, 0);
//     // element_j = get_anchor_by_num(mapped_anc, &anc_j, j_offset, 0);
//     i_offset = element_i->line_offset;
//     // j_offset = element_j->line_offset;
//     // change_struct_el_mmap(mapped_anc, element_i, element_j);
//     printf(" array[%d]: element_i->s[0]:%d element_i->s[1]:%d\n", i, element_i->s[0] + 1, element_i->s[1] + 1); // DEBUG
//     fprintf(anc_sorted_file, "%d\t%d\t%d\t%d\t%d\t%f\n", element_i->s[0] + 1, element_i->s[1] + 1, element_i->b[0], element_i->b[1], element_i->ext, element_i->weight);
//   }
//   fclose(anc_sorted_file);

//   // assemble_list(array, dp, number + 1);

//   int j = 0;
//   for (i = 0; dp->anc_num != array[i]; i++)
//   {
//     j++;
//   }

//   // Swapping index 0 with dp->anc_num to complete sorting
//   change_struct_el(array, 0, j);

//   // printf(" post assemble (frag sort)\n"); //DEBUG
//   //   for (size_t i = 0; i < total_anc_num; i++)
//   //   {
//   //     printf(" array[%d]:%d\n", i, array[i]); // DEBUG
//   //   }

// } /* frag_sort */

void ow_add_mmap(int total_diagonals, mmapped_file *mapped_diags, mmapped_file *mapped_fasta)
{

  int sm1_s0 = 0, zero_idx = 0;
  int sm1_s1 = sm1_s0 + 1;
  int sm2_s0 = sm1_s0 + 1;
  int sm2_s1 = sm2_s0 + 1;
  int diagcnt_idx = 1, b0_idx = 2, b1_idx = 3, sel_idx = 4, ext_idx = 6, trans_idx = 5, weight_idx = 7, ow_idx = 8, sum_idx = 9, cs_idx = 10; // 0-indexed
  diag_value diag_count_sm1_val = get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &zero_idx, &istep, &diagcnt_idx, 0);
  int diag_count_sm1 = (int)diag_count_sm1_val.value - 1;
  // int diag_count_sm1 = (int)diag_count_sm1_val.value;
  size_t diag_offset_sm1 = diag_count_sm1_val.line_offset;
  size_t input_offset_s1 = 0, input_offset_s2 = 0;
  int dia_counter = 0;

  // int total_diagonals = get_seqcount_mmapped_fasta(mapped_diags);

  for (sm1_s0 = 0; sm1_s0 < seqnum; sm1_s0++)
  {
    // if (sm1_s0 == 2)
    // {
    //   exit(0); // DEBUG
    // }
    for (sm1_s1 = sm1_s0 + 1; sm1_s1 < seqnum; sm1_s1++)
    {
      // if (i == j)
      //   continue;

      diag_count_sm1_val = get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &zero_idx, &istep, &diagcnt_idx, diag_offset_sm1);
      diag_count_sm1 = (int)diag_count_sm1_val.value;
      diag_offset_sm1 = diag_count_sm1_val.line_offset;

      printf("===== PROCESS : sm1_s0:%d sm1_s1:%d\n", sm1_s0, sm1_s1); // DEBUG
      // printf("diag count sm1:%d\n", diag_count_sm1);             // DEBUG
      // exit(0);                                                   // DEBUG
      int sm1_trans = (int)(get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &zero_idx, &istep, &trans_idx, 0).value);
      int conslen;

      for (int subdiag_sm1 = 0; subdiag_sm1 < diag_count_sm1; subdiag_sm1++)
      {
        dia_counter++;
        int sm1_b0 = (int)(get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &b0_idx, diag_offset_sm1).value);
        int sm1_b1 = (int)(get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &b1_idx, diag_offset_sm1).value);
        short sm1_ext = (short)(get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &ext_idx, diag_offset_sm1).value);
        double sm1_weight = (double)(get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &weight_idx, diag_offset_sm1).value);
        // printf(" sm1_s0:%d sm1_s1:%d sm1_weight:%g\n", sm1_s0, sm1_s1, sm1_weight); // DEBUG
        // if (sm1_weight <= 0)
        // {
        //   printf("SM 1 WEIGHT 0\n"); // DEBUG
        //   // exit(0);                   // DEBUG
        //   // continue;
        // }
        FILE *fp_st;
        if (print_status)
          if ((dia_counter % 100) == 0)
          {
            fp_st = fopen(pst_name, "w");

            fprintf(fp_st, " dsd  %s \n", input_line);
            fprintf(fp_st, "\n\n\n    Status of the program run:\n");
            fprintf(fp_st, "    ==========================\n\n");
            if (seqnum > 2)
            {
              fprintf(fp_st, "      iteration step %d in ", istep);
              fprintf(fp_st, "multiple alignment\n");
            }
            fprintf(fp_st, "      calculating overlap weight for diagonals\n");
            fprintf(fp_st, "      current diagonal = %d\n\n", dia_counter);
            fprintf(fp_st, "      total number of");
            fprintf(fp_st, " diagonals: %d\n\n\n\n", total_diagonals);
            fclose(fp_st);
          }

        int sm1_s[2] = {sm1_s0, sm1_s1};
        int sm1_b[2] = {sm1_b0, sm1_b1};
        input_offset_s1 = 0;
        for (sm2_s1 = sm1_s1 + 1; sm2_s1 < seqnum; sm2_s1++)
        {
          // sm2_s0 = sm1_s0 + 1;
          sm2_s0 = sm1_s1;

          // if ((sm1_s0 == sm2_s0 && sm1_s1 == sm2_s1) || (sm1_s0 == sm2_s1 && sm1_s1 == sm2_s0))
          // {
          //   continue;
          // }

          diag_value diag_count_sm2_val = get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &zero_idx, &istep, &diagcnt_idx, 0);
          // int diag_count_sm1 = (int)diag_count_sm1_val.value - 1;
          int diag_count_sm2 = (int)diag_count_sm2_val.value;
          size_t diag_offset_sm2 = diag_count_sm2_val.line_offset;

          int sm2_s[2] = {sm2_s0, sm2_s1};

          // printf("(SUB) PROCESS : sm2_s0:%d sm2_s1:%d\n", sm2_s[0], sm2_s[1]); // DEBUG
          // printf("diag count sm2:%d\n", diag_count_sm2);                       // DEBUG

          // size_t diag_offset_sm2 = get_diagval_mmapped(mmapped_diags, &sm2_s0, &sm2_s1, &diag_count, &istep, &diagcnt_idx, diag_offset_sm1).line_offset;

          // for (int subdiag_sm1 = 0; subdiag_sm1 < diag_count_sm1; subdiag_sm1++)
          // {
          // int sm1_b0 = (int)(get_diagval_mmapped(mmapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &b0_idx, diag_offset_sm1).value);
          // int sm1_b1 = (int)(get_diagval_mmapped(mmapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &b1_idx, diag_offset_sm1).value);
          // short sm1_ext = (short)(get_diagval_mmapped(mmapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &ext_idx, diag_offset_sm1).value);
          // double sm1_weight = (double)(get_diagval_mmapped(mmapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &weight_idx, diag_offset_sm1).value);
          // // double sm1_ow = (double)(get_diagval_mmapped(mmapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &ow_idx, diag_offset_sm1).value);
          // int sm1_b[2] = {sm1_b0, sm1_b1};

          // for (int subdiag_sm2 = 0; subdiag_sm2 < diag_count_sm2; subdiag_sm2++)

          int subdiag_val = (subdiag_sm1 * 2) - (subdiag_sm1 % 2);
          // printf(" sm2_s0:%d sm2_s1:%d subdiag_sm1:%d subdiag_val:%d\n", sm2_s0, sm2_s1, subdiag_sm1, subdiag_val); // DEBUG
          // printf("    diag_count_sm1:%d diag_count_sm2:%d\n", diag_count_sm1, diag_count_sm2);                      // DEBUG
          if (diag_count_sm1 > diag_count_sm2)
          {
            // printf("SM 2 : diag_count_sm1 > diag_count_sm2\n"); // DEBUG
            // sm2_weight = sm2_ow;
            // exit(0);                   // DEBUG
            // continue;
            // printf(" sm2_s0:%d sm2_s1:%d subdiag_sm1:%d subdiag_val:%d\n", sm2_s0, sm2_s1, subdiag_sm1, subdiag_val); // DEBUG
            // printf("    diag_count_sm1:%d diag_count_sm2:%d\n", diag_count_sm1, diag_count_sm2);                      // DEBUG
            if (subdiag_val % diag_count_sm2 == 0)
              subdiag_val = subdiag_val % diag_count_sm2;
            else
              subdiag_val = subdiag_val % diag_count_sm2 - 1;
            // printf("subdiag_val:%d\n", subdiag_val); // DEBUG
          }
          // if ((subdiag_val % 2) > 0)
          //   subdiag_val = 0;
          // if (subdiag_val > diag_count_sm2 && subdiag_val <= diag_count_sm1)
          //   subdiag_val = (subdiag_val % diag_count_sm2) - 1;
          // printf("subdiag_sm1:%d\n", subdiag_sm1); // DEBUG
          for (int subdiag_sm2 = subdiag_val; subdiag_sm2 < diag_count_sm2; subdiag_sm2++) //|| subdiag_sm2 < diag_count_sm1
          {
            if (subdiag_sm2 == diag_count_sm2)
              continue;

            // if (subdiag_sm2 > diag_count_sm2 && subdiag_sm2 <= diag_count_sm1)
            //   subdiag_sm2 = (subdiag_sm2 % diag_count_sm2) - (subdiag_sm2 % 2);

            // printf("subdiag_sm2:%d\n", subdiag_sm2); // DEBUG
            double sm1_ow = (double)(get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &ow_idx, diag_offset_sm1).value);
            diag_value sm2_ow_val = get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &ow_idx, diag_offset_sm2);
            double sm2_ow = (double)(sm2_ow_val.value);
            diag_offset_sm2 = sm2_ow_val.line_offset;
            int sm2_b0 = (int)(get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &b0_idx, diag_offset_sm2).value);
            int sm2_b1 = (int)(get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &b1_idx, diag_offset_sm2).value);
            short sm2_ext = (short)(get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &ext_idx, diag_offset_sm2).value);
            double sm2_weight = (double)(get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &weight_idx, diag_offset_sm2).value);
            // double sm2_ow = (double)(get_diagval_mmapped(mmapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &ow_idx, diag_offset_sm2).value);
            int sm2_b[2] = {sm2_b0, sm2_b1};
            // printf(" sm2_s0:%d sm2_s1:%d sm2_weight:%g \n", sm2_s0, sm2_s1, sm2_weight); // DEBUG
            // if (sm2_weight <= 0)
            // {
            //   printf("SM 2 WEIGHT 0\n"); // DEBUG
            //   // sm2_weight = sm2_ow;
            //   // exit(0);                   // DEBUG
            //   // continue;
            //   printf(" sm2_s0:%d sm2_s1:%d sm2_weight:%g sm2_ow:%g subdiag_sm1:%d subdiag_sm2:%d\n", sm2_s0, sm2_s1, sm2_weight, sm2_ow, subdiag_sm1, subdiag_sm2); // DEBUG
            //   printf("    diag_count_sm1:%d diag_count_sm2:%d\n", diag_count_sm1, diag_count_sm2);                                                                  // DEBUG
            // }

            // printf("sm2_ow:%g sm2_weight:%g\n", sm2_ow, sm2_weight); // DEBUG
            input_offset_s2 = 0;
            for (int k = 0; k < 2; k++)
            {
              for (int l = 0; l < 2; l++)
              {
                if (sm1_s[k] == sm2_s[l])
                  if (sm1_s[l] != sm2_s[k])
                    if (sm1_b[k] < sm2_b[l] + sm2_ext &&
                        sm2_b[l] < sm1_b[k] + sm1_ext)
                    {
                      conslen = mini2(sm1_b[k] + sm1_ext, sm2_b[l] + sm2_ext) - maxi2(sm1_b[k], sm2_b[l]);
                      if (
                          (sm1_trans == 0) ||
                          ((conslen % 3) == 0))
                      {
                        int s1 = sm1_s[(k + 1) % 2];
                        int s2 = sm2_s[(l + 1) % 2];

                        // if (!input_offset_s1)
                        //   input_offset_s1 = 0;
                        // printf(" input_offset_s1: %d\n", input_offset_s1); // DEBUG
                        fasta_value *lhs_seq = get_seq_mmapped_fasta(mapped_fasta, &s1, input_offset_s1);
                        // lhs_seq->data++;
                        char *match_lhs = lhs_seq->data;
                        // if (lhs_seq->line_offset)
                        if (input_offset_s1 != 0)
                          input_offset_s1 = lhs_seq->line_offset - strlen(match_lhs) - 1;
                        // else
                        //   input_offset_s1 = 0;

                        // printf(" lhs seq: %s\n", match_lhs); // DEBUG

                        // if (!input_offset_s2)
                        //   input_offset_s2 = 0;
                        // printf(" input_offset_s2: %d\n", input_offset_s2); // DEBUG
                        fasta_value *rhs_seq = get_seq_mmapped_fasta(mapped_fasta, &s2, input_offset_s2);
                        // rhs_seq->data++;
                        char *match_rhs = rhs_seq->data;
                        // if (rhs_seq->line_offset > (size_t)0)
                        if (input_offset_s2 != 0)
                          input_offset_s2 = rhs_seq->line_offset - strlen(match_rhs) - 1;
                        // else
                        //   input_offset_s2 = 0;
                        // printf(" rhs seq: %s\n", match_rhs); // DEBUG

                        int b1 = sm1_b[(k + 1) % 2];
                        int dif = sm2_b[l] - sm1_b[k];
                        if (dif > 0)
                          b1 = b1 + dif;

                        int b2 = sm2_b[(l + 1) % 2];
                        dif = sm1_b[k] - sm2_b[l];
                        if (dif > 0)
                          b2 = b2 + dif;

                        int match = 0;
                        double add_wgt;

                        for (int cons_iter = 0; cons_iter < conslen; cons_iter++)
                        {
                          if (
                              (wgt_type == 0) ||
                              (sm1_trans && ((cons_iter % 3) == 0)))
                            match = match + sim_score[amino[s1][b1 + cons_iter]][amino[s2][b2 + cons_iter]];
                          else
                          { // match = match + ( seq[ s1 ][ b1 + k ] == seq[ s2 ][ b2 + k ] );
                            match = match + (match_lhs[b1 + cons_iter] == match_rhs[b2 + cons_iter]);
                          }
                        }

                        if (wgt_type == 0)
                          add_wgt = wgt_prot[conslen][match];
                        else if (sm1_trans)
                          add_wgt = wgt_trans[conslen / 3][match];
                        else
                          add_wgt = wgt_dna[conslen][match];

                        printf("SETTING (sm1i[%d]),sm1j[%d] sm2i[%d],(sm2j[%d])\n", sm1_s[k], sm1_s[l], sm2_s[k], sm2_s[l]); // DEBUG
                        printf("  (b) s1:%d s2:%d \n", (k + 1) % 2, (l + 1) % 2);                                            // DEBUG
                        printf("  (b) sm1->weight:%f \n", sm1_weight);                                                       // DEBUG
                        printf("  (b) sm2->weight:%f \n", sm2_weight);                                                       // DEBUG
                        printf("  (b) sm1->ow:%g \n", sm1_ow);                                                               // DEBUG
                        printf("  (b) sm2->ow:%g \n", sm2_ow);                                                               // DEBUG

                        sm1_ow = sm1_ow + add_wgt;
                        sm2_ow = sm2_ow + add_wgt;

                        printf("  (a) sm1->ow:%g \n", sm1_ow);                                                                      // DEBUG
                        printf("  (a) sm2->ow:%g \n", sm2_ow);                                                                      // DEBUG
                        printf("  sm1->s[0]:%d sm1->s[1]:%d hv:%d diag:cnt:%d\n", sm1_s[0], sm1_s[1], subdiag_sm1, diag_count_sm1); // DEBUG
                        printf("  sm2->s[0]:%d sm2->s[1]:%d hv:%d diag:cnt:%d\n", sm2_s[0], sm2_s[1], subdiag_sm2, diag_count_sm2); // DEBUG
                        printf("  add_wgt:%g istep:%d\n", add_wgt, istep);                                                          // DEBUG

                        set_diags_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &ow_idx, &sm1_ow, diag_offset_sm1);
                        printf(" SET sm1 ow[%d][%d]:%g\n", sm1_s0, sm1_s1, get_diagval_mmapped(mapped_diags, &sm1_s0, &sm1_s1, &subdiag_sm1, &istep, &ow_idx, diag_offset_sm1).value);

                        set_diags_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &ow_idx, &sm2_ow, diag_offset_sm2);
                        printf(" SET sm2 ow[%d][%d]:%g\n", sm2_s0, sm2_s1, get_diagval_mmapped(mapped_diags, &sm2_s0, &sm2_s1, &subdiag_sm2, &istep, &ow_idx, diag_offset_sm2).value);

                        free(match_lhs);
                        free(match_rhs);
                        // exit(0); // DEBUG
                      }
                    }
              }
            }
          }
        }
        // }

        // exit(0); // DEBUG
      }

      // diag_count_sm1_val = get_diagval_mmapped(mmapped_diags, &sm1_s0, &sm1_s1, &diag_count_sm1, &istep, &diagcnt_idx, diag_offset_sm1);
      // diag_count_sm1 = (int)diag_count_sm1_val.value;
      // diag_offset_sm1 = diag_count_sm1_val.line_offset;
      // exit(0); // DEBUG
    }
    // exit(0); // DEBUG
  }
  // exit(0); // DEBUG

  msync(mapped_diags->mapped_file, mapped_diags->sb.st_size, MS_SYNC);
}

// void seq_shift()
// {
//   int i, hv;

//   for (i = 0; i < seqnum; i++)
//     for (hv = seqlen[i] + 1; hv > 0; hv--)
//       seq[i][hv] = seq[i][hv - 1];
// }

// void filter_all(int *number)
// {
//   char anc_file_name[NAME_LEN];
//   strcpy(anc_file_name, input_name);
//   strcat(anc_file_name, ".sorted.anc");

//   mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
//   // struct stat *anc_sb = &(mmapped_anc->sb);
//   size_t anc_file_size = mmapped_anc->sb.st_size;

//   // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
//   // struct stat *sb_fasta = &(mmapped_fasta->sb);
//   // size_t file_size_fasta = sb_fasta->st_size;

//   int anc_one = 1;
//   struct multi_frag *diag = get_anchor_by_num(mmapped_anc, &anc_one, 0, 0); // its not the first anchor but the last in the diagonal

//   munmap(mmapped_anc->mapped_file, anc_file_size);
//   free(mmapped_anc);

//   // munmap(mmapped_fasta->mapped_file, file_size_fasta);
//   // free(mmapped_fasta);

//   filter(number, diag);
// }

void filter(int *number, int anchor_step, mmapped_file *mmapped_anc, mmapped_file *mmapped_fasta, mmapped_file *mmapped_diags, mmapped_file *mmapped_ow)
{
  /* checks diagonals one by one, if they are consistent with the
     diagonals already included into the alignment. If a new diagonal
     is consistent, it is included into the alignment and the frontiers
     in clos (when GABIOS is used) are changed accordingly */

  int i, j, k, l, sv, hv, ab[2], as[2], ae[2], aext, nv;
  double awgt;

  int test;      /* = 1 if current diagonal consistent; = 0 otherwise */
  int number_bf; /* number of diagonals before filter */

  FILE *fp_st, *fp_cap;

  FILE *open_pos_f;
  //  char *open_pos_fname; // = seq_file;
  int seq_file_length = sizeof(seq_file) / sizeof(char);
  char open_pos_fname[seq_file_length]; // = seq_file;
  snprintf(open_pos_fname, sizeof(seq_file), "%s", seq_file);
  strcat(open_pos_fname, ".op");

  // char all_diag_fname[seq_file_length]; // = seq_file;
  // strcpy(&all_diag_fname, &input_name);
  // strcat(all_diag_fname, ".diags");

  // mmapped_file *all_diags;

  // if (anchor_step == 0)
  //   all_diags = mmap_file(&all_diag_fname, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);

  int anc_idx = 0, diagcnt_idx = 1, b0_idx = 2, b1_idx = 3, ext_idx = 6, sel_idx = 4, trans_idx = 5, weight_idx = 7, ow_idx = 8, sum_idx = 9, cs_idx = 10, added_weights_idx = 11; // 0-indexed - DATA

  // char anc_file_name[NAME_LEN];
  // strcpy(anc_file_name, input_name);
  // strcat(anc_file_name, ".anc");

  // mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
  // struct stat *anc_sb = &(mmapped_anc->sb);
  // size_t anc_file_size = mmapped_anc->sb.st_size;

  mmapped_file *mmapped_op = mmap_file(open_pos_fname, O_RDWR, PROT_READ | PROT_WRITE, MAP_SHARED, 0);
  // struct stat *sb = &(mmapped_op->sb);
  size_t file_size = mmapped_op->sb.st_size;

  // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  // struct stat *sb_fasta = &(mmapped_fasta->sb);
  // size_t file_size_fasta = sb_fasta->st_size;

  size_t input_offset = 0, diag_offset = 0;
  int val_zero = 0;
  double lb, ub;

  struct multi_frag *dia;
  char cap_file_name[NAME_LEN];

  if ((istep > 0) && (iter_cond_prob == 0))
  {

#pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    for (i = 0; i < seqnum; i++)
      for (j = 0; j < seqnum; j++)
        cont_it_p[i][j] = 0;
  }

  // dia = diagonal;
  int line_num = 1;
  if ((dia = (struct multi_frag *)calloc(1, sizeof(struct multi_frag))) == NULL)
  {
    printf(" problems with memory allocation for `diags fragments' %d !  \n \n", line_num);
    exit(1);
  }

  if (anchor_step == 1)
  {
    int zero_idx = 1;
    dia = get_anchor_by_num(mmapped_fasta, mmapped_anc, &zero_idx, 0, 0);
    printf("dia->anc_num:%d\n", dia->anc_num);                                                   // DEBUG
    printf("dia->line_num:%d dia->s[0]:%d dia->s[1]:%d\n", dia->line_num, dia->s[0], dia->s[1]); // DEBUG
  }
  else
  {
    // printf("here 18\n"); // DEBUG
    diag_line line_val = get_diagvals_line_mmapped(mmapped_diags, &line_num, 0);
    diag_offset = line_val.line_offset;
    dia->b[0] = line_val.b0;
    dia->b[1] = line_val.b1;
    dia->s[0] = line_val.s0;
    dia->s[1] = line_val.s1;
    dia->hv = line_val.hv;
    dia->anc_num = line_val.anc_num;
    dia->ext = line_val.ext;
    dia->cs = line_val.cs;
    dia->fasta_offset = 0;
    dia->it = line_val.istep;
    dia->line_num = line_val.line_num;
    dia->line_offset = line_val.line_offset;
    dia->weight = line_val.weight;
    dia->ow = line_val.ow;
    dia->trans = line_val.trans;
    dia->sel = line_val.sel;
    dia->diag_count = line_val.numsubseq;
    // dia->added_weights = line_val.added_weights;
    // printf(">%d,%d,%d,%d:%d,%d,%d,%d,%u,%u,%d,%g,%g,%g,%u\n", line_val.s0, line_val.s1, hv, istep, line_val.anc_num, line_val.numsubseq, line_val.b0, line_val.b1, line_val.sel, line_val.trans, line_val.ext, line_val.weight, line_val.ow, line_val.total_sum, line_val.cs); // DEBUG
  }
  number_bf = *number;
  printf(" number_bf:%d\n", number_bf); // DEBUG

  if ((istep == 0) && anchors && (seqnum > 2))
  {
    strcpy(cap_file_name, input_name);
    strcat(cap_file_name, ".cap");
    fp_cap = fopen(cap_file_name, "w");
  }

  // printf("here 18.1\n"); // DEBUG
  // #pragma omp parallel for shared(clos) num_threads(omp_get_max_threads()) reduction(+ : tot_weight)
  for (nv = 0; nv < number_bf; nv++)
  {
    // printf("here 18.1.0\n"); // DEBUG

    // line_num++;
    // printf(" line_num:%d dia:%p\n", line_num, dia);               // DEBUG
    // printf(" dia->s[0]:%d dia->s[1]:%d\n", dia->s[0], dia->s[1]); // DEBUG
    ab[0] = dia->b[0];        /* begin of n-th diagonal in 1. sequence */
    ab[1] = dia->b[1];        /* begin of n-th diagonal in 2. sequence */
    as[0] = dia->s[0];        /* 1. sequence of n-th diagonal */
    as[1] = dia->s[1];        /* 2. sequence of n-th diagonal */
    aext = dia->ext;          /* length of n-th diagonal */
    awgt = dia->weight;       /* length of n-th diagonal */
    ae[0] = ab[0] + aext - 1; /* end of n-th diagonal in 1. sequence */
    ae[1] = ab[1] + aext - 1; /* end of n-th diagonal in 2. sequence */

    printf(" Filtering sequences anchor s1:%d s2:%d\n", as[0] + 1, as[1] + 1); // DEBUG

    if (print_status)
      if (((nv + 1) % 10) == 0)
      {
        fp_st = fopen(pst_name, "w");

        fprintf(fp_st, "\n\n\n    Status of the program run:\n");
        fprintf(fp_st, "    ==========================\n\n");
        fprintf(fp_st, "      %s \n\n", input_line);
        fprintf(fp_st, "      iteration step %d \n", istep);
        fprintf(fp_st, "      checking diagonal %d for ", nv + 1);
        fprintf(fp_st, "consistency\n\n      total number of");
        fprintf(fp_st, " diagonals = %d \n\n\n\n", number_bf);
        fclose(fp_st);
      }

    test = alignableSegments(clos, as[0], ab[0], as[1], ab[1], aext);

    if (test) /* i.e current diagonal consistent with the diagonals
                 already included into the alignment */
    {
      // #pragma omp critical
      // {
      addAlignedSegments(clos, as[0], ab[0], as[1], ab[1], aext);
      // }

      if (istep)
        for (hv = 0; hv < aext; hv++)
          for (i = 0; i < 2; i++)
          {
            j = (i + 1) % 2;
            //  open_pos[ as[i] ][ as[j] ][ ab[i]+hv ] = 0;
            int ab_hv = ab[i] + hv;
            // #pragma omp critical
            //             {
            // printf("as[i]:%d as[j]:%d ab_hv:%d\n", as[i], as[j], ab_hv); // DEBUG
            openpos_value set_pos_val = set_openpos_mmapped(mmapped_op, &as[i], &as[j], &ab_hv, &val_zero, input_offset);
        // #pragma omp atomic write
        // diag_value set_pos_val = set_diags_mmapped(all_diags, &as[i], &as[j], &ab_hv, &istep, &op_idx, &val_zero, input_offset);
#pragma omp atomic write
            input_offset = set_pos_val.line_offset;
            // }
            // exit(0); // DEBUG
          }

      if (anchor_step)
        dia->sel = 1;
      else
      {
        // exit(0); //DEBUG
        dia->sel = 1;
        double selected = 1.0; // double for consistency because we are writing doubles to file
        // diag_value seq_diag = get_diags_seq_mmapped(all_diags, &(as[i]), &(as[j]), &sel_idx, input_offset);
        printf(" SELECTED s0:%d s1:%d hv:%d istep:%d\n", dia->s[0], dia->s[1], dia->hv, istep); // DEBUG
        // set_diags_mmapped(all_diags, &(dia->s[0]), &(dia->s[1]), &(dia->hv), &istep, &sel_idx, &selected, dia->line_offset);
        set_diags_mmapped(mmapped_diags, &(dia->s[0]), &(dia->s[1]), &(dia->hv), &istep, &sel_idx, &selected, dia->line_offset);
        // exit(0); // DEBUG
      }
      double added_weights = get_ow_mmapped(mmapped_ow, &(dia->s[0]), &(dia->s[1]), &val_zero, dia->line_offset).value + dia->weight;
      // set_diags_mmapped(all_diags, &(dia->s[0]), &(dia->s[1]), &(dia->hv), &istep, &added_weights_idx, &added_weights, dia->line_offset);
      printf("DEBUG dia->s[0]:%d dia->s[1]:%d added_weights:%g\n", dia->s[0], dia->s[1], added_weights); // DEBUG
      set_ow_mmapped(mmapped_ow, &(dia->s[0]), &(dia->s[1]), &val_zero, &added_weights, dia->ow_offset);
      // glob_sim[as[0]][as[1]] =
      //     glob_sim[as[0]][as[1]] + dia->weight;

      if (istep)
      {
#pragma omp atomic update
        tot_weight += dia->weight;
      }
    } /* if test, i.e. current diagonal consistent */
    else /* no consistency */
    {

      // #pragma omp critical
      // {
#pragma omp atomic update
      (*number)--;
      // }
      if (anchor_step)
        dia->sel = 0;
      else
      {
        dia->sel = 0;
        double not_selected = 0.0; // double for consistency because we are writing doubles to file
        // diag_value seq_diag = get_diags_seq_mmapped(all_diags, &(as[i]), &(as[j]), &sel_idx, input_offset);
        printf(" NOT SELECTED s0:%d s1:%d hv:%d istep:%d\n", dia->s[0], dia->s[1], dia->hv, istep); // DEBUG
        // set_diags_mmapped(all_diags, &(dia->s[0]), &(dia->s[1]), &(dia->hv), &istep, &sel_idx, &not_selected, input_offset);
        set_diags_mmapped(mmapped_diags, &(dia->s[0]), &(dia->s[1]), &(dia->hv), &istep, &sel_idx, &not_selected, dia->line_offset);
        // exit(0); // DEBUG
      }

#pragma omp atomic write
      cont_it_p[as[0]][as[1]] = 1;
    }

    if ((istep == 0) && anchors && (seqnum > 2))
    {
      fprintf(fp_cap, " anchor %d %d %d %d %d %f ", as[0] + 1, as[1] + 1, ab[0], ab[1], aext, awgt);
      if (dia->sel == 0)
        fprintf(fp_cap, " inconsistent ");
      fprintf(fp_cap, "\n");
    }

    // #pragma omp critical
    // dia = dia->next;

    // int next_dia = dia->anc_num + 1;
    int next_dia = dia->line_num + 1;
    if (anchor_step)
    {
      if (next_dia <= total_anc_num)
      {
#pragma omp atomic write
        dia = get_anchor_by_num(mmapped_fasta, mmapped_anc, &next_dia, dia->line_offset, 0);
        // dia = dia->next;
      }
    }
    else
    {
      diag_line line_val;
#pragma omp atomic write
      line_val = get_diagvals_line_mmapped(mmapped_diags, &next_dia, diag_offset);
      diag_offset = line_val.line_offset;
      // dia = (struct multi_frag *)calloc(1, sizeof(struct multi_frag) * 1);
      dia = (struct multi_frag *)realloc(dia, sizeof(struct multi_frag) * 1);
      dia->b[0] = line_val.b0;
      dia->b[1] = line_val.b1;
      dia->s[0] = line_val.s0;
      dia->s[1] = line_val.s1;
      dia->anc_num = line_val.anc_num;
      dia->ext = line_val.ext;
      dia->cs = line_val.cs;
      dia->fasta_offset = 0;
      dia->it = line_val.istep;
      dia->line_num = line_val.line_num;
      dia->line_offset = line_val.line_offset;
      dia->weight = line_val.weight;
      dia->ow = line_val.ow;
      dia->trans = line_val.trans;
      dia->sel = line_val.sel;
    }
    printf(" going to filter next s1:%d s2:%d\n", dia->s[0] + 1, dia->s[1] + 1); // DEBUG
  } /*  for(hv = 0 ; hv < number_bf ; hv++ )  */

  if ((istep == 0) && anchors && (seqnum > 2))
    fclose(fp_cap);

  // if (anchor_step == 0)
  // {
  //   msync(all_diags->mapped_file, all_diags->sb.st_size, MS_SYNC);
  //   munmap(all_diags->mapped_file, all_diags->sb.st_size);
  //   free(all_diags->file_name);
  //   free(all_diags);
  // }

  msync(mmapped_op->mapped_file, file_size, MS_SYNC);
  munmap(mmapped_op->mapped_file, file_size);
  free(mmapped_op->file_name);
  free(mmapped_op);

  // msync(mmapped_anc->mapped_file, anc_file_size, MS_SYNC);
  msync(mmapped_anc->mapped_file, mmapped_anc->sb.st_size, MS_SYNC);
  // munmap(mmapped_anc->mapped_file, anc_file_size);
  // free(mmapped_anc);

} /*   filter(  )  */

void sel_test(mmapped_file *mapped_diags)
{
  int hv;
  // struct multi_frag *hp;

  // hp = this_it_dia;
  int line_num = 1;
  int total_diagonals = get_seqcount_mmapped_fasta(mapped_diags);
  diag_line diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, 0);
  size_t diag_offset = diag_data.line_offset;

  // for (hv = 0; hv < num_dia_af[istep]; hv++)
  for (hv = 0; hv < total_diagonals; hv++)
  {
    diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, diag_offset);
    diag_offset = diag_data.line_offset;
    if (diag_data.sel == 0)
    {
      printf("\n \n \n   WARNING:   \n \n \n");
      printf(" sel[%d] = %d \n", hv, diag_data.sel);
      exit(2);
    }
    // hp = hp->next;
    line_num++;
  }
}

void throw_out(double *weight_sum, mmapped_file *mapped_diags)
{
  int nc;
  short consist_found = 0;

  // struct multi_frag *cp; /* current diagonal */
  // struct multi_frag *hp; /* predecedor of cp */

  // // hp = (struct multi_frag *)calloc(1, sizeof(struct multi_frag));
  // hp = (struct multi_frag *)malloc(1 * sizeof(struct multi_frag));
  // cp = this_it_dia;
  // hp = NULL;

  *weight_sum = 0;

  int line_num = 1;
  int total_diagonals = get_seqcount_mmapped_fasta(mapped_diags);
  diag_line diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, 0);
  size_t diag_offset = diag_data.line_offset;
  // for (nc = 0; nc < num_dia_bf[istep]; nc++)
  for (nc = 0; nc < total_diagonals; nc++)
  {
    diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, diag_offset);
    diag_offset = diag_data.line_offset;
    if (diag_data.sel)
    {
      *weight_sum = *weight_sum + diag_data.weight;
      consist_found = 1;

      // // hp = cp;
      // // cp = cp->next;
      // line_num++;
    }

    line_num++;
    //     else
    //     {
    //       // cp = cp->next;
    //       if (consist_found)
    //       {
    // #pragma omp critical
    //         {
    //           free(hp->next);
    //         }
    //         hp->next = cp;
    //       }
    //       else
    //       {
    // #pragma omp critical
    //         {
    //           free(this_it_dia);
    //         }
    //         this_it_dia = cp;
    //       }
    //     }
  }
} /* throw_out */

size_t new_shift_mmap(int s, int p, int dif, mmapped_file *mapped_fasta, size_t input_offset)
/* shifts the elements of sequence s starting with  position p
   for dif elements to the right */
{
  int hv;
  int shift_dif; /* length of a gap (if existing) between position hv
                    and position hv+1. In case of gaps, the function
                    `new_shift' diminishs the lengths of the gaps instead
                    of shifting further sequence elements to the right  */

  fasta_len_value *len_val = get_seqlen_mmapped_fasta(mapped_fasta, &s, input_offset);
  size_t seqlen_s = len_val->len;
  size_t ret_offset = len_val->line_offset;

  // for (hv = p; (hv < seqlen[s] + 1) && (dif > 0); hv++)
  for (hv = p; (hv < seqlen_s + 1) && (dif > 0); hv++)
  {
    shift_dif = shift[s][hv + 1] - shift[s][hv] - 1;
    shift[s][hv] = shift[s][hv] + dif;
    dif = dif - shift_dif;
  }

  free(len_val);
  return ret_offset;
}

wgt_type_count(int num, int e_len, int *plus_cnt, int *minus_cnt,
               int *nuc_cnt, int *frg_inv, struct multi_frag *dia)
{

  int i, dc, pc, s1, pos;

  for (dc = 0; dc < num; dc++)
  {

    for (pc = 0; pc < dia->ext; pc++)
    {
      i = dia->b[0] + pc;
      s1 = dia->s[0];
      pos = shift[s1][i];
      if (dia->trans)
        if (dia->cs)
          minus_cnt[pos] = minus_cnt[pos] + 1;
        else
          plus_cnt[pos] = plus_cnt[pos] + 1;
      else
      {
        nuc_cnt[pos] = nuc_cnt[pos] + 1;
      }
      frg_inv[pos] = frg_inv[pos] + 1;
    }
    dia = dia->next;
  }
}

plot_calc(int num, int e_len, double *w_count, double *pl,
          struct multi_frag *dia, FILE *fp_csc)
{
  int i, dc, pc, s1, pos;
  double max_weight = 0; /* maximum value of `weight_count' */
  double shrink, shrink_csc, hsc;

  for (dc = 0; dc < num; dc++)
  {

    for (pc = 0; pc < dia->ext; pc++)
    {
      i = dia->b[0] + pc;
      s1 = dia->s[0];
      pos = shift[s1][i];
      w_count[pos] = w_count[pos] + dia->weight;
    }
    dia = dia->next;
  }

  for (i = 0; i <= e_len; i++)
    if (max_weight < w_count[i])
      max_weight = w_count[i];

  if (max_weight)
  {
    shrink = plot_num / max_weight;
    shrink_csc = MAX_CSC / max_weight;

    for (i = 0; i <= e_len; i++)
      pl[i] = w_count[i] * shrink;

    if (col_score)
    {
      printf(" e_len = %d \n\n", e_len);
      for (i = 0; i <= e_len; i++)
      {
        hsc = w_count[i] * shrink_csc;
        fprintf(fp_csc, "%5.1f\t0\n", hsc);
      }
    }
  }
  else
  {
    for (i = 0; i <= e_len; i++)
      pl[i] = 0;

    printf(" e_len = %d \n\n", e_len);
    printf(" no max weight\n\n");
  }
}

void av_tree_print_mmap(mmapped_file *mapped_fasta, mmapped_file *mapped_ow)
{
  int i, j, k, connect, max_pair[2], cv, m1, m2;
  struct subtree *all_clades;
  double **clade_similarity, new_similarity;
  double max_sim;
  char *string, l_name[2][20];
  char tree_name[NAME_LEN];
  double max_seq_sim, branch_len[2], depth;

  FILE *t_file;

  if ((all_clades = (struct subtree *)
           calloc(seqnum, sizeof(struct subtree))) == NULL)
  // if ((all_clades = (struct subtree *)
  //  malloc(seqnum * sizeof(struct subtree))) == NULL)
  {
    printf(" problems with memory allocation for `all_clades'\n \n");
    exit(1);
  }

  // if ((clade_similarity = (double **)
  //          calloc(seqnum, sizeof(double *))) == NULL)
  if ((clade_similarity = (double **)
           malloc(seqnum * sizeof(double *))) == NULL)
    exit(1);

  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
    if ((clade_similarity[i] = (double *)
             calloc(seqnum, sizeof(double))) == NULL)
      // if ((clade_similarity[i] = (double *)
      //  malloc(seqnum * sizeof(double))) == NULL)
      exit(1);

  if ((string = (char *)
           calloc(seqnum * 100, sizeof(char))) == NULL)
  // if ((string = (char *)
  //  malloc(seqnum * 100 * sizeof(char))) == NULL)
  {
    printf(" problems with memory allocation for `string'\n \n");
    exit(1);
  }

  size_t input_offset = 0;

  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
  {

    fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &i, input_offset);
    // seq_name_val->data++;
    char *seq_name_i = seq_name_val->data;

#pragma omp atomic write
    input_offset = seq_name_val->line_offset;

    // printf("seq_name_i[%d]: %s\n", i, seq_name_i); // DEBUG

    if ((all_clades[i].member = (int *)
             calloc(seqnum + 1, sizeof(int))) == NULL)
    // if ((all_clades[i].member = (int *)
    //          malloc(seqnum * sizeof(int))) == NULL)
    {
      printf(" problems with memory allocation for `all_clades'\n \n");
      exit(1);
    }

    if ((all_clades[i].name = (char *)
             calloc(strlen(seq_name_val->data) + 2, sizeof(char))) == NULL)
    // if ((all_clades[i].name = (char *)
    //  malloc(seqnum * 100 * sizeof(char))) == NULL)
    {
      printf(" problems with memory allocation for `all_clades'\n \n");
      exit(1);
    }

    // strcpy(all_clades[i].name, seq_name[i]);
    strcpy(all_clades[i].name, seq_name_i);

    all_clades[i].member_num = 1;
    all_clades[i].member[0] = i;
    all_clades[i].valid = 1;
    all_clades[i].depth = 0;

    free(seq_name_i);
    // free(seq_name_val->data);
    free(seq_name_val);
  }

  // printf("\n here 18.1.1 \n"); // DEBUG
  size_t input_offset_i = 0, input_offset_j = 0;
  int added_weights_idx = 11, zero_idx = 0;
#pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
#pragma omp parallel for shared(i) num_threads(omp_get_max_threads())
    for (j = i + 1; j < seqnum; j++)
    {
      // clade_similarity[i][j] = glob_sim[i][j];
      // clade_similarity[j][i] = glob_sim[i][j];
      // ow_value added_weights_ij = get_ow_mmapped(mapped_ow, &i, &j, &added_weights_idx, input_offset_i);
      // ow_value added_weights_ji = get_ow_mmapped(mapped_ow, &j, &i, &added_weights_idx, input_offset_j);
      ow_value added_weights_ij = get_ow_mmapped(mapped_ow, &i, &j, &zero_idx, input_offset_i);
      ow_value added_weights_ji = get_ow_mmapped(mapped_ow, &j, &i, &zero_idx, input_offset_j);
      clade_similarity[i][j] = added_weights_ij.value;
      clade_similarity[j][i] = added_weights_ji.value;
      input_offset_i = added_weights_ij.line_offset;
      input_offset_j = added_weights_ji.line_offset;
    }

  for (connect = 1; connect < seqnum; connect++)
  {
    max_sim = -1;

    for (i = 0; i < seqnum; i++)
      for (j = 0; j < seqnum; j++)
        if (i != j)
          if (all_clades[i].valid && all_clades[j].valid)
            if (clade_similarity[i][j] > max_sim)
            {
              max_sim = clade_similarity[i][j];
              max_pair[0] = i;
              max_pair[1] = j;
            }

    depth = 1 / (max_sim + 1);

    {
      m1 = max_pair[0];
      m2 = max_pair[1];

      for (i = 0; i < seqnum; i++)
        if (all_clades[i].valid)
          if (i != m1)
            if (i != m2)
            {
              if (!strcmp(clust_sim, "av"))
                new_similarity =
                    (clade_similarity[i][m1] * all_clades[m1].member_num +
                     clade_similarity[i][m2] * all_clades[m2].member_num) /
                    (all_clades[m1].member_num + all_clades[m2].member_num);

              if (!strcmp(clust_sim, "max"))
                new_similarity =
                    maxf2(clade_similarity[i][m1], clade_similarity[i][m2]);

              if (!strcmp(clust_sim, "min"))
                new_similarity =
                    minf2(clade_similarity[i][m1], clade_similarity[i][m2]);

              clade_similarity[i][m1] = new_similarity;
              clade_similarity[m1][i] = new_similarity;
            }

      all_clades[m2].valid = 0;

      for (k = 0; k < all_clades[m2].member_num; k++)
        all_clades[m1].member[all_clades[m1].member_num + k] =
            all_clades[m2].member[k];

      all_clades[m1].member_num =
          all_clades[m1].member_num + all_clades[m2].member_num;

      for (k = 0; k < 2; k++)
      {
        branch_len[k] = depth - all_clades[max_pair[k]].depth;
        sprintf(l_name[k], ":%f", branch_len[k]);
      }

      all_clades[m1].depth = depth;

      strcpy(string, "(");
      strcat(string, all_clades[m1].name);
      strcat(string, l_name[0]);
      /*            strcat(string,",\n");   */
      strcat(string, all_clades[m2].name);
      strcat(string, l_name[1]);
      strcat(string, ")");

      strcpy(all_clades[m1].name, string);
    }
  }

  strcat(string, ";");

  i = strlen(string) + 2;
  printf("i:%d\n", i); // DEBUG
  if ((upg_str = (char *)calloc(i, sizeof(char))) == NULL)
  // if ((upg_str = (char *)malloc(i * sizeof(char))) == NULL)
  {
    printf(" problems with memory allocation for `upg_str'\n \n");
    exit(1);
  }

  // printf("\n here 18.1.2 \n"); // DEBUG

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i <= strlen(string); i++)
    upg_str[i] = string[i];

  // free(all_clades->name);
  // free(all_clades->member);
  // free(all_clades);

  // free(string);
  // free(clade_similarity);

} /*  av_tree_print  */

// void print_log(struct multi_frag *d, FILE *fp_l, FILE *fp_fs)
void print_log(FILE *fp_l, FILE *fp_fs, mmapped_file *mapped_fasta, mmapped_file *mapped_diags)
{
  int i, j, mn, pv, percent, this_frag_trans, frg_count = 0;
  // struct multi_frag *diagonal;
  char hc;

  if (long_output)
  {
    fprintf(fp_l, " \n \n  Iteration %d:\n", istep);

    if (istep < 10)
      fprintf(fp_l, "  ------------");
    else
      fprintf(fp_l, "  -------------");
  }

  size_t seqname_offset_i = 0, seqname_offset_j = 0, diag_offset = 0; // d->line_offset;
  long line_num = 1;

  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
  {
    fasta_value *seq_name_i_val = get_seqname_mmapped_fasta(mapped_fasta, &i, seqname_offset_i);
    // seq_name_i_val->data++;
    char *seq_name_i = seq_name_i_val->data;
#pragma omp atomic write
    seqname_offset_i = seq_name_i_val->line_offset + strlen(seq_name_i) + 1;
    fasta_len_value *len_val_i = get_seqlen_mmapped_fasta(mapped_fasta, &i, seqname_offset_i);
    size_t seq_len_i = len_val_i->len;
    fasta_value *seq_i_val = get_seq_mmapped_fasta(mapped_fasta, &i, seqname_offset_i);
    // seq_i_val->data++;
    char *seq_i = seq_i_val->data;
#pragma omp atomic write
    seqname_offset_j = seq_i_val->line_offset + seq_len_i + 1;
    // #pragma omp parallel for shared(i) num_threads(omp_get_max_threads())
    for (j = i + 1; j < seqnum; j++)
    {

      fasta_value *seq_name_j_val = get_seqname_mmapped_fasta(mapped_fasta, &j, seqname_offset_j);
      // seq_name_j_val->data++;
      char *seq_name_j = seq_name_j_val->data;
#pragma omp atomic write
      seqname_offset_j = seq_name_j_val->line_offset + strlen(seq_name_j) + 1;
      fasta_value *seq_j_val = get_seq_mmapped_fasta(mapped_fasta, &j, seqname_offset_j);
      // seq_j_val->data++;
      char *seq_j = seq_j_val->data;
      fasta_len_value *len_val_j = get_seqlen_mmapped_fasta(mapped_fasta, &j, seqname_offset_j);
      size_t seq_len_j = len_val_j->len;
#pragma omp atomic write
      seqname_offset_j = seq_j_val->line_offset + seq_len_j + 1;
      if (long_output)
      {

        if (seqnum > 2)
        {
          fprintf(fp_l, "\n \n \n \n  Pairwise alignment ");
          fprintf(fp_l, "%d/%d", i + 1, j + 1);
          // fprintf(fp_l, " (%s / %s) \n", seq_name[i], seq_name[j]);
          fprintf(fp_l, " (%s / %s) \n", seq_name_i, seq_name_j);
          fprintf(fp_l, "  =========================");
          fprintf(fp_l, "===================== ");
        }
        fprintf(fp_l, " \n \n \n");
      }

      pairalignsum = 0;
      pairalignlen = 0;

      // diagonal = d;
      diag_line diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, diag_offset);
      diag_offset = diag_data.line_offset;
      // diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, diag_offset);
      for (int diag_idx = 0; diag_idx < diag_data.numsubseq; diag_idx++)
      {
        // }
        // while (diagonal != NULL)
        // {
        frg_count++;
        // if (diagonal->s[0] == i && diagonal->s[1] == j)
        if (diag_data.s0 == i && diag_data.s1 == j)
        {
          // if (diagonal->sel)
          if (diag_data.sel)
          {
            if (long_output)
            {
              fprintf(fp_l, "   *");
              // fprintf(fp_l, " (%3d,", diagonal->b[0]);
              fprintf(fp_l, " (%3d,", diag_data.b0);
            }

            // pairalignsum = pairalignsum + diagonal->weight;
            // pairalignlen = pairalignlen + diagonal->ext;
            pairalignsum = pairalignsum + diag_data.weight;
            pairalignlen = pairalignlen + diag_data.ext;
          }
          else if (long_output)
            // fprintf(fp_l, "     (%3d,", diagonal->b[0]);
            fprintf(fp_l, "     (%3d,", diag_data.b0);

          if (long_output)
          {
            // fprintf(fp_l, "%3d)  ", diagonal->b[1]);
            // fprintf(fp_l, " wgt:%7.3f ", diagonal->weight);
            fprintf(fp_l, "%3d)  ", diag_data.b1);
            fprintf(fp_l, " wgt:%7.3f ", diag_data.weight);
            if (seqnum > 2)
              if (overlap_weights)
                // fprintf(fp_l, " olw:%7.3f ", diagonal->ow);
                fprintf(fp_l, " olw:%7.3f ", diag_data.ow);
            // fprintf(fp_l, "len: %2d", diagonal->ext);
            fprintf(fp_l, "len: %2d", diag_data.ext);
            if ((wgt_type == 3) || crick_strand)
            {
              // if (diagonal->trans)
              if (diag_data.trans)
                fprintf(fp_l, "  P-frg");
              else
                fprintf(fp_l, "  N-frg");
            }

            // if (diagonal->trans)
            if (diag_data.trans)
              if (crick_strand)
              {
                // if (diagonal->cs)
                if (diag_data.cs)
                  fprintf(fp_l, ", CRICK strand ");
                else
                  fprintf(fp_l, ", WATSON strand ");
              }
          }

          if (frg_mult_file_v)
          {
            fprintf(fp_fs, "FRG %d ", frg_count);
            // fprintf(fp_fs, "name: %s %s ", seq_name[i], seq_name[j]);
            fprintf(fp_fs, "name: %s %s ", seq_name_i, seq_name_j);

            fprintf(fp_fs, "seq: %d %d ", i + 1, j + 1);
            // fprintf(fp_fs, "beg: %d %d ", diagonal->b[0], diagonal->b[1]);
            // fprintf(fp_fs, "len: %d ", diagonal->ext);
            fprintf(fp_fs, "beg: %d %d ", diag_data.b0, diag_data.b1);
            fprintf(fp_fs, "len: %d ", diag_data.ext);

            // fprintf(fp_fs, "wgt:%7.3f ", diagonal->weight);
            fprintf(fp_fs, "wgt:%7.3f ", diag_data.weight);
            // if (diagonal->sel)
            if (diag_data.sel)
              fprintf(fp_fs, " CONS  ");
            else
              fprintf(fp_fs, " NON-CONS ");
            fprintf(fp_fs, "\n");
            fprintf(fp_fs, "SEG1   ");

            // for (pv = 0; pv < diagonal->ext; pv++)
            //   fprintf(fp_fs, "%c", seq_i[diagonal->b[0] + pv]);
            for (pv = 0; pv < diag_data.ext; pv++)
              fprintf(fp_fs, "%c", seq_i[diag_data.b0 + pv - 1]);
            fprintf(fp_fs, "\n");
            // free(seq_i);
            // // free(seq_i_val->data);
            // free(seq_i_val);

            fprintf(fp_fs, "SEG2   ");

            // for (pv = 0; pv < diagonal->ext; pv++)
            //   fprintf(fp_fs, "%c", seq_j[diagonal->b[1] + pv]);
            for (pv = 0; pv < diag_data.ext; pv++)
              fprintf(fp_fs, "%c", seq_j[diag_data.b1 + pv - 1]);
            fprintf(fp_fs, "\n");
            fprintf(fp_fs, "\n");
            // free(seq_j);
            // // free(seq_j_val->data);
            // free(seq_j_val);
          }
          if (frg_mult_file & !frg_mult_file_v)
          {
            // if (diagonal->sel)
            if (diag_data.sel)
            {
              fprintf(fp_fs, " %d %d ", i + 1, j + 1);
              // fprintf(fp_fs, " %d %d ", diagonal->b[0], diagonal->b[1]);
              // fprintf(fp_fs, " %d \n", diagonal->ext);
              fprintf(fp_fs, " %d %d ", diag_data.b0, diag_data.b1);
              fprintf(fp_fs, " %d \n", diag_data.ext);
            }
          }

          if (long_output)
          {
            fprintf(fp_l, "\n");

            if (
                wgt_type == 2 ||
                // ((wgt_type == 3) && diagonal->trans))
                ((wgt_type == 3) && diag_data.trans))
              this_frag_trans = 1;
            else
              this_frag_trans = 0;

            if (this_frag_trans)
            {
              fprintf(fp_l, "\n           ");
              // for (pv = 0; pv < diagonal->ext; pv++)
              for (pv = 0; pv < diag_data.ext; pv++)
              {
                // hc = amino_acid[amino[i][diagonal->b[0] + pv - 1]];
                hc = amino_acid[amino[i][diag_data.b0 + pv - 1]];
                if (crick_strand)
                  // if (diagonal->cs)
                  if (diag_data.cs)
                    // hc = amino_acid[amino_c[i][diagonal->b[0] + pv - 1]];
                    hc = amino_acid[amino_c[i][diag_data.b0 + pv - 1]];

                if ((pv % 3) == 0)
                  fprintf(fp_l, "/");
                if ((pv % 3) == 1)
                  fprintf(fp_l, "%c", hc);
                if ((pv % 3) == 2)
                  fprintf(fp_l, "\\");
              }
            }

            fprintf(fp_l, "\n           ");
            // #pragma omp parallel for num_threads(omp_get_max_threads())
            // for (pv = 0; pv < diagonal->ext; pv++)
            // fasta_value *seq_i_val = get_seq_mmapped_fasta(mapped_fasta, &i, seqname_offset_i);
            // char *seq_i = seq_i_val->data;
            // printf("seq_i:%s\n", seq_i); // DEBUG
            for (pv = 0; pv < diag_data.ext; pv++)
            {
              // #pragma omp critical
              //               {
              // fprintf(fp_l, "%c", seq_i[diagonal->b[0] + pv]);
              fprintf(fp_l, "%c", seq_i[diag_data.b0 + pv - 1]);
              // }
            }
            fprintf(fp_l, "\n");

            // free(seq_i);
            // // free(seq_i_val->data);
            // free(seq_i_val);

            fprintf(fp_l, "           ");
            // #pragma omp parallel for num_threads(omp_get_max_threads())
            // for (pv = 0; pv < diagonal->ext; pv++)

            // seqname_offset_j = seq_j_val->line_offset + strlen(seq_j) + 1;
            // printf("seq_j:%s\n", seq_j); // DEBUG
            for (pv = 0; pv < diag_data.ext; pv++)
            {

              // #pragma omp critical
              //               {
              // fprintf(fp_l, "%c", seq_j[diagonal->b[1] + pv]);
              fprintf(fp_l, "%c", seq_j[diag_data.b1 + pv - 1]);
              // }
            }
            // free(seq_j);
            // // free(seq_j_val->data);
            // free(seq_j_val);

            if (this_frag_trans)
            {
              fprintf(fp_l, "\n           ");
              // for (pv = 0; pv < diagonal->ext; pv++)
              for (pv = 0; pv < diag_data.ext; pv++)
              {
                // hc = amino_acid[amino[j][diagonal->b[1] + pv - 1]];
                hc = amino_acid[amino[j][diag_data.b1 + pv - 1]];
                if (crick_strand)
                  // if (diagonal->cs)
                  if (diag_data.cs)
                    // hc = amino_acid[amino_c[j][diagonal->b[1] + pv - 1]];
                    hc = amino_acid[amino_c[j][diag_data.b1 + pv - 1]];

                if ((pv % 3) == 0)
                  fprintf(fp_l, "\\");
                if ((pv % 3) == 1)
                  fprintf(fp_l, "%c", hc);
                if ((pv % 3) == 2)
                  fprintf(fp_l, "/");
              }
            }

            fprintf(fp_l, "\n \n");
          }
        } /*  if( diagonal->s[0] == i && diagonal->s[1] == j)  */

        // diagonal = diagonal->next;
        diag_data = get_diagvals_mmapped(mapped_diags, &i, &j, &diag_idx, &(diag_data.istep), diag_data.line_offset);

      } /*  while(diagonal != NULL) */

      // percent = pairalignlen * 100 / mini2(seqlen[i], seqlen[j]);
      percent = pairalignlen * 100 / mini2((int)seq_len_i, (int)seq_len_j);

      if (long_output)
      {
        fprintf(fp_l, "\n      Sum of diagonal scores: %f\n", pairalignsum);
        fprintf(fp_l, "      Aligned residues: %d\n", pairalignlen);
        fprintf(fp_l, "      (%d percent of the shorter", percent); // TODO: Fix percent, fix output logs, fix some bugs in pairwise alignment (some fragments are chosen wrongly, check alignment of seq 1 and 8 in log)
        fprintf(fp_l, " sequence aligned)\n");
      }

      line_num++;

      free(seq_j);
      free(seq_j_val);
      free(seq_name_j);
      // free(seq_name_i_val->data);
      // free(seq_name_j_val->data);
      free(seq_name_j_val);
    } /*   for(j = i + 1 ; j < seqnum ; j++)  */
    free(seq_name_i);
    free(seq_name_i_val);
  } /* for(i = 0     ; i < seqnum ; i++) */

} /*  print_log  */
