

/*******************\
*                   *
*     DIALIGN 2     *
*                   *
*     output.c      *
*                   *
\*******************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "define.h"
#include "dialign.h"
#include "alig_graph_closure.h"

extern char seq_file[NAME_LEN];

extern int cd_gobics, wgt_type_plot, col_score;
extern int ref_seq, anchors, speed_optimized, online;
extern short crick_strand;
extern double sf_mat_thr;
extern int **amino;
extern char amino_acid[22];
extern int quali_num, wgt_plot, mask, lgs_option;
extern char input_line[NAME_LEN];
extern char clust_sim[NAME_LEN];
extern int msf_file, cw_file;
extern int lmax;
extern double threshold;
extern double av_len;
extern int pr_av_max_nd, wgt_type;
extern int num_dia_p, overlap_weights;
extern int fasta_file;
extern char *upg_str;
extern int plot_num;
extern char *seq[MAX_SEQNUM];
extern int *seqlen;
extern int maxlen;
extern char *seq_name[MAX_SEQNUM];
extern char *full_name[MAX_SEQNUM];
extern int **shift;

extern struct multi_frag *pair_dia;
extern struct multi_frag *this_it_dia;
extern struct multi_frag *all_it_dia;

extern CLOSURE *clos;

extern int max_sim_score;
extern int max_dia;
extern int seqnum;
extern int num_all_it_dia;
extern int frg_count;

extern int mini2(int a, int b);
extern size_t new_shift_mmap(int s, int p, int dif, mmapped_file *mapped_fasta, size_t input_offset);
extern void mini(int *a, int b);
extern void maxi(int *a, int b);
extern int int_test(float f);
extern plot_calc(int num, int e_len, double *w_count, double *pl,
                 struct multi_frag *dia, FILE *fp_csc);
extern wgt_type_count(int num, int e_len, int *plus_cnt, int *minus_cnt,
                      int *nuc_cnt, int *frg_inv, struct multi_frag *dia);

extern mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset);
extern size_t get_seqcount_mmapped_fasta(mmapped_file *mapped_fasta);
extern fasta_value *get_seq_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_value *get_seqname_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern fasta_len_value *get_seqlen_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern diag_value get_diagval_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, size_t input_offset);
extern diag_line get_diagvals_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, size_t input_offset);
extern diag_value get_diagval_line_mmapped(mmapped_file *mapped_file, int *line_num, int *val_index, size_t input_offset);
extern diag_line get_diagvals_line_mmapped(mmapped_file *mapped_file, int *line_num, size_t input_offset);
extern size_t find_in_file(char *mapped_file, struct stat *sb, char *search_pattern);
extern size_t find_index_in_fileline(char *line_map, struct stat *sb, int *index, size_t *line_offset);

char *itoa(int num, char *buffer, int base)
{
  int current = 0;
  if (num == 0)
  {
    buffer[current++] = '0';
    buffer[current] = '\0';
    return buffer;
  }
  int num_digits = 0;
  if (num < 0)
  {
    if (base == 10)
    {
      num_digits++;
      buffer[current] = '-';
      current++;
      num *= -1;
    }
    else
      return NULL;
  }
  num_digits += (int)floor(log(num) / log(base)) + 1;
  while (current < num_digits)
  {
    int base_val = (int)pow(base, num_digits - 1 - current);
    int num_val = num / base_val;
    char value = num_val + '0';
    buffer[current] = value;
    current++;
    num -= base_val * num_val;
  }
  buffer[current] = '\0';
  return buffer;
}

static double PRECISION = 0.00000000000001;
static int MAX_NUMBER_STRING_SIZE = 32;

/**
 * Double to ASCII
 */
char *dtoa(char *s, double n, int width, double precision)
{
  // handle special cases
  if (isnan(n))
  {
    strcpy(s, "nan");
  }
  else if (isinf(n))
  {
    strcpy(s, "inf");
  }
  else if (n == 0.0)
  {
    strcpy(s, "0");
  }
  else
  {
    int digit, m, m1;
    char *c = s;
    int neg = (n < 0);
    if (neg)
      n = -n;
    // calculate magnitude
    m = log10(n);
    int useExp = (m >= 14 || (neg && m >= 9) || m <= -9);
    if (neg)
      *(c++) = '-';
    // set up for scientific notation
    if (useExp)
    {
      if (m < 0)
        m -= 1.0;
      n = n / pow(10.0, m);
      m1 = m;
      m = 0;
    }
    if (m < 1.0)
    {
      m = 0;
    }
    // convert the number
    // while (n > PRECISION || m >= 0)
    while (n > precision || m >= 0)
    {
      double weight = pow(10.0, m);
      if (weight > 0 && !isinf(weight))
      {
        digit = floor(n / weight);
        n -= (digit * weight);
        *(c++) = '0' + digit;
      }
      if (m == 0 && n > 0)
        *(c++) = '.';
      m--;
    }
    if (useExp)
    {
      // convert the exponent
      int i, j;
      *(c++) = 'e';
      if (m1 > 0)
      {
        *(c++) = '+';
      }
      else
      {
        *(c++) = '-';
        m1 = -m1;
      }
      m = 0;
      while (m1 > 0)
      {
        *(c++) = '0' + m1 % 10;
        m1 /= 10;
        m++;
      }
      c -= m;
      for (i = 0, j = m - 1; i < j; i++, j--)
      {
        // swap without temporary
        c[i] ^= c[j];
        c[j] ^= c[i];
        c[i] ^= c[j];
      }
      c += m;
    }
    *(c) = '\0';
  }
  return s;
}

int get_substmat_mmapped(mmapped_file *mapped_file, int *i, int *j, int *k, int *l)
{
  char search_pattern[128];
  sprintf(search_pattern, ">%d,%d,%d:", *i, *j, *k); // Construct search pattern

  size_t file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern) + 1; // Move past \n
  char *line_map = &mapped_file[file_offset];
  // printf(" line_map: %s\n", line_map);                                            // DEBUG
  size_t line_offset = find_index_in_fileline(line_map, &(mapped_file->sb), l, &file_offset); //+ 1; // Move past ,
  if (line_offset == -1)
  {
    // printf(" search_pattern: %s\n", search_pattern); // DEBUG
    // printf(" line_offset: %d\n", line_offset);       // DEBUG
    // printf(" file_offset: %d\n", file_offset);       // DEBUG
    // char *line_maptmp = &mapped_file[file_offset];// DEBUG
    // size_t len_tmp = strcspn(line_maptmp, "\n");// DEBUG
    // char *line_tmp = strndup(line_maptmp, len_tmp);// DEBUG
    // printf(" line_tmp: %s\n", line_tmp); // DEBUG

    return -1;
  }
  // printf(" mapped_file[%d]: %c\n", line_offset, mapped_file[line_offset]); // DEBUG

  // size_t file_offset = 0;
  // int w = 0, x = 0, y = 0, z = 0;

  // while (file_offset < sb->st_size)
  // {
  //   char *line_map = &mapped_file[file_offset];

  //   // Find the length of the current line
  //   size_t len = strcspn(line_map, "\n");
  //   char *line = strndup(line_map, len);

  //   // Read i and j index from the line
  //   if (line[0] == '>')
  //   {
  //     sscanf(line, "> %d %d %d", &x, &y, &z);
  //   }
  //   else
  //   {
  //     if (x == *i && y == *j && z == *k)
  //     {
  //       // size_t file_pos = offset - len + *l - 1;
  //       int current_index = 0;
  //       size_t line_offset = file_offset;

  //       while (current_index < *l)
  //       {
  //         if (mapped_file[line_offset] == ',' && mapped_file[line_offset] != '\n')
  //         {
  //           current_index++; // Increment the index when a comma is found
  //         }
  //         line_offset++; // Move to the next character
  //       }

  int ret_value = -1;
// ret_value = line[*l] - '0';
#pragma omp critical
  {
    // ret_value = mapped_file[file_pos] - '0';
    ret_value = mapped_file->mapped_file[line_offset] - '0';
  }
  // free(line); // Free memory
  // printf("ret_value: %d\n", ret_value); // DEBUG
  return ret_value;
  //   }
  // }
  //   file_offset += len + 1;
  //   free(line); // Free the line after processing
  // }

  // return -1; // Return -1 if position not found
}

int add_substmat_mmapped(mmapped_file *mapped_file, int *i, int *j, int *k, int *l)
{

  char search_pattern[128];
  sprintf(search_pattern, ">%d,%d,%d:", *i, *j, *k); // Construct search pattern

  size_t file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern) + 1; // Move past \n
  char *line_map = &mapped_file[file_offset];
  // printf(" line_map: %s\n", line_map);                                            // DEBUG
  size_t line_offset = find_index_in_fileline(line_map, &(mapped_file->sb), l, &file_offset); // + 1; // Move past ,

  int ret_value = -1;

//   size_t file_offset = 0;
//   int x = 0, y = 0, z = 0;

//   while (file_offset < sb->st_size)
//   {
//     char *line_map = NULL;
// #pragma omp critical
//     {
//       line_map = &mapped_file[file_offset];
//     }

//     if (line_map == NULL)
//     {
//       break;
//     }

//     // Find the length of the current line
//     size_t len = strcspn(line_map, "\n");
//     char *line = strndup(line_map, len);

//     if (line[0] == '>')
//     {
//       sscanf(line, "> %d %d %d", &x, &y, &z);
//     }
//     else
//     {
//       if (x == *i && y == *j && z == *k)
//       {
//         // Update the z-th position with the new value
//         // memcpy(&mapped_file[offset - len + z], &value_char, 1);
//         size_t line_offset = file_offset;
//         int current_index = 0;
//         while (current_index < *l)
//         {
//           if (mapped_file[line_offset] == ',' && mapped_file[line_offset] != '\n')
//           {
//             current_index++; // Increment the index when a comma is found
//           }
//           line_offset++; // Move to the next character
//         }
#pragma omp critical
  {
    // ret_value = line[*l] - '0';
    char buffer[256];
    ret_value = mapped_file->mapped_file[line_offset] - '0';
    ret_value++;
    // char value_char = ret_value + '0';
    // char value_char = (char)ret_value;
    // char *value_char = (char *)malloc(sizeof(int) * sizeof(char));
    char *value_char = NULL;
    value_char = itoa(ret_value, buffer, 10);
    // sprintf(value_char, "%d", ret_value);

    // size_t file_pos = offset - len + *l - 1;

    // line[*l]++; // Increment the value in the `line`

    // // printf("offset: %d\n", offset); // DEBUG
    // // printf("len: %d\n", len);       // DEBUG

    // printf("*i:%d *j:%d *k:%d *l:%d\n", *i, *j, *k, *l); // DEBUG
    // printf("line: %s\n", line);                          // DEBUG
    // printf("value: %d\n", ret_value);                    // DEBUG
    // printf("value_char: %c\n", *value_char);             // DEBUG
    // printf("current_index: %d\n", current_index);        // DEBUG
    // printf("*l: %d\n", *l);                              // DEBUG
    // printf("value_char: %c\n", line[*l]); // DEBUG
    // printf("*l: %d\n", *l);               // DEBUG
    // // // // printf("mapped_line_addr: %d\n", &mapped_file[offset + *l]);           // DEBUG
    // // // // printf("mapped_line: %c\n", mapped_file[offset + *l]);                 // DEBUG
    // // printf("mapped_line_addr CPT: %p\n", &mapped_file[offset - len + *l]); // DEBUG
    // // printf("mapped_line CPT: %c\n", mapped_file[offset - len + *l]);       // DEBUG
    // printf("mapped_line_addr CPT: %p\n", &mapped_file[file_pos]); // DEBUG
    // printf("(before) mapped_line[%d] CPT: %c\n", line_offset, mapped_file[line_offset]); // DEBUG

    // // mapped_file[offset - len + *l] = line[*l]; // Update mapped file

    // // // Synchronize the change to the file (optional)
    // // msync(&mapped_file[offset - len + *l], sizeof(char), MS_SYNC);

    // // // memmove(&mapped_file[offset - len + *l], &line[*l], 1); // ChatGPT
    // // // // memmove(&mapped_file[offset + *z], &value_char, 1);

    // // // // Synchronize the change to the file
    // // // // msync(mapped_file, sb->st_size, MS_SYNC);
    // // // msync(&mapped_file[offset - len + *l], sizeof(line[*l]), MS_SYNC);

    // mapped_file[file_pos] = ++line[*l]; // Update mapped file
    // mapped_file[offset] = ret_value + '0'; // Update mapped file
    mapped_file->mapped_file[line_offset] = *value_char; // Update mapped file
    int len_value_char = strlen(value_char);             // sizeof(value_char) / sizeof(value_char[0]);
    // memmove(&mapped_file[file_pos], &value_char, 1); // ChatGPT
    // printf("sizeof value_char: %d\n", sizeof(char) * len_value_char); // DEBUG
    // Synchronize the change to the file (optional)
    // msync(&mapped_file[file_pos], sizeof(char), MS_SYNC);
    msync(&(mapped_file->mapped_file[line_offset]), sizeof(char) * len_value_char, MS_SYNC);
    // printf("(after) mapped_line[%d] CPT: %c\n", line_offset, mapped_file[line_offset]); // DEBUG
    // assert(mapped_file[file_pos] == line[*l]);
    // // printf("mapped_line after replace: %c\n", mapped_file[offset + *l]);           // DEBUG
    // printf("mapped_line after replace CPT: %c\n", mapped_file[offset - len + *l]); // DEBUG
    // printf("line after replace: %s\n", line); // DEBUG

    // ret_value = line[*l] - '0';
    ret_value = mapped_file->mapped_file[line_offset] - '0';
    // printf(" ret_val : %d\n", ret_value); // DEBUG
  }
  // free(line); // Free memory
  return ret_value;
  //     }
  //   }
  //   file_offset += len + 1;

  //   free(line); // Free the line after processing
  // }

  // return -1; // Return -1 if position not found
}

void subst_mat_mmap(char *file_name, int fragno, struct multi_frag *frg, mmapped_file *mapped_fasta)
{

  int s0, s1, i, j, frg_count;
  // short a0, a1;
  int a0, a1;
  // int ****sbsmt;
  struct multi_frag *frag;
  char mat_file_name[NAME_LEN];
  FILE *fp_mat;

  //   if ((sbsmt = (int ****)malloc(seqnum * sizeof(int ***))) == NULL)
  //   {
  //     printf("Problems with memory allocation for sbsmt\n");
  //     exit(1);
  //   }

  // #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic)
  //   for (i = 0; i < seqnum; i++)
  //   {
  //     if ((sbsmt[i] = (int ***)malloc(seqnum * sizeof(int **))) == NULL)
  //     {
  //       printf("Problems with memory allocation for sbsmt\n");
  //       exit(1);
  //     }
  // #pragma omp parallel for shared(i) num_threads(omp_get_max_threads()) schedule(dynamic)
  //     for (j = 0; j < seqnum; j++)
  //     {
  //       if ((sbsmt[i][j] = (int **)malloc(21 * sizeof(int *))) == NULL)
  //       {
  //         printf("Problems with memory allocation for sbsmt\n");
  //         exit(1);
  //       }
  // #pragma omp parallel for shared(i, j) num_threads(omp_get_max_threads()) schedule(dynamic)
  //       for (a0 = 0; a0 < 21; a0++)
  //       {
  //         if ((sbsmt[i][j][a0] = (int *)malloc(21 * sizeof(int))) == NULL)
  //         {
  //           printf("Problems with memory allocation for sbsmt\n");
  //           exit(1);
  //         }
  // #pragma omp parallel for shared(i, j, a0) num_threads(omp_get_max_threads()) schedule(dynamic)
  //         for (a1 = 0; a1 <= 20; a1++)
  //         {
  //           sbsmt[i][j][a0][a1] = 0;
  //         }
  //       }
  //     }
  //   }

  // if ((sbsmt = (int ****)calloc(seqnum, sizeof(int ***))) == NULL)
  // {
  //   printf("Problems with memory allocation for sbsmt\n");
  //   exit(1);
  // }

  // for (i = 0; i < seqnum; i++)
  //   if ((sbsmt[i] = (int ***)calloc(seqnum, sizeof(int **))) == NULL)
  //   {
  //     printf("Problems with memory allocation for sbsmt\n");
  //     exit(1);
  //   }

  // for (i = 0; i < seqnum; i++)
  //   for (j = 0; j < seqnum; j++)
  //     if ((sbsmt[i][j] = (int **)calloc(21, sizeof(int *))) == NULL)
  //     {
  //       printf("Problems with memory allocation for sbsmt\n");
  //       exit(1);
  //     }

  // for (i = 0; i < seqnum; i++)
  //   for (j = 0; j < seqnum; j++)
  //     for (a0 = 0; a0 < 21; a0++)
  //       if ((sbsmt[i][j][a0] = (int *)calloc(21, sizeof(int))) == NULL)
  //       {
  //         printf("Problems with memory allocation for sbsmt\n");
  //         exit(1);
  //       }

  // for (i = 0; i < seqnum; i++)
  //   for (j = 0; j < seqnum; j++)
  //     for (a0 = 0; a0 <= 20; a0++)
  //       for (a1 = 0; a1 <= 20; a1++)
  //         sbsmt[i][j][a0][a1] = 0;

  // printf("\n here 18.2.-2 \n"); // DEBUG

  FILE *sm_f;
  int seq_file_length = sizeof(seq_file) / sizeof(char);
  char sm_fname[seq_file_length]; // = seq_file;
  // printf("\n seq_file: %s \n", seq_file); // DEBUG
  snprintf(sm_fname, sizeof(seq_file), "%s", seq_file);
  strcat(sm_fname, ".sm");
  // printf("\n here 18.2.-1 \n"); // DEBUG
  sm_f = fopen(sm_fname, "w");
  // for (i = 0; i < seqnum; i++)
  // {
  //   for (j = 0; j < seqnum; j++)
  //   {
  //     for (a0 = 0; a0 <= 20; a0++)
  //     {
  //       fprintf(sm_f, ">%d,%d,%d:", i, j, a0);
  //       for (a1 = 0; a1 < 20; a1++)
  //       {
  //         fprintf(sm_f, "0,");
  //       }
  //       fprintf(sm_f, "0");
  //       fprintf(sm_f, "\n");
  //     }
  //   }
  // }

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
  {
#pragma omp parallel for shared(i) num_threads(omp_get_max_threads())
    for (j = 0; j < seqnum; j++)
    {
#pragma omp parallel for shared(i, j) num_threads(omp_get_max_threads())
      for (a0 = 0; a0 <= 20; a0++)
      {
        char *file_line = (char *)malloc(((8 + 1) + (21 * 2)) * sizeof(char));
        strcpy(file_line, "");

        for (a1 = 0; a1 < 20; a1++)
        {
          strcat(file_line, "0,");
        }
        strcat(file_line, "0");
#pragma omp critical
        fprintf(sm_f, ">%d,%d,%d:%s\n", i, j, a0, file_line);
        free(file_line);
      }
    }
  }

  fclose(sm_f);

  // int sm_fd = open(sm_fname, O_RDWR);
  // if (sm_fd == -1)
  // {
  //   perror("Error opening file");
  //   exit(1);
  // }

  // struct stat sm_sb;
  // if (fstat(sm_fd, &sm_sb) == -1)
  // {
  //   perror("Error getting file size");
  //   exit(1);
  // }

  // // mmap the file into memory
  // char *sm_mmap = mmap(NULL, sm_sb.st_size, PROT_WRITE, MAP_SHARED, sm_fd, 0);
  // if (sm_mmap == MAP_FAILED)
  // {
  //   perror("Error mmapping file");
  //   exit(1);
  // }

  // close(sm_fd);

  mmapped_file *mmapped_sm = mmap_file(sm_fname, O_RDWR, PROT_WRITE, MAP_SHARED, 0);
  // struct stat *sm_sb = &(mmapped_sm->sb);
  size_t file_size = mmapped_sm->sb.st_size;

  // int val_zero = 0;
  // int val_one = 1;
  // int val_four = 4;
  // int val_five = 5;

  // add_substmat_mmapped(mmapped_sm->mapped_file, sm_sb, &val_zero, &val_zero, &val_zero, &val_five);
  // add_substmat_mmapped(mmapped_sm->mapped_file, sm_sb, &val_zero, &val_zero, &val_one, &val_four);

  // printf(" get_subst[0,0,0,5]: %d\n", get_substmat_mmapped(mmapped_sm->mapped_file, sm_sb, &val_zero, &val_zero, &val_zero, &val_five)); // DEBUG must be == 1
  // printf(" get_subst[0,0,1,4]: %d\n", get_substmat_mmapped(mmapped_sm->mapped_file, sm_sb, &val_zero, &val_zero, &val_one, &val_four));  // DEBUG must be == 1
  // printf(" get_subst[0,0,4,1]: %d\n", get_substmat_mmapped(mmapped_sm->mapped_file, sm_sb, &val_zero, &val_zero, &val_four, &val_one));  // DEBUG must be == 0

  // msync(mmapped_sm->mapped_file, file_size, MS_SYNC);
  // munmap(mmapped_sm->mapped_file, file_size);
  // exit(0);

  // printf("\n here 18.2.0 \n"); // DEBUG

  strcpy(mat_file_name, file_name);
  strcat(mat_file_name, ".mat");

  fp_mat = fopen(mat_file_name, "w");

  frag = frg;

  // #pragma omp parallel for shared(frag, sf_mat_thr) num_threads(omp_get_max_threads()) schedule(dynamic)
  for (frg_count = 0; frg_count < fragno; frg_count++)
  {
    if (frag->weight > sf_mat_thr)
    {
      // #pragma omp parallel for shared(frag) num_threads(omp_get_max_threads()) schedule(dynamic)
      for (i = 0; i < frag->ext; i++)
      {
        // printf(" i: %d\n", i);                                             // DEBUG
        // printf(" frag->s[0]: %d\n", frag->s[0]);                           // DEBUG
        // printf(" frag->s[1]: %d\n", frag->s[1]);                           // DEBUG
        // printf(" frag->b[0]: %d\n", frag->b[0]);                           // DEBUG
        // printf(" frag->b[1]: %d\n", frag->b[1]);                           // DEBUG
        // printf(" frag->b[0] + i: %d\n", frag->b[0] + i);                   // DEBUG
        // printf(" frag->b[1] + i: %d\n", frag->b[1] + i);                   // DEBUG
        // printf(" amino[0][0]: %d\n", amino[0][0]);                         // DEBUG
        // printf(" amino[s0][b0]: %d\n", amino[frag->s[0]][frag->b[0] + i]); // DEBUG
        // printf(" amino[s1][b1]: %d\n", amino[frag->s[1]][frag->b[1] + i]); // DEBUG

        a0 = amino[frag->s[0]][frag->b[0] + i];
        a1 = amino[frag->s[1]][frag->b[1] + i];
        s0 = frag->s[0];
        s1 = frag->s[1];
        // #pragma omp critical
        //         {
        // sbsmt[s0][s1][a0][a1]++;
        // sbsmt[s1][s0][a1][a0]++;
        add_substmat_mmapped(mmapped_sm, &s0, &s1, &a0, &a1);
        add_substmat_mmapped(mmapped_sm, &s1, &s0, &a1, &a0);
        // }
      }
    }
    // #pragma omp critical
    // {
    frag = frag->next;
    // }
  }

  // printf("\n here 18.2.1 \n"); // DEBUG

  fprintf(fp_mat, "taxanumber: %d ;\n", seqnum);
  fprintf(fp_mat, "description: DIALIGN alignment ;\n");
  fprintf(fp_mat, "description: %s;\n", input_line);

  size_t input_offset = 0;

#pragma omp parallel for num_threads(omp_get_max_threads()) schedule(dynamic)
  for (i = 0; i < seqnum; i++)
  {
    fasta_value *name_val = get_seqname_mmapped_fasta(mapped_fasta, &i, input_offset);

    if (input_offset >= mapped_fasta->sb.st_size)
    {
      printf("Error: Could not find sequence name for sequence %d\n", i);
      exit(1);
    }

    char *seq_name_i = name_val->data;

#pragma omp atomic write
    input_offset = name_val->line_offset;

    // fprintf(fp_mat, "taxon: %.3d  name: %s  ;\n", i + 1, full_name[i]);
    fprintf(fp_mat, "taxon: %.3d  name: %s  ;\n", i + 1, seq_name_i);
    free(seq_name_i);
    // free(name_val->data);
    free(name_val);
  }

  for (s0 = 0; s0 < seqnum; s0++)
    for (s1 = s0 + 1; s1 < seqnum; s1++)
      for (a0 = 1; a0 <= 20; a0++)
        for (a1 = 1; a1 < 21; a1++)
        {
          fprintf(fp_mat, "pair: %.3d %.3d ", s0 + 1, s1 + 1);
          fprintf(fp_mat, " acids: %c%c ", amino_acid[a0], amino_acid[a1]);
          // fprintf(fp_mat, " number: %d ;\n", sbsmt[s0][s1][a0][a1]);
          fprintf(fp_mat, " number: %d ;\n", get_substmat_mmapped(mapped_fasta, &s0, &s1, &a0, &a1));
        }

  msync(mmapped_sm->mapped_file, file_size, MS_SYNC);
  munmap(mmapped_sm->mapped_file, file_size);
  // free(sm_sb);
  free(mmapped_sm->file_name);
  free(mmapped_sm);

} /* subst_mat */

// void print_fragments(struct multi_frag *d, FILE *fp_ff2)
void print_fragments(mmapped_file *mapped_diags, FILE *fp_ff2)
{

  // struct multi_frag *fragment;

  // fragment = d;
  int line_num = 1;
  diag_line diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, 0);
  size_t diag_offset = diag_data.line_offset;
  int diagonal_count = get_seqcount_mmapped_fasta(mapped_diags);
  // while (fragment != NULL)
  while (line_num <= diagonal_count)
  {
    // if (fragment->it)
    if (diag_data.istep)
    {
      frg_count++;
      fprintf(fp_ff2, "%6d) ", frg_count);
      // fprintf(fp_ff2, "seq: %3d %3d  ", fragment->s[0] + 1, fragment->s[1] + 1);
      // fprintf(fp_ff2, "beg: %7d %7d ", fragment->b[0], fragment->b[1]);
      // fprintf(fp_ff2, "len: %3d ", fragment->ext);
      // fprintf(fp_ff2, "wgt: %6.2f ", fragment->weight);
      // fprintf(fp_ff2, "olw: %6.2f ", fragment->ow);
      fprintf(fp_ff2, "seq: %3d %3d  ", diag_data.s0 + 1, diag_data.s1 + 1);
      fprintf(fp_ff2, "beg: %7d %7d ", diag_data.b0, diag_data.b1);
      fprintf(fp_ff2, "len: %3d ", diag_data.ext);
      fprintf(fp_ff2, "wgt: %6.2f ", diag_data.weight);
      fprintf(fp_ff2, "olw: %6.2f ", diag_data.ow);

      // fprintf(fp_ff2, "it: %d ", fragment->it);
      fprintf(fp_ff2, "it: %d ", diag_data.istep);
      // if (fragment->sel)
      if (diag_data.sel)
        fprintf(fp_ff2, "cons   ");
      else
        fprintf(fp_ff2, "incons ");

      if ((wgt_type == 3) || crick_strand)
      {
        // if (fragment->trans)
        if (diag_data.trans)
          fprintf(fp_ff2, " P-frg");
        else
          fprintf(fp_ff2, " N-frg");
        // if (fragment->trans)
        if (diag_data.trans)
          if (crick_strand)
            // if (fragment->cs)
            if (diag_data.cs)
              fprintf(fp_ff2, " -");
            else
              fprintf(fp_ff2, " +");
      }

      fprintf(fp_ff2, "\n");
    }
    // fragment = fragment->next;
    diag_data = get_diagvals_line_mmapped(mapped_diags, &line_num, diag_offset);
  }
}

void weight_print(double **wgt)
{
  int i, j, l, s;
  FILE *fp;

  fp = fopen("weight_table", "w");

  fprintf(fp, "  len1 = %d, len2 = %d\n\n", seqlen[0], seqlen[1]);
  fprintf(fp, "  \n   %s \n\n", input_line);
  for (l = 1; l <= max_dia; l++)
    for (s = 0; s <= l * max_sim_score; s++)
      fprintf(fp, " %d %d %7.8g \n", l, s, wgt[l][s]);

  fclose(fp);

} /* weight_print */

void ali_arrange_mmap(int fragno, struct multi_frag *d, FILE *fp, FILE *fp2, FILE *fp3, FILE *fp4, FILE *fp_col_score, mmapped_file *mapped_fasta)
{
  if (fragno < 0)
  {
    return (0);
  }

  int block_no, char_no;
  int shift_cond, endlen;
  int p, pn, i, j, k, l, hv, bc, lc, max_p;
  int b1, b2, s1, s2, e, dif, sv, lv, add, msf_lines;

  char sim_char;
  float weak_wgt_type_thr = WEAK_WGT_TYPE_THR;
  float strong_wgt_type_thr = STRONG_WGT_TYPE_THR;
  float frac_plus, frac_minus, frac_nuc, f_inv;

  char **endseq;
  char **hseq;
  char *clear_seq;
  double *weight_count;
  int *plus_count;
  int *minus_count;
  int *nuc_count;
  int *frg_involved;
  double *plot; /* plot[i] = sum of weights of fragments involved at
                  position i normalizet such that the maximum value */

  char gap_char = '-';
  char ambi_char = ' ';
  int *begin, *end, *b_len, *first_pos, pl_int;
  int b_size; /* size of fragments */
  struct multi_frag *fragments, *dia;
  int **inv_shift;
  int char_per_line; /* number of residues per line in output file        */
  char aligned;
  char_per_line = ((PAPER_WIDTH - 18) / 11) * 10;

  dia = d;

  if ((endseq = (char **)malloc(seqnum * sizeof(char *))) == NULL)
  {
    printf(" problems with memory allocation for `endseq' !  \n \n");
    exit(1);
  }

  if ((hseq = (char **)malloc(seqnum * sizeof(char *))) == NULL)
  {
    printf(" problems with memory allocation for `hseq' !  \n \n");
    exit(1);
  }

  if ((begin = (int *)malloc(seqnum * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `begin' !  \n \n");
    exit(1);
  }

  if ((end = (int *)malloc(seqnum * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `end' !  \n \n");
    exit(1);
  }

  if ((b_len = (int *)malloc(seqnum * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `b_len' !  \n \n");
    exit(1);
  }

  if ((first_pos = (int *)malloc(seqnum * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `first_pos' !  \n \n");
    exit(1);
  }

  if ((shift = (int **)malloc(seqnum * sizeof(int *))) == NULL)
  {
    printf("not enough memory available for `shift' !!!!\n");
    fprintf(fp, "not enough memory available for `shift' !\n");
    exit(1);
  }
  // printf("\n here 18.3.0 \n"); // DEBUG

  size_t input_offset = 0;

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
  {
    if (input_offset >= mapped_fasta->sb.st_size)
    {
      printf("Error: Could not find sequence length for sequence %d\n", hv);
      exit(1);
    }

    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mapped_fasta, &hv, input_offset);

    size_t seq_len_hv = len_val->len;
#pragma omp atomic write
    input_offset = len_val->line_offset + len_val->len + 1;

    // if ((shift[hv] = (int *)calloc((seqlen[hv] + 2), sizeof(int))) == NULL)
    if ((shift[hv] = (int *)malloc((seq_len_hv + 2) * sizeof(int))) == NULL)
    {
      printf("not enough memory available for `shift' !!!!\n");
      fprintf(fp, "not enough memory available for `shift' !\n");
      exit(1);
    }

    if (fragno >= 0)
    {
      begin[hv] = seq_len_hv;
      end[hv] = 1;
    }
  }

  if (fragno > 0)
  {
    if ((fragments = malloc(fragno * sizeof(struct multi_frag))) == NULL)
    {
      printf("not enough memory available for fragments!\n");
      fprintf(fp, "not enough memory available for fragments!\n");
      exit(1);
    }
  }

  for (hv = 1; hv <= fragno; hv++)
  {
    fragments[hv - 1] = *dia;
    dia = dia->next;
  }

  // printf("\n here 18.3.1 \n"); // DEBUG

  // #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
  for (hv = 0; hv < fragno; hv++)
    for (j = 0; j < 2; j++)
    {
      mini(&begin[fragments[hv].s[j]], fragments[hv].b[j]);
      maxi(&end[fragments[hv].s[j]], fragments[hv].b[j] + fragments[hv].ext);
    }

  // printf("\n here 18.3.2 \n"); // DEBUG

  input_offset = 0;

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
  {
    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mapped_fasta, &hv, input_offset);

    if (input_offset >= mapped_fasta->sb.st_size)
    {
      printf("Error: Could not find sequence length for sequence %d\n", hv);
      exit(1);
    }

    size_t seq_len_hv = len_val->len;

#pragma omp atomic write
    input_offset = len_val->line_offset;

    begin[hv] = 1;
    end[hv] = seq_len_hv + 1;
  }

  b_size = 0;

  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
  {
    b_len[i] = end[i] - begin[i];
#pragma omp critical
    {
      maxi(&b_size, b_len[i]);
    }
  }

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < seqnum; i++)
#pragma omp parallel for shared(i) num_threads(omp_get_max_threads())
    for (hv = 0; hv < b_len[i]; hv++)
      shift[i][begin[i] + hv] = hv;

  shift_cond = 1;

  // printf("\n here 18.3.2.-2 \n"); // DEBUG

  input_offset = 0;

  while (shift_cond)
  {
    shift_cond = 0;

    // #pragma omp parallel for collapse(2) shared(dif) num_threads(omp_get_max_threads())
    for (hv = 0; hv < fragno; hv++)
      for (j = 0; j < 2; j++)
      {
        k = (j + 1) % 2;
        s1 = fragments[hv].s[j];
        s2 = fragments[hv].s[k];
        b1 = fragments[hv].b[j];
        b2 = fragments[hv].b[k];
        e = fragments[hv].ext;

        // #pragma omp parallel for collapse(2) shared(hv, e, j, k, s1, s2, b1, b2, l, dif) num_threads(omp_get_max_threads())
        for (l = e - 1; l >= 0; l--)
        {
          dif = shift[s2][b2 + l] - shift[s1][b1 + l];
          if (dif > 0)
          {
            size_t tmp_offset = new_shift_mmap(s1, b1 + l, dif, mapped_fasta, input_offset);

#pragma omp atomic write
            input_offset = tmp_offset;
#pragma omp atomic write
            shift_cond = 1;
          }
        }
      }
  }

  endlen = 0;

  for (hv = 0; hv < seqnum; hv++)
    maxi(&endlen, shift[hv][end[hv] - 1] + 1);

  if ((inv_shift = (int **)malloc(seqnum * sizeof(int *))) == NULL)
  {
    printf("not enough memory available for `inv_shift' !!!!\n");
    fprintf(fp, "not enough memory available for `inv_shift' !\n");
    exit(1);
  }

  // printf("\n here 18.3.2.-1 \n"); // DEBUG

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
  {
    if ((endseq[hv] = malloc(endlen * sizeof(char))) == NULL)
    {
      printf(" not enough memory available for printing results!\n");
      fprintf(fp, " not enough memory available");
      fprintf(fp, " for printing results!\n");
      exit(1);
    }

    if ((inv_shift[hv] = (int *)malloc((endlen + 2) * sizeof(int))) == NULL)
    {
      printf("not enough memory available for `inv_shift' !!!!\n");
      fprintf(fp, "not enough memory available for `inv_shift' !\n");
      exit(1);
    }

    if ((hseq[hv] = malloc((maxlen + 1) * sizeof(char))) == NULL)
    {
      printf("not enough memory available for printing results! \n");
      fprintf(fp, "not enough memory available");
      fprintf(fp, " for printing results! \n");
      exit(1);
    }
  }

  if ((clear_seq = (char *)malloc((endlen + 1) * sizeof(char))) == NULL)
  {
    printf(" problems with memory allocation for `clear_seq' !  \n \n");
    exit(1);
  }

  if ((weight_count =
           (double *)malloc((endlen + 2) * sizeof(double))) == NULL)
  {
    printf(" problems with memory allocation for `weight_count' !\n \n");
    exit(1);
  }

  if ((plot = (double *)malloc((endlen + 2) * sizeof(double))) == NULL)
  {
    printf(" problems with memory allocation for `plot' ! \n \n");
    exit(1);
  }

  if ((plus_count =
           (int *)malloc((endlen + 2) * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `plus_count' !\n \n");
    exit(1);
  }

  if ((minus_count =
           (int *)malloc((endlen + 2) * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `minus_count' !\n \n");
    exit(1);
  }

  if ((nuc_count =
           (int *)malloc((endlen + 2) * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `nuc_count' !\n \n");
    exit(1);
  }

  if ((frg_involved =
           (int *)malloc((endlen + 2) * sizeof(int))) == NULL)
  {
    printf(" problems with memory allocation for `frg_involved ' !\n \n");
    exit(1);
  }

  input_offset = 0;

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
  {
    if (input_offset >= mapped_fasta->sb.st_size)
    {
      printf("Error: Could not find sequence length for sequence %d\n", hv);
      exit(1);
    }

    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mapped_fasta, &hv, input_offset);

    size_t seq_len_hv = len_val->len;

#pragma omp atomic write
    input_offset = len_val->line_offset + len_val->len + 1;

#pragma omp parallel for shared(hv) num_threads(omp_get_max_threads())
    for (p = 1; p <= seq_len_hv; p++)
      inv_shift[hv][shift[hv][p]] = p;
  }

#pragma omp parallel for collapse(2) shared(hv) num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
    for (i = 0; i < endlen; i++)
      endseq[hv][i] = gap_char;

  // printf("\n here 18.3.2.1 \n"); // DEBUG

  input_offset = 0;

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
  {
    // printf("hv: %d\n", hv);   // DEBUG
    // printf("&hv: %p\n", &hv); // DEBUG

    if (input_offset >= mapped_fasta->sb.st_size)
    {
      printf("Error: Could not find sequence for sequence %d\n", hv);
      exit(1);
    }

    fasta_value *seq_val = get_seq_mmapped_fasta(mapped_fasta, &hv, input_offset);

    char *seq_hv = seq_val->data;

#pragma omp atomic write
    input_offset = seq_val->line_offset + strlen(seq_hv) + 1;

    // printf("seq_hv: %c\n", seq_hv[0]); // DEBUG
#pragma omp parallel for shared(hv, seq_hv) num_threads(omp_get_max_threads())
    for (i = begin[hv]; i < end[hv]; i++)
    {
      hseq[hv][i] = tolower(seq_hv[i]);
    }
    free(seq_hv);
    // free(seq_val->data);
    free(seq_val);
  }

  // printf("\n here 18.3.2.2 \n"); // DEBUG

  input_offset = 0;

#pragma omp parallel for shared(input_offset) num_threads(omp_get_max_threads())
  for (hv = 0; hv < fragno; hv++)
  {
#pragma omp parallel for shared(hv, input_offset) num_threads(omp_get_max_threads())
    for (k = 0; k < 2; k++)
    {

      if (input_offset >= mapped_fasta->sb.st_size)
      {
        printf("Error: Could not find sequence for fragment %d, %d\n", hv, k);
        exit(1);
      }

      int frag_hv_s_k = fragments[hv].s[k];
      int frag_hv_b_k = fragments[hv].b[k];
      int frag_hv_ext = fragments[hv].ext;

      fasta_value *seq_val = get_seq_mmapped_fasta(mapped_fasta, &frag_hv_s_k, input_offset);

      char *seq_frag = seq_val->data;

#pragma omp atomic write
      input_offset = seq_val->line_offset + strlen(seq_frag) + 1;

// #pragma omp parallel for shared(hv, k, seq_frag, frag_hv_s_k, frag_hv_b_k) num_threads(omp_get_max_threads())
//       for (i = fragments[hv].b[k]; i < fragments[hv].b[k] + fragments[hv].ext; i++)
//       {
//         hseq[fragments[hv].s[k]][i] = seq_frag[i];
//       }
#pragma omp parallel for shared(hv, k, seq_frag, frag_hv_s_k, frag_hv_b_k, frag_hv_ext) num_threads(omp_get_max_threads())
      for (i = frag_hv_b_k; i < frag_hv_b_k + frag_hv_ext; i++)
      {
        hseq[frag_hv_s_k][i] = seq_frag[i];
      }
      free(seq_frag);
      // free(seq_val->data);
      free(seq_val);
    }
  }

  // printf("\n here 18.3.2.3 \n"); // DEBUG

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
  {
#pragma omp parallel for shared(hv) num_threads(omp_get_max_threads())
    for (i = begin[hv]; i < end[hv]; i++)
    {
      endseq[hv][shift[hv][i]] = hseq[hv][i];
    }
  }

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (i = 0; i < endlen; i++)
    clear_seq[i] = ' ';

  // printf("\n here 18.3.3 \n"); // DEBUG

  // #pragma omp parallel for num_threads(omp_get_max_threads())
  for (p = 0; p < endlen; p++)
  {
    s1 = 0;
    while (
        (endseq[s1][p] == tolower(endseq[s1][p])) && (s1 < (seqnum - 1)) /* no capital letter */
    )
      s1++;

    if (s1 < (seqnum - 1))
    {
#pragma omp parallel for shared(p, s1, clos) num_threads(omp_get_max_threads())
      for (s2 = s1 + 1; s2 < seqnum; s2++)
      {
        if (endseq[s2][p] != tolower(endseq[s2][p]))
        /* endseq[s2][p] capital letter */
        {
          aligned = alignedPositions(clos, s1, inv_shift[s1][p], s2,
                                     succFrontier(clos, s1, inv_shift[s1][p], s2));

          if (!aligned)
            /* i.e.endseq[s1][p] not aligned with end seq[s2][p]*/
            clear_seq[p] = ambi_char;
        }
      }
    }
  }

  if (mask)
  {
#pragma omp parallel for collapse(2) shared(sv) num_threads(omp_get_max_threads())
    for (sv = 0; sv < seqnum; sv++)
    {
      for (hv = 0; hv < endlen; hv++)
      {
        if (endseq[sv][hv] != gap_char)
        {
          if (endseq[sv][hv] == tolower(endseq[sv][hv]))
            endseq[sv][hv] = '*';
        }
      }
    }
  }

  if (col_score)
  {
    fprintf(fp_col_score, "# 1 %d \n", endlen);
    fprintf(fp_col_score, "# %s \n", upg_str);
  }

  plot_calc(num_all_it_dia, endlen, weight_count, plot, all_it_dia, fp_col_score);

  wgt_type_count(num_all_it_dia, endlen, plus_count, minus_count, nuc_count, frg_involved, all_it_dia);

  lc = (endlen - 1) / char_per_line;

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
    first_pos[hv] = begin[hv];

  input_offset = 0;
  // #pragma omp parallel collapse(2) for num_threads(omp_get_max_threads()) schedule(dynamic) shared(hv, k)
  for (k = 0; k <= lc; k++)
  {
    for (hv = 0; hv < seqnum; hv++)
    {
      if (input_offset >= mapped_fasta->sb.st_size)
      {
        printf("Error: Could not find sequence name for sequence %d\n", hv);
        exit(1);
      }

      fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &hv, input_offset);
      char *seq_name_hv = seq_name_val->data;

#pragma omp atomic write
      input_offset = seq_name_val->line_offset;

      // #pragma omp critical
      // {
      fprintf(fp, "%s", seq_name_hv);
      fprintf(fp, "%8d  ", first_pos[hv]);
      // }

      // #pragma omp parallel for shared(hv) num_threads(omp_get_max_threads()) reduction(+ : first_pos[hv])
      for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
      {
        if (!(i % 10))
        {
          // #pragma omp critical
          // {
          fprintf(fp, " ");
          // }
        }
        // #pragma omp critical
        // {
        fprintf(fp, "%c", endseq[hv][k * char_per_line + i]);
        // }
        if (endseq[hv][k * char_per_line + i] != gap_char)
        {
          // #pragma omp critical
          // {
          first_pos[hv]++;
          // }
        }
      }
      // #pragma omp critical
      // {
      fprintf(fp, " \n");
      // }
      free(seq_name_hv);
      // free(seq_name_val->data);
      free(seq_name_val);
    }

    // #pragma omp critical
    // {
    fprintf(fp, "         ");
    // }

    // #pragma omp parallel for shared(hv, k) num_threads(omp_get_max_threads()) schedule(dynamic)
    for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
    {
      if (!(i % 10))
      {
        // #pragma omp critical
        // {
        fprintf(fp, " ");
        // }
      }
      // #pragma omp critical
      // {
      fprintf(fp, "%c", clear_seq[k * char_per_line + i]);
      // }
    }

    if (plot_num)
    {
      // #pragma omp critical
      // {
      fprintf(fp, " \n");
      // }
    }

    if (quali_num == 0)
    {
      // #pragma omp parallel for shared(hv, k) num_threads(omp_get_max_threads()) schedule(dynamic)
      for (pn = 0; pn < plot_num; pn++)
      {
        fprintf(fp, "                      ");
        for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
        {
          if (!(i % 10))
            fprintf(fp, " ");
          if (plot[k * char_per_line + i] > pn)
            fprintf(fp, "*");
          else
            fprintf(fp, " ");
        }
        fprintf(fp, " \n");

        if (plot_num == 1)
          fprintf(fp, " \n");
      }
    }

    if (quali_num)
    {
      for (i = 0; i < SEQ_NAME_LEN; i++)
      {
        fprintf(fp, " ");
      }

      fprintf(fp, "          ");
      for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
      {
        if (!(i % 10))
          fprintf(fp, " ");
        pl_int = 9 * plot[k * char_per_line + i] / plot_num;
        fprintf(fp, "%d", pl_int);
      }
      fprintf(fp, " \n");
    }

    /***********************************************************************

    fprintf(fp, " \n");
    if( wgt_type > 1 ) {
      for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) {
        fprintf(fp," ");
      }

      fprintf(fp,"  plus    ");
      for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) {
        if( !(i%10) )fprintf(fp, " ");
        fprintf(fp, "%d", plus_count[ k * char_per_line + i ] );
      }
      fprintf(fp, " \n");
    }

    if( wgt_type > 1 ) {
      for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) {
        fprintf(fp," ");
      }

      fprintf(fp,"  minus   ");
      for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) {
        if( !(i%10) )fprintf(fp, " ");
        fprintf(fp, "%d", minus_count[ k * char_per_line + i ] );
      }
      fprintf(fp, " \n");
    }

    if( wgt_type > 1 ) {
      for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) {
        fprintf(fp," ");
      }

      fprintf(fp,"  nuc     ");
      for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) {
        if( !(i%10) )fprintf(fp, " ");
        fprintf(fp, "%d", nuc_count[ k * char_per_line + i ] );
      }
      fprintf(fp, " \n");
      fprintf(fp, " \n");
    }

    ************************************************************************/

    if (wgt_type_plot)
      if (wgt_type == 3)
      {

        fprintf(fp, "sim. level");

        for (i = 0; i < SEQ_NAME_LEN; i++)
        {
          fprintf(fp, " ");
        }

        for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
        {
          if (!(i % 10))
            fprintf(fp, " ");
          sim_char = '.';

          if (frg_involved[k * char_per_line + i])
          {

            f_inv = frg_involved[k * char_per_line + i];
            frac_plus = plus_count[k * char_per_line + i] / f_inv;
            frac_minus = minus_count[k * char_per_line + i] / f_inv;
            frac_nuc = nuc_count[k * char_per_line + i] / f_inv;

            if (frac_plus > weak_wgt_type_thr)
              if (crick_strand)
                sim_char = 'f';
              else
                sim_char = 'p';
            if (frac_plus > strong_wgt_type_thr)
              if (crick_strand)
                sim_char = 'F';
              else
                sim_char = 'P';
            if (frac_minus > weak_wgt_type_thr)
              sim_char = 'r';
            if (frac_minus > strong_wgt_type_thr)
              sim_char = 'R';

            if (frac_nuc > weak_wgt_type_thr)
              sim_char = 'n';
            if (frac_nuc > strong_wgt_type_thr)
              sim_char = 'N';
          }
          fprintf(fp, "%c", sim_char);
        }
        fprintf(fp, " \n");
        fprintf(fp, " \n");
      }

    fprintf(fp, " \n");

  } /*   for(k=0;k<=lc;k++)  */

  input_offset = 0;

  if (fasta_file)
  {
    for (sv = 0; sv < seqnum; sv++)
    {
      if (input_offset >= mapped_fasta->sb.st_size)
      {
        printf("Error: Could not find sequence name for sequence %d\n", sv);
        exit(1);
      }

      fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &sv, input_offset);
      char *seq_name_sv = seq_name_val->data;

#pragma omp atomic write
      input_offset = seq_name_val->line_offset;

      // fprintf(fp2, ">%s", full_name[sv]);
      fprintf(fp2, ">%s", seq_name_sv);
      for (i = 0; i < endlen; i++)
      {
        if (!(i % 50))
          fprintf(fp2, "\n");
        fprintf(fp2, "%c", endseq[sv][i]);
      }

      fprintf(fp2, "\n ");
      if (sv < (seqnum - 1))
        fprintf(fp2, "\n");

      free(seq_name_sv);
      // free(seq_name_val->data);
      free(seq_name_val);
    }
  }

  input_offset = 0;

  if (cw_file)
  {
    block_no = 0;

    fprintf(fp4, "DIALIGN 2.1 multiple sequence alignment \n\n");
    fprintf(fp4, "// \n\n\n");

    while (block_no * 60 < endlen)
    {
      char_no = mini2(60, (endlen - block_no * 60));
      for (sv = 0; sv < seqnum; sv++)
      {
        if (input_offset >= mapped_fasta->sb.st_size)
        {
          printf("Error: Could not find sequence name for sequence %d\n", sv);
          exit(1);
        }

        fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &sv, input_offset);
        char *seq_name_sv = seq_name_val->data;

#pragma omp atomic write
        input_offset = seq_name_val->line_offset;

        // fprintf(fp4, "%s ", seq_name[sv]);
        fprintf(fp4, "%s ", seq_name_sv);
        for (i = 0; i < char_no; i++)
          fprintf(fp4, "%c", endseq[sv][block_no * 60 + i]);
        fprintf(fp4, "\n");
        free(seq_name_sv);
        // free(seq_name_val->data);
        free(seq_name_val);
      }
      fprintf(fp4, "\n\n");
      block_no++;
    }
  }

  input_offset = 0;

  if (msf_file)
  {
    msf_lines = endlen / 50;
    if (endlen % 50)
      msf_lines = msf_lines + 1;

    fprintf(fp3, "DIALIGN 2\n\n\n");
    fprintf(fp3, "   MSF: %d \n\n", endlen);

    for (sv = 0; sv < seqnum; sv++)
    {
      if (input_offset >= mapped_fasta->sb.st_size)
      {
        printf("Error: Could not find sequence name for sequence %d\n", sv);
        exit(1);
      }

      fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &sv, input_offset);
      char *seq_name_sv = seq_name_val->data;

      fasta_len_value *len_val = get_seqlen_mmapped_fasta(mapped_fasta, &sv, input_offset);
      size_t seq_len_sv = len_val->len;

#pragma omp atomic write
      input_offset = seq_name_val->line_offset + seq_len_sv + 1;

      // fprintf(fp3, " Name: %s    Len: %d \n", seq_name[sv], seqlen[sv]);
      fprintf(fp3, " Name: %s    Len: %d \n", seq_name_sv, seq_len_sv);
      free(seq_name_sv);
      // free(seq_name_val->data);
      free(seq_name_val);
    }
    fprintf(fp3, "\n// \n\n");

    input_offset = 0;

    for (lv = 0; lv < msf_lines; lv++)
    {
      add = lv * 50;
      max_p = mini2(endlen - add, 50);

      for (sv = 0; sv < seqnum; sv++)
      {
        if (input_offset >= mapped_fasta->sb.st_size)
        {
          printf("Error: Could not find sequence name for sequence %d\n", sv);
          exit(1);
        }

        fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &sv, input_offset);
        char *seq_name_sv = seq_name_val->data;

#pragma omp atomic write
        input_offset = seq_name_val->line_offset;

        // fprintf(fp3, "%s", seq_name[sv]);
        fprintf(fp3, "%s", seq_name_sv);
        free(seq_name_sv);
        // free(seq_name_val->data);
        free(seq_name_val);

        for (i = 0; i < 4; i++)
          fprintf(fp3, " ");

        for (i = 0; i < max_p; i++)
        {
          if (!(i % 10))
            fprintf(fp3, " ");
          if (endseq[sv][add + i] == '-')
            fprintf(fp3, ".");
          else
            fprintf(fp3, "%c", endseq[sv][add + i]);
        }
        fprintf(fp3, "\n");
      }
      fprintf(fp3, "\n\n");
    }
  }

  if ((seqnum > 2) && (ref_seq == 0))
  {
    fprintf(fp, "\n \n \n   Sequence tree:\n");
    fprintf(fp, "   ==============\n\n");

    if (!strcmp(clust_sim, "av"))
      fprintf(fp, "Tree constructed using UPGMA");
    fprintf(fp, "based on DIALIGN fragment weight scores");

    if (!strcmp(clust_sim, "max"))
      fprintf(fp, "Tree constructed using maximum linkage clustering");

    if (!strcmp(clust_sim, "min"))
      fprintf(fp, "Tree constructed using minimum linkage clustering");

    fprintf(fp, "\n \n%s", upg_str);
  }

  fprintf(fp, "\n \n \n");

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
    free(hseq[hv]);

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
    free(endseq[hv]);

  if (fragno > 0)
    free(fragments);

  free(plot);

  free(weight_count);

  // if (fragno >= 0)
  // {

  //   // for (hv = 0; hv < seqnum; hv++)
  //   // {
  //   //   begin[hv] = seqlen[hv];
  //   //   end[hv] = 1;
  //   // }

  //   // if (fragno > 0)
  //   //   if ((fragments = calloc(fragno, sizeof(struct multi_frag))) == NULL)
  //   //   {
  //   //     printf("not enough memory available for fragments!\n");
  //   //     fprintf(fp, "not enough memory available for fragments!\n");
  //   //     exit(1);
  //   //   }

  //   // for (hv = 1; hv <= fragno; hv++)
  //   // {
  //   //   fragments[hv - 1] = *dia;
  //   //   dia = dia->next;
  //   // }

  //   // for (hv = 0; hv < fragno; hv++)
  //   //   for (j = 0; j < 2; j++)
  //   //   {
  //   //     mini(&begin[fragments[hv].s[j]], fragments[hv].b[j]);
  //   //     maxi(&end[fragments[hv].s[j]], fragments[hv].b[j] + fragments[hv].ext);
  //   //   }

  //   // for (hv = 0; hv < seqnum; hv++)
  //   // {
  //   //   begin[hv] = 1;
  //   //   end[hv] = seqlen[hv] + 1;
  //   // }

  //   // b_size = 0;

  //   // for (i = 0; i < seqnum; i++)
  //   // {
  //   //   b_len[i] = end[i] - begin[i];
  //   //   maxi(&b_size, b_len[i]);
  //   // }

  //   // for (i = 0; i < seqnum; i++)
  //   //   for (hv = 0; hv < b_len[i]; hv++)
  //   //     shift[i][begin[i] + hv] = hv;

  //   // shift_cond = 1;

  //   // while (shift_cond)
  //   // {
  //   //   shift_cond = 0;

  //   //   for (hv = 0; hv < fragno; hv++)
  //   //     for (j = 0; j < 2; j++)
  //   //     {
  //   //       k = (j + 1) % 2;
  //   //       s1 = fragments[hv].s[j];
  //   //       s2 = fragments[hv].s[k];
  //   //       b1 = fragments[hv].b[j];
  //   //       b2 = fragments[hv].b[k];
  //   //       e = fragments[hv].ext;

  //   //       for (l = e - 1; l >= 0; l--)
  //   //       {
  //   //         dif = shift[s2][b2 + l] - shift[s1][b1 + l];
  //   //         if (dif > 0)
  //   //         {
  //   //           new_shift(s1, b1 + l, dif);
  //   //           shift_cond = 1;
  //   //         }
  //   //       }
  //   //     }
  //   // } /*  while (shift_cond)  */

  //   // endlen = 0;

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   maxi(&endlen, shift[hv][end[hv] - 1] + 1);

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   if ((endseq[hv] = calloc(endlen, sizeof(char))) == NULL)
  //   //   {
  //   //     printf(" not enough memory available for printing results!\n");
  //   //     fprintf(fp, " not enough memory available");
  //   //     fprintf(fp, " for printing results!\n");
  //   //     exit(1);
  //   //   }

  //   // if ((inv_shift = (int **)calloc(seqnum, sizeof(int *))) == NULL)
  //   // {
  //   //   printf("not enough memory available for `inv_shift' !!!!\n");
  //   //   fprintf(fp, "not enough memory available for `inv_shift' !\n");
  //   //   exit(1);
  //   // }

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   if ((inv_shift[hv] = (int *)calloc((endlen + 2), sizeof(int))) == NULL)
  //   //   {
  //   //     printf("not enough memory available for `inv_shift' !!!!\n");
  //   //     fprintf(fp, "not enough memory available for `inv_shift' !\n");
  //   //     exit(1);
  //   //   }

  //   // if ((clear_seq = (char *)calloc((endlen + 1), sizeof(char))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `clear_seq' !  \n \n");
  //   //   exit(1);
  //   // }

  //   // if ((weight_count =
  //   //          (float *)calloc((endlen + 2), sizeof(float))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `weight_count' !\n \n");
  //   //   exit(1);
  //   // }

  //   // if ((plot = (float *)calloc((endlen + 2), sizeof(float))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `plot' ! \n \n");
  //   //   exit(1);
  //   // }

  //   // if ((plus_count =
  //   //          (int *)calloc((endlen + 2), sizeof(int))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `plus_count' !\n \n");
  //   //   exit(1);
  //   // }

  //   // if ((minus_count =
  //   //          (int *)calloc((endlen + 2), sizeof(int))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `minus_count' !\n \n");
  //   //   exit(1);
  //   // }

  //   // if ((nuc_count =
  //   //          (int *)calloc((endlen + 2), sizeof(int))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `nuc_count' !\n \n");
  //   //   exit(1);
  //   // }

  //   // if ((frg_involved =
  //   //          (int *)calloc((endlen + 2), sizeof(int))) == NULL)
  //   // {
  //   //   printf(" problems with memory allocation for `frg_involved ' !\n \n");
  //   //   exit(1);
  //   // }

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   for (p = 1; p <= seqlen[hv]; p++)
  //   //     inv_shift[hv][shift[hv][p]] = p;

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   if ((hseq[hv] = calloc((maxlen + 1), sizeof(char))) == NULL)
  //   //   {
  //   //     printf("not enough memory available for printing results! \n");
  //   //     fprintf(fp, "not enough memory available");
  //   //     fprintf(fp, " for printing results! \n");
  //   //     exit(1);
  //   //   }
  //   /*
  //   printf("endlen = %d \n\n", endlen);
  //   */

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   for (i = 0; i < endlen; i++)
  //   //     endseq[hv][i] = gap_char;

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   for (i = begin[hv]; i < end[hv]; i++)
  //   //     hseq[hv][i] = tolower(seq[hv][i]);

  //   // for (hv = 0; hv < fragno; hv++)
  //   //   for (k = 0; k < 2; k++)
  //   //     for (i = fragments[hv].b[k]; i < fragments[hv].b[k] + fragments[hv].ext; i++)
  //   //       hseq[fragments[hv].s[k]][i] = seq[fragments[hv].s[k]][i];

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   for (i = begin[hv]; i < end[hv]; i++)
  //   //     endseq[hv][shift[hv][i]] = hseq[hv][i];

  //   // for (i = 0; i < endlen; i++)
  //   //   clear_seq[i] = ' ';

  //   // for (p = 0; p < endlen; p++)
  //   // {
  //   //   s1 = 0;
  //   //   while (
  //   //       (endseq[s1][p] == tolower(endseq[s1][p])) && (s1 < (seqnum - 1)) /* no capital letter */
  //   //   )
  //   //     s1++;

  //   //   if (s1 < (seqnum - 1))
  //   //   {
  //   //     for (s2 = s1 + 1; s2 < seqnum; s2++)
  //   //     {
  //   //       if (endseq[s2][p] != tolower(endseq[s2][p]))
  //   //       /* endseq[s2][p] capital letter */
  //   //       {
  //   //         aligned = alignedPositions(clos, s1, inv_shift[s1][p], s2,
  //   //                                    succFrontier(clos, s1, inv_shift[s1][p], s2));

  //   //         if (!aligned)
  //   //           /* i.e.endseq[s1][p] not aligned with end seq[s2][p]*/
  //   //           clear_seq[p] = ambi_char;
  //   //       }
  //   //     }
  //   //   }
  //   // }

  //   // if (mask)
  //   //   for (sv = 0; sv < seqnum; sv++)
  //   //     for (hv = 0; hv < endlen; hv++)
  //   //       if (endseq[sv][hv] != gap_char)
  //   //         if (endseq[sv][hv] == tolower(endseq[sv][hv]))
  //   //           endseq[sv][hv] = '*';

  //   // if (col_score)
  //   // {
  //   //   fprintf(fp_col_score, "# 1 %d \n", endlen);
  //   //   fprintf(fp_col_score, "# %s \n", upg_str);
  //   // }

  //   // plot_calc(num_all_it_dia, endlen, weight_count, plot, all_it_dia, fp_col_score);

  //   // wgt_type_count(num_all_it_dia, endlen, plus_count, minus_count, nuc_count, frg_involved, all_it_dia);

  //   // lc = (endlen - 1) / char_per_line;
  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   first_pos[hv] = begin[hv];

  //   // for (k = 0; k <= lc; k++)
  //   // {
  //   //   for (hv = 0; hv < seqnum; hv++)
  //   //   {
  //   //     fprintf(fp, "%s", seq_name[hv]);

  //   //     fprintf(fp, "%8d  ", first_pos[hv]);

  //   //     for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
  //   //     {
  //   //       if (!(i % 10))
  //   //         fprintf(fp, " ");
  //   //       fprintf(fp, "%c", endseq[hv][k * char_per_line + i]);
  //   //       if (endseq[hv][k * char_per_line + i] != gap_char)
  //   //         first_pos[hv]++;
  //   //     }
  //   //     fprintf(fp, " \n");
  //   //   }

  //   //   fprintf(fp, "         ");
  //   //   for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
  //   //   {
  //   //     if (!(i % 10))
  //   //       fprintf(fp, " ");
  //   //     fprintf(fp, "%c", clear_seq[k * char_per_line + i]);
  //   //   }

  //   //   if (plot_num)
  //   //     fprintf(fp, " \n");

  //   //   if (quali_num == 0)
  //   //     for (pn = 0; pn < plot_num; pn++)
  //   //     {
  //   //       fprintf(fp, "                      ");
  //   //       for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
  //   //       {
  //   //         if (!(i % 10))
  //   //           fprintf(fp, " ");
  //   //         if (plot[k * char_per_line + i] > pn)
  //   //           fprintf(fp, "*");
  //   //         else
  //   //           fprintf(fp, " ");
  //   //       }
  //   //       fprintf(fp, " \n");

  //   //       if (plot_num == 1)
  //   //         fprintf(fp, " \n");
  //   //     }

  //   //   if (quali_num)
  //   //   {
  //   //     for (i = 0; i < SEQ_NAME_LEN; i++)
  //   //     {
  //   //       fprintf(fp, " ");
  //   //     }

  //   //     fprintf(fp, "          ");
  //   //     for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
  //   //     {
  //   //       if (!(i % 10))
  //   //         fprintf(fp, " ");
  //   //       pl_int = 9 * plot[k * char_per_line + i] / plot_num;
  //   //       fprintf(fp, "%d", pl_int);
  //   //     }
  //   //     fprintf(fp, " \n");
  //   //   }

  //   //   /***********************************************************************

  //   //   fprintf(fp, " \n");
  //   //   if( wgt_type > 1 ) {
  //   //     for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) {
  //   //       fprintf(fp," ");
  //   //     }

  //   //     fprintf(fp,"  plus    ");
  //   //     for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) {
  //   //       if( !(i%10) )fprintf(fp, " ");
  //   //       fprintf(fp, "%d", plus_count[ k * char_per_line + i ] );
  //   //     }
  //   //     fprintf(fp, " \n");
  //   //   }

  //   //   if( wgt_type > 1 ) {
  //   //     for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) {
  //   //       fprintf(fp," ");
  //   //     }

  //   //     fprintf(fp,"  minus   ");
  //   //     for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) {
  //   //       if( !(i%10) )fprintf(fp, " ");
  //   //       fprintf(fp, "%d", minus_count[ k * char_per_line + i ] );
  //   //     }
  //   //     fprintf(fp, " \n");
  //   //   }

  //   //   if( wgt_type > 1 ) {
  //   //     for( i = 0 ; i < SEQ_NAME_LEN ; i++ ) {
  //   //       fprintf(fp," ");
  //   //     }

  //   //     fprintf(fp,"  nuc     ");
  //   //     for( i = 0 ; i < mini2( char_per_line , endlen-k*char_per_line ) ; i++ ) {
  //   //       if( !(i%10) )fprintf(fp, " ");
  //   //       fprintf(fp, "%d", nuc_count[ k * char_per_line + i ] );
  //   //     }
  //   //     fprintf(fp, " \n");
  //   //     fprintf(fp, " \n");
  //   //   }

  //   //   ************************************************************************/

  //   //   if (wgt_type_plot)
  //   //     if (wgt_type == 3)
  //   //     {

  //   //       fprintf(fp, "sim. level");

  //   //       for (i = 0; i < SEQ_NAME_LEN; i++)
  //   //       {
  //   //         fprintf(fp, " ");
  //   //       }

  //   //       for (i = 0; i < mini2(char_per_line, endlen - k * char_per_line); i++)
  //   //       {
  //   //         if (!(i % 10))
  //   //           fprintf(fp, " ");
  //   //         sim_char = '.';

  //   //         if (frg_involved[k * char_per_line + i])
  //   //         {

  //   //           f_inv = frg_involved[k * char_per_line + i];
  //   //           frac_plus = plus_count[k * char_per_line + i] / f_inv;
  //   //           frac_minus = minus_count[k * char_per_line + i] / f_inv;
  //   //           frac_nuc = nuc_count[k * char_per_line + i] / f_inv;

  //   //           if (frac_plus > weak_wgt_type_thr)
  //   //             if (crick_strand)
  //   //               sim_char = 'f';
  //   //             else
  //   //               sim_char = 'p';
  //   //           if (frac_plus > strong_wgt_type_thr)
  //   //             if (crick_strand)
  //   //               sim_char = 'F';
  //   //             else
  //   //               sim_char = 'P';
  //   //           if (frac_minus > weak_wgt_type_thr)
  //   //             sim_char = 'r';
  //   //           if (frac_minus > strong_wgt_type_thr)
  //   //             sim_char = 'R';

  //   //           if (frac_nuc > weak_wgt_type_thr)
  //   //             sim_char = 'n';
  //   //           if (frac_nuc > strong_wgt_type_thr)
  //   //             sim_char = 'N';
  //   //         }
  //   //         fprintf(fp, "%c", sim_char);
  //   //       }
  //   //       fprintf(fp, " \n");
  //   //       fprintf(fp, " \n");
  //   //     }

  //   //   fprintf(fp, " \n");

  //   // } /*   for(k=0;k<=lc;k++)  */

  //   // if (fasta_file)
  //   // {
  //   //   for (sv = 0; sv < seqnum; sv++)
  //   //   {
  //   //     fprintf(fp2, ">%s", full_name[sv]);
  //   //     for (i = 0; i < endlen; i++)
  //   //     {
  //   //       if (!(i % 50))
  //   //         fprintf(fp2, "\n");
  //   //       fprintf(fp2, "%c", endseq[sv][i]);
  //   //     }

  //   //     fprintf(fp2, "\n ");
  //   //     if (sv < (seqnum - 1))
  //   //       fprintf(fp2, "\n");
  //   //   }
  //   // }

  //   // if (cw_file)
  //   // {
  //   //   block_no = 0;

  //   //   fprintf(fp4, "DIALIGN 2.1 multiple sequence alignment \n\n");
  //   //   fprintf(fp4, "// \n\n\n");

  //   //   while (block_no * 60 < endlen)
  //   //   {
  //   //     char_no = mini2(60, (endlen - block_no * 60));
  //   //     for (sv = 0; sv < seqnum; sv++)
  //   //     {
  //   //       fprintf(fp4, "%s ", seq_name[sv]);
  //   //       for (i = 0; i < char_no; i++)
  //   //         fprintf(fp4, "%c", endseq[sv][block_no * 60 + i]);
  //   //       fprintf(fp4, "\n");
  //   //     }
  //   //     fprintf(fp4, "\n\n");
  //   //     block_no++;
  //   //   }
  //   // }

  //   // if (msf_file)
  //   // {
  //   //   msf_lines = endlen / 50;
  //   //   if (endlen % 50)
  //   //     msf_lines = msf_lines + 1;

  //   //   fprintf(fp3, "DIALIGN 2\n\n\n");
  //   //   fprintf(fp3, "   MSF: %d \n\n", endlen);

  //   //   for (sv = 0; sv < seqnum; sv++)
  //   //     fprintf(fp3, " Name: %s    Len: %d \n", seq_name[sv], seqlen[sv]);
  //   //   fprintf(fp3, "\n// \n\n");

  //   //   for (lv = 0; lv < msf_lines; lv++)
  //   //   {
  //   //     add = lv * 50;
  //   //     max_p = mini2(endlen - add, 50);

  //   //     for (sv = 0; sv < seqnum; sv++)
  //   //     {
  //   //       fprintf(fp3, "%s", seq_name[sv]);
  //   //       for (i = 0; i < 4; i++)
  //   //         fprintf(fp3, " ");

  //   //       for (i = 0; i < max_p; i++)
  //   //       {
  //   //         if (!(i % 10))
  //   //           fprintf(fp3, " ");
  //   //         if (endseq[sv][add + i] == '-')
  //   //           fprintf(fp3, ".");
  //   //         else
  //   //           fprintf(fp3, "%c", endseq[sv][add + i]);
  //   //       }
  //   //       fprintf(fp3, "\n");
  //   //     }
  //   //     fprintf(fp3, "\n\n");
  //   //   }
  //   // }

  //   // if ((seqnum > 2) && (ref_seq == 0))
  //   // {
  //   //   fprintf(fp, "\n \n \n   Sequence tree:\n");
  //   //   fprintf(fp, "   ==============\n\n");

  //   //   if (!strcmp(clust_sim, "av"))
  //   //     fprintf(fp, "Tree constructed using UPGMA");
  //   //   fprintf(fp, "based on DIALIGN fragment weight scores");

  //   //   if (!strcmp(clust_sim, "max"))
  //   //     fprintf(fp, "Tree constructed using maximum linkage clustering");

  //   //   if (!strcmp(clust_sim, "min"))
  //   //     fprintf(fp, "Tree constructed using minimum linkage clustering");

  //   //   fprintf(fp, "\n \n%s", upg_str);
  //   // }

  //   // fprintf(fp, "\n \n \n");

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   free(hseq[hv]);

  //   // for (hv = 0; hv < seqnum; hv++)
  //   //   free(endseq[hv]);

  //   // if (fragno > 0)
  //   //   free(fragments);

  //   // free(plot);

  //   // free(weight_count);

  // } /* for(bc=0;bc<1;bc++) */

#pragma omp parallel for num_threads(omp_get_max_threads())
  for (hv = 0; hv < seqnum; hv++)
    free(shift[hv]);

} /*  ali_arrange  */

void para_print_mmap(char *s_f, FILE *fpi, mmapped_file *mapped_fasta)
{
  int p_count = 1;
  int hv, i;

  {

    if (cd_gobics)
    {
      fprintf(fpi, " \n                          CHAOS / DIALIGN  \n");
      fprintf(fpi, "                          ***************\n \n");

      if (BETA)
        fprintf(fpi, "                           beta version\n\n");

      fprintf(fpi, "          Program code written by \n");
      fprintf(fpi, "          Burkhard Morgenstern, Said Abdeddaim and Michael Brudno \n\n");
      fprintf(fpi, "             e-mail contact: ");
      fprintf(fpi, "dialign (at) gobics (dot) de \n \n");
      fprintf(fpi, "          Published research assisted");
      fprintf(fpi, " by CHAOS / DIALIGN should cite:  \n \n");
      fprintf(fpi, "             Michael Brudno et al.");
      fprintf(fpi, " (2003)\n");
      fprintf(fpi, "             \"Fast and sensitive multiple alignment");
      fprintf(fpi, " of large genomic sequences\" \n");
      fprintf(fpi, "             BMC Bioinformatics 4:66 \n");
      fprintf(fpi, "             http://www.biomedcentral.com/1471-2105/4/66 \n\n");
    }
    else
    {
      fprintf(fpi, " \n                           DIALIGN 2.2.1 \n");
      fprintf(fpi, "                           *************\n \n");

      if (BETA)
        fprintf(fpi, "                           beta version\n\n");

      fprintf(fpi, "          Program code written by Burkhard");
      fprintf(fpi, " Morgenstern and Said Abdeddaim \n");
      fprintf(fpi, "             e-mail contact: ");
      fprintf(fpi, "dialign (at) gobics (dot) de \n \n");
      fprintf(fpi, "          Published research assisted");
      fprintf(fpi, " by DIALIGN 2 should cite:  \n \n");
      fprintf(fpi, "             Burkhard Morgenstern");
      fprintf(fpi, " (1999).\n");

      fprintf(fpi, "             DIALIGN 2: improvement of the");
      fprintf(fpi, " segment-to-segment\n             approach");
      fprintf(fpi, " to multiple sequence alignment.\n");
      fprintf(fpi, "             Bioinformatics 15,");
      fprintf(fpi, " 211 - 218. \n\n");
    }

    fprintf(fpi, "          For more information, please visit");
    fprintf(fpi, " the DIALIGN home page at \n\n             ");
    fprintf(fpi, "http://bibiserv.techfak.uni-bielefeld.de/dialign/");
    fprintf(fpi, " \n \n");

    fprintf(fpi, "         ************************************************************\n \n");
  }

  /*

        fprintf(fpi,"   Options:\n");
        fprintf(fpi,"   ========\n \n");


        if( wgt_type )
          fprintf(fpi,"    %2d) nucleic acid sequences aligned \n", p_count++);
        else
          fprintf(fpi,"    %2d) protein sequences aligned \n", p_count++);

        if( wgt_type == 2 )
          {
            fprintf(fpi,"    %2d) translation",p_count++);
            fprintf(fpi," of nucleotide fragments");
            fprintf(fpi," into peptide fragments\n");
          }

        if( wgt_type == 1 )
          {
            fprintf(fpi,"    %2d) no translation of",p_count++);
            fprintf(fpi," of nucleotide fragments");
            fprintf(fpi," into peptide fragments\n");
          }

        if( wgt_type == 3 )
          {
            fprintf(fpi,"    %2d) mixed alignment consisting", p_count++);
            fprintf(fpi," of P-fragments and N-fragments \n");
          }


        if( seqnum > 2 )
        if( overlap_weights )
          fprintf(fpi,"    %2d) overlap weights used \n", p_count++);
        else
          fprintf(fpi,"    %2d) overlap weights NOT used \n", p_count++);

        if( threshold )
          {
            fprintf(fpi,"    %2d) threshold T =", p_count++);
            fprintf(fpi," %2.2f\n",threshold);
          }

        if( mask )
          {
            fprintf(fpi,"    %2d) non-aligned residues masked", p_count++);
            fprintf(fpi," by `*' \n", p_count++);
          }


        if( lgs_option)
          {
            fprintf(fpi,"    %2d) option for long genomic ", p_count++);
            fprintf(fpi,"sequences used \n");
          }

        if( crick_strand )
          {
            fprintf(fpi,"    %2d) translation of Watson and Crick strand \n", p_count++);
          }



        if( lmax != MAX_DIA )
          {
            fprintf(fpi,"    %2d) maximum length of fragments = %d", p_count++ , lmax );
            if( wgt_type == 0)
              fprintf(fpi," residues ");
            if( wgt_type == 1)
              fprintf(fpi," residues ");
            if( wgt_type == 0)
              fprintf(fpi," codons ");
            if( wgt_type == 0)
              fprintf(fpi," codons / residues ");
            fprintf(fpi," \n");
          }



        if( fasta_file )
          fprintf(fpi,"    %2d) separate file in FASTA format \n", p_count++);

        if( msf_file )
          fprintf(fpi,"    %2d) separate file in msf format \n", p_count++);

        if( cw_file )
          fprintf(fpi,"    %2d) separate file in clustal format \n", p_count++);

       if( plot_num )
         {
           fprintf(fpi,"    %2d) %d \"*\" characters", p_count++, plot_num);
           fprintf(fpi," for regions of maximum similarity\n");
         }
  */

  if (online)
  {
    fprintf(fpi, "\n\n    The following options have been used: \n\n");
    fprintf(fpi, "     - sequences are");
    if (wgt_type == 0)
      fprintf(fpi, " protein sequences \n");
    if (wgt_type == 1)
      fprintf(fpi, " nucleic acid sequences without translation option\n");
    if (wgt_type == 2)
      fprintf(fpi, " nucleic acid sequences with translation option\n");
    if (speed_optimized)
      fprintf(fpi, "     - speed optimized,");
    fprintf(fpi, " see user guide for details \n");
    if (anchors)
      fprintf(fpi, "     - anchor points used\n");
    fprintf(fpi, "\n");
  }
  else
    fprintf(fpi, "\n\n   %s \n\n", input_line);

  fprintf(fpi, " \n");

  fprintf(fpi, "   Aligned sequences:          length:\n");
  fprintf(fpi, "   ==================          =======\n \n");

  size_t input_offset = 0;

  for (hv = 0; hv < seqnum; hv++)
  {

    if (input_offset >= mapped_fasta->sb.st_size)
    {
      fprintf(stderr, "Error: input_offset >= sb->st_size\n");
      exit(1);
    }

    fprintf(fpi, " %3d) ", hv + 1);
    // fprintf(fpi, "%s", seq_name[hv]);
    // fprintf(fpi, "         %9d\n", seqlen[hv]);

    fasta_value *seq_name_val = get_seqname_mmapped_fasta(mapped_fasta, &hv, input_offset);
    char *seq_name_hv = seq_name_val->data;

    fasta_len_value *len_val = get_seqlen_mmapped_fasta(mapped_fasta, &hv, input_offset);
    size_t seq_len_hv = len_val->len;

#pragma omp atomic write
    input_offset = seq_name_val->line_offset + seq_len_hv + 1;

    // printf("seq_name_hv: %s\n", seq_name_hv); //DEBUG
    // printf("seq_len_hv: %d\n", seq_len_hv); //DEBUG
    // exit(0);
    fprintf(fpi, "%s", seq_name_hv);
    fprintf(fpi, "         %9d\n", seq_len_hv);

    free(seq_name_hv);
    // free(seq_name_val->data);
    free(seq_name_val);
  }

  fprintf(fpi, "\n   Average seq. length:");
  fprintf(fpi, "      %9.1f \n", av_len);

  fprintf(fpi, "\n\n   Please note that only upper-case letters are");
  fprintf(fpi, " considered to be aligned. \n");

  fprintf(fpi, "\n\n   Alignment (DIALIGN format):\n");
  fprintf(fpi, "   ===========================\n \n");

} /* para_print */
