
/*******************\
*                   *
*     DIALIGN 2     *
*                   *
*     anchor.c      *
*                   *
\*******************/

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
// #include <regex.h>

extern int total_anc_num, *seqlen;
extern int seqnum;
extern char *seq[MAX_SEQNUM];
extern struct multi_frag *anchor_frg;

extern char seq_file[NAME_LEN];
extern char input_name[NAME_LEN];

// extern ow_value get_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset);
extern fasta_len_value *get_seqlen_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset);
extern mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset);
extern line_value get_offset_of_line(mmapped_file *mapped_file, int *line_num, size_t input_offset);
extern line_value get_line_of_offset(mmapped_file *mapped_file, size_t input_offset);
extern int word_count(char *str);

size_t anchor_check(int s1, int s2, int b1, int b2, int l, double scr, size_t fasta_offset, mmapped_file *fasta_mapped)
{

  if (
      (s1 < 1) ||
      (s1 > seqnum))
  {
    printf(" \n\n  wrong sequence # %d in anchoring file\n\n", s1);
    printf("  data set consists only of %d sequences \n\n", seqnum);
    printf("  PROGRAM TERMINATED \n\n");
    exit(1);
  }
  if (
      (s2 < 1) ||
      (s2 > seqnum))
  {
    printf(" \n\n  wrong sequence # %d in anchoring file\n\n", s2);
    printf("  data set consists only of %d sequences \n\n", seqnum);
    printf("  PROGRAM TERMINATED \n\n");
    exit(1);
  }

  if (s1 == s2)
  {
    printf("\n strange data in anchoring file:\n");
    printf(" sequence # %d anchored with itself.\n\n", s1);
    printf("  PROGRAM TERMINATED \n\n");
    exit(1);
  }

  /*
        if(
          ( b1 < 1 ) ||
          ( b1 + l - 1 > seqlen[ s1 - 1 ] )
        ) {
          printf(" \n\n anchor # %d starts", total_anc_num + 1 ) ;
          printf(" at position %d in sequence %d and has a length of %d.\n", b1, s1, l ) ;
          printf(" This does not fit into sequence # %d " , s1 );
          printf(" (sequence length = %d) \n\n", seqlen[ s1 - 1 ] ) ;
          printf("  PROGRAM TERMINATED \n\n" ) ;
          exit( 1 ) ;
        }
  */

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

  // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  // struct stat *sb = &(mmapped_fasta->sb);
  // size_t file_size = sb->st_size;

  // mmapped_file *fasta_mapped = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  // // struct stat *sb_fasta = &(mmapped_fasta->sb);
  // size_t file_size_fasta = fasta_mapped->sb.st_size;

  // printf(" s1:%d s2:%d\n", s1, s2); // DEBUG

  int s1m1 = s1 - 1;
  int s2m1 = s2 - 1;
  fasta_len_value *seqlen_s1m1_val = get_seqlen_mmapped_fasta(fasta_mapped, &s1m1, fasta_offset);
  size_t seqlen_s1m1 = seqlen_s1m1_val->len;

  fasta_len_value *seqlen_s2m1_val = get_seqlen_mmapped_fasta(fasta_mapped, &s2m1, fasta_offset);
  size_t seqlen_s2m1 = seqlen_s2m1_val->len;

  // munmap(mmapped_fasta->mapped_file, file_size_fasta);
  // free(sb_fasta);
  // free(mmapped_fasta);

  if (
      (b1 < 1) ||
      // ( b1 + l - 1 > seqlen[ s1 - 1 ] )
      (b1 + l - 1 > seqlen_s1m1))
  {
    printf(" \n\n  WARNING:");
    printf(" \n\n  anchor # %d starts", total_anc_num + 1);
    printf(" at position %d in sequence %d\n ", b1, s1);
    printf(" and is %d residues in length.\n", l);
    printf("  However, sequence %d", s1);
    // printf(" is only %d residues in length \n\n", seqlen[ s1 - 1 ] ) ;
    printf(" is only %d residues in length \n\n", seqlen_s1m1);
    printf("  PROGRAM TERMINATED \n\n");
    exit(1);
  }

  if (
      (b2 < 1) ||
      // ( b2 + l - 1 > seqlen[ s2 - 1 ] )
      (b2 + l - 1 > seqlen_s2m1))
  {
    printf(" \n\n  WARNING:");
    printf(" \n\n  anchor # %d starts", total_anc_num + 1);
    printf(" at position %d in sequence %d\n ", b2, s2);
    printf(" and is %d residues in length.\n", l);
    printf("  However, sequence %d", s2);
    // printf(" is only %d residues in length \n\n", seqlen[ s2 - 1 ] ) ;
    printf(" is only %d residues in length \n\n", seqlen_s2m1);
    printf("  PROGRAM TERMINATED \n\n");
    exit(1);
  }

  // munmap(fasta_mapped->mapped_file, file_size_fasta);
  // // free(sb);
  // free(fasta_mapped);
  size_t ret_offset = seqlen_s1m1_val->line_offset;
  free(seqlen_s1m1_val);
  free(seqlen_s2m1_val);

  return ret_offset;
}

// void swap_multifrag_anchor(mmapped_file *mapped_file, int line_num1, int line_num2)
// {
//   size_t line1_len, line2_len;
//   size_t line1_offset, line2_offset;
//   // Get pointers to the start of each line based on their line numbers
//   line1_offset = get_offset_of_line(mapped_file, line_num1, &line1_len);
//   line2_offset = get_offset_of_line(mapped_file, line_num2, &line2_len);
//   char *line1 = &mapped_file[line1_offset];
//   char *line2 = &mapped_file[line2_offset];
//   // Check if both lines are valid and exist
//   if (!line1 || !line2)
//   {
//     fprintf(stderr, "Error: One or both lines not found in the file.\n");
//     return;
//   }

//   // Temporary buffer to hold line content for swapping
//   char *temp1 = strndup(line1, line1_len);
//   char *temp2 = strndup(line2, line2_len);

//   // Ensure the lines are of equal length for an in-place swap, otherwise handle resizing
//   if (line1_len == line2_len)
//   {
//     memcpy(line1, temp2, line2_len);
//     memcpy(line2, temp1, line1_len);
//   }
//   else
//   {
//     // Copy the shorter line to the longer line's memory and pad as needed
//     size_t min_len = line1_len < line2_len ? line1_len : line2_len;
//     memcpy(line1, temp2, min_len);
//     memset(line1 + min_len, ' ', line1_len - min_len);

//     memcpy(line2, temp1, min_len);
//     memset(line2 + min_len, ' ', line2_len - min_len);
//   }

//   // Synchronize changes to the file
//   msync(line1, line1_len, MS_SYNC);
//   msync(line2, line2_len, MS_SYNC);

//   // Free temporary buffers
//   free(temp1);
//   free(temp2);
// }

// void match_pattern_in_mmap(mmapped_file *mapped_file) {
//     regex_t regex;
//     const char *pattern = "^10\\t11"; // Use double backslash for escape sequence
//     int reti;

//     // Compile the regex
//     reti = regcomp(&regex, pattern, REG_EXTENDED);
//     if (reti) {
//         fprintf(stderr, "Could not compile regex\n");
//         exit(1);
//     }

//     // Split the mapped file into lines
//     char *line = strtok(mapped_file->mapped_file, "\n");
//     while (line != NULL) {
//         // Execute the regex
//         reti = regexec(&regex, line, 0, NULL, 0);
//         if (!reti) {
//             printf("Matched line: %s\n", line);
//         }
//         else if (reti == REG_NOMATCH) {
//             // No match
//         }
//         else {
//             char msg[100];
//             regerror(reti, &regex, msg, sizeof(msg));
//             fprintf(stderr, "Regex match failed: %s\n", msg);
//         }
//         line = strtok(NULL, "\n");
//     }

//     // Free compiled regex
//     regfree(&regex);
// }

size_t search_with_memcmp(char *search_pattern, int pattern_length, mmapped_file *mapped_file)
{
  char *line_start = mapped_file->mapped_file;
  char *end = mapped_file->mapped_file + mapped_file->sb.st_size;
  char *newline;

  if (line_start < mapped_file || line_start >= mapped_file->mapped_file + mapped_file->sb.st_size)
  {
    fprintf(stderr, "Invalid line_start pointer\n");
    exit(EXIT_FAILURE);
  }
  if (end < mapped_file || end > mapped_file->mapped_file + mapped_file->sb.st_size)
  {
    fprintf(stderr, "Invalid end pointer\n");
    exit(EXIT_FAILURE);
  }

  size_t length = end - line_start;
  if (length < 0)
  {
    fprintf(stderr, "Invalid length calculation: end is before line_start\n");
    exit(EXIT_FAILURE);
  }

  // while ((newline = memchr(line_start, '\n', end - line_start))) {
  //     size_t line_length = newline - line_start;

  //     // Check if the line is long enough to contain "10\t11"
  //     if (line_length >= pattern_length) {
  //         // Compare the start of the line with "10\t11"
  //         if (memcmp(line_start, search_pattern, pattern_length) == 0) {
  //             // printf("Matched line: %.*s\n", (int)line_length, line_start);
  //             return line_start;
  //         }
  //     }

  //     // Move to the next line
  //     line_start = newline + 1;
  // }

  do
  {
    // Look for a newline character in the remaining part of the file
    newline = memchr(line_start, '\n', end - line_start);
    size_t line_length;

    // If newline was found
    if (newline)
    {
      line_length = newline - line_start;

      // Check if the line is long enough to contain the search pattern
      if (line_length >= pattern_length)
      {
        // Compare the start of the line with the search pattern
        if (memcmp(line_start, search_pattern, pattern_length) == 0)
        {
          return line_start - mapped_file->mapped_file; // Return the offset
        }
      }
      // Move to the next line
      line_start = newline + 1;
    }
    else
    {
      // Handle the last line if no newline was found
      line_length = end - line_start;
      if (line_length >= pattern_length && memcmp(line_start, search_pattern, pattern_length) == 0)
      {
        return line_start - mapped_file->mapped_file; // Return the offset
      }
      break; // Exit the loop if no more lines
    }
  } while (line_start < end); // Continue until the end of the mapped file

  // Check for the last line if it doesn't end with a newline
  if (line_start < end)
  {
    if (memcmp(line_start, search_pattern, pattern_length) == 0)
    {
      // printf("Matched line: %.*s\n", (int)(end - line_start), line_start);
      return line_start;
    }
  }

  return -1;
}

struct multi_frag *get_anchor_by_offset(mmapped_file *mapped_fasta, mmapped_file *mapped_anc, size_t offset_anc, size_t fasta_offset)
{
  unsigned int anc_idx = 0;
  int line_len = 0;
  size_t offset = 0, prev_offset = 0;
  if (offset_anc > 0)
  {
    for (offset = 0; offset < offset_anc; offset++)
    {
      // printf(" mapped_anc[%d]: %c offset_anc:%d\n", offset, mapped_anc[offset], offset_anc); // DEBUG
      line_len++;
      if (mapped_anc->mapped_file[offset] == '\n')
      {
        // printf("offset:%d\n", offset); // DEBUG
        prev_offset = offset - line_len + 1;
        line_len = 0;
        anc_idx++;
        // printf("prev_offset:%d mapped_anc[%d + 1]:%c anc_idx:%d anc_num:%d\n", prev_offset, offset, mapped_anc[offset + 1], anc_idx, *anc_num); // DEBUG
      }
    }
  }

  struct multi_frag *current_frg = NULL;
  char *line_map = &(mapped_anc->mapped_file[offset_anc]);
  int anchor_len = strcspn(line_map, "\n");
  char *line = strndup(line_map, anchor_len);

  char *line_copy = strdup(line);
  // printf(" line:%s\n", line); // DEBUG

  int line_word_count = word_count(line_copy);
  if (line_word_count == 6)
  {
    int beg1 = 0, beg2 = 0, seq1 = 0, seq2 = 0, len = 0;
    double wgt = 0.0f;

    // if ((current_frg = (struct multi_frag *)malloc(sizeof(struct multi_frag))) == NULL)
    if ((current_frg = (struct multi_frag *)calloc(1, sizeof(struct multi_frag))) == NULL)
    {
      printf(" problems with memory allocation for `anchor fragments' at offset %d !  \n \n", offset_anc);
      exit(1);
    }

    sscanf(line, "%d %d %d %d %d %g", &seq1, &seq2, &beg1, &beg2, &len, &wgt);

    current_frg->fasta_offset = anchor_check(seq1, seq2, beg1, beg2, len, wgt, fasta_offset, mapped_fasta);
    current_frg->anc_num = seq1;
    // seq1--;
    // seq2--; // Adjust to zero-indexed
    current_frg->s[0] = seq1 - 1; // Adjust to zero-indexed
    current_frg->s[1] = seq2 - 1; // Adjust to zero-indexed
    current_frg->b[0] = beg1;
    current_frg->b[1] = beg2;
    current_frg->ext = len;
    current_frg->weight = wgt;
    current_frg->prev_offset = prev_offset;
    current_frg->line_num = anc_idx + 1;
    current_frg->line_offset = offset_anc;
    seq1 = seq1 - 1;
    seq2 = seq2 - 1;
    int zero_idx = 0;
    // ow_value ow_val = get_ow_mmapped(mapped_ow, &seq1, &seq2, &zero_idx, 0);
    // if (ow_val.value > 0)
    // {
    //   current_frg->ow = ow_val.value;
    // }
    // else
    // {
    //   current_frg->ow = (double)wgt;
    // }
    // current_frg->ow_offset = ow_val.line_offset;
  }
  else
  {
    if (line_word_count != 0)
    {
      printf("\n\n  Anchor file has wrong format. Each line must contain 6 numbers! \n");
      printf("  Anchor file contains line: %s\n", line_copy);
      printf("  PROGRAM TERMINATED \n\n");
      free(line_copy);
      exit(1);
    }
  }

  free(line);
  free(line_copy);

  return current_frg;
}

struct multi_frag *get_anchor_seqs(mmapped_file *mapped_fasta, mmapped_file *mapped_anc, int *seqnum1, int *seqnum2, size_t fasta_offset, int get_next)
{
  // size_t offset = offset_anc;
  struct multi_frag *current_frg = NULL;
  char search_pattern[120];
  sprintf(search_pattern, "%d\t%d", *seqnum1, *seqnum2);
  // int pattern_len = (sizeof(int) * (sizeof(*seqnum1) + (sizeof(*seqnum2)))) + 1;
  size_t line_offset = search_with_memcmp(search_pattern, strlen(search_pattern), mapped_anc);
  if (line_offset < 0)
  {
    printf("Error: anchors for %d and %d could not be found\n", *seqnum1, *seqnum2);
    exit(0);
  }

  printf(" line_offset:%d\n", line_offset); // DEBUG
  // exit(0); //DEBUG

  current_frg = get_anchor_by_offset(mapped_fasta, mapped_anc, line_offset, fasta_offset);
  current_frg->line_offset = line_offset;
  if (*seqnum1 > 2 && *seqnum2 > 2)
  {
    memset(search_pattern, 0, strlen(search_pattern));
    // pattern_len = (sizeof(int) * (sizeof(*(seqnum1) - 1) + (sizeof(*(seqnum2) - 1)))) + 1;
    sprintf(search_pattern, "%d\t%d", *(seqnum1)-1, *(seqnum1)-1);
    size_t prev_line_offset = search_with_memcmp(search_pattern, strlen(search_pattern), mapped_anc);
    current_frg->prev_offset = prev_line_offset;
  }
  if (get_next >= 1)
  {
    int next_seq1 = *(seqnum1)++;
    int next_seq2 = *(seqnum2)++;
    current_frg->next = get_anchor_seqs(mapped_fasta, mapped_anc, &next_seq1, &next_seq2, current_frg->fasta_offset, 0);
  }

  return current_frg;
}

struct multi_frag *get_anchor_by_num(mmapped_file *mapped_fasta, mmapped_file *mapped_anc, int *anc_num, size_t offset_anc, int get_next)
{
  struct multi_frag *current_frg = NULL;
  unsigned int anc_idx = 0;
  int beg1 = 0, beg2 = 0, seq1 = 0, seq2 = 0, len = 0;
  double wgt = 0.0f;

  size_t offset = 0, prev_offset = 0;
  size_t fasta_offset = 0;
  int found_anchor = 0;
  int line_len = 0;
  // Adjust offset if necessary based on `offset_anc`
  // if (offset_anc > 0)
  // {
  for (size_t tmp_offset = 0; tmp_offset < mapped_anc->sb.st_size; tmp_offset++)
  {
    // printf(" mapped_anc[%d]: %c offset_anc:%d\n", offset, mapped_anc[offset], offset_anc); // DEBUG
    line_len++;
    if (mapped_anc->mapped_file[tmp_offset] == '\n')
    {
      // printf("offset:%d\n", offset); // DEBUG
      prev_offset = offset - line_len + 1;
      offset = tmp_offset + 1;
      line_len = 0;
      anc_idx++;
      // printf("prev_offset:%d mapped_anc[%d + 1]:%c anc_idx:%d anc_num:%d\n", prev_offset, offset, mapped_anc[offset + 1], anc_idx, *anc_num); // DEBUG
      // printf("anc_idx:%d anc_num:%d\n", anc_idx, *anc_num); // DEBUG
    }
    if (anc_idx == *anc_num - 1)
    {
      break;
    }
    if (anc_idx > *anc_num - 1)
    {
      offset = 0;
      break;
    }
  }
  // }

  // Allocate and initialize `current_frg`
  if ((current_frg = (struct multi_frag *)calloc(1, sizeof(struct multi_frag))) == NULL)
  {
    printf(" problems with memory allocation for `anchor fragments' %d !  \n \n", *anc_num - 1);
    exit(1);
  }
  // memset(current_frg, 0, sizeof(struct multi_frag)); // Initialize all fields

  // while (offset < mapped_anc->sb.st_size)
  // {
  // if (anc_idx > *anc_num)
  // {
  //   break;
  // }
  char *line_map = &(mapped_anc->mapped_file[offset]);
  int anchor_len = strcspn(line_map, "\n");
  char *line = strndup(line_map, anchor_len);
  printf(" DEBUG anc_idx:%d anc_num:%d offset:%d\n", anc_idx, *anc_num, offset); // DEBUG
  // Only allocate `line_copy` if needed

  if (anc_idx == *anc_num - 1) //|| *anc_num == 0
  {
    char *line_copy = strdup(line);
    printf(" line_copy:%s\n", line_copy); // DEBUG

    int line_word_count = word_count(line_copy);
    if (line_word_count == 6)
    {
      sscanf(line, "%d %d %d %d %d %lf", &seq1, &seq2, &beg1, &beg2, &len, &wgt);

      fasta_offset = anchor_check(seq1, seq2, beg1, beg2, len, wgt, fasta_offset, mapped_fasta);
      // seq1--;
      // seq2--;
      found_anchor = 1;
      current_frg->s[0] = seq1 - 1; // Adjust to zero-indexed
      current_frg->s[1] = seq2 - 1; // Adjust to zero-indexed
      current_frg->b[0] = beg1;
      current_frg->b[1] = beg2;
      current_frg->ext = len;
      current_frg->weight = wgt;
      printf("GABN - file:%s len:%d wgt:%f\n", mapped_anc->file_name, len, current_frg->weight); // DEBUG
      current_frg->line_offset = offset;
      // current_frg->anc_num = anc_num;
      current_frg->line_num = anc_idx + 1;
      current_frg->anc_num = seq1;
      current_frg->fasta_offset = fasta_offset;
      // if (prev_offset - offset > anchor_len + 1)
      // {
      //   printf("here 1\n"); // DEBUG
      //   prev_offset = offset - (anchor_len + 1);
      // }
      current_frg->prev_offset = prev_offset;
      // current_frg->anc_num = *anc_num;

      current_frg->next = NULL;
      unsigned int next_idx = anc_idx + 1;
      if (get_next == 1 && next_idx <= total_anc_num)
      {
        // current_frg->next = get_anchor_by_num(mapped_anc, &next_idx, current_frg->line_offset, 0);
        current_frg->next = get_anchor_seqs(mapped_fasta, mapped_anc, &seq1, &seq2, current_frg->fasta_offset, 0); // seq1 and seq2 are 1 indexed so no need to add them
      }
      printf("GABN - anc_idx: %d anc_num:%d offset:%d wgt:%f total_anc_num:%d \n", anc_idx, *anc_num, offset, current_frg->weight, total_anc_num); // DEBUG
      printf("%d %d %d %d %d %f\n", seq1, seq2, beg1, beg2, len, wgt);                                                                             // DEBUG
      printf(" line:%s\n", line);                                                                                                                  // DEBUG

      // free(line);
      // free(line_copy);

      // return current_frg;
    }
    else
    {
      if (line_word_count != 0)
      {
        printf("\n\n  Anchor file has wrong format. Each line must contain 6 numbers! \n");
        printf("  Anchor file contains line: %s\n", line_copy);
        printf("  PROGRAM TERMINATED \n\n");
        free(line_copy);
        exit(1);
      }
    }
    free(line_copy); // Ensure line_copy is freed if it was allocated
  }

  // Move to the next line
  // prev_offset = offset;
  // offset += anchor_len + 1;
  // anc_idx++;
  // free(line); // Free `line` once per loop iteration
  // }
  free(line);
  // free(line_copy);

  return current_frg;
  // Handle the end-of-file scenario without infinite recursion
  if (found_anchor == 0 && offset_anc != 0)
  {
    // printf("recursing for anc_num:%d input_offset:%d offset:%d\n", *anc_num, offset_anc, offset); // DEBUG
    return get_anchor_by_num(mapped_fasta, mapped_anc, anc_num, 0, get_next);
  }

  // current_frg->s[0] = 1; // Adjust to zero-indexed
  // current_frg->s[1] = 1; // Adjust to zero-indexed
  // current_frg->b[0] = 0;
  // current_frg->b[1] = 0;
  // current_frg->ext = 0;
  // current_frg->weight = -1;
  // current_frg->line_offset = 0;
  // current_frg->anc_num = 0;
  // current_frg->fasta_offset = 0;

  // printf("returning NULL for anc_num:%d input_offset:%d offset:%d\n", *anc_num, offset_anc, offset); // DEBUG
  return NULL;

  // printf("returning EMPTY for anc_num:%d input_offset:%d offset:%d\n", *anc_num, offset_anc, offset); // DEBUG
  // return current_frg;
}

struct multi_frag *get_anchor_empty()
{
  struct multi_frag *current_frg = NULL;
  if ((current_frg = (struct multi_frag *)calloc(1, sizeof(struct multi_frag))) == NULL)
  {
    printf(" problems with memory allocation for `empty anchor fragment' !  \n \n");
    exit(1);
  }
  current_frg->s[0] = 0;
  current_frg->s[1] = 0;
  current_frg->b[0] = 0;
  current_frg->b[1] = 0;
  current_frg->ext = 0;
  current_frg->weight = 0.0f;
  current_frg->line_offset = 0;
  // current_frg->anc_num = anc_num;
  current_frg->line_num = 0;
  current_frg->anc_num = 0;
  current_frg->fasta_offset = 0;

  return current_frg;
}

// struct multi_frag *get_anchor_by_num(char *mapped_anc, struct stat *sb_anc, char *mapped_fasta, struct stat *sb_fasta, int *anc_num, size_t offset_anc, int get_next)
// {
//   struct multi_frag *current_frg;
//   unsigned int anc_idx = 0;
//   int beg1, beg2, seq1, seq2, len;
//   float wgt;

//   size_t offset = 0;
//   size_t fasta_offset = 0;

//   // mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
//   // struct stat *sb_fasta = &(mmapped_fasta->sb);
//   // size_t file_size_fasta = sb_fasta->st_size;

//   if (offset_anc > 0)
//   {
//     // anc_idx = 0;
//     for (offset = 0; offset < offset_anc; offset++)
//     {
//       if (mapped_anc[offset] == '\n')
//         anc_idx++;
//       if (anc_idx > *anc_num - 1)
//       {
//         offset = 0;
//         break;
//       }
//     }
//   }

//   if ((current_frg = (struct multi_frag *)malloc(1 * sizeof(struct multi_frag))) == NULL)
//   {
//     printf(" problems with memory allocation for `anchor fragments' %d !  \n \n", *anc_num - 1);
//     exit(1);
//   }

//   while (offset < sb_anc->st_size)
//   {
//     // line = &mapped_anc[offset];
//     // anc_idx++;

//     char *line_map = &mapped_anc[offset];
//     int anchor_len = strcspn(line_map, "\n");
//     char *line = strndup(line_map, anchor_len);
//     int line_len = strlen(line);

//     // printf(" mapped_anc[0]: %c\n", mapped_anc[0]);                   // DEBUG
//     // printf(" mapped_anc[%d]: %c\n", offset, mapped_anc[offset]);     // DEBUG
//     // printf(" line_map[0]: %c\n", line_map[0]);               // DEBUG
//     // printf(" line_map[%d]: %c\n", offset, line_map[offset]); // DEBUG
//     // printf(" line: %s\n", line);                             // DEBUG
//     // printf(" strlen(line): %d\n", line_len);                 // DEBUG
//     // printf(" anchor_len: %d\n", anchor_len);                 // DEBUG
//     // printf(" offset: %d\n", offset);                         // DEBUG

//     if (anc_idx == *anc_num - 1)
//     {
//       char *line_copy = (char *)malloc((strlen(line) + 1) * sizeof(char));
//       strcpy(line_copy, line);

//       int line_word_count = word_count(line_copy);

//       // printf("line: %s\n", line_copy);                                 // DEBUG
//       // printf("strlen(line): %d\n", strlen(line));                      // DEBUG
//       // printf("line_len: %d\n", line_word_count);                       // DEBUG
//       // printf("line_offset_len: %d\n", line_word_count * sizeof(char)); // DEBUG
//       // printf("anc_idx: %d anc_num:%d offset:%d \n", anc_idx, *anc_num - 1, offset); // DEBUG

//       if (line_word_count == 6)
//       {
//         sscanf(line, "%d %d %d %d %d %f ", &seq1, &seq2, &beg1, &beg2, &len, &wgt);

//         fasta_offset = anchor_check(mapped_fasta, sb_fasta, seq1, seq2, beg1, beg2, len, wgt, fasta_offset);

//         seq1 = seq1 - 1;
//         seq2 = seq2 - 1;

//         current_frg->s[0] = seq1;
//         current_frg->s[1] = seq2;
//         current_frg->b[0] = beg1;
//         current_frg->b[1] = beg2;
//         current_frg->ext = len;
//         current_frg->weight = wgt;
//         current_frg->line_offset = offset;
//         current_frg->anc_num = *anc_num;

//         // if (offset + line_len + 1 >= sb_anc->st_size)
//         // {
//         //   current_frg->next = NULL;
//         // }
//         // else
//         // {

//         // current_frg->next = (struct multi_frag *)
//         //                  calloc( 1 , sizeof(struct multi_frag) );

//         current_frg->next = NULL;
//         if (get_next == 1)
//         {

//           current_frg->next = (struct multi_frag *)
//               malloc(1 * sizeof(struct multi_frag));

//           unsigned int next_idx = anc_idx + 1;
//           current_frg->next = get_anchor_by_num(mapped_anc, sb_anc, mapped_fasta, sb_fasta, &next_idx, current_frg->line_offset, 0);
//         }

//         free(line);
//         free(line_copy);

//         // munmap(mmapped_fasta->mapped_file, file_size_fasta);
//         // free(mmapped_fasta);

//         return current_frg;
//         // break;
//       }
//       else
//       {
//         if (line_word_count != 0)
//         {
//           printf("\n\n  Anchor file has wrong format. ");
//           printf("\n  Each line must contain 6 numbers! \n");
//           printf("\n  Anchor file contains line \n\n");
//           printf("         %s \n", line_copy);
//           printf("  PROGRAM TERMINATED \n\n");
//           exit(1);
//         }
//       }

//       free(line_copy);
//     }

//     // Move to the next line
//     offset += line_len + 1;
//     anc_idx++;
//     free(line);
//   }

//   if (anc_idx > *anc_num - 1 && offset >= sb_anc->st_size && offset_anc != 0)
//   {
//     // munmap(mmapped_fasta->mapped_file, file_size_fasta);
//     // free(mmapped_fasta);
//     return get_anchor_by_num(mapped_anc, sb_anc, mapped_fasta, sb_fasta, anc_num, 0, 1);
//   }
//   printf(" returning NULL\n"); // DEBUG
//   // munmap(mmapped_fasta->mapped_file, file_size_fasta);
//   // free(mmapped_fasta);
//   return NULL;
// }

int get_anchor_count(char *file_name)
{
  total_anc_num = 0;
  char anc_file_name[NAME_LEN];

  // int fd;
  // struct stat sb;
  // char *mapped;
  // //char *line;
  // //size_t len;

  strcpy(anc_file_name, file_name);
  strcat(anc_file_name, ".anc");

  // // Open the file
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

  mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDONLY, PROT_READ, MAP_SHARED, 0);

  size_t offset = 0;
  int sn = -1;

  // char *mmapped_anc_file = mmapped_anc->mapped_file;
  // printf("\n mmapped_anc->sb ptr: %p \n", mmapped_anc->sb);              // DEBUG
  // printf("\n mmapped_anc->sb ptr: %p \n", &(mmapped_anc->sb));           // DEBUG
  // printf("\n mmapped_anc->sb->st_size: %d \n", mmapped_anc->sb.st_size); // DEBUG
  // struct stat *sb = &(mmapped_anc->sb);
  size_t file_size = mmapped_anc->sb.st_size;

  // printf("\n sb->st_size: %d \n", sb->st_size);                         // DEBUG
  // printf("\n mmapped_anc->sb ptr: %p \n", sb);                          // DEBUG
  // printf("\n mmapped_anc->sb ptr: %p \n", &(mmapped_anc->sb));          // DEBUG
  // printf("\n mmapped_anc->sb.st_size: %d \n", mmapped_anc->sb.st_size); // DEBUG

  while (offset < file_size)
  {
    // line = &mapped[offset];

    char *line_map = &(mmapped_anc->mapped_file[offset]);

    // Find the length of the current line
    int len = strcspn(line_map, "\n");
    char *line = strndup(line_map, len);

    int line_word_count = word_count(line);
    int line_len = strlen(line);

    if (line_word_count == 6)
    {
      printf("anc_num:%d line:%s\n", total_anc_num, line); // DEBUG
      total_anc_num++;
    }
    else
    {
      if (line_word_count != 0)
      {
        printf("\n\n  Anchor file has wrong format. ");
        printf("\n  Each line must contain 6 numbers! \n");
        printf("\n  Anchor file contains line \n\n");
        printf("         %s \n", line);
        printf("  PROGRAM TERMINATED \n\n");
        exit(1);
      }
    }

    // Move to the next line
    offset += line_len + 1;
    free(line);
  }

  // Unmap the file after reading
  munmap(mmapped_anc->mapped_file, file_size);
  // free(sb);
  free(mmapped_anc->file_name);
  free(mmapped_anc);
  return total_anc_num;
}

// int multi_anc_read(char *file_name)
// {

//   char anc_file_name[NAME_LEN];
//   FILE *fp;
//   struct multi_frag *current_frg;
//   char line[10000];
//   int i, len, beg1, beg2, sv = 0, wrdl, hv, word_num;
//   int seq1, seq2;
//   float wgt;

//   strcpy(anc_file_name, file_name);
//   strcat(anc_file_name, ".anc");

//   if ((fp = fopen(anc_file_name, "r")) == NULL)
//     erreur("\n\n cannot find file with anchor points \n\n\n");

//   // if ((anchor_frg = (struct multi_frag *)calloc(1, sizeof(struct multi_frag))) == NULL)
//   if ((anchor_frg = (struct multi_frag *)malloc(1 * sizeof(struct multi_frag))) == NULL)
//   {
//     printf(" problems with memory allocation for `anchor fragments' !  \n \n");
//     exit(1);
//   }

//   current_frg = anchor_frg;

//   while (fgets(line, MLINE, fp) != NULL)
//   {

//     if (word_count(line) == 6)
//     {
//       sscanf(line, "%d %d %d %d %d %f ", &seq1, &seq2, &beg1, &beg2, &len, &wgt);

//       anchor_check(seq1, seq2, beg1, beg2, len, wgt);

//       seq1 = seq1 - 1;
//       seq2 = seq2 - 1;

//       current_frg->s[0] = seq1;
//       current_frg->s[1] = seq2;
//       current_frg->b[0] = beg1;
//       current_frg->b[1] = beg2;
//       current_frg->ext = len;
//       current_frg->weight = wgt;

//       // current_frg->next = (struct multi_frag *)
//       //                  calloc( 1 , sizeof(struct multi_frag) );

//       current_frg->next = (struct multi_frag *)
//           malloc(1 * sizeof(struct multi_frag));

//       current_frg = current_frg->next;
//       total_anc_num++;
//     }
//     else
//     {
//       if (word_count(line) != 0)
//       {
//         printf("\n\n  Anchor file has wrong format. ");
//         printf("\n  Each line must contain 6 numbers! \n");
//         printf("\n  Anchor file contains line \n\n");
//         printf("         %s \n", line);
//         printf("  PROGRAM TERMINATED \n\n");
//         exit(1);
//       }
//     }
//   }
// } /* multi_anc_read  */

int multi_anc_read_mmap(char *file_name)
{

  char anc_file_name[NAME_LEN];
  struct multi_frag *current_frg;
  // char line[ 10000 ] ;
  int i, beg1, beg2, sv = 0, wrdl, hv, word_num;
  int seq1, seq2;
  double wgt;

  // int fd;
  // struct stat sb;
  // char *mapped;
  // //char *line;
  // //size_t len;

  strcpy(anc_file_name, file_name);
  strcat(anc_file_name, ".anc");

  // // Open the file
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

  mmapped_file *mmapped_anc = mmap_file(anc_file_name, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  struct stat *sb = &(mmapped_anc->sb);
  size_t file_size = sb->st_size;

  mmapped_file *mmapped_fasta = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
  // struct stat *sb_fasta = &(mmapped_fasta->sb);
  size_t file_size_fasta = mmapped_fasta->sb.st_size;

  // if( ( anchor_frg = ( struct multi_frag * ) calloc( 1 , sizeof( struct multi_frag ) ))
  //   == NULL)
  if ((anchor_frg = (struct multi_frag *)malloc(1 * sizeof(struct multi_frag))) == NULL)
  {
    printf(" problems with memory allocation for `anchor fragments' !  \n \n");
    exit(1);
  }

  current_frg = anchor_frg;

  size_t offset = 0;
  size_t fasta_offset = 0;

  int sn = -1;

  while (offset < file_size)
  {
    // line = &mapped[offset];

    char *line_map = &(mmapped_anc->mapped_file[offset]);

    // Find the length of the current line
    int len = strcspn(line_map, "\n");
    char *line = strndup(line_map, len);

    // char *line_copy = (char *)malloc((strlen(line)) * sizeof(char));
    // strcpy(line_copy, line);

    int line_word_count = word_count(line);
    int line_len = strlen(line);

    printf("line: %s\n", line);                                      // DEBUG
    printf("strlen(line): %d\n", strlen(line));                      // DEBUG
    printf("line_len: %d\n", line_word_count);                       // DEBUG
    printf("line_offset_len: %d\n", line_word_count * sizeof(char)); // DEBUG
    printf("anc_num: %d\n", total_anc_num);                          // DEBUG

    if (line_word_count == 6)
    {
      sscanf(line, "%d %d %d %d %d %lf", &seq1, &seq2, &beg1, &beg2, &len, &wgt);

      fasta_offset = anchor_check(seq1, seq2, beg1, beg2, len, wgt, fasta_offset, mmapped_fasta);

      // seq1 = seq1 - 1;
      // seq2 = seq2 - 1;

      current_frg->s[0] = seq1 - 1;
      current_frg->s[1] = seq2 - 1;
      current_frg->b[0] = beg1;
      current_frg->b[1] = beg2;
      current_frg->ext = len;
      current_frg->weight = wgt;

      current_frg->next = (struct multi_frag *)
          calloc(1, sizeof(struct multi_frag));

      // current_frg->next = (struct multi_frag *)
      //     malloc(1 * sizeof(struct multi_frag));

      current_frg = current_frg->next;
      total_anc_num++;
    }
    else
    {
      if (line_word_count != 0)
      {
        printf("\n\n  Anchor file has wrong format. ");
        printf("\n  Each line must contain 6 numbers! \n");
        printf("\n  Anchor file contains line \n\n");
        printf("         %s \n", line);
        printf("  PROGRAM TERMINATED \n\n");
        exit(1);
      }
    }

    // Move to the next line
    offset += line_len + 1;
    free(line);
  }

  // Unmap the file after reading
  munmap(mmapped_anc->mapped_file, file_size);
  // free(sb);
  free(mmapped_anc->file_name);
  free(mmapped_anc);

  munmap(mmapped_fasta->mapped_file, file_size_fasta);
  free(mmapped_fasta->file_name);
  free(mmapped_fasta);

  return total_anc_num;
}
