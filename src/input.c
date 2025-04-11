
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

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "define.h"
#include "dialign.h"
#include "alig_graph_closure.h"

extern int max_dia, self_comparison;
extern int sim_score[21][21];
extern int max_sim_score;
extern double av_sim_score_pep;
extern double av_sim_score_nuc;
extern char par_dir[NAME_LEN];
extern double **tp400_prot, **tp400_dna, **tp400_trans;
extern int *seqlen;
extern int seqnum;

extern char *itoa(int num, char *buffer, int base);

extern char *dtoa(char *s, double n, int width, double precision);

int word_count(char *str)
{

  short word = 0;
  int i;
  int word_len = 0;

  for (i = 0; i < strlen(str) - 1; i++)
  {
    if ((str[i] != ' ') && (str[i] != '\t'))
    {
      if (!word)
      {
        word_len++;
        word = 1;
      }
    }
    else
      word = 0;
  }
  return (word_len);
}

void exclude_frg_read(char *file_name, int ***exclude_list)
{

  char exclude_file_name[NAME_LEN];
  FILE *fp;
  char line[10000];
  int i, len, beg1, beg2, seq1, seq2;
  int sv = 0, hv, word_num;

  strcpy(exclude_file_name, file_name);
  strcat(exclude_file_name, ".xfr");

  if ((fp = fopen(exclude_file_name, "r")) == NULL)
    erreur("\n\n cannot find file with excluded fragments \n\n\n");

  while (fgets(line, MLINE, fp) != NULL)
  {
    if (strlen(line) > 4)
    {
      sscanf(line, "%d %d %d %d %d", &seq1, &seq2, &beg1, &beg2, &len);

      if (seq1 > seqnum)
      {
        printf("\n\n exclueded fragment makes no sense!\n\n");
        printf(" wrong sequence no %d in fragment\n\n", seq1);
        printf("%d %d %d %d %d \n\n ", seq1, seq2, beg1, beg2, len);
        exit(1);
      }

      if (seq2 > seqnum)
      {
        printf("\n\n    excluded fragment makes no sense!\n\n");
        printf("    wrong sequence no %d in fragment\n\n", seq2);
        printf("    %d %d %d %d %d \n\n", seq1, seq2, beg1, beg2, len);
        exit(1);
      }

      /*
            seq1 = seq1 - 1;
            seq2 = seq2 - 1;
      */

      if (beg1 + len > seqlen[seq1 - 1] + 1)
      {
        printf("\n\n    excluded fragment makes no sense!\n");
        printf("    fragment");
        printf("     \" %d %d %d %d %d \"\n", seq1, seq2, beg1, beg2, len);
        printf("    doesn't fit into sequence %d:\n", seq1);
        printf("    sequence %d has length =  %d\n\n", seq1, seqlen[seq1 - 1]);
        exit(1);
      }

      for (i = 0; i < len; i++)
      {
        exclude_list[seq1 - 1][seq2 - 1][beg1 + i] = beg2 + i;
      }
    }
  }

  fclose(fp);
} /* excluded_frg_read  */

// void ws_remove( char *str ) {
//   int pv = 0, len;

//   // printf("\n here 3.1 \n"); //DEBUG

//   // Skip over leading spaces, tabs, and newlines
//   while( ( str[ pv ] == ' ' ) || ( str[ pv ] == '\t' ) || ( str[ pv ] == '\n' ) ) {
//     pv++;
//     // printf("\n here 3.2 \n"); //DEBUG
//   }

//   // Use memmove to safely remove leading whitespace
//   memmove( str, str + pv, strlen( str + pv ) + 1 );

//   // printf("\n here 3.3 \n"); //DEBUG

//   // Remove trailing spaces, tabs, and newlines
//   len = strlen( str ) - 1;
//   while( len >= 0 && ( str[ len ] == ' ' || str[ len ] == '\t' || str[ len ] == '\n' ) ) {
//     str[ len ] = '\0';
//     len--;
//     // printf("\n here 3.4 \n"); //DEBUG
//   }

//   // printf("\n here 3.5 \n"); //DEBUG
// }

void ws_remove(char *str)
{
  int pv = 0;

  if (str == NULL)
  {
    printf("Error: str is NULL \n");
    exit(1);
  }

  // printf("\n here 3.1 \n"); //DEBUG

  // Skip over leading spaces and tabs
  while ((str[pv] == ' ') || (str[pv] == '\t'))
  {
    pv++;
    // printf("\n here 3.2 \n"); //DEBUG
  }

  // printf("\n here 3.3 \n"); //DEBUG

  // Use memmove to safely move the string in case of overlapping memory
  memmove(str, str + pv, strlen(str + pv) + 1);

  // printf("\n here 3.4 \n"); //DEBUG
}

// void ws_remove( char *str ) {
//   int pv = 0 ;

//   printf("\n here 3.1 \n"); //DEBUG

//   while( ( str[ pv ] == ' ' ) || ( str[ pv ] == '\t' ) ) {
//     pv++ ;
//     printf("\n here 3.2 \n"); //DEBUG
//   }

//   printf("\n here 3.3 \n"); //DEBUG

//   strcpy( str , str + pv );

//   printf("\n here 3.4 \n"); //DEBUG

// }

void n_clean(char *str)
{
  int pv = 0;
  char *char_ptr;

  // printf("\n here 3.4.1 : %s \n", str); // DEBUG

  while ((str[pv] == ' ') ||
         (str[pv] == '\t') ||
         (str[pv] == '>'))
  {
    pv++;
    // printf("\n here 3.4.2 \n"); // DEBUG
  }
  strcpy(str, str + pv);

  // printf("\n here 3.4.3 \n"); // DEBUG

  if ((char_ptr = strchr(str, ' ')) != NULL)
    *char_ptr = '\0';
  if ((char_ptr = strchr(str, '\t')) != NULL)
    *char_ptr = '\0';
  if ((char_ptr = strchr(str, '\n')) != NULL)
    *char_ptr = '\0';

  // printf("\n here 3.4.4 \n"); // DEBUG
}

void fasta_test(char *seq_file)
{

  int test = 1;
  int pv = 0;
  FILE *fp;

  char line[MAX_INPUT_LINE];

  if ((fp = fopen(seq_file, "r")) == NULL)
  {
    printf("\n\n Cannot find sequence file %s \n\n\n", seq_file);
    exit(1);
  }

  while (test)
  {
    fgets(line, MAX_INPUT_LINE, fp);

    ws_remove(line);

    if (line[0] != '\n')
      if (line[0] == '>')
        test = 0;
      else
        erreur("\n\n  file not in FASTA format  \n\n");
  }

  fclose(fp);
}

// char *get_seq_mmapped_fasta(char* mapped_fasta, int *seqnum, struct stat *sb) {
//     char *line;
//     size_t len;
//     char* sequence;

//     size_t offset = 0;
//     int sn = -1;

//     // printf(" here 8.1.1\n"); //DEBUG
//     // printf("*seqnum: %d\n", *seqnum); //DEBUG

//     while (offset < sb->st_size) {
//         char *line_map = &mapped_fasta[offset];

//         // Find the length of the current line
//         len = strcspn(line_map, "\n");
//         line = strndup(line_map, len);
//         // printf("len: %d\n", len); //DEBUG
//         // printf("line0: %c\n", line[0]); //DEBUG
//         // ws_remove( line );
//         // n_clean( line );
//         // printf("line: %s\n", line); //DEBUG

//         // if(line==NULL){
//         //   break;
//         // }

//         if( line[0] == '>' ) {
//           sn++;
//           if(sn>*seqnum){
//             break;
//           }
//           // printf("line0: %c\n", line[0]); //DEBUG
//         }else{
//           if(sn==*seqnum){
//             sequence = ( char * ) malloc( ( strlen( line ) + 3 ) * sizeof ( char ) );
//             strcpy( sequence , line ) ;
//             //// sequence = line;
//             printf("seq: %s\n", sequence); //DEBUG
//             // break;
//             return sequence;
//           }
//         }
//         offset += len + 1;
//     }
//     return sequence;
// }

// char *get_seqname_mmapped_fasta(char* mapped_fasta, int *seqnum, struct stat *sb) {
//     char *line;
//     size_t len;
//     char* sequence_name;

//     size_t offset = 0;
//     int sn = -1;

//     while (offset < sb->st_size) {
//         char *line_map = &mapped_fasta[offset];

//         // Find the length of the current line
//         len = strcspn(line_map, "\n");
//         line = strndup(line_map, len);
//         // ws_remove( line );
//         // n_clean( line );

//         // if(line==NULL){
//         //   break;
//         // }

//         if( line[0] == '>' ) {
//             sn++;

//             if(sn>*seqnum){
//               break;
//             }

//             if(sn==*seqnum){
//               sequence_name = line;
//               // printf("seq_name: %s\n", sequence_name); //DEBUG
//             break;
//             }
//         }

//         offset += len + 1;
//     }
//     return sequence_name;
// }

fasta_len_value *get_seqlen_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset)
{
  char *line = NULL;
  size_t len;
  size_t sequence_len = 0;

  size_t offset = 0;
  int sn = -1;
  fasta_len_value *ret_value = (fasta_len_value *)calloc(1, sizeof(fasta_len_value));

  ret_value->len = 0;
  ret_value->seq_num = 0;
  ret_value->line_offset = -1;

  int found_seq = 0, found_offset = 0;
  // printf("here 8.1\n"); //DEBUG

  if (input_offset > 0)
  {
    for (size_t tmp_offset = 0; tmp_offset <= input_offset; tmp_offset++)
    {
      // printf(" mapped_fasta[%d]: %c\n", offset, mapped_fasta[offset]); // DEBUG
      if (*seqnum == sn)
      {
        break;
      }
      if (mapped_fasta->mapped_file[tmp_offset] == '>')
      {
        offset = tmp_offset;
        sn++;
      }
    }
  }

  // printf(" (a) mapped_fasta[%d]: %c\n", offset, mapped_fasta[offset]); // DEBUG
  // #pragma omp private(line, len, offset, sn)
  while (offset < mapped_fasta->sb.st_size)
  {
    char *line_map = &(mapped_fasta->mapped_file[offset]);

    // Find the length of the current line
    len = strcspn(line_map, "\n");
    // if (line == NULL)
    // {
    line = (char *)malloc((len + 1) * sizeof(char));
    // }
    // else
    // {
    //   line = (char *)realloc((char *)line, len + 1);
    // }
    line = strndup(line_map, len);

    // ws_remove( line );
    // n_clean( line );

    if (line == NULL)
    {
      break;
    }

    if (line[0] == '>')
    {
      free(line);
      if (sn > *seqnum)
      {
        //   //   // free(line);
        break;
      }
      if (sn == *seqnum && found_seq == 0)
      // if (found_seq == 0)
      {
        ret_value->line_offset = offset;
        ret_value->seq_num = sn;
        found_seq = 1;
      }
      sn++;
    }
    else
    {
      if (sn == *seqnum)
      {
        // printf(" strlen(line): %d\n", strlen(line)); // DEBUG
        // char* sequence = ( char * ) malloc( ( strlen( line ) + 3 ) * sizeof ( char ) );
        // strcpy( sequence , line ) ;
        if (found_offset == 0)
        {
          ret_value->line_offset = offset;
          found_offset = 1;
        }
        sequence_len += strlen(line);
        // printf("*seqnum: %d seq_len: %d\n", *seqnum, sequence_len); // DEBUG

        // #pragma omp critical
        // {
        // if (line != NULL)

        // free(line);

        // }
        // printf(" line: %s\n", line);       // DEBUG
        // printf(" *seqnum: %d\n", *seqnum); // DEBUG

        // return sequence_len;
      }
      free(line);
    }
    // #pragma omp critical
    // {
    // if (line != NULL)
    // free(line);
    // }

    offset += len + 1;
  }
  // return sequence_len;
  // #pragma omp critical
  // {
  // if (line != NULL)
  // free(line);
  // }

  ret_value->len = sequence_len;

  // printf(" *seqnum:%d seqnum:%d len:%d\n", *seqnum, ret_value->seq_num, ret_value->len); // DEBUG

  // printf(" ret_value.len: %d\n", ret_value->len); // DEBUG

  // if (ret_value->len == -1 && input_offset < mapped_fasta->sb.st_size && input_offset != 0)
  if (ret_value->len <= 0 && input_offset < mapped_fasta->sb.st_size && input_offset != 0)
  {
    // printf(" here 8.1.1\n"); // DEBUG
    free(ret_value);
    return get_seqlen_mmapped_fasta(mapped_fasta, seqnum, 0);
  }

  return ret_value;
}

fasta_value *get_seq_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset)
{
  // char *line;
  // size_t len;
  char *sequence = NULL; // Initialize sequence as NULL for safety
  size_t sequence_len = 0;
  size_t offset = 0;
  // size_t ret_offset = -1;
  int current_seqnum = -1;

  // sequence = (char *)malloc(1 * sizeof(char));
  // strcpy(sequence, " ");

  int got_line_offset = 0;
  fasta_value *ret_value = (fasta_value *)calloc(1, sizeof(fasta_value));

  ret_value->data = NULL;
  ret_value->seq_num = -1;
  // ret_value->line_offset = -1;
  ret_value->line_offset = 0;

  if (input_offset > 0)
  {
    for (size_t tmp_offset = 0; tmp_offset <= input_offset; tmp_offset++)
    {
      if (*seqnum == current_seqnum)
      {
        break;
      }
      if (mapped_fasta->mapped_file[tmp_offset] == '>')
      {
        offset = tmp_offset;
        current_seqnum++;
      }
    }
  }

  while (offset < mapped_fasta->sb.st_size)
  {
    char *line_map = &(mapped_fasta->mapped_file[offset]);
    // Find the length of the current line
    size_t len = strcspn(line_map, "\n");
    char *line = strndup(line_map, len); // Duplicates the line

    //   // printf(" here 17.5.1\n"); //DEBUG
    //   // printf(" line_out : %s \n", line); // DEBUG

    //   // Handle sequence headers (lines that start with '>')
    if (line[0] == '>')
    {
      //     current_seqnum++;
      // free(line); // Free the header line
      if (current_seqnum > *seqnum)
      {
        break; // Exit once the desired sequence is processed
      }
      if (current_seqnum == *seqnum && got_line_offset == 0)
      {
        ret_value->line_offset = offset;
        got_line_offset = 1;
      }
      current_seqnum++;
    }
    else
    {
      if (current_seqnum == *seqnum)
      {
        // If sequence is found, accumulate it
        // printf(" here 17.5.2\n"); // DEBUG
        // printf(" line: %s \n", line); // DEBUG
        if (sequence == NULL)
        {
          size_t line_len = strlen(line);
          sequence = (char *)calloc(1, (sequence_len + line_len + 2) * sizeof(char));
          strcat(sequence, " ");
          // strcpy(sequence, "");
          strcpy(sequence + sequence_len + 1, line);
          sequence_len += line_len;
        }
        else
        {
          size_t line_len = strlen(line);
          sequence = realloc(sequence, sequence_len + line_len + 1);
          strcpy(sequence + sequence_len, line); // Append the current line
          sequence_len += line_len;
        }
        // free(line); // Free the non header line
      }
    }
    free(line); // Free the line

    offset += len + 1;
  }

  if (sequence == NULL && input_offset < mapped_fasta->sb.st_size && input_offset != 0)
  // if (sequence == " " && input_offset < mapped_fasta->sb.st_size && input_offset != 0)
  {
    free(ret_value);
    return get_seq_mmapped_fasta(mapped_fasta, seqnum, 0);
  }

  // Null terminate the sequence if it's found
  if (sequence != NULL)
  // if (sequence != " ")
  {
    sequence[sequence_len + 1] = '\0';
    // sequence = &(sequence[1]);
  }

  ret_value->data = sequence;
  ret_value->seq_num = current_seqnum;

  // printf(" seq: %s \n", sequence); // DEBUG
  return ret_value;
}

fasta_value *get_seqname_mmapped_fasta(mmapped_file *mapped_fasta, int *seqnum, size_t input_offset)
{
  // char *line = NULL;
  // size_t len;

  // printf(" seqnum for seqname: %d\n", *seqnum); // DEBUG

  // char *sequence_name = NULL; // Initialize to NULL to avoid dangling pointers

  size_t offset = 0;
  int current_seqnum = -1;

  fasta_value *ret_val = (fasta_value *)calloc(1, sizeof(fasta_value)); // = { NULL, -1, -1 };

  ret_val->data = NULL;
  ret_val->seq_num = -1;
  ret_val->line_offset = -1;

  if (input_offset > 0)
  {
    for (size_t tmp_offset = 0; tmp_offset < input_offset; tmp_offset++)
    {
      if (*seqnum == current_seqnum)
      {
        break;
      }
      if (mapped_fasta->mapped_file[tmp_offset] == '>')
      {
        offset = tmp_offset;
        current_seqnum++;
      }
    }
  }

  while (offset < mapped_fasta->sb.st_size)
  {
    char *line_map = &(mapped_fasta->mapped_file[offset]);

    // Find the length of the current line
    size_t len = strcspn(line_map, "\n");
    char *line = strndup(line_map, len); // Allocate memory for the line

    if (line[0] == '>')
    { // Sequence header
      // current_seqnum++;
      if (current_seqnum == *seqnum)
      {
        ret_val->data = (char *)calloc(len + 1, sizeof(char));
        // ret_val.data = strdup(line + 1); // Copy sequence name
        // *line = line[1];             // Skip the '>' character
        strcpy(ret_val->data, line); // Copy sequence name
        ret_val->seq_num = current_seqnum;
        ret_val->line_offset = offset;
        ret_val->data[len] = '\0';
        // ret_val->data++; // Skip the '>' character
        // ret_val->data = &(ret_val->data[1]);
        // *(ret_val->data) = ret_val->data[1];

        // printf(" num:%d seqname:%s\n", current_seqnum, ret_val.data); // DEBUG

        free(line); // Free the original line

        break; // Exit after finding the sequence name
      }
      current_seqnum++;
    }
    free(line); // Free the memory for the line

    offset += len + 1; // Move to the next line
  }

  if (ret_val->data == NULL && ret_val->seq_num < 0 && input_offset != 0)
  {
    // printf(" retrying with 0 offset\n"); // DEBUG
    free(ret_val);
    return get_seqname_mmapped_fasta(mapped_fasta, seqnum, 0);
  }

  // printf(" num:%d ret_val.seq_num:%d ret_val.line_offset:%d *seqnum:%d input_offset:%d\n", current_seqnum, ret_val->seq_num, ret_val->line_offset, *seqnum, input_offset); // DEBUG
  // printf(" seqname:%s\n", ret_val->data); // DEBUG

  return ret_val; // Return the sequence name
}

size_t get_seqcount_mmapped_fasta(mmapped_file *mapped_fasta)
{
  // char *line;
  // size_t len;

  size_t offset = 0;
  int current_seqnum = 0;

  // printf("mapped_fasta[0]: %c\n", mapped_fasta[0]); // DEBUG
  // printf("sb->st_size: %ld\n", sb->st_size);        // DEBUG

  while (offset < mapped_fasta->sb.st_size)
  {
    char *line_map = &(mapped_fasta->mapped_file[offset]);

    // Find the length of the current line
    size_t len = strcspn(line_map, "\n");
    char *line = strndup(line_map, len); // Allocate memory for the line
    offset += len + 1;                   // Move to the next line

    if (line[0] == '>')
    { // Sequence header
      current_seqnum++;
    }

    free(line); // Free the memory for the line
  }

  // printf("seq_count: %d\n", current_seqnum); // DEBUG
  // free(line);
  return current_seqnum; // Return the sequence name
}

// size_t get_seqlen_mmapped_fasta(char* mapped_fasta, int *seqnum, struct stat *sb) {
//     char *line;
//     size_t len;
//     size_t offset = 0;
//     int current_seqnum = -1;

//     while (offset < sb->st_size) {
//         char *line_map = &mapped_fasta[offset];

//         // Find the length of the current line
//         len = strcspn(line_map, "\n");
//         line = strndup(line_map, len);
//         offset += len + 1;

//         if (line[0] == '>') {  // Sequence header
//             current_seqnum++;
//             if (current_seqnum > *seqnum) {
//                 free(line);
//                 break;
//             }
//         } else {
//             if (current_seqnum == *seqnum) {
//                 size_t seqlen = strlen(line);
//                 free(line);  // Free memory
//                 return seqlen;  // Return sequence length
//             }
//         }
//         free(line);  // Free the line
//     }

//     return 0;  // Return 0 if sequence is not found
// }

size_t find_index_in_fileline(char *line_map, struct stat *sb, int *index, size_t *line_offset)
{

  int current_index = 1;
  // int found_idx = 0;
  size_t offset_val = 0;

  // printf(" line_map[0]: %c line_offset:%d \n", line_map[0], *line_offset); // DEBUG

  size_t header_len = strcspn(line_map, ":");
  size_t len = strcspn(line_map, "\n");

  // printf(" line_map[%d]:%c line_map[%d]:%c \n", header_len, line_map[header_len], header_len + 1, line_map[header_len + 1]); // DEBUG
  // printf(" header_len: %d\n", header_len); // DEBUG
  // printf(" len: %d\n", len);               // DEBUG

  if (current_index == *index)
  {
    return *line_offset + header_len + 1;
  }

  char *line = strndup(line_map, len);

  if (line == NULL || len == 0)
  {
    free(line);
    return -1;
  }

  // printf(" line[0]:%c\n", line[0]);                                   // DEBUG
  // printf(" line:%s\n", line);                                         // DEBUG
  // printf(" len:%d\n", len);                                           // DEBUG
  // printf(" line_map[%d]:%c\n", *line_offset, line_map[*line_offset]); // DEBUG
  // printf("(idx) line[%d]:%c\n", *index, line[*index]);                // DEBUG

  // // #pragma omp parallel for shared(line_map, line_offset, index, offset_val) schedule(dynamic) reduction(+ : current_index)
  // for (offset_val = *line_offset; offset_val < sb->st_size; offset_val++) // Move to the next comma or newline
  // {
  //   // if (found_idx == 1)
  //   // {
  //   //   continue;
  //   // }
  //   if (line_map[*line_offset] == ',' && current_index < *index)
  //   {
  //     // #pragma omp atomic update
  //     current_index += 1; // Found a comma, increment index
  //     // found_idx = 1;
  //   }
  // }

  // #pragma omp parallel for shared(line_map, line_offset, index) schedule(static) reduction(+ : current_index)
  for (offset_val = header_len + 1; offset_val < len && current_index <= *index; offset_val++) // Move to the next comma or newline
  {
    // if (found_idx == 1)
    // {
    //   // continue;
    //   break;
    // }
    if (line[offset_val] != ',')
    {
      // #pragma omp atomic update
      current_index += 1; // Found a comma, increment index
                          // #pragma omp atomic write
      // found_idx = 1;
    }
  }

  // if (*index % 2 != 0)
  // {
  //   offset_val++;
  // }

  // while (current_index < *index && *line_offset < sb->st_size)
  // {
  //   if (line_map[*line_offset] == ',')
  //   {
  //     current_index++; // Found a comma, increment index
  //   }
  //   (*line_offset)++; // Move to the next character
  // }

  // Check if we reached the end of the file or line without finding the element
  // if (*line_offset >= sb->st_size || line_map[*line_offset] == '\n' || line_map[*line_offset] == '\0' || line_map[*line_offset] == NULL) //|| current_index != *index
  // {
  //   return -1;
  // }

  // printf(" line:%s\n", line); // DEBUG

  if (*line_offset + offset_val >= sb->st_size || line[offset_val] == NULL) //|| current_index != *index
  {
    free(line);
    return -1;
  }
  else
  {
    // printf("(before) line_map[%d]:%c\n", *line_offset, line_map[*line_offset]); // DEBUG
    // if (current_index == *index + 1 && line[offset_val] == '\n')
    // {
    //   printf("char \\n found at line[%d]:%c\n", offset_val, line[offset_val]); // DEBUG
    //   *line_offset += offset_val - 1;
    // }
    // else
    // {

    *line_offset += offset_val;
    // printf(" line[0]:%c\n", line[0]);                                          // DEBUG
    // printf(" line:%s\n", line);                                                // DEBUG
    // printf(" len:%d\n", len);                                                  // DEBUG
    // printf("(idx) line[%d]:%c\n", *index, line[*index]);                       // DEBUG
    // printf("(after) line_map[%d]:%c\n", *line_offset, line_map[*line_offset]); // DEBUG
    // }
  }
  free(line);
  // printf("current_index:%d *index:%d line[%d]:%c line[%d + 1]:%c\n", current_index, *index, offset_val, line[offset_val], offset_val + 1, line[offset_val + 1]); // DEBUG
  return *line_offset - 1; // Return the index of the element
}

// mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset)
// {
//   mmapped_file *ret_val = (mmapped_file *)calloc(1, sizeof(mmapped_file));
//   int fd = open(file_name, file_mode);

//   printf("file_name: %s\n", file_name); // DEBUG

//   if (fd == -1)
//   {
//     perror("Error opening file");
//     exit(1);
//   }

//   struct stat sb;
//   if (fstat(fd, &sb) == -1)
//   {
//     perror("Error getting file size");
//     exit(1);
//   }

//   // mmap the file into memory
//   char *mapped_file = mmap(NULL, sb.st_size, protocol, map_mode, fd, offset);
//   if (mapped_file == MAP_FAILED)
//   {
//     perror("Error mmapping file");
//     exit(1);
//   }

//   close(fd);
//   ret_val->mapped_file = mapped_file;
//   ret_val->sb = (struct stat *)calloc(1, sizeof(struct stat));
//   *(ret_val->sb) = sb;

//   printf(" mapped_file[0]: %c\n", ret_val->mapped_file[0]); // DEBUG
//   printf(" sb ptr addr: %p\n", &(ret_val->sb));             // DEBUG
//   // printf(" sb ptr addr (int): %d\n", &(ret_val->sb));       // DEBUG
//   printf(" size: %d\n", ret_val->sb->st_size); // DEBUG

//   return ret_val;
// }

// mmapped_file append_line_to_mmap(mmapped_file *mapped_file, const char *line_to_append)
// {
//   // Get the file descriptor from the file name
//   int fd = open(mapped_file->file_name, O_RDWR);
//   if (fd < 0)
//   {
//     perror("Failed to open file");
//     exit(EXIT_FAILURE);
//   }

//   // Get the current size of the file
//   if (fstat(fd, &mapped_file->sb) < 0)
//   {
//     perror("Failed to get file stats");
//     close(fd);
//     exit(EXIT_FAILURE);
//   }
//   size_t current_size = mapped_file->sb.st_size;
//   size_t new_size = current_size + strlen(line_to_append);

//   // Extend the file size
//   if (ftruncate(fd, new_size) < 0)
//   {
//     perror("Failed to extend file size");
//     close(fd);
//     exit(EXIT_FAILURE);
//   }

//   // Unmap the current mapping
//   if (mapped_file->mapped_file && munmap(mapped_file->mapped_file, mapped_file->sb.st_size) < 0)
//   {
//     perror("Failed to unmap memory");
//     close(fd);
//     exit(EXIT_FAILURE);
//   }

//   // Remap the file with the new size
//   char *new_mapping = mmap(NULL, new_size, PROT_READ | PROT_WRITE, mapped_file->map_mode, fd, 0);
//   if (new_mapping == MAP_FAILED)
//   {
//     perror("Failed to mmap");
//     close(fd);
//     exit(EXIT_FAILURE);
//   }

//   // Append the new line to the mapped region
//   memcpy(new_mapping + current_size, line_to_append, strlen(line_to_append));
//   printf(" line_to_append:%s strlen(line_to_append):%d\n", line_to_append, strlen(line_to_append)); // DEBUG
//   // Sync changes to the file
//   if (msync(new_mapping, new_size, MS_SYNC) < 0)
//   {
//     perror("Failed to sync mmap");
//   }

//   // Update the mmapped_file struct
//   mapped_file->mapped_file = new_mapping;
//   mapped_file->sb.st_size = new_size;

//   // Close the file descriptor
//   close(fd);

//   // exit(0); // DEBUG
//   // Return the updated mmapped_file struct
//   return *mapped_file;
// }

mmapped_file *mmap_file(char *file_name, int file_mode, int protocol, int map_mode, size_t offset)
{
  mmapped_file *ret_val = (mmapped_file *)calloc(1, sizeof(mmapped_file) + 1);
  int fd = open(file_name, file_mode);

  if (fd == -1)
  {
    perror("Error opening file");
    exit(1);
    // fd = open(file_name, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
  }

  // ret_val->fd =

  // Get file size using fstat
  if (fstat(fd, &(ret_val->sb)) == -1) // Directly store in ret_val->sb
  {
    perror("Error getting file size");
    exit(1);
  }

  ret_val->file_mode = file_mode;
  ret_val->file_name = (char *)calloc(1, (strlen(file_name) + 1) * sizeof(char));
  strcpy(ret_val->file_name, file_name);
  ret_val->map_mode = map_mode;
  ret_val->protocol = protocol;
  // mmap the file into memory
  ret_val->mapped_file = mmap(NULL, ret_val->sb.st_size, protocol, map_mode, fd, offset);
  if (ret_val->mapped_file == MAP_FAILED)
  {
    perror("Error mmapping file");
    exit(1);
  }

  close(fd);

  // printf("file_name: %s\n", file_name);                     // DEBUG
  // printf(" mapped_file[0]: %c\n", ret_val->mapped_file[0]); // DEBUG
  // printf(" sb ptr addr: %p\n", &(ret_val->sb));             // DEBUG
  // printf("size: %ld\n", ret_val->sb.st_size);               // DEBUG

  return ret_val;
}

mmapped_file *resize_mmap_file(mmapped_file *mapped_file, size_t new_size)
{
  // mmapped_file *ret_val = (mmapped_file *)calloc(1, sizeof(mmapped_file));
  // Unmap current mapping
  munmap(mapped_file->mapped_file, mapped_file->sb.st_size);

  // printf(" Original size:%d New Size:%d\n", mapped_file->sb.st_size, new_size); // DEBUG

  int fd = open(mapped_file->file_name, mapped_file->file_mode);
  if (fd == -1)
  {
    perror("Error opening file");
    exit(1);
  }

  // Resize the file to new size
  if (ftruncate(fd, new_size) == -1)
  {
    perror("Error resizing the file");
    exit(1);
  }

  // Update stat structure and remap the file
  mapped_file->sb.st_size = new_size;
  // ret_val->file_mode = file_mode;
  // ret_val->file_name = (char *)malloc(1 * strlen(file_name) * sizeof(char));
  // strcpy(ret_val->file_name, file_name);
  // ret_val->map_mode = map_mode;
  // ret_val->protocol = protocol;
  mapped_file->mapped_file = mmap(NULL, new_size, mapped_file->protocol, mapped_file->map_mode, fd, 0);
  if (mapped_file->mapped_file == MAP_FAILED)
  {
    perror("Error remapping the file");
    close(fd);
    exit(1);
  }
  close(fd);

  // printf("resized file:%s size:%d\n", mapped_file->file_name, mapped_file->sb.st_size);
  return mapped_file;
}

size_t get_size_from_offset(mmapped_file *mapped_file, size_t input_offset)
{
  size_t i = 0;
  if (input_offset >= mapped_file->sb.st_size)
  {
    return 0;
  }
  i = input_offset;
  while (i < mapped_file->sb.st_size)
  {
    i++;
  }
  return i;
}

// Function to shift contents within the memory-mapped file
mmapped_file *shift_mmap_contents(mmapped_file *mapped_file, size_t start_offset, size_t shift_amount)
{

  // size_t new_size = mapped_file->sb.st_size + shift_amount;

  // // Resize and remap the file
  // mapped_file = resize_mmap_file(mapped_file, new_size);

  // if (shift_amount > 0)
  // {
  //   // Shift right to make space
  //   memmove(mapped_file->mapped_file + start_offset + shift_amount, mapped_file->mapped_file + start_offset, mapped_file->sb.st_size - start_offset);
  // }
  // else if (shift_amount < 0)
  // {
  //   // Shift left to reduce space
  //   memmove(mapped_file->mapped_file + start_offset + shift_amount, mapped_file->mapped_file + start_offset, mapped_file->sb.st_size - start_offset);
  // }

  if (shift_amount > 0)
  {
    // Shift right to make space
    memmove(mapped_file->mapped_file + start_offset + shift_amount,
            mapped_file->mapped_file + start_offset,
            mapped_file->sb.st_size - start_offset - shift_amount); // Adjusted
    // memset(mapped_file->mapped_file + start_offset, ' ', shift_amount);
  }
  else if (shift_amount < 0)
  {
    // Shift left to reduce space
    memmove(mapped_file->mapped_file + start_offset,
            mapped_file->mapped_file + start_offset - shift_amount,
            mapped_file->sb.st_size - start_offset + shift_amount); // Adjusted
    // memset(mapped_file->mapped_file + mapped_file->sb.st_size + shift_amount, ' ', -shift_amount);
  }

  msync(mapped_file->mapped_file, mapped_file->sb.st_size, MS_SYNC);

  return mapped_file;
}

size_t find_in_file(char *mapped_file, struct stat *sb, char *search_pattern, size_t input_offset)
{

  size_t file_offset = 0, line_len = 0;

  if (input_offset > 0)
  {
    file_offset = input_offset;
  }
  // int found = 0;
  // char *search_result = NULL;
  // #pragma omp parallel for shared(mapped_file, sb, search_pattern, file_offset, line_len) schedule(dynamic) // private(search_result)
  //   for (file_offset = 0; file_offset < sb->st_size; file_offset += line_len + 1)
  //   {
  //     // size_t line_len = 0;
  //     if (found == 1)
  //     {
  //       continue;
  //     }
  // #pragma omp atomic write
  //     line_len = strcspn(&mapped_file[file_offset], "\n");
  //     char *search_result = strstr(&mapped_file[file_offset], search_pattern);
  //     if (search_result != NULL)
  //     {
  //       // printf("Found pattern at position: %ld\n", search_result - mapped_file); //DEBUG
  //       // printf(" mapped_file[%d]:%c\n", file_offset, mapped_file[file_offset]); // DEBUG
  //       file_offset = search_result - mapped_file + strlen(search_pattern) - 1;
  //       // printf(" found mapped_file[%d]:%c\n", file_offset, mapped_file[file_offset]); // DEBUG
  //       // return (file_offset);
  //       found = 1;
  // #pragma omp cancel for
  //     }
  //     // file_offset += line_len + 1;
  //     // file_offset += line_len;
  //   }
  char *search_result = NULL;
  while (file_offset < sb->st_size)
  {
    // printf(" (b) mapped_file[%d]:%c\n", file_offset, mapped_file[file_offset]); // DEBUG
    line_len = strcspn(&mapped_file[file_offset], "\n");
    // printf("Checking line: '%.*s'\n", (int)line_len, &mapped_file[file_offset]);
    // printf(" line_len:%d\n", line_len); // DEBUG
    search_result = strstr(&mapped_file[file_offset], search_pattern);
    if (search_result != NULL)
    {
      // printf("Found pattern at position: %ld\n", search_result - mapped_file); //DEBUG
      // printf("Found %c at position: %d strlen(search_pattern):%d\n", mapped_file[search_result - mapped_file], search_result - mapped_file, strlen(search_pattern)); // DEBUG
      file_offset = (size_t)(search_result - mapped_file); //(search_result - mapped_file);
      // free(search_result);
      // printf(" (a) mapped_file[%d]:%c\n", file_offset, mapped_file[file_offset]); // DEBUG
      break;
    }
    file_offset += line_len + 1;
    // printf(" (a) mapped_file[%d]:%c\n", file_offset, mapped_file[file_offset]); // DEBUG
  }

  if (file_offset >= sb->st_size)
  {
    // printf("%d >= %d : returning -2\n", file_offset, sb->st_size); // DEBUG
    return -2;
  }
  // printf("mapped_file[%d]: %c\n", file_offset, mapped_file[file_offset]); // DEBUG
  return (file_offset);
}

line_value get_offset_of_line(mmapped_file *mapped_file, int *line_num, size_t input_offset)
{
  line_value ret_val;
  size_t offset = 0;
  int current_line = 0;
  size_t line_offset = 0;
  size_t len = 0;

  ret_val.line_num = -1;
  ret_val.line_offset = -1;

  // if (input_offset > 0)
  for (size_t tmp_offset = 0; tmp_offset < input_offset; tmp_offset++)
  {
    if (mapped_file->mapped_file[tmp_offset] == '>')
    {
      offset = tmp_offset;
      current_line++;
    }
    if (current_line == *line_num)
      break;
  }

  // printf(" current_line:%d *line_num:%d\n", current_line, *line_num); // DEBUG

  // while (offset < mapped_file->sb.st_size)
  // {
  //   char *line_map = &(mapped_file->mapped_file[offset]);
  //   len = strcspn(line_map, "\n");
  //   if (current_line == *line_num)
  //   {
  //     line_offset = offset;
  //     break;
  //   }
  //   offset += len + 1;
  //   current_line++;
  // }

  ret_val.line_num = current_line;
  ret_val.line_offset = offset;
  return ret_val;
}

line_value get_line_of_offset(mmapped_file *mapped_file, size_t input_offset)
{
  line_value ret_val;
  size_t offset = 0;
  int current_line = 0;

  ret_val.line_num = -1;
  ret_val.line_offset = -1;

  for (size_t tmp_offset = 0; tmp_offset < input_offset && tmp_offset < mapped_file->sb.st_size; tmp_offset++)
  {
    if (mapped_file->mapped_file[tmp_offset] == '>')
    {
      offset = tmp_offset;
      current_line++;
    }
  }

  ret_val.line_num = current_line;

  line_value offset_line = get_offset_of_line(mapped_file, &ret_val.line_num, input_offset);
  ret_val.line_offset = offset_line.line_offset;
  // printf(" CURR_LINE:%d line_offset:%u line_num:%d\n", current_line, ret_val.line_offset, ret_val.line_num); // DEBUG
  return ret_val;
}

// size_t file_length_from_offset(mmapped_file *mapped_file, size_t input_offset)
// {
//   size_t offset = input_offset;
//   size_t len = 0;
//   for (offset = input_offset; offset < mapped_file->sb.st_size; offset++)
//   {
//     len++;
//   }
//   return len;
// }

openpos_value get_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset)
{
  // size_t file_offset = 0, line_len = 0;
  // int i = 0, j = 0;

  openpos_value ret_val;
  ret_val.x = *x;
  ret_val.y = *y;
  ret_val.z = *z;
  ret_val.value = -1;
  ret_val.line_offset = -1;

  // char *search_result = NULL;
  // while (file_offset < sb->st_size)
  // {
  //   size_t line_len = strcspn(&mapped_file[file_offset], "\n");
  //   // printf("Checking line: '%.*s'\n", (int)line_len, &mapped_file[file_offset]);
  //   // printf(" line_len:%zu\n", line_len); // DEBUG
  //   search_result = strstr(&mapped_file[file_offset], search_pattern);
  //   if (search_result != NULL)
  //   {
  //     // printf("Found pattern at position: %ld\n", search_result - mapped_file); //DEBUG
  //     file_offset = search_result - mapped_file;
  //     break;
  //   }
  //   file_offset += line_len + 1;
  // }

  // #pragma omp parallel for shared(mapped_file, sb, search_pattern, search_result, file_offset, line_len) private(line_len, search_result) schedule(static)
  //   for (file_offset = 0; file_offset < sb->st_size; file_offset += line_len + 1)
  //   {
  //     line_len = strcspn(&mapped_file[file_offset], "\n");
  //     search_result = strstr(&mapped_file[file_offset], search_pattern);
  //     if (search_result != NULL)
  //     {
  //       // printf("Found pattern at position: %ld\n", search_result - mapped_file); //DEBUG
  //       file_offset = search_result - mapped_file;
  //       break;
  //     }
  //   }

  //   if (file_offset > sb->st_size || search_result == NULL)
  //   {
  //     return -1;
  //   }

  size_t file_offset = -1; // = find_in_file(mapped_file, sb, search_pattern, input_offset) + 1;

  // if (input_offset == 0)
  // {
  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d:", *x, *y); // Construct search pattern
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }
  // else
  // {
  //   file_offset = input_offset;
  // }

  if ((int)file_offset < 0)
  {
    // printf(" (p) file_offset: %d\n", file_offset); // DEBUG
    if (input_offset > 0)
    {
      // printf(" (p) input_offset: %d\n", input_offset); // DEBUG
      // free(ret_val);
      return get_openpos_mmapped(mapped_file, x, y, z, 0);
    }
    return ret_val;
  }

  ret_val.line_offset = file_offset;

  char *line_map = &(mapped_file->mapped_file[file_offset]);

  // printf(" *x:%d *y:%d z:%d \n", *x, *y, *z); // DEBUG
  // printf(" input_offset: %d\n", input_offset);                             // DEBUG
  // printf(" file_offset: %d\n", file_offset);                               // DEBUG
  // printf(" line_map[0]: %c\n", line_map[0]);                               // DEBUG
  // printf(" mapped_file[%d]: %c\n", file_offset, mapped_file[file_offset]); // DEBUG

  size_t line_offset = find_index_in_fileline(line_map, &(mapped_file->sb), z, &file_offset); // Move past ,

  // printf(" line_offset: %d\n", line_offset); // DEBUG

  if (line_offset < 0)
  {
    return ret_val;
  }
  // size_t line_offset = file_offset + strlen(search_pattern) * sizeof(char) + 1; // Start after the search pattern
  // size_t len = strcspn(line_map, "\n");
  // char *line = strndup(line_map, len);
  // char *line_tokens = strtok(line, ",");
  // unsigned int line_tokens_offset = &line_tokens[*z] - &line_tokens[0] + 2; // [*z ] + 2 is 1 index
  // // unsigned int line_tokens_offset = &line_tokens[*z] - &line_tokens[0] + 4; // [*z]  + 4 is 0 index

  // printf("line: %s\n", line);                                     // DEBUG
  // printf("strlen(search_pattern): %d\n", strlen(search_pattern)); // DEBUG
  // printf("line_offset: %d\n", line_offset);                       // DEBUG
  // printf("line_map[0]: %c\n", line_map[0]);                            // DEBUG
  // printf("line_offset: %d\n", line_offset);                            // DEBUG
  // printf("file_offset: %d\n", file_offset);                            // DEBUG
  // printf("line_map[%d]: %c\n", line_offset, mapped_file[line_offset]);                                   // DEBUG
  // printf("(f+l) line_map[%d]: %c\n", file_offset + line_offset, mapped_file[file_offset + line_offset]); // DEBUG
  // printf("line[line_offset - 1]: %c\n", line[line_offset - 1]);   // DEBUG
  // printf("line[line_offset + 1]: %c\n", line[line_offset + 1]);   // DEBUG
  // printf("line_tokens[0]: %c\n", *z, line_tokens[0]);             // DEBUG
  // printf("line_tokens[%d]: %c\n", *z, line_tokens[*z]);           // DEBUG
  // printf("line_tokens offset: %d\n", line_tokens_offset);         // DEBUG
  // printf("len: %d\n", len);                                       // DEBUG
  // // int current_index = 1;

  // // // Traverse through the line to find the z-th element
  // // while (current_index < *z && line_offset < sb->st_size)
  // // {
  // //   if (line_map[line_offset] == ',')
  // //   {
  // //     current_index++; // Found a comma, increment index
  // //   }
  // //   line_offset++; // Move to the next character
  // // }

  // // // Check if we reached the end of the file or line without finding the element
  // // if (line_offset >= sb->st_size || line_map[line_offset] == '\n')
  // // {
  // //   return -1;
  // // }
  // size_t file_pos = line_offset + line_tokens_offset;

  // int ret_open_pos = -1;
#pragma omp critical
  {
    // int ret_open_pos = line[*z] - '0';

    // ret_open_pos = mapped_file[line_offset] - '0';

    // printf(" z:%d mapped_file[%d]: %c\n", *z, line_offset, mapped_file[line_offset]); // DEBUG
    ret_val.value = mapped_file->mapped_file[line_offset] - '0';

    // ret_open_pos = line_map[line_offset] - '0';
    // ret_open_pos = line_tokens[*z + 1] - '0';
  }

  // printf(" ret_val: %d\n", ret_val.value);            // DEBUG
  // printf(" ret_val (C) : %c\n", ret_val.value + '0'); // DEBUG

  return ret_val;
}

openpos_value set_openpos_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, int *value, size_t input_offset)
{
  // size_t file_offset = 0, line_len = 0;
  // int i = 0, j = 0;
  // char value_char = *value + '0'; // Convert value to character
  char buffer[256];
  char *value_char = NULL;
  value_char = itoa(*value, buffer, 10);
  int len_value_char = strlen(value_char);
  // printf(" value_char: %c\n", *value_char); // DEBUG

  openpos_value ret_val;

  ret_val.x = *x;
  ret_val.y = *y;
  ret_val.z = *z;
  ret_val.value = -1;
  ret_val.line_offset = -1;

  // printf(" search_pattern: %s\n", search_pattern); // DEBUG
  //  line_len = strcspn(mapped_file, "\n");
  // char *search_result = NULL;
  // while (file_offset < sb->st_size)
  // {
  //   size_t line_len = strcspn(&mapped_file[file_offset], "\n");
  //   // printf("Checking line: '%.*s'\n", (int)line_len, &mapped_file[file_offset]);
  //   // printf(" line_len:%zu\n", line_len); // DEBUG
  //   search_result = strstr(&mapped_file[file_offset], search_pattern);
  //   if (search_result != NULL)
  //   {
  //     // printf("Found pattern at position: %ld\n", search_result - mapped_file);
  //     file_offset = search_result - mapped_file;
  //     break;
  //   }
  //   file_offset += line_len + 1;
  // }

  // #pragma omp parallel for shared(mapped_file, sb, search_pattern, search_result, file_offset, line_len) private(line_len, search_result) schedule(static)
  //   for (file_offset = 0; file_offset < sb->st_size; file_offset += line_len + 1)
  //   {
  //     line_len = strcspn(&mapped_file[file_offset], "\n");
  //     search_result = strstr(&mapped_file[file_offset], search_pattern);
  //     if (search_result != NULL)
  //     {
  //       // printf("Found pattern at position: %ld\n", search_result - mapped_file); //DEBUG
  //       file_offset = search_result - mapped_file;
  //       break;
  //     }
  //   }

  //   if (file_offset > sb->st_size || search_result == NULL)
  //   {
  //     return -1;
  //   }

  // while (file_offset < sb->st_size)
  // {
  //   char *line_map = NULL;
  // #pragma omp critical
  //   {
  //     line_map = &mapped_file[file_offset];
  //   }

  // Find the length of the current line
  // size_t len = strcspn(line_map, "\n");
  // char *line = strndup(line_map, len);
  // printf("line: %s\n", line); // DEBUG

  // if (line[0] == '>')
  // {
  // Read i and j index from the line
  // sscanf(line, "> %d %d", &i, &j);
  // }
  // else
  // {
  // if (i == *x && j == *y)
  // {
  // Update the z-th position with the new value
  // memcpy(&mapped_file[offset - len + z], &value_char, 1);
  // size_t line_offset = file_offset;
  // int current_index = 1;
  // if (line_map[line_offset] == "\n")
  // {
  //   line_offset++;
  // }
  // while (current_index < *z)
  // {
  //   // printf(" line_map[line_offset]: %c\n", line_map[line_offset]); // DEBUG
  //   if (line_map[line_offset] == ',')
  //   {
  //     current_index++; // Increment the index when a comma is found
  //   }
  //   line_offset++; // Move to the next character
  // }

  // if (line_map[line_offset] == "\n")
  // {
  //   return -1;
  // }

  size_t file_offset = -1;

  // if (input_offset == 0)
  // {
  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d:", *x, *y); // Construct search pattern
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }
  // else
  // {
  //   file_offset = input_offset;
  // }

  if ((int)file_offset < 0)
  {
    if (input_offset > 0)
    {
      // free(ret_val);
      return set_openpos_mmapped(mapped_file, x, y, z, value, 0);
    }
    return ret_val;
  }

  ret_val.line_offset = file_offset;

  char *line_map = &mapped_file->mapped_file[file_offset];
  size_t line_offset = find_index_in_fileline(line_map, &(mapped_file->sb), z, &file_offset);

  if (line_offset < 0)
  {
    return ret_val;
  }
  // size_t line_offset = file_offset + strlen(search_pattern) * sizeof(char) + 1; // Start after the search pattern
  // size_t len = strcspn(line_map, "\n");
  // char *line = strndup(line_map, len);
  // char *line_tokens = strtok(line, ",");
  // unsigned int line_tokens_offset = &line_tokens[*z] - &line_tokens[0] + 2; // [*z ] + 2 is 1 index
  // // unsigned int line_tokens_offset = &line_tokens[*z] - &line_tokens[0] + 4; // [*z]  + 4 is 0 index

  // printf("line: %s\n", line);                                     // DEBUG
  // printf("strlen(search_pattern): %d\n", strlen(search_pattern)); // DEBUG
  // printf("line_offset: %d\n", line_offset);                       // DEBUG
  // printf("line_map[0]: %c\n", line_map[0]);                            // DEBUG
  // printf("line_offset: %d\n", line_offset);                            // DEBUG
  // printf("file_offset: %d\n", file_offset);                            // DEBUG
  // printf("(b) line_map[0]: %c\n", line_map[0]);                            // DEBUG
  // printf("(b) line_map[%d]: %c\n", line_offset, mapped_file[line_offset]); // DEBUG
  // printf("(f+l) line_map[%d]: %c\n", file_offset + line_offset, mapped_file[file_offset + line_offset]); // DEBUG

  // printf("line[line_offset - 1]: %c\n", line[line_offset - 1]);   // DEBUG
  // printf("line[line_offset + 1]: %c\n", line[line_offset + 1]);   // DEBUG
  // printf("line_tokens[0]: %c\n", *z, line_tokens[0]);             // DEBUG
  // printf("line_tokens[%d]: %c\n", *z, line_tokens[*z]);           // DEBUG
  // printf("line_tokens offset: %d\n", line_tokens_offset);         // DEBUG
  // printf("len: %d\n", len);                                       // DEBUG

#pragma omp critical
  {
    // size_t file_pos = line_offset - len + *z - 1;
    // size_t file_pos = line_offset + line_tokens_offset;
    // printf("file_pos: %d\n", file_pos); // DEBUG
    // printf("line: %s\n", line);                                            // DEBUG
    // printf("value: %d\n", *value);                                         // DEBUG
    // printf("mapped_line_addr: %d\n", &mapped_file[offset + *z]);           // DEBUG
    // printf("mapped_line: %c\n", mapped_file[line_offset]); // DEBUG
    // printf("mapped_line_addr CPT: %d\n", &mapped_file[offset - len + *z]); // DEBUG
    // printf("mapped_line CPT: %c\n", mapped_file[offset - len + *z]);       // DEBUG
    // memmove(&mapped_file[file_pos], &value_char, 1); // ChatGPT

    // memmove(&mapped_file[line_offset], &value_char, len_value_char); // ChatGPT
    mapped_file->mapped_file[line_offset] = *value_char;
    // line_tokens[*z] = *value_char;

    // mapped_file[offset - len + *z] = value_char; // ChatGPT

    // memmove(&mapped_file[offset + *z], &value_char, 1);

    // Synchronize the change to the file
    // msync(mapped_file, sb->st_size, MS_SYNC);
    // msync(&mapped_file[file_pos], sizeof(value_char), MS_SYNC);
    msync(&(mapped_file->mapped_file[line_offset]), sizeof(char) * len_value_char, MS_SYNC);
    // msync(&mapped_file[offset - len + *z], sizeof(char), MS_SYNC);
    // printf("(a) line_map[%d]: %c\n", line_offset, mapped_file[line_offset]); // DEBUG
    ret_val.value = mapped_file->mapped_file[line_offset] - '0';
  }
  // free(line); // Free memory
  // free(search_result);

  // return *value;
  // free(value_char);
  return ret_val;

  // }
  // }
  // file_offset += len + 1;
  // free(line); // Free the line after processing
  // }

  // return -1; // Return -1 if position not found
}

/*
void set_float_in_csv(char *line, size_t start_offset, int element_index, float new_value) {
    char *current_position = line + start_offset;  // Move to starting offset
    int current_index = 0;

    // Locate the start of the target element
    while (*current_position && current_index < element_index) {
        if (*current_position == ',') {
            current_index++;
        }
        current_position++;
    }

    if (current_index != element_index) {
        // Element index not found in the line, exit
        fprintf(stderr, "Error: Element index %d out of range\n", element_index);
        return;
    }

    // Locate end of the current element
    char *end_position = current_position;
    while (*end_position && *end_position != ',') {
        end_position++;
    }

    // Prepare the new value as a string
    char new_value_str[32];
    snprintf(new_value_str, sizeof(new_value_str), "%.6f", new_value);  // Adjust precision as needed

    // Calculate the length of the new and old values
    size_t old_length = end_position - current_position;
    size_t new_length = strlen(new_value_str);
    size_t length_difference = new_length - old_length;

    // Adjust the string length if needed
    if (length_difference > 0) {
        memmove(end_position + length_difference, end_position, strlen(end_position) + 1);
    } else if (length_difference < 0) {
        memmove(end_position + length_difference, end_position, strlen(end_position) + 1);
    }

    // Copy the new value into the string
    memcpy(current_position, new_value_str, new_length);
}*/

double get_double_in_linedata(char *line, size_t start_offset, int *element_index)
{
  char *current_position = line + start_offset; // Move to starting offset
  int current_index = 0;

  // printf(" here test 0 \n"); // DEBUG

  size_t header_len = strcspn(current_position, ":");

  // printf(" here test 0.1 \n"); // DEBUG
  // printf(" header_len:%d \n", header_len); // DEBUG

  size_t len = strcspn(current_position, "\n");

  // printf(" here test 1 \n"); // DEBUG

  current_position = current_position + header_len + 1;

  // Locate the start of the target element
  // if (*element_index != 0)
  while (*current_position && current_index < *element_index)
  {
    if (*current_position == ',')
    {
      current_index++;
    }
    current_position++;
  }

  // printf(" here test 2 \n"); // DEBUG

  if (current_index != *element_index)
  {
    // Element index not found in the line, exit
    fprintf(stderr, "Error: Element index %d out of range\n", *element_index);
    return;
  }

  // Locate end of the current element
  // char *end_position = current_position;
  // while (*end_position && *end_position != ',')
  // {
  //   end_position++;
  // }

  // printf(" here test 3 \n"); // DEBUG

  char *end_position = current_position + strcspn(current_position, "\n");

  // printf(" here test 4 \n"); // DEBUG

  // Calculate the length of the value
  size_t val_length = end_position - current_position;

  // Create a buffer to store the substring for the value
  char buffer[64]; // Adjust size as needed for large numbers
  if (val_length >= sizeof(buffer))
  {
    fprintf(stderr, "Error: Value is too large to convert\n");
    return 0.0;
  }

  // printf(" here test 5 \n"); // DEBUG

  // Copy the element into the buffer and convert to double
  snprintf(buffer, val_length + 1, "%s", current_position);

  // printf(" val_length + 1:%d \n", val_length + 1);        // DEBUG
  // printf(" current_position:%c \n", current_position[0]); // DEBUG
  // printf(" buffer:%s \n", buffer);                        // DEBUG
  // printf(" here test 6 \n"); // DEBUG

  return atof(buffer);
}

double get_double_in_line(char *line, size_t start_offset, int *element_index)
{
  char *current_position = line + start_offset; // Move to starting offset
  int current_index = 0;

  // printf(" current_position[0]:%c\n", current_position); // DEBUG

  // size_t header_len = strcspn(current_position, ":");

  // // printf(" here test 0.1 \n"); // DEBUG
  // // printf(" header_len:%d \n", header_len); // DEBUG

  // size_t len = strcspn(current_position, "\n");

  // // printf(" here test 1 \n"); // DEBUG
  if (strcspn(current_position, ">") > 0)
    current_position = current_position + strcspn(current_position, ">") + 1;

  // Locate the start of the target element
  // if (*element_index != 0)
  while (*current_position != '\n')
  {
    if (current_index == *element_index)
      break;
    if (*current_position == ',' || *current_position == ':' || *current_position == '>')
    {
      current_index++;
    }
    current_position++;
  }

  // printf(" here test 2 \n"); // DEBUG

  if (current_index != *element_index)
  {
    // Element index not found in the line, exit
    fprintf(stderr, "Error: Element index %d out of range\n", *element_index);
    return;
  }

  // Locate end of the current element
  // char *end_position = current_position;
  // while (*end_position && *end_position != ',')
  // {
  //   end_position++;
  // }

  // printf(" here test 3 \n"); // DEBUG

  char *end_position;

  if (strcspn(current_position, ",") < strcspn(current_position, "\n"))
    end_position = current_position + strcspn(current_position, ",");
  else
    end_position = current_position + strcspn(current_position, "\n");

  // printf(" here test 4 \n"); // DEBUG

  // Calculate the length of the value
  size_t val_length = end_position - current_position;

  // Create a buffer to store the substring for the value
  char buffer[64]; // Adjust size as needed for large numbers
  if (val_length >= sizeof(buffer))
  {
    fprintf(stderr, "Error: Value is too large to convert\n");
    return 0.0;
  }

  // printf(" here test 5 \n"); // DEBUG

  // Copy the element into the buffer and convert to double
  snprintf(buffer, val_length + 1, "%s", current_position);

  // printf(" val_length + 1:%d \n", val_length + 1);        // DEBUG
  // printf(" current_position:%c \n", current_position[0]); // DEBUG
  // printf(" buffer:%s \n", buffer);                        // DEBUG
  // printf(" here test 6 \n"); // DEBUG

  return atof(buffer);
}

ow_value get_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, size_t input_offset)
{
  ow_value ret_val;
  ret_val.x = *x;
  ret_val.y = *y;
  ret_val.z = *z;
  ret_val.value = -1;
  ret_val.line_offset = -1;

  size_t file_offset = -1; // = find_in_file(mapped_file, sb, search_pattern, input_offset) + 1;

  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d:", *x, *y); // Construct search pattern
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

  // if ((int)file_offset < 0)
  // {
  //   // if x,y index cannot be found then we search for y,x index
  //   sprintf(search_pattern, ">%d,%d:", *y, *x); // Construct search pattern
  //   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }

  if ((int)file_offset < 0)
  {
    // printf(" (p) file_offset: %d\n", file_offset); // DEBUG
    if (input_offset > 0)
    {
      // printf(" (p) input_offset: %d\n", input_offset); // DEBUG
      // free(ret_val);
      return get_ow_mmapped(mapped_file, x, y, z, 0);
    }
    return ret_val;
  }

  ret_val.line_offset = file_offset;

  // printf(" get_ow(): x:%d, y:%d, z:%d file_offset:%d \n", *x, *y, *z, file_offset); // DEBUG

#pragma omp critical
  {
    ret_val.value = get_double_in_linedata(mapped_file->mapped_file, file_offset, z);
  }
  // printf("  value:%g\n", ret_val.value); // DEBUG
  return ret_val;
}

int get_double_precision(double value)
{
  char buffer[50];

  // Convert the double to a string with high precision
  // sprintf(buffer, "%f", value);
  // sprintf(buffer, "%.17g", value);
  sprintf(buffer, "%g", value);

  // Find the decimal point in the string
  char *decimal_point = strchr(buffer, '.');
  if (decimal_point == NULL)
  {
    // No decimal point found, so return 0 for integer values
    return 0;
  }

  // printf(" buffer:%s\n", buffer); // DEBUG

  // Count digits after the decimal point
  int decimal_count = 0;
  char *ptr = decimal_point + 1;
  while (*ptr != '\0')
  {
    decimal_count++;
    ptr++;
  }

  return decimal_count;
}

double get_precision_level(double value)
{

  int precision = get_double_precision(value);

  char *buffer = NULL;
  buffer = (char *)malloc((2 + precision) * sizeof(char));
  strcpy(buffer, "");
  strcat(buffer, "0.");
  for (size_t i = 0; i < precision; i++)
  {
    strcat(buffer, "0");
  }
  strcat(buffer, "1");
  strcat(buffer, "\0");
  double ret_val = strtod(buffer, NULL);
  return ret_val;
}

int set_double_in_line(char *line, size_t start_offset, int *element_index, double *new_value)
{
  char *current_position = line + start_offset; // Move to starting offset
  int current_index = 0;

  size_t header_len = strcspn(current_position, ":");
  size_t len = strcspn(current_position, "\n");

  current_position = current_position + header_len + 1;

  // Locate the start of the target element
  // if (*element_index != 0)
  while (*current_position && current_index < *element_index)
  {
    // if (*current_position == '\n')
    // {
    //   break;
    // }
    if (*current_position == ',' || *current_position == ':')
    {
      current_index++;
    }
    current_position++;
  }

  if (current_index != *element_index)
  {
    // Element index not found in the line, exit
    fprintf(stderr, "Error: Element index %d out of range\n", *element_index);
    return;
  }

  // Locate end of the current element
  // char *end_position = current_position;
  // while (*end_position && *end_position != ',')
  // {
  //   end_position++;
  // }

  char *end_position;
  if (strcspn(current_position, ",") < strcspn(current_position, "\n"))
    end_position = current_position + strcspn(current_position, ",") - 1;
  else
    end_position = current_position + strcspn(current_position, "\n") - 1;

  // if (*end_position != '\n')
  // {
  //   end_position--;
  // }

  // printf("value:%g precision:%d\n", *new_value, get_double_precision(*new_value)); // DEBUG
  // Prepare the new value as a string
  int val_precision = get_double_precision(*new_value);
  // val_precision = val_precision <= 1 ? val_precision : val_precision + 1;
  char buffer[val_precision];
  // snprintf(new_value_str, sizeof(new_value_str), "%.6f", *new_value); // Adjust precision as needed
  char *new_value_str = NULL; // dtoa(&buffer, *new_value, val_precision, get_precision_level(*new_value));
  new_value_str = (char *)calloc(1, val_precision * sizeof(char));
  sprintf(new_value_str, "%g", *new_value);
  // Calculate the length of the new and old values
  // size_t old_length = end_position - current_position;
  // size_t new_length = strlen(new_value_str);
  // size_t length_difference = new_length - old_length - 1;
  size_t old_length = end_position - current_position + 1;
  size_t new_length = strlen(new_value_str);
  size_t length_difference = new_length - old_length;

  // printf(" old_len:%d\n", old_length);                  // DEBUG
  // printf(" new_len:%d\n", new_length);                  // DEBUG
  // printf(" length_difference:%d\n", length_difference); // DEBUG

  // Adjust the string length if needed
  if (abs(length_difference) > 0)
  {
    memmove(end_position + length_difference, end_position, strlen(end_position));
    // exit(0); // DEBUG
  }
  // else if (length_difference < 0)
  // {
  //   memmove(end_position + length_difference, end_position, strlen(end_position) + 1);
  // }

  // Copy the new value into the string
  memcpy(current_position, new_value_str, new_length);

  // Calculate the range to msync based on modified bytes
  size_t modified_start = (size_t)(current_position - line);
  size_t modified_end = modified_start + new_length + 1; // +1 to include the comma or end of line
  size_t length_to_sync = modified_end - modified_start;

  // msync to synchronize changes with the file
  // if (msync(line + modified_start, length_to_sync, MS_SYNC) == -1)
  // {
  //   perror("msync failed");
  // }
  if (length_to_sync > 0)
    msync(line + modified_start, length_to_sync, MS_SYNC);
  free(new_value_str);

  // exit(0); // DEBUG

  return length_difference;
}

ow_value set_ow_mmapped(mmapped_file *mapped_file, int *x, int *y, int *z, double *value, size_t input_offset)
{
  // size_t file_offset = 0, line_len = 0;
  // int i = 0, j = 0;
  // char value_char = *value + '0'; // Convert value to character
  char buffer[256];
  char *value_char = NULL;
  int append_mode = FALSE;
  // value_char = dtoa(&buffer, *value, get_double_precision(*value) + 1, get_precision_level(*value));
  value_char = (char *)calloc((get_double_precision(*value) + 1), sizeof(char));
  sprintf(value_char, "%g", *value);
  int len_value_char = strlen(value_char);
  // printf(" get_double_precision(*value): %d\n", get_double_precision(*value)); // DEBUG
  // printf(" get_precision_level(*value): %f\n", get_precision_level(*value));   // DEBUG
  // printf(" value_char: %s\n", value_char);                                     // DEBUG

  ow_value ret_val;

  ret_val.x = *x;
  ret_val.y = *y;
  ret_val.z = *z;
  ret_val.value = -1;
  ret_val.line_offset = -1;

  size_t file_offset = -1;

  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d:", *x, *y); // Construct search pattern
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

  // if ((int)file_offset < 0)
  // {                                             // if we cannot find x,y index then we try to set y,x
  //   sprintf(search_pattern, ">%d,%d:", *y, *x); // Construct search pattern
  //   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }

  if ((int)file_offset < 0)
  {
    if (input_offset > 0)
    {
      return set_ow_mmapped(mapped_file, x, y, z, value, 0);
    }
    // // free(value_char);
    // // return ret_val;
    // append_mode = TRUE;
  }

  //   // if (append_mode == FALSE)
  //   // {
  ret_val.line_offset = file_offset;
  int len_diff = 0;
#pragma omp critical
  {
    len_diff = set_double_in_line(mapped_file->mapped_file, file_offset, z, value);
    if (abs(len_diff) > 0)
    {
      mapped_file = resize_mmap_file(mapped_file, mapped_file->sb.st_size + len_diff);
    }
  }
#pragma omp atomic update
  ret_val.value = *value;
  //   // exit(0); // DEBUG
  //   //   }
  //   //   else
  //   //   {
  //   //     char *append_arr;
  //   //     append_arr = (char *)calloc(sizeof(*x) + sizeof(*y) + sizeof(value_char) + 4, 1);
  //   //     sprintf(append_arr, ">%d,%d:%s\n", *x, *y, value_char);
  //   //     // printf(" append_arr:%s\n", append_arr); // DEBUG
  //   // #pragma omp critical
  //   //     {
  //   //       *mapped_file = append_line_to_mmap(mapped_file, append_arr);
  //   //     }
  //   //     free(append_arr);
  //   //   }
  free(value_char);
  return ret_val;
}

diag_value get_diagval_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, size_t input_offset)
{
  diag_value ret_val;
  ret_val.s0 = *w;
  ret_val.s1 = *x;
  ret_val.hv = *y;
  ret_val.istep = *z;
  ret_val.index = *val_index;
  ret_val.value = -1;
  ret_val.line_offset = -1;
  ret_val.line_num = -1;

  size_t file_offset = -1; // = find_in_file(mapped_file, sb, search_pattern, input_offset) + 1;

  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d,%d,%d:", ret_val.s0, ret_val.s1, ret_val.hv, ret_val.istep); // Construct search pattern
  // printf(" search_pattern:%s:%d\n", search_pattern, *val_index);                               // DEBUG
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

  // if ((int)file_offset < 0)
  // {
  //   // if x,y index cannot be found then we search for y,x index
  //   sprintf(search_pattern, ">%d,%d:", *y, *x); // Construct search pattern
  //   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }

  if ((int)file_offset < 0)
  {
    // printf(" (p) file_offset: %d\n", file_offset); // DEBUG
    if (input_offset > 0)
    {
      // printf(" (p) input_offset: %d\n", input_offset); // DEBUG
      return get_diagval_mmapped(mapped_file, w, x, y, z, val_index, 0);
    }
    return ret_val;
  }

  ret_val.line_offset = file_offset;
  ret_val.line_num = get_line_of_offset(mapped_file, file_offset).line_num;
  // printf(" get_ow(): x:%d, y:%d, z:%d file_offset:%d \n", *x, *y, *z, file_offset); // DEBUG

#pragma omp critical
  {
    ret_val.value = get_double_in_linedata(mapped_file->mapped_file, file_offset, val_index);
  }
  // printf("  value:%g\n", ret_val.value); // DEBUG
  return ret_val;
}

diag_line get_diagvals_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, size_t input_offset)
{
  int anc_idx = 0, diagcnt_idx = 1, b0_idx = 2, b1_idx = 3, ext_idx = 6, sel_idx = 4, trans_idx = 5, weight_idx = 7, ow_idx = 8, sum_idx = 9, cs_idx = 10; // 0-indexed - DATA
  int s0_idx = 1, s1_idx = 2, hv_idx = 3, istep_idx = 4;                                                                                                   // 1-indexed - HEADER

  diag_line ret_val;
  ret_val.s0 = *w;
  ret_val.s1 = *x;
  ret_val.hv = *y;
  ret_val.istep = *z;
  ret_val.line_num = -1;
  ret_val.line_offset = -1;
  size_t file_offset = -1; // = find_in_file(mapped_file, sb, search_pattern, input_offset) + 1;

  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d,%d,%d:", ret_val.s0, ret_val.s1, ret_val.hv, ret_val.istep); // Construct search pattern
  // printf(" search_pattern:%s:%d\n", search_pattern, *val_index);                               // DEBUG
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

  // if ((int)file_offset < 0)
  // {
  //   // if x,y index cannot be found then we search for y,x index
  //   sprintf(search_pattern, ">%d,%d:", *y, *x); // Construct search pattern
  //   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }

  if ((int)file_offset < 0)
  {
    // printf(" (p) file_offset: %d\n", file_offset); // DEBUG
    if (input_offset > 0)
    {
      // printf(" (p) input_offset: %d\n", input_offset); // DEBUG
      return get_diagvals_mmapped(mapped_file, w, x, y, z, val_index, 0);
    }
    return ret_val;
  }

  ret_val.line_offset = file_offset;
  ret_val.line_num = get_line_of_offset(mapped_file, file_offset).line_num;
  // printf(" get_ow(): x:%d, y:%d, z:%d file_offset:%d \n", *x, *y, *z, file_offset); // DEBUG

#pragma omp critical
  {
    ret_val.s0 = get_double_in_line(mapped_file->mapped_file, file_offset, &s0_idx);
    ret_val.s1 = get_double_in_line(mapped_file->mapped_file, file_offset, &s1_idx);
    ret_val.hv = get_double_in_line(mapped_file->mapped_file, file_offset, &hv_idx);
    ret_val.istep = get_double_in_line(mapped_file->mapped_file, file_offset, &istep_idx);

    ret_val.anc_num = (int)get_double_in_linedata(mapped_file->mapped_file, file_offset, &anc_idx);
    ret_val.numsubseq = (int)get_double_in_linedata(mapped_file->mapped_file, file_offset, &diagcnt_idx);
    ret_val.b0 = (int)get_double_in_linedata(mapped_file->mapped_file, file_offset, &b0_idx);
    ret_val.b1 = (int)get_double_in_linedata(mapped_file->mapped_file, file_offset, &b1_idx);
    ret_val.sel = (short)get_double_in_linedata(mapped_file->mapped_file, file_offset, &sel_idx);
    ret_val.trans = (short)get_double_in_linedata(mapped_file->mapped_file, file_offset, &trans_idx);
    ret_val.ext = (int)get_double_in_linedata(mapped_file->mapped_file, file_offset, &ext_idx);
    ret_val.weight = get_double_in_linedata(mapped_file->mapped_file, file_offset, &weight_idx);
    ret_val.ow = get_double_in_linedata(mapped_file->mapped_file, file_offset, &ow_idx);
    ret_val.total_sum = get_double_in_linedata(mapped_file->mapped_file, file_offset, &sum_idx);
    ret_val.cs = (short)get_double_in_linedata(mapped_file->mapped_file, file_offset, &cs_idx);
  }
  // printf("  value:%g\n", ret_val.value); // DEBUG
  return ret_val;
}
// diag_value get_diags_seq_mmapped(mmapped_file *mapped_file, int *x, int *y, int *val_index, size_t input_offset)
// {
//   // int anc_idx = 5, diagcnt_idx = 6, b0_idx = 7, b1_idx = 8, ext_idx = 11, sel_idx = 9, trans_idx = 10, weight_idx = 12, ow_idx = 13, sum_idx = 14, cs_idx = 15, op_idx = 16; // 1-indexed - DATA
//   int s0_idx = 1, s1_idx = 2, hv_idx = 3, istep_idx = 4; // 1-indexed - HEADER
//   diag_value ret_val;
//   ret_val.s0 = *x;
//   ret_val.s1 = *y;
//   ret_val.hv = -1;
//   ret_val.istep = -1;
//   ret_val.index = *val_index;
//   ret_val.value = -1;
//   ret_val.line_offset = -1;
//   ret_val.line_num = -1;

//   size_t file_offset = -1; // = find_in_file(mapped_file, sb, search_pattern, input_offset) + 1;

//   char search_pattern[64];
//   sprintf(search_pattern, ">%d,%d,", ret_val.s0, ret_val.s1); // Construct search pattern
//   // printf(" search_pattern:%s:%d\n", search_pattern, *val_index); // DEBUG
//   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

//   // if ((int)file_offset < 0)
//   // {
//   //   // if x,y index cannot be found then we search for y,x index
//   //   sprintf(search_pattern, ">%d,%d:", *y, *x); // Construct search pattern
//   //   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
//   // }

//   if ((int)file_offset < 0)
//   {
//     // printf(" (p) file_offset: %d\n", file_offset); // DEBUG
//     if (input_offset > 0)
//     {
//       // printf(" (p) input_offset: %d\n", input_offset); // DEBUG
//       return get_diags_seq_mmapped(mapped_file, x, y, val_index, 0);
//     }
//     return ret_val;
//   }

//   ret_val.line_offset = file_offset;
//   ret_val.line_num = get_line_of_offset(mapped_file, file_offset).line_num;
//   // printf(" get_ow(): x:%d, y:%d, z:%d file_offset:%d \n", *x, *y, *z, file_offset); // DEBUG

// #pragma omp critical
//   {
//     ret_val.value = get_double_in_linedata(mapped_file->mapped_file, file_offset, val_index);
//     ret_val.hv = get_double_in_line(mapped_file->mapped_file, file_offset, &hv_idx);
//     ret_val.istep = get_double_in_line(mapped_file->mapped_file, file_offset, &istep_idx);
//   }
//   // printf("  value:%g\n", ret_val.value); // DEBUG
//   return ret_val;
// }

diag_value get_diagval_line_mmapped(mmapped_file *mapped_file, int *line_num, int *val_index, size_t input_offset)
{
  int current_line = 0;
  size_t offset = 0;
  diag_value ret_val;
  int s0_idx = 1, s1_idx = 2, hv_idx = 3, istep_idx = 4; // 1-indexed - HEADER

  ret_val.s0 = -1;
  ret_val.s1 = -1;
  ret_val.hv = -1;
  ret_val.istep = -1;
  ret_val.line_num = -1;
  ret_val.line_offset = 0;
  ret_val.index = *val_index;

  // printf(" line_num:%d\n", *line_num); // DEBUG
  if (input_offset > 0)
  {
    line_value line_val = get_line_of_offset(mapped_file, input_offset);
    current_line = line_val.line_num;
    offset = line_val.line_offset;
    // exit(0); // DEBUG
  }
  else
  {
    for (size_t tmp_offset = 0; tmp_offset < mapped_file->sb.st_size; tmp_offset++)
    {
      if (mapped_file->mapped_file[tmp_offset] == '>')
      {
        offset = tmp_offset;
        current_line++;
      }
      if (current_line == *line_num)
        break;
    }
  }

  // printf(" current_line:%d offset=%u\n", current_line, offset); // DEBUG
  if ((offset >= mapped_file->sb.st_size || current_line != *line_num) && input_offset != 0)
  {
    return get_diagval_line_mmapped(mapped_file, line_num, val_index, 0);
  }

  ret_val.line_num = current_line;
  ret_val.line_offset = offset;

  // char *line_map = &(mapped_file->mapped_file[offset]);
  // // offset = strcspn(line_map, ">") + 1;
  // int line_len = strcspn(line_map, "\n");
  // char *line = strndup(line_map, line_len);
#pragma omp critical
  {
    ret_val.s0 = get_double_in_line(mapped_file->mapped_file, offset, &s0_idx);
    ret_val.s1 = get_double_in_line(mapped_file->mapped_file, offset, &s1_idx);
    ret_val.hv = get_double_in_line(mapped_file->mapped_file, offset, &hv_idx);
    ret_val.istep = get_double_in_line(mapped_file->mapped_file, offset, &istep_idx);
    ret_val.value = get_double_in_linedata(mapped_file->mapped_file, offset, val_index);

    // printf(" ret_val.s0:%d ret_val.s1:%d\n", ret_val.s0, ret_val.s1); // DEBUG
  }
  // free(line);
  // printf(" curr_line:%d line_num:%d input_offset:%d offset:%u ret_val s0:%d s1:%d hv:%d value:%g\n", current_line, *line_num, input_offset, offset, ret_val.s0, ret_val.s1, ret_val.hv, ret_val.value); // DEBUG

  return ret_val;
}

diag_line get_diagvals_line_mmapped(mmapped_file *mapped_file, int *line_num, size_t input_offset)
{
  int current_line = 0;
  size_t offset = 0;
  diag_line ret_val;
  int anc_idx = 0, diagcnt_idx = 1, b0_idx = 2, b1_idx = 3, ext_idx = 6, sel_idx = 4, trans_idx = 5, weight_idx = 7, ow_idx = 8, sum_idx = 9, cs_idx = 10; // 0-indexed - DATA
  int s0_idx = 1, s1_idx = 2, hv_idx = 3, istep_idx = 4;                                                                                                   // 1-indexed - HEADER

  ret_val.s0 = -1;
  ret_val.s1 = -1;
  ret_val.hv = -1;
  ret_val.istep = -1;
  ret_val.line_num = -1;
  ret_val.line_offset = -1;
  // ret_val.added_weights = -1;

  if (input_offset > 0)
  {
    line_value line_val = get_line_of_offset(mapped_file, input_offset);
    current_line = line_val.line_num;
    offset = line_val.line_offset;
  }
  else
  {
    for (size_t tmp_offset = 0; tmp_offset < mapped_file->sb.st_size; tmp_offset++)
    {
      if (mapped_file->mapped_file[tmp_offset] == '>')
      {
        offset = tmp_offset;
        current_line++;
      }
      if (current_line == *line_num)
        break;
    }
  }

  if ((offset >= mapped_file->sb.st_size || current_line != *line_num) && input_offset > 0)
  {
    return get_diagvals_line_mmapped(mapped_file, line_num, 0);
  }

  // printf("get_line_of_offset():%u offset:%u\n", get_line_of_offset(mapped_file, input_offset).line_offset, offset);                                                                                        // DEBUG
  // printf("mapped_file->mapped_file[val]:%c mapped_file->mapped_file[offset]:%c\n", mapped_file->mapped_file[get_line_of_offset(mapped_file, input_offset).line_offset], mapped_file->mapped_file[offset]); // DEBUG

  ret_val.line_num = current_line;
  ret_val.line_offset = offset;

  // printf("first char:%c\n", mapped_file->mapped_file[offset]); // DEBUG
  // if (mapped_file->mapped_file[offset] != '>')                 // DEBUG
  //   exit(0);                                                   // DEBUG
  // char *line_map = &(mapped_file->mapped_file[offset]);
  // int line_len = strcspn(line_map, "\n");
  // char *line = strndup(line_map, line_len);

  // char *line_map = &(mapped_file->mapped_file[offset]);
  // // offset = strcspn(line_map, ">") + 1;
  // int line_len = strcspn(line_map, "\n");
  // char *line = strndup(line_map, line_len);
#pragma omp critical
  {
    ret_val.s0 = get_double_in_line(mapped_file->mapped_file, offset, &s0_idx);
    ret_val.s1 = get_double_in_line(mapped_file->mapped_file, offset, &s1_idx);
    ret_val.hv = get_double_in_line(mapped_file->mapped_file, offset, &hv_idx);
    ret_val.istep = get_double_in_line(mapped_file->mapped_file, offset, &istep_idx);

    ret_val.anc_num = (int)get_double_in_linedata(mapped_file->mapped_file, offset, &anc_idx);
    ret_val.numsubseq = (int)get_double_in_linedata(mapped_file->mapped_file, offset, &diagcnt_idx);
    ret_val.b0 = (int)get_double_in_linedata(mapped_file->mapped_file, offset, &b0_idx);
    ret_val.b1 = (int)get_double_in_linedata(mapped_file->mapped_file, offset, &b1_idx);
    ret_val.sel = (short)get_double_in_linedata(mapped_file->mapped_file, offset, &sel_idx);
    ret_val.trans = (short)get_double_in_linedata(mapped_file->mapped_file, offset, &trans_idx);
    ret_val.ext = (int)get_double_in_linedata(mapped_file->mapped_file, offset, &ext_idx);
    ret_val.weight = get_double_in_linedata(mapped_file->mapped_file, offset, &weight_idx);
    ret_val.ow = get_double_in_linedata(mapped_file->mapped_file, offset, &ow_idx);
    ret_val.total_sum = get_double_in_linedata(mapped_file->mapped_file, offset, &sum_idx);
    ret_val.cs = (short)get_double_in_linedata(mapped_file->mapped_file, offset, &cs_idx);
    // ret_val.added_weights = (double)get_double_in_linedata(mapped_file->mapped_file, offset, &added_weights_idx);

    // sscanf(line, ">%d,%d,%d,%d:%d,%d,%d,%d,%u,%u,%d,%g,%g,%g", &ret_val.s0, &ret_val.s1, &ret_val.hv, &ret_val.istep, &ret_val.anc_num, &ret_val.numsubseq, &ret_val.b0, &ret_val.b1, &ret_val.sel, &ret_val.trans, &ret_val.ext, &ret_val.weight, &ret_val.ow, &ret_val.total_sum);
    // printf("(FULL) ret_val.s0:%d ret_val.s1:%d\n", ret_val.s0, ret_val.s1); // DEBUG
  }
  // free(line);

  return ret_val;
}

diag_value set_diags_mmapped(mmapped_file *mapped_file, int *w, int *x, int *y, int *z, int *val_index, double *value, size_t input_offset)
{
  // size_t file_offset = 0, line_len = 0;
  // int i = 0, j = 0;
  // char value_char = *value + '0'; // Convert value to character
  char buffer[256];
  char *value_char = NULL;
  int append_mode = FALSE;
  // value_char = dtoa(&buffer, *value, get_double_precision(*value) + 1, get_precision_level(*value));
  value_char = (char *)calloc((get_double_precision(*value) + 1), sizeof(char));
  sprintf(value_char, "%g", *value);
  int len_value_char = strlen(value_char);
  // printf(" get_double_precision(*value): %d\n", get_double_precision(*value)); // DEBUG
  // printf(" get_precision_level(*value): %f\n", get_precision_level(*value));   // DEBUG
  // printf(" value_char: %s\n", value_char);                                     // DEBUG

  diag_value ret_val;
  ret_val.s0 = *w;
  ret_val.s1 = *x;
  ret_val.hv = *y;
  ret_val.istep = *z;
  ret_val.index = *val_index;
  ret_val.value = -1;
  ret_val.line_offset = -1;

  size_t file_offset = -1;

  char search_pattern[64];
  sprintf(search_pattern, ">%d,%d,%d,%d:", ret_val.s0, ret_val.s1, ret_val.hv, ret_val.istep); // Construct search pattern
  file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

  // if ((int)file_offset < 0)
  // {                                             // if we cannot find x,y index then we try to set y,x
  //   sprintf(search_pattern, ">%d,%d:", *y, *x); // Construct search pattern
  //   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);
  // }

  if ((int)file_offset < 0)
  {
    if (input_offset > 0)
    {
      return set_diags_mmapped(mapped_file, w, x, y, z, val_index, value, 0);
    }
    // // free(value_char);
    // // return ret_val;
    // append_mode = TRUE;
  }

  // if (append_mode == FALSE)
  // {
  ret_val.line_offset = file_offset;
  int len_diff = 0;
#pragma omp critical
  {
    len_diff = set_double_in_line(mapped_file->mapped_file, file_offset, val_index, value);
    // printf("len_diff:%d\n", len_diff); // DEBUG
    // exit(0);                           // DEBUG
    if (abs(len_diff) > 0)
    {
      mapped_file = resize_mmap_file(mapped_file, mapped_file->sb.st_size + len_diff);
    }
  }
#pragma omp atomic update
  ret_val.value = *value;
  // exit(0); // DEBUG
  //   }
  //   else
  //   {
  //     char *append_arr;
  //     append_arr = (char *)calloc(sizeof(*x) + sizeof(*y) + sizeof(value_char) + 4, 1);
  //     sprintf(append_arr, ">%d,%d:%s\n", *x, *y, value_char);
  //     // printf(" append_arr:%s\n", append_arr); // DEBUG
  // #pragma omp critical
  //     {
  //       *mapped_file = append_line_to_mmap(mapped_file, append_arr);
  //     }
  //     free(append_arr);
  //   }
  free(value_char);
  // exit(0); // DEBUG
  return ret_val;
}

// diag_value set_diags_sel_mmapped(mmapped_file *mapped_file, int *x, int *y, int *val_index, double *value, size_t input_offset){
//   // size_t file_offset = 0, line_len = 0;
//   // int i = 0, j = 0;
//   // char value_char = *value + '0'; // Convert value to character
//   char buffer[256];
//   char *value_char = NULL;
//   int append_mode = FALSE;
//   // value_char = dtoa(&buffer, *value, get_double_precision(*value) + 1, get_precision_level(*value));
//   value_char = (char *)calloc((get_double_precision(*value) + 1), sizeof(char));
//   sprintf(value_char, "%g", *value);
//   int len_value_char = strlen(value_char);
//   // printf(" get_double_precision(*value): %d\n", get_double_precision(*value)); // DEBUG
//   // printf(" get_precision_level(*value): %f\n", get_precision_level(*value));   // DEBUG
//   // printf(" value_char: %s\n", value_char);                                     // DEBUG

//   diag_value ret_val;
//   ret_val.s0 = *x;
//   ret_val.s1 = *y;
//   ret_val.hv = -1;
//   ret_val.istep = -1;
//   ret_val.index = *val_index;
//   ret_val.value = -1;
//   ret_val.line_offset = -1;

//   size_t file_offset = -1;

//   char search_pattern[64];
//   sprintf(search_pattern, ">%d,%d,", ret_val.s0, ret_val.s1); // Construct search pattern
//   file_offset = find_in_file(mapped_file->mapped_file, &(mapped_file->sb), search_pattern, input_offset);

//   if ((int)file_offset < 0)
//   {
//     if (input_offset > 0)
//     {
//       return set_diags_sel_mmapped(mapped_file, x, y, val_index, value, 0);
//     }
//   }

//   ret_val.line_offset = file_offset;

//   ret_val.hv = get_double_in_line();
//   ret_val.istep = -1;

//   int len_diff = 0;

// #pragma omp critical
//   {
//     len_diff = set_double_in_line(mapped_file->mapped_file, file_offset, val_index, value);
//     // printf("len_diff:%d\n", len_diff); // DEBUG
//     // exit(0);                           // DEBUG
//     if (abs(len_diff) > 0)
//     {
//       mapped_file = resize_mmap_file(mapped_file, mapped_file->sb.st_size + len_diff);
//     }
//   }

// #pragma omp atomic update
//   ret_val.value = *value;

//   free(value_char);
//   return ret_val;
// }

// // Function to read sequence file using mmap
// int seq_read_mmap(const char *seq_file, char **sq, char **sqn, char **fsqn)
// {
//   // int fd;
//   // struct stat sb;
//   // char *mapped;
//   // char *line;
//   // size_t len;
//   int i, j, k, pv, crc;
//   int max_char[MAX_SEQNUM];

//   // // Open the file
//   // fd = open(seq_file, O_RDONLY);
//   // if (fd == -1)
//   // {
//   //   perror("Error opening file");
//   //   return -1;
//   // }

//   // // Get the size of the file
//   // if (fstat(fd, &sb) == -1)
//   // {
//   //   perror("Error getting file size");
//   //   close(fd);
//   //   return -1;
//   // }

//   // // Map the file into memory
//   // mapped = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
//   // if (mapped == MAP_FAILED)
//   // {
//   //   perror("Error mapping file");
//   //   close(fd);
//   //   return -1;
//   // }

//   // // Close the file descriptor, as it is no longer needed after mmap
//   // close(fd);

//   mmapped_file *mmapped_seq = mmap_file(seq_file, O_RDONLY, PROT_READ, MAP_SHARED, 0);
//   struct stat *sb = &(mmapped_seq->sb);
//   size_t file_size = sb->st_size;

//   // Start reading the sequences from the mapped file
//   size_t offset = 0;
//   int sn = -1;

//   while (offset < file_size)
//   {
//     char *line_map = &(mmapped_seq->mapped_file[offset]);

//     // Find the length of the current line
//     size_t len = strcspn(line_map, "\n");
//     char *line = strndup(line_map, len);
//     // if (line[0] == '>') {
//     //     // Sequence header line
//     //     sn++;
//     //     fsqn[sn] = strndup(line, len);
//     //     sqn[sn] = strndup(line + 1, len - 1); // Skip '>'
//     // } else {
//     //     // Sequence data line
//     //     printf("Processing sequence %d\n", sn);

//     //     sq[sn] = realloc(sq[sn], (sq[sn] ? strlen(sq[sn]) : 0) + len + 1);
//     //     strncat(sq[sn], line, len);
//     // }

//     ws_remove(line);

//     if (line[0] == '>')
//     {
//       sn++;

//       // if(sn>=20314){
//       //   exit(1);
//       // }

//       printf("Processing sequence %d\n", sn);

//       if (sn >= MAX_SEQNUM - 1)
//       {
//         printf("Error: Number of sequences exceeds the maximum allowed (%d)\n", MAX_SEQNUM);
//         exit(1);
//       }

//       n_clean(line);

//       // fsqn[ sn ] = ( char * ) calloc( strlen( line ) + 3 , sizeof ( char ) );
//       fsqn[sn] = (char *)malloc((strlen(line) + 3) * sizeof(char));

//       // printf("\n here 3.5 \n"); // DEBUG

//       strcpy(fsqn[sn], line);

//       // printf("\n here 3.6 \n"); // DEBUG

//       max_char[sn] = 0;
//       // sqn[ sn ]  = ( char * ) calloc( SEQ_NAME_LEN + 3 , sizeof ( char ) );
//       sqn[sn] = (char *)malloc((SEQ_NAME_LEN + 3) * sizeof(char));
//       if (sqn[sn] == NULL)
//       {
//         printf("Error: Memory allocation failed for sequence %d\n", sn);
//         exit(1);
//       }

//       // printf("\n here 3.7 \n"); // DEBUG

//       for (crc = 0; crc < SEQ_NAME_LEN; crc++)
//         if (crc < strlen(line))
//           sqn[sn][crc] = line[crc];
//         else
//           sqn[sn][crc] = ' ';

//       // printf("\n here 3.8 \n"); // DEBUG

//       sqn[sn][SEQ_NAME_LEN] = '\0';

//       // printf("\n here 3.9 \n"); // DEBUG
//     }
//     else
//     {
//       // if(strcmp(line,' ') != 0 || strcmp(line,"\n") != 0){
//       // printf("\n here 3.10 : %s \n", line); // DEBUG
//       max_char[sn] = max_char[sn] + strlen(line) - 1;
//       // printf("\n here 3.11 \n"); // DEBUG
//       // }
//     }

//     // Move to the next line
//     offset += len + 1;

//     free(line);
//   }

//   // printf("\n here 4 \n"); // DEBUG

//   for (i = 0; i <= sn; i++)
//   {
//     // sq[ i ]  = ( char * ) calloc( max_char[ i ] + 1 , sizeof ( char ) );
//     sq[i] = (char *)malloc((max_char[i] + 1) * sizeof(char));
//     // printf("\n here 4.1 \n"); // DEBUG
//   }

//   // printf("\n here 4.2 \n"); // DEBUG

//   // if( (seqlen = (int *) calloc( ( sn + 1 ) , sizeof(int) )) == NULL)
//   if ((seqlen = (int *)malloc((sn + 1) * sizeof(int))) == NULL)
//   {
//     // printf("\n here 4.3 \n"); // DEBUG
//     erreur("\n\n problems with memory allocation for `seqlen' \n\n");
//   }

//   // printf("\n here 4.4 \n"); // DEBUG

//   // printf("\n here 4.5 \n"); // DEBUG

//   /******************************************/

//   if (self_comparison == 1)
//   {
//     if (sn != 0)
//     {
//       printf("\n\n With option \"self comparison\" input file must contain one single sequence \n\n");
//       exit(1);
//     }

//     // sq[ 1 ]  = ( char * ) calloc( max_char[ 0 ] + 1 , sizeof ( char ) );
//     sq[1] = (char *)malloc((max_char[0] + 1) * sizeof(char));

//     // sqn[ 1 ]  = ( char * ) calloc( strlen( line ) + 3 , sizeof ( char ) );
//     sqn[1] = (char *)malloc((strlen(line) + 3) * sizeof(char));
//     strcpy(sqn[1], sqn[0]);
//   }

//   /******************************************/

//   offset = 0;
//   sn = -1;
//   while (offset < file_size)
//   {
//     // line = &mapped[offset];
//     char *line_map = &(mmapped_seq->mapped_file[offset]);

//     // Find the length of the current line
//     size_t len = strcspn(line_map, "\n");
//     char *line = strndup(line_map, len);

//     // printf("\n here 5.1 : %s \n", line); // DEBUG
//     ws_remove(line);
//     if (line[0] == '>')
//     {
//       sn++;
//       j = 0;
//     }
//     else
//       for (k = 0; k < strlen(line); k++)
//         if (
//             (line[k] >= 65) && (line[k] <= 90) ||
//             (line[k] >= 97) && (line[k] <= 122))
//           sq[sn][j++] = toupper(line[k]);

//     // printf("\n here 5.2 \n"); //DEBUG
//     offset += len + 1;
//   }

//   sn++;

//   // printf("\n here 6 \n"); // DEBUG

//   for (i = 0; i < sn; i++)
//   {
//     seqlen[i] = strlen(sq[i]);
//   }

//   // printf("\n here 7 \n"); // DEBUG

//   if (self_comparison)
//   {
//     seqlen[1] = seqlen[0];
//     for (i = 0; i <= seqlen[0]; i++)
//       sq[1][i] = sq[0][i];
//     sn++;
//   }

//   // printf("\n here 8 \n"); // DEBUG

//   // Unmap the file after reading
//   munmap(mmapped_seq->mapped_file, file_size);
//   // free(sb);
//   free(mmapped_seq);

//   return sn; // + 1; // Return number of sequences read
// }

void matrix_read(FILE *fp_mat)
{
  int i, j;
  char line[MLINE], dummy[MLINE];

  fgets(line, MLINE, fp_mat);
  fgets(line, MLINE, fp_mat);

  for (i = 1; i <= 20; i++)
  {
    for (j = i; j <= 20; j++)
    {
      fscanf(fp_mat, "%d", &sim_score[i][j]);
      sim_score[j][i] = sim_score[i][j];
      if (sim_score[i][j] > max_sim_score)
        max_sim_score = sim_score[i][j];
    }

    fscanf(fp_mat, "%s\n", dummy);
  }

  fclose(fp_mat);

  for (i = 0; i <= 20; i++)
  {
    sim_score[i][0] = 0;
    sim_score[0][i] = 0;
  }

  /*
   sim_score[0][0] = max_sim_score ;
  */
}

void tp400_read(int w_type, double **pr_ptr)
{

  /* reads probabilities from file */
  /* w_type = 0 (protein), 1 (dna w/o transl.), 2 (dna with transl.) */

  char line[MLINE], file_name[MLINE], suffix[10], str[MLINE];
  int sum, len, max_sim, i;
  double pr;

  FILE *fp;

  if (w_type == 0)
  {
    strcpy(suffix, "prot");
  }

  if (w_type == 1)
  {
    strcpy(suffix, "dna");
  }

  if (w_type == 2)
  {
    strcpy(suffix, "trans");
  }

  strcpy(file_name, par_dir);
  strcat(file_name, "/tp400_");
  strcat(file_name, suffix);

  if ((fp = fopen(file_name, "r")) == NULL)
  {
    printf("\n\n Cannot find the file %s \n\n", file_name);
    printf(" Make sure the environment variable DIALIGN2_DIR points\n");
    printf(" to a directory containing the files \n\n");
    printf("   BLOSUM \n   tp400_dna\n   tp400_prot \n   tp400_trans \n\n");
    printf(" These files should be contained in the DIALIGN package \n\n\n");
    exit(1);
  }

  if (fgets(line, MLINE, fp) == NULL)
    erreur("\n\n problem with file %s  \n\n", file_name);
  else if (w_type % 2)
    av_sim_score_nuc = atof(line);
  else
    av_sim_score_pep = atof(line);

  while (fgets(line, MLINE, fp) != NULL)
  {
    sscanf(line, "%d %d %s", &len, &sum, str);

    pr = atof(str);
    pr_ptr[len][sum] = pr;
  }

  fclose(fp);

} /*  tp400_read  */
