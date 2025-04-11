
#pragma once

#define PAPER_WIDTH 80
#define MLINE 1000
#define MAX_REGEX 1000
#define NAME_LEN 1000
#define SEQ_NAME_LEN 50
#define MAX_SEQNUM 50000
#define MAX_ITNUM 3
#define MAX_INPUT_LINE 10000
#define MIN_MOT_WGT 0.1
#define MAX_CSC 10

/**************************\
*                          *
*    default parameters    *
*                          *
\**************************/

#define BETA 0
#define WEB 0
#define OVERLAP_THRESHOLD 35
#define MIN_DIA 1
#define MAX_DIA 40
#define MATNAME "BLOSUM"
#define WEAK_WGT_TYPE_THR 0.5
#define STRONG_WGT_TYPE_THR 0.75

struct pair_frag
{
      int b1, b2, ext;
      double weight;
      short trans, cs;
      struct pair_frag *prec, *last;
      double sum;
};
/*
      fragments in function `pairalign'

      b1, b2:    begin of the diagonal
      ext:       length of the diagonal
      weight:    weight of the diagonal
      prec:      preceding diagonal in dot matrix
      last:      last diagonal ending in the same column
      sum:       sum of weights accumulated
      cs:        crick strand
      trans:     translation
*/

struct multi_frag
{
      int b[2], s[2], ext, it, hv, diag_count;
      double weight;
      double ow;
      short sel, trans;
      short cs;
      struct multi_frag *next, *pred;
      size_t line_offset;
      size_t prev_offset;
      size_t fasta_offset;
      size_t ow_offset;
      size_t trans_offset;
      int anc_num;
      int line_num;
};
/*
      fragments outside function `pairalign'

      b[0], b[1]:  begin of the diagonal
      s[0], s[1]:  sequences, to which diagonal belongs
      hv: hv
      ext:         length of the diagonal
      weight:      individual weight of the diagonal
      ow:          overlap weight of the diagonal
      sel:         1, if accepted in filter proces, 0 else
      trans:       translation
      cs:          crick strand
      it:          iteration step
      *next:       next diagonal
*/

// struct multi_frag_min
// {
//       int idx, s0, s1, hv, istep;
//       double total_sum;
//       int line_num;
//       size_t line_offset;
//       struct multi_frag_min *next;
//       struct multi_frag_min *prev;
// };
// /*
//       fragments outside function `pairalign'

//       idx: Index of current diagonal
//       s0, s1, hv, istep : s[0], s[1], hv, iteration step and sum of fragment weights/ow
//       line_num: line number of next diagonal
//       *next:       next diagonal
//       *prev: previous diagonal
// */

struct leaf
{
      int s1, s2, clade;
};
struct seq_pair
{
      int s1, s2;
      double weight;
};

struct subtree
{
      int member_num, valid;
      int *member;
      char *name;
      float depth;
};

#include <sys/stat.h>

typedef struct
{
      char *mapped_file;
      struct stat sb;
      // char *eof;
      // int fd;
      char *file_name;
      int file_mode;
      int protocol;
      int map_mode;
} mmapped_file;

typedef struct
{
      int x, y, z;
      int value;
      size_t line_offset;
} openpos_value;

typedef struct
{
      int x, y, z;
      double value;
      size_t line_offset;
} ow_value;

typedef struct
{
      int i, j, k, l;
      int value;
      size_t line_offset;
} subst_mat_value;

typedef struct
{
      int seq_num;
      char *data;
      size_t line_offset;
} fasta_value;

typedef struct
{
      int seq_num;
      size_t len;
      size_t line_offset;
} fasta_len_value;

typedef struct
{
      double value;
      size_t line_offset;
} frag_chain_value;

typedef struct
{
      size_t line_offset;
      int line_num;
} line_value;

typedef struct
{
      int s0, s1, hv, istep, index, line_num;
      double value;
      size_t line_offset;
} diag_value;

typedef struct
{
      int s0, s1, hv, istep, line_num, anc_num, numsubseq, b0, b1, ext;
      short sel, trans, cs;
      double weight, ow, total_sum;
      size_t line_offset;
} diag_line;