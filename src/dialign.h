
#define PAPER_WIDTH      80 
#define MLINE 1000 
#define MAX_REGEX 1000 
#define NAME_LEN 1000   
#ifndef SEQ_NAME_LEN
#define SEQ_NAME_LEN 12   /* width of the name column in the alignment output */
#endif
#define MAX_ITNUM 3 

#define MIN_MOT_WGT 0.1 
#define MAX_CSC 10 


         /**************************\
         *                          * 
         *    default parameters    *
         *                          * 
         \**************************/


#define BETA              0
#define WEB               0
#define OVERLAP_THRESHOLD 35
#define MIN_DIA           1
#define MAX_DIA          40
#define MATNAME          "BLOSUM"     
#define WEAK_WGT_TYPE_THR      0.5 
#define STRONG_WGT_TYPE_THR    0.75 


struct pair_frag {int b1, b2, ext; float weight; short trans, cs; 
                   struct pair_frag *prec, *last; float sum; };
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

/* A diagonal ("fragment").  36 bytes (the original: 56).  The list link is a
   32-bit index into the fragment arena (frags.c) instead of two 64-bit
   pointers; sel/trans/cs are bit fields, `it' 16 bit.  Node 0 is reserved,
   an index of 0 means "no next node".  Use NX(p) / SETNX(p,q) for the link. */
struct multi_frag {int b[2], s[2], ext; float weight, ow; unsigned int next;
                   unsigned short it; unsigned char sel:1, trans:1, cs:1;};
     /*
           fragments outside function `pairalign' 

           b[0], b[1]:  begin of the diagonal
           s[0], s[1]:  sequences, to which diagonal belongs
           ext:         length of the diagonal
           weight:      individual weight of the diagonal
           ow:          overlap weight of the diagonal
           sel:         1, if accepted in filter proces, 0 else
           trans:       translation
           cs:          crick strand 
           it:          iteration step 
           next:        index of the next diagonal (0 = end of list) 
     */

struct leaf {int s1, s2, clade;};
struct seq_pair {int s1, s2; float weight;};      

struct subtree { int member_num, valid ; int *member; char *name ;
                 float depth; };         




/* flags "position p of sequence i is already directly aligned with a residue
   of sequence j" (original: open_pos[i][j][p], stored inverted, see dialign.c) */
extern unsigned char **closed_pos_row;
#define CLOSED_POS(i,j) ( closed_pos_row[(i)] + (size_t)(j) * ( (size_t) seqlen[(i)] + 2 ) )


/* ---- fragment arena (frags.c) ---- */
#include <stdint.h>
extern struct multi_frag *frag_base;                 /* node i lives at frag_base + i */
struct multi_frag *frag_new( void );                 /* zero-filled list node        */
void   frag_free( struct multi_frag *p );
static inline unsigned int frag_id_of( const struct multi_frag *p )
  { return p ? (unsigned int)( p - frag_base ) : 0u ; }
#define FRAG_ID(p)   frag_id_of( p )
#define NX(p)        ( (p)->next ? frag_base + (p)->next : (struct multi_frag *) 0 )
#define SETNX(p,q)   ( (p)->next = frag_id_of( q ) )
