/*
 *  frags.c -- arena for the nodes of the diagonal list.
 *
 *  All list nodes live in ONE contiguous, file-backed region so that a node
 *  can be named by a 32-bit index.  Address space for 2^32 nodes is reserved
 *  up front (PROT_NONE, costs nothing); storage is attached in 64 MiB
 *  pieces as the list grows, so memory/disk use follows the real number of
 *  diagonals.  Only the main thread creates or frees nodes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "dialign.h"
#include "mmstore.h"

struct multi_frag *frag_base = NULL;
static size_t   max_slots = 0;            /* usable node indices (1 .. max_slots-1) */
static size_t   used = 1;                 /* next unused index (0 is reserved)      */
static size_t   mapped = 0;               /* indices backed so far                  */
static unsigned free_head = 0;            /* freed nodes, linked through ->next     */
static int      debug = -1;

#define CHUNK_BYTES ((size_t)64 << 20)

static void frag_init( void )
{
#if SIZE_MAX > 0xFFFFFFFFu
    /* Native / any platform with a real reserve-now, commit-later virtual
       memory subsystem and a 64-bit size_t (so `(size_t)1 << 32` below is
       well-defined, not a shift equal to the type's own bit width): try
       for 2^32 nodes of PROT_NONE address space up front (free until
       touched), halving down to a floor if that's refused. */
    size_t       slots = (size_t)1 << 32 ;
    const size_t floor = (size_t)1 << 24 ;
    void *p = NULL ;
    while( slots >= floor ) {
        p = mm_reserve( slots * sizeof(struct multi_frag) );
        if( p ) break ;
        slots >>= 1 ;
    }
    if( !p ) { fprintf(stderr,"dialign: cannot reserve address space for fragments\n"); exit(1); }
    frag_base = (struct multi_frag *) p ;
    max_slots = slots ;
    if( max_slots > 0xFFFFFFFFu ) max_slots = 0xFFFFFFFFu ;
#else
    /* 32-bit size_t (e.g. WASM/Emscripten): there is no real virtual-memory
       subsystem here to reserve-now/commit-later into. Observed directly:
       mm_reserve()'s PROT_NONE mmap can succeed (Emscripten's shim appears
       to treat it as a no-op), but the later mm_map_fixed() that is
       supposed to back part of that same reservation with real memory via
       a second, MAP_FIXED mmap at that exact address does not -- "cannot
       map memory for fragments" on the very first fragment allocation.
       That two-step reserve/commit-at-a-fixed-address pattern has no real
       equivalent in a single growable linear-memory model, and there is
       far less address space to begin with regardless. So on this class
       of platform the arena is just one ordinary, single allocation sized
       to exactly one chunk (CHUNK_BYTES -- the same per-chunk size the
       native path grows by, ~64 MiB) -- plenty for a browser-side job, but
       it does not grow beyond that: by setting mapped == max_slots == per
       up front, a job needing more diagonals than this takes the existing
       "too many diagonals" exit below (see frag_new()) instead of ever
       calling the broken mm_map_fixed() path again. */
    size_t per = ( CHUNK_BYTES / sizeof(struct multi_frag) ) & ~(size_t)16383 ;
    void *p = calloc( per, sizeof(struct multi_frag) );
    if( !p ) { fprintf(stderr,"dialign: cannot allocate memory for fragments\n"); exit(1); }
    frag_base = (struct multi_frag *) p ;
    max_slots = per ;
    mapped    = per ;
#endif
    debug = ( getenv("DIALIGN_MMDEBUG") != NULL );
}

struct multi_frag *frag_new( void )
{
    struct multi_frag *p ;
    if( !frag_base ) frag_init();
    if( free_head ) {
        p = frag_base + free_head ;
        free_head = p->next ;
        memset( p , 0 , sizeof *p );
        return p ;
    }
    if( used >= mapped ) {
        /* per*sizeof must be a multiple of any page size (up to 64 KiB) */
        size_t per = ( CHUNK_BYTES / sizeof(struct multi_frag) ) & ~(size_t)16383 ;
        if( mapped + per > max_slots ) {
            fprintf(stderr,"dialign: too many diagonals (more than %lu)\n",(unsigned long)max_slots);
            exit(1);
        }
        if( mm_map_fixed( (char *) frag_base + mapped * sizeof(struct multi_frag) ,
                          per * sizeof(struct multi_frag) ) == NULL ) {
            fprintf(stderr,"dialign: cannot map memory for fragments\n");
            exit(1);
        }
        mapped += per ;
    }
    return frag_base + ( used++ );
}

void frag_free( struct multi_frag *p )
{
    if( debug ) memset( p , 0xDD , sizeof *p );
    p->next = free_head ;
    free_head = FRAG_ID(p) ;
}
