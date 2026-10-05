/*
 *  mmstore.h  --  file-backed ("mmapped") memory store for dialign-mmap
 *
 *  Every data structure whose size grows with the input (sequences, the
 *  open-position tables, the closure/frontier tables, the fragment lists,
 *  the output buffers ...) is allocated from this store instead of from the
 *  process heap.  The store hands out ordinary, stable pointers, so the
 *  algorithmic code is byte-for-byte the code of the original DIALIGN
 *  2.2.1; but the bytes behind the pointers live in (unlinked, sparse)
 *  temporary files that are mmap()ed MAP_SHARED.  The kernel is therefore
 *  free to evict cold pages to disk, and the resident *anonymous* memory
 *  of the process stays small no matter how large the input gets.
 *
 *  Design
 *  ------
 *   - small blocks (< MM_LARGE bytes) are served from per-thread free
 *     lists that are refilled from 64 MiB file-backed chunks
 *     (size-class allocator, 16 byte alignment, no lock on the fast path);
 *   - large blocks get a dedicated file mapping of their own and are
 *     returned to the OS with munmap() when freed;
 *   - all files are created with O_TMPFILE (or mkstemp + unlink) so nothing
 *     is ever left behind, even if the program is killed;
 *   - fresh file pages read as zero, so mm_calloc() of a fresh mapping
 *     costs nothing (no memset, no page is touched until it is used);
 *   - if the temporary directory is unusable the store degrades to plain
 *     calloc()/free() with a single warning.
 */
#ifndef MMSTORE_H
#define MMSTORE_H

#include <stddef.h>

/* configuration (call before the first allocation; all optional) */
void   mm_set_dir(const char *dir);      /* directory for the temp files   */
void   mm_set_enabled(int on);           /* 0 = plain heap (for testing)   */
int    mm_is_enabled(void);

/* allocation -- the returned memory is always zero-filled for mm_calloc(),
   and *not* guaranteed to be zero for mm_alloc()/mm_realloc() growth      */
void  *mm_alloc(size_t n);
void  *mm_calloc(size_t nmemb, size_t size);
void  *mm_realloc(void *p, size_t n);
void   mm_free(void *p);

/* contiguous 2-D array with a row-pointer table (both file backed).
   Row r starts at  (char*)m[0] + r*cols*elt .  Release with mm_matrix_free */
void **mm_matrix(size_t rows, size_t cols, size_t elt);
void   mm_matrix_free(void **m);

/* address-space reservation + fixed backing (used by the fragment arena) */
void  *mm_reserve(size_t len);
void  *mm_map_fixed(void *addr, size_t len);

/* diagnostics */
size_t mm_bytes_mapped(void);            /* bytes currently mapped         */
size_t mm_bytes_peak(void);              /* peak of the above              */
void   mm_report(void);                  /* one line on stderr             */

#endif
