/*
 *  mmstore.c  --  file-backed ("mmapped") memory store, see mmstore.h
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mmstore.h"

#ifndef O_TMPFILE
#define O_TMPFILE 0
#endif

#define MM_HDR        16                 /* header in front of every block */
#define MM_CHUNK      ((size_t)64 << 20) /* file-backed arena chunk        */
#define MM_LARGE      ((size_t)1  << 20) /* >= this: dedicated mapping     */
#define MM_NCLASS     96

enum { TAG_SMALL = 1, TAG_LARGE = 2, TAG_HEAP = 3 };

typedef struct {
    size_t size;      /* usable bytes (small: class size, large: map length-hdr) */
    size_t tag;
} mm_hdr;

typedef struct { void *next; } mm_free_node;

typedef struct {
    void  *freelist[MM_NCLASS];
    char  *bump, *end;
} mm_thread_heap;

static __thread mm_thread_heap th;

static pthread_mutex_t mm_lock = PTHREAD_MUTEX_INITIALIZER;
static char   mm_dir[4096]     = "";
static int    mm_on            = 1;
static int    mm_warned        = 0;
static int    mm_sig_installed = 0;
static int    mm_debug         = -1;   /* DIALIGN_MMDEBUG=1: poison freed blocks */
static size_t mm_mapped        = 0;
static size_t mm_peak          = 0;
static size_t mm_nchunks       = 0;

/* ------------------------------------------------------------------ */

static void mm_sigbus(int sig)
{
    static const char msg[] =
      "\ndialign: fatal I/O error on the temporary storage (disk full?).\n"
      "         Use  -tmpdir <dir>  to point it at a volume with more room.\n";
    (void)sig;
    if (write(2, msg, sizeof msg - 1) < 0) { /* nothing more to do */ }
    _exit(1);
}

static void mm_install_sigbus(void)
{
    struct sigaction sa;
    if (mm_sig_installed) return;
    mm_sig_installed = 1;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = mm_sigbus;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGBUS, &sa, NULL);
}

void mm_set_dir(const char *dir)
{
    if (dir && *dir) {
        strncpy(mm_dir, dir, sizeof mm_dir - 1);
        mm_dir[sizeof mm_dir - 1] = '\0';
    }
}

void mm_set_enabled(int on) { mm_on = on; }
int  mm_is_enabled(void)    { return mm_on; }

static const char *mm_get_dir(void)
{
    const char *d;
    if (mm_dir[0]) return mm_dir;
    if ((d = getenv("DIALIGN_TMPDIR")) != NULL && *d) return d;
    return ".";
}

/* create an anonymous temp file of `len` bytes and map it (at `fixed` if
   non-NULL, replacing whatever is mapped there); NULL on failure */
static void *mm_map_file_at(void *fixed, size_t len)
{
    const char *dir = mm_get_dir();
    int fd = -1;
    void *p;

#if defined(__linux__)
    if (O_TMPFILE) fd = open(dir, O_TMPFILE | O_RDWR, 0600);
#endif
    if (fd < 0) {
        char path[4200];
        snprintf(path, sizeof path, "%s/.dialign-mm-XXXXXX", dir);
        fd = mkstemp(path);
        if (fd >= 0) unlink(path);
    }
    if (fd < 0) return NULL;
    if (ftruncate(fd, (off_t)len) != 0) { close(fd); return NULL; }
    p = mmap(fixed, len, PROT_READ | PROT_WRITE,
             MAP_SHARED | ( fixed ? MAP_FIXED : 0 ), fd, 0);
    close(fd);                                   /* mapping stays valid */
    if (p == MAP_FAILED) return NULL;
    mm_install_sigbus();
    return p;
}

static void *mm_map_file(size_t len) { return mm_map_file_at(NULL, len); }

static void mm_account(long delta)
{
    if (delta >= 0) {
        size_t v = __atomic_add_fetch(&mm_mapped, (size_t)delta, __ATOMIC_RELAXED);
        size_t pk = __atomic_load_n(&mm_peak, __ATOMIC_RELAXED);
        while (v > pk &&
               !__atomic_compare_exchange_n(&mm_peak, &pk, v, 0,
                                            __ATOMIC_RELAXED, __ATOMIC_RELAXED))
            ;
    } else {
        __atomic_sub_fetch(&mm_mapped, (size_t)(-delta), __ATOMIC_RELAXED);
    }
}

static void mm_warn_fallback(void)
{
    if (!mm_warned) {
        mm_warned = 1;
        fprintf(stderr,
          "dialign: warning: cannot create temporary files in '%s'; "
          "falling back to ordinary heap memory.\n", mm_get_dir());
    }
}

/* ------------------------------------------------------------------ */
/*  size classes                                                       */
/* ------------------------------------------------------------------ */

static size_t class_round(size_t n)
{
    if (n < 16) n = 16;
    if (n <= 256) return (n + 15) & ~(size_t)15;
    {
        int msb = 63 - __builtin_clzl(n - 1);
        size_t step = (size_t)1 << (msb - 2);
        return (n + step - 1) & ~(step - 1);
    }
}

static int class_index(size_t r)          /* r already rounded */
{
    if (r <= 256) return (int)(r / 16) - 1;                 /* 0..15  */
    {
        int msb = 63 - __builtin_clzl(r - 1);
        size_t sub = (r - ((size_t)1 << msb)) / ((size_t)1 << (msb - 2));
        return 16 + (msb - 8) * 4 + (int)sub;               /* 17..   */
    }
}

/* ------------------------------------------------------------------ */

static void *heap_alloc(size_t n, int zero)
{
    mm_hdr *h = zero ? calloc(1, n + MM_HDR) : malloc(n + MM_HDR);
    if (!h) {
        fprintf(stderr, "dialign: out of memory (%zu bytes)\n", n);
        exit(1);
    }
    h->size = n;
    h->tag  = TAG_HEAP;
    return (char *)h + MM_HDR;
}

static void *large_alloc(size_t n)
{
    size_t len, pg = (size_t)sysconf(_SC_PAGESIZE);
    mm_hdr *h;
    len = (n + MM_HDR + pg - 1) & ~(pg - 1);
    h = (mm_hdr *)mm_map_file(len);
    if (!h) return NULL;
    h->size = len - MM_HDR;
    h->tag  = TAG_LARGE;
    mm_account((long)len);
    return (char *)h + MM_HDR;
}

static void *small_alloc(size_t n, int *fresh)
{
    size_t r = class_round(n);
    int ci = class_index(r);
    mm_hdr *h;

    if (ci >= MM_NCLASS) return NULL;
    if (th.freelist[ci]) {
        mm_free_node *f = (mm_free_node *)th.freelist[ci];
        th.freelist[ci] = f->next;
        *fresh = 0;
        return f;
    }
    if (th.bump == NULL || (size_t)(th.end - th.bump) < r + MM_HDR) {
        char *c = (char *)mm_map_file(MM_CHUNK);
        if (!c) return NULL;
        pthread_mutex_lock(&mm_lock);
        mm_nchunks++;
        pthread_mutex_unlock(&mm_lock);
        mm_account((long)MM_CHUNK);
        th.bump = c;
        th.end  = c + MM_CHUNK;
    }
    h = (mm_hdr *)th.bump;
    th.bump += r + MM_HDR;
    h->size = r;
    h->tag  = TAG_SMALL;
    *fresh = 1;
    return (char *)h + MM_HDR;
}

static void *mm_alloc_impl(size_t n, int zero)
{
    void *p;
    int fresh = 1;

    if (n == 0) n = 1;
    if (n > ((size_t)-1) / 2) {
        fprintf(stderr, "dialign: allocation request too large\n");
        exit(1);
    }
    if (!mm_on) return heap_alloc(n, zero);

    if (n >= MM_LARGE) {
        p = large_alloc(n);                         /* fresh => zero */
    } else {
        p = small_alloc(n, &fresh);
        if (p && zero && !fresh) memset(p, 0, ((mm_hdr *)((char *)p - MM_HDR))->size);
    }
    if (!p) {                                       /* degrade gracefully */
        mm_warn_fallback();
        mm_on = 0;
        return heap_alloc(n, zero);
    }
    return p;
}

void *mm_alloc(size_t n)                { return mm_alloc_impl(n, 0); }

void *mm_calloc(size_t nmemb, size_t size)
{
    if (size != 0 && nmemb > ((size_t)-1) / 2 / size) {
        fprintf(stderr, "dialign: allocation request too large\n");
        exit(1);
    }
    return mm_alloc_impl(nmemb * size, 1);
}

void mm_free(void *p)
{
    mm_hdr *h;
    if (!p) return;
    h = (mm_hdr *)((char *)p - MM_HDR);
    switch (h->tag) {
    case TAG_SMALL: {
        int ci = class_index(h->size);
        if( mm_debug < 0 ) mm_debug = ( getenv("DIALIGN_MMDEBUG") != NULL );
        if( mm_debug ) memset( p , 0xDD , h->size );
        mm_free_node *f = (mm_free_node *)p;
        f->next = th.freelist[ci];
        th.freelist[ci] = f;
        break; }
    case TAG_LARGE: {
        size_t len = h->size + MM_HDR;
        munmap(h, len);
        mm_account(-(long)len);
        break; }
    case TAG_HEAP:
        free(h);
        break;
    default:
        fprintf(stderr, "dialign: mm_free() of a pointer not owned by the store\n");
        abort();
    }
}

void *mm_realloc(void *p, size_t n)
{
    mm_hdr *h;
    void *q;
    size_t keep;
    if (!p) return mm_alloc(n);
    h = (mm_hdr *)((char *)p - MM_HDR);
    keep = h->size;
    if (n <= keep && h->tag != TAG_LARGE) return p;   /* small: fits already */
    q = mm_alloc(n);
    memcpy(q, p, keep < n ? keep : n);
    if( n > keep ) memset( (char *)q + keep , 0 , n - keep );   /* grown part reads as zero */
    mm_free(p);
    return q;
}

void **mm_matrix(size_t rows, size_t cols, size_t elt)
{
    void **m;
    char *data;
    size_t r;
    m = (void **)mm_alloc((rows ? rows : 1) * sizeof(void *));
    data = (char *)mm_calloc((rows ? rows : 1) * (cols ? cols : 1), elt ? elt : 1);
    for (r = 0; r < rows; r++) m[r] = data + r * cols * elt;
    if (rows == 0) m[0] = data;
    return m;
}

void mm_matrix_free(void **m)
{
    if (!m) return;
    mm_free(m[0]);          /* contiguous data block */
    mm_free(m);
}

size_t mm_bytes_mapped(void) { return __atomic_load_n(&mm_mapped, __ATOMIC_RELAXED); }
size_t mm_bytes_peak(void)   { return __atomic_load_n(&mm_peak,   __ATOMIC_RELAXED); }

void mm_report(void)
{
    fprintf(stderr, "dialign: mmap store: peak %.1f MiB mapped in %s (%zu chunks)\n",
            (double)mm_bytes_peak() / (1 << 20), mm_get_dir(), mm_nchunks);
}


/* ---- helpers for the fragment arena (frags.c) ---- */
void *mm_reserve(size_t len)          /* address space only, no memory */
{
    void *p = mmap(NULL, len, PROT_NONE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
    return p == MAP_FAILED ? NULL : p;
}

void *mm_map_fixed(void *addr, size_t len)   /* back [addr,addr+len) with storage */
{
    void *p;
    if (mm_on) {
        p = mm_map_file_at(addr, len);
        if (p) { mm_account((long)len); return p; }
        mm_warn_fallback();
        mm_on = 0;
    }
    p = mmap(addr, len, PROT_READ | PROT_WRITE,
             MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE | MAP_FIXED, -1, 0);
    return p == MAP_FAILED ? NULL : p;
}
