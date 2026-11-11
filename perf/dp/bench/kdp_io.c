#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "kdp_kernel.h"
#include "aln_mem.h"
#include "aln_setup.h"

static int read_exact(FILE* f, void* buf, size_t n)
{
        return fread(buf, 1, n, f) == n ? 0 : 1;
}

int kdp_load(const char* path, struct kdp_problem** out)
{
        struct kdp_problem* p = NULL;
        FILE* f = NULL;
        int i;
        size_t na, nb;

        p = calloc(1, sizeof(*p));
        if(!p) return 1;
        snprintf(p->path, sizeof(p->path), "%s", path);
        f = fopen(path, "rb");
        if(!f){
                fprintf(stderr, "kdp: cannot open %s\n", path);
                goto fail;
        }
        if(read_exact(f, &p->h, sizeof(p->h)) || memcmp(p->h.magic, KDP_MAGIC, 8) != 0 ||
           p->h.version != KDP_VERSION){
                fprintf(stderr, "kdp: %s is not a v%d .kdp file\n", path, KDP_VERSION);
                goto fail;
        }
        na = (size_t)p->h.len_a;
        nb = (size_t)p->h.len_b;
        if(p->h.kind == KDP_SS){
                p->seq1 = malloc(na + 1);
                if(!p->seq1 || read_exact(f, p->seq1, na)) goto fail;
        }
        if(p->h.kind == KDP_SS || p->h.kind == KDP_SP){
                p->seq2 = malloc(nb + 1);
                if(!p->seq2 || read_exact(f, p->seq2, nb)) goto fail;
        }
        if(p->h.kind == KDP_SP || p->h.kind == KDP_PP){
                p->prof1 = malloc(sizeof(float) * (na + 2) * 64);
                if(!p->prof1 || read_exact(f, p->prof1, sizeof(float) * (na + 2) * 64)) goto fail;
        }
        if(p->h.kind == KDP_PP){
                p->prof2 = malloc(sizeof(float) * (nb + 2) * 64);
                if(!p->prof2 || read_exact(f, p->prof2, sizeof(float) * (nb + 2) * 64)) goto fail;
        }
        if(p->h.has_path){
                p->ref_path = malloc(sizeof(int32_t) * (na + 2));
                if(!p->ref_path || read_exact(f, p->ref_path, sizeof(int32_t) * (na + 2))) goto fail;
        }
        if(p->h.consistency){
                size_t nb = (size_t)p->h.bonus_rows * (size_t)p->h.bonus_K;
                p->bonus.n_rows = p->h.bonus_rows;
                p->bonus.K = p->h.bonus_K;
                p->bonus.cols = malloc(sizeof(int) * (nb ? nb : 1));
                p->bonus.vals = malloc(sizeof(float) * (nb ? nb : 1));
                if(!p->bonus.cols || !p->bonus.vals) goto fail;
                for(i = 0; i < (int)nb; i++){
                        int32_t v;
                        if(read_exact(f, &v, sizeof(v))) goto fail;
                        p->bonus.cols[i] = v;
                }
                if(read_exact(f, p->bonus.vals, sizeof(float) * nb)) goto fail;
                p->use_bonus = 1;
        }
        fclose(f);
        f = NULL;

        memset(&p->ap, 0, sizeof(p->ap));
        p->ap.nthreads = 1;
        for(i = 0; i < 23; i++){
                p->subm_rows[i] = &p->h.subm[i * 23];
        }
        p->ap.subm = p->subm_rows;
        p->ap.gpo = p->h.gpo;
        p->ap.gpe = p->h.gpe;
        p->ap.tgpe = p->h.tgpe;
        p->ap.subm_offset = p->h.subm_offset;
        *out = p;
        return 0;
fail:
        if(f) fclose(f);
        kdp_free(p);
        return 1;
}

void kdp_free(struct kdp_problem* p)
{
        if(!p) return;
        free(p->seq1);
        free(p->seq2);
        free(p->prof1);
        free(p->prof2);
        free(p->ref_path);
        free(p->bonus.cols);
        free(p->bonus.vals);
        free(p);
}

void kdp_mem_reset(const struct kdp_problem* p, struct aln_mem* m)
{
        int i, g;
        m->starta = 0;
        m->startb = 0;
        m->enda = p->h.len_a;
        m->endb = p->h.len_b;
        m->starta_2 = 0;
        m->enda_2 = p->h.len_a;
        m->len_a = p->h.len_a;
        m->len_b = p->h.len_b;
        m->f[0].a = p->h.f0[0];
        m->f[0].ga = p->h.f0[1];
        m->f[0].gb = p->h.f0[2];
        m->b[0].a = p->h.b0[0];
        m->b[0].ga = p->h.b0[1];
        m->b[0].gb = p->h.b0[2];
        m->margin_sum = 0.0F;
        m->margin_count = 0;
        m->flip_threshold = 0.0F;
        m->flip_counter = 0;
        m->consistency = (p->h.consistency && p->use_bonus) ? (struct sparse_bonus*)&p->bonus : NULL;
        m->mode = ALN_MODE_FULL;
        g = (p->h.len_a > p->h.len_b ? p->h.len_a : p->h.len_b) + 2;
        for(i = 0; i < g; i++){
                m->path[i] = -1;
        }
}

int kdp_mem_create(const struct kdp_problem* p, struct aln_mem** out)
{
        struct aln_mem* m = NULL;
        if(alloc_aln_mem(&m, 256) != 0) return 1;
        m->ap = (struct aln_param*)&p->ap;
        m->seq1 = p->seq1;
        m->seq2 = p->seq2;
        m->prof1 = p->prof1;
        m->prof2 = p->prof2;
        m->sip = p->h.sip;
        m->run_parallel = 0;
        m->len_a = p->h.len_a;
        m->len_b = p->h.len_b;
        if(init_alnmem(m) != 0){
                free_aln_mem(m);
                return 1;
        }
        kdp_mem_reset(p, m);
        *out = m;
        return 0;
}

void kdp_mem_free(struct aln_mem* m)
{
        free_aln_mem(m);
}
