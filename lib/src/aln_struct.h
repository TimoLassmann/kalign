#ifndef ALN_STRUCT_H
#define ALN_STRUCT_H

#include <stdint.h>

#define ALN_MODE_SCORE_ONLY 2
#define ALN_MODE_FULL 1

struct states{
        float a;
        float ga;
        float gb;
        // float x;
};

struct aln_mem{
        const float* prof1;
        const float* prof2;
        const uint8_t* seq1;
        const uint8_t* seq2;
        struct aln_param* ap;
        struct states* f;
        struct states* b;
        int* path;
        int* tmp_path;
        uint8_t run_parallel;
        int alloc_path_len;
        float score;
        float margin_sum;       /* accumulated meetup margins */
        int margin_count;       /* number of meetup calls */

        float flip_threshold;   /* midpoints with margin < this are flip candidates; 0 = no flips */
        int flip_trial;         /* round-robin: current trial (1..K-1, 0 = baseline) */
        int flip_stride;        /* round-robin: number of flip slots */
        int flip_counter;       /* running count of flip candidates encountered */
        uint32_t flip_mask;     /* bitmask of which slots to flip (0 = no flips) */

        float* flip_margins;    /* per-meetup margins recorded during baseline */
        int flip_margin_alloc;  /* allocated size of flip_margins */
        int* flip_bit_map;      /* maps flip_counter → bit index (-1 = not targeted) */
        int flip_n_targets;     /* number of individually targeted midpoints */
        int flip_n_uncertain;   /* total uncertain midpoints (margin < threshold) */

        int starta;
        int starta_2;
        int startb;
        int enda;
        int enda_2;
        int endb;
        int size;
        int len_a;
        int len_b;

        int sip;
        int mode;
};

#endif
