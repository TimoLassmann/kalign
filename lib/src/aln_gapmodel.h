/* Per-member gap model (run option gap_model = 1), see aln_gapmodel.c. */
#ifndef ALN_GAPMODEL_H
#define ALN_GAPMODEL_H

struct aln_mem;
struct aln_param;

/* profile column slots holding the member counts (unused by the default model) */
#define GM_NRES  58   /* members with a residue in the column */
#define GM_LEAD  59   /* members with no residue at or before the column (before their first residue) */
#define GM_TRAIL 60   /* members with no residue at or after the column (after their last residue) */
#define GM_GOPEN 61   /* members with an internal gap here and a residue in the previous column */
#define GM_GEXT  62   /* members with an internal gap here and a gap in the previous column */

int gm_prepare(struct aln_mem* m);
void gm_free(struct aln_mem* m);

int gm_pp_forward(struct aln_mem* m);
int gm_pp_backward(struct aln_mem* m);
int gm_pp_meetup(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score);
int gm_sp_forward(struct aln_mem* m);
int gm_sp_backward(struct aln_mem* m);
int gm_sp_meetup(struct aln_mem* m, int old_cor[], int* meet, int* t, float* score);

void gm_leaf_counts(float* prof, int len);
int update_n_gm(const float* profa, const float* profb, float* newp, struct aln_param* ap, int* path, int sipa, int sipb);

#endif
