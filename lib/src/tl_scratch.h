#ifndef TL_SCRATCH_H
#define TL_SCRATCH_H

/* Scratch buffers kept in _Thread_local variables are not freed when a thread
   ends.  After (re)allocating them, a thread records them with
   tl_scratch_track(slot, ptrs, n) (n <= TL_SCRATCH_PTRS); they are freed when the
   thread exits.  One slot per kind of buffer set. */
enum { TL_SCRATCH_WAVE_SS, TL_SCRATCH_WAVE_SP, TL_SCRATCH_WAVE_PP, TL_SCRATCH_GAPMODEL, TL_SCRATCH_SLOTS };
#define TL_SCRATCH_PTRS 5
void tl_scratch_track(int slot, void* const* ptrs, int n);

#endif
