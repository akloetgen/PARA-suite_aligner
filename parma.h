/*
 * parma.h
 *
 *  Created on: 20.08.2014
 *      Author: akloetgen
 */

#ifndef PARMA_H_
#define PARMA_H_

#include "bwt.h"
#include "bwtaln.h"

typedef struct { // recursion stack
	u_int32_t info; // score<<21 | i
	u_int32_t n_mm:8, n_gapo:8, n_gape:8, state:2, n_seed_mm:6;
	u_int32_t n_ins:16, n_del:16;
	int last_diff_pos;
	double ep_p_val;
	bwtint_t k, l; // (k,l) is the SA region of [i,n-1]
} gap_entry_t_parma;

typedef struct {
	int n_entries, m_entries;
	gap_entry_t_parma *stack;
} gap_stack1_t_parma;

typedef struct {
	int n_stacks, n_entries;
	int best;
	gap_stack1_t_parma *stacks;
} gap_stack_t_parma;

#ifdef __cplusplus
extern "C" {
#endif

	gap_stack_t_parma *gap_init_stack2_parma(int max_score);
	gap_stack_t_parma *gap_init_stack_parma(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt);
	void gap_destroy_stack_parma(gap_stack_t_parma *stack);
	bwt_aln1_t *bwt_match_gap_parma(bwt_t *const bwt, int len, const ubyte_t *seq, bwt_width_t *w,
							  bwt_width_t *seed_w, const gap_opt_t *opt, int *_n_aln, gap_stack_t_parma *stack);

	// DELETE????
	void bwa_aln2seq_ep(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);

#ifdef __cplusplus
}
#endif


#endif /* BWTGAPEP_H_ */
