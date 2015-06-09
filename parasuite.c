/*
 * parasuite.c
 *
 *  Created on: 20.08.2014
 *      Author: akloetgen
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parasuite.h"
#include "parasuitecore.h"
#include "bwtaln.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

#define aln_score_parasuite(m,o,e,p) ((m)*(p)->s_mm + (o)*(p)->s_gapo + (e)*(p)->s_gape)
//#define aln_score_ep(X,i,p) (pow((p)->avg_match, i-X) * pow((p)->avg_mm, X))

gap_stack_t_parasuite *gap_init_stack2_parasuite(int max_score) {
	gap_stack_t_parasuite *stack;
	stack = (gap_stack_t_parasuite*) calloc(1, sizeof(gap_stack_t_parasuite));
	//fprintf(stderr, "max_score=%d\n", max_score);
	stack->n_stacks = max_score;
	stack->stacks = (gap_stack1_t_parasuite*) calloc(stack->n_stacks,
			sizeof(gap_stack1_t_parasuite));
	return stack;
}

gap_stack_t_parasuite *gap_init_stack_parasuite(int max_mm, int max_gapo, int max_gape,
		const gap_opt_t *opt) {
	return gap_init_stack2_parasuite(
			aln_score_parasuite(max_mm + 1, max_gapo + 1, max_gape + 1, opt));
}

void gap_destroy_stack_parasuite(gap_stack_t_parasuite *stack) {
	int i;
	for (i = 0; i != stack->n_stacks; ++i)
		free(stack->stacks[i].stack);
	free(stack->stacks);
	free(stack);
}

static void gap_reset_stack_parasuite(gap_stack_t_parasuite *stack) {
	int i;
	for (i = 0; i != stack->n_stacks; ++i)
		stack->stacks[i].n_entries = 0;
	stack->best = stack->n_stacks;
	stack->n_entries = 0;
}

static inline void gap_push_parasuite(gap_stack_t_parasuite *stack, int i, bwtint_t k,
		bwtint_t l, int n_mm, int n_gapo, int n_gape, int n_ins, int n_del,
		int state, int is_diff, const gap_opt_t *opt, double ep_p_val, int len) {
	//int score;
	gap_entry_t_parasuite *p;
	gap_stack1_t_parasuite *q;
	int score = aln_score_parasuite(n_mm, n_gapo, n_gape, opt);
	//double score = aln_score_ep(opt->X, len - i, opt);
	//double score = pow(opt->avg_match, len - i) * pow(opt->avg_mm, n_mm);
	q = stack->stacks + score;
	if (q->n_entries == q->m_entries) {
		q->m_entries = q->m_entries ? q->m_entries << 1 : 4;
		q->stack = (gap_entry_t_parasuite*) realloc(q->stack,
				sizeof(gap_entry_t_parasuite) * q->m_entries);
	}
	p = q->stack + q->n_entries;
	p->info = (u_int32_t) score << 21 | i;
	p->k = k;
	p->l = l;
	p->n_mm = n_mm;
	p->n_gapo = n_gapo;
	p->n_gape = n_gape;
	p->n_ins = n_ins;
	p->n_del = n_del;
	p->state = state;
	p->last_diff_pos = is_diff ? i : 0;
	p->ep_p_val = ep_p_val;
	++(q->n_entries);
	++(stack->n_entries);
	//fprintf("current score=%d and best=%d\n", score, stack->best);
	if (stack->best > score)
		stack->best = score;
}

static inline void gap_pop_parasuite(gap_stack_t_parasuite *stack, gap_entry_t_parasuite *e) {
	gap_stack1_t_parasuite *q;
	q = stack->stacks + stack->best;
	*e = q->stack[q->n_entries - 1];
	--(q->n_entries);
	--(stack->n_entries);
	if (q->n_entries == 0 && stack->n_entries) { // reset best
		int i;
		for (i = stack->best + 1; i < stack->n_stacks; ++i)
			if (stack->stacks[i].n_entries != 0)
				break;
		stack->best = i;
	} else if (stack->n_entries == 0)
		stack->best = stack->n_stacks;
}

static inline void gap_shadow_parasuite(int x, int len, bwtint_t max,
		int last_diff_pos, bwt_width_t *w) {
	int i, j;
	for (i = j = 0; i < last_diff_pos; ++i) {
		if (w[i].w > x)
			w[i].w -= x;
		else if (w[i].w == x) {
			w[i].bid = 1;
			w[i].w = max - (++j);
		} // else should not happen
	}
}

static inline int int_log2_parasuite(uint32_t v) {
	int c = 0;
	if (v & 0xffff0000u) {
		v >>= 16;
		c |= 16;
	}
	if (v & 0xff00) {
		v >>= 8;
		c |= 8;
	}
	if (v & 0xf0) {
		v >>= 4;
		c |= 4;
	}
	if (v & 0xc) {
		v >>= 2;
		c |= 2;
	}
	if (v & 0x2)
		c |= 1;
	return c;
}

bwt_aln1_t *bwt_match_gap_parasuite(bwt_t * const bwt, int len, const ubyte_t *seq,
		bwt_width_t *width, bwt_width_t *seed_width, const gap_opt_t *opt,
		int *_n_aln, gap_stack_t_parasuite *stack) { // $seq is the reverse complement of the input read

	// p_value for current read mapping step
	double p_threshold = (double) pow((double) opt->avg_match, (len - opt->X))
			* pow((double) opt->avg_mm, opt->X);
	// calc number of maximal mms allowed:

	/*double temp_best_mm = 0.0;
	for (i = 0; i < 4; i++) {
		for (k = 0; k < 4; k++) {
			if (i != k && opt->profile.position_profile[i][k] > temp_best_mm) {
				temp_best_mm = opt->profile.position_profile[i][k];
			}
		}
	}
	fprintf(stderr, "temp_best_mm=%f; best_mm=%f\n", temp_best_mm, opt->best_mm);*/
	double temp_val = 1.0;
	int max_mm = 0;
	while (temp_val > p_threshold) {
		temp_val = (double) pow((double) opt->avg_match, (len - max_mm))
				* pow((double) opt->best_mm, max_mm);
		max_mm++;
	}

	int best_score = aln_score_parasuite(max_mm, opt->max_gapo + 1, opt->max_gape + 1,
			opt);
	int best_diff = max_mm + 1;	//, max_diff = opt->max_diff;
	int max_diff = max_mm;
	int best_cnt = 0;
	int max_entries = 0, j, _j, n_aln, m_aln;
	bwt_aln1_t *aln;

	m_aln = 4;
	n_aln = 0;
	aln = (bwt_aln1_t*) calloc(m_aln, sizeof(bwt_aln1_t));

	// check whether there are too many N
	for (j = _j = 0; j < len; ++j)
		if (seq[j] > 3)
			++_j;
	if (_j > max_diff) {
		*_n_aln = n_aln;
		return aln;
	}

	//for (j = 0; j != len; ++j) printf("#0 %d: [%d,%u]\t[%d,%u]\n", j, w[0][j].bid, w[0][j].w, w[1][j].bid, w[1][j].w);
	gap_reset_stack_parasuite(stack); // reset stack
	/*int bla[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	 int bla2[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	 gap_push_parasuite(stack, len, 0, bwt->seq_len, 0, 0, 0, 0, 0, 0, 0, opt, bla, bla2, 0, -1);*/
	gap_push_parasuite(stack, len, 0, bwt->seq_len, 0, 0, 0, 0, 0, 0, 0, opt, 1.0,
			len);

	//fprintf(stderr, "bla error\n");

	while (stack->n_entries) {

		gap_entry_t_parasuite e;
		int i, m, m_seed = 0, hit_found, allow_diff, allow_M, tmp;
		bwtint_t k, l, cnt_k[4], cnt_l[4], occ;

		//fprintf(stderr, "stack->n_entries=%d\n", stack->n_entries);

		if (max_entries < stack->n_entries)
			max_entries = stack->n_entries;
		if (stack->n_entries > opt->max_entries)
			break;
		gap_pop_parasuite(stack, &e); // get the best entry
		k = e.k;
		l = e.l; // SA interval
		i = e.info & 0xffff; // length
		//fprintf(stderr, "i=%d with char %d\n",i,seq[i]);
		//fprintf(stderr, "info of gap entry: %d\n",i);
		//fprintf(stderr, "k=%ld; l=%ld; i=%d\n",k,l,i);

		// HIER MUSS NEUER VERGLEICH VS. ERRORPROFIL HER UND NICHT VS BEST_SCORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		/*if (!(opt->mode & BWA_MODE_NONSTOP)
				&& e.info >> 21 > best_score + opt->s_mm) {
			fprintf(stderr,"ABBURCH: no need to proceed - VERBESSERN MIT EP\n");
			break; // no need to proceed
		}*/
		if (!(opt->mode & BWA_MODE_NONSTOP)
				&& e.info >> 21 > best_score + opt->s_mm) {
			//fprintf(stderr,"ABBURCH: no need to proceed - VERBESSERN MIT EP\n");
			break; // no need to proceed
		}

		m = max_diff - (e.n_mm + e.n_gapo);
		if (opt->mode & BWA_MODE_GAPE)
			m -= e.n_gape;
		//if (m < 0) continue;

		// CHECK EXCLUDED; MAYBE CHECK IN AGAIN???
		/*if (e.ep_p_val < p_threshold) {
			continue;
		}*/
		if (seed_width) { // apply seeding
			m_seed = opt->max_seed_diff - (e.n_mm + e.n_gapo);
			if (opt->mode & BWA_MODE_GAPE)
				m_seed -= e.n_gape;
		}
		//printf("#1\t[%d,%d,%d,%c]\t[%d,%d,%d]\t[%u,%u]\t[%u,%u]\t%d\n", stack->n_entries, a, i, "MID"[e.state], e.n_mm, e.n_gapo, e.n_gape, width[i-1].bid, width[i-1].w, k, l, e.last_diff_pos);
		//if (i > 0 && m < width[i-1].bid) continue;

		// check whether a hit is found
		hit_found = 0;

		if (i == 0) {
			hit_found = 1;
			//fprintf(stderr, "hit found for i==0 with %d errors\n", e.n_mm);
			//fprintf(stderr, "ended up with p-val=%1.10f with last diff at %d\n", e.ep_p_val, e.last_diff_pos);
			// ONLY IF NOT ALL MMs WERE NEEDED TO FIND MATCH
		} else if (m == 0
				&& (e.state == STATE_M || (opt->mode & BWA_MODE_GAPE)
						|| e.n_gape == opt->max_gape)) { // no diff allowed

			if (bwt_match_exact_alt(bwt, i, seq, &k, &l)) {
				// FINDS HIT IF MAXIMAL NUMBER OF MM IS ALREADY APPLIED
				hit_found = 1;

				//fprintf(stderr, "hit found starting from i=%d and len=%d without error with %d errors\n",i, len, e.n_mm);
				//fprintf(stderr, "ended up with p-val=%1.10f with last diff at %d\n", e.ep_p_val, e.last_diff_pos);
			} else
				continue; // no hit, skip
		}

		if (hit_found) { // action for found hits
			//fprintf(stderr, "number entries passed before hit was found: stack->n_entries=%d\n", stack->n_entries);
			int score = aln_score_parasuite(e.n_mm, e.n_gapo, e.n_gape, opt);
			int do_add = 1;
			//printf("#2 hits found: %d:(%u,%u)\n", e.n_mm+e.n_gapo, k, l);
			if (n_aln == 0) {
				best_score = score;
				best_diff = e.n_mm + e.n_gapo;
				if (opt->mode & BWA_MODE_GAPE)
					best_diff += e.n_gape;
				//if (!(opt->mode & BWA_MODE_NONSTOP))
				//	max_diff = (best_diff + 1 > opt->max_diff)? opt->max_diff : best_diff + 1; // top2 behaviour
			}
			if (score == best_score)
				best_cnt += l - k + 1;
			else if (best_cnt > opt->max_top2) {
				//fprintf(stderr, "BREAK 2nd!!!!!!\n");
				break; // top2b behaviour
			}
			if (e.n_gapo) { // check whether the hit has been found. this may happen when a gap occurs in a tandem repeat
				for (j = 0; j != n_aln; ++j)
					if (aln[j].k == k && aln[j].l == l)
						break;
				if (j < n_aln)
					do_add = 0;
			}
			if (do_add) { // append
				bwt_aln1_t *p;
				gap_shadow_parasuite(l - k + 1, len, bwt->seq_len, e.last_diff_pos,
						width);
				if (n_aln == m_aln) {
					m_aln <<= 1;
					aln = (bwt_aln1_t*) realloc(aln,
							m_aln * sizeof(bwt_aln1_t));
					memset(aln + m_aln / 2, 0, m_aln / 2 * sizeof(bwt_aln1_t));
				}
				p = aln + n_aln;
				p->n_mm = e.n_mm;
				p->n_gapo = e.n_gapo;
				p->n_gape = e.n_gape;
				p->n_ins = e.n_ins;
				p->n_del = e.n_del;
				p->k = k;
				p->l = l;
				p->score = score;
				//fprintf(stderr, "*** n_mm=%d,n_gapo=%d,n_gape=%d,n_ins=%d,n_del=%d\n", e.n_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del);
				++n_aln;
			}
			continue;
		}

		--i;
		bwt_2occ4(bwt, k - 1, l, cnt_k, cnt_l); // retrieve Occ values
		occ = l - k + 1;
		// test whether diff is allowed
		allow_diff = allow_M = 1;
		// VERBIETTET DAS NICHT IN DER REGEL DEN 3. T-C ?!?!?!??!?!?!?!?!?!?!?!?!?
		//if (e.ep_p_val < (p_threshold * opt->avg_mm)) {

		/*
		// OLD CHECK FOR HOW MANY MMs ARE STILL ALLOWED
		if (i > 0) {
			int ii = i - (len - opt->seed_len);
			if (width[i-1].bid > m-1) allow_diff = 0;
			else if (width[i-1].bid == m-1 && width[i].bid == m-1 && width[i-1].w == width[i].w) allow_M = 0;
			if (seed_width && ii > 0) {
				if (seed_width[ii-1].bid > m_seed-1) allow_diff = 0;
				else if (seed_width[ii-1].bid == m_seed-1 && seed_width[ii].bid == m_seed-1
						&& seed_width[ii-1].w == seed_width[ii].w) allow_M = 0;
			}
		}*/
		if (p_threshold >= (e.ep_p_val * opt->best_mm) || p_threshold == 0) {
			//fprintf(stderr, "p-val too small at i=%d\n", i);
			allow_diff = allow_M = 0;
		}

		// indels
		tmp = (opt->mode & BWA_MODE_LOGGAP) ?
				int_log2_parasuite(e.n_gape + e.n_gapo) / 2 + 1 : e.n_gapo + e.n_gape;
		if (allow_diff && i >= opt->indel_end_skip + tmp
				&& len - i >= opt->indel_end_skip + tmp) {
			if (e.state == STATE_M) { // gap open
				//fprintf(stderr, "opt->max_gapo=%d\n",opt->max_gapo);
				//if (e.n_gapo < opt->max_gapo) { // gap open is allowed
				//if (e.n_gapo < opt->X) { // gap open is allowed
					// insertion; check whether threshold would be exceeded
					if ((double) (e.ep_p_val * opt->profile.indels[INSERTION])
							> p_threshold) {
						//	fprintf(stderr, "insert insertion");
						gap_push_parasuite(stack, i, k, l, e.n_mm, e.n_gapo + 1,
								e.n_gape, e.n_ins + 1, e.n_del, STATE_I, 1, opt,
								e.ep_p_val * opt->profile.indels[INSERTION],
								len);
					}
					// deletion; check whether threshold would be exceeded
					/*fprintf(stderr,
							"e.ep_p_val=%.10f, del-value=%.10f p-val with insert/Deletion=%.10f abd p_thres=%.10f\n",
							e.ep_p_val, opt->profile.indels[DELETION],
							(double) (e.ep_p_val * opt->profile.indels[DELETION]),p_threshold);*/
					if ((double) (e.ep_p_val * opt->profile.indels[DELETION])
							> p_threshold) {
						//fprintf(stderr, "insert insertion");
						for (j = 0; j != 4; ++j) {
							k = bwt->L2[j] + cnt_k[j] + 1;
							l = bwt->L2[j] + cnt_l[j];
							if (k <= l)
								gap_push_parasuite(stack, i + 1, k, l, e.n_mm,
										e.n_gapo + 1, e.n_gape, e.n_ins,
										e.n_del + 1, STATE_D, 1, opt,
										e.ep_p_val
												* opt->profile.indels[DELETION],
										len);
						}
					}
				//}
			} else if (e.state == STATE_I) { // extention of an insertion
				if (e.n_gape < opt->max_gape) // gap extention is allowed
					gap_push_parasuite(stack, i, k, l, e.n_mm, e.n_gapo, e.n_gape + 1,
							e.n_ins + 1, e.n_del, STATE_I, 1, opt, e.ep_p_val,
							len);
			} else if (e.state == STATE_D) { // extention of a deletion
				if (e.n_gape < opt->max_gape) { // gap extention is allowed
					if (e.n_gape + e.n_gapo < max_diff
							|| occ < opt->max_del_occ) {
						for (j = 0; j != 4; ++j) {
							k = bwt->L2[j] + cnt_k[j] + 1;
							l = bwt->L2[j] + cnt_l[j];
							if (k <= l)
								gap_push_parasuite(stack, i + 1, k, l, e.n_mm,
										e.n_gapo, e.n_gape + 1, e.n_ins,
										e.n_del + 1, STATE_D, 1, opt,
										e.ep_p_val, len);
						}
					}
				}
			}
		}
		// mismatches
		if (allow_diff && allow_M) { // mismatch is allowed
			//fprintf(stderr, "check for mut at i=%d with char %d\n", i, seq[i]);
			double new_ep_p_val;
			for (j = 1; j <= 4; ++j) {
				int c = (seq[i] + j) & 3;
				// CHECK WHETHER c (mm char) IS APPLICAPLE REGARDING THE p-VALUE!!!! AND UPDATE p!!!!

				/*else {
				 continue;
				 }*/

				int is_mm = (j != 4 || seq[i] > 3);
				k = bwt->L2[c] + cnt_k[c] + 1;
				l = bwt->L2[c] + cnt_l[c];
				if (k <= l) {
					/*if (c == 1 && seq[i] == 3) {
					 fprintf(stderr, "T-C mut at i=%d\n",i);
					 }*/
					//int refmm, readmm, last_err;
					if ((new_ep_p_val = e.ep_p_val * opt->profile.position_profile[3 - c][3 - seq[i]]) < p_threshold) {
						continue;
					}
					if (is_mm) {
						//e.ep_p_val *= opt->profile.position_profile[3-c][3-seq[i]];


						//fprintf(stderr, "at i=%d calc p-val: current=%1.10f with %1.10f results in new=%1.10f\n", i, e.ep_p_val,opt->profile.position_profile[3-c][3-seq[i]], e.ep_p_val * opt->profile.position_profile[3-c][3-seq[i]]);

						/*refmm = c;
						 readmm = seq[i];
						 if (e.all_mm_pos_counter < 10) {
						 e.all_mm_pos[e.all_mm_pos_counter] = i;
						 e.all_mm_ref_char[e.all_mm_pos_counter] = c;
						 fprintf(stderr, "error at ref %d and i=%d\n", c, i);
						 //e.all_mm_pos_counter++;*/
						//} else {
						//fprintf(stderr, "ERROR: TRYING TO SET mm pos counter > 10\n");
						//}
						//last_err = i;
						//fprintf(stderr, "mut from REF=%d to READ=%d found at position i=%d\n",c,seq[i],i);
						//} else {
						//last_err = e.last_err;
					}
					//gap_push_parasuite(stack, i, k, l, e.n_mm + is_mm, e.n_gapo, e.n_gape, e.n_ins, e.n_del, STATE_M, is_mm, opt, e.all_mm_pos, e.all_mm_ref_char, e.all_mm_pos_counter + is_mm, last_err);
					gap_push_parasuite(stack, i, k, l, e.n_mm + is_mm, e.n_gapo,
							e.n_gape, e.n_ins, e.n_del, STATE_M, is_mm, opt,
							new_ep_p_val, len);
				}
			}
		} else if (seq[i] < 4) { // try exact match only
			//fprintf(stderr, "try exact match at bottom at i=%d\n",i);
			int c = seq[i] & 3;
			k = bwt->L2[c] + cnt_k[c] + 1;
			l = bwt->L2[c] + cnt_l[c];
			if (k <= l)
				gap_push_parasuite(stack, i, k, l, e.n_mm, e.n_gapo, e.n_gape, e.n_ins,
						e.n_del, STATE_M, 0, opt, e.ep_p_val, len);
		}
	}

	*_n_aln = n_aln;
	//fprintf(stderr, "max_entries = %d\n", max_entries);
	return aln;
}
