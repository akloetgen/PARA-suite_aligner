/*
 * parmacore.h
 *
 *  Created on: 17.07.2014
 *      Author: akloetgen
 */

#ifndef PARMACORE_H_
#define PARMACORE_H_

#include <stdint.h>
#include "bwt.h"
#include "bwtaln.h"

#ifdef __cplusplus
extern "C" {
#endif

	gap_opt_t *gap_init_opt_parma();
	void parma_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt);

	void parma_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);
	int parma_cal_avgdiff(int length, double avg_err);

#ifdef __cplusplus
}
#endif

#endif /* PARMACORE_H_ */
