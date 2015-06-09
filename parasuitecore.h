/*
 * parasuitecore.h
 *
 *  Created on: 17.07.2014
 *      Author: akloetgen
 */

#ifndef PARASUITECORE_H_
#define PARASUITECORE_H_

#include <stdint.h>
#include "bwt.h"
#include "bwtaln.h"

#ifdef __cplusplus
extern "C" {
#endif

	gap_opt_t *gap_init_opt_parasuite();
	void parasuite_aln_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt);

	void parasuite_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt);
	int parasuite_cal_avgdiff(const gap_opt_t *opt, int length);

#ifdef __cplusplus
}
#endif

#endif /* parasuiteCORE_H_ */
