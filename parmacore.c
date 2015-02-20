/*
 * parmacore.c
 *
 *  Created on: 17.07.2014
 *      Author: akloetgen
 */

#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "parmacore.h"
#include "parma.h"
#include "utils.h"
#include "bwa.h"

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifndef LF
#define LF	10
#endif

int load_indel_profile(gap_opt_t *opt, const char *indel_filename) {
	// load indel rates from file.
	fprintf(stderr, "[%s] Load indel profile from file %s\n", __func__,
			indel_filename);

	int line_length = 255;
	char line[line_length];
	double insertion, deletion;
	FILE *indel_file = fopen(indel_filename, "r");

	if (indel_file != NULL) {
		while (fgets(line, line_length, indel_file)) {
			sscanf(line, "%lf\t%lf", &insertion, &deletion);
			fprintf(stderr, "[%s] insertion rate=%lf\n", __func__, insertion);
			fprintf(stderr, "[%s] deletion rate=%lf\n", __func__, deletion);
			//opt->profile.indels
			opt->profile.indels[INSERTION] = insertion;
			opt->profile.indels[DELETION] = deletion;
		}
		fclose(indel_file);
	} else {
		fprintf(stderr,
				"[%s] ABORT: Couldn't read indel profile from file \"%s\". Probably the"
						" file is missing or rights are not set correctly...\n",
				__func__, indel_filename);
		return 1;
	}

	return 0;
}

int load_error_profile(gap_opt_t *opt, const char *ep_filename) {
	// load error profile from file if already exists, otherwise ep should be calculated on-the-fly
	// should always be loaded from file as java is the main class??!?!

	fprintf(stderr, "[%s] Load error profile from file %s\n", __func__,
			ep_filename);

	int line_length = 255;
	char line[line_length];
	double A, C, T, G;
	FILE *ep_file = fopen(ep_filename, "r");

	if (ep_file != NULL) {
		int k = 0;
		while (fgets(line, line_length, ep_file)) {
			sscanf(line, "%lf\t%lf\t%lf\t%lf", &A, &C, &G, &T);
			fprintf(stderr, "[%s] A=%f; C=%f, G=%f, T=%f\n", __func__, A, C, G,
					T);

			opt->profile.position_profile[k][0] = A;
			opt->profile.position_profile[k][1] = C;
			opt->profile.position_profile[k][2] = G;
			opt->profile.position_profile[k][3] = T;
			k++;
			//fprintf(stderr, "A=%Lf; C=%Lf, G=%Lf, T%Lf\n", opt->profile.position_profile[k][0],opt->profile.position_profile[k][1],opt->profile.position_profile[k][2],opt->profile.position_profile[k][3]);
		}
		fclose(ep_file);
	} else {
		fprintf(stderr,
				"[%s] ABORT: Couldn't read error profile from file \"%s\". Probably the"
						" file is missing or rights are not set correctly...\n",
				__func__, ep_filename);
		//exit(EXIT_FAILURE);
		return 1;
	}

	return 0;
}

int calculateAverageProbabilites(gap_opt_t *opt, int including_indels) {
	// CALCULATE MAX THRESHOLD FOR P-VAL DURING BACKWARD-SEARCH
	int i, j;
	double avg_match = 0.0;
	double avg_mm = 0.0;
	double best_mm = 0.0;
	if (including_indels) {
		double indel = opt->profile.indels[INSERTION]
				+ opt->profile.indels[DELETION];
		double bias = ((double) 1.0 - indel);

		/*	for (i = 0; i < 4; i++) {
		 for (j = 0; j < 4; j++) {
		 if (i == j) {
		 avg_match += (double) opt->profile.position_profile[i][j]
		 * bias;
		 } else {
		 avg_mm += (double) opt->profile.position_profile[i][j]
		 * bias;
		 if ((double) opt->profile.position_profile[i][j] * bias
		 > best_mm) {
		 best_mm = (double) opt->profile.position_profile[i][j]
		 * bias;
		 }
		 }
		 }

		 }*/
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				if (i == j) {
					avg_match += (double) opt->profile.position_profile[i][j];
				} else {
					avg_mm += (double) opt->profile.position_profile[i][j];
					if ((double) opt->profile.position_profile[i][j]
							> best_mm) {
						best_mm = (double) opt->profile.position_profile[i][j];
					}
				}
			}
		}
		avg_mm += indel;
		avg_match += bias;
		opt->avg_match = (double) avg_match / 5;
		opt->avg_mm = (double) avg_mm / 14;
		//opt->avg_mm = (double) avg_mm / 16;
		opt->best_mm = (double) best_mm;
		//fprintf(stderr, "avg_match=%f; avg_mm=%f; best_mm=%f\n", opt->avg_match,
		//		opt->avg_mm, opt->best_mm);
	} else {
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				if (i == j) {
					avg_match += opt->profile.position_profile[i][j];
				} else {
					avg_mm += opt->profile.position_profile[i][j];
					if (opt->profile.position_profile[i][j] > best_mm) {
						best_mm = opt->profile.position_profile[i][j];
					}
				}
			}
		}
		opt->avg_match = (double) avg_match / 4;
		opt->avg_mm = (double) avg_mm / 12;
		opt->best_mm = (double) best_mm;
		//fprintf(stderr, "avg_match=%f; avg_mm=%f; best_mm=%f", opt->avg_match,
		//	opt->avg_mm, opt->best_mm);
	}

	return 0;
}

gap_opt_t *gap_init_opt_parma() {
	gap_opt_t *o;
	o = (gap_opt_t*) calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	 rate. Voilating this requirement will break pairing! */
	o->s_mm = 3;
	o->s_gapo = 11;
	o->s_gape = 4;
	// DISALLOW GAPS!!! MUCH FASTER AND DOESNT HAVE ANY INFLUENCE ON THE ALREADY SHORT READS OF <26BP!!! OTHERWISE ACTIVATE GAP BY SETTING GAPO = 1
	o->max_diff = -1;
	o->max_gapo = 1;
	o->max_gape = -1;
	o->indel_end_skip = 5;
	o->max_del_occ = 10;
	o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32;
	o->max_seed_diff = 2;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			if (i != j) {
				o->profile.position_profile[i][j] = 0.01;
			} else {
				o->profile.position_profile[i][j] = 0.97;
			}
		}
	}
	/*o->profile.indels[INSERTION] = 0.00025;
	 o->profile.indels[DELETION] = 0.0065;*/
	/*o->profile.indels[INSERTION] = 0.000075;
	 o->profile.indels[DELETION] = 0.000075;*/
	o->profile.indels[INSERTION] = 0.001658;
	o->profile.indels[DELETION] = 0.001899;
	o->X = -1;
	return o;
}

int parma_cal_avgdiff(const gap_opt_t *opt, int length) {
	int temp_allowed_mm;
	//if ((temp_allowed_mm = floor((double) length * avg_err * 10)) > 2) {
	if ((temp_allowed_mm = floor((double) length * (opt->avg_mm * 14))) > 2) {
		/*if ((temp_allowed_mm = floor(
		 (double) (4 - opt->profile.position_profile[0][0]
		 - opt->profile.position_profile[1][1]
		 - opt->profile.position_profile[2][2]
		 - opt->profile.position_profile[3][3]) * length)) > 2) {*/
		return temp_allowed_mm;
	}
	return 2;
}

// width must be filled as zero
int bwt_cal_width_parma(const bwt_t *bwt, int len, const ubyte_t *str,
		bwt_width_t *width) {
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0;
	l = bwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = bwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

void parma_cal_sa_reg_gap(int tid, bwt_t * const bwt, int n_seqs,
		bwa_seq_t *seqs, const gap_opt_t *opt) {
	//fprintf(stderr, "bla\n");
	int i, j, max_l = 0, max_len, avg_len = 0;
	gap_stack_t_parma *stack;
	bwt_width_t *w, *seed_w;
	gap_opt_t local_opt = *opt;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i) {
		if (seqs[i].len > max_len) {
			max_len = seqs[i].len;
		}
		avg_len += seqs[i].len;
	}
	avg_len = (int) avg_len / n_seqs;
	// FÜR SPEEDUP KÖNNTE SO ETWAS ÄHNLICHES WIEDER REIN; FÜR DEN CORE ALG IST DAS NICHT NOTWENDIG!

	//if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff_ep(max_len, BWA_AVG_ERR, opt->fnr);
	if (opt->X == -1)
		local_opt.X = parma_cal_avgdiff(opt, avg_len);
	//fprintf(stderr, "[%s] for read lenght %d, average error X = %d\n", __func__, avg_len, local_opt.X);

	if (local_opt.X < local_opt.max_gapo)
		local_opt.max_gapo = local_opt.max_diff;
	//int min_init_value = -1;
	stack = gap_init_stack_parma(local_opt.X + 1, local_opt.max_gapo,
			local_opt.max_gape, &local_opt);
	//fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

	seed_w = (bwt_width_t*) calloc(opt->seed_len + 1, sizeof(bwt_width_t));
	w = 0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
#ifdef HAVE_PTHREAD
		if (i % opt->n_threads != tid) continue;
#endif
		p->sa = 0;
		p->type = BWA_TYPE_NO_MATCH;
		p->c1 = p->c2 = 0;
		p->n_aln = 0;
		p->aln = 0;
		if (max_l < p->len) {
			max_l = p->len;
			w = (bwt_width_t*) realloc(w, (max_l + 1) * sizeof(bwt_width_t));
			memset(w, 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width_parma(bwt, p->len, p->seq, w);

		// HIER DAS SELBE WIE OBEN: ERSETZEN DURCH EP?!?!
		//if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff_ep(p->len, BWA_AVG_ERR, opt->fnr);
		if (opt->X == -1)
			local_opt.X = parma_cal_avgdiff(opt, avg_len);
		local_opt.seed_len =
				opt->seed_len < p->len ? opt->seed_len : 0x7fffffff;
		if (p->len > opt->seed_len)
			bwt_cal_width_parma(bwt, opt->seed_len,
					p->seq + (p->len - opt->seed_len), seed_w);
		// core function
		for (j = 0; j < p->len; ++j) // we need to complement
			p->seq[j] = p->seq[j] > 3 ? 4 : 3 - p->seq[j];

		p->aln = bwt_match_gap_parma(bwt, p->len, p->seq, w,
				p->len <= opt->seed_len ? 0 : seed_w, &local_opt, &p->n_aln,
				stack);

		//fprintf(stderr, "read %d has %i mms\n", p->name, p->aln->n_mm);
		//fprintf(stderr, "mm=%lld,ins=%lld,del=%lld,gapo=%lld\n", p->aln->n_mm, p->aln->n_ins, p->aln->n_del, p->aln->n_gapo);
		// clean up the unused data in the record
		free(p->name);
		free(p->seq);
		free(p->rseq);
		free(p->qual);
		p->name = 0;
		p->seq = p->rseq = p->qual = 0;
	}
	free(seed_w);
	free(w);
	gap_destroy_stack_parma(stack);
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt;
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
}thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	parma_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif

bwa_seqio_t *parma_open_reads(int mode, const char *fn_fa) {
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE)
			which |= 4;
		if (mode & BWA_MODE_BAM_READ1)
			which |= 1;
		if (mode & BWA_MODE_BAM_READ2)
			which |= 2;
		if (which == 0)
			which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else
		ks = bwa_seq_open(fn_fa);
	return ks;
}

void bwa_alnep_core(const char *prefix, const char *fn_fa, const gap_opt_t *opt) {
	int i, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bwt_t *bwt;

	// initialization
	ks = parma_open_reads(opt->mode, fn_fa);

	{ // load BWT
		char *str = (char*) calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix);
		strcat(str, ".bwt");
		bwt = bwt_restore_bwt(str);
		free(str);
	}

	// core loop
	err_fwrite(SAI_MAGIC, 1, 4, stdout);
	err_fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	//while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual))
			!= 0) {
		tot_seqs += n_seqs;
		t = clock();

		fprintf(stderr, "[parma_core] calculate SA coordinate... ");

#ifdef HAVE_PTHREAD
		if (opt->n_threads <= 1) { // no multi-threading at all
			//fprintf(stderr, "\nHAVE_PTHREAD defined: single threading\n");
			parma_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
		} else {
			//fprintf(stderr, "\nHAVE_PTHREAD defined: multi threading\n");
			pthread_t *tid;
			pthread_attr_t attr;
			thread_aux_t *data;
			int j;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
			for (j = 0; j < opt->n_threads; ++j) {
				data[j].tid = j; data[j].bwt = bwt;
				data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = opt;
				pthread_create(&tid[j], &attr, worker, data + j);
			}
			for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
			free(data); free(tid);
		}
#else
		//fprintf(stderr, "\nHAVE_PTHREAD not defined.\n");
		parma_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

		fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		t = clock();
		fprintf(stderr, "[parma_core] write to the disk... ");
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			err_fwrite(&p->n_aln, 4, 1, stdout);
			//fprintf(stderr, "\nk=%lu; l=%lu\n",p->aln->k,p->aln->k);
			if (p->n_aln)
				err_fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
		}
		fprintf(stderr, "%.2f sec\n", (float) (clock() - t) / CLOCKS_PER_SEC);
		t = clock();

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[parma_core] %d sequences have been processed.\n",
				tot_seqs);
	}

	// destroy
	bwt_destroy(bwt);
	bwa_seq_close(ks);
}

int bwa_parma(int argc, char *argv[]) {
	int c, opte = -1;
	gap_opt_t *opt;
	char *prefix;
	int error_profile_set = 0, indel_profile_set = 0;

	opt = gap_init_opt_parma();
	while ((c = getopt(argc, argv,
			"n:p:g:o:e:i:d:l:k:LR:m:t:NM:O:E:q:f:b012YB:I:D:")) >= 0) {
		switch (c) {
		case 'n':
			opt->X = atoi(optarg);
			break;
			/*case 'n':
			 if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			 else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			 break;*/
		case 'p':
			if (load_error_profile(opt, optarg)) {
				free(opt);
				return 1;
			}
			error_profile_set = 1;
			break;
		case 'g':
			if (load_indel_profile(opt, optarg)) {
				free(opt);
				return 1;
			}
			indel_profile_set = 1;
			break;
		case 'o':
			opt->max_gapo = atoi(optarg);
			break;
		case 'e':
			opte = atoi(optarg);
			break;
			//case 'M': opt->s_mm = atoi(optarg); break;
		case 'O':
			opt->s_gapo = atoi(optarg);
			break;
		case 'E':
			opt->s_gape = atoi(optarg);
			break;
		case 'd':
			opt->max_del_occ = atoi(optarg);
			break;
		case 'i':
			opt->indel_end_skip = atoi(optarg);
			break;
		case 'l':
			opt->seed_len = atoi(optarg);
			break;
		case 'k':
			opt->max_seed_diff = atoi(optarg);
			break;
		case 'm':
			opt->max_entries = atoi(optarg);
			break;
		case 't':
			opt->n_threads = atoi(optarg);
			break;
		case 'L':
			opt->mode |= BWA_MODE_LOGGAP;
			break;
		case 'R':
			opt->max_top2 = atoi(optarg);
			break;
		case 'q':
			opt->trim_qual = atoi(optarg);
			break;
		case 'N':
			opt->mode |= BWA_MODE_NONSTOP;
			opt->max_top2 = 0x7fffffff;
			break;
		case 'f':
			xreopen(optarg, "wb", stdout);
			break;
		case 'b':
			opt->mode |= BWA_MODE_BAM;
			break;
		case '0':
			opt->mode |= BWA_MODE_BAM_SE;
			break;
		case '1':
			opt->mode |= BWA_MODE_BAM_READ1;
			break;
		case '2':
			opt->mode |= BWA_MODE_BAM_READ2;
			break;
			//case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y':
			opt->mode |= BWA_MODE_CFY;
			break;
		case 'B':
			opt->mode |= atoi(optarg) << 24;
			break;
		case 'I':
			opt->profile.indels[INSERTION] = atof(optarg);
			break;
		case 'D':
			opt->profile.indels[DELETION] = atof(optarg);
			break;
		default:
			return 1;
		}
	}
	if (error_profile_set && indel_profile_set) {
		fprintf(stderr,
				"[%s] calculate average values for matches and mismatches/indels\n",
				__func__);
		calculateAverageProbabilites(opt, 1);
	} else if (error_profile_set) {
		fprintf(stderr,
				"[%s] calculate average values for matches and mismatches (excluding indels)\n",
				__func__);
		calculateAverageProbabilites(opt, 0);
	}

	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	/*	int q = 15;
	 for (q = 15; q < 40; q++) {
	 fprintf(stderr, "length=%d; X=%d\n", q, parma_cal_avgdiff(opt, q));
	 }*/

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Program: PARMA addon for bwa version 0.7.8\n");
		fprintf(stderr, "Version: 0.5 alpha\n");
		fprintf(stderr,
				"Contact: Andreas Kloetgen <andreas.kloetgen@hhu.de>\n\n");
		fprintf(stderr, "Usage:   bwa parma [options] <prefix> <in.fq>\n\n");
		fprintf(stderr,
				"Options: -n NUM    median #diff (int). Real #diff depends on error profile. [%.d]\n",
				opt->X);
		/*fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
		 BWA_AVG_ERR, opt->fnr);*/
		fprintf(stderr, "         -p FILE   error profile saved to a file\n");
		fprintf(stderr,
				"         -o INT    maximum number or fraction of gap opens [%d]\n",
				opt->max_gapo);
		fprintf(stderr,
				"         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr,
				"         -i INT    do not put an indel within INT bp towards the ends [%d]\n",
				opt->indel_end_skip);
		fprintf(stderr,
				"         -d INT    maximum occurrences for extending a long deletion [%d]\n",
				opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr,
				"         -k INT    maximum differences in the seed [%d]\n",
				opt->max_seed_diff);
		fprintf(stderr,
				"         -m INT    maximum entries in the queue [%d]\n",
				opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n",
				opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n",
				opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n",
				opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n",
				opt->s_gape);
		fprintf(stderr,
				"         -R INT    stop searching when there are >INT equally best hits [%d]\n",
				opt->max_top2);
		fprintf(stderr,
				"         -q INT    quality threshold for read trimming down to %dbp [%d]\n",
				BWA_MIN_RDLEN, opt->trim_qual);
		fprintf(stderr,
				"         -f FILE   file to write output to instead of stdout\n");
		fprintf(stderr, "         -B INT    length of barcode\n");
		fprintf(stderr,
				"         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr,
				"         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		//fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		fprintf(stderr,
				"         -b        the input read file is in the BAM format\n");
		fprintf(stderr,
				"         -0        use single-end reads only (effective with -b)\n");
		fprintf(stderr,
				"         -1        use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr,
				"         -2        use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr,
				"         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "         -I        insertion rate [%f]\n",
				opt->profile.indels[INSERTION]);
		fprintf(stderr, "         -D        insertion rate [%f]\n",
				opt->profile.indels[DELETION]);
		fprintf(stderr, "\n");
		return 1;
	}

	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index\n", __func__);
		free(opt);
		return 1;
	}
	bwa_alnep_core(prefix, argv[optind + 1], opt);
	free(opt);
	free(prefix);
	return 0;
}
