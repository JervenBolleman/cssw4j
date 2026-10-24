/*
 * Throughput benchmark of the C SSW library on the protein test dataset, the counterpart of the Java
 * swiss.sib.cssw4j.Benchmark. Aligns all queries against all targets (BLOSUM50) `rounds` times.
 *
 * Usage: ssw_bench <targets.fa> <queries.fa> <gapO> <gapE> <flag> <rounds>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "ssw.h"

/* Same as ssw_reference.c */
static const int8_t blosum50[] = {
	 5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,
	-2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,
	-1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,
	-2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,
	-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,
	-1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,
	-1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,
	 0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,
	-2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,
	-1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,
	-2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,
	-1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,
	-1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,
	-3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,
	-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,
	 1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,
	 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5,
	-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5,
	-2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5,
	 0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5,
	-2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5,
	-1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5,
	-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
	-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1
};

static const char* AA = "ARNDCQEGHILKMFPSTWYVBZX*";

typedef struct { int8_t* seq; int32_t len; } record;

static record* read_fasta(const char* path, int32_t* count) {
	FILE* f = fopen(path, "r");
	if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
	int8_t table[128];
	for (int i = 0; i < 128; ++i) table[i] = 23;
	for (int i = 0; AA[i]; ++i) { table[(int)AA[i]] = i; table[(int)AA[i] + 32 * (AA[i] >= 'A')] = i; }
	int32_t cap = 64, n = 0;
	record* recs = malloc(cap * sizeof(record));
	char* line = NULL;
	size_t lcap = 0;
	ssize_t l;
	while ((l = getline(&line, &lcap, f)) >= 0) {
		while (l > 0 && (line[l - 1] == '\n' || line[l - 1] == '\r')) --l;
		if (line[0] == '>') {
			if (n == cap) { cap *= 2; recs = realloc(recs, cap * sizeof(record)); }
			recs[n].seq = NULL;
			recs[n].len = 0;
			++n;
		} else if (n > 0 && l > 0) {
			record* r = &recs[n - 1];
			r->seq = realloc(r->seq, r->len + l);
			for (ssize_t i = 0; i < l; ++i) r->seq[r->len + i] = table[line[i] & 127];
			r->len += l;
		}
	}
	free(line);
	fclose(f);
	*count = n;
	return recs;
}

int main(int argc, char* argv[]) {
	if (argc != 7) {
		fprintf(stderr, "Usage: %s <targets.fa> <queries.fa> <gapO> <gapE> <flag> <rounds>\n", argv[0]);
		return 1;
	}
	int32_t nt, nq;
	record* targets = read_fasta(argv[1], &nt);
	record* queries = read_fasta(argv[2], &nq);
	int32_t gapO = atoi(argv[3]), gapE = atoi(argv[4]), flag = atoi(argv[5]), rounds = atoi(argv[6]);
	long long cells = 0, checksum = 0;
	struct timespec t0, t1;
	clock_gettime(CLOCK_MONOTONIC, &t0);
	for (int32_t r = 0; r < rounds; ++r) {
		for (int32_t q = 0; q < nq; ++q) {
			s_profile* p = ssw_init(queries[q].seq, queries[q].len, blosum50, 24, 2);
			for (int32_t t = 0; t < nt; ++t) {
				s_align* a = ssw_align(p, targets[t].seq, targets[t].len, gapO, gapE, flag, 0, 0, queries[q].len / 2);
				checksum += a->score1;
				cells += (long long)queries[q].len * targets[t].len;
				align_destroy(a);
			}
			init_destroy(p);
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &t1);
	double s = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
	printf("C      flag=%d gap=%d/%d: %8.3f s, %6.2f GCUPS (checksum %lld)\n", flag, gapO, gapE, s, cells / s / 1e9, checksum);
	return 0;
}
