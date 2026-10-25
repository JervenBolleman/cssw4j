/*
 * Checks the native library libcssw4j.so (built with GraalVM Native Image) against the
 * reference results of the original C SSW library.
 *
 * Aligns every query against every target through the cssw4j_* C API and prints the first
 * eleven columns of the reference format: query, target, score1, score2, ref_begin1,
 * ref_end1, read_begin1, read_end1, ref_end2, flag, cigar. mark_mismatch is not part of the
 * native API, so the last two reference columns are not produced.
 *
 * Usage: ssw_native_check <protein|dna> <matrix file|blosum50|match:mismatch> <gapO> <gapE> <flag> <targets.fa> <queries.fa>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "libcssw4j.h"

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

static int8_t aa_table[128] = {
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

static int8_t nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct { char* name; char* seq; int32_t len; } record;

static record* read_fasta(const char* path, int32_t* count) {
	FILE* f = fopen(path, "r");
	if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
	int32_t cap = 64, n = 0;
	record* recs = malloc(cap * sizeof(record));
	char* line = NULL;
	size_t lcap = 0;
	ssize_t l;
	while ((l = getline(&line, &lcap, f)) >= 0) {
		while (l > 0 && (line[l - 1] == '\n' || line[l - 1] == '\r')) line[--l] = 0;
		if (line[0] == '>') {
			if (n == cap) { cap *= 2; recs = realloc(recs, cap * sizeof(record)); }
			char* end = strpbrk(line + 1, " \t");
			if (end) *end = 0;
			recs[n].name = strdup(line + 1);
			recs[n].seq = calloc(1, 1);
			recs[n].len = 0;
			++n;
		} else if (n > 0 && l > 0) {
			record* r = &recs[n - 1];
			r->seq = realloc(r->seq, r->len + l + 1);
			memcpy(r->seq + r->len, line, l + 1);
			r->len += l;
		}
	}
	free(line);
	fclose(f);
	*count = n;
	return recs;
}

/* Same parsing rules as main.c of the SSW library. */
static int8_t* read_matrix(const char* path, int32_t* n) {
	FILE* f = fopen(path, "r");
	if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
	char line[128];
	int8_t* mat = calloc(1024, 1);
	int32_t k = 0, m = 0;
	while (fgets(line, 128, f)) {
		if (line[0] == '*' || (line[0] >= 'A' && line[0] <= 'Z')) {
			if (line[0] >= 'A' && line[0] <= 'Z') aa_table[(int)line[0]] = aa_table[(int)line[0] + 32] = m;
			char str[4], *s = str;
			str[0] = '\0';
			for (int32_t l = 1; line[l]; ++l) {
				if ((line[l] >= '0' && line[l] <= '9') || line[l] == '-') *s++ = line[l];
				else if (str[0] != '\0') { *s = '\0'; mat[k++] = (int8_t)atoi(str); s = str; str[0] = '\0'; }
			}
			if (str[0] != '\0') { *s = '\0'; mat[k++] = (int8_t)atoi(str); }
			++m;
		}
	}
	fclose(f);
	*n = m;
	return mat;
}

int main(int argc, char* argv[]) {
	if (argc != 8) {
		fprintf(stderr, "Usage: %s <protein|dna> <matrix file|blosum50|match:mismatch> <gapO> <gapE> <flag> <targets.fa> <queries.fa>\n", argv[0]);
		return 1;
	}
	int protein = strcmp(argv[1], "protein") == 0;
	int32_t n;
	const int8_t* mat;
	int8_t* table;
	if (protein) {
		table = aa_table;
		if (strcmp(argv[2], "blosum50") == 0) { mat = blosum50; n = 24; }
		else mat = read_matrix(argv[2], &n);
	} else {
		int32_t match, mismatch;
		sscanf(argv[2], "%d:%d", &match, &mismatch);
		int8_t* m = calloc(25, 1);
		int32_t k = 0;
		for (int32_t i = 0; i < 4; ++i) {
			for (int32_t j = 0; j < 4; ++j) m[k++] = i == j ? match : -mismatch;
			m[k++] = 0;
		}
		mat = m;
		n = 5;
		table = nt_table;
	}
	int32_t gapO = atoi(argv[3]), gapE = atoi(argv[4]), flag = atoi(argv[5]);
	int32_t nt, nq;
	record* targets = read_fasta(argv[6], &nt);
	record* queries = read_fasta(argv[7], &nq);

	graal_isolate_t* isolate = NULL;
	graal_isolatethread_t* thread = NULL;
	if (graal_create_isolate(NULL, &isolate, &thread) != 0) {
		fprintf(stderr, "could not create the GraalVM isolate\n");
		return 1;
	}
	void* aligner = cssw4j_aligner_new(thread);
	int32_t cigarCap = 1 << 16;
	uint32_t* cigar = malloc(cigarCap * sizeof(uint32_t));

	for (int32_t q = 0; q < nq; ++q) {
		int8_t* qnum = malloc(queries[q].len + 1);
		for (int32_t i = 0; i < queries[q].len; ++i) qnum[i] = table[(int)queries[q].seq[i] & 127];
		void* profile = cssw4j_profile_new(thread, (char*)qnum, queries[q].len, (char*)mat, n, 2);
		int32_t maskLen = queries[q].len / 2;
		for (int32_t t = 0; t < nt; ++t) {
			int8_t* tnum = malloc(targets[t].len + 1);
			for (int32_t i = 0; i < targets[t].len; ++i) tnum[i] = table[(int)targets[t].seq[i] & 127];
			int32_t out[9]; /* score1 score2 refBegin1 refEnd1 readBegin1 readEnd1 refEnd2 flag cigarLen */
			int rc = cssw4j_align(thread, aligner, profile, (char*)tnum, targets[t].len, gapO, gapE, flag, 0, 0,
					maskLen, out, (int*)cigar, cigarCap);
			if (rc != 0) {
				fprintf(stderr, "cssw4j_align failed with %d for %s vs %s\n", rc, queries[q].name, targets[t].name);
				return 1;
			}
			printf("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", queries[q].name, targets[t].name,
				out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7]);
			if (out[8] == 0) printf("*");
			for (int32_t c = 0; c < out[8]; ++c) printf("%u%c", cigar[c] >> 4, "MIDNSHP=X"[(cigar[c] & 0xf) > 8 ? 0 : cigar[c] & 0xf]);
			printf("\n");
			free(tnum);
		}
		cssw4j_profile_free(thread, profile);
		free(qnum);
	}
	cssw4j_aligner_free(thread, aligner);
	graal_tear_down_isolate(thread);
	free(cigar);
	return 0;
}
