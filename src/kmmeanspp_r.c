/**
 * @file kmmeanspp_r.c
 *
 * Purpose: C wrapper for km-kmeans++.
 */

#include <R.h>
#include <Rinternals.h>

#include "kmmeanspp.h"
#include "array.h"

enum RETURN_SLOT {
	WITHIN_SS,
	CENTERS,
	PARTITION,
	NUMBER_SLOTS,
};

char const *return_names[NUMBER_SLOTS] = {
	"criterion",
	"centers",
	"partition"
};


/**
 * This function is called from R.
 *
 * @param data_r	matrix of numeric data (REALSXP pxn matrix)
 * @param K_r		number of clusters (INTSXP scalar)
 * @param ninit_r	number of initializations (INTSXP scalar)
 * @return		list
 */
SEXP kmmeanspp_r(	SEXP data_r,
			SEXP K_r,
			SEXP ninit_r,
			SEXP kmmnsiter_r
	)
{
	int K, n_init, kmmns_iter, nprotect = 0, n, p, **data1, *iptr;
	double *dptr, **data, **data2;
	SEXP data_c = NULL;
	SEXP dim = getAttrib(data_r, R_DimSymbol);
	SEXP return_list = R_NilValue;
	SEXP return_list_names;

	if (!isMatrix(data_r)) {
		Rprintf("Data should be a matrix.\n");
		return R_NilValue;
	}

	if (!isReal(data_r)) {
		PROTECT(data_c = coerceVector(data_r, REALSXP));
		++nprotect;
		dptr = REAL(data_c);
	} else {
		dptr = REAL(data_r);
	}

	if (!isInteger(K_r)) {
		Rprintf("Invalid type for number of clusters.\n");
		return R_NilValue;
	}
	if (!isInteger(ninit_r)) {
		Rprintf("Invalid type for number of initializations.\n");
		return R_NilValue;
	}
	if (!isInteger(kmmnsiter_r)) {
		Rprintf("Invalid type for number of kmmeans iterations.\n");
		return R_NilValue;
 	}

	GetRNGstate();

	K = * INTEGER(K_r);
	n_init = * INTEGER(ninit_r);
 	kmmns_iter = * INTEGER(kmmnsiter_r);
	n = INTEGER(dim)[0];
	p = INTEGER(dim)[1];

	MAKE_MATRIX(data, n, p);
	MAKE_MATRIX(data1, n, p + 3);
	MAKE_MATRIX(data2, K, 2*p + 3);

	for (int i = 0; i < n; ++i)
		for (int j = 0; j < p; ++j) {
			int l = j * n + i;

			/* we intentionally represent missing as infinite */
			if (ISNA(dptr[l]))
				data[i][j] = INFINITY;
			else
				data[i][j] = dptr[l];
		}
	
	double tss = repkmmeanspp(data, K, n, p, data1, data2, n_init, kmmns_iter);
	PROTECT(return_list = allocVector(VECSXP, NUMBER_SLOTS));
	++nprotect;

	SET_VECTOR_ELT(return_list, WITHIN_SS, ScalarReal(tss));
	SET_VECTOR_ELT(return_list, CENTERS, allocMatrix(REALSXP, K, p));

	dptr = REAL(VECTOR_ELT(return_list, CENTERS));
	for (int k = 0; k < K; ++k)
		for (int j = 0; j < p; ++j)
			dptr[K*j + k] = data2[k][j];

	SET_VECTOR_ELT(return_list, PARTITION, allocVector(INTSXP, n));
	iptr = INTEGER(VECTOR_ELT(return_list, PARTITION));
	for (int i = 0; i < n; ++i)
		iptr[i] = data1[i][1];	/* why the hell? */

	FREE_MATRIX(data);
	FREE_MATRIX(data1);
	FREE_MATRIX(data2);

	PROTECT(return_list_names = allocVector(STRSXP, NUMBER_SLOTS));
	++nprotect;

	SET_STRING_ELT(return_list_names, WITHIN_SS, mkChar(return_names[WITHIN_SS]));
	SET_STRING_ELT(return_list_names, CENTERS, mkChar(return_names[CENTERS]));
	SET_STRING_ELT(return_list_names, PARTITION, mkChar(return_names[PARTITION]));

	setAttrib(return_list, R_NamesSymbol, return_list_names);

	PutRNGstate();

	UNPROTECT(nprotect);
	return return_list;
} /* kmmeanspp_r */
