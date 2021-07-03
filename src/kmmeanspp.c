/*Andy Lithio, 5-9-2013, Functions for K-Means Algorithm, 580 Semester Project, See Report*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <R_ext/Print.h> 
#include <R/Rmath.h>
#include <time.h>

#include "array.h"

void makeIV(int n, int m, double **dat, double **iv);
void sumsquares(double **clust, double **dat, int **obs, int k, int n, int m, double *tss, double **iv);
void clustmeans(double **clust, double **iv, int **obs, int k, int n, int m);
void kplusplus(double **clust, double **iv, int k, int n, int m);
void qtrans(double **clust, double **dat, int **obs, double **iv, int k, int n, int m,int *ncp,int *itran,int *indx,double *d);
void otrans(double **clust, double **dat, int **obs, double **iv, int k, int n, int m,int *ncp,int *itran,int *live,int *indx,double *d);
int kmeans2(double **dat,double **clust, int **obs, double *tss, int k, int m, int n, int numiter, int kmnsiter);

/*
  Calculates values for m-by-3n inner variance matrix as described in report. First n columns have IV value for each dimension of each observation (n_ij*sum(y_ij^2)-sum(y_ij)^2), the second n columns have the sum of recorded values in each dimension for each observation (sum(y_ij)), and the third n columns have sum of recorded values squared (sum(y_ij^2))

  Variables
  n: # of dimensions
  m: # of measurements
  dat: matrix of observed data
  iv: inner variance matrix

  Working Variables
  i: row in dat
  j: dimension
  jj: measurement in curobs
  n_i: number of measurements for current observation
  n_ij: number of recorded values for current dimension of current observation
*/
void makeIV(int n, int m, double **dat, double **iv)
{
	int i,j;

	for (i=0; i<m; i++){
		for (j=0; j<n; j++){
			if (isfinite(dat[i][j])){	/* non-missing */
				iv[i][j] = dat[i][j];
				iv[i][j+n] = dat[i][j]*dat[i][j];
			}else{
				iv[i][j] = NAN;
				iv[i][j+n] = NAN;
			}
		}
	}
}

/*
  Updates within cluster sum of squares (WSS) and total sum of squares (TSS). Only called at the end.

  Variables:
  **clust: matrix of cluster centers, etc
  **dat: matrix of observed data
  **obs: matrix of observations
  k: number of clusters
  n: number of dimensions
  m: total number of measurements
  *tss: pointer to total sum of squares
  **iv: inner variance matrix
  */
void sumsquares(double **clust, double **dat, int **obs, int k, int n, int m, double *tss, double **iv)
{
	int s, j;

	/*Sets all WSS to 0*/
	for (s=0; s<k; s++)
		clust[s][2*n+2] = 0;
  
	/*Calculates WSS and stores in clust*/
	for (s=0; s<m; s++)
		for (j=0; j<n; j++)
			if (obs[s][j+3] > 0)
				clust[obs[s][1]][2*n+2] = clust[obs[s][1]][2*n+2] + iv[s][j+n] + clust[obs[s][1]][j]*clust[obs[s][1]][j] - 2.0*clust[obs[s][1]][j]*iv[s][j];

	/*Adds up all WSS and stores in *tss*/
	*tss = 0.;
	for (s=0; s<k; s++)
		*tss = *tss + clust[s][2*n+2];
}


/*
  Calculates Initial Cluster Means, Sizes

  Variables:
  **clust: matrix of cluster centers, etc
  **iv: inner variance matrix
  **obs: matrix of observations
  k: number of clusters
  n: number of dimensions
  m: total number of measurements

  *totals: k-vector for holding sums of observed values
  *cnts: k-vector for holding counts of observed values
  ci: index for cluster that observation belongs to
*/
void clustmeans(double **clust, double **iv, int **obs, int k, int n, int m)
{
	double *totals;
	int *cnts,q,r,s,ci;

	MAKE_VECTOR(totals, k);
	MAKE_VECTOR(cnts, k);
	/*Do the following for each dimension*/
	for (r=0; r<n; r++){
		/*Initialize working vectors to 0*/
		for (q=0; q<k; q++){
			totals[q] = 0.;
			cnts[q] = 0;
		}
		/*For each measurement, if we have a value add it to the total of the cluster it belongs to and increment the count vector*/
		for (s=0; s<m; s++){
			if (obs[s][3+r] > 0){
				ci = obs[s][1];
				totals[ci] = totals[ci] + iv[s][r];
				cnts[ci] = cnts[ci] + 1;
			}
		}
		/*Set the cluster centers to the mean, num of meas per dim, and reset sizes to 0*/
		for (q=0; q<k; q++){
			clust[q][n+r] = (double) cnts[q];
			if (cnts[q] != 0)
				clust[q][r] = totals[q] / (double) cnts[q];
			else
				clust[q][r] = NAN;
			if (r == 0)
				clust[q][2*n] = 0;
		}
	}

	/*In a new loop, for each observation increment the ns (number of obs in each cluster)*/
	for (s=0; s<m; s++)
		clust[obs[s][1]][2*n] = clust[obs[s][1]][2*n] + 1;

	FREE_VECTOR(totals);
	FREE_VECTOR(cnts);
}

/*
  Optimal Transfer Step (OTRANS)

  Variables:
  **clust: matrix of cluster centers, etc
  **dat: matrix of observed data
  **obs: matrix of observations
  **iv: matrix of observation summaries
  k: number of clusters
  n: number of dimensions
  m: total number of measurements
  *ncp:
  *itran:
  *live: 
  *indx:
  *d: vector of reduction in SS if obs is removed from current cluster
  Working Variables:
  curobs: Index of current observation
  curni: n_i of current observation

*/
void otrans(double **clust, double **dat, int **obs, double **iv, int k, int n, int m,int *ncp,int *itran,int *live,int *indx,double *d)
{
	int j, r, l1, l2, minmeas, jj, flag;
	double de, df, da, db, r2, dd, dc, alw, alt, al1, al2;

	/*If transfer in last step, mark cluster as live*/
	for (r=0; r<k; r++)
		if (itran[r] == 1)
			live[r] = m+1;
  
	for (r=0; r<m; r++){
		(*indx) = (*indx) + 1 ;
		l1 = obs[r][1];
		l2 = obs[r][2];
		/*Only consider change if not only obs in clust*/
		minmeas = 1;
		if (((int) clust[l1][2*n]) == 1)
			minmeas = 0;
		if (minmeas == 1){
			/*Refresh reduction in SS if removed from cluster if necessary*/
			if (ncp[l1] != 0) {
				de = 0.;
				for (j=0; j<n; j++){
					if (obs[r][j+3] > 0 && clust[l1][j+n] > 1){
						df = iv[r][j + n] + clust[l1][j]*clust[l1][j] - 2.0*clust[l1][j]*iv[r][j];
						de = de + df*clust[l1][j+n]/(clust[l1][j+n] - (double) obs[r][j+3]);
					}
				}
				d[r] = de;
			}
			/*Calc increase in SS if moved to Clust2*/
			da = 0.;
			for (j=0; j<n; j++){
				if (obs[r][j+3] > 0 && clust[l2][j+n] > 0){
					db =  iv[r][j + n] + clust[l2][j]*clust[l2][j] - 2.0*clust[l2][j]*iv[r][j];
					da = da +db*clust[l2][j+n]/(clust[l2][j+n] + (double) obs[r][j+3]);
				}
			}
			r2 = da;
			/*Any other clusters better than L1 or L2?*/
			for (j=0; j<k; j++){
				if ((r >= live[l1]-1 && r >= live[j]-1) || j == l1 || j == l2) {
				}else{
					dc = 0.;
					flag = 0;
					jj = 0;
					while (flag == 0 && jj < n){
						if (obs[r][jj+3] > 0 && clust[j][n+jj] > 0){
							dd =  iv[r][jj + n] + clust[j][jj]*clust[j][jj] - 2.0*clust[j][jj]*iv[r][jj];
							dc = dc + dd*clust[j][n+jj]/(clust[j][n+jj] + (double) obs[r][jj+3]);
						}
						if (dc >= r2)
							flag = 1;
						jj = jj+1;
					}
					if (flag == 0){
						r2 = dc;
						l2 = j;
					}
				}
			}
			if (r2 >= d[r]) {
				obs[r][2] = l2;
			}else{
				(*indx)=0;
				live[l1] = m+r+1;
				live[l2] = m+r+1;
				ncp[l1] = r+1;
				ncp[l2] = r+1;
				clust[l1][2*n] = clust[l1][2*n] - 1;
				clust[l2][2*n] = clust[l2][2*n] + 1;
				for(j=0; j<n; j++){
					if (obs[r][j+3] > 0){
						al1 = clust[l1][n+j];
						alw = al1 - 1.;
						al2 = clust[l2][n+j];
						alt = al2+1.;
						if (alw > 0)
							clust[l1][j] = (clust[l1][j]*al1 - iv[r][j]) / alw;
						else
							clust[l1][j] = NAN;
						if (al2 > 0)
							clust[l2][j] = (clust[l2][j]*al2 + iv[r][j])/alt;
						else
							clust[l2][j] = iv[r][j];
						clust[l1][n+j] = alw;
						clust[l2][n+j] = alt;
					}
				}

				obs[r][1] = l2;
				obs[r][2] = l1;
			}
		}
		if (*indx == m)
			return;
	}
	for (j=0; j<k; j++){
		itran[j] = 0;
		live[j] = live[j] - m;
	}
	return;
}

/*
  Quick Transfer Step (QTRANS)

  Variables:
  **clust: matrix of cluster centers, etc
  **dat: matrix of observed data
  **obs: matrix of observations
  k: number of clusters
  n: number of dimensions
  m: total number of measurements
  *ncp:
  *itran
  *indx
  *d:

  Working Variables:
  curobs: Index of current observation
  curni: n_i of current observation
*/
void qtrans(double **clust, double **dat, int **obs, double **iv, int k, int n, int m,int *ncp,int *itran,int *indx,double *d)
{
	int r,icoun,istep,iflag,q,minmeas,j,flag,l1,l2;
	double da,db,dd,de,al1,alw,al2,alt,r2;

	icoun = 0;
	istep = 0;
	iflag = 0;

	while (iflag == 0){
		for (r=0; r<m; r++){
			icoun = icoun + 1;
			istep = istep + 1;
			l1 = obs[r][1];
			l2 = obs[r][2];
			/*Only consider change if not only obs in clust*/
			minmeas = 1;
			for (q=0; q<n; q++)
				if (((int) clust[l1][n+q]) == obs[r][3+q])
					minmeas = 0;
			if (minmeas == 1){
				/*If clust has been changed in last m steps, refresh distance*/
				if (istep <= ncp[l1]) {
					da = 0.;
					for (j=0; j<n; j++){
						if (obs[r][j+3] > 0 && clust[l1][j+n] > 1){
							db =  iv[r][j + n] + clust[l1][j]*clust[l1][j] - 2.0*clust[l1][j]*iv[r][j];
							da = da + db*clust[l1][j+n]/(clust[l1][j+n] - ((double) obs[r][j+3]));
						}
					}
					d[r] = da;
				}
				/*If either cluster is live, check to see if we should move point*/
				if (istep<ncp[l1] || istep < ncp[l2]){
					r2 = d[r];
					dd = 0.;
					flag = 0;
					j = 0;
					while (flag == 0 && j < n){
						if (obs[r][j+3] > 0 && clust[l2][n+j] > 0){
							de =  iv[r][j + n] + clust[l2][j]*clust[l2][j] - 2.0*clust[l2][j]*iv[r][j];
							dd = dd + de*clust[l2][n+j]/(clust[l2][n+j] + ((double) obs[r][j+3]));
						}
						if (dd >= r2)
							flag = 1;
						j++;
					}
					/*If reduction in SS greather than increase, reset count, indx, and ncp, mark itran, recalc means*/
					if (flag == 0){
						icoun = 0;
						(*indx) = 0;
						itran[l1] = 1;
						itran[l2] = 1;
						ncp[l1] = istep + m;
						ncp[l2] = istep + m;
						for (j=0; j<n; j++){
							if (obs[r][j+3] > 0){
								al1 = clust[l1][n+j];
								alw = al1 - 1;
								al2 = clust[l2][n+j];
								alt = al2+1;
								if (alw > 0)
									clust[l1][j] = (clust[l1][j]*al1 - iv[r][j])/alw;
								else
									clust[l1][j] = NAN;
								if (al2 > 0)
									clust[l2][j] = (clust[l2][j]*al2 + iv[r][j])/alt;
								else
									clust[l2][j] = iv[r][j];
								clust[l1][n+j] = alw;
								clust[l2][n+j] = alt;
							}
						}
						obs[r][1] = l2;
						obs[r][2] = l1;
					}
				}
			}
			/*If m steps without transfer, done*/
			if (icoun == m) return;
		}
	}
}

/*K++ Initialization (CURRENTLY USED)*/
void kplusplus(double **clust, double **iv, int k, int n, int m)
{

	double myrand, *dists, wtotal, tracker, newdist,tosq;
	int i, j, ii, jj, divnum;

	/*Choose first obs at random and set cluster to its means*/
	myrand = runif(0,1);	/* RUNIF */
	for (i=0; i<m; i++)
		if ((double) (i+1) / m > myrand)
			break;
	for (ii=0; ii<n; ii++)
		clust[0][ii] = iv[i][ii];

	/*Calc Distances between chosen obs and all other obs*/
	MAKE_VECTOR(dists, m);
	for (j=0; j<m; j++){
		dists[j] = 0;
		divnum = 0;
		for (ii=0; ii<n; ii++){
			if (isfinite(iv[j][ii]) == 1 && isfinite(clust[0][ii])){
				divnum = divnum + 1;
				tosq = clust[0][ii] - iv[j][ii];
				dists[j] = dists[j] + tosq*tosq;
			}
		}
		if (divnum > 1)
			dists[j] = dists[j]/((double) divnum);
	}

	/*Do for remaining clusters...*/
	for (jj=1; jj<k; jj++){
		/*Sum distances to use as weights*/
		wtotal = 0;
		for (j=0; j<m; j++)
			wtotal = wtotal + dists[j];
		/*Choose obs proportional to dists and set cluster jj to its means*/
		myrand = runif(0,1);	/* RNG */
		tracker = 0;
		for (i=0; i<m; i++){
			tracker = tracker + dists[i]/wtotal;
			if (tracker >= myrand)
				break;
		}
		for (ii=0; ii<n; ii++)
			clust[jj][ii] = iv[i][ii];

		/*Calc distances between chosen obs and all other, put in dists if less than current entry*/
		if (jj+1 < k){
			for (j=0; j<m; j++){
				newdist = 0;
				divnum = 0;
				for (ii=0; ii<n; ii++){
					if (isfinite(iv[j][ii]) && isfinite(clust[jj][ii])){
						divnum = divnum + 1;
						tosq = clust[jj][ii] - iv[j][ii];
						newdist = newdist + tosq*tosq;
					}
				}
				if (divnum > 1)
					newdist = newdist/((double) divnum);
				if (newdist < dists[j])
					dists[j] = newdist;
			}
		}
	}
	FREE_VECTOR(dists);
}

int kmeans2(double **dat, double **clust, int **obs, double *tss, int k, int m, int n, int numiter, int iter)
{
/*
  K-Means Algorithm

  Variables:
  **clust: matrix of cluster centers, etc
  **dat: matrix of observed data
  **obs: matrix of observations
  k: number of clusters
  n: number of dimensions
  m: total number of measurements
  *tss: pointer to holder for total sum of squares

  Working Variables:
  **iv: matrix of within sum of squares and sums
  means: holds overall mean in each dimension for initialization
  **dists: distances of each point to overall mean with vector of obs index
  stop: if stop=0, no clusters live after OPTRA
  ERR = 0 Success
  1 Initialization Error
  2 Did Not Converge
  3 Input Error
*/

	double **iv, mydist, mydist2, newdist, newtss, **newclust, *d;
	int i, j, initruns, ERR=0, **newobs, *ncp, *itran, *live, indx, *ind/*,code*/, jj;

	if (k <= 1 || m <= k){
		ERR = 3;
#ifdef MATHLIB_STANDALONE
		printf("Requested number of clusters is less than or equal to 1, or greater than or equal to the number of observations\n");
#else
		Rprintf("Requested number of clusters is less than or equal to 1, or greater than or equal to the number of observations\n");
#endif
		return ERR;
	}
	/*Allocate obs matrix, 0 initialize all count columns, then populate count columns*/
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			obs[i][3+j] = 0;
	for (i=0; i<m; i++) {
		obs[i][0] = 1;
		for (j=0; j<n; j++)
			if (isfinite(dat[i][j]))
				obs[i][3+j]++;
	}
    
	/*Allocate and create inner variance matrix iv and a bunch of other vectors*/
	MAKE_MATRIX(iv, m, n*2);
	MAKE_MATRIX(newclust, k, 2*n+3);
	MAKE_MATRIX(newobs, m, 2);
	MAKE_VECTOR(itran, k);
	MAKE_VECTOR(live, k);
	MAKE_VECTOR(ncp, k);
	MAKE_VECTOR(d, m);
	MAKE_VECTOR(ind, k);
	makeIV(n, m, dat, iv);

	/*Do k-means 1000 times and keep the best run*/
	for (initruns=0; initruns<numiter; initruns++){
		kplusplus(newclust, iv, k, n,m);
		/*Initial cluster assignment*/

		/*For each observation, find and store closest two cluster means*/
		for (i=0; i<m; i++){
			mydist = 0.;
			mydist2 = 0.;
			for (j=0; j<k; j++){
				newdist = 0.;
				switch (j) {
				case 0:
					for (jj=0; jj<n; jj++)
						if (obs[i][jj+3])
							mydist = mydist + (iv[i][jj]-newclust[j][jj])*(iv[i][jj]-newclust[j][jj]);
					obs[i][1] = j;
					break;
				case 1:
					for (jj=0; jj<n; jj++)
						if (obs[i][jj+3] > 0)
							newdist = newdist + (iv[i][jj]-newclust[j][jj])*(iv[i][jj]-newclust[j][jj]);
					if (newdist < mydist){
						obs[i][2] = obs[i][1];
						obs[i][1] = j;
						mydist2 = mydist;
						mydist = newdist;
					}else{
						obs[i][2] = j;
						mydist2 = newdist;
					}
					break;
				default:
					newdist = 0.;
					for (jj=0; jj<n; jj++)
						if (obs[i][jj+3] > 0)
							newdist = newdist + (iv[i][jj]-newclust[j][jj])*(iv[i][jj]-newclust[j][jj]);
					if (newdist < mydist){
						obs[i][2] = obs[i][1];
						obs[i][1] = j;
						mydist2 = mydist;
						mydist = newdist;
					} else if (newdist < mydist2){
						obs[i][2] = j;
						mydist2 = newdist;
					}
				}
			}
		}
		/*Update cluster means, n_j's, and within-sum-of-squares*/
		clustmeans(newclust, iv, obs, k, n, m);
		for (i=0; i<k; i++) {
			if (newclust[i][2*n] == 0){
				ERR = 1;
#ifdef MATHLIB_STANDALONE
				//				printf("A cluster is empty after initial assignment. You may want to try a new initialization\n");
#else
				// Rprintf("A cluster is empty after initial assignment. You may want to try a new initialization\n");
#endif
			}
			/*Mark all clusts as needing to be updated*/
			itran[i] = 1;
			ncp[i] = -1;
		}
		indx = 0;
		/*Run initial OPTRA and then do QTRANS,OPTRA until convergence*/
		for (i=0; i<iter; i++){
			otrans(newclust, dat, obs, iv, k, n, m, ncp, itran, live, &indx, d);

			/*If no transfers, done*/
			if (indx == m){
				i = iter;
			}else{
				qtrans(newclust, dat, obs, iv, k, n, m, ncp, itran, &indx, d);
				if (k == 2){
					ERR = 0;
					i = iter;
					break;
				}else{
					/*Reset ncp*/
					for (j=0; j<k; j++)
						ncp[j] = 0;
				}
			}
		}
		/*Calculate sum of squares and keep results if better*/
		sumsquares(newclust, dat, obs, k, n, m, &newtss, iv);
		if (newtss < *tss || initruns == 0){
			for (i=0; i<m; i++)
				for (j=0; j<2; j++)
					newobs[i][j] = obs[i][1+j];
			for (i=0; i<k; i++)
				for (j=0; j<(2*n+3); j++)
					clust[i][j] = newclust[i][j];
			(*tss) = newtss;
		}

	}

	for (i=0; i<m; i++)
		for (j=0; j<2; j++)
			obs[i][1+j] = newobs[i][j];
	if (ERR == 2)
#ifdef MATHLIB_STANDALONE
		printf("WARNING: MAXIMUM NUMBER OF ITERATIONS REACHED!!!\n");
#else
		Rprintf("WARNING: MAXIMUM NUMBER OF ITERATIONS REACHED!!!\n");
#endif

	FREE_MATRIX(iv);
	FREE_MATRIX(newobs);
	FREE_MATRIX(newclust);
	FREE_VECTOR(itran);
	FREE_VECTOR(live);
	FREE_VECTOR(ncp);
	FREE_VECTOR(d);
	FREE_VECTOR(ind);
	return ERR;
}

/**
 * Wrapper to kmeans2, which runs km-kmeans repeatedly with km-kmeans++
 * initialization and returns the best solution found.
 *
 * @param dat		data: m x n matrix with non-finite values for missing data
 * @param k		number of clusters
 * @param m		number of observations
 * @param n		number of coordinates
 * @param obs		m x (n+3) matrix, need not be initialized
 *			closest mean in column 1, second closest mean in column 2,
 *			indicator of non-missing data in columns 3 to 3 + n
 * @param clust		cluster means
 * @param numiter	number of times to initialize (ignore name)
 */
double repkmmeanspp(double **dat, int k, int m, int n, int **obs, double **clust, int numiter, int iter)
{
	double tss = INFINITY;
	int  i, j, ERR, numindim;

	/*Call to kmeans function only necessary if K>1*/
	if (k == 1){			
		for (j=0; j<n; j++){
			clust[0][j] = 0.0;
			numindim = 0;
			for (i=0; i<m; i++){
				if (isfinite(dat[i][j])){
					clust[0][j] = clust[0][j] + dat[i][j];
					numindim = numindim + 1;
				}
			}
			clust[0][j] = clust[0][j] / numindim;	
		}
		for(i=0; i<m; i++)
			obs[i][1] = 0;
		
		tss = 0.0;
		for (i=0; i<m; i++)
			for (j=0; j<n; j++)
				if (isfinite(dat[i][j]))
					tss = tss + (dat[i][j] - clust[0][j]) * (dat[i][j] - clust[0][j]);
	}else{
		ERR = kmeans2(dat, clust, obs, &tss, k, m, n, numiter, iter);
		if (ERR != 0)
#ifdef MATHLIB_STANDALONE
//			printf("Error code: %d\n", ERR);
#else
//			Rprintf("Error code: %d\n", ERR);
#endif
	}
	return tss;
}
