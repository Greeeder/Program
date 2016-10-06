// ---------------------------------------------------------------------------

#ifndef Fit_MRQH
#define Fit_MRQH

#include "Globales.h"
#include "Correlations.h"

/*!
 * \brief	Object for nonlinear least-squares fitting by the Levenberg-Marquardt method, also including 
 * 			the ability to hold specified parameters at fixed, specified values. Call constructor to bind data
 * 			vectors and fitting functions and to input an initial parameter guess. Then call any combination
 * 			of hold, free, and fit as often as desired. fit sets the output quantities a, covar, alpha, and chisq..
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	19/03/2016
 */

struct Fitmrq {
	static const int NDONE = 4;		   //!< The minimum iteration to success
	static const int ITMAX = 2000;	   //!< The number of maxumum iterations

	int ndat;						   //!< The number of data
	int ma;							   //!< The number of coefficients
	int mfit;						   //!< The number of free coefficients

	dVector &x;						   //!< The x array
	dVector &y;						   //!< The y array							
	dVector &sig;					   //!< The standard deviation for each value
	double tol;						   //!< The tolerance

	stCorrelation *funcs;			   //!< The function to be fitted

	bVector ia;						   //!< true if the coefficient is free
	bVector iamax;					   //!< true if the coefficient has a maximum
	bVector iamin;					   //!< true if the coefficient has a minimum

	dVector a;						   //!< The vector of fitted coefficients
	dVector amax;					   //!< The maximum of the coefficients
	dVector amin;					   //!< The minimum of the coefficients

	dMatrix covar;					   //!< The covariance matrix
	dMatrix alpha;					   //!< The curvature matrix
	double chisq;					   //!< The value of \f$\tilde{\chi}^2\f$ for the fit

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	Fitmrq();

	/*!
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	xx   	The x-array input data.
	 * \param [in,out]	yy   	The y-array output data.
	 * \param [in,out]	ssig 	The standard deviation of each point.
	 * \param [in,out]	aa   	The coefficients array.
	 * \param [in,out]	funks	The function to be fitted.
	 * \param	TOL			 	(optional) the tolerance.
	 */

	Fitmrq(dVector &xx, dVector &yy, dVector &ssig, dVector &aa, stCorrelation *funks, const double TOL = 1.e-3) :
		ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig), tol(TOL), ia(ma), iamax(ma), iamin(ma), a(aa), amax(ma),
		amin(ma) {

		funcs = funks;
		DimdMatrix(alpha, ma, ma);
		DimdMatrix(covar, ma, ma);
		for(int i = 0; i < ma; i++)
			ia[i] = true;
	}

	/*!
	 * \brief	Asign data (if the struct is already constructed).
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	xx   	The x-array input data.
	 * \param [in,out]	yy   	The y-array output data.
	 * \param [in,out]	ssig 	The standard deviation of each point.
	 * \param [in,out]	aa   	The coefficients array.
	 * \param [in,out]	funks	The function to be fitted.
	 * \param	TOL			 	(optional) the tolerance.
	 */

	void AsignData(dVector &xx, dVector &yy, dVector &ssig, dVector &aa, stCorrelation *funks, const double TOL = 1.e-3) {
		ndat = xx.size();
		ma = aa.size();
		x = xx;
		y = yy;
		sig = ssig;
		tol = TOL;
		ia.resize(ma);
		iamax.resize(ma);
		iamin.resize(ma);
		a = aa;
		amax.resize(ma);
		amin.resize(ma);

		funcs = funks;
		DimdMatrix(alpha, ma, ma);
		DimdMatrix(covar, ma, ma);
		for(int i = 0; i < ma; i++)
			ia[i] = true;
	}

	/*!
	 * \brief	Function for holding a parameter, identified by a value i in the range 0 ... ma-1,
	 * 			fixed at the value val.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	i  	Index of the coefficient to be held.
	 * \param	val	The value of the coefficient.
	 */

	void hold(const int i, const double val) {
		ia[i] = false;
		a[i] = val;
	}

	/*!
	 * \brief	Function for freeing a parameter that was previously held fixed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	i	Index of the coefficient to be freed.
	 */

	void free(const int i) {
		ia[i] = true;
	}

	/*!
	 * \brief	Fix the maximum value of a coefficient.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	i  	Index of the coefficient.
	 * \param	val	The maximum value.
	 */

	void max_limit(const int i, const double val) {
		iamax[i] = true;
		amax[i] = val;
	}

	/*!
	 * \brief	Fix the minimum value of a coefficient.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	i  	Index of the coefficient.
	 * \param	val	The minimum value.
	 */

	void min_limit(const int i, const double val) {
		iamin[i] = true;
		amin[i] = val;
	}

	/*!
	 * \brief	Check if any coefficient has reach the limits.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	a1	The coefficients.
	 *
	 * \return	true if it succeeds, false if it fails.
	 */

	bool check_limits(dVector& a1) {
		bool retval = false;
		for(int i = 0; i < ma; i++) {
			if(iamax[i] == true && a1[i] > amax[i]) {
				a1[i] = amax[i];
				retval = true;
			} else if(iamin[i] == true && a1[i] < amin[i]) {
				a1[i] = amin[i];
				retval = true;
			}
		}
		return retval;
	}

	/*!
	 * \brief	Used by fit() to evaluate the linearized fitting matrix alpha, vector beta and to calculate \f$\tilde{\chi}^2\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	a	 	The coefficient array.
	 * \param [in,out]	alpha	The curvature matrix.
	 * \param [in,out]	beta 	The beta.
	 */

	void mrqcof(dVector &a, dMatrix &alpha, dVector &beta) {
		int i, j, k, l, m;
		double ymod, wt, sig2i, dy;
		dVector dyda(ma);
		for(j = 0; j < mfit; j++) {
			// Initialize(symmetric)alpha, beta.
			for(k = 0; k <= j; k++)
				alpha[j][k] = 0.0;
			beta[j] = 0.;
		}
		chisq = 0.;
		for(i = 0; i < ndat; i++) {
			// Summation loop over all data.
			funcs->Fun(x[i], a, ymod, dyda);
			sig2i = 1.0 / (sig[i] * sig[i]);
			dy = y[i] - ymod;
			for(j = 0, l = 0; l < ma; l++) {
				if(ia[l]) {
					wt = dyda[l] * sig2i;
					for(k = 0, m = 0; m < l + 1; m++)
						if(ia[m])
							alpha[j][k++] += wt * dyda[m];
					beta[j++] += dy * wt;
				}
			}
			chisq += dy * dy * sig2i;
			// And find  2.
		}
		for(j = 1; j < mfit; j++)
			// Fill in the symmetric side.
			for(k = 0; k < j; k++)
				alpha[k][j] = alpha[j][k];
	}

	/*!
	 * \brief	Expand in storage the covariance matrix covar, so as to take into account parameters that
	 * 			are being held fixed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	covar	The covariance.
	 */

	void covsrt(dMatrix &covar) {
		int i, j, k;
		for(i = mfit; i < ma; i++)
			for(j = 0; j < i + 1; j++)
				covar[i][j] = covar[j][i] = 0.0;
		k = mfit - 1;
		for(j = ma - 1; j >= 0; j--) {
			if(ia[j]) {
				for(i = 0; i < ma; i++)
					Swap(covar[i][k], covar[i][j]);
				for(i = 0; i < ma; i++)
					Swap(covar[k][i], covar[j][i]);
				k--;
			}
		}
	}

	/*!
	 * \brief	Fits the function to the provide data.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void fit() {
		int j, k, l, iter, done = 0;
		double alamda = 0.001, ochisq;
		dVector atry(ma), beta(ma), da(ma);
		mfit = 0;
		for(j = 0; j < ma; j++) {
			if(ia[j])
				mfit++;
		}
		dMatrix oneda, temp;
		DimdMatrix(oneda, mfit, 1);
		DimdMatrix(temp, mfit, mfit);
		mrqcof(a, alpha, beta);
		for(j = 0; j < ma; j++)
			atry[j] = a[j];
		ochisq = chisq;
		for(iter = 0; iter < ITMAX; iter++) {
			if(done == NDONE)
				alamda = 0.;
			for(j = 0; j < mfit; j++) {
				for(k = 0; k < mfit; k++)
					covar[j][k] = alpha[j][k];
				covar[j][j] = alpha[j][j] * (1.0 + alamda);
				for(k = 0; k < mfit; k++) {
					temp[j][k] = covar[j][k];
				}
				oneda[j][0] = beta[j];
			}
			if(!gaussj(temp, oneda))
				return;
			// Matrix solution.
			for(j = 0; j < mfit; j++) {
				for(k = 0; k < mfit; k++)
					covar[j][k] = temp[j][k];
				da[j] = oneda[j][0];
			}
			if(done == NDONE) {
				// Converged.Clean up and return.
				covsrt(covar);
				covsrt(alpha);
				return;
			}
			for(j = 0, l = 0; l < ma; l++) {
				// Did the trial succeed ?
				if(ia[l])
					atry[l] = a[l] + da[j++];
			}

			check_limits(atry);
			mrqcof(atry, covar, da);

			if(fabs(chisq - ochisq) < Max(tol, tol * chisq))
				done++;
			if(chisq < ochisq) {
				// Success, accept the new solution.
				alamda *= 0.1;
				ochisq = chisq;
				for(j = 0; j < mfit; j++) {
					for(k = 0; k < mfit; k++)
						alpha[j][k] = covar[j][k];
					beta[j] = da[j];
				}
				for(l = 0; l < ma; l++)
					a[l] = atry[l];
			} else {
				// Failure, increase alamda.
				alamda *= 10.0;
				chisq = ochisq;
			}
		}
		//throw("Fitmrq too many iterations");
	}

	/*!
	 * \brief	Linear equation solution by Gauss-Jordan elimination. The input matrix
	 * 			is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input containing the m right-hand side vectors.
	 * 			On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of
	 * 			solution vectors.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	a	The a matrix.
	 * \param [in,out]	b	The b matrix.
	 *
	 * \return	true if it succeeds, false if it fails.
	 */

	bool gaussj(dMatrix &a, dMatrix &b) {
		int i, icol, irow, j, k, l, ll, n = a.size(), m = b[0].size();
		double big, dum, pivinv;
		iVector indxc(n), indxr(n), ipiv(n); // These integer arrays are used for bookkeeping on
		for(j = 0; j < n; j++)
			ipiv[j] = 0; // the pivoting.
		for(i = 0; i < n; i++) {     // This is the main loop over the columns to be
			big = 0.0; // reduced.
			for(j = 0; j < n; j++)     // This is the outer loop of the search for a pivot
				if(ipiv[j] != 1)     // element.
					for(k = 0; k < n; k++) {
						if(ipiv[k] == 0) {
							if(fabs(a[j][k]) >= big) {
								big = fabs(a[j][k]);
								irow = j;
								icol = k;
							}
						}
					}
			++(ipiv[icol]);

			if(irow != icol) {
				for(l = 0; l < n; l++)
					Swap(a[irow][l], a[icol][l]);
				for(l = 0; l < m; l++)
					Swap(b[irow][l], b[icol][l]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if(a[icol][icol] == 0.0) {
				//throw("gaussj: Singular Matrix");
				return false;
			}
			pivinv = 1.0 / a[icol][icol];
			a[icol][icol] = 1.0;
			for(l = 0; l < n; l++)
				a[icol][l] *= pivinv;
			for(l = 0; l < m; l++)
				b[icol][l] *= pivinv;
			for(ll = 0; ll < n; ll++)     // Next, we reduce the rows...
				if(ll != icol) {     // ...except for the pivot one, of course.
					dum = a[ll][icol];
					a[ll][icol] = 0.0;
					for(l = 0; l < n; l++)
						a[ll][l] -= a[icol][l] * dum;
					for(l = 0; l < m; l++)
						b[ll][l] -= b[icol][l] * dum;
				}
		}
		for(l = n - 1; l >= 0; l--) {
			if(indxr[l] != indxc[l])
				for(k = 0; k < n; k++)
					Swap(a[k][indxr[l]], a[k][indxc[l]]);
		}
		return true;
	}

};

/*!
* \brief	Object for multiple nonlinear least-squares fitting by the Levenberg-Marquardt method, also including
* 			the ability to hold specified parameters at fixed, specified values. Call constructor to bind data
* 			vectors and fitting functions and to input an initial parameter guess. Then call any combination
* 			of hold, free, and fit as often as desired. fit sets the output quantities a, covar, alpha, and chisq..
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	19/03/2016
 */

struct FitmrqMult {
	static const int NDONE = 4;		   //!< The minimum iteration to success
	static const int ITMAX = 10000;	   //!< The number of maxumum iterations

	int ndat;						   //!< The number of data
	int ma;							   //!< The number of coefficients
	int mfit;						   //!< The number of free coefficients

	dVector &y;						   //!< The y array							
	dVector &sig;					   //!< The standard deviation for each value
	dMatrix &x;						   //!< The dMatrix&amp; to process
	double tol;						   //!< The tolerance

	stCorrelation *funcs;			   //!< The function to be fitted

	bVector ia;						   //!< true if the coefficient is free
	bVector iamax;					   //!< true if the coefficient has a maximum
	bVector iamin;					   //!< true if the coefficient has a minimum

	dVector a;						   //!< The vector of fitted coefficients
	dVector amax;					   //!< The maximum of the coefficients
	dVector amin;					   //!< The minimum of the coefficients

	dMatrix covar;					   //!< The covariance matrix
	dMatrix alpha;					   //!< The curvature matrix
	double chisq;					   //!< The value of \f$\tilde{\chi}^2\f$ for the fit

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	FitmrqMult();

	/*!
	* \brief	Constructor.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	xx   	The x-matrix input data.
	* \param [in,out]	yy   	The y-array output data.
	* \param [in,out]	ssig 	The standard deviation of each point.
	* \param [in,out]	aa   	The coefficients array.
	* \param [in,out]	funks	The function to be fitted.
	* \param	TOL			 	(optional) the tolerance.
	*/

	FitmrqMult(dMatrix &xx, dVector &yy, dVector &ssig, dVector &aa, stCorrelation *funks, const double TOL = 1.e-3) :
		ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig), tol(TOL), ia(ma), iamax(ma), iamin(ma), a(aa), amax(ma),
		amin(ma) {

		funcs = funks;
		DimdMatrix(alpha, ma, ma);
		DimdMatrix(covar, ma, ma);
		for(int i = 0; i < ma; i++)
			ia[i] = true;
	}

	/*!
	* \brief	Asign data (if the struct is already constructed).
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	xx   	The x-array input data.
	* \param [in,out]	yy   	The y-array output data.
	* \param [in,out]	ssig 	The standard deviation of each point.
	* \param [in,out]	aa   	The coefficients array.
	* \param [in,out]	funks	The function to be fitted.
	* \param	TOL			 	(optional) the tolerance.
	*/

	void AsignData(dMatrix &xx, dVector &yy, dVector &ssig, dVector &aa, stCorrelation *funks, const double TOL = 1.e-3) {
		ndat = xx.size();
		ma = aa.size();
		x = xx;
		y = yy;
		sig = ssig;
		tol = TOL;
		ia.resize(ma);
		iamax.resize(ma);
		iamin.resize(ma);
		a = aa;
		amax.resize(ma);
		amin.resize(ma);

		funcs = funks;
		DimdMatrix(alpha, ma, ma);
		DimdMatrix(covar, ma, ma);
		for(int i = 0; i < ma; i++)
			ia[i] = true;
	}

	/*!
	* \brief	Function for holding a parameter, identified by a value i in the range 0 ... ma-1,
	* 			fixed at the value val.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i  	Index of the coefficient to be held.
	* \param	val	The value of the coefficient.
	*/

	void hold(const int i, const double val) {
		ia[i] = false;
		a[i] = val;
	}

	/*!
	* \brief	Function for freeing a parameter that was previously held fixed.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i	Index of the coefficient to be freed.
	*/

	void free(const int i) {
		ia[i] = true;
	}

	/*!
	* \brief	Fix the maximum value of a coefficient.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i  	Index of the coefficient.
	* \param	val	The maximum value.
	*/

	void max_limit(const int i, const double val) {
		iamax[i] = true;
		amax[i] = val;
	}

	/*!
	* \brief	Fix the minimum value of a coefficient.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i  	Index of the coefficient.
	* \param	val	The minimum value.
	*/

	void min_limit(const int i, const double val) {
		iamin[i] = true;
		amin[i] = val;
	}

	/*!
	* \brief	Check if any coefficient has reach the limits.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	a1	The coefficients.
	*
	* \return	true if it succeeds, false if it fails.
	*/

	bool check_limits(dVector& a1) {
		bool retval = false;
		for(int i = 0; i < ma; i++) {
			if(iamax[i] == true && a1[i] > amax[i]) {
				a1[i] = amax[i];
				retval = true;
			} else if(iamin[i] == true && a1[i] < amin[i]) {
				a1[i] = amin[i];
				retval = true;
			}
		}
		return retval;
	}

	/*!
	* \brief	Used by fit() to evaluate the linearized fitting matrix alpha, vector beta and to calculate \f$\tilde{\chi}^2\f$.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	a	 	The coefficient array.
	* \param [in,out]	alpha	The curvature matrix.
	* \param [in,out]	beta 	The beta.
	*/

	void mrqcof(dVector &a, dMatrix &alpha, dVector &beta) {
		int i, j, k, l, m;
		double ymod, wt, sig2i, dy;
		dVector dyda(ma);
		for(j = 0; j < mfit; j++) {
			// Initialize(symmetric)alpha, beta.
			for(k = 0; k <= j; k++)
				alpha[j][k] = 0.0;
			beta[j] = 0.;
		}
		chisq = 0.;
		for(i = 0; i < ndat; i++) {
			// Summation loop over all data.
			funcs->Fun(x[i], a, ymod, dyda);
			sig2i = 1.0 / (sig[i] * sig[i]);
			dy = y[i] - ymod;
			for(j = 0, l = 0; l < ma; l++) {
				if(ia[l]) {
					wt = dyda[l] * sig2i;
					for(k = 0, m = 0; m < l + 1; m++)
						if(ia[m])
							alpha[j][k++] += wt * dyda[m];
					beta[j++] += dy * wt;
				}
			}
			chisq += dy * dy * sig2i;
			// And find  2.
		}
		for(j = 1; j < mfit; j++)
			// Fill in the symmetric side.
			for(k = 0; k < j; k++)
				alpha[k][j] = alpha[j][k];
	}

	/*!
	* \brief	Expand in storage the covariance matrix covar, so as to take into account parameters that
	* 			are being held fixed.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	covar	The covariance.
	*/

	void covsrt(dMatrix &covar) {
		int i, j, k;
		for(i = mfit; i < ma; i++)
			for(j = 0; j < i + 1; j++)
				covar[i][j] = covar[j][i] = 0.0;
		k = mfit - 1;
		for(j = ma - 1; j >= 0; j--) {
			if(ia[j]) {
				for(i = 0; i < ma; i++)
					Swap(covar[i][k], covar[i][j]);
				for(i = 0; i < ma; i++)
					Swap(covar[k][i], covar[j][i]);
				k--;
			}
		}
	}

	/*!
	* \brief	Fits the function to the provide data.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*/

	void fit() {
		int j, k, l, iter, done = 0;
		double alamda = 0.001, ochisq;
		dVector atry(ma), beta(ma), da(ma);
		mfit = 0;
		for(j = 0; j < ma; j++) {
			if(ia[j])
				mfit++;
		}
		dMatrix oneda, temp;
		DimdMatrix(oneda, mfit, 1);
		DimdMatrix(temp, mfit, mfit);
		mrqcof(a, alpha, beta);
		for(j = 0; j < ma; j++)
			atry[j] = a[j];
		ochisq = chisq;
		for(iter = 0; iter < ITMAX; iter++) {
			if(done == NDONE)
				alamda = 0.;
			for(j = 0; j < mfit; j++) {
				for(k = 0; k < mfit; k++)
					covar[j][k] = alpha[j][k];
				covar[j][j] = alpha[j][j] * (1.0 + alamda);
				for(k = 0; k < mfit; k++) {
					temp[j][k] = covar[j][k];
				}
				oneda[j][0] = beta[j];
			}
			if(!gaussj(temp, oneda))
				return;
			// Matrix solution.
			for(j = 0; j < mfit; j++) {
				for(k = 0; k < mfit; k++)
					covar[j][k] = temp[j][k];
				da[j] = oneda[j][0];
			}
			if(done == NDONE) {
				// Converged.Clean up and return.
				covsrt(covar);
				covsrt(alpha);
				return;
			}
			for(j = 0, l = 0; l < ma; l++) {
				// Did the trial succeed ?
				if(ia[l])
					atry[l] = a[l] + da[j++];
			}

			check_limits(atry);
			mrqcof(atry, covar, da);

			if(fabs(chisq - ochisq) < Max(tol, tol * chisq))
				done++;
			if(chisq < ochisq) {
				// Success, accept the new solution.
				alamda *= 0.1;
				ochisq = chisq;
				for(j = 0; j < mfit; j++) {
					for(k = 0; k < mfit; k++)
						alpha[j][k] = covar[j][k];
					beta[j] = da[j];
				}
				for(l = 0; l < ma; l++)
					a[l] = atry[l];
			} else {
				// Failure, increase alamda.
				alamda *= 10.0;
				chisq = ochisq;
			}
		}
		//throw("Fitmrq too many iterations");
	}

	/*!
	* \brief	Linear equation solution by Gauss-Jordan elimination. The input matrix
	* 			is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input containing the m right-hand side vectors.
	* 			On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of
	* 			solution vectors.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	a	The a matrix.
	* \param [in,out]	b	The b matrix.
	*
	* \return	true if it succeeds, false if it fails.
	*/

	bool gaussj(dMatrix &a, dMatrix &b) {
		int i, icol, irow, j, k, l, ll, n = a.size(), m = b[0].size();
		double big, dum, pivinv;
		iVector indxc(n), indxr(n), ipiv(n); // These integer arrays are used for bookkeeping on
		for(j = 0; j < n; j++)
			ipiv[j] = 0; // the pivoting.
		for(i = 0; i < n; i++) {     // This is the main loop over the columns to be
			big = 0.0; // reduced.
			for(j = 0; j < n; j++)     // This is the outer loop of the search for a pivot
				if(ipiv[j] != 1)     // element.
					for(k = 0; k < n; k++) {
						if(ipiv[k] == 0) {
							if(fabs(a[j][k]) >= big) {
								big = fabs(a[j][k]);
								irow = j;
								icol = k;
							}
						}
					}
			++(ipiv[icol]);

			if(irow != icol) {
				for(l = 0; l < n; l++)
					Swap(a[irow][l], a[icol][l]);
				for(l = 0; l < m; l++)
					Swap(b[irow][l], b[icol][l]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if(a[icol][icol] == 0.0) {
				cerr << "gaussj: Singular Matrix" << endl;
				return false;
				//throw("gaussj: Singular Matrix");
			}
			pivinv = 1.0 / a[icol][icol];
			a[icol][icol] = 1.0;
			for(l = 0; l < n; l++)
				a[icol][l] *= pivinv;
			for(l = 0; l < m; l++)
				b[icol][l] *= pivinv;
			for(ll = 0; ll < n; ll++)     // Next, we reduce the rows...
				if(ll != icol) {     // ...except for the pivot one, of course.
					dum = a[ll][icol];
					a[ll][icol] = 0.0;
					for(l = 0; l < n; l++)
						a[ll][l] -= a[icol][l] * dum;
					for(l = 0; l < m; l++)
						b[ll][l] -= b[icol][l] * dum;
				}
		}
		for(l = n - 1; l >= 0; l--) {
			if(indxr[l] != indxc[l])
				for(k = 0; k < n; k++)
					Swap(a[k][indxr[l]], a[k][indxc[l]]);
		}
		return true;
	}

};

/*!
* \brief	Object for multiple nonlinear least-squares fitting by the Levenberg-Marquardt method, also including
* 			the ability to hold specified parameters at fixed, specified values. Call constructor to bind data
* 			vectors and fitting functions and to input an initial parameter guess. Then call any combination
* 			of hold, free, and fit as often as desired. fit sets the output quantities a, covar, alpha, and chisq..
*
* \author	F.J. Arnau (farnau@mot.upv.es)
* \date	19/03/2016
*/

struct FitmrqMult2 {
	static const int NDONE = 4;		   //!< The minimum iteration to success
	static const int ITMAX = 1000;	   //!< The number of maxumum iterations

	int ndat;						   //!< The number of data
	int ma;							   //!< The number of coefficients
	int mfit;						   //!< The number of free coefficients

	dVector &y;						   //!< The y array							
	dVector &sig;					   //!< The standard deviation for each value
	dMatrix &x;						   //!< The dMatrix&amp; to process
	double tol;						   //!< The tolerance

	/*!
	 * \brief	Funcs to be fitted.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	parameter1	The first parameter.
	 * \param [in,out]	parameter2	The second parameter.
	 * \param [in,out]	parameter3	The third parameter.
	 * \param [in,out]	parameter4	The fourth parameter.
	 */

	void (*funcs)(dVector&, dVector&, double&, dVector&);

	bVector ia;						   //!< true if the coefficient is free
	bVector iamax;					   //!< true if the coefficient has a maximum
	bVector iamin;					   //!< true if the coefficient has a minimum

	dVector a;						   //!< The vector of fitted coefficients
	dVector amax;					   //!< The maximum of the coefficients
	dVector amin;					   //!< The minimum of the coefficients

	dMatrix covar;					   //!< The covariance matrix
	dMatrix alpha;					   //!< The curvature matrix
	double chisq;					   //!< The value of \f$\tilde{\chi}^2\f$ for the fit

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	FitmrqMult2();

	/*!
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	xx   	The x-matrix input data.
	 * \param [in,out]	yy   	The y-array output data.
	 * \param [in,out]	ssig 	The standard deviation of each point.
	 * \param [in,out]	aa   	The coefficients array.
	 * \param [in,out]	funks  	The function to be fitted.
	 * \param	TOL				(optional) the tolerance.
	 */

	FitmrqMult2(dMatrix &xx, dVector &yy, dVector &ssig, dVector &aa, void funks(dVector&, dVector&, double &, dVector&),
				const double TOL = 1.e-3) :
		ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig), tol(TOL), funcs(funks), ia(ma), iamax(ma), iamin(ma), a(aa),
		amax(ma), amin(ma) {

		DimdMatrix(alpha, ma, ma);
		DimdMatrix(covar, ma, ma);

		for(int i = 0; i < ma; i++) {
			ia[i] = true;
			iamax[i] = false;
			iamin[i] = false;
		}
	}

	/*!
	* \brief	Function for holding a parameter, identified by a value i in the range 0 ... ma-1,
	* 			fixed at the value val.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i  	Index of the coefficient to be held.
	* \param	val	The value of the coefficient.
	*/

	void hold(const int i, const double val) {
		ia[i] = false;
		a[i] = val;
	}

	/*!
	* \brief	Function for freeing a parameter that was previously held fixed.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i	Index of the coefficient to be freed.
	*/

	void free(const int i) {
		ia[i] = true;
	}

	/*!
	* \brief	Fix the maximum value of a coefficient.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i  	Index of the coefficient.
	* \param	val	The maximum value.
	*/

	void max_limit(const int i, const double val) {
		iamax[i] = true;
		amax[i] = val;
	}

	/*!
	* \brief	Fix the minimum value of a coefficient.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param	i  	Index of the coefficient.
	* \param	val	The minimum value.
	*/

	void min_limit(const int i, const double val) {
		iamin[i] = true;
		amin[i] = val;
	}

	/*!
	* \brief	Check if any coefficient has reach the limits.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	a1	The coefficients.
	*
	* \return	true if it succeeds, false if it fails.
	*/

	void check_limits(dVector& a1) {
		for(int i = 0; i < ma; i++) {
			if(iamax[i] == true && a1[i] > amax[i])
				a1[i] = amax[i];
			else if(iamin[i] == true && a1[i] < amin[i])
				a1[i] = amin[i];
		}
	}

	/*!
	* \brief	Used by fit() to evaluate the linearized fitting matrix alpha, vector beta and to calculate \f$\tilde{\chi}^2\f$.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	a	 	The coefficient array.
	* \param [in,out]	alpha	The curvature matrix.
	* \param [in,out]	beta 	The beta.
	*/

	void mrqcof(dVector &a, dMatrix &alpha, dVector &beta) {
		int i, j, k, l, m;
		double ymod, wt, sig2i, dy;
		dVector dyda(ma);
		for(j = 0; j < mfit; j++) {
			// Initialize(symmetric)alpha, beta.
			for(k = 0; k <= j; k++)
				alpha[j][k] = 0.0;
			beta[j] = 0.;
		}
		chisq = 0.;
		for(i = 0; i < ndat; i++) {
			// Summation loop over all data.
			funcs(x[i], a, ymod, dyda);
			sig2i = 1.0 / (sig[i] * sig[i]);
			dy = y[i] - ymod;
			for(j = 0, l = 0; l < ma; l++) {
				if(ia[l]) {
					wt = dyda[l] * sig2i;
					for(k = 0, m = 0; m < l + 1; m++)
						if(ia[m])
							alpha[j][k++] += wt * dyda[m];
					beta[j++] += dy * wt;
				}
			}
			chisq += dy * dy * sig2i;
			// And find  2.
		}
		for(j = 1; j < mfit; j++)
			// Fill in the symmetric side.
			for(k = 0; k < j; k++)
				alpha[k][j] = alpha[j][k];
	}

	/*!
	* \brief	Expand in storage the covariance matrix covar, so as to take into account parameters that
	* 			are being held fixed.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	covar	The covariance.
	*/

	void covsrt(dMatrix &covar) {
		int i, j, k;
		for(i = mfit; i < ma; i++)
			for(j = 0; j < i + 1; j++)
				covar[i][j] = covar[j][i] = 0.0;
		k = mfit - 1;
		for(j = ma - 1; j >= 0; j--) {
			if(ia[j]) {
				for(i = 0; i < ma; i++)
					Swap(covar[i][k], covar[i][j]);
				for(i = 0; i < ma; i++)
					Swap(covar[k][i], covar[j][i]);
				k--;
			}
		}
	}

	/*!
	* \brief	Fits the function to the provide data.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*/

	void fit() {
		int j, k, l, iter, done = 0;
		double alamda = 0.001, ochisq;
		dVector atry(ma), beta(ma), da(ma);
		mfit = 0;
		for(j = 0; j < ma; j++)
			if(ia[j])
				mfit++;
		dMatrix oneda, temp;
		DimdMatrix(oneda, mfit, 1);
		DimdMatrix(temp, mfit, mfit);
		mrqcof(a, alpha, beta);
		for(j = 0; j < ma; j++)
			atry[j] = a[j];
		ochisq = chisq;
		for(iter = 0; iter < ITMAX; iter++) {
			if(done == NDONE)
				alamda = 0.;
			for(j = 0; j < mfit; j++) {
				for(k = 0; k < mfit; k++)
					covar[j][k] = alpha[j][k];
				covar[j][j] = alpha[j][j] * (1.0 + alamda);
				for(k = 0; k < mfit; k++)
					temp[j][k] = covar[j][k];
				oneda[j][0] = beta[j];
			}
			gaussj(temp, oneda);
			// Matrix solution.
			for(j = 0; j < mfit; j++) {
				for(k = 0; k < mfit; k++)
					covar[j][k] = temp[j][k];
				da[j] = oneda[j][0];
			}
			if(done == NDONE) {
				// Converged.Clean up and return.
				covsrt(covar);
				covsrt(alpha);
				return;
			}
			for(j = 0, l = 0; l < ma; l++)
				// Did the trial succeed ?
				if(ia[l])
					atry[l] = a[l] + da[j++];
			check_limits(atry);
			mrqcof(atry, covar, da);
			if(fabs(chisq - ochisq) < Max(tol, tol * chisq))
				done++;
			if(chisq < ochisq) {
				// Success, accept the new solution.alamda *= 0.1;
				ochisq = chisq;
				for(j = 0; j < mfit; j++) {
					for(k = 0; k < mfit; k++)
						alpha[j][k] = covar[j][k];
					beta[j] = da[j];
				}
				for(l = 0; l < ma; l++)
					a[l] = atry[l];
			} else {
				// Failure, increase alamda.alamda *= 10.0;
				chisq = ochisq;
			}
		}
		throw("Fitmrq too many iterations");
	}

	/*!
	* \brief	Linear equation solution by Gauss-Jordan elimination. The input matrix
	* 			is a[0..n-1][0..n-1]. b[0..n-1][0..m-1] is input containing the m right-hand side vectors.
	* 			On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of
	* 			solution vectors.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	19/03/2016
	*
	* \param [in,out]	a	The a matrix.
	* \param [in,out]	b	The b matrix.
	*
	* \return	true if it succeeds, false if it fails.
	*/

	void gaussj(dMatrix &a, dMatrix &b) {
		int i, icol, irow, j, k, l, ll, n = a.size(), m = b[0].size();
		double big, dum, pivinv;
		iVector indxc(n), indxr(n), ipiv(n); // These integer arrays are used for bookkeeping on
		for(j = 0; j < n; j++)
			ipiv[j] = 0; // the pivoting.
		for(i = 0; i < n; i++) {     // This is the main loop over the columns to be
			big = 0.0; // reduced.
			for(j = 0; j < n; j++)     // This is the outer loop of the search for a pivot
				if(ipiv[j] != 1)     // element.
					for(k = 0; k < n; k++) {
						if(ipiv[k] == 0) {
							if(fabs(a[j][k]) >= big) {
								big = fabs(a[j][k]);
								irow = j;
								icol = k;
							}
						}
					}
			++(ipiv[icol]);

			if(irow != icol) {
				for(l = 0; l < n; l++)
					Swap(a[irow][l], a[icol][l]);
				for(l = 0; l < m; l++)
					Swap(b[irow][l], b[icol][l]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if(a[icol][icol] == 0.0)
				throw("gaussj: Singular Matrix");
			pivinv = 1.0 / a[icol][icol];
			a[icol][icol] = 1.0;
			for(l = 0; l < n; l++)
				a[icol][l] *= pivinv;
			for(l = 0; l < m; l++)
				b[icol][l] *= pivinv;
			for(ll = 0; ll < n; ll++)     // Next, we reduce the rows...
				if(ll != icol) {     // ...except for the pivot one, of course.
					dum = a[ll][icol];
					a[ll][icol] = 0.0;
					for(l = 0; l < n; l++)
						a[ll][l] -= a[icol][l] * dum;
					for(l = 0; l < m; l++)
						b[ll][l] -= b[icol][l] * dum;
				}
		}
		for(l = n - 1; l >= 0; l--) {
			if(indxr[l] != indxc[l])
				for(k = 0; k < n; k++)
					Swap(a[k][indxr[l]], a[k][indxc[l]]);
		}
	}

};
// ---------------------------------------------------------------------------
#endif

