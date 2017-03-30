
#include <qpDUNES.h>
#include <string.h>
#include "./spring_mass_include/spring_mass_data_3.c"
#include "./spring_mass_include/spring_mass_data_4.c"
#include "./spring_mass_include/spring_mass_data_5.c"
#include "./spring_mass_include/spring_mass_data_6.c"
#include "./spring_mass_include/spring_mass_data_7.c"
#include "./spring_mass_include/spring_mass_data_8.c"
#include "./spring_mass_include/spring_mass_data_9.c"
#include "./spring_mass_include/spring_mass_data_10.c"

#define MAXITER 100

// cd .. && make && cd examples/ && ./spring_mass

// TODO(dimitris): Check that we always do Newton and NOT gradient steps (dual_qp.c, line 93)
// TODO(dimitris): if options.nwtnHssnFacAlg ~= QPDUNES_NH_FAC_BAND_REVERSE + blasfeo, throw error

void assign_data(int NM, double **A, double **B, double **c, double **Q, double **P, double **R, 
	double **xLow, double **xUpp, double **xNLow, double **xNUpp, double **uLow, double **uUpp, 
	double **x0, int *NX, int *NU);

void write_timings_to_txt_files(int NM, int nIter, double *tNwtnSetup, double *tNwtnFactor, 
	double *tNwtnSolve, double *tQP, double *tLineSearch, double *tExtra, double *tIter);

int main( )
{
	int N = 20;
	int NX, NU;
	int NM[] = {3, 4, 5, 6, 7, 8, 9, 10};  // NOTE: in ascending order
	int nexp = sizeof(NM)/sizeof(int);
	int NRUNS = 100;
	int i, j;
	boolean_t isLTI;

	double *A;
	double *B;
	double *c;
	double *Q;
	double *P;
	double *R;
	double *xLow;
	double *xUpp;
	double *xNLow;
	double *xNUpp;
	double *uLow;
	double *uUpp;
	double *x0;
	double *S = 0; 
	unsigned int* nD = 0;  // no affine constraints

	// logging
	int nIter;
	double totTime, newtonTime, minTotTime;
	double tNwtnSetup[MAXITER];
	double tNwtnFactor[MAXITER];
	double tNwtnSolve[MAXITER];
	double tQP[MAXITER];
	double tLineSearch[MAXITER];
	double tExtra[MAXITER];
	double tIter[MAXITER];

	qpData_t qpData;

	// NOTE: options need to be set before setup

	qpOptions_t qpOptions = qpDUNES_setupDefaultOptions();
	qpOptions.maxIter = MAXITER;
	qpOptions.printLevel = 1; // 0 = no output, 1 = only errors and success, 2 = additionally iterations and warnings,  3 = debug information, 4 = active set info?
	qpOptions.printIterationTiming = 0;
	qpOptions.logLevel = QPDUNES_LOG_ITERATIONS; // QPDUNES_LOG_OFF, QPDUNES_LOG_ITERATIONS, QPDUNES_LOG_ALL_DATA
	// ...

	for (int iexp = 0; iexp < nexp; iexp++) {

		assign_data(NM[iexp], &A, &B, &c, &Q, &P, &R, &xLow, &xUpp, &xNLow, &xNUpp, &uLow, &uUpp, &x0, &NX, &NU);

		for (j = 0; j < NRUNS; j++) {
		
			// TODO(dimitris): supress prints of qpDUNES_allocInterval inside here
			qpDUNES_setup(&qpData, N, NX, NU, nD, &qpOptions);  // passing 0 in the last argument sets the default QP options

			qpDUNES_setupSimpleBoundedInterval(&qpData, qpData.intervals[0], Q, R, S, A, B, c, x0, x0, uLow, uUpp);
		//   printf("H_TYPE = %d\n", qpData.intervals[0]->H.sparsityType);
			for (i = 1; i < N; ++i) {
				qpDUNES_setupSimpleBoundedInterval(&qpData, qpData.intervals[i], Q, R, S, A, B, c, xLow, xUpp, uLow, uUpp);
			}
			qpDUNES_setupSimpleBoundedInterval(&qpData, qpData.intervals[N], P, 0, 0, 0, 0, 0, xNLow, xNUpp, 0, 0);

			qpDUNES_setupAllLocalQPs(&qpData, isLTI=QPDUNES_TRUE);  // determine local QP solvers and set up auxiliary data

			qpDUNES_solve(&qpData);
			
			// for (i = 0; i < N; ++i) {
			// 	qpDUNES_printMatrixData(qpData.intervals[i]->z.data, 1, NX+NU, "z[%d]:", i);
			// }
			// qpDUNES_printMatrixData(qpData.intervals[N]->z.data, 1, NX, "z[%d]:", i);
			
			nIter = qpData.log.numIter;
			// check if cpu time is better than the current min one
			totTime = 0;
			for (i = 1; i < nIter; i++) {
				totTime += qpData.log.itLog[i].tIt;
			}
			if (j == 0) minTotTime = totTime;

			if (minTotTime >= totTime) { 
				minTotTime = totTime;
				// if yes, store detailed timings
				for (i = 1; i < nIter; i++) {
					tNwtnSetup[i] = qpData.log.itLog[i].tNwtnSetup;
					tNwtnFactor[i] = qpData.log.itLog[i].tNwtnFactor;
					tNwtnSolve[i] = qpData.log.itLog[i].tNwtnSolve;
					tQP[i] = qpData.log.itLog[i].tQP;
					tLineSearch[i] = qpData.log.itLog[i].tLineSearch;
					tExtra[i] = qpData.log.itLog[i].tExtra;
					tIter[i] = qpData.log.itLog[i].tIt;
				}
			}
			
			qpDUNES_cleanup(&qpData);
			printf("end of run %d/%d (%d iterations) \n", j+1, NRUNS, nIter);
		}

		for (i = 1; i < nIter; i++) {
			qpDUNES_printf("\nTimings Iteration %d:", i);
			qpDUNES_printf("Setup of Newton system:         %7.3f ms (%5.2f%%)",
					1e3 * (tNwtnSetup[i]) / 1, (tNwtnSetup[i]) / (tIter[i]) * 100);
			qpDUNES_printf("Factorization of Newton system: %7.3f ms (%5.2f%%)",
					1e3 * (tNwtnFactor[i]) / 1, (tNwtnFactor[i]) / (tIter[i]) * 100);
			qpDUNES_printf("Backsolve of newton system:     %7.3f ms (%5.2f%%)",
					1e3 * (tNwtnSolve[i]) / 1, (tNwtnSolve[i]) / (tIter[i]) * 100);
			qpDUNES_printf("QP solution:                    %7.3f ms (%5.2f%%)",
					1e3 * (tQP[i]) / 1, (tQP[i]) / (tIter[i]) * 100);
			qpDUNES_printf("Line search:                    %7.3f ms (%5.2f%%)",
					1e3 * (tLineSearch[i]) / 1, (tLineSearch[i]) / (tIter[i]) * 100);
			qpDUNES_printf("Overhead (print+log):           %7.3f ms (%5.2f%%)",
					1e3 * (tExtra[i]) / 1, (tExtra[i]) / (tIter[i]) * 100);
			qpDUNES_printf("                               -----------");
			qpDUNES_printf("Full iteration:                 %7.3f ms\n",
					1e3 * (tIter[i]) / 1);
		}

		write_timings_to_txt_files(NM[iexp], nIter, tNwtnSetup, tNwtnFactor, tNwtnSolve, tQP, tLineSearch, tExtra, tIter);

		totTime = 0;
		newtonTime = 0;
		for (i = 1; i < nIter ; i++) {
			totTime += tIter[i];
			newtonTime += tNwtnFactor[i] + tNwtnSolve[i];
		}
		printf("Total time for %d masses:  %7.3f ms\n", NM[iexp], 1000*totTime);
		printf("Newton time for %d masses: %7.3f ms\n", NM[iexp], 1000*newtonTime);

	}

	printf("spring-mass example done.\n");

	return 0;
}


void assign_data(int NM, double **A, double **B, double **c, double **Q, double **P, double **R, 
	double **xLow, double **xUpp, double **xNLow, double **xNUpp, double **uLow, double **uUpp, 
	double **x0, int *NX, int *NU) {

	*NX = 2*NM;
	
	switch (NM) {

		case 3:
			*A = A_3;
			*B = B_3;
			*c = c_3;
			*Q = Q_3;
			*P = P_3;
			*R = R_3;
			*xLow = xLow_3;
			*xUpp = xUpp_3;
			*xNLow = xNLow_3;
			*xNUpp = xNUpp_3;
			*uLow = uLow_3;
			*uUpp = uUpp_3;
			*x0 = x0_3;
			*NU = sizeof(B_3)/sizeof(x0_3);
			break;
		case 4:
			*A = A_4;
			*B = B_4;
			*c = c_4;
			*Q = Q_4;
			*P = P_4;
			*R = R_4;
			*xLow = xLow_4;
			*xUpp = xUpp_4;
			*xNLow = xNLow_4;
			*xNUpp = xNUpp_4;
			*uLow = uLow_4;
			*uUpp = uUpp_4;
			*x0 = x0_4;
			*NU = sizeof(B_4)/sizeof(x0_4);
			break;
		case 5:
			*A = A_5;
			*B = B_5;
			*c = c_5;
			*Q = Q_5;
			*P = P_5;
			*R = R_5;
			*xLow = xLow_5;
			*xUpp = xUpp_5;
			*xNLow = xNLow_5;
			*xNUpp = xNUpp_5;
			*uLow = uLow_5;
			*uUpp = uUpp_5;
			*x0 = x0_5;
			*NU = sizeof(B_5)/sizeof(x0_5);
			break;
		case 6:
			*A = A_6;
			*B = B_6;
			*c = c_6;
			*Q = Q_6;
			*P = P_6;
			*R = R_6;
			*xLow = xLow_6;
			*xUpp = xUpp_6;
			*xNLow = xNLow_6;
			*xNUpp = xNUpp_6;
			*uLow = uLow_6;
			*uUpp = uUpp_6;
			*x0 = x0_6;
			*NU = sizeof(B_6)/sizeof(x0_6);
			break;
		case 7:
			*A = A_7;
			*B = B_7;
			*c = c_7;
			*Q = Q_7;
			*P = P_7;
			*R = R_7;
			*xLow = xLow_7;
			*xUpp = xUpp_7;
			*xNLow = xNLow_7;
			*xNUpp = xNUpp_7;
			*uLow = uLow_7;
			*uUpp = uUpp_7;
			*x0 = x0_7;
			*NU = sizeof(B_7)/sizeof(x0_7);
			break;
		case 8:
			*A = A_8;
			*B = B_8;
			*c = c_8;
			*Q = Q_8;
			*P = P_8;
			*R = R_8;
			*xLow = xLow_8;
			*xUpp = xUpp_8;
			*xNLow = xNLow_8;
			*xNUpp = xNUpp_8;
			*uLow = uLow_8;
			*uUpp = uUpp_8;
			*x0 = x0_8;
			*NU = sizeof(B_8)/sizeof(x0_8);
			break;
		case 9:
			*A = A_9;
			*B = B_9;
			*c = c_9;
			*Q = Q_9;
			*P = P_9;
			*R = R_9;
			*xLow = xLow_9;
			*xUpp = xUpp_9;
			*xNLow = xNLow_9;
			*xNUpp = xNUpp_9;
			*uLow = uLow_9;
			*uUpp = uUpp_9;
			*x0 = x0_9;
			*NU = sizeof(B_9)/sizeof(x0_9);
			break;
		case 10:
			*A = A_10;
			*B = B_10;
			*c = c_10;
			*Q = Q_10;
			*P = P_10;
			*R = R_10;
			*xLow = xLow_10;
			*xUpp = xUpp_10;
			*xNLow = xNLow_10;
			*xNUpp = xNUpp_10;
			*uLow = uLow_10;
			*uUpp = uUpp_10;
			*x0 = x0_10;
			*NU = sizeof(B_10)/sizeof(x0_10);
			break;
	}
}


void write_timings_to_txt_files(int NM, int nIter, double *tNwtnSetup, double *tNwtnFactor, 
	double *tNwtnSolve, double *tQP, double *tLineSearch, double *tExtra, double *tIter) {
		char fname[256];
		char fpath[] = "spring_mass_log/";

		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tNwtnSetup_", NM, ".txt");
		write_double_vector_to_txt(&tNwtnSetup[1], nIter-1, fname);
		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tNwtnFactor_", NM, ".txt");
		write_double_vector_to_txt(&tNwtnFactor[1], nIter-1, fname);
		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tNwtnSolve_", NM, ".txt");
		write_double_vector_to_txt(&tNwtnSolve[1], nIter-1, fname);
		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tQP_", NM, ".txt");
		write_double_vector_to_txt(&tQP[1], nIter-1, fname);
		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tLineSearch_", NM, ".txt");
		write_double_vector_to_txt(&tLineSearch[1], nIter-1, fname);
		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tExtra_", NM, ".txt");
		write_double_vector_to_txt(&tExtra[1], nIter-1, fname);
		snprintf(fname, sizeof(fname), "%s%s%d%s", fpath, "tIter_", NM, ".txt");
		write_double_vector_to_txt(&tIter[1], nIter-1, fname);
}
