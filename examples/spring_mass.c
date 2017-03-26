
#include <qpDUNES.h>
#include "./spring_mass_data.c"

#define INFTY 1.0e12

// cd .. && make && cd examples/ && ./spring_mass

// TODO(dimitris): Check that we always do Newton and NOT gradient steps (dual_qp.c, line 93)
// TODO(dimitris): if options.nwtnHssnFacAlg ~= QPDUNES_NH_FAC_BAND_REVERSE + blasfeo, throw error

int main( )
{
	printf("DIMENSIONS READ FROM FILE: N = %d, NX = %d, NU = %d\n", N, NX, NU);

	int NRUNS = 50;
	int i, j;
	boolean_t isLTI;
	unsigned int* nD = 0;  // no affine constraints
	double *S = 0; 
	
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

	for (j = 0; j < NRUNS; j++) {
	
		// TODO(dimitris): supress prints of qpDUNES_allocInterval inside here
		qpDUNES_setup(&qpData, N, NX, NU, nD, &qpOptions);  // passing 0 in the last argument sets the default QP options

		qpDUNES_setupSimpleBoundedInterval(&qpData, qpData.intervals[0], Q, R, S, A, B, c, x0, x0, uLow, uUpp);
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

		if (minTotTime > totTime) { 
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
				totTime += tIter[i];
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

	totTime = 0;
	newtonTime = 0;
	for (i = 1; i < nIter ; i++) {
		totTime += tIter[i];
		newtonTime += tNwtnFactor[i] + tNwtnSolve[i];
	}

	printf("Total time:  %7.3f ms\n", 1000*totTime);
	printf("Newton time: %7.3f ms\n", 1000*newtonTime);

	printf("spring-mass example done.\n");

	return 0;
}
