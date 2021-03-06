/*
 *    This file was auto-generated by ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2010 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC) under
 *    supervision of Moritz Diehl. All rights reserved.
 *
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

#include <string.h>

#include <stdio.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/time.h>

#include <stdlib.h>

using namespace std;

#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#define NX          ACADO_NX      /* number of differential states  */
#define NU          ACADO_NU      /* number of control inputs       */
#define N          	ACADO_N      /* number of control intervals    */
#define NY			ACADO_NY
#define NYN			ACADO_NYN	
#define NUM_STEPS   200     /* number of real time iterations
// #define NUM_STEPS   2     /* number of real time iterations */
#define VERBOSE     1      /* show iterations: 1, silent: 0  */

ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

extern "C"
{

#if ACADO_QP_SOLVER == ACADO_FORCES

#include "forces.h"
extern forces_info acadoForces_info;
extern forces_params acadoForces_params;
extern forces_output acadoForces_output;

extern void evaluateConstraints();

#elif ACADO_QP_SOLVER == ACADO_QPOASES

extern void condensePrep();

#elif ACADO_QP_SOLVER == ACADO_QPDUNES

#else
#error("Unknown QP solver")
#endif

extern void modelSimulation();
extern void evaluateObjective();
extern int feedbackStep( );

}

bool readDataFromFile( const char* fileName, vector< vector< double > >& data )
{
	ifstream file( fileName );
	string line;
	
	if ( file.is_open() )
	{
		while( getline(file, line) )
		{
			istringstream linestream( line );
			vector< double > linedata;
			double number;
			
			while( linestream >> number )
			{
				linedata.push_back( number );
			}
			
			data.push_back( linedata );
		}
		
		file.close();
	}
	else
		return false;
	
	return true;
}

void printMatrix(real_t* data, const string& name, unsigned nRows, unsigned nCols)
{
	cout << name << ": " << endl;
	for (unsigned i = 0; i < nRows; ++i)
	{
		for (unsigned j = 0; j < nCols; ++j)
			cout << data[i * nCols + j] << "   ";

		cout << endl;
	}
	cout << endl;
}

void isSymmetric(real_t* data, const string& name, unsigned nRows, unsigned nCols)
{
	for (unsigned i = 0; i < nRows; ++i)
	{
		for (unsigned j = i + 1; j < nCols; ++j)
			if (fabs(data[i * nCols + j] - data[j * nCols + i]) > 1e-5)
			{
				cout << "Oh shit!: " << name << " is not symmetric" << endl;

				return;
			}
	}
	cout << name << " is symmetric" << endl;
}

void writeMatrixInMatlabFormat(	ostream& s,
								const string& name,
								const double* data,
								const unsigned m,
								const unsigned n
								)
{
	if (!m && !n)
		return;

	s << name << " = [ ";
	if (m == 1 || n == 1)
	{
		unsigned k = (m == 1) ? n : m;

		for (unsigned i = 0; i < k; ++i)
			s << data[ i ] << " ";
		s << " ]';" << endl << endl;
	}
	else
	{
		for (unsigned i = 0; i < m; ++i)
		{
			for (unsigned j = 0; j < n; ++j)
				s << data[j * m + i] << " ";
			s << "; ";
		}
		s << "];" << endl << endl;
	}
}

int main(int argc, char * const argv[ ])
{
	if (argc != 2)
	{
		cout << "Run index is required!" << endl;
		
		return 1;
	}
	
	const int runIndex = atoi( static_cast< char* >(argv[ 1 ]) );
	
	int i, iter, j;
	real_t measurement[ NX ];

	// Reset all memory, just for safety
	memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
	memset(&acadoVariables, 0, sizeof( acadoVariables ));
	
	int status;
	
	vector< vector< double > > eq;
	vector< vector< double > > dist;
	
	//
	// Initialize logger
	//
	
	stringstream outFile;
	outFile << "simple_compressor_log_N" << N << "_run" << runIndex << ".txt";
	
	cout << outFile.str() << endl;

	vector< vector< double > > log;
	log.resize( NUM_STEPS );
	for(unsigned i = 0; i < NUM_STEPS; ++i)
		log[ i ].resize(NX + NU + 6, 0.0);

	//
	// Timing stuff
	//
	real_t t1, t2, t3, t4;
	real_t fdbSum = 0.0;
	real_t intSum = 0.0;
	real_t condSum = 0.0;
	real_t prepSum = 0.0;

	timer t;

	//
	// Initialize the solver
	//
	initializeSolver();
		
	//
	// Initialize shooting nodes
	//
// 	for (i = 0; i < N + 1; ++i)
// 		for (j = 0; j < NX; ++j)
// 			acadoVariables.x[i * NX + j] = xxx;

	// V2
		
	acadoVariables.x[ 0 ] = 0.941e+000;
	acadoVariables.x[ 1 ] = 1.121e+000;
	acadoVariables.x[ 2 ] = 0.307;
	acadoVariables.x[ 3 ] = 1.191;
	acadoVariables.x[ 4 ] = 0.498e+000;
	acadoVariables.x[ 5 ] = 0.0;
	
	acadoVariables.p[0] = 35.0;
	acadoVariables.p[1] =  1.0;
// 	acadoVariables.p[2] =  0.7;
	acadoVariables.p[2] =  0.15;

	initializeNodesByForwardSimulation();

	//
	// Initialize references
	//
	
	for( i = 0; i < N; ++i )
	{
		acadoVariables.y[i * NY + 0] = 0.941e+000;
		acadoVariables.y[i * NY + 1] = 1.121e+000;
		acadoVariables.y[i * NY + 2] = 0.307;
		acadoVariables.y[i * NY + 3] = 1.191;
		acadoVariables.y[i * NY + 4] = 0.498e-001;
		acadoVariables.y[i * NY + 5] = 0.0;
		acadoVariables.y[i * NY + 6] = 0.0;
		acadoVariables.y[i * NY + 7] = 0.0;
	}
	  
	acadoVariables.yN[ 0 ] = 0.941e+000;
	acadoVariables.yN[ 1 ] = 1.121e+000;
	acadoVariables.yN[ 2 ] = 0.307;
	acadoVariables.yN[ 3 ] = 1.191;
	acadoVariables.yN[ 4 ] = 0.498e-001;
	acadoVariables.yN[ 5 ] = 0.0;
	
	//
	// Initial feedback
	//
	for (i = 0; i < NX; ++i)
		acadoVariables.x0[ i ] = acadoVariables.x[ i ]; // * 1.1;

	// Number of AS changes
	unsigned nIt = 0;
		
	//
	// Main loop
	//
	for( iter = 0; iter < NUM_STEPS; iter++ )
	{
#if ACADO_QP_SOLVER == ACADO_FORCES
		tic( &t );
		modelSimulation();
		evaluateObjective();
		t1 = toc( &t );

		tic( &t );
		evaluateConstraints();
		t2 = toc( &t );

		tic( &t );
		status = feedbackStep( );
		t3 = toc( &t );
		
		t4 = t3;

		nIt = acadoForces_info.it;

#elif ACADO_QP_SOLVER == ACADO_QPOASES
		tic( &t );
		modelSimulation();
		evaluateObjective();
		t1 = toc( &t );

		tic( &t );
		condensePrep(  );
		t2 = toc( &t );

		tic( &t );
		status = feedbackStep( );
		t3 = toc( &t );

		t4 = t3;

		nIt = getNWSR();

#elif ACADO_QP_SOLVER == ACADO_QPDUNES

		tic( &t );
		preparationStep();
		t1 = toc( &t );

		tic( &t );
		t2 = toc( &t );

		tic( &t );
		status = feedbackStep( );
		t3 = toc( &t );

		t4 = t3;

		nIt = 0;

#else
#error "Unknown QP solver option."
#endif

#if ACADO_QP_SOLVER == ACADO_FORCES
		if (status != 1)
#elif ACADO_QP_SOLVER == ACADO_QPOASES
		if ( status )

#elif ACADO_QP_SOLVER == ACADO_QPDUNES
		if (status != 1)

#else
#error "Wrong option for QP solver"
#endif
		{
			cout << "Iteration:" << iter << ", QP problem! QP status: " << status << endl;
			
			break;
		}
		
		for (unsigned ii = 0; ii < (N + 1) * NX; ++ii)
			if (acadoVariables.x[ ii ] != acadoVariables.x[ ii ])
		{
				cout << "Iteration:" << iter << ", NaN problems with diff. variables" << endl;

				exit( 1 );
		}

		for(i = 0; i < NX; i++)
			log[ iter ][ i ] = acadoVariables.x[ i ];
		for(j = 0;  i < NX + NU; i++, j++)
			log[ iter ][ i ] = acadoVariables.u[ j ];

		log[ iter ][ i++ ] = t1;
		log[ iter ][ i++ ] = t2;
		log[ iter ][ i++ ] = t3;
		log[ iter ][ i++ ] = t4;
		log[ iter ][ i++ ] = getKKT();
		log[ iter ][ i++ ] = (real_t)nIt;

		for (i = 0; i < NX; ++i)
			acadoVariables.x0[ i ] = acadoVariables.x[NX + i];

#if ACADO_QP_SOLVER == ACADO_QPDUNES
		shiftQpData();
#endif

 		shiftStates(2, 0, 0);
 		shiftControls( 0 );
		
		prepSum += t1 + t2;
		intSum += t1;
		condSum += t2;
		fdbSum += t3;
	}
	
#if ACADO_QP_SOLVER == ACADO_QPDUNES
	cleanupSolver();
#endif

#if VERBOSE

	cout << "Average fdbTime: " << scientific << fdbSum / NUM_STEPS * 1e6 << "usec" << endl;
	cout << "Average intTime: " << scientific << intSum / NUM_STEPS * 1e6 << "usec" << endl;
	cout << "Average condTime: " << scientific << condSum / NUM_STEPS * 1e6 << "usec" << endl;
	cout << "Average prepTime: " << scientific << prepSum / NUM_STEPS * 1e6 << "usec" << endl;
	
#endif

	ofstream dataLog( outFile.str().c_str() );
	if ( dataLog.is_open() )
	{
		for (i = 0; i < log.size(); i++)
		{
			for (unsigned j = 0; j < log[ i ].size(); j++)
				dataLog << log[ i ][ j ] << " ";
			dataLog << endl;
		}

		dataLog.close();
	}
	else
	{
		cout << "Log file could not be opened" << endl;

		return 1;
	}

    return 0;
}


