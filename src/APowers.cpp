/****************************************************************/
/* Sequential and Parallel Sparse Matrix Multiplication Library */
/* version 2.3 -------------------------------------------------*/
/* date: 01/18/2009 --------------------------------------------*/
/* author: Aydin Buluc (aydin@cs.ucsb.edu) ---------------------*/
/* Customized to generate A3 by Abhinav Sarje <asarje@lbl.gov> -*/
/****************************************************************/

#include <mpi.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>  // Required for stringstreams
#include <ctime>
#include <cmath>
#include <SpParVec.h>
#include <SpTuples.h>
#include <SpDCCols.h>
#include <SpParMat.h>
#include <DenseParMat.h>
#include <DenseParVec.h>
#include <Operations.h>
#include <CombBLAS.h>

using namespace std;


// Simple helper class for declarations: Just the numerical type is templated 
// The index type and the sequential matrix type stays the same for the whole code
// In this case, they are "int" and "SpDCCols"
template <class NT>
class Dist
{ 
public: 
	typedef SpDCCols < int, NT > DCCols;
	typedef SpParMat < int, NT, DCCols > MPI_DCCols;
	typedef DenseParVec < int, NT> MPI_DenseVec;
};

template <class NT>
class PSpMat
{
public:
	typedef SpDCCols < int, NT > DCCols;
	typedef SpParMat < int, NT, DCCols > MPI_DCCols;
};


typedef PlusTimesSRing<double, double> PTDOUBLEDOUBLE;

int main(int argc, char* argv[])
{
	int nprocs, myrank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	typedef PlusTimesSRing<bool, bool> PTBOOLBOOL;
	if(argc < 4)
		{
		if(myrank == 0)
		{	
					cout << "Usage: ./apowers <INPUTADDRESS> <KLIMIT> <NONUM>" << endl;
					cout << "Example: ./apowers Data/input.txt 4 0" << endl;
			cout << "If <NONUM> is set, then the input file is unweighted" << endl;
				}
		MPI_Finalize();
		return -1;
		}

	{
		double klimit = atof(argv[2]);
		string ifilename(argv[1]);		
		ifstream input(ifilename.c_str());
		if( !input ) 
		{
				SpParHelper::Print( "Error opening input stream\n");
				return -1;
  		}
		MPI_Barrier(MPI_COMM_WORLD);
	
		Dist<double>::MPI_DCCols A;	// construct object
		bool nonum = static_cast<bool>(atoi(argv[3]));
		if(myrank ==0)
			cout << "Reading file " << (nonum?string("without"):string("with")) << " numerics" << endl;

		A.ReadDistribute(ifilename.c_str(), 0, nonum);	// read it from file
		if(myrank ==0)
			cout << "Matrix read"<< endl;	
		input.clear();
		input.close();

		int k = 1;
		std::string s;
		std::stringstream out;
		out << k;
		s = out.str();
		A.PrintForPatoh("A"+s+".patoh");

		SpParHelper::Print("Original matrix printed...\n");
		A.PrintInfo();
		if(klimit == 3) {
			Dist<double>::MPI_DCCols B = A;
			A.Square<PTBOOLBOOL>();
			//Dist<double>::MPI_DCCols out = Mult_AnXBn_Synch<SR, NUO, UDERO> (A, B, false, false);
			//Dist<double>::MPI_DCCols out = Mult_AnXBn_Synch<Dist<double>::MPI_DCCols> (A, B, false, false);
			PSpMat<double>::MPI_DCCols out = Mult_AnXBn_Synch<PTDOUBLEDOUBLE, double, PSpMat<double>::DCCols >(A, B);

			A.PrintForPatoh("A2.patoh");
			string message = "Power 2 printed\n"; 
			SpParHelper::Print(message);
			A.PrintInfo();

			out.PrintForPatoh("A3.patoh");
			message = "Power 3 printed\n"; 
			SpParHelper::Print(message);
			out.PrintInfo();
		} else {
			while (k < klimit)
			{
				A.Square<PTBOOLBOOL>();
				k *= 2;
			
				std::string s;
				std::stringstream out;
				out << k;
				s = out.str();
				A.PrintForPatoh("A"+s+".patoh");

				string message = "Power " + s + " printed\n"; 
				SpParHelper::Print(message);
				A.PrintInfo();
			}
		} // if-else
	}	

	// make sure the destructors for all objects are called before MPI::Finalize()
	MPI_Finalize();
	return 0;
}

