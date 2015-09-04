// Distributed two-dimensional Discrete FFT transform
// Guru Das Srinagesh
// ECE 6122 Project 1


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <signal.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "Complex.h"
#include "InputImage.h"

using namespace std;

void Transform1D(Complex* h, int w, Complex* H, int startrow, int rows_per_CPU)
{
  // Implement a simple 1-d DFT using the double summation equation
  // given in the assignment handout.  h is the time-domain input
  // data, w is the width (N), and H is the output array.

  int N = w;    // For convenience's sake
  int row;
  
  Complex W;    // This is the W in the formula
  
  Complex sum;  // Sum of all W[k]*h[k] inner loop products

  /* Build the 1-D transform formula */
  int n, k;
  // cout<<"Rank - "<<(startrow/rows_per_CPU)<<"  1-D transform computation\t\t";
  printf("  (Rank %2d)  1-D transform computation\t\t", (startrow/rows_per_CPU));

  for(row=startrow; row<(startrow + rows_per_CPU); row++)
    for(n=N*row; n<N*(row+1); n++){
      sum.real = 0; sum.imag = 0;
      for(k=0; k<N; k++){
        W.real = cos(2*M_PI*n*k/N); W.imag = -sin(2*M_PI*n*k/N);
        sum = sum + (W * h[(N*row) + k]);
      }
      H[n].real = sum.real;
      H[n].imag = sum.imag;
    } 

  cout<<"[ OK ]"<<endl<<endl;
}

void transpose(Complex* mat, int N){
  Complex temp;

  for(int row=0; row<N; row++)
    for(int column=0; column<N; column++)
        if(row < column){
          temp= mat[N*row + column]; 
          mat[N*row + column]= mat[N*column + row];
          mat[N*column + row] = temp;
        }
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Use MPI to find how many CPUs in total, and which one
  //    this process is
  // 3) Allocate an array of Complex object of sufficient size to
  //    hold the 2d DFT results (size is width * height)
  // 4) Obtain a pointer to the Complex 1d array of input data
  // 5) Do the individual 1D transforms on the rows assigned to your CPU
  // 6) Send the resultant transformed values to the appropriate
  //    other processors for the next phase.
  // 6a) To send and receive columns, you might need a separate
  //     Complex array of the correct size.
  // 7) Receive messages from other processes to collect your columns
  // 8) When all columns received, do the 1D transforms on the columns
  // 9) Send final answers to CPU 0 (unless you are CPU 0)
  //   9a) If you are CPU 0, collect all values from other processors
  //       and print out with SaveImageData().

  /* Step 1:  Read in the image and find out width and height */
  InputImage image(inputFN);  // Create the helper object for reading the image
  int width = image.GetWidth();
  int height = image.GetHeight();

  /* Step 2:  Use MPI to determine how many CPUs there are, etc. */
  int nCPUs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nCPUs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* ----- :  Calculate number of rows per CPU and startrow */

  //rank = 15; /******* TODO TODO TODO Remove this! This is just for testing! *******/

  int rows_per_CPU = height / 16/******* TODO nCPUs ********/;    //cout<<"rows_per_CPU = "<<rows_per_CPU<<endl;
  int startrow = rows_per_CPU * rank;                             //cout<<"startrow = "<<startrow<<endl;
  
  /* Step 3:  Create the H array that contains the 2-D DFT results */
  Complex* H = new Complex[width * height];

  /* Step 4:  Get the data of the input image, namely create the h array */
  Complex* h = image.GetImageData();

  /* Step 5:  Do 1-D transform of input image */
  // If you are rank 0, then no need to send anything, just calculate your bit
  // and sit tight and quietly receive everything from everyone else and put
  // the data in the proper rows of the H array.
  Transform1D(h, 256, H, startrow, rows_per_CPU);

  // // REMOVE THIS!!!!!!!!!!!!!
  // int a[16];
  // for(int i=0; i<16; i++)
  //   a[i] = i; 

  if(rank != 0/*== 5*/){
    /* Advance H pointer to the right place for the purpose of sending to CPU 0 */
    // H = H + (startrow * 256);
    int rc = MPI_Send(&H[startrow*256], rows_per_CPU * 256 * sizeof(Complex), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    if(rc != MPI_SUCCESS){
      cerr<<"Error sending data to CPU 0. Exiting.\n";
      MPI_Finalize(); exit(1);
    }
    // 
    // ---------------------------------------------------------------------------------------------
    // 
    
    // if(rank == 1 || rank == 3){
    //   int rc = MPI_Send(&a[rank], 2*sizeof(int), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    //   if(rc != MPI_SUCCESS){
    //     cerr<<"Error sending data to CPU 0. Exiting.\n";
    //     MPI_Finalize(); exit(1);
    //   }
    // }

  } else if (rank == 0){
    MPI_Status status; int rc;
    MPI_Request request;
    for(int cpu=1; cpu<nCPUs; cpu++){
      // int cpu = 5; // CHAAAAAAAAAAAAAAAANNNNNNNNNNNGEEEEEEEEEEEE this!
      // H = H + (rows_per_CPU * 256 * (cpu-1)); 
      rc = MPI_Irecv(&H[rows_per_CPU*cpu*256], rows_per_CPU * 256 * sizeof(Complex), MPI_CHAR, cpu, 0, MPI_COMM_WORLD, &request);
      if(rc != MPI_SUCCESS){
        cerr<<"Error receiving data from CPU "<<cpu<<". Exiting.\n";
        MPI_Finalize(); exit(1);
      }
      MPI_Wait(&request, &status);
      printf("  (Rank %2d)  Part of 1-D X-form received \t[ OK ]\n\n",cpu);
    }
    /* Step VV: Save 1-D transform image data for debug */
    image.SaveImageData("MyAfter1d.txt", H, width, height);
    cout<<"  ---- File written ----: MyAfter1d.txt\t\t[ OK ]"<<endl<<endl;
    // 
    // -----------------------------------------------------------------------------------
    // 
    // int b[1]; MPI_Status status; MPI_Request request;
    // for(int cpu = 1; cpu<nCPUs; cpu++){
    //   int rc = MPI_Irecv(b, sizeof(double), MPI_INT, cpu, 0, MPI_COMM_WORLD, &request);
    //   if(rc != MPI_SUCCESS){
    //     cerr<<"Error receiving data from CPU "<<cpu<<". Exiting.\n";
    //     MPI_Finalize(); exit(1);
    //   }
    //   MPI_Wait(&request, &status);
    //   cout<<"  Received data: "<<b[0]<<" from CPU "<<cpu<<endl;
    // }
    // 
    // ------------------------------------------------------------------------------------
    // 

    // int *b = (int *) calloc(16,sizeof(int)); b[0] = 0;
    // MPI_Status status; MPI_Request request;
    // for(int cpu = 1; cpu<5; cpu+=2){
    //   // b = b + cpu;
    //   int rc = MPI_Irecv(&b[cpu], 2*sizeof(int), MPI_CHAR, cpu, 0, MPI_COMM_WORLD, &request);
    //   if(rc != MPI_SUCCESS){
    //     cerr<<"Error receiving data from CPU "<<cpu<<". Exiting.\n";
    //     MPI_Finalize(); exit(1);
    //   }
    //   MPI_Wait(&request, &status);
    //   cout<<"  Received data: "<<b[cpu]<<" from CPU "<<cpu<<endl;
    // }

    // cout<<"Final result = ";
    // for(int cpu = 0; cpu<nCPUs; cpu++){
    //   cout<<b[cpu]<<" ";
    // }
    // cout<<endl;

  }

  

  ///* Step WW: Do Transpose of intermediate 1-D transform array */
  //transpose(H, 256);
  //cout<<"  Post 1-D transpose \t\t\t[ OK ]"<<endl<<endl;

  ///* Step XX: Do 1-D transform again on the transpose array */
  //Complex* H_final = new Complex[width * height];
  //Transform1D(H, 256, H_final);

  ///* Step YY: Do Transpose of second intermediate H_final array for final result */
  //transpose(H_final, 256);
  //cout<<"  Post 2-D transpose \t\t\t[ OK ]"<<endl<<endl;

  ///* Step ZZ: Finally, write the 2-D transform values to disk */
  //image.SaveImageData("MyAfter2d.txt", H_final, width, height);
  //cout<<"  File written: MyAfter2d.txt\t\t[ OK ]"<<endl<<endl;
}


int main(int argc, char** argv)
{
  /* Do MPI Init here */
  int rc = MPI_Init(&argc, &argv);

  if(rc != MPI_SUCCESS){
    cerr<<"Error starting MPI program. Exiting.\n";
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
  // cout<<"Rank - "<<rank<<"  Exiting normally.\n"<<endl;
  printf("   Rank %2d   EXIT\t\t\t\t  --\n\n", rank);
  MPI_Finalize();
}  
