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

void Transform1D(Complex* h, int w, Complex* H, int startrow, int rows_per_CPU, int D)
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
  printf("  (Rank %2d)  %d-D transform computation\t\t", (startrow/rows_per_CPU), D);

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
  int rows_per_CPU = height / 16/******* TODO nCPUs ********/;
  int startrow = rows_per_CPU * rank;
  
  /* Step 3:  Create the H array that contains the 2-D DFT results */
  Complex* H = new Complex[width * height];

  /* Step 4:  Get the data of the input image, namely create the h array */
  Complex* h = image.GetImageData();

  /* Step 5:  Do 1-D transform of input image */
  Transform1D(h, 256, H, startrow, rows_per_CPU, 1);

  /* Step GG: Stitch together pieces of 1D X-formed image, gathering from all CPUs */
  if(rank != 0){
    /* Send your section of the 1-D transformed image to CPU 0 */
    int rc = MPI_Send(&H[startrow*256], rows_per_CPU * 256 * sizeof(Complex), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    if(rc != MPI_SUCCESS){
      cerr<<"Error sending data to CPU 0. Exiting.\n";
      MPI_Finalize(); exit(1);
    }

  } else if (rank == 0){
    /* Collect parts of the 1-D transformed image from rest of CPUs and stitch them together */
    MPI_Status status; int rc;
    MPI_Request request;
    for(int cpu=1; cpu<nCPUs; cpu++){
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

    /* Step WW: Do Transpose of intermediate 1-D transform array */
    transpose(H, 256);
    printf("  (Rank %2d)  Transpose of 1-D transformed image [ OK ]\n\n",rank);
  }
  
  /* Step MM: Send entire Transpose of 1-D transformed image to all other CPUs */
  if(rank==0){
    /* Send to all CPUs */
    MPI_Status status; MPI_Request request; int rc;
    for(int cpu=1; cpu<nCPUs; cpu++){
      rc = MPI_Isend(H, 256*256*sizeof(Complex), MPI_CHAR, cpu, 0, MPI_COMM_WORLD, &request);
      if(rc != MPI_SUCCESS){
        cerr<<"Error sending Tranpose of image to CPU "<<cpu<<". Exiting.\n";
        MPI_Finalize(); exit(1);
      }
      MPI_Wait(&request, &status);
    }    

  } else {
    /* Receive from CPU 0 */
    MPI_Status status; MPI_Request request; int rc;
    rc = MPI_Irecv(H, 256*256*sizeof(Complex), MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
    if(rc != MPI_SUCCESS){
      cerr<<"Error receiving Tranpose of image from CPU 0. Exiting.\n";
      MPI_Finalize(); exit(1);
    }
    MPI_Wait(&request, &status);    
  }

  /* Step XX: Do 1-D transform again on the transpose array using MPI */
  Complex* H_final = new Complex[width * height];
  Transform1D(H, 256, H_final, startrow, rows_per_CPU, 2);

  /* Step NN: Stitch together pieces of 2D X-formed image, gathering from all CPUs */
  if(rank != 0){
    /* Send your section of the 1-D transformed image to CPU 0 */
    int rc = MPI_Send(&H_final[startrow*256], rows_per_CPU * 256 * sizeof(Complex), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    if(rc != MPI_SUCCESS){
      cerr<<"Error sending data to CPU 0. Exiting.\n";
      MPI_Finalize(); exit(1);
    }

  } else if (rank == 0){
    /* Collect parts of the 1-D transformed image from rest of CPUs and stitch them together */
    MPI_Status status; MPI_Request request; int rc;
    for(int cpu=1; cpu<nCPUs; cpu++){
      rc = MPI_Irecv(&H_final[rows_per_CPU*cpu*256], rows_per_CPU * 256 * sizeof(Complex), MPI_CHAR, cpu, 0, MPI_COMM_WORLD, &request);
      if(rc != MPI_SUCCESS){
        cerr<<"Error receiving data from CPU "<<cpu<<". Exiting.\n";
        MPI_Finalize(); exit(1);
      }
      MPI_Wait(&request, &status);
      printf("  (Rank %2d)  Part of 1-D X-form received \t[ OK ]\n\n",cpu);
    }

    /* Step WW: Do Transpose of intermediate 1-D transform array */
    transpose(H_final, 256);
    printf("  (Rank %2d)  Transpose of 2-D transformed image [ OK ]\n\n",rank);

    /* Step VV: Save 2-D transform image data for debug */
    image.SaveImageData("MyAfter2d.txt", H_final, width, height);
    cout<<"  ---- File written ----: MyAfter2d.txt\t\t[ OK ]"<<endl<<endl;
  }
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

  string fn("Tower.txt");               // default file name
  if (argc > 1) fn = string(argv[1]);   // if name specified on cmd line
  Transform2D(fn.c_str());              // Perform the transform.

  printf("   Rank %2d   EXIT\t\t\t\t  --\n\n", rank);

  MPI_Finalize();

  return 0;
}  
