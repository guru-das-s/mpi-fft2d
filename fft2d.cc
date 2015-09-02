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

#include "Complex.h"
#include "InputImage.h"

using namespace std;

void Transform1D(Complex* h, int w, Complex* H)
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
  cout<<"  1-D transform computation\t\t";

  for(row=0; row<N; row++)
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
  //          As an initial step, forget about MPI, just get 1-D transform right

  /* Step 3:  Create the H array that contains the 2-D DFT results */
  Complex* H = new Complex[width * height];

  /* Step 4:  Get the data of the input image, namely create the h array */
  Complex* h = image.GetImageData();

  /* Step 5:  Do 1-D transform of input image */
  //          As an initial step, forget about MPI, just get 1-D transform right
  cout<<endl; // For formatting output nicely 
  Transform1D(h, 256, H);

  /* Step VV: Save 1-D transform image data for debug */
  image.SaveImageData("MyAfter1d.txt", H, width, height);
  cout<<"  File written: MyAfter1d.txt\t\t[ OK ]"<<endl<<endl;

  /* Step WW: Do Transpose of intermediate 1-D transform array */
  int N = 256; Complex temp;
  
  for(int row=0; row<N; row++)
    for(int column=0; column<N; column++)
        if(row < column){
          //temp.real = H[N*row + column].real; temp.imag = H[N*row + column].imag; 
          //H[N*row + column].real = H[N*column + row].real; H[N*row + column].imag = H[N*column + row].imag;
          //H[N*column + row].real = temp.real; H[N*column + row].imag = temp.imag;

          temp= H[N*row + column]; 
          H[N*row + column]= H[N*column + row];
          H[N*column + row] = temp;
        }
  cout<<"  Post 1-D transpose \t\t\t[ OK ]"<<endl<<endl;

  /* Step XX: Do 1-D transform again on the transpose array */
  Complex* H_final = new Complex[width * height];
  Transform1D(H, 256, H_final);

  /* Step YY: Do Transpose of second intermediate H_final array for final result */
  for(int row=0; row<N; row++)
    for(int column=0; column<N; column++)
        if(row < column){
          temp= H_final[N*row + column]; 
          H_final[N*row + column]= H_final[N*column + row];
          H_final[N*column + row] = temp;
        }
  cout<<"  Post 2-D transpose \t\t\t[ OK ]"<<endl<<endl;

  /* Step ZZ: Finally, write the 2-D transform values to disk */
  image.SaveImageData("MyAfter2d.txt", H_final, width, height);
  cout<<"  File written: MyAfter2d.txt\t\t[ OK ]"<<endl<<endl;
}


int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
}  
