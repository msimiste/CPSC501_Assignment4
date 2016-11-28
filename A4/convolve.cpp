
/*  Include files  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;


struct wavInfo{
	
unsigned int fileSize;
unsigned int subChunk1Size;
unsigned int numChannels;
unsigned int sampleRate;
unsigned int sampleSize;

float arr[];

	
} xWav,hWav;



/*  Function prototypes  */
void convolve(float x[], int N, float h[], int M, float y[], int P);
void print_vector(char *title, float x[], int N);
short int* readWavFile(char *inputFileName, wavInfo &wav);


/*****************************************************************************
*
*    Function:     main
*
*    Description:  Tests the convolve function with various input signals
*
*****************************************************************************/

int main(int argc, char *argv[])
{
	
	streampos size;
    char *inputFileName = argv[1];
	char *outputFilename = argv[3];
	char *IRFilename = argv[2];

	//short *x_temp = readWavFile(inputFileName,xWav);
	short *h_temp = readWavFile(IRFilename,hWav);

	for(int i = 0; i < hWav.fileSize/2; i++)
	{
		cout << hWav.arr[i] << " ";
	} 

  
  return 0;
}




short * readWavFile(char *inputFileName, wavInfo &wav){
	
    streampos size;
    char * memblock;
    int fileSizeOffset;
  
	ifstream file ( inputFileName, ios::in|ios::binary|ios::ate);
	
	
	
    if (file.is_open())
    {
        size = file.tellg();
        memblock = new char [size];
        file.seekg (0, ios::beg);
        file.read (memblock, size);
        file.close();

        //subChunk1Size = memblock[16];
        memcpy(&wav.subChunk1Size, &memblock[16],4);
    
        fileSizeOffset = 40 + (wav.subChunk1Size - 16);
    
        //cout << ("offset: ") << fileSizeOffset; 
    
        //cout << memblock[40];
        memcpy (&wav.fileSize,&memblock[fileSizeOffset],4);
        short * outArr = new short[wav.fileSize/2];
        int fileSize1 =  memblock[41];
        int fileSize2 = memblock[42];
    
    
    
    
        cout << "filesize :" << (wav.fileSize) << "\n";
        cout << "subChunkSize: " << (wav.subChunk1Size) << "\n";
        //cout << (25) << "\n";
        //cout << "the entire file content is in memory";
        int t = 0;
        for(int i = (fileSizeOffset + 4); i < wav.fileSize; i+=2){	
                
            short left = (memblock[i] << 8);
            short right = memblock[i+1];
            short combo = left & right;
            if(combo == 0){
            wav.arr[t] = combo;
            }            
            else if(combo <= 0){
                /* float tmp = (float)(combo >> 16);
                wav.arr[t] = (((float) (combo))/((float) (-32768)));
                cout << wav.arr[t] << " "; */
            }
            else{
               //wav.arr[t] = (((float) (combo))/((float) (32767)));
               wav.arr[t] = (float) (combo >> 15);
               cout << wav.arr[t] << " ";
            } 
            
            cout << wav.arr[t];
            outArr[t++] = combo;
            //cout << outArr[t++] << " ";
        } 

        return &outArr[1];
    }
    else cout << "Unable to open file";
}
/*****************************************************************************
*
*    Function:     convolve
*
*    Description:  Convolves two signals, producing an output signal.
*                  The convolution is done in the time domain using the
*                  "Input Side Algorithm" (see Smith, p. 112-115).
*
*    Parameters:   x[] is the signal to be convolved
*                  N is the number of samples in the vector x[]
*                  h[] is the impulse response, which is convolved with x[]
*                  M is the number of samples in the vector h[]
*                  y[] is the output signal, the result of the convolution
*                  P is the number of samples in the vector y[].  P must
*                       equal N + M - 1
*
*****************************************************************************/

void convolve(float x[], int N, float h[], int M, float y[], int P)
{
  int n, m;

  /*  Make sure the output buffer is the right size: P = N + M - 1  */
  if (P != (N + M - 1)) {
    printf("Output signal vector is the wrong size\n");
    printf("It is %-d, but should be %-d\n", P, (N + M - 1));
    printf("Aborting convolution\n");
    return;
  }

  /*  Clear the output buffer y[] to all zero values  */  
  for (n = 0; n < P; n++)
    y[n] = 0.0;

  /*  Do the convolution  */
  /*  Outer loop:  process each input value x[n] in turn  */
  for (n = 0; n < N; n++) {
    /*  Inner loop:  process x[n] with each sample of h[]  */
    for (m = 0; m < M; m++)
      y[n+m] += x[n] * h[m];
  }
}


/*****************************************************************************
*
*    Function:     print_vector
*
*    Description:  Prints the vector out to the screen
*
*    Parameters:   title is a string naming the vector
*                  x[] is the vector to be printed out
*                  N is the number of samples in the vector x[]
*
*****************************************************************************/

void print_vector(char *title, float x[], int N)
{
  int i;

  printf("\n%s\n", title);
  printf("Vector size:  %-d\n", N);
  printf("Sample Number \tSample Value\n");
  for (i = 0; i < N; i++)
    printf("%-d\t\t%f\n", i, x[i]);
}
