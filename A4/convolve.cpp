
/*  Include files  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;


struct wavInfo{
	
unsigned int fileSize;
unsigned int ChunkSize;
char * Format;
char * SubChunkID;
unsigned int subChunk1Size;
unsigned int numChannels;
unsigned int sampleRate;
unsigned int sampleSize;

vector<float> arr;

	
} xWav,hWav;



/*  Function prototypes  */
void convolve(float x[], int N, float h[], int M, float y[], int P);
void print_vector(char *title, float x[], int N);
short int* readWavFile(char *inputFileName, wavInfo &wav);
void fillFloatArray(short *	inWav, wavInfo& wav);


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
	short *x_temp = readWavFile(inputFileName,xWav); //143351
	short *h_temp = readWavFile(IRFilename, hWav);

	fillFloatArray(x_temp, xWav);
	fillFloatArray(h_temp, hWav);
	
	for(int i = 0; i < hWav.arr.size(); i++){
	
			cout << hWav.arr.at(i) << " ";
	}
	

  
  return 0;
}


void fillFloatArray(short  * inWav, wavInfo& wav){

	for(int i = 0; i < wav.fileSize/2; i++)
	{
		//hWav.arr.at(i) = 
		if(inWav[i] == 0){ wav.arr.push_back(inWav[i]);}
		else if(inWav[i] > 0){
			wav.arr.push_back(((float)inWav[i]/(float)32767));
		}
		else{
			wav.arr.push_back(((float)inWav[i]/(float)32768));
		}
		
		//cout << wav.arr.at(i) << " ";
		//cout << h_temp[i] <<  " ";
	}
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
        //short * outArr = new short[wav.fileSize/2];
        
        short * outArr = new short[wav.fileSize/2];
        int fileSize1 =  memblock[41];
        int fileSize2 = memblock[42];
        cout << "filesize :" << (wav.fileSize) << "\n";
        cout << "subChunkSize: " << (wav.subChunk1Size) << "\n";
      
        int t = 0;
        for(int i = (fileSizeOffset + 4); i < wav.fileSize; i+=2){	
                
            short left = (memblock[i] << 8);
            short right = memblock[i+1];
            short combo = left & right;
		    outArr[t++] = combo;
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
