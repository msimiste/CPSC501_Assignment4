
/*  Include files  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

/*  CONSTANTS  ***************************************************************/
#define PI                3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY         440.0

/*  Test tone duration in seconds  */
#define DURATION          2.0				

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2

using namespace std;


struct wavInfo{
	
char * ChunkID;
unsigned int ChunkSize;
char * Format;
char * SubChunkID;
unsigned int subChunk1Size;

int_least16_t AudioFormat;
unsigned int numChannels;
unsigned int sampleRate;
unsigned int ByteRate;
unsigned int BlockAlign;
unsigned int sampleSize;

char * SubChunk2ID;
unsigned int Subchunk2Size;


unsigned int fileSize;
vector<float> arr;	
vector<int> arr1;

} xWav,hWav;



/*  Function prototypes  */
void convolve(vector<float> x, int N, vector<float> h, int M, float y[], int P);
void print_vector(char *title, float x[], int N);
short int* readWavFile(char *inputFileName, wavInfo &wav);
void fillFloatArray(short *	inWav, wavInfo& wav);
void fillIntArray(float y[], int *out, int p);
void createTone(double frequency, double duration,
                    int numberOfChannels, int numberOfSamples, int out[], char *filename);
void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile);
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);

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

	for(int i = 0; i < hWav.arr.size(); i++){
		hWav.arr.at(i) = 1.0;
	}
	fillFloatArray(x_temp, xWav);
	fillFloatArray(h_temp, hWav);
	
	//for(int i = 0; i < hWav.arr.size(); i++){
	
	//		cout << hWav.arr.at(i) << " ";
	//}  
	
	int p = xWav.arr.size() + hWav.arr.size() - 1;
	float* y = new float[p];
	
	
	convolve(xWav.arr, xWav.arr.size(), hWav.arr, hWav.arr.size(), y, p);
	
	int outVals[p];
	
	fillIntArray(y,outVals,p);
	
	//float* t = new float[xWav.arr.size()];
	//for(int i = 0; i < xWav.arr.size(); i++){ t[i] = xWav.arr.at(i);}
	
	//fillIntArray(t, outVals, xWav.arr.size());
	
	for(int i = 0; i < p; i++){
		cout << outVals[i] << " ";
	}
	createTone(FREQUENCY, DURATION, MONOPHONIC,p, outVals,outputFilename);
	
  return 0;
}



void fillFloatArray(short  * inWav, wavInfo& wav){

	for(int i = 0; i < wav.fileSize/2; i++)
	{
		if(inWav[i] == 0){ wav.arr.push_back(inWav[i]);}
		else if(inWav[i] > 0){
			wav.arr.push_back(((float)inWav[i]/(float)32767));
		}
		else{
			wav.arr.push_back(((float)inWav[i]/(float)32768));
		}		
	}
}

void fillIntArray(float y[], int *out, int p){
	for(int i = 0; i < p; i++){
		if(y[i] >= 0 ){
			out[i] = (int) (y[i] * 32767);
		}
		else{
				out[i] = (int)(y[i] * 32768);
		}		
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

		//get ChunkId value
        memcpy(&wav.ChunkID, &memblock[0], 4);
        
        //get ChunkSize value
        memcpy(&wav.ChunkSize, &memblock[4],4);      
        
        //get Format value
        memcpy(&wav.Format, &memblock[8],4);
        
        //getSubChunk1ID
        memcpy(&wav.SubChunkID, &memblock[12],4);
        
        //get SubChunk1Size value
        memcpy(&wav.subChunk1Size, &memblock[16],4);
        
        //get AudioFormat
        memcpy(&wav.AudioFormat, &memblock[18], 2);
        
        //get NumChannels
        memcpy(&wav.numChannels, &memblock[20], 2);
        
        //get SampleRate
        memcpy(&wav.sampleRate, &memblock[22],4);    
        
    
		
    
        fileSizeOffset = 40 + (wav.subChunk1Size - 16); 
       
        memcpy (&wav.fileSize,&memblock[fileSizeOffset],4);
       
        
        short * outArr = new short[wav.fileSize/2];
        int fileSize1 =  memblock[41];
        int fileSize2 = memblock[42];
    
        cout << "filesize :" << (wav.fileSize) << "\n";
        cout << "subChunkSize: " << (wav.subChunk1Size) << "\n";
     
        int t = 0;
        for(int i = (fileSizeOffset + 4); i < wav.fileSize; i+=2){	
                
            short left = (memblock[i]) | 0xff00;
				
            short right = (memblock[i+1] << 8 ) | 0x00ff;
            
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

//void convolve(float x[], int N, float h[], int M, float y[], int P)

void convolve(vector<float> x, int N, vector<float> h, int M, float y[], int P)
{
  int n, m;
	float max = 0.0;
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
    for (m = 0; m < M; m++){
      y[n+m] += x.at(n) * h.at(m);      
      if(abs(y[n+m]) > max){
		max = abs(y[n+m]);
		
	   }
     
	}
  }
  
  for(int i = 0; i < P; i++){
	
	y[i] = y[i]/max;
	cout << y[i] << " ";
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


/******************************************************************************
*
*       function:       createTone
*
*       purpose:        Calculates and writes out a sine tone to file
*
*       arguments:      frequency:  frequency of the test tone in Hz
*                       duration:  length of the test tone in seconds
*                       numberOfChannels:  number of audio channels
*                       filename:  name of the file to create
*                       
*       internal
*       functions:      writeWaveFileHeader, fwriteShortLSB
*
*       library
*       functions:      ceil, pow, fopen, fprintf, sin, rint, fclose
*
******************************************************************************/

void createTone(double frequency, double duration,
                    int numberOfChannels, int numberOfSamples, int out[], char *filename)
{
    int i;
		
    /*  Calculate the number of sound samples to create,
        rounding upwards if necessary  */
    //int numberOfSamples = (int)ceil(duration * SAMPLE_RATE);
	
    /*  Calculate the maximum value of a sample  */
    //int maximumValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;
	
    /*  Open a binary output file stream for writing */
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
        fprintf(stderr, "File %s cannot be opened for writing\n", filename);
        return;
    }
	
    /*  Write the WAVE file header  */
    writeWaveFileHeader(numberOfChannels, numberOfSamples,
                        SAMPLE_RATE, outputFileStream);
						
    /*  Create the sine tone and write it to file  */
    /*  Since the frequency is fixed, the angular frequency
        and increment can be precalculated  */
    double angularFrequency = 2.0 * PI * frequency;
    double increment = angularFrequency / SAMPLE_RATE;
    for (i = 0; i < numberOfSamples; i++) {
        /*  Calculate the sine wave in the range -1.0 to + 1.0  */
        double value = sin(i * increment);
		
        /*  Convert the value to a 16-bit integer, with the
            range -maximumValue to + maximumValue.  The calculated
            value is rounded to the nearest integer  */
        //short int sampleValue = rint(value * maximumValue);
		
        /*  Write out the sample as a 16-bit (short) integer
            in little-endian format  */
        fwriteShortLSB(out[i], outputFileStream);
		
        /*  If stereo output, duplicate the sample in the right channel  
        if (numberOfChannels == STEREOPHONIC)
            fwriteShortLSB(sampleValue, outputFileStream);*/
    }
						
    /*  Close the output file stream  */
    fclose(outputFileStream);
}



/******************************************************************************
*
*       function:       writeWaveFileHeader
*
*       purpose:        Writes the header in WAVE format to the output file.
*
*       arguments:      channels:  the number of sound output channels
*                       numberSamples:  the number of sound samples
*                       outputRate:  the sample rate
*                       outputFile:  the output file stream to write to
*                       
*       internal
*       functions:      fwriteIntLSB, fwriteShortLSB
*
*       library
*       functions:      ceil, fputs
*
******************************************************************************/

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}



/******************************************************************************
*
*       function:       fwriteIntLSB
*
*       purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

