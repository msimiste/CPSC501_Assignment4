
/*  Include files  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>



using namespace std;

int getNextPowOf2(int input);

/*****************************************************************************
*
*    Function:     main
*
*    Description:  Tests the convolve function with various input signals
*
*****************************************************************************/

int main(int argc, char *argv[])
{
	int t = getNextPowOf2(63);
	cout << t;
	
  return 0;
}

int getNextPowOf2(int input){

	int i = 1;
	while(i < input)
	{
		i <<=1;
	}
	return i;
}
