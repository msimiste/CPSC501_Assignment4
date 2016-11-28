// reading an entire binary file
#include <iostream>
#include <fstream>
using namespace std;

int main () {
  streampos size;
  char * memblock;

  ifstream file ("somefile.wav", ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();

    //cout << "the entire file content is in memory";
    
   for(int i = 0; i< 100; i++){
		cout << memblock[i];
	} 

    delete[] memblock;
  }
  else cout << "Unable to open file";
  return 0;
}
