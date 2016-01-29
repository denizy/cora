/*Copyright (c) 2015-2016 Deniz Yorukoglu. All rights reserved.*/
#include<iostream>
#include<fstream>
#include<sstream>
#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<zlib.h>
using namespace std;

#define MAX_LINE_LEN 255 //limit for readLength
#define MAX_FRAG_GROUPS 255 //limit for the number of different fragmentation groups (though the total number of fragmentation signals could be larger)

unsigned int numGroups = 0;
unsigned int fragSignalSize = 0;
unsigned char fragHash_2[85][85]; //This is for ACGTN ^ 2 fast access of frag group
unsigned char fragHash_4[85][85][85][85]; // Same with ^ 4
int readLength;

//Convert from integer to DNA base
unsigned char GetNuc(unsigned char c)
{
	switch(c)
	{
		case 0:
			return 'N';
		case 1:
			return 'A';
		case 2:
			return 'C';
		case 3:
			return 'G';
		case 4:
			return 'T';
		default:
			cout << "ERROR: unexpected character: '" << c << "'" << endl;
			exit(8);
			break;
	}
}

//Set fragmentation signal to group id
void SetFragHash4(unsigned char i1, unsigned char i2, unsigned char i3, unsigned char i4, unsigned char val)
{
	unsigned char c1 = GetNuc(i1);
	unsigned char c2 = GetNuc(i2);
	unsigned char c3 = GetNuc(i3);
	unsigned char c4 = GetNuc(i4);

	fragHash_4[c1][c2][c3][c4] = val;
}

//Setup which frag signal belongs to which frag group
void SetupWholeFragmentation(ifstream& fin) //Currently all fragmentation is one character from each end (update this for further fragmentation for incredibly large whole genome runs)
{
	int numLines;
	fin >> numLines >> fragSignalSize;

	if(fragSignalSize == 2)
	{
		string fragModeStr;
		while(fin >> fragModeStr)
		{
			if(fragModeStr[0]=='N')
			{
				fragHash_2['N']['N'] = numGroups;
				fragHash_2['N']['A'] = numGroups;
				fragHash_2['N']['C'] = numGroups;
				fragHash_2['N']['G'] = numGroups;
				fragHash_2['N']['T'] = numGroups;
				fragHash_2['A']['N'] = numGroups;
				fragHash_2['C']['N'] = numGroups;
				fragHash_2['G']['N'] = numGroups;
				fragHash_2['T']['N'] = numGroups;

				fragModeStr = fragModeStr.substr(2, fragModeStr.length()-2);
			}
			
			assert(fragModeStr.length() % 3 == 0);
			int numFrags = fragModeStr.length() / 3;
			for(int i=0; i<numFrags; i++) //All reverse complements are already provided
			{
				fragHash_2[(unsigned char) fragModeStr[3*i]][(unsigned char) fragModeStr[3*i +1]] = numGroups;
			}
				
			numGroups++;
		}
	}
	else if(fragSignalSize == 4)
	{
		string fragModeStr;
		while(fin >> fragModeStr)
		{
			if(fragModeStr[0]=='N')
			{
				for(unsigned char i1=0; i1<=4; i1++) //Efficiency isn't much important in this part
					for(unsigned char i2=0; i2<=4; i2++)
						for(unsigned char i3=0; i3<=4; i3++)
							for(unsigned char i4=0; i4<=4; i4++)
							{
								if(i1 == 0 || i2 == 0 || i3 == 0 || i4 == 0)
								{
									SetFragHash4(i1, i2, i3, i4, numGroups);
								}
							}				
				fragModeStr = fragModeStr.substr(2, fragModeStr.length()-2);
			}
			
			assert(fragModeStr.length() % 5 == 0);
			int numFrags = fragModeStr.length() / 5;
			for(int i=0; i<numFrags; i++) //All reverse complements are already provided
			{
				fragHash_4
					[(unsigned char) fragModeStr[5*i]]
					[(unsigned char) fragModeStr[5*i +1]]
					[(unsigned char) fragModeStr[5*i + 2]]
					[(unsigned char) fragModeStr[5*i + 3]]  = numGroups;
			}
				
			numGroups++;
		}
	}
	else
	{
		cout << "fragsignalsize should be either 2 or 4" << endl;
		cout << "size: " << fragSignalSize << endl;
		exit(7);
	}
}

inline unsigned short GetFragGroupNo_2(unsigned char a, unsigned char b)
{
	return fragHash_2[a][b];
}

inline unsigned short GetFragGroupNo_4(unsigned char a, unsigned char b, unsigned char c, unsigned char d)
{
	return fragHash_4[a][b][c][d];
}

//////////////////////
// Variables & functions related to output buffering
//////////////////////

#define OUTPUT_BUFFER_SIZE 2000000
char* buffer[MAX_FRAG_GROUPS+2]; //buffers for printing output
unsigned int bufferLen[MAX_FRAG_GROUPS + 2]; //buffer lengths
FILE* fout[MAX_FRAG_GROUPS+2]; //output file list

void FlushBuffer(unsigned short group) //flush buffer for given frag group
{
	fwrite(buffer[group], 1, bufferLen[group], fout[group]);
	bufferLen[group] = 0;
}

//add read-mer to buffer for given group
void AddToBuffer(unsigned short group, unsigned long long code, char* seq, unsigned char len) //len is variable depending on whether it's split or not
{
	if(bufferLen[group] > OUTPUT_BUFFER_SIZE)
	{
		FlushBuffer(group);
	}	

	unsigned long long * ptr = (unsigned long long *)  (buffer[group] + bufferLen[group]);
	*ptr = code;

	bufferLen[group] += sizeof(unsigned long long);	

	for(int i=0; i<len; i++) //Copy the entire read sequence
	{
		buffer[group][bufferLen[group]++] = seq[i];
	}
}

//////////////////////
// Variables & functions related to input buffering
//////////////////////

#define INPUT_BUFFER_SIZE 2000000
char firstInBuffer[INPUT_BUFFER_SIZE+5000];
unsigned int firstInBufferSize = 0;
unsigned int firstInBufferOffset = 0;
unsigned int firstInBufferOffset_next = 0;

char secondInBuffer[INPUT_BUFFER_SIZE+5000];
unsigned int secondInBufferSize = 0;
unsigned int secondInBufferOffset = 0;
unsigned int secondInBufferOffset_next = 0;

//functions are self-explanatory
inline unsigned int skipNlines_INT_FirstBuffer(unsigned int curInd,  unsigned int skipCount)
{
	while(skipCount > 0 && curInd < firstInBufferSize)
	{
		if(firstInBuffer[curInd] == '\n')
		{
			skipCount--;
		}
		curInd++;
	}
	return curInd;
}

inline unsigned int skipNlines_INT_SecondBuffer(unsigned int curInd,  unsigned int skipCount)
{
	while(skipCount > 0 && curInd < secondInBufferSize)
	{
		if(secondInBuffer[curInd] == '\n')
		{
			skipCount--;
		}
		curInd++;
	}
	return curInd;
}

inline unsigned int skipUntil_Char_FirstBuffer(unsigned int curInd, char termChar)
{
	while(firstInBuffer[curInd] != termChar)
	{
		curInd++;
	}
	return curInd;
}

inline unsigned int skipUntil_Char_SecondBuffer(unsigned int curInd, char termChar)
{
	while(secondInBuffer[curInd] != termChar)
	{
		curInd++;
	}
	return curInd;
}

void InitFirstBuffer(FILE* fin)
{
	size_t result = fread(firstInBuffer, 1, INPUT_BUFFER_SIZE, fin);	

	firstInBufferOffset = 0;
	firstInBufferSize = result;

	firstInBufferOffset_next = skipNlines_INT_FirstBuffer(firstInBufferOffset, 1); //skips from read seq to name of the next read 
}

void InitSecondBuffer(FILE* fin)
{
	size_t result = fread(secondInBuffer, 1, INPUT_BUFFER_SIZE, fin);	

	secondInBufferOffset = 0;
	secondInBufferSize = result;

	secondInBufferOffset_next = skipNlines_INT_SecondBuffer(secondInBufferOffset, 1); //skips from read seq to name of the next read 
}

unsigned int Min_skip_len;

char* GetNextReadStart_FirstBuffer(FILE* fin)
{
	firstInBufferOffset = firstInBufferOffset_next;  //skips from readName to read seq
	firstInBufferOffset_next = skipNlines_INT_FirstBuffer(firstInBufferOffset + Min_skip_len, 1); //skips from read seq to name of the next read 

	if(firstInBufferOffset_next >= firstInBufferSize)
	{
		if(firstInBufferOffset >= firstInBufferSize)
		{
			return NULL; //there are no more to read
		}

		unsigned int numShiftedChars = firstInBufferSize - firstInBufferOffset;
		for(unsigned int i=0; i<numShiftedChars; i++)
		{
			firstInBuffer[i] = firstInBuffer[firstInBufferOffset + i];
		}
		size_t result = fread(firstInBuffer + numShiftedChars, 1, INPUT_BUFFER_SIZE - numShiftedChars, fin);	

		firstInBufferSize = result + numShiftedChars;

		firstInBufferOffset = 0;
		firstInBufferOffset_next = skipNlines_INT_FirstBuffer(firstInBufferOffset + Min_skip_len, 1); //skips from read seq to name of the next read 
	}

	return firstInBuffer + firstInBufferOffset;
}

char* GetNextReadStart_SecondBuffer(FILE* fin)
{
	secondInBufferOffset = secondInBufferOffset_next;  //skips from readName to read seq
	secondInBufferOffset_next = skipNlines_INT_SecondBuffer(secondInBufferOffset + Min_skip_len, 1); //skips from read seq to name of the next read 
	
	if(secondInBufferOffset_next >= secondInBufferSize)
	{
		if(secondInBufferOffset >= secondInBufferSize)
		{
			return NULL; //there are no more to read
		}

		unsigned int numShiftedChars = secondInBufferSize - secondInBufferOffset;
		for(unsigned int i=0; i<numShiftedChars; i++)
		{
			secondInBuffer[i] = secondInBuffer[secondInBufferOffset + i];
		}
		size_t result = fread(secondInBuffer + numShiftedChars, 1, INPUT_BUFFER_SIZE - numShiftedChars, fin);	

		secondInBufferSize = result + numShiftedChars;

		secondInBufferOffset = 0;
		secondInBufferOffset_next = skipNlines_INT_SecondBuffer(secondInBufferOffset + Min_skip_len, 1); //skips from read seq to name of the next read 
	}

	return secondInBuffer + secondInBufferOffset;
}

//////////////////////////////////////////


int main(int argc, char* argv[])
{
	if(argc != 6)
	{
		cout << "This program splits are given set of fastq files into a list of fragmented binary read-mer files (sequence is still text within binary)." << endl;
		cout << "ARGV[1] = fragmentation signals" << endl;
		cout << "ARGV[2] = readLen" << endl;
		cout << "ARGV[3] = Input FASTQ file list" << endl;
		cout << "ARGV[4] = Output prefix name (*.frag0 ....)" << endl; //Names are assigned as numbers (long long in binary) + Number is same as original --> lineNumber, then mate, then split, then direction
		cout << "ARGV[5] = mapping mode: SINGLE, SPLIT_SINGLE, PAIRED, SPLIT_PAIRED, THREEWAY_SINGLE, THREEWAY_PAIRED" << endl;
		exit(1);
	}
	string mappingModeStr = string(argv[5]);
	
	readLength = atoi(argv[2]);		
	Min_skip_len = 2*readLength + 5; //this is the minimum distance between the beginning of the seqeunce line and the past the beginning of the next read name;

	ifstream finFastqList(argv[3]);

	int numSamples;
	finFastqList >> numSamples;
	
	ifstream finFrags(argv[1]);
	SetupWholeFragmentation(finFrags);

	for(unsigned int i=0; i<numGroups; i++)
	{
		stringstream outputName;
		outputName << argv[4] << ".frag" << i;
		buffer[i] = (char *) calloc (OUTPUT_BUFFER_SIZE + 2000, 1); //2000 for any trailing sequence in the current buffer 
		fout[i] = fopen(outputName.str().c_str(), "wb");
	}
	
	unsigned long long lineCount;

	if(fragSignalSize == 2)
	{

		if(mappingModeStr == "SINGLE")
		{
			char* readLine;
			unsigned long long readCode = 0;
			string curFastqFileName;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName >> lineCount;

				FILE* finRead = fopen(curFastqFileName.c_str(), "rb");
				InitFirstBuffer(finRead);		
				while(true)
				{
					readLine = GetNextReadStart_FirstBuffer(finRead);
					
					if(readLine == NULL)
					{
						break;
					}
		
					unsigned short groupNo = GetFragGroupNo_2(readLine[0], readLine[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine, readLength);

					readCode += 2; //for direction jumps (they're all forward at this point)
				}	

				fclose(finRead);
			}
		}
		else if(mappingModeStr == "SPLIT_SINGLE")
		{
			char* readLine;
		
			unsigned long long readCode = 0;
		
			assert(readLength % 2 == 0);
			unsigned short halfReadLength = readLength / 2;

			string curFastqFileName;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName >> lineCount;
			
				FILE* finRead = fopen(curFastqFileName.c_str(), "rb");
				InitFirstBuffer(finRead);		
				while(true)
				{
					readLine = GetNextReadStart_FirstBuffer(finRead);

					if(readLine == NULL)
					{
						break;
					}

					unsigned short groupNo = GetFragGroupNo_2(readLine[0], readLine[halfReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine, halfReadLength);
					readCode += 2; //switching to second split

					//this is second split
					groupNo = GetFragGroupNo_2(readLine[halfReadLength], readLine[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine + halfReadLength, halfReadLength);
					readCode += 2; //switching to the first split of the next read
				}	
				fclose(finRead);
			}
		}
		else if(mappingModeStr == "THREEWAY_SINGLE")
		{
			char* readLine;
		
			unsigned long long readCode = 0;
		
			assert(readLength % 3 == 0);
			unsigned short oneThirdReadLength = readLength / 3;

			string curFastqFileName;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName >> lineCount;
			
				FILE* finRead = fopen(curFastqFileName.c_str(), "rb");
				InitFirstBuffer(finRead);		
				while(true)
				{
					readLine = GetNextReadStart_FirstBuffer(finRead);

					if(readLine == NULL)
					{
						break;
					}

					unsigned short groupNo = GetFragGroupNo_2(readLine[0], readLine[oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine, oneThirdReadLength);
					readCode += 2; //switching to second split

					//this is second split
					groupNo = GetFragGroupNo_2(readLine[oneThirdReadLength], readLine[readLength-oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine + oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the first split of the next read

					//this is third split
					groupNo = GetFragGroupNo_2(readLine[2 * oneThirdReadLength], readLine[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine + 2 * oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the first split of the next read
				}	
				fclose(finRead);
			}
		}
		else if(mappingModeStr == "PAIRED")
		{
			char *readLine1, *readLine2;
		
			unsigned long long readCode = 0;
			string curFastqFileName1, curFastqFileName2;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName1 >> curFastqFileName2 >> lineCount;
				
				FILE* finRead1 = fopen(curFastqFileName1.c_str(), "rb");
				InitFirstBuffer(finRead1);		
				FILE* finRead2 = fopen(curFastqFileName2.c_str(), "rb");
				InitSecondBuffer(finRead2);		
				while(true)
				{
					readLine1 = GetNextReadStart_FirstBuffer(finRead1);
					readLine2 = GetNextReadStart_SecondBuffer(finRead2);

					if(readLine1 == NULL)
					{
						assert(readLine2 == NULL);
						break;	
					}

					unsigned short groupNo = GetFragGroupNo_2(readLine1[0], readLine1[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine1, readLength);
				
					readCode += 2; //switching to the second mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_2(readLine2[0], readLine2[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine2, readLength);

					readCode += 2; //swithing to the first mate of the next read //for direction jumps (they're all forward at this point)
				}
				fclose(finRead1);
				fclose(finRead2);
			}
		}
		else if(mappingModeStr == "SPLIT_PAIRED")
		{
			char *readLine1, *readLine2;
			unsigned long long readCode = 0;
			assert(readLength % 2 == 0);
			unsigned short halfReadLength = readLength / 2;

			string curFastqFileName1, curFastqFileName2;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName1 >> curFastqFileName2 >> lineCount;
				
				FILE* finRead1 = fopen(curFastqFileName1.c_str(), "rb");
				InitFirstBuffer(finRead1);		
				FILE* finRead2 = fopen(curFastqFileName2.c_str(), "rb");
				InitSecondBuffer(finRead2);		
				while(true)
				{
					readLine1 = GetNextReadStart_FirstBuffer(finRead1);
					readLine2 = GetNextReadStart_SecondBuffer(finRead2);
					
					if(readLine1 == NULL)
					{
						assert(readLine2 == NULL);
						break;	
					}
					
					unsigned short groupNo = GetFragGroupNo_2(readLine1[0], readLine1[halfReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine1, halfReadLength);
					readCode += 2; //switching to the second split of the first mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_2(readLine1[halfReadLength], readLine1[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine1 + halfReadLength, halfReadLength);
					readCode += 2; //switching to the first split of the second mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_2(readLine2[0], readLine2[halfReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine2, halfReadLength);
					readCode += 2; //swithing to the second split of the second mate //for direction jumps (they're all forward at this point)
					
					groupNo = GetFragGroupNo_2(readLine2[halfReadLength], readLine2[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine2 + halfReadLength, halfReadLength);
					readCode +=2 ; //switching to the first split of the first mate of the next read
				}

				fclose(finRead1);
				fclose(finRead2);
			}
		}
		else if(mappingModeStr == "THREEWAY_PAIRED")
		{
			char *readLine1, *readLine2;
			unsigned long long readCode = 0;
			assert(readLength % 3 == 0);
			unsigned short oneThirdReadLength = readLength / 3;

			string curFastqFileName1, curFastqFileName2;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName1 >> curFastqFileName2 >> lineCount;
				
				FILE* finRead1 = fopen(curFastqFileName1.c_str(), "rb");
				InitFirstBuffer(finRead1);		
				FILE* finRead2 = fopen(curFastqFileName2.c_str(), "rb");
				InitSecondBuffer(finRead2);		
				while(true)
				{
					readLine1 = GetNextReadStart_FirstBuffer(finRead1);
					readLine2 = GetNextReadStart_SecondBuffer(finRead2);
					
					if(readLine1 == NULL)
					{
						assert(readLine2 == NULL);
						break;	
					}
					
					unsigned short groupNo = GetFragGroupNo_2(readLine1[0], readLine1[oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine1, oneThirdReadLength);
					readCode += 2; //switching to the second split of the first mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_2(readLine1[oneThirdReadLength], readLine1[readLength-oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine1 + oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the third split of the first mate //for direction jumps (they're all forward at this point)
					
					groupNo = GetFragGroupNo_2(readLine1[2*oneThirdReadLength], readLine1[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine1 + oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the first split of the second mate //for direction jumps (they're all forward at this point)


					groupNo = GetFragGroupNo_2(readLine2[0], readLine2[oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine2, oneThirdReadLength);
					readCode += 2; //swithing to the second split of the second mate //for direction jumps (they're all forward at this point)
					
					groupNo = GetFragGroupNo_2(readLine2[oneThirdReadLength], readLine2[readLength-oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine2 + oneThirdReadLength, oneThirdReadLength);
					readCode +=2 ; //switching to the third split of the second mate 
					
					groupNo = GetFragGroupNo_2(readLine2[2 * oneThirdReadLength], readLine2[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine2 + 2 * oneThirdReadLength, oneThirdReadLength);
					readCode +=2 ; //switching to the first split of the first mate of the next read
				}

				fclose(finRead1);
				fclose(finRead2);
			}
		}



		else
		{
			cout << "ERROR: Mapping mode undefined... " << endl;
			exit(5);
		}
	}
	else if(fragSignalSize == 4)
	{
		if(mappingModeStr == "SINGLE")
		{
			char* readLine;
			unsigned long long readCode = 0;
			string curFastqFileName;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName >> lineCount;

				FILE* finRead = fopen(curFastqFileName.c_str(), "rb");
				InitFirstBuffer(finRead);		
				while(true)
				{
					readLine = GetNextReadStart_FirstBuffer(finRead);
					
					if(readLine == NULL)
					{
						break;
					}
		
					unsigned short groupNo = GetFragGroupNo_4(readLine[0], readLine[1], readLine[readLength-2], readLine[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine, readLength);
					readCode += 2; //for direction jumps (they're all forward at this point)
				}	

				fclose(finRead);
			}
		}
		else if(mappingModeStr == "SPLIT_SINGLE")
		{
			char* readLine;
		
			unsigned long long readCode = 0;
		
			assert(readLength % 2 == 0);
			unsigned short halfReadLength = readLength / 2;

			string curFastqFileName;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName >> lineCount;
			
				FILE* finRead = fopen(curFastqFileName.c_str(), "rb");
				InitFirstBuffer(finRead);		
				while(true)
				{
					readLine = GetNextReadStart_FirstBuffer(finRead);

					if(readLine == NULL)
					{
						break;
					}

					unsigned short groupNo = GetFragGroupNo_4(readLine[0], readLine[1], readLine[halfReadLength-2], readLine[halfReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine, halfReadLength);
					readCode += 2; //switching to second split

					//this is second split
					groupNo = GetFragGroupNo_4(readLine[halfReadLength], readLine[halfReadLength+1], readLine[readLength-2], readLine[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine + halfReadLength, halfReadLength);
					readCode += 2; //switching to the first split of the next read
				}	
				fclose(finRead);
			}
		}
		else if(mappingModeStr == "THREEWAY_SINGLE")
		{
			char* readLine;
		
			unsigned long long readCode = 0;
		
			assert(readLength % 3 == 0);
			unsigned short oneThirdReadLength = readLength / 3;

			string curFastqFileName;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName >> lineCount;
			
				FILE* finRead = fopen(curFastqFileName.c_str(), "rb");
				InitFirstBuffer(finRead);		
				while(true)
				{
					readLine = GetNextReadStart_FirstBuffer(finRead);

					if(readLine == NULL)
					{
						break;
					}

					unsigned short groupNo = GetFragGroupNo_4(readLine[0], readLine[1], readLine[oneThirdReadLength-2], readLine[oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine, oneThirdReadLength);
					readCode += 2; //switching to second split

					//this is second split
					groupNo = GetFragGroupNo_4(readLine[oneThirdReadLength], readLine[oneThirdReadLength+1], readLine[readLength-oneThirdReadLength-2], readLine[readLength-oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine + oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the first split of the next read
					
					//this is third split
					groupNo = GetFragGroupNo_4(readLine[2*oneThirdReadLength], readLine[2*oneThirdReadLength+1], readLine[readLength-2], readLine[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine + 2 * oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the first split of the next read
				}	
				fclose(finRead);
			}


		}
		else if(mappingModeStr == "PAIRED")
		{
			char *readLine1, *readLine2;
		
			unsigned long long readCode = 0;
			string curFastqFileName1, curFastqFileName2;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName1 >> curFastqFileName2 >> lineCount;
				
				FILE* finRead1 = fopen(curFastqFileName1.c_str(), "rb");
				InitFirstBuffer(finRead1);		
				FILE* finRead2 = fopen(curFastqFileName2.c_str(), "rb");
				InitSecondBuffer(finRead2);		
				while(true)
				{
					readLine1 = GetNextReadStart_FirstBuffer(finRead1);
					readLine2 = GetNextReadStart_SecondBuffer(finRead2);

					if(readLine1 == NULL)
					{
						assert(readLine2 == NULL);
						break;	
					}
					
					unsigned short groupNo = GetFragGroupNo_4(readLine1[0], readLine1[1], readLine1[readLength-2], readLine1[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine1, readLength);
					readCode += 2; //switching to the second mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_4(readLine2[0], readLine2[1], readLine2[readLength-2], readLine2[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine2, readLength);
					readCode += 2; //swithing to the first mate of the next read //for direction jumps (they're all forward at this point)
				}
				fclose(finRead1);
				fclose(finRead2);
			}
		}
		else if(mappingModeStr == "SPLIT_PAIRED")
		{
			char *readLine1, *readLine2;
			unsigned long long readCode = 0;
			assert(readLength % 2 == 0);
			unsigned short halfReadLength = readLength / 2;
			
			string curFastqFileName1, curFastqFileName2;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName1 >> curFastqFileName2 >> lineCount;
				
				FILE* finRead1 = fopen(curFastqFileName1.c_str(), "rb");
				InitFirstBuffer(finRead1);		
				FILE* finRead2 = fopen(curFastqFileName2.c_str(), "rb");
				InitSecondBuffer(finRead2);		
				while(true)
				{
					readLine1 = GetNextReadStart_FirstBuffer(finRead1);
					readLine2 = GetNextReadStart_SecondBuffer(finRead2);
					
					if(readLine1 == NULL)
					{
						assert(readLine2 == NULL);
						break;	
					}
				
					unsigned short groupNo = GetFragGroupNo_4(readLine1[0], readLine1[1], readLine1[halfReadLength-2], readLine1[halfReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine1, halfReadLength);
					readCode += 2; //switching to the second split of the first mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_4(readLine1[halfReadLength], readLine1[halfReadLength+1], readLine1[readLength-2], readLine1[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine1 + halfReadLength, halfReadLength);
					readCode += 2; //switching to the first split of the second mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_4(readLine2[0], readLine2[1], readLine2[halfReadLength-2], readLine2[halfReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine2, halfReadLength);
					readCode += 2; //swithing to the second split of the second mate //for direction jumps (they're all forward at this point)
					
					groupNo = GetFragGroupNo_4(readLine2[halfReadLength], readLine2[halfReadLength+1], readLine2[readLength-2], readLine2[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine2 + halfReadLength, halfReadLength);
					readCode +=2 ; //switching to the first split of the first mate of the next read
				}

				fclose(finRead1);
				fclose(finRead2);
			}
		}
		else if(mappingModeStr == "THREEWAY_PAIRED")
		{
			char *readLine1, *readLine2;
			unsigned long long readCode = 0;
			assert(readLength % 3 == 0);
			unsigned short oneThirdReadLength = readLength / 3;
			
			string curFastqFileName1, curFastqFileName2;
			for(int k=0; k<numSamples; k++)
			{
				finFastqList >> curFastqFileName1 >> curFastqFileName2 >> lineCount;
				
				FILE* finRead1 = fopen(curFastqFileName1.c_str(), "rb");
				InitFirstBuffer(finRead1);		
				FILE* finRead2 = fopen(curFastqFileName2.c_str(), "rb");
				InitSecondBuffer(finRead2);		
				while(true)
				{
					readLine1 = GetNextReadStart_FirstBuffer(finRead1);
					readLine2 = GetNextReadStart_SecondBuffer(finRead2);
					
					if(readLine1 == NULL)
					{
						assert(readLine2 == NULL);
						break;	
					}
				
					unsigned short groupNo = GetFragGroupNo_4(readLine1[0], readLine1[1], readLine1[oneThirdReadLength-2], readLine1[oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine1, oneThirdReadLength);
					readCode += 2; //switching to the second split of the first mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_4(readLine1[oneThirdReadLength], readLine1[oneThirdReadLength+1], readLine1[readLength-oneThirdReadLength-2], readLine1[readLength-oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine1 + oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the third split of the first mate //for direction jumps (they're all forward at this point)

					groupNo = GetFragGroupNo_4(readLine1[2*oneThirdReadLength], readLine1[2*oneThirdReadLength+1], readLine1[readLength-2], readLine1[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine1 + 2 * oneThirdReadLength, oneThirdReadLength);
					readCode += 2; //switching to the first split of the second mate //for direction jumps (they're all forward at this point)


					groupNo = GetFragGroupNo_4(readLine2[0], readLine2[1], readLine2[oneThirdReadLength-2], readLine2[oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine2, oneThirdReadLength);
					readCode += 2; //swithing to the second split of the second mate //for direction jumps (they're all forward at this point)
					
					groupNo = GetFragGroupNo_4(readLine2[oneThirdReadLength], readLine2[oneThirdReadLength+1], readLine2[readLength-oneThirdReadLength-2], readLine2[readLength-oneThirdReadLength-1]);
					AddToBuffer(groupNo, readCode, readLine2 + oneThirdReadLength, oneThirdReadLength);
					readCode +=2 ; //switching to the third split of the second mate
					
					groupNo = GetFragGroupNo_4(readLine2[2*oneThirdReadLength], readLine2[2*oneThirdReadLength+1], readLine2[readLength-2], readLine2[readLength-1]);
					AddToBuffer(groupNo, readCode, readLine2 + 2 * oneThirdReadLength, oneThirdReadLength);
					readCode +=2 ; //switching to the first split of the first mate of the next read
				}

				fclose(finRead1);
				fclose(finRead2);
			}
		}



		else
		{
			cout << "ERROR: Mapping mode undefined... " << endl;
			exit(5);

		}
	}
	else
	{
		cout << "unexpected frag signal size: " << fragSignalSize << endl;
		exit(9);
	}

	for(unsigned int i=0; i<numGroups; i++)
	{
		FlushBuffer(i);
	}

	return 0;
}
