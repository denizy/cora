/*Copyright (c) 2015-2016 Deniz Yorukoglu. All rights reserved.*/
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<assert.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>
#include<stdio.h>
using namespace std;

const int FA_LINE_LEN = 50;
const int MAX_NUM_LONG_READ_NAMES = 2000000;
int curLongReadNameArrLen = MAX_NUM_LONG_READ_NAMES;
const int MAX_READ_NAME_LENGTH = 5000000;
int curReadNameArrLen = MAX_READ_NAME_LENGTH;
const int MAX_NUM_CHRS = 250;

#define IDENTITY_ALPHABET_START 33
#define IDENTITY_ALPHABET_END   126
#define IDENTITY_ALPHABET_SIZE  94

int repMapMode;

#define BWA_MODE 101
#define BOWTIE_MODE 102
#define MRSFAST_MODE 103
#define BWA_MEM_MODE 201
#define BOWTIE_2_MODE 202
#define MRSFAST_ULTRA_MODE 203
#define MANUAL_MODE 301

#define MAX_LINE_LEN 200
#define MAX_NUM_CHRS 250
#define MAX_CHR_NAME_LEN 50

#define MAX_NUM_EDITS 2


double getTime()
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return t.tv_sec+t.tv_usec/1000000.0;
}

string ConvertRefEditsToReadEdits(const string& seq, string refEdits);
void PrintCompactLinkFromSamLine_multiChr(const string& genMapLine, ofstream& foutLinks, int  numChr, string chrList[]);

string** longReadNameList; //[MAX_NUM_LONG_READ_NAMES+2];
int numLongReadNames;

unsigned char idDigitLen;

//Modified from Luke Chmella's itoa function implementation
unsigned short itoa10(int value, char* result)
{
        unsigned short len = 0;
        char* ptr = result, *ptr1 = result, tmp_char;
        int tmp_value;

        do
        {
                tmp_value = value;
                value /= 10;
                *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * 10)];
                len++;
        } while ( value );

        // Apply negative sign
        if (tmp_value < 0) *ptr++ = '-';
        *ptr-- = '\0';
        while(ptr1 < ptr)
        {
                tmp_char = *ptr;
                *ptr--= *ptr1;
                *ptr1++ = tmp_char;
        }

        return len;
}

//Data variables for reference
char* fullRef[MAX_NUM_CHRS+2];
unsigned int chrLens[MAX_NUM_CHRS+2];
char chrNames[MAX_NUM_CHRS+2][MAX_CHR_NAME_LEN+2];
char chrCodeNameList[MAX_NUM_CHRS+2][10];
int numChrs;

void InitRef(string refName, string refFaiName, int refLineLen)
{
	//Verifying refFileFai & reference size
        ifstream finFai(refFaiName.c_str());
        if(!finFai.is_open())
        {
                cout << "Could not find fai file for reference file: " << refName << endl;
                exit(16);
        }

        //In MultiChr Version only MaxChrSize is needed
        int chrCode = 0;
        string faiLine;

        while(getline(finFai, faiLine))
        {
		cout << "faiLine: " << faiLine << endl;

                stringstream faiLineSS(faiLine);
                string ctgName;
                int ctgSize;
                faiLineSS >> ctgName >> ctgSize;

                chrCode++;
                fullRef[chrCode] = (char *) calloc (ctgSize+MAX_LINE_LEN, sizeof(char));
                chrLens[chrCode] = ctgSize;
                strcpy(chrNames[chrCode], ctgName.c_str());
       
		stringstream ss;
                ss << chrCode;
                string sss;
                ss >> sss;

		for(unsigned int i=0; i<sss.length(); i++)
                {
                        chrCodeNameList[chrCode][i] = sss[i];
                }
                chrCodeNameList[chrCode][sss.length()] = '\0';
	}
        numChrs = chrCode;

	if(repMapMode == BWA_MEM_MODE)
	{
		//Dummy space allocation for the 0th chromosome (will not be used)
		fullRef[0] = (char *) malloc (500);

		finFai.close();

		FILE* finRef = fopen(refName.c_str(),"r");

		cout << "Reading full reference..." << endl;

		int curChrCode = 0;
		int curChrPos = 1;

		char* curChrPtr = fullRef[0];

		char* ret;
		do
		{
			ret = fgets(curChrPtr + curChrPos, MAX_LINE_LEN, finRef);
			if(curChrPtr[curChrPos] == '>')
			{
				curChrCode++;
				curChrPtr = fullRef[curChrCode];
				curChrPos = 1;
			}
			else
			{
				curChrPos += refLineLen;
			}
		}while(ret!=NULL);

		cout << "Renaming lowercases..." << endl;
		for(int chr=1; chr<=numChrs; chr++)
		{
			char* curPtr = fullRef[chr];
			curPtr[chrLens[chr]] = '\0';    // Puts end of string markers for each chromosome, for easier while-loop search
			int pos = 1;
			while(curPtr[pos] != '\0')
			{
				if(curPtr[pos] > 'Z')
					curPtr[pos] -= ('a'-'A');
				pos++;
			}
		}
		cout << "Finished loading reference." << endl;
		fclose(finRef);
	}
}

inline void IncrementReadId(char arr[], unsigned int rightEndPos)
{
        unsigned int leftEndPos = rightEndPos - idDigitLen + 1;

	arr[rightEndPos]++;
        while(arr[rightEndPos] > IDENTITY_ALPHABET_END)
        {
                arr[rightEndPos] = IDENTITY_ALPHABET_START;
                rightEndPos--;
                if(rightEndPos < leftEndPos)
                        break;
                arr[rightEndPos]++;
        }
}

inline void DecrementReadId(char arr[], unsigned int rightEndPos)
{
	unsigned int leftEndPos = rightEndPos - idDigitLen + 1;


        arr[rightEndPos]--;
        while(arr[rightEndPos] < IDENTITY_ALPHABET_START)
        {
                arr[rightEndPos] = IDENTITY_ALPHABET_END;
                rightEndPos--;
                if(rightEndPos < leftEndPos)
                        break;
                arr[rightEndPos]--;
        }
}

inline int readStr(char*& str, char arr[], char del)
{
	int pos = 0;
	while(str[pos] != del)
	{
		arr[pos] = str[pos];
		pos++;
	}
	arr[pos] = '\0';
	str += pos + 1;	
	return pos;
}


inline int readStr_name(char*& str, char*& arr, char del)
{
	//This also enables resizing of the name arr
	int pos = 0;
	while(str[pos] != del)
	{
		arr[pos] = str[pos];
		pos++;

		if(pos >= curReadNameArrLen)
		{
			int newCurReadNameArrLen = 2 * curReadNameArrLen;

			char* newArr = (char *) calloc (newCurReadNameArrLen + 2, sizeof(char)); 

			for(int i=0; i<curReadNameArrLen; i++)
			{
				newArr[i] = arr[i];	
			}

			free(arr);
			
			arr = newArr;	
			
			curReadNameArrLen = newCurReadNameArrLen;
			
			cout << "curReadNameArrLen resized to: " << curReadNameArrLen << endl;

		}
	}
	arr[pos] = '\0';
	str += pos + 1;	
	return pos;
}


void findFieldWithPrefix(char*& str, char retArr[], unsigned short prefixLen) // retArr contains the prefix,  reads the remaining str and look for the field that starts with prefix (until endl) otherwise put '\0' at beginning of retArr
{
	while(true)
	{
		bool correctPrefix = 1;
		for(int i=0; i<prefixLen; i++)
		{
			if(str[i] != retArr[i])
			{
				correctPrefix = 0;
				break;	
			}
		}

		if(correctPrefix)
		{
			str += prefixLen;
			unsigned short retArrInd = prefixLen;
			while(str[0] > 10)
			{
				retArr[retArrInd++] = str[0];
				str++;
			} 
			retArr[retArrInd] = '\0';
			return;
		}
		else
		{
			while(str[0] > 10) //9 is TAB, 10 is newline
			{
				str++;
			}
			switch(str[0])
			{
				case 9:
					str++;
					break; //just skip to the next field
				case 10:
					return;
				
				default:
					retArr[0] = '\0';
					return;
			}
		}
	}
}

inline void skipNtabs(char*& str, int skipCount)
{
	int numTabs = 0;
	while(numTabs < skipCount)
	{
		if(str[0] == '\t')
		{
			numTabs++;
		}
		str++;
	}
}


void assignEditString(char out[], char chrNo[], char chrPosStr[], char dir[], char seq[], unsigned short readLen)
{
	int numEdits = 0;
	char* ret = out;
	for(int i=1; i<=numChrs; i++)
	{
		if(strcmp(chrNames[i], chrNo) == 0)
		{
			unsigned int chrPos = atoi(chrPosStr);
			unsigned short prevPos = 0;
			for(unsigned short k=1; k<=readLen; k++)
			{
				if(fullRef[i][chrPos + k - 1] != seq[k - 1])
				{
					numEdits++;
					assert(numEdits <= MAX_NUM_EDITS);					

					if(k - prevPos > 1)
					{
						unsigned short len = itoa10(k - prevPos - 1, ret);
						ret += len;
					}
					ret[0] = seq[k-1]; //edit vals will come from seq not ref		 
					ret++;
					prevPos = k;
				}
			}

			if(readLen - prevPos >= 1)
			{
				unsigned short len = itoa10(readLen - prevPos, ret);
                                ret += len;
			}

			ret[0] = '\0';
			return;	
		}
	}
	assert(0);
}

int AddIntegerToString(char* str, int num) //a limited version of atoi for numbers < 999
{
	assert(num < 1000 && num > 0);

	if(num >= 100)
	{
		str[0] = (num / 100) + '0';
		str[1] = ((num % 100) / 10) + '0';
		str[2] = (num % 10) + '0';
		return 3;
	}
	else if(num >= 10)
	{
		str[0] = (num / 10) + '0';
		str[1] = (num % 10) + '0';
		return 2;
	}
	else
	{
		str[0] = num + '0';
		return 1;
	}
}

int main(int argc, char* argv[])
{
	double beginTime = getTime();
	if(argc!=13 || string(argv[7]).length() != 1)
	{
		cout << "This program prints a list of links, given a SAM output file printed by an off-the-shelf mapper" << endl;
		cout << "ARGV[1 - in] genome reference input" << endl;
		cout << "ARGV[2 - in] set of mappings to genome ref" << endl;
		cout << "ARGV[3 - in] inchworm assembly ref (or NOREF)" << endl; // In the future it might be sensible to split this as it gets larger
		cout << "ARGV[4 - in] mappings onto inchworm assembly ref (any string if [3] is NOREF)" << endl;
		cout << "ARGV[5 - in] long read names list for incorporating them" << endl;
		cout << "ARGV[6 - out] read link output - for both template A and B" << endl;
		cout << "ARGV[7 - param] number of errors (edit or hamming) per readMer -- currently only one digit" << endl; 
		cout << "ARGV[8 - param] idDigitLen" << endl;
		cout << "ARGV[9 - param] representative mapping mode" << endl;
		cout << "ARGV[10 - param] readLen" << endl;
		cout << "ARGV[11 - param] refLineLen" << endl;
		cout << "ARGV[12 - string param] temporary directory name" << endl;
		exit(131); 
	}
	string temporaryDirectoryName = argv[12];

	int numErrorsPerReadMer = atoi(argv[7]);

	//representative mapping = coarse mapping
	//Set up off-the-shelf mapper mode
	string representativeMappingModeStr(argv[9]);
	if(representativeMappingModeStr == "BWA")
	{
		repMapMode = BWA_MODE;
	}
	else if(representativeMappingModeStr == "BOWTIE")
	{
		repMapMode = BOWTIE_MODE;
	}
	else if(representativeMappingModeStr == "MRSFAST")
	{
		repMapMode = MRSFAST_MODE;
	}
	else if(representativeMappingModeStr == "BWA_MEM")
	{
		repMapMode = BWA_MEM_MODE;

	}
	else if(representativeMappingModeStr == "BOWTIE_2")
	{
		repMapMode = BOWTIE_2_MODE;
	}
	else if(representativeMappingModeStr == "MRSFAST_ULTRA")
	{
		repMapMode = MRSFAST_ULTRA_MODE;
	}
	else if(representativeMappingModeStr == "MANUAL")
	{
		repMapMode = MANUAL_MODE;
	}
	else
	{
		cout << "Representative Mapping Mode: " << repMapMode << " is not supported" << endl;
		exit(89);
	}

	int readLen = atoi(argv[10]);
	int refLineLen = atoi(argv[11]);

	idDigitLen = (unsigned char) atoi(argv[8]);

	//Initialize chromosome list for the merged ref
	string genRef(argv[1]);
	string genRefFai = genRef + ".fai";
	InitRef(genRef, genRefFai, refLineLen); //Read fai -- reference is also read for BWA-mem (which doesn't print edit strings).

	double faiTime = getTime();
	cout << "FaiTime: " << faiTime - beginTime << endl;

	string iwRef(argv[3]);
	assert(iwRef == "NOREF"); //This part is deactivated for compressive mapping
	
	double refTime = getTime();
	cout << "RefTime: " << refTime - faiTime << endl;
	
	//Set up long read name library
	ifstream finLongReadNames(argv[5]);
	string junkL, longReadName;

	string** longReadNameList = (string**) calloc (MAX_NUM_LONG_READ_NAMES+2, sizeof(string*));

	while(finLongReadNames >> junkL >> longReadName)
	{
		numLongReadNames++;
	
		if(numLongReadNames >= curLongReadNameArrLen)
		{
			int newCurLongReadNameArrLen = 2 * curLongReadNameArrLen;

			string** newLongReadNameList = (string**) calloc (newCurLongReadNameArrLen + 2, sizeof(string*)); 

			for(int i=0; i<curLongReadNameArrLen; i++)
			{
				newLongReadNameList[i] = longReadNameList[i];	
			}

			delete [] longReadNameList;
			
			longReadNameList = newLongReadNameList;	
			
			curLongReadNameArrLen = newCurLongReadNameArrLen;	
			cout << "curLongReadNameArrLen resized to: " << curLongReadNameArrLen << endl;
		}

		longReadNameList[numLongReadNames] = new string(longReadName);
	}

	FILE* finGenMap = fopen(argv[2], "rb");
	
	setvbuf(finGenMap, NULL, _IOFBF, 1048576);

	if(finGenMap == NULL)
	{
		cout << "ERROR: mapping file couldn't be opened: '" << argv[2] << "'" << endl;
		exit(41);
	}

	FILE* foutLinks = fopen(argv[6], "wb");	
	setvbuf(foutLinks, NULL, _IOFBF, 1048576);

	//Generalized header-line skipping for all of the mappers
	string headerCountFile = temporaryDirectoryName + "/__tempHeaderLineCount";
	stringstream headerLineCountCall;
	headerLineCountCall << "head -" << max(100, 2 * (int) numChrs) << " " << argv[2] << " | grep ^@ | wc -l > " << headerCountFile;
	int x = system(headerLineCountCall.str().c_str());
	cout << "headerLineCountCall: " << headerLineCountCall.str() << "\t\texitCode: " << x <<  endl;		
	ifstream finTempHeaderLineCount(headerCountFile.c_str());

	int numHeaderLines;
	finTempHeaderLineCount >> numHeaderLines;
	char* linebuf=0;
	size_t linesiz=0;
	for(int k=0; k<numHeaderLines; k++)
	{
		assert(getline(&linebuf, &linesiz, finGenMap));
		free(linebuf);
		linebuf=NULL;
	}	
	/////////////////////////////////////////////////////////
	

	char* readName = (char*) calloc (MAX_READ_NAME_LENGTH+2, sizeof(char));
	char cigarStr[32];
	char ctgName[32];
	char mapDirStr[8];
	char ctgPosStr[16];
	char seq[256];
	char edits[32];
	char nm[20]; //for BWA mem
	unsigned int filledChars = 0;
        unsigned int bufferSize = 5005000; //This should be larger than the longest link name length (LRN filtering should be lower than this threshold )choosing 5m for LRN
        char *buffer = (char *) malloc (bufferSize + 5);
        buffer[bufferSize] = '\0';

	//[TODO] Make this parametric
	unsigned int minSkip = readLen + idDigitLen + 21 - 15; //This is the minimum distance to skip ahead to find the '\n'
	int readInsertionList[10]; //This holds the insertion positions for a read, the count of them is calculated within the loop for each read
	
	while(!feof(finGenMap)) //A lot of buffer computations below to quickly go through the sam file and print links
	{
		size_t result = fread(buffer+filledChars, 1, bufferSize-filledChars, finGenMap);	

		if(result != bufferSize-filledChars)
                {
                        bufferSize = filledChars + result;
           	}
		buffer[bufferSize] = '\n';

                unsigned int bufferOffset = 0;
		filledChars = 0;

                unsigned int offsetForNextLine = 0;
                while(buffer[bufferOffset] != '\0')
                {
                       //Find sentence size
                        bool fullyContainedFlag = 0;

                        unsigned int searchStartPoint = bufferOffset;
                        searchStartPoint += minSkip;

                        for(unsigned int i=searchStartPoint; i<bufferSize; i++)
                        {
                                if(buffer[i] == '\n')
                                {
                                        fullyContainedFlag = 1;
                                        offsetForNextLine = i+1;
                                        break;
                                }
                        }

                        if(!fullyContainedFlag)
                        {
                                for(unsigned int i=bufferOffset; i<bufferSize; i++)
                                {
                                        buffer[i-bufferOffset] = buffer[i];
                                }
                                filledChars = bufferSize - bufferOffset;
                                break;
                        }

			char* str = buffer + bufferOffset;
			readStr_name(str, readName, '\t');
			readStr(str, mapDirStr, '\t');
			readStr(str, ctgName, '\t');

			if(ctgName[0] == '*')
			{
				bufferOffset = offsetForNextLine;
				continue;
			}	
			
			readStr(str, ctgPosStr, '\t');

			skipNtabs(str, 1);

			//Read Cigar here
			readStr(str, cigarStr, '\t');

			int numInsertionsToRead = 0;

			int totalReadOffset = 0; //this holds the read offset for the current flag;
			int lengthForFlag = 0; //this holds the length of any 'M' 'I' 'D' flag.
			for(int i=0; cigarStr[i] != '\0'; i++)
			{
				if(cigarStr[i] >= '0' && cigarStr[i] <= '9')
				{
					lengthForFlag = lengthForFlag * 10 + cigarStr[i] - '0';	
				}
				else if(cigarStr[i] == 'M')
				{
					totalReadOffset += lengthForFlag;
					lengthForFlag = 0;
				}
				else if(cigarStr[i] == 'D')
				{
					lengthForFlag = 0;
				}
				else if(cigarStr[i] == 'I' || cigarStr[i] == 'S')
				{
					//there is an insertion to the read
					for(int k=1; k<=lengthForFlag; k++)
					{
						readInsertionList[numInsertionsToRead] = totalReadOffset + k;
						numInsertionsToRead++;		
					}
					totalReadOffset += lengthForFlag;
					lengthForFlag = 0;
				}
				else
				{
					cout << "ERROR: Currently on coarse mapping with CIGAR strings with M I D S flags supported" << endl;
					cout << "Current flag: " << cigarStr[i] << endl;
					exit(99);
				}
			}

			skipNtabs(str, 3);

			readStr(str, seq, '\t');


			if(repMapMode == MRSFAST_ULTRA_MODE || repMapMode == BOWTIE_MODE || repMapMode == MRSFAST_MODE)	
			{		
				skipNtabs(str, 2);
				readStr(str, edits, '\n');
			}
			else if(repMapMode == BOWTIE_2_MODE|| repMapMode == MANUAL_MODE)
			{
				edits[0] = 'M';
				edits[1] = 'D';
				edits[2] = ':';
				edits[3] = 'Z';
				edits[4] = ':';

				findFieldWithPrefix(str, edits, 5);
				
				//any field after MD:Z will not be used, automatically skipped to the newline at the end of while loop
				assert(edits[0] != '\0');
			}
			else if(repMapMode == BWA_MODE)
			{
				skipNtabs(str, 2); //these are "*" for quality, and the XT field
				
				int len = readStr(str, nm, '\t');
				if(nm[len-1] > numErrorsPerReadMer + '0')
				{
					bufferOffset = offsetForNextLine;
					continue;
				}

				edits[0] = 'M';
				edits[1] = 'D';
				edits[2] = ':';
				edits[3] = 'Z';
				edits[4] = ':';

				findFieldWithPrefix(str, edits, 5);
				
				//any field after MD:Z will not be used, automatically skipped to the newline at the end of while loop
				assert(edits[0] != '\0');

			}
			else if(repMapMode == BWA_MEM_MODE)
			{
				skipNtabs(str, 1);
				int len = readStr(str, nm, '\t');

				if(nm[len-1] > numErrorsPerReadMer + '0' )
 				{
					bufferOffset = offsetForNextLine;
                                	continue;
				}

				edits[0] = 'M';
                                edits[1] = 'D';
                                edits[2] = ':';
                                edits[3] = 'Z';
                                edits[4] = ':';

				//[TODO] this function might be outdated, check if that's true
				assignEditString(edits + 5, ctgName, ctgPosStr, mapDirStr, seq, readLen);
			}
			else
			{
				cout << "Field number for MD:Z not known... Add representative Mapping mode for: " << representativeMappingModeStr << endl;
                		exit(42);	
			}			


			int nameLen = strlen(readName);
			//long read name conversion as needed
			if(readName[0] == '#')
			{
				if(nameLen % idDigitLen != 0) //LRN's aren't divisible by idDigitLen (this is to prevent some fluke similarity between read identities and LRN signals)
				{
					int code = atoi(readName + 1); //atoi doesn't mind trailing nonnumeric characters
					int longLength = longReadNameList[code]->length(); 
				
					for(int i=0; i<longLength; i++)
					{
						
						readName[i] = longReadNameList[code]->at(i);
					}			
					readName[longLength] = '\0';
					nameLen = longLength;
				}
			}
			
			int mapDir = atoi(mapDirStr);
		
			if(mapDir % 32 >= 16) //means reverse direction alignment
			{
				for(unsigned int i=0; readName[i] != '\0'; i+=idDigitLen)
				{
					if(readName[i + idDigitLen - 1] % 2 == 1) //ONLY works if alphabet size is even
					{
						IncrementReadId(readName, i + idDigitLen - 1);
					}
					else
					{
						DecrementReadId(readName, i + idDigitLen - 1);
					}
				}
			}

			readName[nameLen] = '\t';
			nameLen++;
			
			for(int i=1; i<=numChrs; i++)
			{
				if(strcmp(ctgName, chrNames[i]) == 0)
				{
					for(int k=0; chrCodeNameList[i][k] != '\0'; k++)
					{
						readName[nameLen] = chrCodeNameList[i][k];
						nameLen++;
					}
					readName[nameLen] = '\t';
					nameLen++;
					break;
				}
			}
		
			for(int i=0; ctgPosStr[i]!='\0'; i++)
			{
				readName[nameLen] = ctgPosStr[i];
				nameLen++;
			}	

			readName[nameLen] = '\t';
			nameLen++;

			if(numInsertionsToRead == 0)
			{
				short offset = 0;
				short curOffset = 0;

				for(int i=5; edits[i] != '\0'; i++)
				{
					char ch = edits[i];

					if(ch >= '0' && ch <= '9')
					{
						curOffset = curOffset * 10 + ch - '0';
						readName[nameLen] = ch;
						nameLen++;
					}
					else
					{
						short totalOffset = offset + curOffset;

						if(ch == '^')
						{

							i++;
							totalOffset--; //need to call this only once for the deletion block
							while(edits[i] >= 'A')
							{
								readName[nameLen] = 'D';
								nameLen++;
								i++;

							}
							i--;
						}
						else
						{
							readName[nameLen] = seq[totalOffset];
							nameLen++;
						}

						offset = totalOffset + 1;
						curOffset = 0;
					}
				}
			}
			else //if(numInsertionsToRead > 0) -- Handle insertions simultaneously
			{
				int curInsertionIndex = 0; //this is the index of the insertion in the readInsertionList

				int numericReplacePosition = nameLen; //this is the index of the last position + 1 in the readName char array excluding an integers printed since

				short offset = 0;
				short curOffset = 0;

				for(int i=5; edits[i] != '\0'; i++)
				{
					char ch = edits[i];

					if(ch >= '0' && ch <= '9')
					{
						curOffset = curOffset * 10 + ch - '0';
						readName[nameLen] = ch;
						nameLen++;
					}
					else
					{
						if(curInsertionIndex < numInsertionsToRead && offset < readInsertionList[curInsertionIndex] && readInsertionList[curInsertionIndex] <= offset + curOffset + 1)
						{
							//Rewrite the integer to match the insertions
							int k = curInsertionIndex;
							for(; k<numInsertionsToRead && offset < readInsertionList[k] && readInsertionList[k] <= offset + curOffset + 1; k++)
							{
								if(k == curInsertionIndex)
								{
									if(readInsertionList[k] - offset - 1 > 0)
									{
										numericReplacePosition += AddIntegerToString(readName + numericReplacePosition, readInsertionList[k] - offset - 1);
									}
								}
								else if(readInsertionList[k] - readInsertionList[k-1] != 1) //if it's 0 no need to put a zero in between the insertions
								{
									numericReplacePosition += AddIntegerToString(readName + numericReplacePosition, readInsertionList[k] - readInsertionList[k-1] - 1);
								}
				
								readName[numericReplacePosition] = seq[readInsertionList[k]-1] - 'A' + 'a'; //lowercase letters represent insertions in the read
								numericReplacePosition++;
								curOffset++; //since insertion became part of the block
							}

							curInsertionIndex = k;
			
							if(offset + curOffset - readInsertionList[curInsertionIndex-1] > 0)
							{
								numericReplacePosition += AddIntegerToString(readName + numericReplacePosition, offset + curOffset - readInsertionList[curInsertionIndex-1]);
							}

							nameLen = numericReplacePosition;
						} 

						short totalOffset = offset + curOffset;

						if(ch == '^')
						{
							i++;
							totalOffset--; //need to call this only once for the insertion block
							while(edits[i] >= 'A')
							{
								readName[nameLen] = 'D';
								nameLen++;
								i++;

							}
							i--;
						}
						else
						{
							readName[nameLen] = seq[totalOffset];
							nameLen++;
						}

						numericReplacePosition = nameLen;

						offset = totalOffset + 1;
						curOffset = 0;
					}
				}

				if(curInsertionIndex < numInsertionsToRead && offset < readInsertionList[curInsertionIndex] && readInsertionList[curInsertionIndex] <= offset + curOffset + 1)
				{
					//Rewrite the integer to match the insertions
					int k = curInsertionIndex;
					for(; k<numInsertionsToRead && offset < readInsertionList[k] && readInsertionList[k] <= offset + curOffset + 1; k++)
					{
						if(k == curInsertionIndex)
						{
							if(readInsertionList[k] - offset - 1 > 0)
							{
								numericReplacePosition += AddIntegerToString(readName + numericReplacePosition, readInsertionList[k] - offset - 1);
							}
						}
						else if(readInsertionList[k] - readInsertionList[k-1] != 1) //if it's 0 no need to put a zero in between the insertions
						{
							numericReplacePosition += AddIntegerToString(readName + numericReplacePosition, readInsertionList[k] - readInsertionList[k-1] - 1);
						}
		
						readName[numericReplacePosition] = seq[readInsertionList[k]-1] - 'A' + 'a'; //lowercase letters represent insertions in the read
						numericReplacePosition++;
						curOffset++; //since insertion became part of the block
					}

					curInsertionIndex = k;
	
					if(offset + curOffset - readInsertionList[curInsertionIndex-1] > 0)
					{
						numericReplacePosition += AddIntegerToString(readName + numericReplacePosition, offset + curOffset - readInsertionList[curInsertionIndex-1]);
					}

					nameLen = numericReplacePosition;
				} 
			}
			
			readName[nameLen] = '\n';
			nameLen++;		

			fwrite(readName, sizeof(char), nameLen, foutLinks);  
			bufferOffset = offsetForNextLine;		
		}
	}
	
	double samTime = getTime();
	cout << "SamReadAndPrintTime: " << samTime - refTime << endl; 
	
	return 0;
}
