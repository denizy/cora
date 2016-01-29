/*Copyright (c) 2015-2016 Deniz Yorukoglu. All rights reserved.*/
#include<iostream>
#include<fstream>
#include<string>
#include<assert.h>
#include<sstream>
#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<algorithm>
#include<limits.h>
#include<time.h>
#include<sys/time.h>
#include<tr1/unordered_map>
#include<string.h> //[TODO] see if any library below this is used
#include<iomanip>  
#include<bitset>
#include<unistd.h>
#include<ios>
#include<omp.h>
using namespace std;

#define MAX_NUM_CHR 250 //[TODO] Increase the size of this to 65536
#define MAX_LINE_LEN 250 //[TODO] Increase the size of this to 65536
#define MAX_READ_LEN 250

#define NUM_MAX_SEEDS 4

#define DEBUG_SWITCH_FOR_EXACT_PRINT 0 //If this is turned on, the equivalence class will be printed in text
#define DEBUG_SWITCH_FOR_INEXACT_PRINT 0 //If this is turned on, all inexact homologies will be printed in text 
#define DEBUG_SWITCH_FOR_DETAILED_PRINT 0 //If this is turned on, every homology pair will be printed with all redundancy 

ofstream foutDetailedDebugPrint("DebugDetailedOutput");

#define BUFFER_ITEMS 100000
#define MAX_THREADS 30 //[TODO] Increase the size to 256

//################################################################//

char RevCompChar[256];

void SetupRevCompChar()
{
        RevCompChar['A'] = 'T';
        RevCompChar['T'] = 'A';
        RevCompChar['C'] = 'G';
        RevCompChar['G'] = 'C';
        RevCompChar['N'] = 'N';

//for insertions
        RevCompChar['a'] = 't';
        RevCompChar['t'] = 'a';
        RevCompChar['c'] = 'g';
        RevCompChar['g'] = 'c';
        RevCompChar['n'] = 'n';

//for deletion

        RevCompChar['D'] = 'D';
}

//Reference genome data
int chrLens[MAX_NUM_CHR+2]; //length of each chromosome
string* chrNames[MAX_NUM_CHR+2];
int numChrs; //number of chromosomes

int readLength; //global read-mer length (this is full read length for normal collapsing, but half of the read length for split-collapsing (homology table should be constructed for the read -mer length accordingly)

char* fullRef[MAX_NUM_CHR+2];  //containes nucleotides in the reference genome (1-based)
char* fullRefMarks[MAX_NUM_CHR+2]; //Marks are encoded in this array (the number shows the signal id for easy identification later on)
int refLineLen; //Length of each line in the reference (excluding last line in each chromosome)

double getTime()
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return t.tv_sec+t.tv_usec/1000000.0;
}

void ReadReference(const string& refFileName) //load all nucleotides to memory
{
	FILE* finRef = fopen(refFileName.c_str(),"r");

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
		curPtr[chrLens[chr]+1] = '\0';    // Puts end of string markers for each chromosome, for easier while-loop search
		int pos = 1;
		while(curPtr[pos] != '\0')
		{
			if(curPtr[pos] > 'Z')
				curPtr[pos] -= ('a'-'A');
			pos++;
		}
	}

	cout << "Finished loading reference." << endl;
}

void ConfigureRelatedParameters(const string& refFile)
{
	//Verifying refFile & reference line length
        ifstream finRef(refFile.c_str());
        if(!finRef.is_open())
        {
                cout << "ERROR: Could not find reference file: " << refFile << endl;
                cout << "Exiting.." << endl;
                exit(0);
        }

        string refLine, seqLine;
        finRef >> refLine >> seqLine;
        if(refLine[0] != '>')
        {
                cout << "ERROR: Reference file is not properly formatted..." << endl;
                cout << "Exiting.." << endl;
                exit(0);
        }

        refLineLen = seqLine.length(); //assume all other lines would be the same

        cout << "Reference line length: " << refLineLen << "bp .." << endl;
        finRef.clear();
        finRef.close();

        //Naming fai file
        string refFileFai = refFile + ".fai";

        //Verifying refFileFai & reference size
        ifstream finFai(refFileFai.c_str());
        if(!finFai.is_open())
        {
                cout << "Could not find fai file for reference file: " << refFile << endl;
                exit(0);
        }

       	int chrCode = 0;
	string faiLine;
        while(getline(finFai, faiLine)) //Set up chromosome count and lengths
        {
                stringstream faiLineSS(faiLine);

                string ctgName;
                int ctgSize;
                faiLineSS >> ctgName >> ctgSize;

		chrCode++;
		fullRef[chrCode] = (char *) calloc (ctgSize+MAX_LINE_LEN, sizeof(char));
		fullRefMarks[chrCode] = (char *) calloc (ctgSize+MAX_LINE_LEN, sizeof(char));
        	chrLens[chrCode] = ctgSize;
		chrNames[chrCode] = new string(ctgName);
	}
	numChrs = chrCode;

	//Dummy space allocation for the 0th chromosome (will not be used)
	fullRef[0] = (char *) calloc (500, sizeof(char));
        
	finFai.close();	
}

//###### HOM TABLE EXACT COLLAPSING #######//
struct eqItem
{
	char dir; //dir = 0 from representative, for other items, it is the direction to the repesentative
	char chrCode; //chromosome id
	unsigned int chrPos; //chromosome position
};

struct eqClass //each class has a representative (the first item eqList[0] and the remaining items represent exact homologies to the representative
{
	eqItem* eqList;
	int count; //number of items in the equivalence class together with the representative of the class
};

void DebugPrintEqItem(eqItem x)
{
	cout << "EqItem: chrCode: " << x.chrCode << "  chrPos: " << x.chrPos << "  dir: " << x.dir << endl;
}

bool sort_eqClass_by_representative(eqClass e1, eqClass e2)
{
	assert(e2.eqList != NULL && e1.eqList != NULL);
	return e1.eqList[0].chrCode == e2.eqList[0].chrCode ? (e1.eqList[0].chrPos < e2.eqList[0].chrPos) : e1.eqList[0].chrCode < e2.eqList[0].chrCode;
}

void ClearEquivClass(vector<eqClass>& equivSet)
{
	for(unsigned long long i=0; i<equivSet.size(); i++)
	{
		if(equivSet[i].count != -1)
		{
			free(equivSet[i].eqList);
		}
	}	
	equivSet.clear();
}

//Read the equiv classes 
void FillInEquivClasses(const string& eqInputFile, vector<eqClass>& equivSet)
{
// The input format is a bit tricky: the input will be split in 10 signal related chapters
//	Each chapter has different sections for different chromosomes
//	E[num1] [num2] [num3]
//	num1 represents the equivalence class no (as offset - if doesn't exist, offset is 1)
//	num2 is the direction of the current readmer in the equiv class
//	num3 is the chromosome position (as offset - if doesn't exist, offset is 1) 

	ifstream finEq(eqInputFile.c_str());

	long long numberOfLinesRead = 0;
	long long maxEquivNoSoFar = 0;
	long long prevEquivClassNo = 0; 
	cout << "Pre-scan file to get class sizes" << endl;
	string eqLine;

	double lastTime = getTime();

	while(getline(finEq,eqLine))
        {
		numberOfLinesRead++;

		if(numberOfLinesRead % 10000000 == 0)
		{
			double curTime = getTime();
			cout << numberOfLinesRead / 1000000 << "M\t" << curTime - lastTime << endl;
                	lastTime = curTime;
		}
		stringstream eqLineSS(eqLine);

                string f1 = "",f2 = "",f3 = "";
                eqLineSS >> f1 >> f2 >> f3;
		
		if(f1 == "signals")
		{
			continue;
		}
			
		if(f1 == "chrCode")
		{
			prevEquivClassNo = 0;
			continue;
		}			

		long long eOffset = 1;
		if(f1[0] == 'E')
		{
			eOffset = atoll(f1.substr(1,f1.length()-1).c_str());
			f1 = f2; //shift the fields as if there were no eqfield
			f2 = f3;
		}

		long long tempEquivNoValue = prevEquivClassNo + eOffset;
		long long equivNo = tempEquivNoValue;

		assert(equivNo <= maxEquivNoSoFar + 1);
		assert(equivNo == tempEquivNoValue);

		if(equivNo > maxEquivNoSoFar)
		{
			eqClass es;
			es.eqList = NULL;
			es.count = 1;
			
			equivSet[equivNo-1] = es;
			maxEquivNoSoFar = equivNo;
		}
		else
		{
			equivSet[equivNo-1].count++;
		}

		prevEquivClassNo = equivNo;
	}

	cout << "Finished scanning for counts.. resetting file" << endl;
	
	finEq.clear();
	finEq.close();
	finEq.open(eqInputFile.c_str());

	//Allocate new array
	cout << "Allocating new array" << endl;
	long long numEquivClass = equivSet.size();
	cout << "Resetting count array" << endl;
	//Resetting here in order to start counting again.
	for(long long i=0; i<numEquivClass; i++)
	{
		if(equivSet[i].count > 1)
		{
			equivSet[i].eqList = (eqItem*) malloc (equivSet[i].count * sizeof(eqItem)); //This saves space by not allocating any space for single item equivalence classes
			equivSet[i].count = 0;
		}
		else
		{
			equivSet[i].count = -1; //This is the marker for single item class
		}
	}

	cout << "Scanning files again for filling in equivalence class info into preallocated array.." << endl;
	maxEquivNoSoFar = 0;
	prevEquivClassNo = 0; 
	long long prevChrPos = 0;
	char curChrCode = -1;

	//Scan and place in lists (load everything in a single scan)
	while(getline(finEq,eqLine))
	{
		stringstream eqLineSS(eqLine);
		string f1 = "",f2 = "",f3 = "";
		eqLineSS >> f1 >> f2 >> f3;

		if(f1 == "signals")
		{
			continue;
		}
			
		if(f1 == "chrCode")
		{
			prevEquivClassNo = 0;
			prevChrPos = 0;
			curChrCode = atoi(f2.c_str());		
			continue;
		}			

		long long  eOffset = 1;
		if(f1[0] == 'E')
		{
			eOffset = (long long) atol(f1.substr(1,f1.length()-1).c_str());
			f1 = f2; //shift the fields as if there were no eqfield
			f2 = f3;
		}

		long long pOffset = 1;
		if(f2 != "")
		{
			pOffset = atol(f2.c_str()); 
		}
		
		long long equivNo = prevEquivClassNo + eOffset;
		bool dir = (f1 == "1");
		unsigned int chrPos = prevChrPos + pOffset;

		int curCount = equivSet[equivNo-1].count;

		if(curCount != -1)
		{
			eqItem newEqItem;
			newEqItem.dir = (char) dir;
			newEqItem.chrPos = chrPos;
			
			newEqItem.chrCode = curChrCode; //this reversal is for finding out the end point (it doesn't affect order of equivArr since the first item in each list is guaranteed to have posiitve chrcode
			equivSet[equivNo-1].eqList[curCount] = newEqItem;			
			equivSet[equivNo-1].count++;
		}
	
		prevEquivClassNo = equivNo;
		prevChrPos = chrPos;
	}

	long long totalEquivCount = numEquivClass;

	//Go over and remap the equivalence class values (so that they are smaller and are ordered)
	cout << "Second pass over equivalence class for compacting.." << endl;

	long long finalI = 0; //the position to swap to during selection
	for(long long i=0; i<totalEquivCount; i++)
	{
		int posSize = equivSet[i].count;//equivClassCounts[i];
		if(posSize > 1)
		{
			if(i!=finalI)
			{
				equivSet[finalI]=equivSet[i];
				if(equivSet[finalI].eqList[0].dir == 1) //revert everyone's direction to make base forward
				{
					int curVecLen = equivSet[finalI].count;
					assert(curVecLen > 1);
					for(int j=0; j<curVecLen; j++)
					{
						equivSet[finalI].eqList[j].dir = !equivSet[finalI].eqList[j].dir; 
					}								
				}
			}
			finalI++;
		}
	}
	equivSet.resize(finalI);
	
	for(unsigned long long i=0; i<equivSet.size(); i++)
	{
		assert(equivSet[i].eqList != NULL);
	}


	
	//TODO(denizy) To make this efficient in the future, do a merge sort that is aware of the split intervals (since each split interval is sorted within itself)

	cout << "sorting equivList..." << endl;
	//Need to sort due to split runs
	sort(equivSet.begin(), equivSet.end(), sort_eqClass_by_representative);
	cout << "Equiv Class Loading Finished.." << endl;
}

//Identify whether different equivalence classes can be merged into one
unsigned int FindEquivBlockSize(long long start, vector<eqClass>& equivSet)
{
	unsigned int distance = 1;
	eqItem* baseArr = equivSet[start].eqList;
	int baseArrSize = equivSet[start].count;	

	long long compactEquivClassCount = equivSet.size();
	while(start+distance < compactEquivClassCount && baseArr[0].chrCode == equivSet[start+distance].eqList[0].chrCode && baseArr[0].chrPos + distance == equivSet[start+distance].eqList[0].chrPos)
	{
		//Check incrementation direction (if doesn't fit, break the block)				
		eqItem* compArr = equivSet[start+distance].eqList;
		int compArrSize =  equivSet[start+distance].count;

		if(compArrSize != baseArrSize)
			return distance; //mid-function returns still show the valid block length

		for(int k=1; k<baseArrSize; k++)
		{
			if(baseArr[k].dir != compArr[k].dir)
				return distance;
			if(baseArr[k].chrCode != compArr[k].chrCode)
				return distance;
		
			if(baseArr[k].dir == 0 )
			{
				if(baseArr[k].chrPos + distance != compArr[k].chrPos)
				{
					return distance;
				}	
			}
			else
			{
				if(baseArr[k].chrPos - distance != compArr[k].chrPos)
				{
					return distance;
				}
			}
		}
		distance++;
	}

	return distance;
}

#define NO_CLASS_FLAG INT_MAX
unsigned int* eqLookup[MAX_NUM_CHR+2];  //contains lookups to the equivSet indices 0-based
bool* eqLookupDir[MAX_NUM_CHR+2]; //contains directions for the lookups, 0 if forward, 1 if reverse
vector<eqClass> debugUncompressedEquivSet;

void Debug_SetupEqClassLookupTables_And_Print_E0(vector<eqClass>& equivSet, ofstream& fout)
{
	if(equivSet.size() > INT_MAX)
	{
		cout << "ERROR: Debug for detailed print will fail sinze the eqClass vector is too large" << endl;
	}

	for(int i=0; i<=numChrs; i++)
	{
		eqLookup[i] = (unsigned int *) malloc ((chrLens[i] + 2) * sizeof(unsigned int));
		eqLookupDir[i] = (bool *) malloc ((chrLens[i] + 2) * sizeof(bool));
		for(int k=0; k<=chrLens[i]; k++)
		{
			eqLookup[i][k] = NO_CLASS_FLAG; 
			eqLookupDir[i][k] = 0;
		}
	}

	unsigned int numClass = equivSet.size(); 

	//Here Print Eq Classes
	for(unsigned int i=0; i<numClass; i++)
	{
		int count = equivSet[i].count;
		for(int k1=0; k1<count-1; k1++)
		{
			for(int k2=k1+1; k2<count; k2++)
			{
				fout << *(chrNames[(int) equivSet[i].eqList[k1].chrCode]) << "\t" << equivSet[i].eqList[k1].chrPos << "\t" << 16 * (int) (equivSet[i].eqList[k1].dir ^ equivSet[i].eqList[k2].dir) << "\t" << *(chrNames[ (int) equivSet[i].eqList[k2].chrCode]) << "\t" << equivSet[i].eqList[k2].chrPos << "\t" << 0 << endl;
				fout << *(chrNames[(int) equivSet[i].eqList[k2].chrCode]) << "\t" << equivSet[i].eqList[k2].chrPos << "\t" << 16 * (int) (equivSet[i].eqList[k1].dir ^ equivSet[i].eqList[k2].dir) << "\t" << *(chrNames[ (int) equivSet[i].eqList[k1].chrCode]) << "\t" << equivSet[i].eqList[k1].chrPos << "\t" << 0 << endl;
			}
		}
	}

	//Also save a copy of the equivSet
	debugUncompressedEquivSet = equivSet;

	//need deep copy of the eqList
	for(unsigned int i=0; i<debugUncompressedEquivSet.size(); i++)
	{
		int count = debugUncompressedEquivSet[i].count;
		debugUncompressedEquivSet[i].eqList = new eqItem[count];
		for(int k=0; k<count; k++)
		{
			debugUncompressedEquivSet[i].eqList[k] = equivSet[i].eqList[k];
		}				
	}	

	for(unsigned int i=0; i<numClass; i++)
	{
		//cout << "EQUIVSET: " << equivSet[i].count << endl;

		for(int k=0; k<equivSet[i].count; k++)
		{
			eqLookup[ (int) equivSet[i].eqList[k].chrCode ][ equivSet[i].eqList[k].chrPos ] = i;
			if(equivSet[i].eqList[k].dir == 0)
			{
				eqLookupDir[ (int) equivSet[i].eqList[k].chrCode ][ equivSet[i].eqList[k].chrPos ] = 0;
			}
			else
			{
				eqLookupDir[ (int) equivSet[i].eqList[k].chrCode ][ equivSet[i].eqList[k].chrPos ] = 1;
			}
		}	
	}
}

//Go through each equivalence calss and print corresponding eact homologies within the class
long long PrintBinaryEquivData(FILE* foutBinEquiv, int readLen, vector<eqClass>& equivSet)
{
	//The format is a continuous stream of 
	// [rep_chr_no](byte) [rep_chr_pos](int) [dir1](byte) [chrNo1] [chrPos1] [dir2] [chrNo2] [chrPos2] .........(x3) [-1] [block length]
	
	ofstream foutDebug;
	if(DEBUG_SWITCH_FOR_EXACT_PRINT)
	{
		foutDebug.open("DebugEqClassOutput");
	}

	if(DEBUG_SWITCH_FOR_DETAILED_PRINT)
	{
		Debug_SetupEqClassLookupTables_And_Print_E0(equivSet, foutDetailedDebugPrint);
	}


	long long numClass = equivSet.size();	

	long long numPrintedLines = 0;
	
	for(long long i=0; i<numClass; i++)
	{
		numPrintedLines++;
		//Incrementation will be automatic: has to be +1 for same direction or -1 reverse direction (exceptional palindromic-like cases break blocks - so no need to worry)
		unsigned int blockSize = FindEquivBlockSize(i, equivSet); //Perform block merging on the fly while printing

		eqItem* curArr = equivSet[i].eqList;
		int curArrSize = equivSet[i].count;
		
		assert(curArr[0].chrCode != 0);

		fwrite(&(curArr[0].chrCode), sizeof(char), 1, foutBinEquiv);
		fwrite(&(curArr[0].chrPos), sizeof(unsigned int), 1, foutBinEquiv);

		for(int k=1; k<curArrSize; k++) //The boundaries of this for-loop is intentional (due to the way the representative is stored similary to other items but not printed)
                {
			fwrite(&(curArr[k].dir), sizeof(char), 1, foutBinEquiv);

			assert(curArr[k].chrCode != 0);
			fwrite(&(curArr[k].chrCode), sizeof(char), 1, foutBinEquiv);
			fwrite(&(curArr[k].chrPos), sizeof(unsigned int), 1, foutBinEquiv);
                }

		char termin = -1;
		fwrite(&termin, sizeof(char), 1, foutBinEquiv);
		fwrite(&blockSize, sizeof(unsigned int), 1, foutBinEquiv);	
	
		//If debug is on, print text version as well
		if(DEBUG_SWITCH_FOR_EXACT_PRINT)
		{
			foutDebug << (int) curArr[0].chrCode << "\t" << curArr[0].chrPos << "\t";
			for(int k=1; k <curArrSize; k++)
			{
				foutDebug << (int) curArr[k].dir << "\t" << (int) curArr[k].chrCode << "\t" << curArr[k].chrPos << "\t";
			}
			foutDebug << -1 << "\t" << blockSize << endl;
		}

		i += blockSize -1; //minus 1 since it will be incremented once more 
	}

	return numPrintedLines;
}

//Identify total size of the genome from the fa.fai file
long long GetTotalGenomeSize()
{
	long long totalSize = 0;

	for(int i=1; i<=numChrs; i++)
	{
		totalSize += chrLens[i]; 
	}

	return totalSize;
}

int main_exactCompact(int readLength, string InputBinaryEquivClassFile, string OutputBinaryEquivClassFile, long long genSize)
{
	vector<eqClass> equivSetLocal;
	
	equivSetLocal.resize(genSize);

	FillInEquivClasses(InputBinaryEquivClassFile,  equivSetLocal); //Read the text version of equivalences from input and make them compact	
	
	FILE* foutBinEquiv;
	foutBinEquiv = fopen(OutputBinaryEquivClassFile.c_str(), "wb"); 

	long long numPrintedLines = PrintBinaryEquivData(foutBinEquiv, readLength, equivSetLocal); //Print binary versions of the equivalence from in the output

	ofstream foutClassCount((OutputBinaryEquivClassFile + ".itemCount").c_str());
	foutClassCount << numPrintedLines << endl;
	foutClassCount.close();

	ClearEquivClass(equivSetLocal);

	return 1;
}


//###### HOM TABLE INEXACT COLLAPSING #######//
//This version orders each list in the decreasing order of the blockLen, for more efficient access to the nodes within the homology table map inference step.

ofstream foutDebugInexact; //This is only used if we need to debug print the inexact homologies

//Data structure for reading the binary files
//This stores one instance of inexact homology between two read-mers in the reference
struct homPair
{
	char chrCode1;
	char chrCode2;
	unsigned int chrPos1;
	unsigned int chrPos2;
	char dir; //direction of homology
	
	char err1off; //Offsets of the readmer for the mismatch (this is in the forward direction of the first read-mer)
	char err2off; //Offset2 and offset3 are 0 if no error
	char err3off;
};

void DebugPrintHomPair(homPair x)
{
	cout << "HomPair: chr1: " << (int) x.chrCode1 << "  chr2: " << (int) x.chrCode2 << "  pos1: " << x.chrPos1 << "  pos2: " << x.chrPos2;
	cout << "  dir: " << (bool) x.dir << "  err1off: " << (int) x.err1off << "  err2off: " << (int) x.err2off  << " err3off: " << (int) x.err3off << endl;
}

struct homMapItem
{
	char dir;	
	char chrCode;
	int chrPos;
	char offset1;
	char offset2;
	char offset3;
	short blockLen;
};

#define MAX_BLOCK_LEN 32000

//A special triangle item contains dir = -1 , chrCode = 0 (so it's at the beginning), chrPos = triangleId, offset1 = distance triangle, offset2 = -1;

vector<homMapItem>** homMapList[MAX_NUM_CHR+2]; 
int* faiSplitChrCode; //Determines the boundaries of the splitRun for each iteration [0]; boundaries are inclusive; size is numSplitRuns+1; [0] is 1 & [numSplitRuns] = numChrs is the end chr of the splitRun;

void PrepHomMapList(int& numSplitRuns)
{
	//Initialize homMapListForAllChrSizes (from the fai files)
	//Initialize each vetor as null, so they can be initialized only when needed	

	//Depending on chrLengths, here determine the splitRun boundaries
	faiSplitChrCode = (int *) malloc ((numSplitRuns + 1) * sizeof(int));

	faiSplitChrCode[0] = 0; //note that chrCodes are 1-based (but this actually represents the end chr for the -1 th split (which is 0))

	long long totalChrLength = 0;
	for(int i=1; i<=numChrs; i++)
	{
		totalChrLength += chrLens[i];
	}

	long long aimedSplitLen = totalChrLength / (long long)numSplitRuns;
	
	long long lengthSoFar = 0;
	long long nextThreshold = aimedSplitLen;
	int splitNo = 1;

	for(int i=1; i<=numChrs; i++)
	{
		if(i >= numChrs)
		{
			numSplitRuns = splitNo;
			break;
		}	
		lengthSoFar += chrLens[i];
		if(lengthSoFar >= nextThreshold)
		{
			faiSplitChrCode[splitNo] = i;
			lengthSoFar = 0;
			splitNo++;
		}
	}
	faiSplitChrCode[numSplitRuns] = numChrs;
}

void InitSplitAllocation(int startCtg, int endCtg)
{
	for(int chrCode=startCtg; chrCode<=endCtg; chrCode++)
	{
		homMapList[chrCode] = (vector<homMapItem>**) calloc (chrLens[chrCode]+2, sizeof(vector<homMapItem>*));
	}
}

void ClearSplitAllocation_ContigIntervalAware(int startCtg, int endCtg, int contigStart, int contigEnd)
{
	assert(startCtg == endCtg);

	//Other clearing might need to be done here
	for(int i=startCtg; i<=endCtg; i++)
        {
                for(int k=contigStart; k<= contigEnd; k++)
                {
                        vector<homMapItem>* curVec = homMapList[i][k];
			if(curVec)
			{
				curVec->clear();
			}
			delete curVec;
		}
	}

        for(int chrCode=startCtg; chrCode<=endCtg; chrCode++)
        {
                free(homMapList[chrCode]);
        }
}

void ClearSplitAllocation(int startCtg, int endCtg)
{
	//Other clearing might need to be done here
	for(int i=startCtg; i<=endCtg; i++)
        {
                for(int k=1; k<= chrLens[i]; k++)
                {
                        vector<homMapItem>* curVec = homMapList[i][k];
			if(curVec)
			{
				curVec->clear();
			}
			delete curVec;
		}
	}

        for(int chrCode=startCtg; chrCode<=endCtg; chrCode++)
        {
                free(homMapList[chrCode]);
        }
}

homPair SortHomPair(homPair inPair, int readLen)
{
	homPair outPair = inPair;
	if(inPair.chrCode2 < inPair.chrCode1 || (inPair.chrCode2 == inPair.chrCode1 && inPair.chrPos2 < inPair.chrPos1))
	{
		//Perform reversal here
		outPair.chrCode1 = inPair.chrCode2;
		outPair.chrCode2 = inPair.chrCode1;
		outPair.chrPos1 = inPair.chrPos2;
		outPair.chrPos2 = inPair.chrPos1;
		
		if(inPair.dir)
		{
			assert(inPair.err1off != 0);
			//invert offsets if they are in reverse 
			if(inPair.err2off == 0)
			{
				outPair.err1off = readLen - inPair.err1off + 1;
				outPair.err2off = 0;
				outPair.err3off = 0;
			}
			else if(inPair.err3off == 0)
			{
				outPair.err1off = readLen - inPair.err2off + 1;
				outPair.err2off = readLen - inPair.err1off + 1;
				outPair.err3off = 0;
			}
			else
			{
				outPair.err1off = readLen - inPair.err3off + 1;
				outPair.err2off = readLen - inPair.err2off + 1;
				outPair.err3off = readLen - inPair.err1off + 1;
			}
		}
	}
	return outPair;
}

void DebugDetailedPrint(int sourceChr, int sourcePos, const homMapItem& curItem, ofstream& fout)
{
	int targetChr = curItem.chrCode;
	int targetPos = curItem.chrPos; 		
			
	int mapDir = curItem.dir;

	int numErrors = 0;
	if(curItem.offset1 != 0) 
		numErrors++;
	else 
		assert(curItem.offset2 == 0 && curItem.offset3 == 0);

	if(curItem.offset2 != 0) 
		numErrors++;
	else 
		assert(curItem.offset3 == 0);

	if(curItem.offset3 != 0) 
		numErrors++;
		
	int sourceEqClassId = eqLookup[sourceChr][sourcePos];
	int sourceDir = eqLookupDir[sourceChr][sourcePos];

	int targetEqClassId = eqLookup[targetChr][targetPos];
	int targetDir = eqLookupDir[targetChr][targetPos];
	
	int jointDir = (mapDir ^ sourceDir) ^ targetDir;  

	if(sourceEqClassId == NO_CLASS_FLAG && targetEqClassId == NO_CLASS_FLAG)
	{
		fout << *(chrNames[sourceChr]) << "\t" << sourcePos << "\t" << mapDir * 16 << "\t" << *(chrNames[targetChr]) << "\t" << targetPos << "\t" << numErrors << endl;
	}
	if(sourceEqClassId == NO_CLASS_FLAG && targetEqClassId != NO_CLASS_FLAG)
	{
		for(int j2 = 0; j2< debugUncompressedEquivSet[targetEqClassId].count; j2++)
		{
			int chrSecond = debugUncompressedEquivSet[targetEqClassId].eqList[j2].chrCode;
			int posSecond = debugUncompressedEquivSet[targetEqClassId].eqList[j2].chrPos;
	
			int printDir = jointDir ^ debugUncompressedEquivSet[targetEqClassId].eqList[j2].dir;

			fout << *(chrNames[sourceChr]) << "\t" << sourcePos << "\t" << printDir *  16 << "\t" << *(chrNames[chrSecond]) << "\t" << posSecond << "\t" << numErrors << endl; 			
		}
	}

	if(sourceEqClassId != NO_CLASS_FLAG && targetEqClassId == NO_CLASS_FLAG)
	{
		for(int j1=0; j1< debugUncompressedEquivSet[sourceEqClassId].count; j1++)
		{
			int chrFirst = debugUncompressedEquivSet[sourceEqClassId].eqList[j1].chrCode;
			int posFirst = debugUncompressedEquivSet[sourceEqClassId].eqList[j1].chrPos;

			int printDir = jointDir ^ debugUncompressedEquivSet[sourceEqClassId].eqList[j1].dir;

			fout << *(chrNames[chrFirst]) << "\t" << posFirst << "\t" << printDir * 16 << "\t" << *(chrNames[targetChr]) << "\t" << targetPos << "\t" << numErrors << endl;			
		}
	}

	if(sourceEqClassId != NO_CLASS_FLAG && targetEqClassId != NO_CLASS_FLAG)
	{		
		for(int j1=0; j1< debugUncompressedEquivSet[sourceEqClassId].count; j1++)
		{
			for(int j2 = 0; j2< debugUncompressedEquivSet[targetEqClassId].count; j2++)
			{
				int chrFirst = debugUncompressedEquivSet[sourceEqClassId].eqList[j1].chrCode;
				int chrSecond = debugUncompressedEquivSet[targetEqClassId].eqList[j2].chrCode;
				int posFirst = debugUncompressedEquivSet[sourceEqClassId].eqList[j1].chrPos;
				int posSecond = debugUncompressedEquivSet[targetEqClassId].eqList[j2].chrPos;
		
				int ultimateDir = (debugUncompressedEquivSet[sourceEqClassId].eqList[j1].dir ^ jointDir) ^ debugUncompressedEquivSet[targetEqClassId].eqList[j2].dir;

				fout << *(chrNames[chrFirst]) << "\t" << posFirst << "\t" << ultimateDir *  16 << "\t" << *(chrNames[chrSecond]) << "\t" << posSecond << "\t" << numErrors << endl; 			
			}
		}
	}
}

void LoadHomPair_SplitAware_ContigIntervalAware(homPair curHomPair, int readLength, int splitChrStart, int splitChrEnd, unsigned int contigStart, unsigned int contigEnd) //This is the same function as the LoadHomPair_SplitAware but allows limitation on a subinterval of a single chromosome
{
	//Look at the first position
	//Create the vector if necessary
	//Push back a hom map item, depending on the joint directionality
	//Since this will be printed to a file, can keep a single direction so far (but it might turn out that compacting works better when both are included)
	
	assert(splitChrStart == splitChrEnd); //This function only allows a single chromosome's subintervals.

	if(curHomPair.chrCode1 >= splitChrStart && curHomPair.chrCode1 <= splitChrEnd && curHomPair.chrPos1 >= contigStart && curHomPair.chrPos1 <= contigEnd)
	{
		vector<homMapItem>*& firstToSecVec = homMapList[(int) curHomPair.chrCode1][curHomPair.chrPos1];
		if(firstToSecVec == NULL)
		{
			firstToSecVec = new vector<homMapItem>(0);
		}
		
		homMapItem firstToSecondMap;
		firstToSecondMap.chrCode = curHomPair.chrCode2;
		firstToSecondMap.chrPos = curHomPair.chrPos2;
		firstToSecondMap.offset1 = curHomPair.err1off;
		firstToSecondMap.offset2 = curHomPair.err2off;
		firstToSecondMap.offset3 = curHomPair.err3off;
		firstToSecondMap.dir = curHomPair.dir;	
		firstToSecondMap.blockLen = 1;
		
		firstToSecVec->push_back(firstToSecondMap);

		if(DEBUG_SWITCH_FOR_DETAILED_PRINT)
		{
			DebugDetailedPrint(curHomPair.chrCode1, curHomPair.chrPos1, firstToSecondMap, foutDetailedDebugPrint);
		}		
	}

	if(curHomPair.chrCode2 >= splitChrStart && curHomPair.chrCode2 <= splitChrEnd && curHomPair.chrPos2 >= contigStart && curHomPair.chrPos2 <= contigEnd)
	{
		vector<homMapItem>*& SecToFirstVec = homMapList[(int) curHomPair.chrCode2][curHomPair.chrPos2];

		if(SecToFirstVec == NULL)
		{
			SecToFirstVec = new vector<homMapItem>(0);
		}

		homMapItem secondToFirstMap;
		secondToFirstMap.chrCode = curHomPair.chrCode1;
		secondToFirstMap.chrPos = curHomPair.chrPos1;
		secondToFirstMap.dir = curHomPair.dir;
		secondToFirstMap.blockLen = 1;
		if(curHomPair.dir)
		{
			if(curHomPair.err2off == 0)
			{
				secondToFirstMap.offset1 = readLength - curHomPair.err1off + 1;
				secondToFirstMap.offset2 = 0;
				secondToFirstMap.offset3 = 0;
			}
			else if(curHomPair.err3off == 0)
			{
				secondToFirstMap.offset1 = readLength - curHomPair.err2off + 1;
				secondToFirstMap.offset2 = readLength - curHomPair.err1off + 1;
				secondToFirstMap.offset3 = 0;
			}
			else
			{
				secondToFirstMap.offset1 = readLength - curHomPair.err3off + 1;
				secondToFirstMap.offset2 = readLength - curHomPair.err2off + 1;
				secondToFirstMap.offset3 = readLength - curHomPair.err1off + 1;
			}
		}
		else
		{
			secondToFirstMap.offset1 = curHomPair.err1off;
			secondToFirstMap.offset2 = curHomPair.err2off;
			secondToFirstMap.offset3 = curHomPair.err3off;
		}

		SecToFirstVec->push_back(secondToFirstMap);
		
		if(DEBUG_SWITCH_FOR_DETAILED_PRINT)
		{
			DebugDetailedPrint(curHomPair.chrCode2, curHomPair.chrPos2, secondToFirstMap, foutDetailedDebugPrint);
		}		
	}
}

void LoadHomPair_SplitAware(homPair curHomPair, int readLength, int splitChrStart, int splitChrEnd)
{
	//Look at the first position
	//Create the vector if necessary
	//Push back a hom map item, depending on the joint directionality
	//Since this will be printed to a file, can keep a single direction so far (but it might turn out that compacting works better when both are included)
	
	if(curHomPair.chrCode1 >= splitChrStart && curHomPair.chrCode1 <= splitChrEnd)
	{
		vector<homMapItem>*& firstToSecVec = homMapList[(int) curHomPair.chrCode1][curHomPair.chrPos1];
		if(firstToSecVec == NULL)
		{
			firstToSecVec = new vector<homMapItem>(0);
		}
		
		homMapItem firstToSecondMap;
		firstToSecondMap.chrCode = curHomPair.chrCode2;
		firstToSecondMap.chrPos = curHomPair.chrPos2;
		firstToSecondMap.offset1 = curHomPair.err1off;
		firstToSecondMap.offset2 = curHomPair.err2off;
		firstToSecondMap.offset3 = curHomPair.err3off;
		firstToSecondMap.dir = curHomPair.dir;	
		firstToSecondMap.blockLen = 1;
		
		firstToSecVec->push_back(firstToSecondMap);

		if(DEBUG_SWITCH_FOR_DETAILED_PRINT)
		{
			DebugDetailedPrint(curHomPair.chrCode1, curHomPair.chrPos1, firstToSecondMap, foutDetailedDebugPrint);
		}		
	}

	if(curHomPair.chrCode2 >= splitChrStart && curHomPair.chrCode2 <= splitChrEnd)
	{
		vector<homMapItem>*& SecToFirstVec = homMapList[(int) curHomPair.chrCode2][curHomPair.chrPos2];

		if(SecToFirstVec == NULL)
		{
			SecToFirstVec = new vector<homMapItem>(0);
		}

		homMapItem secondToFirstMap;
		secondToFirstMap.chrCode = curHomPair.chrCode1;
		secondToFirstMap.chrPos = curHomPair.chrPos1;
		secondToFirstMap.dir = curHomPair.dir;
		secondToFirstMap.blockLen = 1;
		if(curHomPair.dir)
		{
			if(curHomPair.err2off == 0)
			{
				secondToFirstMap.offset1 = readLength - curHomPair.err1off + 1;
				secondToFirstMap.offset2 = 0;
				secondToFirstMap.offset3 = 0;
			}
			else if(curHomPair.err3off == 0)
			{
				secondToFirstMap.offset1 = readLength - curHomPair.err2off + 1;
				secondToFirstMap.offset2 = readLength - curHomPair.err1off + 1;
				secondToFirstMap.offset3 = 0;
			}
			else
			{
				secondToFirstMap.offset1 = readLength - curHomPair.err3off + 1;
				secondToFirstMap.offset2 = readLength - curHomPair.err2off + 1;
				secondToFirstMap.offset3 = readLength - curHomPair.err1off + 1;
			}
		}
		else
		{
			secondToFirstMap.offset1 = curHomPair.err1off;
			secondToFirstMap.offset2 = curHomPair.err2off;
			secondToFirstMap.offset3 = curHomPair.err3off;
		}

		SecToFirstVec->push_back(secondToFirstMap);
		
		if(DEBUG_SWITCH_FOR_DETAILED_PRINT)
		{
			DebugDetailedPrint(curHomPair.chrCode2, curHomPair.chrPos2, secondToFirstMap, foutDetailedDebugPrint);
		}		
	}
}

bool sort_by_genomic_pos(homMapItem h1, homMapItem h2)
{
	return h1.chrCode == h2.chrCode ? (h1.chrPos < h2.chrPos) : h1.chrCode < h2.chrCode; //First sort by chrCode and then by chrPos
}

bool sort_by_decreasing_blockLen(homMapItem h1, homMapItem h2)
{
	return h1.blockLen > h2.blockLen;
}

void SortHomMapList_SplitAware_ContigIntervalAware(int splitChrStart, int splitChrEnd, int contigStart, int contigEnd)
{
	assert(splitChrStart == splitChrEnd);

	//go vector by vector and sort each vector (since most lists would be small, this wont take too long)
	for(int i=splitChrStart; i<=splitChrEnd; i++)
	{
		for(int k=contigStart; k<= contigEnd; k++)
		{
			if(homMapList[i][k])
			{
				sort(homMapList[i][k]->begin(), homMapList[i][k]->end(), sort_by_genomic_pos);
			}
		}
	}
}

void SortHomMapList_SplitAware(int splitChrStart, int splitChrEnd)
{
	//go vector by vector and sort each vector (since most lists would be small, this wont take too long)
	for(int i=splitChrStart; i<=splitChrEnd; i++)
	{
		for(int k=1; k<= chrLens[i]; k++)
		{
			if(homMapList[i][k])
			{
				sort(homMapList[i][k]->begin(), homMapList[i][k]->end(), sort_by_genomic_pos);
			}
		}
	}
}

int GetAdjacentIntervalEndPoint(int chrCode, int startPos)
{
	int distance = 1;
	
	while(homMapList[chrCode][startPos + distance])
	{
		distance++;	
	}

	return startPos + distance - 1; //-1 because last try failed
}

int SearchInVec(const vector<homMapItem>& vec, int lastPos, int chrCode, int dir, int chrPos, int offset1, int offset2, int offset3)
{
	int vecSize = vec.size();
	int curInd = lastPos;

	while(curInd < vecSize && (vec[curInd].chrCode < chrCode || (vec[curInd].chrCode == chrCode && vec[curInd].chrPos < chrPos)))
	{
		curInd++;
	}	

	if(curInd >= vecSize)
	{
		curInd = vecSize - 1;
	}
	
	while(curInd >= 0 && (vec[curInd].chrCode > chrCode || (vec[curInd].chrCode == chrCode && vec[curInd].chrPos > chrPos)))
	{
		curInd--;
	}

	if(curInd < 0)
	{
		curInd = 0;
	}

	//At this point if there is a single chrCode/chrPos with proper offsets it is checked
	if(vec[curInd].chrCode != chrCode || vec[curInd].chrPos != chrPos || vec[curInd].dir != dir || vec[curInd].offset1 != offset1 || vec[curInd].offset2 != offset2 || vec[curInd].offset3 != offset3)
	{
		return -1;
	}

	return curInd;
}

void DebugPrintHomMapItemVector(vector<homMapItem>*& vectorToPrint)
{
	vector<homMapItem>& vec = *(vectorToPrint);
	int size = vec.size();
	for(int i=0; i<size; i++)
	{
		cout << i << ": " << (int) vec[i].chrCode << " " << (int) vec[i].dir << " " << (int) vec[i].chrPos << " " << (int) vec[i].offset1 << " " << (int) vec[i].offset2 << " " << (int) vec[i].offset3 << " " << (int) vec[i].blockLen << endl;
	}
}

void CollapseSingleNodeIntervalBlocks(int chrCode, int startPos) //Note that blockLen is limited by 32000 (blockLen -1 means that the node is deactivated)
{
	int distance = 1;
	
	assert(homMapList[chrCode][startPos]);
	vector<homMapItem>& baseList = *(homMapList[chrCode][startPos]);

	int baseListSize = baseList.size();	

	while(homMapList[chrCode][startPos + distance])
	{
		//Check incrementation direction (if doesn't fit, break the block)				
		vector<homMapItem>& compList = *(homMapList[chrCode][startPos + distance]);

		int lastLookedPositionInCompList = 0;
		for(int k=0; k<baseListSize; k++)
		{
			if(baseList[k].blockLen == distance && baseList[k].blockLen <= 65000)
			{
				//search it in current list (binary search could be done - or prev pos could be save and two way search be done on the saved position)
				int chrPosToLookFor = -1;

				if(baseList[k].dir == 0)
				{
					chrPosToLookFor = baseList[k].chrPos + distance;

				}
				else
				{
					chrPosToLookFor = baseList[k].chrPos - distance;

				}

				if(chrPosToLookFor <= 0 || chrPosToLookFor > chrLens[(int) baseList[k].chrCode])
				{
					continue;
				}
				assert(chrPosToLookFor > 0);
				
				int offset1ToLookFor = baseList[k].offset1 - distance;
				int offset2ToLookFor = 0;
				if(baseList[k].offset2 != 0) 
				{
					offset2ToLookFor = baseList[k].offset2 - distance;
				}
				int offset3ToLookFor = 0;
				if(baseList[k].offset3 != 0)
				{
					offset3ToLookFor = baseList[k].offset3 - distance;
				}
				
				if(offset1ToLookFor <= 0 || offset2ToLookFor < 0 || offset3ToLookFor < 0) //Not one if small-equal the other is strictly small (since one error should exist)
				{
					continue;
				}
				
				int compVecPos = SearchInVec(compList, lastLookedPositionInCompList, baseList[k].chrCode, baseList[k].dir, chrPosToLookFor, offset1ToLookFor, offset2ToLookFor, offset3ToLookFor);
				if(compVecPos != -1)
				{
					lastLookedPositionInCompList = compVecPos;
					baseList[k].blockLen++;
					compList[compVecPos].blockLen = -1;
				}
			}
		}
		distance++;
	}
}

unsigned int PrintBinaryHomMapData_SplitAware_ContigIntervalAware(FILE* foutBinMap, int splitChrStart, int splitChrEnd, int numMismatches, int contigStart, int contigEnd) //optionally can compact here
{
	assert(splitChrStart == splitChrEnd);

	unsigned int numberOfLinesPrinted = 0;
	for(int i=splitChrStart; i<=splitChrEnd; i++)
        {
		cout << "Started printing chrCode: " << i << endl;
                for(int k=contigStart; k<=contigEnd; k++)
                {
			vector<homMapItem>* curVec = homMapList[i][k];	
                        
			if(curVec)
			{
				CollapseSingleNodeIntervalBlocks(i, k);
				
				//Sort so that the longer blockLen is at the beginning
				sort(curVec->begin(), curVec->end(), sort_by_decreasing_blockLen);
		
				int noHeaderYetFLAG = 1;
				int size = curVec->size();
				for(int p=0; p<size; p++)
				{
					if((*curVec)[p].blockLen != -1)
					{
						if(noHeaderYetFLAG)
						{
							fwrite(&i, sizeof(char), 1, foutBinMap);
							fwrite(&k, sizeof(unsigned int), 1, foutBinMap);			
						
							if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
							{
								foutDebugInexact << "HEADER: " << i << " " << k << endl; 
							}
	
							noHeaderYetFLAG = 0;
						}
						fwrite(&((*curVec)[p].dir), sizeof(char), 1, foutBinMap);
						fwrite(&((*curVec)[p].chrCode), sizeof(char), 1, foutBinMap);				
						
						assert((*curVec)[p].chrCode != 0);

						fwrite(&((*curVec)[p].chrPos), sizeof(unsigned int), 1, foutBinMap);
						fwrite(&((*curVec)[p].offset1), sizeof(char), 1, foutBinMap);				
						if(numMismatches > 1)
						{
							fwrite(&((*curVec)[p].offset2), sizeof(char), 1, foutBinMap);				
						}
						if(numMismatches > 2)
						{
							fwrite(&((*curVec)[p].offset3), sizeof(char), 1, foutBinMap);				
						}
						fwrite(&((*curVec)[p].blockLen), sizeof(short), 1, foutBinMap);				
				

						if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
						{
							foutDebugInexact << (int) (*curVec)[p].dir << "\t" << (int) (*curVec)[p].chrCode << "\t" << (*curVec)[p].chrPos << "\t" << (int) (*curVec)[p].offset1;
							if(numMismatches > 1)
								foutDebugInexact << "\t" << (int) (*curVec)[p].offset2;
							if(numMismatches > 2)
								foutDebugInexact << "\t" << (int) (*curVec)[p].offset3;
							foutDebugInexact << "\t" << (int) (*curVec)[p].blockLen << endl;;
						}
					}
				}
				
				if(noHeaderYetFLAG == 0)
				{
					numberOfLinesPrinted++;			
					char endValue = -1;
					fwrite(&endValue, sizeof(char), 1, foutBinMap);
					if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
					{
						foutDebugInexact << "End Value: " << (int) endValue << endl; 
					}
				}
			}
		}
	}
	return numberOfLinesPrinted;
}

unsigned int PrintBinaryHomMapData_SplitAware(FILE* foutBinMap, int splitChrStart, int splitChrEnd, int numMismatches) //optionally can compact here
{

	unsigned int numberOfLinesPrinted = 0;
	for(int i=splitChrStart; i<=splitChrEnd; i++)
        {
		cout << "Started printing chrCode: " << i << endl;
                for(int k=1; k<= chrLens[i]; k++)
                {
			vector<homMapItem>* curVec = homMapList[i][k];	
                        
			if(curVec)
			{
				CollapseSingleNodeIntervalBlocks(i, k);
				
				//Sort so that the longer blockLen is at the beginning
				sort(curVec->begin(), curVec->end(), sort_by_decreasing_blockLen);
		
				int noHeaderYetFLAG = 1;
				int size = curVec->size();
				for(int p=0; p<size; p++)
				{
					if((*curVec)[p].blockLen != -1)
					{
						if(noHeaderYetFLAG)
						{
							fwrite(&i, sizeof(char), 1, foutBinMap);
							fwrite(&k, sizeof(unsigned int), 1, foutBinMap);			
						
							if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
							{
								foutDebugInexact << "HEADER: " << i << " " << k << endl; 
							}
	
							noHeaderYetFLAG = 0;
						}
						fwrite(&((*curVec)[p].dir), sizeof(char), 1, foutBinMap);
						fwrite(&((*curVec)[p].chrCode), sizeof(char), 1, foutBinMap);				
						
						assert((*curVec)[p].chrCode != 0);

						fwrite(&((*curVec)[p].chrPos), sizeof(unsigned int), 1, foutBinMap);
						fwrite(&((*curVec)[p].offset1), sizeof(char), 1, foutBinMap);				
						if(numMismatches > 1)
						{
							fwrite(&((*curVec)[p].offset2), sizeof(char), 1, foutBinMap);				
						}
						if(numMismatches > 2)
						{
							fwrite(&((*curVec)[p].offset3), sizeof(char), 1, foutBinMap);				
						}
						fwrite(&((*curVec)[p].blockLen), sizeof(short), 1, foutBinMap);				
				

						if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
						{
							foutDebugInexact << (int) (*curVec)[p].dir << "\t" << (int) (*curVec)[p].chrCode << "\t" << (*curVec)[p].chrPos << "\t" << (int) (*curVec)[p].offset1;
							if(numMismatches > 1)
								foutDebugInexact << "\t" << (int) (*curVec)[p].offset2;
							if(numMismatches > 2)
								foutDebugInexact << "\t" << (int) (*curVec)[p].offset3;
							foutDebugInexact << "\t" << (int) (*curVec)[p].blockLen << endl;;
						}
					}
				}
				
				if(noHeaderYetFLAG == 0)
				{
					numberOfLinesPrinted++;			
					char endValue = -1;
					fwrite(&endValue, sizeof(char), 1, foutBinMap);
					if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
					{
						foutDebugInexact << "End Value: " << (int) endValue << endl; 
					}
				}
			}
		}
	}
	
	return numberOfLinesPrinted;
}

int main_inexactCompact(int numThreads, int numSplitRuns, int readLength, string inputHomMapPrefix, string outputBinaryHomMap, int numMismatches, int maxContigLen) //if maxContigLen is 0, then there is no contigLen based splitting
{
	if(DEBUG_SWITCH_FOR_INEXACT_PRINT)
	{
		foutDebugInexact.open("DebugInexactOutput");
	}

	FILE* foutBinMap;
	foutBinMap = fopen(outputBinaryHomMap.c_str(), "wb");

	//This file will be used to print out auxiliary info to the hom table (such as the total number of class nodes, emHomItems, etc.)
	string binMapitemCountFile = outputBinaryHomMap + ".itemCount";
	ofstream foutBinMap_itemCounts(binMapitemCountFile.c_str());

	if(numSplitRuns < 1)
	{
		cout << "ERROR: number of split runs >= 1" << endl;
		exit(6);
	}

	//Here we setup the chrCodes an chrLengths, as well as depending on the number of split-runs, we determine the staring and ending chromosome of each split (For 0-based xth split, we'll process from faiSplitChrCode[x]+1 till faiSplitChrCode[x+1])
	
	cout << "Before Prep #splitRuns is: " << numSplitRuns << endl;
	PrepHomMapList(numSplitRuns);
	cout << "After Prep #splitRuns is: " << numSplitRuns << endl;

	unsigned int numberOfLinesPrinted = 0;

	for(int splitRunNo = 0; splitRunNo < numSplitRuns; splitRunNo++)
	{
		int splitFaiStart = faiSplitChrCode[splitRunNo] + 1; //These are the chrCodes that will be handled this split run
		int splitFaiEnd = faiSplitChrCode[splitRunNo+1];

		long long totalContigLen = 0;
		for(int k=splitFaiStart; k<=splitFaiEnd; k++)
		{
			totalContigLen = chrLens[k];
		}
		assert(totalContigLen > 0);

		if(maxContigLen != 0 && maxContigLen < totalContigLen)
		{
			//Then go through chromosomes one by one and do subchromosomal splits when needed
			for(int curChrNo=splitFaiStart; curChrNo<=splitFaiEnd; curChrNo++)
			{
				if(chrLens[curChrNo] <= maxContigLen) //process as usual
				{
					cout << "Split Run No: " << splitRunNo << " SubSplit Whole:  ctgStart: " << curChrNo << "  ctgEnd: " << curChrNo << endl;
					
					InitSplitAllocation(curChrNo, curChrNo);

					for(int i=1; i<=numThreads; i++)
					{
						stringstream curFileNameSS;
						curFileNameSS << inputHomMapPrefix << i;
						FILE* finHomMap = fopen(curFileNameSS.str().c_str(), "rb");
					
						if(finHomMap == NULL)
						{
							cout << "File can't be opened: " << curFileNameSS.str() << endl;
							exit(8); 	
						}
				
						homPair curHomPair;
						
						while(true)
						{
							int result = fread(&curHomPair, sizeof(homPair), 1, finHomMap);
							if(result != 1)
							{
								break;
							}
						
							//Add both forward links (first -> second) and reverse links (second -> first)	
							LoadHomPair_SplitAware(curHomPair, readLength, curChrNo, curChrNo);
						}
						cout << "Finished reading file: " << curFileNameSS.str() << "..." << endl;
				
						fclose(finHomMap);
					}

					//Starting sorting
					cout << "Sorting homMapList.." << endl;
					SortHomMapList_SplitAware(curChrNo, curChrNo);

					cout << "Printing binary homology data" << endl;
					numberOfLinesPrinted += PrintBinaryHomMapData_SplitAware(foutBinMap, curChrNo, curChrNo, numMismatches); //optionally can compact here
				
					cout << "Before clearing" << endl;
					ClearSplitAllocation(curChrNo, curChrNo);
					cout << "After clearing" << endl;
				}
				else //then employ splits
				{
					int contigStart = 1;

					while(contigStart <= chrLens[curChrNo])
					{
						int contigEnd = min(chrLens[curChrNo], contigStart + maxContigLen - 1);

						cout << "Split Run No: " << splitRunNo << " SubSplit Partial:  ctgStart: " << curChrNo << "  ctgEnd: " << curChrNo << " contigStart: " << contigStart << " contigEnd: " << contigEnd << endl;
					
						InitSplitAllocation(curChrNo, curChrNo);

						for(int i=1; i<=numThreads; i++)
						{
							stringstream curFileNameSS;
							curFileNameSS << inputHomMapPrefix << i;
							FILE* finHomMap = fopen(curFileNameSS.str().c_str(), "rb");
						
							if(finHomMap == NULL)
							{
								cout << "File can't be opened: " << curFileNameSS.str() << endl;
								exit(8); 	
							}
					
							homPair curHomPair;
							
							while(true)
							{
								int result = fread(&curHomPair, sizeof(homPair), 1, finHomMap);
								if(result != 1)
								{
									break;
								}
							
								//Add both forward links (first -> second) and reverse links (second -> first)	
								LoadHomPair_SplitAware_ContigIntervalAware(curHomPair, readLength, curChrNo, curChrNo, contigStart, contigEnd);
							}
							cout << "Finished reading file: " << curFileNameSS.str() << "..." << endl;
					
							fclose(finHomMap);
						}

						//Starting sorting
						cout << "Sorting homMapList.." << endl;
						SortHomMapList_SplitAware_ContigIntervalAware(curChrNo, curChrNo, contigStart, contigEnd); //There is nothing position specific in sorting

						cout << "Printing binary homology data" << endl;
						numberOfLinesPrinted += PrintBinaryHomMapData_SplitAware_ContigIntervalAware(foutBinMap, curChrNo, curChrNo, numMismatches, contigStart, contigEnd); //optionally can compact here
					
						
						cout << "Before clearing" << endl;
						ClearSplitAllocation_ContigIntervalAware(curChrNo, curChrNo, contigStart, contigEnd);
						cout << "After clearing" << endl;
					
						contigStart += maxContigLen;
					}
				}
			} 
		}
		else
		{
			cout << "Split Run No: " << splitRunNo << " ctgStart: " << splitFaiStart << "  ctgEnd: " << splitFaiEnd << endl;
			
			InitSplitAllocation(splitFaiStart, splitFaiEnd);

			for(int i=1; i<=numThreads; i++)
			{
				stringstream curFileNameSS;
				curFileNameSS << inputHomMapPrefix << i;
				FILE* finHomMap = fopen(curFileNameSS.str().c_str(), "rb");
			
				if(finHomMap == NULL)
				{
					cout << "File can't be opened: " << curFileNameSS.str() << endl;
					exit(8); 	
				}
		
				homPair curHomPair;
				
				while(true)
				{
					int result = fread(&curHomPair, sizeof(homPair), 1, finHomMap);
					if(result != 1)
					{
						break;
					}
				
					//Add both forward links (first -> second) and reverse links (second -> first)	
					LoadHomPair_SplitAware(curHomPair, readLength, splitFaiStart, splitFaiEnd);
				}
				cout << "Finished reading file: " << curFileNameSS.str() << "..." << endl;
		
				fclose(finHomMap);
			}

			//Starting sorting
			cout << "Sorting homMapList.." << endl;
			SortHomMapList_SplitAware(splitFaiStart, splitFaiEnd);

			cout << "Printing binary homology data" << endl;
			numberOfLinesPrinted += PrintBinaryHomMapData_SplitAware(foutBinMap, splitFaiStart, splitFaiEnd, numMismatches); //optionally can compact here
		
			cout << "Before clearing" << endl;
			ClearSplitAllocation(splitFaiStart, splitFaiEnd);
			cout << "After clearing" << endl;
		}
	}

	foutBinMap_itemCounts << numberOfLinesPrinted << endl;
	foutBinMap_itemCounts.close();

	return 1;
}

// Variables for seed (segment) selection
char SegmentTemplate[MAX_READ_LEN+2]; //This template contains shapes for the splits -- useful for simple error checking - Note that the templates are 0-based
char subSegmentLength; // If segment sizes are shorter than 8bp, this handles dummy characters
unsigned short shortRevCompArr[USHRT_MAX + 5]; //Look up table that holds revcomps of all short values

//###################
// Data types for the hash table
//###################

#define numberOfLongs 4

struct nLongs //Holds multiple long longs representing DNA sequences [Make this compiler argument later on]
{
	nLongs()
	{
		for(unsigned char i=0; i<numberOfLongs; i++)
		{
			longList[i] = 0;
		}
	}

	unsigned long long longList[numberOfLongs]; //each one of them stores a 16-mer
};

struct nLongsHash //Hash functor 
{
	long operator() (const nLongs &k) const 
	{ 
		long retLong = 0;
		for(unsigned char i=0; i<numberOfLongs; i++)
		{
			retLong += k.longList[i];
		}
		return retLong;
	}
};

struct nLongsPred //Equality functor
{
	bool operator() (const nLongs &x, const nLongs &y) const 
	{
		for(unsigned char i=0; i<numberOfLongs; i++)
		{
			if(x.longList[i] != y.longList[i])
			{
				return false;
			}
		}
		return true;
	}
};
///////////////////////////////////////

string ReverseComplementString(string seq)
{
        string revCompSeq ((size_t) seq.length(), 'X');
        int len = seq.length();

        for(int i=0; i<len; i++)
        {
                revCompSeq[i] = RevCompChar[seq[len-i-1]];
        }
        return revCompSeq;
}

//printing text and integer versions of nLong 
string nLong2Binary_Debug(nLongs curNLong)
{
	stringstream ss;

	for(unsigned char i=0; i<numberOfLongs; i++)
	{
		std::bitset<64> before(curNLong.longList[i]);
		ss << before << " ";
	}	

	ss << "\nLL: ";

	for(unsigned char i=0; i<numberOfLongs; i++)
	{
		ss << curNLong.longList[i] << " ";
	}

	ss << endl;
	
	return ss.str(); 
}

char GetNucFromCode(int code)
{
	switch(code)
	{	
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		default:
			cout << "ERROR(3) Nucleotide code isn't valid: " << code << endl;
			exit(0);
	}	
}

unsigned int GetValueForNuc(char ch)
{
        switch(ch)
        {
                case 'A':
                        return 0;
                case 'C':
                        return 1;
                case 'G':
                        return 2;
                case 'T':
                        return 3;
               	default:
                        cout << "ERROR(2) character: '" << ch << "' ASCII: " << (int) ch << endl;
                        exit(0);
        }
}

//convert a read-mer sequence into an nLong
nLongs GetNLongsFromString(const char *arr, int readLen)
{
	assert(readLen <= numberOfLongs*32);

	nLongs curNLongs;

	unsigned long long val;
	short llBitLen = 0;
	short llPos = 0; //this is iterated over pos
	for(unsigned char k=0; k<numberOfLongs && llBitLen < readLen; k++)
	{
		val = 0;
		llBitLen = min(readLen, llBitLen+32);

		for(; llPos<llBitLen; llPos++)
		{
			val = val << 2;
			val += GetValueForNuc(arr[llPos]);
		}
	
		if(llBitLen % 32 != 0) //this automatically be the last step
		{
			char amountToShiftHalf = 32 - (llBitLen % 32);
			curNLongs.longList[k] = val << (amountToShiftHalf * 2); //double shift is needed for 2bit per base	
		}
		else
		{
			curNLongs.longList[k] = val;
		}
	}

	return curNLongs;
}

//This created int sized code for smaller segments used for multi-hash inexact homolgoy table creation
//For subSegmentLengths shorter than 8, it puts dummy A-T pairs at the center for place holder that is satisfying reverse complement property
unsigned int GetSegmentCode(char* arr, int runMode, int readLen) //runMode determines which segment to return (This version returns 3 types of 2x8bp segments from outerEdge )
{
	unsigned int val = 0;

	if(readLen >= 16 * NUM_MAX_SEEDS)
	{
		for(int i=0; i<readLen; i++)
		{
			if(SegmentTemplate[i] == runMode)
			{
				val = val << 2;
				val |= GetValueForNuc(arr[i]);
			}
		}
	}
	else
	{
		for(int i=0; i<readLen/2; i++)
		{
			if(SegmentTemplate[i] == runMode)
			{
				val = val << 2;
				val |= GetValueForNuc(arr[i]);
			}
		}
	
		for(int k=0; k<8-subSegmentLength; k++) //this fills the code with dummy A-T pairs
		{
			val = val << 2;
			val |= GetValueForNuc('A');
			val = val << 2;
			val |= GetValueForNuc('T');
		}
			
		for(int i=readLen/2; i<readLen; i++)
		{
			if(SegmentTemplate[i] == runMode)
			{
				val = val << 2;
				val |= GetValueForNuc(arr[i]);
			}
		}
	}

	return val;
}

string DebugGetStringFromSegmentCode(unsigned int code)
{
	string str("XXXXXXXXXXXXXXXX");
	for(int i=15; i>=0; i--)
	{
		str[i] = GetNucFromCode(code % 4);
		code = code >> 2;
	}
	return str;
}

void ShiftNLongLeft(nLongs& curNLong, char shiftCount)
{
	if(shiftCount == 0)
	{
		return;
	}

	//First handle the case of shifting larger than long long size
	unsigned char numCompleteShifts = shiftCount / 64;
	unsigned char effectiveNumLongs = numberOfLongs - numCompleteShifts;
	for(unsigned char k=0; k < effectiveNumLongs; k++)
	{
		curNLong.longList[k] = curNLong.longList[k+numCompleteShifts];
	}
	for(unsigned char k=effectiveNumLongs; k<numberOfLongs; k++)
	{
		curNLong.longList[k] = 0;
	}

	//Now that very large shifts are handled, do the small scale shifting

	unsigned char effectiveShiftCount = shiftCount % 64;

	for(unsigned char k=0; k<effectiveNumLongs-1; k++)
	{
		curNLong.longList[k] = curNLong.longList[k] << effectiveShiftCount;
		curNLong.longList[k] += curNLong.longList[k+1] >> (64 - effectiveShiftCount);
	}
	curNLong.longList[effectiveNumLongs-1] = curNLong.longList[effectiveNumLongs-1] << effectiveShiftCount; 
}

void RollOverNLong(nLongs& curNLong, char charToRollAdd, int readLen) //add the new char to the long long representation
{

	unsigned long long newIncrement = GetValueForNuc(charToRollAdd);
	char amountToShift = (numberOfLongs*32-readLen)*2; //*2 comes from the fact that each base reprents two bits thus shift is twice as much as the required offset
	newIncrement = newIncrement << amountToShift;

	ShiftNLongLeft(curNLong, 2); //This takes care of one base rolling
	curNLong.longList[(readLen-1)/32] += newIncrement; //This is the last block that needs to be modified
}

unsigned short ReverseComplementBitsShort(unsigned short in)
{
	unsigned short out = 0;
	for(int i=0; i<8; i++)
	{
		out <<= 2;
		out += 3-(in%4);
		in >>= 2;
	}
	return out;
}
unsigned int ReverseComplementInt(unsigned int in)
{
	unsigned int revComp = (unsigned int) shortRevCompArr[(unsigned short) in];
	in = in >> 16;
	revComp = revComp << 16;
	revComp += (unsigned int) shortRevCompArr[(unsigned short) in];
	return revComp;
}

void InitShortRevCompArray()
{
	for(int i = 0; i < USHRT_MAX; i++)
	{
		unsigned short revShort = ReverseComplementBitsShort((unsigned short) i);
		shortRevCompArr[(unsigned short) i] = revShort;		
	}  
}

unsigned long long ReverseCompLongLong(unsigned long long num)
{
	unsigned long long revComp = (unsigned long long) shortRevCompArr[(unsigned short) num];
	num = num >> 16;
	revComp = revComp << 16;
	
	revComp += (unsigned long long) shortRevCompArr[(unsigned short) num];
	num = num >> 16;
	revComp = revComp << 16;
	
	revComp += (unsigned long long) shortRevCompArr[(unsigned short) num];
	num = num >> 16;
	revComp = revComp << 16;

	revComp += (unsigned long long) shortRevCompArr[(unsigned short) num];

	return revComp;
}

nLongs ReverseComplementNLong(nLongs curNLong, int readLen)
{
	unsigned char numberOfEffectiveLongs = (readLen-1) / 32 + 1; //don't waste time for empty blocks
	nLongs reverse;
	for(unsigned char k=0; k<numberOfEffectiveLongs; k++)
	{
		reverse.longList[k] = ReverseCompLongLong(curNLong.longList[numberOfEffectiveLongs-k-1]);
	}
	ShiftNLongLeft(reverse, 2*(numberOfEffectiveLongs*32-readLen)); //2* is because each base in a read represents 2 bits
	return reverse;
}

//Skips the current window to a readmer position without N
int GetFirstNonNreadmer(char* chrRef, int pos, int readLen, int refLen) 
{
	int startPos = pos;
	int curPos = startPos;
	while(true)
	{
		if(chrRef[curPos] == 'N')
		{
			curPos++;
			startPos = curPos;
		}
		else
		{
			curPos++;
			if(curPos - startPos >= readLen)
			{
				return startPos;
			}
		}		
		if(curPos >= refLen)
		{
			return -1;
		}
	}

	assert(0);
	return -1;
}

//Checks if the first sequence is lexicographically smaller than the second
int IsLexiSmall_NLongs(const nLongs& first, const nLongs& second)
{
	for(unsigned char k=0; k<numberOfLongs; k++)
	{
		if(first.longList[k] < second.longList[k])
		{
			return 1;
		}
		if(first.longList[k] > second.longList[k])
		{
			return 0;
		}
	}

	//two of them are identical 
	return -1;
}

struct MapSegment
{
	char chrCode; //chromosome ID
	int chrPos; //regardless of the direction of which segment of the read-mer, chrPos shows the first base of the forward strand full read-mer
	bool dir; //direction of the segment w.r.t forward ref
	MapSegment* next; //next segment in the same bin
};

//Generated the template for which seeds (segments) will be sampled from the read-mers
void InitSegmentTemplates(int readLen, int numMismatches)
{
	if(numMismatches < 1 || numMismatches > 3)
	{
		cout << "ERROR: Number of mismatches for inexact homology needs to be between 1 and 3" << endl;
		exit(6);
	}

	//Sample case for 50bp" 12312312312312312312312300321321321321321321321321

	string segmentTemplateStr(readLen, '0');

	int numSeeds = numMismatches + 1;

	//generalized code for multiple mismatches
	if(readLen >= numSeeds * 16)
	{
		for(int i=0; i<numSeeds * 8; i++)
		{
			segmentTemplateStr[i] = '1' + (i%numSeeds);
			segmentTemplateStr[readLen-i-1] = '1' + (i%numSeeds);
		}

		for(int i= numSeeds * 8; i<readLen-numSeeds*8; i++)
		{
			segmentTemplateStr[i] = '0';
		}

		subSegmentLength = 8;
	}
	else
	{
		int coreHalfLen = readLen/2;
		int halfEmptyLen = (readLen % (2*numSeeds)) / 2;	

		for(int i=0; i<coreHalfLen - halfEmptyLen; i++)
		{
			segmentTemplateStr[i] = '1' + (i%numSeeds);
			segmentTemplateStr[readLen-i-1] = '1' + (i%numSeeds);
		}

		subSegmentLength = readLen/(numSeeds * 2); //lower-rounding handles overhead
	}	

	cout << "segmentTemplateStr: " << segmentTemplateStr << endl;

	for(int i=0; i<readLen; i++)
	{	
		SegmentTemplate[i] = segmentTemplateStr[i] - '0';
	}
	
}

//These are the offsets of the parallelization signals used for each segment
void SetSegmentSignalOffsets(int& so, int& eo, int& t_so_1, int& t_eo_1, int& t_so_2, int& t_eo_2, int readLen, int runMode) //returned offsets are 0-based
{
	int count = 0;
	for(int i=0; i<readLen; i++)
	{
		if(SegmentTemplate[i] == runMode)
		{
			if(count == 0)
			{
				so = i;
				count++;
			}
			else if(count == 1)
			{
				t_so_1 = i;
				count++;
			}
			else if(count == 2)
			{
				t_so_2 = i;
				break;
			}
		}	
	}

	count = 0;
	for(int i=readLen-1; i > t_so_2; i--)
	{
		if(SegmentTemplate[i] == runMode) 
		{
			if(count == 0)
			{
				eo = i;
				count++;
			}
			else if(count == 1)
			{
				t_eo_1 = i;
				count++;
			}
			else
			{
				t_eo_2 = i;
				break;
			}
		}
	}
}

//Buffer for printing homologies (in separate threads that don't lock each other)
homPair outBuffer[MAX_THREADS+2][BUFFER_ITEMS+2];
int bufferSize[MAX_THREADS+2];
int bufferLimit[MAX_THREADS+2];

void InitBuffers()
{
	for(int i=1; i<=MAX_THREADS; i++)
	{
		bufferLimit[i] = BUFFER_ITEMS;
	}
}

void FlushBuffer(homPair buff[], int &buffSize, FILE* foutFom)
{
	fwrite(buff, sizeof(homPair), buffSize, foutFom);
	buffSize = 0;
}

//Look at identical segments and report if the read-mers containing them are homologous
void CompareAndReport(MapSegment* seg1, MapSegment* seg2, FILE* foutHomMap, int runMode, int readLen, int threadNo, int errorLimit, bool palindromeFLAG) //runMode will affect which segment from outer edge is matched
{
	//runMode n has nth segment identical => will report any homology unless no segment before n has a perfect match
	assert(runMode >= 1 && runMode <= errorLimit + 1);

	int chrCode1 = seg1->chrCode, chrCode2 = seg2->chrCode;
	int chrPos1 = seg1->chrPos, chrPos2 = seg2->chrPos;

	char* ref1 = fullRef[chrCode1] + chrPos1;
	char* ref2 = fullRef[chrCode2] + chrPos2;

	int numErrors = 0;
        int errorPos[3]; //this is bounded by the maximum number of mismatches
	errorPos[0] = 0; errorPos[1] = 0, errorPos[2] = 0;

	unsigned char segmentErrorCounts[5]; //There could be at most numMismatches + 1 segments
	segmentErrorCounts[1] = 0; segmentErrorCounts[2] = 0; segmentErrorCounts[3] = 0, segmentErrorCounts[4] = 0;

	if(seg1->dir == seg2->dir)
	{
		for(int i=0; i<readLen; i++)
		{
			if(SegmentTemplate[i] != runMode && ref1[i] != ref2[i]) //empty ones need to be checked
			{
				segmentErrorCounts[(short) SegmentTemplate[i]]++; //increment specific error counts for each segment
				errorPos[numErrors] = i + 1;
				numErrors++;
                                if(numErrors > errorLimit)
                                {
					return;
                                }
			}
		}
	}
	else //check in revComp direction
	{
		for(int i=0; i<readLen; i++)
                {
                        if(SegmentTemplate[i] != runMode && ref1[i] != RevCompChar[ref2[readLen-i-1]]) //positions in the template that are different than segment if need to be checked
			{
				segmentErrorCounts[(short) SegmentTemplate[i]]++;
                                errorPos[numErrors] = i + 1;
                                numErrors++;
                                if(numErrors > errorLimit)
                                {
                                        return;
                                }
                        }
		}	
	}

	if(numErrors==0)
	{
		//Create an exception for the case that recoverse complement palindrome is inserted with perfect homology to itself	
		if(palindromeFLAG == 1 && chrCode1 == chrCode2 && chrPos1 == chrPos2)
		{
			return;
		}
		else
		{
			//Otherwise it should have been caught during equivalence class construction
			cout << "ERROR: no substitutions (such a case should have been caught earlier" << endl;
			assert(numErrors > 0); //Otherwise it should have been caught during duplicate collapsing
		}
	}	

	if(segmentErrorCounts[runMode] != 0)
	{
		cout << "ERROR: SegmentErrorCounts for runmode itself cannot be nonzero: " << (int) segmentErrorCounts[runMode] << endl;
		cout << "chrPos1: " << chrPos1 << " chrPos2: " << chrPos2 << " runMode: " << (int) runMode << endl;
	}
	assert(segmentErrorCounts[runMode] == 0);
	assert(errorPos[0] != -1);
	if(numErrors == 2)
	{
		assert(errorPos[1] != -1);
	}
	if(numErrors == 3)
	{
		assert(errorPos[0] != -1);
	}


	//The statement below forces that homology is reported for the current segment (runMode) only if there is no earlier segment that is identically matching	
	switch(runMode)
	{
		case 4:
			if(segmentErrorCounts[3] == 0) return;
		case 3:
			if(segmentErrorCounts[2] == 0) return;
		case 2:
			if(segmentErrorCounts[1] == 0) return;
		default:
			break;
	}

	homPair& outMap = outBuffer[threadNo][bufferSize[threadNo]];
	
	outMap.chrCode1 = chrCode1;
	outMap.chrCode2 = chrCode2;
	outMap.chrPos1 = chrPos1;
	outMap.chrPos2 = chrPos2;
	outMap.dir = !(seg1->dir == seg2->dir);

	// the offset direction will always be in phase with the positive genome direction of the first readmer (corrections are done while reading)
	outMap.err1off = errorPos[0];
	outMap.err2off = 0;
	outMap.err3off = 0;

	if(numErrors >= 2)
	{
		outMap.err2off = errorPos[1]; //name is 1-based but index is 0-based
	}

	if(numErrors >= 3)
	{
		outMap.err3off = errorPos[2];
	}

	assert(outMap.err1off <= readLength && outMap.err2off <= readLength && outMap.err3off <= readLength);

	bufferSize[threadNo]++;
	
	if(bufferSize[threadNo] >= bufferLimit[threadNo])
	{
		FlushBuffer(outBuffer[threadNo], bufferSize[threadNo], foutHomMap);
	}
}

//Recursively clean-up segments
int DestroyMapSegment(MapSegment *item)
{
	if(item)
	{
		int x = DestroyMapSegment(item->next);
		delete item;
		return x + 1;
	}
	return 0;
}

//Destroy segments for the entire hash table
void Clean_Up_Hash(std::tr1::unordered_map<unsigned int, MapSegment*> hash)
{
	for(std::tr1::unordered_map<unsigned int, MapSegment*>::iterator it = hash.begin(); it != hash.end(); ++it)
	{
		MapSegment* curNode = it->second;
		DestroyMapSegment(curNode);
	}
	hash.clear();
}

//2-base or 4-base split signals for inexact homology table parallelization
char splitSignal_4b_hash[25][25][25][25];
char splitSignal_2b_hash[25][25];

//Put proper split signal groups depending on the input...
int FillInSplitSignals(string splitSignalList[], const string& splitSignalListStr)
{
	stringstream SS1(splitSignalListStr);
	int count = 0;
	string junk;
	while(SS1 >> junk)
	{
		count++;
	}

	stringstream SS2(splitSignalListStr);
	count = 0;
	string signal;
	while(SS2 >> signal)
	{
		splitSignalList[count] = signal;
		count++;
	
		string revComp = ReverseComplementString(signal);
		if(revComp != signal)
		{
			splitSignalList[count] = revComp;
			count++;
		}
	}

	assert(count > 0);
	
	return count;
}

//Parse the auxiliary 4-base splitting file
int Read4charSignals(ifstream& fin4charSignal, string splitSignalList[100][150], int signalListCounts[])
{
	string line;
	getline(fin4charSignal, line);
	
	int numLines = atoi(line.c_str());
	for(int i=1; i<=numLines; i++)
	{
		getline(fin4charSignal, line);
		int count = FillInSplitSignals(splitSignalList[i], line);
		signalListCounts[i] = count;
	
		for(int k=0; k<count; k++)
		{
			string str_k = splitSignalList[i][k];
			splitSignal_4b_hash[str_k[0]-'A'][str_k[1]-'A'][str_k[2]-'A'][str_k[3]-'A'] = i; 
		}
	}
	return numLines;
}

//Verify if the characters seen match the 2-base signals
bool DoesSignalMatch_2char_faster(char firstChar, char secChar, int iterNo)
{
        return (splitSignal_2b_hash[firstChar-'A'][secChar-'A'] == iterNo);
}

//Verify if the characters seen match the 4-base signals
bool DoesSignalMatch_4char_faster(char firstChar, char fourthChar, char secChar, char thirdChar, int threadNo)
{
	return (splitSignal_4b_hash[firstChar-'A'][secChar-'A'][thirdChar-'A'][fourthChar-'A'] == threadNo);
}

void AssertRange(int val, int rangeStart, int rangeEnd)
{
	assert(val >= rangeStart && val <= rangeEnd);
}

int main(int argc, char* argv[])
{
	double startTime = getTime();
	if(argc!=10)
	{
		cout << "Given a reference fasta file, this program constructs exact and inexact homology tables in not fully compressed form" << endl;
		cout << "ARGV[1] INPUT: refFile" << endl;
		cout << "ARGV[2] INPUT: readLen [33-64bp]" << endl; //Ask for this separately since the read-mer length could be independent from fastq file
		cout << "ARGV[3] INPUT: 2-char signal list file" << endl; //Only string that will be hashed are the ones that begin/end with the two-char split signal (or revComp) 1-char on ech side
		cout << "ARGV[4] OUTPUT: equalMaps encoded file" << endl;
		cout << "ARGV[5] OUTPUT: rep-to-rep E2 homology file" << endl;
		cout << "ARGV[6] INPUT: 4char splits with number of splits on top" << endl;
		cout << "ARGV[7] MODE: EXACT, BOTH (at this point inexact can't be computed by itself)" << endl;
		cout << "ARGV[8] PARAM: Number of mismatches for satisfying homology" << endl;
		cout << "ARGV[9] PARAM: Partial Run Mode. Options: FULL, ONLY_CONSTRUCT, ONLY_COMPACT, ONLY_EXACT_COMPACT, ONLY_INEXACT_COMPACT" << endl;
		exit(19);
	}

	SetupRevCompChar();

	string runMode(argv[9]);
	readLength = atoi(argv[2]);
	if(readLength < 33 || readLength > 64)
	{
		cout << "This homology table constructor works for 33-64bp reads" << endl;
		cout << "An update is required for other readLengths" << endl;
		exit(65);
	}
	
	string OutputExactHomTable(argv[4]);
	string OutputExactHomTable_PreCompact = OutputExactHomTable + "_preCompact_toDelete";

	int numSignals; //First line of the splits file determines how many splits there are going to be 
	ifstream finSignal(argv[3]); //input file with a list of 2-char signals to divide the genome with
	if(!finSignal.is_open())
	{
		cout << "ERROR: Can't open file: " << argv[3] << endl;
		exit(66);
	}

	string numSignalS;
	getline(finSignal, numSignalS);
		
	numSignals = atoi(numSignalS.c_str());

	string splitSignalList[500][100];
	int signalListCounts[500];
	for(int i=0; i<500; i++)
	{
		signalListCounts[i] = 0;
	}

	ifstream fin4charSignal(argv[6]);
	if(!fin4charSignal.is_open())
	{
		cout << "ERROR: can't open file: " << argv[6] << endl;
		exit(67);
	}

	string splitSignalList_4char[100][150];
	int signalListCounts_4char[100];	
	int num4charSignals = Read4charSignals(fin4charSignal, splitSignalList_4char, signalListCounts_4char);

	//Create multiple output files (so that hom map printing doesn't get locked)
	string OutputInexactHomTable(argv[5]);
	string OutputInexactHomTable_PreCompact = OutputInexactHomTable + "_preCompact_toDelete";	

	int numMismatch = atoi(argv[8]);

	string finRefName(argv[1]);
	ConfigureRelatedParameters(finRefName);

	long long numEquivClasses = -1; //This is used in exact comapct however, the size of it can be optimized during the construct stage

	if(runMode == "FULL" || runMode == "ONLY_CONSTRUCT")
	{
		ReadReference(finRefName);

		InitShortRevCompArray(); //This fills in the lookup table for quick reverse complement of shorts
		InitBuffers(); //Manual buffers for printing output faster

		unsigned int equalMapIndex = 0; //This gives a uniq equivalence class ID to any group of reads that map to the same bin (indices are carried over from previous signal splits);
		
		//Output for equivalence classes (exact homologies)
		FILE* foutEqualMaps = fopen(OutputExactHomTable_PreCompact.c_str(), "w");
		if(foutEqualMaps == NULL)
		{
			cout << "ERROR: Can't write to file: " << OutputExactHomTable_PreCompact << endl;
			exit(68);
		}

		int signalNo = 0;

		cout << "Finding identicals...\t" << getTime() - startTime << endl;

		int numExactPrecompactLinesPrinted = 0;
		string splitSignalListStr;
		//Construct the exact homology table in splits determined by the input (splitting here is not for parallelization but for memory management)
		while(getline(finSignal, splitSignalListStr)) //These are physical splits -- so no single-machine parallelization here
		{
			signalNo++;
			
			signalListCounts[signalNo] = FillInSplitSignals(splitSignalList[signalNo], splitSignalListStr);		

			for(int k=0; k<signalListCounts[signalNo]; k++)
			{
				string str_k = splitSignalList[signalNo][k];
				splitSignal_2b_hash[str_k[0]-'A'][str_k[1]-'A'] = signalNo;
			}

			fprintf(foutEqualMaps,"signals %s\n",splitSignalListStr.c_str());
			numExactPrecompactLinesPrinted++;

			cout << "Signal Level: " << signalNo << "\t" << splitSignalListStr << "\t" << getTime() - startTime << endl;
			
			std::tr1::unordered_map<nLongs, long long, nLongsHash, nLongsPred> collapsed;
			float z = collapsed.max_load_factor();
			collapsed.max_load_factor ( z / 1.13 );
			collapsed.rehash(200000000);

			for(int chrNo = 1; chrNo <= numChrs; chrNo++)
			{
				fprintf(foutEqualMaps,"chrCode %d\n", chrNo);
				numExactPrecompactLinesPrinted++;
				long long prevPrintedEqualMapNo = 0, prevPrintedChrPos = 0;

				char* chrRef = fullRef[chrNo];
				char* chrRefMarks = fullRefMarks[chrNo];		
					
				int startPos = GetFirstNonNreadmer(chrRef, 1, readLength, chrLens[chrNo]); //1 is the position to start searching
				int charToRollAdd = 'N'; // N is assigned as a dummy to fail the assert check if there is a logical error
				bool rollOverFLAG = 0;

				nLongs curNLong;
				for(int curPos = startPos; ; curPos++)
				{
					if(rollOverFLAG == 0)
					{
						curNLong = GetNLongsFromString(chrRef+curPos,readLength);
						rollOverFLAG = 1;
					}
					else
					{
						//roll from previous
						assert(charToRollAdd != 'N');
						RollOverNLong(curNLong, charToRollAdd, readLength);
					}

					if(DoesSignalMatch_2char_faster(chrRef[curPos], chrRef[curPos+readLength-1], signalNo))
					{
						nLongs revCompNLong = ReverseComplementNLong(curNLong, readLength);
						int strandSwitchedFLAG = 0;
						long long *ptr;

						int dirCheck = IsLexiSmall_NLongs(curNLong, revCompNLong); //returns -1 if it's a palindrome
		
						if(dirCheck) //still find this one for palindrome
						{
							ptr = &(collapsed[curNLong]); //collapse forwards
						}
						else
						{
							strandSwitchedFLAG = 1;
							ptr = &(collapsed[revCompNLong]); //collapse backwards
						}

						if(chrRefMarks[curPos] != 0)
						{
							cout << "ERROR 9: Corruption when filling chrRefMarks" << endl;
							exit(0);
						}

						if((*ptr) == 0)
						{
							equalMapIndex++;
							*ptr = equalMapIndex;
							chrRefMarks[curPos] = signalNo; // This line marks the positions in the reference as the representative for the equivalence class (the id of the signal is noted)
						}

						long long curEqualMapNo = *ptr;
						if(curEqualMapNo - prevPrintedEqualMapNo != 1)
						{
							fprintf(foutEqualMaps, "E%lld\t", curEqualMapNo - prevPrintedEqualMapNo);
						}
						fprintf(foutEqualMaps, "%d", strandSwitchedFLAG);
			
						if(curPos - prevPrintedChrPos != 1)
						{
							fprintf(foutEqualMaps, "\t%lld", curPos - prevPrintedChrPos);
						}
						fprintf(foutEqualMaps,"\n");
						numExactPrecompactLinesPrinted++;
						prevPrintedEqualMapNo = curEqualMapNo;
						prevPrintedChrPos = curPos;

						//Process palindrome in addition [Check if this addition of palindromes cause touble downstream]
						if(dirCheck == -1)
						{
							//Forward is already printed, so print the reverse
					
							if(curEqualMapNo - prevPrintedEqualMapNo != 1)
							{
								fprintf(foutEqualMaps, "E%lld\t", curEqualMapNo - prevPrintedEqualMapNo);
							}
							strandSwitchedFLAG = 1;
							fprintf(foutEqualMaps, "%d", strandSwitchedFLAG);

							if(curPos - prevPrintedChrPos != 1)
							{
								fprintf(foutEqualMaps, "\t%lld", curPos - prevPrintedChrPos);
							}
							fprintf(foutEqualMaps,"\n");
							numExactPrecompactLinesPrinted++;
							prevPrintedEqualMapNo = curEqualMapNo;
							prevPrintedChrPos = curPos;
						}
					}

					char newChar= chrRef[curPos+readLength];
					if(newChar == '\0')
					{
						break;
					}	

					if(newChar == 'N')
					{
						curPos = GetFirstNonNreadmer(chrRef, curPos+readLength+1, readLength, chrLens[chrNo]) - 1; //-1 is since it will be incremented in the next loop;
						if(curPos < 0) //Any negative result from the function means that there is no such read-mer left
						{
							break;
						}
						//newChar wont be rolled - reset curLong read
						rollOverFLAG = 0;
					}

					charToRollAdd = newChar;					
				}
			}
			collapsed.clear();
		}
		fclose(foutEqualMaps);

		//If the user only wants to construct exact homology don't go any further
		if(string(argv[7]) == "EXACT")
		{
			//There is nothing more to process (exit)
			return 0;
		}

		numEquivClasses = equalMapIndex;
		int numSegments = numMismatch + 1;

		// Segments are partitioned based on the template provided in the following function (One requirement is that the template should be palindromic - reverse complement will fall into the same segmentation)
		InitSegmentTemplates(readLength, numMismatch); //This fills in the shapes of the 3-way split segments

		assert(numSignals == signalNo);

		finSignal.close();
		
		string homOutputFileNames[100];
		FILE* foutHomMaps[100];
		for(int i=1; i<=num4charSignals; i++)
		{
			stringstream nameSS;
			nameSS << OutputInexactHomTable_PreCompact << "_t" << i; 
			homOutputFileNames[i] = nameSS.str();
			foutHomMaps[i] = fopen(homOutputFileNames[i].c_str(),"wb");
		}	

		//This for loop constructs all inexact homologies (there are 2 layers of splitting here -one is for memory management -- reducing the total size in memory at a given time--- and the other for parallelization
		for(int signalNo = 1; signalNo <= numSignals; signalNo++)
		{
			cout << endl << "SplitNo: " << signalNo << "\t" << getTime() - startTime << endl;
			cout << "Creating Hash Tables\t" << getTime() - startTime << endl;
			
			#pragma omp parallel num_threads(num4charSignals) //Parallelize depending of 4char signals
			{
				int threadSignalNo = omp_get_thread_num() + 1; //returned thread no is originally 0, but here I increment by one to match with signal numbers		
				for(int segmentNo = 1; segmentNo <=numSegments; segmentNo++) //[TODO] Changing this hard-coded 3 here would allow changing 2 error limit
				{
					std::tr1::unordered_map<unsigned int, MapSegment*> collapsed_segment; //This holds the hash table for the current segment		

					float z = collapsed_segment.max_load_factor();
					collapsed_segment.max_load_factor ( z / 1.13 );
					collapsed_segment.rehash(2000000);
				
					int segmentStartSignalOffset = -1, segmentEndSignalOffset = -1;
					int segmentThreadStartSignalOffset_1 = -1, segmentThreadEndSignalOffset_1 = -1; //Thread offsets are one step inside (but still parallel)
					int segmentThreadStartSignalOffset_2 = -1, segmentThreadEndSignalOffset_2 = -1; //second thread offsets are one further inside
				
					SetSegmentSignalOffsets(segmentStartSignalOffset, segmentEndSignalOffset, segmentThreadStartSignalOffset_1, segmentThreadEndSignalOffset_1,
								segmentThreadStartSignalOffset_2, segmentThreadEndSignalOffset_2, readLength, segmentNo); //These return 0-based
				
					for(int chrNo = 1; chrNo <= numChrs; chrNo++)
					{
						char* chrRef = fullRef[chrNo];
						char* chrRefMarks = fullRefMarks[chrNo];
						int lenLimit = chrLens[chrNo]-readLength+1;
						for(int curPos = 1; curPos<=lenLimit; curPos++) //go through all position in the reference
						{
							if(chrRefMarks[curPos]) //marked position -- determined by the exact homology construction -- if marked, it means that it is the representative of its equivalence class
							{
								if(chrRefMarks[curPos] < 0 || chrRefMarks[curPos] > numSignals)
								{
									cout << "ERROR 8: corruption in chrRefMarks" << endl;
									exit(0);
								}
					
								assert(chrRefMarks[curPos] <= numSignals);
								
								//Normal physical split signals
								char s_Start = chrRef[curPos+segmentStartSignalOffset];
								char s_End = chrRef[curPos+segmentEndSignalOffset];					

								//Split signals for threads (4char)
								char sT_Start_1 = chrRef[curPos+segmentThreadStartSignalOffset_1];
								char sT_End_1 = chrRef[curPos+segmentThreadEndSignalOffset_1];
								char sT_Start_2 = chrRef[curPos+segmentThreadStartSignalOffset_2];
								char sT_End_2 = chrRef[curPos+segmentThreadEndSignalOffset_2];
								
								if(DoesSignalMatch_2char_faster(s_Start, s_End, signalNo))
								if(DoesSignalMatch_4char_faster(sT_Start_1, sT_End_1, sT_Start_2, sT_End_2, threadSignalNo))
								{ //This looks at proper segment positions for physical splits					
									unsigned int fwSegmentCode = GetSegmentCode(chrRef+curPos, segmentNo, readLength);
									unsigned int rcSegmentCode = ReverseComplementInt(fwSegmentCode);
						
									MapSegment **ptr;

									bool strandSwitchedFLAG = 0;
									bool palindromeFLAG = 0;			
			
									if(fwSegmentCode <= rcSegmentCode)
									{
										ptr = &(collapsed_segment[fwSegmentCode]); //hashing is done in this line
										if(fwSegmentCode == rcSegmentCode)
										{
											palindromeFLAG = 1;
										}								
									}
									else
									{
										ptr = &(collapsed_segment[rcSegmentCode]); //and the reverse complement here
										strandSwitchedFLAG = 1;
									}

									MapSegment* curMapSegment = new(MapSegment);
									
									curMapSegment->chrCode = chrNo;
									curMapSegment->chrPos = curPos;
									curMapSegment->dir = strandSwitchedFLAG;
									curMapSegment->next = NULL;

									if((*ptr) != 0)
									{
										//Search and compare to every other Map Segment
										MapSegment* otherMapSegment = (*ptr);
										while(otherMapSegment != NULL)
										{
											CompareAndReport(curMapSegment, otherMapSegment, foutHomMaps[threadSignalNo], segmentNo, readLength, threadSignalNo, numMismatch, 0); 
											otherMapSegment = otherMapSegment->next;
										}
									}
									
									//Either way insert MapSegment to the List
									curMapSegment->next = (*ptr);
									(*ptr) = curMapSegment;	
									
									//If the segment is a palindrome, try the other way too (since segments are symmetric both directions would be testable)
									if(palindromeFLAG)
									{
										MapSegment* palindromeMapSegment = new(MapSegment);
										palindromeMapSegment->chrCode = chrNo;
										palindromeMapSegment->chrPos = curPos;
										palindromeMapSegment->dir = 1;
										palindromeMapSegment->next = NULL;

										assert((*ptr) != 0);

										assert(strandSwitchedFLAG == 0);
								
										MapSegment* otherMapSegment = (*ptr);

										while(otherMapSegment != NULL)
										{
											CompareAndReport(palindromeMapSegment, otherMapSegment, foutHomMaps[threadSignalNo], segmentNo, readLength, threadSignalNo, numMismatch, 1); 
											otherMapSegment = otherMapSegment->next;
										}
										palindromeMapSegment->dir = strandSwitchedFLAG;
										//Insert palindrome as well to the list
										palindromeMapSegment->next = (*ptr);
										(*ptr) = palindromeMapSegment;
									}
								}	
							}
						}
					}
					Clean_Up_Hash(collapsed_segment);
				}
			}					
		}

		//This is related to buffer code if will be used later for binary outputs
		for(int i=1; i<=num4charSignals; i++)
		{
			if(bufferSize[i] > 0)
			{
				FlushBuffer(outBuffer[i], bufferSize[i], foutHomMaps[i]);
				fclose(foutHomMaps[i]);
			}
		}

		cout << "Finished finding inexact homologies... " << endl;

		//Free up all the memory used for construction at this point in the code
		assert(fullRefMarks[0] == 0);
		for(int i=0; i<=numChrs; i++)
		{
			free(fullRef[i]);
			free(fullRefMarks[i]);
		}
	}

	//Compaction calls
	if(runMode == "FULL" || runMode == "ONLY_COMPACT" || runMode == "ONLY_EXACT_COMPACT")
	{
		if(runMode == "ONLY_COMPACT" || runMode == "ONLY_EXACT_COMPACT")
		{
			numEquivClasses =  GetTotalGenomeSize();
		}

		assert(numEquivClasses != -1);		

		main_exactCompact(readLength, OutputExactHomTable_PreCompact, OutputExactHomTable, numEquivClasses);
	}
	
	if(runMode == "FULL" || runMode == "ONLY_COMPACT" || runMode == "ONLY_INEXACT_COMPACT")
	{
		if(DEBUG_SWITCH_FOR_DETAILED_PRINT == 1 && runMode == "ONLY_INEXACT_COMPACT")
		{
			cout << "ERROR: Detailed Debug Mode is not compatible with ONLY_INEXACT_COMPACT run mode." << endl;
			cout << "SOLUTION: Use ONLY_COMPACT instead that also runs exact compact" << endl;
			exit(23);
		}

		int numInexactCompactSplits = numChrs;

		// This is the maximum contiguous sequence to compact at a given time for memory handling (sub-chromosomal split)
		int maxContigLen = 0; //0 means there is no minimum contig len 		
		if(numMismatch == 3)
		{
			maxContigLen = 50000000; //For E3 it should be smaller for proper memory handling
		}

		cout << "maxContigLen: " << maxContigLen << endl;

		main_inexactCompact(num4charSignals, numInexactCompactSplits + 2, readLength, OutputInexactHomTable_PreCompact + "_t", OutputInexactHomTable, numMismatch, maxContigLen);
	}

	return 0;
}
