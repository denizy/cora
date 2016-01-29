/*Copyright (c) 2015-2016 Deniz Yorukoglu. All rights reserved.*/

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<assert.h>
#include<stdlib.h>
#include<string>
#include<string.h>
#include<tr1/unordered_map>
#include<time.h>
#include<sys/time.h>
#include<unistd.h>
#include<signal.h>
#include<vector>
#include<limits.h>
using namespace std;

//[FEATURES]: Takes FASTQ input and collapses (no quality values are considered)
//[DONE]: Preserve Identities of each single read while collapsing
//[DONE]: Integrate open closed form printing
//[DONE]: Long read name separation should be done here as well (otherwise it takes time to re-read and print the file).
//[DONE]: Switch to multiple input files (Containg a saparate read fastq file for each individual)

#define MAX_NUM_SAMPLES 65535 //MAX LIMIT ON NUMBER OF SAMPLES is 65535 (Due to MAX_unsigned_short limit) (SAMPLE = number of fastq files in the input dataset)
#define MAX_LINE_LEN 300 //Maximum number of bases in a given line in a read dataset as well as the reference (used for static allocation of strings)
#define N_COUNT_THRESHOLD 2 //Defunct -- Previously used for when a read should be disregarded when it has more than a certain number of unknown bases ('N')
#define MRSFAST_LRN_THRESHOLD 194 //this is the maximum read length that Mrsfast can safely handle (thus anything longer is separated into a different list)
#define BWA_LRN_THRESHOLD 5000000 // ideally 2 billion for BWA but linkConstruct buffer causes trouble //2 billion iis the maximum read length that BWA can safely handle (thus anything longer is separated into a different list)
#define GENERIC_LRN_THRESHOLD 194 //[TODO] Test for Bowtie2 and separate their upper limit as well
#define MAX_ID_DIGIT_LEN 10 //This is the maximum length of a compact read ID (though the actual reasonable limit is ~7 which equals to 192 billion reads!)

////////////////////////////////////////////////
//Different flags used for collapsing
////////////////////////////////////////////////
#define PAIRED_MODE 22 //Input dataset is paired end
#define SINGLE_MODE 11 //Input dataset is single end

#define NO_SPLIT_MODE 33 //Input dataset will be collapsed as whole reads
#define HALF_SPLIT_MODE 44 //Input dataset will be collapsed as half-split reads
#define THREEWAY_SPLIT_MODE 55 //Input dataset will be collapsed as threeway-split reads

#define BWA_MODE 17 //This requires special handling of read names ending with /1 and /2
#define BOWTIE_2_MODE 18 //This requires printing a faux *.fastw files since bowtie2 causes trouble with repeated '>' in read name but not repeated '@'
#define BWA_MEM_MODE 19 //This requires special handling of names ending with /0 through /9 ([TODO] a better fix will be needed for this since there will be a lot of special cases)
#define GENERIC_MODE 34 //Generic collapsing scheme with default CORA read name assignments in FASTA format
////////////////////////////////////////////////

unsigned char inputMode; //stores paired end or single end mode flag
unsigned char splitMode; //stores no split or half split mode flag

inline double getTime() //This is the function that gets timestamps
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return t.tv_sec+t.tv_usec/1000000.0;
}

char idDigitLen; //This stores how many characters each read identity will contain

struct hashItem //A hash item contains locus information + linked list of samples and counts having the read 
{
	char refCode;
	unsigned int refPos;
	char refDir; //This is to store which direction the reference readmer is stored

	//TODO(denizy) For speed improvement, for any name that uses less than a pointer size (store it directly inside) -- which saves a pointer lookup cost
	//TODO(denizy) For both space and speed improvement, store all of these in a large static array - and point to the index (though a proper memory manager has to be added for growing names)
	char* ident; //list of identities within hashed item (representation form is a simple character array) 
	 
	unsigned int identCount; //as identCount grows, the capacity grows as idDigitLen *  0, 1, 2, 4, 8, 16 etc.
};

char RevCompChar[256];
void SetupRevCompChar()
{
        RevCompChar['A'] = 'T';
        RevCompChar['T'] = 'A';
        RevCompChar['C'] = 'G';
        RevCompChar['G'] = 'C';
        RevCompChar['N'] = 'N';
}

//Check if the reverse complement of the given read-mer is lexicographically smaller than the original (if so, store it in revComp array)
inline bool CheckIfComplementIsEarlier(char* refStr, char revCompStr[], int readLen) //this is a slow way to check whether the reverse complement of the seqence is lexicographically smaller than itself (might make this faster)
{
	for(int i=0; i<readLen; i++)
	{
		char compChar = RevCompChar[(int) refStr[readLen-1-i]];
		if(refStr[i] < compChar)
		{
			return 0;
		}
		else 
		{
			revCompStr[i] = compChar;
			if(refStr[i] > compChar)
			{
				for(int k=i+1; k<readLen; k++)
				{
					revCompStr[k] = RevCompChar[(int) refStr[readLen-1-k]];
				}
				revCompStr[readLen] = '\0';
				return 1;
			}
		}
	}
	revCompStr[readLen] = '\0';
	return 0;
}

//Assigns the reverse complements for half-splits of a read -- the same idea as CheckIfComplementIsEarlier is used to prevent excessive base comparisons
inline void AssignHalfSplitComplementsAndDirections(char* readLine, char* compReadLine, int fullReadLength, bool& strandSwitchedFLAG_firstSplit, bool& strandSwitchedFLAG_secondSplit)
{
	int halfLen = fullReadLength / 2;

	char* forwFirstSplit = readLine;
	char* compFirstSplit = compReadLine + halfLen;

	char* forwSecondSplit = readLine + halfLen;
	char* compSecondSplit = compReadLine;

	for(int i=0; i<halfLen; i++)
	{
		compFirstSplit[i] = RevCompChar[(int) forwFirstSplit[halfLen-1-i]];
		if(forwFirstSplit[i] < compFirstSplit[i])
		{
			strandSwitchedFLAG_firstSplit = 0;
			break;
		}	
		else if(forwFirstSplit[i] > compFirstSplit[i])
		{
			strandSwitchedFLAG_firstSplit = 1;
			for(int k=i+1; k<halfLen; k++)
			{
				compFirstSplit[k] = RevCompChar[(int) forwFirstSplit[halfLen-1-k]];
			}
			compFirstSplit[halfLen] = '\0';
			break;
		}
	}

	for(int i=0; i<halfLen; i++)
	{
		compSecondSplit[i] = RevCompChar[(int) forwSecondSplit[halfLen-1-i]];
		if(forwSecondSplit[i] < compSecondSplit[i])
		{
			strandSwitchedFLAG_secondSplit = 0;
			break;
		}
		else if(forwSecondSplit[i] > compSecondSplit[i])
		{
			strandSwitchedFLAG_secondSplit = 1;
			for(int k=i+1; k<halfLen; k++)
			{
				compSecondSplit[k] = RevCompChar[(int) forwSecondSplit[halfLen-1-k]];
			}
			break;
		}
	}
}

//Assigns the reverse complements for third-splits of a read -- the same idea as CheckIfComplementIsEarlier is used to prevent excessive base comparisons
inline void AssignThirdSplitComplementsAndDirections(char* readLine, char* compReadLine, int fullReadLength, bool& strandSwitchedFLAG_firstSplit, bool& strandSwitchedFLAG_secondSplit, bool& strandSwitchedFLAG_thirdSplit)
{
	int thirdLen = fullReadLength / 3;

	char* forwFirstSplit = readLine;
	char* compFirstSplit = compReadLine + thirdLen + thirdLen;

	char* forwSecondSplit = readLine + thirdLen;
	char* compSecondSplit = compReadLine + thirdLen;

	char* forwThirdSplit = readLine + thirdLen + thirdLen;
	char* compThirdSplit = compReadLine;

	for(int i=0; i<thirdLen; i++)
	{
		compFirstSplit[i] = RevCompChar[(int) forwFirstSplit[thirdLen-1-i]];
		if(forwFirstSplit[i] < compFirstSplit[i])
		{
			strandSwitchedFLAG_firstSplit = 0;
			break;
		}	
		else if(forwFirstSplit[i] > compFirstSplit[i])
		{
			strandSwitchedFLAG_firstSplit = 1;
			for(int k=i+1; k<thirdLen; k++)
			{
				compFirstSplit[k] = RevCompChar[(int) forwFirstSplit[thirdLen-1-k]];
			}
			compFirstSplit[thirdLen] = '\0';
			break;
		}
	}

	for(int i=0; i<thirdLen; i++)
	{
		compSecondSplit[i] = RevCompChar[(int) forwSecondSplit[thirdLen-1-i]];
		if(forwSecondSplit[i] < compSecondSplit[i])
		{
			strandSwitchedFLAG_secondSplit = 0;
			break;
		}
		else if(forwSecondSplit[i] > compSecondSplit[i])
		{
			strandSwitchedFLAG_secondSplit = 1;
			for(int k=i+1; k<thirdLen; k++)
			{
				compSecondSplit[k] = RevCompChar[(int) forwSecondSplit[thirdLen-1-k]];
			}
			break;
		}
	}

	for(int i=0; i<thirdLen; i++)
	{
		compThirdSplit[i] = RevCompChar[(int) forwThirdSplit[thirdLen-1-i]];
		if(forwThirdSplit[i] < compThirdSplit[i])
		{
			strandSwitchedFLAG_thirdSplit = 0;
			break;
		}
		else if(forwThirdSplit[i] > compThirdSplit[i])
		{
			strandSwitchedFLAG_thirdSplit = 1;
			for(int k=i+1; k<thirdLen; k++)
			{
				compThirdSplit[k] = RevCompChar[(int) forwThirdSplit[thirdLen-1-k]];
			}
			break;
		}
	}
}

//checks if the current value is a power of 2 -- this check is for resizing without spending extra memory for storing capacity
inline bool IsLog2Integer(unsigned int val) 
{
	while(val % 2 == 0)
	{
		val/=2;
	}

	if(val == 1) //It divided till it became 1, meaning that it's a power of two
	{
		return true;
	}
	
	return false; //if it falls out of the loop it's not a two-power	
}

unsigned int globalLRNCount; //This is the counter for LRN codes, in the LRN list (LRN = Long read Name = Encoded names that are too long for the off-the-shelf mapper used) 	

//Alphabet size of the read name encoding should be text readable otherwise off-the-shelf mappers can't read it (it can be made binary by modifying mappers' source)
#define IDENTITY_ALPHABET_START 33 //This is the starting character of the text readable alphabet "!"
#define IDENTITY_ALPHABET_END   126 //This is the ending character of the text readdable alphabet "~"
#define IDENTITY_ALPHABET_SIZE  94  //This is the full size of the alphabet per character

unsigned char curID[MAX_ID_DIGIT_LEN+2]; //Stores the encoded read ID for the currently processed read or split
//ID is stored as: direction for lowest bit (0 forward, 1 is reverse); if paired-end, mate id for the next lowest bit (0 is first mate, 1 is second); 
//	if half-split mode the split id for the next lowest bit (0 is first split, 1 is second split); remainder is for the order in which the read appears in the original dataset

//Quick incrementation of the read ID
inline void IncrementID(unsigned char prevDirection, unsigned char strandSwitchedFlag) //This incrementation works for both single and paired end modes
{
	unsigned char incrAmount = (2 - prevDirection) + strandSwitchedFlag;

	unsigned char curPos = idDigitLen - 1;
	curID[curPos] += incrAmount;

	while(curPos > 0 && curID[curPos] > IDENTITY_ALPHABET_END) //if curPos becomes 0 and satisfies property then it means that the idDigitLen was not properly created
	{
		curID[curPos] -= IDENTITY_ALPHABET_SIZE;
		curID[curPos -1] ++;
		curPos--;
	}
}

//Prints read id in integer and character form 
void DebugPrintID()
{
	cout << "int: ";
	for(int i=0; i<idDigitLen; i++)
	{
		cout << int(curID[i]) << " "; 
	}
	cout << "\tchar: '";
	for(int i=0; i<idDigitLen; i++)
	{
		cout << curID[i];
	}
	cout << "'" <<endl;
}

int repMapMode; //Which off-the-shelf mapper is going to used for repMap = representative mapping = coarse mapping -- depending on this the collaper will output different types of output 

char globalLRNString[50]; //this stores the encoding of the current LRN name 
	//(#1#, #1143# etc -- the trick is that the length of the LRN is not divisible by idDigitLen (so that mapInfer can understand that it is not the real name but the LRN encoded ID)

inline void AssignGlobalLRN() //This assigns global LRN while making sure that length of the string is not divisible by idDigitLen (for differentiating between LRN and encoded chars)
{
	unsigned char  n = sprintf (globalLRNString+1, "%d", globalLRNCount);
	globalLRNString[n+1] = '#';
	
	if((n+2) % idDigitLen != 0) //Just keep it normal if not dividisble by idDigitLen
	{
		globalLRNString[n+2] = '\0';
	}
	else
	{
		globalLRNString[n+2] = '#';
		globalLRNString[n+3] = '\0';
	}
}

char dummyQualString[256]; //Stores "S" for the entire quality string, used since of Bowtie2 can process fasta files, so a fastq file with dummy quality strings is printed for coarse mapping

//TODO(denizy) Increase this to allow for larger number of chromosomes for mapping
#define MAX_NUM_CHRS 256 //Limit for number of chromosomes in the input -- If longer is needed convert hashItems character refID's to short refID's

char* fullRef[MAX_NUM_CHRS+2]; //char array array storing the reference bases (not compact -- 1 byte per base) -- Positions are 1-based ( fullRef[0] is empty )
unsigned int chrLens[MAX_NUM_CHRS+2]; //Each chromosome's length in fullRef
unsigned short numChrs; //Number of chromosomes in reference

//This function loads reference genome for a given .fa file with a corresponding .fa.fai file
void ReadReference(const string& refFile, int refLineLen)
{
        string refFileFai = refFile + ".fai";
        //Verifying refFileFai & reference size
        ifstream finFai(refFileFai.c_str());

        int chrCode = 0;
        string faiLine;
        while(getline(finFai, faiLine))
        {
                stringstream faiLineSS(faiLine);

                string junk;
                int ctgSize;
                faiLineSS >> junk >> ctgSize;

                chrCode++;
                fullRef[chrCode] = (char *) calloc (ctgSize+MAX_LINE_LEN, sizeof(char));
                chrLens[chrCode] = ctgSize;
        }
        numChrs = chrCode;


        //Dummy space allocation for the 0th chromosome (will not be used)
        fullRef[0] = (char *) malloc (MAX_LINE_LEN);

        FILE* finRef = fopen(refFile.c_str(),"r");

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
        for(unsigned short chr=1; chr<=numChrs; chr++)
        {
                char* curPtr = fullRef[chr];
 	        curPtr[chrLens[chr] + 1] = '\0';    // Puts end of string markers for each chromosome, for easier while-loop search
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

////////////////////////////////////////////
// Functions related ro fragmentation
// 	- For large input datasets that won't fit in the memory for collapsing, is done in fragmented batches.
//	- Splitting hash table into fragments is done by looking at 1 or 2 bases at each end of the read-mers
///////////////////////////////////////////////

char frags[1000][10]; //each of them is a fragmentation signal (AC, TT, etc.)
int numFrags = 0; //total number of fragmentation groups (each group might have multiple signals)
bool nFLAG; //When this is turned on this frag run will also handle strings that start or end with 'N';

void SetupFragmentation_2(string fragModeStr) //fragmentation done by one character from each end (allows splitting into 10 groups)
{
	if(fragModeStr == "NONE")	
	{
		numFrags = 0;
	}
	else
	{
		if(fragModeStr[0]=='N')
		{
			nFLAG = 1;
			fragModeStr = fragModeStr.substr(2, fragModeStr.length()-2);
		}

		assert(fragModeStr.length() % 3 == 0);
		numFrags = fragModeStr.length() / 3;
		for(int i=0; i<numFrags; i++) //All reverse complements are already provided
		{
			frags[i][0] = fragModeStr[3*i];
			frags[i][1] = fragModeStr[3*i +1];
		}
	}
}

void SetupFragmentation_4(string fragModeStr) //fragmentation done by two characters from each end (allows splitting into ~144 groups)
{
	if(fragModeStr == "NONE")	
	{
		numFrags = 0;
	}
	else
	{
		if(fragModeStr[0]=='N')
		{
			nFLAG = 1;
			fragModeStr = fragModeStr.substr(2, fragModeStr.length()-2);
		}

		assert(fragModeStr.length() % 5 == 0);
		numFrags = fragModeStr.length() / 5;
		for(int i=0; i<numFrags; i++) //All reverse complements are already provided
		{
			frags[i][0] = fragModeStr[5*i];
			frags[i][1] = fragModeStr[5*i +1];
			frags[i][2] = fragModeStr[5*i +2];
			frags[i][3] = fragModeStr[5*i +3];
		}
	}
}

//Checks if the current 2-base signal is within the fragment group
bool IsFragIncluded_2(char c1, char c2)
{
	for(int i=0; i<numFrags; i++)
	{
		if(frags[i][0] == c1 && frags[i][1] == c2)
		{
			return 1;
		}
	}
	return 0;
}

//Checks if current 4-base signal is within the fragment group
bool IsFragIncluded_4(char c1, char c2, char c3, char c4)
{
	for(int i=0; i<numFrags; i++)
	{
		if(frags[i][0] == c1 && frags[i][1] == c2 && frags[i][2] == c3 && frags[i][3] == c4)
		{
			return 1;
		}
	}
	return 0;
}
////////////////////////////////////////////
//Create ID from given sequence oter and the direction 
void CreateID(unsigned char* curID, bool dir, unsigned long long seqCode) //This directionality flag now represents which direction the read was compressed against the reference
{
	if(dir)
	{
		seqCode++;
	}

	int ind = idDigitLen - 1;

	while(ind >= 0)
	{
		curID[ind] = (seqCode % IDENTITY_ALPHABET_SIZE) + IDENTITY_ALPHABET_START;
		seqCode /= IDENTITY_ALPHABET_SIZE;
		ind--;
	}
}

int main(int argc, char* argv[])
{
	double runStartTime = getTime();

	if(argc!=17)
	{
		cout << "This program collapses a given set of read sequences (or half splits) and print a compact fasta or fastq file for an off-the-shelf mapper to perform coarse mapping" << endl;
		cout << "It also captures perfect matches to the reference and seperates them for automatic construction of perfect links (coarse mapping will only produce non-perfect links)" << endl;
		cout << "ARGV[1] INPUT:  FASTQ filelist (or the binary file itself if run in Fragmented mode)" << endl; //with read counts as the third field
		cout << "ARGV[2] PARAM: Digit len for representing reads" << endl;
		cout << "ARGV[3] OUTPUT: output simplified read file in non-fasta format (Long read names will be separated with matching ids in ARGV[9]" << endl;
		cout << "ARGV[4] INPUT: reference file input (or NOREF)" << endl;
		cout << "ARGV[5] PARAM: Length per line in Ref (not important if NOREF)" << endl;
		cout << "ARGV[6] PARAM: read length" << endl;
		cout << "ARGV[7] OUTPUT: perfect mapping list output" << endl;
		cout << "ARGV[8] CURRENTLY UNUSED PARAM: Max Chr Size" << endl;
		cout << "ARGV[9] OUTPUT: Long read name list" << endl;
		cout << "ARGV[10] PARAM: Input mode" << endl; //SINGLE or PAIRED
		cout << "ARGv[11] PARAM: Split mode" << endl; //NONE or HALF
		cout << "ARGV[12] PARAM: Representative mapping mode" << endl; //This is needed to BWA integration
		cout << "ARGV[13] PARAM: NONE (no fragmentation) or AA_AC_AG_[as many as the frag signals to be considered in this stage]" << endl; //All the signals for current fragmentation signal group
		cout << "ARGV[14] PARAM: fragmentation signal length" << endl; //2 or 4
		cout << "ARGV[15] OUTPUT: temporary kill signal output" << endl; //This is done due to heavy use of unordered hash tables which takes a while for gcc to deallocate (but OS takes care of it much faster)
		cout << "ARGV[16] PARAM: LRN offset: LRNs will be counted starting from this value" << endl; //Needed if fragmentation is used
		exit(10);
	}

	SetupRevCompChar();

	/////////////////////////////////////////
	//Set parameters from input -- see beginning of code for descriptions
	/////////////////////////////////////////
	
	globalLRNCount = atoi(argv[16]); //this will offset the counts of LRN items
	string fragModeStr = string(argv[13]);

	int fragSignalLen = -1;

	if(fragModeStr != "NONE") //Setup fragmentation signals
	{
		fragSignalLen = atoi(argv[14]);
		if(fragSignalLen == 2)
		{
			SetupFragmentation_2(fragModeStr);
		}
		else if(fragSignalLen == 4)
		{
			SetupFragmentation_4(fragModeStr);
		}
		else
		{
			cout << "fragSignalLen can only be 2 or 4" << endl;
			exit(19);
		}
	}

	//TODO(denizy) Make this flexible for odd length reads (by allowing 1 based overlap or gap at the center when selecting splits)
	unsigned short readLength = atoi(argv[6]); //Length of read-mers ( equal to original read length if No_split mode, equal to half if half-split mode)

	string representativeMappingMode(argv[12]);
	if(representativeMappingMode == "BWA")
	{
		repMapMode = BWA_MODE;
	}
	else if(representativeMappingMode == "BWA_MEM")
	{
		repMapMode = BWA_MEM_MODE;
	}
	else if(representativeMappingMode == "BOWTIE_2")
	{
		repMapMode = BOWTIE_2_MODE;
		for(unsigned short i=0; i<readLength; i++) //Bowtie2 needs fastq output 
		{
			dummyQualString[i] = 'S';
		}
		dummyQualString[readLength] = '\0';
	}
	else
	{
		repMapMode = GENERIC_MODE;
	}

	globalLRNString[0] = '#'; //This is initialization of LRN strings

	string inputModeStr = argv[10];
	if(inputModeStr == "PAIRED")
	{
		inputMode = PAIRED_MODE;
	}
	else if(inputModeStr == "SINGLE")
	{
		inputMode = SINGLE_MODE;
	}
	else
	{
		cout << "ERROR: Input mode unrecognized" << endl;
		exit(8);
	}

	string splitModeStr = argv[11];
	if(splitModeStr == "FULL")
	{
		splitMode = NO_SPLIT_MODE;
	}
	else if(splitModeStr == "HALF")
	{
		splitMode = HALF_SPLIT_MODE;
	}
	else if(splitModeStr == "THREEWAY")
	{
		splitMode = THREEWAY_SPLIT_MODE;
	}
	else
	{
		cout << "ERROR: Split mode unrecognized" << endl;
		exit(9);
	}
	
	idDigitLen = atoi(argv[2]);

	cout << "Total digit length of all read identities is: " << (int) idDigitLen << endl;

	int refLineLen = atoi(argv[5]);
	string finRefName = argv[4];
	//////////////////////////////////////////////////////


	//unordered resizable hash table for read-mer sequences -- optimize load factor and rehash size that determines the starting bucket count
	std::tr1::unordered_map<string, hashItem*> collapsed;
	float z = collapsed.max_load_factor(); 
	collapsed.max_load_factor ( z / 1.13 );
	collapsed.rehash(200000000);
	
	cout << "NUMFRAGS:" << numFrags << endl;
	//This is for the reference hashing
	if(finRefName != "NOREF") //If there is reference is used for collapsing, extract unique k-mers from the reference and add them to the hash table
	{
		char compRefChr[100];
			
		cout << "Reading reference data.." << endl;
		ReadReference(argv[4], refLineLen); //Dumps the entire reference sequence to fullRef[][]

		cout << "Reference read: " << getTime() - runStartTime << endl;

		if(numFrags == 0) //hash all read-mers
		{
			for(unsigned short k=1; k<=numChrs; k++) //go through all chromosomes
			{
				unsigned int curChrLim = chrLens[k] - readLength + 1;	
				unsigned short curNcount = 0;	

				for(unsigned short i=1; i<=readLength; i++) //Counts N's within the first readMer in chromosome (remaining read-mers would be counted dynamically)
				{
					if(fullRef[k][i] == 'N')
					{
						curNcount++;
					}
				}
			
				for(unsigned int i=1; i<=curChrLim; i++) //go through all bases in the reference
				{
					if(curNcount == 0)
					{
						char saveVal = fullRef[k][i+readLength]; //This sets the end of the read-mer to \0 so that hashing function can now it ends there
						fullRef[k][i+readLength] = '\0';
					
						bool strandSwitchedFLAG = CheckIfComplementIsEarlier(fullRef[k] + i, compRefChr, readLength);
		
						hashItem **ptr;
						if(strandSwitchedFLAG == 0) //if forward direction is lexicographicallysmaller, hash the original read-mer
							ptr = &(collapsed[fullRef[k] + i]);
						else
							ptr = &(collapsed[compRefChr]); //otherwise hash the reverse complement

						if((*ptr) == 0) //If read-mer is not hashed before assign chrNo, chrPos, dir (of lexicographically small version), identity info is not used yet (those are for reads)
						{
							//Check revComp
							*(ptr) = new (hashItem); 
							(*ptr)->refCode = k;
							(*ptr)->refPos = i;
							(*ptr)->refDir = (char) strandSwitchedFLAG;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						} //Otherwise k-mer is not unique, don't do anything

						fullRef[k][i+readLength] = saveVal; //restore the character at the end of the read-mer
					}

					//This automatically adjusts the nCount in the current window
					if(fullRef[k][i] == 'N')
					{
						curNcount--;
					}
					if(fullRef[k][i+readLength] == 'N')
					{
						curNcount++;
					}
				}
			}
		}
		else //hash only read-mers fitting the frag signals
		{
			if(fragSignalLen == 2)
			{
				for(unsigned short k=1; k<=numChrs; k++)
				{
					unsigned int curChrLim = chrLens[k] - readLength + 1;	
					unsigned short curNcount = 0;	

					for(unsigned short i=1; i<=readLength; i++)
					{
						if(fullRef[k][i] == 'N')
						{
							curNcount++;
						}
					}
				
					for(unsigned int i=1; i<=curChrLim; i++)
					{
						if(curNcount == 0)
						{
							if(IsFragIncluded_2(fullRef[k][i], fullRef[k][i+readLength-1])) //Check if the 2-base frag signal is within the current frag group (if so conitnue adding the reference read-mer to the hash-table)
							{
								char saveVal = fullRef[k][i+readLength];
								fullRef[k][i+readLength] = '\0';

								bool strandSwitchedFLAG = CheckIfComplementIsEarlier(fullRef[k] + i, compRefChr, readLength);
				
								hashItem **ptr;
								if(strandSwitchedFLAG == 0)
									ptr = &(collapsed[fullRef[k] + i]);
								else
									ptr = &(collapsed[compRefChr]);
								if((*ptr) == 0)
								{
									//Check revComp
									*(ptr) = new (hashItem);
									(*ptr)->refCode = k;
									(*ptr)->refPos = i;
									(*ptr)->refDir = (char) strandSwitchedFLAG;
									(*ptr)->ident = NULL;
									(*ptr)->identCount = 0;
								}

								fullRef[k][i+readLength] = saveVal;
							}
						}

						//This automatically adjusts the nCount in the current window
						if(fullRef[k][i] == 'N')
						{
							curNcount--;
						}
						if(fullRef[k][i+readLength] == 'N')
						{
							curNcount++;
						}
					}
				}
			}
			else if(fragSignalLen == 4)
			{
				for(unsigned short k=1; k<=numChrs; k++)
				{
					unsigned int curChrLim = chrLens[k] - readLength + 1;	
					unsigned short curNcount = 0;	

					for(unsigned short i=1; i<=readLength; i++)
					{
						if(fullRef[k][i] == 'N')
						{
							curNcount++;
						}
					}
				
					for(unsigned int i=1; i<=curChrLim; i++)
					{
						if(curNcount == 0)
						{
							if(IsFragIncluded_4(fullRef[k][i], fullRef[k][i+1], fullRef[k][i+readLength-2], fullRef[k][i+readLength-1])) //Check if the 4-base frag signal is within the current frag group (if so conitnue adding the reference read-mer to the hash-table)
							{
								char saveVal = fullRef[k][i+readLength];
								fullRef[k][i+readLength] = '\0';

								bool strandSwitchedFLAG = CheckIfComplementIsEarlier(fullRef[k] + i, compRefChr, readLength);
				
								hashItem **ptr;
								if(strandSwitchedFLAG == 0)
									ptr = &(collapsed[fullRef[k] + i]);
								else
									ptr = &(collapsed[compRefChr]);
								if((*ptr) == 0)
								{
									//Check revComp
									*(ptr) = new (hashItem);
									(*ptr)->refCode = k;
									(*ptr)->refPos = i;
									(*ptr)->refDir = (char) strandSwitchedFLAG;
									(*ptr)->ident = NULL;
									(*ptr)->identCount = 0;
								}

								fullRef[k][i+readLength] = saveVal;
							}
						}

						//This automatically adjusts the nCount in the current window
						if(fullRef[k][i] == 'N')
						{
							curNcount--;
						}
						if(fullRef[k][i+readLength] == 'N')
						{
							curNcount++;
						}
					}
				}
			}
			else
			{
				cout << "fragSignalLen should be either 2 or 4" << endl;
				cout << "val: " << argv[14] << endl;
				exit(17);		
			}
		}
	}

	double timeAfterFinishingReadingRef = getTime();

	cout << "Time spent for processing reading and hashing ref: " << timeAfterFinishingReadingRef - runStartTime << endl;

	cout << "Reading read data..." << endl;

	//Read length adjustment for read repository creation
	int fullReadLength = readLength;
	int halfReadLength = -1, oneThirdReadLength = -1, twoThirdsReadLength = -1;

	if(splitMode == HALF_SPLIT_MODE)
	{
		fullReadLength = readLength * 2; //If half_split mode, set the full read length to be twice read-mer size (readLength denotes read-mer length)
		halfReadLength = readLength;
	}
	else if(splitMode == THREEWAY_SPLIT_MODE)
	{
		fullReadLength = readLength * 3;
		oneThirdReadLength = readLength;
		twoThirdsReadLength = readLength * 2;
	}

	ifstream finFastqList(argv[1]); //This contains (1st line) number of sample datasets (N+1th line) The nth sample's fastq file name (two tab-delimited names if paired end) and the number of reads within one of them 

	int numSamples; //number of samples (e.g. separate fastq files for individual samples)
	finFastqList >> numSamples;
	if(numSamples > MAX_NUM_SAMPLES)
	{
		cout << "This version will crash if the number of samples is greater than " << MAX_NUM_SAMPLES << endl;
		exit(28);
	}	

	cout << "number of Samples: " << numSamples << endl;

	if(numFrags == 0) //hash all reads as normal -- reading the original fastq [TODO] use the fastq reading version implemented in fastqSplitter
	{
		char readLine[MAX_LINE_LEN+1]; //stores original full read sequence
		char compReadLine[MAX_LINE_LEN+1]; //stores the complement sequence
		char readLine_mate[MAX_LINE_LEN+1]; //stores the original sequence of the mate read (if exists)
		char compReadLine_mate[MAX_LINE_LEN+1]; //stores complement of mate
	
		readLine[fullReadLength] = '\0'; compReadLine[fullReadLength] = '\0'; readLine_mate[fullReadLength] = '\0'; compReadLine_mate[fullReadLength] = '\0';
	
		char junkLine[MAX_LINE_LEN+1]; //for the dummy line between seq and qual as well as the quality line itself

		bool prevDirection = 0;
		for(unsigned char i=0; i<idDigitLen; i++) //initialize the current read-mer id to 0
		{
			curID[i] = IDENTITY_ALPHABET_START;
		}
		curID[idDigitLen - 1] -= 2; //this is manually decremented for being incremented again for the first occurence (does away with if check for first item)

		//This is the loop that goes over different fastq files
		string curFastqFileName;
		string curFastqFileName_mate;
		int readCount_junk; //not used here in this program but is in the input

		while(finFastqList >> curFastqFileName)
		{
			
			FILE* finRead = fopen(curFastqFileName.c_str(),"r");
			FILE* finRead_mate;

			setvbuf ( finRead , NULL , _IOFBF , 16777216);

			if(finRead == NULL)
			{
				cout << "Input fastq file: " << curFastqFileName << " doesn't exist" << endl;
				exit(5);
			}

			//A lot of repetition here, but each is doing slightly different things and I don't want to add a lot if if checks per read 
			if(inputMode == PAIRED_MODE) //for paired end reads
			{
				finFastqList >> curFastqFileName_mate >> readCount_junk;
				finRead_mate = fopen(curFastqFileName_mate.c_str(), "r");

				setvbuf ( finRead_mate , NULL , _IOFBF , 16777216);

				if(finRead_mate == NULL)
				{
					cout << "Read mate input file: " << curFastqFileName_mate << " doesn't exist" << endl;
					exit(6);
				}

				while(fgets(readLine, MAX_LINE_LEN, finRead) != NULL) //read the read header
				{
					assert(fgets(readLine_mate, MAX_LINE_LEN, finRead_mate) != NULL); //read the mate header
					assert(fgets(readLine, MAX_LINE_LEN, finRead) != NULL); //read the actual sequence
					readLine[fullReadLength] = '\0';
					assert(fgets(readLine_mate, MAX_LINE_LEN, finRead_mate) != NULL); //read the mate sequence
					readLine_mate[fullReadLength] = '\0';

					//Read other lines that are not used here (read names and quality scores can be recovered at the end of mapInfer)
					assert(fgets(junkLine, MAX_LINE_LEN, finRead) != NULL);
					assert(fgets(junkLine, MAX_LINE_LEN, finRead) != NULL);
					assert(fgets(junkLine, MAX_LINE_LEN, finRead_mate) != NULL);
					assert(fgets(junkLine, MAX_LINE_LEN, finRead_mate) != NULL);
					
					if(splitMode == NO_SPLIT_MODE) //Collapse each mate as a whole
					{
						//Currently all items are collapsed regardless of the N content (since counting N's in the read is more costly than hashing them, which doesn't need to go through all characters) 
						bool strandSwitchedFLAG = CheckIfComplementIsEarlier(readLine, compReadLine, fullReadLength);
			
						hashItem **ptr;
						if(strandSwitchedFLAG == 0)
							ptr = &(collapsed[readLine]);
						else
							ptr = &(collapsed[compReadLine]);

						if((*ptr) == 0) //If the reference doesn't have this read-mer
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							//ID changes depending on direction (we don't want to computed the id from scratch here, so we increment using previous id and its direction)
							IncrementID(prevDirection, strandSwitchedFLAG); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG; //save the current direction for the next guy
						}
						else //If read-mer exists in reference
						{
							bool readStrandFlag = strandSwitchedFLAG; //then direction represents direction to the forward reference

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag = !(strandSwitchedFLAG == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag;
						}
					
						//This is for adding the new read-mer ID to the hashItem
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else if(splitMode == HALF_SPLIT_MODE) //collapse each mate in halves
					{
						bool strandSwitchedFLAG_firstSplit = 0, strandSwitchedFLAG_secondSplit = 0;
						AssignHalfSplitComplementsAndDirections(readLine, compReadLine, fullReadLength, strandSwitchedFLAG_firstSplit, strandSwitchedFLAG_secondSplit);
	
						hashItem **ptr;
						if(strandSwitchedFLAG_firstSplit == 0)
						{
							char saveChar = readLine[halfReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine[halfReadLength] = '\0';
							ptr = &(collapsed[readLine]);
							readLine[halfReadLength] = saveChar;
						}
						else
						{
							ptr = &(collapsed[compReadLine + halfReadLength]);
						}
						
						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_firstSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_firstSplit;
						}
						else
						{
							bool readStrandFlag_firstSplit = strandSwitchedFLAG_firstSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_firstSplit = !(strandSwitchedFLAG_firstSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_firstSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_firstSplit;
						}
					
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}

						(*ptr)->identCount++;


						//From here on is the second split
						//ptr is already defined above
						if(strandSwitchedFLAG_secondSplit == 0)
						{
							ptr = &(collapsed[readLine + halfReadLength]);
						}
						else
						{
							char saveChar = compReadLine[halfReadLength];
							compReadLine[halfReadLength] = '\0';
							ptr = &(collapsed[compReadLine]);
							compReadLine[halfReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_secondSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_secondSplit;

						}
						else
						{
							bool readStrandFlag_secondSplit = strandSwitchedFLAG_secondSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_secondSplit = !(strandSwitchedFLAG_secondSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_secondSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_secondSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else if(splitMode == THREEWAY_SPLIT_MODE) //collapse each mate in thirds
					{

						bool strandSwitchedFLAG_firstSplit = 0, strandSwitchedFLAG_secondSplit = 0, strandSwitchedFLAG_thirdSplit = 0;
						AssignThirdSplitComplementsAndDirections(readLine, compReadLine, fullReadLength, strandSwitchedFLAG_firstSplit, strandSwitchedFLAG_secondSplit, strandSwitchedFLAG_thirdSplit);
	
						hashItem **ptr;
						if(strandSwitchedFLAG_firstSplit == 0)
						{
							char saveChar = readLine[oneThirdReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine[oneThirdReadLength] = '\0';
							ptr = &(collapsed[readLine]);
							readLine[oneThirdReadLength] = saveChar;
						}
						else
						{
							ptr = &(collapsed[compReadLine + twoThirdsReadLength]);
						}
						
						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_firstSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_firstSplit;
						}
						else
						{
							bool readStrandFlag_firstSplit = strandSwitchedFLAG_firstSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_firstSplit = !(strandSwitchedFLAG_firstSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_firstSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_firstSplit;
						}
					
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}

						(*ptr)->identCount++;


						//From here on is the second split
						//ptr is already defined above
						if(strandSwitchedFLAG_secondSplit == 0)
						{
							char saveChar = readLine[twoThirdsReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine[twoThirdsReadLength] = '\0';
							ptr = &(collapsed[readLine + oneThirdReadLength]);
							readLine[twoThirdsReadLength] = saveChar;
						}
						else
						{
							char saveChar = compReadLine[twoThirdsReadLength];
							compReadLine[twoThirdsReadLength] = '\0';
							ptr = &(collapsed[compReadLine]);
							compReadLine[twoThirdsReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_secondSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_secondSplit;

						}
						else
						{
							bool readStrandFlag_secondSplit = strandSwitchedFLAG_secondSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_secondSplit = !(strandSwitchedFLAG_secondSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_secondSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_secondSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;

						//From here on is the third split
						//ptr is already defined above
						if(strandSwitchedFLAG_thirdSplit == 0)
						{
							ptr = &(collapsed[readLine + twoThirdsReadLength]);
						}
						else
						{
							char saveChar = compReadLine[oneThirdReadLength];
							compReadLine[oneThirdReadLength] = '\0';
							ptr = &(collapsed[compReadLine]);
							compReadLine[oneThirdReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_thirdSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_thirdSplit;

						}
						else
						{
							bool readStrandFlag_thirdSplit = strandSwitchedFLAG_thirdSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_thirdSplit = !(strandSwitchedFLAG_thirdSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_thirdSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_thirdSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else
					{
						assert(0);
					}

					//Do everything so far for the second mate					
					if(splitMode == NO_SPLIT_MODE)
					{
						bool strandSwitchedFLAG = CheckIfComplementIsEarlier(readLine_mate, compReadLine_mate, fullReadLength);
						hashItem **ptr;
						if(strandSwitchedFLAG == 0)
							ptr = &(collapsed[readLine_mate]);
						else
							ptr = &(collapsed[compReadLine_mate]);

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						
							IncrementID(prevDirection, strandSwitchedFLAG); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG;
						}
						else
						{
							bool readStrandFlag = strandSwitchedFLAG;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag = !(strandSwitchedFLAG == (bool) (*ptr)->refDir);
							}				
							IncrementID(prevDirection, readStrandFlag); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag;
						}
					
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identitiy to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else if(splitMode == HALF_SPLIT_MODE)
					{
						//this is for the first split
						bool strandSwitchedFLAG_firstSplit = 0, strandSwitchedFLAG_secondSplit = 0;
						AssignHalfSplitComplementsAndDirections(readLine_mate, compReadLine_mate, fullReadLength, strandSwitchedFLAG_firstSplit, strandSwitchedFLAG_secondSplit);
						hashItem **ptr;
						if(strandSwitchedFLAG_firstSplit == 0)
						{
							char saveChar = readLine_mate[halfReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine_mate[halfReadLength] = '\0';
							ptr = &(collapsed[readLine_mate]);
							readLine_mate[halfReadLength] = saveChar;
						}
						else
						{
							ptr = &(collapsed[compReadLine_mate + halfReadLength]);
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						
							IncrementID(prevDirection, strandSwitchedFLAG_firstSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_firstSplit;
						}
						else
						{
							bool readStrandFlag_firstSplit = strandSwitchedFLAG_firstSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_firstSplit = !(strandSwitchedFLAG_firstSplit == (bool) (*ptr)->refDir);
							}				
							IncrementID(prevDirection, readStrandFlag_firstSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_firstSplit;
						}
					
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identitiy to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;


						//this is for the second split
						//ptr is already defined above
						if(strandSwitchedFLAG_secondSplit == 0)
						{
							ptr = &(collapsed[readLine_mate + halfReadLength]);
						}
						else
						{
							char saveChar = compReadLine_mate[halfReadLength];
							compReadLine_mate[halfReadLength] = '\0';
							ptr = &(collapsed[compReadLine_mate]);
							compReadLine_mate[halfReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						
							IncrementID(prevDirection, strandSwitchedFLAG_secondSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
						
							prevDirection = strandSwitchedFLAG_secondSplit;
						}
						else
						{
							bool readStrandFlag_secondSplit = strandSwitchedFLAG_secondSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_secondSplit = !(strandSwitchedFLAG_secondSplit == (bool) (*ptr)->refDir);
							}				
							IncrementID(prevDirection, readStrandFlag_secondSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_secondSplit;
						}
	
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identitiy to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else if(splitMode == THREEWAY_SPLIT_MODE)
					{
						bool strandSwitchedFLAG_firstSplit = 0, strandSwitchedFLAG_secondSplit = 0, strandSwitchedFLAG_thirdSplit = 0;
						AssignThirdSplitComplementsAndDirections(readLine_mate, compReadLine_mate, fullReadLength, strandSwitchedFLAG_firstSplit, strandSwitchedFLAG_secondSplit, strandSwitchedFLAG_thirdSplit);
	
						hashItem **ptr;
						if(strandSwitchedFLAG_firstSplit == 0)
						{
							char saveChar = readLine_mate[oneThirdReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine_mate[oneThirdReadLength] = '\0';
							ptr = &(collapsed[readLine_mate]);
							readLine_mate[oneThirdReadLength] = saveChar;
						}
						else
						{
							ptr = &(collapsed[compReadLine_mate + twoThirdsReadLength]);
						}
						
						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_firstSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_firstSplit;
						}
						else
						{
							bool readStrandFlag_firstSplit = strandSwitchedFLAG_firstSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_firstSplit = !(strandSwitchedFLAG_firstSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_firstSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_firstSplit;
						}
					
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}

						(*ptr)->identCount++;


						//From here on is the second split
						//ptr is already defined above
						if(strandSwitchedFLAG_secondSplit == 0)
						{
							char saveChar = readLine_mate[twoThirdsReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine_mate[twoThirdsReadLength] = '\0';
							ptr = &(collapsed[readLine_mate + oneThirdReadLength]);
							readLine_mate[twoThirdsReadLength] = saveChar;
						}
						else
						{
							char saveChar = compReadLine_mate[twoThirdsReadLength];
							compReadLine_mate[twoThirdsReadLength] = '\0';
							ptr = &(collapsed[compReadLine_mate]);
							compReadLine_mate[twoThirdsReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_secondSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_secondSplit;

						}
						else
						{
							bool readStrandFlag_secondSplit = strandSwitchedFLAG_secondSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_secondSplit = !(strandSwitchedFLAG_secondSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_secondSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_secondSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;

						//From here on is the third split
						//ptr is already defined above
						if(strandSwitchedFLAG_thirdSplit == 0)
						{
							ptr = &(collapsed[readLine_mate + twoThirdsReadLength]);
						}
						else
						{
							char saveChar = compReadLine_mate[oneThirdReadLength];
							compReadLine_mate[oneThirdReadLength] = '\0';
							ptr = &(collapsed[compReadLine_mate]);
							compReadLine_mate[oneThirdReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_thirdSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_thirdSplit;

						}
						else
						{
							bool readStrandFlag_thirdSplit = strandSwitchedFLAG_thirdSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_thirdSplit = !(strandSwitchedFLAG_thirdSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_thirdSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_thirdSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
						

					}
					else
					{
						assert(0);
					}
				}
			}
			else if(inputMode == SINGLE_MODE) //The input read dataset is single-end (a limited version of the paired end processing -- detailed comments above)
			{
				finFastqList >> readCount_junk;

				while(fgets(readLine, MAX_LINE_LEN, finRead) != NULL)
				{
					assert(fgets(readLine, MAX_LINE_LEN, finRead) != NULL);
					readLine[fullReadLength] = '\0';

					assert(fgets(junkLine, MAX_LINE_LEN, finRead) != NULL);
					assert(fgets(junkLine, MAX_LINE_LEN, finRead) != NULL);

					if(splitMode == NO_SPLIT_MODE)
					{
						bool strandSwitchedFLAG = CheckIfComplementIsEarlier(readLine, compReadLine, fullReadLength);

						hashItem **ptr;
						if(strandSwitchedFLAG == 0)
							ptr = &(collapsed[readLine]);
						else
							ptr = &(collapsed[compReadLine]);

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						
							IncrementID(prevDirection, strandSwitchedFLAG); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG;
						}
						else
						{
							bool readStrandFlag = strandSwitchedFLAG;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{
								readStrandFlag = !(strandSwitchedFLAG == (bool) (*ptr)->refDir);
							}

							IncrementID(prevDirection, readStrandFlag); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag;
						}



						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identitiy to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;

					}
					else if(splitMode == HALF_SPLIT_MODE)
					{
						//this is the first split

						bool strandSwitchedFLAG_firstSplit = 0, strandSwitchedFLAG_secondSplit = 0;
						AssignHalfSplitComplementsAndDirections(readLine, compReadLine, fullReadLength, strandSwitchedFLAG_firstSplit, strandSwitchedFLAG_secondSplit);

						hashItem **ptr;
						if(strandSwitchedFLAG_firstSplit == 0)
						{
							char saveChar = readLine[halfReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine[halfReadLength] = '\0';
							ptr = &(collapsed[readLine]);
							readLine[halfReadLength] = saveChar;
						}
						else
						{
							ptr = &(collapsed[compReadLine + halfReadLength]);
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						
							IncrementID(prevDirection, strandSwitchedFLAG_firstSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_firstSplit;
						}
						else
						{
							bool readStrandFlag_firstSplit = strandSwitchedFLAG_firstSplit;
							
							if((*ptr)->refDir != -1) //which means that the reference is not created
							{
								readStrandFlag_firstSplit = !(strandSwitchedFLAG_firstSplit == (bool) (*ptr)->refDir);

							}

							IncrementID(prevDirection, readStrandFlag_firstSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_firstSplit;
						}

						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identitiy to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;


						//This is the second split

						//ptr is already defined above
						if(strandSwitchedFLAG_secondSplit == 0)
						{
							ptr = &(collapsed[readLine + halfReadLength]);
						}
						else
						{
							char saveChar = compReadLine[halfReadLength];
							compReadLine[halfReadLength] = '\0';
							ptr = &(collapsed[compReadLine]);
							compReadLine[halfReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
						
							IncrementID(prevDirection, strandSwitchedFLAG_secondSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_secondSplit;
						}
						else
						{
							bool readStrandFlag_secondSplit = strandSwitchedFLAG_secondSplit;
							
							if((*ptr)->refDir != -1) //which means that the reference is not created
							{
								readStrandFlag_secondSplit = !(strandSwitchedFLAG_secondSplit == (bool) (*ptr)->refDir);
							}

							IncrementID(prevDirection, readStrandFlag_secondSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_secondSplit;
						}

						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identitiy to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else if(splitMode == THREEWAY_SPLIT_MODE) //collapse each mate in thirds
					{
						bool strandSwitchedFLAG_firstSplit = 0, strandSwitchedFLAG_secondSplit = 0, strandSwitchedFLAG_thirdSplit = 0;
						AssignThirdSplitComplementsAndDirections(readLine, compReadLine, fullReadLength, strandSwitchedFLAG_firstSplit, strandSwitchedFLAG_secondSplit, strandSwitchedFLAG_thirdSplit);
	
						hashItem **ptr;
						if(strandSwitchedFLAG_firstSplit == 0)
						{
							char saveChar = readLine[oneThirdReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine[oneThirdReadLength] = '\0';
							ptr = &(collapsed[readLine]);
							readLine[oneThirdReadLength] = saveChar;
						}
						else
						{
							ptr = &(collapsed[compReadLine + twoThirdsReadLength]);
						}
						
						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_firstSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_firstSplit;
						}
						else
						{
							bool readStrandFlag_firstSplit = strandSwitchedFLAG_firstSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_firstSplit = !(strandSwitchedFLAG_firstSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_firstSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_firstSplit;
						}
					
						unsigned int curIdentCount = (*ptr)->identCount;
						unsigned int curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}

						(*ptr)->identCount++;


						//From here on is the second split
						//ptr is already defined above
						if(strandSwitchedFLAG_secondSplit == 0)
						{
							char saveChar = readLine[twoThirdsReadLength]; //this is a quick fix to make it seem like the first split and restore it again
							readLine[twoThirdsReadLength] = '\0';
							ptr = &(collapsed[readLine + oneThirdReadLength]);
							readLine[twoThirdsReadLength] = saveChar;
						}
						else
						{
							char saveChar = compReadLine[twoThirdsReadLength];
							compReadLine[twoThirdsReadLength] = '\0';
							ptr = &(collapsed[compReadLine]);
							compReadLine[twoThirdsReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_secondSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_secondSplit;

						}
						else
						{
							bool readStrandFlag_secondSplit = strandSwitchedFLAG_secondSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_secondSplit = !(strandSwitchedFLAG_secondSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_secondSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_secondSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;

						//From here on is the third split
						//ptr is already defined above
						if(strandSwitchedFLAG_thirdSplit == 0)
						{
							ptr = &(collapsed[readLine + twoThirdsReadLength]);
						}
						else
						{
							char saveChar = compReadLine[oneThirdReadLength];
							compReadLine[oneThirdReadLength] = '\0';
							ptr = &(collapsed[compReadLine]);
							compReadLine[oneThirdReadLength] = saveChar;
						}

						if((*ptr) == 0)
						{
							*(ptr) = new hashItem();
							(*ptr)->refCode = -1;
							(*ptr)->refPos = -1;
							(*ptr)->refDir = -1;
							(*ptr)->ident = NULL;
							(*ptr)->identCount = 0;
					
							IncrementID(prevDirection, strandSwitchedFLAG_thirdSplit); //Note that this directionality doesn't represent the mapping direction, it represents the direction change through collapsing
							prevDirection = strandSwitchedFLAG_thirdSplit;

						}
						else
						{
							bool readStrandFlag_thirdSplit = strandSwitchedFLAG_thirdSplit;

							if((*ptr)->refDir != -1) //which means that the reference is not created
							{ 
								readStrandFlag_thirdSplit = !(strandSwitchedFLAG_thirdSplit == (bool) (*ptr)->refDir);
							}				

							IncrementID(prevDirection, readStrandFlag_thirdSplit); //This directionality flag now represents which direction the read was compressed against the reference
							prevDirection = readStrandFlag_thirdSplit;
						}
					
						curIdentCount = (*ptr)->identCount;
						curSize = curIdentCount * idDigitLen;
						if(curIdentCount == 0)
						{
							(*ptr)->ident = (char *) malloc (idDigitLen); 
						}
						else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
						{
							char* prevIdent = (*ptr)->ident;
							char* curIdent = (char *) malloc (curSize * 2);
							for(unsigned int k=0; k<curSize; k++)
							{
								curIdent[k] = prevIdent[k];
							}
							free(prevIdent);
							(*ptr)->ident = curIdent;
						}

						//Add identity to hashItem and increment
						for(int k=0; k<idDigitLen; k++)
						{
							(*ptr)->ident[curSize + k] = curID[k];
						}
						(*ptr)->identCount++;
					}
					else
					{
						assert(0);
					}
		

				}
			}
			else
			{
				assert(0);
			}
		}				
	}
	else //fragmentation of the reads is already done before by fastqSplitter, here only difference from above is to read its binary output in [unsigned long long] + [char array (read/split length)]
	{
		char* readLine;//[MAX_LINE_LEN+1];
		char compReadLine[MAX_LINE_LEN+1];

		FILE* finBinary = fopen(argv[1],"rb"); //In the fragmented run mode, argv[1] contains the filename directly		

		if(finBinary == NULL)
		{
			cout << "Binary input file: " << argv[1] << " doesn't exist" << endl;
			exit(5);
		}

		unsigned long long seqCode;					
		unsigned short seqLen = readLength;
	
		size_t sizeofLongLong = sizeof(unsigned long long);
		size_t sizeofSeqItem = sizeofLongLong + seqLen;

		unsigned int BUFFER_ITEM_COUNT = 2000000;
		char* binaryReadBuffer = (char *) malloc (sizeofSeqItem * BUFFER_ITEM_COUNT + 2);

		unsigned int binaryReadBufferOffset = 0;
		unsigned int binaryReadBufferSize = 0;			 

		while(true)
		{
			//cout << "TRUE" << endl;
			if(binaryReadBufferOffset >= binaryReadBufferSize)
			{
				//fill in buffer
				binaryReadBufferSize = fread(binaryReadBuffer, sizeofSeqItem, BUFFER_ITEM_COUNT, finBinary);
				binaryReadBufferSize *= sizeofSeqItem;

				if(binaryReadBufferSize == 0)
				{
					break;
				}
				binaryReadBufferOffset = 0;			
			}

			//Here convert binary data to useful info
			seqCode = *((unsigned long long* ) (binaryReadBuffer + binaryReadBufferOffset));	
		
			readLine = binaryReadBuffer + binaryReadBufferOffset + sizeofLongLong;
			binaryReadBufferOffset += sizeofSeqItem;
	
			bool strandSwitchedFLAG = CheckIfComplementIsEarlier(readLine, compReadLine, seqLen);
			
			hashItem **ptr;
			if(strandSwitchedFLAG == 0)
			{
				char tempChar = readLine[seqLen];
				readLine[seqLen] = 0;
			

				//cout << binaryReadBufferOffset << "\t" << readLine << endl;

				ptr = &(collapsed[readLine]);
				readLine[seqLen] = tempChar;
			}
			else
				ptr = &(collapsed[compReadLine]);

			if((*ptr) == 0)
			{
				*(ptr) = new hashItem();
				(*ptr)->refCode = -1;
				(*ptr)->refPos = -1;
				(*ptr)->refDir = -1;
				(*ptr)->ident = NULL;
				(*ptr)->identCount = 0;
			
				CreateID(curID, strandSwitchedFLAG, seqCode); //This function will just add direction (all other info is already in seqCode)
			}
			else
			{
				bool readStrandFlag = strandSwitchedFLAG;

				if((*ptr)->refDir != -1) //which means that the reference is not created
				{ 
					readStrandFlag = !(strandSwitchedFLAG == (bool) (*ptr)->refDir);
				}				

				CreateID(curID, readStrandFlag, seqCode); //This directionality flag now represents which direction the read was compressed against the reference
			}

				
			unsigned int curIdentCount = (*ptr)->identCount;
			unsigned int curSize = curIdentCount * idDigitLen;
			if(curIdentCount == 0)
			{
				(*ptr)->ident = (char *) malloc (idDigitLen); 
			}
			else if(IsLog2Integer(curIdentCount)) //Deep copy and resize as double
			{
				char* prevIdent = (*ptr)->ident;
				char* curIdent = (char *) malloc (curSize * 2);
				for(unsigned int k=0; k<curSize; k++)
				{
					curIdent[k] = prevIdent[k];
				}
				free(prevIdent);
				(*ptr)->ident = curIdent;
			}

			//Add identity to hashItem and increment
			for(int k=0; k<idDigitLen; k++)
			{
				(*ptr)->ident[curSize + k] = curID[k];
			}
			(*ptr)->identCount++;

		}
	}

	double timeAfterFinishingReadingReadData = getTime();
	cout << "Time spent for processing reading and hashing read data: " << timeAfterFinishingReadingReadData - timeAfterFinishingReadingRef << endl;
	cout << "Printing Collapsed and Perf data..." << endl;

	int numTotal = 0;

	FILE* foutRead = fopen(argv[3],"w");
	FILE* foutPerfect = fopen(argv[7],"w");
	FILE* foutLRN = fopen(argv[9],"w");


	setvbuf(foutRead, NULL, _IOFBF, 16777216);
	setvbuf(foutPerfect, NULL, _IOFBF, 16777216);
	setvbuf(foutLRN, NULL, _IOFBF, 16777216);
	
	//this is the part we go through the hash table and print each uniq read-mer sequence with its encoded read-mer ID
	for(std::tr1::unordered_map<string, hashItem*>::iterator it = collapsed.begin(); it != collapsed.end(); ++it)
	{
		//TODO(denizy) Convert all fprintfs to fwrites in binary for efficiency.
		unsigned int curIdentCount = it->second->identCount;
		char* curIdent = it->second->ident;		

		if(curIdentCount == 0)
		{
			continue;
		}

		int readNameLength = idDigitLen * curIdentCount;
		if((it->second)->refCode == -1)
                {
			if(repMapMode == GENERIC_MODE)
			{
				if(readNameLength > GENERIC_LRN_THRESHOLD) //This is to eliminate unnecessary checks for nodeLists (may not be exact actually -- rewrite this better in the future)
				{
					globalLRNCount++;
					AssignGlobalLRN(); //This assigns global LRN while making sure that length of the string is not divisible by idDigitLen (for differentiating between LRN and encoded chars)
					
					fprintf(foutLRN, "%s\t%.*s\n", globalLRNString, idDigitLen * curIdentCount, curIdent);
					fprintf(foutRead,">%s\n", globalLRNString);
				}
				else
				{
					fprintf(foutRead,">%.*s\n", idDigitLen * curIdentCount, curIdent);
				}
				
				fputs((it->first.c_str()),foutRead);
                        	fprintf(foutRead,"\n");
			}
			else if(repMapMode == BWA_MODE)
			{
				char nextToLastChar = curIdent[idDigitLen * curIdentCount - 2];
				char lastChar = curIdent[idDigitLen * curIdentCount - 1];

				if(readNameLength > BWA_LRN_THRESHOLD || (nextToLastChar == '/' && (lastChar == '1' || lastChar == '2')))
				{
					globalLRNCount++;
					AssignGlobalLRN(); //This assigns global LRN while making sure that length of the string is not divisible by idDigitLen (for differentiating between LRN and encoded chars)
					
					fprintf(foutLRN, "%s\t%.*s\n", globalLRNString, idDigitLen * curIdentCount, curIdent);
					fprintf(foutRead,">%s\n", globalLRNString);
				}
				else
				{
					fprintf(foutRead,">%.*s\n", idDigitLen * curIdentCount, curIdent);
				}
			
				fputs((it->first.c_str()),foutRead);
                        	fprintf(foutRead,"\n");
			}
			else if(repMapMode == BWA_MEM_MODE)
			{
				char nextToLastChar = curIdent[idDigitLen * curIdentCount - 2];
				char lastChar = curIdent[idDigitLen * curIdentCount - 1];

				if(readNameLength > BWA_LRN_THRESHOLD || (nextToLastChar == '/' && (lastChar >= '0' && lastChar <= '9')))
				{
					globalLRNCount++;
					AssignGlobalLRN(); //This assigns global LRN while making sure that length of the string is not divisible by idDigitLen (for differentiating between LRN and encoded chars)
					
					fprintf(foutLRN, "%s\t%.*s\n", globalLRNString, idDigitLen * curIdentCount, curIdent);
					fprintf(foutRead,">%s\n", globalLRNString);
				}
				else
				{
					fprintf(foutRead,">%.*s\n", idDigitLen * curIdentCount, curIdent);
				}
			
				fputs((it->first.c_str()),foutRead);
                        	fprintf(foutRead,"\n");
			}
			else if(repMapMode == BOWTIE_2_MODE) //Need to print fauxq instead of fasta
			{
				if(readNameLength > GENERIC_LRN_THRESHOLD) //This is to eliminate unnecessary checks for nodeLists (may not be exact actually -- rewrite this better in the future)
				{
					globalLRNCount++;
					AssignGlobalLRN(); //This assigns global LRN while making sure that length of the string is not divisible by idDigitLen (for differentiating between LRN and encoded chars)
					
					fprintf(foutLRN, "%s\t%.*s\n", globalLRNString, idDigitLen * curIdentCount, curIdent);
					fprintf(foutRead,"@%s\n", globalLRNString);
				}
				else
				{
					fprintf(foutRead,"@%.*s\n", idDigitLen * curIdentCount, curIdent);
				}
				
				fputs((it->first.c_str()),foutRead);
                        	fprintf(foutRead,"\n+\n%s\n", dummyQualString);
			}
		}
		else
		{
			fprintf(foutPerfect,"%.*s\t", idDigitLen * curIdentCount, curIdent);
	                fprintf(foutPerfect,"%d\t%d\t%d\n",(it->second)->refCode, (it->second)->refPos, readLength); //direction info is encoded in the name
		}

		numTotal++;
	}

	double timeAfterFinishingPrintingData = getTime();
	cout << "Time spent traversing and printing hashtable: " << timeAfterFinishingPrintingData - timeAfterFinishingReadingReadData << endl;
	cout << "numTotal: " << numTotal << endl;
	cout << "FULL Run time (before FLUSH): " << timeAfterFinishingPrintingData - runStartTime << endl;

	fflush(foutRead);
	fflush(foutPerfect);
	fflush(foutLRN);

	fclose(foutRead);
	fclose(foutPerfect);
	fclose(foutLRN);

	double timeAfterFclose = getTime();
	cout << "FULL Run time: " << timeAfterFclose - runStartTime << endl;

	//cout << "killing program on purpose - for faster deallocation " << endl;
	//gcc does a worse job of deallocating the hash table, so let OS take care of it
	ofstream killout(argv[15]);
	killout << "Y" << endl;
	kill(getpid(), SIGINT); //this means intentional killing for faster deallocation

return 0;
}
