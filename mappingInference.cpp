/*Copyright (c) 2015-2016 Deniz Yorukoglu. All rights reserved.*/
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<assert.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<cmath>
#include<stdio.h>
#include<string.h>
#include<tr1/unordered_map>
#include<algorithm>
#include<vector>
#include<limits.h>
using namespace std;

//limits
#ifdef ED_THREE
	#define MAX_EDITS 3
	#define MAX_DOUBLE_EDITS 6
	#define MAX_NUM_HOM_EDITS 3
#elif ED_FOUR
	#define MAX_EDITS 4
	#define MAX_DOUBLE_EDITS 8
	#define MAX_NUM_HOM_EDITS 4
#elif ED_FIVE
	#define MAX_EDITS 5
	#define MAX_DOUBLE_EDITS 10
	#define MAX_NUM_HOM_EDITS 5
#elif ED_SIX
	#define MAX_EDITS 6
	#define MAX_DOUBLE_EDITS 12
	#define MAX_NUM_HOM_EDITS 6
#else
	#define MAX_EDITS 2
	#define MAX_DOUBLE_EDITS 4
	#define MAX_NUM_HOM_EDITS 2
#endif

#define MAX_LINE_LEN 200 
#define MAX_NUM_CHRS 250 
#define MAX_CHR_NAME_LEN 50
#define MAX_READ_LEN 250
#define MAX_RETURNED_LINKSLIST 5000
#define NUM_MAX_SAMPLES 50000
#define COLLAPSED 123
#define UNCOLLAPSED 124
#define MAX_ID_DIGIT_LEN 10
#define MAX_READ_FILE_NAME_LEN 250

//alphabet properties for read name econding
#define IDENTITY_ALPHABET_START 33
#define IDENTITY_ALPHABET_END   126
#define IDENTITY_ALPHABET_SIZE  94

//flag values (arbitrary)
#define MISSING_SPLIT_OFFSET 32 //indicates that the first split is missing (value is high so that 2*DOUBLE_ERROR cannot exceed the value)
#define DOUBLE_MISSING_SPLIT_OFFSET 64 //indicates that the second split is missing

#define FASTA 5
#define SAM 6

#define PAIRED_MODE 22
#define SINGLE_MODE 11

#define NO_SPLIT 66
#define HALF_SPLIT 77

//Mapping Modes

#define ALL_MAPPING_MODE 101

// In the case of best mapping without indels BEST is same as BEST_SENSITIVE
#define BEST_MAPPING_MODE 102
#define BEST_FAST_MAPPING_MODE 106
#define BEST_SENSITIVE_MAPPING_MODE 103
#define UNIQUE_MAPPING_MODE 104
#define STRATUM_MAPPING_MODE 105

// Distance metric encoding
#define HAMMING 444
#define LEVENSHTEIN 555
int distanceMetric = HAMMING;


int mappingMode;
int mappingQualityPrintFlag;

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

string ReverseComplementString(string str) //This is a debug function -- no need to optimize
{
	string retStr = str;

	for(int i=0; i<(int)str.length(); i++)
	{
		retStr[i] = RevCompChar[(unsigned char) str[str.length() - i -1]];
	}

	return retStr;
}

int DebugCountStringDiff(const string& str1, const string& str2)
{
	assert(str1.length() == str2.length());

	int diffCount = 0;

	for(int i = 0; i<(int)str1.length(); i++)
	{
		if(str1[i] != str2[i])
		{
			diffCount++;	
		}
	}		
	
	return diffCount;
}

int numMismatchesPerReadMer; //The actual value of the mismatch per readMer during the run (as opposed to the maximum possible value defined above)
int doubleNumMismatchesPerReadMer;
int numHomTableMismatches; //The number of mismatches used to generate the homology table (as opposed to the maximum vlaue defined above)
int globalMapCountLimit = 0; //If not zero, then all mappings will only report this many mappings
unsigned int minFragmentSize, maxFragmentSize; //size interval for the fragments

double getTime()
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return t.tv_sec+t.tv_usec/1000000.0;
}

//MD is the edit script printed in SAM format
char singleMDstr[50];
char leftMDstr[50];
char rightMDstr[50];

//fast output writing to a buffer with delimiter
inline void writeStr(char*& str, char arr[], char del)
{
	int pos = 0;
	while(arr[pos] != '\0')
	{
		str[pos] = arr[pos];
		pos++;
	}
	str[pos] = del;
	str += pos + 1;
}

//fast output writing to a buffer
inline void writeStr(char*& str, char arr[])
{
	int pos = 0;
	while(arr[pos] != '\0')
	{
		str[pos] = arr[pos];
		pos++;
	}
	str += pos;
}

//Assumes value is between 0 and 1000 exclusive
inline void writeIntToStr(char*& str, int val, char del)
{
	//assert(val > 0);
	if(val < 10)
	{
		str[0] = val + '0';
		str[1] = del;
		str += 2;
	}
	else if(val < 100)
	{
		str[0] = val / 10 + '0';
		str[1] = val % 10 + '0';
		str[2] = del;
		str += 3;
	}
	else if(val < 1000)
	{
		str[0] = (val/100) + '0';
		str[1] = (val%100)/10 + '0';
		str[2] = val % 10 + '0';
		str[3] = del;
		str += 4; 
	}
	else
	{
		cout << "ERROR: SAM flags for read length should be within (0,1000)" << endl;
		exit(96);
	}
}

short itoa10(int value, char* result) 
{
	short len = 0;
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

#define LTOA_BUFSIZE (sizeof(long) * 8 + 1)

//Modified from Robert B. Stout's ltoa function implementation
char *ltoa10(long N, char* str)
{
      register int i = 2;
      long uarg;
      char *tail, *head = str, buf[LTOA_BUFSIZE];

      tail = &buf[LTOA_BUFSIZE - 1];           /* last character position      */
      *tail-- = '\0';

      if (N < 0L)
      {
            *head++ = '-';
            uarg    = -N;
      }
      else  uarg = N;

      if (uarg)
      {
            for (i = 1; uarg; ++i)
            {
                  register ldiv_t r;

                  r       = ldiv(uarg, 10);
                  *tail-- = (char)(r.rem + ((9L < r.rem) ?
                                  ('A' - 10L) : '0'));
                  uarg    = r.quot;
            }
      }
      else  *tail-- = '0';

      memcpy(head, ++tail, i);
      return str;
}

unsigned char inputMode; //paired or single end
unsigned char splitMode; //full read-mer or half-split

//Data variables for reference
char* fullRef[MAX_NUM_CHRS+2];
unsigned int chrLens[MAX_NUM_CHRS+2];
char chrNames[MAX_NUM_CHRS+2][MAX_CHR_NAME_LEN+2];
int numChrs;
int readLen; //readMer len
int finalReadLen; //length of the actual read (twice readmer length if half-split is applied)
unsigned char idDigitLen; //length of the encoded read name (per read)

#define VALSTR_LIMIT 512
char valStr[VALSTR_LIMIT+2][5]; //Valstr is used for temporary char space for integers when printing them to buffer
void SetUpValStr()
{
	for(int i=0; i<=VALSTR_LIMIT; i++)
	{
		itoa10(i, valStr[i]);		
	}
}

//=================================
struct medit //multi-edit: stores a list of substitions with their offsets
{
	unsigned char edits[MAX_EDITS];
	unsigned char pos[MAX_EDITS];
	unsigned char numEdits;
	//Edits: 
	// A C G T -> mismatches
	// D -> deletions in the read (gap in the reference in alignment)
	// a c g t -> insertions in the read (gap in the read in alignment)
};

struct doubleMedit
{
        doubleMedit(){};
        doubleMedit(struct medit& myMedit)
        {
                numEdits = myMedit.numEdits;
                for(int i=0; i<MAX_EDITS; i++)
                {
                        edits[i] = myMedit.edits[i];
                        pos[i] = myMedit.pos[i];
                }
        }

        unsigned char edits[MAX_DOUBLE_EDITS];
        unsigned char pos[MAX_DOUBLE_EDITS];
        unsigned char numEdits;
        //Edits:
        // A C G T -> mismatches
        // D -> deletions in the read (gap in the reference in alignment)
        // a c g t -> insertions in the read (gap in the read in alignment)
};

int GetIndelLengthModifier(const medit& med) //returns #deletions - #insertions in medit, this determines whether the alignment region of the readmer in the reference will get longer or shorter
{
	int offset = 0;
	for(int i=0; i<med.numEdits; i++)
	{
		switch(med.edits[i])
		{
			case 'D':
				offset++;
				break;
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'n':
				offset--;
				break;
			default:
				break;
		}
	} 
	return offset;
}

//Checks whether two medits are equal
bool equals(const medit& m1, const medit& m2)
{
	if(m1.numEdits == m2.numEdits)
	{
		for(unsigned char i=0; i<m1.numEdits; i++)
		{
			if(m1.edits[i] != m2.edits[i] || m1.pos[i] != m2.pos[i])
			{
				return 0;
			}
		}
	}
	else
	{
		return 0;
	}
	return 1;
}
unsigned char meditEditCountLimit; //Different from EDIT_MAX_LENGTH that shows capacity, this represents a variable limit on the size of current medits (which might differ based on splits and length, but mostly will stay static)
medit dummyMedit;

//=====================================
//Data types for storing read links (in the future handle these gracefully by not loading all)

struct link //A link contains the position in the reference it points to with the direction of homology as well as the edits (medit format)
{
	link(){};
	link(unsigned int y, medit z, unsigned char x, bool t): chrPos(y), edit(z), chrCode(x), dir(t) {}

	unsigned int chrPos;	
	medit edit;
	unsigned char chrCode; //if chrCode is 0, chrPos denotes the character position of the nohit array (the exact position in noHitArray that sequence starts)
	char dir;

	bool operator() (const link& a, const link& b) //For sorting
	{
		return ( (unsigned long long) a.chrCode * UINT_MAX + a.chrPos < (unsigned long long) b.chrCode * UINT_MAX + b.chrPos);
	}
} linkObject;

link* links; //each read / mate or split gets one index, directionality doesn't get a new index, but are encoded in bool 
unsigned int linksSize; //total number of links currently stored in memory
link dummyLink;

struct mergedLink //A link contains the position in the reference it points to with the direction of homology as well as the edits (medit format)
{
        mergedLink(){};
        mergedLink(struct link& myLink)
        {
                chrPos = myLink.chrPos; chrCode = myLink.chrCode; dir = myLink.dir; edit = doubleMedit(myLink.edit);
        }
        mergedLink(unsigned int y, medit z, unsigned char x, bool t)
	{
                chrPos = y; chrCode = x; dir = t; edit = doubleMedit(z);
	}
        mergedLink(unsigned int y, doubleMedit z, unsigned char x, bool t): chrPos(y), edit(z), chrCode(x), dir(t) {}

        unsigned int chrPos;
        doubleMedit edit;
        unsigned char chrCode; //if chrCode is 0, chrPos denotes the character position of the nohit array (the exact position in noHitArray that sequence starts)
        char dir;

        bool operator() (const mergedLink& a, const mergedLink& b) //For sorting
        {
                return ( (unsigned long long) a.chrCode * UINT_MAX + a.chrPos < (unsigned long long) b.chrCode * UINT_MAX + b.chrPos);
        }
} mergedLinkObject;

bool sortByErrors(const mergedLink &lhs, const mergedLink &rhs)
{
	return lhs.edit.numEdits < rhs.edit.numEdits;
}

///////////////////////////////////
//Data types for storing no hit sequences (loaded alon the way, and array is overwritten

#define MAX_NO_HIT_ARR_SIZE  500000000 //This is required to be less than DONT_LOAD_NOHIT_FLAG
#define DONT_LOAD_NOHIT_FLAG 4294967290 //This is used for pruning when looking for nohit reads (unmapped reads during coarse mapping) 
short noHitMode;  //SAM or FA (fasta)

char* noHitArr; //All no hit reads are stored in one large array 
unsigned int noHitItemLimit; //number of nohit items in this list
unsigned int noHitSizeLimit; //number of bases in the list

//Allocate proper memory for nohit array
void SetupNoHitArr(unsigned int numTotalReads)
{
	noHitItemLimit = min(numTotalReads, (unsigned int) (MAX_NO_HIT_ARR_SIZE / readLen));	

	//The trigger for reloading nohit arr will be when when we intercept a link with chr = 0 pos = 0, then the following happens 
	// (1) scan links from the that link positions, identify which links to load noHits for
	//	- Any link with more than 1 missing splits should not be loaded (mark those by DONT_LOAD_NOHIT_FLAG)
	//	- Start from 0 in nohitarr and Reserve readLen space in order for each chr=0 link that is not marked DONT_LOAD_NOHIT_FLAG
	//	- While doing this store indices of the reserved space in the chrPos field of the link[splitId]
	//      - Stop when you can't reserve any more space stop and mark the beginning readId and ending readId (not split Id, bt the reduced readId without dir, split or mate)
	// (2) Scan the nohit file from the beginning to the end (depending on SAM or FA)
	//	- If the link[splitID].chrPos is 0 then
	// 	- Place nohit split sequence to their right place in the nohit array as specified by link[splitID].chrPos

	noHitArr = (char *) malloc ((noHitItemLimit + 2) * readLen);
	noHitSizeLimit = noHitItemLimit * readLen;
}

unsigned long long GetReadIdFromStringWithPos(const string& str, unsigned int pos); //Convert from encoded read name of idDigitLen size at offset pos in the given string 
string GetLongReadName(string code); //Extracts the actual encoded read name for long read names (LRN) 

//No hit array is loaded in batches, whenever a new nohit item is needed that is not in memory, the next batch is loaded and so on
void ReloadNoHitArray_SplitSingleEnd(unsigned int startSplitID, char unmappedSplitsFileName[]) //startSplitID is the index in the links table that we start the noHitAssignment (Guaranteed to be the id of the first split of the first mate)
{
	ifstream finUnmappedSplits(unmappedSplitsFileName);
	assert(finUnmappedSplits.is_open());

	unsigned int curPointerInNoHit = 1;
	unsigned int curSplitID = startSplitID;
	unsigned int curIndex = 0;

	while(curPointerInNoHit < noHitSizeLimit && curSplitID < linksSize)
	{
		short numEmpties = 0;
		
		if(links[curSplitID].chrPos == 0)
		{
			numEmpties++;
			curIndex = curSplitID;
			links[curSplitID].chrPos = DONT_LOAD_NOHIT_FLAG; //this is temporary until we're sure of single noHit
		}

		if(links[curSplitID+1].chrPos == 0)
                {
                        numEmpties++;
			links[curSplitID+1].chrPos = DONT_LOAD_NOHIT_FLAG; //this is temporary until we're sure of single noHit

			curIndex = curSplitID+1;
                }

                if(numEmpties != 1)
		{
			curSplitID += 2;
			continue;
		}

		//Now we will be reserving space for curIndex
//		assert(links[curIndex].chrCode == 0);

		links[curIndex].chrPos = curPointerInNoHit;
		curPointerInNoHit += readLen; 	
		curSplitID += 2;
	} 

	unsigned int endingSplitID = curSplitID;

	//now curIndex is the index of the last split that has reserved a place in noHitArray
	//Here scan through the no hit file and retrieve all sequences 

	string name, junk, seq;

	string inputLine;
	while(getline(finUnmappedSplits, inputLine))
	{
		if(noHitMode == SAM)
		{
			stringstream samLineSS(inputLine);
			samLineSS >> name >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> seq;
		}
		else //if(noHitMode == FASTA)
		{
			getline(finUnmappedSplits, seq);			
			name = inputLine.substr(1,inputLine.length()-1);			
		}

		if(name.length() % idDigitLen != 0)
		{
			name = GetLongReadName(name);
		}

		unsigned int nameLength = name.length();
		for(unsigned int i=0; i<nameLength; i += idDigitLen)
		{
			unsigned int splitID = GetReadIdFromStringWithPos(name, i);
			bool curSplitDir = splitID % 2;
			splitID /= 2; //to eliminate direction
	
			if(links[splitID].chrCode != 0 || splitID < startSplitID || splitID >= endingSplitID)
			{
				continue;
			}

			unsigned int curNoHitSeqPos = links[splitID].chrPos; //chrPos stores the index of the reserved space in nohitArr
			if(curNoHitSeqPos != DONT_LOAD_NOHIT_FLAG)
			{
				links[splitID].dir = curSplitDir;
				for(int k=0; k<readLen; k++)
				{
					noHitArr[curNoHitSeqPos + k] = seq[k];
				}
			}
		}	
	}

	finUnmappedSplits.clear();
	finUnmappedSplits.close();
}

//No hit array is loaded in batches, whenever a new nohit item is needed that is not in memory, the next batch is loaded and so on
void ReloadNoHitArray_SplitPaired(unsigned int startSplitID, char unmappedSplitsFileName[]) //startSplitID is the index in the links table that we start the noHitAssignment (Guaranteed to be the id of the first split of the first mate)
{
	ifstream finUnmappedSplits(unmappedSplitsFileName);
	assert(finUnmappedSplits.is_open());

	unsigned int curPointerInNoHit = 1;
	unsigned int curSplitID = startSplitID;
	unsigned int curIndex = 0;

	while(curPointerInNoHit < noHitSizeLimit && curSplitID < linksSize)
	{
		short numEmpties = 0;
		
		if(links[curSplitID].chrPos == 0)
		{
			numEmpties++;
			curIndex = curSplitID;
			links[curSplitID].chrPos = DONT_LOAD_NOHIT_FLAG; //this is temporary until we're sure of single noHit
		}

		if(links[curSplitID+1].chrPos == 0)
                {
                        numEmpties++;
			links[curSplitID+1].chrPos = DONT_LOAD_NOHIT_FLAG; //this is temporary until we're sure of single noHit

			curIndex = curSplitID+1;
                }

		if(links[curSplitID+2].chrPos == 0)
                {
                        numEmpties++;
			links[curSplitID+2].chrPos = DONT_LOAD_NOHIT_FLAG; //this is temporary until we're sure of single noHit

			curIndex = curSplitID+2;
                }

		if(links[curSplitID+3].chrPos == 0)
                {
                        numEmpties++;
			links[curSplitID+3].chrPos = DONT_LOAD_NOHIT_FLAG; //this is temporary until we're sure of single noHit

			curIndex = curSplitID+3;
                }

                if(numEmpties != 1)
		{
			curSplitID += 4;
			continue;
		}

		//Now we will be reserving space for curIndex
//		assert(links[curIndex].chrCode == 0);

		links[curIndex].chrPos = curPointerInNoHit;
		curPointerInNoHit += readLen; 	
		curSplitID += 4;
	} 

	unsigned int endingSplitID = curSplitID;

	//now curIndex is the index of the last split that has reserved a place in noHitArray
	//Here scan through the no hit file and retrieve all sequences 

	string name, junk, seq;

	string inputLine;
	while(getline(finUnmappedSplits, inputLine))
	{
		if(noHitMode == SAM)
		{
			stringstream samLineSS(inputLine);
			samLineSS >> name >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> junk >> seq;
		}
		else //if(noHitMode == FASTA)
		{
			getline(finUnmappedSplits, seq);			
			name = inputLine.substr(1,inputLine.length()-1);			
		}

		if(name.length() % idDigitLen != 0)
		{
			//cout << "name: " << name << endl;
			name = GetLongReadName(name);
		}

		unsigned int nameLength = name.length();
		for(unsigned int i=0; i<nameLength; i += idDigitLen)
		{
			unsigned int splitID = GetReadIdFromStringWithPos(name, i);
			bool curSplitDir = splitID % 2;
			splitID /= 2; //to eliminate direction
	
			if(links[splitID].chrCode != 0 || splitID < startSplitID || splitID >= endingSplitID)
			{
				continue;
			}

			unsigned int curNoHitSeqPos = links[splitID].chrPos; //chrPos stores the index of the reserved space in nohitArr
			if(curNoHitSeqPos != DONT_LOAD_NOHIT_FLAG)
			{
				links[splitID].dir = curSplitDir;
				for(int k=0; k<readLen; k++)
				{
					noHitArr[curNoHitSeqPos + k] = seq[k];
				}
			}
		}	
	}

	finUnmappedSplits.clear();
	finUnmappedSplits.close();
}

////////////////////////////////////
//Data types for storing inexact homology information [extended to multi edit homologies on 11/7/2013]
//Atomic unit of Inexact Homology Table -- format similar to links, but here homologies also have blockLength and edits are not explicitly specified since they are already in fullRef
struct emHomItem //multipleError homologies
{
	unsigned int chrPos; //pos in chr for the target
	unsigned char offsetList[MAX_NUM_HOM_EDITS]; //list of offset positions that represent the substitutions in homology -- offsets are always in the forward direction of the source position
	unsigned int blockLen : 7; //How long the homology is continued as a block: e.g if blockLen is 3, then A ~ B, A+1 ~ B+1, A+2 ~ B+2 for forward homologies, A~B, A+1 ~ B-1, A+2 ~ B-2 for reverseComplement homologies
	unsigned int dir : 1; //direction of homology
	unsigned char chrCode; //chromosome of the target
};

struct emHomClassNode //for each equiv class representative (described below), this stores the list of inexact homology links it has
{
	emHomItem* list; //list of inexact homologies for this position
	unsigned int chrPos; //pos in chr of source
	unsigned int listSize; //number of homologies for this pos
	unsigned char chrCode; //chr of source
	unsigned char maxRevLookup; //max number indices to go back to find all emHoms relevant to this pos (this is required since, due to blockLengths > 1, previous positions might contain homologies relevant to this guy)
};

emHomClassNode* emHomClassList; //1-based list of all equivalence class representatives with their inexact homology info
unsigned int emHomClassListSize; //total number of these

///////////////////////////////////////
//Data types for storing equivalence classes (perfect homology)

//Atomic unit of exact homology table
struct eqItem 
{
	unsigned int chrPos; //position of the target
	char dir; //direction of homology
	unsigned char chrCode; //chromosome of the target
};

//Equivalence class representative that contains exact homology links to the other members of the class
struct equivClassNode
{
	eqItem* list; //list of links other members of the class
	unsigned int blockLen; //length of block for the entire class (unlike inexact homologies for which the block length is per link), each region in the class has perfect homologie for the entire blockLength
	unsigned int chrPos; //position of class rep
	unsigned int listSize; //total number of non-rep members of the class
	unsigned char chrCode; //chromosome of the class rep
};

equivClassNode* equivClassList; //1-based list of all equivalence class representatives with their exact homologies
unsigned int equivClassListSize; //total number of these (this can be different than emHomClassListSize, since for inexact any position in an equivalence class of itself if it has at least one inexact homology, but those are not reprented in this list since they are unnecessary for euqivalence purposes)

/////////////////////////////
std::tr1::unordered_map<unsigned int, unsigned int> equivHub[MAX_NUM_CHRS+2];
std::tr1::unordered_map<unsigned int, unsigned int> homHub[MAX_NUM_CHRS+2];
////////////////////////////

int GetChrCode(const string& chrStr)
{
	for(int i=1; i<=numChrs; i++)
	{
		if(chrStr == chrNames[i])
		{
			return i;
		}
	}
	assert(0);
}

//Updated to extend homologies to multiEdit homologies on 11/7/2013
void LoadInexactHomologies(char* inexactTableFileName, int readLen, int numHomTableMismatches)
{
	cout << "Loading inexact homology maps.." << endl;

	//open itemCounts file that contains the length of the files.
	ifstream finItemCounts((string(inexactTableFileName) + ".itemCount").c_str());

	if(!finItemCounts.is_open())
	{
		cout << "Could not find inexact homology table item counts file: " << string(inexactTableFileName) + ".itemCount" << endl;
		exit(89);
	}

	finItemCounts >> emHomClassListSize;

	emHomClassList = (emHomClassNode *) calloc ((emHomClassListSize+2), sizeof(emHomClassNode)); //Format is a list of main inexact homology blocks pointing to a list of inexact homology links from that position
	int bufferListCap = 10000;
	int bufferListSize = 0;
	emHomItem* bufferList = (emHomItem *) malloc ((bufferListCap+2) * sizeof(emHomItem));	

	FILE* finBin = fopen(inexactTableFileName, "rb");	
	if(!finBin)
	{
		cout << "Could not find inexact homology table: " << inexactTableFileName << endl;
		exit(90);
	}

	for(unsigned int i=1; i<=emHomClassListSize; i++)
	{
		size_t tempSize = fread(&(emHomClassList[i].chrCode), sizeof(unsigned char), 1, finBin);
		
		tempSize &= fread(&(emHomClassList[i].chrPos), sizeof(unsigned int), 1, finBin);

		bufferListSize = 0;
	
		char curDir; //The type of these parameters are important since they are the types in the homology table
		short curBlockLen; //In order to change them, the homology table data types shoudl also be upgraded
		while(true)
		{

			tempSize &= fread(&(curDir), sizeof(unsigned char), 1, finBin);
			if(curDir == -1)
			{
				//End of list -- copy buffer to emHomClassList and quit
				emHomClassList[i].list = (emHomItem *) malloc ((bufferListSize) * sizeof(emHomItem));
				for(int k=0; k<bufferListSize; k++)
				{
					emHomClassList[i].list[k] = bufferList[k];
				}
				emHomClassList[i].listSize = bufferListSize;
				bufferListSize = 0;
				break;
			}

			bufferList[bufferListSize].dir = curDir;

			tempSize &= fread(&(bufferList[bufferListSize].chrCode), sizeof(unsigned char), 1, finBin);
			tempSize &= fread(&(bufferList[bufferListSize].chrPos), sizeof(unsigned int), 1, finBin);
		
			for(unsigned char k=0; k<numHomTableMismatches; k++)
			{
				tempSize &= fread(&(bufferList[bufferListSize].offsetList[k]), sizeof(unsigned char), 1, finBin);
			}
							
			tempSize &= fread(&curBlockLen, sizeof(short), 1, finBin);
			bufferList[bufferListSize].blockLen = curBlockLen;
	
			bufferListSize++;
			if(bufferListSize > bufferListCap) //resize bufferList
			{
				emHomItem* oldBufferList = bufferList;
				bufferListCap *= 2;
				bufferList = (emHomItem *) malloc ((bufferListCap+2)*sizeof(emHomItem));
				for(int k=0; k<bufferListSize; k++)
				{
					bufferList[k] = oldBufferList[k];	
				}	
				free(oldBufferList);
			}
		}
	}
	free(bufferList);

	//Assigning reverse lookups
	//TODO(denizy)  Move these look up values to pre-processing
	unsigned int lookInd = 2;
	for(unsigned int i=1; i <= emHomClassListSize && lookInd <= emHomClassListSize; i++)
	{
		unsigned char curChrCode = emHomClassList[i].chrCode;
		unsigned int curChrPos = emHomClassList[i].chrPos;
		unsigned int lookChrPosLimit = curChrPos + emHomClassList[i].list[0].blockLen - 1; //first item in the list always the longest
		
		while(emHomClassList[lookInd].chrCode == curChrCode && emHomClassList[lookInd].chrPos <= lookChrPosLimit)
		{
			emHomClassList[lookInd].maxRevLookup = lookInd - i;
			lookInd++;
		}
	}
}

//Load the compact exact homology table to memory (format is a list of equivalence class representative linked to other items in the class
void LoadEquivClasses(char* perfTableFileName)
{
	//open itemCounts file that contains the length of the files.
	ifstream finItemCounts((string(perfTableFileName) + ".itemCount").c_str());

	if(!finItemCounts.is_open())
	{
		cout << "Homology Table item count file doesn't exist: " << string(perfTableFileName) + ".itemCount" << endl;
		exit(87);
	}

	finItemCounts >> equivClassListSize;
	equivClassList = (equivClassNode *) malloc ((equivClassListSize+2)*sizeof(equivClassNode));

	int bufferListCap = 10000;
	int bufferListSize = 0;
	eqItem* bufferList = (eqItem *) malloc ((bufferListCap+2)*sizeof(eqItem));

	FILE* finBin = fopen(perfTableFileName, "rb");

	if(!finBin)
	{
		cout << "Homology Table file doesn't exist: " << perfTableFileName << endl;
		exit(88);
	}

	for(unsigned int i=1; i<=equivClassListSize; i++)
	{
		size_t tempSize = fread(&(equivClassList[i].chrCode), sizeof(unsigned char), 1, finBin);
		tempSize &= fread(&(equivClassList[i].chrPos), sizeof(unsigned int), 1, finBin);

		bufferListSize = 0;
		while(true)
		{
			tempSize &= fread(&(bufferList[bufferListSize].dir), sizeof(unsigned char), 1, finBin);
			if(bufferList[bufferListSize].dir == -1)
			{
				//End of list -- read the block length, copy buffer to equivClassList and quit
				tempSize &= fread(&(equivClassList[i].blockLen), sizeof(unsigned int), 1, finBin);
				equivClassList[i].list = (eqItem *) malloc ((bufferListSize) * sizeof(eqItem));
				for(int k=0; k<bufferListSize; k++)
				{
					equivClassList[i].list[k] = bufferList[k];
				}
		
				equivClassList[i].listSize = bufferListSize;

				bufferListSize = 0;
				break;
			}

			tempSize &= fread(&(bufferList[bufferListSize].chrCode), sizeof(unsigned char), 1, finBin);
			tempSize &= fread(&(bufferList[bufferListSize].chrPos), sizeof(unsigned int), 1, finBin);

			bufferListSize++;
			if(bufferListSize > bufferListCap) //resize bufferList
			{
				eqItem* oldBufferList = bufferList;
				bufferListCap *= 2;
				bufferList = (eqItem *) malloc ((bufferListCap+2)*sizeof(eqItem));
				for(int k=0; k<bufferListSize; k++)
				{
					bufferList[k] = oldBufferList[k];	
				}	
				free(oldBufferList);
			}
		}
	}
	free(bufferList);
}

//Convert edit offsets and characters to their reverse complements
void ConvertToRevCompMedit(medit& curMedit, int readLen)
{
	unsigned char midPoint = (curMedit.numEdits + 1) / 2; //+1 since we want mid item to be reverse complemented too
	for(unsigned char k=0; k<midPoint; k++)
	{
		unsigned char secondPos = curMedit.numEdits-k-1;
		unsigned char temppos1 = readLen - curMedit.pos[secondPos] + 1;
		unsigned char tempch1 = RevCompChar[curMedit.edits[secondPos]];	
		
		if(tempch1 >= 'a' && tempch1 <= 't') //Indels should be incremented extra
			temppos1++;
		
		curMedit.edits[secondPos] = RevCompChar[curMedit.edits[k]];	
		curMedit.pos[secondPos] = readLen - curMedit.pos[k] + 1;
		if(curMedit.edits[k] >= 'a' && curMedit.edits[k] <= 't')
		{
			curMedit.pos[secondPos]++;
		}	

		curMedit.pos[k] = temppos1;
		curMedit.edits[k] = tempch1;		
	}
}

void ConvertToRevCompMedit(doubleMedit& curMedit, int readLen)
{
	unsigned char midPoint = (curMedit.numEdits + 1) / 2; //+1 since we want mid item to be reverse complemented too
	for(unsigned char k=0; k<midPoint; k++)
	{
		unsigned char secondPos = curMedit.numEdits-k-1;
		unsigned char temppos1 = readLen - curMedit.pos[secondPos] + 1;
		unsigned char tempch1 = RevCompChar[curMedit.edits[secondPos]];	
		
		if(tempch1 >= 'a' && tempch1 <= 't') //Indels should be incremented extra
			temppos1++;
		
		curMedit.edits[secondPos] = RevCompChar[curMedit.edits[k]];	
		curMedit.pos[secondPos] = readLen - curMedit.pos[k] + 1;
		if(curMedit.edits[k] >= 'a' && curMedit.edits[k] <= 't')
		{
			curMedit.pos[secondPos]++;
		}	

		curMedit.pos[k] = temppos1;
		curMedit.edits[k] = tempch1;		
	}
}

string DebugPrintMedit(const medit& d, double Code)
{
	bool problemMeditFLAG = 0;
	stringstream sout;

	unsigned char numEdits = d.numEdits % MISSING_SPLIT_OFFSET;
	sout << "numEdits: " << (int) d.numEdits << " : ";

	if(numEdits > MAX_EDITS)
	{
		cout << "!!!WARNING: NumEdits higher than allowed" << endl; 
		cout << "original num edits: " << (int) d.numEdits << endl;
		cout << "mod num Edits: " << (int) numEdits << endl;
	}

	for(int i=0; i<numEdits; i++)
	{
		if((int) d.pos[i] == 0)
		{
			cout << "!!!!PROBLEM medit!!!!!" << endl;
			problemMeditFLAG = 1;
		}
		sout << "(" << (int) d.pos[i] << "," << d.edits[i] << ") ";
	}

	sout << " Key: " << Code;
	
	if(problemMeditFLAG == 1)
	{

		cout << "here is problem medit: '" << sout.str() << "'" << endl;
		for(int i=0; i<numEdits; i++)
		{
			cout << "Edit " << i << ": " << (int) d.pos[i] << "\t'" << d.edits[i] << "'   int: " << (int) d.edits[i] << endl;  
		}
		problemMeditFLAG = 0;
		assert(0);
	}

	for(int i=1; i<numEdits; i++)
	{
		if(d.pos[i] < d.pos[i-1])
		{
			cout << "CORRUPT_ORDER: Here is corrupted medit for position order: '" << sout.str() <<  "'" << endl;
		}
		if(d.pos[i] == d.pos[i-1] && d.edits[i-1] < 91)
		{
			cout << "CORRUPT_EQUALITY: Here is the corrupted medit: '" << sout.str() << "'" << endl;
		}
	}

	return sout.str();
}

string DebugPrintMedit(const doubleMedit& d, double Code)
{
	bool problemMeditFLAG = 0;
	stringstream sout;

	unsigned char numEdits = d.numEdits % MISSING_SPLIT_OFFSET;
	sout << "numEdits: " << (int) d.numEdits << " : ";

	if(numEdits > MAX_DOUBLE_EDITS)
	{
		cout << "!!!WARNING: NumEdits higher than allowed" << endl; 
		cout << "original num edits: " << (int) d.numEdits << endl;
		cout << "mod num Edits: " << (int) numEdits << endl;
	}
	
	for(int i=0; i<numEdits; i++)
	{
		if((int) d.pos[i] == 0)
		{
			cout << "!!!!PROBLEM medit!!!!!" << endl;
			problemMeditFLAG = 1;
		}
		sout << "(" << (int) d.pos[i] << "," << d.edits[i] << ") ";
	}

	sout << " Key: " << Code;
	
	if(problemMeditFLAG == 1)
	{

		cout << "here is problem medit: '" << sout.str() << "'" << endl;
		for(int i=0; i<numEdits; i++)
		{
			cout << "Edit " << i << ": " << (int) d.pos[i] << "\t'" << d.edits[i] << "'   int: " << (int) d.edits[i] << endl;  
		}
		problemMeditFLAG = 0;
		assert(0);
	}

	for(int i=1; i<numEdits; i++)
	{
		if(d.pos[i] < d.pos[i-1])
		{
			cout << "CORRUPT_ORDER: Here is corrupted medit for position order: '" << sout.str() <<  "'" << endl;
		}
		if(d.pos[i] == d.pos[i-1] && d.edits[i-1] < 91)
		{
			cout << "CORRUPT_EQUALITY: Here is the corrupted medit: '" << sout.str() << "'" << endl;
		}
	}

	return sout.str();
}

//Create MD edit script from medit while replacing the character for the reference (since medit is the edits in the read and the MD string is the edits in the reference)
void GetMDStrFromMeditWithCharReplacement(const doubleMedit& curMedit, char replacedChars[], int finalReadLen, char outarr[])
{
	//Alignment dir is assumed to be forward
	char* md = outarr;	
	unsigned char numEdits = curMedit.numEdits;
	if(numEdits==0)
	{
		writeStr(md, valStr[finalReadLen], '\0');
	}
	else
	{
		int numInsertionsWithinBlock = 0; //This resets whenever a non-insertion variation is seen 
		int distanceFromVal = 1;

		for(int k=0; k<=numEdits; k++)
		{
			unsigned char curPosVal = finalReadLen + 1;
			if(k < numEdits)
			{
				if(curMedit.edits[k] >= 'a' && curMedit.edits[k] <= 't') //insertions are not reported in MD str
				{
					numInsertionsWithinBlock++;
					continue;
				}
				curPosVal = curMedit.pos[k];
			}

			unsigned char diffMD = curPosVal - distanceFromVal;// + numDeletionsWithinBlock;

			if(k == numEdits)
			{
				diffMD -= numInsertionsWithinBlock;
			}

			distanceFromVal += diffMD + 1; //this will be the position that we start for the next character

			if(diffMD > 0)
			{
				if((md-1)[0] == '0')
				{
					md--;
				}
	
				writeStr(md, valStr[diffMD]);
			}
			
			if(k < numEdits)
			{
				if(curMedit.edits[k] == 'D')
				{
					md[0] = '^'; md++;
					md[0] = replacedChars[k]; md++;
					md[0] = '0'; md++; //to make sure that the next character will not be confused with a mismatch
					numInsertionsWithinBlock--;
				}
				else
				{
					md[0] = replacedChars[k]; md++;
				}
			}
		}			
	}
	md[0] = '\0';
}

//Construct medit from a given editstring (the way it's represented in the links)
medit GetMedit(const string& editStr, bool& successful)
{
	medit curMedit = dummyMedit;

	int fullOffset = 0;
	int editLen = editStr.length();
	for(int i=0; i<editLen; i++)
	{
		int curOffset = 0;
		while(i<editLen && editStr[i]<= '9')
		{
			curOffset *= 10;
			curOffset += editStr[i] - '0';
			i++;
		}		 
	
		if(i>=editLen)
		{
			break;
		}

		fullOffset += curOffset;

		if(curMedit.numEdits >= MAX_EDITS)
		{
			successful = 0;
			break;
		}

		curMedit.pos[curMedit.numEdits] = fullOffset + 1;
		curMedit.edits[curMedit.numEdits] = editStr[i];
		curMedit.numEdits++;

		if(editStr[i] >= 'A' && editStr[i] <= 'T') //Includes Deletions 'D' and all mismatches 'ACGTN'
			fullOffset++;
	}

	return curMedit;
}

//Go through fai file and load all chromosomes to memory (1 char per nucleotide) - 1 based chromosome ids, 1 based position indices
void LoadMultiChrReference(char* refFile, int refLineLen)
{
	string refFileName = refFile;
	string faiFileName = refFileName + ".fai";
	
	cout << "Reading fai file and allocating reference" << endl;

	//Verifying refFileFai & reference size
        ifstream finFai(faiFileName.c_str());
        if(!finFai.is_open())
        {
                cout << "Could not find fai file for reference file: " << refFile << endl;
                exit(16);
        }

        //In MultiChr Version only MaxChrSize is needed
        int chrCode = 0;
        string faiLine;
        while(getline(finFai, faiLine))
        {
                stringstream faiLineSS(faiLine);
                string ctgName;
                int ctgSize;
                faiLineSS >> ctgName >> ctgSize;

                chrCode++;
                fullRef[chrCode] = (char *) calloc (ctgSize+MAX_LINE_LEN, sizeof(char));
		chrLens[chrCode] = ctgSize;
        	strcpy(chrNames[chrCode], ctgName.c_str());
	}
        numChrs = chrCode;

        //Dummy space allocation for the 0th chromosome (will not be used)
        fullRef[0] = (char *) malloc (500);

        finFai.close();

        FILE* finRef = fopen(refFileName.c_str(),"r");
	if(!finRef)
	{
		cout << "Could not find reference file: " << refFile << endl;
		exit(17);
	}

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
	fclose(finRef);
}

//Handles loading both inexact and exact homology tables, as well as setting up the posHub links
void LoadHomologyTables(char* perfTableFileName, char* inexactTableFileName, int readLen, int numHomTableMismatches, int mapMode)
{
	//Homology table is in two pieces, perfect equivalence table (block-blocks) and inexact homology table (node-blocks)
	cout << "Loading perfect homology table.." << endl;
	LoadEquivClasses(perfTableFileName);

	if(mapMode != BEST_MAPPING_MODE && mapMode != BEST_FAST_MAPPING_MODE) //this is the faster best mapping mode
		LoadInexactHomologies(inexactTableFileName, readLen, numHomTableMismatches);

	if(mapMode != BEST_MAPPING_MODE && mapMode != BEST_FAST_MAPPING_MODE)
	{
		//First set up emHoms for all representatives (non-representatives will not be here)	
		for(unsigned int i=1; i<= emHomClassListSize; i++)
		{
			//Assign pointers of posNodes
			homHub[(int) emHomClassList[i].chrCode][emHomClassList[i].chrPos] = i;
		}

		//Second set-up all equivClass posHubs (while carrying over emHomologies)
		for(unsigned int i=1; i<= equivClassListSize; i++)
		{
			//Assign pointers of posNodes
			equivClassNode& curEqNode = equivClassList[i];
			unsigned int curEqBlockLen = curEqNode.blockLen;
			for(unsigned int k=0; k <= curEqBlockLen-1; k++) //dir of the representative is always forward
			{
				equivHub[(short) curEqNode.chrCode][curEqNode.chrPos + k] = i;
			} 

			for(unsigned int k=0; k<curEqNode.listSize; k++)
			{
				eqItem& curHomItem = curEqNode.list[k];
				if(curHomItem.dir == 0)
				{
					for(unsigned int j=0; j<= curEqBlockLen-1; j++)
					{
						equivHub[(short) curHomItem.chrCode][curHomItem.chrPos + j] = i;
					
						//Copy from representative
						homHub[(short) curHomItem.chrCode][curHomItem.chrPos + j] = homHub[(short) curEqNode.chrCode][curEqNode.chrPos + j];
					}
				}
				else
				{
					for(unsigned int j=0; j<= curEqBlockLen-1; j++)
					{
						equivHub[(short) curHomItem.chrCode][curHomItem.chrPos - j] = i;					
	
						//Copy from representative
						homHub[(short) curHomItem.chrCode][curHomItem.chrPos - j] = homHub[(short) curEqNode.chrCode][curEqNode.chrPos + j];
					}
				}
			}
		}
	}
	else
	{
		//Second set-up all equivClass posHubs (without emHomologies)
		for(unsigned int i=1; i<= equivClassListSize; i++)
		{
			//Assign pointers of posNodes
			equivClassNode& curEqNode = equivClassList[i];
			unsigned int curEqBlockLen = curEqNode.blockLen;
			for(unsigned int k=0; k <= curEqBlockLen-1; k++) //dir of the representative is always forward
			{
				equivHub[(short) curEqNode.chrCode][curEqNode.chrPos + k] = i;
			} 
		}
	}
}

////////////////////////////////
/// Variables and functions relevant to LRNs ( long read names that would crash the off-the-shelf mapper so these are handles independently)
string** longList;
unsigned int longListSize; //total number of items in the current LRN list
unsigned int longListCap; //current capacity of the LRN list
#define LONG_LIST_START_SIZE 100000

//Setup the list of LRNs and do some reallocations as needed
void SetupLongNameList(ifstream& finLong)
{
	cout << "Setting up long name list.." << endl;

	//Setup initial array
	longList = (string **) malloc (sizeof(string*) * (LONG_LIST_START_SIZE +2));
	longListCap = LONG_LIST_START_SIZE;

	string num, name;
	while( finLong >> num >> name)
	{
		longListSize++;
		longList[longListSize] = new string(name);
		if(longListSize >= longListCap)
		{
		//TODO(denizy) Look into this function for potential mem leak (valgrind warning)
			string** tempLongList = (string **) malloc (sizeof(string*) * (2 * longListCap+ 2));	
			longListCap *= 2;
			for(unsigned int k=1; k<=longListSize; k++)
			{
				tempLongList[k] = longList[k];
			}		
			delete longList;
			longList = tempLongList;
		}	
	}
	cout << "Number of long names intercepted: " << longListSize << endl;
}

//For a given pseudo-readName (e.g #34###), return the LRN corresponding to the integer inside it
string GetLongReadName(string code)
{
	int secondHash = code.find("#", 1);
	int longCode = atoi(code.substr(1, secondHash-1).c_str());
	return *(longList[longCode]);
}

char readLenStr[4]; //character array for making integer printing to the buffer faster
char replacedQueryCharList[MAX_DOUBLE_EDITS+2]; //list of replacement characters when merging to medits
unsigned long long GetReadIdFromString(const string& str);

unsigned long long GetReadIdFromCharArray(char arr[])
{
	unsigned long long val = 0;
	for(int i=0; i<idDigitLen; i++)
	{
		val *= IDENTITY_ALPHABET_SIZE; 
		val += (arr[i] - IDENTITY_ALPHABET_START);
	}
	return val;
}

//zero and max are the begining and the ending of the current loaded links table (since if links table is very large it is loaded in multiple batches)
unsigned long long curZero;
unsigned long long curMax;
unsigned long long maxRepositorySize = UINT_MAX / 2;

struct memoInd
{
	int begInd; //inclusive beginning location of a memoization block
	int endInd; //inclusive ending location 
	unsigned int chrPos; //The same as chrCode, but the position in the chromosome
	medit edit;
	unsigned char chrCode; //This is for memoization marked links to get their chrCode from (since we don't know which readID will be processed first)
	unsigned char properReadLength; //length of the genome alignment after indels
};

int curTotalMemoSize; 
vector<memoInd> memoInds; //This keeps track of which positions in the memoArr belong to the current memoization block
vector<link> memoArr; //This is the array that stores all the memoized info
int MEMOIZATION_THRESHOLD;

#define MEMO_CHR_CODE 250	 //This the ChrCode flag that tells that the corresponding link is not onto the genome but onto a memoized block

//Loads links table form the unsorted input file (only load the needed batch)
bool LoadLinksRepositoryFromUnsorted(char* linksFileName, unsigned long long numReads, short numTraversal)
{
	curZero = maxRepositorySize * numTraversal;
	curMax = curZero + maxRepositorySize;

	if(numReads > curMax)
	{
		cout << "ERROR: Current memoization version doesn't support large link tables -- update to enable it" << endl;
		exit(0);
	}	

	unsigned long long totalSizeToAllocate = numReads; //later on think about multipass solutions

	if(splitMode == HALF_SPLIT)
		totalSizeToAllocate *= 2;
	if(inputMode == PAIRED_MODE)
		totalSizeToAllocate *= 2;
	
	if(numTraversal == 0)
	{
		if(totalSizeToAllocate > maxRepositorySize)
			totalSizeToAllocate = maxRepositorySize; 

		links = (link *) calloc (totalSizeToAllocate, sizeof(link));
		linksSize = totalSizeToAllocate;
	}
	else
	{
		//Clean Repository and assign new size

		linksSize = min(maxRepositorySize, totalSizeToAllocate - curZero);

		for(unsigned int i=0; i<linksSize; i++)
		{
			//it's enough to just reset chrCode -- it automatically becomes unusable		
			links[i].chrCode = 0;
		}
	}
	
	//TODO(denizy) improve IO speed with faster freads and buffering here

	ifstream finLinks(linksFileName);

	string linkNameCode, linkGenEdit; //links file contains the encoded read names of the link as well as the MD edit scripts
	unsigned int linkChrPos, linkChrCode; //and chromosome and position that the link points towards in the refernece genome

	int nameMemoThresh = MEMOIZATION_THRESHOLD * idDigitLen; //Minimum length for name that will have results memoized 

	//Initialize memoization (here insert a dummy item so that the first index becomes 1 -- othwerwise chrPos = 0 can be confused with other flags
	memoInd dummyMemoInd;
	dummyMemoInd.begInd = -1;
	dummyMemoInd.endInd = -1;
	dummyMemoInd.chrCode = 0;
	dummyMemoInd.chrPos = 0;
	
	memoInds.push_back(dummyMemoInd);
	int numMemoed = 1; //Total number of memoized readMers (each has a construct Maps block)
	
	string linkLine;
	while(getline(finLinks, linkLine))
	{
		stringstream linkLineSS(linkLine);
		
		linkLineSS >> linkNameCode >> linkChrCode >> linkChrPos >> linkGenEdit;

		int nameLen = linkNameCode.length();
		
		if(nameLen % idDigitLen != 0)
		{
			//Assing proper long name 
			linkNameCode = GetLongReadName(linkNameCode);
			nameLen = linkNameCode.length();
		}

		//Make sure that the link provided by coarse mapping is within the errorPerKmer limits, otherwise skip
		bool successfulMeditReturn = 1;
		medit curLinkMedit = GetMedit(linkGenEdit, successfulMeditReturn);
		if(successfulMeditReturn == 0)
		{
			continue;
		}
	
		unsigned long long readId;

		if(nameLen < nameMemoThresh)
		{
			for(int i=0; i<nameLen; i += idDigitLen)
			{
				readId = GetReadIdFromString(linkNameCode.substr(i, idDigitLen));

				bool dirFlag = readId % 2;
				readId /= 2;
				
				if(readId >= curZero && readId < curMax)
				{
					readId -= curZero;

					links[readId].chrCode = linkChrCode;
					links[readId].chrPos = linkChrPos;
					links[readId].edit = curLinkMedit;
					links[readId].dir = dirFlag;
				}
			}
		}
		else //setup memoization
		{
			memoInd newMemoInd;
			newMemoInd.begInd = -1;
			newMemoInd.endInd = -1;
			newMemoInd.chrCode = linkChrCode;
			newMemoInd.chrPos = linkChrPos;
			newMemoInd.edit = curLinkMedit;

			memoInds.push_back(newMemoInd);

			for(int i=0; i<nameLen; i += idDigitLen)
			{
				readId = GetReadIdFromString(linkNameCode.substr(i, idDigitLen));

				bool dirFlag = readId % 2;
				readId /= 2;
				
				if(readId >= curZero && readId < curMax)
				{
					readId -= curZero;
			
					links[readId].chrCode = MEMO_CHR_CODE;
					links[readId].chrPos = numMemoed;
					links[readId].dir = dirFlag;
				}
			}


			numMemoed++;
		}
	}
	
	finLinks.clear();
	finLinks.close();

	if(numReads <= curMax) //All read links traversed in the input, no more rescans needed
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

vector<mergedLink> travArr[6]; //This is the list of links for each type of readend, merged end, etc.
unsigned char travArrReadmerLength[6]; //This the total length of splits so that length modifiying indels can be accomodated;

int SplaySplit(short chrNo, unsigned int chrPos, const doubleMedit& MD, bool dir, char splay[], int properReadLen);

// Debug function for verifying internal state of the medits and check for inconsistencies
void CheckForMeditInconsistencies(unsigned char pushChr,  unsigned int pushPos, medit pushMedit, bool pushDir, unsigned int readId, double code, int properReadLen)
{
	unsigned char linkChr = links[readId].chrCode;
	unsigned int linkPos = links[readId].chrPos;
	doubleMedit linkMedit = links[readId].edit;
	bool linkDir = links[readId].dir;

	if(linkChr == 0)
	{
		//might extract it from nohit array but not really necessary for completeness of debugging
		return;
	}

	//Restore sequnce from linkDir
	char splayedRead[MAX_READ_LEN + 10];
	
	if(SplaySplit(linkChr, linkPos, linkMedit, linkDir, splayedRead, properReadLen))
	{
		char splayedPush[MAX_READ_LEN + 10];
	
		if(!SplaySplit(pushChr, pushPos, pushMedit, pushDir, splayedPush, properReadLen))
		{
			cout << "CORRUPT_SPLAY_ERROR" << endl;

			cout << "Code: " << code << endl;
			cout << "pushChr: " << (int) pushChr << "\tpushPos:" << pushPos << "\tpushMedit: " << DebugPrintMedit(pushMedit, 2) << "\tpushDir: " << pushDir << endl;
			cout << "linkChr: " << (int) linkChr << "\tlinkPos:" << linkPos << "\tlinkMedit: " << DebugPrintMedit(linkMedit, 3) << "\tlinkDir: " << linkDir << endl;
			cout << "Code: " << code << endl;
			cout << "ReadId: " << readId << " readNo: " << readId/4 + 1 << endl;
			cout << "PUSHref: " << string(fullRef[pushChr] + pushPos, properReadLen) << "\t" << endl;
			cout << "LINKref: " << string(fullRef[linkChr] + linkPos, properReadLen) << "\t" << endl;
			cout << "pushPos: " << pushPos << "\tdir: " << pushDir << "\t" << DebugPrintMedit(pushMedit, 4) << endl;
			cout << "linkPos: " << linkPos << "\tdir: " << linkDir << "\t" << DebugPrintMedit(linkMedit, 5) << endl;
			cout << "SLINKrd: " << splayedRead << endl;
		}		
	}	
}

doubleMedit savedEquivRepEdit;

//Get all relevant exact homologies
void TraverseEquivClass(unsigned int readId, unsigned char repChrCode, unsigned int repChrPos, unsigned int curEqClassIndex, int equivBlockOffset, bool dirFromHitToRep, unsigned char arrNo)
{
	int properReadLen = travArrReadmerLength[arrNo]; //considering the indel liength modifications to the read

	//No need to do medit merges, just print the actual medit with depending on the direction (since they are all identical)
	equivClassNode& curEqClass = equivClassList[curEqClassIndex];
		
	medit linkMedit =  links[readId].edit;
	medit revCompLinkMedit = linkMedit;

	bool repIndelOffsetFlag = 0;

	if(dirFromHitToRep == 0)
	{
		ConvertToRevCompMedit(revCompLinkMedit, properReadLen);
	}	
	else
	{
		ConvertToRevCompMedit(linkMedit, properReadLen);
		repChrPos += readLen - properReadLen;
		repIndelOffsetFlag = 1;
	}

	if(links[readId].dir == 1)
	{
		dirFromHitToRep = !dirFromHitToRep;
	}

	travArr[arrNo].push_back(mergedLink(repChrPos + equivBlockOffset, linkMedit, repChrCode, dirFromHitToRep)); 

	if(properReadLen > readLen) //Deletion case --> might need to add more mismatches to the rep
	{
		mergedLink &lastLink = travArr[arrNo][travArr[arrNo].size()-1];	
		
		if(links[readId].dir == dirFromHitToRep)
		{
			//Addition to the end
			for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
			{
				if( fullRef[repChrCode][lastLink.chrPos + nucInd] != fullRef[links[readId].chrCode][links[readId].chrPos + nucInd] )
				{
					if(lastLink.edit.pos[lastLink.edit.numEdits-1] <= nucInd) //NO conflicts when merging so merge directly
					{
						lastLink.edit.pos[lastLink.edit.numEdits] = nucInd+1;
						lastLink.edit.edits[lastLink.edit.numEdits] = fullRef[links[readId].chrCode][links[readId].chrPos + nucInd];
						lastLink.edit.numEdits++;
					}
					else
					{
						int insertionPos = lastLink.edit.numEdits-1;
						while(lastLink.edit.pos[insertionPos] > nucInd + 1)
						{
							lastLink.edit.pos[insertionPos+1] = lastLink.edit.pos[insertionPos];
							lastLink.edit.edits[insertionPos+1] = lastLink.edit.edits[insertionPos];
							insertionPos--;
						}
						
						if(lastLink.edit.pos[insertionPos] == nucInd + 1)
						{
							if(fullRef[lastLink.chrCode][lastLink.chrPos + nucInd] == lastLink.edit.edits[insertionPos])
							{
								//Delete the edit since it negates it -- then shift the remaining edits back
								lastLink.edit.numEdits--;
								for(int i=insertionPos; i<lastLink.edit.numEdits; i++)
								{
									lastLink.edit.pos[insertionPos] = lastLink.edit.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
									lastLink.edit.edits[insertionPos] = lastLink.edit.edits[insertionPos+2];
								}
							}							
							//else do nothing since the same edit would be valid
						}
						else
						{
							lastLink.edit.pos[insertionPos] = nucInd + 1;
							lastLink.edit.edits[insertionPos] = fullRef[links[readId].chrCode][links[readId].chrPos + nucInd];
							lastLink.edit.numEdits++;
						}
					}
				}
			}
		}
		else
		{
			//Addition to the beginning
			for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
			{
				if( fullRef[repChrCode][lastLink.chrPos + properReadLen - nucInd -1] != RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]] )
				{
					if(lastLink.edit.pos[0] > properReadLen - nucInd)
					{
						for(int i=lastLink.edit.numEdits; i>0; i--)
						{
							lastLink.edit.pos[i] = lastLink.edit.pos[i-1];
							lastLink.edit.edits[i] = lastLink.edit.edits[i-1];
						}
						lastLink.edit.pos[0] = properReadLen - nucInd;
						lastLink.edit.edits[0] = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]];
						lastLink.edit.numEdits++;
					}
					else
					{
						int insertionPos = 0;
						while(lastLink.edit.pos[insertionPos] < properReadLen - nucInd)
						{
							insertionPos++;
						}
						
						if(lastLink.edit.pos[insertionPos] == properReadLen - nucInd)
						{
							if(fullRef[lastLink.chrCode][lastLink.chrPos + properReadLen - nucInd - 1] == RevCompChar[lastLink.edit.edits[insertionPos]])
							{
								//Delete the edit since it negates it -- then shift the remaining edits back
								lastLink.edit.numEdits--;
								for(int i=insertionPos; i<lastLink.edit.numEdits; i++)
								{
									lastLink.edit.pos[insertionPos] = lastLink.edit.pos[insertionPos+1];//Since we were pre-shifting, we now have to revert two steps back
									lastLink.edit.edits[insertionPos] = lastLink.edit.edits[insertionPos+1];
								}
							}
							//else do nothing since the same edit would be valid
						}
						else
						{
							for(int i=lastLink.edit.numEdits; i>insertionPos; i--)
							{
								lastLink.edit.pos[i] = lastLink.edit.pos[i-1];
								lastLink.edit.edits[i] = lastLink.edit.edits[i-1];
							}
							lastLink.edit.pos[insertionPos] = properReadLen - nucInd;
							lastLink.edit.edits[insertionPos] = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]];
							lastLink.edit.numEdits++;
						}

					}
				}
			}
		}

		savedEquivRepEdit = travArr[arrNo][travArr[arrNo].size()-1].edit;

		if(travArr[arrNo][travArr[arrNo].size()-1].edit.numEdits > numMismatchesPerReadMer)
		{
			travArr[arrNo].pop_back();
		}
	}

	if(curEqClassIndex != 0)
	{
		for(unsigned int k=0; k<curEqClass.listSize; k++)
		{
			eqItem& curEqItem = curEqClass.list[k];
			if(repIndelOffsetFlag == 0)
			{
				if(curEqItem.dir == 0)
				{
					travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset, linkMedit, curEqItem.chrCode, dirFromHitToRep));
				}
				else
				{
					travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset + (readLen - properReadLen), revCompLinkMedit, curEqItem.chrCode, !dirFromHitToRep));
				}
			}
			else
			{
				if(curEqItem.dir == 0)
				{
					travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset + (readLen - properReadLen), linkMedit, curEqItem.chrCode, dirFromHitToRep));
				}
				else
				{
					travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset, revCompLinkMedit, curEqItem.chrCode, !dirFromHitToRep));
				}
				
			}
			
			//Add additional potential mismatches due to deeltions (note that exact matches may no longer be exact matches in the case of deletions		
			if(properReadLen > readLen) //Deletion case
			{
				if(curEqItem.dir == !(links[readId].dir == dirFromHitToRep)) //Look at additional bases in the ]readLen, properReadLen] range and add any mismatches from the target.
				{
					mergedLink &lastLink = travArr[arrNo][travArr[arrNo].size()-1];

					for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
					{
						if( fullRef[lastLink.chrCode][lastLink.chrPos + nucInd] != fullRef[links[readId].chrCode][links[readId].chrPos + nucInd] )
						{
							if(lastLink.edit.pos[lastLink.edit.numEdits-1] <= nucInd) //NO conflicts when merging so merge directly
							{
								lastLink.edit.pos[lastLink.edit.numEdits] = nucInd+1;
								lastLink.edit.edits[lastLink.edit.numEdits] = fullRef[links[readId].chrCode][links[readId].chrPos + nucInd];
								lastLink.edit.numEdits++;
							}
							else //Need to handle some tricky merging
							{
								int insertionPos = lastLink.edit.numEdits-1;
								while(lastLink.edit.pos[insertionPos] > nucInd + 1)
								{
									lastLink.edit.pos[insertionPos+1] = lastLink.edit.pos[insertionPos];
									lastLink.edit.edits[insertionPos+1] = lastLink.edit.edits[insertionPos];
									insertionPos--;
								}
								
								if(lastLink.edit.pos[insertionPos] == nucInd + 1)
								{
									if(fullRef[lastLink.chrCode][lastLink.chrPos + nucInd] == lastLink.edit.edits[insertionPos])
									{
										//Delete the edit since it negates it -- then shift the remaining edits back
										lastLink.edit.numEdits--;
										for(int i=insertionPos; i<lastLink.edit.numEdits; i++)
										{
											lastLink.edit.pos[insertionPos] = lastLink.edit.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
											lastLink.edit.edits[insertionPos] = lastLink.edit.edits[insertionPos+2];
										}
									}
									//else do nothing since the same edit would be valid							
								}
								else
								{
									lastLink.edit.pos[insertionPos] = nucInd + 1;
									lastLink.edit.edits[insertionPos] = fullRef[links[readId].chrCode][links[readId].chrPos + nucInd];
									lastLink.edit.numEdits++;
								}
							}
						}
					}
				}
				else //Look at additional bases in the ]readLen, properReadLen] range from the source but compare to bases before the beginning position of the target and with the complement characters.
				{
					mergedLink &lastLink = travArr[arrNo][travArr[arrNo].size()-1];	
					for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
					{
						if( fullRef[lastLink.chrCode][lastLink.chrPos + properReadLen - nucInd -1] != RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]])
						{
							if(lastLink.edit.pos[0] > properReadLen - nucInd)
							{
								for(int i=lastLink.edit.numEdits; i>0; i--)
								{
									lastLink.edit.pos[i] = lastLink.edit.pos[i-1];
									lastLink.edit.edits[i] = lastLink.edit.edits[i-1];
								}
								lastLink.edit.pos[0] = properReadLen - nucInd;
								lastLink.edit.edits[0] = RevCompChar[ (unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]];
								lastLink.edit.numEdits++;
							}
							else
							{
								int insertionPos = 0;
								while(lastLink.edit.pos[insertionPos] < properReadLen - nucInd)
								{
									insertionPos++;
								}

								if(lastLink.edit.pos[insertionPos] == properReadLen - nucInd)
								{
									if(fullRef[lastLink.chrCode][lastLink.chrPos + properReadLen - nucInd - 1] == RevCompChar[(unsigned char) lastLink.edit.edits[insertionPos]])
									{
										//Delete the edit since it negates it -- then shift the remaining edits back
										lastLink.edit.numEdits--;
										for(int i=insertionPos; i<lastLink.edit.numEdits; i++)
										{
											lastLink.edit.pos[insertionPos] = lastLink.edit.pos[insertionPos+1];//Since we were pre-shifting, we now have to revert two steps back
											lastLink.edit.edits[insertionPos] = lastLink.edit.edits[insertionPos+1];
										}
									}
									//else do nothing since the same edit would be valid
								}
								else
								{
									for(int i=lastLink.edit.numEdits; i>insertionPos; i--)
									{
										lastLink.edit.pos[i] = lastLink.edit.pos[i-1];
										lastLink.edit.edits[i] = lastLink.edit.edits[i-1];
									}
									lastLink.edit.pos[insertionPos] = properReadLen - nucInd;
									lastLink.edit.edits[insertionPos] = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]];
									lastLink.edit.numEdits++;
								}
							}
						}
					}
				}

				if(travArr[arrNo][travArr[arrNo].size()-1].edit.numEdits > numMismatchesPerReadMer)
				{
					travArr[arrNo].pop_back();
				}
			}
		}
	}
}

unsigned int RetrieveOffsetAndDirInfo(unsigned char linkChrCode, unsigned int linkChrPos, short& equivBlockOffset, bool& dirFromHitToRep);

//Get all relevant inexact homologies
void TraverseInexactHomologyTable(unsigned int readId, unsigned char arrNo, unsigned int linkEquivClassIndex, unsigned short  linkEquivBlockOffset, bool directionFromLinkToEquivRep)
{
	int properReadLen = travArrReadmerLength[arrNo];

	unsigned char linkChrCode = links[readId].chrCode;
	unsigned int linkChrPos = links[readId].chrPos;
	doubleMedit linkMeditOrig(links[readId].edit); //regardless of direction of readLink, the medit should be in forward direction

	bool correctedLinkDir = (links[readId].dir != directionFromLinkToEquivRep); //this corrects for the direction change from link pointer to the actual class representative that we're working with
	
	if(directionFromLinkToEquivRep == 1)
	{
		ConvertToRevCompMedit(linkMeditOrig, properReadLen);
	}	

	if(properReadLen > readLen)
	{
		linkMeditOrig = savedEquivRepEdit; //This is in order the prevent the inconsistencies caused by missing edits from link-to-
	}

	unsigned char chrOfLinkEquivRep;
	unsigned int offsetPosOfTheLinkEquivRep;
	if(linkEquivClassIndex)
	{
		chrOfLinkEquivRep = equivClassList[linkEquivClassIndex].chrCode;
		offsetPosOfTheLinkEquivRep = equivClassList[linkEquivClassIndex].chrPos + linkEquivBlockOffset;
		if(directionFromLinkToEquivRep)
		{
			offsetPosOfTheLinkEquivRep += readLen - properReadLen;
		}
	}
	else
	{
		chrOfLinkEquivRep = linkChrCode;
		offsetPosOfTheLinkEquivRep = linkChrPos;
	}

	//TODO(denizy) Implement faster lookup pointer rather than searching backwards like this
	unsigned int curEmHomIndex = 0;
	for(unsigned int i=0; i<(unsigned short) readLen && i < offsetPosOfTheLinkEquivRep && curEmHomIndex == 0 && i<linkChrPos; i++)
        {
		std::tr1::unordered_map<unsigned int, unsigned int>::iterator ite = homHub[chrOfLinkEquivRep].find(offsetPosOfTheLinkEquivRep-i);
                if(ite != homHub[chrOfLinkEquivRep].end())
                {
                       curEmHomIndex = ite->second;
                }
	}

        unsigned char reverseBackupLim = emHomClassList[curEmHomIndex].maxRevLookup;
	for(unsigned int i=0; i<=reverseBackupLim; i++) //blockLength of inexactHom cannot be larger than readLen
        {
                emHomClassNode& curHomClass = emHomClassList[curEmHomIndex-i];
		unsigned int listSize = curHomClass.listSize;

                for(unsigned int k=0; k<listSize; k++)
                {
			unsigned char offsetFromCurrentBlock = offsetPosOfTheLinkEquivRep - curHomClass.chrPos;
	
			emHomItem& curEmHomItem = curHomClass.list[k];

                        if(curEmHomItem.blockLen > offsetFromCurrentBlock) //then we're interested in this link
                        {
				//This checks if the homology item is also within an equivalence class (in which case we'd need to properly scan this homology class to avoid double printing
				unsigned int curHomClass_offsetChrPos = curHomClass.chrPos + offsetFromCurrentBlock;
				unsigned int curEmHomItem_offsetChrPos = curEmHomItem.chrPos + offsetFromCurrentBlock;
				doubleMedit augmentedLinkMedit = linkMeditOrig;
				
				if(curEmHomItem.dir != 0)
				{
					curEmHomItem_offsetChrPos = curEmHomItem.chrPos - offsetFromCurrentBlock;
				}
				
				//We would like to construct the read link to the curHomItem with proper offsets
				for(unsigned char offInd = 0; offInd < numMismatchesPerReadMer && curEmHomItem.offsetList[offInd] != 0; offInd++) //offset positions are always in the forward direction of the 
				{
					if(curEmHomItem.offsetList[offInd] > offsetFromCurrentBlock)
					{
						unsigned char editOffsetToAdd = curEmHomItem.offsetList[offInd] - offsetFromCurrentBlock;

						char editCharToAddFromHomClass = fullRef[curHomClass.chrCode][curHomClass_offsetChrPos + editOffsetToAdd - 1];  //this is the character in emHomClass that is different from emHomItem (directionality not important (-1 is beacuse offsets are 1 based where are ref pointer directly points at 0
						char charWithinHomItem;
						if(curEmHomItem.dir == 0)
						{
							charWithinHomItem = fullRef[curEmHomItem.chrCode][curEmHomItem_offsetChrPos + editOffsetToAdd - 1]; 
						}
						else
						{
							charWithinHomItem = RevCompChar[(unsigned char) fullRef[curEmHomItem.chrCode][curEmHomItem_offsetChrPos + readLen - editOffsetToAdd]];
						}
						
						bool actionTakenFlag = 0;
						for(unsigned char lind = 0; lind < augmentedLinkMedit.numEdits; lind++)
						{
							if(augmentedLinkMedit.pos[lind] > editOffsetToAdd)
							{
								//Shift and insert
								for(unsigned char lend = augmentedLinkMedit.numEdits; lend > lind; lend--)
								{
									augmentedLinkMedit.pos[lend] = augmentedLinkMedit.pos[lend-1];
									augmentedLinkMedit.edits[lend] = augmentedLinkMedit.edits[lend-1];
								}
								augmentedLinkMedit.pos[lind] = editOffsetToAdd;
								//the character to add
								augmentedLinkMedit.edits[lind] = editCharToAddFromHomClass;
								augmentedLinkMedit.numEdits++;
								actionTakenFlag = 1;
								break;
							}
							else if(augmentedLinkMedit.pos[lind] == editOffsetToAdd)
							{
								//two scenarios: 
								//(1) readLink's editChar is same as the homItem's character at the proper position (in which case the edit will be deleted from medit.
		
								//Shift for all the insertions with the same position index (which could happen due to chaining)
								if(augmentedLinkMedit.edits[lind] >= 'a' && augmentedLinkMedit.edits[lind] <= 't')
								{
									lind++;
									while(lind < augmentedLinkMedit.numEdits && augmentedLinkMedit.pos[lind] == editOffsetToAdd && augmentedLinkMedit.edits[lind] >= 'a' && augmentedLinkMedit.edits[lind] <= 't')
									{
										lind++;
									}
								}
								
								//Now check again if there is going to be delete operation or an insert operation to the medits 
								if(augmentedLinkMedit.pos[lind] == editOffsetToAdd)
								{
									if(augmentedLinkMedit.edits[lind] == charWithinHomItem)
									{
										//deleted edit from medit -- by shifting the remaining edits left
										for(unsigned char ldel = lind + 1; ldel < augmentedLinkMedit.numEdits; ldel++)
										{
											augmentedLinkMedit.pos[ldel-1] = augmentedLinkMedit.pos[ldel];
											augmentedLinkMedit.edits[ldel-1] = augmentedLinkMedit.edits[ldel];
										}								
										augmentedLinkMedit.numEdits--;
									}
									//else do nothing, sinec the original edit will stay
								}
								else
								{
									//Shift and insert
									for(unsigned char lend = augmentedLinkMedit.numEdits; lend > lind; lend--)
									{
										augmentedLinkMedit.pos[lend] = augmentedLinkMedit.pos[lend-1];
										augmentedLinkMedit.edits[lend] = augmentedLinkMedit.edits[lend-1];
									}
									augmentedLinkMedit.pos[lind] = editOffsetToAdd;
									//the character to add
									augmentedLinkMedit.edits[lind] = editCharToAddFromHomClass;
									augmentedLinkMedit.numEdits++;
								}
								actionTakenFlag = 1;						
								break;
								//(2) readLink's editChar is different than homItems's character, in which case (the editCharToAddFromHomClass is irrelevant, nothing changes)
								//This should already handle the case of deletions overlapping with mismatches (in which case deletion should be kept and mismatch should be ignored)
							}
							//else do nothing, let the loop continue to the next index within qugmented link medit
						}
						//couldn't find any place to insert the emHomItem's edit, so append to the end
						if(!actionTakenFlag)
						{
							augmentedLinkMedit.pos[augmentedLinkMedit.numEdits] = editOffsetToAdd;
							augmentedLinkMedit.edits[augmentedLinkMedit.numEdits] = editCharToAddFromHomClass;
							augmentedLinkMedit.numEdits++;
						}
					}
				}

				//Additional processing in the case of indels
				if(properReadLen != readLen)
				{
					if(properReadLen > readLen) //This means that due to some deletions, the medit we should cover (properReadLen - readLen) more bases
					{
						//add additional potential mismatches here
						if(curEmHomItem.dir == 0)
						{
							//Look at additional bases in the ]readLen, properReadLen] range and add any mismatches from the target. 
							for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
							{
								char charToCompare; //This is the link ref character to compare curEmHomItem's edit to
								if(directionFromLinkToEquivRep == 0)
								{
									charToCompare = fullRef[linkChrCode][linkChrPos + nucInd];
								}
								else
								{
									charToCompare = RevCompChar[(unsigned char) fullRef[linkChrCode][linkChrPos + properReadLen - nucInd - 1]];
								}

								if(charToCompare != fullRef[curEmHomItem.chrCode][curEmHomItem_offsetChrPos + nucInd])
								{
									if(augmentedLinkMedit.pos[augmentedLinkMedit.numEdits-1] <= nucInd)
									{

										augmentedLinkMedit.pos[augmentedLinkMedit.numEdits] = nucInd+1;
										augmentedLinkMedit.edits[augmentedLinkMedit.numEdits] = charToCompare;
										augmentedLinkMedit.numEdits++;
									}
									else
									{
										int insertionPos = augmentedLinkMedit.numEdits-1;
										while(augmentedLinkMedit.pos[insertionPos] > nucInd + 1)
										{
											augmentedLinkMedit.pos[insertionPos+1] = augmentedLinkMedit.pos[insertionPos];
											augmentedLinkMedit.edits[insertionPos+1] = augmentedLinkMedit.edits[insertionPos];
											insertionPos--;
										}
										
										if(augmentedLinkMedit.pos[insertionPos] == nucInd + 1)
										{
											if(fullRef[curEmHomItem.chrCode][curEmHomItem_offsetChrPos + nucInd] == augmentedLinkMedit.edits[insertionPos])
											{
												//Delete the edit since it negates it -- then shift the remaining edits back
												augmentedLinkMedit.numEdits--;
												for(int i=insertionPos; i<augmentedLinkMedit.numEdits; i++)
												{
													augmentedLinkMedit.pos[insertionPos] = augmentedLinkMedit.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
													augmentedLinkMedit.edits[insertionPos] = augmentedLinkMedit.edits[insertionPos+2];
												}
											}
											//else do nothing since the same edit would be valid							
										}
										else
										{
											augmentedLinkMedit.pos[insertionPos] = nucInd + 1;
											augmentedLinkMedit.edits[insertionPos] = charToCompare;
											augmentedLinkMedit.numEdits++;
										}
									}
								}	
							}
						}
						else
						{
							for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
							{
								char charToCompare; //This is the link ref character to compare curEmHomItem's edit to
								if(directionFromLinkToEquivRep == 0)
								{
									charToCompare = fullRef[linkChrCode][linkChrPos + nucInd];
								}
								else
								{
									charToCompare = RevCompChar[(unsigned char) fullRef[linkChrCode][linkChrPos + properReadLen - nucInd - 1]];
								}

								if(charToCompare != RevCompChar[(unsigned char) fullRef[curEmHomItem.chrCode][curEmHomItem_offsetChrPos + readLen - nucInd - 1]])
								{
									if(augmentedLinkMedit.pos[augmentedLinkMedit.numEdits-1] <= nucInd)
									{

										augmentedLinkMedit.pos[augmentedLinkMedit.numEdits] = nucInd+1;
										augmentedLinkMedit.edits[augmentedLinkMedit.numEdits] = charToCompare;
										augmentedLinkMedit.numEdits++;
									}
									else
									{
										int insertionPos = augmentedLinkMedit.numEdits-1;
										while(augmentedLinkMedit.pos[insertionPos] > nucInd + 1)
										{
											augmentedLinkMedit.pos[insertionPos+1] = augmentedLinkMedit.pos[insertionPos];
											augmentedLinkMedit.edits[insertionPos+1] = augmentedLinkMedit.edits[insertionPos];
											insertionPos--;
										}
										
										if(augmentedLinkMedit.pos[insertionPos] == nucInd + 1)
										{
											if(RevCompChar[(unsigned char) fullRef[curEmHomItem.chrCode][curEmHomItem_offsetChrPos + readLen - nucInd - 1]] == augmentedLinkMedit.edits[insertionPos])
											{
												//Delete the edit since it negates it -- then shift the remaining edits back
												augmentedLinkMedit.numEdits--;
												for(int i=insertionPos; i<augmentedLinkMedit.numEdits; i++)
												{
													augmentedLinkMedit.pos[insertionPos] = augmentedLinkMedit.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
													augmentedLinkMedit.edits[insertionPos] = augmentedLinkMedit.edits[insertionPos+2];
												}
											}
											//else do nothing since the same edit would be valid
										}
										else
										{
											augmentedLinkMedit.pos[insertionPos] = nucInd + 1;
											augmentedLinkMedit.edits[insertionPos] = charToCompare;
											augmentedLinkMedit.numEdits++;
										}
									}
								}
							}
						}
					}		
					else //This means that due to some insertions, the mismatches at the end [for the last (properReadLen - readLen) bases should be dropped
					{
						while(augmentedLinkMedit.pos[augmentedLinkMedit.numEdits-1] > properReadLen)
						{
							augmentedLinkMedit.numEdits--;
						}
					}
				}
		
				if(augmentedLinkMedit.numEdits > numMismatchesPerReadMer) //To many errors in one split
				{
					continue;
				}

				if(curEmHomItem.dir == 0)
				{
					travArr[arrNo].push_back(mergedLink(curEmHomItem_offsetChrPos, augmentedLinkMedit, curEmHomItem.chrCode, correctedLinkDir));
				}
				else
				{
					//ReverseComplement medit again
					doubleMedit augmentedLinkMedit_copy = augmentedLinkMedit;
					ConvertToRevCompMedit(augmentedLinkMedit_copy, properReadLen);
					travArr[arrNo].push_back(mergedLink(curEmHomItem_offsetChrPos + (readLen - properReadLen), augmentedLinkMedit_copy, curEmHomItem.chrCode, !correctedLinkDir));
				}
				
				unsigned int homItemEquivIndex = 0;
                                std::tr1::unordered_map<unsigned int, unsigned int>::iterator ite = equivHub[curEmHomItem.chrCode].find(curEmHomItem_offsetChrPos);
                                if(ite != equivHub[curEmHomItem.chrCode].end())
                                {
                                        homItemEquivIndex = ite->second;
                                }

				if(homItemEquivIndex != 0)
				{
					//The representative is already printed
					//First identify the representative of the EquivClass
					short equivBlockOffset = 0; //this is for offseting from block-compressed representation
				        bool dirFromHitToRep = 0; //The direction of current position from the representative;
					RetrieveOffsetAndDirInfo(curEmHomItem.chrCode, curEmHomItem_offsetChrPos, equivBlockOffset, dirFromHitToRep);

					doubleMedit linkMeditForEquiv = augmentedLinkMedit;
					if(curEmHomItem.dir != 0)
					{
						ConvertToRevCompMedit(linkMeditForEquiv, properReadLen);	
					}	
				
					equivClassNode& curEqClass = equivClassList[homItemEquivIndex];

					for(unsigned int listk=0; listk<curEqClass.listSize; listk++)
					{
						eqItem& curEqItem = curEqClass.list[listk];

						if(properReadLen == readLen)
						{	
							if(curEqItem.dir  == 0)
							{
								travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset, linkMeditForEquiv, curEqItem.chrCode, (correctedLinkDir != curEmHomItem.dir)));
							}
							else
							{
								doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
								ConvertToRevCompMedit(linkMeditForEquiv_copy, properReadLen);
								travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset, linkMeditForEquiv_copy, curEqItem.chrCode, !(correctedLinkDir != curEmHomItem.dir)));
							}
						}
						else
						{
							if(properReadLen > readLen)//Deletion Case
							{
								if(curEqItem.dir == 0)
								{
									if(directionFromLinkToEquivRep == 0)
									{
										if(curEmHomItem.dir == 0)
										{
											//Position is not modified but medits might be
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
											
											for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
											{
												char linkChar = fullRef[links[readId].chrCode][links[readId].chrPos + nucInd];

												if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + nucInd] != linkChar)
												{
													if(linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits-1] <= nucInd)
													{
														linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits] = nucInd+1;
														linkMeditForEquiv_copy.edits[linkMeditForEquiv_copy.numEdits] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = linkMeditForEquiv_copy.numEdits-1;
														while(linkMeditForEquiv_copy.pos[insertionPos] > nucInd + 1)
														{
															linkMeditForEquiv_copy.pos[insertionPos+1] = linkMeditForEquiv_copy.pos[insertionPos];
															linkMeditForEquiv_copy.edits[insertionPos+1] = linkMeditForEquiv_copy.edits[insertionPos];
															insertionPos--;
														}
														
														if(linkMeditForEquiv_copy.pos[insertionPos] == nucInd + 1)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + nucInd] == linkMeditForEquiv_copy.edits[insertionPos])
															{
																//Delete the edit since it negates it -- then shift the remaining edits back
																linkMeditForEquiv_copy.numEdits--;
																for(int i=insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+2];
																}
															}
															//else do nothing since the same edit would be valid
														}
														else
														{
															linkMeditForEquiv_copy.pos[insertionPos] = nucInd + 1;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
											
											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset, linkMeditForEquiv_copy, curEqItem.chrCode, correctedLinkDir));
											}
										}
										else
										{
											//Position is not modified but medits might be
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
											
											for(int nucInd = readLen; nucInd < properReadLen; nucInd++)
											{
												char linkChar = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]];
												if(linkChar != fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + readLen - nucInd -1]) // readLen is intentional (+/-properReadLens cancel each other out)
												{

													if(linkMeditForEquiv_copy.pos[0] > properReadLen - nucInd)
													{
														for(int i=linkMeditForEquiv_copy.numEdits; i>0; i--)
														{
															linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
															linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
														}
														linkMeditForEquiv_copy.pos[0] = properReadLen - nucInd;
														linkMeditForEquiv_copy.edits[0] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = 0;
														while(linkMeditForEquiv_copy.pos[insertionPos] < properReadLen - nucInd)
															insertionPos++;

														if(linkMeditForEquiv_copy.pos[insertionPos] == properReadLen - nucInd)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + readLen - nucInd -1] == RevCompChar[(unsigned char) linkMeditForEquiv_copy.edits[insertionPos]])
															{
																//Delete the edit
																linkMeditForEquiv_copy.numEdits--;
																for(int i = insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+1];
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+1];
																}
																//else do nothing
															}
														}
														else
														{
															for(int i=linkMeditForEquiv_copy.numEdits; i>insertionPos; i--)
															{
																linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
																linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
															}
															linkMeditForEquiv_copy.pos[insertionPos] = properReadLen - nucInd;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;	
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
											
											int pushOffset = -properReadLen + readLen;
											
											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, !correctedLinkDir));
											}
										}
									}
									else
									{
										if(curEmHomItem.dir == 0)
										{
											//Position is not modified but medits can be
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;

											int localProperReadLen = properReadLen;
											for(int nucInd = readLen; nucInd < localProperReadLen; nucInd++)
											{
												char linkChar = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + localProperReadLen - nucInd - 1]];
												
												if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + nucInd] != linkChar)
												{
													if(linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits-1] <= nucInd)
													{
														linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits] = nucInd+1;
														linkMeditForEquiv_copy.edits[linkMeditForEquiv_copy.numEdits] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = linkMeditForEquiv_copy.numEdits-1;
														while(linkMeditForEquiv_copy.pos[insertionPos] > nucInd + 1)
														{
															linkMeditForEquiv_copy.pos[insertionPos+1] = linkMeditForEquiv_copy.pos[insertionPos];
															linkMeditForEquiv_copy.edits[insertionPos+1] = linkMeditForEquiv_copy.edits[insertionPos];
															insertionPos--;
														}
														
														if(linkMeditForEquiv_copy.pos[insertionPos] == nucInd + 1)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + nucInd] == linkMeditForEquiv_copy.edits[insertionPos])
															{
																//Delete the edit since it negates it -- then shift the remaining edits back
																linkMeditForEquiv_copy.numEdits--;
																for(int i=insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+2];
																}
															}
															//else do nothing since the same edit would be valid
														}
														else
														{
															linkMeditForEquiv_copy.pos[insertionPos] = nucInd + 1;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
										
											int pushOffset = 0;
											if(curEmHomItem.dir)
												pushOffset = -localProperReadLen + readLen;
											
											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, (correctedLinkDir != curEmHomItem.dir)));
											}
										}
										else
										{
											//Position is not modified but medits might be
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;

											int localProperReadLen = properReadLen;
											for(int nucInd = readLen; nucInd < localProperReadLen; nucInd++)
											{
												char linkChar = fullRef[links[readId].chrCode][links[readId].chrPos + localProperReadLen - nucInd - 1];
												if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + readLen - nucInd - 1] != linkChar)
												{
													if(linkMeditForEquiv_copy.pos[0] > properReadLen - nucInd)
													{
														for(int i=linkMeditForEquiv_copy.numEdits; i>0; i--)
														{
															linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
															linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
														}
														linkMeditForEquiv_copy.pos[0] = properReadLen - nucInd;
														linkMeditForEquiv_copy.edits[0] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = 0;
														while(linkMeditForEquiv_copy.pos[insertionPos] < properReadLen - nucInd)
															insertionPos++;

														if(linkMeditForEquiv_copy.pos[insertionPos] == properReadLen - nucInd)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + readLen - nucInd -1] == RevCompChar[(unsigned char) linkMeditForEquiv_copy.edits[insertionPos]])
															{
																//Delete the edit
																linkMeditForEquiv_copy.numEdits--;
																for(int i = insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+1];
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+1];
																}
																//else do nothing
															}
														}
														else
														{
															for(int i=linkMeditForEquiv_copy.numEdits; i>insertionPos; i--)
															{
																linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
																linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
															}
															linkMeditForEquiv_copy.pos[insertionPos] = properReadLen - nucInd;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;	
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
											
											int pushOffset = 0;
											if(curEmHomItem.dir)
												pushOffset = -localProperReadLen + readLen;
											
											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, (correctedLinkDir != curEmHomItem.dir)));
											}
										}
									}
								}
								else
								{
									if(directionFromLinkToEquivRep == 0)
									{
										if(curEmHomItem.dir == 0)
										{
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
											ConvertToRevCompMedit(linkMeditForEquiv_copy, properReadLen);
											int localProperReadLen = properReadLen;
											for(int nucInd = readLen; nucInd < localProperReadLen; nucInd++)
											{
												char linkChar = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + nucInd]];
												
												if( fullRef[curEqItem.chrCode][curEqItem.chrPos - equivBlockOffset + readLen - nucInd -1] != linkChar) // +/- properReadLens cancel out
												{
													if(linkMeditForEquiv_copy.pos[0] > properReadLen - nucInd)
													{
														for(int i=linkMeditForEquiv_copy.numEdits; i>0; i--)
														{
															linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
															linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
														}
														linkMeditForEquiv_copy.pos[0] = properReadLen - nucInd;
														linkMeditForEquiv_copy.edits[0] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = 0;
														while(linkMeditForEquiv_copy.pos[insertionPos] < properReadLen - nucInd)
															insertionPos++;

														if(linkMeditForEquiv_copy.pos[insertionPos] == properReadLen - nucInd)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + readLen - nucInd -1] == RevCompChar[(unsigned char) linkMeditForEquiv_copy.edits[insertionPos]])
															{
																//Delete the edit
																linkMeditForEquiv_copy.numEdits--;
																for(int i = insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+1];
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+1];
																}
																//else do nothing
															}
														}
														else
														{
															for(int i=linkMeditForEquiv_copy.numEdits; i>insertionPos; i--)
															{
																linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
																linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
															}
															linkMeditForEquiv_copy.pos[insertionPos] = properReadLen - nucInd;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;	
															linkMeditForEquiv_copy.numEdits++;
														}
													}
													

												}
											}
																					
											int pushOffset = 0;
											if(curEmHomItem.dir == 0)
												pushOffset = -localProperReadLen + readLen;

											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, (correctedLinkDir == curEmHomItem.dir)));
											}
										}
										else
										{
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
											ConvertToRevCompMedit(linkMeditForEquiv_copy, properReadLen);
											int localProperReadLen = properReadLen;
											
											for(int nucInd = readLen; nucInd < localProperReadLen; nucInd++)
											{
												char linkChar = fullRef[links[readId].chrCode][links[readId].chrPos + nucInd];
			
												if(fullRef[curEqItem.chrCode][curEqItem.chrPos - equivBlockOffset + nucInd] != linkChar)
												{
													if(linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits-1] <= nucInd)
													{
														linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits] = nucInd+1;
														linkMeditForEquiv_copy.edits[linkMeditForEquiv_copy.numEdits] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = linkMeditForEquiv_copy.numEdits-1;
														while(linkMeditForEquiv_copy.pos[insertionPos] > nucInd + 1)
														{
															linkMeditForEquiv_copy.pos[insertionPos+1] = linkMeditForEquiv_copy.pos[insertionPos];
															linkMeditForEquiv_copy.edits[insertionPos+1] = linkMeditForEquiv_copy.edits[insertionPos];
															insertionPos--;
														}
														
														if(linkMeditForEquiv_copy.pos[insertionPos] == nucInd + 1)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + nucInd] == linkMeditForEquiv_copy.edits[insertionPos])
															{
																//Delete the edit since it negates it -- then shift the remaining edits back
																linkMeditForEquiv_copy.numEdits--;
																for(int i=insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+2];
																}
															}
															//else do nothing since the same edit would be valid
														}
														else
														{
															linkMeditForEquiv_copy.pos[insertionPos] = nucInd + 1;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
										
											int pushOffset = 0;
											if(curEmHomItem.dir == 0)
												pushOffset = -localProperReadLen + readLen;

											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, (correctedLinkDir == curEmHomItem.dir)));
											}
										}
									}
									else
									{
										if(curEmHomItem.dir == 0)
										{
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
											ConvertToRevCompMedit(linkMeditForEquiv_copy, properReadLen);
											
											int localProperReadLen = properReadLen;

											for(int nucInd = readLen; nucInd < localProperReadLen; nucInd++)
											{
												char linkChar = fullRef[links[readId].chrCode][links[readId].chrPos + localProperReadLen - nucInd - 1];

												if(fullRef[curEqItem.chrCode][curEqItem.chrPos - equivBlockOffset + readLen - nucInd -1] != linkChar)
												{	
													if(linkMeditForEquiv_copy.pos[0] > properReadLen - nucInd)
													{
														for(int i=linkMeditForEquiv_copy.numEdits; i>0; i--)
														{
															linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
															linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
														}
														linkMeditForEquiv_copy.pos[0] = properReadLen - nucInd;
														linkMeditForEquiv_copy.edits[0] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = 0;
														while(linkMeditForEquiv_copy.pos[insertionPos] < properReadLen - nucInd)
															insertionPos++;

														if(linkMeditForEquiv_copy.pos[insertionPos] == properReadLen - nucInd)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + readLen - nucInd -1] == linkMeditForEquiv_copy.edits[insertionPos])
															{
																//Delete the edit
																linkMeditForEquiv_copy.numEdits--;
																for(int i = insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+1];
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+1];
																}
																//else do nothing
															}
														}
														else
														{
															for(int i=linkMeditForEquiv_copy.numEdits; i>insertionPos; i--)
															{
																linkMeditForEquiv_copy.pos[i] = linkMeditForEquiv_copy.pos[i-1];
																linkMeditForEquiv_copy.edits[i] = linkMeditForEquiv_copy.edits[i-1];
															}
															linkMeditForEquiv_copy.pos[insertionPos] = properReadLen - nucInd;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;	
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
											
											int pushOffset = 0;
											if(curEmHomItem.dir == 0)
												pushOffset = -localProperReadLen + readLen;

											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, (correctedLinkDir == curEmHomItem.dir)));
											}
										}
										else
										{
											doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
											ConvertToRevCompMedit(linkMeditForEquiv_copy, properReadLen);

											int localProperReadLen = properReadLen;

											for(int nucInd = readLen; nucInd < localProperReadLen; nucInd++)
											{
												char linkChar = RevCompChar[(unsigned char) fullRef[links[readId].chrCode][links[readId].chrPos + localProperReadLen - nucInd - 1]];
												if(fullRef[curEqItem.chrCode][curEqItem.chrPos - equivBlockOffset + nucInd] != linkChar)
												{
													if(linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits-1] <= nucInd)
													{
														linkMeditForEquiv_copy.pos[linkMeditForEquiv_copy.numEdits] = nucInd+1;
														linkMeditForEquiv_copy.edits[linkMeditForEquiv_copy.numEdits] = linkChar;
														linkMeditForEquiv_copy.numEdits++;
													}
													else
													{
														int insertionPos = linkMeditForEquiv_copy.numEdits-1;
														while(linkMeditForEquiv_copy.pos[insertionPos] > nucInd + 1)
														{
															linkMeditForEquiv_copy.pos[insertionPos+1] = linkMeditForEquiv_copy.pos[insertionPos];
															linkMeditForEquiv_copy.edits[insertionPos+1] = linkMeditForEquiv_copy.edits[insertionPos];
															insertionPos--;
														}
														
														if(linkMeditForEquiv_copy.pos[insertionPos] == nucInd + 1)
														{
															if(fullRef[curEqItem.chrCode][curEqItem.chrPos + equivBlockOffset + nucInd] == linkMeditForEquiv_copy.edits[insertionPos])
															{
																//Delete the edit since it negates it -- then shift the remaining edits back
																linkMeditForEquiv_copy.numEdits--;
																for(int i=insertionPos; i<linkMeditForEquiv_copy.numEdits; i++)
																{
																	linkMeditForEquiv_copy.pos[insertionPos] = linkMeditForEquiv_copy.pos[insertionPos+2]; //Since we were pre-shifting, we now have to revert two steps back
																	linkMeditForEquiv_copy.edits[insertionPos] = linkMeditForEquiv_copy.edits[insertionPos+2];
																}
															}
															//else do nothing since the same edit would be valid
														}
														else
														{
															linkMeditForEquiv_copy.pos[insertionPos] = nucInd + 1;
															linkMeditForEquiv_copy.edits[insertionPos] = linkChar;
															linkMeditForEquiv_copy.numEdits++;
														}
													}
												}
											}
																					
											int pushOffset = 0;
											if(curEmHomItem.dir == 0)
												pushOffset = -localProperReadLen + readLen;

											if(linkMeditForEquiv_copy.numEdits <= numMismatchesPerReadMer)
											{
												travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, (correctedLinkDir == curEmHomItem.dir)));
											}
										}
									}
								}
							}
							else //Insertion Case
							{
								if(curEqItem.dir == 0)
								{
									int pushOffset = 0;
									if(curEmHomItem.dir == 1)
										pushOffset = -properReadLen + readLen;
									
									//No need for any modification
									travArr[arrNo].push_back(mergedLink(curEqItem.chrPos + equivBlockOffset + pushOffset, linkMeditForEquiv, curEqItem.chrCode, (correctedLinkDir != curEmHomItem.dir)));
								}
								else
								{
									int pushOffset = 0;
									if(curEmHomItem.dir == 0)
										pushOffset = -properReadLen + readLen;

									doubleMedit linkMeditForEquiv_copy = linkMeditForEquiv;
									ConvertToRevCompMedit(linkMeditForEquiv_copy, properReadLen);

									//Position is modified but no need to change medits
									travArr[arrNo].push_back(mergedLink(curEqItem.chrPos - equivBlockOffset + pushOffset, linkMeditForEquiv_copy, curEqItem.chrCode, correctedLinkDir == curEmHomItem.dir));
								}
							}
						}
					}
				}
			}
			else
			{
				break;
			}
		}
	}
}

//Find all relevant maps of the current read-mer
void ConstructMaps(unsigned int readId, unsigned char arrNo)
{
	if(links[readId].chrCode == 0) //there is no link for sequence (though there is a sequence associated with it in the nohits array)
	{
		return;	
	}

	//This part is for exact homology table traversal
	short equivBlockOffset = 0; //this is for offseting from block-compressed representation
	bool dirFromHitToRep = 0; //The direction of current position from the representative;

	unsigned int curEqClassIndex = RetrieveOffsetAndDirInfo(links[readId].chrCode, links[readId].chrPos, equivBlockOffset, dirFromHitToRep);

	unsigned char repChrCode;
	unsigned int repChrPos;
	if(curEqClassIndex == 0)
	{
		repChrCode = links[readId].chrCode;
		repChrPos = links[readId].chrPos;
	}
	else
	{
		repChrCode = equivClassList[curEqClassIndex].chrCode;
		repChrPos = equivClassList[curEqClassIndex].chrPos; //offsetting is done inside the output function
	}
	TraverseEquivClass(readId, repChrCode, repChrPos, curEqClassIndex, equivBlockOffset, dirFromHitToRep, arrNo);

	//This part is inexact homology table traversal
	TraverseInexactHomologyTable(readId, arrNo, curEqClassIndex, equivBlockOffset, dirFromHitToRep);

	sort(travArr[arrNo].begin(), travArr[arrNo].end(), mergedLinkObject); 
}

//Find all relevant maps of the current read-mer
void ConstructMaps_NoSort(unsigned int readId, unsigned char arrNo)
{
	if(links[readId].chrCode == 0) //there is no link for sequence (though there is a sequence associated with it in the nohits array)
	{
		return;	
	}

	//This part is for exact homology table traversal
	short equivBlockOffset = 0; //this is for offseting from block-compressed representation
	bool dirFromHitToRep = 0; //The direction of current position from the representative;

	unsigned int curEqClassIndex = RetrieveOffsetAndDirInfo(links[readId].chrCode, links[readId].chrPos, equivBlockOffset, dirFromHitToRep);

	unsigned char repChrCode;
	unsigned int repChrPos;
	if(curEqClassIndex == 0)
	{
		repChrCode = links[readId].chrCode;
		repChrPos = links[readId].chrPos;
	}
	else
	{
		repChrCode = equivClassList[curEqClassIndex].chrCode;
		repChrPos = equivClassList[curEqClassIndex].chrPos; //offsetting is done inside the output function
	}
	TraverseEquivClass(readId, repChrCode, repChrPos, curEqClassIndex, equivBlockOffset, dirFromHitToRep, arrNo);

	//This part is inexact homology table traversal
	TraverseInexactHomologyTable(readId, arrNo, curEqClassIndex, equivBlockOffset, dirFromHitToRep);
}

//Find all relevant maps of the current read-mer by searching only in the exact hom table
void ConstructMaps_OnlyExact(unsigned int readId, unsigned char arrNo)
{
	if(links[readId].chrCode == 0) //there is no link for sequence (though there is a sequence associated with it in the nohits array)
	{
		return;	
	}

	//This part is for exact homology table traversal
	short equivBlockOffset = 0; //this is for offseting from block-compressed representation
	bool dirFromHitToRep = 0; //The direction of current position from the representative;

	unsigned int curEqClassIndex = RetrieveOffsetAndDirInfo(links[readId].chrCode, links[readId].chrPos, equivBlockOffset, dirFromHitToRep);

	unsigned char repChrCode;
	unsigned int repChrPos;
	if(curEqClassIndex == 0)
	{
		repChrCode = links[readId].chrCode;
		repChrPos = links[readId].chrPos;
	}
	else
	{
		repChrCode = equivClassList[curEqClassIndex].chrCode;
		repChrPos = equivClassList[curEqClassIndex].chrPos; //offsetting is done inside the output function
	}
	TraverseEquivClass(readId, repChrCode, repChrPos, curEqClassIndex, equivBlockOffset, dirFromHitToRep, arrNo);

	sort(travArr[arrNo].begin(), travArr[arrNo].end(), mergedLinkObject); 
}


//Find all relevant maps of the current read-mer by searching only in the exact hom table
void ConstructMaps_OnlyExact_NoSort(unsigned int readId, unsigned char arrNo)
{
	if(links[readId].chrCode == 0) //there is no link for sequence (though there is a sequence associated with it in the nohits array)
	{
		return;	
	}

	//This part is for exact homology table traversal
	short equivBlockOffset = 0; //this is for offseting from block-compressed representation
	bool dirFromHitToRep = 0; //The direction of current position from the representative;

	unsigned int curEqClassIndex = RetrieveOffsetAndDirInfo(links[readId].chrCode, links[readId].chrPos, equivBlockOffset, dirFromHitToRep);

	unsigned char repChrCode;
	unsigned int repChrPos;
	if(curEqClassIndex == 0)
	{
		repChrCode = links[readId].chrCode;
		repChrPos = links[readId].chrPos;
	}
	else
	{
		repChrCode = equivClassList[curEqClassIndex].chrCode;
		repChrPos = equivClassList[curEqClassIndex].chrPos; //offsetting is done inside the output function
	}
	TraverseEquivClass(readId, repChrCode, repChrPos, curEqClassIndex, equivBlockOffset, dirFromHitToRep, arrNo);
}

void MemoizeConstructedMaps(unsigned char arrNo, int memoCode, bool dirFromReadmer) //memoCode is the index of the memoInds 
{
	unsigned curTravArrSize = travArr[arrNo].size();
	memoInds[memoCode].begInd = curTotalMemoSize;
	memoInds[memoCode].endInd = curTotalMemoSize + curTravArrSize - 1;
	memoInds[memoCode].properReadLength = travArrReadmerLength[arrNo]; 

	if(dirFromReadmer == 0)
	{
		for(unsigned int i=0; i < curTravArrSize; i++)
		{
			mergedLink& myMergedLink = travArr[arrNo][i];
			//Construct medit from mergedMedit
			medit myMedit;
			myMedit.numEdits = myMergedLink.edit.numEdits;
			for(int k=0; k<myMedit.numEdits; k++)
			{
				myMedit.edits[k] = myMergedLink.edit.edits[k];
				myMedit.pos[k] = myMergedLink.edit.pos[k];
			}
			memoArr.push_back(link(myMergedLink.chrPos, myMedit, myMergedLink.chrCode, myMergedLink.dir));
		}	
		curTotalMemoSize += curTravArrSize;
	}
	else
	{
		for(unsigned int i=0; i < curTravArrSize; i++)
		{
			mergedLink& myMergedLink = travArr[arrNo][i];
			//Construct medit from mergedMedit
			medit myMedit;
			myMedit.numEdits = myMergedLink.edit.numEdits;
			for(int k=0; k<myMedit.numEdits; k++)
			{
				myMedit.edits[k] = myMergedLink.edit.edits[k];
				myMedit.pos[k] = myMergedLink.edit.pos[k];
			}
		
			memoArr.push_back(link(myMergedLink.chrPos, myMedit, myMergedLink.chrCode, !myMergedLink.dir));
			curTotalMemoSize++;
		}	
	}
}

void RestoreMemoizedBlock(unsigned int readId, unsigned char arrNo)
{
	bool dirFromReadMer = links[readId].dir;	
	int beginInd = memoInds[links[readId].chrPos].begInd;
	int endInd = memoInds[links[readId].chrPos].endInd; //inclusive interval end
	travArrReadmerLength[arrNo] = memoInds[links[readId].chrPos].properReadLength;	
	
	if(dirFromReadMer == 0)
	{
		for(int i=beginInd; i<=endInd; i++)
		{
			travArr[arrNo].push_back(mergedLink(memoArr[i]));
		}
	}
	else
	{
		int travArrSize = 0;
		for(int i=beginInd; i<=endInd; i++)
		{
			travArr[arrNo].push_back(mergedLink(memoArr[i]));
			travArr[arrNo][travArrSize].dir = !memoArr[i].dir;
			travArrSize++;
		}
	}
}

#define FORWARD_ORDER 0
#define REVERSE_ORDER 1

//This function offsets the second dedit by read length and adds to the first medit.
void MergeSplitMeditsIntoFirst(doubleMedit& d1, const doubleMedit& d2, bool mergeDir, int meditLen1, int meditLen2) //second medit is only acessed but first medit might be updated
{
	if(mergeDir == 0)
	{	
		//if mergeDir is 0 (forw), d2 will simply be incremented by readLen and appended to d1
		for(int k=0; k<d2.numEdits; k++)
		{
			d1.pos[d1.numEdits+k] = d2.pos[k] + meditLen1;
			d1.edits[d1.numEdits+k] = d2.edits[k];	
		}
	}
	else
	{
		//if mergeDir is 1 (rev), d1 will need to be inremented and shifted, d2's original positions will be copied to d1
		for(int k=d1.numEdits-1; k>=0; k--)
		{
			d1.pos[d2.numEdits + k] = d1.pos[k] + meditLen2;
			d1.edits[d2.numEdits + k] = d1.edits[k];
		}
		for(int k=0; k<d2.numEdits; k++)
		{
			d1.pos[k] = d2.pos[k];
			d1.edits[k] = d2.edits[k];
		}
	}
	
	d1.numEdits += d2.numEdits;
}

bool firstSeqSplayedFlag, secondSeqSplayedFlag; //Splayed means the sequence is extracted and set up for comparsion purposes (currently naively done as sequence comparison

char firstSplay[MAX_READ_LEN+10];
char secondSplay[MAX_READ_LEN+10];

doubleMedit tempSplayedMedit; //This is the temporary medit obtained from comparing the medit with an extension

//Splay split is obtaining the proper sequence of the missing read from the link and medit 
int SplaySplit(short chrNo, unsigned int chrPos, const doubleMedit& MD, bool dir, char splay[], int properReadLen)
{
	if(chrPos + properReadLen - 1 > chrLens[chrNo])
	{
		return 0;
	}

	if(dir == 0)
	{
		for(int k=0; k<readLen; k++)
		{
			splay[k] = fullRef[chrNo][chrPos+k];
		}
	
		int insertionOffset = 0;
		int deletionOffset = 0;

		for(int k=0; k<MD.numEdits; k++)
		{
			switch(MD.edits[k])
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					splay[MD.pos[k] - 1 + insertionOffset - deletionOffset] = MD.edits[k];
					break;	
				case 'a':
				case 'c':
				case 'g':
				case 't':
				case 'n':
					{
						int insertionBlockLength = 1; //special case for insertion blocks longer than 1
						while(k + insertionBlockLength < MD.numEdits && MD.pos[k+insertionBlockLength] == MD.pos[k] && MD.edits[k + insertionBlockLength] >= 'a' && MD.edits[k + insertionBlockLength] <= 't')
							insertionBlockLength++;

						for(int i=readLen-1; i>=MD.pos[k] + insertionBlockLength - 1 + insertionOffset - deletionOffset; i--) // >= is intentional since edits are 1 based string is 0 based;
							splay[i] = splay[i-insertionBlockLength];

						for(int i = k; i< k + insertionBlockLength; i++)
							splay[MD.pos[i] + i - k - 1 + insertionOffset - deletionOffset] = MD.edits[i] - 'a' + 'A'; 

						k += insertionBlockLength - 1;
						
						insertionOffset += insertionBlockLength;
						break;	
					}
				case 'D':
					for(int i=MD.pos[k] - 1 + insertionOffset - deletionOffset; i<readLen-1; i++)
						splay[i] = splay[i+1];
					
					splay[readLen - 1] = fullRef[chrNo][chrPos + readLen + deletionOffset - insertionOffset]; // It is actually chrPos + (readLen - 1) + 1 but ones cancel out
					deletionOffset++;
					break;

				default:
					cout << "Error UE1: unexpected edit character: '" << MD.edits[k] << "'  (int = " << (int) MD.edits[k] << ")" << endl;
					exit(98);
			}		
		
		}
	}
	else
	{
		for(int k=0; k<readLen; k++)
		{
			splay[readLen-k-1] = RevCompChar[(unsigned char) fullRef[chrNo][chrPos+k]]; 
		}

		int insertionOffset = 0;
		int deletionOffset = 0;

		for(int k=0; k<MD.numEdits; k++)
		{
			switch(MD.edits[k])
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					splay[readLen - MD.pos[k] - insertionOffset + deletionOffset] = RevCompChar[(unsigned char) MD.edits[k]];
					break;	
				case 'a':
				case 'c':
				case 'g':
				case 't':
				case 'n':
					{
						int insertionBlockLength = 1; //special case for insertion blocks longer than 1
						while(k + insertionBlockLength < MD.numEdits && MD.pos[k+insertionBlockLength] == MD.pos[k] &&  MD.edits[k + insertionBlockLength] >= 'a' && MD.edits[k + insertionBlockLength] <= 't')
							insertionBlockLength++;

						for(int i=readLen-1; i>=MD.pos[k] + insertionBlockLength - 1 + insertionOffset - deletionOffset; i--) // >= is intentional since edits are 1 based string is 0 based;
							splay[readLen - i - 1] = splay[readLen - i + insertionBlockLength - 1];

						for(int i = k; i< k + insertionBlockLength; i++)
							splay[readLen - MD.pos[i] - i + k - insertionOffset + deletionOffset] = RevCompChar[(unsigned char) (MD.edits[i] - 'a' + 'A')]; 

						k += insertionBlockLength - 1;
						
						insertionOffset += insertionBlockLength;
						break;	
					}
				case 'D':
					for(int i=MD.pos[k] - 1 + insertionOffset - deletionOffset; i<readLen-1; i++)
						splay[readLen - i - 1] = splay[readLen - i - 2];
					
					splay[0] = RevCompChar[(unsigned char) fullRef[chrNo][chrPos + readLen + deletionOffset - insertionOffset]]; // It is actually chrPos + (readLen - 1) + 1 but ones cancel out
					deletionOffset++;
					break;
				default:
					cout << "Error UE2: unexpected edit character: '" << MD.edits[k] << "'  (int = " << (int) MD.edits[k] << ")" << endl;
					cout << "chrNo: " << (int) chrNo << " chrPos: " << chrPos << endl;
					cout << "MD: " << DebugPrintMedit(MD,444) << endl;
					cout << "dir: " << (int) dir << endl;
					cout << "properReadLen: " << properReadLen << endl;
					exit(98);
			}		
		}
	}

	splay[readLen] = '\0';
	return 1;
}

unsigned char DP[MAX_READ_LEN][MAX_READ_LEN];
char DParrow[MAX_READ_LEN][MAX_READ_LEN];

void InitIndelAlignmentMatrix()
{
	DP[0][0] = 0;
	for(int i=1; i<=MAX_DOUBLE_EDITS; i++)
	{
		DP[i][0] = i;
		DParrow[i][0] = 3;
		DP[0][i] = i;
		DParrow[0][i] = 2;
	}		
}

//TODO(denizy) do path compression for DP arrows for speed.	
//TODO(denizy) Pruning for dead alignments for speed.
bool IndelAlignment(doubleMedit& splayMD, char splay[], char* ref, bool compDir, unsigned char errLimit, bool boundaryDirection, char& modChrPosOffset)
{
	if(compDir == 0)
	{
		if(boundaryDirection == 0)
		{
			for(int i=1; i<=readLen; i++)
			{
				unsigned char minForRow = errLimit + 1;				
				for(int j = max(1, i - errLimit); j <= i+errLimit; j++)
				{
					if(splay[i-1] != ref[j-1])
					{
						DParrow[i][j] = -1; //-1 indicates mismatch
						DP[i][j] = DP[i-1][j-1] + 1;
					}
					else
					{
						DParrow[i][j] = 1;
						DP[i][j] = DP[i-1][j-1];
					}

					if( j > i - errLimit && DP[i][j-1] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i][j-1] + 1;
						DParrow[i][j] = 2; //Deletion from read
					}

					if( j < i + errLimit && DP[i-1][j] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i-1][j] + 1;
						DParrow[i][j] = 3; //Insertion to read;
					}

					minForRow = min(DP[i][j], minForRow);	
				}
	
				if(minForRow > errLimit)
				{
					return 0;
				}
			}
		}
		else
		{
			//Right-genome boundary version
			for(int i=1; i<=readLen; i++)
			{
				unsigned char minForRow = errLimit + 1;				
				for(int j = max(1, i - errLimit); j <= i+errLimit; j++)
				{
					if(splay[readLen-i] != ref[readLen-j])
					{
						DParrow[i][j] = -1; //-1 indicates mismatch
						DP[i][j] = DP[i-1][j-1] + 1;
					}
					else
					{
						DParrow[i][j] = 1;
						DP[i][j] = DP[i-1][j-1];
					}

					if( j > i - errLimit && DP[i][j-1] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i][j-1] + 1;
						DParrow[i][j] = 2; //Deletion from read
					}

					if( j < i + errLimit && DP[i-1][j] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i-1][j] + 1;
						DParrow[i][j] = 3; //Insertion to read;
					}
					
					minForRow = min(DP[i][j], minForRow);	
				}

				if(minForRow > errLimit)
				{
					return 0;
				}
			}
		}
	}
	else
	{
		if(boundaryDirection == 0)
		{
			for(int i=1; i<=readLen; i++)
			{
				unsigned char minForRow = errLimit + 1;				
				for(int j = max(1, i - errLimit); j <= i+errLimit; j++)
				{
					if(RevCompChar[(unsigned char) splay[readLen-i]] != ref[j-1])
					{
						DParrow[i][j] = -1; //-1 indicates mismatch
						DP[i][j] = DP[i-1][j-1] + 1;
					}
					else
					{
						DParrow[i][j] = 1;
						DP[i][j] = DP[i-1][j-1];
					}

					if( j > i - errLimit && DP[i][j-1] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i][j-1] + 1;
						DParrow[i][j] = 2; //Deletion from read
					}

					if( j < i + errLimit && DP[i-1][j] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i-1][j] + 1;
						DParrow[i][j] = 3; //Insertion to read;
					}
					
					minForRow = min(DP[i][j], minForRow);	
				}
				
				if(minForRow > errLimit)
				{
					return 0;
				}
			}
		}
		else
		{
			for(int i=1; i<=readLen; i++)
			{
				unsigned char minForRow = errLimit + 1;				
				for(int j = max(1, i - errLimit); j <= i+errLimit; j++)
				{
					if(RevCompChar[(unsigned char) splay[i-1]] != ref[readLen-j])
					{
						DParrow[i][j] = -1; //-1 indicates mismatch
						DP[i][j] = DP[i-1][j-1] + 1;
					}
					else
					{
						DParrow[i][j] = 1;
						DP[i][j] = DP[i-1][j-1];
					}

					if( j > i - errLimit && DP[i][j-1] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i][j-1] + 1;
						DParrow[i][j] = 2; //Deletion from read
					}

					if( j < i + errLimit && DP[i-1][j] + 1 < DP[i][j])
					{
						DP[i][j] = DP[i-1][j] + 1;
						DParrow[i][j] = 3; //Insertion to read;
					}
					minForRow = min(DP[i][j], minForRow);	
				}

				if(minForRow > errLimit)
				{
					return 0;
				}
				
			}
		}
	}
	
	//Find the best scoring alignment here
	int minLoc = readLen - errLimit;
	for(int j=readLen - errLimit + 1; j <= readLen + errLimit; j++)
		if( DP[readLen][j] < DP[readLen][minLoc])
			minLoc = j;

	if(DP[readLen][minLoc] > errLimit)
	{
		return 0;
	}
	
	//Backtrack DP arrows and assign edits to the splayMD.
	//Note that the number of edits is known from above, so no shifting is required
	splayMD.numEdits = DP[readLen][minLoc];

	if(compDir == 0)
	{
		if(boundaryDirection == 0)
		{
			int curEditPos = splayMD.numEdits - 1;
			int curReadPos = readLen;

			while(curEditPos >= 0)
			{
				switch(DParrow[curReadPos][minLoc])
				{
					case -1:
						splayMD.pos[curEditPos] = minLoc;
						splayMD.edits[curEditPos] = splay[curReadPos-1];
						curEditPos--;
					case 1:
						minLoc--;
						curReadPos--;
						break;
					case 2:
						splayMD.pos[curEditPos] = minLoc;
						splayMD.edits[curEditPos] = 'D';
						curEditPos--; 
						minLoc--;
						break;				
					case 3:
						splayMD.pos[curEditPos] = minLoc + 1;
						splayMD.edits[curEditPos] = splay[curReadPos-1] - 'A' + 'a';
						curReadPos--;
						curEditPos--;
						break;
				}
			}
		}
		else
		{
			modChrPosOffset = char(readLen - minLoc);

			int curEditPos = 0;
			int curRow = readLen;

			while(curEditPos < splayMD.numEdits)
			{
				switch(DParrow[curRow][minLoc])
				{
					case -1:
						splayMD.pos[curEditPos] = readLen - minLoc + 1;
						splayMD.edits[curEditPos] = splay[readLen - curRow];
						curEditPos++;
					case 1:
						minLoc--;
						curRow--;
						break;
					case 2:
						splayMD.pos[curEditPos] = readLen - minLoc + 1;
						splayMD.edits[curEditPos] = 'D';
						curEditPos++; 
						minLoc--;
						break;				
					case 3:
						splayMD.pos[curEditPos] = readLen - minLoc + 1; //Plus +1 due to left shifted formation of the edits in the normal case
						splayMD.edits[curEditPos] = splay[readLen - curRow] - 'A' + 'a';
						curRow--;
						curEditPos++;
						break;
				}
			}
		}
	}
	else
	{
		if(boundaryDirection == 0)
		{
			int curEditPos = splayMD.numEdits - 1;
			int curReadPos = readLen;

			while(curEditPos >= 0)
			{
				switch(DParrow[curReadPos][minLoc])
				{
					case -1:
						splayMD.pos[curEditPos] = minLoc;
						splayMD.edits[curEditPos] = RevCompChar[(unsigned char) splay[readLen-curReadPos]];
						curEditPos--;
					case 1:
						minLoc--;
						curReadPos--;
						break;
					case 2:
						splayMD.pos[curEditPos] = minLoc;
						splayMD.edits[curEditPos] = 'D';
						curEditPos--; 
						minLoc--;
						break;				
					case 3:
						splayMD.pos[curEditPos] = minLoc + 1;
						splayMD.edits[curEditPos] = RevCompChar[(unsigned char) (splay[readLen - curReadPos] - 'A' + 'a')];
						curReadPos--;
						curEditPos--;
						break;
				}
			}
		}
		else
		{
			int curEditPos = 0;
			int curRow = readLen;

			modChrPosOffset = char(readLen - minLoc);

			while(curEditPos < splayMD.numEdits)
			{
				switch(DParrow[curRow][minLoc])
				{
					case -1:
						splayMD.pos[curEditPos] = readLen - minLoc + 1;
						splayMD.edits[curEditPos] = RevCompChar[(unsigned int) splay[curRow-1]];
						curEditPos++;
					case 1:
						minLoc--;
						curRow--;
						break;
					case 2:
						splayMD.pos[curEditPos] = readLen - minLoc + 1;
						splayMD.edits[curEditPos] = 'D';
						curEditPos++; 
						minLoc--;
						break;				
					case 3:
						splayMD.pos[curEditPos] = readLen - minLoc + 1; //-1 + 1
						splayMD.edits[curEditPos] = RevCompChar[(unsigned int) (splay[curRow-1] - 'A' + 'a')];
						curRow--;
						curEditPos++;
						break;
				}
			}
		}
	}
	return 1;
}

//identifies the differences between a reference and a given sequence
//splay with edits in the case that there is no alignment with mismatches
//returns the modified starting chrPos for right bounded alignment
bool CheckSplayError(doubleMedit& splayMD, char splay[], unsigned char refChrNo, unsigned int refChrPos, bool compDir, unsigned char errLimit, bool boundaryDirection, char& modChrPosOffset)  //0 means read is bounded on the left, 1 means ref is bounded on the right
{
	//For indels the reference alignment should be strict on one side (it could be the beginning or the ending), it could be relaxed on the other side
	//The bounded dynamic programming algorithm would start from that side.
	if(refChrPos > chrLens[refChrNo]) //Checks for unsigned int overflow
	{
		return 0;
	} 
	else if(refChrPos + readLen - 1 > chrLens[refChrNo]) //Checks for illegally hanging read ends
	{
		return 0;
	}

	char* ref = fullRef[refChrNo] + refChrPos; 

	int tryIndelFlag = 0;

	if(compDir == 0)
	{
		for(int k=0; k<readLen; k++)
		{
			if(splay[k] != ref[k])
			{
				if(splayMD.numEdits >= errLimit)
				{
					splayMD.numEdits=0;
					tryIndelFlag = 1;
					break;
				}

				splayMD.pos[splayMD.numEdits] = k+1;
				splayMD.edits[splayMD.numEdits] = splay[k]; //note that medits store differences in the read (where as sam MD will store diffs in ref)
				splayMD.numEdits++;				
			}			
		}
	}
	else
	{
		for(int k=0; k<readLen; k++)
		{
			if(RevCompChar[(unsigned char) splay[readLen-k-1]] != ref[k])
			{
				if(splayMD.numEdits >= errLimit)
				{
					splayMD.numEdits=0;
					tryIndelFlag = 1;
					break;
				}

				splayMD.pos[splayMD.numEdits] = k+1;
				splayMD.edits[splayMD.numEdits] = RevCompChar[(unsigned char) splay[readLen-k-1]]; // This are oriented in forwRef
				splayMD.numEdits++;
			}
		}
	}

	if(tryIndelFlag)
	{
		// Hamming distance alignment prevents IndelAlignment
		// Same applies to BEST_FAST Levenshtein mapping mode
		if(distanceMetric == HAMMING || mappingMode == BEST_FAST_MAPPING_MODE)
		{
			return 0;
		}

		bool alignmentFlag = IndelAlignment(splayMD, splay, ref, compDir, errLimit, boundaryDirection, modChrPosOffset);
		return alignmentFlag;
	}
	return 1;
}			

//Identifies if any noHit access will be needed and sets up proper flags for it
void MarkVectorForMissingSplits(unsigned char sourceNo, unsigned char targetNo)
{
	bool splitType = sourceNo % 2; // 0,2-> becomes 0 (left split) 2,4 becomes 1 (right split)

	//Remove all splits with MAX edits and put offset
	int markedIndex = 0;
	unsigned sourceSize = travArr[sourceNo].size();
	for(unsigned int k = 0; k < sourceSize; k++)
	{
		bool dir = travArr[sourceNo][k].dir;
		if(travArr[sourceNo][k].edit.numEdits < numMismatchesPerReadMer)
		{
			travArr[targetNo].push_back(travArr[sourceNo][k]);
			if( dir xor splitType ) //if forward direction and second split, or reverse direction first split 
			{
				if(travArr[targetNo][markedIndex].chrPos <= readLen)
				{
					travArr[targetNo].pop_back();
					//Don't modify marked index here since we popped back what we pushed back
					continue;
				}
				else
				{
					//increment all edits by readLen and decrement chrPos by readLen
					travArr[targetNo][markedIndex].chrPos -= readLen;
					for(int j=0; j < travArr[targetNo][markedIndex].edit.numEdits; j++)
					{
						travArr[targetNo][markedIndex].edit.pos[j] += readLen;
					}
				}
			}
				
			if(splitType) // we have is the first split (this is important, since it doesn't have to be the same as travArr[arrNo1]
			{
				travArr[targetNo][markedIndex].edit.numEdits += MISSING_SPLIT_OFFSET; //Since the second split is missing
			}
			else //split we have is the second split
			{
				travArr[targetNo][markedIndex].edit.numEdits += DOUBLE_MISSING_SPLIT_OFFSET; //Since the first split is missing
			}
			markedIndex++;
		}
	}
}

//After two split's maps are independently constructed, this merges the ones that are in proper offset
//Merging is aware of the length of each readLen alignment after indel length modification
void MergeSplitMaps(unsigned char arrNo1, unsigned char arrNo2)
{
	unsigned char mergedArrNo = (arrNo1/2)+4; // meaning 0->4 3->5

	unsigned int vec1size = travArr[arrNo1].size();
	unsigned int vec2size = travArr[arrNo2].size();

	firstSeqSplayedFlag = 0;
        secondSeqSplayedFlag = 0;

	if(vec1size == 0)
	{
		MarkVectorForMissingSplits(arrNo2, mergedArrNo);
	}
	else if(vec2size == 0)
	{
		MarkVectorForMissingSplits(arrNo1, mergedArrNo);
	}
	else //both mates exist
	{
		//These are redefined due to overwriting that might occur below
		short firstChrCode = travArr[arrNo1][0].chrCode;
        	unsigned int firstChrPos = travArr[arrNo1][0].chrPos;
        	doubleMedit firstMedit = travArr[arrNo1][0].edit;
        	bool firstDir = travArr[arrNo1][0].dir;

		unsigned int posSplit1 = 0, posSplit2_forw = 0, posSplit2_rev = 0;
		
		while(posSplit1 < vec1size)
		{
			bool matchFoundForFirst = 0;
			bool split1dir = travArr[arrNo1][posSplit1].dir;


			if(split1dir == 0)
			{
				while(posSplit2_forw < vec2size) 
				{
					if(travArr[arrNo2][posSplit2_forw].dir == 0)
					{
						if(travArr[arrNo2][posSplit2_forw].chrCode == travArr[arrNo1][posSplit1].chrCode)
						{
							if(travArr[arrNo2][posSplit2_forw].chrPos == travArr[arrNo1][posSplit1].chrPos + travArrReadmerLength[arrNo1]) //This is the offset to be observed
							{
								matchFoundForFirst = 1;
								MergeSplitMeditsIntoFirst(travArr[arrNo1][posSplit1].edit, travArr[arrNo2][posSplit2_forw].edit, split1dir, travArrReadmerLength[arrNo1], travArrReadmerLength[arrNo2]); //0 is forward flag
								travArr[mergedArrNo].push_back(travArr[arrNo1][posSplit1]); 
								//no need to update the chrPos here since it will hold the same value from pos1 
								posSplit2_forw++;
								break;
							}
							else if(travArr[arrNo2][posSplit2_forw].chrPos > travArr[arrNo1][posSplit1].chrPos + travArrReadmerLength[arrNo1]) //need to increment pos1 first
							{
								break;
							}
						}
						else if(travArr[arrNo2][posSplit2_forw].chrCode > travArr[arrNo1][posSplit1].chrCode)
						{
							break;
						}

						if(travArr[arrNo2][posSplit2_forw].edit.numEdits < numMismatchesPerReadMer)
						{
							int didSplay = 1;
							if(!firstSeqSplayedFlag)
							{
								didSplay = SplaySplit(firstChrCode, firstChrPos, firstMedit, firstDir, firstSplay, travArrReadmerLength[arrNo1]); //if returns zero, there is an out of bounds error
								if(didSplay)
								{
									firstSeqSplayedFlag = 1;
								}
							}	

							if(didSplay)
							{
								if(travArr[arrNo2][posSplit2_forw].chrPos > (unsigned int) readLen)
								{
									tempSplayedMedit.numEdits = 0;

									char modChrPosOffset = 0;
									bool validFlag = CheckSplayError(tempSplayedMedit, firstSplay, travArr[arrNo2][posSplit2_forw].chrCode, travArr[arrNo2][posSplit2_forw].chrPos - readLen, 0, doubleNumMismatchesPerReadMer-travArr[arrNo2][posSplit2_forw].edit.numEdits, 1, modChrPosOffset); //1 means right-bounded indel alignment
									if(validFlag)
									{
										MergeSplitMeditsIntoFirst(tempSplayedMedit, travArr[arrNo2][posSplit2_forw].edit, 0, readLen, travArrReadmerLength[arrNo2]); //0 is forward flag

										for(int i=0; i<tempSplayedMedit.numEdits; i++)
											tempSplayedMedit.pos[i] -= modChrPosOffset;

										travArr[mergedArrNo].push_back(travArr[arrNo2][posSplit2_forw]);
										unsigned int mergedSize = travArr[mergedArrNo].size();
										travArr[mergedArrNo][mergedSize - 1].chrPos -= readLen - modChrPosOffset;
										travArr[mergedArrNo][mergedSize - 1].edit = tempSplayedMedit;
									}
								}
							}
						}
					}
					posSplit2_forw++;
				}
			}
			else
			{
				while(posSplit2_rev < vec2size)
				{
					if(travArr[arrNo2][posSplit2_rev].dir == 1)
					{
						if(travArr[arrNo2][posSplit2_rev].chrCode == travArr[arrNo1][posSplit1].chrCode)
						{
							if(travArr[arrNo2][posSplit2_rev].chrPos == travArr[arrNo1][posSplit1].chrPos - travArrReadmerLength[arrNo2])
							{
								matchFoundForFirst = 1;
								MergeSplitMeditsIntoFirst(travArr[arrNo1][posSplit1].edit, travArr[arrNo2][posSplit2_rev].edit, split1dir, travArrReadmerLength[arrNo1], travArrReadmerLength[arrNo2]);
								travArr[mergedArrNo].push_back(travArr[arrNo1][posSplit1]);
								travArr[mergedArrNo][travArr[ mergedArrNo].size() - 1].chrPos = travArr[arrNo2][posSplit2_rev].chrPos; //Since it's reverse, the position will now come from the second split
								posSplit2_rev++;
								break;
							}
							else if(travArr[arrNo2][posSplit2_rev].chrPos > travArr[arrNo1][posSplit1].chrPos - travArrReadmerLength[arrNo2]) //need to increment pos1 first
							{
								break;
							}
						}
						else if(travArr[arrNo2][posSplit2_rev].chrCode > travArr[arrNo1][posSplit1].chrCode)
						{
							break;
						}
					
						if(travArr[arrNo2][posSplit2_rev].edit.numEdits < numMismatchesPerReadMer)
						{
							int didSplay = 1;
							if(!firstSeqSplayedFlag)
							{

								didSplay = SplaySplit(firstChrCode, firstChrPos, firstMedit, firstDir, firstSplay, travArrReadmerLength[arrNo1]);
								if(didSplay)
								{
									firstSeqSplayedFlag = 1;
								}
							}		
		
							if(didSplay)
							{
								tempSplayedMedit.numEdits = 0;

								//here readLen got modified to travArrReadmerLength since the starting positions of the firstSplit to splay is determined by the alignment length of the second split
								char modChrPosOffset = 0;
								bool validFlag = CheckSplayError(tempSplayedMedit, firstSplay, travArr[arrNo2][posSplit2_rev].chrCode, travArr[arrNo2][posSplit2_rev].chrPos + travArrReadmerLength[arrNo2], 1, doubleNumMismatchesPerReadMer-travArr[arrNo2][posSplit2_rev].edit.numEdits, 0, modChrPosOffset); //Indel alignment is left-bounded and reverse
							
								if(validFlag)
								{
									MergeSplitMeditsIntoFirst(tempSplayedMedit, travArr[arrNo2][posSplit2_rev].edit, 1, readLen, travArrReadmerLength[arrNo2]); //1 is reverse flag
									travArr[mergedArrNo].push_back(travArr[arrNo2][posSplit2_rev]);
									unsigned int mergedSize = travArr[mergedArrNo].size();
									travArr[mergedArrNo][mergedSize - 1].edit = tempSplayedMedit;
								}
							}
						}
					}
					posSplit2_rev++;
				}
			}
	
			if(!matchFoundForFirst)
			{
				//No match found for posSplit1
				if(travArr[arrNo1][posSplit1].edit.numEdits < numMismatchesPerReadMer) //then there is a chance that the second split has more errors
				{
					int didSplay = 1;
					if(!secondSeqSplayedFlag)
					{
						didSplay = SplaySplit(travArr[arrNo2][0].chrCode, travArr[arrNo2][0].chrPos, travArr[arrNo2][0].edit, travArr[arrNo2][0].dir, secondSplay, travArrReadmerLength[arrNo2]);
						if(didSplay)
						{
							secondSeqSplayedFlag = 1;
						}
					}
						
					if(didSplay)
					{		
						tempSplayedMedit.numEdits = 0; //fast way to reset medit
						bool validFlag = 0;
						char modChrPosOffset = 0;					
	
						//Here readLen is modified to travArrReadmerLength[arrNo1] only in the forward direction case since, in the forward case the starting position of the split to splay is determined by the indel modified alignment length of the first split, not for the reverse direction case
						if(split1dir == 0)
						{
							validFlag = CheckSplayError(tempSplayedMedit, secondSplay, travArr[arrNo1][posSplit1].chrCode, travArr[arrNo1][posSplit1].chrPos + travArrReadmerLength[arrNo1], split1dir, doubleNumMismatchesPerReadMer-travArr[arrNo1][posSplit1].edit.numEdits, 0, modChrPosOffset);
						}
						else
						{
							if(travArr[arrNo1][posSplit1].chrPos > (unsigned int) readLen)
							{
								validFlag = CheckSplayError(tempSplayedMedit, secondSplay, travArr[arrNo1][posSplit1].chrCode, travArr[arrNo1][posSplit1].chrPos - readLen, split1dir, doubleNumMismatchesPerReadMer-travArr[arrNo1][posSplit1].edit.numEdits, 1, modChrPosOffset);	
							}
						}

						if(validFlag)
						{
							MergeSplitMeditsIntoFirst(travArr[arrNo1][posSplit1].edit, tempSplayedMedit, split1dir, travArrReadmerLength[arrNo1], readLen); //0 is forward flag

							travArr[mergedArrNo].push_back(travArr[arrNo1][posSplit1]);
							if(split1dir == 1)
							{
								unsigned int mergedLastIndex = travArr[mergedArrNo].size() - 1;
								
								for(int i=0; i<travArr[mergedArrNo][mergedLastIndex].edit.numEdits; i++)
									travArr[mergedArrNo][mergedLastIndex].edit.pos[i] -= modChrPosOffset;
								travArr[mergedArrNo][mergedLastIndex].chrPos -= readLen - modChrPosOffset; //Since it's reverse, the position will now come from the second split that is upstream
							}
							//no need to update the chrPos here since it will hold the same value from pos1 
						}
					}
				}
			}
			posSplit1++;
		}

		unsigned int posSplit2_last = max(posSplit2_forw, posSplit2_rev);
		//this is to handle remaining second split -> first extensions after first splits are finished
		while(posSplit2_last < vec2size)
		{
			if(travArr[arrNo2][posSplit2_last].dir == 0)
			{
				if(travArr[arrNo2][posSplit2_last].edit.numEdits < numMismatchesPerReadMer)
				{
					int didSplay = 1;
					if(!firstSeqSplayedFlag)
					{
						didSplay = SplaySplit(firstChrCode, firstChrPos, firstMedit, firstDir, firstSplay, travArrReadmerLength[arrNo1]); //if returns zero, there is an out of bounds error
						if(didSplay)
						{
							firstSeqSplayedFlag = 1;
						}
					}	

					if(didSplay)
					{
						if(travArr[arrNo2][posSplit2_last].chrPos > (unsigned int) readLen)
						{
							tempSplayedMedit.numEdits = 0;
							char modChrPosOffset = 0;
							//generated split is compared to the chromosome code of the second
							bool validFlag = CheckSplayError(tempSplayedMedit, firstSplay, travArr[arrNo2][posSplit2_last].chrCode, travArr[arrNo2][posSplit2_last].chrPos - readLen, 0, doubleNumMismatchesPerReadMer-travArr[arrNo2][posSplit2_last].edit.numEdits, 1, modChrPosOffset);

							if(validFlag)
							{
								MergeSplitMeditsIntoFirst(tempSplayedMedit, travArr[arrNo2][posSplit2_last].edit, 0, readLen, travArrReadmerLength[arrNo2]); //0 is forward flag
								
								for(int i=0; i<tempSplayedMedit.numEdits; i++)
									tempSplayedMedit.pos[i] -= modChrPosOffset;
								
								travArr[mergedArrNo].push_back(travArr[arrNo2][posSplit2_last]);
								unsigned int mergedSize = travArr[mergedArrNo].size();
								travArr[mergedArrNo][mergedSize - 1].chrPos -= readLen - modChrPosOffset;
								travArr[mergedArrNo][mergedSize - 1].edit = tempSplayedMedit;
							}
						}
					}
				}
			}
			else
			{
				if(travArr[arrNo2][posSplit2_last].edit.numEdits < numMismatchesPerReadMer)
				{
					int didSplay = 1;
					if(!firstSeqSplayedFlag)
					{
						didSplay = SplaySplit(firstChrCode, firstChrPos, firstMedit, firstDir, firstSplay, travArrReadmerLength[arrNo1]);
						if(didSplay)
						{
							firstSeqSplayedFlag = 1;
						}
					}		

					if(didSplay)
					{
						tempSplayedMedit.numEdits = 0;
						char modChrPosOffset = 0;
						bool validFlag = CheckSplayError(tempSplayedMedit, firstSplay, travArr[arrNo2][posSplit2_last].chrCode, travArr[arrNo2][posSplit2_last].chrPos + travArrReadmerLength[arrNo2], 1, doubleNumMismatchesPerReadMer-travArr[arrNo2][posSplit2_last].edit.numEdits, 0, modChrPosOffset);
					
						if(validFlag)
						{
							MergeSplitMeditsIntoFirst(tempSplayedMedit, travArr[arrNo2][posSplit2_last].edit, 1, readLen, travArrReadmerLength[arrNo2]); //1 is reverse flag
							travArr[mergedArrNo].push_back(travArr[arrNo2][posSplit2_last]);
							unsigned int mergedSize = travArr[mergedArrNo].size();
							travArr[mergedArrNo][mergedSize - 1].edit = tempSplayedMedit;
						}
					}
				}
			}
			posSplit2_last++;
		}

		sort(travArr[mergedArrNo].begin(), travArr[mergedArrNo].end(), mergedLinkObject);
	}
}

#define MAX_MAP_COUNT_TIERS 2*MAX_DOUBLE_EDITS
unsigned int mapQcounts[MAX_MAP_COUNT_TIERS+7];

unsigned char CalculateMapQ(unsigned char bestDistance)
{
	//Look at an array that has counts for mappings intercepted for best-mapping [0 to 2*MAX_DOUBLE_EDITS]
        int mapQ = 42 / mapQcounts[bestDistance];
        mapQ -= bestDistance;
        mapQ -= min((unsigned int) 5, mapQcounts[bestDistance+1]/4);
        mapQ -= min((unsigned int) 5, mapQcounts[bestDistance+2]/16);
        mapQ -= min((unsigned int) 5, mapQcounts[bestDistance+3]/64);
        mapQ -= min((unsigned int) 5, mapQcounts[bestDistance+4]/256);

        if(mapQ < 0)
        {
                mapQ = 0;
        }

        return (unsigned char) mapQ;
}

int printBufferSize = 5005000;
int bufferOffset = 0;
char printBuffer[9000000]; //additional space is for overflow avoidance (before flushing)

inline void flushBuffer(FILE* fout)
{
	fwrite(printBuffer, sizeof(char), bufferOffset, fout);
	bufferOffset = 0;
}

char singleItemChrPosStr[20];
char leftItemChrPosStr[20];
char rightItemChrPosStr[20];
char tLenStr_plus[20];
char tLenStr_minus[20];
char printReadIdStr[50];
char tempPrintStr[MAX_READ_LEN+3];

char leftQualityScore[MAX_READ_LEN+3];
char rightQualityScore[MAX_READ_LEN+3];

//The following parameters are for readNameRecovery (and quality score)
bool printReadNamesFlag = 0;
bool printQualityScoresFlag = 0;

unsigned long long previousReadCodeLeft = 0;
unsigned long long previousReadCodeRight = 0;

ifstream finFastqLeft;
ifstream finFastqRight;

char** readFileList_left;
char** readFileList_right;

unsigned long long * numReadsInFile;
unsigned long long * lastReadCodeInFile; //The code for the last read in file

int numFastqs;
int curFastqNo;

bool isReadNameCharInvalid[256];

void InitFastqs(ifstream& finRFL)
{
        isReadNameCharInvalid[' '] = 1;
        isReadNameCharInvalid['\0'] = 1;
        isReadNameCharInvalid['\t'] = 1;
        isReadNameCharInvalid['\n'] = 1;
        isReadNameCharInvalid['/'] = 1;

        finRFL >> numFastqs;

        readFileList_left = (char **) malloc (numFastqs * sizeof(char*));
        readFileList_right = (char **) malloc (numFastqs * sizeof(char*));
        for(int i=0; i<numFastqs; i++)
        {
                readFileList_left[i] = (char *) malloc ((MAX_READ_FILE_NAME_LEN+2) * sizeof(char));
                readFileList_right[i] = (char *) malloc ((MAX_READ_FILE_NAME_LEN+2) * sizeof(char));
        }

        numReadsInFile = (unsigned long long *) malloc (numFastqs * sizeof(unsigned long long));
        lastReadCodeInFile = (unsigned long long *) malloc (numFastqs * sizeof(unsigned long long));

        for(int i=0; i<numFastqs; i++)
        {
                finRFL >> readFileList_left[i];
                if(i==0)
                {
                        finFastqLeft.open(readFileList_left[i]);
                        assert(finFastqLeft.is_open());
                }

                if(inputMode == PAIRED_MODE)
                {
                        finRFL >> readFileList_right[i];

			if(i == 0)
                        {
                                finFastqRight.open(readFileList_right[i]);
                                assert(finFastqRight.is_open());
                        }
                }

                finRFL >> numReadsInFile[i];

                if(i==0)
                        lastReadCodeInFile[i] = numReadsInFile[i];
                else
                        lastReadCodeInFile[i] = lastReadCodeInFile[i-1] + numReadsInFile[i];
        }
}

//This function extracts the read name from leftReadFile and also sets the quality score string for left read (doesn't touch the right read's quality scores) 
void SetCurrentReadName(unsigned long long printReadCode) //this always works on the first readFile
{
	if(printReadCode == previousReadCodeLeft)
	{
		return; //If same as before, it is likely a reprint (e.g. all-mapping)
	}

	//If readCode is more the last item in the current file, skip file
        if(printReadCode > lastReadCodeInFile[curFastqNo])
        {
                //If the read no exceeds the current read file, close the ifstream
                finFastqLeft.clear();
                finFastqLeft.close();

		//This skips multiple files if needed
		while(printReadCode > lastReadCodeInFile[curFastqNo])
		{
                	curFastqNo++;
                }

                //And open the file that contains the current read (it's guaranteed that there won't be any read code larger than the total read size in the files.
		finFastqLeft.open(readFileList_left[curFastqNo]);

		//setting previous read code as a distance marker
		previousReadCodeLeft = lastReadCodeInFile[curFastqNo-1];		

                if(printQualityScoresFlag && inputMode == PAIRED_MODE) //Also arrange the second file as well in the paired-end situation
                {
                        finFastqRight.clear();
                        finFastqRight.close();
                        finFastqRight.open(readFileList_right[curFastqNo]);
                        previousReadCodeRight = lastReadCodeInFile[curFastqNo-1];
                }
        }

	//Easy skip for any read that is not relevant
        int numReadsToSkip = printReadCode - previousReadCodeLeft - 1;
        string junk;
        for(int i=0; i<numReadsToSkip; i++)
        {
                //read 4 lines to junk;
                getline(finFastqLeft, junk);
                getline(finFastqLeft, junk);
                getline(finFastqLeft, junk);
                getline(finFastqLeft, junk);
        }

        //Now the ifstream points to the beginning of the read that we'd like to extract
        int printReadIdStrInd = 0;
        char charRead = 'X';
        //Read one junk character whichshould be '@';
        finFastqLeft.get(charRead);

        assert(charRead == '@');
	
        while(!isReadNameCharInvalid[charRead]) //Table check for valid characters in a read name that is printed in SAM file
        {
                finFastqLeft.get(charRead);
                printReadIdStr[printReadIdStrInd++] = charRead;
        }
        printReadIdStr[printReadIdStrInd-1] = '\0';

	//Read the remaining characters on the read line to junk
	if(charRead != '\n')
	{
		getline(finFastqLeft, junk);
	}

	//Now read two junk lines and and the left_quality
        getline(finFastqLeft, junk);
        getline(finFastqLeft, junk);

	//Quality score is always a read length
        for(int i=0; i<finalReadLen; i++)
        {
        	finFastqLeft.get(leftQualityScore[i]);
	}
	leftQualityScore[finalReadLen] = '\0';
	char endLine;
	finFastqLeft.get(endLine);
	assert(endLine == '\n');
}

void PrintQualityScoreLeft(char *& out, bool dir)
{
        if(dir == 0)
        {
                for(int i=0; i<finalReadLen; i++)
                {
                        //out[0] = junk[i];
                        out[0] = leftQualityScore[i];
			out++;
                }
        }
        else
        {
                for(int i=0; i<finalReadLen; i++)
                {
                        //out[0] = junk[finalReadLen-i-1];
                        out[0] = leftQualityScore[finalReadLen-i-1];
                        out++;
                }

        }
}

void PrintQualityScoreRight(char *& out, bool dir, unsigned long long printReadCode)
{
	if(printReadCode != previousReadCodeRight)
	{
	        int numReadsToSkip = printReadCode - previousReadCodeRight - 1;

	        string junk;
	        for(int i=0; i<numReadsToSkip; i++)
	        {
	                //read 4 lines to junk;
	                getline(finFastqRight, junk);
	                getline(finFastqRight, junk);
	                getline(finFastqRight, junk);
	                getline(finFastqRight, junk);
	        }

		//Reached current read, so read three more lines
		getline(finFastqRight, junk);
		getline(finFastqRight, junk);
		getline(finFastqRight, junk);

		//This one is the quality score string we need
		for(int i=0; i<finalReadLen; i++)
		{
			finFastqRight.get(rightQualityScore[i]);
		}
		rightQualityScore[finalReadLen] = '\0';
		char endLine;
		finFastqRight.get(endLine);
		assert(endLine == '\n');
	}

        if(dir == 0)
        {
                for(int i=0; i<finalReadLen; i++)
                {
                        out[0] = rightQualityScore[i];
                        out++;
                }
        }
        else
        {
                for(int i=0; i<finalReadLen; i++)
                {
                        out[0] = rightQualityScore[finalReadLen-i-1];
                        out++;
                }
        }
}

void PrintSingleSamFromLinks(FILE* fout,/* const*/ mergedLink& singleItem, unsigned long long printReadCode, bool flushFlag)
{
        if(printReadNamesFlag)
        {
                SetCurrentReadName(printReadCode);
        }
        else
        {
                printReadIdStr[0] = 'R';
                ltoa10(printReadCode, printReadIdStr+1);
        }

	short singleFlag = (short) singleItem.dir * 16; //64 and 128 come from being the first and second fragments

	unsigned int singleChrPos = singleItem.chrPos;
	itoa10(singleChrPos, singleItemChrPosStr);
	short singleChrCode = singleItem.chrCode;

	//Handle and prep MAPQ (mapping quality)
	unsigned char singleMapQ = 255; //255 represents unspecified mapping quality
	if(mappingQualityPrintFlag)
	{
		singleMapQ = CalculateMapQ(singleItem.edit.numEdits);		
	}

	//Regardless of directiong the tempString stays same
	for(int i=0; i<finalReadLen; i++)
	{
		tempPrintStr[i] = fullRef[singleChrCode][singleChrPos + i];
	} 
	tempPrintStr[finalReadLen] = '\0';

	bool singleIndelFlag = 0; //tells whether there is any form of indel in the singleItem
	int singleInsertionOffset = 0;
	int singleDeletionOffset = 0;

	for(int k=0; k<singleItem.edit.numEdits; k++)
	{
		switch(singleItem.edit.edits[k])
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
				replacedQueryCharList[k] = tempPrintStr[singleItem.edit.pos[k] - 1 + singleInsertionOffset - singleDeletionOffset];
				tempPrintStr[singleItem.edit.pos[k] - 1 + singleInsertionOffset - singleDeletionOffset] = singleItem.edit.edits[k];
				break;	
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'n':
				{
					int insertionBlockLength = 1; //special case for insertion blocks longer than 1
					//Check if there are other insertions ahead to increase the insertion block length
					while(k < singleItem.edit.numEdits - insertionBlockLength && singleItem.edit.pos[k+insertionBlockLength] == singleItem.edit.pos[k] && singleItem.edit.edits[k + insertionBlockLength] >= 'a' && singleItem.edit.edits[k + insertionBlockLength] <= 't')
					{
						insertionBlockLength++;
					}

					for(int i=finalReadLen-1; i>=singleItem.edit.pos[k] + insertionBlockLength - 1 + singleInsertionOffset - singleDeletionOffset; i--) // >= is intentional since edits are 1 based string is 0 based;
					{
						tempPrintStr[i] = tempPrintStr[i-insertionBlockLength];
					}
					
					for(int i = k; i< k + insertionBlockLength; i++)
					{
						tempPrintStr[singleItem.edit.pos[i] + i - k - 1 + singleInsertionOffset - singleDeletionOffset] = singleItem.edit.edits[i] - 'a' + 'A'; 
					}

					k += insertionBlockLength - 1;
					singleIndelFlag = 1;
					
					singleInsertionOffset += insertionBlockLength;
					break;	
				}
			case 'D':
				replacedQueryCharList[k] = tempPrintStr[singleItem.edit.pos[k] - 1 + singleInsertionOffset - singleDeletionOffset];
			
				for(int i=singleItem.edit.pos[k] - 1 + singleInsertionOffset - singleDeletionOffset; i<finalReadLen-1; i++)
					tempPrintStr[i] = tempPrintStr[i+1];
			
				tempPrintStr[finalReadLen - 1] = fullRef[singleChrCode][singleChrPos + finalReadLen + singleDeletionOffset - singleInsertionOffset]; // It is actually singleChrPos + (finalReadLen - 1) + 1 but ones cancel out
			
				singleIndelFlag = 1;
				singleDeletionOffset++;

				break;

			default:
				cout << "Error UE3: unexpected edit character: '" << singleItem.edit.edits[k] << "'  (int = " << (int) singleItem.edit.edits[k] << ")" << endl;
				exit(98);
		}		
	}

	GetMDStrFromMeditWithCharReplacement(singleItem.edit, replacedQueryCharList, finalReadLen, singleMDstr); //direction is assumed to be forward
	
	char* str = printBuffer + bufferOffset;
	writeStr(str, printReadIdStr, '\t');
	writeStr(str, valStr[singleFlag], '\t'); //valstr contains string form of some integer values
	writeStr(str, chrNames[singleItem.chrCode], '\t');
	writeStr(str, singleItemChrPosStr, '\t'); 
	writeStr(str, valStr[singleMapQ], '\t');
	
	//CIGAR print
	if(!singleIndelFlag)
	{ 
		writeStr(str, valStr[finalReadLen]);
		str[0] = 'M'; str++;
	}
	else
	{
		int properReadLen=finalReadLen;

		//walk through the edits and create a cigar string
		int lastOffset = 0;
		for(int k=0; k<singleItem.edit.numEdits; k++)
		{
			switch(singleItem.edit.edits[k])
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					//ignore
					break;
				case 'a':
				case 'c':
				case 'g':
				case 't':
				case 'n':
					{
						if(singleItem.edit.pos[k] - lastOffset > 1)
							writeIntToStr(str, singleItem.edit.pos[k] - lastOffset - 1, 'M');
						
						int numInsertions = 1;
						while(k+1 < singleItem.edit.numEdits && singleItem.edit.pos[k+1] == singleItem.edit.pos[k] && singleItem.edit.edits[k+1] >= 'a' && singleItem.edit.edits[k+1] <= 't')
						{
							numInsertions++;
							k++;
						}
						properReadLen -= numInsertions;
						writeIntToStr(str, numInsertions, 'I'); 	
						lastOffset = singleItem.edit.pos[k] - 1;
						break;
					}					
				case 'D':
					{
						if(singleItem.edit.pos[k] - lastOffset > 1)
							writeIntToStr(str, singleItem.edit.pos[k] - lastOffset - 1, 'M');
		
						int numDeletions = 1;
						while(k+1 < singleItem.edit.numEdits && singleItem.edit.pos[k+1] == singleItem.edit.pos[k]+1 && singleItem.edit.edits[k+1] == 'D')
						{
							numDeletions++;
							k++;
						}
						properReadLen += numDeletions;
						writeIntToStr(str, numDeletions, 'D'); 	
						lastOffset = singleItem.edit.pos[k];
						break;					
					}
				default:
					cout << "Error UE4: unexpected edit character: '" << singleItem.edit.edits[k] << "'  (int = " << (int) singleItem.edit.edits[k] << ")" << endl;
					exit(97);
			}
		}		

		//Also put in the last block if not zero
		if(properReadLen - lastOffset > 0)
		{
			//write the difference
			writeIntToStr(str, properReadLen - lastOffset, 'M');	
		}
	}

	str[0] = '\t'; str++;
	str[0] = '*'; str++;
	str[0] = '\t'; str++;
	str[0] = '0'; str++;
	str[0] = '\t'; str++;
	str[0] = '0'; str++;
	str[0] = '\t'; str++;
	writeStr(str, tempPrintStr, '\t');

        if(printQualityScoresFlag)
        {
		PrintQualityScoreLeft(str, singleItem.dir);
	}
        else
        {
                str[0] = '*'; str++;
        }
	str[0] = '\t'; str++;

	str[0] = 'N'; str++;
	str[0] = 'M'; str++;
	str[0] = ':'; str++;
	str[0] = 'i'; str++;
	str[0] = ':'; str++;
	writeStr(str, valStr[singleItem.edit.numEdits], '\t');
	str[0] = 'M'; str++;
	str[0] = 'D'; str++;
	str[0] = ':'; str++;
	str[0] = 'Z'; str++;
	str[0] = ':'; str++;
	writeStr(str, singleMDstr, '\n');

	bufferOffset = str - printBuffer; //Get how far away current string-write pointer is from the beginning of the buffer

	previousReadCodeLeft = printReadCode;

	if(flushFlag && bufferOffset > printBufferSize)
        {
              flushBuffer(fout);
        }
}

void PrintPairSamFromLinks(FILE* fout, const mergedLink& leftItem, const mergedLink& rightItem, unsigned long long printReadCode, bool flushFlag) //readcode is the joint name for all splits and mates
{
        if(printReadNamesFlag)
        {
                SetCurrentReadName(printReadCode);
        }
        else
        {
                printReadIdStr[0] = 'R';
                ltoa10(printReadCode, printReadIdStr+1);
        }

	short leftFlag = (short) leftItem.dir * 16 + (short) rightItem.dir * 32 + 64 + 1; //64 and 128 come from being the first and second fragments
	short rightFlag = (short) rightItem.dir * 16 + (short) leftItem.dir * 32 + 128 + 1;

	unsigned int leftChrPos = leftItem.chrPos;
	itoa10(leftChrPos, leftItemChrPosStr);	
	short leftChrCode = leftItem.chrCode;
	unsigned int rightChrPos = rightItem.chrPos;
	itoa10(rightChrPos, rightItemChrPosStr);	
	short rightChrCode = rightItem.chrCode;

	//Handle and prep MAPQ (mapping quality)
	unsigned char leftMapQ = 255; //255 represents unspecified mapping quality
	if(mappingQualityPrintFlag)
	{
		leftMapQ = CalculateMapQ(leftItem.edit.numEdits + rightItem.edit.numEdits);		
	}
	unsigned char rightMapQ = leftMapQ; //Both mates have the same mapping quality

	//Here we handle the left mate
	int tLen = abs((int) rightChrPos - (int)leftChrPos) + finalReadLen;
	if(rightChrPos < leftChrPos)
	{
		tLen = -tLen;
	}	
	itoa10(tLen, tLenStr_plus);
	itoa10(-tLen, tLenStr_minus);

	bool isProperlyAligned = 1; // Do the reads have alignments looking toward each other on different strands (No eversion, no inversion of one end)

	if(leftItem.dir == rightItem.dir)
	{
		isProperlyAligned = 0;
	}
	else
	{
		if(leftItem.dir == 0)
		{
			if(rightChrPos < leftChrPos)
			{
				isProperlyAligned = 0;
			}
		}
		else
		{
			if(rightChrPos > leftChrPos)
			{
				isProperlyAligned = 0;
			}
		}
	}

	leftFlag += 2 * (short) isProperlyAligned;
	rightFlag += 2 * (short) isProperlyAligned;	

	//Regardless of directiong the tempString stays same
	for(int i=0; i<finalReadLen; i++)
	{
		tempPrintStr[i] = fullRef[leftChrCode][leftChrPos + i];
	} 
	tempPrintStr[finalReadLen] = '\0';

	bool leftIndelFlag = 0; //tells whether there is any form of indel in the leftItem

	int leftInsertionOffset = 0;
	int leftDeletionOffset = 0;

	for(int k=0; k<leftItem.edit.numEdits; k++)
	{
		switch(leftItem.edit.edits[k])
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
				replacedQueryCharList[k] = tempPrintStr[leftItem.edit.pos[k] - 1 + leftInsertionOffset - leftDeletionOffset];
				tempPrintStr[leftItem.edit.pos[k] - 1 + leftInsertionOffset - leftDeletionOffset] = leftItem.edit.edits[k];
				break;	
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'n':
				{
					int insertionBlockLength = 1; //special case for insertion blocks longer than 1
					//Check if there are other insertions ahead to increase the insertion block length
					while(k < leftItem.edit.numEdits - insertionBlockLength && leftItem.edit.pos[k+insertionBlockLength] == leftItem.edit.pos[k] && leftItem.edit.edits[k + insertionBlockLength] >= 'a' && leftItem.edit.edits[k + insertionBlockLength] <= 't')
					{
						insertionBlockLength++;
					}

					for(int i=finalReadLen-1; i>=leftItem.edit.pos[k] + insertionBlockLength - 1 + leftInsertionOffset - leftDeletionOffset; i--) // >= is intentional since edits are 1 based string is 0 based;
					{
						tempPrintStr[i] = tempPrintStr[i-insertionBlockLength];
					}
					
					for(int i = k; i< k + insertionBlockLength; i++)
					{
						tempPrintStr[leftItem.edit.pos[i] + i - k - 1 + leftInsertionOffset - leftDeletionOffset] = leftItem.edit.edits[i] - 'a' + 'A'; 
					}

					k += insertionBlockLength - 1;
					leftIndelFlag = 1;
					
					leftInsertionOffset += insertionBlockLength;
					break;	
				}
			case 'D':
				replacedQueryCharList[k] = tempPrintStr[leftItem.edit.pos[k] - 1 + leftInsertionOffset - leftDeletionOffset];
			
				for(int i=leftItem.edit.pos[k] - 1 + leftInsertionOffset - leftDeletionOffset; i<finalReadLen-1; i++)
					tempPrintStr[i] = tempPrintStr[i+1];
			
				tempPrintStr[finalReadLen - 1] = fullRef[leftChrCode][leftChrPos + finalReadLen + leftDeletionOffset - leftInsertionOffset]; // It is actually leftChrPos + (finalReadLen - 1) + 1 but ones cancel out
			
				leftIndelFlag = 1;
				leftDeletionOffset++;

				break;

			default:
				cout << "Error UE3: unexpected edit character: '" << leftItem.edit.edits[k] << "'  (int = " << (int) leftItem.edit.edits[k] << ")" << endl;
				exit(98);
		}		
	}

	GetMDStrFromMeditWithCharReplacement(leftItem.edit, replacedQueryCharList, finalReadLen, leftMDstr); //direction is assumed to be forward

	char* str = printBuffer + bufferOffset;
	writeStr(str, printReadIdStr, '\t');
	writeStr(str, valStr[leftFlag], '\t'); //valstr contains string form of some integer values
	writeStr(str, chrNames[leftItem.chrCode], '\t');
	writeStr(str, leftItemChrPosStr, '\t'); 
	writeStr(str, valStr[leftMapQ], '\t');

	if(!leftIndelFlag)
	{ 
		writeStr(str, valStr[finalReadLen]);
		str[0] = 'M'; str++;
	}
	else
	{
		int properReadLen=finalReadLen;

		//walk through the edits and create a cigar string
		int lastOffset = 0;
		for(int k=0; k<leftItem.edit.numEdits; k++)
		{
			switch(leftItem.edit.edits[k])
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					//ignore
					break;
				case 'a':
				case 'c':
				case 'g':
				case 't':
				case 'n':
					{
						if(leftItem.edit.pos[k] - lastOffset > 1)
							writeIntToStr(str, leftItem.edit.pos[k] - lastOffset - 1, 'M');
						
						int numInsertions = 1;
						while(k+1 < leftItem.edit.numEdits && leftItem.edit.pos[k+1] == leftItem.edit.pos[k] && leftItem.edit.edits[k+1] >= 'a' && leftItem.edit.edits[k+1] <= 't')
						{
							numInsertions++;
							k++;
						}
						properReadLen -= numInsertions;
						writeIntToStr(str, numInsertions, 'I'); 	
						lastOffset = leftItem.edit.pos[k] - 1;
						break;
					}					
				case 'D':
					{
						if(leftItem.edit.pos[k] - lastOffset > 1)
							writeIntToStr(str, leftItem.edit.pos[k] - lastOffset - 1, 'M');
		
						int numDeletions = 1;
						while(k+1 < leftItem.edit.numEdits && leftItem.edit.pos[k+1] == leftItem.edit.pos[k]+1 && leftItem.edit.edits[k+1] == 'D')
						{
							numDeletions++;
							k++;
						}
						properReadLen += numDeletions;
						writeIntToStr(str, numDeletions, 'D'); 	
						lastOffset = leftItem.edit.pos[k];
						break;					
					}
				default:
					cout << "Error UE4: unexpected edit character: '" << leftItem.edit.edits[k] << "'  (int = " << (int) leftItem.edit.edits[k] << ")" << endl;
					exit(97);
			}
		}		

		//Also put in the last block if not zero
		if(properReadLen - lastOffset > 0)
		{
			//write the difference
			writeIntToStr(str, properReadLen - lastOffset, 'M');	
		}
	}


	str[0] = '\t'; str++;
	str[0] = '='; str++;
	str[0] = '\t'; str++;
	writeStr(str, rightItemChrPosStr, '\t');
	writeStr(str, tLenStr_plus, '\t');
	writeStr(str, tempPrintStr, '\t');

        if(printQualityScoresFlag)
        {
		PrintQualityScoreLeft(str, leftItem.dir);
	}
        else
        {
                str[0] = '*'; str++;
        }
	str[0] = '\t'; str++;

	str[0] = 'N'; str++;
	str[0] = 'M'; str++;
	str[0] = ':'; str++;
	str[0] = 'i'; str++;
	str[0] = ':'; str++;
	writeStr(str, valStr[leftItem.edit.numEdits], '\t');
	str[0] = 'M'; str++;
	str[0] = 'D'; str++;
	str[0] = ':'; str++;
	str[0] = 'Z'; str++;
	str[0] = ':'; str++;
	writeStr(str, leftMDstr, '\n');

	// Everything until here was about the left mate
	//////////////////////////////////////////////


	//Here we handle the right mate
	for(int i=0; i<finalReadLen; i++)
	{
		tempPrintStr[i] = fullRef[rightChrCode][rightChrPos + i];
	} 
	tempPrintStr[finalReadLen] = '\0';

	bool rightIndelFlag = 0;

	int rightInsertionOffset = 0;
	int rightDeletionOffset = 0;

	for(int k=0; k<rightItem.edit.numEdits; k++)
	{
		switch(rightItem.edit.edits[k])
		{
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
				replacedQueryCharList[k] = tempPrintStr[rightItem.edit.pos[k] - 1 + rightInsertionOffset - rightDeletionOffset];
				tempPrintStr[rightItem.edit.pos[k] - 1 + rightInsertionOffset - rightDeletionOffset] = rightItem.edit.edits[k];
				break;
			case 'a':
			case 'c':
			case 'g':
			case 't':
			case 'n':
				{
					int insertionBlockLength = 1; //special case for insertion blocks longer than 1
					while(k + insertionBlockLength < rightItem.edit.numEdits && rightItem.edit.pos[k+insertionBlockLength] == rightItem.edit.pos[k] && rightItem.edit.edits[k + insertionBlockLength] >= 'a' && rightItem.edit.edits[k + insertionBlockLength] <= 't')
					{
						insertionBlockLength++;
					}

					for(int i=finalReadLen-1; i>=rightItem.edit.pos[k] + insertionBlockLength - 1 + rightInsertionOffset - rightDeletionOffset; i--) // >= is intentional since edits are 1 based string is 0 based;
					{
						tempPrintStr[i] = tempPrintStr[i-insertionBlockLength];
					}
		
					for(int i = k; i< k + insertionBlockLength; i++)
					{
						tempPrintStr[rightItem.edit.pos[i] + i - k - 1 + rightInsertionOffset - rightDeletionOffset] = rightItem.edit.edits[i] - 'a' + 'A'; 
					}

					k += insertionBlockLength - 1;
					rightIndelFlag = 1;
					
					rightInsertionOffset += insertionBlockLength;
					
					break;	
				}
			case 'D':
				replacedQueryCharList[k] = tempPrintStr[rightItem.edit.pos[k] - 1 + rightInsertionOffset - rightDeletionOffset];
				for(int i=rightItem.edit.pos[k] - 1 + rightInsertionOffset - rightDeletionOffset; i<finalReadLen-1; i++)
					tempPrintStr[i] = tempPrintStr[i+1];
				tempPrintStr[finalReadLen - 1] = fullRef[rightChrCode][rightChrPos + finalReadLen + rightDeletionOffset - rightInsertionOffset]; // It is actually righthrPos + (finalReadLen - 1) + 1 but ones cancel out
				rightIndelFlag = 1;
				rightDeletionOffset++;
				break;
			default:
				cout << "Error UE5: unexpected edit character: '" << rightItem.edit.edits[k] << "'  (int = " << (int) rightItem.edit.edits[k] << ")" << endl;
				exit(99);
		}
	}

	GetMDStrFromMeditWithCharReplacement(rightItem.edit, replacedQueryCharList, finalReadLen, rightMDstr); //direction is assumed to be forward

	writeStr(str, printReadIdStr, '\t');
	writeStr(str, valStr[rightFlag], '\t'); //valstr contains string form of some integer values
	writeStr(str, chrNames[rightItem.chrCode], '\t');
	writeStr(str, rightItemChrPosStr, '\t'); 
	writeStr(str, valStr[rightMapQ], '\t');
	
	if(!rightIndelFlag)
	{ 
		writeStr(str, valStr[finalReadLen]);
		str[0] = 'M'; str++;
	}
	else
	{
		int properReadLen=finalReadLen;

		//walk through the edits and create a cigar string
		int lastOffset = 0;
		for(int k=0; k<rightItem.edit.numEdits; k++)
		{
			switch(rightItem.edit.edits[k])
			{
				case 'A':
				case 'C':
				case 'G':
				case 'T':
				case 'N':
					//ignore
					break;
				case 'a':
				case 'c':
				case 'g':
				case 't':
				case 'n':
					{
						if(rightItem.edit.pos[k] - lastOffset > 1)
							writeIntToStr(str, rightItem.edit.pos[k] - lastOffset - 1, 'M');
						
						int numInsertions = 1;
						while(k+1 < rightItem.edit.numEdits && rightItem.edit.pos[k+1] == rightItem.edit.pos[k] && rightItem.edit.edits[k+1] >= 'a' && rightItem.edit.edits[k+1] <= 't')
						{
							numInsertions++;
							k++;
						}
						properReadLen -= numInsertions;
						writeIntToStr(str, numInsertions, 'I'); 	
						lastOffset = rightItem.edit.pos[k] + -1;
						break;
					}					
				case 'D':
					{
						if(rightItem.edit.pos[k] - lastOffset > 1)
							writeIntToStr(str, rightItem.edit.pos[k] - lastOffset - 1, 'M');
		
						int numDeletions = 1;
						while(k+1 < rightItem.edit.numEdits && rightItem.edit.pos[k+1] == rightItem.edit.pos[k]+1 && rightItem.edit.edits[k+1] == 'D')
						{
							numDeletions++;
							k++;
						}
						properReadLen += numDeletions;

						writeIntToStr(str, numDeletions, 'D'); 	
						lastOffset = rightItem.edit.pos[k];	
						break;					
					}
				default:
					cout << "Error UE6: unexpected edit character: '" << rightItem.edit.edits[k] << "'  (int = " << (int) rightItem.edit.edits[k] << ")" << endl;
					exit(97);
			}
		}		

		//Also put in the last block if not zero
		if(properReadLen - lastOffset > 0)
		{
			//write the difference
			writeIntToStr(str, properReadLen - lastOffset, 'M');	
		}
	}

	str[0] = '\t'; str++;
	str[0] = '='; str++;
	str[0] = '\t'; str++;
	writeStr(str, leftItemChrPosStr, '\t');
	writeStr(str, tLenStr_minus, '\t');
	writeStr(str, tempPrintStr, '\t');

        if(printQualityScoresFlag)
        {
                PrintQualityScoreRight(str, rightItem.dir, printReadCode);
                previousReadCodeRight = printReadCode;
        }
        else
        {
                str[0] = '*'; str++;
        }

	str[0] = '\t'; str++;
	str[0] = 'N'; str++;
	str[0] = 'M'; str++;
	str[0] = ':'; str++;
	str[0] = 'i'; str++;
	str[0] = ':'; str++;
	writeStr(str, valStr[rightItem.edit.numEdits], '\t');
	str[0] = 'M'; str++;
	str[0] = 'D'; str++;
	str[0] = ':'; str++;
	str[0] = 'Z'; str++;
	str[0] = ':'; str++;
	writeStr(str, rightMDstr, '\n');

	bufferOffset = str - printBuffer; //Get how far away current string-write pointer is from the beginning of the buffer

	previousReadCodeLeft = printReadCode;
	
	//flushing is done by the function that calls printing
	if(flushFlag && bufferOffset > printBufferSize)
        {
              flushBuffer(fout);
        }
}

doubleMedit tempMedit;

//Merges two medits, provided that they are ordered and no edit goes in between other's
void SimpleOrderedMergeMedits(doubleMedit& toMe, const doubleMedit& fromMe)
{
	unsigned char toNumEdits = toMe.numEdits;
	if(toNumEdits > 0 && toMe.pos[0] > fromMe.pos[0])
	{
		//fromMe will be on the left side
		unsigned char fromNumEdits = fromMe.numEdits;
		for(int k=toNumEdits + fromNumEdits - 1; k>=fromNumEdits; k--)
		{
			toMe.pos[k] = toMe.pos[k - fromNumEdits];
			toMe.edits[k] = toMe.edits[k - fromNumEdits];
		}
		
		for(int k=0; k<fromNumEdits; k++)
		{
			toMe.pos[k] = fromMe.pos[k];	
			toMe.edits[k] = fromMe.edits[k];
		}
	}
	else
	{
		//fromMe will be on the right side
		for(int k=0; k<fromMe.numEdits; k++)
		{
			toMe.pos[toNumEdits + k] = fromMe.pos[k];
			toMe.edits[toNumEdits + k] = fromMe.edits[k]; 
		}
	}
	toMe.numEdits = toNumEdits + fromMe.numEdits;
}


//Merges two medits, provided that they are ordered and no edit goes in between other's
char SimpleOrderedMergeMedits(doubleMedit& toMe, const doubleMedit& fromMe, char posOffset, char posIndelModifier) 
//Overloaded indel version of the same function above -- properly shifts chrPos (returned value) and edit positions depending on the indel offset and which side the split with the indel goes to
{
	unsigned char toNumEdits = toMe.numEdits;
	
	if(posOffset == 0) //fromMe will be inserted to the left, existing positions in toMe will be offset by indelModifier together with the newly inserted fromMe, returned value will be the indelModifier offset for chrPos
	{
		//fromMe will be on the left side
		unsigned char fromNumEdits = fromMe.numEdits;
		for(int k=toNumEdits + fromNumEdits - 1; k>=fromNumEdits; k--)
		{
			toMe.pos[k] = toMe.pos[k - fromNumEdits] - posIndelModifier;
			toMe.edits[k] = toMe.edits[k - fromNumEdits];
		}
		
		for(int k=0; k<fromNumEdits; k++)
		{
			toMe.pos[k] = fromMe.pos[k] - posIndelModifier;	
			toMe.edits[k] = fromMe.edits[k];
		}
	
		toMe.numEdits = toNumEdits + fromMe.numEdits;
		
		//fromMe will be on the left side, so modify chrPos accordingly
		return posIndelModifier;
	}
	else //fromMe will be inserted to the right, existing positions in toMe will not be modified, chrPos will not be modified
	{
		//fromMe will be on the right side
		for(int k=0; k<fromMe.numEdits; k++)
		{
			toMe.pos[toNumEdits + k] = fromMe.pos[k] + posOffset - posIndelModifier;
			toMe.edits[toNumEdits + k] = fromMe.edits[k]; 
		}
		toMe.numEdits = toNumEdits + fromMe.numEdits;
		
		return 0;
	}
}


//Checks if the corresponding noHit read-mer sequence still stays within the allowed error rate for the single-end read
bool IsNoHitExtensionViable_SingleEnd(mergedLink& singleLink, unsigned long long printReadCode)
{
	unsigned long long splitBaseCode = (printReadCode -1)*2; //Times 2 for single-end reads

	tempMedit.numEdits = 0;

	unsigned char extensionErrorLimit;
	bool missDir;
	char* curRef;

	int singleSplitType = singleLink.edit.numEdits / MISSING_SPLIT_OFFSET;
	unsigned char posOffset = 0;

	char* noHitSeq;
	bool noHitDir;

	bool indelBoundaryDir = 0; //This is the boundary for indel alignment
        if(singleSplitType == 1) //first split missing
        {
                //Obtain nohit arr
		noHitSeq = noHitArr + links[splitBaseCode].chrPos;
		noHitDir = links[splitBaseCode].dir;
        	singleLink.edit.numEdits %= MISSING_SPLIT_OFFSET;
		extensionErrorLimit = doubleNumMismatchesPerReadMer - singleLink.edit.numEdits;
		missDir = singleLink.dir;
		curRef = fullRef[singleLink.chrCode] + singleLink.chrPos;
		if(missDir == 1)
		{
			posOffset = travArrReadmerLength[1];
			curRef += travArrReadmerLength[1]; //1 is the indel modified length of the second split
		}
		else
		{
			indelBoundaryDir = 1;
		}
	}
        else// if(singleSplitType == 2)
        {
		noHitSeq = noHitArr + links[splitBaseCode+1].chrPos;
		noHitDir = links[splitBaseCode+1].dir;
		singleLink.edit.numEdits %= MISSING_SPLIT_OFFSET;
        	extensionErrorLimit = doubleNumMismatchesPerReadMer - singleLink.edit.numEdits;
		missDir = singleLink.dir;
		curRef = fullRef[singleLink.chrCode] + singleLink.chrPos;
		if(missDir == 0)
		{
			posOffset = travArrReadmerLength[0];
			curRef += travArrReadmerLength[0];
		}
		else
		{
			indelBoundaryDir = 1;
		}
	}

	bool tryIndelFlag = 0;
	if(noHitDir == missDir)
	{
		for(int k=0; k<readLen; k++)
                {
                        if(curRef[k] != noHitSeq[k])
                        {
                                if(tempMedit.numEdits >= extensionErrorLimit)
                                {
					tempMedit.numEdits = 0;	
					tryIndelFlag=1;
                                        break;
                                }
                                tempMedit.pos[tempMedit.numEdits] = k + 1 + posOffset;
                                tempMedit.edits[tempMedit.numEdits] = noHitSeq[k];
                                tempMedit.numEdits++;
                        }
                }
	}
	else
	{
		for(int k=0; k<readLen; k++)
                {
                        if(curRef[k] != RevCompChar[(unsigned char) noHitSeq[readLen-k-1]])
                        {
                                if(tempMedit.numEdits >= extensionErrorLimit)
                                {
					tempMedit.numEdits = 0;
					tryIndelFlag=1;
					break;
                                }
                                tempMedit.pos[tempMedit.numEdits] = k + 1 + posOffset; //medit positions stay in the forward direction
                            	tempMedit.edits[tempMedit.numEdits] = RevCompChar[(unsigned char) noHitSeq[readLen-k-1]];
                                tempMedit.numEdits++;
                        }
                }
	}

	if(tryIndelFlag)
	{
		if(distanceMetric == HAMMING)
		{
			return 0;
		}

		char posOffsetIndelModifier = 0; //This is shifting of the offset towards left or right based on the right bounded indel alignment 

		if(!IndelAlignment(tempMedit, noHitSeq, curRef, !(noHitDir == missDir), extensionErrorLimit, indelBoundaryDir, posOffsetIndelModifier))
		{
			return 0;
		}

		singleLink.chrPos += SimpleOrderedMergeMedits(singleLink.edit, tempMedit, posOffset, posOffsetIndelModifier); //They are ordered, but which order is not known, they should be merged into the first one
	}
	else
	{
		SimpleOrderedMergeMedits(singleLink.edit, tempMedit); //They are ordered, but which order is not known, they should be merged into the first one
	}
	return 1;
}


//Checks if the corresponding noHit read-mer sequence still stays within the allowed error rate for the whole read
bool IsNoHitExtensionViable(mergedLink& leftLink, mergedLink& rightLink, unsigned long long printReadCode)
{
	unsigned long long splitBaseCode = (printReadCode -1)*4; //

	tempMedit.numEdits = 0;

	unsigned char extensionErrorLimit;
	bool missDir;
	short missMate;
	char* curRef;

	int leftSplitType = leftLink.edit.numEdits / MISSING_SPLIT_OFFSET;
	unsigned char posOffset = 0;

	char* noHitSeq;
	bool noHitDir;

	bool indelBoundaryDir = 0; //This is the boundary for indel alignment
        if(leftSplitType == 1) //first split missing
        {
                //Obtain nohit arr
		missMate = 1;
		noHitSeq = noHitArr + links[splitBaseCode].chrPos;
		noHitDir = links[splitBaseCode].dir;
        	leftLink.edit.numEdits %= MISSING_SPLIT_OFFSET;
		extensionErrorLimit = doubleNumMismatchesPerReadMer - leftLink.edit.numEdits;
		missDir = leftLink.dir;
		curRef = fullRef[leftLink.chrCode] + leftLink.chrPos;
		if(missDir == 1)
		{
			posOffset = travArrReadmerLength[1];
			curRef += travArrReadmerLength[1]; //1 is the indel modified length of the second split
		}
		else
		{
			indelBoundaryDir = 1;
		}
	}
        else if(leftSplitType == 2)
        {
		missMate = 1;
		noHitSeq = noHitArr + links[splitBaseCode+1].chrPos;
		noHitDir = links[splitBaseCode+1].dir;
		leftLink.edit.numEdits %= MISSING_SPLIT_OFFSET;
        	extensionErrorLimit = doubleNumMismatchesPerReadMer - leftLink.edit.numEdits;
		missDir = leftLink.dir;
		curRef = fullRef[leftLink.chrCode] + leftLink.chrPos;
		if(missDir == 0)
		{
			posOffset = travArrReadmerLength[0];
			curRef += travArrReadmerLength[0];
		}
		else
		{
			indelBoundaryDir = 1;
		}
	}
        else
        {
		missMate = 2;
                int rightSplitType = rightLink.edit.numEdits / MISSING_SPLIT_OFFSET;
                if(rightSplitType == 1)
                {
			noHitSeq = noHitArr + links[splitBaseCode+2].chrPos;
			noHitDir = links[splitBaseCode+2].dir;
			rightLink.edit.numEdits %= MISSING_SPLIT_OFFSET;
        		extensionErrorLimit = doubleNumMismatchesPerReadMer - rightLink.edit.numEdits;
			missDir = rightLink.dir;
			curRef = fullRef[rightLink.chrCode] + rightLink.chrPos;
			if(missDir == 1)
			{
				posOffset = travArrReadmerLength[3];
				curRef += travArrReadmerLength[3];
			}
			else
			{
				indelBoundaryDir = 1;
			}
                }
                else
                {
			noHitSeq = noHitArr + links[splitBaseCode+3].chrPos;
			noHitDir = links[splitBaseCode+3].dir;
        		rightLink.edit.numEdits %= MISSING_SPLIT_OFFSET;
			extensionErrorLimit = doubleNumMismatchesPerReadMer - rightLink.edit.numEdits;
			missDir = rightLink.dir;
			curRef = fullRef[rightLink.chrCode] + rightLink.chrPos;
			if(missDir == 0)
			{
				posOffset = travArrReadmerLength[2];
				curRef += travArrReadmerLength[2];
			}
			else
			{
				indelBoundaryDir = 1;
			}
	        }
        }

	bool tryIndelFlag = 0;
	if(noHitDir == missDir)
	{
		for(int k=0; k<readLen; k++)
                {
                        if(curRef[k] != noHitSeq[k])
                        {
                                if(tempMedit.numEdits >= extensionErrorLimit)
                                {
					tempMedit.numEdits = 0;	
					tryIndelFlag=1;
                                        break;
                                }
                                tempMedit.pos[tempMedit.numEdits] = k + 1 + posOffset;
                                tempMedit.edits[tempMedit.numEdits] = noHitSeq[k];
                                tempMedit.numEdits++;
                        }
                }
	}
	else
	{
		for(int k=0; k<readLen; k++)
                {
                        if(curRef[k] != RevCompChar[(unsigned char) noHitSeq[readLen-k-1]])
                        {
                                if(tempMedit.numEdits >= extensionErrorLimit)
                                {
					tempMedit.numEdits = 0;
					tryIndelFlag=1;
					break;
                                }
                                tempMedit.pos[tempMedit.numEdits] = k + 1 + posOffset; //medit positions stay in the forward direction
                            	tempMedit.edits[tempMedit.numEdits] = RevCompChar[(unsigned char) noHitSeq[readLen-k-1]];
                                tempMedit.numEdits++;
                        }
                }
	}

	if(tryIndelFlag)
	{
		if(distanceMetric == HAMMING)
		{
			return 0;
		}

		char posOffsetIndelModifier = 0; //This is shifting of the offset towards left or right based on the right bounded indel alignment 

		if(!IndelAlignment(tempMedit, noHitSeq, curRef, !(noHitDir == missDir), extensionErrorLimit, indelBoundaryDir, posOffsetIndelModifier))
		{
			return 0;
		}

		if(missMate == 1)
		{
			leftLink.chrPos += SimpleOrderedMergeMedits(leftLink.edit, tempMedit, posOffset, posOffsetIndelModifier); //They are ordered, but which order is not known, they should be merged into the first one
		}
		else
		{
			rightLink.chrPos += SimpleOrderedMergeMedits(rightLink.edit, tempMedit, posOffset, posOffsetIndelModifier);
		}
	}
	else
	{
		if(missMate == 1)
		{
			SimpleOrderedMergeMedits(leftLink.edit, tempMedit); //They are ordered, but which order is not known, they should be merged into the first one
		}
		else
		{
			SimpleOrderedMergeMedits(rightLink.edit, tempMedit);
		}
	}
	return 1;
}

//Identifies mate-pairs and prints all that is within specified insert size
void MatchMatesAndPrint(unsigned char arrNo1, unsigned char arrNo2, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* foutPair, unsigned long long printReadCode)
{
	unsigned int leftIndex = 0, rightBegIndex = 0, leftMateListSize = travArr[arrNo1].size(), rightMateListSize = travArr[arrNo2].size();
	while(leftIndex < leftMateListSize)
	{
		while(rightBegIndex < rightMateListSize && (travArr[arrNo2][rightBegIndex].chrCode < travArr[arrNo1][leftIndex].chrCode || (travArr[arrNo2][rightBegIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][rightBegIndex].chrPos < travArr[arrNo1][leftIndex].chrPos - (maxFragmentSize - finalReadLen) )))
		{
			rightBegIndex++; //increment to position that satisfies left-side upper distance threshold
		}			
		unsigned int curRightIndex = rightBegIndex; //copy index since right beg index will be the starting point of the next left guy

		while(curRightIndex < rightMateListSize &&  travArr[arrNo2][curRightIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][curRightIndex].chrPos <= travArr[arrNo1][leftIndex].chrPos + (maxFragmentSize - finalReadLen))
		{
			//So far is the maximum interval that satisfies left-side or right-side pairs
			if(abs((int) travArr[arrNo2][curRightIndex].chrPos - (int) travArr[arrNo1][leftIndex].chrPos) >= (minFragmentSize - finalReadLen))
			{ //this filters out pairs that are too close to each other						
				if(travArr[arrNo1][leftIndex].edit.numEdits + travArr[arrNo2][curRightIndex].edit.numEdits >=  MISSING_SPLIT_OFFSET)
				{
					if(travArr[arrNo1][leftIndex].edit.numEdits < MISSING_SPLIT_OFFSET || travArr[arrNo2][curRightIndex].edit.numEdits < MISSING_SPLIT_OFFSET) //if both are missing don't waste time
					{
						mergedLink leftLink = travArr[arrNo1][leftIndex];
						mergedLink rightLink = travArr[arrNo2][curRightIndex];
						if(IsNoHitExtensionViable(leftLink, rightLink, printReadCode))
						{
							PrintPairSamFromLinks(foutPair, leftLink, rightLink, printReadCode, 1);
						}
					}
				}
				else
				{
					//This is the normal case where the pairs are in proper distance and there is not missing read-mer
					PrintPairSamFromLinks(foutPair, travArr[arrNo1][leftIndex], travArr[arrNo2][curRightIndex], printReadCode, 1);
				}
			}
			curRightIndex++;
		}
		leftIndex++;
	}
	
	if(bufferOffset > printBufferSize)
	{
		flushBuffer(foutPair);
	}
}

//Identifies mate-pairs and prints all that is within specified insert size
void MatchMatesAndPrint_MapCountLimit(unsigned char arrNo1, unsigned char arrNo2, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* foutPair, unsigned long long printReadCode, int curMapCountLimit)
{
	unsigned int leftIndex = 0, rightBegIndex = 0, leftMateListSize = travArr[arrNo1].size(), rightMateListSize = travArr[arrNo2].size();
	while(leftIndex < leftMateListSize)
	{
		while(rightBegIndex < rightMateListSize && (travArr[arrNo2][rightBegIndex].chrCode < travArr[arrNo1][leftIndex].chrCode || (travArr[arrNo2][rightBegIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][rightBegIndex].chrPos < travArr[arrNo1][leftIndex].chrPos - (maxFragmentSize - finalReadLen) )))
		{
			rightBegIndex++; //increment to position that satisfies left-side upper distance threshold
		}			
		unsigned int curRightIndex = rightBegIndex; //copy index since right beg index will be the starting point of the next left guy

		while(curRightIndex < rightMateListSize &&  travArr[arrNo2][curRightIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][curRightIndex].chrPos <= travArr[arrNo1][leftIndex].chrPos + (maxFragmentSize - finalReadLen))
		{
			//So far is the maximum interval that satisfies left-side or right-side pairs
			if(abs((int) travArr[arrNo2][curRightIndex].chrPos - (int) travArr[arrNo1][leftIndex].chrPos) >= (minFragmentSize - finalReadLen))
			{ //this filters out pairs that are too close to each other						
				if(travArr[arrNo1][leftIndex].edit.numEdits + travArr[arrNo2][curRightIndex].edit.numEdits >=  MISSING_SPLIT_OFFSET)
				{
					if(travArr[arrNo1][leftIndex].edit.numEdits < MISSING_SPLIT_OFFSET || travArr[arrNo2][curRightIndex].edit.numEdits < MISSING_SPLIT_OFFSET) //if both are missing don't waste time
					{
						mergedLink leftLink = travArr[arrNo1][leftIndex];
						mergedLink rightLink = travArr[arrNo2][curRightIndex];
						if(IsNoHitExtensionViable(leftLink, rightLink, printReadCode))
						{
							PrintPairSamFromLinks(foutPair, leftLink, rightLink, printReadCode, 1);
							curMapCountLimit--;
							if(!curMapCountLimit)
							{
								if(bufferOffset > printBufferSize)
								{
									flushBuffer(foutPair);
								}
								return;
							}
						}
					}
				}
				else
				{
					//This is the normal case where the pairs are in proper distance and there is not missing read-mer
					PrintPairSamFromLinks(foutPair, travArr[arrNo1][leftIndex], travArr[arrNo2][curRightIndex], printReadCode, 1);
					curMapCountLimit--;
					if(!curMapCountLimit)
					{
						if(bufferOffset > printBufferSize)
						{
							flushBuffer(foutPair);
						}
						return;
					}


				}
			}
			curRightIndex++;
		}
		leftIndex++;
	}
	
	if(bufferOffset > printBufferSize)
	{
		flushBuffer(foutPair);
	}
}


//Stratum mapping version of above (print all mappings within the best tire of errors)
void MatchMatesAndPrint_StratumMapping(unsigned char arrNo1, unsigned char arrNo2, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* foutPair, unsigned long long printReadCode)
{
	short bestErrorSoFar = 999; //this acts as SHORT_MAX
	int startingbufferOffset = bufferOffset; //we save this since after we start filling in the buffer, if we hit a better mapping, we can just restore the original offest in the buffer and start overwriting


	unsigned int leftIndex = 0, rightBegIndex = 0, leftMateListSize = travArr[arrNo1].size(), rightMateListSize = travArr[arrNo2].size();
	while(leftIndex < leftMateListSize)
	{

		while(rightBegIndex < rightMateListSize && (travArr[arrNo2][rightBegIndex].chrCode < travArr[arrNo1][leftIndex].chrCode || (travArr[arrNo2][rightBegIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][rightBegIndex].chrPos < travArr[arrNo1][leftIndex].chrPos - (maxFragmentSize - finalReadLen) )))
		{
			rightBegIndex++; //increment to position that satisfies left-side upper distance threshold
		}			
		unsigned int curRightIndex = rightBegIndex; //copy index since right beg index will be the starting point of the next left guy
		while(curRightIndex < rightMateListSize &&  travArr[arrNo2][curRightIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][curRightIndex].chrPos <= travArr[arrNo1][leftIndex].chrPos + (maxFragmentSize - finalReadLen))
		{
			//So far is the maximum interval that satisfies left-side or right-side pairs
			if(abs((int) travArr[arrNo2][curRightIndex].chrPos - (int) travArr[arrNo1][leftIndex].chrPos) >= (minFragmentSize - finalReadLen))
			{
				//this filters out pairs that are toos close to each other						

				short editCount1 = travArr[arrNo1][leftIndex].edit.numEdits; //didn't want to put this outside the while loop since the missing offset case might cause it to be updated
				short editCount2 = travArr[arrNo2][curRightIndex].edit.numEdits;

				if(editCount1 + editCount2 >=  MISSING_SPLIT_OFFSET)
				{
					if(editCount1 < MISSING_SPLIT_OFFSET || editCount2 < MISSING_SPLIT_OFFSET) //if both are missing don't waste time
					{
						int editNorm1 = editCount1 % MISSING_SPLIT_OFFSET;
						int editNorm2 = editCount2 % MISSING_SPLIT_OFFSET;					

						if(editNorm1 + editNorm2 > bestErrorSoFar)
						{
							if(editNorm1 > bestErrorSoFar)
							{
								break; //this breaks from the innerloop that increases right (since the left guy can't be in a results not worse than the best)
							}
							else
							{
								curRightIndex++;
								continue;
							}
						}
											
						mergedLink leftLink = travArr[arrNo1][leftIndex];
						mergedLink rightLink = travArr[arrNo2][curRightIndex];
						if(IsNoHitExtensionViable(leftLink, rightLink, printReadCode))
						{
							int updatedEditCount1 = leftLink.edit.numEdits;
							int updatedEditCount2 = rightLink.edit.numEdits;
								
							if(updatedEditCount1 + updatedEditCount2 <= bestErrorSoFar)
							{
								if(updatedEditCount1 + updatedEditCount2 < bestErrorSoFar)
								{
									bestErrorSoFar = updatedEditCount1 + updatedEditCount2;
									bufferOffset = startingbufferOffset;							
								}
														
								PrintPairSamFromLinks(foutPair, leftLink, rightLink, printReadCode, 0);
							}
							else
							{
								if(updatedEditCount1 > bestErrorSoFar)
								{
									break;
								}
							}
		
						}
					}
				}
				else
				{
					if(editCount1 + editCount2 <= bestErrorSoFar)
					{
						if(editCount1 + editCount2 < bestErrorSoFar)
						{
							bestErrorSoFar = editCount1 + editCount2;
							bufferOffset = startingbufferOffset;							
						}
						PrintPairSamFromLinks(foutPair, travArr[arrNo1][leftIndex], travArr[arrNo2][curRightIndex], printReadCode, 0);
					}
					else
					{
						if(editCount1 > bestErrorSoFar)
						{
							break; //this breaks from the innerloop that increases right (since the left guy can't be in a results not worse than the best)
						}
					}
				}
			}
			curRightIndex++;
		}
		leftIndex++;
	}
	
	if(bufferOffset > printBufferSize)
	{
		flushBuffer(foutPair);
	}
}

//Unique mapping version -- prints only mappings for reads that don't have more than 1 mapping within the allowed threshold
void MatchMatesAndPrint_UniqueMapping(unsigned char arrNo1, unsigned char arrNo2, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* foutPair, unsigned long long printReadCode)
{
	bool alreadyPrintedFlag = 0;
	int startingbufferOffset = bufferOffset; //we save this since after we start filling in the buffer, if we hit a better mapping, we can just restore the original offest in the buffer and start overwriting

	unsigned int leftIndex = 0, rightBegIndex = 0, leftMateListSize = travArr[arrNo1].size(), rightMateListSize = travArr[arrNo2].size();
	while(leftIndex < leftMateListSize)
	{
		while(rightBegIndex < rightMateListSize && (travArr[arrNo2][rightBegIndex].chrCode < travArr[arrNo1][leftIndex].chrCode || (travArr[arrNo2][rightBegIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][rightBegIndex].chrPos < travArr[arrNo1][leftIndex].chrPos - (maxFragmentSize - finalReadLen) )))
		{
			rightBegIndex++; //increment to position that satisfies left-side upper distance threshold
		}			
		unsigned int curRightIndex = rightBegIndex; //copy index since right beg index will be the starting point of the next left guy
		while(curRightIndex < rightMateListSize &&  travArr[arrNo2][curRightIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][curRightIndex].chrPos <= travArr[arrNo1][leftIndex].chrPos + (maxFragmentSize - finalReadLen))
		{
			//So far is the maximum interval that satisfies left-side or right-side pairs
			if(abs((int) travArr[arrNo2][curRightIndex].chrPos - (int) travArr[arrNo1][leftIndex].chrPos) >= (minFragmentSize - finalReadLen))
			{
				//this filters out pairs that are too close to each other						
				short editCount1 = travArr[arrNo1][leftIndex].edit.numEdits; //didn't want to put this outside the while loop since the missing offset case might cause it to be updated
				short editCount2 = travArr[arrNo2][curRightIndex].edit.numEdits;

				if(editCount1 + editCount2 >=  MISSING_SPLIT_OFFSET)
				{
					if(editCount1 < MISSING_SPLIT_OFFSET || editCount2 < MISSING_SPLIT_OFFSET) //if both are missing don't waste time
					{
						mergedLink leftLink = travArr[arrNo1][leftIndex];
						mergedLink rightLink = travArr[arrNo2][curRightIndex];
						if(IsNoHitExtensionViable(leftLink, rightLink, printReadCode))
						{
							if(alreadyPrintedFlag == 1)
							{
								bufferOffset = startingbufferOffset;
								return;
							}
							else
							{
								alreadyPrintedFlag = 1;
							}
	
							PrintPairSamFromLinks(foutPair, leftLink, rightLink, printReadCode, 0);
						}
					}
				}
				else
				{
					if(alreadyPrintedFlag == 1)
					{
						bufferOffset = startingbufferOffset;
						return;
					}
					else
					{
						alreadyPrintedFlag = 1;
					}

					PrintPairSamFromLinks(foutPair, travArr[arrNo1][leftIndex], travArr[arrNo2][curRightIndex], printReadCode, 0);
				}
			}
			curRightIndex++;
		}
		leftIndex++;
	}

	if(bufferOffset > printBufferSize)
	{
		flushBuffer(foutPair);
	}
}

#define LARGE_BEST_ERROR 999
//Best mapping version -- only one of the mappings within the best tier is reported (selection is arbitrary)
void MatchMatesAndPrint_BestMapping_MapQ(unsigned char arrNo1, unsigned char arrNo2, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* foutPair, unsigned long long printReadCode)
{
	//Clean mapQCounts
	for(int i=0; i<=MAX_MAP_COUNT_TIERS; i++)
	{
		mapQcounts[i]=0;
	}

	short bestErrorSoFar = LARGE_BEST_ERROR; //this acts as SHORT_MAX
	mergedLink bestLeftLink, bestRightLink;

	unsigned int leftIndex = 0, rightBegIndex = 0, leftMateListSize = travArr[arrNo1].size(), rightMateListSize = travArr[arrNo2].size();
	while(leftIndex < leftMateListSize)
	{

		while(rightBegIndex < rightMateListSize && (travArr[arrNo2][rightBegIndex].chrCode < travArr[arrNo1][leftIndex].chrCode || (travArr[arrNo2][rightBegIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][rightBegIndex].chrPos < travArr[arrNo1][leftIndex].chrPos - (maxFragmentSize - finalReadLen) )))
		{
			rightBegIndex++; //increment to position that satisfies left-side upper distance threshold
		}			
		unsigned int curRightIndex = rightBegIndex; //copy index since right beg index will be the starting point of the next left guy
		while(curRightIndex < rightMateListSize &&  travArr[arrNo2][curRightIndex].chrCode == travArr[arrNo1][leftIndex].chrCode && travArr[arrNo2][curRightIndex].chrPos <= travArr[arrNo1][leftIndex].chrPos + (maxFragmentSize - finalReadLen))
		{
			//So far is the maximum interval that satisfies left-side or right-side pairs
			if(abs((int) travArr[arrNo2][curRightIndex].chrPos - (int) travArr[arrNo1][leftIndex].chrPos) >= (minFragmentSize - finalReadLen))
			{
				//this filters out pairs that are toos close to each other						

				short editCount1 = travArr[arrNo1][leftIndex].edit.numEdits; //didn't want to put this outside the while loop since the missing offset case might cause it to be updated
				short editCount2 = travArr[arrNo2][curRightIndex].edit.numEdits;

				if(editCount1 + editCount2 >=  MISSING_SPLIT_OFFSET)
				{
					if(editCount1 < MISSING_SPLIT_OFFSET || editCount2 < MISSING_SPLIT_OFFSET) //if both are missing don't waste time
					{
						int editNorm1 = editCount1 % MISSING_SPLIT_OFFSET;
						int editNorm2 = editCount2 % MISSING_SPLIT_OFFSET;					

						if(editNorm1 + editNorm2 >= bestErrorSoFar)
						{
							if(editNorm1 >= bestErrorSoFar)
							{
								break; //this breaks from the innerloop that increases right (since the left guy can't be in a results not worse than the best)
							}
							else
							{
								curRightIndex++;
								continue;
							}
						}
											
						mergedLink leftLink = travArr[arrNo1][leftIndex];
						mergedLink rightLink = travArr[arrNo2][curRightIndex];
						if(IsNoHitExtensionViable(leftLink, rightLink, printReadCode))
						{
							int updatedEditCount1 = leftLink.edit.numEdits;
							int updatedEditCount2 = rightLink.edit.numEdits;
								
							mapQcounts[updatedEditCount1 + updatedEditCount2]++; 
							if(updatedEditCount1 + updatedEditCount2 < bestErrorSoFar)
							{
								bestErrorSoFar = updatedEditCount1 + updatedEditCount2;
								bestLeftLink = leftLink;
								bestRightLink = rightLink;													
							}
							else
							{
								if(updatedEditCount1 >= bestErrorSoFar)
								{
									break;
								}
							}
						}
					}
				}
				else
				{
					mapQcounts[editCount1 + editCount2]++; 
					if(editCount1 + editCount2 < bestErrorSoFar)
					{
						bestErrorSoFar = editCount1 + editCount2;
						bestLeftLink = travArr[arrNo1][leftIndex];
						bestRightLink = travArr[arrNo2][curRightIndex];													
					}
					else
					{
						if(editCount1 >= bestErrorSoFar)
						{
							break; //this breaks from the innerloop that increases right (since the left guy can't be in a results not worse than the best)
						}
					}
				}
			}
			curRightIndex++;
		}
		leftIndex++;
	}

	if(bestErrorSoFar != LARGE_BEST_ERROR)
	{
		PrintPairSamFromLinks(foutPair, bestLeftLink, bestRightLink, printReadCode, 1);			
	}
}

void ReportSingleEndMappings(unsigned char arrNo, FILE* foutSingle, unsigned long long printReadCode, int mappingMode)
{
	if(splitMode == HALF_SPLIT)
	{
		for(int singleIndex=0; singleIndex<travArr[arrNo].size(); singleIndex++)
		{
			if(travArr[arrNo][singleIndex].edit.numEdits >= MISSING_SPLIT_OFFSET)
			{
				bool viableNoHitExtensionFlag = IsNoHitExtensionViable_SingleEnd(travArr[4][singleIndex], printReadCode);
				if(!viableNoHitExtensionFlag)
				{
					travArr[arrNo][singleIndex].chrCode = 0; //means invalid
					travArr[arrNo][singleIndex].edit.numEdits = 99; //makes it easier for some mapping mode computations below
				}	
			}
		} 

		if(mappingMode == BEST_MAPPING_MODE) //BEST_FAST and BEST_SENSITIVE also arrives at this stage as BEST_MAPPING_MODE, since it doesn't matter when printing
		{
			//Clean mapQCounts
			for(int i=0; i<=MAX_MAP_COUNT_TIERS; i++)
			{
				mapQcounts[i]=0;
			}
			
			//BEST selection computation	
			int lowestErrorSoFar = 99;
			int lowestErrorItemIndex = -1;
			//Search travArr[arrNo] and print the lowest error
			
			for(unsigned int j=0; j<travArr[arrNo].size(); j++)
			{
				mapQcounts[travArr[arrNo][j].edit.numEdits]++;
				if(travArr[arrNo][j].edit.numEdits < lowestErrorSoFar)
				{
					lowestErrorItemIndex = j;
					lowestErrorSoFar = travArr[arrNo][j].edit.numEdits;
					if(lowestErrorSoFar == 0)
					{
						break;
					}
				}
			}

			if(lowestErrorItemIndex != -1)
			{	
				//Print lowest error item here
				PrintSingleSamFromLinks(foutSingle, travArr[arrNo][lowestErrorItemIndex], printReadCode, 1);	
			}
		}
		else if(mappingMode == ALL_MAPPING_MODE)
		{
			int numPrintsRemaining = travArr[arrNo].size();
			if(globalMapCountLimit != 0)
			{
				numPrintsRemaining = globalMapCountLimit;
			} 

			//Print all single items
			for(unsigned int j=0; j<travArr[arrNo].size(); j++)
			{
				if(travArr[arrNo][j].chrCode != 0)
				{
					PrintSingleSamFromLinks(foutSingle, travArr[arrNo][j], printReadCode, 1);
					numPrintsRemaining--;
					if(numPrintsRemaining <= 0)
					{
						break;
					}
				}
			}
		}
		else if(mappingMode == UNIQUE_MAPPING_MODE)
		{
			//Check if the list is larger than 1, but there is only a unique valid mapping
			int uniqueItemIndex = -1; 
			for(unsigned int j=0; j<travArr[arrNo].size(); j++)
			{
				if(travArr[arrNo][j].chrCode != 0)
				{
					if(uniqueItemIndex == -1)
					{
						uniqueItemIndex = j;		
					}
					else
					{
						uniqueItemIndex = -1;
						break;
					}
				}
			} 			

			if(uniqueItemIndex != -1) //Print the unique item if there is only one mapping
			{
				PrintSingleSamFromLinks(foutSingle, travArr[arrNo][uniqueItemIndex], printReadCode, 1);
			}	
		}
		else //mappingMode == STRATUM_MAPPING_MODE
		{
			sort(travArr[arrNo].begin(), travArr[arrNo].end(), sortByErrors); //First sort by errors
			
			if(travArr[arrNo][0].chrCode != 0)
			{
				PrintSingleSamFromLinks(foutSingle, travArr[arrNo][0], printReadCode, 1); //Always print the first item (it's guaranteed that there will be a first item)
				for(unsigned int j=1; j<travArr[arrNo].size(); j++)
				{
					if(travArr[arrNo][j].edit.numEdits != travArr[arrNo][j-1].edit.numEdits) //break whenever the error count changes
					{
						break;
					}
					PrintSingleSamFromLinks(foutSingle, travArr[arrNo][j], printReadCode, 1);
				}	
			}
		}
	}
	else
	{
		if(mappingMode == BEST_MAPPING_MODE) //BEST_FAST and BEST_SENSITIVE also arrives at this stage as BEST_MAPPING_MODE, since it doesn't matter when printing
		{
			//Clean mapQCounts
			for(int i=0; i<=MAX_MAP_COUNT_TIERS; i++)
			{
				mapQcounts[i]=0;
			}
			
			//BEST selection computation	
			int lowestErrorSoFar = 999;
			int lowestErrorItemIndex = 0;
			//Search travArr[arrNo] and print the lowest error
			
			for(unsigned int j=0; j<travArr[arrNo].size(); j++)
			{
				mapQcounts[travArr[arrNo][j].edit.numEdits]++;
				if(travArr[arrNo][j].edit.numEdits < lowestErrorSoFar)
				{
					lowestErrorItemIndex = j;
					lowestErrorSoFar = travArr[arrNo][j].edit.numEdits;
					if(lowestErrorSoFar == 0)
					{
						break;
					}
				}
			}	
			//Print lowest error item here
			PrintSingleSamFromLinks(foutSingle, travArr[arrNo][lowestErrorItemIndex], printReadCode, 1);	
		}
		else if(mappingMode == ALL_MAPPING_MODE)
		{
			int numberOfPrints;
			if(globalMapCountLimit != 0 && globalMapCountLimit < travArr[arrNo].size())
			{
				numberOfPrints = globalMapCountLimit;
			}
			else
			{
				numberOfPrints = travArr[arrNo].size();
			}

			//Print all single items
			for(unsigned int j=0; j<numberOfPrints; j++)
			{
				//REMOVE THIS
				assert(travArr[arrNo][j].edit.numEdits <= 4);
				PrintSingleSamFromLinks(foutSingle, travArr[arrNo][j], printReadCode, 1);
			}
		}
		else if(mappingMode == UNIQUE_MAPPING_MODE)
		{
			if(travArr[arrNo].size() == 1) //Print the unique item if there is only one mapping
			{
				PrintSingleSamFromLinks(foutSingle, travArr[arrNo][0], printReadCode, 1);
			}	
		}
		else //mappingMode == STRATUM_MAPPING_MODE
		{
			sort(travArr[arrNo].begin(), travArr[arrNo].end(), sortByErrors); //First sort by errors
			PrintSingleSamFromLinks(foutSingle, travArr[arrNo][0], printReadCode, 1); //Always print the first item (it's guaranteed that there will be a first item)
			for(unsigned int j=1; j<travArr[arrNo].size(); j++)
			{
				if(travArr[arrNo][j].edit.numEdits != travArr[arrNo][j-1].edit.numEdits) //break whenever the error count changes
				{
					break;
				}
				PrintSingleSamFromLinks(foutSingle, travArr[arrNo][j], printReadCode, 1);
			}	
		}
	}
}
//void MatchMatesAndPrint_BestMapping_MapQ(unsigned char arrNo1, unsigned char arrNo2, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* foutPair, unsigned long long printReadCode)

//Best mapping without inexactSearch
void TraverseLinks_BestMapping_OnlyExactSearch(unsigned int numReads, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* fout, char unmappedSplitsFileName[])
{
	finalReadLen = readLen;
	if(splitMode == HALF_SPLIT)
	{
		finalReadLen = 2 * readLen;
	}

	if(inputMode == SINGLE_MODE)
	{
		if(splitMode == NO_SPLIT)
		{
			travArr[0].resize(1024);
			for(unsigned int i=0; i<linksSize; i++)
			{
				if(links[i].chrCode != 0) //Process if there is a mapping, otherwise move on to the next read since recovery is not possible
				{
					travArr[0].clear();
					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits
					
					//i corresponds to travArr[0] in single end NoSplit mode
					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps_OnlyExact_NoSort(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps_OnlyExact_NoSort(i, 0);
							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
				
					unsigned long long printReadCode = curZero + i + 1; //No division for single-end NoSplitMode 
			
					ReportSingleEndMappings(0, fout, printReadCode, BEST_MAPPING_MODE);
				}
			}	
		}
		else // splitMode == HALF_SPLIT
		{
			travArr[0].resize(1024);
			travArr[1].resize(1024);
			travArr[4].resize(1024);
			
			for(unsigned int i=0; i<linksSize; i+=2)
			{
				short existCount = 0;
			
				if(links[i].chrCode == 0)
				{
					if(links[i].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitSingleEnd(i, unmappedSplitsFileName);
					}
				}
				else
				{
					existCount++;
				}
				
				if(links[i+1].chrCode == 0)
				{
					if(links[i+1].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitSingleEnd(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}

				if(existCount >= 1)
				{
					travArr[0].clear();
					travArr[1].clear();
					travArr[4].clear();
					
					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits

					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps_OnlyExact(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps_OnlyExact(i, 0);

							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
					
					if(links[i+1].chrCode != MEMO_CHR_CODE)
					{	
                        			travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
						ConstructMaps_OnlyExact(i+1, 1);
					}
					else
					{
						int curMemoCode = links[i+1].chrPos;

						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                {
							links[i+1].chrCode = memoInds[curMemoCode].chrCode;
							links[i+1].chrPos = memoInds[curMemoCode].chrPos;
							links[i+1].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
                                                        ConstructMaps_OnlyExact(i+1, 1);
                                                        MemoizeConstructedMaps(1, curMemoCode, links[i+1].dir);
                                               	}
                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
                                                {
                                                        RestoreMemoizedBlock(i+1, 1);
                                                }
					}

					meditEditCountLimit = doubleNumMismatchesPerReadMer; //for merging we're working with full size error limits
					MergeSplitMaps(0, 1);
	
					if(travArr[4].size() != 0)
					{
						unsigned long long printReadCode = (curZero + i)/2 + 1;
						ReportSingleEndMappings(4, fout, printReadCode, BEST_MAPPING_MODE); //Merged travArr[arrNo] is 4
					}
				}
			}
		}
	}
	else //if inputMode == PAIRED_MODE
	{
		if(splitMode == NO_SPLIT)
		{
			travArr[0].resize(1024);
			travArr[2].resize(1024);
		
			for(unsigned int i=0; i<linksSize; i+=2)
			{
				if(links[i].chrCode != 0 && links[i+1].chrCode != 0)
				{
					travArr[0].clear();
					travArr[2].clear();
					
					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits
				
					//i corresponds to travArr[0] in NoSplit mode
					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps_OnlyExact(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps_OnlyExact(i, 0);
							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
					
					//i+1 corresponds to travArr[2] in NoSplit mode
					if(links[i+1].chrCode != MEMO_CHR_CODE)
					{	
                        			travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+1].edit);
						ConstructMaps_OnlyExact(i+1, 2);
					}
					else
					{
						int curMemoCode = links[i+1].chrPos;

						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                {
							links[i+1].chrCode = memoInds[curMemoCode].chrCode;
							links[i+1].chrPos = memoInds[curMemoCode].chrPos;
							links[i+1].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+1].edit);
                                                        ConstructMaps_OnlyExact(i+1, 2);
                                                        MemoizeConstructedMaps(2, curMemoCode, links[i+1].dir);
                                               	}
                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
                                                {
                                                        RestoreMemoizedBlock(i+1, 2);
                                                }
					}

					unsigned long long printReadCode = (curZero + i)/2 + 1; //Divide by 2 for NoSplitMode 
					MatchMatesAndPrint_BestMapping_MapQ(0, 2, minFragmentSize, maxFragmentSize, fout, printReadCode);
				}
			}
		}
		else //splitMode == HALF_SPLIT
		{
			travArr[0].resize(1024);
			travArr[1].resize(1024);
			travArr[2].resize(1024);
			travArr[3].resize(1024);
			travArr[4].resize(1024);
			travArr[5].resize(1024);

			for(unsigned int i=0; i<linksSize; i+=4)
			{
				short existCount = 0;

				if(links[i].chrCode == 0)
				{
					if(links[i].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName);
					}
				}
				else
				{
					existCount++;
				}
				
				if(links[i+1].chrCode == 0)
				{
					if(links[i+1].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}

				if(links[i+2].chrCode == 0)
				{
					if(links[i+2].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}
				
				if(links[i+3].chrCode == 0)
				{
					if(links[i+3].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}

				if(existCount >= 3)
				{
					travArr[0].clear();
					travArr[1].clear();
					travArr[4].clear();

					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits

					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps_OnlyExact(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps_OnlyExact(i, 0);

							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
					
					if(links[i+1].chrCode != MEMO_CHR_CODE)
					{	
                        			travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
						ConstructMaps_OnlyExact(i+1, 1);
					}
					else
					{
						int curMemoCode = links[i+1].chrPos;

						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                {
							links[i+1].chrCode = memoInds[curMemoCode].chrCode;
							links[i+1].chrPos = memoInds[curMemoCode].chrPos;
							links[i+1].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
                                                        ConstructMaps_OnlyExact(i+1, 1);
                                                        MemoizeConstructedMaps(1, curMemoCode, links[i+1].dir);
                                               	}
                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
                                                {
                                                        RestoreMemoizedBlock(i+1, 1);
                                                }
					}

					meditEditCountLimit = doubleNumMismatchesPerReadMer; //for merging we're working with full size error limits
					MergeSplitMaps(0, 1);

					if(travArr[4].size() > 0) //this is the meregd list for the first mate
					{
						travArr[2].clear();
                        			travArr[3].clear();
                        			travArr[5].clear();
					
						meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits
						
						if(links[i+2].chrCode != MEMO_CHR_CODE)
						{
                        				travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+2].edit);
							ConstructMaps_OnlyExact(i+2, 2);
						}
						else
						{
							int curMemoCode = links[i+2].chrPos;

							if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                	{
								links[i+2].chrCode = memoInds[curMemoCode].chrCode;
								links[i+2].chrPos = memoInds[curMemoCode].chrPos;
								links[i+2].edit = memoInds[curMemoCode].edit;
                        					travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+2].edit);
                                                	        ConstructMaps_OnlyExact(i+2, 2);
	                                                        MemoizeConstructedMaps(2, curMemoCode, links[i+2].dir);
							}
	                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
	                                                {
	                                                        RestoreMemoizedBlock(i+2, 2);
	                                                }
						}	

						if(links[i+3].chrCode != MEMO_CHR_CODE)
						{
                        				travArrReadmerLength[3] = readLen + GetIndelLengthModifier(links[i+3].edit);
							ConstructMaps_OnlyExact(i+3, 3);
						}
						else
						{
							int curMemoCode = links[i+3].chrPos;
							
							if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                	{
								links[i+3].chrCode = memoInds[curMemoCode].chrCode;
								links[i+3].chrPos = memoInds[curMemoCode].chrPos;
								links[i+3].edit = memoInds[curMemoCode].edit;
                        					travArrReadmerLength[3] = readLen + GetIndelLengthModifier(links[i+3].edit);
                                                	        ConstructMaps_OnlyExact(i+3, 3);
	                                                        MemoizeConstructedMaps(3, curMemoCode, links[i+3].dir);
							}
	                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
	                                                {
	                                                        RestoreMemoizedBlock(i+3, 3);
	                                                }
						}						
						
						meditEditCountLimit = doubleNumMismatchesPerReadMer; //for merging we're working with full size error limits
						MergeSplitMaps(2,3);

						if(travArr[5].size() > 0) //this the merged list for second mate
						{
							//this is the number of the read that will appear in the output (divided by 4 instead of 8, since different from readCodes, linkIndex doesn't store directionality
							unsigned long long printReadCode = (curZero + i)/4 + 1; 

							MatchMatesAndPrint_BestMapping_MapQ(4, 5, minFragmentSize, maxFragmentSize, fout, printReadCode);
						}
					}
				}
			}	
		}
	} 
}

//All mapping version of traverse links -- where all mapping constructions, split merging as well as pair matching is done
void TraverseLinks(unsigned int numReads, unsigned int minFragmentSize, unsigned int maxFragmentSize, FILE* fout, char unmappedSplitsFileName[])
{
	finalReadLen = readLen;
	if(splitMode == HALF_SPLIT)
	{
		finalReadLen = 2 * readLen;
	}

	if(inputMode == SINGLE_MODE)
	{
		if(splitMode == NO_SPLIT)
		{
			travArr[0].resize(1024);
			for(unsigned int i=0; i<linksSize; i++)
			{
				if(links[i].chrCode != 0) //Process if there is a mapping, otherwise move on to the next read since recovery is not possible
				{
					travArr[0].clear();
					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits
					
					//i corresponds to travArr[0] in single end NoSplit mode
					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps_NoSort(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps_NoSort(i, 0);
							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
				
					unsigned long long printReadCode = curZero + i + 1; //No division for single-end NoSplitMode 
	
					if(mappingMode == BEST_SENSITIVE_MAPPING_MODE)
					{
						ReportSingleEndMappings(0, fout, printReadCode, BEST_MAPPING_MODE);
					}
					else
					{
						ReportSingleEndMappings(0, fout, printReadCode, mappingMode);
					}
				}
			}	
		}
		else // splitMode == HALF_SPLIT
		{
			travArr[0].resize(1024);
			travArr[1].resize(1024);
			travArr[4].resize(1024);
			
			for(unsigned int i=0; i<linksSize; i+=2)
			{
				short existCount = 0;
			
				if(links[i].chrCode == 0)
				{
					if(links[i].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitSingleEnd(i, unmappedSplitsFileName);
					}
				}
				else
				{
					existCount++;
				}
				
				if(links[i+1].chrCode == 0)
				{
					if(links[i+1].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitSingleEnd(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}

				if(existCount >= 1)
				{
					travArr[0].clear();
					travArr[1].clear();
					travArr[4].clear();
					
					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits

					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps(i, 0);

							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
					
					if(links[i+1].chrCode != MEMO_CHR_CODE)
					{	
                        			travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
						ConstructMaps(i+1, 1);
					}
					else
					{
						int curMemoCode = links[i+1].chrPos;

						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                {
							links[i+1].chrCode = memoInds[curMemoCode].chrCode;
							links[i+1].chrPos = memoInds[curMemoCode].chrPos;
							links[i+1].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
                                                        ConstructMaps(i+1, 1);
                                                        MemoizeConstructedMaps(1, curMemoCode, links[i+1].dir);
                                               	}
                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
                                                {
                                                        RestoreMemoizedBlock(i+1, 1);
                                                }
					}

					meditEditCountLimit = doubleNumMismatchesPerReadMer; //for merging we're working with full size error limits
					MergeSplitMaps(0, 1);

					if(travArr[4].size() != 0)
					{
						unsigned long long printReadCode = (curZero + i)/2 + 1;
						if(mappingMode == BEST_SENSITIVE_MAPPING_MODE)
						{
							ReportSingleEndMappings(4, fout, printReadCode, BEST_MAPPING_MODE);
						}
						else
						{
							ReportSingleEndMappings(4, fout, printReadCode, mappingMode);
						}
					}
				}
			}

		}
	}
	else //if inputMode == PAIRED_MODE
	{
		if(splitMode == NO_SPLIT)
		{
			travArr[0].resize(1024);
			travArr[2].resize(1024);
		
			for(unsigned int i=0; i<linksSize; i+=2)
			{
				if(links[i].chrCode != 0 && links[i+1].chrCode != 0)
				{
					travArr[0].clear();
					travArr[2].clear();
					
					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits
				
					//i corresponds to travArr[0] in NoSplit mode
					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps(i, 0);
							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			
					
					//i+1 corresponds to travArr[2] in NoSplit mode
					if(links[i+1].chrCode != MEMO_CHR_CODE)
					{	
                        			travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+1].edit);
						ConstructMaps(i+1, 2);
					}
					else
					{
						int curMemoCode = links[i+1].chrPos;

						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                {
							links[i+1].chrCode = memoInds[curMemoCode].chrCode;
							links[i+1].chrPos = memoInds[curMemoCode].chrPos;
							links[i+1].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+1].edit);
                                                        ConstructMaps(i+1, 2);
                                                        MemoizeConstructedMaps(2, curMemoCode, links[i+1].dir);
                                               	}
                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
                                                {
                                                        RestoreMemoizedBlock(i+1, 2);
                                                }
					}

					unsigned long long printReadCode = (curZero + i)/2 + 1; //Divide by 2 for NoSplitMode 
					
					if(mappingMode == BEST_SENSITIVE_MAPPING_MODE)
					{
						MatchMatesAndPrint_BestMapping_MapQ(0, 2, minFragmentSize, maxFragmentSize, fout, printReadCode);
					}
					else if(mappingMode == ALL_MAPPING_MODE)
					{
						if(globalMapCountLimit == 0)
						{
							MatchMatesAndPrint(0, 2, minFragmentSize, maxFragmentSize, fout, printReadCode);
						}
						else
						{
							MatchMatesAndPrint_MapCountLimit(0, 2, minFragmentSize, maxFragmentSize, fout, printReadCode, globalMapCountLimit);
						}
					}
					else if(mappingMode == STRATUM_MAPPING_MODE)
					{
						MatchMatesAndPrint_StratumMapping(0, 2, minFragmentSize, maxFragmentSize, fout, printReadCode);
					}
					else// if(mappingMode == UNIQUE_MAPPING_MODE)
					{
						MatchMatesAndPrint_UniqueMapping(0, 2, minFragmentSize, maxFragmentSize, fout, printReadCode);
					}
				}
			}
		}
		else //splitMode == HALF_SPLIT
		{
			travArr[0].resize(1024);
			travArr[1].resize(1024);
			travArr[2].resize(1024);
			travArr[3].resize(1024);
			travArr[4].resize(1024);
			travArr[5].resize(1024);

			for(unsigned int i=0; i<linksSize; i+=4)
			{
				short existCount = 0;
			
				if(links[i].chrCode == 0)
				{
					if(links[i].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName);
					}
				}
				else
				{
					existCount++;
				}
				
				if(links[i+1].chrCode == 0)
				{
					if(links[i+1].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}

				if(links[i+2].chrCode == 0)
				{
					if(links[i+2].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}
				
				if(links[i+3].chrCode == 0)
				{
					if(links[i+3].chrPos == 0)
					{
						//Reached to the end of the Nohit list (or maybe didn't start yet)
						ReloadNoHitArray_SplitPaired(i, unmappedSplitsFileName); //intentionally sending the splitID of the first split of the first mate
					}
				}
				else
				{
					existCount++;
				}
			
				if(existCount >= 3)
				{
					travArr[0].clear();
					travArr[1].clear();
                        		travArr[4].clear();

					meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits
				
					if(links[i].chrCode != MEMO_CHR_CODE)
					{
                        			travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
						ConstructMaps(i, 0); //Construct as regular
					}	
					else
					{
						int curMemoCode = links[i].chrPos;
						
						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
						{
							links[i].chrCode = memoInds[curMemoCode].chrCode;
							links[i].chrPos = memoInds[curMemoCode].chrPos;
							links[i].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[0] = readLen + GetIndelLengthModifier(links[i].edit);
							ConstructMaps(i, 0);
							MemoizeConstructedMaps(0, curMemoCode, links[i].dir);
						}
						else //Then this is the following readID for the readMer (copy memoizedBlockDirectly)
						{
							RestoreMemoizedBlock(i, 0);
						}
					}			

					if(links[i+1].chrCode != MEMO_CHR_CODE)
					{	
                        			travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
						ConstructMaps(i+1, 1);

					}
					else
					{
						int curMemoCode = links[i+1].chrPos;

						if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                {
							links[i+1].chrCode = memoInds[curMemoCode].chrCode;
							links[i+1].chrPos = memoInds[curMemoCode].chrPos;
							links[i+1].edit = memoInds[curMemoCode].edit;
                        				travArrReadmerLength[1] = readLen + GetIndelLengthModifier(links[i+1].edit);
                                                        ConstructMaps(i+1, 1);
                                                        MemoizeConstructedMaps(1, curMemoCode, links[i+1].dir);
                                               	}
                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
                                                {
                                                        RestoreMemoizedBlock(i+1, 1);
                                                }
					}

					meditEditCountLimit = doubleNumMismatchesPerReadMer; //for merging we're working with full size error limits
					MergeSplitMaps(0, 1);

					if(travArr[4].size() > 0) //this is the meregd list for the first mate
					{
						travArr[2].clear();
                        			travArr[3].clear();
                        			travArr[5].clear();
					
						meditEditCountLimit = numMismatchesPerReadMer; //for construction we're working with split size error limits

						if(links[i+2].chrCode != MEMO_CHR_CODE)
						{
                        				travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+2].edit);
							ConstructMaps(i+2, 2);
						}
						else
						{
							int curMemoCode = links[i+2].chrPos;

							if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                	{
								links[i+2].chrCode = memoInds[curMemoCode].chrCode;
								links[i+2].chrPos = memoInds[curMemoCode].chrPos;
								links[i+2].edit = memoInds[curMemoCode].edit;
                        					travArrReadmerLength[2] = readLen + GetIndelLengthModifier(links[i+2].edit);
                                                	        ConstructMaps(i+2, 2);
	                                                        MemoizeConstructedMaps(2, curMemoCode, links[i+2].dir);
							}
	                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
	                                                {
	                                                        RestoreMemoizedBlock(i+2, 2);
	                                                }
						}	

						if(links[i+3].chrCode != MEMO_CHR_CODE)
						{
                        				travArrReadmerLength[3] = readLen + GetIndelLengthModifier(links[i+3].edit);
						
							ConstructMaps(i+3, 3);
						}
						else
						{
							int curMemoCode = links[i+3].chrPos;
							
							if(memoInds[curMemoCode].begInd == -1) //Then it hasn't been created yet
                                                	{
								links[i+3].chrCode = memoInds[curMemoCode].chrCode;
								links[i+3].chrPos = memoInds[curMemoCode].chrPos;
								links[i+3].edit = memoInds[curMemoCode].edit;
                        					travArrReadmerLength[3] = readLen + GetIndelLengthModifier(links[i+3].edit);
                                                	        ConstructMaps(i+3, 3);
	                                                        MemoizeConstructedMaps(3, curMemoCode, links[i+3].dir);
							}
	                                                else //Then this is the following readID for the readMer (copy memoizedBlockDirectly
	                                                {
	                                                        RestoreMemoizedBlock(i+3, 3);
	                                                }
						}						

						meditEditCountLimit = doubleNumMismatchesPerReadMer; //for merging we're working with full size error limits
						MergeSplitMaps(2,3);

						if(travArr[5].size() > 0) //this the merged list for second mate
						{
							//this is the number of the read that will appear in the output (divided by 4 instead of 8, since different from readCodes, linkIndex doesn't store directionality
							unsigned long long printReadCode = (curZero + i)/4 + 1; 
							if(mappingMode == BEST_SENSITIVE_MAPPING_MODE)
							{
								MatchMatesAndPrint_BestMapping_MapQ(4, 5, minFragmentSize, maxFragmentSize, fout, printReadCode);
							}
							else if(mappingMode == ALL_MAPPING_MODE)
							{
								if(globalMapCountLimit == 0)
								{
									MatchMatesAndPrint(4, 5, minFragmentSize, maxFragmentSize, fout, printReadCode);
								}
								else
								{
									MatchMatesAndPrint_MapCountLimit(4, 5, minFragmentSize, maxFragmentSize, fout, printReadCode, globalMapCountLimit);
								}
							}
							else if(mappingMode == STRATUM_MAPPING_MODE)
							{
								MatchMatesAndPrint_StratumMapping(4, 5, minFragmentSize, maxFragmentSize, fout, printReadCode);
							}
							else// if(mappingMode == UNIQUE_MAPPING_MODE)
							{
								MatchMatesAndPrint_UniqueMapping(4, 5, minFragmentSize, maxFragmentSize, fout, printReadCode);
							}
						}
					}
				}
			}	
		}
	} 
}

inline bool GetDirectionFromReadCode(unsigned long long code)
{
	return code % 2; //Provided that the alphabet size is even (if the alphabet is changed, this needs to be changed too)
}
				
unsigned long long GetReadIdFromString(const string& str)
{
	unsigned long long val = 0;
	unsigned short len = str.length();
	for(int i=0; i<len; i++)
	{
		val *= IDENTITY_ALPHABET_SIZE; 
		val += (str[i] - IDENTITY_ALPHABET_START);
	}
	return val;
}

unsigned long long GetReadIdFromStringWithPos(const string& str, unsigned int pos)
{
	unsigned long long val = 0;
	unsigned short len = pos + idDigitLen;
	for(unsigned int i=pos; i<len; i++)
	{
		val *= IDENTITY_ALPHABET_SIZE; 
		val += (str[i] - IDENTITY_ALPHABET_START);
	}
	return val;
}

//Obtains equivClass index from posHub, and also determines the proper offset and direction from hit to rep
//Retrieve offset is guaranteed to return correct results since the direction and position is always stored within the list
//BlockOffset is always positive (that's why dir is required as well)
unsigned int RetrieveOffsetAndDirInfo(unsigned char linkChrCode, unsigned int linkChrPos, short& equivBlockOffset, bool& dirFromHitToRep)
{
	unsigned int curEquivIndex = 0;
        std::tr1::unordered_map<unsigned int, unsigned int>::iterator ite = equivHub[(short)linkChrCode].find(linkChrPos);
        if(ite != equivHub[(short)linkChrCode].end())
        {
                curEquivIndex = ite->second;
        }

	if(curEquivIndex == 0)
	{
		return 0; //meaning that there is none;
	}
	
	equivClassNode& curEqNode = equivClassList[curEquivIndex];
	unsigned int curBlockLen = curEqNode.blockLen;
	unsigned int curRepPos = curEqNode.chrPos;
	
	if(linkChrCode == curEqNode.chrCode && linkChrPos >= curRepPos && linkChrPos < curRepPos + curBlockLen)
	{
		dirFromHitToRep = 0;
		equivBlockOffset = linkChrPos - curRepPos;
		return curEquivIndex;
	}
	else
	{
		unsigned int curListSize = curEqNode.listSize;
		for(unsigned int k=0; k<curListSize; k++)
		{
			eqItem& curHomItem = curEqNode.list[k];
			unsigned int curHomPos = curHomItem.chrPos;
			if(linkChrCode == curHomItem.chrCode)
			{
				if(curHomItem.dir == 0)
				{
					if(linkChrPos >= curHomPos && linkChrPos < curHomPos + curBlockLen)
					{
						dirFromHitToRep = 0;
						equivBlockOffset = linkChrPos - curHomPos;
						return curEquivIndex;
					}	
				}
				else
				{
					if(linkChrPos <= curHomPos && linkChrPos > curHomPos - curBlockLen)
					{
						dirFromHitToRep = 1;
						equivBlockOffset = curHomPos - linkChrPos;
						return curEquivIndex;
					}
				}
			}
		}
	}

	assert(0);
}

int main(int argc, char* argv[])
{
	double beginTime = getTime();
	if(argc!=30)
	{
		cout << "This program traverses the links table and homology table simultaneously in order infer final mapping results" << endl;	
		cout << "ARGV[1] homology table perfect map (equivalence class) compact binary file" << endl;
		cout << "ARGV[2] homology table inexact map (E2) compact binary file" << endl;
		cout << "ARGV[3] Read Length" << endl;
		cout << "ARGV[4] UniRef full sequencefile (templateAB), fai should have been created as well" << endl;
		cout << "ARGV[5] Line length of the uniref file" << endl;
		cout << "ARGV[6] UNSORTED combined links table including perfs but without terrible" << endl;
		cout << "ARGV[7] Distance metric (HAMMING or LEVENSHTEIN)" << endl;
		cout << "ARGV[8] mapping domain mode: ALL | STRATUM | BEST | UNIQUE" << endl;
		cout << "ARGV[9] UNUSED" << endl;
		cout << "ARGV[10] Maximum printed mappings (Only used for ALL mapping)" << endl;
		cout << "ARGV[11] output file that represents all mappings with reads (not in SAM format yet)" << endl;
		cout << "ARGV[12] sam output type: COLLAPSED | UNCOLLAPSED" << endl;
		cout << "ARGV[13] UNUSED" << endl;
		cout << "ARGV[14] UNUSED" << endl;
		cout << "ARGV[15] readIdentityDigitLength" << endl;
		cout << "ARGV[16] splitMode" << endl; //curently HALF or NONE
		cout << "ARGV[17] inputMode" << endl;
		cout << "ARGV[18] minFragmentSize" << endl;
		cout << "ARGV[19] maxFragmentSize" << endl;
		cout << "ARGV[20] INPUT: unmapped splits List (SAM format) " << endl;
		cout << "ARGV[21] INPUT: longNameList" << endl;
		cout << "ARGV[22] MODE: Format of the longNameList: FA (FASTA) or SAM (SAM)" << endl;
		cout << "ARGV[23] INFO: Total number of reads in the original dataset (two-splits and mate-pairs are reported only once)" << endl;
		cout << "ARGV[24] PARAM: Memoization Threshold - Any links with at least this many read-mers will be memoized" << endl;
		cout << "ARGV[25] PARAM: Num mismatches per Readmer" << endl;
		cout << "ARGV[26] PARAM: Num mismatches in the inexact homology table" << endl;
		cout << "ARGV[27] FLAG: For printing original read names, 0 or 1, 1 prints original read names, 0 prints the order of the read" << endl;
		cout << "ARGV[28] FLAG: For printing quality scores in the sam output, 0 or 1, 1 prints quality scores, 0 prints '*' in the quality field" << endl;
		cout << "ARGV[29] INPUT: Read File List File" << endl;
		exit(16);
	}	

	SetupRevCompChar();
	InitIndelAlignmentMatrix();

	numMismatchesPerReadMer = atoi(argv[25]);
	doubleNumMismatchesPerReadMer = 2 * numMismatchesPerReadMer;
	numHomTableMismatches = atoi(argv[26]);

	MEMOIZATION_THRESHOLD = atoi(argv[24]);

	if(string(argv[7]) == "LEVENSHTEIN")
	{
		distanceMetric = LEVENSHTEIN;
	}
	else // HAMMING
	{
		distanceMetric = HAMMING;
	}

	ifstream finLongNameList(argv[21]);
	SetupLongNameList(finLongNameList);

	unsigned long long numTotalReads = atol(argv[23]);
	idDigitLen = (unsigned char) atoi(argv[15]);
	readLen = atoi(argv[3]);
	itoa10(readLen, readLenStr);
		
	string noHitModeStr(argv[22]);
	if(noHitModeStr == "FA")
	{
		noHitMode = FASTA;
	}
	else if(noHitModeStr == "SAM")
	{
		noHitMode = SAM;
	}
	else
	{
		cout << "noHitMode should either be FA (fasta) or SAM" << endl;
		exit(23);
	}

	SetupNoHitArr(numTotalReads);
	SetUpValStr(); //this creates char array representations of numbers upto a specific value (it's for faster printing);

	minFragmentSize = atoi(argv[18]);
	maxFragmentSize = atoi(argv[19]);

	string splitModeStr = argv[16];
	if(splitModeStr == "NONE" || splitModeStr == "FULL")
	{
		splitMode = NO_SPLIT;
	}
	else if(splitModeStr == "HALF")
	{
		splitMode = HALF_SPLIT;
	}
	else
	{
		cout << "ERROR: splitMode not recognized" << endl;
		exit(6);
	}

	string inputModeStr = argv[17];
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
		cout << "ERROR: inputMode not recognized" << endl;
		exit(5);
	}

	meditEditCountLimit = numMismatchesPerReadMer;

	int refLineLen = atoi(argv[5]);
	string mappingModeStr(argv[8]);

	if(mappingModeStr == "ALL")
	{
		mappingMode = ALL_MAPPING_MODE;
		globalMapCountLimit = atoi(argv[10]);
	}
	else if(mappingModeStr == "BEST") //This is the faster best mapping mode
	{
		mappingMode = BEST_MAPPING_MODE;
		mappingQualityPrintFlag = 1;
	}
	else if(mappingModeStr == "BEST_SENSITIVE")
	{	
		mappingMode = BEST_SENSITIVE_MAPPING_MODE;
		mappingQualityPrintFlag = 1;
	}
	else if(mappingModeStr == "BEST_FAST")
	{
		mappingMode = BEST_FAST_MAPPING_MODE;
		mappingQualityPrintFlag = 1;
	}
	else if(mappingModeStr == "STRATUM")
	{
		mappingMode = STRATUM_MAPPING_MODE; 
	}
	else if(mappingModeStr == "UNIQUE")
	{
		mappingMode = UNIQUE_MAPPING_MODE;
	}
	else
	{
		cout << "ERROR: Unknown mapping mode: '" << mappingModeStr << "'" << endl;
		exit(898);
	}
	
	// In the case of hamming distance BEST mode is BEST_SENSITIVE MODE that would load the memory into table
	if(distanceMetric == HAMMING && mappingMode == BEST_MAPPING_MODE)
	{
		mappingMode = BEST_SENSITIVE_MAPPING_MODE;
	}

	//Settings for read name and quality score printing
	printReadNamesFlag = atoi(argv[27]);
	printQualityScoresFlag = atoi(argv[28]);
	if(printQualityScoresFlag)
		printReadNamesFlag = 1; //Quality score printing automatically activates read name printing
	if(printReadNamesFlag || printQualityScoresFlag)
	{
		ifstream finRFL(argv[29]); //readFileList
		assert(finRFL.is_open());
		InitFastqs(finRFL);
	}
	
	//Load all needed files
	LoadMultiChrReference(argv[4], refLineLen);
	LoadHomologyTables(argv[1], argv[2], readLen, numHomTableMismatches, mappingMode);

	double refHomTime = getTime();
	cout << "Loaded reference and homology table in " << refHomTime - beginTime << " seconds..." << endl;

	bool allLinksTraversedFLAG = 0;
	unsigned short numTraversal = 0;

	while(allLinksTraversedFLAG == 0)
	{
		allLinksTraversedFLAG = LoadLinksRepositoryFromUnsorted(argv[6], numTotalReads, numTraversal);

		double postLoadTime = getTime();

		cout << "Loading links finished after " << postLoadTime - refHomTime << " seconds..." << endl;

		FILE* fout = fopen(argv[11], "wb");
		setvbuf(fout, NULL, _IOFBF, 4096);

		finalReadLen = readLen;
		if(splitMode == HALF_SPLIT)
		{
			finalReadLen = 2 * readLen;
		}

		if(mappingMode == BEST_MAPPING_MODE || mappingMode == BEST_FAST_MAPPING_MODE)
		{
			TraverseLinks_BestMapping_OnlyExactSearch(numTotalReads, minFragmentSize, maxFragmentSize, fout, argv[20]); //argv[20] is the file name of unmapped splits
		}
		else
		{
			//all - best_sensitive - unique - stratum 
			TraverseLinks(numTotalReads, minFragmentSize, maxFragmentSize, fout, argv[20]); //argv[20] is the file name of unmapped splits
		}
	
		//This prints anything that remains and resets buffer
		flushBuffer(fout);

		double timeAfterTraversal = getTime();
		cout << "Traversed Links: " << timeAfterTraversal - postLoadTime << " seconds..." << endl;

		numTraversal = 0;
	}

	cout << "MemoIndsSize: " << memoInds.size() << endl;
	cout << "MemoArrSize: " << memoArr.size() << endl;

	return 0;
}


