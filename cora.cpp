/*Copyright (c) 2015-2016 Deniz Yorukoglu. All rights reserved.*/
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<assert.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<stdio.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<vector>
using namespace std;

#define VERSION "1.1.5b"

#define HOMINDEX "coraIndex"
#define MAPPERINDEX "mapperIndex"
#define SEARCH "search"
#define READFILEGEN "readFileGen"
#define FAIGENERATE "faiGenerate"

#define NUMCHROMFORSHORT 127 //If the number of chroms in the fai file exceeds this file then large chrom version will be called with "_L"
int numChromInFai; //obtained from Fai file 

double getTime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec+t.tv_usec/1000000.0;
}

//Variables for homology table construction
int homTableConstructionFlag; //This flag determines whether a homology table will be constructed
string homTablePerfectMatchFile = "./Exact"; //This field contains the compact equivalence infomation
string homTableInexactMatchFile = "./Inexact"; //This field contains the prefix for the binary files containing the inexact mappings
int numInexactHomTablePartitions; //Number of partitions (starting from 1) that should be added as a suffix to the prefix above
string PhysicalSplitsFile = "__TEMP__Aux_Physical_Splits_File"; //These splits only contain 2-char signals (they can be grouped) for physical partition of the processes through time for efficient memory management
string ParallelSplitsFile = "__TEMP__Aux_Parallel_Splits_File"; //These splits contain 4-char signals for partitioning into threads
int numPhysicalSplits = 10;
int numParallelSplits = 8;

//Parameters and variables for FASTQ splitting
#define AUTO_FASTQ_SPLIT_MODE 0 // This is the flag to set the number of fastq splits to be ceil( refLen / split_function_interval )
#define AUTO_FASTQ_SPLIT_FUNCTION_INTERVAL 100000000
#define AUTO_FASTQ_SPLIT_FUNCTION_INTERVAL_LOWMEM 66666666
#define MAX_AUTO_FASTQ_SPLIT 144
int numFastqSplits = AUTO_FASTQ_SPLIT_MODE;

// Executable names
string collapseExec("collapse");
string linkConstructExec("linkConstruct");
string homTableSetupExec("homTable_setup");
string mappingInferenceExec("mappingInference");
string fastqSplitterExec("fastqSplit");
///////////////////

// Cora main parameters
string RepresentativeMapperExecutable = "bwa"; //either the mapping executable for present aligners, or the full mapping mode for alignment 
string RepresentativeMappingMode = "BWA";   // MRSFAST, MRSFAST_ULTRA, BWA, BWA_MEM, BOWTIE, BOWTIE_2, MANUAL --> manual requires the exact command to run the mapper
string runMode = "1111";
string fastqInputListFile = "";
string refFile = "", refFileFai = "";
string mappingOutputFile = "CORA_Output.sam";
string mappingReportMode = "ALL"; //ALL or BEST or UNIQUE or STRATUM or BEST_SENSITIVE or BEST_FAST (only for gapped mapping)
string mapperMetric = "HAMMING";
string inputMode = "PAIRED"; //SINGLE or PAIRED
int insertSizeLowerThreshold = 150, insertSizeUpperThreshold = 650;
string splitMode = "HALF"; //FULL or HALF (This might be made more parametric later on) -- Can also be THREEWAY
string collapseMode = "WITHREF"; //NOREF or WITHREF
string readCompressionMode = "OFF"; //OFF or GZIP
string temporaryDirectoryName = "__temporary_CORA_files";
string inputReadGroupString = "";
string killSignalFile = "__TEMP__Kill_Signal_File"; //This file is created if a subprogram is killing itself intentionally (for deallocation speed) so that the parent can ignore it
string collapseFragmentFile = "__TEMP__Aux_Fastq_Splits_File"; //either NONE or a file containing list of fragment groups (AA_GT_, TA_, N_AG_AC_, etc.) -- This is for Fastq Splitting -- It is stuomatically converted to NONE if the number of splits is specified as less than two
int numCoarseMappingThreads = 1;
int memoizationThreshold = 20; //memoization threshold show the lower inclusive bound of joint-readMer reads to be memoized.
int numMismatchesPerReadMer = 2;
int globalMapCountLimit = 0; //If non-zero determines the maximum number of mappings to be printed (only for ALL mapping)

//-=================================================

//Encoded read name (identity) related parameters
#define IDENTITY_ALPHABET_START 33
#define IDENTITY_ALPHABET_END   126
#define IDENTITY_ALPHABET_SIZE  94
long long CapOnTotalNumberOfReads; //This value determines how many chars will be used to identify each read item -- mates counted separately
long long identitySpaceUpperLimit; //This is a derived vale from total number Of reads depending on the extra bits each read will be using (for example directionality, splits etc)

//Length related parameters
char idDigitLen; //This is a derived value from identitySpaceUpperLimit and the 
int refLineLen; //lineLength if reference
int readLen; //full readLength for each read 
int completeReadLen; //regardless of splitMode this is the full length of the read
int kmerLen; //length for the homology table

//Item count related parameters
int numSamples; //number of individual datasets
int maxChrSizeInRef;
unsigned long long numTotalReads;

//These are precomputed LPT (longest processing time) load balancing tables for 2-base and 4-base signals to be used for splitting the k-mer space
string LPT2path("LPT_2.dat");
string LPT4path("LPT_4.dat");

string SplitTwoSignal(int numSplits)
{
    if(numSplits < 1 || numSplits > 10)
    {
        cout << "ERROR: Number of physical splits can be between 1 and 10" << endl;
        exit(89);
    }

    ifstream fin(LPT2path.c_str());

    string twoArr;
    for(int i=1; i<=numSplits; i++)
    {
        getline(fin, twoArr);
    }

    for(int i=0; i<twoArr.length(); i++)
    {
        if(twoArr[i] == 'X')
        {
            twoArr[i] = '\n';
        }
    }

    fin.clear();
    fin.close();

    return twoArr;
}

string SplitFourSignal(int numSplits)
{
    if(numSplits < 1 || numSplits > 144)
    {
        cout << "ERROR: Number of physical splits can be between 1 and 144" << endl;
        exit(89);
    }

    ifstream fin(LPT4path.c_str());

    string fourArr;
    for(int i=1; i<=numSplits; i++)
    {
        getline(fin, fourArr);
    }
    
    for(int i=0; i<fourArr.length(); i++)
    {
        if(fourArr[i] == 'X')
        {
            fourArr[i] = '\n';
        }
    }

    fin.clear();
    fin.close();

    return fourArr;
}


ofstream logOut;

#define IGNORE 17
#define VERIFY 34

//This executes a given command, while tracking runtime and exitcodes
void ExecuteSystemCall(const stringstream& call, char ignoreFlag)
{
    cout << "Command: " << call.str() << endl;
    logOut << "Command: " << call.str() << endl;

    double callBeginTime = getTime();
    int x = system(call.str().c_str());
    double callEndTime = getTime();

    if(ignoreFlag != IGNORE)
    {
        if(x != 256 && x != 0)
        {
            cout << "WARNING: last command returned non-zero exit code: " << x << endl;
            cout << "*** CORA run might have crashed here and the rest of pipeline may not run properly ***" << endl;
        }
    }
    else
    {
        //Check kill signal
        FILE* finKill = fopen(killSignalFile.c_str(), "r");
        
        if(!finKill)
        {
            cout << "WARNING: last command returned non-zero exit code: " << x << endl;
        }
        else // delete killSignalFile
        {
            fclose(finKill);
            stringstream removeKillFile;
            removeKillFile << "rm " << killSignalFile << endl;
            ExecuteSystemCall(removeKillFile, VERIFY);  
        }
    }

    cout << "Exit Code: " << x << "  --- Completed... (in " << callEndTime - callBeginTime << " seconds)" << endl;
    logOut << "runTime: " << callEndTime - callBeginTime << endl;
}

//Fai parsing
int CheckFai_and_ReturnMaxChrSize(string faiFileName)
{
    cout << "FaiFileName: " << faiFileName << endl;
    assert(faiFileName.substr(faiFileName.length()-4,4) == ".fai");
    string actualFileName = faiFileName.substr(0, faiFileName.length()-4);
    
    int maxContigSizeInFai = 0;

    //Verifying faiFile & reference size
    ifstream finFai(faiFileName.c_str());
    if(!finFai.is_open())
    {
        finFai.clear();
        finFai.close();
        cout << "ERROR: Could not find fai file: " << faiFileName << endl;
        cout << "Please run 'cora " << FAIGENERATE << " " << faiFileName.substr(0,faiFileName.length()-3) << "' or 'samtools faidx " << faiFileName.substr(0,faiFileName.length()-3) << "'";
        exit(99);
    }

    //In MultiChr Version only MaxChrSize is needed
    string faiLine;
    numChromInFai = 0;
    while(getline(finFai, faiLine))
    {
        stringstream faiLineSS(faiLine);

        string ctgName;
        int ctgSize;
        faiLineSS >> ctgName >> ctgSize;

        if(ctgSize > maxContigSizeInFai)
        {
            maxContigSizeInFai = ctgSize;
        }

        numChromInFai++;   
    }
    finFai.close();
   
    if(numChromInFai > 32000)
    {
        cout << "ERROR: Current verison of CORA only supports at most 32000 chromosomes within a reference file." << endl;
        cout << "Let denizy@mit.edu know if a version that supports more chromosomes is needed." << endl;
        exit(100);
    }
 
    return maxContigSizeInFai;
}

long long CheckFai_and_ReturnTotalGenomeSize(string faiFileName)
{
    long long totalGenomeSize = 0;

    //Verifying faiFile & reference size
    ifstream finFai(faiFileName.c_str());
    if(!finFai.is_open())
    {
        finFai.clear();
        finFai.close();
        cout << "ERROR: Could not find fai file: " << faiFileName << endl;
        cout << "Please run 'cora " << FAIGENERATE << " " << faiFileName.substr(0,faiFileName.length()-3) << "' or 'samtools faidx " << faiFileName.substr(0,faiFileName.length()-3) << "'";
        exit(98);
    }

    string faiLine;
    while(getline(finFai, faiLine))
    {
        stringstream faiLineSS(faiLine);

        string ctgName;
        int ctgSize;
        faiLineSS >> ctgName >> ctgSize;

        totalGenomeSize += ctgSize;
    }
    finFai.close();

    return totalGenomeSize;
}

void SetupAuxiliaryFiles()
{
    if(homTableConstructionFlag)
    {   
        string phyStr = SplitTwoSignal(numPhysicalSplits);
        
        //TODO(denizy) Check if homtable file paths exist and if they are writable 

        //Put auxiliary files in the temporary directory
        ofstream foutPhy(PhysicalSplitsFile.c_str());
        foutPhy << numPhysicalSplits << endl << phyStr << endl;
        foutPhy.clear();
        foutPhy.close();
        
        string parStr = SplitFourSignal(numParallelSplits);
        ofstream foutPar(ParallelSplitsFile.c_str());
        foutPar << numParallelSplits << endl << parStr << endl;
        foutPar.clear();
        foutPar.close();
    }
    else if(runMode[0] != '0') //Fastq split is only needed for the first stage
    {
        if(numFastqSplits == AUTO_FASTQ_SPLIT_MODE)
        {
            long long totalGenomeSize = CheckFai_and_ReturnTotalGenomeSize(refFile + ".fai");

            numFastqSplits = 1 + totalGenomeSize / AUTO_FASTQ_SPLIT_FUNCTION_INTERVAL;      
            if((mappingReportMode == "BEST" || mappingReportMode == "BEST_FAST"))
            {
                numFastqSplits = 1 + totalGenomeSize / AUTO_FASTQ_SPLIT_FUNCTION_INTERVAL_LOWMEM;   
            }

            if(numFastqSplits > MAX_AUTO_FASTQ_SPLIT)
            {
                numFastqSplits = MAX_AUTO_FASTQ_SPLIT;
            } 
        }

        if(numFastqSplits > 1)
        {
            string fastqSplitStr = SplitFourSignal(numFastqSplits);
            
            //Edit string here to make it fit fastq splitting -- This is because fastqSplit executable requires a particular input
            ofstream foutFastqSplit(collapseFragmentFile.c_str());
            foutFastqSplit << numFastqSplits << " " << 4 << endl;       

            stringstream fastqSplitStrSS(fastqSplitStr);

            int lineCount = 0;
            string line;
            while(getline(fastqSplitStrSS, line))
            {
                lineCount++;
                
                string lineRC(line);
                for(int i=0; i<(int)line.length(); i++)
                {
                    lineRC[i] = line[line.length() - i - 1];

                    if(lineRC[i] == 'A')
                    {
                        lineRC[i] = 'T';
                    }
                    else if(lineRC[i] == 'C')
                    {
                        lineRC[i] = 'G';
                    }
                    else if(lineRC[i] == 'G')
                    {
                        lineRC[i] = 'C';
                    }
                    else if(lineRC[i] == 'T')
                    {
                        lineRC[i] = 'A';
                    }
                }

                string stringToPrint = line + " " + lineRC + " ";

                if(lineCount == numFastqSplits)
                {
                    stringToPrint = "N " + stringToPrint;
                }

                //Replace all spaces with underscores + add N at the end
                for(int i=0; i<(int)stringToPrint.length(); i++)
                {
                    if(stringToPrint[i] == ' ')
                    {
                        stringToPrint[i] = '_';
                    }
                }   
                    
                foutFastqSplit << stringToPrint << endl; //Addition of the N string to the last split
                //cout << stringToPrint << endl;
            }
        }
    }
}

void VerifyRefFile(string givenRef)
{
    ifstream finRef(givenRef.c_str());
    if(!finRef.is_open())
    {
        cout << "ERROR: Could not find reference file: " << refFile << endl;
        cout << "Exiting.." << endl;
        exit(10);
    }

    string refLine, seqLine;
    getline(finRef, refLine);
    getline(finRef, seqLine);
    
    if(refLine[0] != '>')
    {
        cout << "ERROR: Reference file is corrupted (not multi-fasta format)..." << endl;
        cout << "Exiting.." << endl;
        exit(11);
    }

    finRef.clear();
        finRef.close();
}

void VerifyRefFileAndFai()
{
    //Verifying refFile & reference line length
    ifstream finRef(refFile.c_str());
    if(!finRef.is_open())
    {
        cout << "ERROR: Could not find reference file: " << refFile << endl;
        cout << "Exiting.." << endl;
        exit(10);
    }

    string refLine, seqLine;
    getline(finRef, refLine);
    getline(finRef, seqLine);
    
    if(refLine[0] != '>')
    {
        cout << "ERROR: Reference file is corrupted (not multi-fasta format)..." << endl;
        cout << "Exiting.." << endl;
        exit(11);
    }
        
    refLineLen = seqLine.length(); //assume all other lines would be the same

    cout << "Reference line length: " << refLineLen << "bp .." << endl;
    finRef.clear();
    finRef.close();

    //Check fai file and get the maximum chromosome size
    refFileFai = refFile + ".fai";
    maxChrSizeInRef = CheckFai_and_ReturnMaxChrSize(refFileFai);
        
    cout << "Max contig size in reference: " << maxChrSizeInRef << endl;
}

void CreateDirectoryIfDoesNotExist(string name)
{
    //Create temporary files folder if doesn't exist
    struct stat st = {0};
    if(stat(name.c_str(), &st) == -1)
    {
        mkdir(name.c_str(), 0700);
    }
}

//Set up auxiliary parameters 
void ConfigureRelatedParameters(string coraExecPath)
{
    //Create temporary files folder if doesn't exist
    CreateDirectoryIfDoesNotExist(temporaryDirectoryName);
    killSignalFile = temporaryDirectoryName + "/" + killSignalFile;
    collapseFragmentFile = temporaryDirectoryName + "/" + collapseFragmentFile;
    PhysicalSplitsFile = temporaryDirectoryName + "/" + PhysicalSplitsFile;
    ParallelSplitsFile = temporaryDirectoryName + "/" + ParallelSplitsFile;

    size_t slashPos = coraExecPath.rfind("/");
    string execPathPrefix;  
    if(slashPos == string::npos)
        execPathPrefix = "";
    else
        execPathPrefix = coraExecPath.substr(0,slashPos+1);  

    collapseExec = execPathPrefix + collapseExec;
    linkConstructExec = execPathPrefix + linkConstructExec;
    homTableSetupExec = execPathPrefix + homTableSetupExec;
    mappingInferenceExec = execPathPrefix + mappingInferenceExec;
    fastqSplitterExec = execPathPrefix + fastqSplitterExec;

    if(numChromInFai > NUMCHROMFORSHORT)
    {
        homTableSetupExec += "_L"; //This is the homTableSetup version that supports larger chromosomes
        mappingInferenceExec += "_L"; //This is the mappingInference version that supports larger chromosomes
        collapseExec += "_L"; //This is the collapse version that supports larger chromosomes
    }

    LPT2path = execPathPrefix + LPT2path;
    LPT4path = execPathPrefix + LPT4path;

    assert(killSignalFile != ""); //this is the file that an intentional kill will report in
    
    //Remove any existing kill signal files
    FILE* killIn = fopen(killSignalFile.c_str(), "r");
    if(killIn)
    {
        fclose(killIn);
        stringstream removeKillFile;
        removeKillFile << "rm " << killSignalFile << endl;
        ExecuteSystemCall(removeKillFile, VERIFY);  
    }

    completeReadLen = readLen;

    if(splitMode == "HALF")
    {
        readLen /= 2; //If full read length is odd, it will be (fullLen-1)/2    
    }
    else if(splitMode == "THREEWAY")
    {
        readLen /= 3; //If full read length is not divisible by three, k-mer len will be interger division
    }

    //Check whether config file is proper
    if(runMode.length() != 4)
    {
        cout << "ERROR: RunMode should have 4 stages defined" << endl;
        exit(5);
    }

    //This generates the splitting scheme for both physical and parallel splits
    SetupAuxiliaryFiles();

    if(homTableConstructionFlag)
    {
        //Check whether splits files are prepared (physical splits are due to repeated smaller runs for memory management -- parallel splits are simultaenous)
        ifstream finPhySplit(PhysicalSplitsFile.c_str());
        ifstream finParSplit(ParallelSplitsFile.c_str());
        if(!finPhySplit.is_open() || !finParSplit.is_open())
        {
            cout << "ERROR: AUX split files do not exist" << endl;
            cout << "Exiting.." << endl;
            exit(7);
        }   

        //There will be as many inexact homology table partitions as the number of parallel splits (they are merged during collapsing)
        finParSplit >> numInexactHomTablePartitions;
        
        finPhySplit.clear();
        finPhySplit.close();
        finParSplit.clear();
        finParSplit.close();

        //Check whether we are overwriting the hom table files
        if(homTableConstructionFlag)
        {
            ifstream finHomPerf(homTablePerfectMatchFile.c_str(), ios::binary);     
            if(finHomPerf.is_open())
            {
                cout << "ERROR: HomTable Perf file already exists: '" << homTablePerfectMatchFile.c_str() << "'";
                cout << "Exiting.." << endl;
                exit(8);
            }
            finHomPerf.close();
            ifstream finHomInexact(homTableInexactMatchFile.c_str(), ios::binary);  
            if(finHomInexact.is_open())
            {
                cout << "ERROR: HomTable Inexact file already exists: '" << homTableInexactMatchFile.c_str() << "'";
                cout << "Exiting.." << endl;
                exit(9);
            }
            finHomInexact.close();
        } //if flag is false, there is no reason to check files until search part of island
    }   

    if(runMode != "0000") //If there is any non-zero runMode option, still read through the fai
    {
        //Set read length -- read lengths should be uniform throughout the dataset
        ifstream finReadFile(fastqInputListFile.c_str());
        if(!finReadFile.is_open())
        {
            cout << "ERROR: Read fastq file could not be opened: " << fastqInputListFile << endl;
            cout << "Exiting.." << endl;
            exit(12);
        }
        else
        {
            string firstLine;
            getline(finReadFile, firstLine);
            numSamples = atoi(firstLine.c_str());
            
            string listLine;
            for(int i=0; i<numSamples; i++)
            {
                if(!getline(finReadFile, listLine))
                {
                    cout << "Mismatch between number of samples and number of lines in input list" << endl;
                    exit(45);
                }

                stringstream listLineSS(listLine);

                if(inputMode.substr(0,6) == "PAIRED")
                {
                    string fileName_left, fileName_right;
                    unsigned long long count;
                
                    if(listLineSS >> fileName_left >> fileName_right >> count)
                    {
                        ifstream left(fileName_left.c_str());
                        ifstream right(fileName_right.c_str());
                        if(!left.is_open() || !right.is_open())
                        {
                            cout << "One of the following files is missing: " << endl;
                            cout << fileName_left << endl;
                            cout << fileName_right << endl;
                            exit(46);
                        }       
                        
                        left.clear();
                        left.close();
                        right.clear();
                        right.close();  
                    
                        numTotalReads += count;
                    }
                    else
                    {
                        cout << "The input list line: " << listLine << endl;
                        cout << "It should contain [FIRST READ FILE] [SECOND READ FILE] [COUNT]" << endl;
                        exit(47);
                    }
                }
                else //if(inputMode.substr(0,6) == "SINGLE")
                {
                    string fileName_single;
                    unsigned long long count;
                    if(listLineSS >> fileName_single >> count)
                    {
                        ifstream single(fileName_single.c_str());
                        if(!single.is_open())
                        {
                            cout << "One of the following files is missing: " << endl;
                            cout << fileName_single << endl;
                            exit(46);
                        }       
                        
                        single.clear();
                        single.close();
                    
                        numTotalReads += count;
                    }
                    else
                    {
                        cout << "The input list line: " << listLine << endl;
                        cout << "It should contain [SINGLE END FILE] [COUNT]" << endl;
                        exit(47);
                    }
                }
            }
        }
        finReadFile.clear();
        finReadFile.close();

        CapOnTotalNumberOfReads = numTotalReads;
        
        //Configure identity related parameters
        identitySpaceUpperLimit = CapOnTotalNumberOfReads * 2 + 5; //Times four is for storing the directionality of compression and mate information (will grow in the future with splits) [There is 5 offset to allow the highest 5 idDigitLen strings to be special cases]

        if(inputMode.substr(0,6) == "PAIRED")
        {
            identitySpaceUpperLimit *= 2;
        }

        if(splitMode == "HALF")
        {
            identitySpaceUpperLimit *= 2;
        }
        else if(splitMode == "THREEWAY")
        {
            identitySpaceUpperLimit *= 3;
        }
        
        idDigitLen = 0;
        while(identitySpaceUpperLimit > 0)
        {
            idDigitLen++;
            identitySpaceUpperLimit /= IDENTITY_ALPHABET_SIZE;
        }

        idDigitLen = (unsigned char) max(3, (int) idDigitLen); //Minimum idDigitLen is 3 for division purposes of long read names (unlikely but can be a faulty corner case)
    }
}

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void PrintHomManual(string error)
{
    cout << "====================================================================================================" << endl;
    cout << "Usage: cora " << HOMINDEX << " [options] <Reference> <ExactHom> <InexactHom>" << endl ;
    cout << "====================================================================================================" << endl;
    
    cout << "Notes: Reference should be in fasta or multi-fasta format and indexed" << endl;
    cout << "       You can use 'cora " << FAIGENERATE << "' or 'samtools faidx' to index the reference" << endl;
    cout << "       ExactHom and InexactHom are file names to be used in the search pipeline of CORA" << endl << endl;

    cout << "Important Options: " << endl << endl;
    cout << "    -K     [INT 33-64] k-mer length [no default]" << endl;
    cout << "    -H     [INT 1-3]   Hamming-distance threshold per k-mer for inexact homology table. [2] " << endl << endl;

    cout << "Performance Options: " << endl << endl;
    cout << "    -p     [INT 1-10]  Number of physical file splits for computing the homology tables [10]" << endl;
    cout << "                       Higher values substantially reduce memory use, sacrificing some run time." << endl;
    cout << "    -t     [INT 1-24]  Number of parallel threads to use for computing the inexact homology table [8]" << endl << endl;

    if(error != "")
    {
        cout << "====================================================================================================" << endl;
        cout << "Input Error: "; printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", error.c_str());
        cout << "Please see manual above... " << endl << endl;
    }

    exit(0);
}

void PrintFaiGenerateManual(string error)
{
    cout << "====================================================================================================" << endl;
    cout << "Usage: cora " << FAIGENERATE << " <Reference>" << endl ;
    cout << "====================================================================================================" << endl << endl;
    cout << "Notes: This command is a replacement for 'samtools faidx' in case samtools is not installed" << endl;
    cout << "       Reference should be in fasta or multi-fasta format" << endl << endl;

    if(error != "")
    {
        cout << "====================================================================================================" << endl;
        cout << "Input Error: "; printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", error.c_str());
        cout << "Please see manual above... " << endl << endl;
    }

    exit(0);
}

void PrintReadFileGenManual(string error)
{
    cout << "====================================================================================================" << endl;
    cout << "Usage: To generate SINGLE-end read file list: " << endl << endl;
    cout << "        cora " << READFILEGEN << " [option] <FileName> -S <FASTQ-A> <FASTQ-B> <...>" << endl << endl;
    cout << "       To generate PAIRED-end read file list: " << endl << endl;
    cout << "        cora " << READFILEGEN << " [option] <FileName> -P <FASTQ-A1> <FASTQ-A2> <FASTQ-B1> <FASTQ-B2> <...>" << endl << endl;
    cout << "====================================================================================================" << endl << endl;
    cout << "Notes: <FileName> is the output read file list name for CORA, which stores read dataset info" << endl;
    cout << "       FASTQ files are input read datasets. Paired and single-end read datasets cannot be mixed" << endl << endl;
 
    cout << "Options   : --ReadComp  [STRING] Compression format of the input reads (GZIP or OFF) [default OFF]" << endl;    
    cout << "                                     OFF  -> input datasets are in plain FASTQ format" << endl;
    cout << "                                     GZIP -> input datasets are compressed with gzip" << endl << endl; 

    if(error != "")
    {
        cout << "====================================================================================================" << endl;
        cout << "Input Error: "; printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", error.c_str());
        cout << "Please see manual above... " << endl << endl;
    }

    exit(0);    
}

void PrintMapperIndexManual(string error)
{
    cout << "====================================================================================================" << endl;
    cout << "Usage: cora " << MAPPERINDEX << " [options] <Reference>" << endl ;
    cout << "====================================================================================================" << endl << endl;
    
    cout << "Important Options: " << endl << endl;
    cout << "    --Map        Off-the shelf mapper to be used for coarse mapping (stages 2-3) [default BWA]" << endl;
    cout << "                 current valid options: BWA, BWA_MEM, BOWTIE, BOWTIE_2, MRSFAST, MRSFAST_ULTRA" << endl;
    cout << "    --Exec       Name of the executable path for mapper [default bwa]" << endl;
    cout << "                 Can be full path (e.g. /home/folder/bwa) or just executable if installed (e.g. bwa)" << endl;
    cout << "    --Index      Full path for the index to be constructed if [default is same as <Reference>]" << endl << endl;
    
    cout << "Miscellaneous Options: " << endl << endl;
    cout << "    --opt        Additional options to be passed to the mapper for indexing" << endl;
    cout << "                 Warning: some of the non-default indexing options may not be compatible with CORA" << endl << endl;
    
    if(error != "")
    {
        cout << "====================================================================================================" << endl;
        cout << "Input Error: "; printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", error.c_str());
        cout << "Please see manual above... " << endl << endl;
    }
    
    exit(0);
}


void PrintSearchManual(string error)
{
    cout << "====================================================================================================" << endl;
    cout << "Usage: cora " << SEARCH << "[options] <Read_File_List> <Reference> <ExactHom> <InexactHom>" << endl;
    cout << "====================================================================================================" << endl << endl;
    
    cout << "Notes: <Reference> should be in fasta or multi-fasta format and indexed." << endl;
    cout << "       <ExactHom> and <InexactHom> are files generated in the " << HOMINDEX << " run" << endl;
    cout << "       <Read_File_List> can be generated using " << READFILEGEN << " command or manually." << endl << endl;

    cout << "Important Options: " << endl << endl;
    cout << "    -C            Determines which stages of CORA pipeline will be executed. [default 1111]" << endl;
    cout << "                  Other valid options are: 1000, 1100, 1110, 0100, 0110, 0111, 0010, 0011, 0001" << endl;
    cout << "                  These indicate 4 flags for running following 4 stages: " << endl;
    cout << "                      1) Compressing reads into unique k-mers." << endl;
    cout << "                      2) Coarse mapping k-mers to the reference using an off-the-shelf aligner." << endl;
    cout << "                      3) Converting coarse mapping results to read links." << endl;
    cout << "                      4) Traverse homology index for generating final mappings." << endl;

    cout << "    --Mode        Mapping mode: ALL, BEST, BEST_SENSITIVE, BEST_FAST, STRATUM, or UNIQUE [ALL]." << endl;

    cout << "    --Map         Off-the shelf mapper to be used for coarse mapping [default BWA]" << endl;
    cout << "                  Can be BWA, BWA_MEM, BOWTIE, BOWTIE_2, MRSFAST, MRSFAST_ULTRA, MANUAL" << endl;

    cout << "    --Exec        Name of the executable file to be run for coarse mapping" << endl;
    cout << "                  Either full path (e.g. /home/folder/bwa) or executable if installed (e.g. bwa)" << endl;

    cout << "    -R            Read input mode, SINGLE for single-end or PAIRED for paired-end reads [PAIRED]" << endl;

    cout << "    --MinI        [INT] for PAIRED (paired-end) mode, min insert length (|TLEN| in SAM) [150]" << endl;
    cout << "    --MaxI        [INT] for PAIRED (paired-end) mode, max insert length (|TLEN| in SAM) [650]" << endl;

    cout << "    -O            Output SAM file to print the final mapping output [default = CORA_Output.sam]" << endl;

    cout << "    -L            [INT] Read length for CORA mapping (per read end) [no default]. " << endl;
    
    cout << "    -K            K-mer compression mode, FULL, HALF or THREEWAY [default = HALF]." << endl;
    cout << "                  FULL, HALF and THREEWAY denote 1, 2, or 3 k-mers per read-end, respectively." << endl; 
    cout << "                  K-mer length in " << HOMINDEX << " stage should concordant with -K (see below): " << endl;
    cout << "                      FULL -> index k-mer should be same as read length" << endl;
    cout << "                      HALF -> k-mer length should be floor(read_length/2)" << endl;
    cout << "                      THREEWAY -> k-mer length should be floor(read_length/3)" << endl; 

    cout << "    --Metric      Distance metric to be used for mapping (HAMMING or EDIT) [default HAMMING]" << endl;
    cout << "                      HAMMING -> distance is substitution only" << endl;
    cout << "                      EDIT -> (Levenshtein) distance allows indels (requires --Map BWA or BOWTIE_2)" << endl;
    
    cout << "    -E            [INT] Distance threshold per k-mer for inexact homology table. (1 to 3) [default 2] " << endl;
    cout << "                  Use of Hamming or Edit distance is determined by --Metric option." << endl;
    cout << "                  Should be the same as -H used in " << HOMINDEX << " stage. " << endl;

    cout << "Performance Options: " << endl << endl;
    cout << "    --memoi       [INT] Activate memoization for k-mers appearing at least as many times as INT [20]." << endl;
    cout << "                  Must be larger than 1. Lower values save runtime during Stage 4 using more memory." << endl;
    
    cout << "    --fs          [INT or AUTO] Number of physical file splits for Input FASTQ file for compression" << endl;
    cout << "                  Valid parameters are either AUTO or an integer between 1 to 144 [default is AUTO] " << endl;
    cout << "                  Higher values reduce memory use during Stage 1 sacrificing some run time. " << endl;
    cout << "                  AUTO mode determines the number of FASTQ splits as: 1 + genome_size / 10M." << endl << endl;    

    cout << "    --cm          [STRING] Determines if the reference should be used for compressing reads or not." << endl;
    cout << "                  Options: WITHREF or NOREF [default WITHREF]" << endl << endl;

    cout << "Paralelization Options: " << endl << endl;
    cout << "    --coarseP     [INT] The number of parallel threads to be used for coarse mapping. [default is 1]" << endl;
    cout << "                  Capped at maximum number of threads allowed by coarse-mapper" << endl << endl;

    cout << "Miscellaneous Options: " << endl << endl;
    
    cout << "    --RG          [\"STRING\"] All read group data for the read datasets being mapped [default NONE]" << endl;
    cout << "                  This string will appear in the header and IDs will be attached to each mapping line" << endl;
    cout << "                  The whole string should be double quoted, the first identifier should be a unique ID" << endl;
    cout << "                  Multiple read datasets's RG data should be comma-separated in the order of Read_File_List" << endl;
    cout << "                  For tab-delimiting identifiers within the command line you can use Ctrl+V -> Ctrl+I" << endl; 
    cout << "                  The following is an example --RG value for 3 single-end or paired-end read datasets" << endl;
    cout << "                  \"ID:xx1\tCN:yya\tDS:za zb\tDT:ta,ID:xx2\tCN:yya,ID:xx3\tCN:yyb\tDS:zc zd\"" << endl << endl;

    cout << "    --TempDir     [STRING] The directory to be used for temp CORA files. [__temporary_CORA_files]." << endl;
    cout << "                  This enables running two CORA jobs in the same folder with different TempDir." << endl;

    cout << "    --ReadComp    [STRING] Compression format of the input reads (GZIP or OFF) [default OFF]" << endl;    
    cout << "                  OFF means that the input datasets are in plain FASTQ format" << endl << endl;

    cout << "    --MaxMapCount [INT] Maximum number of mapping to print in ALL mapping mode [default is infinite]" << endl << endl;

    if(error != "")
    {
        cout << "====================================================================================================" << endl;
        cout << "Input Error: "; printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", error.c_str());
        cout << "Please see manual above... " << endl << endl;
    }

    exit(0);
}

void ParseCommandLineArguments(int argc, char* argv[])
{
    cout << endl;

    if(argc == 1 || string(argv[1]) == "-v" || string(argv[1]) == "--help" || string(argv[1]) == "help")
    {
        cout << "====================================================================================================" << endl;
        printf("=======   " ANSI_COLOR_YELLOW "  CORA v" VERSION " by Deniz Yorukoglu (denizy@mit.edu, http://cora.csail.mit.edu/) " ANSI_COLOR_RESET  "   ========\n");
        cout << "====================================================================================================" << endl << endl;

        cout << "Usage: cora <command> [options]" << endl << endl;
        cout << "Command: " << HOMINDEX << "\t Create a homology table for CORA." << endl;
        cout << "         " << FAIGENERATE << "\t Generate .fai file for reference." << endl;
        cout << "         " << MAPPERINDEX << "\t Index the reference genome for coarse mapper." << endl;
        cout << "         " << READFILEGEN << "\t Automatically generate CORA's read file list input." << endl; 
        cout << "         " << SEARCH << "\t\t Run CORA's compressive read alignment pipeline." << endl;
        cout << endl << "To print the manual for each command, run './cora <command>' without any options." << endl << endl;
    
        exit(0);
    }

    string firstArg(argv[1]);
    if(firstArg == HOMINDEX)
    {
        // k-mer length [no default, INT] 
        // Num Physical Splits [default 1]
        // Num Parallel Splits [default 1]
        // Mismatches per k-mer [default 2] // There are some constraints regarding read-mer length

        homTableConstructionFlag = 1;
        runMode = "0000";
    
        if(argc <= 2 || string(argv[2]) == "--help" || string(argv[2]) == "help")
            PrintHomManual("");
        else
        {
            for(int i=2; i<argc;)
            {
                string argMark = string(argv[i]);

                if(argMark[0] == '-' && argc < i+4)
                {
                    PrintHomManual("Either missing parameter for argument " + argMark + " or last three arguments, <Reference> <ExactHom> <InexactHom>, are not set");              
                }           
    
                if(argMark == "-K")
                {
                    kmerLen = atoi(argv[i+1]);
                    if(kmerLen <= 32)
                        PrintHomManual("-K argument needs to be an integer larger than 32");
                    i+=2;
                }               
                else if(argMark == "-H")
                {
                    numMismatchesPerReadMer = atoi(argv[i+1]);
                    if(numMismatchesPerReadMer < 1 || numMismatchesPerReadMer > 3)
                        PrintHomManual("-H argument parameter should be 1 <= INT <= 3");
                    i+=2;
                }
                else if(argMark == "-p")
                {
                    numPhysicalSplits = atoi(argv[i+1]);        
                    if(numPhysicalSplits < 1 || numPhysicalSplits > 10)
                         PrintHomManual("-p argument parameter should be: 1 <= INT <= 10");
                    i+=2;
                }               
                else if(argMark == "-t")
                {
                    numParallelSplits = atoi(argv[i+1]);
                    if(numParallelSplits < 1 || numParallelSplits > 24)
                        PrintHomManual("-t argument parameter should be: 1 <= INT <= 24");
                    i+=2;
                }
                else
                {
                    if(argMark[0] == '-') //False command
                        PrintHomManual("Unknown command line argument: " + argMark);
                
                    //i is the reference
                    //i+1 is the exectHom
                    //i+2 is the inexactHom
                    if(argc != i+3)
                        PrintHomManual("<Reference> <ExactHom> <InexactHom> should be the last three arguments");
                    
                    refFile = argv[i];
                    homTablePerfectMatchFile = argv[i+1];
                    homTableInexactMatchFile = argv[i+2];   

                    cout << "refFile: " << refFile << "\thomTablePerfectMatchFile: " << homTablePerfectMatchFile << "\thomTableInexactMatchFile: " << homTableInexactMatchFile << endl;

                    if(homTablePerfectMatchFile == homTableInexactMatchFile)
                        PrintHomManual("<ExactHom> should be different than <InexactHom>");

                    VerifyRefFileAndFai();  

                    i+=3;
                }
            }
        
            if(refFile == "")
            {
                PrintHomManual("Missing <Reference> argument.");
            }   
    
            if(kmerLen == 0)
                PrintHomManual("-K argument is required."); 
        
            CreateDirectoryIfDoesNotExist(temporaryDirectoryName);
            logOut.open((string(temporaryDirectoryName + "/OUTPUT_run_log_file.") + HOMINDEX).c_str());
        }
    }
    else if(firstArg == FAIGENERATE)
    {
        if(argc != 3 || string(argv[2]) == "--help" || string(argv[2]) == "help")
        {
            PrintFaiGenerateManual("");
        }
        else
        {
            CreateDirectoryIfDoesNotExist(temporaryDirectoryName);
            logOut.open((string(temporaryDirectoryName + "/OUTPUT_run_log_file.") + FAIGENERATE).c_str());
            
            string refPath(argv[2]);
            string faiPath = refPath + ".fai";

            VerifyRefFile(refPath);
            
            stringstream faiGenerateCall;
            faiGenerateCall << "LC_ALL=C awk 'substr($1,1,1) == \">\" {if(curChr != 0) {printf \"%s\\t%ld\\t%.0f\\t%ld\\t%ld\\n\", curChrName, chrLen, chrOffset, lineBases, lineBases+1} curChrName = substr($1,2,length($1)-1); curChr+=1; prevLineIsChrName=1; totalOffset = totalOffset + length($0) + 1; chrOffset=totalOffset; chrLen=0} substr($1,1,1) != \">\" {if(prevLineIsChrName) {prevLineIsChrName=0; lineBases = length($1)} totalOffset+=length($0)+1; chrLen+=length($1)} END{printf \"%s\\t%ld\\t%0.f\\t%ld\\t%ld\\n\", curChrName, chrLen, chrOffset, lineBases, lineBases+1}' ";
        
            faiGenerateCall << refPath << " > " << faiPath;

            ExecuteSystemCall(faiGenerateCall, VERIFY);
        
            logOut.close();
            exit(0);
        }
    }
    else if(firstArg == READFILEGEN)
    {
        //This command should generate the read file list by taking in the list of read files the user wants to map and then count their number of reads (possibly by wc -l)
        
        if(argc <= 2 || string(argv[2]) == "--help" || string(argv[2]) == "help")
        {
            PrintReadFileGenManual("");
        }
        else
        {
            CreateDirectoryIfDoesNotExist(temporaryDirectoryName);
            logOut.open((string(temporaryDirectoryName + "/OUTPUT_run_log_file.") + READFILEGEN).c_str());

            if(argc < 5)
                PrintReadFileGenManual("Single-end (-S) or paired-end (-P) should be specified, and read FASTQ files should be listed");
        
            int curArgv = 2;
            string readCompressionMode = "OFF";
            if(string(argv[curArgv]) == "--ReadComp")
            {
                readCompressionMode = argv[curArgv+1];
                if(readCompressionMode != "OFF" && readCompressionMode != "GZIP")
                {
                    PrintSearchManual("Valid options for --ReadComp argument are OFF or GZIP");             
                }
                curArgv+=2;
            }
            
            string fileListName = argv[curArgv];
            ofstream foutFile(fileListName.c_str());
            curArgv++;          

            string readType = argv[curArgv];
            curArgv++;
            string tempFileLenFile = temporaryDirectoryName + "/__TEMP__Read_File_Len";

            if(readType == "-S")
            {
                int numDatasets = argc - curArgv; 

                foutFile << numDatasets << endl;

                for(int i=curArgv; i<argc; i++)
                {
                    stringstream fileLenCall;
                    if(readCompressionMode == "GZIP")
                        fileLenCall << "gunzip -c " << argv[i] << " | wc -l > " << tempFileLenFile;
                    else // OFF
                        fileLenCall << "wc -l " << argv[i] << " > " << tempFileLenFile;
                    ExecuteSystemCall(fileLenCall, VERIFY);
                
                    ifstream fileLenFile(tempFileLenFile.c_str());
                    int fileLength = 0;
                    fileLenFile >> fileLength;
            
                    if(fileLength == 0)
                    {
                        PrintReadFileGenManual(string("Read file \"") + argv[i] + "\" doesn't exist or is empty");
                    }

                    if(fileLength % 4 != 0)
                    {
                        PrintReadFileGenManual(string("Read file \"") + argv[i] + "\" isn't in FASTQ format or is possibly corrupted");
                    }

                    foutFile << argv[i] << "\t" << fileLength / 4 << endl;

                    fileLenFile.close();
                }
            }
            else if(readType == "-P")
            {
                if(argc % 2 == 1)
                {
                    PrintReadFileGenManual("In paired-end mode there should even number of read files listed (alternating between first mate and second mate)"); 
                }

                int numDatasets = (argc - curArgv)/2; 

                foutFile << numDatasets << endl; 

                for(int i=curArgv; i<argc; i+=2)
                {
                    stringstream fileLenCall;
                    if(readCompressionMode == "GZIP")
                        fileLenCall << "gunzip -c " << argv[i] << " | wc -l > " << tempFileLenFile;
                    else // OFF
                        fileLenCall << "wc -l " << argv[i] << " > " << tempFileLenFile;
                    ExecuteSystemCall(fileLenCall, VERIFY);
                    
                    ifstream fileLenFile(tempFileLenFile.c_str());
                    int fileLength = 0;
                    fileLenFile >> fileLength;
            
                    if(fileLength == 0)
                    {
                        PrintReadFileGenManual(string("Read file \"") + argv[i] + "\" doesn't exist or is empty");
                    }

                    if(fileLength % 4 != 0)
                    {
                        PrintReadFileGenManual(string("Read file \"") + argv[i] + "\" isn't in FASTQ format or is possibly corrupted");
                    }

                    foutFile << argv[i] << "\t" << argv[i+1] << "\t" << fileLength / 4 << endl;

                    fileLenFile.close();
                }               
            }
            else
            {
                PrintReadFileGenManual("The argument after <FileListName> should be either '-S' or '-P'; but is instead: " + readType);
            }

            foutFile.close();
        }

        logOut.close();
        exit(0);

    }
    else if(firstArg == MAPPERINDEX)
    {
        //Call reference indexing using the provided mapper 

        //Required parameters
        //Mapper type: [no default]] BWA or BOWTIE2 or MRSFAST or MRSFAST-ULTRA or MANUAL
        //Mapper executable: [no default]
        //Modified Indexing options in quotes '': Default empty
        //Reference file in fasta format: [no default]

        if(argc <= 2 || string(argv[2]) == "--help" || string(argv[2]) == "help")
        {
            PrintMapperIndexManual("");
        }
        else
        {
            string indexPath = "";
            string referencePath = "";
            string mapperType = "BWA";
            string mapperExecPath = "bwa";
            string addOptions = "";
            
            for(int i=2; i<argc; )
            {
                string argMark = string(argv[i]);
                
                if(argMark[0] == '-' && argc < i+2)
                {
                    PrintMapperIndexManual("Either missing parameter for argument " + argMark + " or <Reference> is not set");              
                }           
                                
                if(argMark == "--Map" )
                {
                    mapperType = argv[i+1];

                    if(mapperType == "BWA")
                        mapperExecPath = "bwa";
                    else if(mapperType == "BWA_MEM")
                        mapperExecPath = "bwa";
                    else if(mapperType == "BOWTIE")
                        mapperExecPath = "bowtie2";
                    else if(mapperType == "BOWTIE_2")
                        mapperExecPath = "bowtie";  
                    else if(mapperType == "MRSFAST")
                        mapperExecPath = "mrsfast";
                    else if(mapperType == "MRSFAST_ULTRA")
                        mapperExecPath = "mrsfast";
                    else
                    { 
                        PrintMapperIndexManual("Supported mappers for coarse mapping are: BWA, BWA_MEM, BOWTIE, BOWTIE_2, MRSFAST, MRSFAST_ULTRA, MANUAL\nIf the coarse mapping method of your choice is not available, please use MANUAL coarse mapping guide to manually pass coarse mapping commands to CORA.\nIf MANUAL mode fails or if you would like better integration of a mapper to CORA pipeline, please send an e-mail to denizy@mit.edu");  
                    }
                    i+=2;
                }
                else if(argMark == "--Exec")
                {
                    mapperExecPath = argv[i+1];
                    i+=2;
                }
                else if(argMark == "--Index")
                {
                    indexPath = argv[i+1];
                    i+=2;
                }   
                else if(argMark == "--opt")
                {
                    addOptions = argv[i+1];
                    i+=2;   
                }   
                else
                {
                    if(argMark[0] == '-') //False command
                        PrintMapperIndexManual("Unknown command line argument: " + argMark);
                    
                    if(argc != i+1)
                        PrintMapperIndexManual("<Reference> should be the last argument");
                
                    referencePath = argv[i];
                    i++;
                }
            }

            if(referencePath == "")
            {
                PrintMapperIndexManual("Missing <Reference> argument.");
            }
        
            CreateDirectoryIfDoesNotExist(temporaryDirectoryName);
            logOut.open((string(temporaryDirectoryName + "/OUTPUT_run_log_file.") + MAPPERINDEX).c_str());


            if(indexPath == "")
                indexPath = referencePath;

            if(mapperType == "BWA" || mapperType == "BWA_MEM")
            {
                stringstream bwaIndexCall;
                bwaIndexCall << mapperExecPath << " index " << addOptions << " -p " << indexPath << " " << referencePath;
                ExecuteSystemCall(bwaIndexCall, VERIFY);
            }   
            else if(mapperType == "MRSFAST" || mapperType == "MRSFAST_ULTRA")
            {
                stringstream mrsFastIndexCall;
                mrsFastIndexCall << mapperExecPath << " --index " << referencePath << " " << addOptions;
                ExecuteSystemCall(mrsFastIndexCall, VERIFY);
            }
            else if(mapperType == "BOWTIE" || mapperType == "BOWTIE_2")
            {
                stringstream bowtieIndexCall;
                bowtieIndexCall << mapperExecPath << "-build " << addOptions << " " << referencePath << " " << indexPath;
                ExecuteSystemCall(bowtieIndexCall, VERIFY);
            }
            else
            {
                assert(0);
            }   
        }

        exit(0);
    }
    else if(firstArg == SEARCH)
    {
        //Call remainder of cora pipeline

        //Required parameters
        // RunMode = default 1111, other valid options 1000, 1100, 1110, 0100, 0110, 0010, 0011, 0001
            // This indicates 4 flags for running 4 stages of CORA's Compressive Search
            // 1) Compressing reads into unique k-mers
            // 2) Coarse mapping k-mers to the reference
            // 3) Converting mapping results to links
            // 4) Traverse links and homology tables simultaneously for generating final mappings

        //Number of physical splits for FASTQ read compression: [default 1] (valid >1) Required only 
                
        //Min insert size: (default 150) Required only if Read dataset type is PAIRED and flag 4 is turned on
        //Max insert size: (default 650) Required only if Read dataset type is PAIRED and flag 4 is turned on
        
        //Mapping Mode: [default ALL] [Other valid options: BEST, STRATUM, UNIQUE (required only if 4th flag is on)
        //Memoziation Threshold: default 20 (valid >1) [required only if 4th flag is on]
        //Collapse Mode

        //Output file name [default CORA_Output.sam]  (required only if 4th flag is on)

        homTableConstructionFlag = 0;

        if(argc <= 2 || string(argv[2]) == "--help" || string(argv[2]) == "help")
        {
            PrintSearchManual("");
        }
        else
        {
            for(int i=2; i<argc; )
            {
                string argMark = string(argv[i]);
            

                if(argMark[0] == '-' && argc < i+5)
                {
                    PrintSearchManual("Either missing parameter for argument " + argMark + " or last four arguments, <Read_File_List> <Reference> <ExactHom> <InexactHom>, are not set");               
                }           

                if(argMark == "-K")
                {
                    splitMode = argv[i+1];
                    if(splitMode != "HALF" && splitMode != "FULL" && splitMode != "THREEWAY") //need parameter value
                        PrintSearchManual("Valid options for -K argument are: HALF or FULL");
                    i+=2;
                }               
                else if(argMark == "-L")
                {
                    readLen = atoi(argv[i+1]);
                    if(readLen < 32)
                        PrintSearchManual("-L argument requires a parameter: INT > 32");    
                    i+=2;   
                }
                else if(argMark == "-E")
                {
                    numMismatchesPerReadMer = atoi(argv[i+1]);
                    if(numMismatchesPerReadMer < 1 || numMismatchesPerReadMer > 3)
                        PrintSearchManual("-H argument requires a parameter: 1 <= INT <= 3");
                    i+=2;
                }
                else if(argMark == "-C")
                {
                    runMode = argv[i+1];        

                    if(runMode.length() != 4)
                        PrintSearchManual("-C argument requires 4 flags (set as 0 or 1)");
                    
                    for(int k=0; k<4; k++)
                        if(runMode[k] != '0' && runMode[k] != '1')
                            PrintSearchManual("-C argument requires all values to be set to 0 or 1 (0 indicates skipped stage, 1 indicates executed stage)");
                    
                    if(runMode == "0000")
                        PrintSearchManual("-C argument requires at least one stage's flag to be turned on");
                    
                    if(runMode[1] == '0')
                        if(runMode[2] + runMode[3] != 0)
                            if(runMode[0] == 1)
                                PrintSearchManual("Valid -C argument parameters are: 1000, 1100, 1110, 1111, 0100, 0110, 0111, 0010, 0011, 0001"); 
                    if(runMode[2] == '0')
                        if(runMode[0] + runMode[1] != 0)
                            if(runMode[4] != 0)
                                PrintSearchManual("Valid -C argument parameters are: 1000, 1100, 1110, 1111, 0100, 0110, 0111, 0010, 0011, 0001");
                    i+=2;
                }   
                else if(argMark == "--Mode")    
                {
                    mappingReportMode = argv[i+1];
                    if(mappingReportMode != "ALL" && mappingReportMode != "BEST" && mappingReportMode != "STRATUM" && mappingReportMode != "UNIQUE" && mappingReportMode != "BEST_SENSITIVE" && mappingReportMode != "BEST_FAST")
                        PrintSearchManual("Supported mapping report modes are: ALL, BEST, STRATUM, and UNIQUE\nPlease let us know about suggestions for additional useful mapping modes to CORA at denizy@mit.edu");
                    i+=2;
                }
                else if(argMark == "--Map")
                {
                    RepresentativeMappingMode = argv[i+1];
                    
                    if(RepresentativeMappingMode == "BWA")
                        RepresentativeMapperExecutable = "bwa";
                    else if(RepresentativeMappingMode == "BWA_MEM")
                        RepresentativeMapperExecutable = "bwa";
                    else if(RepresentativeMappingMode == "BOWTIE")
                        RepresentativeMapperExecutable = "bowtie2";
                    else if(RepresentativeMappingMode == "BOWTIE_2")
                        RepresentativeMapperExecutable = "bowtie";  
                    else if(RepresentativeMappingMode == "MRSFAST")
                        RepresentativeMapperExecutable = "mrsfast";
                    else if(RepresentativeMappingMode == "MRSFAST_ULTRA")
                        RepresentativeMapperExecutable = "mrsfast";
                    else if(RepresentativeMappingMode == "MANUAL")
                    {
                        //Executable specification is required
                    }
                    else
                    { 
                        PrintSearchManual("Supported mappers for coarse mapping are: BWA, BWA_MEM, BOWTIE, BOWTIE_2, MRSFAST, MRSFAST_ULTRA, MANUAL\nIf the coarse mapping method of your choice is not available, please use MANUAL coarse mapping guide to manually pass coarse mapping commands to CORA.\nIf MANUAL mode fails or if you would like better integration of a mapper to CORA pipeline, please send an e-mail to denizy@mit.edu");  
                    }
                    i+=2;
                }
                else if(argMark == "--Metric")
                {
                    mapperMetric = argv[i+1];
                    
                    if(mapperMetric != "HAMMING" && mapperMetric != "EDIT")
                    {
                        PrintSearchManual("Supported mapping distance metrics are HAMMING and EDIT. HAMMING distance is substitution only. EDIT distance can allow indels.");   
                    }
                    i+=2;
                }
                else if(argMark == "--Exec")
                {
                    RepresentativeMapperExecutable = argv[i+1];
                    i+=2;
                }
                else if(argMark == "-R")
                {
                    inputMode = argv[i+1];
                    if(inputMode != "SINGLE" && inputMode != "PAIRED")
                        PrintSearchManual("Valid options for -R argument is SINGLE or PAIRED");
                    i+=2;
                }       
                else if(argMark == "--MinI")
                {
                    insertSizeLowerThreshold = atoi(argv[i+1]);
                    if(insertSizeLowerThreshold <= 0)
                        PrintSearchManual("--MinI should be larger than 0");
                    i+=2;
                }
                else if(argMark == "--MaxI")
                {
                    insertSizeUpperThreshold = atoi(argv[i+1]);
                    if(insertSizeUpperThreshold <= 0)
                        PrintSearchManual("--MaxI should be larger than 0");
                    i+=2;
                }
                else if(argMark == "-O")
                {
                    mappingOutputFile = argv[i+1];
                    i+=2;   
                }
                else if(argMark == "--memoi")
                {
                    memoizationThreshold = atoi(argv[i+1]); 
                    if(memoizationThreshold <= 1)
                        PrintSearchManual("--memoi parameter should be larger than 1");
                    i+=2;
                }
                else if(argMark == "--fs")
                {
                    if(string(argv[i+1]) == "AUTO")
                    {
                        numFastqSplits = AUTO_FASTQ_SPLIT_MODE;
                    }
                    else
                    {
                        numFastqSplits = atoi(argv[i+1]);
                        if(numFastqSplits < 1 || numFastqSplits > 144)
                            PrintSearchManual("--fs parameter should be 1 <= INT <= 144 or AUTO");
                    }
                    i+=2;           
                }
                else if(argMark == "--cm")
                {
                    collapseMode = argv[i+1];
                    if(collapseMode != "WITHREF" && collapseMode != "NOREF")
                        PrintSearchManual("Valid options for --cm argument are WITHREF or NOREF");
                    i+=2;   
                }
                else if(argMark == "--coarseP")
                {
                    numCoarseMappingThreads = atoi(argv[i+1]);
                    if(numCoarseMappingThreads < 1)
                        PrintSearchManual("Valid options for --coarseP are positive integers");
                    i+=2;
                }
                else if(argMark == "--RG")
                {
                    inputReadGroupString = string(argv[i+1]);
                    i+=2;
                }
                else if(argMark == "--TempDir")
                {
                    temporaryDirectoryName = argv[i+1];
                    i+=2;
                }
                else if(argMark == "--ReadComp")
                {
                    readCompressionMode = argv[i+1];
                    if(readCompressionMode != "OFF" && readCompressionMode != "GZIP")
                        PrintSearchManual("Valid options for --ReadComp argument are OFF or GZIP");
                    i+=2;
                }
                else if(argMark == "--MaxMapCount")
                {
                    string field(argv[i+1]);
                    for(int k=0; k<(int) field.length(); k++)
                    {
                        if(field[k] < '0' || field[k] > '9')
                        {
                            PrintSearchManual("Valid options for --MaxMapCount is a positive integer");
                        }
                    }
                    
                    globalMapCountLimit = atoi(argv[i+1]);
                    if(globalMapCountLimit <= 0)
                    {
                        PrintSearchManual("Valid options for --MaxMapCount is a positive integer");
                    }
                    

                    i+=2;
                }
                else
                {
                    // i = <Read_File_List>      i+1 = <Reference>     i+2 = <ExactHom>     i+3 = <InexactHom>

                    if(argMark[0] == '-') //False command
                        PrintSearchManual("Unknown command line argument: " + argMark);
                    
                    if(argc != i+4)
                        PrintSearchManual("<Read_File_List> <Reference> <ExactHom> <InexactHom> should be the last four arguments");
                
                    fastqInputListFile = argv[i];   
                    refFile = argv[i+1];
                    homTablePerfectMatchFile = argv[i+2];
                    homTableInexactMatchFile = argv[i+3];   

                    //cout << "fastqInputListFile: " << fastqInputListFile << "\trefFile: " << refFile << "\thomTablePerfectMatchFile: " << homTablePerfectMatchFile << "\thomTableInexactMatchFile: " << homTableInexactMatchFile << endl;

                    if(homTablePerfectMatchFile == homTableInexactMatchFile)
                        PrintSearchManual("<ExactHom> should be different than <InexactHom>");

                    VerifyRefFileAndFai();  

                    i+=4;
                }   
            }
    
            if(mapperMetric == "EDIT" && RepresentativeMappingMode != "BWA" && RepresentativeMappingMode != "BOWTIE_2")
            {
                PrintSearchManual("Currently the only supported coarse mapper for EDIT distance metric is BWA or BOWTIE_2. Other mappers will be added in the upcoming updates.");
            } 

            if(fastqInputListFile == "")
            {
                PrintSearchManual("Missing <Read_File_List> argument.");
            }
    
            if(readLen == 0)
            {
                PrintSearchManual("-L argument is required."); 
            }
            
            if(RepresentativeMappingMode == "MANUAL" && RepresentativeMapperExecutable == "")
            {
                PrintSearchManual("If MANUAL mode is used for --Map argument, --Exec argument is required.");
            }   

            if(inputMode == "PAIRED")
            {
                if(insertSizeLowerThreshold > insertSizeUpperThreshold)
                {
                    PrintSearchManual("--MinI cannot be larger than --MaxI");
                }
            }   
        }
        
        CreateDirectoryIfDoesNotExist(temporaryDirectoryName);
        logOut.open((string(temporaryDirectoryName + "/OUTPUT_run_log_file.") + SEARCH).c_str());
    }
    else
    {
        cout << "The first argument should be '" << HOMINDEX << "', '" << FAIGENERATE << "', '" << MAPPERINDEX << "', '" << READFILEGEN << "' , or '" << SEARCH << "'. For help, run \"cora help\"" << endl << endl;
        exit(98);
    }
}

int main(int argc, char* argv[])
{
    double beginTime = getTime();

    ParseCommandLineArguments(argc, argv);
    cout << "Command line arguments are parsed successfully..." << endl;

    ConfigureRelatedParameters(string(argv[0]));
    cout << "Related parameters are configured..." << endl;

    //Outside of runMode, with hom table flag is on, it will be construced and compacted
    if(homTableConstructionFlag == 1)
    {
        cout << "Running Step H: Constructing Homology Table" << endl;

        stringstream homSetupCall;
        homSetupCall << homTableSetupExec << " " << refFile << " " << kmerLen << " " << PhysicalSplitsFile << " " << homTablePerfectMatchFile << " " << homTableInexactMatchFile << " " << ParallelSplitsFile << " BOTH " << numMismatchesPerReadMer << " FULL";
        ExecuteSystemCall(homSetupCall, VERIFY);

        //Delete large temporary files generated during the construction of Final homTables (Be careful that, this should be deactivated for easier debugging)
        string homTablePerfectMatchFile_beforeCompacting = homTablePerfectMatchFile + "_preCompact_toDelete";
        string homTableInexactMatchFilePrefix_beforeCompacting = homTableInexactMatchFile + "_preCompact_toDelete";
        
        stringstream deleteTempHomTableFilesCall;
        deleteTempHomTableFilesCall << "rm " << homTablePerfectMatchFile_beforeCompacting << " " << homTableInexactMatchFilePrefix_beforeCompacting << "_t*";
        ExecuteSystemCall(deleteTempHomTableFilesCall, VERIFY);

        if(runMode == "0000")
        {
            exit(0);
        }
    }

    //FileNames defined in runMode.at(0) -- or Step 1
    string fragmentedFilePrefix = temporaryDirectoryName + "/__TEMP__fragmentedInput";
    string collapsedFastqFile_longReadsList = temporaryDirectoryName + "/Step1outA_CollapseOutput_DifficultReadNames";
    string collapsedFastqFile_longReadsSeperated = temporaryDirectoryName + "/Step1outB_CollapseOutput_CompactKmersToCoarseMap";
    string perfMapsFile = temporaryDirectoryName + "/Step1outC_CollapseOutput_PerfectKmerMatches";
    //FileNames defined in runMode.at(1) -- or Step 2
    string msMapFile = temporaryDirectoryName + "/Step2outA_CoarseMappingResult";
    //FileNames defined in runMode.at(2) -- or Step 3
    string mergedReadEditLinksFile = temporaryDirectoryName + "/Step3outA_KmerEditLinks";
    string pibLinks = temporaryDirectoryName + "/Step3outB_MergedKmerlinks"; //TODO(denizy) Eliminate the merging need, just read from two separate files

    //This stage runs collapsing, depending on whether reference is used for collapsing or fragmentation is done different runs will be made
    if(runMode[0] >= '1')
    {
        if(readCompressionMode == "GZIP")
        {
            //Unzip input files
            ifstream finReadFile(fastqInputListFile.c_str());
            string uncompressedFastqInputListFile = temporaryDirectoryName + "/__UncompressedReadFileList"; 
            ofstream foutModifiedReadFile(uncompressedFastqInputListFile.c_str());

            string line;
            getline(finReadFile,line); //first line is not used and printed without modification
            foutModifiedReadFile << line << endl;

            for(int lineCounter = 1; getline(finReadFile,line); lineCounter++)
            { 
                //Each field apart from the last is an input file in gzipped format
                stringstream lineSS(line);
                vector<string> fields;
                string field;
                while(lineSS >> field)
                {
                    fields.push_back(field);
                }
            
                for(int fieldNo = 0; fieldNo < (int) fields.size(); fieldNo++)
                {
                    if(fieldNo < (int) fields.size() - 1)
                    {
                        stringstream unzippedFileName;
                        unzippedFileName << temporaryDirectoryName << "/__uncompressedReadDataset" << lineCounter << ".field" << fieldNo; 
                    
                        //call gunzip on file               
                        stringstream gunzipCall;
                        gunzipCall << "gunzip -c " << fields[fieldNo] << " > " << unzippedFileName.str();
                        ExecuteSystemCall(gunzipCall, VERIFY);      
    
                        foutModifiedReadFile << unzippedFileName.str() << "\t";
                    }
                    else
                    {
                        foutModifiedReadFile << fields[fieldNo] << endl;
                    }
                }
            }
            finReadFile.close();
            foutModifiedReadFile.close();
        
            //fastqInputListFile is updated with the modified file
            fastqInputListFile = uncompressedFastqInputListFile;
        }

        cout << "Running Step 0: Collapsing Read Sequences" << endl;

        ifstream finFastqList(fastqInputListFile.c_str());

        int numLRNs = 0; //This keeps track of how many long read names have been printed (in order to keep track of which id they have across frag runs)

        //[This tool is updated by Aug 4 2013 version, now it automatically does the compact renaming as well as long read name separation 
        //Perform read collapsing and perfMap generation and mapping (with a workaround for long read names)
        if(numFastqSplits == 1)
        {
            stringstream CollapseCall;
            CollapseCall << collapseExec << " " << fastqInputListFile << " " << (int) idDigitLen << " " << collapsedFastqFile_longReadsSeperated << " ";

            if(collapseMode == "NOREF")
            {
                CollapseCall << "NOREF";
            }
            else if(collapseMode == "WITHREF")
            {
                CollapseCall << refFile;
            }

            CollapseCall << " " << refLineLen << " " << completeReadLen << " " << perfMapsFile << " " << maxChrSizeInRef << " " << collapsedFastqFile_longReadsList << " " << inputMode << " " << splitMode << " " << RepresentativeMappingMode;
            CollapseCall << " NONE 0 " << killSignalFile << " " << numLRNs;
        
            ExecuteSystemCall(CollapseCall, IGNORE); //ignore the kill signal since it is for faster deallocation
        }
        else
        {
            //Here call fastq splitter tool
            stringstream fastqSplitCall;
            fastqSplitCall << fastqSplitterExec << " " << collapseFragmentFile << " " << completeReadLen << " " << fastqInputListFile << " " << fragmentedFilePrefix << " ";

            //Determine whether splitting mode is SINGLE, SPLIT_SINGLE, PAIRED, SPLIT_PAIRED
            if(inputMode == "PAIRED")
            {
                if(splitMode == "HALF")
                {
                    fastqSplitCall << "SPLIT_PAIRED";
                }
                else if(splitMode == "THREEWAY")
                {
                    fastqSplitCall << "THREEWAY_PAIRED";
                }
                else
                {
                    fastqSplitCall << "PAIRED";
                }
            }
            else
            {
                if(splitMode == "HALF")
                {
                    fastqSplitCall << "SPLIT_SINGLE";                   
                }
                else if(splitMode == "THREEWAY")
                {
                    fastqSplitCall << "THREEWAY_SINGLE";
                }
                else
                {
                    fastqSplitCall << "SINGLE";
                }
            }
            ExecuteSystemCall(fastqSplitCall, VERIFY);

            //Here call collapseExec once for each split
            unsigned short numFrags = 0;
            ifstream fragFile(collapseFragmentFile.c_str());
            int numLines, fragSignalSize;
            fragFile >> numLines >> fragSignalSize;
            
            string tempWCLfileName = temporaryDirectoryName + "/__tempWCLfile";         

            string fragLine;
            while(fragFile >> fragLine)
            {
                stringstream fragmentedPerfMapsFileSS, fragmentedLongReadsListSS, fragmentedLongReadsSeparatedSS, fragmentedInputReadsFileSS;
                fragmentedPerfMapsFileSS << perfMapsFile << ".split" << numFrags;
                fragmentedLongReadsListSS << collapsedFastqFile_longReadsList << ".split" << numFrags;
                fragmentedLongReadsSeparatedSS << collapsedFastqFile_longReadsSeperated << ".split" << numFrags;
                fragmentedInputReadsFileSS << fragmentedFilePrefix << ".frag" << numFrags;
                
                stringstream fragmentedCollapseCall;
                fragmentedCollapseCall << collapseExec << " " << fragmentedInputReadsFileSS.str() << " " << (int) idDigitLen << " " << fragmentedLongReadsSeparatedSS.str() << " ";

                if(collapseMode == "NOREF")
                {
                    fragmentedCollapseCall << "NOREF";
                }
                else if(collapseMode == "WITHREF")
                {
                    fragmentedCollapseCall << refFile;
                }
                
                fragmentedCollapseCall << " " << refLineLen << " " << completeReadLen << " " << fragmentedPerfMapsFileSS.str() << " " << maxChrSizeInRef << " " << fragmentedLongReadsListSS.str() << " " << inputMode << " " << splitMode << " " << RepresentativeMappingMode << " " << fragLine << " " << fragSignalSize << " " << killSignalFile << " " << numLRNs;

                ExecuteSystemCall(fragmentedCollapseCall, IGNORE);
                
                stringstream wclCall;
                wclCall << "LC_ALL=C wc -l " << fragmentedLongReadsListSS.str() << " > " << tempWCLfileName;
                ExecuteSystemCall(wclCall, VERIFY);

                ifstream finWCL(tempWCLfileName.c_str());
                int curNumLRNs = -1;
                finWCL >> curNumLRNs;

                assert(curNumLRNs != -1);

                numLRNs += curNumLRNs; 

                finWCL.clear();
                finWCL.close();
                numFrags++;
            }

            //Cat-merge fragmented files
            stringstream catLongsTogether;
            catLongsTogether << "LC_ALL=C cat ";
            for(int i=0; i<numFrags; i++)
            {
                catLongsTogether << " " << collapsedFastqFile_longReadsList << ".split" << i;
            }
            catLongsTogether << " > " << collapsedFastqFile_longReadsList;          
            ExecuteSystemCall(catLongsTogether, VERIFY);

            stringstream catPerfsTogether;
            catPerfsTogether << "LC_ALL=C cat ";
            for(int i=0; i<numFrags; i++)
            {
                catPerfsTogether << " " << perfMapsFile << ".split" << i;
            }
            catPerfsTogether << " > " << perfMapsFile;
            ExecuteSystemCall(catPerfsTogether, VERIFY);

            stringstream catInexactsTogether;
            catInexactsTogether << "LC_ALL=C cat ";
            for(int i=0; i<numFrags; i++)
            {
                catInexactsTogether << " " << collapsedFastqFile_longReadsSeperated << ".split" << i; 
            }
            catInexactsTogether << " > " << collapsedFastqFile_longReadsSeperated;
            ExecuteSystemCall(catInexactsTogether, VERIFY);

            //TODO(denizy) Provide option clean up fragmented files
        }
    }

    int subErrorThreshold = numMismatchesPerReadMer;
    
    //coarse (representative) mapping for links table construction
    if(runMode[1] >= '1')
    {
        //[TODO: Add auto indexing for all files, but do it before collapsing, so that preprocessing doesn't get confused
        double representativeMapBeginTime = getTime();

        if(RepresentativeMappingMode == "MANUAL")
        {
            while(true)//RepresentativeMapperExecutable.find(";") != string::npos)
            {
                cout << "Rep: " << RepresentativeMapperExecutable << endl;
                
                size_t semiColonPos = RepresentativeMapperExecutable.find(";");

                string curSentence = RepresentativeMapperExecutable.substr(0, semiColonPos);
            
                while(curSentence.find("#INPUT#") != string::npos)
                {
                    int location = curSentence.find("#INPUT#");

                    string newCurSentence = curSentence.substr(0,location) + collapsedFastqFile_longReadsSeperated + curSentence.substr(location + 7, curSentence.length() - location - 7); 

                    curSentence = newCurSentence;
                }

                while(curSentence.find("#OUTPUT#") != string::npos)
                {
                    int location = curSentence.find("#OUTPUT#");

                    string newCurSentence = curSentence.substr(0,location) + msMapFile + curSentence.substr(location + 8, curSentence.length() - location - 8); 

                    curSentence = newCurSentence;   
                }
            
                while(curSentence.find("#REF#") != string::npos)
                {
                    int location = curSentence.find("#REF#");

                    string newCurSentence = curSentence.substr(0,location) + refFile + curSentence.substr(location + 5, curSentence.length() - location - 5); 

                    curSentence = newCurSentence;   
                }

                stringstream RepMapCall;
                RepMapCall << curSentence;

                ExecuteSystemCall(RepMapCall, VERIFY);

                if(semiColonPos == string::npos)
                {
                    break;
                }
                else
                {
                    RepresentativeMapperExecutable = RepresentativeMapperExecutable.substr(semiColonPos+1, RepresentativeMapperExecutable.length() - semiColonPos - 1);
                }
            }

            //some tools would add unnecessary extensions to the mapped file name, so check for them here and ask the user to give better commands for renaming things, such as a mv
            ifstream fin_check(msMapFile.c_str());
        
            if(!fin_check.is_open())
            {
                cout << "ERROR: the coarse mapping file is not properly generated or named" << endl;
                cout << "If the mapping file exists but need renaming (due to a file extension that mapper arbitrarily adds), please add another manual command sentence after ',', with the proper mv command, such as: mv #OUTPUT#.sam #OUTPUT#" << endl; 
                exit(123);
            }
            
            fin_check.clear();
            fin_check.close();
        }
        else if(RepresentativeMappingMode == "MRSFAST_ULTRA") //This is for the mrsfast-ULTRA 3.1
        {
            stringstream RepMapCall; //currently best map format --> normally it should either be any rep or iterative -e 1 -n 1 & -e 2 -n 1 calls
            RepMapCall << RepresentativeMapperExecutable << " --search " << refFile << " --seq " << collapsedFastqFile_longReadsSeperated << " --threads " << numCoarseMappingThreads << " --best " << subErrorThreshold << " -e " << subErrorThreshold << " -o " << msMapFile;
            ExecuteSystemCall(RepMapCall, VERIFY);
            //SAMPLE call: ../../MRSFAST_3_new/mrsfast --search ../REF/chr21.fa --seq Step0outA_CollapsedReads.fastq.OCform.longReadSeperated --threads 16 --mem 40 --best 2 -o Step0outC_CollapsedReads.fastq.OCform.longReadSeperated.Best2Chr21.sam

            stringstream renameCall; //this is for fixing extra ".sam" addition at the end
            renameCall << "mv " << msMapFile << ".sam " << msMapFile;
            ExecuteSystemCall(renameCall, VERIFY);
        }
        else if(RepresentativeMappingMode == "BOWTIE")
        {
            stringstream BowtieRepMapCall;
            BowtieRepMapCall << RepresentativeMapperExecutable << " --quiet -S -v " << subErrorThreshold << " --best -k 1 -f --sam-nohead " << refFile << " " << collapsedFastqFile_longReadsSeperated << " " << msMapFile;
            ExecuteSystemCall(BowtieRepMapCall, VERIFY);
            //SAMPLE call: ../../BOWTIE-TOOL/bowtie-0.12.8/bowtie --quiet -S -v 2 --best -k 1 -q --sam-nohead ../REF_BOWTIE/chr20.fa QUERY_Mapping/IBS6_QUERY/IBS6.samples.50bp.renamedAll
        }
        else if(RepresentativeMappingMode == "BOWTIE_2")
        {
            
            //double lvalue = -((1.5 + (subErrorThreshold + 6))  / double(readLen)); //1.5 is just for rounding safety  
            int maxPenalty = subErrorThreshold * 6;
            stringstream bowtie2call;

            
            if(mapperMetric == "HAMMING")
            {
                bowtie2call << RepresentativeMapperExecutable << " --end-to-end --ignore-quals -p " << numCoarseMappingThreads << " --rdg 1000,1000 --rfg 1000,1000 --np 6 --score-min L," << -maxPenalty << "," << 0 << " -x " << refFile << " -U " << collapsedFastqFile_longReadsSeperated << " -S " << msMapFile;
            }
            else //EDIT
            {
                bowtie2call << RepresentativeMapperExecutable << " --end-to-end --ignore-quals -p " << numCoarseMappingThreads << " --rdg 0,6 --rfg 0,6 --np 6 --score-min L," << -maxPenalty << "," << 0  << " -x " << refFile << " -U " << collapsedFastqFile_longReadsSeperated << " -S " << msMapFile;
            }

            ExecuteSystemCall(bowtie2call, VERIFY);
            //SAMPLE call for readLen=54  --> ../../BOWTIE2-TOOL/bowtie2-2.1.0/bowtie2 -q --end-to-end --very-sensitive --ignore-quals --rdg 100,100 --rfg 100,100 --np 6 --score-min L,-0.25,-0.25 -x ../REF_BOWTIE_2/chr20.fa -U Step0outA_ColReads.fq.OC.LRNSep.fauxq -S temp.moreParams3.best.sam
        }
        else if(RepresentativeMappingMode == "MRSFAST") //This is the best-mapping version of original mrsFAST (earlier than 3.0)
        {
            if(subErrorThreshold == 1) //If error threshold is 1, then no need for iterative running
            {
                stringstream SingleMrsfastCall;
                SingleMrsfastCall << RepresentativeMapperExecutable << " --search " << refFile << " --seq " << collapsedFastqFile_longReadsSeperated << " -e 1 -n 1 -o " << msMapFile;
            }
            else
            {
                //Need to do iterative alignment, either E1+E2 or E1+E2+E3

                string msMapFile_E1 =  msMapFile + "_E1";

                stringstream RepMapCall_E1;
                RepMapCall_E1 << RepresentativeMapperExecutable << " --search " << refFile << " --seq " << collapsedFastqFile_longReadsSeperated << " -e 1 -n 1 -o " << msMapFile_E1;
                ExecuteSystemCall(RepMapCall_E1, VERIFY);
        
                string E1_nohit_file = msMapFile_E1 + ".nohit";
                string msMapFile_E2 =  msMapFile + "_E2";

                stringstream RepMapCall_E2;
                RepMapCall_E2 << RepresentativeMapperExecutable << " --search " << refFile << " --seq " << E1_nohit_file << " -e 2 -n 1 -o " << msMapFile_E2;
                ExecuteSystemCall(RepMapCall_E2, VERIFY);
            
                string msMapFile_E3;
    
                if(subErrorThreshold > 2)
                {
                    string E2_nohit_file = msMapFile_E2 + ".nohit";
                    msMapFile_E3 =  msMapFile + "_E3";

                    stringstream RepMapCall_E3;
                    RepMapCall_E3 << RepresentativeMapperExecutable << " --search " << refFile << " --seq " << E2_nohit_file << " -e 3 -n 1 -o " << msMapFile_E3;
                    ExecuteSystemCall(RepMapCall_E3, VERIFY);
                }

                if(subErrorThreshold == 2)
                {
                    stringstream RepMapMergeCall;
                    RepMapMergeCall << "cat " << msMapFile_E1 << " " << msMapFile_E2 << " > " << msMapFile;
                    ExecuteSystemCall(RepMapMergeCall, VERIFY);
                }
                else if(subErrorThreshold == 3)
                {
                    stringstream RepMapMergeCall;
                    RepMapMergeCall << "cat " << msMapFile_E1 << " " << msMapFile_E2 << " " << msMapFile_E3 << " > " << msMapFile;
                    ExecuteSystemCall(RepMapMergeCall, VERIFY);
                }
                else
                {
                    cout << "ERROR: Currently CORA doesn't support error threshold higher than 3 per readMer. Your input was: " << subErrorThreshold << endl;
                    exit(99);
                }
            }
        }
        else if(RepresentativeMappingMode == "BWA")
        {
            string sai_file = msMapFile + "sai"; 
            
            if(mapperMetric == "HAMMING")
            {
                stringstream bwaRepMapCall;
                bwaRepMapCall << RepresentativeMapperExecutable << " aln -t " << numCoarseMappingThreads << " -n " << subErrorThreshold << " -o 0 " << refFile << " " << collapsedFastqFile_longReadsSeperated << " > " << sai_file;
                ExecuteSystemCall(bwaRepMapCall, VERIFY);
                //SAMPLE call: ../../../BWA-TOOL/bwa-0.7.5a/bwa aln -n 2 -o 0 ../../REF_BWA/chr20.bwa-0.7.5a.fa Step0outA_ColReads.fq.OC.LRNSep > Step0outA_ColReads.fq.OC.LRNSep.n2o0.sai

                stringstream samseCall;
                samseCall << RepresentativeMapperExecutable << " samse " << refFile << " " << sai_file << " " << collapsedFastqFile_longReadsSeperated << " > " << msMapFile; 
                ExecuteSystemCall(samseCall, VERIFY);
                //SAMPLE call: ../../../BWA-TOOL/bwa-0.7.5a/bwa samse ../../REF_BWA/chr20.bwa-0.7.5a.fa Step0outA_ColReads.fq.OC.LRNSep.n2o0.sai Step0outA_ColReads.fq.OC.LRNSep > Step0outA_ColReads.fq.OC.LRNSep.n2o0.sam
            }
            else
            {
                stringstream bwaRepMapCall;
                bwaRepMapCall << RepresentativeMapperExecutable << " aln -t " << numCoarseMappingThreads << " -n " << subErrorThreshold << " -o " << subErrorThreshold << " -d 0 -i 0 " << refFile << " " << collapsedFastqFile_longReadsSeperated << " > " << sai_file;
                ExecuteSystemCall(bwaRepMapCall, VERIFY);
                //SAMPLE call: ../../../BWA-TOOL/bwa-0.7.5a/bwa aln -n 2 -o 0 ../../REF_BWA/chr20.bwa-0.7.5a.fa Step0outA_ColReads.fq.OC.LRNSep > Step0outA_ColReads.fq.OC.LRNSep.n2o0.sai

                stringstream samseCall;
                samseCall << RepresentativeMapperExecutable << " samse " << refFile << " " << sai_file << " " << collapsedFastqFile_longReadsSeperated << " > " << msMapFile; 
                ExecuteSystemCall(samseCall, VERIFY);
                //SAMPLE call: ../../../BWA-TOOL/bwa-0.7.5a/bwa samse ../../REF_BWA/chr20.bwa-0.7.5a.fa Step0outA_ColReads.fq.OC.LRNSep.n2o0.sai Step0outA_ColReads.fq.OC.LRNSep > Step0outA_ColReads.fq.OC.LRNSep.n2o0.sam
            }
        }
        else if(RepresentativeMappingMode == "BWA_MEM")
        {
            stringstream bwaMemCall;
            int tValue = readLen - 2 * subErrorThreshold; //this is for two errors
            bwaMemCall << RepresentativeMapperExecutable << " mem -t " << numCoarseMappingThreads << " -A 1 -B 1 -O 1000 -E 1000 -L 1000 -U 1 -T " << tValue << " " << refFile << " " << collapsedFastqFile_longReadsSeperated << " > " << msMapFile;
            ExecuteSystemCall(bwaMemCall, VERIFY);
            //SAMPLE call for 54bp: /mnt/work/denizy/BWA-TOOL/bwa-0.7.5a/bwa mem -A 1 -B 6 -O 1000 -E 1000 -L 1000 -U 6 -T 40 /mnt/work/denizy/ISLAND_clone_folder/REF_BWA/chr20.bwa-0.7.5a.fa Step0outA_ColReads.fq.OC.LRNSep > Step0outA_ColReads.fq.OC.LRNSep.BWAMEMmap
        }
        else
        {
            cout << "There is no representative mapping mode named '" << RepresentativeMappingMode << endl;
            exit(6);
        }
        
        double representativeMapEndTime = getTime();
        cout << "RepMapTime: " << representativeMapEndTime - representativeMapBeginTime << endl;
    }

    //Construct links from sam outputs of coarse mapping
    if(runMode[2] >= '1')
    {
        if(RepresentativeMappingMode == "MRSFAST_ULTRA_NO_MEM_SPEC")
        {
            RepresentativeMappingMode = "MRSFAST_ULTRA";
        }
        
        //Repository Construction -- (Options are NOTER and NOREF for simpler variants)
        stringstream linkConstructCall;
        linkConstructCall << linkConstructExec << " " << refFile << " " << msMapFile << " NOREF NOREF ";
        linkConstructCall << collapsedFastqFile_longReadsList << " " <<  mergedReadEditLinksFile << " " << subErrorThreshold << " " << (int) idDigitLen << " " << RepresentativeMappingMode << " " << readLen << " " << refLineLen << " " << temporaryDirectoryName;
        ExecuteSystemCall(linkConstructCall, VERIFY);
        //This is changed ../ISLAND_Code_Rep/1_RepoLinks/linkConstruct_handleLongReadNames ../REF/chr21.fa Step0outC_CollapsedReads.fastq.OCform.longReadSeperated.Best2Chr21.sam NOREF NOREF NOTER Step1_MergedReference_OnlyA.fa Step1_ReadLinks_OnlyA Step1_ContigOffsetFile_OnlyA Step0outA_CollapsedReads.fastq.OCform.longReadNamesList
    }

    double constructTimeEnd = getTime();
    logOut << "Construction Time: " << constructTimeEnd - beginTime << endl;

    //Mapping inference stage
    if(runMode[3] >= '1')
    {
        stringstream mergePerfAndInexactLinksCall;
        mergePerfAndInexactLinksCall << "cat " << perfMapsFile << " " << mergedReadEditLinksFile << " > " << pibLinks;
        ExecuteSystemCall(mergePerfAndInexactLinksCall, VERIFY);

        stringstream mapAllInferCall;

        int meditSizeNeeded = numMismatchesPerReadMer;
        if(splitMode == "THREEWAY")
        {
            meditSizeNeeded *= 3;
        }
        else if(splitMode == "HALF")
        {
            meditSizeNeeded *= 2; 
        }


        if(meditSizeNeeded <= 4)
        {
            mapAllInferCall << mappingInferenceExec << " ";
        }
        else if(meditSizeNeeded <= 6)
        {
            mapAllInferCall << mappingInferenceExec << "_E3 ";
        }
        else if(meditSizeNeeded <= 8)
        {
            mapAllInferCall << mappingInferenceExec << "_E4 ";
        }
        else if(meditSizeNeeded <= 10)
        {
            mapAllInferCall << mappingInferenceExec << "_E5 ";
        }
        else if(meditSizeNeeded <= 12)
        {
            mapAllInferCall << mappingInferenceExec << "_E6 ";
        }
        else
        {
            cout << "Current version only supports up to 12 errors per read end" << endl;
            exit(8);
        }

        mapAllInferCall << homTablePerfectMatchFile << " " << homTableInexactMatchFile << " " << completeReadLen << " " << refFile << " " << refLineLen << " " << pibLinks << " ";

        if(mapperMetric == "EDIT")
        {
            mapAllInferCall << "LEVENSHTEIN ";
        }
        else
        {
            mapAllInferCall << "HAMMING ";
        }
    
        mapAllInferCall << mappingReportMode << " " << fastqInputListFile << " " << globalMapCountLimit << " " << mappingOutputFile << " UNCOLLAPSED";

        if(inputReadGroupString != "")
        {
            mapAllInferCall << " \"" << inputReadGroupString << "\"";
        }
        else
        {
            mapAllInferCall << " NULL";
        }
        
        mapAllInferCall << " NULL " << (int) idDigitLen << " " << splitMode << " " << inputMode << " ";

        if(inputMode == "PAIRED")
        {
            //specify the insert size constraints
            mapAllInferCall << insertSizeLowerThreshold << " " << insertSizeUpperThreshold << " ";              
        }
        else
        {
            mapAllInferCall << "NULL NULL" << " ";
        }
        mapAllInferCall << "NULL " << collapsedFastqFile_longReadsList << " NULL " << numTotalReads << " " << memoizationThreshold << " " << numMismatchesPerReadMer << " " << numMismatchesPerReadMer;

        ExecuteSystemCall(mapAllInferCall, VERIFY);
    }

    double endTime = getTime();
    logOut << "Search time: " << endTime - constructTimeEnd << endl;
    logOut << "Time between BEGIN-END: " << endTime - beginTime << endl;
}
