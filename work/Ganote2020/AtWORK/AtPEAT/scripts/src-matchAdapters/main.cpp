/*
 * main.cpp
 *
 *  Created on: Aug 12, 2011
 *      Author: Andrea Gossett
 *      Purpose: Try to recreate the IdentifyPairsFaster.java file
 *
 */

// Load header files
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <map>
#include <new>
#include <cstdlib>
//#include "PairedRead.h"

// Declare std usage
using std::cout;
using std::endl;
using std::string;

// Declare functions
void Usage();
int matchAdapter(string str1, string adapter, int adapterLength, int maxMismatches, 
		int minReportLength, int maxReportLength, int counter);

// Declare the main function
int main(int argc, char *argv[]) {

//	set_new_handler(NoMemory);

	// Declare variables
	char input1[1000];
	char input2[1000];
	char outputMatch[1000];
	char outputNonMatch[1000];

	string format = "fa";
	string headerFormat = "new";
	string adapterFivePrime = "GTTGGACTCGAGCGTACATCGTTAGAAGCTTGAG";
	string adapterThreePrime = "GTCGGATAGGCCGTCTTCAGCCGC";
	int minReportLength = 18;
	int maxReportLength = 200;
	int maxMismatches = 1;
	int adapterLength = 12;

	//typedef map<string, PairedRead> readInformation;

	// Variables for file reading
	string lineOne1;
	string lineOne2;
	string lineOne3;
	string lineOne4;
	string lineTwo1;
	string lineTwo2;
	string lineTwo3;
	string lineTwo4;
	int lcount = 0;
	int nMatch = 0;
	int nNoMatch = 0;
	int nLineOneMatch = 0;
	int nLineTwoMatch = 0;
	int nFivePrimeMatch = 0;
	int nThreePrimeMatch = 0;
	int nFivePrimeTooShort = 0;
	int nThreePrimeTooShort = 0;
	int nLineOneTooShort = 0;
	int nLineTwoTooShort = 0;
	int nLineOneNoMatch = 0;
	int nLineTwoNoMatch = 0;

//	int lineLimit = 1000000;

	std::fstream fin1, fin2;
	std::fstream fout1, fout0;

	if(argc < 1) {
		Usage();	// print usage statement
	} else {
		for(int i = 0; i < argc; i++) {
			if(strcmp(argv[i], "-i") == 0) {
				strcpy(input1, argv[i+1]);
				strcpy(input2, argv[i+2]);
			} else if(strcmp(argv[i], "-oM") == 0) {
				strcpy(outputMatch, argv[i+1]);
			} else if(strcmp(argv[i], "-oN") == 0) {
				strcpy(outputNonMatch, argv[i+1]);
			} else if(strcmp(argv[i], "-aF") == 0) {
				adapterFivePrime = string(argv[i+1]);
			} else if(strcmp(argv[i], "-aT") == 0) {
				adapterThreePrime = string(argv[i+1]);
			} else if(strcmp(argv[i], "-mm") == 0) {
				maxMismatches = atoi(argv[i+1]);
			} else if(strcmp(argv[i], "--fa") == 0) {
				format = string("fa");
			} else if(strcmp(argv[i], "--fq") == 0) {
				format = string("fq");
			} else if(strcmp(argv[i], "--old") == 0) {
				headerFormat = string("old");
			} else if(strcmp(argv[i], "--new") == 0) {
				headerFormat = string("new");
			} else if(strcmp(argv[i], "-l") == 0) {
				adapterLength = atoi(argv[i+1]);
			} else if(strcmp(argv[i], "-mr") == 0) {
				minReportLength = atoi(argv[i+1]);
			}
		}
	}

	if((strlen(input1) == 0) || (strlen(input2) == 0) || (strlen(outputMatch) == 0)) {
		cout << "One or more file names was not specified.\n";
		Usage();
		return 0;
	}

	// open input/output files; if one fails, close
	fin1.open(input1, std::ios::in);
	if(!fin1.is_open()) {
		cout << "Error opening input1:\t" << input1 << "\n";
		return 0;
	}
	fin2.open(input2, std::ios::in);
	if(!fin2.is_open()) {
		cout << "Error opening input2:\t" << input2 << "\n";
		fin1.close();
		return 0;
	}
	fout1.open(outputMatch, std::ios::out);
	if(!fout1.is_open()) {
		cout << "Error opening output1:\t" << outputMatch << "\n";
		fin1.close();
		fin2.close();
		return 0;
	}
	if(strlen(outputNonMatch) > 0) {
		fout0.open(outputNonMatch, std::ios::out);
		if(!fout0.is_open()) {
			cout << "Error opening output1:\t" << outputNonMatch << "\n";
			fin1.close();
			fin2.close();
			fout1.close();
			return 0;
		}
	}

	// Read in lines one by 1 from fin1 and fin2
 	while(getline(fin1, lineOne1)) {

		if(lcount%10000000 == 0) {
			cout << "Currently matching line " << lcount << "\n";
		}

		// Read in lines 1-2 no matter what from file 1
		getline(fin1, lineOne2);

		// Read in lines 1-2 no matter what from file 2
		getline(fin2, lineTwo1);
		getline(fin2, lineTwo2);

		// If the input format was fastq, read in another 2 lines
		if(format.find("fq") != string::npos) {
			getline(fin1, lineOne3);
			getline(fin1, lineOne4);
			getline(fin2, lineTwo3);
			getline(fin2, lineTwo4);
		}

		int lineOneFive = 0;
		int lineOneThree = 0;
		int lineTwoFive = 0;
		int lineTwoThree = 0;

		// SeqID is line 1
		// Sequence is line 2
		// SeqID is line 3
		// Qual is line 4

		// Check to see that line 1 is the same in file1 and file2
		if((headerFormat.compare("old") == 0) && (lineOne1.substr(0, (lineOne1.length()-2)).compare(
				lineTwo1.substr(0, (lineTwo1.length()-2))) != 0)) {
			cout << "File 1 and 2's sequence IDs did not match at line " << lcount << "\n";
			return 0;
		} else if ((headerFormat.compare("new") == 0) && (lineOne1.substr(0,(lineOne1.find_first_of(" "))).compare(
				lineTwo1.substr(0, (lineTwo1.find_first_of(" ")))))) {
			cout << "File 1 and 2's sequence IDs did not match at line " << lcount << "\n";
			return 0;
		} else {
		// If the string IDs are the same
			
			// We need to see if either sequence contains the adapter
			lineOneFive = matchAdapter(lineOne2, adapterFivePrime, adapterLength,
						maxMismatches, minReportLength, maxReportLength, lcount);
			if(lineOneFive == 0) {
				lineOneThree = matchAdapter(lineOne2, adapterThreePrime, adapterLength,
						maxMismatches, minReportLength, maxReportLength, lcount);
			}
	/*
			if(lineOneFive != 0){
				cout << "Found 5' adapter in file 1 (" << lineOneFive << ").\n";
			} else if (lineOneThree != 0) {
				cout << "Found 3' adapter in file 1 (" << lineOneThree << ").\n";
			}
	*/
			lineTwoFive = matchAdapter(lineTwo2, adapterFivePrime, adapterLength,
					maxMismatches, minReportLength, maxReportLength, lcount);
			if(lineTwoFive == 0) {
				lineTwoThree = matchAdapter(lineTwo2, adapterThreePrime, adapterLength,
									maxMismatches, minReportLength, maxReportLength, lcount);
			}
		/*
			if(lineTwoFive != 0){
				cout << "Found 5' adapter in file 2 (" << lineTwoFive << ").\n";
			} else if (lineTwoThree != 0) {
				cout << "Found 3' adapter in file 2 (" << lineTwoThree << ").\n";
			}

			cout << "-----\n";
		*/
			if((lineOneThree != 0) || (lineOneFive != 0)) {
				nLineOneMatch++;
			} else {
				nLineOneNoMatch++;
			}
			if((lineTwoThree != 0) || (lineTwoFive != 0)) {
				nLineTwoMatch++;
			} else {
				nLineTwoNoMatch++;
			}

			if((lineOneThree != 0) || (lineTwoThree != 0)) {
				nThreePrimeMatch++;
			}
			if((lineOneFive != 0) || (lineTwoFive != 0)) {
				nFivePrimeMatch++;
			}

			// If either side's read is too short
			if(((lineTwoThree != 0) && (lineTwoThree < minReportLength)) ||
				((lineOneThree != 0) && (lineOneThree < minReportLength))) {
				nThreePrimeTooShort++;
			}
			if(((lineTwoFive != 0) && (lineTwoFive < minReportLength)) ||
				((lineOneFive != 0) && (lineOneFive < minReportLength))) {
				nFivePrimeTooShort++;
			}

			if(((lineOneThree != 0) && (lineOneThree < minReportLength)) ||
				((lineOneFive != 0) && (lineOneFive < minReportLength))) {
				nLineOneTooShort++;
			}
			if(((lineTwoThree != 0) && (lineTwoThree < minReportLength)) ||
				((lineTwoFive != 0) && (lineTwoFive < minReportLength))) {
				nLineTwoTooShort++;
			}


			// If the format is "old", then only subtract 2 bp
			int subtractLength = 4;
			if(headerFormat.compare("old") == 0) {
				subtractLength = 2;
			}

			// If both reads had adapters
			if(((lineOneFive >= minReportLength) || (lineOneThree >= minReportLength)) && 
				((lineTwoFive >= minReportLength) || (lineTwoThree >= minReportLength))) {
				// Output whichever reads
				if(lineOneFive >= minReportLength) {
					fout1 << lineOne1.substr(0, lineOne1.length()-subtractLength) << ":F\n";
					fout1 << lineOne2.substr(0, lineOneFive) << "\n";
					if(format.find("fq") != string::npos) {
						fout1 << lineOne3.substr(0, lineOne1.length()-subtractLength) << ":F\n";
						fout1 << lineOne4.substr(0, lineOneFive) << "\n";
					}
				} else if(lineOneThree >= minReportLength) {
					fout1 << lineOne1.substr(0, lineOne1.length()-subtractLength) << ":T\n";
					fout1 << lineOne2.substr(0, lineOneThree) << "\n";
					if(format.find("fq") != string::npos) {
						fout1 << lineOne3.substr(0, lineOne1.length()-subtractLength) << ":T\n";
						fout1 << lineOne4.substr(0, lineOneThree) << "\n";
					}
				}
				
				if(lineTwoFive >= minReportLength) {
					fout1 << lineTwo1.substr(0, lineTwo1.length() - subtractLength) << ":F\n";
					fout1 << lineTwo2.substr(0, lineTwoFive) << "\n";
					if(format.find("fq") != string::npos) {
						fout1 << lineTwo3.substr(0, lineTwo1.length() - subtractLength) << ":F\n";
						fout1 << lineTwo4.substr(0, lineTwoFive) << "\n";
					}
				} else if(lineTwoThree >= minReportLength) {
					fout1 << lineTwo1.substr(0, lineTwo1.length() - subtractLength) << ":T\n";
					fout1 << lineTwo2.substr(0, lineTwoThree) << "\n";
					if(format.find("fq") != string::npos) {
						fout1 << lineTwo3.substr(0, lineTwo1.length() - subtractLength) << ":T\n";
						fout1 << lineTwo4.substr(0, lineTwoThree) << "\n";
					}
				}
				nMatch++;
			// If they don't both contain the adapter...
			} else {
				if(fout0.is_open()) {
					if(lineOneFive >= minReportLength) {
						fout0 << lineOne1.substr(0, lineOne1.length()-subtractLength) << ":F\n";
						fout0 << lineOne2.substr(0, lineOneFive) << "\n";
					} else if (lineOneThree >= minReportLength) {
						fout0 << lineOne1.substr(0, lineOne1.length()-subtractLength) << ":T\n";
						fout0 << lineOne2.substr(0, lineOneThree) << "\n";
					} else {
						fout0 << lineOne1.substr(0, lineOne1.length()-subtractLength) << ":N\n";
						fout0 << lineOne2 << "\n";
					}
					if(lineTwoFive >= minReportLength) {
						fout0 << lineTwo1.substr(0, lineTwo1.length()-subtractLength) << ":F\n";
						fout0 << lineTwo2.substr(0, lineTwoFive) << "\n";
					} else if (lineTwoThree >= minReportLength) {
						fout0 << lineTwo1.substr(0, lineTwo1.length()-subtractLength) << ":T\n";
						fout0 << lineTwo2.substr(0, lineTwoThree) << "\n";
					} else {
						fout0 << lineTwo1.substr(0, lineTwo1.length()-subtractLength) << ":N\n";
						fout0 << lineTwo2 << "\n";
					}
				}
				nNoMatch++;
			}

			lcount++;
		}
	}

	// Close the files that we opened
	fin1.close();
	fin2.close();
	fout1.close();
	if(fout0.is_open()) {
		fout0.close();
	}

	cout << "Matched " << nMatch << " pairs.\n";
	cout << nNoMatch << " pairs did not contain both a 5' and 3' read and a sufficiently long trimmed sequence.\n";
	cout << "File breakdown:\n";
	cout << "\tMatched an adapter in " << nLineOneMatch << " reads from the first file.\n";
	cout << "\tMatched an adapter in " << nLineTwoMatch << " reads from the second file.\n";
	cout << "Side breakdown:\n";
	cout << "\tMatched a 5' adapter in " << nFivePrimeMatch << " reads.\n";
	cout << "\tMatched a 3' adapter in " << nThreePrimeMatch << " reads.\n";
	cout << "Length breakdown:\n";
	cout << "\t" << nFivePrimeTooShort << " reads were too short on 5' side.\n";
	cout << "\t" << nThreePrimeTooShort << " reads were too short on 3' side.\n";
	cout << "\t" << nLineOneTooShort << " reads were too short in file 1.\n";
	cout << "\t" << nLineTwoTooShort << " reads were too short in file 2.\n";

	return 0;
}

// Define a usage function to provide help to the user
void Usage() {
	cout << "Parameters for getReadsWithAdapters:\n";
	cout << "Last modified 2012/03/21\n";
	cout << "\t-i <file1> <file2>\tFile names for the two PE read files.\n";
	cout << "\t--fa/--fq\t\tInput format (FASTA, default, or FASTQ).\n";
	cout << "\t--old/--new\tHeader format (old = #0/File; new = spaceFile:Y:0).\n";
	cout << "\t-oM <file>\t\tOutput filenames for matched reads.\n";
	cout << "\t-oN <file>\t\tOutput filenames for non-matched reads.\n";
	cout << "\t-aF <string>\t\t5' adapter sequence (default: GTTGGACTCGAGCGTACATCGTTAGAAGCTTGA).\n";
	cout << "\t-tF <string>\t\t3' adapter sequence (default: GTCGGATAGGCCGTCTTCAGCCGCTCAAGCTTC).\n";
	cout << "\t-mm <int>\t\tMax mismatches allowed (default: 1).\n";
	cout << "\t-l <int>\t\tLength of adapter to use (default: 12).\n";
	cout << "\t-mr <int>\t\tMinimum length of string to count as a match (default: 18).\n";
}

int matchAdapter(string str1, string adapter, int adapterLength, int maxMismatches, 
		int minReportLength, int maxReportLength, int counter) {
	int badMatch = 0;

	if(str1.find(adapter.substr(0, adapterLength)) != string::npos) {
			return str1.find(adapter.substr(0, adapterLength));
	} else {
		// If it doesn't, we're going to have to check individually
		for(unsigned int i = 0; i < (str1.length() - adapterLength); i++) {
			// create a test string of same length as adapter
			string teststr = str1.substr(i, adapterLength);

			badMatch = 0;

			// Test base by base
			for(int j = 0; j < (adapterLength); j++) {
				if(!teststr.substr(j,1).compare(adapter.substr(j,1)) == 0){
					if(badMatch <= maxMismatches+1) {
						badMatch++;
					} else {
						j = adapterLength;
					}
				}
			}

			if(badMatch <= maxMismatches) {
				return i;
			}
		}
	}
	
	return 0;
}

