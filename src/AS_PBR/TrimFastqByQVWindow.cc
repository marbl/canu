#include "AS_UTL_fileIO.H"
#include "AS_UTL_fasta.H"

#include <iostream>
#include <fstream>
#include <vector>

#define FASTA_LINE_LENGTH 60

using namespace std;

static string readRecord(ifstream &file, string &header,
        vector<int32> *qualValues = NULL) {
    string toReturn;
    long position = file.tellg();
    string line;

    while (file.good()) {
        position = file.tellg();
        file >> line;
        if (line.at(0) == '>') { // record start
            if (toReturn.length() != 0) {
                file.seekg(position, ios::beg);
                break;
            }
            header = line.substr(1);
        } else {
            toReturn += line;
            if (qualValues != NULL) {
                toReturn += " ";
                qualValues->push_back(atol(line.c_str()));
            }
        }
    }

    return toReturn;
}

static void processRecord(string &header, string& fasta,
        vector<int32> &qualValues, FILE *fastaOut, FILE *qualOut,  double qvCut, uint32 minSeqLen = 50, uint32 offset = 33) {
    string qual;
    string qualCAEncoding;

    double maxSum = 0;
    uint32 maxStartIndex = 0;
    uint32 maxEndIndex = 0;
    double currSum = (double) qualValues.front() - qvCut;
    uint32 currStartIndex = 0;
    qual += ((char) (qualValues[0] + offset));
    qualCAEncoding += ((char) (qualValues[0] + '0'));

    for (size_t i = 1; i < qualValues.size(); i++) {
        qual += ((char) (qualValues[i] + offset));
        qualCAEncoding += ((char) (qualValues[i] + '0'));

        if (currSum > 0) {
            currSum = currSum + ((double) qualValues[i] - qvCut);
        } else {
            currSum = ((double) qualValues[i] - qvCut);
            currStartIndex = i;
        }
        if (currSum > maxSum) {
            maxSum = currSum;
            maxStartIndex = currStartIndex;
            maxEndIndex = i;
        }
    }
    if (maxEndIndex - maxStartIndex < minSeqLen) {
        return;
    }

   AS_UTL_writeFastA(fastaOut, (char *)fasta.substr(maxStartIndex, maxEndIndex - maxStartIndex + 1).c_str(), maxEndIndex - maxStartIndex + 1, FASTA_LINE_LENGTH, ">%s/"F_U32"_"F_U32"\n", (char *)header.c_str(), maxStartIndex, maxEndIndex);
   AS_UTL_writeQVFastA(qualOut, (char *)qualCAEncoding.substr(maxStartIndex, maxEndIndex - maxStartIndex + 1).c_str(), maxEndIndex - maxStartIndex + 1, (FASTA_LINE_LENGTH/3), ">%s/"F_U32"_"F_U32"\n", (char *)header.c_str(), maxStartIndex, maxEndIndex);

    cout << "@" << header << "/" << maxStartIndex << "_" << maxEndIndex << endl;
    cout << fasta.substr(maxStartIndex, maxEndIndex - maxStartIndex + 1)
            << endl;
    cout << "+" << header << "/" << maxStartIndex << "_" << maxEndIndex << endl;
    cout << qual.substr(maxStartIndex, maxEndIndex - maxStartIndex + 1) << endl;
}

int main(int argc, char * argv[]) {

    bool debug = false;
    double qvCut = 59.5;
    int32 minSeqLen = 50;
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <fasta file> <0-60 PHRED qual file> <fasta output> <qual output> [optional qv cutoff = " << qvCut << "]" << "[ optional min length = " << minSeqLen << "]"
                << endl;
        exit(1);
    }
    if (argc > 5) {
        qvCut = atof(argv[5]);
    }
    if (argc > 6) {
        minSeqLen = atoi(argv[6]);
    }

    ifstream fastaFile(argv[1], ifstream::in);
    ifstream qualFile(argv[2], ifstream::in);
    FILE *fastaOut = fopen(argv[3], "w");
    FILE *qualOut = fopen(argv[4], "w");

    string fasta;
    string qual;
    vector<int32> qualValues;
    string header;
    string qualHeader;

    while (fastaFile.good()) {
        fasta = readRecord(fastaFile, header);
        qual = readRecord(qualFile, qualHeader, &qualValues);
        if (header != qualHeader) {
            cerr << "Error, header files are not equal" << endl;
            exit(1);
        }

        processRecord(header, fasta, qualValues, fastaOut, qualOut, qvCut, minSeqLen);
        qualValues.clear();
    }

    fastaFile.close();
    qualFile.close();
    fclose(fastaOut);
    fclose(qualOut);
}

