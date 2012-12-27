#include "AS_UTL_fileIO.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

static string readRecord(ifstream &file, string &header,
        vector<int32> *qualValues = NULL) {
    string toReturn;
    int position = file.tellg();
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
        vector<int32> &qualValues, double qvCut, uint32 offset = 33) {
    string qual;
    uint32 minSeqLen = 50;

    double maxSum = 0;
    uint32 maxStartIndex = 0;
    uint32 maxEndIndex = 0;
    double currSum = (double) qualValues.front() - qvCut;
    uint32 currStartIndex = 0;

    for (size_t i = 0; i < qualValues.size(); i++) {
        qual += ((char) (qualValues[i] + offset));

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
    cout << "@" << header << "/" << maxStartIndex << "_" << maxEndIndex << endl;
    cout << fasta.substr(maxStartIndex, maxEndIndex - maxStartIndex + 1)
            << endl;
    cout << "+" << header << "/" << maxStartIndex << "_" << maxEndIndex << endl;
    cout << qual.substr(maxStartIndex, maxEndIndex - maxStartIndex + 1) << endl;
}

int main(int argc, char * argv[]) {

    bool debug = false;
    double qvCut = 59.5;
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <fasta file> <0-60 PHRED qual file>"
                << endl;
        exit(1);
    }
    if (argc > 4) {
        qvCut = atof(argv[3]);
    }

    ifstream fastaFile(argv[1], ifstream::in);
    ifstream qualFile(argv[2], ifstream::in);

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

        processRecord(header, fasta, qualValues, qvCut);
        qualValues.clear();
    }

    fastaFile.close();
    qualFile.close();
}

