#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <map>
#include <utility>
#include <algorithm>

using namespace std;

bool loadSampleMap(ifstream& ifs, const string& infilename, map<string, unsigned int>& mapOut);
bool loadSampleMap(ifstream& ifs, const string& infilename, map<string, unsigned int>& mapOut)
{
  const int BUFSIZE(1000);
  char buf[BUFSIZE];
  unsigned int linenum(0);

  while (ifs.getline(buf, BUFSIZE))
    {
      linenum++;
      if (strchr(buf,'\t'))
	{
	  cerr << "Error:  Tab-delimited data found on line " << linenum << " of " << infilename
	       << ".\nContents should be one column of IDs whose order defines the order of the output columns."
	       << endl << endl;
	  return false;
	}
      mapOut[string(buf)] = linenum;
    }
  return true;
}

bool loadChunkIDmap(ifstream& ifs, const string& infilename, map<string, unsigned long>& mapOut);
bool loadChunkIDmap(ifstream& ifs, const string& infilename, map<string, unsigned long>& mapOut)
{
  const int BUFSIZE(1000);
  char buf[BUFSIZE];
  unsigned int linenum(0);

  while (ifs.getline(buf, BUFSIZE))
    {
      linenum++;
      if (strchr(buf,'\t'))
	{
	  cerr << "Error:  Tab-delimited data found on line " << linenum << " of " << infilename
	       << ".\nContents should be one column of IDs whose order defines the order of the output rows."
	       << endl << endl;
	  return false;
	}
      mapOut[string(buf)] = linenum;
    }
  return true;
}

//   map<string, unsigned long> mapFromChunkIDtoOutputOrder;
bool parseChunk(ifstream& ifs, const string& infilename, const map<string, unsigned int>& sampToColNum,
		const map<string, unsigned long>& chunkIDtoRowOrder, map<unsigned long, map<unsigned int, string> >& outData);
bool parseChunk(ifstream& ifs, const string& infilename, const map<string, unsigned int>& sampToColNum,
		const map<string, unsigned long>& chunkIDtoRowOrder, map<unsigned long, map<unsigned int, string> >& outData)
{
  const int BUFSIZE(1000);
  char *p, buf[BUFSIZE];
  string sampName, signal, chunkID;
  map<string, unsigned long>::const_iterator cit;
  map<unsigned long, map<unsigned int, string> >::iterator odit;
  map<string, unsigned int>::const_iterator sit; // iterator over map from sample ID to output column
  map<unsigned int, string> emptyMap;
  map<unsigned int, string>::iterator signal_it; // iterator over (ID,signal) values observed within a DHS (i.e. within a chunkID)
  unsigned long linenum(0);
  unsigned int fieldnum(0);

  while (ifs.getline(buf, BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf,"\t")))
	{
	MissingField:
	  cerr << "Error:  Failed to find required field (column) " << fieldnum
	       << "on line " << linenum << " of " << infilename << '.' << endl << endl;
	  return false;
	}
      // ignore column 1, chromosome
      while ((p = strtok(NULL,"\t")))
	{
	  fieldnum++;
	  switch (fieldnum) {
	  case 4:
	    sampName = string(p);
	    break;
	  case 5:
	    signal = string(p);
	    break;
	  case 8:
	    chunkID = string(p);
	    break;
	  default:
	    break;
	  }
	}
      if (fieldnum < 8)
	goto MissingField;
      if (fieldnum > 8)
	{
	  cerr << "Error:  Found additional fields (columns) on line " << linenum
	       << " of " << infilename << "; expected exactly 8." << endl << endl;
	  return false;
	}
      cit = chunkIDtoRowOrder.find(chunkID);
      if (chunkIDtoRowOrder.end() == cit)
	{
	  cerr << "Error:  Chunk ID " << chunkID << ", on line " << linenum
	       << " of " << infilename << ", was not found in the list of ordered chunkIDs."
	       << endl << endl;
	  return false;
	}
      sit = sampToColNum.find(sampName);
      if (sampToColNum.end() == sit)
	{
	  cerr << "Error:  Sample name \"" << sampName << "\", in column 4 on line " << linenum
	       << " of " << infilename << ", was not found in the list of sample IDs." << endl << endl;
	  return false;
	}
      if (outData.empty())
	{
	  outData[cit->second] = emptyMap;
	  outData[cit->second][sit->second] = signal;
	  continue;
	}
      odit = outData.lower_bound(cit->second); // odit points to chunkID's output order, or the closest entry greater than it, or end()
      if (outData.end() == odit || odit->first != cit->second)
	{
	  // Insert an entry for cit->second, which is the output order for this chunkID. Use a "hint" for max-efficiency insertion.
	  if (odit != outData.begin())
	    odit--;
	  odit = outData.insert(odit, pair<unsigned long, map<unsigned int, string> >(cit->second, emptyMap));
	  odit->second[sit->second] = signal;
	  continue;
	}
      // Occasionally, a sample contributes multiple peaks to a DHS.
      // In such cases, report the maximum signal observed among the sample's peaks within the DHS.
      signal_it = odit->second.find(sit->second);
      if (odit->second.end() == signal_it)
	odit->second[sit->second] = signal;
      else
	{
	  if (atof(signal.c_str()) > atof(signal_it->second.c_str()))
	    signal_it->second = signal;
	}
    }

  return true;
}

void writeRows(ofstream& ofsSignal, ofstream& ofsBinary, const map<unsigned long, map<unsigned int, string> >& rowData, const unsigned int& nCols);
void writeRows(ofstream& ofsSignal, ofstream& ofsBinary, const map<unsigned long, map<unsigned int, string> >& rowData, const unsigned int& nCols)
{
  for (map<unsigned long, map<unsigned int, string> >::const_iterator rowData_it = rowData.begin();
       rowData_it != rowData.end(); rowData_it++)
    {
      unsigned int curCol(2); // curCol = current column in the output
      map<unsigned int, string>::const_iterator it = rowData_it->second.begin();
      // Write the first column separately, then precede each subsequent column with a tab character.
      if (1 == it->first)
	{
	  ofsSignal << (it++)->second;
	  ofsBinary << 1;
	}
      else
	{
	  ofsSignal << 0;
	  ofsBinary << 0;
	}
      while (it != rowData_it->second.end())
	{
	  while (curCol < it->first)
	    {
	      ofsSignal << "\t0";
	      ofsBinary << "\t0";
	      curCol++;
	    }
	  ofsSignal << '\t' << (it++)->second;
	  ofsBinary << "\t1";
	  curCol++;
	}
      while (curCol <= nCols)
	{
	  ofsSignal << "\t0";
	  ofsBinary << "\t0";
	  curCol++;
	}
      ofsSignal << endl;
      ofsBinary << endl;
    }
}

int main(int argc, const char* argv[])
{
  if (argc != 6)
    {
      cerr << "Usage:  " << argv[0] << " fileOfSampleIDsInOutputColumnOrder fileOfAllChunkIDsInRowOrder fileOfChunkDataWith8columns outfileSignal outfileBinary\n"
	   << "\twhere outputSignal will a (subset of a) signal matrix with no row or column labels (i.e., not a bed file)\n"
	   << "\tand outputBinary will be the same as outputSignal but with each nonzero entry replaced with 1.\n"
	   << "\tThe first input file must contain a single column of unique sample IDs that encompass all IDs\n"
	   << "\t\tthat can possibly occur in column 4 of fileOfChunkDataWith8columns.\n"
	   << "\tThe second input file must contain the 4th column of the appropriate *_chunkIDs.bed Index file\n"
	   << "\t\t(the chunkIDs in final row order), as a single column with no additional data.\n"
	   << "\tfileOfChunkDataWith8columns must be a chunk*bed file from the appropriate peaks_* subdirectory\n"
	   << "\t\tof the directory that received all the temporary files while the Index was being built."
	   << endl << endl;
      return -1;
    }
  ifstream infileSampleNames(argv[1]), infileChunkIDorder(argv[2]), infileChunkData(argv[3]);
  ofstream outfileSignal(argv[4]), outfileBinary(argv[5]);
  string fnameSampleNames(argv[1]), fnameChunkIDorder(argv[2]), fnameChunkData(argv[3]),
    fnameOutputSignal(argv[4]), fnameOutputBinary(argv[5]);
  if (!infileSampleNames)
    {
      cerr << "Error:  Unable to open file \"" << fnameSampleNames << "\" for read." << endl << endl;
      return -1;
    }
  if (!infileChunkIDorder)
    {
      cerr << "Error:  Unable to open file \"" << fnameChunkIDorder << "\" for read." << endl << endl;
      return -1;
    }
  if (!infileChunkData)
    {
      cerr << "Error:  Unable to open file \"" << fnameChunkData << "\" for read." << endl << endl;
      return -1;
    }
  if (!outfileSignal)
    {
      cerr << "Error:  Unable to open file \"" << fnameOutputSignal << "\" for write." << endl << endl;
      return -1;
    }
  if (!outfileBinary)
    {
      cerr << "Error:  Unable to open file \"" << fnameOutputBinary << "\" for write." << endl << endl;
      return -1;
    }    
  map<string, unsigned int> mapFromSampleNameToColumnNum;
  if (!loadSampleMap(infileSampleNames, fnameSampleNames, mapFromSampleNameToColumnNum))
    return -1;
  map<string, unsigned long> mapFromChunkIDtoOutputOrder;
  if (!loadChunkIDmap(infileChunkIDorder, fnameChunkIDorder, mapFromChunkIDtoOutputOrder))
    return -1;
  map<unsigned long, map<unsigned int, string> > outputRowData;
  if (!parseChunk(infileChunkData, fnameChunkData, mapFromSampleNameToColumnNum,
		  mapFromChunkIDtoOutputOrder, outputRowData))
    return -1;
  writeRows(outfileSignal, outfileBinary, outputRowData, mapFromSampleNameToColumnNum.size());

  return 0;
}
