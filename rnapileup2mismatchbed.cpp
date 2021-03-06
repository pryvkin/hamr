//  Copyright (c) 2013 University of Pennsylvania
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "hamr.h"

using namespace std;

int rnapileup2mismatchbed_main( const vector<string> &args ) {
  if (args.size() < 2) {
    cerr << "USAGE: " << args[0] << " in.rnapileup\n";
    return(1);
  }

  // arg! have to do these annoying tricks to get an input stream
  // to come from stdin or a file
  istream* p_infile;
  ifstream* p_in;
  if (args[1] != "-") {
    p_in = new ifstream(args[1].c_str());
    if (!p_in->is_open()) {
      cerr << "Could not open file " << args[1] << "\n";
      return(1);
    }
    p_infile = p_in;
  }
  else
    p_infile = &cin;

  istream& infile = (*p_infile);

  vector<char> relevant_nucs;
  relevant_nucs.push_back(',');
  relevant_nucs.push_back('.');
  relevant_nucs.push_back('A');
  relevant_nucs.push_back('C');
  relevant_nucs.push_back('G');
  relevant_nucs.push_back('T');
  relevant_nucs.push_back('N');
  relevant_nucs.push_back('a');
  relevant_nucs.push_back('c');
  relevant_nucs.push_back('g');
  relevant_nucs.push_back('t');
  relevant_nucs.push_back('n');

  bool rev_strand[256];
  memset(rev_strand, 0, sizeof(bool)*256);
  rev_strand['a'] = rev_strand['c'] = rev_strand['g'] = rev_strand['t'] = 
    rev_strand['n'] = rev_strand[','] = true;

  char complement[256];
  memset(complement, 0, sizeof(char)*256);
  complement['A'] = complement['a'] = 'T';
  complement['C'] = complement['c'] = 'G';
  complement['G'] = complement['g'] = 'C';
  complement['T'] = complement['t'] = 'A';
  complement['N'] = complement['n'] = 'N';
  complement['.'] = complement[','] = '.';


  string line;
  while (getline(infile, line)) {
    istringstream linestr(line);
    string chr, posstr, refstr, nstr;
    getline(linestr, chr, '\t');
    getline(linestr, posstr, '\t');
    unsigned int pos = atoi(posstr.c_str());

    getline(linestr, refstr, '\t');
    char ref=refstr[0];

    getline(linestr, nstr, '\t');
    //int n = atoi(nstr.c_str());

    string nucstr;
    getline(linestr, nucstr, '\t');

    string qualstr;
    getline(linestr, qualstr, '\t');

    string readposstr;
    getline(linestr, readposstr, '\t');

    int counts[256];
    memset(counts, 0, sizeof(int)*256);

    string readposstr_by_nuc[256];

    int read_idx=0;
    for(unsigned int i=0; i < nucstr.size(); ++i) {
      if (nucstr[i] == '^')
	i += 2; // skip ^ and mapq
      else if (nucstr[i] == '$')
	i += 1; // skip $
      ++ counts[(unsigned int)nucstr[i]];
      readposstr_by_nuc[ (unsigned int)nucstr[i] ] += readposstr[read_idx];
      ++ read_idx;
    }

    for(unsigned int i=0; i < relevant_nucs.size(); ++i) {
      char strand = '+';
      char new_ref = ref;
      char nuc = relevant_nucs[i];
      char orig_nuc = nuc;
      unsigned int nuc_idx = (unsigned int)nuc;
      unsigned int ref_nuc_idx = (unsigned int)ref;
      unsigned int orig_nuc_idx = (unsigned int)orig_nuc;

      int count = counts[nuc_idx];
      if (count > 0) {
	if (rev_strand[nuc_idx]) {
	  strand = '-';
	  new_ref = complement[ref_nuc_idx];
	  nuc = complement[nuc_idx];
	}
	// condense read position string (which is Sanger encoded)
	// into a histogram of the form x:count,x:count,...
	int readpos_counts[256];
	memset(readpos_counts, 0, sizeof(int)*256);

	// track max_readpos to speed up the loop below
	int max_readpos=0;
	for(unsigned int i=0; i < readposstr_by_nuc[orig_nuc_idx].size(); ++i) {
	  int readpos = int(readposstr_by_nuc[orig_nuc_idx][i])-33;
	  ++ readpos_counts[readpos];
	  max_readpos = (readpos > max_readpos) ? readpos : max_readpos;
	}

	// output BED format
	cout << chr << "\t" << (pos-1) << "\t" << pos << "\t"
	     << new_ref << ">" << nuc << "\t" << count << ";";
	bool first=true;
	for(int i=0; i <= max_readpos; ++i) {
	  if (readpos_counts[i] > 0) {
	    if (!first)
	      cout << ",";
	    cout << i << ":" << readpos_counts[i];
	    first=false;
	  }
	}
      
	cout << "\t" << strand << "\n";
      }
    }
  }

  if (p_infile != &cin)
    delete p_infile;

  return(0);
}

