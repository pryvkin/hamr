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
#include <algorithm>

#include "hamr.h"

using namespace std;

struct pileup_pos {
  bool start, end;
  char mapq, nuc, q, p;
  pileup_pos() : start(false), end(false), mapq('~'), nuc(' '), q(' '), p(' ')
   { }
};

struct ThresholdQuality {
  int minq;
  ThresholdQuality(int mq) : minq(mq) { }
  bool operator() (const pileup_pos &p) const {
    return ( int((p.q)-33) < minq );
  }
};

// remove positions corresponding to 5' or 3' ends of reads
// which can be error-prone
struct FilterEnds {
  FilterEnds() { }
  bool operator() (const pileup_pos &p) const {
    return ( p.start || p.end );
  }
};

int filter_pileup_main( const vector<string> &args ) {
  if (args.size() < 4) {
    cerr << "USAGE: " << args[0] << " in.rnapileup  minQ min_coverage [remove_ends]\n";
    return(1);
  }

  // arg have to do these tricks to get "-"=>stdin to work
  istream *p_infile;
  ifstream *p_in;

  if (args[1] != "-") {
    p_in = new ifstream(args[1].c_str());
    if (!p_in->is_open()) {
      cerr << "Could not open file " << args[1] << "\n";
      return(1);
    }
    p_infile = p_in;
  } else
    p_infile = &cin;

  bool do_remove_ends = false;
  if (args.size() == 4 && args[4] == "1") {
    cerr << "removing ends\n";
    do_remove_ends = true;
  }

  istream &infile = *p_infile;

  int minq = atoi(args[2].c_str());
  unsigned int min_cov = atoi(args[3].c_str());

  ThresholdQuality threshold_quality(minq);
  FilterEnds filter_ends;

  string line;
  while (getline(infile, line)) {
    istringstream linestr(line);
    string chr, pos, ref, nstr;
    getline(linestr, chr, '\t');
    getline(linestr, pos, '\t');
    getline(linestr, ref, '\t');
    getline(linestr, nstr, '\t');
    unsigned int n = atoi(nstr.c_str());

    if (n < min_cov)
      continue;

    vector<pileup_pos> data(n);
    string nucstr, qstr, pstr;
    getline(linestr, nucstr, '\t');
    getline(linestr, qstr, '\t');
    getline(linestr, pstr, '\t');
    
    int read = 0;  // 0..n
    for(unsigned int i=0; i < nucstr.size(); ++i) {
      if (nucstr[i] == '^') {
	data[read].start = true;
	data[read].mapq = nucstr[++i];
	data[read].nuc = nucstr[++i];
      } else if (nucstr[i] == '$') {
	data[read].end = true;
	data[read].nuc = nucstr[++i];
      } else {
	data[read].nuc = nucstr[i];
      }

      data[read].q = qstr[read];
      data[read].p = pstr[read];
      ++read;
    }

    data.erase(remove_if(data.begin(), data.end(), threshold_quality),
	       data.end());
    
    if (do_remove_ends)
      data.erase(remove_if(data.begin(), data.end(), filter_ends),
		 data.end());

    if (data.size() < min_cov)
      continue;

    cout << chr << "\t" << pos << "\t" << ref << "\t"
	 << data.size() << "\t";

    string new_quals(data.size(), ' ');
    string new_pos(data.size(), ' ');
    for(unsigned int i=0; i < data.size(); ++i) {
      new_quals[i] = data[i].q;
      new_pos[i] = data[i].p;
      if (data[i].start)
	cout << "^" << data[i].mapq;
      else if (data[i].end)
	cout << "$";
      cout << data[i].nuc;
    }

    cout << "\t" << new_quals << "\t" << new_pos << "\n";
  }

  if (p_infile != &cin)
    delete p_infile;

  return(0);
}

