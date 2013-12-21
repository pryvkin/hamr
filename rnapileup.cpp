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

////  rnapileup
// Generates pileup for RNAseq data:
//      Requires a sorted BAM file
//      Genomic strands are treated separately unless --no-strand-specific
//      Does not support indels

//  2.0 - Added another column with position-along-read data
//        (Sanger encoded just like base quals; X-33 = position along
//          read from 5' end, starting at 0)
//  2.1 - Now supports soft-clipping and discards reads with indels
//  2.2 - Integrated into HAMR
//        Added several options for filtering input, yielding
//          a smaller output file
    
// #define DEBUGMODE

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <deque>
#include <cctype>

#include "sam.h"
#include "faidx.h"
#include "hamr.h"

using namespace std;

///////////////////////

// represents pileup data at one site
struct Pileup {
  int pos;
  char ref;
  int nreads;
  string pileup;
  string quals;
  string readpos;

  // for debugging
#ifdef DEBUGMODE
  vector<string> read_ids;
#endif
  Pileup() : pos(0), ref('N'), nreads(0), pileup(), quals() { }
  Pileup(int p) : pos(p), ref('N'), nreads(0), pileup(), quals() { }
};

void process_queue(deque<Pileup> &q, int upto_pos, bool process_all,
		   const string &ref_id,
		   int min_coverage,
		   unsigned long &sites_excluded_cov,
		   unsigned long &sites_encountered) {
  // output one-based coords
  while( (!q.empty()) &&
	 (process_all || (q.front().pos < upto_pos))) {

    ++sites_encountered;

    // exclude sites with not enough reads covering
    if (q.front().nreads < min_coverage) {
      ++sites_excluded_cov;
      q.pop_front();
      continue;
    }

    cout << ref_id << "\t"
	 << 1+(q.front().pos) << "\t"
	 << q.front().ref << "\t"
	 << q.front().nreads << "\t"
	 << q.front().pileup << "\t"
	 << q.front().quals << "\t" 
	 << q.front().readpos;
#ifdef DEBUGMODE
    cout << "\t";
    for (int i=0; i < q.front().read_ids.size(); ++i) {
      cout << q.front().read_ids[i] << ",";
    }
#endif
      cout << "\n";
    q.pop_front();
  }
}

/////////////////////
/*class DNAComplementer {
  char c[256];
public:
  DNAComplementer() { 
    for(int i=0; i<256; ++i)
      c[i] = 'n';
    c['a'] = 't'; c['c'] = 'g';  c['g'] = 'c'; c['t'] = 'a';
    c['m'] = 'k'; c['r'] = 'y';  c['k'] = 'm'; c['y'] = 'r';
    c['s'] = 's'; c['w'] = 'w';
    c['v'] = 'b'; c['h'] = 'd'; c['b'] = 'v'; c['d'] = 'h';
    for(int i='a'; i<='z'; ++i)
      c[toupper(i)] = char(toupper(c[int(i)]));
  }
  char operator () (char x) const { return c[int(x)]; }
  };*/

/////////////////////

void print_usage(const vector<string> &args, bool options_only=false) {
  if (!options_only) {
    cerr << "USAGE: " << args[0] << " [OPTIONS] reads.bam genome.fasta\n\n"
	 << "    OPTIONS:\n";
  }
  cerr   << "      --exclude-ends         Exclude 5' and 3' ends of reads\n"
         << "      --min-q=N              Exclude bases with Q score < N (15)\n"
         << "      --min-coverage=N       Exclude sites with < N reads covering (10)\n"
	 << "      --not-strand-specific  Library not strand-specific (convert everything to +)\n";

}

int rnapileup_main(const vector<string> &args) {
  arg_collection value_args;
  vector<string> positional_args;

  parse_arguments(args, value_args, positional_args);

  // library not stand specific
  bool no_ss = false;
  bool exclude_ends = false;
  int min_coverage = 10;
  int min_q = 15;

  // collect and validate command line arguments
  for (arg_collection::iterator it = value_args.begin();
       it != value_args.end(); ++it) {
    string key = it->first;
    string value = it->second;
    bool conv_success = false;

    if (key == "--not-strand-specific") {
      no_ss = true;

    } else if (key == "--exclude-ends") {
      exclude_ends = true;

    } else if (key == "--min-q") {
      min_q = from_s<int>(value, conv_success);
      if (!conv_success || (min_q < 0)) {
	cerr << "Invalid value for --min-q: " << value << "; must be a non-negative integer\n";
	return(1);
      }

    } else if (key == "--min-coverage") {
      min_coverage = from_s<int>(value, conv_success);
      if (!conv_success || (min_coverage < 0)) {
	cerr << "Invalid value for --min-coverage: "
	     << value << "; must be a non-negative integer\n";
	return(1);
      }
    } else if (key == "--list-options") {
      print_usage(args, true);
      return(0);

    } else {
      /*cerr << "Invalid option: " << key << "\n";
      print_usage(args);
      return(1);*/
    }
  }

  if (positional_args.size() < 3) {
    print_usage(args);
    return(1);
  }

  string bam_fn( positional_args[1] );
  string fas_fn( positional_args[2] );

  // output supplied arguments
  cerr << "  Processing BAM file " << bam_fn << "\n";
  cerr << "  Using genome fasta file " << fas_fn << "\n";
  if (no_ss)
    cerr << "  Treating library as non-stand-specific\n";
  if (exclude_ends)
    cerr << "  Excluding ends of reads\n";
  cerr << "  Requiring Q-score >= " << min_q << "\n";
  cerr << "  Requiring " << min_coverage << " coverage at a site\n";

  // index the fasta file by finding out where each chr starts
  ifstream file(fas_fn.c_str());
  if (!file.is_open()) {
    cerr << "Failed to open FASTA file " << fas_fn << "\n";
    return 1;
  }

  faidx_t *fai = fai_load(fas_fn.c_str());

  bamFile bam_file;
  bam_header_t *bam_hdr;

  if ((bam_file = bam_open(bam_fn.c_str(), "r")) == 0) {
    cerr << "Failed to open BAM file " << bam_fn << "\n";
    return 1;
  }

  // read BAM header
  bam_hdr = bam_header_read(bam_file);
  bam1_t *bam = bam_init1();

  string curr_ref;
  string prev_ref;
  char *ref_seq = NULL;
  int ref_len(0);
  bool changed_ref = false;

  string bases(16, 'X');
  bases[1] = 'A';
  bases[2] = 'C';
  bases[4] = 'G';
  bases[8] = 'T';
  bases[15] = 'N';

  // track numbers for filtered bases
  unsigned long bases_excluded_end = 0;
  unsigned long bases_excluded_q = 0;
  unsigned long bases_encountered = 0;
  unsigned long sites_excluded_cov = 0;
  unsigned long sites_encountered = 0;


  // maintain a queue of pileup data and output sites (process_queue)
  // when we encounter a read that starts after them
  deque<Pileup> q;

  while( bam_read1(bam_file, bam) > 0 ) {
    changed_ref = false;

    // skip non-unique reads
    //int num_hits = bam_aux2i( bam_aux_get(bam, "NH") );
    //if (num_hits > 1)
    //  continue;

    // skip unmapped reads
    if ((bam->core.flag & 0x4) > 0)
      continue;

    // check for soft-clipped reads (clipped seq is present in SEQ field but missing from reference)
    // on both ends
    bool has_indels(false);
    int nclipstart(0);
    int nclipend(0);
    
    // check the CIGAR string for various operations
    for(int i=0; i < bam->core.n_cigar; ++i) {
      // discards reads with indels
      if (bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CINS || 
	  bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CDEL )  {
	has_indels = true;
	break;
      }

      // note soft-clipping so we can ignore those bases
      if (bam_cigar_op(bam1_cigar(bam)[i]) == BAM_CSOFT_CLIP) {
	if (i == 0)
	  nclipstart = bam_cigar_oplen(bam1_cigar(bam)[i]);
	else
	  nclipend = bam_cigar_oplen(bam1_cigar(bam)[i]);
      }
    }

    if (has_indels)
      continue;

    // get chr name for this bam line
    string ref( bam_hdr->target_name[bam->core.tid] );

    // load genomic sequence for this chromosome
    if (ref != curr_ref) {
      cerr << "Loading sequence for " << ref << " ...";

      if (ref_seq) {
	// deallocate previous chr sequence
	free(ref_seq);
	changed_ref = true;
	prev_ref = curr_ref;
      } else {
	// initialize prev_ref to ref
	prev_ref = ref;
      }

      curr_ref = ref;
      ref_seq = fai_fetch(fai, ref.c_str(), &ref_len);

      // convert sequence to uppercase
      for(int i=0; i < ref_len; ++i)
	ref_seq[i] = toupper(ref_seq[i]);

      cerr << "done.\n";
      cerr.flush();
    }
    
    int read_pos(bam->core.pos);
    int read_len(bam->core.l_qseq);
    string read_qual(read_len, '#');
    string read_seq(read_len, 'N');

#ifdef DEBUGMODE
    string read_id(bam1_qname(bam));
#endif

#ifdef DEBUGMODE
    if (nclipstart != 0 || nclipend != 0)
      cerr << "Soft-clipped: (" << nclipstart << ", " << nclipend << "\n";
#endif

    // remove all the soft-clipped bases from the start
    read_seq.erase(0, nclipstart);
    read_qual.erase(0, nclipstart);
    read_len -= nclipstart;
    
    // process queue
    process_queue(q, read_pos, changed_ref, 
		  changed_ref ? prev_ref : curr_ref,
		  min_coverage,
		  sites_excluded_cov,
		  sites_encountered);

    // build read seq
    for(int i=0; i < (read_len-nclipend); ++i) {
      read_seq[i] = bases[bam1_seqi( bam1_seq(bam), nclipstart + i)];
      read_qual[i] = 33+bam1_qual(bam)[nclipstart + i];
    }
    //cout << "Read: " << read_seq << "\n";

    bool rev_strand = bam1_strand(bam);

    // for this read, simultaneously loop through its sequence
    // and the pileup queue
    deque<Pileup>::iterator q_it( q.begin() );
    for(int i=0; i < (read_len-nclipend); ++i, ++q_it) {
      int g(read_pos + i);  // genomic position

      // reported POS in bam is actually the first MATCHING base, so adjust it by the starting soft clip
      //g -= nclipstart;

      if (q_it == q.end()) {
	q.push_back(Pileup(g));
	q_it = q.end()-1;
      }

      // skip clipped bases
      if ( i < nclipstart || i > ((read_len-nclipend)-1) )
	continue;

      ++bases_encountered;

      // exclude read-ends
      if (exclude_ends && 
	  ((i == 0) || (i == (read_len - 1))) ) {
	++bases_excluded_end;
	continue;
      }
      // exclude low-quality bases
      if ((int(read_qual[i])-33) < min_q) {
	++bases_excluded_q;
	continue;
      }

      if (q_it->pos != g) {
	cerr << "ERROR: read pos " << g << " != queue pos " << q_it->pos << "\n";
#ifdef DEBUGMODE
	cerr << read_pos << " " << read_id << "\n";
#endif
	return 1;
      }

      // cout << read_seq[i] << " vs " << ref_seq[g] << "\n";
     
      if (i == 0)
	q_it->pileup += (rev_strand && !no_ss) ? "$" : "^~";
      if (i == ((read_len-nclipend) - 1))
	q_it->pileup += (rev_strand && !no_ss) ? "^~" : "$";

      // make sure we don't go past end of ref seq
      if (g >= ref_len) {
	cerr << "ERROR: genomic pos " << g << " >= chr length (" << ref_len << ")\n";
	return 1;
      }

      if (ref_seq[g] == read_seq[i])
	q_it->pileup += (rev_strand && !no_ss) ? ',' : '.';
      else
	q_it->pileup += (rev_strand && !no_ss) ? tolower(read_seq[i]) : read_seq[i];

      q_it->readpos += char(33 + ((rev_strand && !no_ss) ? (read_len-(1+i)) : i));
      q_it->ref = ref_seq[g];
      q_it->quals += read_qual[i];

#ifdef DEBUGMODE
      q_it->read_ids.push_back(read_id);
#endif

      ++(q_it->nreads);
    }
  }

  // process queue
  process_queue(q, 0, true, curr_ref, min_coverage,
		sites_excluded_cov, sites_encountered);

  // output statistics
  double bases_excluded_end_pct = 100.0 * double(bases_excluded_end) / 
    double(bases_encountered);
  double bases_excluded_q_pct = 100.0 * double(bases_excluded_q) / 
    double(bases_encountered);
  double sites_excluded_cov_pct = 100.0 * double(sites_excluded_cov) / 
    double(sites_encountered);

  cerr << "Bases encountered: " << bases_encountered << "\n"
       << "Bases excluded due to being on read-end: " << setw(3) << bases_excluded_end_pct << "%\n"
       << "Bases excluded due to low Q: " << setw(3) << bases_excluded_q_pct << "%\n"
       << "Sites encountered: " << sites_encountered << "\n"
       << "Sites excluded due to low coverage: " << setw(3) << sites_excluded_cov_pct << "%\n";

  return 0;
}
