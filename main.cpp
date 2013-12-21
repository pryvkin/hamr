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
#include <string>
#include <vector>

#include "hamr.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "USAGE: " << argv[0] << " cmd\n" 
	 << "    where cmd is rnapileup|filter_pileup|rnapileup2mismatchbed\n";
    return(1);
  }

  vector<string> args;

  // make args[0] be "hamr_cmd <cmd>"
  string cmd = argv[1];
  args.push_back(string(argv[0]) + " " + cmd);
  for(int i=2; i < argc; ++i)
    args.push_back(argv[i]);

  if (cmd == "rnapileup")
    return (rnapileup_main(args));
  else if (cmd == "filter_pileup")
    return (filter_pileup_main(args));
  else if (cmd == "rnapileup2mismatchbed")
    return (rnapileup2mismatchbed_main(args));
  else {
    cerr << "Invalid command: " << cmd << "\n";
    return(1);
  }
  

  return(0);
}
