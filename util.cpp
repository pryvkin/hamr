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

#include <string>
#include <vector>
#include <map>

#include "hamr.h"

using namespace std;

// parses a list of command line arguments where each argument is either
// <value> or --switch or --key=value
// the first type, positional arguments, go into positional_args
// the rest go into a map of key=value
void parse_arguments(const vector<string> &args,
		     arg_collection &value_args,
		     vector<string> &positional_args) {
  for(unsigned int i=0; i < args.size(); ++i) {
    if (args[i].empty())
      continue;
      
    if (args[i][0] == '-') {
      unsigned int equals_at = args[i].find('=');
      string key = args[i].substr(0, equals_at);
      string value = "";

      if (equals_at != string::npos &&
	  equals_at < (args[i].size() - 1))
	  value = args[i].substr(1+equals_at, string::npos );

      value_args[key] = value;

    } else {
      positional_args.push_back(args[i]);
    }
  }
}
