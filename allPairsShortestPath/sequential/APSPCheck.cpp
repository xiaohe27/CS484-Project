//
// Created by xiaohe on 4/18/16.
//
// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "utils.h"
#include "apsp.h"
#include <assert.h>

using namespace std;


int main(int argc, char* argv[]) {
    string iFileStr = "randLocalGraph_WE_5_3";
    const char *iFile = iFileStr.c_str();
    string oFileStr = "randLocalGraph_WE_5_3.tmp";
    const char *oFile = oFileStr.c_str();

    string expectedFilePath = "randLocalGraph_WE_5_3.expected";

    write_allPairsShortestPath_2_file(strdup(iFile), strdup(oFile));
    bool eq = utils::areTwoFilesEqual(oFileStr, expectedFilePath);
    assert(eq);
}
