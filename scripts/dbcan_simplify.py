#!/usr/bin/env python

import os
import sys

if len(sys.argv) == 1:
    print("usage: script.py cazy_annotations.tsv")
    sys.exit(1)

with open(sys.argv[1]) as fin:
    next(fin)
    for line in fin:
        cazy = set()
        line = line.rstrip()
        cols = line.split("\t")
        for col in cols[1:-1]:
            if col != "-":
                col = col.split("+")
                for enz in col:
                    enz = enz.split("(")[0].split("_")[0]
                    if not enz[0].isdigit():
                        cazy.add(enz)
        print("%s\t%s" % (cols[0], "\t".join(list(cazy))))
                   
                
        
