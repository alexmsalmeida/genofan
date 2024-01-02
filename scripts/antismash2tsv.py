#!/usr/bin/env python

import os
import glob
import sys
from Bio import SeqIO

genome_name = sys.argv[2]
genome_to_bgc = {}
for input_file in glob.glob(os.path.join(sys.argv[1], "*.gbk")):
    genome_name = "_".join(os.path.basename(input_file).split("_")[:-1])
    with open(input_file) as file_in:
       for record in SeqIO.parse(file_in, "genbank"):
           for feature in record.features:
               if feature.type == "cand_cluster":
                   for cluster in feature.qualifiers["product"]:
                       print("%s\t%s" % (genome_name, cluster))
    
