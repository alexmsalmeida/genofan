#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
import logging

def run_hmmscan(threads, output, db, query):
    cmd_hmmscan = ["hmmsearch",
                   "--cpu", threads,
                   "--noali",
                   "--cut_ga",
                   "--domtblout", output,
                   db, query]
    subprocess.call(cmd_hmmscan)

def parse_hmmscan(hmm_file):
    kos = {}
    with open(hmm_file) as fin:
        for line in fin:
            if line[0] != "#":
                cols = line.rstrip().split()
                if cols[0] not in kos:
                    kos[cols[0]] = [cols[3]]
                else:
                    kos[cols[0]].append(cols[3])
    return kos

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Query proteins against KOfam database')
    parser.add_argument('-t', '--threads', metavar='', required=True, \
                        help='Number of threads')
    parser.add_argument('-q', '--query', metavar='', required=True,\
                        help='Protein FASTA file to query')
    parser.add_argument('-d', '--database', metavar='', required=False,\
                        help='KOfam database files (.hmm)')
    parser.add_argument('-o', '--output', metavar='', required=True,\
                        help='Directory name to store output files')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        logging.basicConfig(format='%(asctime)s\t%(message)s', 
                            datefmt='%d/%m/%Y %I:%M:%S %p', level=logging.INFO)
        logging.info("Starting KOfamScan script ...")
        if not os.path.exists(args.output):
            os.makedirs(args.output)
        logging.info("Running hmmscan ...")
        query_basename = os.path.basename(args.query).split(".fa")[0]
        hmm_file = args.output+"/kofam_raw.tsv"
        run_hmmscan(args.threads, hmm_file, args.database, args.query)
        logging.info("Parsing results ...")
        results = parse_hmmscan(hmm_file)
        with open(args.output+"/kegg_orthologs.tsv", "w") as fout:
            for gene in results:
                fout.write("%s\t%s\n" % (gene, "\t".join(set(results[gene]))))
