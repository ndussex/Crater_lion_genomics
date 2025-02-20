#!/usr/bin/python
from __future__ import division
from sys import argv
import sys
import subprocess
from optparse import OptionParser
from itertools import izip, islice
import random
import operator
import gzip
import math
from Bio import SeqIO

"""
author: Tom van der Valk
email: tom.vandervalk@nrm.se
usage: python vcf_replace_ref_allel.py --gzvcf filename.vcf.gz --outgroup outgroup.fasta --output DIR/output_name
"""

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("--gzvcf", action="store", type="string", dest="gzvcf")
parser.add_option("--outgroup-fasta", action="store", type="string", dest="outgroup")
parser.add_option("--output", action="store", type="string", dest="output")
(options, args) = parser.parse_args()

print "loading outgroup fasta...."
input_file = options.outgroup
chrom_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

print "iterating through vcf file"

outputfile = gzip.open(options.output +".vcf.gz","w")

counter = 0
total_sites = 0
with gzip.open(options.gzvcf) as f1:
        for line in f1:
            if line.startswith("#"):
                outputfile.write(line)
            else:
                counter += 1
                if counter % 100000 == 0:
                    print "processed variants....: ", counter
                splitted = line.strip().split("\t")
                chrom,pos,ref_allel,alt_allel,quality,alleles = splitted[0],int(splitted[1]),splitted[3],splitted[4],splitted[5],splitted[9:]
                if "2" in alleles:
                    print "more than 2 alleles found at site:",chrom,pos,"skipping variant"
                else:
                    if chrom not in chrom_dict:
                        print "chromosome:", chrom,"not found in fasta file"
                    else:
                        ancesteral_allel = chrom_dict[chrom][pos-1]
                        if ancesteral_allel != "N":
                            total_sites += 1
                            if ancesteral_allel == ref_allel:
                                outputfile.write(line)
                            elif ancesteral_allel == alt_allel:
                                parsed_alleles = []
                                for i in alleles:
                                    i = i.split(":")[0]
				    if i == "0/0":
                                        parsed_alleles += ["1/1"]
                                    elif i == "0/1":
                                        parsed_alleles += ["0/1"]
                                    elif i == "1/1":
                                        parsed_alleles += ["0/0"]
                                outputfile.write(chrom + "\t" + str(pos) + "\t" + "." + "\t" +  alt_allel + "\t" + ref_allel + "\t" + quality + "\t" + "." + "\t" + \
                                                splitted[7] + "\t" + splitted[8] + "\t" + "\t".join(parsed_alleles) + "\n")

print "total variants:",counter
print "sites without ancesteral state (filtered): ",counter-total_sites
outputfile.close()
