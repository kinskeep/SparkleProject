#!/usr/bin/python
#Remove individuals with missingness above given threshold
#MMD 3/10/15

import sys, re, os;

#Print error/usage message if correct arguments not provided.
if len(sys.argv) != 3:
        print "ERROR: include 2 arguments.";
        print "Usage: python RemoveMissing.py <input> <threshold>";
        sys.exit(0);

#Declare variables from arguments input
infile = sys.argv[1];
thresh = sys.argv[2];
cutoff = float(thresh)/100;
outfile = infile.rstrip("str") + thresh + ".str";

print cutoff;

input = open(infile);
output = open(outfile, "w");

#Loop over file calculating %missing
line = input.readline().strip();
while line != "":
	parts = line.split();
        if parts[0] == "SNP_1":
                output.write(line + "\n");
                line = input.readline().strip();
                continue;
	length = len(parts);
	totSNPs = 0.0;
	missSNPs = 0.0; 
	for i in xrange(2,length):
		totSNPs = totSNPs + 1;
		if parts[i] == "-9":
			missSNPs = missSNPs + 1;
		else:
			missSNPs = missSNPs;
        if missSNPs/totSNPs <= cutoff:
                output.write(line + "\n");
                line = input.readline().strip();
        else:
		print missSNPs/totSNPs;
                line = input.readline().strip();

input.close();
output.close();
