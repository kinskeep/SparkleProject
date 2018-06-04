#!/usr/bin/python
#4/7/2014

#Find overlapping SNPs in 2 files
import sys;

file1 = sys.argv[1];
file2 = sys.argv[2];
outfile = sys.argv[3];

output = open(outfile, "w");

goodSNPs = {};

input = open(file1);
line = input.readline().strip();
while line != "":
	parts = line.split();
	if parts[1] == "pVal":
		line = input.readline().strip();
		continue;
	if parts[0] in goodSNPs:
		goodSNPs[parts[0]] = goodSNPs[parts[0]] + 1;
	else:
		goodSNPs[parts[0]] = 1;
	line = input.readline().strip();
input.close();

input = open(file2);
line = input.readline().strip();
count = 0;
while line != "":
        parts = line.split();
        if parts[0] in goodSNPs:
                output.write(parts[0] + "\n");
		count = count + 1;
	line = input.readline().strip();
input.close();
output.close();
print count;
