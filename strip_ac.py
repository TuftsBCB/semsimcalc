#!/usr/bin/python

import argparse

def main():
	""" Parses files from Swiss-Prot, found here:
		http://www.uniprot.org/downloads
		Tested on the Swiss-Prot file text file.
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--infile", help="Input AC file", required=True)
	parser.add_argument("-o", "--outfile", help="Output stripped AC file", required=True)

	options = parser.parse_args()

	infile = open(options.infile)
	outfile = open(options.outfile, 'w')

	for line in infile:
		if line.startswith("AC"):
			name = line.strip().strip(';').split()[1]
			outfile.write(name)
			outfile.write('\n')
			outfile.flush()
		if (line.startswith("DR") and "GO:" in line):
			cols = line.strip().strip(';').split(';')
			term = [s for s in cols if ("GO:" in s)][0]
			outfile.write(term.strip())
			outfile.write('\n')
			outfile.flush()
		if line.startswith("ID"):
			outfile.write("-\n")
			outfile.flush

	outfile.write('-\n')
	outfile.flush()

	infile.close()
	outfile.close()

if __name__ == "__main__":
	main()
