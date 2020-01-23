#!/usr/bin/python

import sys, argparse, os

def process(inputfile, outputfile):
    print(inputfile)
    f_in = open(inputfile, 'r')
    f_out = open(outputfile, 'w')
    f_out.write("chr\tpos\tN\tX\n")	

    for line in f_in:
        cols = line.split('\t')
        f_out.write(cols[0] + '\t' + cols[1] + '\t' + cols[4] + '\t' + str(int(0.5 + float(cols[3]) * float(cols[4]))) + '\n')
    f_in.close()
    f_out.close()
        #print 'Only CpGs with a minimum number of reads are concidered (before: ' + str(cpgs_before) + ' CpGs, now: ' + str(cpgs_after) + ' CpGs)'

def main(argv):
        parser = argparse.ArgumentParser(description='The Script takes a bed file as input and converts it to the format chr\tpos\tN\tX\t(N = total number of reads, X = number of reads showing methylation).')
        parser.add_argument('-i', required=True, help='Input bed file', metavar='Input bed file', dest='inputfile')
        parser.add_argument('-o', required=True, help='Output file', metavar='Output file', dest='outputfile')
        args = parser.parse_args()

        if not os.path.isfile(args.inputfile):
            print('Error: Input file is not valid.\tExecution will be stopped.')
            sys.exit()
        process(args.inputfile, args.outputfile)
        
main(sys.argv)
