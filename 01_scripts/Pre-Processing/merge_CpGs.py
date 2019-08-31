#!/usr/bin/python

import sys, argparse

def main(argv):

        parser = argparse.ArgumentParser(description='Script takes a bed input and merges the 2 strands of a CpG to one')
        parser.add_argument('-i', required=True, help='Input bed file', metavar='Input bed file', dest='inputfile')
        parser.add_argument('-o', required=True, help='Output file', metavar='Output file', dest='outputfile')
        args = parser.parse_args()

        f_out = open(args.outputfile, 'w')
        f_in =  open(args.inputfile, 'r')
	
        line1_again = 0
        number_of_lines_before = 0
        number_of_lines_after = 0
	
	# Read Header
        f_in.readline()

        while 1:
        	# Set line 1
                if (not line1_again):
                    line1 = f_in.readline().split('\t')
                    if len(line1) == 1:
                        break
                    number_of_lines_before += 1
                line1_again = 0

		# Set line 2
                line2 = f_in.readline().split('\t')
                if len(line2) == 1:
                    assert(int(line1[1]) + 1 == int(line1[2]))
                    f_out.write(line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + str(float(line1[3])/100) + '\t' + line1[4] + '\t' + line1[5] + '\t' + line1[6] + '\t' + line1[7] + '\t' + line1[8] + '\t' + line1[9] + '\t' + line1[10])
                    number_of_lines_after += 1
                    break
                number_of_lines_before += 1

		# Merge lines if possible
		# (1) Is one position in the + and one on the # strand? (2) Are the positions subsequent and on the same chromosome?
                if (line1[5] == '+' and line2[5] == '-') and (int(line1[2]) == int(line2[1])) and (line1[0] == line2[0]):
                    methylation = ((float(line1[3])/100) * float(line1[4]) + (float(line2[3])/100) * float(line2[4]))/(float(line1[4]) + float(line2[4]))
                    assert(int(line1[2]) == int(line2[1]))
                    f_out.write(line1[0] + '\t' + line1[1] + '\t' + line2[2] + '\t' + str(methylation) + '\t' + str(int(line1[4]) + int(line2[4])) + '\t' + '+-' + '\t' + line1[6] + '\t' + line1[7] + '\t' + line1[8] + '\t' + line1[9] + '\t' + line1[10])
                    number_of_lines_after += 1
		# 2 lines are not one CpG
                else:
                    assert(int(line1[1]) + 1 == int(line1[2]))
                    if (line1[5] == '+'):
                        f_out.write(line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + str(float(line1[3])/100.0) + '\t' + line1[4] + '\t' + line1[5] + '\t' + line1[6] + '\t' + line1[7] + '\t' + line1[8] + '\t' + line1[9] + '\t' + line1[10])
                    elif (line1[5] == '-'):
                        f_out.write(line1[0] + '\t' + str(int(line1[1])-1) + '\t' + str(int(line1[2])-1) + '\t' + str(float(line1[3])/100.0) + '\t' + line1[4] + '\t' + line1[5] + '\t' + line1[6] + '\t' + line1[7] + '\t' + line1[8] + '\t' + line1[9] + '\t' + line1[10])            
                    number_of_lines_after += 1
                line1 = line2
                line1_again = 1
		
	#print '\nMerging \t\t\t done \t --> Reduced from ' + str(number_of_lines_before) + ' lines to ' + str(number_of_lines_after) + '\n'
        f_in.close()
        f_out.close()		

if __name__ == "__main__":
   main(sys.argv[1:])
