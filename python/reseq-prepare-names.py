#!/usr/bin/python
# Changes read names, so that pairs have identical names and tiles are not stripped of during the mapping

import getopt
import gzip
import string
import sys

# Ignore broken pipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def modName(name, num_spaces, num_colons):
    split_line = name[:-1].split(' ')
    split_part = split_line[num_spaces].split(':')
    part = ':'.join(split_part[:num_colons+1])
    start = '_'.join(split_line[:num_spaces])
    return start+'_'+part

def prepareNames(file1,file2):
    with gzip.open(file2, 'rb') if 'gz' == file2.split('.')[-1] else open(file2, 'r') as f2:
        first_line2 = f2.readline()
        pass

    with gzip.open(file1, 'rb') if 'gz' == file1.split('.')[-1] else open(file1, 'r') as f1:
        # Find the space and the semicolon where to separate
        first_line1 = f1.readline()

        num_spaces = 0;
        num_colons = 0;
        last_num_colons = 0;
        read_pos = 0;
    
        while(read_pos < min(len(first_line1),len(first_line2)) and first_line1[read_pos] == first_line2[read_pos] and (2 > num_colons or ' ' != first_line1[read_pos]) ):
            if ' ' == first_line1[read_pos]:
                num_spaces += 1
                last_num_colons = num_colons
                num_colons = 0
                pass
            if ':' == first_line1[read_pos]:
                num_colons += 1
                pass
                
            read_pos += 1
            pass
            
        # In case we ended with a difference, reduce it by one so the difference is not included
        if read_pos < min(len(first_line1),len(first_line2)) and first_line1[read_pos] != first_line2[read_pos]:
            if num_colons:
                num_colons -= 1
                pass
            elif num_spaces:
                num_spaces -= 1
                num_colons = last_num_colons
                pass
            pass
            
        # Pass file1 to stdout
        print modName(first_line1, num_spaces, num_colons)
        n_lines = 1 # We already have the first line
        for line in f1:
            if 0 == n_lines % 4:
                print modName(line, num_spaces, num_colons)
                pass
            else:
                print line[:-1]
                pass
                
            n_lines += 1
            pass
        pass

    return

def usage():
    print "Usage: python reseq-prepare-names.py [OPTIONS] File1 File2"
    print "Returns File1 in stdout with changed read names, so that pairs have identical names and tiles are not stripped of during the mapping"
    print "  -h, --help            display this help and exit"
    return

def main(argv):
    try:
        optlist, args = getopt.getopt(argv, 'h', ['help'])
        pass
    except getopt.GetoptError:
	print "Unknown option\n"
        usage()
        sys.exit(2)
        pass


    for opt, par in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
            pass
        pass

    if 2 != len(args):
        print "Wrong number of files. Exactly two are required.\n"
        usage()
        sys.exit(2)
        pass

    prepareNames(args[0],args[1])
    return

if __name__ == "__main__":
    main(sys.argv[1:])
    pass
