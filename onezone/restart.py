#!/usr/bin/python

from onezone import zone
import sys, getopt


def main(argv):

    inputfile = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o")
    except getopt.GetoptError:
        print('restart.py -i <pickled_output_filename>')
        sys.exit(2)

    for opt, arg in opts:

        if opt == '-h':
            print('Restart onezone model from a picked dump file')
            print('restart.py -i <pickled_output_filename>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        else:
            print("option not understood")

    print("Restarting from file ", inputfile)

    zone.restart(inputfile)    


if __name__ == "__main__":
    main(sys.argv[1:])

