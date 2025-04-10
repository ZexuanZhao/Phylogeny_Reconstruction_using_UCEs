#!/usr/bin/env python3
"""
 Gblocks.py - Python wrapper script to run Gblocks from BioLegato.

 Synopsis:
   
        Gblocks.py infile [options]
   

@modified: March 2 2021
@author: Brian Fristensky
@contact: brian.fristensky@umanitoba.ca
"""

import argparse
import os
import subprocess
import sys


PROGRAM = "Gblocks.py: "
USAGE = "\n\t USAGE: Gblocks.py infile [options] "

DEBUG = True

NUMSEQS=0 # Number of sequences in infile

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"Wrapper class for command line parameters"
class Parameters:

    def __init__(self):
        """
                Initializes arguments:
                Then calls read_args() to fill in their values from command line
                """

        self.IFN = "" 
        self.TYPE = "p"
        self.B3 = 8 
        self.B4 = 10
        self.B5 = "n"
        self.B6 = "y"
        self.V = "60"
        self.NUMSEQS = 0
        self.read_args()

    def read_args(self):
        """
                Read command line arguments into a Paramters object
                """
        parser = argparse.ArgumentParser()
        parser.add_argument("infile", action="store", default="", help="input file")
        parser.add_argument("--type", action="store", default="p", help="[p|d|c] protein,DNA,Codons")
        parser.add_argument("--b1", action="store", default="51", help="Min. % of Sequences For A Conserved Position")
        parser.add_argument("--b2", action="store", default="85", help="Min. % of Sequences for a flank Posn.")
        parser.add_argument("--b3", action="store", default="8", help="Max. No. of Contiguous Nonconserved Positions")
        parser.add_argument("--b4", action="store", default="10", help="Min. length of a block")
        parser.add_argument("--b5", action="store", default="n", help="Allowed gap posns.")
        parser.add_argument("--b6", action="store", default="8", help="")
        parser.add_argument("--v", action="store", default="50", help="")

        try:
            args = parser.parse_args()
            
            self.IFN = args.infile
            if args.type in ["p","d","c"] :
                self.TYPE = args.type
            if not args.b1 == "" :
                self.B1 = int(args.b1)
            if not args.b2 == "" :
                self.B2 = int(args.b2)
            if not args.b3 == "" :
                self.B3 = int(args.b3)
            if not args.b4 == "" :
                tempB4=int(args.b4)
                if tempB4 > 2 :
                    self.B4 = tempB4
                else :
                    self.B4 = 10
            if args.b5 in ["n","h","a"] :
                self.B5 = args.b5
            if args.b6 in ["y","n"] :
                self.B6 = args.b6
            tempV = int(args.v)
            if tempV >= 50 :
                self.V = tempV
            self.NUMSEQS = FastaCount(self.IFN)

            # Sanity checking.
            if self.B1 < 51 :
                self.B1 = 51
            if self.B1 > 100 :
                self.B1 = 100
            if self.B2 < 85 :
                self.B2 = 85
            if self.B2 > 100 :
                self.B2 = 100
            # From Gblocks paper, "conserved (>=IS and < FS)"
            if self.B1 > self.B2 :
                self.B2 = self.B1

        except ValueError:
            print(USAGE)

        if DEBUG :
            print("IFN: " + self.IFN)
            print("TYPE: " + self.TYPE)
            print("B1: " + str(self.B1))
            print("B2: " + str(self.B2))
            print("B3: " + str(self.B3))
            print("B4: " + str(self.B4))
            print("B5: " + self.B5)
            print("B6: " + self.B6)
            print("V: " + str(self.V))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"Return the number of sequences in a fasta file"
def FastaCount(FN):
    N = 0
    fafile = open(FN,"r")
    lines = fafile.readlines()
    fafile.close()
    for l in lines :
        if l.startswith(">") :
            N +=1
    return N

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"Run Gblocks"
def RunGblocks(P):

        COMMAND=["Gblocks",P.IFN, "-e=.fsa", "-t=" + P.TYPE]
        MinConSeqs = round((P.B1/100) * P.NUMSEQS)
        MinFlank = round((P.B2/100) * P.NUMSEQS)
        COMMAND.extend(["-b1=" +str(MinConSeqs), "-b2=" +str(MinFlank)])
        COMMAND.extend(["-b3=" + str(P.B3), "-b4=" + str(P.B4), "-b5=" + P.B5, "-b6=" + P.B6, "-v=" + str(P.V)])
        print(COMMAND)
        p = subprocess.Popen(COMMAND)
        p.wait()

#======================== MAIN PROCEDURE ==========================
def main():
    """
        Called when not in documentation mode.
        """
    P = Parameters ()

    RunGblocks(P)        

if ("-test" in sys.argv):
    pass
else:
    main()
