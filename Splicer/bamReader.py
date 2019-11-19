import sqlite3
import os
import bamnostic as bn

#get connection to the sqlite database
conn = sqlite3.connect("E:\speedSplice" + os.path.sep + 'splice.sqlite', isolation_level=None)
c = conn.cursor()


samfile = bn.AlignmentFile("hg19test.bam", "rb")

i = 0
for read in samfile:
    cigar = read.cigarstring
    start = read.reference_start+read.cigar[0][1]
    stop = read.reference_start+read.cigar[0][1]+read.cigar[1][1]
    if "N" in cigar:
        i+=1
        c.execute("SELECT * FROM splice WHERE from_pos='"+start+"' AND to_pos='"+stop+"'")
        print (cigar, read.cigartuples, start, stop, read.reference_name)
    if i > 50:
        break


