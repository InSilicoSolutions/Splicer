import os
import zipfile
import _pickle as cPickle
import time
import gzip
import gc
import argparse

import numpy as np
import pysam as ps

#count the matching sequence in the rest of a read starting at the 'start' element in the cigar string
def remaining_match(cigartups, start):
    match = 0
    for i in range(start, len(cigartups)):
        tup = cigartups[i]
       
        if tup[0] in [0, 7]:
            match += tup[1]
    return match        

#For a give read, update the running count of splice junctions observed.
def tally_splices(read):
    #skip secondary alignments, poor quality reads, and unmapped reads
    if not read.is_secondary and not read.is_qcfail and not read.is_unmapped:
        
        cigar = read.cigarstring
        #If this is a junction read
        if "N" in cigar:
        
            #print (bn.utils.flag_decode(read.flag))
 
            match_bases = 0
            pos = read.reference_start
            chrom = read.reference_name
                
            for i, tup in enumerate(read.cigartuples):
                #read the 'N' items in the cigar string but require at least 5 base match on each side of the splice
                if tup[0] == 3:

                    if match_bases >= 5 and remaining_match(read.cigartuples, i) >= 5:
                        splice_5prime = pos
                        splice_3prime = pos+tup[1]+1
                        

                        #print(read.query_name, chrom + ' ' + str(splice_5prime) + '-' + str(splice_3prime), cigar, pos)
                        if chrom not in splice_counts:
                            splice_counts[chrom] = {}
                        if splice_5prime not in splice_counts[chrom]:   
                            splice_counts[chrom][splice_5prime] = {}
                        if splice_3prime not in splice_counts[chrom][splice_5prime]:
                            splice_counts[chrom][splice_5prime][splice_3prime] = 1
                        else:
                            splice_counts[chrom][splice_5prime][splice_3prime] = splice_counts[chrom][splice_5prime][splice_3prime] + 1
                    match_bases = 0
                else:   
                    match_bases = 0
                
                #Cigar String letters M D N = X that advance the reference position
                if tup[0] in [0,2,3,7,8]:
                    pos = pos + tup[1]
            
                if tup[0] in [0, 7]:
                    match_bases = match_bases + tup[1]
    
    
#For a give read, update the running count of read depth at each 
#location on the current chromosome where the read aligns. 
def tally_exons(read):
    
    cigar = read.cigarstring        

    if not read.is_qcfail and not read.is_unmapped:

        #Not sure what to do with secondary alignments - keep for now if < 4 mismatches.
        #This means ambiguous reads will increase counts in all locations where they align.
        mismatch = read.get_tag('NM')

        #Require moderately high identity and low mismatch/clipping.
        if mismatch <= 7 and read.query_length - (read.query_alignment_end - read.query_alignment_start) <= 7:

            pos = read.reference_start+1
            chrom = read.reference_name

            for i, tup in enumerate(read.cigartuples):
                #For matching sections of the cigar string, increase counts for that section of reference.
                if tup[0] == 0 or tup[0] == 7:
                    end = pos + tup[1]
                    countsArray[pos : end] += np.ones(tup[1], dtype=np.uint16)

                #Cigar String letters M D N = X that advance the reference position
                if tup[0] in [0,2,3,7,8]:
                    pos = pos + tup[1]
                

def estimateReadLen(bamfile):
    length = 0
    counter = 0
    for read in bamfile:
        counter += 1
        length += read.query_length
        if counter == 100:
            break
    length = length / 100
    return length
    
    

parser = argparse.ArgumentParser(description='Takes the given .bam file and looks through all the reads to construct a count of all exons and splices. Produces dictionary of splices and vector / chromosome of read counts', usage='splice_counter <bam_file>')
parser.add_argument("bam", help='The .bam file to be processed.  Note a companion .bai file is expected to be present')
parser.add_argument("-s", '--sample', nargs='?', help='Optional - sample identifier. Used to name output .zip file if present')
args = parser.parse_args()

start = time.time()
chr_time = time.time()

# Nested dictionary structure keyed by chromosome, then splice 5' end, then splice 3' end
# with a running count of how many times the splice is observed.
splice_counts = {}


bamfile = ps.AlignmentFile(args.bam, "rb")
#Assume sample ID is encoded in bam file name before the first '.'
if args.sample is None:
    sample = os. path. basename(args.bam)
    sample = sample[0:sample.find('.')]
else:
    sample = args.sample


#Open a zip file archive for output
zfile = zipfile.ZipFile(sample + '.zip', 'w', zipfile.ZIP_DEFLATED)

print("Procesing ", args.bam, "  Sample: " + sample)
print("   Total Aligned Reads: ", bamfile.mapped)



count = 0

#For speed and memory size, data is process one chromosome at a time.  A np array is used
#to aggregate read counts for each position in the genome.  So, iterate on chromosome (ref).
for pos, ref in enumerate(bamfile.references):
    if 'chrUn' in ref or 'random' in ref:
        continue

    #Allocate array to accumulate read counts per position
    countsArray = np.zeros(bamfile.lengths[pos]+1, dtype=np.int32)
    print("   Starting: ", ref)

    for read in bamfile.fetch(ref):
        count+=1

        tally_splices(read)
        
        tally_exons(read)
    
        if count % 1000000 == 0:
            print ("  Processed: ", count, " reads.")
   
    #write the chromsome read depth array to the output zip file.
    zfile.writestr(ref, cPickle.dumps(countsArray, -1)) 
    
    #free up memory used by the array
    del countsArray
    gc.collect()
    
    print ("  Completed: ", ref, " Time: ", str(time.time() - chr_time))
    chr_time = time.time()

print ("Total Alignments: " + str(count))  

#write splice counts to output file.
zfile.writestr('splice_count', cPickle.dumps(splice_counts, -1))

meta = {}
meta['Sample'] = sample
meta['Mapped Reads'] = bamfile.mapped
meta['Read length'] = estimateReadLen(bamfile)
meta['Alignments'] = count

zfile.writestr('metadata', cPickle.dumps(meta, -1))
zfile.close()

bamfile.close()

end = time.time()
print("Finished Elapsed: " + str(end - start))

