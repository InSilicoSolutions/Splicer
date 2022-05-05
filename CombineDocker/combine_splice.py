import zipfile
import _pickle as cPickle
import numpy as np
import os
import sys
import gc

sample_files = []
sample_ids = []

if len(sys.argv) < 3:
    print("Usage combine_splice directory tissue")
    exit(-1)

directory = sys.argv[1]
tissue = sys.argv[2]

if os.path.exists(os.path.join(directory, tissue + '.zip')):
    os.remove(os.path.join(directory, tissue + '.zip'))

#Build a list of the samples in the directory
for (root,dirs,files) in os.walk(directory):
    for file in files:
        if file.endswith(".zip"):
            sample_files.append(os.path.join(root,file))
            if '_' in file:
                id = file[file.rindex('_')+1:]
            else:
                id = file
            if id.startswith("GTEX"):
                sample_ids.append(id.rpartition('.')[0])
            else:
                sample_ids.append(id.rpartition('-')[0])
                
num_samples = len(sample_files)


#Make on pass through the sample files to get a list of chromosomes,
#number of aligned reads per sample, and verify the sample id.
num_reads = np.zeros(num_samples, dtype=np.int32)
chrom_list = []
for position, sample_file in enumerate(sample_files):

    zf = zipfile.ZipFile(os.path.join(directory, sample_file), 'r') 
    if len(chrom_list) == 0:  
        #Get a list of the chromosomes we need to process. 
        splices = cPickle.load(zf.open('splice_count'))
        chrom_list = splices.keys()
    meta = cPickle.load(zf.open('metadata'))
    num_reads[position] = meta['Mapped Reads']  
    
    #Check that the sample file is the one we are expecting.  
    sample_id = meta['Sample']
    if sample_id.endswith('-TP'):
        sample_id = sample_id[0:-3]
    if (sample_ids[position] != sample_id):
        print("Error Sample name mismatch " + sample_ids[position] + " " + sample_id)
        exit(-1)
              
    zf.close()



#Open output file for combined data pickles
zOut = zipfile.ZipFile(os.path.join(directory, tissue + '.zip'), 'w', zipfile.ZIP_DEFLATED) 
#store the sample identifiers in output.  Position indicates relative position of sample in splice read counts.
zOut.writestr('samples', cPickle.dumps(sample_ids, -1))
#Store the total number of reads per sample
zOut.writestr('num_reads', cPickle.dumps(num_reads, -1))


#Combine the counts for splices across all samples.  A dictionary per chromosome is created
#and for each observed splice, a vector is put into the dictionary with a count of reads 
#for that splice.  A samples array identifies the relative position of each sample in the 
#read count arrays.  A num_reads array identifies the total mapped reads for each sample.


#loop one chromosome at a time to keep memory requirements for construction and use reasonable
for chr in chrom_list:
    combine_splices = {}
    print("Processing " + chr + " for " + str(num_samples) + " splice count sample files.")
    for position, sample_file in enumerate(sample_files):
     
        zf = zipfile.ZipFile(os.path.join(directory, sample_file), 'r')    
    
        splices = cPickle.load(zf.open('splice_count'))
        
        if chr in splices:
            chr_dict = splices[chr]
            for splice_from in chr_dict:
                splice_from_dict = chr_dict[splice_from]
                for splice_to in splice_from_dict:
                    splice_str = str(splice_from) + '-' + str(splice_to)
                    if splice_str not in combine_splices:
                        combine_splices[splice_str] = np.zeros(num_samples, dtype=np.int32)
                    combine_splices[splice_str][position] = splice_from_dict[splice_to]
      
        zf.close()



    print ("Writing chromosome " + chr)
    zOut.writestr(chr, cPickle.dumps(combine_splices, -1))
    del combine_splices
    gc.collect()
    

zOut.close()
print('Done')


