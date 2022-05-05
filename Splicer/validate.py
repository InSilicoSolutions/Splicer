import zipfile
import _pickle as cPickle
import numpy as np

#zf = zipfile.ZipFile('C:/Splicer/Data/TCGA_OV.zip', 'r')    
#zf = zipfile.ZipFile('C:/Splicer/Data/GTEX_ovary_GTEX-13OVI-0726-SM-5L3DD.zip', 'r')    
 
#print('Reading')
#ov_splices = cPickle.load(zf.open('combined_splices'))

#ov_meta = cPickle.load(zf.open('metadata'))

#zf.close()

#print('looking')
#print(ov_splices['chr16']['766490-766668'])
#print('iso')
#print(ov_splices['chr16']['766490-766644'])
#print('done')

zf = zipfile.ZipFile('C:/Splicer/combined/TCGA_OV.zip', 'r')    
np.set_printoptions(precision=2, suppress=True)

print('Reading')
splices = cPickle.load(zf.open('chr16'))
num_reads = cPickle.load(zf.open('num_reads'))
reads_per_m = num_reads/1000000
zf.close()
print('meso')
meso = splices['766490-766668']
print(meso)
print('iso meso')
iso_meso = splices['766490-766644']
print(iso_meso)
print()
print()
print('normalized meso')
print(meso/reads_per_m)
print('normalized iso meso')
print(iso_meso/reads_per_m)
print('done')


#with gzip.open('C:/Splicer/Validate/notebooks_GTEX-1117F-0426-SM-5EGHI.splice.pickle.gz', 'rb') as handle:    
#    splices=pickle.load(handle)
    
#with gzip.open('C:/Splicer/Validate/notebooks_GTEX-1117F-0426-SM-5EGHI.exon_counts.chr1.pickle.gz', 'rb') as handle:    
#chr1_ex=cPickle.load(zf.open('chr1'))
    
i = 1   
        
#with gzip.open('C:/Splicer/Validate/notebooks_GTEX-1117F-0426-SM-5EGHI.splice.pickle.gz', 'rb') as handle:    
#    splices=pickle.load(handle)
    
#with gzip.open('C:/Splicer/Validate/notebooks_GTEX-1117F-0426-SM-5EGHI.exon_counts.chr1.pickle.gz', 'rb') as handle:    
#chr1_ex=cPickle.load(zf.open('chr1'))
    
  
#zf.close()
    


