import sqlite3
import os
import bamnostic as bn
import argparse
import bisect
import time

parser = argparse.ArgumentParser(description='takes the given .bam file and looks through all the reads to construct a count of all exons and splices in the reference splice graphs', usage='geneModelImport database_directory  GTF|UCSC GTF_File|knownGene.txt kgTxInfo.txt kgXRef.txt')
parser.add_argument("Database", help='The path to where you want to store the database file.')
parser.add_argument("Bam", help='The .bam file to count the reads from')
parser.add_argument("sampleName", help='GTF input file or UCSC input files.')
args = parser.parse_args()

#get connection to the sqlite database
conn = sqlite3.connect(args.Database + os.path.sep + 'splice.sqlite', isolation_level=None)
c = conn.cursor()


c.execute('''CREATE TABLE IF NOT EXISTS Sample
             (Sample_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Sample_Name varchar(40) NOT NULL DEFAULT NULL,
              Total_Reads INTEGER NOT NULL DEFAULT NULL, 
              Transcript_Reads INTEGER NOT NULL DEFAULT NULL);''')
c.execute("CREATE INDEX IF NOT EXISTS idx_sample_name ON Sample(sample_name);")

c.execute('''CREATE TABLE IF NOT EXISTS Sample_Attribs
             (Sample_ID INTEGER NOT NULL DEFAULT NULL,
              Attribute varchar(255) NOT NULL DEFAULT NULL,
              Value varchar(255) NOT NULL DEFAULT NULL);''')

c.execute('''CREATE TABLE IF NOT EXISTS Exon_Counts
             (Sample_ID INTEGER  NOT NULL DEFAULT NULL,
              SG_Exon_ID INTEGER NOT NULL DEFAULT NULL,
              Count INTEGER NOT NULL DEFAULT NULL);''')
c.execute("CREATE INDEX IF NOT EXISTS idx_ec_sample_id ON Exon_Counts(sample_id);")
c.execute("CREATE INDEX IF NOT EXISTS idx_ec_exon_id ON Exon_Counts(SG_Exon_ID);")

c.execute('''CREATE TABLE IF NOT EXISTS Splice_Counts
             (Sample_ID INTEGER  NOT NULL DEFAULT NULL,
              SG_Splice_ID INTEGER NOT NULL DEFAULT NULL,
              Count INTEGER NOT NULL DEFAULT NULL);''')
c.execute("CREATE INDEX IF NOT EXISTS idx_sc_sample_id ON Splice_Counts(sample_id);")
c.execute("CREATE INDEX IF NOT EXISTS idx_sc_splice_id ON Splice_Counts(SG_Splice_ID);")


#Find out the next assignable ID for this sample
c.execute("Select MAX(Sample_ID) FROM Sample")
ret = c.fetchone()
prevId = ret[0]
if prevId:
    Sample_Id = int(prevId)+1
else:
    Sample_Id = 1


#initialize the splice dictionary
sDict = {}
eDict = {}
epDict = {}
ecDict = {}
scDict = {}

start_time = time.time()

c.execute("SELECT SG_Splice_ID, Start_Position, Stop_Position, Chromosome FROM SG_Splice")
ret = c.fetchall()
#load the splice dictionary keyed on start-stop with the sg id as the value
for y in range(len(ret)):
    key = str(ret[y][1])+'-'+str(ret[y][2])
    chrom = ret[y][3]
    if not chrom.startswith("chr"):
        chrom = "chr"+chrom
    if chrom == "chrMT":
        chrom = "chrM"
    if chrom not in sDict:
        sDict[chrom] = {}
    sDict[chrom][key] = ret[y][0]
    
    
    
c.execute("SELECT SG_Exon_ID, Start_Position, Stop_Position, Chromosome FROM SG_Exon")
ret = c.fetchall()
#load the exon dictionary keyed on the start and the stop with the sg id as the value
for y in range(len(ret)):
    chrom = ret[y][3]
    if not chrom.startswith("chr"):
        chrom = "chr"+chrom
    if chrom == "chrMT":
        chrom = "chrM"
    if chrom not in eDict:
        eDict[chrom] = {}
        epDict[chrom] = []
    
    #add start    
    eDict[chrom][ret[y][1]] = ret[y][0]
    #add stop
    eDict[chrom][ret[y][2]] = ret[y][0]
    #add to tuple exon positions list (flip start and stop to correct if negative strand)
    if ret[y][1] < ret[y][2]:
        epDict[chrom].append((ret[y][1], ret[y][2]))
    else:
        epDict[chrom].append((ret[y][2], ret[y][1]))

#sorted list of all exon start stop tuples keyed on chromosome
for key in epDict:
    epDict[key] = sorted(epDict[key])

#"hg19test.bam"
samfile = bn.AlignmentFile(args.Bam, "rb")

def exonIncrement(start, stop, chro):
    exonIds = []
    try:
        pList = epDict[chro]
        #flip start in stop to correct for negative strand
        if start > stop:
            temp = start
            start = stop
            stop = temp
        #find the index that start belongs at
        idx = bisect.bisect(pList, (start,stop))
        
        i = idx
        if i == len(pList):
            return([])
        #move down the exons adding the ids of those included in the read
        while (i > -1  and start <= pList[i][1]):
            exonIds.append(eDict[chro][pList[i][0]])
            i-=1
        #ISSUE IF NEVER LOOP****************
        #if it goes off the end of a known exon add none and scrap the read
        if start < pList[i+1][0]:
            return([])
        
        i = idx
        looped = False
        #move up the exons adding ids of those included in the read
        while (i < len(pList) and i > -1 and stop >= pList[i][0]):
            exonIds.append(eDict[chro][pList[i][1]])
            i+=1
            looped = True
        #if it goes of the end of a known exon add none and scrap the read
        if looped and stop > pList[i-1][1]:
            return([])
    except:
        fe.write(chro+'\n')
    
    return(exonIds)

f = open('diagnostic.txt', 'w')    
fe = open('exceptions.txt', 'w')
i = 0
totalCount = 0
tranCount = 0
totalDupeCount = 0
tranDupeCount = 0
totalJRcount = 0
tranJRcount = 0
totalNJRcount = 0
tranNJRcount = 0
exceptionCount = 0
prevRead = ""
prevExons = ""
prevSplices = ""
for read in samfile:
    dupeTag = False
    exonSet = set()
    spliceSet = set()
    readR_S = read.reference_start
    readR_E = read.reference_end
    i+=1
    totalCount += 1
    #print(f"{totalCount:,d}")
    tranBool = False
    cigarString = read.cigarstring
    cigar = read.cigar
    chro = read.reference_name
    if str(readR_S)+"-"+str(readR_E)+"-"+cigarString+"-"+chro == prevRead:
        dupeTag = True
        totalDupeCount += 1
        for exon in prevExons:
            tranBool = True
            ecDict[exon] += 1
        for splice in prevSplices:
            tranBool = True
            scDict[splice] += 1
        if tranBool == True:
            tranDupeCount += 1
    elif "N" in cigarString:
        totalJRcount += 1
        #initialize the start and stop based on the first junction 
        start = readR_S+cigar[0][1]
        stop = start+cigar[1][1]+1
        #exon check from the start of the read to the start of the first splice
        exonSet.update(exonIncrement(readR_S+1, start, chro))
        for x in range(int(len(cigar)/2)):
            #if this is not the first junction adjust the start and stop
            if x != 0:
                start = stop+cigar[x*2][1]-1
                #exon check from the end of the last splice to the beginning of this splice 
                exonSet.update(exonIncrement(stop, start, chro))
                stop = start+cigar[x*2+1][1]+1
            
            #check if the splice is known and count it if so
            try:
                if str(start)+"-"+str(stop) in sDict[chro]:
                    spliceID = sDict[chro][str(start)+"-"+str(stop)]
                    spliceSet.add(spliceID)
                    tranBool = True
                    if spliceID in scDict:
                        scDict[spliceID] += 1
                    else:
                        scDict[spliceID] = 1
            except Exception as e:
                exceptionCount += 1
                fe.write(str(type(e))+': '+str(e)+"\n")
            exonID = ""
            
        exonSet.update(exonIncrement(stop, readR_E, chro))
        if (tranBool or len(exonSet) != 0):
            tranJRcount += 1        
    else:
        totalNJRcount += 1
        start = readR_S+1
        stop = start+cigar[0][1]
        exonSet.update(exonIncrement(start, stop, chro))
        if (len(exonSet) != 0):
            tranNJRcount += 1
            
    
    #add in all the sets
    for exon in exonSet:
        tranBool = True
        #print("exon: "+str(exon))
        if exon in ecDict:
            ecDict[exon] += 1
        else:
            ecDict[exon] = 1
            
    if tranBool == True:
        tranCount += 1
    elif dupeTag == False:
        f.write("read start:"+str(readR_S+1)+"  read stop:"+str(readR_E)+ "  negative strand?: "+ str(read.is_reverse)+"  "+chro+"\n")
    
    #set this line to prevRead
    if dupeTag == False:
        prevRead = str(readR_S)+"-"+str(readR_E)+"-"+cigarString+"-"+chro
        prevExons = exonSet
        prevSplices = spliceSet
    
    #if i == 5000000:
     #   break
     
c.execute('begin')
for key in scDict:
    c.execute("INSERT INTO Splice_Counts VALUES("+str(Sample_Id)+", "+str(key)+", "+str(scDict[key])+")")

for key in ecDict:
    c.execute("INSERT INTO Exon_Counts VALUES("+str(Sample_Id)+", "+str(key)+", "+str(ecDict[key])+")")

#add this sample to the sample table    
c.execute("INSERT INTO Sample VALUES("+str(Sample_Id)+", "+args.sampleName+", "+str(totalCount)+", "+str(tranCount)+")")
c.execute('commit')
f.close()
fe.close()
print("transcript junction reads: "+str(tranJRcount))
print("total junction reads: "+str(totalJRcount))
print("transcript duplicate reads: "+str(tranDupeCount))
print("total duplicate reads: "+str(totalDupeCount))
print("transcript Non junction reads: "+str(tranNJRcount))
print("total Non junction reads: "+str(totalNJRcount))

print("number of exceptions caught: "+ str(exceptionCount))

print("--- %.2f seconds---" % (time.time() - start_time))
print('Done')

