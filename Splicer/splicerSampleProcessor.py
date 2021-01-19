import sqlite3
import os
import bamnostic_mod as bn
import argparse
import bisect
import time

#Downloads\Lung\47b982b3-c7ce-4ca7-8c86-c71c15979620\G28588.NCI-H1915.1.bam
#Downloads\Lung\98a0206b-29f5-42d3-957b-6480e2fde185\G20483.HCC-15.2.bam
#Downloads\Lung\18004fb1-89a2-4ba1-a321-a0aa854e98c3\G25210.NCI-H510.1.bam
#Downloads\Lung\47030e40-acbd-4905-939c-d57441d7042e\G25222.NCI-H2171.1.bam
#Downloads\Lung\1357785f-f84b-4688-9b4c-0c2b5472ef51\G27281.RERF-LC-MS.1.bam
#Downloads\Lung\e48ea2ee-1dda-4061-a199-6e22fd2df382\G25212.NCI-H661.1.bam
#Downloads\Lung\f03dbfee-a523-438f-8459-f47f2ff1880f\G25224.NCI-H2066.1.bam

#Downloads\HeadAndNeck\0e67231f-97be-447c-b3b0-a656fc30a62d\G27454.PE_CA-PJ15.2.bam
#Downloads\HeadAndNeck\1acf65a0-0268-4288-9904-33bff618a31d\G27515.PE_CA-PJ41__clone_D2_.2.bam
#Downloads\HeadAndNeck\1f290458-df28-4c78-b73d-0202fb53bb0e\G27220.SCC-4.1.bam
#Downloads\HeadAndNeck\2b507086-977b-4cb7-abd9-83ee4ce9a893\G27489.PE_CA-PJ34__clone_C12_.2.bam
#Downloads\HeadAndNeck\7ed3e895-6826-430d-a39d-338111f16083\G27512.SNU-1214.2.bam
#Downloads\HeadAndNeck\c11aa745-72ea-44ca-b70d-7811c2f244b7\G27533.SNU-1066.2.bam
#Downloads\HeadAndNeck\dc8393c0-7d9e-4040-a91a-5783544cac35\G28853.HSC-4.3.bam

parser = argparse.ArgumentParser(description='takes the given .bam file and looks through all the reads to construct a count of all exons and splices in the reference splice graphs', usage='splicerSampleProcessor database_directory bam_file sample_name novelSplicesToggle(True|False)')
parser.add_argument("Database", help='The path to where you want to store the database file.')
parser.add_argument("Bam", help='The .bam file to count the reads from')
parser.add_argument("sampleName", help='Name for current sample in the sample table')
parser.add_argument("novelSplices", choices=['True', 'False'], help='Controls whether the program tries to find new splices')
args = parser.parse_args()

os.chdir("D:\\")
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

novelSplices = args.novelSplices

#initialize the splice dictionary
sDict = {}
eDict = {}
epDict = {}
ecDict = {}
scDict = {}
discoverySplices = {}

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
        if novelSplices == 'True':
            discoverySplices[chrom] = {}
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
    except Exception:
        pass
    
    return(exonIds)

fns = open ('novelSplices.txt', 'w')
i = 0

totalCount = 0
missingAttrCount = 0
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
    #read does not have an alignment
    if (
        not hasattr(read, 'reference_name')  or read.reference_name == None  or
        not hasattr(read, 'reference_start') or read.reference_start == None or
        not hasattr(read, 'reference_end')   or read.reference_end == None   or
        not hasattr(read, 'cigarstring')     or read.cigarstring == None     or
        not hasattr(read, 'cigar')           or read.cigar == None
       ):
        missingAttrCount += 1
        continue
    
    dupeTag = False
    exonSet = set()
    spliceSet = set()
    readR_S = read.reference_start
    readR_E = read.reference_end
    i+=1
    totalCount += 1
    if totalCount % 1000000 == 0:
        print(f"{totalCount:,d}")
        break
    tranBool = False
    cigarString = read.cigarstring
    cigar = read.cigar
    chro = read.reference_name
    if not chro.startswith("chr"):
        chro = "chr"+chro
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
                if str(start)+"-"+str(stop) in sDict[chro] or str(stop)+"-"+str(start) in sDict[chro]:
                    spliceID = sDict[chro][str(start)+"-"+str(stop)]
                    spliceSet.add(spliceID)
                    tranBool = True
                    if spliceID in scDict:
                        scDict[spliceID] += 1
                    else:
                        scDict[spliceID] = 1
                elif novelSplices == 'True':
                    if start in eDict[chro] and stop in eDict[chro]:
                        if str(start)+"-"+str(stop) in discoverySplices[chro]:
                            discoverySplices[chro][str(start)+"-"+str(stop)]+=1
                        else:
                            discoverySplices[chro][str(start)+"-"+str(stop)]=1
                            
                        experiSplicect = 1
            except Exception as e:
                exceptionCount += 1
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
        else:
            print("Missing: " + chro + ' ' + str(start) + ' ' + str(stop))    
    
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
for chromkey in discoverySplices:
    for skey in discoverySplices[chromkey]:
        fns.write(skey + ", Count: " + str(discoverySplices[chromkey][skey])+'\n')
fns.close()
print("missing attribute reads: " + str(missingAttrCount))
print("transcript junction reads: "+str(tranJRcount))
print("total junction reads: "+str(totalJRcount))
print("transcript duplicate reads: "+str(tranDupeCount))
print("total duplicate reads: "+str(totalDupeCount))
print("transcript Non junction reads: "+str(tranNJRcount))
print("total Non junction reads: "+str(totalNJRcount))

print("number of exceptions caught: "+ str(exceptionCount))

print("--- %.2f seconds---" % (time.time() - start_time))
print('Done')

