import sqlite3
import os
import argparse

parser = argparse.ArgumentParser(description='Process initial database tables into refined SG tables.')
parser.add_argument("Database", help='The path to where you want to store the database file')
args = parser.parse_args()

#get connection to the sqlite database
conn = sqlite3.connect(args.Database + os.path.sep + 'splice.sqlite', isolation_level=None)
c = conn.cursor()

c.execute("DROP TABLE IF EXISTS SG_Exon;")
c.execute('''CREATE TABLE SG_Exon
             (Gene_ID INTEGER NOT NULL DEFAULT NULL,
              SG_Exon_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Exon_Number varchar(30) NOT NULL DEFAULT NULL,
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL);''')

c.execute("DROP TABLE IF EXISTS SG_Splice;")
c.execute('''CREATE TABLE SG_Splice
             (Gene_ID INTEGER NOT NULL DEFAULT NULL,
              SG_Splice_ID INTEGER PRIMARY KEY NOT NULL DEFAULT NULL,
              Chromosome varchar(30) NOT NULL DEFAULT NULL,
              Start_Position INTEGER NOT NULL DEFAULT NULL,
              Stop_Position INTEGER NOT NULL DEFAULT NULL);''')

#little object to hold on to the exons that are buing built and processed in exonBuilder
class SubExon:
    def __init__(self, name, start, stop):
        self.name = name
        self.start = start
        self.stop = stop

#builds and names exons from all the starts and stops in one gene
#then in populates SG_Exon table with these new exons    
def exonBuilder(starts, stops, strand, Gene_ID, chromosome, SG_Exon_ID):
    #add all the positions to a list and sort them
    posList = list(starts)
    posList.extend(stops)
    posList.sort()
    #if negative strand reverse
    if strand == "-":
        temp = starts
        starts = stops
        stops = temp
        posList.reverse()
    #boolean used to keep track of if the loop has seen a stop position
    boolStop = False
    current = SubExon("1.1", posList[0], None)
    exonCount = 1
    subExonCount = 1
    exonList = []
    #loop through the positions making appropriate subExon objects
    for i in range(1,len(posList)):
        
        #note the positions here are splice starts and stops so where a splice stops is where an exon starts
        if posList[i] in stops or boolStop == False:
            #set the stop marker if this is a stop position
            if posList[i] in stops:
                boolStop = True
            #update the current sub exon with its end position
            current.stop = posList[i]
            #if this is the first time the exon has been split change the name of the exon to have a .1
            subExonCount += 1
            #add the previous exon to the list before writing over the current variable
            exonList.append(current)
            #make the new sub exon object
            current = SubExon(str(exonCount)+"."+str(subExonCount), posList[i], None)
            
        else:
            boolStop = False
            #if the exon has no subexons remove the '.1' from the name
            if current.name.split('.')[1] == "2":
                exonList[-1].name = exonList[-1].name.split(".")[0]
            exonCount += 1
            subExonCount = 1
            current = SubExon(str(exonCount)+"."+str(subExonCount), posList[i], None)
    if current.name.split('.')[1] == "2":
        exonList[-1].name = exonList[-1].name.split(".")[0]
        #for exon in exonList:
        #add to database
    for exon in exonList:
        c.execute("INSERT INTO SG_Exon VALUES("+str(Gene_ID)+", "+str(SG_Exon_ID)+", "+exon.name+", '"+chromosome+"', "+str(exon.start)+", "+str(exon.stop)+")")  
        print("Exon: "+str(SG_Exon_ID))
        SG_Exon_ID += 1 
    return SG_Exon_ID


c.execute("select count(DISTINCT Gene_ID) From Transcript")
ret = c.fetchone()
numGenes = ret[0]
SG_Splice_ID = 1
SG_Exon_ID = 1
c.execute("begin")
for i in range(numGenes):
    c.execute("SELECT Transcript.Gene_ID, Exon.Chromosome, Exon.Start_Position, Exon.Stop_Position, Exon.Transcript_ID, Gene.Strand FROM Exon Join Transcript on Transcript.Transcript_ID = Exon.Transcript_ID JOIN Gene ON Transcript.Gene_ID = Gene.Gene_ID  WHERE Transcript.Gene_ID = "+str(i+1))
    ret = c.fetchall()
    previousTrans = 0
    previousStop = 0
    starts = set()
    stops = set()
    print(ret[0][0])
    for y in range(len(ret)):
        starts.add(ret[y][2])
        stops.add(ret[y][3])
        if previousTrans == ret[y][4]:
            Gene_ID = ret[y][0]
            thisStart = ret[y][2]
            chromosome = ret[y][1]
            c.execute("INSERT INTO SG_Splice VALUES("+str(Gene_ID)+", "+str(SG_Splice_ID)+", '"+chromosome+"', "+str(previousStop)+", "+str(thisStart)+")")
            SG_Splice_ID += 1
        previousStop = ret[y][3]
        previousTrans = ret[y][4]
        
    #run exon builder on the starts and stops gathered from the original exon entries for this gene
    SG_Exon_ID = exonBuilder(starts, stops, ret[0][5], ret[0][0], ret[0][1], SG_Exon_ID)
c.execute("commit")
print("done")
    
