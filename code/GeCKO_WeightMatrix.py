#! /usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import datetime
starttime = datetime.datetime.now()

#forward primer matrix
m1 = np.zeros((4,141),dtype='uint64')
#reverse primer matrix
m2 = np.zeros((4,141),dtype='uint64')

line_num=0
file = open("data_aligned_primer.sam")
file.readline()
file.readline()
file.readline()
file.readline()
line=file.readline()

f=open("record.txt","w+")

while line !="EOF":
    if len(line)<=1:
        break
    #total line

    line_num=line_num+1
    print(line_num)
    f.writelines(str(line_num))
    f.writelines("\n")
    #preprocess
    newline = line.split("\t")
    #print(newline)
    flag =newline[1]
    if flag=='4':
        line = file.readline()
        continue

    seqid1 = newline[0]
    print(seqid1)
    f.writelines(seqid1 + "\n")
    line2 = file.readline()
    newline2 = line2.split("\t")
    seqid2 = newline2[0]
    #print(seqid2)
    if seqid1 == seqid2:
        print("match!")
        f.writelines("match!" + "\n")
        refseq = newline[2]
        flag2 = newline[8]
        cigar = newline[5]
        if len(flag2)>1:
            seq = flag2[1:142]
        else: seq=newline[9]


        #print(len(newseq))

        newseq = list(seq)
        #print(newseq)
        i=0
        j=0
        while len(cigar)>0:
            if '0'<=cigar[i]<='9':
                if i<len(cigar)-1:
                    i=i+1
                continue
            if cigar[i]=='M' or cigar[i]=='D':
                j=j+int(cigar[0:i])
                cigar=cigar[i+1:]
                #print(cigar)
                i=0
            elif cigar[i]=='I' or cigar[i]=='S':
                if int(cigar[0:i])>40:
                    break
                temp=1
                while temp<=int(cigar[0:i]) and j<len(newseq):
                    #print(j)
                    newseq[j]='X'
                    j = j + 1
                    temp=temp+1
                cigar=cigar[i+1:]
                #print(cigar)
                i=0

        #if len(flag2)>1:
        #    newseq = newseq[0:141]
        #print(newseq)
        #read next



        cigar2 = newline2[5]

        i=0
        l=140
        while len(cigar2)>0:
            if '0'<=cigar2[i]<='9':
                if i<len(cigar2)-1:
                    i=i+1
                continue
            if cigar2[i]=='M' or cigar2[i]=='D':
                l=l-int(cigar2[0:i])
                cigar2=cigar2[i+1:]
                #print(cigar)
                i=0
            elif cigar2[i]=='I' or cigar2[i]=='S':
                if int(cigar2[0:i])>40:
                    break
                temp=1
                while temp<=int(cigar2[0:i]) and l<len(newseq):
                    #print(j)
                    newseq[l]='X'
                    l = l - 1
                    temp=temp+1
                cigar2=cigar2[i+1:]
                #print(cigar)
                i=0

        while j<=l and j<len(newseq):
            newseq[j]='X'
            j=j+1

        #from list to string
        finseq = ''.join(newseq)
        print(finseq)
        if len(finseq)>141 or len(finseq)<140:
            line = file.readline()
            continue

        f.writelines(seqid1 + "\n")
        f.writelines(finseq + "\n")

        k=0
        if refseq.find("P5")==0:
           while k<141 and k<len(finseq):
               if finseq[k]=='A':
                 m1[0][k]=m1[0][k]+1
               elif finseq[k]=='C':
                 m1[1][k]=m1[1][k]+1
               elif finseq[k]=='G':
                 m1[2][k]=m1[2][k]+1
               elif finseq[k]=='T':
                 m1[3][k]=m1[3][k]+1
               k=k+1
        elif refseq.find("P3")==0:
           while k<141 and k<len(finseq):
               if finseq[k]=='A':
                 m2[0][k]=m2[0][k]+1
               elif finseq[k]=='C':
                 m2[1][k]=m2[1][k]+1
               elif finseq[k]=='G':
                 m2[2][k]=m2[2][k]+1
               elif finseq[k]=='T':
                 m2[3][k]=m2[3][k]+1
               k=k+1

        line = file.readline()

    else:
        line = line2
file.close()

#print(m1)
#print(m2)

#save to txt
np.savetxt('forward_primer',m1)
np.savetxt('reverse_primer',m2)

print("Done!")

f.writelines("Done!" + "\n")
f.close()

endtime = datetime.datetime.now()
print(endtime - starttime)
