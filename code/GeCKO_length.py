#! /usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import datetime
starttime = datetime.datetime.now()

# length seq: from 15 to 25
length_seq=np.zeros(21)
# forward begin bp: from 65 to 115
forward_begin=np.zeros(61)
# forward_end bp: from 65 to 115
forward_end = np.zeros(61)
#reverse_begin bp: from 65 to 115
reverse_begin = np.zeros(61)
#reverse_end bp: from 65 to 115
reverse_end = np.zeros(61)

line_num=0
file = open("data_aligned_primer.sam")
file.readline()
file.readline()
file.readline()
file.readline()
line=file.readline()

f=open("record2.txt","w+")

while line !="EOF":
    if len(line)<=1:
        break
    #total line

    line_num=line_num+1
    print(line_num)
    #f.writelines(str(line_num))
    #f.writelines("\n")
    #preprocess
    newline = line.split("\t")
    #print(newline)
    flag =newline[1]
    if flag=='4':
        line = file.readline()
        continue

    seqid1 = newline[0]
    print(seqid1)
    #f.writelines(seqid1 + "\n")
    line2 = file.readline()
    newline2 = line2.split("\t")
    seqid2 = newline2[0]
    #print(seqid2)
    if seqid1 == seqid2:
        print("match!")
        #f.writelines("match!" + "\n")
        refseq = newline[2]
        flag2 = newline[8]
        cigar = newline[5]
        if len(flag2)>1:
            seq = flag2[1:142]
        else: seq=newline[9]


        #print(len(newseq))
        print(cigar)
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


        cigar2 = newline2[5]
        print(cigar2)
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
        #print(finseq)
        #print(len(finseq))
        if len(finseq)>141 or len(finseq)<140:
            line = file.readline()
            continue

        f.writelines(seqid1 + "\n")
        f.writelines(finseq + "\n")

        #
        k=64
        sum_=0
        while sum_<6 and k<125:
            if finseq[k]=='X':
                sum_=sum_+1
            else:
                sum_=0
            k=k+1

        while k<125 and finseq[k]=='X':
            sum_=sum_+1
            k=k+1

        #print(sum_)
        end=k-1
        begin=end-sum_+1
        #print(end)
        #print(begin)

        length_seq_pos=sum_-15
        if 0<=length_seq_pos<=20:
            length_seq[length_seq_pos]=length_seq[length_seq_pos]+1
            f.writelines("right"+"\n")
        else:
            f.writelines("wrong"+"\n")
        #print(length_seq_pos)

        begin_pos=begin-65+1
        #print(begin_pos)
        #seq end_pos next
        end_pos=end-65+2
        #print(end_pos)

        if refseq.find("P5")==0:
            if 0 <= begin_pos <= 60 and 0 <= end_pos <= 60:
                forward_begin[begin_pos]=forward_begin[begin_pos]+1
                forward_end[end_pos]=forward_end[end_pos]+1
            #print(forward_begin)
            #print(forward_end)
        elif refseq.find("P3")==0:
            if 0<=begin_pos<=60 and 0<=end_pos<=60:
                reverse_begin[begin_pos]=reverse_begin[begin_pos]+1
                reverse_end[end_pos]=reverse_end[end_pos]+1
            #print(reverse_begin)
            #print(reverse_end)

        line = file.readline()

    else:
        line = line2
file.close()

"""
print(length_seq)
print(forward_begin)
print(forward_end)
print(reverse_begin)
print(reverse_end)
"""


#save to txt
np.savetxt('length_seq',length_seq)
np.savetxt('forward_begin',forward_begin)
np.savetxt('forward_end',forward_end)
np.savetxt('reverse_begin',reverse_begin)
np.savetxt('reverse_end',reverse_end)


print("Done!")
endtime = datetime.datetime.now()
print(endtime - starttime)
#f.writelines("Done!" + "\n")
#f.close()
