
''' This code generates Set cover problem patterns and overlap ratio is calculated 
w.r.t transactions as items. Prune transactions based on minRF'''

import time
import sys
import numpy as np
import copy
import itertools
from bitarray import bitarray
from PIL import Image

tot_CPs=[]
freqList=[]
#freqEle=[]
class cmine():
    def __init__(self, minRF, minCS, maxOR, inpfile, outfile,pathStr,datasetname):
        self.minRF = minRF
        self.minCS = minCS
        self.maxOR = maxOR
        self.inpfile = inpfile
        self.outfile = outfile
        self.nofs = self.getlines(inpfile)
        #self.fout = open(outfile,'w')
        self.noofCTP=0
        self.noof_Candi_CTP=0
        self.NOk = []
        self.Candi_CTP=[]
        self.final_CTPs=[]
        self.nofs,self.bitpattern = self.dbscanSCP(inpfile)
        #[self.items, self.bitpattern,self.TidKey_Dict] = self.dbscan(inpfile)
        
        "Checking feature coverage of a graph greaterthan minFC or not"
        print("No.of 1 size candidates:",len(self.bitpattern))        
        for key, value in self.bitpattern.items():
            self.noof_Candi_CTP=self.noof_Candi_CTP+1
            if (1.0*value.count()/self.nofs) >= minCS:
                tt=[[key],1.0*value.count()/self.nofs]
                tt.append(0)
                tot_CPs.append(tt)
            elif (1.0*value.count()/self.nofs) >= minRF:
                temp=[1.0*value.count()/self.nofs,key]
                self.NOk.append(temp)
                
        print("NO.of 1-size candidate patterns, 1-size GCP=",self.noof_Candi_CTP,len(tot_CPs))
        #self.NOk.sort(reverse=True)
        sorteditems = sorted(self.NOk, key = lambda a: (-a[0],a[1]))
        #print(len(sorteditems))        
        self.NOk=[]
        for i in sorteditems:
            self.NOk.append([i[1]])
            
            
            
        print("1-size nonOverlap transactions",len(self.NOk))
        freqItemcnt=[]
        for key, value in self.bitpattern.items():
            #print(key,value)
            freqItemcnt.append(value.count())
        #for ii in freqItemcnt:
         #   print(ii)    
        one_size_coverage=[]       
        
      
        
    def get_overlapratio_cs(self, new_tp):  
        cov_set=self.nofs*bitarray('0')
        last_tra_cov=self.nofs*bitarray('0')
        for tid in new_tp[:-1]:
            cov_set = cov_set | self.bitpattern[tid]
        #for tid in new_tp[-1]:
        last_tra_cov = self.bitpattern[new_tp[-1]]
            
        cov_sup = 1.0*(cov_set | last_tra_cov).count()/self.nofs
        
        ov_ratio= 1.0*(cov_set & last_tra_cov).count()/(last_tra_cov.count())
        return cov_sup,ov_ratio
   


    def expand(self,pathStr,datasetname):
        cnt = 0
        cnt1 = 0
        length = 1;
        while len(self.NOk)>0:
            #print("length",length,len(self.NOk))
            l_size_GCP=0
            l_size_Candi_pat=0
            temp_NOk = self.NOk
            self.NOk = []
            #print(temp_NOk)
            for i in range(len(temp_NOk)):
                for j in range(i+1, len(temp_NOk)):
                    cnt += 1
                    if temp_NOk[i][:-1] == temp_NOk[j][:-1]:
                        cnt1 += 1
                        newpattern = temp_NOk[i] + [temp_NOk[j][-1]]
                        cs,ov_ra=self.get_overlapratio_cs(newpattern)
                        l_size_Candi_pat=l_size_Candi_pat+1
                        self.noof_Candi_CTP=self.noof_Candi_CTP+1
                        if ov_ra <= maxOR:                            
                            if cs >= minCS:
                                self.noofCTP=self.noofCTP+1
                                temp=[newpattern,cs,ov_ra]
                                l_size_GCP=l_size_GCP+1
                                tot_CPs.append(temp)
                            else:
                                self.NOk.append(newpattern)
                            
                                
                    else:
                        break
            length += 1
            print(length, "length Candi patterns and GCP=",l_size_Candi_pat,l_size_GCP)
            print("pattern length=",length,len(self.NOk))
#        self.fout.close()
        return self.noof_Candi_CTP,self.noofCTP
    
   
    def getlines(self,filename):
        with open(filename,"r") as f:
            return sum(1 for _ in f)
            
    def dbscanSCP(self,db):
        # x=[]  
        noofbits = 0
        Tran_feat_Dict={}
        fidlist=[]
        g= open(db, 'r')
        for row in g:
            #if row[0]=='x':
            row = row.rstrip('\n')
            r = row.split(' ')
            for fid in r:
                if fid not in fidlist and fid != '':
                    fidlist.append(fid)
        self.nofs=len(fidlist)
        print("No.of features",len(fidlist))            
        f = open(db, 'r')
        Tidbitpat ={}        
        fid=0
        count=0
        for row in f:
            #if row[0]=='x':
            row = row.rstrip('\n')
            r = row.split(' ')
            bit_tra=self.nofs*bitarray('0')
            for fid in r:
                if fid !='':
                    bit_tra[int(fid)] = 'True'
            Tran_feat_Dict[count]=bit_tra
                   
            count=count+1
        return self.nofs,Tran_feat_Dict
              
              

start_time = time.time()
minRF = float(sys.argv[1])
minCS = float(sys.argv[2])
maxOR = float(sys.argv[3])
datasetname = sys.argv[4]
writePatterns = sys.argv[5]
pathStr='./Dataset/'
inpfile = pathStr+datasetname+".txt"
outfile = pathStr+str(datasetname)+"SetCover_Results.txt"
obj=cmine(minRF, minCS, maxOR, inpfile, outfile,pathStr,datasetname)
candidate_patterns,CTPs = obj.expand(pathStr,datasetname)
CTP_time = time.time()
print(str(datasetname)+", TC="+str(minRF)+", TPC="+str(minCS)+",overlap ratio"+str(maxOR)+", Exex Time="+str(CTP_time-start_time)+", No.of Candidate Patterns="+str(candidate_patterns)+",No.of GTCP="+str(len(tot_CPs)))
outStr=str(datasetname)+", TC="+str(minRF)+", TPC="+str(minCS)+",overlap ratio"+str(maxOR)+", Exex Time="+str(CTP_time-start_time)+", No.of Candidate Patterns="+str(candidate_patterns)+",No.of GTCP="+str(len(tot_CPs))
#outStr=str(datasetname)+","+str(minRF)+","+str(minCS)+","+str(maxOR)+","+str(len(tot_CPs))+","+str(candidate_patterns)+","+str(CTP_time-start_time)+"\n"
with open(str(pathStr)+"/"+str(datasetname)+"_Results.txt", 'a', newline = '') as q:
    q.write(str(outStr)+"\n")    
    q.close()      
if writePatterns=='True':
    with open(str(pathStr)+"/"+str(datasetname)+"_GTCPs.txt", 'a', newline = '') as q:
        for pat in tot_CPs:
            q.write(str(pat)+"\n")    
        q.close()

