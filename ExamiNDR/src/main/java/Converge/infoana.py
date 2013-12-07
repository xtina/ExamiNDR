
import sys
import glob
import re
import string
from TAMO import MotifTools
import ConvergeMotifTools

'''
[10-28-02] Noticed problem that not all "similar" sequences are removed
           from info files.  Example: 9ant

[11-20-02] Moved WeightMask from MotifTools to this file
'''


def main():
    if len(sys.argv) < 2:
        print "Usage: %s <info_file> [mask_threshold (single digit)]"%(re.sub('^.*/','',sys.argv[0]))
        sys.exit(1)
    infofile   = sys.argv[1]
    maskthresh = 0
    if len(sys.argv) > 2:
        maskthresh = int(sys.argv[2])
    TAMO(infofile,maskthresh) # "Then a miracle occurs..."

################################################################################
# BSite                                                                        #
# Represent a binding site from transfac or DNA footprint                      #
#                                                                              #
# VAR: type - TF (Transfac) or FP (Footprint)                                  #
#      qual - 0 (FP) or 1-6 (Transfac) or 9 (? or YTFD or SCPD)                #
#      seq  - sequence of binding site                                         #
#                                                                              #
################################################################################
class BSite:
    '''Class to store transfac (or other) binding site information'''
    def __init__(self,text):
        self.type = ''
        self.qual = 6
        self.seq  = ''
        (info,self.seq) = text.split(':')
        if info[0:2] == 'TF':
            self.type = 'TF'
            if info[2] == '?':
                self.qual = 9
            else:
                self.qual = int(info[2])
        elif info[0:2] == 'FP':
            self.type = 'FP'
            self.qual = 0
    def cleantxt(self):
        if self.type == 'FP':
            return(re.sub('/.+','',self.seq).upper())
        else:
            _t = re.sub('[actg]','N',self.seq)
            _t = re.sub('^N+|N+$','',_t)
            return(_t)
    def __repr__(self):
        return "%s%s:%s"%(self.type,self.qual,self.seq)

################################################################################
# Pcp:                                                                         #
# Store a PCP string                                                           #
#                                                                              #
# VAR: seq  - PCP sequence                                                     #
#                                                                              #
################################################################################
class Pcp:
    '''Class to store a PCP string in PCP.seq'''
    def __init__(self,text):
        self.seq = text
    def __repr__(self):
        return(self.seq)

################################################################################
# Infofile                                                                     #
# Load and process information from a ".info" file                             #
#                                                                              #
# VAR: query      - DICT info of SITE                                          #
#      query_name - STR  full name of query                                    #
#      query_syns - LIST synonyms of query_name                                #
#      weights    - STR  % Representation of PCPs among structure (1-10, "X"=10#
#      neighbors  - DICT of DICTS information of sequences w similar PCPs      #
#                                                                              #
################################################################################
class Infofile:
    '''
    Class to load and manipulate an ".info" file, containing
    specificity and PCP proximity information to a query sequence
    '''
    def __init__(self,filename,no_rm_querys=''):
        self.query      = {}
        self.query_name = ''
        self.query_syns = []
        self.weights    = ''
        self.maskthresh = 0
        self.masktxt    = ''
        self.neighbors  = {}
        self._load_file(filename)
        self._query_syns()
        if not no_rm_querys:
            self._rm_querys()

    def __repr_dict(self,D):
        s = ''
        sites = ''
        i = 0
        for site in D['bsites']:
            i = i + 1
            if (i%4 == 0): sites = sites + '\n\t\t\t\t\t'
            sites = sites + " " + site.__repr__()
        s = s + "%3d %s \t%s" %(D['pcnt'],D['pcp'],sites)
        return s
    def __repr__(self):
        s = ''
        s = s + "# %s\n"%self.query_name
        if self.query_syns:
            s = s + "# Synonyms: "
            for syn in self.query_syns:
                s = s + "%s "%syn
            s = s + "\n"
        if self.query:
            sites =  " ".join(map(lambda x: x.__repr__(),self.query['bsites']))
            s = s + "#Q %-22s %s\n"%(self.query_name, self.__repr_dict(self.query))
        if self.weights:
            s = s + "#M %22s     %s\n"%("Weights",self.weights)
            s = s + "#M %22s     %s\n"%("Used for comparison:",self.masktxt)
        if self.neighbors:
            for Nkey in self.keys():
                N = self.neighbors[Nkey]
                s = s + "   %-22s %s\n" %(Nkey,self.__repr_dict(N))
        return(s)
               
    def keys(self):
        Nkeys = self.neighbors.keys()
        Nkeys.sort(lambda a,b,N=self.neighbors: cmp(N[b]['pcnt'],N[a]['pcnt']))
        return(Nkeys)

    def has_key(self,key):
        return(self.neighbors.has_key(key))
        
    def __getitem__(self,key):
        return self.neighbors[key]
                
    def _rm_querys(self):
        '''
        Remove any entry that might be the same as the query.  I
        wouldnt want to cheat by using the answer when making
        the prediction!
        '''
        if (len(self.query_syns) == 0):
            self._query_syns()
        #First remove Query line (to "self.query")
        if self.neighbors.has_key(self.query_name):
            self.query = self.neighbors[self.query_name]
            del(self.neighbors[self.query_name])
        for Nkey in self.neighbors.keys():
            for query_syn in self.query_syns:
                if Nkey.find(query_syn) >= 0:   #Similar to query
                    del self.neighbors[Nkey]
                    break

    def _query_syns(self):
        '''What names is the query sequence known by?'''
        self.query_syns.append(self.query_name)
        if (self.query_name.find(':') > 1):
            names = self.query_name.split(':')
            for name in names:
                self.query_syns.append(name)
                if re.search('_\d+$',name):
                    truncated_name = re.sub('_\d+$','',name)
                    self.query_syns.append(truncated_name)
    def _load_file(self,filename):
        '''Load the file into the "Infofile" data structure'''
        FH = open(filename,'r')
        for line in FH.readlines():
            toks = line.split()
            if toks[1] == 'Mask':
                self.weights = toks[-1]
                continue
            if (not self.query):
                self.query_name = toks[0]
            D = {}
            D['pcnt']  = int(toks[3])
            D['pcp']   = Pcp(toks[4])
            D['bsites'] = []
            for bs in toks[5:]:
                D['bsites'].append(BSite(bs))
            self.neighbors[toks[2]] = D
        FH.close()
    def fasta_by_key(self, filename, Nkey):
        if not self.has_key(Nkey):
            print "Error - No such key in info file: %s"%Nkey
            return
        FH = open(filename,'w')
        N = self[Nkey]
        i = 1
        for site in N['bsites']:
            if (i<10): i_txt = '0' + `i`
            else:      i_txt = `i`
            i = i + 1
            if site.seq.find('/') >= 0:
                seq = re.sub('/\w+','',site.seq)
            else:
                seq = re.sub('^[actg]+|[actg]+$','',site.seq)
            FH.write(">%s_%s:%s %3d\n%s\n"%(i_txt,Nkey,site.__repr__()[0:3].replace(':',''),
                                            N['pcnt'],seq))
        
    def bsites2seqs(self,thresh=100):
        seqs = []
        for Nkey in self.keys():
            N = self[Nkey]
            if N['pcnt'] >= thresh:
                for site in N['bsites']:
                    seqs.append(site.cleantxt())
        return(seqs)
                

    def fasta_threshold(self,filename,thresh=100):
        FH = open(filename,'w')
        for Nkey in self.keys():
            N = self[Nkey]
            if N['pcnt'] >= thresh:
                i = 1
                for site in N['bsites']:
                    if (i<10): i_txt = '0' + `i`
                    else:      i_txt = `i`
                    i = i + 1
                    if site.seq.find('/') >= 0:
                        seq = re.sub('/\w+','',site.seq)
                    else:
                        seq = re.sub('^[actg]+|[actg]+$','',site.seq)
                    FH.write(">%s_%s:%s %3d\n%s\n"%(i_txt,Nkey,site.__repr__()[0:3].replace(':',''),
                                                 N['pcnt'],seq))
    def _masktxt(self):
        maskA = []
        for w in self.weights:
            if ((w.isdigit() and int(w) > self.maskthresh) or w == 'X'):
                maskA.append('X')
            else:
                maskA.append('-')
        self.masktxt = ''.join(maskA)

    def remask(self,threshold=0):
        self.maskthresh = threshold
        self._masktxt()
        for Nkey in self.keys():
            N = self[Nkey]
            count = 0
            tot   = 0
            for i in range(len(N['pcp'].seq)):
                w = self.weights[i]
                if ((w.isdigit() and int(w) > threshold) or w == 'X'):
                    tot = tot + 1
                    if N['pcp'].seq[i] == self.query['pcp'].seq[i]:
                        count = count + 1
            N['pcnt'] = count * 100 / tot
        
################################################################################
# WeightMask                                                                   #
# Mask "hot areas" of probes by occurenct of top ranking n-mers                #
#                                                                              #
# ________________________________11122334456665544443311111__________________ #
# ATCGATGCTGATGCGTATAGCTATATGCGATATCCGCGCGCGCGCCCCGCCGCTCTCGATATCGCGATCTCGATAG #
#                                                                              #
# Then output these "hot" regions.  Intended to escape prior knowledge of      #
# motif width .                                                                #
#                                                                              #
################################################################################
class WeightMask:
    '''
    WeightMask                                                                   
    Mask "hot areas" of probes by occurenct of top ranking n-mers                
    
    ________________________________11122334456665544443311111__________________
    ATCGATGCTGATGCGTATAGCTATATGCGATATCCGCGCGCGCGCCCCGCCGCTCTCGATATCGCGATCTCGATAG
    '''
    def __init__(self, seqs, nmer_size=5):
        self.seqs      = seqs
        self.vecs      = []
        self.nmer_size = nmer_size
        self.nmersP    = []
        self.redseqs   = []
        self.mseqs     = []
        self._do_weights()
        self._reduce()

    def revcomp(self,seq):
        _t = map(lambda x,D=revcomp: D[x], seq)
        _t.reverse()
        return ''.join(_t)

    def _do_weights(self):
        self.nmersP = top_nmers(self.nmer_size,self.seqs,1) # 1 for with counts
        for seq in self.seqs:
            self.vecs.append(map(lambda x:0, seq))
        for (nmer,weight) in self.nmersP:
            nmerRC = self.revcomp(nmer)
            if weight < 2: continue
            for (seq,vec) in zip(self.seqs,self.vecs):
                for offset in range(len(seq)-self.nmer_size+1):
                    subseq = seq[offset:offset+len(nmer)]
                    if subseq == nmer or subseq == nmerRC:
                        for j in range(offset,offset+self.nmer_size):
                            vec[j] = vec[j] + weight

    def _reduce(self):
        redseqs = []
        mseqs   = []
        for (seq,vec) in zip(self.seqs,self.vecs):
            redseq = []
            maximum = max(vec)
            for (letter,count) in zip(seq,vec):
                if count < maximum/2:
                    redseq.append(' ')
                else:
                    redseq.append(letter)
            #print "%s\n%s\n%s\n"%(seq,''.join(redseq),''.join(redseq).split())
            redseqs.extend(''.join(redseq).split())
        redseqs.sort(lambda x,y: cmp(len(y),len(x)))
        self.redseqs = redseqs[:]
        for (nmer,count) in self.nmersP:
            done = 0
            for mseq in mseqs:          #Have we already covered this nmer?
                if mseq.find(nmer) >= 0:
                    #print "%s: Already found %s"%(nmer,mseq)
                    done = 1
            if not done:
                for redseq in self.redseqs:
                    if redseq.find(nmer) >= 0:
                        mseqs.append(redseq)
                        #print "%s\t%3d  %s"%(nmer,count,redseq)
                        break
        self.mseqs = mseqs

            

    def print_weights(self):
        for i in range(4):
            (s,c) = self.nmersP[i]
            print "%s\t%d"%(s,c)
        for i in range(len(self.seqs)):
            seq = self.seqs[i]
            vec = self.vecs[i]
            s = map(lambda x: '%-3s'%x, seq)
            v = map(lambda x: '%-3d'%x, vec)
            print ' '.join(s)
            print ' '.join(v)
            print


def ReduceInfo2seqs(Info,Threshold=50, C_FCN = lambda L: ConvergeMotifTools.compress_seqs_AA(L,50)):
    mseqs = []
    for key in Info.keys():
        if Info[key]['pcnt'] < Threshold:
            continue
        print "Working on %s :"%key
        cleanseqs = map(lambda x:x.cleantxt(), Info[key]['bsites'])
        if not cleanseqs:
            print "no sequences for %s"%key
            continue
        if len(Info[key]['bsites']) > 3:
            newseqs = C_FCN(cleanseqs)
        else:
            newseqs = cleanseqs
        mseqs.extend(newseqs)
        print "Adding %-22s %3d%% (%s)"%(key,Info[key]['pcnt'],' '.join(newseqs))
    return(mseqs)



def applyWeightMask(seqs):
    W = WeightMask(seqs,5)
    #W.print_weights()
    return(W.mseqs)

def Reduce_WeightMask(Info):
    print 'Computing weightmasks'
    mseqs  = ReduceInfo2seqs(Info, 80, applyWeightMask)
    mmseqs = applyWeightMask(mseqs[0:2])
    print mmseqs

def Reduce_Nmers(Info):
    print 'COMPUTING Nmers ....'
    mseqs = ReduceInfo2seqs(Info,70, lambda L: MotifTools.top_nmers(6,L)[0:3])
    print "Combining representative sequences...: "
    for i in range(len(mseqs)):
        i = i + 1
        print '\t%s'%mseqs[i-1],
        if (i%5 == 0): print
    print 

    top_seq_pairs = MotifTools.top_nmers(5,mseqs,1)
    total_nmers = 0
    for (mner,count) in top_seq_pairs:
        total_nmers = total_nmers + count
    for (nmer,count) in top_seq_pairs[0:8]:
        print "RESULT: %s\t%2d (%5.2f%%) occurences:  "%(nmer,count,
                                                         100*float(count)/total_nmers),
        for bsite in Info.query['bsites']:
            seq = bsite.cleantxt()
            (max,s1,s2) = MotifTools.compare_seqs(nmer,seq)
            print '   %s vs %s %4.2f correct'%(s1,s2,max)


def Reduce_AA(Info):
    mseqs  = ReduceInfo2seqs(Info,100)
    print "Combining representative sequences...: "
    for i in range(len(mseqs)):
        i = i + 1
        print '\t%s'%mseqs[i-1]
        if (i%3 == 0): print
    print 

        
    motifs = ConvergeMotifTools.compress_seqs_AA(mseqs, 50, 1) #1 Returns motif objects

    for motif in motifs:
        print "RESULT: %s\t%5.1f"%( motif,motif.MAP)
        for bsite in Info.query['bsites']:
            seq = bsite.cleantxt()
            #seq = 'NNNNNN%sNNNNNN'%seq
            (matches,ends,scores) = motif.scan(seq)
            for i in range(len(matches)):
                print "Good match: %s (%f)"%(matches[i],scores[i])


def TAMO(infofile,maskthreshold=0):
    '''
    Then a miracle occurs
    '''
    I = Infofile(infofile)
    I.remask(maskthreshold)
    print I
    #Reduce_WeightMask(I)
    Reduce_Nmers(I)


    #Reduce_AA(I)
    #reduce(I)
    #I.fasta_threshold('t.fsa',60)


if (__name__ == '__main__'): main()
