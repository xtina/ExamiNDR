import sys, re, os, math, string, tempfile, copy, pickle
from   random import random,shuffle
COPY = copy

#Wishlist: Better ambiguity codes (based on information content)
#          Extract flanking subsequences (covering multiple matches)
#          Print bona-fide logo

one2two = {  'W':'AT',    'M':'AC',   'R':'AG',
             'S':'CG',    'Y':'CT',   'K':'GT'}
two2one = { 'AT': 'W',   'AC': 'M',  'AG': 'R',
            'CG': 'S',   'CT': 'Y',  'GT': 'K'}
revcomp = { 'A':'T',      'T':'A',    'C':'G',   'G':'C',
            'W':'W',      'S':'S',    'K':'M',   'M':'K',
            'Y':'R',      'R':'Y',    'N':'N',
            'B':'N', 'D':'N', 'H':'N', 'V':'N', ' ':'N'}  #[12-11-02] Needs fixing

ACGT = list('ACGT')
YEAST_BG = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31} #Yeast Default


revcomplement_memo = {'A':'T'}
revcompTBL = string.maketrans("AGCTagctWSKMYRnN", "TCGAtcgaWSMKTYnN")
def revcomplement(seq):
    global revcomplement_memo
    try:
        rc = revcomplement_memo[seq]
    except KeyError:
        #_t = map(lambda x,D=revcomp: D[x], seq)
        #get = revcomp.get
        #_t = map(get, seq)
        _t = list(seq.translate(revcompTBL))
        _t.reverse()
        rc = ''.join(_t)
        revcomplement_memo[seq] = rc
        revcomplement_memo[rc]  = seq
    return(rc)


def Motif_from_ll(ll):
    m = Motif(None,None)
    m.compute_from_ll(ll)
    return(m)

def Motif_from_counts(countmat,beta=0.01,bg={'A':.25,'C':.25,'G':.25,'T':.25}):
    m = Motif('',bg)
    m.compute_from_counts(countmat,beta)
    return m

def Motif_from_text(text,beta=0.05,source='',bg=None):
    if not bg: bg={'A':.25,'C':.25,'G':.25,'T':.25}
    m = Motif('',bg)
    m.compute_from_text(text,beta)
    m.source = source
    return(m)

def copy(motif):
    a = COPY.deepcopy(motif)
    #a.__dict__ = motif.__dict__.copy()
    return a



class Motif:
    def __init__(self,list_of_seqs_or_text=[],backgroundD=None):
        self.MAP       = 0
        self.evalue    = None
        self.oneletter = ''
        self.nseqs     = 0
        self.counts    = []
        self.closematches = []
        self.width     = 0
        self.fracs     = []
        self.logP      = []
        self.ll        = []
        self.bits      = []
        self.totalbits = 0
        self.maxscore  = 0
        self.minscore  = 0
        self.min_kscore  = 0
        self.pvalue      = 1
        self.pvalue_rank = 1
        self.church      = 1
        self.church_rank = 1
        self.Cpvalue     = 1
        self.Cpvalue_rank= 1
        self.Cchurch     = 1
        self.Cchurch_rank= 1
        self.E_seq       = None
        self.frac        = None
        self.E_site      = None
        self.E_chi2      = None
        self.kellis      = None
        self.MNCP        = None
        self.ROC_auc     = None
        self.realpvalue  = None
        self.Cfrac       = None
        self.CRA         = None
        self.valid     = None
        self.seeddist  = 0
        self.seednum   = -1
        self.seedtxt   = None
        self.family    = None
        self.source    = None
        self.threshold = None
        self._bestseqs = None
        self.bgscale   = 1
        self.best_pvalue = None
        self.best_factor = None
        self.gamma     = None
        self.nbound    = 0
        self.matchids  = []
        self.overlap   = None
        self.cumP      = []
        self.numbound      = 0
        self.nummotif      = 0
        self.numboundmotif = 0
        self.thetas    = []
        if backgroundD:
            self.background = backgroundD
        else:
            self.background = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31} #Yeast Default

        if type(list_of_seqs_or_text) == type(''):
            self.seqs      = []
            text           = list_of_seqs_or_text
            self.compute_from_text(text)
        else:
            self.seqs      = list_of_seqs_or_text
        if self.seqs:
            self._parse_seqs(list_of_seqs_or_text)
            self._compute_ll()
            self._compute_oneletter()
            #self._compute_threshold(2.0)

    def __repr__(self):
        return "%s (%d)"%(self.oneletter, self.nseqs)
    def summary(self):
        m = self
        txt = "%-34s (Bits: %5.2f  MAP: %7.2f   D: %5.3f  %3d)  E: %7.3f"%(
            m, m.totalbits, m.MAP, m.seeddist, m.seednum, nlog10(m.pvalue))
        if m.church != None:  txt = txt + '  ch: %6.2f'%(nlog10(m.church))
        if m.frac   != None:  txt = txt + '  f: %5.3f'%(m.frac)
        if m.E_site != None:  txt = txt + '  Es: %6.2f'%(nlog10(m.E_site))
        if m.E_seq  != None:  txt = txt + '  Eq: %6.2f'%(nlog10(m.E_seq))
        if m.MNCP   != None:  txt = txt + '  mn: %6.2f'%(m.MNCP)
        if m.ROC_auc!= None:  txt = txt + '  Ra: %6.4f'%(m.ROC_auc)
        if m.E_chi2 != None:
            if m.E_chi2 == 0: m.E_chi2=1e-20
            txt = txt + ' x2: %5.2f'%(nlog10(m.E_chi2))
        if m.CRA    != None:  txt = txt + '  cR: %6.4f'%(m.CRA)
        if m.Cfrac  != None:  txt = txt + '  Cf: %5.3f'%(m.Cfrac)
        if m.realpvalue != None: txt = txt + '  P: %6.4e'%(m.realpvalue)
        if m.kellis != None:  txt = txt +  '  k: %6.2f'%(m.kellis)
        if m.numbound      :  txt = txt +  '  b: %3d'%(m.numbound)
        if m.nummotif      :  txt = txt +  '  nG: %3d'%(m.nummotif)
        if m.numboundmotif :  txt = txt +  '  bn: %3d'%(m.numboundmotif)

        return txt

    def minimal_raw_seqs(self):
        ''' Return minimal list of seqs that represent consensus '''
        seqs = [[], []]
        for letter in self.oneletter:
            if one2two.has_key(letter):
                seqs[0].append(one2two[letter][0])
                seqs[1].append(one2two[letter][1])
            else:
                seqs[0].append(letter)
                seqs[1].append(letter)
        if ''.join(seqs[0]) == ''.join(seqs[1]):
            return( [''.join(seqs[0])] )
        else:
            return( [''.join(seqs[0]), ''.join(seqs[0])] )
    def _compute_oneletter(self):
        letters = []
        for i in range(self.width):
            downcase = None
            if self.bits[i] < 0.25:
                letters.append('.')
                continue
            if self.bits[i] < 1.0: downcase = 'True'
            tups = [(self.ll[i][x],x) for x in ACGT if self.ll[i][x] > 0.0]
            if not tups:  #Kludge if all values are negative (can this really happen?)
                tups = [(self.ll[i][x],x) for x in ACGT]
                tups.sort()
                tups.reverse()
                tups = [tups[0]]
                downcase = 'True'
            tups.sort()      #Rank by LL
            tups.reverse()
            bases = [x[1] for x in tups[0:2]]
            bases.sort()
            if len(bases) == 2: L = two2one[''.join(bases)]
            else:               L = bases[0]
            if downcase: L = L.lower()
            letters.append(L)
        self.oneletter = ''.join(letters)
    def _compute_oneletter_old(self):
        letters = []
        for i in range(self.width):
            #if self.bits[i] < .2:
            #    letters.append('N')
            #    continue
            list = []
            keys = ['A', 'C', 'T', 'G']
            keys.sort(lambda x,y,D=self.ll[i]: cmp(D[y],D[x]))
            for key in keys:
                if self.ll[i][key] > 0.0:
                    list.append(key)
                else: break
            if len(list) == 1:
                letters.append(list[0])
            elif len(list) == 2:
                slist = list[0:2]
                slist.sort()
                pair = ''.join(slist)
                letters.append( two2one[pair] )
            else:
                letters.append('N')
        self.oneletter = ''.join(letters)

    def _parse_seqs(self, LOS):
        self.nseqs = len(LOS)
        self.width = len(LOS[0])
        for i in range(self.width):
            Dc = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N': 0}
            for seq in LOS:
                key = seq[i]
                Dc[key] = Dc[key] + 1
            del(Dc['N'])
            self.counts.append(Dc)

    def _compute_ll(self):
        self.fracs = []
        self.logP  = []
        self.ll    = []
        for i in range(self.width):
            Dll  = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            Df   = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            DlogP= {'A': 0, 'C': 0, 'T': 0, 'G': 0}
            for key in self.counts[i].keys():
                #print i,key,self.counts[i][key],self.nseqs
                Pij = self.counts[i][key] / float(self.nseqs)
                Df [key] = Pij
                Dll[key] = (math.log( (self.counts[i][key] + self.bgscale*self.background[key] ) / 
                                      ((self.nseqs + self.bgscale) * self.background[key])     ) /
                            math.log(2))
                if Pij > 0:
                    DlogP[key]  = math.log(Pij)/math.log(2)
                else:
                    DlogP[key]  = -100  #Near zero
            self.fracs.append(Df)
            self.logP.append (DlogP)
            self.ll.append   (Dll)
        self.P = self.fracs
        self._compute_bits()
        self._compute_ambig_ll()
        self._maxscore()


    def compute_from_ll(self,ll):
        self.ll    = ll
        self.closematches = []
        self.width = len(ll)
        self._compute_bg_from_ll()
        self._compute_logP_from_ll()
        self._compute_ambig_ll()
        self._compute_bits()
        self._compute_oneletter()
        self._maxscore()

    def _computeP(self):
        P = []
        for i in range(self.width):
            #print i,
            _p = {}
            for L in ACGT: _p[L] = math.pow(2.0,self.logP[i][L])
            P.append(_p)
        #print
        self.P = P

    def _compute_bits(self):
        bits = []
        totbits = 0
        bgbits  = 0
        bg      = self.background
        UNCERT  = lambda x: x*math.log(x)/math.log(2.0)
        for letter in ACGT:
            bgbits = bgbits + UNCERT(bg[letter])
        for i in range(self.width):
            tot = 0
            for letter in ACGT:
                Pij = pow(2.0, self.logP[i][letter])
                tot = tot + UNCERT(Pij)
                #bit = Pij * self.ll[i][letter]
                #if bit > 0:
                #    tot = tot + bit
            #print tot, bgbits, tot-bgbits
            bits.append(max(0,tot-bgbits))
            totbits = totbits + max(0,tot-bgbits)
            #added november 8th, 2004 by KM
            #bits.append(tot)
            #totbits = totbits+tot
        self.bits = bits
        self.totalbits = totbits

    def _cut_ll(self,val=-1.0):
        for Dll in self.ll:
            for key in Dll.keys():
                if Dll[key] < val:
                    Dll[key] = val
        
    def denoise(self,bitthresh=0.5):
        for i in range(self.width):
            tot = 0
            for letter in ACGT:
                if self.logP:
                    Pij = pow(2.0, self.logP[i][letter])
                else:
                    Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                if Pij > 0.01:
                    bit = Pij * self.ll[i][letter]
                    tot = tot + bit
            if tot < bitthresh:  #Zero Column
                for letter in ACGT:
                    self.ll[i][letter] = 0.0
        self.compute_from_ll(self.ll)

    def giflogo(self,id,title=None,scale=0.8,info_str=''):
        return giflogo(self,id,title,scale,info_str)

    def _print_bits(self,norm=2.3, height=8.0):
        bits   = []
        tots   = []
        str    = []
        for i in range(self.width):
            D = {}
            tot = 0
            for letter in ['A', 'C', 'T', 'G']:
                if self.logP:
                    Pij = pow(2.0, self.logP[i][letter])
                else:
                    Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                if Pij > 0.01:
                    '''Old'''
                    D[letter] = Pij * self.ll[i][letter]
                    #'''new'''
                    #Q = self.background[letter]
                    #D[letter] = ( Pij * math.log(Pij) - Pij * math.log(Q) ) / math.log(2.0)
                    '''for both old and new'''
                    tot = tot + D[letter]
            bits.append(D)
            tots.append(tot)
        for i in range(self.width):
            s = []
            _l = bits[i].keys()
            _l.sort(lambda x,y,D=bits[i]: cmp(D[y],D[x]))
            for key in _l:
                for j in range(int(bits[i][key] / norm * height)):
                    s.append(key)
            str.append(''.join(s))
        fmt = '%%%ds'%height
        print '#  %s'%('-'*self.width)
        for h in range(height):
            sys.stdout.write("#  ")
            for i in range(self.width):
                sys.stdout.write((fmt%str[i])[h])
            if h == 0:
                sys.stdout.write(' -- %4.2f bits\n'%norm)
            elif h == height-1:
                sys.stdout.write(' -- %4.2f bits\n'%(norm/height))
            else:
                sys.stdout.write('\n')
        print '#  %s'%('-'*self.width)
        print '#  %s'%self.oneletter

    def print_logo(self,norm=2.0, height=8.0):
        bits   = []
        tots   = []
        str    = []
        for i in range(self.width):
            D = {}
            tot = 0
            for letter in ['A', 'C', 'T', 'G']:
                if self.logP:
                    Pij = pow(2.0, self.logP[i][letter])
                else:
                    Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                if Pij > 0.01:
                    '''Old'''
                    D[letter] = Pij * self.bits[i]
                    #'''new'''
                    #Q = self.background[letter]
                    #D[letter] = ( Pij * math.log(Pij) - Pij * math.log(Q) ) / math.log(2.0)
                    '''for both old and new'''
                    tot = tot + D[letter]
            bits.append(D)
            tots.append(tot)
        for i in range(self.width):
            s = []
            _l = bits[i].keys()
            _l.sort(lambda x,y,D=bits[i]: cmp(D[y],D[x]))
            for key in _l:
                for j in range(int(bits[i][key] / norm * height)):
                    s.append(key)
            str.append(''.join(s))
        fmt = '%%%ds'%height
        print '#  %s'%('-'*self.width)
        for h in range(height):
            sys.stdout.write("#  ")
            for i in range(self.width):
                sys.stdout.write((fmt%str[i])[h])
            if h == 0:
                sys.stdout.write(' -- %4.2f bits\n'%norm)
            elif h == height-1:
                sys.stdout.write(' -- %4.2f bits\n'%(norm/height))
            else:
                sys.stdout.write('\n')
        print '#  %s'%('-'*self.width)
        print '#  %s'%self.oneletter
                    

    def _compute_ambig_ll(self):
        for Dll in self.ll:
            for L in one2two.keys():
                Dll[L] = max(Dll[one2two[L][0]],  Dll[one2two[L][1]] )
            Dll['N'] = 0.0
            Dll['B'] = 0.0

    def compute_from_nmer(self,nmer,beta=0.001):  #For reverse compatibility
        self.compute_from_text(nmer,beta)

    def compute_from_text(self,text,beta=0.001):
        prevlett = {'B':'A', 'D':'C', 'V':'T', 'H':'G'}
        countmat = []
        text = re.sub('[\.\-]','N',text.upper())

        for i in range(len(text)):
            D = {'A': 0, 'C': 0, 'T':0, 'G':0}
            letter = text[i]
            if letter in ['B', 'D', 'V', 'H']:  #B == no "A", etc...
                _omit = prevlett[letter]
                for L in ACGT:
                    if L != _omit: D[L] = 0.3333
            elif one2two.has_key(letter):  #Covers WSMYRK
                for L in list(one2two[letter]):
                    D[L] = 0.5
            elif letter == 'N':
                for L in D.keys():
                    D[L] = self.background[L]
            elif letter == '@':
                for L in D.keys():
                    D[L] = self.background[L]-(0.0001)
                D['A'] = D['A'] + 0.0004
            else:
                D[letter] = 1.0
            countmat.append(D)
        self.compute_from_counts(countmat,beta)

    def new_bg(self,bg):
        counts = []
        for pos in self.logP:
            D = {}
            for L,lp in pos.items():
                D[L] = math.pow(2.0,lp)
            counts.append(D)
        self.background = bg
        self.compute_from_counts(counts,0)
        
    def addpseudocounts(self,beta=0):
        self.compute_from_counts(self.counts,beta)
    
    def compute_from_counts(self,countmat,beta=0):
        self.counts  = countmat
        self.width   = len(countmat)
        self.bgscale = 0

        maxcount = 0
        #Determine Biggest column
        for col in countmat:
            tot = 0
            for v in col.values():
                tot = tot + v
            if tot > maxcount: maxcount = tot

        #Pad counts of remaining columns
        for col in countmat:
            tot = 0
            for c in col.values():
                tot = tot + c
            pad = maxcount - tot
            for L in col.keys():
                col[L] = col[L] + pad * self.background[L]
                
        self.nseqs = maxcount
        nseqs      = maxcount
        
        #Add pseudocounts        
        if beta > 0:  
            multfactor = {}
            bgprob = self.background
            pcounts= {}
            for L in bgprob.keys():
                pcounts[L] = beta*bgprob[L]*nseqs 
            for i in range(self.width):
                for L in countmat[i].keys():
                    _t = (countmat[i][L] + pcounts[L]) #Add pseudo
                    _t = _t / (1.0 + beta)    #Renomalize
                    countmat[i][L] = _t

        #Build Motif
        self.counts = countmat
        self._compute_ll()
        self._compute_oneletter()
        self._maxscore()


    def compute_from_counts_old(self,countmat,beta=0):
        self.counts = countmat
        self.width  = len(countmat)
        self.bgscale= 0
        
        _T = []
        _t = 0
        for i in range(self.width):
            _T.append(0)
            for L in countmat[i].keys():
                _T[i] = _T[i] + countmat[i][L]
                _t    = _t    + countmat[i][L]
        if (math.fabs( float(_t)/self.width - _T[0]) > 0.001):
            print "Problem with counts: %f != %f * %f\n%s"%(
                _t, _T[0], self.width, _T)

        self.nseqs = float(_T[0])
        nseqs = self.nseqs
        if beta > 0:  #Add pseudocounts
            multfactor = {}
            bgprob = self.background
            pcounts= {}
            for L in bgprob.keys():
                pcounts[L] = beta*bgprob[L]*nseqs 
            for i in range(self.width):
                for L in countmat[i].keys():
                    _t = (countmat[i][L] + pcounts[L]) #Add pseudo
                    _t = _t / (1.0 + beta)    #Renomalize
                    countmat[i][L] = _t
        self.counts = countmat
        self._compute_ll()
        self._compute_oneletter()
        self._maxscore()

    def _compute_bg_from_ll(self):
        '''Compute background model from log-likelihood matrix
        by noting that:   pA  + pT  + pC  + pG  = 1
                  and     bgA + bgT + bgC + bgG = 1
                  and     bgA = bgT,   bgC = bgG
                  and so  bgA = 0.5 - bgC
                  and     pA = lA / bgA,  etc for T, C, G
                  so...
                         (lA + lT)bgA + (lC + lG)bgC          =  1
                         (lA + lT)bgA + (lC + lG)(0.5 - bgA)  =  1
                         (lA + lT - lC - lG)bgA +(lC +lG)*0.5 =  1
                          bgA                                 =  {1 - 0.5(lC + lG)} / (lA + lT - lC - lG)
        + Gain accuracy by taking average of bgA over all positions of PSSM
        '''
        pow = math.pow
        bgATtot = 0
        nocount = 0
        near0   = lambda x:(-0.01 < x and x < 0.01)
        for i in range(self.width):
            _D = self.ll[i]
            ATtot = pow(2,_D['A']) + pow(2,_D['T'])
            GCtot = pow(2,_D['C']) + pow(2,_D['G'])
            if near0(_D['A']) and near0(_D['T']) and near0(_D['G']) and near0(_D['C']):
                nocount = nocount + 1
                continue
            if near0(ATtot-GCtot):     #Kludge to deal with indeterminate case
                nocount = nocount + 1
                continue
            bgAT   = (1.0 - 0.5*GCtot)/(ATtot - GCtot)
            if (bgAT < 0.1) or (bgAT > 1.1):
                nocount = nocount + 1
                continue
            bgATtot = bgATtot + bgAT
        if nocount == self.width:  #Kludge to deal with different indeterminate case
            self.background = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
            return
        bgAT = bgATtot / (self.width - nocount)
        bgGC = 0.5 - bgAT
        self.background = {'A':bgAT, 'C':bgGC, 'G':bgGC, 'T':bgAT}            
        
    def _compute_logP_from_ll(self):
        log = math.log
        logP = []
        for i in range(self.width):
            D = {}
            for L in ACGT:
                ''' if   ll = log(p/b) then
                       2^ll = p/b
                  and    ll = log(p) - log(b)
                  so log(p) = ll + log(b)'''
                #Pij = pow(2.0, self.ll[i][letter]) * self.background[letter]
                D[L] = self.ll[i][L] + log(self.background[L])/log(2.)
            logP.append(D)
        self.logP = logP

    def _print_ll(self):
        print "#  ",
        for i in range(self.width):
            print "  %4d   "%i,
        print
        for L in ['A', 'C', 'T', 'G']:
            print "#%s "%L,
            for i in range(self.width):
                print  "%8.3f "%self.ll[i][L],
            print
    def _print_p(self):
        print "#  ",
        for i in range(self.width):
            print "  %4d   "%i,
        print
        for L in ['A', 'C', 'T', 'G']:
            print "#%s "%L,
            for i in range(self.width):
                print  "%8.3f "%math.pow(2,self.logP[i][L]),
            print
    def _print_counts(self):
        print "#  ",
        for i in range(self.width):
            print "  %4d   "%i,
        print
        for L in ['A', 'C', 'T', 'G']:
            print "#%s "%L,
            for i in range(self.width):
                print  "%8.3f "%self.counts[i][L],
            print
        
    def _maxscore(self):
        total = 0
        lowtot= 0
        for lli in self.ll:
            total = total + max(lli.values())
            lowtot= lowtot+ min(lli.values())
        self.maxscore = total
        self.minscore = lowtot
        self.kscore()

    def kscore(self):
        KLs = []
        KLsum = 0.0
        kmin = 0.0    #best possible match
        for i in range(self.width):
            dvg = 0
            for L in 'ACGT':
                dvg = dvg + math.pow(2,self.logP[i][L])*self.ll[i][L]
            KLs.append(dvg)
        for logPi,KL in zip(self.logP,KLs):
            KLsum = KLsum + KL
            kmin = kmin + KL*(1 - math.pow(2,max(logPi.values())))
        if (kmin==0.0):
            self.min_kscore = 0.0
        else:
            self.min_kscore = kmin/KLsum
        
    def _compute_threshold(self,z=2.0):
        scoretally = []
        for seq in self.seqs:
            matches,endpoints,scores = self.scan(seq,-100)
            scoretally.append(scores[0])
        ave,std = avestd(scoretally)
        self.threshold = ave - z *std
        #print '#%s: threshold %5.2f = %5.2f - %4.1f * %5.2f'%(
        #    self, self.threshold, ave, z, std)

    def bestscanseq(self,seq):
        matches,endpoints,scores = self.scan(seq,-100)
        t = zip(scores,matches)
        t.sort()
        bestseq   = t[-1][1]
        bestscore = t[-1][0]
        return bestscore, bestseq
    
    def bestscan(self,seq):
        matches,endpoints,scores = self.scan(seq,-100)
        if not scores: return -100
        scores.sort()
        best = scores[-1]
        return best

    def matchstartorient(self,seq, factor=0.7):
        ans = []
        txts,endpoints,scores = self.scan(seq,factor=factor)
        for txt, startstop in zip(txts,endpoints):
            start, stop = startstop
            rctxt  = revcomplement(txt)
            orient = (self.bestscore(txt,1) >= self.bestscore(rctxt,1))
            ans.append((start,orient))
        return ans

    def scan(self, seq, threshold = '', factor=0.7):
        if len(seq) < self.width:
            return(self._scan_smaller(seq,threshold))
        else:
            return(self._scan(seq,threshold,factor=factor))

    def scansum(self,seq,threshold = -1000):
        ll = self.ll
        sum = 0
        width        = self.width
        width_r      = range(width)
        width_rcr    = range(width-1,-1,-1)
        width_ranges = zip(width_r,width_rcr)
        seqcomp      = seq.translate(revcompTBL)

        total = 0
        hits  = 0
        etotal= 0
        for offset in range(len(seq)-width+1):
            total_f = 0
            total_r = 0
            for i,ir in width_ranges:
                pos = offset+i
                total_f = total_f + ll[i][    seq[pos]]
                total_r = total_r + ll[i][seqcomp[pos]]
            total_max = max(total_f,total_r)
            if total_max >= threshold:
                total = total + total_max
                etotal = etotal + math.exp(total)
                hits  = hits + 1
            if not hits:
                ave = 0
            else:
                ave = float(total)/float(hits)
        return(total,hits,ave,math.log(etotal))
    def score(self, seq, fwd='Y'):
        matches, endpoints, scores = self._scan(seq,threshold=-100000,forw_only=fwd)
        return scores[0]
    def bestscore(self,seq, fwd=''):
        matches, endpoints, scores = self._scan(seq,threshold=-100000,forw_only=fwd)
        if scores: return max(scores)
        else:      return -1000

    def _scan(self, seq,threshold='',forw_only='',factor=0.7):
        ll = self.ll #Shortcut for Log-likelihood matrix
        if not threshold: threshold = factor * self.maxscore
        
        #print '%5.3f'%(threshold/self.maxscore)
        matches       = []
        endpoints     = []
        scores        = []
        width         = self.width
        width_r       = range(width)
        width_rcr     = range(width-1,-1,-1)
        width_ranges  = zip(width_r,width_rcr)

        seqcomp = seq.translate(revcompTBL)
        

        for offset in range(len(seq)-self.width+1):    #Check if +/-1 needed
            total_f = 0
            total_r = 0
            for i,ir in width_ranges:
                pos = offset+i
                total_f = total_f + ll[i ][    seq[pos]]
                total_r = total_r + ll[ir][seqcomp[pos]]

            if 0 and total_f > 1:
                for i,ir in width_ranges:
                    print seq[offset+i],'%6.3f'%ll[i ][        seq[offset+i] ],'   ',
                print '= %7.3f'%total_f
                
            if 0:
                print "\t\t%s vs %s: F=%6.2f R=%6.2f %6.2f %4.2f"%(seq[offset:offset+self.width],
                                                                   self.oneletter,total_f,total_r,
                                                                   self.maxscore,
                                                                   max([total_f,total_r])/self.maxscore)
            if total_f > threshold and ((total_f > total_r) or forw_only):
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_f)
                matches.append(seq[offset:offset+self.width])
            elif total_r > threshold:
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_r)
                matches.append(seq[offset:offset+self.width])
        return(matches,endpoints,scores)
    def _scan_smaller(self, seq, threshold=''):
        '''The sequence is smaller than the PSSM.  Are there
           good matches to regions of the PSSM?
        '''
        ll = self.ll #Shortcut for Log-likelihood matrix
        matches   = []
        endpoints = []
        scores    = []
        w         = self.width
        for offset in range(self.width-len(seq)+1):    #Check if +/-1 needed
            maximum = 0
            for i in range(len(seq)):
                maximum = maximum + max(ll[i+offset].values())
            if not threshold: threshold = 0.8 * maximum
            total_f = 0
            total_r = 0
            for i in range(len(seq)):
                total_f = total_f + ll[i+offset      ][        seq[i] ]
                total_r = total_r + ll[w-(i+offset)-1][revcomp[seq[i]]]
            if 0:
                print "\t\t%s vs %s: F=%6.2f R=%6.2f %6.2f %4.2f"%(seq, self.oneletter[offset:offset+len(seq)],
                                                                   total_f, total_r,  maximum,
                                                                   max([total_f,total_r])/self.maxscore)
            if total_f > threshold and total_f > total_r:
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_f)
                matches.append(seq[offset:offset+self.width])
            elif total_r > threshold:
                endpoints.append( (offset,offset+self.width-1) )
                scores.append(total_r)
                matches.append(seq[offset:offset+self.width])
        return(matches,endpoints,scores)                

    def mask_seq(self,seq):
        masked = ''
        matches, endpoints, scores = self.scan(seq)
        cursor = 0
        for start, stop in endpoints:
            masked = masked + seq[cursor:start] + 'N'*self.width
            cursor = stop+1
        masked = masked + seq[cursor:]
        return masked

    def masked_neighborhoods(self,seq,flanksize):
        ns = self.seq_neighborhoods(seq,flanksize)
        return [self.mask_seq(n) for n in ns]

    def seq_neighborhoods(self,seq,flanksize):
        subseqs = []
        matches, endpoints, scores = self.scan(seq)
        laststart, laststop = -1, -1
        for start, stop in endpoints:
            curstart, curstop = max(0,start-flanksize), min(stop+flanksize,len(seq))
            if curstart > laststop:
                if laststop != -1:
                    subseqs.append(seq[laststart:laststop])
                laststart, laststop = curstart, curstop
            else:
                laststop = curstop
        if endpoints: subseqs.append(seq[laststart:laststop])
        return subseqs
        

    def maskdiff(self,other):
        return maskdiff(self,other)

    def __sub__(self,other):
        '''Return Euclidean distance between this pssm and another'''
        if type(other) != type(self):
            print "computing distance of unlike pssms (types %s, %s)"%(
                type(other),type(self))
            print 'First: %s'%other
            print 'Self:  %s'%self
            sys.exit(1)
        if other.width != self.width:
            print "computing distance of unlike pssms (width %d != %d)"%(
                other.width,self.width)
            sys.exit(1)
        D = 0
        FABS = math.fabs
        POW  = math.pow
        for L in self.logP[0].keys():
            for i in range(self.width):
                D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
                #D = D + FABS( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]))
                #D = D + FABS(self.logP[i][L] - other.logP[i][L])
        return(math.sqrt(D))
    def maxdiff(self):
        '''Compute maximum possible Euclidean distance to another motif.'''
        POW  = math.pow
        D = 0
        for i in range(self.width):
            _min = 100
            _max = -100
            for L in ACGT:
                val = POW(2,self.logP[i][L])
                if   val > _max:
                    _max  = val
                    _maxL = L
                elif val < _min:
                    _min  = val
                    _minL = L
            for L in ACGT:
                if L == _minL:
                    delta = 1-POW(2,self.logP[i][L])           #1-val
                    D = D + delta*delta
                else:
                    D = D + POW( POW(2,self.logP[i][L]), 2)    #0-val
        return(math.sqrt(D))
                
    def revcomp(self):
        return revcompmotif(self)
    def trimmed(self,thresh=0.1):
        for start in range(0,self.width-1):
            if self.bits[start]>=thresh: break
        for stop  in range(self.width,1,-1):
            if self.bits[stop-1]>=thresh: break
        m = self[start,stop]
        return m
    def bestseqs(self,thresh=None):
        if not thresh:
            if self._bestseqs:
                return self._bestseqs
        if not thresh: thresh = 0.8 * self.maxscore
        self._bestseqs = bestseqs(self,thresh)
        return self._bestseqs
    def emit(self,prob_min=0.0,prob_max=1.0):
        if not self.cumP:
            for logcol in self.logP:
                tups = []
                for L in ACGT:
                    p = math.pow(2,logcol[L])
                    tups.append((p,L))
                tups.sort()
                cumu = []
                tot  = 0
                for p,L in tups:
                    tot = tot + p
                    cumu.append((tot,L))
                self.cumP.append(cumu)
        s = []
        #u = random()+0.01 #Can make higher for more consistent motifs
        u = (prob_max-prob_min)*random() + prob_min
        for cumu in self.cumP:
            #u = random()+0.01 #Can make higher for more consistent motifs
            last = 0
            for p,L in cumu:
                if last < u and u <= p:
                    letter = L
                    break
                else: last = p
#           print L,'%8.4f'%u,cumu
            s.append(L)
        #print ''.join(s)
        return ''.join(s)
            
                
    def random_kmer(self):
        if not self._bestseqs: self._bestseqs = self.bestseqs()
        seqs   = self._bestseqs
        pos = int(random() * len(seqs))
        print 'Random: ',self.oneletter,seqs[pos][1]
        return(seqs[pos][1])
    def __getitem__(self,tup):
        '''Interface to submotif.  Less pythonish, but more reliable'''
        if len(tup) != 2:
            print "Motif[i,j] requires two arguments, not ",tup
        else:
            beg, end = tup[0], tup[1]
            return(submotif(self,beg,end))
    def __getslice__(self,beg,end):
        if beg >= end:
            #Probably python converted negative idx.  Undo
            beg = beg - self.width
        return(submotif(self,beg,end))
    def __add__(self,other):
        return merge(self,other,0)
    def __len__(self):
        return(self.width)
    def shuffledP(self):
        return shuffledP(self)
    def copy(self):
        a = Motif()
        a.__dict__ = self.__dict__.copy()
        return a
    def random_diff_avestd(self,iters=5000):
        return(random_diff_avestd(self,iters))
    def bogus_kmers(self,count=200):
        '''For making logos with logo program'''
        POW  = math.pow
        #Build p-value inspired matrix
        #Make totals cummulative:
        # A: 0.1 C: 0.4 T:0.2 G:0.3
        #                            ->  A:0.0 C:0.1 T:0.5 G:0.7  0.0
        
        #Take bg into account:
        # We want to pick P' for each letter such that:
        #     P'/0.25  = P/Q
        # so  P'       = 0.25*P/Q
        
        m = []
        for i in range(self.width):
            _col = []
            tot   = 0.0
            for L in ACGT:
                _col.append( tot )
                tot = tot + POW(2,self.logP[i][L]) * 0.25 / self.background[L]
            _col.append(tot)
            #Renormalize
            for idx in range(len(_col)):
                _col[idx] = _col[idx] / _col[-1]
            m.append(_col)

        for p in range(0): #Was 5
            for i in range(len(m)):
                print '%6.4f  '%m[i][p],
            print

        seqs=[]
        for seqnum in range(count+1):
            f = float(seqnum)/(count+1)
            s = []
            for i in range(self.width):
                for j in range(4):
                    if (m[i][j] <= f and f < m[i][j+1]):
                        s.append(ACGT[j])
                        break
            seqs.append(''.join(s))

        del(seqs[0])
        #for i in range(count):
        #    print ">%3d\n%s"%(i,seqs[i])

        return(seqs)


def minwindowdiff(M1,M2,overlap=5,diffmethod='diff'):
    #Alternate method: maskdiff, infomaskdiff
    if type(M1) != type(M2):
        print "Error: Attempted to compute alignment of objects that are not both Motifs"
        print "       types %s: %s  and %s: %s"%(M1,type(M1),M2,type(M2))
        sys.exit(1)

    if M1.width <= M2.width: A = M1; Borig = M2
    else:                    A = M2; Borig = M1
    wA = A.width
    wB = Borig.width
    O  = overlap

    if   diffmethod == 'diff':
        diff_fcn = diff
    elif diffmethod == 'maskdiff':
        diff_fcn = maskdiff
    elif diffmethod == 'infomaskdiff':
        diff_fcn = infomaskdiff
        
    mindiff = 1000
    #print 'minwindodebug    wA ', wA, 'wB ', wB, 'O ', O, 'wA-0', wA-O, 'wB-O', wB-O
    for Astart in range(wA-O+1):
        subA = A[Astart:Astart+O]
        for B in [Borig, Borig.revcomp()]:
            for Bstart in range(wB-O+1):
                subB = B[Bstart:Bstart+O]
                mindiff = min(mindiff, diff_fcn(subA,subB))
                #print 'minwindodebug     ',subA, subB, diff_fcn(subA,subB)
    return(mindiff)
    

def minaligndiff(M1,M2,overlap=5,diffmethod='diff'):
    #Alternate method: maskdiff, infomaskdiff
    if type(M1) != type(M2):
        print "Error: Attempted to compute alignment of objects that are not both Motifs"
        print "       types %s: %s  and %s: %s"%(M1,type(M1),M2,type(M2))
        sys.exit(1)

    if M1.width <= M2.width:
        A = M1; Borig = M2
        switch = 0
    else:
        A = M2; Borig = M1
        switch = 1
    wA = A.width
    wB = Borig.width
    O  = overlap

    '''
    Here is the figure to imagine:
       012345678901234567890   wA: 6  Bstart: 6-3     = 3
         A         (A)         wB: 11 Bstop:  6+11-3-1= 13
       ------     %%%%%%        O: 3  lastA:  6+11-3-3= 11
          -----------
          |O|  B
    '''

    if   diffmethod == 'diff':
        diff_fcn = diff
    elif diffmethod == 'maskdiff':
        diff_fcn = maskdiff
    elif diffmethod == 'infomaskdiff':
        diff_fcn = infomaskdiff
    
    Bstart = wA-O
    Bstop  = wA+wB-O-1
    lastA  = wA+wB-O-O
    Dmin = 1000
    Dmins=[]
    #print A
    #print '%s%s'%(' '*Bstart,Borig)
    for B in [Borig, Borig.revcomp()]:
        for start in range(0,lastA+1):
            Bpos = []
            Apos = []
            for offset in range(wA):
                abs = start+offset
                if abs >= Bstart and abs <= Bstop:
                    Apos.append(offset)
                    Bpos.append(abs-Bstart)
            subA = A[min(Apos),max(Apos)+1]
            subB = B[min(Bpos),max(Bpos)+1]
            #print '%s%s\n%s%s  %f'%(
            #    ' '*start, subA,
            #    ' '*start, subB,   diff_fcn(subA,subB))
            if switch: _diff = diff_fcn(subB,subA)
            else:      _diff = diff_fcn(subA,subB)
            Dmin = min(Dmin, _diff)
    return(Dmin)
    
'''
To compare 2 motifs of the same width, there are these five functions:

m1 - m2            - Euclidean Distance (sqrt(sum_col(sum_row)))
diff(m1,m2)        - psuedo-Euclidean (sum_col(sqrt(norm(sum_row)))/#col
maskdiff(m1,m2)    - diff, but excluding positions with "N" in m2
infomaskdiff(m1,m2)- diff, but scaling distance by normalized
     information content at each position in m2.
diverge(m1,m2)     - Mutual information sum[p log (p/q)]

**Note that maskdiff, infomaskdiff, and diverge are not symmetric functions

To compare 2 motifs of different widths, there is the function:

minaligndiff(M1,M2,overlap=5,diffmethod='diff')

this does a 'sliding' comparison of two motifs and reports the minimum
distance over all alignments.  overlap refers to the minumum overlap
required while sliding.  Below, overlap is '2'.  The default is '5'.

      ------
          -----------

You can optionally specify the distance metric as a text string.
The default is 'diff'.

'''


def diff(self,other):
    '''Return pseudo-Euclidean distance'''
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)
    POW     = math.pow
    Dtot    = 0
    for i in range(self.width):
        '''Computes distance'''
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
        Dtot = Dtot + math.sqrt(D)/math.sqrt(2.0)
    return(Dtot/self.width)
    

def maskdiff(self,other):
    '''Return pseudo-Euclidean distance, but only include columns that are not background'''
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)

    Dtot = 0
    POW  = math.pow
    NEAR0= lambda x:(-0.01 < x and x < 0.01)
    divisor = 0
    for i in range(self.width):
        nearcount = 0

        '''Implements mask'''
        for L in ACGT:
            diff = POW(2,other.logP[i][L]) - other.background[L]
            if NEAR0(diff): nearcount = nearcount + 1
        if nearcount == 4:
            #print 'Skipping position %d :'%i,other.logP[i]
            continue

        '''Computes distance'''
        divisor = divisor + 1
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
        Dtot = Dtot + math.sqrt(D)/math.sqrt(2.0)
    return(Dtot/divisor)

def infomaskdiff(self,other):
    '''
    Return pseudo-Euclidean distance, but scale column distance by
    information content of "other".
    '''
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)

    maxbits = math.log( 1.0/min(other.background.values()) ) / math.log(2.0)
    '''or... alternatively'''
    #print maxbits, max(other.bits)
    #print other.bits
    maxbits = max(other.bits)
    if maxbits < 0.1:  #'''There is nothing important here'''
        return 1
    
    Dtot    = 0
    POW     = math.pow
    divisor = 0
    '''Computes distance'''
    for i in range(self.width):
        D = 0
        for L in ACGT:
            D = D + POW( POW(2,self.logP[i][L]) - POW(2,other.logP[i][L]), 2 )
        col_dist  = math.sqrt(D)/math.sqrt(2.0)
        col_scale = other.bits[i]/maxbits
        divisor = divisor + col_scale
        Dtot = Dtot + col_dist*col_scale
    return(Dtot/divisor)

def diverge(self,other):
    '''Return Euclidean distance, but only include columns that are not background'''
    if type(other) != type(self):
        print "computing distance of unlike pssms (types %s, %s)"%(
            type(other),type(self))
        print 'First: %s'%other
        print 'Self:  %s'%self
        sys.exit(1)
    if other.width != self.width:
        print "computing distance of unlike pssms (width %d != %d)"%(
            other.width,self.width)
        sys.exit(1)

    Dtot = 0
    POW  = math.pow
    LOG2 = lambda x:math.log(x)/math.log(2.0)
    NEAR0= lambda x:(-0.01 < x and x < 0.01)
    divisor = 0
    for i in range(self.width):
        nearcount = 0

        '''Implements mask'''
        for L in ACGT:
            diff = POW(2,other.logP[i][L]) - self.background[L]
            if NEAR0(diff): nearcount = nearcount + 1
        if nearcount == 4:
            #print 'Skipping position %d :'%i,other.logP[i]
            continue

        '''Computes distance'''
        divisor = divisor + 1
        D = 0
        for L in ACGT:
            Pself = POW(2, self.logP[i][L])
            Pother= POW(2,other.logP[i][L])
            D = D + Pself * LOG2(Pself/Pother)
        Dtot = Dtot + D
    return(Dtot/divisor)



def bestseqs(motif,thresh,
            seq='',score=0,depth=0,bestcomplete=None,SEQS=[]):
    '''
    This function returns all sequences that a motif could
    match match with a sum(log-odds) score greater than thresh
    Aborts at 2000 sequences
    '''
    if depth == 0:
        SEQS = []  #Must be a python bug. I shouldn't have to do this
    if not bestcomplete:
        M = motif
        maxs = []
        for i in range(M.width):
            bestj = 'A'
            for j in ['C', 'G', 'T']:
                if M.ll[i][j] > M.ll[i][bestj]:
                    bestj = j
            maxs.append(M.ll[i][bestj])
        bestcomplete = []
        for i in range(M.width):
            tot = 0
            for j in range(i,M.width):
                tot = tot + maxs[j]
            bestcomplete.append(tot)
    if depth == motif.width:
        if score > thresh:
            SEQS.append((score,seq))
        if len(SEQS) > 2000:
            thresh = 1000.0 # Return Early, You don't really want all these sequences, do you?
        return
    if depth==-1:
        print '# %-10s %6.3f %6.3f %2d'%(seq, score, bestcomplete[depth], depth)
    if score + bestcomplete[depth] < thresh: return
    if depth > 0 and len(SEQS) > 2000:
        return
    for L in ACGT:
        newseq   = seq + L
        newscore = score + motif.ll[depth][L]
        bestseqs(motif,thresh,newseq,newscore,depth+1,bestcomplete,SEQS)
    if depth == 0:
        SEQS.sort()
        SEQS.reverse()
        return SEQS

def seqs2fasta(seqs,fasta_file = ''):
    if not fasta_file:
        fasta_file = tempfile.mktemp()
    FH = open(fasta_file,'w')
    for i in range(len(seqs)):
        FH.write(">%d\n%s\n"%(i,seqs[i]))
    FH.close()
    return(fasta_file)

def compress_seqs_AA(seqs,AA_iters=10,returnMotifs=0):
    compressed_seqs = []
    fasta_file = seqs2fasta(seqs)
    MA = MetaAce(fasta_file,4,AA_iters)
    os.remove(fasta_file)
    for motif in MA.best(2):
        for seq in motif.minimal_raw_seqs():
            compressed_seqs.append(seq)
    if returnMotifs:
        return(MA.best(10))
    else:
        return(compressed_seqs)

def compress_seqs_nmer(width,seqs):
    return(top_nmers(width,seqs)[0:2])

def top_nmers(N,seqs,with_counts = 0,purge_Ns = '',gap = 0):
    Nmers = {}
    revcompTBL = string.maketrans("AGCTagctnN", "TCGAtcganN")
    if gap:
        flank_length = N/3
        gap_flank = N - (N/3)
    for seq in seqs:
        for i in range(len(seq)-N+1):
            if gap: Nmer = seq[i:i+flank_length] + seq[(i+gap_flank):i+N]
            else: Nmer = seq[i:i+N]
            if purge_Ns:
                if Nmer.find('N') >= 0: continue
            _t = list(Nmer.translate(revcompTBL))
            _t.reverse()
            NmerRC = ''.join(_t)   # _t used until here to revese comp seq
            _t = [Nmer, NmerRC]
            _t.sort()
            NmerKey = _t[0]        # _t used until here to get alphabetically first seq
            if Nmers.has_key(NmerKey):
                Nmers[NmerKey] = Nmers[NmerKey] + 1
            else:
                Nmers[NmerKey] = 1
    sorted = Nmers.keys()
    sorted.sort(lambda x,y,D=Nmers:cmp(D[y],D[x]) or cmp(x,y))
    #for i in range(10):
    #    print "# %2d  %s %d"%(i,sorted[i],Nmers[sorted[i]])
    if with_counts:
        return(zip(sorted,map(lambda x,N=Nmers:N[x], sorted)))
    else:
        return(sorted)

def m_matches(seqs,wmer,m):
    '''
    list = m_matches(seqs,wmer,m)
    m represents minimum number of identitical positions to wmer
    '''
    matches = []
    width = len(wmer)
    for (nmer, count) in top_nmers(width,seqs,'with counts'):
        match = 0
        for i in range(width):
            if nmer[i] == wmer[i]:
                match = match+1
        if match >= m:
            for i in range(count):
                matches.append(nmer)
    return(matches)

def compare_seqs(s1, s2):
    if len(s1) > len(s2):
        long  = s1
        short = s2
    else:
        long  = s2
        short = s1
    (maxcount,max_i) = (0,0)
    for i in range(len(long)-len(short)+1):
        idcount_f = 0
        idcount_r = 0
        for j in range(len(short)):
            if short[j] == long[i+j]:
                idcount_f = idcount_f + 1
            if short[-(j+1)] == revcomp[long[i+j]]:
                idcount_r = idcount_r + 1
        if (idcount_f > maxcount and idcount_f >= idcount_r):
            maxcount = idcount_f
            max_i    = i
        elif (idcount_r > maxcount):
            maxcount = idcount_r
            max_i    = i
        #print i,j,idcount_f,idcount_r,maxcount
    maxfrac = float(maxcount) / len(short)
    print maxfrac,maxcount,len(short)
    return(maxfrac,short,long[max_i:max_i+len(short)])

class WeightMask:
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


def random_diff_avestd(motif,iters=5000):
    w = motif.width
    vals = []
    for i in range(iters):
        vals.append(motif - Random_motif(w))
    return(avestd(vals))

def Random_motif(w):
    C = []
    for i in range(w):
        D = {}
        tot = 0
        p = int(random.random() * 4)
        Lup = ACGT[p]
        for L in ACGT:
            D[L] = 0.1
            tot = tot + 0.001
        D[Lup] = D[Lup] + 1
        for L in ACGT:
            D[L] = D[L]/tot
        C.append(D)
    m = Motif()
    m.compute_from_counts(C)
    return(m)

def toDict(M):
    '''Return a list of dictionaries, assuming M entries
       are in alphabetical order (ACGT)'''
    if type(M[0]) == type(0.0):
        return(toDictVect(M))
    else:
        a = []
        for i in range(len(M)):
            a.append(toDictVect(M[i]))
        return(a)
        
def toDictVect(V):
    '''Return a dictionarie, assuming V entries
       are in alphabetical order (ACGT)'''
    D = {}
    for L,i in (('A',0), ('C',1), ('G',2), ('T',3)):
        D[L]=V[i]
    return(D)

def main():
    if len(sys.argv) < 2:
        print "Usage: %s <fasta_file>"%(re.sub('^.*/','',sys.argv[0]))
        sys.exit(1)
    
    J = AlignAce(sys.argv[1])
    for motif in J.motifs:
        if motif.MAP > 10:
            print motif
    print J

                
def submotif(self,beg,end):
    bg = self.background.copy()
    P = []

    #Determine if any 'zeros' should be added at begining
    #because the user has specified a negative beg index
    for i in range(beg,0):
        P.append(bg.copy())

    #Copy relevant content of motif
    start = max(beg,0)
    stop  = min(end,self.width)
    for i in range(start,stop):
        D = {}
        for L in ACGT:
            D[L] = math.pow(2.,self.logP[i][L])
        P.append(D)

    #Determine if any 'zeros' should be added at the end
    #because the user has specified a width too large
    for i in range(self.width,end):
        P.append(bg.copy())

    #print "BEG, END", beg,end
    #for i in range(beg,end):
    #    print i,P[i]

    #Build the Motif
    M = COPY.deepcopy(self)
    #M = Motif(None,bg.copy())
    M.compute_from_counts(P)
    M.source = self.source
    return(M)
                
def shuffledP(self):
    bg = self.background.copy()
    P = []

    #Copy relevant content of motif
    for i in range(0,self.width):
        D = {}
        _s = ACGT[:]
        shuffle(_s)
        for L,_L in zip(ACGT,_s):
            D[L] = math.pow(2.,self.logP[i][_L])
        P.append(D)

    #Build the Motif
    M = COPY.deepcopy(self)
    #M = Motif(None,bg.copy())
    M.compute_from_counts(P)
    M.source = self.source
    return(M)

def revcompmotif(self):
    bg = self.background.copy()
    P = []

    for i in range(self.width):
        D = {}
        for L in ACGT:
            D[L] = math.pow(2.,self.logP[self.width-i-1][revcomp[L]])
        P.append(D)

    #Build the Motif
    M = COPY.deepcopy(self)
    M.compute_from_counts(P)
    return(M)
        

def sum(motifs,weights=[]):
    if not weights:
        weights = [1.0] * len(motifs)
    tot = 0.0
    for w in weights: tot=tot+float(w)
    weights = [(w/tot) for w in weights]
    C = []
    for c in motifs[0].fracs:
        D = {}
        for L in ACGT: D[L] = 0.0
        C.append(D)
    for m,w in zip(motifs,weights):
        for i in range(m.width):
            for L in ACGT:
                C[i][L] = C[i][L] + m.fracs[i][L]*w
    motif = Motif_from_counts(C,0.0,bg=motifs[0].background)
    return motif.trimmed()
                
'''
Bob- If this is an annoyance, dont go through the pain of porting it again.
'''
def giflogo(motif,id,title=None,scale=0.8,info_str=''):
    SEQLOGO = '/cluster/bgordon/biobin/seqlogo'
    kmers   = motif.bogus_kmers(100)
    width   = float(len(kmers[0]) )
    height  = float(4)
    m       = motif
    width, height = width*scale, height*scale
    tmp     = tempfile.mktemp() + '.fsa'
    if title==None: title = ' '

    '''
    if not info_str:
        info_str = "%s b: %5.2f  D: %5.3f s: %2d  E: %6.3f S: %6.3f"%(
            m.oneletter, m.totalbits, m.seeddist, m.seednum,
            -math.log(m.pvalue)/math.log(10.), -math.log(m.church)/math.log(10.))
        #print m, m.source, m.seedtxt, m.seednum
        if m.source:
            info_str = '%s  %s'%(m.source,info_str)
    '''
    if not info_str: info_str = ' '
    seqs2fasta(kmers,tmp)
    cmd = '%s -F GIF -acpY -w%d -h%d -k 1 -F GIF -o %s -M -f %s -t "%s" -Z "%s"'%(
        SEQLOGO, width, height, id, tmp, title, info_str)

    GIMP = '/usr/bin/gimp'
    convcmd = "%s -c -d -i -s -b \'(let* ((img  (car (file-gif-load 1 \"%s.gif\" \"%s.gif\" ))) "%(GIMP,id,id) + \
              ' (draw (car (gimp-image-active-drawable img)))) ' + \
              ' (file-gif-save 1 img draw  \"%s.gif\"  \"%s.gif\" 1 0 0 0))\' '%(id,id)  + \
              " \'(gimp-quit 0)\' < /dev/null > /dev/null 2>&1"

    #print convcmd 
    os.system(cmd)
    os.system(convcmd)
    os.unlink(tmp)
    return "%s.gif"%id

def merge(A,B,overlap=0):
    '''Merge motifs A and B into a new motif, with a specified overlap'''
    if (overlap < 0 or overlap > A.width or overlap >B.width):
        print 'Cannot overlap %s with %s by %d bases'%(A.oneletter,B.oneletter,overlap)
        return(None)

    #Build Probability matrix.  Width will be A.width + B.width - overlap
    w = A.width + B.width - overlap

    P = []
    #Make a copy of A's probabilities into P
    for i in range(A.width):
        D = {}
        logP = A.logP[i]
        for L in logP.keys():
            D[L] = math.pow(2,logP[L])
        P.append(D)
    #Add B's first 'overlap' probabilities to last 'overlap' probabilities of P
    for i in range(overlap):
        logP = B.logP[i]
        Pidx = len(P)-overlap+i
        _tot = 0
        for L in logP.keys():
            P[Pidx][L] = (P[Pidx][L] + math.pow(2,logP[L])) / 2.0
            P[Pidx][L] = max(P[Pidx][L],math.pow(2,logP[L]))
            _tot = _tot + P[Pidx][L]
        for L in logP.keys():
            P[Pidx][L] = P[Pidx][L] / _tot
    #Append B's remaining probabilites to P
    for i in range(overlap,B.width):
        D = {}
        logP = B.logP[i]
        for L in logP.keys():
            D[L] = math.pow(2,logP[L])
        P.append(D)
        
    #Build a motif
    M = Motif(None,A.background.copy())
    M.source = A.source,B.source
    M.compute_from_counts(P)
    return(M)

def avestd(vals):
    (sum, sum2) = (0.,0.)
    N = float(len(vals))
    for val in vals:
        sum  = sum  + float(val)
        sum2 = sum2 + float(val)*float(val)
    if N == 1:
        ave = sum
        std = 0
    else:
        ave = sum /  N
        std = math.sqrt( (sum2-(N*ave*ave)) / (N-1.0) )
    return(ave,std)


def load(filename):
    FID = open(filename,'r')
    lines = FID.readlines()
    FID.close()
    motifs   = []
    seedD    = {}
    seedfile = ''
    for i in range(len(lines)):
        if lines[i][0:10] == 'Log-odds matrix'[0:10]:
            w = len(lines[i+1].split())-1
            ll = []
            for pos in range(w):
                ll.append({})
            for j in range(0,4):
                toks = lines[i+j+2].split()
                L = toks[0][1]
                for pos in range(w):
                    ll[pos][L] = float(toks[pos+1])
            m = Motif_from_ll(ll)
            motifs.append(m)
        if lines[i][0:6] == 'Motif '[0:6]:
            toks =  lines[i].split()
            motifs[-1].nseqs    = float(re.sub('[\(\)]','',toks[3]))
            motifs[-1].totalbits= float(toks[5])
            motifs[-1].MAP      = float(toks[7])
            motifs[-1].seeddist = float(toks[9])
            motifs[-1].seednum  = int(toks[10][0:-1])
            motifs[-1].pvalue   = math.pow(10,-float(toks[12]))

            if 'ch:' in toks:
                _idx = toks.index('ch:')
                motifs[-1].church = math.pow(10,-float(toks[_idx+1]))
            if 'Es:' in toks:
                _idx = toks.index('Es:')
                motifs[-1].E_site = math.pow(10,-float(toks[_idx+1]))
            if 'x2:' in toks:
                _idx = toks.index('x2:')
                motifs[-1].E_chi2 = math.pow(10,-float(toks[_idx+1]))
            if 'Eq:' in toks:
                _idx = toks.index('Eq:')
                motifs[-1].E_seq = math.pow(10,-float(toks[_idx+1]))
            if 'mn:' in toks:
                _idx = toks.index('mn:')
                motifs[-1].MNCP = float(toks[_idx+1])
            if 'f:' in toks:
                _idx = toks.index('f:')
                motifs[-1].frac = float(toks[_idx+1])
            if 'Ra:' in toks:
                _idx = toks.index('Ra:')
                motifs[-1].ROC_auc = float(toks[_idx+1])
            if 'cR:' in toks:
                _idx = toks.index('cR:')
                motifs[-1].CRA     = float(toks[_idx+1])
            if 'Cf:' in toks:
                _idx = toks.index('Cf:')
                motifs[-1].Cfrac   = float(toks[_idx+1])
            if 'k:' in toks:
                _idx = toks.index('k:')
                motifs[-1].kellis  = float(toks[_idx+1])

            if 'b:' in toks:
                _idx = toks.index('b:')
                motifs[-1].numbound = int(toks[_idx+1])
            if 'nG:' in toks:
                _idx = toks.index('nG:')
                motifs[-1].nummotif = int(toks[_idx+1])
            if 'bn:' in toks:
                _idx = toks.index('bn:')
                motifs[-1].numboundmotif = int(toks[_idx+1])



        if lines[i][0:10] == 'Threshold: '[0:10]:
            toks =  lines[i].split()
            motifs[-1].threshold= float(toks[1])
        if lines[i][0:5] == 'Seed '[0:5]:
            toks = lines[i].split()
            id = int(toks[1][0:-1])  #'10:' -> '10'
            seedD[id] = toks[2]
        if lines[i][0:7] == 'Source: '[0:7]:
            motifs[-1].source = lines[i][7:].strip()
        if lines[i][0:6] == 'Gamma: '[0:6]:
            motifs[-1].gamma = float(lines[i][6:])
        if lines[i][0:6] == 'Evalue: '[0:6]:
            motifs[-1].evalue = float(lines[i][7:].strip())
        if lines[i].find('Using')>=0 and lines[i].find('as seeds')>=0:
            '''#Using all (132) motifs in SLT_081503.seeds as seeds:'''
            seedfile = lines[i].split()[-3]
        if lines[i][0:6] == 'thetas: '[0:6]:
            thetas = lines[i][7:].split()
            for theta in thetas:
                motifs[-1].thetas.append(float(theta))
    for i in range(len(motifs)):
        if seedfile: motifs[i].seedfile = seedfile
        seednum = motifs[i].seednum
        if seedD.has_key(seednum):
            motifs[i].seedtxt = seedD[seednum]
    return(motifs)
    
def save_motifs(motifs,filename,kmer_count=20):
    _old_stdout = sys.stdout  #Cache REAL stdout
    try:
        sys.stdout  = open(filename,'w')
        print_motifs(motifs,kmer_count)
        sys.stdout.close()
    except:
        print '!-- Error saving motifs to %s'%filename
        raise
    sys.stdout  = _old_stdout
    
def print_motif(motif,kmer_count=20,istart=0):
    print_motifs([motif],kmer_count,istart)
    sys.stdout.flush()

def print_motifs(motifs,kmer_count=20,istart=0):
    i = istart-1
    for m in motifs:
        i = i + 1
        print "Log-odds matrix for Motif %3d %s"%(i,m)
        m._print_ll()
        print "Sequence Logo"
        m._print_bits()
        for newprop in ['gamma', 'church', 'E_site', 'E_seq', 'E_chi2', 'realpvalue',
                        'kellis', 'MNCP', 'ROC_auc', 'CRA', 'Cfrac', 'frac']:
            if not m.__dict__.has_key(newprop):   #Kludge to deal w/ old shelves
                m.__dict__[newprop] = None
        if m.seedtxt:  print "Seed: %3d %s"%(i,m.seedtxt)
        if m.gamma:    print "Gamma: %7.5f"%m.gamma
        if m.evalue != None: print 'Evalue: %6.3e'%m.evalue
        if m.family:   print "Family: ",m.family
        if m.source:   print "Source: ",m.source
        #Motif   0 NGAGGGGGNN (0)            (Bits:   8.24   MAP:   6.53   D:  0.21  0)  Enr: 54.000 
        print "Motif %3d %-25s (Bits: %5.2f  MAP: %5.2f   D: %5.3f  %2d) E: %6.3f"%(
            i, m, m.totalbits, m.MAP, m.seeddist, m.seednum, nlog10(m.pvalue)),
        if m.church != None:  print ' ch: %5.2f'%nlog10(m.church),
        if m.frac   != None:  print ' f: %5.2f'%(m.frac),
        if m.E_site != None:  print ' Es: %5.2f'%nlog10(m.E_site),
        if m.E_seq != None:  print ' Eq: %5.2f'%(nlog10(m.E_seq)),
        if m.MNCP   != None:  print ' mn: %5.2f'%(m.MNCP),
        if m.ROC_auc!= None:  print ' Ra: %6.4f'%(m.ROC_auc),
        if m.E_chi2 != None:
            if m.E_chi2 == 0: m.E_chi2=1e-20
            print ' x2: %5.2f'%(nlog10(m.E_chi2)),
        if m.CRA    != None:  print ' cR: %6.4f'%(m.CRA),
        if m.Cfrac  != None:  print ' Cf: %6.4f'%(m.Cfrac),
        if m.realpvalue != None: print ' P: %6.4e'%(m.realpvalue),
        if m.kellis != None:  print ' k: %5.2f'%(m.kellis),
        try:
            if m.numbound      :  print ' b: %3d'%(m.numbound),
            if m.nummotif      :  print ' nG: %3d'%(m.nummotif),
            if m.numboundmotif :  print ' bn: %3d'%(m.numboundmotif),
        except: pass
        print
        if m.thetas != []:
            tstr = "thetas:"
            for theta in m.thetas:
                tstr = tstr + " " + str(theta)
            print tstr
        _max = m.maxscore
        m.maxscore = -100
        if kmer_count >= 0:  seqs = m.bogus_kmers(kmer_count)
        else:                seqs = m.seqs
        for seq in seqs:     print seq,i,m.scan(seq)[2][0]
                
        m.maxscore = _max
        print '*'*m.width
        print "MAP Score: %f"%(m.MAP)
        sys.stdout.flush()

def nlog10(x,min=1e-323):
    if x < min: x=min
    try:
        return math.fabs(math.log(x)/math.log(10))
    except:
        return 0

def txt2motifs(txt,VERBOSE=1):
    motifs = []
    exists = os.path.exists
    toks   = txt.split(':')
    if exists(toks[0]):               #It's a file!!
        fname = toks[0]
        if fname.find('.pickle') > 0: #It's a pickle!!
            return pickletxt2motifs(toks)
        else:                         #It's a "Motif" file!!
            if VERBOSE:
                print "# Loading motif from %s"%fname
            allmotifs = load(fname)
        if len(toks) == 1: motifs = allmotifs
        else:
            idxs   = [int(x) for x in toks[1].split(',')]
            motifs = [allmotifs[x] for x in idxs]
    else:                             #It's a text string!!
        fname = 'TXT'
        for t in txt.split(','):
            motifs.append(Motif_from_text(t))
    for i in range(len(motifs)): motifs[i].index = i
    for i in range(len(motifs)): motifs[i].file = fname
    return motifs

def pickletxt2motifs(toks):
    fname = toks[0]
    print "# Loading motif pickle from %s"%fname
    F = open(fname,'r')
    DA = pickle.load(F)
    F.close()
    ans = []
    if type(DA) == type({}):
        if len(toks) > 1:
            keys = [x.replace('%',' ') for x in toks[1].split(',')]
            for k in keys: ans.append(DA[k])
        else:
            for k in DA.keys(): DA[k].key = k
            ans = DA.values()
    else: #Assuming DA is a list
        if len(toks) > 1:
            idxs = [int(x) for x in toks[1].split(',')]
            ans  = [DA[x] for x in idxs]
        else:
            ans  = DA
    return ans
    

if __name__ == '__main__': main()


