
import sys, re, os, math, time, string, tempfile, pickle, random

import Numeric
import time
import CVGsupport                #From local directory
import ConvergeMotifTools
from forKenzie import KenzieSequences

'''
This file contains class definitions for 5 classes
class EM:                  Implementation of an EM

class MarkovBackground:    Implementation of a n-th order Markov Background model

class Probe:               A probe (string with additional baggage)

class Wmer:                A wmer  (string with additional baggage)

class MDscanCandidate:     Essentially a "super" Motif, that can compute its MAP score
'''
ConvergeDataDir='/home/romer/romerweb/ConvergeSupport/data/'
theMarkovBackground = None
MarkBG = {}

class MDscanCandidate:
    '''
    A candidate for MDscan consists of the following:
       A set of segments "wmers"
       A pssm computed from these segments
       A score (probably the MAP score described in the literature)
       General propeties:
         # segments
         Width

       A candidate must also be able to:
         Evaluate its own score
         Modify itself:
           Add wmer
           Remove wmer      
    '''
    
    def __init__(self, wmers=''):
        self.MAP  = 0
        self.log2 = math.log(2)
        self.Q = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31} #Yeast defaults
        self.last_bgprob = 0
        self.deltaMAP = 0.0

        global theMarkovBackground
        if theMarkovBackground:
            self.bgprob = theMarkovBackground.zeroth()

        if wmers:
            self.wmers    = wmers
            self._update()
            self.needs_update = 1

    def find_wmers(self,seqs):
        self.wmers = []
        for seq in seqs:
            (bestscore,bestmatch) = (0,'')
            (matches,endpoints,scores) = self.pssm._scan(seq,-10)
            for match,score in zip(matches,scores):
                if score > bestscore:
                    bestscore = score
                    bestmatch = match
            if bestmatch:
                self.wmers.append(bestmatch)
        self._update()

    def _plogp_new(self,wmers):
        nwmers = float(len(wmers))
        toAscii = {'A':97, 'C':99, 'G':103, 'T':116}
        sumsD = {}
        AsciiSeq = Numeric.array(wmers).astype(type(1))
        plogp = 0
        for key in toAscii.keys():
            sumsD[key] = Numeric.sum(Numeric.equal(AsciiSeq,toAscii[key])).astype(type(1.))
            f = sumsD[key]/nwmers + 0.0000000001
            plogp = plogp + Numeric.sum(Numeric.log(f) * f / Numeric.log(2))
        return(plogp)

    def _plogp(self,wmers):
        nwmers = float(len(wmers))
        plogp  = 0
        for i in range(len(wmers[0])):
            C = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'N':0}
            for wmer in wmers:
                C[wmer[i]] = C[wmer[i]] + 1
            for key in C.keys():
                if key == 'N': continue
                f = C[key]/nwmers
                if f > 0:
                    plogp = plogp + f * math.log(f)/self.log2
        return(plogp)

    def _plogp_memo(self, wmer, count):
        nwmers = float(len(self.wmers) + count)
        plogp = 0
        for i in range(len(wmer)):
            letter = wmer[i]
            for key in ['A', 'C', 'T', 'G']:
                lettercount = self.pssm.counts[i][key]
                if key == letter:
                    lettercount = lettercount + count
                #print lettercount, nwmers
                f = float(lettercount)/nwmers
                if f > 0:
                    plogp = plogp + f * math.log(f)/self.log2
        return(plogp)


    def computeMAP_memo(self, wmer, count):
        wmers = self.wmers
        nwmers = float(len(wmers))  + count
        width  = float(len(wmers[0]))
        plogp  = self._plogp_memo(wmer,count)
        bgprob = self.last_bgprob + count * self.logPbackground(wmer)
        #print "MEMO: %d %f %f"%(nwmers,plogp,bgprob)
        MAP = math.log(nwmers)/ (self.log2 * width) \
              * (plogp - (1.0/nwmers) * bgprob)
        #MAP = nwmers/ width * (plogp - (1.0/nwmers) * bgprob - math.log(75)/math.log(2))
        return(MAP,bgprob)
        
    def _update(self):
        self.nwmers   = len(self.wmers)
        self.pssm     = ConvergeMotifTools.Motif(self.wmers,self.Q)
        self.MAP      = self.computeMAP() 
        self.pssm.MAP = self.MAP
        
    def _recompute(self,_seqs):
        self.nwmers       = len(_seqs)
        self.wmers        = _seqs
        self.pssm         = ConvergeMotifTools.Motif(_seqs,self.Q)
        #self.last_bgprob  = _bgprob
        self.MAP          = self.computeMAP()
        self.pssm.MAP     = self.MAP
        self.needs_update = 1

    def computeMAP(self,in_wmers=''):
        if in_wmers:
            wmers = in_wmers
        else:
            wmers = self.wmers
        nwmers = float(len(wmers))
        width  = float(len(wmers[0]))
        plogp  = self._plogp(wmers)
        bgprob = 0
        for wmer in wmers:
            bgprob = bgprob + self.logPbackground(wmer)
        #print "ORIG: %d %f %f"%(nwmers,plogp,bgprob)
        MAP = math.log(nwmers)/ (self.log2 * width) * (plogp - (1.0/nwmers) * bgprob)
        if not in_wmers:
            self.last_bgprob = bgprob
        return(MAP)
    
    def logPbackground(self,wmer):
        if theMarkovBackground:
            Ptot = theMarkovBackground.logbackground(wmer)
        else:
            Ptot = math.log(0.25)/self.log2 * len(wmer)
        return(Ptot)
    def __repr__(self):
        s = ''
        s = s + '%-10s  (Bits: %6.2f   MAP: %6.2f   D: %5.2f %2d)'%\
                 (self.pssm.__repr__(), self.pssm.totalbits,self.MAP,
                  self.pssm.seeddist,self.pssm.seednum)
        return s

    def check_and_update(self,wmer,count,verbose=''):
        if self.has_wmer(wmer):
            return
        scanscore,maxscore = self.pssm._scan(wmer,-1000,"FORW_ONLY")[2][0],self.pssm.maxscore
        if (scanscore/maxscore) < 0:  #This is a shortcut: Use Scan to decide
            return()                  #Whether the full MAP calculation is worthwhile
        
        # Does adding one copy help?  (we don't want to overwhelm with
        # all copies.)
        #_wmers = self.wmers[:]
        #_wmers.append(wmer)
        #_MAP  = self.computeMAP(_wmers)
        (_MAP,bgprob) = self.computeMAP_memo(wmer,1)
        #print ('%8.4f %8.4f %8.4f %8.4f ')%(
        #    self.MAP,_MAP,self.computeMAP(_wmers), \
        #    _MAP-self.computeMAP(_wmers))
        delta = _MAP - self.MAP
        if 0:  #Used for printing scanscore vs MAP score stats
            print '## %10.5f %10.5f %10.5f   %10.5f %10.5f %10.5f'%(
                scanscore,maxscore,scanscore/maxscore,_MAP,self.MAP,_MAP/self.MAP),
            if _MAP > self.MAP:
                print 1,wmer
            else:
                print 0,wmer
            return()
        
        if _MAP > self.MAP + self.deltaMAP:
            # The one wmer helps the map score.  So add all its copies
            # This seems like the right thing to do... doesn't it?  How
            # could we justify only adding _some_ of the copies?
            _wmers = self.wmers[:]
            for i in range(count):
                _wmers.append(wmer)
            _pssm = ConvergeMotifTools.Motif(_wmers,self.Q)
            (_MAP,self.last_bgprob)  = self.computeMAP_memo(wmer,count)
            if verbose: print 'Adding %d copies of %s: O:%f N:%f (delta:%f)'% \
               (count,wmer,self.MAP,_MAP,delta)
            sys.stdout.flush()
            self.nwmers  = len(_wmers)
            self.wmers   = _wmers
            self.pssm    = _pssm
            self.MAP     = _MAP
            self.pssm.MAP= _MAP
            self.needs_update = 1
            return("ACCEPTED")

    def MAPscan(self,nmers):
        (maxMAP,maxnmer) = (0,'')
        for nmer in nmers:
            scanscore = self.pssm._scan(nmer,-1000,"Forw_Only")[2][0]
            if scanscore < 0: continue
            (_MAP, bgprob) = self.computeMAP_memo(nmer,1)
            if _MAP > maxMAP:
                maxMAP  = _MAP
                maxnmer = nmer
        return(maxnmer)

    def MAPpurge(self,verbose=''):
        Tot_Removed = 0
        while 1:
            wmers    = self.wmers
            new_list = self._MAPpurge_list()
            delta = len(wmers) - len(new_list)
            Tot_Removed = Tot_Removed + delta
            if delta == 0:  #Nothing purged
                break
            if verbose:
                print "\tPurged %4d (Total %4d) sequences (leaving %4d) %s"% \
                      (delta,Tot_Removed,len(new_list),self.pssm.oneletter)
                sys.stdout.flush()
            self._recompute(new_list)
            #print self

    def _MAPpurge_list(self):
        omit_list = []
        for wmer in self.wmers:
            count = 0
            for _wmer in self.wmers:
                if wmer == _wmer:
                    count = count + 1
            if count == len(self.wmers):  #Only one type of wmer in this thing
                return(self.wmers)        #Don't purge it!!!
            (_MAP,_bgprob) = self.computeMAP_memo(wmer,-count)
            if _MAP > self.MAP:
                omit_list.append(wmer)
        new_list  = []
        for wmer in self.wmers:
            if not (wmer in omit_list):
                new_list.append(wmer)
        return(new_list)

    def purge(self,verbose = ''):
        if verbose: print "Purging %s"%self
        #Build list of wmers
        wmersD = {}
        for wmer in self.wmers:
            if not wmersD.has_key(wmer):
                wmersD[wmer] = 1
        for wmer_omit in wmersD.keys():
            #Build temporary sequence list
            _seqs = []
            for wmer in self.wmers:
                if wmer != wmer_omit:
                    _seqs.append(wmer)
            count = len(self.wmers) - len(_seqs)
            if len(_seqs) == 0: continue
            (_MAP,_bgprob)  = self.computeMAP_memo(wmer_omit,-count)
            #print ('%8.4f %8.4f %8.4f %8.4f ')%(
            #    self.MAP,_MAP,self.computeMAP(_seqs), \
            #    _MAP-self.computeMAP(_seqs))
            #_MAP  = self.computeMAP(_seqs)
            if _MAP > self.MAP:
                if verbose: print 'Purging %d copies of %s: O:%f N:%f (delta:%f)'% \
                   (count,wmer_omit,self.MAP,_MAP,_MAP-self.MAP)
                sys.stdout.flush()
                self.nwmers  = len(_seqs)
                self.wmers   = _seqs
                self.pssm    = ConvergeMotifTools.Motif(_seqs,self.Q)
                self.last_bgprob = _bgprob
                self.MAP     = _MAP
                self.pssm.MAP= _MAP
                self.needs_update = 1
        
    def has_wmer(self,wmer):
        rc = ConvergeMotifTools.revcomplement(wmer)
        if (wmer in self.wmers) or (rc in self.wmers):
            return(1)
        else:
            return(0)

class EM:
    def __init__(self,seed_seqs, all_seqs, probes, width = 6, verbose = ''):
        self.seed_seqs  = seed_seqs #Sequences to be scanned for seeds
        self.seqs       = all_seqs
	self.probeids	= probes
        self.candidates = []
        self.models     = []      #Set directly or computed from seed_seqs
        self.width      = width
	self.probes 	= []
        self.verbose    = verbose
        if width:
            self.goodwmersT = ConvergeMotifTools.top_nmers(self.width,self.seed_seqs,1,"")
        else:
            self.goodwmersT = zip(self.seed_seqs,range(len(self.seed_seqs)))
        self.bgprob     = {'A': 0.31, 'C': .19, 'G': .19, 'T': .31}
        self.beta       = 0.001
        self.cbeta      = 0.001
        self.deltamin   = 1e-3
        self.method     = "ZOOPS" # OOPS or ZOOPS )
        self.param      = {}
        self.gapflank   = 0
        self.gapweight  = 0.2
        self.seedbeta   = 0.02
        self.joint      = 1
        self.numgenomes = 0

        global theMarkovBackground
        global MarkBG
        if theMarkovBackground:
            self.bgprob = theMarkovBackground.zeroth()

	if theMarkovBackground.species == 'YEAST':
        	self.thetas     = [1.0, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50]
                #self.thetas     = [1.0, 0.77, 0.67, 0.61, 0.64, 0.45, 0.45]
                self.zeta1s     = [0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                self.zeta2s     = [0.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
	elif theMarkovBackground.species == 'human':
		self.thetas	= [1.0, 0.37, 0.35]
                self.zeta1s     = [0.0, 0.1, 0.1]
                self.zeta2s     = [0.0, 0.5, 0.5]
        elif theMarkovBackground.species == 'Ciona':
            self.thetas = [1.0, 0.5]
            self.zeta1s = [0.0, 0.1]
            self.zeta2s = [0.0, 0.5]
	else:
		self.thetas     = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        '''DELETE
        if all_seqs:  #Should we Go?
            self.EM()
        '''

    def report(self):
        if self.width:
            print "# A total of %5d %d-mer seeds specified or found in top %3d sequences"% \
                  (len(self.goodwmersT),self.width,len(self.seed_seqs))
        else:
            print "# A total of %5d seeds specified or found in top %3d sequences"% \
                  (len(self.goodwmersT),len(self.seed_seqs))
            
        print "# Using beta = %f   and  Converging once d(theta) < %f"% \
              (self.beta, self.deltamin)
        if theMarkovBackground:
            print "# Using High-order (%d-th)  Markov Background model"% \
                  (theMarkovBackground.highestorder -1)
        if self.gapflank > 0:
            print "# Using Gapped model: flank %d, weight %3.1f"%(self.gapflank, self.gapweight)
        sys.stdout.flush()

    def seed_models(self):
        V = self.verbose
        beta = self.seedbeta
        for nmer,count in self.goodwmersT:
            #print "nmer: %s beta: %f"%(nmer,beta)
            M = ConvergeMotifTools.Motif(None,self.bgprob) #No automatic computation
            M.compute_from_nmer(nmer,beta)
            self.models.append(M)

    def calcmask(self,width):
        self.mask = map(lambda x:1,range(width))
        if self.gapflank > 0:
            gf = self.gapflank
            gw = self.gapweight
            for pos in range(gf,width-gf):
                self.mask[pos] = gw
            self.mask[0]  = gw
            self.mask[-1] = gw

    def EM_Cstart(self):
        verbose = self.verbose
        if verbose:
            print "Seeding models..."
            sys.stdout.flush()
        self.seed_models()

        #Initialize parameters
        if not self.param.has_key('gamma'): self.param['gamma'] = 0.2
        timings = {'Probes':0, 'Background':0, 'C EM':0, 'Post':0}
        _time = time.time()

	#USE THIS FOR MAMMALS!!!! LIKE HUMANS ETC.....
	if theMarkovBackground.species == 'human':
                self.numgenomes = 3
                '''
		af = open(ConvergeDataDir+'upstream2000.maf') 
    		seqDict = {}
    		line = af.readline()
		while line != "":
			while line[0] == 'a':
				line = af.readline()
				if line[0] == 's':
					line = line.split()
					tag = line[1]
					seqH = line[len(line)-1]
					line = af.readline()
					line = line.split()
					seqM = line[len(line)-1]
					line = af.readline()
					line = line.split()
					seqR = line[len(line)-1]
                                	#This block of code is used to expunge gaps in the human sequence
					seqH = seqH.replace('.','-')
                                	dash = seqH.find('-')
                                	while dash != -1:
                                	    seqH = seqH[:dash]+seqH[dash+1:]
                                	    seqM = seqM[:dash]+seqM[dash+1:]
                                	    seqR = seqR[:dash]+seqR[dash+1:]
                                	    dash = seqH.find('-')
                                	#end of gap elimination stuff
				seqDict[tag] = [('human', seqH),('mouse',seqM.replace('.','-')),('rat',seqR.replace('.','-'))]
				line = af.readline()
                        line = af.readline()
                af.close()
		pf = open('humanmouserat.pickle', 'w')
		pickle.dump(seqDict, pf)
		pf.close()
                '''
		pf = open('humanmouserat.pickle', 'r')
		seqDict = pickle.load(pf)
	elif theMarkovBackground.species == 'YEAST':
		sequence_repository=KenzieSequences(ConvergeDataDir+'ManoliProbeAlignments_050105.pickle')
                self.numgenomes = 4
        elif theMarkovBackground.species == 'Ciona':
            self.numgenomes = 2
            seq_file = open(ConvergeDataDir+'CionaSeq.pickle','r')
            sequence_repository = pickle.load(seq_file)
            seq_file.close()
	discard = 0

	for probe_name in self.probeids:
            #try:
            if theMarkovBackground.species == 'human':
                try:
                    mam_seqs = seqDict[probe_name]
                    ms = []
                    ms.append(mam_seqs[0])
                    #print mam_seqs[0][0]
		    seq_list  = ms
                except:
                    discard = discard + 1
                #print "%s\n"%seq
            elif theMarkovBackground.species == 'YEAST':
                seq_list= sequence_repository.returnSequencesForProbes(probe_name)
            elif theMarkovBackground.species == 'Ciona':
                seq_list = sequence_repository[probe_name]
            #print seq_list
            P = Probe(seq_list)
            #except:
                #discard = discard + 1   
            if (len(P)>0):
                self.probes.append(P)
            else:
                discard = discard + 1

	print "Number of discarded Probes: %i"%discard
	print "number of probes: %i"%len(self.probes)

	#for seq in self.seqs:
        #    P = Probe(seq)
        #    self.probes.append(P)              
       
        _time2 = time.time(); timings['Probes'] = _time2-_time; _time = _time2

        if verbose: print "Optimizing candidates by EM."
        if verbose: sys.stdout.flush()

        c_logZ_sets = {}
	#consbg_set = {}
        for Model,i in zip(self.models,range(len(self.models))):
            width = Model.width
	    self.calcmask(width)              
            if not c_logZ_sets.has_key(width):
                	c_logZs_set = []
			if verbose: print "#%s   |%s|"%(' '*28,'-'*len(self.seqs))
			if verbose: sys.stdout.flush()
			if verbose: print "Computing background (width %2d)  "%width,
			for P in self.probes:
				if verbose: sys.stdout.write('.')
				if verbose: sys.stdout.flush()	
				logZs = []
				for seq,label in zip(P.seqs,P.labels):
					logZs[len(logZs):len(logZs)] = self.all_Wmers(width,seq,label)
				for seq, label in zip(P.seqs,P.labels):
					s = seq.upper()
                			s = MarkBG[label].replace_N(s)
					#s = s.replace('-','')
					logZs.append(MarkBG[label].logbackground(s))
				c_logZs = CVGsupport.list2double(logZs)
				c_logZs_set.append(c_logZs)
			c_logZ_sets[width] = c_logZs_set
			if verbose: print

	    c_logZ_set = c_logZ_sets[width]

	    for P,c_logZs in zip(self.probes,c_logZ_set):
               	P.c_wmerbgs = c_logZs
            
            _time2 = time.time()
            timings['Background'] = timings['Background'] +_time2-_time
            _time = _time2

            '''Perform EM'''
            _time  = time.time()
            if (len(self.probes)>0):
                newModel = self.EM_C(Model, self.probes)
            _time2 = time.time(); timings['C EM'] = timings['C EM'] + _time2-_time; _time = _time2

            #print "cLL: ",newModel.joint
            #print "pLL: ",self.compute_joint(newModel,Wmers_by_seq)

            '''Was there a problem?'''
            if newModel == None:
                continue


            '''Set various things in PSSM'''
            #Distance(s)
            seeddist = ConvergeMotifTools.infomaskdiff(newModel,Model)
            print '%s ----> %s'%(Model,newModel)
            print "Seed %2d: %s  -->  %s  mask:%9.5f  infoMask:%9.5f d:%9.5f"%(
                i, Model, newModel,
                ConvergeMotifTools.maskdiff(newModel,Model),
                ConvergeMotifTools.infomaskdiff(newModel,Model), #order is important
                Model-newModel)
            #Seed
            if Model.seedtxt: newModel.seedtxt = Model.seedtxt
            if Model.source:  newModel.source  = Model.source
           
            #newModel.denoise()
            newModel.seeddist = seeddist
            newModel.seednum  = i
            print newModel
            newModel._print_p()
            newModel._print_ll()

            '''Set various things in Candidate (like a wrapper for PSSM)'''
            C = MDscanCandidate()
            C.pssm = newModel.copy()
            #C.wmers = self.best_by_Z(Wmers_by_seq)
            C.wmers  = [newModel.emit() for junk in range(20)]
            #C._update()  #MAJOR REMOVAL????????? DBG 10-14-03
            #C.MAPpurge()
            C.pssm = newModel.copy()  
            self.candidates.append(C)
            _time2 = time.time(); timings['Post'] = timings['Post']+_time2-_time;_time = _time2

        '''Print Timing Information'''
        if verbose:
            print "# Timing Information"
            _t = 0
            for timing in timings.keys():
            	_t = _t + timings[timing]
            for timing in timings.keys():
            	print "# %12s %f  %f%%"%(timing,timings[timing],timings[timing]*100/_t)

    '''
    def calc_cons_prob(self, seq):
	logP = 0
	for i in range(len(seq)):
		logP = logP + math.log(self.cbgprob[seq[i]])/math.log(2.0)
	return logP
    '''

    def compute_joint(self,model,Wmers_by_seq):
        sum = 0
        #param_gamma = self.param['gamma']
        gamma = model.gamma
        for Wmer_set in Wmers_by_seq:
            print len(Wmer_set)
            param_lambda = float(gamma) / len(Wmer_set[0].src)
            Qi = 0 

            for wmer in Wmer_set:
                sum = sum + wmer.Z * wmer.score * wmer.count
                Qi  = Qi + wmer.Z
                #if wmer.Z > 0.5:
                #    print "Z:%f s:%f c:%f\t%s\t%f"%(wmer.Z,wmer.score,wmer.count,wmer,sum)
            sum = sum + (1 - Qi) * Wmer_set[0].srcQ
            sum = sum + (1 - Qi) * math.log(gamma)/math.log(2)
            sum = sum +      Qi  * math.log(param_lambda)         /math.log(2)

        return(sum)
                

    def EM_C(self, Model, probes, store_Zs=''):
        width = Model.width
        bg    = self.bgprob

        '''Build Probelist'''
        c_Probelist = CVGsupport.Probelist()
        for m_probe in probes:
            c_Probelist.append(m_probe.c_intA, len(m_probe), m_probe.genomes, m_probe.logP, m_probe.c_wmerbgs, m_probe.seqsbgs)            
        '''Build C++ PSSM'''
        c_Model          = CVGsupport.SeqMat(width, self.numgenomes)
        c_Model.gamma    = self.param['gamma']
        c_Model.gammawt  = 0.8
        c_Model.deltamin = self.deltamin
        c_Model.beta     = self.beta
        #c_Model.num_genomes    = self.num_genomes
        c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'],0.0,0)
        if (theMarkovBackground.species=='YEAST'):
            indexes = {'Scer': 0, 'Spar': 1, 'Smik': 2, 'Sbay': 3, 'Skud': 4}
            #'Scas': 5, 'Sklu': 6}
        for label in MarkBG.keys():
            bg = MarkBG[label].zeroth()
            index = indexes[label]
            if (label == 'Scer'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], 0, index)
            elif (label == 'Smik'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], bg['-'], index)
            elif (label == 'Sbay'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], bg['-'], index)
            elif (label == 'Spar'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], bg['-'], index)
            elif (label  == 'Skud'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], bg['-'], index)
                #c_Model.setBg(0.251, 0.150, 0.150, 0.251, 0.197, index)
            elif (label  == 'Scas'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], bg['-'], index)                
                #c_Model.setBg(0.256, 0.131, 0.131, 0.256, 0.226, index)
            elif (label  == 'Sklu'):
                c_Model.setBg(bg['A'], bg['C'], bg['G'], bg['T'], bg['-'], index)                
                #c_Model.setBg(0.247, 0.154, 0.154, 0.247, 0.198, index)
            elif (label == 'savignyi'):
                c_Model.setBg(0.32, 0.18, 0.18, 0.32, 0.0, index)
            elif (label == 'intestinalis'):
                c_Model.setBg(0.249, 0.149, 0.149, 0.249, 0.204, index)
	    elif (label == 'human'):
		c_Model.setBg(0.2352, 0.2648, 0.2648, 0.2352, 0, index)
	    elif (label == 'mouse'):
		c_Model.setBg(0.1388, 0.1537, 0.1537, 0.1388, 0.4149, index)
  	    elif (label == 'rat'):
		c_Model.setBg(0.1352, 0.1481, 0.1481, 0.1352, 0.4334, index)
        
        #Set the value of theta and zetas for each aligned genome
	i = 0
	for val in self.thetas:
        	c_Model.setTheta(val, i)
		i = i+1
                
        i=0
        for zeta1, zeta2 in zip(self.zeta1s, self.zeta2s):
            c_Model.setZeta(zeta1, zeta2, i)
            i = i+1
        
        '''Gap maddness'''
        for i in range(width):
            c_Model.setmask(i,self.mask[i])

        '''Copy python PSSM data into C++ PSSM'''
        Ljs = zip(['A','C','G','T'],[0,1,2,3])
        for i in range(width):
            for L,j in Ljs:
                c_Model.set(i,j, Model.logP[i][L])

        '''Do EM until convergence (c_Model.deltamin)'''
        c_Model.EMstep(c_Probelist, self.param['gamma'])

        '''Store the Zs with each Wmer, which is useful for doing'''
        '''some post-processing later.'''
        if 0 and store_Zs:
            print '#Re-organizing Z data...',;sys.stdout.flush()
            for Wmer_set,i in zip(Wmers_by_seq, range(len(Wmers_by_seq))):
                for wmer,j in zip(Wmer_set, range(len(Wmer_set))):
                    wmer.Z = c_Probelist.get_Z(i,j)
            print 'Done.'; sys.stdout.flush()

        '''Build new python Motif from data in updated C++ Motif'''
        newW = []
        Ljs = zip(['A','C','G','T'],[0,1,2,3])
        _t = 0
        for i in range(width):
            d = {}
            for L,j in Ljs:
                d[L] = c_Model.get_c(i,j)
                _t   = _t + d[L]
            newW.append(d)

        if _t == 0:
            return None

        Mtmp = ConvergeMotifTools.Motif(None,self.bgprob)
        Mtmp.compute_from_counts(newW,self.beta)
        Mtmp.gamma = c_Model.gamma
        for theta, i in zip(self.thetas,range(len(self.thetas))):
            Mtmp.thetas.append(c_Model.get_theta(i))
        Mtmp.joint = c_Model.joint
        '''Pass it back'''
        return(Mtmp)

    def best_by_Z(self,Wmers_by_seq):
        ans = []
        for Wmer_set in Wmers_by_seq:
            bestwmer,bestZ = '', 0
            for wmer in Wmer_set:
                if wmer.Z > bestZ:
                    bestZ = wmer.Z
                    bestwmer = wmer.seq
            if bestwmer: ans.append(bestwmer)
            if 0 and self.verbose: print "Best: %s %f"%(bestwmer,math.log(bestZ))
        return(ans)

    def all_Wmers(self,N,seq,label):
        forw = []
        rev  = []
        #seqrc = ConvergeMotifTools.revcomplement(seq)
        '''Mlh = theMarkovBackground.highestorder
        Mlb = theMarkovBackground.logbackground
        MCP = theMarkovBackground.CP'''
        Mlh = MarkBG[label].highestorder
        Mlb = MarkBG[label].logbackground
        MCP = MarkBG[label].CP

	Mr = MarkBG[label].replace_N
	s1   = seq.replace('.','-')
	#s1 = s1.replace('.','')
	s2 = ''
        Fbg = Mlb(s1)
        #Rbg = Mlb(seqrc)
        nmask = map(lambda x:1-x, self.mask)

        '''
        ?? QUESTION: Is it sensible to compute the background probabilities
        this way?
        
        1) BG of complementary strand is taken as equal to primary strand.
        2) Letters inside the motif window are not used for conditional probabilities.
           As a result, the calculation essentially breaks down to the log probability the
           background emits the sequence to the left of the window plus the log probability
           the background emits the sequence to the right.
        3) I\'ve worked out an efficient way to compute this by
           a) Compute the background probability for the entire probe/sequence
           b) (Quick) Compute logQdiff below
           c) Subtract
        '''

        for i in range(len(seq)-N+1):
            subseq = seq[i:i+N]
            '''Build Wmer information'''
            #Wtmp        = Wmer(subseq)
            left        = seq[0:i]
            right       = seq[i+N:]
            #Wtmp.lflank = left
            #Wtmp.rflank = right
            #if i==0: Wtmp.src    = seq
            #Wtmp.srcQ   = Fbg
            #Wtmp.i      = i

            '''This is the fast way'''
	    subseq = Mr(subseq)
	    subseq = subseq.replace('.','-')
	    #subseq = subseq.replace('.','')
	    s1 = left[-Mlh:] + subseq + right[0:Mlh]
	    #s1 = s1.replace('-','')
	    s1 = s1.replace('.','-')
	    s2 = left[-Mlh:]
	    #s2 = s2.replace('-','')
	    s2 = s2.replace('.','-')
	    s3 = right[0:Mlh]
	    #s3 = s3.replace('-','')
	    s3 = s3.replace('.','-')
	    logQdiff = Mlb(s1) - Mlb(s2) - Mlb(s3)
            #logQdiff = Mlb(left[-Mlh:] + subseq + right[0:Mlh]) - Mlb(left[-Mlh:]) - Mlb(right[0:Mlh])
            logQtot = Fbg - logQdiff

            '''Add a bit back for intervening bases in the "gap" '''
            gapbg = 0
            for p in range(len(subseq)):
                gapbg = gapbg + MCP[subseq[p]] * nmask[p]
            logQtot = logQtot + gapbg
            
            '''Build Wmer-reverse complement information'''
            #Wtmprc = Wmer(Wtmp.rc)
            #Wtmprc.lflank = seqrc[0:-(i+N)]  #Check this in case it is ever necessary
            #if i!=0:
            #    Wtmprc.rflank = seqrc[-i:]   #Necessary [11-12-02]
            #else:
            #    Wtmprc.rflank = ''
            #Wtmprc.logQtot = Wtmp.logQtot
            #Wtmprc.srcQ    = Wtmp.srcQ
            #Wtmprc.i       = i
            forw.append(logQtot)
            rev.append(logQtot)
        W = []
        W.extend(forw)
        W.extend(rev)
        #seq.c_wmerbgs = CVGsupport.list2double(map(lambda x: x.logQtot, W))
        #CVGsupport.printdouble(seq.c_wmerbgs,len(W))
        return(W)

        
class MarkovBackground:
    def __init__(self,file,numprobes,species='YEAST',seqs=''):
        self.highestorder = 4
        if   species[0:5] == 'YEAST':
            #self.sourcefile = ConvergeDataDir+'yeast.nc.6.freq'
            if numprobes<500:
                self.sourcefile = ConvergeDataDir+'Scer.bg'
            else:
                self.sourcefile = file.split('.')[0] + '.Scer' + '.bg'
        elif species[0:5] == 'HUMAN':
            if numprobes<500:
                self.sourcefile = ConvergeDataDir+'human_promoter.6.freq'
            else:
                self.sourcefile = file.split('.')[0] + '.' + species + '.bg'
        elif species == 'BAC_ORF':
            self.sourcefile = ConvergeDataDir+'SLR14.2.freq'
        elif species == 'Smik':
            if numprobes<500:
                self.sourcefile = ConvergeDataDir+'Smik.bg'
            else:
                self.sourcefile = file.split('.')[0] + '.' + species + '.bg'
        elif species == 'Spar':
            if numprobes<500:
                self.sourcefile = ConvergeDataDir+'Spar.bg'
            else:
                self.sourcefile = file.split('.')[0] + '.' + species + '.bg'
        elif species == 'Sbay':
            if numprobes<500:
                self.sourcefile = ConvergeDataDir+'Sbay.bg'
            else:
                self.sourcefile = file.split('.')[0] + '.' + species + '.bg'
	elif species == 'Skud':
	    self.sourcefile = ConvergeDataDir+'Skud.bg'
	elif species == 'Scas':
	    self.sourcefile = ConvergeDataDir+'Scas.bg'
	elif species == 'Sklu':
	    self.sourcefile = ConvergeDataDir+'Sklu.bg'
	elif species == 'human':
	    self.sourcefile = ConvergeDataDir+'human.bg'
	elif species == 'mouse':
	    self.sourcefile = ConvergeDataDir+'mouse.bg'
	elif species == 'rat':
	    self.sourcefile = ConvergeDataDir+'rat.bg'
        elif species == 'Ciona':
            self.sourcefile = ConvergeDataDir+'savignyi.bg'
            self.highestorder = 3
        elif species == 'intestinalis':
            self.sourcefile = ConvergeDataDir+'intestinalis.bg'
            self.highestorder = 3
        else:
            print 'EM.MarkovBackground: Unknown species %s, using Yeast'
            self.sourcefile = ConvergeDataDir+'yeast.nc.6.freq'
        self.species = species
        self.D  = {}
        self.F  = {}         #Frequencies
        self.CP = {}         #log2(Conditional Probabilities)  CP['ACTG'] = p( G | ACT ) 
        self.nmers_by_size = map(lambda x:[],range(0,10))
        if seqs:
            self.highestorder = 0
            self.freq_from_seqs(seqs)
        else:
            self.freq_from_file()
        self.compute_conditional()
        self.totD = {}
        #self.totD = shelve.open('MarkovBackgroundTotals'+`self.highestorder`+'.shelve')

    def zeroth(self):
        D = {}
        for L in ['A', 'C', 'G', 'T', '-']:
            D[L] = self.F[L]
        return(D)

    def permute(self,letters, depth, seqs=[''],curdepth=0):
        newseqs = []
        for seq in seqs:
            for letter in letters:
                newseqs.append(seq + letter)
        if depth > curdepth:
            return(self.permute(letters,depth,newseqs,curdepth + 1))
        else:
            return(seqs)

    def compute_conditional(self):
        #print self.permute(['A','C','T','G'], 3)
        #print '>>>'
        for lett in ['A','C','T','G','-']:
            self.CP[lett]= math.log(self.F[lett]) / math.log(2)
        for depth in range(2,self.highestorder+1):
            for leading_seq in self.nmers_by_size[depth-1]:
                sub_total = self.F[leading_seq]
                sub_tot   = 0
                for trailing_lett in ['A','C','T','G','-']:
                    sub_key = "%s%s"%(leading_seq,trailing_lett)
                    if self.F.has_key(sub_key):
                        sub_tot = sub_tot + self.F[sub_key]
                for trailing_lett in ['A','C','T','G','-']:
                    sub_key = "%s%s"%(leading_seq,trailing_lett)
                    if self.F.has_key(sub_key):
                        self.CP[sub_key] = math.log(self.F[sub_key]/sub_tot) / math.log(2)
                    #else:
                    #    self.CP[sub_key] = math.log(.01*pow(10.,-depth)) / math.log(2)
                    #    self.F[sub_key]  = 0
                        #print "%10s %8.4f"%(sub_key,pow(2,self.CP[sub_key]))
                #print '-'*20, (sub_total-sub_tot)/sub_tot*100,'%'
    def freq_from_seqs(self,seqs):
        for depth in range(1,6):
            nmersT = ConvergeMotifTools.top_nmers(depth, seqs, "TUPLES")
            self.nmers_by_size[depth] = map(lambda x:x[0],nmersT)
            total = 0
            for nmer,count in nmersT:
                total = total + count
            for nmer,count in nmersT:
                rc = ConvergeMotifTools.revcomplement(nmer)
                if nmer == rc:                       #correct top_nmers 
                    f   = float(count)/total         #palindrome count
                else:
                    f   = float(count)/total/2
                self.F[nmer] = f
                self.F[rc]   = f
        for depth in range(0):                       #For debugging
            total = 0
            for k in self.F.keys():
                if len(k) == depth:
                    total = total + self.F[k]
                    print k, self.F[k]
            print depth,total
                
    def study_seqs(self,seqs):
        for depth in range(1,6):
            nmersT = ConvergeMotifTools.top_nmers(depth, seqs, "TUPLES")
            total = 0
            for nmer,count in nmersT:
                total = total + count
                rc = ConvergeMotifTools.revcomplement(nmer)
            for nmer,count in nmersT:
                f   = math.log(float(count)/total)/math.log(2)
                f_2 = math.log(0.5 * float(count)/total)/math.log(2)
                rc = ConvergeMotifTools.revcomplement(nmer)
                if rc != nmer:
                    self.D[nmer] = f_2
                    self.D[rc]   = f_2
                else:
                    self.D[nmer] = f
        for depth in range(0):
            total = 0
            for k in self.D.keys():
                if len(k) == depth:
                    total = total + pow(2,self.D[k])
                    print k, pow(2,self.D[k])
            print depth,total
        self.highestorder = 5

    def freq_from_file(self):
        print "# Loading High-order Markov background model from:\n#\t\t%s"%self.sourcefile
        FH = open(self.sourcefile,'r')
        lines = FH.readlines()
        FH.close()
        for line in lines:
            if line[0] == '#':  #Comment line
                continue
            (key,freq) = line.split()
            self.F[key.upper()] = float(freq)
        self.highestorder = -1 + max(map(len,self.F.keys()))
        for depth in range(1,self.highestorder):
            a = []
            for key in self.F.keys():
                if len(key) == depth:
                    a.append(key)
            #print depth,len(a),self.nmers_by_size
            self.nmers_by_size[depth] = a
    def logbackground(self,s):
	seq = self.replace_N(s)
	_T = 0
        if  len(seq) == 0:
		_T = 0
        elif self.totD.has_key(seq):   #Memo-ize frequent lookups
            _T = self.totD[seq]
        else:
            for i in range(len(seq)):
                start  = max(0,i+1-self.highestorder)
                subseq = seq[start:i+1]
                #print start,i+1-self.highestorder,subseq,self.D[subseq]
                #print "subSEQ: ",start,i+1,subseq
                _T = _T + self.CP[subseq]
            self.totD[seq] = _T
        return(_T)
    def replace_N(self,seq):
	  s = seq
          tot = self.CP['A']+self.CP['C']+self.CP['G']+self.CP['T']
	  r = [self.CP['A']/tot, (self.CP['A']+self.CP['C'])/tot, (self.CP['A']+self.CP['C']+self.CP['G'])/tot, tot]
	  i = s.find('N')
	  while i != -1:
	  	rnum = random.random()
	  	if rnum<r[0]:
			s = s[:i]+'A'+s[i+1:]
	  	elif rnum<r[1]:
			s = s[:i]+'C'+s[i+1:]
	  	elif rnum<r[2]:
			s = s[:i]+'G'+s[i+1:]
		else:
			s = s[:i]+'T'+s[i+1:]
		i = s.find('N')
	  return(s)
    def _keys(self):
        return self.D.keys()
    def ___getitem__(self,key):
        return self.D[key]
    def _has_key(self):
        return D.has_key(key)

class Probe:
    '''
    Probe object: Extension of string class containing
    extra information, including logP and other information
    '''
    def __init__(self,seqs=''):
        global MarkBG
        global theMarkovBackground

	self.seqs = []
        self.labels = []
	self.c_intA = []
        self.genomes = []
        self.info   = ''
        self.logP   = 0
	self.seqbgs = None
        self.c_wmerbgs= None
        intseq = ''

        genomes = '00000000000'
        yeast = {'Scer': 0, 'Spar': 1, 'Smik': 2, 'Sbay': 3}
        #'Skud': 4, 'Scas': 5, 'Sklu': 6}
        human = {'human': 0, 'mouse': 1, 'rat': 2}
        ciona = {'savignyi': 0, 'intestinalis': 1}
        
	if (len(seqs)==0):
            return
        for i in range(7):
            for seq in seqs:
                if (yeast[seq[0]]==i):
                    self.labels.append(seq[0])
                    s = seq[1].upper()
                    s = MarkBG[seq[0]].replace_N(s)
                    self.seqs.append(s)
        #self.numgenomes = len(self.seqs)
        for seq in self.seqs:
            intseq = intseq + seq
        self.c_intA = CVGsupport.seq2int(intseq)
        self.compute_logP()

        if theMarkovBackground.species == 'YEAST':
            for label in self.labels:
                genomes = genomes[:yeast[label]] + '1' + genomes[(yeast[label]+1):]
        elif theMarkovBackground.species == 'human':
            for label in self.labels:
                genomes = genomes[:human[label]] + '1' + genomes[(human[label]+1):]
        elif theMarkovBackground.species == 'Ciona':
            for label in self.labels:
                genomes = genomes[:ciona[label]] + '1' + genomes[(ciona[label]+1):]
        self.genomes = CVGsupport.seq2int(genomes)
        
    def compute_logP(self):
        global MarkBG
	seqsbgs = []
	self.logP = 0
	for seq,label in zip(self.seqs,self.labels):
		#s = seq.replace('-', '')
		bg = MarkBG[label].logbackground(seq)
		seqsbgs.append(bg)
        	self.logP  = self.logP + bg
	self.seqsbgs = CVGsupport.list2double(seqsbgs)
    def __cmp__(self, other):
        if type(other) == type(self):
            seq = other.seq
        else:
            seq = other
        return cmp(self.seq,seq)
    def __repr__(self):
        return self.seqs
    def __len__(self):
        if len(self.seqs):
            return len(self.seqs[0])
        else:
            return 0
    def __getitem__(self,n):
        return self.seqs[n]
    def __getslice__(self,i,j):
        return self.seqs[i:j]
    def __hash__(self):
        return(hash(self.seqs[0]))
    def __commented_del__(self):
        print "PROBE: del (len %d)"%len(self.seqs)
        

class Wmer:
    '''
    This class supports an extensible way to store information about a W-mer,
    including such things as source sequence, flanking residues, count, etc.
    '''
    def __init__(self, seq, count=1):
        self.seq    = seq
        self.rc     = ConvergeMotifTools.revcomplement(seq)
        self.count  = count
        self.src    = ''
        self.srcQ   = 0
        self.lflank = ''
        self.rflank = ''
        self.logQtot= 0
        self.bgprob = 0
        self.score  = 0
        self.Z      = 0
    def __cmp__(self, other):
        if type(other) == type(self):
            seq = other.seq
        else:
            seq = other
        #print seq,'vs.',self.seq, self.rc
        if (seq == self.seq or seq == self.rc):
            return(0)
        else:
            return(1)
    def __repr__(self):
        return self.seq + "(x%d)"%self.count
    def __len__(self):
        return len(self.seq)
    def __getitem__(self,n):
        return self.seq[n]
    def __getslice__(self,i,j):
        return self.seq[i:j]
    def __hash__(self):
        return(hash(self.seq))
    
def loadMarkovBackground(file,numprobes,species='YEAST'):
    global theMarkovBackground
    global MarkBG
    theMarkovBackground = MarkovBackground(file,numprobes,species)
    if species == 'YEAST':
        MarkBG = {'Scer': theMarkovBackground, 'Spar': MarkovBackground(file,numprobes,'Spar'), 'Smik': MarkovBackground(file,numprobes,'Smik'), 'Sbay': MarkovBackground(file,numprobes,'Sbay')}
        #, 'Skud': MarkovBackground(file,numprobes,'Skud'), 'Scas': MarkovBackground(file,numprobes,'Scas'), 'Sklu': MarkovBackground(file,numprobes,'Sklu')}
    elif species == 'human':
	MarkBG = {'human': theMarkovBackground, 'mouse': MarkovBackground('mouse'), 'rat': MarkovBackground('rat')}
    elif species == 'Ciona':
        MarkBG = {'savignyi': theMarkovBackground, 'intestinalis': MarkovBackground(file,numprobes,'intestinalis')}
        
def permute(depth, letters=['A','C','G','T','-'], seqs=[''],curdepth=0):
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(permute(depth,letters,newseqs,curdepth + 1))
    else:
        return(seqs)

def fasta2dict(filename):
    return fasta2seqs(filename,'WANT DICT')

def fasta2seqs(filename,want_dict=''):
    '''
    Quick and dirty Fasta Loader
    '''
    seqs  = []
    seq   = ''
    key   = ''
    keys  = []
    FH    = open(filename,'r')
    lines = FH.readlines()
    FH.close()
    for line in lines:
        if line[0] == '>':
            if key: oldkey = key
            key = line[1:].split()[0]
            if seq:
                seqs.append(seq)
                keys.append(oldkey)
                seq = ''
        else:
            seq = seq + line.strip().upper()
    if seq:
        seqs.append(seq)
        keys.append(key)
    if want_dict:
        D = {}
        for key,seq in zip(keys,seqs):
            D[key]=seq
        return(D)
    else:
        return(seqs)

def main():
    if len(sys.argv) < 2:
        print "Usage: %s <fasta_file> [width] [-random_background]"%(re.sub('^.*/','',sys.argv[0]))
        sys.exit(1)
    filename = sys.argv[1]

    if len(sys.argv) >2:
        width = int(sys.argv[2])
    else: width = 5

    if not ('-random_background' in sys.argv or '-nomarkov' in sys.argv):
        loadMarkovBackground()

    algorithm = ''

    good_count = 2   #Default: Take the top 2
    good_s     = []  #Initialize seq array
    for tok,i in zip(sys.argv,xrange(10)):
        if tok == '-top':
            good_count = int(sys.argv[i+1])
        if tok == '-prior':
            good_s.append(sys.argv[i+1])
        if tok == '-greedy':
            algorithm = "GREEDY"
            
    seqs = fasta2seqs(filename)
    good_s.extend(seqs[0:min(good_count,len(seqs))])
    other_s = seqs[min(good_count,len(seqs)):]

    #MD = MDscan(good_s,other_s,width,"VERBOSE")
    if algorithm == "GREEDY":
        MD = MDscan(good_s,other_s,width,"VERBOSE")
    else:
        MD = EM(good_s,other_s,width,"VERBOSE")

    for C,i in zip(MD.candidates,xrange(1,1000)):
        C.pssm.maxscore = -100
        print "Motif %3d %s"%(i,C)
        for seq in C.pssm.seqs:
            print seq,i,C.pssm.scan(seq)[2][0]
        print '*'*len(seq)
        print "MAP Score: %f"%C.MAP


def log2_sum(logx, logy):
    #ans = logx + math.log(1 + math.exp(logy - logx))
    if (logx > logy):
        (loga,logb) = (logx,logy)
    else:
        (loga,logb) = (logy,logx)
    if  logb < -1e100:
        ans = loga
    else:
        ans = loga + math.log(1 + math.pow(2,logb - loga))/math.log(2)
    return(ans)

log2_sum = CVGsupport.log2_sum    

def tamofile2motifs(filename):
    FID = open(filename,'r')
    lines = FID.readlines()
    FID.close()
    motifs = []
    for i in range(len(lines)):
        if lines[i][0:10] == 'Log-odds matrix'[0:10]:
            m = ConvergeMotifTools.Motif()
            w = len(lines[i+1].split())-1
            for pos in range(w):
                m.ll.append({})
            m.width = w
            for j in range(0,4):
                toks = lines[i+j+2].split()
                L = toks[0][1]
                for pos in range(w):
                    m.ll[pos][L] = float(toks[pos+1])
            m._compute_ambig_ll()
            m._compute_oneletter()
            m._maxscore()
            motifs.append(m)
    return(motifs)

if __name__ == '__main__':
    main()
    #profile.run('main()','md.prof')
    #stats=pstats.Stats('md.prof')
    #stats.sort_stats('cumulative').print_stats()

    #main()
