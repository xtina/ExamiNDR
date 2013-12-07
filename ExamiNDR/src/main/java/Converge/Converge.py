
import sys, re, os, math, time, string, tempfile, pickle
sys.path.append('/home/romer/romerweb/ConvergeSupport')

import ConvergeMotifTools  #Motif, top_nmers, oneletter, etc...
import Cvg as EM
import infoana     #Load .info files for easy manipulation
import Validate    #Check answers agains Transfac
from TAMO.seq import Fasta
from TAMO.util import Arith
import cluster3 as cluster
import profile
from forKenzie import KenzieSequences

if __name__ == '__main__':
    from ConvergeMotifMetrics import ProbeSet, gapped_motifs

import profile
import pstats
from count2motif import tfmats

import gc
if __name__ == '__main__':
    gc.disable()
    print '#|'
flush = sys.stdout.flush

def main():
    if len(sys.argv) < 3:
        print "Usage: %s <fasta_file> [width = None ] [options]"%(re.sub('^.*/','',sys.argv[0]))
        print "Options include:"
        print "                  -valid  <tf_name> Check answers against Transfac"
        print " EM Parameters:"
        print "                  -beta    [0.01]   Beta for pseudocounts"
        print "                  -seedbeta[0.02]   Beta for pseudocounts for seeds from text"
        print "                  -gamma   [0.2]    Gamma (fraction of sequences)"
        print "                  -delta   [0.001]  Convergence criteria"
        print " "
        print " Seeds (not actually proper priors)"
        print "                  -prior            Seqences or motifs for seeds (may be repeated)"
        print "                  -top N   [0]      Include w-mers in top N probes"
        print "                  -gap    string    sample gapped motifs"
        print "                  -TF               Seed with (all) TRANSFAC PSSMs (buggy)"
        print "                  -info <file.info> for structural priors"
        print "                  -pad              add NN..NN to seed"
        print " "
        print " Genome / Background model "
        print "                  -human (250,1000) Use Human Background model"
        print "                  -Y2K, -Y5C        Use Yeast Upstream Intergenic regions (2000, 500)"
        print "                  -B                Use Bacterial Orfs"
        print " " 
        print "Examples:"
        print " %s t.fsa 5 -prior GGGTA -prior AAAAAC "%(sys.argv[0].split('/')[-1])
        print "   will start an EM with 3 seeds: GGGTA, AAAAA, and AAAAC"
        print 
        print " %s t.fsa 5 -info CUP9.info -gamma 0.5 "%(sys.argv[0].split('/')[-1])
        print "   will start an EM with Enriched seeds in CUP9.info, with"
        print "   Gamma expectation of 50% of all probes"
        print 
        print " %s t.fsa -prior MCM1_5.tamo:0 "%(sys.argv[0].split('/')[-1])
        print "   will start an EM with 0th motif of the file MCM1_5.tamo"
        print "   as a seed"
        print 
        sys.exit(1)
    fastafile = sys.argv[1] 

    #Echo the command line
    print "#" + ' '.join(map(lambda x: re.sub(' ','\ ',x), sys.argv))

    if sys.argv[2].isdigit():
        width = sys.argv[2]
    else: width = None
    
    algorithm = ''
    beta      = ''
    seedbeta  = ''
    cbeta     = ''
    deltamin  = ''
    gamma     = 0.2
    infofile  = ''
    seedmodels= []
    species   = 'YEAST'
    valid_tfs = []
    gapped_syl= None
    gapflank  = 0
    gapweight = 0.2
    enrichfact= 0.7
    pmax      = 0  #False
    TFSEEDS   = 0
    TFMids    = []
    pad       = None
    padlen    = 0
    thetas    = []
    seed_count = 0   #Default: Take the top 0
    seed_s     = []  #Initialize seq array
    sp_seed   = 0
    
    '''Parse command-line arguments'''
    for tok,i in zip(sys.argv,xrange(len(sys.argv))):
        if   tok == '-top'   :   seed_count = int(sys.argv[i+1])
        elif tok == '-greedy':   algorithm  = "GREEDY"
        elif tok == '-prior' :   seed_s.append(sys.argv[i+1])
        elif tok == 'sp'     :   sp_seed    = 1
        elif tok == '-beta'  :   beta       = float(sys.argv[i+1])
        elif tok == '-beta'  :   seedbeta   = float(sys.argv[i+1])
        elif tok == '-cbeta' :   cbeta      = float(sys.argv[i+1])
	elif tok == '-thetas':   
		for j in range(int(sys.argv[i+1])):
			thetas.append(float(sys.argv[i+j+2]))
        elif tok == '-gamma' :   gamma      = float(sys.argv[i+1])
        elif tok == '-delta' :   deltamin   = float(sys.argv[i+1])
        elif tok == '-info'  :   infofile   = sys.argv[i+1]
        elif tok == '-valid' :   valid_tfs.append(sys.argv[i+1])
        elif tok == '-w'     :   width      = sys.argv[i+1]
        elif tok == '-width' :   width      = sys.argv[i+1]
        elif tok == '-gap'   :   gapped_syl = sys.argv[i+1]
        elif tok == '-gapflank' :gapflank   = int(sys.argv[i+1])
        elif tok == '-gapweight':gapweight  = float(sys.argv[i+1])
        elif tok == '-enrichfact':enrichfact= float(sys.argv[i+1])
        elif tok == '-pmax'  :   pmax       = 1
        elif tok == '-Y2K'   :   species    = "YEAST_2000_UP"
        elif tok == '-Y5C'   :   species    = "YEAST_500_UP"
        elif tok == '-B'     :   species    = "BAC_ORF"
        elif tok == '-Ch22'  :   species    = "Ch22"
        elif tok == '-genome':   species    = sys.argv[i+1]
        elif tok == '-pad'   :
            pad        = sys.argv[i+1]
            padlen     = sys.argv[i+2]     
        elif tok == '-TF'    :
            TFSEEDS = 1
            for j in range(i+1,len(sys.argv)):
                if re.match('M0',sys.argv[j]):
                    TFMids.append(sys.argv[j])
                else:
                    break
        elif tok == '-human' :
            _s = ''
            if sys.argv[i+1].isdigit(): _s = '_'+sys.argv[i+1]
            else:                       _s = ''
            species    = 'HUMAN'+_s

    seqs = []
    fsaD     = Fasta.load(fastafile)
    probes = fsaD.keys()
    '''
    for probeid in fsaD.keys():
        seqs.append  (fsaD    [probeid])
    '''
    numprobes = len(probes)
    #print "numprobes: %i"%numprobes
    if not ('-random_background' in sys.argv or '-nomarkov' in sys.argv):
        EM.loadMarkovBackground(fastafile,numprobes,species)

    #seqs     = EM.fasta2seqs(fastafile)
    all_seqs = seqs
    seed_s.extend(seqs[0:min(seed_count,len(seqs))])
    #not necessary --- seed_c.extend(c_seqs[0:min(seed_count,len(seqs))])

    if infofile and width=='info':
        width = info2width(infofile)
    elif width != None:
        width = int(width)

    #Alternate source of seeds
    if infofile:
        if 1 or width:
            seedmodels.extend(info2seeds(width,infofile,fastafile,species))
        else:
            print 'Error: need to specify motif width w/ .info file'
    
    #Any -prior pointers to motifs in other files?
    (seed_s, motifs) = parse_priors(seed_s)
    seedmodels.extend(motifs)

    #Should we get seeds from TRANSFAC?
    if TFSEEDS:
        tf = []
        D  = tfmats()
        if not TFMids:
            keys = D.keys()
        else:
            keys = []
            for TFMid in TFMids:
                for key in D.keys():
                    if key[0:6] == TFMid:
                        keys.append(key)
                        break
        for key in keys:
            m = D[key]
            m.seednum = int(re.sub('M0*','',key.split()[0]))
            m.seedtxt = '%-24s %s'%(m,key)
            tf.append(m)
        tf.sort(lambda x,y: cmp(x.seednum,y.seednum))
        seedmodels.extend(tf)
        #seedmodels.append(tf[33])

    if gapped_syl:
        gapped_priors = gapped_motifs(gapped_syl)
        gapped_priors = map(lambda x:'N'+x+'N', gapped_priors)
        seed_s.extend(gapped_priors)

    if pad:
        print '# Padding models with NN-m-NN'
        newmodels = []
        for m in seedmodels:
            newmodels.append(m[-2,m.width+2])
        seedmodels = newmodels

    '''
    Set everything up and GO!!
    '''
    global theEM
    theEM = EM.EM(seed_s,[],[],width,"VERBOSE")
    if beta:     theEM.beta     = beta
    if cbeta:    theEM.cbeta    = cbeta
    if deltamin: theEM.deltamin = deltamin
    if seedbeta: theEM.seedbeta = seedbeta
    if thetas:   theEM.thetas = thetas
    theEM.param['gamma']        = gamma
    theEM.probeids.extend(probes)
    theEM.seqs.extend(all_seqs)
    #theEM.cons_seqs.extend(c_seqs)
    theEM.models    = seedmodels
    theEM.gapflank  = gapflank
    theEM.gapweight = gapweight
    theEM.report()
    theEM.EM_Cstart()    #GO!!

    #print "#Sorting candidates"
    #sys.stdout.flush()
    #EM.candidates.sort(lambda x,y: cmp(y.MAP,x.MAP))

    #sys.exit(0)
    
    '''
    Compute some metrics
    '''
    print "#Loading Genome %s"%species ; sys.stdout.flush()
    if species == 'human':
	Genome = ProbeSet('HUMAN',enrichfact)
    else:
    	Genome = ProbeSet(species,enrichfact)
    ids    = Genome.ids_from_file(fastafile)

    #fsaDict = Fasta.load(fastafile)
    #probes = fsaDict.keys()
    #cons_pickle = fastafile.split('.')[0] + '.cpickle'
    for C in theEM.candidates:
        #p_cons = conservation_pvalue(C.pssm,probes,fsaDict,ConsDict,4)
        #print p_cons
        if not pmax:
            w_dict = Genome.w_dict
            for key,i in zip(w_dict.keys(),range(len(C.pssm.thetas))):
                w_dict[key] = C.pssm.thetas[i]
            Genome.w_dict = w_dict
            C.pssm.pvalue = Genome.p_value(C.pssm,ids,'verbose')
            #print "P-VAL: %f"%(Genome.p_value(C.pssm,ids,'verbose')*p_cons)
            C.pssm.church = Genome.church(C.pssm,ids)
        else:
            (p,frac) = Genome.best_p_value(C.pssm,ids)
            C.pssm.pvalue    = p
            C.pssm.threshold = frac * C.pssm.maxscore
            print "Bests:",p,frac

    for valid_tf in valid_tfs:
        C.pssm.valid = Validate.validate(C.pssm,valid_tf,'Verbose',"Want Tuple")
    
    '''
    Print out all motifs (sorted by Enrichment) in an AlignACE-like form
    '''

    theEM.candidates.sort(lambda x,y: cmp(x.pssm.pvalue,y.pssm.pvalue))
    for C,i in zip(theEM.candidates,range(len(theEM.candidates))):
        C.pssm.maxscore = -100  #May have side effects.  Recompute when done
        if C.pssm.valid:
            _t = C.pssm.valid
            if not _t[0]:
                vstring = "(--- %8.4f %8.4f %s)"%(_t[1],_t[2],_t[3])
            else:
                vstring = "(HIT %8.4f %8.4f %s)"%(_t[1],_t[2],_t[3])
        else:
            vstring = ''
        C.pssm._maxscore()     #Recomputed
        print "Log-odds matrix for Motif %3d %s"%(i,C)
        C.pssm._print_ll()
        print "Sequence Logo"
        C.pssm._print_bits()
        flush()
        #print '# %3d matching sequences at 90%%'%len(C.pssm.bestseqs(C.pssm.maxscore * 0.9))
        flush()
        m = C.pssm
        if not m.__dict__.has_key('gamma'):  m.gamma = None #Kludge to deal w/ old shelves
        if m.seedtxt:     print "Seed: %3d %s"%(i,m.seedtxt)
        if m.source:      print "Source: ",m.source
        if m.gamma:       print "Gamma: %7.5f"%m.gamma
        if m.threshold:   print "Threshold: %5.2f"%m.threshold
        if m.thetas != []:
            tstr = "thetas:"
            for theta in m.thetas:
                tstr = tstr + " " + str(theta)
            print tstr                                    
        #if C.pssm.seedtxt:
        #    print 'Seed  %3d %-25s'%(i,C.pssm.seedtxt)
        if C.pssm.church != None: vstring = 'ch: %5.2f  %s'%(
            math.fabs(math.log(C.pssm.church)/math.log(10)), vstring)
        print "Motif %3d %-25s  nlog(p): %6.3f  %s"%(i,C,-math.log(C.pssm.pvalue)/math.log(10),vstring)
        if C.pssm.threshold:
            print "Threshold: %6.3f  %4.1f%%"%(
                C.pssm.threshold, 100.0*C.pssm.threshold/C.pssm.maxscore)
            

        C.pssm.maxscore = -1e100  #May have side effects.  Recompute when done
        for seq in C.wmers:
            print seq,i,C.pssm.scan(seq)[2][0]
        C.pssm._maxscore()      #Recomputed
        print '*'*len(seq)
        print "MAP Score: %f"%C.MAP
        sys.stdout.flush()
    sys.stdout.flush()
    sys.exit(0) #Avoid ridiculous python cleanup times

def info2width(infofile): 
    width = 8
    if re.search('.info$',infofile):
        I    = infoana.Infofile(infofile,'DONT REMOVE QUERY')
        lmax = 0
        for key in I.keys():
            for site in I[key]['bsites']:
                if site.type == 'FP':
                    l = len(site.cleantxt())
                    if l > lmax: lmax = l
        if lmax < 7:
            width = 8
        else:
            width = lmax+1
        if re.search('zf\-C2H2',infofile):     #Usually multimeric 
            width = 9
        if re.search('Zn_clus',infofile):      #Usually dimer with (long?) gap
            width = 18
        if re.search('SRF-TF',infofile):       #Usually dimer
            width = width * 2
        if re.search('HSF_DNA-bind',infofile): #HSF binds as trimer
            width = width * 3                
    elif re.search('.fsa$',infofile):
        seqs = theEM.fasta2seqs(infofile)
        lenmax = max(map(len,seqs))
        width  = lenmax+1
    return width

def info2seeds(N,infofile,probefile,species='YEAST'):
    if species == 'human':
	species = 'HUMAN'
    G    = ProbeSet(species)
    IDs  = G.ids_from_file(probefile)
    Q    = EM.theMarkovBackground.zeroth()
    seqs = []
    
    if re.search('.info$',infofile):
        #I    = infoana.Infofile(infofile,'DONT REMOVE QUERY')
        I    = infoana.Infofile(infofile)
        print "# Loading infofile: %s"%infofile
        print I
        seqs = map(lambda x: 'NNNN%sNNNN'%x, I.bsites2seqs(50.0))
    elif re.search('.fsa$',infofile):
        fsaDict = Fasta.load(infofile)
        probes = fsaDict.keys()
        #sequence_repository = KenzieSequences()
        cons_pickle = infofile.split('.')[0] + '.cpickle'
        try:
            CFH = open(cons_pickle, 'r')
            ConsDict = pickle.load(CFH)
            CFH.close()
        except:
            ConsDict = {}
            for probe in probes:
                seqs = []
                cons = []
                try:
                    seq_list = G.alignments[probe]
                except:
                    continue
                if (seq_list!=[]):
                    cer_seq = seq_list[0][1]
                else:
                    cer_seq = ''
                cer_seq = cer_seq.upper()
                numg = len(seq_list) - 1
                for i in range(1,4):
                    try:
                        seqs.append(seq_list[i][1].upper())
                    except:
                        seqs.append('')
                    cons.append([])
                for position in range(len(cer_seq)):
                    ref = cer_seq[position]
                    for i in range(3):
                        if (seqs[i]==''): continue 
                        if (seqs[i][position]!=ref):
                            cons[i].append(1)
                        else:
                            cons[i].append(0)
                ConsDict[probe] = cons
            CFH = open(cons_pickle, 'w')
            pickle.dump(ConsDict, CFH)
            CFH.close()
       
        for probe in probes:
            superseq = ''
            try:
                seq_list = G.alignments[probe]
            except:
                continue
            for seq in seq_list:
                subseq = seq[1].replace('-','')
                subseq = subseq.replace('.','')
                seqs.append(subseq)
    
    if not N:
        nmers = seqs
    else:
        if (N<11):
            nmers= ConvergeMotifTools.top_nmers(N,seqs)
        else:
            gaplen = N - 2*(N/3)
            gr = ''
            for i in range(gaplen):
                gr = gr + 'N'
            nmers = ConvergeMotifTools.top_nmers(N,seqs,0,'',1)
            gnmers = []
            for nmer in nmers:
                gnmers.append(nmer[0:(N/3)]+gr+nmer[(N/3):2*(N/3)])
            nmers = gnmers
        if len(nmers) > 201: nmers = nmers[0:200]
        
    print "Scoring enrichment of %d nmers from .info file"%len(nmers)

    nmers_scoresT = []
    for nmer in nmers:
        if nmer[0:(N/3)].isalpha():
            p = G.p_value(nmer,IDs,'verbose')
            #if (species=='Ciona'): ng = 2
            #else: ng = 4
            #p_cons = conservation_pvalue(nmer,IDs,fsaDict,ConsDict,ng)
            #if (p_cons<0.1):
            nmers_scoresT.append((nmer,p))
    nmers_scoresT.sort(lambda x,y: cmp(x[1],y[1]))
    #for tup in nmers_scoresT:
    #    print tup
    last = min(20,len(nmers_scoresT))
    models = []
    for i in range(last):
        seq = nmers_scoresT[i][0]
        m = ConvergeMotifTools.Motif('',Q)
        m.compute_from_text(seq,0.1)
        models.append(m)
    return(models)
    
def parse_priors(seeds):
    'Some of these seed might be pointers into a file'
    'Such as ACE2.tamo:4 or ACE2.tamo'
    txt_seeds = []
    mat_seeds = []
    for seed in seeds:
        if re.search('\.(tamo)',seed) or re.search('\.seeds*',seed):
            toks = seed.split(':')
            filename = toks[0]
            del toks[0]
            motifs = tamofile2motifs(filename)
            if toks:  #Are any particular motif number specified?
                print '#Motifs from %s to be used as seeds:'%(filename)
                for tok in toks:
                    if tok.isdigit():
                        m = motifs[int(tok)]
                        mat_seeds.append(m)
                        print '# %3d %s'%(int(tok),m)
            else:     #No specifics, so try them all
                print '#Using all (%d) motifs in %s as seeds:'%(len(motifs),filename)
                mat_seeds.extend(motifs)
                for m,idx in zip(motifs,range(len(motifs))):
                    print '# %3d %s'%(idx,m)
        else:
            txt_seeds.append(seed)
    return(txt_seeds,mat_seeds)

def conservation_pvalue(nmer, IDs, fsaDict, ConsDict, num_alignments):
    width = len(nmer)
    total = []
    unconserved = []
    thetas = []
    for i in range(num_alignments-1):
        total.append(0)
        unconserved.append(0)
        tot_positions = 0
        tot_unconserved = 0
        for ID in IDs:
            try:
                cons = ConsDict[ID]
            except:
                continue
            if (cons[i]==[]): continue
            try:
                tot_positions = tot_positions + len(cons[i])
                tot_unconserved = tot_unconserved + cons[i].count(1)
            except:
                print "cons: %s"%cons
        try:
            thetas.append(float(tot_unconserved)/tot_positions)
        except:
            thetas.append(1.0)
    #print thetas
    for ID in IDs:
        seq = fsaDict[ID]
        seqrc = ConvergeMotifTools.revcomplement(seq)
        try:
            cons = ConsDict[ID]
        except:
            continue
        if (cons==[]):
            continue
        hits = []
        hitsr = []
        if (type(nmer)==type(ConvergeMotifTools.Motif())):
            hits = cluster.matches_old(nmer, seq, 0.7)
            #print len(hits)
        else:
            site_re = re.compile(AmbigToRegExp(nmer))
            hit = site_re.search(seq)
            hitr = site_re.search(seqrc)
            while (hit!=None):
                hits.append(hit.start())
                hit = site_re.search(seq,hit.end())
            while (hitr!=None):
                if (hits.count(len(seq)-hitr.end()-1)==0):
                    hits.append(len(seq) - hitr.end()-1)
                hitr = site_re.search(seqrc,hitr.end())
        for hit in hits:
            for i in range(width):
                for j in range(num_alignments-1):
                    if (cons[j]!=[]):
                        total[j] = total[j] + 1
                        unconserved[j] = unconserved[j] + cons[j][hit+i]
    mean = 0
    var = 0
    uc = 0
    for i in range(num_alignments-1):
        m = total[i]*thetas[i]
        #nq = total[i]*(1-thetas[i])
        #if (min(nq,m)<5): return(0.5)
        mean = mean + m
        var = var + total[i]*thetas[i]*(1-thetas[i])
        uc = uc + unconserved[i]
    stdev = math.sqrt(var)
    if (stdev>0):
        Z = ((uc + 0.5)-mean)/stdev
    else:
        Z = 0
    if (Z>=0):
        p = Arith.lzprob(Z)
    else:
        p = 1.0 - Arith.lzprob(-Z)
    return(p)

def n_k(n,k):
    if ((k==0)|(k==n)):
        answer = 1
    else:
        answer = n_k(n-1,k-1) + n_k(n-1,k)
    return(answer)

def remove_hits(hits,width):
    for i in range(1,len(hits)):
        if ((hits[i]-width)<hits[i-1]):
            hits.remove(hits[i])
            return(1)
    return(0)

def save_motifs(motifs,filename,kmer_count=20):
    _old_stdout = sys.stdout  #Cache REAL stdout
    sys.stdout  = open(filename,'w')
    print_motifs(motifs,kmer_count)
    sys.stdout.close()
    sys.stdout  = _old_stdout
    
def print_motif(motif,kmer_count=20,istart=0):
    print_motifs([motif],kmer_count,istart)

def print_motifs(motifs,kmer_count=20,istart=0):
    i = istart-1
    for m in motifs:
        i = i + 1
        print "Log-odds matrix for Motif %3d %s"%(i,m)
        m._print_ll()
        print "Sequence Logo"
        m._print_bits()
        if not m.__dict__.has_key('gamma'):  m.gamma = None #Kludge to deal w/ old shelves
        if not m.__dict__.has_key('church'):  m.church = None #Kludge to deal w/ old shelves
        if not m.__dict__.has_key('realpvalue'): m.realpvalue = None #Kludge to deal w/ old shelves
        if m.seedtxt:  print "Seed: %3d %s"%(i,m.seedtxt)
        if m.gamma:    print "Gamma: %7.5f"%m.gamma
        if m.evalue != None: print 'Evalue: %6.3e'%m.evalue
        if m.source:   print "Source: ",m.source
        #Motif   0 NGAGGGGGNN (0)            (Bits:   8.24   MAP:   6.53   D:  0.21  0)  Enr: 54.000 
        print "Motif %3d %-25s (Bits: %5.2f  MAP: %5.2f   D: %5.3f  %2d) Enr: %6.3f"%(
            i, m, m.totalbits, m.MAP, m.seeddist, m.seednum, math.fabs(math.log(m.pvalue)/math.log(10.))),
        if m.church != None:  print ' ch: %5.2f'%(
            math.fabs(math.log(m.church)/math.log(10.))),
        if m.realpvalue != None: print ' P: %6.4e'%(m.realpvalue),
        print

        _max = m.maxscore
        m.maxscore = -100
        if kmer_count >= 0:  seqs = m.bogus_kmers(kmer_count)
        else:                seqs = m.seqs
        for seq in seqs:     print seq,i,m.scan(seq)[2][0]
                
        m.maxscore = _max
        print '*'*m.width
        print "MAP Score: %f"%(m.MAP)
        sys.stdout.flush()

def load(filename): return tamofile2motifs(filename)
def tamofile2motifs(filename):
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
            m = ConvergeMotifTools.Motif_from_ll(ll)
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
                motifs[-1].church = math.pow(10,-float(toks[14]))
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
    for i in range(len(motifs)):
        if seedfile: motifs[i].seedfile = seedfile
        seednum = motifs[i].seednum
        if seedD.has_key(seednum):
            motifs[i].seedtxt = seedD[seednum]
    return(motifs)

def AmbigToRegExp(expression):
    expression=expression.replace('R','[RAG]')
    expression=expression.replace('Y','[YTC]')
    expression=expression.replace('W','[WTA]')
    expression=expression.replace('S','[SCG]')
    expression=expression.replace('M','[MAC]')
    expression=expression.replace('K','[KGT]')
    expression=expression.replace('H','[HATC]')
    expression=expression.replace('B','[BGCT]')
    expression=expression.replace('V','[VGAC]')
    expression=expression.replace('D','[DGAT]')
    expression=expression.replace('N','[NATCGRYWSMKHBVD-]')
    return expression
                                                

if __name__ == '__main__':
    if 1:
        main()
    else:
        profile.run('main()','md.prof')
        stats=pstats.Stats('md.prof')
        stats.sort_stats('cumulative').print_stats()
