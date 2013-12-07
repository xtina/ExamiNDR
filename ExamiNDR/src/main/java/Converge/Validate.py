
import sys, re, os, math, time, string, tempfile
from TAMO import MotifTools  #Motif, top_nmers, oneletter, etc...

#########################################################################
# THESE LINES ARE NEEDED TO ACCESS TRANSFAC                             #
#                                                                       #
import Transfac                                                         #
transfacdir = 'Transfac/'                       #
#TF_FACTORS  = Transfac.openTfactorDB(transfacdir + "ALLfactor.pydb")   #
#TF_SITES    = Transfac.openTfSiteDB (transfacdir + "ALLsite.pydb")     #
#SWP2TF      = Transfac.openAuxDB    (transfacdir + "ALLSwissprot.pydb")#
#                                                                       #
#                                                                       #
#########################################################################
TF_FACTORS = None
TF_SITES   = None
SWP2TF     = None

def load_Transfac():
    'Load shelves of Transfac Data'
    global TF_FACTORS
    global TF_SITES
    global SWP2TF
    TF_FACTORS  = Transfac.openTfactorDB(transfacdir + "ALLfactor.pydb")   
    TF_SITES    = Transfac.openTfSiteDB (transfacdir + "ALLsite.pydb")     
    SWP2TF      = Transfac.openAuxDB    (transfacdir + "ALLSwissprot.pydb")
    

def validate(motif, tf_name, verbose='',want_tuple=''):
    'validate(tf_name, motif): Determine if motif is similar to Transfac entries'
    if not TF_FACTORS: load_Transfac
    real_name = get_transfac_id(tf_name)
    
    if real_name:  #It was possible to find a transfac entry
        Tfactor   = TF_FACTORS[real_name]
        Sites     = get_sites(Tfactor)
        if Sites:
            match = check_validity(motif,Sites,verbose,want_tuple)
        else:
            match = None
            print 'Could not find binding sites in Transfac for %s'%tf_name
    else:
        print "Could not find transfac entry for %s"%tf_name
        match     = None
    return(match)

def check_validity(motif,Sites,verbose='',want_tuple=''):
    maxscore  = motif.maxscore
    resultsT  = []
    matchq    = None
    motif._cut_ll(0.0)
    ICOUNT = 0
    for Site in Sites:
        if len(Site.seqs) > 0:
            maxsitelen = 0
            for seq in Site.seqs:
                if len(seq) > maxsitelen: maxsitelen = len(seq)
            for seq in Site.seqs:
                seq = re.sub('^[actgN]', '' , seq)      #Remove leading garbage
                seq = re.sub('[actgN]$', '' , seq)      #Remove tailing garbage
                seq = re.sub('[acgt]'  , 'N', seq)      #Mask   middle  garbage
                if     (maxsitelen >= 6 and len(seq) < 6) or \
                       (not seq.isalpha()): #or re.search('[SWKRMY]',seq):
                    continue
                (matches,endpoints,scores) = motif.scan("NNNNN%sNNNNN"%seq.upper(),-10000.)
                _T = zip(scores,matches)
                #for i in range(len(_T)): print '  %8.3f %s %2d %s %7.2f'%(_T[i][0],_T[i][1],i,"NNNNN%sNNNNN"%seq.upper(),maxscore)
                _T.sort()
                _T.reverse()
                score = _T[0][0]
                match = _T[0][1]
                #print " >>",score,match
                ICOUNT = ICOUNT + 1
                #print resultsT, ICOUNT
                resultsT.append( (score,score/maxscore,seq) )
                if score/maxscore > 0.75:
                    matchq = 1
    #for i in range(len(resultsT)): print resultsT[i]
    resultsT.sort(lambda x,y: cmp(x[0],y[0]) or cmp(len(y[2]),len(x[2])))
    resultsT.reverse()
    if verbose:
        for score,frac,seq in resultsT:
            print "%8.3f  %8.3f  %s"%(score,frac,seq)
    if want_tuple:
        return(matchq,resultsT[0][0],resultsT[0][1],resultsT[0][2])
    else:
        return(matchq)

def get_transfac_id(tf_name):
    if not TF_FACTORS: load_Transfac()
    if TF_FACTORS.has_key(tf_name):
        name = tf_name
    elif SWP2TF.has_key(tf_name):
        name = SWP2TF[tf_name]
    elif SWP2TF.has_key(tf_name + '_YEAST'):
        name = SWP2TF[tf_name + '_YEAST']
    elif SWP2TF.has_key(tf_name[0:2]+tf_name[3:]+'_YEAST'): #STE12 -> ST12
        name = SWP2TF[tf_name[0:2]+tf_name[3:]+"_YEAST"]
    elif SWP2TF.has_key(tf_name + '_HUMAN'):
        name = SWP2TF[tf_name + '_HUMAN']
    elif SWP2TF.has_key(tf_name[0:2]+tf_name[3:]+'_HUMAN'): #STE12 -> ST12
        name = SWP2TF[tf_name[0:2]+tf_name[3:]+"_HUMAN"]
    else:
        name = ''
    return name
        
def get_sites(Tfactor):
    if not TF_FACTORS: load_Transfac()
    sites = []
    for site_id in Tfactor.get_site_ids():
        if TF_SITES.has_key(site_id):
            sites.append(TF_SITES[site_id])
    return(sites)


def test():
    motifs = []
    betalist =  [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0]
    for beta in [1.0]:
        m = MotifTools.Motif()
        m.compute_from_text('GGTTTCAT', beta) #STE12 binding site
        print m
        m._print_ll()
        print "Against Ste12:"
        match  = validate(m,"STE12",'V','T')
        print "Against Fkh2:"
        fmatch = validate(m,"FKH2", 'V','T')
        print beta, match, fmatch


if __name__ == '__main__':
    test()
    
