
import glob, os, sys, math, random
import ConvergeMotifTools

def find_matches(sequences, motif, threshhold, gen_weights = [], metric = 'weighted'):
   locs = []    #motif locations
   weights = []    #genome weightings
   wsum = 0    #sum of weightings (used for normalization purposes)
   m = ConvergeMotifTools.Motif()
   m.compute_from_ll(motif.ll)
   minscore = m.min_kscore    #best possible distance score obtainable
   maxscore = 0    #worst possible overall score obtainable
   width = m.width    #determine motif width
   num_genomes = len(sequences)   #get the number of alignments for this sequence
   num_positions = len(sequences[0])-width    #determine number of possible motif positions

   #setup genome weightings
   if (gen_weights==[]):
      for i in range(num_genomes):
         weights.append(1)
         wsum = wsum + 1
   else:
      for i in range(num_genomes):
         weights.append(gen_weights[i])
         wsum = wsum + gen_weights[i]
   
   #determine best possible score
   for i in range(num_genomes):
      maxscore = maxscore + weights[i]*(1 - minscore)

   mini = 1000
   maxi = 0
   
   #calculate distance score at each position in the alignment
   for i in range(num_positions):
      distf = 0    #forward distance score
      distr = 0    #reverse distance score
      bseq = []
      for j in range(num_genomes):
         seq = sequences[j][i:(i+width)]
         bseq.append(seq)
         distf = distf + weights[j]*(calc_dist(m, [], seq, metric, 'f')-minscore)
         distr = distr + weights[j]*(calc_dist(m, [], seq, metric, 'r')-minscore)
      score = min(distf,distr)/maxscore
      print score
      '''
      if (score<mini):
         mini = score
         best_seq = bseq
      if (score>maxi): maxi = score
      if (score < threshhold): locs.append([i, score])
   print "min: %f max: %f"%(mini,maxi)
   for seq in best_seq:
      print seq
      '''
   return(locs)

def find_matches_old(sequence, motif, threshold):
   loc = []
   width = len(motif)
   cutoff = threshold*motif.maxscore
   best = 0
   #print cutoff
   for i in range(len(sequence)-width+1):
      seq = sequence[i:(i+width)]
      seqrc = rc(seq)
      #if (len(seq)==width): print "OK"
      f = 0
      r = 0
      for j in range(width):
         try:
            f = f + motif.ll[j][seq[j]]
         except:
            for L in 'ACGT':
               f = f + 0.25*motif.ll[j][L]
         try:
            r = r + motif.ll[width-j-1][seqrc[j]]
         except:
            for L in 'ACGT':
               r = r + 0.25*motif.ll[j][L]
      #print f
      if ((f>cutoff)|(r>cutoff)):
         loc.append(i)
   return(loc)

def matches_old(motif, sequence, threshold):
   loc = []
   width = len(motif)
   cutoff = threshold*motif.maxscore
   best = 0
   #print cutoff
   for i in range(len(sequence)-width):
      seq = sequence[i:(i+width)]
      seqrc = rc(seq)
      #if (len(seq)==width): print "OK"
      f = 0
      r = 0
      for j in range(width):
         try:
            f = f + motif.ll[j][seq[j]]
         except:
            for L in 'ACGT':
               f = f + 0.25*motif.ll[j][L]
         try:
            r = r + motif.ll[width-j-1][seqrc[j]]
         except:
            for L in 'ACGT':
               r = r + 0.25*motif.ll[j][L]
      #print f
      #if ((f>cutoff)|(r>cutoff)):
      #   loc.append(i)
      if (f>best):
         best = f
         best_loc = i
      if (r>best):
         best = r
         best_loc = i
   if (best>cutoff): loc.append(best_loc)
   return(loc)
         
def calc_dist(motif1, motif2=[], seq=[], metric = 'weighted', direction = 'both', bg = [0.31, 0.19, 0.19, 0.31]):
   PSSM1 = []
   PSSM2 = []
   KL1 = []	#KL-divergence with background
   KL2 = []
   vec = []
   dvg = 0
   width1 = motif1.width
   if (motif2 != []): width2 = motif2.width
   else: width2 = 0
   maxw = max(width1,width2)
   minw = min(width1,width2)
   minolp = min(maxw/2+1, minw)
   #if (minw==minolp):
   ohng_pen = 1.0/(maxw-minolp)
   m1 = None
   m2 = None
   if (width1>width2):
      m1 = motif1
      m2 = motif2
   else:
      m1 = motif2
      m2 = motif1
      width1 = width2
      width2 = m2.width

   #Set up the motif PSSMs
   for i in range(width1):
        dvg = 0
	vec = []
        for L,j in zip('ACGT', range(4)):
	   prob = math.pow(2,m1.logP[i][L])
	   vec.append(prob)
	   dvg = dvg + prob*m1.ll[i][L]
	PSSM1.append(vec)
	KL1.append(dvg)
   if (m2 != []):
      for i in range(width2):
         dvg = 0
         vec = []
         for L,j in zip ('ACGT', range(4)):
            prob = math.pow(2,m2.logP[i][L])
            vec.append(prob)
            dvg = dvg + prob*m2.ll[i][L]
         PSSM2.append(vec)
         KL2.append(dvg)
   #calculate the distance score for the 2 motif case
   dist = 0.0
   mindist = 1
   #print "width1: %i width2: %i"%(width1, width2)
   if (m2 != []):
      overlap = 1
      sl = 0
      left_ohng = 0
      for i in range((width1+width2)-1):
         usKLsum = 0.0    #unscaled KLsum
         KLsum = 0.0      #scaled KLsum
         olp_KLsum = 0.0  #overlapping KLsum
         dist = 0.0
         if (sl==0):
            #add overhang penalties
            for j in range(overlap,width1):
               dist = dist + ohng_pen*KL1[j]*distance(PSSM1[j],bg)
               KLsum = KLsum + ohng_pen*KL1[j]
               usKLsum = usKLsum + KL1[j]
            for j in range(width2-overlap):
               dist = dist + ohng_pen*KL2[j]*distance(PSSM2[j],bg)
               KLsum = KLsum + ohng_pen*KL2[j]
               usKLsum = usKLsum + KL2[j]
            #add overlap distances
            i1 = 0
            i2 = width2-overlap
            for j in range(overlap):
               dist = dist + max(KL1[i1],KL2[i2])*distance(PSSM1[i1],PSSM2[i2])
               KLsum = KLsum + max(KL1[i1],KL2[i2])
               olp_KLsum = olp_KLsum + max(KL1[i1],KL2[i2])
               usKLsum = usKLsum + max(KL1[i1],KL2[i2])
               i1 = i1 + 1
               i2 = i2 + 1
            overlap = overlap + 1
            if (overlap > width2):
               left_ohng = 1
               if (width1>width2):
                  sl = 1
                  overlap = width2
               else:
                  sl = 2
                  overlap = width2 - 1
         elif (sl==1):
            #add overhang penalties
            for j in range(left_ohng):
               dist = dist + ohng_pen*KL1[j]*distance(PSSM1[j],bg)
               KLsum = KLsum + ohng_pen*KL1[j]
               usKLsum = usKLsum + KL1[j]
            for j in range(left_ohng+overlap, width1):
               dist = dist + ohng_pen*KL1[j]*distance(PSSM1[j],bg)
               KLsum = KLsum + ohng_pen*KL1[j]
               usKLsum = usKLsum + KL1[j]
            #add overlap distances
            i1 = left_ohng
            i2 = 0
            for j in range(overlap):
               dist = dist + max(KL1[i1],KL2[i2])*distance(PSSM1[i1],PSSM2[i2])
               KLsum = KLsum + max(KL1[i1],KL2[i2])
               olp_KLsum = olp_KLsum + max(KL1[i1],KL2[i2])
               usKLsum = usKLsum +  max(KL1[i1],KL2[i2])
               i1 = i1 + 1
               i2 = i2 + 1
            left_ohng = left_ohng + 1
            if ((left_ohng+overlap)>width1):
               overlap = overlap - 1
               sl = 2
         elif (sl==2):
            #add overhang penalties
            for j in range(left_ohng):
               dist = dist + ohng_pen*KL1[j]*distance(PSSM1[j],bg)
               KLsum = KLsum + ohng_pen*KL1[j]
               usKLsum = usKLsum + KL1[j]
            for j in range(overlap,width2):
               dist = dist + ohng_pen*KL2[j]*distance(PSSM2[j],bg)
               KLsum = KLsum + ohng_pen*KL2[j]
               usKLsum = usKLsum + KL2[j]
            #add overlap distances
            i1 = left_ohng
            i2 = 0
            for j in range(overlap):
               dist = dist + max(KL1[i1],KL2[i2])*distance(PSSM1[i1],PSSM2[i2])
               KLsum = KLsum + max(KL1[i1],KL2[i2])
               olp_KLsum = olp_KLsum + max(KL1[i1],KL2[i2])
               usKLsum = usKLsum + max(KL1[i1],KL2[i2])
               i1 = i1 + 1
               i2 = i2 + 1
            left_ohng = left_ohng + 1
            #print "i: %i left overhang: %i"%(i,left_ohng)
            overlap = overlap - 1
         dist = dist/KLsum
         olp_bits = olp_KLsum/usKLsum
         if ((dist < mindist)&(olp_bits>(2.0/3.0))&((sl>0)&((overlap+1)>=minolp)|(sl==0)&(overlap>minolp))):
            bKLs = KLsum
            mindist = dist
            
   else:	#motif distance to sequence
	Ldict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        KLsum = 0
        seqrc = rc(seq)
        distrc = 0
	for index1 in range(width1):
           if (metric=='weighted'):
	      if (seq[index1]=='-'):
	          dist = dist + KL1[index1]
              elif (seq[index1]=='N'):
                  Nval = 0.0
                  for key in Ldict.keys():
                     Nval = Nval + bg[Ldict[key]]*PSSM1[index1][Ldict[key]]
                  dist = dist + KL1[index1]*(1-Nval)
              else:
                  dist = dist + KL1[index1]*(1 - PSSM1[index1][Ldict[seq[index1]]])
              if (seqrc[index1]=='-'):
                 distrc = distrc + KL1[index1]
              elif (seqrc[index1]=='N'):
                  Nval = 0.0
                  for key in Ldict.keys():
                     Nval = Nval + bg[Ldict[key]]*PSSM1[index1][Ldict[key]]
                  distrc = distrc + KL1[index1]*(1-Nval)
              else:
                 distrc = distrc + KL1[index1]*(1-PSSM1[index1][Ldict[seqrc[index1]]])
              KLsum = KLsum + KL1[index1]
           else:
 	      if (seq[index1]=='-'):
	          dist = dist + 1
              elif (seq[index1]=='N'):
                  Nval = 0.0
                  for key in Ldict.keys():
                     Nval = Nval + bg[Ldict[key]]*PSSM1[index1][Ldict[key]]
                  dist = dist + (1-Nval)
              else:
                  dist = dist + (1 - PSSM1[index1][seq[index1]])
              if (seqrc[index1]=='-'):
                 distrc = distrc + 1
              elif (seqrc[index1]=='N'):
                  Nval = 0.0
                  for key in Ldict.keys():
                     Nval = Nval + bg[Ldict[key]]*PSSM1[index1][Ldict[key]]
                  distrc = distrc + (1-Nval)
              else:
                 distrc = distrc + (1-PSSM1[index1][seqrc[index1]])
              KLsum = KLsum + 1
        if (direction=='both'):
           mindist = min(dist,distrc)/KLsum
        elif (direction=='f'):
           mindist = dist/KLsum
        elif (direction=='r'):
           mindist = distrc/KLsum
   #print "best overlap: %i"%besto
   #print "KLsum: %f"%(bKLs)
   return(mindist)

def distance(probs1, probs2):
    exp_mismatch = 1
    for i in range(len(probs1)):
	exp_mismatch = exp_mismatch - probs1[i]*probs2[i]
    return(exp_mismatch)

def rc(seq):
   seqrc = ''
   width = len(seq)
   RCDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
   for i in range(width):
      seqrc = seqrc + RCDict[seq[width-1-i]]
   return(seqrc)
      
#---------------------------------------------------------
#PAMSIL clustering of motifs
#---------------------------------------------------------
def cluster(Motifs, max_clusters = []):
    #generate distance matrix
    num_motifs = len(Motifs)	#number of motifs to cluster
    m_dist = []			#distance matrix between motifs
    labels = []			#cluster label for each motif
    S = []			#silhouette value for each motif
    col = []
    mindif = -0.00001
    best_sum = -100000
    if max_clusters:
       Kmax = max_clusters
    else: Kmax = 5
    #print "Kmax: %i"%Kmax

    for i in range(num_motifs):
	m1 = Motifs[i]
	col = []
	for j in range(num_motifs):
	    m2 = Motifs[j]
	    dist1 = calc_dist(m1, m2)
	    m2 = ConvergeMotifTools.revcompmotif(m2)
	    dist2 = calc_dist(m1, m2)
	    col.append(min(dist1, dist2))
	m_dist.append(col)
        #print col
	S.append(0)
        
    #do PAMSIL clustering for K = 2
    K = 2
    while (K<=Kmax):
	labels = []
	for i in range(K):
	    labels.append([])
        #initialization
	for i in range(num_motifs):
	    labels[random.randint(0,(K-1))].append(i)	    
	#calculate average silhouette
	S = calc_sil(labels, m_dist, S)
	sumS = 0
	for s in S:
	    sumS = sumS + s

	cvgdif = -100000
        new_S = S
	oldsum = sumS
	while (cvgdif < mindif):    
	      #cycle through swaps and find the best one
	      best = (0, 0, 0)
	      bestdif = 0
              temp_labels = labels
              i = 0
	      for label in labels:
		  for motif in label:
		      label.remove(motif)
                      j = 0
		      for l in labels:
			  if (l!=label):
			     l.append(motif)
			     new_S = calc_sil(labels, m_dist, new_S)
                             dif = 0
                             for s in new_S:
                                dif = dif + s
			     dif = dif - sumS
			     if (dif < bestdif):
				best_dif = dif
				best = (i, j, motif)
			     l.remove(motif)
                          j = j + 1
		      label.append(motif)
                  i = i + 1
	      #apply the best swap
              if (best[0] != 0):    #if there is a good swap available
                 label = labels[best[0]]
                 label.remove(best[2])
                 label = labels[best[1]]
                 label.append(best[2])

	         #calculate average silhouette
                 S = calc_sil(labels, m_dist, S)
                 sumS = 0
                 for s in S:
                    sumS = sumS + s
                 cvgdif = sumS - oldsum
              else: cvgdif = 0
              
              #print "K: %i cvgdif = %f"%(K,cvgdif)
              #print labels
	      oldsum = sumS

	if (sumS < best_sum):
	   bestK = K
	   best_sum = sumS
           print "best sum: %f"%best_sum
	   best_S = S
	   best_labels = labels
	
	K = K + 1

	#repeat clustering until K = Kmax
 
    #return clusters with maximum avg silhouette
    return(labels)

#---------------------------------------------
#calculate the silhouette
#---------------------------------------------
def calc_sil(labels, m_dist, _S):
    a = 0
    b = 0
    S = _S
    for label in labels:
	for motif in label:
	    a = 0
	    b = 0
	    bmin = 100000
	    for l in labels:
		if (label==l):
		   for m in l:
		       a = a + m_dist[motif][m]
		   a = a/len(l)
		else:
		   for m in l:
		       b = b + m_dist[motif][m]
                   if (len(l)>0):
                      b = b/len(l)
                   else:
                      b = 100000
		   if (b<bmin): bmin = b
                   b = 0
	    S[motif] = (a-bmin)/(max(a,bmin))
    return(S)
