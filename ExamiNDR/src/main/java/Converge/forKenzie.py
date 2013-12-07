import pickle,glob,re
import ConvergeFasta


record=ConvergeFasta.Record()




class KenzieSequences:
	def __init__(self,ProbeAlignmentPickleFilename):
		#the following pickle excludes any probe that is predicted to amplify > 1 region
		#pickle_file=open('/cluster/fraenkel/ef/Young/04.01.22.dir/ProbeAlignments_040126.pickle')
		#the following pickle excludes any probe that is predicted to amplify > 2 regions
		#pickle_file=open('/cluster/fraenkel/ef/Young/04.01.22.dir/ProbeAlignments_040604.pickle')
		pickle_file=open(ProbeAlignmentPickleFilename)
		self.ProbeDict= pickle.load(pickle_file)
		#self.genomes=['Scer', 'Smik', 'Spar', 'Sbay']
		self.genomes=['Scer', 'Sbay', 'Scas', 'Sklu', 'Skud', 'Smik', 'Spar']

	
	def returnSequencesForProbes(self,probe):
		sequence_list=[]
		msa = self.ProbeDict.get(probe, None)
		if not msa:
			return []
		if not len(msa.get_gapped_seq('Scer')):
				return []
		for genome in self.genomes:
			seq=msa.data.get(genome)
			if not seq:
				continue
			sequence_list.append((genome, seq))
		return sequence_list

			
		

	
	

if __name__=='__main__':
	from Fraenkel.utilities import averageList

	'''
	probes=['iYAL068C-0',\
'iYAL068C-1',\
'iYAL060W',\
'iYBL019W',\
'iYPR201W-0',\
'iYPR201W-1']
	for probe in probes:

		
		k=KenzieSequences()
		print k.returnSequencesForProbes(probe)
	'''
	num_changes={}
	#k=KenzieSequences('/cluster/fraenkel/ef/Kenzie/05.01.04.dir/ManoliProbeAlignments_050105.pickle')
	k=KenzieSequences('/cluster/fraenkel/ef/Kenzie/05.01.04.dir/USCProbeAlignments_050105.pickle')
	for g in k.genomes:
		num_changes[g]=[0,0]
	for probe in k.ProbeDict.keys():
		num_changes_probe=k.ProbeDict[probe].countChanges('Scer',skip_all_gaps=1)
		for g in [ x for x in k.genomes if x in num_changes_probe.keys()]:
			#print probe, g, num_changes_probe[g]
			num_changes[g][0]+=num_changes_probe[g][0]
			num_changes[g][1]+=num_changes_probe[g][1]

	print '%-10s %10s'%('genome','1-mutation rate')
	for g in k.genomes:
			if num_changes[g][1]:
				f=num_changes[g][0]*1.0/num_changes[g][1]
				print '%-10s %10.2f'%(g,1.-f)
			
		

