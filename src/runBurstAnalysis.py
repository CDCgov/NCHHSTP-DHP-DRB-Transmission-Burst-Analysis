#!/usr/bin/python

# THIS PROGRAM DESCRIBES BURSTS OF HIV TRANSMISSION ACROSS A CLOCK-TREE  
import argparse
import re
import os
import sys

"""
EXAMPLE COMMAND:
>python runBurstAnalysis.py -i Jurisdiction1_Burst_Input.txt -s 2014 -e 2016 -b 3 -d 2019

"""

"""
REQUIRED DEPENDENCIES
ete3 from http://etetoolkit.org/download/
FastTree from http://www.microbesonline.org/fasttree/#Install
TreeTime from https://treetime.readthedocs.io/en/latest/installation.html
"""

def getArgs():
	parser = argparse.ArgumentParser(description='ALERT: Must open Anaconda ETE via "export PATH=~/anaconda_ete/bin:$PATH" or "conda activate ete3" command in Terminal before executing program')
	parser.add_argument('-i', '--input-data-file', help='tab-delimited data file', required=True)
	parser.add_argument('-s', '--epoch-start', help='year epoch begins', required=True)
	parser.add_argument('-e', '--epoch-end', help='year epoch ends', required=True)
	parser.add_argument('-b', '--burst-size', help='minumum transmission events in burst', required=True)
	parser.add_argument('-d', '--end-year', help='year after which diagnoses are ignored', required=True)
	parser.add_argument('--skip-alignment', help='skip sequence alignment step (must have been completed previously)', required=False, action='store_true')
	parser.add_argument('--skip-tree', help='skip constructing phylogenetic tree (must have been completed previously)', required=False, action='store_true')
	parser.add_argument('--skip-clock', help='skip molecular clock inference (must have been completed previously)', required=False, action='store_true')
	args = parser.parse_args()
	return args

arg = getArgs()

epochStart = int(arg.epoch_start)
epochEnd = int(arg.epoch_end)
events = int(arg.burst_size)
upperDxLimit = int(arg.end_year)
dataFile = arg.input_data_file
rootPath = str(arg.input_data_file.split('.')[0])+'.EpochStart'+str(epochStart)+'.EpochEnd'+str(epochEnd)+'.BurstSize'+str(events)+'.EndYear'+str(upperDxLimit)

from ete3 import Tree

"""
_hiv_aids_dx_dt_num: YYYYMMDD; Missing <-YYYY0000
sample_dt: YYYYMMDD
"""


#HXB2 SEQUENCE INCLUDES FULL PR/RT/IN SEQUENCE WITH NO NNNs

HXB2seq = '>HXB2\n\
CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTA\
GATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTAT\
CAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACA\
TAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCA\
GGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGA\
AGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAAT\
TAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAA\
AAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACC\
TAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTA\
GCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGAC\
TTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCA\
GAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCT\
GGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAA\
CTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAA\
AGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAA\
TTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAG\
GCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAAC\
ATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAG\
AGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACT\
AATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTC\
GGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCA\
ATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTA\
GATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAA\
TTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAG\
AAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTT\
CATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAG\
ATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGC\
AGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGA\
GATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGC\
AGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTT\
ATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGT\
GACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACA\
GGATGAGGAT\n'


"""
DEFINE FUNCTIONS
"""

def makeMetaDict(dataFile,rootPath,epochStart):
	metaDict = {}
	for line in open(dataFile):
		p = line.split('\t')
		if 'mpv_uid' in line:
			varDict = {}
			for i, v in enumerate(line.rstrip().split('\t')):
				if v == 'mpv_uid': uidPos = i
				if v == '_hiv_aids_dxdt_num': dxdtPos = i
				if v == 'sample_dt': sampledtPos = i
				if v == 'clean_seq': seqPos = i
		else:
			dxdt = int(str(p[dxdtPos])[0:4])
			sampledt = p[sampledtPos]
			if dxdt >= epochStart and sampledt.isdecimal() == True and len(re.sub('[N-]','',p[seqPos])) > 500 and dxdt < upperDxLimit + 1:
				metaDict[p[uidPos]] = {'_hiv_aids_dx_dt_num':p[dxdtPos],'sample_dt':p[sampledtPos],'clean_seq':p[seqPos]}

	return metaDict



def alignSeqs(rootPath,HXB2seq):
	# WRITE OUT UNALIGNED FASTA FILE
	seqFile = open(rootPath+'.unaligned.fa','w')

	for uid in metaDict: seqFile.write('>'+uid+'\n'+metaDict[uid]['clean_seq']+'\n')

	seqFile.close()

	# WRITE OUT HXB2 REFERENCE FILE
	hxb2File = open(rootPath+'.HXB2ref.fa','w')
	hxb2File.write(HXB2seq)
	hxb2File.close()

	# ALIGN SEQUENCES
	rootPath = re.sub(' ','\ ',rootPath)
	os.system('bealign -r '+rootPath+'.HXB2ref.fa -m HIV_BETWEEN_F -R '+rootPath+'.unaligned.fa '+rootPath+'.bam')
	os.system('bam2msa '+rootPath+'.bam '+rootPath+'.FullPol.fa')

	prrtDict = {}
	for line in open(rootPath+'.FullPol.fa'):
		if '>' in line:
			uid = line.rstrip().strip('>')
			prrtDict[uid] = ''

		else: prrtDict[uid] += line.rstrip()

	prrtFile = open(rootPath+'.fa','w')

	for uid in metaDict: prrtFile.write('>'+uid+'\n'+prrtDict[uid][0:1497]+'\n')

	prrtFile.close()

	# CLEAN UP DIRECTORY
	os.system('rm '+rootPath+'.HXB2ref.fa')
	os.system('rm '+rootPath+'.bam')
	os.system('rm '+rootPath+'.bam.bai')
	os.system('rm '+rootPath+'.FullPol.fa')
	os.system('rm '+rootPath+'.unaligned.fa')

def makeTree(rootPath):
	# FASTTREE
	rootPath = re.sub(' ','\ ',rootPath)
	os.system('FastTree -nt -gtr -cat 20 < '+rootPath+'.fa > '+rootPath+'.GTRCAT20.tree')

def clockTree(rootPath,epochStart):
	dateFile = open(rootPath+'.Dates.csv','w')
	dateFile.write('name,date\n')

	for uid in metaDict:
		if int(metaDict[uid]['_hiv_aids_dx_dt_num'][0:4]) >= epochStart: dateFile.write(uid+','+metaDict[uid]['sample_dt'][0:4]+'-'+metaDict[uid]['sample_dt'][4:6]+'-'+metaDict[uid]['sample_dt'][6:8]+'\n')

	dateFile.close()

	#treetime
	rootPath = re.sub(' ','\ ',rootPath)
	os.system('treetime --clock-rate 0.00122 --aln '+rootPath+'.fa --tree '+rootPath+'.GTRCAT20.tree --dates '+rootPath+'.Dates.csv --keep-polytomies --clock-filter 0 --outdir '+rootPath+'_TreeTime/')

def readTree(rootPath):
	for line in open(rootPath+'_TreeTime/timetree.nexus','r'):
		if 'Tree tree1=' in line:
			nexTree = line.rstrip()

			break

	nwkTree = re.sub('\[[^]]*\]','',nexTree)
	nwkTree = nwkTree.split('=')[1]

	t = Tree(nwkTree,format=1)

	rootDate = 2000
	youngestTip = 1900

	for line in open(rootPath+'_TreeTime/dates.tsv'):
		if 'date' not in line:
			d = float(line.split('\t')[2].rstrip()) 
			if d < rootDate: rootDate = d

			if d > youngestTip: youngestTip = d

	return [t,rootDate,youngestTip]



def makeAgeDict(t,rootName):
	ageDict = {}

	for node in t.traverse():
		if node != rootName:
			nodeAge = rootDate + t.get_distance(node,rootName)
			ageDict[node] = nodeAge

		else: ageDict[node] = rootDate

	return ageDict



"""
DO BURSTS GIVE RISE TO DISPROPORTIONATELY MORE TRANSMISSION EVENTS IN EPOCH #2 COMPARED WITH ALL OTHER BRANCHES IN THE TREE?
"""

def describeBurstRatio(window,ageDict,rootPath):
	burstDict = {'EventCount':0,'EndEpoch1Nodes':{},'Epoch1Tips':set([]),'Epoch2Nodes':{},'BurstNodes':set([]),'BurstTips':set([])}

	for node in ageDict:
		nodeAge = ageDict[node]
		if node.up != None:

			# IDENTIFY ALL BRANCHES SPANNING EPOCH 1 AND EPOCH 2
			# STORE PARENT NODE OF THESE BRANCHES AND THEIR DESCENDANT NODES
			parent = node.up
			parentAge = ageDict[parent]

			if nodeAge >= epochEnd + 1 and parentAge < epochEnd + 1 and parent not in burstDict['EndEpoch1Nodes']:
				burstDict['EndEpoch1Nodes'][parent] = set([])
				for child in parent.get_children():
					childAge = ageDict[child]
					if childAge >= epochEnd + 1 or child.is_leaf() == True:
						burstDict['EndEpoch1Nodes'][parent].add(child)

			# TIPS THAT ARE SAMPLED DURING EPOCH 1
			if nodeAge >= epochStart + 3 - window and nodeAge < epochEnd + 1 and node.is_leaf() == True:
				burstDict['Epoch1Tips'].add(node)

			# NODES THAT ARISE DURING EPOCH 2
			if nodeAge >= epochEnd + 1 and node.is_leaf() == False:
				burstDict['Epoch2Nodes'][node] = node.children

	for node in ageDict:
		nodeAge = ageDict[node]

		if nodeAge >= epochEnd + 1 - window and nodeAge < epochEnd + 1:
			eventCount = len(node.children) - 1
			burstTips = set([])
			burstNodes = set([node])

			for child in node.iter_descendants():
				childAge = ageDict[child]

				# CHILD IS LESS THAN THE END OF THE EPOCH (I.E., END OF 2016; epochEnd)
				if child.is_leaf() == False and childAge < epochEnd + 1:
					eventCount += len(child.children) - 1
					burstNodes.add(child)

				elif child.is_leaf() == True: burstTips.add(child)

			if eventCount >= events:
				if len(burstTips.intersection(burstDict['BurstTips'])) == 0:
					burstDict['EventCount'] += 1

				burstDict['BurstTips'] = burstTips.union(burstDict['BurstTips'])
				burstDict['BurstNodes'] = burstNodes.union(burstDict['BurstNodes'])

	Epoch1BurstBranch = set([])
	Epoch1NonBurstBranch = set([])
	Epoch2BurstNodeCount = 0
	Epoch2NonBurstNodeCount = 0

	for node in burstDict['EndEpoch1Nodes']:
		tempTips = set([])
		tempNodes = set([])

		for d in node.iter_descendants():
			if d.is_leaf() == True: tempTips.add(d)
			else: tempNodes.add(d)

		burst = len(tempTips.intersection(burstDict['BurstTips'])) > 0

		for D in burstDict['EndEpoch1Nodes'][node]:
			if D not in burstDict['Epoch1Tips']:
				if burst == True: Epoch1BurstBranch.add(D)
				else: Epoch1NonBurstBranch.add(D)

	for tip in burstDict['Epoch1Tips']:
		if tip in burstDict['BurstTips']: Epoch1BurstBranch.add(tip)
		else: Epoch1NonBurstBranch.add(tip)

	for node in burstDict['Epoch2Nodes']:
		tempTips = set([])
		for d in node.iter_descendants():
			if d.is_leaf() == True: tempTips.add(d)

		burst = len(tempTips.intersection(burstDict['BurstTips'])) > 0

		if burst == True: Epoch2BurstNodeCount += len(node.children)-1
		else: Epoch2NonBurstNodeCount += len(node.children)-1

	Epoch1BurstBranchCount= len(Epoch1BurstBranch)
	Epoch1NonBurstBranchCount= len(Epoch1NonBurstBranch)

	if Epoch1BurstBranchCount != 0: BurstNodeToBranch = float(Epoch2BurstNodeCount) / Epoch1BurstBranchCount
	else: BurstNodeToBranch = 0

	if Epoch1NonBurstBranchCount != 0: NonBurstNodeToBranch = float(Epoch2NonBurstNodeCount) / Epoch1NonBurstBranchCount
	else: NonBurstNodeToBranch = 0

	if NonBurstNodeToBranch != 0: BurstToNonBurstRatio = BurstNodeToBranch / NonBurstNodeToBranch
	else: BurstToNonBurstRatio = 0

	epoch1BurstTips = set([])
	epoch1NonBurstTips = set([])

	tipPath = open(rootPath+'.DescendentsofBurstSize'+str(window)+'.csv','w')
	tipPath.write('mpv_uid,BurstDescendent\n')

	for node in ageDict:
		if node.is_leaf() == True:
			if node.name in burstDict['BurstTips']:
				BD = 'Y'
			else:
				BD = 'N'

			tipPath.write(str(node.name.split('|')[0])+','+str(BD)+'\n')

	tipPath.close()			

	return str(window)+','+str(Epoch1BurstBranchCount)+','+str(Epoch1NonBurstBranchCount)+','+str(Epoch2BurstNodeCount)+','+str(Epoch2NonBurstNodeCount)+','+str(round(BurstNodeToBranch,4))+','+str(round(NonBurstNodeToBranch,4))+','+str(round(BurstToNonBurstRatio,4))


"""
ORIGINAL PERIOD: 2014-2016 (inclusive), Follow on 2017-2019 (inclusive)
WHAT PROPORTION OF PEOPLE DIAGNOSED FROM 2014-2019 DESCEND OR ARE PART OF A BURST OCCURING POST-2014 BURST A 1-, 2-, OR 3-YEAR BURST?
START IN 2014 AND GO FORWARD IN TIME
	PROVIDE START DATE OF BURSTS
"""
def recentBurstTipsFixed(ageDict,epochStart,epochEnd,upperDxLimit,burstTime,events):
	dxFile = open(rootPath+'.DiagnosesFollowingAnyBurst.'+str(burstTime)+'Years.tsv','w')
	dxFile.write('mpv_uid,GenoDate,DxYear,BurstDuration,BurstStart,BurstEnd,TotalDescendents\n')

	for node in ageDict:
		nodeAge = ageDict[node]
		if node.is_leaf() == True and nodeAge >= epochStart:
			parent = node.up
			parentAgeList = [[parent, round(ageDict[parent],2)]]

			while ageDict[parent.up] >= epochStart:
				parent = parent.up
				parentAgeList.append([parent, round(ageDict[parent],2)])

			parentAgeList.sort(key=lambda x:x[1], reverse=False)

			if parentAgeList[0][1] >= epochStart:
				for i, refNode in enumerate(parentAgeList):
					burstCheck = 0
					refNodeName = refNode[0]
					refNodeAge = refNode[1]
					children = [[refNodeName,refNodeAge,len(refNodeName.children)-1]]

					for child in refNode[0].get_descendants():
						if child.is_leaf() == False and refNodeAge - ageDict[child] <= burstTime:
							children.append([child,ageDict[child],len(child.children)-1])		

					eventCount = 0

					children.sort(key=lambda x:x[1], reverse=False)

					for AB in children:
						if AB[1] - refNodeAge <= burstTime: eventCount += AB[2]
						if eventCount >= events:
							dxFile.write(node.name+','+str(nodeAge)+','+metaDict[node.name]['_hiv_aids_dx_dt_num'][0:4]+','+str(round(AB[1] - refNodeAge,2))+','+str(round(refNodeAge,2))+','+str(round(AB[1],2))+','+str(len(refNodeName.get_leaf_names()))+'\n')
							burstCheck = 1
							break

					if burstCheck == 1: break

	dxFile.close()


"""
RUNNING ACTUAL ANALYSES
"""

print ("\nREADING IN DATA FILE\n")

metaDict = makeMetaDict(dataFile,rootPath,epochStart)

if arg.skip_alignment == False:
	print ("BUILDING SEQUENCE ALIGNMENT\n")
	alignSeqs(rootPath,HXB2seq)

else: print ("SKIPPING SEQUENCE ALIGNMENT\n")



if arg.skip_tree == False: 
	print ("BUILDING PHYLOGENETIC TREE\n")
	makeTree(rootPath)

else: print ("SKIPPING PHYLOGENETIC TREE BUILDING\n")

if arg.skip_clock == False:
	print ("DATING PHYLOGENETIC TREE\n")
	clockTree(rootPath,epochStart)
else: print ("SKIPPING MOLECULAR CLOCK DATING")

print ("\nREADING IN TREE FILE\n")

treeData = readTree(rootPath)
t = treeData[0]
rootDate = treeData[1]
youngestTip = treeData[2]
rootName = t.get_tree_root()

print ("PRE-COMPUTING INTERNAL NODE AGES\n")
ageDict = makeAgeDict(t,rootName)

print ("ANALYZING TRANSMISSION BURSTS\n")
burstSummary = open(rootPath+'.BurstSummaries.csv','w')
burstSummary.write('WindowSize,BranchesInBurst,BranchesNotInBurst,BurstDecendentNodes,NonBurstDescendentNodes,RatioOfBurstNodesToBranches,RatioOfNonBurstNodesToBranches,RatioOfBurstToNonBurst\n')

for window in [0.5,1.0,1.5,2.0,2.5,3.0]: burstSummary.write(describeBurstRatio(window,ageDict,rootPath)+'\n')
burstSummary.close()

print ("FINDING DESCENDENTS OF RECENT BURSTS\n")
recentBurstTipsFixed(ageDict,epochStart,epochEnd,upperDxLimit,1,events)
recentBurstTipsFixed(ageDict,epochStart,epochEnd,upperDxLimit,2,events)
recentBurstTipsFixed(ageDict,epochStart,epochEnd,upperDxLimit,3,events)

print("FINISHED ANALYZING BURSTS\n")