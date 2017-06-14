import os
import random

###Program starts on line 588

def PopulateTempArr(TempArr, UpRegFile, DownRegFile, CombinedFile):
	###When this function is done, Temparr should have five elements
	###[study, dataset, upregdata, downregdata, combineddata]
	###each of the three last indices should be lists of lists of three 
	###elements that look like this [GeneSym, logFC, adj. P Val]
	UpRegDict = {}
	header = UpRegFile.readline().strip().split('\t')
	if len(header) > 1:
		PIndex = header.index("adj.P.Val")
		FCIndex = header.index("logFC")
		GeneSymIndex = header.index("ID")
		for line in UpRegFile:
			arr = line.strip().split('\t')

			WeirdGene = arr[GeneSymIndex].split(" /// ")
			if len(WeirdGene) > 1:
				for symbol in WeirdGene:
					if symbol in UpRegDict:
						UpRegDict[symbol].append(arr[FCIndex])
						UpRegDict[symbol].append(arr[PIndex])
					else:
						UpRegDict[symbol] = [arr[FCIndex], arr[PIndex]]
			else:
				if arr[GeneSymIndex] in UpRegDict:
					UpRegDict[arr[GeneSymIndex]].append(arr[FCIndex])
					UpRegDict[arr[GeneSymIndex]].append(arr[PIndex])
				else:
					UpRegDict[arr[GeneSymIndex]] =  [arr[FCIndex], arr[PIndex]]
	TempArr.append(UpRegDict)

	DownRegDict = {}
	header = DownRegFile.readline().strip().split('\t')
	if len(header) > 1:
		PIndex = header.index("adj.P.Val")
		FCIndex = header.index("logFC")
		GeneSymIndex = header.index("ID")
		for line in DownRegFile:
			arr = line.strip().split('\t')
			
			WeirdGene = arr[GeneSymIndex].split(" /// ")
			if len(WeirdGene) > 1:
				for symbol in WeirdGene:
					if symbol in DownRegDict:
						DownRegDict[symbol].append(arr[FCIndex])
						DownRegDict[symbol].append(arr[PIndex])
					else:
						DownRegDict[symbol] = [arr[FCIndex], arr[PIndex]]
			else:
				if arr[GeneSymIndex] in DownRegDict:
					DownRegDict[arr[GeneSymIndex]].append(arr[FCIndex])
					DownRegDict[arr[GeneSymIndex]].append(arr[PIndex])
				else:
					DownRegDict[arr[GeneSymIndex]] =  [arr[FCIndex], arr[PIndex]]
	TempArr.append(DownRegDict)

	CombinedDict = {}
	header = CombinedFile.readline().strip().split('\t')
	if len(header) > 1:
		PIndex = header.index("adj.P.Val")
		FCIndex = header.index("logFC")
		GeneSymIndex = header.index("ID")
		for line in CombinedFile:
			arr = line.strip().split('\t')
			
			WeirdGene = arr[GeneSymIndex].split(" /// ")
			if len(WeirdGene) > 1:
				for symbol in WeirdGene:
					if symbol in CombinedDict:
						CombinedDict[symbol].append(arr[FCIndex])
						CombinedDict[symbol].append(arr[PIndex])
					else:
						CombinedDict[symbol] = [arr[FCIndex], arr[PIndex]]
			else:
				if arr[GeneSymIndex] in CombinedDict:
					CombinedDict[arr[GeneSymIndex]].append(arr[FCIndex])
					CombinedDict[arr[GeneSymIndex]].append(arr[PIndex])
				else:
					CombinedDict[arr[GeneSymIndex]] =  [arr[FCIndex], arr[PIndex]]

	TempArr.append(CombinedDict)

def BuildRosettaStone():
	#constructs the Rosetta Stone
	RosettaStone = []
	RSFile = open("RosettaStone.txt", 'r')
	for line in RSFile:
		arr = line.split('\t')
		if len(arr) > 1:
			RosettaStone.append(line.strip().split("\t"))
	
	###Generate dictionaries that connect each human gene to species orthologues
	###and dictionaries that connect each species orthologue to its human gene
	PSDCounts = []
	PairwiseSpeciesDicts = []
	for x in range(len(RosettaStone[0])):
		PairwiseSpeciesDicts.append([RosettaStone[0][0]+"_"+RosettaStone[0][x], {}])
		PSDCounts.append(0)

	####Example table...
	#### Human ||  Dog || Zebrafish
	#### ABC   || ABCD || abc, abcl
	####link ABC to ABCD

	SpeciesIndices = RosettaStone[0]
	###dictionary of all human genes in the RosettaStone
	for x in range(1, len(RosettaStone)):
		for y in range(0, len(RosettaStone[x])):
			SpeciesGenes = list(set(RosettaStone[x][y].strip('"').split(",")))
			###for each y index, map the human gene to the first gene in the list
			if RosettaStone[x][y] == "NO_ORTHO":
				continue
			SpeciesToHumanIndex = FindPSDIndex(SpeciesIndices[0], SpeciesIndices[y], PairwiseSpeciesDicts)
				###this should set the species gene = to the set of human orthologues
				###SpeciesToHumanIndex is the correct species comparison to add the gene to
				###RosettaStone[x][0] is the human gene
			for Gene in SpeciesGenes:
				if Gene not in PairwiseSpeciesDicts[SpeciesToHumanIndex][1]:
					PairwiseSpeciesDicts[SpeciesToHumanIndex][1][Gene] = set()
					PairwiseSpeciesDicts[SpeciesToHumanIndex][1][Gene].add(RosettaStone[x][0])
				else:
					PairwiseSpeciesDicts[SpeciesToHumanIndex][1][Gene].add(RosettaStone[x][0])

	return PairwiseSpeciesDicts, SpeciesIndices

def CalcMaximumSharedGenes(PSD):
	###the purpose of this function is to calculate the maximum number of
	###shared genes between datasets for Monte Carlo Analysis
	###this involves 1) counting the number of Human genes in the Rosetta 
	###stone that are mapped by the dataset and 2) counting the number of 
	###potential shared genes in the Rosetta stone based on the mapped
	###human genes

	StudiesAndMappedGenes = []
	for study in os.listdir(os.getcwd()):
		if  "." in study or study == "Stratifications" or study == "TestData" or study == "OtherStuff":
			continue
		if study == "OLD_SHARMA":
			continue
		Species = study.split('_')[1]

		###first see how many genes in each "CLEAN_FINAL_EXPRESSION_SET.txt"
		###file map to the human_species comparison in PSD
		#SpeciesPSDIndex AKA SPI is the index where the dictionary of human_species is found
		SPI= FindPSDIndex("human", Species, PSD)

		###read in all genes from the study's expression file 
		###and try to map them to PSD[Species1PSDIndex]
		StudyExpressionFile = open("./"+study+"/CLEAN_FINAL_EXPRESSION_SET.txt", 'r')
		###MissedGenesFile contains all the genes that didn't map using the Rosetta Stone
		MissedGenesFile = open("./"+study+"/MissedGenes.txt", 'w')
		MissedGenes = {}
		###this will go through and find all occurences where
		###a species gene is found in the Rosetta Stone and 
		###write these into a new dictionary called found_dict
		found_dict = {}
		for line in StudyExpressionFile:
			arr = line.strip().split('\t')
			gene = arr[-1].strip('"')
			if "///" in gene:
				MultiGeneArr = gene.split(" /// ")
				for x in MultiGeneArr:
					if x in PSD[SPI][1]:
						if x not in found_dict:
							found_dict[x] = set(PSD[SPI][1][x])
						else:
							for y in PSD[SPI][1][x]:
								found_dict[x].add(y)
					else:
						MissedGenes[x] = None
			elif gene in PSD[SPI][1]:
				if gene not in found_dict:
					found_dict[gene] = set(PSD[SPI][1][gene])
				else:
					for y in PSD[SPI][1][gene]:
						found_dict[gene].add(y)
			else:
				MissedGenes[gene] = None
		MissedGenesFile.write("Missed Genes: " + str(len(MissedGenes))+'\n')
		for x in MissedGenes:
			MissedGenesFile.write(x+'\n')
		MissedGenesFile.close()

		###found_dict now has all the species genes from the study that mapped with human genes 
		###in the Rosetta stone
		StudiesAndMappedGenes.append([study, found_dict])

	###now do comparisons between each study in StudiesAndMappedGenes to see
	###how many Rosetta Stone genes they potentially share
	for Study in StudiesAndMappedGenes:
		TempSet = set()
		for gene in Study[1]:
			TempList = list(Study[1][gene])
			for HumanGeneSym in TempList:
				TempSet.add(HumanGeneSym)
		Study.append(TempSet)
	###StudiesAndMappedGenes now contains 3 elements
	###element 0 is the study name
	###element 1 is a dictionary of species genes = set of human genes
	###element 2 is a set of all unique human genes that were found in study

	PairwiseMaxGenes = []
	for x in StudiesAndMappedGenes:
		for y in StudiesAndMappedGenes:
			TempArr = [x[0],y[0]]
			CommonGenes = x[2].intersection(y[2])
			TempArr.append(CommonGenes)
			PairwiseMaxGenes.append([x[0],y[0],len(CommonGenes), CommonGenes])

	return PairwiseMaxGenes

def FindPSDIndex(species1, species2, PSD):
	###given two species and the list of PairwiseSpeciesDictionarys (PSD)
	###this function will find the correct PSD index for the comparison 
	###between the two species
	species1 = species1.lower()
	species2 = species2.lower()

	for x in range(len(PSD)):
		if PSD[x][0].lower() == species1 + "_" + species2:
			return x
		elif PSD[x][0].lower() == species2 + "_" + species1:
			print("UNEXPECTED RESULT IN FindPSDIndex")
			print(species1, species2)
			return x
	return -1

def BuildRobsTable(PMG, StratArr, SpeciesIndices, strat, PSD):

	UpFile = open("./Stratifications/" +strat+ "/"+strat+ "_UpGenes.txt", 'w')
	DownFile = open("./Stratifications/" +strat+ "/"+strat+ "_DownGenes.txt", 'w')
	CombinedFile = open("./Stratifications/" +strat+ "/"+strat+ "_CombinedGenes.txt", 'w')

	###All the human genes that need to be in the Rosetta Stone
	HumanGenes = PSD[0][1]

	###put headers on all the files
	header = "HumanGene"
	for Dataset in StratArr:
		header += "\t" + Dataset[0].split("_")[0] + "_" + Dataset[1]
	header += '\tDataSetOccurences\n'
	UpFile.write(header)
	DownFile.write(header)
	CombinedFile.write(header)

	HumanFirstDatasets = []
	for dataset in StratArr:
		###dataset[0] is the study name Author_Species
		###dataset[1] is the dataset
		###dataset[2] upregulated data
		###dataset[3] downregulated data
		###dataset[4] combined data    
		###all data look like dict[SpeciesGene] = [FC, p-val]
		###first put all up data in terms of the human gene.  
		###it should go from this... dict[SpeciesGene] = [FC, p-val]
		###to this... dict[HumanGene] = ["SpeciesGene_FC_p-val"]
		species = dataset[0].split("_")[1]
		index = FindPSDIndex("human", species, PSD)
		TempArr = [dataset[0],dataset[1]]

		UpDict = {}
		for gene in dataset[2]:
			#print(gene)
			if gene in PSD[index][1]:
				HGList = list(PSD[index][1][gene])
				for HG in HGList:
					if HG in UpDict:
						UpDict[HG].append(gene+"_"+dataset[2][gene][0]+"_"+dataset[2][gene][1])
					else:
						UpDict[HG] = [gene+"_"+dataset[2][gene][0]+"_"+dataset[2][gene][1]]
		TempArr.append(UpDict)

		DownDict = {}
		for gene in dataset[3]:
			if gene in PSD[index][1]:
				HGList = list(PSD[index][1][gene])
				for HG in HGList:
					if HG in DownDict:
						DownDict[HG].append(gene+"_"+dataset[3][gene][0]+"_"+dataset[3][gene][1])
					else:
						DownDict[HG] = [gene+"_"+dataset[3][gene][0]+"_"+dataset[3][gene][1]]
		TempArr.append(DownDict)

		CombinedDict = {}
		for gene in dataset[4]:
			if gene in PSD[index][1]:
				HGList = list(PSD[index][1][gene])
				for HG in HGList:
					if HG in CombinedDict:
						CombinedDict[HG].append(gene+"_"+dataset[4][gene][0]+"_"+dataset[4][gene][1])
					else:
						CombinedDict[HG] = [gene+"_"+dataset[4][gene][0]+"_"+dataset[4][gene][1]]
		TempArr.append(CombinedDict)

		HumanFirstDatasets.append(TempArr)

	for HumanGene in HumanGenes:
		if HumanGene == "":
			continue
		UpLine = HumanGene + '\t'
		DownLine = HumanGene + '\t'
		CombinedLine = HumanGene + '\t'
		###count the number of dataset occurences
		UpCount = 0
		DownCount = 0
		CombinedCount = 0

		for Dataset in HumanFirstDatasets:
			if HumanGene in Dataset[2]:
				UpCount += 1
				for x in Dataset[2][HumanGene]:
					UpLine += x + ","
				UpLine.strip(',')
			UpLine += '\t'
			if HumanGene in Dataset[3]:
				DownCount += 1
				for x in Dataset[3][HumanGene]:
					DownLine += x + ','
				DownLine.strip(',')
			DownLine += '\t'
			if HumanGene in Dataset[4]:
				CombinedCount += 1
				for x in Dataset[4][HumanGene]:
					CombinedLine += x + ','
				CombinedLine.strip(",")
			CombinedLine += '\t'
		UpFile.write(UpLine+str(UpCount)+'\n')
		DownFile.write(DownLine+str(DownCount)+'\n')
		CombinedFile.write(CombinedLine+str(CombinedCount)+'\n')


	UpFile.close()
	DownFile.close()
	CombinedFile.close()
	SuccessFile = open("./Stratifications/" +strat+ "/RTMade.txt", 'w')
	SuccessFile.close()

def SetUpMonteCarlo(ExpressionFile, PMG, strat, genegroup, PSD):
	###ExpressionFile is the RT file
	header = ExpressionFile.readline().strip().split('\t')
	datasets = header[1:-1]

	#identify the directory that contains the expresion data for each dataset
	#dir ends up being a string that looks like "Study_species"
	DatasetToRSMaps = []
	for x in datasets:
		DatasetRoot = x.split("_")[0]
		for dir in os.listdir(os.getcwd()):
			if DatasetRoot in dir:
				DatasetToRSMaps.append([dir])

	#this section sets up the the ability to randomly sample genes from each dataset 
	#by mapping genes found in each dataset to DatasetToRSMap via the Rosetta Stone
	for x in DatasetToRSMaps:
		DatasetToRSMap = {}
		DatasetFile = open("./"+x[0]+"/CLEAN_FINAL_EXPRESSION_SET.txt",'r')
		species = x[0].split("_")[1]
		PSDIndex = FindPSDIndex("human", species, PSD)

		#for each line in the expression file for the dataset, if the gene maps to the PSD
		#Add the mapped gene to DatasetToRSMap, 
		for line in DatasetFile:
			Gene = line.strip().split('\t')[-1].strip('"')
			if Gene in PSD[PSDIndex][1]:
				#print(PSD[PSDIndex][1][Gene])
				DatasetToRSMap[Gene] = PSD[PSDIndex][1][Gene]
			else:
				DatasetToRSMap[Gene] = None
		x.append(DatasetToRSMap)
	###we can now randomly pick genes from the DatasetToRSMaps

	###now figure out how many DE genes are in each dataset and how many DE genes are shared between datasets
	###HumanDEGeneLists will keep track of which RS genes were DE in each dataset
	HumanDEGeneLists = []
	###SpeciesDEGeneLists will keep a set of all species genes that were DE in the dataset
	###the number of unique DE genes mapped to RS will be used as the number of genes to pull 
	###randomly from DatasetToRSMaps
	SpeciesDEGeneLists = []
	for d in datasets:
		HumanDEGeneLists.append([d, set()])
		SpeciesDEGeneLists.append([d, set()])

	###we need to go to the limmaout files to create SpeciesDEGeneLists
	###this seems to work
	###The limmaout file contains all DE genes for a given dataset.  
	for x in range(len(DatasetToRSMaps)):
		filenamearr = datasets[x].split('_')[1:]
		filename = ""
		for zzz in filenamearr:
			filename += zzz + "_"
		filename = filename.strip('_')
		if genegroup == "Down":
			limmafile = open("./"+DatasetToRSMaps[x][0]+"/Limmaout/"+filename+"_down.txt", 'r')
			limmafile.readline()
		elif genegroup == "Up":
			limmafile = open("./"+DatasetToRSMaps[x][0]+"/Limmaout/"+filename+"_up.txt", 'r')
			limmafile.readline()
		for line in limmafile:
			arr = line.split('\t')
			if len(arr) > 1:
				SpeciesDEGeneLists[x][1].add(arr[1])
		limmafile.close()

	count = 0
	for line in ExpressionFile:
		#Count tracks our progress through the RT file and should add up to the number of Human Genes 
		count += 1
		#Split each line of the file into an array
		arr = line.split('\t')[1:-1]
		for ind in range(len(arr)):
			#if an entry != "" then there is a DE gene for that dataset
			#therefore, add the count index to that datasets list of DE genes
			if arr[ind] != "":
				HumanDEGeneLists[ind][1].add(count)
	ExpressionFile.close()

	#now we can pick len(SpeciesDEGeneLists[ind][1]) genes from DatasetToRSMaps[ind]
	print("Pairwise")
	#this will perform pairwise monte carlo simulations to assess if the number of shared DE genes between datasets is significant
	PairwisePvalDict = PairwiseMonteCarlo(strat, datasets, DatasetToRSMaps, HumanDEGeneLists, SpeciesDEGeneLists, genegroup)
	print("Distribution")
	#we run this to test if the distribution of shared DE genes between all datasets is significant.  We didn't include this in the manuscript and it is not guaranteed to work
	DistributionPValArr, ActualDistribution, AverageRandomDistribution = DistributionMonteCarlo(strat, datasets, DatasetToRSMaps, HumanDEGeneLists, SpeciesDEGeneLists, genegroup)

	return PairwisePvalDict, DistributionPValArr, ActualDistribution, AverageRandomDistribution, datasets

def PairwiseMonteCarlo(strat, datasets, DatasetToRSMaps, HumanDEGenesLists, SpeciesDEGeneLists, genegroup):
	#store the results of the simulations in SimulationResultsDict
	SimulationResultsDict = {}

	#this should track the number of simulations that need to be run
	compcount = 0
	for x in range(len(datasets)):
		compcount += x

	#the PairwiseDistributionFile is a file that tracks the distribution of outcomes from the Monte Carlo simulation
	#we used this for debugging purposes and to create examples of how the analysis works
	PairwiseDistributionFile = open("./Stratifications/"+strat+"/"+strat+genegroup+"_PairwiseDistributions.txt", 'w')

	DatasetCompCount = 0
	for x in range(len(datasets)):
		for y in range(x+1, len(datasets)):
			DatasetCompCount += 1
			comparison = datasets[x] + "_" + datasets[y]
			#find how many RS genes are shared between the two datasets
			NumberOfOverlappingRSGenes = len(HumanDEGenesLists[x][1].intersection(HumanDEGenesLists[y][1]))
			if NumberOfOverlappingRSGenes == 0:
				SimulationResultsDict[datasets[x],datasets[y]] = 1
				continue
			#find how many species genes are in each SpeciesDEGeneList
			#this is the number of genes we will randomly select in the MC simulation.  
			XSpeciesGenes = len(SpeciesDEGeneLists[x][1])
			YSpeciesGenes = len(SpeciesDEGeneLists[y][1])
			
			###run the MC simulations
			#GreaterThanCount tracks the number of simulations that produce a number of shared genes than is empirically observed
			GreaterThanCount = 0
			IterCount = 10000
			DistributionTracker = {}
			for iter in range(IterCount):
				#this carriage return line seems to work when the script is run on Cygwin but not from the python shell
				print "Comparison ", DatasetCompCount, " of ", compcount, " MC Iterations Completed ", iter, "                                 \r",
				#randomly select genes from DatasetToRSMaps
				XGenes = random.sample(DatasetToRSMaps[x][1],XSpeciesGenes)
				YGenes = random.sample(DatasetToRSMaps[y][1],YSpeciesGenes)
				XHumanGenes = set()
				YHumanGenes = set()
				#map the selected species genes to the Rosetta Stone
				for gene in XGenes:
					if DatasetToRSMaps[x][1][gene] != None:
						XHumanGenes.update(DatasetToRSMaps[x][1][gene])
				for gene in YGenes:
					if DatasetToRSMaps[y][1][gene] != None:
						YHumanGenes.update(DatasetToRSMaps[y][1][gene])
				#find the randomly selected genes that overlap
				SharedGenes = XHumanGenes.intersection(YHumanGenes)
				NumSharedGenes = len(SharedGenes)
				
				#add the number of randomly shared genes to the distribution tracker
				if NumSharedGenes in DistributionTracker:
					DistributionTracker[NumSharedGenes] += 1
				else:
					DistributionTracker[NumSharedGenes] = 1

				if NumSharedGenes >= NumberOfOverlappingRSGenes:
					GreaterThanCount += 1

			PairwiseDistributionFile.write(datasets[x] + " vs " + datasets[y] +'\n')
			RandomSharedGenesLine = "Number of Overlapping Genes in MC Simulation:"
			SimulationsNumLine = "Number of MC Simulations:"
			#The program is nice enough to access the elements of DistributionTracker in numerical order
			#I didn't realize that dictionaries worked that way...
			for z in DistributionTracker:
				RandomSharedGenesLine += '\t'+ str(z)
				SimulationsNumLine += '\t' + str(DistributionTracker[z])
			PairwiseDistributionFile.write(RandomSharedGenesLine +'\n')
			PairwiseDistributionFile.write(SimulationsNumLine+'\n')
			PairwiseDistributionFile.write("Actual Number Overlapping Genes:\t"+str(NumberOfOverlappingRSGenes)+'\n')

			#The pval is calculated as the proportion of simulations that give a number of randomly 
			#shared genes that is greater than or equal to the number of empirically shared genes
			pval = float(GreaterThanCount)/float(IterCount)
			SimulationResultsDict[datasets[x],datasets[y]] = pval

			PairwiseDistributionFile.write("P-Value:\t"+str(pval)+'\n\n')


	PairwiseDistributionFile.close()
	#for x in SimulationResultsDict:
	#   print(x, SimulationResultsDict[x])
	return SimulationResultsDict

def DistributionMonteCarlo(strat, datasets, DatasetToRSMaps, HumanDEGeneLists, SpeciesDEGeneLists, genegroup):
	AllGenesWithOccurrences = {}
	for x in HumanDEGeneLists:
		for gene in x[1]:
			if gene in AllGenesWithOccurrences:
				AllGenesWithOccurrences[gene] += 1
			else:
				AllGenesWithOccurrences[gene] = 1
	ActualDistribution = [0]
	RandomDistCopy = [0]
	GreaterThanDist = [0]
	for x in datasets:
		ActualDistribution.append(0)
		RandomDistCopy.append(0)
		GreaterThanDist.append(0)
	for gene in AllGenesWithOccurrences:
		ActualDistribution[AllGenesWithOccurrences[gene]] += 1
	actualsums = list(RandomDistCopy)
	for x in range(len(ActualDistribution)):
		for y in range(x,len(ActualDistribution)):
			actualsums[x] += ActualDistribution[y]
	ARD = list(RandomDistCopy)

	itercount = 10000
	for iter in range(itercount):
		print " MC Iterations Completed ", iter, "              \r",
		RandomDist = list(RandomDistCopy)
		AllRandomGenes = {}
		for x in range(len(datasets)): 
			DatasetGenes = random.sample(DatasetToRSMaps[x][1], len(SpeciesDEGeneLists[x][1]))
			TempDatasetGeneList = {}
			for gene in DatasetGenes:
				if DatasetToRSMaps[x][1][gene] != None:
					for Humgene in DatasetToRSMaps[x][1][gene]:
						TempDatasetGeneList[Humgene] = None 
			for x in TempDatasetGeneList:
				if x in AllRandomGenes:
					AllRandomGenes[x] += 1
				else:
					AllRandomGenes[x] = 1   
		
		for gene in AllRandomGenes:
			RandomDist[AllRandomGenes[gene]] += 1

		###we need to count how many times an index plus all indexes that are greater than the index have a number of shared genes that is greater than the observed value  
		randomsums = list(RandomDistCopy)
		for x in range(len(RandomDist)):
			ARD[x]+= RandomDist[x]
			for y in range(x,len(RandomDist)):
				randomsums[x] += RandomDist[y]
		
		for x in range(len(randomsums)):
			if randomsums[x] > actualsums[x]:
				GreaterThanDist[x] += 1

	PValArr = list(RandomDistCopy)
	for x in range(len(GreaterThanDist)):
		PValArr[x] = float(GreaterThanDist[x])/float(itercount)
		ARD[x] = float(ARD[x])/float(itercount)

	return PValArr, ActualDistribution, ARD

"""The purpose of this script is to take normalized expression data from multiple studies and compare the Differentially Expressed (DE) genes
found in each.  To compare between species we use a comparison table called the Rosetta Stone.  Statistical significance of the number of
overlapping DE genes between datasets is assessed using monte carlo simulations. For more info, or to express comments/concerns feel free
to contact me at rjdonahue@wisc.edu"""

###There are multiple lines of code that may seem extraneous and these are generally used for debugging or for tracking the progress of the program

#PSD = Pairwise Species Dicts is a list of dictionaries that link human genes to genes from each species 
#and genes from each species to their human genes.  
PSD, SpeciesIndices = BuildRosettaStone()

#Calculate the maximum number of shared genes between each species
#an return Pairwise Maximum shared Genes (PMG) which is an array 
#that details the maximum number of genes two datasets may share
PMG = CalcMaximumSharedGenes(PSD)

###for each stratification of the data we want to construct a Rob's Table (RT) 
###and calculate the number of overlapping genes between datasets
for strat in os.listdir("./Stratifications/"):
	###StratArr will hold all the pertinent information about a dataset
	###including title, species and DE genes
	StratArr = []
	if ".txt" in strat or".xlsx" in strat or ".py" in strat or "OLD" in strat or "pptm" in strat or "docx" in strat or "FIGURES" in strat or "Done" in strat:
		continue
	#if an RT has already been made, don't make it again
	if "RTMade.txt" not in os.listdir("./Stratifications/"+strat):
		DatasetsFile = open("./Stratifications/"+strat+"/Datasets.txt",'r')
		###each line should be "author_species dataset"
		for line in DatasetsFile:
			if line == "\r\n":
				continue
			arr = line.split(" ")
			TempArr = []
			study = arr[0]
			dataset = arr[1].strip()
			TempArr.append(study)
			TempArr.append(dataset)

			###now go find the correct Limmaout folder and get the information
			###the downregulated genes
			DownRegFile = open(study+'/Limmaout/'+dataset+"_down.txt", 'r')
			###the upregulated genes
			UpRegFile = open(study+'/Limmaout/'+dataset+"_up.txt", 'r')
			###the combined down and up regulated genes
			CombinedFile = open(study+'/Limmaout/'+dataset+"_updown.txt", 'r')

			PopulateTempArr(TempArr, UpRegFile, DownRegFile, CombinedFile)
			StratArr.append(TempArr)
			DownRegFile.close()
			UpRegFile.close()
			CombinedFile.close()

	else:
		#print(strat, " already done")
		continue

	###we now have a built RosettaStone named PSD
	###and the expression data for all pertinent datasets
	###in  StratArr... we need to organize all the data
	###in Strat Arr into Rob's Table. 
	BuildRobsTable(PMG, StratArr, SpeciesIndices, strat, PSD)

###now count the # of pairwise shared genes and do MonteCarlo Analysis on those
###there are functions to perform PairwiseMonteCarlo and Distribution Monte Carlo
###note that we never used the Distribution Monte Carlo analysis in our manuscript
for strat in os.listdir("./Stratifications"):
	if ".txt" in strat or".xlsx" in strat or ".py" in strat or "OLD" in strat or "pptm" in strat or "docx" in strat or "FIGURES" in strat or "Done" in strat:
		continue

	UpRTFile = open("./Stratifications/"+strat+"/"+strat+"_UpGenes.txt", 'r')
	DownRTFile = open("./Stratifications/"+strat+"/"+strat+"_DownGenes.txt", 'r')

	###upregulated stuff 
	print(strat)
	print("UpReg")
	UpPairwisePValDict, UpDistributionPValArr, UpDistribution, UpARD, UpDatasets = SetUpMonteCarlo(UpRTFile, PMG, strat, "Up", PSD)
	###downregulated stuff
	print("DownReg")
	DownPairwisePValDict, DownDistributionPValArr, DownDistribution, DownARD, DownDatasets = SetUpMonteCarlo(DownRTFile, PMG, strat, "Down", PSD)

	UpRTFile.close()
	DownRTFile.close()

	###write data into a table
	outfile = open("./Stratifications/"+strat+"/"+strat+"_MonteCarloResults.txt", 'w')
	outfile.write("Results of Pairwise Monte Carlo Analysis of Upregulated Genes\n\n")
	
	StratDatasets = UpDatasets
	header = "Datasets"
	for x in StratDatasets: 
		header += '\t'+ x 
	header += '\n'
	outfile.write(header)

	count = 0
	for x in range(len(StratDatasets)):
		line = StratDatasets[x] + '\t'
		for number in range(count):
			line += '\t' 
		line += "NA\t"
		for y in range(x+1, len(StratDatasets)):
			if (StratDatasets[x], StratDatasets[y]) in UpPairwisePValDict:
				line += str(UpPairwisePValDict[StratDatasets[x], StratDatasets[y]]) + '\t'
			elif (StratDatasets[y], StratDatasets[x]) in UpPairwisePValDict:
				line += str(UpPairwisePValDict[StratDatasets[y], StratDatasets[x]]) + '\t'
			else:
				print("KEY ERROR UP")
				print(StratDatasets[y], StratDatasets[x])
		outfile.write(line.strip('\t') + '\n')
		count += 1

	outfile.write("Results of Distribution Monte Carlo Analysis of UpRegulated Genes\n")
	HeaderLine = "Number of Datasets In Which A Gene Occurs"
	RealDistLine = "Actual Gene Distribution"
	ARDLine = "Average of the Random Distributions"
	PValline = "P-Value From 10000 Simulations"
	for x in range(len(UpDistribution)):
		HeaderLine += '\t'+str(x) 
		RealDistLine += '\t' + str(UpDistribution[x])
		ARDLine += '\t' + str(UpARD[x])
		PValline += '\t' + str(UpDistributionPValArr[x])
	outfile.write(HeaderLine+'\n')
	outfile.write(RealDistLine+'\n')
	outfile.write(ARDLine+'\n')
	outfile.write(PValline+'\n')

	outfile.write("\n\nResults of Pairwise Monte Carlo Analysis of Downregulated Genes\n\n")

	#DownPairwisePValDict, DownDistributionPValArr, DownDistribution, DownARD, DownDatasets

	StratDatasets = DownDatasets
	header = "Datasets"
	for x in StratDatasets: 
		header += '\t'+ x 
	header += '\n'
	outfile.write(header)

	count = 0
	for x in range(len(StratDatasets)):
		line = StratDatasets[x] + '\t'
		for number in range(count):
			line += '\t' 
		line += "NA\t"
		for y in range(x+1,len(StratDatasets)):
			if (StratDatasets[x], StratDatasets[y]) in DownPairwisePValDict:
				line += str(DownPairwisePValDict[StratDatasets[x], StratDatasets[y]]) + '\t'
			elif (StratDatasets[y], StratDatasets[x]) in DownPairwisePValDict:
				line += str(DownPairwisePValDict[StratDatasets[y], StratDatasets[x]]) + '\t'
			else:
				print("KEY ERROR DOWN")
				print(StratDatasets[y], StratDatasets[x])
		outfile.write(line.strip('\t') + '\n')
		count += 1

	outfile.write("Results of Distribution Monte Carlo Analysis of DownRegulated Genes\n")
	HeaderLine = "Number of Datasets In Which A Gene Occurs"
	RealDistLine = "Actual Gene Distribution"
	ARDLine = "Average of the Random Distributions"
	PValline = "P-Value From 10000 Simulations"
	for x in range(len(DownDistribution)):
		HeaderLine += '\t'+str(x) 
		ARDLine += '\t' + str(DownARD[x])
		RealDistLine += '\t' + str(DownDistribution[x])
		PValline += '\t' + str(DownDistributionPValArr[x])
	outfile.write(HeaderLine+'\n')
	outfile.write(RealDistLine+'\n')
	outfile.write(ARDLine +'\n')
	outfile.write(PValline+'\n')

	outfile.close()
