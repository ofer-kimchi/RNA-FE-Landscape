This code refers to the paper "RNA structure prediction including pseudoknots through direct enumeration of states" by Kimchi et al. All files should be run in MatLab. It has been tested in version R2017b but is expected to work in previous versions as well.

The main file is RNALandscape_main.m. The inputs and outputs of that file are described in comments at the top of the file. The two .mat files saved are the results of running ennumerateDisallowedPermutationsWithPseudoknots. This function need only be run once (not once per sequence).

This code takes as input an RNA or DNA sequence, or multiple sequences together, and computes the entire free energy landscape -- all possible secondary structures and their probabilities. It also coarse grains over different topologies (if the boolean variable storeGraphs is set to true) and returns the possible topologies the sequence can fold into and the probability of each topology forming.

At this point, we only allow up to three strands of RNA/DNA interacting. We also only consider topologies which can be decomposed into the minimal subgraphs listed in the paper.

Sample code:

sequences = {'UGUGUCUUGGAUCGCGCGGGUCAAAUGUAUAUGGUUCAUAUACAUCCGCAGGCACGUAAUAAAGCGA'}; T = 300; vs = 0.02; duplexEntropyPenaltyInKB = 17; minBPInRegion = 4; numParForLoops = 0; pairwiseEnergies = true; storeGraphs = true; makingFigures = true; printProgressUpdate = true; allowParallelStrands = false; allowPseudoknots = true; minNumCompatible = 0; substems = 'all';

[sortedProbs, indexSortedProbs, sortedFE, indexSortedFE,graphProbs, weightMatrixList, numBondsMatrixList, startAndPermuTime, checkFxnTime, totalTime, STable, possiblePermutations, numRegions, expFEeff] = RNALandscape_main(sequences, T, vs, duplexEntropyPenaltyInKB, minBPInRegion, numParForLoops, pairwiseEnergies, storeGraphs, makingFigures, printProgressUpdate, allowParallelStrands, allowPseudoknots, minNumCompatible, substems);