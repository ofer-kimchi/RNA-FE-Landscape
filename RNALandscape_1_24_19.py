#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 14:37:24 2018

@author: Ofer
"""
# =============================================================================
# Things left to do (not urgent so I haven't done them yet):
# make structure/graph histograms, parallelizing, specialLoops, different vs for the two strands,
# have inputs corruptE_Std, corruptS_Std, include corruptFE_Std
# =============================================================================

import numpy as np
import copy
import time
import scipy
import pickle #for saving and loading
import glob, os #for saving and loading
import cProfile #for profiling
import itertools
import networkx as nx #for graphs
import networkx.algorithms.isomorphism as iso #for graph isomorphism comparison
import matplotlib.pyplot as plt #for plotting

# =============================================================================
# # ===========================================================================
# # # =========================================================================
# # # 
# # #                                   TO RUN THE CODE:
# # #
# # # a = RNALandscape(['GCGCAAAAGCGC'],storeGraphs=True) #define the object with the parameters you want
# # # sol = a.calculateFELandscape() #calculate the FE Lanscape
# # # # The results are stored in, e.g., sol.sortedFEs (the free energies, sorted from lowest to highest)
# # #
# # # =========================================================================
# # ===========================================================================
# =============================================================================

class RNALandscape:
    
    def __init__(self, sequences, DNA = [False,False], sequenceFrequencies = [1,1], b = 0.8/0.33, T = 310.15, 
                  vs=0.02, duplexEntropyPenaltyInKB = 0, minBPInStem = 2, numParForLoops = 1, 
                  storeGraphs = False, makeFigures = True, printProgressUpdate = True, 
                  allowPseudoknots = True, minNumStemsInStructure = 0, substems = 'all', frozenBPs = [],
                  toSave = True, tryToLoadVariables = True, maxSizeSTable = 10**4, 
                  onlyAllowSubsetsOfLongestStems = False, onlyConsiderSubstemsFromEdges = False,
                  onlyConsiderBondedStrands = False, considerTabulatedLoops = False, 
                  unmatchedBPPenalty = True, minNtsInHairpin = 3, considerAllAsTerminalMismatches = False,
                  includeTerminalMismatches = True, includeFlushCoaxialStacks = True, 
                  includeDanglingEnds = True, includeTerminalAUATPenalties = True, 
                  considerC3andC4 = True,  corruptFESeed = False, unboundButCouldBindPenalties = [0,0]):
        
        
#From the MatLab, we removed the following functionalities: 
# pairwiseEnergies is now forced to be True;
# allowParallelStrands is now forced to be False
# =============================================================================
# Description of inputs:
#         
# This function takes as input an array of sequences (eg:
#     ['AUGCGC','GCGCAU']) and outputs the entire free energy landscape of these sequences interacting.
#     These sequences can be DNA or RNA. The code can be run with only one sequence e.g.
#     ['GCGAAAACGC']. 
 
# DNA is a list the same length as sequences. For each sequence, it is True if 
#     that sequence is DNA and False if it is RNA

# sequenceFrequencies has no effect (it is a placeholder; in the future we should be able to give
#        More precise entropic costs of binding based on the concentrations of each sequence)

# b is the persistence length of single stranded RNA, measured in units of nucleotides (nts) (~3.3 angstroms) 
#     This assumes a persistence length of RNA of ~0.8 nm found in The Snakelike Chain Character of Unstructured RNA
#     by  David R. Jacobson, Dustin B. McIntosh, and Omar A. Saleh
#     As a plus, this agrees well with Siggia's formula b=2.5a (where a is the size of one ntd).
#     This isn't very different for s.s. DNA. See Ionic strength-dependent persistence lengths of single-stranded RNA and DNA
#     by Huimin Chen, Steve P. Meisburger, Suzette A. Pabit, Julie L. Sutton, Watt W. Webb, and Lois Pollack
#
#     We use g = 3/(2*b) = gamma in place of b for practical reasons 
#     (that's the way b enters into the entropy formulae).
        
# T sets the temperature, which affects the free energy calculation FE=E-TS.
#     It is measured in Kelvin

# vs is the volume (in ntds^3) within which two nucleotides (ntds or nts, I use 
#     them interchangably) can bond.

# duplexEntropyPenaltyInKB only affects runs with more than one sequence --
#     it defines the entropy penalty (in units of Boltzmann's constant) of two
#     strands being bound. This will in general depend on the size of the container and the
#     relative concentrations of the two sequences. Future iterations of this
#     code will take the latter into account as well.

# minBPInStem defines the minimum length of a stem (it is equal to the
#     parameter m in the paper). It is an integer > 0. min value it can take is 
#     1; Pipas & McMahon set it to 3. 
# 
# numParForLoops affects the run itself. Set it to 0 to use a for loop
#     rather than a parfor loop (parallel computing) to calculate free energies.
#     Otherwise set it to any integer >0 to define how many iterations of parfor
#     will be run.
# 
# pairwiseEnergies is a boolean. It is true if you want to use the Turner
#     energy parameters. Otherwise, we have defined a simple model that uses
#     only the energies of single base pair (bp) interactions based on work by
#     Tristan Cragnolini. This functionality was removed in the Python code 
#     (it is now forced to be True)
# 
# storeGraphs is a boolean. Set it to True if you'd like to store the
#     topology of each structure and all the topologies. The code runs faster if
#     it is set to False, but this disallows coarse-graining by topology.
# 
# makeFigures defines whether you'd like the code to output figures at the
#     end of the calculation. It is a boolean
# 
# printProgressUpdate is a boolean. Set it to true to print messages
#     updating you on how the algorithm is progressing (e.g. how long it's spent
#     on each subfunction).
# 
# allowParallelStrands defines whether to allow parallel stems in the
#     computation. If pairwiseEnergies is True, it should be set to False;
#     otherwise, you can set it to True if you'd like. This functionality was removed in the Python code 
#     (it is now forced to be False)
        
# allowPseudoknots defines whether to allow pseudoknots in the computation.
#     It is a boolean. It does not affect ``pseudoknots" created by the binding of two strands.

# minNumStemsInStructure (integer) defines the minimum number of compatible stems that can
#     comprise a structure. Set it to 0 to get the most general results.
#     Otherwise set it to an integer > 0. If it is zero, it also includes the fully unfolded chain
#     If it is 1 it includes structures created by single stems; etc.
#     In the MatLab code, this parameter was called minNumCompatible.

# substems defines whether to allow all possible substems. If a stem of
#     length l is found to be possible, but minBPInStem = m is < l, subsets of this
#     stem of any length between l & m are also allowed. If substems is set to 'all', all of
#     these possible stems will be considered as well. If substems is set to an
#     integer >= 0, the code will only consider subsets of the full stem of
#     length l-substems.
        
# frozenBPs is a list of base pairs that are forced to be in each returned structure.
#     An example: frozenBPs = [[2, 9], [3, 8], [37, 46], [38, 45]]
        
# toSave is a boolean which tells youdetermines whether to save the results

# tryToLoadVariables is a boolean which tells the code whether to try to load
#     the results from a saved file (it can also work even if the saved file has a different
#     T, vs, or duplexEntropyInKB).
        
# maxSizeSTable is for the preallocation of space for STable. It tells the code how much space to preallocate
#     (in other words, what we expect the maximum total number of stems will be).
        
# onlyAllowSubsetsOfLongestStems is a boolean which should be set to True if you want to only
#     include the longest possible stems and ignore the rest. By default, the code also includes
#     any stem of length >= 16 bps. 

# onlyConsiderSubstemsFromEdges is a boolean. It should be set to True if you want to consider,
#     only substems that include the edges of the full stem (in order to simulate zipperings).
#     For example, for a stem made up of base pairs [[1,100], [2,99], [3,98], [4,97]],
#     the substems [[1,100], [2,99], [3,98]], [[2,99], [3,98], [4,97]], 
#     [[1,100], [2,99]], and [[3,98], [4,97]] will be considered, but
#     the substem [[2,99], [3,98]] will not be since it doesn't include the edge of the full stem.
        
# onlyConsiderBondedStrands is a boolean which should be set to True if you want all returned structures
#     to include at least one bond between the two strands.
        
# considerTabulatedLoops is a boolean. Set it to True to use tabulated energy/entropy parameters
#     from the Turner/Mathews groups for special hairpins and internal loops.
        
# unmatchedBPPenalty is a boolean. The standard Turner/Mathews parameters have a rule that if 
#     two nts are unpaired in the specific structure you're considering, but could pair in another structure,
#     the purine should be replaced with A and the pyrimidine with C for the purposes of the energy calculation
#     (i.e. for terminal mismatch purposes, if the terminal mismatch could in fact have paired, there's a penalty).
#     Set this to True to include that penalty, and False to ignore this rule.
        
# minNtsInHairpin is an integer. It specifies the minimum number of nts in a hairpin. It should be set to 3.

# considerAllAsTerminalMismatches is a boolean. Set it to True to ignore both
#     dangling ends and flush coaxial stacks, and to consider all ends of stems 
#     as terminal mismatches (even if the next nts over are both bound as part of different stems)
#     If this is True, includeTerminalMismatches = True; includeFlushCoaxialStacks = False;
#     and includeDanglingEnds = False.

# includeTerminalMismatches is a boolean. Set it to True to consider terminal mismatches
#     when tabulating the bondFEs. If you change this the whole CHECK function needs to be rerun.

# includeFlushCoaxialStacks is a boolean. Set it to True to consider that if two stems are adjacent,
#     the stems will stack for the purposes of considering terminal mismatches (or rather, not 
#     mismatches, in this case). See the main code for how our use of flush coaxial stacks differs
#     from the NNDB's (Turner/Mathews). If you change this the whole CHECK function needs to be rerun.
        
# includeTerminalAUATPenalties is a boolean. Set it to True to include the Turner/Mathews
#     energy and entropy penalties for AU/GU/AT base pairs at the ends of a stem. We tabulate
#     these closures regardless since it's very inexpensive, but only include them in the energy/entropy 
#     calculations if this is set to True. Therefore, this can be changed with only running
#     the postCalculationFxn.

# includeDanglingEnds is a boolean. Set it to True to use dangling ends when performing the 
#     energy/entropy calculations. We tabulate these regardless since it's relatively inexpensive,
#     but only include them in the energy/entropy calculations if this is set to True. 
#     Therefore, this can be changed with only running the postCalculationFxn.
        
# considerC3andC4 is a boolean. Set it to True if you want to calculate the three- and four-way 
#     compatibility tensors C3 and C4 in order to get rid of many structures with high-order pseudoknots.
#     This is most useful (speeds up the code most) if considerOnlyBondedStrands = False, 
#     since these tensors are not used to deduce compatibilities including stems that define 
#     intermolecular binding.

# corruptFESeed is a float. Set it to False (0) to use the Turner/Mathews parameters for
#     energies/entropies. Set it to a (nonzero) float to use that float as a seed for corrupting
#     those energy/entropy parameters.
        
# unboundButCouldBindPenalties is a list [x,y]. If two nts at the edge of a stem are unbound but could bind 
#     we may treat them (following the NNDB) as an AC pair. However, we also introduce 
#     corrections to this since that decision appears to have been arbitrary. 
#     x is the energy penalty we add, and y is the entropy penalty.
#     Even if unmatchedBPPenalty = False, this penalty is introduced if unboundButCouldBindPenalties ~=[0,0]
# =============================================================================
        
        
# =============================================================================
#         Save variables to self
# =============================================================================        
        self.version = '1_24_19'
        
        sequences, DNA, sequenceFrequencies = self.basicCorrectionsToInput(sequences, DNA, sequenceFrequencies) 
#            perform basic corrections to input like sorting sequences, making sure they're in a list, etc.
        
        self.sequences = sequences            
        self.DNA = DNA
        self.sequenceFrequencies = sequenceFrequencies
        self.b = b #in units of nucleotides (nts or ntds -- I use both interchangably)
        self.g = 3 / (2*self.b)
        self.T = T #in units of Kelvin
        self.vs = vs #in units of ntds^3
        self.duplexEntropyPenaltyInKB = duplexEntropyPenaltyInKB
        self.minBPInStem = minBPInStem
        self.numParForLoops = numParForLoops
        self.pairwiseEnergies = True
        self.storeGraphs = storeGraphs
        self.makeFigures = makeFigures
        self.printProgressUpdate = printProgressUpdate
        self.tryToLoadVariables = tryToLoadVariables
        self.allowParallelStrands = False
        self.allowPseudoknots = allowPseudoknots
        self.minNumStemsInStructure = minNumStemsInStructure
        self.substems = substems
        self.frozenBPs = frozenBPs
        self.maxSizeSTable = maxSizeSTable
        self.toSave = toSave
        self.initiationPenalty = False
        self.onlyAllowSubsetsOfLongestStems = onlyAllowSubsetsOfLongestStems
        self.onlyConsiderSubstemsFromEdges = onlyConsiderSubstemsFromEdges
        self.onlyConsiderBondedStrands = onlyConsiderBondedStrands
        self.considerTabulatedLoops = considerTabulatedLoops
        self.unmatchedBPPenalty = unmatchedBPPenalty
        self.kB = 0.0019872 #units of kcal/(mol*Kelvin)
        self.minNtsInHairpin = minNtsInHairpin #min number of nts in hairpin
        self.considerAllAsTerminalMismatches = considerAllAsTerminalMismatches
        
        if considerAllAsTerminalMismatches:
            includeTerminalMismatches = True
            includeFlushCoaxialStacks = False
            includeDanglingEnds = False
        
        self.includeTerminalMismatches = includeTerminalMismatches
        self.includeFlushCoaxialStacks = includeFlushCoaxialStacks
        self.includeDanglingEnds = includeDanglingEnds
        self.includeTerminalAUATPenalties = includeTerminalAUATPenalties
        self.considerC3andC4 = considerC3andC4
        
        if not self.allowPseudoknots:
            self.considerC3andC4 = False #since all those will already be not considered
            #However, considerC3andC4 should still be True if numSequences > 1 since we automatically only 
            #consider only intramolecular pseudoknots for C3 and C4
        
        self.corruptFESeed = corruptFESeed
        self.unboundButCouldBindPenalties = unboundButCouldBindPenalties
        self.completed = False
# =============================================================================
#         Make the list of sequences into one long sequence. 
#         Also convert the sequence to numbers
# =============================================================================
        self.sequence,self.sequenceInNumbers,self.numSequences,self.numNt,\
            self.numNt1,self.numNt2,self.linkerPos = self.multipleStrandsSetup()

        if self.numSequences == 1:
            self.onlyConsiderBondedStrands = False #doesn't make sense otherwise
        
        if self.onlyConsiderBondedStrands: #if only consider bonded strands, don't include completely unfolded structure
            self.minNumStemsInStructure = max(1, self.minNumStemsInStructure) 
        
# =============================================================================
#         Make the filename to save the data
# =============================================================================
        if self.toSave or self.tryToLoadVariables:
            self.fileVars, self.fileTxt, self.folderName, self.fileFigs = self.makeFileName()

    def basicCorrectionsToInput(self, sequences, DNA, sequenceFrequencies):
# =============================================================================
# perform basic corrections to input like sorting sequences, making sure they're in a list, etc.
# =============================================================================
        
        #If we have only one sequence and we forgot to put it in a list, 
        #and instead just used a string, we should correct it
        if not isinstance(sequences,list) and isinstance(sequences,str):
            sequences = [sequences] 

        if len(sequences) > 2:
            print('ERROR: We only accept at most two sequences at once for now') 
            #not myPrintFxn since the files haven't been created yet
            return(None)
        
        if len(DNA) < len(sequences):
            print('ERROR: len(DNA) < len(sequences)')
            return(None)
        else:
            DNA = DNA[:len(sequences)]

#       Make sure sequences is sorted so that we load the results if we tried it with opposite order
        if not sequences == sorted(sequences):
            DNA = [DNA[1], DNA[0]]
            sequences = [sequences[1], sequences[0]]
            sequenceFrequencies = [sequenceFrequencies[1], sequenceFrequencies[0]]
            
#        If we have DNA, change all U's to T's, and vice versa if we have RNA
        for i in range(len(sequences)):
            if DNA[i]:
                sequences[i].replace('U','T')
            else:
                sequences[i].replace('T','U')
            
        return(sequences, DNA, sequenceFrequencies)
        
    def calculateFELandscape(self):
# =============================================================================
#         Try to load the variables
# =============================================================================
        if self.printProgressUpdate:
            myPrintFxn('sequence = ' + self.sequence, self.fileTxt)
            myPrintFxn('minBPInStem = ' + str(self.minBPInStem), self.fileTxt)

        startTime = time.time()
        self.startTime = startTime
        
        if self.tryToLoadVariables:
            loadedVars = self.loadVars(self.folderName)
            if loadedVars:
                return(loadedVars) #don't need returnFxn here since loadVars already returns the returnFxn

# =============================================================================
#         If we couldn't load the variables (otherwise this would all be skipped)
#         Start enumerating the structures. First, the START function:
#            Enumerate the stems and define their compatibilities
# =============================================================================   
                
        self.numStems, self.STableStructure, self.STableBPs = self.createSTable()
        
        if self.printProgressUpdate:
            myPrintFxn('numStems = ' + str(self.numStems), self.fileTxt)
    
        self.frozenStems = self.frozenStemsFromFrozenBPs()
        # =============================================================================
        # given a list of base pairs that we want fixed, what regions do they correspond to?
        # frozenRegions is a cell array, where frozenRegions{i} is a vector of
        # regions containing all the base pairs of the i'th stem in frozenBPs. Thus,
        # in order to keep all of frozenBPs fixed, each possible structure must have
        # one stem from each vector in each element of frozenRegions.
        # =============================================================================
        
        #Make a list of which stems correpsond to the two strands being bound
        self.linkedStems = self.findLinkedStems() #a boolean np array of length numStems.
#           linkedStems[i] == True if stem i corresponds to the two strands being bound, and False otherwise
        
        #Make compatibility matrices
        self.C = self.makeCompatibilityMatrix()
        
        if not all(np.diagonal(self.C)): #i.e. if there are some disallowed stems (because of frozenBPs)
            self.numStems, self.frozenStems, self.STableBPs, self.STableStructure, self.C, self.linkedStems = self.removeDisallowedStems()
            if self.printProgressUpdate:
                myPrintFxn('modified numStems = ' + str(self.numStems), self.fileTxt)
        
        if self.considerC3andC4:
            self.C3 = self.makeThreewayCompatibilityTensor()
            self.C4 = self.makeFourwayCompatibilityTensor()
        else:
            self.C3 = True; self.C4 = True #so the variables exist
            
# =============================================================================
#         I had thought to calculate bond energy and entropies at the level of stems, but it's inefficient and slow.
# =============================================================================
        #Calculate the bond energy and entropy of each stem, not including what to do at the ends of the stems 
        #(i.e. dangling ends, terminal mismatches, and the like)
#        if self.pairwiseEnergies:
#            self.stemFECounts = self.calculateStemFreeEnergiesPairwise()
#            self.stemEnergies, self.stemEntropies = self.calculateStemFreeEnergiesPairwise()
#        else:
#            print('non-pairwise FE has not yet been created')
#            self.stemEnergies, self.stemEntropies = self.calculateStemFreeEnergiesSingle()
            
        if self.printProgressUpdate:
            myPrintFxn('Time for START function = ' + str(time.time() - startTime), self.fileTxt)
        
# =============================================================================
#         Next, the PERMU function: Find all possibile allowed structures by 
#         ennumerating all allowed permutations of the stems (regions) found in STable
# =============================================================================
#            We no longer do the following (explained below):
#            Taking a leaf from TT2NE's book, we'll start by defining helipoints -- 
#            These are combinations of helices separated by zero or one nts. The benefit of these
#            is that their bond energies and entropies are entirely well-defined (with the exception of
#            if we choose to include 2x1 or 2x2 internal loop parameters, but I don't 
#            think we should since those seem so wishy-washy).
# #            Update: making the helipoints is fairly fast, but making the helipoint compatibility matrix is 
# #            very very slow. The reason for that is that the number of helipoints grows quickly
# #            (almost as fast as the total number of structures) and so the helipoint compatibility matrix
# #            is extremely large. Here is the code which now we don't use anymore.
# #        helipointStartTime = time.time()
# #        self.numHelipoints, self.helipoints, self.adjStems = self.makeHelipoints()
# #        self.CHeli = self.makeHelipointCompatibilityMatrix()
# #        #self.helipointEnergies, self.helipointEntropies = self.calculateHelipointFreeEnergies()
# #        
# #        #Define helipoint energies
# #        helipointTime = time.time() - helipointStartTime
# #        if self.printProgressUpdate:
# #            print('Time for helipoint function = ' + str(helipointTime))
# =============================================================================
# =============================================================================
        
        permuFxnStartTime = time.time()
        self.numStructures, self.structures, self.linkedStructures = self.enumerateStructures()
        
        if self.printProgressUpdate:
            myPrintFxn('numStructures = ' + str(self.numStructures), self.fileTxt)
            myPrintFxn('Time for PERMU function = ' + str(time.time() - permuFxnStartTime), self.fileTxt)

        checkFxnStartTime = time.time()
        self.checkFxnStartTime = checkFxnStartTime
        (self.bondFECounts, self.dangling5Count, self.dangling3Count, 
         self.allNumVs, self.allComponentGraphs, self.structureGraphList, 
         self.allWhichStructureGraph,self.unboundButCouldBindCounts) = self.calculateFreeEnergies()
        
        if self.printProgressUpdate:
            myPrintFxn('Time for CHECK function = ' + str(time.time() - checkFxnStartTime), self.fileTxt)

        return(self.postCalculationFxn())
    
    
    def makeThreewayCompatibilityTensor(self):
# =============================================================================
# In order to avoid having any numerical integration, only consider pseduoknots of max order 2,
# meaning only containing maximum two stems (for example don't consider intramolecular kissing hairpins)
# If we have three (antiparallel) stems each defined by a base pair: A bound to a, B to b, C to c (A<a, B<b, C<c)
# such that A<B<C, what matters for whether or not they can all coexist is their respective order.
# The first has to be A. Then, either a or B can be next. But if it's a, then that makes a hairpin
# and whatever happens next will result in a pseudoknot of max order 2.
        
# So the game is to come up with an ordering for A,a,B,b,C,c -- or, more intuitively, using parentheses, of
# the six characters ([{}]) such that ( always comes before [ which always comes before {,
# such that no subset of consecutive brackets correspond to a real substructure. 
# For example, in the sequence ([{]}), the subset '[{]}' corresponds to a real substructure (it is
# a complete set of opening and closing brackets).
        
# So ( has to be first, then ) or [, but as discussed, in coming up with disallowed structures, next has to be [.
# After ([ next can be either ) or {  -- if it's ] that closes that bracket and the pair can be removed.
# Let's consider these one by one. 
#    First, consider ([). Next can't be ] since that would be a complete structure on its own.
#    So it has to go ([){. Next can't come }, it has to be ].
#    So we have our first disallowed set: ([){]}, or ABaCbc. The order of A,a,B,b,C,c is [0, 2, 1, 4, 3, 5]
#    
#    Next, consider ([{. After that can be either ] or ), so we consider both.
#        So we have ([{]. Next can't be } since [{]} is a complete substructure, 
#        so we have our next disallowed set: ([{])}, or ABCbac. The order of A,a,B,b,C,c is [0, 2, 4, 3, 1, 5]
#        
#        Next we have ([{). Either possibilities after are acceptable, so our final two disallowed sets are:
#        ([{)]} or ABCabc -- The order of A,a,B,b,C,c is [0, 2, 4, 1, 3, 5]
#        and ([{)}] or ABCacb -- The order of A,a,B,b,C,c is [0, 2, 4, 1, 5, 3]
#        
# This argument is where we get threewayDisallowedPermutations, the four sorted orders of 
#        A,a,B,b,C,c for the four disallowed sets.
# =============================================================================
        numStems = self.numStems
        STableBPs = self.STableBPs
        CMat = self.C
        linkedStems = self.linkedStems #a boolean np array of length numStems.
#        linkedStems[i] == True iff stem i corresponds to a bond between two strands
        
        threewayDisallowedPermutations = [[0, 2, 1, 4, 3, 5], [0, 2, 4, 3, 1, 5], [0, 2, 4, 1, 3, 5], [0, 2, 4, 1, 5, 3]]
        
        C3 = np.zeros((numStems,numStems,numStems),dtype=bool)
#        Go through each triplet of stems
        for i in range(numStems):
            if not linkedStems[i]: #things that appear to be high-order pseudoknots because of the linker of 'OOO's
#                may in reality be pseudoknots of much lower order
                for j in range(i+1,numStems):
                    if CMat[i,j] and not linkedStems[j]: #they need to be pairwise compatible
                        for k in range(j+1,numStems):
                            if CMat[j,k] and not linkedStems[k]:
                                A = STableBPs[i][0][0] #first nt of stem
                                a = STableBPs[i][0][1] #what it's bound to. a>A
                                B = STableBPs[j][0][0]
                                b = STableBPs[j][0][1] #b>B
                                C = STableBPs[k][0][0]
                                c = STableBPs[k][0][1] #c>C
                                
                                #put the stems into order such that a<c<e
                                l = sorted([[A,a],[B,b],[C,c]]) #python automatically sorts by first element in each sublist
                                A, a = l[0]
                                B, b = l[1]
                                C, c = l[2]
                                
                                #What matters for determining whether this triplet is allowed is the order of a,b,c,d,e,f
                                if list(np.argsort([A,a,B,b,C,c])) not in threewayDisallowedPermutations:
                                    C3[i,j,k] = True; C3[i,k,j] = True; C3[j,i,k] = True 
                                    C3[j,k,i] = True; C3[k,i,j] = True; C3[k,j,i] = True
                        
        return(C3)
        
    def makeFourwayCompatibilityTensor(self):
# =============================================================================
# Similarly to three-way compatibility, we also consider fourway compatibility.
# There are far more than four disallowed ways of having four stems, and we'll show how you can enumerate them here.
# We have A,B,C,D,a,b,c,d, such that A<a, B<b, etc. and A<B<C<D. We want to find how we can arrange the 8 bps
# (within these rules) such that no subsequence forms a complete structure (contains all its own opening and closing nts
# where a is the closing nt of A, etc.)
#        
# We have to start with AB. Next can be a or C.
#   ABa. Next can't be b so it has to be C
#   ABaC. Next can be b or D. 
#       ABaCb. Next can only be D. Then has to come c then d. --> ABaCbDcd
#       ABaCD. Next can be b or c
#          ABaCDbcd, ABaCDbdc
#          ABaCDcbd, ABaCDcdb
#   ABC. Next can be a, b, or D. If it's D, we can have any combination of the closing pairs as long as d isn't next and a isn't last.
#        ABCDabcd, ABCDabdc, ABCDacbd, ABCDacdb, ABCDadbc, ABCDadcb
#        ABCDbcad, ABCDbdac, ABCDbacd, ABCDbadc
#        ABCDcbad, ABCDcdab, ABCDcabd, ABCDcadb
#    
#       ABCa. Next can be b,c, or D
#           ABCabDcd
#           ABCacDbd
#           ABCaD. Next can be b or c
#           ABCaDbcd, ABCaDbdc, ABCaDcbd, ABCaDcdb
#       ABCb. Next can be a, c, or D (this is as above, but a can't be the last one)
#           ABCbaDcd, ABCbcDad, ABCbDacd, ABCbDadc, ABCbDcad
# I checked this against the MatLab code -- besides parallel stems, the difference is that the MatLab code also
# includes those sets of four stems that lead to a pseduoknot of order 3 (we're only considering those of order 4 here)
# =============================================================================
        numStems = self.numStems
        STableBPs = self.STableBPs
        CMat = self.C
        C3 = self.C3
        printProgressUpdate = self.printProgressUpdate
        linkedStems = self.linkedStems #a boolean np array of length numStems.
#        linkedStems[i] == True iff stem i corresponds to a bond between two strands
        
        fourwayDisallowedPermutationsStr = ['ABaCbDcd','ABaCDbcd','ABaCDbdc','ABaCDcbd','ABaCDcdb','ABCDabcd', 
                                            'ABCDabdc', 'ABCDacbd', 'ABCDacdb', 'ABCDadbc', 'ABCDadcb','ABCDbcad', 
                                            'ABCDbdac', 'ABCDbacd', 'ABCDbadc','ABCDcbad', 'ABCDcdab','ABCDcabd', 
                                            'ABCDcadb','ABCabDcd','ABCacDbd','ABCaDbcd', 'ABCaDbdc', 'ABCaDcbd', 
                                            'ABCaDcdb','ABCbaDcd', 'ABCbcDad', 'ABCbDacd', 'ABCbDadc', 'ABCbDcad']
        disallowedPermuStr2Num = {'A':0, 'a':1, 'B':2, 'b':3, 'C':4, 'c':5, 'D':6, 'd':7}
        
        fourwayDisallowedPermutations = [[disallowedPermuStr2Num[i] for i in s] for s in fourwayDisallowedPermutationsStr]
        
        C4 = np.zeros((numStems,numStems,numStems,numStems),dtype=bool)
#        Go through each quadruplet of stems
        counter = 0
        
        for i in range(numStems):
            if not linkedStems[i]:
                for j in range(i+1,numStems):
                    if CMat[i,j] and not linkedStems[j]: #they need to be pairwise compatible
                        for k in range(j+1,numStems):
                            if C3[i,j,k] and not linkedStems[k]: #they need to be pairwise compatible and three-way compatible
                                for q in range(k+1,numStems):
                                    if C3[i,j,q] and C3[i,k,q] and C3[j,k,q] and not linkedStems[q]: 
                                            #C3 includes pairwise compatibility as well as three-way compatibility
                                        A = STableBPs[i][0][0] #first nt of stem
                                        a = STableBPs[i][0][1] #what it's bound to. a>A
                                        B = STableBPs[j][0][0]
                                        b = STableBPs[j][0][1] #b>B
                                        C = STableBPs[k][0][0]
                                        c = STableBPs[k][0][1] #c>C
                                        D = STableBPs[q][0][0]
                                        d = STableBPs[q][0][1]
                                        
                                        #put the stems into order such that a<c<e
                                        l = sorted([[A,a],[B,b],[C,c],[D,d]]) #python automatically sorts by first element in each sublist
                                        A, a = l[0]
                                        B, b = l[1]
                                        C, c = l[2]
                                        D, d = l[3]
                                        
                                        #What matters for determining whether this triplet is allowed is the order of a,b,c,d,e,f
                                        if list(np.argsort([A,a,B,b,C,c,D,d])) not in fourwayDisallowedPermutations:
                                            C4[i,j,k,q] = True; C4[i,j,q,k] = True; C4[i,k,j,q] = True; C4[i,k,q,j] = True
                                            C4[i,q,j,k] = True; C4[i,q,k,j] = True; C4[j,i,k,q] = True; C4[j,i,q,k] = True
                                            C4[j,k,i,q] = True; C4[j,k,q,i] = True; C4[j,q,i,k] = True; C4[j,q,k,i] = True
                                            C4[k,i,j,q] = True; C4[k,i,q,j] = True; C4[k,j,i,q] = True; C4[k,j,q,i] = True
                                            C4[k,q,i,j] = True; C4[k,q,j,i] = True; C4[q,i,j,k] = True; C4[q,i,k,j] = True
                                            C4[q,j,i,k] = True; C4[q,j,k,i] = True; C4[q,k,i,j] = True; C4[q,k,j,i] = True
                                            
                                            counter += 1
                                            if counter % 5e5 == 0 and printProgressUpdate:
                                                myPrintFxn('Setting up C4: i = ' + str(i) + '; length c4 indices = ' + str(counter), self.fileTxt)
                            
        return(C4)
    
    def postCalculationFxn(self):
# =============================================================================
#         Take the energies and entropies that have been calculated, add the 
#         relevant vs factors to the entropies, calculate free energies
#         calculate probabilities, sort the structures, calculate prob of bonded/unbonded, etc.
# =============================================================================
        postCalculationFxnTime = time.time()
        
        numStructures = self.numStructures
        numSequences = self.numSequences
        linkedStructures = self.linkedStructures
        DNA = self.DNA
        corruptFESeed = self.corruptFESeed
        g = self.g #gamma
        T = self.T
        duplexEntropyPenaltyInKB = self.duplexEntropyPenaltyInKB
        allComponentGraphs = self.allComponentGraphs
        allNumVs = self.allNumVs
        bondFECounts = self.bondFECounts
        dangling5Count = self.dangling5Count
        dangling3Count = self.dangling3Count
        unboundButCouldBindCounts = self.unboundButCouldBindCounts
        
        logVs = np.log(self.vs)
        kB = self.kB
        printProgressUpdate = self.printProgressUpdate
        storeGraphs = self.storeGraphs
        unboundButCouldBindPenalties = self.unboundButCouldBindPenalties
        
# =============================================================================
# ********* Perform the loop entropy calculation given the component graphs ********
#     each componentGraphs is a dictionary, where each entry of the dictionary yields a vector of vectors
#     where each element of the vector is a different parameter set for a different instance of that net.
#     The parameter sets are defined as vectors in the order [s1, s2, s3, s4, l1, l2].
#     If the parameter is not defined for the graph its value is set to be -1 (like s4 for an open-net-2a)
#     Example usage: componentGraphs['closedNet2a'][0][3] gives the value of s4 for the first
#     instance of closedNet2a found in the structure. componentGraphs has these keys: 
#    ['openNet0','closedNet0','openNet1','closedNet1','openOpenNet2','openNet2a',
#    'closedNet2a','openNet2b','closedNet2b','openNet2c','closedNet2c','tooComplex']
# =============================================================================

        compGraphKeys = ['openNet0','closedNet0','openNet1','closedNet1','openOpenNet2','openNet2a',
                         'closedNet2a','openNet2b','closedNet2b','openNet2c','closedNet2c','tooComplex']
        allLoopEntropies = np.zeros(numStructures)
        
        for whichStruct in range(numStructures):
            loopEntropy = 0
            componentGraphs = allComponentGraphs[whichStruct]
            for netType in compGraphKeys:
                for slVec in componentGraphs[netType]:
                    loopEntropy += entropyOfNet(netType, slVec, g)
            allLoopEntropies[whichStruct] = kB*(loopEntropy + allNumVs[whichStruct]*logVs)
            
            
# =============================================================================
#   Perform bond energy and bond entropy calculations (including dangling ends)
# =============================================================================
        if not corruptFESeed:
            energyMatrices, entropyMatrices = bondFreeEnergies()
        else:
            energyMatrices, entropyMatrices = corruptBondFreeEnergies(corruptFESeed)
        
        allBondEnergies = np.zeros(numStructures)
        allBondEntropies = np.zeros(numStructures)
#        The fact that these matrices are sparse slows down the code (for ~3000 structures
#       for non-sparse matrices postCalculationFxn took 0.17 s which changed to 0.8 s after introducing sparse matrices)
#        We cut this in half by only considering RNA(DNA) for systems that only have RNA(DNA)
        bpTypeList = []
        if any([not i for i in DNA]):
            bpTypeList.append(0) #RNA/RNA bonds
        if any(DNA):
            bpTypeList.append(1) #DNA/DNA bonds
        if bpTypeList == [0,1]:
            bpTypeList.append(2) #RNA/DNA bonds
        if self.includeTerminalAUATPenalties:
            bpTypeList.append(3) #terminal AU/GU/AT penalties
            
        for bpType in bpTypeList: #RNARNA, DNADNA, RNADNA, terminalAUATPenalty (don't need unknown since it just gives 0)
            allBondEnergies += np.sum(np.sum(np.multiply(energyMatrices[bpType], 
                                np.array([bondFECounts[whichStruct][bpType].todense() 
                                for whichStruct in range(numStructures)])),axis=1),axis=1)
            allBondEntropies += np.sum(np.sum(np.multiply(entropyMatrices[bpType], 
                                 np.array([bondFECounts[whichStruct][bpType].todense() 
                                 for whichStruct in range(numStructures)])),axis=1),axis=1)
        
        bpTypeList = []
        if any([not i for i in DNA]):
            bpTypeList.append(0) #RNA
        if any(DNA):
            bpTypeList.append(1) #DNA
            
        if self.includeDanglingEnds:
            if not corruptFESeed:
                dangling5Energies, dangling5Entropies, dangling3Energies, dangling3Entropies = danglingEndMatrices()
            else:
                dangling5Energies, dangling5Entropies, dangling3Energies, dangling3Entropies = corruptDanglingEndMatrices(corruptFESeed)
            
            for bpType in bpTypeList: #RNA, DNA, (don't need unknown)
                allBondEnergies += np.sum(np.sum(np.multiply(dangling5Energies[bpType], 
                                np.array([dangling5Count[whichStruct][bpType].todense() 
                                for whichStruct in range(numStructures)])),axis=1),axis=1)
                allBondEnergies += np.sum(np.sum(np.multiply(dangling3Energies[bpType], 
                                np.array([dangling3Count[whichStruct][bpType].todense() 
                                for whichStruct in range(numStructures)])),axis=1),axis=1)
                allBondEntropies += np.sum(np.sum(np.multiply(dangling5Entropies[bpType], 
                                 np.array([dangling5Count[whichStruct][bpType].todense() 
                                 for whichStruct in range(numStructures)])),axis=1),axis=1)
                allBondEntropies += np.sum(np.sum(np.multiply(dangling3Entropies[bpType], 
                                 np.array([dangling3Count[whichStruct][bpType].todense() 
                                 for whichStruct in range(numStructures)])),axis=1),axis=1)
        
        if numSequences > 1:
            self.allDuplexEntropies = linkedStructures * kB * duplexEntropyPenaltyInKB
        else:
            self.allDuplexEntropies = np.zeros(numStructures)
        
        if any(unboundButCouldBindPenalties):
            allBondEnergies += unboundButCouldBindCounts * unboundButCouldBindPenalties[0]
            allBondEntropies += unboundButCouldBindCounts * unboundButCouldBindPenalties[1]
            
        self.allLoopEntropies = allLoopEntropies
        self.allBondEnergies = allBondEnergies
        self.allBondEntropies = allBondEntropies
        
        allFreeEnergies = allBondEnergies - T * (allBondEntropies + allLoopEntropies + self.allDuplexEntropies)
        
        probabilities = np.exp(-allFreeEnergies / (kB*T))
        partitionFxn = np.sum(probabilities)
        probabilities /= partitionFxn
        
        self.indexSort = np.flip(np.argsort(probabilities)) #flip so that probabilities are decreasing
#        The minFE structure will be structures[indexSort[0]]
#        More generally, the structure with the n^th lowest FE will be structures[indexSort[n]]
#        and its FE will be given by sortedFEs[n]
        self.sortedProbs = probabilities[self.indexSort]
        self.sortedFEs = allFreeEnergies[self.indexSort]
        
        if storeGraphs:
            structureGraphList = self.structureGraphList
            allWhichStructureGraph = self.allWhichStructureGraph
            numGraphs = len(structureGraphList)
            graphProbs = np.array([np.sum([probabilities[i] for i in range(numStructures) if 
                                           allWhichStructureGraph[i] == j]) for j in range(numGraphs)])
            self.indexSortedGraphProbs = np.flip(np.argsort(graphProbs))
            self.sortedGraphProbs = graphProbs[self.indexSortedGraphProbs]
        else:
            self.indexSortedGraphProbs = []
            self.sortedGraphProbs = []
            
        if numSequences > 1:
            probBonded = np.sum([probabilities[i] for i in range(numStructures) if linkedStructures[i]])
            if printProgressUpdate:
                myPrintFxn('probability of sequences being bound = ' + str(probBonded), self.fileTxt)
        else:
            probBonded = 0
        self.probSeparateStrands = 1 - probBonded
        
        if printProgressUpdate:
            myPrintFxn('Time for postCalculationFxn function = ' + str(time.time() - postCalculationFxnTime), self.fileTxt)
        
        self.completed = True
        
        if self.toSave: #saving can take some time
            save(self,self.fileVars)
            
        return(self.returnFxn())   
        
        
    
    def calculateFreeEnergies(self):
        numStructures = self.numStructures
        numParForLoops = self.numParForLoops
        considerTabulatedLoops = self.considerTabulatedLoops
        storeGraphs = self.storeGraphs
        STableStructure = self.STableStructure
        structures = self.structures
        
# =============================================================================
#         Initialize the tables we're going to return
# =============================================================================
        bondFECounts = [None] * numStructures
        dangling5Count = [None] * numStructures
        dangling3Count = [None] * numStructures
        allComponentGraphs = [None] * numStructures
        allNumVs = np.zeros(numStructures,dtype=int)
        unboundButCouldBindCounts = np.zeros(numStructures,dtype=int)
        
        if considerTabulatedLoops:
            allSpecialLoopEnergies = np.zeros(numStructures)
            allSpecialLoopEntropies = np.zeros(numStructures)
        else:
            allSpecialLoopEnergies = []
            allSpecialLoopEntropies = []
                
        #keep track of weight matrices and numBonds matrices, and which of 
        #the graphs in weightMatrixList,numBondsMatrixList each structure corresponds to
        structureGraphList = []
        if storeGraphs:
            allWhichStructureGraph = np.zeros(numStructures,dtype=int) 
        else:
            allWhichStructureGraph = []
        
        if numParForLoops > 0:
            listOfAllNewStructureGraphs = [None]*numStructures
        
# =============================================================================
#         Create various arrays once so we don't need to do it for each structure over again
# =============================================================================
        allPerms = []
        for i in range(2,6):
            allPerms.append(list(itertools.permutations(range(1,i))))
        
#        anyDNA = any(DNA)
        self.refNets = createRefNets()
#        [specialHairpins, IL11FE, IL11E, IL12FE, IL12E, IL22FE, IL22E] = makeVars(anyDNA,considerTabulatedLoops);


# =============================================================================
#         Iterate over structures, calculating the free energy for each
# =============================================================================
        #Parallelize the next part: Also, Maybe use map?
        
        storedGraphEigs = []
        #Potentially get rid of graphs if they are showing that they don't have a high probability of occurring (add cutoff as input var)
        for whichStruct in range(numStructures):
            structure = [None]*len(structures[whichStruct])
            for i in range(len(structure)):
                structure[i] = STableStructure[structures[whichStruct][i]]
            
            allNumVs[whichStruct], allComponentGraphs[whichStruct], structureGraph = self.calculateLoopEntropy(structure)
            (bondFECounts[whichStruct], dangling5Count[whichStruct], dangling3Count[whichStruct],
             specialLoopEnergies, specialLoopEntropies,unboundButCouldBindCounts[whichStruct]) = \
                 self.calculateBondEnergy(structure, structureGraph)
            
            if considerTabulatedLoops:
                allSpecialLoopEnergies[whichStruct] = specialLoopEnergies
                allSpecialLoopEntropies[whichStruct] = specialLoopEntropies
                
            if storeGraphs: #faster_could_be_isomorphic barely adds time to 
#                case without storeGraphs, but makes a lot of mistakes (thinks graphs are isomorphic when they're not)

                eigGraph = np.sort(np.linalg.eigvals(nx.convert_matrix.to_numpy_matrix(structureGraph)))
                GisoTest = [i for i in range(len(structureGraphList)) 
                    if nx.faster_could_be_isomorphic(structureGraph,structureGraphList[i]) and
                    all(np.abs(eigGraph - storedGraphEigs[i]) < 1e-6)]
                Giso = []
                for i in GisoTest: #break makes this ~2x faster than list comprehension, 
#                    but still much slower than not having storeGraphs
#                   Considering eigenvalue equality (above) cuts time in half again
                    if nx.is_isomorphic(structureGraph,structureGraphList[i]):
                        Giso = [i]
                        break
                if not len(Giso):
                    structureGraphList.append(structureGraph)
                    storedGraphEigs.append(eigGraph)
                    allWhichStructureGraph[whichStruct] = len(structureGraphList) - 1
                else: #which structureGraph does it correspond to?
                    allWhichStructureGraph[whichStruct] = Giso[0]
            
            if whichStruct % 100000 == 0 and whichStruct > 0:
                myPrintFxn('Calculating free energies: whichStruct = ' + str(whichStruct), self.fileTxt)
                myPrintFxn('Time elapsed since start of CHECK = ' + str(time.time() - self.checkFxnStartTime), self.fileTxt)
                #myPrintFxn('', self.fileTxt) #add an extra line for clarity
        return(bondFECounts, dangling5Count, dangling3Count, allNumVs, allComponentGraphs, 
               structureGraphList, allWhichStructureGraph,unboundButCouldBindCounts)

    
    def calculateBondEnergy(self, structure, G):
        sequenceInNumbers = self.sequenceInNumbers
        numNt = self.numNt
        linkerPos = self.linkerPos
        unmatchedBPPenalty = self.unmatchedBPPenalty

        if len(linkerPos):
            firstNts = [0, linkerPos[-1] + 1]
            lastNts = [linkerPos[0] - 1, numNt - 1]
        else:
            firstNts = [0]
            lastNts = [numNt - 1]

        boundNts = structure2boundNt(structure,numNt)
        bpsList = structure2bpsList(structure)
        
        RNARNACount = scipy.sparse.lil_matrix((6,16),dtype=int) #np.zeros((6,4,4),dtype = int)
        #First index tells you if the first bp of the set is AU (0) CG (1) GC (2) UA (3) GU (4) or UG (5)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3) (which row of table 1 in Serra & Turner).
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3) (which column of table 1 in Serra & Turner).
            
#    Had to make this and the others 2D array to be able to use scipy sparse functions, so 
#    second/third index is replaced by (4*second index) + third index. That's why it's
#        np.zeros(6,16) and not np.zeros(6,4,4).
#        lil_matrix was chosen because, from scipy reference page, it is fast for constructing
#        sparse matrices incrementally. For operations (later we'll multiply it) it should be converted to
#        another form.

        DNADNACount = scipy.sparse.lil_matrix((4,16),dtype=int) #np.zeros((4,16),dtype = int)
        #First index tells you if the first bp of the set is AT (0) CG (1) GC (2) or TA (3)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or T(3)
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or T(3)
            
        RNADNACount = scipy.sparse.lil_matrix((8,16),dtype=int) #np.zeros((8,16),dtype = int)
        #First index tells you if the set is (putting RNA first) 
        #5'AX3'/3'TY5' (0) 5CX3/3GY5 (1) 5GX3/3CY5 (2) 5UX3/3AY5 (3) 5XA3/3YT5 (4) 5XC3/3YG5 (5) 5XG3/3YC5 (6) or 5XU3/3YA5 (7).
        #Second index tells you if X is A (0) C (1) G (2) or U (3)
        #Third index tells you if Y is A (0) C (1) G (2) or T (3)
        
        terminalAUATCount = scipy.sparse.lil_matrix((1,2),dtype=int)
#        Number of terminal AU (or GU) pairs, number of terminal AT pairs
        
        unknownCount = scipy.sparse.lil_matrix((1,1),dtype=int) #np.zeros((1,1), dtype = int)
        
        bondFECounts = [RNARNACount, DNADNACount, RNADNACount, terminalAUATCount, unknownCount]
        
        dangling5RNARNACount = scipy.sparse.lil_matrix((4,4),dtype=int) #np.zeros((4,4),dtype = int)
        dangling5DNADNACount = scipy.sparse.lil_matrix((4,4),dtype=int) #np.zeros((4,4),dtype = int)
        dangling3RNARNACount = scipy.sparse.lil_matrix((4,4),dtype=int) #np.zeros((4,4),dtype = int)
        dangling3DNADNACount = scipy.sparse.lil_matrix((4,4),dtype=int) #np.zeros((4,4),dtype = int)
        unknownDanglingCount = scipy.sparse.lil_matrix((1,1),dtype=int) #np.zeros((1,1),dtype = int)
        
        dangling5Count = [dangling5RNARNACount, dangling5DNADNACount, unknownDanglingCount]
        dangling3Count = [dangling3RNARNACount, dangling3DNADNACount, unknownDanglingCount]
#        first index tells you if the paired ntd is an A(1) C(2) G(3) U/T(4) 
#        second index tells you if the dangling ntd is A(1) C(2) G(3) or U/T(4).
#        Uses dangling nt to determine if it's RNA/RNA or DNA/DNA

        specialLoopEnergies = 0; specialLoopEntropies = 0
        
        unboundButCouldBindCounts = 0
        
        if self.considerTabulatedLoops: #XXX NOT CREATED YET
            #check if any hairpin loops are special
            for edge in G.selfloop_edges(): 
                closingBPs = G.nodes[edge[0]]['bp'] #edge[0] and edge[1] are equal by definition
                hairpinLoopSeq = sequenceInNumbers[closingBPs[0]:closingBPs[1]+1]
                if (not any(hairpinLoopSeq == 0) #can't have unknown ntds or multiple strands bonded
                        and len(hairpinLoopSeq) <= 8): #%none of the special hairpin loops have more than 6 unpaired ntds
                    additionalLoopEnergy, additionalLoopEntropy = specialHairpinFxn(hairpinLoopSeq)
                    specialLoopEnergies += additionalLoopEnergy
                    specialLoopEntropies += additionalLoopEntropy
                    
#        from Zhang & Chen (2006) which cites Turner ,... Freier (1988) 
#        the energy gained from a single bp is zero (in the nearest neighbor model)
        stemCounter = 0 #remove stems of length 1 for the nearest-neighbor model energy calculation
        while stemCounter < len(structure): 
            if len(structure[stemCounter]) == 2:
                del structure[stemCounter]
            else:
                stemCounter += 1
        
#        #If we previously (at the level of stems) made the stemFECounts matrix
#        for stemIndex in self.structures[whichStruct]: #bondFE from the stems
#            for bpType in range(5): #RNA/RNA, DNA/DNA, RNA/DNA, terminalAU/AT, unknown
#                bondFECounts[bpType] += self.stemFECounts[stemIndex][bpType]
                
        for stem in structure: #base pairs in stem
            numBonds = int(len(stem)/2)
            for j in range(numBonds - 1):
                firstNtIndex = stem[j]
                firstBPIndex = stem[j+numBonds]
                secondNtIndex = stem[j+1]
                secondBPIndex = stem[j+1+numBonds]
                index, bpType = freeEnergyMatrixIndices(sequenceInNumbers, firstNtIndex,firstBPIndex,
                                                            secondNtIndex,secondBPIndex, bound = [1,1],
                                                            unmatchedBPPenalty = unmatchedBPPenalty)
                bondFECounts[bpType][index[0], index[1]] += 1
            
#            Terminal AU/GU/AT penalties (bpType = 3). Tabulate them no matter what,
#             but we only include them in the FE calculation if includeTerminalAUATPenalties==True.
            for j in [0, numBonds - 1]: #penalties apply to ends of helices
                if sequenceInNumbers[stem[j]] == 4 or sequenceInNumbers[stem[j + numBonds]] == 4:
                    bondFECounts[3][0,0] += 1 #terminal AU/GU was found
                elif sequenceInNumbers[stem[j]] == 8 or sequenceInNumbers[stem[j + numBonds]] == 8:
                    bondFECounts[3][0,1] += 1 #terminal AT was found
                
                

        for stem in structure: #two potential terminal mismatches for each stem
            numBonds = int(len(stem)/2)
            
            for c, (firstNtIndex, firstBPIndex) in enumerate(zip(
                    [stem[numBonds], stem[numBonds-1]], [stem[0], stem[-1]])):

                secondNtIndex = firstNtIndex + 1 
                secondBPIndex = firstBPIndex - 1
                bound = [1, 0] #secondNt is not actually bound to secondBP

# =============================================================================
#                   Terminal mismatches
#   For the terminal mismatch just after stem ends, firstNt = stem[numBonds - 1], firstBP = stem[-1]
#   For the mismatch just before stem starts, secondNt = stem[0] and secondBP = stem[numBonds]
#   but if you flip it upside down, you get firstNt = stem[numBonds], firstBP = stem[0]
#   
#   Alternatively:
#    For the terminal mismatch just before stem starts, secondNt = stem[0] and secondBP = stem[numBonds]
#    For the mismatch just after stem ends, firstNt = stem[numBonds - 1], firstBP = stem[-1]
#    but if you flip it upside down, you get secondBP = stem[numBonds - 1], secondNt = stem[-1]
#
#    don't need to worry if secondNtIndex or secondBPIndex are not real nts since freeEnergyMatrixIndices takes care of that
# =============================================================================
                
                if ((secondNtIndex not in boundNts and secondBPIndex not in boundNts
                     and self.includeTerminalMismatches) or self.considerAllAsTerminalMismatches):
                    index, bpType = freeEnergyMatrixIndices(sequenceInNumbers, firstNtIndex,firstBPIndex,
                                                            secondNtIndex,secondBPIndex, bound = bound,
                                                            unmatchedBPPenalty = unmatchedBPPenalty)
                    bondFECounts[bpType][index[0], index[1]] += 1
                    
                    #If the two nts are complementary and could be bound, add one to unboundButCouldBindCounts
                    if not (secondNtIndex < 0 or secondNtIndex > numNt - 1 or secondBPIndex < 0 or secondBPIndex > numNt - 1):
                        if isComplementary(sequenceInNumbers[secondNtIndex], sequenceInNumbers[secondBPIndex]):
                            unboundButCouldBindCounts += 1
                    
# =============================================================================
#                        Flush coaxial stacks
# if there's a bulge loop of length 1, we also need to consider the stacking energy of adjacent bp
# see footnote to Table 3 in A set of nearest neighbor parameters for predicting the enthalpy change 
# of RNA secondary structure formation by Zhi John Lu  Douglas H. Turner  David H. Mathews
#                        
# We want to actually consider any bulge loop to have a flush coaxial stack.
# The reason for that is that the Turner parameters have a lower free energy 
# cost of forming bulge loops than hairpin loops, but we treat them the same
# (as closed nets 0). Therefore, while the Turner parameters don't give bulge loops 
# the benefitial energy of, say, a flush coaxial stack, we do allow that.
#                        
# This doesn't agree with Turner methodology though in that if we have 
# a 3-way junction, where each stem is flush with the next,
# we'll consider two flush coaxial stacks here, while Turner methodology
# says you should only consider the best one (and treat the other as a dangling end)
# =============================================================================
                elif secondNtIndex in boundNts and self.includeFlushCoaxialStacks: 
                    #This considers all possible flush coaxial stacks
#                        Don't want to double count these so only consider them from one side 
#                        (i.e. we don't consider if secondBPIndex in boundNts)
                
                    #find what secondNtIndex is bound to
                    secondBPIndex = [bp[1-x] for bp in bpsList for x in [0,1] if bp[x] == secondNtIndex][0]
                    
                    index, bpType = freeEnergyMatrixIndices(sequenceInNumbers, firstNtIndex,firstBPIndex,
                                                            secondNtIndex,secondBPIndex, bound = bound,
                                                            unmatchedBPPenalty = unmatchedBPPenalty)
                    bondFECounts[bpType][index[0], index[1]] += 1
                    
# =============================================================================
#                     Dangling ends
#    We only consider dangling ends at edges of sequence.
#    This doesn't agree with Turner methodology as discussed above.
#    There are four possible dangling ends (we use N1 to mean firstNt, and B1 to mean firstBP):
#                        
#   5' dangling end       3' dangling end      5' dangling end       3' dangling end
# 5'----B2--B1----3'              B1----3'    5'----N1              5'----N1--N2----3'
#           N1----5'    3'----N2--N1----5'    3'----B1--B2----3'    3'----B1
#    N1 is last nt          B1 is first nt      N1 is last nt          B1 is first nt
                    
#  So 5' dangling ends have N1 as the last nt and 3' dangling ends have B1 as the first nt
# (When you orient yourself so that N1 and B1 are always bound to one another)
                    
# We tabulate these regardless, but only include them in the FE calculation (in the postCalculationFxn)
# if includeDanglingEnds = True
# =============================================================================
                if firstNtIndex in lastNts: #Don't need to worry if lastNt is bound to firstNt
#                        since that'll be taken care of by danglingMatrixIndices
                    index, bpType = danglingMatrixIndices(sequenceInNumbers, firstNtIndex, secondNtIndex)
                    dangling5Count[bpType][index[0], index[1]] += 1
                if firstBPIndex in firstNts: 
                    index, bpType = danglingMatrixIndices(sequenceInNumbers, firstBPIndex, secondBPIndex)
                    dangling3Count[bpType][index[0], index[1]] += 1

        return(bondFECounts,dangling5Count,dangling3Count,specialLoopEnergies,specialLoopEntropies,
               unboundButCouldBindCounts)
    
    def calculateLoopEntropy(self,structure):
        G, numVs, bridgeList = createGraphFromStructure(structure, self.numNt, self.linkerPos)
        componentGraphs = calculateEntropyFromGraph(G, self.refNets, bridgeList)        
        
        return(numVs, componentGraphs, G)
        
    def enumerateStructures(self):
# =============================================================================
#        Make the list of structures. Each structure is defined by the list of stems that comprise it
#        (previously helipoints)
# =============================================================================
        numStems = self.numStems
        C = self.C; C3 = self.C3; C4 = self.C4
        printProgressUpdate = self.printProgressUpdate
        considerC3andC4 = self.considerC3andC4
        minNumStemsInStructure = self.minNumStemsInStructure
        linkedStems = self.linkedStems
        numSequences = self.numSequences
        onlyConsiderBondedStrands = self.onlyConsiderBondedStrands
        
        numStructures = 0
        structures = []
        
        prevNumStructures = -1 #keep track of this just so that we don't print the same statement multiple times
 
        for i in range(numStems):
            for j in range(i+1,numStems):
                if C[i,j]:
                    currStructure = [i,j]                    
                    lenCurrStructure = 2
                    k = j #the next stem we'll try adding (we're about to add one so that's why it's j and not j+1)
                    while lenCurrStructure >= 2:
                        while k < numStems - 1:
                            k += 1
                            
                            mutuallyCompatible = True
                            #Check mutual 2-way compatibility between the stem we want to add and all stems 
                            #present in the current structure.
                            
                            for l in currStructure[:lenCurrStructure]: #the same as just "in currStructure" but written
                                #this way so that if we want to preallocate space it'll be easier
                                if not C[k,l]:
                                    mutuallyCompatible = False
                                    break
                            
                            if considerC3andC4 and mutuallyCompatible:
#                                Check 3-way compatibility. Iterate over all pairs in current structure
                                for l,m in itertools.combinations(currStructure[:lenCurrStructure],2):
                                    if not C3[l,m,k]:
                                        mutuallyCompatible = False;
                                        break
                                    
#                                Check 4-way compatibility
                                if mutuallyCompatible and lenCurrStructure > 2: 
                                    #don't actually need to specify and lenCurrStructure > 2 since itertools.combinations
#                                       would just return an empty set if lenCurrStructure == 2
#                                    Iterate over all triplets in current structure
                                    for l,m,n in itertools.combinations(currStructure[:lenCurrStructure],3):
                                        if not C4[l,m,n,k]:
                                            mutuallyCompatible = False;
                                            break
                            
                            if mutuallyCompatible:
                                lenCurrStructure += 1
                                currStructure.append(k)
                                
                        #Add the structure to our list
                        hasEnoughBondedStrands = True
                        if onlyConsiderBondedStrands:
                            if not any([linkedStems[l] for l in currStructure[:lenCurrStructure]]):
                                hasEnoughBondedStrands = False
                                
                        if minNumStemsInStructure <= lenCurrStructure and hasEnoughBondedStrands:
                            structures.append(copy.copy(currStructure)) #even better, make it numpy array.                        
                            numStructures += 1

                        if numStructures % 5e5 == 0 and printProgressUpdate and prevNumStructures != numStructures:
                            myPrintFxn('Setting up structures: i = ' + str(i) + '; numStructures = ' + str(numStructures), self.fileTxt)
                            prevNumStructures = numStructures #so we don't print the same thing multiple times
                            
                        k = currStructure[lenCurrStructure - 1]
                        currStructure.pop()
                        lenCurrStructure -= 1
                        
        if minNumStemsInStructure <= 1:
            for i in range(numStems): #each stem can be its own helipoint
                structures.append([i])
            numStructures += numStems
            
        if minNumStemsInStructure == 0: #then a completely free structure is allowed as well
            structures.append([])
            numStructures += 1
        
        #for each structure, determine if it has a bond between the two strands or not
        linkedStructures = []
        if numSequences > 1:
            linkedStructures = np.zeros(numStructures,dtype=bool)
            for i in range(numStructures):
                if any([linkedStems[j] for j in structures[i]]):
                    linkedStructures[i] = 1

        
        return(numStructures, structures, linkedStructures)

    def calculateHelipointFreeEnergies(self):
# =============================================================================
#         For each helipoint, calculate its energy and entropy
# =============================================================================
        pass

    def makeHelipointCompatibilityMatrix(self):
# =============================================================================
#         #A helipoint is incompatible with another that contains any stems incompatible with any of the
#         #stems in the first helipoint, or any of the same stems.
#         #A helipoint is also incompatible with one that contains a stem that the first helipoint
#         #could have contained but didn't. For example, if stem i is separated by a bulge loop from stem j
#         #then stems i and j could be contained in a single helipoint. But there would also be a helipoint
#         #that only includes stem i and doesn't include stem j, and that helipoint is incompatible with 
#         #any that include stem j.
# =============================================================================
        numHelipoints = self.numHelipoints
        helipoints = self.helipoints
        adjStems = self.adjStems
        C = self.C
        
        def helipointPairs(helipointsI, helipointsJ):
            """Produce pairs of indexes"""
            for i in helipointsI:
                for j in helipointsJ:
                    yield i, j
        
        CHeli = np.zeros((numHelipoints,numHelipoints),dtype = bool)
        for i in range(numHelipoints):
            for j in range(i,numHelipoints):
                if i == j:
                    CHeli[i,j] = 1
                else:# all([C[x,y] and not adjStems[x,y] for x in helipoints[i] for y in helipoints[j]]):
                    possiblyCompatible = True
                    for x,y in helipointPairs(helipoints[i],helipoints[j]):
                        if adjStems[x,y] or not C[x,y]:
                            possiblyCompatible = False
                            break
                    if possiblyCompatible:       
                        #C[x,y] checks that all stems are mutually compatible. 
                        #not adjStems[x,y] checks that none of the stems are the same or are adjacent
                        #if any of the stems of helipoint i are adjacent to any stem in helipoint j, the two
                        #helipoints can't be compatible. Since in that case, there is another helipoint that includes
                        #both of the adjacent stems.
                        CHeli[i,j] = 1
                        CHeli[j,i] = 1
        return(CHeli)        
    
    def makeHelipoints(self):
# =============================================================================
#        Make the list of helipoints. Here we use the term helipoints slightly differently from TT2NE's usage.
#        They were concerned with coarse-graining over very similar structures, that differ by the presence of 
#        a 1x1 internal loop or the placement of a size-1 bulge loop. However, we're interested not in 
#        coarse-graining over structures, but in doing all of our bond enthalpy/entropy calculations
#        up front at the lowest possible level, to avoid repeating work. We can't do it completely at the level
#        of stems, since at that level it's unclear whether to include terminal mismatches, or dangling ends, etc.
#        We therefore construct our helipoints, combinations of stems that are separated by bulge loops of any size
#        (in that case, the ends of stems include dangling mismatch terms for bulge loops of size > 1 and
#        terminal mismatches for bulge loops of size 1) or 1 x anything internal loops, or in general,
#        stems separated by zero or one intervening unpaired nts. We can then perform the bond energy calculation
#        on these (including terminal mismatch terms at the ends of helipoints).
# =============================================================================
        numStems = self.numStems
        STableBPs = self.STableBPs
        C = self.C
        printProgressUpdate = self.printProgressUpdate
        
# =============================================================================
#        For each stem, make a list of the stems from which it is separated by at most one unpaired nt
#        Perhaps at this step we should also calculate the bond free energy from including both of those stems together.
#        Do this through creating a matrix adjStems where adjStems[i,j] = 1 iff the two stems
#        are separated by at most one unpaired nt. Set diagonal to one by convention. This property is used to make compatibility matrix
# =============================================================================
        adjStems = np.zeros((numStems,numStems),dtype = bool)
        for i in range(numStems):
            edgesOfI = [STableBPs[i][0][0],STableBPs[i][0][1],STableBPs[i][-1][0],STableBPs[i][-1][1]]
            adjStems[i,i] = 1
            for j in range(i+1,numStems):
                edgesOfJ = [STableBPs[j][0][0],STableBPs[j][0][1],STableBPs[j][-1][0],STableBPs[j][-1][1]]
                if C[i,j]:
                    separationBetweenStems = [abs(x - y) for x in edgesOfI for y in edgesOfJ]
                    #We're looking for there to be at least one 1 (meaning stem i is 
                    #directly adjacent to j) or 2 (meaning there is an intervening
                    #unpaired nt.) There can't be any 0's since then they wouldn't be mutually compatible.
                    if any([x < 3 for x in separationBetweenStems]):
                        adjStems[i,j] = 1
                        adjStems[j,i] = 1
               
# =============================================================================
#         Now construct the list of helipoints
# =============================================================================
        numHelipoints = 0
        helipoints = []
        
        for i in range(numStems): #each stem can be its own helipoint
            helipoints.append([i])
        numHelipoints = numStems
        
        prevNumHelipoints = -1 #keep track of this just so that we don't print the same statement multiple times
 
        for i in range(numStems):
            for j in range(i+1,numStems):
                if C[i,j] and adjStems[i,j]:
                    currHelipoint = [i,j]                    
                    lenCurrHelipoint = 2
                    k = j #the next stem we'll try adding (we're about to add one so that's why it's j and not j+1)
                    while lenCurrHelipoint >= 2:
                        while k < numStems - 1:
                            k += 1
                            
                            mutuallyCompatible = True
                            #Check mutual 2-way compatibility between the stem we want to add and all stems 
                            #present in the current helipoint.
                            
                            for l in currHelipoint[:lenCurrHelipoint]: #the same as just "in currHelipoint" but written
                                #this way so that if we want to preallocate space it'll be easier
                                if not C[k,l]:
                                    mutuallyCompatible = False
                                    break
                            
                            if mutuallyCompatible:
                                #Is stem k separated from any of the current stems by at most one unpaired nt?
                                kAdjacentToStemInCurrHelipoint = [True for x in currHelipoint if adjStems[x,k]]#k in listOfAdjacentStems[x]]
                                if not kAdjacentToStemInCurrHelipoint:
                                    mutuallyCompatible = False
                            
                            #Include 3-way and 4-way compatibilities here as well.
                            
                            if mutuallyCompatible:
                                lenCurrHelipoint += 1
                                currHelipoint.append(k)
                                
                        #Add the helipoint to our list
                        helipoints.append(copy.copy(currHelipoint)) #even better, make it numpy array.
                        
                        numHelipoints += 1

                        if numHelipoints % 5e5 == 0 and printProgressUpdate and prevNumHelipoints != numHelipoints:
                            myPrintFxn('Setting up helipoints: i = ' + str(i) + '; numHelipoints = ' + str(numHelipoints), self.fileTxt)
                            prevNumHelipoints = numHelipoints #so we don't print the same thing multiple times
                            
                        k = currHelipoint[lenCurrHelipoint - 1]
                        currHelipoint.pop()
                        lenCurrHelipoint -= 1
        
        return(numHelipoints, helipoints, adjStems)
               
    def calculateStemFreeEnergiesPairwise(self):
# =============================================================================
#         Calculate the bond energy and entropy of each stem, not including what to do at the 
#        ends of the stems (i.e. dangling ends, terminal mismatches, and the like)
# =============================================================================
        numStems = self.numStems
        STableStructure = self.STableStructure
        sequenceInNumbers = self.sequenceInNumbers
        unmatchedBPPenalty = self.unmatchedBPPenalty
        
#        #Define energy (enthalpy; deltaH; units of kcal/mol) and entropy (deltaS; units of kcal/(mol*K)) matrices for bonds.
#        energyMatrices, entropyMatrices = bondFreeEnergies()
#        
#        stemEnergies = np.zeros(numStems)
#        stemEntropies = np.zeros(numStems)
#        
#        if self.includeTerminalAUATPenalties:
#            (terminal_AU_penalty_energy, terminal_AU_penalty_entropy, 
#             terminal_AT_penalty_energy, terminal_AT_penalty_entropy) = terminalATAUPenalties()
#        else:
#            (terminal_AU_penalty_energy, terminal_AU_penalty_entropy, 
#             terminal_AT_penalty_energy, terminal_AT_penalty_entropy) = [0,0,0,0]
        
        
        RNARNACount = scipy.sparse.lil_matrix((6,16),dtype=int) #np.zeros((6,4,4),dtype = int)
        #First index tells you if the first bp of the set is AU (0) CG (1) GC (2) UA (3) GU (4) or UG (5)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3) (which row of table 1 in Serra & Turner).
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3) (which column of table 1 in Serra & Turner).
            
#    Had to make this and the others 2D array to be able to use scipy sparse functions, so 
#    second/third index is replaced by (4*second index) + third index. That's why it's
#        np.zeros(6,16) and not np.zeros(6,4,4).
#        lil_matrix was chosen because, from scipy reference page, it is fast for constructing
#        sparse matrices incrementally. For operations (later we'll multiply it) it should be converted to
#        another form.

        DNADNACount = scipy.sparse.lil_matrix((4,16),dtype=int) #np.zeros((4,16),dtype = int)
        #First index tells you if the first bp of the set is AT (0) CG (1) GC (2) or TA (3)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or T(3)
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or T(3)
            
        RNADNACount = scipy.sparse.lil_matrix((8,16),dtype=int) #np.zeros((8,16),dtype = int)
        #First index tells you if the set is (putting RNA first) 
        #5'AX3'/3'TY5' (0) 5CX3/3GY5 (1) 5GX3/3CY5 (2) 5UX3/3AY5 (3) 5XA3/3YT5 (4) 5XC3/3YG5 (5) 5XG3/3YC5 (6) or 5XU3/3YA5 (7).
        #Second index tells you if X is A (0) C (1) G (2) or U (3)
        #Third index tells you if Y is A (0) C (1) G (2) or T (3)
        
        terminalAUATCount = scipy.sparse.lil_matrix((1,2),dtype=int)
#        Number of terminal AU (or GU) pairs, number of terminal AT pairs
        
        unknownCount = scipy.sparse.lil_matrix((1,1),dtype=int) #np.zeros((1,1), dtype = int)
        
        bondFECounts = [[RNARNACount, DNADNACount, RNADNACount, terminalAUATCount, unknownCount]]
        
        stemFECounts = bondFECounts*numStems
        
        
        for stemNumber, stem in enumerate(STableStructure):
            numBonds = int(len(stem)/2)
            for j in range(numBonds - 1):
                firstNtIndex = stem[j] #the 5' nt (its position in the sequence)
                firstBPIndex = stem[j+numBonds] #what it's bonded to
                secondNtIndex = firstNtIndex + 1 #the 3' nt
                secondBPIndex = firstBPIndex - 1 #what it's bonded to.
                #we're here assuming antiparallel stem. Otherwise, change to + 2*isParallel - 1
                
                index, bpType = freeEnergyMatrixIndices(sequenceInNumbers,firstNtIndex,firstBPIndex,
                                                        secondNtIndex,secondBPIndex, bound = [1,1], 
                                                        unmatchedBPPenalty = unmatchedBPPenalty)
                
                stemFECounts[stemNumber][bpType][index[0], index[1]] += 1
                
                
#                stemEnergies[stemNumber] += energyMatrices[bpType][index[0], index[1]]
#                stemEntropies[stemNumber] += entropyMatrices[bpType][index[0], index[1]]
                
                
# =============================================================================
#    Include terminal penalties for AU, GU, or AT pairs at the ends of helices 
#    (from Xia et al. 1998 for RNA, or SantaLucia and Hicks (2004) for DNA):
#    GU penalty is assumed to be same as for AU (Xia paper, or NNDB)
#    From Xia, Mathews, Turner paper (from Soll, Nishimura, Moore book):
#    "Note that when dangling ends or terminal mismatches follow terminal AU 
#    or GU pairs, the penalty for a terminal AU is still applied"
#                
#    Guess that  this penalty also applies for DNA/RNA interactions since
#    it's physically motivated by the number of H-bonds AU/AT pairs have
#     if (bpType == 0 and index[0] == 0) or (bpType == 2 and (index[0] == 3 or index[0] == 7)):
#         terminalAUPenaltyCounter += 1
#     elif (bpType == 1 and index[0] == 0) or (bpType == 2 and (index[0] == 0 or index[0] == 4)):
#         terminalATPenaltyCounter += 1
# =============================================================================
            for j in [0, numBonds - 1]: #penalties apply to ends of helices
                if sequenceInNumbers[stem[j]] == 4 or sequenceInNumbers[stem[j + numBonds]] == 4:
                    stemFECounts[stemNumber][3][0,0] += 1 #terminal AU/GU was found
                elif sequenceInNumbers[stem[j]] == 8 or sequenceInNumbers[stem[j + numBonds]] == 8:
                    stemFECounts[stemNumber][3][0,1] += 1 #terminal AT was found
            
#            if not self.considerAllAsTerminalMismatches:
#                for j in [0, numBonds - 1]: #penalties apply to ends of helices
#                    if sequenceInNumbers[stem[j]] == 4 or sequenceInNumbers[stem[j + numBonds]] == 4:
#                        stemEnergies[stemNumber] += terminal_AU_penalty_energy
#                        stemEntropies[stemNumber] += terminal_AU_penalty_entropy
#                    elif sequenceInNumbers[stem[j]] == 8 or sequenceInNumbers[stem[j + numBonds]] == 8:
#                        stemEnergies[stemNumber] += terminal_AT_penalty_energy
#                        stemEntropies[stemNumber] += terminal_AT_penalty_entropy
        
        return(stemFECounts)
        
        
        
#    def calculateStemFreeEnergiesSingle(self):
## =============================================================================
##         Calculate the bond energy of each stem assuming each base pair gives an energy independent of its context.
## =============================================================================       
#        numStems = self.numStems
#        STableStructure = self.STableStructure
#        sequenceInNumbers = self.sequenceInNumbers
#        
#        bondEnergyMatrix = np.array([[-2.06209, -0.182809, -2.06209, -2.06209],
#                                     [-0.182809, 0, -4.86521, -2.06209],
#                                     [-2.06209, -4.86521, 0, -2.06209],
#                                     [-2.06209, -2.06209, -2.06209, 0]])
#        #    BondEnergyMatrix[i,j] gives energy for bond between nucelotides i and j
#        #    where A=0, C=1, G=2, U=3. This comes from Cragnolini's paper
#        #    "Coarse-grained HiRE-RNA model for ab initio RNA folding beyond simple 
#        #    molecules, including noncanonical and multiple base pairings"
#        #   There is no functionality currently for DNA/DNA or RNA/DNA pairs, only for RNA/RNA
#        
#        stemEnergies = np.zeros(numStems)
#        stemEntropies = np.zeros(numStems)
#        
#        for stemNumber, stem in enumerate(STableStructure):
#            numBonds = int(len(stem)/2)
#            for j in range(numBonds):
#                firstNtIndex = stem[j] #the nt (its position in the sequence)
#                firstBPIndex = stem[j+numBonds] #what it's bonded to
#                
#                firstNt = sequenceInNumbers[firstNtIndex] #the identity of the nt (A is 1, C is 2, etc.)
#                firstBP = sequenceInNumbers[firstBPIndex]
#                
#                if firstNt and firstBP: #if either of the nts is unknown, don't add anything to the free energy
#                    #this part isn't actually necessary, since if they're unknown, they couldn't have been bound.
#                
#                    if firstNt > 4: #There is no functionality currently for DNA/DNA or RNA/DNA pairs, only for RNA/RNA
#                        firstNt -= 4
#                    if firstBP > 4:
#                        firstBP -= 4
#                    
#                    stemEnergies[stemNumber] += bondEnergyMatrix[firstNt - 1, firstBP - 1]
#                    #need -1 for RNA nts since nts are numbered 1-4 but indices are 0-3
#                
#        return(stemEnergies,stemEntropies)
                    
    
    def removeDisallowedStems(self):
# =============================================================================
#         If there are frozenBPs, there may be some stems that are incompatible with
#        stems that need to be included to include the frozenBPs. We therefore want to remove these
#        stems. 
#        This means modifying the C matrix, the frozenStems, STableBPs, and STableStructure lists, and numStems.
# =============================================================================
        numStems = self.numStems
        C = self.C
        STableBPs = self.STableBPs
        STableStructure = self.STableStructure
        frozenStems = self.frozenStems
        linkedStems = self.linkedStems
        
        disallowedStems = [i for i in range(numStems) if not C[i,i]]
        #There should be exactly no overlap between disallowedStems and frozenStems. Otherwise, something is weird.
#        I guess I could imagine such a scenario (a stem required for freezing one BP incompatible with another) but I can't
#        come up with an example off the top of my head
        if len([i for i in disallowedStems for j in frozenStems if i in j]):
            myPrintFxn('Something is strange -- a disallowed stem is also a frozen stem', self.fileTxt)
        
        #correct the various arrays/lists   
        correctedNumStems = numStems - len(disallowedStems)
        
        semiCorrectedC = np.delete(C,disallowedStems,0) #delete the rows of C corresponding to disallowed stems
        correctedC = np.delete(semiCorrectedC,disallowedStems,1) #and delete the columns
        
        correctedSTableBPs = [STableBPs[i] for i in range(numStems) if i not in disallowedStems] #mind that this is the old numStems
        correctedSTableStructure = [STableStructure[i] for i in range(numStems) if i not in disallowedStems]
        
        correctedLinkedStems = np.delete(linkedStems,disallowedStems) #delete the stems that are entirely disallowed
        
        #modifying frozenStems is a bit harder because it doesn't involve just deleting stems, but rather changing their numbering.
        #If stem 4 used to be frozen, but now stem 3 has been deleted, the frozen stem needs to be renamed 3.
        correctedFrozenStems = copy.copy(frozenStems)
        for i in range(len(correctedFrozenStems)):
            j = 0
            while j < len(correctedFrozenStems[i]):
                if correctedFrozenStems[i][j] in disallowedStems: #don't expect this to ever happen
                    del correctedFrozenStems[i][j] 
                else: 
                    correctedFrozenStems[i][j] -= len([True for x in disallowedStems if x < correctedFrozenStems[i][j]])
                    j += 1
        
        correctedFrozenStems = [i for i in correctedFrozenStems if i] #remove empty lists                     
        
        return(correctedNumStems,correctedFrozenStems,correctedSTableBPs,correctedSTableStructure,correctedC,correctedLinkedStems)     
    

    def makeCompatibilityMatrix(self):
# =============================================================================
#         Defines a matrix C of numStems x numStems where the ij'th element is 1
#        if stems i and j can coexist in the same structure and 0 otherwise.
#        We set the diagonal elements to 1.
# =============================================================================
        numStems = self.numStems
        STableStructure = self.STableStructure
        minNtsInHairpin = self.minNtsInHairpin
        frozenStems = self.frozenStems
        STableBPs = self.STableBPs
        linkedStems = self.linkedStems
        numSequences = self.numSequences
        allowPseudoknots = self.allowPseudoknots
        
# =============================================================================
#         First, start with the basic compatibility matrix. Stems i and j are compatible
#        if they don't share any of the same nts, they don't make hairpins less than minNtsInHairpin
#        and, if we disallow pseudoknots, they don't make a pseudoknot.
# =============================================================================
        
        C = np.zeros((numStems,numStems),dtype = bool)
        for i in range(numStems):
            for j in range(i,numStems):
                if i == j:
                    C[i,j] = 1
                    C[j,i] = 1
                elif not twoListsShareElement(STableStructure[i],STableStructure[j]):
                    C[i,j] = 1
                    C[j,i] = 1
                    disallowPseudoknotsIJ = not allowPseudoknots
                    if numSequences > 1:
                        if linkedStems[i] or linkedStems[j]:
                            disallowPseudoknotsIJ = False
                    if disallowPseudoknotsIJ:
# =============================================================================
#                 next part lets us disallow pseudoknots
#                 From Mauri 2005 (substitute i,j,k,l for a,b,c,d):
#                 If ntd a is paired with ntd b>a, and c>a is paired with d>c,
#                 then either a<b<c<d or a<c<d<b, but a<c<b<d is a pseudoknot.
#                   This only applies if we're considering a single RNA moleucle, not two bonded molecules.
# =============================================================================
                        a = STableStructure[i][0]
                        b = STableStructure[i][int(len(STableStructure[i]/2))]
                        c = STableStructure[j][0]
                        d = STableStructure[j][int(len(STableStructure[j]/2))]
                        
                        if c < a: #switch around labels so we have c>a
                            a = STableStructure[j][0]
                            b = STableStructure[j][int(len(STableStructure[j]/2))]
                            c = STableStructure[i][0]
                            d = STableStructure[i][int(len(STableStructure[i]/2))]
                        
                        if (a<c and c<b and b<d):
                            C[i,j] = 0
                            C[j,i] = 0
                    
                    #make sure each hairpin has at least minBPInHairpin unpaired bps in it. 
                    #Right now, this only constrains pairwise compatible regions, but it's a start
                    iHairpin = list(range(STableStructure[i][int(len(STableStructure[i])/2) - 1] + 1, 
                                         STableStructure[i][int(len(STableStructure[i])/2)]))
                    #The unpaired nts between the start and end of stem i
                    jHairpin = list(range(STableStructure[j][int(len(STableStructure[j])/2) - 1] + 1, 
                                         STableStructure[j][int(len(STableStructure[j])/2)]))
                    #same for stem j
                    
                    #if the number of unpaired nts is less than minNtsInHairpin
                    if (len(np.setdiff1d(iHairpin,STableStructure[j])) < minNtsInHairpin or 
                        len(np.setdiff1d(jHairpin,STableStructure[i])) < minNtsInHairpin):
                        C[i,j] = 0
                        C[j,i] = 0
        
        #That does it for the basic compatibility matrix.
        
# =============================================================================
#         #Next, if a stem isn't compatible with a stem we need to include, it can't be
#         #present in any structure so don't include it. We designate these stems by setting their diagonal elements to zero.
# =============================================================================
        if frozenStems:
            for i in range(numStems): #for each stem
                if C[i,i]: #not important here, but useful when we repeat this process after the next block.
                    for j in range(len(frozenStems)): #for each set of stems that need to be included
                        compatibleWithFrozen = False #is stem i compatible with at least one of the stems 
                        #out of the j'th list in frozenStems? It needs to be (for all j) to be included in any structure.
                        for k in range(len(frozenStems[j])):
                            if C[i,frozenStems[j][k]]:
                                compatibleWithFrozen = True #it's compatible with at least one of the stems
                                break
    
                        if not compatibleWithFrozen: #if for any set of regions one of which needs to be included,
                            #region i isn't compatible with any of them, we can't include it in our structure
                            C[i,:] = 0
                            C[:,i] = 0

# =============================================================================
#        We don't want to allow two substems to be compatible if including
#        both of them just looks like one longer stem. For example, a substem
#        of ntd 1 bonded to 10 shouldn't be compatible with a substem of ntd 2
#        bonded to 9, even if they came from different full stems initially. 
#        Really, we just need to check if the combined stems form one single stem.
# =============================================================================
        for i in range(numStems):
            for j in range(i+1,numStems):
                if C[i,j]: #not worth going through this code if we already know the regions aren't compatible
                    combinedBPs = STableBPs[i] + STableBPs[j]
                    combinedStem = bpsList2structure(combinedBPs)
                    if len(combinedStem) == 1: #then the combined stems form one single stem
                        C[i,j] = 0
                        C[j,i] = 0
        
# =============================================================================
#         Finally, we repeat the process of removing stems incompatible with frozenStems.
#           We could of course just do it here and delete the one above, but it's faster to do it
#           earlier as well, so that the previous step is fast if many stems are disallowed then.
#           This is a precise copy of what we did two blocks ago.
# =============================================================================
        if frozenStems:
            for i in range(numStems): #for each stem
                if C[i,i]:
                    for j in range(len(frozenStems)): #for each set of stems that need to be included
                        compatibleWithFrozen = False #is stem i compatible with at least one of the stems 
                        #out of the j'th list in frozenStems? It needs to be (for all j) to be included in any structure.
                        for k in range(len(frozenStems[j])):
                            if C[i,frozenStems[j][k]]:
                                compatibleWithFrozen = True #it's compatible with at least one of the stems
                                break
    
                        if not compatibleWithFrozen: #if for any set of regions one of which needs to be included,
                            #region i isn't compatible with any of them, we can't include it in our structure
                            C[i,:] = 0
                            C[:,i] = 0
            
        return(C)
    

    def frozenStemsFromFrozenBPs(self): #still need to check more carefully
# =============================================================================
#        given a list of base pairs that we want fixed, what regions do you need to have to include all of them?
#        frozenBPs is a list of base pairs we want to be in the final structures
#        The question is, in order to constrain e.g. base pair number 1, what is the list of stems
#        that contain base pair number 1? You need to have one of those stems in your structures
#        Make a mega-list of lists of stems, such that each structure needs to have 
#        one stem from each list in the mega-list.
#        
#       You might think you could just do:
#       frozenStems = [[i for i in range(numStems) if j in STableBPs[i]] for j in frozenBPs]
        #This gives, for each bp in frozenBPs, the list of stems that contain it.
        #However, this is not the minimal list. For example, one of the elements of the first list
        #could be incompatible with all of the elements of the second list; therefore, 
        #that element actually can't be chosen at all and should be excluded.
        #An example is if two bps are adjacent.
# =============================================================================
        
        frozenBPs = self.frozenBPs
        if not frozenBPs:
            return([])
        
        numStems = self.numStems
        STableBPs = self.STableBPs
        #if two or more base pairs of frozenBPs are adjacent, then any stem that includes them
        #must include all of them to be considered (otherwise there'd be no way to include the others).
        frozenStructBPs = bpsList2structBPs(frozenBPs) #a list of frozen stems where each stem is in base-pair (bp) format
        
        frozenStems = [[i for i in range(numStems) if all([k in STableBPs[i] for k in j])] for j in frozenStructBPs]
        
        for j in range(len(frozenStems)):
            if not frozenStems[j]: #if we couldn't find any regions in STable that'll satisfy the j'th stem in frozenStructBPs
                myPrintFxn('In frozenStemsFromFrozenBPs, couldn''t find a match in STable for (and removed) the stem ' 
                      + str(frozenStructBPs[j]), self.fileTxt)
        
        frozenStems = [i for i in frozenStems if i] #remove empty lists 
        return(frozenStems)
        
    
    def findLinkedStems(self):
# =============================================================================
#         Given the STable, find which stems correspond to the two strands being bound
# =============================================================================
        STableStructure = self.STableStructure
        numStems = self.numStems
        numNt1 = self.numNt1
        
        linkedStems = np.zeros(numStems,dtype=bool)
        if numStems > 1: #I'm pretty sure the next part is fine even if numStems == 1
            linkedStemsList = [i for i in range(numStems) if (any([j < numNt1 for j in STableStructure[i][:int(len(STableStructure[i])/2)]])
                           and any([j > numNt1 for j in STableStructure[i][int(len(STableStructure[i])/2):]]))]
            #The stem has to have the first strand before the second since j>i in the main loop of createSTable
            for i in range(numStems):
                if i in linkedStemsList:
                    linkedStems[i] = 1            
        
        return(linkedStems)
    
    
    def createSTable(self):
# =============================================================================
#         Make the STable. It is a made of two lists.
#         The first list is the stem written in "structure" format.
#         This is a vector of length 2*(length of stem) which defines the base pairs
#         in the stem. For example, [1 2 3 51 50 49] means ntd 1 is bonded to 51, 2
#         to 50, 3 to 49. The second list is the same data in another
#         format: a (length of stem)x2 matrix, where each row defines a base pair.
#         The same stem would be given by [[1, 51],[2, 50],[3 49]].
# =============================================================================
        numNt = self.numNt
        numNt1 = self.numNt1
        seqInNum = self.sequenceInNumbers
        minBPInStem = self.minBPInStem
        minNtsInHairpin = self.minNtsInHairpin
        substems = self.substems
        
        maxNumStems = self.maxSizeSTable #to preallocate space for STable.
        STableStructure = [None]*maxNumStems #first column of STable
        STableBPs = [None]*maxNumStems #second column of STable
        
        B = np.zeros((numNt,numNt)) #matrix describing which bases are allowed to bond to which others
        for i in range(numNt):
            for j in range(numNt):
                if isComplementary(seqInNum[i], seqInNum[j]):
                    B[i,j] = 1
                    B[j,i] = 1
#                if seqInNum[i]==1 and seqInNum[j]==4: #for Watson-Crick base pairs, Bij=1
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 5 and seqInNum[j] == 8:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 1 and seqInNum[j] == 8:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 4 and seqInNum[j] == 5:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 2 and seqInNum[j] == 3:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 6 and seqInNum[j] == 7:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 2 and seqInNum[j] == 7:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 3 and seqInNum[j] == 6:
#                    B[i,j] = 1
#                    B[j,i] = 1
#                elif seqInNum[i] == 3 and seqInNum[j] == 4:
#                    #can set Bij=2 if we want to distinguish non-Watson-Crick base
#                    #pairs, GU. This could allow us, for example, to ensure GU base
#                    #pairs don't close any stem.
#                    B[i,j] = 1
#                    B[j,i] = 1
#                #Can't have DNA G bind to RNA U
        
# =============================================================================
#        a stem is started when we find two BPs, i and j, that can bind. Then we
#        look if i+1 can bind to j-1, etc., and in the next for loop, if i+1 can
#        bind to j+1, etc.
#        For now, we don't include substems -- i.e. we only consider the longest possible stem
#        for each starting bp.
# =============================================================================
 
        numStems = 0
        
# =============================================================================
#       Allow for parallel stems as well if you'd like here
#        if self.allowParallelStrands:
#            isParallelVec = [0,1]
#        else:
#            isParallelVec = [0]
# =============================================================================
        
        BPsAlreadyConsidered = [] #so we avoid subsets of stems until we consider them explicitly later
        maxI = numNt - 2*minBPInStem - minNtsInHairpin + 1 #we need at least minBPInStem consecutive base pairs binding
#        if onlyConsiderBondedStrands: %Should consider other stems regardless, since even if onlyConsiderBondedStrands,
#        we want to include structures that include both inter- and intra-molecular bonds
#            maxI = numNt1 #don't go beyond the end of the first sequence for i if we only consider bonded strands
        for i in range(maxI): 
            minJ = i + 2*minBPInStem + minNtsInHairpin - 2
#            if onlyConsiderBondedStrands:
#                minJ = numNt1 + len(self.linkerPos) - 1 #don't go at all into the first sequence for j 
##                   if we only consider bonded strands. The -1 isn't necessary -it's just to be confident we're not 
##                   making an off-by-one error (and since the B[i,j] == 0 if j is in the linker, it isn't an issue).
            for j in range(numNt - 1, minJ, -1): #range(start,stop,step)
                #so for example, if we have numNt = 7, minBPInStem = 2, minNtsInHairpin = 3, the only possibility
                #is a hairpin with 0 bonded to 6, and 1 to 5
                if B[i,j] == 1: #if they can bond
                    
                    if [i,j] not in BPsAlreadyConsidered:
#                        BPsAlreadyConsidered.append([i,j]) #unnecessary since we definitely won't consider this bp again
                         
                        currentStemI = [i]
                        currentStemJ = [j]
                        listOfPairs = [[i,j]]
                        
# =============================================================================
#                         If we do include parallel stems, add stems of length 1 in after
#                        adding in all the substems (and change the substem addition to not include stems of length 1)
#                       Code:
#                        if minBPInStem == 1: #Then even just ntd i bonded to j is a full stem
#                            STableStructure[numStems] = currentStemI + currentStemJ
#                            STableBPs[numStems] = [currentStemI + currentStemJ]
#                            numStems += 1
# =============================================================================
                        
                        #now try to lengthen the stem, to include nts i + lenStem and j - lenStem
                        lenStem = 0 #lenStem is one less than the number of bps in a stem.
                        endOfStem = False
                        while not endOfStem:
                            lenStem += 1
                            newI = i + lenStem
                            newJ = j - lenStem
                            if (newI > numNt - 1 or newJ < 0 or #if we've gone beyond the edges of the RNA
                                newJ - newI <= minNtsInHairpin or #or the bps are too close together 
                                B[newI, newJ] == 0): #or the bps can't actually bind
                                
                                endOfStem = True
                                lenStem -= 1 #to correct for adding one earlier
                            else:
                                currentStemI.append(newI)
                                currentStemJ.append(newJ)
                                listOfPairs.append([newI,newJ])
                                BPsAlreadyConsidered.append([newI, newJ])
                        #now that we've finished making the longest possible stem starting with the base pair i,j
                        #add that stem to STable
                        if len(currentStemI) >= minBPInStem:
                            STableStructure[numStems] = currentStemI + currentStemJ
                            STableBPs[numStems] = listOfPairs
                            numStems += 1
                        
        if self.onlyAllowSubsetsOfLongestStems:
            #remove all stems but the longest stems so that in the next step we add
            #only the subsets of the longest stems
            
            maxLengthStem = minBPInStem
            for i in range(numStems):
                if len(STableStructure[i])/2 > maxLengthStem:
                    maxLengthStem = len(STableStructure[i])/2
                    
            minMaxLengthStem = 16 #don't chop off stems of length more than minMaxLengthStem
            maxLengthStem = min(maxLengthStem,minMaxLengthStem); 
            STableBPs = [STableBPs[i] for i in range(numStems) if len(STableStructure[i])/2 >= maxLengthStem]
            STableStructure = [STableStructure[i] for i in range(numStems) if len(STableStructure[i])/2 >= maxLengthStem]
            numStems = len(STableStructure)
            
            #since we preallocate space, need to make more empty room since we just removed all extra space
            STableBPs += [None]*maxNumStems
            STableStructure += [None]*maxNumStems
            
# =============================================================================
# Based on a comment made by Zuker & Sankoff (1984) regarding the
# Pipas-McMahon algorithm, we also want to include substems -- i.e.
# fractions of the whole stem. Thus, we allow the possibility that two
# full stems may be incompatible, but that their compatible substems are
# more desirable than a single full stem. 
# This process of course takes zero time if substems == 0
# =============================================================================
        for i in range(numStems):
            fullStemI = STableStructure[i] #the full stem we're considering
            lenStem = int(len(fullStemI)/2)
            
            #What substems should we consider? This is given by the substems argument to the code
            if substems == 'all':
                minBPInStemSub = minBPInStem
            else: #then substems is an integer
                minBPInStemSub = max(lenStem - substems,minBPInStem)
            
            #we can make substems of length lengthStem-1 till minBPInStem. 
            for j in range(1,lenStem-minBPInStemSub+1): #There are j possible lengths.
                #j also tells us how much shorter the substem is compared to the full stem.
                possSubstemCounters = np.arange(j+1)
                if self.onlyConsiderSubstemsFromEdges:
                    possSubstemCounters = [0,j]
                
                for k in possSubstemCounters: #substems come from getting rid of either edge.
                    truncatedStemI = fullStemI[k:lenStem-j+k] + fullStemI[lenStem+k:len(fullStemI)-j+k]
                    #truncatedStemI is the truncated stem. Add it to STable
                    currentStemI = truncatedStemI[:int(len(truncatedStemI)/2)]
                    currentStemJ = truncatedStemI[int(len(truncatedStemI)/2):]
                    STableStructure[numStems] = currentStemI + currentStemJ
                    
                    listOfPairs = [[currentStemI[0],currentStemJ[0]]]
                    for l in range(1,len(currentStemI)):
                        listOfPairs.append([currentStemI[l],currentStemJ[l]])
                    STableBPs[numStems] = listOfPairs
                    numStems += 1
    
        #remove preallocated space
        STableBPs = [STableBPs[i] for i in range(numStems)]
        STableStructure = [STableStructure[i] for i in range(numStems)]
        
        return(numStems, STableStructure, STableBPs)
    
    
    def multipleStrandsSetup(self):
    
        sequences = self.sequences
        DNA = self.DNA
        minNtsInHairpin = self.minNtsInHairpin
        
        numSequences = len(sequences);
    
        sequence1 = sequences[0]
        numNt1 = len(sequence1)
        sequence2 = ''
        
        if numSequences == 1:
            linkerPos = np.array([],dtype = int)
            sequence = sequence1
            numNt2 = 0 
            DNA = DNA[0]
            isDNA = [DNA]*len(sequence) #which ntds are DNA (rather than RNA). 1 for ntds which are DNA, 0 for RNA.
                
        elif numSequences == 2:
            sequence2 = sequences[1]
            numNt2 = len(sequence2)
            
            linker = 'O'*(2*minNtsInHairpin) #want linker to not be allowed to bind, and long enough not to constrain anything
            sequence = sequence1 + linker + sequence2
            linkerPos = np.arange(numNt1, numNt1 + len(linker)) #ntds of sequence which are actually the linker
            
            #which ntds are DNA (rather than RNA). 1 for ntds which are DNA, 0 for RNA.
            isDNA = [DNA[0]]*numNt1
            isDNA += [0]*len(linker)
            isDNA += [DNA[1]]*numNt2
        
        numNt = len(sequence)
        #convert sequence into numbers 
        #write the sequence as a vector where for RNA each instance of A is replaced with 1, 
        #C with 2, G with 3, U with 4, and other with 0.
        #For DNA, ntds become 5, 6, 7, 8 for A, C, G, T.
        
        sequenceInNumbers = np.array([1 if (sequence[i] == 'A' and not isDNA[i]) else 
                                      2 if (sequence[i] == 'C' and not isDNA[i]) else
                                      3 if (sequence[i] == 'G' and not isDNA[i]) else 
                                      4 if (sequence[i] == 'U' and not isDNA[i]) else
                                      5 if (sequence[i] == 'A' and isDNA[i]) else 
                                      6 if (sequence[i] == 'C' and isDNA[i]) else
                                      7 if (sequence[i] == 'G' and isDNA[i]) else 
                                      8 if (sequence[i] == 'T' and isDNA[i]) else
                                      0 for i in range(numNt)])
    
        return(sequence,sequenceInNumbers,numSequences,numNt,numNt1,numNt2,linkerPos)
        #sequence is a string, sequenceInNumbers is a np array, numSequences is an int
        #numNt is an int, numNt1 is an int, numNt2 is an int, linkerPos is a np array    
    
    
    def returnFxn(self):
        #return whatever needs to be returned
        
        if self.makeFigures:
            self.figureFxn()

        if self.printProgressUpdate:
            myPrintFxn('Total time = ' + str(time.time() - self.startTime), self.fileTxt)
            myPrintFxn('sortedProbs = ' + str(self.sortedProbs[:min(len(self.sortedProbs),5)]), self.fileTxt)
            myPrintFxn('sortedFEs = ' + str(self.sortedFEs[:min(len(self.sortedFEs),5)]), self.fileTxt)
            myPrintFxn('indexSort = ' + str(self.indexSort[:min(len(self.indexSort),5)]), self.fileTxt)
            
        return(self)
        
    def figureFxn(self):
        fileFigs = self.fileFigs
        #MFE structure plot
        if self.toSave:
            if not os.path.exists(fileFigs):
                os.makedirs(fileFigs)
        
        f = drawRNAStructure(self.structureWithNthLowestFE(0), self.sequence, 
                             linkerPos = self.linkerPos, showFig = True)
        if self.toSave:
            f.savefig(fileFigs + 'MFE_RNA_Structure.png')
        
        #Histogram of FEs
        displayedFE = [self.sortedFEs[i] for i in range(self.numStructures) if self.sortedFEs[i] < 30]
        if len(displayedFE) > 20:
            n, bins, patches = plt.hist(displayedFE, bins=min(500,int(len(displayedFE)/10)))
            plt.xlabel('Free Energy (kcal/mol)')
            plt.ylabel('Structures')
            if self.toSave:
                plt.savefig(fileFigs + 'FE_Histogram') # bbox_inches='tight' should get rid of extra whitespace
            plt.show()
            
#        Histogram of structure probabilities
            
#        Histogram of graph probabilities
        

    def loadVars(self,folderName):
# =============================================================================
# #       if we ran the same code previously, we can just load the results
# =============================================================================
        realStartTime = self.startTime
        realMakeFigures = self.makeFigures
        try: 
            loadedVars = load(self.fileVars)
            loadedVars.startTime = realStartTime
            loadedVars.makeFigures = realMakeFigures
            return(loadedVars.returnFxn()) #return variables and make figures
            
        except: 
#          otherwise, we might have ran a similar enough version where we can easily return the results
#          The simplest case is if we stored graphs in the past but don't need to now. 
            if not self.storeGraphs: 
                self.storeGraphs = True
                fileVars, fileTxt, folderName, fileFigs = self.makeFileName()
                try:
                    loadedVars = load(fileVars)
                    self.storeGraphs = False
                    loadedVars.storeGraphs = False
                    loadedVars.returnFxn() #the return function shouldn't re-save the results (waste of time, but also could be misleading here)
                except: #then we need to undo the setting of storeGraphs to True
                    self.storeGraphs = False
                    fileVars, fileTxt, folderName, fileFigs = self.makeFileName()
            
#           Otherwise, we can check if we have the results at a different T, vs, b (gamma),
#               duplexEntropy penalty, pairwiseEnergies, corruptFESeed, 
#               includeTerminalAUATPenalties, and/or includeDanglingEnds.
                    
            trueT = self.T
            trueVs = self.vs
            trueDuplexEntropyPenalty = self.duplexEntropyPenaltyInKB            
            truePairwiseEnergies = self.pairwiseEnergies
            trueStoreGraphs = self.storeGraphs
            trueG = self.g
            trueB = self.b
            trueCorruptFESeed = self.corruptFESeed
            trueIncludeTerminalAUATPenalties = self.includeTerminalAUATPenalties
            trueIncludeDanglingEnds = self.includeDanglingEnds
            trueUnboundButCouldBindPenaltiesE = self.unboundButCouldBindPenalties[0]
            trueUnboundButCouldBindPenaltiesS = self.unboundButCouldBindPenalties[1]
            
            tryFiles = folderName + self.fileNameFxn(self.numSequences,str(self.pairwiseEnergies), 
                str(self.storeGraphs),strT = '*', strVs = '*', strG = '*', 
                strDuplexEntropyPenaltyInKB = '*', strCorruptFESeed = '*', 
                strIncludeTerminalAUATPenalties = '*',strIncludeDanglingEnds = '*',
                strUnboundButCouldBindPenaltiesE = '*', strUnboundButCouldBindPenaltiesS = '*') + 'vars.*'
                
            for file in glob.glob(tryFiles):
                loadedVars = load(file)
                if loadedVars.storeGraphs == trueStoreGraphs or trueStoreGraphs == False:
#                    Otherwise, if you want to get the graphProbs but we didn't store graphs
#                    on the previous run, you won't be able to get the graphProbs.
                    loadedVars.T = trueT
                    loadedVars.vs = trueVs
                    loadedVars.duplexEntropyPenaltyInKB = trueDuplexEntropyPenalty
                    loadedVars.pairwiseEnergies = truePairwiseEnergies
                    loadedVars.storeGraphs = trueStoreGraphs
                    loadedVars.g = trueG
                    loadedVars.b = trueB
                    loadedVars.corruptFESeed = trueCorruptFESeed
                    loadedVars.includeTerminalAUATPenalties = trueIncludeTerminalAUATPenalties
                    loadedVars.includeDanglingEnds = trueIncludeDanglingEnds
                    loadedVars.unboundButCouldBindPenalties = [trueUnboundButCouldBindPenaltiesE,trueUnboundButCouldBindPenaltiesS]
                    
                    loadedVars.startTime = realStartTime
                    loadedVars.makeFigures = realMakeFigures
                    loadedVars.fileVars, loadedVars.fileTxt, loadedVars.folderName, loadedVars.fileFigs = loadedVars.makeFileName()
                    return(loadedVars.postCalculationFxn())
                
        return(False) #if we weren't able to load anything
            
        
    def makeFileName(self):
# =============================================================================
#         #Make the folder and the file names
#         
#         #Example: for sequences == ['ABCDEF','GHIJ']
#         #File is saved in the folder RNALandscape_#version/numSeq_2/numNt_6/A/AB/ABC/ABCD/ABCDE/ABCDEF/numNt_4/G/GH/GHI/GHIJ/
# 
#         #if runningOnCluster:
#         #    folderName = strcat('/n/home05/okimchi/RNALandscape_multipleSequences_frozenBonds_07_07_18','_outputs/');
# =============================================================================
                
        folderName = 'RNALandscape_' + self.version + '/numSeq_' + str(self.numSequences) + '/'
        numNtVec = [self.numNt1,self.numNt2]
        
        maxNumNtsNested = 10 #gets to be too many nested folders if the sequence is too long
        maxNumCharsInFolderName = 100 #Can't/shouldn't have folder name of more than 100 characters:
        
        for k in range(self.numSequences):
            numNtK = numNtVec[k]
            nNtNested = min(numNtK + 1,maxNumNtsNested);
            folderName += 'numNt_' + str(numNtK) + '/'
            
            for i in range(1,nNtNested): 
                folderName += self.sequences[k][:i] + '/'
                
            if numNtK >= maxNumNtsNested:
                for i in range(numNtK // maxNumCharsInFolderName + 1):
                    folderName += self.sequences[k][i*maxNumCharsInFolderName:(i+1)*maxNumCharsInFolderName] + '/'
        
        folderName += 'minBPInStem_' + str(self.minBPInStem) + '/'
        
        #then make sure we're deliminating the right frozenBPs, by keeping track of the start and stop of each stem in frozenBPs
        frozenBPs = self.frozenBPs
        if frozenBPs: 
            #for example, if we only are freezing the stem [2 3 4 5 33 32 31 30]
            #frozenBPsAsStr will be '2-33_5-30' and will have a '__' before the next stem if it exists
            frozenStructBPs = bpsList2structBPs(frozenBPs)
            frozenBPsAsStr = ''
            for i in range(len(frozenStructBPs)): #for each stem
                frozenBPsAsStr += str(frozenStructBPs[i][0][0]) + '-' + str(frozenStructBPs[i][0][1])
                if len(frozenStructBPs[i]) > 1: #more than one bp in stem
                    frozenBPsAsStr += '_' + str(frozenStructBPs[i][-1][0]) + '-' + str(frozenStructBPs[i][-1][1])
                if i < len(frozenStructBPs) - 1:
                    frozenBPsAsStr += '__'
             
#            prevFrozenBP = frozenBPs[0]
#            frozenBPsAsStr += str(prevFrozenBP[0]) + '_' + str(prevFrozenBP[1]) + '__'
#            for i in range(1,np.shape(frozenBPs)[0]): #double counts stems of length 1, but that's fine as long as it's consistent
#                currFrozenBP = frozenBPs[i]
#                if not all(currFrozenBP - prevFrozenBP == [1,-1]): #i.e. if the next bp is not the continuation of a stem
#                    frozenBPsAsStr += str(prevFrozenBP[0]) + '_' + str(prevFrozenBP[1]) + '__'
#                    frozenBPsAsStr += str(currFrozenBP[0]) + '_' + str(currFrozenBP[1]) + '__'
#                prevFrozenBP = currFrozenBP #don't need to copy
#
#            frozenBPsAsStr += str(prevFrozenBP[0]) + '_' + str(prevFrozenBP[1])
            #for example, if we only are freezing the stem [2 3 4 5 33 32 31 30]
            #frozenBPsAsStr will be '2_33__5_30'.
            
            folderName += 'frozenBPs_' + frozenBPsAsStr + '/'

        if self.considerTabulatedLoops:
            folderName += 'considerTabulatedLoops_' + str(self.considerTabulatedLoops) + '/'
        if self.unmatchedBPPenalty:
            folderName += 'unmatchedBPPenalty_' + str(self.unmatchedBPPenalty) + '/'
        if self.onlyConsiderBondedStrands:
            folderName += 'onlyConsiderBondedStrands_' + str(self.onlyConsiderBondedStrands) + '/'
        if self.onlyAllowSubsetsOfLongestStems:
            folderName += 'onlyAllowSubsetsOfLongestStems_' + str(self.onlyAllowSubsetsOfLongestStems) + '/'
        if self.onlyConsiderSubstemsFromEdges:
            folderName += 'onlyConsiderSubstemsFromEdges_' + str(self.onlyConsiderSubstemsFromEdges) + '/'
        if not self.substems == 'all':
            folderName += 'substems_' + str(self.substems) + '/'
        if self.minNumStemsInStructure > 0:
            folderName += 'minNumStemsInStructure_' + str(self.minNumStemsInStructure) + '/'
        if not self.allowPseudoknots:
            folderName += 'allowPseudoknots_' + str(self.allowPseudoknots) + '/'
        if self.allowParallelStrands:
            folderName += 'allowParallel_' + str(self.allowParallelStrands) + '/'
        if not self.minNtsInHairpin == 3:
            folderName += 'minNtsInHairpin_' + str(self.minNtsInHairpin) + '/'
        if self.considerAllAsTerminalMismatches:
            folderName += 'considerAllAsTerminalMismatches_' + str(self.considerAllAsTerminalMismatches) + '/'
        if not self.includeTerminalMismatches:
            folderName += 'includeTerminalMismatches' + str(self.includeTerminalMismatches) + '/'
        if not self.includeFlushCoaxialStacks:
            folderName += 'includeFlushCoaxialStacks' + str(self.includeFlushCoaxialStacks) + '/'

        if not os.path.exists(folderName):
            os.makedirs(folderName)
        
# =============================================================================
#         Finally, now that we have the folder name, append the file name to it.
#         The choice of what to put in the file name rather than in the folders
#         is made so that we can load another file from the same folder and 
#         ideally only need to make minor changes in order to get what we want.
#         for example, changing the temperature is trivial once we've solved for
#         the energies and entropies.
# =============================================================================
        
        fileName = folderName + self.fileNameFxn(self.numSequences, str(self.pairwiseEnergies), 
                str(self.storeGraphs),str(self.T), str(self.vs), str(self.g), 
                str(self.duplexEntropyPenaltyInKB), str(self.corruptFESeed), 
                str(self.includeTerminalAUATPenalties), str(self.includeDanglingEnds),
                str(self.unboundButCouldBindPenalties[0]),str(self.unboundButCouldBindPenalties[1]))    

        fileTxt = fileName + 'output.txt' #make the textfile name to store the progress updates
        fileVars = fileName + 'vars' #make the variable file name to store the variables
        fileFigs = fileName + 'figs/' #make the folder name to store the figures
        
        return(fileVars,fileTxt,folderName,fileFigs)
        
    
    def fileNameFxn(self, numSequences, strPairwiseEnergies, strStoreGraphs, strT, 
                    strVs, strG, strDuplexEntropyPenaltyInKB, strCorruptFESeed,
                    strIncludeTerminalAUATPenalties, strIncludeDanglingEnds,
                    strUnboundButCouldBindPenaltiesE,strUnboundButCouldBindPenaltiesS):
        #Make the file name -- i.e. just the name of the file within the folder
        
        fileName = 'pairwiseEnergies_' + strPairwiseEnergies + '_storeGraphs_' + strStoreGraphs + \
                        '_T_' + strT + '_vs_' + strVs + '_g_' + strG + \
                        '_corrFESeed_' + strCorruptFESeed + \
                        '_incTermAUAT_' + strIncludeTerminalAUATPenalties + \
                        '_incDanglingEnds_' + strIncludeDanglingEnds + \
                        '_unboundButCouldBindPenaltiesE_' + strUnboundButCouldBindPenaltiesE + \
                        '_unboundButCouldBindPenaltiesS_' + strUnboundButCouldBindPenaltiesS + '_'
                        
        if numSequences > 1: #then duplexEntropyPenaltyInKB is important
            fileName += 'duplexEntropyPenalty_' + strDuplexEntropyPenaltyInKB + '_'  
            
        return(fileName)
        
    def profileFxn(self):
        cProfile.runctx('self.calculateFELandscape()', globals(), locals())
        
    def structureWithNthLowestFE(self,n): #Can't use this until after calculation is finished
#        the structure with the n^th lowest FE will be structures[indexSort[n]]
#        its FE is given by self.sortedFEs[n]
        whichStruct = self.indexSort[n]
        structure = [None]*len(self.structures[whichStruct])
        for i in range(len(structure)):
            structure[i] = self.STableStructure[self.structures[whichStruct][i]]
        return(structure)
        
    def whichStructure(self,structure):
#        Given a structure, find its index and position in self.sortedFEs
        structStems = []
        for stem in structure:
            stemI = [stemIndex for stemIndex in range(self.numStems) if self.STableStructure[stemIndex] == stem]
            if stemI:
                structStems.append(stemI[0])
            else:
                return([],[])
        structIndex = [i for i in range(self.numStructures) if all([x in self.structures[i] for x in structStems]) \
                       and all([x in structStems for x in self.structures[i]])]
        if structIndex:
            return(structIndex[0], np.where(self.indexSort == structIndex[0])[0][0]) #index in self.structures, index in self.sortedFEs
        else:
            return([],[])

# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================

# =============================================================================
# GENERAL FUNCTIONS
# =============================================================================
        
def twoListsShareElement(a,b):
# =============================================================================
#     Check if two lists, a and b, share any items. This is supposedly the fastest way to test this.
#    https://stackoverflow.com/questions/3170055/test-if-lists-share-any-items-in-python
#    Returns True if the lists share any elements, False otherwise
# =============================================================================
    return(not set(a).isdisjoint(b))

# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================

def isComplementary(x,y):
#    x,y define two ntds with the code RNA: A=1,C=2,G=3,U=4; DNA: A=5,C=6,G=7,T=8
#    unknown = 0
    x, y = sorted([x,y])
    if x == 1 and y == 4:
        return(True)
    elif x == 5 and y == 8:
        return(True)
    elif x == 1 and y == 8:
        return(True)
    elif x == 4 and y == 5:
        return(True)
    elif x == 2 and y == 3:
        return(True)
    elif x == 6 and y == 7:
        return(True)
    elif x == 2 and y == 7:
        return(True)
    elif x == 3 and y == 6:
        return(True)
    elif x == 3 and y == 4:
        return(True)
    #Can't have DNA G bind to RNA U
    
    return(False)
    
    
# =============================================================================
# FUNCTIONS NEEDED FOR FREE ENERGY CALCULATION
# =============================================================================
        
def bondFreeEnergiesRNARNA():
    #define the RNA/RNA bond enthalpy and entropy arrays
    
    #Sources: Table 4 of Xia et al Biochemistry '98
        #Table 4 of Mathews et al. JMB '99
        #Table 3 of Xia, Mathews, Turner "Thermodynamics of RNA secondary structure formation" in book by Soll, Nishmura, Moore
    #First index tells you if the first bp of the set is AU (0) CG (1) GC (2) UA (3) GU (4) or UG (5)
    #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3) (which row of table 1 in Serra & Turner).
    #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3) (which column of table 1 in Serra & Turner).
    
#    Had to make this 2D array to be able to use scipy sparse functions, so second/third index is replaced
#    by (4*second index) + third index
    bondEnergyMatrixRNARNA = np.array([[-3.9, 2, -3.5, -6.82, -2.3, 6, -11.4, -0.3, -3.1, -10.48, -3.5, -3.21, -9.38, 4.6, -8.81, -1.7],
                                       [-9.1, -5.6, -5.6, -10.44, -5.7, -3.4, -13.39, -2.7, -8.2, -10.64, -9.2, -5.61, -10.48, -5.3, -12.11, -8.6],
                                       [-5.2, -4, -5.6, -12.44, -7.2, 0.5, -14.88, -4.2, -7.1, -13.39, -6.2, -8.33, -11.4, -0.3, -12.59, -5],
                                       [-4, -6.3, -8.9, -7.69, -4.3, -5.1, -12.44, -1.8, -3.8, -10.44, -8.9, -6.99, -6.82, -1.4, -12.83, 1.4],
                                       [3.4, 2, -3.5, -12.83, -2.3, 6, -12.59, -0.3, 0.6, -12.11, -3.5, -13.47, -8.81, 4.6, -14.59, -1.7],
                                       [-4.8, -6.3, -8.9, -6.99, -4.3, -5.1, -8.33, -1.8, -3.1, -5.61, 1.5, -9.26, -3.21, -1.4, -13.47, 1.4]])
    #matrix of enthalpies (deltaH). Units are kcal/mol.
    
    bondEntropyMatrixRNARNA = np.array([[-10.2, 9.6, -8.7, -19, -5.3, 21.6, -29.5, 1.5, -7.3, -27.1, -8.7, -8.6, -26.7, 17.4, -24, -2.7],
                                        [-24.5, -13.5, -13.4, -26.9, -15.2, -7.6, -32.7, -6.3, -21.8, -26.7, -24.6, -13.5, -27.1, -12.6, -32.2, -23.9],
                                        [-13.2, -8.2, -13.9, -32.5, -19.6, 3.9, -36.9, -12.2, -17.8, -32.7, -15.1, -21.9, -29.5, -2.1, -32.5, -14],
                                        [-9.7, -17.1, -25.2, -20.5, -11.6, -14.6, -32.5, -4.2, -8.5, -26.9, -25, -19.3, -19, -2.5, -37.3, 6],
                                        [10, 9.6, -8.7, -37.3, -5.3, 21.6, -32.5, 1.5, 0, -32.2, -8.7, -44.9, -24, 17.4, -51.2, -2.7],
                                        [-12.1, -17.7, -25.2, -19.3, -11.6, -14.6, -21.9, -4.2, -11.2, -13.5, 2.1, -30.8, -8.6, -2.5, -44.9, 6]])
    #matrix of entropies (deltaS). Units are initially eu, but then converted to kcal/(mol*K).
     
    bondEntropyMatrixRNARNA /= 1000 #to convert from eu (entropy units) to kcal/(mol*K)

#    bondFreeEnergyMatrixRNARNA = np.zeros((6,4,4)) #matrix of stacking energies (deltaG). Units are kcal/mol.    
#    bondFreeEnergyMatrixRNARNA[0,:,:] = [[-0.8, -1.0, -0.8, -0.93], [-0.6, -0.7, -2.24, -0.7], [-0.8, -2.08, -0.8, -0.55], [-1.10, -0.8, -1.36, -0.8]]
#    bondFreeEnergyMatrixRNARNA[1,:,:] = [[-1.5, -1.5, -1.4, -2.11], [-1.0, -1.1, -3.26, -0.8], [-1.4, -2.36, -1.6, -1.41], [-2.08, -1.4, -2.11, -1.2]]
#    bondFreeEnergyMatrixRNARNA[2,:,:] = [[-1.1, -1.5, -1.3, -2.35], [-1.1, -0.7, -3.42, -0.5], [-1.6, -3.26, -1.4, -1.53], [-2.24, -1.0, -2.51, -0.7]]
#    bondFreeEnergyMatrixRNARNA[3,:,:] = [[-1.0, -0.8, -1.1, -1.33], [-0.7, -0.6, -2.35, -0.5], [-1.1, -2.11, -1.2, -1.00], [-0.93, -0.6, -1.27, -0.5]]
#    bondFreeEnergyMatrixRNARNA[4,:,:] = [[0.3, -1.0, -0.8, -1.27], [-0.6, -0.7, -2.51, -0.7], [0.6, -2.11, -0.8, -0.500], [-1.36, -0.8, 1.29, -0.8]] 
#    bondFreeEnergyMatrixRNARNA[5,:,:] = [[-1.0, -0.8, -1.1, -1.00], [-0.7, -0.6, -1.53, -0.5], [0.5, -1.41, 0.8, 0.30], [-0.55, -0.6, -0.500, -0.5]]
#    #the -0.500 was actually measured at +0.47 but the authors claim -0.5 is a better estimate.
    
    return(bondEnergyMatrixRNARNA, bondEntropyMatrixRNARNA)
    
def bondFreeEnergiesDNADNA():
    #define the DNA/DNA bond enthalpy and entropy arrays
    
    #Data is from From: Thermodynamics and NMR of Internal G·T Mismatches in DNA and other papers by Allawi and 
    #various SantaLucia publications (cited as 28-32 in Mfold web server for nucleic acid folding and hybridization prediction)
    #First index tells you if the first bp of the set is AT (0) CG (1) GC (2) or TA (3)
    #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or T(3)
    #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or T(3) 
    
#    Had to make this 2D array to be able to use scipy sparse functions, so second/third index is replaced
#    by (4*second index) + third index 
    
    bondFreeEnergyMatrixDNADNA = np.array([[0.61, 0.88, 0.14, -1.0, 0.77, 1.33, -1.44, 0.64, 0.02, -1.28, -0.13, 0.71, -0.88, 0.73, 0.07, 0.69],
                                           [0.43, 0.75, 0.03, -1.45, 0.79, 0.7, -1.84, 0.62, 0.11, -2.17, -0.11, -0.47, -1.28, 0.4, -0.32, -0.21],
                                           [0.17, 0.81, -0.25, -1.3, 0.47, 0.79, -2.24, 0.62, -0.52, -1.84, -1.11, 0.08, -1.44, 0.98, -0.59, 0.45],
                                           [0.69, 0.92, 0.42, -0.58, 1.33, 1.05, -1.3, 0.97, 0.74, -1.45, 0.44, 0.43, -1.0, 0.75, 0.34, 0.68]])
    #matrix of enthalpies (deltaH). Units are kcal/mol.
    
    bondEnergyMatrixDNADNA = np.array([[1.2, 2.3, -0.6, -7.9, 5.3, 0.0, -8.4, 0.7, -0.7, -7.8, -3.1, 1.0, -7.2, -1.2, -2.5, -2.7],
                                       [-0.9, 1.9, -0.7, -8.5, 0.6, -1.5, -8, -0.8, -4, -10.6, -4.9, -4.1, -7.8, -1.5, -2.8, -5],
                                       [-2.9, 5.2, -0.6, -8.2, -0.7, 3.6, -9.8, 2.3, 0.5, -8, -6, 3.3, -8.4, 5.2, -4.4, -2.2],
                                       [4.7, 3.4, 0.7, -7.2, 7.6, 6.1, -8.2, 1.2, 3, -8.5, 1.6, -0.1, -7.9, 1.0, -1.3, 0.2]])
    #matrix of free energies (deltaG). Units are kcal/mol.
    
    bondEntropyMatrixDNADNA = -(bondFreeEnergyMatrixDNADNA - bondEnergyMatrixDNADNA) / (273.15+37) 
    #matrix of entropies (deltaS). Units are converted to kcal/(mol*K) by dividing by 310.15 since the free energies were measured at 310.15 K
    
    return(bondEnergyMatrixDNADNA, bondEntropyMatrixDNADNA)


def bondFreeEnergiesRNADNA():
    #define the RNA/DNA bond enthalpy and entropy arrays
    
    #First index tells you if the set is (putting RNA first) 
    #5'AX3'/3'TY5' (0) 5CX3/3GY5 (1) 5GX3/3CY5 (2) 5UX3/3AY5 (3) 5XA3/3YT5 (4) 5XC3/3YG5 (5) 5XG3/3YC5 (6) or 5XU3/3YA5 (7).
    #Second index tells you if X is A (0) C (1) G (2) or U (3)
    #Third index tells you if Y is A (0) C (1) G (2) or T (3)
    #Data is from From: Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid Duplexes & 
    #Thermodynamic contributions of single internal rA·dA, rC·dC, rG·dG and rU·dT mismatches in RNA/DNA duplexes
    
#    Had to make this 2D array to be able to use scipy sparse functions, so second/third index is replaced
#    by  (4*second index) + third index
    
    #100 is put in place of elements for which there aren't published parameters
    bondFreeEnergyMatrixRNADNA = np.array([[1.07, 100, 100, -1.0, 100, 1.64, -2.1, 100, 100, -1.8, 0.31, 100, -0.9, 100, 100, 0.63],
                                           [0.90, 100, 100, -0.9, 100, 1.04, -2.1, 100, 100, -1.7, 0.14, 100, -0.9, 100, 100, 0.49],
                                           [0.51, 100, 100, -1.3, 100, 0.96, -2.7, 100, 100, -2.9, -0.58, 100, -1.1, 100, 100, 0.18],
                                           [1.13, 100, 100, -0.6, 100, 1.15, -1.5, 100, 100, -1.6, 0.44, 100, -0.2, 100, 100, 1.07],
                                           [1.36, 100, 100, -1.0, 100, 1.70, -0.9, 100, 100, -1.3, 0.50, 100, -0.6, 100, 100, 1.21],
                                           [0.19, 100, 100, -2.1, 100, 0.73, -2.1, 100, 100, -2.7, -0.83, 100, -1.5, 100, 100, -0.02],
                                           [0.21, 100, 100, -1.8, 100, 0.46, -1.7, 100, 100, -2.9, -0.33, 100, -1.6, 100, 100, 0.14],
                                           [1.85, 100, 100, -0.9, 100, 1.88, -0.9, 100, 100, -1.1, 0.97, 100, -0.2, 100, 100, 1.03]])
    #matrix of stacking free energies (deltaG). Units are kcal/mol.
    
    bondEnergyMatrixRNADNA = np.array([[-4.3, 100, 100, -7.8, 100, -8.8, -5.9, 100, 100, -9.1, -3.3, 100, -8.3, 100, 100, 0.6],
                                       [5.5, 100, 100, -9.0, 100, 10.5, -9.3, 100, 100, -16.3, -8.9, 100, -7.0, 100, 100, -0.4],
                                       [-1.9, 100, 100, -5.5, 100, -0.1, -8.0, 100, 100, -12.8, -8.0, 100, -7.8, 100, 100, -11.6],
                                       [-1.7, 100, 100, -7.8, 100, -3.3, -8.6, 100, 100, -10.4, -5.8, 100, -11.5, 100, 100, -2.2],
                                       [3.0, 100, 100, -7.8, 100, -0.3, -9.0, 100, 100, -5.5, 1.1, 100, -7.8, 100, 100, -3.3],
                                       [-6.0, 100, 100, -5.9, 100, 9.3, -9.3, 100, 100, -8.0, -7.0, 100, -8.6, 100, 100, 0.1],
                                       [-10.5, 100, 100, -9.1, 100, -11.5, -16.3, 100, 100, -12.8, -16.5, 100, -10.4, 100, 100, -13.4],
                                       [5.6, 100, 100, -8.3, 100, 0.8, -7.0, 100, 100, -7.8, -3.7, 100, -11.5, 100, 100, 3.0]])
    #matrix of enthalpies (deltaH). Units are kcal/mol.
    
    bondEntropyMatrixRNADNA = -(bondFreeEnergyMatrixRNADNA - bondEnergyMatrixRNADNA) / (273.15+37)
    
    #for elements for which there aren't published parameters, average the values for RNARNA and DNADNA bonds.
    bondEnergyMatrixDNADNA, bondEntropyMatrixDNADNA = bondFreeEnergiesDNADNA()
    bondEnergyMatrixRNARNA, bondEntropyMatrixRNARNA = bondFreeEnergiesRNARNA()
    
    for i in range(8):
        for j in range(4):
            for k in range(4):
                if bondEntropyMatrixRNADNA[i,k] == 0:
                    if i < 4:
                        #bondFreeEnergyMatrixRNADNA[i,j,k] = (bondFreeEnergyMatrixDNADNA[i,j,k] + bondFreeEnergyMatrixRNARNA[i,j,k]) / 2;
                        bondEnergyMatrixRNADNA[i,4*j+k] = (bondEnergyMatrixDNADNA[i,4*j+k] + bondEnergyMatrixRNARNA[i,4*j+k]) / 2;
                        bondEntropyMatrixRNADNA[i,4*j+k] = (bondEntropyMatrixDNADNA[i,4*j+k] + bondEntropyMatrixRNARNA[i,4*j+k]) / 2;
                    else: #flip the stem upside down to put the X,Y pair second. The order is then TA, GC, CG, AU 
                        #(i.e. the opposite order from the i<=4 case). Also, then need to switch j and k (i.e. switch X and Y).
                        
                        #bondFreeEnergyMatrixRNADNA[i,4*j+k] = (bondFreeEnergyMatrixDNADNA[7-i,4*k+j] + bondFreeEnergyMatrixRNARNA[7-i,4*k+j]) / 2;
                        bondEnergyMatrixRNADNA[i,4*j+k] = (bondEnergyMatrixDNADNA[7-i,4*k+j] + bondEnergyMatrixRNARNA[7-i,4*k+j]) / 2;
                        bondEntropyMatrixRNADNA[i,4*j+k] = (bondEntropyMatrixDNADNA[7-i,4*k+j] + bondEntropyMatrixRNARNA[7-i,4*k+j]) / 2;

    return(bondEnergyMatrixRNADNA, bondEntropyMatrixRNADNA)
    
def bondFreeEnergies():
    bondEnergyMatrixDNADNA, bondEntropyMatrixDNADNA = bondFreeEnergiesDNADNA()
    bondEnergyMatrixRNARNA, bondEntropyMatrixRNARNA = bondFreeEnergiesRNARNA() 
    bondEnergyMatrixRNADNA, bondEntropyMatrixRNADNA = bondFreeEnergiesRNADNA() 
    terminal_AU_penalty_energy, terminal_AU_penalty_entropy, \
        terminal_AT_penalty_energy, terminal_AT_penalty_entropy = terminalATAUPenalties()
    unknownEnergyMatrix = np.zeros((1,1)); unknownEntropyMatrix = np.zeros((1,1))
    
    energyMatrices = [bondEnergyMatrixRNARNA, bondEnergyMatrixDNADNA, bondEnergyMatrixRNADNA, 
                      np.array([[terminal_AU_penalty_energy, terminal_AT_penalty_energy]]), unknownEnergyMatrix]        
    entropyMatrices = [bondEntropyMatrixRNARNA, bondEntropyMatrixDNADNA, bondEntropyMatrixRNADNA, 
                       np.array([[terminal_AU_penalty_entropy, terminal_AT_penalty_entropy]]), unknownEntropyMatrix]

    return(energyMatrices, entropyMatrices)

def corruptBondFreeEnergies(corruptFESeed, engStd = 0.065, entStd = 0.073, FEStd = 0.024):
# =============================================================================
# Use a corrupted free energy table for the bond free energies. Keep track of which
#    corrupted version we used with corruptFESeed.
# Corrupt the FEs by changing \DeltaH, \DeltaS == dH, dS by random amounts drawn from Gaussians
#    with standard deviations engStd, entStd. Do this such that dG_37 = dH - (310.15)*dS
#    changes by a random amount with standard deviation FEStd.
    
# From Xia paper: (replaced \Delta H, etc. with dH, etc., and use dG to mean dG_37).
#    TM is the melting temperature.
#    
#    Methods for estimating errors in thermodynamic parameters have been described in detail.
# Instrumental fluctuations contribute negligibly to the uncertainties. Standard deviations of parameters for
# single measurements are typically about 6.5%, 7.3%, and 2.4% for dH, dS, and dG, respectively.
# Standard deviations for parameters calculated from Equation (8} based on 7-10 measurements are
# typically 2.9%, 3.3%, and 1.0% for dH, dS. and dG, respectively, The relative uncertainty in dS
# is usually about 13% larger than that in dH, because dS depends on more experimental parameters
# than dH. Uncertainty in TM is normally about 1.6°C. The errors in dG and TM are less than those in
# dH and dS because errors in dH and dS are highly correlated, with average observed correlation
# coefficients being greater than 0.999. Thus, dG and TM are more accurate parameters than either
# dG or dS individually.
    
# For the first sentence, Xia et al. cite two papers:
# Functional group substitutions as probes of hydrogen bonding between GA mismatches in RNA internal loops
#    by John SantaLucia Jr., Ryszard Kierzek, and Douglas H. Turner (1991 J. AM. Chem. Soc.)
# Thermodynamic parameters for an expanded nearest-neighbor model for formation of RNA duplexes with Watson-Crick base pairs.
#    by Xia T1, SantaLucia J Jr, Burkard ME, Kierzek R, Schroeder SJ, Jiao X, Cox C, Turner DH (1998 Biochemistry)
# =============================================================================
    
# =============================================================================
#     We corrupt the energy/entropy matrices by multiplying them by numbers drawn from a 2D Gaussian
#    For a multivariate normal with mean [mx, my], cov = [[sx**2, rho*sx*sy],[rho*sx*sy, sy**2]]
#    p(Ax + By) = \int p(x',y') \delta(Ax' + By' - (Ax + By)) = normal distribution
#    with mean A*mx + B*my, and variance (A*sx)**2 + (B*sy)**2 + 2*A*B*sx*sy*rho.
#    This means that a parameter with energy E, entropy S can be corrupted by multiplying
#    E by a random number drawn from a Gaussian with std sE, and S by a random number drawn 
#    from a Gaussian with std. sS. If we want G=E-TS to be changed by multiplication by a random number
#    drawn from a Gaussian with std. sG, we have that:
#    G*sG = sqrt((E*sE)**2 + (S*T*sS)**2 - 2*E*S*T*sE*sS*rho). Thus, 
#    rho = ((E*sE)**2 + (S*T*sS)**2 -(G*sG)**2) / (2*E*S*T*sE*sS).
#    Using the findRhos function we wrote below, we find that in order to get as low sG 
#    as Xia et al. claim, you need rho to be maxed out at 1 -- perfect correlation between
#    the corruption of E and that of S.
    
#    We therefore haven't yet included FEStd, though doing so is relatively straightforward
#    (we'd just need to make rho parameter-specific)
# =============================================================================
    
    np.random.seed(corruptFESeed)
    
    corruptEMatrices, corruptSMatrices = bondFreeEnergies()
    
    sE = engStd #0.065 according to Xia et al. for RNARNA
    sS = entStd #0.073 according to Xia et al. for RNARNA
    rho = 1 #between -1 and 1 to make cov positive semidefinite
#    Using the findRhos function, we find that, at least for RNARNA, the energy/entropy measurement
#    errors have to be perfectly correlated to get as low free energy measurement error as Xia et al. claim.
#   If you want rho to be parameter specific, move this into the for loop.

    sES = sE*sS*rho #covariance in x,y 
    cov = [[sE**2,sES],[sES,sS**2]] #must be symmetric and positive semidefinite
    
    for bpType in range(len(corruptEMatrices)):
        rs = np.random.multivariate_normal([1,1], cov, size = corruptEMatrices[bpType].shape)
#            means are [1,1] since we're going to multiply the random numbers by the energy/entropy matrices
        corruptEMatrices[bpType] *= rs[:,:,0]
        corruptSMatrices[bpType] *= rs[:,:,1]
        
    return(corruptEMatrices, corruptSMatrices)
    
def findRhos(bondEnergyMatrix, bondEntropyMatrix, engStd, entStd, FEStd):
# =============================================================================
#     for each tabulated energy, entropy, if the error in the energy measurement = engStd
#     and the error in the entropy measurement = entStd, what is the required correlation coefficient
#     to get the error in the Free Energy (at 37C) to be FEStd?
# =============================================================================
    Es = [bondEnergyMatrix[i,j] for i in range(bondEnergyMatrix.shape[0]) 
        for j in range(bondEnergyMatrix.shape[1])]
    Ss = [bondEntropyMatrix[i,j] for i in range(bondEnergyMatrix.shape[0]) 
        for j in range(bondEnergyMatrix.shape[1])]
    Gs = [Es[i] - 310.15*Ss[i] for i in range(len(Es))]
    sE = engStd #0.065 according to Xia et al. for RNARNA
    sS = entStd #0.073 according to Xia et al. for RNARNA
    sG = FEStd #0.024 according to Xia et al. for RNARNA

    #what rhos would we need to achieve this sG value for all the different parameters?
    rhos = [(E**2 * sE**2 + S**2 * 310.15**2 * sS**2 - G**2 * sG**2) / (2*E*S*310.15*sE*sS) for E,S,G in zip(Es, Ss, Gs)]
    return(rhos)
#    as you can see, the rho values are all >= 0.999 except for two which are -0.90, for the RNARNA parameters
#    rhos < -1 or >1 are unphysical.
    
    
def danglingEndMatrices():
    #Dangling ends also have an associated free energy. 

    #first index tells you if the paired ntd is an A(1) C(2) G(3) or U(4);
    #by paired ntd, we mean the paired ntd either 1 before or 1 after the
    #dangling ntd (on the same strand as the dangling end)
    #second index tells you if the dangling ntd is A(1) C(2) G(3) or U(4). 
    
    #Data taken from Xia, Mathews, Turner "Thermodynamics of RNA secondary structure 
    #formation" in book by Soll, Nishmura, Moore.
    #It is the same as the data from Serra & Turner Table 2 (second half of the table).
    
    #GU base pairs are treated here as if they were AU (so G is replaced with A
    #if the base pair is GU).
    
    #For a 5' dangling end:
#    dangling5FreeEnergyRNARNA = np.array([[-0.3, -0.3, -0.4, -0.2], [-0.5, -0.3, -0.2, -0.1], [-0.2, -0.3, -0.0, -0.0], [-0.3, -0.1, -0.2, -0.2]])
    dangling5EnergyRNARNA = np.array([[1.6, 2.2, 0.7, 3.1], [-2.4, 3.3, 0.8, -1.4], [-1.6, 0.7, -4.6, -0.4], [-0.5, 6.9, 0.6, 0.6]])
    dangling5EntropyRNARNA = np.array([[6.1, 7.9, 3.4, 10.6], [-6.0, 11.8, 3.4, -4.3], [-4.5, 3.1, -14.8, -1.2], [-0.7, 22.8, 2.7, 2.7]]) / 1000 
        #divide by 1000 to convert from eu to kcal/(mol*K)
    
    #For a 3' dangling end:
#    dangling3FreeEnergyRNARNA = np.array([[-0.8, -0.5, -0.8, -0.6], [-1.7, -0.8, -1.7, -1.2], [-1.1, -0.4, -1.3, -0.6], [-0.7, -0.1, -0.7, -0.1]])
    dangling3EnergyRNARNA = np.array([[-4.9, -0.9, -5.5, -2.3], [-9.0, -4.1, -8.6, -7.5], [-7.4, -2.8, -6.4, -3.6], [-5.7, -0.7, -5.8, -2.2]])
    dangling3EntropyRNARNA = np.array([[-13.2, -1.2, -15.0, -5.4], [-23.4, -10.7, -22.2, -20.4], [-20.0, -7.9, -16.6, -9.7], [-16.4, -1.8, -16.4, -6.8]]) /1000 
        #divide by 1000 to convert from eu to kcal/(mol*K)
    
    
    #DNA dangling end parameters taken from "Thermodynamic parameters for DNA
    #sequences with dangling ends". All energies in kcal/mol, so entropies are
    #in kcal/(mol*K). All parameters were checked with The Thermodynamics of
    #DNA Structural Motifs by SantaLucia and Hicks
    
    #first index tells you if the paired ntd is an A(1) C(2) G(3) or T(4);
    #second index tells you if the dangling ntd is A(1) C(2) G(3) or T(4).
    
    #For a 5' dangling end:
    dangling5FreeEnergyDNADNA = np.transpose(np.array([[-0.51, -0.96, -0.58, -0.5], [-0.42, -0.52, -0.34, -0.02], [-0.62, -0.72, -0.56, 0.48], [-0.71, -0.58, -0.61, -0.10]]))
    dangling5EnergyDNADNA = np.transpose(np.array([[0.2, -0.62, -3.7, -2.9], [0.6, -4.4, -4.0, -4.1], [-1.1, -5.1, -3.9, -4.2], [-6.9, -4.0, -4.9, -0.2]])) 
        #transpose because I copied this from excel doc which had matrix transposed
        #in  ... -6.2 is replaced with -6.3 but I believe it's a mistake because
        #they cite the same paper I used and it's the only difference
    dangling5EntropyDNADNA = -(dangling5FreeEnergyDNADNA - dangling5EnergyDNADNA) / (273.15+37)
    
    #For a 3' dangling end:
    dangling3FreeEnergyDNADNA = np.transpose(np.array([[-0.12, -0.82, -0.92, -0.48], [0.28, -0.31, -0.23, -0.19], [-0.01, -0.01, -0.44, -0.50], [0.13, -0.52, -0.35, -0.29]]))
    dangling3EnergyDNADNA = np.transpose(np.array([[-0.5, -5.9, -2.1, -0.7], [4.7, -2.6, -0.2, 4.4], [-4.1, -3.2, -3.9, -1.6], [-3.8, -5.2, -4.4, 2.9]]))
    dangling3EntropyDNADNA = -(dangling3FreeEnergyDNADNA - dangling3EnergyDNADNA) / (273.15+37)

    unknownDanglingMatrix = np.zeros((1,1))
    
    dangling5Energies = [dangling5EnergyRNARNA, dangling5EnergyDNADNA, unknownDanglingMatrix]
    dangling5Entropies = [dangling5EntropyRNARNA, dangling5EntropyDNADNA, unknownDanglingMatrix]
    dangling3Energies = [dangling3EnergyRNARNA, dangling3EnergyDNADNA, unknownDanglingMatrix]
    dangling3Entropies = [dangling3EntropyRNARNA, dangling3EntropyDNADNA, unknownDanglingMatrix]
    
    return(dangling5Energies, dangling5Entropies, dangling3Energies, dangling3Entropies)


def corruptDanglingEndMatrices(corruptFESeed, engStd = 0.065, entStd = 0.073, FEStd = 0.024):

    np.random.seed(corruptFESeed)
    
    dangling5Energies, dangling5Entropies, dangling3Energies, dangling3Entropies = danglingEndMatrices()
    
    sE = engStd #0.065 according to Xia et al. for RNARNA
    sS = entStd #0.073 according to Xia et al. for RNARNA
    rho = 1 #between -1 and 1 to make cov positive semidefinite
#    Using the findRhos function, we find that, at least for RNARNA, the energy/entropy measurement
#    errors have to be perfectly correlated to get as low free energy measurement error as Xia et al. claim.
#   If you want rho to be parameter specific, move this into the for loop.

    sES = sE*sS*rho #covariance in x,y 
    cov = [[sE**2,sES],[sES,sS**2]] #must be symmetric and positive semidefinite
    
    for bpType in range(len(dangling5Energies)):
        rs = np.random.multivariate_normal([1,1], cov, size = dangling5Energies[bpType].shape)
    #            means are [1,1] since we're going to multiply the random numbers by the energy/entropy matrices
        dangling5Energies[bpType] *= rs[:,:,0]
        dangling5Entropies[bpType] *= rs[:,:,1]
        
        rs = np.random.multivariate_normal([1,1], cov, size = dangling3Energies[bpType].shape)
    #            means are [1,1] since we're going to multiply the random numbers by the energy/entropy matrices
        dangling3Energies[bpType] *= rs[:,:,0]
        dangling3Entropies[bpType] *= rs[:,:,1]
    
    return(dangling5Energies, dangling5Entropies, dangling3Energies, dangling3Entropies)
    
def terminalATAUPenalties():
    terminal_AT_penalty_energy = 2.2; terminal_AT_penalty_entropy = (2.2-0.05)/310.15
    #from SantaLucia and Hicks review The Termodynamicso of DNA Structural Motifs 
    
    terminal_AU_penalty_energy = 3.72; terminal_AU_penalty_entropy = (3.72 - 0.45)/310.15
    #from Xia et al. 1998
    
    return(terminal_AU_penalty_energy, terminal_AU_penalty_entropy, 
             terminal_AT_penalty_energy, terminal_AT_penalty_entropy)
    
    
def danglingMatrixIndices(sequenceInNumbers, pairedNtIndex, danglingNtIndex):
# =============================================================================
#     Get the indices of the energy/entropy matrix to use for a dangling end#    
#    where pairedNt is the bonded nt on the same strand as the dangling end, and danglingNt is the dangling end
# =============================================================================
    numNt = len(sequenceInNumbers)
    
    if danglingNtIndex < 0 or danglingNtIndex > numNt - 1: #pairedNt must be all right
        pairedNt = 0; danglingNt = 0 #if anything is disallowed, return zero
    else:
        pairedNt = sequenceInNumbers[pairedNtIndex] #the identity of the nt (A is 1, C is 2, etc.)
        danglingNt = sequenceInNumbers[danglingNtIndex]

    if pairedNt == 0 or danglingNt == 0:
        return([0,0,0], 2) #if one of the sequence elements is unknown, return bpType = 2
#        This also means we don't need to worry about whether or not the base pairs sent
#       correspond to the linker or something.
        
    if danglingNt <= 4: 
        bpType = 0 #RNARNA
    elif danglingNt > 4:
        bpType = 1 #DNADNA
#    This doesn't deal well with RNA/DNA hybrids (it uses the parameters of whichever
#    ntd is dangling) because I couldn't find experimental estimates of the relevant parameters.

        #need -1 for RNA nts since nts are numbered 1-4 but indices are 0-3
        #need -5 for DNA nts since nts are numbered 5-8 but indices are 0-3    
    return([pairedNt -4*bpType - 1, danglingNt -4*bpType - 1], bpType)


#def deltaGforBPsAndTerminalMismatches_Fxn(bondEnergyMatrix,bondEntropyMatrix,firstNt,firstBP,secondNt,secondBP,bpType):
## =============================================================================
##     Calculate bond energy and bond entropy for internal base pairs and terminal mismatches
## =============================================================================
#    
#    bondEnergy = 0
#    bondEntropy = 0
#    
#    if firstNt == 0 or firstBP == 0 or secondNt == 0 or secondBP == 0:
#        return(bondEnergy,bondEntropy)
#        #if one of the sequence elements is unknown, don't do anything
#        
#    elif bpType == 0: #we're dealing with an RNA/RNA bond
#        #For the bondEnergyMatrix,
#        #First index tells you if the first bp of the set is AU (0) CG (1) GC (2) UA (3) GU (4) or UG (5)
#        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3)
#        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3)
#        if firstNt == 1 and firstBP == 4:
#            basePair = 0
#        elif firstNt == 2 and firstBP == 3:
#            basePair = 1
#        elif firstNt == 3 and firstBP == 2:
#            basePair = 2
#        elif firstNt == 4 and firstBP == 1:
#            basePair = 3
#        elif firstNt == 3 and firstBP == 4:
#            basePair = 4
#        elif firstNt == 4 and firstBP == 3:
#            basePair = 5
#        else:
#            print('The deltaG function has a problem RNA/RNA')
#            
#        bondEnergy += bondEnergyMatrix[basePair, secondNt - 1, secondBP - 1]
#        bondEntropy += bondEntropyMatrix[basePair, secondNt - 1, secondBP - 1]
#        #need -1 for RNA nts since nts are numbered 1-4 but indices are 0-3
#        
#    elif bpType == 1: #we're dealing with a DNA/DNA pair
#        #For the bondEnergyMatrix. We're dealing with DNA/DNA pair
#        #First index tells you if the first bp of the set is AT (0) CG (1) GC (2) or TA (3)
#        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or T(3)
#        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or T(3)
#        if firstNt == 5 and firstBP == 8:
#            basePair = 0
#        elif firstNt == 6 and firstBP == 7:
#            basePair = 1
#        elif firstNt == 7 and firstBP == 6:
#            basePair = 2
#        elif firstNt == 8 and firstBP == 5:
#            basePair = 3
#        else:
#            print('The deltaG function has a problem DNA/DNA')
#        bondEnergy += bondEnergyMatrix[basePair, secondNt - 5, secondBP - 5]
#        bondEntropy += bondEntropyMatrix[basePair,secondNt - 5,secondBP - 5]
#        #need -5 for DNA nts since nts are numbered 5-8 but indices are 0-3
#    
#    elif bpType == 2: #we're dealing with an RNA/DNA bond
#        #For the bondEnergyMatrix,
#        #First index tells you if the first bp of the set is AT (0) CG (1) GC (2) UA (3) or any of the reverse. 
#        #firstNt and secondNt are fixed to be RNA while firstBP and secondBP are DNA. Therefore the order of which
#        #base pair is first or second is in some cases switched.
#        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3)
#        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3)
#        if firstNt == 1 and firstBP == 8:
#            basePair = 0
#        elif firstNt == 2 and firstBP == 7:
#            basePair = 1
#        elif firstNt == 3 and firstBP == 6:
#            basePair = 2
#        elif firstNt == 4 and firstBP == 5:
#            basePair = 3
#        elif secondNt == 1 and secondBP == 8:
#            basePair = 4
#        elif secondNt == 2 and secondBP == 7:
#            basePair = 5
#        elif secondNt == 3 and secondBP == 6:
#            basePair = 6
#        elif secondNt == 4 and secondBP == 5:
#            basePair = 7
#        else:
#            print('The deltaG function has a problem RNA/DNA')
#
#        bondEnergy += bondEnergyMatrix[basePair,secondNt - 1,secondBP - 5]
#        bondEntropy += bondEntropyMatrix[basePair,secondNt - 1,secondBP - 5]
#        #need -1 for RNA nts since nts are numbered 1-4 but indices are 0-3
#        #need -5 for DNA nts since nts are numbered 5-8 but indices are 0-3
#        
#    
#    return(bondEnergy, bondEntropy)
    
def freeEnergyMatrixIndices(sequenceInNumbers, firstNtIndex,firstBPIndex, secondNtIndex,
                            secondBPIndex, bound = [1,1], unmatchedBPPenalty = True):
# =============================================================================
#     Get the indices of the energy/entropy matrix to use for
#    the base pair stack:
#    
#    5' firstNt, secondNt 3'
#    3' firstBP, secondBP 5'
#    
#    where firstNt is sequenceInNumbers[firstNtIndex], etc.
#    
#    Assumes firstNt and firstBP are actually bound, but secondNt and secondBP need not be
#    bound = [are_firstNt_and_firstBP_bound, are_secondNt_and_secondBP_bound]
#    unmatchedBPPenalty = True if two nts that could be bound but aren't should be treated as
#    an A-C pair. On top of whether or not we treat them as an A-C pair, we introduce an 
#    [energy, entropy] penalty given by unboundButCouldBindPenalties
# =============================================================================
    if not bound[0]: #if we accidentally have firstNt and firstBP not bound, flip the stem
        realSecondBPIndex = copy.copy(firstNtIndex); realFirstBPIndex = copy.copy(secondNtIndex)
        realSecondNtIndex = copy.copy(firstBPIndex); realFirstNtIndex = copy.copy(secondBPIndex)
        firstNtIndex = realFirstNtIndex; firstBPIndex = realFirstBPIndex
        secondNtIndex = realSecondNtIndex; secondBPIndex = realSecondBPIndex
    
    numNt = len(sequenceInNumbers)
    
    if (firstNtIndex < 0 or firstNtIndex > numNt - 1 or secondNtIndex < 0 or secondNtIndex > numNt - 1 or 
            firstBPIndex < 0 or firstBPIndex > numNt - 1 or secondBPIndex < 0 or secondBPIndex > numNt - 1):
        firstNt = 0; firstBP = 0; secondNt = 0; secondBP = 0 #if anything is disallowed, return zero
    else:
        firstNt = sequenceInNumbers[firstNtIndex] #the identity of the nt (A is 1, C is 2, etc.)
        firstBP = sequenceInNumbers[firstBPIndex]
        secondNt = sequenceInNumbers[secondNtIndex]
        secondBP = sequenceInNumbers[secondBPIndex]

    if firstNt == 0 or firstBP == 0 or secondNt == 0 or secondBP == 0:
        return([0,0,0], 4) #if one of the sequence elements is unknown, return bpType = 4
#        This also means we don't need to worry about whether or not the base pairs sent
#       correspond to the linker or something.
    
#    if base pair could bind but aren't bound in this structure, Lu, Turner, Mathews (NAR, 2006) say 
#    that we should treat them as an AC pair (where A replaces the purine and C the pyrimidine)
    if unmatchedBPPenalty and not bound[1]:
        if ((secondNt == 1 and secondBP == 4) or (secondNt == 4 and secondBP == 1) or
                (secondNt == 2 and secondBP == 3) or (secondNt == 3 and secondBP == 2) or
                (secondNt == 3 and secondBP == 4) or (secondNt == 4 and secondBP == 3)):
            if secondNt == 1 or secondNt == 3: #if secondNt is a purine
                secondNt = 1 #the purine is replaced with A 
                secondBP = 2 #and the pyrimidine with C
            else:
                secondNt = 2
                secondBP = 1

# =============================================================================
#     Find type of bond
# =============================================================================
    if firstNt <=4 and firstBP <= 4 and secondNt <=4 and secondBP <=4: 
        bpType = 0 #RNARNA

    elif firstNt > 4 and firstBP > 4 and secondNt > 4 and secondBP > 4:
        bpType = 1 #DNADNA
                
    elif firstNt <=4 and firstBP > 4 and secondNt <=4 and secondBP > 4:
        bpType = 2 #RNADNA
        
    elif firstNt > 4 and firstBP <=4 and secondNt > 4 and secondBP <=4:
        bpType = 2 #RNADNA
#       need to flip the stack upside down (i.e. just look at it as if the second strand were the first.)
        realSecondBP = copy.copy(firstNt); realFirstBP = copy.copy(secondNt)
        realSecondNt = copy.copy(firstBP); realFirstNt = copy.copy(secondBP)
        firstNt = realFirstNt; firstBP = realFirstBP; secondNt = realSecondNt; secondBP = realSecondBP
        
#    elif firstNt > 4 and firstBP <=4 and isParallel and secondNt > 4 and secondBP <=4:
#        bpType = 1 #RNADNA
#        bondEnergy,bondEntropy =  deltaGforBPsAndTerminalMismatches_Fxn(
#               bondEnergyMatrixRNADNA,bondEntropyMatrixRNADNA,firstBP,firstNt,secondBP,secondNt,bpType)
    
# =============================================================================
#     Get indices of matrix for specific pair of bps
# =============================================================================
    
    if bpType == 0: #we're dealing with an RNA/RNA bond
        #For the bondEnergyMatrix,
        #First index tells you if the first bp of the set is AU (0) CG (1) GC (2) UA (3) GU (4) or UG (5)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3)
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3)
        if firstNt == 1 and firstBP == 4:
            basePair = 0
        elif firstNt == 2 and firstBP == 3:
            basePair = 1
        elif firstNt == 3 and firstBP == 2:
            basePair = 2
        elif firstNt == 4 and firstBP == 1:
            basePair = 3
        elif firstNt == 3 and firstBP == 4:
            basePair = 4
        elif firstNt == 4 and firstBP == 3:
            basePair = 5
        else:
            print('The deltaG function has a problem RNA/RNA')
        
        return([basePair, 4*(secondNt - 1) + (secondBP - 1)], bpType)
        #need -1 for RNA nts since nts are numbered 1-4 but indices are 0-3
#        Had to make the energy/entropy matrices 2D which is why secondNt and secondBP are
#        not being returned as separate indices, but as one index (4*secondNt + secondBP)
        
    elif bpType == 1: #we're dealing with a DNA/DNA pair
        #For the bondEnergyMatrix. We're dealing with DNA/DNA pair
        #First index tells you if the first bp of the set is AT (0) CG (1) GC (2) or TA (3)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or T(3)
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or T(3)
        if firstNt == 5 and firstBP == 8:
            basePair = 0
        elif firstNt == 6 and firstBP == 7:
            basePair = 1
        elif firstNt == 7 and firstBP == 6:
            basePair = 2
        elif firstNt == 8 and firstBP == 5:
            basePair = 3
        else:
            print('The deltaG function has a problem DNA/DNA')
            
        return([basePair, 4*(secondNt - 5) + (secondBP - 5)], bpType)
        #need -5 for DNA nts since nts are numbered 5-8 but indices are 0-3
#        Had to make the energy/entropy matrices 2D which is why secondNt and secondBP are
#        not being returned as separate indices, but as one index (4*secondNt + secondBP)

    elif bpType == 2: #we're dealing with an RNA/DNA bond
        #For the bondEnergyMatrix,
#        First index tells you if the set is (putting RNA first) (0) 5'AX3'/3'TY5' (1) 5CX3/3GY5 
#       (2) 5GX3/3CY5 (3) 5UX3/3AY5 (4) 5XA3/3YT5 (5) 5XC3/3YG5 (6) 5XG3/3YC5 or (7) 5XU3/3YA5.

        #firstNt and secondNt are fixed to be RNA while firstBP and secondBP are DNA. Therefore the order of which
        #base pair is first or second is in some cases switched (which is why we may have had to flip the stem
#        upside down a minute ago)
        #Second index tells you if the 3' ntd of the second bp is A (0) C (1) G(2) or U(3)
        #Third index tells you if the 5' ntd of the second bp is A (0) C (1) G(2) or U(3)

        if firstNt == 1 and firstBP == 8:
            basePair = 0
        elif firstNt == 2 and firstBP == 7:
            basePair = 1
        elif firstNt == 3 and firstBP == 6:
            basePair = 2
        elif firstNt == 4 and firstBP == 5:
            basePair = 3
        elif secondNt == 1 and secondBP == 8:
            basePair = 4
        elif secondNt == 2 and secondBP == 7:
            basePair = 5
        elif secondNt == 3 and secondBP == 6:
            basePair = 6
        elif secondNt == 4 and secondBP == 5:
            basePair = 7
        else:
            print('The deltaG function has a problem RNA/DNA')

        return([basePair,4*(secondNt - 1) + (secondBP - 5)], bpType)
        #need -1 for RNA nts since nts are numbered 1-4 but indices are 0-3
        #need -5 for DNA nts since nts are numbered 5-8 but indices are 0-3
#        Had to make the energy/entropy matrices 2D which is why secondNt and secondBP are
#        not being returned as separate indices, but as one index (4*secondNt + secondBP)

# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================


# =============================================================================
# Converting between different ways to express secondary structure.
# We define the following, and give for each an example:
# (The example is for a structure composed of two stems: the first stem has ntds 
# 1 bound to 10, 2 to 9, and 3 to 8; the second has 15 bound to 25, 16 to 24, 17 to 23, and 18 to 22.):
#    structure = list of stems, where each stem is represented by a single list giving first the 
#       first strand, and then in order, the ntds they're bound to.
#       Example: [[1,2,3,10,9,8],[15,16,17,18,25,24,23,22]]
#    structBPs = list of stems, where each stem is represented by a list of base pairs in that stem,
#       and each base pair is represented as a list of the two ntds comprising it 
#       (in increasing order along the chain).
#       Example: [[[1,10],[2,9],[3,8]],[[15,25],[16,24],[17,23],[18,22]]]
#    bpsList = list of base pairs, where each base pair is represented as a list of the two nts
#       comprising it (in increasing order along the chain)
#    Finally, for the graph entropy calculation, it is useful to convert from a structure
#    (defined by any of the above) to a list of bound nts. From there, it is useful to make 
#    a list of the consecutive unpaired nts
# =============================================================================
    
def bpsList2structure(bpsList):
# =============================================================================
#     Given a list of base pairs, make them in "structure" format.
# =============================================================================
    bpsList = sorted(bpsList)
    structure = []
    startOfStem = False
    stem1 = [] #first nts of stem
    stem2 = [] #complementary nts of stem
    for i in range(len(bpsList)): #for each base pair
        firstBP = bpsList[i][0]
        secondBP = bpsList[i][1]
        
        if i: #i.e. if it's not the first base pair being considered
            #check if it's a continuation of the previous stem.
            if firstBP - stem1[-1] == 1: 
                if secondBP - stem2[-1] == 1: #we're dealing with a parallel stem. 
                    startOfStem = False
                elif secondBP - stem2[-1] == -1: #antiparallel (normal) stem
                    startOfStem = False
                else:
                    startOfStem = True
            else:
                startOfStem = True
        
        if not startOfStem:
            stem1.append(firstBP) 
            stem2.append(secondBP)
        else:
            structure.append(stem1 + stem2)
            stem1 = [firstBP]
            stem2 = [secondBP]
    
        #if we've reached the end, make sure to put in the last base pair
        if i == len(bpsList) - 1:
            structure.append(stem1 + stem2)
    return(structure)
    
def bpsList2structBPs(bpsList):
# =============================================================================
#     Given a list of base pairs, put them in structBPs format, meaning take the list of base pairs
#     and make a list of stems, where each stem is itself a list of base pairs in that stem.
#     For exmample, a bpList would be [[1,10],[2,9],[3,8],[15,25],[16,24],[18,22],[17,23]]
#     The corresponding structBPs would be [[[1,10],[2,9],[3,8]],[[15,25],[16,24],[17,23],[18,22]]]
#     and the corresponding structure would be [[1,2,3,10,9,8],[15,16,17,18,25,24,23,22]]
# =============================================================================
    
    bpsList = sorted(bpsList)
    structBPs = [[]]
    startOfStem = False
    for i in range(len(bpsList)): #for each base pair
        firstBP = bpsList[i][0]
        secondBP = bpsList[i][1]
        
        if i: #i.e. if it's not the first base pair being considered
            #check if it's a continuation of the previous stem.
            if firstBP - structBPs[-1][-1][0] == 1: 
                if secondBP - structBPs[-1][-1][1] == 1: #we're dealing with a parallel stem. 
                    startOfStem = False
                elif secondBP - structBPs[-1][-1][1] == -1: #antiparallel (normal) stem
                    startOfStem = False
                else:
                    startOfStem = True
            else:
                startOfStem = True
                
        if not startOfStem:
            structBPs[-1].append([firstBP,secondBP])
        else:
            structBPs.append([[firstBP,secondBP]])
    return(structBPs)

def structBPs2structure(structBPs):
    structure = []
    for i in range(len(structBPs)): #for each stem
        stem1 = [] #first nts of stem
        stem2 = [] #complementary nts of stem
        for j in range(len(structBPs[i])): #for each base pair
            stem1.append(structBPs[i][j][0])
            stem2.append(structBPs[i][j][1])
        structure.append(stem1 + stem2)
    return(structure)

def structBPs2bpsList(structBPs):
    return([bp for stem in structBPs for bp in stem])
    
def structure2structBPs(structure):
# =============================================================================
#     Take a structure format and make it into a structBPs format. This is a mega-list of
#    sub-lists where each sub-list is a list of base pairs that belong to a single stem
#    and the mega-list is a list of the stems (again, where each stem is expressed as a list of base pairs)
# =============================================================================
    structBPs = [] #a list of stems where each stem is in base-pair (bp) format
    for i in range(len(structure)): #for each stem
        structBPs.append([])
        for j in range(int(len(structure[i])/2)): #for each base pair in the stem
            structBPs[i].append([structure[i][j],structure[i][j + int(len(structure[i])/2)]])
            
    return(structBPs)


def structure2bpsList(structure):
    return(structBPs2bpsList(structure2structBPs(structure)))


def structure2boundNt(structure, numNt):
# =============================================================================
#     Return a list of the nts that are bound where the total length of the sequence is numNts
# =============================================================================
    return([i for i in range(numNt) for j in structure if i in j])

def structBPs2boundNt(structBPs, numNt):
    return([i for i in range(numNt) for j in structBPs for k in j if i in k])
    
def bpsList2boundNt(bpsList, numNt):
    return([i for i in range(numNt) for j in bpsList if i in j])
    




# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================

# =============================================================================
# Graph functions
# =============================================================================

def myNumConnCompFromMat(M): #M is adjacency matrix of a graph in np.array form
# =============================================================================
#  From stackexchange discussion: If you put all 1 on the diagonal of your adjacency matrix A, 
#    and all edge weights are positive then when you multiply A2=A*A you get a non-zero entry 
#    aij in A2 if and only if there exist non-zero aik and akj in A for some k, i.e. there 
#    is a path of length 2 between i and j if k≠j and k≠i and there is a path of length 1 if k=j or k=i. 
#    So the non-zero entries in A2 tell you all pairs of nodes that are connected by a path of length 1 or 2. 
#    Similarly the entries in Ak tell you all pairs of nodes that are connected by a path of length k or less. 
#    So if you start with A and keep squaring until you get Ak where k≥n where n is the number of nodes, 
#    then the non-zero entries in row i tell you all the nodes that are connected to node i 
#    (since two connected nodes must be connected by a path of length n or less). 
#    So if you have a row in Ak that is all non-zero, then the graph is connected. 
#    If the graph is not connected, you can similarly tell the connected components from the rows of Ak.
    
#    M should have ones on the diagonals (or any nonzero number really)
# =============================================================================
    n = len(M)
    
    power = int(2**np.ceil(np.log2(n)))
    
    C = np.linalg.matrix_power(M,power) #takes roughly half the time of the code
    
#    connectedNodes = []
#    i = 0
#    while i < n:
#        if not any(i in connComp for connComp in connectedNodes):
#            connectedNodes.append([j for j in range(n) if C[i,j]])
#        i += 1
#    return(len(connectedNodes))

    return(np.round(np.sum(1/np.sum(C,axis=0)))) 
    
def myNumConnCompFromMatC(M, edges): #same principle as above
# =============================================================================
# Given a matrix M and a set of edges to try, try removing each edge, and if doing so
# disconnects M, then add that to the list of bridges. 
#    This code was written to be simple to translate into C
#    Requires M to be an array of ones and zeros.
# =============================================================================

    n = len(M)
    nEdges = edges.shape[0] #or 1, not sure
    power = int(2**np.ceil(np.log2(n))) #want to raise it to a power >= n, and it's faster to raise
#       to a power that's a multiple of 2, since the code raises to a power by successive self-multiplications.
    
    bridgeList = np.zeros(nEdges,dtype=int) - 1 #
    nBridges = 0
    
    C = np.linalg.matrix_power(M,power)
    nCC = np.round(np.sum(1/np.sum(C,axis=0)))
    for i in range(nEdges):
        edge0 = edges[i][0]
        edge1 = edges[i][1]
        M2 = copy.copy(M)
        M2[edge0, edge1] = 0; M2[edge1, edge0] = 0
        C = np.linalg.matrix_power(M2,power) 
        
        if np.round(np.sum(1/np.sum(C,axis=0))) > nCC:
            bridgeList[nBridges] = i
            nBridges += 1
            
    return(bridgeList) 
    
    
    
def drawGraph(G,figsize = (7,7), drawTopology=True, drawEdgeLengths = True):
# =============================================================================
#     Assumes edges can either have between 0 and 2 single bonds and 0 to 1 double bonds
#    Edges with one single bond are blue; one double bond are red; two single bonds are green;
#    one single one double are magenta; two single one double are black
#    Nodes with a self-loop are drawn in blue; others in red.
#    Actually draws the graph twice -- once just to show topology, once including lengths in edges
# =============================================================================
    GGraph = nx.Graph(G)
    GGraphEdgesColors = []
    nodesWithSelfLoops = []
    for edge in GGraph.edges(): #with parentheses after edges, this doesn't include the key (len(edge)==2)
        numSBonds = 0
        numDBonds = 0
        for key in G[edge[0]][edge[1]].keys():
            if G[edge[0]][edge[1]][key]['singleBond']:
                numSBonds += 1
            else:
                numDBonds += 1
        if numSBonds == 1 and numDBonds == 0:
            GGraphEdgesColors.append('b')
        elif numSBonds == 2 and numDBonds == 0:
            GGraphEdgesColors.append('g')
        elif numSBonds == 0 and numDBonds == 1:
            GGraphEdgesColors.append('r')
        elif numSBonds == 1 and numDBonds == 1:
            GGraphEdgesColors.append('m')
        elif numSBonds == 2 and numDBonds == 1:
            GGraphEdgesColors.append('k')
        
        if edge[0] == edge[1]:
            nodesWithSelfLoops.append(edge[0])
    GGraphNodeColors = ['xkcd:blue' if node in nodesWithSelfLoops else 'r' for node in GGraph.nodes()]                        
    
    pos = nx.spring_layout(GGraph)
    
    if drawTopology:
        plt.figure(figsize=figsize)  
        nx.draw_networkx(GGraph, pos, edge_color=GGraphEdgesColors, node_color = GGraphNodeColors, alpha = 0.7)
        plt.axis('off')
        plt.show()

    if drawEdgeLengths:
        plt.figure(figsize=figsize)    
#        nx.draw(GGraph,pos,edge_color='black',width=1,linewidths=1,\
#        node_size=500,node_color='white',alpha=0.3,\
#        labels={node:node for node in GGraph.nodes()})
        nx.draw_networkx(GGraph, pos, edge_color=GGraphEdgesColors, node_color = GGraphNodeColors, alpha = 0.2)
        elabels = dict()
        for edge in G.edges:
            elabels[(edge[0], edge[1])] = ''
        for edge in G.edges: #without parentheses after edges, this includes the key as edge[2]
            if G[edge[0]][edge[1]][edge[2]]['singleBond']:
                elabels[(edge[0], edge[1])] += 's' + str(G[edge[0]][edge[1]][edge[2]]['length']) + ' '
            else:
                elabels[(edge[0], edge[1])] += 'd' + str(G[edge[0]][edge[1]][edge[2]]['length']) + ' '
        nx.draw_networkx_edge_labels(GGraph,pos,edge_labels=elabels,font_color='black',font_size=10)
        plt.axis('off')
        plt.show()


def drawRNAStructure(structure,sequence,linkerPos = [], figsize = (7,7), showFig = True):
# =============================================================================
    
# =============================================================================
    f = plt.figure(figsize=figsize)  
    numNt = len(sequence)
    G = nx.Graph()
    nodeLabels = {}
    
    for i in range(numNt): #make a node for all nts (including those in linker)
#        We'll remove linker nodes later, but otherwise it makes nodelist confusing
        G.add_node(sequence[i] + str(i))
        nodeLabels[sequence[i] + str(i)] = sequence[i]  + str(i)
    
    edgeColors = []
    nodeColors = ['xkcd:salmon']*numNt
    nodelist = list(G.nodes)
    
#    Make the backbone edges
    for i in range(numNt - 1):
        if i not in linkerPos and i+1 not in linkerPos:
            G.add_edge(nodelist[i],nodelist[i+1], backbone = True, weight = 1.5) #higher weight means more distance between bonded nodes
    
#    Make the bond edges
    for stem in structure:
        numBonds = int(len(stem)/2)
        for i in range(numBonds):
            G.add_edge(nodelist[stem[i]], nodelist[stem[i+numBonds]], weight = 1, backbone = False)
            nodeColors[stem[i]] = 'xkcd:sky blue'; nodeColors[stem[i + numBonds]] = 'xkcd:sky blue'
    
    nodeColors = [nodeColors[i] for i in range(numNt) if i not in linkerPos]
    for i in linkerPos:
        G.remove_node(sequence[i] + str(i))
        del(nodeLabels[sequence[i] + str(i)])
    
    for edge in G.edges:
        if int(edge[1][1:]) - int(edge[0][1:]) == 1:
            edgeColors.append('black')
        else:
            edgeColors.append('xkcd:sky blue')
        
    pos = nx.kamada_kawai_layout(G,weight = 'weight')
    nx.draw_networkx(G, pos, edge_color=edgeColors, node_color = nodeColors, labels=nodeLabels, alpha = 0.7, font_size=9)
    plt.axis('off')
    if showFig:
        plt.show()
    return(f)
#    plt.show()

def createRefNets():
# =============================================================================
#     Create the list of reference graphs (nets) for which we've calculated the entropy
#    (this is Fig. S2 of the supplement, plus the open-open-net-2).
#    If we choose to define stems of length 0 as single nodes, we need to add here
#    openNet2Single, closedNet2Single, openNet2cL1Zero, closedNet2cL1Zero
# =============================================================================
    openNet0 = nx.MultiGraph()
    openNet0.add_edge(0, 1, singleBond=True, length=1)
    
    closedNet0 = nx.MultiGraph()
    closedNet0.add_edge(0, 0, singleBond=True, length=1)
    
    openNet1 = nx.MultiGraph()
    openNet1.add_edge(0, 1, singleBond=False, length=1)
    openNet1.add_edge(0, 1, singleBond=True, length=1)
    
    closedNet1 = openNet1.copy()
    closedNet1.add_edge(0, 1, singleBond=True, length=1)
    
    openNet2a = nx.MultiGraph()
    openNet2a.add_edge(0, 1, singleBond=False, length=1)
    openNet2a.add_edge(0, 2, singleBond=True, length=1)
    openNet2a.add_edge(1, 2, singleBond=True, length=1)
    openNet2a.add_edge(1, 3, singleBond=True, length=1)
    openNet2a.add_edge(2, 3, singleBond=False, length=1)
    
    closedNet2a = openNet2a.copy()
    closedNet2a.add_edge(0, 3, singleBond=True, length=1)
    
    openNet2b = nx.MultiGraph()
    openNet2b.add_edge(0, 1, singleBond = False, length = 1)
    openNet2b.add_edge(0, 2, singleBond = True, length = 1)
    openNet2b.add_edge(1, 3, singleBond = True, length = 1)
    openNet2b.add_edge(2, 3, singleBond = True, length = 1)
    openNet2b.add_edge(2, 3, singleBond = False, length = 1)
    
    closedNet2b = openNet2b.copy()
    closedNet2b.add_edge(0, 1, singleBond = True, length = 1)
    
    openNet2c = nx.MultiGraph()
    openNet2c.add_edge(0, 1, singleBond = False, length = 1)
    openNet2c.add_edge(0, 2, singleBond = True, length = 1)
    openNet2c.add_edge(0, 2, singleBond = True, length = 1)
    openNet2c.add_edge(1, 3, singleBond = True, length = 1)
    openNet2c.add_edge(2, 3, singleBond = False, length = 1)
    
    closedNet2c = openNet2c.copy()
    closedNet2c.add_edge(1, 3, singleBond = True, length = 1)
    
    openOpenNet2 = nx.MultiGraph()
    openOpenNet2.add_edge(0, 1, singleBond=False, length=1)
    openOpenNet2.add_edge(0, 2, singleBond=True, length=1)
    openOpenNet2.add_edge(1, 3, singleBond=True, length=1)
    openOpenNet2.add_edge(2, 3, singleBond=False, length=1)
    
    refNets = dict(openNet0 = openNet0, closedNet0 = closedNet0, 
                   openNet1 = openNet1, closedNet1 = closedNet1, 
                   openNet2a = openNet2a, closedNet2a = closedNet2a,
                   openNet2b = openNet2b, closedNet2b = closedNet2b, 
                   openNet2c = openNet2c, closedNet2c = closedNet2c, 
                   openOpenNet2 = openOpenNet2)
    
    return(refNets)
    
    
def createGraphFromStructure(structure, numNt, linkerPos): 
    
    G = nx.MultiGraph()
    numStems = len(structure)
    numVs = 0 #each double bond comes with a factor of v_s in the entropy calculation
    
    if not numStems:
        return(G, numVs, [])
# =============================================================================
#   first set up nodes and double bonded edges. It is possible do do this at the step of creating the
#   stems (STable), and then to use the networkx functions union_all, compose_all, or disjoint_union_all
# =============================================================================
    nodeCounter = 0
    edgeCounter = 0
##    keep track of the nodes that start/end the stem. If the base pair that defines
##    the node is [i,j], the former have a single bond connection to: 1) the node for which one of the bps is
##    i-x for the smallest possible x, and 2) the node for which one of the bps is j+y for the smallest possible y
##    The latter have connections to i+z for the smallest z, and j-w for the smallest w
    
#    make a dictionary of the base pairs which define the nodes where for each bp, we save the following:
#    [corresponding node, need to find previous nt to which it's connected by a single bond, 
#        need to find next nt to which it's connected by a single bond]
    
#    The graph decomposition process is actually much easier if we treat 
#    nodes that correspond to stems of length 1 the same way we treat other stems,
#    that is, create two new nodes (in this case they'll have the same bps)
#    and create a double-bond edge between them (in this case it'll have length 0)
    
#    #From when we wanted to treat stems of length 1 as just a single node:
#    bondedNtNodeDict = {}
#    for i in range(numStems):
#        lenStem = int(len(structure[i])/2)
#        if lenStem > 1: #if we're not dealing with a single bp
#            #the stem corresponds to two nodes (for the two edges of the stem) connected by a double bond edge
#            firstNodeBPs = [structure[i][0], structure[i][lenStem]]
#            G.add_node(nodeCounter, bp = firstNodeBPs)
#            bondedNtNodeDict[firstNodeBPs[0]] = [nodeCounter, True, False] 
#            bondedNtNodeDict[firstNodeBPs[1]] = [nodeCounter, False, True] 
#                #[node, single bond connects to previous node, single bond connects to next node]
#            nodeCounter += 1
#            
#            secondNodeBPs = [structure[i][lenStem - 1], structure[i][-1]]
#            G.add_node(nodeCounter, bp = secondNodeBPs)
#            bondedNtNodeDict[secondNodeBPs[0]] = [nodeCounter, False, True] 
#            bondedNtNodeDict[secondNodeBPs[1]] = [nodeCounter, True, False] 
#            
#            G.add_edge(nodeCounter - 1, nodeCounter, singleBond = False, 
#                       length = lenStem - 1, key = edgeCounter) #make the double bond edge
#            nodeCounter += 1
#            edgeCounter += 1
#        else: #If I want to treat stems of length 1 with just one node:
#            firstAndSecondNodeBPs = [structure[i][0], structure[i][lenStem]]
#            G.add_node(nodeCounter,bp = firstAndSecondNodeBPs)
##            firstNodeOfStemList[i] = firstAndSecondNodeBPs
##            secondNodeOfStemList[i] = firstAndSecondNodeBPs
#            for bondedNt in firstAndSecondNodeBPs:
#                bondedNtNodeDict[bondedNt] = [nodeCounter, True, True]
#            nodeCounter += 1
    
    #Create the adjacency matrix for the graph to make it fast to find bridges --
#    edges which when you remove them disconnect the graph.
    M = np.zeros((100,100),dtype=int) #assume we'll have fewer than 100 nodes
    edgeList = []
    
    bondedNtNodeDict = {}
    for i in range(numStems):
        lenStem = int(len(structure[i])/2)
        
        firstNodeBPs = [structure[i][0], structure[i][lenStem]]
        secondNodeBPs = [structure[i][lenStem - 1], structure[i][-1]]
        #if lenStem == 1, firstNodeBPs and secondNodeBPs are identical
        for j in firstNodeBPs + secondNodeBPs: #initialize dictionary entries
#            This is useful because if lenStem == 1, the dictionary keys are identical 
#            and so we want to use appending functions here rather than equalities
            bondedNtNodeDict[j] = []
            
        G.add_node(nodeCounter, bp = firstNodeBPs)
        bondedNtNodeDict[firstNodeBPs[0]].append([nodeCounter, True, False])
        bondedNtNodeDict[firstNodeBPs[1]].append([nodeCounter, False, True])
            #[node, single bond connects to previous node, single bond connects to next node]
        M[nodeCounter,nodeCounter] = 1
        nodeCounter += 1
        
        G.add_node(nodeCounter, bp = secondNodeBPs)
        bondedNtNodeDict[secondNodeBPs[0]].append([nodeCounter, False, True])
        bondedNtNodeDict[secondNodeBPs[1]].append([nodeCounter, True, False])
        M[nodeCounter,nodeCounter] = 1
        
        G.add_edge(nodeCounter - 1, nodeCounter, singleBond = False, 
                   length = lenStem - 1)#, key = edgeCounter) #make the double bond edge
        M[nodeCounter - 1, nodeCounter] += 1; M[nodeCounter, nodeCounter - 1] += 1; 
        edgeList.append((nodeCounter - 1,nodeCounter))
        nodeCounter += 1
        edgeCounter += 1
        numVs += 1
    
# =============================================================================
#     #set up nodes for the first and last nts of each sequence (if they haven't yet been 
#     #put in nodes). If we're just using the graph for the entropy calculation, 
#     #this in principle isn't necessary since whatever (single) bonds these nodes
#     #will be connected to will be immediately discarded and won't affect the entropy
#     #However, it makes it easier since including this makes it so we don't need to 
#     #explicitly check when we make a single bond edge if it corresponds to the linker
# =============================================================================
    boundNts = structure2boundNt(structure, numNt)
    if boundNts[0] != 0: #then nt zero hasn't yet been put into a node
        boundNt = 0
        G.add_node(nodeCounter, bp = [boundNt, boundNt])
        bondedNtNodeDict[boundNt] = [[nodeCounter, False, True]]
        M[nodeCounter,nodeCounter] = 1
        nodeCounter += 1
    else: #if it has been put into a node, correct it so that we don't look for previous nt
        boundNt = 0
        for dictEntry in bondedNtNodeDict[boundNt]:
            dictEntry[1] = False
        
    if boundNts[-1] != numNt - 1: #then the last nt hasn't been put into a node
        boundNt = numNt - 1
        G.add_node(nodeCounter, bp = [boundNt, boundNt])
        bondedNtNodeDict[boundNt] = [[nodeCounter, True, False]]
        M[nodeCounter,nodeCounter] = 1
        nodeCounter += 1
    else: #if it has been put into a node, correct it so that we don't look for next nt
        boundNt = numNt - 1
        for dictEntry in bondedNtNodeDict[boundNt]:
            dictEntry[2] = False
        
    if len(linkerPos): #then we also need to consider the last nt of first strand and first nt of second strand
        if linkerPos[-1] + 1 not in boundNts: #first nt of second strand
            boundNt = linkerPos[-1] + 1
            G.add_node(nodeCounter, bp = [boundNt, boundNt])
            bondedNtNodeDict[boundNt] = [[nodeCounter, False, True]]
            M[nodeCounter,nodeCounter] = 1
            nodeCounter += 1
        else: #if it has been put into a node, correct it so that we don't look for previous nt
            boundNt = linkerPos[-1] + 1
            for dictEntry in bondedNtNodeDict[boundNt]:
                dictEntry[1] = False
            
        if linkerPos[0] - 1 not in boundNts: #last nt of first strand
            boundNt = linkerPos[0] - 1
            G.add_node(nodeCounter, bp = [boundNt, boundNt])
            bondedNtNodeDict[boundNt] = [[nodeCounter, True, False]] 
            M[nodeCounter,nodeCounter] = 1
            nodeCounter += 1
        else: #if it has been put into a node, correct it so that we don't look for next nt
            boundNt = linkerPos[0] - 1
            for dictEntry in bondedNtNodeDict[boundNt]:
                dictEntry[2] = False
            
            
# =============================================================================
#   now that we've created the nodes and double-bond edges, create the single-bond edges.
#   Before making a single-bond edge, make sure that it doesn't correspond to the linker.
#   This step is unnecessary if we include the first and last nts of each strand as nodes
# =============================================================================
    for bondedNt in bondedNtNodeDict.keys():
        for dictEntry in bondedNtNodeDict[bondedNt]: #added when we started treating stems of length 1 as two nodes
            node = dictEntry[0] #node we're considering
            
            if dictEntry[1]: #need to find previous bondedNt to which it's connected by a single stem
                previousBondedNts = [i for i in bondedNtNodeDict.keys() if i < bondedNt]
                if previousBondedNts: #meaning there is a single bond to be made between the nodes.
    #                This could be empty if bondedNt is the first nt of a strand or 
    #                is connected to the first nt (which hasn't yet been made into a node) 
                    connectingNt = max(previousBondedNts)
                    if len(bondedNtNodeDict[connectingNt]) == 1: #there's only one entry for this node
                        connectingDictEntry = bondedNtNodeDict[connectingNt][0]
                    elif bondedNtNodeDict[connectingNt][0][2]:
                        connectingDictEntry = bondedNtNodeDict[connectingNt][0]
                    elif bondedNtNodeDict[connectingNt][1][2]:
                        connectingDictEntry = bondedNtNodeDict[connectingNt][1]
                    else:
                        print('Weird! We weren''t supposed to make the edge between nts ' + str(bondedNt) + ', ' + str(connectingNt))
                    
                    connectingNode = connectingDictEntry[0]
                        
                    #create the edge
                    G.add_edge(node, connectingNode, singleBond = True, 
                           length = bondedNt - connectingNt)#, key = edgeCounter)
                    M[node,connectingNode] += 1; M[connectingNode,node] += 1
                    edgeList.append((node,connectingNode))
                    edgeCounter += 1
                    
                    connectingDictEntry[2] = False #don't want to make the edge again when we get to connectingNt
                    dictEntry[1] = False #since we just found its connecting single edge
                    
                else: #The next few lines are really just for checking purposes. Can just erase them and make the previous line unindented by one
                    if bondedNt == 0: 
                        dictEntry[1] = False #the first nt won't be connected
                    elif len(linkerPos):
                        if bondedNt == linkerPos[-1] + 1: #then it's the first nt of the second sequence
                            dictEntry[1] = False
                    if dictEntry[1]:
                        print('Could not find the previous bonded nt for nt ' + str(bondedNt))
    
                    
            if dictEntry[2]: #need to find next bondedNt to which it's connected by a single stem
                nextBondedNts = [i for i in bondedNtNodeDict.keys() if i > bondedNt]
                if nextBondedNts: #meaning there is a single bond to be made between the nodes
    #                This could be empty if bondedNt is the last nt of a strand or 
    #                is connected to the last nt (if we don't make it into a node) 
                    connectingNt = min(nextBondedNts)
                    if len(bondedNtNodeDict[connectingNt]) == 1: #there's only one entry for this node
                        connectingDictEntry = bondedNtNodeDict[connectingNt][0]
                    elif bondedNtNodeDict[connectingNt][0][1]:
                        connectingDictEntry = bondedNtNodeDict[connectingNt][0]
                    elif bondedNtNodeDict[connectingNt][1][1]:
                        connectingDictEntry = bondedNtNodeDict[connectingNt][1]
                    else:
                        print('Weird (2)! We weren''t supposed to make the edge between nts ' + str(bondedNt) + ', ' + str(connectingNt))
                    
                    connectingNode = connectingDictEntry[0]
                    
                    #create the edge
                    G.add_edge(node, connectingNode, singleBond = True, 
                           length = connectingNt - bondedNt)#, key = edgeCounter)
                    M[node, connectingNode] += 1; M[connectingNode, node] += 1
                    edgeList.append((node,connectingNode))
                    edgeCounter += 1
                    
                    connectingDictEntry[1] = False #don't want to make the edge again when we get to connectingNt
                    dictEntry[2] = False #since we just found its connecting single stem
                
                else: #The next few lines are really just for checking purposes. Can just erase them and make the previous line unindented by one
                    if bondedNt == numNt - 1: 
                        dictEntry[2] = False #the last nt won't be connected
                    elif len(linkerPos):
                        if bondedNt == linkerPos[0] - 1: #then it's the last nt of the first sequence
                            dictEntry[2] = False
                    if dictEntry[2]:
                        print('Could not find the next bonded nt for nt ' + str(bondedNt))
        
    #The graph should now be properly created. We can check that we didn't miss single edges:
    for bondedNt in bondedNtNodeDict.keys():
        for dictEntry in bondedNtNodeDict[bondedNt]:
            if dictEntry[1]:
                print('Never found the preceding nt to make single edge for nt ' + str(bondedNt))
            if dictEntry[2]:
                print('Never found the next nt to make single edge for nt ' + str(bondedNt))
    
# =============================================================================
#     Find the bridges
# =============================================================================
    edgeList = list(G.edges)
    edgeList = np.array([(edge[0], edge[1]) for edge in edgeList if M[edge[0]][edge[1]] == 1], dtype = int)
    #don't want to try to remove edges that are duplicated (e.g. two edges between the same two nodes for internal loop)
    whichEdgesToRemove = myNumConnCompFromMatC(M[:nodeCounter,:nodeCounter].astype(bool), edgeList)
    bridgeList = [(edgeList[edge][0], edgeList[edge][1]) for edge in whichEdgesToRemove if edge > -1]
    
    return(G, numVs, bridgeList)
       
#    G.add_edge(1,2,singleBond=False,length=2,key=2) #how graph edge addition should look. 
#        #each different edge between the same two nodes must have a different key. It's unclear 
#        #at this point if creating the keys will be important (I think it won't be)
#    G.number_of_edges(1,2) #returns number of edges between nodes 1 and 2
    
def graphDecomposition(G, componentGraphs=[], bridgeList = []):
# =============================================================================
#     Perform the graph decomposition on a graph G. Keep track of component graphs
#     found (and removed) during the graph decomposition process
#    
#     The graph decomposition process proceeds as follows:
#     We check each bond if removing it will disconnect the graph
#     into more components. If so, we remove it (logging the fact that we did so)
#     We then check each node if it is not connected to anything (if so we remove it)
#     or if it is connected only to two single bonds, in which case we can concatenate those two bonds
#     and remove the node (i.e. we remove the node and connect the nodes to which it was connected to each other).
#    
#    If we were to define a stem of length one as a single node, we would need to also do the following:
#    For practical reasons, it's useful to also check if a node is connected to itself, and if so,
#    to create a new node connected to itself and delete that one. Otherwise, if the node that's connected
#    to itself comes in the middle of a single-bond region of a larger structure (say a pseudoknot), it's tricky
#    to tell that it needs to be disconnected from that structure. This is only an issue if minBPInStem == 1.
# =============================================================================

#    Gstart = G.copy()
    if not len(componentGraphs): #initialize componentGraphs if it isn't inputted
        componentGraphs = {}     
        netTypes = ['openNet0','closedNet0','openNet1','closedNet1','openOpenNet2',
        'openNet2a','closedNet2a','openNet2b','closedNet2b','openNet2c','closedNet2c','tooComplex']
        for net in netTypes:
            componentGraphs[net] = []
    
    #completedGraphDecomposition = False
            
#        #This is slower:
#    nConnectedComponents = nx.number_connected_components(G) #will be 1 here at the start
#    edgeList = list(G.edges) #need this so it doesn't change within the for loop
#    for edge in edgeList:
#        if len(G[edge[0]][edge[1]]) == 1: #because we don't want to get rid of edges with two bonds
#            H = nx.Graph(G) #Things are faster for graphs than multigraphs, and for checking connectedness
##            a graph is fine (as long as we previously check that we're not getting rid of edges with two bonds)
#            H.remove_edge(edge[0], edge[1])
#            if nx.number_connected_components(H) > nConnectedComponents:
#                if G.edges[edge[0], edge[1], edge[2]]['singleBond']: #key must be zero since there can't be another edge here
#                    #keep track of single edges deleted
#                    componentGraphs['openNet0'].append([G.edges[edge[0], edge[1], edge[2]]['length'], -1, -1, -1, -1, -1])
#                G.remove_edge(edge[0], edge[1], edge[2])
#                nConnectedComponents = nx.number_connected_components(H)
    
#    graphDecompCounter = 0
    
#    Having bridgeList computed before using our built-in code (that can be translated to C
#       to make it faster) makes the code take basically exactly the same amount of time
#       as using nx.bridges() (it appears ever so slightly faster).
    for edge in bridgeList:
        if G.edges[edge[0], edge[1], 0]['singleBond']: #key must be zero since there can't be another edge here
            #keep track of single edges deleted
            componentGraphs['openNet0'].append([G.edges[edge[0], edge[1], 0]['length'], -1, -1, -1, -1, -1])
        G.remove_edge(edge[0], edge[1])
    
#        if graphDecompCounter == 3:
#            print('need to put graph decomposition in loop for G = ')
#            print('numBonds adjacency matrix = ')
#            print(nx.to_numpy_matrix(Gstart))
#            print('singleBonds adjacency matrix = ')
#            print(nx.to_numpy_matrix(Gstart, weight='singleBond'))
        
#        bridgeList = nx.bridges(nx.Graph(G)) #faster than checking each edge individually
##        Making bridgeList is ~10x faster than the minimum_edge_cut function
##        However, bridgeList is very slow compared to the rest of the code (definitely rate-limiting).
##        It's slow to make the generator, and it's slow to use the generator to make the edges (removing edges is fast)
#        for edge in bridgeList:
#            if len(G[edge[0]][edge[1]]) == 1: #because we had to convert the multigraph to a regular graph
#    #                and we don't want to get rid of edges with two bonds
#                if G.edges[edge[0], edge[1], 0]['singleBond']: #key must be zero since there can't be another edge here
#                    #keep track of single edges deleted
#                    componentGraphs['openNet0'].append([G.edges[edge[0], edge[1], 0]['length'], -1, -1, -1, -1, -1])
#                G.remove_edge(edge[0], edge[1])
#                completedGraphDecomposition = False

    nodeList = list(G.nodes)
        
        
#        This is also slower than bridgeList (by about 2x). Besides myNumConnCompFromMatC, making M is slow
#        M = nx.to_numpy_array(G, dtype = bool) #doesn't appear to matter for speed if it's bool or int or default
#        if len(M):
#            M += np.eye(len(M), dtype=bool) #put ones on diagonals. It's ok if elements are 2, etc.
#            edgeList = list(G.edges)
#            edgeList = np.array([(nodeList.index(edge[0]), nodeList.index(edge[1])) for edge in edgeList
#                                  if len(G[edge[0]][edge[1]]) == 1 and not (edge[0] == edge[1])], dtype = int)
#            
#            whichEdgesToRemove = myNumConnCompFromMatC(M,edgeList)
#            
#            edgesToRemove = [(nodeList[edgeList[edge][0]], nodeList[edgeList[edge][1]]) for edge in whichEdgesToRemove if edge > -1]
#            for edge in edgesToRemove:
#                if G.edges[edge[0], edge[1], 0]['singleBond']: #key must be zero since there can't be another edge here
#                    #keep track of single edges deleted
#                    componentGraphs['openNet0'].append([G.edges[edge[0], edge[1], 0]['length'],
#                                                       -1, -1, -1, -1, -1])
#                completedGraphDecomposition = False
#            G.remove_edges_from(edgesToRemove)
        
# This is also slower than bridgeList
#        M = nx.to_numpy_array(G, dtype = bool) #doesn't appear to matter for speed if it's bool or int or default
#        M += np.eye(len(M), dtype=bool) #put ones on diagonals. It's ok if elements are 2, etc.
#        nConnectedComponents = nx.number_connected_components(G) #will be 1 here at the start
#        edgeList = list(G.edges) #need this so it doesn't change within the for loop
        
#        edgesToRemove = []
#        for edge in edgeList:
#            if len(G[edge[0]][edge[1]]) == 1 and not (edge[0] == edge[1]):
#    #                since we don't want to get rid of edges with two bonds or self-loops
#                M2 = copy.copy(M)
#                mEdge = (nodeList.index(edge[0]), nodeList.index(edge[1]))
#                M2[mEdge[0], mEdge[1]] = 0; M2[mEdge[1], mEdge[0]] = 0
#                M2NCC = myNumConnCompFromMat(M2)
#                
#                if M2NCC > nConnectedComponents:
#                    if G.edges[edge[0], edge[1], 0]['singleBond']: #key must be zero since there can't be another edge here
#                        #keep track of single edges deleted
#                        componentGraphs['openNet0'].append([G.edges[edge[0], edge[1], 0]['length'], -1, -1, -1, -1, -1])
#                    edgesToRemove.append(edge)
#                    completedGraphDecomposition = False
#                    nConnectedComponents = M2NCC
#                    M = M2
#        G.remove_edges_from(edgesToRemove)
        
    for node in nodeList:
        if G.degree(node) == 0:
            G.remove_node(node)
            #completedGraphDecomposition = False
#        if node is connected only to two single bonds, delete it and connect the 
#        two nodes it's connected to to one another
        elif G.degree(node) == 2 and G.degree(node, weight='singleBond') == 2:
            nodeNeigh = list(G.neighbors(node))
            if not nodeNeigh[0] == node: #if it's a hairpin loop we don't want to do anything
#                The following will work both if node is connected to two different nodes or twice to the same node
                lengths = []
                for neigh in nodeNeigh:
                    for key in G.get_edge_data(node,neigh).keys():
                        lengths.append(G.get_edge_data(node,neigh)[key]['length'])
                G.add_edge(nodeNeigh[0],nodeNeigh[-1], singleBond = True, length = sum(lengths))#, key = edgeCounter)
#                    edgeCounter += 1
                G.remove_node(node)
#                completedGraphDecomposition = False
    
    #    if not G.number_of_nodes():
    #        break
    
    return(G, componentGraphs)
    

def calculateEntropyFromGraph(G, refNets, bridgeList):
# =============================================================================
#     Given a graph representing an RNA structure, return its entropy. Consider
#    the contribution of vs to the entropy separately (simply keep track of the 
#    number of factors of vs needed to be included, to make it easy to change vs 
#    without rerunning this).
# =============================================================================
    em = iso.numerical_multiedge_match('singleBond',0) 

    componentGraphs = {} #which of the possible subgraphs do we get, and what parameters for each?
# =============================================================================
#     componentGraphs is a dictionary, where each entry of the dictionary yields a vector of vectors
#     where each element of the vector is a different parameter set for a different instance of that net.
#     The parameter sets are defined as vectors in the order [s1, s2, s3, s4, l1, l2].
#     If the parameter is not defined for the graph its value is set to be -1 (like s4 for an open-net-2a)
#     Example usage: componentGraphs['closedNet2a'][0][3] gives the value of s4 for the first
#     instance of closedNet2a found in the structure. componentGraphs has these keys: 
#    ['openNet0','closedNet0','openNet1','closedNet1','openOpenNet2','openNet2a',
#    'closedNet2a','openNet2b','closedNet2b','openNet2c','closedNet2c','tooComplex']
# =============================================================================
    
    netTypes = ['openNet0','closedNet0','openNet1','closedNet1','openOpenNet2',
    'openNet2a','closedNet2a','openNet2b','closedNet2b','openNet2c','closedNet2c','tooComplex']
        #'closedNets0_2', need to also be accounted for if we consider stems of length zero as a single node
        #We decided here not to consider 'openNet2Single','closedNet2Single','openNet2cL1Zero','closedNet2cL1Zero'
        #as separate types of graphs (these all are created when at least one stem has length 0)
        
    for net in netTypes:
        componentGraphs[net] = []
            
    #Perform the graph decomposition
    G, componentGraphs = graphDecomposition(G, componentGraphs, bridgeList)
    
    #Find for each subgraph what net it's isomorphic to
    for c in nx.connected_components(G):
        Gsub = G.subgraph(c).copy() #create subgraph of connected components
        
        numNodesGsub = Gsub.number_of_nodes()
        numEdgesGsub = Gsub.number_of_edges()
        
        netTooComplex = False
        couldntMatchNet = False
        foundIsomorphism = False
        
        if numNodesGsub > 4: #then we know none of our evaluated nets will be isomorphic to gSub
            netTooComplex = True
            foundIsomorphism = True
        
        while not foundIsomorphism: #use a while loop so we can break out of it quickly
#            We expect is_isomorphic to be slow (there are faster heuristics we may be able to use)
#            and so start by checking if gSub even has the right number of edges, nodes before 
#            checking isomorphism.
#            #We loop over potential keys since we don't know for sure what the key will be, but really,
#                that's not quite true -- if we let the keys be automatically assigned, it should be zero
#                The only way it wouldn't be zero is if there were already another edge between these two nodes
#                in which case which key was zero and which was one would be arbitrary
            
            if numEdgesGsub == 1 and numNodesGsub == 1:
                #closedNet0
                GM = nx.isomorphism.GraphMatcher(refNets['closedNet0'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('closedNet0')
                    for key in Gsub[GM.mapping[0]][GM.mapping[0]].keys():
                        s1 = Gsub[GM.mapping[0]][GM.mapping[0]][key]['length']
                    
                    componentGraphs['closedNet0'].append([s1, -1, -1, -1, -1, -1])
                    break
            
            elif numEdgesGsub == 2 and numNodesGsub == 2:
                #openNet1
                GM = nx.isomorphism.GraphMatcher(refNets['openNet1'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('openNet1')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        if Gsub[GM.mapping[0]][GM.mapping[1]][key]['singleBond']:
                            s1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                        else:
                            l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    
                    componentGraphs['openNet1'].append([s1, -1, -1, -1, l1, -1])
                    break
                
            elif numEdgesGsub == 3 and numNodesGsub == 2:
                #closedNet1
                GM = nx.isomorphism.GraphMatcher(refNets['closedNet1'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('closedNet1')
                    s12 = [] #doesn't matter which single bond is s1 and which s2
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        if Gsub[GM.mapping[0]][GM.mapping[1]][key]['singleBond']:
                            s12.append(Gsub[GM.mapping[0]][GM.mapping[1]][key]['length'])
                        else:
                            l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    
                    componentGraphs['closedNet1'].append([s12[0], s12[1], -1, -1, l1, -1])
                    break
            
            elif numEdgesGsub == 5 and numNodesGsub == 4:
                #openNet2a
                GM = nx.isomorphism.GraphMatcher(refNets['openNet2a'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('openNet2a')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s3 = Gsub[GM.mapping[0]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[2]].keys():
                        s2 = Gsub[GM.mapping[1]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s1 = Gsub[GM.mapping[1]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                    
                    componentGraphs['openNet2a'].append([s1, s2, s3, -1, l1, l2])
                    break
                
                #openNet2b
                GM = nx.isomorphism.GraphMatcher(refNets['openNet2b'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('openNet2b')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s1 = Gsub[GM.mapping[0]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s3 = Gsub[GM.mapping[1]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        if Gsub[GM.mapping[2]][GM.mapping[3]][key]['singleBond']:
                            s2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                        else:
                            l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                    componentGraphs['openNet2b'].append([s1, s2, s3, -1, l1, l2])
                    break
                
                #openNet2c
                GM = nx.isomorphism.GraphMatcher(refNets['openNet2c'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('openNet2c')
                    s12 = [] #doesn't matter which single bond is s1 and which s2
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s12.append(Gsub[GM.mapping[0]][GM.mapping[2]][key]['length'])
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s3 = Gsub[GM.mapping[1]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                    
                    componentGraphs['openNet2c'].append([s12[0], s12[1], s3, -1, l1, l2])
                    break
                
            elif numEdgesGsub == 6 and numNodesGsub == 4:
                #closedNet2a
                GM = nx.isomorphism.GraphMatcher(refNets['closedNet2a'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('closedNet2a')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s3 = Gsub[GM.mapping[0]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[3]].keys():
                        s4 = Gsub[GM.mapping[0]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[2]].keys():
                        s2 = Gsub[GM.mapping[1]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s1 = Gsub[GM.mapping[1]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                    
                    componentGraphs['closedNet2a'].append([s1, s2, s3, s4, l1, l2])
                    break
                
                #closedNet2b
                GM = nx.isomorphism.GraphMatcher(refNets['closedNet2b'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('closedNet2b')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        if Gsub[GM.mapping[0]][GM.mapping[1]][key]['singleBond']:
                            s4 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                        else:
                            l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s1 = Gsub[GM.mapping[0]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s3 = Gsub[GM.mapping[1]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        if Gsub[GM.mapping[2]][GM.mapping[3]][key]['singleBond']:
                            s2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                        else:
                            l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']
                    
                    componentGraphs['closedNet2b'].append([s1, s2, s3, s4, l1, l2])
                    break
                
                #closedNet2c
                GM = nx.isomorphism.GraphMatcher(refNets['closedNet2c'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('closedNet2c')
                    s12 = [] #doesn't matter which single bond is s1 and which s2
                    s34 = []
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s12.append(Gsub[GM.mapping[0]][GM.mapping[2]][key]['length'])
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s34.append(Gsub[GM.mapping[1]][GM.mapping[3]][key]['length'])
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length'] 
                    
                    componentGraphs['closedNet2c'].append([s12[0], s12[1], s34[0], s34[1], l1, l2])
                    break
                
            elif numEdgesGsub == 4 and numNodesGsub == 4:
                #openOpenNet2 (also known as the "very open net 2" in our paper)
                GM = nx.isomorphism.GraphMatcher(refNets['openOpenNet2'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('openOpenNet2')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys():
                        l1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    for key in Gsub[GM.mapping[0]][GM.mapping[2]].keys():
                        s1 = Gsub[GM.mapping[0]][GM.mapping[2]][key]['length']
                    for key in Gsub[GM.mapping[1]][GM.mapping[3]].keys():
                        s2 = Gsub[GM.mapping[1]][GM.mapping[3]][key]['length']
                    for key in Gsub[GM.mapping[2]][GM.mapping[3]].keys():
                        l2 = Gsub[GM.mapping[2]][GM.mapping[3]][key]['length']    
                    
                    componentGraphs['openOpenNet2'].append([s1, s2, -1, -1, l1, l2])
                    break
                
            elif numEdgesGsub == 1 and numNodesGsub == 2:
                #openNet0. Should already have been removed by graph decomposition process
                GM = nx.isomorphism.GraphMatcher(refNets['openNet0'],Gsub,edge_match=em)
                if GM.is_isomorphic():
                    #print('openNet0')
                    for key in Gsub[GM.mapping[0]][GM.mapping[1]].keys(): #nodes connected by edge
                        s1 = Gsub[GM.mapping[0]][GM.mapping[1]][key]['length']
                    
                    componentGraphs['openNet0'].append([s1, -1, -1, -1, -1, -1])
                    break
                
            print('couldn''t match net')
            print('numBonds adjacency matrix = ')
            print(nx.to_numpy_matrix(Gsub))
            print('singleBonds adjacency matrix = ')
            print(nx.to_numpy_matrix(Gsub, weight='singleBond'))
            couldntMatchNet = True
            foundIsomorphism = True
        
        if couldntMatchNet or netTooComplex:
            componentGraphs['tooComplex'].append([-1,-1,-1,-1,-1,-1])

    return(componentGraphs)


# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================

# =============================================================================
# Entropy calculation given graph (formale found in Supplementary Fig. S2), and limits thereof
# =============================================================================    
def openNet0Entropy(s1, g):
    return(0) #or np.log(1)

def closedNet0Entropy(s1, g):
    expEnt = (g / (np.pi * s1))**(3/2)
    return(np.log(expEnt))

def openNet2SingleEntropy(s1, s2, s3, g):
    D = s1*s2 + s1*s3 + s2*s3
    expEnt = g**3 / (np.pi**3 * D**(3/2))
    return(np.log(expEnt))

def closedNet2SingleEntropy(s1, s2, s3, s4, g):
    D = s1*s2*s3 + s1*s2*s4 + s1*s3*s4 + s2*s3*s4
    expEnt = g**(9/2) / (np.pi**(9/2) * D**(3/2))
    return(np.log(expEnt))

def openNet1Entropy(s1, l1, g):
    if l1 == 0:
        return(closedNet0Entropy(s1 = s1, g = g))
    
    expEnt = (g / (np.pi * s1))**(3/2) * np.exp(-g * l1**2 / s1)
    return(np.log(expEnt))

def closedNet1Entropy(s1, s2, l1, g):
    if l1 == 0:
        return(closedNet0Entropy(s1 = s1, g = g) + closedNet0Entropy(s1 = s2, g = g))
    expEnt = (((g / (np.pi * s1)) * (g / (np.pi * s2)))**(3/2) * 
              np.exp(-g * (s1 + s2) * l1**2 / (s1 * s2)))
    return(np.log(expEnt))

def openNet2cL1ZeroEntropy(s1, s2, s3, l1, g):
#    From the MatLab, s0->s1, s1->s2, s2->s3
    if l1 == 0:
        return(openNet2SingleEntropy(s1 = s1, s2 = s2, s3 = s3, g = g))
    
    D = s1*s2 + s1*s3 + s2*s3
    expEnt = g**3 / (np.pi**3 * D**(3/2)) * np.exp(-g * l1**2 * (s1 + s2) / D)
    return(np.log(expEnt))

def closedNet2cL1ZeroEntropy(s1, s2, s3, s4, l1, g):
#    From the MatLab, s0->s1, s1->s2, s2->s3, s3->s4
    if l1 == 0:
        return(closedNet2SingleEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, g = g))
    
    D = s1*s2*s3 + s1*s2*s4 + s1*s3*s4 + s2*s3*s4
    expEnt = g**(9/2) / (np.pi**(9/2) * D**(3/2)) * np.exp(-g * l1**2 * (s1 + s2)*(s3 + s4) / D)
    return(np.log(expEnt))    

def openOpenNet2Entropy(s1, s2, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(closedNet0Entropy(s1 = s1 + s2, g = g))
    elif l1 == 0:
        return(openNet1Entropy(s1 = s1 + s2, l1 = l2, g = g))
    elif l2 == 0:
        return(openNet1Entropy(s1 = s1 + s2, l1 = l1, g = g))
    
    expEnt = (g**(1/2) * np.sinh(2 * g * l1 * l2 / (s1 + s2)) / 
              (2 * np.pi**(3/2) * l1 * l2 * np.sqrt(s1 + s2)) *
              np.exp(-g * (l1**2 + l2**2) / (s1 + s2)))
    return(np.log(expEnt))

def openNet2aEntropy(s1, s2, s3, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(openNet2SingleEntropy(s1 = s1, s2 = s2, s3 = s3, g = g))
    elif l1 == 0: #s1 and s2 should be connected to the same node so s1 and s3 need to be switched here
        return(openNet2cL1ZeroEntropy(s1 = s3, s2 = s2, s3 = s1, l1 = l2, g = g))
    elif l2 == 0: #s1 and s2 should be connected to the same node so it's fine as is
        return(openNet2cL1ZeroEntropy(s1 = s1, s2 = s2, s3 = s3, l1 = l1, g = g))
    
    D = s1*s2 + s1*s3 + s2*s3
    expEnt = (g**2 * np.sinh(2*g*l1*l2*s2 / D) / 
              (2*np.pi**3*l1*l2*s2*np.sqrt(D)) * 
              np.exp(-g*(l1**2 * (s1 + s2) + l2**2 * (s2 + s3)) / D))
    return(np.log(expEnt))

def closedNet2aEntropy(s1, s2, s3, s4, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(closedNet2SingleEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, g = g))
    elif l1 == 0: #s1 and s2 should be equivalent, as should s3 and s4, so we need to switch things around
        return(closedNet2cL1ZeroEntropy(s1 = s3, s2 = s2, s3 = s1, s4 = s4, l1 = l2, g = g))
    elif l2 == 0: #s1 and s2 should be equivalent, as should s3 and s4, so things work out
        return(closedNet2cL1ZeroEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, l1 = l1, g = g))
    
    D = s1*s2*s3 + s1*s2*s4 + s1*s3*s4 + s2*s3*s4
    if s1*s3 == s2*s4: #need to take limit where s1*s3-s2*s4 == 0
        expEnt = (g**(9/2) / (np.pi**(9/2) * D**(3/2)) *
                  np.exp(-g*(l1**2 * (s1+s2)*(s3+s4) + l2**2 * (s1+s4)*(s2+s3)) / D))
    else:
        expEnt = (g**(7/2) * np.sinh(2*g*l1*l2*(s1*s3 - s2*s4) / D) / 
                  (2*np.pi**(9/2) * l1*l2*(s1*s3 - s2*s4)*np.sqrt(D)) * 
                  np.exp(-g*(l1**2 * (s1+s2)*(s3+s4) + l2**2 * (s1+s4)*(s2+s3)) / D))
    return(np.log(expEnt))

def openNet2bEntropy(s1, s2, s3, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(closedNet0Entropy(s1 = s1 + s3, g = g) + closedNet0Entropy(s1 = s2, g = g))
    elif l1 == 0: 
        return(closedNet1Entropy(s1 = s1 + s3, s2 = s2, l1 = l2, g = g))
    elif l2 == 0: 
        return(openNet1Entropy(s1 = s1 + s3, l1 = l1, g = g) + closedNet0Entropy(s1 = s2, g = g))
    
    expEnt = (g**2 * np.sinh(2*g*l1*l2 / (s1 + s3)) / 
              (2 * np.pi**3 * l1*l2 * s2**(3/2) * np.sqrt(s1 + s3)) * 
              np.exp(-g*(l1**2 * s2 + l2**2 * (s1 + s2 + s3)) / (s2 * (s1 + s3))))
    return(np.log(expEnt))

def closedNet2bEntropy(s1, s2, s3, s4, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(closedNet0Entropy(s1 + s3, g = g) + closedNet0Entropy(s2, g = g) + closedNet0Entropy(s4, g = g))
    elif l1 == 0: 
        return(closedNet1Entropy(s1 = s1 + s3, s2 = s2, l1 = l2, g = g) + closedNet0Entropy(s1 = s4, g = g))
    elif l2 == 0: 
        return(closedNet1Entropy(s1 = s1 + s3, s2 = s4, l1 = l1, g = g) + closedNet0Entropy(s1 = s2, g = g))
    
    expEnt = (g**(7/2) * np.sinh(2*g*l1*l2 / (s1 + s3)) / 
              (2 * np.pi**(9/2) * l1*l2 * (s2*s4)**(3/2) * np.sqrt(s1 + s3)) * 
              np.exp(-g*(l1**2 * s2*(s1+s3+s4) + l2**2 * s4*(s1+s2+s3)) / (s2 * s4 * (s1 + s3))))
    return(np.log(expEnt))
    
def openNet2cEntropy(s1, s2, s3, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(openNet2SingleEntropy(s1 = s1, s2 = s2, s3 = s3, g = g))
    elif l1 == 0: #s1 and s2 should be connected to the same node so it's fine as is
        return(openNet2cL1ZeroEntropy(s1 = s1, s2 = s2, s3 = s3, l1 = l2, g = g))
    elif l2 == 0: #s1 and s2 should be connected to the same node so it's fine as is
        return(openNet2cL1ZeroEntropy(s1 = s1, s2 = s2, s3 = s3, l1 = l1, g = g))
    
    D = s1*s2 + s1*s3 + s2*s3
    expEnt = (g**2 * np.sinh(2 * g * l1 * l2 * (s1 + s2) / D) / 
              (2 * np.pi**3 * l1*l2*(s1+s2) * np.sqrt(D)) * 
              np.exp(-g * (l1**2 + l2**2) * (s1 + s2) / D))
    return(np.log(expEnt))
    

def closedNet2cEntropy(s1, s2, s3, s4, l1, l2, g):
    if l1 == 0 and l2 == 0:
        return(closedNet2SingleEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, g = g))
    elif l1 == 0: #s1 and s2 should be equivalent, as should s3 and s4, so things work out
        return(closedNet2cL1ZeroEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, l1 = l2, g = g))
    elif l2 == 0: #s1 and s2 should be equivalent, as should s3 and s4, so things work out
        return(closedNet2cL1ZeroEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, l1 = l1, g = g))
    
    D = s1*s2*s3 + s1*s2*s4 + s1*s3*s4 + s2*s3*s4
    expEnt = (g**(7/2) * np.sinh(2*g*l1*l2*(s1 + s2)*(s3 + s4) / D) / 
              (2 * np.pi**(9/2) * l1*l2*(s1 + s2)*(s3 + s4) * np.sqrt(D)) * 
              np.exp(-g * (l1**2 + l2**2) * (s1 + s2) * (s3 + s4) / D))
    return(np.log(expEnt))
 
def entropyOfNet(netType, slVec, g):
    s1, s2, s3, s4, l1, l2 = slVec
    if netType == 'openNet0':
        return(openNet0Entropy(s1 = s1, g = g))
    if netType == 'closedNet0':
        return(closedNet0Entropy(s1 = s1, g = g))
    if netType == 'openNet1':
        return(openNet1Entropy(s1 = s1, l1 = l1, g = g))
    if netType == 'closedNet1':
        return(closedNet1Entropy(s1 = s1, s2 = s2, l1 = l1, g = g))
    if netType == 'openOpenNet2':
        return(openOpenNet2Entropy(s1 = s1, s2 = s2, l1 = l1, l2 = l2, g = g))
    if netType == 'openNet2a':
        return(openNet2aEntropy(s1 = s1, s2 = s2, s3 = s3, l1 = l1, l2 = l2, g = g))
    if netType == 'closedNet2a':
        return(closedNet2aEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, l1 = l1, l2 = l2, g = g))
    if netType == 'openNet2b':
        return(openNet2bEntropy(s1 = s1, s2 = s2, s3 = s3, l1 = l1, l2 = l2, g = g))
    if netType == 'closedNet2b':
        return(closedNet2bEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, l1 = l1, l2 = l2, g = g))
    if netType == 'openNet2c':
        return(openNet2cEntropy(s1 = s1, s2 = s2, s3 = s3, l1 = l1, l2 = l2, g = g))
    if netType == 'closedNet2c':
        return(closedNet2cEntropy(s1 = s1, s2 = s2, s3 = s3, s4 = s4, l1 = l1, l2 = l2, g = g))
    if netType == 'tooComplex':
        return(-10000)



# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================
        
def save(obj, filename):
# =============================================================================
#     Save an object to a file
# =============================================================================
    if not filename[-7:] == '.pickle':
        filename = filename + '.pickle'
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=4) #protocol 4 came out with Python version 3.4
    
def load(filename):
# =============================================================================
#     Load an object from a file
# =============================================================================
    if not filename[-7:] == '.pickle':
        filename = filename + '.pickle'
    with open(filename, 'rb') as handle:
        data = pickle.load(handle)
    return data

def writeToTxt(filename,textToWrite):
# =============================================================================
#     Output the progress update to the text file
# =============================================================================
    with open(filename, "a") as f: 
        f.write(textToWrite)
        f.write('\n')

def myPrintFxn(txtToPrint, fileTxt):
# =============================================================================
#     Print statement in editor and print it in the text file
# =============================================================================
    writeToTxt(fileTxt, txtToPrint)
    print(txtToPrint)
    
# =============================================================================
# *****************************************************************************
# *****************************************************************************
# =============================================================================

#Note for myself: For debugging, use ! before variable while debugging to not have it be used as a shortcut
#CGCGCAAGACACUGGUCAAGCGCUCC
#Example code:
q = RNALandscape(['AUCUGAUACUGUGCUAUGUCUGAGAUAGC'],storeGraphs=True, tryToLoadVariables = False) #define the object with the parameters you want
q.calculateFELandscape() #calculate the FE Lanscape
#The results are stored in, e.g., q.sortedFEs (the free energies, sorted from lowest to highest)
