function [listMyMFEBPs, MFEProbs, tp, fp, fn, sensitivity,ppv,whichStructureList,probStructureList,minBPInRegionList,probTopology,probTopologyPerBase,fracTopologyPerBase,whichTopology,totalNumStructures,totalNumTopologies,probPseudoknot,startAndPermuTimeVec,checkTimeVec,totalTimeVec] = ...
    compareToExperiments_Fxn(sequences,listStructBPs,maxNumRegions,vs,pairwiseEnergies,T,changeSubstems)
%used in conjunction with compareToExperiments.m
%takes a list of sequences and experimentally determined structures and
%finds the predicted structure and compares it to the experimentally
%determined structure.

numSeqs = length(sequences); %number of sequences/structures to check

listMyMFEBPs = cell(1,length(sequences)); %for each structure, what's the predicted MFE structure?
MFEProbs = zeros(1,numSeqs); %what's the probability of folding into that structure?

whichStructureList = zeros(1,numSeqs); %1 if the MFE structure is the same as the experimental, 2 if the predicted second lowest FE structure is the experimental, etc; 10^8 if the experimental structure isn't in the list of structures allowed by our chosen constraints on the algorithm.
probStructureList = zeros(1,numSeqs); %predicted probability of folding into the experimentally determined structure
fracBPPredictedInMFE = zeros(1,numSeqs); %same as tp divided by number of base pairs in experimental structure, so this is the same as sensitivity
fracBPMissedInMFE = zeros(1,numSeqs); %same as fn divided by number of base pairs in experimental structure
fracBPRoughlyPredictedInMFE = zeros(1,numSeqs); %same as fracBPPredictedInMFE but calling a bp i,j correctly predicted if experimentally there's a bp i,j+1 or i-1,j, etc.
fracBPRoughlyMissedInMFE = zeros(1,numSeqs); %same as fracBPMissedInMFE but calling a bp i,j correctly predicted if experimentally there's a bp i,j+1 or i-1,j, etc.
probBPRoughlyPredicted = zeros(1,numSeqs); %total prob of bp being roughly predicted
probBPRoughlyMissed = zeros(1,numSeqs); %total prob of bp being roughly missed
probBPExactlyPredicted = zeros(1,numSeqs); %total prob of bp being exactly predicted
probBPExactlyMissed = zeros(1,numSeqs); %total prob of bp being exactly missed
minBPInRegionList = zeros(1,numSeqs); %list of chosen constraints m = minimum number of base pairs in a stem for each sequence.

probTopology = zeros(1,numSeqs); %predicted probability of folding into the experimentally determined topology
probTopologyPerBase = zeros(1,numSeqs); %for each base pair, the predicted probability it'll be part of the right topology, averaged over the base pairs
fracTopologyPerBase = zeros(1,numSeqs); %actually never calculated in this function -- so it just stays a vector of zeros
whichTopology = zeros(1,numSeqs); %returns 1 if the MFE topology is the experimental one, 2 if it's the second MFE, etc.
totalNumStructures = zeros(1,numSeqs); %total number of possible structures tested for each sequence
totalNumTopologies = zeros(1,numSeqs); %total number of potential topologies predicted for each sequence
probPseudoknot = zeros(1,numSeqs); %predicted probability of each sequence folding into a pseudoknot

tp = zeros(1,numSeqs); %true positives -- base pairs in predicted MFE structure also in experimental structure
%tn = zeros(1,numPDBEntries); %meaningless in this context
fp = zeros(1,numSeqs); %false positives -- base pairs in predicted MFE structure not in experimental structure
fn = zeros(1,numSeqs); %false negatives -- base pairs not in predicted MFE structure which are in experimental structure

startAndPermuTimeVec = zeros(1,numSeqs);
checkTimeVec = zeros(1,numSeqs);
totalTimeVec = zeros(1,numSeqs);

%sensitivity is defined as the proportion of the number of correctly predicted bps to the number of bps of the known structure
%ppv is defined as the proportion of the number of correctly predicted bps to the number of bps predicted by the algorithms
%as in "prediction of RNA secondary structure with pseudoknots using
%integer programming" by Poolsap, Kato and Akutsu (BMC Bioinformatics
%2009).
sensitivity = zeros(1,numSeqs);
ppv = zeros(1,numSeqs);

for i = 1:numSeqs
    if rem(i,25) ==0
        disp(['Finished ', num2str(i),' out of ', num2str(numSeqs), ' sequences']);
    end
    
    sequence = sequences{i};
    strucBPs = listStructBPs{i};
    
    c1 = 1; c2 = 1; storeGraphs = true;
    makingFigures = false; printProgressUpdate = false; 
    minBPInHairpin = 3; numParForLoops = 5; allowParallelStrands = false; allowPseudoknots = true; minNumCompatible = 0; substems = 'all';
    
    [~,sequenceInNumbers,~,~,~,~,~,~,~,~,~] = multipleStrandsSetup_Fxn({sequence});

    minBPInRegion = 1;
    [numRegions, ~] = createSTable_Fxn(sequenceInNumbers,minBPInRegion,minBPInHairpin,allowParallelStrands,substems);
    while numRegions > maxNumRegions %make minBPInRegion the smallest it can be such that numRegions is not > maxNumRegions
        if changeSubstems %set to false in analyses
            if minBPInRegion >= 3 && strcmp(substems,'all')
                substems = 2;
            else
                minBPInRegion = minBPInRegion + 1; substems = 'all';
            end
        else
            minBPInRegion = minBPInRegion + 1;
        end
        [numRegions, ~] = createSTable_Fxn(sequenceInNumbers,minBPInRegion,minBPInHairpin,allowParallelStrands,substems);
    end
    if numRegions == 0
        minBPInRegion = minBPInRegion - 1;
    end
         
    minBPInRegionList(i) = minBPInRegion;
    [sortedProbs, indexSortedProbs,~, ~,graphProbs,weightMatrixList,numBondsMatrixList,startAndPermuTime,checkFxnTime,totalTime,STable,possiblePermutations,~,~] = ...
        RNALandscape_main({sequence},1,c1,c2,T,vs,0,minBPInRegion,numParForLoops,pairwiseEnergies,storeGraphs,makingFigures,...
        printProgressUpdate,allowParallelStrands, allowPseudoknots, minNumCompatible, substems);
    
    totalTimeVec(i) = totalTime;
    startAndPermuTimeVec(i) = startAndPermuTime;
    checkTimeVec(i) = checkFxnTime;
    totalNumStructures(i) = length(possiblePermutations);
    totalNumTopologies(i) = length(graphProbs);
    %%
    %what is the probability we calculate for folding into "structure"?
    foundStructure = false;
    counter = 0;
    while ~foundStructure
        counter = counter + 1;
        whichStruct = indexSortedProbs(counter);
        myStrucBPs = [];
        
        for j = 1:length(possiblePermutations{whichStruct})
            myStrucBPs = [myStrucBPs; STable{possiblePermutations{whichStruct}(j),3}];
        end
        
        myStrucBPs = sortrows(myStrucBPs);
        
        if isequal(strucBPs,myStrucBPs)
            foundStructure = true;
        end
        if counter == length(indexSortedProbs)
            break
        end
    end
    if foundStructure
        whichStructureList(i) = counter;
        probStructureList(i) = sortedProbs(counter);
    else 
        whichStructureList(i) = 10^8;
        probStructureList(i) = 0;
    end

    
    %now calculate how well we predicted structure.
    %first, consider the predicted MFE structure. How many of the base
    %pairs did it accurately predict?
    whichStruct = indexSortedProbs(1);
    myMFEBPs = [];
    for j = 1:length(possiblePermutations{whichStruct})
        myMFEBPs = [myMFEBPs; STable{possiblePermutations{whichStruct}(j),3}];
    end
    
    listMyMFEBPs{i} = myMFEBPs;
    MFEProbs(i) = sortedProbs(1);
    
    numBP = size(strucBPs,1);
    for k = 1:numBP
        bp = strucBPs(k,:); %base pair in structure
        
        %is this base pair in myMFEStructure?
        if ~isempty(myMFEBPs) && any(all(myMFEBPs' == bp'))
            bpFound = true;
        else
            bpFound = false;
        end
        
        if bpFound
            fracBPPredictedInMFE(i) = fracBPPredictedInMFE(i) + 1;
            tp(i) = tp(i) + 1;
        else
            fracBPMissedInMFE(i) = fracBPMissedInMFE(i) + 1;
            fn(i) = fn(i) + 1;
        end
    end
    
    numBP = size(myMFEBPs,1);
    for k = 1:numBP
        bp = myMFEBPs(k,:); %base pair in structure
        
        if ~isempty(strucBPs) && any(all(strucBPs' == bp'))
            bpFound = true;
        else
            bpFound = false;
        end
        if bpFound
            %tp(i) = tp(i) + 1; %already took this into account
        else
            fp(i) = fp(i) + 1;
        end
    end
    
    numBPsInMyStructure = size(myMFEBPs,1);
    numBPsInExpStructure = size(strucBPs,1);
    if numBPsInExpStructure > 0
        sensitivity(i) = tp(i)/numBPsInExpStructure;
    else
        sensitivity(i) = 0;
    end
    if numBPsInMyStructure > 0 
        ppv(i) = tp(i)/numBPsInMyStructure;
    else
        ppv(i) = 0;
    end
    
    
    
    %now, consider how many bps it predicted when we use the standard
    %approximation that "Known base-pairs are considered correctly
    %predicted if they occur in the predicted structure or if the predicted
    %structure contains a base-pair shifted by at most one nucleotide on
    %one side of the pair. For example, a known base-pair between nucleotides
    %i and j is considered correctly predicted by base-pairs of i to j,
    %i to j - 1,i to j + 1, i - 1 to j, or i + 1 to j. A predicted pair of
    %i + 1 to j - 1, however, is not considered correct."
    
    numBP = size(strucBPs,1);
    for k = 1:numBP
        bp = strucBPs(k,:); %base pair in structure
        
        %is this base pair in myMFEStructure?
        bpFound = false;
        bpTryMatrix = [bp;bp(1),bp(2)-1;bp(1),bp(2)+1;bp(1)-1,bp(2);bp(1)+1,bp(2)];
        for bpTry = 1:5
            if ~isempty(myMFEBPs) && any(all(myMFEBPs' == bpTryMatrix(bpTry,:)')) %if any of the modified bps are found in myMFEBPs
                bpFound = true;
                break
            end
        end
        
        if bpFound
            fracBPRoughlyPredictedInMFE(i) = fracBPRoughlyPredictedInMFE(i) + 1;
        else
            fracBPRoughlyMissedInMFE(i) = fracBPRoughlyMissedInMFE(i) + 1;
        end
    end
    
    %%
    %now calculate the total probability of predicting each bp found in
    %the MFE structure, and the probability of not calculating each bp not
    %found in the MFE structure. Could use this to generate a ROC curve
    
    numBP = size(strucBPs,1);
    totalNumBP = numBP;
    for k = 1:numBP
        bp = strucBPs(k,:); %base pair in structure
        
        %with what probability is this base pair predicted?
        whichStructCounter = 1;
        whichStruct = indexSortedProbs(1);
        while sortedProbs(whichStructCounter) > 1e-7 && whichStructCounter < length(sortedProbs) %arbitrary cutoff. Too low and the code takes too long; too high and probabilities aren't normalized
            myBPs = [];
            for l = 1:length(possiblePermutations{whichStruct})
                myBPs = [myBPs; STable{possiblePermutations{whichStruct}(l),3}]; %#ok<*AGROW>
            end
            
            if ~isempty(myBPs)
                bpFound = false;
                bpTryMatrix = [bp;bp(1),bp(2)-1;bp(1),bp(2)+1;bp(1)-1,bp(2);bp(1)+1,bp(2)];
                for bpTry = 1:5
                    if any(all(myBPs' == bpTryMatrix(bpTry,:)')) %if any of the modified bps are found in myMFEBPs
                        bpFound = true;
                        break
                    end
                end
                
                if bpFound && bpTry == 1
                    probBPExactlyPredicted(i) = probBPExactlyPredicted(i) + sortedProbs(whichStructCounter);
                    probBPRoughlyPredicted(i) = probBPRoughlyPredicted(i) + sortedProbs(whichStructCounter);
                elseif bpFound
                    probBPRoughlyPredicted(i) = probBPRoughlyPredicted(i) + sortedProbs(whichStructCounter);
                    probBPExactlyMissed(i) = probBPExactlyMissed(i) + sortedProbs(whichStructCounter);
                else
                    probBPRoughlyMissed(i) = probBPRoughlyMissed(i) + sortedProbs(whichStructCounter);
                    probBPExactlyMissed(i) = probBPExactlyMissed(i) + sortedProbs(whichStructCounter);
                end
            else
                probBPRoughlyMissed(i) = probBPRoughlyMissed(i) + sortedProbs(whichStructCounter);
                probBPExactlyMissed(i) = probBPExactlyMissed(i) + sortedProbs(whichStructCounter);
            end
            whichStructCounter = whichStructCounter + 1;
            whichStruct = indexSortedProbs(whichStructCounter);
        end
    end
    probBPRoughlyPredicted(i) = probBPRoughlyPredicted(i)/totalNumBP;
    probBPRoughlyMissed(i) = probBPRoughlyMissed(i)/totalNumBP;
    probBPExactlyPredicted(i) = probBPExactlyPredicted(i)/totalNumBP;
    probBPExactlyMissed(i) = probBPExactlyMissed(i)/totalNumBP;
    
    fracBPRoughlyPredictedInMFE(i) = fracBPRoughlyPredictedInMFE(i)/totalNumBP;
    fracBPPredictedInMFE(i) = fracBPPredictedInMFE(i)/totalNumBP;
    fracBPMissedInMFE(i) = fracBPMissedInMFE(i)/totalNumBP;
    fracBPRoughlyMissedInMFE(i) = fracBPRoughlyMissedInMFE(i)/totalNumBP;
    
    
    
    
    %calculate the probability for folding into the experimental topology
    numNtds = length(sequence);
    expStructureString = structBPs2structureString(numNtds,strucBPs);   
    [expWM,expNBM] = structureString2graph(expStructureString);
    
    if storeGraphs
        [sortedGraphProbs, sortOrder] = sort(graphProbs,'descend'); %sortedGraphProbs = graphProbs(sortOrder).
        sortedWMList = weightMatrixList(sortOrder);
        sortedNBMList = numBondsMatrixList(sortOrder);
        j = 1;
        while sortedGraphProbs(j) > 1e-8 && j < length(sortedGraphProbs) %otherwise not worth it to go through this code
            if myPseudoIsIsomorphic(expWM,sortedWMList{j},expNBM,sortedNBMList{j})
                probTopology(i) = sortedGraphProbs(j);
                whichTopology(i) = j;
                break
            end
            j = j + 1;
        end    
    end
    
    %calculate the total probability of predicting each ntds's topology found in
    %the experimental structure.
    [perBaseGraphExp] = structureString2perBaseGraph(expStructureString);
    
    whichStructCounter = 1;
    whichStruct = indexSortedProbs(1);
    while sortedProbs(whichStructCounter) > 2e-10 && whichStructCounter < length(sortedProbs) %arbitrary cutoff. Too low and the code takes too long; too high and probabilities aren't normalized
        myBPs = [];
        for l = 1:length(possiblePermutations{whichStruct})
            myBPs = [myBPs; STable{possiblePermutations{whichStruct}(l),3}]; %#ok<*AGROW>
        end
        predStructureString = structBPs2structureString(numNtds,myBPs);
        if ~storeGraphs
            [predWM,predNBM] = structureString2graph(predStructureString);
            if myPseudoIsIsomorphic(expWM,predWM,expNBM,predNBM)
                probTopology(i) = probTopology(i) + sortedProbs(whichStructCounter);
                if whichTopology(i) == 0
                    whichTopology(i) = whichStructCounter;
                end
            end
            
        end
        
        [perBaseGraphPred] = structureString2perBaseGraph(predStructureString);
        predPseudoknot = false;
        for k = 1:numNtds-1
            if strcmp(perBaseGraphPred{k},perBaseGraphExp{k})
                %fracTopologyPerBase(i) = fracTopologyPerBase(i) + 1/(numNtds-1);
                probTopologyPerBase(i) = probTopologyPerBase(i) + sortedProbs(whichStructCounter)/(numNtds-1);
            elseif (strcmp(perBaseGraphPred{k},'open_net_2_single') && strcmp(perBaseGraphExp{k},'open_net_2a_s1')) ||...
                    (strcmp(perBaseGraphPred{k},'open_net_2_single') && strcmp(perBaseGraphExp{k},'open_net_2a_s02')) ||... %open_net_2_single is equivalent to open nets 2a and 2b with l1 and l2 = 0
                    (strcmp(perBaseGraphPred{k},'open_net_2_single') && strcmp(perBaseGraphExp{k},'open_net_2c_s2')) ||...
                    (strcmp(perBaseGraphPred{k},'open_net_2_single') && strcmp(perBaseGraphExp{k},'open_net_2c_s01')) ||...
                    ...
                    (strcmp(perBaseGraphExp{k},'open_net_2_single') && strcmp(perBaseGraphPred{k},'open_net_2a_s1')) ||...
                    (strcmp(perBaseGraphExp{k},'open_net_2_single') && strcmp(perBaseGraphPred{k},'open_net_2a_s02')) ||... %open_net_2_single is equivalent to open nets 2a and 2b with l1 and l2 = 0
                    (strcmp(perBaseGraphExp{k},'open_net_2_single') && strcmp(perBaseGraphPred{k},'open_net_2c_s2')) ||...
                    (strcmp(perBaseGraphExp{k},'open_net_2_single') && strcmp(perBaseGraphPred{k},'open_net_2c_s01')) ||...
                    ...
                    (strcmp(perBaseGraphPred{k},'closed_net_2_single') && strcmp(perBaseGraphExp{k},'closed_net_2a_s12')) ||...
                    (strcmp(perBaseGraphPred{k},'closed_net_2_single') && strcmp(perBaseGraphExp{k},'closed_net_2a_s03')) ||... %closed_net_2_single is equivalent to closed nets 2a and 2b with l1 and l2 = 0
                    (strcmp(perBaseGraphPred{k},'closed_net_2_single') && strcmp(perBaseGraphExp{k},'closed_net_2c_s23')) ||...
                    (strcmp(perBaseGraphPred{k},'closed_net_2_single') && strcmp(perBaseGraphExp{k},'closed_net_2c_s01')) ||...
                    ...
                    (strcmp(perBaseGraphExp{k},'closed_net_2_single') && strcmp(perBaseGraphPred{k},'closed_net_2a_s12')) ||...
                    (strcmp(perBaseGraphExp{k},'closed_net_2_single') && strcmp(perBaseGraphPred{k},'closed_net_2a_s03')) ||... %closed_net_2_single is equivalent to closed nets 2a and 2b with l1 and l2 = 0
                    (strcmp(perBaseGraphExp{k},'closed_net_2_single') && strcmp(perBaseGraphPred{k},'closed_net_2c_s23')) ||...
                    (strcmp(perBaseGraphExp{k},'closed_net_2_single') && strcmp(perBaseGraphPred{k},'closed_net_2c_s01'))
                probTopologyPerBase(i) = probTopologyPerBase(i) + sortedProbs(whichStructCounter)/(numNtds-1);
            end
            if ~predPseudoknot && (strcmp(perBaseGraphPred{k},'closed_net_2c_s01')|| strcmp(perBaseGraphPred{k},'closed_net_2c_s23')|| ...
                    strcmp(perBaseGraphPred{k},'open_net_2c_s01')|| strcmp(perBaseGraphPred{k},'closed_net_2c_l12')|| ...
                    strcmp(perBaseGraphPred{k},'open_net_2c_s2')|| strcmp(perBaseGraphPred{k},'open_net_2c_l12')|| ...
                    strcmp(perBaseGraphPred{k},'closed_net_2b_l12')|| strcmp(perBaseGraphPred{k},'closed_net_2b_s13')|| ...
                    strcmp(perBaseGraphPred{k},'closed_net_2b_s02')|| strcmp(perBaseGraphPred{k},'open_net_2b_l2')|| ...
                    strcmp(perBaseGraphPred{k},'open_net_2b_l1')|| strcmp(perBaseGraphPred{k},'open_net_2b_s1')|| ...
                    strcmp(perBaseGraphPred{k},'open_net_2b_s02')|| strcmp(perBaseGraphPred{k},'closed_net_2a_l12')|| ...
                    strcmp(perBaseGraphPred{k},'closed_net_2a_s12')|| strcmp(perBaseGraphPred{k},'closed_net_2a_s03')|| ...
                    strcmp(perBaseGraphPred{k},'open_net_2a_l12')|| strcmp(perBaseGraphPred{k},'open_net_2a_s1') ||...
                    strcmp(perBaseGraphPred{k},'open_net_2a_s02')|| strcmp(perBaseGraphPred{k},'open_net_2_single') ||...
                    strcmp(perBaseGraphPred{k},'closed_net_2_single')|| strcmp(perBaseGraphPred{k},'closed_net_1_l1')||...
                    strcmp(perBaseGraphPred{k},'closed_net_1_s0')|| strcmp(perBaseGraphPred{k},'open_net_1_l1')||...
                    strcmp(perBaseGraphPred{k},'open_net_1_s0'))
                predPseudoknot = true;
            end
        end
        if predPseudoknot
            probPseudoknot(i) = probPseudoknot(i) + sortedProbs(whichStructCounter);
        end
        whichStructCounter = whichStructCounter + 1;
        whichStruct = indexSortedProbs(whichStructCounter);
    end
end

figure;
histogram(probStructureList,0:0.1:1)

disp(['mean sensitivity = ',num2str(mean(sensitivity))])
disp(['mean specificity = ',num2str(mean(ppv))])


disp(['max total time = ',num2str(max(totalTimeVec))]);