%compare topologies to experiment
%reads as input the table of prediction results, and generates Fig. 4 or Fig. S3 along
%with other figures. In order to generate Fig. S3, set ignoreNonCanonical
%and ignoreUnknown to true (to generate Fig. 4, set them to false). 
%In order to see results with N_stems^min = x, set maxNumRegions to x (x can be 150 or 200).

numNonPseudoknotSequences = 186; %number of nonpseudoknotted sequences
numPseudoknotSequences = 235; %number of pseudoknotted sequences
maxNumRegions = 200;
tab = readtable(['prediction_tools_comparison_NStemsMax_',num2str(maxNumRegions),'.xlsx']);

sequences = tab.sequence;
numSeqs = length(sequences);
indicesOfStructures = [4,26,32,44:6:size(tab,2)]; %simply gives the column number of each prediction tool we test
lengthOfSeqs = tab.length;

whichTopology = tab.pred_WhichTopology; %whichTopology considers the total probability of forming a structure with the right topology, not just the MFE structure

numI = length(indicesOfStructures);
correctlyPredictedTopology = zeros(numSeqs,numI);
fractionPredictedTopology = zeros(numSeqs,numI);
predictedParallel = zeros(numSeqs,numI); %did the predicted structure include any parallel loops?
avgPredictedParallelPerBase = zeros(numSeqs,numI); %what fraction of the base pairs were part of topologies which included parallel loops?
fracTopologyPerBase = zeros(numSeqs,numI);
predPseudoknot = zeros(numSeqs,numI);
expPseudoknot = zeros(numSeqs,1);
hasNonCanonicalPair = zeros(numSeqs,1);

tp = zeros(numSeqs,numI); %although this is already found in the table, we recalculate tp (true positives) etc. here, just to demonstrate how we calculated it
fp = zeros(numSeqs,numI);
fn = zeros(numSeqs,numI);
ppv = zeros(numSeqs,numI);
sensitivity = zeros(numSeqs,numI);
sensDiff = zeros(numSeqs,numI); %the fact that this vector ends up being close to zero just demonstrates that we correctly calculated the sensitivity in the table.
ppvDiff = zeros(numSeqs,numI);
tpDiff = zeros(numSeqs,numI);
fpDiff = zeros(numSeqs,numI);
fnDiff = zeros(numSeqs,numI);

top1Topologies = zeros(numSeqs,1);
top5Topologies = zeros(numSeqs,1);
top10Topologies = zeros(numSeqs,1);

for i = 1:numSeqs %calculate sensitivity and specificity
    %really the information for these is already in the table but here we just show how we calculated it (and
    %using the sensDiff etc. vectors show that the info in the table is accurate)
    expStrucBPs = structureString2structBPs(tab.experimental_Structure{i});
    
    for j = 1:size(expStrucBPs,1)
        ntd1 = sequences{i}(expStrucBPs(j,1));
        ntd2 = sequences{i}(expStrucBPs(j,2));
        if strcmp(ntd1,'A')
            if ~strcmp(ntd2,'U')
                hasNonCanonicalPair(i) = 1;
            end
        elseif strcmp(ntd1,'C')
            if ~strcmp(ntd2,'G')
                hasNonCanonicalPair(i) = 1;
            end
        elseif strcmp(ntd1,'G')
            if ~strcmp(ntd2,'C') && ~strcmp(ntd2,'U')
                hasNonCanonicalPair(i) = 1;
            end
        elseif strcmp(ntd1,'U')
            if ~strcmp(ntd2,'G') && ~strcmp(ntd2,'A')
                hasNonCanonicalPair(i) = 1;
            end
        end
    end
    
    
    numBPsInExpStructure = size(expStrucBPs,1);
    for j = 1:numI
        index = indicesOfStructures(j);
        predStrucBPs = structureString2structBPs(tab{i,index}{1});

        numBPExp = size(expStrucBPs,1);
        for k = 1:numBPExp
            bp = expStrucBPs(k,:); %base pair in experimental structure
            %is this base pair in predicted structure?
            if ~isempty(predStrucBPs) && any(all(predStrucBPs' == bp')) 
                bpFound = true;
            else
                bpFound = false;
            end

            if bpFound
                tp(i,j) = tp(i,j) + 1;
            else
                fn(i,j) = fn(i,j) + 1;
            end
        end

        numBPPred = size(predStrucBPs,1);
        for k = 1:numBPPred
            bp = predStrucBPs(k,:); %base pair in structure

            if ~isempty(expStrucBPs) && any(all(expStrucBPs' == bp'))
                bpFound = true;
            else
                bpFound = false;
            end
            if bpFound
                %tp(i) = tp(i) + 1; %already took this into account
            else
                fp(i,j) = fp(i,j) + 1;
            end
        end
    
        if numBPsInExpStructure > 0
            sensitivity(i,j) = tp(i,j)/numBPsInExpStructure;
        else
            sensitivity(i,j) = 0;
        end
        if numBPPred > 0
            ppv(i,j) = tp(i,j)/numBPPred;
        else
            ppv(i,j) = 0;
        end
        
        if j>1
            tpDiff(i,j) = tp(i,j) - tab{i,index+1};
            fpDiff(i,j) = fp(i,j) - tab{i,index+2};
            fnDiff(i,j) = fn(i,j) - tab{i,index+3};
            sensDiff(i,j) = sensitivity(i,j) - tab{i,index+4};
            ppvDiff(i,j) = ppv(i,j) - tab{i,index+5};
        end
    end
end

numNPtopologies = 0; %find number of experimental topologies for nonpseudoknotted sequences
numPtopologies = 0; %same for pseudoknotted
NPtopologiesWM = {}; %keep track of the topologies found via Weight Matrix and NumBonds Matrix
NPtopologiesNBM = {};
PtopologiesWM = {};
PtopologiesNBM = {};
unknown = zeros(1,numSeqs); %one if the topology formed by the experimental structure is more complex than the topologies we chose to allow in our algorithm
for i = 1:numSeqs
    [expWM,expNBM] = structureString2graph(tab.experimental_Structure{i});
    %find if this is a new experimental topology
    if i ==1 
        numNPtopologies = numNPtopologies+1;
        NPtopologiesWM{numNPtopologies} = expWM;
        NPtopologiesNBM{numNPtopologies} = expNBM;
    elseif i<=numNonPseudoknotSequences %ie for nonpseudoknots
        newTop = true;
        for j = 1:numNPtopologies
            if myPseudoIsIsomorphic(expWM,NPtopologiesWM{j},expNBM,NPtopologiesNBM{j})
                newTop = false;
            end
        end
        if newTop
            numNPtopologies = numNPtopologies+1;
            NPtopologiesWM{numNPtopologies} = expWM;
            NPtopologiesNBM{numNPtopologies} = expNBM;
        end
    elseif i == numNonPseudoknotSequences+1 %now for pseudoknots
        numPtopologies = numPtopologies+1;
        PtopologiesWM{numPtopologies} = expWM;
        PtopologiesNBM{numPtopologies} = expNBM;
    else
        newTop = true;
        for j = 1:numPtopologies
            if myPseudoIsIsomorphic(expWM,PtopologiesWM{j},expNBM,PtopologiesNBM{j})
                newTop = false;
            end
        end
        if newTop
            numPtopologies = numPtopologies+1;
            PtopologiesWM{numPtopologies} = expWM;
            PtopologiesNBM{numPtopologies} = expNBM;
        end 
    end
    
    
    for j = 1:numI
        [predWM,predNBM] = structureString2graph(tab{i,indicesOfStructures(j)}{1});
        correctlyPredictedTopology(i,j) = myPseudoIsIsomorphic(expWM,predWM,expNBM,predNBM);
    end
    [perBaseGraphExp] = structureString2perBaseGraph(tab.experimental_Structure{i});
    for k = 1:length(perBaseGraphExp)
        if strcmp(perBaseGraphExp{k},'unknown')
            unknown(i) = 1;
        end
    end
    for j = 1:numI
        [perBaseGraphPred] = structureString2perBaseGraph(tab{i,indicesOfStructures(j)}{1});
        %perBaseCorrectPrediction = strcmp(perBaseGraphExp,perBaseGraphPred);
        for k = 1:length(perBaseGraphExp)
%             if strcmp(perBaseGraphExp{k},'open_net_2_single') || strcmp(perBaseGraphPred{k},'open_net_2_single')...
%                     || strcmp(perBaseGraphExp{k},'closed_net_2_single') || strcmp(perBaseGraphPred{k},'closed_net_2_single')
%                 disp('Need to deal with net_2_single')
%             end
%             if contains(string(perBaseGraphPred{k}),"open_net_2c") || contains(string(perBaseGraphPred{k}),"closed_net_2c")...
%                     ||contains(string(perBaseGraphPred{k}),"open_net_2b") || contains(string(perBaseGraphPred{k}),"closed_net_2b")...
%                     ||contains(string(perBaseGraphPred{k}),"open_net_1") || contains(string(perBaseGraphPred{k}),"closed_net_1")...
%                 
%                 predictedParallel(i,j) = 1;
%                 avgPredictedParallelPerBase(i,j) = avgPredictedParallelPerBase(i,j) + 1;
%             end
            if strcmp(perBaseGraphPred{k},perBaseGraphExp{k})
                fracTopologyPerBase(i,j) = fracTopologyPerBase(i,j) + 1/(length(perBaseGraphExp));
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
                fracTopologyPerBase(i,j) = fracTopologyPerBase(i,j) + 1/(length(perBaseGraphExp));
            end
            if ~predPseudoknot(i,j) && (strcmp(perBaseGraphPred{k},'closed_net_2c_s01')|| strcmp(perBaseGraphPred{k},'closed_net_2c_s23')|| ...
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
                    strcmp(perBaseGraphPred{k},'open_net_1_s0')|| strcmp(perBaseGraphPred{k},'unknown'))
                predPseudoknot(i,j) = true;
            end
            if ~expPseudoknot(i) && (strcmp(perBaseGraphExp{k},'closed_net_2c_s01')|| strcmp(perBaseGraphExp{k},'closed_net_2c_s23')|| ...
                    strcmp(perBaseGraphExp{k},'open_net_2c_s01')|| strcmp(perBaseGraphExp{k},'closed_net_2c_l12')|| ...
                    strcmp(perBaseGraphExp{k},'open_net_2c_s2')|| strcmp(perBaseGraphExp{k},'open_net_2c_l12')|| ...
                    strcmp(perBaseGraphExp{k},'closed_net_2b_l12')|| strcmp(perBaseGraphExp{k},'closed_net_2b_s13')|| ...
                    strcmp(perBaseGraphExp{k},'closed_net_2b_s02')|| strcmp(perBaseGraphExp{k},'open_net_2b_l2')|| ...
                    strcmp(perBaseGraphExp{k},'open_net_2b_l1')|| strcmp(perBaseGraphExp{k},'open_net_2b_s1')|| ...
                    strcmp(perBaseGraphExp{k},'open_net_2b_s02')|| strcmp(perBaseGraphExp{k},'closed_net_2a_l12')|| ...
                    strcmp(perBaseGraphExp{k},'closed_net_2a_s12')|| strcmp(perBaseGraphExp{k},'closed_net_2a_s03')|| ...
                    strcmp(perBaseGraphExp{k},'open_net_2a_l12')|| strcmp(perBaseGraphExp{k},'open_net_2a_s1') ||...
                    strcmp(perBaseGraphExp{k},'open_net_2a_s02')|| strcmp(perBaseGraphExp{k},'open_net_2_single') ||...
                    strcmp(perBaseGraphExp{k},'closed_net_2_single')|| strcmp(perBaseGraphExp{k},'closed_net_1_l1')||...
                    strcmp(perBaseGraphExp{k},'closed_net_1_s0')|| strcmp(perBaseGraphExp{k},'open_net_1_l1')||...
                    strcmp(perBaseGraphExp{k},'open_net_1_s0')|| strcmp(perBaseGraphExp{k},'unknown'))
                expPseudoknot(i) = true;
            end
        end
        avgPredictedParallelPerBase(i,j) = avgPredictedParallelPerBase(i,j)/length(perBaseGraphExp);
        fractionPredictedTopology(i,j) = fracTopologyPerBase(i,j);%mean(strcmp(perBaseGraphExp,perBaseGraphPred));
    end
    
    if whichTopology(i) <= 1 && whichTopology(i) > 0 
        top1Topologies(i) = top1Topologies(i) + 1;
    end
    if whichTopology(i) <= 5 && whichTopology(i) > 0 
        top5Topologies(i) = top5Topologies(i) + 1;
    end
    if whichTopology(i) <= 10 && whichTopology(i) > 0 
        top10Topologies(i) = top10Topologies(i) + 1;
    end
end

disp(['average topologies predicted within 1 (non-pseudoknotted) = ',num2str(mean(top1Topologies(1:numNonPseudoknotSequences,:)))])
disp(['average topologies predicted within 1 (pseudoknotted) = ',num2str(mean(top1Topologies(numNonPseudoknotSequences+1:end,:)))])

disp(['average topologies predicted within 5 (non-pseudoknotted) = ',num2str(mean(top5Topologies(1:numNonPseudoknotSequences,:)))])
disp(['average topologies predicted within 5 (pseudoknotted) = ',num2str(mean(top5Topologies(numNonPseudoknotSequences+1:end,:)))])

disp(['average topologies predicted within 10 (non-pseudoknotted) = ',num2str(mean(top10Topologies(1:numNonPseudoknotSequences,:)))])
disp(['average topologies predicted within 10 (pseudoknotted) = ',num2str(mean(top10Topologies(numNonPseudoknotSequences+1:end,:)))])

ignoreUnknown = false; ignoreNonCanonical = false; %Set both to true to get supplementary figure S3
%ignoreUnknown says to not consider sequences which experimentally fold into topologies more complex than those
%we've solved for analytically; ignoreNonCanonical says to ignore sequences
%that have experimentally predicted base pairs that aren't AU, CG, or GU.
if ignoreUnknown
    sensitivity = sensitivity(~logical(unknown),:);  %ok to do this way because the first components of sensitivity, etc are the 
    %non-pseudoknotted (which are never unknown), and we always call 1:numNonPseudoknotSequences and numNonPseudoknotSequences+1:end
    correctlyPredictedTopology = correctlyPredictedTopology(~logical(unknown),:); 
    fractionPredictedTopology = fractionPredictedTopology(~logical(unknown),:); 
    ppv = ppv(~logical(unknown),:); 
    predPseudoknot = predPseudoknot(~logical(unknown),:); 
    
    top1Topologies = top1Topologies(~logical(unknown),:); 
    top5Topologies = top5Topologies(~logical(unknown),:); 
    top10Topologies = top10Topologies(~logical(unknown),:); 
    
    disp(['average topologies predicted within 1 (non-pseudoknotted) = ',num2str(mean(top1Topologies(1:numNonPseudoknotSequences,:)))])
    disp(['average topologies predicted within 1 (pseudoknotted) = ',num2str(mean(top1Topologies(numNonPseudoknotSequences+1:end,:)))])

    disp(['average topologies predicted within 5 (non-pseudoknotted) = ',num2str(mean(top5Topologies(1:numNonPseudoknotSequences,:)))])
    disp(['average topologies predicted within 5 (pseudoknotted) = ',num2str(mean(top5Topologies(numNonPseudoknotSequences+1:end,:)))])

    disp(['average topologies predicted within 10 (non-pseudoknotted) = ',num2str(mean(top10Topologies(1:numNonPseudoknotSequences,:)))])
    disp(['average topologies predicted within 10 (pseudoknotted) = ',num2str(mean(top10Topologies(numNonPseudoknotSequences+1:end,:)))])
end

if ignoreNonCanonical
    numNonPseudoknotSequences = sum(~expPseudoknot(~logical(hasNonCanonicalPair)));
    if ignoreUnknown
        hasNonCanonicalPair = hasNonCanonicalPair(~logical(unknown));
    end
    sensitivity = sensitivity(~logical(hasNonCanonicalPair),:); 
    correctlyPredictedTopology = correctlyPredictedTopology(~logical(hasNonCanonicalPair),:); 
    fractionPredictedTopology = fractionPredictedTopology(~logical(hasNonCanonicalPair),:); 
    ppv = ppv(~logical(hasNonCanonicalPair),:); 
    predPseudoknot = predPseudoknot(~logical(hasNonCanonicalPair),:); 
    
    top1Topologies = top1Topologies(~logical(hasNonCanonicalPair),:); 
    top5Topologies = top5Topologies(~logical(hasNonCanonicalPair),:); 
    top10Topologies = top10Topologies(~logical(hasNonCanonicalPair),:); 
        
    disp(['average topologies predicted within 1 (non-pseudoknotted) = ',num2str(mean(top1Topologies(1:numNonPseudoknotSequences,:)))])
    disp(['average topologies predicted within 1 (pseudoknotted) = ',num2str(mean(top1Topologies(numNonPseudoknotSequences+1:end,:)))])

    disp(['average topologies predicted within 5 (non-pseudoknotted) = ',num2str(mean(top5Topologies(1:numNonPseudoknotSequences,:)))])
    disp(['average topologies predicted within 5 (pseudoknotted) = ',num2str(mean(top5Topologies(numNonPseudoknotSequences+1:end,:)))])

    disp(['average topologies predicted within 10 (non-pseudoknotted) = ',num2str(mean(top10Topologies(1:numNonPseudoknotSequences,:)))])
    disp(['average topologies predicted within 10 (pseudoknotted) = ',num2str(mean(top10Topologies(numNonPseudoknotSequences+1:end,:)))])
end

disp(['average correctly predicted nonpseudoknotted topology = ',num2str(mean(correctlyPredictedTopology(1:numNonPseudoknotSequences,:)))])
disp(['average correctly predicted pseudoknotted topology = ',num2str(mean(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,:)))])
disp(['average correctly predicted nonpseudoknotted per-base topology = ',num2str(mean(fractionPredictedTopology(1:numNonPseudoknotSequences,:)))])
disp(['average correctly predicted pseudoknotted per-base topology = ',num2str(mean(fractionPredictedTopology(numNonPseudoknotSequences+1:end,:)))])
disp(['average sensitivity for nonpseudoknotted = ',num2str(mean(sensitivity(1:numNonPseudoknotSequences,:)))])
disp(['average sensitivity for pseudoknotted = ',num2str(mean(sensitivity(numNonPseudoknotSequences+1:end,:)))])
disp(['average ppv for nonpseudoknotted = ',num2str(mean(ppv(1:numNonPseudoknotSequences,:)))])
disp(['average ppv for pseudoknotted = ',num2str(mean(ppv(numNonPseudoknotSequences+1:end,:)))])

if numI > 1
[~,i] = max(mean(correctlyPredictedTopology(1:numNonPseudoknotSequences,2:end)));
disp(['best algorithm for nonpseudoknotted correctlyPredictedTopology was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(correctlyPredictedTopology(1:numNonPseudoknotSequences,i+1)))])

[~,i] = max(mean(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,2:end)));
disp(['best algorithm for pseudoknotted correctlyPredictedTopology was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,i+1)))])

[~,i] = max(mean(fractionPredictedTopology(1:numNonPseudoknotSequences,2:end)));
disp(['best algorithm for nonpseudoknotted per-base topology was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(fractionPredictedTopology(1:numNonPseudoknotSequences,i+1)))])

[~,i] = max(mean(fractionPredictedTopology(numNonPseudoknotSequences+1:end,2:end)));
disp(['best algorithm for pseudoknotted per-base topology was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(fractionPredictedTopology(numNonPseudoknotSequences+1:end,i+1)))])

[~,i] = max(mean(sensitivity(1:numNonPseudoknotSequences,2:end)));
disp(['best algorithm for nonpseudoknotted sensitivity was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(sensitivity(1:numNonPseudoknotSequences,i+1)))])

[~,i] = max(mean(sensitivity(numNonPseudoknotSequences+1:end,2:end)));
disp(['best algorithm for pseudoknotted sensitivity was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(sensitivity(numNonPseudoknotSequences+1:end,i+1)))])

[~,i] = max(mean(ppv(1:numNonPseudoknotSequences,2:end)));
disp(['best algorithm for nonpseudoknotted ppv was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(ppv(1:numNonPseudoknotSequences,i+1)))])

[~,i] = max(mean(ppv(numNonPseudoknotSequences+1:end,2:end)));
disp(['best algorithm for pseudoknotted ppv was ', tab.Properties.VariableNames{indicesOfStructures(i+1)},' = ', num2str(mean(ppv(numNonPseudoknotSequences+1:end,i+1)))])
end

%% Generate figures showing comparisons
%DEPENDENCIES: Uses Chad Greene's "Intuitive RGB color values from xkcd"
%function to generate colors.
colorVec = rand(numI,3);
colorNames = {'bright blue','kelly green','grass','moss','olive','yellowish green','leaf green','green','maroon','rose','pink','rust','scarlet','brick red','light red','fuchsia'};
for i = 1:numI
    colorVec(i,:) = rgb(colorNames{i});
end

labels = cell(1,numI);
labels{1} = 'Our algorithm';
emp = cell(1,numI); emp{1} = '';
for i= 2:numI
    labels{i} = tab.Properties.VariableNames{indicesOfStructures(i)}; 
    emp{i} = '';
end

labels = {'Our algorithm','RNAFold','Vienna (Andronescu)',...
    'MFold','CONTRAFold','PPfold','CentroidFold','ContextFold','HotKnots (DP)','HotKnots (RE)','HotKnots (CC)','ProbKnot','pknots','RNAPKplex','ILM'};

labels= labels(1:numI);
a = [];
figure;
hold on
for i = numI:-1:1
    a=[a,scatter(mean(sensitivity(1:numNonPseudoknotSequences,i)),mean(ppv(1:numNonPseudoknotSequences,i)),2000,colorVec(i,:),'filled','MarkerFaceAlpha',0.5)];
end

xlim([0.3 1]);
ylim([0.3 1]);
xlabel('nonpseudoknotted sensitivity')
ylabel('nonpseudoknotted ppv')
set(gca,'fontsize',18)
a = fliplr(a);
l = legend(a,labels);
set(l,'location','northwest')


a = [];
figure;
%subplot(1,2,2);
hold on
for i = numI:-1:1
    a=[a,scatter(mean(sensitivity(numNonPseudoknotSequences+1:end,i)),mean(ppv(numNonPseudoknotSequences+1:end,i)),2000,colorVec(i,:),'filled','MarkerFaceAlpha',0.5)];
end
xlim([0.3 1]);
ylim([0.3 1]);
xlabel('pseudoknotted sensitivity')
ylabel('pseudoknotted ppv')
set(gca,'fontsize',18)


%%
colorVec = rand(numI,3);
colorNames = {'bright blue','forest green','burnt sienna','kelly green','mint','royal blue','dull yellow','pink','maroon','aqua','lavender','sky blue','apple green','gold','light grey','faded purple'};
for i = 1:numI
    colorVec(i,:) = rgb(colorNames{i});
end

figure;
set(gcf, 'Position', [300, 500, 1000, 540])


subplot(2,6,1)
b1= scatter(1:15,100.*mean(correctlyPredictedTopology(1:numNonPseudoknotSequences,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(correctlyPredictedTopology(1:numNonPseudoknotSequences,:)),std(100.*correctlyPredictedTopology(1:numNonPseudoknotSequences,:))./sqrt(numNonPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])

ylim([0,100])
xlim([0,16])
xticklabels(emp)
newYTickLables = yticklabels;
for k = 2:length(newYTickLables)
    newYTickLables{k} = [newYTickLables{k},'%'];
end
yticklabels(newYTickLables)
ylabel('Predicted topology (acc.)')
set(gca,'fontsize',16)


subplot(2,6,2)
b2= scatter(1:15,100.*mean(fractionPredictedTopology(1:numNonPseudoknotSequences,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(fractionPredictedTopology(1:numNonPseudoknotSequences,:)),std(100.*fractionPredictedTopology(1:numNonPseudoknotSequences,:))./sqrt(numNonPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('Per-base topology (acc.)')
set(gca,'fontsize',16)


subplot(2,6,3)
b3= scatter(1:15,100.*mean(sensitivity(1:numNonPseudoknotSequences,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(sensitivity(1:numNonPseudoknotSequences,:)),std(100.*sensitivity(1:numNonPseudoknotSequences,:))./sqrt(numNonPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('Sensitivity')
set(gca,'fontsize',16)

subplot(2,6,4)
b4= scatter(1:15,100.*mean(ppv(1:numNonPseudoknotSequences,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(ppv(1:numNonPseudoknotSequences,:)),std(100.*ppv(1:numNonPseudoknotSequences,:))./sqrt(numNonPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('PPV')
set(gca,'fontsize',16)

subplot(2,6,5)
b4= scatter(1:15,100.*mean(predPseudoknot(1:numNonPseudoknotSequences,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(predPseudoknot(1:numNonPseudoknotSequences,:)),std(100.*predPseudoknot(1:numNonPseudoknotSequences,:))./sqrt(numNonPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('Predicted pseudoknots')
set(gca,'fontsize',16)

subplot(2,6,7)
b5= scatter(1:15,100.*mean(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,:)),std(100.*correctlyPredictedTopology(numNonPseudoknotSequences+1:end,:))./sqrt(numPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(newYTickLables)
ylabel('Predicted topology (acc.)')
set(gca,'fontsize',16)

subplot(2,6,8)
b6= scatter(1:15,100.*mean(fractionPredictedTopology(numNonPseudoknotSequences+1:end,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(fractionPredictedTopology(numNonPseudoknotSequences+1:end,:)),std(100.*fractionPredictedTopology(numNonPseudoknotSequences+1:end,:))./sqrt(numPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('Per-base topology (acc.)')
set(gca,'fontsize',16)

subplot(2,6,9)
b= scatter(1:15,100.*mean(sensitivity(numNonPseudoknotSequences+1:end,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(sensitivity(numNonPseudoknotSequences+1:end,:)),std(100.*sensitivity(numNonPseudoknotSequences+1:end,:))./sqrt(numPseudoknotSequences+1),'.','color',[0.3 0.3 0.3]) 
ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('Sensitivity')
set(gca,'fontsize',16)

subplot(2,6,10)
b= scatter(1:15,100.*mean(ppv(numNonPseudoknotSequences+1:end,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(ppv(numNonPseudoknotSequences+1:end,:)),std(100.*ppv(numNonPseudoknotSequences+1:end,:))./sqrt(numPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])


ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('PPV')
set(gca,'fontsize',16)
xlim([-0.2 numI+1.2])
yticks(yticks)


subplot(2,6,11)
b= scatter(1:15,100.*mean(predPseudoknot(numNonPseudoknotSequences+1:end,:)),40,colorVec,'filled');
hold on
errorbar(100.*mean(predPseudoknot(numNonPseudoknotSequences+1:end,:)),std(100.*predPseudoknot(numNonPseudoknotSequences+1:end,:))./sqrt(numPseudoknotSequences+1),'.','color',[0.3 0.3 0.3])

ylim([0,100])
xlim([0,16])

xticklabels(emp)
yticklabels(emp)
ylabel('Predicted pseudoknots')
set(gca,'fontsize',16)
xlim([-0.2 numI+1.2])
yticks(yticks)

[~,h3]=suplabel('No Pseudoknots','y',[.075 .562 .65 .4]);
set(h3,'fontsize',26,'color',rgb('black'),'fontweight','bold')
[~,h3]=suplabel('Pseudoknots','y',[.075 .08 .65 .4]);
set(h3,'fontsize',26,'color',rgb('black'),'fontweight','bold')

a = [];
for i = 1:numI
    a = [a,bar(10000+i,10000+i,'facecolor',colorVec(i,:))];%,'markersize',14)];
end

subplot(2,6,[6 12])
axis off
set(gca,'color','none')
legend(a,labels,'Position',[0.855 0.263 0.05 0.5])

%%

if ~ignoreUnknown && ~ignoreNonCanonical
    %Estimate how much our heuristics hurt us
    %takes several minutes to run
    minBPInRegionTooHigh = zeros(1, numSeqs); %was the minimum length of stems we allowed too large compared to the minimum length of the experimental stem?
    minBPinRegionDist = zeros(1, numSeqs); %how too high was it
    singleBPs = zeros(1,numSeqs); %how many times do we have single bps in the experimental dataset?
    index = indicesOfStructures(1);
    for i = 1:numSeqs
        expStruc = structBPs2structure(structureString2structBPs(tab.experimental_Structure{i}));
        predStruc = structBPs2structure(structureString2structBPs(tab{i,index}{1}));
        minLengthStemExp = 100;
        for j = 1:length(expStruc) %how long is the shortest stem in experimental structure?
            if length(expStruc{j})/2 < minLengthStemExp
                minLengthStemExp = length(expStruc{j})/2;
            end
        end

        %find the minBPInRegion we allowed
        sequence = sequences{i};
        [~,sequenceInNumbers,~,~,~,~,~,~,~,~,~] = multipleStrandsSetup_Fxn({sequence});

        minBPInHairpin = 3; allowParallelStrands = false; substems = 'all';

        minBPInRegion = 1; 
        [numRegions, ~] = createSTable_Fxn(sequenceInNumbers,minBPInRegion,minBPInHairpin,allowParallelStrands,substems);
        while numRegions > maxNumRegions 
            minBPInRegion = minBPInRegion + 1;
            [numRegions, ~] = createSTable_Fxn(sequenceInNumbers,minBPInRegion,minBPInHairpin,allowParallelStrands,substems);
        end
        if numRegions == 0
            minBPInRegion = minBPInRegion - 1;
        end

        if minBPInRegion > minLengthStemExp
            minBPInRegionTooHigh(i) = 1;
            minBPinRegionDist(i) = minBPInRegion - minLengthStemExp;
        end
        if minLengthStemExp == 1
            singleBPs(i) = 1;
        end
    end
    sum(minBPInRegionTooHigh(correctlyPredictedTopology(1:numNonPseudoknotSequences,1)~=1))

    sum(minBPInRegionTooHigh(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,1) ~= 1))




    sum(unknown(1,:)' & ~correctlyPredictedTopology(:,1) & expPseudoknot(:,1))
    sum(~correctlyPredictedTopology & expPseudoknot)
    numNPtopologies
    numPtopologies
    sum(hasNonCanonicalPair & expPseudoknot)
    sum(hasNonCanonicalPair & ~expPseudoknot)

    sum(minBPInRegionTooHigh(1:numNonPseudoknotSequences))
    sum(minBPInRegionTooHigh(numNonPseudoknotSequences+1:end))
end