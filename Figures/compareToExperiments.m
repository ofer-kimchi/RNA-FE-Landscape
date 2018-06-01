minDist = 0.2; %get rid of sequences within a Jukes-Cantor distance of minDist
changeSubstems = false; %do we want to allow the parameter "substems" to deviate from "all"?
maxNumRegions = 200; %was 150 to make figure that's in paper currently (Trying to change it on 5/8/18 and will update paper if necessary)

%% First, make a list of the sequences that we're going to test
%load sequences from Pseudobase++ database
[pseudoknotSequences, pseudoknotListStructBPs,pseudoknotsStructuresToPrint] = readBPSEQFile('pseudobasePPSequences80Ntds.txt',0);
%don't get rid of sequences because they're too similar at this stage,
%since we're going to do that soon.
%pseudoknotSequences = {}; pseudoknotListStructBPs = {}; pseudoknotsStructuresToPrint = {};

vs = 0.02;
pairwiseEnergies = true;
T = 300; 

correctStructureNonCanonical = false; %false means don't remove base pairs that aren't AU, GC, or GU
correctStructureHairpin = true; %true means remove base pairs until the minimum length of any hairpin is 3 ntds. We ultimately remove any sequences we had to modify in this way from our analysis.

%load sequences from RNAStrand database
folderName = 'RNAStrand80Ntds/';

allFiles = dir(folderName); %all files for this sequence
allFiles = struct2cell(allFiles);
fileNames = allFiles(1,:); %names of all files in this folder
realFileIndices = regexp(fileNames,strcat('.*.txt'));
realFileIndices = find(~cellfun(@isempty,realFileIndices));
fileNames = fileNames(realFileIndices); %cell array with all the files in folderName ending in .txt

%make a list of sequences and structures 
%get rid of files with the same sequence and structure (duplicates) and
%make a list of the files that have the same sequence but different structures
sequences = cell(1,length(fileNames));
nonPseudoknotStructures = cell(1,length(fileNames));
duplicates = cell(1);
numDuplicates = 0;
sameSequences = cell(1); %list of files that have the same sequence but different structures
numSameSequences = 0;
numCorrectedHairpins = 0;
correctedHairpins = cell(1); %list of files that have too short hairpins. We don't include these sequences in our analysis.

i = 1;
while i <= length(fileNames)

    fileName = strcat(folderName,fileNames{i});
    isDuplicate = false;
    isSameSequence = false;
    
    [sequence, structure,correctedNC,correctedHairpin] = readCTFile(fileName,correctStructureNonCanonical,correctStructureHairpin);
    
    for j = 1:i-1
        if strcmp(sequences{j},sequence) && isequal(nonPseudoknotStructures{j},structure)
            isDuplicate = true; 
            break
        end
    end
    
    if ~isDuplicate %want to make sure we get rid of duplicates before checking 
                %if we have some with the same sequence and different structures
        for j = 1:i-1 
            if strcmp(sequences{j},sequence)
                isSameSequence = true; 
                break
            end
        end
    end
    
    nonAUCG = 0;
%     for j = 1:length(sequence)
%         if ~strcmp(sequence(j),'A') && ~strcmp(sequence(j),'C') && ~strcmp(sequence(j),'G') &&...
%                 ~strcmp(sequence(j),'U')
%             nonAUCG = 1;
%         end
%     end
    
    if ~isDuplicate && ~isSameSequence && ~nonAUCG && ~correctedHairpin
        sequences{i} = sequence;
        nonPseudoknotStructures{i} = structure;
        i=i+1;
    elseif isDuplicate
        numDuplicates = numDuplicates + 1;
        duplicates{numDuplicates} = {fileNames{j},fileNames{i}};
        fileNames(i) = [];
    elseif isSameSequence
        numSameSequences = numSameSequences + 1;
        sameSequences{numSameSequences} = {fileNames{j},fileNames{i}};
        fileNames(i) = []; %don't want to test the same sequence for two different structures.
%     elseif nonAUCG
%         fileNames(i) = [];
    elseif correctedHairpin 
        numCorrectedHairpins = numCorrectedHairpins + 1;
        correctedHairpins{numCorrectedHairpins} = {fileNames{i}};
        fileNames(i) = []; 
    end
end

sequences(i:end) = [];
nonPseudoknotStructures(i:end) = [];


listStructBPs = {};
for i = 1:length(nonPseudoknotStructures)
    listStructBPs{length(listStructBPs)+1} = structure2structBPs(nonPseudoknotStructures{i}); %#ok<*SAGROW>
end


%separate pseudoknotted and not pseudoknotted structures
i=1;
while i <= length(nonPseudoknotStructures)
    if hasPseudoknot(nonPseudoknotStructures{i})
        pseudoknotListStructBPs = [pseudoknotListStructBPs,listStructBPs{i}]; %#ok<*AGROW>
        pseudoknotSequences = [pseudoknotSequences,sequences{i}];
        sequences(i) = [];
        listStructBPs(i) = []; 
        nonPseudoknotStructures(i) = [];
    else 
        i = i+1;
    end
end


listStructBPs = [listStructBPs,pseudoknotListStructBPs]; %listStructBPs is a cell array where the i'th element lists the experimentally determined base pairs in structure i
sequences = [sequences, pseudoknotSequences]; %list of sequences
%%
%load sequences from CompraRNA database
[newSequences, newSeqListStructBPs,newStructuresToPrint] = readBPSEQFile('CompraRNA_PDB.txt',0);
i = 1;
while i <= length(newSequences) 
    if length(newSequences{i})>80
        newSequences(i) = []; newSeqListStructBPs(i) = []; newStructuresToPrint(i) = [];
    else
        i=i+1;
    end
end

%get rid of sequences that we already have
i = 1;
while i <= length(newSequences) 
    isDup = false;
    for j = 1:length(sequences)
        if strcmp(sequences{j},newSequences{i})
            isDup = true;
            if ~isequal(newSeqListStructBPs{i},listStructBPs{j}) %meaning that CompraRNA database has the same sequence as found in RNAStrand but a different experimentally determined structure
                %disp('CompraRNA database isn''t giving same results as RNAStrand or Pseudobase++')
                %we trust CompraRNA database more than RNAStrand (I checked
                %the cases where they disagree by hand).
                isDup = false;
                sequences(j) = []; listStructBPs(j) = [];
            end
            break
        end
    end
    if isDup
        newSequences(i) = []; newSeqListStructBPs(i) = []; newStructuresToPrint(i) = [];
    else
        i=i+1;
    end
end

listStructBPs = [listStructBPs,newSeqListStructBPs];
sequences = [sequences, newSequences];
%%

%get rid of sequences which have some ntds unspecified
i = 1;
while i <= length(sequences)
    unspecifiedNtds = false;
    for j = 1:length(sequences{i})
        if sequences{i}(j) ~= 'A' && sequences{i}(j) ~= 'C' && ...
                sequences{i}(j) ~= 'G' && sequences{i}(j) ~= 'U'
            sequences(i) = [];
            listStructBPs(i) = [];
            unspecifiedNtds = true;
            break
        end
    end
    if ~unspecifiedNtds
        i = i+1;
    end
end


%make sure we got rid of non-canonical base pairs if that's what we wanted
%to do (but we have this set to false)
if correctStructureNonCanonical
    for i = 1:length(listStructBPs)
        for j = 1:size(listStructBPs{i},1)
            ntd1 = sequences{i}(listStructBPs{i}(j,1));
            ntd2 = sequences{i}(listStructBPs{i}(j,2));
            if strcmp(ntd1,'A')
                if ~strcmp(ntd2,'U')
                    disp(['STRUCTURE ',num2str(i),' HAS NON-CANONICAL BPS'])
                end
            elseif strcmp(ntd1,'C')
                if ~strcmp(ntd2,'G')
                    disp(['STRUCTURE ',num2str(i),' HAS NON-CANONICAL BPS'])
                end 
            elseif strcmp(ntd1,'G')
                if ~strcmp(ntd2,'C') && ~strcmp(ntd2,'U')
                    disp(['STRUCTURE ',num2str(i),' HAS NON-CANONICAL BPS'])
                end 
            elseif strcmp(ntd1,'U')
                if ~strcmp(ntd2,'G') && ~strcmp(ntd2,'A')
                    disp(['STRUCTURE ',num2str(i),' HAS NON-CANONICAL BPS'])
                end 
            end
        end
    end
end


%Now that we have all the sequences from all three databases, get rid of sequences within a Jukes-Cantor distance of minDist

seqD = seqpdist(sequences,'Alphabet','NT'); %measures pairwise Jukes-Cantor distances between sequences
seqDMat = zeros(length(sequences)); %matrix of pairwise distances between sequences
counter = 0;
for j = 1:length(sequences)
    for i = j+1:length(sequences)
        counter = counter + 1;
        seqDMat(i,j) = seqD(counter);
        seqDMat(j,i) = seqD(counter); %to make dMat symmetric
    end
end

%now we want to generate a new set of sequences (and structures) that
%doesn't include sequences that are very near each other in sequence space
%and structure space

indexCulledSequences = []; %sequences to remove
for i = 1:length(sequences)
    seqDMat(i,i) = 1e6; %so we don't worry about diagonal distances being < minDist
end

for i = 1:length(sequences)
    if ~any(indexCulledSequences == i) && any(seqDMat(i,:)<minDist) 
        nearbySequences = setdiff(find(seqDMat(i,:)<minDist),indexCulledSequences); %setdiff(A,B) gives elements of A not in B
        %disp([i,nearbySequences])
        indexCulledSequences = [indexCulledSequences, nearbySequences];
    end
end

sequences(indexCulledSequences) = [];
listStructBPs(indexCulledSequences) = [];

hasNCBPs = zeros(1,length(sequences)); %for each sequence, does it include base pairs outside AU, GU, and GC?
for i = 1:length(listStructBPs)
    for j = 1:size(listStructBPs{i},1)
        ntd1 = sequences{i}(listStructBPs{i}(j,1));
        ntd2 = sequences{i}(listStructBPs{i}(j,2));
        if strcmp(ntd1,'A')
            if ~strcmp(ntd2,'U')
                hasNCBPs(i) = 1;
            end
        elseif strcmp(ntd1,'C')
            if ~strcmp(ntd2,'G')
                hasNCBPs(i) = 1;
            end
        elseif strcmp(ntd1,'G')
            if ~strcmp(ntd2,'C') && ~strcmp(ntd2,'U')
                hasNCBPs(i) = 1;
            end
        elseif strcmp(ntd1,'U')
            if ~strcmp(ntd2,'G') && ~strcmp(ntd2,'A')
                hasNCBPs(i) = 1;
            end
        end
    end
end




%find number of pseudoknotted and non-pseudoknotted sequences
for i = 1:length(pseudoknotSequences)
    IndexC = strfind(sequences, pseudoknotSequences{i});
    indexFirstPseudoknot = find(not(cellfun('isempty', IndexC)));
    if ~isempty(indexFirstPseudoknot)
        numNonPseudoknotSequences = indexFirstPseudoknot - 1;
        numPseudoknotSequences = length(sequences) - numNonPseudoknotSequences;
        break
    end
end

NPSequences = {}; NPListStructBPs = {}; PSequences = {}; PListStructBPs = {};
numPseudoknotSequences = 0; numNonPseudoknotSequences = 0;
oldSequences = sequences; oldListStructBPs = listStructBPs;
for i = 1:length(sequences)
    if hasPseudoknotStructBPs(listStructBPs{i})
        numPseudoknotSequences = numPseudoknotSequences + 1;
        PSequences{numPseudoknotSequences} = sequences{i};
        PListStructBPs{numPseudoknotSequences} = listStructBPs{i};
    else
        numNonPseudoknotSequences = numNonPseudoknotSequences + 1;
        NPSequences{numNonPseudoknotSequences} = sequences{i};
        NPListStructBPs{numNonPseudoknotSequences} = listStructBPs{i};
    end
end
sequences = [NPSequences,PSequences];
listStructBPs = [NPListStructBPs,PListStructBPs];


disp([num2str(numNonPseudoknotSequences), ' non-pseudoknotted structures']);
disp([num2str(numPseudoknotSequences), ' pseudoknotted structures']);


%% take our list of sequences and experimentally determined structures and find the predicted structures. Then compare them to the experimentally determined structure.
%sequences = fliplr(sequences); listStructBPs = fliplr(listStructBPs);
[listMyMFEBPs, MFEProbs, tp, fp, fn, sensitivity,ppv,whichStructureList, probStructureList,minBPInRegionList,probTopology,probTopologyPerBase,fracTopologyPerBase,whichTopology,totalNumStructures,totalNumTopologies,probPseudoknot,startAndPermuTimeVec,checkTimeVec,totalTimeVec] = ...
   compareToExperimentsFxn4(sequences,listStructBPs,maxNumRegions,vs,pairwiseEnergies,T,changeSubstems);

%%
minBPInRegionHist = hist(minBPInRegionList,1:max(minBPInRegionList));
disp(['minBPInRegion fractions: ',num2str(minBPInRegionHist./sum(minBPInRegionHist))])


disp([num2str(numNonPseudoknotSequences), ' non-pseudoknotted structures']);
disp([num2str(numPseudoknotSequences), ' pseudoknotted structures']);


structureMyMFEStrings = {};
for i = 1:length(listMyMFEBPs)
    sequence = sequences{i};
    structBPs = listMyMFEBPs{i};
    
    printedStructCell = cell(1,length(sequence)); %way to print structure. For now, a cell array where the j'th element is equal to k if ntd j is paired with k, and 0 otherwise
    numBPs = size(structBPs,1);
    for j = 1:length(printedStructCell)
        paired = find(structBPs == j);
        if isempty(paired)
            printedStructCell{j} = '0 ';
        else
            if paired <= numBPs
                printedStructCell{j} = [num2str(structBPs(paired+numBPs)),' ']; %don't use strcat to preserve spaces
            else
                printedStructCell{j} = [num2str(structBPs(paired-numBPs)),' '];
            end
        end
    end
    printedStructString = '';
    for j = 1:length(printedStructCell)
        printedStructString = [printedStructString,printedStructCell{j}]; 
    end
    structureMyMFEStrings{length(structureMyMFEStrings)+1} = printedStructString; 
end

structureExpStrings = {};
for i = 1:length(listStructBPs)
    sequence = sequences{i};
    structBPs = listStructBPs{i};
    
    printedStructCell = cell(1,length(sequence)); %way to print structure. For now, a cell array where the j'th element is equal to k if ntd j is paired with k, and 0 otherwise
    numBPs = size(structBPs,1);
    for j = 1:length(printedStructCell)
        paired = find(structBPs == j);
        if isempty(paired)
            printedStructCell{j} = '0 ';
        else
            if paired <= numBPs
                printedStructCell{j} = [num2str(structBPs(paired+numBPs)),' ']; %don't use strcat to preserve spaces
            else
                printedStructCell{j} = [num2str(structBPs(paired-numBPs)),' '];
            end
        end
    end
    printedStructString = '';
    for j = 1:length(printedStructCell)
        printedStructString = [printedStructString,printedStructCell{j}]; 
    end
    structureExpStrings{length(structureExpStrings)+1} = printedStructString;
end



lengthSequences = zeros(1,length(sequences));
for i = 1:length(lengthSequences)
    lengthSequences(i) = length(sequences{i});
end

resultsTable = table(sequences',lengthSequences',structureExpStrings',structureMyMFEStrings',...
    MFEProbs', tp', fp', fn', sensitivity',ppv',whichStructureList', probStructureList',probTopology',probTopologyPerBase',whichTopology',totalNumStructures',totalNumTopologies',probPseudoknot');
resultsTable.Properties.VariableNames = {'sequence', 'length','experimental_Structure','pred_Structure',...
    'pred_Prob','pred_TP','pred_FP','pred_FN','pred_Sensitivity','pred_PPV','pred_WhichStructExp','pred_ProbExpStruct','pred_ProbTopology','pred_ProbTopologyPerBase','pred_WhichTopology','pred_totalNumStructures','pred_totalNumTopologies','pred_probPseudoknot'};


writetable(resultsTable,['our_prediction_',...
    'minDist_', num2str(minDist),'_maxNumRegions_',num2str(maxNumRegions),'_vs_',num2str(vs),'_pairwiseEnergies_',num2str(pairwiseEnergies),'_corStrucNC_',...
    num2str(correctStructureNonCanonical),'_corStrucHairpin_',num2str(correctStructureHairpin),'_changeSubstems_',num2str(changeSubstems),'.txt'],'Delimiter',',');

%% %compare topologies to experiment. Generates Fig. 5 of the paper
%This section actually doesn't affect anything -- it's just for comparison
%to the real analysis which comes in a different file. 
numSeqs = length(sequences);
correctlyPredictedTopology = zeros(numSeqs,1);
fractionPredictedTopology = zeros(numSeqs,1);
predictedParallel = zeros(numSeqs,1); %did the predicted structure include any parallel loops?
avgPredictedParallelPerBase = zeros(numSeqs,1); %what fraction of the base pairs were part of topologies which included parallel loops?

top1Topologies = zeros(numSeqs,1); %did the prediction algorithm predict the correct topology?
top5Topologies = zeros(numSeqs,1); %did the prediction algorithm predict the correct topology within the top 5?
top10Topologies = zeros(numSeqs,1); %did the prediction algorithm predict the correct topology within the top 10?

top1PercentTopologies = zeros(numSeqs,1); %did the prediction algorithm predict the correct topology within the top 1%?
top5PercentTopologies = zeros(numSeqs,1); %did the prediction algorithm predict the correct topology within the top 5%?
top10PercentTopologies = zeros(numSeqs,1); %did the prediction algorithm predict the correct topology within the top 10%?

unknown = zeros(numSeqs,1); %is the pseudoknot too complicated for our algorithm to predict?
for i = 1:numSeqs
    [expWM,expNBM] = structureString2graph(resultsTable.experimental_Structure{i});
    
    [predWM,predNBM] = structureString2graph(resultsTable{i,4}{1});
    correctlyPredictedTopology(i) = myPseudoIsIsomorphic(expWM,predWM,expNBM,predNBM);
    
    [perBaseGraphExp] = structureString2perBaseGraph(resultsTable.experimental_Structure{i});
    
    [perBaseGraphPred] = structureString2perBaseGraph(resultsTable{i,4}{1});
    %perBaseCorrectPrediction = strcmp(perBaseGraphExp,perBaseGraphPred);
    for k = 1:length(perBaseGraphExp)
%         if strcmp(perBaseGraphExp{k},'open_net_2_single') || strcmp(perBaseGraphPred{k},'open_net_2_single')...
%                 || strcmp(perBaseGraphExp{k},'closed_net_2_single') || strcmp(perBaseGraphPred{k},'closed_net_2_single')
%             disp('Need to deal with net_2_single')
%         end
%         if contains(string(perBaseGraphPred{k}),"open_net_2c") || contains(string(perBaseGraphPred{k}),"closed_net_2c")...
%                 ||contains(string(perBaseGraphPred{k}),"open_net_2b") || contains(string(perBaseGraphPred{k}),"closed_net_2b")...
%                 ||contains(string(perBaseGraphPred{k}),"open_net_1") || contains(string(perBaseGraphPred{k}),"closed_net_1")
%             predictedParallel(i) = 1;
%             avgPredictedParallelPerBase(i) = avgPredictedParallelPerBase(i) + 1;
%         end
        if strcmp(perBaseGraphExp{k},'unknown')
            unknown(i) = 1;
        end
    end
    avgPredictedParallelPerBase(i) = avgPredictedParallelPerBase(i)/length(perBaseGraphExp);
    fractionPredictedTopology(i) = mean(strcmp(perBaseGraphExp,perBaseGraphPred));
    
    if whichTopology(i) <= 1 && whichTopology(i) > 0 
        top1Topologies(i) = top1Topologies(i) + 1;
    end
    if whichTopology(i) <= 5 && whichTopology(i) > 0 
        top5Topologies(i) = top5Topologies(i) + 1;
    end
    if whichTopology(i) <= 10 && whichTopology(i) > 0 
        top10Topologies(i) = top10Topologies(i) + 1;
    end
    if whichTopology(i) <= ceil(totalNumTopologies(i)*0.01) && whichTopology(i) > 0 
        top1PercentTopologies(i) = top1PercentTopologies(i) + 1;
    end
    if whichTopology(i) <= ceil(totalNumTopologies(i)*0.05) && whichTopology(i) > 0 
        top5PercentTopologies(i) = top5PercentTopologies(i) + 1;
    end
    if whichTopology(i) <= ceil(totalNumTopologies(i)*0.1) && whichTopology(i) > 0 
        top10PercentTopologies(i) = top10PercentTopologies(i) + 1;
    end
end
unknown = logical(unknown);

%disp(['average correctly predicted topology = ',num2str(mean(correctlyPredictedTopology))])
disp(['fraction of MFE structures which had the same topology as experimental (non-pseudoknotted) = ',num2str(mean(correctlyPredictedTopology(1:numNonPseudoknotSequences,:)))])
disp(['fraction of MFE structures which had the same topology as experimental (pseudoknotted) = ',num2str(mean(correctlyPredictedTopology(numNonPseudoknotSequences+1:end,:)))])

disp(['average topologies predicted within 1 (non-pseudoknotted) = ',num2str(mean(top1Topologies(1:numNonPseudoknotSequences,:)))])
disp(['average topologies predicted within 1 (pseudoknotted) = ',num2str(mean(top1Topologies(numNonPseudoknotSequences+1:end,:)))])

disp(['average topologies predicted within 5 (non-pseudoknotted) = ',num2str(mean(top5Topologies(1:numNonPseudoknotSequences,:)))])
disp(['average topologies predicted within 5 (pseudoknotted) = ',num2str(mean(top5Topologies(numNonPseudoknotSequences+1:end,:)))])

disp(['average topologies predicted within 10 (non-pseudoknotted) = ',num2str(mean(top10Topologies(1:numNonPseudoknotSequences,:)))])
disp(['average topologies predicted within 10 (pseudoknotted) = ',num2str(mean(top10Topologies(numNonPseudoknotSequences+1:end,:)))])
% 
% disp(['average topologies predicted within 1% (non-pseudoknotted) = ',num2str(mean(top1PercentTopologies(1:54,:)))])
% disp(['average topologies predicted within 1% (pseudoknotted) = ',num2str(mean(top1PercentTopologies(55:197,:)))])
% 
% disp(['average topologies predicted within 5% (non-pseudoknotted) = ',num2str(mean(top5PercentTopologies(1:54,:)))])
% disp(['average topologies predicted within 5% (pseudoknotted) = ',num2str(mean(top5PercentTopologies(55:197,:)))])
% 
% disp(['average topologies predicted within 10% (non-pseudoknotted) = ',num2str(mean(top10PercentTopologies(1:54,:)))])
% disp(['average topologies predicted within 10% (pseudoknotted) = ',num2str(mean(top10PercentTopologies(55:197,:)))])


disp(['average predicted probability of folding into experimental topology (non-pseudoknotted) = ',num2str(mean(probTopology(1:numNonPseudoknotSequences)))])
disp(['average predicted probability of folding into experimental topology (pseudoknotted) = ',num2str(mean(probTopology(numNonPseudoknotSequences+1:end)))])

figure;
histogram(probTopology(1:numNonPseudoknotSequences),0:0.1:1)
xlabel('predicted probability of folding into experimental topology')
ylabel('frequency')
title('non-pseudoknotted structures')
set(gca,'fontsize',16)
figure;
histogram(probTopology(numNonPseudoknotSequences+1:end),0:0.1:1)
xlabel('predicted probability of folding into experimental topology')
ylabel('frequency')
title('pseudoknotted structures')
set(gca,'fontsize',16)
% 
% disp(['average base pairs part of topologies which included parallel loops (non-pseudoknotted) = ',num2str(mean(avgPredictedParallelPerBase(1:54,:)))])
% disp(['average base pairs part of topologies which included parallel loops (pseudoknotted) = ',num2str(mean(avgPredictedParallelPerBase(55:197,:)))])
% 

%disp(['average correctly predicted per-base topology = ',num2str(mean(fractionPredictedTopology))])
disp(['average correctly predicted non-pseudoknotted per-base topology = ',num2str(mean(fractionPredictedTopology(1:numNonPseudoknotSequences,:)))])
disp(['average correctly predicted pseudoknotted per-base topology = ',num2str(mean(fractionPredictedTopology(numNonPseudoknotSequences+1:end,:)))])


%disp(['average predicted topologies included parallel loops = ',num2str(mean(predictedParallel))])
% disp(['average predicted topologies included parallel loops (non-pseudoknotted) = ',num2str(mean(predictedParallel(1:numNonPseudoknotSequences,:)))])
% disp(['average predicted topologies included parallel loops (pseudoknotted) = ',num2str(mean(predictedParallel(numNonPseudoknotSequences+1:end,:)))])
fracPredTopPseud = fractionPredictedTopology(numNonPseudoknotSequences+1:end,:)';
mean(fracPredTopPseud(~unknown(numNonPseudoknotSequences+1:end)))



%% Generate Fig. 5
figure;
hold on
%plot the probability of pseudoknot
probPseudoknot(probPseudoknot<=1e-9) = 2e-10;
nonPseudoProbPseudo = probPseudoknot(1:numNonPseudoknotSequences);
pseudoProbPseudo = probPseudoknot(numNonPseudoknotSequences+1:end);
sortedNonPseudoProbPseudo = sort(probPseudoknot(1:numNonPseudoknotSequences));
sortedPseudoProbPseudo = sort(probPseudoknot(numNonPseudoknotSequences+1:end));
a=scatter(1:numNonPseudoknotSequences,sortedNonPseudoProbPseudo,[],rgb('bright blue'),'o','MarkerEdgeAlpha',0.5);
gapPseudo = 20;
b=scatter(gapPseudo+numNonPseudoknotSequences+1:numNonPseudoknotSequences/numPseudoknotSequences:gapPseudo+numNonPseudoknotSequences*2+1-numNonPseudoknotSequences/numPseudoknotSequences,sortedPseudoProbPseudo,[],rgb('salmon'),'o','MarkerEdgeAlpha',0.5);
set(gca,'yscale','log','fontsize',16)
set(gca,'xtick',[numNonPseudoknotSequences/2,gapPseudo+numNonPseudoknotSequences*3/2],'xticklabels',{'\color[rgb]{0.0039, 0.3961, 0.9882} Non-Pseudoknotted','\color[rgb]{1.0000, 0.4745, 0.4235} Pseudoknotted'})
ax = gca;
ax.XAxis.TickLength = [0.00 0.0];
ylabel('Probability of forming pseudoknot')
%legend([a,b],{'Non-Pseudoknotted','Pseudoknotted'})

% % get the current tick labeks
% ticklabels = get(gca,'XTickLabel');
% % prepend a color for each tick label
% ticklabels_new = cell(size(ticklabels));
% for i = 1:length(ticklabels)
%     ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
% end
% % set the tick labels
% set(gca, 'YTickLabel', ticklabels_new);