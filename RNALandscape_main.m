function[sortedProbs, indexSortedProbs,sortedFE, indexSortedFE,graphProbs,weightMatrixList,numBondsMatrixList,startAndPermuTime,checkFxnTime,totalTime,STable,possiblePermutations,numRegions,expFEeff] = ...
    RNALandscape_main(sequences,T,vs,duplexEntropyPenaltyInKB,minBPInRegion,numParForLoops,pairwiseEnergies,storeGraphs,makingFigures,...
    printProgressUpdate,allowParallelStrands, allowPseudoknots, minNumCompatible, substems)
%This function takes as input a cell array of sequences (ex:
%{'AUGCGC','GCGCAU'}) and outputs the entire free energy landscape of these sequences interacting.
%These sequences can be DNA or RNA -- the code recognizes DNA sequences if
%they have a T. The code can be run with only one sequence e.g.
%{'GCGAAAACGC'}. 

%T sets the temperature, which affects the free energy calculation FE=E-TS.
%It is measured in Kelvin

%vs is the volume (in ntds^3) within which two nucleotides (ntds) can bond.

%duplexEntropyPenaltyInKB only affects runs with more than one sequence --
%it defines the entropy penalty (in units of Boltzmann's constant) of two
%strands being bound. This will in general depend on the size of the container and the
%relative concentrations of the two sequences. Future iterations of this
%code will take the latter into account as well.

%minBPInRegion defines the minimum length of a stem (it is equal to the
%parameter m in the paper). It is an integer > 0. min value it can take is 
%1; Pipas & McMahon set it to 3. 

%numParForLoops affects the run itself. Set it to 0 to use a for loop
%rather than a parfor loop (parallel computing) to calculate free energies.
%Otherwise set it to any integer >0 to define how many iterations of parfor
%will be run.

%pairwiseEnergies is a boolean. It is true if you want to use the Turner
%energy parameters. Otherwise, we have defined a simple model that uses
%only the energies of single base pair (bp) interactions based on work by
%Tristan Cragnolini.

%storeGraphs is a boolean. Set it to true if you'd like to store the
%topology of each structure and all the topologies. The code runs faster if
%it is set to false, but this disallows coarse-graining by topology.

%makingFigures defines whether you'd like the code to output figures at the
%end of the calculation. It is a boolean

%printProgressUpdate is a boolean. Set it to true to print messages
%updating you on how the algorithm is progressing (e.g. how long it's spent
%on each subfunction).

%allowParallelStrands defines whether to allow parallel stems in the
%computation. If pairwiseEnergies is true, it should be set to false;
%otherwise, you can set it to true.

%allowPseudoknots defines whether to allow pseudoknots in the computation

%minNumCompatible defines the minimum number of compatible stems that can
%comprise a structure. Set it to 0 to get the most general results.
%Otherwise set it to an integer >0.

%substems defines whether to allow all possible substems. If a stem of
%length l is found to be possible, but minBPInRegion = m<l, subsets of this
%stem of any length between l & m are also allowed. If substems is set to 'all', all of
%these possible stems will be considered as well. If substems is set to an
%integer >= 0, the code will only consider subsets of the full stem of
%length l-substems.

%----------------------------------------------------------------------
%sortedProbs is a vector which gives the probability for each structure
%forming. It is sorted so that sortedProbs(1)>=sortedProbs(2), etc.
%indexSortedProbs maps from each element of sortedProbs which element of
%possiblePermutations (i.e. which structure) it refers to.

%sortedFE is similar -- it gives the free energy of each structure, in a
%sorted fashion. indexSortedFE is the same

%graphProbs is also a vector which gives for each topology its probability of forming. 
%Topology i is defined by weightMatrixList{i} and numBondsMatrixList{i}
%which are two matrices of size numNodes x numNodes (i.e. the number of
%nodes in the graph defining that topology). weightMatrixList{i}(j,k) gives
%the "weight" of all the bonds between nodes j and k. Bonds defining stems
%(i.e. double-stranded regions) count as weight 2; bonds defining
%single-stranded regions count as weight 1. numBondsMatrixList{i}(j,k) is
%similar but simply gives the number of bonds between nodes j and k. For
%example if we have an open-net 0, nodes 1 & 2 are connected by one double
%bond and one single bond. Thus the 2x2 weight matrix has 3's on the
%off-diagonal, and 0's on the diagonal and the 2x2 numBonds matrix has 2's
%on the off-diagonal and 0's on the diagonal.

%startAndPermuTime,checkFxnTime, and totalTime give the time it took to
%perform each of these calculations (the start and permu functions
%together, i.e. enumerating all possible structures; the check function,
%i.e. calculating all the free energies; and the total time, respectively).
%The names of the functions comes from Pipas & McMahon's excellent paper.

%STable is a numStems x 3 cell array. The first column is just the number
%of the stem. The second column is the stem written in "structure" format.
%This is a vector of length 2x(length of stem) which defines the base pairs
%in the stem. For example, [1 2 3 51 50 49] means ntd 1 is bonded to 51, 2
%to 50, 3 to 49. The third column of STable is the same data in another
%format: a (length of stem) x 2 matrix, where each row defines a base pair.
%The same stem would be given by [1 51; 2 50; 3 49].

%possiblePermutations defines the structures. It is a numStructures x 1
%cell array. Each element of the cell array is a vector which is a list of
%which stems (from STable) are found in that given structure. For example,
%if possiblePermutations{i} = [1 3 20] that means that the i'th structure
%has stems 1, 3, and 20 (defined as the first, third, and twentieth rows of
%STable).

%numRegions is the number of possible stems, given as the number of rows of
%STable. Note: I use the word "region" to mean "stem" throughout.

%expFEeff only applies if you are considering multiple strands. It is the
%probability of the strands being bound vs separate and is given as a
%vector

%%
if nargin < 20
    substems = 'all';
end
if nargin < 19
    minNumCompatible = 0;
end
if nargin < 18
    allowPseudoknots = true; %should we allow pseudoknots to form? 
end
if nargin < 17
    allowParallelStrands = false;
end

sequences = sort(sequences); %sorts alphabetically

[sequence,sequenceInNumbers,numSequences,numBP,~,numBP1,linkerPos,~,numBP2,~,numBP3] = ...
    multipleStrandsSetup_Fxn(sequences); %deals with multiple sequences
%by concatenating them into one long sequence with a linker of O's which
%can't bind to anything in between. We don't consider the entropy
%contribution of these O's, of course.

%% START function
startTime = tic;
    
omega = ones(1,numBP+1); %omega(s) term from paper. Set to 1 to get Delta S as described in paper.

kB = 0.0019872; %units of kcal/(mol*Kelvin)

minBPInHairpin = 3; %minimum number of nucleotides in a hairpin.

disp(['sequence = ',sequence]);
if printProgressUpdate
    disp(['minBPInRegion = ',num2str(minBPInRegion)]);
end
diary off %need this to write everything till this point to the diary
diary on
%%

[numRegions, STable] = createSTable_Fxn(sequenceInNumbers,minBPInRegion,minBPInHairpin,allowParallelStrands,substems);

if printProgressUpdate
    disp(['numRegions = ',num2str(numRegions)]);
    diary off %need this to write everything till this point to the diary
    diary on
end

C = zeros(numRegions,numRegions); %matrix describing which regions are compatible
%if regions i and j are compatible, Cij = 1; else, Cij = 0.
for i = 1:numRegions
    for j = 1:numRegions
        if i ==j
            C(i,j) = 1; %diagonal elements of C matrix are set to 1.
        elseif all(ismember(STable{i,2},STable{j,2}) == 0) %checks if regions i and j share any of the same bases
            C(i,j) = 1;
            if ~allowPseudoknots
                %-----------------next part lets us disallow pseudoknots
                %From Mauri et al. 2005 (substitute i,j,k,l for a,b,c,d):
                %If ntd a is paired with ntd b>a, and c>a is paired with d>c,
                %then either a<b<c<d or a<c<d<b, but a<c<b<d is a pseudoknot.
                
                a = STable{i,2}(1);
                b = STable{i,2}(1+length(STable{i,2})/2);
                c = STable{j,2}(1);
                d = STable{j,2}(1+length(STable{j,2})/2);
                
                if c<a %switch around labels so we have c>a
                    a = STable{j,2}(1);
                    b = STable{j,2}(end);
                    c = STable{i,2}(1);
                    d = STable{i,2}(end);
                end
                
                allOnOneStrand = true; %are all ntds here on the same strand?
                for k  = 1:length(linkerPos)
                    if sign(a - linkerPos{k}(1)) == sign(b - linkerPos{k}(1)) && ...
                            sign(a - linkerPos{k}(1)) == sign(c - linkerPos{k}(1)) &&...
                            sign(a - linkerPos{k}(1)) == sign(d - linkerPos{k}(1))
                        %allOnOneStrand = true;
                    else
                        allOnOneStrand = false;
                    end
                end
                
                if (a<c && c<b && b<d) && allOnOneStrand
                    C(i,j) = 0;
                end
                %---------------------
            end
        end
    end
end

%But: we don't want to allow two subregions to be compatible if including
%both of them just looks like one longer region. For example, a subregion
%of ntd 1 bonded to 10 shouldn't be compatible with a subregion of ntd 2
%bonded to 9, even if they came from different full regions initially. 

C2 = C; %C2 is the corrected version of the C matrix which we will now use.

for i = 1:numRegions
    for j = 1:numRegions
        if C2(i,j) && i ~=j %not worth going through this code if C2(i,j) ==0
            regionI = STable{i,2};
            lenRegionI = length(regionI);
            regionJ = STable{j,2};
            lenRegionJ = length(regionJ);
            
            if STable{j,2}(1)<STable{i,2}(1) %switch label of the regions so regionI starts before regionJ
                regionSwap = regionJ;
                regionJ = regionI;
                regionI = regionSwap;
                lenRegionI = length(regionI);
                lenRegionJ = length(regionJ);
            end
            combinedRegion = [regionI(1:lenRegionI/2), regionJ(1:lenRegionJ/2),...
                regionI(lenRegionI/2+1:lenRegionI),regionJ(lenRegionJ/2+1:lenRegionJ)];
            lenCombined = length(combinedRegion);
            isExistingRegion = false; %we want to know if the combined region is an already-existing region.
            for k = 1:numRegions 
                if isequal(combinedRegion,STable{k,2}) 
                    isExistingRegion = true;
                    break
                end
            end
            
            if isExistingRegion
                C2(i,j) = 0;
                C2(j,i) = 0;
            
            %or if it's a region that could have formed and didn't, get rid
            %of it (this could be e.g. because it leads to a hairpin that's
            %too short)
            elseif isequal(combinedRegion(1:lenCombined/2),combinedRegion(1):combinedRegion(lenCombined/2)) &&...
                    (isequal(combinedRegion(lenCombined/2+1:end),combinedRegion(lenCombined/2+1):combinedRegion(end)) ||...
                    isequal(combinedRegion(lenCombined/2+1:end),combinedRegion(end):-1:combinedRegion(lenCombined/2+1)))
                C2(i,j) = 0;
                C2(j,i) = 0;
            end
        end
    end
end

%%

C3 = []; %matrix specifying threeway compatibility
%regions i,j, and k are compatible if they are pairwise compatible and if
%they do not form a pseudoknot of higher order than what we can calculate
%(i.e. of order 3 or above).

if numSequences >1
    numRegionsNecessaryforC3 = 10000; %don't make a C3 matrix if we have numSequences>1
else 
    numRegionsNecessaryforC3 = 0;%how many regions do we need to have for using C3 to be worth it?
end
if numRegions>=numRegionsNecessaryforC3
    C3 = false(numRegions,numRegions,numRegions);
    try %there is a one-time cost (not once per sequence but one time total) to run the code ennumerateDisallowedPermutationsWithPseudoknots(3). 
        %So we run it once, and any time after that we run code we can just
        %load it.
        
        %load the results of ennumerateDisallowedPermutations(3) since it takes
        %over ten minutes to run
        newData1 = load('-mat', '~/disallowedPermutationsWithPseudoknots3.mat');
        vars = fieldnames(newData1);  %for some reason, the data is saved as a struct
        threewayDisallowedPermutations = newData1.(vars{1});
    catch
        threewayDisallowedPermutations = ennumerateDisallowedPermutationsWithPseudoknots(3);
        save('~/disallowedPermutationsWithPseudoknots3.mat','threewayDisallowedPermutations');
    end


    for i = 1:numRegions
        for j = i+1:numRegions
            if C2(i,j)
            for k = j+1:numRegions
                if C2(i,k) && C2(j,k) %they need to be pairwise compatible
                    regionI = STable{i,2};
                    regionJ = STable{j,2};
                    regionK = STable{k,2};
                    a = regionI(1);
                    b = regionI(length(regionI)/2+1);
                    c = regionJ(1);
                    d = regionJ(length(regionJ)/2+1);
                    e = regionK(1);
                    f = regionK(length(regionK)/2+1);
                    
                    %given a,b,c,d,e,f, (more importantly, their relations
                    %to each other, like which is greater than which) we
                    %can already rule out a structure in some cases as
                    %being a pseudoknot of too high order for us to
                    %compute.

                    %check if the three regions are otherwise not compatible
                    %put the regions into order such that a<c<e
                    [sortThree,index] = sortrows([a b; c d; e f]);
                    a = sortThree(1); c = sortThree(2); e = sortThree(3);
                    b = sortThree(4); d = sortThree(5); f = sortThree(6);
                    regionIKList = [{regionI},{regionJ},{regionK}];
                    regionI = regionIKList{index(1)}; regionJ = regionIKList{index(2)}; regionK = regionIKList{index(3)}; 

                    abPseudoknot = false;
                    abNonPseudoknot = false;
                    if length(regionI)>2 
                        if regionI(length(regionI)/2+2) - regionI(length(regionI)/2+1) == 1
                            abPseudoknot = true;
                        else
                            abNonPseudoknot = true;
                        end
                    end
                    cdPseudoknot = false;
                    cdNonPseudoknot = false;
                    if length(regionJ)>2 
                        if regionJ(length(regionJ)/2+2) - regionJ(length(regionJ)/2+1) == 1
                            cdPseudoknot = true;
                        else
                            cdNonPseudoknot = true;
                        end
                    end
                    efPseudoknot = false;
                    efNonPseudoknot = false;
                    if length(regionK)>2  
                        if regionK(length(regionK)/2+2) - regionK(length(regionK)/2+1) == 1
                            efPseudoknot = true;
                        else
                            efNonPseudoknot = true;
                        end
                    end


                    [~,index] = sort([a,b,c,d,e,f]);
                    checkDisallowedPermu = [index, abPseudoknot,abNonPseudoknot,...
                        cdPseudoknot,cdNonPseudoknot,efPseudoknot,efNonPseudoknot];
                    disallowedPseudoknot = false;
                    if any(all((threewayDisallowedPermutations == checkDisallowedPermu)'))
                        %this line (above) is equivalent to: ismember(checkDisallowedPermu,threewayDisallowedPermutations,'rows')
                        disallowedPseudoknot = true;
                    end

                    if ~disallowedPseudoknot
                        C3(i,j,k) = true;
                        C3(i,k,j) = true;
                        C3(j,i,k) = true;
                        C3(j,k,i) = true;
                        C3(k,i,j) = true;
                        C3(k,j,i) = true;
                    end
                end
            end
            end
        end
    end
end

%now make a matrix of fourway compatibility in the same way.
if numSequences >1
    numRegionsNecessaryforC4 = 10000; %don't make a C4 matrix if we have numSequences>1
else 
    numRegionsNecessaryforC4 = 0; %how many regions do we need to have for using C4 to be worth it?
end

if numRegions>=numRegionsNecessaryforC4 
    C4 = false(numRegions,numRegions,numRegions,numRegions);
    counter = 1;
    prevCounter = 0;
    %regions i,j, k, and l are compatible if they are pairwise & threewise compatible and if
    %they do not form a pseudoknot of higher order than what we can calculate
    %(i.e. of order 3 or above).
    try 
        %load the results of ennumerateDisallowedPermutations(4) since it takes
        %over ten minutes to run
        newData1 = load('-mat', '~/disallowedPermutationsWithPseudoknots4.mat');
        vars = fieldnames(newData1); %for some reason, the data is saved as a struct
        fourwayDisallowedPermutations = newData1.(vars{1});
    catch 
        fourwayDisallowedPermutations =ennumerateDisallowedPermutationsWithPseudoknots(4);
        save('~/disallowedPermutationsWithPseudoknots4.mat','fourwayDisallowedPermutations');
    end

    startTimeWithoutC4 = toc(startTime);
    if printProgressUpdate
        disp(['Time until using C4 = ',num2str(startTimeWithoutC4)]);
        diary off %need this to write everything till this point to the diary
        diary on
    end
    

    for i = 1:numRegions
        for j = i+1:numRegions
            if C2(i,j)
            for k = j+1:numRegions
                if C3(i,j,k) %they need to be pairwise compatible in order to be three-way compatible
                for l = k+1:numRegions
                    if C3(i,j,l) && C3(i,k,l) && C3(j,k,l) %they need to be pairwise compatible and threewise compatible
                        regionI = STable{i,2};
                        regionJ = STable{j,2};
                        regionK = STable{k,2};
                        regionL = STable{l,2};
                        a = regionI(1);
                        b = regionI(length(regionI)/2+1);
                        c = regionJ(1);
                        d = regionJ(length(regionJ)/2+1);
                        e = regionK(1);
                        f = regionK(length(regionK)/2+1);
                        q = regionL(1);
                        r = regionL(length(regionL)/2+1);

                        %put the regions into order such that a<c<e<g
                        [sortFour,index] = sortrows([a b; c d; e f; q r]);
                        a = sortFour(1); c = sortFour(2); e = sortFour(3);q= sortFour(4);
                        b = sortFour(5); d = sortFour(6); f = sortFour(7);r= sortFour(8);
                        regionILList = [{regionI},{regionJ},{regionK},{regionL}];
                        regionI = regionILList{index(1)}; regionJ = regionILList{index(2)};
                        regionK = regionILList{index(3)}; regionL = regionILList{index(4)};

                        abPseudoknot = false;
                        abNonPseudoknot = false;
                        if length(regionI)>2
                            if regionI(length(regionI)/2+2) - regionI(length(regionI)/2+1) == 1
                                abPseudoknot = true;
                            else
                                abNonPseudoknot = true;
                            end
                        end
                        cdPseudoknot = false;
                        cdNonPseudoknot = false;
                        if length(regionJ)>2
                            if regionJ(length(regionJ)/2+2) - regionJ(length(regionJ)/2+1) == 1
                                cdPseudoknot = true;
                            else
                                cdNonPseudoknot = true;
                            end
                        end
                        efPseudoknot = false;
                        efNonPseudoknot = false;
                        if length(regionK)>2
                            if regionK(length(regionK)/2+2) - regionK(length(regionK)/2+1) == 1
                                efPseudoknot = true;
                            else
                                efNonPseudoknot = true;
                            end
                        end
                        qrPseudoknot = false;
                        qrNonPseudoknot = false;
                        if length(regionL)>2
                            if regionL(length(regionL)/2+2) - regionL(length(regionL)/2+1) == 1
                                qrPseudoknot = true;
                            else
                                qrNonPseudoknot = true;
                            end
                        end



                        [~,index] = sort([a,b,c,d,e,f,q,r]);
                        checkDisallowedPermu = [index, abPseudoknot,abNonPseudoknot,...
                        cdPseudoknot,cdNonPseudoknot,efPseudoknot,efNonPseudoknot,qrPseudoknot,qrNonPseudoknot];
                        disallowedPseudoknot = false;
                        if any(all((fourwayDisallowedPermutations == checkDisallowedPermu)'))
                            %this line (above) is equivalent to: ismember(checkDisallowedPermu,fourwayDisallowedPermutations,'rows')
                            disallowedPseudoknot = true;
                        end
                        if ~disallowedPseudoknot
                            counter = counter + 1;
                            C4(i,j,k,l) = true; C4(i,j,l,k) = true; C4(i,k,j,l) = true; C4(i,k,l,j) = true; C4(i,l,j,k) = true; C4(i,l,k,j) = true;
                            C4(j,i,k,l) = true; C4(j,i,l,k) = true; C4(j,k,i,l) = true; C4(j,k,l,i) = true; C4(j,l,i,k) = true; C4(j,l,k,i) = true;
                            C4(k,i,j,l) = true; C4(k,i,l,j) = true; C4(k,j,i,l) = true; C4(k,j,l,i) = true; C4(k,l,i,j) = true; C4(k,l,j,i) = true;
                            C4(l,i,j,k) = true; C4(l,i,k,j) = true; C4(l,j,i,k) = true; C4(l,j,k,i) = true; C4(l,k,i,j) = true; C4(l,k,j,i) = true;
                        end
                        if rem(counter,500000 )==0 && printProgressUpdate && counter ~= prevCounter
                            disp(['Setting up C4: i = ',num2str(i)])
                            disp(['length c4 indices = ',num2str(counter)])
                            prevCounter = counter; %so we don't save the same thing multiple times
                            diary off %need this to write everything till this point to the diary
                            diary on
                        end
                    end
                end
                end
            end
            end
        end
    end
else
    C4 = []; %to save space and so that if numRegions<numRegionsNecessaryforC4 we can check if C4 is empty to know whether to use it
end

startFxnTime = toc(startTime);
if printProgressUpdate
    if numRegions>=numRegionsNecessaryforC4 
        disp(['length c4 indices = ',num2str(counter)])
    else
        disp('numRegions < numRegionsNecessaryforC4')
    end
    disp(['Time for START function = ',num2str(startFxnTime)]);
    diary off %need this to write everything till this point to the diary
    diary on
end
%% PERMU function
permuTime = tic;
possiblePermutations = PERMU_Fxn(numRegions,C2,C3,C4,minNumCompatible,...
    numRegionsNecessaryforC3,numRegionsNecessaryforC4,printProgressUpdate);
permuFxnTime = toc(permuTime);
startAndPermuTime = toc(startTime);
if printProgressUpdate
    disp(['Time for PERMU function = ',num2str(permuFxnTime)]);
    diary off %need this to write everything till this point to the diary
    diary on
end


clear('C3', 'C4') %to save space and since they can relatively quickly be regenerated


%% CHECK function
%Assign free energy to each possible structure

numStructures = length(possiblePermutations);
disp(['numStructures = ',num2str(numStructures)]);
diary off %need this to write everything till this point to the diary
diary on
allFreeEnergies = zeros(1,numStructures); %table of all deltaGs of the structures we find.
allBondEnergies = zeros(1,numStructures);
allBondEntropies = zeros(1,numStructures);
allLoopEntropies = zeros(1,numStructures);
if numSequences>1
    numStrandsLinkedVec = zeros(1,numStructures); %for each structure, how many strands are bonded to one another in that structure?
end
if storeGraphs
    %list of distinct graphs our structures correspond to
    weightMatrixList = cell(0); %list of distinct weight matrix and numBond matrix pairs our structures correspond to
    %(in other words, there might be repeats in either list, but these
    %still correspond to distinct structures because the same weight matrix
    %with a different numBonds matrix is a different graph).
    numBondsMatrixList = cell(0);
    allGraphs = zeros(1,numStructures); %which of the graphs in graphList does each structure correspond to?
    if numParForLoops > 0
        listOfAllNewWeightMatrices = cell(1,numStructures);
        listOfAllNewNumBondsMatrices = cell(1,numStructures);
    end
else
    allGraphs = 0;
    weightMatrixList = 0; numBondsMatrixList = 0;
end



checkTime = tic;

whichStructLeftToCompute = 1:numStructures; %structures whose FE we have yet to compute

allPerms = cell(1,4); %Used in entropy calculation and in determining graph isomorphisms. run this once so that we don't need to run it ever again
for i = 1:4
    allPerms{i} = perms(1:i);
end

if numParForLoops > 0
    whichStructSeries = cell(1,numParForLoops);
    numStructuresToCompute = length(whichStructLeftToCompute);
    for i = 1:numParForLoops
        whichStructSeries{i} = whichStructLeftToCompute(floor(numStructuresToCompute*(i-1)/numParForLoops)+1:...
            floor(numStructuresToCompute*i/numParForLoops));
    end

    for parForCounter = 1:numParForLoops
        disp(['parForCounter = ',num2str(parForCounter)])
        
%         if numRegions > 198 && parForCounter > 1
%             %get rid of most of graphs if they don't have enough probability of
%             %occuring to save on space
%             probabilities = exp(-allFreeEnergies./(kB*T));
%             partitionFxn = sum(probabilities);
%             probabilities = probabilities./partitionFxn;
%             numGraphs = length(weightMatrixList);
%             graphProbs = zeros(1,numGraphs);
%             for i = 1:whichStruct
%                 if allGraphs(i)>0
%                     graphProbs(allGraphs(i)) = graphProbs(allGraphs(i)) + ...
%                         probabilities(i);
%                 end
%             end
%             graphProbCutoff = 1e-7;
%             weightMatrixList(graphProbs<graphProbCutoff) = [];
%             numBondsMatrixList(graphProbs<graphProbCutoff) = [];
%             %now need to update allGraphs
%             for i = 1:whichStruct
%                 index = allGraphs(i);
%                 numDeleted = sum(graphProbs(1:index)<graphProbCutoff);
%                 allGraphs(i) = max(0,allGraphs(i) - numDeleted);
%             end
%             
%         else
%             graphProbCutoff = -1;
%         end
        
        diary off %need this to write everything till this point to the diary
        diary on
        
        parfor whichStruct = whichStructSeries{parForCounter} %#ok<PFRNG> %which structure do you want to calculate free energies for?
            structure = cell(1,length(possiblePermutations{whichStruct})); %list of bonded bps in structure    
            %%
            strucMatrix = zeros(numBP); %connectivity matrix of structure (strucMatrix_i,j = 1 if ntds i and j are bonded)

            %set up structure and strucMatrix
            for i = 1:length(structure)
                structure{i} = STable{possiblePermutations{whichStruct}(i),2};
                numBonds = length(structure{i})/2;
                for j = 1:numBonds
                    k = structure{i}(j);
                    l = structure{i}(j+numBonds);
                    strucMatrix(k,l) = 1;
                    strucMatrix(l,k) = 1;
                end
            end

            %calculate energy of bonds
            if pairwiseEnergies
                [bondEnergy,bondEntropy] = pairwiseBondFreeEnergies_Fxn(structure,sequenceInNumbers,linkerPos,T);
            else
                bondEnergy = singleBondEnergies(structure,sequenceInNumbers);
                bondEntropy = 0;
            end
            allBondEnergies(whichStruct) = bondEnergy;
            allBondEntropies(whichStruct) = bondEntropy;
            
            bondFreeEnergy = bondEnergy -T*bondEntropy;

            [loopEntropy,initialWeightMatrix,initialNumBondsMatrix] = ...
                calculateEntropy_Fxn(structure,strucMatrix,numBP,omega,minBPInHairpin,minBPInRegion,whichStruct,allPerms,vs,linkerPos,allowParallelStrands);
            
            if storeGraphs 
                listOfAllNewWeightMatrices{whichStruct} = initialWeightMatrix;
                listOfAllNewNumBondsMatrices{whichStruct} = initialNumBondsMatrix;                
            end
            
            %calculate the number of strands linked together in this structure
            if numSequences > 1
                numStrandsLinked = 1;
                if any(any(strucMatrix(1:numBP1,linkerPos{1}(end)+1:linkerPos{1}(end)+numBP2)))
                    numStrandsLinked = numStrandsLinked + 1;
                end
                if numSequences == 3 
                    if any(any(strucMatrix(1:numBP1,linkerPos{2}(end)+1:linkerPos{2}(end)+numBP3)))
                        numStrandsLinked = numStrandsLinked + 1;
                    end
                    if numStrandsLinked < 3 %then we still haven't determined if all 3 sequences are bound to one another
                        if any(any(strucMatrix(linkerPos{1}(end)+1:linkerPos{1}(end)+numBP2,linkerPos{2}(end)+1:linkerPos{2}(end)+numBP3)))
                            numStrandsLinked = numStrandsLinked + 1;
                        end
                    end
                end
            else
                numStrandsLinked = 1;
            end
            if numSequences > 1
                numStrandsLinkedVec(whichStruct) = numStrandsLinked;
            end
  
            loopEntropy = loopEntropy - kB*duplexEntropyPenaltyInKB*(numStrandsLinked-1); %give a translational entropy penalty if strands are bonded to one another
            
            
            allLoopEntropies(whichStruct) = loopEntropy;

            loopFreeEnergy = -T*loopEntropy;
            allFreeEnergies(whichStruct) = bondFreeEnergy+loopFreeEnergy;

        end
        whichStruct = whichStructSeries{parForCounter}(end);
        oldLengthWeightMatrixList = length(weightMatrixList); 
        if oldLengthWeightMatrixList == 0 %as it is for first time through (at least if the code hadn't been started before)
            oldLengthWeightMatrixList = 1;
        end
        if storeGraphs
            for whichStruct = whichStructSeries{parForCounter}
                %for each structure, check if it makes a new graph topology for this sequence
                newGraph = true;
                initialWeightMatrix = listOfAllNewWeightMatrices{whichStruct};
                if ~isempty(initialWeightMatrix)
                    initialNumBondsMatrix = listOfAllNewNumBondsMatrices{whichStruct};
                    if isempty(weightMatrixList) %meaning we need to initialize allGraphs, and the two matrix lists
                        allGraphs(whichStruct) = 1;
                        weightMatrixList = {initialWeightMatrix};
                        numBondsMatrixList = {initialNumBondsMatrix};
                    else
                        for i = length(weightMatrixList):-1:1 %if in parForCounter >1 we have a new weightMatrix,
                            %chances are that it's replicated in other
                            %structures.
                            if myPseudoIsIsomorphic(initialWeightMatrix,...
                                    weightMatrixList{i},initialNumBondsMatrix,numBondsMatrixList{i})
                                allGraphs(whichStruct) = i;
                                newGraph = false;
                                break
                            end
                        end
                        if newGraph
                            weightMatrixList{length(weightMatrixList)+1} = initialWeightMatrix;
                            allGraphs(whichStruct) = length(weightMatrixList);
                            numBondsMatrixList{length(numBondsMatrixList)+1} = initialNumBondsMatrix;
                        end
                    end
                    listOfAllNewWeightMatrices{whichStruct} = {}; %to free up space
                    listOfAllNewNumBondsMatrices{whichStruct} = {};
                end
            end
        end
        disp(['time elapsed since start of check = ',num2str(toc(checkTime))])
    end

else %if we don't want to use parallel for loops
    for whichStruct = whichStructLeftToCompute %which structure do you want to calculate free energies for?

        if rem(whichStruct,100000)==0  && ~redoEnergies
            disp(whichStruct)
            disp(['time elapsed since start of check = ',num2str(toc(checkTime))])
            diary off %need this to write everything till this point to the diary
            diary on
        end
        structure = cell(1,length(possiblePermutations{whichStruct})); %list of bonded bps in structure    
        %%
        strucMatrix = zeros(numBP); %connectivity matrix of structure (strucMatrix_i,j = 1 if ntds i and j are bonded)

        %set up structure and strucMatrix
        for i = 1:length(structure)
            structure{i} = STable{possiblePermutations{whichStruct}(i),2};
            numBonds = length(structure{i})/2;
            for j = 1:numBonds
                k = structure{i}(j);
                l = structure{i}(j+numBonds);
                strucMatrix(k,l) = 1;
                strucMatrix(l,k) = 1;
            end
        end

        if pairwiseEnergies
            [bondEnergy,bondEntropy] = pairwiseBondFreeEnergies_Fxn(structure,sequenceInNumbers,linkerPos,T);
        else
            bondEnergy = singleBondEnergies(structure,sequenceInNumbers);
            bondEntropy = 0;
        end
        
        allBondEnergies(whichStruct) = bondEnergy;
        allBondEntropies(whichStruct) = bondEntropy;
            
        bondFreeEnergy = bondEnergy -T*bondEntropy;
            
        
        %calculate the number of strands linked together in this structure
        if numSequences > 1
            numStrandsLinked = 1;
            if any(any(strucMatrix(1:numBP1,linkerPos{1}(end)+1:linkerPos{1}(end)+numBP2)))
                numStrandsLinked = numStrandsLinked + 1;
            end
            if numSequences == 3 
                if any(any(strucMatrix(1:numBP1,linkerPos{2}(end)+1:linkerPos{2}(end)+numBP3)))
                    numStrandsLinked = numStrandsLinked + 1;
                end
                if numStrandsLinked < 3 %then we still haven't determined if all 3 sequences are bound to one another
                    if any(any(strucMatrix(linkerPos{1}(end)+1:linkerPos{1}(end)+numBP2,linkerPos{2}(end)+1:linkerPos{2}(end)+numBP3)))
                        numStrandsLinked = numStrandsLinked + 1;
                    end
                end
            end
        else 
            numStrandsLinked = 1;
        end
        
        
        [loopEntropy,initialWeightMatrix,initialNumBondsMatrix] = ...
            calculateEntropy_Fxn(structure,strucMatrix,numBP,omega,minBPInHairpin,minBPInRegion,whichStruct,allPerms,vs,linkerPos,allowParallelStrands);
        
        loopEntropy = loopEntropy - kB*duplexEntropyPenaltyInKB*(numStrandsLinked-1); %give a translational entropy penalty if strands are bonded to one another
        
        if storeGraphs
            %check if this is a new graph topology for this sequence
            newGraph = true;
            if isempty(weightMatrixList)
                allGraphs(whichStruct) = 1;
                weightMatrixList = {initialWeightMatrix};
                numBondsMatrixList = {initialNumBondsMatrix};
            else
                for i = 1:length(weightMatrixList)
                    if myPseudoIsIsomorphic(initialWeightMatrix,...
                            weightMatrixList{i},initialNumBondsMatrix,numBondsMatrixList{i})
                        allGraphs(whichStruct) = i;
                        newGraph = false;
                        break
                    end
                end
                if newGraph
                    weightMatrixList{length(weightMatrixList)+1} = initialWeightMatrix;
                    allGraphs(whichStruct) = length(weightMatrixList);
                    numBondsMatrixList{length(numBondsMatrixList)+1} = initialNumBondsMatrix;
                end
            end
        end
        
        allLoopEntropies(whichStruct) = loopEntropy;
        
        loopFreeEnergy = -T*loopEntropy;
        allFreeEnergies(whichStruct) = bondFreeEnergy-T*allLoopEntropies(whichStruct); 
        if numSequences > 1
            numStrandsLinkedVec(whichStruct) = numStrandsLinked;
        end
    end
end

expFEeff = zeros(1,numSequences); %what is the effective free energy of having each strand separate (index 1), two strands bonded (index 2) or three strands bonded (index 3)?
%effective free energy is defined here as = -kB*T*log(\sum(exp(-FE/kB*T)))
%where the sum is over all structures with 1) each strand separate; 2) two
%strands bonded; or 3) all three strands bonded.
if numSequences>1
    for whichStruct = 1:numStructures
        numStrandsLinked = numStrandsLinkedVec(whichStruct);
        expFEeff(numStrandsLinked) = expFEeff(numStrandsLinked) + exp(-allFreeEnergies(whichStruct)/(kB*T));
    end
end

if numSequences > 1
    disp(['probability of having each sequence separate = ', num2str(expFEeff(1)/sum(expFEeff))])
    if numSequences == 2
        disp(['probability of having sequences bonded = ', num2str(expFEeff(2)/sum(expFEeff))])
    elseif numSequences == 3
        disp(['probability of having 2 sequences bonded = ', num2str(expFEeff(2)/sum(expFEeff))])
        disp(['probability of having each all 3 sequences bonded = ', num2str(expFEeff(3)/sum(expFEeff))])
    end
end

checkFxnTime = toc(checkTime);
if printProgressUpdate
    disp(['Time for CHECK function = ',num2str(checkFxnTime)]);
    diary off %need this to write everything till this point to the diary
    diary on
end


probabilities = exp(-allFreeEnergies./(kB*T));
partitionFxn = sum(probabilities);
probabilities = probabilities./partitionFxn;

[sortedProbs,indexSortedProbs] =sort(probabilities,'descend'); %structure number indexSortedProbs(i) occurs with probability sortedProbs(i)
[sortedFE,indexSortedFE] = sort(allFreeEnergies,'ascend');
[minFE,whichStructMinFE] = min(allFreeEnergies);

if printProgressUpdate
    disp(['sortedProbs = ', num2str(sortedProbs(1:min(length(sortedProbs),5)))])
    disp(['indexSortedProbs = ', num2str(indexSortedProbs(1:min(length(indexSortedProbs),5)))])
    disp(['sortedFE = ', num2str(sortedFE(1:min(length(sortedFE),5)))])
    disp(['indexSortedFE = ', num2str(indexSortedFE(1:min(length(sortedFE),5)))]) %of course, indexSortedFE is the same as indexSortedProbs, but it doesn't hurt to check
    disp(['minFE = ',num2str(minFE)]);
    disp(['whichStructMinFE = ',num2str(whichStructMinFE)]);
    diary off %need this to write everything till this point to the diary
    diary on
end
    
if storeGraphs
    numGraphs = length(weightMatrixList);
    graphProbs = zeros(1,numGraphs);
    for i = 1:numStructures
        if allGraphs(indexSortedProbs(i))>0
            graphProbs(allGraphs(indexSortedProbs(i))) = graphProbs(allGraphs(indexSortedProbs(i))) + ...
                sortedProbs(i);
        end
    end

    [sortedGraphProbs,indexSortedGraphProbs] =sort(graphProbs,'descend'); %graph number indexSortedGraphProbs(i) occurs with probability sortedGraphProbs(i)
    if printProgressUpdate
        disp(['sortedGraphProbs = ', num2str(sortedGraphProbs(1:min(length(sortedGraphProbs),5)))])
        disp(['indexSortedProbs = ', num2str(indexSortedGraphProbs(1:min(length(indexSortedGraphProbs),5)))])
        diary off %need this to write everything till this point to the diary
        diary on
    end
else
    graphProbs = 0; sortedGraphProbs = 0; indexSortedGraphProbs = 0; numGraphs = 0;
end


%% Generating figures

if makingFigures
    figureFxnTime = RNALandscape_makeFiguresFxn(allFreeEnergies,minFE,...
        storeGraphs,sortedProbs,weightMatrixList,numBondsMatrixList,allGraphs,indexSortedProbs,...
        possiblePermutations,sequence,indexSortedGraphProbs,sortedGraphProbs,STable,numBP,whichStructMinFE);
end

totalTime = toc(startTime);
disp(['Total time= ',num2str(totalTime)]);

disp('____________________________________') %to separate entries in the diary
 
diary off %need this for the diary to display everything till this point