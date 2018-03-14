function [numRegions, STable] = createSTable_Fxn(sequenceInNumbers,minBPInRegion,minBPInHairpin,allowParallel,substem)

%sequenceInNumbers is the sequence. RNA is represented by 1(A) 2(C) 3(G) or
%4(U); DNA by 5(A) 6(C) 7(G) or 8(T).

numBP = length(sequenceInNumbers);
if nargin == 3
    allowParallel = true; %should we allow pseudoknots to form?
end
allowOnlyCanonical = false; %should we allow GU bonds or other non-canonical base pairings? 
%(if we should allow GU pairs, then set parameter to false).


B = zeros(numBP,numBP); %matrix describing which bases can bond to which others
for i = 1:numBP
    for j = 1:numBP
        if sequenceInNumbers(i)==1 && sequenceInNumbers(j)==4 %for Watson-Crick base pairs, Bij=1
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 5 && sequenceInNumbers(j) == 8
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 1 && sequenceInNumbers(j) == 8
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 4 && sequenceInNumbers(j) == 5
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 2 && sequenceInNumbers(j) == 3
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 6 && sequenceInNumbers(j) == 7
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 2 && sequenceInNumbers(j) == 7
            B(i,j) = 1;
            B(j,i) = 1;
        elseif sequenceInNumbers(i) == 3 && sequenceInNumbers(j) == 6
            B(i,j) = 1;
            B(j,i) = 1;
        elseif ~allowOnlyCanonical && sequenceInNumbers(i) == 3 && sequenceInNumbers(j) == 4
            %can set Bij=2 if we want to distinguish non-Watson-Crick base
            %pairs, GU. This could allow us, for example, to ensure GU base
            %pairs don't close any region.
            B(i,j) = 1;
            B(j,i) = 1;
        end
    end
end


STable = cell(10^5,3); %creates empty cell array with preallocated space
%STable lists all possible regions, where a region is a consecutive set of
%bases that binds to another conseccutive set of bases.
%The first column of STable is the region number (a counter basically of
%the number of regions). The second column is a vector (e.g. [1 2 3 10 9
%8]) which is read as ntd 1 bonded to ntd 10; 2 to 9; 3 to 8. The third
%column is a matrix, (eg. [1 10; 2 9; 3 8]) with the same meaning.


%initialize some variables for use in constructing the STable
numRegions = 0; %should be zero, since we increase by one when we find another region.

%a region is started when we find two BPs, i and j, that can bind. Then we
%look if i+1 can bind to j-1, etc., and in the next for loop, if i+1 can
%bind to j+1, etc.

for i = 1:numBP-2*minBPInRegion+1 %because we need at least minBPInRegion consecutive base pairs binding
    for j = numBP:-1:i+2*minBPInRegion
        if B(i,j) == 1 && j-i-1>=minBPInHairpin
            currentRegionI = i;
            currentRegionJ = j;
            if minBPInRegion == 1 %Then even just ntd i bonded to j is a full region.
                numRegions = numRegions + 1;
                STable(numRegions,:) = {numRegions,[currentRegionI,currentRegionJ],[currentRegionI,currentRegionJ]};
            end
            
            firstTwoPairs = [[i,j];[i+1,j-1]];
            isSubset = false; %checks if the region we're starting to make is a subset of a previously-made region.
            for k = 1:numRegions
                if any(all((STable{k,3} == firstTwoPairs(1,:))')) && any(all((STable{k,3} == firstTwoPairs(2,:))'))
                        %this line (above) is equivalent to: all(ismember(firstTwoPairs,STable{k,3},'rows'))
                    isSubset = true;
                    break
                end
            end
            
            if ~isSubset
                numRegions = numRegions + 1;
                lenRegion = 0; %lenRegion is one less than the number of bps in a region.
                endOfRegion = false;
                while endOfRegion == false
                    lenRegion = lenRegion + 1;
                    if i+lenRegion>numBP || j-lenRegion<1 || j-lenRegion - (i+lenRegion) -1 < 1
                        endOfRegion = true;
                        lenRegion = lenRegion - 1; %to correct for adding one earlier
                    elseif B(i+lenRegion,j-lenRegion) ~=0
                        currentRegionI = [currentRegionI, i+lenRegion]; %#ok<*AGROW>
                        currentRegionJ = [currentRegionJ, j-lenRegion];
                    else
                        endOfRegion = true;
                        lenRegion = lenRegion - 1; %to correct for adding one earlier
                    end
                end

                %Pruning the regions we've added by mistake
                while j-lenRegion - (i+lenRegion) -1  < minBPInHairpin... %We cannot have hairpin loops of less than minBPInHairpin bases
                        || B(i+lenRegion,j-lenRegion)==2 %if we want to impose that GU cannot occur at the end of a region.

                    currentRegionI = currentRegionI(1:end-1);
                    currentRegionJ = currentRegionJ(1:end-1);
                    lenRegion = lenRegion - 1; %to correct for adding one by mistake before.
                    if lenRegion <1
                        break
                    end
                end

                if length(currentRegionI) < max(minBPInRegion,2) %because we need at least minBPInRegion consecutive bps in a region
                    %and we already included all regions of length 1.
                    numRegions = numRegions - 1; %to correct for adding one by mistake before.
                    %elseif  %it's important to have any other pruning constraints as elseif so that we don't oversubtract from numRegions
                    
                else
                    listOfPairs = [currentRegionI(1),currentRegionJ(1)]; %construct the third column of this entry of STable
                    for k = 2:length(currentRegionI)
                        listOfPairs = [listOfPairs;[currentRegionI(k),currentRegionJ(k)]];
                    end
                    STable(numRegions,:) = {numRegions,[currentRegionI,currentRegionJ],listOfPairs};
                end
            end
        end
    end
end

if allowParallel %if we allow pseudoknots, then we can have i and j both increasing through a region.
    %a region is started when we find two BPs, i and j, that can bind. Then we
    %look if i+1 can bind to j+1, etc.
    for i = 1:numBP-minBPInRegion+1 %because we need at least minBPInRegion consecutive base pairs binding
        %and we already considered regions consisting of single base pairs.
        for j = i+2:numBP-minBPInRegion+1
            %disp([i,j])
            if B(i,j) == 1 && j-i-1>=minBPInHairpin +minBPInRegion -1
                currentRegionI = i;
                currentRegionJ = j;
                
                firstTwoPairs = [[i,j];[i+1,j+1]];
                isSubset = false; %checks if the region we're starting to make is a subset of a previously-made region.
                for k = 1:numRegions
                    if any(all((STable{k,3} == firstTwoPairs(1,:))')) && any(all((STable{k,3} == firstTwoPairs(2,:))'))
                        %this line (above) is equivalent to: all(ismember(firstTwoPairs,STable{k,3},'rows'))
                        isSubset = true;
                        break
                    end
                end
                
                if ~isSubset
                    numRegions = numRegions + 1;
                    lenRegion = 0;
                    endOfRegion = false;
                    while endOfRegion == false
                        lenRegion = lenRegion + 1;
                        if i+lenRegion>numBP || j+lenRegion>numBP || j - (i+lenRegion) <= lenRegion+minBPInHairpin %1 
                            endOfRegion = true;
                            lenRegion = lenRegion - 1; %to correct for adding one earlier
                        elseif B(i+lenRegion,j+lenRegion) ~=0
                            currentRegionI = [currentRegionI, i+lenRegion];
                            currentRegionJ = [currentRegionJ, j+lenRegion];
                        else
                            endOfRegion = true;
                            lenRegion = lenRegion - 1; %to correct for adding one earlier
                        end
                    end
                    
                    %Pruning the regions we've added by mistake
                    while j - (i+lenRegion)  <= lenRegion + minBPInHairpin... %We cannot have hairpin loops of less than minBPInHairpin bases
                            || B(i+lenRegion,j+lenRegion)==2 %if we want to impose that GU cannot occur at the end of a region.
                        
                        currentRegionI = currentRegionI(1:end-1);
                        currentRegionJ = currentRegionJ(1:end-1);
                        lenRegion = lenRegion - 1; %to correct for adding one by mistake before.
                        if lenRegion <1
                            break
                        end
                    end
                    
                    if length(currentRegionI) < max(minBPInRegion,2) %because we need at least minBPInRegion consecutive bps in a region
                        %and we already included all regions of length 1.
                        numRegions = numRegions - 1; %to correct for adding one by mistake before.
                        %elseif  %it's important to have any other pruning constraints as elseif so that we don't oversubtract from numRegions
                        
                    else
                        listOfPairs = [currentRegionI(1),currentRegionJ(1)];
                        for k = 2:length(currentRegionI)
                            listOfPairs = [listOfPairs;[currentRegionI(k),currentRegionJ(k)]];
                        end
                        STable(numRegions,:) = {numRegions,[currentRegionI,currentRegionJ],listOfPairs};
                    end
                end
            end
        end
    end
end

%Based on a comment made by Zuker & Sankoff (1984) regarding the
%Pipas-McMahon algorithm, we also want to include subregions -- i.e.
%fractions of the whole region. Thus, we allow the possibility that two
%full regions may be incompatible, but that their compatible subregions are
%more desirable than a single full region.

if strcmp(substem,'all')
    for i = 1:numRegions
        regionI = STable{i,2}; %full region we're considering
        lengthRegion = length(regionI)/2;
        for j = 1:lengthRegion-minBPInRegion %we can make subregions of length lengthRegion-1 till minBPInRegion. There are j possible lengths.
            for k = 0:j %subregions come from getting rid of either edge.
                regionJ = [regionI(1+k:lengthRegion-j+k),regionI(lengthRegion+1+k:end-j+k)];
                %regionJ is the truncated region.
                %check if this region already exists. If not, add it.
                novelRegion = true;
                for l = 1:numRegions
                    if isequal(regionJ, STable{l,2})
                        novelRegion = false;
                        break
                    end
                end
                if novelRegion
                    numRegions = numRegions+1;

                    %add the region to STable
                    currentRegionI = regionJ(1:length(regionJ)/2);
                    currentRegionJ = regionJ(length(regionJ)/2+1:length(regionJ));
                    listOfPairs = [currentRegionI(1),currentRegionJ(1)];
                    for l = 2:length(currentRegionI)
                        listOfPairs = [listOfPairs;[currentRegionI(l),currentRegionJ(l)]];
                    end

                    STable(numRegions,:) = {numRegions,regionJ,listOfPairs}; %j tells us how long the subregion is compared to the full region.
                end
            end
        end
    end
else
    %first, remove all subregions that are already present (you can
    %probably be more efficient than this but the code is basically
    %instantaneous to run so it's not worth optimizing)
    for i = 1:numRegions
        regionI = STable{i,2}; %full region we're considering
        lengthRegion = length(regionI)/2;
        for j = 1:lengthRegion-minBPInRegion %we can make subregions of length lengthRegion-1 till minBPInRegion. There are j possible lengths.
            for k = 0:j %subregions come from getting rid of either edge.
                regionJ = [regionI(1+k:lengthRegion-j+k),regionI(lengthRegion+1+k:end-j+k)];
                %regionJ is the truncated region.
                %check if this region already exists. If not, add it.
                novelRegion = false;
                while ~novelRegion
                    novelRegion = true;
                    for l = 1:numRegions
                        if isequal(regionJ, STable{l,2})
                            novelRegion = false;
                            break
                        end
                    end
                    if ~novelRegion
                        STable(l,:) = [];
                        for q = l:size(STable,1)
                            STable{q,1} = STable{q,1} - 1;
                        end
                        numRegions = numRegions - 1;
                    end
                end 
            end
        end
    end

    
    %then, only add substems of length up to its maximum - "substem"
    for i = 1:numRegions
        regionI = STable{i,2}; %full region we're considering
        lengthRegion = length(regionI)/2;
        minBPInRegionSub = max(lengthRegion - substem,minBPInRegion);
        for j = 1:lengthRegion-minBPInRegionSub %we can make subregions of length lengthRegion-1 till minBPInRegion. There are j possible lengths.
            for k = 0:j %subregions come from getting rid of either edge.
                regionJ = [regionI(1+k:lengthRegion-j+k),regionI(lengthRegion+1+k:end-j+k)];
                %regionJ is the truncated region.
                %don't need to check if this region already exists. If not, add it.
                novelRegion = true;
                for l = 1:numRegions
                    if isequal(regionJ, STable{l,2})
                        novelRegion = false;
                        break
                    end
                end
                if novelRegion
                    numRegions = numRegions+1;

                    %add the region to STable
                    currentRegionI = regionJ(1:length(regionJ)/2);
                    currentRegionJ = regionJ(length(regionJ)/2+1:length(regionJ));
                    listOfPairs = [currentRegionI(1),currentRegionJ(1)];
                    for l = 2:length(currentRegionI)
                        listOfPairs = [listOfPairs;[currentRegionI(l),currentRegionJ(l)]];
                    end

                    STable(numRegions,:) = {numRegions,regionJ,listOfPairs}; %j tells us how long the subregion is compared to the full region.
                end
            end
        end
    end
end

STable(numRegions+1:end,:) = [];