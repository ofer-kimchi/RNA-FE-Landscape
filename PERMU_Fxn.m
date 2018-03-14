function possiblePermutations = PERMU_Fxn(numRegions,C2,C3,C4,minNumCompatible,...
    numRegionsNecessaryforC3,numRegionsNecessaryforC4,printProgressUpdate)

%minNumCompatible = 0; %How many compatible regions make up a structure?
%set to 0 to include a structure that is completely unfolded as well
if numRegions > 247
    possiblePermutations = cell(1.7*10^9,1); %cell array of all possible structures.
%initially much larger than it needs to be, but unnecessary elements will be truncated
else
    possiblePermutations = cell(4.5*10^8,1);
end
numStructures = 1;
prevNumStructures = 0;

for i = 1:numRegions
    for j = i+1:numRegions
        if C2(i,j) ==1
            currPermu = zeros(1,numRegions);
            currPermu(1:2) = [i,j]; %the current set of regions we're considering as a structure
            lengthCurrPermu = 2;
            k = j; %the next region we'll try adding (we're about to add one so that's why it's j and not j+1).
            while lengthCurrPermu >= 2
                while k<numRegions
                    k = k+1;
                    mutuallyCompatible = true;
                    for l = currPermu(1:lengthCurrPermu)
                        if C2(k,l) == 0
                            mutuallyCompatible = false;
                            break
                        end
                    end
                    
                    %also check threeway compatibility
                    if numRegions>=numRegionsNecessaryforC3
                        if mutuallyCompatible %no reason to check this if we've already decided they aren't mutually compatible
                            for l = 1:lengthCurrPermu-1
                                for m = l+1:lengthCurrPermu
                                    if C3(currPermu(l),currPermu(m),k) == 0
                                        mutuallyCompatible = false;
                                        break
                                    end
                                end
                            end
                        end
                    end
                    
                    if numRegions>=numRegionsNecessaryforC4 %also check fourway compatibility
                        if mutuallyCompatible && lengthCurrPermu>2 %no reason to check this if we've already decided they aren't mutually compatible
                            if ~isempty(C4)
                                for l = 1:lengthCurrPermu-2
                                    for m = l+1:lengthCurrPermu-1
                                        for n = m+1:lengthCurrPermu
                                            if ~C4(currPermu(l),currPermu(m),currPermu(n),k)
                                                mutuallyCompatible = false;
                                                break
                                            end
                                        end
                                    end
                                end
                            else 
                                mutuallyCompatible = false;
                            end
                        end
                    end
                    
                    
                    if mutuallyCompatible
                        lengthCurrPermu = lengthCurrPermu+1;
                        currPermu(lengthCurrPermu) = k;
                    end
                end
                
                %do we have enough compatible regions to make a structure?
                enoughRegions = true; 
                if lengthCurrPermu < minNumCompatible
                    enoughRegions = false;
                end
                   
                if enoughRegions %if the permutation is long enough, add it to the list.
                    possiblePermutations{numStructures} = currPermu(1:lengthCurrPermu);
                    numStructures = numStructures+1;
                    if rem(numStructures,5*10^5)==0 && printProgressUpdate && prevNumStructures ~= numStructures
                        disp(['Setting up possiblePermutations: i = ',num2str(i)])
                        disp(['numStructures = ',num2str(numStructures)])
                        prevNumStructures = numStructures; %so we don't save the same thing multiple times
                        diary off %need this to write everything till this point to the diary
                        diary on
                    end
                end
                
                k = currPermu(lengthCurrPermu); %the next region we'll try adding
                currPermu(lengthCurrPermu) = 0; %are there alternatives to the last region we added?
                lengthCurrPermu = lengthCurrPermu-1;
            end
        end
    end
end


if minNumCompatible <= 1 %also add each region individually as a possible structure
    for i = 1:numRegions
        possiblePermutations{numStructures} = i;
        numStructures = numStructures+1;
    end
end

if minNumCompatible <= 0 %add another structure with no bonded bps
    possiblePermutations(numStructures+1:end) = []; %truncate unnecessary space we preallocated for possiblePermutations.
else
    possiblePermutations(numStructures:end) = [];
end

