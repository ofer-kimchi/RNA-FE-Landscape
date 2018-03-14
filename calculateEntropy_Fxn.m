function [loopEntropy,initialWeightMatrix,initialNumBondsMatrix] = ...
    calculateEntropy_Fxn(structure,strucMatrix,numBP,omega,minBPInHairpin,minBPInRegion,whichStruct,allPerms,vs,linkers,allowParallelStrands)
%Now we need to calculate the free energy contributions from the
%various loops in the structure. We start by enumerating the
%unbound bps, and use graphs to categorize the loops, and finally
%calculate the associated free energies.

%This function converts the structure being considered to a graph

%set up cell array for lists of unpaired regions
structUnpaired = cell(1,1); %list of unpaired regions
numStems = length(structure);

numUnpairedRegions = 0; %keeps track of number of unpaired regions
i = 1;
while i <= numBP
    %check if i is unpaired
    iUnpaired = 1-any(strucMatrix(i,:)==1); %0 if i is paired, 1 if it's unpaired
    if iUnpaired && ~any([linkers{:}]==i) %~ismember(i,[linkers{:}]) %don't put linker in this cell because this cell is used to construct singleEdges
        if numUnpairedRegions ==0 || ~any(structUnpaired{numUnpairedRegions}==i-1) %if i-1 is paired
            numUnpairedRegions = numUnpairedRegions+1;
            structUnpaired{numUnpairedRegions} = i;
        else
            structUnpaired{numUnpairedRegions} = [structUnpaired{numUnpairedRegions},i];
        end
    end
    i=i+1;
end

%set up graph of structure

%first set up nodes and double bonded edges
nodes = zeros(2,length(structure)); %each node is actually two bonded bps, so it is defined by two numbers (the positions of the two ntds in the sequence)
nodeCounter = 0;
stemCounter = 0;
numDoubleEdges = 0;

doubleEdges = zeros(2,numStems); %each doubleEdge denotes a double stranded region connecting two nodes.
%Each edge is defined by two numbers -- the numbers of the nodes it connects.
lengthDoubleEdges = zeros(1,numStems); %this array gives for each double the length of that edge.
%importantly, it is in the same ordering as doubleEdges so we can
%use the two arrays in conjunction.

for j = 0:length(linkers) %check if first ntd is part of any stem
    if j == 0
        linkerPos = 0;
    else
        linkerPos = linkers{j};
    end
    firstNtdInStem = false;
    for i = 1:numStems
        if any(structure{i}==linkerPos(end)+1) %i.e. the first ntd of any sequence
            firstNtdInStem = true;
        end
    end
    if ~firstNtdInStem %then the first ntd is not part of a stem, and we need to make a node for it (otherwise the next routine will do that automatically)
        nodeCounter = nodeCounter+1;
        nodes(:,nodeCounter) = [linkerPos(end),linkerPos(end)];
    end
end

numNodes = size(nodes,2);

pseudoSingleEdges = zeros(2,2*numNodes); %happens when two consecutive nucleotides are double bonded,
%but not part of the same strand of the same stem
numPseudoSingleEdges = 0;
while stemCounter < numStems
    nodeCounter = nodeCounter+1;
    stemCounter = stemCounter+1;
    stemLengthTot = length(structure{stemCounter});
    if stemLengthTot >2
        nodes(:,nodeCounter) = [structure{stemCounter}(1),structure{stemCounter}(stemLengthTot/2+1)];
        nodeCounter = nodeCounter+1;
        nodes(:,nodeCounter) = [structure{stemCounter}(stemLengthTot/2),structure{stemCounter}(stemLengthTot)];
        
        numDoubleEdges = numDoubleEdges + 1;
        doubleEdges(:,numDoubleEdges) = [nodeCounter-1, nodeCounter];
        lengthDoubleEdges(numDoubleEdges) = abs(structure{stemCounter}(stemLengthTot/2)-structure{stemCounter}(1));
        %if a double-stranded region contains 4 ntds, its length is actually 3.
        
        if structure{stemCounter}(stemLengthTot/2)+1== structure{stemCounter}(stemLengthTot/2+1)
            %then there is a pseudoSingleBond in addition to a double
            %bond between these nodes.
            numPseudoSingleEdges = numPseudoSingleEdges+1;
            pseudoSingleEdges(:,numPseudoSingleEdges) = [nodeCounter-1;nodeCounter];
        end
    elseif stemLengthTot == 2
        nodes(:,nodeCounter) = [structure{stemCounter}(1),structure{stemCounter}(stemLengthTot/2+1)];
    else
        disp('something strange happened')
    end
end
doubleEdges(:,numDoubleEdges+1:numStems) = [];

for j = 0:length(linkers) %check if last ntd is part of any stem
    if j == 0
        linkerPos = 0;
    else
        linkerPos = linkers{j};
    end
    lastNtdInStem = false;
    for i = 1:numStems
        if any(structure{i}==mod(linkerPos(1)-1,numBP+1)) %i.e. last element of any sequence
            lastNtdInStem = true;
        end
    end
    if ~lastNtdInStem %then the last ntd is not part of a stem, and we need to make a node for it
        nodeCounter = nodeCounter+1;
        nodes(:,nodeCounter) = [mod(linkerPos(1)-1,numBP+1)+1,mod(linkerPos(1)-1,numBP+1)+1]; %returns linkerPos(1) unless it's 0, in which case it returns numBP+1
    end
end


numNodes = size(nodes,2);

%now do pseudoSingleEdges
for i = 1:numNodes-1
    for j = i+1:numNodes
        numPseudoSingleEdgesToCreate = 0;
        if min(abs(nodes(:,i)-nodes(:,j)))==1  %meaning that the nodes are one ntd apart and thus
            %there are no unpaired ntds in between them, though
            %there should be an edge.
            numPseudoSingleEdgesToCreate = numPseudoSingleEdgesToCreate+1;
        end
        if min(abs(flipud(nodes(:,i))-nodes(:,j)))==1
            numPseudoSingleEdgesToCreate = numPseudoSingleEdgesToCreate+1;
        end
        for k =1:numPseudoSingleEdgesToCreate %because there may be more than one pseudoSingleEdge between the same two nodes
            %if the stem itself is only two ntds long, ignore it.
            if ~any(all((doubleEdges' == [i;j]')'))
                %this line (above) is equivalent to: ~ismember([i;j]',doubleEdges','rows')
                numPseudoSingleEdges = numPseudoSingleEdges + 1;
                pseudoSingleEdges(:,numPseudoSingleEdges) = [i;j];
            else
                if minBPInRegion >2
                    disp('tried to make an incorrect pseudoSingleEdge')
                end
            end
        end
    end
end
if minBPInHairpin == 0
    for i = 1:length(structure)
        if structure{i}(length(structure{i})/2) == structure{i}(end) -1
            if find(nodes(1,:)==structure{i}(length(structure{i})/2)) == find(nodes(2,:)== structure{i}(end))
                numPseudoSingleEdges = numPseudoSingleEdges + 1;
                psNode = find(nodes(1,:)==structure{i}(length(structure{i})/2));
                pseudoSingleEdges(:,numPseudoSingleEdges) = [psNode;psNode];
            end
        end
    end
end

pseudoSingleEdges(:,numPseudoSingleEdges+1:end) = [];


%find single-strand edges of graph
numSingleEdges = length(structUnpaired);
if isempty(structUnpaired{1}) %meaning actually there are no unpaired ntds
    numSingleEdges = 0;
end
singleEdges = zeros(2,numSingleEdges); %each singleEdge denotes a single stranded region connecting two nodes.
%Each edge is defined by two numbers -- the numbers of the nodes it connects.

singleEdgesSeq = zeros(2,numSingleEdges); %this array gives for each singleEdge the locations of the first and last ntds.
%importantly, it is in the same ordering as singleEdges so we can
%use the two arrays in conjunction.
for i = 1:numSingleEdges
    ssRegion = structUnpaired{i}; %vector of consecutive unbonded ntds.
    [~,firstNode] = find(nodes==ssRegion(1)-1,1); %only returns the first time it finds it,
    %so that if we're dealing with nodes at the very beginning
    %or end of the sequence, (which might just be the first or
    %last sequence element repeated twice) the node doesn't get
    %returned twice
    [~,secondNode] = find(nodes==ssRegion(end)+1,1);
    singleEdges(:,i) = [min(firstNode,secondNode),max(firstNode,secondNode)];
    singleEdgesSeq(:,i) = [min(ssRegion(1),ssRegion(end))-1,max(ssRegion(1),ssRegion(end))+1];
    for j = 0:length(linkers)
        if j == 0
            linkerPos = 0;
        else
            linkerPos = linkers{j};
        end
        if singleEdgesSeq(1,i) == linkerPos(end) %then it shouldn't actually be 0 it should be 1.
            singleEdgesSeq(1,i) = linkerPos(end)+1;
        end
        if singleEdgesSeq(2,i) == mod(linkerPos(1)-1,numBP+1)+1 %then it shouldn't actually be numBP+1 it should be numBP.
            singleEdgesSeq(2,i) = mod(linkerPos(1)-1,numBP+1)+1;
        end
    end
end



%create EdgeTable to define graph 
gEdgesMatrix = unique([doubleEdges,singleEdges,pseudoSingleEdges]','rows');
%list of all edges (no specifications yet as to weights or numBonds or lengths)
numEdges = length(gEdgesMatrix(:,1));
gEdgesMatrix = [gEdgesMatrix,zeros(numEdges,7)]; %columns of the matrix correspond to:
%(1&2) nodes between which there is an edge, whose properties are listed in the rest of the columns
%(3) weight of the edge -- sum of the number of bonds composing the edge, with double bonds counted as 2 and single bonds as 1.
%(4) number of bonds in the given edge
%(5) the length of a double bond between the nodes (if one exists) (and 0 otherwise)
%(6) the length of a single bond between the nodes (if one exists) (and 0 otherwise)
%(7) the length of a second single bond between the nodes (if one exists) (and 0 otherwise)
%(8) the length of a third single bond between the nodes (if one exists) (and 0 otherwise)
%(9) the length of a fourth single bond between the nodes (if one exists) (and 0 otherwise)

initialWeightMatrix = zeros(numNodes);
initialNumBondsMatrix = zeros(numNodes);

for i = 1:numSingleEdges
    %this line (below) is equivalent to: [~, index] = ismember(singleEdges(:,i)',gEdgesMatrix(:,1:2),'rows'); 
    index = find(all((gEdgesMatrix(:,1:2) == singleEdges(:,i)')')); %EndNodes 
    gEdgesMatrix(index,3) = gEdgesMatrix(index,3)+1; %update weight by adding 1
    gEdgesMatrix(index,4) = gEdgesMatrix(index,4)+1; %update numBonds by adding 1
    initialWeightMatrix(singleEdges(1,i),singleEdges(2,i)) = initialWeightMatrix(singleEdges(1,i),singleEdges(2,i)) + 1;
    initialNumBondsMatrix(singleEdges(1,i),singleEdges(2,i)) = initialNumBondsMatrix(singleEdges(1,i),singleEdges(2,i)) + 1;
    if singleEdges(1,i) ~= singleEdges(2,i)
        initialWeightMatrix(singleEdges(2,i),singleEdges(1,i)) = initialWeightMatrix(singleEdges(2,i),singleEdges(1,i)) + 1;
        initialNumBondsMatrix(singleEdges(2,i),singleEdges(1,i)) = initialNumBondsMatrix(singleEdges(2,i),singleEdges(1,i)) + 1;
    end
    
    if gEdgesMatrix(index,6) == 0 %singleBond1Length
        gEdgesMatrix(index,6) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i); %if a single stranded region contains 4 ntds, its length is actually 5
    elseif gEdgesMatrix(index,7) == 0 %singleBond2Length
        gEdgesMatrix(index,7) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
    elseif gEdgesMatrix(index,8) == 0 %singleBond3Length
        gEdgesMatrix(index,8) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
    elseif gEdgesMatrix(index,9) == 0 %singleBond4Length
        gEdgesMatrix(index,9) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
    else
        disp('something is off')
    end
end

for i = 1:numDoubleEdges
    %this line (below) is equivalent to: [~, index] = ismember(doubleEdges(:,i)',gEdgesMatrix(:,1:2),'rows');
    index = find(all((gEdgesMatrix(:,1:2) == doubleEdges(:,i)')')); %EndNodes 
    gEdgesMatrix(index,3) = gEdgesMatrix(index,3)+2; %update weight by adding 2
    gEdgesMatrix(index,4) = gEdgesMatrix(index,4)+1; %update numBonds by adding 1
    initialWeightMatrix(doubleEdges(1,i),doubleEdges(2,i)) = initialWeightMatrix(doubleEdges(1,i),doubleEdges(2,i)) + 2;
    initialNumBondsMatrix(doubleEdges(1,i),doubleEdges(2,i)) = initialNumBondsMatrix(doubleEdges(1,i),doubleEdges(2,i)) + 1;
    if doubleEdges(1,i) ~= doubleEdges(2,i)
        initialWeightMatrix(doubleEdges(2,i),doubleEdges(1,i)) = initialWeightMatrix(doubleEdges(2,i),doubleEdges(1,i)) + 2;
        initialNumBondsMatrix(doubleEdges(2,i),doubleEdges(1,i)) = initialNumBondsMatrix(doubleEdges(2,i),doubleEdges(1,i)) + 1;
    end
    
    if gEdgesMatrix(index,5) == 0 %doubleBondLength
        gEdgesMatrix(index,5) = lengthDoubleEdges(i);
    else
        disp('something is off 2')
    end
end

for i = 1:numPseudoSingleEdges
    %this line (below) is equivalent to: [~, index] = ismember(pseudoSingleEdges(:,i)',gEdgesMatrix(:,1:2),'rows');
    index = find(all((gEdgesMatrix(:,1:2) == pseudoSingleEdges(:,i)')')); %EndNodes %we know that this edge is a member, just want to find its index.
    gEdgesMatrix(index,3) = gEdgesMatrix(index,3)+1; %update weight by adding 1
    gEdgesMatrix(index,4) = gEdgesMatrix(index,4)+1; %update numBonds by adding 1
    initialWeightMatrix(pseudoSingleEdges(1,i),pseudoSingleEdges(2,i)) = initialWeightMatrix(pseudoSingleEdges(1,i),pseudoSingleEdges(2,i)) + 1;
    initialNumBondsMatrix(pseudoSingleEdges(1,i),pseudoSingleEdges(2,i)) = initialNumBondsMatrix(pseudoSingleEdges(1,i),pseudoSingleEdges(2,i)) + 1;
    if pseudoSingleEdges(1,i)~=pseudoSingleEdges(2,i)
        initialWeightMatrix(pseudoSingleEdges(2,i),pseudoSingleEdges(1,i)) = initialWeightMatrix(pseudoSingleEdges(2,i),pseudoSingleEdges(1,i)) + 1;
        initialNumBondsMatrix(pseudoSingleEdges(2,i),pseudoSingleEdges(1,i)) = initialNumBondsMatrix(pseudoSingleEdges(2,i),pseudoSingleEdges(1,i)) + 1;
    end
    
    if gEdgesMatrix(index,6) == 0 %singleBond1Length
        gEdgesMatrix(index,6) = 1;
    elseif gEdgesMatrix(index,7) == 0 %singleBond2Length
        gEdgesMatrix(index,7) = 1;
    elseif gEdgesMatrix(index,8) == 0 %singleBond3Length
        gEdgesMatrix(index,8) = 1;
    elseif gEdgesMatrix(index,9) == 0 %singleBond4Length
        gEdgesMatrix(index,9) = 1;
    else
        disp('something is off 3')
    end
end


%create graph
weightMatrix = initialWeightMatrix;
numBondsMatrix = initialNumBondsMatrix;

%%
[loopEntropy,componentGraphs] = calculateEntropyFromGraph_Fxn(gEdgesMatrix,omega,minBPInHairpin,weightMatrix,numBondsMatrix,allPerms,vs);

%% The code from here on can be commented out with no effect
if nargin <= 10 
    allowParallelStrands = true;
end

%To be doubly sure we got rid of all parallel strands, penalize any
%parallel strands that seemed to form
if ~allowParallelStrands && ~isempty(linkers) %you can have e.g. open nets 1 not be parallel if it's a duplex.
    %penalize parallel strands
    parallelStrandsEntropyPenalty = 10000; %entropy penalty for parallel strands to get rid of parallel strands
    if componentGraphs(4)>0 || componentGraphs(5)>0 || componentGraphs(10)>0 || ...
            componentGraphs(11)>0 || componentGraphs(12)>0 ||componentGraphs(13)>0 ||componentGraphs(14)>0 || componentGraphs(15)>0
        loopEntropy = loopEntropy - parallelStrandsEntropyPenalty; %doesn't include open & closed nets 2-single
    end
end
