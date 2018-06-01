function [perBaseGraph] = structureString2perBaseGraph(stringStructure) 
%%
%takes as input a string '10 9 8 0 0 0 3 2 1' which means ntd 1 is bonded
%to 10, 2 to 9, 3 to 8, 4-7 are unbonded.
%Returns a cell vector perBaseGraph of length numBP-1 where perBaseGraph(i) is
%a string representing what topology the bp i is a part of (e.g. stem of
%hairpin loop, etc.)

minBPInRegion = 1;
componentGraphsVec = zeros(1,15);

strucBPs = structreString2structBPs(stringStructure);
structure = structBPs2structure(strucBPs);

numBP = length(stringStructure(stringStructure == ' '));
if stringStructure(end) ~= ' '
    numBP = numBP + 1;
end

perBaseGraph = cell(1,numBP-1);

strucMatrix = zeros(numBP);
for i = 1:size(strucBPs,1)
    strucMatrix(strucBPs(i,1),strucBPs(i,2)) = 1;
    strucMatrix(strucBPs(i,2),strucBPs(i,1)) = 1;
end

allPerms = cell(1,4);
for i = 1:4
    allPerms{i} = perms(1:i);
end

minBPInHairpin = 3;

%%
%taken from calculateEntropy_10_25_17_noGraph
structUnpaired = cell(1,1); %list of unpaired regions
numStems = length(structure);

numUnpairedRegions = 0; %keeps track of number of unpaired regions
i = 1;
while i <= numBP
    %check if i is unpaired
    iUnpaired = 1-any(strucMatrix(i,:)==1); %0 if i is paired, 1 if it's unpaired
    if iUnpaired
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
doubleEdgesSeq = zeros(2,0); %this array gives for each singleEdge the locations of the first and last ntds.

%check if first ntd is part of any stem
firstNtdInStem = false;
for i = 1:numStems
    if any(structure{i}==1)
        firstNtdInStem = true;
    end
end
if ~firstNtdInStem %then the first ntd is not part of a stem, and we need to make a node for it (otherwise the next routine will do that automatically)
    nodeCounter = nodeCounter+1;
    nodes(:,nodeCounter) = [0,0];
end
numNodes = size(nodes,2);

pseudoSingleEdges = zeros(2,2*numNodes); %happens when two consecutive nucleotides are double bonded,
%but not part of the same strand of the same stem
numPseudoSingleEdges = 0;

pseudoSingleEdgesSeq = zeros(2,numPseudoSingleEdges);
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

        doubleEdgesSeq(:,size(doubleEdgesSeq,2)+1) = sort([structure{stemCounter}(stemLengthTot/2),structure{stemCounter}(1)]);
        %doubleEdgesSeq(:,size(doubleEdgesSeq,2)+1) = sort([structure{stemCounter}(stemLengthTot/2+1),structure{stemCounter}(end)]);

        if structure{stemCounter}(stemLengthTot/2)+1== structure{stemCounter}(stemLengthTot/2+1)
            %then there is a pseudoSingleBond in addition to a double
            %bond between these nodes.
            numPseudoSingleEdges = numPseudoSingleEdges+1;
            pseudoSingleEdges(:,numPseudoSingleEdges) = [nodeCounter-1;nodeCounter];
            pseudoSingleEdgesSeq(:,numPseudoSingleEdges) = sort([structure{stemCounter}(stemLengthTot/2),structure{stemCounter}(1)]);
        end
    elseif stemLengthTot == 2
        nodes(:,nodeCounter) = [structure{stemCounter}(1),structure{stemCounter}(stemLengthTot/2+1)];
    else
        disp('something strange happened')
    end
end
doubleEdges(:,numDoubleEdges+1:numStems) = [];
%check if last ntd is part of any stem
lastNtdInStem = false;
for i = 1:numStems
    if any(structure{i}==numBP)
        lastNtdInStem = true;
    end
end
if ~lastNtdInStem %then the last ntd is not part of a stem, and we need to make a node for it
    nodeCounter = nodeCounter+1;
    nodes(:,nodeCounter) = [numBP+1,numBP+1];
end

%     [~,nodeWithLastNtd]=find(nodes==numBP+1);
%     if isempty(nodeWithLastNtd)
%         [~,nodeWithLastNtd]=find(nodes==numBP);
%     end
%     nodeWithLastNtd = nodeWithLastNtd(1);
%     [~,nodeWithFirstNtd]=find(nodes==1);
%     if isempty(nodeWithFirstNtd)
%         [~,nodeWithFirstNtd]=find(nodes==0);
%     end
%     nodeWithFirstNtd = nodeWithFirstNtd(1);

numNodes = size(nodes,2);

%now do pseudoSingleEdges
for i = 1:numNodes-1
    for j = i+1:numNodes
        numPseudoSingleEdgesToCreate = 0;
        pseudoSingleEdgesSeqToCreate = [[0 0];[0 0]];
        if min(abs(nodes(:,i)-nodes(:,j)))==1  %meaning that the nodes are one ntd apart and thus
            %there are no unpaired ntds in between them, though
            %there should be an edge.
            numPseudoSingleEdgesToCreate = numPseudoSingleEdgesToCreate+1;
            if ~any(all((doubleEdges' == [i;j]')'))
                nodesI = nodes(:,i); nodesI = nodesI(abs(nodes(:,i)-nodes(:,j))==1);
                nodesJ = nodes(:,j); nodesJ = nodesJ(abs(nodes(:,i)-nodes(:,j))==1);
                pseudoSingleEdgesSeqToCreate(numPseudoSingleEdgesToCreate,:) = sort([nodesI nodesJ]);
            end
        end
        if min(abs(flipud(nodes(:,i))-nodes(:,j)))==1
            numPseudoSingleEdgesToCreate = numPseudoSingleEdgesToCreate+1;
            if ~any(all((doubleEdges' == [i;j]')'))
                nodesI = flipud(nodes(:,i)); nodesI = nodesI(abs(flipud(nodes(:,i))-nodes(:,j))==1);
                nodesJ = nodes(:,j); nodesJ = nodesJ(abs(flipud(nodes(:,i))-nodes(:,j))==1);
                pseudoSingleEdgesSeqToCreate(numPseudoSingleEdgesToCreate,:) = sort([nodesI nodesJ]);
            end
        end
        for k =1:numPseudoSingleEdgesToCreate %because there may be more than one pseudoSingleEdge between the same two nodes
            %if the stem itself is only two ntds long, ignore it.
            %I made a substitution for ismember to make the code run
            %faster. Check that this substitution is correct:
%             if  ~any(all((doubleEdges' == [i;j]')')) ~= ~ismember([i;j]',doubleEdges','rows')
%                 disp('I made a mistake in substituting ismember')
%             end
            if ~any(all((doubleEdges' == [i;j]')'))
                %this line (above) is equivalent to: ~ismember([i;j]',doubleEdges','rows')
                numPseudoSingleEdges = numPseudoSingleEdges + 1;
                pseudoSingleEdges(:,numPseudoSingleEdges) = [i;j];
                pseudoSingleEdgesSeq(:,numPseudoSingleEdges) = pseudoSingleEdgesSeqToCreate(k,:);
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
    if singleEdgesSeq(1,i) == 0 %then it shouldn't actually be 0 it should be 1.
        singleEdgesSeq(1,i) = 1;
    end
    if singleEdgesSeq(2,i) == numBP+1 %then it shouldn't actually be numBP+1 it should be numBP.
        singleEdgesSeq(2,i) = numBP;
    end
end


%%
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

gEdgesMatrixWithNtds  = cell(size(gEdgesMatrix)); %for each edge in gEdgesMatrix, instead of writing the length, 
    %write which ntds contribute to the double bond, single bond, etc.


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
        gEdgesMatrixWithNtds{index,6} = singleEdgesSeq(1,i):singleEdgesSeq(2,i);
    elseif gEdgesMatrix(index,7) == 0 %singleBond2Length
        gEdgesMatrix(index,7) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
        gEdgesMatrixWithNtds{index,7} = singleEdgesSeq(1,i):singleEdgesSeq(2,i);
    elseif gEdgesMatrix(index,8) == 0 %singleBond3Length
        gEdgesMatrix(index,8) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
        gEdgesMatrixWithNtds{index,8} = singleEdgesSeq(1,i):singleEdgesSeq(2,i);
    elseif gEdgesMatrix(index,9) == 0 %singleBond4Length
        gEdgesMatrix(index,9) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
        gEdgesMatrixWithNtds{index,9} = singleEdgesSeq(1,i):singleEdgesSeq(2,i);
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
        gEdgesMatrixWithNtds{index,5} = doubleEdgesSeq(1,i):doubleEdgesSeq(2,i);
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
        gEdgesMatrix(index,6) = 1; %NOT SURE IF THIS IS ACCURATE.
        gEdgesMatrixWithNtds{index,6} = pseudoSingleEdgesSeq(1,i):pseudoSingleEdgesSeq(2,i);
    elseif gEdgesMatrix(index,7) == 0 %singleBond2Length
        gEdgesMatrix(index,7) = 1;
        gEdgesMatrixWithNtds{index,7} = pseudoSingleEdgesSeq(1,i):pseudoSingleEdgesSeq(2,i);
    elseif gEdgesMatrix(index,8) == 0 %singleBond3Length
        gEdgesMatrix(index,8) = 1;
        gEdgesMatrixWithNtds{index,8} = pseudoSingleEdgesSeq(1,i):pseudoSingleEdgesSeq(2,i);
    elseif gEdgesMatrix(index,9) == 0 %singleBond4Length
        gEdgesMatrix(index,9) = 1;
        gEdgesMatrixWithNtds{index,9} = pseudoSingleEdgesSeq(1,i):pseudoSingleEdgesSeq(2,i);
    else
        disp('something is off 3')
    end
end


%create graph
weightMatrix = initialWeightMatrix;
numBondsMatrix = initialNumBondsMatrix;



%%
%from calculateEntropyFromGraph_10_25_17_noGraph

b = 0.8/(0.33); %in units of nucleotides (approximate size of ntd is 3.3 angstroms)
    %This assumes a persistence length of RNA of ~0.8 nm found in The Snakelike Chain Character of Unstructured RNA
    %by  David R. Jacobson, Dustin B. McIntosh, and Omar A. Saleh
    %As a plus, this agrees well with Siggia's formula b=2.5a (where a is the size of one ntd).
%vs = 0.01378; %in units of ntds^3. Calculated using the avg bond length (rho) in Tristan's simulations
    %which was 0.504385 nm. More precisely, vs = 14.9566 ntds^3.
    %0.01378 is used to fit deltaS to experimental results in e.g. Fig. 6 of
    %Discrete state model and accurate estimation of loop entropy of RNA secondary structures
    %by Jian Zhang, Ming Lin, Rong Chen, Wei Wang, and Jie Liang
beta = 3/(2*b);

kB = 0.0019872; %units of kcal/(mol*Kelvin)

loopEntropy = 0;
omega = ones(numBP,1);
vs = 0.02;

%define graphs for the nets we know of

openNet0WM = [[0,1];[1,0]]; %weight matrix for openNet0. This is the simplest version of open-net-0. Just two connected nodes.
openNet0NBM = [[0,1];[1,0]]; %numBonds matrix for openNet0

closedNet0WM = 1;
closedNet0NBM = 1;

closedNets0_2WM = 2;
closedNets0_2NBM = 2;

openNet1WM = [[0,3];[3,0]];
openNet1NBM = [[0,2];[2,0]];
 %like closed net 0 and open net 0, an openNet1 has two nodes connected by an edge
    %To distinguish openNet1 from closedNet0 and openNet0, we have the
    %edge have weight 3 with numBonds =2 because there is one double bond
    %and one single bond.

closedNet1WM = [[0,4];[4,0]];
closedNet1NBM = [[0,3];[3,0]];
%closedNet1 has two single bonds and a double bond. 
    %Hence, the edge between the two nodes has weight 4 and numBonds =3.

openNet2aWM = [[0 2 1 0];[2 0 1 1];[1 1 0 2];[0 1 2 0]];
openNet2aNBM = [[0 1 1 0];[1 0 1 1];[1 1 0 1];[0 1 1 0]];

closedNet2aWM = [[0 2 1 1];[2 0 1 1];[1 1 0 2];[1 1 2 0]];
closedNet2aNBM = [[0 1 1 1];[1 0 1 1];[1 1 0 1];[1 1 1 0]];

openNet2bWM = [[0 3 1 0];[3 0 0 1];[1 0 0 2];[0 1 2 0]];
openNet2bNBM = [[0 2 1 0];[2 0 0 1];[1 0 0 1];[0 1 1 0]];

closedNet2bWM = [[0 3 1 0];[3 0 0 1];[1 0 0 3];[0 1 3 0]];
closedNet2bNBM = [[0 2 1 0];[2 0 0 1];[1 0 0 2];[0 1 2 0]];

openNet2cWM = [[0 2 2 0];[2 0 0 1];[2 0 0 2];[0 1 2 0]];
openNet2cNBM = [[0 1 2 0];[1 0 0 1];[2 0 0 1];[0 1 1 0]];

closedNet2cWM = [[0 2 2 0];[2 0 0 2];[2 0 0 2];[0 2 2 0]];
closedNet2cNBM = [[0 1 2 0];[1 0 0 2];[2 0 0 1];[0 2 1 0]];

openNet2SingleWM = [[0,3];[3,0]];
openNet2SingleNBM = [[0,3];[3,0]];
%like openNet1 but the double bond is replaced by two single bonds
%more relevantly, like an openNet2a or openNet2c (the expressions yield the same result) with l1 = l2 = 0.

closedNet2SingleWM = [[0,4];[4,0]];
closedNet2SingleNBM = [[0,4];[4,0]];
%like closedNet1 but the double bond is replaced by two single bonds
%more relevantly, like a closedNet2a or closedNed2c (the expressions yield the same result) with l1 = l2 = 0.

openNet2cL1ZeroWM = [[0 2 1];[2 0 2];[1 2 0]];
openNet2cL1ZeroNBM = [[0 2 1];[2 0 1];[1 1 0]];

closedNet2cL1ZeroWM = [[0 2 2];[2 0 2];[2 2 0]];
closedNet2cL1ZeroNBM = [[0 2 2];[2 0 1];[2 1 0]];
%like openNet2c or openNet2a with l1 = 0 (and l2 nonzero).


%if you take open(closed)Net2b and set l1 = l2 = 0, you end up with 2(3)
%closedNets0. The two expressions match up, so all is good.

%a check:
%openNet2bIsomorphism = graph([[0 1 0 2];[1 0 3 0];[0 3 0 1];[2 0 1 0]]);
%openNet2bIsomorphism.Edges.numBonds = [1 1 2 1]'; 

%------------------------------------------------------------------------
%%

numEdges = size(gEdgesMatrix,1);
listEdges = 1:numEdges;
useEffectiveHelixLengths = false; %should we correct helix lengths using the equation from Xayaphoummine et al.?
loopEntropy = 0;

numDisconnectedGraphs = max(myConnComp2(weightMatrix)); %number of disconnected graphs g can be separated into with the deletion of double stranded regions.
%will be equal to 1 since this is the start of the code.

%try to remove each bond , and if you can -- by
%removing it -- separate the graph into two disconnected graphs,
%then remove the bond. By doing this, we're trying to
%decompose the graph into the smallest nets possible.

oneBondList = listEdges(gEdgesMatrix(:,4) == 1); %list of edges corresponding to numBonds = 1
%don't want to include bonds which connect a node to itself since we deal
%with those separately.
for i = length(oneBondList):-1:1 %go backwards so you can delete edges from gEdges without affecting oneBondList
    edgeIndex = oneBondList(i); 
    hWeightMatrix = weightMatrix; 
    hWeightMatrix(gEdgesMatrix(edgeIndex,1),gEdgesMatrix(edgeIndex,2)) = 0;
    hWeightMatrix(gEdgesMatrix(edgeIndex,2),gEdgesMatrix(edgeIndex,1)) = 0;
    %hWeightMatrix is the weightMatrix when you remove bond oneBondList(i) from weightMatrix.
    potentialNumDisconnectedGraphs = max(myConnComp2(hWeightMatrix));
    if potentialNumDisconnectedGraphs>numDisconnectedGraphs %if h can be separated into more disconnected graphs than g can
        if gEdgesMatrix(edgeIndex,6)>=1 %singleBond1Length
            s0 = gEdgesMatrix(edgeIndex,6);
            loopEntropy = loopEntropy + kB*log(omega(s0));
            componentGraphsVec(3) = componentGraphsVec(3) + 1;
            for j = 1:length(gEdgesMatrixWithNtds{edgeIndex,6})-1
                perBaseGraph{gEdgesMatrixWithNtds{edgeIndex,6}(j)} = 'open_net_0';
            end
        elseif gEdgesMatrix(edgeIndex,5)>=1 %then the bond deleted was a double bond
            for j = 1:length(gEdgesMatrixWithNtds{edgeIndex,5})-1
                perBaseGraph{gEdgesMatrixWithNtds{edgeIndex,5}(j)} = 'double bond'; %meaning just a double bond that's not part of a larger structure
            end
        else
            for j = 1:100
                disp('?????????')
            end
        end
        
        weightMatrix = hWeightMatrix; %remove bond 
        numBondsMatrix(gEdgesMatrix(edgeIndex,1),gEdgesMatrix(edgeIndex,2)) = 0;
        numBondsMatrix(gEdgesMatrix(edgeIndex,2),gEdgesMatrix(edgeIndex,1)) = 0;
        gEdgesMatrix(edgeIndex,:) = []; %update the edge table as well
        gEdgesMatrixWithNtds(edgeIndex,:) = [];
        numDisconnectedGraphs = potentialNumDisconnectedGraphs;
    end
end


%If any ntd was only connected to the rest of a graph by open-nets-0,
%then the node is unimportant and can (and should) be
%deleted. We therefore remove all nodes with no edges.
nodesWithNoEdges = find(sum(weightMatrix)==0);
% if isempty(weightMatrix) %then find(sum(weightMatrix)) will return 0 instead of empty array
%     if nargout > 1
%         componentGraphs = componentGraphsVec;
%     end
%     return
% end

%Update the names of the nodes in gEdgesMatrix and update weightMatrix and numBondsMatrix:
for noEdgesCounter = length(nodesWithNoEdges):-1:1 %important for it to descend here
    gEdgesMatrix(gEdgesMatrix(:,1)>nodesWithNoEdges(noEdgesCounter),1) = ...
        gEdgesMatrix(gEdgesMatrix(:,1)>nodesWithNoEdges(noEdgesCounter),1) - 1;
    gEdgesMatrix(gEdgesMatrix(:,2)>nodesWithNoEdges(noEdgesCounter),2) = ...
        gEdgesMatrix(gEdgesMatrix(:,2)>nodesWithNoEdges(noEdgesCounter),2) - 1;
    weightMatrix(:,nodesWithNoEdges(noEdgesCounter)) = [];
    weightMatrix(nodesWithNoEdges(noEdgesCounter),:) = [];
    numBondsMatrix(:,nodesWithNoEdges(noEdgesCounter)) = [];
    numBondsMatrix(nodesWithNoEdges(noEdgesCounter),:) = [];
end


%for testing:
%         threeClosedNet0 = graph([[0 1 1];[1 0 1];[1 1 0]])
%         threeClosedNet0.Edges.Weight = [];
%         threeClosedNet0.Edges.Weight = [1 1 1]';
%         threeClosedNet0.Edges.numBonds = [1 1 1]';
%         threeClosedNet0.Edges.singleBond1Length = [3 4 5]'
%         threeClosedNet0.Edges.singleBond2Length = [0 0 0]';
%         threeClosedNet0.Edges.doubleBondLength = [0 0 0]';
%
%         fourClosedNet0 = graph([[0 1 0 1];[1 0 1 0];[0 1 0 1];[1 0 1 0]])
%         fourClosedNet0.Edges.Weight = [];
%         fourClosedNet0.Edges.Weight = [1 1 1 1]';
%         fourClosedNet0.Edges.numBonds = [1 1 1 1]';
%         fourClosedNet0.Edges.singleBond1Length = [3 4 5 6]'
%         fourClosedNet0.Edges.singleBond2Length = [0 0 0 0]';
%         fourClosedNet0.Edges.doubleBondLength = [0 0 0 0]';

%%
isUnchanged = false; %make sure to do the following loop at least once
while ~isUnchanged
    %at the end of this loop, we want to check if g has changed -- if it has, we might need to
    %rerun this loop
    hWeightMatrix = weightMatrix;
    
    %if a node is connected to itself, we know that it's a closed-net-0.
    %Therefore, we make a new node that is connected to itself with the
    %same length, and get rid of the closed-net-0 edge (this deals with
    %errors of the form
    
    %     EndNodes    Weight    numBonds    doubleBondLength    singleBond1Length    singleBond2Length    singleBond3Length    singleBond4Length
    %     ________    ______    ________    ________________    _________________    _________________    _________________    _________________
    %
    %     1    1      1         1           0                   3                    0                    0                    0
    %     1    2      1         1           0                   1                    0                    0                    0
    %     1    3      1         1           0                   2                    0                    0                    0
    %     2    3      2         1           1                   0                    0                    0                    0
    %     2    5      1         1           0                   6                    0                    0                    0
    %     3    4      1         1           0                   1                    0                    0                    0
    %     4    5      3         2           1                   3                    0                    0                    0
    
    index = find(gEdgesMatrix(:,1) == gEdgesMatrix(:,2)); %index of edge connecting a node to itself
    i = 1;
    while i <= length(index)
        if gEdgesMatrix(index(i),7) ~= 0 %singleBond2Length %don't want to get rid of double hairpin loops
            index(i) = [];
        else 
            i = i+1;
        end
    end
    counter = 1;
    while counter <= length(index)
        if nnz(weightMatrix(gEdgesMatrix(index(counter),1),:))==1 
            index(counter) = []; %also don't need to mess with nodes that aren't connected to anything but themselves
        else
            counter = counter + 1;
        end
    end
    %         if any(gEdges.doubleBondLength(index) ~= zeros(length(index),1)) || ...
    %                 any(gEdges.singleBond3Length(index) ~= zeros(length(index),1)) || ...
    %                 any(gEdges.singleBond4Length(index) ~= zeros(length(index),1))
    %             disp('not sure why a node is connected to itself in this way')
    %         end
    
    indexCounter = 0;
    lengthIndexInitial = length(index);
    while ~isempty(index)
        indexCounter = indexCounter + 1;
        if indexCounter > lengthIndexInitial
            for i = 1:100 %make sure I really see it
                disp(strcat('problem with whichStruct = ',num2str(whichStruct)));
            end
            if indexCounter > length(index)*100
                if nargout > 1
                    componentGraphs = componentGraphsVec;
                end
                return %so you don't get stuck in this loop and so the message will be printed
            end
        end
        
        %make a new node that's connected to itself with this single bond
        edgeToRemove = index(1);
        numNodesG = length(weightMatrix)+1;
        edgeToAdd = [numNodesG, numNodesG,1,1,0,gEdgesMatrix(edgeToRemove,6),0,0,0];
        edgeToAddCell = {[],[],[],[],[],gEdgesMatrixWithNtds{edgeToRemove,6},[],[],[]};
        %add the new edge to gEdges
        [gEdgesMatrix,sortOrder] = sortrows([gEdgesMatrix;edgeToAdd]); %sortrows doesn't do anything
        gEdgesMatrixWithNtds = [gEdgesMatrixWithNtds; edgeToAddCell]; gEdgesMatrixWithNtds = gEdgesMatrixWithNtds(sortOrder,:); %#ok<AGROW>
        
        
        weightMatrix(numNodesG, numNodesG) = 1; %add the new node and edge to weightMatrix and numBondsMatrix
        numBondsMatrix(numNodesG, numNodesG) = 1;
        
        %now get rid of the old edge connecting the node to itself
        weightMatrix(gEdgesMatrix(edgeToRemove,1),gEdgesMatrix(edgeToRemove,1)) = 0; 
            %this syntax works because the node is connected to itself
        numBondsMatrix(gEdgesMatrix(edgeToRemove,1),gEdgesMatrix(edgeToRemove,1)) = 0; 
        gEdgesMatrix(edgeToRemove,:) = [];  %make sure to do this after deleting the node from weightMatrix and numBondsMatrix!
        for j = 1:length(gEdgesMatrixWithNtds{edgeToRemove,6})-1
            perBaseGraph{gEdgesMatrixWithNtds{edgeToRemove,6}(j)} = 'closed_net_0';
        end
        gEdgesMatrixWithNtds(edgeToRemove,:) = [];
        
        %repeat the process
        index = find(gEdgesMatrix(:,1) == gEdgesMatrix(:,2)); %index of edge connecting a node to itself
        i = 1;
        while i <= length(index)
            if gEdgesMatrix(index(i),7) ~= 0 %singleBond2Length %don't want to get rid of double hairpin loops
                index(i) = [];
            else
                i = i+1;
            end
        end
        counter = 1;
        while counter <= length(index)
            if nnz(weightMatrix(gEdgesMatrix(index(counter),1),:))==1
                index(counter) = []; %also don't need to mess with nodes that aren't connected to anything but themselves
            else
                counter = counter + 1;
            end
        end
    end
    
    
    %If any ntd was only connected to the rest of a graph by open-nets-0,
    %then the node is unimportant and can (and should) be
    %deleted. We therefore remove all nodes with no edges.
    nodesWithNoEdges = find(sum(weightMatrix)==0);
    if isempty(weightMatrix) %then find(sum(weightMatrix)) will return 0 instead of empty array
        if nargout > 1
            componentGraphs = componentGraphsVec;
        end
        return
    end    
    
    %Update the names of the nodes in gEdgesMatrix and update weightMatrix and numBondsMatrix:
    for noEdgesCounter = length(nodesWithNoEdges):-1:1 %important for it to descend here
        gEdgesMatrix(gEdgesMatrix(:,1)>nodesWithNoEdges(noEdgesCounter),1) = ...
            gEdgesMatrix(gEdgesMatrix(:,1)>nodesWithNoEdges(noEdgesCounter),1) - 1;
        gEdgesMatrix(gEdgesMatrix(:,2)>nodesWithNoEdges(noEdgesCounter),2) = ...
            gEdgesMatrix(gEdgesMatrix(:,2)>nodesWithNoEdges(noEdgesCounter),2) - 1;
        weightMatrix(:,nodesWithNoEdges(noEdgesCounter)) = [];
        weightMatrix(nodesWithNoEdges(noEdgesCounter),:) = [];
        numBondsMatrix(:,nodesWithNoEdges(noEdgesCounter)) = [];
        numBondsMatrix(nodesWithNoEdges(noEdgesCounter),:) = [];
    end
    
    
    %now, if a node is connected only by two single bonds (as long as they're not
    %connected to the same node which I can't deal with as easily), I can delete
    %it for the purposes of determining which net the graph is, because
    %it means that the node was connected to another double bond before
    %which got deleted. So I delete all nodes with weight and degree equal to 2
    %(one by one, so as to not over-delete).
    %Then, I need to make sure to transfer the edges from the nodes
    %I delete to the nodes on either side.
    [gEdgesMatrix,weightMatrix,numBondsMatrix,gEdgesMatrixWithNtds,perBaseGraph] = ...
        deleteSuperfluousNodesWithNtds_12_26_17_noGraph(gEdgesMatrix,weightMatrix,numBondsMatrix,gEdgesMatrixWithNtds,perBaseGraph);
    
    if isequal(weightMatrix,hWeightMatrix) 
        isUnchanged = true;
    end
end

%%
netTooComplex = false;
conncompG = myConnComp2(weightMatrix);

for i = 1:max(conncompG)
    %get subgraph that is disconnected from the rest
    %check to which net the subgraph is isomorphic
    weightMatrixSubgraph = weightMatrix(conncompG==i,conncompG==i);
    numBondsMatrixSubgraph = numBondsMatrix(conncompG==i,conncompG==i);
    %find edges which contain nodes in conncompG
    %gEdgesMatrixSubgraph = gEdgesMatrix(ismember(gEdgesMatrix(:,1),find(conncompG==i)) | ...
    %    ismember(gEdgesMatrix(:,2),find(conncompG==i)),:);
    gEdgesMatrixSubgraph = gEdgesMatrix(any(ismember(gEdgesMatrix(:,1:2),find(conncompG==i))')',:);
    gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtds(any(ismember(gEdgesMatrix(:,1:2),find(conncompG==i))')',:);
    numNodesGSubgraph = length(weightMatrixSubgraph);
    numEdgesGSubgraph =  size(gEdgesMatrixSubgraph,1);
    
    currentNodeNames = unique(gEdgesMatrixSubgraph(:,1:2)); %returns the current node names in sorted order
    for j = 1:numNodesGSubgraph %renumber nodes in gEdgesMatrixSubgraph to be 1 to numNodesGSubgraph
        gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2)==currentNodeNames(j)) = j;
    end
    
    if max(max(gEdgesMatrixSubgraph(:,1:2))) ~= numNodesGSubgraph || nnz(triu(weightMatrixSubgraph)) ~= numEdgesGSubgraph
        disp(strcat('I seem to have made a mistake in whichStruct = ',num2str(whichStruct)))
    end
    
    foundIsomorphism = false;
    if numNodesGSubgraph>4 %then we know none of our evaluated nets will be isomorphic to gSubgraph
        netTooComplex = true;
        foundIsomorphism = true;
        for k = 5:9 %includes 5 and 9
            for edgeIndex = 1:size(gEdgesMatrixWithNtdsSubgraph,1)
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{edgeIndex,k})-1
                    if gEdgesMatrixWithNtdsSubgraph{edgeIndex,k}(j+1) - gEdgesMatrixWithNtdsSubgraph{edgeIndex,k}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{edgeIndex,k}(j)} = 'unknown';
                        %the reason for this if statement is that when we
                        %concatenate multiple s-s regions to to form a new
                        %edge, we can get things like a closed-net-0 loop
                        %that looks like [29 30 31 32 6 7]. which
                        %represents a bulge loop.
                    end
                end
            end
        end
    end
    while ~foundIsomorphism
        %Getting the edges of a graph is slow, and isisomorphic is
        %slow, so we'll first check if gSubgraph even has the right
        %number of edges and nodes to perhaps be isomorphic, and then if it does,
        %we'll check isisomorphic
        
        %if gSubgraph is isomorphic to the net, reordernodes(gSubgraph,p) has the same structure as net.
        %else, p is an empty vector.
        %Using isisomorphic and only then doing the isomorphism doesn't
        %have any effect -- isisomorphic just calls isomorphism and
        %sees if it returns an empty array
        
        if numEdgesGSubgraph == 1 && numNodesGSubgraph == 1
            % if already took care of closedNet0 earlier in the code, can
            % comment out this part
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet0WM,numBondsMatrixSubgraph,closedNet0NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                if s0 <= minBPInHairpin-1
                    loopEntropy = loopEntropy - 10000;
                else
                    loopEntropy = loopEntropy + kB*(log(omega(s0)*vs)-(3/2)*log(s0/(beta/pi)));
                end
                componentGraphsVec(1) = componentGraphsVec(1) + 1;
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'closed_net_0';
                    end
                end
                %disp('closedNet0')
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNets0_2WM,numBondsMatrixSubgraph,closedNets0_2NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(7); %singleBond2Length;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*vs)-(3/2)*log(s0/(beta/pi)));
                loopEntropy = loopEntropy + kB*(log(omega(s1)*vs)-(3/2)*log(s1/(beta/pi)));
                if s0 <= minBPInHairpin-1 || s1 <= minBPInHairpin-1
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(2) = componentGraphsVec(2) + 1;
                for k = 6:7 %since s0 and s1 are symmetric
                    for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,k})-1
                        if gEdgesMatrixWithNtdsSubgraph{1,k}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,k}(j) == 1
                            perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,k}(j)} = 'closed_net_0';
                        end
                    end
                end
                %disp('closedNet0_2')
                break
            end
        elseif numEdgesGSubgraph == 1 && numNodesGSubgraph == 2
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet0WM,numBondsMatrixSubgraph,openNet0NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                loopEntropy = loopEntropy + kB*log(omega(s0));
                %disp('openNet0')
                componentGraphsVec(3) = componentGraphsVec(3) + 1;
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'open_net_0';
                    end
                end
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet1WM,numBondsMatrixSubgraph,openNet1NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(5); %doubleBond1Length;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1); %#ok<*UNRCH>
                end
                loopEntropy = loopEntropy + kB*(log(omega(s0)*vs)-(3/2)*log(s0/(beta/pi))-(beta*l1^2)/s0);
                if s0<=l1+ minBPInHairpin+1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(4) = componentGraphsVec(4) + 1;
                
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'open_net_1_s0';
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,5}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,5}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'open_net_1_l1';
                    end
                end
                %disp('openNet1')
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet1WM,numBondsMatrixSubgraph,closedNet1NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(7); %singleBond2Length;
                l1 = gEdgesMatrixSubgraph(5); %doubleBond1Length;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                end
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*vs^2)-(3/2)*log(s0*s1/((beta/pi)^2))-(beta*(s0 + s1)*l1^2)/(s0*s1));
                if s0<=l1+ minBPInHairpin+1 || s1 <= l1 + minBPInHairpin+1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(5) = componentGraphsVec(5) + 1;
                for k = 6:7 %since s0 and s1 are symmetric
                    for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                        if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                            perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'closed_net_1_s0';
                        end
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1 %double bond edges can't get concatenated so we don't actually need an if statement here
                    if gEdgesMatrixWithNtdsSubgraph{1,5}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,5}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'closed_net_1_l1';
                    end
                end
                %disp('closedNet1')
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet2SingleWM,numBondsMatrixSubgraph,openNet2SingleNBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(8); %singleBond3Length;
                D = s0*s1+s1*s2+s0*s2;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)-(3/2)*log(D/((beta/pi)^2)));
                %disp('openNet2Single')
                componentGraphsVec(6) = componentGraphsVec(6) + 1;
                for k = 6:8 %since s0, s1, and s2 are all symmetric
                    for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,k})-1
                        if gEdgesMatrixWithNtdsSubgraph{1,k}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,k}(j) == 1
                            perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,k}(j)} = 'open_net_2_single';
                        end
                    end
                end
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet2SingleWM,numBondsMatrixSubgraph,closedNet2SingleNBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(8); %singleBond3Length;
                s3 = gEdgesMatrixSubgraph(9); %singleBond4Length;
                D = s0*s1*s2+s0*s1*s3+s0*s2*s3+s1*s2*s3;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)-(3/2)*log(D/((beta/pi)^3)));
                %disp('closedNet2Single')
                componentGraphsVec(7) = componentGraphsVec(7) + 1;
                for k = 6:9 %since s0, s1, s2, and s3 are all symmetric
                    for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,k})-1
                        if gEdgesMatrixWithNtdsSubgraph{1,k}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,k}(j) == 1
                            perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,k}(j)} = 'closed_net_2_single';
                        end
                    end
                end
                break
            end
        elseif numEdgesGSubgraph == 5 && numNodesGSubgraph == 4
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet2aWM,numBondsMatrixSubgraph,openNet2aNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph; 
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
        
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                s2 = gEdgesMatrixSubgraph(4,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBond1Length;
                l2 = gEdgesMatrixSubgraph(5,5); %doubleBond1Length;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                    l2 = effHelixLength(l2);
                end
                D = s0*s1 + s0*s2 + s1*s2;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*(s1+s2)*l1^2)/D - (beta*(s0+s1)*l2^2)/D...
                    -log(2*beta*s1*l1*l2/D) + log(sinh(2*beta*s1*l1*l2/D)));
                %disp('openNet2a')
                componentGraphsVec(8) = componentGraphsVec(8) + 1;
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'open_net_2a_s02'; %since s0 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{4,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{4,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,6}(j)} = 'open_net_2a_s02'; %since s0 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,6}(j)} = 'open_net_2a_s1'; 
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1 %double bond edges can't get concatenated so we don't need an if statement here
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'open_net_2a_l12'; %since l1 and l2 are symmetric
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{5,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{5,5}(j)} = 'open_net_2a_l12'; %since l1 and l2 are symmetric
                end
                break
            end
        elseif numEdgesGSubgraph == 6 && numNodesGSubgraph == 4
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet2aWM,numBondsMatrixSubgraph,closedNet2aNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(4,6); %singleBond1Length;
                s2 = gEdgesMatrixSubgraph(5,6); %singleBond1Length;
                s3 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(6,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                    l2 = effHelixLength(l2);
                end
                D = s0*s3*(s1+s2) + s1*s2*(s0+s3);
                if s0*s2 ~= s1*s3
                    loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                        -(3/2)*log(D/((beta/pi)^3))-(beta*(s0+s3)*(s1+s2)*l1^2)/D - (beta*(s2+s3)*(s0+s1)*l2^2)/D...
                        -log(2*beta*(s0*s2-s1*s3)*l1*l2/D) + log(sinh(2*beta*(s0*s2-s1*s3)*l1*l2/D)));
                else %then we need to take lim x->0 of sinh(x)/x=1
                    loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                        -(3/2)*log(D/((beta/pi)^3))-(beta*(s0+s3)*(s1+s2)*l1^2)/D - (beta*(s2+s3)*(s0+s1)*l2^2)/D);
                end
                componentGraphsVec(9) = componentGraphsVec(9) + 1;
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'closed_net_2a_s03'; %since s0 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,6}(j)} = 'closed_net_2a_s03'; %since s0 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{4,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{4,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,6}(j)} = 'closed_net_2a_s12'; %since s1 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{5,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{5,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{5,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{5,6}(j)} = 'closed_net_2a_s12'; %since s1 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'closed_net_2a_l12'; %since l1 and l2 are symmetric
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{6,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{6,5}(j)} = 'closed_net_2a_l12'; %since l1 and l2 are symmetric
                end
                %disp('closedNet2a')
                break
            end
        elseif numEdgesGSubgraph == 4 && numNodesGSubgraph == 4
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet2bWM,numBondsMatrixSubgraph,openNet2bNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(1,6); %singleBond1Length;
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                    l2 = effHelixLength(l2);
                end
                D = s1*(s0+s2);
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*s1*l1^2)/D - (beta*(s0+s1+s2)*l2^2)/D...
                    -log(2*beta*s1*l1*l2/D) + log(sinh(2*beta*s1*l1*l2/D)));
                if s1<=l2 + minBPInHairpin+ 1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(10) = componentGraphsVec(10) + 1;
                %disp('openNet2b')
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'open_net_2b_s02'; %since s0 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,6}(j)} = 'open_net_2b_s02'; %since s0 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'open_net_2b_s1'; 
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,5}(j)} = 'open_net_2b_l1';
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'open_net_2b_l2';
                end
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet2bWM,numBondsMatrixSubgraph,closedNet2bNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(1,6); %singleBond1Length;
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                s3 = gEdgesMatrixSubgraph(4,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                    l2 = effHelixLength(l2);
                end
                D = s1*s3*(s0+s2);
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                    -(3/2)*log(D/((beta/pi)^3))-(beta*(s0+s2+s3)*s1*l1^2)/D - (beta*(s0+s1+s2)*s3*l2^2)/D...
                    -log(2*beta*s1*s3*l1*l2/D) + log(sinh(2*beta*s1*s3*l1*l2/D)));
                if s1<=l2 + minBPInHairpin + 1 || s3 <= l1 + minBPInHairpin + 1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(11) = componentGraphsVec(11) + 1;
                %disp('closedNet2b')
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'closed_net_2b_s02'; %since s0 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,6}(j)} = 'closed_net_2b_s02'; %since s0 and s2 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'closed_net_2b_s13'; %since s1 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{4,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{4,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,6}(j)} = 'closed_net_2b_s13'; %since s1 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,5}(j)} = 'closed_net_2b_l12'; %since l1 and l2 are symmetric
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'closed_net_2b_l12'; %since l1 and l2 are symmetric
                end
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet2cWM,numBondsMatrixSubgraph,openNet2cNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(2,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                    l2 = effHelixLength(l2);
                end
                D = s0*s1+s1*s2+s0*s2;
                A = (s0+s1)/D;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*A*l1^2) - (beta*A*l2^2)...
                    -log(2*beta*A*l1*l2) + log(sinh(2*beta*A*l1*l2)));
                componentGraphsVec(12) = componentGraphsVec(12) + 1;
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'open_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,7})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,7}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,7}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,7}(j)} = 'open_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,6}(j)} = 'open_net_2c_s2'; 
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'open_net_2c_l12'; %since l1 and l2 are symmetric
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,5}(j)} = 'open_net_2c_l12'; %since l1 and l2 are symmetric
                end
                %disp('openNet2c')
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet2cWM,numBondsMatrixSubgraph,closedNet2cNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(2,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                s3 = gEdgesMatrixSubgraph(3,7); %singleBond2Length;
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                    l2 = effHelixLength(l2);
                end
                D = s0*s1*s2+s0*s1*s3+s0*s2*s3+s1*s2*s3;
                A = (s0+s1)*(s2+s3)/D;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                    -(3/2)*log(D/((beta/pi)^3))-(beta*A*l1^2) - (beta*A*l2^2)...
                    -log(2*beta*A*l1*l2) + log(sinh(2*beta*A*l1*l2)));
                componentGraphsVec(13) = componentGraphsVec(13) + 1;
                %disp('closedNet2c')
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'closed_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,7})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,7}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,7}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,7}(j)} = 'closed_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,6}(j)} = 'closed_net_2c_s23'; %since s2 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,7})-1
                    if gEdgesMatrixWithNtdsSubgraph{3,7}(j+1) - gEdgesMatrixWithNtdsSubgraph{3,7}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,7}(j)} = 'closed_net_2c_s23'; %since s2 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,5}(j)} = 'closed_net_2c_l12'; %since l1 and l2 are symmetric
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{4,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{4,5}(j)} = 'closed_net_2c_l12'; %since l1 and l2 are symmetric
                end
                break
            end
        elseif numEdgesGSubgraph == 3 && numNodesGSubgraph == 3
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet2cL1ZeroWM,numBondsMatrixSubgraph,openNet2cL1ZeroNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(1,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(1,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(3,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                end
                D = s0*s1+s1*s2+s0*s2;
                A = (s0+s1)/D;
                if s0+s2 <= minBPInHairpin-1 || s1 + s2 <= minBPInHairpin-1
                    loopEntropy = loopEntropy - 10000;
                end
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*A*l1^2));
                componentGraphsVec(14) = componentGraphsVec(14) + 1;
                %disp('openNet2cL1Zero')
                
                %keep this as if it was an open net 2c since really it is
                %(just with l2 = 0)
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'open_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,7})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,7}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,7}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,7}(j)} = 'open_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'open_net_2c_s2'; 
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,5}(j)} = 'open_net_2c_l12';
                end
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet2cL1ZeroWM,numBondsMatrixSubgraph,closedNet2cL1ZeroNBM);
            if ~isempty(p)
                for nodeOrdering = 1:numNodesGSubgraph %here we need to actually reorder the nodes because we have more than 1 edge.
                    gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2) == p(nodeOrdering)) = nodeOrdering+numNodesGSubgraph;
                end
                gEdgesMatrixSubgraph(:,1:2) = gEdgesMatrixSubgraph(:,1:2) - numNodesGSubgraph;
                %after reordering nodes, we need to sort the matrix.
                gEdgesMatrixSubgraph(:,1:2) = sort(gEdgesMatrixSubgraph(:,1:2),2); %makes each node have the lower node in column 1 and the higher in column 2
                [gEdgesMatrixSubgraph,sortOrder] = sortrows(gEdgesMatrixSubgraph); %sort the rows
                gEdgesMatrixWithNtdsSubgraph = gEdgesMatrixWithNtdsSubgraph(sortOrder,:); 
                
                s0 = gEdgesMatrixSubgraph(1,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(1,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s3 = gEdgesMatrixSubgraph(2,7); %singleBond2Length;
                l1 = gEdgesMatrixSubgraph(3,5); %doubleBondLength;
                
                if useEffectiveHelixLengths
                    l1 = effHelixLength(l1);
                end
                D = s0*s1*s2+s0*s1*s3+s0*s2*s3+s1*s2*s3;
                A = (s0+s1)*(s2+s3)/D;
                if s0+s2 <= minBPInHairpin-1 || s0 + s3 <= minBPInHairpin -1 || ...
                        s1+s2 <= minBPInHairpin-1 || s1 + s3 <= minBPInHairpin-1
                    loopEntropy = loopEntropy - 10000;
                end
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                    -(3/2)*log(D/((beta/pi)^3))-(beta*A*l1^2));
                componentGraphsVec(15) = componentGraphsVec(15) + 1;
                %disp('closedNet2cL1Zero')
                %keep this as if it was an closed net 2c since really it is
                %(just with l2 = 0)
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,6}(j)} = 'closed_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{1,7})-1
                    if gEdgesMatrixWithNtdsSubgraph{1,7}(j+1) - gEdgesMatrixWithNtdsSubgraph{1,7}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{1,7}(j)} = 'closed_net_2c_s01'; %since s0 and s1 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,6})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,6}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,6}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,6}(j)} = 'closed_net_2c_s23'; %since s2 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{2,7})-1
                    if gEdgesMatrixWithNtdsSubgraph{2,7}(j+1) - gEdgesMatrixWithNtdsSubgraph{2,7}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{2,7}(j)} = 'closed_net_2c_s23'; %since s2 and s3 are symmetric
                    end
                end
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{3,5})-1
                    perBaseGraph{gEdgesMatrixWithNtdsSubgraph{3,5}(j)} = 'closed_net_2c_l12'; %since l1 and l2 are symmetric
                end
                break
            end
        end
        foundIsomorphism = true;
        netTooComplex = true;
    end
    if netTooComplex
        loopEntropy = loopEntropy - 10000;
        for l = 1:numEdgesGSubgraph
            for k = 5:9 
                for j = 1:length(gEdgesMatrixWithNtdsSubgraph{l,k})-1
                    if gEdgesMatrixWithNtdsSubgraph{l,k}(j+1) - gEdgesMatrixWithNtdsSubgraph{l,k}(j) == 1
                        perBaseGraph{gEdgesMatrixWithNtdsSubgraph{l,k}(j)} = 'unknown';
                    end
                end
            end
        end
        %             gSubgraph.Edges;
        %disp(strcat('net too complex for structure ',num2str(whichStruct)))
        
        %             if nnz(gSubgraph.Edges.doubleBondLength)>1 || ... %if there is more than one double bond in the subgraph
        %                  length(gSubgraph.Edges.doubleBondLength)<3
        %                 disp(gSubgraph.Edges)
        %             end
        %             allLoopEntropies(whichStruct) = loopEntropy;
        %             loopEnergy = -T*loopEntropy;
        %             allFreeEnergies(whichStruct) = bondEnergy+loopEnergy;
        break %no need to evaluate the rest of the structure.
    end
end


%for double bonds, need to make both strands have the same entry in perBaseGraph
for i = 1:length(structure)
    nBP = length(structure{i})/2;
    for j = 1:nBP-1
        if ~isempty(perBaseGraph{min(structure{i}(j),structure{i}(j+1))})
            perBaseGraph{min(structure{i}(j+nBP),structure{i}(j+nBP+1))} =...
                perBaseGraph{min(structure{i}(j),structure{i}(j+1))};
        elseif ~isempty(perBaseGraph{min(structure{i}(j+nBP),structure{i}(j+nBP+1))})
            perBaseGraph{min(structure{i}(j),structure{i}(j+1))} =...
                perBaseGraph{min(structure{i}(j+nBP),structure{i}(j+nBP+1))};
        else
            disp('error in perBaseGraph')
        end
    end
end

