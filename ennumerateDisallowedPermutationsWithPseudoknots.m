function[disallowedPermutations] = ennumerateDisallowedPermutationsWithPseudoknots(threeOrFour)
%threeOrFour can be 3 if we're making the three-way compatibility matrix,
%or 4 if we're making the four-way compatibility matrix.

if threeOrFour == 4
    sequence = ['a' 'b' 'c' 'd' 'e' 'f' 'q' 'r']; %this is our structure. 
    %Each index defines a position on the RNA strand. Index a is bonded to b, 
    %c to d, e to f, and q to r. These pairs are ordered such that a<c<e<q.
    %(we don't use g and h to not confuse with the graphs). 
    %Also, a<b, c<d, e<f, and q<r.
    %In this code, we allow any combination of the threeOrFour nodes to
    %represent a pseudoknot.
    %We want to find which permutations of this set of four base pairs is not
    %allowed by our rules of no order-3 or higher pseudoknots. 
    %(By permutation, we mean the relations between b, d, f, and r (which
    %is greater than which, etc.)
    %basically, we go through the calculateEntropyFromGraph_Fxn code and
    %try to decompose all possible permutations of these bps into minimal
    %subgraphs. If we can't decompose a permutation into a graph that we
    %can calculate the topology for without numerical integration, we say
    %that this permutation isn't allowed.
elseif threeOrFour == 3
    sequence = ['a' 'b' 'c' 'd' 'e' 'f'];
end

%tic

%define graphs for the nets we know of

openNet0 = graph([[0,1];[1,0]]); %simplest version of open-net-0. Just two connected nodes.
openNet0.Edges.numBonds = ones(length(openNet0.Edges.Weight),1); 

closedNet0 = graph(1); %simplest version of closed-net-0. Just a node connected to itself.
closedNet0.Edges.numBonds = ones(length(closedNet0.Edges.Weight),1); 
 
closedNets0_2 = graph(2); %2 closed nets zero on the same node
closedNets0_2.Edges.numBonds = 2*ones(length(closedNet0.Edges.Weight),1); 

openNet1 = graph([[0,3];[3,0]]); %like closed net 0 and open net 0, an openNet1 has two nodes connected by an edge
    %To distinguish openNet1 from closedNet0 and openNet0, we have the
    %edge have weight 3 with numBonds =2 because there is one double bond
    %and one single bond.
openNet1.Edges.numBonds = 2.*ones(length(openNet1.Edges.Weight),1); 

closedNet1 = graph([[0,4];[4,0]]); %closedNet1 has two single bonds and a double bond. 
    %Hence, the edge between the two nodes has weight 4 and numBonds =3.
closedNet1.Edges.numBonds = 3.*ones(length(closedNet1.Edges.Weight),1); 

openNet2a = graph([[0 2 1 0];[2 0 1 1];[1 1 0 2];[0 1 2 0]]); 
openNet2a.Edges.numBonds = ones(length(openNet2a.Edges.Weight),1); 

closedNet2a = graph([[0 2 1 1];[2 0 1 1];[1 1 0 2];[1 1 2 0]]); 
closedNet2a.Edges.numBonds = ones(length(closedNet2a.Edges.Weight),1); 

openNet2b = graph([[0 3 1 0];[3 0 0 1];[1 0 0 2];[0 1 2 0]]); 
openNet2b.Edges.numBonds = [2 1 1 1]'; 

closedNet2b = graph([[0 3 1 0];[3 0 0 1];[1 0 0 3];[0 1 3 0]]); 
closedNet2b.Edges.numBonds = [2 1 1 2]'; 

openNet2c = graph([[0 2 2 0];[2 0 0 1];[2 0 0 2];[0 1 2 0]]);
openNet2c.Edges.numBonds = [1 2 1 1]';

closedNet2c = graph([[0 2 2 0];[2 0 0 2];[2 0 0 2];[0 2 2 0]]);
closedNet2c.Edges.numBonds = [1 2 2 1]';

openNet2Single = graph([[0,3];[3,0]]); %like openNet1 but the double bond is replaced by two single bonds
%more relevantly, like an openNet2a or openNet2c (the expressions yield the same result) with l1 = l2 = 0.
openNet2Single.Edges.numBonds = 3.*ones(length(openNet2Single.Edges.Weight),1); 

closedNet2Single = graph([[0,4];[4,0]]); %like closedNet1 but the double bond is replaced by two single bonds
%more relevantly, like a closedNet2a or closedNed2c (the expressions yield the same result) with l1 = l2 = 0.
closedNet2Single.Edges.numBonds = 4.*ones(length(closedNet2Single.Edges.Weight),1);

openNet2cL1Zero = graph([[0 2 1];[2 0 2];[1 2 0]]); %like openNet2c or openNet2a with l1 = 0 (and l2 nonzero).
openNet2cL1Zero.Edges.numBonds = [2 1 1]';

closedNet2cL1Zero = graph([[0 2 2];[2 0 2];[2 2 0]]); %like openNet2c or openNet2a with l1 = 0 (and l2 nonzero).
closedNet2cL1Zero.Edges.numBonds = [2 2 1]';


complexNets = zeros(1,factorial(length(sequence))); %list of permutations
%which lead to high order pseudoknots.
disallowedPermutations = zeros(factorial(length(sequence))*27,4*threeOrFour);
netCounter = 1;

if threeOrFour == 4
    qrPseudoknotPossibilities = [0,1];
else
    qrPseudoknotPossibilities = 0;
end

seqPerms = perms(sequence);
for perm = 1:factorial(length(sequence))
for abPseudoknot = [0,1]
for cdPseudoknot = [0,1]
for efPseudoknot = [0,1]
for qrPseudoknot = qrPseudoknotPossibilities
for abNonPseudoknot = [0,1]
for cdNonPseudoknot = [0,1]
for efNonPseudoknot = [0,1]
for qrNonPseudoknot = qrPseudoknotPossibilities
if ~(abPseudoknot && abNonPseudoknot)
if ~(cdPseudoknot && cdNonPseudoknot)
if ~(efPseudoknot && efNonPseudoknot)
if ~(qrPseudoknot && qrNonPseudoknot)
orderedSeq = seqPerms(perm,:);
%make the sequence into single bonds
if threeOrFour == 4
    orderedSeq = strcat(orderedSeq(1),'xxxxx',orderedSeq(2),'xxxxx',orderedSeq(3),'xxxxx',...
        orderedSeq(4),'xxxxx',orderedSeq(5),'xxxxx',orderedSeq(6),'xxxxx',orderedSeq(7),'xxxxx',orderedSeq(8),'xxxxx'); 
elseif threeOrFour == 3
    orderedSeq = strcat(orderedSeq(1),'xxxxx',orderedSeq(2),'xxxxx',orderedSeq(3),'xxxxx',...
        orderedSeq(4),'xxxxx',orderedSeq(5),'xxxxx',orderedSeq(6),'xxxxx'); 
end
a = strfind(orderedSeq,'a'); c = strfind(orderedSeq,'c'); 
if a < c
e = strfind(orderedSeq,'e'); 
if c < e
if threeOrFour == 4
    q = strfind(orderedSeq,'q'); 
elseif threeOrFour == 3
    q = length(orderedSeq) + 10;
end
if e < q
b = strfind(orderedSeq,'b'); 
if a < b
d = strfind(orderedSeq,'d');
if c < d
f = strfind(orderedSeq,'f'); 
if e < f
if threeOrFour == 4
    r = strfind(orderedSeq,'r');
elseif threeOrFour == 3
    r = length(orderedSeq) + 20;
end
if q < r
    if abPseudoknot
        structureAB = [a, a+1, b, b+1];
    elseif abNonPseudoknot
        structureAB = [a, a+1, b, b-1];
    else 
        structureAB = [a,b];
    end
    if cdPseudoknot
        structureCD = [c, c+1, d, d+1];
    elseif cdNonPseudoknot
        structureCD = [c, c+1, d, d-1];
    else 
        structureCD = [c,d];
    end
    if efPseudoknot
        structureEF = [e, e+1, f, f+1];
    elseif efNonPseudoknot
        structureEF = [e, e+1, f, f-1];
    else 
        structureEF = [e, f];
    end
    if qrPseudoknot
        structureQR = [q, q+1, r, r+1];
    elseif qrNonPseudoknot
        structureQR = [q, q+1, r, r-1];
    else
        structureQR = [q,r];
    end
    structure = {structureAB,structureCD,structureEF,structureQR};
    numBP = max([a,b,c,d,e,f,q,r]);
    strucMatrix = zeros(numBP); %connectivity matrix of structure (strucMatrix_i,j = 1 if ntds i and j are bonded)
    
    %set up structure and strucMatrix
    for i = 1:length(structure)
        numBonds = length(structure{i})/2;
        for j = 1:numBonds
            k = structure{i}(j);
            l = structure{i}(j+numBonds);
            strucMatrix(k,l) = 1;
            strucMatrix(l,k) = 1;
        end
    end
    %set up cell array for lists of unpaired regions
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
                if ~ismember([i;j]',doubleEdges','rows') 
                    numPseudoSingleEdges = numPseudoSingleEdges + 1;
                    pseudoSingleEdges(:,numPseudoSingleEdges) = [i;j];
                end
            end
        end
    end
    pseudoSingleEdges(:,numPseudoSingleEdges+1:end) = [];
    
    
    %find single-strand edges of graph
    numSingleEdges = length(structUnpaired);
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
        if singleEdgesSeq(2,i) == numBP+1 %then it shouldn't actually be 0 it should be 1.
            singleEdgesSeq(2,i) = numBP;
        end
    end
    
    
    
    %create EdgeTable to define graph
    
    EdgeTable = unique(table([doubleEdges,singleEdges,pseudoSingleEdges]','VariableNames',{'EndNodes'})); 
    %list of all edges (no specifications yet as to weights or numBonds or lengths)
    numEdges = length(EdgeTable.EndNodes(:,1));
    EdgeTable.Weight = zeros(numEdges,1); %weight of an edge -- sum of the number of bonds composing an edge, with double bonds counted as 2 and single bonds as 1.
    EdgeTable.numBonds = zeros(numEdges,1); %number of bonds in a given edge
    EdgeTable.doubleBondLength = zeros(numEdges,1); %the length of a double bond between the nodes (if one exists)
    EdgeTable.singleBond1Length = zeros(numEdges,1); %the length of a single bond between the nodes (if one exists)
    EdgeTable.singleBond2Length = zeros(numEdges,1); %the length of a second single bond between the nodes (if one exists)
    EdgeTable.singleBond3Length = zeros(numEdges,1); %the length of a third single bond between the nodes (if one exists)
    EdgeTable.singleBond4Length = zeros(numEdges,1); %the length of a fourth single bond between the nodes (if one exists)
    
    for i = 1:numSingleEdges
        [~, index] = ismember(singleEdges(:,i)',EdgeTable.EndNodes,'rows'); %we know that this edge is a member, just want to find its index.
        EdgeTable.Weight(index) = EdgeTable.Weight(index) + 1;
        EdgeTable.numBonds(index) = EdgeTable.numBonds(index) + 1;
        if EdgeTable.singleBond1Length(index) == 0
            EdgeTable.singleBond1Length(index) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i); %if a sinlge stranded region contains 4 ntds, its length is actually 5
        elseif EdgeTable.singleBond2Length(index) == 0
            EdgeTable.singleBond2Length(index) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
        elseif EdgeTable.singleBond3Length(index) == 0
            EdgeTable.singleBond3Length(index) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
        elseif EdgeTable.singleBond4Length(index) == 0
            EdgeTable.singleBond4Length(index) = singleEdgesSeq(2,i) - singleEdgesSeq(1,i);
        else
            disp('something is off')
        end
    end
    
    for i = 1:numDoubleEdges
        [~, index] = ismember(doubleEdges(:,i)',EdgeTable.EndNodes,'rows'); %we know that this edge is a member, just want to find its index.
        EdgeTable.Weight(index) = EdgeTable.Weight(index) + 2;
        EdgeTable.numBonds(index) = EdgeTable.numBonds(index) + 1;
        if EdgeTable.doubleBondLength(index) == 0
            EdgeTable.doubleBondLength(index) = lengthDoubleEdges(i);
        else
            disp('something is off 2')
        end
    end
    
    for i = 1:numPseudoSingleEdges
        [~, index] = ismember(pseudoSingleEdges(:,i)',EdgeTable.EndNodes,'rows'); %we know that this edge is a member, just want to find its index.
        EdgeTable.Weight(index) = EdgeTable.Weight(index) + 1;
        EdgeTable.numBonds(index) = EdgeTable.numBonds(index) + 1;
        if EdgeTable.singleBond1Length(index) == 0
            EdgeTable.singleBond1Length(index) = 1; 
        elseif EdgeTable.singleBond2Length(index) == 0
            EdgeTable.singleBond2Length(index) = 1;
        elseif EdgeTable.singleBond3Length(index) == 0
            EdgeTable.singleBond3Length(index) = 1;
        elseif EdgeTable.singleBond4Length(index) == 0
            EdgeTable.singleBond4Length(index) = 1;
        else
            disp('something is off 3')
        end
    end
    
    
    %create graph
    g = graph(EdgeTable); 

    gEdges = EdgeTable;
    numDisconnectedGraphs = 1; %number of disconnected graphs g can be separated into with the deletion of double stranded regions.
    
    %try to remove each bond , and if you can -- by
    %removing it -- separate the graph into two disconnected graphs,
    %then remove the bond. By doing this, we're trying to
    %decompose the graph into the smallest nets possible.
    
    %first try all double bonds
    for i = 1:numDoubleEdges
        edgeIndex = findedge(g,doubleEdges(1,i),doubleEdges(2,i));
        if gEdges.Weight(edgeIndex) == 2 %meaning there is only a double bond connecting these nodes
            h = rmedge(g,edgeIndex); %h is the graph when you remove double bond i from g.
            potentialNumDisconnectedGraphs = max(conncomp(h));
            if potentialNumDisconnectedGraphs>numDisconnectedGraphs %if h can be separated into more disconnected graphs than g can
                g = h;
                gEdges(edgeIndex,:) = []; %update the edge table as well
                numDisconnectedGraphs = potentialNumDisconnectedGraphs;
            end
        end
    end
    
    %then all single bonds. In this case,  we
    %will account for the entropy cost of these (open net 0)
    for i = 1:numSingleEdges
        edgeIndex = findedge(g,singleEdges(1,i),singleEdges(2,i));
        if gEdges.Weight(edgeIndex) == 1 %meaning there is only one single bond connecting these nodes
            h = rmedge(g,edgeIndex); %h is the graph when you remove single bond i from g.
            potentialNumDisconnectedGraphs = max(conncomp(h));
            if potentialNumDisconnectedGraphs>numDisconnectedGraphs %if h can be separated into more disconnected graphs than g can
                g = h;
                gEdges(edgeIndex,:) = []; %update the edge table as well
                numDisconnectedGraphs = potentialNumDisconnectedGraphs;
            end
        end
    end
    
    %and finally all pseudoSingle bonds. I don't think we should
    %include an entropy cost to doing this.
    for i = 1:numPseudoSingleEdges
        edgeIndex = findedge(g,pseudoSingleEdges(1,i),pseudoSingleEdges(2,i));
        if gEdges.Weight(edgeIndex) == 1 %meaning there is only one pseudoSingle bond connecting these nodes
            h = rmedge(g,edgeIndex); %h is the graph when you remove pseudoSingle bond i from g.
            potentialNumDisconnectedGraphs = max(conncomp(h));
            if potentialNumDisconnectedGraphs>numDisconnectedGraphs %if h can be separated into more disconnected graphs than g can
                g = h;
                gEdges(edgeIndex,:) = []; %update the edge table as well
                numDisconnectedGraphs = potentialNumDisconnectedGraphs;
            end
        end
    end
      
    
    
    isIsomorphic = false; %make sure to do the following loop at least once
    while ~isIsomorphic
        h=g; %at the end of this loop, we want to check if g has changed -- if it has, we might need to
        %rerun this loop
        
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
        

        index = find(gEdges.EndNodes(:,1) == gEdges.EndNodes(:,2)); %index of edge connecting a node to itself
        for i = 1:length(index)
            if gEdges.singleBond2Length(index) ~=0 %don't want to get rid of double hairpin loops
                index(i) = [];
            end
        end
        counter = 1;
        while counter <= length(index)
            if isempty(nearest(g,gEdges.EndNodes(index(counter),1),1,'Method','unweighted'))
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
        %make a new node that's connected to itself with this single bond
        for i = 1:length(index) %only repeat this the number of times we found hairpin loops (otherwise it'll never stop)
            edgeToRemove = index(1);
            g = addnode(g,1);
            edgeToAdd = table([numnodes(g) numnodes(g)],'VariableNames',{'EndNodes'});
            edgeToAdd.Weight = 1;
            edgeToAdd.numBonds = 1;
            edgeToAdd.doubleBondLength = 0;
            edgeToAdd.singleBond1Length = gEdges.singleBond1Length(edgeToRemove);
            edgeToAdd.singleBond2Length = 0;
            edgeToAdd.singleBond3Length = 0;
            edgeToAdd.singleBond4Length = 0;
            g = addedge(g,edgeToAdd);
            gEdges = sortrows([gEdges;edgeToAdd]);
            
            %now get rid of the old edge connecting the node to itself
            g = rmedge(g,edgeToRemove);
            gEdges(edgeToRemove,:) = [];
            index = find(gEdges.EndNodes(:,1) == gEdges.EndNodes(:,2));
        end
        
        
        %If any ntd was only connected to the rest of a graph by open-nets-0,
        %then the node is unimportant and can (and should) be
        %deleted. We therefore remove all nodes with no edges.
        nodesWithNoEdges = find(centrality(g,'degree','Importance',gEdges.Weight)==0);
        g = rmnode(g,nodesWithNoEdges);
        %The function centrality has the property that a self-loop counts as two edges connecting to the node (as we'd like it to).
        
        %now, if a node is connected only by two single bonds (as long as they're not
        %connected to the same node which I can't deal with as easily), I can delete
        %it for the purposes of determining which net the graph is, because
        %it means that the node was connected to another double bond before
        %which got deleted. So I delete all nodes with weight and degree equal to 2
        %(one by one, so as to not over-delete).
        %Then, I need to make sure to transfer the edges from the nodes
        %I delete to the nodes on either side.
        
        %gEdges is still updated here since even though we just removed a
        %node, it was a node with no edges. We just need to update the
        %names of the nodes in the edges.
        for noEdgesCounter = length(nodesWithNoEdges):-1:1 %important for it to descend here
            gEdges.EndNodes(gEdges.EndNodes>nodesWithNoEdges(noEdgesCounter)) = ...
                gEdges.EndNodes(gEdges.EndNodes>nodesWithNoEdges(noEdgesCounter))-1;
        end
        [g,gEdges]=deleteSuperfluousNodes_withGraph_Fxn(g,gEdges); %deleting these nodes takes more code time than any other part.
        
        if max(conncomp(g)) ~= max(conncomp(h)) %MatLab has a bug in 
            %versions before 2017a which make isisomorphic throw an error if
            %one of the input graphs has exactly one connected component, and the other has more.
          isIsomorphic = false;
        else
          isIsomorphic = isisomorphic(g, h);
        end
        
    end
    g=h;
    
    
    netTooComplex = false;
    netTooLarge = false;
    for i = 1:max(conncomp(g)) 
        gSubgraph = subgraph(g,find(conncomp(g)==i)); %subgraph that is disconnected from the rest
        %check which net gSubgraph is isomorphic to 
        foundIsomorphism = false;
        while ~foundIsomorphism
            %Getting the edges of a graph is slow, and isisomorphic is
            %slow, so we'll first check if gSubgraph even has the right
            %number of edgesand nodes to perhaps be isomorphic, and then if it does,
            %we'll check isisomorphic
            
            %if gSubgraph is isomorphic to the net, reordernodes(gSubgraph,p) has the same structure as net.
            %else, p is an empty vector.
            %Using isisomorphic and only then doing the isomorphism doesn't
            %have any effect -- isisomorphic just calls isomorphism and
            %sees if it returns an empty array
            
            if numedges(gSubgraph) == 1 && numnodes(gSubgraph) == 1
                % if already took care of closedNet0 earlier in the code, can
                % comment out this part
                p = isomorphism(closedNet0,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(closedNets0_2,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
            elseif numedges(gSubgraph) == 1 && numnodes(gSubgraph) == 2
                p = isomorphism(openNet0,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(openNet1,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(closedNet1,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(openNet2Single,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(closedNet2Single,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
            elseif numedges(gSubgraph) == 5 && numnodes(gSubgraph) == 4
                p = isomorphism(openNet2a,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
            elseif numedges(gSubgraph) == 6 && numnodes(gSubgraph) == 4
                p = isomorphism(closedNet2a,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
            elseif numedges(gSubgraph) == 4 && numnodes(gSubgraph) == 4    
                p = isomorphism(openNet2b,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(closedNet2b,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(openNet2c,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(closedNet2c,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
            elseif numedges(gSubgraph) == 3 && numnodes(gSubgraph) == 3
                p = isomorphism(openNet2cL1Zero,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
                p = isomorphism(closedNet2cL1Zero,gSubgraph,'EdgeVariables',{'numBonds','Weight'});
                if ~isempty(p)
                    break
                end
            end
            foundIsomorphism = true;
            netTooComplex = true;
        end
        if netTooComplex
            complexNets(netCounter) = perm;
            if threeOrFour == 4
                [~,index] = sort([a,b,c,d,e,f,q,r]);
                pseudoknotList = [abPseudoknot,abNonPseudoknot,...
                    cdPseudoknot,cdNonPseudoknot,efPseudoknot,efNonPseudoknot,qrPseudoknot,qrNonPseudoknot];
            elseif threeOrFour == 3
                [~,index] = sort([a,b,c,d,e,f]);
                pseudoknotList = [abPseudoknot,abNonPseudoknot,...
                    cdPseudoknot,cdNonPseudoknot,efPseudoknot,efNonPseudoknot];
            end
            disallowedPermutations(netCounter,:) = [index,pseudoknotList];
            netCounter = netCounter + 1;
            break %no need to evaluate the rest of the structure.
        end
    end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end
end


complexNets(netCounter:end) = [];
disallowedPermutations(netCounter:end,:) = [];
netCounter = netCounter - 1;
%toc