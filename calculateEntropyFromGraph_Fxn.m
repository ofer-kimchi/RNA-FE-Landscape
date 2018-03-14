function [loopEntropy,componentGraphs] = calculateEntropyFromGraph_Fxn(gEdgesMatrix,omega,minBPInHairpin,weightMatrix,numBondsMatrix,allPerms,vs)
%%
%given a graph of a structure, calculate the entropy. To do so, we
%decompose the graph into the minimal sub-graphs of which it's comprised.

b = 0.8/(0.33); %in units of nucleotides (approximate size of ntd is 3.3 angstroms)
    %This assumes a persistence length of RNA of ~0.8 nm found in The Snakelike Chain Character of Unstructured RNA
    %by  David R. Jacobson, Dustin B. McIntosh, and Omar A. Saleh
    %This agrees well with Siggia's formula b=2.5a (where a is the size of one ntd).
beta = 3/(2*b);

kB = 0.0019872; %units of kcal/(mol*Kelvin)


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

componentGraphsVec = zeros(1,15); %a vector telling you how many of each type of net are found in this graph
%------------------------------------------------------------------------
%%

numEdges = size(gEdgesMatrix,1);
listEdges = 1:numEdges;
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
        end
        weightMatrix = hWeightMatrix; %remove bond 
        numBondsMatrix(gEdgesMatrix(edgeIndex,1),gEdgesMatrix(edgeIndex,2)) = 0;
        numBondsMatrix(gEdgesMatrix(edgeIndex,2),gEdgesMatrix(edgeIndex,1)) = 0;
        gEdgesMatrix(edgeIndex,:) = []; %update the edge table as well
        numDisconnectedGraphs = potentialNumDisconnectedGraphs;
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
        %add the new edge to gEdges
        gEdgesMatrix = sortrows([gEdgesMatrix;edgeToAdd]);
        weightMatrix(numNodesG, numNodesG) = 1; %add the new node and edge to weightMatrix and numBondsMatrix
        numBondsMatrix(numNodesG, numNodesG) = 1;
        
        %now get rid of the old edge connecting the node to itself
        weightMatrix(gEdgesMatrix(edgeToRemove,1),gEdgesMatrix(edgeToRemove,1)) = 0; 
            %this syntax works because the node is connected to itself
        numBondsMatrix(gEdgesMatrix(edgeToRemove,1),gEdgesMatrix(edgeToRemove,1)) = 0; 
        gEdgesMatrix(edgeToRemove,:) = [];  %make sure to do this after deleting the node from weightMatrix and numBondsMatrix!
        
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
    [gEdgesMatrix,weightMatrix,numBondsMatrix]=deleteSuperfluousNodes_Fxn(gEdgesMatrix,weightMatrix,numBondsMatrix);
    
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
    gEdgesMatrixSubgraph = gEdgesMatrix(any(ismember(gEdgesMatrix(:,1:2),find(conncompG==i))')',:);
    
    numNodesGSubgraph = length(weightMatrixSubgraph);
    numEdgesGSubgraph =  size(gEdgesMatrixSubgraph,1);
    
    currentNodeNames = unique(gEdgesMatrixSubgraph(:,1:2)); %returns the current node names in sorted order
    for j = 1:numNodesGSubgraph %renumber nodes in gEdgesMatrixSubgraph to be 1 to numNodesGSubgraph
        gEdgesMatrixSubgraph(gEdgesMatrixSubgraph(:,1:2)==currentNodeNames(j)) = j;
    end
    
    
    foundIsomorphism = false;
    if numNodesGSubgraph>4 %then we know none of our evaluated nets will be isomorphic to gSubgraph
        netTooComplex = true;
        foundIsomorphism = true;
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
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,openNet1WM,numBondsMatrixSubgraph,openNet1NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(5); %doubleBond1Length;
                
                loopEntropy = loopEntropy + kB*(log(omega(s0)*vs)-(3/2)*log(s0/(beta/pi))-(beta*l1^2)/s0);
                if s0<=l1+ minBPInHairpin+1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(4) = componentGraphsVec(4) + 1;
                %disp('openNet1')
                break
            end
            p = myPseudoIsomorphism(allPerms{numNodesGSubgraph},weightMatrixSubgraph,closedNet1WM,numBondsMatrixSubgraph,closedNet1NBM);
            if ~isempty(p)
                s0 = gEdgesMatrixSubgraph(6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(7); %singleBond2Length;
                l1 = gEdgesMatrixSubgraph(5); %doubleBond1Length;
                
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*vs^2)-(3/2)*log(s0*s1/((beta/pi)^2))-(beta*(s0 + s1)*l1^2)/(s0*s1));
                if s0<=l1+ minBPInHairpin+1 || s1 <= l1 + minBPInHairpin+1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(5) = componentGraphsVec(5) + 1;
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length. connects nodes 0 & 2. corresponds to our s3
                s1 = gEdgesMatrixSubgraph(3,6); %singleBond1Length. connects nodes 1 & 2. corresponds to our s2
                s2 = gEdgesMatrixSubgraph(4,6); %singleBond1Length. connects nodes 1 & 3. corresponds to our s1
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBond1Length;
                l2 = gEdgesMatrixSubgraph(5,5); %doubleBond1Length;
                
                D = s0*s1 + s0*s2 + s1*s2;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*(s1+s2)*l1^2)/D - (beta*(s0+s1)*l2^2)/D...
                    -log(2*beta*s1*l1*l2/D) + log(sinh(2*beta*s1*l1*l2/D)));
                %disp('openNet2a')
                componentGraphsVec(8) = componentGraphsVec(8) + 1;
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length. connects nodes 0 & 2, corresponds to our s3
                s1 = gEdgesMatrixSubgraph(4,6); %singleBond1Length. connects nodes 1 & 2, corresponds to our s2
                s2 = gEdgesMatrixSubgraph(5,6); %singleBond1Length. connects nodes 1 & 3, corresponds to our s1
                s3 = gEdgesMatrixSubgraph(3,6); %singleBond1Length. connects nodes 0 & 3, corresponds to our s4
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(6,5); %doubleBondLength;
                
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length; connects nodes 0 & 2, corresponds to our s1
                s1 = gEdgesMatrixSubgraph(1,6); %singleBond1Length; connects nodes 2 & 3, corresponds to our s2
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length; connects nodes 1 & 3, corresponds to our s3
                l1 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                
                D = s1*(s0+s2);
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*s1*l1^2)/D - (beta*(s0+s1+s2)*l2^2)/D...
                    -log(2*beta*s1*l1*l2/D) + log(sinh(2*beta*s1*l1*l2/D)));
                if s1<=l2 + minBPInHairpin+ 1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(10) = componentGraphsVec(10) + 1;
                %disp('openNet2b')
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length; connects nodes 0 & 2, corresponds to our s1
                s1 = gEdgesMatrixSubgraph(1,6); %singleBond1Length; connects nodes 2 & 3, corresponds to our s2
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length; connects nodes 1 & 3, corresponds to our s3
                s3 = gEdgesMatrixSubgraph(4,6); %singleBond1Length; connects nodes 0 & 1, corresponds to our s4
                l1 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                
                D = s1*s3*(s0+s2);
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                    -(3/2)*log(D/((beta/pi)^3))-(beta*(s0+s2+s3)*s1*l1^2)/D - (beta*(s0+s1+s2)*s3*l2^2)/D...
                    -log(2*beta*s1*s3*l1*l2/D) + log(sinh(2*beta*s1*s3*l1*l2/D)));
                if s1<=l2 + minBPInHairpin + 1 || s3 <= l1 + minBPInHairpin + 1 %guess that in order to make a parallel loop it needs to be this long.
                    loopEntropy = loopEntropy - 10000;
                end
                componentGraphsVec(11) = componentGraphsVec(11) + 1;
                %disp('closedNet2b')
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(2,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;
                
                D = s0*s1+s1*s2+s0*s2;
                A = (s0+s1)/D;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*A*l1^2) - (beta*A*l2^2)...
                    -log(2*beta*A*l1*l2) + log(sinh(2*beta*A*l1*l2)));
                componentGraphsVec(12) = componentGraphsVec(12) + 1;
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(2,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(3,6); %singleBond1Length;
                s3 = gEdgesMatrixSubgraph(3,7); %singleBond2Length;
                l1 = gEdgesMatrixSubgraph(1,5); %doubleBondLength;
                l2 = gEdgesMatrixSubgraph(4,5); %doubleBondLength;

                D = s0*s1*s2+s0*s1*s3+s0*s2*s3+s1*s2*s3;
                A = (s0+s1)*(s2+s3)/D;
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*omega(s3)*vs^3)...
                    -(3/2)*log(D/((beta/pi)^3))-(beta*A*l1^2) - (beta*A*l2^2)...
                    -log(2*beta*A*l1*l2) + log(sinh(2*beta*A*l1*l2)));
                componentGraphsVec(13) = componentGraphsVec(13) + 1;
                %disp('closedNet2c')
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(1,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(1,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                l1 = gEdgesMatrixSubgraph(3,5); %doubleBondLength;
                
                D = s0*s1+s1*s2+s0*s2;
                A = (s0+s1)/D;
                if s0+s2 <= minBPInHairpin-1 || s1 + s2 <= minBPInHairpin-1
                    loopEntropy = loopEntropy - 10000;
                end
                loopEntropy = loopEntropy + kB*(log(omega(s0)*omega(s1)*omega(s2)*vs^2)...
                    -(3/2)*log(D/((beta/pi)^2))-(beta*A*l1^2));
                componentGraphsVec(14) = componentGraphsVec(14) + 1;
                %disp('openNet2cL1Zero')
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
                gEdgesMatrixSubgraph = sortrows(gEdgesMatrixSubgraph); %sort the rows
                
                s0 = gEdgesMatrixSubgraph(1,6); %singleBond1Length;
                s1 = gEdgesMatrixSubgraph(1,7); %singleBond2Length;
                s2 = gEdgesMatrixSubgraph(2,6); %singleBond1Length;
                s3 = gEdgesMatrixSubgraph(2,7); %singleBond2Length;
                l1 = gEdgesMatrixSubgraph(3,5); %doubleBondLength;
                
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
                break
            end
        end
        foundIsomorphism = true;
        netTooComplex = true;
    end
    if netTooComplex
        loopEntropy = loopEntropy - 10000; %disallow all structures with topologies more complex than we can deal with by including this entropy penalty
        break %no need to evaluate the rest of the structure.
    end
end

if nargout > 1
    componentGraphs = componentGraphsVec;
end