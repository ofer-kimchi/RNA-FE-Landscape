function[gEdgesMatrix,weightMatrix,numBondsMatrix] = ...
    deleteSuperfluousNodes_Fxn(gEdgesMatrix,weightMatrix,numBondsMatrix)
%if a node is connected only by two single bonds (as long as they're not
    %connected to the same node which I can't deal with as easily), I can delete
    %it for the purposes of determining which net the graph is, because
    %it means that the node was connected to another double bond before
    %which got deleted. So I delete all nodes with weight and degree equal to 2
    %(one by one, so as to not over-delete).
    %Then, I need to make sure to transfer the edges from the nodes
    %I delete to the nodes on either side.


%first we consider an edge case -- nodes which have four single bonds
%attached. There are two cases with these nodes -- 1) either all their
%neighbors are connected to one another (without them); or, 2) two of their
%neighbors are connected to each other, and the other two are as well. This
%is because if three of their neighbors are connected to each other and the
%other is not connected, then the bond between that neighbor and the node
%in question would have already been broken in order to make a disconnected
%graph. We are interested in case 2: in this case, we can delete the node,
%and make the graph disconnected.

%%
possibleNodesToDelete = find(sum(weightMatrix)==4 & sum(numBondsMatrix)==4)';
while ~isempty(possibleNodesToDelete)
    for i = 1:length(possibleNodesToDelete)
        node = possibleNodesToDelete(i);
        hWeightMatrix = weightMatrix; 
        hWeightMatrix(node,:) = [];
        hWeightMatrix(:,node) = [];
        if max(myConnComp2(hWeightMatrix)) > max(myConnComp2(weightMatrix)) %then deleting that node disconnects the graph
            %Delete this node. If the node is connected to nodes
            %a, b, c, d, and when it's deleted the graph splits into a,b and c,d, then
            %make a connection between a and b whose length is the sum of the lengths
            %of the connections between a and node, and b and node (and same for c,d).
            allConnectedNodes = find(weightMatrix(node,:)~=0)';
            
            if length(allConnectedNodes)<3 %don't need to consider these cases because they're already handled by the 
                %``index'' part of the function in the for loop preceding
                %the call of deleteSuperfluousNodes()
                possibleNodesToDelete(i) = [];
                break
            end
            
            if length(allConnectedNodes) == 3
                %then it must be connected to one node twice (if it was
                %connected to one node only once, and by removing it you
                %disconnect the graph, then the removal of that edge would
                %have disconnected the graph and we would have gotten rid
                %of it previously)
                
                %we first find which node that is, then take care of it,
                %then take care of the other two.
                
                for j = 1:3 %to find which node is connected twice
                    edgesOfNode = [allConnectedNodes(j) node];
                    %this line (below) is equivalent to: [~,index1] = ismember(edgesOfNode(1,:),gEdgesMatrix(:,1:2),'rows');
                    index1 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(1,:))'));
                    if ~isempty(index1) && gEdgesMatrix(index1,7) ~= 0 %singleBond2Length %then this is the node that's connected twice
                        break
                    end
                    %this line (below) is equivalent to: [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdgesMatrix(:,1:2),'rows');
                    index1 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(1,:)))'));
                    if ~isempty(index1) && gEdgesMatrix(index1,7) ~= 0 %singleBond2Length %then this is the node that's connected twice
                        break
                    end
                    if j == 3 %if we've reached this point at j ==3 then our assumptions are wrong
                        disp('our assumptions are wrong')
                    end
                end
                
                connectedNodes = allConnectedNodes(j); %this is the node that's connected twice
                edgeToAddLength = gEdgesMatrix(index1,6) + gEdgesMatrix(index1,7)... %singleBond1Length + singleBond2Length
                    + gEdgesMatrix(index1,8) + gEdgesMatrix(index1,9); %singleBond3Length + singleBond4Length
                if ~any(all((gEdgesMatrix(:,1:2) == [connectedNodes connectedNodes])')) 
                    %this line (above) is equivalent to: ~ismember([connectedNodes connectedNodes], gEdgesMatrix(:,1:2),'rows')
                    edgeToAdd = [connectedNodes, connectedNodes, 1, 1, 0, edgeToAddLength, 0, 0, 0];
                    gEdgesMatrix = sortrows([gEdgesMatrix;edgeToAdd]);
                    weightMatrix(connectedNodes,connectedNodes) = 1;
                    numBondsMatrix(connectedNodes,connectedNodes) = 1;
                else
                    %this line (below) is equivalent to: [~,index1] = ismember([connectedNodes connectedNodes],gEdgesMatrix(:,1:2),'rows');
                    index1 = find(all((gEdgesMatrix(:,1:2) == [connectedNodes connectedNodes])'));
                    gEdgesMatrix(index1,3) = gEdgesMatrix(index1,3)+1; %update weight by adding 1
                    gEdgesMatrix(index1,4) = gEdgesMatrix(index1,4)+1; %update numBonds by adding 1
                    weightMatrix(connectedNodes,connectedNodes) = weightMatrix(connectedNodes,connectedNodes) + 1;
                    numBondsMatrix(connectedNodes,connectedNodes) = numBondsMatrix(connectedNodes,connectedNodes) + 1;
                    if gEdgesMatrix(index1,6) == 0 %singleBond1Length
                        gEdgesMatrix(index1,6) = edgeToAddLength;
                    elseif gEdgesMatrix(index1,7) == 0 %singleBond2Length
                        gEdgesMatrix(index1,7) = edgeToAddLength;
                    elseif gEdgesMatrix(index1,8) == 0 %singleBond3Length
                        gEdgesMatrix(index1,8) = edgeToAddLength;
                    elseif gEdgesMatrix(index1,9) == 0 %singleBond4Length
                        gEdgesMatrix(index1,9) = edgeToAddLength;
                    else
                        disp('something went strangely wrong')
                    end
                end

                
                allConnectedNodes(j) = []; %so now allConnectedNodes is the list of the other two singly connected nodes.
                connectedNodes = allConnectedNodes;
                edgesOfNode = [connectedNodes, [node; node]];
                %figure out how long (how many ntds) the bond we're adding is.
                %this line (below) is equivalent to: [~,index1] = ismember(edgesOfNode(1,:),gEdgesMatrix(:,1:2),'rows');
                index1 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(1,:))'));
                if isempty(index1) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                    %this line (below) is equivalent to: [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdgesMatrix(:,1:2),'rows');
                    index1 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(1,:)))'));
                end
                %this line (below) is equivalent to: [~,index2] = ismember(edgesOfNode(2,:),gEdgesMatrix(:,1:2),'rows');
                index2 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(2,:))'));
                if isempty(index2) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                    %this line (below) is equivalent to: [~,index2] = ismember(fliplr(edgesOfNode(2,:)),gEdgesMatrix(:,1:2),'rows');
                    index2 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(2,:)))'));
                end
                edgeToAddLength = gEdgesMatrix(index1,6) + gEdgesMatrix(index2,6); %singleBond1Length + singleBond1Length
                
                %make edge between connectedNodes, if an edge doesn't already exist
                %this line (below) is equivalent to: [edgeExists,edgeIndex] = ismember(connectedNodes',gEdgesMatrix(:,1:2),'rows');
                edgeIndex = find(all((gEdgesMatrix(:,1:2) == connectedNodes')'));
                edgeExists = ~isempty(edgeIndex);
                if ~edgeExists %make sure edge doesn't already exist or else program will produce an error message
                    edgeToAdd = [connectedNodes', 1, 1, 0, edgeToAddLength, 0, 0, 0];
                    
                    gEdgesMatrix = sortrows([gEdgesMatrix;edgeToAdd]);
                    weightMatrix(connectedNodes(1),connectedNodes(2)) = weightMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                    weightMatrix(connectedNodes(2),connectedNodes(1)) = weightMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                    numBondsMatrix(connectedNodes(1),connectedNodes(2)) = numBondsMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                    numBondsMatrix(connectedNodes(2),connectedNodes(1)) = numBondsMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                else %if edge does already exist, then modify it so numBonds = numBonds + 1, and weight = weight+1, and add a single bond with the proper length.
                    gEdgesMatrix(edgeIndex,3) = gEdgesMatrix(edgeIndex,3)+1; %update weight by adding 1
                    gEdgesMatrix(edgeIndex,4) = gEdgesMatrix(edgeIndex,4)+1; %update numBonds by adding 1
                    weightMatrix(connectedNodes(1),connectedNodes(2)) = weightMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                    weightMatrix(connectedNodes(2),connectedNodes(1)) = weightMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                    numBondsMatrix(connectedNodes(1),connectedNodes(2)) = numBondsMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                    numBondsMatrix(connectedNodes(2),connectedNodes(1)) = numBondsMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                    
                    if gEdgesMatrix(edgeIndex,6) == 0 %singleBond1Length
                        gEdgesMatrix(edgeIndex,6) = edgeToAddLength;
                    elseif gEdgesMatrix(edgeIndex,7) == 0 %singleBond2Length
                        gEdgesMatrix(edgeIndex,7) = edgeToAddLength;
                    elseif gEdgesMatrix(edgeIndex,8) == 0 %singleBond3Length
                        gEdgesMatrix(edgeIndex,8) = edgeToAddLength;
                    elseif gEdgesMatrix(edgeIndex,9) == 0 %singleBond4Length
                        gEdgesMatrix(edgeIndex,9) = edgeToAddLength;
                    else
                        disp('something is going on')
                    end
                end
                
                
            elseif length(allConnectedNodes) == 4
                %find indices of connectedNodes in the graph h, which doesn't
                %have the node ``node''. We'll use that to find out how the
                %graph disconnects with the removal of ``node''.
                connectedNodesH = allConnectedNodes;
                connectedNodesH(connectedNodesH>node) = connectedNodesH(connectedNodesH>node) - 1;
                conncompH = myConnComp2(hWeightMatrix);
                aa = conncompH(connectedNodesH(1));
                bb = conncompH(connectedNodesH(2));
                cc = conncompH(connectedNodesH(3));
                dd = conncompH(connectedNodesH(4));

                %rearrange order in connectedNodes (which are the real values of the nodes) so that aa == bb
                if aa == cc
                    allConnectedNodes = [allConnectedNodes(1), allConnectedNodes(3), allConnectedNodes(2), allConnectedNodes(4)]';
                elseif aa == dd
                    allConnectedNodes = [allConnectedNodes(1), allConnectedNodes(4), allConnectedNodes(2), allConnectedNodes(3)]';
                end

                for j = [1,3]
                    connectedNodes = allConnectedNodes(j:j+1);

                    edgesOfNode = [connectedNodes, [node; node]];
                    %figure out how long (how many ntds) the bond we're adding is.
                    %this line (below) is equivalent to: [~,index1] = ismember(edgesOfNode(1,:),gEdgesMatrix(:,1:2),'rows'); %EndNodes
                    index1 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(1,:))'));
                    if isempty(index1) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                        %this line (below) is equivalent to: [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdgesMatrix(:,1:2),'rows'); %EndNodes
                        index1 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(1,:)))'));
                    end
                    %this line (below) is equivalent to: [~,index2] = ismember(edgesOfNode(2,:),gEdgesMatrix(:,1:2),'rows'); %EndNodes
                    index2 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(2,:))'));
                    if isempty(index2) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                        %this line (below) is equivalent to: [~,index2] = ismember(fliplr(edgesOfNode(2,:)),gEdgesMatrix(:,1:2),'rows');
                        index2 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(2,:)))'));
                    end
                    edgeToAddLength = gEdgesMatrix(index1,6) + gEdgesMatrix(index2,6); %singleBond1Length + singleBond1Length

                    %make edge between connectedNodes, if an edge doesn't already exist
                    %this line (below) is equivalent to: [edgeExists,edgeIndex] = ismember(connectedNodes',gEdgesMatrix(:,1:2),'rows');
                    edgeIndex = find(all((gEdgesMatrix(:,1:2) == connectedNodes')'));
                    edgeExists = ~isempty(edgeIndex);
                    if ~edgeExists %make sure edge doesn't already exist or else program will produce an error message
                        edgeToAdd = [connectedNodes', 1, 1, 0, edgeToAddLength, 0, 0, 0];
                        gEdgesMatrix = sortrows([gEdgesMatrix;edgeToAdd]);
                        weightMatrix(connectedNodes(1),connectedNodes(2)) = weightMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                        weightMatrix(connectedNodes(2),connectedNodes(1)) = weightMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                        numBondsMatrix(connectedNodes(1),connectedNodes(2)) = numBondsMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                        numBondsMatrix(connectedNodes(2),connectedNodes(1)) = numBondsMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                    else %if edge does already exist, then modify it so numBonds = numBonds + 1, and weight = weight+1, and add a single bond with the proper length.
                        gEdgesMatrix(edgeIndex,3) = gEdgesMatrix(edgeIndex,3)+1; %update weight by adding 1
                        gEdgesMatrix(edgeIndex,4) = gEdgesMatrix(edgeIndex,4)+1; %update numBonds by adding 1
                        weightMatrix(connectedNodes(1),connectedNodes(2)) = weightMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                        weightMatrix(connectedNodes(2),connectedNodes(1)) = weightMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                        numBondsMatrix(connectedNodes(1),connectedNodes(2)) = numBondsMatrix(connectedNodes(1),connectedNodes(2)) + 1;
                        numBondsMatrix(connectedNodes(2),connectedNodes(1)) = numBondsMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                        
                        
                        if gEdgesMatrix(edgeIndex,6) == 0 %singleBond1Length
                            gEdgesMatrix(edgeIndex,6) = edgeToAddLength;
                        elseif gEdgesMatrix(edgeIndex,7) == 0 %singleBond2Length
                            gEdgesMatrix(edgeIndex,7) = edgeToAddLength;
                        elseif gEdgesMatrix(edgeIndex,8) == 0 %singleBond3Length
                            gEdgesMatrix(edgeIndex,8) = edgeToAddLength;
                        elseif gEdgesMatrix(edgeIndex,9) == 0 %singleBond4Length
                            gEdgesMatrix(edgeIndex,9) = edgeToAddLength;
                        else
                            disp('something is going on')
                        end
                    end
                end
            
            else
                disp('huh thats a mistake')
            end
            
            
            %find edges that include node node and delete them. 
            gEdgesMatrix(gEdgesMatrix(:,1)==node | gEdgesMatrix(:,2)==node ,:) = [];

            %After deleting these edges, we need to update the names of the
            %nodes.
            gEdgesMatrix(gEdgesMatrix(:,1)>node,1) = ...
                gEdgesMatrix(gEdgesMatrix(:,1)>node,1) - 1;
            gEdgesMatrix(gEdgesMatrix(:,2)>node,2) = ...
                gEdgesMatrix(gEdgesMatrix(:,2)>node,2) - 1;

            weightMatrix(:,node) = [];
            weightMatrix(node,:) = [];
            numBondsMatrix(:,node) = [];
            numBondsMatrix(node,:) = [];
            
            possibleNodesToDelete = find(sum(weightMatrix)==4 & sum(numBondsMatrix)==4)';
            break %don't delete any more nodes before checking if we still need to after having deleted this one.
        else 
            possibleNodesToDelete(i) = [];
            break %because the index i needs to be changed.
        end
    end
end




%We don't need to consider the other edge case -- nodes which have three 
%single bonds attached. There are two cases with these nodes -- 1) either 
%all their neighbors are connected to one another (without them) (in which 
%case we don't want to delete them); or, 2) by deleting the node, we can 
%make the graph disconnected. There are three possibilities for case 2: 
%Either the node in question is connected to two nodes, and to one of them 
%twice -- in that case, the free single edge would have already been deleted; 
%or it's connected to three nodes, two of which are connected to one another 
%-- also in this case the free single edge would have been deleted; or it's 
%connected to one node with three single bonds, in which case, we don't want to delete it.


%%

%now do the main thing which is deleting nodes whose only connections are
%two single bonds (but which aren't connected to themselves because those are
%closed nets 0).

nodesToDelete = find(sum(weightMatrix)==2 & sum(numBondsMatrix)==2)';
%note that the indicies of these will change as we delete nodes.

i = 1; %don't include nodes that are connected to themselves. 
%If we don't include this, nodes with two single bonds connecting them to themselves will be in nodesToDelete
while i <= length(nodesToDelete)
    if weightMatrix(nodesToDelete(i),nodesToDelete(i))~=0
        nodesToDelete(i) = [];
    else
        i = i+1;
    end
end

while ~isempty(nodesToDelete) %delete nodes one by one so that you don't accidentally delete all three nodes in a cycle, for example.
    node = nodesToDelete(1); %make sure not to call it nodes to distinguish it from the vector of nodes
    connectedNodes = find(weightMatrix(node,:)~=0)';
    %find all nodes within a distance 1 of the given node,
    %treating all weights as 1 (so we don't have to separately deal with the times when numBonds = Weight = 2)
    
    if isempty(connectedNodes) %then we have a node connected to itself
        %we know that the node is a closed net zero so we might as well
        %compute its entropy and delete it
        
        nodesToDelete(1)=[];
    elseif length(connectedNodes)==1 %then we have two nodes connected to each other with two single bonds
        edgesOfNode = [node, connectedNodes];
        %this line (below) is equivalent to: [~,index1] = ismember(edgesOfNode(1,:),gEdgesMatrix(:,1:2),'rows');
        index1 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(1,:))'));
        if isempty(index1) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
            %this line (below) is equivalent to: [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdgesMatrix(:,1:2),'rows');
            index1 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(1,:)))'));
        end
        edgeToAddLength = gEdgesMatrix(index1,6) + gEdgesMatrix(index1,7)... %singleBond1Length + singleBond2Length
                    + gEdgesMatrix(index1,8) + gEdgesMatrix(index1,9); %singleBond3Length + singleBond4Length
        if ~any(all((gEdgesMatrix(:,1:2) == [connectedNodes connectedNodes])')) 
            %this line (above) is equivalent to: ~ismember([connectedNodes connectedNodes], gEdgesMatrix(:,1:2),'rows')
            edgeToAdd = [connectedNodes, connectedNodes, 1, 1, 0, edgeToAddLength, 0, 0, 0];%note that connectedNodes has length 1
            gEdgesMatrix = sortrows([gEdgesMatrix;edgeToAdd]);
            weightMatrix(connectedNodes,connectedNodes) = weightMatrix(connectedNodes,connectedNodes) + 1;
            numBondsMatrix(connectedNodes,connectedNodes) = numBondsMatrix(connectedNodes,connectedNodes) + 1;
        else     
            %this line (below) is equivalent to: [~,index1] = ismember([connectedNodes connectedNodes],gEdgesMatrix(:,1:2),'rows');
            index1 = find(all((gEdgesMatrix(:,1:2) == [connectedNodes connectedNodes])'));
            gEdgesMatrix(index1,3) = gEdgesMatrix(index1,3)+1; %update weight by adding 1
            gEdgesMatrix(index1,4) = gEdgesMatrix(index1,4)+1; %update numBonds by adding 1
            weightMatrix(connectedNodes,connectedNodes) = weightMatrix(connectedNodes,connectedNodes) + 1;
            numBondsMatrix(connectedNodes,connectedNodes) = numBondsMatrix(connectedNodes,connectedNodes) + 1;
            if gEdgesMatrix(index1,6) == 0 %singleBond1Length
                gEdgesMatrix(index1,6) = edgeToAddLength;
            elseif gEdgesMatrix(index1,7) == 0 %singleBond2Length
                gEdgesMatrix(index1,7) = edgeToAddLength;
            elseif gEdgesMatrix(index1,8) == 0 %singleBond3Length
                gEdgesMatrix(index1,8) = edgeToAddLength;
            elseif gEdgesMatrix(index1,9) == 0 %singleBond4Length
                gEdgesMatrix(index1,9) = edgeToAddLength;
            else
                disp('something went strangely wrong')
            end
        end
        
        %delete edges which include node node
        gEdgesMatrix(gEdgesMatrix(:,1)==node | gEdgesMatrix(:,2)==node ,:) = [];
        %After deleting these edges, we need to update the names of the
        %nodes.
        gEdgesMatrix(gEdgesMatrix(:,1)>node,1) = gEdgesMatrix(gEdgesMatrix(:,1)>node,1) - 1;
        gEdgesMatrix(gEdgesMatrix(:,2)>node,2) = gEdgesMatrix(gEdgesMatrix(:,2)>node,2) - 1;
        weightMatrix(:,node) = [];
        weightMatrix(node,:) = [];
        numBondsMatrix(:,node) = [];
        numBondsMatrix(node,:) = [];
        
        %take connectedNodes off of list of nodesToDelete
        nodesToDelete(nodesToDelete==connectedNodes)=[];
        
        nodesToDelete = nodesToDelete-1;
        nodesToDelete(1)=[];
    elseif length(connectedNodes) == 2
        edgesOfNode = [connectedNodes, [node; node]];
        %figure out how long (how many ntds) the bond we're adding is.
        %this line (below) is equivalent to: [~,index1] = ismember(edgesOfNode(1,:),gEdgesMatrix(:,1:2),'rows');
        index1 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(1,:))'));
        if isempty(index1) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
            %this line (below) is equivalent to: [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdgesMatrix(:,1:2),'rows');
            index1 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(1,:)))'));
        end
        %this line (below) is equivalent to: [~,index2] = ismember(edgesOfNode(2,:),gEdgesMatrix(:,1:2),'rows');
        index2 = find(all((gEdgesMatrix(:,1:2) == edgesOfNode(2,:))'));
        if isempty(index2) %it's because we need to switch the order of the edges in edgesOfNode(:,1)
            %this line (below) is equivalent to: [~,index2] = ismember(fliplr(edgesOfNode(2,:)),gEdgesMatrix(:,1:2),'rows');
            index2 = find(all((gEdgesMatrix(:,1:2) == fliplr(edgesOfNode(2,:)))'));
        end
       
        edgeToAddLength = gEdgesMatrix(index1,6) + gEdgesMatrix(index2,6); %singleBond1Length + singleBond1Length
        
        %make edge between connectedNodes, if an edge doesn't
        %already exist
        %this line (below) is equivalent to: [edgeExists,edgeIndex] = ismember(connectedNodes',gEdgesMatrix(:,1:2),'rows');
        edgeIndex = find(all((gEdgesMatrix(:,1:2) == connectedNodes')'));
        edgeExists = ~isempty(edgeIndex);
        if ~edgeExists %make sure edge doesn't already exist or else program will produce an error message
            edgeToAdd = [connectedNodes', 1, 1, 0, edgeToAddLength, 0, 0, 0];
            
            gEdgesMatrix = sortrows([gEdgesMatrix;edgeToAdd]);
            weightMatrix(connectedNodes(1),connectedNodes(2)) = weightMatrix(connectedNodes(1),connectedNodes(2)) + 1;
            weightMatrix(connectedNodes(2),connectedNodes(1)) = weightMatrix(connectedNodes(2),connectedNodes(1)) + 1;
            numBondsMatrix(connectedNodes(1),connectedNodes(2)) = numBondsMatrix(connectedNodes(1),connectedNodes(2)) + 1;
            numBondsMatrix(connectedNodes(2),connectedNodes(1)) = numBondsMatrix(connectedNodes(2),connectedNodes(1)) + 1;
                        
        else %if edge does already exist, then modify it so numBonds = numBonds + 1, and weight = weight+1, and add a single bond with the proper length.
            gEdgesMatrix(edgeIndex,3) = gEdgesMatrix(edgeIndex,3)+1; %update weight by adding 1
            gEdgesMatrix(edgeIndex,4) = gEdgesMatrix(edgeIndex,4)+1; %update numBonds by adding 1
            weightMatrix(connectedNodes(1),connectedNodes(2)) = weightMatrix(connectedNodes(1),connectedNodes(2)) + 1;
            weightMatrix(connectedNodes(2),connectedNodes(1)) = weightMatrix(connectedNodes(2),connectedNodes(1)) + 1;
            numBondsMatrix(connectedNodes(1),connectedNodes(2)) = numBondsMatrix(connectedNodes(1),connectedNodes(2)) + 1;
            numBondsMatrix(connectedNodes(2),connectedNodes(1)) = numBondsMatrix(connectedNodes(2),connectedNodes(1)) + 1;
            
            if gEdgesMatrix(edgeIndex,6) == 0 %singleBond1Length
                gEdgesMatrix(edgeIndex,6) = edgeToAddLength;
            elseif gEdgesMatrix(edgeIndex,7) == 0 %singleBond2Length
                gEdgesMatrix(edgeIndex,7) = edgeToAddLength;
            elseif gEdgesMatrix(edgeIndex,8) == 0 %singleBond3Length
                gEdgesMatrix(edgeIndex,8) = edgeToAddLength;
            elseif gEdgesMatrix(edgeIndex,9) == 0 %singleBond4Length
                gEdgesMatrix(edgeIndex,9) = edgeToAddLength;
            else
                disp('something is going on')
            end
            
        end

        
        %delete edges which include node node
        gEdgesMatrix(gEdgesMatrix(:,1)==node | gEdgesMatrix(:,2)==node ,:) = [];
        %After deleting these edges, we need to update the names of the
        %nodes.
        gEdgesMatrix(gEdgesMatrix(:,1)>node,1) = gEdgesMatrix(gEdgesMatrix(:,1)>node,1) - 1;
        gEdgesMatrix(gEdgesMatrix(:,2)>node,2) = gEdgesMatrix(gEdgesMatrix(:,2)>node,2) - 1;
        
        weightMatrix(:,node) = [];
        weightMatrix(node,:) = [];
        numBondsMatrix(:,node) = [];
        numBondsMatrix(node,:) = [];
        
        nodesToDelete = nodesToDelete-1;
        nodesToDelete(1)=[];
    else
        disp('not sure why we have more than two connected nodes')
    end
end