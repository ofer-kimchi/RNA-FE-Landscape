function[g,gEdges] = deleteSuperfluousNodes_withGraph_Fxn(g,gEdges)
%Essentially the same as deleteSuperfluousNodes_Fxn but considers graphs
%instead of matrices

numNodes = numnodes(g);
nodeList = 1:numNodes;

%first we consider an edge case -- nodes which have four single bonds
%attached. There are two cases with these nodes -- 1) either all their
%neighbors are connected to one another (without them); or, 2) two of their
%neighbors are connected to each other, and the other two are as well. This
%is because if three of their neighbors are connected to each other and the
%other is not connected, then the bond between that neighbor and the node
%in question would have already been broken in order to make a disconnected
%graph. We are interested in case 2: in this case, we can delete the node,
%and make the graph disconnected.


possibleNodesToDelete = nodeList(centrality(g,'degree','Importance',gEdges.Weight)==4 & ...
    centrality(g,'degree','Importance',gEdges.numBonds)==4)';
while ~isempty(possibleNodesToDelete)
    for i = 1:length(possibleNodesToDelete)
        node = possibleNodesToDelete(i);
        h = rmnode(g,node);
        if max(conncomp(h)) > max(conncomp(g)) %then deleting that node disconnects the graph
            %Delete this node. If the node is connected to nodes
            %a, b, c, d, and when it's deleted the graph splits into a,b and c,d, then
            %make a connection between a and b whose length is the sum of the lengths
            %of the connections between a and node, and b and node (and same for c,d).
            allConnectedNodes = nearest(g,node,1,'Method','unweighted');
            
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
                    [~,index1] = ismember(edgesOfNode(1,:),gEdges.EndNodes,'rows');
                    if index1 ~= 0 && gEdges.singleBond2Length(index1) ~= 0 %then this is the node that's connected twice
                        break
                    end
                    [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdges.EndNodes,'rows');
                    if index1 ~= 0 && gEdges.singleBond2Length(index1) ~= 0 %then this is the node that's connected twice
                        break
                    end
                    if j == 3 %if we've reached this point at j ==3 then our assumptions are wrong
                        disp('our assumptions are wrong')
                    end
                end
                
                connectedNodes = allConnectedNodes(j); %this is the node that's connected twice
                edgeToAddLength = gEdges.singleBond1Length(index1) + gEdges.singleBond2Length(index1)...
                    + gEdges.singleBond3Length(index1) + gEdges.singleBond4Length(index1);
                if ~ismember([connectedNodes connectedNodes], gEdges.EndNodes,'rows')
                    edgeToAdd = table([connectedNodes connectedNodes],'VariableNames',{'EndNodes'});
                    edgeToAdd.Weight = 1;
                    edgeToAdd.numBonds = 1;
                    edgeToAdd.doubleBondLength = 0;
                    edgeToAdd.singleBond1Length = edgeToAddLength;
                    edgeToAdd.singleBond2Length = 0;
                    edgeToAdd.singleBond3Length = 0;
                    edgeToAdd.singleBond4Length = 0;
                    g = addedge(g,edgeToAdd);
                    gEdges = sortrows([gEdges;edgeToAdd]);
                else
                    [~,index1] = ismember([connectedNodes connectedNodes],gEdges.EndNodes,'rows');
                    g.Edges.Weight(index1) = g.Edges.Weight(index1)+1;
                    gEdges.Weight(index1) = gEdges.Weight(index1)+1; %first change the graph, then the table which is separate from the graph
                    g.Edges.numBonds(index1) = g.Edges.numBonds(index1)+1;
                    gEdges.numBonds(index1) = gEdges.numBonds(index1)+1;
                    if gEdges.singleBond1Length(index1) == 0
                        g.Edges.singleBond1Length(index1) = edgeToAddLength;
                        gEdges.singleBond1Length(index1) = edgeToAddLength;
                    elseif gEdges.singleBond2Length(index1) == 0
                        g.Edges.singleBond2Length(index1) = edgeToAddLength;
                        gEdges.singleBond2Length(index1) = edgeToAddLength;
                    elseif gEdges.singleBond3Length(index1) == 0
                        g.Edges.singleBond3Length(index1) = edgeToAddLength;
                        gEdges.singleBond3Length(index1) = edgeToAddLength;
                    elseif gEdges.singleBond4Length(index1) == 0
                        g.Edges.singleBond4Length(index1) = edgeToAddLength;
                        gEdges.singleBond4Length(index1) = edgeToAddLength;
                    else
                        disp('something went strangely wrong')
                    end
                end

                
                allConnectedNodes(j) = []; %so now allConnectedNodes is the list of the other two singly connected nodes.
                connectedNodes = allConnectedNodes;
                edgesOfNode = [connectedNodes, [node; node]];
                %figure out how long (how many ntds) the bond we're adding is.
                [~,index1] = ismember(edgesOfNode(1,:),gEdges.EndNodes,'rows');
                if index1 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                    [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdges.EndNodes,'rows');
                end
                [~,index2] = ismember(edgesOfNode(2,:),gEdges.EndNodes,'rows');
                if index2 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                    [~,index2] = ismember(fliplr(edgesOfNode(2,:)),gEdges.EndNodes,'rows');
                end
                edgeToAddLength = gEdges.singleBond1Length(index1) + gEdges.singleBond1Length(index2);
                
                %make edge between connectedNodes, if an edge doesn't
                %already exist
                [edgeExists,edgeIndex] = ismember(connectedNodes',gEdges.EndNodes,'rows');
                if ~edgeExists %make sure edge doesn't already exist or else program will produce an error message
                    edgeToAdd = table(connectedNodes','VariableNames',{'EndNodes'});
                    edgeToAdd.Weight = 1;
                    edgeToAdd.numBonds = 1;
                    edgeToAdd.doubleBondLength = 0;
                    edgeToAdd.singleBond1Length = edgeToAddLength;
                    edgeToAdd.singleBond2Length = 0;
                    edgeToAdd.singleBond3Length = 0;
                    edgeToAdd.singleBond4Length = 0;
                    
                    g = addedge(g,edgeToAdd);
                    gEdges = sortrows([gEdges;edgeToAdd]);
                else %if edge does already exist, then modify it so numBonds = numBonds + 1, and weight = weight+1, and add a single bond with the proper length.
                    g.Edges.numBonds(edgeIndex)=g.Edges.numBonds(edgeIndex)+1;
                    gEdges.numBonds(edgeIndex)=gEdges.numBonds(edgeIndex)+1;
                    g.Edges.Weight(edgeIndex)=g.Edges.Weight(edgeIndex)+1;
                    gEdges.Weight(edgeIndex)=gEdges.Weight(edgeIndex)+1;
                    if gEdges.singleBond1Length(edgeIndex) == 0
                        g.Edges.singleBond1Length(edgeIndex) = edgeToAddLength;
                        gEdges.singleBond1Length(edgeIndex) = edgeToAddLength;
                    elseif gEdges.singleBond2Length(edgeIndex) == 0
                        g.Edges.singleBond2Length(edgeIndex) = edgeToAddLength;
                        gEdges.singleBond2Length(edgeIndex) = edgeToAddLength;
                    elseif gEdges.singleBond3Length(edgeIndex) == 0
                        g.Edges.singleBond3Length(edgeIndex) = edgeToAddLength;
                        gEdges.singleBond3Length(edgeIndex) = edgeToAddLength;
                    elseif gEdges.singleBond4Length(edgeIndex) == 0
                        g.Edges.singleBond4Length(edgeIndex) = edgeToAddLength;
                        gEdges.singleBond4Length(edgeIndex) = edgeToAddLength;
                    else
                        disp('something is going on');
                    end
                end
                
                
            elseif length(allConnectedNodes) == 4
                %find indices of connectedNodes in the graph h, which doesn't
                %have the node ``node''. We'll use that to find out how the
                %graph disconnects with the removal of ``node''.
                connectedNodesH = allConnectedNodes;
                connectedNodesH(connectedNodesH>node) = connectedNodesH(connectedNodesH>node) - 1;
                conncompH = conncomp(h);
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
                    [~,index1] = ismember(edgesOfNode(1,:),gEdges.EndNodes,'rows');
                    if index1 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                        [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdges.EndNodes,'rows');
                    end
                    [~,index2] = ismember(edgesOfNode(2,:),gEdges.EndNodes,'rows');
                    if index2 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
                        [~,index2] = ismember(fliplr(edgesOfNode(2,:)),gEdges.EndNodes,'rows');
                    end
                    edgeToAddLength = gEdges.singleBond1Length(index1) + gEdges.singleBond1Length(index2);

                    %make edge between connectedNodes, if an edge doesn't
                    %already exist
                    [edgeExists,edgeIndex] = ismember(connectedNodes',gEdges.EndNodes,'rows');
                    if ~edgeExists %make sure edge doesn't already exist or else program will produce an error message
                        edgeToAdd = table(connectedNodes','VariableNames',{'EndNodes'});
                        edgeToAdd.Weight = 1;
                        edgeToAdd.numBonds = 1;
                        edgeToAdd.doubleBondLength = 0;
                        edgeToAdd.singleBond1Length = edgeToAddLength;
                        edgeToAdd.singleBond2Length = 0;
                        edgeToAdd.singleBond3Length = 0;
                        edgeToAdd.singleBond4Length = 0;

                        g = addedge(g,edgeToAdd);
                        gEdges = sortrows([gEdges;edgeToAdd]);
                    else %if edge does already exist, then modify it so numBonds = numBonds + 1, and weight = weight+1, and add a single bond with the proper length.
                        g.Edges.numBonds(edgeIndex)=g.Edges.numBonds(edgeIndex)+1;
                        gEdges.numBonds(edgeIndex)=gEdges.numBonds(edgeIndex)+1;
                        g.Edges.Weight(edgeIndex)=g.Edges.Weight(edgeIndex)+1;
                        gEdges.Weight(edgeIndex)=gEdges.Weight(edgeIndex)+1;
                        if gEdges.singleBond1Length(edgeIndex) == 0
                            g.Edges.singleBond1Length(edgeIndex) = edgeToAddLength;
                            gEdges.singleBond1Length(edgeIndex) = edgeToAddLength;
                        elseif gEdges.singleBond2Length(edgeIndex) == 0
                            g.Edges.singleBond2Length(edgeIndex) = edgeToAddLength;
                            gEdges.singleBond2Length(edgeIndex) = edgeToAddLength;
                        elseif gEdges.singleBond3Length(edgeIndex) == 0
                            g.Edges.singleBond3Length(edgeIndex) = edgeToAddLength;
                            gEdges.singleBond3Length(edgeIndex) = edgeToAddLength;
                        elseif gEdges.singleBond4Length(edgeIndex) == 0
                            g.Edges.singleBond4Length(edgeIndex) = edgeToAddLength;
                            gEdges.singleBond4Length(edgeIndex) = edgeToAddLength;
                        else
                            disp('something is going on');
                        end
                    end
                end
            
            else
                disp('huh thats a mistake')
            end
            
            g = rmnode(g,node); 
            gEdges = g.Edges; %after removing a node, all the edges containing that node get deleted, and the nodes get reordered.
            
            possibleNodesToDelete = nodeList(centrality(g,'degree','Importance',gEdges.Weight)==4 & ...
                centrality(g,'degree','Importance',gEdges.numBonds)==4)';
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



%now do the main thing which is deleting nodes whose only connections are
%two single bonds (but which aren't connected to themselves because those are
%closed nets 0).

nodesToDelete = nodeList(centrality(g,'degree','Importance',gEdges.Weight)==2 & ...
    centrality(g,'degree','Importance',gEdges.numBonds)==2)'; %note that the indicies of these will change as we delete nodes.
        
while ~isempty(nodesToDelete) %delete nodes one by one so that you don't accidentally delete all three nodes in a cycle, for example.
    node = nodesToDelete(1); %make sure not to call it nodes to distinguish it from the vector of nodes
    connectedNodes = nearest(g,node,1,'Method','unweighted'); %find all nodes within a distance 1 of the given node,
    %treating all weights as 1 (so we don't have to separately deal with the times when numBonds = Weight = 2)
    
    if isempty(connectedNodes) %then we have a node connected to itself
        %we know that the node is a closed net zero so we might as well
        %compute its entropy and delete it
        
        %uncomment the following six lines (with % at their end) if we are
        %having isisomorphic check for closed net 0 as well
        %             edgeIndex = findedge(g,node,node);% %find the index of the edge connecting node to itself
        %             s0 = g.Edges.singleBond1Length(edgeIndex);%
        %             loopEntropy = loopEntropy + kB*(log(alpha)-(3/2)*log(s0/beta));%
        %
        %             g = rmnode(g,node);%
        %             numNodes = numNodes-1;%
        %             nodesToDelete = nodesToDelete-1;%
        
        nodesToDelete(1)=[];
    elseif length(connectedNodes)==1 %then we have two nodes connected to each other with two single bonds
        edgesOfNode = [node, connectedNodes];
        [~,index1] = ismember(edgesOfNode(1,:),gEdges.EndNodes,'rows');
        if index1 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
            [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdges.EndNodes,'rows');
        end
        edgeToAddLength = gEdges.singleBond1Length(index1) + gEdges.singleBond2Length(index1)...
            + gEdges.singleBond3Length(index1) + gEdges.singleBond4Length(index1);
        if ~ismember([connectedNodes connectedNodes], gEdges.EndNodes,'rows')
            edgeToAdd = table([connectedNodes connectedNodes],'VariableNames',{'EndNodes'});
            edgeToAdd.Weight = 1;
            edgeToAdd.numBonds = 1;
            edgeToAdd.doubleBondLength = 0;
            edgeToAdd.singleBond1Length = edgeToAddLength;
            edgeToAdd.singleBond2Length = 0;
            edgeToAdd.singleBond3Length = 0;
            edgeToAdd.singleBond4Length = 0;
            g = addedge(g,edgeToAdd);
            gEdges = sortrows([gEdges;edgeToAdd]);
        else
            [~,index1] = ismember([connectedNodes connectedNodes],gEdges.EndNodes,'rows');
            g.Edges.Weight(index1) = g.Edges.Weight(index1)+1;
            gEdges.Weight(index1) = gEdges.Weight(index1)+1;
            g.Edges.numBonds(index1) = g.Edges.numBonds(index1)+1;
            gEdges.numBonds(index1) = gEdges.numBonds(index1)+1;
            if gEdges.singleBond1Length(index1) == 0
                g.Edges.singleBond1Length(index1) = edgeToAddLength;
                gEdges.singleBond1Length(index1) = edgeToAddLength;
            elseif gEdges.singleBond2Length(index1) == 0
                g.Edges.singleBond2Length(index1) = edgeToAddLength;
                gEdges.singleBond2Length(index1) = edgeToAddLength;
            elseif gEdges.singleBond3Length(index1) == 0
                g.Edges.singleBond3Length(index1) = edgeToAddLength;
                gEdges.singleBond3Length(index1) = edgeToAddLength;
            elseif gEdges.singleBond4Length(index1) == 0
                g.Edges.singleBond4Length(index1) = edgeToAddLength;
                gEdges.singleBond4Length(index1) = edgeToAddLength;
            else
                disp('something went strangely wrong')
            end
        end
        g = rmnode(g,node);
        gEdges = g.Edges; %after removing a node, all the edges containing that node get deleted, and the nodes get reordered.
        
        %take connectedNodes off of list of nodesToDelete
        nodesToDelete(nodesToDelete==connectedNodes)=[];
        
        nodesToDelete = nodesToDelete-1;
        nodesToDelete(1)=[];
    elseif length(connectedNodes) == 2
        edgesOfNode = [connectedNodes, [node; node]];
        %figure out how long (how many ntds) the bond we're adding is.
        [~,index1] = ismember(edgesOfNode(1,:),gEdges.EndNodes,'rows');
        if index1 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
            [~,index1] = ismember(fliplr(edgesOfNode(1,:)),gEdges.EndNodes,'rows');
        end
        [~,index2] = ismember(edgesOfNode(2,:),gEdges.EndNodes,'rows');
        if index2 == 0 %it's because we need to switch the order of the edges in edgesOfNode(:,1)
            [~,index2] = ismember(fliplr(edgesOfNode(2,:)),gEdges.EndNodes,'rows');
        end
        edgeToAddLength = gEdges.singleBond1Length(index1) + gEdges.singleBond1Length(index2);
        
        %make edge between connectedNodes, if an edge doesn't
        %already exist
        [edgeExists,edgeIndex] = ismember(connectedNodes',gEdges.EndNodes,'rows');
        if ~edgeExists %make sure edge doesn't already exist or else program will produce an error message
            edgeToAdd = table(connectedNodes','VariableNames',{'EndNodes'});
            edgeToAdd.Weight = 1;
            edgeToAdd.numBonds = 1;
            edgeToAdd.doubleBondLength = 0;
            edgeToAdd.singleBond1Length = edgeToAddLength;
            edgeToAdd.singleBond2Length = 0;
            edgeToAdd.singleBond3Length = 0;
            edgeToAdd.singleBond4Length = 0;
            
            g = addedge(g,edgeToAdd);
            gEdges = sortrows([gEdges;edgeToAdd]);
        else %if edge does already exist, then modify it so numBonds = numBonds + 1, and weight = weight+1, and add a single bond with the proper length.
            g.Edges.numBonds(edgeIndex)=g.Edges.numBonds(edgeIndex)+1;
            gEdges.numBonds(edgeIndex)=gEdges.numBonds(edgeIndex)+1;
            g.Edges.Weight(edgeIndex)=g.Edges.Weight(edgeIndex)+1;
            gEdges.Weight(edgeIndex)=gEdges.Weight(edgeIndex)+1;
            if gEdges.singleBond1Length(edgeIndex) == 0
                g.Edges.singleBond1Length(edgeIndex) = edgeToAddLength;
                gEdges.singleBond1Length(edgeIndex) = edgeToAddLength;
            elseif gEdges.singleBond2Length(edgeIndex) == 0
                g.Edges.singleBond2Length(edgeIndex) = edgeToAddLength;
                gEdges.singleBond2Length(edgeIndex) = edgeToAddLength;
            elseif gEdges.singleBond3Length(edgeIndex) == 0
                g.Edges.singleBond3Length(edgeIndex) = edgeToAddLength;
                gEdges.singleBond3Length(edgeIndex) = edgeToAddLength;
            elseif gEdges.singleBond4Length(edgeIndex) == 0
                g.Edges.singleBond4Length(edgeIndex) = edgeToAddLength;
                gEdges.singleBond4Length(edgeIndex) = edgeToAddLength;
            else
                disp('something is going on');
            end
        end
        g = rmnode(g,node);
        gEdges = g.Edges; %after removing a node, all the edges containing that node get deleted, and the nodes get reordered.
        nodesToDelete = nodesToDelete-1;
        nodesToDelete(1)=[];
    else
        disp('not sure why we have more than two connected nodes')
    end
end