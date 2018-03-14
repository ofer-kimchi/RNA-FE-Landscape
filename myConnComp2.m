function [nodeComponent] = myConnComp2(weightMatrix)
%return a vector of length numNodes where if elements i & j of the vector are
%equal, the nodes i & j belong to the same connected component of the graph
%created by the adjacency (weighted or unweighted) matrix given as an
%input.

a = weightMatrix~=0;
numNodes = length(a);
a(1:numNodes+1:numNodes^2) = ones(1,numNodes); %put ones on diagonal elements

[p,~,r,~] = dmperm(a);
nodeComponent = zeros(1,numNodes);
for i = 1:length(r)-1
    nodeComponent(p(r(i):r(i+1)-1)) = i;
end