function [p] = myPseudoIsomorphism(permList,weightMatrix1,weightMatrix2,numBondsMatrix1,numBondsMatrix2)
%returns p such that if the two graphs are isomorphic, you can set: 
%wm1 = weightMatrix1(:,p); wm1 = wm1(p,:); and 
%nbm1 = numBondsMatrix1(:,p); nbm1 = nbm1(p,:);
%and then wm1 is weightMatrix2 and nbm1 is numBondsMatrix2.
%if p(i)=j, that means that the node which in weightMatrix1 is called j is
%called i in weightMatrix2. After reordering the nodes (see wm1 above), the
%node that was called j will be called i.

%permList is equivalent to allPerms{length(weightMatrix1)} in previous
%functions. We previously defined allPerms to be
% allPerms = cell(1,4); 
% for i = 1:4
%     allPerms{i} = perms(1:i);
% end

%first, check if the graphs are isomorphic
if nnz(weightMatrix1) == nnz(weightMatrix2) &&...
        isequal(sum(sum(weightMatrix1)),sum(sum(weightMatrix2))) &&...
        isequal(sum(sum(numBondsMatrix1)),sum(sum(numBondsMatrix2))) &&...
        all(abs(eig(weightMatrix1)-eig(weightMatrix2))<1e-8) &&...
        all(abs(eig(numBondsMatrix1)-eig(numBondsMatrix2))<1e-8) 
    if isequal(weightMatrix1,weightMatrix2) && isequal(numBondsMatrix1,numBondsMatrix2) %if the graphs are equivalent, we're done
        p = 1:length(weightMatrix1);
    else %if you know the two are isomorphic just go through all possible permutations till you find the right one.
        for i = 1:size(permList,1) 
            p = permList(i,:);
            weightMatrix1Perm = weightMatrix1(p,:);
            weightMatrix1Perm = weightMatrix1Perm(:,p);
            numBondsMatrix1Perm = numBondsMatrix1(:,p);
            numBondsMatrix1Perm = numBondsMatrix1Perm(p,:);
            if isequal(weightMatrix1Perm,weightMatrix2) && isequal(numBondsMatrix1Perm,numBondsMatrix2)
                break %we've found p
            end
        end
        if isequal(p,1:length(weightMatrix1)) %if we reach this point, that means we erroneously claimed the graphs were isomorphic
            %and that our function therefore needs retooling.
            disp('theres a problem in myPsuedoIsomorphism')
        end
    end
else
    p = [];
end