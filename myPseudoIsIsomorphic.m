function [yn] = myPseudoIsIsomorphic(weightMatrix1,weightMatrix2,numBondsMatrix1,numBondsMatrix2)
%returns yn=true if the two graphs are isomorphic.

%Actually, since we're just checking whether the eigenvalues are the same,
%we haven't proved that this works for all possible graphs; we only proved 
%(by exhaustive enumeration) that it works for all graphs of up to 4 nodes 
%(which represent topologies we can deal with in our calculations). However,
%we have not found evidence of any two graphs that this function says are 
%isomorphic but aren't in reality. 

%More precisely, this function returns true if the two graphs are probably 
%isomorphic, and false if they definitely aren't isomorphic. If it returns 
%false, the graphs definitely aren't isomorphic; if it returns true, they still
%might not be isomorphic but it would be slow to check.


%the code seems to run faster if we put everything in separate if statements
%instead of using &&. 

if length(weightMatrix1)==length(weightMatrix2) %check that their sizes are the same
    if nnz(weightMatrix1) == nnz(weightMatrix2)
        if isequal(sum(sum(weightMatrix1)),sum(sum(weightMatrix2)))
            if isequal(sum(sum(numBondsMatrix1)),sum(sum(numBondsMatrix2)))
                if all(abs(eig(weightMatrix1)-eig(weightMatrix2))<1e-8)
                    if all(abs(eig(numBondsMatrix1)-eig(numBondsMatrix2))<1e-8)
                        yn = true;
                    else
                        yn = false;
                    end
                else
                    yn = false;
                end
            else
                yn = false;
            end
        else
            yn = false;
        end
    else
        yn = false;
    end
else
    yn = false;
end