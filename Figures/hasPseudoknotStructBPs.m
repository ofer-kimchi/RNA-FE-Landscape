function tf = hasPseudoknotStructBPs(strucBPs)

%returns true if the structure inputted contains a pseudoknot and false if
%it contains only hairpins.

strucBPs = sortrows(strucBPs); 
%I think sortrows is unnecessary because it should already be sorted but why not make sure

%From Mauri 2005 (substitute i,j,k,l for a,b,c,d):
%If ntd a is paired with ntd b>a, and c>a is paired with d>c,
%then either a<b<c<d or a<c<d<b, but a<c<b<d is a pseudoknot.
if size(strucBPs,1) <= 1
    tf = false;
else
    tf = false;
    for i = 1:size(strucBPs,1)-1
        for j = i+1:size(strucBPs,1)
            a = strucBPs(i,1);
            b = strucBPs(i,2);
            c = strucBPs(j,1);
            d = strucBPs(j,2);

            if b<a || d < c
                disp('hasPsuedoknot huh')
            end
            if c < a
                disp('hasPseudoknot double huh')
            end

            if c < b && b < d %then it contains a pseudoknot
                tf = true;
            end
        end
    end
end