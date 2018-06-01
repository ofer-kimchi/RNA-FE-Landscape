function strucBPs = dotBracket2structBPs(dotBracket)
strucBPs = [];
dotBracket(dotBracket == ':') = '.';
allDots = repmat('.',1,length(dotBracket));
while ~strcmp(dotBracket,allDots)
    for i = length(dotBracket)-1:-1:1
        if dotBracket(i) == '(' && dotBracket(i+1) ~= '('
            foundClose = false;
            j=i;
            while ~foundClose
                j = j+1;
                if dotBracket(j) == ')'
                    foundClose = true;
                end
            end
            strucBPs = [strucBPs;[i,j]];
            dotBracket(i) = '.';
            dotBracket(j) = '.';
            break
        end
        if dotBracket(i) == '[' && dotBracket(i+1) ~= '['
            foundClose = false;
            j=i;
            while ~foundClose
                j = j+1;
                if dotBracket(j) == ']'
                    foundClose = true;
                end
            end
            strucBPs = [strucBPs;[i,j]];
            dotBracket(i) = '.';
            dotBracket(j) = '.';
        end
        if dotBracket(i) == '{' && dotBracket(i+1) ~= '{'
            foundClose = false;
            j=i;
            while ~foundClose
                j = j+1;
                if dotBracket(j) == '}'
                    foundClose = true;
                end
            end
            strucBPs = [strucBPs;[i,j]];
            dotBracket(i) = '.';
            dotBracket(j) = '.';
        end
    end
end
strucBPs = sortrows(strucBPs);
