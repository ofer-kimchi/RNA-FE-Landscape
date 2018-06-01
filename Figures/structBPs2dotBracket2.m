function dotBracket = structBPs2dotBracket2(numBPs,structure)
% tested on [1 13;2 12; 3 11; 4 16; 5 15; 6 14; 7 20; 8 19; 9 18];
dotBracketCell = cell(1,numBPs);
for i = 1:length(structure)
    for j = 1:length(structure{i})/2
        dotBracketCell{structure{i}(j)} = [i,structure{i}(j+length(structure{i})/2)];
    end
end

dotBracket = repmat('.',1,numBPs);
indexClosedParen = []; %index of closed parenthesis
indexClosedBracket = [];
indexClosedBrace = [];
indexOpenParen = []; %index of first open parenthesis
indexOpenBracket = [];
indexOpenBrace = [];
for i = 1:numBPs
    if ~isempty(dotBracketCell{i}) %meaning that ntd i is bonded
        j = dotBracketCell{i}(2); %the ntd that i is bonded to
        ijCanBeParen = true;
        for k = 1:length(indexOpenParen)
            if (i>indexOpenParen(k) && i<indexClosedParen(k)) && ~(j>indexOpenParen(k) && j<indexClosedParen(k)) ||...
                    ~(i>indexOpenParen(k) && i<indexClosedParen(k)) && (j>indexOpenParen(k) && j<indexClosedParen(k))
                ijCanBeParen = false;
            end
        end
        if ijCanBeParen
            dotBracket(i) = '(';
            dotBracket(j) = ')';
            indexOpenParen = [indexOpenParen,i];
            indexClosedParen = [indexClosedParen,j];
        else
            ijCanBeBracket = true;
            for k = 1:length(indexOpenBracket)
                if (i>indexOpenBracket(k) && i<indexClosedBracket(k)) && ~(j>indexOpenBracket(k) && j<indexClosedBracket(k)) ||...
                        ~(i>indexOpenBracket(k) && i<indexClosedBracket(k)) && (j>indexOpenBracket(k) && j<indexClosedBracket(k))
                    ijCanBeBracket = false;
                end
            end
            if ijCanBeBracket
                dotBracket(i) = '[';
                dotBracket(j) = ']';
                indexOpenBracket = [indexOpenBracket,i];
                indexClosedBracket = [indexClosedBracket,j];
            else
                ijCanBeBrace = true;
                for k = 1:length(indexOpenBrace)
                    if (i>indexOpenBrace(k) && i<indexClosedBrace(k)) && ~(j>indexOpenBrace(k) && j<indexClosedBrace(k))||...
                            ~(i>indexOpenBrace(k) && i<indexClosedBrace(k)) && (j>indexOpenBrace(k) && j<indexClosedBrace(k))
                        ijCanBeBrace = false;
                    end
                end
                if ijCanBeBrace
                    dotBracket(i) = '{';
                    dotBracket(j) = '}';
                    indexOpenBrace = [indexOpenBrace,i];
                    indexClosedBrace = [indexClosedBrace,j];
                else
                    disp('problem with structBPs2dotBracket2')
                end
            end
        end
    end
end

