function [ha, H] = rnaplotEdited(s, hf,ha,varargin)
%RNAPLOT draws RNA secondary structure
%
%   RNAPLOT(RNA2NDSTRUCT) draws the secondary structure specified in RNA2NDSTRUCT
%   in circle format. The secondary structure RNA2NDSTRUCT can be either in
%   bracket notation or in matrix notation.
%
%   RNAPLOT(RNA2NDSTRUCT, 'SEQUENCE', SEQ) uses SEQ to label each residue
%   position in the plot. SEQ can be a string of characters or a structure
%   with a field named Sequence.
%
%   RNAPLOT(RNA2NDSTRUCT, 'FORMAT', PLOTFORMAT) draws the secondary structure in a
%   representation specified by PLOTFORMAT. The valid formats are:
%
%     'circle'     - each base is represented by a dot on the circumference
%     (default)      of a circle of arbitrary size. Lines connect bases that
%                    pair with each other.
%     'diagram'    - two dimensional representation of RNA secondary
%                    structure. Each base is identified by a label
%                    corresponding to the residue position and type.
%                    Backbone and hydrogen bonds between pairing bases are
%                    represented by lines. SEQUENCE information is
%                    mandatory for this type of plot.
%     'dotdiagram' - two dimensional representation of RNA secondary
%                    structure. Each base is represented by a dot. Backbone
%                    and hydrogen bonds between pairing bases are
%                    represented by lines.
%     'graph'      - bases are displayed in their sequence position order
%                    on the abscissa of a graph. Semi-elliptical lines
%                    connect bases that pair with each other.
%     'mountain'   - two-dimensional plot where the base position is in the
%                    abscissa and the number of base pairs enclosing a given
%                    base is in the ordinate.
%     'tree'       - rooted tree where unpaired bases are represented by
%                    leaf nodes and base pairs are represented by internal
%                    nodes. The tree root is a fictitious node, not
%                    associated with any base in the secondary structure.
%
%    RNAPLOT(RNA2NDSTRUCT, ... , 'SELECTION', SEL) highlights a subset of
%    residues specified by SEL. SEL can be a numeric array representing
%    valid positions within the structure or any of the following choices:
%
%      'paired'     - include all residues that are paired.
%      'unpaired'   - include all residues that are unpaired.
%      'AU' or 'UA' - include all AU or UA base pairs.
%      'GC' or 'CG' - include all GC or CG base pairs.
%      'GU' or 'UG' - include all GU or UG base pairs.
%
%    If SEL is 'AU', 'UA', 'GC', 'CG', 'GU', or 'UG' the RNA sequence must be
%    specified through the option SEQUENCE.
%
%   RNAPLOT(RNA2NDSTRUCT, ... , 'COLORBY', COLORSCHEME) applies the specified color
%   scheme to the plot. Valid options are:
%
%      'state'   - color residues according to their pair state.
%      'residue' - color residues according to their type.
%      'pair'    - color residues according to their base pair type.
%
%   If COLORSCHEME is 'residue' or 'pair', the RNA sequence must be
%   provided with the SEQUENCE property.
%   Because internal nodes of a tree correspond to two residues, the
%   'residue' option is not allowed in the 'tree' format.
%
%   ha = RNAPLOT(...) returns the handle of the figure axis.
%
%   [ha, H] = RNAPLOT(...) returns a structure with the following fields,
%   based on the options chosen for SELECTION and COLORBY:
%
%     'Paired'   - handles to all paired residues.
%     'Unpaired' - handles to all unpaired residues.
%     'A'        - handles to all A residues.
%     'C'        - handles to all C residues.
%     'G'        - handles to all G residues.
%     'U'        - handles to all U residues.
%     'AU'       - handles to all AU/UA base pairs.
%     'GC'       - handles to all GC/CG base pairs.
%     'GU'       - handles to all GU/UG base pairs.
%     'Selected' - handles to all selected residues.
%
%   Examples:
%
%     % Predict the secondary structure and plot it using the circle representation
%     seq = 'GCGCCCGUAGCUCAAUUGGAUAGAGCGUUUGACUACGGAUCAAAAGGUUAGGGGUUCGACUCCUCUCGGGCGCG';
%     ss  = rnafold(seq);
%     rnaplot(ss);
%
%     % Plot the structure as graph and color residues by pair type
%     rnaplot(ss, 'sequence', seq, 'format', 'graph', 'colorby', 'pair');
%
%     % Plot the structure as mountain and color residues by type
%     ha = rnaplot(ss, 'sequence', seq, 'format', 'mountain', 'colorby', 'residue');
%     title(ha, 'Bacillus halodurans, tRNA Arg');
%
%     % Mutate the first six positions in the sequence and observe the
%     % effect the change has on the secondary structure
%     seqMut = seq;
%     seqMut(1:6) = 'AAAAAA';
%     ssMut = rnafold(seqMut);
%     rnaplot(ss, 'sequence', seq, 'format', 'diagram', 'selection', 1:6);
%     rnaplot(ssMut, 'sequence', seq, 'format', 'diagram', 'selection', 1:6);
%
%   See also RNACONVERT, RNADEMO, RNAFOLD.

%   Copyright 2007-2012 The MathWorks, Inc.


bioinfochecknargin(nargin,1,mfilename);

%=== default values
selmode = 0;     % state (paired vs. unpaired)
selected = [];
seq = '';
typePlot = 1;    % circle
colorScheme = 1; % state (paired vs. unpaired)
colors.green = [0 0.7 0];
colors.blue = [0 0 1];
colors.red = [1 0 0];
colors.cyan = [0 0.9 1];
colors.yellow = [1 0.8 0];

%=== input parameter parsing and error checking
nvarargin = numel(varargin);
if nvarargin
    if rem(nvarargin,2) == 1
        error(message('bioinfo:rnaplot:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'format', 'selection', 'sequence', 'colorby'};
    oksels = {'AU', 'UA', 'GU', 'UG', 'GC', 'CG', 'paired', 'unpaired'};
    oktypes = {'circle', 'diagram', 'dotdiagram', 'graph', 'mountain', 'tree'};
    okcolors = {'state', 'residue', 'pair'};

    for j=1:2:nvarargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error(message('bioinfo:rnaplot:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:rnaplot:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1 % format
                    typePlot = find(strncmpi(pval, oktypes, numel(pval)));
                    if isempty(typePlot)
                        error(message('bioinfo:rnaplot:InvalidOutputFormat'));
                    elseif length(typePlot)>1
                        error(message('bioinfo:rnaplot:AmbiguousPlotTypeName', pval));
                    end
                case 2 % selection
                    if ischar(pval)
                        selmode = find(strcmpi(pval, oksels));
                        if isempty(selmode)
                            error(message('bioinfo:rnaplot:InvalidSelectionMode'));
                        end
                    elseif isnumeric(pval) && isvector(pval) && max(pval) <= length(s) && min(pval) > 0
                        selmode = 0;
                        if size(pval,1) > 1
                            selected = pval; % make sure the numeric selection is  n x 1
                        else
                            selected = pval';
                        end
                    else
                        error(message('bioinfo:rnaplot:InvalidSelection'))
                    end
                case 3 % sequence
                    if isstruct(pval)
                        if any(size(pval) > 1)
                            error(message('bioinfo:rnaplot:InvalidLabelsStruct'));
                        else
                            pval = upper(bioinfoprivate.seqfromstruct(pval));
                        end
                    end
                    if ischar(pval) && isvector(pval)
                        % make sure it is a 1 x n char array
                        if size(pval,1) > 1
                            pval = pval';
                        end
                        % make sure the labels are right for the structure
                        if (length(pval) == length(s))
                            seq = upper(pval);
                        else
                            error(message('bioinfo:rnaplot:InvalidLabelsLength'))
                        end
                    else
                        error(message('bioinfo:rnaplot:InvalidLabels'))
                    end
                case 4 % colorby
                    if ischar(pval)
                        colorScheme = find(strncmpi(pval, okcolors, numel(pval)));
                        if isempty(colorScheme)
                            error(message('bioinfo:rnaplot:InvalidColorSchemeChoice'));
                        end
                    else
                        error(message('bioinfo:rnaplot:InvalidColorScheme'))

                    end
            end
        end
    end
end

%=== make sure sequence info is provided when needed
% seq info is required by: plot = diagram; selection = AU, GC, GU; colorScheme = residue or pair
if isempty(seq) && (typePlot == 2 || ((selmode < 7) && (selmode > 0)) || colorScheme > 1)
    error(message('bioinfo:rnaplot:MissingLabels'));
end

%=== color by residue is not allowed with tree plot
if typePlot == 6 && colorScheme == 2
    error(message('bioinfo:rnaplot:InvalidcolorSchemeCombination'))
end

%==========================================================================
% MAIN
%==========================================================================

%=== convert secondary structure into needed formats (s1 = vector, s2 = matrix)
if size(s,1) > size(s,2)
    s = s'; % if bracket notation, must be a vector 1 x n
end

try
    s2 = rnaconvert(s);
catch theErr 
    error(message('bioinfo:rnaplot:InvalidInput'));
end

if isvector(s)
    s1 = s;
else
    s1 = s2;
    s2 = s;
end

%=== define selected set
[b1, b2, basepaired, baseunpaired] = localGroupByStatus(s2);
if selmode > 0
    if selmode < 7
        pairs = localGroupByPair(seq, s2);
    end
    switch selmode
        case {1,2} % AU or UA
            sel = pairs{1,1};
            selected = unique([sel(:,1); sel(:,2)]);
        case {3,4} % GU or UG
            sel = pairs{1,2};
            selected = unique([sel(:,1); sel(:,2)]);
        case {5,6} % GC or CG
            sel = pairs{1,3};
            selected = unique([sel(:,1); sel(:,2)]);
        case 7 % Paired
            selected = basepaired;
        case 8 % Unpaired
            selected = baseunpaired;
    end

end

%=== plot
switch typePlot
    case 1 % circle
        [ha, H] = localMatrix2circle(s2, selected, seq, colorScheme, colors);
    case 2 % diagram
        [ha, H] = localMatrix2diagram(s2, selected, seq, colorScheme, colors, 0,hf,ha);
    case 3 % dotdiagram
        [ha, H] = localMatrix2diagram(s2, selected, seq, colorScheme, colors, 1,hf,ha);
    case 4 % graph
        [ha H] = localMatrix2graph(s2, selected, seq, colorScheme, colors);
    case 5 % mountain
        [ha, H] = localBracket2mountain(s1, s2, selected, seq, colorScheme, colors);
    case 6 % tree
        [ha, H] = localBracket2tree(s1, s2, selected, seq, colorScheme, colors);
end

if nargout == 0
    clear ha
end


%==========================================================================
% SUBFUNCTIONS
%==========================================================================

function [ha H] = localMatrix2circle(matrix, selected, seq, colorScheme, colors)
% Draw circle plot for the input structure in matrix notation. Residues
% are placed around a full circle and edges are drawn between residues that
% pair. Selected positions are highlighted.

N = size(matrix,1);
t = 2 * pi / N;
theta = (1:N) * t;

labelRadius = 1.07; % radius where numbering labels are placed

%=== determine numbering labels ticks
step = floor(N/100) * 5;
if step == 0
    step = 5;
end
ticks = 0:step:N;
ticks(1) = 1;
if N - ticks(end) < step/4
    ticks(end) = [];
end
tickLabels = cellfun(@num2str, num2cell(ticks), 'UniformOutput', false);

%=== prepare figure
hf = figure('Color','w');
ha = axes('Tag','Bioinfo:rnaplot:circle');
hold on; axis square; axis off;

%=== plot lines
[i, j, bp, bu] = localGroupByStatus(matrix);
hlines = plot([cos(theta(i)); cos(theta(j))], [sin(theta(i)); sin(theta(j))], ...
    '-', 'Color', 'black', 'clipping', 'off');

%=== color according to the scheme
switch colorScheme
    case 1 % color by state
        [bp_sel, bu_sel, bp_nosel, bu_nosel] = localGroupBySelState(bp, bu, selected);

        hbp_nosel = localPlotOnCircle(theta(bp_nosel), colors.blue);
        hbu_nosel = localPlotOnCircle(theta(bu_nosel), colors.red);
        hbp_sel = localPlotOnCircle(theta(bp_sel), colors.green);
        hbu_sel = localPlotOnCircle(theta(bu_sel), colors.green);

        H.Paired = [hbp_sel; hbp_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hbp_sel; hbu_sel];

        hnot = localSetLegendAll(hbp_nosel, hbu_nosel, [], [], H.Selected, colorScheme, colors, []);

    case 2 % color by residue type
        [A_sel, C_sel, G_sel, U_sel, A_nosel, C_nosel, G_nosel, U_nosel] = ...
            localGroupBySelResidue(seq, selected);

        hA_nosel = localPlotOnCircle(theta(A_nosel), colors.cyan);
        hC_nosel = localPlotOnCircle(theta(C_nosel), colors.blue);
        hG_nosel = localPlotOnCircle(theta(G_nosel), colors.yellow);
        hU_nosel = localPlotOnCircle(theta(U_nosel), colors.red);

        hA_sel = localPlotOnCircle(theta(A_sel), colors.green);
        hC_sel = localPlotOnCircle(theta(C_sel), colors.green);
        hG_sel = localPlotOnCircle(theta(G_sel), colors.green);
        hU_sel = localPlotOnCircle(theta(U_sel), colors.green);

        H.A = [hA_sel; hA_nosel];
        H.C = [hC_sel; hC_nosel];
        H.G = [hG_sel; hG_nosel];
        H.U = [hU_sel; hU_nosel];

        H.Selected = [hA_sel; hC_sel; hG_sel; hU_sel];
        hnot = localSetLegendAll(hA_nosel, hC_nosel, hG_nosel, hU_nosel, H.Selected, colorScheme, colors, []);

    case 3 % color by pair type
        pairs = localGroupByPair(seq, matrix);
        [AU_sel, GU_sel, GC_sel, bu_sel, AU_nosel, GU_nosel, GC_nosel, bu_nosel] = ...
            localGroupBySelPair(pairs, selected, bu);

        hAU_nosel = localPlotOnCircle(theta(AU_nosel), colors.cyan);
        hGU_nosel = localPlotOnCircle(theta(GU_nosel), colors.yellow);
        hGC_nosel = localPlotOnCircle(theta(GC_nosel), colors.blue);
        hbu_nosel = localPlotOnCircle(theta(bu_nosel), colors.red);

        hAU_sel = localPlotOnCircle(theta(AU_sel), colors.green);
        hGU_sel = localPlotOnCircle(theta(GU_sel), colors.green);
        hGC_sel = localPlotOnCircle(theta(GC_sel), colors.green);
        hbu_sel = localPlotOnCircle(theta(bu_sel), colors.green);

        H.AU = [hAU_sel; hAU_nosel];
        H.GU = [hGU_sel; hGU_nosel];
        H.GC = [hGC_sel; hGC_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hAU_sel; hGU_sel; hGC_sel; hbu_sel];

        hnot = localSetLegendAll(hAU_nosel, hGU_nosel, hGC_nosel, hbu_nosel, H.Selected, colorScheme, colors, []);

end

%=== complete figure details (numbering and axis limits)
text(labelRadius * cos(theta(ticks)), labelRadius * sin(theta(ticks)), tickLabels, ...
    'HorizontalAlignment', 'center','pickableparts','none'); % numbering
set(ha, 'yLim', [-1 1.1]);
%set(ha, 'xLim', [-1 1.1]);

%=== Datatips and callbacks
datatipLabels = localCreateDatatipLabels(N, seq);
localCreateDatatips(hf, cos(theta), sin(theta), datatipLabels);

%=== Remove unwanted objects from legend
for i = 1:numel(hlines)
    hasbehavior(hlines(i),'legend', false)
end
for i = 1:numel(hnot)
    hasbehavior(hnot(i),'legend', false)
end

legend('location', 'NorthEastOutside') % default location overlaps with the plot itself

%==========================================================================

function [ha H] = localMatrix2diagram(matrix, selected,  seq,  colorScheme, colors, dot,hf,ha)
% Draw a secondary structure diagram for the structure in matrix notation.
% Each residue is represented by a labeled node, whereas backbone bonds and
% hydrogen bonds are represented by edges. The optimal layout of the
% structure is determined by the java-based FADE layout algorithm. Selected
% residue positions are highlighted.

if dot
    dotSize = 16; % size of marker for dot
end

if ismac
    fontSize = 10;
else
    fontSize = 8;
end

%=== set FADE algorithm parameter constants
numSteps = 3000;
radius =  size(matrix,1) * 20 / 6.28; % assume 20pt for each marker
canvasSize = max([50000 radius]);
%maxStepSize = 5;
springConstant = 4e-05;
edgeLength = 10;
x0 = canvasSize/2; y0 = -canvasSize/2;

%=== add two fictitious nodes for the termini (5' and 3')
M = zeros(size(matrix,1)+2);
M(2:end-1, 2:end-1) = matrix;
seq = ['5' seq '3']; %XXXXX


%=== add connectivity representing the RNA backbone
M = diag(ones(1,size(M,1)-1),1) | M;

%=== create nodes and lay them out in a circle to lessen knotting
nNodes = size(M,1);
nodes = javaArray('com.mathworks.bde.graphLayout.Node',nNodes);

angles = (1:nNodes) * 2 * pi/nNodes;
for i = 1:nNodes
    nodes(i) = com.mathworks.bde.graphLayout.Node(x0 + radius * cos(angles(i)), y0 + radius * sin(angles(i)));
end

%=== create edges according to the connectivity matrix
[m,n] = find(M);
nEdges = length(m);
edges = javaArray('com.mathworks.bde.graphLayout.Edge',nEdges);
for i = 1:nEdges
    edges(i) = com.mathworks.bde.graphLayout.Edge(nodes(m(i)),nodes(n(i)));
end

%=== perform the layout
layout = com.mathworks.bde.graphLayout.FadeLayout;
layout.configure(nodes,edges, canvasSize, canvasSize);
%layout.setMaximumStepSize(maxStepSize);
layout.setSpringConstant(springConstant);
layout.setEdgeLength(edgeLength);
for step = 1:numSteps
    layout.fadeStep
end

%=== get the coordinates of all nodes
nNodes = length(nodes);
x = zeros(nNodes,1);
y = zeros(nNodes,1);
for i = 1:nNodes
    x(i) = nodes(i).bounds.x;
    y(i) = nodes(i).bounds.y;
end

%=== prepare figure plot
%hf = figure('Color','white');
%ha = axes('Tag','Bioinfo:rnaplot:alldiagram');
hold on; axis off;

%=== plot backbone and hydrogen bonds
hbk = plot([x(m),x(n)]',[y(m),y(n)]', 'k.-', 'MarkerEdgeColor', 'white', ...
    'MarkerSize', 20, 'Clipping', 'off');

%=== setup different values for dotdiagram and diagram
if dot
    sz = dotSize;
    isdiagram = 0;
    diagCoord = [];
else
    sz = fontSize;
    isdiagram = 1;
    diagCoord = [x(1) y(1)]; % used to create legend for diagram only
end

%=== plot according to the color scheme
[b1, b2, bp, bu] = localGroupByStatus(matrix);
switch colorScheme

    case 1 % color by state
        [bp_sel, bu_sel, bp_nosel, bu_nosel] = localGroupBySelState(bp + 1, bu + 1, selected + 1);

        hbp_nosel = localPlotOnDiagram(x, y, bp_nosel, seq', colors.blue, sz, isdiagram);
        hbu_nosel = localPlotOnDiagram(x, y, bu_nosel, seq', colors.red, sz, isdiagram);
        hbp_sel = localPlotOnDiagram(x, y, bp_sel, seq', colors.green, sz, isdiagram);
        hbu_sel = localPlotOnDiagram(x, y, bu_sel, seq', colors.green, sz, isdiagram);

        H.Paired = [hbp_sel; hbp_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hbp_sel; hbu_sel];

        hnot = localSetLegendAll(hbp_nosel, hbu_nosel, [], [], H.Selected, colorScheme, colors, diagCoord);

    case 2 % color by residue
        [A_sel, C_sel, G_sel, U_sel, A_nosel, C_nosel, G_nosel, U_nosel] = ...
            localGroupBySelResidue(seq(2:end-1), selected);

        hA_nosel = localPlotOnDiagram(x, y, A_nosel + 1, seq', colors.cyan, sz, isdiagram);
        hC_nosel = localPlotOnDiagram(x, y, C_nosel + 1, seq', colors.blue, sz, isdiagram);
        hG_nosel = localPlotOnDiagram(x, y, G_nosel + 1, seq', colors.yellow, sz, isdiagram);
        hU_nosel = localPlotOnDiagram(x, y, U_nosel + 1, seq', colors.red, sz, isdiagram);

        hA_sel = localPlotOnDiagram(x, y, A_sel + 1, seq', colors.green, sz, isdiagram);
        hC_sel = localPlotOnDiagram(x, y, C_sel + 1, seq', colors.green, sz, isdiagram);
        hG_sel = localPlotOnDiagram(x, y, G_sel + 1, seq', colors.green, sz, isdiagram);
        hU_sel = localPlotOnDiagram(x, y, U_sel + 1, seq', colors.green, sz, isdiagram);
        % add one because the selection is relative to the original
        % sequence, while the node coordinates are relative to the extended
        % sequences (with the two terminal nodes)

        H.A = [hA_sel; hA_nosel];
        H.C = [hC_sel; hC_nosel];
        H.G = [hG_sel; hG_nosel];
        H.U = [hU_sel; hU_nosel];
        H.Selected = [hA_sel; hC_sel; hG_sel; hU_sel];

        hnot = localSetLegendAll(hA_nosel, hC_nosel, hG_nosel, hU_nosel, H.Selected, colorScheme, colors, diagCoord);

    case 3 % color by pair
        pairs = localGroupByPair(seq(2:end-1), matrix);
        [AU_sel, GU_sel, GC_sel, bu_sel, AU_nosel, GU_nosel, GC_nosel, bu_nosel] = ...
            localGroupBySelPair(pairs, selected, bu);

        hAU_nosel = localPlotOnDiagram(x, y, AU_nosel + 1, seq', colors.cyan, sz, isdiagram);
        hGU_nosel = localPlotOnDiagram(x, y, GU_nosel + 1, seq', colors.yellow, sz, isdiagram);
        hGC_nosel = localPlotOnDiagram(x, y, GC_nosel + 1, seq', colors.blue, sz, isdiagram);
        hbu_nosel = localPlotOnDiagram(x, y, bu_nosel + 1, seq', colors.red, sz, isdiagram);

        hAU_sel = localPlotOnDiagram(x, y, AU_sel + 1, seq', colors.green, sz, isdiagram);
        hGU_sel = localPlotOnDiagram(x, y, GU_sel + 1, seq', colors.green, sz, isdiagram);
        hGC_sel = localPlotOnDiagram(x, y, GC_sel + 1, seq', colors.green, sz, isdiagram);
        hbu_sel = localPlotOnDiagram(x, y, bu_sel + 1, seq', colors.green, sz, isdiagram);

        H.AU = [hAU_sel; hAU_nosel];
        H.GU = [hGU_sel; hGU_nosel];
        H.GC = [hGC_sel; hGC_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hAU_sel; hGU_sel; hGC_sel; hbu_sel];

        hnot = localSetLegendAll(hAU_nosel, hGU_nosel, hGC_nosel, hbu_nosel, H.Selected, colorScheme, colors, diagCoord);

end

%=== add terminal labels
text(x(1), y(1), seq(1), 'HorizontalAlignment','center', 'Color', [0 0 0], ...
    'FontSize', fontSize, 'BackgroundColor', [1 1 1],'pickableparts','none');
text(x(end), y(end), seq(end), 'HorizontalAlignment','center', 'Color', [0 0 0], ...
    'FontSize', fontSize, 'BackgroundColor', [1 1 1],'pickableparts','none');

%=== datatips and callbacks
datatipLabels = localCreateDatatipLabels(nNodes-2, seq(2:end-1));
localCreateDatatips(hf, x(2:end-1), y(2:end-1), datatipLabels);

%=== remove unwanted objects from legend
for i = 1:numel(hbk)
    hasbehavior(hbk(i),'legend', false)
end
for i = 1:numel(hnot)
    hasbehavior(hnot(i),'legend', false)
end


%==========================================================================
function [ha H] = localMatrix2graph(matrix, selected, seq, colorScheme, colors)
% Draw a graph plot for the input structure in matrix notation. Each
% residues is placed on a horizontal line according to their sequence
% position. A semi-elliptical line is drawn between residues that pair.
% Selected residue positions are highlighted.

%=== determine radius and center of ellipses
[i, j, bp, bu] = localGroupByStatus(matrix);
r = (j - i) / 2; % radii
xo = i + r;      % center of ellipse
theta = linspace(0, pi, max([max(r) 100])); % angles

%=== set up figure plot
hf = figure('Color','White');
ha = axes('Tag','Bioinfo:rnaplot:graph');
hold on; axis equal;

%=== plot semi-ellipses
he = zeros(1, numel(r));
for k = 1:numel(r)
    he(k) = plot(xo(k) + r(k) * cos(theta), 0.5 * r(k) * sin(theta), 'Color', 'black');
end

%=== color according to the scheme
switch colorScheme

    case 1 % color by state
        [bp_sel, bu_sel, bp_nosel, bu_nosel] = localGroupBySelState(bp, bu, selected);
        hbp_nosel = plot(bp_nosel, zeros(1, numel(bp_nosel)), '.', 'Color', colors.blue);
        hbu_nosel = plot(bu_nosel, zeros(1, numel(bu_nosel)), '.', 'Color', colors.red);
        hbp_sel = plot(bp_sel, zeros(1, numel(bp_sel)), '.', 'Color', colors.green);
        hbu_sel = plot(bu_sel, zeros(1, numel(bu_sel)), '.', 'Color', colors.green);

        H.Paired = [hbp_sel; hbp_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hbp_sel; hbu_sel];

        hnot = localSetLegendAll(hbp_nosel, hbu_nosel, [], [], H.Selected, colorScheme, colors, []);

    case 2 % color by residue
        [A_sel, C_sel, G_sel, U_sel, A_nosel, C_nosel, G_nosel, U_nosel] = ...
            localGroupBySelResidue(seq, selected);

        hA_nosel = plot(A_nosel, zeros(1,numel(A_nosel)), '.', 'Color', colors.cyan);
        hC_nosel = plot(C_nosel, zeros(1,numel(C_nosel)), '.', 'Color', colors.blue);
        hG_nosel = plot(G_nosel, zeros(1,numel(G_nosel)), '.', 'Color', colors.yellow);
        hU_nosel = plot(U_nosel, zeros(1,numel(U_nosel)), '.', 'Color', colors.red);

        hA_sel = plot(A_sel, zeros(1,numel(A_sel)), '.', 'Color', colors.green);
        hC_sel = plot(C_sel, zeros(1,numel(C_sel)), '.', 'Color', colors.green);
        hG_sel = plot(G_sel, zeros(1,numel(G_sel)), '.', 'Color', colors.green);
        hU_sel = plot(U_sel, zeros(1,numel(U_sel)), '.', 'Color', colors.green);

        H.A = [hA_sel; hA_nosel];
        H.C = [hC_sel; hC_nosel];
        H.G = [hG_sel; hG_nosel];
        H.U = [hU_sel; hU_nosel];
        H.Selected = [hA_sel; hC_sel; hG_sel; hU_sel];

        hnot = localSetLegendAll(hA_nosel, hC_nosel, hG_nosel, hU_nosel, H.Selected, colorScheme, colors, []);

    case 3 % color by pair
        pairs = localGroupByPair(seq, matrix);
        [AU_sel, GU_sel, GC_sel, bu_sel, AU_nosel, GU_nosel, GC_nosel, bu_nosel] = ...
            localGroupBySelPair(pairs, selected, bu);

        hAU_nosel = plot(AU_nosel, zeros(1,numel(AU_nosel)), '.', 'Color', colors.cyan);
        hGU_nosel = plot(GU_nosel, zeros(1,numel(GU_nosel)), '.', 'Color', colors.yellow);
        hGC_nosel = plot(GC_nosel, zeros(1,numel(GC_nosel)), '.', 'Color', colors.blue);
        hbu_nosel = plot(bu_nosel, zeros(1,numel(bu_nosel)), '.', 'Color', colors.red);

        hAU_sel = plot(AU_sel, zeros(1,numel(AU_sel)), '.', 'Color', colors.green);
        hGU_sel = plot(GU_sel, zeros(1,numel(GU_sel)), '.', 'Color', colors.green);
        hGC_sel = plot(GC_sel, zeros(1,numel(GC_sel)), '.', 'Color', colors.green);
        hbu_sel = plot(bu_sel, zeros(1,numel(bu_sel)), '.', 'Color', colors.green);

        H.AU = [hAU_sel; hAU_nosel];
        H.GU = [hGU_sel; hGU_nosel];
        H.GC = [hGC_sel; hGC_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hAU_sel; hGU_sel; hGC_sel; hbu_sel];

        hnot = localSetLegendAll(hAU_nosel, hGU_nosel, hGC_nosel, hbu_nosel, H.Selected, colorScheme, colors, []);

end

%=== datatips and callbacks
N = size(matrix, 1);
datatipLabels = localCreateDatatipLabels(N, seq);
localCreateDatatips(hf, 1:N, zeros(1,N), datatipLabels);

%=== remove unwanted objects from legend
for i = 1:numel(he)
    hasbehavior(he(i),'legend', false)
end
for i = 1:numel(hnot)
    hasbehavior(hnot(i),'legend', false)
end

%=== Complete figure details
xlabel('Sequence position');
ylabel('Relative distance between paired bases');
set(ha, 'Ylim',[0 max(r)]);
set(ha, 'Xlim',[0 size(matrix,1)+1]);
hold off

%==========================================================================
function [ha H] = localBracket2mountain(bracket, matrix, selected, seq, colorScheme, colors)
% Draw a two-dimensional plot where the base position is in the abscissa
% and the number of base pairs enclosing a given base is in the ordinate.
% Selected residue positions are highlighted.

openBrackets = bracket == '(';
closeBrackets = bracket == ')';
mountain = cumsum(openBrackets - closeBrackets);
mountain(openBrackets) = mountain(openBrackets)-1;
N = numel(mountain);

%=== prepare figure
hf = figure('color','w');
ha = axes('Tag','Bioinfo:rnaplot:mountain');
hold on;
hlines = plot(1:N, mountain, 'k-');

%=== determine paired and unpaired residues
[i, j, bp, bu] = localGroupByStatus(matrix);

%=== plot according to scheme
switch colorScheme

    case 1 % color by state
        [bp_sel, bu_sel, bp_nosel, bu_nosel] = localGroupBySelState(bp, bu, selected);

        hbp_nosel = plot(bp_nosel, mountain(bp_nosel), '.', 'Color', colors.blue);
        hbu_nosel = plot(bu_nosel, mountain(bu_nosel), '.', 'Color', colors.red);
        hbp_sel = plot(bp_sel, mountain(bp_sel), '.', 'Color', colors.green);
        hbu_sel = plot(bu_sel, mountain(bu_sel), '.', 'Color', colors.green);

        H.Paired = [hbp_sel; hbp_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hbp_sel; hbu_sel];

        hnot = localSetLegendAll(hbp_nosel, hbu_nosel, [], [], H.Selected, colorScheme, colors, []);

    case 2 % color by residue
        [A_sel, C_sel, G_sel, U_sel, A_nosel, C_nosel, G_nosel, U_nosel] = ...
            localGroupBySelResidue(seq, selected);

        hA_nosel = plot(A_nosel, mountain(A_nosel), '.', 'Color', colors.cyan);
        hC_nosel = plot(C_nosel, mountain(C_nosel), '.', 'Color', colors.blue);
        hG_nosel = plot(G_nosel, mountain(G_nosel), '.', 'Color', colors.yellow);
        hU_nosel = plot(U_nosel, mountain(U_nosel), '.', 'Color', colors.red);

        hA_sel = plot(A_sel, mountain(A_sel), '.', 'Color', colors.green);
        hC_sel = plot(C_sel, mountain(C_sel), '.', 'Color', colors.green);
        hG_sel = plot(G_sel, mountain(G_sel), '.', 'Color', colors.green);
        hU_sel = plot(U_sel, mountain(U_sel), '.', 'Color', colors.green);

        H.A = [hA_sel; hA_nosel];
        H.C = [hC_sel; hC_nosel];
        H.G = [hG_sel; hG_nosel];
        H.U = [hU_sel; hU_nosel];
        H.Selected = [hA_sel; hC_sel; hG_sel; hU_sel];

        hnot = localSetLegendAll(hA_nosel, hC_nosel, hG_nosel, hU_nosel, H.Selected, colorScheme, colors, []);

    case 3 % color by pair
        pairs = localGroupByPair(seq, matrix);
        [AU_sel, GU_sel, GC_sel, bu_sel, AU_nosel, GU_nosel, GC_nosel, bu_nosel] = ...
            localGroupBySelPair(pairs, selected, bu);

        hAU_nosel = plot(AU_nosel, mountain(AU_nosel), '.', 'Color', colors.cyan);
        hGU_nosel = plot(GU_nosel, mountain(GU_nosel), '.', 'Color', colors.yellow);
        hGC_nosel = plot(GC_nosel, mountain(GC_nosel), '.', 'Color', colors.blue);
        hbu_nosel = plot(bu_nosel, mountain(bu_nosel), '.', 'Color', colors.red);

        hAU_sel = plot(AU_sel, mountain(AU_sel), '.', 'Color', colors.green);
        hGU_sel = plot(GU_sel, mountain(GU_sel), '.', 'Color', colors.green);
        hGC_sel = plot(GC_sel, mountain(GC_sel), '.', 'Color', colors.green);
        hbu_sel = plot(bu_sel, mountain(bu_sel), '.', 'Color', colors.green);

        H.AU = [hAU_sel; hAU_nosel];
        H.GU = [hGU_sel; hGU_nosel];
        H.GC = [hGC_sel; hGC_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hAU_sel; hGU_sel; hGC_sel; hbu_sel];

        hnot = localSetLegendAll(hAU_nosel, hGU_nosel, hGC_nosel, hbu_nosel, H.Selected, colorScheme, colors, []);
end

%=== datatips and callbacks
datatipLabels = localCreateDatatipLabels(N, seq);
localCreateDatatips(hf, 1:N, mountain, datatipLabels);

%=== remove unwanted objects from legend
for i = 1:numel(hlines)
    hasbehavior(hlines(i),'legend', false)
end
for i = 1:numel(hnot)
    hasbehavior(hnot(i),'legend', false)
end

%=== complete figure details
xlabel('Sequence position');
ylabel('Number of enclosing base pairs');
set(ha, 'Xlim', [0 N+1], 'Ylim', [0 max(mountain)+1])
hold off;


%==========================================================================
function [ha H] = localBracket2tree(bracket, matrix, selected, seq, colorScheme, colors)
% Draw the tree representation for the structure in bracket notation.
% Leaves represent unpaired residues, while internal nodes are associated
% with basepairs. The root is a fictitious node, not associated with any
% basepair. The layout of the structure can be inferred by traversing the
% tree in pre-order. Selected residue positions are highlighted.

nodeMarkerSize = 14;

%=== add fictitious pair that will become the root
extBracket = ['(' bracket ')'];
N = numel(extBracket);
extMatrix = zeros(N,N);
extMatrix(2:end-1, 2:end-1) = matrix;
extMatrix(1,end) = 1;
selected = selected + 1;

%=== determine parent-child matrix for the tree
last = 1;
CM = zeros(N,N);      % matrix parent-child
parent = zeros(N,1);  % parent of the node i
closing = zeros(N,1); % true if i is the second nt of a bp (not represented)

for i = 1:N
    if(extBracket(i)=='(')
        parent(i) = last;
        CM(parent(i),i) = 1;
        last = i;
    elseif(extBracket(i)=='.')
        parent(i) = last;
        CM(parent(i),i) = 1;
    else
        parent(i) = parent(last);
        last = parent(i);
        closing(i) = 1;
    end
end

closing = logical(closing);

%=== remove nodes corresponding to second base in a bp
CM(closing,:) = [];
CM(:,closing) = [];
CM(1,1) = 0; % avoid biograph warning on self-connected nodes

%=== find nodes connected by edges
[m,n] = find(CM);

%=== map seq position into node ids and viceversa
[b1, b2, bp, bu] = localGroupByStatus(extMatrix);

pairMap = zeros(N,2); % map containing each res and its partner
pairMap(b1, 1:2) = [b1 b2];
pairMap(b2, 1:2) = [b2 b1];
pairMap(bu, 1:2) = [bu bu];

nodeIds = sort([b1; bu]);

pos2node = zeros(N,1);
pos2node(nodeIds) = nodeIds;
pos2node(b2) = b1;
node2pos = pairMap(nodeIds,:);

%=== create a biograph
ids = cellfun(@num2str, num2cell(nodeIds), 'UniformOutput', false);
bg = biograph(CM, ids, 'Showarrows', 'off', 'Showtextinnodes', 'none');
dolayout(bg);
xy = zeros(numel(bg.Nodes),2);
for i = 1:numel(bg.Nodes)
    xy(i,1:2) = bg.Nodes(i).Position;
end

%=== prepare figure
hf = figure('Color','w');
ha = axes('Tag','Bioinfo:rnaplot:tree');
hold on; axis off;

%=== plot tree edges
he = plot([xy(m,1),xy(n,1)]', [xy(m,2),xy(n,2)]', 'k-', 'Clipping', 'off');

%=== plot according to color scheme
switch colorScheme

    case 1 % color by state
        [bp_sel, bu_sel, bp_nosel, bu_nosel] = localGroupBySelState(bp, bu, selected);

        %=== remove root from set of basepairs
        %bp_nosel((bp_nosel == 1 | bp_nosel == N)) = [];

        %=== map sequence positions into node ids
        bu_noselIds = cellfun(@num2str, num2cell(pos2node(bu_nosel)), 'UniformOutput', false);
        bp_noselIds = cellfun(@num2str, num2cell(pos2node(bp_nosel)), 'UniformOutput', false);
        bu_selIds = cellfun(@num2str, num2cell(pos2node(bu_sel)), 'UniformOutput', false);
        bp_selIds = cellfun(@num2str, num2cell(pos2node(bp_sel)), 'UniformOutput', false);

        buNodes_nosel = getnodesbyid(bg, bu_noselIds);
        bpNodes_nosel = getnodesbyid(bg, bp_noselIds);
        buNodes_sel = getnodesbyid(bg, bu_selIds);
        bpNodes_sel = getnodesbyid(bg, bp_selIds);

        %=== get node coordinates and plot
        hbu_nosel = localPlotOnTree(buNodes_nosel, colors.red, nodeMarkerSize);
        hbp_nosel = localPlotOnTree(bpNodes_nosel, colors.blue, nodeMarkerSize);
        hbu_sel = localPlotOnTree(buNodes_sel, colors.green, nodeMarkerSize);
        hbp_sel = localPlotOnTree(bpNodes_sel, colors.green, nodeMarkerSize);

        %=== return handles
        H.Paired = [hbp_sel; hbp_nosel];
        H.Unpaired = [hbu_sel; hbu_nosel];
        H.Selected = [hbp_sel; hbu_sel];

        hnot = localSetLegendAll(hbp_nosel, hbu_nosel, [], [], H.Selected, colorScheme, colors, []);

    case 2 % color by residues
        [A_sel, C_sel, G_sel, U_sel, A_nosel, C_nosel, G_nosel, U_nosel] = ...
            localGroupBySelResidue(seq, selected);

        %=== map sequence positions into node ids
        A_selIds = cellfun(@num2str, num2cell(pos2node(A_sel)), 'UniformOutput', false);
        C_selIds = cellfun(@num2str, num2cell(pos2node(C_sel)), 'UniformOutput', false);
        G_selIds = cellfun(@num2str, num2cell(pos2node(G_sel)), 'UniformOutput', false);
        U_selIds = cellfun(@num2str, num2cell(pos2node(U_sel)), 'UniformOutput', false);

        A_noselIds = cellfun(@num2str, num2cell(pos2node(A_nosel)), 'UniformOutput', false);
        C_noselIds = cellfun(@num2str, num2cell(pos2node(C_nosel)), 'UniformOutput', false);
        G_noselIds = cellfun(@num2str, num2cell(pos2node(G_nosel)), 'UniformOutput', false);
        U_noselIds = cellfun(@num2str, num2cell(pos2node(U_nosel)), 'UniformOutput', false);

        ANodes_sel = getnodesbyid(bg, A_selIds);
        CNodes_sel = getnodesbyid(bg, C_selIds);
        GNodes_sel = getnodesbyid(bg, G_selIds);
        UNodes_sel = getnodesbyid(bg, U_selIds);

        ANodes_nosel = getnodesbyid(bg, A_noselIds);
        CNodes_nosel = getnodesbyid(bg, C_noselIds);
        GNodes_nosel = getnodesbyid(bg, G_noselIds);
        UNodes_nosel = getnodesbyid(bg, U_noselIds);

        %=== get node coordinates and plot
        hA_nosel = localPlotOnTree(ANodes_nosel, colors.cyan, nodeMarkerSize);
        hC_nosel = localPlotOnTree(CNodes_nosel, colors.blue, nodeMarkerSize);
        hG_nosel = localPlotOnTree(GNodes_nosel, colors.yellow, nodeMarkerSize);
        hU_nosel = localPlotOnTree(UNodes_nosel, colors.red, nodeMarkerSize);

        hA_sel = localPlotOnTree(ANodes_sel, colors.green, nodeMarkerSize);
        hC_sel = localPlotOnTree(CNodes_sel, colors.green, nodeMarkerSize);
        hG_sel = localPlotOnTree(GNodes_sel, colors.green, nodeMarkerSize);
        hU_sel = localPlotOnTree(UNodes_sel, colors.green, nodeMarkerSize);

        H.A = [hA_sel; hA_nosel];
        H.C = [hC_sel; hC_nosel];
        H.G = [hG_sel; hG_nosel];
        H.U = [hU_sel; hU_nosel];
        H.Selected = [hA_sel; hC_sel; hG_sel; hU_sel];

        hnot = localSetLegendAll(hA_nosel, hC_nosel, hG_nosel, hU_nosel, H.Selected, colorScheme, colors, []);

    case 3 % color by pair
        pairs = localGroupByPair(seq, matrix);
        pairs{1} = pairs{1} + 1; % add one because of the root
        pairs{2} = pairs{2} + 1;
        pairs{3} = pairs{3} + 1;

        [AU_sel, GU_sel, GC_sel, bu_sel, AU_nosel, GU_nosel, GC_nosel, bu_nosel] = ...
            localGroupBySelPair(pairs, selected, bu);

        %=== map sequence positions into node ids
        AU_selIds = cellfun(@num2str, num2cell(pos2node(AU_sel)), 'UniformOutput', false);
        GU_selIds = cellfun(@num2str, num2cell(pos2node(GU_sel)), 'UniformOutput', false);
        GC_selIds = cellfun(@num2str, num2cell(pos2node(GC_sel)), 'UniformOutput', false);
        bu_selIds = cellfun(@num2str, num2cell(pos2node(bu_sel)), 'UniformOutput', false); % bu comes from extMatrix (root is included)

        AU_noselIds = cellfun(@num2str, num2cell(pos2node(AU_nosel)), 'UniformOutput', false);
        GU_noselIds = cellfun(@num2str, num2cell(pos2node(GU_nosel)), 'UniformOutput', false);
        GC_noselIds = cellfun(@num2str, num2cell(pos2node(GC_nosel)), 'UniformOutput', false);
        bu_noselIds = cellfun(@num2str, num2cell(pos2node(bu_nosel)), 'UniformOutput', false);

        AUNodes_sel = getnodesbyid(bg, AU_selIds);
        GUNodes_sel = getnodesbyid(bg, GU_selIds);
        GCNodes_sel = getnodesbyid(bg, GC_selIds);
        buNodes_sel = getnodesbyid(bg, bu_selIds);

        AUNodes_nosel = getnodesbyid(bg, AU_noselIds);
        GUNodes_nosel = getnodesbyid(bg, GU_noselIds);
        GCNodes_nosel = getnodesbyid(bg, GC_noselIds);
        buNodes_nosel = getnodesbyid(bg, bu_noselIds);

        %=== get node coordinates and plot
        hAU_nosel = localPlotOnTree(AUNodes_nosel, colors.cyan, nodeMarkerSize);
        hGU_nosel = localPlotOnTree(GUNodes_nosel, colors.yellow, nodeMarkerSize);
        hGC_nosel = localPlotOnTree(GCNodes_nosel, colors.blue, nodeMarkerSize);
        hbu_nosel = localPlotOnTree(buNodes_nosel, colors.red, nodeMarkerSize);

        hAU_sel = localPlotOnTree(AUNodes_sel, colors.green, nodeMarkerSize);
        hGU_sel = localPlotOnTree(GUNodes_sel, colors.green, nodeMarkerSize);
        hGC_sel = localPlotOnTree(GCNodes_sel, colors.green, nodeMarkerSize);
        hbu_sel = localPlotOnTree(buNodes_sel, colors.green, nodeMarkerSize);

        H.AU = [hAU_sel; hAU_nosel];
        H.GU = [hGU_sel; hGU_nosel];
        H.GC = [hGC_sel; hGC_nosel];
        H.bu = [hbu_sel; hbu_nosel];
        H.Selected = [hAU_sel; hGU_sel; hGC_sel; hbu_sel];

        hnot = localSetLegendAll(hAU_nosel, hGU_nosel, hGC_nosel, hbu_nosel, H.Selected, colorScheme, colors, []);
end

%=== datatips and callbacks
datatipLabels = cell(1, numel(bg.Nodes));
if isempty(seq) % labels are just sequence positions
    for k = 2:numel(bg.Nodes)
        if node2pos(k,1) ==  node2pos(k,2) % unpaired base
            datatipLabels{k} = num2str(node2pos(k,1) - 1);
        else % paired base
            datatipLabels{k} = [num2str(node2pos(k,1) - 1 ) ' - ' num2str(node2pos(k,2) - 1)];
        end
    end

else
    seq = [' ' seq ' '];
    for k = 2:numel(bg.Nodes)
        if  node2pos(k,1) ==   node2pos(k,2) % unpaired base
            datatipLabels{k} = [seq(node2pos(k,1)) num2str( node2pos(k,1) - 1)];
        else % paired base
            datatipLabels{k} = [seq( node2pos(k,1)) num2str( node2pos(k,1) - 1 ) ...
                ' - ' seq( node2pos(k,2)) num2str( node2pos(k,2) - 1)];
        end
    end
end

localCreateDatatips(hf, xy(:,1), xy(:,2), datatipLabels);

%=== remove unwanted objects from legend
for i = 1:numel(he)
    hasbehavior(he(i),'legend', false)
end
for i = 1:numel(hnot)
    hasbehavior(hnot(i),'legend', false)
end

hold off;

%==========================================================================
% AUXILIARY SUBFUNCTIONS
%==========================================================================

%==========================================================================
function [a, c, g, u] = localGroupByType(labels)
% Given a sequence, group the residues by type (A,C,G,T) and return
% their positions.

numLabels = nt2int(labels);
a = find(numLabels == 1);
c = find(numLabels == 2);
g = find(numLabels == 3);
u = find(numLabels == 4);

%==========================================================================
function [i,j,p,u] = localGroupByStatus(str)
% Given a matrix representing the secondary structure, group the residues
% by pair status (paired vs. unpaired) and return their positions.

[i,j] = find(triu(str));
p = sort([i; j]);                    % paired bases
u = setdiff(1:size(str,1),p)';       % unpaired bases

%==========================================================================
function pairs = localGroupByPair(seq, matrix)
% Given a sequence and the secondary structure in a matrix, group the
% residues by pair type (AU/UA, GC/CG, GU/UG) and return their positions.

seqn = nt2int(seq);
cseqn = 5 - seqn;
[X,Y]= meshgrid(seqn, cseqn); 

AUs = (X == 4 & matrix & Y == 4);
CGs = (X == 3 & matrix & Y == 3);
GCs = (X == 2 & matrix);
UAs = (X == 1 & matrix);

matrix = triu(matrix & (X + X' == 7));
GUs = X == 4 & matrix & Y == 2;
UGs = X == 3 & matrix & Y == 1;

[AU1 AU2] = find(AUs | UAs);
[GU1 GU2] = find(GUs | UGs);
[GC1 GC2] = find(GCs | CGs);

pairs{1} = [AU1 AU2];
pairs{2} = [GU1 GU2];
pairs{3} = [GC1 GC2];

%==========================================================================
function datatipLabels = localCreateDatatipLabels(N, seq)
% Create N datatip labels. If sequence is provided, each residue is labeled
% by the residue type and the residue position, otherwise only the position
% is used.

if isempty(seq)  % no seq info provided
    datatipLabels = cellfun(@num2str, num2cell((1:N)'), 'UniformOutput', false);
else
    positions = cellfun(@num2str, num2cell((1:N)'), 'UniformOutput', false);
    residues = cellstr(seq');
    datatipLabels = strcat(residues, positions);
end

%==========================================================================
function localCreateDatatips(hf, x, y, labels)
% Create datatip in axes ha at positions given by coordinates x and y,
% using the provided labels.

appdata.x = x;
appdata.y = y;
appdata.labels = labels;
appdata.hAxis = hf.CurrentAxes;
appdata.hTip = text(0,0,'','BackgroundColor', [1 1 0.933333],'Color', [0 0 0],...
    'EdgeColor', [0.8 0.8 0.8], 'Visible','off',...
    'Tag','Bioinfo:rnaplot:dataTip','Interpreter','none',...
    'fontsize',9,'fontname','Helvetica', 'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center');
set(findall(appdata.hAxis,'Marker','.'),'ButtonDownFcn',@clickOnPlot);
set(hf,'WindowButtonUpFcn',@(h,e) set(appdata.hTip,'Visible','off'));
setappdata(hf,'RNAPlot',appdata)

%==========================================================================
function clickOnPlot(varargin)
% callback function highlights selected element and displays label.
hf = gcbf;

appdata = getappdata(hf,'RNAPlot');
point = get(appdata.hAxis,'CurrentPoint');

% find the closest point
axesExtentInPoints = hgconvertunits(hf,appdata.hAxis.Position,appdata.hAxis.Units,'Points',hf);
scaleX = axesExtentInPoints(1)./ diff(appdata.hAxis.XLim);
scaleY = axesExtentInPoints(2)./ diff(appdata.hAxis.YLim);

[~, index] = min( ( scaleX.*(appdata.x - point(1,1)) ).^2 +...
                  ( scaleY.*(appdata.y - point(1,2)) ).^2 );
			  
set(appdata.hTip,'String',appdata.labels{index},'visible','on')
set(appdata.hTip,'Position',point(1,:))

%==========================================================================
function [bp_sel, bu_sel, bp_nosel, bu_nosel] = localGroupBySelState(bp, bu, selected)
% Group residue positions into: 1) basepaired in the selection; 2) unpaired
% in the selection; 3) basepaired not in the selection; 4) unpaired not in
% the selection

if ~isempty(selected)
    bp_sel = intersect(bp, selected); % basepaired bases in selection
    bu_sel = intersect(bu, selected); % unpaired bases in selection
    bp_nosel = setdiff(bp, selected); % basepaired bases not in selection
    bu_nosel = setdiff(bu, selected); % unpaired bases not in selection
else
    bp_sel = [];
    bu_sel = [];
    bp_nosel = bp;
    bu_nosel = bu;
end
%==========================================================================
function  [A_sel, C_sel, G_sel, U_sel, A_nosel, C_nosel, G_nosel, U_nosel] = ...
    localGroupBySelResidue(seq, selected)
% Group residues by type and depending whether they are in the selection or not.

[As, Cs, Gs, Us] = localGroupByType(seq);

if ~isempty(selected)
    A_sel = intersect(As, selected);
    C_sel = intersect(Cs, selected);
    G_sel = intersect(Gs, selected);
    U_sel = intersect(Us, selected);

    A_nosel = setdiff(As, selected);
    C_nosel = setdiff(Cs, selected);
    G_nosel = setdiff(Gs, selected);
    U_nosel = setdiff(Us, selected);
else
    A_sel = [];
    C_sel = [];
    G_sel = [];
    U_sel = [];

    A_nosel = As;
    C_nosel = Cs;
    G_nosel = Gs;
    U_nosel = Us;
end

%==========================================================================
function [AU_sel, GU_sel, GC_sel, bu_sel, AU_nosel, GU_nosel, GC_nosel, bu_nosel] = ...
    localGroupBySelPair(pairs, selected, bu)
% Group residues by pair type and whether they are in the selection or not.

AUs = pairs{1};
GUs = pairs{2};
GCs = pairs{3};

if ~isempty(selected)
    AU_sel = intersect([AUs(:,1); AUs(:,2)], selected);
    GU_sel = intersect([GUs(:,1); GUs(:,2)], selected);
    GC_sel = intersect([GCs(:,1); GCs(:,2)], selected);
    bu_sel = intersect(bu, selected);

    AU_nosel = setdiff([AUs(:,1); AUs(:,2)], selected);
    GU_nosel = setdiff([GUs(:,1); GUs(:,2)], selected);
    GC_nosel = setdiff([GCs(:,1); GCs(:,2)], selected);
    bu_nosel = setdiff(bu, selected);
else
    AU_sel = [];
    GU_sel = [];
    GC_sel = [];
    bu_sel = [];

    AU_nosel = [AUs(:,1); AUs(:,2)];
    GU_nosel = [GUs(:,1); GUs(:,2)];
    GC_nosel = [GCs(:,1); GCs(:,2)];
    bu_nosel = bu;
end

%==========================================================================
function q = localPlotOnTree(nodes, mycolor, nodeMarkerSize)
% Plot tree using the coordinates from nodes (bg object)

if isempty(nodes)
    q = [];
else
    xy = zeros(numel(nodes, 2));
    for i = 1:numel(nodes)
        xy(i,1:2) = nodes(i).Position;
    end
    q = plot(xy(:,1),xy(:,2), '.', 'Color', mycolor, ...
        'MarkerSize', nodeMarkerSize, 'Clipping', 'off');
end
%==========================================================================
function q = localPlotOnCircle(x, mycolor)
% Plot on a circle using angles given by x

q = plot(cos(x), sin(x), '.', 'Color', mycolor,'Clipping', 'off');

%==========================================================================
function q = localPlotOnDiagram(x, y, z, seq, mycolor, mysize, isdiagram)
% Plot diagram or dotdiagram using coordinates x and y of subset z

if isdiagram
    q = text(x(z), y(z), seq(z), 'HorizontalAlignment', 'Center', ...
        'FontSize', mysize, 'Color', mycolor, 'Clipping', 'off','pickableparts','none');
else % dot
    q = plot(x(z), y(z), '.', 'MarkerSize', mysize, 'Color', mycolor, 'Clipping', 'off');
end

%==========================================================================
function hnot = localSetLegendAll(h1, h2, h3, h4, h5, colorScheme, colors, diagCoord)
% Create legend for handles h1...h5, using the appropriate color scheme. If
% the plot is a diagram, an artificial legend is created using objects
% plotted behind the 5' text with coordinates diagCoord. Return handle to
% objects that will have to be removed from the legend.

switch colorScheme
    case 1
        legendStrings = {'Paired', 'Unpaired', '', '', 'Selected'};
        legendColors = {colors.blue, colors.red, '', '', colors.green};
    case 2
        legendStrings = {'A', 'C', 'G', 'U', 'Selected'};
        legendColors = {colors.cyan, colors.blue, colors.yellow, colors.red, colors.green};
    case 3
        legendStrings = {'AU/UA', 'GU/UG', 'GC/CG', 'Unpaired', 'Selected'};
        legendColors = {colors.cyan, colors.yellow, colors.blue, colors.red, colors.green};
end

hnot = []; % in case h5 is empty
hLeg = zeros(1,5);

if ~isempty(h1) % Paired | A | AU
    if isempty(diagCoord)
        hLeg(1) = h1(1);
    else
        hLeg(1) = plot(diagCoord(1), diagCoord(2), 'o', 'Color', legendColors{1}, 'Visible', 'on');
    end
end

if ~isempty(h2) % Unpaired | C | GU
    if isempty(diagCoord)
        hLeg(2) = h2(1);
    else
        hLeg(2) = plot(diagCoord(1), diagCoord(2), 'o', 'Color', legendColors{2}, 'Visible', 'on');
    end
end

if ~isempty(h3) % [] | G | GC
    if isempty(diagCoord)
        hLeg(3) = h3(1);
    else
        hLeg(3) = plot(diagCoord(1), diagCoord(2), 'o', 'Color', legendColors{3}, 'Visible', 'on');
    end
end

if ~isempty(h4) % [] | U | Unpaired
    if isempty(diagCoord)
        hLeg(4) = h4(1);
    else
        hLeg(4) = plot(diagCoord(1), diagCoord(2), 'o', 'Color', legendColors{4}, 'Visible', 'on');
    end
end

if ~isempty(h5) % Selected
    if isempty(diagCoord)
        hLeg(5) = h5(1);
    else
        hLeg(5) = plot(diagCoord(1), diagCoord(2), 'o', 'Color', legendColors{5}, 'Visible', 'on');
    end
    if length(h5) > 1
        hnot = h5(2:end);
    end
end

isEmptyHandle = hLeg ~= 0;
legend(hLeg(isEmptyHandle), legendStrings(isEmptyHandle));
legend off; %XXXX
