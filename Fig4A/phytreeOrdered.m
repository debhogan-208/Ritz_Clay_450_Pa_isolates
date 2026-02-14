% %% David Ritz
% PhD candidate | Schultz lab
% Created: Summer 2025
% Description: Write a phytree to Newick format with custom leaf order,
% bypassing Matlab's internal traversal/optimization
%
%   phytreeOrdered(tree, leaves, filename)
%
%   tree     : MATLAB phytree object
%   leaves   : cell array of leaf names, in desired order
%   filename : output .nwk file
%
% Example:
%   phytreeOrdered(smallTree, leaves, 'custom_tree.nwk')

function phytreeOrdered(tree, leaves, filename)

    % Extract tree data
    names     = get(tree, 'LeafNames');
    distances = get(tree, 'Distances');
    pointers  = get(tree, 'Pointers');

    % Root node index (last internal node)
    nLeaves   = length(names);
    nNodes    = size(pointers,1) + nLeaves;
    root      = nNodes;  

    % Build Newick string
    newick = [toNewick(root, pointers, names, distances, leaves, nLeaves), ';'];

    % Write to file
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', newick);
    fclose(fid);
end

%% Recursive function to traverse tree
function newickStr = toNewick(node, pointers, names, distances, leaves, nLeaves)
    if node <= nLeaves
        % Leaf node
        newickStr = sprintf('%s:%.6f', names{node}, distances(node));
    else
        % Internal node: its row in pointers is node - nLeaves
        children = pointers(node - nLeaves, :);

        % Enforce custom leaf order
        children = orderChildren(children, pointers, names, leaves, nLeaves);

        % Build strings for each child
        childStrs = cell(1, numel(children));
        for k = 1:numel(children)
            childStrs{k} = toNewick(children(k), pointers, names, distances, leaves, nLeaves);
        end

        newickStr = sprintf('(%s):%.6f', strjoin(childStrs, ','), distances(node));
    end
end

%% Helper to reorder children according to `leaves`
function ordered = orderChildren(children, pointers, names, leaves, nLeaves)
    orderVals = zeros(size(children));

    for i = 1:numel(children)
        descLeaves = getDescendantLeaves(children(i), pointers, names, nLeaves);
        [~, pos]   = ismember(descLeaves, leaves);
        orderVals(i) = min(pos(pos > 0)); % choose earliest matching leaf index
    end

    [~, idx] = sort(orderVals);
    ordered  = children(idx);
end

%% Helper: recursively collect all descendant leaves for a node
function desc = getDescendantLeaves(node, pointers, names, nLeaves)
    if node <= nLeaves
        % Leaf node
        desc = {names{node}};
    else
        children = pointers(node - nLeaves, :);
        desc = {};
        for k = 1:numel(children)
            desc = [desc, getDescendantLeaves(children(k), pointers, names, nLeaves)]; %#ok<AGROW>
        end
    end
end