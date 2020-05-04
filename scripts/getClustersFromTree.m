function [cluster_member] = getClustersFromTree(tree)
%Gets the different clusters, based on introductions from a single tree
tree_tmp = sprintf(strrep(tree,')', ')n%d'), [1:length(strfind(tree, ')'))]);
tree_tmp = regexprep(tree_tmp, 'E-(\d*)','');

tree_tmp = regexprep(tree_tmp,'\[&type="(\w*)",reaction="Introduction",time=(\d*).(\d*)\]','_outside');
tree_tmp = regexprep(tree_tmp,'\[&type="(\w)",reaction="Infection",time=(\d*).(\d*)\]','');
tree_tmp = regexprep(tree_tmp,'\[&type="(\w)",reaction="Sampling",time=(\d*).(\d*)\]','');
% tree_tmp = regexprep(tree_tmp,'\[&type="(\w)",time=(\d*).(\d*)\]','');
% tree_tmp = regexprep(tree_tmp,'\[&type="(\w*)",location="(\d*)",reaction="(\w*)",time=(\d*).(\d*)\]','_$2');
% tree_tmp = regexprep(tree_tmp,'\[&type="Istart",reaction="Start",time=(\d*).(\d*)\]','_1');

% read in as matlab tree object
ptree = phytreeread(tree_tmp);
nodenames = get(ptree,'nodenames');
matrix = getmatrix(ptree);


% get all nodes not in 0
make_subtrees = false(1,length(nodenames));
for i = 1:length(nodenames)
    if contains(nodenames{i}, '_outside')
        children = find(matrix(i,:));
        if ~contains(nodenames{children(1)}, '_outside')
            make_subtrees(children(1)) = true;
        end
        if ~contains(nodenames{children(2)}, '_outside')
            make_subtrees(children(2)) = true;
        end
    end
end
% get the indices for the different subtrees
indices = find(make_subtrees);

leafs = cell(0,0);
cluster_member = zeros(0,0);

for i = 1 : length(indices)
    if indices(i) <= (length(nodenames)-1)/2 +1
        cluster_member(str2double(nodenames{indices(i)})) = i;
    else
        sub = subtree(ptree, indices(i));
        cl_leafs = get(sub,'leafnames');
        for j = 1 : length(cl_leafs)
            cluster_member(str2double(cl_leafs(j))) = i;
        end
    end
end
end

