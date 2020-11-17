% get the clade mapping 
subset = 0.01:0.01:0.99


f = fopen('../../ncov-severity/across-states/strain_to_clade.tsv');fgets(f);

c = 1;
strain.id = cell(0,0);strain.strain = cell(0,0);
while ~feof(f)
    line = strsplit(strtrim(fgets(f)));
    strain.strain{c} = line{1};
    strain.id{c} = line{2};
    c = c+1;
end
fclose(f);h

% exclude WA1 
exlude = {'USA/WA1/2020'};
sample_cutoff = '2020-07-23';
clade = {'all', 'D', 'G', 'Yakima'};

end_date = '2020-01-31';

s = fopen('../results/mrsi.tsv', 'w');
cls = fopen('../results/cluster_size.tsv', 'w');
st = fopen('../results/sampling_times.tsv', 'w');

fprintf(s,'filename\tmrsi\tclade\n');
fprintf(cls,'filename\tnumber\tsize\tclade\n');
fprintf(st,'Date\tnumber\n');

%% distance based clustering
date_cutoff = datenum(sample_cutoff);
% read in iq tree from nextstrain pipeline
full_tree = phytreeread('../../ncov/results/wa_state/tree.nwk');
% get all WA sequences
f = fopen('../data/combined_meta.tsv');
line = strsplit(fgets(f), '\t');
div_id = find(ismember(line,'division'));
date_id = find(ismember(line,'date'));
originating_id = find(ismember(line,'originating_lab'));
submitting_id = find(ismember(line,'submitting_lab'));
location_id = find(ismember(line,'location'));
c=1;
id = cell(0,0);
clade_membership = cell(0,0);
division = cell(0,0);
lab = cell(0,0);
sub = cell(0,0);
date = cell(0,0);
date_val = zeros(0,0);
while ~feof(f)
    line = strsplit(fgets(f), '\t','CollapseDelimiters', false);
    id{c,1} = line{1};
    strain_id = find(ismember(strain.id,line{1}));
    if ~isempty(strain_id)
        clade_membership{c,1} = strain.strain{strain_id};
    end
    date{c,1} = line{date_id};
    lab{c,1} = line{originating_id};
    sub{c,1} = line{submitting_id};


    if ~isempty(find(ismember(id{c,1}, exlude)))
        division{c,1} = 'NA';
    elseif ~contains(line{date_id}, 'X') && sum(date{c,1}=='-')==2 && datenum(date{c,1})<=date_cutoff
        date_val(c,1) = datenum(date{c,1});
        division{c,1} = line{div_id};
    else
        division{c,1} = 'NA';
    end
    c=c+1;
end
fclose(f);

for rep = 1 : length(subset)
    

    % get all leafnames
    leafs = get(full_tree, 'leafnames');

    %%
    is_wa = find(ismember(division, 'Washington'));

    is_wa_leaf = false(length(leafs), 1);
    for i = 1 : length(leafs)
        if sum(ismember(id(is_wa), leafs{i}))==1
            is_wa_leaf(i) = true;
        end
    end
    
    keep_nodes = true(length(leafs), 1);
    keep_nodes(randsample(find(is_wa_leaf==0), round(sum(is_wa_leaf==0)*(1-subset(rep))))) = false;

    tree = prune(full_tree, ~keep_nodes);
    
    leafs = get(tree, 'leafnames');
    is_wa = find(ismember(division, 'Washington'));

    is_wa_leaf = false(length(leafs), 1);
    for i = 1 : length(leafs)
        if sum(ismember(id(is_wa), leafs{i}))==1
            is_wa_leaf(i) = true;
        end
    end

    
    nodenames = get(tree, 'nodenames');    

    % get the connectivity matrix and the distances between nodes
    mat = getmatrix(tree);
    dist = pdist(tree, 'Squareform',true,'Nodes','all');

    % initialize the location vector
    location = cell(length(nodenames),1);
    visited = false(length(nodenames),1);

    for j = 1 : length(is_wa_leaf)           
        visited(j) = true;
        if is_wa_leaf(j)
            location{j} = "WA";         
        else
            location{j} = "Outside";
        end
    end

    % upwards path of the parsimony calling    
    not_visited = find(~visited);
    while ~isempty(not_visited)
        disp(length(not_visited))
        not_visited = find(~visited);
        for j = 1 : length(not_visited)

            children = find(mat(not_visited(j),:));
            if ~isempty(location{children(1)}) && ~isempty(location{children(2)})
                % get the distance of the new sample from a basel sequence
                int = intersect(location{children(1)}, location{children(2)});
                if isempty(int)
                    location{not_visited(j)} = [location{children(1)}, location{children(2)}];
                else
                    location{not_visited(j)} = int;
                end
                visited(not_visited(j)) = true; 
            end
        end
    end

    % downwards calling 
    visited(length(is_wa_leaf)+1:end) = false;
    visited(end) = true;

    not_visited = find(~visited);
    while ~isempty(not_visited)
        not_visited = find(~visited);
        for j = length(not_visited) : -1 : 1
            parent = find(mat(:,not_visited(j)));            
            if sum(ismember(not_visited, parent))==0
                int = intersect(location{not_visited(j)}, location{parent});
                if ~isempty(int)
                    location{not_visited(j)} = int;
                end
                visited(not_visited(j)) = true;
            end
        end
    end


    %%

    % get all nodes that are only in WA
    onlyWA = false(length(visited),1);   
    uncertainCOunt = 0;
    for j = 1 : length(onlyWA)
        if length(location{j})==1 && location{j}=="WA"
            onlyWA(j) = true;
        elseif location{j}=="WA"
            uncertainCOunt = uncertainCOunt+1;
        end
    end


    isInWA = find(onlyWA(1:length(is_wa_leaf)));

    wa_names = str2double(nodenames(isInWA));    

    clustered_node = false(size(nodenames));
    clustered_node(isInWA) = true;


    % for each isInWA leaf, get all parent nodes in WA
    waParent = cell(length(isInWA),1);
    for a = 1 : length(isInWA)
        parent = find(mat(:,isInWA(a)));
        while onlyWA(parent)
            clustered_node(parent) = true;
            waParent{a} = [waParent{a} parent];
            parent = find(mat(:,parent));
        end
    end

    in_cluster_mat = zeros(length(isInWA), length(isInWA));

    new_mat = zeros(size(in_cluster_mat));
    % get for each pair of leaves if they are in the same cluster
    for a = 1 : length(isInWA)
        disp(a)
        if ~isempty(waParent{a})
            for b = a+1 : length(isInWA)
                if ~isempty(waParent{b})
                    int = intersect(waParent{a}, waParent{b});
                    if ~isempty(int)
                        new_mat(a,b) = 1;
                    end
                end
            end     
        end
    end

    %%
    wa_leafs = leafs(is_wa_leaf);

    connected = new_mat;

    stop=false;
    while ~stop
        [a,b] = find(connected==1,1,'first');
        if isempty(a)
            stop = true;
        else
            % combine a and b
            wa_leafs{a} = [wa_leafs{a} ',' wa_leafs{b}];
            % combine the rows of a and b
            for i = 1 : length(wa_leafs)
                connected(a,i) = max([connected(a,i),connected(b,i)]);
            end
            % remove b
            connected(b,:) = 0;
            connected(:,b) = 0;
            wa_leafs{b} = 'NA';
        end
    end

    % get the cluster
    cl_ind = find(~ismember(wa_leafs, 'NA'));
    wa_clusters = wa_leafs(cl_ind);


    f = fopen(['../results/cluster_assignment_' num2str(rep) '.tsv'], 'w');
    fprintf(f,'strain\tcluster\n')
    for a = 1 : length(wa_clusters)
        seqs = strsplit(wa_clusters{a}, ',');
        for b = 1 : length(seqs)
            if length(seqs{b})>2
                fprintf(f, '%s\t%d\n', seqs{b}, a)
            end
        end
    end
    fclose(f);
end

