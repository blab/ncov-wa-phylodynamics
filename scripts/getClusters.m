threshold = 5;

% get the clade mapping 
f = fopen('../../ncov-severity/across-states/strain_to_clade.tsv');fgets(f);

c = 1;
strain.id = cell(0,0);strain.strain = cell(0,0);
while ~feof(f)
    line = strsplit(strtrim(fgets(f)));
    strain.strain{c} = line{1};
    strain.id{c} = line{2};
    c = c+1;
end
fclose(f);

% some random keyboard hits for the random number generator
rng(108978544);
% define how many samples to take that are not from UW
max_nr_samples = [1200 750 750 750];

% exclude WA1 
exlude = {'USA/WA1/2020'};
sample_cutoff = '2020-07-03';
clade = {'all', 'D', 'G', 'Yakima'};

% read in the number of positive tests per county
f = fopen('../data/PUBLIC_Tests_by_Specimen_Collection.csv');fgets(f);
while ~feof(f)
    line = strsplit(fgets(f), ',');
    if contains(line{1}, 'County')
        
    end
end


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
tree = phytreeread('../../ncov/results/wa_state/tree_raw.nwk');
% get all WA sequences
f = fopen('../../ncov/data/example_metadata.tsv');
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
location = cell(0,0);
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
        location{c,1} = 'NA';
    elseif ~contains(line{date_id}, 'X') && sum(date{c,1}=='-')==2 && datenum(date{c,1})<=date_cutoff
        date_val(c,1) = datenum(date{c,1});
        division{c,1} = line{div_id};
        location{c,1} = line{location_id};
    else
        division{c,1} = 'NA';
        location{c,1} = 'NA';
    end
    c=c+1;
end
fclose(f);

% get all leafnames
leafs = get(tree, 'leafnames');

pariwise_distances = pdist(tree, 'Squareform', true)*29000;


%%
is_wa = find(ismember(division, 'Washington'));

is_wa_leaf = false(length(leafs), 1);
for i = 1 : length(leafs)
    if sum(ismember(id(is_wa), leafs{i}))==1
        is_wa_leaf(i) = true;
    end
end

wa_leafs = leafs(is_wa_leaf);
% get the pairwise genetic distance between all samples (times the 
% approximate number of nucleotides in the virus
% take only the wa leafs
wa_distances = pariwise_distances(is_wa_leaf,is_wa_leaf);
% consider two sequences linked if they ahve less than 3 mutations
% difference
% set all diagonal entries to large
for i = 1 : length(wa_leafs)
    wa_distances(i,i) = threshold*2;
end
below_threshold = wa_distances<threshold;
below_threshold = triu(wa_distances<threshold);

stop=false;
while ~stop
    [a,b] = find(below_threshold==1,1,'first');
    if isempty(a)
        stop = true;
    else
        % combine a and b
        wa_leafs{a} = [wa_leafs{a} ',' wa_leafs{b}];
        % combine the rows of a and b
        for i = 1 : length(wa_leafs)
            below_threshold(a,i) = max([below_threshold(a,i),below_threshold(b,i)]);
        end
        % remove b
        below_threshold(b,:) = 0;
        below_threshold(:,b) = 0;
        wa_leafs{b} = 'NA';
    end
end
% get the cluster
cl_ind = find(~ismember(wa_leafs, 'NA'));
wa_clusters = wa_leafs(cl_ind);


% remove clusters if they are not in the correct clade
clade_label = cell(length(wa_clusters),1);
for a = length(wa_clusters):-1:1
    seqs = strsplit(wa_clusters{a}, ',');
    is_in_clade = cell(0,0);
    for b = 1 : length(seqs)
        ind = find(ismember(id,seqs{b}));
        is_in_clade{b} = clade_membership{ind};
    end
    u_cl = unique(is_in_clade);
    clade_label{a} = u_cl{1};
    if strcmp(clade{sc}, 'D') || strcmp(clade{sc}, 'G')
        % if cluster is not of clade, remove cluster
        if isempty(find(ismember(u_cl, clade{sc})))
            wa_clusters(a) = [];
            clade_label{a} = [];
        end
    end
end 
if sc ==1
    f = fopen('../results/cluster_assignment.tsv', 'w');
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
