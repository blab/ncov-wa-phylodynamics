threshold = 1;

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

exlude = {''};

date_cutoff = datenum('2020-12-31');

% distance based clustering
% read in iq tree from nextstrain pipeline
tree = phytreeread('../../ncov/results/wa_state/tree_raw.nwk');


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

leafs = get(tree, 'leafnames');

pariwise_distances = pdist(tree, 'Squareform', true)*29000;


%% get all leafnames


is_wa = find(ismember(division, 'Washington') & ~ismember(location, '')& ~ismember(location, 'NA'));
is_wa_leaf = false(length(leafs), 1);
for i = 1 : length(leafs)
    if sum(ismember(id(is_wa), leafs{i}))==1
        is_wa_leaf(i) = true;
        county{i} = location{find(ismember(id, leafs{i}))};
        date_val_wa(i) = date_val(find(ismember(id, leafs{i})));
    end
end

wa_leafs = leafs(is_wa_leaf);
wa_county = county(is_wa_leaf);
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

below_threshold = wa_distances < threshold;
below_threshold = triu(wa_distances < threshold);


%% count the number of pairs
use_pairs = [true, false];
f = fopen('../results/connections.tsv','w');
fprintf(f, 'from_county\tto_county');

fprintf(f, '\tfrom_lat\tto_lat\tfrom_long\tto_long');
fprintf(f, '\tprob\tpairs\tpairwise\n');


for j = 1 : 2
    if sum(ismember(wa_county,'Seattle'))>0
        wa_county{ismember(wa_county,'Seattle')} = 'King County';
    end
    unique_counties = unique(wa_county);
    pairs = zeros(length(unique_counties), length(unique_counties));

    wa_county_number = zeros(1,length(wa_county));
    for i = 1 : length(wa_county_number)
        wa_county_number(i) = find(ismember(unique_counties,wa_county(i)));
    end

    below_threshold_tmp = below_threshold;
    wa_county_number_tmp = wa_county_number;
    % only keep one individual from the same county per cluster
    remove = zeros(0,0);
    for i = 1 : size(below_threshold_tmp)
        related = find(below_threshold_tmp(i,:));
        ind = find(wa_county_number_tmp(related)==wa_county_number_tmp(i));
        if ~use_pairs(j)
            remove = [remove, related(ind)];
        end
    end

    remove = unique(remove);
    wa_county_number_tmp(remove) = [];
    below_threshold_tmp(:,remove) = [];
    below_threshold_tmp(remove,:) = [];


    [a,b] = find(below_threshold_tmp);
    for i = 1 : length(a)
        min_val = min(wa_county_number_tmp(a(i)), wa_county_number_tmp(b(i)));
        max_val = max(wa_county_number_tmp(a(i)), wa_county_number_tmp(b(i)));
        pairs(min_val,max_val) = pairs(min_val,max_val)+1;
    end

    nr_reps = 10000;

    prob_pairs_larger = zeros(size(pairs));
    prob_pairs_smaller = zeros(size(pairs));

    for r = 1 : nr_reps
        disp(r)
        wa_county_number_random = wa_county_number_tmp(randsample(length(wa_county_number_tmp),length(wa_county_number_tmp)));
        random_pairs = zeros(length(unique_counties), length(unique_counties));


        for i = 1 : length(a)
            min_val = min(wa_county_number_random(a(i)), wa_county_number_random(b(i)));
            max_val = max(wa_county_number_random(a(i)), wa_county_number_random(b(i)));
            random_pairs(min_val,max_val) = random_pairs(min_val,max_val)+1;
        end

        prob_pairs_larger = prob_pairs_larger+(pairs>random_pairs);  
        prob_pairs_smaller = prob_pairs_smaller+(pairs<random_pairs);  
    end

    for a = 1 : size(prob_pairs_larger,1)-1
        for b = a+1 : size(prob_pairs_larger,1)
            prob_pairs_larger(b,a)= prob_pairs_larger(a,b);
            prob_pairs_smaller(b,a)= prob_pairs_smaller(a,b);
        end
    end

    prob_pairs_larger = prob_pairs_larger/nr_reps;
    % prob_pairs_larger(prob_pairs_larger<0.95) = 0;

    prob_pairs_smaller = prob_pairs_smaller/nr_reps;
    prob_pairs_smaller(prob_pairs_smaller<0.95) = 0;

    for a = 1 : size(prob_pairs_larger,1)-1
        for b = a+1 : size(prob_pairs_larger,1)
            prob_pairs_larger(b,a) = prob_pairs_larger(a,b);
            prob_pairs_smaller(b,a) = prob_pairs_smaller(a,b);
        end
    end


    % tbl = array2table(prob_pairs_larger, 'VariableNames', unique_counties, 'RowNames',unique_counties);
    % heatmap(tbl)

    % read in lat and longs
    coordsdata = importdata('../../ncov/defaults/lat_longs.tsv');
    coords = zeros(0,2);
    for i = 1 : length(unique_counties)
        ind = find(ismember(coordsdata.textdata(:,2), unique_counties{i}));
        length(ind)
        coords(i, :) = coordsdata.data(ind,:);
    end

    seg_x = zeros(0,2);
    seg_y = zeros(0,2);
    width = zeros(0,1);
    c=1;
    for a = 1 : size(prob_pairs_larger,1)-1
        for b = a+1 : size(prob_pairs_larger,1)
            fprintf(f, '%s\t%s', unique_counties{a},unique_counties{b});
            fprintf(f, '\t%f\t%f\t%f\t%f', coords(a,1),coords(b,1), coords(a,2),coords(b,2));
            fprintf(f, '\t%f\t%f\t%d\n', prob_pairs_larger(a,b),pairs(a,b), use_pairs(j));
            
            if prob_pairs_larger(a,b)>0.95
                seg_x(c,1) = coords(a,2);
                seg_x(c,2) = coords(b,2);
                seg_y(c,1) = coords(a,1);
                seg_y(c,2) = coords(b,1);
                width(c,1) = prob_pairs_larger(a,b);
                pairs(c,1) = pairs(a,b);
                c = c+1;
            end
        end
    end
end
fclose(f);

% figure()
% for i = 1 : size(seg_x,1)
%     plot(seg_x(i,:),seg_y(i,:), 'LineWidth', width(i),'Color','k'); hold on
% end
% 
% figure()
% h = heatmap(prob_pairs_larger);
% h.XData = unique_counties;
% h.YData = unique_counties;



%% get the sampling times of pairs from different counties
if sum(ismember(wa_county,'Seattle'))>0
    wa_county{ismember(wa_county,'Seattle')} = 'King County';
end
unique_counties = unique(wa_county);
pairs = zeros(length(unique_counties), length(unique_counties));

wa_county_number = zeros(1,length(wa_county));
for i = 1 : length(wa_county_number)
    wa_county_number(i) = find(ismember(unique_counties,wa_county(i)));
end

below_threshold_tmp = below_threshold;
wa_county_number_tmp = wa_county_number;
% only keep one individual from the same county per cluster
remove = zeros(0,0);
relatives=cell(0,0);
for i = 1 : size(below_threshold_tmp)
    related = find(below_threshold_tmp(i,:));
    ind = find(wa_county_number_tmp(related)==wa_county_number_tmp(i));
    relatives{i,1} = [i related(ind)];
    remove = [remove, related(ind)];        
end


remove = unique(remove);
wa_county_number_tmp(remove) = [];
below_threshold_tmp(:,remove) = [];
below_threshold_tmp(remove,:) = [];
relatives(remove) = [];

% get the sampling times
samp_time = date_val_wa(is_wa_leaf);


[a,b] = find(below_threshold_tmp);
fopen('../results/pairs.tsv', 'w');
for i = 1 : length(a)
    fprintf(f, '%s',datestr(min(samp_time(relatives{a(i)})), 'yyyy-mm-dd'));
    fprintf(f, '\t%s',datestr(min(samp_time(relatives{b(i)})), 'yyyy-mm-dd'));
    fprintf(f, '\t%s',unique_counties{wa_county_number_tmp(a(i))});
    fprintf(f, '\t%s',unique_counties{wa_county_number_tmp(b(i))});
    fprintf(f, '\n');
end
fclose(f);
