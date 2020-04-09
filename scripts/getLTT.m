% function [] = getLTT()
tree_files = dir('../out/multicoal_skygrid_2*.trees');

max_day = 70;

skip_nr = 10;

clear treeDat;

for i = 1 : length(tree_files)
    f = fopen([tree_files(i).folder '/' tree_files(i).name]);
    trees = cell(0,0);
    c = 1;
    all_vals = zeros(1,max_day);
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'tree STATE')
            if c>=skip_nr
                tmp = strsplit(line);
                ptree = phytreeread(tmp{4});
                dist = pdist(ptree, 'SquareForm',true, 'Nodes', 'all');
                % only take distances relative to the root
                dist = dist(:,end);
                % convert to heights in days
                clear height
                height(:,1) = -(dist-max(dist))*366;
                % get the connectivity
                mat = getmatrix(ptree);
                [a,b] = find(mat);
                height(b,2) = height(a,1);
                vals = zeros(1,max_day);
                for j = 1 : length(vals)
                    vals(j) = sum(height(:,1)<=j & height(:,2)>j);
                end
                
                
                all_vals(c-skip_nr+1,:) = vals;
            end
            c = c+1;
        end
    end
    fclose(f);
    % get the mean values;
    treeDat(i).mean_vals = mean(all_vals);
end

%%
% print to file
g = fopen('../results/ltt_mean.tsv','w');
fprintf(g, 'day')
for i = 1 : length(tree_files)
    tmp = strsplit(tree_files(i).name, '.');
    fprintf(g, '\t%s', strrep(tmp{2}, 'lc_',''));
end
fprintf(g, '\n');

for j = 1 : length(treeDat(1).mean_vals)
    fprintf(g, '%d', j-1);
    for i = 1 : length(tree_files)
        fprintf(g, '\t%f', treeDat(i).mean_vals(j));
    end 
    fprintf(g, '\n');
end
fclose(f);

