% function [] = getLTT()
entry=1;
clear treeDat;

for nr =1:4
    tree_files = dir(sprintf('../out/multicoal_skygrowth_%d.*.trees', nr));

    max_day = 119;

    skip_nr = 1000;


    for i = 1 : length(tree_files)
        disp(tree_files(i).name)
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
        treeDat(entry).mean_vals = mean(all_vals);
        treeDat(entry).nr = nr;
        entry = entry+1;
    end
end

%%
% print to file
g = fopen('../results/ltt_mean.tsv','w');
c = 1;
fprintf(g, 'dataset\tcluster');
for i = 1 : length(treeDat(1).mean_vals)
     fprintf(g, '\t%d', i);
end
fprintf(g, '\n');

for nr = 1:4
    
    tree_files = dir(sprintf('../out/multicoal_skygrowth_%d.*.trees', nr));
    
    for i = 1 : length(tree_files)
        tmp = strsplit(tree_files(i).name, '.');
        fprintf(g, '%d\t%s',treeDat(c).nr, strrep(tmp{2}, 'lc_',''));
        for j = 1 : length(treeDat(c).mean_vals)
            fprintf(g, '\t%f', treeDat(c).mean_vals(j));
        end
        c = c+1;
        fprintf(g, '\n');

    end
end
fclose(g);

