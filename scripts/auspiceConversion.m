clear
% converts an auspice file to include, exlude different keys and change
% colors
counties = {'King County', 'Snohomish County', 'Island County', 'Whatcom County', 'Skagit County', 'San Juan County', 'Pierce County', 'Mason County', 'Kitsap County', 'Jefferson County WA', 'Clallam County', 'Thurston County', 'Lewis County', 'Pacific County', 'Grays Harbor County', 'Wahkiakum County', 'Cowlitz County', 'Skamania County', 'Clark County', 'Yakima County', 'Kittitas County', 'Benton County', 'Klickitat County', 'Franklin County WA', 'Walla Walla County', 'Adams County WA', 'Whitman County', 'Garfield County', 'Columbia County', 'Asotin County', 'Lincoln County', 'Spokane County', 'Ferry County', 'Stevens County', 'Pend Oreille County', 'Grant County WA', 'Okanogan County', 'Chelan County', 'Douglas County'};
county_colors ={'#5E1D9D','#4F1FAB','#492AB5','#4235C0','#4042C7','#3F50CC','#3E5DD0','#416ACF','#4377CD','#4682C9','#4A8CC2','#4F96BC','#549DB2','#5AA4A8','#60AA9E','#67AF93','#6EB389','#77B67F','#80B974','#89BB6B','#92BC63','#9CBE5B','#A6BE55','#B0BD4F','#B9BC4A','#C2BA46','#CBB742','#D2B240','#D9AD3D','#DEA63B','#E29D39','#E69337','#E68634','#E67932','#E56A2F','#E25A2C','#E04929','#DE3926','#DB2823'};
regions = {'Washington State', 'Europe', 'North America', 'Asia', 'Oceania', 'South America', 'Africa', 'Not WA'};
regions_colors = {'#268457', '#14309A', '#BF0B30', '#3CB3DF', '#030303', '#FBE000', '#F8A93D', '#A5A5A5'};

f = fopen('../results/node_locations.tsv');c=1;id=cell(0,0);loc=cell(0,0);
while ~feof(f)
    line = strsplit(strtrim(fgets(f)), '\t');
    id{c} = line{1};loc{c} = line{2};c=c+1;
end
fclose(f);

f = fopen('../../ncov-severity/data/WA_df.tsv');c=1;id_uw=cell(0,0);loc_uw=cell(0,0);
while ~feof(f)
    line =  strsplit(strtrim(fgets(f)), '\t');
    id_uw{c} = line{1};loc_uw{c} = line{5};c=c+1;
end
fclose(f);

f = fopen('../results/cluster_assignment.tsv');c=1;id_cl=cell(0,0);cluster=zeros(0,0);
while ~feof(f)
    line = strsplit(strtrim(fgets(f)), '\t');
    id_cl{c} = line{1};cluster(c) = str2double(line{2});c=c+1;
end
fclose(f);



f = fopen('../../ncov/auspice/ncov_wa_state.json');
g = fopen('../results/wa.json', 'w');

first_name = true;
lala=false;

while ~feof(f)
    line = fgets(f);    
    if contains(line, '"key": "recency",')% removes recency
        %%
        line = fgets(f);
        while ~contains(line, '{')
            line = fgets(f);
        end
    elseif contains(line, '"filters": [')
        fprintf(g,'%s',line);fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);
        fprintf(g,'      "region",\n');
        fprintf(g,'      "location",\n');
        fprintf(g,'      "cluster_size",\n');
        fprintf(g,'      "region",\n');
        fprintf(g,'      "author"\n');
   elseif contains (line, '"key": "country",')% removes country
        %%
        line = fgets(f);
        while ~contains(line, '{')
            line = fgets(f);
        end
    elseif contains (line, '"key": "location",')% only keeps WA counties
%         fprintf(g, '        "key": "location",\n');
%         fprintf(g, '        "title": "Location",\n');
%         fprintf(g, '        "type": "categorical"\n');
%         fprintf(g, '      },\n');
% 
%         line = fgets(f);
%         while ~contains(line, '{')
%             line = fgets(f);
%         end
%         fprintf(g, '%s',line);

        %%
        fprintf(g, '%s', line);line = fgets(f);
        fprintf(g, '%s', line);line = fgets(f);        
        first = true;
        while ~contains(line, 'title')            
            if contains(line, '",')
                tmp = strsplit(line, '"');
                if sum(ismember(counties,tmp{2}))==1
                    if first
                        first = false;
                    else
                        fprintf(g, '          ],\n');
                    end
                    fprintf(g,'          [\n');
                    fprintf(g,'            "%s",\n', tmp{2});line = fgets(f);     
                    fprintf(g,'            "%s"\n', county_colors{find(ismember(counties,tmp{2}))});
                    line = fgets(f);
                end
            end
            line = fgets(f);
        end
        fprintf(g, '          ]\n');
        fprintf(g, '	],\n');
        fprintf(g, '%s', line);line = fgets(f);
        fprintf(g, '%s', line);line = fgets(f);
        fprintf(g, '%s', line);
        fprintf(g, '      {\n');
        fprintf(g, '        "key": "cluster_size",\n');
        fprintf(g, '        "title": "Cluster size",\n');
        fprintf(g, '        "type": "categorical"\n');
        fprintf(g, '      },\n');

    elseif contains (line, '"key": "division",')% removes divison
        %%
        line = fgets(f);
        while ~contains(line, '{')
            line = fgets(f);
        end
    elseif contains(line, '"description": "Hi!')
        while ~contains(line, '},')
            line = fgets(f);
        end
        fprintf(g,'    "description": "This is a build for Washington State sequences alongside sequences from around the world",\n');
        fprintf(g,'    "display_defaults": {\n');
        fprintf(g,'      "branch_label": "clade",\n');
        fprintf(g,'      "color_by": "region",\n');
        fprintf(g,'      "distance_measure": "num_date",\n');
        fprintf(g,'      "geo_resolution": "location",\n');
        fprintf(g,'      "map_triplicate": true\n');
        fprintf(g,'    },\n');
        


    elseif contains (line, '"key": "region",')% removes and re-adds region
        %%
        line = fgets(f);
        while ~contains(line, '{')
            line = fgets(f);
        end
        
        fprintf(g,'        "key": "region",\n');
        fprintf(g,'        "scale": [\n');
        
        for i = 1 : length(regions)
            fprintf(g,'          [\n');
            fprintf(g,'            "%s",\n', regions{i});
            fprintf(g,'            "%s"\n',regions_colors{i});
            if i==length(regions)
                fprintf(g,'          ]\n');
            else
                fprintf(g,'          ],\n');
            end
        end
        fprintf(g,'        ],\n');
        fprintf(g,'        "title": "Region",\n');
        fprintf(g,'        "type": "categorical"\n');
        fprintf(g,'      },\n');
        fprintf(g,'      {\n');

    elseif contains(line, '"name": "USA/WA') && ~contains(line, 'travel')
        tmp = line;
        name_val = strsplit(tmp, '"');
        cl_ind = find(ismember(id_cl, name_val{4}));
        if isempty(cl_ind)
             fprintf(g, '%s', line);
        else
            has_loc = false;
            is_wa = true;
            while ~contains(line, '"region": {')
                if contains(line, 'location')
                    lines = cell(0,0);
                    lines{1} = line;
                    lines{2} = fgets(f);
                    lines{3} = fgets(f);
                    
                    tmp = strsplit(lines{2}, '"');
                    has_loc = true;
                    if sum(ismember(counties, tmp{4}))==1
                        fprintf(g, '%s', lines{1});
                        fprintf(g, '%s', lines{2});
                        fprintf(g, '%s', lines{3});
                    else                        
                        if strcmp(tmp{4}, 'Seattle')
                            fprintf(g, '%s', lines{1});
                            fprintf(g, '%s', strrep(lines{2}, 'Seattle','King County'));
                            fprintf(g, '%s', lines{3});
                        elseif strcmp(tmp{4}, 'Kirkland')
                            fprintf(g, '%s', lines{1});
                            fprintf(g, '%s', strrep(lines{2}, 'Seattle','King County'));
                            fprintf(g, '%s', lines{3});
                        else                            
                            is_wa = false;
                            has_loc = false;
                        end
                    end
                    has_loc = true;
                    line = fgets(f);
                else
                    fprintf(g, '%s', line);line = fgets(f);
                end
            end
            fprintf(g, '%s', line);line = fgets(f);
            if is_wa
                fprintf(g, '%s', strrep(line, 'North America', 'Washington State'));
            else
                fprintf(g, '%s', line);
            end
            line = fgets(f);

            fprintf(g, '%s', line);

            whitespace = regexp(line, '(\s*)', 'match');

            fprintf(g,'%s"cluster_size": {\n',whitespace{1});
            fprintf(g,'%s  "value": "%d"\n',whitespace{1}, sum(ismember(cluster, cluster(cl_ind))));
            fprintf(g,'%s},\n',whitespace{1});       
            if ~has_loc && is_wa
                ind = find(ismember(id_uw, name_val{4}));
                if ~isempty(ind) && ~strcmp(loc_uw{ind},'NA')
                    whitespace = regexp(line, '(\s*)', 'match');
                    fprintf(g, '%s"location": {\n',whitespace{1});
                    fprintf(g, '%s  "value": "%s County"\n',whitespace{1},loc_uw{ind});
                    fprintf(g, '%s },\n',whitespace{1});                           
                end
            end       
        end
    elseif contains(line, '"geo_resolutions": [')
        h = fopen('../data/geo_resolution.txt');
        while ~feof(h)
            line_h = fgets(h);
            fprintf(g, line_h);
        end
        while ~contains(line,'],')
            line = fgets(f);
        end
    elseif contains(line, '"name": "NODE')
        tmp = strsplit(line,'"');
        fprintf(g, '%s', line);line = fgets(f);
        fprintf(g, '%s', line);
        whitespace = regexp(line, '(\s*)', 'match');
        ind = find(ismember(id, tmp{4}));
        fprintf(g,'%s  "region": {\n',whitespace{1});
        fprintf(g,'%s  "value": "%s"\n',whitespace{1}, loc{ind});
        fprintf(g,'%s  },\n',whitespace{1});
    elseif contains(line, '"name": "')
        tmp = strsplit(line,'"');
        if first_name
            first_name = false;
            fprintf(g, '%s', strrep(line, 'the Nextstrain team', 'Nicola Felix Mueller'));line = fgets(f);
            fprintf(g, '%s', strrep(line, 'https://nextstrain.org/', 'https://bedford.io/team/nicola-mueller/'));
        else
            breakout=false;
            while ~contains(line, '"location": {')
                fprintf(g, '%s', line);line = fgets(f);
                if contains(line, '"url"')                    
                    break
                end
            end
            if contains(line, '"url"')
                fprintf(g, '%s', line);
            else
                line = fgets(f);line = fgets(f);
            end
        end
    else
        fprintf(g, '%s', line);
    end
end
fclose('all');
disp('done')



%%
f = fopen('../results/wa.json');
g = fopen('../results/ncov_wa-phylodynamics.json', 'w');

first_name = true;
lala=false;

while ~feof(f)
    line = fgets(f);
    if contains(line, 'location')
        fprintf(g, '%s', strrep(line, 'location','county'));
    elseif contains(line, 'Location')
        fprintf(g, '%s', strrep(line, 'Location','County'));
    else
        fprintf(g, '%s', line);
    end
end
fclose('all');
