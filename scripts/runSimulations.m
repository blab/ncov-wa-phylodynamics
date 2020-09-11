% run sims
sample_cutoff = {'2020-06-31'};
end_date = '2020-01-25';

start_tree = '2020-02-06';

clock_rate = 0.0011;

k = 1;

syms p r;
reporting_delay=0;

rate_shift = datenum(start_tree, 'yyyy-mm-dd'):2:(datenum(sample_cutoff, 'yyyy-mm-dd')+14);
rate_shift_year = (rate_shift - datenum(start_tree, 'yyyy-mm-dd'))/366;

R0 = [repmat(2.5,1,8) 2.5:-(1.4/9):0.6 repmat(0.6,1,15) 0.6:0.05:1];
R0(end+1:length(rate_shift_year)) = R0(end);

eqn1 = p*r/(1-p);
eqn2 = p*r/(1-p)^2;


for i = 1 : length(R0)
    sol = solve([eqn1==R0(i), eqn2 == R0(i)+R0(i)^2/k], [p,r]);
    rvalue(i) = double(sol.r);
    pvalue(i) = 1-double(sol.p);
end
x = [1:100];

unin_rate = 52.2857;


y = zeros(length(x), length(R0));
for a = 1 : size(y,1)
    for b = 1 : size(y,2)
        y(a,b) = nbinpdf(x(a),rvalue(b),pvalue(b))*unin_rate;
    end
end


for rep = 0 : 8
    f = fopen('../xml_templates/simulation_template.xml');
    g = fopen(sprintf('../simulations/eir_%d.xml',rep),'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_transmission')
%             rate_string = sprintf('%f', R0(1)*unin_rate);
%             for b = 2 : length(R0)
%                 rate_string = sprintf('%s,%f:%f', rate_string, R0(b)*unin_rate,rate_shift_year(b));
%             end
%             fprintf(g, '\t\t<reactionGroup spec=''ReactionGroup'' reactionGroupName=''Infection''>\n');
%             fprintf(g, '\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
%             fprintf(g, '\t\t\t\tI:1 -> %dI:1\n',2);
%             fprintf(g, '\t\t\t</reaction>\n');
%             fprintf(g, '\t\t</reactionGroup>\n');

            for i = 1 : length(x)
                rate_string = sprintf('%f', y(i,1));
                for b = 2 : length(R0)
                    rate_string = sprintf('%s,%f:%f', rate_string, y(i,b),rate_shift_year(b-1));
                end
                fprintf(g, '\t\t<reactionGroup spec=''ReactionGroup'' reactionGroupName=''Infection''>\n');
                fprintf(g, '\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
%                 fprintf(g, '\t\t\t\tI:1 -> I:1 + %dI:1\n',x(i));
                fprintf(g, '\t\t\t\tI:1 -> %dI:1\n',x(i)+1);
                fprintf(g, '\t\t\t</reaction>\n');
                fprintf(g, '\t\t</reactionGroup>\n');
%                 fprintf(g, '\t\t<reactionGroup spec=''ReactionGroup'' reactionGroupName=''Infection''>\n');
%                 fprintf(g, '\t\t\t<reaction spec=''Reaction'' rate="%s">\n',rate_string);
%                 fprintf(g, '\t\t\t\tE:1 -> %dE:1\n',x(i)+1);
%                 fprintf(g, '\t\t\t</reaction>\n');
%                 fprintf(g, '\t\t</reactionGroup>\n');
            end
        elseif contains(line, 'insert_simulation_time')
             fprintf(g, strrep(line, 'insert_simulation_time', num2str(max(rate_shift_year))));
        else
            fprintf(g, line);
        end
    end
    fclose('all');
    system(sprintf('/Applications/BEAST\\ 2.6.0/bin/beast -overwrite ../simulations/eir_%d.xml', rep));
    
    origin = datenum(start_tree);

    t = fopen(sprintf('../simulations/eir_%d.tree', rep));

    tree = fgets(t);fclose(t);
    g = fopen(sprintf('../simulations/eir_%d.conv.tree', rep),'w');
    sampling_times = zeros(0,0);

    % add the reporting delay on them

    name = regexp(tree, '\((\d*)\:(\d*).(\d*)?(E-(\d*)|)','match');
    name = [name, regexp(tree, ',(\d*):(\d*).(\d*)?(E-(\d*)|)','match')];
    for i = 1 : length(name)
        tmp = strsplit(name{i}, ':');
        new_length = [tmp{1} ':' sprintf('%.12f',(str2double(tmp{2})+reporting_delay/366))];
%         fprintf('%s %s\n',tmp{2}, new_length)
        tree = strrep(tree, name{i}, new_length);

    end

    % rename tips
    tree = regexprep(tree, '((\d*):','(i$1:');
    tree = regexprep(tree, ',(\d*):',',i$1:');

    f = fopen(sprintf('../simulations/eir_%d.nexus', rep));fgets(f);fgets(f);fgets(f);
    line = fgets(f);
    % get all the sampling times
    st = regexp(line, '(\d*)\[&type="I",reaction="Sampling",time=(\d*)\.(\d*)\]','match');
    % get the cluster membership of leafs
    tree_str = strsplit(strtrim(line));
    cluster_member = getClustersFromTree(tree_str{end});
    
    
    m = fopen(sprintf('../simulations/eir_%d.tsv', rep),'w');
    for i = 1 : length(st)
        tmp = strsplit(st{i}, '[');
        tmp2 = strsplit(st{i}, '=');
        day = origin+floor((str2double(strrep(tmp2{end}, ']','')))*366);
        day = datestr(day+reporting_delay, 'yyyy-mm-dd');
        fprintf(m,'i%s\t%s\t%d\n',tmp{1}, day, cluster_member(str2double(tmp{1})));
    end

    ptree = phytreeread(tree);
    fprintf(g, getnewickstr(ptree));
    fclose('all');

    delete(sprintf('../simulations/eir_%d.fasta', rep));

    command = sprintf('%s -mHKY -t6.5 -a 0.05 %s -l 29000 -s %.12f < ../simulations/eir_%d.conv.tree > %s',...
    '../Software/seq-gen1.3.3',...
        '-f0.3,0.2,0.2,0.3',clock_rate, rep, sprintf('../simulations/eir_%d.fasta', rep));
    system(command);
end