% run sims
sample_cutoff = {'2020-03-24'};
end_date = '2020-01-25';

for rep = 0 : 8
    f = fopen('../xml_templates/simulation_template.xml');
    g = fopen(sprintf('../simulations/eir_%d.xml',rep),'w');
    while ~feof(f)
        line = fgets(f);
        if contains(line, 'insert_transmission')
            x = [1:200];
            k = 0.2;
            unin_rate = 5;
            for i = 1 : length(x)
                y(1) = nbinpdf(x(i),1,k)*unin_rate*5;
                y(2) = nbinpdf(x(i),1,k)*unin_rate*2;
                y(3) = nbinpdf(x(i),1,k)*unin_rate;

                fprintf(g, '\t\t<reactionGroup spec=''ReactionGroup'' reactionGroupName=''Infection''>\n');
                fprintf(g, '\t\t\t<reaction spec=''Reaction'' rate="%f,%f:0.0405,%f:0.0787 ">\n',y(1),y(2),y(3));
                fprintf(g, '\t\t\t\tI:1 -> I:1 + %dE:1\n',x(i));
                fprintf(g, '\t\t\t</reaction>\n');
                fprintf(g, '\t\t</reactionGroup>\n');
                fprintf(g, '\t\t<reactionGroup spec=''ReactionGroup'' reactionGroupName=''Infection''>\n');
                fprintf(g, '\t\t\t<reaction spec=''Reaction'' rate="%f,%f:0.0405,%f:0.0787 ">\n',y(1),y(2),y(3));
                fprintf(g, '\t\t\t\tE:1 -> %dE:1\n',x(i)+1);
                fprintf(g, '\t\t\t</reaction>\n');
                fprintf(g, '\t\t</reactionGroup>\n');
            end
        else
            fprintf(g, line);
        end
    end
    fclose('all');
    system(sprintf('/Applications/BEAST\\ 2.6.0/bin/beast -overwrite ../simulations/eir_%d.xml', rep));
    
    origin = datenum(sample_cutoff{1})-0.117*365;

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
        fprintf('%s %s\n',tmp{2}, new_length)
        tree = strrep(tree, name{i}, new_length);

    end

    % rename tips
    tree = regexprep(tree, '((\d*):','(i$1:');
    tree = regexprep(tree, ',(\d*):',',i$1:');

    f = fopen(sprintf('../simulations/eir_%d.nexus', rep));fgets(f);fgets(f);fgets(f);
    line = fgets(f);
    % get all the sampling times
    st = regexp(line, '(\d*)\[&type="I",reaction="Sampling",time=(\d*)\.(\d*)\]','match');

    m = fopen(sprintf('../simulations/eir_%d.tsv', rep),'w');
    for i = 1 : length(st)
        tmp = strsplit(st{i}, '[');
        tmp2 = strsplit(st{i}, '=');
        day = origin+floor((str2double(strrep(tmp2{end}, ']','')))*366);
        day = datestr(day+reporting_delay, 'yyyy-mm-dd');
        fprintf(m,'i%s\t%s\n',tmp{1}, day);
    end

    ptree = phytreeread(tree);
    fprintf(g, getnewickstr(ptree));
    fclose('all');

    delete(sprintf('../simulations/eir_%d.fasta', rep));

    command = sprintf('%s -mHKY -t6.5 -a 0.05 %s -l 29000 -s %.12f < ../simulations/eir_%d.conv.tree > %s',...
    '../Software/seq-gen1.3.3',...
        '-f0.3,0.2,0.2,0.3',0.0008, rep, sprintf('../simulations/eir_%d.fasta', rep));
    system(command);
end