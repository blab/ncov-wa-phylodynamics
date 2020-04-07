% function [] = getWACluster()
threshold = 5;

% exclude WA1 
exlude = {'USA/WA1/2020'};
sample_cutoff = {'2020-03-10', '2020-03-17', '2020-03-24'};

end_date = '2020-01-10';

s = fopen('../config/mrsi.tsv', 'w');
cls = fopen('../config/cluster_size.tsv', 'w');

fprintf(s,'filename\tmrsi\n');
fprintf(cls,'filename\tnumber\tsize\n');

for sc = 1 : length(sample_cutoff)
    rate_shifts = [3/366:3/366:(datenum(sample_cutoff(sc))-datenum(end_date))/366];
    
    date_cutoff = datenum(sample_cutoff{sc});
    % read in iq tree from nextstrain pipeline
    tree = phytreeread('../../ncov/results/tree_raw.nwk');
    % get all WA sequences
    f = fopen('../../ncov/data/metadata.tsv');
    line = strsplit(fgets(f), '\t');
    div_id = find(ismember(line,'division'));
    date_id = find(ismember(line,'date'));
    c=1;
    id = cell(0,0);
    division = cell(0,0);
    date = cell(0,0);
    date_val = zeros(0,0);
    while ~feof(f)
        line = strsplit(fgets(f), '\t');
        id{c,1} = line{1};
        date{c,1} = line{date_id};
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
    % get all leafnames
    leafs = get(tree, 'leafnames');
    % get all WA samples
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
    pariwise_distances = pdist(tree, 'Squareform', true)*29000;
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


    %% get the sampling times of each cluster
    sampling_times = cell(length(wa_clusters),1);
    max_sampling_times = zeros(length(wa_clusters),1);
    for a = 1 : length(sampling_times)
        sampling_times{a} = zeros(0,0);
        seqs = strsplit(wa_clusters{a}, ',');
        for b = 1: length(seqs)
            % find the sequence index
            ind = find(ismember(id, seqs{b}));
            sampling_times{a}(b) = date_val(ind);
        end
        max_sampling_times(a) = max(sampling_times{a});

    end

    %% build the mutlti coal xml
    method = {'skygrid', 'independent', 'skygrowth'};
    for sp = 1 : length(method)
        f = fopen('../xml_templates/multicoal_template.xml');
        g = fopen(sprintf('../xmls/multicoal_%s_%d.xml', method{sp}, sc), 'w');
        fprintf(s, 'multicoal_%s_%d\t%s\n',method{sp}, sc, datestr(max(max_sampling_times), 'yyyy-mm-dd'));
       
        % read in the sequence data
        fasta = fastaread('../../ncov/results/masked.fasta');
        seq_id = cell(0,0);
        for i = 1 : length(fasta)
            seq_id{i} = fasta(i).Header;
        end

        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_data')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    fprintf(cls, 'multicoal_%s_%d\t%d\t%d\n',method{sp}, sc,a,length(seqs));                    
                    fprintf(g, '\t<data id="sequences_meta:lc_%d" spec="Alignment">\n',a);
                    for b = 1: length(seqs)
                        % find the sequence index
                        ind = find(ismember(seq_id, seqs{b}));
                        fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',seqs{b},seqs{b},...
                        fasta(ind).Sequence);    
                    end
                    fprintf(g, '\t</data>\n');
                end
            elseif contains(line, 'insert_tree')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    date_vals = 'rem';
                    for b = 1: length(seqs)
                        % find the sequence index
                        ind = find(ismember(seq_id, seqs{b}));
                        ind_date = find(ismember(id, seqs{b}));
                        date_vals = [date_vals ',' seqs{b} '=' date{ind_date}];
                    end
                    date_vals = strrep(date_vals, 'rem,','');

                    fprintf(g, '\t\t\t<tree id="Tree.t:lc_%d" spec="beast.evolution.tree.Tree" name="stateNode">\n', a);
                    fprintf(g, '\t\t\t\t<trait id="dateTrait.t:sequences_meta:lc_%d_1" spec="beast.evolution.tree.TraitSet" dateFormat="yyyy-M-dd" traitname="date" value="%s">\n', a, date_vals);
                    fprintf(g, '\t\t\t\t\t<taxa id="TaxonSet.sequences_meta:lc_%d_1" spec="TaxonSet">\n', a);
                    fprintf(g, '\t\t\t\t\t\t<alignment id="sequences_meta:lc_%d_1" spec="FilteredAlignment" filter="1::1">\n', a);
                    fprintf(g, '\t\t\t\t\t\t\t<data idref="sequences_meta:lc_%d"/>\n', a);
                    fprintf(g, '\t\t\t\t\t\t</alignment>\n');
                    fprintf(g, '\t\t\t\t\t</taxa>\n');
                    fprintf(g, '\t\t\t\t</trait>\n');
                    fprintf(g, '\t\t\t\t<taxonset idref="TaxonSet.sequences_meta:lc_%d_1"/>\n',a);
                    fprintf(g, '\t\t\t</tree>\n');
                end
            elseif contains(line, 'id="sigma.Ne" name="stateNode" dimension="1">') && sp==3
                fprintf(g, strrep(line,'>2<','>100<'));

            elseif contains(line, 'insert_rate_shifts')
                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<parameter id="rootLength:lc_%d" name="stateNode" upper="1.0" dimension="1">0.1</parameter>\n',a);
                end

                fprintf(g,'\t\t\t<parameter id="rateShifts" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts));
            elseif contains(line, 'insert_init_tree')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    fprintf(g, '\t\t<init id="RandomTree.t:sequences_meta_%d" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:lc_%d" taxa="@sequences_meta:lc_%d">\n', a,a,a);
                    fprintf(g, '\t\t\t<populationModel id="ConstantPopulation0.%d" spec="ConstantPopulation">\n', a);
                    fprintf(g, '\t\t\t\t<parameter id="randomPopSize.t:%d" spec="parameter.RealParameter" name="popSize">0.001</parameter>\n', a);
                    fprintf(g, '\t\t\t</populationModel>\n');
                    fprintf(g, '\t\t</init>\n');
                end
            elseif contains(line, 'insert_priors')

                fprintf(g,'\t\t\t\t<prior id="Sigmaprior2" name="distribution" x="@sigma.immi">\n');
                fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.4" name="distr" meanInRealSpace="true" M="0.5" S="0.25"/>\n');
                fprintf(g,'\t\t\t\t</prior>\n');

                if sp==1
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" M="0" S="1"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');

                    fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@Ne">\n');
                    fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                    fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t\t\t</distr>\n');
%                     fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="0" sigma="10"/>\n');
                    fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="-10" sigma="1"/>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');

                elseif sp==2
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@Ne">\n');
                    fprintf(g,'\t\t\t\t\t<Normal id="Uniform.3" name="distr" mean="0" sigma="2"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');
                elseif sp==3
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="100" S="0.5"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');

                    fprintf(g,'\t\t\t\t<distribution spec=''nab.skygrid.GrowthRateSmoothingPriorRealParam'' NeLog="@Ne" rateShifts="@rateShifts">\n');
                    fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                    fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t\t\t</distr>\n');
                    fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="-10" sigma="1"/>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');

                end
                fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@immigrationRate">\n');
                fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                fprintf(g,'\t\t\t\t\t<sigma idref="sigma.immi"/>\n');
                fprintf(g,'\t\t\t\t\t</distr>\n');
                fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="0" sigma="10"/>\n');
    %             fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="0" sigma="10"/>\n');
                fprintf(g,'\t\t\t\t</distribution>\n');
                fprintf(g,'\t\t\t\t<distribution id="CoalescentConstant.t" spec="nab.multitree.MultiTreeCoalescent" rateIsBackwards="true">\n');
                fprintf(g,'\t\t\t\t\t<populationModel id="Skygrid" spec="nab.skygrid.Skygrowth" logNe="@Ne" rateShifts="@rateShifts"/>\n');
                fprintf(g,'\t\t\t\t\t<immigrationRate id="timeVaryingMigrationRates" spec="nab.skygrid.TimeVaryingRates" rate="@immigrationRate" rateShifts="@rateShifts"/>\n');
                fprintf(g,'\t\t\t\t\t<multiTreeIntervals id="TreeIntervals.t" spec="nab.multitree.MultiTreeIntervals">\n');
                for a = 1 : length(wa_clusters)
                    offset = (max(max_sampling_times)-max_sampling_times(a))/365;
                    fprintf(g,'\t\t\t\t\t\t<tree idref="Tree.t:lc_%d"/>\n', a);
                    fprintf(g,'\t\t\t\t\t\t<parameter id="offset:lc_%d" estimate="true" name="offset">%f</parameter>\n', a, offset);
                    fprintf(g,'\t\t\t\t\t\t<rootLength idref="rootLength:lc_%d"/>\n',a);
                end
                fprintf(g,'\t\t\t\t\t</multiTreeIntervals>\n');
                fprintf(g,'\t\t\t\t</distribution>\n');
            elseif contains(line, 'insert_likelihood')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    if length(seqs)>1
                        fprintf(g,'\t\t\t\t<distribution id="treeLikelihood:lc_%d" spec="ThreadedTreeLikelihood" data="@sequences_meta:lc_%d" tree="@Tree.t:lc_%d" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n',a,a,a);
                    end
                end

            elseif contains(line, 'insert_operators')
                max_size = 0;
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    max_size = max([max_size, length(seqs)]);
                end


                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    if length(seqs)>1
                        rel_weight=length(seqs)/max_size;
                        fprintf(g, '\t\t<operator id="CoalescentExponentialTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialTreeRootScaler.t:lc_%d" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialUniformOperator.t:lc_%d" spec="Uniform" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,30*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialSubtreeSlide.t:lc_%d" spec="SubtreeSlide" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialNarrow.t:lc_%d" spec="Exchange" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialWide.t:lc_%d" spec="Exchange" isNarrow="false" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialWilsonBalding.t:lc_%d" spec="WilsonBalding" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                    end

                    fprintf(g,'\t\t<operator id="CoalescentRootLengthTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:lc_%d" weight="1.0"/>\n',a,a);

                end
                fprintf(g,'\t\t<operator id="AMVGoperator1" spec="AdaptableVarianceMultivariateNormalOperator" every="100" beta="0.1" scaleFactor="0.1" weight="10.0">\n');
                fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@Ne"/>\n');
                if sp==1
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.Ne"/>\n');
                end
                fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@immigrationRate"/>\n');
                fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.immi"/>\n');
                fprintf(g,'\t\t</operator>\n');
                fprintf(g,'\t\t<operator id="RMW" spec="RealRandomWalkOperator" windowSize="0.5" parameter="@Ne" weight="10.0"/>\n');

            elseif contains(line, 'insert_logs')
                fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="sigma.immi"/>\n');

                fprintf(g,'\t\t\t<log idref="Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="immigrationRate"/>\n');
                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<log idref="rootLength:lc_%d"/>\n',a);
                end

                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:lc_%d"/>\n',a,a);
                    fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger1:%d" spec="nab.util.MultiTreeStatLogger" heightOnly="true" tree="@Tree.t:lc_%d" offset="@offset:lc_%d" rootLength="@rootLength:lc_%d" />\n',a,a,a,a);
                    fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger2:%d" spec="nab.util.MultiTreeStatLogger" originOnly="true" tree="@Tree.t:lc_%d" offset="@offset:lc_%d" rootLength="@rootLength:lc_%d" />\n',a,a,a,a);
                end
            elseif contains(line, 'insert_logtree')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    if length(seqs)>1
                        fprintf(g,'\t\t<logger id="treelog.t:sequences_meta_%d" spec="Logger" fileName="$(filebase).lc_%d.trees" logEvery="100000" mode="tree">\n',a,a);
                        fprintf(g,'\t\t\t<log id="TreeWithMetaDataLogger.t:sequences_meta_%d" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:lc_%d"/>\n',a,a);
                        fprintf(g,'\t\t</logger>\n');
                    end
                end
            else
                fprintf(g, line);
            end
        end
        fclose(g);
        fclose(f);
    end
    disp('done')
    
    %% build the mutlti bdsky xml
    method = {'skygrid'};
    for sp = 1 : length(method)
        f = fopen('../xml_templates/multibd_template.xml');
        g = fopen(sprintf('../xmls/multibd_%s_%d.xml', method{sp}, sc), 'w');
        fprintf(s, 'multibd_%s_%d\t%s\n',method{sp}, sc, datestr(max(max_sampling_times), 'yyyy-mm-dd'));
       
        % read in the sequence data
        fasta = fastaread('../../ncov/results/masked.fasta');
        seq_id = cell(0,0);
        for i = 1 : length(fasta)
            seq_id{i} = fasta(i).Header;
        end

        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_data')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    fprintf(cls, 'multibd_%s_%d\t%d\t%d\n',method{sp}, sc,a,length(seqs));                    
                    fprintf(g, '\t<data id="sequences_meta:lc_%d" spec="Alignment">\n',a);
                    for b = 1: length(seqs)
                        % find the sequence index
                        ind = find(ismember(seq_id, seqs{b}));
                        fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',seqs{b},seqs{b},...
                        fasta(ind).Sequence);    
                    end
                    fprintf(g, '\t</data>\n');
                end
            elseif contains(line, 'insert_tree')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    date_vals = 'rem';
                    for b = 1: length(seqs)
                        % find the sequence index
                        ind = find(ismember(seq_id, seqs{b}));
                        ind_date = find(ismember(id, seqs{b}));
                        date_vals = [date_vals ',' seqs{b} '=' date{ind_date}];
                    end
                    date_vals = strrep(date_vals, 'rem,','');

                    fprintf(g, '\t\t\t<tree id="Tree.t:lc_%d" spec="beast.evolution.tree.Tree" name="stateNode">\n', a);
                    fprintf(g, '\t\t\t\t<trait id="dateTrait.t:sequences_meta:lc_%d_1" spec="beast.evolution.tree.TraitSet" dateFormat="yyyy-M-dd" traitname="date" value="%s">\n', a, date_vals);
                    fprintf(g, '\t\t\t\t\t<taxa id="TaxonSet.sequences_meta:lc_%d_1" spec="TaxonSet">\n', a);
                    fprintf(g, '\t\t\t\t\t\t<alignment id="sequences_meta:lc_%d_1" spec="FilteredAlignment" filter="1::1">\n', a);
                    fprintf(g, '\t\t\t\t\t\t\t<data idref="sequences_meta:lc_%d"/>\n', a);
                    fprintf(g, '\t\t\t\t\t\t</alignment>\n');
                    fprintf(g, '\t\t\t\t\t</taxa>\n');
                    fprintf(g, '\t\t\t\t</trait>\n');
                    fprintf(g, '\t\t\t\t<taxonset idref="TaxonSet.sequences_meta:lc_%d_1"/>\n',a);
                    fprintf(g, '\t\t\t</tree>\n');
                end
            elseif contains(line, 'id="sigma.Ne" name="stateNode" dimension="1">') && sp==3
                fprintf(g, strrep(line,'>2<','>100<'));
 
            elseif contains(line, 'insert_rate_shifts')
                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<parameter id="rootLength:lc_%d" name="stateNode" upper="0.1" dimension="1">0.01</parameter>\n',a);
                end
                fprintf(g,'\t\t\t<parameter id="logReproductiveNumber" dimension="%d" name="stateNode">0</parameter>\n', length(rate_shifts)+1);
                fprintf(g,'\t\t\t<parameter id="samplingProportion" dimension="%d" name="stateNode">-1</parameter>\n', length(rate_shifts)+1);
                fprintf(g,'\t\t\t<parameter id="sigma.Sampling" dimension="1" name="stateNode">1</parameter>\n');
            elseif contains(line, 'insert_init_tree')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    fprintf(g, '\t\t<init id="RandomTree.t:sequences_meta_%d" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:lc_%d" taxa="@sequences_meta:lc_%d">\n', a,a,a);
                    fprintf(g, '\t\t\t<populationModel id="ConstantPopulation0.%d" spec="ConstantPopulation">\n', a);
                    fprintf(g, '\t\t\t\t<parameter id="randomPopSize.t:%d" spec="parameter.RealParameter" name="popSize">0.001</parameter>\n', a);
                    fprintf(g, '\t\t\t</populationModel>\n');
                    fprintf(g, '\t\t</init>\n');
                end
            elseif contains(line, 'insert_priors')
                
                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t\t<prior id="RootLengthPrior:lc_%d" name="distribution" x="@rootLength:lc_%d">\n',a,a);
                    fprintf(g,'\t\t\t\t\t<LogNormal id="RootLengthPriorLogNormal:lc_%d" name="distr" meanInRealSpace="true" M="0.025" S="0.25"/>\n', a);
                    fprintf(g,'\t\t\t\t</prior>\n');
                end


                if sp==1
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" M="0" S="1"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');

                    fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@logReproductiveNumber">\n');
                    fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                    fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t\t\t</distr>\n');
                    fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="0" sigma="1"/>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');
                    
                    fprintf(g,'\t\t\t\t<prior id="Sigmapriodsr1" name="distribution" x="@sigma.Sampling">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Unidsform.3" name="distr" M="0" S="1"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');

                    fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@samplingProportion">\n');
                    fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                    fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Sampling"/>\n');
                    fprintf(g,'\t\t\t\t\t</distr>\n');
                    fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="-3" sigma="1"/>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');


                elseif sp==2
                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@Ne">\n');
                    fprintf(g,'\t\t\t\t\t<Normal id="Uniform.3" name="distr" mean="0" sigma="2"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');
                end
                
                for a = 1 : length(wa_clusters)
                    offset = (max(max_sampling_times)-max_sampling_times(a))/365;
                    rate_shifts_offsetted = rate_shifts-offset;
                    offset_val = find(rate_shifts_offsetted>0);
                    offset_val = offset_val(1)-1;
                    rate_shifts_offsetted(rate_shifts_offsetted<=0)=[];
                    rate_shifts_offsetted = [0 rate_shifts_offsetted];
                    fprintf(g,'\t\t\t\t<distribution id="BitrhDeathSkySerial.t:lc_%d" spec="beast.evolution.speciation.BirthDeathSkylineModel" origin="@rootLength:lc_%d" absoluteReproductiveNumber="@absoluteReproductiveNumber" tree="@Tree.t:lc_%d" originIsRootEdge="true">\n', a, a, a);
                    fprintf(g,'\t\t\t\t\t<reverseTimeArrays id="reverse.lc_%d" estimate="false" spec="parameter.BooleanParameter">true true true true</reverseTimeArrays>\n', a);
                    fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes1.lc_%d" estimate="false" name="birthRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', rate_shifts_offsetted));
                    fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes2.lc_%d" estimate="false" name="deathRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', 0));
                    fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes3.lc_%d" estimate="false" name="samplingRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', rate_shifts_offsetted));
                    fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="logReproductiveNumber" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam1.lc_%d" parameter="@logReproductiveNumber" offset="%d"/>\n', a, offset_val);
                    fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="becomeUninfectiousRate" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam2.lc_%d" parameter="@becomeUninfectiousRate" offset="%d"/>\n', a, 0);
                    fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="samplingProportion" spec="beast.evolution.speciation.TruncatedRealParameter" isLog="true" id="truncatedParam3.lc_%d" parameter="@samplingProportion" offset="%d"/>\n', a, offset_val);
                    fprintf(g,'\t\t\t\t</distribution>\n');
                
                end
            elseif contains(line, 'insert_likelihood')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    if length(seqs)>1
                        fprintf(g,'\t\t\t\t<distribution id="treeLikelihood:lc_%d" spec="ThreadedTreeLikelihood" data="@sequences_meta:lc_%d" tree="@Tree.t:lc_%d" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n',a,a,a);
                    end
                end

            elseif contains(line, 'insert_operators')
                max_size = 0;
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    max_size = max([max_size, length(seqs)]);
                end


                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    if length(seqs)>1
                        rel_weight=length(seqs)/max_size;
                        fprintf(g, '\t\t<operator id="CoalescentExponentialTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialTreeRootScaler.t:lc_%d" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialUniformOperator.t:lc_%d" spec="Uniform" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,30*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialSubtreeSlide.t:lc_%d" spec="SubtreeSlide" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialNarrow.t:lc_%d" spec="Exchange" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialWide.t:lc_%d" spec="Exchange" isNarrow="false" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        fprintf(g, '\t\t<operator id="CoalescentExponentialWilsonBalding.t:lc_%d" spec="WilsonBalding" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                    end

                    fprintf(g,'\t\t<operator id="CoalescentRootLengthTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:lc_%d" weight="1.0"/>\n',a,a);

                end
                fprintf(g,'\t\t<operator id="AMVGoperator1" spec="AdaptableVarianceMultivariateNormalOperator" every="100" beta="0.1" scaleFactor="0.1" weight="10.0">\n');
                fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@logReproductiveNumber"/>\n');
                if sp==1
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.Ne"/>\n');
                end
                fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.Sampling"/>\n');
                fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@samplingProportion"/>\n');
                fprintf(g,'\t\t</operator>\n');

            elseif contains(line, 'insert_logs')
                fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                fprintf(g,'\t\t\t<log idref="sigma.Sampling"/>\n');

                fprintf(g,'\t\t\t<log idref="logReproductiveNumber"/>\n');
                fprintf(g,'\t\t\t<log idref="samplingProportion"/>\n');
                
                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<log idref="rootLength:lc_%d"/>\n',a);
                end

                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:lc_%d"/>\n',a,a);
                end
                for a = 1 : length(wa_clusters)
                    fprintf(g,'\t\t\t<log idref="BitrhDeathSkySerial.t:lc_%d"/>\n',a);
                end

            elseif contains(line, 'insert_logtree')
                for a = 1 : length(wa_clusters)
                    seqs = strsplit(wa_clusters{a}, ',');
                    if length(seqs)>1
                        fprintf(g,'\t\t<logger id="treelog.t:sequences_meta_%d" spec="Logger" fileName="$(filebase).lc_%d.trees" logEvery="100000" mode="tree">\n',a,a);
                        fprintf(g,'\t\t\t<log id="TreeWithMetaDataLogger.t:sequences_meta_%d" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:lc_%d"/>\n',a,a);
                        fprintf(g,'\t\t</logger>\n');
                    end
                end
            else
                fprintf(g, line);
            end
        end
        fclose(g);
        fclose(f);
    end
    disp('done')
    

end
fclose('all')
