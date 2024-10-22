% function [] = getWACluster()
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
max_nr_samples = [1500 500 1000 750];

% exclude WA1 
exlude = {'USA/WA1/2020'};
sample_cutoff = '2020-07-01';
clade = {'all', 'D', 'G', 'Yakima'};


end_date = '2020-01-31';


%% distance based clustering
date_cutoff = datenum(sample_cutoff);
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
    clade_found=false;
    if ~isempty(strain_id)
        clade_membership{c,1} = strain.strain{strain_id};
        clade_found = true;
    end
    date{c,1} = line{date_id};
    lab{c,1} = line{originating_id};
    sub{c,1} = line{submitting_id};


    if ~isempty(find(ismember(id{c,1}, exlude))) || ~clade_found
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

s = fopen('../results/mrsi.tsv', 'w');
cls = fopen('../results/cluster_size.tsv', 'w');
st = fopen('../results/sampling_times.tsv', 'w');

fprintf(s,'filename\tmrsi\tclade\n');
fprintf(cls,'filename\tnumber\tsize\tclade\n');
fprintf(st,'Date\tnumber\n');


for sc = 1 : length(clade)
    
% get all WA samples
    if strcmp(clade{sc}, 'Yakima')
        is_wa = find(ismember(location, 'Yakima County'));
    else
        is_wa = find(ismember(division, 'Washington') & ~ismember(location, 'Yakima County') & contains(location, 'County'));
    end
    
    % read in the cluster membership
	f = fopen('../results/cluster_assignment.tsv');fgets(f);c=1;
    cl_nr = zeros(0,0);
    while ~feof(f)
        tmp = strsplit(fgets(f));
        ind = find(ismember(id, tmp{1}));
        
        cl_nr(ind) = str2double(tmp{2});
    end
    fclose(f);
    
    
    
    wa_cl_nr = cl_nr(is_wa);
    cl_id = id(is_wa);
    uni_cl = unique(wa_cl_nr);
    uni_cl(uni_cl==0) = [];
    wa_clusters = cell(0,0);
    for i = 1 : length(uni_cl)
        ind = find(ismember(wa_cl_nr, uni_cl(i)));
        wa_clusters{i} = cl_id{ind(1)};
        if length(ind)>1
            for a = 2 : length(ind)
                wa_clusters{i} = [ wa_clusters{i} ',' cl_id{ind(a)}];
            end
        end
    end
        
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
        vals = 1;
        if length(u_cl)>1
            max_rep = sum(ismember(is_in_clade,u_cl{1}));            
            for u=2:length(u_cl)
                if sum(ismember(is_in_clade,u_cl{u}))>max_rep
                    vals=u;
                    max_rep=sum(ismember(is_in_clade,u_cl{u}));
                end
            end            
        end
        clade_label{a} = u_cl{vals};
        if strcmp(clade{sc}, 'D') || strcmp(clade{sc}, 'G')
            % if cluster is not of clade, remove cluster
            if isempty(find(ismember(u_cl, clade{sc})))
                wa_clusters(a) = [];
                clade_label{a} = [];
            end
        end
    end 

    %% randomly subsample was clusters   
    all_seqs = cell(0,0);
    wa_ind = zeros(0,0);
    for a = 1 : length(wa_clusters)
        seqs = strsplit(wa_clusters{a}, ',');
        for b = 1 : length(seqs)
            all_seqs{end+1,1} = seqs{b};
            wa_ind(end+1,1) = find(ismember(id,seqs{b}));
        end
    end
    
    % subsample based on sequence dates
    sample_date = date_val(wa_ind);
    use_sample = zeros(0,0);
    
    max_samples = min(length(wa_ind),max_nr_samples(sc));

    potential_samples = wa_ind;
    use_sample = randsample(potential_samples,max_samples);
%     while length(use_sample)<max_samples
%         % sample a date
%         uni_dates = unique(sample_date);
%         this_date = randsample(uni_dates,1);
%         if this_date ~= -1
%             potential_samples = wa_ind(sample_date==this_date);
%             if length(potential_samples)>1
%                 use_sample(end+1,1) = randsample(potential_samples,1);
%             else
%                 use_sample(end+1,1) = potential_samples;
%             end
%             sample_date(wa_ind==use_sample(end)) = -1;
%         end
%     end
    
    wa_clusters_curr = wa_clusters;
    clade_label_curr = clade_label;
    
    keep_samples = id(use_sample);
    for a = 1 : length(wa_clusters_curr)
        seqs = strsplit(wa_clusters_curr{a}, ',');
        for b = 1 : length(seqs)
            if sum(ismember(keep_samples, seqs{b}))==0
                wa_clusters_curr{a}=strrep(wa_clusters_curr{a}, seqs{b}, '');
                wa_clusters_curr{a}=strrep(wa_clusters_curr{a}, ',,', ',');
            end
        end
    end
    
    for a = length(wa_clusters_curr):-1:1
        if isempty(wa_clusters_curr{a}) || strcmp(wa_clusters_curr{a}, ',')
            wa_clusters_curr(a) = [];
            clade_label_curr(a) = [];
        else
            if strcmp(wa_clusters_curr{a}(1), ',')
                wa_clusters_curr{a}(1)=[];
            end
            if strcmp(wa_clusters_curr{a}(end), ',')
                wa_clusters_curr{a}(end)=[];
            end
        end
    end



    %% get the sampling times of each cluster
    sampling_times = cell(length(wa_clusters_curr),1);
    max_sampling_times = zeros(length(wa_clusters_curr),1);
    all_sampling = zeros(0,0);

    for a = 1 : length(sampling_times)
        sampling_times{a} = zeros(0,0);
        seqs = strsplit(wa_clusters_curr{a}, ',');
        for b = 1: length(seqs)
            % find the sequence index
            ind = find(ismember(id, seqs{b}));
            sampling_times{a}(b) = date_val(ind);
            all_sampling(end+1,1) =  date_val(ind);
        end
        max_sampling_times(a) = max(sampling_times{a});

    end
    disp(length(all_sampling))
    
    if sc==length(sample_cutoff)
        uni_sampling = unique(all_sampling);
        for us = 1 : length(uni_sampling)
            fprintf(st, '%s\t%d\n', datestr(uni_sampling(us), 'yyyy-mm-dd'), sum(all_sampling==uni_sampling(us)));
        end
    end
    
    disp(datestr(max(max_sampling_times)))
    
    rate_shifts = [3.5/366:3.5/366:(max(max_sampling_times)-datenum(end_date))/366 0.5 1];
    rate_shifts_immi = [7/366:7/366:(max(max_sampling_times)-datenum(end_date))/366 0.5 1];


    %% build the mutlti coal xml
    method = {'skygrid', 'independent', 'skygrowth', 'skysampling'};
    for sp = 1 : length(method)
        disp(sp)
        if sp == 2 || sp == 4
            continue
        end
        disp(sp)
        fprintf(s, 'multicoal_%s_%d\t%s\t%s\n',method{sp}, sc, datestr(max(max_sampling_times), 'yyyy-mm-dd'), clade{sc});
        for rep = 0 : 2
            f = fopen('../xml_templates/multicoal_template.xml');
            g = fopen(sprintf('../xmls/multicoal_%s_%d_rep%d.xml', method{sp}, sc, rep), 'w');

            % read in the sequence data
            fasta = fastaread('../../ncov/results/wa_state/sample-division.fasta');
            seq_id = cell(0,0);
            for i = 1 : length(fasta)
                seq_id{i} = fasta(i).Header;
            end

            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_data')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        fprintf(cls, 'multicoal_%s_%d\t%d\t%d\t%s\n',method{sp}, sc,a,length(seqs), clade_label_curr{a});                    
                        fprintf(g, '\t<data id="sequences_meta:lc_%d" spec="Alignment">\n',a);
                        for b = 1: length(seqs)
                            % find the sequence index
                            ind = find(ismember(seq_id, seqs{b}));
                            if isempty(ind)
                                error('sequence not found')
                            end
                            fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',seqs{b},seqs{b},...
                            fasta(ind).Sequence);    
                        end
                        fprintf(g, '\t</data>\n');
                    end
                elseif contains(line, 'insert_tree')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        date_vals = 'rem';
                        for b = 1 : length(seqs)
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
                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<parameter id="rootLength:lc_%d" name="stateNode" upper="1.0" dimension="1">0.1</parameter>\n',a);
                    end

                    fprintf(g,'\t\t\t<parameter id="rateShifts" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts));
                    fprintf(g,'\t\t\t<parameter id="rateShifts.immi" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts_immi));
                    if sp==4
                        fprintf(g,'\t\t\t<parameter id="samplingRate" name="stateNode">0</parameter>\n');
                        fprintf(g,'\t\t\t<parameter id="rateShifts.samp" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts_immi));
                    end

                elseif contains(line, 'insert_init_tree')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        
                        
                        fprintf(g, '\t\t<init spec="beast.util.ClusterTree" id="RandomTree.t:sequences_meta_%d" initial="@Tree.t:lc_%d" clusterType="upgma" taxa="@sequences_meta:lc_%d"/>\n', a,a,a);

                        
%                         fprintf(g, '\t\t<init id="RandomTree.t:sequences_meta_%d" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:lc_%d" taxa="@sequences_meta:lc_%d">\n', a,a,a);
%                         fprintf(g, '\t\t\t<populationModel id="ConstantPopulation0.%d" spec="ConstantPopulation">\n', a);
%                         fprintf(g, '\t\t\t\t<parameter id="randomPopSize.t:%d" spec="parameter.RealParameter" name="popSize">0.001</parameter>\n', a);
%                         fprintf(g, '\t\t\t</populationModel>\n');
%                         fprintf(g, '\t\t</init>\n');
                    end
                elseif contains(line, 'insert_priors')

                    fprintf(g,'\t\t\t\t<prior id="Sigmaprior2" name="distribution" x="@sigma.immi">\n');
                    fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.4" name="distr" meanInRealSpace="true" M="0.5" S="0.25"/>\n');
                    fprintf(g,'\t\t\t\t</prior>\n');

                    if sp==1
                        fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                        fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="0.1" S="0.5"/>\n');
                        fprintf(g,'\t\t\t\t</prior>\n');

                        fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@Ne">\n');
                        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                        fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                        fprintf(g,'\t\t\t\t\t</distr>\n');
    %                     fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="0" sigma="10"/>\n');
                        fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="-5" sigma="5"/>\n');
                        fprintf(g,'\t\t\t\t</distribution>\n');

                    elseif sp==2
                        fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@Ne">\n');
                        fprintf(g,'\t\t\t\t\t<Normal id="Uniform.3" name="distr" mean="0" sigma="2"/>\n');
                        fprintf(g,'\t\t\t\t</prior>\n');
                    elseif sp==3 || sp==4
                        fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                        fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="20" S="0.5"/>\n');
                        fprintf(g,'\t\t\t\t</prior>\n');

                        fprintf(g,'\t\t\t\t<distribution spec=''nab.skygrid.GrowthRateSmoothingPriorRealParam'' NeLog="@Ne" rateShifts="@rateShifts">\n');
                        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                        fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                        fprintf(g,'\t\t\t\t\t</distr>\n');
                        fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="-5" sigma="5"/>\n');
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
                    fprintf(g,'\t\t\t\t\t<immigrationRate id="timeVaryingMigrationRates" spec="nab.skygrid.TimeVaryingRates" rate="@immigrationRate" rateShifts="@rateShifts.immi"/>\n');
                    if sp==4
                        fprintf(g,'\t\t\t\t\t<samplingRate id="timeVaryingMigrationRates2" spec="nab.skygrid.TimeVaryingRates" rate="@samplingRate" rateShifts="@rateShifts.samp"/>\n');
                    end

                    fprintf(g,'\t\t\t\t\t<multiTreeIntervals id="TreeIntervals.t" spec="nab.multitree.MultiTreeIntervals">\n');
                    for a = 1 : length(wa_clusters_curr)
                        offset = (max(max_sampling_times)-max_sampling_times(a))/365;
                        fprintf(g,'\t\t\t\t\t\t<tree idref="Tree.t:lc_%d"/>\n', a);
                        fprintf(g,'\t\t\t\t\t\t<parameter id="offset:lc_%d" estimate="true" name="offset">%f</parameter>\n', a, offset);
                        fprintf(g,'\t\t\t\t\t\t<rootLength idref="rootLength:lc_%d"/>\n',a);
                    end
                    fprintf(g,'\t\t\t\t\t</multiTreeIntervals>\n');
                    fprintf(g,'\t\t\t\t</distribution>\n');
                elseif contains(line, 'insert_likelihood')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        if length(seqs)>1
                            fprintf(g,'\t\t\t\t<distribution id="treeLikelihood:lc_%d" spec="ThreadedTreeLikelihood" data="@sequences_meta:lc_%d" tree="@Tree.t:lc_%d" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n',a,a,a);
                        end
                    end

                elseif contains(line, 'insert_operators')
                    max_size = 0;
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        max_size = max([max_size, length(seqs)]);
                    end


                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        if length(seqs)>1
                            rel_weight=sqrt(length(seqs))/sqrt(max_size);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialTreeRootScaler.t:lc_%d" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialUniformOperator.t:lc_%d" spec="Uniform" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,30*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialSubtreeSlide.t:lc_%d" spec="SubtreeSlide" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialNarrow.t:lc_%d" spec="Exchange" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialWide.t:lc_%d" spec="Exchange" isNarrow="false" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialWilsonBalding.t:lc_%d" spec="WilsonBalding" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        end

                        fprintf(g,'\t\t<operator id="CoalescentRootLengthTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:lc_%d" weight="0.1"/>\n',a,a);

                    end
                    fprintf(g,'\t\t<operator id="AMVGoperator1" spec="AdaptableVarianceMultivariateNormalOperator" every="100" beta="0.1" scaleFactor="0.1" weight="10.0">\n');
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@Ne"/>\n');
                    if sp~=2
                        fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.Ne"/>\n');
                    end
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@immigrationRate"/>\n');
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.immi"/>\n');
                    if sp==4
                        fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@samplingRate"/>\n');
                    end

                    fprintf(g,'\t\t</operator>\n');
                    fprintf(g,'\t\t<operator id="RMW" spec="RealRandomWalkOperator" windowSize="0.5" parameter="@Ne" weight="5.0"/>\n');

                elseif contains(line, 'insert_logs')
                    fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t<log idref="sigma.immi"/>\n');

                    fprintf(g,'\t\t\t<log idref="Ne"/>\n');
                    fprintf(g,'\t\t\t<log idref="immigrationRate"/>\n');
                    if sp==4
                        fprintf(g,'\t\t\t<log idref="samplingRate"/>\n');
                    end

                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<log idref="rootLength:lc_%d"/>\n',a);
                    end

                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:lc_%d"/>\n',a,a);
                        fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger1:%d" spec="nab.util.MultiTreeStatLogger" heightOnly="true" tree="@Tree.t:lc_%d" offset="@offset:lc_%d" rootLength="@rootLength:lc_%d" />\n',a,a,a,a);
                        fprintf(g,'\t\t\t<log id="MultiTreeStatsLogger2:%d" spec="nab.util.MultiTreeStatLogger" originOnly="true" tree="@Tree.t:lc_%d" offset="@offset:lc_%d" rootLength="@rootLength:lc_%d" />\n',a,a,a,a);
                    end
                elseif contains(line, 'insert_logtree')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
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
    end
    disp('done')
    
    %% build the mutlti bdsky xml
    method = {'skygrid'};
    for sp = 1 : length(method)
        fprintf(s, 'multibd_%s_%d\t%s\t%s\n',method{sp}, sc, datestr(max(max_sampling_times), 'yyyy-mm-dd'), clade{sc});

        for rep = 0 : 2
            f = fopen('../xml_templates/multibd_template.xml');
            g = fopen(sprintf('../xmls/multibd_%s_%d_rep%d.xml', method{sp}, sc, rep), 'w');

            % read in the sequence data
            fasta = fastaread('../../ncov/results/wa_state/sample-division.fasta');
            seq_id = cell(0,0);
            for i = 1 : length(fasta)
                seq_id{i} = fasta(i).Header;
            end

            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_data')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        fprintf(cls, 'multibd_%s_%d\t%d\t%d\t%s\n',method{sp}, sc,a,length(seqs), clade_label_curr{a});                    
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
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
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
                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<parameter id="rootLength:lc_%d" name="stateNode" upper="0.1" dimension="1">0.01</parameter>\n',a);
                    end
                    fprintf(g,'\t\t\t<parameter id="logReproductiveNumber" dimension="%d" name="stateNode">0</parameter>\n', length(rate_shifts)+1);
                    fprintf(g,'\t\t\t<parameter id="samplingProportion" dimension="%d" name="stateNode">-2</parameter>\n',length(rate_shifts_immi)+1);
                    fprintf(g,'\t\t\t<parameter id="sigma.Sampling" dimension="1" name="stateNode">1000</parameter>\n');
                elseif contains(line, 'insert_init_tree')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        fprintf(g, '\t\t<init spec="beast.util.ClusterTree" id="RandomTree.t:sequences_meta_%d" initial="@Tree.t:lc_%d" clusterType="upgma" taxa="@sequences_meta:lc_%d"/>\n', a,a,a);

%                         fprintf(g, '\t\t<init id="RandomTree.t:sequences_meta_%d" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:lc_%d" taxa="@sequences_meta:lc_%d">\n', a,a,a);
%                         fprintf(g, '\t\t\t<populationModel id="ConstantPopulation0.%d" spec="ConstantPopulation">\n', a);
%                         fprintf(g, '\t\t\t\t<parameter id="randomPopSize.t:%d" spec="parameter.RealParameter" name="popSize">0.001</parameter>\n', a);
%                         fprintf(g, '\t\t\t</populationModel>\n');
%                         fprintf(g, '\t\t</init>\n');
                    end
                elseif contains(line, 'insert_priors')

                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t\t<prior id="RootLengthPrior:lc_%d" name="distribution" x="@rootLength:lc_%d">\n',a,a);
                        fprintf(g,'\t\t\t\t\t<LogNormal id="RootLengthPriorLogNormal:lc_%d" name="distr" meanInRealSpace="true" M="0.025" S="0.25"/>\n', a);
                        fprintf(g,'\t\t\t\t</prior>\n');
                    end


                    if sp==1
                        fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@sigma.Ne">\n');
                        fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="0.1" S="0.25"/>\n');
                        fprintf(g,'\t\t\t\t</prior>\n');

                        fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@logReproductiveNumber">\n');
                        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                        fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                        fprintf(g,'\t\t\t\t\t</distr>\n');
                        fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="1.3863" sigma="1"/>\n');
                        fprintf(g,'\t\t\t\t</distribution>\n');

    %                     fprintf(g,'\t\t\t\t<prior id="Sigmapriodsr1" name="distribution" x="@sigma.Sampling">\n');
    %                     fprintf(g,'\t\t\t\t\t<LogNormal id="Unidsform.3" name="distr" M="0" S="1"/>\n');
    %                     fprintf(g,'\t\t\t\t</prior>\n');

                        fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@samplingProportion">\n');
                        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                        fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Sampling"/>\n');
                        fprintf(g,'\t\t\t\t\t</distr>\n');
%                         fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Uniform" lower="-100000" upper="100000"/>\n');
                        fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="-3" sigma="1"/>\n');
%                         fprintf(g,'\t\t\t\t\t<finalDistr spec="beast.math.distributions.Normal"  mean="-3" sigma="1"/>\n');
                        fprintf(g,'\t\t\t\t</distribution>\n');


                    elseif sp==2
                        fprintf(g,'\t\t\t\t<prior id="Sigmaprior1" name="distribution" x="@Ne">\n');
                        fprintf(g,'\t\t\t\t\t<Normal id="Uniform.3" name="distr" mean="0" sigma="2"/>\n');
                        fprintf(g,'\t\t\t\t</prior>\n');
                    end

                    for a = 1 : length(wa_clusters_curr)
                        offset = (max(max_sampling_times)-max_sampling_times(a))/365;
                        
                        % compute rate shifts for R0
                        rate_shifts_offsetted = rate_shifts-offset;
                        offset_val = find(rate_shifts_offsetted>0);
                        offset_val = offset_val(1)-1;
                        rate_shifts_offsetted(rate_shifts_offsetted<=0)=[];
                        rate_shifts_offsetted = [0 rate_shifts_offsetted];
                        
                        % compute rate shifts for sampling
                        rate_shifts_offsetted_samp = rate_shifts_immi-offset;
                        offset_val_sampling = find(rate_shifts_offsetted_samp>0);
                        offset_val_sampling = offset_val_sampling(1)-1;
                        rate_shifts_offsetted_samp(rate_shifts_offsetted_samp<=0)=[];
                        rate_shifts_offsetted_samp = [0 rate_shifts_offsetted_samp];

                        fprintf(g,'\t\t\t\t<distribution id="BitrhDeathSkySerial.t:lc_%d" spec="beast.evolution.speciation.BirthDeathSkylineModel" origin="@rootLength:lc_%d" absoluteReproductiveNumber="@absoluteReproductiveNumber" tree="@Tree.t:lc_%d" originIsRootEdge="true">\n', a, a, a);
                        fprintf(g,'\t\t\t\t\t<reverseTimeArrays id="reverse.lc_%d" estimate="false" spec="parameter.BooleanParameter">true true true true</reverseTimeArrays>\n', a);
                        fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes1.lc_%d" estimate="false" name="birthRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', rate_shifts_offsetted));
                        fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes2.lc_%d" estimate="false" name="deathRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', 0));
                        fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes3.lc_%d" estimate="false" name="samplingRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', rate_shifts_offsetted_samp));
                        fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="logReproductiveNumber" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam1.lc_%d" parameter="@logReproductiveNumber" offset="%d"/>\n', a, offset_val);
                        fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="becomeUninfectiousRate" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam2.lc_%d" parameter="@becomeUninfectiousRate" offset="%d"/>\n', a, 0);
                        fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="samplingProportion" spec="beast.evolution.speciation.TruncatedRealParameter" isLog="true" id="truncatedParam3.lc_%d" parameter="@samplingProportion" offset="%d"/>\n', a, offset_val_sampling);
                        fprintf(g,'\t\t\t\t</distribution>\n');

                    end
                elseif contains(line, 'insert_likelihood')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        if length(seqs)>1
                            fprintf(g,'\t\t\t\t<distribution id="treeLikelihood:lc_%d" spec="ThreadedTreeLikelihood" data="@sequences_meta:lc_%d" tree="@Tree.t:lc_%d" siteModel="@SiteModel" branchRateModel="@ClockModel"/>\n',a,a,a);
                        end
                    end

                elseif contains(line, 'insert_operators')
                    max_size = 0;
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        max_size = max([max_size, length(seqs)]);
                    end


                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
                        if length(seqs)>1
                            rel_weight=sqrt(length(seqs))/sqrt(max_size);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialTreeRootScaler.t:lc_%d" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialUniformOperator.t:lc_%d" spec="Uniform" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,30*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialSubtreeSlide.t:lc_%d" spec="SubtreeSlide" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialNarrow.t:lc_%d" spec="Exchange" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a,15*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialWide.t:lc_%d" spec="Exchange" isNarrow="false" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                            fprintf(g, '\t\t<operator id="CoalescentExponentialWilsonBalding.t:lc_%d" spec="WilsonBalding" tree="@Tree.t:lc_%d" weight="%f"/>\n',a,a, 3*rel_weight);
                        end

                        fprintf(g,'\t\t<operator id="CoalescentRootLengthTreeScaler.t:lc_%d" spec="ScaleOperator" scaleFactor="0.5" parameter="@rootLength:lc_%d" weight="0.1"/>\n',a,a);

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

                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<log idref="rootLength:lc_%d"/>\n',a);
                    end

                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<log id="TreeStatsLogger:%d" spec="beast.evolution.tree.TreeStatLogger" tree="@Tree.t:lc_%d"/>\n',a,a);
                    end
                    for a = 1 : length(wa_clusters_curr)
                        fprintf(g,'\t\t\t<log idref="BitrhDeathSkySerial.t:lc_%d"/>\n',a);
                    end

                elseif contains(line, 'insert_logtree')
                    for a = 1 : length(wa_clusters_curr)
                        seqs = strsplit(wa_clusters_curr{a}, ',');
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
    end
    disp('done')
    

end
fclose('all');
