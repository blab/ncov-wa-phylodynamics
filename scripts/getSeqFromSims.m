% getSeqFromSims

% defines the start of the simulation
sample_cutoff = {'2020-06-31'};
end_date = '2020-01-25';

% define the reporting delay in days
reporting_delay = 0;

for rep = 0:8
    for sc = 1 : length(sample_cutoff)
        rate_shifts = [2/366:2/366:(datenum(sample_cutoff(sc))-datenum(end_date))/366 0.5];
        rate_shifts_immi = [14/366:14/366:(datenum(sample_cutoff(sc))-datenum(end_date))/366 0.5];

        f = fopen(sprintf('../simulations/eir_%d.tsv', rep));
        c=1;
        id = cell(0,0);
        date = cell(0,0);
        date_val = zeros(0,0);
        wa_clusters = cell(1000,1);
        
        while ~feof(f)
            line = strsplit(fgets(f), '\t');
            date_num = datenum(line{2});
            if date_num<=datenum(sample_cutoff(sc))
                id{c,1} = line{1};
                date{c,1} = line{2};
                
                wa_clusters{str2double(line{3})} = [wa_clusters{str2double(line{3})} ',' id{c,1}];
                c=c+1;
            end
        end
        wa_clusters = wa_clusters(~cellfun('isempty',wa_clusters));
        for i = 1:length(wa_clusters)
            wa_clusters{i} = wa_clusters{i}(2:end);
        end
        fclose(f);
        
        %% get the sampling times of each cluster
        sampling_times = cell(length(wa_clusters),1);
        max_sampling_times = zeros(length(wa_clusters),1);
        all_sampling = zeros(0,0);

        for a = 1 : length(sampling_times)
            sampling_times{a} = zeros(0,0);
            seqs = strsplit(wa_clusters{a}, ',');
            for b = 1: length(seqs)
                % find the sequence index
                ind = find(ismember(id, seqs{b}));
                sampling_times{a}(b) = datenum(date(ind));
                all_sampling(end+1,1) =  datenum(date(ind));
            end
            max_sampling_times(a) = max(sampling_times{a});
        end



        %% build the mutlti coal xml
        method = {'skygrid', 'independent', 'skygrowth', 'skysampling'};
        for sp = 1 : length(method)
            if sp == 2 || sp == 4
                continue;
            end

            f = fopen('../xml_templates/multicoal_template.xml');
            g = fopen(sprintf('../xmls/simmulticoal_%s_%d.xml', method{sp}, rep), 'w');

            s = fopen(sprintf('../simulations/eir_%d.fasta', rep));
            fgets(s);c=1;seq_id=cell(0,0);
            while ~feof(s)
                line = strsplit(strtrim(fgets(s)));
                fasta(c).Header = line{1};
                seq_id{c} = line{1};
                fasta(c).Sequence = line{2};
                c = c+1;
            end
            fclose(s);


            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_data')
                    for a = 1 : length(wa_clusters)
                        seqs = strsplit(wa_clusters{a}, ',');
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
                    fprintf(g,'\t\t\t<parameter id="rateShifts.immi" name="stateNode">%s</parameter>\n', sprintf('%f ', rate_shifts_immi));
                    if sp==4
                        fprintf(g,'\t\t\t<parameter id="samplingRate" name="stateNode">0</parameter>\n');
                        fprintf(g,'\t\t\t<parameter id="rateShifts.samp" name="stateNode">0.5</parameter>\n');
                    end
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
                    if sp~=2
                        fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.Ne"/>\n');
                    end
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@immigrationRate"/>\n');
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$LogTransform" f="@sigma.immi"/>\n');
                    fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@immigrationRate"/>\n');
                    if sp==4
                        fprintf(g,'\t\t\t<transformations spec="beast.util.Transform$NoTransform" f="@samplingRate"/>\n');
                    end
                    fprintf(g,'\t\t</operator>\n');
                    fprintf(g,'\t\t<operator id="RMW" spec="RealRandomWalkOperator" windowSize="0.5" parameter="@Ne" weight="10.0"/>\n');

                elseif contains(line, 'insert_logs')
                    fprintf(g,'\t\t\t<log idref="sigma.Ne"/>\n');
                    fprintf(g,'\t\t\t<log idref="sigma.immi"/>\n');

                    fprintf(g,'\t\t\t<log idref="Ne"/>\n');
                    fprintf(g,'\t\t\t<log idref="immigrationRate"/>\n');
                    if sp==4
                        fprintf(g,'\t\t\t<log idref="samplingRate"/>\n');
                    end
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
            g = fopen(sprintf('../xmls/simmultibd_%s_%d.xml', method{sp}, rep), 'w');


            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_data')
                    for a = 1 : length(wa_clusters)
                        seqs = strsplit(wa_clusters{a}, ',');
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
                    fprintf(g,'\t\t\t<parameter id="samplingProportion" dimension="%d" name="stateNode">-1</parameter>\n', 1);
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
                        fprintf(g,'\t\t\t\t\t<LogNormal id="Uniform.3" name="distr" meanInRealSpace="true" M="0.1" S="0.25"/>\n');
                        fprintf(g,'\t\t\t\t</prior>\n');

                        fprintf(g,'\t\t\t\t<distribution spec=''beast.mascotskyline.skyline.LogSmoothingPrior'' NeLog="@logReproductiveNumber">\n');
                        fprintf(g,'\t\t\t\t\t<distr spec="beast.math.distributions.Normal"  mean="0">\n');
                        fprintf(g,'\t\t\t\t\t<sigma idref="sigma.Ne"/>\n');
                        fprintf(g,'\t\t\t\t\t</distr>\n');
                        fprintf(g,'\t\t\t\t\t<initialDistr spec="beast.math.distributions.Normal"  mean="1.3863" sigma="1"/>\n');
                        fprintf(g,'\t\t\t\t</distribution>\n');
                    
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
                        fprintf(g,'\t\t\t\t\t<parameter id="intervalTimes3.lc_%d" estimate="false" name="samplingRateChangeTimes">%s</parameter>\n', a, sprintf('%f ', 0));
                        fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="logReproductiveNumber" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam1.lc_%d" parameter="@logReproductiveNumber" offset="%d"/>\n', a, offset_val);
                        fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="becomeUninfectiousRate" spec="beast.evolution.speciation.TruncatedRealParameter" id="truncatedParam2.lc_%d" parameter="@becomeUninfectiousRate" offset="%d"/>\n', a, 0);
                        fprintf(g,'\t\t\t\t\t<truncatedRealParameter name="samplingProportion" spec="beast.evolution.speciation.TruncatedRealParameter" isLog="true" id="truncatedParam3.lc_%d" parameter="@samplingProportion" offset="%d"/>\n', a, 0);
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
end
