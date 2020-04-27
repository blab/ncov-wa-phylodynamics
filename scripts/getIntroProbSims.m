clear

intro_perc = [0.1, 0.01, 0.001];
samples = round(logspace(0,log10(915),20));
f = fopen('../results/sim_cluster_size.tsv', 'w');
fprintf(f, 'nr_samples\tintro_prob\trun\tintro_percentage\n');

for i = 1 : length(intro_perc)
    disp(i)
    loc = 1;
    next_state = 1;
    for rep = 1 : 10
        while length(loc) < 10000
            if rand < intro_perc(i) % is intro
                next_state = next_state+1;
                loc(end+1) = next_state;
            else % is local
                new_loc = randsample(loc,1);
                loc(end+1) = new_loc;
            end    
        end

        for s = 1:length(samples)
            prob_curr = 0;
            nr_reps = 10000;
            for r = 1: nr_reps
                has_locs = randsample(loc, samples(s));
                has_locs_wo = has_locs;
                has_locs_wo(end) = [];
                if length(unique(has_locs))>length(unique(has_locs_wo))
                    prob_curr = prob_curr+1;
                end
            end
            fprintf(f, '%d\t%f\t%d\t%f\n', samples(s), prob_curr/nr_reps, rep, intro_perc(i));
        end
    end
end
fclose(f)
