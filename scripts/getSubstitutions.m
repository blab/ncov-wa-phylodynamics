% 
clear
% read in the sequence data
fasta = fastaread('../../ncov/results/wa_state/sample-division.fasta');
% 
fasta_ref = fastaread('../../ncov/results/wa_state/sample-global.fasta');
for i=1:length(fasta_ref)
    if contains(fasta_ref(i).Header, 'Wuhan/WH01/2019')
        reference = fasta_ref(i).Sequence;
    end    
end
gb_ref = genbankread('../../ncov/defaults/reference_seq.gb');



f = fopen('../data/combined_meta.tsv');
line = strsplit(fgets(f), '\t');
date_id = find(ismember(line,'date'));
c=1;
id = cell(0,0);
date = cell(0,0);
while ~feof(f)
    line = strsplit(fgets(f), '\t','CollapseDelimiters', false);
    id{c,1} = line{1};
    date{c,1} = line{date_id};
    c=c+1;
end
fclose(f);

%
for j = 1 : length(gb_ref.CDS)
    protein{j} = nt2aa(reference(gb_ref.CDS(j).indices(1):gb_ref.CDS(j).indices(2)), 'ACGTOnly', false);
end

f = fopen('../results/subs_and_muts.tsv', 'w');
fprintf(f, 'Sequence\tDate');
for j = 1 : length(gb_ref.CDS)
    fprintf(f, '\t%s', gb_ref.CDS(j).gene);
end
fprintf(f, '\tAA\tMutations\n');

for i = 1 : length(fasta)
    ind = find(ismember(id,fasta(i).Header));
    if isempty(ind)
        date_str = 'NA';
    else
        date_str = date{ind};
    end

    fprintf(f,'%s\t%s', fasta(i).Header,date_str);

    AA = 0;
    % compare AA
    fasta(i).Sequence(fasta(i).Sequence=='-')='N';
    for j = 1 : length(gb_ref.CDS)
        AA_vals = nt2aa(fasta(i).Sequence(gb_ref.CDS(j).indices(1):gb_ref.CDS(j).indices(2)), 'ACGTOnly', false);
        diffs = protein{j}~=AA_vals;
        diffs(AA_vals=='X')=0;
        diffs(protein{j}=='X')=0;
        fprintf(f, '\t%d', sum(diffs));
        AA=AA+sum(diffs);
    end
    % compare RNA
    diffs = reference~=fasta(i).Sequence;
    diffs(reference=='N') = 0;
    diffs(fasta(i).Sequence=='N') = 0;
    fprintf(f,'\t%d\t%d\n', AA, sum(diffs));
end
fclose(f);



