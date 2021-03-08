clear
%% distance based clustering
% get all WA sequences
f = fopen('../data/combined_meta.tsv');

line = strsplit(fgets(f), '\t');
div_id = find(ismember(line,'division'));
date_id = find(ismember(line,'date'));
location_id = find(ismember(line,'location'));

c=1;
id = cell(0,0);
clade_membership = cell(0,0);
division = cell(0,0);
sub = cell(0,0);
location = cell(0,0);
date = cell(0,0);
date_val = zeros(0,0);
while ~feof(f)
    line = strsplit(fgets(f), '\t','CollapseDelimiters', false);
    id{c,1} = line{1};
    date{c,1} = line{date_id};


    if ~contains(line{date_id}, 'X') && sum(date{c,1}=='-')==2
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



mask = [4221,8886,10554,13687,14222,14223,14225,17178,17179,17182,...
21149,21151,21209,21212,27658,27660,28184,29594];

%%
fasta = fastaread('../../ncov/results/wa_state/sample-division.fasta');
clear Sequence
seq_id = cell(0,0);
for i = 1 : length(fasta)
    Sequence(i,:) = fasta(i).Sequence(mask);
    seq_id{i} = fasta(i).Header;
end
%%
is_problem = zeros(0,0);
for i = 1 : size(Sequence,2)
    u = unique(Sequence(:,i));
    u(u=='N') = [];
    u(u=='-') = [];
%     disp(u')

    freqs = zeros(0,0);
    for j = 1 : length(u)
        freqs(j) = sum(Sequence(:,i)==u(j));
    end
%     disp(freqs)
    freqs(freqs==max(freqs))=-1;
    [~,max_ind]=max(freqs);
    disp(u(max_ind))
    disp(freqs(max_ind))
    is_problem = [is_problem find(Sequence(:,i)'==u(max_ind))];
end

uni = unique(is_problem);
%%

dates = zeros(0,0);
loc=cell(0,0);
emp = false(0,0);
for i = 1 : length(uni)
    ind = find(ismember(id, seq_id{uni(i)}));    
    dates(i) = date_val(ind); 
    loc{i} = location{ind};
    emp(i) = isempty(location{ind});
end

cutoff = datenum('2020-07-01'); 
fig=figure; 
hist(dates); hold on
line([cutoff cutoff])

sum(~emp)

