% main_poissonbeta.m
% PoissonBeta: a Poisson-Beta model for RNA-Seq data analysis
%
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: October 29 2011

clc; clear;
dir_data = 'C:\research\project\TangCellStemCell\data';
dir_poissonbeta = 'C:\research\research_program\PoissonBeta\source\PoissonBeta_Gene';

% read data
cd(dir_data);
[num, ~, raw] = xlsread('mmc3');
[num_row, num_column] = size(raw);

% select ESC
esc_index = zeros(num_column-1,1);
for ii=1:num_column
    if ~isempty(strfind(raw{1, ii}, 'ESC_A'))
        esc_index(ii-1) = 1;
    end;
end;

transcript_length = num(:, 6);
read_count = num(:, logical(esc_index));
read_count(:,13) = []; % remove a technical replicate

cd(dir_poissonbeta);
size_factor = PoissonBetaSizeFactor(read_count);

% select expressed genes

expressed_index = (sum(read_count,2)~=0);
read_count = read_count(expressed_index,:);
transcript_length = transcript_length(expressed_index);
gene_id = raw(2:end,1);
gene_id = gene_id(expressed_index);

cd(dir_poissonbeta);
Tmax = 1000;
clear num; clear raw;
[samples_Pij, samples_Si, samples_Koni, samples_Koffi, Qiter] = ...
    PoissonBeta(read_count, size_factor, transcript_length, Tmax);

cd(dir_data)
save mESC_genewise_T10000 samples_Pij samples_Si samples_Koni samples_Koffi Qiter ...
    read_count size_factor transcript_length gene_id;
