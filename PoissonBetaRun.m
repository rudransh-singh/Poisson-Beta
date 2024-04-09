function fid_out = PoissonBetaRun(input_file, output_file, Tmax)

% PoissonBeta
% 
% Input: 
% 1. input_file: input file for count data
% 2. output_file: output file for results
% 3. Tmax: a maximum number of Gibbs sampling interations
% 
% By Jong Kyoung Kim 
% jkkim@ebi.ac.uk
% Last Update: 26 December 2012

read_count = importdata(input_file, '\t', 1);
transcript_length = read_count.data(:,1);
gene_id = read_count.textdata(2:end,1);
read_count = read_count.data(:, 2:end);

if ischar(Tmax)
    Tmax = str2double(Tmax);
end;

[num_gene, num_cell] = size(read_count);

expressed_index = (sum(read_count,2)~=0);
read_count_expressed = read_count(expressed_index,:);

size_factor = PoissonBetaSizeFactor(read_count_expressed);

[samples_Pij, samples_Si, samples_Koni, samples_Koffi, Qiter] = ...
    PoissonBeta(read_count_expressed, size_factor, transcript_length(expressed_index), Tmax);

mean_Si = mean(samples_Si, 2);
mean_Koni = mean(samples_Koni, 2);
mean_Koffi = mean(samples_Koffi, 2);
mean_Pij = 1-samples_Pij;
mean_SKoffi = mean(samples_Si./samples_Koffi,2);
mean_Exi = mean(samples_Si.*(samples_Koni./(samples_Koffi+samples_Koni)),2);

fid_out = fopen(output_file, 'w');
fprintf(fid_out, '%s\t%s\t%s\t%s\t%s\t%s\t', 'Gene', 'Si', 'Koni', 'Koffi', 'SKoffi', 'Exi');
for jj=1:num_cell
    fprintf(fid_out, '%s\t', strcat('Pij', '(Cell', num2str(jj),')'));
end;
fprintf(fid_out, '\n');
count = 1;
for ii=1:num_gene
    fprintf(fid_out, '%s\t', gene_id{ii});
    if expressed_index(ii) == 1
        fprintf(fid_out, '%f\t%f\t%f\t%f\t%f\t', mean_Si(count), mean_Koni(count), ...
            mean_Koffi(count), mean_SKoffi(count), mean_Exi(count));
        for jj=1:num_cell
            fprintf(fid_out, '%f\t', mean_Pij(count,jj));
        end;
        count = count + 1;
    else
        fprintf(fid_out, '%f\t%f\t%f\t%f\t%f\t', 0, 0, ...
            0, 0, 0);
        for jj=1:num_cell
            fprintf(fid_out, '%f\t', 0);
        end;  
    end;
    fprintf(fid_out, '\n');
end;
fclose(fid_out);