function size_factor = PoissonBetaSizeFactor(read_count)
% PoissonBetaSizeFactor
% 
% by Jong Kyoung Kim
% Last update: 17 October 2011

[temp,n] = size(read_count);
nonzero_gene = (sum(log(read_count),2)~=-Inf);
normalized_read_count = read_count(nonzero_gene,:);
normalized_read_count = normalized_read_count./repmat(geomean(normalized_read_count,2), 1, n);
size_factor = zeros(1,n);
for ii=1:n
    size_factor(ii) = median(normalized_read_count(normalized_read_count(:,ii)~=0,ii));
end;
