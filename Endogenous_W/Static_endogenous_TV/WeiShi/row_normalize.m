function W = row_normalize(W0)

n = length(W0(:,1));
W = zeros(n,n);

W_sum = zeros(n,1);
for i = 1:n
    for j = 1:n
        W_sum(i) = W_sum(i) + W0(i,j);
    end
end

for i = 1:n
    if W_sum(i) ~= 0
        for j = 1:n
            W(i,j) = W0(i,j)/W_sum(i);
        end
    end
end