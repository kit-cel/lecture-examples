function H = getIrregularH(N,lambda,rho)
% H = getRegularH(N,lambda,rho)
%
% lambda is the degree distribution polynomial from edge perspective and
% the first element is the coefficient corresponding to the largest degree,
% i.e., lambda = [lambda_N, lambda_{N-1}, ..., lambda_2, lambda_1]
%
% For example, consider
% lambda = [0.5 0 0 0.2 0.3 0]
% This would correspond to a polynomial
%                  5        2        
% Lambda(Z) = 0.5 Z  + 0.2 Z  + 0.3 Z
%
% rho is similarly defined
%
% For example, to generate the second irregular code example of the lecture
% with n = 10000, we can call
% H = getIrregularH(10000, [0.1151 0.1971 0 0 0.0768 0.202 0.409 0], [1 0 0 0 0 0])
% 
% we can verify that d_{c,avg} = 6 using dcavg = mean(sum(H,2)) and that
% d_{v,avg} = 3 using dvacg = mean(sum(H,1))

[t_dv,t_a_dv] = edge_to_node(lambda);
[t_dc,t_a_dc] = edge_to_node(rho);

[dv, n_dv, dc, n_dc] = degdist_to_integer(t_dv, t_a_dv, t_dc, t_a_dc, N);

n_dv = round(n_dv);
n_dc = round(n_dc);

E = dv(:)'*n_dv(:);


idx_i = [];
offset = 0;
for k=1:numel(dv)
    temp = reshape(repmat([1:n_dv(k)],dv(k),1),[],1);
    idx_i = [idx_i; temp+offset];
    offset = offset + n_dv(k);
end
idx_j = [];
offset = 0;
for k=1:numel(dc)
    temp = reshape(repmat([1:n_dc(k)],dc(k),1),[],1);
    idx_j = [idx_j; temp+offset];
    offset = offset + n_dc(k);
end


% find permutation
permutation = randperm(numel(idx_j));

idx_j = idx_j(permutation);

abort = false;
while ~abort

    H = sparse(idx_j, idx_i, ones(size(idx_i)));
    if isempty(find(H==2))
        abort = true;        
    else             
        % eliminate double edges
        for k = 1:N
            ti = find(idx_i == k);
            if length(unique(idx_j(ti))) ~= length(idx_j(ti))
                [nu,~] = histc(idx_j(ti),1:sum(n_dc));
                largeval = find(nu > 1,1);
                largeidx = find(idx_j(ti)==largeval,1);
                temp = idx_j(ti(largeidx));
                rv = randi(numel(idx_j),1);
                idx_j(ti(largeidx)) = idx_j(rv);
                idx_j(rv) = temp;
            end
        end
    end
end


end

function [dv,a] = edge_to_node(poly)
% [dv,a] = node_to_edge(poly)
%
% converts degree distribution from edge perspective as in density evolution for example
% to distribution from node perspective, 

int_pol = polyint(poly);
int_pol = int_pol ./ polyval(int_pol, 1);

idx = find(int_pol > 0);
a = int_pol(idx);
dv = length(int_pol) - idx;

a = a(end:-1:1);
dv = dv(end:-1:1);

end