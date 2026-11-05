function H = getRegularH(N,dv,dc)
% H = getRegularH(N,dv,dc)
%
% get parity-check matrix of a regular [dv,dc] LDPC code using the socket
% and permutation technique

% number of edges
E = N*dv;

% number of check nodes
M = E/dc;

if M ~= floor(M)
    error('Code construction with given size not possible, as it would lead to a non-integer number of check nodes');
end

% nodes connected to sockets
idx_j = ceil([1:E]/dc);
idx_i = ceil([1:E]/dv);

% random permutation
idx_j = idx_j(randperm(numel(idx_j)));


% try to construct matrix and eliminate double edges
abort = false;
while ~abort

    H = sparse(idx_j, idx_i, ones(size(idx_i)));
    % check for double edges
    if isempty(find(H==2))
        abort = true;        
    else             
        % eliminate double edges
        for k = 1:N
            ti = find(idx_i == k);
            if length(unique(idx_j(ti))) ~= length(idx_j(ti))
                [nu,~] = histc(idx_j(ti),1:M);
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
