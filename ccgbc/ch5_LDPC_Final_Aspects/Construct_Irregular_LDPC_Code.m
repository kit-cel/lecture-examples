function H = Construct_Irregular_LDPC_Code(variable_nodes, check_nodes)
% H = Construct_Irregular_LDPC_Code(variable_nodes, check_nodes)
%
% construct an irregular parity-check matrix using the socket method
% described in the lecture. The parameters are the following:
%
% variable_nodes: an nx2 array containing in each row the variable node
%                 degree and the number of variable node of that degree
% check_nodes:    similarly, an nx2 array containing in each row the check
%                 node degree and the number of check nodes of that degree 
%
% Example: H = Construct_Irregular_LDPC_Code([2 3; 3 3], [3, 5]) constructs
%          the code from the lecture with 3 VNs of degree 2, 3 VNs of
%          degree 3 and 3 CNs of degree 5
% ATTENTION: the code is for teaching purposes only and does *not* remove
% 4-cycles, this needs to be done in a separate step.

dv = variable_nodes(:,1);
n_dv = variable_nodes(:,2);

dc = check_nodes(:,1);
n_dc = check_nodes(:,2);

% number of variable nodes
N = sum(n_dv);

% number of edges
E = dv(:)'*n_dv(:);

% sanity check
if E ~= dc(:)'*n_dc(:)
    error('incompatible code definition ... number of edges do not match');
end

% construct variable node sockets
% idx_i tell for each socket, to which VN it is attached
idx_i = [];
offset = 0;
for k=1:numel(dv)
    temp = reshape(repmat([1:n_dv(k)],dv(k),1),[],1);
    idx_i = [idx_i; temp+offset];
    offset = offset + n_dv(k);
end

% construct check node sockets
% idx_j tell for each socket, to which CN it is attached
idx_j = [];
offset = 0;
for k=1:numel(dc)
    temp = reshape(repmat([1:n_dc(k)],dc(k),1),[],1);
    idx_j = [idx_j; temp+offset];
    offset = offset + n_dc(k);
end

% generate random permutation between sockets
permutation = randperm(numel(idx_j));

% permute check node sockets
idx_j = idx_j(permutation);

% eliminate parallel edges. A parallel edge will have an entry > 1 in the parity-check matrix due to the fact that the sparse function adds entries 
abort = false;
while ~abort
    % construct the parity-check matrix
    H = sparse(idx_j, idx_i, ones(size(idx_i)));
    % are there parallel edges?
    if isempty(find(H >= 2))
        abort = true;        
    else             
        % eliminate double edges, loop over every variable nodes
        for k = 1:N
            % find socket indices that connect to VN k
            ti = find(idx_i == k);
            % check if the connected CNs contain duplicates (unique)            
            if length(unique(idx_j(ti))) ~= length(idx_j(ti))
                % count connected check nodes
                [nu,~] = histc(idx_j(ti),1:sum(n_dc));
                % find largest one (multiple connections to that check node)
                largeval = find(nu > 1,1);
                largeidx = find(idx_j(ti)==largeval,1);
                
                % find a random edge
                rv = randi(numel(idx_j),1);
                temp = idx_j(ti(largeidx));
                idx_j(ti(largeidx)) = idx_j(rv);
                idx_j(rv) = temp;
            end
        end
    end
end