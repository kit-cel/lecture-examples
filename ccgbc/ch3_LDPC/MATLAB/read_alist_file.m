function H = read_alist_file(filename)
% H = read_alist_file(filename)

[fileid,message] = fopen(filename, 'rt');
if (fileid == -1)
    error(message)
end

t = fscanf(fileid, '%i',[1,2]);
N = t(1);
M = t(2);

t = fscanf(fileid,'%i', [1,2]);
max_sum_dv = t(1);
max_sum_dc = t(2);

sum_dv = fscanf(fileid, '%i\n', N);
sum_dc = fscanf(fileid, '%i\n', M);

row_idx = reshape(fscanf(fileid, '%i\n', [max_sum_dv,N]).', 1, N*max_sum_dv);
col_idx = reshape(repmat([1:N]',1,max_sum_dv), 1, N*max_sum_dv);

fclose(fileid);

idx = find(row_idx > 0);

H = sparse(row_idx(idx), col_idx(idx), ones(1,numel(idx)), M, N);