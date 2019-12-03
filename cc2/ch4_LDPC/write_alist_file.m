function write_alist_file(H, filename)
% write_alist_file(H, filename)

if((exist(filename)==2))
    over=input(['Warning file: ' filename ' exists!!! Overwrite y/n [n]?'],'s');
    if(isempty(over)|(over=='n')|(over=='N'))
        return;
    end
end

[fileid,message]=fopen(filename,'wt');
if (fileid==-1)
    error(message)
end

[M N]=size(H);

% write dimensions of H
fprintf(fileid,'%i %i\n',N,M);

sum_dv = full(sum(H,1));
sum_dc = full(sum(H,2));

% write max weights
%fprintf(fileid,'%i %i\n',max(sum2),max(sum1));
fprintf(fileid,'%i %i\n',max(sum_dv),max(sum_dc));

fprintf(fileid, '%i ', sum_dv);
fprintf(fileid, '\n');
fprintf(fileid, '%i ', sum_dc);
fprintf(fileid, '\n');
fclose(fileid);

[ind, col] = find(H);
col = (col-1) * max(sum_dv) + 1;
for k=1:(max(sum_dv)-1)
    t = [1;    col(2:end)-col(1:end-1)];
    col(t==0) = col(t==0) + 1;
end
dvidx = zeros(max(sum_dv), N);
dvidx(col) = ind;   

dlmwrite(filename, dvidx',  '-append', 'delimiter', ' ', 'precision','%d');




[ind, row] = find(H');
row = (row-1) * max(sum_dc) + 1;
for k=1:(max(sum_dc)-1)
    t = [1;    row(2:end)-row(1:end-1)];
    row(t==0) = row(t==0) + 1;
end
dcidx = zeros(max(sum_dc), M);
dcidx(row) = ind;   

dlmwrite(filename, dcidx',  '-append', 'delimiter', ' ', 'precision','%d');


