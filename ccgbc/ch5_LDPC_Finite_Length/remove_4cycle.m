function [Hout,success] = remove_4cycle(H)
% removes cycles of length 4
success = true;
N = size(H,2);
M = size(H,1);

abort = false;
trials = 0;
while ~abort
    T = H'*H - diag(sum(H,1));
    [v1s,v2s] = find(triu(T) >= 2);
    %numel(v1s)
    if isempty(v1s)
        abort = true;
    end
    for k=1:numel(v1s)
        v1 = v1s(k);
        v2 = v2s(k);
        
        % common check node
        cncs = find(H(:,v1) + H(:,v2) == 2);
        if numel(cncs) >= 2
            cn = find(H(:,v1)+H(:,v2)==2,1);
            find_replace = false;
            while ~find_replace
                vni = myrandsample(setdiff(1:N,[v1 v2]),1);
                if all(H(:,v2).*H(:,vni) == 0)
                    cnt = myrandsample(find(H(:,vni)),1);
                    % swap edges
                    temp = H([cnt,cn],v2);
                    H([cnt,cn],v2) = H([cnt,cn],vni);
                    H([cnt,cn],vni) = temp;
                    find_replace = true;
                end
            end
        end
    end
    trials = trials + 1;
    if trials > 65
        success = false;
        abort = true;
    end
end
Hout = H;
end

function u = myrandsample(v,num)
idx = randi(numel(v),1,num);
u = v(idx);
end
