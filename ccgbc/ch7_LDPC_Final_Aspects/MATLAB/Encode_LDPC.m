function x = Encode_LDPC(H, u)
% function x = encode_LDPC(H, u)
%
% H is sparse parity check matrix
% u is input vector
% x is output vector

    persistent Hs;
    
    persistent T Ti A B E C D colintl phi_inv;

    if ~isequal(Hs,H)
        disp('Precomputing step ... will be avoided next time you encode using the same matrix');
        Hs = H;
        [T,A,B,E,C,D,colintl] = greedy_upper_triangulate_forenc(H);
        Ti = invert_upper_diagonal(T);
        temp = E*Ti;
       
        
        while rank(gf(full(mod(C + temp*A,2)),1)) < size(C,1)
            rc1 = randi(size(C,2));
            rc2 = randi(size(D,2));     
            colintl(size(E,2)+[rc1 size(A,2)+rc2]) = colintl(size(E,2)+[size(A,2)+rc2 rc1]);
            tt = A(:,rc1); A(:,rc1) = B(:,rc2); B(:,rc2) = tt;
            tt = C(:,rc1); C(:,rc1) = D(:,rc2); D(:,rc2) = tt;
        end                
        
        Hp = mod([eye(size(temp,2)), zeros(size(temp,2), size(temp,1)); temp, eye(size(temp,1))] * [T,A,B;E,C,D],2);
        
        phi = mod(C + temp*A,2);
        phi_inv_gf = inv(gf(full(phi)));
        phi_inv = double(phi_inv_gf.x);
    end
    
    Heq = [T, A, B; E C D];
        
    p2 = mod(phi_inv*mod(mod(D + E*Ti*B,2)*u(:),2),2);
    p1 = mod(Ti*(A*p2+B*u(:)),2);
    x = [p1;p2;u(:)];
    x(colintl) = x; 
end

function Ti = invert_upper_diagonal(T)
% function Ti = invert_upper_diagonal(T)
% T should be a sparse binary matrix (not gf!)
if size(T,1) < 5 
    tTi = inv(gf(full(T)));
    [i,j] = find(tTi == gf(1));
    % back to sparse
    Ti = sparse(i,j, ones(length(i),1), size(T,1), size(T,2));
else
    sht1 = ceil(size(T,1)/2);
    sht2 = size(T,1)-sht1;
    T11 = T(1:sht1,1:sht1);
    T22 = T((sht1+1):end,(sht1+1):end);
    T12 = T(1:sht1, (sht1+1):end);
    
    T11i = invert_upper_diagonal(T11);
    T22i = invert_upper_diagonal(T22);
    
    Tz = sparse(sht2,sht1);
    Tur = mod(T11i*T12*T22i,2);
    Ti = [T11i, Tur; Tz, T22i];
end
end

function [T,A,B,E,C,D,colintl] = greedy_upper_triangulate_forenc(H)
% carries out greedy upper triangulation as described in Appendix A.2 of Richardson, Urbanke, "Moder Coding Theory"   
resdegrees = full(sum(H,1));
resdegrees(resdegrees==0) = 9999;
n = size(H,2);
k = size(H,2) - size(H,1);

Hg = H;
t=0;
g=0;

colintl = 1:size(H,2);
while t ~= n-k-g
    % EXTEND
    idx = find(resdegrees==min(resdegrees));
    rancol = idx(randi(length(idx)));
    rs = find(Hg(:,rancol))';
    rs = intersect(rs, (t+1):(n-k-g));
    Hg(:,[t+1,rancol]) = Hg(:,[rancol,t+1]);    
    Hg([t+1;rs(1)],:) = Hg([rs(1);t+1],:);
    
    colintl([t+1,rancol]) = colintl([rancol,t+1]);
    
    t = t+1;
    if min(resdegrees) ~= 1
        % CHOOSE
        % increase gap
        Hg([n-k-g-[0:(length(rs)-2)], rs(2:end)],:) = Hg([rs(2:end), n-k-g-[0:(length(rs)-2)]],:);        
        g = g + length(rs) - 1;        
    end
    % update residual degreess
    resdegrees = full(sum(Hg((t+1):(n-k-g),:),1));
    resdegrees(1:(t+1)) = 9999;
    resdegrees(resdegrees == 0) = 9999;
end


Ts = size(H,1)-g;
T = Hg(1:Ts,1:Ts);
A = Hg(1:Ts,(Ts+1):(Ts+g));
B = Hg(1:Ts, (Ts+g+1):end);
E = Hg((Ts+1):end,1:Ts);
C = Hg((Ts+1):end, (Ts+1):(Ts+g));
D = Hg((Ts+1):end, (Ts+g+1):end);

end
