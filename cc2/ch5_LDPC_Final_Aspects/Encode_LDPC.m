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