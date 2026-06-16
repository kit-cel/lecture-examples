function E_vals = forney(L,S,Xs,t)
%FORNEY Find the error values given syndromes S and an error locator
%polynomial L, at error locators Xs_i for an RS code with distance 2t+1
% Using the Forney formula E(X_i) = - Omega(X_i^{-1})/L'(X_i^{-1}) with
% Omega(Z) = L(Z)S(Z) mod Z^{2t}
arguments (Input)
    L = gf([2 4 1],3)
    S = gf([3 6 3 5],3)
    Xs = gf([3 7],3)
    t = 2
end

arguments (Output)
    E_vals
end
    % Compute Omega
    Omega = conv(L,S);
    Omega = Omega(end-2*t+1:end);

    % Compute the coefficients of the formal derivative L'(Z)
    d_L = real_mult(L(1:end-1),flip(1:length(L)-1));

    E_vals = 0*Xs; % Shorthand initialisation

    % Use Forney's formula
    for k = 1:length(Xs)
        E_vals(k) = -polyval(Omega,1/Xs(k))/polyval(d_L,1/Xs(k));
    end
    
end

% Compute entrywise real multiplication of two vectors n and x: 
% n_i*x_i = x_i + x_i + ... + x_i, n_i times. For binary
% extension fields this is x_i for odd n_i and 0 for even n_i
function out = real_mult(x_gf,n)
    out = 0*x_gf;
    for i = 1:length(x_gf)
        for j = 1:n(i)
            out(i) = out(i) + x_gf(i);
        end
    end
end