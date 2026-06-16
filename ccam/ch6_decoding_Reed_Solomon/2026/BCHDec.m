function xhat = BCHDec(r,m,t)
%BCHDEC Decode a received word r = c + e where c is BCH encoded with base field GF(2^m)
%   We implicitly assume that the BCH code is primitive and narrow sense,
%   i.e. the generator polynomial is g(Z) = (Z-a)(Z-a^2)...(Z-a^{2t})*'all conjugate zeros', and
%   that a = gf(2,m)
%   This is mostly a copy 
arguments (Input)
    r (1,:) {mustBeA(r,'gf')} = gf([0,1,0,1,0,0,0,0,1,0,0,0,0,0,0],4)
    m = 4
    t = 6
end

arguments (Output)
    xhat
end
    n = 2^m - 1;
    a = gf(2,m);

    % Compute Syndrome (only for a^1,...,a^{2t})
    S = gf(zeros(1,t),m);
    for i = 1:t
        S(i) = polyval(flip(r),a^i);
    end
    
    % Choose you favorite algorithm to compute the error locator polynomial
    L = bm_LFSR(S,m);
    %L = sugiyama(flip(S),2^m);

    % Find the zeros of L(Z) via chien search
    E_L = chienSearch(L,m);
    E_L = 1./E_L;
    E_pos = log(E_L); % Get the exponent
    E = zeros(1,n); 
    E(E_pos+1) = 1; % Reconstruct error
    xhat = r + E; % Correct the error
end