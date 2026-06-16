function c_hat = RSdec(r,q,t)
%RSDEC Decode a received word r = c + e where c is RS encoded over GF(q)
%   We implicitly assume that the RS code is primitive and narrow sense,
%   i.e. the generator polynomial is g(Z) = (Z-a)(Z-a^2)...(Z-a^{2t}), and
%   that a = gf(2,log2(q))
    arguments (Input)
        r = gf([6 15 10 3 8 6 2 2 5 2 6 8 13 4 5],4);
        q = 2^4
        t = 3
    end
    
    arguments (Output)
        % Decoded codeword
        c_hat
    end

    m = log2(q);

    % Compute Syndromes
    a = gf(2,m);
    S = gf(zeros(1,2*t),m);
    for i = 1:2*t
        S(i) = polyval(r,a^i);
    end
    

    % Error Locator Polynomial
    % Choose your favourite algorithm
    %L = bm_LFSR(S,m);
    L = sugiyama(flip(S),2^m);

    % Find Zeros and Error Positions
    E_L = chienSearch(L,m);
    E_L = 1./E_L; % invert
    E_pos = log(E_L)+1; % get exponent
    Xs = E_L;
    
    % Find Error Values (Forneys Algorithm)
    E_vals = forney(flip(L),flip(S),Xs,t);

    % Reconstruct the error vector
    e = 0*r; % shorthand for e = gf(zeros(1,length(r)),m)
    e(E_pos) = E_vals;

    % Correct the error
    c_hat = r+flip(e);
end