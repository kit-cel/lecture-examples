function c = RSenc(message,g,fieldSize)
%RSENC Encode a sequence of q-ary symbols with the RS code generator
%polynomial g
%   message: vector of symbols from GF(2^fieldSize)
%   g: generator polynomial, also a vector over GF(2^fieldSize)
%   representing g_rZ^r + g_{g-1}Z^{r-1} + ... + g_1Z + g_0

arguments (Input)
    % R.h.s. sets default values so the function can be tested in isolation
    message     = gf([5 2 1 6 8 3 10 15 4],4)
    g           = gf([12 10 12 3 9 7 1],4)
    fieldSize   = 2^4
end

arguments (Output)
    % Encoded vector
    c
end

r = length(g)-1; % Degree of g

% Systematic Encoding
message_shifted = lshift(flip(message),r);              % Z^r m(Z)
[~,parity]      = longDiv(message_shifted,g,fieldSize); % Z^r m(Z) mod g(Z)
c               = message_shifted - parity;             % c(Z) = Z^r m(Z) - (Z^r m(Z))
end

% Compute the leftshift of a vector padding with zeros, e.g. lshift([0 0 1 2 3],1) = [0 1 2 3 0] 
% If a non-zero entry is shifted out of range, the vector is expanded, e.g.
% lshift([0 1 2 3],2) = [1 2 3 0 0]
function out = lshift(b,k)
    idxb = find(b~=0,1,'first');
    
    % If shifted out of range, enlager b
    if idxb - k <= 0
        b = [b(1)*zeros(1,k-idxb+1),b];
        idxb = find(b~=0,1,'first');
    end
    
    % Pad with zeros
    b(idxb-k:end-k) = b(idxb:end);
    b(end-k+1:end) = 0;
    out = b;
end