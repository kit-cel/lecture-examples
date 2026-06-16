function [q,r] = longDiv(a,b, fieldSize)
%LONGDIV Long division of a(Z) by b(Z) over GF(2^fieldSize)
%   Detailed explanation goes here
arguments (Input)
    a = gf([1 0 0 0 0],4)
    b = gf([11 3 8 0],4)
    fieldSize = 16
end

arguments (Output)
    % a(Z) = q(Z) b(Z) + r(Z) with deg(r)<deg(b)
    q
    r
end

% Make sure b has the same length as a
m = log2(fieldSize);
b = [gf(zeros(1,length(a)-length(b)),m),b];


% Euclidean Algorithm
r = a; % Initialize
q = gf(zeros(size(a)),log2(fieldSize));

while degree(r)>= degree(b)
    % Divide highest order terms
    idxr = find(r~=0,1,'first');
    idxb = find(b~=0,1,'first');
    tmp = r(idxr)/b(idxb);

    % Subtract
    k = degree(r)-degree(b);
    q(end-k) = tmp;
    r = r - tmp*lshift(b,k);
end

end

function out = lshift(b,k)
    idxb = find(b~=0,1,'first');
    
    % If shifted out of range, enlager b
    if idxb - k <= 0
        b = [b(1)*zeros(1,k-idxb+1),b];
        idxb = find(b~=0,1,'first');
    end
    
    b(idxb-k:end-k) = b(idxb:end);
    b(end-k+1:end) = 0;
    out = b;
end

function out = degree(a)
    out = find(a~=0,1,'first');
    out = length(a) - out;
end