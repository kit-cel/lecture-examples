function L = sugiyama(a,fieldSize)
%SUGIYAMA Given a sequence [a_0,a_1,...,a_{2t}] find the smallest LSFR filter
%[L_0,L_1,..,L_nu] such that a_k = sum_i=1^nu L(i)a(k-i) for k > nu
%   Detailed explanation goes here
arguments (Input)
    a = [11 3 8 0]
    fieldSize = 16
end

arguments (Output)
    L
end
    n = length(a); % n = 2t
    x = zeros(1, n + 1); % Initialize output array
    x(1) = 1; % Representing x^{2t}
    [~,~,L] = euclidean_poly(x,a,fieldSize); % Run Extended Euclidean Algorithm
    L = L/L(end); % Normalize the polynomial to have L_0 = 1

    % Remove leading zeros and reverse the order to be compatible with the
    % rest of the code
    idx = find(L~=0,1,"first");
    L = L(idx:end);
    L = flip(L);
end