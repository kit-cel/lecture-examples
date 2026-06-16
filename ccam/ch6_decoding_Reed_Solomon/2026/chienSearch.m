function pol_zeros = chienSearch(poly,fieldSize)
%CHIENSEARCH Find the zeros of poly(Z) by evaluating all field elements
%of GF(2^fieldSize) using chien search
%   Detailed explanation goes here
arguments (Input)
    poly = gf([1,1,0,6],4)
    fieldSize = 4
end

arguments (Output)
    % Vector Z_i such that poly(Z_i) = 0
    pol_zeros
end
    % Initialize registers
    pol_zeros = [];
    a = gf(2,fieldSize);
    A = a.^(0:length(poly)-1);
    v_ones = gf(ones(size(poly)),fieldSize);
    base = poly;

    % Test Z = 1
    if v_ones*base'==0
        pol_zeros = [pol_zeros,0];
    end

    % Chien Search
    for i = 1:2^fieldSize-1
        base = base.*A;
        val = v_ones*base'; % sum of base coefficients
        if val==0
            pol_zeros = [pol_zeros,i];
        end
    end
    pol_zeros = a.^pol_zeros;

    % Brute-force evaluation would look like this:
    %vals = polyval(flip(poly),0:2^fieldSize-1);
    %pol_zeros = gf(find(vals == 0)-1,fieldSize);
end