function [A, jb] = gfrref(A)
%GFRREF   Reduced row echelon form for Galois field matrices
%   R = GFRREF(A) produces the reduced row echelon form of A.
%
%   [R,jb] = GFRREF(A) also returns a vector, jb, so that:
%       r = length(jb) is this algorithm's idea of the rank of A,
%       x(jb) are the bound variables in a linear system, Ax = b,
%       A(:,jb) is a basis for the range of A,
%       R(1:r,jb) is the r-by-r identity matrix.
%
%   Based on rref.m, copyright 1984-2017 The MathWorks, Inc. 
%   modified to work with matrices over GF(2^m) and removed many of the
%   checks

[m,n] = size(A);

% do not use any tests as we know what we are inputting

% Loop over the entire matrix.
i = 1;
j = 1;


jb = zeros(1,0);


while i <= m && j <= n
   % Find value and index of first non-zero element in the remainder of column j.
   temp = A(i:m, j);
   k = find(temp.x ~= 0, 1);
   if isempty(k)
       % no element found, move on
       j = j + 1;
   else
      k = k+i-1;
      % Remember column index
      jb(end+1) = j; 
      % Swap i-th and k-th rows.
      A([i k],j:n) = A([k i],j:n);
      % Divide the pivot row by the pivot element.
      A(i,j:n) = A(i,j:n)./A(i,j);
      % Subtract multiples of the pivot row from all the other rows.
      for k = [1:i-1 i+1:m]
         A(k,j:n) = A(k,j:n) - A(k,j).*A(i,j:n);
      end
      i = i + 1;
      j = j + 1;
   end
end


