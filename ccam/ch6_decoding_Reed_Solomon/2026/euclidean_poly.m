function [g,s,t] = euclidean_poly(a,b, fieldSize)
%ECULIDEAN_INT Compute g,s, and t such that a(x)s(x) + b(x)t(x) = g(x)
%   Detailed explanation goes here
arguments (Input)
    % Assume length(a)>=length(b)
    a = [1 0 0 0 0 0 0]
    b = [3 2 2 4 3 2]
    fieldSize = 5
end

arguments (Output)
    g % gcd such that a(x)s(x) + b(x)t(x) = g(x)
    s  
    t
end
    m = log2(fieldSize);
    a = gf(a,m);
    b = gf(b,m);

    % Equalize lengths
    n = length(a);
    bz = gf(zeros(size(a)),m);
    bz(n-length(b)+1:n) = b;
    b = bz;
    s_old = gf(zeros(size(a)),m);
    t_old = gf(zeros(size(a)),m);

    % Initialize
    s = s_old;
    t = t_old;
    s_old(end) = 1;
    t_old(end) = 0;
    s(end) = 0;
    t(end) = 1;
    
    r_old = a;
    r = b;
    
    while ~isempty(degree(r))
        
        
        r_temp = r;
        s_temp = s;
        t_temp = t;
        
        [q,r] = longDiv(r_old,r,fieldSize);
        if isprime(fieldSize)
            s = mod(s_old - longMultGF(q,s,fieldSize,n),fieldSize);
            t = mod(t_old - longMultGF(q,t,fieldSize,n),fieldSize);
        else
            s = s_old - longMultGF(q,s,fieldSize,n);
            t = t_old - longMultGF(q,t,fieldSize,n);
        end
        
        r_old = r_temp;
        t_old = t_temp;
        s_old = s_temp;    

        if degree(r)<(n-1)/2 % n = 2t + 1
            break
        end
    end
    % Normalise gcd
    idx = find(r~=0,1,'first');
    leading = r(idx);


    g = r/leading;
    s = s/leading;
    t = t/leading;


    % See if indeed a(Z)s(Z) + b(Z)t(Z) = g(Z)
    %longMultGF(a,s,fieldSize,n) + longMultGF(b,t,fieldSize,n)
end

% For q prime, for q=p^m replace by GF operations
function out = longMultGF(a,b,q,n)
    out = conv(gf(a,log2(q)),gf(b,log2(q)));

    % Trim/expand zeros to have length n
    z = gf(zeros(1,n),log2(q));
    if length(out) > n
        z = out(end-n+1:end);
    else
        z(n-length(out)+1:end) = out;
    end

    out = z;
end


function out = degree(a)
    out = find(a~=0,1,'first');
    out = length(a) - out;
end