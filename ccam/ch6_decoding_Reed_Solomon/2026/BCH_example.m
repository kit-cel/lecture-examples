function BCH_example()
%BCH_example Example run. Create a BCH codes, encode a message, introduce
%errors, decode
    
    % Example over GF(16)
    fieldSize = 8;

    % Primitive Element (2 is always primitive in matlabs GF implementation)
    a = gf(2,fieldSize);
    t = 3;

    % Find conjugacy classes
    classes = {};
    used = [];
    for i = 1:2*t
        beta = a^i;
        if ~isempty(find(used==beta,1))
            continue
        end
        class = [];
        while isempty(find(class==beta,1))
            class   = [class,beta];
            used = [used,beta];
            beta    = beta^2;
        end
        classes{i} = class;
    end
    
    % Build generator polynomial
    M = 1;
    for i = 1:length(classes)
        elements = classes{i};
        for j = 1:length(elements)
            M = conv(M,[-elements(j) 1]);
        end
    end
    g = M;

    % Print g in nice form
    %coeffsToPoly(g.x,'Z')
    
    % Create message
    n = 2^fieldSize - 1;
    k = n - length(g)+1; % degree of g is length(g) - 1
    m = randi([0,1],1,k);
    c = BCHEncode(g,m,n);
    

    %Create n_mc random Error Patterns. If n_e <= t all should be corrected
    n_e = 3; % No. of Errors
    n_mc = 50; % No. of patterns
    E_pos = zeros(n_mc,n_e);
    for nn = 1:n_mc
        E_pos(nn,:) = randperm(n,n_e);
    end

    % Simulation loop to count the number of decoding errors
    bl_errs = 0;
    for i = 1:size(E_pos,1)
        e = zeros(1,n);
        e(E_pos(i,:)) = 1; % set error pattern
        r = c + e; % apply error
        c_hat = BCHDec(r,fieldSize,2*t); % decode
        berrs = sum(c_hat.x~=c.x); % count bit errors
        if berrs>0
            bl_errs = bl_errs + 1; % count block errors
        end
    end
    fprintf("A total of %d block errors in %d runs \n",bl_errs,n_mc);
end

% Simple non-systematic encoding by polynomial multiplication
function out = BCHEncode(g,m,n)
    g = [g,zeros(1,n-length(g))];
    c = zeros(1,n);
    for i = 1:length(m)
        c = c + m(i)*circshift(g,i-1);
    end
    out = c;
end

% Function for nice printing of polynomial vectors coeffsToPoly([1 2
% 3],'Z') prints 1 + 2*Z + 3*Z^{2}
function s = coeffsToPoly(c, varname)
% coeffsToPoly  Return polynomial string from coefficient vector
%   s = coeffsToPoly(c)         % default variable 'x'
%   s = coeffsToPoly(c,'X')     % custom variable name
if nargin < 2 || isempty(varname), varname = 'x'; end
c = c(:).';                % ensure row
n = numel(c) - 1;
terms = {};

for k = 0:n
    a = c(k+1);
    if a == 0
        continue
    end
    % Coefficient string
    if k == 0
        coefStr = sprintf('%g', a);
    else
        if a == 1
            coefStr = '';
        elseif a == -1
            coefStr = '-';
        else
            coefStr = sprintf('%g*', a);
        end
    end
    % Power string
    if k == 0
        powStr = '';
    elseif k == 1
        powStr = varname;
    else
        powStr = sprintf('%s^{%d}', varname, k);
    end
    terms{end+1} = [coefStr powStr]; %#ok<AGROW>
end

if isempty(terms)
    s = '0';
    return
end

% Join terms handling signs
s = terms{1};
for i = 2:numel(terms)
    term = terms{i};
    % If term starts with '-', it's negative; otherwise prepend ' + '
    if startsWith(term, '-')
        s = [s ' - ' term(2:end)];
    else
        s = [s ' + ' term];
    end
end
end