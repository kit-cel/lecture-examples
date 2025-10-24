% This code is provided as supplementary material of the lecture Channel Coding: Graph-based Codes.
%
% This code illustrates
%
%     Calculating the channel capacity using convex optimization
%     ATTENTION: in order to run this code, it is advised that you to have CVX installed.
%     You can get it (Linux/Windows/macOS) at http://cvxr.com


% first example, weakly symmetric channel
fprintf('Weakly symmetric channel:\n');
P = [1/3, 1/3;
     1/2, 1/6;
     1/6, 1/2];

[C, px] = channel_capacity(P);
fprintf('Capacity achieving input distribution:\n');
px
fprintf('Capacity:\n');
C

% second example, non-symmetric channel from lecture
fprintf('Arbitrary channel:\n');
P = [1/2, 1/8;
     1/3, 5/8;
     1/6, 1/4];
 
[C, px] = channel_capacity(P);
fprintf('Capacity achieving input distribution:\n');
px
fprintf('Capacity:\n');
C


% Example: Z-channel
q = 0.1;
P = [1, q;
     0, 1-q];
[C, px] = channel_capacity(P);
fprintf('Capacity achieving input distribution:\n');
px
h = @(x)(-x.*log2(x) - (1-x).*log2(1-x));
u = ((1-q)+(1-q)*2^(h(q)/(1-q)))^(-1);
fprintf('From analytical formula: [%1.4f %1.4f]\n', 1-u, u);
fprintf('Capacity:\n');
C
fprintf('From analytical formula: %1.4f\n',h(u*(1-q))-u*h(q));
 

% check if channel is weakly symmetric
function retval = is_weakly_symmetric(P)
    V = size(P,2);
    W = size(P,1);

    retval = true;
    % first check if all columns are permutations of each other
    col1 = sort(P(:,1));
    for k=2:size(P,2)
        if any(col1 ~= sort(P(:,k)))
            retval = false;
        end
    end
    
    % now check if row sums are equal
    row_sums = sum(P,2);
    if ~all(row_sums == row_sums(1))
        retval = false;
    end    
    
end

% calculate the capacity and the capacity-achieving input distribution for
% a DMC specified by a transition matrix P
function [C, px] = channel_capacity(P)    
    V = size(P,2);
    W = size(P,1);
    
    if is_weakly_symmetric(P)
    % use the formula for weakly symmetric channels
        px = ones(1,V)/V;
        col1 = P(:,1);
        Hc = -sum(col1 .* log2(col1));
        C = log2(W) - Hc;
    else
        P_tilde = P .* log2(P);
        P_tilde(isnan(P_tilde)) = 0;  % 0 * log2(0) = 0
                
        
        if exist('cvx_begin') > 0
        
            % the optimization magic starts here
            cvx_expert true
            cvx_begin quiet
                variable px(V)
                py = P*px;    % output probability
                maximize (sum(P_tilde, 1)*px  + sum(entr(P*px))/log(2))
                subject to
                    px >= 0
                    sum(px) == 1
            cvx_end
            C = cvx_optval;        
        else % cvx is not available, use built in fmincon function from MATLAB. The interface is a little bit less beautiful than the cvx interface
            myentr = @(x)(max([-x.*log(x)],0));
            [px,C] = fmincon(@(px)(-sum(P_tilde, 1)*px  - sum(myentr(P*px))/log(2)), ones(V,1)/V, [], [], ones(1,V), 1, zeros(V,1), ones(V,1), [], optimoptions(@fmincon, 'Display', 'off'));
            C = -C;
        end
    end
    
end


