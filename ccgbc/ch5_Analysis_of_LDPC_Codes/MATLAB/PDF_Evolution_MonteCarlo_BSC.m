% dewnsity evolution of a simple regular LDPC code with the particle method
% (population dynamics)
function PDF_Evolution_MonteCarlo_BSC
dv = 3;
dc = 6;
delta = 0.05;

% number of decoding iterations
iterations = 10;

% slows down quite a lot! Use with care and only for debugging
display_histogram = true;

% population size
N = 1000000;


% all-zero codeword (+1) with errors
z = 1-2*(rand(N,1) < delta);

Lc = log((1-delta)/delta);

LLRs = z * Lc;



CN_messages = zeros(N,1);
for i=1:iterations
    % compute VN to CN messages
    if i == 1
        VN_messages = LLRs;
    else
        % VN operation
        VN_messages = LLRs(randi(N,N,1)) + sum(CN_messages(randi(N,N,dv-1)),2);
    end
    
    if display_histogram
        % histograms for display/debug purposes
        if i <= 2
            um = unique(VN_messages);
            umc = arrayfun(@(x)sum(VN_messages == x), um, 'UniformOutput', true);
            VNh{i} = umc/sum(umc);
            VXh{i} = um;
        else
            [H,X] = hist(VN_messages,-10:0.5:30.5);
            H = H./sum(H);
            VNh{i} = H;
            VXh{i} = X;        
        end
    end
    
    if ~display_histogram && all(VN_messages > 0)
        % no more errors, abort
        fprintf('No errors after %d iterations', i);
        break;
    end
    
    CN_messages = zeros(N,1);
    % here, implement the sum product rule
    idx = randi(N,N,dc-1);
    CN_messages = phi(sum(phi(abs(VN_messages(idx))),2)) .* prod(sign_corr(VN_messages(idx)),2);
    
    if display_histogram        
        if i <= 2
            um = unique(CN_messages);
            umc = arrayfun(@(x)sum(CN_messages == x), um, 'UniformOutput', true);
            CNh{i} = umc/sum(umc);
            CXh{i} = um;
        else    
            [H,X] = hist(CN_messages,-10:0.5:30.5);
            H = H./sum(H);
            CNh{i} = H;
            CXh{i} = X;        
        end
    end
    
    if all(VN_messages > 0)
        % no more errors, abort
        fprintf('No errors');
        break;
    end
    
end
 if display_histogram
     figure(1);
     subplot(1,2,1);
     for k=1:i
         if k <= 2
             stem(VXh{k}, VNh{k});
         else
             plot(VXh{k},VNh{k});
         end
         hold all;
     end
     hold off;
     
     subplot(1,2,2);
     for k=1:i
        if k <= 2
            stem(CXh{k}, CNh{k});
        else
            plot(CXh{k},CNh{k});
        end
        hold all;
     end
     hold off;
     shg    
 end
end

function s = sign_corr(x)
s = sign(x);
s(x == 0) = (-1).^randi([0 1], size(s(x==0)));
end


%INF_CONDITION = [0, 692];
function y = phi(x)
y = log(coth(x/2));
y(x  <= 1e-300) = 692;
end
%phi = @(x)((x > 1e-300)*log(coth(x/2)) + INF_CONDITION(1+(x <= 1e-300)));
