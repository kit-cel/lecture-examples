% dewnsity evolution of a simple regular LDPC code with the particle method
% (population dynamics)
function PDF_Evolution_MonteCarlo_BiAWGN
dv = 3;
dc = 6;
EbN0_dB = 2.1;

% design rate
rd = 1 - dv/dc;

EsN0_dB = EbN0_dB + 10*log10(rd);

% noise standard deviation
sigma_n = sqrt(0.5 * 10.^(-EsN0_dB/10));

% number of decoding iterations
iterations = 11;

% slows down quite a lot! Use with care and only for debugging
display_histogram = true;

% population size
N = 1000000;


% all-zero codeword (+1) with errors
z = 1 + sigma_n*randn(N,1);
 

LLRs = z * 2 / (sigma_n^2);



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
        [H,X] = hist(VN_messages,-10:0.5:30.5);
        H = H./sum(H);
        VNh{i} = H;
        VXh{i} = X;        
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
        [H,X] = hist(CN_messages,-10:0.1:30.5);
        H = H./sum(H);
        CNh{i} = H;
        CXh{i} = X;        
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
         plot(VXh{k},VNh{k});         
         hold all;
     end
     hold off;
     
     subplot(1,2,2);
     for k=1:i
        plot(CXh{k},CNh{k});
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
