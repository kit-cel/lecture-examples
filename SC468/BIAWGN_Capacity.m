function BIAWGN_Capacity

esno_range = -5:0.5:10;

C_shannon = 0.5*log2(1 + 2*10.^((esno_range)/10.0));
C_bpsk_soft = channelCapacity_BIAWGN(esno_range);

% compute for W = 2
esno_lin = 10.^(esno_range/10);
% equivalent error rate of BSC
delta = 0.5*erfc(sqrt(esno_lin));

% BSC capacity
hb = @(x)(-x.*log2(x)-(1-x).*log2(1-x));

figure(2);
plot(esno_range, C_bpsk_soft);
hold on;
plot(esno_range, 1-hb(delta));

% get capacity of quantized channel
C4 = quantized_capacity(esno_range, 4);
C8 = quantized_capacity(esno_range, 8);
C16 =  quantized_capacity(esno_range, 16);

plot(esno_range, C4, '--');
plot(esno_range, C8, '--');
plot(esno_range, C16, '--');

plot(esno_range, C_shannon, 'k-');

legend('$W \to \infty$ (soft dec.)', '$W = 2$ (hard dec.)', '$W = 4$ (2 bit)', '$W = 8$ (3 bit)', '$W = 16$ (4 bit)','Location','SouthEast','Interpreter','latex');
grid on;
hold off;
ylim([0 1.2]);
end



% helper function, compute capacity of BIAWGN channel using Gauss-Hermite
% quadrature
function cap = channelCapacity_BIAWGN(esno_range)
% Gauss-Hermite polynomial taken from Python's numpy library np.polynomial.hermite.hermgauss(40)
x_GH = [-8.09876114 -7.41158253 -6.84023731 -6.32825535 -5.85409506 -5.40665425 -4.97926098 -4.56750207 -4.16825707 -3.77920675 -3.39855827 -3.02487988  -2.656996   -2.29391714 -1.93479147 -1.57886989 -1.22548011 -0.87400661  -0.52387471 -0.17453721  0.17453721  0.52387471  0.87400661  1.22548011   1.57886989  1.93479147  2.29391714  2.656996    3.02487988  3.39855827   3.77920675  4.16825707  4.56750207  4.97926098  5.40665425  5.85409506   6.32825535  6.84023731  7.41158253  8.09876114];
w_GH = [2.59104371e-29 8.54405696e-25 2.56759337e-21 1.98918101e-18 6.00835879e-16 8.80570765e-14 7.15652805e-12 3.52562079e-10 1.12123608e-08 2.41114416e-07 3.63157615e-06 3.93693398e-05 3.13853595e-04 1.87149683e-03 8.46088801e-03 2.93125655e-02 7.84746059e-02 1.63378733e-01 2.65728252e-01 3.38643277e-01 3.38643277e-01 2.65728252e-01 1.63378733e-01 7.84746059e-02 2.93125655e-02 8.46088801e-03 1.87149683e-03 3.13853595e-04 3.93693398e-05 3.63157615e-06 2.41114416e-07 1.12123608e-08 3.52562079e-10 7.15652805e-12 8.80570765e-14 6.00835879e-16 1.98918101e-18 2.56759337e-21 8.54405696e-25 2.59104371e-29];

% channel transition probability
f_YgivenX = @(y,x,sigman)(exp(-((y-x).^2)/(2*sigman.^2))./sqrt(2*pi)./sigman);
% probability of channel output assuming equiprobable channel inputs
f_Y = @(y,sigman)(0.5*(f_YgivenX(y,+1,sigman)+f_YgivenX(y,-1,sigman)));


sigman_range = sqrt(0.5*10.^(-esno_range/10));
cap = zeros(size(esno_range));
for si = 1:numel(sigman_range)    
    integral_xplus1 = sum(w_GH .* [log2(f_Y(sqrt(2)*sigman_range(si)*x_GH + 1, sigman_range(si)))]);    
    integral_xminus1 = sum(w_GH .* [log2(f_Y(sqrt(2)*sigman_range(si)*x_GH - 1, sigman_range(si)))]);    

    integral = (integral_xplus1 + integral_xminus1)/2/sqrt(pi);
    cap(si) =  -integral - 0.5*log2(2*pi*exp(1)*sigman_range(si)^2);
end
end

% compure capacity of quantized channel. Assuming that the channel is
% (weakly) symmetric, we can use a simplified version
function C = quantized_capacity(esno_range, W)
C = zeros(size(esno_range));

if mod(W,2) ~= 0 || W < 4
    error('only works for even W > 2');
end

for k=1:numel(esno_range)
    esnolin = 10.^(esno_range(k)/10);
    sigman =sqrt(1/2/esnolin);
    
    % test a bunch of quantization intervals. Quantize channel output uniformly between
    % -max_x and +max_x for 1000 different values and select the one
    % maximizing the capacity
    max_x = linspace(0.5,3,1000);
    
    Ctemp = zeros(size(max_x));
    for mi = 1:numel(max_x)
        delta = max_x(mi)/((W-2)/2);
        P = zeros(W,2);
        start = -max_x(mi);
        for j=1:W
            if j == 1
                P(j,1) = qfunc((-10-(+1))/sigman) - qfunc((start-(+1))/sigman);
                P(j,2) = qfunc((-10-(-1))/sigman) - qfunc((start-(-1))/sigman);
            elseif j == W
                P(j,1) = qfunc((start-(+1))/sigman) - qfunc((start+10-(+1))/sigman);
                P(j,2) = qfunc((start-(-1))/sigman) - qfunc((start+10-(-1))/sigman);
            else
                P(j,1) = qfunc((start-(+1))/sigman) - qfunc((start+delta-(+1))/sigman);
                P(j,2) = qfunc((start-(-1))/sigman) - qfunc((start+delta-(-1))/sigman);
                start = start+delta;
            end
        end
        Ctemp(mi) = symmetric_channel_capacity(P);
    end
    C(k) = max(Ctemp);
end

end

% helper function that computes the capacity of a symmetric channel
function C = symmetric_channel_capacity(P)
n = size(P,2); % input 
m = size(P,1); % output

P = P./repmat(sum(P),m,1); 
c = P.*log2(P);
c(isnan(c)) = 0;
c = sum(c)';

x = ones(n,1)/n;

C = c'*x + sum(-(P*x).*log2(P*x));
end
