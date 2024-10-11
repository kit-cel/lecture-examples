% Demo Fehlerwahrscheinlichkeit QPSK
% Vorlesung Nachrichtentechnik II
% Institut fuer Nachrichtentechnik (CEL)

clc
clear;
close all;

%Signalenergie
Es=1;

%Signal-zu-Rauschverhaeltnis (SNR) in Dezibel (dB)
SNR=5;

%Rauschenergie
N0=Es/10^(SNR./10);

A=sqrt(Es);

mu=[A 0; 0 A; -A 0; 0 -A];

K=1/pi/N0;

pts = -4:.1:4;
N = length(pts);
X = reshape(repmat(pts,1,N),N,N);
Y = reshape(repmat(pts,N,1),N,N);
Z = K*exp((-(X-mu(1,1)).^2-(Y-mu(1,2)).^2)./N0)+K*exp((-(X-mu(2,1)).^2-(Y-mu(2,2)).^2)./N0) +K* exp((-(X-mu(3,1)).^2-(Y-mu(3,2)).^2)./N0)+K*exp((-(X-mu(4,1)).^2-(Y-mu(4,2)).^2)./N0);
figure;surf(X,Y,Z)

% Symbol 1
Z1 = K*exp((-(X-mu(1,1)).^2-(Y-mu(1,2)).^2)./N0);
figure;surf(X,Y,Z1)
B = zeros(N,N);
for i=1+(N-1)/2:N
    for k=1+(N-1)/2:N-1
        if (i==k)
            B(i,k)=max(Z1(:));
            B(i,41-k+1+40)=max(Z1(:));
        end
    end
end

hold on
surf(X,Y,B)
colormap('jet')