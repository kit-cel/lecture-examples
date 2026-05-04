function simulate

remove_4cycles = true;

dv = 3;
dc = 6;
n = 200;

eps_range = linspace(0.1,0.43, 15);

epochs = linspace(200, 200, numel(eps_range));


frames = 3000;

errors = [];

epsi = 1;
for epsilon = eps_range
    epsilon
    num_err = 0;
    for ei = 1:epochs(epsi)
        H = getRegularH(n, dv, dc);
        if remove_4cycles
            [H,~] = remove_4cycle(H);
        end
        for fi = 1:frames
            erasures = find(rand(1,n) < epsilon);
            E = Peeling_Decoder(H, erasures);
            if E > 0
                num_err = num_err + 1;
            end
         end
    end
    errors(end+1) = num_err / epochs(epsi) / frames;
    epsi = epsi + 1;
end
figure(1);semilogy(eps_range, errors);
format long
[eps_range(:), errors(:)]
end

function E = Peeling_Decoder(Hcl, erasures)
% function E = Peeling_Decoder(Hcl, erasures)
%
% returns reisudla number of erasures after peeling decoding, erasures is a
% binary (logical) vector indicating at which vns erasures are present
% note that the performance of the peeling decoder is equivalent to the one
% of the BEC message passing decoder, as has been proven in M. Stinner, L.
% Barletta, P. Olmos, "Finite-length scaling based on belief propagation
% for spatially coupled LDPC codes", https://arxiv.org/pdf/1604.05111.pdf

H = Hcl(:,erasures);
% find degree-1 checks
dc1checks = sum(H,2) == 1;
[~,assoc_vns] = find(H(dc1checks,:));

assoc_vns = assoc_vns(:).';

E = size(H,2);

while numel(assoc_vns) > 0
    H(:,assoc_vns) = [];
    
    dc1checks = sum(H,2)==1;
    [~,assoc_vns] = find(H(dc1checks,:));
    assoc_vns = assoc_vns(:).';
    E = size(H,2);
end
end




