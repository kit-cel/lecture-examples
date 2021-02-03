% construct an irregular matrix
H = Construct_Irregular_LDPC_Code([2 6; 3 12; 6 12], [8 15]);
spy(H);

% generate random information word
u = randi(2,1,15)-1;

% uses a slightly different representation using an upper triangular matrix
% T than in the lecture. Conceptually identical though
x = Encode_LDPC(H, u);

% verify that x is a codeword
mod(H*x,2)