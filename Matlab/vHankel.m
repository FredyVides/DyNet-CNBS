% Code by: Fredy Vides
% For paper, Sparse autoregressive network model identification 
% for inflation dynamics forecasting
function H=vHankel(y,L)
    N = size(y,1);
    H = [];
    for k = 1:(N-L+1)
        y0 = y(k:(k+L-1),:).';
        H = [H y0(:)];
    end
end