function [sys1,sys2] = InflationPredictor(S,N,T,Tv,n1,n2,alpha,type)
% Code by: Fredy Vides
% For paper, Sparse autoregressive network model identification 
% for inflation dynamics forecasting
% Examples:  
% [sys1,sys2] = InflationPredictor(49,28,10,2,12,[18 2]);
% [sys1,sys2] = InflationPredictor(49,28,10,2,12,[18 2:4]);
% [sys1,sys2] = InflationPredictor(49,39,10,2,16,[22 1]);
% [sys1,sys2] = InflationPredictor(49,39,10,2,16,[22 1:3]);

if nargin <=6
    alpha = abs(randn);
    type = 0;
elseif nargin < 8
    type = 0;
end

format long g;
%Outputs = csvread("InflationOutputs.csv");
Outputs = csvread("HNInflationData.csv");
% Outputs := hn
%Inputs = csvread("InflationInputs.csv");
Inputs = csvread("RCInflationData.csv");
% Inputs := [us gt sv ca cr]

scale = norminv(1-.01/2);

L = length(Outputs);

S = min(L,S);

L = S;

Ntrain = min(N,L-1);

T = min(S-Ntrain,T);

t = (Ntrain:(Ntrain+T))+1979;

data_id = iddata(Outputs(1:Ntrain,1),[],1,'TimeUnit', 'years');

sys1 = arx(data_id,n1);

H = vHankel(Outputs(1:Ntrain,1),n1);

H0 = H(:,1:(end-1));

H1 = H(end,2:end);

E00 = eye(size(H0,1));

if type == 0
    Wp = SpSolver(H0*H0.'+alpha*E00,H0*H1.',20,1e-6,1,1e-3).';
    sys1.A(2:end) = -fliplr(Wp);
end

[yF, ~, ~, yFSD]  = forecast(sys1, data_id, T);

figure(1),subplot(211),plot((1:(Ntrain+T))+1979,Outputs(1:(Ntrain+T)),'b',t, [Outputs(Ntrain);yF.OutputData], 'r.-', t, [Outputs(Ntrain);yF.OutputData]+[0;scale*yFSD], 'r--', t, [Outputs(Ntrain);yF.OutputData]-[0;scale*yFSD], 'r--')
axis tight
grid on,xlabel('t (years)'),ylabel('Inflation rate (%)'),legend('Reference data','Predictions')
mv = min([Outputs(Ntrain);yF.OutputData]-[0;scale*yFSD]);
Mv = max(Outputs);


Y = Outputs;

ln2 = length(n2)-1;

for k = 1:ln2
    Y = [Y Inputs(:,n2(k+1))];
end

Yt = Y(1:Ntrain,:);

sys2 = varm(ln2+1,n2(1));
H = vHankel(Yt,n2(1));

H0 = H(:,1:(end-1));

H1 = H((end-ln2):end,2:end);

E0 = eye((ln2+1)*n2(1));

if type == 0
    W = SpSolver(H0*H0.'+alpha*E0,H0*H1.',prod(size(H1)),1e-3,1,5e-3).';
else
    W = H1*H0.'/(H0*H0.'+alpha*E0);
end

Cw = sparse((ln2+1)*n2(1),(ln2+1)*n2(1));

Cw(1:(end-ln2-1),(ln2+2):end) = speye((ln2+1)*n2(1)-ln2-1);

Cw((end-ln2):end,:) = W;

Cs = H0(:,1);

M = min(ceil(Ntrain/2),size(H0,2)-1);

for k = 1:M
    Cs = [Cs Cw*Cs(:,k)];
end

Cs = Cs((end-ln2):end,:);
Cs0 = H0((end-ln2):end,1:(M+1))-Cs((end-ln2):end,:);

sys2.Covariance = cov(Cs0.');

sys2.Constant = zeros(ln2+1,1);

for k = 1:n2(1), sys2.AR{k} = W(:,k:(k+ln2));end

[YF,YFSD] = forecast(sys2,T,Yt);

Tv = min(Tv,T-1);

q0 = [Cs';YF(1:Tv,:)];
nq = size(q0,2);
q1 = [H0((end-ln2):end,1:(M+1)).';Y((Ntrain+1):(Ntrain+Tv),:)];
options.tol = 1e-6;
q = lsqnonneg([(q0.'*q0+alpha*eye(nq));ones(1,nq)],[(q0'*q1(:,1));1],options)

sum(q)

extractMSE = @(x)diag(x)';
MSE = cellfun(extractMSE,YFSD,UniformOutput=false);
SE = sqrt(cell2mat(MSE));
SE = SE(:,1);

YFp = YF*q;

figure(1),subplot(212),plot((1:(Ntrain+T))+1979,Outputs(1:(Ntrain+T)),'b',t, [Outputs(Ntrain);YFp], 'r.-', t, [Outputs(Ntrain);YFp]+[0;scale*SE], 'r--', t, [Outputs(Ntrain);YFp]-[0;scale*SE], 'r--')
xlim([1980 1980+length(Outputs)-1]),ylim([mv Mv]), grid on,xlabel('t (years)'),ylabel('Inflation rate (%)'),legend('Reference data','Predictions')

figure(2),subplot(121),G = CGraphGen(Outputs(1:Ntrain,:),Inputs(1:Ntrain,:),0.2);axis square,axis tight

figure(2),subplot(122),G = CGraphGen(Yt(:,1),Yt(:,2:end),0.22);axis square,axis tight

figure(3),subplot(121),spy(Cw),axis tight,grid on

figure(3),subplot(122),plot(exp(2*pi*i*(0:1/100:1))),hold on,plot(eig(full(Cw)),'r.'),hold off,axis equal,axis tight,grid on