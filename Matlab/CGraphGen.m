% Code by: Fredy Vides
% For paper, Sparse autoregressive network model identification 
% for inflation dynamics forecasting
function G = CGraphGen(Outputs,Inputs,th)
Cm = corr([Outputs Inputs]);
Cm = abs(Cm)>th;
sCm = sum(Cm);
s=[];
t=[];
for k = 1:length(Cm) 
    s = [s k*ones(1,sCm(k))];
    t=[t find(Cm(k,:))];
end
G = digraph(s,t);
plot(G)
end