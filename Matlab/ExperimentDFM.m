lt = 1;

n = 5;

if lt == 0 
    N  = 160;
    M = N + 20;
else
    N = 140;
    M = N + 40;
end

th = 0.4;



FechasDFM
GenDatosDFM

time = DatesDFM.Variables;

EDTable = DatosDFM.Variables;
EDTable = EDTable.double;

SEDTable = EDTable(:,[1 [3:2:17 20:2:32 33 34 36 76 78:80]+2 end-2 end]);
TI_data = csvread('TIData.csv');
SEDTable = [SEDTable TI_data];
f = find(isnan(SEDTable));
SEDTable(f) = 0;
if lt ==1
    f0 = find(sum(SEDTable>0)==241);
    SEDTable(:,f0) = log(SEDTable(:,f0));
end
SEDTable0 = SEDTable;

ts = size(SEDTable,2);

c = ones(M,2);
c(:,1) = 0:(M-1);
t = 0:(length(SEDTable)-1);

for k = 1:ts
    d = c\SEDTable(1:M,k);
    mt(k) = max(abs(d(1)*t+d(2)));
end

SEDTable = SEDTable*diag(1./mt);

SEDTableSample = SEDTable(1:N,:);

cv = cov(SEDTableSample);

[u,sc,v] = svd(cv,'econ');
sc = diag(sc);

disp('Porcentaje de variación explicada por el modelo DFM con respecto a los datos de entrenamiento:');
disp(100*sum(sc(1:n))/sum(sc));

[u,~,v] = svd(SEDTableSample,'econ');
p = v(:,1:n);

fact = SEDTable(1:M,:)*p;
w0p = cov(fact);
w0p = (w0p+w0p.')/2;
[w0p,~] = eig(w0p);
fact = fact*w0p;

t0 = smoothdata(fact(:,1),"sgolay",24,"Degree",3);
for k = 2:n
    t0 = [t0 smoothdata(fact(:,k),"sgolay",24,"Degree",3)];
end

W = SEDTable(1:M,:).'*t0;
[uw,sw,vw]= svd(W,'econ');
W = uw*vw.';

fact_w = SEDTable*W;

t1w = smoothdata(fact_w(:,1),"sgolay",24,"Degree",3);
for k = 2:n
    t1w = [t1w smoothdata(fact_w(:,k),"sgolay",24,"Degree",3)];
end

w0w = cov(t1w);
w0w = (w0w+w0w.')/2;
[w0w,~] = eig(w0w);
W = W*w0w;
fact_w = fact_w*w0w;
t1w = t1w*w0w;

dtt1w = [];
for k = 1:n
    d = c\t0(:,k);
    dtt1w = [dtt1w t1w(:,k)-d(1)*t.'-d(2)];
end

fact_p = SEDTable*p;

t1p = smoothdata(fact_p(:,1),"sgolay",24,"Degree",3);
for k = 2:n
    t1p = [t1p smoothdata(fact_p(:,k),"sgolay",24,"Degree",3)];
end

w0p = cov(t1p);
w0p = (w0p+w0p.')/2;
[w0p,~] = eig(w0p);
p = p*w0p;
fact_p = fact_p*w0p;
t1p = t1p*w0p;

dtt1p = [];
for k = 1:n
    d = c\t0(:,k);
    dtt1p = [dtt1p t1p(:,k)-d(1)*t.'-d(2)];
end

dmt = diag(mt);
RSEDTable_w = t1w*W.'*dmt;
RSEDTable_p = t1p*p.'*dmt;
RSEDTable = SEDTable*dmt;

if lt == 0
    Titles = ["L_T(IMAE)" "L_T(IPM(Exp))" "L_T(PIB(Real))" "L_T(IREM)" "L_T(TCP(IPC))" "L_T(TIAR)"];
else
    Titles = ["N_T(IMAE)" "N_T(IPM(Exp))" "N_T(PIB(Real))" "N_T(IREM)" "N_T(TCP(IPC))" "N_T(TIAR)"];
end

k = [1 27 24 15 25 20];

for j=1:6
    figure(1),subplot(2,3,j),plot(time,RSEDTable(:,k(j)),time,RSEDTable_w(:,k(j)),'r.-',time,RSEDTable_p(:,k(j)),'k.-'),
    axis tight,grid on,title(Titles(j)),recessionplot
end

figure(2),
for k = 1:n,subplot(n,1,k),plot(time,fact_w(:,k),time,t1w(:,k)),end

figure(3),plotmatrix(fact_w-t1w)

c0 = corr(SEDTable(1:N,:));
c0 = c0.*(abs(c0)>=th);
c0 = (c0+c0.')/2;
distC0 =  sqrt(2*(1-c0));
distC0 =  (distC0+distC0.')/2;
G = graph(c0);
T = graph(distC0);
T = minspantree(T);
figure(4),subplot(221),plot(G),axis square,title("G_p(t)"),
subplot(222),plot(T),axis square,title("T_p(t)")

c0 = corr(SEDTable(1:M,:));
c0 = c0.*(abs(c0)>=th);
c0 = (c0+c0.')/2;
distC0 =  sqrt(2*(1-c0));
distC0 =  (distC0+distC0.')/2;
G = graph(c0);
T = graph(distC0);
T = minspantree(T);
figure(4),subplot(223),plot(G),axis square,title("G_p(t)"),
subplot(224),plot(T),axis square,title("T_p(t)")

Z = zeros(size(SEDTable,2));
Z(:,20) = c0(:,20);
Z(20,:) = c0(20,:);
G0 = graph(Z);
figure(5),plot(G0)