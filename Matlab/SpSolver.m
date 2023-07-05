% Code by: Fredy Vides
% For paper, Sparse Autoregressive Reservoir Computers for 
% Dynamic Financial Processes Identification
%%
%%
function X=SpSolver(A,Y,L,tol,sl,delta)
if nargin<=4 
	sl=1;delta=tol;
elseif nargin<=5 
	delta=tol;
end
N=size(Y,2);
%if length(L)==1, L=L*ones(1,N);end
[u,s,~]=svd(A,0);
rk=sum(diag(s)>tol);
A=u(:,1:rk)'*A;
Y=u(:,1:rk)'*Y;
for k=1:N
	X(:,k)=SpLSSolver(A,Y(:,k),L,tol,sl,delta);
end
end
