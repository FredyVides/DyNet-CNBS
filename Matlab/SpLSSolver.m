% Code by: Fredy Vides
% For paper, Sparse Autoregressive Reservoir Computers for 
% Dynamic Financial Processes Identification
%%
%%
function x=SpLSSolver(A,y,L,tol,sl,delta)
	if nargin<=4 
		sl=1;delta=tol;
	elseif nargin<=5 
		delta=tol;
	end
N=size(A,2);
K=1;
Error = abs(1+tol);
w=sparse(N,1);
x=w;
[u,s,v]=svd(A,0);
rk=sum(diag(s)>tol);
u=u(:,1:rk);
s=s(1:rk,1:rk);
v=v(:,1:rk);
c=v*(s\(u'*y));
ac=abs(c);
[~,f]=sort(-ac);
L=min(max(sum(ac(f)>delta),1),L);
while K<=L && Error>tol
	ff=sort(f(1:K));
	x=w;
	[u,s,v]=svd(A(:,ff),0);
	rk=sum(diag(s)>tol);
	u=u(:,1:rk);
	s=s(1:rk,1:rk);
	v=v(:,1:rk);
	c=v*(s\((u'*y)));
	x(ff,1)=c;
	Error = norm(A*x-y);
	K=K+sl;
end
end
