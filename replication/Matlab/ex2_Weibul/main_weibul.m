%  RE-parameterize
clear all; clc
rng default;

diary off
lbl1=mfilename;
logfile_name=strcat(lbl1,'-Diary.log');
delete(logfile_name);
diary(logfile_name);


rho=0.0; %must be positive in program
r=sqrt(rho/(1-rho));
k=3;
n=10^7;
% x1, x2, x3 are teh x's in the different time periods
x1=randn(n,k);
a=randn(n,1);
x1(:,2)=(x1(:,2)+a*r)/sqrt(1+r*r);
x1(:,3)=(x1(:,3)+a*r)/sqrt(1+r*r);
x1(:,1)=ones(n,1);

x2=x1;
x2(:,3)=(x2(:,3)+randn(n,1))*sqrt(0.5);
x3=x2;
x3(:,3)=(x3(:,3)+randn(n,1))*sqrt(0.5);

beta=ones(k,1)*sqrt(0.5);
beta(1)=-1;
alpha=2.0;
eta=exp(randn(n,1)*beta(2));

E=-log(rand(n,1));


censor=300+zeros(n,1);
% converged=zeros(n,1);

t1=zeros(n,1);
t2=censor;
Z1=zeros(n,1);
[Z2]=fun.Z(t2,x1,x2,x3,eta,[beta;alpha]);
if length(find(E>Z2))>0
    'censoring'
    stop
end



toler=1.0d-12;
iter=0;

while max(abs(t2-t1)) > toler
    iter=iter+1;
    tm=(t1+t2)*0.5;
    [Zm]=fun.Z(tm,x1,x2,x3,eta,[beta;alpha]);
    ii=(Zm>E);
    t2=ii.*tm+(1-ii).*t2;
    t1=ii.*t1+(1-ii).*tm;
end

T=(t1+t2)/2;

fun.print([mean(T<1) mean((T<2).*(T>1)) mean((T>2))]) 


% stop

[ZZ]=fun.Z(T,x1,x2,x3,eta,[beta;alpha]);

%  REPARAMETERIZING 
theta=[beta/alpha;alpha];

mean(log(ZZ)+double(eulergamma))
mean(log(ZZ)+double(eulergamma))/(std(log(ZZ))/sqrt(n))

[ZZ0]=fun.Z(T,x1,x2,x3,1+0*eta,[beta;alpha]);

mean(log(ZZ0)+double(eulergamma))
mean(log(ZZ0)+double(eulergamma))/(std(log(ZZ0))/sqrt(n))



[m,S]=fun.MoM(theta,T,x1,x2,x3);




% calculate the analytic S

V1=pi^2/6+beta(2)^2;
Sanalytic=V1*[1 0 0 0 0; ...
           0 1 0 0 0; ...
           0 0 1 sqrt(1/2) 1/2; ...
           0 0 sqrt(1/2) 1 sqrt(1/2); ...
           0 0 1/2 sqrt(1/2) 1];



kk=length(m);
% W=inv(Sanalytic);
W=inv(S);
[SE,grad,sens,Avar] = fun.sensitivity(theta,T,x1,x2,x3,S,W);

diary on

fprintf('OPTIMAL\n\n')

b=sqrt(diag(Avar));
fun.print(Avar./(b*b'))

lbl={ '$\beta_0$' '$\beta_1$' '$\beta_2$' '$\alpha$'};

hdr={'$E\left[ e \right] $' 
    '$E\left[ e x_{11}\right] $' 
    '$E\left[ e x_{21}\right] $' 
    '$E\left[ e x_{22}\right] $ '
    '$E\left[ e x_{23}\right] $ '
};


meas_str = {'M_1','\mathcal{E}_2','\mathcal{E}_3','\mathcal{E}_4','\mathcal{E}_5','\mathcal{E}_6'};
meas = {sens.M1,sens.M2e,sens.M3e,sens.M4e,sens.M5e,sens.M6e};
fun.write_res(meas,meas_str,lbl)


W=inv(diag(diag(S)));

[SE,grad,sens,Avar] = fun.sensitivity(theta,T,x1,x2,x3,S,W);

diary on

fprintf('NON-OPTIMAL\n\n')

b=sqrt(diag(Avar));
fun.print(Avar./(b*b'))


meas = {sens.M1,sens.M2e,sens.M3e,sens.M4e,sens.M5e,sens.M6e};
fun.write_res(meas,meas_str,lbl)

    
    
diary off
