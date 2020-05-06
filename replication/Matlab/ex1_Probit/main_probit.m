clear all; clc
rng default;

diary off
lbl1=mfilename;
logfile_name=strcat(lbl1,'-Diary.log');
delete(logfile_name);
diary(logfile_name);
diary on

rho=0.5; %must be positive in program
r=sqrt(rho/(1-rho));
k=3;
n=10000000;
x=randn(n,k);
a=randn(n,1);
x(:,2)=(x(:,2)+a*r)/sqrt(1+r*r);
x(:,3)=(x(:,3)+a*r)/sqrt(1+r*r);
x(:,1)=ones(n,1);

beta=ones(k,1)/sqrt(2+2*rho);

y=(x*beta+randn(n,1)>0);

[m,S]=fun.MoM(beta,y,x);


kk=k*(k+1)/2;


fprintf('Optimal')

W=inv(S);
[SE,grad,sens,Avar] = fun.sensitivity(beta,y,x,S,W);

lbl={'$\beta_0$' '$\beta_1$' '$\beta_1$' };

hdr={'$E\left[ e \right] $' 
    '$E\left[ e x_{1}\right] $' 
    '$E\left[ e x_{2}\right] $ '
    '$E\left[ e x_{1}^{2} \right] $'
    '$E\left[ e x_{1}x_{2}\right] $' 
    '$E\left[ e x_{2}^{2}\right] $'};

titl='Method of Moments Estimation of Probit Model. M1';
fun.write_res(sens.M1,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M2e';
fun.write_res(sens.M2e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M3e';
fun.write_res(sens.M3e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M4e';
fun.write_res(sens.M4e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M5e';
fun.write_res(sens.M5e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M6e';
fun.write_res(sens.M6e,titl,lbl,hdr)



fprintf('Non-Optimal')

W=inv(diag(diag(S)));
[SE,grad,sens,Avar] = fun.sensitivity(beta,y,x,S,W);

lbl={'$\beta_0$' '$\beta_1$' '$\beta_1$' };

hdr={'$E\left[ e \right] $' 
    '$E\left[ e x_{1}\right] $' 
    '$E\left[ e x_{2}\right] $ '
    '$E\left[ e x_{1}^{2} \right] $'
    '$E\left[ e x_{1}x_{2}\right] $' 
    '$E\left[ e x_{2}^{2}\right] $'};

titl='Method of Moments Estimation of Probit Model. M1';
fun.write_res(sens.M1,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M2e';
fun.write_res(sens.M2e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M3e';
fun.write_res(sens.M3e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M4e';
fun.write_res(sens.M4e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M5e';
fun.write_res(sens.M5e,titl,lbl,hdr)

titl='Method of Moments Estimation of Probit Model. M6e';
fun.write_res(sens.M6e,titl,lbl,hdr)

diary off



    
    

