mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

m=5*msun;
e=1e51;
r=0.2*rsun;

kappa=0.2;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;


theta=linspace(0,pi/2,2048);


radius=linspace(0,2,2048)*rconv;
r=linspace(0,2,2048).*rconv;
for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
gradp=zeros(2048,2048);

cita=0;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end


time=load([path '/processeddata/timesteps_1024.mat']);
syms x
f=@(x) 1./(x.^(3/2).*(1-x.^-4))

for t=200:200
    load([path '/processeddata/ic/x0_' num2str(t-1) '.mat'], 'x0')
    for k=2:length(x0)
        y = logspace(-5,4,50);
        q=zeros(1,50);
        for i=1:length(y)
            q(i) = integral(f,1+y(i),x0(k))
        end
        g= fit(q',log10(y)','poly7');
        save([path '/processeddata/ic/func_ic_' num2str(t-1) '.mat'], 'g')
    end
end
