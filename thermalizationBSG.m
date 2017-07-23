mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=15*msun;
e=1e51;
r=49*rsun;
Z=0.005;
a=7.566e-15;
X=0.7;
kappa_T=0.2*(1+X);

c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;

pconv=econv/(rconv^3);
rhoconv=mconv/(rconv^3);
vconv=sqrt(econv/mconv);

theta=linspace(pi/2,0,2048);
radius=linspace(0,2,2048)*rconv;
x=linspace(0,2,2048)*rconv;
y=linspace(0,2,2048)*rconv;


xp=zeros(2048,2048);
yp=zeros(2048,2048);
%theta=linspace(0,pi/2,2048);
for i=1:2048
    for j=1:2048
       xp(i,j)=radius(i)*sin(theta(j));
       yp(i,j)=radius(i)*cos(theta(j)); 
    end
end

xc=reshape(xp,2048*2048,1);
yc=reshape(yp,2048*2048,1);
[xg,yg] = meshgrid(x,y);
tau=zeros(2048,2048);


for t=28:28
    t
    name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
    name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
    velr= csvread(name).*vconv;
    name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
    velt= csvread(name).*vconv;    
    kappa_ff=3.7e22*(1+X)*(1-Z)*dens.*((3*pres/a).^-(7/8));
    kappa=kappa_T+kappa_ff;
    kappa_c=reshape(kappa,2048*2048,1);
    dens_c=reshape(dens,2048*2048,1);
    F_kappa = scatteredInterpolant(xc,yc,kappa_c);
    F_dens =scatteredInterpolant(xc,yc,dens_c);
    kappa_c=F_kappa(xg,yg);
    dens_c=F_dens(xg,yg);
    F=kappa_c.*dens_c*2/2048*rconv;
    for i=1:2048
       % xx=zeros(2048,2048-i);
       % xx=xg(:,i:2048)*rconv;
        tau(:,i)=trapz(F(:,i:2048),2);
    end
    name=[ '/home/nilou/Data/processeddata/BSG/tauBSG_' int2str(t-1) '.csv'];
    csvwrite(name,tau);
end