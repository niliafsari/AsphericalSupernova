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
radius=linspace(0,2,2048);
x=linspace(0,2,2048);
y=linspace(0,2,2048);


xp=zeros(2048,2048);
yp=zeros(2048,2048);
%theta=linspace(0,pi/2,2048);
for i=1:2048
    for j=1:2048
       xp(i,j)=radius(i)*sin(theta(j));
       yp(i,j)=radius(i)*cos(theta(j)); 
    end
end

radius=zeros(1,2048);

for t=1:50        
    name1=['/home/nilou/Data/rawdata/dcodeunit/dcodeunit_' int2str(t-1) '.csv'] ;
    d=csvread(name1);
    indexes=find(isnan(d) | isinf(d));
    indexes_tot=indexes((xp(indexes).^2+yp(indexes).^2)>0.5^2);
    name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
    dens(indexes_tot)=0;
    name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
    velr= csvread(name).*vconv;
    name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
    velt= csvread(name).*vconv;    
    kappa_ff=3.7e22*(1+X)*(1-Z)*dens.*((3*pres/a).^-(7/8));
    kappa=kappa_T+kappa_ff;
    dtau_ff=kappa_ff.*dens*2/2048*rconv;
    dtau=kappa.*dens*2/2048*rconv;
    for i=1:2048
        tau(i,:)=trapz(dtau(i:2048,:),1);
        tau_ff(i,:)=trapz(dtau_ff(i:2048,:),1); 
    end
    name=[ '/home/nilou/Data/processeddata/BSG/tauBSGradial_' int2str(t-1) '.csv'];
    nameff=[ '/home/nilou/Data/processeddata/BSG/tauBSGradialFF_' int2str(t-1) '.csv'];
    csvwrite(name,tau);
    csvwrite(nameff,tau_ff);
end