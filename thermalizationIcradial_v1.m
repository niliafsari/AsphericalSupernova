mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
c=3e10;
m_bsg=15*msun;
e_bsg=1e51;
r_bsg=49*rsun;
rconv_bsg=r_bsg/rtot;
econv_bsg=e_bsg/etot;
mconv_bsg=m_bsg/mtot;
pconv_bsg=econv_bsg/(rconv_bsg^3);
rhoconv_bsg=mconv_bsg/(rconv_bsg^3);
vconv_bsg=sqrt(econv_bsg/mconv_bsg);

m=5*msun;
e=1e51;
r=0.2*rsun;
Z=0.005;
a=7.566e-15;
X=0;
kappa_T=0.2*(1+X);
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

for t=2:50      
    t
    name=[ '/home/nilou/Data/processeddata/BSG/gradp_' int2str(t-2) '.csv'];
    gradp=csvread(name)*(rconv_bsg/pconv_bsg)*pconv/rconv;
%     name1=['/home/nilou/Data/rawdata/dcodeunit/dcodeunit_' int2str(t-1) '.csv'] ;
%     d=csvread(name1);
%     indexes=find(isnan(d) | isinf(d));
%     indexes_tot=indexes((xp(indexes).^2+yp(indexes).^2)>0.5^2);
    name=['/home/nilou/Data/rawdata/pressure/pres2048_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
%     dens(indexes_tot)=0;
%     name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
%     velr= csvread(name).*vconv;
%     name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
%     velt= csvread(name).*vconv;    
    kappa_ff=3.7e22*(1+X)*(1-Z)*dens.*((3*pres/a).^-(7/8));
    kappa=kappa_T+kappa_ff;
    tau_ff=kappa_ff.*dens.*pres./gradp;
    tau=kappa.*dens.*pres./gradp;
    name=[ '/home/nilou/Data/processeddata/ic/tauIcradial_v2_' int2str(t-1) '.csv'];
    nameff=[ '/home/nilou/Data/processeddata/ic/tauIcradialFF_v2_' int2str(t-1) '.csv'];
    csvwrite(name,tau);
    csvwrite(nameff,tau_ff);
end