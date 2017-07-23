mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;
m=15*msun;
e=1e51;
r=49*rsun;
kappa=0.34;
c=3e10;
rconv=r/rtot;
econv=e/etot;
mconv=m/mtot;


pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;
times=(time.time1*tconv);
load('luminosity.mat', 'luminosity')

luminosity=luminosity*kappa*rhoconv/(pconv*rconv);

save('luminosityCode.mat', 'luminosity')