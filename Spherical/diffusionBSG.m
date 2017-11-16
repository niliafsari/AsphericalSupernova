mtot=0.0048187313*2;
etot=9.1997e-4*2;
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

pconv=econv/(rconv^3);
rhoconv=mconv/(rconv^3);
vconv=sqrt(econv/mconv);

r=linspace(0,8,16384).*rconv;

for t=1:200
    name=['/home/nilou/bgq/raw_data/pres16384_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/bgq/raw_data/dens16384_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
    name=['/home/nilou/bgq/raw_data/velr16384_' int2str(t-1) '.csv'] ;
    velr= csvread(name).*vconv;
    d=zeros(16384,1);
    dr=8*rconv/16384;
    for i=1:16384 
        if (i==16384)
            delvel=(velr(i)-velr(i-1))/dr ;
            gradp2=((pres(i)-pres(i-1))/dr)^2;
        else
            delvel=(velr(i+1)-velr(i))/dr ;
            gradp2=((pres(i+1)-pres(i))/dr)^2;
        end
        d(i,1)=(3*kappa*dens(i)*abs(delvel)*(pres(i)^2))/(c*gradp2);
    end
    name=[ '/home/nilou/bgq/processeddata/BSG/dparamBSG_' int2str(t-1) '.csv'];
    csvwrite(name,d);
end