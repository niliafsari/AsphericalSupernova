a=7.566e-15;
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

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);
tconv=rconv/vconv;

r=0.5*rconv;

radius=linspace(0,8,16384)*rconv;
diffrsg=zeros(200,1);
time=load('/home/nilou/bgq/output/timesteps.mat')
time_rsg=time.time1*tconv;
diff=load('/home/nilou/bgq/processeddata/BSG/diffrsg.mat');
diff_front=diff.diffbsg;
ltot=zeros(200,1);
t_loc_all=zeros(267,1);
for t=1:200
    name=['/home/nilou/bgq/raw_data/pres16384_' int2str(t-1) '.csv'] ;
    pres= csvread(name).*pconv;
    name=['/home/nilou/bgq/raw_data/dens16384_' int2str(t-1) '.csv'] ;
    dens= csvread(name).*rhoconv;
    name=['/home/nilou/bgq/raw_data/velr16384_' int2str(t-1) '.csv'] ;
    velr= csvread(name).*vconv;
    gradp=zeros(16384,1);
    dr=8*rconv/16384;
    for i=1:16384 
        if (i==16384)
            delvel=(velr(i)-velr(i-1))/dr ;
            gradp(i)=(pres(i)-pres(i-1))/dr;
        else
            delvel=(velr(i+1)-velr(i))/dr ;
            gradp(i)=(pres(i+1)-pres(i))/dr;
        end
    end
    [u,i]=min(abs(radius-diff_front(t)));
    ltot(t)=-4*pi*diff_front(t)^2*c*gradp(i)/(kappa*dens(i));
    t_loc_all(t)=(3*pres(i)/a).^0.25;
    tdiff_all=(3*kappa*dens(i)).*(pres(i).^2)./(c*gradp(i).^2);
    etha_all(t)=(7e5./min(time_rsg(t),tdiff_all)).*((dens(i)/1e-10).^-2).*(t_loc_all(t)/(100*11604.52)).^3.5;
    
end
save('/home/nilou/Data/processeddata/BSG/colortemp_BSG_spherical.mat','t_loc_all','time_rsg','etha_all');
