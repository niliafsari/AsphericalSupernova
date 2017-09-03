mtot=4.82611e-3*2;
etot=1.3402e-2*2;
rtot=0.5;
msun=1.989e33;
rsun=6.955e10;

kappa=1;
c=1;
rconv=1;
econv=1;
mconv=1;

pconv=econv/rconv^3;
rhoconv=mconv/rconv^3;
vconv=sqrt(econv/mconv);

r=linspace(0,2,2048).*rconv;
theta=linspace(0,pi/2,2048);

cita=1;

if cita==1
    path='/mnt/scratch-lustre/nafsari/Data2048';
else
    path='/home/nilou/Data';
end

for t=161:170
    t
    d=zeros(2048,2048);
    if (t<=200 && t~=176)
        name=[path '/rawdata/pressure/pres1024_' int2str(t-1) '.csv'] ;
        pres= csvread(name).*pconv;
        name=[path '/rawdata/density/dens1024_' int2str(t-1) '.csv'] ;
        dens= csvread(name).*rhoconv;
        name=[path '/rawdata/velocity/velr1024_' int2str(t-1) '.csv'] ;
        velr= csvread(name).*vconv;
        name=[path '/rawdata/velocity/velth1024_' int2str(t-1) '.csv'] ;
        velt= csvread(name).*vconv;        
    else
        name=[path '/rawdata/pressure/pres1024_' int2str(t-1) '.mat'] ;
        load(name,'pres_data');
        pres=pres_data*pconv;
        name=[path '/rawdata/density/dens1024_' int2str(t-1) '.mat'] ;
        load(name,'dens_data');
        dens=dens_data*rhoconv;
        name=[path '/rawdata/velocity/velr1024_' int2str(t-1) '.mat'] ;
        load(name,'velx_data');
        velr=velx_data*vconv;
        name=[path '/rawdata/velocity/velth1024_' int2str(t-1) '.mat'] ;
        load(name,'vely_data');
        velt=vely_data*vconv;  
    end
    for i=1:2048
        for j=1:2048
            if (i==1)
                pr=r(i); 
                pindex_r=i;
            else
                pr=r(i-1);
                pindex_r=i-1;
            end
            if (j==1)
                pt=theta(j);
                pindex_t=j;
            else
                pt=theta(j-1);
                pindex_t=j-1;
            end
            if (i==2048)
                nr=r(i);
                nindex_r=i;
            else
                nr=r(i+1);
                nindex_r=i+1;
            end
            if (j==2048)
                nt=theta(j);
                nindex_t=j;
            else
                nt=theta(j+1);
                nindex_t=j+1;
            end
            dr=nr-pr;
            cr=pr+ dr/2;
            dt=nt-pt;
            ct=pt+ dt/2;
            delvel= ((nr^2 * velr(nindex_r,j) - pr^2 * velr(pindex_r,j))/ (cr^2 * dr) )+ (( sin(nt)*velt(i,nindex_t) - sin(pt) * velt(i,pindex_t))/(cr * sin(ct)* dt));
            gradp2=((pres(nindex_r,j)-pres(pindex_r,j))/dr)^2+ ((pres(i,nindex_t)-pres(i,pindex_t))/(cr*dt))^2;
            d(i,j)=(3*kappa*dens(i,j)*abs(delvel)*(pres(i,j)^2))/(c*gradp2);
        end
    end
    name=[ path '/rawdata/dcodeunit/dcodeunit1024_' int2str(t-1) '.mat'];
    save(name,'d');
    name=[ path '/rawdata/gradp2/gradp21024_' int2str(t-1) '.mat'];
    save(name,'gradp2');    
end