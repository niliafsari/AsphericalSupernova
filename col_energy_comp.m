mtot=0.0048187313;
etot=0.0130949665;
rtot=0.5;


ptot=etot/rtot^3;
rhotot=mtot/rtot^3;
vtot=sqrt(etot/mtot);
ttot=rtot/vtot;

theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048)/rtot;

xx=zeros(2048,2048);
yy=zeros(2048,2048);

thet=zeros(2048,2048);


for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
       thet(i,j)=theta(j);
       rad(i,j)=radius(i);
    end
end
time=load('/home/nilou/Data/timesteps.mat');
cone=[0 0 0 0 0 0 0.02 0.02 0.1 0.1 0.1 0.15 0.17 0.17 0.2 0.2 0.22 0.23 0.23 0.24 0.25 0.27 0.3 0.32 ...
    0.33 0.34 0.34 0.34 0.34 0.34 0.34 0.35 0.36 0.36 0.36 0.37 0.37 0.38 0.39  0.4 0.4 0.41 0.42 0.42 0.42];

for t=33:71
    t
    time_current=time.time1(t)/ttot;
    name=['/home/nilou/Data/rawdata/velocity/velr2048_' int2str(t-1) '.csv'] ;
    velr= csvread(name)/vtot;
    name=['/home/nilou/Data/rawdata/velocity/velth2048_' int2str(t-1) '.csv'] ;
    velt= csvread(name)/vtot;
    v_y=-velt.*sin(thet) + velr.*cos(thet);
    v_x=velt.*cos(thet) + velr.*sin(thet);
    t_y=(yy(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3))-((xx(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3))-1)*cone(t-26)/3))./-v_y(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3));
    r_x=1+v_x(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3)).*t_y;
    t_y(isinf(t_y))=median(t_y)+10^5;
    angleth=acos(1./r_x)*180/pi; 
    
    name=['/home/nilou/Data/rawdata/density/dens2048_' int2str(t-1) '.csv'] ;
    dens= csvread(name)/rhotot;
    dv=(rad.^2).*sin(thet)*2*pi*(4/2048)*(pi/(2*2048));
    e_z=(dens(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3)).*dv(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3))).*((v_y(v_y<0 & rad>1& yy>((xx-1)*cone(t-26)/3)).^2 +v_x(v_y<0 & rad>1 & yy>((xx-1)*cone(t-26)/3)).^2));    
    
    angleth_next=zeros(0,1);
    if (t<71)
        dt=(time.time1(t+1)/ttot-time_current);
        rnext=sqrt((xx+v_x*dt).^2 + (yy+v_y*dt).^2);
        rnext_c=rnext(v_y<0 &  rad>1 & yy>((xx-1)*cone(t-26)/3));

        t_colnext =t_y(t_y< dt | rnext_c>=4 )+time_current;
        r_xnext =r_x(t_y< dt | rnext_c>=4);
        e_znext =e_z(t_y< dt | rnext_c>=4);
        if (~isempty(r_xnext ))
            angleth_next =acos(1./r_xnext )*180/pi; 
        end
    else
        t_colnext =t_y+time_current;
        r_xnext =r_x;
        e_znext =e_z; 
        
        if (~isempty(r_xnext))
            angleth_next =acos(1./r_xnext )*180/pi; 
        end
    end
    size(e_znext )
    if (~isempty(r_xnext))
        save(['/home/nilou/Data/processeddata/Collision/collision_comp_v2_' num2str(t) '.mat'],'e_znext','t_colnext','angleth_next','r_xnext'); 
    end
    
    
    
    
    clear e_znext t_colnext angleth_next r_xnext
end


% grid off
% axis([0 90 -0.5 2])
% %caxis([log10(median(e_z)/((pi/4096)*(4/2048)))-3,log10(max(e_z)/((pi/4096)*(4/2048)))])
% caxis([-10 -7])
% colormap jet
% c=colorbar('Location','eastoutside');
% name=['/home/nilou/Data/plot/collision/col_comp_energy.png'];
% set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
% print(gcf, '-dpng', '-r900', name)
