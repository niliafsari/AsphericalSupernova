
theta=linspace(0,pi/2,2048);
radius=linspace(0,2,2048);

xx=zeros(2048,2048);
yy=zeros(2048,2048);

for i=1:2048
    for j=1:2048
       xx(i,j)=radius(i)*sin(theta(j));
       yy(i,j)=radius(i)*cos(theta(j)); 
    end
end
    a=get(gcf,'Position');
    close all

     myFigure = figure ;    
     set(myFigure,'units','points','position',a,'Color','white')  
surf(xx,yy, log10(dens_data));    
axis equal
grid off
view(2)