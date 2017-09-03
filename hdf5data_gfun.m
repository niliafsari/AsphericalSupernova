name='sedov_sphhdf5_chk_00';
cell=1024;
dens_data=zeros(cell,cell);
gfun_data=zeros(cell,cell);
pres_data=zeros(cell,cell);
velx_data=zeros(cell,cell);
vely_data=zeros(cell,cell);
for i=175:175
    i
    if i<10
        name='sedov_sphhdf5_chk_000';
        name=[name int2str(i)];
    elseif i<100
        name='sedov_sphhdf5_chk_00';
        name=[name int2str(i)];  
    else
        name='sedov_sphhdf5_chk_0';
        name=[name int2str(i)];          
    end
    hinfo = hdf5info(['/home/nilou/gpc/' name]);
    coord = hdf5read(hinfo.GroupHierarchy.Datasets(3)); 
    dens = hdf5read(hinfo.GroupHierarchy.Datasets(4));   
    pres = hdf5read(hinfo.GroupHierarchy.Datasets(16));
    gfun = hdf5read(hinfo.GroupHierarchy.Datasets(9));
    velx = hdf5read(hinfo.GroupHierarchy.Datasets(26)); 
    vely = hdf5read(hinfo.GroupHierarchy.Datasets(27)); 
    dens=reshape(dens,[64,64,cell]);    
    pres=reshape(pres,[64,64,cell]);    
    velx=reshape(velx,[64,64,cell]);  
    vely=reshape(vely,[64,64,cell]); 
    gfun=reshape(gfun,[64,64,cell]); 
    for j=1:32
        ind=coord(2,((j-1)*32)+1);
        same=find(coord(2,:)==ind);
        %size(dset1(:,:,same'))
        for k=1:length(same)
            dens_data((k-1)*64+1:k*64,(j-1)*64+1:j*64)=dens(1:64,1:64,same(k));
            gfun_data((k-1)*64+1:k*64,(j-1)*64+1:j*64)=gfun(1:64,1:64,same(k));
            pres_data((k-1)*64+1:k*64,(j-1)*64+1:j*64)=pres(1:64,1:64,same(k));
            velx_data((k-1)*64+1:k*64,(j-1)*64+1:j*64)=velx(1:64,1:64,same(k));
            vely_data((k-1)*64+1:k*64,(j-1)*64+1:j*64)=vely(1:64,1:64,same(k));
        end
    end
    save(['/home/nilou/Data/rawdata/density/dens' int2str(cell) '_' int2str(i) '.mat'],'dens_data');
    save(['/home/nilou/Data/rawdata/pressure/pres' int2str(cell) '_' int2str(i) '.mat'],'pres_data');
    save(['/home/nilou/Data/rawdata/velocity/velr' int2str(cell) '_' int2str(i) '.mat'],'velx_data');
    save(['/home/nilou/Data/rawdata/gfunction/gfun' int2str(cell) '_' int2str(i) '.mat'],'gfun_data');
    save(['/home/nilou/Data/rawdata/velocity/velth' int2str(cell) '_' int2str(i) '.mat'],'vely_data');
end 
