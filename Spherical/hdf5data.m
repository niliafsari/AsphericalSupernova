name='sedov_sphhdf5_chk_00';
cell=16384;
dens_data=zeros(cell,1);
eint_data=zeros(cell,1);
ener_data=zeros(cell,1);
pres_data=zeros(cell,1);
velx_data=zeros(cell,1);
vely_data=zeros(cell,1);
for i=143:267
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
    hinfo = hdf5info(['/home/nilou/bgq/output/' name]);
    coord = hdf5read(hinfo.GroupHierarchy.Datasets(3)); 
    dens = hdf5read(hinfo.GroupHierarchy.Datasets(4));   
    eint = hdf5read(hinfo.GroupHierarchy.Datasets(5));
    ener = hdf5read(hinfo.GroupHierarchy.Datasets(6));  
    pres = hdf5read(hinfo.GroupHierarchy.Datasets(19));
    velx = hdf5read(hinfo.GroupHierarchy.Datasets(29)); 
    vely = hdf5read(hinfo.GroupHierarchy.Datasets(30)); 
    dens=reshape(dens,[32,2,512]);    
    eint=reshape(eint,[32,2,512]);  
    ener=reshape(ener,[32,2,512]); 
    pres=reshape(pres,[32,2,512]);    
    velx=reshape(velx,[32,2,512]);  
    vely=reshape(vely,[32,2,512]); 

    for k=1:512
        dens_data((k-1)*32+1:k*32)=dens(1:32,1,k);
        eint_data((k-1)*32+1:k*32)=eint(1:32,1,k);
        ener_data((k-1)*32+1:k*32)=ener(1:32,1,k);
        pres_data((k-1)*32+1:k*32)=pres(1:32,1,k);
        velx_data((k-1)*32+1:k*32)=velx(1:32,1,k);
        vely_data((k-1)*32+1:k*32)=vely(1:32,1,k);
    end

    xlswrite(['/home/nilou/bgq/raw_data/dens' int2str(cell) '_' int2str(i) '.xlsx'],dens_data);
    xlswrite(['/home/nilou/bgq/raw_data/eint' int2str(cell) '_' int2str(i) '.xlsx'],eint_data);
    xlswrite(['/home/nilou/bgq/raw_data/ener' int2str(cell) '_' int2str(i) '.xlsx'],ener_data);
    xlswrite(['/home/nilou/bgq/raw_data/pres' int2str(cell) '_' int2str(i) '.xlsx'],pres_data);
    xlswrite(['/home/nilou/bgq/raw_data/velr' int2str(cell) '_' int2str(i) '.xlsx'],velx_data);
    xlswrite(['/home/nilou/bgq/raw_data/velth' int2str(cell) '_' int2str(i) '.xlsx'],vely_data);
end 
