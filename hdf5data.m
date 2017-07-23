name='sedov_sphhdf5_chk_00';
cell=2048;
dens_data=zeros(cell,cell);
eint_data=zeros(cell,cell);
ener_data=zeros(cell,cell);
pres_data=zeros(cell,cell);
velx_data=zeros(cell,cell);
vely_data=zeros(cell,cell);
for i=0:75
    if i<10
        name='sedov_sphhdf5_chk_000';
        name=[name int2str(i)];
    else
        name='sedov_sphhdf5_chk_00';
        name=[name int2str(i)];    
    end
    hinfo = hdf5info(['/mnt/scratch-lustre/nafsari/Data2048/' name]);
    coord = hdf5read(hinfo.GroupHierarchy.Datasets(3)); 
    dens = hdf5read(hinfo.GroupHierarchy.Datasets(4));   
    eint = hdf5read(hinfo.GroupHierarchy.Datasets(5));
    ener = hdf5read(hinfo.GroupHierarchy.Datasets(6));  
    pres = hdf5read(hinfo.GroupHierarchy.Datasets(23));
    velx = hdf5read(hinfo.GroupHierarchy.Datasets(33)); 
    vely = hdf5read(hinfo.GroupHierarchy.Datasets(34)); 
    dens=reshape(dens,[32,32,cell*2]);    
    eint=reshape(eint,[32,32,cell*2]);  
    ener=reshape(ener,[32,32,cell*2]); 
    pres=reshape(pres,[32,32,cell*2]);    
    velx=reshape(velx,[32,32,cell*2]);  
    vely=reshape(vely,[32,32,cell*2]); 
    for j=1:64
        ind=coord(2,j*64);
        same=find(coord(2,:)==ind);
        %size(dset1(:,:,same'))
        for k=1:length(same)
            dens_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=dens(1:32,1:32,same(k));
            eint_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=eint(1:32,1:32,same(k));
            ener_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=ener(1:32,1:32,same(k));
            pres_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=pres(1:32,1:32,same(k));
            velx_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=velx(1:32,1:32,same(k));
            vely_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=vely(1:32,1:32,same(k));
        end
    end
    xlswrite(['dens' int2str(cell) '_' int2str(i) '.xlsx'],dens_data);
    xlswrite(['eint' int2str(cell) '_' int2str(i) '.xlsx'],eint_data);
    xlswrite(['ener' int2str(cell) '_' int2str(i) '.xlsx'],ener_data);
    xlswrite(['pres' int2str(cell) '_' int2str(i) '.xlsx'],pres_data);
    xlswrite(['velr' int2str(cell) '_' int2str(i) '.xlsx'],velx_data);
    xlswrite(['velth' int2str(cell) '_' int2str(i) '.xlsx'],vely_data);
end 
