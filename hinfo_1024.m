name='sedov_sphhdf5_chk_00';
cell=1024;
dens_data=zeros(cell,cell);

pres_data=zeros(cell,cell);

for i=24:24
    if i<10
        name='sedov_sphhdf5_chk_000';
        name=[name int2str(i)];
    else
        name='sedov_sphhdf5_chk_00';
        name=[name int2str(i)];    
    end
    hinfo = hdf5info(['/home/nilou/cita/nxb1024nyb1024v2/' name]);
    coord = hdf5read(hinfo.GroupHierarchy.Datasets(3)); 
    dens = hdf5read(hinfo.GroupHierarchy.Datasets(4));   
  
    pres = hdf5read(hinfo.GroupHierarchy.Datasets(23));
 
    dens=reshape(dens,[32,32,cell]);    
    pres=reshape(pres,[32,32,cell]);    

    for j=1:32
        ind=coord(2,j*32);
        same=find(coord(2,:)==ind);
        %size(dset1(:,:,same'))
        for k=1:length(same)
            dens_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=dens(1:32,1:32,same(k));
            pres_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=pres(1:32,1:32,same(k));
        end
    end
    xlswrite(['dens' int2str(cell) '_' int2str(i) '.xlsx'],dens_data);
 
    xlswrite(['pres' int2str(cell) '_' int2str(i) '.xlsx'],pres_data);

end 
