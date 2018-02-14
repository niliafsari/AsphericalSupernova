name='sedov_sphhdf5_chk_00';
cell=256;
dens_data=zeros(cell,cell);
eint_data=zeros(cell,cell);
ener_data=zeros(cell,cell);
pres_data=zeros(cell,cell);
velx_data=zeros(cell,cell);
vely_data=zeros(cell,cell);
for i=216:216
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
    hinfo = hdf5info(['/home/nilou/gpc256/' name]);
    coord = hdf5read(hinfo.GroupHierarchy.Datasets(3)); 
    dens = hdf5read(hinfo.GroupHierarchy.Datasets(4));   
    %eint = hdf5read(hinfo.GroupHierarchy.Datasets(5));
    %ener = hdf5read(hinfo.GroupHierarchy.Datasets(6));  
    pres = hdf5read(hinfo.GroupHierarchy.Datasets(16));
    %velx = hdf5read(hinfo.GroupHierarchy.Datasets(33)); 
    %vely = hdf5read(hinfo.GroupHierarchy.Datasets(34)); 
    dens=reshape(dens,[cell,cell]);    
  
    pres=reshape(pres,[cell,cell]);    


    xlswrite(['dens' int2str(cell) '_' int2str(i) '.xlsx'],dens);
    %xlswrite(['eint' int2str(cell) '_' int2str(i) '.xlsx'],eint_data);
    %xlswrite(['ener' int2str(cell) '_' int2str(i) '.xlsx'],ener_data);
    xlswrite(['pres' int2str(cell) '_' int2str(i) '.xlsx'],pres);
    %xlswrite(['velr' int2str(cell) '_' int2str(i) '.xlsx'],velx_data);
    %xlswrite(['velth' int2str(cell) '_' int2str(i) '.xlsx'],vely_data);
end 
