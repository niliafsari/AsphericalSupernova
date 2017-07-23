name='sedov_sphhdf5_chk_00';
cell=2048;
time_data=zeros(cell,cell);

for i=21:21
    i
    if i<10
        name='sedov_sphhdf5_chk_000';
        name=[name int2str(i)];
    else
        name='sedov_sphhdf5_chk_00';
        name=[name int2str(i)];    
    end
    hinfo = hdf5info(['/home/nilou/Data/' name]);
    coord = hdf5read(hinfo.GroupHierarchy.Datasets(3)); 
    mtime = hdf5read(hinfo.GroupHierarchy.Datasets(4));   
    mtime=reshape(mtime,[32,32,cell*2]);     
    for j=1:64
        ind=coord(2,j*64);
        same=find(coord(2,:)==ind);
        %size(dset1(:,:,same'))
        for k=1:length(same)
            time_data((k-1)*32+1:k*32,(j-1)*32+1:j*32)=mtime(1:32,1:32,same(k));
        end
    end
    csvwrite(['/home/nilou/Data/densmtime' int2str(cell) '_' int2str(i) '.csv'],time_data);
end