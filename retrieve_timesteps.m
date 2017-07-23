time1=zeros(75,1);
for i=0:75
    if i<10
        name='sedov_sphhdf5_chk_000';
        name=[name int2str(i)];
    else
        name='sedov_sphhdf5_chk_00';
        name=[name int2str(i)];    
    end
    hinfo = hdf5info(['/home/nilou/Data/' name]);
    a=hdf5read(hinfo.GroupHierarchy.Datasets(26));
    time1(i+1,1)=(a(1,1).Data{2})
end
save('timesteps.mat',time1);