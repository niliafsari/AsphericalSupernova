time1=zeros(142,1);
for i=0:266
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
    a=hdf5read(hinfo.GroupHierarchy.Datasets(22));
    time1(i+1,1)=(a(1,1).Data{2})
end
save('/home/nilou/bgq/output/timesteps.mat','time1');