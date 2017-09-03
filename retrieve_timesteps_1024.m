time1=zeros(550,1);
for i=0:550
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
    a=hdf5read(hinfo.GroupHierarchy.Datasets(19));
    time1(i+1,1)=(a(1,1).Data{2});
end
save('/home/nilou/Data/processeddata/timesteps_1024.mat','time1');