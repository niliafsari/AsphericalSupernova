load('diff_bsg.mat','diff_bsg');

for i=1:50
    xdiff=diff_bsg{i}(1,:);
    ydiff=diff_bsg{i}(2,:);
    if sum(xdiff<0)
        i
        sum(xdiff<0)
    end
    if sum(ydiff<0)
        i
        sum(ydiff<0)
    end
    
end