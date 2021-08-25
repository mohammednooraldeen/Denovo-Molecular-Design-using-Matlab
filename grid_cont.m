
s=size(dddddd);c=1
for x=1:100:14000
    y=x+100;c=c-0.005;
 scatter3 (dddddd(x:y,1),dddddd(x:y,2),dddddd(x:y,3),50,[c 0 0],'fill'); hold on; 
end

scatter3(receptor(:,2),receptor(:,3),receptor(:,4),10,'fill'); hold off;