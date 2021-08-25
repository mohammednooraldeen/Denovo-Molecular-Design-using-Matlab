x= 10; scatter3 (ddd(1:x:end,1),ddd(1:x:end,2),ddd(1:x:end,3),20,data_norm(1:x:end,1:3)); hold on; scatter3(receptor(:,2),receptor(:,3),receptor(:,4),10,'fill'); hold off


x= 10; scatter3 (ddd(1:x:end,1),ddd(1:x:end,2),ddd(1:x:end,3),20,[1 0 0]); hold on; scatter3(receptor(:,2),receptor(:,3),receptor(:,4),10,'fill'); hold off
