function gene_out = hybrid_simul( gene,receptor_map,lb,ub )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

 [x,feval]=simulannealbnd(@(A) simulated_pos_quat(A,gene,'receptor.C.map'),...
     [45.686 46 18.443 ],lb(1:3),ub(1:3), saoptimset('InitialTemperature',[50 50 50],...
     'TemperatureFcn',@temperaturefast,'PlotFcns',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping}))
gene_out=gene;
gene_out{2,1}=x;

end

