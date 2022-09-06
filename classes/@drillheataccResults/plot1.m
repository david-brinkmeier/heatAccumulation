function plot1(obj)
% call results.plot1(results)
% get heatsource lenght / drilldepths
hslens = cellfun(@max,obj.sim.Heatsource2Dlen); % accounts for weber edgecase

fig = genORselectfigbyname('Overview [metrics]',obj.UUID);
fig.Color = [1 1 1]; fig.Position = [810 600 560 420];
plot([1 obj.sim.TotalNoPulses],1e6*[obj.drilldepth_max obj.drilldepth_max],'--k','DisplayName','Förster')
hold on
plot(1:obj.sim.TotalNoPulses,1e6*hslens,'-k','DisplayName','Weber')
ylim([0 round(1.05*1e6*obj.drilldepth_max,-1)])
title({strcat('P_{avg,max}=',32,num2str(round(obj.MeanmaxavgLaserPower{2},2)),'W,',32,'Ep_{max}=',32,num2str(round(obj.laser.Ep*1e6)),'µJ,',32,...
    'frep_{max}=',32,num2str(max(obj.sim.reprates{1}*1e-3)),'kHz,',32,'w0 =',32,num2str(obj.laser.w0*1e6),'µm'),...
    strcat('P_{avg}=',32,num2str(round(obj.MeanmaxavgLaserPower{1},2)),'W,',32,'Q_{in,tot}=',32,num2str(round(obj.Qtotal_in,2)),'J,',32,...
    '"drill_{efficiency}"=',32,num2str(round(obj.MillimeterPerJoule,2)),32,'mm/J')})
xlabel('Number of pulses','FontSize',12);
ylabel('Depth in µm','FontSize',12);
yyaxis right
plot(1:obj.sim.TotalNoPulses,obj.gouffe,'--r','DisplayName','Gouffé')
hold on
plot(1:obj.sim.TotalNoPulses,obj.sim.EnergyRampFull{1},'-r','DisplayName','Ep')
plot(1:obj.sim.TotalNoPulses,obj.sim.reprates{1}./max(obj.sim.reprates{1}(:)),'-b','DisplayName','frep')
plot(1:obj.sim.TotalNoPulses,obj.AvgLaserPower./obj.MeanmaxavgLaserPower{2},'-.','Color','#0072BD','DisplayName','Pavg')
ylabel('Gouffé / Energy / Power / Reprate [0-1]')
legend('Location','Southeast','NumColumns',2,'AutoUpdate','off')
xlim([1 obj.sim.TotalNoPulses]),
ylim([min([obj.gouffe(:); obj.sim.EnergyRampFull{1}(:)]), ...
    round(max([obj.gouffe(:); obj.sim.EnergyRampFull{1}(:)])*1.05,2)])

% make y axes black
ax1 = gca;
ax1.LineWidth = 1;
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = 'k';
drawnow
end