load allOutputsDecrease;
figure
subplot(121)
ref = out1323000; % reference is 30x samplerate
upsamp = length(ref) / length(out44100);
range = [1:10,15, 20, 25];
error = zeros(length(range), 4410);
j = 1;
for i = range
    eval("error(j, :) = sqrt((out" + num2str(i * 44100)+"(1:"+i+":end) - ref(1:"+ upsamp +":end)).^2);")
    j = j + 1;
end
loglog(range, sum(error, 2)'/length(out44100), 'k', 'Linewidth', 2)
xLab = xlabel('$f_s$ (in Hz) (x44100)', 'interpreter', 'latex');
yLab = ylabel('MSE for 0.1s', 'FontName', 'times');
title("Decreasing $\mathcal{N}$", 'interpreter', 'latex');
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', [1, 5:5:25], ...
    'xticklabel', [1, 5:5:25], 'TickLabelInterpreter', 'latex', 'Position', [0.0949 0.1776 0.3931 0.7240]);
set(gcf, 'Position', [477 582 592 275]);
grid on
ylim([0.0004, 0.06])
xLab.Position = [5.0000 2.3932e-04 -1];
yLab.Position = [0.5810 0.0049 -1.0000];
%%
subplot(122)
load allOutputsIncrease
ref = outTest1323000; % reference is 10x samplerate
upsamp = length(ref) / length(outTest44100);
range = [1:10, 15, 20, 25];
error = zeros(length(range), 4410);
j = 1;
for i = range
    eval("error(j, :) = sqrt((outTest" + num2str(i * 44100)+"(1:"+i+":end) - ref(1:"+ upsamp +":end)).^2);")
    j = j + 1;
end
loglog(range, sum(error, 2)'/length(out44100), 'k', 'Linewidth', 2)
xLab = xlabel('$f_s$ (in Hz) (x44100)', 'interpreter', 'latex');
% ylabel('MSE for 0.1s', 'FontName', 'times');
title("Increasing $\mathcal{N}$", 'interpreter', 'latex');
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', [1, 5:5:25], ...
    'xticklabel', [1, 5:5:25], 'TickLabelInterpreter', 'latex', 'Position', [0.5649 0.1776 0.3931 0.7240]);
grid on
xLab.Position = [5.0742 2.4545e-04 -1];
% yLab.Position = [0.611,0.0037,-1];
ylim([0.0004, 0.06])

%%
% figure;
% subplot(211);
% plot(out44100, '-.', 'Linewidth', 1, 'color', 'r');
% hold on
% plot(out220500(1:5:end), '--', 'Linewidth', 0.75, 'color', 'b')
% plot(out441000(1:10:end), 'Linewidth', 0.5, 'color', 'g')
% legend({'$f_s = 44100$', '$f_s = 220500$', '$f_s = 441000$'}, 'interpreter', 'latex')
% subplot(212);
% plot(out44100, '-.', 'Linewidth', 1, 'color', 'r');
% hold on
% plot(out220500(1:5:end), '--', 'Linewidth', 0.75, 'color', 'b')
% plot(out441000(1:10:end), 'Linewidth', 0.5, 'color', 'g')
% legend({'$f_s = 44100$', '$f_s = 220500$', '$f_s = 441000$'}, 'interpreter', 'latex')
% xlim([(length(out44100) - 100), length(out44100)])



%% Create data
% for fs = (44100*20):44100*5:(44100*30)
%     outTest = convergenceCheck(fs, false);
%     eval("outTest" + num2str(fs)+"= outTest;");
% end
