uSave = load('/Users/SilvinW/repositories/RealTimeDynamic/Builds/MacOSX/DerivedData/uSaveDynamic.csv');

for n = 1:size(uSave, 1)
    hold off
    plot (0:length(uSave(n,:))-3, uSave(n, 1:end-2))
    hold on;
    plot (length(uSave(n,:))-3:length(uSave(n,:))-2, uSave(n, end-1:end))

    pause(0.05)
    drawnow;
end