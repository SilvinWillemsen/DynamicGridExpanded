clear all;
close all;
clc;

membraneStates = load("/Users/SilvinW/repositories/ModularVST/Builds/MacOSX/build/Debug/statesSaveMembrane.csv");
NxNy = load("/Users/SilvinW/repositories/ModularVST/Builds/MacOSX/build/Debug/NxNy.csv");
NxC = NxNy(1) - 3; % Nx clamped
NyC = NxNy(2) - 3; % Ny clamped

lengthSound = length(membraneStates / NyC);
for n = 1:lengthSound
    imagesc(reshape(membraneStates(((n-1)*NyC + 1):(n*NyC), :), NxC, NyC))
    pause(0.5);
    max(max(abs(membraneStates(((n-1)*NyC + 1):(n*NyC), :))))
    drawnow;
end