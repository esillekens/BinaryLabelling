clear
close all

folder = 'Binary-Labeling-for-2D-and-4D-constellations/coordinates&labeling4D/';
% mod = 'b4_64';
% mod = 'l4_64';
% mod = 'w4_256';
mod = 'a4_256';

% mod = 'SO-PM-QPSK4_16';

% mod = '7b2A8PSK';

% mod = '4D-32';
% mod = '4D-64PRS';
% mod = '4D-OS128';
% 

% 
% mod = 'rand4_64';
% mod = 'randn4_64';

if strncmp(mod,'4D-',3) || strcmp(mod(2:end),'b2A8PSK')
    fist = dir(fullfile(folder,mod,'*.mat'));
    matfile = load(fullfile(folder,mod,fist(1).name));
    X = matfile.X;
    label = bi2de(matfile.L);
elseif strncmp(mod,'rand',4)
    shape = fliplr(str2double(strsplit(mod(5:end),'_')));
    X = sqrt(3)*(2*rand(shape)-1);
    label = [];
elseif strncmp(mod,'randn',5)
    shape = fliplr(str2double(strsplit(mod(6:end),'_')));
    X = randn(shape);
    label = [];
else
    X = readmatrix(fullfile(folder,mod,[mod,'_X.txt']));
    label = getfield(load(fullfile(folder,mod,[mod,'_Lbin.mat'])),'best_L_BSA_full');
end
gmifun = @(snr,x) arrayfun(@(s)GMI_withgrad_nD(s,x),snr );

X = X./sqrt(mean(sum(abs(X).^2,2)));
M = size(X,1);

SNR = 4:2:20;
figure, hold on, legend('show')

plot(SNR,gmifun(SNR,X),...
    'DisplayName','Original')

if ~isempty(label)
    plot(SNR,gmifun(SNR,relabel(X,label)),...
        'DisplayName','BSA')
end

plot(SNR,gmifun(SNR,relabel(X,SpectralLabel(X,1))),...
    'DisplayName','Spectral')


plot(SNR,gmifun(SNR,RelabelSquare(X)),...
    'DisplayName','QAM')

