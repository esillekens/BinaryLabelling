clear
close all

folder = 'Binary-Labeling-for-2D-and-4D-constellations/coordinates&labeling2D/';
% mod = 'cross2_64';
% mod = 'c2_64';
% mod = 'honeycomb2_64';
% mod = 'DSQ2_32';

mod = '2D-GS-64';

% mod = 'opt12dB_64';

% mod = 'rand2_64';
% mod = 'randn2_64';

if strncmp(mod,'opt',3)
    matfile = load(fullfile(folder,mod,[mod,'.mat']));
    X = matfile.X;
    label = bi2de(matfile.Lbin);
elseif strncmp(mod,'2D-GS',4)
    matfile = load(fullfile(folder,mod,[mod,'.mat']));
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

SNR = 0:20;
figure, hold on, legend('show')

plot(SNR,gmifun(SNR,X),...
    'DisplayName','Original')

if ~isempty(label)
    plot(SNR,gmifun(SNR,relabel(X,label)),...
        'DisplayName','BSA')
end

plot(SNR,gmifun(SNR,relabel(X,SpectralLabel(X,1))),...
    'DisplayName','Spectral')


plot(SNR,gmifun(SNR,RelabelAPSK(X)),...
    'DisplayName','APSK')


plot(SNR,gmifun(SNR,RelabelAPSK_spectral(X)),...
    'DisplayName','APSK_{spectral}')


plot(SNR,gmifun(SNR,RelabelSquare(X)),...
    'DisplayName','QAM')

