%%
clear all;

file = load('AdhesiveNail.mat');
AN = file.TDM;
AdhesiveNail = AN(:,15);

file = load('AdhesiveNonObject.mat');
AdhesiveOnly = file.non_obj;

file = load('Air.mat');
Air = file.mid_air;

file = load('FoamAdhesive.mat');
FAdhesive = file.FDM;


file = load('FoamAir.mat');
FAir = file.non_obj;


file = load('FoamNonObj.mat');
FoamOnly = file.non_obj;


file = load('FoamNail.mat');
FN = file.FNM;
FoamNail = FN(:,15);

NailOnly = FoamNail - FoamOnly;

%% Original Data

figure;
TOTAL = 4;
subplot(TOTAL,1,1), 
plot(transpose(db(AN(:,13))) - min(transpose(db(AN(:,13))))), title('Signal of AdhesiveNail');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');
subplot(TOTAL,1,2), plot(db(Air) - min(db(Air))), title('Signal of Air');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')
subplot(TOTAL,1,3), plot(db(AdhesiveOnly) - min(db(AdhesiveOnly))),  title('Signal of AdhesiveOnly');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')
subplot(TOTAL,1,3), plot(db(AdhesiveOnly) - min(db(AdhesiveOnly))),  title('Signal of AdhesiveOnly');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')
subplot(TOTAL,1,4),
plot(db(FoamOnly)),  title('Nail');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')

%% Independent Components Analysis using FastICA
% addpath('/Users/qizhu/Documents/MATLAB/Matlab Code-2/FastICA_25');
PCNo = 4;
ICNo = 4;

% 12:16
% AirRep = repmat(db(Air), 1, 5);
% X = transpose(db(AN(:,12:16) - AirRep));
X = db(transpose(db(AN(:,12:16))));

[icasig, A, W] = fastica(X,'numOfIC',ICNo,'displayMode','off','firstEig',1,'lastEig',PCNo); % fast ICA

% [icasig, A, W] = fastica(X,'numOfIC',ICNo,'displayMode','off','interactivePCA','on'); % fast ICA
fasticag(X)  % used the correct number of PCA can always produce good results

figure;
for i=1:size(icasig,1)
    subplot(size(icasig,1),1,i),plot(icasig(i,:));
end

%% reconstruct signal
rejectICA=4;
% reconstruct the signal
A2=A;
icasig2=icasig;
A2(:,rejectICA)=[];
icasig2(rejectICA,:)=[];

newX=(A2*icasig2);

figure, 
for i=1:size(newX,1)
    subplot(size(newX,1),1,i)
%     plot(X,newX(i,:)),xlim([X(1) X(end)])
    plot(newX(i,:))  ;
end

%% decomposed signal of air
% FoamNail - FNonObj (Nail) -- 1 : Nail (0.280)
% Air -- 4 -- (0.372)
% AdhesiveOnly -- 1 -- (0.310)

%               1       2       3       4
%   Air        0.24    0.23    0.26    0.37        Air      -- 4 : 0.37
%   Adhesive   0.23    0.23    0.26    0.34        Adhesive -- 1 : 0.31
%   Nail       0.21    0.20    0.22    0.31        Nail     -- 1 : 0.27

%               1       2       3       4
%   Air        0.24    0.29    0.26    0.26        Air      -- 4 : 0.37
%   Adhesive   0.23    0.28    0.25    0.24        Adhesive -- 1 : 0.31
%   Nail       0.21    0.25    0.24    0.22        Nail     -- 1 : 0.27
No_Air = 1;
Comparison = Air;
% Comparison = FoamOnly;
ComparisonScale = (max(db(Comparison)) - min(db(Comparison)));
DecomposedScale = max(icasig(No_Air,:)) - min(icasig(No_Air,:));

Comparison = db(Comparison) - min(db(Comparison));
DecomposedNormalized = (icasig(No_Air,:) - min(icasig(No_Air,:))) / DecomposedScale * ComparisonScale;
% [C2, lag] = xcorr(DecomposedNormalized, Comparison, 'coeff');
[C2, lag] = xcorr((icasig(No_Air,:) - min(icasig(No_Air,:))) / (max(icasig(No_Air,:)) - min(icasig(No_Air,:))) * (max(db(Comparison)) - min(db(Comparison))), Comparison, 'coeff');
abs(corr(Comparison, transpose(DecomposedNormalized)))
max(abs(C2))

figure;
X = 1:201;
h1 = plot(X, DecomposedNormalized, 'r');
hold on;
h2 = plot(X, Comparison, '-.b');
hold off;
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');
legend([h1,h2],'Signal After ICA', 'Original Signal');

%% 
% [C2, lag] = xcorr((icasig(No_Air,:) - min(icasig(No_Air,:))) / (max(icasig(No_Air,:)) - min(icasig(No_Air,:))) , Comparison, 'coeff');
figure;
plot(lag, abs(C2));
title('Correlation(Air, Decomposed_Air)')
ylabel('Amplitude')
grid on
xtitle = 'Time';
xtitle = [xtitle newline newline 'Figure 4.7 Correlation(Air, Decomposed_air)'];
xlabel(xtitle, 'FontSize', 14);
ylabel('Correlation', 'FontSize', 14);

%% 