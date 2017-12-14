% 12/04 => 11:30 AM
clear all;

file = load('AdhesiveNail.mat');
AN = file.TDM;
AdhesiveNail = AN(:,13);


file = load('AdhesiveNonObject.mat');
AdhesiveOnly = file.non_obj;

file = load('Air.mat');
Air = file.mid_air;

file = load('FoamAdhesive.mat');
FAdhesive = file.FDM;


file = load('FoamAir.mat');
FAir = file.non_obj;


file = load('FoamNonObj.mat');
FNonObj = file.non_obj;


file = load('FoamNail.mat');
FN = file.FNM;
FoamNail = FN(:,13);

%% ifft data
AN = ifft(AN);
FN = ifft(FN);
AdhesiveNail = ifft(AdhesiveNail);
AdhesiveOnly = ifft(AdhesiveOnly);
Air = ifft(Air);
FAdhesive = ifft(FAdhesive);
FAir = ifft(FAir);
FNonObj = ifft(FNonObj);
FoamNail = ifft(FoamNail);

file = load('FoamOnlyExtraction.mat');
FoamOnlyExtraction = file.FoamOnlyExtraction;

file = load('NailOnlyExtraction.mat');
NailOnlyExtraction = file.NailOnlyExtraction;

file = load('NailOnlyExtraction1.mat');
NailOnlyExtraction1 = file.NailOnlyExtraction1;

file = load('NailOnlyExtraction2.mat');
NailOnlyExtraction2 = file.NailOnlyExtraction2;

%% Show input signal component
figure;
TOTAL = 4;
subplot(TOTAL,1,1), 
plot(transpose(db(FN(:,13))) - min(transpose(db(FN(:,13))))), title('Signal of FoamNail');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');
subplot(TOTAL,1,2), plot(db(Air) - min(db(Air))), title('Signal of Air');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')
subplot(TOTAL,1,3), plot(db(FNonObj) - min(db(FNonObj))),  title('Signal of FoamOnly');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')
subplot(TOTAL,1,4),
plot(db(FoamNail - FNonObj)),  title('Nail');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points')
%% Independent Components Analysis using FastICA
% addpath('/Users/qizhu/Documents/MATLAB/Matlab Code-2/FastICA_25');
PCNo = 4;
ICNo = 4;

X = db(transpose(db(AN(:,12:16))));
[icasig, A, W] = fastica(X,'numOfIC',ICNo,'displayMode','off','firstEig',1,'lastEig',PCNo); % fast ICA

% [icasig, A, W] = fastica(X,'numOfIC',ICNo,'displayMode','off','interactivePCA','on'); % fast ICA
fasticag(X)  % used the correct number of PCA can always produce good results

figure;
for i=1:size(icasig,1)
    subplot(size(icasig,1),1,i),plot(icasig(i,:));
end

% return
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

%%
% figure;
% subplot(2,1,1)
% plot(db(Air)), title('Original Air Signal');
% subplot(2,1,2)
% plot(icasig(3,:)), title('Signal Component After ICA Extraction');

%%
figure;
% for i=1:size(icasig,1)
%     subplot(size(icasig,1),1,i),plot((icasig(3,:) - min(icasig(3,:))) / (max(icasig(3,:)) - min(icasig(3,:))) * 40);
%     title('Decomposed Signal 1');
%     ylabel('Amplitude')
%     grid on 
%     xlabel('Data points');
% end
subplot(3,1,1)
plot((icasig(1,:) - min(icasig(1,:))) / (max(icasig(1,:)) - min(icasig(1,:))) * 40);
title('Decomposed Signal 1');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');

subplot(3,1,2)
plot((icasig(2,:) - min(icasig(2,:))) / (max(icasig(2,:)) - min(icasig(2,:))) * 40);
title('Decomposed Signal 2');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');

subplot(3,1,3)
plot((icasig(3,:) - min(icasig(3,:))) / (max(icasig(3,:)) - min(icasig(3,:))) * 40);
title('Decomposed Signal 3');
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');
%% decomposed signal of air
% FoamNail - FNonObj (Nail) -- 2 : Nail (0.2432)
% Air -- 4 -- (0.3718)
% AdhesiveOnly -- 1 -- (0.3135)
% transpose(NailOnlyExtraction)

No_Air = 4;
Comparison = transpose(NailOnlyExtraction2);

ComparisonScale = (max(db(Comparison)) - min(db(Comparison)));
DecomposedScale = (max(icasig(No_Air,:)) - min(icasig(No_Air,:)));

% NormalizedDecomposedSignal = (icasig(No_Air,:) - min(icasig(No_Air,:))) / DecomposedScale * ComparisonScale;

% NormalizedDecomposedSignal = (db(icasig(No_Air,:)) - min(db(icasig(No_Air,:)))) / DecomposedScale;
% NormalizedComparison = transpose((Comparison - min(db(Comparison))) / ComparisonScale);

% NormalizedDecomposedSignal = (icasig(No_Air,:) - min(icasig(No_Air,:))) / (max(icasig(No_Air,:)) - min(icasig(No_Air,:))) * (max(db(Comparison)) - min(db(Comparison)));
% NormalizedDecomposedSignal = transpose(NormalizedDecomposedSignal);
% NormalizedComparison = db(Comparison) - min(db(Comparison));

NormalizedDecomposedSignal = (icasig(No_Air,:) - min(icasig(No_Air,:))) / (max(icasig(No_Air,:)) - min(icasig(No_Air,:)));
NormalizedDecomposedSignal = transpose(NormalizedDecomposedSignal);
NormalizedComparison = (db(Comparison) - min(db(Comparison))) / (max(db(Comparison)) - min(db(Comparison)));

[C1, lag] = xcorr(NormalizedDecomposedSignal, NormalizedComparison);
[C2, lag] = xcorr(NormalizedComparison);
corr(NormalizedDecomposedSignal, NormalizedComparison)
% [~,I] = max(abs(C1));
% lagDiff = lag(I)
% max(abs(C1))
% max(abs(C1) ./ max(abs(C2)))

figure;
X = 1:201;
h1 = plot(X, (icasig(No_Air,:) - min(icasig(No_Air,:))) / (max(icasig(No_Air,:)) - min(icasig(No_Air,:))) * (max(db(Comparison)) - min(db(Comparison))), 'r');
hold on;
h2 = plot(X, db(Comparison) - min(db(Comparison)), 'k');
grid on
hold off;
ylabel('Amplitude(db)')
grid on 
xlabel('Data points');
legend([h1,h2],'Signal 3 of FoamNail after ICA', 'Original Signal of FoamOnly');

figure;
subplot(2,1,1);
bar(lag, C1);
title('cross correlation between FoamOnly and decomposed signal 3 of FoamNail after ICA')

subplot(2,1,2);
bar(lag, C2);
title('cross correlation between FoamOnly and FoamOnly')