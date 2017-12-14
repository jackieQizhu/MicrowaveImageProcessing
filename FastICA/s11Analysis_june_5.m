% clear all;
% close all;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %               the following code is for bottle with adhesive and nail
% %               the nail is inside of the adhesive in the bottle
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename = 'Adhesive_nail_01.s2p';
% TDM = zeros(201, 25);
% TScale = linspace(0, 200/16, 201);
% ExScale = 1:25;
% for i = 1:9
%     filename(16) = num2str(i,2);
%     Fnames(i,:) = filename;
% end
% for i = 10:19
%     filename(15) = num2str(1,2);
%     filename(16) = num2str(i-10,2);
%     Fnames(i,:) = filename;
% end
% for i = 20:25
%     filename(15) = num2str(2,2);
%     filename(16) = num2str(i-20,2);
%     Fnames(i,:) = filename;
% end    
% for i = 1 : 25
%     name = Fnames(i,:);
%     dat = sparameters(name);
%     
%     s11 = dat.Parameters(1,1,:);
%     x = s11;
%     s11 = reshape(x,size(x,3),size(x,1)*size(x,2));
%     
%     td = s11;
%     td = ifft(s11);
%     TDM(:,i) = td;
% end

% dat = sparameters('mid_anb.s2p');
% s11 = dat.Parameters(1,1,:);
% x = s11;
% s11 = reshape(x,size(x,3),size(x,1)*size(x,2));
% mid_air = ifft(s11);
% save('MidAir.mat', 'mid_air');
% 
% TDM_File = matfile('AdhesiveNail.mat');
% TDM = TDM_File.TDM;
% figure
% imagesc(ExScale, TScale, db(TDM));
% colorbar;
% xtitle = 'Position';
% xtitle = [xtitle newline newline 'Figure 4.1 bottle-with-adhesive-nail-time-position'];
% xlabel(xtitle, 'FontSize', 14);
% caxis([-70 -30]);
% save('AdhesiveNail.mat', 'TDM');
% % title('Figure 4.1 bottle-with-adhesive-nail-time-position')
% grid;
% % 
% % sub = db((TDM) - (TDM(:,1)) * ones(1,25));
% % % figure;
% % % imagesc(ExScale, TScale, sub);
% % % colorbar;
% % % title('botthe-with-adhesive-nail-time-position-substract-center')
% % % caxis([-55 -35]);
% % % 
% % 
% % % read non-object data
dat = sparameters('Adhesive_no_object.s2p');
non_obj = dat.Parameters(1,1,:);
x = non_obj;
non_obj = reshape(x,size(x,3),size(x,1)*size(x,2));
non_obj = ifft(non_obj);
save('AdhesiveNonObject', 'non_obj');

% % cor_res = zeros(25,1);
% % for i = 1:25
% %     cor_res(i) = db(corr(TDM(:,i), TDM(:,i) - non_obj));
% % end
% % 
% % cor_res = - (cor_res / min(cor_res(:))) * 0.5;
% % 
% % figure;
% % stem(cor_res);
% % % title('Figure 4.4 Adhesive-Nail-substract-air-time-position-correlation-non-object')
% % xtitle = 'Position';
% % xtitle = [xtitle newline newline 'Figure 4.2 Adhesive-Nail-time-position-correlation-Adhesive-Nail'];
% % xlabel(xtitle, 'FontSize', 14);
% % ylabel('Correlation', 'FontSize', 14);
% 
filename = 'foam_adhensive_01.s2p';
FDM = zeros(201, 25);
TScale = linspace(0, 200/16, 201);
ExScale = 1:25;

for i = 1:9
    filename(17) = num2str(i,2);
    Fnames1(i,:) = filename;
end
for i = 10:19
    filename(16) = num2str(1,2);
    filename(17) = num2str(i-10,2);
    Fnames1(i,:) = filename;
end
for i = 20:25
    filename(16) = num2str(2,2);
    filename(17) = num2str(i-20,2);
    Fnames1(i,:) = filename;
end
for i = 1 : 25
    name = Fnames1(i,:);
    dat = sparameters(name);
    
    s21 = dat.Parameters(2,1,:);
    x = s21;
    s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
    
    td = ifft(s21);
    FDM(:,i) = td;
end

save('FoamAdhesive.mat', 'FDM');
% 
% [C, lag] = xcorr(FDM(:,1), TDM(:,1) - non_obj);
% figure;
% plot(lag,100*C,'k')
% ylabel('Amplitude')
% grid on
% % title('Figure 4.4 Adhesive-Nail-substract-air-time-position-correlation-non-object')
% xtitle = 'Time';
% xtitle = [xtitle newline newline 'Figure 4.8 Adhesive-Nail-time-position-correlation-Adhesive-Nail'];
% xlabel(xtitle, 'FontSize', 14);
% ylabel('Correlation', 'FontSize', 14);
% 
% % 
sub = db((TDM) - (mid_air) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
% pcolor(ExScale, TScale, sub);
%shading interp;
colorbar;
% title('Figure 4.3 botthe-with-adhesive-nail-time-position-substract-center-air')
xtitle = 'Position';
xtitle = [xtitle newline newline 'Figure 4.3 botthe-with-adhesive-nail-time-position-substract-center-air'];
xlabel(xtitle, 'FontSize', 14);
ylabel('Time Slot', 'FontSize', 14);
caxis([-60 -35]);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % %               the following code is for 
% % %               (adhesive substract air)
% % %               non object
% % %               correlation
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % % S(an~) * S(an(i))
% % cor_res = zeros(25,1);
% % for i = 1:25
% %     cor_res(i) = corr(sub(:,i), sub(:,i) - non_obj);
% % end
% 
% % figure;
% % stem(cor_res);
% % % title('Figure 4.4 Adhesive-Nail-substract-air-time-position-correlation-non-object')
% % xtitle = 'Position';
% % xtitle = [xtitle newline newline 'Figure 4.4 Adhesive-Nail-substract-air-time-position-correlation-Adhesive-Nail'];
% % xlabel(xtitle, 'FontSize', 14);
% % ylabel('Correlation', 'FontSize', 14);
% %
% 
% [C, lag] = xcorr(TDM(:,1), TDM(:,1) - non_obj - mid_air);
% figure;
% plot(lag,10*C,'k')
% ylabel('Amplitude')
% grid on
% % title('Figure 4.4 Adhesive-Nail-substract-air-time-position-correlation-non-object')
% xtitle = 'Time';
% xtitle = [xtitle newline newline 'Figure 4.4 Adhesive-Nail-time-position-correlation-Adhesive-Nail'];
% xlabel(xtitle, 'FontSize', 14);
% ylabel('Correlation', 'FontSize', 14);
% 
% 
% cor_res = zeros(25,1);
% for i = 1:25
%     cor_res(i) = corr(non_obj, TDM(:,i) - non_obj);
% end
% 
% figure;
% stem(cor_res);
% xtitle = 'Position';
% xtitle = [xtitle newline newline 'Figure 4.5 Adhesive-Nail-substract-air-time-position-correlation-Non-Object'];
% xlabel(xtitle, 'FontSize', 14);
% ylabel('Correlation', 'FontSize', 14);
% title('botthe-with-adhesive-nail-correlation-with-non-object');
% 
% % cor_res = zeros(25,1);
% % for i = 1:25
% %     cor_res(i) = db(corr(TDM(:,1), TDM(:,i)));
% % end
% % 
% % figure;
% % stem(cor_res);
% % title('botthe-with-adhesive-nail-correlation-with-bottle-edge')
% % 
% % cor_res = zeros(201,1);
% % non_obj_extend = zeros(201,25);
% % for i = 1:25
% %     non_obj_extend(:,i) = non_obj;
% % end
% % for i = 1:201
% %     cor_res(i) = db(corr(transpose(TDM(i,:)), transpose(non_obj_extend(i,:))));
% % end
% % 
% % figure;
% % title('botthe-with-adhesive-nail-correlation-with-non-object')
% % stem(cor_res);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %               the following code is for foam adhensive
% %               the adhesive fluid is poured into the foam
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% filename = 'foam_adhensive_01.s2p';
% TDM = zeros(201, 25);
% TScale = linspace(0, 200/16, 201);
% ExScale = 1:25;
% 
% for i = 1:9
%     filename(17) = num2str(i,2);
%     Fnames1(i,:) = filename;
% end
% for i = 10:19
%     filename(16) = num2str(1,2);
%     filename(17) = num2str(i-10,2);
%     Fnames1(i,:) = filename;
% end
% for i = 20:25
%     filename(16) = num2str(2,2);
%     filename(17) = num2str(i-20,2);
%     Fnames1(i,:) = filename;
% end
% for i = 1 : 25
%     name = Fnames1(i,:);
%     dat = sparameters(name);
%     
%     s21 = dat.Parameters(2,1,:);
%     x = s21;
%     s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
%     
%     td = ifft(s21);
% %     td = td ./ abs(td);
%     TDM(:,i) = td;
% end
% FAD = TDM;
% % figure
% % imagesc(ExScale, TScale, db(TDM));
% % colorbar;
% % caxis([-70 -30]);
% % title('foam-adsv-time-position')
% % grid;
% % 
% % sub = db(TDM) - db(TDM(:,13)) * ones(1,25);
% % figure;
% % imagesc(ExScale, TScale, sub);
% % colorbar;
% % title('foam-adsv-position-substract-center')
% % caxis([-30 20]);
% % 
% sub = db(TDM) - db(mid_air) * ones(1,25);
% figure;
% imagesc(ExScale, TScale, sub);
% colorbar;
% xtitle = 'Position';
% xtitle = [xtitle newline newline 'Figure 4.6 foam-adsv-position-substract-center-air'];
% xlabel(xtitle, 'FontSize', 14);
% ylabel('Reflection Time', 'FontSize', 14);
% % title('foam-adsv-position-substract-center-air')
% caxis([-30 20]);
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %
% % %               the following code is for foam nail
% % %               the nail is place inside of the foam
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
filename = 'foam_nail_01.s2p';
FNM = zeros(201, 25);
TScale = linspace(0, 200/16, 201);
ExScale = 1:25;

for i = 1:9
    filename(12) = num2str(i,2);
    Fnames2(i,:) = filename;
end
for i = 10:19
    filename(11) = num2str(1,2);
    filename(12) = num2str(i-10,2);
    Fnames2(i,:) = filename;
end
for i = 20:25
    filename(11) = num2str(2,2);
    filename(12) = num2str(i-20,2);
    Fnames2(i,:) = filename;
end  

for i = 1 : 25
    name = Fnames2(i,:);
    dat = sparameters(name);
    
    s21 = dat.Parameters(2,1,:);
    x = s21;
    s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
    
    td = ifft(s21);
    FNM(:,i) = td;
end

save('FoamNail.mat', 'FNM')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %               the following code is for adhesive non object
% %               correlation
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% read non-object data
dat = sparameters('foam_no_object.s2p');
non_obj = dat.Parameters(1,1,:);
x = non_obj;
non_obj = reshape(x,size(x,3),size(x,1)*size(x,2));
non_obj = ifft(non_obj);
save('FoamNonObj.mat', 'non_obj');

dat = sparameters('foam_air_2.s2p');
foam_air = dat.Parameters(1,1,:);
x = foam_air;
foam_air = reshape(x,size(x,3),size(x,1)*size(x,2));
foam_air = ifft(foam_air);
save('FoamAir.mat', 'non_obj');

sub = db(TDM) - db(foam_air) * ones(1,25);
% figure;
% imagesc(ExScale, TScale, sub);
% colorbar;
% title('Figure 4.5 foam-adsv-position-substract-center-air')
% caxis([-30 20]);

% S(an~) * S(an(i))
cor_res = zeros(25,1);
for i = 1:25
    cor_res(i) = corr(TDM(:,i), TDM(:,i) - non_obj - foam_air);
end

figure;
stem(cor_res);
xtitle = 'Position';
xtitle = [xtitle newline newline 'Figure 4.8 botthe-with-adhesive-foam-substract-foam-and-air-correlation-with-foam-nail'];
xlabel(xtitle, 'FontSize', 14);
ylabel('Correlation', 'FontSize', 14);
% title('Figure 4.6 botthe-with-adhesive-foam-substract-foam-and-air-correlation-with-foam');
% 
% cor_res = zeros(25,1);
% for i = 1:25
%     cor_res(i) = corr(foam_air, TDM(:,i) - non_obj - foam_air);
% end
% 
% % figure;
% % stem(cor_res);
% % % title('botthe-with-adhesive-foam-substract-foam-and-air-correlation-with-non-object')
% % xtitle = 'Position';
% % xtitle = [xtitle newline newline 'Figure 4.7 botthe-with-adhesive-foam-substract-foam-and-air-correlation-with-foam-air'];
% % xlabel(xtitle, 'FontSize', 14);
% % ylabel('Correlation', 'FontSize', 14);
% 
[C, lag] = xcorr(TDM(:,1), TDM(:,1) - foam_air);
figure;
plot(lag,10*C,'k')
ylabel('Amplitude')
grid on
% title('Figure 4.4 Adhesive-Nail-substract-air-time-position-correlation-non-object')
xtitle = 'Time';
xtitle = [xtitle newline newline 'Figure 4.7 Adhesive-Nail-time-position-correlation-Adhesive-Nail'];
xlabel(xtitle, 'FontSize', 14);
ylabel('Correlation', 'FontSize', 14);
% 
% 
% % cor_res = zeros(25,1);
% % for i = 1:25
% %     cor_res(i) = db(corr(TDM(:,1), TDM(:,i)));
% % end
% % 
% % figure;
% % title('botthe-with-adhesive-nail-correlation-with-edge-position')
% % stem(cor_res);
% 
% % figure
% % imagesc(ExScale, TScale, db(TDM));
% % colorbar;
% % caxis([-70 -30]);
% % title('foam-nail-time-position')
% % grid;
% 
% % sub = db(TDM) - db(TDM(:,13)) * ones(1,25);
% % figure;
% % imagesc(ExScale, TScale, sub);
% % colorbar;
% % title('foam-nail-position-substract-center')
% % caxis([-30 20]);
% % 
% % sub = db(TDM) - db(exp_generator) * ones(1,25);
% % figure;
% % imagesc(ExScale, TScale, sub);
% % colorbar;
% % title('foam-nail-position-substract-center-air')
% % caxis([-30 20]);
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % % 
% % % 
% % % 
% % 
% % 
% % 
% % How to signal processing about the features
% % 
% % 