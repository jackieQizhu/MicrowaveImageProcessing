clear all;
close all
filename = 'strips_01.s2p';
TDM = zeros(201, 25);
TScale = linspace(0, 200/16, 201);
ExScale = 1:25;
for i = 1:9
    filename(9) = num2str(i,2);
    Fnames(i,:) = filename;
end
for i = 10:19
    filename(8) = num2str(1,2);
    filename(9) = num2str(i-10,2);
    Fnames(i,:) = filename;
end
for i = 20:25
    filename(8) = num2str(2,2);
    filename(9) = num2str(i-20,2);
    Fnames(i,:) = filename;
end 


for i = 1 : 25
    name = Fnames(i,:);
    dat = sparameters(name);
    
    s21 = dat.Parameters(2,1,:);
    x = s21;
    s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
    
    td = ifft(s21);
    TDM(:,i) = td;
end

dat = sparameters('mid_air.s2p');
s21 = dat.Parameters(2,1,:);
x = s21;
s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
mid_air = ifft(s21);

figure
imagesc(ExScale, TScale, db(TDM));
colorbar;
caxis([-70 -30]);
title('strips-time-position')
grid;

sub = db(TDM) - db(TDM(:,13)) * ones(1,25);
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('strips-time-position-substract-center')
caxis([-30 20]);

sub = db(TDM) - db(mid_air) * ones(1,25);
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('strips-time-position-substract-center-air')
caxis([-30 20]);