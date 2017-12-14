clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               the following code is for strips
%               the metal strips attached on the foam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
dat = sparameters('mid_anb.s2p');
%dat = sparameters('mid_air.s2p');
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

sub = db((TDM) - (TDM(:,1)) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('strips-time-position-substract-first')
caxis([-30 20]);

sub = db((TDM) - (mid_air) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('strips-time-position-substract-center-air')
caxis([-30 20]);











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               the following code is for foam adhensive
%               the adhesive fluid is poured into the foam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename = 'foam_adhensive_01.s2p';
TDM = zeros(201, 25);
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
    TDM(:,i) = td;
end
figure
imagesc(ExScale, TScale, db(TDM));
colorbar;
caxis([-70 -30]);
title('foam-adsv-time-position')
grid;

sub = db((TDM) -(TDM(:,1)) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('foam-adsv-position-substract-first')
caxis([-60 -35]);

sub = db((TDM) -(mid_air) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('foam-adsv-position-substract-center-air')
caxis([-60 -40]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               the following code is for foam nail
%               the nail is place inside of the foam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename = 'foam_nail_01.s2p';
TDM = zeros(201, 25);
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
    TDM(:,i) = td;
end
figure
imagesc(ExScale, TScale, db(TDM));
colorbar;
caxis([-70 -30]);
title('foam-nail-time-position')
grid;

sub = db((TDM) -(TDM(:,11)) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('foam-nail-position-substract-first')
caxis([-30 20]);

sub = db((TDM) - (mid_air) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('foam-nail-position-substract-center-air')
caxis([-30 20]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               the following code is for bottle with adhesive and nail
%               the nail is inside of the adhesive in the bottle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Adhesive_nail_01.s2p';
TDM = zeros(201, 25);
TScale = linspace(0, 200/16, 201);
ExScale = 1:25;
for i = 1:9
    filename(16) = num2str(i,2);
    Fnames3(i,:) = filename;
end
for i = 10:19
    filename(15) = num2str(1,2);
    filename(16) = num2str(i-10,2);
    Fnames3(i,:) = filename;
end
for i = 20:25
    filename(15) = num2str(2,2);
    filename(16) = num2str(i-20,2);
    Fnames3(i,:) = filename;
end    
for i = 1 : 25
    name = Fnames3(i,:);
    dat = sparameters(name);
    
    s21 = dat.Parameters(2,1,:);
    x = s21;
    s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
    
    td = ifft(s21);
    TDM(:,i) = td;
end

figure
imagesc(ExScale, TScale, db(TDM));
colorbar;
caxis([-70 -30]);
title('bottle-with-adhesive-nail-time-position')
grid;

sub = db((TDM) - (TDM(:,1)) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
colorbar;
title('botthe-with-adhesive-nail-time-position-substract-first')
caxis([-50 -30]);

sub = db((TDM) - (mid_air) * ones(1,25));
figure;
imagesc(ExScale, TScale, sub);
% pcolor(ExScale, TScale, sub);
%shading interp;
colorbar;
title('botthe-with-adhesive-nail-time-position-substract-center-air')
caxis([-60 -35]);

