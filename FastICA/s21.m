clear all;
close all
filename = 's21_00.s2p';
TDM = zeros(201, 19);
TScale = linspace(0, 200/16, 201);
ExScale = 1:19;

for i = 9:-1:0
    filename(6) = num2str(i,2);
    Fnames(10-i,:) = filename;
end

filename = 's21_10.s2p';
for i = 1:9
    filename(6) = num2str(i,2);
    Fnames(10+i,:) = filename;
end
    
for i = 1 : 19
    name = Fnames(i,:);
    dat = sparameters(name);
    
    s21 = dat.Parameters(2,1,:);
    x = s21;
    s21 = reshape(x,size(x,3),size(x,1)*size(x,2));
    
    td = ifft(s21);
    TDM(:,i) = td;
end

imagesc(ExScale, TScale, db(TDM));
grid;

sub = TDM - TDM(:,1) * ones(1,19);
figure;
imagesc(ExScale, TScale, db(sub));