clear 
close all
clc
%%
% Load data
% I_REF=imread('dc.tif');
% I_REF=double(I_REF);
% I_REF=(I_REF(861:861+299,1:1+299,:));
% I_REF=double(I_REF);
% I_REF(I_REF<0)=0;
% I_REF = I_REF/max(I_REF(:));
% [Nn,Mm,Ll] = size(I_REF);

load 'C:\Users\dayan\Documents\dataset_UIS\Jorge\BGUDatabase\Cleaned\BGU_4.mat'
I_REF = rad(50:305,50:305,:);
I_REF = I_REF/max(I_REF(:));
[Nn,Mm,Ll] = size(I_REF);

% Number of grooves per mm
N = [200 300 400 500 600 700];

% sensor 
pix_sensor = (4.65)*1e-6;

%% Littrow configuration
G = N(end)/1e-3; 
m = -1;  % mode

lambda = linspace(400,700,Ll)';
mean_lambda = (lambda(1) + lambda(end))/2; 
mean_lambda = mean_lambda*1e-9;

angle_inp = rad2deg(asin((m*mean_lambda*G)/(2)));

% Calculate focal distance
LD = Ll.*pix_sensor;
d_focal = ((Ll-2).*pix_sensor).*cos(deg2rad(-angle_inp))./(G*((lambda(end)-lambda(1)).*1e-9));

%% check
d = 1e-3/N(end);
ref =  rad2deg(asin((sin(deg2rad(angle_inp)) - (m*(lambda).*1e-9)/d)));

%% grating
size_g  = 16;
rep     = 16; 

spatial_dis = round(rand(rep,rep)*(length(N)-1) + 1);
D_matrix_idx = kron(spatial_dis,ones(size_g, size_g));
figure,imagesc(D_matrix_idx)

%%
resd = zeros(1,6);
for n = 1:6
    
d = 1e-3/N(n);
angle_out = rad2deg(asin((sin(deg2rad(angle_inp)) - (m*(lambda*1e-9))/d)));

% gamm = angle_inp -  rad2deg(asin((sin(deg2rad(angle_inp)) - (m*(lambda*1e-9))/d)));
gamm = angle_inp -  rad2deg(asin((sin(deg2rad(angle_inp)) - (m*(mean_lambda))/d)));
gamm = gamm/2;
angle_out = angle_out - (2*angle_inp) + (2*gamm);

alt = d_focal*tan(deg2rad(angle_out));
resd(n) = min(floor(alt/pix_sensor));
alt = alt - d_focal*tan(deg2rad(angle_out(1)));

ind_pix = floor(alt/pix_sensor);
% ind_pix(end) = ind_pix(end) - 1;

figure;
subplot(1,2,1),plot(lambda, alt),title('Wavelength vs heights');
subplot(1,2,2),plot(lambda, ind_pix),title('Wavelength vs sensor element');

% grating matrix
D_aux = zeros(Ll, max(ind_pix)+1);
for i = 1:max(ind_pix)+1
    zeros_aux = zeros(Ll,1);
    posc = find(ind_pix == i-1);
    zeros_aux(posc) = 1;
    D_aux(:,i) = zeros_aux;
end
sum_D = 1./sum(D_aux);
D_aux = D_aux.*repmat(sum_D, Ll,1);
D_aux(isnan(D_aux)) = 0;
imagesc(D_aux)
pause(1)
D{n} = D_aux;
end
%%
%% Sensor
elements_sensor = Nn+Ll-1;
dim_sensor = elements_sensor*pix_sensor;
resd = resd - min(resd);

%%
sensor = zeros(Nn, elements_sensor);

% center = round((d_focal*tan(deg2rad(angle_out_s(1))))/pix_sensor);
figure;
for k = 1:6
    [row,col] = find(D_matrix_idx==k);
    ext = size(D{k},2);
    
    indx = [0:ext-1] + resd(k);
    indx_acc{k} = indx;
    % indx = [0:ext-1] + resd(k);
    for i=1:length(row)
        sensor(row(i),indx + col(i)) = sensor(row(i),indx + col(i)) + squeeze(I_REF(row(i),col(i),:))'*D{k};
    end
    c=gray;

    imagesc(sensor),title('Measurements');
    colormap(c);
end
%%
figure;
subplot(1,2,1),imagesc(I_REF(:,:,10)),title('Image');
subplot(1,2,2),imagesc(sensor),title('Measurements');
colormap(c);
%%
%% Vectorización de la imagen
t = 0;
for i  = 1:Nn
    for j = 1:Mm
        for k = 1:Ll
            t = t + 1;
            f(t) = I_REF(i,j,k);
        end
        
    end
end
%%
[Ns,Ms] = size(sensor); 
H = sparse(Ns*Ms,Nn*Mm*Ll);
D_ind_vect = D_matrix_idx' ;
D_ind_vect = D_ind_vect(:) ;
 %clear H
rowin = 0;
rowfin = 0;
y = [];
t = 0;
dif = 0;
for i = 1:Mm*Nn
    %condicional para pasar de fila
    if mod(i-1,Mm) == 0 && i>1
    t = t+1;
    dif = dif + t*Ms - rowin+min(indx_acc{D_ind_vect(i-1)});
    end
    
    rowin = dif + min(indx_acc{D_ind_vect(i)})+i; 
    rowfin = rowin + size(D{D_ind_vect(i)},2)-1;
    
    H(rowin:rowfin,Ll*(i-1)+1:i*Ll) =  D{D_ind_vect(i)}';

    disp(['Progreso' num2str(i/(Mm*Nn)*100) '%'])
end
y = H*f';


%%

%%
y = H*f';
t = 0;
for i = 1:Ns
    for j = 1:Ms
        t = t+1;
        Y(i,j) = y(t);
    end
end
%y = reshape(y(1:end-119),Ns,[]);
figure;
subplot(1,2,1),imshow(sensor./max(sensor(:))),title('Simulación');
subplot(1,2,2),imshow(Y./max(Y(:))),title('Producto matricial');
colormap(c)


save('H.mat','H')


